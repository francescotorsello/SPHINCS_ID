! File:         submodule_bssn_formulation_io.f90
! Authors:      Francesco Torsello (FT)
!************************************************************************
! Copyright (C) 2020-2023 Francesco Torsello                            *
!                                                                       *
! This file is part of SPHINCS_ID                                       *
!                                                                       *
! SPHINCS_ID is free software: you can redistribute it and/or modify    *
! it under the terms of the GNU General Public License as published by  *
! the Free Software Foundation, either version 3 of the License, or     *
! (at your option) any later version.                                   *
!                                                                       *
! SPHINCS_ID is distributed in the hope that it will be useful,         *
! but WITHOUT ANY WARRANTY; without even the implied warranty of        *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          *
! GNU General Public License for more details.                          *
!                                                                       *
! You should have received a copy of the GNU General Public License     *
! along with SPHINCS_ID. If not, see <https://www.gnu.org/licenses/>.   *
! The copy of the GNU General Public License should be in the file      *
! 'COPYING'.                                                            *
!************************************************************************

SUBMODULE (bssn_formulation) io

  !********************************************
  !
  !# This submodule contains the implementation
  !  of the methods of TYPE [[bssn]] that
  !  handle I/O (input/output)
  !
  !  FT 05.07.2022
  !
  !********************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE read_bssn_dump_print_formatted

    !************************************************
    !
    !# Read the BSSN ID from the binary file output
    !  by write_BSSN_dump, and print it to a
    !  formatted file
    !
    !  FT 08.02.2021
    !
    !************************************************

    USE mesh_refinement,  ONLY: levels, nlevels
    USE tensor,           ONLY: jxx, jxy, jxz, jyy, jyz, jzz, jx, jy, jz
    USE ADM_refine,       ONLY: lapse, shift_u, &
                                allocate_ADM, deallocate_ADM
    USE BSSN_refine,  ONLY: allocate_BSSN, deallocate_BSSN, &
                            Gamma_u,        & ! Conformal connection
                            phi,            & ! Conformal factor
                            trK,            & ! Trace of extrinsic curvature
                            A_BSSN3_ll,     & ! Conformal traceless
                                              ! extrinsic curvature
                            g_BSSN3_ll,     & ! Conformal metric
                            !Theta_Z4,       & ! Vector in the CCZ4 formulation.
                                              ! Loaded here because ADM_TO_BSSN
                                              ! calls SUBROUTINES that need it
                                              ! as input; however, it is not
                                              ! evolved in BSSN
                            !lapse_A_BSSN,   & ! Time derivative of lapse
                            !shift_B_BSSN_u, & ! Time derivativeof shift
                            read_BSSN_dump

    IMPLICIT NONE

    INTEGER:: i, j, k, l, min_ix_y, min_iy_y, min_iz_y, &
              min_ix_z, min_iy_z, min_iz_z

    DOUBLE PRECISION:: min_abs_y, min_abs_z

    LOGICAL:: exist

    CHARACTER(LEN= :), ALLOCATABLE:: finalnamefile

    PRINT *, "** Executing the read_bssn_dump_print_formatted subroutine..."

    levels = this% levels
    nlevels= this% nlevels

    CALL allocate_ADM()
    CALL allocate_BSSN()

    CALL read_BSSN_dump( 00000, namefile_bin )

    IF( this% call_flag == 0 )THEN
      PRINT *, "** The SUBROUTINE print_formatted_id_bssn_variables ", &
        " must be called after compute_and_print_bssn_variables, otherwise", &
        " there are no bssn fields to export to the formatted file."
      PRINT *, "   Aborting."
      PRINT *
      STOP
    ENDIF

    IF( PRESENT(namefile) )THEN
      finalnamefile= namefile
    ELSE
      finalnamefile= "bssn_vars.dat"
    ENDIF

    INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

    IF( exist )THEN
      OPEN( UNIT= 20, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
            FORM= "FORMATTED", &
            POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
            IOMSG= err_msg )
    ELSE
      OPEN( UNIT= 20, FILE= TRIM(finalnamefile), STATUS= "NEW", &
      FORM= "FORMATTED", &
            ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening ", TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when opening " &
    !         // TRIM(finalnamefile) )

    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id
    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Values of the fields (including coordinates) exported by LORENE "&
    // "on each grid point"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 1 in ", TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 1 in "&
    !         // TRIM(finalnamefile) )
    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# column:      1        2       3       4       5", &
    "       6       7       8", &
    "       9       10      11", &
    "       12      13      14", &
    "       15      16      17      18      19", &
    "       20      21      22", &
    "       23      24      25"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 2 in ", TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 2 in "&
    !        // TRIM(finalnamefile) )
    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "#      refinement level    x [km]       y [km]       z [km]       lapse", &
    "       shift_x [c]    shift_y [c]    shift_z [c]", &
    "       conformal factor phi        trace of extr. curv. trK", &
    "       g_BSSN_xx       g_BSSN_xy      g_BSSN_xz", &
    "       g_BSSN_yy       g_BSSN_yz      g_BSSN_zz", &
    "       A_BSSN_xx       A_BSSN_xy      A_BSSN_xz", &
    "       A_BSSN_yy       A_BSSN_yz      A_BSSN_zz", &
    "       Gamma_u_x       Gamma_u_y      Gamma_u_z"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 3 in ", TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 3 in "&
    !        // TRIM(finalnamefile) )

    DO l= 1, this% nlevels, 1

      ASSOCIATE( lapse      => lapse% levels(l)% var, &
                 shift_u    => shift_u% levels(l)% var, &
                 g_BSSN3_ll => g_BSSN3_ll% levels(l)% var, &
                 A_BSSN3_ll => A_BSSN3_ll% levels(l)% var, &
                 phi        => phi% levels(l)% var, &
                 trK        => trK% levels(l)% var, &
                 Gamma_u    => Gamma_u% levels(l)% var &
      )

        min_abs_y= 1D+20
        min_abs_z= 1D+20
        DO k= 1, this% get_ngrid_z(l), 1
          DO j= 1, this% get_ngrid_y(l), 1
            DO i= 1, this% get_ngrid_x(l), 1

              IF( ABS( this% coords% levels(l)% var( i, j, k, jy ) ) &
                  < min_abs_y )THEN
                min_abs_y= ABS( this% coords% levels(l)% var( i, j, k, jy ) )
                min_ix_y= i
                min_iy_y= j
                min_iz_y= k
              ENDIF

              IF( ABS( this% coords% levels(l)% var( i, j, k, jz ) ) &
                  < min_abs_z )THEN
                min_abs_z= ABS( this% coords% levels(l)% var( i, j, k, jz ) )
                min_ix_z= i
                min_iy_z= j
                min_iz_z= k
              ENDIF

            ENDDO
          ENDDO
        ENDDO

        DO k= 1, this% get_ngrid_z(l), 1

          DO j= 1, this% get_ngrid_y(l), 1

            DO i= 1, this% get_ngrid_x(l), 1

              IF( this% export_form_xy .AND. &
                  ( this% coords% levels(l)% var( i, j, k, jz ) /= &
                    this% coords% levels(l)% var( min_ix_z, min_iy_z, &
                                                  min_iz_z, jz ) ) )THEN
                CYCLE
              ENDIF
              IF( this% export_form_x .AND. &
                  ( this% coords% levels(l)% var( i, j, k, jz ) /= &
                    this% coords% levels(l)% var( min_ix_z, min_iy_z, &
                                                  min_iz_z, jz ) &
                    .OR. &
                    this% coords% levels(l)% var( i, j, k, jy ) /= &
                    this% coords% levels(l)% var( min_ix_y, min_iy_y, &
                                                  min_iz_y, jy ) ) )THEN
                CYCLE
              ENDIF

              WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * )&
                  l, &                                              ! 1
                  this% coords% levels(l)% var( i, j, k, jx ), &    ! 2
                  this% coords% levels(l)% var( i, j, k, jy ), &    ! 3
                  this% coords% levels(l)% var( i, j, k, jz ), &    ! 4
                  lapse( i, j, k ), &                               ! 5
                  shift_u( i, j, k, jx ), &                         ! 6
                  shift_u( i, j, k, jy ), &                         ! 7
                  shift_u( i, j, k, jz ), &                         ! 8
                  phi( i, j, k ), &                                 ! 9
                  trK( i, j, k ), &                                 ! 10
                  g_BSSN3_ll( i, j, k, jxx ), &                     ! 11
                  g_BSSN3_ll( i, j, k, jxy ), &                     ! 12
                  g_BSSN3_ll( i, j, k, jxz ), &                     ! 13
                  g_BSSN3_ll( i, j, k, jyy ), &                     ! 14
                  g_BSSN3_ll( i, j, k, jyz ), &                     ! 15
                  g_BSSN3_ll( i, j, k, jzz ), &                     ! 16
                  A_BSSN3_ll( i, j, k, jxx ), &                     ! 17
                  A_BSSN3_ll( i, j, k, jxy ), &                     ! 18
                  A_BSSN3_ll( i, j, k, jxz ), &                     ! 19
                  A_BSSN3_ll( i, j, k, jyy ), &                     ! 20
                  A_BSSN3_ll( i, j, k, jyz ), &                     ! 21
                  A_BSSN3_ll( i, j, k, jzz ), &                     ! 22
                  Gamma_u( i, j, k, jx ), &                         ! 23
                  Gamma_u( i, j, k, jy ), &                         ! 24
                  Gamma_u( i, j, k, jz )                            ! 25

              IF( ios > 0 )THEN
                PRINT *, "...error when writing the arrays in ", &
                         TRIM(finalnamefile), ". The error message is", err_msg
                STOP
              ENDIF
              !CALL test_status( ios, err_msg, "...error when writing " &
              !                  // "the arrays in " // TRIM(namefile) )

            ENDDO
          ENDDO
        ENDDO
      END ASSOCIATE
    ENDDO

    CLOSE( UNIT= 20 )

    !
    !-- Deallocate MODULE variables
    !
    CALL deallocate_ADM()
    CALL deallocate_BSSN()

    PRINT *, " * BSSN ID on the refined mesh, to be supplied to ", &
             "SPHINCS_BSSN, printed to formatted file ", TRIM(namefile)

    PRINT *, "** Subroutine read_bssn_dump_print_formatted " &
             // "executed."
    PRINT *

  END PROCEDURE read_bssn_dump_print_formatted


  MODULE PROCEDURE print_formatted_id_bssn_variables

    !************************************************
    !
    !# Print the BSSN ID, computed on the gravity
    !  grid, to a formatted file
    !
    !  FT 26.10.2020
    !
    !************************************************

    USE tensor,              ONLY: jxx, jxy, jxz, jyy, jyz, jzz, jx, jy, jz

    IMPLICIT NONE

    INTEGER:: i, j, k, l, min_ix_y, min_iy_y, min_iz_y, &
              min_ix_z, min_iy_z, min_iz_z

    DOUBLE PRECISION:: min_abs_y, min_abs_z

    LOGICAL:: exist

    CHARACTER(LEN= :), ALLOCATABLE:: finalnamefile

    ! Being abs_grid a local array, it is good practice to allocate it on the
    ! heap, otherwise it will be stored on the stack which has a very limited
    ! size. This results in a segmentation fault.
    !ALLOCATE( abs_grid( 3, this% ngrid_x, this% ngrid_y, this% ngrid_z ) )

    PRINT *, "** Executing the print_formatted_id_BSSN_variables " &
             // "subroutine..."

    IF( this% call_flag == 0 )THEN
      PRINT *, "** The SUBROUTINE print_formatted_id_bssn_variables ", &
        " must be called after compute_and_print_bssn_variables, otherwise", &
        " there are no bssn fields to export to the formatted file."
      PRINT *, "   Aborting."
      PRINT *
      STOP
    ENDIF

    IF( PRESENT(namefile) )THEN
      finalnamefile= namefile
    ELSE
      finalnamefile= "bssn-id.dat"
    ENDIF

    INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

    IF( exist )THEN
      OPEN( UNIT= 20, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
            FORM= "FORMATTED", &
            POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
            IOMSG= err_msg )
    ELSE
      OPEN( UNIT= 20, FILE= TRIM(finalnamefile), STATUS= "NEW", &
      FORM= "FORMATTED", &
            ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening ", TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when opening " &
    !         // TRIM(finalnamefile) )

    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id
    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Values of the fields (including coordinates) exported by LORENE "&
    // "on each grid point"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 1 in ", TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 1 in "&
    !         // TRIM(finalnamefile) )
    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# column:      1        2       3       4       5", &
    "       6       7       8", &
    "       9       10      11", &
    "       12      13      14", &
    "       15      16      17      18      19", &
    "       20      21      22", &
    "       23      24      25"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 2 in ", TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 2 in "&
    !        // TRIM(finalnamefile) )
    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "#      refinement level    x [km]       y [km]       z [km]       lapse", &
    "       shift_x [c]    shift_y [c]    shift_z [c]", &
    "       conformal factor phi        trace of extr. curv. trK", &
    "       g_BSSN_xx       g_BSSN_xy      g_BSSN_xz", &
    "       g_BSSN_yy       g_BSSN_yz      g_BSSN_zz", &
    "       A_BSSN_xx       A_BSSN_xy      A_BSSN_xz    ", &
    "       A_BSSN_yy       A_BSSN_yz      A_BSSN_zz", &
    "       Gamma_u_x       Gamma_u_y      Gamma_u_z"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 3 in ", TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 3 in "&
    !        // TRIM(finalnamefile) )

    DO l= 1, this% nlevels, 1

      ASSOCIATE( lapse           => this% lapse% levels(l)% var, &
                 shift_u         => this% shift_u% levels(l)% var, &
                 phi             => this% phi% levels(l)% var, &
                 trK             => this% trK% levels(l)% var, &
                 g_BSSN3_ll      => this% g_BSSN3_ll% levels(l)% var, &
                 A_BSSN3_ll      => this% A_BSSN3_ll% levels(l)% var, &
                 Gamma_u         => this% Gamma_u% levels(l)% var &
      )


        !DO iz= 1, this% ngrid_z, 1
        !  DO iy= 1, this% ngrid_y, 1
        !    DO ix= 1, this% ngrid_x, 1
        !      abs_grid( 1, ix, iy, iz )= ABS( this% grid( 1, ix, iy, iz ) )
        !      abs_grid( 2, ix, iy, iz )= ABS( this% grid( 2, ix, iy, iz ) )
        !      abs_grid( 3, ix, iy, iz )= ABS( this% grid( 3, ix, iy, iz ) )
        !    ENDDO
        !  ENDDO
        !ENDDO

        min_abs_y= HUGE(1.0D0)
        min_abs_z= HUGE(1.0D0)
        DO k= 1, this% get_ngrid_z(l), 1
          DO j= 1, this% get_ngrid_y(l), 1
            DO i= 1, this% get_ngrid_x(l), 1

              IF( ABS( this% coords% levels(l)% var( i, j, k, jy ) ) &
                  < min_abs_y )THEN
                min_abs_y= ABS( this% coords% levels(l)% var( i, j, k, jy ) )
                min_ix_y= i
                min_iy_y= j
                min_iz_y= k
              ENDIF

              IF( ABS( this% coords% levels(l)% var( i, j, k, jz ) ) &
                  < min_abs_z )THEN
                min_abs_z= ABS( this% coords% levels(l)% var( i, j, k, jz ) )
                min_ix_z= i
                min_iy_z= j
                min_iz_z= k
              ENDIF

            ENDDO
          ENDDO
        ENDDO

        DO k= 1, this% get_ngrid_z(l), 1

          DO j= 1, this% get_ngrid_y(l), 1

            DO i= 1, this% get_ngrid_x(l), 1

              IF( this% export_form_xy .AND. &
                  ( this% coords% levels(l)% var( i, j, k, jz ) /= &
                    this% coords% levels(l)% var( min_ix_z, min_iy_z, &
                                                  min_iz_z, jz ) ) )THEN
                CYCLE
              ENDIF
              IF( this% export_form_x .AND. &
                  ( this% coords% levels(l)% var( i, j, k, jz ) /= &
                    this% coords% levels(l)% var( min_ix_z, min_iy_z, &
                                                  min_iz_z, jz ) &
                    .OR. &
                    this% coords% levels(l)% var( i, j, k, jy ) /= &
                    this% coords% levels(l)% var( min_ix_y, min_iy_y, &
                                                  min_iz_y, jy ) ) )THEN
                CYCLE
              ENDIF

              WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * )&
                l, &                                              ! 1 !
                this% coords% levels(l)% var( i, j, k, jx ), &    ! 2 !
                this% coords% levels(l)% var( i, j, k, jy ), &    ! 3
                this% coords% levels(l)% var( i, j, k, jz ), &    ! 4
                lapse( i, j, k ), &                               ! 5
                shift_u( i, j, k, jx ), &                         ! 6
                shift_u( i, j, k, jy ), &                         ! 7
                shift_u( i, j, k, jz ), &                         ! 8
                phi( i, j, k ), &                                 ! 9
                trK( i, j, k ), &                                 ! 10
                g_BSSN3_ll( i, j, k, jxx ), &                     ! 11
                g_BSSN3_ll( i, j, k, jxy ), &                     ! 12
                g_BSSN3_ll( i, j, k, jxz ), &                     ! 13
                g_BSSN3_ll( i, j, k, jyy ), &                     ! 14
                g_BSSN3_ll( i, j, k, jyz ), &                     ! 15
                g_BSSN3_ll( i, j, k, jzz ), &                     ! 16
                A_BSSN3_ll( i, j, k, jxx ), &                     ! 17
                A_BSSN3_ll( i, j, k, jxy ), &                     ! 18
                A_BSSN3_ll( i, j, k, jxz ), &                     ! 19
                A_BSSN3_ll( i, j, k, jyy ), &                     ! 20
                A_BSSN3_ll( i, j, k, jyz ), &                     ! 21
                A_BSSN3_ll( i, j, k, jzz ), &                     ! 22
                Gamma_u( i, j, k, jx ), &                         ! 23
                Gamma_u( i, j, k, jy ), &                         ! 24
                Gamma_u( i, j, k, jz )                            ! 25

              IF( ios > 0 )THEN
                PRINT *, "...error when writing the arrays in ", &
                         TRIM(finalnamefile), ". The error message is", err_msg
                STOP
              ENDIF
              !CALL test_status( ios, err_msg, "...error when writing " &
              !                  // "the arrays in " // TRIM(finalnamefile) )

            ENDDO
          ENDDO
        ENDDO
      END ASSOCIATE
    ENDDO

    CLOSE( UNIT= 20 )

    PRINT *, " * BSSN ID on the gravity grid saved to formatted " &
             // "file ", TRIM(finalnamefile)

    PRINT *, "** Subroutine print_formatted_id_BSSN_variables " &
             // "executed."
    PRINT *

  END PROCEDURE print_formatted_id_bssn_variables


END SUBMODULE io
