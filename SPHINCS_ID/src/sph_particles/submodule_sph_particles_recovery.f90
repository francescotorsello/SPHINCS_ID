! File:         submodule_sph_particles_recovery.f90
! Authors:      Francesco Torsello (FT)
!************************************************************************
! Copyright (C) 2020, 2021, 2022 Francesco Torsello                     *
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

SUBMODULE (sph_particles) recovery

  !************************************************
  !
  !# This SUBMODULE contains the implementation
  !  of the method test_recovery of TYPE particles.
  !
  !  FT 18.02.2020
  !
  !************************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE test_recovery

    !************************************************
    !
    !# Tests the recovery. Computes the conserved
    !  variables from the physical ones, and then the
    !  physical ones from the conserved ones. It then
    !  compares the variables computed with the
    !  recovery PROCEDURES, with those computed with
    !  |sphincsid|. @todo add reference for recovery
    !
    !  FT 18.02.2020
    !
    !************************************************

    USE recovery,  ONLY: phys_2_cons, cons_2_phys
    USE tensor,    ONLY: jx, jy, jz
    USE constants, ONLY: zero

    IMPLICIT NONE

    INTEGER, PARAMETER:: unit_recovery= 34956

    DOUBLE PRECISION, DIMENSION(npart)  :: nlrf_rec
    DOUBLE PRECISION, DIMENSION(npart)  :: u_rec
    DOUBLE PRECISION, DIMENSION(npart)  :: pr_rec
    DOUBLE PRECISION, DIMENSION(3,npart):: vel_u_rec
    DOUBLE PRECISION, DIMENSION(npart)  :: theta_rec
    DOUBLE PRECISION, DIMENSION(npart)  :: nstar_rec
    DOUBLE PRECISION, DIMENSION(3,npart):: s_l_rec
    DOUBLE PRECISION, DIMENSION(npart)  :: e_hat_rec

    LOGICAL:: exist

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    ! Initialize local arrays
    nlrf_rec = zero
    u_rec    = zero
    pr_rec   = zero
    vel_u_rec= zero
    theta_rec= zero
    nstar_rec= zero
    s_l_rec  = zero
    e_hat_rec= zero

    !
    !-- Compute conserved fields from physical fields
    !
    CALL phys_2_cons( npart, nlrf, u, pr, vel_u, &
                      ! following is output
                      nstar_rec, s_l_rec, e_hat_rec )

    !
    !-- Recover physical fields from conserved fields
    !
  !  pr_rec= pr
  !
  !  CALL cons_2_phys( npart, nstar_rec, s_l_rec, e_hat_rec, &
  !                    ! following is output (pressure is INOUT)
  !                    nlrf_rec, vel_u_rec, u_rec, pr_rec, theta_rec )

    !
    !-- Print the original and recovered fields to formatted file
    !
    IF( PRESENT(namefile) )THEN
      finalnamefile= namefile
    ELSE
      finalnamefile= "recovery_test.dat"
    ENDIF

    INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

    IF( exist )THEN
      OPEN( UNIT= unit_recovery, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
            FORM= "FORMATTED", &
            POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
            IOMSG= err_msg )
    ELSE
      OPEN( UNIT= unit_recovery, FILE= TRIM(finalnamefile), STATUS= "NEW", &
            FORM= "FORMATTED", &
            ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening " // TRIM(finalnamefile), &
              ". The error message is", err_msg
      STOP
    ENDIF

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Values of the hydro fields computed by SPHINCS_ID and by the " &
    //"recovery, on the particles"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 1 in " // TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# column:      1        2       3       4       5", &
    "       6       7       8", &
    "       9       10      11", &
    "       12      13      14", &
    "       15      16      17", &
    "       18      19      20"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 2 in " // TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "#      particle      x [Msun_geo]       y [Msun_geo]       z [Msun_geo]", &
    "       local rest frame proper baryon density     " &
         //"recovered local rest frame proper baryon density", &
    "       local rest frame baryon density     " &
         //"recovered local rest frame baryon density", &
    "       specific internal energy     recovered specific internal energy", &
    "       pressure     recovered pressure", &
    "       x component of the computing frame velocity " &
         //"x component of the recovered computing frame velocity", &
    "       y component of the computing frame velocity " &
         //"y component of the recovered computing frame velocity", &
    "       z component of the computing frame velocity " &
         //"z component of the recovered computing frame velocity", &
    "       generalized Lorentz factor " &
         //"recovered generalized Lorentz factor"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 3 in " // TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF

    print_data_loop: DO itr = 1, npart, 1

      IF( THIS% export_form_xy .AND. &
          ( pos( 3, itr ) >=  0.5D0 .OR. &
            pos( 3, itr ) <= -0.5D0 ) &
      )THEN
        CYCLE
      ENDIF
      IF( THIS% export_form_x .AND. &
          ( pos( 3, itr ) >=  0.5D0 .OR. &
            pos( 3, itr ) <= -0.5D0 .OR. &
            pos( 2, itr ) >=  0.5D0 .OR. &
            pos( 2, itr ) <= -0.5D0 ) &
      )THEN
        CYCLE
      ENDIF
      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        itr, &                                                        ! 1
        pos      ( jx, itr ), &                                       ! 2
        pos      ( jy, itr ), &                                       ! 3
        pos      ( jz, itr ), &                                       ! 4
        nstar    ( itr ),     &                                       ! 5
        nstar_rec( itr ),     &                                       ! 6
        nlrf     ( itr ),     &                                       ! 7
        nlrf_rec ( itr ),     &                                       ! 8
        u        ( itr ),     &                                       ! 9
        u_rec    ( itr ),     &                                       ! 10
        pr       ( itr ),     &                                       ! 11
        pr_rec   ( itr ),     &                                       ! 12
        vel_u    ( jx, itr ), &                                       ! 13
        vel_u_rec( jx, itr ), &                                       ! 14
        vel_u    ( jy, itr ), &                                       ! 15
        vel_u_rec( jy, itr ), &                                       ! 16
        vel_u    ( jz, itr ), &                                       ! 17
        vel_u_rec( jz, itr ), &                                       ! 18
        theta    ( itr ),     &                                       ! 19
        theta_rec( itr )                                              ! 20

      IF( ios > 0 )THEN
        PRINT *, "...error when writing the arrays in " &
                 // TRIM(finalnamefile), ". The error message is", err_msg
        STOP
      ENDIF

    ENDDO print_data_loop

    PRINT *, " * Results from the recovery test printed to file ", finalnamefile
    PRINT *

    STOP

  END PROCEDURE test_recovery


END SUBMODULE recovery
