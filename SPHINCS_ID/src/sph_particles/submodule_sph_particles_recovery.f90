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

    USE recovery,             ONLY: phys_2_cons, cons_2_phys
    USE tensor,               ONLY: jx, jy, jz, n_sym4x4
    USE constants,            ONLY: zero, one
    USE deactivate_particles, ONLY: nlrf_fb, u_fb, pr_fb, vel_u_fb, theta_fb, &
                                    cs_fb
    USE metric_on_particles,  ONLY: allocate_metric_on_particles, &
                                    deallocate_metric_on_particles, &
                                    g4_ll
    USE utility,              ONLY: compute_g4, determinant_sym4x4
    USE tmp,                  ONLY: fill_arrays

    IMPLICIT NONE

    INTEGER, PARAMETER:: unit_recovery= 34956

    INTEGER:: i_matter, a, a_max

    DOUBLE PRECISION:: det, p_max

    DOUBLE PRECISION, DIMENSION(npart)  :: nlrf_rec
    DOUBLE PRECISION, DIMENSION(npart)  :: u_rec
    DOUBLE PRECISION, DIMENSION(npart)  :: pr_rec
    DOUBLE PRECISION, DIMENSION(3,npart):: vel_u_rec
    DOUBLE PRECISION, DIMENSION(npart)  :: theta_rec
    DOUBLE PRECISION, DIMENSION(npart)  :: nstar_rec
    DOUBLE PRECISION, DIMENSION(3,npart):: s_l_rec
    DOUBLE PRECISION, DIMENSION(npart)  :: e_hat_rec

    LOGICAL:: exist

    CHARACTER( LEN= 2 ):: i_mat
    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    LOGICAL, PARAMETER:: debug= .FALSE.

    PRINT *, " * Testing recovery..."
    PRINT *

    IF( .NOT.ALLOCATED(g4_ll) )THEN
      CALL allocate_metric_on_particles(this% npart)
    ENDIF

    DO a= 1, this% npart, 1

      CALL compute_g4( this% lapse(a), &
            [this% shift_x(a), this% shift_y(a), this% shift_z(a)], &
            [this% g_xx(a), this% g_xy(a), this% g_xz(a), &
             this% g_yy(a), this% g_yz(a), this% g_zz(a)], &
             g4_ll(1:n_sym4x4,a) )

      CALL determinant_sym4x4( g4_ll(1:n_sym4x4,a), det )
      IF( ABS(det) < 1D-10 )THEN
          PRINT *, "** ERROR! The determinant of the spacetime metric is " &
                   // "effectively 0 at particle ", a
          STOP
      ELSEIF( det > 0 )THEN
          PRINT *, "** ERROR! The determinant of the spacetime metric is " &
                   // "positive at particle ", a
          STOP
      ENDIF

    ENDDO

    ALLOCATE( nlrf_fb (npart) )
    ALLOCATE( u_fb    (npart) )
    ALLOCATE( pr_fb   (npart) )
    ALLOCATE( vel_u_fb(3,npart) )
    ALLOCATE( theta_fb(npart) )
    ALLOCATE( cs_fb(npart) )
    nlrf_fb = nlrf
    u_fb    = u
    pr_fb   = pr
    vel_u_fb= vel_u(1:3,:)
    theta_fb= theta
    ! TODO: set the sound speed properly. From pwp_eos MODULE:
    ! enth= 1.0D0 + u + rho_rest/P
    ! cs=   SQRT((Gamma*P_cold + Gamma_th*P_th)/(rho_rest*enth))
    cs_fb   = one


    ! Initialize local arrays
    nlrf_rec = zero
    u_rec    = zero
    pr_rec   = zero
    vel_u_rec= zero
    theta_rec= zero
    nstar_rec= zero
    s_l_rec  = zero
    e_hat_rec= zero

    IF( debug ) PRINT *, "1"

    !
    !-- Compute conserved fields from physical fields
    !
    CALL phys_2_cons( npart, nlrf, u, pr, vel_u, &
                      ! following is output
                      nstar_rec, s_l_rec, e_hat_rec )

    IF( debug ) PRINT *, "2"

    !-------------------------------------------------------
    !-----DEBUGGING
    !-------------------------------------------------------
    p_max= 0.0D0
    DO a= 1, npart, 1
      IF( pr(a) > p_max )THEN
        p_max= pr(a)
        a_max= a
      ENDIF
    ENDDO
    CALL fill_arrays( npart, nstar, s_l_rec, e_hat_rec, nlrf, u, &
                      pr, vel_u, theta )
    !-------------------------------------------------------
    !-------------------------------------------------------

    !
    !-- Recover physical fields from conserved fields
    !
    pr_rec= pr

    CALL cons_2_phys( npart, nstar_rec, s_l_rec, e_hat_rec, &
                      ! following is output (pressure is INOUT)
                      nlrf_rec, vel_u_rec, u_rec, pr_rec, theta_rec, a_max )

    IF( debug ) PRINT *, "3"

    DEALLOCATE( nlrf_fb )
    DEALLOCATE( u_fb )
    DEALLOCATE( pr_fb )
    DEALLOCATE( vel_u_fb )
    DEALLOCATE( theta_fb )
    DEALLOCATE( cs_fb )

    CALL deallocate_metric_on_particles

    IF( debug ) PRINT *, "4"

    matter_objects_loop: DO i_matter= 1, this% n_matter, 1

      !PRINT *, " * Testing recovery on matter object", i_matter, "..."
      !PRINT *

      IF( i_matter > 9 )THEN
        WRITE( i_mat, "(I2)" ) i_matter
      ELSE
        WRITE( i_mat, "(I1)" ) i_matter
      ENDIF
      finalnamefile= "recovery_test-"//TRIM(i_mat)//".dat"

      ASSOCIATE( npart_in   => this% npart_i(i_matter-1) + 1, &
                 npart_fin  => this% npart_i(i_matter-1) +    &
                               this% npart_i(i_matter) )

      !
      !-- Print the original and recovered fields to formatted file
      !
      !IF( PRESENT(namefile) )THEN
      !  finalnamefile= namefile
      !ELSE
      !  finalnamefile= "recovery_test.dat"
      !ENDIF

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

      WRITE( UNIT = unit_recovery, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id

      WRITE( UNIT = unit_recovery, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "# Values of the hydro fields computed by SPHINCS_ID and by the " &
      //"recovery, on the particles"
      IF( ios > 0 )THEN
        PRINT *, "...error when writing line 1 in " // TRIM(finalnamefile), &
                 ". The error message is", err_msg
        STOP
      ENDIF

      WRITE( UNIT = unit_recovery, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
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

      WRITE( UNIT = unit_recovery, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
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

      print_data_loop: DO a = npart_in, npart_fin, 1

        IF( THIS% export_form_xy .AND. &
            ( pos( 3, a ) >=  0.5D0 .OR. &
              pos( 3, a ) <= -0.5D0 ) &
        )THEN
          CYCLE
        ENDIF
        IF( THIS% export_form_x .AND. &
            ( pos( 3, a ) >=  0.5D0 .OR. &
              pos( 3, a ) <= -0.5D0 .OR. &
              pos( 2, a ) >=  0.5D0 .OR. &
              pos( 2, a ) <= -0.5D0 ) &
        )THEN
          CYCLE
        ENDIF

        WRITE( UNIT = unit_recovery, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
          a, &                                                        ! 1
          pos      ( jx, a ) - this% barycenter(i_matter, jx), &      ! 2
          pos      ( jy, a ) - this% barycenter(i_matter, jy), &      ! 3
          pos      ( jz, a ) - this% barycenter(i_matter, jz), &      ! 4
          nstar    ( a ),     &                                       ! 5
          nstar_rec( a ),     &                                       ! 6
          nlrf     ( a ),     &                                       ! 7
          nlrf_rec ( a ),     &                                       ! 8
          u        ( a ),     &                                       ! 9
          u_rec    ( a ),     &                                       ! 10
          pr       ( a ),     &                                       ! 11
          pr_rec   ( a ),     &                                       ! 12
          vel_u    ( jx, a ), &                                       ! 13
          vel_u_rec( jx, a ), &                                       ! 14
          vel_u    ( jy, a ), &                                       ! 15
          vel_u_rec( jy, a ), &                                       ! 16
          vel_u    ( jz, a ), &                                       ! 17
          vel_u_rec( jz, a ), &                                       ! 18
          theta    ( a ),     &                                       ! 19
          theta_rec( a )                                              ! 20

        IF( ios > 0 )THEN
          PRINT *, "...error when writing the arrays in " &
                   // TRIM(finalnamefile), ". The error message is", err_msg
          STOP
        ENDIF

      ENDDO print_data_loop

      CLOSE( unit_recovery )

      END ASSOCIATE

      PRINT *, " * Results from the recovery test on matter object ", &
               i_matter, " printed to file ", finalnamefile
      PRINT *

    ENDDO matter_objects_loop

    !STOP

  END PROCEDURE test_recovery


END SUBMODULE recovery
