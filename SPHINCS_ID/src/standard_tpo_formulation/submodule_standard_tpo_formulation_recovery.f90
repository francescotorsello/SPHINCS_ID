! File:         submodule_standard_tpo_formulation_recovery.f90
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

SUBMODULE (standard_tpo_formulation) recovery_m2p

  !************************************************
  !
  !# This SUBMODULE contains the implementation
  !  of the method test_recovery of TYPE particles.
  !
  !  FT 25.02.2020
  !
  !************************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE test_recovery_m2p

    !************************************************
    !
    !# Tests the recovery using the metric mapped
    !  from the emsh to the particles. Computes the
    !  conserved variables from the physical ones,
    !  and then the physical ones from the conserved
    !  ones. It then compares the variables computed
    !  with the recovery PROCEDURES, with those
    !  computed with |sphincsid|.
    !  @todo add reference for recovery
    !
    !  FT 25.02.2020
    !
    !************************************************

    USE recovery,             ONLY: phys_2_cons, cons_2_phys
    USE tensor,               ONLY: jx, jy, jz, n_sym4x4
    USE constants,            ONLY: zero, one
    USE deactivate_particles, ONLY: nlrf_fb, u_fb, pr_fb, vel_u_fb, theta_fb, &
                                    cs_fb
    USE metric_on_particles,  ONLY: allocate_metric_on_particles, &
                                    deallocate_metric_on_particles, &
                                    get_metric_on_particles, g4_ll
    !USE map_metric_2_particles_refine, &
    !                          ONLY: update_ADM_metric_on_particles
    USE ADM_refine,           ONLY: lapse, shift_u, g_phys3_ll, &
                                    allocate_ADM, deallocate_ADM
    USE BSSN_refine,          ONLY: allocate_BSSN, deallocate_BSSN
    USE mesh_refinement,      ONLY: nlevels, levels, rad_coord, coords, &
                                    allocate_grid_function, &
                                    deallocate_grid_function
    USE gradient,             ONLY: allocate_gradient, deallocate_gradient
    USE McLachlan_refine,     ONLY: allocate_Ztmp, deallocate_Ztmp
    USE alive_flag,           ONLY: alive
    USE input_output,         ONLY: read_options
    USE units,                ONLY: set_units
    USE sph_variables,        ONLY: allocate_SPH_memory, &
                                    deallocate_SPH_memory
    USE utility,              ONLY: compute_g4, determinant_sym4x4

    IMPLICIT NONE

    INTEGER, PARAMETER:: unit_recovery= 34156

    INTEGER:: npart, i_matter, a, l

    !DOUBLE PRECISION:: det

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nlrf_rec
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: u_rec
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pr_rec
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: vel_u_rec
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: theta_rec
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar_rec
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: s_l_rec
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: e_hat_rec

    !DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: lapse_parts
    !DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: shift_parts
    !DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: g3_parts

    LOGICAL:: exist

    CHARACTER( LEN= 2 ):: i_mat
    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    LOGICAL, PARAMETER:: debug= .FALSE.

    npart= parts% get_npart()
    pos  = parts% get_pos()

    ALLOCATE ( levels( THIS% nlevels ), STAT=ios )
    IF( ios > 0 )THEN
     PRINT*,'...allocation error for levels'
     STOP
    ENDIF
    nlevels= THIS% nlevels
    levels = THIS% levels
    coords = THIS% coords

    CALL allocate_ADM()
    CALL allocate_BSSN()

    ! Allocate temporary memory for time integration
    CALL allocate_Ztmp()

    ! Allocate memory for the derivatives of the ADM variables
   ! CALL allocate_GravityAcceleration()

    CALL allocate_grid_function( rad_coord, 'rad_coord', 1 )

    ! Initialize the stress-energy tensor to 0
    DO l= 1, THIS% nlevels, 1
      rad_coord%  levels(l)% var= THIS% rad_coord%  levels(l)% var
      g_phys3_ll% levels(l)% var= THIS% g_phys3_ll% levels(l)% var
      shift_u%    levels(l)% var= THIS% shift_u%    levels(l)% var
      lapse%      levels(l)% var= THIS% lapse%      levels(l)% var
    ENDDO

    CALL set_units('NSM')
    CALL read_options

    CALL allocate_SPH_memory

    ALLOCATE( nlrf_rec (npart)   )
    ALLOCATE( u_rec    (npart)   )
    ALLOCATE( pr_rec   (npart)   )
    ALLOCATE( vel_u_rec(3,npart) )
    ALLOCATE( theta_rec(npart)   )
    ALLOCATE( nstar_rec(npart)   )
    ALLOCATE( s_l_rec  (3,npart) )
    ALLOCATE( e_hat_rec(npart)   )

    IF( debug ) PRINT *, "0"

    ! flag that particles are 'alive'
    IF( .NOT.ALLOCATED( alive ) ) ALLOCATE( alive( npart ) )
    alive( 1:npart )= 1

    CALL allocate_gradient( npart )

    IF( .NOT.ALLOCATED(g4_ll) )THEN
      CALL allocate_metric_on_particles(npart)
    ENDIF

    IF( debug ) PRINT *, "0.25"

    !
    !-- Uncomment the following lines to use he metric from the particles
    !-- This is to compare with he SUBROUTINE recovery_test in TYPE particles,
    !-- and check that both give the same results when using the same data.
    !-- They do on 25.02.2022
    !

    !ALLOCATE( lapse_parts(npart) )
    !ALLOCATE( shift_parts(3,npart) )
    !ALLOCATE( g3_parts   (6,npart) )
    !lapse_parts= parts% get_lapse()
    !shift_parts= parts% get_shift()
    !g3_parts   = parts% get_g3()

    !DO a= 1, npart, 1
    !
    !  CALL compute_g4( lapse_parts(a), &
    !        [shift_parts(1,a), shift_parts(2,a), shift_parts(3,a)], &
    !        [g3_parts(1,a), g3_parts(2,a), g3_parts(3,a), &
    !         g3_parts(4,a), g3_parts(5,a), g3_parts(6,a)], &
    !         g4_ll(1:n_sym4x4,a) )
    !
    !  CALL determinant_sym4x4( g4_ll(1:n_sym4x4,a), det )
    !  IF( ABS(det) < 1D-10 )THEN
    !      PRINT *, "** ERROR! The determinant of the spacetime metric is " &
    !               // "effectively 0 at particle ", a
    !      STOP
    !  ELSEIF( det > 0 )THEN
    !      PRINT *, "** ERROR! The determinant of the spacetime metric is " &
    !               // "positive at particle ", a
    !      STOP
    !  ENDIF
    !
    !ENDDO

    ! Uncomment the following line to use the metric mapped from the mesh
    ! to the particles

    CALL get_metric_on_particles( npart, pos )

    IF( debug ) PRINT *, "0.5"

    ALLOCATE( nlrf_fb (npart) )
    ALLOCATE( u_fb    (npart) )
    ALLOCATE( pr_fb   (npart) )
    ALLOCATE( vel_u_fb(3,npart) )
    ALLOCATE( theta_fb(npart) )
    ALLOCATE( cs_fb   (npart) )

    nstar   = parts% get_nstar_sph()
    nlrf_fb = parts% get_nlrf_sph()
    u_fb    = parts% get_u_sph()
    pr_fb   = parts% get_pressure_cu()
    vel_u_fb= parts% get_vel()
    theta_fb= parts% get_theta()
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
    CALL phys_2_cons( npart, nlrf_fb, u_fb, pr_fb, vel_u_fb, &
                      ! following is output
                      nstar_rec, s_l_rec, e_hat_rec )

    IF( debug ) PRINT *, "2"

    !
    !-- Recover physical fields from conserved fields
    !
    pr_rec= pr_fb

    CALL cons_2_phys( npart, nstar_rec, s_l_rec, e_hat_rec, &
                      ! following is output (pressure is INOUT)
                      nlrf_rec, vel_u_rec, u_rec, pr_rec, theta_rec )

    IF( debug ) PRINT *, "3"

    CALL deallocate_metric_on_particles
    CALL deallocate_gradient
    DEALLOCATE(alive)
    CALL deallocate_sph_memory

    CALL deallocate_grid_function ( rad_coord, 'rad_coord' )
    CALL deallocate_ADM()
    CALL deallocate_Ztmp()
    !CALL deallocate_GravityAcceleration()
    CALL deallocate_BSSN()
    DEALLOCATE( levels )

    IF( debug ) PRINT *, "4"

    matter_objects_loop: DO i_matter= 1, parts% get_n_matter(), 1

      IF( i_matter > 9 )THEN
        WRITE( i_mat, "(I2)" ) i_matter
      ELSE
        WRITE( i_mat, "(I1)" ) i_matter
      ENDIF
      finalnamefile= TRIM(namefile)//"-"//TRIM(i_mat)//".dat"

      !PRINT *, "------------"
      !PRINT *, TRIM(namefile)
      !PRINT *, "------------"
      !PRINT *, finalnamefile
      !PRINT *, "------------"

      ASSOCIATE( npart_in   => parts% get_npart_i(i_matter-1) + 1, &
                 npart_fin  => parts% get_npart_i(i_matter-1) +    &
                               parts% get_npart_i(i_matter) )

      !
      !-- Print the original and recovered fields to formatted file
      !

      INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

      IF( exist )THEN
        OPEN( UNIT= unit_recovery, FILE= TRIM(finalnamefile), &
              STATUS= "REPLACE", FORM= "FORMATTED", &
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

       ! TODO: make the variables export_fomr-xy and export_form_x member
       !       of tpo, not bssn
       ! IF( this% export_form_xy .AND. &
       !     ( pos( 3, a ) >=  0.5D0 .OR. &
       !       pos( 3, a ) <= -0.5D0 ) &
       ! )THEN
       !   CYCLE
       ! ENDIF
       ! IF( this% export_form_x .AND. &
       !     ( pos( 3, a ) >=  0.5D0 .OR. &
       !       pos( 3, a ) <= -0.5D0 .OR. &
       !       pos( 2, a ) >=  0.5D0 .OR. &
       !       pos( 2, a ) <= -0.5D0 ) &
       ! )THEN
       !   CYCLE
       ! ENDIF

        WRITE( UNIT = unit_recovery, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
          a, &                                                        ! 1
          pos      ( jx, a ), &                                       ! 2
          pos      ( jy, a ), &                                       ! 3
          pos      ( jz, a ), &                                       ! 4
          nstar    ( a ),     &                                       ! 5
          nstar_rec( a ),     &                                       ! 6
          nlrf_fb  ( a ),     &                                       ! 7
          nlrf_rec ( a ),     &                                       ! 8
          u_fb     ( a ),     &                                       ! 9
          u_rec    ( a ),     &                                       ! 10
          pr_fb    ( a ),     &                                       ! 11
          pr_rec   ( a ),     &                                       ! 12
          vel_u_fb ( jx, a ), &                                       ! 13
          vel_u_rec( jx, a ), &                                       ! 14
          vel_u_fb ( jy, a ), &                                       ! 15
          vel_u_rec( jy, a ), &                                       ! 16
          vel_u_fb ( jz, a ), &                                       ! 17
          vel_u_rec( jz, a ), &                                       ! 18
          theta_fb ( a ),     &                                       ! 19
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
               i_matter, ", using the metric mapped from the emsh to the ", &
               "particles printed to file ", finalnamefile
      PRINT *

    ENDDO matter_objects_loop

    DEALLOCATE( nlrf_fb )
    DEALLOCATE( u_fb )
    DEALLOCATE( pr_fb )
    DEALLOCATE( vel_u_fb )
    DEALLOCATE( theta_fb )
    DEALLOCATE( cs_fb )

    DEALLOCATE( pos       )
    DEALLOCATE( nstar     )
    DEALLOCATE( nlrf_rec  )
    DEALLOCATE( u_rec     )
    DEALLOCATE( pr_rec    )
    DEALLOCATE( vel_u_rec )
    DEALLOCATE( theta_rec )
    DEALLOCATE( nstar_rec )
    DEALLOCATE( s_l_rec   )
    DEALLOCATE( e_hat_rec )

  END PROCEDURE test_recovery_m2p


END SUBMODULE recovery_m2p
