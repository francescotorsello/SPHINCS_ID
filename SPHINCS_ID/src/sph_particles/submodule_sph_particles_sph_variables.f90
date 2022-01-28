! File:         submodule_sph_particles_sph_variables.f90
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

SUBMODULE (sph_particles) sph_variables

  !****************************************************
  !
  !# THIS SUBMODULE contains the implementation of
  !  the method of TYPE sph_particles
  !  that computes the |sph| variables.
  !
  !  FT 16.10.2020
  !
  !  Renamed from particles_methods to
  !  particles_sph_variables upon improving modularity
  !
  !  FT 12.07.2021
  !
  !****************************************************


  USE constants, ONLY: zero, half, one, two, three


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE compute_and_export_SPH_variables

    !************************************************
    !
    !# Compute the |sph| variables from the
    !  |id|, and print them to a binary file to be
    !  read by |sphincsbssn|, and to a formatted file
    !
    !  FT 18.09.2020
    !
    !************************************************

    USE constants,           ONLY: km2m, m2cm, amu, MSun_geo, &
                                   third, Msun, k_lorene2hydrobase, &
                                   Msun
    USE units,               ONLY: m0c2_cu, set_units
    USE matrix,              ONLY: determinant_4x4_matrix
    USE sph_variables,       ONLY: npart, &  ! particle number
                                   n1,    &  ! particle number for star 1
                                   n2,    &  ! particle number for star 2
                                   pos_u, &  ! particle positions
                                   vel_u, &  ! particle velocities in
                                             ! coordinate frame
                                   nlrf,  &  ! baryon number density in
                                             ! local rest frame
                                   !ehat,  &  ! canonical energy per baryon
                                   nu,    &  ! canonical baryon number per
                                             ! particle
                                   Theta, &  ! Generalized Lorentz factor
                                   h,     &  ! Smoothing length
                                   Pr,    &  ! Pressure
                                   u,     &  ! Internal energy in local rest
                                             ! frame (no kinetic energy)
                                   temp,  &  ! Temperature
                                   av,    &  ! Dissipation
                                   ye,    &  ! Electron fraction
                                   divv,  &  ! Divergence of velocity vel_u
                                   allocate_SPH_memory, &
                                   deallocate_SPH_memory
    USE metric_on_particles, ONLY: allocate_metric_on_particles, &
                                   deallocate_metric_on_particles, &
                                   sq_det_g4
    USE options,             ONLY: basename
    USE input_output,        ONLY: dcount, write_SPHINCS_dump, read_options
    USE NR,                  ONLY: indexx

    !USE RCB_tree_3D,         ONLY: allocate_RCB_tree_memory_3D,&
    !                               deallocate_RCB_tree_memory_3D, iorig
    USE APM,                 ONLY: density_loop
    USE kernel_table,        ONLY: ktable!, dWdv_no_norm, &
                                   !dv_table, dv_table_1, &
                                   !W_no_norm, n_tab_entry
    USE options,             ONLY: ndes
    USE set_h,               ONLY: exact_nei_tree_update
    USE gradient,            ONLY: allocate_gradient, deallocate_gradient
    USE sphincs_sph,         ONLY: density, ncand!, &
                                   !all_clists!, flag_dead_ll_cells
    USE alive_flag,          ONLY: alive
    USE APM,                 ONLY: assign_h
    USE pwp_EOS,             ONLY: select_EOS_parameters, gen_pwp_eos_all, &
                                   get_u_pwp, shorten_eos_name, Gamma_th_1
    USE RCB_tree_3D,         ONLY: iorig, nic, nfinal, nprev, lpart, &
                                   rpart, allocate_RCB_tree_memory_3D, &
                                   deallocate_RCB_tree_memory_3D
    USE matrix,              ONLY: invert_3x3_matrix
    USE analyze,             ONLY: COM
    USE tensor,              ONLY: n_sym4x4, &
                                   itt, itx, ity, itz, &
                                   ixx, ixy, ixz, iyy, iyz, izz

    IMPLICIT NONE

    ! The flag call_flag is set different than 0 if the SUBROUTINE
    ! compute_and_export_SPH_variables is called
    INTEGER, SAVE:: call_flag= 0

    !INTEGER, PARAMETER:: max_it_h= 1

    ! Spacetime indices \mu and \nu
    INTEGER:: nus, mus, cnt1, a, i_matter!, itr2, inde, index1!, cnt2
    INTEGER:: n_problematic_h
    INTEGER:: itot, l, ill!, b, k

    DOUBLE PRECISION:: g4(0:3,0:3)
    DOUBLE PRECISION:: gg4(THIS% npart,n_sym4x4)
    DOUBLE PRECISION:: sq_detg4(THIS% npart)
    DOUBLE PRECISION:: det, sq_g, Theta_a!, &!nu_max1, nu_max2, &
                       !nu_tmp, nu_thres1, nu_thres2
    DOUBLE PRECISION:: com_x_newt, com_y_newt, com_z_newt, com_d_newt
    DOUBLE PRECISION:: com_x_1pn, com_y_1pn, com_z_1pn, com_d_1pn
    DOUBLE PRECISION:: px, py, pz

    !DOUBLE PRECISION:: ha, ha_1, ha_3, va, mat(3,3), mat_1(3,3), xa, ya, za
    !DOUBLE PRECISION:: mat_xx, mat_xy, mat_xz, mat_yy
    !DOUBLE PRECISION:: mat_yz, mat_zz, Wdx, Wdy, Wdz, dx, dy, dz, Wab, &
    !                   Wab_ha, Wi, Wi1, dvv

    LOGICAL:: few_ncand!, invertible_matrix

    LOGICAL, PARAMETER:: debug= .FALSE.

    CHARACTER( LEN= : ), ALLOCATABLE:: compose_namefile
    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    TYPE(timer):: find_h_bruteforce_timer

    find_h_bruteforce_timer= timer( "find_h_bruteforce_timer" )

    PRINT *, "** Executing the compute_and_export_SPH_variables " &
             // "subroutine..."
    PRINT *

    !
    !-- Set up the MODULE variables in MODULE sph_variables
    !-- (used by write_SPHINCS_dump)
    !
    npart= THIS% npart
    n1= THIS% npart_i(1)
    IF( THIS% n_matter == 2 ) n2= THIS% npart_i(2)

    CALL set_units('NSM')
    CALL read_options

    CALL allocate_SPH_memory
    CALL allocate_metric_on_particles( THIS% npart )

    IF( debug ) PRINT *, "1"

    IF(.NOT.ALLOCATED( THIS% nu ))THEN
      ALLOCATE( THIS% nu( THIS% npart ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array nu ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array nu" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% nlrf ))THEN
      ALLOCATE( THIS% nlrf( THIS% npart ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array nlrf ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array nlrf" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% Theta ))THEN
      ALLOCATE( THIS% Theta( THIS% npart ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array Theta ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array Theta" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% v ))THEN
      ALLOCATE( THIS% v( 0:3, THIS% npart ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array v ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array v" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% h ))THEN
      ALLOCATE( THIS% h( THIS% npart ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array h ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array h" )
    ENDIF

    IF( debug ) PRINT *, "2"

    !
    !-- Compute SPH quantities
    !

    CALL THIS% sph_computer_timer% start_timer()
!    !$OMP PARALLEL DO DEFAULT( NONE ) &
!    !$OMP             SHARED( THIS, pos_u, vel_u, sq_det_g4, Theta, nlrf, &
!    !$OMP                     Pr, u, temp, av, divv ) &
!    !$OMP             PRIVATE( itr, g4, gg4, det, sq_g, Theta_a, &
!    !$OMP                      nus, mus )
    compute_SPH_variables_on_particles: DO itr= 1, THIS% npart, 1

      ! Particle positions [Msun_geo]
      pos_u(1,itr)= THIS% pos(1,itr)
      pos_u(2,itr)= THIS% pos(2,itr)
      pos_u(3,itr)= THIS% pos(3,itr)

      ! Coordinate velocity of the fluid [c]
      THIS% v(0,itr)= one
      THIS% v(1,itr)= THIS% lapse_parts(itr)*THIS% v_euler_parts_x(itr) &
                    - THIS% shift_parts_x(itr)
      THIS% v(2,itr)= THIS% lapse_parts(itr)*THIS% v_euler_parts_y(itr) &
                    - THIS% shift_parts_y(itr)
      THIS% v(3,itr)= THIS% lapse_parts(itr)*THIS% v_euler_parts_z(itr) &
                    - THIS% shift_parts_z(itr)
      vel_u(1,itr)  = THIS% v(1,itr)
      vel_u(2,itr)  = THIS% v(2,itr)
      vel_u(3,itr)  = THIS% v(3,itr)

      !
      !-- Metric as matrix for easy manipulation
      !
      g4(0,0)= - THIS% lapse_parts(itr)**2 &
             + THIS% g_xx_parts(itr)*THIS% shift_parts_x(itr) &
              *THIS% shift_parts_x(itr)&
             + 2 * THIS% g_xy_parts(itr)*THIS% shift_parts_x(itr) &
              *THIS% shift_parts_y(itr)&
             + 2 * THIS% g_xz_parts(itr)*THIS% shift_parts_x(itr) &
              *THIS% shift_parts_z(itr)&
             + THIS% g_yy_parts(itr)*THIS% shift_parts_y(itr) &
              *THIS% shift_parts_y(itr)&
             + 2 * THIS% g_yz_parts(itr)*THIS% shift_parts_y(itr) &
              *THIS% shift_parts_z(itr)&
             + THIS% g_zz_parts(itr)*THIS% shift_parts_z(itr) &
              *THIS% shift_parts_z(itr)
      g4(0,1)= THIS% g_xx_parts(itr)*THIS% shift_parts_x(itr) &
             + THIS% g_xy_parts(itr)*THIS% shift_parts_y(itr) &
             + THIS% g_xz_parts(itr)*THIS% shift_parts_z(itr)
      g4(0,2)= THIS% g_xy_parts(itr)*THIS% shift_parts_x(itr) &
             + THIS% g_yy_parts(itr)*THIS% shift_parts_y(itr) &
             + THIS% g_yz_parts(itr)*THIS% shift_parts_z(itr)
      g4(0,3)= THIS% g_xz_parts(itr)*THIS% shift_parts_x(itr) &
             + THIS% g_yz_parts(itr)*THIS% shift_parts_y(itr) &
             + THIS% g_zz_parts(itr)*THIS% shift_parts_z(itr)

      g4(1,0)= THIS% g_xx_parts(itr)*THIS% shift_parts_x(itr) &
             + THIS% g_xy_parts(itr)*THIS% shift_parts_y(itr) &
             + THIS% g_xz_parts(itr)*THIS% shift_parts_z(itr)
      g4(1,1)= THIS% g_xx_parts(itr)
      g4(1,2)= THIS% g_xy_parts(itr)
      g4(1,3)= THIS% g_xz_parts(itr)

      g4(2,0)= THIS% g_xy_parts(itr)*THIS% shift_parts_x(itr) &
             + THIS% g_yy_parts(itr)*THIS% shift_parts_y(itr) &
             + THIS% g_yz_parts(itr)*THIS% shift_parts_z(itr)
      g4(2,1)= THIS% g_xy_parts(itr)
      g4(2,2)= THIS% g_yy_parts(itr)
      g4(2,3)= THIS% g_yz_parts(itr)

      g4(3,0)= THIS% g_xz_parts(itr)*THIS% shift_parts_x(itr) &
             + THIS% g_yz_parts(itr)*THIS% shift_parts_y(itr) &
             + THIS% g_zz_parts(itr)*THIS% shift_parts_z(itr)
      g4(3,1)= THIS% g_xz_parts(itr)
      g4(3,2)= THIS% g_yz_parts(itr)
      g4(3,3)= THIS% g_zz_parts(itr)

      gg4(itr,1)= g4(0,0)
      gg4(itr,2)= g4(0,1)
      gg4(itr,3)= g4(0,2)
      gg4(itr,4)= g4(0,3)
      gg4(itr,5)= g4(1,1)
      gg4(itr,6)= g4(1,2)
      gg4(itr,7)= g4(1,3)
      gg4(itr,8)= g4(2,2)
      gg4(itr,9)= g4(2,3)
      gg4(itr,10)= g4(3,3)

      ! sqrt(-det(g4))
      CALL determinant_4x4_matrix(g4,det)
      IF( ABS(det) < 1D-10 )THEN
          PRINT *, "The determinant of the spacetime metric is " &
                   // "effectively 0 at particle ", itr
          STOP
      ELSEIF( det > 0 )THEN
          PRINT *, "The determinant of the spacetime metric is " &
                   // "positive at particle ", itr
          STOP
      ENDIF
      sq_g= SQRT(-det)
      sq_det_g4(itr)= sq_g
      sq_detg4(itr)= sq_g

      !
      !-- Generalized Lorentz factor
      !
      Theta_a= zero
      DO nus=0,3
        DO mus=0,3
          Theta_a= Theta_a &
                   + g4(mus,nus)*THIS% v(mus,itr)*THIS% v(nus,itr)
        ENDDO
      ENDDO
      Theta_a         = one/SQRT(-Theta_a)
      Theta(itr)      = Theta_a
      THIS% Theta(itr)= Theta_a

      ! Baryon density in the local rest frame [baryon (Msun_geo)^{-3}]
      ! Computed from the LORENE baryon mass density in [kg/m^3]
      nlrf(itr)= THIS% baryon_density_parts(itr)
      THIS% nlrf(itr)= THIS% baryon_density_parts(itr)

      ! Specific internal energy [c^2]
      u(itr)= THIS% specific_energy_parts(itr)

      ! Pressure [amu*c**2/(Msun_geo**3)]
      !          dimensions: [(M/L**3)*L**2/T**2]= [M/(L*T**2)], same as
      !                      energy density
      Pr(itr)= THIS% pressure_parts(itr)
      THIS% pressure_parts_cu(itr)= Pr(itr)

    ENDDO compute_SPH_variables_on_particles
!    !$OMP END PARALLEL DO

    ! Temperature: here dummy
    temp= zero

    ! Dissipation parameter
    av= one

    ! Velocity divergence
    divv= zero

    IF( debug ) PRINT *, "3"


    ! Compute nstar (proper baryon number density) from the ID
    THIS% nstar= ( THIS% nlrf*THIS% Theta )*sq_det_g4

    !
    !-- Compute the particle proper mass, if not computed yet
    !
    IF( .NOT.ALLOCATED( THIS% pvol ) )THEN
      PRINT *, "** ERROR! The array pvol is not allocated. ", &
               "Stopping..."
      PRINT *
      STOP
    ENDIF
    IF( .NOT.( THIS% distribution_id == id_particles_on_spherical_surfaces &
        .OR. &
        ( THIS% distribution_id == id_particles_from_file &
          .AND. THIS% read_nu ) ) )THEN

      THIS% pmass= THIS% nstar*THIS% pvol

    ENDIF

    ! Compute particle number density from the ID
    THIS% particle_density= ( THIS% nstar )/( THIS% pmass )

    IF( debug ) PRINT *, "4"

    !
    !-- Compute the first guess for the smoothing length, if the APM was not
    !-- used
    !
    DO i_matter= 1, THIS% n_matter, 1

      IF( .NOT.THIS% apm_iterate(i_matter) )THEN

        IF( debug ) PRINT *, "Compute first guess for the smoothing length ", &
                             "h, for particles on matter object", itr,"..."

        compute_h_guess: &
        DO itr= THIS% npart_i(i_matter-1) + 1, &
                THIS% npart_i(i_matter-1) + THIS% npart_i(i_matter), 1

          THIS% h(itr)= three*(THIS% pvol(itr))**third
          h(itr)= THIS% h(itr)
          ! /(Msun_geo**3)
          IF( debug .AND. THIS% h(itr) <= zero )THEN
            PRINT *, "** ERROR! h(", itr, ")=", THIS% h(itr)
            PRINT *, "Stopping..."
            PRINT *
            STOP
          ENDIF

        ENDDO compute_h_guess

      ENDIF

    ENDDO

    IF( debug ) PRINT *, "1"

    !-------------------------------------!
    !--  Assignment of nu on the stars. --!
    !-------------------------------------!

!    IF( THIS% redistribute_nu )THEN
!
!      !---------------------------------------------------------------------!
!      !--  Assignment of nu on the stars, with the purpose                --!
!      !--  of having a more uniform nu over the particles without losing  --!
!      !--  baryon mass. This is used only on the lattice, optionally.     --!
!      !---------------------------------------------------------------------!
!
!      IF( THIS% distribution_id == id_particles_on_spherical_surfaces )THEN
!        PRINT *, "** ERROR! Particle placer ", THIS% distribution_id, &
!                 " is not compatible with redistribute_nu= .TRUE."
!        PRINT *, " * Check the parameter file lorene_bns_id_particles.par. ", &
!                 "Stopping..."
!        PRINT *
!        STOP
!      ENDIF
!
!      nu_max1= nlrf( THIS% baryon_density_index( THIS% npart1 ) )&
!              *THIS% pvol( THIS% npart1 ) &
!              *Theta( THIS% baryon_density_index( THIS% npart1 ) )&
!              *sq_det_g4( THIS% baryon_density_index( THIS% npart1 ) )
!      nu_max2= nlrf( THIS% baryon_density_index( THIS% npart ) )&
!              *THIS% pvol( THIS% npart ) &
!              *Theta( THIS% baryon_density_index( THIS% npart ) )&
!              *sq_det_g4( THIS% baryon_density_index( THIS% npart ) )
!
!      nu_thres1= nu_max1/THIS% nu_ratio
!      nu_thres2= nu_max2/THIS% nu_ratio
!
!      ! Reset the total baryon number to 0 (necessary), and nu to an arbitrary
!      ! value (to make debugging easier)
!
!      nu= one
!      THIS% nu= one
!      THIS% nbar_tot= zero
!      THIS% nbar1= zero
!      THIS% nbar2= zero
!
!      cnt1= 0
!      compute_nu_on_particles_star1: DO itr= THIS% npart1, 1, -1
!
!        cnt1= cnt1 + 1
!
!        nu_tmp= nlrf( THIS% baryon_density_index( itr ) ) &
!                *THIS% pvol(itr) &
!                *Theta( THIS% baryon_density_index( itr ) )&
!                *sq_det_g4( THIS% baryon_density_index( itr ) )
!
!        !IF( itr == THIS% npart1 ) nu_max= nu_tmp ! move this out of the loop
!
!        IF( nu_tmp > nu_thres1 )THEN
!          nu( THIS% baryon_density_index( itr ) )      = nu_tmp
!          THIS% nu( THIS% baryon_density_index( itr ) )= nu_tmp
!        ELSE
!          nu( THIS% baryon_density_index( itr ) )      = nu_thres1
!          THIS% nu( THIS% baryon_density_index( itr ) )= nu_thres1
!        ENDIF
!
!        THIS% nbar1= THIS% nbar1 + &
!                     THIS% nu( THIS% baryon_density_index( itr ) )
!
!        IF( THIS% nbar1*amu/MSun > THIS% masses(1) )THEN
!          EXIT
!        ENDIF
!
!      ENDDO compute_nu_on_particles_star1
!
!      cnt2= 0
!      compute_nu_on_particles_star2: DO itr= THIS% npart, THIS% npart1 + 1, -1
!
!        cnt2= cnt2 + 1
!
!        nu_tmp= nlrf( THIS% baryon_density_index( itr ) ) &
!                *THIS% pvol(itr) &
!                *Theta( THIS% baryon_density_index( itr ) ) &
!                *sq_det_g4( THIS% baryon_density_index( itr ) )
!
!        !IF( itr == THIS% npart ) nu_max= nu_tmp
!
!        IF( nu_tmp > nu_thres2 )THEN
!          nu( THIS% baryon_density_index( itr ) )      = nu_tmp
!          THIS% nu( THIS% baryon_density_index( itr ) )= nu_tmp
!        ELSE
!          nu( THIS% baryon_density_index( itr ) )      = nu_thres2
!          THIS% nu( THIS% baryon_density_index( itr ) )= nu_thres2
!        ENDIF
!
!        THIS% nbar2= THIS% nbar2 + &
!                     THIS% nu( THIS% baryon_density_index( itr ) )
!
!        IF( THIS% nbar2*amu/MSun > THIS% masses(2) )THEN
!          EXIT
!        ENDIF
!
!      ENDDO compute_nu_on_particles_star2
!      THIS% nbar_tot= THIS% nbar1 + THIS% nbar2
!
!      !
!      !-- Reshape MODULE variables
!      !
!
!      CALL THIS% reshape_sph_field( pos_u, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( vel_u, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( Theta, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( h, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( nlrf, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( u, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( Pr, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( nu, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( temp, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( av, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( divv, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      !
!      !-- Reshape TYPE member SPH variables
!      !
!
!      CALL THIS% reshape_sph_field( THIS% pos, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% v, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% v_euler_parts_x, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% v_euler_parts_y, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% v_euler_parts_z, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% Theta, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% h, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% baryon_density_parts, cnt1, &
!                                    cnt2, THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% nlrf, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% energy_density_parts, cnt1, &
!                                    cnt2, THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% specific_energy_parts, cnt1, &
!                                    cnt2, THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% pressure_parts, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% pressure_parts_cu, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% nu, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% pvol, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      !
!      !-- Reshape TYPE member spacetime variables
!      !
!
!      CALL THIS% reshape_sph_field( THIS% lapse_parts, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% shift_parts_x, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% shift_parts_y, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% shift_parts_z, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% g_xx_parts, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% g_xy_parts, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% g_xz_parts, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% g_yy_parts, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% g_yz_parts, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      CALL THIS% reshape_sph_field( THIS% g_zz_parts, cnt1, cnt2, &
!                                    THIS% baryon_density_index )
!
!      !
!      !-- Reassign particle numbers
!      !
!
!      npart= cnt1 + cnt2
!      THIS% npart= npart
!      THIS% npart1= cnt1
!      THIS% npart2= cnt2
!      n1= THIS% npart1
!      n2= THIS% npart2
!
!      PRINT *, " * Particles replaced after reassigning nu."
!      PRINT *, " * New number of particles=", THIS% npart
!      PRINT *
!      PRINT *, " * Number of particles on NS 1=", THIS% npart1
!      PRINT *, " * Number of particles on NS 2=", THIS% npart2
!      PRINT *

    !----------------------------------------------!
    !--  Assignment of nu on the matter objects. --!
    !----------------------------------------------!

    ! TODO: The code would be much cleaner if nu and h were assigned in the
    !       constructor for all the particle placers (now only the lattices
    !       and the particles ad from file without nu don't assign them)
    !       If that was the case, then this SUBROUTINE could start out by
    !       estimating the SPH density, and compute everything from there.

    DO i_matter= 1, THIS% n_matter, 1

      ASSOCIATE( npart_in   => THIS% npart_i(i_matter-1) + 1, &
                 npart_fin  => THIS% npart_i(i_matter-1) +    &
                               THIS% npart_i(i_matter) )

        IF( THIS% apm_iterate(i_matter) )THEN

          ! If the APM was used...

          ! Do nothing, nu is already computed and reflected in the constructor
          nu( npart_in : npart_fin )= THIS% nu( npart_in : npart_fin )
          THIS% nbar_i(i_matter)= SUM( THIS% nu( npart_in : npart_fin ), DIM=1 )

        ELSEIF( THIS% distribution_id == id_particles_from_file &
                .AND. THIS% read_nu )THEN

          ! If the particle positions and nu were read from formatted file...

          ! Do nothing, nu is already read from file
          nu( npart_in : npart_fin )= THIS% nu( npart_in : npart_fin )
          THIS% nbar_i(i_matter)= SUM( THIS% nu( npart_in : npart_fin ), DIM=1 )

        ELSEIF( THIS% distribution_id== id_particles_on_spherical_surfaces )THEN
        !ELSE

          ! If the APM was not used...

          ! Set nu based on the particle mass...

          DO itr= npart_in, npart_fin, 1
            nu(itr)= THIS% pmass(itr)*MSun/amu
            THIS% nu(itr)= nu(itr)
            THIS% nbar_i(i_matter)= THIS% nbar_i(i_matter) + nu(itr)
          ENDDO

        ELSE

          ! If the APM was not used and the particles are on lattices...

          DO itr= npart_in, npart_fin, 1
            nu(itr)= nlrf(itr)*THIS% pvol(itr)*Theta( itr )*sq_det_g4( itr )
            THIS% nu(itr)= nu(itr)
            THIS% nbar_i(i_matter)= THIS% nbar_i(i_matter) + nu(itr)
          ENDDO

        ENDIF
        THIS% nbar_tot= THIS% nbar_tot + THIS% nbar_i(i_matter)

        equal_mass_binary: &
        IF( i_matter == 1 .AND. THIS% n_matter == 2 )THEN

          IF( ABS(THIS% mass_ratios(1) - THIS% mass_ratios(2)) &
            /THIS% mass_ratios(2) <= 0.005 .AND. THIS% reflect_particles_x &
            .AND. &
            THIS% distribution_id /= id_particles_from_file )THEN

            ! Consistency check
            IF( THIS% npart_i(i_matter) +      &
            THIS% npart_i(i_matter+1) /= THIS% npart )THEN
              PRINT *, "** ERROR! npart_next /= THIS% npart! "
              PRINT *, "   npart_next=", THIS% npart_i(i_matter) +      &
              THIS% npart_i(i_matter+1)
              PRINT *, "   THIS% npart=", THIS% npart
              PRINT *, "   Stopping..."
              PRINT *
              STOP
            ENDIF

            nu( npart_fin + 1:THIS% npart_i(i_matter) +      &
            THIS% npart_i(i_matter+1) )= nu( npart_in : npart_fin )
            THIS% nu( npart_fin + 1:THIS% npart_i(i_matter) +      &
            THIS% npart_i(i_matter+1) )= nu( npart_in : npart_fin )
            THIS% nbar_i(i_matter + 1)= THIS% nbar_i(i_matter)

            THIS% nbar_tot= THIS% nbar_tot + THIS% nbar_i(i_matter + 1)

            ! Consistency check
            IF( THIS% nbar_tot /= 2*THIS% nbar_i(i_matter + 1) )THEN
              PRINT *, "** ERROR! THIS% nbar_tot /= 2*THIS% nbar(i_matter + 1) ",&
                       "   when reflecting particles or a binary system"
              PRINT *, "   THIS% nbar_tot=", THIS% nbar_tot
              PRINT *, "   2*THIS% nbar(", i_matter + 1, ")=", &
                           2*THIS% nbar_i(i_matter + 1)
              PRINT *, "   Stopping..."
              PRINT *
              STOP
            ENDIF

            EXIT

          ENDIF

        ENDIF equal_mass_binary

      END ASSOCIATE

    ENDDO

    !------------------------------------------------------------------------!
    ! Compute SPH density estimate (nu has to be assigned before this step)  !
    !------------------------------------------------------------------------!

    CALL allocate_RCB_tree_memory_3D(npart)
    iorig(1:npart)= (/ (a,a=1,npart) /)

    CALL allocate_gradient( npart )

    IF( debug ) PRINT *, "-1"

    PRINT *, " * Assigning h..."
    PRINT *

    ! Determine smoothing length so that each particle has exactly
    ! 300 neighbours inside 2h
    CALL assign_h( ndes, &
                   THIS% npart, &
                   THIS% pos, THIS% h, & ! Input
                   h )             ! Output

    IF( debug ) PRINT *, "101.5"

    CALL find_h_bruteforce_timer% start_timer()
    n_problematic_h= 0
    check_h: DO a= 1, THIS% npart, 1

      IF( ISNAN( h(a) ) .OR. h(a) <= zero )THEN

        n_problematic_h= n_problematic_h + 1

        h(a)= find_h_backup( a, THIS% npart, THIS% pos, ndes )

        IF( ISNAN( h(a) ) .OR. h(a) <= zero )THEN
          PRINT *, "** ERROR! h=0 on particle ", a, "even with the brute", &
                   " force method."
          PRINT *, "   Particle position: ", THIS% pos(:,a)
          STOP
        ENDIF

      ENDIF

    ENDDO check_h
    CALL find_h_bruteforce_timer% stop_timer()
    CALL find_h_bruteforce_timer% print_timer( 2 )

    PRINT *, " * The smoothing length was found brute-force for ", &
             n_problematic_h, " particles."
    PRINT *

    PRINT *, " * Computing neighbours' tree..."
    PRINT *

    IF( SUM(THIS% nu, DIM= 1)/SIZE(THIS% nu) == zero )THEN
      PRINT *, "** ERROR! Average nu= 0. Are you assigning values to the ", &
               "TYPE member array?"
      PRINT *, "Stopping..."
      STOP
    ENDIF

    cnt1= 0
    DO

      few_ncand= .FALSE.

      ! Redo the previous step slightly different (it's built-in;
      ! exact_nei_tree_update does not work if I don't call assign_h first),
      ! then update the neighbour-tree and fill the neighbour-data
      CALL exact_nei_tree_update( ndes,        &
                                  THIS% npart, &
                                  THIS% pos,   &
                                  THIS% nu )

      EXIT

      ll_cell_loop2: DO ill= 1, nfinal

        itot= nprev + ill
        IF( nic(itot) == 0 ) CYCLE

        IF( ncand(ill) < ndes - 1 )THEN

          ! Increase the smoothing lengths of the paricles inside the cell,
          ! and rebuild the tree
          few_ncand= .TRUE.

          particle_in_cell_loop: DO l= lpart(itot), rpart(itot)

            h(l)= three*h(l)

          ENDDO particle_in_cell_loop

        ELSE

          few_ncand= .FALSE.

        ENDIF

      ENDDO ll_cell_loop2

      cnt1= cnt1 + 1

      IF( .NOT.few_ncand .OR. cnt1 >= max_it_tree )THEN
        PRINT *, " * Smoothing lengths assigned and neighbours' tree is built."
        EXIT
      ENDIF

    ENDDO

    CALL find_h_bruteforce_timer% start_timer()
    n_problematic_h= 0
    check_h2: DO a= 1, THIS% npart, 1

      IF( ISNAN( h(a) ) .OR. h(a) <= zero )THEN

        n_problematic_h= n_problematic_h + 1
        h(a)= find_h_backup( a, THIS% npart, THIS% pos, ndes )
        PRINT *, h(a)
        IF( ISNAN( h(a) ) .OR. h(a) <= zero )THEN
          PRINT *, "** ERROR! h=0 on particle ", a, "even with the brute", &
                   " force method."
          PRINT *, "   Particle position: ", THIS% pos(:,a)
          STOP
        ENDIF

      ENDIF

    ENDDO check_h2
    CALL find_h_bruteforce_timer% stop_timer()
    CALL find_h_bruteforce_timer% print_timer( 2 )

    PRINT *, " * The smoothing length was found brute-force for ", &
             n_problematic_h, " particles."
    PRINT *

    ! Update the member variables storing smoothing length and particle volume
    THIS% h= h
    THIS% pvol= ( THIS% h/three )**three

 !   !PRINT *
 !   !PRINT *, "nfinal= ", nfinal
 !   ll_cell_loop: DO ill= 1, nfinal
 !
 !     itot= nprev + ill
 !     IF( nic(itot) == 0 ) CYCLE
 !
 !     particle_loop: DO l= lpart(itot), rpart(itot)
 !
 !       a=         iorig(l)
 !
 !       ha=        h(a)
 !       ha_1=      one/ha
 !       ha_3=      ha_1*ha_1*ha_1
 !
 !       xa=        pos_u(1,a)
 !       ya=        pos_u(2,a)
 !       za=        pos_u(3,a)
 !
 !       ! initialize correction matrix
 !       mat_xx=    zero
 !       mat_xy=    zero
 !       mat_xz=    zero
 !       mat_yy=    zero
 !       mat_yz=    zero
 !       mat_zz=    zero
 !
 !       cnt1= 0
 !       cnt2= 0
 !       cand_loop: DO k= 1, ncand(ill)
 !
 !         b=      all_clists(ill)%list(k)
 !
 !         IF( b == a )THEN
 !           cnt1= cnt1 + 1
 !         ENDIF
 !         IF( xa == pos_u(1,b) .AND. ya == pos_u(2,b) .AND. za == pos_u(3,b) &
 !         )THEN
 !           cnt2= cnt2 + 1
 !         ENDIF
 !
 !         ! Distances (ATTENTION: flatspace version !!!)
 !         dx=     xa - pos_u(1,b)
 !         dy=     ya - pos_u(2,b)
 !         dz=     za - pos_u(3,b)
 !         va=     SQRT(dx*dx + dy*dy + dz*dz)*ha_1
 !
 !         !IF( dx == 0 .AND. dy == 0 .AND. dz == 0 )THEN
 !         !  PRINT *, "va=", va
 !         !  PRINT *, "dz=", dx
 !         !  PRINT *, "dy=", dy
 !         !  PRINT *, "dz=", dz
 !         !  PRINT *, "xa=", xa
 !         !  PRINT *, "ya=", ya
 !         !  PRINT *, "za=", za
 !         !  PRINT *, "pos_u(1,b)", pos_u(1,b)
 !         !  PRINT *, "pos_u(2,b)", pos_u(2,b)
 !         !  PRINT *, "pos_u(3,b)", pos_u(3,b)
 !         !  STOP
 !         !ENDIF
 !
 !         ! get interpolation indices
 !         inde=  MIN(INT(va*dv_table_1),n_tab_entry)
 !         index1= MIN(inde + 1,n_tab_entry)
 !
 !         ! get tabulated values
 !         Wi=     W_no_norm(inde)
 !         Wi1=    W_no_norm(index1)
 !
 !         ! interpolate
 !         dvv=    (va - DBLE(inde)*dv_table)*dv_table_1
 !         Wab_ha= Wi + (Wi1 - Wi)*dvv
 !
 !         ! "correction matrix" for derivs
 !         Wdx=    Wab_ha*dx
 !         Wdy=    Wab_ha*dy
 !         Wdz=    Wab_ha*dz
 !         mat_xx= mat_xx + Wdx*dx
 !         mat_xy= mat_xy + Wdx*dy
 !         mat_xz= mat_xz + Wdx*dz
 !         mat_yy= mat_yy + Wdy*dy
 !         mat_yz= mat_yz + Wdy*dz
 !         mat_zz= mat_zz + Wdz*dz
 !
 !       ENDDO cand_loop
 !
 !       ! correction matrix
 !       mat(1,1)= mat_xx
 !       mat(2,1)= mat_xy
 !       mat(3,1)= mat_xz
 !
 !       mat(1,2)= mat_xy
 !       mat(2,2)= mat_yy
 !       mat(3,2)= mat_yz
 !
 !       mat(1,3)= mat_xz
 !       mat(2,3)= mat_yz
 !       mat(3,3)= mat_zz
 !
 !       ! invert it
 !       CALL invert_3x3_matrix(mat,mat_1,invertible_matrix)
 !
 !      ! IF( .NOT.invertible_matrix )THEN
 !      !   PRINT *, "a= ", a
 !      !   PRINT *, "h(a)= ", h(a)
 !      !   PRINT *, "pos_u= ", pos_u(1,b), pos_u(2,b), pos_u(3,b)
 !      !   PRINT *, "nprev= ", nprev
 !      !   PRINT *, "ill= ", ill
 !      !   PRINT *, "itot= ", itot
 !      !   PRINT *, "ncand(ill)= ", ncand(ill)
 !      !   PRINT *, "cnt1= ", cnt1
 !      !   PRINT *, "cnt2= ", cnt2
 !      !   PRINT *
 !      !   STOP
 !      ! ENDIF
 !
 !     ENDDO particle_loop
 !
 !   ENDDO ll_cell_loop

    !
    !-- Compute the proper baryon number density with kernel interpolation
    !

    PRINT *, " * Computing SPH proper baryon number density with kernel", &
             " interpolation..."
    PRINT *
    ! density calls dens_ll_cell, which computes nstar on particle a as
    ! Na=     Na + nu(b)*Wab_ha, so this is nstar= nlrf*sq_g*Theta
    ! It has to be compared with nstar= nlrf*sq_g*Theta
    CALL density( THIS% npart, &
                  THIS% pos,   &
                  THIS% nstar_int )

    !-------------------------------------------------------------------------!
    !-------------------------------------------------------------------------!
    ! This point here is CRUCIAL: the particle distribution may NOT resolve   !
    ! properly the steep density gradient at the surface, even if the APM     !
    ! is used. This implies that the kernel interpolated nstar_int will be    !
    ! different than nstar close to the surface.                              !
    ! The error can be such that the recovery fails in SPHINCS_BSSN, and      !
    ! this is of course a problem. Now, the density used in SPHINCS_BSSN      !
    ! during the evolution is not the one given in the ID. Hence, once        !
    ! nstar_int is computed, nlrf should be recomputed from it, so that the   !
    ! density on the particles corresponds to the density that "they see",    !
    ! that is, the kernel interpolated density that uses these values of nu.  !
    !-------------------------------------------------------------------------!
    !-------------------------------------------------------------------------!

    THIS% nlrf_int= ( THIS% nstar_int/THIS% Theta )/sq_det_g4
    nlrf= THIS% nlrf_int

    !-----------------------------------------------------------------------!
    ! For single and piecewise polytropes, do not use the pressure          !
    ! and specific internal energy from the ID.                             !
    ! Compute them using the exact formulas for piecewise                   !
    ! polytropic instead, starting from the kernel interpolated density     !
    !-----------------------------------------------------------------------!

    matter_objects_loop: DO i_matter= 1, THIS% n_matter, 1

      ASSOCIATE( npart_in  => THIS% npart_i(i_matter-1) + 1, &
                 npart_fin => THIS% npart_i(i_matter-1) +    &
                              THIS% npart_i(i_matter) )

      IF( THIS% all_eos(i_matter)% eos_parameters(1) == DBLE(1) )THEN
      ! If the |eos| is polytropic

        PRINT *, " * Computing pressure and specific internal energy from", &
                 " the baryon mass density, using the exact formulas for", &
                 " single polytropic EOS..."
        PRINT *

        ! Formulas from Read et al. (2009), https://arxiv.org/abs/0812.2163

        IF( THIS% cold_system )THEN
          ! If the system is cold, compute pressure and specific energy
          ! exactly using the polytropic EOS
          Pr(npart_in:npart_fin)= THIS% all_eos(i_matter)% eos_parameters(3) &
                               *( THIS% nlrf_int(npart_in:npart_fin)*m0c2_cu ) &
                                **THIS% all_eos(i_matter)% eos_parameters(2)

          u(npart_in:npart_fin)= ( Pr(npart_in:npart_fin) &
                    /(THIS% nlrf_int(npart_in:npart_fin)*m0c2_cu &
                    *( THIS% all_eos(i_matter)% eos_parameters(2) - one ) ) )

          Pr(npart_in:npart_fin)= Pr(npart_in:npart_fin)/m0c2_cu
          THIS% pressure_parts_cu(npart_in:npart_fin)= Pr(npart_in:npart_fin)
          THIS% u_pwp(npart_in:npart_fin)= u(npart_in:npart_fin)
        ELSE
          ! If the system is hot, that is, has a thermal component, then
          ! the density and the specific energy (the latter including both
          ! cold and thermal part) should be supplied in the ID.
          ! The pressure is computed using them (see pwp_EOS MODULE).
          u(npart_in:npart_fin)= THIS% specific_energy_parts(npart_in:npart_fin)

          DO a= npart_in, npart_fin, 1

            Pr(a)= &
            ! cold pressure
            THIS% all_eos(i_matter)% eos_parameters(3) &
              *( THIS% nlrf_int(a)*m0c2_cu ) &
              **THIS% all_eos(i_matter)% eos_parameters(2) &
            + &
            ! thermal pressure
            Gamma_th_1*( THIS% nlrf_int(a)*m0c2_cu )* &
              MAX(u(a) - ( Pr(a)/(THIS% nlrf_int(a)*m0c2_cu &
                *( THIS% all_eos(i_matter)% eos_parameters(2) - one ) ) ), zero)

          ENDDO
          THIS% pressure_parts_cu(npart_in:npart_fin)= Pr(npart_in:npart_fin)
          THIS% u_pwp(npart_in:npart_fin)= u(npart_in:npart_fin)

        ENDIF

      ELSEIF( THIS% all_eos(i_matter)% eos_parameters(1) == DBLE(110) )THEN
      ! If the |eos| is piecewise polytropic

        PRINT *, " * Computing pressure and specific internal energy from", &
                 " the baryon mass density, using the exact formulas for", &
                 " piecewise polytropic EOS..."
        PRINT *

        CALL select_EOS_parameters( &
                shorten_eos_name(THIS% all_eos(i_matter)% eos_name) )

        CALL gen_pwp_eos_all( THIS% npart_i(i_matter), &
                              THIS% nlrf_int(npart_in:npart_fin)*m0c2_cu, &
                              u(npart_in:npart_fin) )

        THIS% pressure_parts_cu(npart_in:npart_fin)= Pr(npart_in:npart_fin)

        IF( THIS% cold_system )THEN
          ! If the system is cold, get the specific energy computed
          ! exactly using the piecewise polytropic EOS
          THIS% u_pwp(npart_in:npart_fin)= get_u_pwp()
          u(npart_in:npart_fin)= get_u_pwp()
        ENDIF

      ENDIF

      END ASSOCIATE

    ENDDO matter_objects_loop

    !-------------------!
    ! Assignment of Ye  !
    !-------------------!

    IF(.NOT.ALLOCATED( THIS% Ye ))THEN
      ALLOCATE( THIS% Ye( THIS% npart ), STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array Ye ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array lapse_parts" )
    ENDIF

    assign_ye_on_particles: IF( THIS% compose_eos )THEN

      PRINT *, "Assigning electron fraction using the CompOSE file ", &
               TRIM(THIS% compose_path)//TRIM(THIS% compose_filename)

      compose_namefile= TRIM(THIS% compose_path)//TRIM(THIS% compose_filename)
      CALL THIS% read_compose_composition( compose_namefile )
      CALL THIS% compute_Ye()

      PRINT *, "Electron fraction assigned."
      PRINT *

    ELSE

      THIS% Ye= zero

    ENDIF assign_ye_on_particles
    Ye= THIS% Ye

    CALL THIS% sph_computer_timer% stop_timer()

    !
    !-- Printouts
    !
    !"(A28,E15.8,A10)"
    DO i_matter= 1, THIS% n_matter, 1

      ASSOCIATE( npart_in   => THIS% npart_i(i_matter-1) + 1, &
                 npart_fin  => THIS% npart_i(i_matter-1) +    &
                               THIS% npart_i(i_matter) )

      PRINT *, " * Maximum baryon density on object", i_matter, "=", &
                MAXVAL(THIS% baryon_density_parts(npart_in:npart_fin), DIM=1) &
                /((Msun_geo*km2m*m2cm)**3)*(amu), " g cm^{-3}"
                !*amu/(m2cm**3), " amu cm^{-3} (TODO: CHECK UNITS)"
      PRINT *, " * Minimum baryon density on object", i_matter, "=", &
                MINVAL( THIS% baryon_density_parts(npart_in:npart_fin), DIM=1) &
                /((Msun_geo*km2m*m2cm)**3)*(amu), " g cm^{-3}"
                !*amu/(m2cm**3), " amu cm^{-3} (TODO: CHECK UNITS)"
      PRINT *, " * Ratio between the two=", &
               MAXVAL(THIS% baryon_density_parts(npart_in:npart_fin), DIM=1)/ &
               MINVAL(THIS% baryon_density_parts(npart_in:npart_fin), DIM=1)
      PRINT *

      PRINT *, " * Maximum interpolated nlrf on object", i_matter, "=", &
               MAXVAL( THIS% nlrf_int(npart_in:npart_fin), DIM= 1 ) &
               /((Msun_geo*km2m*m2cm)**3), " baryon cm^{-3}"
               !" m0c2_cu (TODO: CHECK UNITS)"!, &
               !"baryon Msun_geo^{-3}"
      PRINT *, " * Minimum interpolated nlrf on object", i_matter, "=", &
               MINVAL( THIS% nlrf_int(npart_in:npart_fin), DIM= 1 ) &
               /((Msun_geo*km2m*m2cm)**3), " baryon cm^{-3}"
               !" m0c2_cu (TODO: CHECK UNITS)"!, &
               !"baryon Msun_geo^{-3}"
      PRINT *, " * Ratio between the two=", &
               MAXVAL( THIS% nlrf_int(npart_in:npart_fin), DIM= 1 )/ &
               MINVAL( THIS% nlrf_int(npart_in:npart_fin), DIM= 1 )
      PRINT *

      THIS% nuratio_i(i_matter)= MAXVAL( THIS% nu(npart_in:npart_fin), DIM= 1 )&
                                /MINVAL( THIS% nu(npart_in:npart_fin), DIM= 1 )
      PRINT *, " * Maximum n. baryon per particle (nu) on object", i_matter, &
                          "=", MAXVAL( THIS% nu(npart_in:npart_fin), DIM= 1 )
      PRINT *, " * Minimum n. baryon per particle (nu) on object", i_matter, &
                          "=", MINVAL( THIS% nu(npart_in:npart_fin), DIM= 1 )
      PRINT *, " * Ratio between the two=", THIS% nuratio_i(i_matter)
      PRINT *

      PRINT *, " * Number of baryons on object", i_matter, "=", &
               THIS% nbar_i(i_matter)
      PRINT *, " * Total mass of the baryons on object", i_matter, "=", &
               THIS% nbar_i(i_matter)*amu/Msun, "Msun =", &
               THIS% nbar_i(i_matter)*amu/Msun/THIS% masses(i_matter), &
               "of the baryon mass of object", i_matter, "."
      PRINT *

      END ASSOCIATE

    ENDDO

    THIS% nuratio= MAXVAL( THIS% nu, DIM= 1 )/MINVAL( THIS% nu, DIM= 1 )
    PRINT *, " * Baryon number ratio across the stars=", THIS% nuratio
    PRINT *
    PRINT *, " * Total mass of the baryons=", &
             THIS% nbar_tot*amu/Msun, "Msun =", &
             THIS% nbar_tot*amu/Msun/(SUM(THIS% masses, DIM=1)), &
             "of the total baryon mass."
    PRINT *

    !
    !-- Adjusting the baryon number per particle uniformly so that
    !-- the baryon mass is correct, but the ratio between nu_max and nu_min
    !-- does not change.
    !-- nlrf is not to be rescaled, according to Stephan, since:
    !--   (i)  it is directly computed from the LORENE ID and should therefore
    !--        be consistent with it
    !--   (ii) it is anyway immediately recomputed in SPHINCS_BSSN
    !
    IF( THIS% correct_nu )THEN

      THIS% nbar_tot= zero
      DO i_matter= 1, THIS% n_matter, 1

        ASSOCIATE( npart_in   => THIS% npart_i(i_matter-1) + 1, &
                   npart_fin  => THIS% npart_i(i_matter-1) +    &
                                 THIS% npart_i(i_matter) )

        THIS% nu( npart_in:npart_fin )= THIS% nu( npart_in:npart_fin ) &
                      /(THIS% nbar_i(i_matter)*amu/Msun/THIS% masses(i_matter))
        nu( npart_in:npart_fin )= THIS% nu( npart_in:npart_fin )

        THIS% nbar_i(i_matter)= THIS% nbar_i(i_matter) &
                      /(THIS% nbar_i(i_matter)*amu/Msun/THIS% masses(i_matter))

        THIS% nbar_tot= THIS% nbar_tot + THIS% nbar_i(i_matter)

        PRINT *, " * Number of corrected baryons on object", i_matter, "=", &
                 THIS% nbar_i(i_matter)
        PRINT *, " * Total mass of the corrected baryons object", i_matter, &
                 "=", THIS% nbar_i(i_matter)*amu/Msun, "Msun =", &
                 THIS% nbar_i(i_matter)*amu/Msun/THIS% masses(i_matter), &
                 "of the baryon mass of object", i_matter, "."

        END ASSOCIATE

      ENDDO

      PRINT *, " * Total number of corrected baryons=", THIS% nbar_tot
      PRINT *, " * Total mass of the corrected baryons=", &
               THIS% nbar_tot*amu/Msun, "Msun =", &
               THIS% nbar_tot*amu/Msun/(SUM(THIS% masses, DIM=1)), &
               "of the total baryon mass."
      PRINT *

    ENDIF

    !
    !-- Exporting the SPH ID to a binary file, for evolution
    !
    IF( THIS% export_bin )THEN

      IF( PRESENT(namefile) )THEN

        finalnamefile= TRIM( namefile ) // "00000"
        dcount= -1 ! since it is increased before writing
        CALL write_SPHINCS_dump( finalnamefile )

      ELSE

        basename= "NSNS."
        dcount= -1 ! since it is increased before writing
        CALL write_SPHINCS_dump()

      ENDIF

    ENDIF

    PRINT *, " * Computing particle number density by kernel interpolation..."
    PRINT *
    nu= one
    CALL density_loop( THIS% npart, THIS% pos, nu, h, &
                       THIS% particle_density_int )

    IF( debug ) PRINT *, "100"

    !CALL COM( THIS% npart, THIS% pos, THIS% nu, &
    !          com_x_newt, com_y_newt, com_z_newt, com_d_newt )
    CALL COM( THIS% npart_i(1), THIS% pos(:,1:THIS% npart_i(1)), &
              THIS% nu(1:THIS% npart_i(1)), &
              com_x_newt, com_y_newt, com_z_newt, com_d_newt )

    !CALL COM_1PN( THIS% npart, THIS% pos, THIS% v, &
    !              !THIS% v_euler_parts_x, &
    !              THIS% nu, THIS% baryon_density_parts, &
    !              THIS% specific_energy_parts, THIS% nstar_int, sq_detg4, gg4, &
    !              com_x_1pn, com_y_1pn, com_z_1pn, com_d_1pn )
    CALL COM_1PN( THIS% npart_i(1), THIS% pos(:,1:THIS% npart_i(1)), &
                  THIS% v(:,1:THIS% npart_i(1)), &
                  !THIS% v_euler_parts_x, &
                  THIS% nu(1:THIS% npart_i(1)), &
                  !THIS% baryon_density_parts(1:THIS% npart_i(1)), &
                  THIS% nlrf_int(1:THIS% npart_i(1)), &
                  THIS% u_pwp(1:THIS% npart_i(1)), &
                  THIS% nstar_int(1:THIS% npart_i(1)), sq_detg4, gg4, &
                  com_x_1pn, com_y_1pn, com_z_1pn, com_d_1pn )

    IF( debug ) PRINT *, "101"

    !CALL momentum_1pn( THIS% npart, THIS% pos, &
    !                   THIS% v, &
    !                   !THIS% v_euler_parts_x, &
    !                   THIS% nu, THIS% baryon_density_parts, &
    !                   THIS% specific_energy_parts, THIS% pressure_parts_cu, &
    !                   THIS% nstar_int, THIS% h, &
    !                   sq_detg4, gg4, px, py, pz )
    !CALL momentum_1pn( THIS% npart_i(1), THIS% pos(:,1:THIS% npart_i(1)), &
    !                   THIS% v(:,1:THIS% npart_i(1)), &
    !                   !THIS% v_euler_parts_x, &
    !                   THIS% nu(1:THIS% npart_i(1)), &
    !                   THIS% baryon_density_parts(1:THIS% npart_i(1)), &
    !                   THIS% specific_energy_parts(1:THIS% npart_i(1)), &
    !                   THIS% pressure_parts_cu(1:THIS% npart_i(1)), &
    !                   THIS% nstar_int(1:THIS% npart_i(1)), &
    !                   THIS% h(1:THIS% npart_i(1)), &
    !                   sq_detg4, gg4, px, py, pz )

    PRINT *, "LORENE COM:            ", &
             THIS% barycenter(1,:)! + THIS% barycenter(2,:)
    PRINT *, "Newtonian COM:         ", &
             com_x_newt, com_y_newt, com_z_newt!, com_d_newt
    PRINT *, "1PN COM:               ", &
             com_x_1pn, com_y_1pn, com_z_1pn!, com_d_1pn
    !PRINT *, "1PN spacetime momentum:", px, py, pz

    PRINT *, " * Deallocating MODULE variables..."
    PRINT *
    CALL deallocate_metric_on_particles
    CALL deallocate_gradient
    DEALLOCATE( alive )
    CALL deallocate_RCB_tree_memory_3D
    CALL deallocate_SPH_memory

    STOP

    call_flag= call_flag + 1
    THIS% call_flag= call_flag

    PRINT *, "** Subroutine compute_and_export_SPH_variables executed."
    PRINT *

  END PROCEDURE compute_and_export_SPH_variables


END SUBMODULE sph_variables

