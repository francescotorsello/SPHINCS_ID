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


  USE constants,  ONLY: zero, half, one, two, three, c_light2
  USE options,    ONLY: eos_str


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE compute_and_print_sph_variables

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
                                   Msun, zero, one
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
                                   u,     &  ! Specific internal energy in local
                                             ! rest frame (no kinetic energy)
                                   temp,  &  ! Temperature
                                   av,    &  ! Dissipation
                                   ye,    &  ! Electron fraction
                                   divv,  &  ! Divergence of velocity vel_u
                                   cs,    &  ! Sound speed
                                   allocate_SPH_memory, &
                                   deallocate_SPH_memory
    USE metric_on_particles, ONLY: allocate_metric_on_particles, &
                                   deallocate_metric_on_particles, &
                                   sq_det_g4, g4_ll
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
                                   gen_pwp_eos, gen_pwp_cold_eos, &
                                   get_u_pwp, shorten_eos_name, Gamma_th_1
    USE RCB_tree_3D,         ONLY: iorig, nic, nfinal, nprev, lpart, &
                                   rpart, allocate_RCB_tree_memory_3D, &
                                   deallocate_RCB_tree_memory_3D
    USE matrix,              ONLY: invert_3x3_matrix
    USE analyze,             ONLY: COM, lin_mom
    USE tensor,              ONLY: n_sym4x4, &
                                   itt, itx, ity, itz, &
                                   ixx, ixy, ixz, iyy, iyz, izz
    USE utility,             ONLY: compute_g4, determinant_sym4x4, &
                                   spacetime_vector_norm_sym4x4

    IMPLICIT NONE

    ! The flag call_flag is set different than 0 if the SUBROUTINE
    ! compute_and_print_sph_variables is called
    INTEGER, SAVE:: call_flag= 0

    !INTEGER, PARAMETER:: max_it_h= 1

    ! Spacetime indices \mu and \nu
    INTEGER:: nus, mus, cnt1, a, i_matter!, itr2, inde, index1!, cnt2
    INTEGER:: n_problematic_h
    INTEGER:: itot, l, ill!, b, k

    INTEGER:: pwp_gamma, pwp_kappa

    !DOUBLE PRECISION:: g4(0:3,0:3)
    DOUBLE PRECISION:: g4(n_sym4x4)
    DOUBLE PRECISION:: gg4(n_sym4x4,this% npart)
    DOUBLE PRECISION:: sq_detg4(this% npart)
    DOUBLE PRECISION:: det, sq_g, Theta_a!, &!nu_max1, nu_max2, &
                       !nu_tmp, nu_thres1, nu_thres2
    DOUBLE PRECISION:: com_x_newt, com_y_newt, com_z_newt, com_d_newt, mass_newt
    DOUBLE PRECISION:: com_x_1pn, com_y_1pn, com_z_1pn, com_d_1pn, mass_1pn
    DOUBLE PRECISION:: px_newt, py_newt, pz_newt, pnorm_newt
    DOUBLE PRECISION:: px, py, pz, pnorm, tmp

    !DOUBLE PRECISION:: ha, ha_1, ha_3, va, mat(3,3), mat_1(3,3), xa, ya, za
    !DOUBLE PRECISION:: mat_xx, mat_xy, mat_xz, mat_yy
    !DOUBLE PRECISION:: mat_yz, mat_zz, Wdx, Wdy, Wdz, dx, dy, dz, Wab, &
    !                   Wab_ha, Wi, Wi1, dvv

    LOGICAL:: few_ncand!, invertible_matrix

    LOGICAL, PARAMETER:: debug= .FALSE.

    CHARACTER( LEN= 2 ):: i_mat
    CHARACTER( LEN= : ), ALLOCATABLE:: compose_namefile
    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    TYPE(timer):: find_h_bruteforce_timer

    find_h_bruteforce_timer= timer( "find_h_bruteforce_timer" )

    PRINT *, "** Executing the compute_and_print_sph_variables " &
             // "subroutine..."
    PRINT *

    !
    !-- Set up the MODULE variables in MODULE sph_variables
    !-- (used by write_SPHINCS_dump)
    !
    npart= this% npart
    n1= this% npart_i(1)
    IF( this% n_matter == 2 ) n2= this% npart_i(2)

    CALL set_units('NSM')
    CALL read_options

    CALL allocate_SPH_memory
    CALL allocate_metric_on_particles( this% npart )

    IF( debug ) PRINT *, "1"

    IF(.NOT.ALLOCATED( this% nu ))THEN
      ALLOCATE( this% nu( this% npart ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array nu ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array nu" )
    ENDIF
    IF(.NOT.ALLOCATED( this% nlrf ))THEN
      ALLOCATE( this% nlrf( this% npart ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array nlrf ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array nlrf" )
    ENDIF
    IF(.NOT.ALLOCATED( this% Theta ))THEN
      ALLOCATE( this% Theta( this% npart ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array Theta ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array Theta" )
    ENDIF
    IF(.NOT.ALLOCATED( this% v ))THEN
      ALLOCATE( this% v( 0:3, this% npart ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array v ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array v" )
    ENDIF
    IF(.NOT.ALLOCATED( this% h ))THEN
      ALLOCATE( this% h( this% npart ), STAT= ios, &
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

    CALL this% sph_computer_timer% start_timer()
!    !$OMP PARALLEL DO DEFAULT( NONE ) &
!    !$OMP             SHARED( THIS, pos_u, vel_u, sq_det_g4, Theta, nlrf, &
!    !$OMP                     Pr, u, temp, av, divv ) &
!    !$OMP             PRIVATE( itr, g4, gg4, det, sq_g, Theta_a, &
!    !$OMP                      nus, mus )
    compute_SPH_variables_on_particles: DO itr= 1, this% npart, 1

      ! Particle positions [Msun_geo]
      pos_u(1,itr)= this% pos(1,itr)
      pos_u(2,itr)= this% pos(2,itr)
      pos_u(3,itr)= this% pos(3,itr)

      ! Coordinate velocity of the fluid [c]
      this% v(0,itr)= one
      this% v(1,itr)= this% lapse(itr)*this% v_euler_x(itr) &
                    - this% shift_x(itr)
      this% v(2,itr)= this% lapse(itr)*this% v_euler_y(itr) &
                    - this% shift_y(itr)
      this% v(3,itr)= this% lapse(itr)*this% v_euler_z(itr) &
                    - this% shift_z(itr)
      vel_u(1,itr)  = this% v(1,itr)
      vel_u(2,itr)  = this% v(2,itr)
      vel_u(3,itr)  = this% v(3,itr)

      CALL compute_g4( this% lapse(itr), &
            [this% shift_x(itr), this% shift_y(itr), this% shift_z(itr)], &
            [this% g_xx(itr), this% g_xy(itr), this% g_xz(itr), &
             this% g_yy(itr), this% g_yz(itr), this% g_zz(itr)], g4 )

      CALL determinant_sym4x4( g4, det )
      IF( ABS(det) < 1D-10 )THEN
          PRINT *, "The determinant of the spacetime metric is " &
                   // "effectively 0 at particle ", itr
          PRINT *
          STOP
      ELSEIF( det > 0 )THEN
          PRINT *, "The determinant of the spacetime metric is " &
                   // "positive at particle ", itr
          PRINT *
          STOP
      ENDIF
      sq_g= SQRT(-det)
      sq_det_g4(itr)= sq_g
      sq_detg4(itr)= sq_g

      !
      !-- Generalized Lorentz factor
      !
      Theta_a= zero
      CALL spacetime_vector_norm_sym4x4( g4, this% v(:,itr), Theta_a )
      IF( Theta_a > zero )THEN
        PRINT *, "** ERROR! The computing frame particle 4-velocity is ", &
                 "spacelike at particle ", itr
        PRINT *, " * Its norm is ", Theta_a
        PRINT *, " * Stopping.."
        PRINT *
        STOP
      ELSEIF( Theta_a == zero )THEN
        PRINT *, "** ERROR! The computing frame particle 4-velocity is ", &
                 "null at particle ", itr
        PRINT *, " * Its norm is ", Theta_a
        PRINT *, " * Stopping.."
        PRINT *
        STOP
      ENDIF
      Theta_a         = one/SQRT(-Theta_a)
      Theta(itr)      = Theta_a
      this% Theta(itr)= Theta_a

      ! Baryon density in the local rest frame [baryon (Msun_geo)^{-3}]
      ! Computed from the LORENE baryon mass density in [kg/m^3]
      nlrf(itr)= this% baryon_density(itr)
      this% nlrf(itr)= this% baryon_density(itr)

      ! Specific internal energy [c^2]
      u(itr)= this% specific_energy(itr)

      ! Pressure [amu*c**2/(Msun_geo**3)]
      !          dimensions: [(M/L**3)*L**2/T**2]= [M/(L*T**2)], same as
      !                      energy density
      Pr(itr)= this% pressure(itr)
      this% pressure_cu(itr)= Pr(itr)

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
    this% nstar= ( this% nlrf*this% Theta )*sq_det_g4

    !
    !-- Compute the particle proper mass, if not computed yet
    !
    IF( .NOT.ALLOCATED( this% pvol ) )THEN
      PRINT *, "** ERROR! The array pvol is not allocated. ", &
               "Stopping..."
      PRINT *
      STOP
    ENDIF
    IF( .NOT.( this% distribution_id == id_particles_on_spherical_surfaces &
        .OR. &
        ( this% distribution_id == id_particles_from_file &
          .AND. this% read_nu ) ) )THEN

      this% pmass= this% nstar*this% pvol

    ENDIF

    ! Compute particle number density from the ID
    this% particle_density= ( this% nstar )/( this% pmass )

    IF( debug ) PRINT *, "4"

    !
    !-- Compute the first guess for the smoothing length, if the APM was not
    !-- used
    !
    DO i_matter= 1, this% n_matter, 1

      IF( .NOT.this% apm_iterate(i_matter) )THEN

        IF( debug ) PRINT *, "Compute first guess for the smoothing length ", &
                             "h, for particles on matter object", itr,"..."

        compute_h_guess: &
        DO itr= this% npart_i(i_matter-1) + 1, &
                this% npart_i(i_matter-1) + this% npart_i(i_matter), 1

          this% h(itr)= three*(this% pvol(itr))**third
          h(itr)= this% h(itr)
          ! /(Msun_geo**3)
          IF( debug .AND. this% h(itr) <= zero )THEN
            PRINT *, "** ERROR! h(", itr, ")=", this% h(itr)
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

!    IF( this% redistribute_nu )THEN
!
!      !---------------------------------------------------------------------!
!      !--  Assignment of nu on the stars, with the purpose                --!
!      !--  of having a more uniform nu over the particles without losing  --!
!      !--  baryon mass. This is used only on the lattice, optionally.     --!
!      !---------------------------------------------------------------------!
!
!      IF( this% distribution_id == id_particles_on_spherical_surfaces )THEN
!        PRINT *, "** ERROR! Particle placer ", this% distribution_id, &
!                 " is not compatible with redistribute_nu= .TRUE."
!        PRINT *, " * Check the parameter file lorene_bns_id_particles.par. ", &
!                 "Stopping..."
!        PRINT *
!        STOP
!      ENDIF
!
!      nu_max1= nlrf( this% baryon_density_index( this% npart1 ) )&
!              *this% pvol( this% npart1 ) &
!              *Theta( this% baryon_density_index( this% npart1 ) )&
!              *sq_det_g4( this% baryon_density_index( this% npart1 ) )
!      nu_max2= nlrf( this% baryon_density_index( this% npart ) )&
!              *this% pvol( this% npart ) &
!              *Theta( this% baryon_density_index( this% npart ) )&
!              *sq_det_g4( this% baryon_density_index( this% npart ) )
!
!      nu_thres1= nu_max1/this% nu_ratio
!      nu_thres2= nu_max2/this% nu_ratio
!
!      ! Reset the total baryon number to 0 (necessary), and nu to an arbitrary
!      ! value (to make debugging easier)
!
!      nu= one
!      this% nu= one
!      this% nbar_tot= zero
!      this% nbar1= zero
!      this% nbar2= zero
!
!      cnt1= 0
!      compute_nu_on_particles_star1: DO itr= this% npart1, 1, -1
!
!        cnt1= cnt1 + 1
!
!        nu_tmp= nlrf( this% baryon_density_index( itr ) ) &
!                *this% pvol(itr) &
!                *Theta( this% baryon_density_index( itr ) )&
!                *sq_det_g4( this% baryon_density_index( itr ) )
!
!        !IF( itr == this% npart1 ) nu_max= nu_tmp ! move this out of the loop
!
!        IF( nu_tmp > nu_thres1 )THEN
!          nu( this% baryon_density_index( itr ) )      = nu_tmp
!          this% nu( this% baryon_density_index( itr ) )= nu_tmp
!        ELSE
!          nu( this% baryon_density_index( itr ) )      = nu_thres1
!          this% nu( this% baryon_density_index( itr ) )= nu_thres1
!        ENDIF
!
!        this% nbar1= this% nbar1 + &
!                     this% nu( this% baryon_density_index( itr ) )
!
!        IF( this% nbar1*amu/MSun > this% masses(1) )THEN
!          EXIT
!        ENDIF
!
!      ENDDO compute_nu_on_particles_star1
!
!      cnt2= 0
!      compute_nu_on_particles_star2: DO itr= this% npart, this% npart1 + 1, -1
!
!        cnt2= cnt2 + 1
!
!        nu_tmp= nlrf( this% baryon_density_index( itr ) ) &
!                *this% pvol(itr) &
!                *Theta( this% baryon_density_index( itr ) ) &
!                *sq_det_g4( this% baryon_density_index( itr ) )
!
!        !IF( itr == this% npart ) nu_max= nu_tmp
!
!        IF( nu_tmp > nu_thres2 )THEN
!          nu( this% baryon_density_index( itr ) )      = nu_tmp
!          this% nu( this% baryon_density_index( itr ) )= nu_tmp
!        ELSE
!          nu( this% baryon_density_index( itr ) )      = nu_thres2
!          this% nu( this% baryon_density_index( itr ) )= nu_thres2
!        ENDIF
!
!        this% nbar2= this% nbar2 + &
!                     this% nu( this% baryon_density_index( itr ) )
!
!        IF( this% nbar2*amu/MSun > this% masses(2) )THEN
!          EXIT
!        ENDIF
!
!      ENDDO compute_nu_on_particles_star2
!      this% nbar_tot= this% nbar1 + this% nbar2
!
!      !
!      !-- Reshape MODULE variables
!      !
!
!      CALL this% reshape_sph_field( pos_u, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( vel_u, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( Theta, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( h, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( nlrf, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( u, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( Pr, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( nu, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( temp, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( av, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( divv, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      !
!      !-- Reshape TYPE member SPH variables
!      !
!
!      CALL this% reshape_sph_field( this% pos, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% v, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% v_euler_x, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% v_euler_y, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% v_euler_z, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% Theta, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% h, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% baryon_density, cnt1, &
!                                    cnt2, this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% nlrf, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% energy_density, cnt1, &
!                                    cnt2, this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% specific_energy, cnt1, &
!                                    cnt2, this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% pressure, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% pressure_cu, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% nu, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% pvol, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      !
!      !-- Reshape TYPE member spacetime variables
!      !
!
!      CALL this% reshape_sph_field( this% lapse, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% shift_x, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% shift_y, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% shift_z, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% g_xx, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% g_xy, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% g_xz, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% g_yy, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% g_yz, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% g_zz, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      !
!      !-- Reassign particle numbers
!      !
!
!      npart= cnt1 + cnt2
!      this% npart= npart
!      this% npart1= cnt1
!      this% npart2= cnt2
!      n1= this% npart1
!      n2= this% npart2
!
!      PRINT *, " * Particles replaced after reassigning nu."
!      PRINT *, " * New number of particles=", this% npart
!      PRINT *
!      PRINT *, " * Number of particles on NS 1=", this% npart1
!      PRINT *, " * Number of particles on NS 2=", this% npart2
!      PRINT *

    !----------------------------------------------!
    !--  Assignment of nu on the matter objects. --!
    !----------------------------------------------!

    ! TODO: The code would be much cleaner if nu and h were assigned in the
    !       constructor for all the particle placers (now only the lattices
    !       and the particles ad from file without nu don't assign them)
    !       If that was the case, then this SUBROUTINE could start out by
    !       estimating the SPH density, and compute everything from there.

    DO i_matter= 1, this% n_matter, 1

      ASSOCIATE( npart_in   => this% npart_i(i_matter-1) + 1, &
                 npart_fin  => this% npart_i(i_matter-1) +    &
                               this% npart_i(i_matter) )

        IF( this% apm_iterate(i_matter) )THEN

          ! If the APM was used...

          ! Do nothing, nu is already computed and reflected in the constructor
          nu( npart_in : npart_fin )= this% nu( npart_in : npart_fin )
          this% nbar_i(i_matter)= SUM( this% nu( npart_in : npart_fin ), DIM=1 )

        ELSEIF( this% distribution_id == id_particles_from_file &
                .AND. this% read_nu )THEN

          ! If the particle positions and nu were read from formatted file...

          ! Do nothing, nu is already read from file
          nu( npart_in : npart_fin )= this% nu( npart_in : npart_fin )
          this% nbar_i(i_matter)= SUM( this% nu( npart_in : npart_fin ), DIM=1 )

        ELSEIF( this% distribution_id== id_particles_on_spherical_surfaces )THEN
        !ELSE

          ! If the APM was not used...

          ! Set nu based on the particle mass...

          DO itr= npart_in, npart_fin, 1
            nu(itr)= this% pmass(itr)*MSun/amu
            this% nu(itr)= nu(itr)
            this% nbar_i(i_matter)= this% nbar_i(i_matter) + nu(itr)
          ENDDO

        ELSE

          ! If the APM was not used and the particles are on lattices...

          DO itr= npart_in, npart_fin, 1
            nu(itr)= nlrf(itr)*this% pvol(itr)*Theta( itr )*sq_det_g4( itr )
            this% nu(itr)= nu(itr)
            this% nbar_i(i_matter)= this% nbar_i(i_matter) + nu(itr)
          ENDDO

        ENDIF
        this% nbar_tot= this% nbar_tot + this% nbar_i(i_matter)

        equal_mass_binary: &
        IF( i_matter == 1 .AND. this% n_matter == 2 )THEN

          IF( ABS(this% mass_ratios(1) - this% mass_ratios(2)) &
            /this% mass_ratios(2) <= 0.005 .AND. this% reflect_particles_x &
            .AND. &
            this% distribution_id /= id_particles_from_file )THEN

            ! Consistency check
            IF( this% npart_i(i_matter) +      &
            this% npart_i(i_matter+1) /= this% npart )THEN
              PRINT *, "** ERROR! npart_next /= this% npart! "
              PRINT *, "   npart_next=", this% npart_i(i_matter) +      &
              this% npart_i(i_matter+1)
              PRINT *, "   this% npart=", this% npart
              PRINT *, "   Stopping..."
              PRINT *
              STOP
            ENDIF

            nu( npart_fin + 1:this% npart_i(i_matter) +      &
            this% npart_i(i_matter+1) )= nu( npart_in : npart_fin )
            this% nu( npart_fin + 1:this% npart_i(i_matter) +      &
            this% npart_i(i_matter+1) )= nu( npart_in : npart_fin )
            this% nbar_i(i_matter + 1)= this% nbar_i(i_matter)

            this% nbar_tot= this% nbar_tot + this% nbar_i(i_matter + 1)

            ! Consistency check
            IF( this% nbar_tot /= 2*this% nbar_i(i_matter + 1) )THEN
              PRINT *, "** ERROR! this% nbar_tot /= 2*this% nbar(i_matter + 1) ",&
                       "   when reflecting particles or a binary system"
              PRINT *, "   this% nbar_tot=", this% nbar_tot
              PRINT *, "   2*this% nbar(", i_matter + 1, ")=", &
                           2*this% nbar_i(i_matter + 1)
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
                   this% npart, &
                   this% pos, this% h, & ! Input
                   h )             ! Output

    IF( debug ) PRINT *, "101.5"

    CALL find_h_bruteforce_timer% start_timer()
    n_problematic_h= 0
    check_h: DO a= 1, this% npart, 1

      IF( ISNAN( h(a) ) .OR. h(a) <= zero )THEN

        n_problematic_h= n_problematic_h + 1

        h(a)= find_h_backup( a, this% npart, this% pos, ndes )

        IF( ISNAN( h(a) ) .OR. h(a) <= zero )THEN
          PRINT *, "** ERROR! h=0 on particle ", a, "even with the brute", &
                   " force method."
          PRINT *, "   Particle position: ", this% pos(:,a)
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

    IF( SUM(this% nu, DIM= 1)/SIZE(this% nu) == zero )THEN
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
                                  this% npart, &
                                  this% pos,   &
                                  this% nu )

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
    check_h2: DO a= 1, this% npart, 1

      IF( ISNAN( h(a) ) .OR. h(a) <= zero )THEN

        n_problematic_h= n_problematic_h + 1
        h(a)= find_h_backup( a, this% npart, this% pos, ndes )
        PRINT *, h(a)
        IF( ISNAN( h(a) ) .OR. h(a) <= zero )THEN
          PRINT *, "** ERROR! h=0 on particle ", a, "even with the brute", &
                   " force method."
          PRINT *, "   Particle position: ", this% pos(:,a)
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
    this% h= h
    this% pvol= ( this% h/three )**three

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
    CALL density( this% npart, &
                  this% pos,   &
                  this% nstar_int )

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

    this% nlrf_int= ( this% nstar_int/this% Theta )/sq_det_g4
    nlrf= this% nlrf_int

    !-----------------------------------------------------------------------!
    ! For single and piecewise polytropes, do not use the pressure          !
    ! and specific internal energy from the ID.                             !
    ! Compute them using the exact formulas for piecewise                   !
    ! polytropic instead, starting from the kernel interpolated density     !
    !-----------------------------------------------------------------------!

    matter_objects_loop: DO i_matter= 1, this% n_matter, 1

      ASSOCIATE( npart_in  => this% npart_i(i_matter-1) + 1, &
                 npart_fin => this% npart_i(i_matter-1) +    &
                              this% npart_i(i_matter) )

      IF( this% all_eos(i_matter)% eos_parameters(1) == DBLE(1) )THEN
      ! If the |eos| is polytropic

        PRINT *, " * Computing pressure and specific internal energy from", &
                 " the baryon mass density, using the exact formulas for", &
                 " single polytropic EOS, on matter object", i_matter,"..."

        ! Formulas from Read et al. (2009), https://arxiv.org/abs/0812.2163

        IF( this% cold_system )THEN
        ! If the system is cold, compute pressure and specific energy
        ! exactly using the polytropic EOS

          PRINT *, " * Assuming a cold system: no thermal component considered."
          PRINT *

          Pr(npart_in:npart_fin)= &
            this% all_eos(i_matter)% eos_parameters($poly_kappa) &
            *( this% nlrf_int(npart_in:npart_fin)*m0c2_cu ) &
            **this% all_eos(i_matter)% eos_parameters($poly_gamma)

          ! Using this internal energy gives machine-precision relative errors
          ! after the recovery, since it is computed from nlrf_int
          ! Using the internal energy from the ID gives largr errors
          ! For the piecewise polytropes, we have (?) to use the energy
          ! from the ID
          u(npart_in:npart_fin)= ( Pr(npart_in:npart_fin) &
            /(this% nlrf_int(npart_in:npart_fin)*m0c2_cu &
            *( this% all_eos(i_matter)% eos_parameters($poly_gamma) - one ) ) )

          this% enthalpy(npart_in:npart_fin)= one + u(npart_in:npart_fin) &
            + this% nlrf_int(npart_in:npart_fin)*m0c2_cu/Pr(npart_in:npart_fin)

          cs(npart_in:npart_fin)= SQRT( &
            this% all_eos(i_matter)% eos_parameters($poly_gamma) &
              *Pr(npart_in:npart_fin)/ &
            (this% nlrf_int(npart_in:npart_fin)*m0c2_cu &
            *this% enthalpy(npart_in:npart_fin)) )

          !
          !-- Leaving the following code here, commented, because it allows
          !-- to test the pwp_eos MODULE using single polytropes
          !-- All tests were passed on 23.02.2022
          !
    !      CALL select_EOS_parameters( 'soft' )
    !
    !      DO a= npart_in, npart_fin, 1
    !
    !        CALL gen_pwp_cold_eos( this% nlrf_int(a)*m0c2_cu, &
    !                               Pr(a), u(a), cs(a) )
    !
    !        !CALL gen_pwp_eos( this% nlrf_int(a)*m0c2_cu, &
    !        !                  this% u_pwp(a), tmp, &
    !        !                  u(a), &
    !        !                  Pr(a), cs(a) )
    !      ENDDO

          Pr(npart_in:npart_fin)= Pr(npart_in:npart_fin)/m0c2_cu
          this% pressure_cu(npart_in:npart_fin)= Pr(npart_in:npart_fin)
          this% u_pwp(npart_in:npart_fin)= u(npart_in:npart_fin)

        ELSE
        ! If the system is hot, that is, has a thermal component, then
        ! the density and the specific energy (the latter including both
        ! cold and thermal part) should be supplied in the ID.
        ! The pressure is computed using them (see pwp_EOS MODULE).

          PRINT *, " * Assuming a hot system: thermal component considered."
          PRINT *

          u(npart_in:npart_fin)= this% specific_energy(npart_in:npart_fin)

          DO a= npart_in, npart_fin, 1

            Pr(a)= &
            ! cold pressure
            this% all_eos(i_matter)% eos_parameters($poly_kappa) &
              *( this% nlrf_int(a)*m0c2_cu ) &
              **this% all_eos(i_matter)% eos_parameters($poly_gamma) &
            + &
            ! thermal pressure
            Gamma_th_1*( this% nlrf_int(a)*m0c2_cu )* &
              MAX(u(a) - ( Pr(a)/(this% nlrf_int(a)*m0c2_cu &
                *( this% all_eos(i_matter)% eos_parameters($poly_gamma) &
                   - one ) ) ), zero)

          ENDDO
          this% enthalpy(npart_in:npart_fin)= one + u(npart_in:npart_fin) &
            + this% nlrf_int(npart_in:npart_fin)*m0c2_cu/Pr(npart_in:npart_fin)

          cs(npart_in:npart_fin)= SQRT( &
            this% all_eos(i_matter)% eos_parameters($poly_gamma) &
              *Pr(npart_in:npart_fin)/ &
            (this% nlrf_int(npart_in:npart_fin)*m0c2_cu &
            *this% enthalpy(npart_in:npart_fin)) )

          Pr(npart_in:npart_fin)= Pr(npart_in:npart_fin)/m0c2_cu
          this% pressure_cu(npart_in:npart_fin)= Pr(npart_in:npart_fin)
          this% u_pwp(npart_in:npart_fin)= u(npart_in:npart_fin)

        ENDIF

      ELSEIF( this% all_eos(i_matter)% eos_parameters(1) == DBLE(110) )THEN
      ! If the |eos| is piecewise polytropic

        PRINT *, " * Computing pressure and specific internal energy from", &
                 " the baryon mass density, using the exact formulas for", &
                 " piecewise polytropic EOS..."
        PRINT *

        IF( this% cold_system )THEN
        ! If the system is cold, compute pressure and specific energy
        ! exactly using the piecewise polytropic EOS

          PRINT *, " * Assuming a cold system: no thermal component considered."
          PRINT *

          CALL select_EOS_parameters( &
                        shorten_eos_name(this% all_eos(i_matter)% eos_name) )

          PRINT *, this% all_eos(i_matter)% eos_name
          PRINT *, shorten_eos_name(this% all_eos(i_matter)% eos_name)
          PRINT *, eos_str
          STOP

          DO a= npart_in, npart_fin, 1

            CALL gen_pwp_cold_eos( this% nlrf_int(a)*m0c2_cu, &
                                   Pr(a), u(a), cs(a) )

          ENDDO
          Pr(npart_in:npart_fin)= Pr(npart_in:npart_fin)/m0c2_cu
          this% pressure_cu(npart_in:npart_fin)= Pr(npart_in:npart_fin)
          this% u_pwp(npart_in:npart_fin)= u(npart_in:npart_fin)

        ELSE
        ! If the system is hot, that is, has a thermal component, then
        ! the density and the specific energy (the latter including both
        ! cold and thermal part) should be supplied in the ID.
        ! The pressure is computed using them (see pwp_EOS MODULE).

          PRINT *, " * Assuming a hot system: thermal component considered."
          PRINT *

          u(npart_in:npart_fin)= this% specific_energy(npart_in:npart_fin)

          CALL select_EOS_parameters( &
                        shorten_eos_name(this% all_eos(i_matter)% eos_name) )

          DO a= npart_in, npart_fin, 1

            CALL gen_pwp_eos( this% nlrf_int(a)*m0c2_cu, &
                              this% u_pwp(a), tmp, &
                              u(a), &
                              Pr(a), cs(a) )

          ENDDO
          Pr(npart_in:npart_fin)= Pr(npart_in:npart_fin)/m0c2_cu
          this% pressure_cu(npart_in:npart_fin)= Pr(npart_in:npart_fin)
          this% u_pwp(npart_in:npart_fin)= u(npart_in:npart_fin)

        ENDIF

     !   CALL select_EOS_parameters( &
     !           shorten_eos_name(this% all_eos(i_matter)% eos_name) )
     !
     !   !CALL gen_pwp_eos_all( this% npart_i(i_matter), &
     !   !                      this% nlrf_int(npart_in:npart_fin)*m0c2_cu, &
     !   !                      u(npart_in:npart_fin) )
     !
     !   DO a= npart_in, npart_fin, 1
     !
     !     CALL gen_pwp_eos( this% nlrf_int(a)*m0c2_cu, &
     !                       this% u_pwp(a), tmp, &
     !                       this% specific_energy(a), Pr(a), cs(a) )
     !
     !   ENDDO
     !
     !   Pr(npart_in:npart_fin)= Pr(npart_in:npart_fin)/m0c2_cu
     !   this% pressure_cu(npart_in:npart_fin)= Pr(npart_in:npart_fin)
     !
     !   IF( this% cold_system )THEN
     !   ! If the system is cold, get the specific internal energy computed
     !   ! exactly using the piecewise polytropic EOS
     !     PRINT *, " * Assuming a cold system: no thermal component."
     !     PRINT *
     !     u(npart_in:npart_fin)= this% u_pwp(npart_in:npart_fin)
     !   ELSE
     !   ! Otherwise, get the specific nternal nergy from the ID
     !     PRINT *, " * Assuming a hot system: thermal component added."
     !     PRINT *
     !     u(npart_in:npart_fin)= this% specific_energy(npart_in:npart_fin)
     !   ENDIF

      ENDIF

      END ASSOCIATE

    ENDDO matter_objects_loop

    !-------------------!
    ! Assignment of Ye  !
    !-------------------!

    IF(.NOT.ALLOCATED( this% Ye ))THEN
      ALLOCATE( this% Ye( this% npart ), STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array Ye ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array lapse" )
    ENDIF

    assign_ye_on_particles: IF( this% compose_eos )THEN

      PRINT *, "Assigning electron fraction using the CompOSE file ", &
               TRIM(this% compose_path)//TRIM(this% compose_filename)

      compose_namefile= TRIM(this% compose_path)//TRIM(this% compose_filename)
      CALL this% read_compose_composition( compose_namefile )
      CALL this% compute_Ye()

      PRINT *, "Electron fraction assigned."
      PRINT *

    ELSE

      this% Ye= zero

    ENDIF assign_ye_on_particles
    Ye= this% Ye

    CALL this% sph_computer_timer% stop_timer()

    !
    !-- Printouts
    !
    !"(A28,E15.8,A10)"
    DO i_matter= 1, this% n_matter, 1

      ASSOCIATE( npart_in   => this% npart_i(i_matter-1) + 1, &
                 npart_fin  => this% npart_i(i_matter-1) +    &
                               this% npart_i(i_matter) )

      PRINT *, " * Maximum baryon density on object", i_matter, "=", &
                MAXVAL(this% baryon_density(npart_in:npart_fin), DIM=1) &
                /((Msun_geo*km2m*m2cm)**3)*(amu), " g cm^{-3}"
                !*amu/(m2cm**3), " amu cm^{-3} (TODO: CHECK UNITS)"
      PRINT *, " * Minimum baryon density on object", i_matter, "=", &
                MINVAL( this% baryon_density(npart_in:npart_fin), DIM=1) &
                /((Msun_geo*km2m*m2cm)**3)*(amu), " g cm^{-3}"
                !*amu/(m2cm**3), " amu cm^{-3} (TODO: CHECK UNITS)"
      PRINT *, " * Ratio between the two=", &
               MAXVAL(this% baryon_density(npart_in:npart_fin), DIM=1)/ &
               MINVAL(this% baryon_density(npart_in:npart_fin), DIM=1)
      PRINT *

      PRINT *, " * Maximum interpolated nlrf on object", i_matter, "=", &
               MAXVAL( this% nlrf_int(npart_in:npart_fin), DIM= 1 ) &
               /((Msun_geo*km2m*m2cm)**3), " baryon cm^{-3}"
               !" m0c2_cu (TODO: CHECK UNITS)"!, &
               !"baryon Msun_geo^{-3}"
      PRINT *, " * Minimum interpolated nlrf on object", i_matter, "=", &
               MINVAL( this% nlrf_int(npart_in:npart_fin), DIM= 1 ) &
               /((Msun_geo*km2m*m2cm)**3), " baryon cm^{-3}"
               !" m0c2_cu (TODO: CHECK UNITS)"!, &
               !"baryon Msun_geo^{-3}"
      PRINT *, " * Ratio between the two=", &
               MAXVAL( this% nlrf_int(npart_in:npart_fin), DIM= 1 )/ &
               MINVAL( this% nlrf_int(npart_in:npart_fin), DIM= 1 )
      PRINT *

      PRINT *, " * Maximum pressure", i_matter, "=", &
               MAXVAL( this% pressure_cu(npart_in:npart_fin), DIM= 1 ) &
               *amu*c_light2/((Msun_geo*km2m*m2cm)**3), " Ba"
               !" m0c2_cu (TODO: CHECK UNITS)"!, &
               !"baryon Msun_geo^{-3}"
      PRINT *, " * Minimum pressure", i_matter, "=", &
               MINVAL( this% pressure_cu(npart_in:npart_fin), DIM= 1 ) &
               *amu*c_light2/((Msun_geo*km2m*m2cm)**3), " Ba"
               !" m0c2_cu (TODO: CHECK UNITS)"!, &
               !"baryon Msun_geo^{-3}"
      PRINT *, " * Ratio between the two=", &
               MAXVAL( this% pressure_cu(npart_in:npart_fin), DIM= 1 )/ &
               MINVAL( this% pressure_cu(npart_in:npart_fin), DIM= 1 )
      PRINT *

      PRINT *, " * Maximum specific internal energy", i_matter, "=", &
               MAXVAL( this% u_pwp(npart_in:npart_fin), DIM= 1 ), " c^2"
               !" m0c2_cu (TODO: CHECK UNITS)"!, &
               !"baryon Msun_geo^{-3}"
      PRINT *, " * Minimum specific internal energy", i_matter, "=", &
               MINVAL( this% u_pwp(npart_in:npart_fin), DIM= 1 ), " c^2"
               !" m0c2_cu (TODO: CHECK UNITS)"!, &
               !"baryon Msun_geo^{-3}"
      PRINT *, " * Ratio between the two=", &
               MAXVAL( this% u_pwp(npart_in:npart_fin), DIM= 1 )/ &
               MINVAL( this% u_pwp(npart_in:npart_fin), DIM= 1 )
      PRINT *

      this% nuratio_i(i_matter)= MAXVAL( this% nu(npart_in:npart_fin), DIM= 1 )&
                                /MINVAL( this% nu(npart_in:npart_fin), DIM= 1 )
      PRINT *, " * Maximum n. baryon per particle (nu) on object", i_matter, &
                          "=", MAXVAL( this% nu(npart_in:npart_fin), DIM= 1 )
      PRINT *, " * Minimum n. baryon per particle (nu) on object", i_matter, &
                          "=", MINVAL( this% nu(npart_in:npart_fin), DIM= 1 )
      PRINT *, " * Ratio between the two=", this% nuratio_i(i_matter)
      PRINT *

      PRINT *, " * Number of baryons on object", i_matter, "=", &
               this% nbar_i(i_matter)
      PRINT *, " * Total mass of the baryons on object", i_matter, "=", &
               this% nbar_i(i_matter)*amu/Msun, "Msun =", &
               this% nbar_i(i_matter)*amu/Msun/this% masses(i_matter), &
               "of the baryon mass of object", i_matter, "."
      PRINT *

      END ASSOCIATE

    ENDDO

    this% nuratio= MAXVAL( this% nu, DIM= 1 )/MINVAL( this% nu, DIM= 1 )
    PRINT *, " * Baryon number ratio across the stars=", this% nuratio
    PRINT *
    PRINT *, " * Total mass of the baryons=", &
             this% nbar_tot*amu/Msun, "Msun =", &
             this% nbar_tot*amu/Msun/(SUM(this% masses, DIM=1)), &
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
    IF( this% correct_nu )THEN

      this% nbar_tot= zero
      DO i_matter= 1, this% n_matter, 1

        ASSOCIATE( npart_in   => this% npart_i(i_matter-1) + 1, &
                   npart_fin  => this% npart_i(i_matter-1) +    &
                                 this% npart_i(i_matter) )

        this% nu( npart_in:npart_fin )= this% nu( npart_in:npart_fin ) &
                      /(this% nbar_i(i_matter)*amu/Msun/this% masses(i_matter))
        nu( npart_in:npart_fin )= this% nu( npart_in:npart_fin )

        this% nbar_i(i_matter)= this% nbar_i(i_matter) &
                      /(this% nbar_i(i_matter)*amu/Msun/this% masses(i_matter))

        this% nbar_tot= this% nbar_tot + this% nbar_i(i_matter)

        PRINT *, " * Number of corrected baryons on object", i_matter, "=", &
                 this% nbar_i(i_matter)
        PRINT *, " * Total mass of the corrected baryons object", i_matter, &
                 "=", this% nbar_i(i_matter)*amu/Msun, "Msun =", &
                 this% nbar_i(i_matter)*amu/Msun/this% masses(i_matter), &
                 "of the baryon mass of object", i_matter, "."

        END ASSOCIATE

      ENDDO

      PRINT *, " * Total number of corrected baryons=", this% nbar_tot
      PRINT *, " * Total mass of the corrected baryons=", &
               this% nbar_tot*amu/Msun, "Msun =", &
               this% nbar_tot*amu/Msun/(SUM(this% masses, DIM=1)), &
               "of the total baryon mass."
      PRINT *

    ENDIF

    !
    !-- Exporting the SPH ID to a binary file, for evolution
    !
    IF( this% export_bin )THEN

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

    !
    !-- Test the recovery
    !

    ! Assign values to the MODULE variable
    !  DO itr= 1, this% npart, 1
    !
    !    CALL compute_g4( this% lapse(itr), &
    !          [this% shift_x(itr), this% shift_y(itr), this% shift_z(itr)], &
    !          [this% g_xx(itr), this% g_xy(itr), this% g_xz(itr), &
    !           this% g_yy(itr), this% g_yz(itr), this% g_zz(itr)], &
    !           g4_ll(1:n_sym4x4,itr) )
    !
    !  ENDDO

    CALL this% test_recovery( this% npart,       &
                              this% pos,         &
                              this% nlrf_int,    &
                              this% u_pwp,       &
                              this% pressure_cu, &
                              this% v(1:3,:),    &
                              this% theta,       &
                              this% nstar_int )

    ! Test the recovery on ech matter object separately
    ! DO i_matter= 1, this% n_matter, 1
    !
    !   PRINT *, " * Testing recovery on matter object", i_matter, "..."
    !   PRINT *
    !
    !   IF( i_matter > 9 )THEN
    !     WRITE( i_mat, "(I2)" ) i_matter
    !   ELSE
    !     WRITE( i_mat, "(I1)" ) i_matter
    !   ENDIF
    !   finalnamefile= "recovery_test-"//TRIM(i_mat)//".dat"
    !
    !   ASSOCIATE( npart_in   => this% npart_i(i_matter-1) + 1, &
    !              npart_fin  => this% npart_i(i_matter-1) +    &
    !                            this% npart_i(i_matter) )
    !
    !   CALL this% test_recovery( this% npart_i    (i_matter),               &
    !                             this% pos        (:,npart_in:npart_fin),   &
    !                             this% nlrf_int   (npart_in:npart_fin),     &
    !                             this% u_pwp      (npart_in:npart_fin),     &
    !                             this% pressure_cu(npart_in:npart_fin),     &
    !                             this% v          (1:3,npart_in:npart_fin), &
    !                             this% theta      (npart_in:npart_fin),     &
    !                             this% nstar_int  (npart_in:npart_fin),     &
    !                             finalnamefile )
    !
    !   END ASSOCIATE
    !
    ! ENDDO

    PRINT *, " * Computing particle number density by kernel interpolation..."
    PRINT *
    nu= one
    CALL density_loop( this% npart, this% pos, nu, h, &
                       this% particle_density_int )

    IF( debug ) PRINT *, "100"

    !CALL COM( this% npart, this% pos, this% nu, &
    !          com_x_newt, com_y_newt, com_z_newt, com_d_newt )
    !CALL COM( this% npart_i(1), this% pos(:,1:this% npart_i(1)), &
    !          this% nu(1:this% npart_i(1)), &
    !          com_x_newt, com_y_newt, com_z_newt, com_d_newt )

    !CALL COM_1PN( this% npart, this% pos, this% v, &
    !              !this% v_euler_x, &
    !              this% nu, this% baryon_density, &
    !              this% specific_energy, this% nstar_int, sq_detg4, gg4, &
    !              com_x_1pn, com_y_1pn, com_z_1pn, com_d_1pn )
    !CALL COM_1PN( this% npart_i(1), this% pos(:,1:this% npart_i(1)), &
    !              this% v(1:3,1:this% npart_i(1)), &
    !              !this% v_euler_x, &
    !              this% nu(1:this% npart_i(1)), &!/sq_det_g4(1:this% npart_i(1))/this% Theta(1:this% npart_i(1)), &
    !              !this% baryon_density(1:this% npart_i(1)), &
    !              this% nlrf_int(1:this% npart_i(1)), &
    !              this% u_pwp(1:this% npart_i(1)), &
    !              this% nstar_int(1:this% npart_i(1)), &
    !              sq_det_g4(1:this% npart_i(1)), &
    !              gg4(:,1:this% npart_i(1)), &
    !              com_x_1pn, com_y_1pn, com_z_1pn, com_d_1pn, mass_1pn )

    IF( debug ) PRINT *, "101"

    !nu   = this% nu
    !vel_u= this% v(1:3,:)
    !CALL lin_mom( pnorm_newt, px_newt, py_newt, pz_newt )
    !CALL momentum_1pn( this% npart, this% pos, &
    !                   this% v, &
    !                   !this% v_euler_x, &
    !                   this% nu, this% baryon_density, &
    !                   this% specific_energy, this% pressure_cu, &
    !                   this% nstar_int, this% h, &
    !                   sq_detg4, gg4, px, py, pz )
    !CALL momentum_1pn( this% npart, this% pos, &
    !                   this% v, &
    !                   !this% v_euler_x, &
    !                   this% nu, &
    !                   this% nlrf_int, &
    !                   this% u_pwp, &
    !                   this% pressure_cu, &
    !                   this% nstar_int, &
    !                   this% h, &
    !                   sq_detg4, gg4, px, py, pz )

    !PRINT *, "LORENE mass:            ", this% masses(1), mass_1pn*amu/Msun
    !PRINT *, "LORENE COM:            ", &
    !         this% barycenter(1,:)! + this% barycenter(2,:)
    !PRINT *, "Newtonian COM:         ", &
    !         com_x_newt, com_y_newt, com_z_newt!, com_d_newt
    !PRINT *, "1PN COM:               ", &
    !         com_x_1pn, com_y_1pn, com_z_1pn!, com_d_1pn
    !PRINT *, "Newtonian spacetime momentum:", px_newt/SUM(nu,DIM=1), py_newt/SUM(nu,DIM=1), pz_newt/SUM(nu,DIM=1), pnorm_newt/SUM(nu,DIM=1)
    !PRINT *, "1PN spacetime momentum:", px/mass_1pn, py/mass_1pn, pz/mass_1pn, pnorm/mass_1pn
    !PRINT *
    !PRINT *

    !
    !-- Deallocate MODULE variables
    !
    PRINT *, " * Deallocating MODULE variables..."
    PRINT *
    IF( ALLOCATED(g4_ll) ) CALL deallocate_metric_on_particles
    CALL deallocate_gradient
    DEALLOCATE( alive )
    CALL deallocate_RCB_tree_memory_3D
    CALL deallocate_SPH_memory

    !STOP

    call_flag= call_flag + 1
    this% call_flag= call_flag

    PRINT *, "** Subroutine compute_and_print_sph_variables executed."
    PRINT *

  END PROCEDURE compute_and_print_sph_variables


END SUBMODULE sph_variables

