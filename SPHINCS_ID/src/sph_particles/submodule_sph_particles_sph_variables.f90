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
  !# This SUBMODULE contains the implementation of
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


  USE constants,  ONLY: half, c_light2
  USE utility,    ONLY: zero, one, two, three
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

    USE constants,           ONLY: amu, third, Msun
    USE utility,             ONLY: MSun_geo, k_lorene2hydrobase, zero, one, &
                                   compute_g4, determinant_sym4x4, &
                                   spacetime_vector_norm_sym4x4, km2m, m2cm, &
                                   sph_path
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
    USE APM,                 ONLY: density_loop
    USE kernel_table,        ONLY: ktable
    USE options,             ONLY: ndes
    USE set_h,               ONLY: exact_nei_tree_update
    USE gradient,            ONLY: allocate_gradient, deallocate_gradient
    USE sphincs_sph,         ONLY: ncand
    USE alive_flag,          ONLY: alive
    USE APM,                 ONLY: assign_h
    USE pwp_EOS,             ONLY: select_EOS_parameters, gen_pwp_cold_eos, &
                                   gen_pwp_eos, Gamma_th_1
    USE RCB_tree_3D,         ONLY: iorig, nic, nfinal, nprev, lpart, &
                                   rpart, allocate_RCB_tree_memory_3D, &
                                   deallocate_RCB_tree_memory_3D
    USE matrix,              ONLY: invert_3x3_matrix
    USE analyze,             ONLY: COM, lin_mom
    USE tensor,              ONLY: n_sym4x4, &
                                   itt, itx, ity, itz, &
                                   ixx, ixy, ixz, iyy, iyz, izz

    IMPLICIT NONE

    ! The flag call_flag is set different than 0 if the SUBROUTINE
    ! compute_and_print_sph_variables is called
    INTEGER, SAVE:: call_flag= 0

    INTEGER:: cnt1, a, i_matter
    INTEGER:: n_problematic_h
    INTEGER:: itot, l, ill

    DOUBLE PRECISION:: g4(n_sym4x4)
    DOUBLE PRECISION:: det, Theta_a, tmp

    LOGICAL:: few_ncand, exist

    LOGICAL, PARAMETER:: debug= .FALSE.

    CHARACTER( LEN= : ), ALLOCATABLE:: compose_namefile
    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    TYPE(timer):: find_h_bruteforce_timer

    find_h_bruteforce_timer= timer( "find_h_bruteforce_timer" )

    CALL this% sph_computer_timer% start_timer()

    PRINT *, "** Executing the compute_and_print_sph_variables " &
             // "subroutine..."
    PRINT *

    !
    !-- Set up the MODULE variables from MODULE sph_variables
    !-- (used by write_SPHINCS_dump).
    !-- If write_SPHINCS_dump had arguments, it would not be needed
    !-- to assign values to the MODULE variables.
    !-- The TYPE-bound variables could be passed as arguments.
    !
    npart= this% npart
    n1   = this% npart_i(1)
    IF( this% n_matter == 2 ) n2= this% npart_i(2)

    CALL set_units('NSM')
    CALL read_options

    CALL allocate_SPH_memory
    CALL allocate_metric_on_particles( this% npart )

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( this, pos_u, h, nu ) &
    !$OMP             PRIVATE( a )
    assign_module_variables: DO a= 1, this% npart, 1
      pos_u(:,a)= this% pos(:,a)
      h    (a)  = this% h  (a)
      nu   (a)  = this% nu (a)
    ENDDO assign_module_variables
    !$OMP END PARALLEL DO

    IF( debug ) PRINT *, "1"

    !---------------------------------------------------------------!
    ! Compute SPH density estimate after updating the neighbor tree !
    ! (nu has to be assigned before this step)  !                   !
    !---------------------------------------------------------------!

    CALL allocate_RCB_tree_memory_3D(npart)
    iorig(1:npart)= (/ (a,a=1,npart) /)

    CALL allocate_gradient( npart )

    IF( debug ) PRINT *, "-1"

    PRINT *, " * Assigning h..."
    PRINT *

    ! Determine smoothing length so that each particle has exactly
    ! ndes neighbours
    CALL assign_h( ndes, &
                   this% npart, &
                   this% pos, this% h, & ! Input
                   h )                   ! Output

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
               "TYPE-bound array?"
      PRINT *, "Stopping..."
      STOP
    ENDIF

    CALL exact_nei_tree_update( ndes,        &
                                this% npart, &
                                this% pos,   &
                                this% nu )

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

    ! Update the member variable storing the smoothing length
    this% h= h

    !
    !-- Compute the SPH estimate of the relativistic density variable nstar
    !
    PRINT *, " * Computing the SPH estimate of the relativistic density ", &
             "variable nstar..."
    PRINT *

    CALL density_loop( this% npart,     &
                       this% pos,       &
                       this% nu,        &
                       this% h,         &
                       this% nstar_sph )

    IF( debug ) PRINT *, "20"

    !
    !-- Compute square root of minus the metric determinant,
    !-- and the generalized Lorentz factor.
    !-- They are needed to compute nlrf from the SPH estimate of nstar
    !
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( this, sq_det_g4 ) &
    !$OMP             PRIVATE( a, g4, det )
    compute_sph_variables_on_particles: DO a= 1, this% npart, 1

      ! Coordinate velocity of the fluid [c]
      this% v(0,a)= one
      this% v(1,a)= this% lapse(a)*this% v_euler_x(a) &
                    - this% shift_x(a)
      this% v(2,a)= this% lapse(a)*this% v_euler_y(a) &
                    - this% shift_y(a)
      this% v(3,a)= this% lapse(a)*this% v_euler_z(a) &
                    - this% shift_z(a)

      CALL compute_g4( this% lapse(a), &
                       [this% shift_x(a), this% shift_y(a), this% shift_z(a)], &
                       [this% g_xx(a), this% g_xy(a), this% g_xz(a), &
                        this% g_yy(a), this% g_yz(a), this% g_zz(a)], g4 )

      CALL determinant_sym4x4( g4, det )
      IF( ABS(det) < 1D-10 )THEN
        PRINT *, "The determinant of the spacetime metric is " &
                 // "effectively 0 at particle ", a
        PRINT *
        STOP
      ELSEIF( det > 0 )THEN
        PRINT *, "The determinant of the spacetime metric is " &
                 // "positive at particle ", a
        PRINT *
        STOP
      ENDIF
      sq_det_g4(a)= SQRT(-det)

      !
      !-- Generalized Lorentz factor
      !
      this% Theta(a)= zero
      CALL spacetime_vector_norm_sym4x4( g4, this% v(:,a), this% Theta(a) )
      IF( this% Theta(a) > zero )THEN
        PRINT *, "** ERROR! The computing frame particle 4-velocity is ", &
                 "spacelike at particle ", a
        PRINT *, " * Its norm is ",this% Theta(a)
        PRINT *, " * Stopping.."
        PRINT *
        STOP
      ELSEIF( this% Theta(a) == zero )THEN
        PRINT *, "** ERROR! The computing frame particle 4-velocity is ", &
                 "null at particle ", a
        PRINT *, " * Its norm is ", this% Theta(a)
        PRINT *, " * Stopping.."
        PRINT *
        STOP
      ENDIF
      this% Theta(a)= one/SQRT(-this% Theta(a))

    ENDDO compute_sph_variables_on_particles
    !$OMP END PARALLEL DO

    IF( debug ) PRINT *, "21"

    !-------------------------------------------------------------------------!
    !-------------------------------------------------------------------------!
    ! This point here is CRUCIAL: the particle distribution may NOT resolve   !
    ! properly the steep density gradient at the surface, even if the APM     !
    ! is used. This implies that the kernel interpolated nstar_sph will be    !
    ! different than nstar close to the surface.                              !
    ! The error can be such that the recovery fails in SPHINCS_BSSN, and      !
    ! this is of course a problem. Now, the density used in SPHINCS_BSSN      !
    ! during the evolution is not the one given in the ID. Hence, once        !
    ! nstar_sph is computed, nlrf should be recomputed from it, so that the   !
    ! density on the particles corresponds to the density that "they see",    !
    ! that is, the kernel interpolated density that uses these values of nu.  !
    !-------------------------------------------------------------------------!
    !-------------------------------------------------------------------------!

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( this, nlrf, sq_det_g4 ) &
    !$OMP             PRIVATE( a )
    compute_nlrf: DO a= 1, this% npart, 1

      ! Baryon mass density in the local rest frame
      this% nlrf_sph(a)= (this% nstar_sph(a)/this% Theta(a))/sq_det_g4(a)
      nlrf(a)          = this% nlrf_sph(a)
      ! Particle volume
      this% pvol(a)            = this% nu(a)/this% nstar_sph(a)
      ! Particle number density
      this% particle_density(a)= this% nstar_sph(a)/this% nu(a)

    ENDDO compute_nlrf
    !$OMP END PARALLEL DO

    IF( debug ) PRINT *, "22"

    !
    !-- Assign remaining particle properties to MODULE variables
    !
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( this, vel_u, Theta, temp, av, divv ) &
    !$OMP             PRIVATE( a )
    assign_module_variables2: DO a= 1, this% npart, 1
      ! Computing frame spatial velocity
      vel_u(:,a)= this% v(1:3,a)
      ! Generalized Lorentz factor
      Theta(a)  = this% Theta(a)
      ! Temperature
      temp(a)   = zero
      ! Disspation parameter
      av(a)     = zero
      ! Velocity divergence
      divv(a)   = zero
    ENDDO assign_module_variables2
    !$OMP END PARALLEL DO

    IF( debug ) PRINT *, "3"

    !
    !-- Compute nlrf and nstar from the ID (for diagnostics)
    !-- pressure and specific internal energy are already assigned
    !
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( this, sq_det_g4 ) &
    !$OMP             PRIVATE( a )
    compute_densities_from_id: DO a= 1, this% npart, 1
      this% nlrf (a)= this% baryon_density(a)
      this% nstar(a)= ( this% nlrf(a)*this% Theta(a) )*sq_det_g4(a)
    ENDDO compute_densities_from_id
    !$OMP END PARALLEL DO

    IF( debug ) PRINT *, "4"

    !-----------------------------------------------------------------------!
    ! Computation of the pressure.                                          !
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
            this% all_eos(i_matter)% eos_parameters(poly$kappa) &
            *( this% nlrf_sph(npart_in:npart_fin)*m0c2_cu ) &
            **this% all_eos(i_matter)% eos_parameters(poly$gamma)

          ! Using this internal energy gives machine-precision relative errors
          ! after the recovery, since it is computed from nlrf_sph
          ! Using the internal energy from the ID gives larger errors
          u(npart_in:npart_fin)= ( Pr(npart_in:npart_fin) &
            /(this% nlrf_sph(npart_in:npart_fin)*m0c2_cu &
            *( this% all_eos(i_matter)% eos_parameters(poly$gamma) - one ) ) )

          this% enthalpy(npart_in:npart_fin)= one + u(npart_in:npart_fin) &
           + Pr(npart_in:npart_fin)/(this% nlrf_sph(npart_in:npart_fin)*m0c2_cu)

          cs(npart_in:npart_fin)= SQRT( &
            this% all_eos(i_matter)% eos_parameters(poly$gamma) &
              *Pr(npart_in:npart_fin)/ &
            (this% nlrf_sph(npart_in:npart_fin)*m0c2_cu &
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
    !        CALL gen_pwp_cold_eos( this% nlrf_sph(a)*m0c2_cu, &
    !                               Pr(a), u(a), cs(a) )
    !
    !        !CALL gen_pwp_eos( this% nlrf_sph(a)*m0c2_cu, &
    !        !                  this% u_sph(a), tmp, &
    !        !                  u(a), &
    !        !                  Pr(a), cs(a) )
    !      ENDDO

          Pr(npart_in:npart_fin)= Pr(npart_in:npart_fin)/m0c2_cu
          this% pressure_sph(npart_in:npart_fin)= Pr(npart_in:npart_fin)
          this% u_sph(npart_in:npart_fin)= u(npart_in:npart_fin)

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
            this% all_eos(i_matter)% eos_parameters(poly$kappa) &
              *( this% nlrf_sph(a)*m0c2_cu ) &
              **this% all_eos(i_matter)% eos_parameters(poly$gamma) &
            + &
            ! thermal pressure
            Gamma_th_1*( this% nlrf_sph(a)*m0c2_cu )* &
              MAX(u(a) - ( Pr(a)/(this% nlrf_sph(a)*m0c2_cu &
                *( this% all_eos(i_matter)% eos_parameters(poly$gamma) &
                   - one ) ) ), zero)

          ENDDO
          this% enthalpy(npart_in:npart_fin)= one + u(npart_in:npart_fin) &
           + Pr(npart_in:npart_fin)/(this% nlrf_sph(npart_in:npart_fin)*m0c2_cu)

          cs(npart_in:npart_fin)= SQRT( &
            this% all_eos(i_matter)% eos_parameters(poly$gamma) &
              *Pr(npart_in:npart_fin)/ &
            (this% nlrf_sph(npart_in:npart_fin)*m0c2_cu &
            *this% enthalpy(npart_in:npart_fin)) )

          Pr(npart_in:npart_fin)= Pr(npart_in:npart_fin)/m0c2_cu
          this% pressure_sph(npart_in:npart_fin)= Pr(npart_in:npart_fin)
          this% u_sph(npart_in:npart_fin)= u(npart_in:npart_fin)

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

          CALL select_EOS_parameters( this% all_eos(i_matter)% eos_name )

          DO a= npart_in, npart_fin, 1

            CALL gen_pwp_cold_eos( this% nlrf_sph(a)*m0c2_cu, &
                                   Pr(a), u(a), cs(a) )

          ENDDO
          Pr(npart_in:npart_fin)= Pr(npart_in:npart_fin)/m0c2_cu
          this% pressure_sph(npart_in:npart_fin)= Pr(npart_in:npart_fin)
          this% u_sph(npart_in:npart_fin)= u(npart_in:npart_fin)

        ELSE
        ! If the system is hot, that is, has a thermal component, then
        ! the density and the specific energy (the latter including both
        ! cold and thermal part) should be supplied in the ID.
        ! The pressure is computed using them (see pwp_EOS MODULE).

          PRINT *, " * Assuming a hot system: thermal component considered."
          PRINT *

          u(npart_in:npart_fin)= this% specific_energy(npart_in:npart_fin)

          CALL select_EOS_parameters( this% all_eos(i_matter)% eos_name )

          DO a= npart_in, npart_fin, 1

            CALL gen_pwp_eos( this% nlrf_sph(a)*m0c2_cu, &
                              this% u_sph(a), tmp, &
                              u(a), &
                              Pr(a), cs(a) )

          ENDDO
          Pr(npart_in:npart_fin)= Pr(npart_in:npart_fin)/m0c2_cu
          this% pressure_sph(npart_in:npart_fin)= Pr(npart_in:npart_fin)
          this% u_sph(npart_in:npart_fin)= u(npart_in:npart_fin)

        ENDIF

      ENDIF

      END ASSOCIATE

    ENDDO matter_objects_loop

    !------------------!
    ! Assignment of Ye !
    !------------------!

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
               MAXVAL( this% nlrf_sph(npart_in:npart_fin), DIM= 1 ) &
               /((Msun_geo*km2m*m2cm)**3), " baryon cm^{-3}"
               !" m0c2_cu (TODO: CHECK UNITS)"!, &
               !"baryon Msun_geo^{-3}"
      PRINT *, " * Minimum interpolated nlrf on object", i_matter, "=", &
               MINVAL( this% nlrf_sph(npart_in:npart_fin), DIM= 1 ) &
               /((Msun_geo*km2m*m2cm)**3), " baryon cm^{-3}"
               !" m0c2_cu (TODO: CHECK UNITS)"!, &
               !"baryon Msun_geo^{-3}"
      PRINT *, " * Ratio between the two=", &
               MAXVAL( this% nlrf_sph(npart_in:npart_fin), DIM= 1 )/ &
               MINVAL( this% nlrf_sph(npart_in:npart_fin), DIM= 1 )
      PRINT *

      PRINT *, " * Maximum pressure", i_matter, "=", &
               MAXVAL( this% pressure_sph(npart_in:npart_fin), DIM= 1 ) &
               *amu*c_light2/((Msun_geo*km2m*m2cm)**3), " Ba"
               !" m0c2_cu (TODO: CHECK UNITS)"!, &
               !"baryon Msun_geo^{-3}"
      PRINT *, " * Minimum pressure", i_matter, "=", &
               MINVAL( this% pressure_sph(npart_in:npart_fin), DIM= 1 ) &
               *amu*c_light2/((Msun_geo*km2m*m2cm)**3), " Ba"
               !" m0c2_cu (TODO: CHECK UNITS)"!, &
               !"baryon Msun_geo^{-3}"
      PRINT *, " * Ratio between the two=", &
               MAXVAL( this% pressure_sph(npart_in:npart_fin), DIM= 1 )/ &
               MINVAL( this% pressure_sph(npart_in:npart_fin), DIM= 1 )
      PRINT *

      PRINT *, " * Maximum specific internal energy", i_matter, "=", &
               MAXVAL( this% u_sph(npart_in:npart_fin), DIM= 1 ), " c^2"
               !" m0c2_cu (TODO: CHECK UNITS)"!, &
               !"baryon Msun_geo^{-3}"
      PRINT *, " * Minimum specific internal energy", i_matter, "=", &
               MINVAL( this% u_sph(npart_in:npart_fin), DIM= 1 ), " c^2"
               !" m0c2_cu (TODO: CHECK UNITS)"!, &
               !"baryon Msun_geo^{-3}"
      PRINT *, " * Ratio between the two=", &
               MAXVAL( this% u_sph(npart_in:npart_fin), DIM= 1 )/ &
               MINVAL( this% u_sph(npart_in:npart_fin), DIM= 1 )
      PRINT *

      END ASSOCIATE

    ENDDO

    !
    !-- Printing the SPH ID to a binary file, for SPHINCS_BSSN
    !
    IF( this% export_bin )THEN

      IF( PRESENT(namefile) )THEN

        finalnamefile= TRIM(namefile)//"00000"
        dcount= -1 ! since it is increased before writing
        CALL write_SPHINCS_dump(finalnamefile)

      ELSE

        basename= "NSNS."
        dcount= -1 ! since it is increased before writing
        CALL write_SPHINCS_dump()

      ENDIF

    ENDIF


    !-------------!
    ! Diagnostics !
    !-------------!

    !
    !-- Test the recovery
    !
    CALL this% test_recovery( this% npart,       &
                              this% pos,         &
                              this% nlrf_sph,    &
                              this% u_sph,       &
                              this% pressure_sph, &
                              this% v(1:3,:),    &
                              this% theta,       &
                              this% nstar_sph )

    !CALL this% test_recovery( npart,       &
    !                          pos_u,         &
    !                          nlrf,    &
    !                          u,       &
    !                          Pr, &
    !                          vel_u,    &
    !                          Theta,       &
    !                          this% nstar_sph )


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
    !                             this% nlrf_sph   (npart_in:npart_fin),     &
    !                             this% u_sph      (npart_in:npart_fin),     &
    !                             this% pressure_sph(npart_in:npart_fin),     &
    !                             this% v          (1:3,npart_in:npart_fin), &
    !                             this% theta      (npart_in:npart_fin),     &
    !                             this% nstar_sph  (npart_in:npart_fin),     &
    !                             finalnamefile )
    !
    !   END ASSOCIATE
    !
    ! ENDDO

    !CALL exact_nei_tree_update( ndes,        &
    !                            this% npart, &
    !                            this% pos,   &
    !                            this% nu )
    !this% h= h

    CALL compute_and_print_quality_indicators &
      (this% npart, this% pos, this% h, this% nu, this% nstar_sph, &
       TRIM(sph_path))


    ALLOCATE( this% adm_linear_momentum_i( this% n_matter, 3 ) )
    this% adm_linear_momentum_fluid= zero
    DO i_matter= 1, this% n_matter, 1

      PRINT *, " * Estimating the ADM linear momentum using the canonical ", &
               "SPH momentum per baryon on the particles, ", &
               "on matter object ", i_matter, "..."
      PRINT *

      ASSOCIATE( npart_in   => this% npart_i(i_matter-1) + 1, &
                 npart_fin  => this% npart_i(i_matter-1) +    &
                               this% npart_i(i_matter) )

      CALL compute_adm_momentum_fluid_fields(                             &
                                  npart_fin - npart_in + 1,               &
                                  this% g_xx(npart_in:npart_fin),         &
                                  this% g_xy(npart_in:npart_fin),         &
                                  this% g_xz(npart_in:npart_fin),         &
                                  this% g_yy(npart_in:npart_fin),         &
                                  this% g_yz(npart_in:npart_fin),         &
                                  this% g_zz(npart_in:npart_fin),         &
                                  this% lapse(npart_in:npart_fin),        &
                                  this% shift_x(npart_in:npart_fin),      &
                                  this% shift_y(npart_in:npart_fin),      &
                                  this% shift_z(npart_in:npart_fin),      &
                                  this% nu(npart_in:npart_fin),           &
                                  this% Theta(npart_in:npart_fin),        &
                                  this% nlrf_sph(npart_in:npart_fin),     &
                                  this% pressure_sph(npart_in:npart_fin),  &
                                  this% u_sph(npart_in:npart_fin),        &
                                  this% v(1:3,npart_in:npart_fin),        &
                                  this% adm_linear_momentum_i(i_matter,:) )

      PRINT *, "   SPH estimate of the ADM linear momentum computed using ", &
               "the canonical momentum per baryon, on matter object", &
               i_matter,"= "
      PRINT *, "   (", this% adm_linear_momentum_i(i_matter, 1), ","
      PRINT *, "    ", this% adm_linear_momentum_i(i_matter, 2), ","
      PRINT *, "    ", this% adm_linear_momentum_i(i_matter, 3), ") Msun*c"
      PRINT *
      this% adm_linear_momentum_fluid= this% adm_linear_momentum_fluid + &
                                       this% adm_linear_momentum_i(i_matter,:)

      END ASSOCIATE

    ENDDO
    PRINT *, "   SPH estimate of the ADM momentum of the fluid ", &
             "computed using the canonical momentum per baryon= "
    PRINT *, "   (", this% adm_linear_momentum_fluid(1), ","
    PRINT *, "    ", this% adm_linear_momentum_fluid(2), ","
    PRINT *, "    ", this% adm_linear_momentum_fluid(3), ") Msun*c"
    PRINT *


    IF( ASSOCIATED(this% post_process_sph_id) )THEN

      CALL this% post_process_sph_id( this% npart, this% pos, &
                                    this% nlrf_sph, &
                                    this% u_sph, &
                                    this% pressure_sph, this% v(1:3,:), &
                                    this% theta, this% nstar_sph, this% nu, &
                                    this% g_xx,      &
                                    this% g_xy,      &
                                    this% g_xz,      &
                                    this% g_yy,      &
                                    this% g_yz,      &
                                    this% g_zz,      &
                                    this% lapse,     &
                                    this% shift_x,   &
                                    this% shift_y,   &
                                    this% shift_z,   &
                                    this% adm_linear_momentum_fluid, &
                                    this% adm_mass )

    ELSE

      PRINT *, "** ERROR! The PROCEDURE POINTER post_process_sph_id ", &
               "is not associated with any PROCEDURE!"
      PRINT *, " * Stopping..."
      PRINT *
      STOP

    ENDIF

    vel_u= this% v(1:3,:)

    this% adm_linear_momentum_fluid= zero
    DO i_matter= 1, this% n_matter, 1

      PRINT *, " * Estimating the ADM linear momentum using the canonical ", &
               "SPH momentum per baryon on the particles, ", &
               "on matter object ", i_matter, "..."
      PRINT *

      ASSOCIATE( npart_in   => this% npart_i(i_matter-1) + 1, &
                 npart_fin  => this% npart_i(i_matter-1) +    &
                               this% npart_i(i_matter) )

      CALL compute_adm_momentum_fluid_fields(                             &
                                  npart_fin - npart_in + 1,               &
                                  this% g_xx(npart_in:npart_fin),         &
                                  this% g_xy(npart_in:npart_fin),         &
                                  this% g_xz(npart_in:npart_fin),         &
                                  this% g_yy(npart_in:npart_fin),         &
                                  this% g_yz(npart_in:npart_fin),         &
                                  this% g_zz(npart_in:npart_fin),         &
                                  this% lapse(npart_in:npart_fin),        &
                                  this% shift_x(npart_in:npart_fin),      &
                                  this% shift_y(npart_in:npart_fin),      &
                                  this% shift_z(npart_in:npart_fin),      &
                                  this% nu(npart_in:npart_fin),           &
                                  this% Theta(npart_in:npart_fin),        &
                                  this% nlrf_sph(npart_in:npart_fin),     &
                                  this% pressure_sph(npart_in:npart_fin),  &
                                  this% u_sph(npart_in:npart_fin),        &
                                  vel_u(1:3,npart_in:npart_fin),        &
                                  this% adm_linear_momentum_i(i_matter,:) )

      PRINT *, "   SPH estimate of the ADM linear momentum computed using ", &
               "the canonical momentum per baryon, on matter object", &
               i_matter,"= "
      PRINT *, "   (", this% adm_linear_momentum_i(i_matter, 1), ","
      PRINT *, "    ", this% adm_linear_momentum_i(i_matter, 2), ","
      PRINT *, "    ", this% adm_linear_momentum_i(i_matter, 3), ") Msun*c"
      PRINT *
      this% adm_linear_momentum_fluid= this% adm_linear_momentum_fluid + &
                                       this% adm_linear_momentum_i(i_matter,:)

      END ASSOCIATE

    ENDDO
    PRINT *, "   SPH estimate of the ADM momentum of the fluid ", &
             "computed using the canonical momentum per baryon= "
    PRINT *, "   (", this% adm_linear_momentum_fluid(1), ","
    PRINT *, "    ", this% adm_linear_momentum_fluid(2), ","
    PRINT *, "    ", this% adm_linear_momentum_fluid(3), ") Msun*c"
    PRINT *

    !
    !-- Exporting the SPH ID to a binary file, for SPHINCS_BSSN
    !
 !   IF( this% export_bin )THEN
 !
 !     IF( PRESENT(namefile) )THEN
 !
 !       finalnamefile= TRIM( namefile ) // "00000"
 !       dcount= -1 ! since it is increased before writing
 !       CALL write_SPHINCS_dump( finalnamefile )
 !
 !     ELSE
 !
 !       basename= "NSNS."
 !       dcount= -1 ! since it is increased before writing
 !       CALL write_SPHINCS_dump()
 !
 !     ENDIF
 !
 !   ENDIF

    !
    !-- Compute particle number density
    !

    PRINT *, " * Computing particle number density by kernel interpolation..."
    PRINT *
    nu= one
    CALL density_loop( this% npart, this% pos, nu, h, &
                       this% particle_density_sph )

    IF( debug ) PRINT *, "100"

    !CALL COM( this% npart, this% pos, this% nu, &
    !          com_x_newt, com_y_newt, com_z_newt, com_d_newt )
    !CALL COM( this% npart_i(1), this% pos(:,1:this% npart_i(1)), &
    !          this% nu(1:this% npart_i(1)), &
    !          com_x_newt, com_y_newt, com_z_newt, com_d_newt )

    !CALL COM_1PN( this% npart, this% pos, this% v, &
    !              !this% v_euler_x, &
    !              this% nu, this% baryon_density, &
    !              this% specific_energy, this% nstar_sph, sq_detg4, gg4, &
    !              com_x_1pn, com_y_1pn, com_z_1pn, com_d_1pn )
    !CALL COM_1PN( this% npart_i(1), this% pos(:,1:this% npart_i(1)), &
    !              this% v(1:3,1:this% npart_i(1)), &
    !              !this% v_euler_x, &
    !              this% nu(1:this% npart_i(1)), &!/sq_det_g4(1:this% npart_i(1))/this% Theta(1:this% npart_i(1)), &
    !              !this% baryon_density(1:this% npart_i(1)), &
    !              this% nlrf_sph(1:this% npart_i(1)), &
    !              this% u_sph(1:this% npart_i(1)), &
    !              this% nstar_sph(1:this% npart_i(1)), &
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
    !                   this% specific_energy, this% pressure_sph, &
    !                   this% nstar_sph, this% h, &
    !                   sq_detg4, gg4, px, py, pz )
    !CALL momentum_1pn( this% npart, this% pos, &
    !                   this% v, &
    !                   !this% v_euler_x, &
    !                   this% nu, &
    !                   this% nlrf_sph, &
    !                   this% u_sph, &
    !                   this% pressure_sph, &
    !                   this% nstar_sph, &
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
!      CALL this% reshape_sph_field( this% pressure_sph, cnt1, cnt2, &
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

