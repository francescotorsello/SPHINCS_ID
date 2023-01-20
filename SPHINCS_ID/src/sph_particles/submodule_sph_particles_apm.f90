!# File:         submodule_sph_particles_apm.f90
!  Authors:      Francesco Torsello (FT)

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

SUBMODULE (sph_particles) apm

  !***********************************
  !
  !# This SUBMODULE contains the
  !  implementation of the method
  !  perform_apm of TYPE particles.
  !
  !  FT 04.06.2021
  !
  !***********************************

  USE constants,  ONLY: quarter
  USE utility,    ONLY: zero, one, two, three, four, five, ten, sph_path, &
                        Msun_geo


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE perform_apm

    !*****************************************************
    !
    !# Compute the particle positions as follows:
    !
    !   1. Take initial particle distribution as input
    !   2. Assume that the particles have the same mass
    !   3. Do the APM iteration so that the final
    !      SPH kernel estimate of the baryon mass
    !      density matches the baryon density in the
    !      given |id|
    !   4. Correct the particle masses ONCE, in order
    !      to match the density even better. Since we
    !      don't want a large mass ratio, we impose a
    !      maximum mass ratio when performing this
    !      correction.
    !
    !  After this procedure, the resulting particle
    !  distribution has positions and baryon numbers
    !  that kernel-estimate very well the mass density
    !  of the star, and has a low mass ratio.
    !
    !  This procedure assigns positions, smoothing
    !  lengths \(h\), and \(\nu\).
    !
    !  @warning
    !  If the outer layers of a star have a very low density
    !  compared to the core, it can happen that, irrespective
    !  of the initial particle distribution and the APM
    !  parameters, the particle distribution output by the
    !  APM does not have a smooth surface. In this case,
    !  the only solution (that has been found as of 20.10.2021)
    !  is to increase the particle number.
    !
    !  As of 20.10.2021, this has only happened with the
    !  CompOSE tabulated EOS (all of them), but not with
    !  any piecewise polytropic or polytropic EOS.
    !
    !  This problem can manifest itself with error messages
    !  concerning particles occupying the sameposition,
    !  or some smoothing lengths being 0, or a matrix not being
    !  invertible.
    !
    !  FT 20.10.2021
    !  @endwarning
    !
    !  FT 04.06.2021
    !
    !*****************************************************

    USE utility,             ONLY: cnt, spherical_from_cartesian, &
                                   cartesian_from_spherical
    USE constants,           ONLY: half, third, Msun, amu, pi

    USE sph_variables,       ONLY: allocate_sph_memory, deallocate_sph_memory, &
                                   npart, h, nu
    USE metric_on_particles, ONLY: allocate_metric_on_particles, &
                                   deallocate_metric_on_particles
    USE gradient,            ONLY: allocate_gradient, deallocate_gradient
    USE set_h,               ONLY: exact_nei_tree_update, posmash
    USE units,               ONLY: umass, m0c2_cu

    USE APM,                 ONLY: density_loop, position_correction, assign_h
    USE analyze,             ONLY: COM
    USE matrix,              ONLY: determinant_4x4_matrix

    USE sphincs_sph,         ONLY: density, ncand, all_clists
    USE matrix,              ONLY: invert_3x3_matrix
    USE kernel_table,        ONLY: &!dWdv_no_norm, &
                                   dv_table, dv_table_1,&
                                   W_no_norm, n_tab_entry, &
                                   interp_gradW_table,interp_W_gradW_table
    USE RCB_tree_3D,         ONLY: iorig, &
                                   !nic, nfinal, nprev, lpart, rpart, &
                                   allocate_RCB_tree_memory_3D, &
                                   deallocate_RCB_tree_memory_3D

    IMPLICIT NONE

    INTEGER,          PARAMETER:: max_npart        = 20D+6
    INTEGER,          PARAMETER:: nn_des           = 301
    INTEGER,          PARAMETER:: m_max_it         = 50
    INTEGER,          PARAMETER:: search_pos       = 10
    !INTEGER,          PARAMETER:: print_step       = 15
    INTEGER,          PARAMETER:: nuratio_max_steps= 200
    INTEGER,          PARAMETER:: nuratio_min_it   = 100

    DOUBLE PRECISION, PARAMETER:: tol               = 1.0D-3
    !DOUBLE PRECISION, PARAMETER:: iter_tol          = 2.0D-2
    !DOUBLE PRECISION, PARAMETER:: backup_h          = 0.25D0
    DOUBLE PRECISION, PARAMETER:: max_art_pr_ghost  = 1.0D+10
    DOUBLE PRECISION, PARAMETER:: tiny_real         = 1.0D-10
    DOUBLE PRECISION, PARAMETER:: nuratio_tol       = 2.5D-3

    INTEGER:: a, itr, itr2, i_shell, n_inc, cnt1, b, inde, index1   ! iterators
    INTEGER:: npart_real, npart_real_half, npart_ghost, npart_all
    INTEGER:: nx, ny, nz, i, j, k
    INTEGER:: a_numin, a_numin2, a_numax, a_numax2
    INTEGER:: nuratio_cnt
    INTEGER:: dim_seed, rel_sign
    INTEGER:: n_problematic_h, ill, l, itot, cnt_push_ghost, max_push_ghost
    INTEGER, DIMENSION(:), ALLOCATABLE:: cnt_move

    DOUBLE PRECISION:: ghost_displacement, min_radius, max_radius, radius_part
    DOUBLE PRECISION:: ellipse_thickness
    DOUBLE PRECISION:: radius_x, radius_y, radius_z
    DOUBLE PRECISION:: h_max, h_av, tmp, dens_min, atmosphere_density!, delta
    DOUBLE PRECISION:: xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, &
                       rad_x, rad_y, rad_z, com_x, com_y, com_z, com_d
    DOUBLE PRECISION:: max_z_real
    DOUBLE PRECISION:: xtemp, ytemp, ztemp, x_ell, y_ell, z_ell
    DOUBLE PRECISION:: min_nu, max_nu, min_nu2, max_nu2
    ! The value of nu equal for all the particles, used during the APM iteration
    DOUBLE PRECISION:: nu_all, nu_ghost
    DOUBLE PRECISION:: err_N_mean_min, err_N_mean_min_old, err_N_mean, &
                       err_mean_old, err_n_min, err_N_max, dN, &!dNstar, &
                       nstar_id_err, nstar_sph_err, dN_max, dN_av
    DOUBLE PRECISION:: art_pr_max
    DOUBLE PRECISION:: nu_tot, nu_ratio, nu_tmp2, nuratio_tmp, nuratio_tmp_prev
    DOUBLE PRECISION:: variance_nu, stddev_nu, mean_nu
    DOUBLE PRECISION:: variance_dN, stddev_dN
    DOUBLE PRECISION:: rand_num, rand_num2
    DOUBLE PRECISION:: r, theta, phi
    DOUBLE PRECISION:: r_ell, theta_ell, phi_ell
    DOUBLE PRECISION:: dS_norm_av
    DOUBLE PRECISION:: nstar_id_av, nstar_sph_av, nlrf_id_av, pressure_id_av
    DOUBLE PRECISION:: l2norm_displacement, average_displacement

    INTEGER, DIMENSION(:), ALLOCATABLE:: neighbors_lists
    INTEGER, DIMENSION(:), ALLOCATABLE:: n_neighbors
    INTEGER, DIMENSION(:), ALLOCATABLE:: seed

    DOUBLE PRECISION, DIMENSION(3):: pos_corr_tmp
    DOUBLE PRECISION, DIMENSION(3):: pos_maxerr
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos_tmp
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: ghost_pos
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: all_pos
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: all_pos_prev
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: correction_pos

    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: h_guess
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: h_tmp
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: h_guess_tmp
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: cnt_array

    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: rho_tmp
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar_id
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nlrf_id
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_id
    !DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar_eul_id
    !DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nu_eul
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar_sph
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nlrf_sph
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_sph
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: sqg
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: dS
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: dNstar
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: art_pr

    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nu_tmp
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pvol_tmp
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nu_one

    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar_int
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: particle_density_final

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: nearest_neighbors

    LOGICAL:: exist

    CHARACTER(LEN=:), ALLOCATABLE:: finalnamefile

    LOGICAL, PARAMETER:: debug= .FALSE.
    !LOGICAL:: few_ncand!, invertible_matrix
    LOGICAL:: push_away_ghosts

    TYPE(timer):: find_h_bruteforce_timer

    find_h_bruteforce_timer= timer( "find_h_bruteforce_timer" )

    CALL RANDOM_SEED( SIZE= dim_seed )
    ALLOCATE( seed( dim_seed ) )
    seed( 1 )= 2
    seed( 2 )= 1
    DO itr= 3, dim_seed
      seed( itr )= seed( itr - 1 ) + seed( itr - 2 )
    ENDDO
    CALL RANDOM_SEED( PUT= seed )

    IF( debug ) PRINT *, "0"

    npart_real= SIZE( pos_input(1,:) )

    IF( debug ) PRINT *, "npart_real= ", npart_real

    PRINT *, "** Steering parameters for the APM iteration:"
    PRINT *
    PRINT *, "   mass_it= ",           mass_it
    PRINT *, "   use_pressure= ",      use_pressure
    PRINT *, "   adapt_ghosts= ",      adapt_ghosts
    PRINT *, "   move_away_ghosts= ",  move_away_ghosts
    PRINT *, "   use_atmosphere= ",    use_atmosphere
    PRINT *, "   remove_atmosphere= ", remove_atmosphere
    PRINT *

    IF(use_pressure)THEN
      PRINT *, "** Using the physical pressure rather than the density to ", &
               "compute the artificial pressure."
      PRINT *
    ENDIF

    !---------------------------------------------------------------!
    !-- Tune the displacement the ghosts depending on the desired --!
    !-- EXIT condition for the APM iteration                      --!
    !---------------------------------------------------------------!

    IF( nuratio_des > zero )THEN

      ghost_displacement= third/DBLE(nuratio_max_steps)/ten
      max_push_ghost    = two*nuratio_max_steps

    ELSE

      ghost_displacement= third/DBLE(max_inc)/ten
      max_push_ghost    = max_inc

    ENDIF

    !------------------------------------------------!
    !-- If desired, compute the atmosphere density --!
    !------------------------------------------------!

    IF( use_atmosphere )THEN

      dens_min= HUGE(one)
      DO a= 1, npart_real, 1

        tmp= get_density( pos_input(1,a), pos_input(2,a), pos_input(3,a) )

        IF( tmp < dens_min )THEN

          dens_min= tmp

        ENDIF

      ENDDO
      atmosphere_density= zero*dens_min*1.0D-30

    ENDIF

    !---------------------------------------!
    !-- Allocate, assign and test h_guess --!
    !---------------------------------------!

    IF(.NOT.ALLOCATED( h_guess ))THEN
      ALLOCATE( h_guess( max_npart ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array h_guess in SUBROUTINE ", &
                  "perform_apm. The error message is",&
                  err_msg
         STOP
      ENDIF
    ENDIF

    h_guess= zero
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( h_guess, pvol, npart_real ) &
    !$OMP             PRIVATE( a )
    DO a= 1, npart_real, 1
      h_guess(a)= three*(pvol(a)**third)
    ENDDO
    !$OMP END PARALLEL DO

    find_h_bruteforce_timer= timer( "find_h_bruteforce_timer" )
    CALL find_h_bruteforce_timer% start_timer()
    n_problematic_h= 0
    check_h_guess: DO a= 1, npart_real, 1
    ! find_h_backup, called below, is OMP parallelized, so this loop
    ! should not be parallelized as well

      IF( .NOT.is_finite_number( h_guess(a) ) .OR. h_guess(a) <= zero )THEN

        n_problematic_h= n_problematic_h + 1
        h_guess(a)= find_h_backup( a, npart_real, pos_input, nn_des )
        IF( .NOT.is_finite_number( h_guess(a) ) .OR. h_guess(a) <= zero )THEN
          PRINT *, "** ERROR! h=0 on particle ", a, "even with the brute", &
                   " force method."
          PRINT *, "   Particle position: ", pos_input(:,a)
          STOP
        ENDIF

      ENDIF

    ENDDO check_h_guess
    CALL find_h_bruteforce_timer% stop_timer()
    CALL find_h_bruteforce_timer% print_timer( 2 )

    PRINT *, " * The smoothing length was found brute-force for ", &
             n_problematic_h, " particles."
    PRINT *

    IF( debug ) PRINT *, "0.5"

    !--------------------------------------------------------------------!
    !-- Store particles above xy plane as the first half of the array, --!
    !-- and mirror them to the second half                             --!
    !--------------------------------------------------------------------!

    pos_tmp= pos_input
    h_tmp= h_guess
    itr= 0

    DO a= 1, npart_real, 1
      IF( pos_tmp(3,a) > zero )THEN
        itr= itr + 1
        pos_input(1,itr)= pos_tmp(1,a)
        pos_input(2,itr)= pos_tmp(2,a)
        pos_input(3,itr)= pos_tmp(3,a)
        h_guess( itr )     = h_tmp( itr )
      ENDIF
    ENDDO
    npart_real_half= itr

    DO a= 1, npart_real_half, 1
      pos_input(1, npart_real_half + a)=   pos_input(1,a)
      pos_input(2, npart_real_half + a)=   pos_input(2,a)
      pos_input(3, npart_real_half + a)= - pos_input(3,a)
      h_guess( npart_real_half + a )     =   h_guess( a )
    ENDDO
    npart_real= 2*npart_real_half

    IF( debug ) PRINT *, "1"
    IF( debug ) PRINT *, "npart_real= ", npart_real

    !--------------------------------------------------------------------!
    !-- Find the maximum and the average smoothing length of the       --!
    !-- particles whose distance from the center is higher than        --!
    !-- radius_z, and use them to place ghost particles a little more  --!
    !-- outside than the surface of the particles.                     --!
    !--------------------------------------------------------------------!

    radius_x= MAX( sizes(1), sizes(2) )
    radius_y= MAX( sizes(3), sizes(4) )
    radius_z= MAX( sizes(5), sizes(6) )
    min_radius = MINVAL([radius_x, radius_y, radius_z])
    max_radius = MAXVAL([radius_x, radius_y, radius_z])

    h_max= zero
    h_av = zero
    itr  = 0
    max_z_real= ABS( MAXVAL( ABS(pos_input(3,:) - center(3)), DIM= 1 ) )
    ALLOCATE(cnt_array(npart_real))
    cnt_array= zero
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( pos_input, npart_real, center, h_guess, &
    !$OMP                     cnt_array, max_z_real ) &
    !$OMP             PRIVATE( a ) &
    !$OMP             REDUCTION( MAX: h_max ) &
    !$OMP             REDUCTION( +: h_av )
    DO a= 1, npart_real, 1

      IF( SQRT( ( pos_input(1,a) - center(1) )**two &
              + ( pos_input(2,a) - center(2) )**two &
              + ( pos_input(3,a) - center(3) )**two ) &
                > (one - one/(ten*ten))*max_z_real )THEN

        cnt_array(a)= one
        !IF( h_guess(a) > h_max )THEN
        !  h_max= h_guess(a)
        !ENDIF
        h_max= MAX(h_max, h_guess(a))
        h_av = h_av + h_guess(a)

      ENDIF

    ENDDO
    !$OMP END PARALLEL DO
    h_av= h_av/SUM(cnt_array, DIM=1)
    IF( debug ) PRINT *, "h_av=", h_av
    IF( debug ) PRINT *

    IF( debug ) PRINT *, "2"

    !-------------------------------!
    !--  Placing ghost particles  --!
    !-------------------------------!

    CALL place_and_print_ghost_particles()

    ! Assign values to all_pos
    IF(.NOT.ALLOCATED( all_pos ))THEN
      ALLOCATE( all_pos( 3, npart_all ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array ghost_pos in SUBROUTINE ", &
                  "perform_apm. The error message is",&
                  err_msg
         STOP
      ENDIF
    ENDIF
    all_pos(:, 1:npart_real)          = pos_input
    all_pos(:, npart_real+1:npart_all)= ghost_pos

    CALL allocate_apm_fields( npart_real, npart_ghost )

    IF(debug)THEN
      !
      !-- Estimate ghost density considering also the real particles
      !

      npart= npart_all

      CALL allocate_SPH_memory

      CALL allocate_RCB_tree_memory_3D(npart)
      iorig(1:npart)= (/ (a,a=1,npart) /)

      IF( debug ) PRINT *, "10"

      CALL allocate_gradient( npart )
      CALL allocate_metric_on_particles( npart )

      ! Determine smoothing length so that each particle has exactly
      ! 300 neighbours inside 2h
      CALL assign_h( nn_des, &
                     npart_all, &
                     all_pos, &
                     h_guess, &
                     h )

      find_h_bruteforce_timer= timer( "find_h_bruteforce_timer" )
      CALL find_h_bruteforce_timer% start_timer()
      n_problematic_h= 0
      DO a= 1, npart_all, 1
      ! find_h_backup, called below, is OMP parallelized, so this loop
      ! should not be parallelized as well

        IF( .NOT.is_finite_number( h(a) ) .OR. h(a) <= zero )THEN

          n_problematic_h= n_problematic_h + 1
          h(a)= find_h_backup( a, npart_all, all_pos(:,:), nn_des)
          IF( .NOT.is_finite_number( h(a) ) .OR. h(a) <= zero )THEN
            PRINT *, "** ERROR! h=0 on particle ", a, "even with the brute", &
                     " force method."
            PRINT *, "   Particle position: ", all_pos(:,a)
            STOP
          ENDIF

        ENDIF

      ENDDO
      CALL find_h_bruteforce_timer% stop_timer()
      CALL find_h_bruteforce_timer% print_timer( 2 )

      PRINT *, " * The smoothing length was found brute-force for ", &
               n_problematic_h, " particles."
      PRINT *

      PRINT *, " * Measure SPH particle number density..."
      PRINT *

      nu(1:npart_real)          = (mass/DBLE(npart_real))*umass/amu
      nu(npart_real+1:npart_all)= nu_ghost
      PRINT *, "npart_all         =", npart_all
      PRINT *, "SIZE(all_pos(1,:))=", SIZE(all_pos(1,:))
      PRINT *, "SIZE(nu)          =", SIZE(nu)
      PRINT *, "SIZE(h)           =", SIZE(h)
      PRINT *, "SIZE(nstar_sph)   =", SIZE(nstar_sph)
      PRINT *
      CALL density_loop( npart_all, all_pos, &
                         nu, h, nstar_sph )

      !nstar_sph_ghost_av= SUM(nstar_sph_ghost)/npart_ghost
      PRINT *, "nstar_sph_ghost_max=", &
               MAXVAL(nstar_sph(npart_real+1:npart_all), DIM=1)
      PRINT *, "nstar_sph_ghost_av=", &
               SUM(nstar_sph(npart_real+1:npart_all))/npart_ghost
      PRINT *, "nstar_sph_ghost_min=", &
               MINVAL(nstar_sph(npart_real+1:npart_all), DIM=1)

      CALL deallocate_metric_on_particles()
      CALL deallocate_gradient()
      CALL deallocate_RCB_tree_memory_3D()
      CALL deallocate_SPH_memory()

    ENDIF

    !STOP

    !----------------------------!
    !-- Allocate needed memory --!
    !----------------------------!

    npart= npart_all

    CALL allocate_SPH_memory

    CALL allocate_RCB_tree_memory_3D(npart)
    iorig(1:npart)= (/ (a,a=1,npart) /)

    IF( debug ) PRINT *, "10"

    CALL allocate_gradient( npart )
    CALL allocate_metric_on_particles( npart )

    !-------------------------------------!
    !-- Setting up ID for APM iteration --!
    !-------------------------------------!

    PRINT *, "** Setting up ID for APM iteration..."
    PRINT *

    PRINT *, " * Assigning h..."
    PRINT *

    ! Determine smoothing length so that each particle has exactly
    ! 300 neighbours inside 2h
    CALL assign_h( nn_des, &
                   npart_all, &
                   all_pos, h_guess, &
                   h )

    find_h_bruteforce_timer= timer( "find_h_bruteforce_timer" )
    CALL find_h_bruteforce_timer% start_timer()
    n_problematic_h= 0
    check_h1: DO a= 1, npart_real, 1
    ! find_h_backup, called below, is OMP parallelized, so this loop
    ! should not be parallelized as well

      IF( .NOT.is_finite_number( h(a) ) .OR. h(a) <= zero )THEN

        n_problematic_h= n_problematic_h + 1
        h(a)= find_h_backup( a, npart_real, all_pos(:,1:npart_real), nn_des )
        IF( .NOT.is_finite_number( h(a) ) .OR. h(a) <= zero )THEN
          PRINT *, "** ERROR! h=0 on particle ", a, "even with the brute", &
                   " force method."
          PRINT *, "   Particle position: ", all_pos(:,a)
          STOP
        ENDIF

      ENDIF

    ENDDO check_h1
    CALL find_h_bruteforce_timer% stop_timer()
    CALL find_h_bruteforce_timer% print_timer( 2 )

    PRINT *, " * The smoothing length was found brute-force for ", &
             n_problematic_h, " particles."
    PRINT *

    PRINT *, " * Measure SPH particle number density..."
    PRINT *

    nu= one
    CALL density_loop( npart_all, all_pos, &    ! input
                       nu, h, nstar_sph )      ! output

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( all_pos, npart_all, nstar_sph, h, nu, &
    !$OMP                     center ) &
    !$OMP             PRIVATE( a )
    check_nstar_sph: DO a= 1, npart_all, 1

      IF( .NOT.is_finite_number( nstar_sph( a ) ) )THEN

        PRINT *, "** WARNING! nstar_sph(", a, ") is a not a finite number ",&
                 "in SUBROUTINE perform_apm!"
        IF( debug ) PRINT *, " * h(", a, ")=", h(a)
        IF( debug ) PRINT *, " * nu(", a, ")=", nu(a)
        IF( debug ) PRINT *, " * all_pos(", a, ")=", all_pos(:,a)
        IF( debug ) PRINT *, " * r(", a, ")=", &
                              SQRT( ( all_pos(1,a) - center )**two &
                              + all_pos(2,a)**two + all_pos(3,a)**two )
        PRINT *
        STOP

      ENDIF

    ENDDO check_nstar_sph
    !$OMP END PARALLEL DO

    IF( debug ) PRINT *, "4"

    max_nu= zero
    min_nu= HUGE(one)

    IF( debug ) PRINT *, "7"

    CALL get_nstar_id_atm( npart_real, all_pos(1,1:npart_real), &
                           all_pos(2,1:npart_real), &
                           all_pos(3,1:npart_real), &
                           nstar_sph, nstar_id, nlrf_sph, sqg, &
                           use_atmosphere )

  ! The following test is done inside get_nstar_id_atm. Kept here for paranoia
  !  !$OMP PARALLEL DO DEFAULT( NONE ) &
  !  !$OMP             SHARED( all_pos, npart_all, nstar_id, h, nu, &
  !  !$OMP                     center, dNstar ) &
  !  !$OMP             PRIVATE( a )
  !  check_nstar_id: DO a= 1, npart_all, 1
  !
  !    IF( .NOT.is_finite_number( nstar_id( a ) ) )THEN
  !
  !      PRINT *, "** WARNING! nstar_id(", a, ") is a not a finite number ", &
  !               "in SUBROUTINE perform_apm!"
  !      PRINT *, "   nstar_id(", a, ")=", nstar_id(a)
  !      PRINT *, "   dNstar(", a, ")=", dNstar(a)
  !      PRINT *, "   rho(", a, ")=", get_density( all_pos(1,a), &
  !                                                all_pos(2,a), &
  !                                                all_pos(3,a) )
  !      IF( debug ) PRINT *, " * h(", a, ")=", h(a)
  !      IF( debug ) PRINT *, " * nu(", a, ")=", nu(a)
  !      IF( debug ) PRINT *, " * all_pos(", a, ")=", all_pos(:,a)
  !      IF( debug ) PRINT *, " * r(", a, ")=", &
  !                            SQRT( ( all_pos(1,a) - center )**two &
  !                            + all_pos(2,a)**two + all_pos(3,a)**two )
  !      PRINT *
  !      STOP
  !
  !    ENDIF
  !
  !  ENDDO check_nstar_id
  !  !$OMP END PARALLEL DO

    IF( debug ) PRINT *, "8"

    !----------------------------------------------------!
    !-- enforce centre of mass after having changed nu --!
    !----------------------------------------------------!

    IF( com_star(1) == zero &
        .AND. com_star(2) == zero .AND. com_star(3) == zero )THEN

      CALL COM( npart_real, all_pos(:,1:npart_real), nu(1:npart_real), & ! input
                com_star(1), com_star(2), com_star(3), com_d ) ! output

    ENDIF

    CALL correct_center_of_mass( npart_real, all_pos(:,1:npart_real), &
                                 nu(1:npart_real), get_density, &
                                 validate_position_final, com_star, &
                                 verbose= .TRUE. )

    !-----------------------------------------------------------------------!
    !-- Mirror the positions after having repositioned the center of mass --!
    !-----------------------------------------------------------------------!

    CALL impose_equatorial_plane_symmetry( npart_real, &
                                           all_pos(:,1:npart_real), &
                                           nu(1:npart_real) )

    PRINT *, " * ID set up for the APM iteration."
    PRINT *

    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !--                           APM iteration                           --!
    !-- Assume equal mass particles, and move them so that the SPH kernel --!
    !-- estimate of the mass density matches the star mass density as     --!
    !-- well as reasonably possible.                                      --!
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!



    PRINT *, " * Performing APM iteration..."
    PRINT *

    cnt_move= 0

    ! Set the particles to be equal-mass
    nu_all= (mass/DBLE(npart_real))*umass/amu
    nu(1:npart_real)          = nu_all
    nu(npart_real+1:npart_all)= nu_ghost
    ! nu_ghost is set in place_and_print_ghost_particles

    CALL correct_center_of_mass( npart_real, all_pos(:,1:npart_real), &
                                 nu(1:npart_real), get_density, &
                                 validate_position_final, com_star, &
                                 verbose= .TRUE. )

    all_pos_prev= -one
    PRINT *, " * The APM iteration starts here."
    PRINT *

    n_inc           = 0
    err_N_mean_min  = HUGE(one)
    nuratio_cnt     = 0
    nuratio_tmp_prev= 1.D-8
    cnt_push_ghost  = 0
    apm_iteration: DO itr= 1, apm_max_it, 1

      PRINT *, "------------------------------------------"
      PRINT *, " * Starting with APM step #: ", itr
      PRINT *

      IF( print_step /= 0 &
          .AND. &
          MOD( itr, print_step ) == 0 )THEN

     ! DEBUGGING
     !
     !   DO a= 1, npart_real, 1
     !     IF( check_particle_position( a - 1, &
     !                                  all_pos(:,1:a-1), &
     !                                  all_pos(:,a) ) > 0 &
     !         .AND. &
     !         check_particle_position( npart_real - a, &
     !                                  all_pos(:,a+1:npart_real), &
     !                                  all_pos(:,a) ) > 0 &
     !     )THEN
     !
     !       CALL RANDOM_NUMBER( rand_num )
     !       CALL RANDOM_NUMBER( rand_num2 )
     !
     !       IF( rand_num2 < half )  rel_sign= - 1
     !       IF( rand_num2 >= half ) rel_sign=   1
     !
     !       all_pos(:,a)= all_pos(:,a)*( one + &
     !                                    DBLE(rel_sign)*rand_num*half*third )
     !
     !     ENDIF
     !   ENDDO
     !
     ! END DEBUGGING

        PRINT *, " * Printing temporary APM data to file..."
        CALL dump_apm_pos()
        PRINT *, "   ...done."
        PRINT *

      ENDIF

      IF( debug ) PRINT *, "enforcing center of mass..."

      CALL correct_center_of_mass( npart_real, all_pos(:,1:npart_real), &
                                   nu(1:npart_real), get_density, &
                                   validate_position_final, com_star )

      IF( debug ) PRINT *, "mirroring particles..."

      CALL impose_equatorial_plane_symmetry( npart_real, &
                                             all_pos(:,1:npart_real), &
                                             nu(1:npart_real) )

      IF( debug )THEN

        CALL COM( npart_real, all_pos(:,1:npart_real), &  !
                  nu(1:npart_real), &                     ! input
                  com_x, com_y, com_z, com_d )            ! output

        PRINT *, "** After center of mass correction:"
        PRINT *, " * x coordinate of the center of mass of the star, ", &
                 "from ID: com_star= ", com_star, "Msun_geo"
        PRINT *, " * x coordinate of the center of mass of the particle ", &
                 "distribution: com_x= ", com_x, "Msun_geo"
        PRINT *, " * y coordinate of the center of mass of the particle ", &
                 "distribution: com_y= ", com_y, "Msun_geo"
        PRINT *, " * z coordinate of the center of mass of the particle ", &
                 "distribution: com_z= ", com_z, "Msun_geo"
        PRINT *, " * Distance of the center of mass of the particle ", &
                 "distribution from the  origin: com_d= ", com_d
        PRINT *, " * |com_x-com_star(1)|/|com_star(1)+1|=", &
                 ABS( com_x-com_star(1) )/ABS( com_star(1) + 1 )
        PRINT *, " * |com_y-com_star(2)|/|com_star(2)+1|=", &
                 ABS( com_y-com_star(2) )/ABS( com_star(2) + 1 )
        PRINT *, " * |com_z-com_star(3)|/|com_star(3)+1|=", &
                 ABS( com_z-com_star(3) )/ABS( com_star(3) + 1 )
        PRINT *

        finalnamefile= TRIM(sph_path)//"dbg-pos.dat"

        INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

        IF( exist )THEN
            OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
                  FORM= "FORMATTED", &
                  POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
                  IOMSG= err_msg )
        ELSE
            OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "NEW", &
                  FORM= "FORMATTED", &
                  ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
        ENDIF
        IF( ios > 0 )THEN
          PRINT *, "...error when opening " // TRIM(finalnamefile), &
                   ". The error message is", err_msg
          STOP
        ENDIF

        DO a= 1, npart_real, 1
          WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
            1, a, &
            all_pos(1,a), &
            all_pos(2,a), &
            all_pos(3,a)
        ENDDO

        DO a= npart_real+1, npart_all, 1
          WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
            2, a, &
            all_pos(1,a), &
            all_pos(2,a), &
            all_pos(3,a)
        ENDDO

        CLOSE( UNIT= 2 )

      ENDIF

      IF( debug ) PRINT *, "assign h..."

      h_guess(1:npart_real)= h(1:npart_real)
      !h_guess(npart_real+1:npart_all)= dx*dy*dz

      CALL assign_h( nn_des, &
                     npart_all, &
                     all_pos, h_guess, h )

      find_h_bruteforce_timer= timer( "find_h_bruteforce_timer" )
      CALL find_h_bruteforce_timer% start_timer()
      n_problematic_h= 0
      check_h2: DO a= 1, npart_real, 1
      ! find_h_backup, called below, is OMP parallelized, so this loop
      ! should not be parallelized as well

        IF( .NOT.is_finite_number( h(a) ) .OR. h(a) <= zero )THEN

          n_problematic_h= n_problematic_h + 1
          h(a)= find_h_backup( a, npart_real, all_pos(:,1:npart_real), nn_des )
          IF( .NOT.is_finite_number( h(a) ) .OR. h(a) <= zero )THEN
            PRINT *, "** ERROR! h=0 on particle ", a, "even with the brute", &
                     " force method."
            PRINT *, "   Particle position: ", all_pos(:,a)
            STOP
          ENDIF

        ENDIF

      ENDDO check_h2
      CALL find_h_bruteforce_timer% stop_timer()
      CALL find_h_bruteforce_timer% print_timer( 2 )

      PRINT *, " * The smoothing length was found brute-force for ", &
               n_problematic_h, " particles."
      PRINT *

      IF( debug ) PRINT *, "density_loop..."

      CALL density_loop( npart_all, all_pos, &    ! input
                         nu, h, nstar_sph )      ! output

      IF( debug ) PRINT *, "npart_real= ", npart_real
      IF( debug ) PRINT *, "npart_all= ", npart_all
      IF( debug ) PRINT *

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( all_pos, npart_all, nstar_sph, h, nu, &
      !$OMP                     center ) &
      !$OMP             PRIVATE( a )
      check_nstar_sph2: DO a= 1, npart_all, 1

        IF( .NOT.is_finite_number( nstar_sph( a ) ) )THEN

          PRINT *, "** WARNING! nstar_sph(", a, ") is a not a finite number ",&
                   "in SUBROUTINE perform_apm!"
          IF( debug ) PRINT *, " * h(", a, ")=", h(a)
          IF( debug ) PRINT *, " * nu(", a, ")=", nu(a)
          IF( debug ) PRINT *, " * all_pos(", a, ")=", all_pos(:,a)
          IF( debug ) PRINT *, " * r(", a, ")=", &
                                SQRT( ( all_pos(1,a) - center )**two &
                                + all_pos(2,a)**two + all_pos(3,a)**two )
          PRINT *
          STOP

        ENDIF

      ENDDO check_nstar_sph2
      !$OMP END PARALLEL DO

      CALL get_nstar_id_atm( npart_real, all_pos(1,1:npart_real), &
                             all_pos(2,1:npart_real), &
                             all_pos(3,1:npart_real), &
                             nstar_sph, nstar_id, nlrf_sph, sqg, &
                             use_atmosphere )

! The following test is done inside get_nstar_id_atm. Kept here for paranoia
!      !$OMP PARALLEL DO DEFAULT( NONE ) &
!      !$OMP             SHARED( all_pos, npart_all, nstar_id, h, nu, &
!      !$OMP                     center, dNstar ) &
!      !$OMP             PRIVATE( a )
!      check_nstar_id2: DO a= 1, npart_all, 1
!
!        IF( .NOT.is_finite_number( nstar_id( a ) ) )THEN
!
!          PRINT *, "** WARNING! nstar_id(", a, ") is a not a finite number ", &
!                   "in SUBROUTINE perform_apm!"
!          PRINT *, "   nstar_id(", a, ")=", nstar_id(a)
!          PRINT *, "   dNstar(", a, ")=", dNstar(a)
!          PRINT *, "   rho(", a, ")=", get_density( all_pos(1,a), &
!                                                    all_pos(2,a), &
!                                                    all_pos(3,a) )
!          IF( debug ) PRINT *, " * h(", a, ")=", h(a)
!          IF( debug ) PRINT *, " * nu(", a, ")=", nu(a)
!          IF( debug ) PRINT *, " * all_pos(", a, ")=", all_pos(:,a)
!          IF( debug ) PRINT *, " * r(", a, ")=", &
!                                SQRT( ( all_pos(1,a) - center )**two &
!                                + all_pos(2,a)**two + all_pos(3,a)**two )
!          PRINT *
!          STOP
!
!        ENDIF
!
!      ENDDO check_nstar_id2
!      !$OMP END PARALLEL DO

      IF(use_pressure)THEN

        CALL read_pressure_id( npart_real, &
                               all_pos(1,1:npart_real), &
                               all_pos(2,1:npart_real), &
                               all_pos(3,1:npart_real), &
                               pressure_id(1:npart_real) )

        CALL compute_pressure( npart_all, &
                               all_pos(1,:), all_pos(2,:), all_pos(3,:), &
                               nlrf_sph, eqos, pressure_sph, .FALSE. )

      ENDIF

      !PRINT *, pressure_sph(10)
      !PRINT *, pressure_id(10)
      !STOP

      !
      !-- Assign artificial pressure to the real particles
      !
      art_pr_max= zero
      err_N_max=  zero
      err_N_min=  HUGE(one)
      err_N_mean= zero

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( npart_real, nstar_sph, nstar_id, &
      !$OMP                     pressure_sph, pressure_id, &
      !$OMP                     dNstar, art_pr, use_pressure ) &
      !$OMP             PRIVATE( a )
      assign_artificial_pressure_on_real_particles: DO a= 1, npart_real, 1

        IF( nstar_id(a) <= zero )THEN

          dNstar(a)= zero

        ELSEIF(use_pressure)THEN

          dNstar(a)= ( pressure_sph(a) - pressure_id(a) )/pressure_id(a)

        ELSE

          dNstar(a)= ( nstar_sph(a) - nstar_id(a) )/nstar_id(a)

        ENDIF
        art_pr(a) = MAX( one + dNstar(a), one/ten )

      ENDDO assign_artificial_pressure_on_real_particles
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( npart_all, npart_real, art_pr, nstar_sph, &
      !$OMP                     nstar_id, dNstar, all_pos ) &
      !$OMP             PRIVATE( a )
      find_nan_in_art_pr: DO a= 1, npart_real, 1

        IF( .NOT.is_finite_number(art_pr(a)) )THEN
          PRINT *, "** ERROR! art_pr(", a, ")= ", art_pr(a), &
                   " is not a finite number on a real particle!"
          PRINT *, "   nstar_sph(", a, ")=", nstar_sph(a)
          PRINT *, "   nstar_id(", a, ")=", nstar_id(a)
          PRINT *, "   dNstar(", a, ")=", dNstar(a)
          PRINT *, "   rho(", a, ")=", get_density( all_pos(1,a), &
                                                    all_pos(2,a), &
                                                    all_pos(3,a) )
          PRINT *, " * Stopping..."
          PRINT *
          STOP
        ENDIF

      ENDDO find_nan_in_art_pr
      !$OMP END PARALLEL DO

      art_pr_max= - HUGE(one)
      err_N_max = - HUGE(one)
      err_N_min =   HUGE(one)
      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( npart_real, art_pr, dNstar ) &
      !$OMP             PRIVATE( a ) &
      !$OMP             REDUCTION( MAX: art_pr_max, err_N_max ) &
      !$OMP             REDUCTION( MIN: err_N_min )
      DO a= 1, npart_real, 1
        art_pr_max= MAX( art_pr_max, art_pr(a) )
        err_N_max = MAX( err_N_max, ABS(dNstar(a)) )
        err_N_min = MIN( err_N_min, ABS(dNstar(a)) )
      ENDDO
      !$OMP END PARALLEL DO
      IF( .NOT.is_finite_number( art_pr_max ) )THEN
        PRINT *, "** ERROR! art_pr_max is not a finite number!", &
                 " Stopping.."
        PRINT *
        STOP
      ENDIF

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( npart_real, dNstar, nstar_id ) &
      !$OMP             PRIVATE( a ) &
      !$OMP             REDUCTION( +: err_N_mean )
      DO a= 1, npart_real, 1
        err_N_mean= err_N_mean + ABS(dNstar(a))
      ENDDO
      !$OMP END PARALLEL DO
      err_N_mean= err_N_mean/npart_real
      err_N_mean_min= MIN( err_N_mean, err_N_mean_min )

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( all_pos, npart_real, nstar_sph, nstar_id, &
      !$OMP                     dNstar, err_N_max, &
      !$OMP                     pos_maxerr, nstar_sph_err, nstar_id_err ) &
      !$OMP             PRIVATE( a )
      DO a= 1, npart_real, 1

        IF( dNstar(a) == err_N_max )THEN
          pos_maxerr   = all_pos(:,a)
          nstar_sph_err= nstar_sph(a)
          nstar_id_err = nstar_id(a)
        ENDIF

        IF( .NOT.is_finite_number(dNstar(a)) )THEN
          PRINT *, "** ERROR! dNstar(", a, ")= ", dNstar(a), &
                   " is not a finite number on a real particle!"
          PRINT *, "   nstar_sph= ", nstar_sph(a)
          PRINT *, "   nstar_id= ", nstar_id(a)
          STOP
        ENDIF

      ENDDO
      !$OMP END PARALLEL DO

      IF( .NOT.is_finite_number( nu_all ) )THEN
        PRINT *, "** ERROR! nu_all is not a finite number!", &
                 " * Stopping.."
        PRINT *
        STOP
      ENDIF

      !
      !-- Assign artificial pressure to the ghost particles
      !
      IF(adapt_ghosts)THEN

        IF(debug) PRINT *, "adapt_ghosts=", adapt_ghosts
        IF(debug) PRINT *

        nstar_id( npart_real+1:npart_all )= nstar_id_av
        IF(debug) PRINT*, "nstar_id_av=", nstar_id_av

        !$OMP PARALLEL DO DEFAULT( NONE ) &
        !$OMP             SHARED( npart_real, npart_all, nstar_sph, nstar_id, &
        !$OMP                     dNstar, art_pr, nstar_id_av, pressure_id_av, &
        !$OMP                     pressure_sph, use_pressure ) &
        !$OMP             PRIVATE( a )
        assign_artificial_pressure_on_ghost_particles_adapt: &
        DO a= npart_real + 1, npart_all, 1

          IF(use_pressure)THEN

            art_pr(a)= MAX(one + (pressure_sph(a) - pressure_id_av) &
                                  /pressure_id_av, zero)

          ELSE

            art_pr(a)= MAX(one + (nstar_sph(a) - nstar_id_av)/nstar_id_av, zero)

          ENDIF

        ENDDO assign_artificial_pressure_on_ghost_particles_adapt
        !$OMP END PARALLEL DO

      ELSE

        nstar_id( npart_real+1:npart_all )= zero
        !art_pr ( npart_real+1:npart_all )= 6.0D0*art_pr_max

        !$OMP PARALLEL DO DEFAULT( NONE ) &
        !$OMP             SHARED( all_pos, npart_all, npart_real, center, &
        !$OMP                     dNstar, art_pr, rad_x, rad_y, rad_z, &
        !$OMP                     art_pr_max, itr, ellipse_thickness ) &
        !$OMP             PRIVATE( a, r, theta, phi, x_ell, y_ell, z_ell, r_ell, &
        !$OMP                      i_shell )
        assign_artificial_pressure_on_ghost_particles: &
        !
        !-- Assign a pressure to the ghosts, that grows linearly with the
        !-- distance from the center of the matter object
        !
        DO a= npart_real + 1, npart_all, 1

          CALL spherical_from_cartesian( &
            all_pos(1,a), all_pos(2,a), all_pos(3,a), &
            center(1), center(2), center(3), &
            r, theta, phi )

          CALL cartesian_from_spherical( &
            rad_x, theta, phi, &
            center(1), center(2), center(3), &
            x_ell, y_ell, z_ell, rad_y/rad_x, rad_z/rad_x )

          r_ell= SQRT( ( x_ell - center(1) )**two &
                     + ( y_ell - center(2) )**two &
                     + ( z_ell - center(3) )**two )

          shell_loop: DO i_shell= 1, 10, 1

            IF( r <= (one + (ellipse_thickness - one)*DBLE(i_shell)/ten)*r_ell &
                .AND. &
            r >= (one + (ellipse_thickness - one)*DBLE(i_shell - 1)/ten)*r_ell &
            )THEN

              art_pr(a)= DBLE(3*i_shell)*art_pr_max
              IF( .NOT.is_finite_number(art_pr(a)) &
                  .OR. &
                  art_pr(a) > max_art_pr_ghost &
              )THEN
                art_pr(a)= max_art_pr_ghost
              ENDIF
              EXIT

            ENDIF

          ENDDO shell_loop

        ENDDO assign_artificial_pressure_on_ghost_particles
        !$OMP END PARALLEL DO

      ENDIF

      !STOP

      IF( debug ) PRINT *, "Before calling position_correction"

      IF( debug ) PRINT *, npart_all
      IF( debug ) PRINT *, SIZE(all_pos(1,:))
      IF( debug ) PRINT *, SIZE(h)
      IF( debug ) PRINT *, SIZE(art_pr)
      IF( debug ) PRINT *, SIZE(nstar_sph)
      IF( debug ) PRINT *, SIZE(correction_pos(1,:))

      !
      !-- The following loop shouldn't be needed, but apparently
      !-- the test in the previous loop is not enough to remove
      !-- random NaNs from the artificial pressure on the ghost particles.
      !-- Which, btw, shouldn't be there at all since all the quantities are
      !-- tested, and no errors are detected...
      !-- Also, this loop is needed only when compiling with gfortran.
      !
      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( npart_all, npart_real, art_pr, nstar_sph, &
      !$OMP                     nstar_id, all_pos ) &
      !$OMP             PRIVATE( a )
      fix_nan_in_art_pr_ghost: DO a= npart_real + 1, npart_all, 1

        IF( .NOT.is_finite_number(art_pr(a)) )THEN
          art_pr(a)= max_art_pr_ghost
          !PRINT *, "** ERROR! art_pr(", a, ")= ", art_pr(a), &
          !         " is not a finite number on a ghost particle!"
          !PRINT *, "   nstar_sph(", a, ")=", nstar_sph(a)
          !PRINT *, "   nstar_id(", a, ")=", nstar_id(a)
          !PRINT *, "   rho(", a, ")=", get_density( all_pos(1,a), &
          !                                          all_pos(2,a), &
          !                                          all_pos(3,a) )
          !PRINT *, " * Stopping.."
          !PRINT *
          !STOP
        ENDIF

      ENDDO fix_nan_in_art_pr_ghost
      !$OMP END PARALLEL DO

      PRINT *, " * Maximum relative error between the star density profile", &
               " and its SPH estimate: err_N_max= ", err_N_max
      PRINT *, "     at position: x=", pos_maxerr(1), ", y=", pos_maxerr(2), &
               ", z=", pos_maxerr(3)
      PRINT *, "     with r/(system size)= ", SQRT( &
                              ( ABS(pos_maxerr(1)) - ABS(center(1)) )**two &
                            + ( ABS(pos_maxerr(2)) - ABS(center(2)) )**two &
                            + ( ABS(pos_maxerr(3)) - ABS(center(3)) )**two ) &
                            /sizes(1)
      PRINT *, "   The ID density is   = ", nstar_id_err
      PRINT *, "   The SPH estimate is= ", nstar_sph_err
      PRINT *
      PRINT *, " * Minimum relative error between the star density profile", &
               " and its SPH estimate: ", err_N_min
      PRINT *, " * Average relative error between the star density profile", &
               " and its SPH estimate: ", err_N_mean
      PRINT *, " * Minimum of the average relative error between the star", &
               " density profile and its SPH estimate: ", err_N_mean_min
      PRINT *

      !
      !-- Compute what would be the baryon number at this step
      !

      ! Compute particle number density
      nu= one
      CALL density_loop( npart_all, all_pos, &    ! input
                         nu, h, nstar_sph )      ! output
      ! Reset nu to its values immediately
      nu(1:npart_real)          = nu_all
      nu(npart_real+1:npart_all)= nu_ghost

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( nu_tmp, nu, nstar_id, nstar_sph, &
      !$OMP                     nuratio_thres, npart_real ) &
      !$OMP             PRIVATE( nu_tmp2, a )
      cap_nu: DO a= 1, npart_real, 1

        nu_tmp2= nu(a)
        nu_tmp(a)= nstar_id(a)/nstar_sph(a)

          IF( nu_tmp(a) > nu_tmp2*SQRT(nuratio_thres) ) nu_tmp(a)= &
                                            nu_tmp2*SQRT(nuratio_thres)
          IF( nu_tmp(a) < nu_tmp2/SQRT(nuratio_thres) ) nu_tmp(a)= &
                                            nu_tmp2/SQRT(nuratio_thres)

      ENDDO cap_nu
      !$OMP END PARALLEL DO
      nuratio_tmp= MAXVAL( nu_tmp(1:npart_real), DIM= 1 )&
                  /MINVAL( nu_tmp(1:npart_real), DIM= 1 )

      PRINT *, " * Stopping the APM iteration at this step, with a threshold", &
               " for the baryon number ratio equal to "
      PRINT *, "   nu_thres=", nuratio_thres, ","
      PRINT *, "   the baryon number ratio would be equal to the following."
      PRINT *, " * Temporary CORRECTED maximum baryon number at this step=", &
               MAXVAL( nu_tmp(1:npart_real), DIM= 1 )
      PRINT *, " * Temporary CORRECTED minimum baryon number at this step=", &
               MINVAL( nu_tmp(1:npart_real), DIM= 1 )
      PRINT *, " * Temporary CORRECTED baryon number ratio at this step=", &
               nuratio_tmp
      PRINT *

      !
      !-- Setting variables to steer the APM iteration
      !
      push_away_ghosts= .FALSE.
      IF( err_N_mean > err_mean_old )THEN
        n_inc= n_inc + 1
        !IF( n_inc >= max_inc*half )THEN
        !  push_away_ghosts= .TRUE.
        !ENDIF
      ENDIF
      IF( itr > nuratio_min_it .AND. nuratio_tmp /= nuratio_thres .AND. &
          ABS(nuratio_tmp - nuratio_tmp_prev)/nuratio_tmp_prev &
          <= nuratio_tol )THEN
        nuratio_cnt= nuratio_cnt + 1
        !push_away_ghosts= .TRUE.
      ELSE
        nuratio_cnt= 0
      ENDIF

      radius_part= zero
      h_av       = zero
      cnt_array  = zero
      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( all_pos, npart_real, center, max_radius, h, &
      !$OMP                     cnt_array ) &
      !$OMP             PRIVATE( a, r ) &
      !$OMP             REDUCTION( MAX: radius_part ) &
      !$OMP             REDUCTION( +: h_av )
      find_radius_part: DO a= 1, npart_real, 1

        r= SQRT((all_pos(1,a) - center(1))**2 + (all_pos(2,a) - center(2))**2 &
              + (all_pos(3,a) - center(3))**2)

        radius_part= MAX(radius_part, r)

        IF(r > (one - five/(ten*ten))*max_radius)THEN
          h_av= h_av + h(a)
          cnt_array(a)= one
        ENDIF

      ENDDO find_radius_part
      !$OMP END PARALLEL DO
      h_av= h_av/SUM(cnt_array, DIM=1)
      PRINT *, " * Larger radius among the particles=", radius_part
      PRINT *, " * Larger radius of the star=", max_radius
      PRINT *, " * Their difference: radius_part - max_radius=", &
               radius_part - max_radius
      PRINT *, " * Their relative difference: ", &
               "ABS(radius_part - max_radius)/max_radius=", &
               ABS(radius_part - max_radius)/max_radius
      PRINT *, " * Average smoothing length over the outer layers ", &
               "(r>95% of the star radius)= ", h_av, "Msun_geo="
      PRINT *, "   ", h_av*Msun_geo, "km=", &
               h_av/max_radius*ten*ten, "% of the larger radius of the star"
      PRINT *, " * Ghosts are moved if ", &
               "ABS(radius_part - max_radius) > h_av/four=", h_av/four
      PRINT *

      !IF( ABS(radius_part - max_radius)/max_radius > four*ten*tol )THEN
      IF( ABS(radius_part - max_radius) > h_av/four )THEN
        push_away_ghosts= .TRUE.
        IF( radius_part > max_radius )THEN
          ghost_displacement= - ghost_displacement
        ENDIF
      ENDIF

      ! POSSIBLE EXIT CONDITION. DEPRECATED?
      !
      !IF( ABS( err_N_mean - err_mean_old )/ABS( err_mean_old ) < iter_tol &
      !    .AND. &
      !    err_N_max < ten &
      !)THEN
      !  n_inc= n_inc + 1
      !  PRINT *, "n_inc/max_inc= ", n_inc, "/", max_inc
      !  PRINT *, "ABS( err_N_mean - err_mean_old )/ABS(err_mean_old)= ", &
      !           ABS( err_N_mean - err_mean_old )/ABS(err_mean_old)
      !ENDIF
      !IF( ABS(err_N_mean_min - err_N_mean_min_old)/ABS(err_N_mean_min_old) &
      !      < ten*iter_tol &
      !    .AND. &
      !    ABS(err_N_mean - err_N_mean_min)/ABS(err_N_mean_min) < iter_tol &
      !)THEN
      !  n_inc= n_inc + 1
      !  PRINT *, "n_inc/max_inc= ", n_inc, "/", max_inc
      !  PRINT *, err_N_mean, "err_N_mean"
      !  PRINT *, err_N_mean_min, "err_N_mean_min"
      !  PRINT *, "ABS(err_N_mean - err_N_mean_min)/ABS(err_N_mean_min)= ", &
      !           ABS(err_N_mean - err_N_mean_min)/ABS(err_N_mean_min), " < ", &
      !           iter_tol
      !ELSE
      !  n_inc= 0
      !ENDIF

      !
      !-- Estimating the particle contribution to the SPH momentum equation
      !-- Note that it is much easier to compute the momentum in
      !-- SPHINCS_BSSN, so this feature will most likely be deleted
      !

      !IF( MOD( itr, 5 ) == 0 )THEN
      !
      !  !PRINT *, "Before calling exact_nei_tree_update..."
      !  !PRINT *
      !
      !  CALL exact_nei_tree_update( nn_des, &
      !                              npart_real, &
      !                              all_pos(:,1:npart_real), &
      !                              nu_tmp(1:npart_real) )
      !
      !  CALL compute_hydro_momentum()
      !
      !ENDIF

      !STOP

      !
      !-- Exit conditions
      !
      IF(itr == apm_max_it)THEN

        PRINT *, " * Exit condition satisfied: the APM reached the maximum ", &
                 "number of iterations apm_max_it=", apm_max_it
        PRINT *
        EXIT

      ENDIF
      IF(nuratio_des > zero)THEN
      ! If there is a desired baryon number ratio...

        IF( (nuratio_tmp >= nuratio_des*(one - quarter/ten) .AND. &
             nuratio_tmp <= nuratio_des*(one + quarter/ten) .AND. &
             nuratio_tmp /= nuratio_thres) )THEN
        ! ...and the actual baryon number ratio differs from it by 2.5% at the
        ! most, but it's not equal to the baryon number ratio threshold
        ! (the 'cap' on nu), then EXIT the APM iteration. Also, EXIT if
        ! the hard limit on the number of APM steps (equal to apm_max_it)
        ! is reached.

          PRINT *, " * Exit condition satisfied: the baryon number ratio is ", &
                   nuratio_tmp, " <= ", nuratio_des, "*1.025= ", &
                   nuratio_des*(one + quarter/ten)
          PRINT *
          EXIT

        ENDIF

        IF( (nuratio_cnt >= nuratio_max_steps .AND. &
             nuratio_tmp /= nuratio_thres) &
        )THEN
        ! ...and the actual baryon number ratio did not change significantly
        ! for more than nuratio_max_steps steps, but it's not equal to the
        ! baryon number ratio threshold (the 'cap' on nu), then EXIT the APM
        ! iteration. Also, EXIT if the hard limit on the number of APM steps
        ! (equal to apm_max_it) is reached.

          PRINT *, " * Exit condition satisfied: the baryon number ratio ", &
                   "did not change by more than ", nuratio_tol*ten*ten, &
                   "% for ", nuratio_max_steps, "steps."
          PRINT *
          EXIT

        ELSE
        ! Print the counter

          PRINT *, " * nuratio_cnt/nuratio_max_steps= ", &
                   nuratio_cnt, "/", nuratio_max_steps
          PRINT *

        ENDIF

      ELSE
      ! else if there is NOT a desired baryon number ratio, EXIT the APM
      ! iteration if the relative error on the density increases max_inc times.
      ! Also, EXIT if the hard limit on the number of APM steps
      ! (equal to apm_max_it) is reached.

        PRINT *, " * n_inc= ", n_inc
        PRINT *
        IF(n_inc == max_inc)THEN

          PRINT *, " * Exit condition satisfied: the average over the ", &
                   "particles of the relative difference between the ID  ", &
                   "baryon mass density and its SPH estimate, grew ", &
                   "for ", max_inc, "steps."
          PRINT *
          EXIT

        ENDIF

      ENDIF
      err_mean_old      = err_N_mean
      err_N_mean_min_old= err_N_mean_min
      nuratio_tmp_prev  = nuratio_tmp

      !
      !-- If the EXIT conditions are not satisfied, update the particle
      !-- distribution
      !
      PRINT *, " * Updating positions..."

      all_pos_prev= all_pos

      CALL density_loop( npart_all, all_pos, &    ! input
                         nu, h, nstar_sph )      ! output

      CALL position_correction( npart_all, &
                                all_pos, h, nu_all, art_pr, nstar_sph, &
                                correction_pos )

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( npart_all, correction_pos, all_pos, center, &
      !$OMP                     radius_x ) &
      !$OMP             PRIVATE( a, itr2, r, theta, phi )
      find_nan_in_correction_pos: DO a= 1, npart_all, 1

        loop_over_spatial_components: DO itr2= 1, 3, 1

          IF( .NOT.is_finite_number( correction_pos( itr2, a ) ) )THEN

            CALL spherical_from_cartesian( all_pos(1,a), all_pos(2,a), &
                                           all_pos(3,a), &
                                           center(1), center(2), center(3), &
                                           r, theta, phi )

          !  correction_pos(1,a)= - one/(two*ten)*SIN(theta)*COS(phi)
          !  correction_pos(2,a)= - one/(two*ten)*SIN(theta)*SIN(phi)
          !  correction_pos(3,a)= - one/(two*ten)*COS(theta)
            !correction_pos( itr2, a )= zero

            PRINT *, "** ERROR! correction_pos(", itr2, ",", a, ") is a NaN!"
            PRINT *, " *        correction_pos: x=", correction_pos(1,a), &
                     ", y=", correction_pos(2,a), ", z=", correction_pos(3,a)
            PRINT *, " *        Particle position: x=", all_pos(1,a), &
                     ", y=", all_pos(2,a), ", z=", all_pos(3,a)

            CALL spherical_from_cartesian( &
                              all_pos(1,a), all_pos(2,a), all_pos(3,a), &
                              center(1), center(2), center(3), &
                              r, theta, phi )

            !r_tmp= SQRT( ( all_pos(1,a) - center(1) )**two + &
            !             ( all_pos(2,a) - center(2) )**two + &
            !             ( all_pos(3,a) - center(3) )**two )

            PRINT *, " *        Particle radius: r=", r, &
                     "=", r/radius_x*ten*ten, &
                     "% of the larger radius of the star."
            PRINT *, " *        Particle colatitude: theta=", theta/pi," pi"
            PRINT *, " *        Particle longitude: phi=", phi/pi, " pi"
            PRINT *, " * Stopping.."
            PRINT *
            STOP

          ENDIF

        ENDDO loop_over_spatial_components

      ENDDO find_nan_in_correction_pos
      !$OMP END PARALLEL DO

      IF( debug ) PRINT *, "After calling position_correction"

      cnt_move= 0
      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( use_atmosphere, all_pos, correction_pos, &
      !$OMP                     dNstar, npart_real, nstar_id, cnt_move ) &
      !$OMP             PRIVATE( pos_corr_tmp, a, cnt, rand_num, rand_num2, &
      !$OMP                      rel_sign )
      real_particle_loop: DO a= 1, npart_real, 1

        adapt_displacement_to_error: &
        IF( dNstar(a) >= ten*ten &
            .AND. &
            validate_position_final( &
              all_pos(1,a) + ten*correction_pos(1,a), &
              all_pos(2,a) + ten*correction_pos(2,a), &
              all_pos(3,a) + ten*correction_pos(3,a) ) )THEN

          pos_corr_tmp= all_pos(:,a) + ten*correction_pos(:,a) ! 10


        ELSEIF( dNstar(a) >= ten &
                .AND. &
                validate_position_final( &
                  all_pos(1,a) + three*correction_pos(1,a), &
                  all_pos(2,a) + three*correction_pos(2,a), &
                  all_pos(3,a) + three*correction_pos(3,a) ) )THEN

          pos_corr_tmp= all_pos(:,a) + three*correction_pos(:,a) ! 3


        ELSE

          pos_corr_tmp= all_pos(:,a) + correction_pos(:,a) ! 1

        ENDIF adapt_displacement_to_error

        if_atmosphere: IF( use_atmosphere )THEN
        ! If the atmosphere is used...

          all_pos(:,a)= pos_corr_tmp
          cnt_move(a)= 1
          !...move the particle without any validation

        ELSE

          cnt= 0
          determine_new_position: DO

            test_position: IF( get_density( &
                pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) > zero &
                .AND. &
                validate_position_final( &
                    pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) &
                !.AND. &
                !check_particle_position( a - 1, &
                !                         all_pos(:,1:a-1), &
                !                         pos_corr_tmp ) == 0 &
                !.AND. &
                !check_particle_position( npart_real - a, &
                !                         all_pos(:,a+1:npart_real), &
                !                         pos_corr_tmp ) == 0 &
            )THEN
            ! If the new position is valid...

              all_pos(:,a)= pos_corr_tmp
              cnt_move(a)= 1
              EXIT
              !...move the particle, and exit the 'determine_new_position' loop

            ELSEIF( cnt <= search_pos )THEN
            ! ...else if the new position is invalid,
            !    and the current step is lower than search_pos, change the new
            !    position randomly, independently in x, y, z

              cnt= cnt + 1

              !
              !-- Add random noise to x coordinate
              !
              CALL RANDOM_NUMBER( rand_num )
              CALL RANDOM_NUMBER( rand_num2 )

              IF( rand_num2 < half )  rel_sign= - 1
              IF( rand_num2 >= half ) rel_sign=   1

              pos_corr_tmp(1)= all_pos(1,a) + &
                correction_pos(1,a)*( one + DBLE(rel_sign)*rand_num*half )

              !
              !-- Add random noise to y coordinate
              !
              CALL RANDOM_NUMBER( rand_num )
              CALL RANDOM_NUMBER( rand_num2 )

              IF( rand_num2 < half )  rel_sign= - 1
              IF( rand_num2 >= half ) rel_sign=   1

              pos_corr_tmp(2)= all_pos(2,a) + &
                correction_pos(2,a)*( one + DBLE(rel_sign)*rand_num*half )

              !
              !-- Add random noise to z coordinate
              !
              CALL RANDOM_NUMBER( rand_num )
              CALL RANDOM_NUMBER( rand_num2 )

              IF( rand_num2 < half )  rel_sign= - 1
              IF( rand_num2 >= half ) rel_sign=   1

              pos_corr_tmp(3)= all_pos(3,a) + &
                correction_pos(3,a)*( one + DBLE(rel_sign)*rand_num*half )

              !pos_corr_tmp*( one + DBLE(rel_sign)*rand_num*half*third )

            ELSEIF( cnt == search_pos + 1 )THEN
            ! ...else if the new position was changed randomly search_pos
            !    times, do not move the particle at this step,
            !    and exit the 'determine_new_position' loop

              cnt_move(a)= 0

              ! cnt= cnt + 1
              ! CALL RANDOM_NUMBER( rand_num )
              ! CALL RANDOM_NUMBER( rand_num2 )
              !
              ! IF( rand_num2 < half )  rel_sign= - 1
              ! IF( rand_num2 >= half ) rel_sign=   1
              ! all_pos(:,a)= all_pos(:,a)*( one -rand_num*half*third )

              EXIT

            ENDIF test_position

          ENDDO determine_new_position

        ENDIF if_atmosphere

      ENDDO real_particle_loop
      !$OMP END PARALLEL DO
      PRINT *, " * The fraction of particles that moved at this step is", &
               DBLE(SUM(cnt_move))/DBLE(npart_real)
      PRINT *

      !
      !-- Compute some measures of the displacement vector field over the
      !-- particles
      !
      l2norm_displacement = zero
      average_displacement= zero
      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( npart_real, all_pos, all_pos_prev ) &
      !$OMP             PRIVATE( a ) &
      !$OMP             REDUCTION( +: l2norm_displacement, average_displacement)
      l2norm_disp: DO a= 1, npart_real, 1

        ! l2 norm
        l2norm_displacement= l2norm_displacement &
                           + NORM2(all_pos(:,a) - all_pos_prev(:,a))**2

        ! Arithmetic average of the Euclidean norm of the displacement
        average_displacement= average_displacement &
                            + NORM2(all_pos(:,a) - all_pos_prev(:,a))

      ENDDO l2norm_disp
      !$OMP END PARALLEL DO
      l2norm_displacement= SQRT(l2norm_displacement)
      average_displacement= l2norm_displacement/npart_real
      PRINT *, " * l_2 norm of the displacement of the particles= ", &
               l2norm_displacement
      PRINT *, "   (note that the l2 norm grows with the number of particles)"
      PRINT *, " * Arithmetic average of the Euclidean norm of the ", &
               "displacement of the particles= ", average_displacement
      PRINT *

      !
      !-- Displace ghosts, if needed
      !
      IF(debug) PRINT *, "push_away_ghosts:", push_away_ghosts
      IF(debug) PRINT *, "move_away_ghosts:", move_away_ghosts
      IF(debug) PRINT *, "push_away_ghosts .AND. move_away_ghosts:", &
                         push_away_ghosts .AND. move_away_ghosts
      IF(debug) PRINT *, "cnt_push_ghost=", cnt_push_ghost
      IF(debug) PRINT *, "max_push_ghost=", max_push_ghost

      IF( push_away_ghosts .AND. move_away_ghosts &
      !    .AND. cnt_push_ghost <= max_push_ghost &
      )THEN

        PRINT *, " * Displacing ghosts..."

        !max_r_ghost= (one + third)*MAXVAL([radius_x, radius_y, radius_z])

        !$OMP PARALLEL DO DEFAULT( NONE ) &
        !$OMP             SHARED( all_pos, npart_real, npart_all, center, &
        !$OMP                     ghost_displacement, min_radius ) &
        !$OMP             PRIVATE( a, r, theta, phi )
        ghost_particle_loop: DO a= npart_real + 1, npart_all, 1

          CALL spherical_from_cartesian( &
            all_pos(1,a), all_pos(2,a), all_pos(3,a), &
            center(1), center(2), center(3), &
            r, theta, phi )

          r= r + min_radius*ghost_displacement

          CALL cartesian_from_spherical( &
            r, theta, phi, &
            center(1), center(2), center(3), &
            all_pos(1,a), all_pos(2,a), all_pos(3,a) )

          !all_pos(1,a)= center(1) + r*SIN(theta)*COS(phi)
          !all_pos(2,a)= center(2) + r*SIN(theta)*SIN(phi)
          !all_pos(3,a)= center(3) + r*COS(theta)

         ! all_pos(:,a)= center + ghost_displacement*(all_pos(:,a) - center)

         ! This will break the ghost ellipsoid into eight disconnected pieces,
         ! one per quadrant in the 3D space
         ! all_pos(1,a)= all_pos(1,a) &
         !    + SIGN(min_radius*ghost_displacement, all_pos(1,a) - center(1))
         ! all_pos(2,a)= all_pos(2,a) &
         !    + SIGN(min_radius*ghost_displacement, all_pos(2,a) - center(2))
         ! all_pos(3,a)= all_pos(3,a) &
         !    + SIGN(min_radius*ghost_displacement, all_pos(3,a) - center(3))

        ENDDO ghost_particle_loop
        !$OMP END PARALLEL DO

        cnt_push_ghost= cnt_push_ghost + 1

        PRINT *, "   ...done."
        PRINT *

      ENDIF

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( npart_all, all_pos ) &
      !$OMP             PRIVATE( a, itr2 )
      find_nan_in_all_pos: DO a= 1, npart_all, 1

        DO itr2= 1, 3, 1
          IF( .NOT.is_finite_number( all_pos( itr2, a ) ) )THEN
            PRINT *, "** ERROR! all_pos(", itr2, ",", a, ") is a NaN!", &
                     " Stopping.."
            PRINT *
            STOP
          ENDIF
        ENDDO

      ENDDO find_nan_in_all_pos
      !$OMP END PARALLEL DO

      ! If some of the particles crossed the xy plane in the
      ! last step, reflect them back above the xy plane

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( all_pos, all_pos_prev, npart_real ) &
      !$OMP             PRIVATE( a )
      DO a= 1, npart_real, 1

        IF( (all_pos_prev(3,a) > zero .AND. all_pos(3,a) <= zero) &
            .OR. &
            (all_pos_prev(3,a) < zero .AND. all_pos(3,a) >= zero) &
        )THEN

          all_pos(3,a)= all_pos_prev(3,a)

        ENDIF

      ENDDO
      !$OMP END PARALLEL DO

      IF( debug ) PRINT *, "After correcting positions"

    ENDDO apm_iteration

    PRINT *, "** APM iteration completed."
    PRINT *



    !--------------------------!
    !--------------------------!
    !-- END OF APM ITERATION --!
    !--------------------------!
    !--------------------------!



    !-----------------------------!
    !-- Discard ghost particles --!
    !-----------------------------!

    pos= all_pos(:, 1:npart_real)
    IF( debug ) PRINT *, npart

    h      = h(1:npart_real)
    h_guess= h_guess(1:npart_real)
    nu     = nu(1:npart_real)

    !-----------------------------------------------!
    !-- Remove atmosphere, if present and desired --!
    !-----------------------------------------------!

    IF( use_atmosphere .AND. remove_atmosphere )THEN

      CALL discard_atmosphere( npart_real )

    ENDIF

    !---------------!
    !-- Set npart --!
    !---------------!

    npart= npart_real

    !----------------------------!
    !-- enforce centre of mass --!
    !----------------------------!

    CALL correct_center_of_mass( npart_real, pos, nu, get_density, &
                                 validate_position_final, com_star, &
                                 verbose= .TRUE. )

    !-----------------------------------------------------------------------!
    !-- Mirror the positions after having repositioned the center of mass --!
    !-----------------------------------------------------------------------!

    CALL impose_equatorial_plane_symmetry( npart_real, pos, nu )

    !-----------------------------!
    !-- Print positions to file --!
    !-----------------------------!

    IF( PRESENT(namefile_pos) )THEN
      finalnamefile= namefile_pos
    ELSE
      finalnamefile= TRIM(sph_path)//"apm_pos.dat"
    ENDIF

    INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

    IF( exist )THEN
        OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
              FORM= "FORMATTED", &
              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
              IOMSG= err_msg )
    ELSE
        OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "NEW", &
              FORM= "FORMATTED", &
              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening " // TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF

    DO a= 1, npart_real, 1
      tmp= get_density( pos(1,a), pos(2,a), pos(3,a) )
      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        1, a, &
        pos(1,a), &
        pos(2,a), &
        pos(3,a), &
        tmp, cnt_move(a)
    ENDDO

    DO a= npart_real + 1, npart_all, 1
      tmp= get_density( all_pos(1,a), all_pos(2,a), all_pos(3,a) )
      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        2, a, &
        all_pos(1,a), &
        all_pos(2,a), &
        all_pos(3,a), &
        tmp
    ENDDO

    CLOSE( UNIT= 2 )

    !---------------------------------------------------------------!
    !-- Assign baryon number to match profile as good as possible --!
    !---------------------------------------------------------------!

    PRINT *, " * Assign baryon number..."
    PRINT *

    IF( debug ) PRINT *, "1"

    CALL assign_h( nn_des, &           !
                   npart_real, &        !
                   pos, h_guess, & ! Input
                   h )                 ! Output

    find_h_bruteforce_timer= timer( "find_h_bruteforce_timer" )
    CALL find_h_bruteforce_timer% start_timer()
    n_problematic_h= 0
    check_h3: DO a= 1, npart_real, 1
    ! find_h_backup, called below, is OMP parallelized, so this loop
    ! should not be parallelized as well

      IF( .NOT.is_finite_number( h(a) ) .OR. h(a) <= zero )THEN

        n_problematic_h= n_problematic_h + 1
        h(a)= find_h_backup( a, npart_real, pos, nn_des )
        IF( .NOT.is_finite_number( h(a) ) .OR. h(a) <= zero )THEN
          PRINT *, "** ERROR! h=0 on particle ", a, "even with the brute", &
                   " force method."
          PRINT *, "   Particle position: ", pos(:,a)
          STOP
        ENDIF

      ENDIF

    ENDDO check_h3
    CALL find_h_bruteforce_timer% stop_timer()
    CALL find_h_bruteforce_timer% print_timer( 2 )

    IF( debug ) PRINT *, "2"

    ! Measure SPH particle number density
    nu= one
    CALL density_loop( npart_real, pos, &    ! input
                       nu, h, nstar_sph )      ! output

    IF( debug ) PRINT *, "3"

    CALL get_nstar_id_atm( npart_real, pos(1,:), pos(2,:), pos(3,:), &
                           nstar_sph, nstar_id, nlrf_sph, sqg, &
                           use_atmosphere )

    nu(1:npart_real)= nu_all
    PRINT *, " * Baryon number on all particles before correction nu_all= ", &
             nu_all

    nu_tot= zero
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( nu, npart_real ) &
    !$OMP             PRIVATE( a ) &
    !$OMP             REDUCTION(+: nu_tot)
    DO a= 1, npart_real, 1
      nu_tot= nu_tot + nu(a)
    ENDDO
    !$OMP END PARALLEL DO

    PRINT *, " * Total baryon number nu_tot=", nu_tot
    PRINT *, " * Total baryon mass= ", nu_tot*amu/MSun, "=", &
             ten*ten*nu_tot*amu/MSun/mass, "% of the ID baryon mass"
    PRINT *

    IF( debug ) PRINT *, "4"

    IF( debug ) PRINT *, "npart_real= ", npart_real
    IF( debug ) PRINT *, "SIZE(nu)= ", SIZE(nu)
    IF( debug ) PRINT *

    IF( debug ) nu_ratio= MAXVAL( nu, DIM= 1 )/MINVAL( nu, DIM= 1 )
    IF( debug ) PRINT *, " * nu_ratio before correction = ", nu_ratio
    IF( debug ) PRINT *

    !-----------------------------------------!
    !-- Cap nu by the desired nuratio_thres --!
    !-----------------------------------------!

    PRINT *, " * Correcting nu..."
    PRINT *

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( nu_tmp, nu, nstar_id, nstar_sph, &
    !$OMP                     nuratio_thres, npart_real ) &
    !$OMP             PRIVATE( nu_tmp2, a )
    DO a= 1, npart_real, 1

      nu_tmp2= nu(a)
      nu(a)= nstar_id(a)/nstar_sph(a)

        IF( nu(a) > nu_tmp2*SQRT(nuratio_thres) ) nu(a)= &
                                          nu_tmp2*SQRT(nuratio_thres)
        IF( nu(a) < nu_tmp2/SQRT(nuratio_thres) ) nu(a)= &
                                          nu_tmp2/SQRT(nuratio_thres)

    ENDDO
    !$OMP END PARALLEL DO

    !
    !-- Check that nu is acceptable
    !
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( nu, nstar_id, nstar_sph, npart_real ) &
    !$OMP             PRIVATE( a )
    DO a= 1, npart_real, 1

      IF( .NOT.is_finite_number( nu(a) ) )THEN
        PRINT *, " * ERROR! nu(", a, ") is a NaN."
        PRINT *, " nstar_sph(a)= ", nstar_sph(a)
        PRINT *, " nstar_id(a)= ", nstar_id(a)
        PRINT *, " Stopping..."
        PRINT *
        STOP
      ENDIF
      IF( nu(a) < zero )THEN
        PRINT *, " * ERROR! nu(", a, ") is negative."
        PRINT *, " nu(a)= ", nu(a)
        PRINT *, " nstar_sph(a)= ", nstar_sph(a)
        PRINT *, " nstar_id(a)= ", nstar_id(a)
        PRINT *, " Stopping..."
        PRINT *
        STOP
      ENDIF

    ENDDO
    !$OMP END PARALLEL DO

    nu_ratio= MAXVAL( nu, DIM= 1 )/MINVAL( nu, DIM= 1 )
    PRINT *, " * nu_ratio after correction = ", nu_ratio
    PRINT *

    max_nu= zero
    min_nu= HUGE(one)
    DO a= 1, npart_real, 1
      IF( nu(a) > max_nu )THEN
        max_nu= nu(a)
        a_numax= a
      ENDIF
      IF( nu(a) < min_nu )THEN
        min_nu= nu(a)
        a_numin= a
      ENDIF
    ENDDO

    PRINT *, " * Baryon number assigned."
    PRINT *

    !
    !-- Optionally change baryon number without moving the particles
    !
    IF( mass_it )THEN

      ! just a few iterations to NOT get the nu-ratio too large
      mass_iteration: DO itr= 1, m_max_it, 1

         ! measure density
         CALL density_loop( npart_real, pos, &    ! input
                            nu, h, nstar_sph )      ! output


         CALL get_nstar_id_atm( npart_real, pos(1,:), &
                                pos(2,:), &
                                pos(3,:), nstar_sph, nstar_id, nlrf_sph, sqg, &
                                use_atmosphere )

         !nstar_id( npart_real+1:npart_all )= zero

         ! get RELATIVE nu's right
         dN_av= zero
         max_nu= zero
         min_nu= HUGE(one)
         DO a= 1, npart_real, 1
            dN=    (nstar_sph(a)-nstar_id(a))/nstar_id(a)
            nu(a)= nu(a)*(one - dN)
            dN_av= dN_av + dN
            IF( nu(a) > max_nu )THEN
              max_nu= nu(a)
              a_numax= a
            ENDIF
            IF( nu(a) < min_nu )THEN
              min_nu= nu(a)
              a_numin= a
            ENDIF
         ENDDO
         dN_av= dN_av/DBLE(npart_real)

         ! exit condition
         IF( dN_av < tol ) EXIT

      ENDDO mass_iteration

    ENDIF

    PRINT *, " * CORRECTED maximum baryon number at this step=", &
             max_nu
    PRINT *, " * CORRECTED minimum baryon number at this step=", &
             min_nu
    PRINT *, " * CORRECTED baryon number ratio at this step=", &
             max_nu/min_nu
    PRINT *

    max_nu2= zero
    min_nu2= HUGE(one)
    DO a= 1, npart_real, 1
       IF( nu(a) > max_nu2 .AND. a /= a_numax )THEN
         max_nu2= nu(a)
         a_numax2= a
       ENDIF
       IF( nu(a) < min_nu2 .AND. a /= a_numin )THEN
         min_nu2= nu(a)
         a_numin2= a
       ENDIF
    ENDDO

    PRINT *, " * Excluding the absolute max and min of nu:"
    PRINT *
    PRINT *, "   max_nu=", max_nu2
    PRINT *, "        at ", pos(:, a_numax2), " r= ", &
             NORM2( pos(:, a_numax2) )/radius_x
    PRINT *, "   min_nu=", min_nu2
    PRINT *, "        at ", pos(:, a_numin2), " r= ", &
             NORM2( pos(:, a_numin2) )/radius_x
    PRINT *, "   max_nu/min_nu=", max_nu2/min_nu2
    PRINT *

    nu_tot= zero
    DO a= 1, npart_real, 1
      nu_tot= nu_tot + nu(a)
    ENDDO
    mean_nu= nu_tot/npart_real

    variance_nu = zero                       ! compute variance
    DO a = 1, npart_real, 1
      variance_nu = variance_nu + (nu(a) - mean_nu)**two
    END DO
    variance_nu = variance_nu / DBLE(npart_real - 1)
    stddev_nu   = SQRT(variance_nu)            ! compute standard deviation

    PRINT *, "mean_nu=", mean_nu
    PRINT *, "variance_nu=", variance_nu
    PRINT *, "stddev_nu=", stddev_nu
    PRINT *, "stddev_nu/mean_nu=", stddev_nu/mean_nu
    PRINT *

    PRINT *, "Before correcting nu to match the mass of the star..."
    PRINT *
    PRINT *, "nu_tot=", nu_tot
    PRINT *, "mass estimate= ", nu_tot*amu/MSun, "=", &
             ten*ten*nu_tot*amu/MSun/mass, "% of the ID baryon mass"
    PRINT *

    IF( correct_nu )THEN

      nu= nu/(nu_tot*amu/Msun/mass)
      nu_tot= zero
      DO a= 1, npart_real, 1
        nu_tot= nu_tot + nu(a)
      ENDDO

      PRINT *, "After correcting nu to match the mass of the star..."
      PRINT *
      PRINT *, "nu_tot=", nu_tot
      PRINT *, "mass estimate= ", nu_tot*amu/MSun, "=", &
               ten*ten*nu_tot*amu/MSun/mass, "% of the ID baryon mass"
      PRINT *

    ENDIF


  !  finalnamefile= TRIM(sph_path)//"dbg_apm_pos1.dat"
  !
  !  INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )
  !
  !  IF( exist )THEN
  !    OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
  !          FORM= "FORMATTED", &
  !          POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
  !          IOMSG= err_msg )
  !  ELSE
  !    OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "NEW", &
  !          FORM= "FORMATTED", &
  !          ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
  !  ENDIF
  !  IF( ios > 0 )THEN
  !  PRINT *, "...error when opening " // TRIM(finalnamefile), &
  !           ". The error message is", err_msg
  !  STOP
  !  ENDIF
  !
  !  DO a= 1, npart_real, 1
  !  tmp= get_density( pos(1,a), pos(2,a), pos(3,a) )
  !  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
  !    a, &
  !    pos(1,a), &
  !    pos(2,a), &
  !    pos(3,a), &
  !    tmp
  !  ENDDO
  !
  !  CLOSE( UNIT= 2 )

    !----------------------------!
    !-- enforce centre of mass --!
    !----------------------------!

    CALL correct_center_of_mass( npart_real, pos, nu, get_density, &
                                 validate_position_final, com_star, &
                                 verbose= .TRUE. )

    !-----------------------------------------------------------------------!
    !-- Mirror the positions after having repositioned the center of mass --!
    !-----------------------------------------------------------------------!

    CALL impose_equatorial_plane_symmetry( npart_real, pos, nu )


  !  finalnamefile= TRIM(sph_path)//"dbg_apm_pos2.dat"
  !
  !  INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )
  !
  !  IF( exist )THEN
  !    OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
  !          FORM= "FORMATTED", &
  !          POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
  !          IOMSG= err_msg )
  !  ELSE
  !    OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "NEW", &
  !          FORM= "FORMATTED", &
  !          ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
  !  ENDIF
  !  IF( ios > 0 )THEN
  !  PRINT *, "...error when opening " // TRIM(finalnamefile), &
  !           ". The error message is", err_msg
  !  STOP
  !  ENDIF
  !
  !  DO a= 1, npart_real, 1
  !    tmp= get_density( pos(1,a), pos(2,a), pos(3,a) )
  !    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
  !      a, &
  !      pos(1,a), &
  !      pos(2,a), &
  !      pos(3,a), &
  !      tmp
  !  ENDDO
  !
  !  CLOSE( UNIT= 2 )

    !-------------------!
    !-- monitoring... --!
    !-------------------!

    ! measure density
    CALL density_loop( npart_real, pos, &    ! input
                       nu, h, nstar_sph )      ! output

    CALL get_nstar_id_atm( npart_real, pos(1,:), &
                           pos(2,:), &
                           pos(3,:), nstar_sph, nstar_id, nlrf_sph, sqg, &
                           use_atmosphere )

    dN_av = zero
    dN_max= zero
    cnt1= 0
    DO a= 1, npart_real, 1
      IF( get_density( pos(1,a), pos(2,a), pos(3,a) ) > zero )THEN
        dN= ABS(nstar_sph(a)-nstar_id(a))/nstar_id(a)
        dN_av=  dN_av + dN
        dN_max= MAX(dN_max,dN)
        cnt1= cnt1 + 1
      ENDIF
    ENDDO
    dN_av= dN_av/DBLE(cnt1)

    variance_dN = zero                       ! compute variance
    DO a = 1, npart_real, 1
      IF( get_density( pos(1,a), pos(2,a), pos(3,a) ) > zero )THEN
        dN= ABS(nstar_sph(a)-nstar_id(a))/nstar_id(a)
        variance_dN = variance_dN + (dN - dN_av)**two
        cnt1= cnt1 + 1
      ENDIF
    END DO
    variance_dN = variance_dN / DBLE(cnt1)
    stddev_dN   = SQRT(variance_dN)            ! compute standard deviation

    PRINT *, " * Final maximum relative error between the density from the ", &
             "ID and the SPH density estimate: dN_max=", dN_max
    PRINT *, " * Final average relative error between the density from the ", &
             "ID and the SPH density estimate: dN_max=", dN_av
    PRINT *, " * Final variance of the relative error between the density ", &
             "from the ID and the SPH density estimate: variance_dN=", &
             variance_dN
    PRINT *, " * Final standard deviation of the relative error between the ", &
             "density from the ID and the SPH density estimate: stddev_dN=", &
             stddev_dN
    PRINT *, " * Final normalized standard deviation of the relative error ", &
             "between the density from the ID and the SPH density ", &
             "estimate: stddev_dN/dN_a=", stddev_dN/dN_av
    PRINT *

    IF( debug ) PRINT *, "100"

    IF( .NOT.ALLOCATED( nstar_int ) ) ALLOCATE( nstar_int( npart_real ) )

    IF( debug ) PRINT *, "101"

    PRINT *, "** Assigning final smoothing length..."
    PRINT *

    ! Determine smoothing length so that each particle has exactly
    ! 300 neighbours inside 2h
    CALL assign_h( nn_des, &
                   npart_real, &
                   pos, h_guess, & ! Input
                   h )             ! Output

    IF( debug ) PRINT *, "101.5"

    find_h_bruteforce_timer= timer( "find_h_bruteforce_timer" )
    CALL find_h_bruteforce_timer% start_timer()
    n_problematic_h= 0
    check_h4: DO a= 1, npart_real, 1
    ! find_h_backup, called below, is OMP parallelized, so this loop
    ! should not be parallelized as well

      IF( .NOT.is_finite_number( h(a) ) .OR. h(a) <= zero )THEN

        n_problematic_h= n_problematic_h + 1
        h(a)= find_h_backup( a, npart_real, pos, nn_des )
        IF( .NOT.is_finite_number( h(a) ) .OR. h(a) <= zero )THEN
          PRINT *, "** ERROR! h=0 on particle ", a, "even with the brute", &
                   " force method."
          PRINT *, "   Particle position: ", pos(:,a)
          STOP
        ENDIF

      ENDIF

    ENDDO check_h4
    CALL find_h_bruteforce_timer% stop_timer()
    CALL find_h_bruteforce_timer% print_timer( 2 )

    PRINT *, " * The smoothing length was found brute-force for ", &
             n_problematic_h, " particles."
    PRINT *

    IF( SUM(nu, DIM= 1)/SIZE(nu) == zero )THEN
      PRINT *, "** ERROR! Average nu= 0. Are you assigning values to the ", &
               "TYPE member array?"
      PRINT *, "Stopping..."
      STOP
    ENDIF

    PRINT *, "** Building neighbors tree..."
    PRINT *
    cnt1= 0
 !   DO
 !
 !     few_ncand= .FALSE.
 !
 !     ! Redo the previous step slightly different (it's built-in;
 !     ! exact_nei_tree_update does not work if I don't call assign_h first),
 !     ! then update the neighbour-tree and fill the neighbour-data
      CALL exact_nei_tree_update( nn_des, &
                                  npart_real, &
                                  pos, nu )

  !    ll_cell_loop: DO ill= 1, nfinal
  !
  !      itot= nprev + ill
  !      IF( nic(itot) == 0 ) CYCLE
  !
  !      IF( ncand(ill) < nn_des - 1 )THEN
  !
  !        ! Increase the smoothing lengths of the paricles inside the cell,
  !        ! and rebuild the tree
  !        few_ncand= .TRUE.
  !
  !        particle_in_cell_loop: DO l= lpart(itot), rpart(itot)
  !
  !          h(l)= three*h(l)
  !
  !        ENDDO particle_in_cell_loop
  !
  !      ELSE
  !
  !        few_ncand= .FALSE.
  !
  !      ENDIF
  !
  !    ENDDO ll_cell_loop
  !
  !    cnt1= cnt1 + 1
  !
  !    IF( debug ) PRINT *, cnt1
  !
  !    IF( .NOT.few_ncand .OR. cnt1 >= max_it_tree )THEN
  !      EXIT
  !    ENDIF
  !
  !  ENDDO

    find_h_bruteforce_timer= timer( "find_h_bruteforce_timer" )
    CALL find_h_bruteforce_timer% start_timer()
    n_problematic_h= 0
    check_h5: DO a= 1, npart_real, 1
    ! find_h_backup, called below, is OMP parallelized, so this loop
    ! should not be parallelized as well

      IF( .NOT.is_finite_number( h(a) ) .OR. h(a) <= zero )THEN

        n_problematic_h= n_problematic_h + 1
        h(a)= find_h_backup( a, npart_real, pos, nn_des )
        IF( .NOT.is_finite_number( h(a) ) .OR. h(a) <= zero )THEN
          PRINT *, "** ERROR! h=0 on particle ", a, "even with the brute", &
                   " force method."
          PRINT *, "   Particle position: ", pos(:,a)
          STOP
        ENDIF

      ENDIF

    ENDDO check_h5
    CALL find_h_bruteforce_timer% stop_timer()
    CALL find_h_bruteforce_timer% print_timer( 2 )

    PRINT *, " * The smoothing length was found brute-force for ", &
             n_problematic_h, " particles."
    PRINT *

    PRINT *, " * Smoothing lengths assigned and tree is built."
    PRINT *

    IF( debug ) PRINT *, "102"

    CALL density( npart_real, pos, nstar_int )
    !nstar_int= zero

    IF( debug ) PRINT *, "103"

  !  PRINT *, "** Finding nearest neighbors..."

    ALLOCATE( neighbors_lists( npart_real ) )
    ALLOCATE( n_neighbors( npart_real ) )
    ALLOCATE( nearest_neighbors( 2, npart_real ) )

    neighbors_lists= 0
    n_neighbors= 0
    nearest_neighbors(1,:)= 0
    nearest_neighbors(2,:)= HUGE(one)

 !   find_neighbors: DO a= 1, npart_real, 1
 !
 !     CALL get_neighbours_bf( a, npart_real, pos, h, 3, &           ! Input
 !                             n_neighbors(a), neighbors_lists(:) )  ! Output
 !
 !     DO itr= 1, n_neighbors(a), 1
 !
 !       dist= NORM2( pos(:,a) - pos(:,neighbors_lists(itr)) )
 !
 !       !PRINT *, "dist= ", dist
 !       !PRINT *, "nearest_neighbors(2,a)= ", nearest_neighbors(2,a)
 !       !PRINT *, "dist < nearest_neighbors(2,a)= ", dist < nearest_neighbors(2,a)
 !
 !       IF( dist < nearest_neighbors(2,a) )THEN
 !
 !         nearest_neighbors(1,a)= neighbors_lists(itr)
 !         nearest_neighbors(2,a)= dist
 !
 !       ENDIF
 !
 !     ENDDO
 !
 !     !PRINT *, "dist= ", dist
 !     !PRINT *, "nearest_neighbors(2,a)= ", nearest_neighbors(2,a)
 !     !PRINT *
 !
 !   ENDDO find_neighbors
 !
 !   PRINT *, " * Nearest neighbors found. "
 !   PRINT *, " * Average number of neighbors= ", DBLE(SUM(n_neighbors))/DBLE(npart_real)
 !   PRINT *

    IF( debug ) PRINT *, "0"

    !
    !-- Print final results of the APM to file
    !
    IF( PRESENT(namefile_results) )THEN
      finalnamefile= namefile_results
    ELSE
      finalnamefile= TRIM(sph_path)//"apm_results.dat"
    ENDIF

    PRINT *, "** Printing results to file ", TRIM(namefile_results), "..."
    PRINT *

    INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

    IF( exist )THEN
        OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
              FORM= "FORMATTED", &
              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
              IOMSG= err_msg )
    ELSE
        OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "NEW", &
              FORM= "FORMATTED", &
              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening " // TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF

    IF( debug ) PRINT *, "1"

  !  IF( .NOT.ALLOCATED( nu_one ) ) ALLOCATE( nu_one( npart_real ) )
  !  IF( .NOT.ALLOCATED( particle_density_final ) ) &
  !    ALLOCATE( particle_density_final( npart_real ) )
    nu_one= one
    CALL density_loop( npart_real, pos, &    ! input
                       nu_one, h, particle_density_final )      ! output

    DO a= 1, npart_real, 1
      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        a, &
        pos(1,a), pos(2,a), pos(3,a), &
        nstar_id( a ), &
        nstar_int( a ), &
        particle_density_final( a ), &
        particle_density_final( a )*nstar_id( 1 )/particle_density_final( 1 ), &
        ABS(nstar_sph(a)-nstar_id(a))/nstar_id(a), &
        nu(a), &
        nearest_neighbors(2,a)
    ENDDO

    CLOSE( UNIT= 2 )

    IF( debug ) PRINT *, "2"

    !
    !-- Finalize output of the APM iteration
    !
    CALL reallocate_output_fields( npart_real )

    pos_input(:,1:npart_real)= pos(:,1:npart_real)
    h_output(1:npart_real)   = h(1:npart_real)
    nu_output(1:npart_real)  = nu(1:npart_real)

    npart_output= npart_real

    IF( debug ) PRINT *, "3"

    !
    !-- Deallocate global variables
    !
    IF( ALLOCATED( posmash ) ) DEALLOCATE( posmash )
    CALL deallocate_metric_on_particles()
    IF( debug ) PRINT *, "4"
    CALL deallocate_gradient()
    IF( debug ) PRINT *, "5"
    CALL deallocate_RCB_tree_memory_3D()
    IF( debug ) PRINT *, "6"
    CALL deallocate_SPH_memory()

    !
    !-- Check that there aren't multiple particles at the same position
    !
    PRINT *, "** Checking that there aren't particles with the same position..."
    PRINT *

    CALL check_particle_positions( npart_real, pos )

    !STOP



    CONTAINS



    FUNCTION validate_position_final( x, y, z ) RESULT( answer )

      !*******************************************************
      !
      !# Returns validate_position( x, y, z ) if the latter
      !  is present, 0 otherwise
      !
      !  FT 22.09.2021
      !
      !*******************************************************

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: x
      !! \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN):: y
      !! \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN):: z
      !! \(z\) coordinate of the desired point
      LOGICAL:: answer
      !! validate_position( x, y, z ) if the latter is present, 0 otherwise

      IF( PRESENT(validate_position) )THEN

        answer= validate_position( x, y, z )
        !IF( validate_position( x, y, z ) == 1 ) answer= .TRUE.
        !IF( validate_position( x, y, z ) == 0 ) answer= .FALSE.

      ELSE

        answer= .TRUE.

      ENDIF

    END FUNCTION validate_position_final


    SUBROUTINE get_nstar_id_atm &
    ( npart_real, x, y, z, nstar_sph, nstar_id, nlrf_sph, sqg, use_atmosphere )
    !, nstar_eul_id, use_atmosphere )

      !*******************************************************
      !
      !#
      !
      !
      !  FT 5.12.2021
      !
      !*******************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN):: npart_real
      !! Number of real particles (i.e., no ghost particles included here)
      DOUBLE PRECISION, INTENT(IN):: x(npart_real)
      !! Array of \(x\) coordinates
      DOUBLE PRECISION, INTENT(IN):: y(npart_real)
      !! Array of \(y\) coordinates
      DOUBLE PRECISION, INTENT(IN):: z(npart_real)
      !! Array of \(z\) coordinates
      DOUBLE PRECISION, INTENT(IN):: nstar_sph(npart)
      !! |sph| proper baryon density
      DOUBLE PRECISION, INTENT(OUT):: nstar_id(npart)
      !! Array to store the computed proper baryon number density
      DOUBLE PRECISION, INTENT(OUT):: nlrf_sph(npart)
      !# Array to store the local rest frame baryon density computed from
      !  the |sph| proper baryon density
      DOUBLE PRECISION, INTENT(OUT):: sqg(npart)
      !# Square root of minus the determinant of the sacetime metric
      !DOUBLE PRECISION, INTENT(OUT):: nstar_eul_id(npart_real)
      !# Array to store the computed proper baryon number density seen
      !  by the Eulerian observer
      LOGICAL,  INTENT( IN ):: use_atmosphere
      !# `.TRUE.` if an atmosphere should be used during the APM, to allow
      !  the real aprticles more freedom to move around and adjust;
      !  `.FALSE.` otherwise

      CALL get_nstar_id(npart_real, x, y, z, nstar_sph, nstar_id, nlrf_sph, sqg)

      IF( use_atmosphere .EQV. .TRUE. )THEN

        !$OMP PARALLEL DO DEFAULT( NONE ) &
        !$OMP             SHARED( npart_real, nstar_id, &
        !$OMP                     atmosphere_density ) &
        !$OMP             PRIVATE( a )
        DO a= 1, npart_real, 1
          IF( nstar_id(a) <= atmosphere_density )THEN
            nstar_id(a)= atmosphere_density
          ENDIF
          !IF( nstar_eul_id(a) <= atmosphere_density )THEN
          !  nstar_eul_id(a)= atmosphere_density
          !ENDIF
        ENDDO
        !$OMP END PARALLEL DO

      ELSE

        !$OMP PARALLEL DO DEFAULT( NONE ) &
        !$OMP             SHARED( npart_real, nstar_id ) &
        !$OMP             PRIVATE( a )
        DO a= 1, npart_real, 1

          IF( nstar_id( a ) < tiny_real )THEN
            PRINT *, "** ERROR! nstar_id(", a, ")=", nstar_id( a ), &
                     " in SUBROUTINE get_nstar_id_atm."
            PRINT *, " * Stopping.."
            PRINT *
            STOP
          ENDIF
          !IF( nstar_eul_id( a ) < tiny_real )THEN
          !  PRINT *, "** ERROR! nstar_eul_id(", a, ")=", nstar_eul_id( a ), &
          !           " in SUBROUTINE get_nstar_id_atm."
          !  PRINT *, " * Stopping.."
          !  PRINT *
          !  STOP
          !ENDIF

        ENDDO
        !$OMP END PARALLEL DO

      ENDIF

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( npart_real, nstar_id ) &
      !$OMP             PRIVATE( a )
      DO a= 1, npart_real, 1

        IF( .NOT.is_finite_number( nstar_id( a ) ) )THEN
          PRINT *, "** ERROR! nstar_id(", a, ")= ", nstar_id( a ), &
                   "is a not a finite number!", &
                   " in SUBROUTINE get_nstar_id_atm."
          PRINT *, " * Stopping.."
          PRINT *
          STOP
        ENDIF
        !IF( .NOT.is_finite_number( nstar_eul_id( a ) ) )THEN
        !  PRINT *, "** ERROR! nstar_eul_id(", a, ")= ", nstar_eul_id( a ), &
        !           "is a not a finite number!", &
        !           " in SUBROUTINE get_nstar_id_atm."
        !  PRINT *, " * Stopping.."
        !  PRINT *
        !  STOP
        !ENDIF

      ENDDO
      !$OMP END PARALLEL DO

    END SUBROUTINE get_nstar_id_atm


    SUBROUTINE read_pressure_id( npart_real, x, y, z, pressure_id )

      !*******************************************************
      !
      !#
      !
      !
      !  FT 6.12.2022
      !
      !*******************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN):: npart_real
      !! Number of real particles (i.e., no ghost particles included here)
      DOUBLE PRECISION, INTENT(IN):: x(npart_real)
      !! Array of \(x\) coordinates
      DOUBLE PRECISION, INTENT(IN):: y(npart_real)
      !! Array of \(y\) coordinates
      DOUBLE PRECISION, INTENT(IN):: z(npart_real)
      !! Array of \(z\) coordinates
      DOUBLE PRECISION, INTENT(OUT):: pressure_id(npart)
      !! Array to store the pressure read from the |id|

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( npart_real, pressure_id, x, y, z ) &
      !$OMP             PRIVATE( a )
      DO a= 1, npart_real, 1
        pressure_id(a)= get_pressure_id(x(a),y(a),z(a))
      ENDDO
      !$OMP END PARALLEL DO

    END SUBROUTINE read_pressure_id


    SUBROUTINE allocate_apm_fields( npart_real, npart_ghost )

      !*******************************************************
      !
      !# Allocate the fields used during the APM iteration
      !
      !  FT 20.04.2022
      !
      !*******************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN):: npart_real
      INTEGER, INTENT(IN):: npart_ghost

      INTEGER:: npart_all

      npart_all= npart_real + npart_ghost

      IF(.NOT.ALLOCATED( nstar_sph ))THEN
        ALLOCATE( nstar_sph( npart_all ), STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array nstar_sph in SUBROUTINE ", &
                    "allocate_apm_fields. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( nlrf_sph ))THEN
        ALLOCATE( nlrf_sph( npart_all ), STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array nlrf_real in SUBROUTINE ", &
                    "allocate_apm_fields. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( pressure_sph ))THEN
        ALLOCATE( pressure_sph( npart_all ), STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array pr_real in SUBROUTINE ", &
                    "allocate_apm_fields. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( sqg ))THEN
        ALLOCATE( sqg( npart_real ), STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array sqg in SUBROUTINE ", &
                    "allocate_apm_fields. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( dS ))THEN
        ALLOCATE( dS( 3, npart_real ), STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array dS in SUBROUTINE ", &
                    "allocate_apm_fields. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( nstar_id ))THEN
        ALLOCATE( nstar_id( npart_all ), STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array nstar_id in SUBROUTINE ", &
                    "allocate_apm_fields. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( nlrf_id ))THEN
        ALLOCATE( nlrf_id( npart_all ), STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array nlrf_id in SUBROUTINE ", &
                    "allocate_apm_fields. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( pressure_id ))THEN
        ALLOCATE( pressure_id( npart_all ), STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array pressure_id in SUBROUTINE ", &
                    "allocate_apm_fields. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( art_pr ))THEN
        ALLOCATE( art_pr( npart_all ), STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array art_pr in SUBROUTINE ", &
                    "allocate_apm_fields. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( correction_pos ))THEN
        ALLOCATE( correction_pos( 3, npart_all ) )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array correction_pos in ", &
                    "SUBROUTINE allocate_apm_fields. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( all_pos_prev ))THEN
        ALLOCATE( all_pos_prev( 3, npart_all ) )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array all_pos_prev in ", &
                    "SUBROUTINE allocate_apm_fields. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( cnt_move ))THEN
        ALLOCATE( cnt_move( npart_real ) )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array cnt_move in ", &
                    "SUBROUTINE allocate_apm_fields. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( dNstar ))THEN
        ALLOCATE( dNstar( npart_real ), STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array dNstar in SUBROUTINE ", &
                    "allocate_apm_fields. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( nu_tmp ))THEN
        ALLOCATE( nu_tmp( npart_real ), STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array nu_tmp in SUBROUTINE ", &
                    "allocate_apm_fields. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( pos ))THEN
        ALLOCATE( pos( 3, npart_real ), STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array pos in SUBROUTINE ", &
                    "allocate_apm_fields. The error message is", &
                    err_msg
           STOP
        ENDIF
      ENDIF
      IF( .NOT.ALLOCATED( nu_one ) )THEN
        ALLOCATE( nu_one( npart_real ), STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array nu_one in SUBROUTINE ", &
                    "allocate_apm_fields. The error message is", &
                    err_msg
           STOP
        ENDIF
      ENDIF
      IF( .NOT.ALLOCATED( particle_density_final ) )THEN
        ALLOCATE( particle_density_final( npart_real ), &
                  STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array particle_density_final in ", &
                    "SUBROUTINE allocate_apm_fields. The error message is", &
                    err_msg
           STOP
        ENDIF
      ENDIF
      IF( .NOT.ALLOCATED( cnt_array ) )THEN
        ALLOCATE( cnt_array( npart_real ), STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array cnt_array in ", &
                    "SUBROUTINE allocate_apm_fields. The error message is", &
                    err_msg
           STOP
        ENDIF
      ENDIF

    END SUBROUTINE allocate_apm_fields


    SUBROUTINE discard_atmosphere( npart_real )

      !*******************************************************
      !
      !# Remove the atmosphere
      !
      !  FT 20.04.2022
      !
      !*******************************************************

      IMPLICIT NONE

      INTEGER, INTENT(INOUT):: npart_real
      INTEGER:: a

      ALLOCATE( rho_tmp( npart_real ) )
      IF(ALLOCATED(pos_tmp)) DEALLOCATE(pos_tmp)
      ALLOCATE( pos_tmp( 3, npart_real ) )
      IF(ALLOCATED(h_tmp)) DEALLOCATE(h_tmp)
      ALLOCATE( h_tmp( npart_real ) )
      IF(ALLOCATED(h_guess_tmp)) DEALLOCATE(h_guess_tmp)
      ALLOCATE( h_guess_tmp( npart_real ) )
      IF(ALLOCATED(nu_tmp)) DEALLOCATE(nu_tmp)
      ALLOCATE( nu_tmp( npart_real ) )

      pos_tmp    = HUGE(one)
      h_tmp      = HUGE(one)
      h_guess_tmp= HUGE(one)
      nu_tmp     = HUGE(one)

      npart= 0
      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( pos, pos_tmp, h, h_tmp, rho_tmp, npart_real, &
      !$OMP                 h_guess, h_guess_tmp, nu, nu_tmp, pvol_tmp, pvol ) &
      !$OMP             PRIVATE( a ) &
      !$OMP             REDUCTION( +: npart )
      DO a= 1, npart_real, 1
        rho_tmp(a)= get_density( pos(1,a), pos(2,a), pos(3,a) )
        IF( rho_tmp(a) > zero )THEN
          npart= npart + 1
          pos_tmp(:,a)  = pos(:,a)
          h_tmp(a)      = h(a)
          h_guess_tmp(a)= h_guess(a)
          nu_tmp(a)     = nu(a)
          !pvol_tmp(a)   = pvol(a)
        ENDIF
      ENDDO
      !$OMP END PARALLEL DO

      IF(ALLOCATED(pos)) DEALLOCATE(pos)
      ALLOCATE( pos( 3, npart ) )
      IF(ALLOCATED(h)) DEALLOCATE(h)
      ALLOCATE( h( npart ) )
      IF(ALLOCATED(h_guess)) DEALLOCATE(h_guess)
      ALLOCATE( h_guess( npart ) )
      IF(ALLOCATED(nu)) DEALLOCATE(nu)
      ALLOCATE( nu( npart ) )
      !IF(ALLOCATED(pvol)) DEALLOCATE(pvol)
      !ALLOCATE( pvol( npart ) )

  !   !$OMP PARALLEL DO DEFAULT( NONE ) &
  !   !$OMP             SHARED( pos, pos_tmp, h, h_tmp, rho_tmp, npart_real, &
  !   !$OMP                     h_guess, h_guess_tmp, nu, nu_tmp ) &
  !   !$OMP             PRIVATE( a )
      cnt1= 0
      DO a= 1, npart_real, 1
        IF( h_tmp(a) < HUGE(one) )THEN
          cnt1= cnt1 + 1
          pos(:,cnt1)  = pos_tmp(:,a)
          h(cnt1)      = h_tmp(a)
          h_guess(cnt1)= h_guess_tmp(a)
          nu(cnt1)     = nu_tmp(a)
          !pvol(cnt1)   = pvol_tmp(a)
        ENDIF
      ENDDO
  !   !$OMP END PARALLEL DO

      npart_real= npart


    END SUBROUTINE discard_atmosphere


    SUBROUTINE reallocate_output_fields( npart_real )

      !*******************************************************
      !
      !# Reallocate the fields to be returned by perform_apm
      !
      !  FT 20.04.2022
      !
      !*******************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN):: npart_real

      IF( ALLOCATED( pos_input ) ) DEALLOCATE( pos_input )
      ALLOCATE( pos_input( 3, npart_real ), STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array pos_input in ", &
                 "SUBROUTINE reallocate_output_fields. The error message is", &
                 err_msg
        STOP
      ENDIF

      IF( ALLOCATED( h_output ) ) DEALLOCATE( h_output )
      ALLOCATE( h_output( npart_real ), STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array h_output in ", &
                 "SUBROUTINE reallocate_output_fields. The error message is", &
                 err_msg
        STOP
      ENDIF

      IF( ALLOCATED( nu_output ) ) DEALLOCATE( nu_output )
      ALLOCATE( nu_output( npart_real ), STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array nu_output in ", &
                 "SUBROUTINE reallocate_output_fields. The error message is", &
                 err_msg
        STOP
      ENDIF


    END SUBROUTINE reallocate_output_fields


    SUBROUTINE place_and_print_ghost_particles()

      !*******************************************************
      !
      !# Place ghost particles around the matter object,
      !  and print their positions together with the positions
      !  of the real particles to a formatted file
      !
      !  FT 20.04.2022
      !
      !*******************************************************

      USE utility, ONLY: g2kg, m2cm, lorene2hydrobase

      IMPLICIT NONE

      INTEGER:: a, a_nu_min

      DOUBLE PRECISION, PARAMETER:: eps          = 1.D0 !2.0D-1 !5.0D0
      !DOUBLE PRECISION, PARAMETER:: multiple_h_av= 1.0D0

      DOUBLE PRECISION:: nu_av, max_r_real, nstar_sph_ghost_av
      DOUBLE PRECISION, DIMENSION(npart_real):: tmp, tmp2, tmp3, &
        nstar_id_arr, nstar_sph_arr, nlrf_sph_arr, nlrf_id_arr, pressure_id_arr
      DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: ghost_pos_tmp
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: nu_ghost_arr, &
        nstar_sph_ghost, h_ghost

      IF(adapt_ghosts)THEN
        ellipse_thickness= 1.4D0
      ELSE
        ellipse_thickness= 1.1D0
      ENDIF

      npart= npart_real

      CALL allocate_SPH_memory

      CALL allocate_RCB_tree_memory_3D(npart)
      iorig(1:npart)= (/ (a,a=1,npart) /)

      IF( debug ) PRINT *, "10"

      CALL allocate_gradient( npart )
      CALL allocate_metric_on_particles( npart )

      IF(.NOT.ALLOCATED( ghost_pos ))THEN
        ALLOCATE( ghost_pos( 3, max_npart ), STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
          PRINT *, "...allocation error for array ghost_pos in SUBROUTINE ", &
                  "perform_apm. The error message is",&
                  err_msg
          STOP
        ENDIF
      ENDIF

      ghost_pos= zero

      PRINT *, " * Placing ghost particles on a lattice between ellipsodial ", &
               "surfaces..."
      PRINT *

      IF( debug ) PRINT *, "npart_real= ", npart_real

      ! Find the maximum radial coordinate over the real particles
   !   max_r_real= zero
   !   !$OMP PARALLEL DO DEFAULT( NONE ) &
   !   !$OMP             SHARED( npart_real, pos_input, center ) &
   !   !$OMP             PRIVATE( a, r_real ) &
   !   !$OMP             REDUCTION( MAX: max_r_real )
   !   DO a= 1, npart_real, 1
   !
   !     r_real= SQRT( ( pos_input(1,a) - center(1) )**two &
   !                 + ( pos_input(2,a) - center(2) )**two &
   !                 + ( pos_input(3,a) - center(3) )**two )
   !     !IF( r_real > max_r_real ) max_r_real= r_real
   !     max_r_real= MAX( max_r_real, r_real )
   !
   !   ENDDO
   !   !$OMP END PARALLEL DO

      IF(debug) PRINT *, "max_z_real =", max_z_real

   !   CALL get_nstar_id( npart_real, pos_input(1,1:npart_real), &
   !                                  pos_input(2,1:npart_real), &
   !                                  pos_input(3,1:npart_real), tmp, &
   !                                  nstar_id_arr, tmp2, tmp3 )

      ! Determine smoothing length so that each particle has exactly
      ! 300 neighbours inside 2h
      CALL assign_h( nn_des, &
                     npart_real, &
                     pos_input(:,1:npart_real), h_guess(1:npart_real), &
                     h(1:npart_real) )

      find_h_bruteforce_timer= timer( "find_h_bruteforce_timer" )
      CALL find_h_bruteforce_timer% start_timer()
      n_problematic_h= 0
      check_h1: DO a= 1, npart_real, 1
      ! find_h_backup, called below, is OMP parallelized, so this loop
      ! should not be parallelized as well

        IF( .NOT.is_finite_number( h(a) ) .OR. h(a) <= zero )THEN

          n_problematic_h= n_problematic_h + 1
          h(a)= find_h_backup( a, npart_real, pos_input(:,1:npart_real), nn_des)
          IF( .NOT.is_finite_number( h(a) ) .OR. h(a) <= zero )THEN
            PRINT *, "** ERROR! h=0 on particle ", a, "even with the brute", &
                     " force method."
            PRINT *, "   Particle position: ", pos_input(:,a)
            STOP
          ENDIF

        ENDIF

      ENDDO check_h1
      CALL find_h_bruteforce_timer% stop_timer()
      CALL find_h_bruteforce_timer% print_timer( 2 )

      PRINT *, " * The smoothing length was found brute-force for ", &
               n_problematic_h, " particles."
      PRINT *

      PRINT *, " * Compute SPH density..."
      PRINT *

      CALL density_loop( npart_real, pos_input(:,1:npart_real), &
        nu_output(1:npart_real), h(1:npart_real), nstar_sph_arr(1:npart_real) )

      PRINT *, " * Read ID density..."
      PRINT *

      CALL get_nstar_id( npart_real, pos_input(1,1:npart_real), &
                                     pos_input(2,1:npart_real), &
                                     pos_input(3,1:npart_real), &
                                     nstar_sph_arr, &
                                     nstar_id_arr, nlrf_sph_arr, tmp3 )

      IF(use_pressure)THEN

        !$OMP PARALLEL DO DEFAULT( NONE ) &
        !$OMP             SHARED( npart_real, nlrf_id_arr, nstar_id_arr, &
        !$OMP                     nstar_sph_arr, nlrf_sph_arr ) &
        !$OMP             PRIVATE( a )
        DO a= 1, npart_real, 1
          nlrf_id_arr(a)= nstar_id_arr(a)*nlrf_sph_arr(a)/nstar_sph_arr(a)
        ENDDO
        !$OMP END PARALLEL DO

        CALL compute_pressure( npart_real, &
                               pos_input(1,1:npart_real), &
                               pos_input(2,1:npart_real), &
                               pos_input(3,1:npart_real), &
                               nlrf_id_arr, eqos, pressure_id_arr )

      ENDIF

      itr           = 0
      nu_av         = zero
      nstar_id_av   = zero
      nlrf_id_av    = zero
      nstar_sph_av  = zero
      pressure_id_av= zero
      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( npart_real, nlrf_id_arr, nstar_id_arr, &
      !$OMP                     nstar_sph_arr, pos_input, pressure_id_arr, &
      !$OMP                     max_z_real, center, nu_output ) &
      !$OMP             PRIVATE( a ) &
      !$OMP             REDUCTION(+: itr, nu_av, nstar_id_av, nlrf_id_av, &
      !$OMP                          nstar_sph_av, pressure_id_av )
      DO a= 1, npart_real, 1

        IF( SQRT( ( pos_input(1,a) - center(1) )**two &
                + ( pos_input(2,a) - center(2) )**two &
                + ( pos_input(3,a) - center(3) )**two ) &
                  > (one - one/(ten*ten))*max_z_real )THEN

          itr           = itr + 1
          nu_av         = nu_av          + nu_output(a)
          nstar_id_av   = nstar_id_av    + nstar_id_arr(a)
          nlrf_id_av    = nlrf_id_av     + nlrf_id_arr(a)
          nstar_sph_av  = nstar_sph_av   + nstar_sph_arr(a)
          pressure_id_av= pressure_id_av + pressure_id_arr(a)

        ENDIF

      ENDDO
      !$OMP END PARALLEL DO
      nu_av         = nu_av/itr
      nstar_id_av   = nstar_id_av/itr
      nlrf_id_av    = nlrf_id_av/itr
      nstar_sph_av  = nstar_sph_av/itr
      pressure_id_av= pressure_id_av/itr

      xmin= center(1) - sizes(1)*( one + eps )
      xmax= center(1) + sizes(2)*( one + eps )
      ymin= center(2) - sizes(3)*( one + eps )
      ymax= center(2) + sizes(4)*( one + eps )
      zmin= center(3) - sizes(5)*( one + eps )
      zmax= center(3) + sizes(6)*( one + eps )

      IF(adapt_ghosts)THEN

        !
        !-- Determine dx,dy,dz from desired density from the ID, and compute
        !-- nx,ny,nz,ghost_dist,nu_ghost from them
        !
        dx        = ( nu_av/nstar_id_av )**third
        dy        = dx
        dz        = dx
        ghost_dist= dx
        nu_ghost  = dx*dy*dz*nstar_id_av
        nx= NINT( ABS( xmax - xmin )/dx )
        ny= NINT( ABS( ymax - ymin )/dy )
        nz= NINT( ABS( zmax - zmin )/dz )

        IF(debug) PRINT*, "# particles over which the averages are taken=", itr
        IF(debug) PRINT*, "nstar_id_av      =", nstar_id_av
        IF(debug) PRINT*, "nstar_sph_av     =", nstar_sph_av
        IF(debug) PRINT*, "(1.D+12g/cm**3)  =", &
          1.D+12*(g2kg*m2cm**3)*lorene2hydrobase*umass/amu
        IF(debug) PRINT*, "nu_all           =",(mass/DBLE(npart_real))*umass/amu
        IF(debug) PRINT*, "nu_av            =", nu_av
        IF(debug) PRINT*, "nu_av/nstar_id_av=", nu_av/nstar_id_av
        !STOP

      ELSE

        !
        !-- Read nx,ny,nz from parameter file, and compute dx,dy,dz from them.
        !-- The ghosts have the same mass as the real particles
        !
        nx= nx_gh
        ny= ny_gh
        nz= nz_gh
        dx= ABS( xmax - xmin )/DBLE( nx )
        dy= ABS( ymax - ymin )/DBLE( ny )
        dz= ABS( zmax - zmin )/DBLE( nz )
        nu_ghost= (mass/DBLE(npart_real))*umass/amu

      ENDIF

      IF(debug) PRINT *, "dx =", dx
      IF(debug) PRINT *, "dy =", dy
      IF(debug) PRINT *, "dz =", dz
      IF(debug) PRINT *, "nx =", nx
      IF(debug) PRINT *, "ny =", ny
      IF(debug) PRINT *, "nz =", nz
      IF(debug) PRINT *, "nu_ghost =", nu_ghost

      IF(.NOT.ALLOCATED( ghost_pos_tmp ))THEN
        ALLOCATE( ghost_pos_tmp( 3, nx, ny, nz ), STAT= ios, &
            ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array ghost_pos in SUBROUTINE ", &
                    "perform_apm. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF

      rad_x= radius_x + ghost_dist !+ multiple_h_av*h_av
      rad_y= radius_y + ghost_dist !+ multiple_h_av*h_av
      rad_z= radius_z + ghost_dist !+ multiple_h_av*h_av

      IF(adapt_ghosts)THEN

        rad_x= ABS( MAXVAL( ABS(pos_input(1,:) - center(1)), DIM= 1 ) )
        rad_y= ABS( MAXVAL( ABS(pos_input(2,:) - center(2)), DIM= 1 ) )
        rad_z= ABS( MAXVAL( ABS(pos_input(3,:) - center(3)), DIM= 1 ) )

      ENDIF

      PRINT *, "** Distance between the size of the object and the ghost ", &
               "particles: ghost_dist =", ghost_dist
      PRINT *

      IF( debug ) PRINT *, "radius_x= ", radius_x
      IF( debug ) PRINT *, "radius_y= ", radius_y
      IF( debug ) PRINT *, "radius_z= ", radius_z
      IF( debug ) PRINT *, "rad_x= ", rad_x
      IF( debug ) PRINT *, "rad_y= ", rad_y
      IF( debug ) PRINT *, "rad_z= ", rad_z
      IF( debug ) PRINT *

      ghost_pos_tmp= HUGE(zero)

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( nx, ny, nz, xmin, ymin, zmin, dx, dy, dz, &
      !$OMP                     ghost_pos_tmp, ellipse_thickness, &
      !$OMP                     center, rad_x, rad_y, rad_z ) &
      !$OMP             PRIVATE( i, j, k, xtemp, ytemp, ztemp, &
      !$OMP                      x_ell, y_ell, z_ell, r, theta, phi, &
      !$OMP                      r_ell, theta_ell, phi_ell )
      DO k= 1, nz, 1

        ztemp= zmin + dz/two + DBLE(k - 1)*dz

        DO j= 1, ny, 1

          ytemp= ymin + dy/two + DBLE(j - 1)*dy

          DO i= 1, nx, 1

            xtemp= xmin + dx/two + DBLE(i - 1)*dx

            ! Compute te spherical polar coordinates of the point, to get
            ! its angular coordinates
            CALL spherical_from_cartesian( xtemp, ytemp, ztemp, &
                                           center(1), center(2), center(3), &
                                           r, theta, phi )

            !x_ell= center(1) + rad_x*SIN(theta)*COS(phi)
            !y_ell= center(2) + rad_y*SIN(theta)*SIN(phi)
            !z_ell= center(3) + rad_z*COS(theta)

            ! Use the angular coordinates of the point and the ellipse semiaxes,
            ! to obtain the Cartesian coordinates of the point on the ellipsoid
            ! having the same angular cordinates of the input point
            CALL cartesian_from_spherical( rad_x, theta, phi, &
              center(1), center(2), center(3), &
              x_ell, y_ell, z_ell, rad_y/rad_x, rad_z/rad_x )

            ! Compute the spherical polar coordinates of the point on the
            ! ellipsoid
            CALL spherical_from_cartesian( x_ell, y_ell, z_ell, &
                                           center(1), center(2), center(3), &
                                           r_ell, theta_ell, phi_ell )

            ! Place a ghost particle if: (i) its radial coordinate is larger
            ! than the radial coordinate the ellipsoid and lower than
            ! ellipse_thickness times it; (ii) the density is <=0
            IF( ( r <= ellipse_thickness*r_ell .AND. r >= r_ell &
                  .AND. &
                  get_density(xtemp, ytemp, ztemp) <= zero ) &
            )THEN

              ghost_pos_tmp(1, i, j, k)= xtemp
              ghost_pos_tmp(2, i, j, k)= ytemp
              ghost_pos_tmp(3, i, j, k)= ztemp

            ENDIF

           ENDDO
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

      a= 0
      DO k= 1, nz, 1

        DO j= 1, ny, 1

          DO i= 1, nx, 1

            IF( ghost_pos_tmp(1, i, j, k) < HUGE(zero) )THEN

              a= a + 1
              ghost_pos(1,a)= ghost_pos_tmp( 1, i, j, k )
              ghost_pos(2,a)= ghost_pos_tmp( 2, i, j, k )
              ghost_pos(3,a)= ghost_pos_tmp( 3, i, j, k )

            ENDIF

           ENDDO
        ENDDO
      ENDDO
      npart_ghost= a
      IF( npart_ghost == 0 )THEN
        PRINT *, "** ERROR: No ghost particles were placed. Is the ", &
                 "PARAMETER 'ghost_dist' appropriate for the physical system?"
        PRINT *, " * Stopping.."
        PRINT *
        STOP
      ENDIF
      IF( npart_ghost > max_npart )THEN
        PRINT *
        PRINT *, "** ERROR! Too many ghost particles placed!"
        PRINT *, " * npart_ghost= ", npart_ghost
        PRINT *, " * npart_ghost should not be larger than max_npart= ", &
                 max_npart
        PRINT *, " * How are the ghost particles placed? Does the algorithm ", &
                 "make reasonable sense? Perhaps try to set a smaller ", &
                 "parameter eps."
        PRINT *, " * Stopping..."
        PRINT *
        STOP
      ENDIF
      IF( npart_ghost > max_npart )THEN
        PRINT *
        PRINT *, "** ERROR! Too many ghost particles placed!"
        PRINT *, " * npart_ghost= ", npart_ghost
        PRINT *, " * npart_ghost should not be larger than max_npart= ", &
                 max_npart
        PRINT *, " * How are the ghost particles placed? Does the algorithm ", &
                 "make reasonable sense? Perhaps try to set a smaller ", &
                 "parameter eps."
        PRINT *, " * Stopping..."
        PRINT *
        STOP
      ENDIF
      ghost_pos = ghost_pos( :, 1:npart_ghost )

      PRINT *, " * ", npart_ghost, " ghost particles placed around ", &
               npart_real, "real particles."
      PRINT *

      DEALLOCATE( ghost_pos_tmp )

      npart_all= npart_real + npart_ghost

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( h_guess, npart_real, npart_all, dx, dy, dz ) &
      !$OMP             PRIVATE( a )
      DO a= npart_real + 1, npart_all, 1
        h_guess(a)= ( dx*dy*dz )**third
      ENDDO
      !$OMP END PARALLEL DO

      CALL deallocate_metric_on_particles()
      CALL deallocate_gradient()
      CALL deallocate_RCB_tree_memory_3D()
      CALL deallocate_SPH_memory()

      IF(debug)THEN
        !
        !-- Estimate ghost density
        !

        npart= npart_ghost

        CALL allocate_SPH_memory

        CALL allocate_RCB_tree_memory_3D(npart)
        iorig(1:npart)= (/ (a,a=1,npart) /)

        IF( debug ) PRINT *, "10"

        CALL allocate_gradient( npart )
        CALL allocate_metric_on_particles( npart )

        ALLOCATE( nu_ghost_arr(npart_ghost) )
        ALLOCATE( nstar_sph_ghost(npart_ghost) )
        ALLOCATE( h_ghost(npart_ghost) )

        ! Determine smoothing length so that each particle has exactly
        ! 300 neighbours inside 2h
        CALL assign_h( nn_des, &
                       npart_ghost, &
                       ghost_pos(:,1:npart_real), &
                       h_guess(npart_real+1:npart_all), &
                       h_ghost )

        find_h_bruteforce_timer= timer( "find_h_bruteforce_timer" )
        CALL find_h_bruteforce_timer% start_timer()
        n_problematic_h= 0
        DO a= 1, npart_ghost, 1
        ! find_h_backup, called below, is OMP parallelized, so this loop
        ! should not be parallelized as well

          IF( .NOT.is_finite_number( h(a) ) .OR. h_ghost(a) <= zero )THEN

            n_problematic_h= n_problematic_h + 1
            h_ghost(a)= find_h_backup( a, npart_ghost, ghost_pos(:,1:npart_ghost), nn_des)
            IF( .NOT.is_finite_number( h_ghost(a) ) .OR. h_ghost(a) <= zero )THEN
              PRINT *, "** ERROR! h=0 on particle ", a, "even with the brute", &
                       " force method."
              PRINT *, "   Particle position: ", ghost_pos(:,a)
              STOP
            ENDIF

          ENDIF

        ENDDO
        CALL find_h_bruteforce_timer% stop_timer()
        CALL find_h_bruteforce_timer% print_timer( 2 )

        PRINT *, " * The smoothing length was found brute-force for ", &
                 n_problematic_h, " particles."
        PRINT *

        PRINT *, " * Measure SPH particle number density..."
        PRINT *

        nu_ghost_arr= nu_ghost
        CALL density_loop( npart_ghost, ghost_pos, &
                           nu_ghost_arr, h_ghost, nstar_sph_ghost )

        !nstar_sph_ghost_av= SUM(nstar_sph_ghost)/npart_ghost
        PRINT *, "nstar_sph_ghost_max=", MAXVAL(nstar_sph_ghost, DIM=1)
        PRINT *, "nstar_sph_ghost_av=", SUM(nstar_sph_ghost)/npart_ghost
        PRINT *, "nstar_sph_ghost_min=", MINVAL(nstar_sph_ghost, DIM=1)

        CALL deallocate_metric_on_particles()
        CALL deallocate_gradient()
        CALL deallocate_RCB_tree_memory_3D()
        CALL deallocate_SPH_memory()

        npart= npart_real

      ENDIF

      !--------------------------------------------------!
      !--  Printing ghost particles to formatted file  --!
      !--------------------------------------------------!

      PRINT *, " * Printing ghost particles to file..."

      IF( PRESENT(namefile_pos_id) )THEN
        finalnamefile= namefile_pos_id
      ELSE
        finalnamefile= TRIM(sph_path)//"apm_pos_id.dat"
      ENDIF

      INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

      IF( exist )THEN
        OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
              FORM= "FORMATTED", &
              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
              IOMSG= err_msg )
      ELSE
        OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "NEW", &
              FORM= "FORMATTED", &
              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
      ENDIF
      IF( ios > 0 )THEN
        PRINT *, "...error when opening " // TRIM(finalnamefile), &
                 ". The error message is", err_msg
        STOP
      ENDIF

      DO a= 1, npart_real, 1
        WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
          1, a, &
          pos_input(1,a), &
          pos_input(2,a), &
          pos_input(3,a)
      ENDDO

      DO a= 1, npart_ghost, 1
        WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
          2, a, &
          ghost_pos(1,a), &
          ghost_pos(2,a), &
          ghost_pos(3,a)
      ENDDO

      CLOSE( UNIT= 2 )

      PRINT *, " * Positions of ghost and real particles printed to ", &
               TRIM(finalnamefile), " ."

      !STOP


    END SUBROUTINE place_and_print_ghost_particles


    SUBROUTINE dump_apm_pos()

      !*******************************************************
      !
      !# Print the positions of the real and ghost particles
      !  to a formatted file
      !
      !  FT 20.04.2022
      !
      !*******************************************************

      IMPLICIT NONE

      INTEGER, PARAMETER:: unit_dump= 7314891

      IF( debug ) PRINT *, "printing positions to file..."

      IF( PRESENT(namefile_pos) )THEN
        finalnamefile= namefile_pos
      ELSE
        finalnamefile= TRIM(sph_path)//"apm_pos.dat"
      ENDIF

      INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

      IF( exist )THEN
          OPEN( UNIT= unit_dump, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
                FORM= "FORMATTED", &
                POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
                IOMSG= err_msg )
      ELSE
          OPEN( UNIT= unit_dump, FILE= TRIM(finalnamefile), STATUS= "NEW", &
                FORM= "FORMATTED", &
                ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
      ENDIF
      IF( ios > 0 )THEN
        PRINT *, "...error when opening " // TRIM(finalnamefile), &
                 ". The error message is", err_msg
        STOP
      ENDIF

      DO a= 1, npart_real, 1
        tmp= get_density( all_pos(1,a), all_pos(2,a), all_pos(3,a) )
        WRITE( UNIT = unit_dump, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
          1, a, &
          all_pos(1,a), &
          all_pos(2,a), &
          all_pos(3,a), &
          nu_tmp(a), &
          tmp, cnt_move(a)
      ENDDO

      DO a= npart_real + 1, npart_all, 1
        tmp= get_density( all_pos(1,a), all_pos(2,a), all_pos(3,a) )
        WRITE( UNIT = unit_dump, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
          2, a, &
          all_pos(1,a), &
          all_pos(2,a), &
          all_pos(3,a), &
          tmp, &
          art_pr(a), nstar_sph(a), nstar_id_av
      ENDDO

      CLOSE( UNIT= unit_dump )


    END SUBROUTINE dump_apm_pos


    SUBROUTINE compute_hydro_momentum()

      !*******************************************************
      !
      !#
      !
      !  FT 17.06.2022
      !
      !*******************************************************


      USE RCB_tree_3D,         ONLY: iorig, nic, nfinal, nprev, lpart, &
                                     rpart

      IMPLICIT NONE

      DOUBLE PRECISION:: ha, ha_1, ha_3, ha_4, va, xa, ya, za, &
                         hb, hb_1, hb_3, hb_4, xb, yb, zb, rab, rab2, rab_1, &
                         ha2, ha2_4, hb2_4
      !DOUBLE PRECISION:: mat_xx, mat_xy, mat_xz, mat_yy, mat_yz, mat_zz
      DOUBLE PRECISION:: &!Wdx, Wdy, Wdz, ddx, ddy, ddz, Wab, &
                         Wab_ha, Wi, Wi1, dvv, &
                         grW_ha_x, grW_ha_y, grW_ha_z, &
                         grW_hb_x, grW_hb_y, grW_hb_z, eab(3), &
                         Wa, grW, grWa, grWb, vb, prgNa, prgNb


  !    !$OMP PARALLEL DO DEFAULT( NONE ) &
  !    !$OMP             SHARED( this, nlrf_sph, pressure_sph, npart_real, &
  !    !$OMP                     m0c2_cu ) &
  !    !$OMP             PRIVATE( a )
  !    compute_pressure: DO a= 1, npart_real, 1
  !
  !      pressure_sph(a)= this% all_eos(1)% eos_parameters(poly$kappa)/m0c2_cu &
  !                 *(nlrf_sph(a)*m0c2_cu) &
  !                 **this% all_eos(1)% eos_parameters(poly$gamma)
  !
  !    ENDDO compute_pressure
  !    !$OMP END PARALLEL DO

      CALL compute_pressure( npart_real, &
                             all_pos(1,1:npart_real), &
                             all_pos(2,1:npart_real), &
                             all_pos(3,1:npart_real), &
                             nlrf_sph(1:npart_real), eqos, &
                             pressure_sph(1:npart_real) )

      !PRINT *, "Before calling ll_cell_loop..."
      !PRINT *

      dS= zero
      dS_norm_av= zero
      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( nfinal, nprev, iorig, lpart, rpart, nic, &
      !$OMP                     ncand, h, all_pos, all_clists, nlrf_sph, &
      !$OMP                     pressure_sph, dS, sqg, nstar_sph, nu_tmp, &
      !$OMP                     m0c2_cu ) &
      !$OMP             PRIVATE( ill, itot, a, b, l, &
      !$OMP                      ha, ha_1, ha_3, ha_4, ha2, ha2_4, &
      !$OMP                      hb, hb_1, hb_3, hb_4, hb2_4, &
      !$OMP                      xa, ya, za, xb, yb, zb, dx, dy, dz, rab2, &
      !$OMP                      inde, va, dv_table_1, index1, Wi, W_no_norm, &
      !$OMP                      dvv, dv_table, Wab_ha, Wa, grW, Wi1, &
      !$OMP                      grW_ha_x, grW_ha_y, grW_ha_z, grWa, grWb, &
      !$OMP                      grW_hb_x, grW_hb_y, grW_hb_z, eab, rab, vb, &
      !$OMP                      prgNa, prgNb, rab_1 )
      ll_cell_loop2: DO ill= 1, nfinal

        itot= nprev + ill
        IF( nic(itot) == 0 ) CYCLE

        particle_loop2: DO l= lpart(itot), rpart(itot)

          a=         iorig(l)

          ha=        h(a)
          ha_1=      one/ha
          ha_3=      ha_1*ha_1*ha_1
          ha_4=      ha_3*ha_1
          ha2=       ha*ha
          ha2_4=     two*two*ha2

          xa=        all_pos(1,a)
          ya=        all_pos(2,a)
          za=        all_pos(3,a)

          prgNa= (pressure_sph(a)*sqg(a)/(nstar_sph(a)/m0c2_cu)**2)

          cand_loop2: DO k= 1, ncand(ill)

            b= all_clists(ill)%list(k)

            hb=   h(b)
            hb_1= one/hb
            hb_3= hb_1*hb_1*hb_1
            hb_4= hb_3*hb_1
            hb2_4=     two*two*hb*hb

            xb=   all_pos(1,b)  ! CONTRA-variant
            yb=   all_pos(2,b)
            zb=   all_pos(3,b)

            ! potentially bail out
            dx= xa  - xb
            dy= ya  - yb
            dz= za  - zb

            rab2 = dx*dx + dy*dy + dz*dz
            rab  = SQRT(rab2) + 1.D-30
            rab_1= 1.D0/rab
            IF( rab2 > ha2_4 .AND. rab2 > hb2_4 ) CYCLE

            !--------!
            !-- ha --!
            !--------!
            va= rab*ha_1

            ! get interpolation indices
            inde  = MIN(INT(va*dv_table_1),n_tab_entry)
            index1= MIN(inde + 1,n_tab_entry)

            ! get tabulated values
            Wi=     W_no_norm(inde)
            Wi1=    W_no_norm(index1)

            ! interpolate
            dvv=    (va - DBLE(inde)*dv_table)*dv_table_1
            Wab_ha= Wi + (Wi1 - Wi)*dvv

            ! unit vector ("a-b" --> from b to a)
            eab(1)=    dx*rab_1
            eab(2)=    dy*rab_1
            eab(3)=    dz*rab_1

            ! kernel and its gradient
            !DIR$ INLINE
            CALL interp_W_gradW_table( va, Wa, grW )
            Wa=       Wa*ha_3
            grW=      grW*ha_4

            ! nabla_a Wab(ha)
            grW_ha_x= grW*eab(1)
            grW_ha_y= grW*eab(2)
            grW_ha_z= grW*eab(3)

            grWa=     grW_ha_x*eab(1) + &
                      grW_ha_y*eab(2) + &
                      grW_ha_z*eab(3)

            !--------!
            !-- hb --!
            !--------!
            vb=       rab*hb_1

            ! kernel and its gradient
            !DIR$ INLINE
            CALL interp_gradW_table(vb,grW)
            grW=      grW*hb_4

            ! nabla_a Wab(hb)
            grW_hb_x= grW*eab(1)
            grW_hb_y= grW*eab(2)
            grW_hb_z= grW*eab(3)

            grWb= grW_hb_x*eab(1) + &
                  grW_hb_y*eab(2) + &
                  grW_hb_z*eab(3)

            prgNb= pressure_sph(b)*sqg(b)/((nstar_sph(b)/m0c2_cu)**2)

            dS(1,a)= dS(1,a) - nu_tmp(b)*( prgNa*grW_ha_x + prgNb*grW_hb_x )
            dS(2,a)= dS(2,a) - nu_tmp(b)*( prgNa*grW_ha_y + prgNb*grW_hb_y )
            dS(3,a)= dS(3,a) - nu_tmp(b)*( prgNa*grW_ha_z + prgNb*grW_hb_z )

          ENDDO cand_loop2

        ENDDO particle_loop2

      ENDDO ll_cell_loop2
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( dS, npart_real ) &
      !$OMP             PRIVATE( a ) &
      !$OMP             REDUCTION(+: dS_norm_av)
      particle_loop3: DO a= 1, npart_real, 1

         dS_norm_av= dS_norm_av + SQRT(dS(1,a)**2 + dS(2,a)**2 + dS(3,a)**2)

      END DO particle_loop3
      !$OMP END PARALLEL DO
      dS_norm_av= dS_norm_av/npart_real

      PRINT *, " * The average over the particles, of the norm of the", &
               " canonical momentum due to the particle distribution only is:",&
               dS_norm_av
      PRINT *

      !PRINT *, "After calling ll_cell_loop..."
      !PRINT *

      INQUIRE( FILE= TRIM("momentum.dat"), EXIST= exist )

      IF( exist )THEN
          OPEN( UNIT= 23, FILE= TRIM("momentum.dat"), STATUS= "REPLACE", &
                FORM= "FORMATTED", &
                POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
                IOMSG= err_msg )
      ELSE
          OPEN( UNIT= 23, FILE= TRIM("momentum.dat"), STATUS= "NEW", &
                FORM= "FORMATTED", &
                ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
      ENDIF
      IF( ios > 0 )THEN
        PRINT *, "...error when opening " // TRIM("momentum.dat"), &
                 ". The error message is", err_msg
        STOP
      ENDIF

      DO a= 1, npart_real, 1
        WRITE( UNIT = 23, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
          a, &
          all_pos(1,a), &
          all_pos(2,a), &
          all_pos(3,a), &
          dS(1,a), &
          dS(2,a), &
          dS(3,a), &
          pressure_sph(a)!, &
          !nu_tmp(a), &
          !nstar_sph(a)/m0c2_cu
      ENDDO

      CLOSE( UNIT= 23 )


      !STOP


    END SUBROUTINE compute_hydro_momentum


  END PROCEDURE perform_apm

  
END SUBMODULE apm
