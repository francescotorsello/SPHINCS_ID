!& File:         submodule_sph_particles_apm.f90
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


  USE constants,            ONLY: zero, quarter, one, two, three, ten
  USE utility,              ONLY: is_finite_number
  USE sph_variables,        ONLY: allocate_sph_memory, deallocate_sph_memory, &
                                  npart
  USE metric_on_particles,  ONLY: allocate_metric_on_particles, &
                                  deallocate_metric_on_particles
  USE gradient,             ONLY: allocate_gradient, deallocate_gradient

  USE RCB_tree_3D,          ONLY: allocate_RCB_tree_memory_3D, &
                                  deallocate_RCB_tree_memory_3D, iorig
  USE set_h,                ONLY: posmash


  IMPLICIT NONE


  TYPE apm_fields
  !! Data structure containing the fields used during the APM iteration


    INTEGER, DIMENSION(:), ALLOCATABLE:: cnt_move
    !# Array storing \(0\) if a particle does not move at an APM step,
    !  \(1\) if it moves

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: lapse,                      &
                                                  shift_x, shift_y, shift_z,  &
                                                  g_xx, g_xy, g_xz,           &
                                                  g_yy, g_yz, g_zz,           &
                                                  baryon_density,             &
                                                  energy_density,             &
                                                  specific_energy,            &
                                                  pressure,                   &
                                                  v_euler_x,v_euler_y,v_euler_z

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: all_pos
    !! Array storing the real particles and the ghost particles, in this order
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: all_pos_prev_dump
    !# Array storing the real particles and the ghost particles, in this order,
    !  from the previous dump (for plotting purposes only)
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: correction_pos
    !# Array storing the correction to the positions of the real particles and
    !  the ghost particles. Thoe for the ghost particles are zero

    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: h_guess
    !! First guesses for the smoothing lengths

    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar_id
    !! Proper baryon number density computed from the |id|
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar_sph
    !! |sph| estimate of the proper baryon number density
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: dNstar
    !# Relative difference between the proper baryon number density computed
    !  from the |id| and its |sph| estimate
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: art_pr
    !! Artificial pressure
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nu_tmp
    !! Cache-array used to temporarily store the baryon number nu

    LOGICAL, DIMENSION(:), ALLOCATABLE:: good_rho
    !! Array storing `.TRUE.` if the density is acceptable, `.FALSE.` otherwise


    CONTAINS


    PROCEDURE:: allocate_apm_fields
    !! Allocates the APM fields member of TYPE [[apm_fields]]

    PROCEDURE:: deallocate_apm_fields
    !! Deallocates the APM fields member of TYPE [[apm_fields]]

    PROCEDURE:: reallocate_apm_fields
    !! Reallocates the APM fields member of TYPE [[apm_fields]]


  END TYPE


  ! TODO: VERY IMPORTANT! These fields here are global!
  !       Maybe make an object called 'apm_fields' and include them as members,
  !       so you handle their allocation and deallocatio locally and SAFELY
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: lapse, &
                                             shift_x, shift_y, shift_z, &
                                             g_xx, g_xy, g_xz, &
                                             g_yy, g_yz, g_zz, &
                                             baryon_density, &
                                             energy_density, &
                                             specific_energy, &
                                             pressure, &
                                             v_euler_x, v_euler_y, v_euler_z
  LOGICAL, DIMENSION(:),   ALLOCATABLE:: good_rho

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: all_pos
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: all_pos_prev_dump
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: correction_pos

  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: h_guess

  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar_p
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar_real
  !DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: dN
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: dNstar
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: art_pr
  !DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: freeze
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nu_tmp
  INTEGER, DIMENSION(:), ALLOCATABLE:: cnt_move


  CONTAINS


  SUBROUTINE deallocate_apm_fields( this )

    !*****************************************************
    !
    !# Deallocates the APM fields member of TYPE [[apm_fields]]
    !
    !  FT 29.03.2022
    !
    !*****************************************************

    IMPLICIT NONE

    CLASS(apm_fields), INTENT(INOUT):: this
    !! The [[apm_fields]] object that this PROCEDURE is a member of

    IF( ALLOCATED( this% nstar_sph )         ) DEALLOCATE( this% nstar_sph )
    IF( ALLOCATED( this% nstar_id )          ) DEALLOCATE( this% nstar_id )
    IF( ALLOCATED( this% art_pr )            ) DEALLOCATE( this% art_pr )
    IF( ALLOCATED( this% dNstar )            ) DEALLOCATE( this% dNstar )
    IF( ALLOCATED( this% nu_tmp )            ) DEALLOCATE( this% nu_tmp )
    IF( ALLOCATED( this% good_rho )          ) DEALLOCATE( this% good_rho )
    IF( ALLOCATED( this% correction_pos )    ) &
      DEALLOCATE( this% correction_pos )
    IF( ALLOCATED( this% all_pos_prev_dump ) ) &
      DEALLOCATE( this% all_pos_prev_dump )
    IF( ALLOCATED( this% cnt_move )          ) DEALLOCATE( this% cnt_move )

    IF( ALLOCATED( this% lapse            )) DEALLOCATE( this% lapse           )
    IF( ALLOCATED( this% shift_x          )) DEALLOCATE( this% shift_x         )
    IF( ALLOCATED( this% shift_y          )) DEALLOCATE( this% shift_y         )
    IF( ALLOCATED( this% shift_z          )) DEALLOCATE( this% shift_z         )
    IF( ALLOCATED( this% g_xx             )) DEALLOCATE( this% g_xx            )
    IF( ALLOCATED( this% g_xy             )) DEALLOCATE( this% g_xy            )
    IF( ALLOCATED( this% g_xz             )) DEALLOCATE( this% g_xz            )
    IF( ALLOCATED( this% g_yy             )) DEALLOCATE( this% g_yy            )
    IF( ALLOCATED( this% g_yz             )) DEALLOCATE( this% g_yz            )
    IF( ALLOCATED( this% g_zz             )) DEALLOCATE( this% g_zz            )
    IF( ALLOCATED( this% baryon_density   )) DEALLOCATE( this% baryon_density  )
    IF( ALLOCATED( this% energy_density   )) DEALLOCATE( this% energy_density  )
    IF( ALLOCATED( this% specific_energy  )) DEALLOCATE( this% specific_energy )
    IF( ALLOCATED( this% pressure         )) DEALLOCATE( this% pressure        )
    IF( ALLOCATED( this% v_euler_x        )) DEALLOCATE( this% v_euler_x       )
    IF( ALLOCATED( this% v_euler_y        )) DEALLOCATE( this% v_euler_y       )
    IF( ALLOCATED( this% v_euler_z        )) DEALLOCATE( this% v_euler_z       )

  END SUBROUTINE deallocate_apm_fields


  SUBROUTINE allocate_apm_fields( this, npart_real, npart_ghost )

    !*****************************************************
    !
    !# Allocates the APM fields member of TYPE [[apm_fields]]
    !
    !  FT 29.03.2022
    !
    !*****************************************************

    IMPLICIT NONE

    CLASS(apm_fields), INTENT(INOUT):: this
    !! The [[apm_fields]] object that this PROCEDURE is a member of
    INTEGER, INTENT(IN):: npart_real
    !! The number of real particles
    INTEGER, INTENT(IN):: npart_ghost
    !! The number of ghost particles

    ALLOCATE( this% nstar_sph        (npart_real + npart_ghost) )
    ALLOCATE( this% nstar_id         (npart_real + npart_ghost) )
    ALLOCATE( this% art_pr           (npart_real + npart_ghost) )
    ALLOCATE( this% correction_pos   (3,npart_real + npart_ghost) )
    ALLOCATE( this% all_pos_prev_dump(3,npart_real + npart_ghost) )

    ALLOCATE( this% cnt_move         (npart_real) )
    ALLOCATE( this% dNstar           (npart_real) )
    ALLOCATE( this% nu_tmp           (npart_real) )
    ALLOCATE( this% good_rho         (npart_real) )

    ALLOCATE( this% lapse            (npart_real) )
    ALLOCATE( this% shift_x          (npart_real) )
    ALLOCATE( this% shift_y          (npart_real) )
    ALLOCATE( this% shift_z          (npart_real) )
    ALLOCATE( this% g_xx             (npart_real) )
    ALLOCATE( this% g_xy             (npart_real) )
    ALLOCATE( this% g_xz             (npart_real) )
    ALLOCATE( this% g_yy             (npart_real) )
    ALLOCATE( this% g_yz             (npart_real) )
    ALLOCATE( this% g_zz             (npart_real) )
    ALLOCATE( this% baryon_density   (npart_real) )
    ALLOCATE( this% energy_density   (npart_real) )
    ALLOCATE( this% specific_energy  (npart_real) )
    ALLOCATE( this% pressure         (npart_real) )
    ALLOCATE( this% v_euler_x        (npart_real) )
    ALLOCATE( this% v_euler_y        (npart_real) )
    ALLOCATE( this% v_euler_z        (npart_real) )

  END SUBROUTINE allocate_apm_fields


  SUBROUTINE reallocate_apm_fields( this, npart_real, npart_ghost )

    !*****************************************************
    !
    !# Reallocates the APM fields member of TYPE [[apm_fields]]
    !
    !  FT 29.03.2022
    !
    !*****************************************************

    IMPLICIT NONE

    CLASS(apm_fields), INTENT(INOUT):: this
    !! The [[apm_fields]] object that this PROCEDURE is a member of
    INTEGER, INTENT(IN):: npart_real
    !! The number of real particles
    INTEGER, INTENT(IN):: npart_ghost
    !! The number of ghost particles

    CALL this% deallocate_apm_fields
    CALL this% allocate_apm_fields( npart_real, npart_ghost )

    this% all_pos_prev_dump= zero
    this% cnt_move= 0

  END SUBROUTINE reallocate_apm_fields


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

    USE utility,             ONLY: cnt, spherical_from_cartesian
    USE constants,           ONLY: half, third, Msun, amu, pi

    USE sph_variables,       ONLY: h, nu
    USE set_h,               ONLY: exact_nei_tree_update
    USE units,               ONLY: umass

    USE APM,                 ONLY: density_loop, position_correction, assign_h
    USE analyze,             ONLY: COM
    USE matrix,              ONLY: determinant_4x4_matrix

    USE sphincs_sph,         ONLY: density, ncand!, all_clists
    USE RCB_tree_3D,         ONLY: nic, nfinal, nprev, lpart, &
                                   rpart
    USE matrix,              ONLY: invert_3x3_matrix

    IMPLICIT NONE

    INTEGER,          PARAMETER:: max_npart        = 10D+6
    INTEGER,          PARAMETER:: nn_des           = 301
    INTEGER,          PARAMETER:: m_max_it         = 50
    INTEGER,          PARAMETER:: search_pos       = 10
    !INTEGER,          PARAMETER:: print_step       = 15
    DOUBLE PRECISION, PARAMETER:: eps              = 5.0D-1
    DOUBLE PRECISION, PARAMETER:: ellipse_thickness= 1.1D0
    !DOUBLE PRECISION, PARAMETER:: ghost_dist       = 0.25D0!0.25D0 !30.0D0
    DOUBLE PRECISION, PARAMETER:: multiple_h_av    = 1.0D0
    DOUBLE PRECISION, PARAMETER:: tol              = 1.0D-3
    !DOUBLE PRECISION, PARAMETER:: iter_tol         = 2.0D-2
    !DOUBLE PRECISION, PARAMETER:: backup_h         = 0.25D0

    INTEGER:: a, a2, itr, itr2, n_inc, cnt1!, inde, index1   ! iterators
    INTEGER:: npart_real, npart_real_half, npart_ghost, npart_all
    INTEGER:: nx, ny, nz, i, j, k
    INTEGER:: a_numin, a_numin2, a_numax, a_numax2
    INTEGER:: dim_seed, rel_sign
    INTEGER:: n_problematic_h, ill, l, itot

    DOUBLE PRECISION:: smaller_radius, larger_radius, radius_y, radius_z
    DOUBLE PRECISION:: h_max, h_av, tmp, dens_min, atmosphere_density!, delta
    DOUBLE PRECISION:: xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, &
                       rad_x, rad_y, rad_z, com_x, com_y, com_z, com_d
    DOUBLE PRECISION:: max_r_real, r_real, max_z_real
    DOUBLE PRECISION:: xtemp, ytemp, ztemp, x_ell, y_ell, z_ell
    DOUBLE PRECISION:: min_nu, max_nu, min_nu2, max_nu2
    ! The value of nu equal for all the particles, used during the APM iteration
    DOUBLE PRECISION:: nu_all
    DOUBLE PRECISION:: err_N_mean_min, err_N_mean_min_old, err_N_mean, &
                       err_mean_old, err_n_min, err_N_max, &!dNstar, &
                       nstar_p_err, nstar_real_err, dN_max, dN_av, &
                       r_tmp, theta_tmp, phi_tmp
    DOUBLE PRECISION:: art_pr_max
    DOUBLE PRECISION:: nu_tot, nu_ratio, nu_tmp2, nuratio_tmp
    DOUBLE PRECISION:: variance_nu, stddev_nu, mean_nu
    DOUBLE PRECISION:: variance_dN, stddev_dN
    DOUBLE PRECISION:: rand_num, rand_num2
    DOUBLE PRECISION:: r, theta, phi
    DOUBLE PRECISION:: r_ell, theta_ell, phi_ell

    INTEGER, DIMENSION(:), ALLOCATABLE:: neighbors_lists
    INTEGER, DIMENSION(:), ALLOCATABLE:: n_neighbors
    INTEGER, DIMENSION(:), ALLOCATABLE:: seed

    !DOUBLE PRECISION:: ha, ha_1, ha_3, va, mat(3,3), mat_1(3,3), xa, ya, za
    !DOUBLE PRECISION:: mat_xx, mat_xy, mat_xz, mat_yy
    !DOUBLE PRECISION:: mat_yz, mat_zz, Wdx, Wdy, Wdz, ddx, ddy, ddz, Wab, &
    !                   Wab_ha, Wi, Wi1, dvv

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos_tmp
    DOUBLE PRECISION, DIMENSION(3):: pos_corr_tmp
    DOUBLE PRECISION, DIMENSION(3):: pos_maxerr

    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pvol_tmp
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nu_one

    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar_int
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: particle_density_final

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: nearest_neighbors

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: ghost_pos
    !# Array storing the positions of the ghost particles. Only used when
    !  placing the ghost particles.
    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: ghost_pos_tmp
    !! Cache-array used when placing ghost particles

    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: h_tmp
    !! Cache-array used to temporarily store the smoothing lengths
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: h_guess_tmp
    !# Cache-array used to temporarily store the first guesses for the
    !  smoothing lengths

    LOGICAL:: exist

    !CHARACTER:: it_n
    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    LOGICAL, PARAMETER:: debug= .FALSE.
    LOGICAL:: few_ncand!, invertible_matrix

    TYPE(apm_fields):: fld

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
    DO a= 1, npart_real, 1
      h_guess(a)= three*(pvol(a)**third)
    ENDDO

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
      IF( pos_tmp( 3, a ) > zero )THEN
        itr= itr + 1
        pos_input( 1, itr )= pos_tmp( 1, a )
        pos_input( 2, itr )= pos_tmp( 2, a )
        pos_input( 3, itr )= pos_tmp( 3, a )
        h_guess( itr )     = h_tmp( itr )
      ENDIF
    ENDDO
    npart_real_half= itr

    DO a= 1, npart_real_half, 1
      pos_input( 1, npart_real_half + a )=   pos_input( 1, a )
      pos_input( 2, npart_real_half + a )=   pos_input( 2, a )
      pos_input( 3, npart_real_half + a )= - pos_input( 3, a )
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

    !smaller_radius= ABS( MINVAL( pos_input( 1, : ), DIM= 1 ) - center )
    !larger_radius = ABS( center - MAXVAL( pos_input( 1, : ), DIM= 1 ) )
    !radius_y= ABS( MAXVAL( pos_input( 2, : ), DIM= 1 ) )
    !radius_z= ABS( MAXVAL( pos_input( 3, : ), DIM= 1 ) )

  !  IF( pos_input( 1, 10 ) < 0 )THEN
  !
  !    smaller_radius= MIN( binary% get_radius1_x_comp(), &
  !                         binary% get_radius1_x_opp() )
  !    larger_radius = MAX( binary% get_radius1_x_comp(), &
  !                         binary% get_radius1_x_opp() )
  !    radius_y= binary% get_radius1_y()
  !    radius_z= binary% get_radius1_z()
  !
  !  ELSE
  !
  !    smaller_radius= MIN( binary% get_radius2_x_comp(), &
  !                         binary% get_radius2_x_opp() )
  !    larger_radius = MAX( binary% get_radius2_x_comp(), &
  !                         binary% get_radius2_x_opp() )
  !    radius_y= binary% get_radius2_y()
  !    radius_z= binary% get_radius2_z()
  !
  !  ENDIF
    smaller_radius= MIN( sizes(1), sizes(2) )
    larger_radius = MAX( sizes(1), sizes(2) )
    radius_y= sizes(3)
    radius_z= sizes(5)

    h_max= zero
    h_av = zero
    itr  = 0
    max_z_real= ABS( MAXVAL( pos_input( 3, : ), DIM= 1 ) )
    DO a= 1, npart_real, 1

      IF( SQRT( ( pos_input( 1, a ) - center(1) )**two &
              + ( pos_input( 2, a ) - center(2) )**two &
              + ( pos_input( 3, a ) - center(3) )**two ) &
                > 0.99D0*max_z_real )THEN

        itr= itr + 1
        IF( h_guess(a) > h_max )THEN
          h_max= h_guess(a)
        ENDIF
        h_av= h_av + h_guess(a)

      ENDIF

    ENDDO
    h_av= h_av/itr
    IF( debug ) PRINT *, "h_av=", h_av
    IF( debug ) PRINT *

    IF( debug ) PRINT *, "2"

    !-------------------------------!
    !--  Placing ghost particles  --!
    !-------------------------------!

    IF(.NOT.ALLOCATED( ghost_pos ))THEN
      ALLOCATE( ghost_pos( 3, max_npart ), STAT= ios, &
          ERRMSG= err_msg )
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

   ! max_r_real= zero
   ! DO itr= 1, npart_real, 1
   !
   !   r_real= SQRT( ( pos_input( 1, itr ) - center(1) )**two &
   !               + ( pos_input( 2, itr ) - center(2) )**two &
   !               + ( pos_input( 3, itr ) - center(3) )**two )
   !   IF( r_real > max_r_real ) max_r_real= r_real
   !
   ! ENDDO

    nx= nx_gh
    ny= ny_gh
    nz= nz_gh
    xmin= center(1) - sizes(1)*( one + eps )
    xmax= center(1) + sizes(2)*( one + eps )
    ymin= center(2) - sizes(3)*( one + eps )
    ymax= center(2) + sizes(4)*( one + eps )
    zmin= center(3) - sizes(5)*( one + eps )
    zmax= center(3) + sizes(6)*( one + eps )
    dx= ABS( xmax - xmin )/DBLE( nx )
    dy= ABS( ymax - ymin )/DBLE( ny )
    dz= ABS( zmax - zmin )/DBLE( nz )

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

    rad_x= larger_radius + ghost_dist !+ multiple_h_av*h_av
    rad_y= radius_y      + ghost_dist !+ multiple_h_av*h_av
    rad_z= radius_z      + ghost_dist !+ multiple_h_av*h_av

    PRINT *, "** Distance between the size of the object and the ghost ", &
             "particles: ghost_dist =", ghost_dist
    PRINT *

    IF( debug ) PRINT *, "larger_radius= ", larger_radius
    IF( debug ) PRINT *, "radius_y= ", radius_y
    IF( debug ) PRINT *, "radius_z= ", radius_z
    IF( debug ) PRINT *, "rad_x= ", rad_x
    IF( debug ) PRINT *, "rad_y= ", rad_y
    IF( debug ) PRINT *, "rad_z= ", rad_z
    IF( debug ) PRINT *

    ghost_pos_tmp= HUGE(zero)

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( nx, ny, nz, xmin, ymin, zmin, dx, dy, dz, &
    !$OMP                     ghost_pos_tmp, &
    !$OMP                     center, rad_x, rad_y, rad_z ) &
    !$OMP             PRIVATE( i, j, k, xtemp, ytemp, ztemp, &
    !$OMP                      x_ell, y_ell, z_ell, r, theta, phi, &
    !$OMP                      r_ell, theta_ell, phi_ell )
    DO k= 1, nz, 1

      ztemp= zmin + dz/two + DBLE( k - 1 )*dz

      DO j= 1, ny, 1

        ytemp= ymin + dy/two + DBLE( j - 1 )*dy

        DO i= 1, nx, 1

          xtemp= xmin + dx/two + DBLE( i - 1 )*dx

          CALL spherical_from_cartesian( xtemp, ytemp, ztemp, &
                                         center(1), center(2), center(3), &
                                         r, theta, phi )

          x_ell= center(1) + rad_x*SIN(theta)*COS(phi)
          y_ell= center(2) + rad_y*SIN(theta)*SIN(phi)
          z_ell= center(3) + rad_z*COS(theta)

       !  x_ell= center(1) + rad_x*COS(ATAN( ( ytemp - center(2) )/( xtemp - center(1) ) )) &
       !         *SIN(ACOS(( ztemp - center(3) )/SQRT( ( xtemp - center(1) )**two &
       !         + ( ytemp - center(2) )**two &
       !         + ( ztemp - center(3) )**two )))
       !
       !  y_ell= center(2) + rad_y*SIN(ATAN( ( ytemp - center(2) )/( xtemp - center(1) ) )) &
       !         *SIN(ACOS(( ztemp - center(3) )/SQRT( ( xtemp - center(1) )**two &
       !         + ( ytemp - center(2) )**two &
       !         + ( ztemp - center(3) )**two )))
       !
       !  z_ell= center(3) + rad_z*( ( ztemp - center(3) )/SQRT( ( xtemp - center(1) )**two &
       !         + ( ytemp - center(2) )**two &
       !         + ( ztemp - center(3) )**two ))

          CALL spherical_from_cartesian( x_ell, y_ell, z_ell, &
                                         center(1), center(2), center(3), &
                                         r_ell, theta_ell, phi_ell )

          IF( ( r <= ellipse_thickness*r_ell .AND. r >= r_ell &
                .AND. &
                get_density( xtemp, ytemp, ztemp ) <= zero ) &
           !   .OR. &
           !   ( SQRT( ( xtemp - center )**two + ytemp**two &
           !   + ztemp**two ) <= 55zero&
           !   !ellipse_thickness*SQRT( ( x_ell - center )**two &
           !   !            + y_ell**two + z_ell**two ) &
           !   .AND. &
           !   SQRT( ( xtemp - center )**two + ytemp**two &
           !         + ztemp**two ) >= 50zero ) & !
           !   !SQRT( ( x_ell - center )**two + y_ell**two &
           !   !      + z_ell**two ) )
          )THEN

            ghost_pos_tmp( 1, i, j, k )= xtemp
            ghost_pos_tmp( 2, i, j, k )= ytemp
            ghost_pos_tmp( 3, i, j, k )= ztemp

        !  ELSE
        !
        !    PRINT *, SQRT( ( xtemp - center )**two + ytemp**two &
        !    + ztemp**two ) <= & !
        !    ellipse_thickness*SQRT( ( x_ell - center )**two &
        !                + y_ell**two + z_ell**two )
        !    PRINT *, SQRT( ( xtemp - center )**two + ytemp**two &
        !    + ztemp**two ) >= &
        !    SQRT( ( x_ell - center )**two + y_ell**two &
        !    + z_ell**two )
        !    PRINT *, get_density( xtemp, ytemp, ztemp ) <= zero

          ENDIF

         ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    itr= 0
    DO k= 1, nz, 1

      DO j= 1, ny, 1

        DO i= 1, nx, 1

          IF( ghost_pos_tmp( 1, i, j, k ) < HUGE(zero) )THEN

            itr= itr + 1
            ghost_pos( 1, itr )= ghost_pos_tmp( 1, i, j, k )
            ghost_pos( 2, itr )= ghost_pos_tmp( 2, i, j, k )
            ghost_pos( 3, itr )= ghost_pos_tmp( 3, i, j, k )

          ENDIF

         ENDDO
      ENDDO
    ENDDO
    npart_ghost= itr
    IF( npart_ghost == 0 )THEN
      PRINT *, "** ERROR: No ghost particles were placed. Is the PARAMETER ", &
               "'ghost_dist' appropriate for the physical system?"
      PRINT *, "Stopping.."
      PRINT *
      STOP
    ENDIF
    ghost_pos = ghost_pos( :, 1:npart_ghost )

    PRINT *, " * ", npart_ghost, " ghost particles placed around ", &
             npart_real, "real particles."
    PRINT *

    DEALLOCATE( ghost_pos_tmp )

    PRINT *, " * Printing ghost particles to file..."

    IF( PRESENT(namefile_pos_id) )THEN
      finalnamefile= namefile_pos_id
    ELSE
      finalnamefile= "apm_pos_id.dat"
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
        pos_input( 1, a ), &
        pos_input( 2, a ), &
        pos_input( 3, a )
    ENDDO

    DO a= 1, npart_ghost, 1
      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        2, a, &
        ghost_pos( 1, a ), &
        ghost_pos( 2, a ), &
        ghost_pos( 3, a )
    ENDDO

    CLOSE( UNIT= 2 )

    PRINT *, " * Positions of ghost and real particles printed to ", &
             TRIM(finalnamefile), " ."

    !STOP

    npart_all= npart_real + npart_ghost

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

    all_pos( :, 1:npart_real )          = pos_input
    all_pos( :, npart_real+1:npart_all )= ghost_pos

    h_guess= h_guess(1:npart_all)
    h_guess(npart_real+1:npart_all)= ( dx*dy*dz )**third

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
    check_h1: DO a= 1, npart_all, 1
    ! find_h_backup, called below, is OMP parallelized, so this loop
    ! should not be parallelized as well

      IF( .NOT.is_finite_number( h(a) ) .OR. h(a) <= zero )THEN

        n_problematic_h= n_problematic_h + 1
        h(a)= find_h_backup( a, npart_all, all_pos, nn_des )
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

    IF(.NOT.ALLOCATED( nstar_real ))THEN
      ALLOCATE( nstar_real( npart_all ), STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array nstar_real in SUBROUTINE ", &
                  "perform_apm. The error message is",&
                  err_msg
         STOP
      ENDIF
    ENDIF
    !CALL allocate_apm_fields( npart_real, npart_ghost )

    nu= one
    CALL density_loop( npart_all, all_pos, &    ! input
                       nu, h, nstar_real )      ! output

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( all_pos, npart_real, nstar_real, h, nu, &
    !$OMP                     center ) &
    !$OMP             PRIVATE( a )
    DO a= 1, npart_real, 1

      IF( .NOT.is_finite_number( nstar_real( a ) ) )THEN

        PRINT *, "** ERROR! nstar_real(", a, ") is not a finite number on ", &
                 "a real particle!"
        IF( debug ) PRINT *, " * h(", a, ")=", h(a)
        IF( debug ) PRINT *, " * nu(", a, ")=", nu(a)
        IF( debug ) PRINT *, " * all_pos(", a, ")=", all_pos(:,a)
        IF( debug ) PRINT *, " * r(", a, ")=", &
                              SQRT( ( all_pos(1,a) - center )**two &
                              + all_pos(2,a)**two + all_pos(3,a)**two )
        PRINT *, " * Check if the smoothing length is 0 for some particles,", &
                 "   and if so, make its initial guess, h_guess, a bit larger."
        PRINT *
        STOP

      ENDIF

    ENDDO
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( all_pos, npart_real, nstar_real, h, nu, &
    !$OMP                     center, npart_all ) &
    !$OMP             PRIVATE( a )
    DO a= npart_real + 1, npart_all, 1

      IF( .NOT.is_finite_number( nstar_real( a ) ) )THEN

        PRINT *, "** ERROR! nstar_real(", a, ") is not a finite number on ", &
                 "a ghost particle!"
        IF( debug ) PRINT *, " * h(", a, ")=", h(a)
        IF( debug ) PRINT *, " * nu(", a, ")=", nu(a)
        IF( debug ) PRINT *, " * all_pos(", a, ")=", all_pos(:,a)
        IF( debug ) PRINT *, " * r(", a, ")=", &
                              SQRT( ( all_pos(1,a) - center )**two &
                              + all_pos(2,a)**two + all_pos(3,a)**two )
        PRINT *, " * Check if the smoothing length is 0 for some particles,", &
                 "   and if so, make its initial guess, h_guess, a bit larger."
        PRINT *
        STOP

      ENDIF

    ENDDO
    !$OMP END PARALLEL DO

    IF( debug ) PRINT *, "4"

    IF(.NOT.ALLOCATED( nstar_p ))THEN
      ALLOCATE( nstar_p( npart_all ), STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array nstar_p in SUBROUTINE ", &
                  "perform_apm. The error message is",&
                  err_msg
         STOP
      ENDIF
    ENDIF

    IF( debug ) PRINT *, "5"

    IF(.NOT.ALLOCATED( lapse           ))ALLOCATE( lapse          (npart_real) )
    IF(.NOT.ALLOCATED( shift_x         ))ALLOCATE( shift_x        (npart_real) )
    IF(.NOT.ALLOCATED( shift_y         ))ALLOCATE( shift_y        (npart_real) )
    IF(.NOT.ALLOCATED( shift_z         ))ALLOCATE( shift_z        (npart_real) )
    IF(.NOT.ALLOCATED( g_xx            ))ALLOCATE( g_xx           (npart_real) )
    IF(.NOT.ALLOCATED( g_xy            ))ALLOCATE( g_xy           (npart_real) )
    IF(.NOT.ALLOCATED( g_xz            ))ALLOCATE( g_xz           (npart_real) )
    IF(.NOT.ALLOCATED( g_yy            ))ALLOCATE( g_yy           (npart_real) )
    IF(.NOT.ALLOCATED( g_yz            ))ALLOCATE( g_yz           (npart_real) )
    IF(.NOT.ALLOCATED( g_zz            ))ALLOCATE( g_zz           (npart_real) )
    IF(.NOT.ALLOCATED( baryon_density  ))ALLOCATE( baryon_density (npart_real) )
    IF(.NOT.ALLOCATED( energy_density  ))ALLOCATE( energy_density (npart_real) )
    IF(.NOT.ALLOCATED( specific_energy ))ALLOCATE( specific_energy(npart_real) )
    IF(.NOT.ALLOCATED( pressure        ))ALLOCATE( pressure       (npart_real) )
    IF(.NOT.ALLOCATED( v_euler_x       ))ALLOCATE( v_euler_x      (npart_real) )
    IF(.NOT.ALLOCATED( v_euler_y       ))ALLOCATE( v_euler_y      (npart_real) )
    IF(.NOT.ALLOCATED( v_euler_z       ))ALLOCATE( v_euler_z      (npart_real) )

    IF( debug ) PRINT *, "6"

    IF( debug ) PRINT *, "7"

    max_nu= zero
    min_nu= HUGE(one)!1.0D60

    IF(.NOT.ALLOCATED( art_pr ))THEN
      ALLOCATE( art_pr( npart_all ), STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array art_pr in SUBROUTINE ", &
                  "perform_apm. The error message is",&
                  err_msg
         STOP
      ENDIF
    ENDIF

    IF( debug ) PRINT *, "7"

    CALL get_nstar_p_atm( npart_real, all_pos(1,1:npart_real), &
                                      all_pos(2,1:npart_real), &
                                      all_pos(3,1:npart_real), &
                                      nstar_p, &
                                      use_atmosphere )

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

    CALL impose_equatorial_plane_symmetry_apm( npart_all, npart_real, &
                                               npart_ghost, &
                                               all_pos, nu= nu )

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

    !IF( .NOT.ALLOCATED(freeze) ) ALLOCATE( freeze( npart_all ) )
    IF( .NOT.ALLOCATED(correction_pos) ) ALLOCATE( correction_pos( 3, npart_all ) )
    !IF( .NOT.ALLOCATED(all_pos_tmp) ) ALLOCATE( all_pos_tmp( 3, npart_all ) )
    !IF( .NOT.ALLOCATED(all_pos_prev) ) ALLOCATE( all_pos_prev( 3, npart_all ) )
    IF( .NOT.ALLOCATED(all_pos_prev_dump) ) ALLOCATE( all_pos_prev_dump( 3, npart_all ) )
    IF( .NOT.ALLOCATED(cnt_move) ) ALLOCATE( cnt_move( npart_real ) )
    cnt_move= 0

    ! Set the particles to be equal-mass
    nu_all= (mass/DBLE(npart_real))*umass/amu
    nu= nu_all
!    !$OMP PARALLEL DO DEFAULT( NONE ) &
!    !$OMP             SHARED( npart_real, freeze ) &
!    !$OMP             PRIVATE( a )
!    DO a= 1, npart_real, 1
!      freeze(a)= 0
!    ENDDO
!    !$OMP END PARALLEL DO
!    !$OMP PARALLEL DO DEFAULT( NONE ) &
!    !$OMP             SHARED( npart_real, npart_all, freeze ) &
!    !$OMP             PRIVATE( a )
!    DO a= npart_real + 1, npart_all ,1
!      freeze(a)= 1
!    ENDDO
!    !$OMP END PARALLEL DO

    CALL correct_center_of_mass( npart_real, all_pos(:,1:npart_real), &
                                 nu(1:npart_real), get_density, &
                                 validate_position_final, com_star, &
                                 verbose= .TRUE. )

    IF(.NOT.ALLOCATED( dNstar ))THEN
      ALLOCATE( dNstar( npart_real ), STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array dNstar in SUBROUTINE ", &
                  "perform_apm. The error message is",&
                  err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( nu_tmp ))THEN
      ALLOCATE( nu_tmp( npart_real ), STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array nu_tmp in SUBROUTINE ", &
                  "perform_apm. The error message is",&
                  err_msg
         STOP
      ENDIF
    ENDIF
    !IF(.NOT.ALLOCATED( dN ))THEN
    !  ALLOCATE( dN( npart_real ), STAT= ios, ERRMSG= err_msg )
    !  IF( ios > 0 )THEN
    !     PRINT *, "...allocation error for array dN in SUBROUTINE ", &
    !              "perform_apm. The error message is",&
    !              err_msg
    !     STOP
    !  ENDIF
    !ENDIF
    IF(.NOT.ALLOCATED( good_rho ))THEN
      ALLOCATE( good_rho( npart_real ), STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array good_rho in SUBROUTINE ", &
                  "perform_apm. The error message is",&
                  err_msg
         STOP
      ENDIF
    ENDIF

    !all_pos_prev     = -one
    all_pos_prev_dump= zero
    PRINT *, " * The APM iteration starts here."
    PRINT *

    n_inc= 0
    err_N_mean_min= HUGE(one)
    apm_iteration: DO itr= 1, apm_max_it, 1

      PRINT *, "---------------------------------------------------------------"
      PRINT *, " * Starting with APM step #: ", itr
      PRINT *

      dump_pos: IF( (print_step /= 0) &
                    .AND. &
                    (MOD( itr, print_step ) == 0) )THEN

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

        IF( debug ) PRINT *, "printing positions to file..."

        IF( PRESENT(namefile_pos) )THEN
          finalnamefile= namefile_pos
        ELSE
          finalnamefile= "apm_pos.dat"
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
          tmp= get_density( all_pos( 1, a ), all_pos( 2, a ), all_pos( 3, a ) )
          WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
            1, a, &
            all_pos( 1, a ), &
            all_pos( 2, a ), &
            all_pos( 3, a ), &
            tmp, cnt_move(a), &
            all_pos_prev_dump( 1, a ), &
            all_pos_prev_dump( 2, a ), &
            all_pos_prev_dump( 3, a )
        ENDDO

        DO a= npart_real + 1, npart_all, 1
          tmp= get_density( all_pos( 1, a ), all_pos( 2, a ), all_pos( 3, a ) )
          WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
            2, a, &
            all_pos( 1, a ), &
            all_pos( 2, a ), &
            all_pos( 3, a ), &
            tmp
        ENDDO

        CLOSE( UNIT= 2 )

        all_pos_prev_dump= all_pos

      ENDIF dump_pos

      IF( debug ) PRINT *, "enforcing center of mass..."

      CALL correct_center_of_mass( npart_real, all_pos(:,1:npart_real), &
                                   nu(1:npart_real), get_density, &
                                   validate_position_final, com_star )

      IF( debug ) PRINT *, "mirroring particles..."

      !CALL impose_equatorial_plane_symmetry_apm( npart_real, npart_ghost, &
      !                                all_pos(:,1:npart_real), &
      !                                nu= nu(1:npart_real), &
      !                                pos_prev= all_pos_prev(:,1:npart_real) )

      CALL impose_equatorial_plane_symmetry_apm( npart_all, npart_real, &
                                                 npart_ghost, all_pos, nu= nu )

      IF( debug ) PRINT *, "assign h..."

      h_guess(1:npart_real)= h(1:npart_real)

      CALL assign_h( nn_des, &
                     npart_all, &
                     all_pos, h_guess, h )

      find_h_bruteforce_timer= timer( "find_h_bruteforce_timer" )
      CALL find_h_bruteforce_timer% start_timer()
      n_problematic_h= 0
      check_h2: DO a= 1, npart_all, 1
      ! find_h_backup, called below, is OMP parallelized, so this loop
      ! should not be parallelized as well

        IF( .NOT.is_finite_number( h(a) ) .OR. h(a) <= zero )THEN

          n_problematic_h= n_problematic_h + 1
          h(a)= find_h_backup( a, npart_all, all_pos, nn_des )
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
                         nu, h, nstar_real )      ! output

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( all_pos, npart_real, nstar_real, h, nu, &
      !$OMP                     center ) &
      !$OMP             PRIVATE( a )
      DO a= 1, npart_real, 1

        IF( .NOT.is_finite_number( nstar_real( a ) ) )THEN

          PRINT *, "** ERROR! nstar_real(", a, ") is not a finite number on ", &
                   "a real particle!"
          IF( debug ) PRINT *, " * h(", a, ")=", h(a)
          IF( debug ) PRINT *, " * nu(", a, ")=", nu(a)
          IF( debug ) PRINT *, " * all_pos(", a, ")=", all_pos(:,a)
          IF( debug ) PRINT *, " * r(", a, ")=", &
                                SQRT( ( all_pos(1,a) - center )**two &
                                + all_pos(2,a)**two + all_pos(3,a)**two )
          PRINT *, " * Check if the smoothing length is 0 for some particles,", &
                   "   and if so, make its initial guess, h_guess, a bit larger."
          PRINT *
          STOP

        ENDIF

      ENDDO
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( all_pos, npart_real, nstar_real, h, nu, &
      !$OMP                     center, npart_all ) &
      !$OMP             PRIVATE( a )
      DO a= npart_real + 1, npart_all, 1

        IF( .NOT.is_finite_number( nstar_real( a ) ) )THEN

          PRINT *, "** ERROR! nstar_real(", a, ") is not a finite number on ", &
                   "a ghost particle!"
          IF( debug ) PRINT *, " * h(", a, ")=", h(a)
          IF( debug ) PRINT *, " * nu(", a, ")=", nu(a)
          IF( debug ) PRINT *, " * all_pos(", a, ")=", all_pos(:,a)
          IF( debug ) PRINT *, " * r(", a, ")=", &
                                SQRT( ( all_pos(1,a) - center )**two &
                                + all_pos(2,a)**two + all_pos(3,a)**two )
          PRINT *, " * Check if the smoothing length is 0 for some particles,", &
                   "   and if so, make its initial guess, h_guess, a bit larger."
          PRINT *
          STOP

        ENDIF

      ENDDO
      !$OMP END PARALLEL DO

      CALL get_nstar_p_atm( npart_real, all_pos(1,1:npart_real), &
                                        all_pos(2,1:npart_real), &
                                        all_pos(3,1:npart_real), &
                                        nstar_p, &
                                        use_atmosphere )

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( all_pos, npart_all, nstar_p, h, nu, &
      !$OMP                     center, dNstar ) &
      !$OMP             PRIVATE( a )
      check_nstar_p: DO a= 1, npart_all, 1

        IF( .NOT.is_finite_number( nstar_p( a ) ) )THEN

          PRINT *, "** ERROR! nstar_p(", a, ") is a not a finite number ", &
                   "in SUBROUTINE perform_apm!"
          PRINT *, "   nstar_p(", a, ")=", nstar_p(a)
          PRINT *, "   dNstar(", a, ")=", dNstar(a)
          PRINT *, "   rho(", a, ")=", get_density( all_pos(1,a), &
                                                    all_pos(2,a), &
                                                    all_pos(3,a) )
          IF( debug ) PRINT *, " * h(", a, ")=", h(a)
          IF( debug ) PRINT *, " * nu(", a, ")=", nu(a)
          IF( debug ) PRINT *, " * all_pos(", a, ")=", all_pos(:,a)
          IF( debug ) PRINT *, " * r(", a, ")=", &
                                SQRT( ( all_pos(1,a) - center )**two &
                                + all_pos(2,a)**two + all_pos(3,a)**two )
          PRINT *
          STOP

        ENDIF

      ENDDO check_nstar_p
      !$OMP END PARALLEL DO

      art_pr_max= zero
      err_N_max=  zero
      err_N_min=  HUGE(one)
      err_N_mean= zero

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( npart_real, nstar_real, nstar_p, &
      !$OMP                     dNstar, art_pr ) &
      !$OMP             PRIVATE( a )
      assign_artificial_pressure_on_real_particles: DO a= 1, npart_real, 1

        IF( nstar_p(a) <= zero )THEN

          dNstar(a)= zero

        ELSE

          dNstar(a)= ( nstar_real(a) - nstar_p(a) )/nstar_p(a)

        ENDIF
        art_pr(a) = MAX( one + dNstar(a), one/ten )

      ENDDO assign_artificial_pressure_on_real_particles
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( npart_all, npart_real, art_pr, nstar_real, &
      !$OMP                     nstar_p, dNstar, all_pos ) &
      !$OMP             PRIVATE( a )
      find_nan_in_art_pr: DO a= 1, npart_real, 1

        IF( .NOT.is_finite_number(art_pr(a)) )THEN
          PRINT *, "** ERROR! art_pr(", a, ")= ", art_pr(a), &
                   " is not a finite number on a real particle!"
          PRINT *, "   nstar_real(", a, ")=", nstar_real(a)
          PRINT *, "   nstar_p(", a, ")=", nstar_p(a)
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

      art_pr_max= MAXVAL( art_pr, DIM= 1 )

      IF( .NOT.is_finite_number(art_pr_max) )THEN
        PRINT *, "** ERROR! art_pr_max =", art_pr_max
        PRINT *, " * Stopping..."
        PRINT *
        STOP
      ENDIF

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( nstar_p, npart_real, good_rho ) &
      !$OMP             PRIVATE( a )
      DO a= 1, npart_real, 1
        good_rho(a)= nstar_p(a) > zero
      ENDDO
      !$OMP END PARALLEL DO

      err_N_max     = MAXVAL( ABS(dNstar), MASK= good_rho )
      err_N_min     = MINVAL( ABS(dNstar), MASK= good_rho )
      err_N_mean    = SUM( ABS(dNstar), DIM= 1 )/SIZE( dNstar )
      err_N_mean_min= MIN( err_N_mean, err_N_mean_min )

      !IF( err_N_mean_min == err_N_mean )THEN
      !  all_pos_best= all_pos
      !ENDIF

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( all_pos, npart_real, nstar_real, nstar_p, &
      !$OMP                     dNstar, err_N_max, &
      !$OMP                     pos_maxerr, nstar_real_err, nstar_p_err ) &
      !$OMP             PRIVATE( a )
      DO a= 1, npart_real, 1

        !IF( ABS(dNstar(a)) > err_N_max &
        !    .AND. &
        !    get_density( all_pos(1,a), &
        !                 all_pos(2,a), &
        !                 all_pos(3,a) ) > zero )THEN
        !
        !  err_N_max     = ABS(dNstar(a))
        !  pos_maxerr    = all_pos(:,a)
        !  nstar_real_err= nstar_real(a)
        !  nstar_p_err   = nstar_p(a)
        !
        !ENDIF

        IF( dNstar(a) == err_N_max )THEN
          pos_maxerr    = all_pos(:,a)
          nstar_real_err= nstar_real(a)
          nstar_p_err   = nstar_p(a)
        ENDIF

        !err_N_max = MAX( err_N_max, ABS(dNstar) )
        !IF( get_density( all_pos(1,a), &
        !                 all_pos(2,a), &
        !                 all_pos(3,a) ) > zero )THEN
        !
        !  err_N_min = MIN( err_N_min, ABS(dNstar(a)) )
        !  err_N_mean= err_N_mean + ABS(dNstar(a))
        !
        !ENDIF

        IF( ISNAN(dNstar(a)) )THEN
          PRINT *, "dNstar= ", dNstar(a), " at particle ", a
          PRINT *, "nstar_real= ", nstar_real(a)
          PRINT *, "nstar_p= ", nstar_p(a)
          STOP
        ENDIF

      ENDDO
      !$OMP END PARALLEL DO

      IF( .NOT.is_finite_number( art_pr_max ) )THEN
        PRINT *, "** ERROR! art_pr_max is not a finite number!", &
                 " Stopping.."
        PRINT *
        STOP
      ENDIF
      IF( .NOT.is_finite_number( nu_all ) )THEN
        PRINT *, "** ERROR! nu_all is a not a finite number!", &
                 " Stopping.."
        PRINT *
        STOP
      ENDIF

      !
      !-- Assign artificial pressure to the ghost particles
      !

      nstar_p( npart_real+1:npart_all )= zero
      !art_pr ( npart_real+1:npart_all )= 6.0D0*art_pr_max

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( all_pos, npart_all, npart_real, center, &
      !$OMP                     dNstar, art_pr, rad_x, rad_y, rad_z, &
      !$OMP                     art_pr_max ) &
      !$OMP             PRIVATE( a, r, theta, phi, x_ell, y_ell, z_ell, r_ell, &
      !$OMP                      itr, itr2 )
      assign_artificial_pressure_on_ghost_particles: &
      DO a= npart_real + 1, npart_all, 1

        CALL spherical_from_cartesian( &
                              all_pos(1,a), all_pos(2,a), all_pos(3,a), &
                              center(1), center(2), center(3), &
                              r, theta, phi )

        x_ell= center(1) + rad_x*SIN(theta)*COS(phi)
        y_ell= center(2) + rad_y*SIN(theta)*SIN(phi)
        z_ell= center(3) + rad_z*COS(theta)

        r_ell= SQRT( ( x_ell - center(1) )**two &
                   + ( y_ell - center(2) )**two &
                   + ( z_ell - center(3) )**two )

        shell_loop: DO itr2= 0, 20, 1

          IF( r <= ( one + ( ellipse_thickness - one )*DBLE(itr2)/ten )*r_ell &
              .AND. &
              r >= ( one + ( ellipse_thickness - one )*DBLE(itr2-1)/ten )*r_ell&
          ! If the ghost particle is contained within the i-th spherical shell..

              !r <= ( 500.D0 + 50.D0*DBLE(itr2)/ten ) &
              !.AND. &
              !r > ( 500.D0 + 50.D0*DBLE(itr2-1)/ten ) &

          )THEN

            ! ..assign a pressure that increases with i, to build a pressure
            !   gradient
            art_pr(a)= ten*three*three*art_pr_max
            !art_pr(a)= DBLE(itr+1)*art_pr_max

          !ELSE
          !
          !  art_pr( a )= art_pr_max

          ENDIF

        ENDDO shell_loop

      ENDDO assign_artificial_pressure_on_ghost_particles
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( npart_all, npart_real, art_pr, nstar_real, &
      !$OMP                     nstar_p, all_pos ) &
      !$OMP             PRIVATE( a )
      find_nan_in_art_pr_ghost: DO a= npart_real + 1, npart_all, 1

        IF( .NOT.is_finite_number(art_pr(a)) )THEN
          PRINT *, "** ERROR! art_pr(", a, ")= ", art_pr(a), &
                   " is not a finite number on a ghost particle!"
          PRINT *, "   nstar_real(", a, ")=", nstar_real(a)
          PRINT *, "   nstar_p(", a, ")=", nstar_p(a)
          PRINT *, "   rho(", a, ")=", get_density( all_pos(1,a), &
                                                    all_pos(2,a), &
                                                    all_pos(3,a) )
          PRINT *, " * Stopping.."
          PRINT *
          STOP
        ENDIF

      ENDDO find_nan_in_art_pr_ghost
      !$OMP END PARALLEL DO

      IF( debug ) PRINT *, "Before calling position_correction"

      IF( debug ) PRINT *, npart_all
      IF( debug ) PRINT *, SIZE(all_pos(1,:))
      IF( debug ) PRINT *, SIZE(h)
      IF( debug ) PRINT *, SIZE(art_pr)
      IF( debug ) PRINT *, SIZE(nstar_real)
      IF( debug ) PRINT *, SIZE(correction_pos(1,:))

      PRINT *, " * Maximum relative error between the star density profile", &
               " and its SPH estimate: err_N_max= ", err_N_max
      PRINT *, "     at position: x=", pos_maxerr(1), ", y=", pos_maxerr(2), &
               ", z=", pos_maxerr(3)
      PRINT *, "     with r/(system size)= ", SQRT( &
                              ( ABS(pos_maxerr(1)) - ABS(center(1)) )**two &
                            + ( ABS(pos_maxerr(2)) - ABS(center(2)) )**two &
                            + ( ABS(pos_maxerr(3)) - ABS(center(3)) )**two ) &
                            /sizes(1)
      PRINT *, "   The ID density is   = ", nstar_p_err
      PRINT *, "   The SPH estimate is= ", nstar_real_err
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
                         nu, h, nstar_real )      ! output
      nu= nu_all

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( nu_tmp, nu, nstar_p, nstar_real, &
      !$OMP                     nuratio_thres, npart_real ) &
      !$OMP             PRIVATE( nu_tmp2, a )
      cap_nu: DO a= 1, npart_real, 1

        nu_tmp2= nu(a)
        nu_tmp(a)= nstar_p(a)/nstar_real(a)

          IF( nu_tmp(a) > nu_tmp2*SQRT(nuratio_thres) ) nu_tmp(a)= &
                                            nu_tmp2*SQRT(nuratio_thres)
          IF( nu_tmp(a) < nu_tmp2/SQRT(nuratio_thres) ) nu_tmp(a)= &
                                            nu_tmp2/SQRT(nuratio_thres)

      ENDDO cap_nu
      !$OMP END PARALLEL DO
      nuratio_tmp= MAXVAL( nu_tmp(1:npart_real), DIM= 1 )&
                  /MINVAL( nu_tmp(1:npart_real), DIM= 1 )

      PRINT *, " * If the APM iteration was stopped at this step, ", &
               "being the threshold for the baryon number ratio equal to "
      PRINT *, "   nuratio_thres=", nuratio_thres, ","
      PRINT *, "   the baryon number ratio would be equal to the following:"
      PRINT *, " * Temporary CORRECTED maximum baryon number at this step=", &
               MAXVAL( nu_tmp(1:npart_real), DIM= 1 )
      PRINT *, " * Temporary CORRECTED minimum baryon number at this step=", &
               MINVAL( nu_tmp(1:npart_real), DIM= 1 )
      PRINT *, " * Temporary CORRECTED baryon number ratio at this step=", &
               nuratio_tmp

      IF( nuratio_des /= zero )THEN

        PRINT *, " * The APM iteration will stop when the baryon number ", &
                 "ratio is between ", nuratio_des*(one - quarter/ten), &
                 " and ", nuratio_des*(one + quarter/ten)
        PRINT *

      ENDIF

      IF( err_N_mean > err_mean_old )THEN
        n_inc= n_inc + 1
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

      PRINT *, " * npart_real= ",  npart_real
      PRINT *, " * npart_ghost= ", npart_ghost
      PRINT *, " * npart_all= ",   npart_all
      PRINT *

      !
      !-- EXIT conditions
      !
      IF( nuratio_des > zero )THEN

        IF( ( nuratio_tmp >= nuratio_des*(one - quarter/ten) .AND. &
              nuratio_tmp <= nuratio_des*(one + quarter/ten) .AND. &
              nuratio_tmp /= nuratio_thres ) .OR. itr == apm_max_it ) EXIT

      ELSE

        PRINT *, " * n_inc= ", n_inc
        PRINT *
        IF( n_inc == max_inc .OR. itr == apm_max_it ) EXIT

      ENDIF
      err_mean_old      = err_N_mean
      err_N_mean_min_old= err_N_mean_min

      !
      !-- If the particle distribution is not yet good enough, update it
      !
      PRINT *, " * Updating positions..."

      !all_pos_prev= all_pos

      IF( debug ) PRINT *, "21"

      CALL position_correction( npart_all, &
                                all_pos, h, nu_all, art_pr, nstar_real, &
                                correction_pos )

      IF( debug ) PRINT *, "22"

      correction_pos(:,npart_real+1:npart_all)= zero

      IF( debug ) PRINT *, "23"

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( npart_all, correction_pos, all_pos, center, &
      !$OMP                     larger_radius ) &
      !$OMP             PRIVATE( a, itr2, r, theta, phi )
      find_nan_in_correction_pos: DO a= 1, npart_all, 1

        loop_over_spatial_components: DO itr2= 1, 3, 1

          IF( .NOT.is_finite_number( correction_pos( itr2, a ) ) )THEN

            CALL spherical_from_cartesian( all_pos(1,a), all_pos(2,a), &
                                           all_pos(3,a), &
                                           center(1), center(2), center(3), &
                                           r, theta, phi )

            correction_pos( 1, a )= - one/(two*ten)*SIN(theta)*COS(phi)
            correction_pos( 2, a )= - one/(two*ten)*SIN(theta)*SIN(phi)
            correction_pos( 3, a )= - one/(two*ten)*COS(theta)
            !correction_pos( itr2, a )= zero

           ! PRINT *, "** ERROR! correction_pos(", itr2, ",", a, ") is a NaN!"
           ! PRINT *, " *        correction_pos: x=", correction_pos(1,a), &
           !          ", y=", correction_pos(2,a), ", z=", correction_pos(3,a)
           ! PRINT *, " *        Particle position: x=", all_pos(1,a), &
           !          ", y=", all_pos(2,a), ", z=", all_pos(3,a)
           !
           ! CALL spherical_from_cartesian( &
           !                   all_pos(1,a), all_pos(2,a), all_pos(3,a), &
           !                   center(1), center(2), center(3), &
           !                   r_tmp, theta_tmp, phi_tmp )
           !
           ! !r_tmp= SQRT( ( all_pos(1,a) - center(1) )**two + &
           ! !             ( all_pos(2,a) - center(2) )**two + &
           ! !             ( all_pos(3,a) - center(3) )**two )
           !
           ! PRINT *, " *        Particle radius: r=", r_tmp, &
           !          "=", r_tmp/larger_radius*ten*ten, &
           !          "% of the larger radius of the star."
           ! PRINT *, " *        Particle colatitude: theta=", theta_tmp/pi," pi"
           ! PRINT *, " *        Particle longitude: phi=", phi_tmp/pi, " pi"
           ! PRINT *, " * Stopping.."
           ! PRINT *
           ! STOP

          ENDIF

        ENDDO loop_over_spatial_components

      ENDDO find_nan_in_correction_pos
      !$OMP END PARALLEL DO

      IF( debug ) PRINT *, "After calling position_correction"

      cnt_move= 0
      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( use_atmosphere, all_pos, correction_pos, &
      !$OMP                     dNstar, npart_real, nstar_p, cnt_move, center )&
      !$OMP             PRIVATE( pos_corr_tmp, a, cnt, rand_num, rand_num2, &
      !$OMP                      rel_sign, r, theta, phi, theta_tmp )
      particle_loop: DO a= 1, npart_real, 1

        adapt_displacement_to_error: &
        IF( dNstar(a) >= ten*ten &
            .AND. &
            validate_position_final( &
              all_pos(1,a) + ten*correction_pos(1,a), &
              all_pos(2,a) + ten*correction_pos(2,a), &
              all_pos(3,a) + ten*correction_pos(3,a) ) == 1 )THEN

          pos_corr_tmp= all_pos(:,a) + ten*correction_pos(:,a) ! 10


        ELSEIF( dNstar(a) >= ten &
                .AND. &
                validate_position_final( &
                  all_pos(1,a) + three*correction_pos(1,a), &
                  all_pos(2,a) + three*correction_pos(2,a), &
                  all_pos(3,a) + three*correction_pos(3,a) ) == 1 )THEN

          pos_corr_tmp= all_pos(:,a) + three*correction_pos(:,a) ! 3


        ELSE

          pos_corr_tmp= all_pos(:,a) + correction_pos(:,a) ! 1

        ENDIF adapt_displacement_to_error

        if_atmosphere: IF( use_atmosphere )THEN
        ! If the atmosphere is used...

          all_pos(:,a)= pos_corr_tmp
          cnt_move(a)= 1
          !...move the particle without any validation, and exit the loop

        ELSE

          cnt= 0
          determine_new_position: DO

            test_position: IF( get_density( &
                pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) > zero &
                .AND. &
                validate_position_final( &
                    pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) == 0 &
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
            !    and the current step is lower than search_pos

              cnt= cnt + 1

              CALL RANDOM_NUMBER( rand_num )
              CALL RANDOM_NUMBER( rand_num2 )

              IF( rand_num2 < half )  rel_sign= - 1
              IF( rand_num2 >= half ) rel_sign=   1

              pos_corr_tmp(1)= all_pos(1,a) + &
                correction_pos(1,a)*( one + DBLE(rel_sign)*rand_num*half )

              CALL RANDOM_NUMBER( rand_num )
              CALL RANDOM_NUMBER( rand_num2 )

              IF( rand_num2 < half )  rel_sign= - 1
              IF( rand_num2 >= half ) rel_sign=   1

              pos_corr_tmp(2)= all_pos(2,a) + &
                correction_pos(2,a)*( one + DBLE(rel_sign)*rand_num*half )

              CALL RANDOM_NUMBER( rand_num )
              CALL RANDOM_NUMBER( rand_num2 )

              IF( rand_num2 < half )  rel_sign= - 1
              IF( rand_num2 >= half ) rel_sign=   1

              pos_corr_tmp(3)= all_pos(3,a) + &
                correction_pos(3,a)*( one + DBLE(rel_sign)*rand_num*half )

              !pos_corr_tmp*( one + DBLE(rel_sign)*rand_num*half*third )

              ! ...change the new position randomly, independently in x, y, z,
              !    and repeat the test

            ELSEIF( cnt == search_pos + 1 )THEN
            ! ...else if the new position was changed randomly search_pos
            !    times, do not move the particle at this step,
            !    and exit the 'determine_new_position' loop

              ! cnt= cnt + 1
              ! CALL RANDOM_NUMBER( rand_num )
              ! CALL RANDOM_NUMBER( rand_num2 )
              !
              ! IF( rand_num2 < half )  rel_sign= - 1
              ! IF( rand_num2 >= half ) rel_sign=   1
              ! all_pos(:,a)= all_pos(:,a)*( one -rand_num*half*third )

              CALL spherical_from_cartesian( all_pos(1,a), all_pos(2,a), &
                                             all_pos(3,a), &
                                             center(1), center(2), center(3), &
                                             r, theta, phi )

              !all_pos(1,a)= all_pos(1,a) &
              !      - one/(two*ten)*SIN(theta*(one + one/(ten*ten)))*COS(phi)
              !all_pos(2,a)= all_pos(2,a) &
              !      - one/(two*ten)*SIN(theta*(one + one/(ten*ten)))*SIN(phi)
              !all_pos(3,a)= all_pos(3,a) &
              !      - one/(two*ten)*COS(theta*(one + one/(ten*ten)))

              ! ...push the particle closer to the origin close to its radial
              !    line, and place it to the first place having an acceptable
              !    density
              ! TODO: what if the particles is placed on top of another
              !       particle?
              scan_radial_line: DO itr2= (ten - three)*ten, 1, -1

               ! CALL RANDOM_NUMBER( rand_num )
               ! CALL RANDOM_NUMBER( rand_num2 )
               !
               ! IF( rand_num2 < half )  rel_sign= - 1
               ! IF( rand_num2 >= half ) rel_sign=   1
               !
               ! theta_tmp= theta*( one + DBLE(rel_sign)*rand_num*half )
               !
               ! IF( theta_tmp < pi .AND. theta_tmp > pi/two ) theta= theta_tmp
               !
               ! CALL RANDOM_NUMBER( rand_num )
               ! CALL RANDOM_NUMBER( rand_num2 )
               !
               ! IF( rand_num2 < half )  rel_sign= - 1
               ! IF( rand_num2 >= half ) rel_sign=   1
               !
               ! phi= phi*( one + DBLE(rel_sign)*rand_num*half )

                pos_corr_tmp(1)= center(1) +r*itr2/(ten*ten)*SIN(theta)*COS(phi)
                pos_corr_tmp(2)= center(2) +r*itr2/(ten*ten)*SIN(theta)*SIN(phi)
                pos_corr_tmp(3)= center(3) +r*itr2/(ten*ten)*COS(theta)

                IF( get_density( &
                    pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) > zero &
                    .AND. &
                    validate_position_final( &
                    pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) == 0 &
                )THEN

                  all_pos(:,a)= pos_corr_tmp
                  cnt_move(a)= 1
                  EXIT

                ENDIF

              ENDDO scan_radial_line

              EXIT

            ENDIF test_position

          ENDDO determine_new_position

        ENDIF if_atmosphere

      ENDDO particle_loop
      !$OMP END PARALLEL DO

      PRINT *, " * The fraction of particles that moved at this step is", &
               DBLE(SUM(cnt_move))/DBLE(npart_real)
      PRINT *

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( npart_all, all_pos ) &
      !$OMP             PRIVATE( a, itr2 )
      find_nan_in_all_pos: DO a= 1, npart_all, 1

        DO itr2= 1, 3, 1
          IF( ISNAN( all_pos( itr2, a ) ) )THEN
            PRINT *, "** ERROR! all_pos(", itr2, ",", a, ") is a NaN!", &
                     " Stopping.."
            PRINT *
            STOP
          ENDIF
        ENDDO

      ENDDO find_nan_in_all_pos
      !$OMP END PARALLEL DO

      ! If some of the particles crossed the xy plane top-down in the
      ! last step, reflect them back above the xy plane

   !   !$OMP PARALLEL DO DEFAULT( NONE ) &
   !   !$OMP             SHARED( all_pos, all_pos_prev, npart_real ) &
   !   !$OMP             PRIVATE( a )
   !   DO a= 1, npart_real, 1
   !
   !     IF( all_pos_prev( 3, a ) > 0 .AND. &
   !         all_pos( 3, a ) <= 0 )THEN
   !
   !       all_pos( 3, a )= all_pos_prev( 3, a )
   !
   !     ENDIF
   !
   !   ENDDO
   !   !$OMP END PARALLEL DO

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

    IF(ALLOCATED( pos )) DEALLOCATE( pos )
    ALLOCATE( pos( 3, npart_real ), STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...allocation error for array pos in SUBROUTINE ", &
                "perform_apm. The error message is", &
                err_msg
       STOP
    ENDIF

    pos= all_pos( :, 1:npart_real )
    IF( debug ) PRINT *, npart

    h      = h(1:npart_real)
    h_guess= h_guess(1:npart_real)
    nu     = nu(1:npart_real)

    !------------------------------------------------!
    !-- Discard atmosphere, if present and desired --!
    !------------------------------------------------!

    IF( use_atmosphere .AND. remove_atmosphere )THEN

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
      !$OMP             SHARED( pos, pos_tmp, h, h_tmp, npart_real, &
      !$OMP                     nstar_p, h_guess, h_guess_tmp, nu, nu_tmp, &
      !$OMP                     pvol_tmp, pvol ) &
      !$OMP             PRIVATE( a ) &
      !$OMP             REDUCTION( +: npart )
      DO a= 1, npart_real, 1

        IF( nstar_p(a) > zero )THEN
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
   !   !$OMP             SHARED( pos, pos_tmp, h, h_tmp, npart_real, &
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

    ENDIF

    !---------------!
    !-- Set npart --!
    !---------------!

    npart= npart_real

    !----------------------------!
    !-- enforce centre of mass --!
    !----------------------------!

    PRINT *, " * Correcting center of mass..."
    PRINT *

    CALL correct_center_of_mass( npart_real, pos, nu, get_density, &
                                 validate_position_final, com_star, &
                                 verbose= .TRUE. )

    !-----------------------------------------------------------------------!
    !-- Mirror the positions after having repositioned the center of mass --!
    !-----------------------------------------------------------------------!

    PRINT *, " * Imposing equatorial plane symmetry..."
    PRINT *

    !CALL impose_equatorial_plane_symmetry_apm( npart_all, npart_real, &
    !                                           npart_ghost, all_pos )

    CALL impose_equatorial_plane_symmetry( npart_real, all_pos, nu= nu )

    !-----------------------------!
    !-- Print positions to file --!
    !-----------------------------!

    IF( PRESENT(namefile_pos) )THEN
      finalnamefile= namefile_pos
    ELSE
      finalnamefile= "apm_pos.dat"
    ENDIF

    PRINT *, " * Printing APM positions to formatted file ", &
             TRIM(finalnamefile) ,"..."
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

    DO a= 1, npart_real, 1
      tmp= get_density( pos( 1, a ), pos( 2, a ), pos( 3, a ) )
      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        1, a, &
        pos( 1, a ), &
        pos( 2, a ), &
        pos( 3, a ), &
        tmp, cnt_move(a)
    ENDDO

    DO a= npart_real + 1, npart_all, 1
      tmp= get_density( all_pos( 1, a ), all_pos( 2, a ), all_pos( 3, a ) )
      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        2, a, &
        all_pos( 1, a ), &
        all_pos( 2, a ), &
        all_pos( 3, a ), &
        tmp
    ENDDO

    CLOSE( UNIT= 2 )

    !-------------------------------------------------------------------!
    !-- now assign baryon number to match profile as good as possible --!
    !-------------------------------------------------------------------!

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

   ! find_problem_in_h_2: DO a= 1, npart_real, 1
   !
   !   IF( .NOT.is_finite_number( h( a ) ) )THEN
   !     PRINT *, "** ERROR! h(", a, ") is a NaN!"
   !     PRINT *, " * h_guess(", a, ")= ", h_guess(a)
   !     PRINT *, " * all_pos(:,", a, ")= ", all_pos(:,a)
   !     PRINT *, " Stopping..."
   !     PRINT *
   !     STOP
   !   ENDIF
   !   IF( h( a ) <= zero )THEN
   !    ! PRINT *, "** ERROR! h(", a, ") is zero or negative!"
   !    ! PRINT *, " * h_guess(", a, ")= ", h_guess(a)
   !    ! PRINT *, " * all_pos(:,", a, ")= ", all_pos(:,a)
   !    ! PRINT *, " * h(", a, ")= ", h(a)
   !    ! PRINT *, " Stopping..."
   !    ! PRINT *
   !    ! STOP
   !     IF( a == 1 )THEN
   !       DO a2= 2, npart_real, 1
   !         IF( h( a2 ) > zero )THEN
   !           h( a )= h( a2 )
   !           EXIT
   !         ENDIF
   !       ENDDO
   !     ELSE
   !       h(a) = h(a - 1)
   !     ENDIF
   !   ENDIF
   !
   ! ENDDO find_problem_in_h_2

    IF( debug ) PRINT *, "2"

    ! Measure SPH particle number density
    nu= one
    CALL density_loop( npart_real, pos, &    ! input
                       nu, h, nstar_real )      ! output                   

    IF( debug ) PRINT *, "3"

    CALL get_nstar_p_atm( npart_real, pos(1,:), pos(2,:), pos(3,:), &
                          nstar_p, use_atmosphere )

    nu= nu_all
    PRINT *, " * Baryon number on all particles before correction nu_all= ", &
             nu_all

    !nu_tot= zero
    !DO a= 1, npart_real, 1
    !  nu_tot= nu_tot + nu(a)
    !ENDDO
    nu_tot= SUM( nu, DIM= 1 )

    PRINT *, " * Total baryon number nu_tot=", nu_tot
    PRINT *, " * Total baryon mass= ", nu_tot*amu/MSun, "=", &
             ten*ten*nu_tot*amu/MSun/mass, "% of the LORENE baryon mass"
    PRINT *

    IF( debug ) PRINT *, "4"

    IF( debug ) PRINT *, "npart_real= ", npart_real
    IF( debug ) PRINT *, "SIZE(nu)= ", SIZE(nu)
    IF( debug ) PRINT *

    IF( debug ) nu_ratio= MAXVAL( nu, DIM= 1 )/MINVAL( nu, DIM= 1 )
    IF( debug ) PRINT *, " * nu_ratio before correction = ", nu_ratio
    IF( debug ) PRINT *

    !----------------------!
    !-- Correcting nu... --!
    !----------------------!

    PRINT *, " * Correcting nu..."
    PRINT *

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( nu_tmp, nu, nstar_p, nstar_real, &
    !$OMP                     nuratio_thres, npart_real ) &
    !$OMP             PRIVATE( nu_tmp2, a )
    DO a= 1, npart_real, 1

      nu_tmp2= nu(a)
      nu(a)= nstar_p(a)/nstar_real(a)

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
    !$OMP             SHARED( npart_real, nu, nstar_real, nstar_p ) &
    !$OMP             PRIVATE( a )
    DO a= 1, npart_real, 1

      IF( ISNAN( nu(a) ) )THEN
        PRINT *, " * ERROR! nu(", a, ") is a NaN."
        PRINT *, " nstar_real(a)= ", nstar_real(a)
        PRINT *, " nstar_p(a)= ", nstar_p(a)
        PRINT *, " Stopping..."
        PRINT *
        STOP
      ENDIF
      IF( nu(a) < zero )THEN
        PRINT *, " * ERROR! nu(", a, ") is negative."
        PRINT *, " nu(a)= ", nu(a)
        PRINT *, " nstar_real(a)= ", nstar_real(a)
        PRINT *, " nstar_p(a)= ", nstar_p(a)
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
!    !$OMP PARALLEL DO DEFAULT( NONE ) &
!    !$OMP             SHARED( nu, max_nu, npart_real, &
!    !$OMP                     nuratio_thres, npart_real ) &
!    !$OMP             PRIVATE( a )
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
!    !$OMP END PARALLEL DO

    PRINT *, " * Baryon number assigned."
    PRINT *

    IF( mass_it )THEN

      ! just a few iterations to NOT get the nu-ratio too large
      mass_iteration: DO itr= 1, m_max_it, 1

         ! measure density
         CALL density_loop( npart_real, pos, &    ! input
                            nu, h, nstar_real )      ! output


         CALL get_nstar_p_atm( npart_real, pos(1,:), &
                                           pos(2,:), &
                                           pos(3,:), nstar_p, use_atmosphere )

         !nstar_p( npart_real+1:npart_all )= zero

         ! get RELATIVE nu's right
         dN_av= zero
         max_nu= zero
         min_nu= HUGE(one)
         DO a= 1, npart_real, 1

          dNstar(a)= (nstar_real(a)-nstar_p(a))/nstar_p(a)
          nu(a)= nu(a)*(one - dNstar(a))
          dN_av= dN_av + dNstar(a)
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

  !  PRINT *, "max_nu=", max_nu
  !  PRINT *, "        at ", pos(:, a_numax), " r= ", &
  !           NORM2( pos(:, a_numax) )/larger_radius
  !  PRINT *, "min_nu=", min_nu
  !  PRINT *, "        at ", pos(:, a_numin), " r= ", &
  !           NORM2( pos(:, a_numin) )/larger_radius
  !  PRINT *, "max_nu/min_nu=", max_nu/min_nu
  !  PRINT *
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

    PRINT *, "Excluding the absolute max and min of nu:"
    PRINT *
    PRINT *, "max_nu=", max_nu2
    PRINT *, "        at ", pos(:, a_numax2), " r= ", &
             NORM2( pos(:, a_numax2) )/larger_radius
    PRINT *, "min_nu=", min_nu2
    PRINT *, "        at ", pos(:, a_numin2), " r= ", &
             NORM2( pos(:, a_numin2) )/larger_radius
    PRINT *, "max_nu/min_nu=", max_nu2/min_nu2
    PRINT *

    !nu_tot= zero
    !DO a= 1, npart_real, 1
    !  nu_tot= nu_tot + nu(a)
    !ENDDO
    nu_tot= SUM(nu, DIM= 1)
    mean_nu= nu_tot/npart_real

    variance_nu = zero                       ! compute variance
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( npart_real, nu, mean_nu ) &
    !$OMP             PRIVATE( a ) &
    !$OMP             REDUCTION( +: variance_nu )
    DO a = 1, npart_real, 1
      variance_nu = variance_nu + (nu(a) - mean_nu)**two
    ENDDO
    !$OMP END PARALLEL DO
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
             ten*ten*nu_tot*amu/MSun/mass, "% of the LORENE baryon mass"
    PRINT *

    IF( correct_nu )THEN

      nu= nu/(nu_tot*amu/Msun/mass)
      !nu_tot= zero
      !DO a= 1, npart_real, 1
      !  nu_tot= nu_tot + nu(a)
      !ENDDO
      nu_tot= SUM(nu, DIM= 1)

      PRINT *, "After correcting nu to match the mass of the star..."
      PRINT *
      PRINT *, "nu_tot=", nu_tot
      PRINT *, "mass estimate= ", nu_tot*amu/MSun, "=", &
               ten*ten*nu_tot*amu/MSun/mass, "% of the LORENE baryon mass"
      PRINT *

    ENDIF


  !  finalnamefile= "dbg_apm_pos1.dat"
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
  !  tmp= get_density( pos( 1, a ), pos( 2, a ), pos( 3, a ) )
  !  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
  !    a, &
  !    pos( 1, a ), &
  !    pos( 2, a ), &
  !    pos( 3, a ), &
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

    !CALL impose_equatorial_plane_symmetry_apm( npart_all, npart_real, &
    !                                           npart_ghost, all_pos )

    CALL impose_equatorial_plane_symmetry( npart_real, all_pos, nu= nu )


  !  finalnamefile= "dbg_apm_pos2.dat"
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
  !    tmp= get_density( pos( 1, a ), pos( 2, a ), pos( 3, a ) )
  !    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
  !      a, &
  !      pos( 1, a ), &
  !      pos( 2, a ), &
  !      pos( 3, a ), &
  !      tmp
  !  ENDDO
  !
  !  CLOSE( UNIT= 2 )

    !-------------------!
    !-- monitoring... --!
    !-------------------!

    ! measure density
    CALL density_loop( npart_real, pos, &    ! input
                       nu, h, nstar_real )      ! output

    CALL get_nstar_p_atm( npart_real, pos(1,:), &
                                      pos(2,:), &
                                      pos(3,:), nstar_p, use_atmosphere )

    dNstar= zero
    dN_av = zero
    dN_max= zero
    cnt1= 0
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( npart_real, pos, nstar_real, nstar_p, dNstar ) &
    !$OMP             PRIVATE( a ) &
    !$OMP             REDUCTION( +: dN_av, cnt1 )
    DO a= 1, npart_real, 1

      !IF( get_density( pos(1,a), pos(2,a), pos(3,a) ) > zero )THEN
      IF( nstar_p(a) > zero )THEN

        dNstar(a)= ABS(nstar_real(a)-nstar_p(a))/nstar_p(a)
        dN_av=  dN_av + dNstar(a)
        cnt1= cnt1 + 1

      ENDIF

    ENDDO
    !$OMP END PARALLEL DO
    dN_av= dN_av/DBLE(cnt1)
    dN_max= MAXVAL(dNstar, DIM= 1)

    variance_dN = zero                       ! compute variance
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( npart_real,pos,nstar_real,nstar_p,dNstar,dN_av )&
    !$OMP             PRIVATE( a ) &
    !$OMP             REDUCTION( +: variance_dN, cnt1 )
    DO a= 1, npart_real, 1

      IF( nstar_p(a) > zero )THEN

        variance_dN = variance_dN + (dNstar(a) - dN_av)**two
        cnt1= cnt1 + 1

      ENDIF

    ENDDO
    !$OMP END PARALLEL DO
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
    DO

      few_ncand= .FALSE.

      ! Redo the previous step slightly different (it's built-in;
      ! exact_nei_tree_update does not work if I don't call assign_h first),
      ! then update the neighbour-tree and fill the neighbour-data
      CALL exact_nei_tree_update( nn_des, &
                                  npart_real, &
                                  pos, nu )

      !ll_cell_loop: DO ill= 1, nfinal
      !
      !  itot= nprev + ill
      !  IF( nic(itot) == 0 ) CYCLE
      !
      !  IF( ncand(ill) < nn_des - 1 )THEN
      !
      !    ! Increase the smoothing lengths of the paricles inside the cell,
      !    ! and rebuild the tree
      !    few_ncand= .TRUE.
      !
      !    particle_in_cell_loop: DO l= lpart(itot), rpart(itot)
      !
      !      h(l)= three*h(l)
      !
      !    ENDDO particle_in_cell_loop
      !
      !  ELSE
      !
      !    few_ncand= .FALSE.
      !
      !  ENDIF
      !
      !ENDDO ll_cell_loop

      cnt1= cnt1 + 1

      IF( debug ) PRINT *, cnt1

      IF( .NOT.few_ncand .OR. cnt1 >= max_it_tree )THEN
        EXIT
      ENDIF

    ENDDO

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

    !
    !-- Check that the smoothing length is acceptable
    !
   ! check_h: DO a= 1, npart_real, 1
   !
   !   IF( .NOT.is_finite_number( h(a) ) )THEN
   !     PRINT *, "** ERROR! h(", a, ") is a NaN"
   !     !PRINT *, "Stopping..."
   !    ! PRINT *
   !     !STOP
   !     IF( a > npart_real/2 )THEN
   !       DO itr= CEILING(DBLE(npart_real/2)) - 1, 1, -1
   !         IF( h(itr) >= backup_h )THEN
   !           h(a) = h(itr)
   !           EXIT
   !         ENDIF
   !       ENDDO
   !     ELSE
   !       !h(a) = h(a - 1)
   !       DO itr= a + 1, npart_real, 1
   !         IF( h(itr) >= backup_h )THEN
   !           h(a) = h(itr)
   !           EXIT
   !         ENDIF
   !       ENDDO
   !     ENDIF
   !     !PRINT *, "** ERROR! h(", a, ")=", h(a)
   !     !PRINT *
   !   ENDIF
   !   IF( h(a) <= zero )THEN
   !     PRINT *, "** ERROR! h(", a, ")=", h(a)
   !     !PRINT *, "Stopping..."
   !     !PRINT *
   !     !STOP
   !     IF( a > npart_real/2 )THEN
   !       DO itr= CEILING(DBLE(npart_real/2)) - 1, 1, -1
   !         IF( h(itr) >= backup_h )THEN
   !           h(a) = h(itr)
   !           EXIT
   !         ENDIF
   !       ENDDO
   !     ELSE
   !       !h(a) = h(a - 1)
   !       DO itr= a + 1, npart_real, 1
   !         IF( h(itr) >= backup_h )THEN
   !           h(a) = h(itr)
   !           EXIT
   !         ENDIF
   !       ENDDO
   !     ENDIF
   !     !PRINT *, "** ERROR! h(", a, ")=", h(a)
   !     !PRINT *
   !   ENDIF
   !
   ! ENDDO check_h

!    PRINT *
!    PRINT *, "nfinal= ", nfinal
!    ll_cell_loop: DO ill= 1, nfinal
!
!      itot= nprev + ill
!      IF( nic(itot) == 0 ) CYCLE
!
!      particle_loop: DO l= lpart(itot), rpart(itot)
!
!        a=         iorig(l)
!
!        ha=        h(a)
!        ha_1=      one/ha
!        ha_3=      ha_1*ha_1*ha_1
!
!        xa=        pos(1,a)
!        ya=        pos(2,a)
!        za=        pos(3,a)
!
!        ! initialize correction matrix
!        mat_xx=    0.D0
!        mat_xy=    0.D0
!        mat_xz=    0.D0
!        mat_yy=    0.D0
!        mat_yz=    0.D0
!        mat_zz=    0.D0
!
!        cnt1= 0
!        cnt2= 0
!
!        !finalnamefile= "mat.dat"
!        !
!        !INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )
!        !
!        !IF( exist )THEN
!        !  OPEN( UNIT= 20, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
!        !        FORM= "FORMATTED", &
!        !        POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
!        !        IOMSG= err_msg )
!        !ELSE
!        !  OPEN( UNIT= 20, FILE= TRIM(finalnamefile), STATUS= "NEW", &
!        !  FORM= "FORMATTED", &
!        !        ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
!        !ENDIF
!        !IF( ios > 0 )THEN
!        !  PRINT *, "...error when opening ", TRIM(finalnamefile), &
!        !           ". The error message is", err_msg
!        !  STOP
!        !ENDIF
!        cand_loop: DO k= 1, ncand(ill)
!
!          b=      all_clists(ill)%list(k)
!
!          IF( b == a )THEN
!            cnt1= cnt1 + 1
!          ENDIF
!          IF( xa == pos(1,b) .AND. ya == pos(2,b) .AND. za == pos(3,b) &
!          )THEN
!            cnt2= cnt2 + 1
!          ENDIF
!
!          ! Distances (ATTENTION: flatspace version !!!)
!          ddx=     xa - pos(1,b)
!          ddy=     ya - pos(2,b)
!          ddz=     za - pos(3,b)
!          va=     SQRT(ddx*ddx + ddy*ddy + ddz*ddz)*ha_1
!
!          !IF( dx == 0 .AND. dy == 0 .AND. dz == 0 )THEN
!          !  PRINT *, "va=", va
!          !  PRINT *, "dz=", dx
!          !  PRINT *, "dy=", dy
!          !  PRINT *, "dz=", dz
!          !  PRINT *, "xa=", xa
!          !  PRINT *, "ya=", ya
!          !  PRINT *, "za=", za
!          !  PRINT *, "pos_u(1,b)", pos_u(1,b)
!          !  PRINT *, "pos_u(2,b)", pos_u(2,b)
!          !  PRINT *, "pos_u(3,b)", pos_u(3,b)
!          !  STOP
!          !ENDIF
!
!          ! get interpolation indices
!          inde=  MIN(INT(va*dv_table_1),n_tab_entry)
!          index1= MIN(inde + 1,n_tab_entry)
!
!          ! get tabulated values
!          Wi=     W_no_norm(inde)
!          Wi1=    W_no_norm(index1)
!
!          ! interpolate
!          dvv=    (va - DBLE(inde)*dv_table)*dv_table_1
!          Wab_ha= Wi + (Wi1 - Wi)*dvv
!
!          ! "correction matrix" for derivs
!          Wdx=    Wab_ha*ddx
!          Wdy=    Wab_ha*ddy
!          Wdz=    Wab_ha*ddz
!          mat_xx= mat_xx + Wdx*ddx
!          mat_xy= mat_xy + Wdx*ddy
!          mat_xz= mat_xz + Wdx*ddz
!          mat_yy= mat_yy + Wdy*ddy
!          mat_yz= mat_yz + Wdy*ddz
!          mat_zz= mat_zz + Wdz*ddz
!
!          !WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
!          !  mat_xx, mat_xy, mat_xz, mat_yy, mat_yz, mat_zz
!
!        ENDDO cand_loop
!
!        !CLOSE( UNIT= 20 )
!
!        ! correction matrix
!        mat(1,1)= mat_xx
!        mat(2,1)= mat_xy
!        mat(3,1)= mat_xz
!
!        mat(1,2)= mat_xy
!        mat(2,2)= mat_yy
!        mat(3,2)= mat_yz
!
!        mat(1,3)= mat_xz
!        mat(2,3)= mat_yz
!        mat(3,3)= mat_zz
!
!        ! invert it
!        CALL invert_3x3_matrix( mat, mat_1, invertible_matrix )
!        !PRINT *, invertible_matrix
!
!        IF( .NOT.invertible_matrix )THEN
!          PRINT *, "a= ", a
!          PRINT *, "h(a)= ", h(a)
!          PRINT *, "pos_u= ", pos(1,b), pos(2,b), pos(3,b)
!          PRINT *, "nprev= ", nprev
!          PRINT *, "ill= ", ill
!          PRINT *, "itot= ", itot
!          PRINT *, "nfinal= ", nfinal
!          PRINT *, "ncand(ill)= ", ncand(ill)
!          PRINT *, "cnt1= ", cnt1
!          PRINT *, "cnt2= ", cnt2
!          PRINT *
!          STOP
!        ENDIF
!
!      ENDDO particle_loop
!
!    ENDDO ll_cell_loop

    PRINT *, " * Smoothing lengths assigned and tree is built."
    PRINT *, " * The smoothing length was found brute-force for ", &
             n_problematic_h, " particles."
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

    IF( PRESENT(namefile_results) )THEN
      finalnamefile= namefile_results
    ELSE
      finalnamefile= "apm_results.dat"
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

    IF( .NOT.ALLOCATED( nu_one ) ) ALLOCATE( nu_one( npart_real ) )
    IF( .NOT.ALLOCATED( particle_density_final ) ) &
      ALLOCATE( particle_density_final( npart_real ) )
    nu_one= one
    CALL density_loop( npart_real, pos, &    ! input
                       nu_one, h, particle_density_final )      ! output

    DO a= 1, npart_real, 1
      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        a, &
        pos( 1, a ), pos( 2, a ), pos( 3, a ), &
        nstar_p( a ), &
        nstar_int( a ), &
        particle_density_final( a ), &
        particle_density_final( a )*nstar_p( 1 )/particle_density_final( 1 ), &
        ABS(nstar_real(a)-nstar_p(a))/nstar_p(a), &
        nu(a), &
        nearest_neighbors(2,a)
    ENDDO

    CLOSE( UNIT= 2 )

    IF( debug ) PRINT *, "2"

    IF( ALLOCATED( pos_input ) ) DEALLOCATE( pos_input )
    ALLOCATE( pos_input( 3, npart_real ) )
    pos_input(:,1:npart_real)= pos(:,1:npart_real)

    IF( debug ) PRINT *, "2.5"

    IF( ALLOCATED( h_output ) ) DEALLOCATE( h_output )
    ALLOCATE( h_output( npart_real ) )
    h_output(1:npart_real)= h(1:npart_real)

    IF( debug ) PRINT *, "2.6"

    IF( ALLOCATED( nu_output ) ) DEALLOCATE( nu_output )
    ALLOCATE( nu_output( npart_real ) )
    nu_output(1:npart_real)= nu(1:npart_real)

    npart_output= npart_real

    IF( debug ) PRINT *, "3"

    IF( ALLOCATED(posmash) ) DEALLOCATE(posmash)
    CALL deallocate_metric_on_particles()
    IF( debug ) PRINT *, "4"
    CALL deallocate_gradient()
    IF( debug ) PRINT *, "5"
    CALL deallocate_RCB_tree_memory_3D()
    IF( debug ) PRINT *, "6"
    CALL deallocate_SPH_memory()

    !
    !-- Check that there aren't particles with the same positions
    !
    !IF( debug ) finalnamefile= "negative_hydro.dat"
    !IF( debug ) CALL THIS% analyze_hydro( finalnamefile )

    PRINT *, "** Checking that there aren't particles with the same position..."
    PRINT *

    CALL check_particle_positions( npart_real, pos )



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
      INTEGER:: answer
      !! validate_position( x, y, z ) if the latter is present, 0 otherwise

      IF( PRESENT(validate_position) )THEN

        answer= validate_position( x, y, z )
        !IF( validate_position( x, y, z ) == 1 ) answer= .TRUE.
        !IF( validate_position( x, y, z ) == 0 ) answer= .FALSE.

      ELSE

        answer= 0

      ENDIF

    END FUNCTION validate_position_final


    SUBROUTINE get_nstar_p_atm( npart_real, x, y, z, nstar_p, use_atmosphere )

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
      DOUBLE PRECISION, INTENT(OUT):: nstar_p(npart_real)
      !! Array to store the computed proper baryon number density
      LOGICAL,  INTENT( IN ):: use_atmosphere
      !# `.TRUE.` if an atmosphere should be used during the APM, to allow
      !  the real aprticles more freedom to move around and adjust;
      !  `.FALSE.` otherwise

      CALL get_nstar_p( npart_real, x, y, z, nstar_p )

      IF( use_atmosphere .EQV. .TRUE. )THEN

        !$OMP PARALLEL DO DEFAULT( NONE ) &
        !$OMP             SHARED( npart_real, nstar_p, atmosphere_density ) &
        !$OMP             PRIVATE( a )
        DO a= 1, npart_real, 1
          IF( nstar_p(a) <= atmosphere_density )THEN
            nstar_p(a)= atmosphere_density
          ENDIF
        ENDDO
        !$OMP END PARALLEL DO

      ELSE

        !$OMP PARALLEL DO DEFAULT( NONE ) &
        !$OMP             SHARED( npart_real, nstar_p ) &
        !$OMP             PRIVATE( a )
        DO a= 1, npart_real, 1

          IF( nstar_p( a ) < 1.0D-10 )THEN
            PRINT *, "** ERROR! nstar_p(", a, ")=", nstar_p( a ), &
                     " in SUBROUTINE get_nstar_p_atm."
            PRINT *, " * Stopping.."
            PRINT *
            STOP
          ENDIF

        ENDDO
        !$OMP END PARALLEL DO

      ENDIF

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( npart_real, nstar_p ) &
      !$OMP             PRIVATE( a )
      DO a= 1, npart_real, 1

        IF( .NOT.is_finite_number( nstar_p( a ) ) )THEN
          PRINT *, "** ERROR! nstar_p(", a, ") is a not a finite number!", &
                   " in SUBROUTINE get_nstar_p_atm."
          PRINT *, " * Stopping.."
          PRINT *
          STOP
        ENDIF

      ENDDO
      !$OMP END PARALLEL DO

    END SUBROUTINE get_nstar_p_atm


  END PROCEDURE perform_apm


  SUBROUTINE impose_equatorial_plane_symmetry_apm( npart_all, npart_real, &
                                                   npart_ghost, pos, nu, &
                                                   com_star, verbose )

    !*************************************************************
    !
    !# Mirror the particle with z>0 with respect to the xy plane,
    !  to impose the equatorial-plane symmetry
    !
    !  FT 1.09.2021
    !
    !*************************************************************

    USE constants,      ONLY: zero
    USE analyze,        ONLY: COM
    USE sph_variables,  ONLY: h

    IMPLICIT NONE

    INTEGER, INTENT(INOUT):: npart_all
    INTEGER, INTENT(INOUT):: npart_real
    INTEGER, INTENT(IN)   :: npart_ghost

    INTEGER:: npart_real_old

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT):: pos
    !TYPE(apm_fields), INTENT(INOUT):: fld

    DOUBLE PRECISION,                     INTENT(IN),    OPTIONAL:: com_star
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE,   INTENT(INOUT), OPTIONAL:: nu
    LOGICAL,                              INTENT(IN),    OPTIONAL:: verbose

    INTEGER:: a, npart_above_xy
    DOUBLE PRECISION:: com_x, com_y, com_z, com_d

    INTEGER, DIMENSION(:), ALLOCATABLE:: above_xy_plane_a

    DOUBLE PRECISION, DIMENSION(3,npart_all):: postmp
    DOUBLE PRECISION, DIMENSION(npart_all)  :: nutmp
    DOUBLE PRECISION, DIMENSION(npart_all)  :: htmp
    DOUBLE PRECISION, DIMENSION(npart_all)  :: nstarptmp
    DOUBLE PRECISION, DIMENSION(3,npart_all):: all_pos_prev_dump_tmp

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos_below
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nu_below

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile
    LOGICAL:: exist

    IF( MOD(npart_real,2) /= 0 )THEN
      PRINT *, "** ERROR! If the equatorial symmetry has to be imposed, ", &
               "the particle number must be even."
      PRINT *, " * npart_real= ", npart_real
      PRINT *, " * Stopping..."
      PRINT *
    ENDIF

    npart_real_old= npart_real
    postmp        = pos
    IF( PRESENT(nu) ) nutmp = nu
    htmp= h
    nstarptmp= nstar_p
    IF( ALLOCATED(all_pos_prev_dump) ) all_pos_prev_dump_tmp= all_pos_prev_dump

    CALL find_particles_above_xy_plane( npart_real, postmp(:,1:npart_real), &
                                        npart_above_xy, above_xy_plane_a )

    ALLOCATE( pos_below(3,npart_above_xy) )

    IF( PRESENT(nu) )THEN

      ALLOCATE( nu_below(npart_above_xy) )
      CALL reflect_particles_xy_plane( npart_real, postmp(:,1:npart_real), &
                                       pos_below, npart_above_xy, &
                                       above_xy_plane_a, nutmp(1:npart_real), &
                                       nu_below )

    ELSE

      CALL reflect_particles_xy_plane( npart_real, postmp(:,1:npart_real), &
                                       pos_below, npart_above_xy, &
                                       above_xy_plane_a )

    ENDIF

    IF( npart_real/2 /= npart_above_xy )THEN

      npart_real= 2*npart_above_xy
      npart_all= npart_real + npart_ghost

      DEALLOCATE(pos)
      ALLOCATE( pos(3,npart_all) )

      IF( PRESENT(nu) )THEN
        DEALLOCATE(nu)
        ALLOCATE( nu(npart_all) )
      ENDIF

      !CALL fld% reallocate_apm_fields( npart_real, npart_ghost )

      IF( ALLOCATED( nstar_real )        ) DEALLOCATE( nstar_real )
      IF( ALLOCATED( nstar_p )           ) DEALLOCATE( nstar_p )
      IF( ALLOCATED( art_pr )            ) DEALLOCATE( art_pr )
      IF( ALLOCATED( dNstar )            ) DEALLOCATE( dNstar )
      IF( ALLOCATED( nu_tmp )            ) DEALLOCATE( nu_tmp )
      !IF( ALLOCATED( dN )                ) DEALLOCATE( dN )
      IF( ALLOCATED( good_rho )          ) DEALLOCATE( good_rho )
      !IF( ALLOCATED( freeze )            ) DEALLOCATE( freeze )
      IF( ALLOCATED( correction_pos )    ) DEALLOCATE( correction_pos )
      !IF( ALLOCATED( all_pos_tmp )       ) DEALLOCATE( all_pos_tmp )
      !IF( ALLOCATED( all_pos_prev )      ) DEALLOCATE( all_pos_prev )
      IF( ALLOCATED( all_pos_prev_dump ) ) DEALLOCATE( all_pos_prev_dump )
      IF( ALLOCATED( cnt_move )          ) DEALLOCATE( cnt_move )

      ALLOCATE( nstar_real(npart_all) )
      ALLOCATE( nstar_p(npart_all) )
      ALLOCATE( art_pr(npart_all) )
      !ALLOCATE( freeze( npart_all ) )
      ALLOCATE( correction_pos( 3, npart_all ) )
      !ALLOCATE( all_pos_tmp( 3, npart_all ) )
      !ALLOCATE( all_pos_prev( 3, npart_all ) )
      ALLOCATE( all_pos_prev_dump( 3, npart_all ) )

      ALLOCATE( cnt_move( npart_real ) )
      ALLOCATE( dNstar( npart_real ) )
      ALLOCATE( nu_tmp( npart_real ) )
      !ALLOCATE( dN( npart_real ) )
      ALLOCATE( good_rho( npart_real ) )
      !all_pos_prev     = -one
      all_pos_prev_dump= zero
      cnt_move= 0
 !     !$OMP PARALLEL DO DEFAULT( NONE ) &
 !     !$OMP             SHARED( npart_real, freeze ) &
 !     !$OMP             PRIVATE( a )
 !     DO a= 1, npart_real, 1
 !       freeze(a)= 0
 !     ENDDO
 !     !$OMP END PARALLEL DO
 !     !$OMP PARALLEL DO DEFAULT( NONE ) &
 !     !$OMP             SHARED( npart_real, npart_all, freeze ) &
 !     !$OMP             PRIVATE( a )
 !     DO a= npart_real + 1, npart_all , 1
 !       freeze(a)= 1
 !     ENDDO
 !     !$OMP END PARALLEL DO

      DEALLOCATE( lapse           )
      DEALLOCATE( shift_x         )
      DEALLOCATE( shift_y         )
      DEALLOCATE( shift_z         )
      DEALLOCATE( g_xx            )
      DEALLOCATE( g_xy            )
      DEALLOCATE( g_xz            )
      DEALLOCATE( g_yy            )
      DEALLOCATE( g_yz            )
      DEALLOCATE( g_zz            )
      DEALLOCATE( baryon_density  )
      DEALLOCATE( energy_density  )
      DEALLOCATE( specific_energy )
      DEALLOCATE( pressure        )
      DEALLOCATE( v_euler_x       )
      DEALLOCATE( v_euler_y       )
      DEALLOCATE( v_euler_z       )
      ALLOCATE( lapse          (npart_real) )
      ALLOCATE( shift_x        (npart_real) )
      ALLOCATE( shift_y        (npart_real) )
      ALLOCATE( shift_z        (npart_real) )
      ALLOCATE( g_xx           (npart_real) )
      ALLOCATE( g_xy           (npart_real) )
      ALLOCATE( g_xz           (npart_real) )
      ALLOCATE( g_yy           (npart_real) )
      ALLOCATE( g_yz           (npart_real) )
      ALLOCATE( g_zz           (npart_real) )
      ALLOCATE( baryon_density (npart_real) )
      ALLOCATE( energy_density (npart_real) )
      ALLOCATE( specific_energy(npart_real) )
      ALLOCATE( pressure       (npart_real) )
      ALLOCATE( v_euler_x      (npart_real) )
      ALLOCATE( v_euler_y      (npart_real) )
      ALLOCATE( v_euler_z      (npart_real) )

      npart= npart_all
      IF( ALLOCATED(posmash) ) DEALLOCATE(posmash)
      CALL deallocate_metric_on_particles
      CALL deallocate_gradient
      CALL deallocate_RCB_tree_memory_3D
      CALL deallocate_SPH_memory

      CALL allocate_SPH_memory
      CALL allocate_RCB_tree_memory_3D(npart_all)
      iorig(1:npart_all)= (/ (a,a=1,npart_all) /)
      CALL allocate_gradient(npart_all)
      CALL allocate_metric_on_particles(npart_all)

    ENDIF

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( pos, npart_above_xy, nu, postmp, nutmp, &
    !$OMP                     pos_below, nu_below, h, htmp, above_xy_plane_a, &
    !$OMP                     all_pos_prev_dump, all_pos_prev_dump_tmp, &
    !$OMP                     nstar_p, nstarptmp ) &
    !$OMP             PRIVATE( a )
    DO a= 1, npart_above_xy, 1
      pos( :, a )= postmp( :, a )
      pos( :, npart_above_xy + a )= pos_below( :, a )
      h( a )= htmp( above_xy_plane_a(a) )
      h( npart_above_xy + a )= htmp( above_xy_plane_a(a) )
      nstar_p( a )= nstarptmp( above_xy_plane_a(a) )
      nstar_p( npart_above_xy + a )= nstarptmp( above_xy_plane_a(a) )
      IF( ALLOCATED(all_pos_prev_dump) )THEN
        all_pos_prev_dump(1,a)= all_pos_prev_dump_tmp(1,above_xy_plane_a(a))
        all_pos_prev_dump(2,a)= all_pos_prev_dump_tmp(2,above_xy_plane_a(a))
        all_pos_prev_dump(3,a)= all_pos_prev_dump_tmp(3,above_xy_plane_a(a))
        all_pos_prev_dump(1,npart_above_xy + a)= &
                                all_pos_prev_dump_tmp(1,above_xy_plane_a(a))
        all_pos_prev_dump(2,npart_above_xy + a)= &
                                all_pos_prev_dump_tmp(2,above_xy_plane_a(a))
        all_pos_prev_dump(3,npart_above_xy + a)= &
                              - all_pos_prev_dump_tmp(3,above_xy_plane_a(a))
      ENDIF
      IF( PRESENT(nu) )THEN
        nu( a )= nutmp( a )
        nu( npart_above_xy + a )= nu_below( a )
      ENDIF
    ENDDO
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( pos, npart_real, npart_all, nu, postmp, nutmp, &
    !$OMP                     npart_real_old, npart_ghost, h, htmp, &
    !$OMP                     all_pos_prev_dump, all_pos_prev_dump_tmp, &
    !$OMP                     nstar_p, nstarptmp ) &
    !$OMP             PRIVATE( a )
    DO a= 1, npart_ghost, 1
      pos( :, npart_real + a ) = postmp( :, npart_real_old + a )
      IF( PRESENT(nu) ) nu( a )= nutmp( npart_real_old + a )
      h( a )= htmp( npart_real_old + a )
      nstar_p( a )= nstarptmp( npart_real_old + a )
      IF( ALLOCATED(all_pos_prev_dump) )THEN
        all_pos_prev_dump( :, npart_real + a )= &
                                all_pos_prev_dump_tmp( :, npart_real_old + a )
      ENDIF
    ENDDO
    !$OMP END PARALLEL DO

    IF( PRESENT(verbose) .AND. verbose .EQV. .TRUE. )THEN

      CALL COM( npart_real, pos(:,1:npart_real), nu(1:npart_real), & ! input
                com_x, com_y, com_z, com_d ) ! output

      PRINT *, "** After mirroring particles:"
      IF( PRESENT(com_star) ) PRINT *, &
               " * x coordinate of the center of mass of the star, ", &
               "from LORENE: com_star= ", com_star, "Msun_geo"
      PRINT *, " * x coordinate of the center of mass of the particle ", &
               "distribution: com_x= ", com_x, "Msun_geo"
      PRINT *, " * y coordinate of the center of mass of the particle ", &
               "distribution: com_y= ", com_y, "Msun_geo"
      PRINT *, " * z coordinate of the center of mass of the particle ", &
               "distribution: com_z= ", com_z, "Msun_geo"
      PRINT *, " * Distance of the center of mass of the particle ", &
               "distribution from the  origin: com_d= ", com_d
      IF( PRESENT(com_star) ) PRINT *, " * |com_x-com_star/com_star|=", &
               ABS( com_x-com_star )/ABS( com_star + 1 )
      PRINT *

    ENDIF

  END SUBROUTINE impose_equatorial_plane_symmetry_apm


  SUBROUTINE correct_center_of_mass( npart_real, pos, nu, get_density, &
                                     validate_pos, com_star, verbose )

    !***********************************************************
    !
    !# Translate the particles so that their center of mass
    !  coincides with the center of mass of the star, given by
    !  LORENE
    !
    !  FT 1.09.2021
    !
    !***********************************************************

    USE analyze, ONLY: COM

    IMPLICIT NONE

    INTEGER, INTENT(IN):: npart_real
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: com_star
    LOGICAL, INTENT(IN), OPTIONAL:: verbose

    INTERFACE
      FUNCTION get_density( x, y, z ) RESULT( density )
        DOUBLE PRECISION, INTENT(IN):: x
        DOUBLE PRECISION, INTENT(IN):: y
        DOUBLE PRECISION, INTENT(IN):: z
        DOUBLE PRECISION:: density
      END FUNCTION
    END INTERFACE
    INTERFACE
      FUNCTION validate_pos( x, y, z ) RESULT( answer )
        DOUBLE PRECISION, INTENT(IN):: x
        DOUBLE PRECISION, INTENT(IN):: y
        DOUBLE PRECISION, INTENT(IN):: z
        INTEGER:: answer
      END FUNCTION
    END INTERFACE

    DOUBLE PRECISION, DIMENSION(3,npart_real), INTENT(INOUT):: pos
    DOUBLE PRECISION, DIMENSION(npart_real),   INTENT(INOUT):: nu

    INTEGER:: a
    DOUBLE PRECISION:: com_x, com_y, com_z, com_d
    DOUBLE PRECISION, DIMENSION(3):: pos_corr_tmp

    CALL COM( npart_real, pos, nu, &       ! input
              com_x, com_y, com_z, com_d ) ! output

    IF( PRESENT(verbose) .AND. verbose .EQV. .TRUE. )THEN
      PRINT *, "** Before center of mass correction:"
      PRINT *, " * x coordinate of the center of mass of the star, ", &
               "from the ID: com_star= ", com_star(1), "Msun_geo"
      PRINT *, " * y coordinate of the center of mass of the star, ", &
               "from the ID: com_star= ", com_star(2), "Msun_geo"
      PRINT *, " * z coordinate of the center of mass of the star, ", &
               "from the ID: com_star= ", com_star(3), "Msun_geo"
      PRINT *, " * x coordinate of the center of mass of the particle ", &
               "distribution: com_x= ", com_x, "Msun_geo"
      PRINT *, " * y coordinate of the center of mass of the particle ", &
               "distribution: com_y= ", com_y, "Msun_geo"
      PRINT *, " * z coordinate of the center of mass of the particle ", &
               "distribution: com_z= ", com_z, "Msun_geo"
      PRINT *, " * Distance of the center of mass of the particle ", &
               "distribution from the  origin: com_d= ", com_d
      PRINT *, " * |com_x-com_star_x/com_star_x|=", &
               ABS( com_x-com_star(1) )/ABS( com_star(1) + 1 )
      PRINT *
    ENDIF

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( pos, com_star, &
    !$OMP                     com_x, com_y, com_z, npart_real ) &
    !$OMP             PRIVATE( pos_corr_tmp, a )
    DO a= 1, npart_real, 1

      pos_corr_tmp(1)= pos(1,a) - ( com_x - com_star(1) )
      pos_corr_tmp(2)= pos(2,a) - ( com_y - com_star(2) )
      pos_corr_tmp(3)= pos(3,a) - ( com_z - com_star(3) )

      IF( get_density( &
                  pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) > zero &
          .AND. &
          validate_pos( &
                  pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) == 0 &
      )THEN

        pos(:,a)= pos_corr_tmp

      ENDIF

    ENDDO
    !$OMP END PARALLEL DO

    CALL COM( npart_real, pos, nu, & ! input
              com_x, com_y, com_z, com_d ) ! output

    IF( PRESENT(verbose) .AND. verbose .EQV. .TRUE. )THEN
      PRINT *, "** After center of mass correction:"
      PRINT *, " * x coordinate of the center of mass of the star, ", &
               "from the ID: com_star= ", com_star(1), "Msun_geo"
      PRINT *, " * y coordinate of the center of mass of the star, ", &
               "from the ID: com_star= ", com_star(2), "Msun_geo"
      PRINT *, " * z coordinate of the center of mass of the star, ", &
               "from the ID: com_star= ", com_star(3), "Msun_geo"
      PRINT *, " * x coordinate of the center of mass of the particle ", &
               "distribution: com_x= ", com_x, "Msun_geo"
      PRINT *, " * y coordinate of the center of mass of the particle ", &
               "distribution: com_y= ", com_y, "Msun_geo"
      PRINT *, " * z coordinate of the center of mass of the particle ", &
               "distribution: com_z= ", com_z, "Msun_geo"
      PRINT *, " * Distance of the center of mass of the particle ", &
               "distribution from the  origin: com_d= ", com_d
      PRINT *, " * |com_x-com_star_x/com_star_x|=", &
               ABS( com_x-com_star(1) )/ABS( com_star(1) + 1 )
      PRINT *
    ENDIF

  END SUBROUTINE correct_center_of_mass


  SUBROUTINE get_neighbours_bf(ipart,npart,pos,h,dimensions,nnei,neilist)

    !**************************************************************
    !
    !# just for test purposes: get neighbours of particle ipart in
    !  a "brute force" way; ipart is ALSO on the neighbour list;
    !  SKR 8.2.2010
    !
    !  Removed ipart from its own neighbors' list
    !  FT 04.06.2021
    !
    !**************************************************************

    IMPLICIT NONE

    INTEGER,INTENT(IN)::          ipart,npart,dimensions
    DOUBLE PRECISION,INTENT(IN):: pos(dimensions,npart),h(npart)
    INTEGER,INTENT(OUT)::         nnei,neilist(npart)
    INTEGER a
    DOUBLE PRECISION diff(dimensions),d2,r_int2

    ! square of interaction radius
    r_int2= (two*h(ipart))**2

    nnei= 0
    !$OMP PARALLEL DO SHARED(pos,dimensions,ipart,npart,r_int2,nnei,neilist)&
    !$OMP             PRIVATE(a,diff,d2)
    DO a= 1, npart, 1

      IF( a /= ipart )THEN

        diff= pos(1:dimensions,a)-pos(1:dimensions,ipart)
        d2= DOT_PRODUCT(diff,diff)

        ! neighbour?
        IF(d2 < r_int2)THEN
          nnei= nnei + 1
          neilist(nnei)= a
        ENDIF

      ENDIF

    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE get_neighbours_bf

  
END SUBMODULE apm
