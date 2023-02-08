!# File:         submodule_sph_particles_ellipsoidal_surfaces.f90
!  Authors:      Francesco Torsello (FT)

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

SUBMODULE (sph_particles) ellipsoidal_surfaces

  !************************************************
  !
  !# This SUBMODULE contains the implementation
  !  of the method of TYPE sph_particles that places
  !  particles on ellipsoidal surfaces inside
  !  a star
  !
  !  FT 19.04.2021
  !
  !************************************************


  IMPLICIT NONE

  ! Be careful! if you define quantities here, they will be global
  ! If you call the SUBROUTINES multiple times, they will use the SAME variables

  !PRIVATE


  CONTAINS


  MODULE PROCEDURE place_particles_ellipsoidal_surfaces

    !**********************************************
    !
    !# Places particles on spherical surfaces
    !  inside a star
    !
    !  FT 19.04.2021
    !
    !  Upgraded to ellipsoidal surfaces. The user can
    !  choose to place prticles on ellipsoidal or
    !  spherical surfaces
    !
    !  FT 15.11.2022
    !
    !**********************************************

    USE constants,  ONLY: pi, half, amu, Msun
    USE utility,    ONLY: zero, one, two, three, five, seven, ten, &
                          cartesian_from_spherical, is_finite_number, &
                          g2kg, m2cm, lorene2hydrobase
    USE matrix,     ONLY: determinant_4x4_matrix
    USE NR,         ONLY: indexx
    USE APM,        ONLY: assign_h
    !USE options,    ONLY: ndes
    !USE units,      ONLY: umass

    IMPLICIT NONE

    INTEGER:: n_surfaces, itr2, cnt, &
              r, th, phi, i_shell, npart_test, npart_shell_tmp, &
              cnt2, rel_sign, dim_seed, r_cnt, prev_shell, first_r, &
              npart_discard, npart_shell_cnt, size_pos_shell

    INTEGER, DIMENSION(:), ALLOCATABLE:: mass_profile_idx, seed
    INTEGER, DIMENSION(:), ALLOCATABLE:: npart_shell, npart_shelleq

    DOUBLE PRECISION:: a_x, a_y, a_z, max_radius, max_center, &
                       xtemp, ytemp, ztemp, m_p, &
                       dr, dth, dphi, phase, phase_th, mass, &
                       dr_shells, dth_shells, dphi_shells, col, long, rad, &
                       proper_volume, mass_test, mass_test2,&
                       proper_volume_test, npart_shell_kept, &
                       rand_num, rand_num2, delta_r, shell_thickness, &
                       upper_bound_tmp, lower_bound_tmp, col_tmp!, &
                       !rho_to_be_resolved

    DOUBLE PRECISION, PARAMETER:: huge_real= 1.0D30!ABS( HUGE(0.0D0) )

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: mass_profile
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: surface_radii, surface_masses, &
                                                  alpha, m_parts, surface_vol, &
                                                  surface_vol2, mass_surface, &
                                                  mass_surface2, shell_scales

    LOGICAL:: exist, high_mass, low_mass, kept_all

    CHARACTER(LEN=:), ALLOCATABLE:: finalnamefile, finalnamefile2

    TYPE:: colatitude_pos_shell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: colatitudes
    END TYPE

    TYPE(colatitude_pos_shell), DIMENSION(:), ALLOCATABLE:: colatitude_pos

    TYPE:: pos_on_surfaces
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos_shell
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pos_th
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pos_phi
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: psurface_vol
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: psurface_vol2
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: sqdetg
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: baryon_density
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: gamma_euler
    END TYPE

    TYPE(pos_on_surfaces), DIMENSION(:), ALLOCATABLE:: pos_surfaces

    INTEGER,          DIMENSION(:,:),   ALLOCATABLE:: npart_surface_tmp
    INTEGER,          DIMENSION(:,:),   ALLOCATABLE:: npart_discarded
    DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE:: pos_shell_tmp
    DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE:: sqdetg_tmp
    DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE:: bar_density_tmp
    DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE:: gam_euler_tmp
    DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE:: pvol_tmp

    CHARACTER(LEN=:), ALLOCATABLE:: surface_type

    LOGICAL, PARAMETER:: debug= .FALSE.

    PRINT *, "** Executing the place_particles_ellipsoidal_surfaces SUBROUTINE..."
    PRINT *

    CALL RANDOM_SEED( SIZE= dim_seed )
    ALLOCATE( seed( dim_seed ) )
    seed( 1 )= 0
    seed( 2 )= 1
    DO itr= 3, dim_seed
      seed( itr )= seed( itr - 1 ) + seed( itr - 2 )
    ENDDO
    CALL RANDOM_SEED( PUT= seed )

    IF(debug) PRINT *, dim_seed
    IF(debug) PRINT *, seed

    IF(PRESENT(radii))THEN

      surface_type="ellipsoidal"

      max_radius= MAXVAL([radius,radii(1),radii(2)])

      IF(max_radius == radius)THEN

        a_x= one
        a_y= radii(1)/max_radius
        a_z= radii(2)/max_radius
        max_center= center(1)

      ELSEIF(max_radius == radii(1))THEN

        a_x= radius/max_radius
        a_y= one
        a_z= radii(2)/max_radius
        max_center= center(2)

      ELSEIF(max_radius == radii(2))THEN

        a_x= radius/max_radius
        a_y= radii(1)/max_radius
        a_z= one
        max_center= center(3)

      ENDIF

    ELSE

      surface_type="spherical"

      max_radius= radius
      a_x= one
      a_y= one
      a_z= one
      max_center= center(1)

    ENDIF

    !-----------------------------------!
    !-- Compute desired particle mass --!
    !-----------------------------------!

    !IF( PRESENT(pmass_des) )THEN
    !  m_p= pmass_des
    !ELSE
      m_p= mass_star/DBLE(npart_des)
    !ENDIF

    !PRINT *, "std m_p             =", m_p
    !PRINT *, "std npart           =", npart_des
    !rho_to_be_resolved= 1.D+13*(g2kg*m2cm**3)*lorene2hydrobase!1.0D-5
    !!rho_to_be_resolved= 1.138065390333111E-004
    !m_p= rho_to_be_resolved/(ndes/two)
    !!PRINT *, "average rho_id on the outer layers for 1.3 MPA1=", &
    !!         1.363246651849556D+53*amu/umass !=1.138065390333111E-004
    !PRINT *, "rho to be resolved  =", rho_to_be_resolved
    !PRINT *, "m_p to resolve rho  =", m_p
    !PRINT *, "npart to resolve rho=", NINT(mass_star/m_p)
    !STOP

    !--------------------------------!
    !-- Compute number of surfaces --!
    !--------------------------------!

    n_surfaces= number_surfaces( m_p, max_center, max_radius, get_density )

    !--------------------------------------!
    !-- Allocate memory for the surfaces --!
    !--------------------------------------!

    CALL allocate_surface_memory()

    !--------------------------------------------------------!
    !-- Place surfaces based on mass density at that point --!
    !--------------------------------------------------------!

    CALL place_surfaces( central_density, max_center, max_radius, m_p, &
                         n_surfaces, &
                         surface_radii, last_r, get_density )

    ! Printout
    PRINT *, " * Number of the ", surface_type, " surfaces= ", n_surfaces
    PRINT *, " * Radii of the surfaces in units of the equatorial radius", &
             " of the star, towards the companion= "
    PRINT *, surface_radii/max_radius
    PRINT *

    !---------------------------------!
    !-- Compute radial mass profile --!
    !---------------------------------!

    PRINT *, " * Integrating the baryon mass density to get the mass profile..."
    PRINT *

    dr  = max_radius/(three*ten*ten)
    dth = pi/two/(two*ten*ten)
    dphi= two*pi/(two*ten*ten)

    ALLOCATE( mass_profile( 3, 0:NINT(max_radius/dr) ), &
              STAT= ios, ERRMSG= err_msg )
    ALLOCATE( mass_profile_idx( 0:NINT(max_radius/dr) ),&
              STAT= ios, ERRMSG= err_msg )

    CALL integrate_density( center, max_radius, &
                            central_density, &
                            dr, dth, dphi, &
                            mass, mass_profile, &
                            mass_profile_idx, radii )

    mass_profile( 2:3, : )= mass_profile( 2:3, : )*mass_star/mass

    !-----------------------------------!
    !-- Assign masses to each surface --!
    !-----------------------------------!

    CALL assign_surfaces_mass( surface_masses, surface_radii, max_radius, dr, &
                               n_surfaces, mass_profile_idx, mass_profile, &
                               mass_star )

    !----------------------------------------------------!
    !-- Print mass profile and surfaces' radii to file --!
    !----------------------------------------------------!

    IF( PRESENT(filename_mass_profile) )THEN
      finalnamefile= filename_mass_profile
    ELSE
      finalnamefile= "mass_profile.dat"
    ENDIF

    IF( PRESENT(filename_shells_radii) )THEN
      finalnamefile2= filename_shells_radii
    ELSE
      finalnamefile2= "surface_radii.dat"
    ENDIF

    CALL print_mass_profile_surface_radii( mass_profile, mass_profile_idx, &
                                           surface_radii, max_radius, dr, &
                                           n_surfaces, &
                                           filename_mass_profile, &
                                           filename_shells_radii )

    !---------------------------------------------------------!
    !-- Initialize quantities before starting the iteration --!
    !---------------------------------------------------------!

    PRINT *, " * Initializing quantities before starting the iteration..."
    PRINT *

    CALL initialize_surfaces()

    !----------------------------------------!
    !--  Main iteration over the surfaces  --!
    !----------------------------------------!

    PRINT *, " * Assigning first half of particle positions..."
    PRINT *

    place_particles_on_northern_emispheres: DO

      ! Correct npart_shelleq to be divisible by 4
      IF( MOD( npart_shelleq(r), 2 ) /= 0 )THEN
        CALL RANDOM_NUMBER( rand_num2 )
        IF( rand_num2 >= half ) rel_sign=  1
        IF( rand_num2 < half )  rel_sign= -1
        npart_shelleq(r)= npart_shelleq(r) + rel_sign
      ENDIF
      IF( MOD( npart_shelleq(r)/2, 2 ) /= 0 )THEN
        CALL RANDOM_NUMBER( rand_num2 )
        IF( rand_num2 >= half ) rel_sign=  1
        IF( rand_num2 < half )  rel_sign= -1
        npart_shelleq(r)= 2*( npart_shelleq(r)/2 + rel_sign )
      ENDIF
      IF( MOD( npart_shelleq(r), 4 ) /= 0 )THEN
        PRINT *, " * ERROR! npart_shelleq(", r, ")=", npart_shelleq(r), &
                 " is not divisible by 4 in the main iteration. ", &
                 " Check the algorithm. Stopping..."
        STOP
      ENDIF

      ! Compute number of particles on the surface
      npart_shell(r)= NINT(( npart_shelleq(r)**two )/two)

      ! Compute angular step in azimuth phi (constant on each shell)
      IF( npart_shelleq(r) == 0 )THEN
        alpha(r)= zero
      ELSE
        alpha(r)= two*pi/DBLE(npart_shelleq(r))
      ENDIF

      ! Compute angular positions in colatitude theta,
      ! according to https://mathworld.wolfram.com/SpherePointPicking.html
      IF( ALLOCATED( colatitude_pos(r)% colatitudes ) ) &
        DEALLOCATE( colatitude_pos(r)% colatitudes )
      ALLOCATE( colatitude_pos(r)% colatitudes( npart_shelleq(r)/4 ) )

      IF( surface_radii(r) < (one - five/(ten*ten))*last_r*max_radius )THEN

        CALL compute_colatitudes_uniformly_in( pi/two, (ten - five/ten)/ten*pi,&
                                      colatitude_pos(r)% colatitudes( : ) )

      ELSE

        CALL compute_colatitudes_uniformly_in( pi/two, (ten - five/ten)/ten*pi,&
                                      colatitude_pos(r)% colatitudes( : ) )
        !CALL compute_colatitudes_uniformly_in( pi/two, two/3.0D0*pi, &
        !                              colatitude_pos(r)% colatitudes( : ) )

      ENDIF

          !            alpha(r)*one/two + ( itr2 - 1 )*alpha(r)

        !  ACOS( two*( one - COS( pi/3.0D0*( two/3.0D0 + DBLE(itr2 - 1)*DBLE(npart_shelleq(r)/4 + one -(one/two)-(two/3.0D0) )/DBLE(npart_shelleq(r)/4 - one ) ) &
        !                   /DBLE(npart_shelleq(r)/4 + one ) ) ) &
        !      - one )
              !5.0D0/1two

        !colatitude_pos(r)% colatitudes( itr2 )= &
        !              colatitude_pos(r)% colatitudes( itr2 ) &
        !              *( 1 + rel_sign*0.05D0*phase_th )
      DO itr2= 1, npart_shelleq(r)/4, 1

        IF( colatitude_pos(r)% colatitudes( itr2 ) <= pi/two .OR. &
            colatitude_pos(r)% colatitudes( itr2 ) >= pi &
        )THEN
          PRINT *, "** ERROR! ", &
                   "The colatitudes are not in the OPEN interval (pi/2,pi). ", &
                   "Stopping..."
          STOP
        ENDIF

      ENDDO

      npart_discard    = 0
      npart_shell_cnt  = 0
      npart_shell_tmp  = npart_shell(r)
      ! Initialize te  mporary arrays
      pos_shell_tmp    = huge_real
      sqdetg_tmp       = zero
      bar_density_tmp  = zero
      gam_euler_tmp    = zero
      pvol_tmp         = zero
      npart_discarded  = zero
      npart_surface_tmp= zero

      IF( debug ) PRINT *, "Right before OMP, shell ", r, "iteration ", cnt2 + 1

      dphi_shells= alpha(r)

      !$OMP PARALLEL DO DEFAULT(NONE), &
      !$OMP             PRIVATE( phase, col, col_tmp, xtemp, ytemp, ztemp, &
      !$OMP                      dth_shells, delta_r, long, &
      !$OMP                      th, phi, rand_num2, phase_th, rel_sign ), &
      !$OMP             SHARED( r, npart_shelleq, center, rad, alpha, &
      !$OMP                     pos_surfaces, colatitude_pos, n_surfaces, &
      !$OMP                     dr_shells, surface_radii, shell_thickness, &
      !$OMP                     a_x, a_y, a_z, this, &
      !$OMP                     sqdetg_tmp, bar_density_tmp, gam_euler_tmp, &
      !$OMP                  pos_shell_tmp, pvol_tmp, dphi_shells, max_radius, &
      !$OMP                     npart_discarded, npart_surface_tmp, last_r )
      DO phi= 1, npart_shelleq(r), 1

        IF( debug ) PRINT *, "Right before loop over phi"

        DO th= 1, npart_shelleq(r)/4, 1 !npart_shelleq(r) is even, see above

          !
          !-- Randomize positions, if specified by the user in the
          !-- parameter file lorene_bns_id_particles.par
          !
          IF( this% randomize_phi )THEN

            CALL RANDOM_NUMBER( phase )
            phase= phase*alpha(r)

          ENDIF

        !  IF( surface_radii(r) < 0.95D0*last_r*max_radius )THEN
        !
        !    long= phase + phi*alpha(r)
        !
        !  ELSE
        !
        !    long= phase + phi*alpha(r)/3.0D0 - pi/3.0D0
        !
        !  ENDIF
          long= phase + phi*alpha(r)


          col= colatitude_pos(r)% colatitudes(th)
          IF( this% randomize_theta )THEN

            CALL RANDOM_NUMBER( phase_th )
            CALL RANDOM_NUMBER( rand_num2 )
            IF( rand_num2 >= half ) rel_sign=  1
            IF( rand_num2 < half )  rel_sign= -1

            col_tmp= col*( one + rel_sign*five/ten/ten*phase_th )

            IF( col_tmp < pi .AND. col_tmp > pi/two )THEN

              col= col_tmp

            ENDIF

          ENDIF

          rad= surface_radii(r)
          IF( this% randomize_r )THEN

            CALL RANDOM_NUMBER( delta_r )
            CALL RANDOM_NUMBER( rand_num2 )
            IF( rand_num2 >= half ) rel_sign=  1
            IF( rand_num2 < half )  rel_sign= -1

            IF( r/n_surfaces < (one - five/(ten*ten)) )THEN

              rad= rad + rel_sign*delta_r*(seven*five/(ten*ten))*dr_shells

            ELSE

              !rad= rad - (one + delta_r)*0.35D0*dr_shells

              rad= rad - (delta_r*(seven*five/(ten*ten)) + five/ten)*dr_shells

              !rad= rad - delta_r*(seven*five/(ten*ten))*dr_shells

            ENDIF

          ENDIF

          IF( rad < 0 )THEN
            PRINT *, " * ERROR! rad < 0. Check the computation of the radial", &
                     " coordinates of the particles. Stopping.."
            STOP
          ENDIF

          !
          !-- Compute Cartesian coordinates of the candidate particle positions
          !
          CALL cartesian_from_spherical( &
            a_x*rad, col, long, &
            center(1), center(2), center(3), &
            xtemp, ytemp, ztemp, a_y/a_x, a_z/a_x )

          IF( .NOT.is_finite_number( xtemp ) )THEN
            PRINT *, "** ERROR when placing first half of the particles! ", &
                     "xtemp is a NaN. Stopping.."
            STOP
          ENDIF
          IF( .NOT.is_finite_number( ytemp ) )THEN
            PRINT *, "** ERROR when placing first half of the particles! ", &
                     "ytemp is a NaN. Stopping.."
            STOP
          ENDIF
          IF( .NOT.is_finite_number( ztemp ) )THEN
            PRINT *, "** ERROR when placing first half of the particles! ", &
                     "ztemp is a NaN. Stopping.."
            STOP
          ENDIF

          ! Import ID needed to compute the particle masses
          CALL get_id( xtemp, ytemp, ztemp, &
                       sqdetg_tmp( th, phi ), &
                       bar_density_tmp( th, phi ), &
                       gam_euler_tmp( th, phi ) )

          ! Place a particle at a given position only if the hydro
          ! is acceptable
          IF( validate_position_final( xtemp, ytemp, ztemp ) )THEN

            npart_surface_tmp( th, phi )= 1
            pos_shell_tmp( 1, th, phi )= xtemp
            pos_shell_tmp( 2, th, phi )= ytemp
            pos_shell_tmp( 3, th, phi )= ztemp

            ! Compute particle volume
            pvol_tmp( th, phi )= particle_volume( rad, col, dr_shells, &
                                               dth_shells, dphi_shells, th, &
                                               colatitude_pos(r)% colatitudes, &
                                               npart_shelleq(r) )

            ! Safety check
            IF( pvol_tmp( th, phi ) <= 0 )THEN
                ! pos_surfaces(r)% psurface_vol2( itr + 1 ) <= 0 )THEN
              PRINT *, "When placing first half of particles"
              PRINT *, "pvol_tmp( ", r, ",", th, ",", phi, " ) =", &
                       pvol_tmp( th, phi )
              PRINT *, "dr_shells=", dr_shells, &
                       "dth_shells=", dth_shells, &
                       "dphi_shells=", dphi_shells
              STOP
            ENDIF

          ELSE

            ! If the hydro is not positive, or the position is outside the star,
            ! discard the position and count the number of discarded positions
            !npart_discard= npart_discard + 2
            npart_discarded( th, phi )= 2

          ENDIF

        ENDDO

      ENDDO
      !$OMP END PARALLEL DO

      npart_discard   = SUM( SUM( npart_discarded, DIM= 1 ), DIM= 1 )
      npart_shell_cnt = SUM( SUM( npart_surface_tmp, DIM= 1 ), DIM= 1 )
      npart_shell(r)  = MAX( npart_shell(r) - npart_discard, 0 )
      npart_out       = npart_out + npart_shell(r)/2

    !  PRINT *, 'r=', r
    !  PRINT *, "npart_discard=", npart_discard
    !  PRINT *, "npart_shell_cnt=", npart_shell_cnt
    !  PRINT *, 'npart_shell(r)=', npart_shell(r)
    !  PRINT *
    !  STOP

      IF( debug ) PRINT *, "Right after OMP"

      ! Safety check
      IF( npart_shell_cnt /= npart_shell(r)/2 )THEN
        PRINT *, "** ERROR! Mismatch in the particle counters on shell ", r
        PRINT *, " * npart_shell_cnt=", npart_shell_cnt, &
                 "npart_shell(r)/2=", npart_shell(r)/2
        PRINT *, " * npart_shell_cnt should be equal to npart_shell(r)/2. " &
                 // "Stopping..."
        PRINT *
        STOP
      ENDIF

      ! Set up the next step in pathological cases
      IF( npart_shell(r) < 0 ) npart_shell(r)= 0
      IF( npart_shell(r) == 0 )THEN
        IF( r == first_r )THEN
          PRINT *, " ** ERROR! No particles were placed on the first ", &
                   "", surface_type, &
                   " surface! Maybe the system has a geometry such ", &
                   "that there is no matter on the ", surface_type, &
                   " surface with ", &
                   "radius r ~ R/2, with R its radial size?"
          PRINT *, "    If so, the first ", surface_type, &
                   " surface to be populated ", &
                   "has to be changed by changing the variable first_r ", &
                   "in the (internal) SUBROUTINE initialize_surfaces in ", &
                   "SUBMODULE sph_particles@ellipsoidal_surfaces."
          PRINT *, "    Another possibility is that the FUNCTION ", &
                   "validate_position given as argument to the SUBROUTINE ", &
                   "place_particles_ellipsoidal_surfaces, does not return ", &
                   ".TRUE. when a position is acceptable, and .FALSE. when ", &
                   "it is not."
          PRINT *, "    Yet another possibility is that the (internal) ", &
                   "FUNCTION validate_position_final is not set up ", &
                   "properly in SUBMODULE sph_particles@ellipsoidal_surfaces."
          PRINT *, " * Stopping..."
          PRINT *
          STOP
        ENDIF
        m_parts(r)= m_parts( prev_shell )
        PRINT *, " * Placed", npart_shell(r)/2, &
                 " particles on one emisphere of ", surface_type, " surface ", r, &
                 " out of ", n_surfaces
        IF( r == 1 )THEN
          EXIT
        ELSEIF( r < first_r )THEN
          !PRINT *, "r=", r
          r= r - 1
          cnt2 = 0
          upper_bound_tmp= upper_bound
          lower_bound_tmp= lower_bound
          CYCLE
        ELSEIF( r == n_surfaces )THEN
          !PRINT *, "r=", r
          r= first_r - 1
          r_cnt= r_cnt + 1
          cnt2 = 0
          upper_bound_tmp= upper_bound
          lower_bound_tmp= lower_bound
          CYCLE
        ELSEIF( r >= first_r )THEN
          !PRINT *, "r=", r
          r= r + 1
          cnt2 = 0
          upper_bound_tmp= upper_bound
          lower_bound_tmp= lower_bound
          CYCLE
        ENDIF
      ELSE
        m_parts(r)= surface_masses(r)/DBLE(npart_shell(r))
      ENDIF

      IF( debug ) PRINT *, " * Before storing the particles"
      IF( debug ) PRINT *, "11"

      IF( debug ) PRINT *, "Right before correction of particle number"
      IF( debug ) PRINT *, "npart_out=", npart_out

      ! If it's not the first populated surface
      not_first_populated_surface: IF( r /= first_r )THEN

        ! Identify the previous surface
        IF( r < first_r )THEN
          prev_shell= r + 1
        ELSEIF( r > first_r )THEN
          prev_shell= r - 1
        ELSEIF( r == 1 )THEN
          EXIT
        ENDIF

        ! Logical variables that steer the iteration

        ! This speeds up the iteration considerably, since the inner layers
        ! seem to always want a larger tolerance
        IF( r <= (five*three*three/(ten*ten))*n_surfaces )THEN
          upper_bound_tmp= upper_bound_tmp*(one + one/ten)
          lower_bound_tmp= lower_bound_tmp*(one - one/ten)
        ENDIF

        ! Is the particle mass too high?
        high_mass= m_parts(r)/m_parts( prev_shell ) > upper_bound_tmp
        ! Is the particle mass too low?
        low_mass = m_parts(r)/m_parts( prev_shell ) < lower_bound_tmp
        ! How many positions were kept, placing a particle on them?
        npart_shell_kept= DBLE(npart_shell(r))/DBLE(npart_shell_tmp)
        ! Were all the positions kept?
        kept_all = npart_shell_kept == one

        ! If the particle mass is too high and all positions were kept
        adjust_particle_number_surface: IF( high_mass .AND. kept_all )THEN

          IF( debug ) PRINT *, "Case 1"

          cnt2= cnt2 + 1

          ! If this is the (max_steps + 1)st step
          too_many_steps: IF( cnt2 > max_steps )THEN

            ! Allow for a bit more different particle mass
            upper_bound_tmp= upper_bound_tmp*upper_factor
            lower_bound_tmp= lower_bound_tmp*lower_factor

            !
            !-- Special treatment for the positions near the surface
            !

            ! If the range of particle masses is getting too generous
            ! near the surface
            IF( r > (one - two/ten)*n_surfaces .AND. &
              m_parts(r)/m_parts( prev_shell ) > (one + one/ten)*upper_bound &
            )THEN

              ! Increase the number of positions on this surface
              ! by a factor between 5 and 10
              CALL RANDOM_NUMBER( rand_num2 )
              rand_num= NINT( five*( rand_num2 + one ) )
              npart_shelleq(r)= NINT( rand_num*npart_shelleq(r) )
                                !  npart_shelleq( r - 1 ) &
                                !+ rel_sign*NINT( 1 + rand_num )

              ! Reset the particle mass tolerance range
              upper_bound_tmp= upper_bound
              lower_bound_tmp= lower_bound

              ! Reset total particle number
              npart_out= npart_out - npart_shell(r)/2

              ! Reset counter
              cnt2= 1

              ! Replace particles on this surface
              CYCLE

            ENDIF

            ! Reset counter
            cnt2= 1

          ENDIF too_many_steps

          ! If this is not yet the (max_steps + 1)st step

          ! Reset total particle number
          npart_out= npart_out - npart_shell(r)/2

          ! Increase randomly particle number on the equator which determines
          ! the particle number on the surface
          ! More particles = lower particle mass, with fixed surface mass
          CALL RANDOM_NUMBER( rand_num )
          CALL RANDOM_NUMBER( rand_num2 )
          npart_shelleq(r)= npart_shelleq(r) + 1*NINT( one + one*rand_num )&
                                                 + 1*NINT( one + one*rand_num2 )

          ! Treat pathological cases
          IF( npart_shelleq(r) == 0 .OR. npart_shell(r) == 0 )THEN
            CALL RANDOM_NUMBER( rand_num )
            CALL RANDOM_NUMBER( rand_num2 )
            IF( rand_num2 < half )  rel_sign= - 1
            IF( rand_num2 >= half ) rel_sign=   1
            npart_shelleq(r)= npart_shelleq( r - 1 ) &
                              + rel_sign*NINT( 1 + rand_num )
          ENDIF

          ! Replace particles on this surace
          CYCLE

        ! The cases below do similar things to the one above, so there are no
        ! comments on each line. See the case above for explanations
        ELSEIF( low_mass .AND. kept_all )THEN

          IF( debug ) PRINT *, "Case 2"

          cnt2= cnt2 + 1
          IF( cnt2 > max_steps )THEN
            upper_bound_tmp= upper_bound_tmp*upper_factor
            lower_bound_tmp= lower_bound_tmp*lower_factor
            IF( r > (one - two/ten)*n_surfaces .AND. &
              m_parts(r)/m_parts( prev_shell ) < (one - one/ten)*lower_bound &
            )THEN
              CALL RANDOM_NUMBER( rand_num2 )
              rand_num= NINT( five*( rand_num2 + one ) )
              npart_shelleq(r)= NINT( npart_shelleq(r)/rand_num )
                                !  npart_shelleq( r - 1 ) &
                                !+ rel_sign*NINT( 1 + rand_num )
              upper_bound_tmp= upper_bound
              lower_bound_tmp= lower_bound
              npart_out= npart_out - npart_shell(r)/2
              cnt2= 1
              CYCLE
            ENDIF
            cnt2= 1
          ENDIF

          npart_out= npart_out - npart_shell(r)/2

          CALL RANDOM_NUMBER( rand_num )
          CALL RANDOM_NUMBER( rand_num2 )
          npart_shelleq(r)= npart_shelleq(r) - 1*NINT( one + one*rand_num )&
                                                 - 1*NINT( one + one*rand_num2 )

          IF( npart_shelleq(r) == 0 .OR. npart_shell(r) == 0 )THEN
            CALL RANDOM_NUMBER( rand_num )
            CALL RANDOM_NUMBER( rand_num2 )
            IF( rand_num2 < half )  rel_sign= - 1
            IF( rand_num2 >= half ) rel_sign=   1
            npart_shelleq(r)= npart_shelleq( r - 1 ) &
                              + rel_sign*NINT( 1 + rand_num )
          ENDIF

          CYCLE

        ELSEIF( high_mass .AND. .NOT.kept_all ) THEN

          IF( debug ) PRINT *, "Case 3"

          cnt2= cnt2 + 1
          IF( cnt2 > max_steps )THEN
            upper_bound_tmp= upper_bound_tmp*upper_factor
            lower_bound_tmp= lower_bound_tmp*lower_factor
            IF( r > (one - two/ten)*n_surfaces .AND. &
              m_parts(r)/m_parts( prev_shell ) > (one + one/ten)*upper_bound &
                !upper_bound_tmp > 1.1D0*upper_bound &
            )THEN
              CALL RANDOM_NUMBER( rand_num2 )
              rand_num= NINT( five*( rand_num2 + one ) )
              npart_shelleq(r)= NINT( rand_num*npart_shelleq(r) )
                                !  npart_shelleq( r - 1 ) &
                                !+ rel_sign*NINT( 1 + rand_num )
              upper_bound_tmp= upper_bound
              lower_bound_tmp= lower_bound
              npart_out= npart_out - npart_shell(r)/2
              cnt2= 1
              CYCLE
            ENDIF
            cnt2= 1
          ENDIF

          npart_out= npart_out - npart_shell(r)/2

          CALL RANDOM_NUMBER( rand_num )
          CALL RANDOM_NUMBER( rand_num2 )
          IF( rand_num2 < half )  rel_sign= - 1
          IF( rand_num2 >= half ) rel_sign=   1

          ! If x% of the positions were kept, divide the old particle number
          ! on the equator by x, and adjust with some other random and
          ! non-random factors which turn out to work well a posteriori
          npart_shelleq(r)= CEILING( SQRT( &
                                2*(surface_masses(r)/m_parts( prev_shell )) &
                                /npart_shell_kept &
                              ) ) + rel_sign*NINT( 1 + rand_num )

          IF( npart_shelleq(r) == 0 .OR. npart_shell(r) == 0 )THEN
            CALL RANDOM_NUMBER( rand_num )
            CALL RANDOM_NUMBER( rand_num2 )
            IF( rand_num2 < half )  rel_sign= - 1
            IF( rand_num2 >= half ) rel_sign=   1
            npart_shelleq(r)= npart_shelleq( prev_shell ) &
                              + rel_sign*NINT( 1 + rand_num )
          ENDIF

          CYCLE

        ELSEIF( low_mass .AND. .NOT.kept_all ) THEN

          IF( debug ) PRINT *, "Case 4"

          cnt2= cnt2 + 1
          IF( cnt2 > max_steps )THEN
            upper_bound_tmp= upper_bound_tmp*upper_factor
            lower_bound_tmp= lower_bound_tmp*lower_factor
            IF( r > (one - two/ten)*n_surfaces .AND. &
              m_parts(r)/m_parts( prev_shell ) < (one - one/ten)*lower_bound &
                !lower_bound_tmp < 0.9D0*lower_bound &
            )THEN
              CALL RANDOM_NUMBER( rand_num2 )
              rand_num= NINT( five*( rand_num2 + one ) )
              npart_shelleq(r)= NINT( npart_shelleq(r)/rand_num )
                                !  npart_shelleq( r - 1 ) &
                                !+ rel_sign*NINT( 1 + rand_num )
              upper_bound_tmp= upper_bound
              lower_bound_tmp= lower_bound
              npart_out= npart_out - npart_shell(r)/2
              cnt2= 1
              CYCLE
            ENDIF
            cnt2= 1
          ENDIF

          npart_out= npart_out - npart_shell(r)/2

          CALL RANDOM_NUMBER( rand_num )
          CALL RANDOM_NUMBER( rand_num2 )
          IF( rand_num2 < half )  rel_sign= - 1
          IF( rand_num2 >= half ) rel_sign=   1

          ! If x% of the positions were kept, divide the old particle number
          ! on the equator by x, and adjust with some other random and
          ! non-random factors which turn out to work well a posteriori
          npart_shelleq(r)= CEILING( SQRT( &
                                2*(surface_masses(r)/m_parts( prev_shell )) &
                                /npart_shell_kept &
                              ) ) + rel_sign*NINT( 1 + rand_num )

          IF( npart_shelleq(r) == 0 .OR. npart_shell(r) == 0 )THEN
            CALL RANDOM_NUMBER( rand_num )
            CALL RANDOM_NUMBER( rand_num2 )
            IF( rand_num2 < half )  rel_sign= - 1
            IF( rand_num2 >= half ) rel_sign=   1
            npart_shelleq(r)= npart_shelleq( prev_shell ) &
                              + rel_sign*NINT( 1 + rand_num )
          ENDIF

          IF( debug ) PRINT *, "Right after correction of particle number"

          ! Replace particle on this surface
          CYCLE

        ENDIF adjust_particle_number_surface

      ENDIF not_first_populated_surface

      ! >TODO: Safety check to be updated
      !IF( r_cnt > 1 )THEN
      !  npart_test= 0
      !  DO itr= 1, r, 1
      !    npart_test= npart_test + npart_shell( itr )
      !  ENDDO
      !  IF( npart_test/2 /= npart_out )THEN
      !    PRINT *, "** ERROR! The sum of the particles on the shells is not ", &
      !             "equal to the total number of particles. Stopping.."
      !    PRINT *, "npart_test=", npart_test/2, ", npart_out=", npart_out
      !    PRINT *
      !    STOP
      !  ENDIF
      !ENDIF

      IF( debug ) PRINT *, "10"

      ! At this point, the particles are placed on this surface
      ! Print out the result
      PRINT *, " * Placed", npart_shell(r)/2, &
               " particles on one emisphere of ", surface_type, " surface ", r, &
               " out of ", n_surfaces
      PRINT *, "   Surface radius= ", surface_radii(r)/max_radius*ten*ten, &
              "% of the smaller radius of the matter object"
      PRINT *, "   Placed", npart_out, " particles overall, so far."
      IF( r /= first_r ) PRINT *, &
               "   Ratio of particle masses on last 2 surfaces: ", &
               "   m_parts(", r, ")/m_parts(", prev_shell, ")= ",  &
               m_parts(r)/m_parts( prev_shell )

      ! Save particles to non-temporary variables
      size_pos_shell= SIZE( pos_surfaces(r)% pos_shell( 1, : ) )

      IF( npart_shelleq(r)*npart_shelleq(r)/4 > size_pos_shell )THEN

        CALL reallocate_array_2d( pos_surfaces(r)% pos_shell, 3, &
                      npart_shelleq(r)*npart_shelleq(r)/4 )

        CALL reallocate_array_1d( pos_surfaces(r)% sqdetg, &
                      npart_shelleq(r)*npart_shelleq(r)/4 )

        CALL reallocate_array_1d( pos_surfaces(r)% baryon_density, &
                      npart_shelleq(r)*npart_shelleq(r)/4 )

        CALL reallocate_array_1d( pos_surfaces(r)% gamma_euler, &
                      npart_shelleq(r)*npart_shelleq(r)/4 )

        CALL reallocate_array_1d( pos_surfaces(r)% psurface_vol2, &
                      npart_shelleq(r)*npart_shelleq(r)/4 )

      ENDIF

      itr= 0
      DO th= 1, npart_shelleq(r)/4, 1
        DO phi= 1, npart_shelleq(r), 1

          IF( pos_shell_tmp( 1, th, phi ) /= huge_real &
              !pos_shell_tmp( 1, th, phi ) < center + 1.2D0*max_radius &
              !.AND. &
              !pos_shell_tmp( 1, th, phi ) > center - 1.2D0*max_radius
          )THEN

            itr= itr + 1

            pos_surfaces(r)% pos_shell( 1, itr )= pos_shell_tmp( 1, th, phi )
            pos_surfaces(r)% pos_shell( 2, itr )= pos_shell_tmp( 2, th, phi )
            pos_surfaces(r)% pos_shell( 3, itr )= pos_shell_tmp( 3, th, phi )

            pos_surfaces(r)% sqdetg( itr )        = sqdetg_tmp( th, phi )
            pos_surfaces(r)% baryon_density( itr )= bar_density_tmp( th, phi )
            pos_surfaces(r)% gamma_euler( itr )   = gam_euler_tmp( th, phi )
            pos_surfaces(r)% psurface_vol2( itr )   = pvol_tmp( th, phi )


          ENDIF

        ENDDO
      ENDDO
      ! Safety check
      IF( npart_shell_cnt /= itr )THEN
        PRINT *, "** ERROR! Mismatch in the particle counters on shell ", r
        PRINT *, " * npart_shell_cnt=", npart_shell_cnt, &
                 ", itr=", itr, ". npart_shell_cnt should be equal to itr. "
        PRINT *, " * npart_shell(r)/2=", npart_shell(r)/2
        STOP
      ENDIF

      ! Set up next step
      IF( r == n_surfaces )THEN
        r= first_r - 1
        r_cnt= r_cnt + 1
        cnt2 = 0
        upper_bound_tmp= upper_bound
        lower_bound_tmp= lower_bound
        IF( debug ) PRINT *, "last shell"
      ELSEIF( r == 1 )THEN
        IF( debug ) PRINT *, "exit"
        EXIT
      ELSEIF( r < first_r )THEN
        r= r - 1
        r_cnt= r_cnt + 1
        cnt2 = 0
        upper_bound_tmp= upper_bound
        lower_bound_tmp= lower_bound
        IF( debug ) PRINT *, "inner layers"
      ELSEIF( r >= first_r )THEN
        r= r + 1
        r_cnt= r_cnt + 1
        cnt2 = 0
        upper_bound_tmp= upper_bound
        lower_bound_tmp= lower_bound
        IF( debug ) PRINT *, "outer layers"
      ENDIF

      IF( debug ) PRINT *, "12"

    ENDDO place_particles_on_northern_emispheres

    !-----------------------------!
    !--  End of main iteration  --!
    !-----------------------------!

    ! Print out the total number of particles on the northern emispheres,
    ! and the final mass ratio
    PRINT *, " * Particles on the northern emispheres=", npart_out
    PRINT *, " * Particle mass ratio= ", MAXVAL(m_parts)/MINVAL(m_parts)
    PRINT *

    ! Safety check
    npart_test= SUM( npart_shell, DIM= 1 )
    IF( npart_test/2 /= npart_out )THEN
      PRINT *, "** ERROR! The sum of the particles on the shells is not ", &
               "equal to the total number of particles. Stopping.."
      PRINT *, " * npart_test/2=", npart_test/2, ", npart_out=", npart_out
      PRINT *, " * Array npart_shell=", npart_shell
      PRINT *
      STOP
    ENDIF

    IF( debug ) PRINT *, "13"

    ! Deallocate temporary arrays
    IF( ALLOCATED(pos_shell_tmp) )THEN
      DEALLOCATE( pos_shell_tmp  , STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array m_parts in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED(sqdetg_tmp) )THEN
      DEALLOCATE( sqdetg_tmp       , STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array m_parts in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED(bar_density_tmp) )THEN
      DEALLOCATE( bar_density_tmp, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array m_parts in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED(gam_euler_tmp) )THEN
      DEALLOCATE( gam_euler_tmp  , STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array m_parts in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED(pvol_tmp) )THEN
      DEALLOCATE( pvol_tmp       , STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array m_parts in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED(npart_discarded) )THEN
      DEALLOCATE( npart_discarded       , STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array npart_discarded in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED(npart_surface_tmp) )THEN
      DEALLOCATE( npart_surface_tmp       , STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array m_parts in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF

    IF( debug ) PRINT *, "14"

    !
    !-- Mirror particles from the northern emispheres to the southern ones
    !
    PRINT *, " * Mirroring particles from the northern emispheres to the", &
             " southern ones..."
    PRINT *

    IF( debug ) PRINT *, " * npart/2=", npart_out

    DO r= 1, n_surfaces, 1

      DO itr= 1, npart_shell(r)/2, 1

        npart_out= npart_out + 1

        pos_surfaces(r)% pos_shell( 1, npart_shell(r)/2 + itr )= &
                                          pos_surfaces(r)% pos_shell( 1, itr )
        pos_surfaces(r)% pos_shell( 2, npart_shell(r)/2 + itr )= &
                                          pos_surfaces(r)% pos_shell( 2, itr )
        pos_surfaces(r)% pos_shell( 3, npart_shell(r)/2 + itr )= &
                                        - pos_surfaces(r)% pos_shell( 3, itr )
        pos_surfaces(r)% sqdetg( npart_shell(r)/2 + itr )= &
                                                 pos_surfaces(r)% sqdetg( itr )
        pos_surfaces(r)% baryon_density( npart_shell(r)/2 + itr )= &
                                        pos_surfaces(r)% baryon_density( itr )
        pos_surfaces(r)% gamma_euler( npart_shell(r)/2 + itr )= &
                                          pos_surfaces(r)% gamma_euler( itr )
        pos_surfaces(r)% psurface_vol2( npart_shell(r)/2 + itr )= &
                                           pos_surfaces(r)% psurface_vol2( itr )

        ! Safety checks
        IF( pos_surfaces(r)% baryon_density( itr ) == 0 )THEN
          PRINT *, "When mirroring particles"
          PRINT *, r, itr, pos_surfaces(r)% pos_shell( 1, itr ), &
                   pos_surfaces(r)% pos_shell( 2, itr ), &
                   pos_surfaces(r)% pos_shell( 3, itr ), &
                   pos_surfaces(r)% baryon_density( itr )
          STOP
        ENDIF
        IF( pos_surfaces(r)% psurface_vol2( itr ) < 0 )THEN
          PRINT *, "When mirroring particles"
          PRINT *, "pos_surfaces(", r, ")% psurface_vol2( ", itr, " ) =", &
                   pos_surfaces(r)% psurface_vol2( itr )
          STOP
        ENDIF

      ENDDO

    ENDDO
    PRINT *, " * Final number of particles=", npart_out
    PRINT *

    ! Safety checks (maybe redundant at this point, but better to be paranoid)
    npart_test= 0
    DO r= 1, n_surfaces, 1
      npart_test= npart_test + npart_shell(r)
    ENDDO
    IF( npart_test /= npart_out )THEN
      PRINT *, "** ERROR! The sum of the particles on the shells is not ", &
               "equal to the total number of particles. Stopping.."
      PRINT *, " * npart_test", npart_test, ", npart_out=", npart_out
      PRINT *
      STOP
    ENDIF
    IF( SUM( npart_shell, DIM=1 ) /= npart_out )THEN
      PRINT *, "** ERROR! The sum of the particles on the shells is not ", &
               "equal to the total number of particles. Stopping.."
      PRINT *, " * SUM( npart_shell )", SUM( npart_shell, DIM=1 ), &
               ", npart_out=", npart_out
      PRINT *
      STOP
    ENDIF

    debug_if: IF( debug )THEN

      mass_test= zero
      mass_test2= zero
      proper_volume_test= zero
      proper_volume= zero
      i_shell= 1
      DO r= 1, n_surfaces, 1
        !DO itr= i_shell, (i_shell - 1) + (npart_shelleq(r)**two)/two, 1
        !  CALL bns_obj% import_id( &
        !           pos( 1, itr ), pos( 2, itr ), pos( 3, itr ), &
        !           g_xx, baryon_density, gamma_euler )
        !
        !  pvol( itr )= m_parts(r)/( baryon_density*g_xx*SQRT(g_xx)*gamma_euler )
        !
        !  proper_volume_test= proper_volume_test + two*pvol( itr )*g_xx*SQRT(g_xx)
        !  mass_test= mass_test + &
        !              two*baryon_density*pvol( itr )*g_xx*SQRT(g_xx)*gamma_euler
        !ENDDO
        !i_shell= i_shell + (npart_shelleq(r)**two)/two
        DO itr= 1, npart_shell(r), 1

          IF( pos_surfaces(r)% baryon_density( itr ) == 0 )THEN
            PRINT *, "When computing particle volume"
            PRINT *, r, itr, pos_surfaces(r)% pos_shell( 1, itr ), &
                     pos_surfaces(r)% pos_shell( 2, itr ), &
                     pos_surfaces(r)% pos_shell( 3, itr ), &
                     pos_surfaces(r)% baryon_density( itr )
          ENDIF
          !IF( pos_surfaces(r)% psurface_vol2( itr ) <= zero )THEN
          !  PRINT *, "When computing particle volume"
          !  PRINT *, "pos_surfaces(", r, ")% psurface_vol2( ", itr, " ) =", &
          !           pos_surfaces(r)% psurface_vol2( itr )
          !  STOP
          !ENDIF

          pos_surfaces(r)% psurface_vol( itr )= m_parts(r) &
                            /( pos_surfaces(r)% baryon_density( itr ) &
                              *pos_surfaces(r)% sqdetg( itr ) &
                              *pos_surfaces(r)% gamma_euler( itr ) )

          proper_volume_test= proper_volume_test + &
                              pos_surfaces(r)% psurface_vol( itr )

          proper_volume= proper_volume + pos_surfaces(r)% psurface_vol2( itr )

          mass_test= mass_test + pos_surfaces(r)% baryon_density( itr ) &
                    *pos_surfaces(r)% psurface_vol( itr ) &
                    *pos_surfaces(r)% sqdetg( itr ) &
                    *pos_surfaces(r)% gamma_euler( itr )
          mass_test2= mass_test2 + pos_surfaces(r)% baryon_density( itr ) &
                    *pos_surfaces(r)% psurface_vol2( itr ) &
                    *pos_surfaces(r)% sqdetg( itr ) &
                    *pos_surfaces(r)% gamma_euler( itr )

        ENDDO
      ENDDO

      PRINT *, mass_test, mass_test2, mass_star, proper_volume_test, proper_volume
      PRINT *

      mass_test = zero
      mass_test2= zero
      proper_volume_test = zero
      proper_volume= zero
      surface_vol(r)  = zero
      surface_vol2(r) = zero
      mass_surface(r) = zero
      mass_surface2(r)= zero
      DO r= 1, n_surfaces, 1
        DO itr= 1, npart_shell(r), 1
          !IF( pos_surfaces(r)% psurface_vol2( itr ) <= 0 )THEN
          !  PRINT *, "When computing shell volumes and masses"
          !  PRINT *, "pos_surfaces(", r, ")% psurface_vol2( ", itr, " ) =", &
          !           pos_surfaces(r)% psurface_vol2( itr )
          !  STOP
          !ENDIF
          surface_vol(r)  = surface_vol(r)  + pos_surfaces(r)% psurface_vol( itr )
          surface_vol2(r) = surface_vol2(r) + pos_surfaces(r)% psurface_vol2( itr )
          mass_surface(r) = mass_surface(r) + &
                            pos_surfaces(r)% baryon_density( itr ) &
                           *pos_surfaces(r)% psurface_vol( itr ) &
                           *pos_surfaces(r)% sqdetg( itr ) &
                           *pos_surfaces(r)% gamma_euler( itr )
          mass_surface2(r)= mass_surface2(r) + &
                            pos_surfaces(r)% baryon_density( itr ) &
                           *pos_surfaces(r)% psurface_vol2( itr ) &
                           *pos_surfaces(r)% sqdetg( itr ) &
                           *pos_surfaces(r)% gamma_euler( itr )
        ENDDO
        mass_test= mass_test + mass_surface(r)
        mass_test2= mass_test2 + mass_surface2(r)
        proper_volume= proper_volume + surface_vol(r)
        proper_volume_test= proper_volume_test + surface_vol2(r)
        IF( r > 1 )THEN
          PRINT *, "shell", r
          PRINT *, "  shell volumes:", surface_vol(r), surface_vol2(r), &
                   4.0D0/3.0D0*pi* &
                   ( surface_radii(r)**3.0D0 - surface_radii( r - 1 )**3.0D0 )
          PRINT *, "  shell masses:", mass_surface(r), mass_surface2(r), &
                   surface_masses(r)
          PRINT *
        ENDIF
      ENDDO
      PRINT *
      PRINT *, "masses of the star:", mass_test, mass_test2, mass_star
      PRINT *, "volumes of the star:", proper_volume, proper_volume_test, &
                                       4.0D0/3.0D0*pi*max_radius**3.0D0
      PRINT *

      !STOP

      !DO r= 1, n_surfaces, 1
      !  DO itr= 1, npart_shell(r), 1
      !
      !    PRINT*, (m_parts(r)*MSun/amu)/pos_surfaces(r)% psurface_vol( itr ) &
      !            /(pos_surfaces(r)% g_xx( itr ) &
      !              *SQRT(pos_surfaces(r)% g_xx( itr )) &
      !              *pos_surfaces(r)% gamma_euler( itr )), &
      !            pos_surfaces(r)% baryon_density( itr )*MSun/amu
      !  ENDDO
      !ENDDO
      !STOP

    ENDIF debug_if

    !---------------------------------------!
    !--  Save particles to output arrays  --!
    !---------------------------------------!

    IF(.NOT.ALLOCATED( pos ))THEN
      ALLOCATE( pos( 3, npart_out ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array pos in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( pvol ))THEN
      ALLOCATE( pvol( npart_out ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array pvol in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( nu ))THEN
      ALLOCATE( nu( npart_out ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array nu in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF

    cnt= 0
    npart_shell_tmp= 0
    DO r= 1, n_surfaces, 1
      DO itr= 1, npart_shell(r), 1
        pos( 1, itr + npart_shell_tmp )= pos_surfaces(r)% pos_shell( 1, itr )
        pos( 2, itr + npart_shell_tmp )= pos_surfaces(r)% pos_shell( 2, itr )
        pos( 3, itr + npart_shell_tmp )= pos_surfaces(r)% pos_shell( 3, itr )
        pvol(   itr + npart_shell_tmp )= pos_surfaces(r)% psurface_vol2( itr )
        nu(     itr + npart_shell_tmp )= m_parts(r)/amu*Msun
        cnt= cnt + 1
      ENDDO
      npart_shell_tmp= cnt
    ENDDO
    ! Safety check
    IF( cnt /= npart_out )THEN
      PRINT *, "** ERROR! The sum of the particles on the shells is not ", &
               "equal to the total number of particles. Stopping.."
      PRINT *, "cnt", cnt, ", npart_out=", npart_out
      PRINT *
      STOP
    ENDIF

    !-------------------------------------------------------------------!
    !--  Print particle positions to file (TODO: make this optional)  --!
    !-------------------------------------------------------------------!

    PRINT *, " * Printing particle positions to file..."
    PRINT *

    IF( PRESENT(filename_shells_pos) )THEN
      finalnamefile= filename_shells_pos
    ELSE
      finalnamefile= surface_type//"_surfaces_pos.dat"
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

    !DO itr = 1, npart_out, 1
    !
    !  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    !    pos( 1, itr ), pos( 2, itr ), pos( 3, itr )
    !
    !  IF( ios > 0 )THEN
    !    PRINT *, "...error when writing the arrays in " &
    !             // TRIM(finalnamefile), ". The error message is", err_msg
    !    STOP
    !  ENDIF
    !
    !ENDDO

    DO r= 1, n_surfaces, 1
      DO itr= 1, npart_shell(r), 1

        WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
          r, pos_surfaces(r)% pos_shell( 1, itr ), &
          pos_surfaces(r)% pos_shell( 2, itr ), &
          pos_surfaces(r)% pos_shell( 3, itr )

        IF( ios > 0 )THEN
          PRINT *, "...error when writing the arrays in " &
                   // TRIM(finalnamefile), ". The error message is", err_msg
          STOP
        ENDIF

      ENDDO
    ENDDO

    CLOSE( UNIT= 2 )

    PRINT *, " * SUBROUTINE place_particles_ellipsoidal_surfaces executed."
    PRINT *

    IF( debug ) PRINT *, "20"



    CONTAINS



    FUNCTION validate_position_final( x, y, z ) RESULT( answer )

      !*******************************************************
      !
      !# Returns validate_position( x, y, z ) if the latter
      !  is present, `.TRUE.` otherwise
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
      !# validate_position( x, y, z ) if the latter is present,
      !  `.TRUE.` otherwise

      IF( PRESENT(validate_position) )THEN

        answer= validate_position( x, y, z )

      ELSE

        answer= .TRUE.

      ENDIF

    END FUNCTION validate_position_final


    SUBROUTINE allocate_surface_memory()

      !**************************************
      !
      !# Allocates memory for the surfaces
      !
      !  FT 21.04.2022
      !
      !**************************************

      IMPLICIT NONE

      IF(.NOT.ALLOCATED( surface_radii ))THEN
        ALLOCATE( surface_radii( n_surfaces ), STAT= ios, &
                  ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array surface_radii in SUBROUTINE" &
                    // "allocate_surface_memory. ", &
                    "The error message is", err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( surface_masses ))THEN
        ALLOCATE( surface_masses( n_surfaces ), STAT= ios, &
                  ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array surface_masses in SUBROUTINE" &
                    // "allocate_surface_memory. ", &
                    "The error message is", err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( shell_scales ))THEN
        ALLOCATE( shell_scales( n_surfaces ), STAT= ios, &
                  ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array shell_scales in SUBROUTINE" &
                    // "allocate_surface_memory. ", &
                    "The error message is", err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( surface_vol ))THEN
        ALLOCATE( surface_vol( n_surfaces ), STAT= ios, &
                  ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array surface_vol in SUBROUTINE" &
                    // "allocate_surface_memory. ", &
                    "The error message is", err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( surface_vol2 ))THEN
        ALLOCATE( surface_vol2( n_surfaces ), STAT= ios, &
                  ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array surface_vol2 in SUBROUTINE" &
                    // "allocate_surface_memory. ", &
                    "The error message is", err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( mass_surface ))THEN
        ALLOCATE( mass_surface( n_surfaces ), STAT= ios, &
                  ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array surface_vol in SUBROUTINE" &
                    // "allocate_surface_memory. ", &
                    "The error message is", err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( mass_surface2 ))THEN
        ALLOCATE( mass_surface2( n_surfaces ), STAT= ios, &
                  ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array surface_vol2 in SUBROUTINE" &
                    // "allocate_surface_memory. ", &
                    "The error message is", err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( m_parts ))THEN
        ALLOCATE( m_parts( 1:n_surfaces ), STAT= ios, &
                  ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array m_parts in SUBROUTINE" &
                    // "allocate_surface_memory. ", &
                    "The error message is", err_msg
           STOP
        ENDIF
      ENDIF

    END SUBROUTINE allocate_surface_memory


    SUBROUTINE initialize_surfaces()

      !**************************************
      !
      !# Initializes the fields before the
      !  iteration over the surfaces
      !
      !  FT 21.04.2022
      !
      !**************************************

      IMPLICIT NONE

      ALLOCATE( npart_shell   ( n_surfaces ) )
      ALLOCATE( npart_shelleq ( n_surfaces ) )
      ALLOCATE( alpha         ( n_surfaces ) )
      ALLOCATE( colatitude_pos( n_surfaces ) )
      ALLOCATE( pos_surfaces  ( n_surfaces ) )

      initialization: DO r= 1, n_surfaces, 1

        IF( ALLOCATED( pos_surfaces(r)% pos_shell ) ) &
          DEALLOCATE( pos_surfaces(r)% pos_shell )

        IF( ALLOCATED( pos_surfaces(r)% psurface_vol ) ) &
          DEALLOCATE( pos_surfaces(r)% psurface_vol )

        IF( ALLOCATED( pos_surfaces(r)% psurface_vol2 ) ) &
          DEALLOCATE( pos_surfaces(r)% psurface_vol2 )

        IF( ALLOCATED( pos_surfaces(r)% sqdetg ) )&
          DEALLOCATE( pos_surfaces(r)% sqdetg )

        IF( ALLOCATED( pos_surfaces(r)% baryon_density ) ) &
          DEALLOCATE( pos_surfaces(r)% baryon_density )

        IF( ALLOCATED( pos_surfaces(r)% gamma_euler ) ) &
          DEALLOCATE( pos_surfaces(r)% gamma_euler )

        IF( ALLOCATED( pos_surfaces(r)% pos_th ) ) &
          DEALLOCATE( pos_surfaces(r)% pos_th )

        IF( ALLOCATED( pos_surfaces(r)% pos_phi ) ) &
          DEALLOCATE( pos_surfaces(r)% pos_phi )

        ALLOCATE( pos_surfaces(r)% pos_shell  ( 3, npart_des ) )
        ALLOCATE( pos_surfaces(r)% psurface_vol  ( npart_des ) )
        ALLOCATE( pos_surfaces(r)% psurface_vol2 ( npart_des ) )
        ALLOCATE( pos_surfaces(r)% sqdetg        ( npart_des ) )
        ALLOCATE( pos_surfaces(r)% baryon_density( npart_des ) )
        ALLOCATE( pos_surfaces(r)% gamma_euler   ( npart_des ) )
        ALLOCATE( pos_surfaces(r)% pos_th        ( npart_des ) )
        ALLOCATE( pos_surfaces(r)% pos_phi       ( npart_des ) )

        pos_surfaces(r)% pos_shell     =   zero
        pos_surfaces(r)% pos_phi       = - one
        pos_surfaces(r)% pos_th        = - one
        pos_surfaces(r)% psurface_vol  =   zero
        pos_surfaces(r)% psurface_vol2 =   zero
        pos_surfaces(r)% sqdetg        =   zero
        pos_surfaces(r)% baryon_density=   zero
        pos_surfaces(r)% gamma_euler   =   zero
        m_parts(r)                     =   m_p
        npart_shelleq(r)= CEILING(SQRT(DBLE(2*surface_masses(r)/m_parts(r))))

      ENDDO initialization

      IF( ALLOCATED(pos) )THEN
        DEALLOCATE(pos)
        ALLOCATE( pos( 3, 2*npart_des ) )
      ENDIF
      IF( ALLOCATED(pvol) )THEN
        DEALLOCATE(pvol)
        ALLOCATE( pvol( 2*npart_des ) )
      ENDIF
      IF( ALLOCATED(nu) )THEN
        DEALLOCATE(nu)
        ALLOCATE( nu( 2*npart_des ) )
      ENDIF

      pos            = zero
      nu             = zero
      phase          = zero
      proper_volume  = zero
      surface_vol    = zero
      surface_vol2   = zero
      dr_shells      = max_radius/n_surfaces
      npart_out      = 0
      upper_bound_tmp= upper_bound
      lower_bound_tmp= lower_bound
      first_r        = CEILING(DBLE(n_surfaces)/two)
      r              = first_r
      cnt2 = 0
      r_cnt= 1

      ! These array are needed to be able to parallelize the loops on each surface
      ALLOCATE( pos_shell_tmp  ( 3, 5*CEILING(SQRT(DBLE(2*npart_des))), &
                                    5*CEILING(SQRT(DBLE(2*npart_des))) ) )
      ALLOCATE( sqdetg_tmp     (    5*CEILING(SQRT(DBLE(2*npart_des))), &
                                    5*CEILING(SQRT(DBLE(2*npart_des))) ) )
      ALLOCATE( bar_density_tmp(    5*CEILING(SQRT(DBLE(2*npart_des))), &
                                    5*CEILING(SQRT(DBLE(2*npart_des))) ) )
      ALLOCATE( gam_euler_tmp  (    5*CEILING(SQRT(DBLE(2*npart_des))), &
                                    5*CEILING(SQRT(DBLE(2*npart_des))) ) )
      ALLOCATE( pvol_tmp       (    5*CEILING(SQRT(DBLE(2*npart_des))), &
                                    5*CEILING(SQRT(DBLE(2*npart_des))) ) )
      ALLOCATE( npart_discarded(    5*CEILING(SQRT(DBLE(2*npart_des))), &
                                    5*CEILING(SQRT(DBLE(2*npart_des))) ) )
      ALLOCATE( npart_surface_tmp(  5*CEILING(SQRT(DBLE(2*npart_des))), &
                                    5*CEILING(SQRT(DBLE(2*npart_des))) ) )

    END SUBROUTINE initialize_surfaces


  END PROCEDURE place_particles_ellipsoidal_surfaces


  FUNCTION number_surfaces( m_p, center, radius, get_dens ) &
           RESULT( n_surfaces )

    !************************************************
    !
    !# Compute the number of surfaces
    !  by integrating the linear particle density
    !  along the larger equatorial radius
    !
    !  FT 22.07.2021
    !
    !************************************************

    USE constants, ONLY: third
    USE utility,   ONLY: zero, is_finite_number

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT( IN ):: m_p, center, radius
    INTERFACE
      FUNCTION get_dens( x, y, z ) RESULT( density )
        !! Returns the baryon mass density at the desired point
        DOUBLE PRECISION, INTENT(IN):: x
        !! \(x\) coordinate of the desired point
        DOUBLE PRECISION, INTENT(IN):: y
        !! \(y\) coordinate of the desired point
        DOUBLE PRECISION, INTENT(IN):: z
        !! \(z\) coordinate of the desired point
        DOUBLE PRECISION:: density
        !> Baryon mass density at \((x,y,z)\)
      END FUNCTION get_dens
    END INTERFACE


    INTEGER:: n_surfaces

    INTEGER:: r
    DOUBLE PRECISION:: n_surfaces_tmp, x_tmp, rho_tmp
  !  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: particle_profile
  !
  !  IF(.NOT.ALLOCATED( particle_profile ))THEN
  !    ALLOCATE( particle_profile( 2, 500 ), STAT= ios, &
  !              ERRMSG= err_msg )
  !    IF( ios > 0 )THEN
  !       PRINT *, "...allocation error for array particle_profile in" &
  !                // "FUNCTION number_surfaces. ", &
  !                "The error message is", err_msg
  !       STOP
  !    ENDIF
  !  ENDIF

    n_surfaces_tmp= zero

    DO r= 1, 500, 1

      x_tmp  = center + r*radius/DBLE(500)
      rho_tmp= get_dens( x_tmp, zero, zero )

      IF( .NOT.is_finite_number(( rho_tmp/m_p )**third) )THEN
        CYCLE
      ENDIF

      n_surfaces_tmp= n_surfaces_tmp + &
                      radius/DBLE(500)*(( rho_tmp/m_p )**third)

    ENDDO

    IF( .NOT.is_finite_number(n_surfaces_tmp) )THEN
      PRINT *, "** ERROR in SUBROUTINE number_surfaces!", &
               " n_surfaces_tmp = ", n_surfaces_tmp, " is not a finite number!"
      PRINT *
      STOP
    ENDIF
    IF( n_surfaces_tmp <= zero )THEN
      PRINT *, "** ERROR in SUBROUTINE number_surfaces!", &
               " n_surfaces_tmp = ", n_surfaces_tmp, " is nonpositive!"
      PRINT *
      STOP
    ENDIF

    n_surfaces= NINT( n_surfaces_tmp )


  END FUNCTION number_surfaces


  SUBROUTINE reallocate_array_1d( array, new_dim )

    !************************************
    !
    !# Reallocate a 1-dimensional array
    !
    !  FT 22.07.2021
    !
    !************************************

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(IN OUT):: array
    INTEGER, INTENT(IN):: new_dim

    DEALLOCATE( array, STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...deallocation error in SUBROUTINE" &
                // "reallocate_tmp_variable. ", &
                "The error message is", err_msg, ", and IOSTAT= ", ios
       STOP
    ENDIF

    ALLOCATE( array( new_dim ), STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...allocation error in SUBROUTINE" &
                // "reallocate_tmp_variable. ", &
                "The error message is", err_msg, ", and IOSTAT= ", ios
       STOP
    ENDIF

  END SUBROUTINE reallocate_array_1d


  SUBROUTINE reallocate_array_2d( array, new_dim, new_dim2 )

    !************************************
    !
    !# Reallocate a 2-dimensional array
    !
    !  FT 22.07.2021
    !
    !************************************

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(IN OUT):: array
    INTEGER, INTENT(IN):: new_dim, new_dim2

    DEALLOCATE( array, STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...deallocation error in SUBROUTINE" &
                // "reallocate_tmp_variable. ", &
                "The error message is", err_msg, ", and IOSTAT= ", ios
       STOP
    ENDIF

    ALLOCATE( array( new_dim, new_dim2 ), STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...allocation error in SUBROUTINE" &
                // "reallocate_tmp_variable. ", &
                "The error message is", err_msg, ", and IOSTAT= ", ios
       STOP
    ENDIF

  END SUBROUTINE reallocate_array_2d


  SUBROUTINE place_surfaces( central_dens, center, radius, m_p, n_surfaces, &
                             surface_radii, last_r, get_dens )

    !************************************************
    !
    !# Place the surfaces, according to
    !  the baryon mass density of the star
    !  along the larger equatorial radius
    !
    !  FT 23.07.2021
    !
    !************************************************

    USE constants,  ONLY: third, pi
    USE utility,    ONLY: zero, two, ten

    IMPLICIT NONE

    INTEGER,          INTENT( IN ):: n_surfaces
    DOUBLE PRECISION, INTENT( IN ):: central_dens, center, radius, m_p, last_r
    DOUBLE PRECISION, DIMENSION( n_surfaces ), INTENT( IN OUT ):: surface_radii
    INTERFACE
      FUNCTION get_dens( x, y, z ) RESULT( density )
        !! Returns the baryon mass density at the desired point
        DOUBLE PRECISION, INTENT(IN):: x
        !! \(x\) coordinate of the desired point
        DOUBLE PRECISION, INTENT(IN):: y
        !! \(y\) coordinate of the desired point
        DOUBLE PRECISION, INTENT(IN):: z
        !! \(z\) coordinate of the desired point
        DOUBLE PRECISION:: density
        !> Baryon mass density at \((x,y,z)\)
      END FUNCTION get_dens
    END INTERFACE

    INTEGER:: i, th, phi
    DOUBLE PRECISION:: rho_tmp, long, lat, rad

    surface_radii= zero

    IF( central_dens > zero )THEN
      surface_radii(1)= ( central_dens/m_p )**(-third)
    ELSE
      DO i= 1, 1000, 1
        IF( get_dens( center + radius*i/(ten*ten*ten), zero, zero ) > zero )THEN
          surface_radii(1)= &
      ( get_dens( center + radius*i/(ten*ten*ten), zero, zero )/m_p )**(-third)
          EXIT
        ENDIF
      ENDDO
    ENDIF

    DO itr= 2, n_surfaces, 1

      rho_tmp= get_dens( center + surface_radii( itr - 1 ), zero, zero )

      IF( rho_tmp <= 1.0D-13 )THEN

        rad_loop: DO i= 1, 1000, 1
          DO th= 0, 10, 1
            DO phi= 0, 20, 1

              long= phi/(two*ten)*(two*pi)
              lat = th/ten*pi
              rad = (center + surface_radii( itr - 1 ) + radius*i/(ten*ten*ten))

              IF( get_dens( rad*SIN(lat)*COS(long),&
                            rad*SIN(lat)*SIN(long), &
                            rad*COS(lat) ) > 1.0D-13 )THEN

                rho_tmp= get_dens( rad*SIN(lat)*COS(long),&
                                   rad*SIN(lat)*SIN(long), &
                                   rad*COS(lat) )

                EXIT rad_loop

              ENDIF

            ENDDO
          ENDDO
        ENDDO rad_loop

      ENDIF

      IF( rho_tmp == zero )THEN
        surface_radii= surface_radii*itr/n_surfaces
      ENDIF

      surface_radii( itr )= surface_radii( itr - 1 ) + ( rho_tmp/m_p )**(-third)

    ENDDO
    surface_radii= surface_radii*(radius*last_r/surface_radii(n_surfaces))

  END SUBROUTINE place_surfaces


  SUBROUTINE assign_surfaces_mass( surface_masses, surface_radii, radius, dr, &
                                   n_surfaces, mass_profile_idx, mass_profile, &
                                   mass_star )

    !*************************************************
    !
    !# Assign a mass to each surface,
    !  based on the radial mass profile of the star
    !  (computed along the larger equatorial radius)
    !
    !  FT 23.07.2021
    !
    !*************************************************

    USE utility,  ONLY: zero

    IMPLICIT NONE

    INTEGER,          INTENT( IN ):: n_surfaces
    DOUBLE PRECISION, INTENT( IN ):: radius, dr, mass_star

    INTEGER, DIMENSION( : ),                 INTENT( IN ):: mass_profile_idx
    DOUBLE PRECISION, DIMENSION( n_surfaces ), INTENT( IN ):: surface_radii
    DOUBLE PRECISION, DIMENSION( :, : ),     INTENT( IN ):: mass_profile
    DOUBLE PRECISION, DIMENSION( n_surfaces ), INTENT( IN OUT ):: surface_masses

    INTEGER shell_index, itr2

    shell_index= 1
    itr2= 1
    surface_masses= zero
    assign_masses_to_surfaces: DO itr= 1, NINT(radius/dr), 1

      IF( shell_index == n_surfaces )THEN

        surface_masses( shell_index )= SUM( mass_profile( 2, &
         mass_profile_idx(itr2):mass_profile_idx(NINT(radius/dr)-1) ), DIM= 1 )

        EXIT

      ENDIF

      IF( mass_profile( 1, mass_profile_idx(itr) ) &
          >= surface_radii( shell_index ) &!+ radius/DBLE(2*n_surfaces)
      )THEN

       surface_masses( shell_index )= SUM( mass_profile( 2, &
                     mass_profile_idx(itr2):mass_profile_idx(itr) ), DIM= 1 )

       itr2= itr + 1
       shell_index= shell_index + 1

      ENDIF

    ENDDO assign_masses_to_surfaces

    ! Safety check
    IF( ABS( SUM( surface_masses, DIM= 1 ) - mass_star )/mass_star > 5.0D-3 )THEN
      PRINT *, " ** The masses of the shells do not add up to the ", &
               "mass of the star. Stopping..."
      PRINT *, " * SUM( surface_masses )= ", SUM( surface_masses, DIM=1 )
      PRINT *, " * Baryon mass of the star= ", mass_star
      PRINT *, " * Array surface_masses=", surface_masses
      PRINT *
      STOP
    ENDIF

  END SUBROUTINE assign_surfaces_mass


  SUBROUTINE print_mass_profile_surface_radii( mass_profile, mass_profile_idx, &
                                               surface_radii, radius, dr, &
                                               n_surfaces, &
                                               filename_mass_profile, &
                                               filename_shells_radii )

    !*************************************************
    !
    !# Print star's radial mass profile and radii of
    !  surfaces to different ASCII files
    !
    !  FT 23.07.2021
    !
    !*************************************************

    !USE constants,  ONLY: third

    IMPLICIT NONE

    INTEGER,          INTENT( IN ):: n_surfaces
    DOUBLE PRECISION, INTENT( IN ):: radius, dr

    INTEGER, DIMENSION( : ),                 INTENT( IN ):: mass_profile_idx
    DOUBLE PRECISION, DIMENSION( n_surfaces ), INTENT( IN ):: surface_radii
    DOUBLE PRECISION, DIMENSION( :, : ),     INTENT( IN ):: mass_profile

    CHARACTER( LEN= * ), INTENT( IN ):: filename_mass_profile, &
                                        filename_shells_radii

    LOGICAL:: exist

    PRINT *, " * Print mass profile to file..."
    PRINT *

    INQUIRE( FILE= TRIM(filename_mass_profile), EXIST= exist )

    IF( exist )THEN
      OPEN( UNIT= 2, FILE= TRIM(filename_mass_profile), STATUS= "REPLACE", &
            FORM= "FORMATTED", &
            POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
            IOMSG= err_msg )
    ELSE
      OPEN( UNIT= 2, FILE= TRIM(filename_mass_profile), STATUS= "NEW", &
            FORM= "FORMATTED", &
            ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening " // TRIM(filename_mass_profile), &
               ". The error message is", err_msg
      STOP
    ENDIF

    write_data_loop: DO itr = 1, NINT(radius/dr) - 1, 1

      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        mass_profile( 1, mass_profile_idx(itr) ), &
        mass_profile( 2, mass_profile_idx(itr) ), &
        mass_profile( 3, mass_profile_idx(itr) )

      IF( ios > 0 )THEN
        PRINT *, "...error when writing the arrays in " &
                 // TRIM(filename_mass_profile), ". The error message is", &
                 err_msg
        STOP
      ENDIF

    ENDDO write_data_loop

    CLOSE( UNIT= 2 )

    PRINT *, " * Print surfaces' radii to file..."
    PRINT *

    INQUIRE( FILE= TRIM(filename_shells_radii), EXIST= exist )

    IF( exist )THEN
      OPEN( UNIT= 2, FILE= TRIM(filename_shells_radii), STATUS= "REPLACE", &
            FORM= "FORMATTED", &
            POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
            IOMSG= err_msg )
    ELSE
      OPEN( UNIT= 2, FILE= TRIM(filename_shells_radii), STATUS= "NEW", &
            FORM= "FORMATTED", &
            ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening " // TRIM(filename_shells_radii), &
              ". The error message is", err_msg
      STOP
    ENDIF

    DO itr = 1, n_surfaces, 1

      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        surface_radii( itr )

      IF( ios > 0 )THEN
        PRINT *, "...error when writing the arrays in " &
                 // TRIM(filename_shells_radii), ". The error message is", &
                 err_msg
        STOP
      ENDIF

    ENDDO

    CLOSE( UNIT= 2 )

  END SUBROUTINE print_mass_profile_surface_radii


  FUNCTION particle_volume( rad, col, dr_shells, dth_shells, dphi_shells, th, &
                            colatitudes, npart_equator ) RESULT( pvol )

    !*******************************************
    !
    !# Compute the geometrical particle volume
    !  not the proper particle volume.
    !
    !  FT 23.07.2021
    !
    !*******************************************

    USE constants,  ONLY: pi
    USE utility,    ONLY: two

    IMPLICIT NONE

    INTEGER,          INTENT( IN ):: th, npart_equator
    DOUBLE PRECISION, INTENT( IN ):: rad, col, dr_shells, dphi_shells
    DOUBLE PRECISION, INTENT( IN OUT ):: dth_shells
    DOUBLE PRECISION, DIMENSION(:), INTENT( IN ):: colatitudes

    DOUBLE PRECISION:: pvol

    IF( th == 1 )THEN

    !dth_shells= pi - ( col + colatitude_pos(r)% colatitudes(th+1) )/two
      IF( npart_equator == 4 )THEN

        dth_shells= pi

      ELSE

        dth_shells= two*ABS( col - &
                ( col + colatitudes(th + 1) )/two )
      ENDIF

    ELSEIF( th == npart_equator/4 )THEN

    !dth_shells= ( colatitude_pos(r)% colatitudes(th-1) + col - pi )/two
      dth_shells= two*ABS( ( colatitudes(th - 1) &
                + col )/two - col )

    ELSE

      dth_shells= ABS( &
              ( colatitudes(th + 1) + col )/two &
            - ( col + colatitudes(th - 1) )/two )

    ENDIF

    pvol= rad**two*SIN(col)*dr_shells*dth_shells*dphi_shells! &

  END FUNCTION particle_volume


  SUBROUTINE compute_colatitudes_uniformly_in( alpha, beta, colatitudes )

    !**************************************************
    !
    !# Compute the colatitudes according to a
    !  uniform distribution over a
    !  surface, between alpha and beta, with
    !  pi/2 < alpha < beta < pi.
    !  The values are stored in the array colatitudes
    !  See https://mathworld.wolfram.com/SpherePointPicking.html
    !
    !  FT 6.10.2021
    !
    !**************************************************

    USE constants,  ONLY: pi
    USE utility,    ONLY: one, two, four

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: alpha, beta
    DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: colatitudes

    INTEGER:: n, size_col, i

    IF( alpha < pi/two .OR. alpha > pi )THEN
      PRINT *, "ERROR in SUBROUTINE compute_colatitudes_uniformly_in!", &
               " Argument alpha should lie in [pi/2,pi]. Stopping..."
      PRINT *
      STOP
    ENDIF
    IF( beta < pi/two .OR. beta > pi )THEN
      PRINT *, "** ERROR in SUBROUTINE compute_colatitudes_uniformly_in!", &
               " Argument beta should lie in [pi/2,pi]. Stopping..."
      PRINT *
      STOP
    ENDIF
    IF( alpha > beta )THEN
      PRINT *, "** ERROR in SUBROUTINE compute_colatitudes_uniformly_in!", &
               " Argument alpha should be less than argument beta. Stopping..."
      PRINT *
      STOP
    ENDIF

    size_col= SIZE(colatitudes)
    n= 4*size_col

    DO i= 1, size_col, 1

      colatitudes(i)= ACOS( two*( DBLE( i + 1 )*( COS(alpha) + one )/two &
            + ( ( COS(beta) + one )/two &
              - (DBLE(n + 1)/four + one)*( COS(alpha) + one )/two ) &
                      *four*DBLE(i)/DBLE(n + 1) ) - one )

      !PRINT *, "colatitudes(", i, ")", colatitudes(i)/pi, "pi"

      IF( .NOT.is_finite_number( colatitudes(i) ) )THEN
        PRINT *, "** ERROR in SUBROUTINE compute_colatitudes_uniformly_in! ", &
                 "colatitudes(", i, ") is a NaN! Stopping.."
        PRINT *, DBLE( i + 1 )*( COS(alpha) + one )/two
        PRINT *, ( COS(beta) + one )/two
        PRINT *, DBLE(size_col) + one
        PRINT *, DBLE(i)/DBLE(size_col)
        PRINT *
        STOP
      ENDIF

    ENDDO

  END SUBROUTINE compute_colatitudes_uniformly_in


END SUBMODULE ellipsoidal_surfaces
