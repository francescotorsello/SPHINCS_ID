! File:         construct_newtonian_binary.f90
! Author:       Francesco Torsello (FT)
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

PROGRAM construct_newtonian_binary

  !*****************************************************
  !
  !# Read two |tov| star TOV and SPH ID, and construct
  !  an Newtonian binary system based on the
  !  Newtonian 2-body problem.
  !
  !  See Goldstein, Poole, Safko, "Classical mechanics",
  !  Chapter 3, and Landau, Lifshitz, "Mechanics",
  !  Chapter III
  !
  !  FT 12.12.2022
  !
  !*****************************************************

  !
  !-- SPHINCS_fix_metric MODULES
  !
  USE sph_variables,  ONLY: npart,  &  ! particle number
                            n1, n2, &  ! particle numbers for each star
                            pos_u,  &  ! particle positions
                            vel_u,  &  ! particle velocities in
                                       ! coordinate frame
                            nu,     &  ! canonical baryon number per
                                       ! particle
                            Theta,  &  ! Generalized Lorentz factor
                            deallocate_sph_memory
  USE input_output,   ONLY: write_sphincs_dump
  USE units,          ONLY: m0c2_cu
  USE tensor,         ONLY: n_sym3x3

  !
  !-- BSSN MODULES
  !
  USE ADM_refine,                 ONLY: deallocate_ADM
  USE BSSN_refine,                ONLY: allocate_BSSN, deallocate_BSSN, &
                                        write_BSSN_dump
  USE Tmunu_refine,               ONLY: deallocate_Tmunu
  USE GravityAcceleration_refine, ONLY: allocate_GravityAcceleration, &
                                        deallocate_GravityAcceleration
  USE McLachlan_refine,           ONLY: allocate_Ztmp, deallocate_Ztmp, &
                                        ADM_to_BSSN

  !
  !-- SPHINCS_ID MODULES
  !
  USE utility,        ONLY: zero, one, two, Msun_geo, &
                            spacetime_vector_norm_sym4x4
  USE lorentz_group,  ONLY: eta, lorentz_boost


  IMPLICIT NONE


  DOUBLE PRECISION, PARAMETER:: periastron_parameter= 1.D0
  !DOUBLE PRECISION, PARAMETER:: periastron_km= 10.D0
  DOUBLE PRECISION, PARAMETER:: distance_km  = 50000.D0
  DOUBLE PRECISION, PARAMETER:: energy       = zero

  INTEGER:: a
  DOUBLE PRECISION:: periastron, mass1, mass2, radius1, radius2, x1, x2, &
                     angular_momentum, distance
  DOUBLE PRECISION, DIMENSION(3):: v1, v2
  CHARACTER(LEN=:), ALLOCATABLE:: filename1, filename2

  ! Convert periastron and initial distance to code units
  !periastron= periastron_km/Msun_geo
  distance  = distance_km/Msun_geo

  !PRINT *, " * Chosen periastron=", periastron_km, "km=", periastron, "Msun_geo"
  PRINT *, " * Chosen distance between the stars=", distance_km, "km=", &
           distance, "Msun_geo"
  PRINT *


  !--------------!
  !--  SPH ID  --!
  !--------------!

  !
  !-- Read the two TOV SPH ID
  !-- The first star will be displaced to negative x,
  !-- the second star to positive x, depending on the value of the periastron
  !
  filename1 = 'tov-id-files/NSxx.00000'
  filename2 = 'tov-id-files/NSx2.00000'
  CALL read_tov_sph_id(filename1, filename2)

  !
  !-- Find the radii of the stars, as the maximum radial coordinate
  !-- of a particle
  !
  radius1= zero
  !$OMP PARALLEL DO DEFAULT( NONE ) &
  !$OMP             SHARED( n1, pos_u ) &
  !$OMP             PRIVATE( a ) &
  !$OMP             REDUCTION( MAX: radius1 )
  find_radius_star1: DO a= 1, n1, 1

    radius1= MAX( radius1, SQRT(pos_u(1,a)**2 + pos_u(2,a)**2 + pos_u(3,a)**2) )

  ENDDO find_radius_star1
  !$OMP END PARALLEL DO
  radius2= zero
  !$OMP PARALLEL DO DEFAULT( NONE ) &
  !$OMP             SHARED( npart, n1, pos_u ) &
  !$OMP             PRIVATE( a ) &
  !$OMP             REDUCTION( MAX: radius2 )
  find_radius_star2: DO a= n1 + 1, npart, 1

    radius2= MAX( radius2, SQRT(pos_u(1,a)**2 + pos_u(2,a)**2 + pos_u(3,a)**2) )

  ENDDO find_radius_star2
  !$OMP END PARALLEL DO
  PRINT *, " * Radius of star 1=", radius1, "Msun=", radius1*Msun_geo, "km"
  PRINT *, " * Radius of star 2=", radius2, "Msun=", radius2*Msun_geo, "km"
  PRINT *

  !
  !-- Set periastron between the stars
  !
  periastron= periastron_parameter*(radius1 + radius2)
  PRINT *, " * Chosen periastron_parameter=", periastron_parameter
  PRINT *, " * Periastron = periastron_parameter*(radius1 + radius2) =", &
           periastron, "Msun_geo", periastron*Msun_geo, "km="

  !
  !-- Compute masses of the stars
  !
  mass1= SUM(nu(1:n1), DIM=1)*m0c2_cu
  mass2= SUM(nu(n1+1:npart), DIM=1)*m0c2_cu
  PRINT *, " * Mass of star 1=", mass1, "Msun"
  PRINT *, " * Mass of star 2=", mass2, "Msun"
  PRINT *

  !
  !-- Translate the stars from the origin, along the x axis, so that the
  !-- center of mass of the system is at the origin
  !
  x1= - mass2*distance/(mass1 + mass2)
  x2=   mass1*distance/(mass1 + mass2)
  PRINT *, " * x coordinate of the center of mass of star 1=", x1, "Msun"
  PRINT *, " * x coordinate of the center of mass of star 2=", x2, "Msun"
  PRINT *, " * x coordinate of the center of mass of the system=", &
           (mass1*x1 + mass2*x2)/(mass1 + mass2), "Msun"
  PRINT *

  pos_u(1,1:n1)         = pos_u(1,1:n1) + x1
  pos_u(1,n1 + 1: npart)= pos_u(1,n1 + 1: npart) + x2

  !
  !-- Compute total, Newtonian angular momentum of the system
  !
  angular_momentum= newtonian_angular_momentum(energy,periastron)
  PRINT *, " * Angular_momentum of the system=", angular_momentum, "Msun**2"
  PRINT *

  !
  !-- Assign Newtonian velocity and generalized Lorentz factor to the particles
  !
  CALL newtonian_speeds &
    (mass1, mass2, energy, angular_momentum, distance, v1, v2)
  IF(NORM2(v1) > one)THEN
    PRINT *, "** ERROR! The Newtonian speed for star 1 is larger than the ", &
             "speed of light!"
    PRINT *, " * Newtonian speed=", NORM2(v1), "c"
    PRINT *, " * Newtonian velocity=", v1, "c"
    PRINT *, " * Stopping..."
    STOP
  ENDIF
  IF(NORM2(v2) > one)THEN
    PRINT *, "** ERROR! The Newtonian speed for star 2 is larger than the ", &
             "speed of light!"
    PRINT *, " * Newtonian speed=", NORM2(v2), "c"
    PRINT *, " * Newtonian velocity=", v2, "c"
    PRINT *, " * Stopping..."
    STOP
  ENDIF
  PRINT *, " * Newtonian velocity for star 1=", v1, "c"
  PRINT *, " * Newtonian velocity for star 2=", v2, "c"
  PRINT *
  PRINT *, " * Newtonian speed for star 1=", NORM2(v1), "c"
  PRINT *, " * Newtonian speed for star 2=", NORM2(v2), "c"
  PRINT *

  !$OMP PARALLEL DO DEFAULT( NONE ) &
  !$OMP             SHARED( npart, n1, vel_u, Theta, v1, v2 ) &
  !$OMP             PRIVATE( a )
  compute_vel_and_theta_on_particles: DO a= 1, npart, 1

    IF(a <= n1) vel_u(:,a)= v1
    IF(a >  n1) vel_u(:,a)= v2

    CALL spacetime_vector_norm_sym4x4( eta, &
                                       [one,vel_u(1,a),vel_u(2,a),vel_u(3,a)], &
                                       Theta(a) )
    IF( Theta(a) > zero )THEN
      PRINT *, "** ERROR! The computing frame particle 4-velocity is ", &
               "spacelike at particle ", a
      PRINT *, " * Its norm is ", Theta(a)
      PRINT *, " * Stopping.."
      PRINT *
      STOP
    ELSEIF( Theta(a) == zero )THEN
      PRINT *, "** ERROR! The computing frame particle 4-velocity is ", &
               "null at particle ", a
      PRINT *, " * Its norm is ", Theta(a)
      PRINT *, " * Stopping.."
      PRINT *
      STOP
    ENDIF
    Theta(a)= one/SQRT(-Theta(a))

  ENDDO compute_vel_and_theta_on_particles
  !$OMP END PARALLEL DO

  !
  !-- Print the SPH ID
  !
  PRINT *, " * Printing SPH ID to file..."
  filename1= 'eccentric-binary-id-files/NSNS.00000'
  CALL write_sphincs_dump(filename1)
  PRINT *, "...done."
  PRINT *

  !
  !-- Deallocate SPH memory
  !
  CALL deallocate_sph_memory()


  !---------------!
  !--  BSSN ID  --!
  !---------------!

  filename1 = 'tov-id-files/TOV.00000'
  filename2 = 'tov-id-files/TO2.00000'
  CALL read_boost_superimpose_tov_adm_id(filename1, filename2, x1, x2)

  !
  !-- Compute BSSN ID
  !
  CALL allocate_BSSN()
  CALL allocate_Ztmp()
  CALL allocate_GravityAcceleration()

  CALL ADM_to_BSSN()

  CALL deallocate_Ztmp()
  CALL deallocate_Tmunu()
  CALL deallocate_GravityAcceleration()

  !
  !-- Print the BSSN ID
  !
  PRINT *, " * Printing BSSN ID to file..."
  filename1= 'eccentric-binary-id-files/BSSN_vars.00000'
  CALL write_BSSN_dump(filename1)
  PRINT *, "...done."
  PRINT *

  !
  !-- Deallocate ADM and BSSN memory
  !
  CALL deallocate_ADM()
  CALL deallocate_BSSN()



  CONTAINS



  SUBROUTINE read_tov_sph_id(filename1, filename2)

    !***********************************************************
    !
    !# Read the two SPH TOV ID files produced with setup_TOV.x,
    !  and place them symmetrically on the \(x\) axis so that
    !  their distance is equal to the periastron given as input
    !
    !  FT 13.12.2022
    !
    !***********************************************************

    USE sph_variables,  ONLY: npart,  &  ! particle number
                              n1, n2, &  ! particle numbers for each star
                              pos_u,  &  ! particle positions
                              vel_u,  &  ! particle velocities in
                                         ! coordinate frame
                              nlrf,   &  ! baryon number density in
                                         ! local rest frame
                              nu,     &  ! canonical baryon number per
                                         ! particle
                              Theta,  &  ! Generalized Lorentz factor
                              h,      &  ! Smoothing length
                              Pr,     &  ! Pressure
                              u,      &  ! Internal energy in local rest
                                         ! frame (no kinetic energy)
                              Ye,     &  ! Electron fraction
                              allocate_sph_memory, deallocate_sph_memory
    USE input_output,   ONLY: set_units, read_sphincs_dump

    IMPLICIT NONE

    CHARACTER(LEN=:), ALLOCATABLE, INTENT(INOUT):: filename1, filename2

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos_u1, vel_u1, &
                                                    pos_u2, vel_u2
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: u1, nu1, h1, nlrf1, &
                                                    Pr1, Ye1, Theta1, &
                                                    u2, nu2, h2, nlrf2, &
                                                    Pr2, Ye2, Theta2

    CALL set_units('NSM')

    !
    !-- Read just the particle number, to be able to allocate needed memory
    !
    OPEN(10, file= filename1, form='UNFORMATTED')
    READ(10) npart
    CLOSE(10)

    CALL allocate_sph_memory()

    CALL read_sphincs_dump(filename1)

    !PRINT *, "npart=", npart

    ALLOCATE(pos_u1(3,npart))
    ALLOCATE(vel_u1(3,npart))
    ALLOCATE(u1    (npart))
    ALLOCATE(nu1   (npart))
    ALLOCATE(h1    (npart))
    ALLOCATE(nlrf1 (npart))
    ALLOCATE(Pr1   (npart))
    ALLOCATE(Ye1   (npart))
    ALLOCATE(Theta1(npart))

    n1    = npart
    pos_u1= pos_u(:,1:npart)
    vel_u1= vel_u(:,1:npart)
    u1    = u    (1:npart)
    nu1   = nu   (1:npart)
    h1    = h    (1:npart)
    nlrf1 = nlrf (1:npart)
    Pr1   = Pr   (1:npart)
    Ye1   = Ye   (1:npart)
    Theta1= Theta(1:npart)

    !PRINT *, "SIZE(nlrf)=", SIZE(nlrf1)

    CALL deallocate_sph_memory()

    !
    !-- Read just the particle number, to be able to allocate needed memory
    !
    OPEN(10, file= filename2, form='UNFORMATTED')
    READ(10) npart
    CLOSE(10)

    CALL allocate_sph_memory()

    CALL read_sphincs_dump(filename2)

    ALLOCATE(pos_u2(3,npart))
    ALLOCATE(vel_u2(3,npart))
    ALLOCATE(u2    (npart))
    ALLOCATE(nu2   (npart))
    ALLOCATE(h2    (npart))
    ALLOCATE(nlrf2 (npart))
    ALLOCATE(Pr2   (npart))
    ALLOCATE(Ye2   (npart))
    ALLOCATE(Theta2(npart))

    n2    = npart
    pos_u2= pos_u(:,1:npart)
    vel_u2= vel_u(:,1:npart)
    u2    = u    (1:npart)
    nu2   = nu   (1:npart)
    h2    = h    (1:npart)
    nlrf2 = nlrf (1:npart)
    Pr2   = Pr   (1:npart)
    Ye2   = Ye   (1:npart)
    Theta2= Theta(1:npart)

    CALL deallocate_sph_memory()

    npart= n1 + n2

    CALL allocate_sph_memory()

    pos_u(1,1:n1)  = pos_u1(1,:)
    pos_u(2:3,1:n1)= pos_u1(2:3,:)
    vel_u(:,1:n1)  = vel_u1
    u    (1:n1)    = u1
    nu   (1:n1)    = nu1
    h    (1:n1)    = h1
    nlrf (1:n1)    = nlrf1
    Pr   (1:n1)    = Pr1
    Ye   (1:n1)    = Ye1
    Theta(1:n1)    = Theta1

    pos_u(1,n1 + 1: npart)  = pos_u2(1,:)
    pos_u(2:3,n1 + 1: npart)= pos_u2(2:3,:)
    vel_u(:,n1 + 1: npart)  = vel_u2
    u    (n1 + 1: npart)    = u2
    nu   (n1 + 1: npart)    = nu2
    h    (n1 + 1: npart)    = h2
    nlrf (n1 + 1: npart)    = nlrf2
    Pr   (n1 + 1: npart)    = Pr2
    Ye   (n1 + 1: npart)    = Ye2
    Theta(n1 + 1: npart)    = Theta2

    DEALLOCATE(pos_u1)
    DEALLOCATE(vel_u1)
    DEALLOCATE(u1)
    DEALLOCATE(nu1)
    DEALLOCATE(h1)
    DEALLOCATE(nlrf1)
    DEALLOCATE(Pr1)
    DEALLOCATE(Ye1)
    DEALLOCATE(Theta1)
    DEALLOCATE(pos_u2)
    DEALLOCATE(vel_u2)
    DEALLOCATE(u2)
    DEALLOCATE(nu2)
    DEALLOCATE(h2)
    DEALLOCATE(nlrf2)
    DEALLOCATE(Pr2)
    DEALLOCATE(Ye2)
    DEALLOCATE(Theta2)

  END SUBROUTINE read_tov_sph_id


  SUBROUTINE read_boost_superimpose_tov_adm_id(filename1, filename2, x1, x2)

    !***********************************************************
    !
    !# Read the two BSSN TOV ID files produced with setup_TOV.x,
    !  and place them symmetrically on the \(x\) axis so that
    !  their distance is equal to the periastron given as input
    !
    !  FT 13.12.2022
    !
    !***********************************************************

    USE tensor,          ONLY: n_sym4x4
    USE mesh_refinement, ONLY: nlevels, levels, initialize_grid, &
                               grid_function_scalar, grid_function, &
                               read_grid_params, coords, &
                               allocate_grid_function, deallocate_grid_function
    USE ADM_refine,      ONLY: allocate_ADM, lapse, shift_u, &
                               g_phys3_ll, K_phys3_ll, dt_lapse, dt_shift_u
    USE Tmunu_refine,    ONLY: Tmunu_ll, allocate_Tmunu, deallocate_Tmunu
    USE BSSN_refine,     ONLY: phi, trK, Theta_Z4, lapse_A_BSSN, &
                               shift_B_BSSN_u, Gamma_u, g_BSSN3_ll, A_BSSN3_ll
    USE TOV_refine,      ONLY: read_TOV_dump, allocate_tov, deallocate_tov, &
                               get_tov_metric
    USE utility,         ONLY: compute_tpo_metric

    IMPLICIT NONE

    DOUBLE PRECISION,              INTENT(IN)   :: x1, x2
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(INOUT):: filename1, filename2

    INTEGER, PARAMETER:: tov_np= 100001
    INTEGER :: io_error, allocation_status, i, j, k, l
    INTEGER, DIMENSION(3) :: array_shape

    DOUBLE PRECISION:: tmp, tmp2, tmp3, tmp4, &
                       g00, g01, g02, g03, g11, g12, g13, g22, g23, g33

    DOUBLE PRECISION, DIMENSION(4,4):: g(n_sym4x4)

    TYPE(grid_function_scalar):: lapse1, phi1, trK1, Theta_Z41, lapse_A_BSSN1, &
                                 lapse2, phi2, trK2, Theta_Z42, lapse_A_BSSN2, &
                                 dt_lapse1, dt_lapse2
    TYPE(grid_function):: shift_u1, shift_B_BSSN_u1, Gamma_u1, &
                          g_phys3_ll1, g_BSSN3_ll1, A_BSSN3_ll1, &
                          shift_u2, shift_B_BSSN_u2, Gamma_u2, &
                          g_phys3_ll2, g_BSSN3_ll2, A_BSSN3_ll2, &
                          dt_shift_u1, dt_shift_u2, &
                          K_phys3_ll1, K_phys3_ll2, &
                          Tmunu_ll1, Tmunu_ll2

    TYPE(lorentz_boost):: boost1, boost2

    CALL read_grid_params()
    CALL initialize_grid()

    CALL allocate_tov(tov_np)

    PRINT *
    PRINT *, " * Reading ID for first TOV star..."
    CALL read_tov_dump(filename1)

    CALL allocate_grid_function(lapse1,          'lapse1')
    CALL allocate_grid_function(shift_u1,        'shift_u1', 3)
    CALL allocate_grid_function(g_phys3_ll1,     'g_phys3_ll1', n_sym3x3)
    CALL allocate_grid_function(K_phys3_ll1,     'K_phys3_ll1', n_sym3x3)
    CALL allocate_grid_function(Tmunu_ll1,       'Tmunu_ll1', n_sym3x3)
    CALL allocate_grid_function(dt_lapse1,       'dt_lapse1', n_sym3x3)
    CALL allocate_grid_function(dt_shift_u1,     'dt_shift_u1', n_sym3x3)

    CALL allocate_grid_function(Gamma_u1,        'Gamma_u1', 3)
    CALL allocate_grid_function(phi1,            'phi1')
    CALL allocate_grid_function(trK1,            'trK1')
    CALL allocate_grid_function(A_BSSN3_ll1,     'A_BSSN3_ll1', n_sym3x3)
    CALL allocate_grid_function(g_BSSN3_ll1,     'g_BSSN3_ll1', n_sym3x3)
    CALL allocate_grid_function(lapse_A_BSSN1,   'lapse_A_BSSN1')
    CALL allocate_grid_function(shift_B_BSSN_u1, 'shift_B_BSSN_u1', 3)
    CALL allocate_grid_function(Theta_Z41,       'Theta_Z41')

    boost1= lorentz_boost(v1)

    read_tov1_id_on_the_mesh: DO l= 1, nlevels, 1
      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( levels, l, coords, lapse1, shift_u1, &
      !$OMP                     g_phys3_ll1, dt_lapse1, dt_shift_u1, &
      !$OMP                     K_phys3_ll1, Tmunu_ll1, x1, boost1 ) &
      !$OMP             PRIVATE( i, j, k, tmp, tmp2, tmp3, &
      !$OMP                      g00,g01,g02,g03,g11,g12,g13,g22,g23,g33,g )
      DO k= 1, levels(l)% ngrid_z, 1
        DO j= 1, levels(l)% ngrid_y, 1
          DO i= 1, levels(l)% ngrid_x, 1

            CALL get_tov_metric(coords% levels(l)% var(i,j,k,1) - x1, &
                                coords% levels(l)% var(i,j,k,2), &
                                coords% levels(l)% var(i,j,k,3), &
                                tmp, tmp2, tmp3, &
                                g00,g01,g02,g03,g11,g12,g13,g22,g23,g33 )

            g= boost1% &
               apply_as_congruence([g00,g01,g02,g03,g11,g12,g13,g22,g23,g33])

            CALL compute_tpo_metric( g, &
                                     lapse1% levels(l)% var(i,j,k), &
                                     shift_u1% levels(l)% var(i,j,k,:), &
                                     g_phys3_ll1% levels(l)% var(i,j,k,:) )

            dt_lapse1%   levels(l)% var(i,j,k)  = zero
            dt_shift_u1% levels(l)% var(i,j,k,:)= zero
            K_phys3_ll1% levels(l)% var(i,j,k,:)= zero
            Tmunu_ll1%   levels(l)% var(i,j,k,:)= zero

          ENDDO
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
    ENDDO read_tov1_id_on_the_mesh
    PRINT *, "...done"

    PRINT *
    PRINT *, " * Reading ID for second TOV star..."
    CALL read_tov_dump(filename2)

    CALL allocate_grid_function(lapse2,          'lapse2')
    CALL allocate_grid_function(shift_u2,        'shift_u2', 3)
    CALL allocate_grid_function(g_phys3_ll2,     'g_phys3_ll2', n_sym3x3)
    CALL allocate_grid_function(K_phys3_ll2,     'K_phys3_ll2', n_sym3x3)
    CALL allocate_grid_function(Tmunu_ll2,       'Tmunu_ll2', n_sym3x3)
    CALL allocate_grid_function(dt_lapse2,       'dt_lapse2', n_sym3x3)
    CALL allocate_grid_function(dt_shift_u2,     'dt_shift_u2', n_sym3x3)

    CALL allocate_grid_function(Gamma_u2,        'Gamma_u2', 3)
    CALL allocate_grid_function(phi2,            'phi2')
    CALL allocate_grid_function(trK2,            'trK2')
    CALL allocate_grid_function(A_BSSN3_ll2,     'A_BSSN3_ll2', n_sym3x3)
    CALL allocate_grid_function(g_BSSN3_ll2,     'g_BSSN3_ll2', n_sym3x3)
    CALL allocate_grid_function(lapse_A_BSSN2,   'lapse_A_BSSN2')
    CALL allocate_grid_function(shift_B_BSSN_u2, 'shift_B_BSSN_u2', 3)
    CALL allocate_grid_function(Theta_Z42,       'Theta_Z42')

    boost2= lorentz_boost(v2)

    read_tov2_id_on_the_mesh: DO l= 1, nlevels, 1
      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( levels, l, coords, lapse2, shift_u2, &
      !$OMP                     g_phys3_ll2, dt_lapse2, dt_shift_u2, &
      !$OMP                     K_phys3_ll2, Tmunu_ll2, x2, boost2 ) &
      !$OMP             PRIVATE( i, j, k, tmp, tmp2, tmp3, &
      !$OMP                      g00,g01,g02,g03,g11,g12,g13,g22,g23,g33,g )
      DO k= 1, levels(l)% ngrid_z, 1
        DO j= 1, levels(l)% ngrid_y, 1
          DO i= 1, levels(l)% ngrid_x, 1

            CALL get_tov_metric(coords% levels(l)% var(i,j,k,1) - x2, &
                                coords% levels(l)% var(i,j,k,2), &
                                coords% levels(l)% var(i,j,k,3), &
                                tmp, tmp2, tmp3, &
                                g00,g01,g02,g03,g11,g12,g13,g22,g23,g33 )

            g= boost2% &
               apply_as_congruence([g00,g01,g02,g03,g11,g12,g13,g22,g23,g33])

            CALL compute_tpo_metric( g, &
                                     lapse2% levels(l)% var(i,j,k), &
                                     shift_u2% levels(l)% var(i,j,k,:), &
                                     g_phys3_ll2% levels(l)% var(i,j,k,:) )

            dt_lapse2%   levels(l)% var(i,j,k)  = zero
            dt_shift_u2% levels(l)% var(i,j,k,:)= zero
            K_phys3_ll2% levels(l)% var(i,j,k,:)= zero
            Tmunu_ll2%   levels(l)% var(i,j,k,:)= zero

          ENDDO
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
    ENDDO read_tov2_id_on_the_mesh
    PRINT *, "...done"

    CALL allocate_ADM()
    CALL allocate_Tmunu()

    !
    !-- Sum the translated and boosted TOV ID
    !
    PRINT *
    PRINT *, " * Summing the two TOV ID..."
    sum_tov_id: DO l= 1, nlevels, 1
      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( levels, l, coords, lapse1, shift_u1, &
      !$OMP                     g_phys3_ll1, dt_lapse1, dt_shift_u1, &
      !$OMP                     K_phys3_ll1, Tmunu_ll1, lapse2, shift_u2, &
      !$OMP                     g_phys3_ll2, dt_lapse2, dt_shift_u2, &
      !$OMP                     K_phys3_ll2, Tmunu_ll2, g_phys3_ll, &
      !$OMP                     K_phys3_ll, dt_lapse, dt_shift_u, Tmunu_ll, &
      !$OMP                     lapse, shift_u ) &
      !$OMP             PRIVATE( i, j, k, tmp, tmp2, tmp3, &
      !$OMP                      g00,g01,g02,g03,g11,g12,g13,g22,g23,g33 )
      DO k= 1, levels(l)% ngrid_z, 1
        DO j= 1, levels(l)% ngrid_y, 1
          DO i= 1, levels(l)% ngrid_x, 1

            g_phys3_ll% levels(l)% var(i,j,k,1)= one + &
              (g_phys3_ll1% levels(l)% var(i,j,k,1) - one) + &
              (g_phys3_ll2% levels(l)% var(i,j,k,1) - one)

            g_phys3_ll% levels(l)% var(i,j,k,2)= &
              g_phys3_ll1% levels(l)% var(i,j,k,2) + &
              g_phys3_ll2% levels(l)% var(i,j,k,2)

            g_phys3_ll% levels(l)% var(i,j,k,3)= &
              g_phys3_ll1% levels(l)% var(i,j,k,3) + &
              g_phys3_ll2% levels(l)% var(i,j,k,3)

            g_phys3_ll% levels(l)% var(i,j,k,4)= one + &
              (g_phys3_ll1% levels(l)% var(i,j,k,4) - one) + &
              (g_phys3_ll2% levels(l)% var(i,j,k,4) - one)

            g_phys3_ll% levels(l)% var(i,j,k,5)= &
              g_phys3_ll1% levels(l)% var(i,j,k,5) + &
              g_phys3_ll2% levels(l)% var(i,j,k,5)

            g_phys3_ll% levels(l)% var(i,j,k,6)= one + &
              (g_phys3_ll1% levels(l)% var(i,j,k,6) - one) + &
              (g_phys3_ll2% levels(l)% var(i,j,k,6) - one)

            lapse%      levels(l)% var(i,j,k)  = - one + &
              (lapse1% levels(l)% var(i,j,k) + one) + &
              (lapse2% levels(l)% var(i,j,k) + one)

            shift_u%    levels(l)% var(i,j,k,:)= &
              shift_u1% levels(l)% var(i,j,k,:) + &
              shift_u2% levels(l)% var(i,j,k,:)

            dt_lapse%   levels(l)% var(i,j,k)  = zero
            dt_shift_u% levels(l)% var(i,j,k,:)= zero
            K_phys3_ll% levels(l)% var(i,j,k,:)= zero
            Tmunu_ll%   levels(l)% var(i,j,k,:)= zero

          ENDDO
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
    ENDDO sum_tov_id
    PRINT *, "...done"
    PRINT *


    !
    !-- Deallocate temporary memory
    !
    CALL deallocate_grid_function(lapse1,          'lapse1')
    CALL deallocate_grid_function(shift_u1,        'shift_u1')
    CALL deallocate_grid_function(g_phys3_ll1,     'g_phys3_ll1')
    CALL deallocate_grid_function(K_phys3_ll1,     'K_phys3_ll1')
    CALL deallocate_grid_function(Tmunu_ll1,       'Tmunu_ll1')
    CALL deallocate_grid_function(dt_lapse1,       'dt_lapse1')
    CALL deallocate_grid_function(dt_shift_u1,     'dt_shift_u1')

    CALL deallocate_grid_function(Gamma_u1,        'Gamma_u1')
    CALL deallocate_grid_function(phi1,            'phi1')
    CALL deallocate_grid_function(trK1,            'trK1')
    CALL deallocate_grid_function(A_BSSN3_ll1,     'A_BSSN3_ll1')
    CALL deallocate_grid_function(g_BSSN3_ll1,     'g_BSSN3_ll1')
    CALL deallocate_grid_function(lapse_A_BSSN1,   'lapse_A_BSSN1')
    CALL deallocate_grid_function(shift_B_BSSN_u1, 'shift_B_BSSN_u1')
    CALL deallocate_grid_function(Theta_Z41,       'Theta_Z41')

    CALL deallocate_grid_function(lapse2,          'lapse2')
    CALL deallocate_grid_function(shift_u2,        'shift_u2')
    CALL deallocate_grid_function(g_phys3_ll2,     'g_phys3_ll2')
    CALL deallocate_grid_function(K_phys3_ll2,     'K_phys3_ll2')
    CALL deallocate_grid_function(Tmunu_ll2,       'Tmunu_ll2')
    CALL deallocate_grid_function(dt_lapse2,       'dt_lapse2')
    CALL deallocate_grid_function(dt_shift_u2,     'dt_shift_u2')

    CALL deallocate_grid_function(Gamma_u2,        'Gamma_u2')
    CALL deallocate_grid_function(phi2,            'phi2')
    CALL deallocate_grid_function(trK2,            'trK2')
    CALL deallocate_grid_function(A_BSSN3_ll2,     'A_BSSN3_ll2')
    CALL deallocate_grid_function(g_BSSN3_ll2,     'g_BSSN3_ll2')
    CALL deallocate_grid_function(lapse_A_BSSN2,   'lapse_A_BSSN2')
    CALL deallocate_grid_function(shift_B_BSSN_u2, 'shift_B_BSSN_u2')
    CALL deallocate_grid_function(Theta_Z42,       'Theta_Z42')

  END SUBROUTINE read_boost_superimpose_tov_adm_id


  PURE SUBROUTINE newtonian_speeds &
    (mass1, mass2, energy, angular_momentum, distance, v1, v2)

    !***********************************************************
    !
    !# Compute Newtonian speeds for two stars on parabolic orbits
    !  at a given distance, applying conservation of energy
    !  and momentum.
    !  For a 2-body problem with parabolic orbit, the total
    !  energy is zero.
    !
    !  See Goldstein, Poole, Safko, "Classical mechanics",
    !  Sec.3.2, eq. (3.16); Sec.3.3, eq.(3.21); Sec.3.7
    !  See Landau, Lifshitz, "Mechanics", Chapter III
    !
    !  FT 13.12.2022
    !
    !***********************************************************

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: mass1, mass2, energy, angular_momentum, &
                                   distance
    DOUBLE PRECISION, DIMENSION(3), INTENT(OUT):: v1, v2

    DOUBLE PRECISION:: mu, radial_speed_fictitious, total_speed_fictitious, &
                       angular_speed_fictitious

    DOUBLE PRECISION, DIMENSION(3):: v_fictitious

    mu= mass1*mass2/(mass1 + mass2)

    radial_speed_fictitious = SQRT(two/mu*(energy + mass1*mass2/distance &
                                   - angular_momentum**2/(two*mu*distance**2)))

    total_speed_fictitious  = SQRT(two/mu*(energy + mass1*mass2/distance))

    angular_speed_fictitious= SQRT((total_speed_fictitious**2 &
                                    - radial_speed_fictitious**2)/distance**2)

    v_fictitious(1)= radial_speed_fictitious
    v_fictitious(2)= distance*angular_speed_fictitious
    v_fictitious(3)= zero

    v1(1)= mass2*v_fictitious(1)/(mass1 + mass2)
    v1(2)= mass2*v_fictitious(2)/(mass1 + mass2)
    v1(3)= zero

    v2(1)= - mass1*v_fictitious(1)/(mass1 + mass2)
    v2(2)= - mass1*v_fictitious(2)/(mass1 + mass2)
    v2(3)=   zero

  END SUBROUTINE newtonian_speeds


  PURE FUNCTION newtonian_angular_momentum(energy, periastron) &
    RESULT(angular_momentum)

    !***********************************************************
    !
    !# Compute the classical angular momentum of the system,
    !  imposing that the radial velocity of the fictitious
    !  body is 0 at the desired periastron.
    !
    !  See Goldstein, Poole, Safko, "Classical mechanics",
    !  Sec.3.2, eq.(3.16) with \(\dot(r)=0\)
    !  See Landau, Lifshitz, "Mechanics", Chapter III,
    !  eq.(14.5) with \(\dot(r)=0\)
    !
    !  FT 16.12.2022
    !
    !***********************************************************

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: energy, periastron

    DOUBLE PRECISION:: angular_momentum

    DOUBLE PRECISION:: mu

    mu= mass1*mass2/(mass1 + mass2)

    angular_momentum= &
      SQRT(two*mu*periastron**2*(energy + mass1*mass2/periastron))

  END FUNCTION newtonian_angular_momentum


END PROGRAM construct_newtonian_binary
