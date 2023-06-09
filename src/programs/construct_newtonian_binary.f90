! File:         construct_newtonian_binary.f90
! Author:       Francesco Torsello (FT)
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

PROGRAM construct_newtonian_binary

  !*****************************************************
  !
  !# Read two |tov| and |sph| |id|, and construct
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
                            n1,     &  ! particle number for star 1
                            pos_u,  &  ! particle positions
                            vel_u,  &  ! particle velocities in
                                       ! coordinate frame
                            nu,     &  ! canonical baryon number per
                                       ! particle
                            Theta,  &  ! Generalized Lorentz factor
                            deallocate_sph_memory
  USE input_output,   ONLY: write_sphincs_dump
  USE units,          ONLY: m0c2_cu
  USE tensor,         ONLY: n_sym3x3, jx, jy, jz, &
                            jxx, jxy, jxz, jyy, jyz, jzz, &
                            itt, itx, ity, itz, ixx, ixy, ixz, iyy, iyz, izz

  !
  !-- BSSN MODULES
  !
  USE ADM_refine,                 ONLY: deallocate_ADM
  USE BSSN_refine,                ONLY: allocate_BSSN, deallocate_BSSN, &
                                        write_BSSN_dump
  USE Tmunu_refine,               ONLY: deallocate_Tmunu
  USE GravityAcceleration_refine, ONLY: allocate_GravityAcceleration, &
                                        deallocate_GravityAcceleration
  USE mesh_refinement,            ONLY: output_1d, output_2d
  USE McLachlan_refine,           ONLY: allocate_Ztmp, deallocate_Ztmp, &
                                        ADM_to_BSSN, BSSN_Constraints
  USE BSSN_refine,                ONLY: write_BSSN_constraints

  !
  !-- SPHINCS_ID MODULES
  !
  USE utility,        ONLY: zero, one, two, Msun_geo, &
                            spacetime_vector_norm_sym4x4, is_finite_number
  USE lorentz_group,  ONLY: eta, lorentz_boost


  IMPLICIT NONE


  INTEGER:: a
  DOUBLE PRECISION:: periastron, mass1, mass2, radius1, radius2, x1, x2, &
                     energy, angular_momentum, distance, &
                     semimajor_axis, semiminor_axis, apoastron
  DOUBLE PRECISION, DIMENSION(3):: v1, v2
  CHARACTER(LEN=:), ALLOCATABLE:: filename1, filename2

  INTEGER, PARAMETER:: parameters_unit= 17
  INTEGER, PARAMETER:: max_length= 100
  INTEGER:: stat
  DOUBLE PRECISION:: periastron_parameter, distance_km, eccentricity
  CHARACTER(LEN=:), ALLOCATABLE:: parameters_namefile
  CHARACTER(LEN=max_length):: &
    common_path, filename_sph1, filename_sph2, filename_tov1, filename_tov2, &
    output_directory, sph_output_file, bssn_output_file
  CHARACTER(LEN=100):: msg
  LOGICAL:: file_exists

  NAMELIST /newtonian_binary_parameters/ &
            periastron_parameter, distance_km, eccentricity, &
            common_path, filename_sph1, filename_sph2, &
            filename_tov1, filename_tov2, output_directory, &
            sph_output_file, bssn_output_file

  !
  !-- Read parameters
  !
  parameters_namefile= 'newtonian_binary_parameters.dat'

  INQUIRE( FILE= parameters_namefile, EXIST= file_exists )
  IF( file_exists )THEN

    OPEN( parameters_unit, FILE= parameters_namefile, STATUS= 'OLD' )

  ELSE

    PRINT*
    PRINT*,'** ERROR: ', parameters_namefile, " file not found!"
    PRINT*
    STOP

  ENDIF

  READ( UNIT= parameters_unit, NML= newtonian_binary_parameters, IOSTAT= stat, &
        IOMSG= msg )

  IF( stat /= 0 )THEN
    PRINT *, "** ERROR: Error in reading ", parameters_namefile, &
             ". The IOSTAT variable is ", stat, &
             "The error message is", msg
    STOP
  ENDIF

  CLOSE( UNIT= parameters_unit )

  !
  !-- Check that the parameters are reasonable
  !
  IF(eccentricity < zero)THEN

    PRINT *, "** ERROR! The value for the eccentricity in the parameter", &
             " file newtonian_binary_parameters.dat is negative!"
    PRINT *, "   eccentricity= ", eccentricity
    PRINT *, " * Stopping..."
    PRINT *
    STOP

  ENDIF
  IF(periastron_parameter <= zero)THEN

    PRINT *, "** ERROR! The value for the periastron_parameter in the", &
             " parameter file newtonian_binary_parameters.dat is nonpositive!"
    PRINT *, "   periastron_parameter= ", periastron_parameter
    PRINT *, " * Stopping..."
    PRINT *
    STOP

  ENDIF
  IF(distance_km <= zero)THEN

    PRINT *, "** ERROR! The value for the initial distance in the", &
             " parameter file newtonian_binary_parameters.dat is nonpositive!"
    PRINT *, "   distance_km= ", distance_km
    PRINT *, " * Stopping..."
    PRINT *
    STOP

  ENDIF

  ! Convert initial distance to code units
  distance= distance_km/Msun_geo


  !--------------!
  !--  SPH ID  --!
  !--------------!

  !
  !-- Read the two TOV SPH ID
  !-- The first star will be displaced to negative x,
  !-- the second star to positive x, depending on the value of the periastron
  !
  filename1= TRIM(common_path)//TRIM(filename_sph1)
  filename2= TRIM(common_path)//TRIM(filename_sph2)
  CALL read_tov_sph_id(filename1,filename2)

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

    radius1= MAX(radius1, SQRT(pos_u(1,a)**2 + pos_u(2,a)**2 + pos_u(3,a)**2))

  ENDDO find_radius_star1
  !$OMP END PARALLEL DO
  radius2= zero
  !$OMP PARALLEL DO DEFAULT( NONE ) &
  !$OMP             SHARED( npart, n1, pos_u ) &
  !$OMP             PRIVATE( a ) &
  !$OMP             REDUCTION( MAX: radius2 )
  find_radius_star2: DO a= n1 + 1, npart, 1

    radius2= MAX(radius2, SQRT(pos_u(1,a)**2 + pos_u(2,a)**2 + pos_u(3,a)**2))

  ENDDO find_radius_star2
  !$OMP END PARALLEL DO
  PRINT *, " * Radius of star 1=", radius1, "Msun=", radius1*Msun_geo, "km"
  PRINT *, " * Radius of star 2=", radius2, "Msun=", radius2*Msun_geo, "km"
  PRINT *

  !
  !-- Set periastron between the stars
  !
  periastron= (radius1 + radius2)/periastron_parameter
  PRINT *, " * Chosen periastron_parameter=", periastron_parameter
  PRINT *, " * Periastron= (radius1 + radius2)/periastron_parameter =", &
           periastron, "Msun_geo=", periastron*Msun_geo, "km"
  PRINT *

  ! Check that the requested initial distance is equal to, or larger than,
  ! the periastron
  IF(distance < periastron)THEN

    PRINT *
    PRINT *, "** ERROR! The chosen initial distance is strictly smaller than", &
             " the chosen periastron!"
    PRINT *, " * Initial distance= ", distance, "Msun_geo=", distance_km, "km"
    PRINT *, " * Periastron= (radius1 + radius2)/periastron_parameter =", &
             periastron, "Msun_geo=", periastron*Msun_geo, "km"
    PRINT *, " * Stopping..."
    PRINT *
    STOP

  ENDIF
  ! Check that the requested initial distance is equal to, or larger than,
  ! the sum of the two radii
  IF(distance < radius1 + radius2)THEN

    PRINT *
    PRINT *, "** ERROR! The chosen initial distance is strictly smaller than", &
             " the sum of the radii of the stars!"
    PRINT *, " * Initial distance= ", distance, "Msun_geo=", distance_km, "km"
    PRINT *, " * radius1 + radius2 =", &
             radius1 + radius2, "Msun_geo=", (radius1 + radius2)*Msun_geo, "km"
    PRINT *, " * Stopping..."
    PRINT *
    STOP

  ENDIF

  PRINT *, " * Chosen initial distance between the stars=", distance, &
           "Msun_geo=", distance_km, "km"
  PRINT *

  PRINT *, " * Chosen eccentricity=", eccentricity
  IF(eccentricity == zero)THEN
  ! Circle

     PRINT *, " * The orbit is a circle."

  ELSEIF(eccentricity < one)THEN
  ! Ellipse

    ! Compute ellipse parameters
    semimajor_axis= periastron/(one - eccentricity)
    apoastron     = (one + eccentricity)*semimajor_axis
    semiminor_axis= SQRT(periastron*apoastron)

    PRINT *, " * The orbit is an ellipse."
    PRINT *, " * Apoastron= ", apoastron, "Msun_geo=", apoastron*Msun_geo, "km"
    PRINT *, " * Semi-major axis= ", semimajor_axis, "Msun_geo=", &
             semimajor_axis*Msun_geo, "km"
    PRINT *, " * Semi-minor axis= ", semiminor_axis, "Msun_geo=", &
             semiminor_axis*Msun_geo, "km"

    IF(apoastron < radius1 + radius2)THEN

      PRINT *
      PRINT *, "** ERROR! The apoastron is strictly smaller than", &
               " the sum of the radii of the stars!"
      PRINT *, " * Apoastron= ", apoastron, "Msun_geo=", apoastron*Msun_geo,"km"
      PRINT *, " * radius1 + radius2 =", &
               radius1 + radius2, "Msun_geo=", (radius1 + radius2)*Msun_geo,"km"
      PRINT *, " * Stopping..."
      PRINT *
      STOP

    ENDIF

    IF(distance > apoastron)THEN

      PRINT *
      PRINT *, "** ERROR! The chosen initial distance is strictly larger than",&
               " the apoastron!"
      PRINT *, "   Initial distance= ", distance, "Msun_geo=", distance_km, "km"
      PRINT *, " * Apoastron= ", apoastron, "Msun_geo=", apoastron*Msun_geo,"km"
      PRINT *, " * Stopping..."
      PRINT *
      STOP

    ENDIF

  ELSEIF(eccentricity == one)THEN
  ! Parabola (straight line is not considered here; it would have zero
  ! angular momentum)

    PRINT *, " * The orbit is a parabola."

  ELSEIF(eccentricity > one)THEN
  ! Hyperbola

    PRINT *, " * The orbit is a hyperbola."

  ENDIF
  PRINT *

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
  !-- Compute total, Newtonian, energy and angular momentum of the system
  !
  CALL newtonian_energy_angular_momentum &
    (eccentricity, periastron, mass1, mass2, energy, angular_momentum)

  PRINT *, " * Energy of the system=", energy, "Msun"
  PRINT *, " * Angular_momentum of the system=", angular_momentum, "Msun**2"
  PRINT *

  !
  !-- Compute Newtonian velocities and generalized Lorentz factors,
  !-- and assign them to the particles
  !
  CALL newtonian_speeds &
    (mass1, mass2, energy, angular_momentum, distance, v1, v2)
  ! Check that the velocities are acceptable
  IF(.NOT.is_finite_number(NORM2(v1)))THEN
    PRINT *, "** ERROR! The Newtonian speed for star 1 has some NaN ", &
             "components!"
    PRINT *, " * Newtonian speed=", NORM2(v1), "c"
    PRINT *, " * Newtonian velocity=", v1, "c"
    PRINT *, " * Stopping..."
    STOP
  ENDIF
  IF(.NOT.is_finite_number(NORM2(v2)))THEN
    PRINT *, "** ERROR! The Newtonian speed for star 2 has some NaN ", &
             "components!"
    PRINT *, " * Newtonian speed=", NORM2(v2), "c"
    PRINT *, " * Newtonian velocity=", v2, "c"
    PRINT *, " * Stopping..."
    STOP
  ENDIF
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
  filename1= TRIM(output_directory)//TRIM(sph_output_file)
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

  filename1= TRIM(common_path)//TRIM(filename_tov1)
  filename2= TRIM(common_path)//TRIM(filename_tov2)
  CALL read_boost_superimpose_tov_adm_id &
    (filename1,filename2, x1, x2, v1, v2, radius1, radius2)

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
  filename1= TRIM(output_directory)//TRIM(bssn_output_file)
  CALL write_BSSN_dump(filename1)
  PRINT *, "...done."
  PRINT *

  !
  !-- Print the BSSN constraints
  !-- TODO: To do this, the stress-energy tensor is needed
  !
!  CALL BSSN_Constraints
!  CALL output_2d( Tmunu_ll, 3,   dcount, ivar=itt, output_ghosts=.TRUE. )
!  CALL output_2d( Tmunu_ll, 3,   dcount, ivar=itx, output_ghosts=.TRUE. )
!  CALL output_2d( Tmunu_ll, 3,   dcount, ivar=ixx, output_ghosts=.TRUE. )
!  CALL output_2d( lapse, 3,      dcount, output_ghosts=.TRUE. )
!  CALL output_2d( shift_u, 3,    dcount, ivar=jx, output_ghosts=.TRUE. )
!  CALL output_2d( Gamma_u, 3,    dcount, ivar=jx, output_ghosts=.TRUE. )
!  CALL output_2d( phi, 3,        dcount, output_ghosts=.TRUE. )
!  CALL output_2d( trK, 3,        dcount, output_ghosts=.TRUE. )
!  CALL output_2d( A_BSSN3_ll, 3, dcount, ivar=jxx, output_ghosts=.TRUE. )
!  CALL output_2d( g_BSSN3_ll, 3, dcount, ivar=jxx, output_ghosts=.TRUE. )
!  CALL output_2d( Ham, 3,        dcount, output_ghosts=.TRUE. )
!  CALL output_2d( M_l, 3,        dcount, ivar=jx, output_ghosts=.TRUE. )
!  CALL output_2d( M_l, 3,        dcount, ivar=jy, output_ghosts=.TRUE. )
!  CALL output_2d( M_l, 3,        dcount, ivar=jz, output_ghosts=.TRUE. )
!  CALL output_1d( Ham, 1,        dcount, output_ghosts=.TRUE. )
!  CALL output_1d( M_l, 1,        dcount, ivar=jx, output_ghosts=.TRUE. )
!  CALL output_1d( M_l, 1,        dcount, ivar=jy, output_ghosts=.TRUE. )
!  CALL output_1d( M_l, 1,        dcount, ivar=jz, output_ghosts=.TRUE. )

  !
  !-- Deallocate ADM and BSSN memory
  !
  CALL deallocate_ADM()
  CALL deallocate_BSSN()


  PRINT *, "** End of execution. ID files are:"
  PRINT *, "   SPH ID: ", TRIM(output_directory)//TRIM(sph_output_file)
  PRINT *, "   BSSN ID: ", filename1
  PRINT *



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
    USE utility,        ONLY: scan_1d_array_for_nans

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(INOUT):: filename1, filename2

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

    !
    !-- Ensure that the ID does not contain NaNs or infinities
    !
    PRINT *, "** Ensuring that the SPH ID does not have any NaNs or", &
             "infinities..."

    CALL scan_1d_array_for_nans( npart, pos_u(1,:), "pos_u(:,1)" )
    CALL scan_1d_array_for_nans( npart, pos_u(2,:), "pos_u(:,2)" )
    CALL scan_1d_array_for_nans( npart, pos_u(3,:), "pos_u(:,3)" )

    CALL scan_1d_array_for_nans( npart, nlrf, "nlrf" )
    CALL scan_1d_array_for_nans( npart, nu, "nu" )
    CALL scan_1d_array_for_nans( npart, u, "u" )
    CALL scan_1d_array_for_nans( npart, h, "h" )
    CALL scan_1d_array_for_nans( npart, Pr, "Pr" )
    CALL scan_1d_array_for_nans( npart, Ye, "Ye" )
    CALL scan_1d_array_for_nans( npart, Theta, "Theta" )

    CALL scan_1d_array_for_nans( npart, vel_u(1,:), "vel_u(:,1)" )
    CALL scan_1d_array_for_nans( npart, vel_u(2,:), "vel_u(:,2)" )
    CALL scan_1d_array_for_nans( npart, vel_u(3,:), "vel_u(:,3)" )

    PRINT *, " * the SPH ID does not have NaNs or infinities."
    PRINT *

  END SUBROUTINE read_tov_sph_id


  SUBROUTINE read_boost_superimpose_tov_adm_id &
    (filename1, filename2, x1, x2, v1, v2, radius1, radius2)

    !***********************************************************
    !
    !# Read the two BSSN TOV ID files produced with setup_TOV.x,
    !  and place them symmetrically on the \(x\) axis so that
    !  their distance is equal to the periastron given as input
    !
    !  FT 13.12.2022
    !
    !***********************************************************

    USE tensor,          ONLY: jx, jy, jz, jxx, jxy, jxz, &
                               jyy, jyz, jzz, n_sym3x3, n_sym4x4
    USE mesh_refinement, ONLY: nlevels, levels, initialize_grid, &
                               grid_function_scalar, grid_function, &
                               read_grid_params, coords, &
                               allocate_grid_function, deallocate_grid_function
    USE ADM_refine,      ONLY: allocate_ADM, lapse, shift_u, &
                               g_phys3_ll, K_phys3_ll, dt_lapse, dt_shift_u
    USE Tmunu_refine,    ONLY: Tmunu_ll, allocate_Tmunu, deallocate_Tmunu
    USE TOV_refine,      ONLY: read_TOV_dump, allocate_tov, deallocate_tov, &
                               get_tov_metric
    USE utility,         ONLY: zero, compute_tpo_metric, determinant_sym3x3, &
                               scan_3d_array_for_nans, one, two, four

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: x1, x2, radius1, radius2
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: v1, v2
    CHARACTER(LEN=*), INTENT(INOUT):: filename1, filename2

    INTEGER,          PARAMETER:: tov_np= 100001
    DOUBLE PRECISION, PARAMETER:: eps   = 1.75D-1

    DOUBLE PRECISION:: min_abs_z, distance, sigma1, sigma2!= ABS(x1) + ABS(x2)

    INTEGER :: i, j, k, l, unit_att_out, ios

    DOUBLE PRECISION:: tmp, tmp2, tmp3, xtmp, ytmp, ztmp, &
                       g00, g01, g02, g03, g11, g12, g13, g22, g23, g33, &
                       gamma1, gamma2, detg

    DOUBLE PRECISION, DIMENSION(4,4):: g(n_sym4x4)

    CHARACTER(LEN=:), ALLOCATABLE:: attfunc_namefile, err_msg

    LOGICAL:: exist

    TYPE(grid_function_scalar):: lapse1, phi1, trK1, Theta_Z41, lapse_A_BSSN1, &
                                 lapse2, phi2, trK2, Theta_Z42, lapse_A_BSSN2, &
                                 dt_lapse1, dt_lapse2
    TYPE(grid_function_scalar):: attenuating_function1, attenuating_function2
    TYPE(grid_function):: shift_u1, shift_B_BSSN_u1, Gamma_u1, &
                          g_phys3_ll1, g_BSSN3_ll1, A_BSSN3_ll1, &
                          shift_u2, shift_B_BSSN_u2, Gamma_u2, &
                          g_phys3_ll2, g_BSSN3_ll2, A_BSSN3_ll2, &
                          dt_shift_u1, dt_shift_u2, &
                          K_phys3_ll1, K_phys3_ll2, &
                          Tmunu_ll1, Tmunu_ll2

    TYPE(lorentz_boost):: boost1, boost2

    distance= ABS(x1) + ABS(x2)
    !sigma= ABS(x1) + ABS(x2)! - radius1 - radius2

    CALL read_grid_params()
    CALL initialize_grid()

    CALL allocate_tov(tov_np)

    PRINT *
    PRINT *, " * Reading ID for first TOV star..."
    CALL read_tov_dump(filename1)

    !
    !-- Construct boosts and get their Lorentz factors
    !
    boost1= lorentz_boost(v1)
    boost2= lorentz_boost(v2)
    gamma1= boost1% get_lambda()
    gamma2= boost2% get_lambda()

    !sigma1= radius2
    !sigma2= radius1
    sigma1= gamma2*radius2
    sigma2= gamma1*radius1
    !sigma1= ( gamma2*radius2 + gamma2*(distance-radius1) )/two
    !sigma2= ( gamma1*radius1 + gamma1*(distance-radius2) )/two
    !sigma1= gamma2*radius2/((LOG(one/(one - eps)))**(one/four))
    !sigma2= gamma1*radius1/((LOG(one/(one - eps)))**(one/four))

    PRINT *
    PRINT *, "gamma1=", gamma1
    PRINT *, "gamma2=", gamma2
    PRINT *, "sigma1=", sigma1
    PRINT *, "sigma2=", sigma2
    PRINT *

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

    CALL allocate_grid_function(attenuating_function1, 'att_func1')

    read_tov1_id_on_the_mesh: DO l= 1, nlevels, 1
      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( levels, l, coords, lapse1, shift_u1, &
      !$OMP                     g_phys3_ll1, dt_lapse1, dt_shift_u1, &
      !$OMP                     K_phys3_ll1, Tmunu_ll1, x1, x2, boost1, &
      !$OMP                     gamma2, sigma1, attenuating_function1 ) &
      !$OMP             PRIVATE( i, j, k, tmp, tmp2, tmp3, xtmp, ytmp, ztmp, &
      !$OMP                      g00,g01,g02,g03,g11,g12,g13,g22,g23,g33,g )
      DO k= 1, levels(l)% ngrid_z, 1
        DO j= 1, levels(l)% ngrid_y, 1
          DO i= 1, levels(l)% ngrid_x, 1

            xtmp= coords% levels(l)% var(i,j,k,1)! - x1
            ytmp= coords% levels(l)% var(i,j,k,2)
            ztmp= coords% levels(l)% var(i,j,k,3)

            CALL get_tov_metric(xtmp - x1, ytmp, ztmp, &
                                tmp, tmp2, tmp3, &
                                g00,g01,g02,g03,g11,g12,g13,g22,g23,g33)

            g= boost1% &
               apply_as_congruence([g00,g01,g02,g03,g11,g12,g13,g22,g23,g33])

            CALL compute_tpo_metric(g, &
                                    lapse1% levels(l)% var(i,j,k), &
                                    shift_u1% levels(l)% var(i,j,k,:), &
                                    g_phys3_ll1% levels(l)% var(i,j,k,:))

            dt_lapse1%   levels(l)% var(i,j,k)  = zero
            dt_shift_u1% levels(l)% var(i,j,k,:)= zero
            K_phys3_ll1% levels(l)% var(i,j,k,:)= zero
            Tmunu_ll1%   levels(l)% var(i,j,k,:)= zero

            tmp= (gamma2*(xtmp - x2))**2 + ytmp**2 + ztmp**2

            attenuating_function1% levels(l)% var(i,j,k)= &
              one !- EXP( -(tmp**2)/(sigma1**4) )

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

    CALL allocate_grid_function(attenuating_function2, 'att_func2')

    read_tov2_id_on_the_mesh: DO l= 1, nlevels, 1
      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( levels, l, coords, lapse2, shift_u2, &
      !$OMP                     g_phys3_ll2, dt_lapse2, dt_shift_u2, &
      !$OMP                     K_phys3_ll2, Tmunu_ll2, x1, x2, boost2, &
      !$OMP                     gamma1, sigma2, attenuating_function2 ) &
      !$OMP             PRIVATE( i, j, k, tmp, tmp2, tmp3, xtmp, ytmp, ztmp, &
      !$OMP                      g00,g01,g02,g03,g11,g12,g13,g22,g23,g33,g )
      DO k= 1, levels(l)% ngrid_z, 1
        DO j= 1, levels(l)% ngrid_y, 1
          DO i= 1, levels(l)% ngrid_x, 1

            xtmp= coords% levels(l)% var(i,j,k,1)! - x2
            ytmp= coords% levels(l)% var(i,j,k,2)
            ztmp= coords% levels(l)% var(i,j,k,3)

            CALL get_tov_metric(xtmp -x2, ytmp, ztmp, &
                                tmp, tmp2, tmp3, &
                                g00,g01,g02,g03,g11,g12,g13,g22,g23,g33)

            g= boost2% &
               apply_as_congruence([g00,g01,g02,g03,g11,g12,g13,g22,g23,g33])

            CALL compute_tpo_metric(g, &
                                    lapse2% levels(l)% var(i,j,k), &
                                    shift_u2% levels(l)% var(i,j,k,:), &
                                    g_phys3_ll2% levels(l)% var(i,j,k,:))

            dt_lapse2%   levels(l)% var(i,j,k)  = zero
            dt_shift_u2% levels(l)% var(i,j,k,:)= zero
            K_phys3_ll2% levels(l)% var(i,j,k,:)= zero
            Tmunu_ll2%   levels(l)% var(i,j,k,:)= zero

            tmp= (gamma1*(xtmp - x1))**2 + ytmp**2 + ztmp**2

            attenuating_function2% levels(l)% var(i,j,k)= &
              one !- EXP( -(tmp**2)/(sigma2**4) )

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
      !$OMP                     lapse, shift_u, attenuating_function1, &
      !$OMP                     attenuating_function2 ) &
      !$OMP             PRIVATE( i, j, k, tmp, tmp2, tmp3, &
      !$OMP                      g00,g01,g02,g03,g11,g12,g13,g22,g23,g33 )
      DO k= 1, levels(l)% ngrid_z, 1
        DO j= 1, levels(l)% ngrid_y, 1
          DO i= 1, levels(l)% ngrid_x, 1

            g_phys3_ll% levels(l)% var(i,j,k,1)= one + &
              attenuating_function1% levels(l)% var(i,j,k) &
              *(g_phys3_ll1% levels(l)% var(i,j,k,1) - one) + &
              attenuating_function2% levels(l)% var(i,j,k) &
              *(g_phys3_ll2% levels(l)% var(i,j,k,1) - one)

            g_phys3_ll% levels(l)% var(i,j,k,2)= &
             attenuating_function1% levels(l)% var(i,j,k) &
             *g_phys3_ll1% levels(l)% var(i,j,k,2) + &
             attenuating_function2% levels(l)% var(i,j,k) &
             *g_phys3_ll2% levels(l)% var(i,j,k,2)

            g_phys3_ll% levels(l)% var(i,j,k,3)= &
             attenuating_function1% levels(l)% var(i,j,k) &
             *g_phys3_ll1% levels(l)% var(i,j,k,3) + &
             attenuating_function2% levels(l)% var(i,j,k) &
             *g_phys3_ll2% levels(l)% var(i,j,k,3)

            g_phys3_ll% levels(l)% var(i,j,k,4)= one + &
             attenuating_function1% levels(l)% var(i,j,k) &
             *(g_phys3_ll1% levels(l)% var(i,j,k,4) - one) + &
             attenuating_function2% levels(l)% var(i,j,k) &
             *(g_phys3_ll2% levels(l)% var(i,j,k,4) - one)

            g_phys3_ll% levels(l)% var(i,j,k,5)= &
             attenuating_function1% levels(l)% var(i,j,k) &
             *g_phys3_ll1% levels(l)% var(i,j,k,5) + &
             attenuating_function2% levels(l)% var(i,j,k) &
             *g_phys3_ll2% levels(l)% var(i,j,k,5)

            g_phys3_ll% levels(l)% var(i,j,k,6)= one + &
             attenuating_function1% levels(l)% var(i,j,k) &
             *(g_phys3_ll1% levels(l)% var(i,j,k,6) - one) + &
             attenuating_function2% levels(l)% var(i,j,k) &
             *(g_phys3_ll2% levels(l)% var(i,j,k,6) - one)

            lapse%     levels(l)% var(i,j,k)= - one + &
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

    attfunc_namefile= "attenutating_functions.dat"

    PRINT *, "** Printing attenuating functions to file ", &
             TRIM(attfunc_namefile), "..."

    INQUIRE( FILE= TRIM(attfunc_namefile), EXIST= exist )

    IF( exist )THEN
        OPEN( UNIT= unit_att_out, &
              FILE= TRIM(attfunc_namefile), &
              STATUS= "REPLACE", FORM= "FORMATTED", &
              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
              IOMSG= err_msg )
    ELSE
        OPEN( UNIT= unit_att_out, &
              FILE= TRIM(attfunc_namefile), &
              STATUS= "NEW", FORM= "FORMATTED", &
              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening " // &
               TRIM(attfunc_namefile), &
               ". The error message is", err_msg
      STOP
    ENDIF

    DO l= 1, nlevels, 1

      min_abs_z= HUGE(one)
      DO k= 1, levels(l)% ngrid_z, 1
        IF( ABS( coords% levels(l)% var(1,1,k,3) ) < ABS( min_abs_z ) )THEN
          min_abs_z= coords% levels(l)% var(1,1,k,3)
        ENDIF
      ENDDO

      DO k= 1, levels(l)% ngrid_z, 1
        DO j= 1, levels(l)% ngrid_y, 1
          DO i= 1, levels(l)% ngrid_x, 1

             IF( ABS(coords% levels(l)% var(i,j,k,3) - min_abs_z) &
                 /ABS(min_abs_z) > 1.D-5 &
             )THEN

               CYCLE

             ENDIF

             WRITE( UNIT = unit_att_out, IOSTAT = ios, IOMSG = err_msg, &
                    FMT = * ) &
               l, &
               coords% levels(l)% var(i,j,k,1), &
               coords% levels(l)% var(i,j,k,2), &
               coords% levels(l)% var(i,j,k,3), &
               attenuating_function1% levels(l)% var(i,j,k), &
               attenuating_function2% levels(l)% var(i,j,k)
          ENDDO
        ENDDO
      ENDDO

    ENDDO

    CLOSE( UNIT= unit_att_out )

    PRINT *, " * attenuating functions printed to file ", &
             TRIM(attfunc_namefile)
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

    !
    !-- Deallocate attenuating functions
    !
    CALL deallocate_grid_function(attenuating_function1, 'att_func1')
    CALL deallocate_grid_function(attenuating_function2, 'att_func2')

    !
    !-- Ensure that the standard 3+1 ID does not contain NaNs,
    !-- and that the determinant of the spatial metric is
    !-- strictly positive
    !
    PRINT *, "** Ensuring that the ID does not have any NaNs or infinities, ", &
             "and that the determinant of the spatial metric is strictly ", &
             "positive..."

    DO l= 1, nlevels, 1

      ASSOCIATE( nx     => levels(l)% ngrid_x, &
                 ny     => levels(l)% ngrid_y, &
                 nz     => levels(l)% ngrid_z, &
                 lapse  => lapse%      levels(l)% var, &
                 shift  => shift_u%    levels(l)% var, &
                 g      => g_phys3_ll% levels(l)% var, &
                 eK     => K_phys3_ll% levels(l)% var )

      CALL scan_3d_array_for_nans( nx, ny, nz, lapse, "lapse" )

      CALL scan_3d_array_for_nans( nx, ny, nz, shift(:,:,:,jx), &
                                   "shift(:,:,:,jx)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, shift(:,:,:,jy), &
                                   "shift(:,:,:,jy)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, shift(:,:,:,jz), &
                                   "shift(:,:,:,jz)" )

      CALL scan_3d_array_for_nans( nx, ny, nz, g(:,:,:,jxx), &
                                   "g_phys3_ll(:,:,:,jxx)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, g(:,:,:,jxy), &
                                   "g_phys3_ll(:,:,:,jxy)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, g(:,:,:,jxz), &
                                   "g_phys3_ll(:,:,:,jxz)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, g(:,:,:,jyy), &
                                   "g_phys3_ll(:,:,:,jyy)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, g(:,:,:,jyz), &
                                   "g_phys3_ll(:,:,:,jyz)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, g(:,:,:,jzz), &
                                   "g_phys3_ll(:,:,:,jzz)" )

      CALL scan_3d_array_for_nans( nx, ny, nz, eK(:,:,:,jxx), &
                                   "K_phys3_ll(:,:,:,jxx)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, eK(:,:,:,jxy), &
                                   "K_phys3_ll(:,:,:,jxy)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, eK(:,:,:,jxz), &
                                   "K_phys3_ll(:,:,:,jxz)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, eK(:,:,:,jyy), &
                                   "K_phys3_ll(:,:,:,jyy)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, eK(:,:,:,jyz), &
                                   "K_phys3_ll(:,:,:,jyz)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, eK(:,:,:,jzz), &
                                   "K_phys3_ll(:,:,:,jzz)" )

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( l, levels, g_phys3_ll, coords ) &
      !$OMP             PRIVATE( i, j, k, detg )
      DO k= 1, levels(l)% ngrid_z, 1
        DO j= 1, levels(l)% ngrid_y, 1
          DO i= 1, levels(l)% ngrid_x, 1

            CALL determinant_sym3x3(g_phys3_ll% levels(l)% var(i,j,k,:), detg)

            IF( detg < 1.D-10 )THEN

              PRINT *, "** ERROR! The " &
                       // "determinant of the spatial metric is " &
                       // "effectively 0 at the grid point " &
                       // "(i,j,k)= (", i, ",", j,",",k, "), " &
                       // "(x,y,z)= ", "(", &
                       coords% levels(l)% var(i, j, k, 1), ",", &
                       coords% levels(l)% var(i, j, k, 2), ",", &
                       coords% levels(l)% var(i, j, k, 3), ")."
              PRINT *
              PRINT *, "   nx, ny, nz =", &
                levels(l)% ngrid_x, levels(l)% ngrid_y, levels(l)% ngrid_z
              PRINT *
              PRINT *, "   detg=", detg
              PRINT *
              PRINT *, "   g_xx=", g_phys3_ll% levels(l)% var(i,j,k,jxx)
              PRINT *, "   g_xy=", g_phys3_ll% levels(l)% var(i,j,k,jxy)
              PRINT *, "   g_xz=", g_phys3_ll% levels(l)% var(i,j,k,jxz)
              PRINT *, "   g_yy=", g_phys3_ll% levels(l)% var(i,j,k,jyy)
              PRINT *, "   g_yz=", g_phys3_ll% levels(l)% var(i,j,k,jyz)
              PRINT *, "   g_zz=", g_phys3_ll% levels(l)% var(i,j,k,jzz)
              PRINT *
              STOP

            ELSEIF( detg < zero )THEN

              PRINT *, "** ERROR! The " &
                       // "determinant of the spatial metric is " &
                       // "negative at the grid point " &
                       // "(i,j,k)= (", i, ",", j,",",k, "), " &
                       // "(x,y,z)= ", "(", &
                       coords% levels(l)% var(i, j, k, 1), ",", &
                       coords% levels(l)% var(i, j, k, 2), ",", &
                       coords% levels(l)% var(i, j, k, 3), ")."
              PRINT *
              PRINT *, "   nx, ny, nz =", &
                levels(l)% ngrid_x, levels(l)% ngrid_y, levels(l)% ngrid_z
              PRINT *
              PRINT *, "   detg=", detg
              PRINT *
              PRINT *, "   g_xx=", g_phys3_ll% levels(l)% var(i,j,k,jxx)
              PRINT *, "   g_xy=", g_phys3_ll% levels(l)% var(i,j,k,jxy)
              PRINT *, "   g_xz=", g_phys3_ll% levels(l)% var(i,j,k,jxz)
              PRINT *, "   g_yy=", g_phys3_ll% levels(l)% var(i,j,k,jyy)
              PRINT *, "   g_yz=", g_phys3_ll% levels(l)% var(i,j,k,jyz)
              PRINT *, "   g_zz=", g_phys3_ll% levels(l)% var(i,j,k,jzz)
              PRINT *
              STOP

            ENDIF

          ENDDO
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

      END ASSOCIATE

    ENDDO

    PRINT *, " * the standard 3+1 ID does not contain NaNs or infinites, ", &
             "and the determinant of the spatial metric is strictly positive."
    PRINT *

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

    !PRINT *, "radial_speed_fictitious=", radial_speed_fictitious

    total_speed_fictitious  = SQRT(two/mu*(energy + mass1*mass2/distance))

    !PRINT *, "total_speed_fictitious=", total_speed_fictitious

    angular_speed_fictitious= SQRT((total_speed_fictitious**2 &
                                    - radial_speed_fictitious**2)/distance**2)

    !PRINT *, "angular_speed_fictitious=", angular_speed_fictitious

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


  SUBROUTINE newtonian_energy_angular_momentum &
    (eccentricity, periastron, mass1, mass2, energy, angular_momentum)

    !***********************************************************
    !
    !# Compute the Newtonian energy and angular momentum of the system,
    !  imposing that the radial velocity of the fictitious
    !  body is 0 at the desired periastron, with the desired eccentricity.
    !
    !  The formulas used here are found by solving the equations
    !  that can be found in:
    !
    !  Goldstein, Poole, Safko, "Classical mechanics",
    !  Sec.3.2, eq.(3.16) with \(\dot(r)=0\), and Sec.3.7, eq.(3.57)
    !  See Landau, Lifshitz, "Mechanics", Chapter III,
    !  eq.(14.5) with \(\dot(r)=0\)
    !
    !  for the energy and the angular momentum.
    !
    !  FT 16.12.2022
    !
    !***********************************************************

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: eccentricity, periastron, mass1, mass2

    DOUBLE PRECISION, INTENT(OUT):: energy, angular_momentum

    DOUBLE PRECISION:: mu

    mu= mass1*mass2/(mass1 + mass2)

    IF(eccentricity == zero)THEN
    ! Circle

      angular_momentum= SQRT(mu*mass1*mass2*periastron)

    ELSEIF(eccentricity == one)THEN
    ! Parabola (straight line is not considered here; it would have zero
    ! angular momentum)

      angular_momentum= SQRT(two*mu*mass1*mass2*periastron)

    ELSEIF(eccentricity > one)THEN
    ! Hyperbola

      angular_momentum= SQRT((one + eccentricity)*mu*mass1*mass2*periastron)

    ELSEIF(zero < eccentricity .AND. eccentricity < one)THEN
    ! Ellipse [SQRT((one - eccentricity)*mu*mass1*mass2*periastron) would be
    ! for an ellispe having apoastron equal to our value of the periastron]

      angular_momentum= SQRT((one + eccentricity)*mu*mass1*mass2*periastron)

    ELSE

      PRINT *, "** ERROR in SUBROUTINE newtonian_energy_angular_momentum!"
      PRINT *, " * The value for the eccentricity is negative!"
      PRINT *, "   eccentricity= ", eccentricity
      PRINT *, " * Stopping..."
      PRINT *
      STOP

    ENDIF

    energy= &
      mu*(mass1*mass2)**2*(eccentricity**2 - one)/(two*angular_momentum**2)

  END SUBROUTINE newtonian_energy_angular_momentum


END PROGRAM construct_newtonian_binary
