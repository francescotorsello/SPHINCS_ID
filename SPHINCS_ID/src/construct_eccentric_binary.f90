! File:         construct_eccentric_binary.f90
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

PROGRAM construct_eccentric_binary

  !*****************************************************
  !
  !# Read two |tov| star BSSN and SPH ID, and construct
  !  an approximate eccentric binary system with a
  !  velocity determined with Newtonian mechanics and
  !  Newtonian gravitational law
  !
  !  FT 12.12.2022
  !
  !*****************************************************

  !
  !-- SPHINCS_fix_metric MODULES
  !
  USE sph_variables, ONLY: npart, &  ! particle number
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
                           !temp,  &  ! Temperature
                           !av,    &  ! Dissipation
                           ye,    &  ! Electron fraction
                           !divv,  &  ! Divergence of velocity vel_u
                           allocate_SPH_memory, &
                           deallocate_SPH_memory, &
                           n1, n2
  USE input_output,  ONLY: set_units, read_sphincs_dump, write_sphincs_dump
  USE units,         ONLY: m0c2_cu

  !
  !-- BSSN MODULES
  !
  USE mesh_refinement,  ONLY: nlevels, levels, initialize_grid
  USE ADM_refine,       ONLY: allocate_ADM
  USE BSSN_refine,      ONLY: allocate_BSSN, read_BSSN_dump

  !
  !-- SPHINCS_ID MODULES
  !
  USE utility,       ONLY: zero, one, Msun_geo, spacetime_vector_norm_sym4x4
  USE lorentz_group, ONLY: eta


  IMPLICIT NONE


  INTEGER,          PARAMETER:: max_npart    = 5D+7
  DOUBLE PRECISION, PARAMETER:: periastron_km= 45.D0
  INTEGER, PARAMETER :: my_unit = 42

  INTEGER:: a
  DOUBLE PRECISION:: periastron, mass1, mass2, speed
  CHARACTER(LEN=:), ALLOCATABLE:: filename1, filename2

  INTEGER :: io_error, allocation_status
  INTEGER, DIMENSION(3) :: array_shape
  INTEGER :: l

  !--------------!
  !--  SPH ID  --!
  !--------------!

  !
  !-- Allocate SPH memory (reallocated inside read_tov_id)
  !
  npart= max_npart
  CALL set_units('NSM')
  CALL allocate_SPH_memory()

  !
  !-- Read the two TOV SPH ID
  !-- The first star will be displaced to negative x,
  !-- the second star to positive x, depending on the value of the periastron
  !
  filename1 = 'tov-id-files/NSxx.00000'
  filename2 = 'tov-id-files/NSx2.00000'
  periastron= periastron_km/Msun_geo
  CALL read_tov_id(filename1, filename2, periastron)

  !
  !-- Assign Newtonian velocity and generalized Lorentz factor to the particles
  !
  mass1= SUM(nu(1:n1), DIM=1)*m0c2_cu
  mass2= SUM(nu(n1+1:npart), DIM=1)*m0c2_cu
  speed= newtonian_speed(mass1, mass2, periastron)
  IF(speed > one)THEN
    PRINT *, "** ERROR! The Newtonian speed is larger than the speed of light!"
    PRINT *, " * Newtonian speed=", speed
    PRINT *, " * Stopping..."
    STOP
  ENDIF
  PRINT *, " * Newtonian speed=", speed

  ! At periastron, only the \(y\) component of the velocity is nonzero,
  ! assuming that the motion happens in the \(xy\) plane
  !$OMP PARALLEL DO DEFAULT( NONE ) &
  !$OMP             SHARED( npart, n1, vel_u, Theta, speed ) &
  !$OMP             PRIVATE( a )
  compute_vel_and_theta_on_particles: DO a= 1, npart, 1

    vel_u(1,a)= zero
    IF(a <= n1) vel_u(2,a)=   speed
    IF(a >  n1) vel_u(2,a)= - speed
    vel_u(3,a)= zero

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
  filename1= 'eccentric-binary-id-files/NSNS.00000'
  CALL write_sphincs_dump(filename1)

  !
  !-- Deallocate SPH memory
  !
  CALL deallocate_SPH_memory()


  !---------------!
  !--  BSSN ID  --!
  !---------------!

  filename1 = 'tov-id-files/BSSN_vars.00000'
  OPEN(UNIT=my_unit, FILE=filename1, FORM='unformatted', ACTION='read', &
       STATUS='old', IOSTAT=io_error)
  IF (io_error/=0) THEN
    PRINT*,'Error opening ', filename1, ' for reading'
    STOP
  ENDIF

  READ(my_unit) nlevels

  ALLOCATE (levels(nlevels),STAT=allocation_status)
  IF(allocation_status > 0)THEN
     PRINT*,'...allocation error for levels'
     STOP
  ENDIF

  DO l= 1, nlevels

    READ(my_unit) array_shape
    levels(l)% ngrid_x= array_shape(1)
    levels(l)% ngrid_y= array_shape(2)
    levels(l)% ngrid_z= array_shape(3)

  ENDDO

  CLOSE(my_unit)

 ! CALL initialize_grid()

  CALL allocate_ADM()
  CALL allocate_BSSN()

  filename1 = 'tov-id-files/BSSN_vars.00000'
  CALL read_BSSN_dump( 00000, filename1 )

  !
  !-- Print the BSSN ID
  !



  CONTAINS



  SUBROUTINE read_tov_id(filename1, filename2, periastron)

    !***********************************************************
    !
    !# Read the two TOV ID produced with setup_TOV.x, and
    !  place them symmetrically on the \(x\) axis so that
    !  their distance is equal to the periastron given as input
    !
    !  FT 13.12.2022
    !
    !***********************************************************

    IMPLICIT NONE

    DOUBLE PRECISION,              INTENT(IN)   :: periastron
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(INOUT):: filename1, filename2

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos_u1, vel_u1, &
                                                    pos_u2, vel_u2
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: u1, nu1, h1, nlrf1, &
                                                    Pr1, Ye1, Theta1, &
                                                    u2, nu2, h2, nlrf2, &
                                                    Pr2, Ye2, Theta2

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

    CALL deallocate_SPH_memory()

    npart= n1 + n2

    CALL allocate_SPH_memory()

    pos_u(1,1:n1)  = pos_u1(1,:) - periastron/2.D0
    pos_u(2:3,1:n1)= pos_u1(2:3,:)
    vel_u(:,1:n1)  = vel_u1
    u    (1:n1)    = u1
    nu   (1:n1)    = nu1
    h    (1:n1)    = h1
    nlrf (1:n1)    = nlrf1
    Pr   (1:n1)    = Pr1
    Ye   (1:n1)    = Ye1
    Theta(1:n1)    = Theta1

    pos_u(1,n1 + 1: npart)  = pos_u2(1,:) + periastron/2.D0
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

  END SUBROUTINE read_tov_id


  FUNCTION newtonian_speed(mass1, mass2, periastron) RESULT(v)

    !***********************************************************
    !
    !# Compute Newtonian speed at the periastron from the
    !  conservatn of total energy
    !
    !  FT 13.12.2022
    !
    !***********************************************************

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: mass1, mass2, periastron

    DOUBLE PRECISION:: v

    v= SQRT(4*mass1*mass2/(mass1 + mass2))/periastron

  END FUNCTION newtonian_speed


END PROGRAM construct_eccentric_binary
