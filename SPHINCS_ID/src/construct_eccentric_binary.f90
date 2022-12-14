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
  USE tensor,        ONLY: n_sym3x3

  !
  !-- BSSN MODULES
  !


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
  DOUBLE PRECISION:: periastron, mass1, mass2, v1, v2
  CHARACTER(LEN=:), ALLOCATABLE:: filename1, filename2


  ! Convert periastron to code units
  periastron= periastron_km/Msun_geo


  !--------------!
  !--  SPH ID  --!
  !--------------!

  !
  !-- Allocate SPH memory (reallocated inside read_tov_sph_id)
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
  CALL read_tov_sph_id(filename1, filename2, periastron)

  !
  !-- Assign Newtonian velocity and generalized Lorentz factor to the particles
  !
  mass1= SUM(nu(1:n1), DIM=1)*m0c2_cu
  mass2= SUM(nu(n1+1:npart), DIM=1)*m0c2_cu
  !CALL newtonian_speeds(mass1, mass2, periastron, v1, v2)
  v1= zero
  v2= zero
  IF(v1 > one)THEN
    PRINT *, "** ERROR! The Newtonian speed for star 1 is larger than the ", &
             "speed of light!"
    PRINT *, " * Newtonian speed=", v1
    PRINT *, " * Stopping..."
    STOP
  ENDIF
  IF(v2 > one)THEN
    PRINT *, "** ERROR! The Newtonian speed for star 2 is larger than the ", &
             "speed of light!"
    PRINT *, " * Newtonian speed=", v2
    PRINT *, " * Stopping..."
    STOP
  ENDIF
  PRINT *, " * Newtonian speed for star 1=", v1
  PRINT *, " * Newtonian speed for star 2=", v2
  PRINT *

  ! At periastron, only the \(y\) component of the velocity is nonzero,
  ! assuming that the motion happens in the \(xy\) plane
  !$OMP PARALLEL DO DEFAULT( NONE ) &
  !$OMP             SHARED( npart, n1, vel_u, Theta, v1, v2 ) &
  !$OMP             PRIVATE( a )
  compute_vel_and_theta_on_particles: DO a= 1, npart, 1

    vel_u(1,a)= zero
    IF(a <= n1) vel_u(2,a)=   v1
    IF(a >  n1) vel_u(2,a)= - v2
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

  filename1 = 'tov-id-files/TOV.00000'
  filename2 = 'tov-id-files/TO2.00000'
  CALL read_tov_bssn_id(filename1, filename2, periastron)

  !
  !-- Print the BSSN ID
  !

  !
  !-- Deallocate ADM and BSSN memory
  !
  !CALL deallocate_BSSN()
  !CALL deallocate_ADM()



  CONTAINS



  SUBROUTINE read_tov_sph_id(filename1, filename2, periastron)

    !***********************************************************
    !
    !# Read the two SPH TOV ID files produced with setup_TOV.x,
    !  and place them symmetrically on the \(x\) axis so that
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

  END SUBROUTINE read_tov_sph_id


  SUBROUTINE read_tov_bssn_id(filename1, filename2, periastron)

    !***********************************************************
    !
    !# Read the two BSSN TOV ID files produced with setup_TOV.x,
    !  and place them symmetrically on the \(x\) axis so that
    !  their distance is equal to the periastron given as input
    !
    !  FT 13.12.2022
    !
    !***********************************************************

    USE mesh_refinement, ONLY: nlevels, levels, initialize_grid, &
                               grid_function_scalar, grid_function, &
                               read_grid_params, coords, &
                               allocate_grid_function, deallocate_grid_function
    USE ADM_refine,      ONLY: allocate_ADM, deallocate_ADM, lapse, shift_u
    USE BSSN_refine,     ONLY: allocate_BSSN, deallocate_BSSN, write_BSSN_dump, &
                               phi, trK, Theta_Z4, lapse_A_BSSN, shift_B_BSSN_u,&
                               Gamma_u, g_BSSN3_ll, A_BSSN3_ll
    USE TOV_refine,      ONLY: read_TOV_dump, allocate_tov, deallocate_tov, &
                               get_tov_metric

    USE utility,         ONLY: compute_tpo_metric

    IMPLICIT NONE

    DOUBLE PRECISION,              INTENT(IN)   :: periastron
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(INOUT):: filename1, filename2

    INTEGER, PARAMETER:: tov_np= 100001
    INTEGER :: io_error, allocation_status, i, j, k, l
    INTEGER, DIMENSION(3) :: array_shape

    DOUBLE PRECISION:: tmp, tmp2, tmp3, tmp4, &
                       g00, g01, g02, g03, g11, g12, g13, g22, g23, g33

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

  !  OPEN(UNIT=my_unit, FILE=filename1, FORM='unformatted', ACTION='read', &
  !       STATUS='old', IOSTAT=io_error)
  !  IF (io_error/=0) THEN
  !    PRINT*,'Error opening ', filename1, ' for reading'
  !    STOP
  !  ENDIF
  !
  !  READ(my_unit) nlevels
  !
  !  ALLOCATE (levels(nlevels),STAT=allocation_status)
  !  IF(allocation_status > 0)THEN
  !     PRINT*,'...allocation error for levels'
  !     STOP
  !  ENDIF
  !
  !  DO l= 1, nlevels, 1
  !
  !    READ(my_unit) array_shape
  !    levels(l)% ngrid_x= array_shape(1)
  !    levels(l)% ngrid_y= array_shape(2)
  !    levels(l)% ngrid_z= array_shape(3)
  !
  !  ENDDO
  !
  !  CLOSE(my_unit)
  !
  !  ! CALL initialize_grid()
  !
  !  CALL allocate_ADM()
  !  CALL allocate_BSSN()

    CALL read_grid_params()
    CALL initialize_grid()

    CALL allocate_tov(tov_np)

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

    read_tov1_id_on_the_mesh: DO l= 1, nlevels, 1
      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( levels, l, coords, lapse1, shift_u1, &
      !$OMP                     g_phys3_ll1, dt_lapse1, dt_shift_u1, &
      !$OMP                     K_phys3_ll1, Tmunu_ll1 ) &
      !$OMP             PRIVATE( i, j, k, tmp, tmp2, tmp3, &
      !$OMP                      g00,g01,g02,g03,g11,g12,g13,g22,g23,g33 )
      DO k= 1, levels(l)% ngrid_z, 1
        DO j= 1, levels(l)% ngrid_y, 1
          DO i= 1, levels(l)% ngrid_x, 1

            CALL get_tov_metric(coords% levels(l)% var(i,j,k,1), &
                                coords% levels(l)% var(i,j,k,2), &
                                coords% levels(l)% var(i,j,k,3), &
                                tmp, tmp2, tmp3, &
                                g00,g01,g02,g03,g11,g12,g13,g22,g23,g33 )

            CALL compute_tpo_metric &
              ( [g00, g01, g02, g03, g11, g12, g13, g22, g23, g33] , &
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

  !  DO l= 1, nlevels, 1
  !
  !    lapse1% levels(l)% var         = lapse% levels(l)% var
  !    shift_u1% levels(l)% var       = shift_u% levels(l)% var
  !    Gamma_u1% levels(l)% var       = Gamma_u% levels(l)% var
  !    phi1% levels(l)% var           = phi% levels(l)% var
  !    trK1% levels(l)% var           = trK% levels(l)% var
  !    A_BSSN3_ll1% levels(l)% var    = A_BSSN3_ll% levels(l)% var
  !    g_BSSN3_ll1% levels(l)% var    = g_BSSN3_ll% levels(l)% var
  !    lapse_A_BSSN1% levels(l)% var  = lapse_A_BSSN% levels(l)% var
  !    shift_B_BSSN_u1% levels(l)% var= shift_B_BSSN_u% levels(l)% var
  !    Theta_Z41% levels(l)% var      = Theta_Z4% levels(l)% var
  !
  !  ENDDO
  !
  !  CALL deallocate_BSSN()
  !  CALL deallocate_ADM()
  !  DEALLOCATE(levels)

  !  OPEN(UNIT=my_unit, FILE=filename2, FORM='unformatted', ACTION='read', &
  !       STATUS='old', IOSTAT=io_error)
  !  IF (io_error/=0) THEN
  !    PRINT*,'Error opening ', filename2, ' for reading'
  !    STOP
  !  ENDIF
  !
  !  READ(my_unit) nlevels
  !
  !  ALLOCATE (levels(nlevels),STAT=allocation_status)
  !  IF(allocation_status > 0)THEN
  !     PRINT*,'...allocation error for levels'
  !     STOP
  !  ENDIF
  !
  !  DO l= 1, nlevels, 1
  !
  !    READ(my_unit) array_shape
  !    levels(l)% ngrid_x= array_shape(1)
  !    levels(l)% ngrid_y= array_shape(2)
  !    levels(l)% ngrid_z= array_shape(3)
  !
  !  ENDDO
  !
  !  CLOSE(my_unit)
  !
  !  ! CALL initialize_grid()
  !
  !  CALL allocate_ADM()
  !  CALL allocate_BSSN()

    !filename2 = 'tov-id-files/BSSN_var2.00000'
    !CALL read_BSSN_dump( 00000, filename2 )

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

    read_tov2_id_on_the_mesh: DO l= 1, nlevels, 1
      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( levels, l, coords, lapse2, shift_u2, &
      !$OMP                     g_phys3_ll2, dt_lapse2, dt_shift_u2, &
      !$OMP                     K_phys3_ll2, Tmunu_ll2 ) &
      !$OMP             PRIVATE( i, j, k, tmp, tmp2, tmp3, &
      !$OMP                      g00,g01,g02,g03,g11,g12,g13,g22,g23,g33 )
      DO k= 1, levels(l)% ngrid_z, 1
        DO j= 1, levels(l)% ngrid_y, 1
          DO i= 1, levels(l)% ngrid_x, 1

            CALL get_tov_metric(coords% levels(l)% var(i,j,k,1), &
                                coords% levels(l)% var(i,j,k,2), &
                                coords% levels(l)% var(i,j,k,3), &
                                tmp, tmp2, tmp3, &
                                g00,g01,g02,g03,g11,g12,g13,g22,g23,g33 )

            CALL compute_tpo_metric &
              ( [g00, g01, g02, g03, g11, g12, g13, g22, g23, g33] , &
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
    PRINT *

 !   DO l= 1, nlevels, 1
 !
 !     lapse2% levels(l)% var         = lapse% levels(l)% var
 !     shift_u2% levels(l)% var       = shift_u% levels(l)% var
 !     Gamma_u2% levels(l)% var       = Gamma_u% levels(l)% var
 !     phi2% levels(l)% var           = phi% levels(l)% var
 !     trK2% levels(l)% var           = trK% levels(l)% var
 !     A_BSSN3_ll2% levels(l)% var    = A_BSSN3_ll% levels(l)% var
 !     g_BSSN3_ll2% levels(l)% var    = g_BSSN3_ll% levels(l)% var
 !     lapse_A_BSSN2% levels(l)% var  = lapse_A_BSSN% levels(l)% var
 !     shift_B_BSSN_u2% levels(l)% var= shift_B_BSSN_u% levels(l)% var
 !     Theta_Z42% levels(l)% var      = Theta_Z4% levels(l)% var
 !
 !   ENDDO
 !
 !   CALL deallocate_BSSN()
 !   CALL deallocate_ADM()
 !   DEALLOCATE(levels)

    !
    !-- Place the ID according to the periastron
    !



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

  END SUBROUTINE read_tov_bssn_id


  PURE SUBROUTINE newtonian_speeds(energy, mass1, mass2, r, v1, v2)

    !***********************************************************
    !
    !# Compute Newtonian speeds for two stars at a given distance,
    !  applying conservation of energy and momentum
    !
    !  FT 13.12.2022
    !
    !***********************************************************

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: energy, mass1, mass2, r
    DOUBLE PRECISION, INTENT(OUT):: v1, v2

    v1= SQRT( 2.D0*mass2*(energy + mass1*mass2/(r**2))/(mass1*(mass1 + mass2)) )

    v2= mass1/mass2*v1

  END SUBROUTINE newtonian_speeds


  PURE FUNCTION total_energy(mass1, mass2, r) RESULT(energy)

    !***********************************************************
    !
    !# Compute Newtonian speeds for two stars at a given distance,
    !  applying conservation of energy and momentum
    !
    !  FT 13.12.2022
    !
    !***********************************************************

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: mass1, mass2, r
    DOUBLE PRECISION:: energy

  END FUNCTION total_energy


END PROGRAM construct_eccentric_binary
