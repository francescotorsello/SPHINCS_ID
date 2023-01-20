! File:         submodule_bnsfuka_constructor.f90
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

SUBMODULE (bns_fuka) constructor

  !*********************************************************
  !
  !# Implementation of the constructor and
  !  destructor of TYPE [[bnsfuka]], and of the
  !  [[bnsfuka]]-member
  !  PROCEDURES that call the C-bound PROCEDURES
  !  constructig and destructing the |fuka|
  !  |bnsexp| object
  !
  !  FT 23.10.2020
  !
  !*********************************************************


  IMPLICIT NONE


  CONTAINS


  !
  !-- Implementation of the constructor of the bns object
  !
  MODULE PROCEDURE construct_bnsfuka

    !****************************************************
    !
    !# Constructs an object of TYPE [[bnsfuka]]
    !
    !  Created:     FT 09.02.2022
    !  Last update: FT 27.05.2022
    !
    !****************************************************

    USE utility,  ONLY: ten, Msun_geo, flag$tpo

    IMPLICIT NONE

    INTEGER, SAVE:: bns_counter= 1

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: length_scale_pressure


    IF(max$tpo < flag$tpo)THEN
      PRINT *, "** ERROR in MODULE bns_fuka! It must hold max$tpo > flag$tpo,",&
               "where both variables are negative integers. Instead:"
      PRINT *, " * max$tpo=", max$tpo
      PRINT *, " * flag$tpo=", flag$tpo
      PRINT *, " * Stopping..."
      PRINT *
      STOP
    ENDIF

    CALL derived_type% set_n_matter(2)
    CALL derived_type% set_cold_system(.TRUE.)

    derived_type% construction_timer= timer( "binary_construction_timer" )

    ! Construct |fuka| |bnsexp| object
    CALL derived_type% construct_binary( filename )
    derived_type% filename= filename

    ! Read the parameters of the binary system
    CALL read_fuka_id_params( derived_type )

    ! Assign a unique identifier to the bns_fuka object
    derived_type% bns_identifier= bns_counter
    bns_counter= bns_counter + 1

    ! Do not use the geodesic gauge by default
    CALL derived_type% set_one_lapse ( .FALSE. )
    CALL derived_type% set_zero_shift( .FALSE. )

    ! Since Kadath is not thread-safe, we cannot parallelize it using OMP
    ! within SPHINCS_ID. Hence, we chose to make a system call to a program
    ! within Kadath that reads the ID from the FUKA output file and prints it
    ! on a lattice. The ID on the particles will be interplated from this fine
    ! lattice.
    !CALL set_up_lattices_around_stars( filename )

    ! Compute typical length scales of the system using the pressure
    IF( derived_type% get_estimate_length_scale() )THEN

      ALLOCATE( length_scale_pressure(derived_type% get_n_matter()) )
      length_scale_pressure= derived_type% estimate_lengthscale_field( &
                                                get_pressure, &
                                                derived_type% get_n_matter() )

      PRINT *, " * Minimum length scale to resolve on star 1, based on ", &
               "pressure= ", length_scale_pressure(1)*Msun_geo*ten*ten*ten, "m"
      PRINT *, " * Minimum length scale to resolve on star 2, based on ", &
               "pressure= ", length_scale_pressure(1)*Msun_geo*ten*ten*ten, "m"
      PRINT *

    ENDIF

    ! Assign PROCEDURE POINTER to the desired PROCEDURE
    derived_type% finalize_sph_id_ptr => finalize

    CONTAINS

    FUNCTION get_pressure( x, y, z ) RESULT( val )
    !! Returns the value of the pressure at the desired point

      DOUBLE PRECISION, INTENT(IN):: x
      !! \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN):: y
      !! \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN):: z
      !! \(z\) coordinate of the desired point
      DOUBLE PRECISION:: val
      !! Pressure at \((x,y,z)\)

      val= derived_type% read_pressure( x, y, z )

    END FUNCTION get_pressure


  END PROCEDURE construct_bnsfuka


  MODULE PROCEDURE finalize

    !***********************************************
    !
    !#
    !
    !  FT 14.04.2022
    !
    !***********************************************

    IMPLICIT NONE

    ! Temporary implementation, to avoid warnings about unused variables

    pos  = pos
    nlrf = nlrf
    nu   = nu
    pr   = pr
    vel_u= vel_u
    theta= theta
    nstar= nstar
    u    = u

  END PROCEDURE finalize


  MODULE PROCEDURE initialize_id_bnsfuka

    !***********************************************
    !
    !# Initialize the |fuka| |bns| |id|.
    !
    !  - If `flag`= [[utility:flag$sph]], set up the
    !    lattices around the stars for the |bns|
    !    produced with |fuka|.
    !  - If `flag`= [[utility:flag$tpo]], allocate
    !    memory for the hydro grid functions.
    !  - If `flag` > 0, assign its value to [[bnsfuka:l_curr]].
    !  - If [[utility:flag$tpo]] < `flag` < 0,
    !    assign its value to [[bnsfuka:tpo_curr]].
    !
    !  FT 16.09.2022
    !
    !***********************************************

    USE utility,          ONLY: flag$sph, flag$tpo
    USE mesh_refinement,  ONLY: allocate_grid_function

    IMPLICIT NONE

    INTEGER, SAVE:: tpo_counter= 1
    !! Counts how many times the PROCEDURE construct_particles_idase is called

    INTEGER:: i

    CHARACTER(LEN= 3):: cnt_i
    CHARACTER(LEN= :), ALLOCATABLE:: name_mass_density
    CHARACTER(LEN= :), ALLOCATABLE:: name_specific_energy
    CHARACTER(LEN= :), ALLOCATABLE:: name_pressure
    CHARACTER(LEN= :), ALLOCATABLE:: name_v_euler_x
    CHARACTER(LEN= :), ALLOCATABLE:: name_v_euler_y
    CHARACTER(LEN= :), ALLOCATABLE:: name_v_euler_z

    LOGICAL:: wanted_tpo

    IF( PRESENT(switch) )THEN
      IF( switch .EQV. .TRUE. )THEN

        this% tpo$log(tpo_counter)= flag
        this% tpo_curr= tpo_counter

        tpo_counter= tpo_counter + 1

        CALL initialize_id_bnsfuka(this, flag$tpo)

      ENDIF
    ENDIF

    IF( flag /= flag$sph .AND. flag /= flag$tpo .AND. flag < -max$tpo )THEN

      PRINT *, "** ERROR in SUBROUTINE initialize_id_bnsfuka! The INTEGER ", &
               "argument 'flag' should be in the set [1,nlevels], or ", &
               "it should be equal to either flag$sph or flag$tpo, defined ", &
               "in MODULE utility."
      PRINT *, " * Stopping..."
      PRINT *
      STOP

    ENDIF

    IF( flag == flag$sph )THEN

      CALL this% set_up_lattices_around_stars()
      ! Find the surfaces of the stars and print them to a formatted file
      CALL this% find_print_surfaces()

    ELSEIF( flag >= flag$tpo .AND. flag <= -1 )THEN

      wanted_tpo= .FALSE.
      DO i= 1, max$tpo, 1

        IF(this% tpo$log(i) == flag)THEN

          wanted_tpo= .TRUE.
          this% tpo_curr= i

        ENDIF

      ENDDO

      IF(flag == flag$tpo)THEN

        IF( this% tpo_curr <= 9 ) WRITE( cnt_i, "(I1)" ) this% tpo_curr
        IF( this% tpo_curr >= 10 .AND. flag <= 99 ) WRITE( cnt_i, "(I2)" ) &
          this% tpo_curr
        IF( this% tpo_curr >= 100 .AND. flag <= 999 ) WRITE( cnt_i, "(I3)" ) &
          this% tpo_curr

        PRINT *
        PRINT *, "cnt_i=", TRIM(cnt_i)
        PRINT *

        name_mass_density   = "mass_density_fuka-"//TRIM(cnt_i)
        name_specific_energy= "specific_energy_fuka-"//TRIM(cnt_i)
        name_pressure       = "pressure_fuka-"//TRIM(cnt_i)
        name_v_euler_x      = "v_euler_x_fuka-"//TRIM(cnt_i)
        name_v_euler_y      = "v_euler_y_fuka-"//TRIM(cnt_i)
        name_v_euler_z      = "v_euler_z_fuka-"//TRIM(cnt_i)

        PRINT *
        PRINT *, "name_mass_density=", name_mass_density
        PRINT *

        CALL allocate_grid_function( this% mass_density(this% tpo_curr), &
                                     TRIM(name_mass_density), 1 )
        CALL allocate_grid_function( this% specific_energy(this% tpo_curr), &
                                     TRIM(name_specific_energy), 1 )
        CALL allocate_grid_function( this% pressure(this% tpo_curr), &
                                     TRIM(name_pressure), 1 )
        CALL allocate_grid_function( this% v_euler_x(this% tpo_curr), &
                                     TRIM(name_v_euler_x), 1 )
        CALL allocate_grid_function( this% v_euler_y(this% tpo_curr), &
                                     TRIM(name_v_euler_y), 1 )
        CALL allocate_grid_function( this% v_euler_z(this% tpo_curr), &
                                     TRIM(name_v_euler_z), 1 )

      ELSE

        IF(.NOT.wanted_tpo)THEN
          PRINT *, "** ERROR! Mismatch between bns_fuka and bssn objects!"
          PRINT *, " * This should never happen: there is most likely a bug ", &
                   "in SUBROUTINE initialize_id_bnsfuka, or a bug in the ", &
                   "places where it is called."
          PRINT *, " * Stopping..."
          PRINT *
          STOP
        ENDIF

      ENDIF

    ELSE

      this% l_curr= flag

    ENDIF

  END PROCEDURE initialize_id_bnsfuka


  !
  !-- Implementation of the destructor of the bns object
  !
  MODULE PROCEDURE destruct_bnsfuka

    !***********************************************
    !
    !# Destructs an object of TYPE [[bnsfuka]]
    !
    !  Created:     FT 09.02.2022
    !  Last update: FT 30.06.2022
    !
    !***********************************************

    IMPLICIT NONE

    INTEGER:: i_star

    ! Deallocate memory
    CALL this% deallocate_bnsfuka_memory()
    DO i_star=1, 2, 1
      CALL this% star_lattice(i_star)% deallocate_lattice_memory()
    ENDDO

  END PROCEDURE destruct_bnsfuka


  MODULE PROCEDURE construct_binary

    !***********************************************
    !
    !# Construct the |fuka| |bnsexp| object
    !
    !  Created:     FT 09.02.2022
    !  Last update: FT 27.05.2022
    !
    !***********************************************

    IMPLICIT NONE

    LOGICAL:: exist

#ifdef __INTEL_COMPILER

    IF ( C_ASSOCIATED( this% bns_ptr ) ) THEN

      CALL destruct_bns_fuka( this% bns_ptr )

    ENDIF

#endif

    INQUIRE( FILE= fukafile, EXIST= exist )

    IF( exist )THEN

      CALL this% construction_timer% start_timer()
      this% bns_ptr = construct_bns_fuka( fukafile//C_NULL_CHAR )
      CALL this% construction_timer% stop_timer()

    ELSE

      PRINT *, "** ERROR in SUBROUTINE construct_binary: file ", &
               fukafile, " cannot be found!"
      PRINT *
      STOP

    ENDIF


  END PROCEDURE construct_binary


  MODULE PROCEDURE destruct_binary

    !************************************************
    !
    !# Destructs the |fuka| |bnsexp| object and frees
    !  the pointer [[bns:bns_ptr]] pointing to it
    !
    !  Created:     FT 09.02.2022
    !  Last update: FT 27.05.2022
    !
    !************************************************

    IMPLICIT NONE


    IF ( C_ASSOCIATED( this% bns_ptr ) ) THEN

      CALL destruct_bns_fuka( this% bns_ptr )
      this% bns_ptr = C_NULL_PTR

    ENDIF

  END PROCEDURE destruct_binary


END SUBMODULE constructor
