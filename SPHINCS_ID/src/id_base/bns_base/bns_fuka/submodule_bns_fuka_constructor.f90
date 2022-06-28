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
  !  |binns| object
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

    USE utility,  ONLY: ten, Msun_geo

    IMPLICIT NONE

    INTEGER, SAVE:: bns_counter= 1

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: length_scale_pressure

    CALL derived_type% set_n_matter(2)
    CALL derived_type% set_cold_system(.TRUE.)

    derived_type% construction_timer= timer( "binary_construction_timer" )

    ! Construct |fuka| bns_export object
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

      val= derived_type% read_fuka_pressure( x, y, z )

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


  !
  !-- Implementation of the destructor of the bns object
  !
  MODULE PROCEDURE destruct_bnsfuka

    !***********************************************
    !
    !# Destructs an object of TYPE [[bnsfuka]]
    !
    !  Created:     FT 09.02.2022
    !  Last update: FT 27.05.2022
    !
    !***********************************************

    IMPLICIT NONE

    ! Deallocate memory
    CALL this% deallocate_bnsfuka_memory()

  END PROCEDURE destruct_bnsfuka


  MODULE PROCEDURE construct_binary

    !***********************************************
    !
    !# Construct the |fuka| bns_export object
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
    !# Destructs the |fuka| bns_export object and frees
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
