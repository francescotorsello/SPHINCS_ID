! File:         submodule_bnslorene_constructor.f90
! Authors:      Francesco Torsello (FT)
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

SUBMODULE (bns_lorene) constructor

  !*********************************************************
  !
  !# Implementation of the constructor and
  !  destructor of TYPE [[bnslorene]], and of the
  !  [[bnslorene]]-member
  !  PROCEDURES that call the C-bound PROCEDURES
  !  constructig and destructing the |lorene|
  !  |binns| object
  !
  !  FT 23.10.2020
  !
  !*********************************************************


  IMPLICIT NONE


  CONTAINS


  !
  !-- Implementation of the constructor of the bnslorene object
  !
  MODULE PROCEDURE construct_bnslorene

    !****************************************************
    !
    !# Constructs an object of TYPE [[bnslorene]]
    !
    !  FT
    !
    !****************************************************

    USE utility,  ONLY: ten, Msun_geo, use_eos_from_id, common_eos_path

    IMPLICIT NONE

    INTEGER, SAVE:: bns_counter= 1

    INTEGER, PARAMETER:: n_matter= 2

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: length_scale_pressure

    CALL derived_type% set_n_matter(n_matter)
    CALL derived_type% set_cold_system(.TRUE.)

    derived_type% eos_filenames(1)= &
      TRIM(common_eos_path)//TRIM(eos_filenames(1))
    derived_type% eos_filenames(2)= &
      TRIM(common_eos_path)//TRIM(eos_filenames(2))

    derived_type% construction_timer= timer("binary_construction_timer")

    ! Construct |lorene| |binns| object
    IF( PRESENT(filename) .AND. use_eos_from_id )THEN

      CALL derived_type% construct_binary(filename, ["use_id","use_id"])

    ELSEIF( PRESENT(filename) .AND. .NOT.use_eos_from_id )THEN

      CALL derived_type% construct_binary(filename, &
                                          derived_type% eos_filenames)

    ELSE

      !CALL derived_type% construct_binary()
      STOP

    ENDIF
    ! Import the properties of the BNS
    CALL read_bns_properties(derived_type)

    ! Assign a unique identifier to the bnslorene object
    derived_type% bns_identifier= bns_counter
    bns_counter= bns_counter + 1

    ! Do not use the geodesic gauge by default
    CALL derived_type% set_one_lapse (.FALSE.)
    CALL derived_type% set_zero_shift(.FALSE.)

    ! Find the surfaces of the stars and print them to a formatted file
    CALL derived_type% find_print_surfaces()

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
    derived_type% finalize_sph_id_ptr => correct_adm_linear_momentum


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

  END PROCEDURE construct_bnslorene


  MODULE PROCEDURE nothing

    !***********************************************
    !
    !# Procedure that does nothing. It is used to instantiate a deferred
    !  idbase procedure which is not needed in TYPE [[bnslorene]].
    !  It also serves as a placeholder in case the idbase procedure
    !  will be needed in the future.
    !
    !  FT 15.09.2022
    !
    !***********************************************

    IMPLICIT NONE

  END PROCEDURE nothing


  !
  !-- Implementation of the destructor of the bnslorene object
  !
  MODULE PROCEDURE destruct_bnslorene

    !***********************************************
    !
    !# Destructs an object of TYPE [[bnslorene]]
    !
    !  FT
    !
    !***********************************************

    IMPLICIT NONE

    !PRINT *, "Inside destructor of bns."
    !PRINT *

    ! Deallocate memory
    CALL this% deallocate_bnslorene_memory()

  END PROCEDURE destruct_bnslorene


  MODULE PROCEDURE construct_binary

    !***********************************************
    !
    !# Construct the |lorene| |binns| object
    !
    !  FT
    !
    !***********************************************

    IMPLICIT NONE

    !CHARACTER(KIND=C_CHAR, LEN= 7):: default_case
    CHARACTER(KIND=C_CHAR, LEN= :), ALLOCATABLE:: eos1, eos2
    LOGICAL:: exist

    !PRINT *, "** Executing the construct_binary subroutine..."

    eos1= TRIM(eos_filenames(1))//C_NULL_CHAR
    eos2= TRIM(eos_filenames(2))//C_NULL_CHAR

#ifdef __INTEL_COMPILER

    IF ( C_ASSOCIATED( this% bns_ptr ) ) THEN

      CALL destruct_bin_ns( this% bns_ptr )

    ENDIF

#endif

    !
    !-- If the name of the |lorene| binary file id_file is given as argument to
    !-- construct_binary, use it. Otherwise, give the string "read_it"
    !-- to construct_bin_ns as argument, which makes |lorene| read the name of
    !-- the file from the parameter file read_bin_ns.par
    !
    IF( PRESENT( id_file ) )THEN

      INQUIRE( FILE= id_file, EXIST= exist )

      IF( exist )THEN

        CALL this% construction_timer% start_timer()
        this% bns_ptr = construct_bin_ns(id_file//C_NULL_CHAR, eos1, eos2)
        CALL this% construction_timer% stop_timer()

      ELSE

        PRINT *, "** ERROR in SUBROUTINE construct_binary: file ", &
                 id_file, " cannot be found!"
        PRINT *
        STOP

      ENDIF

    ELSE

      !default_case= "read_it"
      !CALL this% construction_timer% start_timer()
      !this% bns_ptr = construct_bin_ns( default_case//C_NULL_CHAR )
      !CALL this% construction_timer% stop_timer()
      STOP

    ENDIF

    !PRINT *, "** Subroutine construct_binary executed."
    !PRINT *

  END PROCEDURE construct_binary


  MODULE PROCEDURE destruct_binary

    !************************************************
    !
    !# Destructs the |lorene| |binns| object and frees
    !  the pointer [[bnslorene:bns_ptr]] pointing to it
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    !PRINT *, "** Executing the destruct_binary subroutine."

    IF ( C_ASSOCIATED( this% bns_ptr ) ) THEN

      CALL destruct_bin_ns( this% bns_ptr )
      this% bns_ptr = C_NULL_PTR

    ENDIF

    !PRINT *, "** Subroutine destruct_binary executed."
    !PRINT *

  END PROCEDURE destruct_binary


END SUBMODULE constructor
