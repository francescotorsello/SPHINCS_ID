! File:         submodule_diffstar_lorene_constructor.f90
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

SUBMODULE (diffstar_lorene) constructor

  !*********************************************************
  !
  !# Implementation of the constructor and
  !  destructor of TYPE [[diffstarlorene]], and of the
  !  [[diffstarlorene]]-member
  !  PROCEDURES that call the C-bound PROCEDURES
  !  constructig and destructing the |lorene|
  !  |etdiffrot| object
  !
  !  FT 25.10.2021
  !
  !*********************************************************


  USE, INTRINSIC:: ISO_C_BINDING, ONLY: C_ASSOCIATED, C_NULL_CHAR, C_NULL_PTR


  IMPLICIT NONE


  CONTAINS


  !
  !-- Implementation of the constructor of the bns object
  !
  MODULE PROCEDURE construct_diffstarlorene

    !****************************************************
    !
    !# Constructs an object of TYPE [[diffstarlorene]]
    !
    !  FT 25.10.2021
    !
    !****************************************************

    USE utility,  ONLY: use_eos_from_id, common_eos_path

    IMPLICIT NONE

    INTEGER, SAVE:: diffstar_counter= 1

    CALL derived_type% set_n_matter(1)
    CALL derived_type% set_cold_system(.TRUE.)

    derived_type% eos_filename(1)= &
      TRIM(common_eos_path)//TRIM(eos_filenames(1))

    derived_type% construction_timer= timer("drs_construction_timer")

    ! Construct |lorene| |etdiffrot| object
    IF( PRESENT(filename) .AND. use_eos_from_id )THEN

      CALL derived_type% construct_drs(filename, "use_id")

    ELSEIF( PRESENT(filename) .AND. .NOT.use_eos_from_id )THEN

      CALL derived_type% construct_drs(filename, derived_type% eos_filename(1))

    ELSE

      !CALL derived_type% construct_drs()
      STOP

    ENDIF

    ! Import the properties of the differentially rotating star
    CALL read_diffstar_properties(derived_type)

    ! Assign a unique identifier to the bns object
    derived_type% diffstar_identifier= diffstar_counter
    diffstar_counter= diffstar_counter + 1

    ! Do not use the geodesic gauge by default
    CALL derived_type% set_one_lapse (.FALSE.)
    CALL derived_type% set_zero_shift(.FALSE.)

    derived_type% finalize_sph_id_ptr => finalize

  END PROCEDURE construct_diffstarlorene


  MODULE PROCEDURE finalize

    !***********************************************
    !
    !# This SUBROUTINE is curretly just a placeholder.
    !  It could be used, for exmaple, to correct
    !  the ADM momentum at the end of the execution,
    !  or correct for residual eccentricity, etc.
    !
    !  @note Temporary implementation, to avoid warnings
    !        about unused variables.
    !
    !  FT 14.04.2022
    !
    !***********************************************

    IMPLICIT NONE

    ! Temporary implementation, to avoid warnings about unuesed variables

    pos  = pos
    nlrf = nlrf
    nu   = nu
    pr   = pr
    vel_u= vel_u
    theta= theta
    nstar= nstar
    u    = u


  END PROCEDURE finalize


  MODULE PROCEDURE nothing

    !***********************************************
    !
    !# Procedure that does nothing. It is used to instantiate a deferred
    !  idbase procedure which s not needed in TYPE [[ejecta_generic]].
    !  It also serves as a placeholder in case the idbase procedure
    !  will be needed in the future.
    !
    !  FT 15.09.2022
    !
    !***********************************************

    IMPLICIT NONE

  END PROCEDURE nothing


  !
  !-- Implementation of the destructor of the bns object
  !
  MODULE PROCEDURE destruct_diffstarlorene

    !***********************************************
    !
    !# Destructs an object of TYPE [[diffstarlorene]]
    !
    !  FT 25.10.2021
    !
    !***********************************************

    IMPLICIT NONE

    ! Deallocate memory
    CALL this% deallocate_diffstar_memory()

  END PROCEDURE destruct_diffstarlorene


  MODULE PROCEDURE construct_drs

    !***********************************************
    !
    !# Construct the |lorene| |etdiffrot| object
    !
    !  FT 25.10.2021
    !
    !***********************************************

    IMPLICIT NONE

    !CHARACTER(KIND= C_CHAR, LEN= 7):: default_case
    CHARACTER(KIND=C_CHAR, LEN= :), ALLOCATABLE:: eos
    LOGICAL:: exist

    !PRINT *, "** Executing the construct_binary subroutine..."

    eos= TRIM(eos_filename)//C_NULL_CHAR

#ifdef __INTEL_COMPILER

    IF ( C_ASSOCIATED( this% diffstar_ptr ) ) THEN

      CALL destruct_etdiffrot( this% diffstar_ptr )

    ENDIF

#endif

    !
    !-- If the name of the |lorene| binary file id_file is given as argument to
    !-- construct_binary, use it. Otherwise, give the string "read_it"
    !-- to construct_drs as argument, which makes |lorene| read the name of
    !-- the file from the parameter file read_bin_ns.par
    !
    IF( PRESENT( id_file ) )THEN

      INQUIRE( FILE= id_file, EXIST= exist )

      IF( exist )THEN

        CALL this% construction_timer% start_timer()
        this% diffstar_ptr = construct_etdiffrot( id_file//C_NULL_CHAR, eos )
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
      !this% diffstar_ptr = construct_etdiffrot( default_case//C_NULL_CHAR, &
      !                                TRIM(eos_filename)//C_NULL_CHAR ) )
      !CALL this% construction_timer% stop_timer()
      STOP

    ENDIF

  END PROCEDURE construct_drs


  MODULE PROCEDURE destruct_drs

    !************************************************
    !
    !# Destructs the |lorene| |etdiffrot| object and frees
    !  the pointer [[diffstarlorene:diffstar_ptr]] pointing to it
    !
    !  FT 25.10.2021
    !
    !************************************************

    IMPLICIT NONE

    IF ( C_ASSOCIATED( this% diffstar_ptr ) ) THEN

      CALL destruct_etdiffrot( this% diffstar_ptr )
      this% diffstar_ptr = C_NULL_PTR

    ENDIF

  END PROCEDURE destruct_drs


END SUBMODULE constructor
