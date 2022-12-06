! File:         submodule_diffstar_lorene_constructor.f90
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

    IMPLICIT NONE

    INTEGER, SAVE:: diffstar_counter= 1

    CALL derived_type% set_n_matter(1)
    CALL derived_type% set_cold_system(.TRUE.)

    derived_type% construction_timer= timer( "drs_construction_timer" )

    ! Construct |lorene| |etdiffrot| object
    IF( PRESENT( filename ) )THEN
        CALL derived_type% construct_drs( filename )
    ELSE
        CALL derived_type% construct_drs()
    ENDIF

    ! Import the parameters of the binary system
    CALL read_diffstar_params( derived_type )

    ! Assign a unique identifier to the bns object
    derived_type% diffstar_identifier= diffstar_counter
    diffstar_counter= diffstar_counter + 1

    ! Do not use the geodesic gauge by default
    CALL derived_type% set_one_lapse ( .FALSE. )
    CALL derived_type% set_zero_shift( .FALSE. )

    derived_type% finalize_sph_id_ptr => finalize

  END PROCEDURE construct_diffstarlorene


  MODULE PROCEDURE finalize

    !***********************************************
    !
    !#
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

    CHARACTER(KIND= C_CHAR, LEN= 7):: default_case
    LOGICAL:: exist

    !PRINT *, "** Executing the construct_binary subroutine..."

#ifdef __INTEL_COMPILER

    IF ( C_ASSOCIATED( this% diffstar_ptr ) ) THEN

      CALL destruct_etdiffrot( this% diffstar_ptr )

    ENDIF

#endif

    !
    !-- If the name of the |lorene| binary file resu_file is given as argument to
    !-- construct_binary, use it. Otherwise, give the string "read_it"
    !-- to construct_drs as argument, which makes |lorene| read the name of
    !-- the file from the parameter file read_bin_ns.par
    !
    IF( PRESENT( resu_file ) )THEN

      INQUIRE( FILE= resu_file, EXIST= exist )

      IF( exist )THEN

        CALL this% construction_timer% start_timer()
        this% diffstar_ptr = construct_etdiffrot( resu_file//C_NULL_CHAR )
        CALL this% construction_timer% stop_timer()

      ELSE

        PRINT *, "** ERROR in SUBROUTINE construct_binary: file ", &
                 resu_file, " cannot be found!"
        PRINT *
        STOP

      ENDIF

    ELSE

      default_case= "read_it"
      CALL this% construction_timer% start_timer()
      this% diffstar_ptr = construct_etdiffrot( default_case//C_NULL_CHAR )
      CALL this% construction_timer% stop_timer()

    ENDIF

  END PROCEDURE construct_drs


  MODULE PROCEDURE destruct_drs

    !************************************************
    !
    !# Destructs the |lorene| |etdiffrot| object and frees
    !  the pointer [[diffstar:diffstar_ptr]] pointing to it
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
