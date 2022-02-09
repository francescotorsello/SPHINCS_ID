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


 ! MODULE PROCEDURE construct_bnsfuka2
 !
 !   !****************************************************
 !   !
 !   !# Constructs an object of TYPE [[bnsfuka]]
 !   !
 !   !  FT
 !   !
 !   !****************************************************
 !
 !   IMPLICIT NONE
 !
 !   CHARACTER(LEN=10) :: resu_file
 !
 !   derived_type => construct_bnsfuka( resu_file )
 !
 ! END PROCEDURE construct_bnsfuka2


  !
  !-- Implementation of the constructor of the bns object
  !
  MODULE PROCEDURE construct_bnsfuka

    !****************************************************
    !
    !# Constructs an object of TYPE [[bnsfuka]]
    !
    !  FT 09.02.2022
    !
    !****************************************************

    IMPLICIT NONE

!    INTEGER, SAVE:: bns_counter= 1
!
!    CALL derived_type% set_n_matter(2)
!    CALL derived_type% set_cold_system(.TRUE.)
!
!    derived_type% construction_timer= timer( "binary_construction_timer" )
!
!    ! Construct |fuka| |binns| object
!    IF( PRESENT( filename ) )THEN
!        CALL derived_type% construct_binary( filename )
!    ELSE
!        CALL derived_type% construct_binary()
!    ENDIF
!
!    ! Import the parameters of the binary system
!    CALL import_id_params( derived_type )
!
!    ! Assign a unique identifier to the bns object
!    derived_type% bns_identifier= bns_counter
!    bns_counter= bns_counter + 1
!
!    ! Do not use the geodesic gauge by default
!    CALL derived_type% set_one_lapse ( .FALSE. )
!    CALL derived_type% set_zero_shift( .FALSE. )

  END PROCEDURE construct_bnsfuka


  !
  !-- Implementation of the destructor of the bns object
  !
  MODULE PROCEDURE destruct_bnsfuka

    !***********************************************
    !
    !# Destructs an object of TYPE [[bnsfuka]]
    !
    !  FT 09.02.2022
    !
    !***********************************************

    IMPLICIT NONE

    ! Deallocate memory
  !  CALL THIS% deallocate_bnsfuka_memory()

  END PROCEDURE destruct_bnsfuka


  MODULE PROCEDURE construct_binary

    !***********************************************
    !
    !# Construct the |fuka| ?? object
    !
    !  FT 09.02.2022
    !
    !***********************************************

    IMPLICIT NONE

!    CHARACTER(KIND= C_CHAR, LEN= 7):: default_case
!    LOGICAL:: exist
!
!    !PRINT *, "** Executing the construct_binary subroutine..."
!
!#ifdef __INTEL_COMPILER
!
!    IF ( C_ASSOCIATED( THIS% bns_ptr ) ) THEN
!
!      CALL destruct_bin_ns( THIS% bns_ptr )
!
!    ENDIF
!
!#endif
!
!    !
!    !-- If the name of the |fuka| binary file resu_file is given as argument to
!    !-- construct_binary, use it. Otherwise, give the string "read_it"
!    !-- to construct_bin_ns as argument, which makes |fuka| read the name of
!    !-- the file from the parameter file read_bin_ns.par
!    !
!    IF( PRESENT( resu_file ) )THEN
!
!      INQUIRE( FILE= resu_file, EXIST= exist )
!
!      IF( exist )THEN
!
!        CALL THIS% construction_timer% start_timer()
!        THIS% bns_ptr = construct_bin_ns( resu_file//C_NULL_CHAR )
!        CALL THIS% construction_timer% stop_timer()
!
!      ELSE
!
!        PRINT *, "** ERROR in SUBROUTINE construct_binary: file ", &
!                 resu_file, " cannot be found!"
!        PRINT *
!        STOP
!
!      ENDIF
!
!    ELSE
!
!      default_case= "read_it"
!      CALL THIS% construction_timer% start_timer()
!      THIS% bns_ptr = construct_bin_ns( default_case//C_NULL_CHAR )
!      CALL THIS% construction_timer% stop_timer()
!
!    ENDIF
!

  END PROCEDURE construct_binary


  MODULE PROCEDURE destruct_binary

    !************************************************
    !
    !# Destructs the |fuka| ?? object and frees
    !  the pointer [[bns:bns_ptr]] pointing to it
    !
    !  FT 09.02.2022
    !
    !************************************************

    IMPLICIT NONE


 !   IF ( C_ASSOCIATED( THIS% bns_ptr ) ) THEN
 !
 !     CALL destruct_bin_ns( THIS% bns_ptr )
 !     THIS% bns_ptr = C_NULL_PTR
 !
 !   ENDIF

  END PROCEDURE destruct_binary


END SUBMODULE constructor
