! File:         module_sphincs_id_interpolate.f90
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

MODULE sphincs_id_interpolate


  !*********************************************
  !
  !# This module contains data and PROCEDURES
  !  needed to set up |id| prepared  on a grid in
  !  PROGRAM [[sphincs_id]]
  !
  !  FT 19.11.2020
  !
  !*********************************************


  USE id_base,         ONLY: idbase
  USE ejecta_generic,  ONLY: ejecta


  IMPLICIT NONE


  CHARACTER( LEN= 5 ), PARAMETER:: ejecta_grid= "EJECT"
  !# String that identifies an ejecta prepared on a uniform Cartesian grid


  CONTAINS


  SUBROUTINE allocate_idbase( id, filename, system, system_name )

    !*********************************************
    !
    !# This SUBROUTINE allocates a polymorphic
    !  object of class idbase to its dynamic type.
    !  The dynamic type is one among those that
    !  interpolate the ID from a grid
    !
    !  FT 19.11.2020
    !
    !*********************************************

    IMPLICIT NONE

    CLASS( idbase ), ALLOCATABLE, INTENT( IN OUT ):: id
    CHARACTER(LEN=*), INTENT( IN ) :: filename
    CHARACTER(LEN=5), INTENT( IN OUT ):: system
    CHARACTER(LEN=5), INTENT( IN OUT ):: system_name

    IF( ALLOCATED(id) )THEN

      PRINT *, "** ERROR in allocate_idbase! ", &
               " The polymorphic allocatable argument 'id' ",&
               " is already allocated. This SUBROUTINE allocates and", &
               " initializes a polymorphic object of CLASS idbase, hence ", &
               " its argument of CLASS idbase should not be already allocated."
      PRINT *, "   Stopping..."
      PRINT *
      STOP

    ENDIF

    IF( filename(1:5) == ejecta_grid )THEN

      ALLOCATE( ejecta:: id )
      system= ejecta_grid
      system_name= "EJEC."

    ELSE

      PRINT *, "** ERROR! Unrecognized physical system ", system
      PRINT *
      PRINT *, "   Please specify the type of physical system in the first 5",&
               " characters of the name of the file containing the initial", &
               " data."
      PRINT *
      PRINT *, "   The 5-character names, and associated physical systems,", &
              " supported by this version of SPHINCS_ID, with flavour = 4, are:"
      PRINT *
      PRINT *, "   EJECT: Ejecta data on a uniform Cartesian grid"
      PRINT *
      STOP

    ENDIF

  END SUBROUTINE allocate_idbase


END MODULE sphincs_id_interpolate
