! File:         submodule_ejecta_generic_memory.f90
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

SUBMODULE (ejecta_generic) memory

  !***********************************************
  !
  !# Implementation of the methods of TYPE ejecta
  !  that (de)allocate memory
  !
  ! FT 14.01.2022
  !
  !***********************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE allocate_gridid_memory

    !***********************************************
    !
    !# Allocate the memory to store the ID
    !  in the member arrays
    !
    !  FT 14.01.2022
    !
    !***********************************************

    IMPLICIT NONE

    IF(.NOT.ALLOCATED( this% grid ))THEN
      ALLOCATE( this% grid( this% nx_grid, this% ny_grid, this% nz_grid, 3 ) )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array grid in SUBROUTINE ", &
                  "allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( this% baryon_mass_density ))THEN
      ALLOCATE( this% baryon_mass_density( this% nx_grid, this% ny_grid, &
                                           this% nz_grid ) )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array baryon_mass_density in ", &
                  "SUBROUTINE allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( this% specific_energy ))THEN
      ALLOCATE( this% specific_energy( this% nx_grid, this% ny_grid, &
                                       this% nz_grid ) )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array specific_energy in ", &
                  "SUBROUTINE allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( this% vel ))THEN
      ALLOCATE( this% vel( this% nx_grid, this% ny_grid, this% nz_grid, 3 ) )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array vel in SUBROUTINE ", &
                  "allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( this% masses ))THEN
      ALLOCATE( this% masses( n_matter ) )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array masses in SUBROUTINE ", &
                  "allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( this% sizes ))THEN
      ALLOCATE( this% sizes( n_matter, 6 ) )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array sizes in SUBROUTINE ", &
                  "allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( this% centers ))THEN
      ALLOCATE( this% centers( n_matter, 3 ) )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array centers in SUBROUTINE ", &
                  "allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( this% barycenters ))THEN
      ALLOCATE( this% barycenters( n_matter, 3 ) )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array barycenters in SUBROUTINE ", &
                  "allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF

  END PROCEDURE allocate_gridid_memory


  MODULE PROCEDURE deallocate_gridid_memory

    !***********************************************
    !
    !# Deallocate the memory to store the ID
    !  in the member arrays
    !
    !  FT 14.01.2022
    !
    !***********************************************

    IMPLICIT NONE

    IF(ALLOCATED( this% grid ))THEN
      DEALLOCATE( this% grid )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array grid in SUBROUTINE ", &
                  "allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(ALLOCATED( this% baryon_mass_density ))THEN
      DEALLOCATE( this% baryon_mass_density )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array baryon_mass_density in ", &
                  "SUBROUTINE allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(ALLOCATED( this% specific_energy ))THEN
      DEALLOCATE( this% specific_energy )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array specific_energy in ", &
                  "SUBROUTINE allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(ALLOCATED( this% vel ))THEN
      DEALLOCATE( this% vel )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array vel in SUBROUTINE ", &
                  "allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(ALLOCATED( this% masses ))THEN
      DEALLOCATE( this% masses )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array masses in SUBROUTINE ", &
                  "allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(ALLOCATED( this% sizes ))THEN
      DEALLOCATE( this% sizes )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array sizes in ", &
                  "SUBROUTINE allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(ALLOCATED( this% centers ))THEN
      DEALLOCATE( this% centers )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array centers in ", &
                  "SUBROUTINE allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(ALLOCATED( this% barycenters ))THEN
      DEALLOCATE( this% barycenters )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array barycenters in SUBROUTINE ", &
                  "allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF

  END PROCEDURE deallocate_gridid_memory


END SUBMODULE memory
