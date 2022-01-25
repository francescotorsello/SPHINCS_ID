! File:         submodule_ejecta_generic_memory.f90
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

    IF(.NOT.ALLOCATED( THIS% grid ))THEN
      ALLOCATE( THIS% grid( THIS% nx_grid, THIS% ny_grid, THIS% nz_grid, 3 ) )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array grid in SUBROUTINE ", &
                  "allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( THIS% baryon_mass_density ))THEN
      ALLOCATE( THIS% baryon_mass_density( THIS% nx_grid, THIS% ny_grid, &
                                           THIS% nz_grid ) )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array baryon_mass_density in ", &
                  "SUBROUTINE allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( THIS% specific_energy ))THEN
      ALLOCATE( THIS% specific_energy( THIS% nx_grid, THIS% ny_grid, &
                                       THIS% nz_grid ) )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array specific_energy in ", &
                  "SUBROUTINE allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( THIS% vel ))THEN
      ALLOCATE( THIS% vel( THIS% nx_grid, THIS% ny_grid, THIS% nz_grid, 3 ) )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array vel in SUBROUTINE ", &
                  "allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( THIS% masses ))THEN
      ALLOCATE( THIS% masses( n_matter ) )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array masses in SUBROUTINE ", &
                  "allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( THIS% sizes ))THEN
      ALLOCATE( THIS% sizes( n_matter, 6 ) )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array sizes in SUBROUTINE ", &
                  "allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( THIS% centers ))THEN
      ALLOCATE( THIS% centers( n_matter, 3 ) )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array centers in SUBROUTINE ", &
                  "allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( THIS% barycenters ))THEN
      ALLOCATE( THIS% barycenters( n_matter, 3 ) )
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

    IF(ALLOCATED( THIS% grid ))THEN
      DEALLOCATE( THIS% grid )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array grid in SUBROUTINE ", &
                  "allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(ALLOCATED( THIS% baryon_mass_density ))THEN
      DEALLOCATE( THIS% baryon_mass_density )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array baryon_mass_density in ", &
                  "SUBROUTINE allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(ALLOCATED( THIS% specific_energy ))THEN
      DEALLOCATE( THIS% specific_energy )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array specific_energy in ", &
                  "SUBROUTINE allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(ALLOCATED( THIS% vel ))THEN
      DEALLOCATE( THIS% vel )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array vel in SUBROUTINE ", &
                  "allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(ALLOCATED( THIS% masses ))THEN
      DEALLOCATE( THIS% masses )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array masses in SUBROUTINE ", &
                  "allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(ALLOCATED( THIS% sizes ))THEN
      DEALLOCATE( THIS% sizes )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array sizes in ", &
                  "SUBROUTINE allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(ALLOCATED( THIS% centers ))THEN
      DEALLOCATE( THIS% centers )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array centers in ", &
                  "SUBROUTINE allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(ALLOCATED( THIS% barycenters ))THEN
      DEALLOCATE( THIS% barycenters )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array barycenters in SUBROUTINE ", &
                  "allocate_gridid_memory.", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF

  END PROCEDURE deallocate_gridid_memory


END SUBMODULE memory
