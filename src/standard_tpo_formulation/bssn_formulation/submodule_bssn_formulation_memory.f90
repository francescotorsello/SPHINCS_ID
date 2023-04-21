! File:         submodule_bssn_formulation_memory.f90
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

SUBMODULE (bssn_formulation) memory

  !************************************************
  !
  !# Implementation of the methods of TYPE [[bssn]]
  !  that (de)allocate memory
  !
  !  FT 9.07.2021
  !
  !************************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE allocate_bssn_fields

    !***********************************************
    !
    !# Allocate memory for the BSSN variables.
    !
    !  FT 23.10.2020
    !
    !  Updated to support mesh refinement
    !
    !  FT 26.03.2021
    !
    !***********************************************

    USE mesh_refinement,  ONLY: allocate_grid_function

    IMPLICIT NONE

    CHARACTER(LEN= 2):: tpo_id

    WRITE( tpo_id, "(I2)" ) this% tpo_id_number

    IF( .NOT.ALLOCATED( this% Gamma_u% levels ) )THEN
      CALL allocate_grid_function( this% Gamma_u, "Gamma_u_id"//tpo_id, 3 )
    ENDIF

    IF( .NOT.ALLOCATED( this% phi% levels ) )THEN
      CALL allocate_grid_function( this% phi, "phi_id"//tpo_id, 1 )
    ENDIF

    IF( .NOT.ALLOCATED( this% trK% levels ) )THEN
      CALL allocate_grid_function( this% trK, "trK_id"//tpo_id, 1 )
    ENDIF

    IF( .NOT.ALLOCATED( this% A_BSSN3_ll% levels ) )THEN
      CALL allocate_grid_function( this% A_BSSN3_ll, "A_BSSN3_ll_id" &
        //tpo_id, 6 )
    ENDIF

    IF( .NOT.ALLOCATED( this% g_BSSN3_ll% levels ) )THEN
      CALL allocate_grid_function( this% g_BSSN3_ll, "g_BSSN3_ll_id" &
        //tpo_id, 6 )
    ENDIF

  END PROCEDURE allocate_bssn_fields


  MODULE PROCEDURE deallocate_bssn_fields

    !**************************************************
    !
    !# Deallocate BSSN memory
    !
    !  FT
    !
    !**************************************************

    USE mesh_refinement, ONLY: deallocate_grid_function

    IMPLICIT NONE

    CHARACTER(LEN= 2):: tpo_id

    WRITE( tpo_id, "(I2)" ) this% tpo_id_number

    IF( ALLOCATED( this% Gamma_u% levels ) )THEN
      CALL deallocate_grid_function( this% Gamma_u, "Gamma_u_id"//tpo_id )
    ENDIF

    IF( ALLOCATED( this% phi% levels ) )THEN
      CALL deallocate_grid_function( this% phi, "phi_id"//tpo_id )
    ENDIF

    IF( ALLOCATED( this% trK% levels ) )THEN
      CALL deallocate_grid_function( this% trK, "trK_id"//tpo_id )
    ENDIF

    IF( ALLOCATED( this% A_BSSN3_ll% levels ) )THEN
      CALL deallocate_grid_function( this% A_BSSN3_ll, "A_BSSN3_ll_id"//tpo_id )
    ENDIF

    IF( ALLOCATED( this% g_BSSN3_ll% levels ) )THEN
      CALL deallocate_grid_function( this% g_BSSN3_ll, "g_BSSN3_ll_id"//tpo_id )
    ENDIF

    IF( ALLOCATED( this% GC% levels ) )THEN
      CALL deallocate_grid_function( this% GC, "GC_id"//tpo_id )
    ENDIF

    IF( ALLOCATED( this% GC_parts% levels ) )THEN
      CALL deallocate_grid_function( this% GC, "GC_parts_id"//tpo_id )
    ENDIF

    IF( ALLOCATED( this% Ricci_ll% levels ) )THEN
      CALL deallocate_grid_function( this% Ricci_ll, "Ricci_ll_id"//tpo_id )
    ENDIF

    IF( ALLOCATED( this% Ricci_scalar% levels ) )THEN
      CALL deallocate_grid_function( this% Ricci_scalar, &
                                     "Ricci_scalar_id"//tpo_id )
    ENDIF

    !IF( ALLOCATED( this% rho% levels ) )THEN
    !  CALL deallocate_grid_function( this% rho, "rho" )
    !ENDIF
    !
    !IF( ALLOCATED( this% S% levels ) )THEN
    !  CALL deallocate_grid_function( this% S, "S" )
    !ENDIF

    !IF( ALLOCATED( this% rho_parts% levels ) )THEN
    !  CALL deallocate_grid_function( this% rho_parts, "rho_parts" )
    !ENDIF
    !
    !IF( ALLOCATED( this% S_parts% levels ) )THEN
    !  CALL deallocate_grid_function( this% S_parts, "S_parts" )
    !ENDIF

  END PROCEDURE deallocate_bssn_fields


END SUBMODULE memory
