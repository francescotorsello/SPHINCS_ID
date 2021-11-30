! File:         submodule_bssn_id_memory.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (formul_bssn_id) bssn_id_memory

  !************************************************
  !
  !# Implementation of the methods of TYPE bssn_id
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

    IF( .NOT.ALLOCATED( THIS% Gamma_u% levels ) )THEN
      CALL allocate_grid_function( THIS% Gamma_u, "Gamma_u_id", 3 )
    ENDIF

    IF( .NOT.ALLOCATED( THIS% phi% levels ) )THEN
      CALL allocate_grid_function( THIS% phi, "phi_id", 1 )
    ENDIF

    IF( .NOT.ALLOCATED( THIS% trK% levels ) )THEN
      CALL allocate_grid_function( THIS% trK, "trK_id", 1 )
    ENDIF

    IF( .NOT.ALLOCATED( THIS% A_BSSN3_ll% levels ) )THEN
      CALL allocate_grid_function( THIS% A_BSSN3_ll, "A_BSSN3_ll_id", 6 )
    ENDIF

    IF( .NOT.ALLOCATED( THIS% g_BSSN3_ll% levels ) )THEN
      CALL allocate_grid_function( THIS% g_BSSN3_ll, "g_BSSN3_ll_id", 6 )
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

    IF( ALLOCATED( THIS% Gamma_u% levels ) )THEN
      CALL deallocate_grid_function( THIS% Gamma_u, "Gamma_u_id" )
    ENDIF

    IF( ALLOCATED( THIS% phi% levels ) )THEN
      CALL deallocate_grid_function( THIS% phi, "phi_id" )
    ENDIF

    IF( ALLOCATED( THIS% trK% levels ) )THEN
      CALL deallocate_grid_function( THIS% trK, "trK_id" )
    ENDIF

    IF( ALLOCATED( THIS% A_BSSN3_ll% levels ) )THEN
      CALL deallocate_grid_function( THIS% A_BSSN3_ll, "A_BSSN3_ll_id" )
    ENDIF

    IF( ALLOCATED( THIS% g_BSSN3_ll% levels ) )THEN
      CALL deallocate_grid_function( THIS% g_BSSN3_ll, "g_BSSN3_ll_id" )
    ENDIF

    IF( ALLOCATED( THIS% GC% levels ) )THEN
      CALL deallocate_grid_function( THIS% GC, "GC_id" )
    ENDIF

  END PROCEDURE deallocate_bssn_fields


END SUBMODULE bssn_id_memory
