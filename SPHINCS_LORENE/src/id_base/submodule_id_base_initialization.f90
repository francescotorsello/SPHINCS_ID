! File:         submodule_id_base_initialization.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (id_base) id_base_initialization

  !********************************************
  !
  !# Implementation of the methods of TYPE
  !  [[idbase]] that initialize objects of
  !  [[idbase]]-extended TYPE
  !
  !  FT 8.11.2021
  !
  !********************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE sanity_check

    !************************************************
    !
    !# Checks that [[idbase:n_matter]] and the sizes
    !  returned by [[idbase:return_spatial_extent]] and
    !  [[idbase:return_total_spatial_extent]]
    !  are acceptable. It is called by initialize,
    !  after the constructor of the derived type.
    !
    !  FT 8.11.2021
    !
    !************************************************

    IMPLICIT NONE

  END PROCEDURE sanity_check


  MODULE PROCEDURE initialize_idbase

    !************************************************
    !
    !# This PROCEDURE calls the constructor of the
    !  [[idbase]]-extended type and the SUBROUTINE
    !  [[idbase:sanity_check]] afterwards. It is recommended
    !  to use this SUBROUTINE to construct objects of
    !  [[idbase]]-extended type since the sanity check is
    !  performed automatically.
    !
    !  FT 8.11.2021
    !
    !************************************************

    USE bns_lorene,      ONLY: bnslorene
    USE diffstar_lorene, ONLY: diffstarlorene
    USE utility,         ONLY: bnslo, drslo

    IMPLICIT NONE

    IF( ALLOCATED(derived_type) )THEN

      PRINT *, "** ERRORin initialize_idbase! ", &
               " The polymorphic allocatable argument 'derived_type' ",&
               " is already allocated. This SUBROUTINE allocates and", &
               " initializes a polymorphic object of CLASS idbase, hence ", &
               " its argument of CLASS idbase should not be already allocated."
      PRINT *, "   Stopping..."
      PRINT *
      STOP

    ENDIF

    IF( filename(1:5) == bnslo )THEN

      ALLOCATE( bnslorene:: derived_type )

    ELSEIF( filename(1:5) == drslo )THEN

      ALLOCATE( diffstarlorene:: derived_type )

    ELSE
      PRINT *, "** ERROR! Unknown name for the physical system: ", filename(1:5)
      PRINT *, "   Set the variable 'system' in the parameter file ", &
               "sphincs_lorene_parameters.par to one of the values listed there."
      PRINT *, "   Stopping..."
      PRINT *
      STOP
    ENDIF

    CALL derived_type% derived_type_constructor( path//filename )

    CALL derived_type% sanity_check()


  END PROCEDURE initialize_idbase


END SUBMODULE id_base_initialization
