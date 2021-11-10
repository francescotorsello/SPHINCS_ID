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
    !  [[idbase:get_total_spatial_extent]]
    !  are acceptable. It is called by initialize,
    !  after the constructor of the derived type.
    !
    !  FT 8.11.2021
    !
    !************************************************

    IMPLICIT NONE

  END PROCEDURE sanity_check


  MODULE PROCEDURE initialize

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


    IMPLICIT NONE


    CALL derived_type% derived_type_constructor( filename )

    CALL derived_type% sanity_check()

    !foo= derived_type


  END PROCEDURE initialize


END SUBMODULE id_base_initialization
