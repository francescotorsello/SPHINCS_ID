! File:         submodule_ejecta_generic_constructor.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (ejecta_generic) ejecta_generic_constructor

  !*********************************************************
  !
  !# Implementation of the constructor and
  !  destructor of TYPE [[diffstarlorene]], and of the
  !  [[diffstarlorene]]-member
  !  PROCEDURES that call the C-bound PROCEDURES
  !  constructig and destructing the |lorene|
  !  |etdiffrot| object
  !
  !  FT 19.11.2021
  !
  !*********************************************************


  IMPLICIT NONE


  CONTAINS


  !
  !-- Implementation of the constructor of the bns object
  !
  MODULE PROCEDURE construct_ejecta

    !****************************************************
    !
    !# Constructs an object of TYPE [[diffstarlorene]]
    !
    !  FT 19.11.2021
    !
    !****************************************************

    IMPLICIT NONE


  END PROCEDURE construct_ejecta


  !
  !-- Implementation of the destructor of the bns object
  !
  MODULE PROCEDURE destruct_ejecta

    !****************************************************
    !
    !# Destructs an object of TYPE [[diffstarlorene]]
    !
    !  FT 19.11.2021
    !
    !****************************************************

    IMPLICIT NONE


  END PROCEDURE destruct_ejecta


END SUBMODULE ejecta_generic_constructor
