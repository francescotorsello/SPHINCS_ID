! File:         module_bns_base.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

MODULE bns_base

  !*******************************************************
  !
  !#
  !
  !   FT 24.09.2021
  !
  !*******************************************************


  USE id_base, ONLY: idbase


  IMPLICIT NONE


  !*******************************************************
  !                                                      *
  !     Definition of TYPE bns (binary neutron star)     *
  !                                                      *
  !   This class imports and stores the LORENE BNS ID    *
  !                                                      *
  !*******************************************************

  TYPE, ABSTRACT, EXTENDS(idbase):: bnsbase

  END TYPE bnsbase

END MODULE bns_base

