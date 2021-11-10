! File:         module_sphincs_lorene.f90
! Author:       Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

MODULE sphincs_lorene


  !*********************************************
  !
  !# This module contains data and PROCEDURES
  !  needed to set up the |lorene| |id| in
  !  PROGRAM [[sphincs_id_lorene]]
  !
  !  FT 23.10.2020
  !
  !*********************************************


  USE id_base,         ONLY: idbase
  USE bns_lorene,      ONLY: bnslorene
  USE diffstar_lorene, ONLY: diffstarlorene


  IMPLICIT NONE


  CHARACTER( LEN= 5 ), PARAMETER:: bnslo= "BNSLO"
  !# String that identifies a binary system of neutron stars computed
  !  with LORENE
  CHARACTER( LEN= 5 ), PARAMETER:: drslo= "DRSLO"
  !# String that identifies a differentially rotating star computed
  !  with LORENE


  CONTAINS


  SUBROUTINE allocate_idbase( id, filename )

    !*********************************************
    !
    !# This SUBROUTINE allocates a polymorphic
    !  object of class idbase to its dynamic type.
    !  The dynamic type is one among those that
    !  use the |lorene| |id|
    !
    !  FT 9.11.2020
    !
    !*********************************************

    IMPLICIT NONE

    CLASS( idbase ), ALLOCATABLE, INTENT( IN OUT ):: id
    CHARACTER(LEN=*), INTENT( IN ) :: filename

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

    IF( filename(1:5) == bnslo )THEN

      ALLOCATE( bnslorene:: id )

    ELSEIF( filename(1:5) == drslo )THEN

      ALLOCATE( diffstarlorene:: id )

    ELSE

      PRINT *, "** ERROR! Unknown name for the physical system: ", filename(1:5)
      PRINT *, "   Set the variable 'system' in the parameter file ", &
               "sphincs_lorene_parameters.par to one of the values ", &
               "listed there."
      PRINT *, "   Stopping..."
      PRINT *
      STOP

    ENDIF

  END SUBROUTINE allocate_idbase


END MODULE sphincs_lorene
