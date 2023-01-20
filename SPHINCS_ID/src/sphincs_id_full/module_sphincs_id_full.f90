! File:         module_sphincs_id_full.f90
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

MODULE sphincs_id_full


  !*********************************************
  !
  !# This module contains data and PROCEDURES
  !  needed to set up all the supported
  !  |id| in PROGRAM [[sphincs_id]]
  !
  !  FT 09.02.2022
  !
  !*********************************************


  USE id_base,         ONLY: idbase
  USE bns_lorene,      ONLY: bnslorene
  USE diffstar_lorene, ONLY: diffstarlorene
  USE bns_fuka,        ONLY: bnsfuka
  USE ejecta_generic,  ONLY: ejecta


  IMPLICIT NONE


  CHARACTER( LEN= 5 ), PARAMETER:: bnslo= "BNSLO"
  !# String that identifies a binary system of neutron stars computed
  !  with |lorene|
  CHARACTER( LEN= 5 ), PARAMETER:: drslo= "DRSLO"
  !# String that identifies a differentially rotating star computed
  !  with |lorene|
  CHARACTER( LEN= 5 ), PARAMETER:: bnsfu= "BNSFU"
  !# String that identifies a binary system of neutron stars computed
  !  with |fuka|
  CHARACTER( LEN= 5 ), PARAMETER:: ejecta_grid= "EJECT"
  !# String that identifies an ejecta prepared on a uniform Cartesian grid


  CONTAINS


  SUBROUTINE allocate_idbase( id, filename, system, system_name )

    !*********************************************
    !
    !# This SUBROUTINE allocates a polymorphic
    !  object of class idbase to its dynamic type.
    !  The dynamic type is one among  all the
    !  supported |id|
    !
    !  FT 9.11.2020
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

    IF( filename(1:5) == bnslo )THEN

      ALLOCATE( bnslorene:: id )
      system= bnslo
      system_name= "NSNS."

    ELSEIF( filename(1:5) == drslo )THEN

      ALLOCATE( diffstarlorene:: id )
      system= drslo
      system_name= "DRSx."

    ELSEIF( filename(1:5) == bnsfu )THEN

      ALLOCATE( bnsfuka:: id )
      system= bnsfu
      system_name= "NSNS."

    ELSEIF( filename(1:5) == ejecta_grid )THEN

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
              " supported by this version of SPHINCS_ID, with flavour = 1, are:"
      PRINT *
      PRINT *, "   BNSLO: Binary Neutron Stars produced with LORENE"
      PRINT *, "   DRSLO: Differentially Rotating Star produced with LORENE"
      PRINT *, "   BNSFU: Binary Neutron Stars produced with FUKA"
      PRINT *, "   EJECT: Ejecta data on a uniform Cartesian grid"
      PRINT *
      STOP

    ENDIF

  END SUBROUTINE allocate_idbase


  FUNCTION shorten_eos_name( eos_long ) RESULT( eos_str )

    !*********************************************
    !                                            *
    ! Converts the longer names from |lorene|    *
    ! to the 4-character strings needed by       *
    ! select_EOS_parameters                      *
    !                                            *
    ! FT 11.08.2021                              *
    !                                            *
    !*********************************************

    IMPLICIT NONE

    CHARACTER(4):: eos_str
    CHARACTER(LEN=*), INTENT(IN):: eos_long

    SELECT CASE( eos_long )

      CASE( 'Multipolytropic SLy  EOS' )

        eos_str= 'SLy '

      CASE( 'Multipolytropic ALF2 EOS' )

        eos_str= 'ALF2'

      CASE( 'Multipolytropic ALF4 EOS' )

        eos_str= 'ALF4'

      CASE( 'Multipolytropic ENG  EOS' )

        eos_str= 'ENG '

      CASE( 'Multipolytropic H4  EOS' )

        eos_str= 'H4  '

      CASE( 'Multipolytropic MPA1 EOS' )

        eos_str= 'MPA1'

      CASE( 'Star 1  Multipolytropic MPA1 EOS' )

        eos_str= 'MPA1'

      CASE( 'Multipolytropic MS1  EOS' )

        eos_str= 'MS1 '

      CASE( 'Multipolytropic MS1b EOS' )

        eos_str= 'MS1b'

      CASE( 'Multipolytropic APR3 EOS' )

        eos_str= 'AP3 '

      CASE( 'Star 1  Multipolytropic APR3 EOS' )

        eos_str= 'AP3 '

      CASE( 'Multipolytropic AP3  EOS' )

        eos_str= 'AP3 '

      CASE( 'Multipolytropic WFF1 EOS' )

        eos_str= 'WFF1'

      CASE( 'Multipolytropic WFF2 EOS' )

        eos_str= 'WFF2'

      CASE( 'Multipolytropic GNH3 EOS' )

        eos_str= 'GNH3'

      CASE( 'Multipolytropic APR4 EOS' )

        eos_str= 'APR4'

      CASE( 'Multipolytropic AP4  EOS' )

        eos_str= 'APR4'

      CASE( 'Multipolytropic haso EOS' )

        eos_str= 'haso'

      CASE( 'Polytropic EOS' )

        eos_str= '   '

      CASE DEFAULT

        PRINT *, "** ERROR! Unknown EOS name: ", eos_long
        PRINT *, " * Please add the name to SUBROUTINE shorten_eos_name. ", &
                 "   This SUBROUTINE exists mostly for backward compatibility."
        PRINT *, " * Stopping..."
        PRINT *
        STOP

    END SELECT

  END FUNCTION shorten_eos_name


  FUNCTION shorten_eos_name_fuka( eos_long ) RESULT( eos_str )

    !*********************************************
    !                                            *
    ! Converts the longer names from |fuka|      *
    ! to the 4-character strings needed by       *
    ! select_EOS_parameters                      *
    !                                            *
    ! FT 26.08.2022                              *
    !                                            *
    !*********************************************

    IMPLICIT NONE

    CHARACTER(4):: eos_str
    CHARACTER(LEN=*), INTENT(IN):: eos_long

    SELECT CASE( TRIM(eos_long) )

      CASE( 'sly4.pwpoly' )

        eos_str= 'SLy '

      CASE( 'alf2.pwpoly' )

        eos_str= 'ALF2'

      CASE( 'alf4.pwpoly' )

        eos_str= 'ALF4'

      CASE( 'eng.pwpoly' )

        eos_str= 'ENG '

      CASE( 'h4.pwpoly' )

        eos_str= 'H4  '

      CASE( 'mpa1.pwpoly' )

        eos_str= 'MPA1'

      CASE( 'ms1.pwpoly' )

        eos_str= 'MS1 '

      CASE( 'ms1b.pwpoly' )

        eos_str= 'MS1b'

      CASE( 'apr3.pwpoly' )

        eos_str= 'AP3 '

      CASE( 'ap3.pwpoly' )

        eos_str= 'AP3 '

      CASE( 'wff1.pwpoly' )

        eos_str= 'WFF1'

      CASE( 'wff2.pwpoly' )

        eos_str= 'WFF2'

      CASE( 'gnh3.pwpoly' )

        eos_str= 'GNH3'

      CASE( 'apr4.pwpoly' )

        eos_str= 'APR4'

      CASE( 'ap4.pwpoly' )

        eos_str= 'APR4'

      CASE( 'apr4.polytrope' )

        eos_str= 'APR4'

      CASE( 'haso.pwpoly' )

        eos_str= 'haso'

      CASE( 'gam2.polytr' )

        eos_str= '   '

      CASE( 'gam2.polytrope' )

        eos_str= '   '

      CASE DEFAULT

        PRINT *, "** ERROR! Unknown EOS name: ", TRIM(eos_long)
        PRINT *, " * Please add the name to SUBROUTINE shorten_eos_name_fuka. "
        PRINT *, " * Stopping..."
        PRINT *
        STOP

    END SELECT

  END FUNCTION shorten_eos_name_fuka


END MODULE sphincs_id_full
