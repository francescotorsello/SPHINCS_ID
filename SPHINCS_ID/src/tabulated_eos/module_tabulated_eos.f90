! File:         module_tabulated_eos.f90
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

MODULE tabulated_eos


  !***********************************************************
  !
  !# This module contains the PROCEDURES needed to deal
  !  with tabulated |eos|.
  !
  !  FT 10.03.2023
  !
  !***********************************************************


  USE utility,    ONLY: ios, err_msg, km2Msun_geo, MeV2amuc2
  USE constants,  ONLY: fm2cm, cm2km, MeV2erg, amu, Msun, c_light2_si
  !USE units, ONLY: m0c2_cu, set_units ! needed for debugging


  IMPLICIT NONE


  DOUBLE PRECISION, PARAMETER:: fm2Msun_geo= fm2cm*cm2km*km2Msun_geo


  CONTAINS


  SUBROUTINE read_compose_beta_equilibrated_eos(namefile, table_eos)

    !**********************************************
    !
    !# Read the |compose| files *.thermo.ns,
    !  *.t, *.nb, *.yq,
    !  containing the \(\beta\)-equilibrated |eos|
    !  computed by the |compose| software.
    !
    !  See [the |compose| Manual](https://compose.obspm.fr/manual){:target="_blank"} for more details.
    !
    !  FT 10.03.2023
    !
    !**********************************************


    IMPLICIT NONE


    CHARACTER(LEN=*), INTENT(IN):: namefile
    !! Path of the |eos| file without '.thermo.ns' extension
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT):: table_eos
    !#  The |eos| table in \compose| format contains:
    !
    ! - First table index; temperature
    !
    ! - Second table index; baryon number density
    !
    ! - Third table index; charge fraction of strongly interacting particles
    !
    ! - Pressure divided by the baryon number density:
    !   \(\dfrac{p}{|nb|}\mathrm{(MeV)}\). After reading it,
    !   it is converted into pressure \(M_\odot c^2 L_\odot^{-3}\)
    !
    ! - Entropy density divided by the baryon number density
    !   (entropy per baryon): \(\dfrac{s}{|nb|} \mathrm{(dimensionless)}\)
    !
    ! - Scaled and shifted baryon chemical potential:
    !   \(\dfrac{\mu_\mathrm{b}}{m_\mathrm{n}} - 1 \mathrm{(dimensionless)}\)
    !
    ! - Scaled charge chemical potential:
    !   \(\dfrac{\mu_\mathrm{q}}{m_\mathrm{n}} \mathrm{(dimensionless)}\)
    !
    ! - Scaled effective lepton chemical potential:
    !   \(\dfrac{\mu_\mathrm{l}}{m_\mathrm{n}} \mathrm{(dimensionless)}\)
    !
    ! - Scaled free energy per baryon:
    !   \(\dfrac{f}{|nb|m_\mathrm{n}} - 1 \mathrm{(dimensionless)}\)
    !
    ! - Scaled internal energy per baryon:
    !   \(\dfrac{e}{|nb|m_\mathrm{n}} \mathrm{(dimensionless)}\)
    !
    ! - Number of additional quantities (this should be 2, as per CompOSE
    !   Manual v3.00, 13.03.2023)
    !
    ! - Electron fraction in beta equilibrium
    !
    ! - Enthalpy density without the temperature term \(Ts\):
    !   \(h= e + p \mathrm{MeV}\,\mathrm{fm}^{-3}\)
    !
    ! - Baryon number density \((\mathrm{fm}^{-3})\). After reading it,
    !   it is converted in baryon mass density \(M_\odot L_\odot^{-3}\)
    !

    INTEGER, PARAMETER:: unit_compose= 876
    INTEGER, PARAMETER:: max_length_table= 10000

    INTEGER:: itr, n_lines, tmp

    !
    !-- Data contained in the *.thermo.ns file
    !
    DOUBLE PRECISION:: m_n
    !! Mass of the neutron \(m_\mathrm{n} \mathrm{(MeV)}\)
    DOUBLE PRECISION:: m_p
    !! Mass of the proton \(m_\mathrm{p} \mathrm{(MeV)}\)
    INTEGER:: i_leptons
    !# Index that says if the |eos| considers leptons.
    !
    !  - `i_leptons`= 1 if leptons are considered
    !  - `i_leptons`= anything else if leptons are not considered
!    INTEGER:: i_t
!    !! First table index; temperature
!    INTEGER:: i_nb
!    !! Second table index; baryon number density
!    INTEGER:: i_q
!    !! Third table index; charge fraction of strongly interacting particles
!    DOUBLE PRECISION:: p_nb
!    !# Pressure divided by the baryon number density:
!    !  \(\dfrac{p}{|nb|}\mathrm{(MeV)}\)
!    DOUBLE PRECISION:: s_nb
!    !# Entropy density divided by the baryon number density
!    !  (entropy per baryon): \(\dfrac{s}{|nb|} \mathrm{(dimensionless)}\)
!    DOUBLE PRECISION:: mub_mn
!    !# Scaled and shifted baryon chemical potential:
!    !  \(\dfrac{\mu_\mathrm{b}}{m_\mathrm{n}} - 1 \mathrm{(dimensionless)}\)
!    DOUBLE PRECISION:: muq_mn
!    !# Scaled charge chemical potential:
!    !  \(\dfrac{\mu_\mathrm{q}}{m_\mathrm{n}} \mathrm{(dimensionless)}\)
!    DOUBLE PRECISION:: mul_mn
!    !# Scaled effective lepton chemical potential:
!    !  \(\dfrac{\mu_\mathrm{l}}{m_\mathrm{n}} \mathrm{(dimensionless)}\)
!    DOUBLE PRECISION:: f_nbmn
!    !# Scaled free energy per baryon:
!    !  \(\dfrac{f}{|nb|m_\mathrm{n}} - 1 \mathrm{(dimensionless)}\)
!    DOUBLE PRECISION:: e_nbmn
!    !# Scaled internal energy per baryon:
!    !  \(\dfrac{e}{|nb|m_\mathrm{n}} \mathrm{(dimensionless)}\)
!    INTEGER:: n_add
!    !# Number of additional quantities
!    DOUBLE PRECISION:: Y_e
!    !# Electron fraction in beta equilibrium
!    DOUBLE PRECISION:: h
!    !# Enthalpy density without the temperature term \(Ts\):
!    !  \(h= e + p \mathrm{MeV}\,\mathrm{fm}^{-3}\)

    LOGICAL:: exist

    CHARACTER(LEN=:), ALLOCATABLE:: finalnamefile


    !PRINT *
    !PRINT *, "** Executing the read_compose_eos subroutine..."
    PRINT *

    ALLOCATE( table_eos(14,max_length_table) )
    table_eos= 0.D0

    finalnamefile= TRIM(namefile)//"eos.thermo.ns"

    INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

    IF( exist )THEN
      OPEN( UNIT= unit_compose, FILE= TRIM(finalnamefile), &
            FORM= "FORMATTED", ACTION= "READ", IOSTAT= ios, &
            IOMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...error when opening " // TRIM(finalnamefile), &
                ". The error message is", err_msg
        PRINT *
        STOP
      ENDIF
    ELSE
      PRINT *, "** ERROR! Unable to find file " // TRIM(finalnamefile)
      PRINT *
      STOP
    ENDIF

    READ( UNIT= unit_compose, FMT= *, IOSTAT = ios, IOMSG= err_msg ) &
                                                        m_n, m_p, i_leptons

    PRINT *, " * Reading file " // TRIM(finalnamefile) // "..."
    n_lines= 0
    read_compose_table: DO itr= 1, max_length_table, 1

      READ( UNIT= unit_compose, FMT= *, IOSTAT = ios, IOMSG= err_msg ) &
                        ! Variables for .thermo.ns file
                        table_eos(1,itr),  table_eos(2,itr), table_eos(3,itr), &
                        table_eos(4,itr),  table_eos(5,itr), table_eos(6,itr), &
                        table_eos(7,itr),  table_eos(8,itr), table_eos(9,itr), &
                        table_eos(10,itr), table_eos(11,itr), &
                        table_eos(12,itr), table_eos(13,itr)
      IF( ios > 0 )THEN
        PRINT *, "...error when reading " // TRIM(finalnamefile), &
                ". The error message is", err_msg
        PRINT *
        STOP
      ENDIF
      IF( ios < 0 )THEN
        PRINT *, " * Reached end of file " // TRIM(finalnamefile)
        PRINT *
        EXIT read_compose_table
      ENDIF
      n_lines= n_lines + 1

    ENDDO read_compose_table

    ! Reallocate arrays to delete unnecessary elements
    table_eos= table_eos(:,1:n_lines)

    CLOSE( unit_compose )

    ! Read baryon number density from *.nb.ns file
    finalnamefile= TRIM(namefile)//"eos.nb.ns"

    INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

    IF( exist )THEN
      OPEN( UNIT= unit_compose, FILE= TRIM(finalnamefile), &
            FORM= "FORMATTED", ACTION= "READ", IOSTAT= ios, &
            IOMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...error when opening " // TRIM(finalnamefile), &
                ". The error message is", err_msg
        PRINT *
        STOP
      ENDIF
    ELSE
      PRINT *, "** ERROR! Unable to find file " // TRIM(finalnamefile)
      PRINT *
      STOP
    ENDIF

    READ( UNIT= unit_compose, FMT= *, IOSTAT= ios, IOMSG= err_msg ) tmp, n_lines

    PRINT *, " * Reading file " // TRIM(finalnamefile) // "..."

    read_compose_nb: DO itr= 1, n_lines, 1
      READ( UNIT= unit_compose, FMT= *, IOSTAT = ios, IOMSG= err_msg ) &
                        ! Variable for .nb file
                        ! (baryon number density, independent variable)
                        table_eos(14,itr)
      IF( ios > 0 )THEN
        PRINT *, "...error when reading " // TRIM(finalnamefile), &
                ". The error message is", err_msg
        PRINT *
        STOP
      ENDIF

    ENDDO read_compose_nb

    CLOSE( unit_compose )

    PRINT *, "...done."
    PRINT *

    !
    !-- Convert (some of the) physical quantities to desired quantities
    !-- in the desiredunits
    !

    ! Pressure: MeV fm^{-3} to Msun_geo c^2 Msun_geo^{-3}
    ! (remember that c=1, so we do not have to divide by it)
    table_eos(4,:)= table_eos(4,:)*table_eos(14,:) &
                    *(MeV2amuc2*amu/Msun)/(fm2Msun_geo**3)

    ! Baryon mass density: MeV fm^{-3} to Msun_geo Msun_geo^{-3}
    ! (remember that c=1, so we do not have to divide by it)
    table_eos(14,:)= m_n*table_eos(14,:) &
                     *(MeV2amuc2*amu/Msun)/(fm2Msun_geo**3)

    !CALL set_units('NSM')
    !
    !PRINT *, "n_lines:", n_lines
    !PRINT *, "SIZE(pressure):", SIZE(table_eos(4,:))
    !PRINT *, "SIZE(baryon mass density):", SIZE(table_eos(14,:))
    !PRINT *, "SIZE(specific internal energy):", SIZE(table_eos(10,:))
    !PRINT *, "pressure:", table_eos(4,210)/lorene2hydrobase
    !PRINT *, "baryon mass density:", table_eos(14,210)/lorene2hydrobase
    !PRINT *, "specific internal energy:", table_eos(10,210)
    !PRINT *, "m0c2_cu:", m0c2_cu
    !PRINT *, "pressure*m0c2_cu:", table_eos(4,210)/m0c2_cu
    !PRINT *, "baryon mass density*m0c2_cu:", table_eos(14,210)/m0c2_cu
    !PRINT *
    !STOP


    !PRINT *, "** Subroutine read_compose_eos executed."
    !PRINT *


  END SUBROUTINE read_compose_beta_equilibrated_eos


END MODULE tabulated_eos


!    ! Read charge fraction of strongly interacting particles from *.t file
!    ! (we don't have charge fraction of strongly interacting particles for the
!    ! beta equilibrated data)
!    finalnamefile= TRIM(namefile)//".yq"
!
!    INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )
!
!    IF( exist )THEN
!      OPEN( UNIT= unit_compose, FILE= TRIM(finalnamefile), &
!            FORM= "FORMATTED", ACTION= "READ", IOSTAT= ios, &
!            IOMSG= err_msg )
!      IF( ios > 0 )THEN
!        PRINT *, "...error when opening " // TRIM(finalnamefile), &
!                ". The error message is", err_msg
!        PRINT *
!        STOP
!      ENDIF
!    ELSE
!      PRINT *, "** ERROR! Unable to find file " // TRIM(finalnamefile)
!      PRINT *
!      STOP
!    ENDIF
!
!    READ( UNIT= unit_compose, FMT= *, IOSTAT= ios, IOMSG= err_msg )tmp, n_lines
!
!    PRINT *, " * Reading file " // TRIM(finalnamefile) // "..."
!
!    read_compose_nb: DO itr= 1, n_lines, 1
!      READ( UNIT= unit_compose, FMT= *, IOSTAT = ios, IOMSG= err_msg ) &
!                        ! Variable for .nb file
!                        ! (charge fraction of strongly
!                        ! interacting particles, independent variable)
!                        table_eos(6,itr)
!      IF( ios > 0 )THEN
!        PRINT *, "...error when reading " // TRIM(finalnamefile), &
!                ". The error message is", err_msg
!        PRINT *
!        STOP
!      ENDIF
!
!    ENDDO read_compose_nb
!
!    CLOSE( unit_compose )


!    ! Read temperature from *.t file (we don't have temperature for the
!    ! beta equilibrated data)
!    finalnamefile= TRIM(namefile)//".t"
!
!    INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )
!
!    IF( exist )THEN
!      OPEN( UNIT= unit_compose, FILE= TRIM(finalnamefile), &
!            FORM= "FORMATTED", ACTION= "READ", IOSTAT= ios, &
!            IOMSG= err_msg )
!      IF( ios > 0 )THEN
!        PRINT *, "...error when opening " // TRIM(finalnamefile), &
!                ". The error message is", err_msg
!        PRINT *
!        STOP
!      ENDIF
!    ELSE
!      PRINT *, "** ERROR! Unable to find file " // TRIM(finalnamefile)
!      PRINT *
!      STOP
!    ENDIF
!
!    READ( UNIT= unit_compose, FMT= *, IOSTAT= ios, IOMSG= err_msg )tmp, n_lines
!
!    PRINT *, " * Reading file " // TRIM(finalnamefile) // "..."
!
!    read_compose_t: DO itr= 1, n_lines, 1
!      READ( UNIT= unit_compose, FMT= *, IOSTAT = ios, IOMSG= err_msg ) &
!                        ! Variable for .t file (temperature, independent variable)
!                        table_eos(4,itr)
!      IF( ios > 0 )THEN
!        PRINT *, "...error when reading " // TRIM(finalnamefile), &
!                ". The error message is", err_msg
!        PRINT *
!        STOP
!      ENDIF
!
!    ENDDO read_compose_t
!
!    CLOSE( unit_compose )
