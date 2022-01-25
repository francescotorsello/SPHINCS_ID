! File:         submodule_bns_base_io.f90
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

SUBMODULE (bns_base) io

  !***************************************************
  !
  !# This submodule contains the implementation of the
  !  methods of TYPE bnsbase that handle I/O (input/output)
  !
  !  FT 5.11.2021
  !
  !***************************************************


  IMPLICIT NONE


  CONTAINS


  !------------------------------!
  !--  OVERRIDING SUBROUTINES  --!
  !------------------------------!


  MODULE PROCEDURE print_summary_bns

    !************************************************
    !
    !# Prints a summary of the physical properties the system
    !  to the standard output and, optionally, to a formatted
    !  file whose name is given as the optional argument `filename`
    !
    !  FT 5.11.2021
    !
    !************************************************

    USE constants,      ONLY: lorene2hydrobase, kg2g, m2cm

    IMPLICIT NONE

    PRINT *, " * Binary system of neutron stars:"
    PRINT *
    PRINT *, "   Baryon mass of star 1=", THIS% get_mass1(), "Msun"
    PRINT *, "   Baryon mass of star 2=", THIS% get_mass2(), "Msun"
    PRINT *, "   Gravitational mass of star 1=", THIS% get_grav_mass1(), "Msun"
    PRINT *, "   Gravitational mass of star 2=", THIS% get_grav_mass2(), "Msun"
    PRINT *
    PRINT *, "   Equatorial (not areal) radius of neutron star 1 towards " &
             // "companion= ", &
                 THIS% get_radius1_x_comp(), "Msun_geo"
    PRINT *, "   Equatorial (not areal) radius of neutron star 2 towards " &
             // "companion= ", &
                 THIS% get_radius2_x_comp(), "Msun_geo"
    PRINT *, "   Equatorial (not areal) radius of neutron star 1 opposite to " &
             // "companion= ", &
                 THIS% get_radius1_x_opp(), "Msun_geo"
    PRINT *, "   Equatorial (not areal) radius of neutron star 2 opposite to " &
             // "companion= ", &
                 THIS% get_radius2_x_opp(), "Msun_geo"
    PRINT *, "   Radius (not areal) along y of neutron star 1= ", &
                 THIS% get_radius1_y(), "Msun_geo"
    PRINT *, "   Radius (not areal) along y of neutron star 2= ", &
                 THIS% get_radius2_y(), "Msun_geo"
    PRINT *, "   Radius (not areal) along y of neutron star 1= ", &
                 THIS% get_radius1_z(), "Msun_geo"
    PRINT *, "   Radius (not areal) along y of neutron star 2= ", &
                 THIS% get_radius2_z(), "Msun_geo"
    PRINT *
    PRINT *, "   EOS for neutron star 1= ", &
                 THIS% get_eos1()
    PRINT *, "   EOS for neutron star 2= ", &
                 THIS% get_eos2()
    PRINT *
    PRINT *, "   Central baryon mass density for star 1= ", &
                 THIS% get_rho_center1(), "Msun/Msun_geo**3= ", &
                 THIS% get_rho_center1() &
                 /lorene2hydrobase*kg2g/(m2cm**3), "g cm^{-3}"
    PRINT *, "   Central baryon mass density for star 2= ", &
                 THIS% get_rho_center2(), "Msun/Msun_geo**3= ", &
                 THIS% get_rho_center2() &
                 /lorene2hydrobase*kg2g/(m2cm**3), "g cm^{-3}"
    PRINT *


  END PROCEDURE print_summary_bns


END SUBMODULE io
