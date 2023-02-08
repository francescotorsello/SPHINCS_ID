! File:         submodule_diffstar_base_io.f90
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

SUBMODULE (diffstar_base) io

  !***************************************************
  !
  !# This submodule contains the implementation of the
  !  methods of TYPE diffstarbase that handle I/O (input/output)
  !
  !  FT 5.11.2021
  !
  !***************************************************


  IMPLICIT NONE


  CONTAINS


  !------------------------------!
  !--  OVERRIDING SUBROUTINES  --!
  !------------------------------!


  MODULE PROCEDURE print_summary_drs

    !************************************************
    !
    !# Prints a summary of the physical properties the system
    !  to the standard output and, optionally, to a formatted
    !  file whose name is given as the optional argument `filename`
    !
    !  FT 5.11.2021
    !
    !************************************************

    USE utility,  ONLY: lorene2hydrobase, kg2g, m2cm

    IMPLICIT NONE

    PRINT *, " * Differentially rotating star (DRS):"
    PRINT *
    PRINT *, "   Baryon mass of the DRS=", this% mass, "Msun"
    PRINT *, "   Gravitational mass of the DRS=", this% mass_grav, "Msun"
    PRINT *
    PRINT *, "   Equatorial (not areal) radius of the DRS at phi=0 = ", &
                 this% r_eq, "Msun_geo"
    PRINT *, "   Equatorial (not areal) radius of the DRS at phi=pi/2 = ", &
                 this% r_eq_pi2, "Msun_geo"
    PRINT *, "   Equatorial (not areal) radius of the DRS at phi=pi = ", &
                 this% r_eq_pi, "Msun_geo"
    PRINT *, "   Equatorial (not areal) radius of the DRS at phi=3pi/2 = ", &
                 this% r_eq_3pi2, "Msun_geo"
    PRINT *, "   Polar radius of the DRS= ", this% r_pole, "Msun_geo"
    PRINT *, "   Ratio between polar radius at equatiorial radius at phi=0= ", &
                 this% r_ratio, "Msun_geo"
    PRINT *
    PRINT *, "   EOS for the DRS= ", this% eos
    PRINT *
    PRINT *, "   Central baryon mass density for the DRS= ", &
                 this% rho_center, "Msun/Msun_geo**3= ", &
                 this% rho_center/lorene2hydrobase*kg2g/(m2cm**3), "g cm^{-3}"
    PRINT *
    PRINT *, "   Ratio T/|W| between the rotaional kinetic energy and ", &
             "the gravitational binding energy: ", this% tsw
    PRINT *, "   For axisymmetric configurations as this one, the ", &
             "threshold for dynamical bar-mode instability is T/|W|~0.25 ", &
             " [Masaru Shibata et al 2000 ApJ 542 453, ", &
             "https://arxiv.org/pdf/astro-ph/0005378.pdf]. See also ", &
             "[Manca et al., Classical and Quantum Gravity, 24, 171, ", &
             "https://arxiv.org/abs/0705.1826], ", &
             "Sec.3.3 in [Galeazzi et al., Astron Astrophys 541:A156, ", &
             "arXiv:1101.2664], and Sec.5.1.3 in ", &
             "[Paschalidis, V., Stergioulas, N., Rotating stars in ", &
             "relativity. Living Rev Relativ 20, 7 (2017), ", &
             "https://link.springer.com/article/10.1007%2Fs41114-017-0008-x]."
    PRINT *


  END PROCEDURE print_summary_drs


END SUBMODULE io
