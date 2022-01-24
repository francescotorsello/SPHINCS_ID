! File:         submodule_diffstar_lorene_io.f90
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

SUBMODULE (diffstar_lorene) io

  !********************************************
  !
  !# This submodule contains the implementation of the
  !  methods of TYPE diffstarlorene that handle I/O (input/output)
  !
  !  FT 05.11.2021
  !
  !********************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE print_diffstar_params

    !****************************************************
    !
    !# Print the parameters of the binary neutron
    !  stars' initial data computed by |lorene|
    !
    !  FT 8.10.2020
    !
    !****************************************************

    USE constants, ONLY: k_lorene2hydrobase, Msun_geo, km2m, m2cm, kg2g, &
                         lorene2hydrobase, zero

    IMPLICIT NONE

    IF( THIS% angular_momentum == zero )THEN

      PRINT *
      PRINT *, " ** The parameters have not ben read yet. ", &
          "Call the SUBROUTINE import_diffstar_params to read them."
      PRINT *

    ELSE

      PRINT *
      PRINT *, " ** The parameters of the differentially rotating star are:"
      PRINT *
      PRINT *, " Baryonic mass = ", THIS% mass, " M_sun"
      PRINT *, " Gravitational mass = ", THIS% mass_grav, " M_sun"
      PRINT *, " Angular momentum = ", THIS% angular_momentum, " G M_sun^2 /c"
      PRINT *, " Surface area = ", THIS% surface_area, " M_sun^2", &
                                   THIS% surface_area*Msun_geo**2.0D0, " km^2"
      PRINT *
      PRINT *, " Radii: "
      PRINT *, "  Areal (or circumferential) radius for the star in the", &
               "  binary system [the one used in the", &
               "  (gravitational)mass-(areal)radius diagrams",&
               "  is for a TOV star], x direction:", &
               THIS% area_radius, " M_sun^geo = ", &
               THIS% area_radius*Msun_geo, " km", &
               THIS% r_circ, " M_sun^geo = ", &
               THIS% r_circ*Msun_geo, " km"
      PRINT *, "  Mean radius = ", THIS% r_mean, " M_sun^geo"
      PRINT *, "  Equatorial radius at phi=0 = ",  THIS% r_eq, " M_sun^geo"
      PRINT *, "  Equatorial radius at phi=pi/2 = ",  THIS% r_eq_pi2, &
               " M_sun^geo"
      PRINT *, "  Equatorial radius at phi=pi = ",  THIS% r_eq_pi, " M_sun^geo"
      PRINT *, "  Equatorial radius at phi=3pi/2 = ",  THIS% r_eq_3pi2, &
               " M_sun^geo"
      PRINT *, "  Polar radius = ",  THIS% r_pole, " M_sun^geo"
      PRINT *, "  Polar radius/(Equatorial radius at phi=0) = ", &
               THIS% r_ratio
      PRINT *, " Radius of the Innermost Stable Circular Orbit (ISCO) = ", &
               THIS% r_isco, " M_sun^geo"
      PRINT *, " Orbital frequency of the Innermost Stable Circular Orbit ", &
               "(ISCO) = ", THIS% f_isco, " M_sun^geo"
      PRINT *, " Specific energy of a test particle at the Innermost Stable ", &
               "Circular Orbit (ISCO) ", THIS% specific_energy_isco, &
               " c^2"
      PRINT *, " Specific angular momentum of a test particle at the ", &
               "Innermost Stable Circular Orbit (ISCO) ", &
               THIS% specific_angular_momentum_isco, " G M_sun /c"
      PRINT *
      PRINT *, " Hydro quantities at the center of the star: "
      PRINT *, "  Central enthalpy = ", THIS% ent_center, " c^2"
      PRINT *, "  Central baryon number density = ", THIS% nbar_center, &
               " (M_sun^geo)^{-3} =", &
               THIS% nbar_center/(MSun_geo*km2m*m2cm)**3, "cm^{-3}"
      PRINT *, "  Central baryon mass density = ", THIS% rho_center, &
               " M_sun^geo (M_sun^geo)^{-3} =", &
               THIS% rho_center/lorene2hydrobase*kg2g/(m2cm**3), "g cm^{-3}"
      PRINT *, "  Central energy density = ", THIS% energy_density_center, &
               " M_sun^geo c^2 (M_sun^geo)^{-3}", &
               THIS% energy_density_center/lorene2hydrobase*kg2g/(m2cm**3), &
               "g c^2 cm^{-3}"
      PRINT *, "  Central specific energy = ", THIS% specific_energy_center, &
               " c^2"
      PRINT *, "  Central pressure = ", THIS% pressure_center, &
               " M_sun^geo c^2 (M_sun^geo)^{-3}", &
               THIS% pressure_center/lorene2hydrobase*kg2g/(m2cm**3), &
               "g c^2 cm^{-3}"
      PRINT *
      PRINT *, " Ratio T/|W| between the rotational kinetic energy and ", &
               "the gravitational binding energy: ", THIS% tsw
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

      PRINT *, " Equations of state for star 1 (EOS1) = ", TRIM(THIS% eos)

      IF( THIS% eos_loreneid == 1 )THEN ! If the EOS is polytropic

        PRINT *, " Parameters for EOS: "
        PRINT *, "  Polytopic index gamma = ", THIS% gamma
        PRINT *, "  Pressure coefficient = ",&
                 THIS% kappa/k_lorene2hydrobase( THIS% gamma ), &
                 "rho_nuc c^2 / n_nuc^gamma = ", THIS% kappa, &
                 "[pure number]"
        PRINT *

      ELSEIF( THIS% gamma0 /= 0 )THEN ! If the EOS is piecewise polytropic

        PRINT *, " Parameters for EOS1: "
        PRINT *, "  Number of polytropic indexes = ", THIS% npeos
        PRINT *, "  Polytopic index gamma0 = ", THIS% gamma0
        PRINT *, "  Polytopic index gamma1 = ", THIS% gamma1
        PRINT *, "  Polytopic index gamma2 = ", THIS% gamma2
        PRINT *, "  Polytopic index gamma3 = ", THIS% gamma3
        PRINT *, "  Pressure coefficient for the crust (here from SLy) = ",&
                 THIS% kappa0/k_lorene2hydrobase( THIS% gamma0 ), &
                 "rho_nuc c^2 / n_nuc^gamma0 = ", THIS% kappa0, &
                 "[pure number]"
        PRINT *, "  Pressure coefficient for the first polytrope = ",&
                 THIS% kappa1/k_lorene2hydrobase( THIS% gamma1 ), &
                 "rho_nuc c^2 / n_nuc^gamma1", THIS% kappa1, &
                 "[pure number]"
        PRINT *, "  Pressure coefficient for the second polytrope = ",&
                 THIS% kappa2/k_lorene2hydrobase( THIS% gamma2 ), &
                 "rho_nuc c^2 / n_nuc^gamma2", THIS% kappa2, &
                 "[pure number]"
        PRINT *, "  Pressure coefficient for the third polytrope = ",&
                 THIS% kappa3/k_lorene2hydrobase( THIS% gamma3 ), &
                 "rho_nuc c^2 / n_nuc^gamma3", THIS% kappa3, &
                 "[pure number]"
        PRINT *, "  Base 10 exponent of the pressure at the first fiducial " &
                 // "density (between gamma_0 and gamma_1) (dyne/cm^2)= ", &
                 THIS% logP1
        PRINT *, "  Base 10 exponent of first fiducial density (g/cm^3) = ", &
                 THIS% logRho0
        PRINT *, "  Base 10 exponent of second fiducial density (g/cm^3) = ",&
                 THIS% logRho1
        PRINT *, "  Base 10 exponent of third fiducial density (g/cm^3) = ", &
                 THIS% logRho2
        PRINT *

      ELSEIF( THIS% eos_loreneid == 17 .OR. THIS% eos_loreneid == 20 )THEN
      ! If the EOS is tabulated

      ELSE

        PRINT *, "** ERROR in SUBROUTINE import_lorene_id_params!", &
                 " The equation of state is unknown!"
        STOP

      ENDIF

    ENDIF

  END PROCEDURE print_diffstar_params


END SUBMODULE io
