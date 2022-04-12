! File:         submodule_bns_lorene_io.f90
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

SUBMODULE (bns_lorene) io

  !********************************************
  !
  !# This submodule contains the implementation of the
  !  methods of TYPE bnslorene that handle I/O (input/output)
  !
  !  FT 05.11.2021
  !
  !********************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE print_summary_bnslorene

    !************************************************
    !
    !# Prints a summary of the physical properties the |bns| system
    !  produced by | lorene to the standard output and, optionally,
    !  to a formatted file whose name is given as the optional
    !  argument `filename`
    !
    !  FT 4.02.2022
    !
    !************************************************

    IMPLICIT NONE

    PRINT *, "   * Binary system of neutron stars produced by LORENE:"
    PRINT *
    PRINT *, "     x coordinate of the center of mass of the system= ", &
             !"weighted with the baryonic mass= ", &
            (this% barycenter1_x*this% mass1 + this% barycenter2_x*this% mass2) &
             /(this% mass1 + this% mass2)
    !PRINT *, "   Center of mass of the system, weighted with the ", &
    !         "gravitational mass= ", (this% barycenter1_x*this% mass_grav1 &
    !         + this% barycenter2_x*this% mass_grav2) &
    !         /(this% mass_grav1 + this% mass_grav2)
    PRINT *
    PRINT *, "     ADM mass of the system= ", this% adm_mass, "MSun"
    PRINT *
    PRINT *, "     ADM linear momentum of the system=(", &
             this% linear_momentum_x, ", "
    PRINT *, "                                        ", &
             this% linear_momentum_y, ", "
    PRINT *, "                                        ", &
             this% linear_momentum_z, ") Msun*c"
    PRINT *
    PRINT *, "     Velocity of the center of mass of the system="
    PRINT *, "     ADM linear momentum / ADM mass =(", &
             this% linear_momentum_x/this% adm_mass, ", "
    PRINT *, "                                      ", &
             this% linear_momentum_y/this% adm_mass, ", "
    PRINT *, "                                      ", &
             this% linear_momentum_z/this% adm_mass, ") c"
    PRINT *
    PRINT *, "     Bowen-York angular momentum of the system= (", &
             this% angular_momentum_x, ", "
    PRINT *, "                                                 ", &
             this% angular_momentum_y, ", "
    PRINT *, "                                                 ", &
             this% angular_momentum_z, ") G*Msun^2/c"
    PRINT *


  END PROCEDURE print_summary_bnslorene


  MODULE PROCEDURE print_id_params

    !****************************************************
    !
    !# Print the parameters of the binary neutron
    !  stars' initial data computed by |lorene|
    !
    !  FT 8.10.2020
    !
    !****************************************************

    USE constants, ONLY: k_lorene2hydrobase, Msun_geo, km2m, m2cm, kg2g, &
                         lorene2hydrobase, c_light, cm2km

    IMPLICIT NONE

    IF( this% angular_momentum_z == 0.0D0 )THEN

      PRINT *
      PRINT *, " ** The parameters have not ben read yet. ", &
          "Call the SUBROUTINE import_lorene_id_params to read them."
      PRINT *

    ELSE

      PRINT *
      PRINT *, " ** The parameters of the binary system are:"
      PRINT *
      PRINT *, " Distance between the points of highest density = ",&
               this% distance, " M_sun^geo = ", this% distance*Msun_geo, " km"
      PRINT *, " Distance between the barycenters = ", &
               this% distance_com, " M_sun^geo", this% distance_com*Msun_geo, &
               " km"
      PRINT *
      PRINT *, " Baryonic mass of NS 1 = ", this% mass1, " M_sun"
      PRINT *, " Baryonic mass of NS 2 = ", this% mass2, " M_sun"
      PRINT *, " Gravitational mass of NS 1 = ", this% mass_grav1, " M_sun"
      PRINT *, " Gravitational mass of NS 2 = ", this% mass_grav2, " M_sun"
      PRINT *, " ADM mass = ", this% adm_mass, " M_sun"
      PRINT *
      PRINT *, " Stellar center of NS 1 = ", this% center1_x, " M_sun^geo"
      PRINT *, " Stellar center of NS 2 = ", this% center2_x, " M_sun^geo"
      PRINT *, " Barycenter of NS 1 = ", this% barycenter1_x, " M_sun^geo"
      PRINT *, " Barycenter of NS 2 = ", this% barycenter2_x, " M_sun^geo"
      PRINT *, " Orbital angular velocity Omega_0 = ", &
               this% angular_vel, " rad/s = ", &
               this% angular_vel/(c_light*cm2km), "km^{-1}"
      PRINT *, " mOmega = ", &
               "Omega_0[km^{-1}]*(mass_grav1[km] + mass_grav2[km]) = ",&
               this% mOmega, "[pure number]"
      PRINT *, " Bowen-York angular momentum of the system, x component = ", &
               this% angular_momentum_x, " G M_sun^2 /c"
      PRINT *, " Bowen-York angular momentum of the system, y component = ", &
               this% angular_momentum_y, " G M_sun^2 /c"
      PRINT *, " Bowen-York angular momentum of the system, z component = ", &
               this% angular_momentum_z, " G M_sun^2 /c"
      PRINT *, " Estimated time of the merger t_merger = ", this% t_merger, &
               " M_sun^geo = ", this% t_merger*MSun_geo/(c_light*cm2km)*1000.0,&
               " ms, from Peters_PR_136_B1224_1964, eq. (5.10)"
      PRINT *
      PRINT *, " Estimated separation to have the merger at t_merger = 2000", &
               " Msun_geo = ", 2000.0D0*MSun_geo/(c_light*cm2km)*1000.0, &
               " ms :", &
               ( 2000.0D0*( this% mass_grav1*this% mass_grav2* &
                  ( this% mass_grav1 + this% mass_grav2 ) )/(5.0D0/256.0D0) ) &
                **(1.0D0/4.0D0), "M_sun^geo = ", &
                ( 2000.0D0*( this% mass_grav1*this% mass_grav2* &
                ( this% mass_grav1 + this% mass_grav2 ) )/(5.0D0/256.0D0) ) &
                **(1.0D0/4.0D0)*Msun_geo, &
                "km, from Peters_PR_136_B1224_1964, eq. (5.10)"
      PRINT *
      PRINT *, " Radii of star 1: "
      PRINT *, "  Areal (or circumferential) radius for the star in the", &
               "  binary system [the one used in the", &
               "  (gravitational)mass-(areal)radius diagrams",&
               "  is for a TOV star], x direction:", &
               this% area_radius1, " M_sun^geo = ", &
               this% area_radius1*Msun_geo, " km"
      PRINT *, "  x direction, towards companion = ", &
               this% radius1_x_comp, " M_sun^geo"
      PRINT *, "  x direction, opposite to companion = ", &
               this% radius1_x_opp, " M_sun^geo"
      PRINT *, "  y direction = ", this% radius1_y, " M_sun^geo"
      PRINT *, "  z direction = ", this% radius1_z, " M_sun^geo"
      PRINT *, " Radii of star 2 :"
      PRINT *, "  Areal (or circumferential) radius for the star in the", &
               "  binary system [the one used in the", &
               "  (gravitational)mass-(areal)radius diagrams",&
               "  is for a TOV star], x direction:", &
               this% area_radius2, " M_sun^geo", &
               this% area_radius2*Msun_geo, " km"
      PRINT *, "  x direction, towards companion = ", &
               this% radius2_x_comp, " M_sun^geo"
      PRINT *, "  x direction, opposite to companion = ", &
               this% radius2_x_opp, " M_sun^geo"
      PRINT *, "  y direction = ", this% radius2_y, " M_sun^geo"
      PRINT *, "  z direction = ", this% radius2_z, " M_sun^geo"
      PRINT *
      PRINT *, " Hydro quantities at the center of star 1: "
      PRINT *, "  Central enthalpy = ", this% ent_center1, " c^2"
      PRINT *, "  Central baryon number density = ", this% nbar_center1, &
               " (M_sun^geo)^{-3} =", &
               this% nbar_center1/(MSun_geo*km2m*m2cm)**3, "cm^{-3}"
      PRINT *, "  Central baryon mass density = ", this% rho_center1, &
               " M_sun^geo (M_sun^geo)^{-3} =", &
               this% rho_center1/lorene2hydrobase*kg2g/(m2cm**3), "g cm^{-3}"
      PRINT *, "  Central energy density = ", this% energy_density_center1, &
               " M_sun^geo c^2 (M_sun^geo)^{-3}", &
               this% energy_density_center1/lorene2hydrobase*kg2g/(m2cm**3), &
               "g c^2 cm^{-3}"
      PRINT *, "  Central specific energy = ", this% specific_energy_center1, &
               " c^2"
      PRINT *, "  Central pressure = ", this% pressure_center1, &
               " M_sun^geo c^2 (M_sun^geo)^{-3}", &
               this% pressure_center1/lorene2hydrobase*kg2g/(m2cm**3), &
               "g c^2 cm^{-3}"
      PRINT *, " Hydro quantities at the center of star 2: "
      PRINT *, "  Central enthalpy = ", this% ent_center2, " c^2"
      PRINT *, "  Central baryon number density = ", this% nbar_center2, &
               " (M_sun^geo)^{-3} =", &
               this% nbar_center2/(MSun_geo*km2m*m2cm)**3, "cm^{-3}"
      PRINT *, "  Central baryon mass density = ", this% rho_center2, &
               " M_sun^geo (M_sun^geo)^{-3} =", &
               this% rho_center2/lorene2hydrobase*kg2g/(m2cm**3), "g cm^{-3}"
      PRINT *, "  Central energy density = ", this% energy_density_center2, &
               " M_sun^geo c^2 (M_sun^geo)^{-3}", &
               this% energy_density_center2/lorene2hydrobase*kg2g/(m2cm**3), &
               "g c^2 cm^{-3}"
      PRINT *, "  Central specific energy = ", this% specific_energy_center2, &
               " c^2"
      PRINT *, "  Central pressure = ", this% pressure_center2, &
               " M_sun^geo c^2 (M_sun^geo)^{-3}", &
               this% pressure_center2/lorene2hydrobase*kg2g/(m2cm**3), &
               "g c^2 cm^{-3}"
      PRINT *
      !IF( show_progress ) &
        PRINT *, " Equations of state for star 1 (EOS1) = ", TRIM(this% eos1)
      !IF( show_progress ) &
        PRINT *, " Equations of state for star 2 (EOS2) = ", TRIM(this% eos2)
      !IF( show_progress ) PRINT *

      IF( this% eos1_loreneid == 1 )THEN ! If the EOS is polytropic

        PRINT *, " Parameters for EOS1: "
        PRINT *, "  Polytopic index gamma_1 = ", this% gamma_1
        PRINT *, "  Pressure coefficient = ",&
                 this% kappa_1/k_lorene2hydrobase( this% gamma_1 ), &
                 "rho_nuc c^2 / n_nuc^gamma_1 = ", this% kappa_1, &
                 "[pure number]"
        PRINT *, " Parameters for EOS2: "
        PRINT *, "  Polytopic index gamma_2 = ", this% gamma_2
        PRINT *, "  Pressure coefficient = ",&
                 this% kappa_2/k_lorene2hydrobase( this% gamma_2 ), &
                 "rho_nuc c^2 / n_nuc^gamma_2 = ", this% kappa_2, &
                 "[pure number]"
        PRINT *

      ELSEIF( this% gamma0_1 /= 0 )THEN ! If the EOS is piecewise polytropic

        PRINT *, " Parameters for EOS1: "
        PRINT *, "  Number of polytropic indexes = ", this% npeos_1
        PRINT *, "  Polytopic index gamma0_1 = ", this% gamma0_1
        PRINT *, "  Polytopic index gamma1_1 = ", this% gamma1_1
        PRINT *, "  Polytopic index gamma2_1 = ", this% gamma2_1
        PRINT *, "  Polytopic index gamma3_1 = ", this% gamma3_1
        PRINT *, "  Pressure coefficient for the crust (here from SLy) = ",&
                 this% kappa0_1/k_lorene2hydrobase( this% gamma0_1 ), &
                 "rho_nuc c^2 / n_nuc^gamma0_1 = ", this% kappa0_1, &
                 "[pure number]"
        PRINT *, "  Pressure coefficient for the first polytrope = ",&
                 this% kappa1_1/k_lorene2hydrobase( this% gamma1_1 ), &
                 "rho_nuc c^2 / n_nuc^gamma1_1", this% kappa1_1, &
                 "[pure number]"
        PRINT *, "  Pressure coefficient for the second polytrope = ",&
                 this% kappa2_1/k_lorene2hydrobase( this% gamma2_1 ), &
                 "rho_nuc c^2 / n_nuc^gamma2_1", this% kappa2_1, &
                 "[pure number]"
        PRINT *, "  Pressure coefficient for the third polytrope = ",&
                 this% kappa3_1/k_lorene2hydrobase( this% gamma3_1 ), &
                 "rho_nuc c^2 / n_nuc^gamma3_1", this% kappa3_1, &
                 "[pure number]"
        PRINT *, "  Base 10 exponent of the pressure at the first fiducial " &
                 // "density (between gamma_0 and gamma_1) (dyne/cm^2)= ", &
                 this% logP1_1
        PRINT *, "  Base 10 exponent of first fiducial density (g/cm^3) = ", &
                 this% logRho0_1
        PRINT *, "  Base 10 exponent of second fiducial density (g/cm^3) = ",&
                 this% logRho1_1
        PRINT *, "  Base 10 exponent of third fiducial density (g/cm^3) = ", &
                 this% logRho2_1
        PRINT *
        PRINT *, " Parameters for EOS2: "
        PRINT *, "  Number of polytropic indexes = ", this% npeos_2
        PRINT *, "  Polytopic index gamma0_2 = ", this% gamma0_2
        PRINT *, "  Polytopic index gamma1_2 = ", this% gamma1_2
        PRINT *, "  Polytopic index gamma2_2 = ", this% gamma2_2
        PRINT *, "  Polytopic index gamma3_2 = ", this% gamma3_2
        PRINT *, "  Pressure coefficient for the crust (here from SLy) = ",&
                 this% kappa0_2/k_lorene2hydrobase( this% gamma0_2 ), &
                 "rho_nuc c^2 / n_nuc^gamma0_2 = ", this% kappa0_2, &
                 "[pure number]"
        PRINT *, "  Pressure coefficient for the first polytrope = ",&
                 this% kappa1_2/k_lorene2hydrobase( this% gamma1_2 ), &
                 "rho_nuc c^2 / n_nuc^gamma1_2", this% kappa1_2, &
                 "[pure number]"
        PRINT *, "  Pressure coefficient for the second polytrope = ",&
                 this% kappa2_2/k_lorene2hydrobase( this% gamma2_2 ), &
                 "rho_nuc c^2 / n_nuc^gamma2_2", this% kappa2_2, &
                 "[pure number]"
        PRINT *, "  Pressure coefficient for the third polytrope = ",&
                 this% kappa3_2/k_lorene2hydrobase( this% gamma3_2 ), &
                 "rho_nuc c^2 / n_nuc^gamma3_2", this% kappa3_2, &
                 "[pure number]"
        PRINT *, "  Base 10 exponent of the pressure at the first fiducial " &
                 // "density (between gamma_0 and gamma_1) (dyne/cm^2)= ", &
                 this% logP1_2
        PRINT *, "  Base 10 exponent of first fiducial density (g/cm^3) = ", &
                 this% logRho0_2
        PRINT *, "  Base 10 exponent of second fiducial density (g/cm^3) = ",&
                 this% logRho1_2
        PRINT *, "  Base 10 exponent of third fiducial density (g/cm^3) = ", &
                 this% logRho2_2
        PRINT *

      ELSEIF( this% eos1_loreneid == 17 .OR. this% eos1_loreneid == 20 )THEN
      ! If the EOS is tabulated

      ELSE

        PRINT *, "** ERROR in SUBROUTINE import_lorene_id_params!", &
                 " The equation of state is unknown!"
        STOP

      ENDIF

    ENDIF

  END PROCEDURE print_id_params


END SUBMODULE io
