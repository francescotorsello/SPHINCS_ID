! File:         submodule_diffstar_base_access.f90
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

SUBMODULE (diffstar_base) access

  !***************************************************
  !
  !# The module contains the implementation of the
  !  methods of TYPE diffstarbase that allow to
  !  access PRIVATE members.
  !
  !  FT 22.10.2021
  !
  !***************************************************


  IMPLICIT NONE


  CONTAINS


  !----------------------------!
  !--  OVERRIDING FUNCTIONS  --!
  !----------------------------!


  MODULE PROCEDURE get_mass

    !************************************************
    !
    !# Returns the baryon mass of the DRS [\(M_\odot\)]
    !
    !  FT 2.11.2021
    !
    !************************************************

    IMPLICIT NONE

    CALL this% check_i_matter(i_matter)

    get_mass= this% mass

  END PROCEDURE get_mass


  MODULE PROCEDURE get_adm_mass

    !************************************************
    !
    !# Returns 0 (the ADM mass is not necessarily
    !  known for this TYPE)
    !
    !  FT 03.11.2022
    !
    !************************************************

    USE utility, ONLY: zero

    IMPLICIT NONE

    get_adm_mass= zero

  END PROCEDURE get_adm_mass


  MODULE PROCEDURE get_radii

    !************************************************
    !
    !# Returns the radii of the DRS [\(L_\odot\)]
    !
    !  FT 2.11.2021
    !
    !************************************************

    IMPLICIT NONE

    CALL this% check_i_matter(i_matter)

    get_radii= this% radii(:)

  END PROCEDURE get_radii


  MODULE PROCEDURE get_center

    !************************************************
    !
    !# Returns the center of the DRS [\(L_\odot\)]
    !
    !  FT 2.11.2021
    !
    !************************************************

    IMPLICIT NONE

    CALL this% check_i_matter(i_matter)

    get_center= this% center(:)

  END PROCEDURE get_center


  MODULE PROCEDURE get_barycenter

    !************************************************
    !
    !# Returns the barycenter of the DRS [\(L_\odot\)]
    !
    !  FT 2.11.2021
    !
    !************************************************

    IMPLICIT NONE

    CALL this% check_i_matter(i_matter)

    get_barycenter= this% barycenter(:)

  END PROCEDURE get_barycenter


  MODULE PROCEDURE get_eos

    !************************************************
    !
    !# Returns the |eos| name of the DRS
    !
    !  FT 2.11.2021
    !
    !************************************************

    IMPLICIT NONE

    CALL this% check_i_matter(i_matter)

    get_eos= this% eos

  END PROCEDURE get_eos


  !-----------------!
  !--  FUNCTIONS  --!
  !-----------------!


  MODULE PROCEDURE get_eos_id

    !************************************************
    !
    !# Returns the |eos| identifier of the
    !  `i_matter`-th star
    !
    !  FT 2.12.2022
    !
    !************************************************

    IMPLICIT NONE

    CALL this% check_i_matter(i_matter)

    get_eos_id= this% eos_id

  END PROCEDURE get_eos_id


  MODULE PROCEDURE get_gamma

    !************************************************
    !
    !# Returns the value of [[gamma]], the
    !  polytropic index for polytropic EOS,
    !  not piecewise polytropic EOS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_gamma= this% gamma

  END PROCEDURE get_gamma


  MODULE PROCEDURE get_kappa

    !************************************************
    !
    !# Returns the value of [[kappa]], the
    !  polytropic constant for polytropic
    !  EOS, not piecewise polytropic EOS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_kappa= this% kappa

  END PROCEDURE get_kappa


  MODULE PROCEDURE get_omega_c

    !************************************************
    !
    !# Returns [[omega_c]], the central  angular
    !  velocity of the system
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_omega_c= this% omega_c

  END PROCEDURE get_omega_c


  MODULE PROCEDURE get_mass_grav

    !************************************************
    !
    !# Returns the gravitational mass of the DRS [\(M_\odot\)]
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_mass_grav= this% mass_grav

  END PROCEDURE get_mass_grav


  MODULE PROCEDURE get_angular_momentum

    !************************************************
    !
    !# Returns the angular momentum of the DRS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_angular_momentum= this% angular_momentum

  END PROCEDURE get_angular_momentum


  MODULE PROCEDURE get_tsw

    !************************************************
    !
    !# Returns [[diffstarbase:tsw]], the ratio \(T/W\)
    !  between the kinetic and gravitational potential
    !  energy of the DRS
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_tsw= this% tsw

  END PROCEDURE get_tsw


  MODULE PROCEDURE get_grv2

    !************************************************
    !
    !# Returns [[diffstarbase:grv2]], the error on the
    !  virial identity  \({\rm GRV2}\).
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_grv2= this% grv2

  END PROCEDURE get_grv2


  MODULE PROCEDURE get_grv3

    !************************************************
    !
    !# Returns [[diffstarbase:grv3]], the error on the
    !  virial identity  \({\rm GRV3}\).
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_grv3= this% grv3

  END PROCEDURE get_grv3


  MODULE PROCEDURE get_r_circ

    !************************************************
    !
    !# Returns [[diffstarbase:r_circ]], the
    !  circumferential radius of the DRS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_r_circ= this% r_circ

  END PROCEDURE get_r_circ


  MODULE PROCEDURE get_r_mean

    !************************************************
    !
    !# Returns [[diffstarbase:r_mean]], the
    !  circumferential radius of the DRS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_r_mean= this% r_mean

  END PROCEDURE get_r_mean


  MODULE PROCEDURE get_r_eq

    !************************************************
    !
    !# Returns [[diffstarbase:r_eq]], the
    !  equatorial radius of the DRS at \(\phi=0\)
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_r_eq= this% r_eq

  END PROCEDURE get_r_eq


  MODULE PROCEDURE get_r_eq_pi2

    !************************************************
    !
    !# Returns [[diffstarbase:r_eq_pi2]], the
    !  equatorial radius of the DRS at \(\phi=\dfrac{\pi}{2}\)
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_r_eq_pi2= this% r_eq_pi2

  END PROCEDURE get_r_eq_pi2


  MODULE PROCEDURE get_r_eq_pi

    !************************************************
    !
    !# Returns [[diffstarbase:r_eq_pi]], the
    !  equatorial radius of the DRS at \(\phi=\pi\)
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_r_eq_pi= this% r_eq_pi

  END PROCEDURE get_r_eq_pi


  MODULE PROCEDURE get_r_eq_3pi2

    !************************************************
    !
    !# Returns [[diffstarbase:r_eq_3pi2]], the
    !  equatorial radius of the DRS at \(\phi=\dfrac{3\pi}{2}\)
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_r_eq_3pi2= this% r_eq_3pi2

  END PROCEDURE get_r_eq_3pi2


  MODULE PROCEDURE get_r_pole

    !************************************************
    !
    !# Returns [[diffstarbase:r_pole]], the
    !  polar radius of the DRS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_r_pole= this% r_pole

  END PROCEDURE get_r_pole


  MODULE PROCEDURE get_r_ratio

    !************************************************
    !
    !# Returns [[diffstarbase:r_ratio]], the
    !  Ratio [[diffstarbase:r_pole]]/[[diffstarbase:r_eq]]
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_r_ratio= this% r_ratio

  END PROCEDURE get_r_ratio


  MODULE PROCEDURE get_r_isco

    !************************************************
    !
    !# Returns [[diffstarbase:r_isco]], the
    !  radius of the Innermost Stable Circular Orbit (ISCO)
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_r_isco= this% r_isco

  END PROCEDURE get_r_isco


  MODULE PROCEDURE get_f_isco

    !************************************************
    !
    !# Returns [[diffstarbase:f_isco]], the orbital
    !  frequency of the Innermost Stable Circular Orbit
    !  (ISCO)
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_f_isco= this% f_isco

  END PROCEDURE get_f_isco


  MODULE PROCEDURE get_specific_energy_isco

    !************************************************
    !
    !# Returns [[diffstarbase:specific_energy_isco]],
    !  the specific energy of a test particle at the
    !  Innermost Stable Circular Orbit (ISCO)
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_specific_energy_isco= this% specific_energy_isco

  END PROCEDURE get_specific_energy_isco


  MODULE PROCEDURE get_specific_angular_momentum_isco

    !************************************************
    !
    !# Returns [[diffstarbase:specific_angular_momentum_isco]],
    !  the specific angular momentum of a test particle
    !  at the Innermost Stable Circular Orbit (ISCO)
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_specific_angular_momentum_isco= this% specific_angular_momentum_isco

  END PROCEDURE get_specific_angular_momentum_isco


  MODULE PROCEDURE get_surface_area

    !************************************************
    !
    !# Returns [[diffstarbase:surface_area]], the
    !  surface area of the DRS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_surface_area= this% surface_area

  END PROCEDURE get_surface_area


  MODULE PROCEDURE get_area_radius

    !************************************************
    !
    !# Returns [[diffstarbase:area_radius]], the
    !  areal radius of the DRS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_area_radius= this% area_radius

  END PROCEDURE get_area_radius


  MODULE PROCEDURE get_ent_center

    !************************************************
    !
    !# Returns the central enthalpy of the DRS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_ent_center= this% ent_center

  END PROCEDURE get_ent_center


  MODULE PROCEDURE get_nbar_center

    !************************************************
    !
    !# Returns the central baryon number density
    !  of the DRS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_nbar_center= this% nbar_center

  END PROCEDURE get_nbar_center


  MODULE PROCEDURE get_rho_center

    !************************************************
    !
    !# Returns the central baryon mass density
    !  of the DRS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_rho_center= this% rho_center

  END PROCEDURE get_rho_center


  MODULE PROCEDURE get_energy_density_center

    !************************************************
    !
    !# Returns the central energy density of the DRS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_energy_density_center= this% energy_density_center

  END PROCEDURE get_energy_density_center


  MODULE PROCEDURE get_specific_energy_center

    !************************************************
    !
    !# Returns the central specific energy of the DRS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_specific_energy_center= this% specific_energy_center

  END PROCEDURE get_specific_energy_center


  MODULE PROCEDURE get_pressure_center

    !************************************************
    !
    !# Returns the central pressure of the DRS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_pressure_center= this% pressure_center

  END PROCEDURE get_pressure_center


  MODULE PROCEDURE get_npeos

    !************************************************
    !
    !# Returns the identifier of the EOS for the DRS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_npeos= this% npeos

  END PROCEDURE get_npeos


  MODULE PROCEDURE get_gamma0

    !************************************************
    !
    !# Returns the value of [[gamma0]], the crust's
    !  polytropic index for the DRS with piecewise
    !  polytropic EOS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_gamma0= this% gamma0

  END PROCEDURE get_gamma0


  MODULE PROCEDURE get_gamma1

    !************************************************
    !
    !# Returns the value of [[gamma1]], the first
    !  polytropic index for the DRS with piecewise
    !  polytropic EOS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_gamma1= this% gamma1

  END PROCEDURE get_gamma1


  MODULE PROCEDURE get_gamma2

    !************************************************
    !
    !# Returns the value of [[gamma2]], the second
    !  polytropic index for NS 2 with piecewise
    !  polytropic EOS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_gamma2= this% gamma2

  END PROCEDURE get_gamma2


  MODULE PROCEDURE get_gamma3

    !************************************************
    !
    !# Returns the value of [[gamma3]], the third
    !  polytropic index for the DRS with piecewise
    !  polytropic EOS (innermost index)
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_gamma3= this% gamma3

  END PROCEDURE get_gamma3


  MODULE PROCEDURE get_kappa0

    !************************************************
    !
    !# Returns the value of [[kappa0]], the crust's
    !  polytropic constant for the DRS with piecewise
    !  polytropic EOS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_kappa0= this% kappa0

  END PROCEDURE get_kappa0


  MODULE PROCEDURE get_kappa1

    !************************************************
    !
    !# Returns the value of [[kappa1]], the first
    !  polytropic constant for the DRS with piecewise
    !  polytropic EOS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_kappa1= this% kappa1

  END PROCEDURE get_kappa1


  MODULE PROCEDURE get_kappa2

    !************************************************
    !
    !# Returns the value of [[kappa2]], the second
    !  polytropic constant for the DRS with piecewise
    !  polytropic EOS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_kappa2= this% kappa2

  END PROCEDURE get_kappa2


  MODULE PROCEDURE get_kappa3

    !************************************************
    !
    !# Returns the value of [[kappa3]], the third
    !  polytropic constant for the DRS with piecewise
    !  polytropic EOS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_kappa3= this% kappa3

  END PROCEDURE get_kappa3


  MODULE PROCEDURE get_logp1

    !************************************************
    !
    !# Returns the value of [[logp1]], the base 10
    !  logarithm of the pressure where the gamma1
    !  polytrope starts, for the DRS with piecewise
    !  polytropic EOS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_logp1= this% logp1

  END PROCEDURE get_logp1



  MODULE PROCEDURE get_logRho0

    !************************************************
    !
    !# Returns the value of [[logRho0]], the base 10
    !  logarithm of the mass density where the
    !  gamma1 polytrope starts, for the DRS with
    !  piecewise polytropic EOS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_logRho0= this% logRho0

  END PROCEDURE get_logRho0


  MODULE PROCEDURE get_logRho1

    !************************************************
    !
    !# Returns the value of [[logRho1]], the base 10
    !  logarithm of the mass density where the
    !  gamma2 polytrope starts, for the DRS with
    !  piecewise polytropic EOS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_logRho1= this% logRho1

  END PROCEDURE get_logRho1


  MODULE PROCEDURE get_logRho2

    !************************************************
    !
    !# Returns the value of [[logRho2]], the base 10
    !  logarithm of the mass density where the
    !  gamma3 polytrope starts, for the DRS with
    !  piecewise polytropic EOS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_logRho2= this% logRho2

  END PROCEDURE get_logRho2


END SUBMODULE access
