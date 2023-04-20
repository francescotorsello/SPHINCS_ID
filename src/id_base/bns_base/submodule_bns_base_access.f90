! File:         submodule_bns_base_access.f90
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

SUBMODULE (bns_base) access

  !***************************************************
  !
  !# The module contains the implementation of the
  !  methods of TYPE [[bnsbase]] that allow to access
  !  PRIVATE members.
  !
  !  FT 12.07.2021
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
    !# Returns the baryon mass of NS `i_matter`-th star
    !  [\(M_\odot\)]
    !
    !  FT 27.10.2021
    !
    !************************************************

    IMPLICIT NONE

    CALL this% check_i_matter(i_matter)

    get_mass= this% mass(i_matter)

  END PROCEDURE get_mass


  MODULE PROCEDURE get_radii

    !************************************************
    !
    !# Returns the radii of the `i_matter`-th star
    !  [\(L_\odot\)]
    !
    !  FT 27.10.2021
    !
    !************************************************

    IMPLICIT NONE

    CALL this% check_i_matter(i_matter)

    get_radii= this% radii(i_matter,:)

  END PROCEDURE get_radii


  MODULE PROCEDURE get_center

    !************************************************
    !
    !# Returns the center of the `i_matter`-th star
    !  [\(L_\odot\)]
    !
    !  FT 27.10.2021
    !
    !************************************************

    IMPLICIT NONE

    CALL this% check_i_matter(i_matter)

    get_center= this% center(i_matter,:)

  END PROCEDURE get_center


  MODULE PROCEDURE get_barycenter

    !************************************************
    !
    !# Returns the barycenter of the `i_matter`-th star
    !  [\(L_\odot\)]
    !
    !  FT 27.10.2021
    !
    !************************************************

    IMPLICIT NONE

    CALL this% check_i_matter(i_matter)

    get_barycenter= this% barycenter(i_matter,:)

  END PROCEDURE get_barycenter


  MODULE PROCEDURE get_eos

    !************************************************
    !
    !# Returns the |eos| name of the `i_matter`-th star
    !  [\(L_\odot\)]
    !
    !  FT 27.10.2021
    !
    !************************************************

    IMPLICIT NONE

    CALL this% check_i_matter(i_matter)

    IF( i_matter == 1 ) get_eos= this% eos1
    IF( i_matter == 2 ) get_eos= this% eos2

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

    IF( i_matter == 1 ) get_eos_id= this% eos1_id
    IF( i_matter == 2 ) get_eos_id= this% eos2_id

  END PROCEDURE get_eos_id


  MODULE PROCEDURE get_gamma_1

    !************************************************
    !
    !# Returns the value of [[bnsbase:gamma_1]], the
    !  polytropic index for NS 1 with polytropic |eos|,
    !  not piecewise polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_gamma_1= this% gamma_1

  END PROCEDURE get_gamma_1


  MODULE PROCEDURE get_gamma_2

    !************************************************
    !
    !# Returns the value of [[bnsbase:gamma_2]], the
    !  polytropic index for NS 2 with polytropic |eos|,
    !  not piecewise polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_gamma_2= this% gamma_2

  END PROCEDURE get_gamma_2


  MODULE PROCEDURE get_kappa_1

    !************************************************
    !
    !# Returns the value of [[bnsbase:kappa_1]], the
    !  polytropic constant for NS 1 with polytropic
    !  |eos|, not piecewise polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_kappa_1= this% kappa_1

  END PROCEDURE get_kappa_1


  MODULE PROCEDURE get_kappa_2

    !************************************************
    !
    !# Returns the value of [[bnsbase:kappa_2]], the
    !  polytropic constant for NS 2 with polytropic
    !  |eos|, not piecewise polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_kappa_2= this% kappa_2

  END PROCEDURE get_kappa_2


  MODULE PROCEDURE get_angular_vel

    !************************************************
    !
    !# Returns the angular velocity of the system
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_angular_vel= this% angular_vel

  END PROCEDURE get_angular_vel


  MODULE PROCEDURE get_distance

    !************************************************
    !
    !# Returns the distance between the NSs
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_distance= this% distance

  END PROCEDURE get_distance


  MODULE PROCEDURE get_distance_com

    !************************************************
    !
    !# Returns the distance between the centers of
    !  mass of the NSs
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_distance_com= this% distance_com

  END PROCEDURE get_distance_com


  MODULE PROCEDURE get_mass1

    !************************************************
    !
    !# Returns the baryon mass of NS 1 [\(M_\odot\)]
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_mass1= this% mass(1)

  END PROCEDURE get_mass1


  MODULE PROCEDURE get_mass2

    !************************************************
    !
    !# Returns the baryon mass of NS 2 [\(M_\odot\)]
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_mass2= this% mass(2)

  END PROCEDURE get_mass2


  MODULE PROCEDURE get_grav_mass1

    !************************************************
    !
    !# Returns the gravitational mass of NS 1 [\(M_\odot\)]
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_grav_mass1= this% mass_grav(1)

  END PROCEDURE get_grav_mass1


  MODULE PROCEDURE get_grav_mass2

    !************************************************
    !
    !# Returns the gravitational mass of NS 2 [\(M_\odot\)]
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_grav_mass2= this% mass_grav(2)

  END PROCEDURE get_grav_mass2


  MODULE PROCEDURE get_adm_mass

    !************************************************
    !
    !# Returns the ADM mass of the system
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_adm_mass= this% adm_mass

  END PROCEDURE get_adm_mass


  MODULE PROCEDURE get_linear_momentum

    !************************************************
    !
    !# Returns the linear momentum of the system
    !
    !  FT 04.02.2022
    !
    !************************************************

    IMPLICIT NONE

    get_linear_momentum= [ this% linear_momentum_x, &
                           this% linear_momentum_y, &
                           this% linear_momentum_z ]

  END PROCEDURE get_linear_momentum


  MODULE PROCEDURE get_angular_momentum

    !************************************************
    !
    !# Returns the angular momentum of the system
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_angular_momentum= [ this% angular_momentum_x, &
                            this% angular_momentum_y, &
                            this% angular_momentum_z ]

  END PROCEDURE get_angular_momentum


  MODULE PROCEDURE get_radius1_x_comp

    !************************************************
    !
    !# Returns the radius of NS 1 along the x axis
    !  on the side of the companion
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_radius1_x_comp= this% radius1_x_comp

  END PROCEDURE get_radius1_x_comp


  MODULE PROCEDURE get_radius1_y

    !************************************************
    !
    !# Returns the radius of NS 1 along the y axis
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_radius1_y= this% radius1_y

  END PROCEDURE get_radius1_y


  MODULE PROCEDURE get_radius1_z

    !************************************************
    !
    !# Returns the radius of NS 1 along the z axis
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_radius1_z= this% radius1_z

  END PROCEDURE get_radius1_z


  MODULE PROCEDURE get_radius1_x_opp

    !************************************************
    !
    !# Returns the radius of NS 1 along the x axis
    !  on the side opposite to the companion
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_radius1_x_opp= this% radius1_x_opp

  END PROCEDURE get_radius1_x_opp


  MODULE PROCEDURE get_center1_x

    !************************************************
    !
    !# Returns the stellar center of NS 1, i.e., the
    !  origin of the LORENE chart centered on NS 1
    !
    !  FT 09.02.2021
    !
    !************************************************

    IMPLICIT NONE

    get_center1_x= this% center1_x

  END PROCEDURE get_center1_x


  MODULE PROCEDURE get_barycenter1_x

    !************************************************
    !
    !# Returns the barycenter of NS 1
    !
    !  FT 09.02.2021
    !
    !************************************************

    IMPLICIT NONE

    get_barycenter1_x= this% barycenter1_x

  END PROCEDURE get_barycenter1_x


  MODULE PROCEDURE get_radius2_x_comp

    !************************************************
    !
    !# Returns the radius of NS 2 along the x axis
    !  on the side of the companion
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_radius2_x_comp= this% radius2_x_comp

  END PROCEDURE get_radius2_x_comp


  MODULE PROCEDURE get_radius2_y

    !************************************************
    !
    !# Returns the radius of NS 2 along the y axis
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_radius2_y= this% radius2_y

  END PROCEDURE get_radius2_y


  MODULE PROCEDURE get_radius2_z

    !************************************************
    !
    !# Returns the radius of NS 2 along the z axis
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_radius2_z= this% radius2_z

  END PROCEDURE get_radius2_z


  MODULE PROCEDURE get_radius2_x_opp

    !************************************************
    !
    !# Returns the radius of NS 2 along the x axis
    !  on the side opposite to the companion
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_radius2_x_opp= this% radius2_x_opp

  END PROCEDURE get_radius2_x_opp


  MODULE PROCEDURE get_center2_x

    !************************************************
    !
    !# Returns the stellar center of NS 2, i.e., the
    !  origin of the LORENE chart centered on NS 2
    !
    !  FT 09.02.2021
    !
    !************************************************

    IMPLICIT NONE

    get_center2_x= this% center2_x

  END PROCEDURE get_center2_x


  MODULE PROCEDURE get_barycenter2_x

    !************************************************
    !
    !# Returns the barycenter of NS 2
    !
    !  FT 09.02.2021
    !
    !************************************************

    IMPLICIT NONE

    get_barycenter2_x= this% barycenter2_x

  END PROCEDURE get_barycenter2_x


  MODULE PROCEDURE get_ent_center1

    !************************************************
    !
    !# Returns the central enthalpy of NS 1
    !
    !  FT 12.02.2021
    !
    !************************************************

    IMPLICIT NONE

    get_ent_center1= this% ent_center1

  END PROCEDURE get_ent_center1


  MODULE PROCEDURE get_nbar_center1

    !************************************************
    !
    !# Returns the central baryon number density
    !  of NS 1
    !
    !  FT 12.02.2021
    !
    !************************************************

    IMPLICIT NONE

    get_nbar_center1= this% nbar_center1

  END PROCEDURE get_nbar_center1


  MODULE PROCEDURE get_rho_center1

    !************************************************
    !
    !# Returns the central baryon mass density
    !  of NS 1
    !
    !  FT 12.02.2021
    !
    !************************************************

    IMPLICIT NONE

    get_rho_center1= this% rho_center1

  END PROCEDURE get_rho_center1


  MODULE PROCEDURE get_energy_density_center1

    !************************************************
    !
    !# Returns the central energy density of NS 1
    !
    !  FT 12.02.2021
    !
    !************************************************

    IMPLICIT NONE

    get_energy_density_center1= this% energy_density_center1

  END PROCEDURE get_energy_density_center1


  MODULE PROCEDURE get_specific_energy_center1

    !************************************************
    !
    !# Returns the central specific energy of NS 1
    !
    !  FT 12.02.2021
    !
    !************************************************

    IMPLICIT NONE

    get_specific_energy_center1= this% specific_energy_center1

  END PROCEDURE get_specific_energy_center1


  MODULE PROCEDURE get_pressure_center1

    !************************************************
    !
    !# Returns the central pressure of NS 1
    !
    !  FT 12.02.2021
    !
    !************************************************

    IMPLICIT NONE

    get_pressure_center1= this% pressure_center1

  END PROCEDURE get_pressure_center1


  MODULE PROCEDURE get_ent_center2

    !************************************************
    !
    !# Returns the central enthalpy of NS 2
    !
    !  FT 12.02.2021
    !
    !************************************************

    IMPLICIT NONE

    get_ent_center2= this% ent_center2

  END PROCEDURE get_ent_center2


  MODULE PROCEDURE get_nbar_center2

    !************************************************
    !
    !# Returns the central baryon number density
    !  of NS 2
    !
    !  FT 12.02.2021
    !
    !************************************************

    IMPLICIT NONE

    get_nbar_center2= this% nbar_center2

  END PROCEDURE get_nbar_center2


  MODULE PROCEDURE get_rho_center2

    !************************************************
    !
    !# Returns the central baryon mass density
    !  of NS 2
    !
    !  FT 12.02.2021
    !
    !************************************************

    IMPLICIT NONE

    get_rho_center2= this% rho_center2

  END PROCEDURE get_rho_center2


  MODULE PROCEDURE get_energy_density_center2

    !************************************************
    !
    !# Returns the central energy density of NS 2
    !
    !  FT 12.02.2021
    !
    !************************************************

    IMPLICIT NONE

    get_energy_density_center2= this% energy_density_center2

  END PROCEDURE get_energy_density_center2


  MODULE PROCEDURE get_specific_energy_center2

    !************************************************
    !
    !# Returns the central specific energy of NS 2
    !
    !  FT 12.02.2021
    !
    !************************************************

    IMPLICIT NONE

    get_specific_energy_center2= this% specific_energy_center2

  END PROCEDURE get_specific_energy_center2


  MODULE PROCEDURE get_pressure_center2

    !************************************************
    !
    !# Returns the central pressure of NS 2
    !
    !  FT 12.02.2021
    !
    !************************************************

    IMPLICIT NONE

    get_pressure_center2= this% pressure_center2

  END PROCEDURE get_pressure_center2


  MODULE PROCEDURE get_eos1

    !************************************************
    !
    !# Returns the name of the |eos| for NS 1
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_eos1= this% eos1

  END PROCEDURE get_eos1


  MODULE PROCEDURE get_eos2

    !************************************************
    !
    !# Returns the name of the |eos| for NS 2
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_eos2= this% eos2

  END PROCEDURE get_eos2


  MODULE PROCEDURE get_npeos_1

    !************************************************
    !
    !# Returns the identifier of the |eos| for NS 1
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_npeos_1= this% npeos_1

  END PROCEDURE get_npeos_1


  MODULE PROCEDURE get_npeos_2

    !************************************************
    !
    !# Returns the identifier of the |eos| for NS 2
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_npeos_2= this% npeos_2

  END PROCEDURE get_npeos_2


  MODULE PROCEDURE get_gamma0_1

    !************************************************
    !
    !# Returns the value of [[bnsbase:gamma0_1]], the crust's
    !  polytropic index for NS 1 with piecewise
    !  polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_gamma0_1= this% gamma0_1

  END PROCEDURE get_gamma0_1


  MODULE PROCEDURE get_gamma0_2

    !************************************************
    !
    !# Returns the value of [[bnsbase:gamma0_2]], the crust's
    !  polytropic index for NS 2 with piecewise
    !  polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_gamma0_2= this% gamma0_2

  END PROCEDURE get_gamma0_2


  MODULE PROCEDURE get_gamma1_1

    !************************************************
    !
    !# Returns the value of [[bnsbase:gamma1_1]], the first
    !  polytropic index for NS 1 with piecewise
    !  polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_gamma1_1= this% gamma1_1

  END PROCEDURE get_gamma1_1


  MODULE PROCEDURE get_gamma1_2

    !************************************************
    !
    !# Returns the value of [[bnsbase:gamma1_2]], the first
    !  polytropic index for NS 2 with piecewise
    !  polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_gamma1_2= this% gamma1_2

  END PROCEDURE get_gamma1_2


  MODULE PROCEDURE get_gamma2_1

    !************************************************
    !
    !# Returns the value of [[bnsbase:gamma2_1]], the second
    !  polytropic index for NS 2 with piecewise
    !  polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_gamma2_1= this% gamma2_1

  END PROCEDURE get_gamma2_1


  MODULE PROCEDURE get_gamma2_2

    !************************************************
    !
    !# Returns the value of [[bnsbase:gamma2_2]], the second
    !  polytropic index for NS 2 with piecewise
    !  polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_gamma2_2= this% gamma2_2

  END PROCEDURE get_gamma2_2


  MODULE PROCEDURE get_gamma3_1

    !************************************************
    !
    !# Returns the value of [[bnsbase:gamma3_1]], the third
    !  polytropic index for NS 1 with piecewise
    !  polytropic |eos| (innermost index)
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_gamma3_1= this% gamma3_1

  END PROCEDURE get_gamma3_1


  MODULE PROCEDURE get_gamma3_2

    !************************************************
    !
    !# Returns the value of [[bnsbase:gamma3_2]], the third
    !  polytropic index for NS 2 with piecewise
    !  polytropic |eos| (innermost index)
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_gamma3_2= this% gamma3_2

  END PROCEDURE get_gamma3_2


  MODULE PROCEDURE get_kappa0_1

    !************************************************
    !
    !# Returns the value of [[bnsbase:kappa0_1]], the crust's
    !  polytropic constant for NS 1 with piecewise
    !  polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_kappa0_1= this% kappa0_1

  END PROCEDURE get_kappa0_1


  MODULE PROCEDURE get_kappa1_1

    !************************************************
    !
    !# Returns the value of [[bnsbase:kappa1_1]], the first
    !  polytropic constant for NS 1 with piecewise
    !  polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_kappa1_1= this% kappa1_1

  END PROCEDURE get_kappa1_1


  MODULE PROCEDURE get_kappa2_1

    !************************************************
    !
    !# Returns the value of [[bnsbase:kappa2_1]], the second
    !  polytropic constant for NS 1 with piecewise
    !  polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_kappa2_1= this% kappa2_1

  END PROCEDURE get_kappa2_1


  MODULE PROCEDURE get_kappa3_1

    !************************************************
    !
    !# Returns the value of [[bnsbase:kappa3_1]], the third
    !  polytropic constant for NS 1 with piecewise
    !  polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_kappa3_1= this% kappa3_1

  END PROCEDURE get_kappa3_1


  MODULE PROCEDURE get_kappa0_2

    !************************************************
    !
    !# Returns the value of [[bnsbase:kappa0_2]], the crust's
    !  polytropic constant for NS 2 with piecewise
    !  polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_kappa0_2= this% kappa0_2

  END PROCEDURE get_kappa0_2


  MODULE PROCEDURE get_kappa1_2

    !************************************************
    !
    !# Returns the value of [[bnsbase:kappa1_2]], the first
    !  polytropic constant for NS 2 with piecewise
    !  polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_kappa1_2= this% kappa1_2

  END PROCEDURE get_kappa1_2


  MODULE PROCEDURE get_kappa2_2

    !************************************************
    !
    !# Returns the value of [[bnsbase:kappa2_2]], the second
    !  polytropic constant for NS 2 with piecewise
    !  polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_kappa2_2= this% kappa2_2

  END PROCEDURE get_kappa2_2


  MODULE PROCEDURE get_kappa3_2

    !************************************************
    !
    !# Returns the value of [[bnsbase:kappa3_2]], the third
    !  polytropic constant for NS 2 with piecewise
    !  polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_kappa3_2= this% kappa3_2

  END PROCEDURE get_kappa3_2


  MODULE PROCEDURE get_logp1_1

    !************************************************
    !
    !# Returns the value of [[bnsbase:logp1_1]], the base 10
    !  logarithm of the pressure where the gamma1_1
    !  polytrope starts, for NS 1 with piecewise
    !  polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_logp1_1= this% logp1_1

  END PROCEDURE get_logp1_1


  MODULE PROCEDURE get_logp1_2

    !************************************************
    !
    !# Returns the value of [[bnsbase:logp1_2]], the base 10
    !  logarithm of the pressure where the gamma1_2
    !  polytrope starts, for NS 2 with piecewise
    !  polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_logp1_2= this% logp1_2

  END PROCEDURE get_logp1_2


  MODULE PROCEDURE get_logRho0_1

    !************************************************
    !
    !# Returns the value of [[bnsbase:logRho0_1]], the base 10
    !  logarithm of the mass density where the
    !  gamma1_1 polytrope starts, for NS 1 with
    !  piecewise polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_logRho0_1= this% logRho0_1

  END PROCEDURE get_logRho0_1


  MODULE PROCEDURE get_logRho0_2

    !************************************************
    !
    !# Returns the value of [[bnsbase:logRho0_2]], the base 10
    !  logarithm of the mass density where the
    !  gamma1_2 polytrope starts, for NS 2 with
    !  piecewise polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_logRho0_2= this% logRho0_2

  END PROCEDURE get_logRho0_2


  MODULE PROCEDURE get_logRho1_1

    !************************************************
    !
    !# Returns the value of [[bnsbase:logRho1_1]], the base 10
    !  logarithm of the mass density where the
    !  gamma2_1 polytrope starts, for NS 1 with
    !  piecewise polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_logRho1_1= this% logRho1_1

  END PROCEDURE get_logRho1_1


  MODULE PROCEDURE get_logRho1_2

    !************************************************
    !
    !# Returns the value of [[bnsbase:logRho1_2]], the base 10
    !  logarithm of the mass density where the
    !  gamma2_2 polytrope starts, for NS 2 with
    !  piecewise polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_logRho1_2= this% logRho1_2

  END PROCEDURE get_logRho1_2


  MODULE PROCEDURE get_logRho2_1

    !************************************************
    !
    !# Returns the value of [[bnsbase:logRho2_1]], the base 10
    !  logarithm of the mass density where the
    !  gamma3_1 polytrope starts, for NS 1 with
    !  piecewise polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_logRho2_1= this% logRho2_1

  END PROCEDURE get_logRho2_1


  MODULE PROCEDURE get_logRho2_2

    !************************************************
    !
    !# Returns the value of [[bnsbase:logRho2_2]]]], the base 10
    !  logarithm of the mass density where the
    !  gamma3_2 polytrope starts, for NS 2 with
    !  piecewise polytropic |eos|
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_logRho2_2= this% logRho2_2

  END PROCEDURE get_logRho2_2


END SUBMODULE access
