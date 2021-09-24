! File:         submodule_bns_access.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (bns_base) bns_access

  !***************************************************
  !
  !# The module contains the implementation of the
  !  methods of TYPE bns that allow to access PRIVATE
  !  members.
  !
  !  FT 12.07.2021
  !
  !***************************************************


  IMPLICIT NONE


  CONTAINS


  !-----------------!
  !--  FUNCTIONS  --!
  !-----------------!


  MODULE PROCEDURE get_gamma_1

    !************************************************
    !
    !# Returns the value of [[gamma_1]], the
    !  polytropic index for NS 1 with polytropic EOS,
    !  not piecewise polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_gamma_1= THIS% gamma_1

  END PROCEDURE get_gamma_1


  MODULE PROCEDURE get_gamma_2

    !************************************************
    !
    !# Returns the value of [[gamma_2]], the
    !  polytropic index for NS 2 with polytropic EOS,
    !  not piecewise polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_gamma_2= THIS% gamma_2

  END PROCEDURE get_gamma_2


  MODULE PROCEDURE get_kappa_1

    !************************************************
    !
    !# Returns the value of [[kappa_1]], the
    !  polytropic constant for NS 1 with polytropic
    !  EOS, not piecewise polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_kappa_1= THIS% kappa_1

  END PROCEDURE get_kappa_1


  MODULE PROCEDURE get_kappa_2

    !************************************************
    !
    !# Returns the value of [[kappa_2]], the
    !  polytropic constant for NS 2 with polytropic
    !  EOS, not piecewise polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_kappa_2= THIS% kappa_2

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

    get_angular_vel= THIS% angular_vel

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

    get_distance= THIS% distance

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

    get_distance_com= THIS% distance_com

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

    get_mass1= THIS% mass1

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

    get_mass2= THIS% mass2

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

    get_grav_mass1= THIS% mass_grav1

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

    get_grav_mass2= THIS% mass_grav2

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

    get_adm_mass= THIS% adm_mass

  END PROCEDURE get_adm_mass


  MODULE PROCEDURE get_angular_momentum

    !************************************************
    !
    !# Returns the angular momentum of the system
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_angular_momentum= THIS% angular_momentum

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

    get_radius1_x_comp= THIS% radius1_x_comp

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

    get_radius1_y= THIS% radius1_y

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

    get_radius1_z= THIS% radius1_z

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

    get_radius1_x_opp= THIS% radius1_x_opp

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

    get_center1_x= THIS% center1_x

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

    get_barycenter1_x= THIS% barycenter1_x

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

    get_radius2_x_comp= THIS% radius2_x_comp

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

    get_radius2_y= THIS% radius2_y

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

    get_radius2_z= THIS% radius2_z

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

    get_radius2_x_opp= THIS% radius2_x_opp

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

    get_center2_x= THIS% center2_x

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

    get_barycenter2_x= THIS% barycenter2_x

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

    get_ent_center1= THIS% ent_center1

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

    get_nbar_center1= THIS% nbar_center1

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

    get_rho_center1= THIS% rho_center1

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

    get_energy_density_center1= THIS% energy_density_center1

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

    get_specific_energy_center1= THIS% specific_energy_center1

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

    get_pressure_center1= THIS% pressure_center1

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

    get_ent_center2= THIS% ent_center2

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

    get_nbar_center2= THIS% nbar_center2

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

    get_rho_center2= THIS% rho_center2

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

    get_energy_density_center2= THIS% energy_density_center2

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

    get_specific_energy_center2= THIS% specific_energy_center2

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

    get_pressure_center2= THIS% pressure_center2

  END PROCEDURE get_pressure_center2


  MODULE PROCEDURE get_eos1

    !************************************************
    !
    !# Returns the name of the EOS for NS 1
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_eos1= THIS% eos1

  END PROCEDURE get_eos1


  MODULE PROCEDURE get_eos2

    !************************************************
    !
    !# Returns the name of the EOS for NS 2
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_eos2= THIS% eos2

  END PROCEDURE get_eos2


  MODULE PROCEDURE get_eos1_id

    !**************************************************
    !
    !# Returns the LORENE ID-number of the EOS for NS 1
    !
    !  FT
    !
    !**************************************************

    IMPLICIT NONE

    get_eos1_id= THIS% eos1_id

  END PROCEDURE get_eos1_id


  MODULE PROCEDURE get_eos2_id

    !**************************************************
    !
    !# Returns the LORENE ID-number of the EOS for NS 2
    !
    !  FT
    !
    !**************************************************

    IMPLICIT NONE

    get_eos2_id= THIS% eos2_id

  END PROCEDURE get_eos2_id


  MODULE PROCEDURE get_npeos_1

    !************************************************
    !
    !# Returns the identifier of the EOS for NS 1
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_npeos_1= THIS% npeos_1

  END PROCEDURE get_npeos_1


  MODULE PROCEDURE get_npeos_2

    !************************************************
    !
    !# Returns the identifier of the EOS for NS 2
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_npeos_2= THIS% npeos_2

  END PROCEDURE get_npeos_2


  MODULE PROCEDURE get_gamma0_1

    !************************************************
    !
    !# Returns the value of [[gamma0_1]], the crust's
    !  polytropic index for NS 1 with piecewise
    !  polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_gamma0_1= THIS% gamma0_1

  END PROCEDURE get_gamma0_1


  MODULE PROCEDURE get_gamma0_2

    !************************************************
    !
    !# Returns the value of [[gamma0_2]], the crust's
    !  polytropic index for NS 2 with piecewise
    !  polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_gamma0_2= THIS% gamma0_2

  END PROCEDURE get_gamma0_2


  MODULE PROCEDURE get_gamma1_1

    !************************************************
    !
    !# Returns the value of [[gamma1_1]], the first
    !  polytropic index for NS 1 with piecewise
    !  polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_gamma1_1= THIS% gamma1_1

  END PROCEDURE get_gamma1_1


  MODULE PROCEDURE get_gamma1_2

    !************************************************
    !
    !# Returns the value of [[gamma1_2]], the first
    !  polytropic index for NS 2 with piecewise
    !  polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_gamma1_2= THIS% gamma1_2

  END PROCEDURE get_gamma1_2


  MODULE PROCEDURE get_gamma2_1

    !************************************************
    !
    !# Returns the value of [[gamma2_1]], the second
    !  polytropic index for NS 2 with piecewise
    !  polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_gamma2_1= THIS% gamma2_1

  END PROCEDURE get_gamma2_1


  MODULE PROCEDURE get_gamma2_2

    !************************************************
    !
    !# Returns the value of [[gamma2_2]], the second
    !  polytropic index for NS 2 with piecewise
    !  polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_gamma2_2= THIS% gamma2_2

  END PROCEDURE get_gamma2_2


  MODULE PROCEDURE get_gamma3_1

    !************************************************
    !
    !# Returns the value of [[gamma3_1]], the third
    !  polytropic index for NS 1 with piecewise
    !  polytropic EOS (innermost index)
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_gamma3_1= THIS% gamma3_1

  END PROCEDURE get_gamma3_1


  MODULE PROCEDURE get_gamma3_2

    !************************************************
    !
    !# Returns the value of [[gamma3_2]], the third
    !  polytropic index for NS 2 with piecewise
    !  polytropic EOS (innermost index)
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_gamma3_2= THIS% gamma3_2

  END PROCEDURE get_gamma3_2


  MODULE PROCEDURE get_kappa0_1

    !************************************************
    !
    !# Returns the value of [[kappa0_1]], the crust's
    !  polytropic constant for NS 1 with piecewise
    !  polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_kappa0_1= THIS% kappa0_1

  END PROCEDURE get_kappa0_1


  MODULE PROCEDURE get_kappa1_1

    !************************************************
    !
    !# Returns the value of [[kappa1_1]], the first
    !  polytropic constant for NS 1 with piecewise
    !  polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_kappa1_1= THIS% kappa1_1

  END PROCEDURE get_kappa1_1


  MODULE PROCEDURE get_kappa2_1

    !************************************************
    !
    !# Returns the value of [[kappa2_1]], the second
    !  polytropic constant for NS 1 with piecewise
    !  polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_kappa2_1= THIS% kappa2_1

  END PROCEDURE get_kappa2_1


  MODULE PROCEDURE get_kappa3_1

    !************************************************
    !
    !# Returns the value of [[kappa3_1]], the third
    !  polytropic constant for NS 1 with piecewise
    !  polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_kappa3_1= THIS% kappa3_1

  END PROCEDURE get_kappa3_1


  MODULE PROCEDURE get_kappa0_2

    !************************************************
    !
    !# Returns the value of [[kappa0_2]], the crust's
    !  polytropic constant for NS 2 with piecewise
    !  polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_kappa0_2= THIS% kappa0_2

  END PROCEDURE get_kappa0_2


  MODULE PROCEDURE get_kappa1_2

    !************************************************
    !
    !# Returns the value of [[kappa1_2]], the first
    !  polytropic constant for NS 2 with piecewise
    !  polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_kappa1_2= THIS% kappa1_2

  END PROCEDURE get_kappa1_2


  MODULE PROCEDURE get_kappa2_2

    !************************************************
    !
    !# Returns the value of [[kappa2_2]], the second
    !  polytropic constant for NS 2 with piecewise
    !  polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_kappa2_2= THIS% kappa2_2

  END PROCEDURE get_kappa2_2


  MODULE PROCEDURE get_kappa3_2

    !************************************************
    !
    !# Returns the value of [[kappa3_2]], the third
    !  polytropic constant for NS 2 with piecewise
    !  polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_kappa3_2= THIS% kappa3_2

  END PROCEDURE get_kappa3_2


  MODULE PROCEDURE get_logp1_1

    !************************************************
    !
    !# Returns the value of [[logp1_1]], the base 10
    !  logarithm of the pressure where the gamma1_1
    !  polytrope starts, for NS 1 with piecewise
    !  polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_logp1_1= THIS% logp1_1

  END PROCEDURE get_logp1_1


  MODULE PROCEDURE get_logp1_2

    !************************************************
    !
    !# Returns the value of [[logp1_2]], the base 10
    !  logarithm of the pressure where the gamma1_2
    !  polytrope starts, for NS 2 with piecewise
    !  polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_logp1_2= THIS% logp1_2

  END PROCEDURE get_logp1_2


  MODULE PROCEDURE get_logRho0_1

    !************************************************
    !
    !# Returns the value of [[logRho0_1]], the base 10
    !  logarithm of the mass density where the
    !  gamma1_1 polytrope starts, for NS 1 with
    !  piecewise polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_logRho0_1= THIS% logRho0_1

  END PROCEDURE get_logRho0_1


  MODULE PROCEDURE get_logRho0_2

    !************************************************
    !
    !# Returns the value of [[logRho0_2]], the base 10
    !  logarithm of the mass density where the
    !  gamma1_2 polytrope starts, for NS 2 with
    !  piecewise polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_logRho0_2= THIS% logRho0_2

  END PROCEDURE get_logRho0_2


  MODULE PROCEDURE get_logRho1_1

    !************************************************
    !
    !# Returns the value of [[logRho1_1]], the base 10
    !  logarithm of the mass density where the
    !  gamma2_1 polytrope starts, for NS 1 with
    !  piecewise polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_logRho1_1= THIS% logRho1_1

  END PROCEDURE get_logRho1_1


  MODULE PROCEDURE get_logRho1_2

    !************************************************
    !
    !# Returns the value of [[logRho1_2]], the base 10
    !  logarithm of the mass density where the
    !  gamma2_2 polytrope starts, for NS 2 with
    !  piecewise polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_logRho1_2= THIS% logRho1_2

  END PROCEDURE get_logRho1_2


  MODULE PROCEDURE get_logRho2_1

    !************************************************
    !
    !# Returns the value of [[logRho2_1]], the base 10
    !  logarithm of the mass density where the
    !  gamma3_1 polytrope starts, for NS 1 with
    !  piecewise polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_logRho2_1= THIS% logRho2_1

  END PROCEDURE get_logRho2_1


  MODULE PROCEDURE get_logRho2_2

    !************************************************
    !
    !# Returns the value of [[logRho2_2]]]], the base 10
    !  logarithm of the mass density where the
    !  gamma3_2 polytrope starts, for NS 2 with
    !  piecewise polytropic EOS
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_logRho2_2= THIS% logRho2_2

  END PROCEDURE get_logRho2_2


END SUBMODULE bns_access
