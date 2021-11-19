! File:         module_ejecta_generic_access.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (ejecta_generic) ejecta_generic_access

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

    CALL THIS% check_i_matter(i_matter)

  END PROCEDURE get_mass


  MODULE PROCEDURE get_radii

    !************************************************
    !
    !# Returns the radii of the DRS [\(L_\odot\)]
    !
    !  FT 2.11.2021
    !
    !************************************************

    IMPLICIT NONE

    CALL THIS% check_i_matter(i_matter)

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

    CALL THIS% check_i_matter(i_matter)

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

    CALL THIS% check_i_matter(i_matter)

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

    CALL THIS% check_i_matter(i_matter)

  END PROCEDURE get_eos


  !-----------------!
  !--  FUNCTIONS  --!
  !-----------------!


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

    get_gamma= THIS% gamma

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

    get_kappa= THIS% kappa

  END PROCEDURE get_kappa


  MODULE PROCEDURE get_npeos

    !************************************************
    !
    !# Returns the identifier of the EOS for the DRS
    !
    !  FT 22.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_npeos= THIS% npeos

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

    get_gamma0= THIS% gamma0

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

    get_gamma1= THIS% gamma1

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

    get_gamma2= THIS% gamma2

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

    get_gamma3= THIS% gamma3

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

    get_kappa0= THIS% kappa0

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

    get_kappa1= THIS% kappa1

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

    get_kappa2= THIS% kappa2

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

    get_kappa3= THIS% kappa3

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

    get_logp1= THIS% logp1

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

    get_logRho0= THIS% logRho0

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

    get_logRho1= THIS% logRho1

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

    get_logRho2= THIS% logRho2

  END PROCEDURE get_logRho2


  MODULE PROCEDURE get_eos_ejectaid

    !**************************************************
    !
    !# Returns the |lorene| ID-number of the EOS of the DRS
    !
    !  FT 19.11.2021
    !
    !**************************************************

    IMPLICIT NONE

  END PROCEDURE get_eos_ejectaid


  MODULE PROCEDURE get_eos_parameters

    !**************************************************
    !
    !# Returns the |eos| parameters of the DRS
    !
    !  FT 19.11.2021
    !
    !**************************************************

    IMPLICIT NONE

  END PROCEDURE get_eos_parameters


END SUBMODULE ejecta_generic_access
