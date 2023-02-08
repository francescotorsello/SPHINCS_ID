! File:         module_ejecta_generic_access.f90
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

SUBMODULE (ejecta_generic) access

  !***************************************************
  !
  !# The module contains the implementation of the
  !  methods of TYPE ejecta that allow to
  !  access PRIVATE members.
  !
  !  FT xx.11.2021
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
    !# Returns the baryon mass of the system [\(M_\odot\)]
    !
    !  FT xx.11.2021
    !
    !************************************************

    IMPLICIT NONE

    CALL this% check_i_matter(i_matter)

    get_mass= this% masses(i_matter)

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
    !# Returns the radii of the system [\(L_\odot\)]
    !
    !  FT xx.11.2021
    !
    !************************************************

    IMPLICIT NONE

    CALL this% check_i_matter(i_matter)

    get_radii= this% sizes(i_matter,:)

  END PROCEDURE get_radii


  MODULE PROCEDURE get_center

    !************************************************
    !
    !# Returns the center of the system [\(L_\odot\)]
    !
    !  FT xx.11.2021
    !
    !************************************************

    IMPLICIT NONE

    CALL this% check_i_matter(i_matter)

    get_center= this% centers(i_matter,:)

  END PROCEDURE get_center


  MODULE PROCEDURE get_barycenter

    !************************************************
    !
    !# Returns the barycenter of the system [\(L_\odot\)]
    !
    !  FT xx.11.2021
    !
    !************************************************

    IMPLICIT NONE

    CALL this% check_i_matter(i_matter)

    get_barycenter= this% barycenters(i_matter,:)

  END PROCEDURE get_barycenter


  MODULE PROCEDURE get_eos

    !************************************************
    !
    !# Returns the |eos| name of the system
    !
    !  FT xx.11.2021
    !
    !************************************************

    IMPLICIT NONE

    CALL this% check_i_matter(i_matter)

    get_eos= "APR4"

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
    !  FT xx.11.2021
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
    !  FT xx.11.2021
    !
    !************************************************

    IMPLICIT NONE

    get_kappa= this% kappa

  END PROCEDURE get_kappa


  MODULE PROCEDURE get_npeos

    !************************************************
    !
    !# Returns the identifier of the EOS for the system
    !
    !  FT xx.11.2021
    !
    !************************************************

    IMPLICIT NONE

    get_npeos= this% npeos

  END PROCEDURE get_npeos


  MODULE PROCEDURE get_gamma0

    !************************************************
    !
    !# Returns the value of [[ejecta:gamma0]], the crust's
    !  polytropic index for the system with piecewise
    !  polytropic EOS
    !
    !  FT xx.11.2021
    !
    !************************************************

    IMPLICIT NONE

    get_gamma0= this% gamma0

  END PROCEDURE get_gamma0


  MODULE PROCEDURE get_gamma1

    !************************************************
    !
    !# Returns the value of [[ejecta:gamma1]], the first
    !  polytropic index for the system with piecewise
    !  polytropic EOS
    !
    !  FT xx.11.2021
    !
    !************************************************

    IMPLICIT NONE

    get_gamma1= this% gamma1

  END PROCEDURE get_gamma1


  MODULE PROCEDURE get_gamma2

    !************************************************
    !
    !# Returns the value of [[ejecta:gamma2]], the second
    !  polytropic index for NS 2 with piecewise
    !  polytropic EOS
    !
    !  FT xx.11.2021
    !
    !************************************************

    IMPLICIT NONE

    get_gamma2= this% gamma2

  END PROCEDURE get_gamma2


  MODULE PROCEDURE get_gamma3

    !************************************************
    !
    !# Returns the value of [[ejecta:gamma3]], the third
    !  polytropic index for the system with piecewise
    !  polytropic EOS (innermost index)
    !
    !  FT xx.11.2021
    !
    !************************************************

    IMPLICIT NONE

    get_gamma3= this% gamma3

  END PROCEDURE get_gamma3


  MODULE PROCEDURE get_kappa0

    !************************************************
    !
    !# Returns the value of [[ejecta:kappa0]], the crust's
    !  polytropic constant for the system with piecewise
    !  polytropic EOS
    !
    !  FT xx.11.2021
    !
    !************************************************

    IMPLICIT NONE

    get_kappa0= this% kappa0

  END PROCEDURE get_kappa0


  MODULE PROCEDURE get_kappa1

    !************************************************
    !
    !# Returns the value of [[ejecta:kappa1]], the first
    !  polytropic constant for the system with piecewise
    !  polytropic EOS
    !
    !  FT xx.11.2021
    !
    !************************************************

    IMPLICIT NONE

    get_kappa1= this% kappa1

  END PROCEDURE get_kappa1


  MODULE PROCEDURE get_kappa2

    !************************************************
    !
    !# Returns the value of [[ejecta:kappa2]], the second
    !  polytropic constant for the system with piecewise
    !  polytropic EOS
    !
    !  FT xx.11.2021
    !
    !************************************************

    IMPLICIT NONE

    get_kappa2= this% kappa2

  END PROCEDURE get_kappa2


  MODULE PROCEDURE get_kappa3

    !************************************************
    !
    !# Returns the value of [[ejecta:kappa3]], the third
    !  polytropic constant for the system with piecewise
    !  polytropic EOS
    !
    !  FT xx.11.2021
    !
    !************************************************

    IMPLICIT NONE

    get_kappa3= this% kappa3

  END PROCEDURE get_kappa3


  MODULE PROCEDURE get_logp1

    !************************************************
    !
    !# Returns the value of [[ejecta:logp1]], the base 10
    !  logarithm of the pressure where the gamma1
    !  polytrope starts, for the system with piecewise
    !  polytropic EOS
    !
    !  FT xx.11.2021
    !
    !************************************************

    IMPLICIT NONE

    get_logp1= this% logp1

  END PROCEDURE get_logp1



  MODULE PROCEDURE get_logRho0

    !************************************************
    !
    !# Returns the value of [[ejecta:logRho0]], the base 10
    !  logarithm of the mass density where the
    !  \(\gamma_1\) polytrope starts, for the system with
    !  piecewise polytropic EOS
    !
    !  FT xx.11.2021
    !
    !************************************************

    IMPLICIT NONE

    get_logRho0= this% logRho0

  END PROCEDURE get_logRho0


  MODULE PROCEDURE get_logRho1

    !************************************************
    !
    !# Returns the value of [[ejecta:logRho1]], the base 10
    !  logarithm of the mass density where the
    !  \(\gamma_2\) polytrope starts, for the system with
    !  piecewise polytropic EOS
    !
    !  FT xx.11.2021
    !
    !************************************************

    IMPLICIT NONE

    get_logRho1= this% logRho1

  END PROCEDURE get_logRho1


  MODULE PROCEDURE get_logRho2

    !************************************************
    !
    !# Returns the value of [[ejecta:logRho2]], the base 10
    !  logarithm of the mass density where the
    !  \(\gamma_3\) polytrope starts, for the system with
    !  piecewise polytropic EOS
    !
    !  FT xx.11.2021
    !
    !************************************************

    IMPLICIT NONE

    get_logRho2= this% logRho2

  END PROCEDURE get_logRho2


  MODULE PROCEDURE get_eos_id

    !**************************************************
    !
    !# Returns the |sphincsid| identifier of the |eos|
    !  of the system
    !
    !  FT 19.11.2021
    !
    !**************************************************

    IMPLICIT NONE

    get_eos_id= this% eos_id

  END PROCEDURE get_eos_id


  MODULE PROCEDURE get_eos_parameters

    !**************************************************
    !
    !# Returns the |eos| parameters of the system
    !  @todo extend to handle single polytropes and tabulated
    !
    !  FT 19.11.2021
    !
    !**************************************************

    IMPLICIT NONE

    ! TODO: Implement single polytropes. Only piecewise polytropes are
    !       supported so far.
    eos_params= [ DBLE(this% eos_id), DBLE(this% npeos), &
          this% gamma0, this% gamma1, this% gamma2, this% gamma3, &
          this% kappa0, this% kappa1, this% kappa2, this% kappa3, &
          this% logP1, &
          this% logRho0, this% logRho1, this% logRho2 ]

  END PROCEDURE get_eos_parameters


END SUBMODULE access
