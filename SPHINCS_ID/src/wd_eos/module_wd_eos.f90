! File:         module_wd_eos.f90
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

MODULE wd_eos

  !********************************************
  !
  !# This MODULE implements the Chandrasekhar's
  !  degenerate |eos| for white dwarfs.
  !
  !  See Benz W., Bowers R.L., Cameron A.G.W.,
  !  Press W.H., 1990, APJ, 348, 647.
  !  doi:10.1086/168273
  !  eqs.(2.4)-(2.5)
  !
  !  FT 19.12.2022
  !
  !********************************************


  USE constants,  ONLY: third, pa2barye, c_light2
  USE units,      ONLY: m0c2_cu
  USE utility,    ONLY: zero, one, two, three, four, ten, lorene2hydrobase, &
                        kg2g, m2cm, c_light2_SI

  IMPLICIT NONE


  PRIVATE


  DOUBLE PRECISION, PARAMETER:: mu_e= two
  !! Mean molecular weight per electron

  DOUBLE PRECISION, PARAMETER:: a_wd_cgs= 6.02D+21*pa2barye
  DOUBLE PRECISION, PARAMETER:: b_wd_cgs= mu_e*9.82D+8*kg2g/(m2cm**3)

  DOUBLE PRECISION, PARAMETER:: a_wd= 6.02D+21/c_light2_SI*lorene2hydrobase
  !! Constant with dimensions of a pressure in code units
  DOUBLE PRECISION, PARAMETER:: b_wd= mu_e*9.82D+8*lorene2hydrobase
  !! Constant with dimensions of a density in code units

 ! DOUBLE PRECISION:: a_wd_m0c2= &
 !                               6.02D+21/c_light2_SI*lorene2hydrobase/m0c2_cu
 ! !! Constant with dimensions of a pressure in units of \(m_0 c^2\)
 ! DOUBLE PRECISION:: b_wd_m0c2= mu_e*9.82D+8*lorene2hydrobase/m0c2_cu
 ! !! Constant with dimensions of a density in units of \(m_0 c^2\)


  PUBLIC:: pr_wd, u_wd, rho_wd, test_wd_eos_cgs


  CONTAINS


  !--------------------------!
  !--  PRIVATE PROCEDURES  --!
  !--------------------------!


  DOUBLE PRECISION PURE FUNCTION f_wd(x)

    !********************************************
    !
    !# Function used in the computation of
    !  the degenerate pressure
    !
    !  FT 19.12.2022
    !
    !********************************************

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: x

    f_wd= x*(two*x**2 - three)*SQRT(x**2 + one) + three*ASINH(x)

  END FUNCTION f_wd


  DOUBLE PRECISION PURE FUNCTION g_wd(x)

    !********************************************
    !
    !# Function used in the computation of
    !  the degenerate specific internal energy
    !
    !  FT 19.12.2022
    !
    !********************************************

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: x

    g_wd= two*four*x**3*(SQRT(x**2 + one) - one) - f_wd(x)

  END FUNCTION g_wd


  !-------------------------!
  !--  PUBLIC PROCEDURES  --!
  !-------------------------!


  DOUBLE PRECISION PURE FUNCTION pr_wd(rho)

    !********************************************
    !
    !# Degenerate pressure as a function of density
    !
    !  See Benz W., Bowers R.L., Cameron A.G.W.,
    !  Press W.H., 1990, APJ, 348, 647.
    !  doi:10.1086/168273
    !  eqs.(2.4)
    !
    !
    !  FT 19.12.2022
    !
    !********************************************

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: rho

    DOUBLE PRECISION:: x

    x= (rho/b_wd)**third

    pr_wd= a_wd*f_wd(x)

  END FUNCTION pr_wd


  DOUBLE PRECISION PURE FUNCTION pr_wd_cgs(rho)

    !********************************************
    !
    !# Degenerate pressure as a function of density,
    !  all in CGS units
    !
    !  See Benz W., Bowers R.L., Cameron A.G.W.,
    !  Press W.H., 1990, APJ, 348, 647.
    !  doi:10.1086/168273
    !  eqs.(2.4)
    !
    !
    !  FT 19.12.2022
    !
    !********************************************

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: rho

    DOUBLE PRECISION:: x

    x= (rho/b_wd_cgs)**third

    pr_wd_cgs= a_wd_cgs*f_wd(x)

  END FUNCTION pr_wd_cgs


  DOUBLE PRECISION PURE FUNCTION u_wd(rho)

    !********************************************
    !
    !# Degenerate specific internal energy as a
    !  function of density
    !
    !  See Benz W., Bowers R.L., Cameron A.G.W.,
    !  Press W.H., 1990, APJ, 348, 647.
    !  doi:10.1086/168273
    !  eqs.(2.5)
    !
    !  FT 19.12.2022
    !
    !********************************************

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: rho

    DOUBLE PRECISION:: x

    x= (rho/b_wd)**third

    u_wd= a_wd*g_wd(x)/rho

  END FUNCTION u_wd


  DOUBLE PRECISION FUNCTION rho_wd(pr, rho_left_bracket, rho_right_bracket)

    !********************************************
    !
    !# Degenerate density as a function of pressure
    !
    !  Use bisection to find the value of the density,
    !  given the pressure
    !
    !  FT 19.12.2022
    !
    !********************************************

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: pr
    DOUBLE PRECISION, INTENT(IN):: rho_left_bracket
    DOUBLE PRECISION, INTENT(IN):: rho_right_bracket

    INTEGER,          PARAMETER:: tolerance_magnitude= 3
    DOUBLE PRECISION, PARAMETER:: tolerance_pr       = 1.D-6
    DOUBLE PRECISION, PARAMETER:: tolerance_rho      = 1.D-10

    INTEGER:: cnt
    DOUBLE PRECISION:: rho_left, rho_right, rho_mean, &
                       pr_left, pr_right, pr_mean, pr_cgs

    IF(pr == zero)THEN
      rho_wd= zero
      RETURN
    ENDIF

    rho_left = FLOOR(LOG10(rho_left_bracket))
    rho_right= CEILING(LOG10(rho_right_bracket))
    ! Bisection in logarithmic scale, to find the approximate order of magnitude
    DO

      pr_left = pr_wd(ten**rho_left)
      pr_right= pr_wd(ten**rho_right)
      IF( pr_left <= pr .AND. pr_right > pr )THEN

        rho_mean= FLOOR((rho_left + rho_right)/two)
        pr_mean = pr_wd(ten**rho_mean)
        !PRINT *, FLOOR(ABS(LOG10(pr_right) - LOG10(pr_left)))
        !PRINT *, LOG10(pr_right)
        !PRINT *, LOG10(pr_left)
        IF( FLOOR(ABS(LOG10(pr_right) - LOG10(pr_left))) &
            <= tolerance_magnitude )THEN
          EXIT
        ENDIF

        IF(pr_mean < pr)THEN
          rho_left = rho_mean
        ELSE
          rho_right= rho_mean
        ENDIF

      ELSEIF( pr_left < pr .AND. pr_right < pr )THEN

        rho_right= rho_right + two

      ELSEIF( pr_left > pr .AND. pr_right > pr )THEN

        rho_left= rho_left - two

      ELSE

        PRINT *
        PRINT *, "** ERROR in SUBROUTINE rho_wd in MODULE wd_eos!"
        PRINT *, "rho_left=", rho_left
        PRINT *, "rho_right=", rho_right
        PRINT *, "pr_left=", pr_left
        PRINT *, "pr=", pr
        PRINT *, "pr_right=", pr_right
        PRINT *
        STOP

      ENDIF

    ENDDO
    ! Bisection in linear scale, to find the precise value
    !rho_left = (ten**rho_left )/lorene2hydrobase*kg2g/(m2cm**3)
    !rho_right= (ten**rho_right)/lorene2hydrobase*kg2g/(m2cm**3)
    pr_cgs   = pr!/lorene2hydrobase*kg2g/(m2cm**3)*c_light2
    rho_left = (ten**rho_left )
    rho_right= (ten**rho_right)
    !PRINT *, "rho_left =", rho_left
    !PRINT *, "rho_right=", rho_right
    !PRINT *, "pr_cgs   =", pr_cgs
    cnt= 0
    DO

      !pr_left = pr_wd_cgs(rho_left)
      !pr_right= pr_wd_cgs(rho_right)
      pr_left = pr_wd(rho_left)
      pr_right= pr_wd(rho_right)
      IF( pr_left <= pr_cgs .AND. pr_right > pr_cgs )THEN

        rho_mean= (rho_left + rho_right)/two
        !pr_mean = pr_wd_cgs(rho_mean)
        pr_mean = pr_wd(rho_mean)
        !PRINT *, "pr_left  =", pr_left
        !PRINT *, "pr_right =", pr_right
        !PRINT *, "rho_mean =", rho_mean
        !PRINT *, "pr_mean  =", pr_mean
        !STOP
        !IF(cnt > 100) PRINT *, "ABS((pr_left - pr_right)/pr_right)=", &
        !                       ABS((pr_left - pr_right)/pr_right)
        !IF(cnt > 100) PRINT *, "pr_left  =", pr_left
        !IF(cnt > 100) PRINT *, "pr_right =", pr_right
        !IF(cnt > 100) PRINT *, "pr_mean  =", pr_mean
        !IF(cnt > 100) PRINT *, "rho_left =", rho_left
        !IF(cnt > 100) PRINT *, "rho_right=", rho_right
        !IF(cnt > 100) PRINT *, "rho_mean =", rho_mean
        !IF(cnt > 100) PRINT *, "pr_cgs   =", pr_cgs
        IF( ABS((pr_left - pr_right)/pr_right) < tolerance_pr &
            .OR. &
            ABS((rho_left - rho_right)/rho_right) < tolerance_rho )THEN
          EXIT
        ENDIF

        IF(pr_mean < pr_cgs)THEN
          rho_left = rho_mean
        ELSE
          rho_right= rho_mean
        ENDIF

      ELSEIF( pr_left < pr_cgs .AND. pr_right < pr_cgs )THEN

        rho_right= two*rho_right

      ELSEIF( pr_left > pr_cgs .AND. pr_right > pr_cgs )THEN

        rho_left= rho_left/two

      ELSE

        PRINT *
        PRINT *, "** ERROR in SUBROUTINE rho_wd in MODULE wd_eos!"
        PRINT *, "rho_left=", rho_left
        PRINT *, "rho_right=", rho_right
        PRINT *, "pr_left=", pr_left
        PRINT *, "pr_cgs=", pr_cgs
        PRINT *, "pr_right=", pr_right
        PRINT *
        STOP

      ENDIF

      cnt= cnt + 1
      IF(cnt > 150)THEN
        STOP
      ENDIF

    ENDDO

    !rho_wd= rho_mean*lorene2hydrobase/kg2g*(m2cm**3)
    rho_wd= rho_mean
    !PRINT *, "completed"
    !PRINT *

  END FUNCTION rho_wd


  SUBROUTINE test_wd_eos_cgs(rho_input, rho, pr, u)

    !********************************************
    !
    !# Test that the implementation is correct.
    !  Takes the density in code units as input,
    !  and assigns CGS values to the pressure,
    !  and the density recomputed from the pressure.
    !  It also computed the dimensionless
    !  specific internal energy.
    !
    !  FT 20.12.2022
    !
    !********************************************

    USE utility,   ONLY: lorene2hydrobase, kg2g, m2cm

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: rho_input
    DOUBLE PRECISION, INTENT(OUT):: rho, pr, u

    PRINT *, "** Testing that the implementation of the Chandrasekhar EOS", &
             " for white dwarfs is correct."
    PRINT *

    PRINT *, "   a_wd in code units=", a_wd
    PRINT *, "   b_wd in code units=", b_wd
    PRINT *, "   a_wd in CGS units=", a_wd_cgs
    PRINT *, "   b_wd in CGS units=", b_wd_cgs
    PRINT *
    PRINT *, "   rho_input in code units=", rho_input
    PRINT *, "   rho_input in CGS units=", &
             rho_input/lorene2hydrobase*kg2g/(m2cm**3)
    PRINT *

    pr = pr_wd(rho_input)
    u  = u_wd(rho_input)
    rho= rho_wd(pr,zero,one)

    PRINT *, "   pr in code units=", pr
    pr= pr/lorene2hydrobase*kg2g/(m2cm**3)*c_light2
    PRINT *, "   pr in CGS units=", pr
    PRINT *

    PRINT *, "   u (dimensionless)=", u
    PRINT *

    PRINT *, "   Recomputed rho in code units=", rho
    rho= rho/lorene2hydrobase*kg2g/(m2cm**3)
    PRINT *, "   Recomputed rho in CGS units=", rho
    PRINT *

  END SUBROUTINE test_wd_eos_cgs


END MODULE wd_eos
