! File:         write_par_eos.f90
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

PROGRAM write_par_eos

  !*****************************************************
  !
  !# Write the \(\mathrm{par\_eos.d}\) parameter file
  !  to specify the |eos| when computing |bns| or |drs|
  !  with |lorene|
  !
  !  FT
  !
  !*****************************************************

  USE utility,    ONLY: density_si2cu, k_lorene2cu, k_lorene2cu_pwp, one
  USE constants,  ONLY: kg2g, m2cm
  USE pwp_EOS,    ONLY: Gamma0, K0, select_EOS_parameters, gen_pwp_eos, &
                        get_rho_0, get_rho_1, get_rho_2, &
                        get_Gamma1, get_Gamma2, get_Gamma3, &
                        get_K1, get_K2, get_K3

  IMPLICIT NONE

  INTEGER, PARAMETER:: npoly= 4
  INTEGER:: ios

  !DOUBLE PRECISION:: gamma1_lorene, gamma2_lorene, gamma3_lorene
  DOUBLE PRECISION:: kappa0_lorene, log10_p0_lorene, log10_rho0_lorene, &
                     log10_rho1_lorene, log10_rho2_lorene

  DOUBLE PRECISION:: kappa_poly, gamma_poly

  LOGICAL:: exist
  CHARACTER(LEN= :), ALLOCATABLE:: namefile
  CHARACTER(LEN= :), ALLOCATABLE:: err_msg
  CHARACTER(LEN= 4):: eos_type
  CHARACTER(LEN= 4):: eos

  namefile= "par_eos.d"

  WRITE(*,'(A,/,A)') "** Do you want to produce a parameter file for a polytropic EOS or piecewise polytropic EOS?", "   Please type `poly` for the first, `pwp` for the second: "
  READ(*,'(A)') eos_type

  select_eos_type: SELECT CASE( eos_type )

  CASE( "poly" )

#ifdef __INTEL_COMPILER
  WRITE(*,'("** Please write the DOUBLE PRECISION value of the polytropic exponent gamma: ",\)')
  READ(*,'(F)') gamma_poly
#endif
#ifdef __GFORTRAN__
  PRINT *, "** Please write the DOUBLE PRECISION value of the polytropic constant gamma, with 4 total digits and 3 decimal digits: "
  READ(*,'(F4.3)') gamma_poly
#endif

#ifdef __INTEL_COMPILER
  WRITE(*,'("** Please write the DOUBLE PRECISION value of the polytropic constant K in SPHINCS units: ",\)')
  READ(*,'(F)') kappa_poly
#endif
#ifdef __GFORTRAN__
  PRINT *, "** Please write the DOUBLE PRECISION value of the polytropic constant K in SPHINCS units, with 6 total digits and 2 decimal digits: "
  READ(*,'(F6.2)') kappa_poly
#endif

    kappa_poly= kappa_poly/k_lorene2cu(gamma_poly)

    INQUIRE( FILE= TRIM(namefile), EXIST= exist )

    IF( exist )THEN
        OPEN( UNIT= 2, FILE= TRIM(namefile), STATUS= "REPLACE", &
              FORM= "FORMATTED", &
              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
              IOMSG= err_msg )
    ELSE
        OPEN( UNIT= 2, FILE= TRIM(namefile), STATUS= "NEW", &
              FORM= "FORMATTED", &
              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening " // TRIM(namefile), &
               ". The error message is", err_msg
      STOP
    ENDIF

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = "(A67)" ) &
      "1          Polytropic EOS (cf. documentation of Eos::eos_from_file)"

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = "(A14)" ) &
      "Polytropic EOS"

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = "(F11.9,A18)" ) &
      gamma_poly, "gamma"!, polytropic exponent"

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = "(F19.17,A38)" ) &
      kappa_poly, "kappa [rho_nuc c^2 / n_nuc^gamma]"!, polytropic constant"
    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = "(F11.9,A37)" ) &
      one,        "mean particle mass [m_b]"

    CLOSE( UNIT= 2 )

  CASE( "pwp " )

#ifdef __INTEL_COMPILER
  WRITE(*,'("** Please write an up-to-4-character string containing the name " &
            "of the piecewise polytropic EOS: ",\)')
#endif
#ifdef __GFORTRAN__
  WRITE(*,'("** Please write an up-to-4-character string containing the name " &
            "of the piecewise polytropic EOS: ")')
#endif
  READ(*,'(A)') eos

    CALL select_EOS_parameters(eos)

    kappa0_lorene= K0/k_lorene2cu_pwp(Gamma0)

    log10_p0_lorene= LOG10( K0/k_lorene2cu_pwp(Gamma0) &
                *( get_rho_0()/density_si2cu*kg2g/(m2cm**3.0D0) )**(Gamma0) )

    log10_rho0_lorene= LOG10( get_rho_0()/density_si2cu*kg2g/(m2cm**3.0D0) )

    log10_rho1_lorene= LOG10( get_rho_1()/density_si2cu*kg2g/(m2cm**3.0D0) )

    log10_rho2_lorene= LOG10( get_rho_2()/density_si2cu*kg2g/(m2cm**3.0D0) )

    INQUIRE( FILE= TRIM(namefile), EXIST= exist )

    IF( exist )THEN
        OPEN( UNIT= 2, FILE= TRIM(namefile), STATUS= "REPLACE", &
              FORM= "FORMATTED", &
              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
              IOMSG= err_msg )
    ELSE
        OPEN( UNIT= 2, FILE= TRIM(namefile), STATUS= "NEW", &
              FORM= "FORMATTED", &
              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening " // TRIM(namefile), &
               ". The error message is", err_msg
      STOP
    ENDIF

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = "(A69)" ) &
      "110        Multipolytropic EOS (cf. document. of Eos::eos_multi_poly)"

    !WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = "(A16,A4,A4)" ) eos
    !  "Multipolytropic ", eos, " EOS"
    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = "(A4)" ) &
      ADJUSTL(TRIM(eos))

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = "(I1,A52)" ) &
      npoly, "npoly,         number of polytropes"

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = "(F11.9,A43)" ) &
      Gamma0, "gamma0,        crust (here from SLy)"
    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = "(F11.9,A55)" ) &
      get_Gamma1(), "gamma1,        adiabatic indexes (crust to core)"
    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = "(F11.9,A13)" ) &
      get_Gamma2(), "gamma2"
    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = "(F11.9,A13)" ) &
      get_Gamma3(), "gamma3"
    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = "(ES16.9E3,A55)" ) &
      kappa0_lorene, &
      "kappa0,        K for gamma0 polytrope (here from SLy)"
    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = "(F12.9,A73)" ) &
      log10_p0_lorene, &
      "log10(P0/c^2), log(p) between gamma0 and gamma1 (dyne/cm^2/c_cgs^2)"

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = "(F12.9,A58)" ) &
      log10_rho0_lorene, &
      "log10(rho0),   logs of transition densities (g/cm^3)"
    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = "(F12.9,A17)" ) &
      log10_rho1_lorene, "log10(rho1)"
    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = "(F12.9,A17)" ) &
      log10_rho2_lorene, "log10(rho2)"

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = "(F3.1,A77)" ) &
      0.0, "decInc,       array (size npoly-1) of percentages which change"
    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = "(F3.1,A69)" ) &
      0.0, "the transition densities by changing the"
    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = "(F3.1,A73)" ) &
      0.0, "transition enthalpies (set to 0. to disable)"

    CLOSE( UNIT= 2 )

  END SELECT select_eos_type

  PRINT *
  PRINT *, "** Parameter file ", namefile," written."
  PRINT *

END PROGRAM write_par_eos
