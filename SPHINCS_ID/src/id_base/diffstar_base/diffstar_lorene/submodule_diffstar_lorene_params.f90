! File:         submodule_diffstar_lorene_params.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (diffstar_lorene) params

  !********************************************
  !
  !# Implementation of the methods of TYPE diffstar
  !  that import from |lorene| the
  !  parameters of the binary system,
  !  and print them to the standard output.
  !
  !  FT 09.07.2021
  !
  !********************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE import_diffstar_params

    !***************************************************
    !
    !# Store the parameters of the binary neutron
    !  stars' |lorene| ID into member variables
    !
    !  FT 5.10.2020
    !
    !***************************************************

    USE, INTRINSIC :: ISO_C_BINDING,  ONLY: C_CHAR, C_NULL_CHAR
    USE constants, ONLY: Msun_geo, km2m, lorene2hydrobase, k_lorene2hydrobase, &
                         k_lorene2hydrobase_piecewisepolytrope

    IMPLICIT NONE

    INTEGER:: i, nchars
    INTEGER, PARAMETER:: str_length= 100

    CHARACTER(KIND= C_CHAR), DIMENSION(str_length):: eos_tmp_c

    PRINT *
    PRINT *, "** Executing the import_diffstar_params subroutine..."

    CALL get_diffstar_params( THIS% diffstar_ptr,                   &
                              THIS% omega_c,                        &
                              THIS% mass,                           &
                              THIS% mass_grav,                      &
                              THIS% angular_momentum,               &
                              THIS% tsw,                            &
                              THIS% grv2,                           &
                              THIS% grv3,                           &
                              THIS% r_circ,                         &
                              THIS% surface_area,                   &
                              THIS% r_mean,                         &
                              THIS% r_eq,                           &
                              THIS% r_eq_pi2,                       &
                              THIS% r_eq_pi,                        &
                              THIS% r_eq_3pi2,                      &
                              THIS% r_pole,                         &
                              THIS% r_ratio,                        &
                              THIS% r_isco,                         &
                              THIS% f_isco,                         &
                              THIS% specific_energy_isco,           &
                              THIS% specific_angular_momentum_isco, &
                              THIS% area_radius,                    &
                              THIS% ent_center,                     &
                              THIS% nbar_center,                    &
                              THIS% rho_center,                     &
                              THIS% energy_density_center,          &
                              THIS% specific_energy_center,         &
                              THIS% pressure_center,                &
                              THIS% redshift_eqf,                   &
                              THIS% redshift_eqb,                   &
                              THIS% redshift_pole,                  &
                              eos_tmp_c,                            &
                              THIS% eos_loreneid,                   &
                              THIS% gamma,                          &
                              THIS% kappa,                          &
                              THIS% npeos,                          &
                              THIS% gamma0,                         &
                              THIS% gamma1,                         &
                              THIS% gamma2,                         &
                              THIS% gamma3,                         &
                              THIS% kappa0,                         &
                              THIS% kappa1,                         &
                              THIS% kappa2,                         &
                              THIS% kappa3,                         &
                              THIS% logP1,                          &
                              THIS% logRho0,                        &
                              THIS% logRho1,                        &
                              THIS% logRho2 )

    ! Convert distances from |lorene| units (km) to SPHINCS units (Msun_geo)
    ! See MODULE constants for the definition of Msun_geo
    THIS% r_circ      = THIS% r_circ/Msun_geo
    THIS% r_mean      = THIS% r_mean/Msun_geo
    THIS% area_radius = THIS% area_radius/Msun_geo
    THIS% r_eq        = THIS% r_eq/Msun_geo
    THIS% r_eq_pi2    = THIS% r_eq_pi2/Msun_geo
    THIS% r_eq_pi     = THIS% r_eq_pi/Msun_geo
    THIS% r_eq_3pi2   = THIS% r_eq_3pi2/Msun_geo
    THIS% r_pole      = THIS% r_pole/Msun_geo
    THIS% r_isco      = THIS% r_isco/Msun_geo
    THIS% surface_area= THIS% surface_area/(Msun_geo**2.0D0)

    THIS% radii(:)= [THIS% r_eq_pi, THIS% r_eq, &
                     THIS% r_eq_pi2, THIS% r_eq_pi2, &
                     THIS% r_pole, THIS% r_pole]

    THIS% center(:)= [0.0D0, 0.0D0, 0.0D0]

    THIS% barycenter(:)= [0.0D0, 0.0D0, 0.0D0]

    ! Convert hydro quantities from |lorene| units to SPHINCS units
    THIS% nbar_center           = THIS% nbar_center*(MSun_geo*km2m)**3
    THIS% rho_center            = THIS% rho_center*lorene2hydrobase
    THIS% energy_density_center = THIS% energy_density_center*lorene2hydrobase
    THIS% pressure_center       = THIS% pressure_center*lorene2hydrobase

    ! Convert polytropic constants from |lorene| units to SPHINCS units
    IF( THIS% eos_loreneid == 1 )THEN ! If the EOS is polytropic

      THIS% kappa= THIS% kappa*k_lorene2hydrobase( THIS% gamma )

    ELSEIF( THIS% gamma0 /= 0 )THEN ! If the EOS is piecewise polytropic

      THIS% kappa0= THIS% kappa0 &
                      *k_lorene2hydrobase_piecewisepolytrope( THIS% gamma0 )
      THIS% kappa1= THIS% kappa1 &
                      *k_lorene2hydrobase_piecewisepolytrope( THIS% gamma1 )
      THIS% kappa2= THIS% kappa2 &
                      *k_lorene2hydrobase_piecewisepolytrope( THIS% gamma2 )
      THIS% kappa3= THIS% kappa3 &
                      *k_lorene2hydrobase_piecewisepolytrope( THIS% gamma3 )

    ELSEIF( THIS% eos_loreneid == 17 .OR. THIS% eos_loreneid == 20 )THEN
    ! If the EOS is tabulated

    ELSE

      PRINT *, "** ERROR in SUBROUTINE import_lorene_id_params!", &
               " The equation of state is unknown! LORENE EOS IDs=", &
               THIS% eos_loreneid, ", ", THIS% eos_loreneid
      STOP

    ENDIF

    ! Convert C++ strings to FORTRAN strings
    i= 1
    DO
      IF( eos_tmp_c(i) == C_NULL_CHAR .OR. i == str_length ) EXIT
      i= i + 1
    ENDDO
    nchars = i - 1

    ALLOCATE( CHARACTER(nchars):: THIS% eos, STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...allocation error for string eos. ", &
                "The error message is ", err_msg
       PRINT *, "The STAT variable is ", ios
       STOP
    ENDIF
    THIS% eos= TRANSFER( eos_tmp_c(1:nchars), THIS% eos )

    CALL print_diffstar_params( THIS )

    PRINT *, "** Subroutine import_diffstar_params executed."
    PRINT *

  END PROCEDURE import_diffstar_params


END SUBMODULE params
