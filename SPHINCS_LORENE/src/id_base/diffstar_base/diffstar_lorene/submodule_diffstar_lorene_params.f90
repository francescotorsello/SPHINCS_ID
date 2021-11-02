! File:         submodule_diffstar_params.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (diffstar_lorene) diffstar_lorene_params

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
                         c_light, cm2km, k_lorene2hydrobase_piecewisepolytrope

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
                         lorene2hydrobase, c_light, cm2km

    IMPLICIT NONE

    IF( THIS% angular_momentum == 0.0D0 )THEN

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


END SUBMODULE diffstar_lorene_params
