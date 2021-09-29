! File:         submodule_bns_params.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (bns_lorene) bns_params

  !********************************************
  !
  !# Implementation of the methods of TYPE bns
  !  that import from LORENE the
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


  MODULE PROCEDURE import_id_params

    !***************************************************
    !
    !# Store the parameters of the binary neutron
    !  stars' LORENE ID into member variables
    !
    !  FT 5.10.2020
    !
    !***************************************************

    USE, INTRINSIC :: ISO_C_BINDING,  ONLY: C_CHAR

    USE constants, ONLY: Msun_geo, km2m, lorene2hydrobase, k_lorene2hydrobase, &
                         c_light, cm2km, k_lorene2hydrobase_piecewisepolytrope

    IMPLICIT NONE

    INTEGER:: i, nchars
    INTEGER, PARAMETER:: str_length= 100

    CHARACTER(KIND= C_CHAR), DIMENSION(str_length):: eos1_tmp_c
    CHARACTER(KIND= C_CHAR), DIMENSION(str_length):: eos2_tmp_c
    !CHARACTER, DIMENSION(:), ALLOCATABLE:: eos1_tmp
    !CHARACTER, DIMENSION(:), ALLOCATABLE:: eos2_tmp

    PRINT *, "** Executing the import_lorene_id_params subroutine..."

    CALL get_lorene_id_params( THIS% bns_ptr, &
                               THIS% angular_vel, &
                               THIS% distance, &
                               THIS% distance_com, &
                               THIS% mass1, &
                               THIS% mass2, &
                               THIS% mass_grav1, &
                               THIS% mass_grav2, &
                               THIS% adm_mass, &
                               THIS% angular_momentum, &
                               THIS% area_radius1, &
                               THIS% radius1_x_comp, &
                               THIS% radius1_y, &
                               THIS% radius1_z, &
                               THIS% radius1_x_opp, &
                               THIS% center1_x, &
                               THIS% barycenter1_x, &
                               THIS% area_radius2, &
                               THIS% radius2_x_comp, &
                               THIS% radius2_y, &
                               THIS% radius2_z, &
                               THIS% radius2_x_opp, &
                               THIS% center2_x, &
                               THIS% barycenter2_x, &
                               THIS% ent_center1, &
                               THIS% nbar_center1, &
                               THIS% rho_center1, &
                               THIS% energy_density_center1, &
                               THIS% specific_energy_center1, &
                               THIS% pressure_center1, &
                               THIS% ent_center2, &
                               THIS% nbar_center2, &
                               THIS% rho_center2, &
                               THIS% energy_density_center2, &
                               THIS% specific_energy_center2, &
                               THIS% pressure_center2, &
                               eos1_tmp_c, &
                               eos2_tmp_c, &
                               THIS% eos1_loreneid, &
                               THIS% eos2_loreneid, &
                               THIS% gamma_1, &
                               THIS% kappa_1, &
                               THIS% gamma_2, &
                               THIS% kappa_2, &
                               THIS% npeos_1, &
                               THIS% gamma0_1, &
                               THIS% gamma1_1, &
                               THIS% gamma2_1, &
                               THIS% gamma3_1, &
                               THIS% kappa0_1, &
                               THIS% kappa1_1, &
                               THIS% kappa2_1, &
                               THIS% kappa3_1, &
                               THIS% logP1_1,  &
                               THIS% logRho0_1,&
                               THIS% logRho1_1,&
                               THIS% logRho2_1,&
                               THIS% npeos_2,  &
                               THIS% gamma0_2, &
                               THIS% gamma1_2, &
                               THIS% gamma2_2, &
                               THIS% gamma3_2, &
                               THIS% kappa0_2, &
                               THIS% kappa1_2, &
                               THIS% kappa2_2, &
                               THIS% kappa3_2, &
                               THIS% logP1_2,  &
                               THIS% logRho0_2,&
                               THIS% logRho1_2,&
                               THIS% logRho2_2)

    ! Convert distances from LORENE units (km) to SPHINCS units (Msun_geo)
    ! See MODULE constants for the definition of Msun_geo
    THIS% distance      = THIS% distance/Msun_geo
    THIS% distance_com  = THIS% distance_com/Msun_geo
    THIS% area_radius1  = THIS% area_radius1/Msun_geo
    THIS% radius1_x_comp= THIS% radius1_x_comp/Msun_geo
    THIS% radius1_y     = THIS% radius1_y/Msun_geo
    THIS% radius1_z     = THIS% radius1_z/Msun_geo
    THIS% radius1_x_opp = THIS% radius1_x_opp/Msun_geo
    THIS% center1_x     = THIS% center1_x/Msun_geo
    THIS% barycenter1_x = THIS% barycenter1_x/Msun_geo
    THIS% area_radius2  = THIS% area_radius2/Msun_geo
    THIS% radius2_x_comp= THIS% radius2_x_comp/Msun_geo
    THIS% radius2_y     = THIS% radius2_y/Msun_geo
    THIS% radius2_z     = THIS% radius2_z/Msun_geo
    THIS% radius2_x_opp = THIS% radius2_x_opp/Msun_geo
    THIS% center2_x     = THIS% center2_x/Msun_geo
    THIS% barycenter2_x = THIS% barycenter2_x/Msun_geo

    ! Convert hydro quantities from LORENE units to SPHINCS units
    THIS% nbar_center1           = THIS% nbar_center1*(MSun_geo*km2m)**3
    THIS% rho_center1            = THIS% rho_center1*lorene2hydrobase
    THIS% energy_density_center1 = THIS% energy_density_center1*lorene2hydrobase
    THIS% pressure_center1       = THIS% pressure_center1*lorene2hydrobase
    THIS% nbar_center2           = THIS% nbar_center2*(MSun_geo*km2m)**3
    THIS% rho_center2            = THIS% rho_center2*lorene2hydrobase
    THIS% energy_density_center2 = THIS% energy_density_center2*lorene2hydrobase
    THIS% pressure_center2       = THIS% pressure_center2*lorene2hydrobase

    ! Convert polytropic constants from LORENE units to SPHINCS units
    IF( THIS% eos1_loreneid == 1 )THEN ! If the EOS is polytropic

      THIS% kappa_1= THIS% kappa_1*k_lorene2hydrobase( THIS% gamma_1 )
      THIS% kappa_2= THIS% kappa_2*k_lorene2hydrobase( THIS% gamma_2 )

    ELSEIF( THIS% gamma0_1 /= 0 )THEN ! If the EOS is piecewise polytropic

      THIS% kappa0_1= THIS% kappa0_1 &
                      *k_lorene2hydrobase_piecewisepolytrope( THIS% gamma0_1 )
      THIS% kappa1_1= THIS% kappa1_2 &
                      *k_lorene2hydrobase_piecewisepolytrope( THIS% gamma1_1 )
      THIS% kappa2_1= THIS% kappa2_1 &
                      *k_lorene2hydrobase_piecewisepolytrope( THIS% gamma2_1 )
      THIS% kappa3_1= THIS% kappa3_2 &
                      *k_lorene2hydrobase_piecewisepolytrope( THIS% gamma3_1 )
      THIS% kappa0_2= THIS% kappa0_1 &
                      *k_lorene2hydrobase_piecewisepolytrope( THIS% gamma0_2 )
      THIS% kappa1_2= THIS% kappa1_2 &
                      *k_lorene2hydrobase_piecewisepolytrope( THIS% gamma1_2 )
      THIS% kappa2_2= THIS% kappa2_1 &
                      *k_lorene2hydrobase_piecewisepolytrope( THIS% gamma2_2 )
      THIS% kappa3_2= THIS% kappa3_2 &
                      *k_lorene2hydrobase_piecewisepolytrope( THIS% gamma3_2 )

    ELSEIF( THIS% eos1_loreneid == 17 )THEN ! If the EOS is tabulated

    ELSE

      PRINT *, "** ERROR in SUBROUTINE import_lorene_id_params!", &
               " The equation of state is unknown! LORENE EOS IDs=", &
               THIS% eos1_loreneid, ", ", THIS% eos2_loreneid
      STOP

    ENDIF

    ! Compute mOmega

    THIS% mOmega= THIS% angular_vel/(c_light*cm2km) &
                  *(THIS% mass_grav1 + THIS% mass_grav2)*Msun_geo

    ! Compute t_merger

    THIS% t_merger= 5.0D0/256.0D0*(THIS% distance**4.0D0) &
                    /( THIS% mass_grav1*THIS% mass_grav2* &
                       ( THIS% mass_grav1 + THIS% mass_grav2 ) )

    ! Convert C++ strings to FORTRAN strings
    i= 1
    DO
      IF( eos1_tmp_c(i) == C_NULL_CHAR .OR. i == str_length ) EXIT
      i= i + 1
    ENDDO
    nchars = i - 1

    !ALLOCATE( eos1_tmp( nchars ), STAT= ios, ERRMSG= err_msg )
    !IF( ios > 0 )THEN
    !   PRINT *, "...allocation error for array eos1_tmp. ", &
    !            "The error message is ", err_msg
    !   PRINT *, "The STAT variable is ", ios
    !   STOP
    !ENDIF
    !eos1_tmp = TRANSFER( eos1_tmp_c(1:nchars), eos1_tmp )

    ALLOCATE( CHARACTER(nchars):: THIS% eos1, STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...allocation error for string eos1. ", &
                "The error message is ", err_msg
       PRINT *, "The STAT variable is ", ios
       STOP
    ENDIF
    THIS% eos1= TRANSFER( eos1_tmp_c(1:nchars), THIS% eos1 )
    !DO i= 1, nchars, 1
    !  THIS% eos1(i:i)= eos1_tmp(i)
    !ENDDO

    i= 1
    DO
      IF( eos2_tmp_c(i) == C_NULL_CHAR .OR. i == str_length ) EXIT
      i= i + 1
    ENDDO
    nchars = i - 1

    !ALLOCATE( eos2_tmp( nchars ), STAT= ios, ERRMSG= err_msg )
    !IF( ios > 0 )THEN
    !   PRINT *, "...allocation error for array eos2_tmp. ", &
    !            "The error message is ", err_msg
    !   PRINT *, "The STAT variable is ", ios
    !   STOP
    !ENDIF
    !eos2_tmp = TRANSFER( eos2_tmp_c(1:nchars), eos2_tmp )

    ALLOCATE( CHARACTER(nchars):: THIS% eos2, STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...allocation error for string eos2. ", &
                "The error message is ", err_msg
       PRINT *, "The STAT variable is ", ios
       STOP
    ENDIF
    THIS% eos2= TRANSFER( eos2_tmp_c(1:nchars), THIS% eos2 )
    !DO i= 1, nchars, 1
    !  THIS% eos2(i:i)= eos2_tmp(i)
    !ENDDO

    CALL print_id_params( THIS )

    PRINT *, "** Subroutine import_lorene_id_params executed."
    PRINT *

  END PROCEDURE import_id_params


  MODULE PROCEDURE print_id_params

    !****************************************************
    !
    !# Print the parameters of the binary neutron
    !  stars' initial data computed by LORENE
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
          "Call the SUBROUTINE import_lorene_id_params to read them."
      PRINT *

    ELSE

      PRINT *
      PRINT *, " ** The parameters of the binary system are:"
      PRINT *
      PRINT *, " Distance between the points of highest density = ",&
               THIS% distance, " M_sun^geo = ", THIS% distance*Msun_geo, " km"
      PRINT *, " Distance between the barycenters = ", &
               THIS% distance_com, " M_sun^geo", THIS% distance_com*Msun_geo, &
               " km"
      PRINT *
      PRINT *, " Baryonic mass of NS 1 = ", THIS% mass1, " M_sun"
      PRINT *, " Baryonic mass of NS 2 = ", THIS% mass2, " M_sun"
      PRINT *, " Gravitational mass of NS 1 = ", THIS% mass_grav1, " M_sun"
      PRINT *, " Gravitational mass of NS 2 = ", THIS% mass_grav2, " M_sun"
      PRINT *, " ADM mass = ", THIS% adm_mass, " M_sun"
      PRINT *
      PRINT *, " Stellar center of NS 1 = ", THIS% center1_x, " M_sun^geo"
      PRINT *, " Stellar center of NS 2 = ", THIS% center2_x, " M_sun^geo"
      PRINT *, " Barycenter of NS 1 = ", THIS% barycenter1_x, " M_sun^geo"
      PRINT *, " Barycenter of NS 2 = ", THIS% barycenter2_x, " M_sun^geo"
      PRINT *, " Angular velocity Omega_0 = ", THIS% angular_vel, " rad/s = ", &
               THIS% angular_vel/(c_light*cm2km), "km^{-1}"
      PRINT *, " mOmega = ", &
               "Omega_0[km^{-1}]*(mass_grav1[km] + mass_grav2[km]) = ",&
               THIS% mOmega, "[pure number]"
      PRINT *, " Angular momentum of the system = ", &
               THIS% angular_momentum, " G M_sun^2 /c"
      PRINT *, " Estimated time of the merger t_merger = ", THIS% t_merger, &
               " M_sun^geo = ", THIS% t_merger*MSun_geo/(c_light*cm2km)*1000.0,&
               " ms, from Peters_PR_136_B1224_1964, eq. (5.10)"
      PRINT *
      PRINT *, " Estimated separation to have the merger at t_merger = 2000", &
               " Msun_geo = ", 2000.0D0*MSun_geo/(c_light*cm2km)*1000.0, &
               " ms :", &
               ( 2000.0D0*( THIS% mass_grav1*THIS% mass_grav2* &
                  ( THIS% mass_grav1 + THIS% mass_grav2 ) )/(5.0D0/256.0D0) ) &
                **(1.0D0/4.0D0), "M_sun^geo = ", &
                ( 2000.0D0*( THIS% mass_grav1*THIS% mass_grav2* &
                ( THIS% mass_grav1 + THIS% mass_grav2 ) )/(5.0D0/256.0D0) ) &
                **(1.0D0/4.0D0)*Msun_geo, &
                "km, from Peters_PR_136_B1224_1964, eq. (5.10)"
      PRINT *
      PRINT *, " Radii of star 1: "
      PRINT *, "  Areal (or circumferential) radius for the star in the", &
               "  binary system [the one used in the", &
               "  (gravitational)mass-(areal)radius diagrams",&
               "  is for a TOV star], x direction:", &
               THIS% area_radius1, " M_sun^geo = ", &
               THIS% area_radius1*Msun_geo, " km"
      PRINT *, "  x direction, towards companion = ", &
               THIS% radius1_x_comp, " M_sun^geo"
      PRINT *, "  x direction, opposite to companion = ", &
               THIS% radius1_x_opp, " M_sun^geo"
      PRINT *, "  y direction = ", THIS% radius1_y, " M_sun^geo"
      PRINT *, "  z direction = ", THIS% radius1_z, " M_sun^geo"
      PRINT *, " Radii of star 2 :"
      PRINT *, "  Areal (or circumferential) radius for the star in the", &
               "  binary system [the one used in the", &
               "  (gravitational)mass-(areal)radius diagrams",&
               "  is for a TOV star], x direction:", &
               THIS% area_radius2, " M_sun^geo", &
               THIS% area_radius2*Msun_geo, " km"
      PRINT *, "  x direction, towards companion = ", &
               THIS% radius2_x_comp, " M_sun^geo"
      PRINT *, "  x direction, opposite to companion = ", &
               THIS% radius2_x_opp, " M_sun^geo"
      PRINT *, "  y direction = ", THIS% radius2_y, " M_sun^geo"
      PRINT *, "  z direction = ", THIS% radius2_z, " M_sun^geo"
      PRINT *
      PRINT *, " Hydro quantities at the center of star 1: "
      PRINT *, "  Central enthalpy = ", THIS% ent_center1, " c^2"
      PRINT *, "  Central baryon number density = ", THIS% nbar_center1, &
               " (M_sun^geo)^{-3} =", &
               THIS% nbar_center1/(MSun_geo*km2m*m2cm)**3, "cm^{-3}"
      PRINT *, "  Central baryon mass density = ", THIS% rho_center1, &
               " M_sun^geo (M_sun^geo)^{-3} =", &
               THIS% rho_center1/lorene2hydrobase*kg2g/(m2cm**3), "g cm^{-3}"
      PRINT *, "  Central energy density = ", THIS% energy_density_center1, &
               " M_sun^geo c^2 (M_sun^geo)^{-3}", &
               THIS% energy_density_center1/lorene2hydrobase*kg2g/(m2cm**3), &
               "g c^2 cm^{-3}"
      PRINT *, "  Central specific energy = ", THIS% specific_energy_center1, &
               " c^2"
      PRINT *, "  Central pressure = ", THIS% pressure_center1, &
               " M_sun^geo c^2 (M_sun^geo)^{-3}", &
               THIS% pressure_center1/lorene2hydrobase*kg2g/(m2cm**3), &
               "g c^2 cm^{-3}"
      PRINT *, " Hydro quantities at the center of star 2: "
      PRINT *, "  Central enthalpy = ", THIS% ent_center2, " c^2"
      PRINT *, "  Central baryon number density = ", THIS% nbar_center2, &
               " (M_sun^geo)^{-3} =", &
               THIS% nbar_center2/(MSun_geo*km2m*m2cm)**3, "cm^{-3}"
      PRINT *, "  Central baryon mass density = ", THIS% rho_center2, &
               " M_sun^geo (M_sun^geo)^{-3} =", &
               THIS% rho_center2/lorene2hydrobase*kg2g/(m2cm**3), "g cm^{-3}"
      PRINT *, "  Central energy density = ", THIS% energy_density_center2, &
               " M_sun^geo c^2 (M_sun^geo)^{-3}", &
               THIS% energy_density_center2/lorene2hydrobase*kg2g/(m2cm**3), &
               "g c^2 cm^{-3}"
      PRINT *, "  Central specific energy = ", THIS% specific_energy_center2, &
               " c^2"
      PRINT *, "  Central pressure = ", THIS% pressure_center2, &
               " M_sun^geo c^2 (M_sun^geo)^{-3}", &
               THIS% pressure_center2/lorene2hydrobase*kg2g/(m2cm**3), &
               "g c^2 cm^{-3}"
      PRINT *
      !IF( show_progress ) &
        PRINT *, " Equations of state for star 1 (EOS1) = ", TRIM(THIS% eos1)
      !IF( show_progress ) &
        PRINT *, " Equations of state for star 2 (EOS2) = ", TRIM(THIS% eos2)
      !IF( show_progress ) PRINT *

      IF( THIS% gamma0_1 == 0 )THEN ! If the EOS is polytropic

        PRINT *, " Parameters for EOS1: "
        PRINT *, "  Polytopic index gamma_1 = ", THIS% gamma_1
        PRINT *, "  Pressure coefficient = ",&
                 THIS% kappa_1/k_lorene2hydrobase( THIS% gamma_1 ), &
                 "rho_nuc c^2 / n_nuc^gamma_1 = ", THIS% kappa_1, &
                 "[pure number]"
        PRINT *, " Parameters for EOS2: "
        PRINT *, "  Polytopic index gamma_2 = ", THIS% gamma_2
        PRINT *, "  Pressure coefficient = ",&
                 THIS% kappa_2/k_lorene2hydrobase( THIS% gamma_2 ), &
                 "rho_nuc c^2 / n_nuc^gamma_2 = ", THIS% kappa_2, &
                 "[pure number]"
        PRINT *

      ELSEIF( THIS% gamma0_1 /= 0 )THEN ! If the EOS is piecewise polytropic

        PRINT *, " Parameters for EOS1: "
        PRINT *, "  Number of polytropic indexes = ", THIS% npeos_1
        PRINT *, "  Polytopic index gamma0_1 = ", THIS% gamma0_1
        PRINT *, "  Polytopic index gamma1_1 = ", THIS% gamma1_1
        PRINT *, "  Polytopic index gamma2_1 = ", THIS% gamma2_1
        PRINT *, "  Polytopic index gamma3_1 = ", THIS% gamma3_1
        PRINT *, "  Pressure coefficient for the crust (here from SLy) = ",&
                 THIS% kappa0_1/k_lorene2hydrobase( THIS% gamma0_1 ), &
                 "rho_nuc c^2 / n_nuc^gamma0_1 = ", THIS% kappa0_1, &
                 "[pure number]"
        PRINT *, "  Pressure coefficient for the first polytrope = ",&
                 THIS% kappa1_1/k_lorene2hydrobase( THIS% gamma1_1 ), &
                 "rho_nuc c^2 / n_nuc^gamma1_1", THIS% kappa1_1, &
                 "[pure number]"
        PRINT *, "  Pressure coefficient for the second polytrope = ",&
                 THIS% kappa2_1/k_lorene2hydrobase( THIS% gamma2_1 ), &
                 "rho_nuc c^2 / n_nuc^gamma2_1", THIS% kappa2_1, &
                 "[pure number]"
        PRINT *, "  Pressure coefficient for the third polytrope = ",&
                 THIS% kappa3_1/k_lorene2hydrobase( THIS% gamma3_1 ), &
                 "rho_nuc c^2 / n_nuc^gamma3_1", THIS% kappa3_1, &
                 "[pure number]"
        PRINT *, "  Base 10 exponent of the pressure at the first fiducial " &
                 // "density (between gamma_0 and gamma_1) (dyne/cm^2)= ", &
                 THIS% logP1_1
        PRINT *, "  Base 10 exponent of first fiducial density (g/cm^3) = ", &
                 THIS% logRho0_1
        PRINT *, "  Base 10 exponent of second fiducial density (g/cm^3) = ",&
                 THIS% logRho1_1
        PRINT *, "  Base 10 exponent of third fiducial density (g/cm^3) = ", &
                 THIS% logRho2_1
        PRINT *
        PRINT *, " Parameters for EOS2: "
        PRINT *, "  Number of polytropic indexes = ", THIS% npeos_2
        PRINT *, "  Polytopic index gamma0_2 = ", THIS% gamma0_2
        PRINT *, "  Polytopic index gamma1_2 = ", THIS% gamma1_2
        PRINT *, "  Polytopic index gamma2_2 = ", THIS% gamma2_2
        PRINT *, "  Polytopic index gamma3_2 = ", THIS% gamma3_2
        PRINT *, "  Pressure coefficient for the crust (here from SLy) = ",&
                 THIS% kappa0_2/k_lorene2hydrobase( THIS% gamma0_2 ), &
                 "rho_nuc c^2 / n_nuc^gamma0_2 = ", THIS% kappa0_2, &
                 "[pure number]"
        PRINT *, "  Pressure coefficient for the first polytrope = ",&
                 THIS% kappa1_2/k_lorene2hydrobase( THIS% gamma1_2 ), &
                 "rho_nuc c^2 / n_nuc^gamma1_2", THIS% kappa1_2, &
                 "[pure number]"
        PRINT *, "  Pressure coefficient for the second polytrope = ",&
                 THIS% kappa2_2/k_lorene2hydrobase( THIS% gamma2_2 ), &
                 "rho_nuc c^2 / n_nuc^gamma2_2", THIS% kappa2_2, &
                 "[pure number]"
        PRINT *, "  Pressure coefficient for the third polytrope = ",&
                 THIS% kappa3_2/k_lorene2hydrobase( THIS% gamma3_2 ), &
                 "rho_nuc c^2 / n_nuc^gamma3_2", THIS% kappa3_2, &
                 "[pure number]"
        PRINT *, "  Base 10 exponent of the pressure at the first fiducial " &
                 // "density (between gamma_0 and gamma_1) (dyne/cm^2)= ", &
                 THIS% logP1_2
        PRINT *, "  Base 10 exponent of first fiducial density (g/cm^3) = ", &
                 THIS% logRho0_2
        PRINT *, "  Base 10 exponent of second fiducial density (g/cm^3) = ",&
                 THIS% logRho1_2
        PRINT *, "  Base 10 exponent of third fiducial density (g/cm^3) = ", &
                 THIS% logRho2_2
        PRINT *

      ELSE

        PRINT *, "** ERROR in SUBROUTINE import_lorene_id_params!", &
                 " The equation of state is unknown!"
        STOP

      ENDIF

    ENDIF

  END PROCEDURE print_id_params


END SUBMODULE bns_params
