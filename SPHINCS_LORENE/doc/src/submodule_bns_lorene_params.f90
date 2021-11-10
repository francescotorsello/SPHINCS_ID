! File:         submodule_bns_lorene_params.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (bns_lorene) bns_lorene_params

  !********************************************
  !
  !# Implementation of the methods of TYPE bns
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


  MODULE PROCEDURE import_id_params

    !***************************************************
    !
    !# Store the parameters of the binary neutron
    !  stars' |lorene| ID into member variables
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
                               THIS% logRho2_2 )

    ! Convert distances from |lorene| units (km) to SPHINCS units (Msun_geo)
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

    THIS% mass(1)= THIS% mass1
    THIS% mass(2)= THIS% mass2

    THIS% radii(1,:)= [THIS% radius1_x_opp, THIS% radius1_x_comp, &
                       THIS% radius1_y, THIS% radius1_y, &
                       THIS% radius1_z, THIS% radius1_z]
    THIS% radii(2,:)= [THIS% radius2_x_opp, THIS% radius2_x_comp, &
                       THIS% radius2_y, THIS% radius2_y, &
                       THIS% radius2_z, THIS% radius2_z]

    THIS% center(1,:)= [THIS% center1_x, 0.0D0, 0.0D0]
    THIS% center(2,:)= [THIS% center2_x, 0.0D0, 0.0D0]

    THIS% barycenter(1,:)= [THIS% barycenter1_x, 0.0D0, 0.0D0]
    THIS% barycenter(2,:)= [THIS% barycenter2_x, 0.0D0, 0.0D0]

    ! Convert hydro quantities from |lorene| units to SPHINCS units
    THIS% nbar_center1           = THIS% nbar_center1*(MSun_geo*km2m)**3
    THIS% rho_center1            = THIS% rho_center1*lorene2hydrobase
    THIS% energy_density_center1 = THIS% energy_density_center1*lorene2hydrobase
    THIS% pressure_center1       = THIS% pressure_center1*lorene2hydrobase
    THIS% nbar_center2           = THIS% nbar_center2*(MSun_geo*km2m)**3
    THIS% rho_center2            = THIS% rho_center2*lorene2hydrobase
    THIS% energy_density_center2 = THIS% energy_density_center2*lorene2hydrobase
    THIS% pressure_center2       = THIS% pressure_center2*lorene2hydrobase

    ! Convert polytropic constants from |lorene| units to SPHINCS units
    IF( THIS% eos1_loreneid == 1 )THEN ! If the EOS is polytropic

      THIS% kappa_1= THIS% kappa_1*k_lorene2hydrobase( THIS% gamma_1 )
      THIS% kappa_2= THIS% kappa_2*k_lorene2hydrobase( THIS% gamma_2 )

    ELSEIF( THIS% gamma0_1 /= 0 )THEN ! If the EOS is piecewise polytropic

      THIS% kappa0_1= THIS% kappa0_1 &
                      *k_lorene2hydrobase_piecewisepolytrope( THIS% gamma0_1 )
      THIS% kappa1_1= THIS% kappa1_1 &
                      *k_lorene2hydrobase_piecewisepolytrope( THIS% gamma1_1 )
      THIS% kappa2_1= THIS% kappa2_1 &
                      *k_lorene2hydrobase_piecewisepolytrope( THIS% gamma2_1 )
      THIS% kappa3_1= THIS% kappa3_1 &
                      *k_lorene2hydrobase_piecewisepolytrope( THIS% gamma3_1 )
      THIS% kappa0_2= THIS% kappa0_2 &
                      *k_lorene2hydrobase_piecewisepolytrope( THIS% gamma0_2 )
      THIS% kappa1_2= THIS% kappa1_2 &
                      *k_lorene2hydrobase_piecewisepolytrope( THIS% gamma1_2 )
      THIS% kappa2_2= THIS% kappa2_2 &
                      *k_lorene2hydrobase_piecewisepolytrope( THIS% gamma2_2 )
      THIS% kappa3_2= THIS% kappa3_2 &
                      *k_lorene2hydrobase_piecewisepolytrope( THIS% gamma3_2 )

    ELSEIF( THIS% eos1_loreneid == 17 .OR. THIS% eos1_loreneid == 20 )THEN
    ! If the EOS is tabulated

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


END SUBMODULE bns_lorene_params
