! File:         submodule_bns_lorene_properties.f90
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

SUBMODULE (bns_lorene) properties

  !********************************************
  !
  !# Implementation of the methods of TYPE [[bnslorene]]
  !  that read from |lorene| the
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


  MODULE PROCEDURE read_bns_properties

    !***************************************************
    !
    !# Store the parameters of the binary neutron
    !  stars' |lorene| ID into member variables
    !
    !  FT 5.10.2020
    !
    !***************************************************

    USE, INTRINSIC :: ISO_C_BINDING,  ONLY: C_CHAR

    USE tabulated_eos,  ONLY: read_compose_beta_equilibrated_eos
    USE constants,      ONLY: c_light, cm2km
    USE utility,        ONLY: Msun_geo, km2m, &
                              density_si2cu, k_lorene2cu, &
                              k_lorene2cu_pwp, &
                              zero, two, four, five, &
                              eos$poly, eos$pwpoly, eos$tabu$compose

#if flavour == 1

  USE sphincs_id_full,    ONLY: shorten_eos_name

#elif flavour == 2

  USE sphincs_id_lorene,  ONLY: shorten_eos_name

#endif

    IMPLICIT NONE

    INTEGER:: i, nchars
    INTEGER, PARAMETER:: str_length = 100
    INTEGER, PARAMETER:: str_length2= 500

    CHARACTER(KIND= C_CHAR), DIMENSION(str_length) :: eos1_tmp_c
    CHARACTER(KIND= C_CHAR), DIMENSION(str_length) :: eos2_tmp_c
    CHARACTER(KIND= C_CHAR), DIMENSION(str_length2):: eostable1_tmp_c
    CHARACTER(KIND= C_CHAR), DIMENSION(str_length2):: eostable2_tmp_c

    PRINT *, "** Executing the read_bns_properties subroutine..."

    CALL get_lorene_bns_params( this% bns_ptr,                 &
                                this% angular_vel,             &
                                this% distance,                &
                                this% distance_com,            &
                                this% mass1,                   &
                                this% mass2,                   &
                                this% mass_grav1,              &
                                this% mass_grav2,              &
                                this% adm_mass,                &
                                this% linear_momentum_x,       &
                                this% linear_momentum_y,       &
                                this% linear_momentum_z,       &
                                this% angular_momentum_x,      &
                                this% angular_momentum_y,      &
                                this% angular_momentum_z,      &
                                this% area_radius1,            &
                                this% radius1_x_comp,          &
                                this% radius1_y,               &
                                this% radius1_z,               &
                                this% radius1_x_opp,           &
                                this% center1_x,               &
                                this% barycenter1_x,           &
                                this% area_radius2,            &
                                this% radius2_x_comp,          &
                                this% radius2_y,               &
                                this% radius2_z,               &
                                this% radius2_x_opp,           &
                                this% center2_x,               &
                                this% barycenter2_x,           &
                                this% ent_center1,             &
                                this% nbar_center1,            &
                                this% rho_center1,             &
                                this% energy_density_center1,  &
                                this% specific_energy_center1, &
                                this% pressure_center1,        &
                                this% ent_center2,             &
                                this% nbar_center2,            &
                                this% rho_center2,             &
                                this% energy_density_center2,  &
                                this% specific_energy_center2, &
                                this% pressure_center2,        &
                                eos1_tmp_c,                    &
                                eos2_tmp_c,                    &
                                this% eos1_loreneid,           &
                                this% eos2_loreneid,           &
                                this% gamma_1,                 &
                                this% kappa_1,                 &
                                this% gamma_2,                 &
                                this% kappa_2,                 &
                                this% npeos_1,                 &
                                this% gamma0_1,                &
                                this% gamma1_1,                &
                                this% gamma2_1,                &
                                this% gamma3_1,                &
                                this% kappa0_1,                &
                                this% kappa1_1,                &
                                this% kappa2_1,                &
                                this% kappa3_1,                &
                                this% logP1_1,                 &
                                this% logRho0_1,               &
                                this% logRho1_1,               &
                                this% logRho2_1,               &
                                this% npeos_2,                 &
                                this% gamma0_2,                &
                                this% gamma1_2,                &
                                this% gamma2_2,                &
                                this% gamma3_2,                &
                                this% kappa0_2,                &
                                this% kappa1_2,                &
                                this% kappa2_2,                &
                                this% kappa3_2,                &
                                this% logP1_2,                 &
                                this% logRho0_2,               &
                                this% logRho1_2,               &
                                this% logRho2_2,               &
                                eostable1_tmp_c,               &
                                eostable2_tmp_c )


    ! Convert distances from |lorene| units (km) to SPHINCS units (Msun_geo)
    this% distance      = this% distance/Msun_geo
    this% distance_com  = this% distance_com/Msun_geo
    this% area_radius1  = this% area_radius1/Msun_geo
    this% radius1_x_comp= this% radius1_x_comp/Msun_geo
    this% radius1_y     = this% radius1_y/Msun_geo
    this% radius1_z     = this% radius1_z/Msun_geo
    this% radius1_x_opp = this% radius1_x_opp/Msun_geo
    this% center1_x     = this% center1_x/Msun_geo
    this% barycenter1_x = this% barycenter1_x/Msun_geo
    this% area_radius2  = this% area_radius2/Msun_geo
    this% radius2_x_comp= this% radius2_x_comp/Msun_geo
    this% radius2_y     = this% radius2_y/Msun_geo
    this% radius2_z     = this% radius2_z/Msun_geo
    this% radius2_x_opp = this% radius2_x_opp/Msun_geo
    this% center2_x     = this% center2_x/Msun_geo
    this% barycenter2_x = this% barycenter2_x/Msun_geo

    this% mass(1)= this% mass1
    this% mass(2)= this% mass2

    this% mass_grav(1)= this% mass_grav1
    this% mass_grav(2)= this% mass_grav2

    this% radii(1,:)= [this% radius1_x_opp, this% radius1_x_comp, &
                       this% radius1_y, this% radius1_y, &
                       this% radius1_z, this% radius1_z]
    this% radii(2,:)= [this% radius2_x_comp, this% radius2_x_opp, &
                       this% radius2_y, this% radius2_y, &
                       this% radius2_z, this% radius2_z]

    this% center(1,:)= [this% center1_x, zero, zero]
    this% center(2,:)= [this% center2_x, zero, zero]

    this% barycenter(1,:)= [this% barycenter1_x, zero, zero]
    this% barycenter(2,:)= [this% barycenter2_x, zero, zero]

    ! Convert hydro quantities from |lorene| units to SPHINCS units
    this% nbar_center1           = this% nbar_center1*(MSun_geo*km2m)**3
    this% rho_center1            = this% rho_center1*density_si2cu
    this% energy_density_center1 = this% energy_density_center1*density_si2cu
    this% pressure_center1       = this% pressure_center1*density_si2cu
    this% nbar_center2           = this% nbar_center2*(MSun_geo*km2m)**3
    this% rho_center2            = this% rho_center2*density_si2cu
    this% energy_density_center2 = this% energy_density_center2*density_si2cu
    this% pressure_center2       = this% pressure_center2*density_si2cu

    !
    !-- Convert C++ strings to Fortran strings
    !

    ! Name of EOS for star 1
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

    ALLOCATE( CHARACTER(nchars):: this% eos1, STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...allocation error for string eos1 in SUBROUTINE", &
                " read_bns_properties in SUBMODULE bns_lorene@properties.",&
                "The error message is ", err_msg
       PRINT *, "The STAT variable is ", ios
       PRINT *
       STOP
    ENDIF
    this% eos1= TRANSFER( eos1_tmp_c(1:nchars), this% eos1 )
    !DO i= 1, nchars, 1
    !  this% eos1(i:i)= eos1_tmp(i)
    !ENDDO

    ! Name of EOS for star 2
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

    ALLOCATE( CHARACTER(nchars):: this% eos2, STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...allocation error for string eos2 in SUBROUTINE", &
                " read_bns_properties in SUBMODULE bns_lorene@properties.", &
                "The error message is ", err_msg
       PRINT *, "The STAT variable is ", ios
       PRINT *
       STOP
    ENDIF
    this% eos2= TRANSFER( eos2_tmp_c(1:nchars), this% eos2 )
    !DO i= 1, nchars, 1
    !  this% eos2(i:i)= eos2_tmp(i)
    !ENDDO

    this% eos1= shorten_eos_name(this% eos1)
    this% eos2= shorten_eos_name(this% eos2)

    ! Name of file containing the EOS table for star 1
    i= 1
    DO
      IF( eostable1_tmp_c(i) == C_NULL_CHAR .OR. i == str_length ) EXIT
      i= i + 1
    ENDDO
    nchars = i - 1

    ALLOCATE( CHARACTER(nchars):: this% eos_table1, &
              STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...allocation error for string eos_tables1 in SUBROUTINE", &
                " read_bns_properties in SUBMODULE bns_lorene@properties.", &
                "The error message is ", err_msg
       PRINT *, "The STAT variable is ", ios
       PRINT *
       STOP
    ENDIF
    this% eos_table1= &
      TRANSFER( eostable1_tmp_c(1:nchars), this% eos_table1 )

    ! Name of file containing the EOS table for star 2
    i= 1
    DO
      IF( eostable2_tmp_c(i) == C_NULL_CHAR .OR. i == str_length ) EXIT
      i= i + 1
    ENDDO
    nchars = i - 1

    ALLOCATE( CHARACTER(nchars):: this% eos_table2, &
              STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...allocation error for string eos_tables2 in SUBROUTINE", &
                " read_bns_properties in SUBMODULE bns_lorene@properties.", &
                "The error message is ", err_msg
       PRINT *, "The STAT variable is ", ios
       PRINT *
       STOP
    ENDIF
    this% eos_table2= &
      TRANSFER( eostable2_tmp_c(1:nchars), this% eos_table2 )


    !
    !-- Convert polytropic constants from |lorene| units to SPHINCS units,
    !-- and assign |sphincsid| identifiers to the |eos|
    !

    ! Star 1
    IF( this% eos1_loreneid == 1 )THEN
    ! If the EOS is polytropic

      this% eos1_id= eos$poly
      this% kappa_1= this% kappa_1*k_lorene2cu( this% gamma_1 )

    ELSEIF( this% eos1_loreneid == 110 )THEN
    ! If the EOS is piecewise polytropic

      this% eos1_id = eos$pwpoly
      this% kappa0_1= this% kappa0_1 &
                      *k_lorene2cu_pwp( this% gamma0_1 )
      this% kappa1_1= this% kappa1_1 &
                      *k_lorene2cu_pwp( this% gamma1_1 )
      this% kappa2_1= this% kappa2_1 &
                      *k_lorene2cu_pwp( this% gamma2_1 )
      this% kappa3_1= this% kappa3_1 &
                      *k_lorene2cu_pwp( this% gamma3_1 )

    ELSEIF( this% eos1_loreneid == 17 .OR. this% eos1_loreneid == 20 )THEN
    ! If the EOS is tabulated, in CompOSE format

      IF(.NOT.ALLOCATED(this% tab_eos)) ALLOCATE(this% tab_eos(2))

      this% eos1_id= eos$tabu$compose

      CALL read_compose_beta_equilibrated_eos &
        (this% eos_table1, this% tab_eos(1)% table_eos)

    ELSE

      PRINT *, "** ERROR in SUBROUTINE read_bns_properties!", &
               " The equation of state on star 1 is unknown! LORENE EOS ID=", &
               this% eos1_loreneid
      PRINT *
      STOP

    ENDIF

    ! Star 2
    IF( this% eos2_loreneid == 1 )THEN
    ! If the EOS is polytropic

      this% eos2_id= eos$poly
      this% kappa_2= this% kappa_2*k_lorene2cu( this% gamma_2 )

    ELSEIF( this% eos2_loreneid == 110 )THEN
    ! If the EOS is piecewise polytropic

      this% eos2_id = eos$pwpoly
      this% kappa0_2= this% kappa0_2 &
                      *k_lorene2cu_pwp( this% gamma0_2 )
      this% kappa1_2= this% kappa1_2 &
                      *k_lorene2cu_pwp( this% gamma1_2 )
      this% kappa2_2= this% kappa2_2 &
                      *k_lorene2cu_pwp( this% gamma2_2 )
      this% kappa3_2= this% kappa3_2 &
                      *k_lorene2cu_pwp( this% gamma3_2 )

    ELSEIF( this% eos2_loreneid == 17 .OR. this% eos2_loreneid == 20 )THEN
    ! If the EOS is tabulated, in CompOSE format

      IF(.NOT.ALLOCATED(this% tab_eos)) ALLOCATE(this% tab_eos(2))

      this% eos2_id= eos$tabu$compose

      CALL read_compose_beta_equilibrated_eos &
        (this% eos_table2, this% tab_eos(2)% table_eos)

    ELSE

      PRINT *, "** ERROR in SUBROUTINE read_bns_properties!", &
               " The equation of state on star 2 is unknown! LORENE EOS ID=", &
               this% eos2_loreneid
      PRINT *
      STOP

    ENDIF

    ! Compute mOmega (see documentation for details)
    this% mOmega= this% angular_vel/(c_light*cm2km) &
                  *(this% mass_grav1 + this% mass_grav2)*Msun_geo

    ! Compute t_merger (see documentation for details)
    this% t_merger= five/(two**8)*( this% distance**4 ) &
                    /( this% mass_grav1*this% mass_grav2* &
                       ( this% mass_grav1 + this% mass_grav2 ) )


    CALL print_bns_properties(this)


  END PROCEDURE read_bns_properties


END SUBMODULE properties
