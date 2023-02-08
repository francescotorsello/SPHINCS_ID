! File:         submodule_diffstar_lorene_params.f90
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


  MODULE PROCEDURE read_diffstar_params

    !***************************************************
    !
    !# Store the parameters of the binary neutron
    !  stars' |lorene| ID into member variables
    !
    !  FT 5.10.2020
    !
    !***************************************************

    USE, INTRINSIC :: ISO_C_BINDING,  ONLY: C_CHAR, C_NULL_CHAR
    USE utility, ONLY: Msun_geo, km2m, lorene2hydrobase, k_lorene2hydrobase, &
                       k_lorene2hydrobase_piecewisepolytrope, zero, two, &
                       eos$poly, eos$pwpoly

#if flavour == 1

  USE sphincs_id_full,    ONLY: shorten_eos_name

#elif flavour == 2

  USE sphincs_id_lorene,  ONLY: shorten_eos_name

#endif

    IMPLICIT NONE

    INTEGER:: i, nchars
    INTEGER, PARAMETER:: str_length= 100

    CHARACTER(KIND= C_CHAR), DIMENSION(str_length):: eos_tmp_c

    PRINT *
    PRINT *, "** Executing the read_diffstar_params subroutine..."

    CALL get_diffstar_params( this% diffstar_ptr,                   &
                              this% omega_c,                        &
                              this% mass,                           &
                              this% mass_grav,                      &
                              this% angular_momentum,               &
                              this% tsw,                            &
                              this% grv2,                           &
                              this% grv3,                           &
                              this% r_circ,                         &
                              this% surface_area,                   &
                              this% r_mean,                         &
                              this% r_eq,                           &
                              this% r_eq_pi2,                       &
                              this% r_eq_pi,                        &
                              this% r_eq_3pi2,                      &
                              this% r_pole,                         &
                              this% r_ratio,                        &
                              this% r_isco,                         &
                              this% f_isco,                         &
                              this% specific_energy_isco,           &
                              this% specific_angular_momentum_isco, &
                              this% area_radius,                    &
                              this% ent_center,                     &
                              this% nbar_center,                    &
                              this% rho_center,                     &
                              this% energy_density_center,          &
                              this% specific_energy_center,         &
                              this% pressure_center,                &
                              this% redshift_eqf,                   &
                              this% redshift_eqb,                   &
                              this% redshift_pole,                  &
                              eos_tmp_c,                            &
                              this% eos_loreneid,                   &
                              this% gamma,                          &
                              this% kappa,                          &
                              this% npeos,                          &
                              this% gamma0,                         &
                              this% gamma1,                         &
                              this% gamma2,                         &
                              this% gamma3,                         &
                              this% kappa0,                         &
                              this% kappa1,                         &
                              this% kappa2,                         &
                              this% kappa3,                         &
                              this% logP1,                          &
                              this% logRho0,                        &
                              this% logRho1,                        &
                              this% logRho2 )

    ! Convert distances from |lorene| units (km) to SPHINCS units (Msun_geo)
    ! See MODULE constants for the definition of Msun_geo
    this% r_circ      = this% r_circ/Msun_geo
    this% r_mean      = this% r_mean/Msun_geo
    this% area_radius = this% area_radius/Msun_geo
    this% r_eq        = this% r_eq/Msun_geo
    this% r_eq_pi2    = this% r_eq_pi2/Msun_geo
    this% r_eq_pi     = this% r_eq_pi/Msun_geo
    this% r_eq_3pi2   = this% r_eq_3pi2/Msun_geo
    this% r_pole      = this% r_pole/Msun_geo
    this% r_isco      = this% r_isco/Msun_geo
    this% surface_area= this% surface_area/(Msun_geo**two)

    ! Note that here the radii of the star along the z axis are set equal to
    ! the equatorial ones. This is because the star can have a toroidal shape
    ! and a box with a z-size determined by the polar radius won't enclose
    ! the entire star. This would be problematic both in the particle and in the
    ! tpo objects.
    this% radii(:)= [this% r_eq_pi, this% r_eq, &
                     this% r_eq_pi2, this% r_eq_pi2, &
                     this% r_eq_pi2, this% r_eq_pi2]
                     !this% r_pole, this% r_pole]

    this% center(:)= [zero, zero, zero]

    this% barycenter(:)= [zero, zero, zero]

    ! Convert hydro quantities from |lorene| units to SPHINCS units
    this% nbar_center           = this% nbar_center*(MSun_geo*km2m)**3
    this% rho_center            = this% rho_center*lorene2hydrobase
    this% energy_density_center = this% energy_density_center*lorene2hydrobase
    this% pressure_center       = this% pressure_center*lorene2hydrobase

    ! Convert polytropic constants from |lorene| units to SPHINCS units
    IF( this% eos_loreneid == 1 )THEN
    ! If the EOS is polytropic

      this% eos_id= eos$poly
      this% kappa= this% kappa*k_lorene2hydrobase( this% gamma )

    ELSEIF( this% gamma0 /= 0 )THEN
    ! If the EOS is piecewise polytropic

      this% eos_id = eos$pwpoly
      this% kappa0= this% kappa0 &
                      *k_lorene2hydrobase_piecewisepolytrope( this% gamma0 )
      this% kappa1= this% kappa1 &
                      *k_lorene2hydrobase_piecewisepolytrope( this% gamma1 )
      this% kappa2= this% kappa2 &
                      *k_lorene2hydrobase_piecewisepolytrope( this% gamma2 )
      this% kappa3= this% kappa3 &
                      *k_lorene2hydrobase_piecewisepolytrope( this% gamma3 )

    !ELSEIF( this% eos_loreneid == 17 .OR. this% eos_loreneid == 20 )THEN
    ! If the EOS is tabulated

    ELSE

      PRINT *, "** ERROR in SUBROUTINE read_lorene_id_params!", &
               " The equation of state is unknown! LORENE EOS IDs=", &
               this% eos_loreneid, ", ", this% eos_loreneid
      STOP

    ENDIF

    ! Convert C++ strings to FORTRAN strings
    i= 1
    DO
      IF( eos_tmp_c(i) == C_NULL_CHAR .OR. i == str_length ) EXIT
      i= i + 1
    ENDDO
    nchars = i - 1

    ALLOCATE( CHARACTER(nchars):: this% eos, STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...allocation error for string eos. ", &
                "The error message is ", err_msg
       PRINT *, "The STAT variable is ", ios
       STOP
    ENDIF
    this% eos= TRANSFER( eos_tmp_c(1:nchars), this% eos )

    this% eos= shorten_eos_name(this% eos)

    CALL print_diffstar_params( this )

    PRINT *, "** Subroutine read_diffstar_params executed."
    PRINT *

  END PROCEDURE read_diffstar_params


END SUBMODULE params
