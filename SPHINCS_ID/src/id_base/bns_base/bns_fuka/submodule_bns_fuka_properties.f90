! File:         submodule_bns_fuka_properties.f90
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

SUBMODULE (bns_fuka) properties

  !********************************************
  !
  !# Implementation of the methods of TYPE [[bnsfuka]]
  !  that import from |fuka| the
  !  parameters of the binary system,
  !  and print them to the standard output.
  !
  !  FT 09.02.2022
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
    !  stars' |fuka| |id| into member variables
    !
    !  FT 09.02.2022
    !
    !***************************************************

    USE, INTRINSIC :: ISO_C_BINDING,  ONLY: C_CHAR

    USE constants,  ONLY: c_light, cm2km
    USE utility,    ONLY: Msun_geo, zero, one, two, four, five, &
                          eos$poly, eos$pwpoly

#if flavour == 1

  USE sphincs_id_full,  ONLY: shorten_eos_name_fuka

#elif flavour == 3

  USE sphincs_id_fuka,  ONLY: shorten_eos_name_fuka

#endif

    IMPLICIT NONE

    INTEGER, PARAMETER:: str_length= 100
    LOGICAL, PARAMETER:: debug= .FALSE.

    INTEGER:: i, nchars

    CHARACTER(KIND= C_CHAR), DIMENSION(str_length):: eos_type_1_tmp_c
    CHARACTER(KIND= C_CHAR), DIMENSION(str_length):: eos_file_1_tmp_c
    CHARACTER(KIND= C_CHAR), DIMENSION(str_length):: eos_type_2_tmp_c
    CHARACTER(KIND= C_CHAR), DIMENSION(str_length):: eos_file_2_tmp_c


    CALL get_fuka_id_params( this% bns_ptr               , &
                             this% angular_vel           , &
                             this% distance              , &
                             this% mass(1)               , &
                             this% mass(2)               , &
                             this% mass_grav(1)          , &
                             this% mass_grav(2)          , &
                             this% radius1_x_opp         , &
                             this% radius1_x_comp        , &
                             this% radius2_x_opp         , &
                             this% radius2_x_comp        , &
                             this% adm_mass              , &
                             this% komar_mass            , &
                             this% linear_momentum_x     , &
                             this% linear_momentum_y     , &
                             this% linear_momentum_z     , &
                             this% angular_momentum_z    , &
                             this% barycenter_system(1)  , &
                             this% barycenter_system(2)  , &
                             this% barycenter_system(3)  , &
                             this% area_radius1          , &
                             this% center1_x             , &
                             this% area_radius2          , &
                             this% center2_x             , &
                             this% ent_center1           , &
                             this% rho_center1           , &
                             this% energy_density_center1, &
                             this% ent_center2           , &
                             this% rho_center2           , &
                             this% energy_density_center2, &
                             eos_type_1_tmp_c            , &
                             eos_file_1_tmp_c            , &
                             eos_type_2_tmp_c            , &
                             eos_file_2_tmp_c            , &
                             this% gamma_1               , &
                             this% kappa_1               , &
                             this% npeos_1               , &
                             this% gamma0_1              , &
                             this% gamma1_1              , &
                             this% gamma2_1              , &
                             this% gamma3_1              , &
                             this% kappa0_1              , &
                             this% kappa1_1              , &
                             this% kappa2_1              , &
                             this% kappa3_1              , &
                             this% logP1_1               , &
                             this% logRho0_1             , &
                             this% logRho1_1             , &
                             this% logRho2_1 )

    this% gamma_2  = this% gamma_1
    this% kappa_2  = this% kappa_1
    this% npeos_2  = this% npeos_1
    this% gamma0_2 = this% gamma0_1
    this% gamma1_2 = this% gamma1_1
    this% gamma2_2 = this% gamma2_1
    this% gamma3_2 = this% gamma3_1
    this% kappa0_2 = this% kappa0_1
    this% kappa1_2 = this% kappa1_1
    this% kappa2_2 = this% kappa2_1
    this% kappa3_2 = this% kappa3_1
    this% logP1_2  = this% logP1_1
    this% logRho0_2= this% logRho0_1
    this% logRho1_2= this% logRho1_1
    this% logRho2_2= this% logRho2_1

    PRINT *, "** Finding centers and radii of the stars..."
    PRINT *
    IF(debug) PRINT *, "max radius 1=", this% radius1_x_comp
    IF(debug) PRINT *, "min radius 1=", this% radius1_x_opp
    IF(debug) PRINT *, "max radius 2=", this% radius2_x_comp
    IF(debug) PRINT *, "min radius 2=", this% radius2_x_opp
    IF(debug) PRINT *

    ! Find the centers of the stars
    this% center1_x= this% find_center(this% distance, -one, read_density)
    this% center2_x= this% find_center(this% distance, one, read_density)

    ! Find the radii of the stars
    this% radius1_x_comp= this% find_radius([this% center1_x,zero,zero], &
                                       [one,zero,zero], read_density)
    this% radius1_x_opp= this% find_radius([this% center1_x,zero,zero], &
                                       [-one,zero,zero], read_density)
    this% radius1_y= this% find_radius([this% center1_x,zero,zero], &
                                       [zero,one,zero], read_density)
    this% radius1_z= this% find_radius([this% center1_x,zero,zero], &
                                       [zero,zero,one], read_density)

    this% radius2_x_comp= this% find_radius([this% center2_x,zero,zero], &
                                       [-one,zero,zero], read_density)
    this% radius2_x_opp= this% find_radius([this% center2_x,zero,zero], &
                                       [one,zero,zero], read_density)
    this% radius2_y= this% find_radius([this% center2_x,zero,zero], &
                                       [zero,one,zero], read_density)
    this% radius2_z= this% find_radius([this% center2_x,zero,zero], &
                                       [zero,zero,one], read_density)

    IF(debug) PRINT *
    IF(debug) PRINT *, "x radius 1+=", this% radius1_x_comp
    IF(debug) PRINT *, "x radius 1-=", this% radius1_x_opp
    IF(debug) PRINT *, "y radius 1=", this% radius1_y
    IF(debug) PRINT *, "z radius 1=", this% radius1_z
    IF(debug) PRINT *, "x radius 2-=", this% radius2_x_comp
    IF(debug) PRINT *, "x radius 2+=", this% radius2_x_opp
    IF(debug) PRINT *, "y radius 2=", this% radius2_y
    IF(debug) PRINT *, "z radius 2=", this% radius2_z
    IF(debug) PRINT *
    IF(debug) PRINT *, "this% center1_x=", this% center1_x
    IF(debug) PRINT *, "this% center2_x=", this% center2_x
    IF(debug) PRINT *, "this% distance=", this% distance
    IF(debug) PRINT *
    IF(debug) PRINT *, "this% rho_center1=", this% rho_center1
    IF(debug) PRINT *, "this% rho_center2=", this% rho_center2
    IF(debug) PRINT *
    IF(debug) STOP

    this% radii(1,:)= [this% radius1_x_opp, this% radius1_x_comp, &
                       this% radius1_y, this% radius1_y, &
                       this% radius1_z, this% radius1_z]
    this% radii(2,:)= [this% radius2_x_comp, this% radius2_x_opp, &
                       this% radius2_y, this% radius2_y, &
                       this% radius2_z, this% radius2_z]

    this% center(1,:)= [this% center1_x, zero, zero]
    this% center(2,:)= [this% center2_x, zero, zero]

  !  this% barycenter(1,:)= [this% barycenter1_x, zero, zero]
  !  this% barycenter(2,:)= [this% barycenter2_x, zero, zero]

    ! Compute mOmega (see documentation for details)
    this% mOmega= this% angular_vel/(c_light*cm2km) &
                  *(this% mass_grav1 + this% mass_grav2)*Msun_geo

    ! Compute t_merger (see documentation for details)
    this% t_merger= five/(two**(two*four))*(this% distance**four) &
                    /( this% mass_grav1*this% mass_grav2* &
                       ( this% mass_grav1 + this% mass_grav2 ) )

    !
    !-- Convert C++ strings to FORTRAN strings
    !

    ! Star 1
    i= 1
    DO
      IF( eos_type_1_tmp_c(i) == C_NULL_CHAR .OR. i == str_length ) EXIT
      i= i + 1
    ENDDO
    nchars = i - 1

    ALLOCATE( CHARACTER(nchars):: this% eos_type_1, STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
      PRINT *, "...allocation error for string eos_type_1. ", &
              "The error message is ", err_msg
      PRINT *, "The STAT variable is ", ios
      PRINT *, "nchars= ", nchars
      STOP
    ENDIF
    this% eos_type_1= TRANSFER( eos_type_1_tmp_c(1:nchars), this% eos_type_1 )

    i= 1
    DO
      IF( eos_file_1_tmp_c(i) == C_NULL_CHAR .OR. i == str_length ) EXIT
      i= i + 1
    ENDDO
    nchars = i - 1

    ALLOCATE( CHARACTER(nchars):: this% eos_file_1, STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
      PRINT *, "...allocation error for string eos_file_1. ", &
              "The error message is ", err_msg
      PRINT *, "The STAT variable is ", ios
      PRINT *, "nchars= ", nchars
      STOP
    ENDIF
    this% eos_file_1= TRANSFER( eos_file_1_tmp_c(1:nchars), this% eos_file_1 )
    IF(debug) PRINT *, "this% eos_file_1=", this% eos_file_1
    this% eos1= shorten_eos_name_fuka(this% eos_file_1)

    ! Star 2
    i= 1
    DO
      IF( eos_type_2_tmp_c(i) == C_NULL_CHAR .OR. i == str_length ) EXIT
      i= i + 1
    ENDDO
    nchars = i - 1

    ALLOCATE( CHARACTER(nchars):: this% eos_type_2, STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
      PRINT *, "...allocation error for string eos_type_2. ", &
              "The error message is ", err_msg
      PRINT *, "The STAT variable is ", ios
      PRINT *, "nchars= ", nchars
      STOP
    ENDIF
    this% eos_type_2= TRANSFER( eos_type_2_tmp_c(1:nchars), this% eos_type_2 )
    this% eos_type_2= this% eos_type_1
    ! this% eos_type_1 is used because eos_type_2_tmp_c is empty...check if
    ! FUKA actually prints it

    i= 1
    DO
      IF( eos_file_2_tmp_c(i) == C_NULL_CHAR .OR. i == str_length ) EXIT
      i= i + 1
    ENDDO
    nchars = i - 1

    ALLOCATE( CHARACTER(nchars):: this% eos_file_2, STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
      PRINT *, "...allocation error for string eos_file_2. ", &
               "The error message is ", err_msg
      PRINT *, "The STAT variable is ", ios
      PRINT *, "nchars= ", nchars
      STOP
    ENDIF
    this% eos_file_2= TRANSFER( eos_file_2_tmp_c(1:nchars), this% eos_file_2 )
    this% eos_file_2= this% eos_file_1
    ! this% eos_file_1 is used because eos_file_2_tmp_c is empty...check if
    ! FUKA actually prints it
    IF(debug) PRINT *, "this% eos_file_2=", this% eos_file_2
    this% eos2= shorten_eos_name_fuka(this% eos_file_2)

    !
    !-- Assign |sphincsid| identifiers to the |eos|
    !

    ! Star 1
    IF( this% eos_type_1 == "Cold_PWPoly" )THEN

      IF( this% npeos_1 == 1 )THEN

        this% eos1_id= eos$poly

      ELSEIF( this% npeos_1 > 1 )THEN

        this% eos1_id = eos$pwpoly

      ENDIF

    ELSE

      PRINT *, "** ERROR in SUBROUTINE read_bns_properties!", &
               " The EOS on star 1 is unknown! FUKA EOS type=", &
               this% eos_type_1
      STOP

    ENDIF

    ! Star 2
    IF( this% eos_type_2 == "Cold_PWPoly" )THEN

      IF( this% npeos_2 == 1 )THEN

        this% eos2_id= eos$poly

      ELSEIF( this% npeos_2 > 1 )THEN

        this% eos2_id = eos$pwpoly

      ENDIF

    ELSE

      PRINT *, "** ERROR in SUBROUTINE read_bns_properties!", &
               " The EOS on star 2 is unknown! FUKA EOS type=", &
               this% eos_type_2
      STOP

    ENDIF


    CALL print_bns_properties(this)


    CONTAINS


    FUNCTION read_density(x, y, z) RESULT(rho)

      DOUBLE PRECISION, INTENT(IN):: x
      DOUBLE PRECISION, INTENT(IN):: y
      DOUBLE PRECISION, INTENT(IN):: z
      DOUBLE PRECISION:: rho

      rho= this% read_fuka_mass_density(x, y, z)

    END FUNCTION read_density


  END PROCEDURE read_bns_properties


END SUBMODULE properties
