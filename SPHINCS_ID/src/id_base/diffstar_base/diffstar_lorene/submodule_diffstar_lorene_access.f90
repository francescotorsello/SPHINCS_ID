! File:         submodule_diffstar_lorene_access.f90
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

SUBMODULE (diffstar_lorene) access

  !***************************************************
  !
  !# The module contains the implementation of the
  !  methods of TYPE diffstar that allow to access PRIVATE
  !  members.
  !
  !  FT 25.10.2021
  !
  !***************************************************


  IMPLICIT NONE


  CONTAINS


  !-----------------!
  !--  FUNCTIONS  --!
  !-----------------!


  MODULE PROCEDURE get_field_array

    !***********************************************
    !
    !# Returns one of the member arrays, selected
    !  with the string input.
    !
    !  FT 25.10.2021
    !
    !***********************************************

    IMPLICIT NONE

    select_field: SELECT CASE( field )

    CASE( "lapse" )

      field_array= this% lapse

    CASE( "shift_x" )

      field_array= this% shift_x

    CASE( "shift_y" )

      field_array= this% shift_y

    CASE( "shift_z" )

      field_array= this% shift_z

    CASE( "g_xx" )

      field_array= this% g_xx

    CASE( "g_xy" )

      field_array= this% g_xy

    CASE( "g_xz" )

      field_array= this% g_xz

    CASE( "g_yy" )

      field_array= this% g_yy

    CASE( "g_yz" )

      field_array= this% g_yz

    CASE( "g_zz" )

      field_array= this% g_zz

    CASE( "k_xx" )

      field_array= this% k_xx

    CASE( "k_xy" )

      field_array= this% k_xy

    CASE( "k_xz" )

      field_array= this% k_xz

    CASE( "k_yy" )

      field_array= this% k_yy

    CASE( "k_yz" )

      field_array= this% k_yz

    CASE( "k_zz" )

      field_array= this% k_zz

    CASE( "baryon_density" )

      field_array= this% baryon_density

    CASE( "energy_density" )

      field_array= this% energy_density

    CASE( "specific_energy" )

      field_array= this% specific_energy

    CASE( "v_euler_x" )

      field_array= this% v_euler_x

    CASE( "v_euler_y" )

      field_array= this% v_euler_y

    CASE( "v_euler_z" )

      field_array= this% v_euler_z

    CASE DEFAULT

      PRINT *, "** There is no field named ", field, "in TYPE diffstar."
      STOP

    END SELECT select_field

  END PROCEDURE get_field_array


  MODULE PROCEDURE get_field_value

    !************************************************
    !
    !# Returns the value of one of the member arrays,
    !  selected with the string input, at the point
    !  given as argument.
    !
    !  FT 25.10.2021
    !
    !************************************************

    IMPLICIT NONE

    select_field: SELECT CASE( field )

    CASE( "lapse" )

      field_value= this% lapse( n )

    CASE( "shift_x" )

      field_value= this% shift_x( n )

    CASE( "shift_y" )

      field_value= this% shift_y( n )

    CASE( "shift_z" )

      field_value= this% shift_z( n )

    CASE( "g_xx" )

      field_value= this% g_xx( n )

    CASE( "g_xy" )

      field_value= this% g_xy( n )

    CASE( "g_xz" )

      field_value= this% g_xz( n )

    CASE( "g_yy" )

      field_value= this% g_yy( n )

    CASE( "g_yz" )

      field_value= this% g_yz( n )

    CASE( "g_zz" )

      field_value= this% g_zz( n )

    CASE( "k_xx" )

      field_value= this% k_xx( n )

    CASE( "k_xy" )

      field_value= this% k_xy( n )

    CASE( "k_xz" )

      field_value= this% k_xz( n )

    CASE( "k_yy" )

      field_value= this% k_yy( n )

    CASE( "k_yz" )

      field_value= this% k_yz( n )

    CASE( "k_zz" )

      field_value= this% k_zz( n )

    CASE( "baryon_density" )

      field_value= this% baryon_density( n )

    CASE( "energy_density" )

      field_value= this% energy_density( n )

    CASE( "specific_energy" )

      field_value= this% specific_energy( n )

    CASE( "v_euler_x" )

      field_value= this% v_euler_x( n )

    CASE( "v_euler_y" )

      field_value= this% v_euler_y( n )

    CASE( "v_euler_z" )

      field_value= this% v_euler_z( n )

    CASE DEFAULT

      PRINT *, "** There is no field named ", field, "in TYPE diffstar."
      STOP

    END SELECT select_field

  END PROCEDURE get_field_value


  MODULE PROCEDURE get_diffstar_identifier

    !************************************************
    !
    !# Returns the value of [[diffstar_identifier]], the
    !  integer identifier of the diffstar object
    !
    !  FT 25.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_diffstar_identifier= this% diffstar_identifier

  END PROCEDURE get_diffstar_identifier


  !MODULE PROCEDURE get_diffstar_ptr
  !
  !  !************************************************
  !  !
  !  !# Returns the value of [[diffstar_ptr]], the C pointer
  !  ! to the |lorene|'s |etdiffrot| object
  !  ! N.B. This variable is global. The pointer
  !  !      to the second |lorene| |etdiffrot| object will
  !  !      overwrite the first one, and so on.
  !  !      This variable stores the pointer to
  !  !      the last defined |lorene| |etdiffrot| object.
  !  !      That's why it is not freed in the
  !  !      destructor of a diffstar object. Presently, it
  !  !      has to be freed by the user at the end of
  !  !      the PROGRAM. See the last part of the
  !  !      PROGRAM in setup_lorene_id.f90, for
  !  !      example.
  !  !
  !  !  FT 25.10.2021
  !  !
  !  !************************************************
  !
  !  IMPLICIT NONE
  !
  !  get_diffstar_ptr= this% diffstar_ptr
  !
  !END PROCEDURE get_diffstar_ptr


  MODULE PROCEDURE get_eos_loreneid

    !**************************************************
    !
    !# Returns the |lorene| ID-number of the EOS of the DRS
    !
    !  FT
    !
    !**************************************************

    IMPLICIT NONE

    get_eos_loreneid= this% eos_loreneid

  END PROCEDURE get_eos_loreneid


  MODULE PROCEDURE get_eos_parameters

    !**************************************************
    !
    !# Returns the |eos| parameters of the DRS
    !
    !  FT 2.11.2021
    !
    !**************************************************

    USE utility,  ONLY: eos$poly, eos$pwpoly

    IMPLICIT NONE

    CALL this% check_i_matter(i_matter)

    IF( this% eos_id == eos$poly )THEN

      eos_params= [ DBLE(this% eos_id), this% gamma, this% kappa ]

    ELSEIF( this% eos_id == eos$pwpoly )THEN

      eos_params= [ DBLE(this% eos_id), DBLE(this% npeos), &
            this% gamma0, this% gamma1, this% gamma2, this% gamma3, &
            this% kappa0, this% kappa1, this% kappa2, this% kappa3, &
            this% logP1, &
            this% logRho0, this% logRho1, this% logRho2 ]

    !ELSEIF( this% eos_loreneid == 17 .OR. this% eos_loreneid == 20 )THEN
    !
    !  eos_params= [ DBLE(this% eos_loreneid) ]

    ELSE

      PRINT *, "** ERROR in SUBROUTINE get_eos_parameters!", &
               " The EOS on the DRS is unknown! SPHINCS_ID EOS ID=", &
               this% eos_id
      STOP

    ENDIF

  END PROCEDURE get_eos_parameters


END SUBMODULE access
