! File:         submodule_bnsfuka_access.f90
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

SUBMODULE (bns_fuka) access

  !***************************************************
  !
  !# The module contains the implementation of the
  !  methods of TYPE [[bnsfuka]] that allow to access
  !  PRIVATE members.
  !
  !  FT 09.02.2022
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
    !  Created:     FT 09.02.2022
    !  Last update: FT 27.05.2022
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

    !CASE( "baryon_density" )
    !
    !  !field_array= this% mass_density
    !
    !CASE( "pressure" )
    !
    !  field_array= this% pressure
    !
    !CASE( "specific_energy" )
    !
    !  field_array= this% specific_energy
    !
    !CASE( "v_euler_x" )
    !
    !  field_array= this% v_euler_x
    !
    !CASE( "v_euler_y" )
    !
    !  field_array= this% v_euler_y
    !
    !CASE( "v_euler_z" )
    !
    !  field_array= this% v_euler_z

    CASE DEFAULT

      PRINT *, "** There is no field named ", field, "in TYPE bnsfuka."
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
    !  Created:     FT 09.02.2022
    !  Last update: FT 27.05.2022
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

    !CASE( "baryon_density" )
    !
    !  !field_value= this% mass_density( n )
    !
    !CASE( "pressure" )
    !
    !  field_value= this% pressure( n )
    !
    !CASE( "specific_energy" )
    !
    !  field_value= this% specific_energy( n )
    !
    !CASE( "v_euler_x" )
    !
    !  field_value= this% v_euler_x( n )
    !
    !CASE( "v_euler_y" )
    !
    !  field_value= this% v_euler_y( n )
    !
    !CASE( "v_euler_z" )
    !
    !  field_value= this% v_euler_z( n )

    CASE DEFAULT

      PRINT *, "** There is no field named ", field, "in TYPE bnsfuka."
      STOP

    END SELECT select_field

  END PROCEDURE get_field_value


  MODULE PROCEDURE get_bns_identifier

    !************************************************
    !
    !# Returns the value of [[bns_identifier]], the
    !  integer identifier of the bns object
    !
    !  Created:     FT 09.02.2022
    !  Last update: FT 27.05.2022
    !
    !************************************************

    IMPLICIT NONE

    get_bns_identifier= this% bns_identifier

  END PROCEDURE get_bns_identifier


  !MODULE PROCEDURE get_bns_ptr
  !
  !  !************************************************
  !  !
  !  !# Returns the value of [[bns_ptr]], the C pointer
  !  ! to the |fuka|'s Bin_NS object
  !  ! N.B. This variable is global. The pointer
  !  !      to the second |fuka| Bin_NS object will
  !  !      overwrite the first one, and so on.
  !  !      This variable stores the pointer to
  !  !      the last defined |fuka| Bin_NS object.
  !  !      That's why it is not freed in the
  !  !      destructor of a bns object. Presently, it
  !  !      has to be freed by the user at the end of
  !  !      the PROGRAM. See the last part of the
  !  !      PROGRAM in sphincs_id.f90, for
  !  !      example.
  !  !
  !  !  FT
  !  !
  !  !************************************************
  !
  !  IMPLICIT NONE
  !
  !  get_bns_ptr= this% bns_ptr
  !
  !END PROCEDURE get_bns_ptr


  MODULE PROCEDURE get_eos1_fukaid

    !**************************************************
    !
    !# Returns the |fuka| ID-number of the EOS for NS 1
    !
    !  FT
    !
    !**************************************************

    IMPLICIT NONE

  !  get_eos1_fukaid= this% eos1_fukaid

  END PROCEDURE get_eos1_fukaid


  MODULE PROCEDURE get_eos2_fukaid

    !**************************************************
    !
    !# Returns the |fuka| ID-number of the EOS for NS 2
    !
    !  FT
    !
    !**************************************************

    IMPLICIT NONE

  !  get_eos2_fukaid= this% eos2_fukaid

  END PROCEDURE get_eos2_fukaid


  MODULE PROCEDURE get_eos_parameters

    !**************************************************
    !
    !# Returns the |eos| parameters of the
    !  `i_matter`-s star
    !
    !  FT
    !
    !**************************************************

    USE utility,  ONLY: eos$poly, eos$pwpoly

    IMPLICIT NONE

    IF( i_matter == 1 )THEN

      IF( this% eos1_id == eos$poly )THEN

        eos_params= [ DBLE(eos$poly), this% gamma_1, this% kappa_1 ]

      ELSEIF( this% eos1_id == eos$pwpoly )THEN

        eos_params= [ DBLE(eos$pwpoly), DBLE(this% npeos_1), &
              this% gamma0_1, this% gamma1_1, this% gamma2_1, this% gamma3_1, &
              this% kappa0_1, this% kappa1_1, this% kappa2_1, this% kappa3_1, &
              this% logP1_1, &
              this% logRho0_1, this% logRho1_1, this% logRho2_1 ]

      ELSE

        PRINT *, "** ERROR in SUBROUTINE get_eos_parameters!", &
                 " The EOS on star 1 is unknown! FUKA EOS type=", &
                 this% eos_type_1
        STOP

      ENDIF

    ELSEIF( i_matter == 2 )THEN

      IF( this% eos1_id == eos$poly )THEN

        eos_params= [ DBLE(eos$poly), this% gamma_1, this% kappa_1 ]

      ELSEIF( this% eos1_id == eos$pwpoly )THEN

        eos_params= [ DBLE(eos$pwpoly), DBLE(this% npeos_1), &
              this% gamma0_1, this% gamma1_1, this% gamma2_1, this% gamma3_1, &
              this% kappa0_1, this% kappa1_1, this% kappa2_1, this% kappa3_1, &
              this% logP1_1, &
              this% logRho0_1, this% logRho1_1, this% logRho2_1 ]

      ELSE

        PRINT *, "** ERROR in SUBROUTINE get_eos_parameters!", &
                 " The EOS on star 1 is unknown! FUKA EOS type=", &
                 this% eos_type_2
        STOP

      ENDIF

    ENDIF

  END PROCEDURE get_eos_parameters


END SUBMODULE access
