! File:         submodule_bnslorene_access.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (bns_lorene) bns_lorene_access

  !***************************************************
  !
  !# The module contains the implementation of the
  !  methods of TYPE bns that allow to access PRIVATE
  !  members.
  !
  !  FT 12.07.2021
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
    !  FT
    !
    !***********************************************

    IMPLICIT NONE

    select_field: SELECT CASE( field )

    CASE( "lapse" )

      field_array= THIS% lapse

    CASE( "shift_x" )

      field_array= THIS% shift_x

    CASE( "shift_y" )

      field_array= THIS% shift_y

    CASE( "shift_z" )

      field_array= THIS% shift_z

    CASE( "g_xx" )

      field_array= THIS% g_xx

    CASE( "g_xy" )

      field_array= THIS% g_xy

    CASE( "g_xz" )

      field_array= THIS% g_xz

    CASE( "g_yy" )

      field_array= THIS% g_yy

    CASE( "g_yz" )

      field_array= THIS% g_yz

    CASE( "g_zz" )

      field_array= THIS% g_zz

    CASE( "k_xx" )

      field_array= THIS% k_xx

    CASE( "k_xy" )

      field_array= THIS% k_xy

    CASE( "k_xz" )

      field_array= THIS% k_xz

    CASE( "k_yy" )

      field_array= THIS% k_yy

    CASE( "k_yz" )

      field_array= THIS% k_yz

    CASE( "k_zz" )

      field_array= THIS% k_zz

    CASE( "baryon_density" )

      field_array= THIS% baryon_density

    CASE( "energy_density" )

      field_array= THIS% energy_density

    CASE( "specific_energy" )

      field_array= THIS% specific_energy

    CASE( "v_euler_x" )

      field_array= THIS% v_euler_x

    CASE( "v_euler_y" )

      field_array= THIS% v_euler_y

    CASE( "v_euler_z" )

      field_array= THIS% v_euler_z

    CASE DEFAULT

      PRINT *, "** There is no field named ", field, "in TYPE bns."
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
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    select_field: SELECT CASE( field )

    CASE( "lapse" )

      field_value= THIS% lapse( n )

    CASE( "shift_x" )

      field_value= THIS% shift_x( n )

    CASE( "shift_y" )

      field_value= THIS% shift_y( n )

    CASE( "shift_z" )

      field_value= THIS% shift_z( n )

    CASE( "g_xx" )

      field_value= THIS% g_xx( n )

    CASE( "g_xy" )

      field_value= THIS% g_xy( n )

    CASE( "g_xz" )

      field_value= THIS% g_xz( n )

    CASE( "g_yy" )

      field_value= THIS% g_yy( n )

    CASE( "g_yz" )

      field_value= THIS% g_yz( n )

    CASE( "g_zz" )

      field_value= THIS% g_zz( n )

    CASE( "k_xx" )

      field_value= THIS% k_xx( n )

    CASE( "k_xy" )

      field_value= THIS% k_xy( n )

    CASE( "k_xz" )

      field_value= THIS% k_xz( n )

    CASE( "k_yy" )

      field_value= THIS% k_yy( n )

    CASE( "k_yz" )

      field_value= THIS% k_yz( n )

    CASE( "k_zz" )

      field_value= THIS% k_zz( n )

    CASE( "baryon_density" )

      field_value= THIS% baryon_density( n )

    CASE( "energy_density" )

      field_value= THIS% energy_density( n )

    CASE( "specific_energy" )

      field_value= THIS% specific_energy( n )

    CASE( "v_euler_x" )

      field_value= THIS% v_euler_x( n )

    CASE( "v_euler_y" )

      field_value= THIS% v_euler_y( n )

    CASE( "v_euler_z" )

      field_value= THIS% v_euler_z( n )

    CASE DEFAULT

      PRINT *, "** There is no field named ", field, "in TYPE bns."
      STOP

    END SELECT select_field

  END PROCEDURE get_field_value


  MODULE PROCEDURE get_bns_identifier

    !************************************************
    !
    !# Returns the value of [[bns_identifier]], the
    !  integer identifier of the bns object
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    get_bns_identifier= THIS% bns_identifier

  END PROCEDURE get_bns_identifier


  !MODULE PROCEDURE get_bns_ptr
  !
  !  !************************************************
  !  !
  !  !# Returns the value of [[bns_ptr]], the C pointer
  !  ! to the |lorene|'s Bin_NS object
  !  ! N.B. This variable is global. The pointer
  !  !      to the second |lorene| Bin_NS object will
  !  !      overwrite the first one, and so on.
  !  !      This variable stores the pointer to
  !  !      the last defined |lorene| Bin_NS object.
  !  !      That's why it is not freed in the
  !  !      destructor of a bns object. Presently, it
  !  !      has to be freed by the user at the end of
  !  !      the PROGRAM. See the last part of the
  !  !      PROGRAM in setup_lorene_id.f90, for
  !  !      example.
  !  !
  !  !  FT
  !  !
  !  !************************************************
  !
  !  IMPLICIT NONE
  !
  !  get_bns_ptr= THIS% bns_ptr
  !
  !END PROCEDURE get_bns_ptr


  MODULE PROCEDURE get_eos1_loreneid

    !**************************************************
    !
    !# Returns the |lorene| ID-number of the EOS for NS 1
    !
    !  FT
    !
    !**************************************************

    IMPLICIT NONE

    get_eos1_loreneid= THIS% eos1_loreneid

  END PROCEDURE get_eos1_loreneid


  MODULE PROCEDURE get_eos2_loreneid

    !**************************************************
    !
    !# Returns the |lorene| ID-number of the EOS for NS 2
    !
    !  FT
    !
    !**************************************************

    IMPLICIT NONE

    get_eos2_loreneid= THIS% eos2_loreneid

  END PROCEDURE get_eos2_loreneid


END SUBMODULE bns_lorene_access
