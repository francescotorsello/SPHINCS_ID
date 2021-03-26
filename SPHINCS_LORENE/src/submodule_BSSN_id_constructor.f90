! File:         submodule_BSSN_id_constructor.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (formul_bssn_id) bssn_id_constructor

  !************************************************
  !                                               *
  ! Implementation of the constructor of TYPE     *
  ! bssn_id, and the methods it calls             *
  !                                               *
  ! FT 23.10.2020                                 *
  !                                              *
  ! Updated to mesh refinement                   *
  !                                              *
  ! FT 26.03.2021                                *
  !                                               *
  !************************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE construct_bssn_id_bns

    !***************************************************
    !                                                  *
    ! This constructor of TYPE bssn_id calls the       *
    ! SUBROUTINES that rely on an bns object, and      *
    ! allocates memory. It constructs the grid         *
    ! using the number of grid points along each axis. *
    !                                                  *
    ! FT 23.10.2020                                    *
    !                                                  *
    !***************************************************

    USE McLachlan_refine, ONLY: initialize_BSSN, deallocate_BSSN

    IMPLICIT NONE

    ! Initialize the timer
    bssn_obj% bssn_computer_timer= timer( "bssn_computer_timer" )

    ! Construct the gravity grid and import the LORENE ID on it,
    ! in standard 3+1 formulation
    CALL bssn_obj% construct_formul_3p1( bns_obj )

    ! Read and store the BSSN parameters
    !CALL bssn_obj% set_up_bssn()
    CALL initialize_BSSN()
    CALL deallocate_BSSN()

    ! The construct_formul_3p1 SUBROUTINE constructs the grid,
    ! hence the dimensions of the arrays imported from the module BSSN
    ! are know and the arrays can be allocated
    CALL allocate_bssn_fields( bssn_obj )

    PRINT *
    PRINT *, " * Ready to compute BSSN variables."
    PRINT *

  END PROCEDURE construct_bssn_id_bns


 ! MODULE PROCEDURE construct_bssn_id_bns_spacings
 !
 !   !************************************************
 !   !                                               *
 !   ! This constructor of TYPE bssn_id calls the    *
 !   ! SUBROUTINES that rely on an bns object, and   *
 !   ! allocates memory. It constructs the grid      *
 !   ! using the grid spacings.                      *
 !   !                                               *
 !   ! FT                                            *
 !   !                                               *
 !   !************************************************
 !
 !   USE McLachlan_refine, ONLY: initialize_BSSN
 !
 !   IMPLICIT NONE
 !
 !   ! Initialize the timer
 !   bssn_obj% bssn_computer_timer= timer( "bssn_computer_timer" )
 !
 !   ! Construct the gravity grid and import the LORENE ID on it,
 !   ! in standard 3+1 formulation
 !   CALL bssn_obj% construct_formul_3p1( bns_obj, dx, dy, dz )
 !
 !   ! The construct_formul_3p1 SUBROUTINE constructs the grid,
 !   ! hence the dimensions of the arrays imported from the module BSSN
 !   ! are know and the arrays can be allocated
 !   CALL allocate_bssn_fields( bssn_obj )  ! NOPASS to this SUBROUTINE?
 !
 !   ! Read and store the BSSN parameters
 !   !CALL bssn_obj% set_up_bssn()
 !   CALL initialize_BSSN()
 !
 !   PRINT *
 !   PRINT *, " * Ready to compute BSSN variables."
 !   PRINT *
 !
 ! END PROCEDURE construct_bssn_id_bns_spacings


  MODULE PROCEDURE allocate_bssn_fields

    !***********************************************
    !                                              *
    ! Allocate memory for the BSSN variables.      *
    !                                              *
    ! FT 23.10.2020                                *
    !                                              *
    ! Updated to mesh refinement                   *
    !                                              *
    ! FT 26.03.2021                                *
    !                                              *
    !***********************************************

    USE mesh_refinement,  ONLY: allocate_grid_function

    IMPLICIT NONE

    IF( .NOT.ALLOCATED( THIS% Gamma_u% levels ) )THEN
      CALL allocate_grid_function( THIS% Gamma_u, "Gamma_u_id", 3 )
    ENDIF

    IF( .NOT.ALLOCATED( THIS% phi% levels ) )THEN
      CALL allocate_grid_function( THIS% phi, "phi_id", 1 )
    ENDIF

    IF( .NOT.ALLOCATED( THIS% trK% levels ) )THEN
      CALL allocate_grid_function( THIS% trK, "trK_id", 1 )
    ENDIF

    IF( .NOT.ALLOCATED( THIS% A_BSSN3_ll% levels ) )THEN
      CALL allocate_grid_function( THIS% A_BSSN3_ll, "A_BSSN3_ll_id", 6 )
    ENDIF

    IF( .NOT.ALLOCATED( THIS% g_BSSN3_ll% levels ) )THEN
      CALL allocate_grid_function( THIS% g_BSSN3_ll, "g_BSSN3_ll_id", 6 )
    ENDIF

  END PROCEDURE allocate_bssn_fields


END SUBMODULE bssn_id_constructor
