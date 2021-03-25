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

    USE McLachlan_refine, ONLY: initialize_BSSN

    IMPLICIT NONE

    ! Initialize the timer
    bssn_obj% bssn_computer_timer= timer( "bssn_computer_timer" )

    ! Construct the gravity grid and import the LORENE ID on it,
    ! in standard 3+1 formulation
    CALL bssn_obj% construct_formul_3p1( bns_obj )

    ! The construct_formul_3p1 SUBROUTINE constructs the grid,
    ! hence the dimensions of the arrays imported from the module BSSN
    ! are know and the arrays can be allocated
    CALL allocate_bssn_fields( bssn_obj )  ! NOPASS to this SUBROUTINE?

    ! Read and store the BSSN parameters
    !CALL bssn_obj% set_up_bssn()
    CALL initialize_BSSN()

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


!  MODULE PROCEDURE set_up_bssn
!
!    !************************************************
!    !                                               *
!    ! Read and store the BSSN parameters            *
!    !                                               *
!    ! FT                                            *
!    !                                               *
!    !************************************************
!
!    USE BSSN_parameters,  ONLY: Read_BSSN_Parameters, Set_BSSN_Parameters, &
!                                fdorder, formulation, conformal_method, &
!                                harmonic_n, harmonic_f, &
!                                advect_lapse, advect_shift, eps_diss, &
!                                eps_diss_max, &
!                                evolve_a, evolve_b, alpha_driver, &
!                                damp_k1, damp_k2, fix_advection_terms, &
!                                gamma_shift, use_spatial_beta_driver, &
!                                beta_driver, spatial_beta_driver_radius, &
!                                use_spatial_shift_gamma_coeff, &
!                                shift_gamma_coeff, &
!                                spatial_shift_gamma_coeff_radius, &
!                                shift_formulation, shift_alpha_power, &
!                                minimum_lapse, lapse_threshold, &
!                                fermi_width, use_spatial_eps
!    !USE McLachlan_refine, ONLY: ghost_size
!
!    IMPLICIT NONE
!
!    LOGICAL, SAVE :: first = .TRUE.
!
!    ! Read the BSSN parameters from the parameter file, set the C++
!    ! version of the parameters, set the ghost
!    ! size depending on the finite-differencing order
!    PRINT *, " * Reading BSSN parameters..."
!    PRINT *
!
!    IF( first )THEN
!      CALL Read_BSSN_Parameters()
!      CALL Set_BSSN_Parameters( fdorder, formulation, conformal_method, &
!                                harmonic_n, harmonic_f, &
!                                advect_lapse, advect_shift, eps_diss, &
!                                eps_diss_max, &
!                                evolve_a, evolve_b, alpha_driver, &
!                                damp_k1, damp_k2, fix_advection_terms, &
!                                gamma_shift, use_spatial_beta_driver, &
!                                beta_driver, spatial_beta_driver_radius, &
!                                use_spatial_shift_gamma_coeff, &
!                                shift_gamma_coeff, &
!                                spatial_shift_gamma_coeff_radius, &
!                                shift_formulation, shift_alpha_power, &
!                                minimum_lapse, lapse_threshold, &
!                                fermi_width, use_spatial_eps )
!      first = .FALSE.
!    ENDIF
!
!    !SELECT CASE( fdorder )
!    !CASE(4)
!    !    THIS% levels% ghost_size= 3
!    !CASE(6)
!    !    THIS% levels% ghost_size= 4
!    !CASE(8)
!    !    THIS% levels% ghost_size= 5
!    !CASE DEFAULT
!    !    PRINT *, 'Invalid value chosen for fdorder chosen.'
!    !    PRINT *, 'Valued values are: 4, 6, 8'
!    !    STOP
!    !END SELECT
!
!  END PROCEDURE set_up_bssn


  MODULE PROCEDURE allocate_bssn_fields

    !***********************************************
    !                                              *
    ! Allocate memory for the BSSN variables.      *
    !                                              *
    ! FT 23.10.2020                                *
    !                                              *
    !***********************************************

    !USE tensor,    ONLY: n_sym3x3
    USE mesh_refinement,  ONLY: allocate_grid_function

    IMPLICIT NONE

    CALL allocate_grid_function( THIS% Gamma_u, "myGamma_u", 3 )

    CALL allocate_grid_function( THIS% phi, "myphi", 1 )

    CALL allocate_grid_function( THIS% trK, "mytrK", 1 )

    CALL allocate_grid_function( THIS% A_BSSN3_ll, "myA_BSSN3_ll", 6 )

    CALL allocate_grid_function( THIS% g_BSSN3_ll, "myg_BSSN3_ll", 6 )

    !IF(.NOT.ALLOCATED( THIS% Gamma_u ))THEN
    !  ALLOCATE( THIS% Gamma_u( THIS% ngrid_x, THIS% ngrid_y, &
    !                           THIS% ngrid_z, 3 ), &
    !                           STAT= ios, ERRMSG= err_msg )
    !  IF( ios > 0 )THEN
    !    PRINT *, "...allocation error for array Gamma_u ", &
    !             ". The error message is", err_msg
    !    STOP
    !  ENDIF
    !  !CALL test_status( ios, err_msg, "...allocation error for array Gamma_u" )
    !ENDIF
    !IF(.NOT.ALLOCATED( THIS% phi ))THEN
    !  ALLOCATE( THIS% phi( THIS% ngrid_x, THIS% ngrid_y, &
    !                       THIS% ngrid_z ), &
    !                       STAT= ios, ERRMSG= err_msg )
    !  IF( ios > 0 )THEN
    !    PRINT *, "...allocation error for array phi ", &
    !             ". The error message is", err_msg
    !    STOP
    !  ENDIF
    !  !CALL test_status( ios, err_msg, "...allocation error for array phi" )
    !ENDIF
    !IF(.NOT.ALLOCATED( THIS% trK ))THEN
    !  ALLOCATE( THIS% trK( THIS% ngrid_x, THIS% ngrid_y, &
    !                       THIS% ngrid_z ), &
    !                       STAT= ios, ERRMSG= err_msg )
    !  IF( ios > 0 )THEN
    !    PRINT *, "...allocation error for array trK ", &
    !             ". The error message is", err_msg
    !    STOP
    !  ENDIF
    !  !CALL test_status( ios, err_msg, "...allocation error for array trK" )
    !ENDIF
    !IF(.NOT.ALLOCATED( THIS% A_BSSN3_ll ))THEN
    !  ALLOCATE( THIS% A_BSSN3_ll( THIS% ngrid_x, THIS% ngrid_y, &
    !                              THIS% ngrid_z, n_sym3x3 ), &
    !                              STAT= ios, ERRMSG= err_msg )
    !  IF( ios > 0 )THEN
    !    PRINT *, "...allocation error for array A_BSSN3_ll ", &
    !             ". The error message is", err_msg
    !    STOP
    !  ENDIF
    !  !CALL test_status( ios, err_msg, &
    !  !                "...allocation error for array A_BSSN3_ll" )
    !ENDIF
    !IF(.NOT.ALLOCATED( THIS% g_BSSN3_ll ))THEN
    !  ALLOCATE( THIS% g_BSSN3_ll( THIS% ngrid_x, THIS% ngrid_y, &
    !                              THIS% ngrid_z, n_sym3x3 ), &
    !                              STAT= ios, ERRMSG= err_msg )
    !  IF( ios > 0 )THEN
    !    PRINT *, "...allocation error for array g_BSSN3_ll ", &
    !             ". The error message is", err_msg
    !    STOP
    !  ENDIF
    !  !CALL test_status( ios, err_msg, &
    !  !                "...allocation error for array g_BSSN3_ll" )
    !ENDIF

  END PROCEDURE allocate_bssn_fields


END SUBMODULE bssn_id_constructor
