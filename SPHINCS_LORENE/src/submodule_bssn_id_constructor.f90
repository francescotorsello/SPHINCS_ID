! File:         submodule_BSSN_id_constructor.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (formul_bssn_id) bssn_id_constructor

  !************************************************
  !                                               *
  ! Implementation of the constructor and         *
  ! destructor of TYPE bssn_id                    *
  !                                               *
  ! FT 23.10.2020                                 *
  !                                               *
  ! Updated to mesh refinement                    *
  !                                               *
  ! FT 26.03.2021                                 *
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
    USE mesh_refinement,  ONLY: levels
    USE Extract_Mass,     ONLY: radius2

    IMPLICIT NONE

    ! Initialize the timer
    bssn_obj% bssn_computer_timer= timer( "bssn_computer_timer" )

    ! Construct the gravity grid and import the LORENE ID on it,
    ! in standard 3+1 formulation
    IF( PRESENT(dx) .AND. PRESENT(dy) .AND. PRESENT(dz) )THEN

      CALL bssn_obj% construct_formul_3p1( bns_obj, dx, dy, dz )

    ELSE

      CALL bssn_obj% construct_formul_3p1( bns_obj )

    ENDIF

    ! Read and store the BSSN parameters
    !CALL bssn_obj% set_up_bssn()
    CALL initialize_BSSN()
    CALL deallocate_BSSN()

    ! The construct_formul_3p1 SUBROUTINE constructs the grid,
    ! hence the dimensions of the arrays imported from the module BSSN
    ! are know and the arrays can be allocated
    CALL allocate_bssn_fields( bssn_obj )

    DEALLOCATE( levels )

    ! radius2 is the extraction radius. If not set here, then it is 0 by default
    ! and the metric is not interpolate on the particle in
    ! get_metric_on_particles
    radius2= HUGE(DBLE(1.0D0))

    PRINT *
    PRINT *, " * Ready to compute BSSN variables."
    PRINT *

  END PROCEDURE construct_bssn_id_bns


!  MODULE PROCEDURE construct_bssn_id_bns_spacings
!
!    !************************************************
!    !                                               *
!    ! This constructor of TYPE bssn_id calls the    *
!    ! SUBROUTINES that rely on an bns object, and   *
!    ! allocates memory. It constructs the grid      *
!    ! using the grid spacings.                      *
!    !                                               *
!    ! FT 21.06.2021                                 *
!    !                                               *
!    !************************************************
!
!    USE McLachlan_refine, ONLY: initialize_BSSN
!
!    IMPLICIT NONE
!
!    ! Initialize the timer
!    bssn_obj% bssn_computer_timer= timer( "bssn_computer_timer" )
!
!    ! Construct the gravity grid and import the LORENE ID on it,
!    ! in standard 3+1 formulation
!    CALL bssn_obj% construct_formul_3p1( bns_obj, dx, dy, dz )
!
!    ! Read and store the BSSN parameters
!    !CALL bssn_obj% set_up_bssn()
!    CALL initialize_BSSN()
!    CALL deallocate_BSSN()
!
!    ! The construct_formul_3p1 SUBROUTINE constructs the grid,
!    ! hence the dimensions of the arrays imported from the module BSSN
!    ! are know and the arrays can be allocated
!    CALL allocate_bssn_fields( bssn_obj )
!
!    PRINT *
!    PRINT *, " * Ready to compute BSSN variables."
!    PRINT *
!
!  END PROCEDURE construct_bssn_id_bns_spacings

  !
  !-- Keeping the following two SUBROUTINES separate in case it is needed
  !-- to add other PROCEDURES to the destructor (probably superfluous...)
  !
  MODULE PROCEDURE destruct_bssn_id

    !**************************************************
    !                                                 *
    ! Finalizer for members of the extended class     *
    ! bssn_id, not the primitive class formul_3p1     *
    !                                                 *
    ! FT                                              *
    !                                                 *
    !**************************************************

    IMPLICIT NONE

    !PRINT *, "Inside destruct_BSSN_id"
    CALL deallocate_bssn_fields( THIS )

  END PROCEDURE destruct_bssn_id


  MODULE PROCEDURE destructor

    !**************************************************
    !                                                 *
    ! Destructor of TYPE bssn_id                      *
    !                                                 *
    ! FT                                              *
    !                                                 *
    !**************************************************

    IMPLICIT NONE

    !PRINT *, "Inside destructor"
    CALL destruct_bssn_id( THIS )

#ifdef __INTEL_COMPILER
  CALL destruct_formul_3p1( THIS )
#endif

#ifdef __GFORTRAN__
  CALL THIS% destruct_formul_3p1
#endif

  END PROCEDURE destructor


END SUBMODULE bssn_id_constructor
