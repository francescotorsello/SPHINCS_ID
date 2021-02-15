! File:         module_bssn_id.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

MODULE formul_bssn_id

  !***********************************************************
  !                                                          *
  !   This module contains the definition of TYPE bssn_id    *
  !                                                          *
  !***********************************************************


  USE utility,       ONLY: ios, err_msg, perc, creturn, run_id, test_status, &
                           compute_g4, determinant_sym4x4_grid, show_progress
  USE bns_id,        ONLY: bns
  USE formul_3p1_id, ONLY: formul_3p1
  USE particles_id,  ONLY: particles
  USE timing,        ONLY: timer


  IMPLICIT NONE


  !********************************************************
  !                                                       *
  !              Definition of TYPE bssn_id               *
  !                                                       *
  ! This class extends the ABSTRACT TYPE formul_3p1 by    *
  ! implementing its deferred methods such that the BSSN  *
  ! variable sare computed on the grid for the LORENE ID, *
  ! stored, exported to a binary file for evolution and   *
  ! to a formatted file. The BSSN constraints can also    *
  ! be computed in different ways, analyzed, and exported *
  ! in different ways.                                    *
  !                                                       *
  !********************************************************

  TYPE, EXTENDS( formul_3p1 ):: bssn_id

    ! The flag call_flag is set different than 0 if the SUBROUTINE
    ! compute_and_export_bssn_variables is called
    INTEGER:: call_flag= 0

    !
    !-- Arrays storing the BSSN variables for the LORENE ID on the grid
    !
    ! Conformal connection
    DOUBLE PRECISION, ALLOCATABLE :: Gamma_u(:,:,:,:)
    ! Conformal factor
    DOUBLE PRECISION, ALLOCATABLE :: phi(:,:,:)
    ! Trace of extrinsic curvature
    DOUBLE PRECISION, ALLOCATABLE :: trK(:,:,:)
    ! Conformal traceless extrinsic curvature
    DOUBLE PRECISION, ALLOCATABLE :: A_BSSN3_ll(:,:,:,:)
    ! Conformal metric
    DOUBLE PRECISION, ALLOCATABLE :: g_BSSN3_ll(:,:,:,:)

    !
    !-- Connection constraints and its l2 norm
    !
    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: GC
    DOUBLE PRECISION, DIMENSION(3):: GC_l2
    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: GC_parts
    DOUBLE PRECISION, DIMENSION(3):: GC_parts_l2

    LOGICAL, PUBLIC:: export_bin
    LOGICAL, PUBLIC:: export_form_xy, export_form_x
    LOGICAL, PUBLIC:: export_constraints_xy, export_constraints_x

    TYPE(timer):: bssn_computer_timer


    CONTAINS


    !-------------------!
    !--  SUBROUTINES  --!
    !-------------------!

    PROCEDURE :: set_up_bssn

    PROCEDURE :: define_allocate_fields => allocate_bssn_fields

    PROCEDURE :: compute_and_export_3p1_variables &
                    => compute_and_export_bssn_variables

    PROCEDURE, PUBLIC :: read_bssn_dump_print_formatted

    PROCEDURE :: print_formatted_lorene_id_3p1_variables &
                    => print_formatted_lorene_id_bssn_variables

    PROCEDURE :: compute_and_export_3p1_constraints_grid &
                    => compute_and_export_bssn_constraints_grid

    PROCEDURE :: compute_and_export_3p1_constraints_particles &
                    => compute_and_export_bssn_constraints_particles

    PROCEDURE :: deallocate_fields => deallocate_bssn_fields

    ! Finalizer for members of the extended class bssn_id, not the
    ! primitive class formul_3p1
    PROCEDURE :: destruct_bssn_id

    ! Destructor; finalizes members from both CLASSES formul_3p1, and bssn_id,
    ! by calling destruct_formul_3p1 and destruct_bssn_id
    FINAL     :: destructor

  END TYPE bssn_id

  !
  !-- Interface of the TYPE bssn_id
  !-- (i.e., declaration of the overloaded constructor)
  !
  INTERFACE bssn_id

    ! Constructs the bssn_id object from the number of grid points
    ! along each axis
    MODULE PROCEDURE:: construct_bssn_id_bns
    ! Constructs the bssn_id object from the grid spacings
    MODULE PROCEDURE:: construct_bssn_id_bns_spacings

  END INTERFACE bssn_id

  !
  !-- Interface of the overloaded constructor of TYPE bssn_id
  !-- Its implementation is in submodule_BSSN_id_constructor.f90
  !
  INTERFACE

    MODULE FUNCTION construct_bssn_id_bns( bns_obj ) RESULT ( bssn_obj )

      CLASS(bns), INTENT( IN OUT ):: bns_obj
      TYPE(bssn_id)               :: bssn_obj

    END FUNCTION construct_bssn_id_bns

    MODULE FUNCTION construct_bssn_id_bns_spacings( bns_obj, dx, dy, dz ) &
                    RESULT ( bssn_obj )

      CLASS(bns), INTENT( IN OUT )  :: bns_obj
      TYPE(bssn_id)                 :: bssn_obj
      DOUBLE PRECISION, INTENT( IN ):: dx, dy, dz

    END FUNCTION construct_bssn_id_bns_spacings

  END INTERFACE

  !
  !-- Interfaces of the methods of TYPE bssn_id
  !-- Their implementations are in submodule_BSSN_id_methods.f90
  !
  INTERFACE

    MODULE SUBROUTINE set_up_bssn( THIS )

      CLASS(bssn_id), INTENT( IN OUT ):: THIS

    END SUBROUTINE set_up_bssn

    MODULE SUBROUTINE allocate_bssn_fields( THIS )

      CLASS(bssn_id), INTENT( IN OUT ):: THIS

    END SUBROUTINE allocate_bssn_fields

    MODULE SUBROUTINE compute_and_export_bssn_variables( THIS, namefile )

      CLASS(bssn_id),      INTENT( IN OUT )           :: THIS
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE compute_and_export_bssn_variables

    MODULE SUBROUTINE read_bssn_dump_print_formatted( THIS, namefile_bin, &
                                                            namefile )

      CLASS(bssn_id),      INTENT( IN OUT )           :: THIS
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile_bin
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE read_bssn_dump_print_formatted

    MODULE SUBROUTINE print_formatted_lorene_id_bssn_variables( THIS, &
                                                                namefile )

      CLASS(bssn_id),      INTENT( IN OUT )           :: THIS
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE print_formatted_lorene_id_bssn_variables

    MODULE SUBROUTINE compute_and_export_bssn_constraints_grid( THIS, &
                                                           bns_obj, &
                                                           namefile, &
                                                           name_logfile )

      CLASS(bssn_id),      INTENT( IN OUT ):: THIS
      CLASS(bns),          INTENT( IN OUT ):: bns_obj
      CHARACTER( LEN= * ), INTENT( IN OUT ):: namefile
      CHARACTER( LEN= * ), INTENT( IN OUT ):: name_logfile

    END SUBROUTINE compute_and_export_bssn_constraints_grid

    MODULE SUBROUTINE compute_and_export_bssn_constraints_particles( THIS, &
                                                           parts_obj, &
                                                           namefile, &
                                                           namefile_sph, &
                                                           name_logfile )

      CLASS(bssn_id),      INTENT( IN OUT ):: THIS
      CLASS(particles),    INTENT( IN OUT ):: parts_obj
      CHARACTER( LEN= * ), INTENT( IN OUT ):: namefile
      CHARACTER( LEN= * ), INTENT( IN OUT ):: namefile_sph
      CHARACTER( LEN= * ), INTENT( IN OUT ):: name_logfile

    END SUBROUTINE compute_and_export_bssn_constraints_particles

    MODULE SUBROUTINE deallocate_bssn_fields( THIS )

      CLASS(bssn_id), INTENT( IN OUT ):: THIS

    END SUBROUTINE deallocate_bssn_fields

    MODULE SUBROUTINE destruct_bssn_id( THIS )

      CLASS(bssn_id), INTENT( IN OUT ):: THIS

    END SUBROUTINE destruct_bssn_id

    MODULE SUBROUTINE destructor( THIS )

      TYPE(bssn_id), INTENT( IN OUT ):: THIS

    END SUBROUTINE destructor

  END INTERFACE


END MODULE formul_bssn_id
