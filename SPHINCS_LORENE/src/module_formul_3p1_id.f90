! File:         module_formul_3p1_id.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

MODULE formul_3p1_id

  !**********************************************************************
  !                                                                     *
  !   This module contains the definition of ABSTRACT TYPE formul_3p1   *
  !                                                                     *
  !**********************************************************************


  USE utility,          ONLY: ios, err_msg, perc, creturn, run_id, &
                              test_status, show_progress
  USE bns_id,           ONLY: bns
  USE particles_id,     ONLY: particles
  USE timing,           ONLY: timer
  USE mesh_refinement,  ONLY: grid_function_scalar, grid_function, level


  IMPLICIT NONE


  !*******************************************************
  !                                                      *
  !      Definition of abstract TYPE formul_3p1          *
  !                                                      *
  ! Abstract class for a 3+1 formulation of the Einstein *
  ! equations. It imports the LORENE ID on the gravity   *
  ! grid, in the standard 3+1 formulation, and defines   *
  ! DEFERRED PROCEDURES to be implemented in the derived *
  ! TYPES of the actual 3+1 formulations (for example,   *
  ! the BSSN formulation).                               *
  !                                                      *
  !*******************************************************

  TYPE, ABSTRACT:: formul_3p1

    INTEGER, PUBLIC:: cons_step

    !INTEGER:: ngrid_x, ngrid_y, ngrid_z

    !DOUBLE PRECISION:: dx, dy, dz, dx_1, dy_1, dz_1

    ! Array storing the Cartesian coordinates of the grid points
    !DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: grid

    !
    !-- Arrays storing the 3+1 variables for the LORENE ID on the grid
    !
    !DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE:: lapse
    !DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: shift_u
    !DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: g_phys3_ll
    !DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: K_phys3_ll
    !
    !! Constraints and their l2 norm
    !DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE:: HC
    !DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: MC
    !DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE:: HC_parts
    !DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: MC_parts

    ! Variables to decide if and how to export the constraints
    LOGICAL, PUBLIC:: export_constraints
    LOGICAL, PUBLIC:: export_constraints_details

    !
    !-- New variables for the refined mesh
    !
    INTEGER:: nlevels

    TYPE(timer):: grid_timer
    TYPE(timer):: importer_timer

    TYPE(level), DIMENSION(:), ALLOCATABLE :: levels
    TYPE(grid_function):: coords

    TYPE(grid_function_scalar):: lapse
    TYPE(grid_function):: shift_u
    TYPE(grid_function):: g_phys3_ll
    TYPE(grid_function):: K_phys3_ll

    TYPE(grid_function_scalar):: HC
    TYPE(grid_function_scalar):: HC_parts
    TYPE(grid_function):: MC
    TYPE(grid_function):: MC_parts

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: HC_l2
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: MC_l2
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: HC_parts_l2
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: MC_parts_l2


    CONTAINS


    !-------------------!
    !--  SUBROUTINES  --!
    !-------------------!

    GENERIC, PUBLIC:: construct_formul_3p1 => &
                                        construct_formul_3p1_bns_ptr, &
                                        construct_formul_3p1_bns_spacings_ptr

    PROCEDURE::       construct_formul_3p1_bns_ptr      => &
                                        construct_formul_3p1_bns

    PROCEDURE::       construct_formul_3p1_bns_spacings_ptr => &
                                        construct_formul_3p1_bns_spacings

    PROCEDURE:: analyze_constraint

    PROCEDURE(define_allocate_fields_interface), DEFERRED:: &
                            define_allocate_fields

    PROCEDURE(compute_and_export_3p1_variables_interface), PUBLIC, &
                            DEFERRED:: compute_and_export_3p1_variables

    PROCEDURE(print_formatted_lorene_id_3p1_variables_interface), PUBLIC, &
                            DEFERRED:: print_formatted_lorene_id_3p1_variables

    GENERIC, PUBLIC:: compute_and_export_3p1_constraints => &
                      compute_and_export_3p1_constraints_grid, &
                      compute_and_export_3p1_constraints_particles

    PROCEDURE(compute_and_export_3p1_constraints_grid_interface), &
              DEFERRED:: compute_and_export_3p1_constraints_grid

    PROCEDURE(compute_and_export_3p1_constraints_particles_interface), &
              DEFERRED:: compute_and_export_3p1_constraints_particles

    PROCEDURE(deallocate_fields_interface), DEFERRED:: deallocate_fields

    PROCEDURE:: destruct_formul_3p1

    !-----------------!
    !--  FUNCTIONS  --!
    !-----------------!

    PROCEDURE:: abs_values_in

    PROCEDURE, PUBLIC:: get_grid_point

    PROCEDURE, PUBLIC:: get_x_spacing

    PROCEDURE, PUBLIC:: get_ngrid_x

    PROCEDURE, PUBLIC:: get_ngrid_y

    PROCEDURE, PUBLIC:: get_ngrid_z

    PROCEDURE, PUBLIC:: get_HC

    PROCEDURE, PUBLIC:: get_MC

    PROCEDURE, PUBLIC:: get_HC_parts

    PROCEDURE, PUBLIC:: get_MC_parts

  END TYPE formul_3p1

  !
  !-- Interface of the cores of the constructors and destructos of TYPES
  !-- derived from formul_3p1
  !-- Their implementations are in submodule formul_3p1_methods.f90
  !
  INTERFACE

    MODULE SUBROUTINE construct_formul_3p1_bns( f3p1_obj, bns_obj )

      CLASS(bns),        INTENT( IN OUT ):: bns_obj
      CLASS(formul_3p1), INTENT( IN OUT ):: f3p1_obj

    END SUBROUTINE construct_formul_3p1_bns

    MODULE SUBROUTINE construct_formul_3p1_bns_spacings( f3p1_obj, bns_obj, &
                                                         dx, dy, dz )

      CLASS(bns),        INTENT( IN OUT ):: bns_obj
      CLASS(formul_3p1), INTENT( IN OUT ):: f3p1_obj
      DOUBLE PRECISION,  INTENT( IN )    :: dx, dy, dz

    END SUBROUTINE construct_formul_3p1_bns_spacings

    MODULE SUBROUTINE destruct_formul_3p1( f3p1_obj )

      CLASS(formul_3p1), INTENT( IN OUT ):: f3p1_obj

    END SUBROUTINE destruct_formul_3p1

  END INTERFACE

  !
  !-- Interface of the methods of TYPES derived from formul_3p1
  !-- Their implementations are in submodule formul_3p1_methods.f90
  !
  INTERFACE

    MODULE SUBROUTINE analyze_constraint( THIS, &
                                          nx, ny, nz, &
                                          constraint, &
                                          name_constraint, &
                                          name_analysis, &
                                          unit_logfile, &
                                          l2_norm, &
                                          loo_norm )

      CLASS(formul_3p1),                  INTENT( IN OUT ):: THIS
      INTEGER,                            INTENT( IN )    :: nx, ny, nz
      DOUBLE PRECISION, DIMENSION(:,:,:), INTENT( IN )    :: constraint
      CHARACTER( LEN= * ),                INTENT( IN )    :: name_constraint
      CHARACTER( LEN= * ),                INTENT( IN )    :: name_analysis
      INTEGER,                            INTENT( IN )    :: unit_logfile
      DOUBLE PRECISION,                   INTENT( OUT )   :: l2_norm
      DOUBLE PRECISION,                   INTENT( OUT )   :: loo_norm

    END SUBROUTINE analyze_constraint

    MODULE FUNCTION abs_values_in( THIS, lower_bound, upper_bound, &
                                   constraint, nx, ny, nz, &
                                   unit_analysis ) RESULT( cnt )

      ! Arguments
      CLASS(formul_3p1),                  INTENT( IN OUT ):: THIS
      DOUBLE PRECISION                                    :: lower_bound, &
                                                             upper_bound
      DOUBLE PRECISION, DIMENSION(:,:,:), INTENT( IN )    :: constraint
      INTEGER,                            INTENT( IN )    :: nx, ny, nz
      INTEGER,                            INTENT( IN )    :: unit_analysis
      ! Result
      INTEGER                                             :: cnt

    END FUNCTION abs_values_in

    MODULE FUNCTION get_grid_point( THIS, i, j, k, l ) RESULT( grid_point )

      ! Arguments
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS
      INTEGER,           INTENT( IN )    :: i, j, k, l
      ! Result
      DOUBLE PRECISION, DIMENSION(3)     :: grid_point

    END FUNCTION get_grid_point

    MODULE FUNCTION get_x_spacing( THIS, l ) RESULT( dx )

      ! Arguments
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS
      INTEGER,           INTENT( IN )    :: l
      ! Result
      DOUBLE PRECISION:: dx

    END FUNCTION get_x_spacing

    MODULE FUNCTION get_ngrid_x( THIS, l ) RESULT( ngrid_x )

      ! Arguments
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS
      INTEGER,           INTENT( IN )    :: l
      ! Result
      DOUBLE PRECISION:: ngrid_x

    END FUNCTION get_ngrid_x

    MODULE FUNCTION get_ngrid_y( THIS, l ) RESULT( ngrid_y )

      ! Arguments
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS
      INTEGER,           INTENT( IN )    :: l
      ! Result
      DOUBLE PRECISION:: ngrid_y

    END FUNCTION get_ngrid_y

    MODULE FUNCTION get_ngrid_z( THIS, l ) RESULT( ngrid_z )

      ! Arguments
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS
      INTEGER,           INTENT( IN )    :: l
      ! Result
      DOUBLE PRECISION:: ngrid_z

    END FUNCTION get_ngrid_z

    MODULE FUNCTION get_HC( THIS, i, j, k, l ) RESULT( HC_value )

      ! Arguments
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS
      INTEGER,           INTENT( IN )    :: i, j, k, l
      ! Result
      DOUBLE PRECISION                   :: HC_value

    END FUNCTION get_HC

    MODULE FUNCTION get_MC( THIS, i, j, k, l ) RESULT( MC_value )

      ! Arguments
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS
      INTEGER,           INTENT( IN )    :: i, j, k, l
      ! Result
      DOUBLE PRECISION, DIMENSION(3)     :: MC_value

    END FUNCTION get_MC

    MODULE FUNCTION get_HC_parts( THIS, i, j, k, l ) RESULT( HC_value )

      ! Arguments
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS
      INTEGER,           INTENT( IN )    :: i, j, k, l
      ! Result
      DOUBLE PRECISION                   :: HC_value

    END FUNCTION get_HC_parts

    MODULE FUNCTION get_MC_parts( THIS, i, j, k, l ) RESULT( MC_value )

      ! Arguments
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS
      INTEGER,           INTENT( IN )    :: i, j, k, l
      ! Result
      DOUBLE PRECISION, DIMENSION(3)     :: MC_value

    END FUNCTION get_MC_parts

  END INTERFACE

  !
  !-- Interfaces of the deferred methods of TYPE formul_3p1
  !-- Their implementations are deferred to derived TYPES
  !
  ABSTRACT INTERFACE

    SUBROUTINE define_allocate_fields_interface( THIS )

      IMPORT:: formul_3p1
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS

    END SUBROUTINE define_allocate_fields_interface

    SUBROUTINE compute_and_export_3p1_variables_interface( THIS, &
                                                                  namefile )

      IMPORT:: formul_3p1
      CLASS(formul_3p1),   INTENT( IN OUT )           :: THIS
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE compute_and_export_3p1_variables_interface

    SUBROUTINE print_formatted_lorene_id_3p1_variables_interface &
                                                    ( THIS, namefile )

      IMPORT:: formul_3p1
      CLASS(formul_3p1),   INTENT( IN OUT )           :: THIS
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE print_formatted_lorene_id_3p1_variables_interface

    SUBROUTINE compute_and_export_3p1_constraints_grid_interface( THIS, &
                                                             bns_obj, &
                                                             namefile, &
                                                             name_logfile )

      IMPORT:: formul_3p1
      IMPORT:: bns
      CLASS(formul_3p1),   INTENT( IN OUT ):: THIS
      CLASS(bns),          INTENT( IN OUT ):: bns_obj
      CHARACTER( LEN= * ), INTENT( IN OUT ):: namefile
      CHARACTER( LEN= * ), INTENT( IN OUT ):: name_logfile

    END SUBROUTINE compute_and_export_3p1_constraints_grid_interface

    SUBROUTINE compute_and_export_3p1_constraints_particles_interface( THIS, &
                                                             parts_obj, &
                                                             namefile, &
                                                             namefile_sph, &
                                                             name_logfile )

      IMPORT:: formul_3p1
      IMPORT:: particles
      CLASS(formul_3p1),   INTENT( IN OUT ):: THIS
      CLASS(particles),    INTENT( IN OUT ):: parts_obj
      CHARACTER( LEN= * ), INTENT( IN OUT ):: namefile
      CHARACTER( LEN= * ), INTENT( IN OUT ):: namefile_sph
      CHARACTER( LEN= * ), INTENT( IN OUT ):: name_logfile

    END SUBROUTINE compute_and_export_3p1_constraints_particles_interface

    SUBROUTINE deallocate_fields_interface( THIS )

      IMPORT:: formul_3p1
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS

    END SUBROUTINE deallocate_fields_interface

  END INTERFACE


END MODULE formul_3p1_id
