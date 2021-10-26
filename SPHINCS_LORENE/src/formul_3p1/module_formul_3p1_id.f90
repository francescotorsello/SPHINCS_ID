! File:         module_formul_3p1_id.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

MODULE formul_3p1_id

  !**********************************************************************
  !
  !#   This module contains the definition of ABSTRACT TYPE formul_3p1
  !
  !**********************************************************************


  USE utility,          ONLY: ios, err_msg, perc, creturn, run_id, &
                              test_status, show_progress
  USE id_base,          ONLY: idbase
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
  !# ABSTRACT TYPE representing a generic 3+1 formulation of the Einsten
  !  equations

    !
    !-- Mesh variables
    !

    INTEGER:: nlevels
    !! Number of refinement levels

    TYPE(level), DIMENSION(:), ALLOCATABLE :: levels
    !! Array containing the information on each refinement level
    TYPE(grid_function):: coords
    !! Grid function storing the Cartesian coordinates
    TYPE(grid_function_scalar):: rad_coord
    !! Grid scalar function storing the radial coordinates of each grid point

    !
    !-- ADM fields
    !

    TYPE(grid_function_scalar):: lapse
    !! Grid scalar function storing the lapse function
    TYPE(grid_function):: shift_u
    !! Grid function storing the shift vector
    TYPE(grid_function):: g_phys3_ll
    !! Grid function storing the spatial metric
    TYPE(grid_function):: K_phys3_ll
    !! Grid function storing the extrinsic curvature

    !
    !-- Constraint violations
    !

    TYPE(grid_function_scalar):: HC
    !# Grid scalar function storing the Hamiltonian constraint (violations)
    !  computed using the full |lorene| ID on the mesh
    TYPE(grid_function_scalar):: HC_parts
    !# Grid scalar function storing the Hamiltonian constraint (violations)
    !  computed using the stress-energy tensor mapped from the particles to the
    !  mesh
    TYPE(grid_function):: MC
    !# Grid function storing the momentum constraint (violations)
    !  computed using the full |lorene| ID on the mesh
    TYPE(grid_function):: MC_parts
    !# Grid function storing the momentum constraint (violations)
    !  computed using the stress-energy tensor mapped from the particles to the
    !  mesh

    !
    !-- Norms of constraint violations
    !

    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: HC_l2
    !! \(\ell_2\) norm of the Hamiltonian constraint computed on the mesh
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: HC_parts_l2
    !# \(\ell_2\) norm of the Hamiltonian constraint computed on the mesh, using
    !  the stress-energy tensor mapped from the particles
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: HC_loo
    !# \(\ell_\infty\) norm of the Hamiltonian constraint computed on the mesh
    !  (i.e., its maximum)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: HC_parts_loo
    !# \(\ell_\infty\) norm of the Hamiltonian constraint computed on the mesh,
    !  using the stress-energy tensor mapped from the particles
    !  (i.e., its maximum)
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: MC_l2
    !! \(\ell_2\) norm of the momentum constraint computed on the mesh
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: MC_parts_l2
    !# \(\ell_2\) norm of the Hamiltonian constraint computed on the mesh, using
    !  the stress-energy tensor mapped from the particles
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: MC_loo
    !# \(\ell_\infty\) norm of the momentum constraint computed on the mesh
    !  (i.e., its maximum)
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: MC_parts_loo
    !# \(\ell_\infty\) norm of the momentum constraint computed on the mesh,
    !  using the stress-energy tensor mapped from the particles
    !  (i.e., its maximum)

    !
    !-- Steering variables
    !

    ! Variables to decide if and how to export the constraints

    INTEGER, PUBLIC:: cons_step
    !! Constraint violations are printed to file every cons_step-th grid point
    !  along each Cartesian direction
    LOGICAL, PUBLIC:: export_constraints
    !# `.TRUE.` if the constraint violations  are to be printed to file,
    !  `.FALSE.` otherwise
    LOGICAL, PUBLIC:: export_constraints_details
    !# `.TRUE.` if the points at which the constraints violations are within the
    !  intervals \([0,10^{-7}],[10^3,\infty],[10^n,10^m]\), with
    !  \(n\in\{-7,2\}\) and \(m=n+1\), are to be printed to file;
    !  `.FALSE.` otherwise
    LOGICAL, PUBLIC:: export_constraints_xy
    !! `.TRUE.` if only the constrain violations on the \(xy\) plane are to be
    !  printed to file, `.FALSE.` otherwise
    LOGICAL, PUBLIC:: export_constraints_x
    !! `.TRUE.` if only the constrain violations on the \(x\) axis are to be
    !  printed to file, `.FALSE.` otherwise

    !
    !-- Timers
    !

    TYPE(timer):: grid_timer
    !# Timer that times how long it takes to set up the grid and allocate
    !  the grid functions
    TYPE(timer):: importer_timer
    !# Timer that times how long it takes to import the |lorene| ID
    !  on the mesh


    CONTAINS


    !-------------------!
    !--  SUBROUTINES  --!
    !-------------------!

    GENERIC, PUBLIC:: construct_formul_3p1 => &
                                        construct_formul_3p1_bns_ptr!, &
                                        !construct_formul_3p1_bns_spacings_ptr

    PROCEDURE::       construct_formul_3p1_bns_ptr => &
                                        construct_formul_3p1_bns

    !PROCEDURE::       construct_formul_3p1_bns_spacings_ptr => &
    !                                    construct_formul_3p1_bns_spacings

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

    PROCEDURE, PUBLIC:: get_nlevels

    PROCEDURE, PUBLIC:: get_levels

    PROCEDURE, PUBLIC:: get_dx

    PROCEDURE, PUBLIC:: get_dy

    PROCEDURE, PUBLIC:: get_dz

    PROCEDURE, PUBLIC:: get_ngrid_x

    PROCEDURE, PUBLIC:: get_ngrid_y

    PROCEDURE, PUBLIC:: get_ngrid_z

    PROCEDURE, PUBLIC:: get_xR

    PROCEDURE, PUBLIC:: get_yR

    PROCEDURE, PUBLIC:: get_zR

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

    MODULE SUBROUTINE construct_formul_3p1_bns( f3p1, id, dx, dy, dz )

      CLASS(idbase),    INTENT( IN OUT ):: id
      CLASS(formul_3p1), INTENT( IN OUT ):: f3p1
      DOUBLE PRECISION, OPTIONAL         :: dx, dy, dz

    END SUBROUTINE construct_formul_3p1_bns

 !   MODULE SUBROUTINE construct_formul_3p1_bns_spacings( f3p1, id, &
 !                                                        dx, dy, dz )
 !
 !     CLASS(bns),        INTENT( IN OUT ):: id
 !     CLASS(formul_3p1), INTENT( IN OUT ):: f3p1
 !     DOUBLE PRECISION,  INTENT( IN )    :: dx, dy, dz
 !
 !   END SUBROUTINE construct_formul_3p1_bns_spacings

    MODULE SUBROUTINE destruct_formul_3p1( f3p1 )

      CLASS(formul_3p1), INTENT( IN OUT ):: f3p1

    END SUBROUTINE destruct_formul_3p1

  END INTERFACE

  !
  !-- Interface of the methods of TYPES derived from formul_3p1
  !-- Their implementations are in submodule formul_3p1_methods.f90
  !
  INTERFACE

    MODULE SUBROUTINE analyze_constraint( THIS, &
                                          l, &
                                          constraint, &
                                          name_constraint, &
                                          unit_logfile, &
                                          name_analysis, &
                                          l2_norm, &
                                          loo_norm )

      CLASS(formul_3p1),                  INTENT( IN OUT ):: THIS
      INTEGER,                            INTENT( IN )    :: l
      DOUBLE PRECISION, DIMENSION(:,:,:), INTENT( IN )    :: constraint
      CHARACTER( LEN= * ),                INTENT( IN )    :: name_constraint
      CHARACTER( LEN= * ),                INTENT( IN )    :: name_analysis
      INTEGER,                            INTENT( IN )    :: unit_logfile
      DOUBLE PRECISION,                   INTENT( OUT )   :: l2_norm
      DOUBLE PRECISION,                   INTENT( OUT )   :: loo_norm

    END SUBROUTINE analyze_constraint

    MODULE SUBROUTINE abs_values_in( THIS, lower_bound, upper_bound, &
                                     constraint, l, &
                                     export, unit_analysis, cnt )

      CLASS(formul_3p1),                  INTENT( IN OUT ):: THIS
      DOUBLE PRECISION                                    :: lower_bound, &
                                                             upper_bound
      DOUBLE PRECISION, DIMENSION(:,:,:), INTENT( IN )    :: constraint
      INTEGER,                            INTENT( IN )    :: l
      LOGICAL,                            INTENT( IN )    :: export
      INTEGER,                            INTENT( IN )    :: unit_analysis
      INTEGER,                            INTENT( OUT )   :: cnt

    END SUBROUTINE abs_values_in

    MODULE FUNCTION get_grid_point( THIS, i, j, k, l ) RESULT( grid_point )

      ! Arguments
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS
      INTEGER,           INTENT( IN )    :: i, j, k, l
      ! Result
      DOUBLE PRECISION, DIMENSION(3)     :: grid_point

    END FUNCTION get_grid_point

    MODULE FUNCTION get_nlevels( THIS ) RESULT( nlevels )

      ! Arguments
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS
      ! Result
      DOUBLE PRECISION:: nlevels

    END FUNCTION get_nlevels

    MODULE FUNCTION get_levels( THIS, l ) RESULT( levels )

      ! Arguments
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS
      INTEGER,           INTENT( IN )    :: l
      ! Result
      TYPE(level), DIMENSION(:), ALLOCATABLE:: levels

    END FUNCTION get_levels

    MODULE FUNCTION get_dx( THIS, l ) RESULT( dx )

      ! Arguments
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS
      INTEGER,           INTENT( IN )    :: l
      ! Result
      DOUBLE PRECISION:: dx

    END FUNCTION get_dx

    MODULE FUNCTION get_dy( THIS, l ) RESULT( dy )

      ! Arguments
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS
      INTEGER,           INTENT( IN )    :: l
      ! Result
      DOUBLE PRECISION:: dy

    END FUNCTION get_dy

    MODULE FUNCTION get_dz( THIS, l ) RESULT( dz )

      ! Arguments
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS
      INTEGER,           INTENT( IN )    :: l
      ! Result
      DOUBLE PRECISION:: dz

    END FUNCTION get_dz

    MODULE FUNCTION get_ngrid_x( THIS, l ) RESULT( ngrid_x )

      ! Arguments
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS
      INTEGER,           INTENT( IN )    :: l
      ! Result
      INTEGER:: ngrid_x

    END FUNCTION get_ngrid_x

    MODULE FUNCTION get_ngrid_y( THIS, l ) RESULT( ngrid_y )

      ! Arguments
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS
      INTEGER,           INTENT( IN )    :: l
      ! Result
      INTEGER:: ngrid_y

    END FUNCTION get_ngrid_y

    MODULE FUNCTION get_ngrid_z( THIS, l ) RESULT( ngrid_z )

      ! Arguments
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS
      INTEGER,           INTENT( IN )    :: l
      ! Result
      INTEGER:: ngrid_z

    END FUNCTION get_ngrid_z

    MODULE FUNCTION get_xR( THIS, l ) RESULT( xR )

      ! Arguments
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS
      INTEGER,           INTENT( IN )    :: l
      ! Result
      INTEGER:: xR

    END FUNCTION get_xR

    MODULE FUNCTION get_yR( THIS, l ) RESULT( yR )

      ! Arguments
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS
      INTEGER,           INTENT( IN )    :: l
      ! Result
      INTEGER:: yR

    END FUNCTION get_yR

    MODULE FUNCTION get_zR( THIS, l ) RESULT( zR )

      ! Arguments
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS
      INTEGER,           INTENT( IN )    :: l
      ! Result
      INTEGER:: zR

    END FUNCTION get_zR

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

    SUBROUTINE compute_and_export_3p1_variables_interface( THIS, namefile )

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
                                                             id, &
                                                             namefile, &
                                                             name_logfile )

      IMPORT:: formul_3p1
      IMPORT:: idbase
      CLASS(formul_3p1),   INTENT( IN OUT ):: THIS
      CLASS(idbase),          INTENT( IN OUT ):: id
      CHARACTER( LEN= * ), INTENT( IN OUT ):: namefile
      CHARACTER( LEN= * ), INTENT( IN OUT ):: name_logfile

    END SUBROUTINE compute_and_export_3p1_constraints_grid_interface

    SUBROUTINE compute_and_export_3p1_constraints_particles_interface( THIS, &
                                                             parts_obj, &
                                                             namefile, &
                                                             name_logfile )

      IMPORT:: formul_3p1
      IMPORT:: particles
      CLASS(formul_3p1),   INTENT( IN OUT ):: THIS
      CLASS(particles),    INTENT( IN OUT ):: parts_obj
      CHARACTER( LEN= * ), INTENT( IN OUT ):: namefile
      CHARACTER( LEN= * ), INTENT( IN OUT ):: name_logfile

    END SUBROUTINE compute_and_export_3p1_constraints_particles_interface

    SUBROUTINE deallocate_fields_interface( THIS )

      IMPORT:: formul_3p1
      CLASS(formul_3p1), INTENT( IN OUT ):: THIS

    END SUBROUTINE deallocate_fields_interface

  END INTERFACE


END MODULE formul_3p1_id
