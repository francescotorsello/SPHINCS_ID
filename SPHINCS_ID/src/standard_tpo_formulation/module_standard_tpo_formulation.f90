! File:         module_standard_tpo_formulation.f90
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

MODULE standard_tpo_formulation

  !***********************************************************************
  !
  !# This module contains the definition of ABSTRACT TYPE tpo
  !
  !***********************************************************************


  USE id_base,          ONLY: idbase
  USE sph_particles,    ONLY: particles
  USE mesh_refinement,  ONLY: grid_function_scalar, grid_function, level
  USE timing,           ONLY: timer
  USE utility,          ONLY: ios, err_msg, perc, creturn, run_id, &
                              test_status, show_progress


  IMPLICIT NONE


  !*******************************************************
  !                                                      *
  !     Definition of abstract TYPE tpo                  *
  !                                                      *
  ! Abstract TYPE for a \(3+1\) formulation of the       *
  ! Einstein equations. It imports the ID on the gravity *
  ! grid, in the standard 3+1 formulation, and defines   *
  ! DEFERRED PROCEDURES to be implemented in the derived *
  ! TYPES of the actual 3+1 formulations (for example,   *
  ! the BSSN formulation).                               *
  !                                                      *
  !*******************************************************

  TYPE, ABSTRACT:: tpo
  !! ABSTRACT TYPE representing a generic \(3+1\) formulation of the |ee|

    INTEGER:: n_matter
    !! Number of matter objects in the physical system

    INTEGER, DIMENSION(:), ALLOCATABLE:: npoints_xaxis
    !# Array containing the number of mesh points of the highest-resolution
    !  refinement level across the x-axis-size of the matter objects

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
    !  computed using the ID on the mesh
    TYPE(grid_function_scalar):: HC_parts
    !# Grid scalar function storing the Hamiltonian constraint (violations)
    !  computed using the stress-energy tensor mapped from the particles to the
    !  mesh
    TYPE(grid_function_scalar):: rho
    !# Grid scalar function storing the matter source in the Hamiltonian
    !  constraint computed using the ID on the mesh, multiplied by \(16\pi\)
    TYPE(grid_function):: S
    !# Grid function storing the matter source in the momentum
    !  constraint computed using the ID on the mesh, multiplied by \(8\pi\)
    TYPE(grid_function):: MC
    !# Grid function storing the momentum constraint (violations)
    !  computed using the ID on the mesh
    TYPE(grid_function):: MC_parts
    !# Grid function storing the momentum constraint (violations)
    !  computed using the stress-energy tensor mapped from the particles to the
    !  mesh
    TYPE(grid_function_scalar):: rho_parts
    !# Grid scalar function storing the matter source in the Hamiltonian
    !  constraint computed using the stress-energy tensor mapped from the
    !  particles to the mesh, multiplied by \(16\pi\)
    TYPE(grid_function):: S_parts
    !# Grid function storing the matter source in the momentum
    !  constraint computed using the stress-energy tensor mapped from the
    !  particles to the mesh, multiplied by \(8\pi\)

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
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: HC_int
    !# Integral the Hamiltonian constraint computed using the stress-energy
    !  tensor mapped from the particles
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: HC_parts_int
    !! Integral the Hamiltonian constraint computed on the mesh
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: MC_int
    !! Integral of the momentum constraint computed on the mesh
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: MC_parts_int
    !# Integral of the momentum constraint computed using the stress-energy
    !  tensor mapped from the particles

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


    !-----------------------------------!
    !--  NON_OVERRIDABLE SUBROUTINES  --!
    !-----------------------------------!


    PROCEDURE, NON_OVERRIDABLE:: setup_standard_tpo_variables
    !# Set up the refined mesh by reading the `gravity_grid_parameter.dat`
    !  parameter file, and read the standard 3+1 |id| using the given
    !  [[idbase]] object, on the refined mesh

    PROCEDURE, NON_OVERRIDABLE:: deallocate_standard_tpo_variables
    !! Deallocates memory for the standard 3+1 fields

    PROCEDURE, NON_OVERRIDABLE:: analyze_constraint
    !# Analyze a constraint (or an arbitrary scalar grid function) by
    !  examining its values at the refined mesh. Computes the number of mesh
    !  points at which the scalar grid function has values lying within
    !  predefined (hard-coded) intervals
    !  @todo complete this documentation entries with more details

    PROCEDURE, PUBLIC, NON_OVERRIDABLE:: test_recovery_m2p
    !# Computes the conserved variables from the physical ones, and vice versa,
    !  to test that the recovered physical variables are the same to those
    !  computed from the |id|. Uses the metric mapped from the grid to the
    !  particles. @todo add reference for recovery

    PROCEDURE, PUBLIC, NON_OVERRIDABLE:: compute_adm_momentum_fluid_m2p
    !# Computes an estimate of the \(\mathrm{ADM}\) linear momentum of the
    !  fluid using the |sph| hydro fields, and the spacetime metric mapped
    !  from the mesh
    !  @todo add reference

    PROCEDURE, NON_OVERRIDABLE:: print_summary
    !# Prints a summary about the features of the refined mesh


    !----------------------------!
    !--  DEFERRED SUBROUTINES  --!
    !----------------------------!


    PROCEDURE(define_allocate_fields_interface), DEFERRED:: &
                      define_allocate_fields
    !# Allocates memory for the fields specific to the formulation identified
    !  by an EXTENDED TYPE

    PROCEDURE(deallocate_fields_interface), DEFERRED:: deallocate_fields
    !# Deallocates memory for the fields specific to the formulation identified
    !  by an EXTENDED TYPE

    GENERIC, PUBLIC:: compute_and_print_tpo_constraints => &
                      compute_and_print_tpo_constraints_grid, &
                      compute_and_print_tpo_constraints_particles
    !# Overloaded PROCEDURE to compute the constraints using only the |id|
    !  on the refined mesh, or the spacetime |id| on the refined mesh and
    !  the hydrodynamical |id| mapped from the particles to the refined mesh

    PROCEDURE(compute_and_print_tpo_constraints_grid_interface), &
              DEFERRED:: compute_and_print_tpo_constraints_grid
    !# Computes the constraints specific to the formulation identified by an
    !  EXTENDED TYPE, using the full |id| on the refined mesh

    PROCEDURE(compute_and_print_tpo_constraints_particles_interface), &
              DEFERRED:: compute_and_print_tpo_constraints_particles
    !# Computes the constraints specific to the formulation identified by an
    !  EXTENDED TYPE, using the |bssn| |id| on the refined
    !  mesh and the hydrodynamical |id| mapped from the particles to the mesh

    PROCEDURE(compute_and_print_tpo_variables_interface), PUBLIC, &
                            DEFERRED:: compute_and_print_tpo_variables
    !# Compute the fields specific to the formulation identified by an
    !  EXTENDED TYPE, starting from the standard 3+1 fields

    PROCEDURE(print_formatted_id_tpo_variables_interface), PUBLIC, &
                            DEFERRED:: print_formatted_id_tpo_variables
    !! Prints the spacetime |id| to a formatted file


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


  END TYPE tpo

  !
  !-- Interface of the cores of the constructors and destructos of TYPES
  !-- derived from tpo
  !-- Their implementations are in submodule tpo_methods.f90
  !
  INTERFACE

    MODULE SUBROUTINE setup_standard_tpo_variables( ftpo, id, dx, dy, dz )

      CLASS(idbase),    INTENT( IN OUT ):: id
      CLASS(tpo), INTENT( IN OUT ):: ftpo
      DOUBLE PRECISION, OPTIONAL         :: dx, dy, dz

    END SUBROUTINE setup_standard_tpo_variables

 !   MODULE SUBROUTINE construct_tpo_bns_spacings( ftpo, id, &
 !                                                        dx, dy, dz )
 !
 !     CLASS(bns),        INTENT( IN OUT ):: id
 !     CLASS(tpo), INTENT( IN OUT ):: ftpo
 !     DOUBLE PRECISION,  INTENT( IN )    :: dx, dy, dz
 !
 !   END SUBROUTINE construct_tpo_bns_spacings

    MODULE SUBROUTINE deallocate_standard_tpo_variables( ftpo )

      CLASS(tpo), INTENT( IN OUT ):: ftpo

    END SUBROUTINE deallocate_standard_tpo_variables

  END INTERFACE

  !
  !-- Interface of the methods of TYPES derived from tpo
  !-- Their implementations are in submodule tpo_methods.f90
  !
  INTERFACE


    MODULE SUBROUTINE analyze_constraint( this, &
                                          l, &
                                          constraint, &
                                          name_constraint, &
                                          unit_logfile, &
                                          name_analysis, &
                                          l2_norm, &
                                          loo_norm, &
                                          integral, &
                                          source )

      CLASS(tpo),                  INTENT( IN OUT ):: this
      INTEGER,                            INTENT( IN ) :: l
      DOUBLE PRECISION, DIMENSION(:,:,:), INTENT( IN ) :: constraint
      CHARACTER( LEN= * ),                INTENT( IN ) :: name_constraint
      CHARACTER( LEN= * ),                INTENT( IN ) :: name_analysis
      INTEGER,                            INTENT( IN ) :: unit_logfile
      DOUBLE PRECISION,                   INTENT( OUT ):: l2_norm
      DOUBLE PRECISION,                   INTENT( OUT ):: loo_norm
      DOUBLE PRECISION,                   INTENT( OUT ):: integral
      DOUBLE PRECISION, DIMENSION(:,:,:), INTENT( IN ), OPTIONAL:: source

    END SUBROUTINE analyze_constraint


    MODULE SUBROUTINE test_recovery_m2p( this, parts, namefile )
    !# Tests the recovery. Computes the conserved variables from the physical
    !  ones, and then the physical ones from the conserved ones. It then
    !  compares the variables computed with the recovery PROCEDURES, with
    !  those computed with |sphincsid|. Uses the mesh-2-particle.
    !  @todo add reference for recovery

      CLASS(tpo),          INTENT(IN)   :: this
      !! [[tpo]] object which this PROCEDURE is a member of
      CLASS(particles),    INTENT(IN)   :: parts
      !# [[sph_particles:particles]] object used to map the metric
      !  from the mesh to the particles, and to call the recovery procedures
      CHARACTER( LEN= * ), INTENT(INOUT):: namefile
      !! Name of the formatted file where the data is printed

    END SUBROUTINE test_recovery_m2p


    MODULE SUBROUTINE compute_adm_momentum_fluid_m2p( this, parts, adm_mom )
    !# Computes an estimate of the \(\mathrm{ADM}\) linear momentum of the
    !  fluid using the |sph| hydro fields, and the spacetime metric mapped
    !  from the mesh
    !  @todo add reference

      CLASS(tpo),                     INTENT(IN) :: this
      !! [[tpo]] object which this PROCEDURE is a member of
      CLASS(particles),               INTENT(IN) :: parts
      !# [[sph_particles:particles]] object used to map the metric
      !  from the mesh to the particles, and to call the recovery procedures
      DOUBLE PRECISION, DIMENSION(3), INTENT(OUT):: adm_mom
      !# ADM linear momentum of the fluid computed using the metric mapped
      !  with the mesh-to-particle mapping

    END SUBROUTINE compute_adm_momentum_fluid_m2p


    MODULE SUBROUTINE print_summary( this )!, filename )
    !# Prints a summary of the properties of the refined mesh,
    !  and optionally, to a formatted file whose name
    !  is given as the optional argument `filename`


      CLASS(tpo), INTENT( IN OUT ):: this
      !CHARACTER( LEN= * ), INTENT( INOUT ), OPTIONAL:: filename
      !! Name of the formatted file to print the summary to

    END SUBROUTINE print_summary


    MODULE SUBROUTINE abs_values_in( this, lower_bound, upper_bound, &
                                     constraint, l, &
                                     export, unit_analysis, cnt )

      CLASS(tpo),                  INTENT( IN OUT ):: this
      DOUBLE PRECISION                                    :: lower_bound, &
                                                             upper_bound
      DOUBLE PRECISION, DIMENSION(:,:,:), INTENT( IN )    :: constraint
      INTEGER,                            INTENT( IN )    :: l
      LOGICAL,                            INTENT( IN )    :: export
      INTEGER,                            INTENT( IN )    :: unit_analysis
      INTEGER,                            INTENT( OUT )   :: cnt

    END SUBROUTINE abs_values_in


    MODULE FUNCTION get_grid_point( this, i, j, k, l ) RESULT( grid_point )

      ! Arguments
      CLASS(tpo), INTENT( IN OUT ):: this
      INTEGER,    INTENT( IN )    :: i, j, k, l
      ! Result
      DOUBLE PRECISION, DIMENSION(3)     :: grid_point

    END FUNCTION get_grid_point


    MODULE FUNCTION get_nlevels( this ) RESULT( nlevels )

      ! Arguments
      CLASS(tpo), INTENT( IN OUT ):: this
      ! Result
      INTEGER:: nlevels

    END FUNCTION get_nlevels


    MODULE FUNCTION get_levels( this ) RESULT( levels )

      ! Arguments
      CLASS(tpo), INTENT( IN OUT ):: this
      ! Result
      TYPE(level), DIMENSION(:), ALLOCATABLE:: levels

    END FUNCTION get_levels


    MODULE FUNCTION get_dx( this, l ) RESULT( dx )

      ! Arguments
      CLASS(tpo), INTENT( IN OUT ):: this
      INTEGER,    INTENT( IN )    :: l
      ! Result
      DOUBLE PRECISION:: dx

    END FUNCTION get_dx


    MODULE FUNCTION get_dy( this, l ) RESULT( dy )

      ! Arguments
      CLASS(tpo), INTENT( IN OUT ):: this
      INTEGER,    INTENT( IN )    :: l
      ! Result
      DOUBLE PRECISION:: dy

    END FUNCTION get_dy


    MODULE FUNCTION get_dz( this, l ) RESULT( dz )

      ! Arguments
      CLASS(tpo), INTENT( IN OUT ):: this
      INTEGER,    INTENT( IN )    :: l
      ! Result
      DOUBLE PRECISION:: dz

    END FUNCTION get_dz


    MODULE FUNCTION get_ngrid_x( this, l ) RESULT( ngrid_x )

      ! Arguments
      CLASS(tpo), INTENT( IN OUT ):: this
      INTEGER,    INTENT( IN )    :: l
      ! Result
      INTEGER:: ngrid_x

    END FUNCTION get_ngrid_x


    MODULE FUNCTION get_ngrid_y( this, l ) RESULT( ngrid_y )

      ! Arguments
      CLASS(tpo), INTENT( IN OUT ):: this
      INTEGER,    INTENT( IN )    :: l
      ! Result
      INTEGER:: ngrid_y

    END FUNCTION get_ngrid_y


    MODULE FUNCTION get_ngrid_z( this, l ) RESULT( ngrid_z )

      ! Arguments
      CLASS(tpo), INTENT( IN OUT ):: this
      INTEGER,    INTENT( IN )    :: l
      ! Result
      INTEGER:: ngrid_z

    END FUNCTION get_ngrid_z


    MODULE FUNCTION get_xR( this, l ) RESULT( xR )

      ! Arguments
      CLASS(tpo), INTENT( IN OUT ):: this
      INTEGER,    INTENT( IN )    :: l
      ! Result
      DOUBLE PRECISION:: xR

    END FUNCTION get_xR


    MODULE FUNCTION get_yR( this, l ) RESULT( yR )

      ! Arguments
      CLASS(tpo), INTENT( IN OUT ):: this
      INTEGER,    INTENT( IN )    :: l
      ! Result
      DOUBLE PRECISION:: yR

    END FUNCTION get_yR


    MODULE FUNCTION get_zR( this, l ) RESULT( zR )

      ! Arguments
      CLASS(tpo), INTENT( IN OUT ):: this
      INTEGER,    INTENT( IN )    :: l
      ! Result
      DOUBLE PRECISION:: zR

    END FUNCTION get_zR


    MODULE FUNCTION get_HC( this, i, j, k, l ) RESULT( HC_value )

      ! Arguments
      CLASS(tpo), INTENT( IN OUT ):: this
      INTEGER,    INTENT( IN )    :: i, j, k, l
      ! Result
      DOUBLE PRECISION                   :: HC_value

    END FUNCTION get_HC


    MODULE FUNCTION get_MC( this, i, j, k, l ) RESULT( MC_value )

      ! Arguments
      CLASS(tpo), INTENT( IN OUT ):: this
      INTEGER,    INTENT( IN )    :: i, j, k, l
      ! Result
      DOUBLE PRECISION, DIMENSION(3)     :: MC_value

    END FUNCTION get_MC


    MODULE FUNCTION get_HC_parts( this, i, j, k, l ) RESULT( HC_value )

      ! Arguments
      CLASS(tpo), INTENT( IN OUT ):: this
      INTEGER,    INTENT( IN )    :: i, j, k, l
      ! Result
      DOUBLE PRECISION                   :: HC_value

    END FUNCTION get_HC_parts


    MODULE FUNCTION get_MC_parts( this, i, j, k, l ) RESULT( MC_value )

      ! Arguments
      CLASS(tpo), INTENT( IN OUT ):: this
      INTEGER,    INTENT( IN )    :: i, j, k, l
      ! Result
      DOUBLE PRECISION, DIMENSION(3)     :: MC_value

    END FUNCTION get_MC_parts


  END INTERFACE

  !
  !-- Interfaces of the deferred methods of TYPE tpo
  !-- Their implementations are deferred to derived TYPES
  !
  ABSTRACT INTERFACE

    SUBROUTINE define_allocate_fields_interface( this )

      IMPORT:: tpo
      CLASS(tpo), INTENT( IN OUT ):: this

    END SUBROUTINE define_allocate_fields_interface

    SUBROUTINE compute_and_print_tpo_variables_interface( this, namefile )

      IMPORT:: tpo
      CLASS(tpo),   INTENT( IN OUT )           :: this
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE compute_and_print_tpo_variables_interface

    SUBROUTINE print_formatted_id_tpo_variables_interface &
                                                    ( this, namefile )

      IMPORT:: tpo
      CLASS(tpo),   INTENT( IN OUT )           :: this
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE print_formatted_id_tpo_variables_interface

    SUBROUTINE compute_and_print_tpo_constraints_grid_interface( this, &
                                                             id, &
                                                             namefile, &
                                                             name_logfile )

      IMPORT:: tpo
      IMPORT:: idbase
      CLASS(tpo),   INTENT( IN OUT ):: this
      CLASS(idbase),       INTENT( IN OUT ):: id
      CHARACTER( LEN= * ), INTENT( IN OUT ):: namefile
      CHARACTER( LEN= * ), INTENT( IN OUT ):: name_logfile

    END SUBROUTINE compute_and_print_tpo_constraints_grid_interface

    SUBROUTINE compute_and_print_tpo_constraints_particles_interface( this, &
                                                             parts_obj, &
                                                             namefile, &
                                                             name_logfile )

      IMPORT:: tpo
      IMPORT:: particles
      CLASS(tpo),   INTENT( IN OUT ):: this
      CLASS(particles),    INTENT( IN OUT ):: parts_obj
      CHARACTER( LEN= * ), INTENT( IN OUT ):: namefile
      CHARACTER( LEN= * ), INTENT( IN OUT ):: name_logfile

    END SUBROUTINE compute_and_print_tpo_constraints_particles_interface

    SUBROUTINE deallocate_fields_interface( this )

      IMPORT:: tpo
      CLASS(tpo), INTENT( IN OUT ):: this

    END SUBROUTINE deallocate_fields_interface

  END INTERFACE


END MODULE standard_tpo_formulation
