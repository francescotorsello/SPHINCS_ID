! File:         module_diffstar_lorene.f90
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

MODULE diffstar_lorene

  !***********************************************************
  !
  !#  This module contains the definition of TYPE diffstarlorene,
  !   and the SUBROUTINES that bind to the methods
  !   of |lorene|'s class |etrotdiff|, defined in
  !   Lorene/Export/C++/Source/Etoile
  !
  !   [|lorene| official repository](https://lorene.obspm.fr/index.html){:target="_blank"}
  !
  !***********************************************************


  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_CHAR, C_PTR
  USE diffstar_base,               ONLY: diffstarbase
  USE id_base,                     ONLY: idbase
  USE utility,                     ONLY: itr, ios, err_msg, test_status, &
                                         perc, creturn, compute_g4, &
                                         determinant_sym4x4, show_progress
  USE timing,                      ONLY: timer


  IMPLICIT NONE


  !*******************************************************
  !                                                      *
  !            Definition of TYPE diffstarlorene         *
  !                                                      *
  !   This class imports and stores the |lorene| DRS ID  *
  !                                                      *
  !*******************************************************

  TYPE, EXTENDS(diffstarbase):: diffstarlorene
  !! TYPE representing a differentially rotating star (DRS)


    PRIVATE


    INTEGER:: diffstar_identifier= 0
    !! Identifier of the diffstarlorene object
    INTEGER:: eos_loreneid
    !! |lorene| identifier for the EoS

    !& C pointer to the |lorene|'s Etdiffrot object
    ! N.B. This variable is global. The pointer to the second |lorene| Etdiffrot
    !      object will overwrite the first one, and so on.
    !      This variable stores the pointer to the last defined |lorene| Etdiffrot
    !      object. That's why it is not freed in the destructor of a bns object.
    !      Presently, it has to be freed by the user at the end of the PROGRAM.
    !      See the last part of the PROGRAM in setup_diffstar.f90, for example.
    TYPE(C_PTR):: diffstar_ptr


    CONTAINS


    !-------------------!
    !--  SUBROUTINES  --!
    !-------------------!

    PROCEDURE:: derived_type_constructor => construct_diffstarlorene

    PROCEDURE:: construct_drs
    !! Constructs the |lorene| Etdiffrot object

    PROCEDURE:: destruct_drs
    !! Destructs the |lorene| Etdiffrot object

    PROCEDURE:: allocate_diffstar_memory
    !! Allocates memory for the [[diffstarlorene]] member arrays

    PROCEDURE:: deallocate_diffstar_memory
    !! Deallocates memory for the [[diffstarlorene]] member arrays

    PROCEDURE:: import_diffstar_params
    !! Imports the parameters of the DRS from |lorene|

    PROCEDURE, PUBLIC:: print_diffstar_params
    !! Prints the parameters of the DRS to the standard output

    PROCEDURE:: import_id_int
    !! Stores the ID in the [[diffstarlorene]] member arrays

    PROCEDURE:: read_id_full      => import_id_full
    PROCEDURE:: read_id_spacetime => import_id_spacetime
    PROCEDURE:: read_id_particles => import_id_particles
    PROCEDURE:: read_id_hydro     => import_id_hydro
    PROCEDURE:: read_id_mass_b    => import_id_mass_b
    PROCEDURE:: read_id_k         => import_id_k

    PROCEDURE:: nothing
    !# Procedure that does nothing. It is used to instantiate a deferred
    !  idbase procedure which is not needed in TYPE [[diffstarlorene]].
    !  It also serves as a placeholder in case the idbase procedure
    !  will be needed in the future.

    PROCEDURE:: initialize_id => nothing

    !-----------------!
    !--  FUNCTIONS  --!
    !-----------------!

    PROCEDURE:: read_mass_density => import_mass_density
    !! Returns the |lorene|'s mass density at the given point

    PROCEDURE:: import_spatial_metric
    !! Returns the |lorene|'s conformally flat spatial ADM metric

    PROCEDURE:: test_position => is_hydro_positive
    !# Returns `.TRUE.` if the energy density or the specific energy or the
    !  pressure are positive

    !
    !-- Overloaded FUNCTION to access the fields as arrays and as values
    !

    GENERIC, PUBLIC:: get_field => get_fa, get_fv
    !# GENERIC PROCEDURE, overloded to access the [[diffstarlorene]]-member
    !  variables as arrays and as values
    PROCEDURE::       get_fa    => get_field_array
    !! Access the [[diffstarlorene]]-member arrays
    PROCEDURE::       get_fv    => get_field_value
    !! Access the components of the [[diffstarlorene]]-member arrays

    !
    !-- FUNCTIONS that access member variables
    !

    PROCEDURE:: get_eos_id => get_eos_loreneid
    !! Returns the |lorene| identifier for the EOS

    PROCEDURE:: return_eos_parameters => get_eos_parameters

    PROCEDURE, PUBLIC:: get_eos_loreneid
    !! Returns [[diffstarlorene:eos_loreneid]]

    PROCEDURE, PUBLIC:: get_diffstar_identifier
    !! Returns [[diffstarlorene:diffstar_identifier]]]

    !PROCEDURE, PUBLIC:: get_diffstar_ptr

    !PROCEDURE:: derived_type_destructor => destruct_diffstarlorene

    FINAL:: destruct_diffstarlorene
    !! Finalizer (Destructor) of a [[diffstarlorene]] object

    PROCEDURE, NOPASS:: finalize
    !# Corrects the |sph| |id| so that the linear \(\mathrm{ADM}\) momentum
    !  is zero

  END TYPE diffstarlorene

  !
  !-- Interface of the TYPE diffstarlorene
  !-- (i.e., declaration of the constructor)
  !
  !INTERFACE diffstarlorene
  !!! Interface of TYPE [[diffstarlorene]]
  !
  !  MODULE PROCEDURE:: construct_diffstarlorene
  !  !! Constructs a [[diffstarlorene]] object
  !
  !END INTERFACE diffstarlorene

  !
  !-- Interfaces of the constructor and destructor of the TYPE diffstarlorene
  !
  INTERFACE

   ! MODULE FUNCTION construct_diffstarlorene( &!derived_type,
   ! filename ) RESULT( foo )
   ! !! Constructs a [[diffstarlorene]] object
   !
   !   CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: filename
   !   !! |lorene| binary file containing the spectral DRS ID
   !   !CLASS(diffstarlorene):: derived_type
   !   CLASS(idbase), ALLOCATABLE:: foo
   !   !! Constructed [[diffstarlorene]] object
   !
   ! END FUNCTION construct_diffstarlorene

    MODULE SUBROUTINE destruct_diffstarlorene( this )
    !! Destruct a [[diffstarlorene]] object

      TYPE(diffstarlorene), INTENT( IN OUT ):: this
      !! [[diffstarlorene]] object to be destructed

    END SUBROUTINE destruct_diffstarlorene

  END INTERFACE

  !
  !-- Interfaces of the methods of the TYPE diffstarlorene
  !-- Their implementations are in submodule_diffstarlorene_methods.f90
  !
  INTERFACE


    !
    !-- SUBROUTINES
    !

    MODULE SUBROUTINE construct_diffstarlorene( derived_type, filename )
    !! Constructs a [[diffstarlorene]] object
    !# Prints a summary of the physical properties the system
    !  to the standard output and, optionally, to a formatted file whose name
    !  is given as the optional argument `filename`

      CHARACTER(LEN=*), INTENT( IN ), OPTIONAL :: filename
      !! |lorene| binary file containing the spectral DRS ID
      CLASS(diffstarlorene), INTENT( OUT ):: derived_type
      !! Constructed [[diffstarlorene]] object

    END SUBROUTINE construct_diffstarlorene


    MODULE SUBROUTINE construct_drs( this, resu_file )
    !! Interface of the subroutine that constructs the |lorene| Etdiffrot object

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),                     INTENT( IN OUT )      :: this
      !> |lorene| binary file containing the spectral DRS ID
      CHARACTER(KIND= C_CHAR, LEN=*), INTENT( IN ), OPTIONAL:: resu_file

    END SUBROUTINE construct_drs


    MODULE SUBROUTINE destruct_drs( this )
    !! Destructs a |lorene| Etdiffrot object

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene), INTENT( IN OUT ):: this

    END SUBROUTINE destruct_drs


    MODULE SUBROUTINE allocate_diffstar_memory( this, d )
    !! Allocates allocatable arrays member of a [[diffstarlorene]] object

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene), INTENT( IN OUT ):: this
      !> Dimension of the arrays
      INTEGER,    INTENT( IN )    :: d

    END SUBROUTINE allocate_diffstar_memory


    MODULE SUBROUTINE deallocate_diffstar_memory( this )
    !! Deallocates allocatable arrays member of a [[diffstarlorene]] object

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene), INTENT( IN OUT ):: this

    END SUBROUTINE deallocate_diffstar_memory


    MODULE SUBROUTINE import_diffstar_params( this )
    !! Imports the DRS parameters from |lorene|

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene), INTENT( IN OUT ):: this

    END SUBROUTINE import_diffstar_params


    MODULE SUBROUTINE print_diffstar_params( this )
    !! Prints the DRS parameters to the standard output

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene), INTENT( IN OUT ):: this

    END SUBROUTINE print_diffstar_params


    MODULE SUBROUTINE import_id_int( this, n, x, y, z )
    !! Stores the ID in the [[diffstarlorene]] member arrays

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),          INTENT( IN OUT ):: this
      INTEGER, INTENT( IN ):: n
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: x
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: y
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: z

    END SUBROUTINE import_id_int


    MODULE SUBROUTINE import_id_full( this, n, x, y, z,&
                                      lapse, &
                                      shift_x, shift_y, shift_z, &
                                      g_xx, g_xy, g_xz, &
                                      g_yy, g_yz, g_zz, &
                                      k_xx, k_xy, k_xz, &
                                      k_yy, k_yz, k_zz, &
                                      baryon_density, &
                                      energy_density, &
                                      specific_energy, &
                                      pressure, &
                                      u_euler_x, u_euler_y, u_euler_z )
    !# Stores the ID in non [[diffstarlorene]]-member arrays with the same
    !  shape as the [[diffstarlorene]] member arrays

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),          INTENT( IN OUT ):: this
      INTEGER,                        INTENT( IN )    :: n
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: x
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: y
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: z
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: lapse
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: shift_x
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: shift_y
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: shift_z
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_xx
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_xy
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_xz
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_yy
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_yz
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_zz
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_xx
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_xy
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_xz
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_yy
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_yz
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_zz
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: baryon_density
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: energy_density
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: specific_energy
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: pressure
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: u_euler_x
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: u_euler_y
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: u_euler_z

    END SUBROUTINE import_id_full


    MODULE SUBROUTINE import_id_spacetime( this, nx, ny, nz, &
                                              pos, &
                                              lapse, &
                                              shift, &
                                              g, &
                                              ek )
    !# Stores the spacetime ID in multi-dimensional arrays needed to compute
    !  the BSSN variables and constraints

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),                INTENT( IN OUT ):: this
      INTEGER,                              INTENT( IN )    :: nx
      INTEGER,                              INTENT( IN )    :: ny
      INTEGER,                              INTENT( IN )    :: nz
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN )    :: pos
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: lapse
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: shift
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: g
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: ek

    END SUBROUTINE import_id_spacetime


    MODULE SUBROUTINE import_id_hydro( this, nx, ny, nz, &
                                             pos, &
                                             baryon_density, &
                                             energy_density, &
                                             specific_energy, &
                                             pressure, &
                                             u_euler )
    !# Stores the hydro ID in the arrays needed to compute the constraints
    !  on the refined mesh

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),                INTENT( IN OUT ):: this
      INTEGER,                              INTENT( IN )    :: nx
      INTEGER,                              INTENT( IN )    :: ny
      INTEGER,                              INTENT( IN )    :: nz
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN )    :: pos
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: baryon_density
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: energy_density
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: specific_energy
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: pressure
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: u_euler

    END SUBROUTINE import_id_hydro


    MODULE SUBROUTINE import_id_particles( this, n, x, y, z, &
                                           lapse, &
                                           shift_x, shift_y, shift_z, &
                                           g_xx, g_xy, g_xz, &
                                           g_yy, g_yz, g_zz, &
                                           baryon_density, &
                                           energy_density, &
                                           specific_energy, &
                                           pressure, &
                                           u_euler_x, u_euler_y, u_euler_z )
    !! Stores the hydro ID in the arrays needed to compute the SPH ID

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),          INTENT( IN OUT ):: this
      INTEGER,                        INTENT( IN )    :: n
      REAL(C_DOUBLE),   DIMENSION(:), INTENT( IN )    :: x
      REAL(C_DOUBLE),   DIMENSION(:), INTENT( IN )    :: y
      REAL(C_DOUBLE),   DIMENSION(:), INTENT( IN )    :: z
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: lapse
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: shift_x
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: shift_y
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: shift_z
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_xx
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_xy
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_xz
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_yy
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_yz
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: g_zz
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: baryon_density
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: energy_density
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: specific_energy
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: pressure
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: u_euler_x
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: u_euler_y
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: u_euler_z

    END SUBROUTINE import_id_particles


    MODULE SUBROUTINE import_id_mass_b( this, x, y, z, &
                                        g, &
                                        baryon_density, &
                                        gamma_euler )
    !! Stores the hydro ID in the arrays needed to compute the baryon mass

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),          INTENT( IN OUT ):: this
      DOUBLE PRECISION,               INTENT( IN )    :: x
      DOUBLE PRECISION,               INTENT( IN )    :: y
      DOUBLE PRECISION,               INTENT( IN)     :: z
      DOUBLE PRECISION, DIMENSION(6), INTENT( OUT )   :: g
      DOUBLE PRECISION,               INTENT( OUT )   :: baryon_density
      DOUBLE PRECISION,               INTENT( OUT )   :: gamma_euler

    END SUBROUTINE import_id_mass_b


    MODULE SUBROUTINE import_id_k( this, n, x, y, z,&
                                         k_xx, k_xy, k_xz, &
                                         k_yy, k_yz, k_zz )
   !! Stores the components of the extrinsic curvature in arrays

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),          INTENT( IN OUT ):: this
      INTEGER,                        INTENT( IN )    :: n
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: x
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: y
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: z
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_xx
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_xy
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_xz
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_yy
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_yz
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: k_zz

    END SUBROUTINE import_id_k


    !
    !-- FUNCTIONS
    !
    MODULE FUNCTION import_mass_density( this, x, y, z ) RESULT( res )
    !! Returns the |lorene| baryon mass density at a point \((x,y,z)\)

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene), INTENT( IN ):: this
      !> \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: x
      !> \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: y
      !> \(z\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: z
      !> Baryon mass density at \((x,y,z)\)
      DOUBLE PRECISION:: res

    END FUNCTION import_mass_density


    MODULE FUNCTION import_spatial_metric( this, x, y, z ) RESULT( res )
    !# Returns the |lorene| conformally flat spatial metric component
    !  \(g_{xx}=g_{yy}=g_{zz}\) at a point \((x,y,z)\)

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),     INTENT( IN )       :: this
      !> \(x\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT( IN ), VALUE:: x
      !> \(y\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT( IN ), VALUE:: y
      !> \(z\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT( IN ), VALUE:: z
      !> \(g_{xx}=g_{yy}=g_{zz}\) at \((x,y,z)\)
      REAL(C_DOUBLE):: res

    END FUNCTION import_spatial_metric


    MODULE FUNCTION is_hydro_positive( this, x, y, z ) RESULT( res )
    !# Returns .TRUE. if the energy density or the specific energy or the
    !  pressure are positive, .FALSE. otherwise

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),     INTENT( IN )       :: this
      !> \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: x
      !> \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: y
      !> \(z\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: z
      !& 1 if the energy density or the specific energy or the pressure
      !  are negative, 0 otherwise
      LOGICAL:: res

    END FUNCTION is_hydro_positive


    MODULE FUNCTION get_field_array( this, field ) RESULT( field_array )
    !! Returns the [[diffstarlorene]] member arrays named field

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),          INTENT( IN )             :: this
      !> Name of the desired [[diffstarlorene]] member array
      CHARACTER( LEN= : ), INTENT( IN ), ALLOCATABLE:: field
      !> Desired [[diffstarlorene]] member array
      DOUBLE PRECISION, DIMENSION(:),    ALLOCATABLE:: field_array

    END FUNCTION get_field_array


    MODULE FUNCTION get_field_value( this, field, n ) RESULT( field_value )
    !! Returns the component n of the [[diffstarlorene]] member arrays named field

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),          INTENT( IN )             :: this
      !> Name of the desired [[diffstarlorene]] member array
      CHARACTER( LEN= : ), INTENT( IN ), ALLOCATABLE:: field
      !> Component of the desired [[diffstarlorene]] member array
      INTEGER,             INTENT( IN )             :: n
      !> Component n of the desired [[diffstarlorene]] member array
      DOUBLE PRECISION                              :: field_value

    END FUNCTION get_field_value


    MODULE FUNCTION get_diffstar_identifier( this )

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene), INTENT( IN ):: this
      ! Result
      DOUBLE PRECISION:: get_diffstar_identifier

    END FUNCTION get_diffstar_identifier


    MODULE FUNCTION get_eos_loreneid( this )

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene), INTENT( IN ):: this
      ! Result
      INTEGER:: get_eos_loreneid

    END FUNCTION get_eos_loreneid


    MODULE SUBROUTINE get_eos_parameters( this, i_matter, eos_params )

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene), INTENT( IN ):: this
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT):: eos_params
      !# Array containing the parameters of the |eos| for the DRS

    END SUBROUTINE get_eos_parameters


    MODULE SUBROUTINE finalize &
      ( npart, pos, nlrf, u, pr, vel_u, theta, nstar, nu )
    !# Post-process the |sph| |id|; for example, correct for the residual
    !  ADM linear momentum.

      !IMPORT:: idbase
      !CLASS(idbase),                        INTENT(IN)   :: this
      INTEGER,                              INTENT(IN)   :: npart
      !! Particle number
      DOUBLE PRECISION, DIMENSION(3,npart), INTENT(INOUT):: pos
      !! Particle positions
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(INOUT):: nlrf
      !! Baryon density in the local rest frame on the particles
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(INOUT):: u
      !! Specific internal energy on the particles
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(INOUT):: pr
      !! Pressure on the particles
      DOUBLE PRECISION, DIMENSION(3,npart), INTENT(INOUT):: vel_u
      !! Spatial velocity in the computing frame on the particles
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(INOUT):: theta
      !! Generalized Lorentz factor on the particles
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(INOUT):: nstar
      !! Proper baryon density in the local rest frame on the particles
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(INOUT):: nu
      !! Baryon number per particle

    END SUBROUTINE finalize


    MODULE SUBROUTINE nothing( this, flag, switch )
    !# Procedure that does nothing. It is used to instantiate a deferred
    !  idbase procedure which is not needed in TYPE [[diffstarlorene]].
    !  It also serves as a placeholder in case the idbase procedure
    !  will be needed in the future.

      CLASS(diffstarlorene), INTENT(INOUT)       :: this
      INTEGER,               INTENT(IN)          :: flag
      !! Identifies what kind of initialization has to be done
      LOGICAL,               INTENT(IN), OPTIONAL:: switch
      !! If `.TRUE.`, switch to a different initialization

    END SUBROUTINE nothing


    !MODULE FUNCTION get_diffstar_ptr( this )
    !
    !  ! Argument
    !  CLASS(diffstarlorene), INTENT( IN ):: this
    !  ! Result
    !  TYPE(C_PTR):: get_diffstar_ptr
    !
    !END FUNCTION get_diffstar_ptr


  END INTERFACE


  PRIVATE:: construct_etdiffrot, get_diffstar_full, get_diffstar_spacetime, &
            get_diffstar_particles, get_diffstar_mass_b, &
            get_diffstar_hydro, get_diffstar_mass_density, &
            get_diffstar_spatial_metric, positive_hydro, get_diffstar_params, &
            destruct_etdiffrot


  !-----------------------------------------------------------------!
  !--  Interfaces to the methods of |lorene|'s class |etdiffrot|  --!
  !-----------------------------------------------------------------!


  INTERFACE


    FUNCTION construct_etdiffrot( c_resu_file ) RESULT( optr ) &
      BIND(C, NAME= "construct_et_diffrot")

      !***********************************************
      !
      !# Interface to the |lorene| method of class
      !  |etdiffrot| with the same name, that constructs
      !  the |lorene| |etdiffrot| object
      !
      !  FT 24.10.2021
      !
      !***********************************************

      IMPORT :: C_PTR, C_CHAR

      IMPLICIT NONE

      !& C string of the name of the |lorene| binary file storing the spectral
      !  DRS ID
      CHARACTER(KIND= C_CHAR), DIMENSION(*), INTENT(IN), OPTIONAL :: &
                                                              c_resu_file
      !> C pointer pointing to the constructed |lorene| |etdiffrot| object
      TYPE(C_PTR) :: optr

    END FUNCTION construct_etdiffrot


    SUBROUTINE get_diffstar_full( optr, &
                                  x, y, z, &
                                  lapse, &
                                  shift_x, shift_y, shift_z, &
                                  g_rr, g_tt, g_pp, &
                                  k_xx, k_xy, k_xz, &
                                  k_yy, k_yz, k_zz, &
                                  baryon_density, &
                                  energy_density, &
                                  specific_energy, &
                                  pressure, &
                                  v_euler_x, v_euler_y, v_euler_z ) &
      BIND(C, NAME= "get_rotdiff_id")

      !*************************************************
      !
      !# Interface to the |lorene| method of class
      !  |etdiffrot| with the same name, that reads the full
      !  |lorene| ID at the specified point.
      !  That is, imports the metric fields, the
      !  components of the extrinsic curvature [c/km],
      !  and the hydro fields.
      !
      !  - shift vector [c]
      !  - baryon mass density [kg m^{-3}]
      !  - energy density [kg c^2 m^{-3}]
      !  - pressure [kg c^2 m^{-3}]
      !  - specific internal energy [c^2]
      !  - fluid 3-velocity with respect to the
      !    Eulerian observer [c]
      !
      !  FT 24.10.2021
      !
      !*************************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      !> C pointer pointing to a |lorene| |etdiffrot| object
      TYPE(C_PTR),    INTENT(IN), VALUE :: optr
      !> \(x\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: x
      !> \(y\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: y
      !> \(z\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: z
      REAL(C_DOUBLE), INTENT(OUT)       :: lapse
      REAL(C_DOUBLE), INTENT(OUT)       :: shift_x
      REAL(C_DOUBLE), INTENT(OUT)       :: shift_y
      REAL(C_DOUBLE), INTENT(OUT)       :: shift_z
      REAL(C_DOUBLE), INTENT(OUT)       :: g_rr, g_tt, g_pp
      REAL(C_DOUBLE), INTENT(OUT)       :: k_xx
      REAL(C_DOUBLE), INTENT(OUT)       :: k_xy
      REAL(C_DOUBLE), INTENT(OUT)       :: k_xz
      REAL(C_DOUBLE), INTENT(OUT)       :: k_yy
      REAL(C_DOUBLE), INTENT(OUT)       :: k_yz
      REAL(C_DOUBLE), INTENT(OUT)       :: k_zz
      REAL(C_DOUBLE), INTENT(OUT)       :: baryon_density
      REAL(C_DOUBLE), INTENT(OUT)       :: energy_density
      REAL(C_DOUBLE), INTENT(OUT)       :: specific_energy
      REAL(C_DOUBLE), INTENT(OUT)       :: pressure
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_x
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_y
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_z

    END SUBROUTINE get_diffstar_full


    SUBROUTINE get_diffstar_spacetime( optr, &
                                       x, y, z, &
                                       lapse, &
                                       shift_x, shift_y, shift_z, &
                                       g_rr, g_tt, g_pp, &
                                       k_xx, k_xy, k_xz, &
                                       k_yy, k_yz, k_zz ) &
      BIND(C, NAME= "get_rotdiff_spacetime")

      !*************************************************
      !
      !# Interface to the |lorene| method of class
      !  |etdiffrot| with the same name, that reads the
      !  metric fields and the components
      !  of the extrinsic curvature [c/km] from |lorene|,
      !  at the specified point
      !
      !  FT 24.10.2021
      !
      !*************************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      !> C pointer pointing to a |lorene| |etdiffrot| object
      TYPE(C_PTR),    INTENT(IN), VALUE :: optr
      !> \(x\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: x
      !> \(y\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: y
      !> \(z\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: z
      REAL(C_DOUBLE), INTENT(OUT)       :: lapse
      REAL(C_DOUBLE), INTENT(OUT)       :: shift_x
      REAL(C_DOUBLE), INTENT(OUT)       :: shift_y
      REAL(C_DOUBLE), INTENT(OUT)       :: shift_z
      REAL(C_DOUBLE), INTENT(OUT)       :: g_rr, g_tt, g_pp
      REAL(C_DOUBLE), INTENT(OUT)       :: k_xx
      REAL(C_DOUBLE), INTENT(OUT)       :: k_xy
      REAL(C_DOUBLE), INTENT(OUT)       :: k_xz
      REAL(C_DOUBLE), INTENT(OUT)       :: k_yy
      REAL(C_DOUBLE), INTENT(OUT)       :: k_yz
      REAL(C_DOUBLE), INTENT(OUT)       :: k_zz

    END SUBROUTINE get_diffstar_spacetime


    SUBROUTINE get_diffstar_particles( optr, &
                                       x, y, z, &
                                       lapse, &
                                       shift_x, shift_y, shift_z, &
                                       g_rr, g_tt, g_pp, &
                                       baryon_density, &
                                       energy_density, &
                                       specific_energy, &
                                       pressure, &
                                       v_euler_x, v_euler_y, v_euler_z ) &
      BIND(C, NAME= "get_rotdiff_particles")

      !**********************************************
      !
      !# Interface to the |lorene| method of class
      !  |etdiffrot| with the same name, that reads the
      !  hydro fields and the metric fields *
      !  from |lorene|, at the specified point
      !
      !  - shift vector [c]
      !  - baryon mass density [kg m^{-3}]
      !  - energy density [kg c^2 m^{-3}]
      !  - pressure [kg c^2 m^{-3}]
      !  - specific internal energy [c^2]
      !  - fluid 3-velocity with respect to the
      !    Eulerian observer [c]
      !
      !  FT 24.10.2021
      !
      !**********************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      !> C pointer pointing to a |lorene| |etdiffrot| object
      TYPE(C_PTR),    INTENT(IN), VALUE :: optr
      !> \(x\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: x
      !> \(y\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: y
      !> \(z\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: z
      REAL(C_DOUBLE), INTENT(OUT)       :: lapse
      REAL(C_DOUBLE), INTENT(OUT)       :: shift_x
      REAL(C_DOUBLE), INTENT(OUT)       :: shift_y
      REAL(C_DOUBLE), INTENT(OUT)       :: shift_z
      REAL(C_DOUBLE), INTENT(OUT)       :: g_rr, g_tt, g_pp
      REAL(C_DOUBLE), INTENT(OUT)       :: baryon_density
      REAL(C_DOUBLE), INTENT(OUT)       :: energy_density
      REAL(C_DOUBLE), INTENT(OUT)       :: specific_energy
      REAL(C_DOUBLE), INTENT(OUT)       :: pressure
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_x
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_y
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_z

    END SUBROUTINE get_diffstar_particles


    SUBROUTINE get_diffstar_mass_b( optr, &
                                    x, y, z, &
                                    g_rr, g_tt, g_pp, &
                                    baryon_density, &
                                    gamma_euler ) &
      BIND(C, NAME= "get_rotdiff_mass_b")

      !************************************************
      !
      !# Interface to the |lorene| method of class
      !  |etdiffrot| with the same name, that reads the
      !  hydro fields and the metric fields
      !  from |lorene|, at the specified point,
      !  needed to compute the baryon mass.
      !
      !  - shift vector [c]
      !  - baryon mass density [kg m^{-3}]
      !  - fluid 3-velocity with respect to the
      !    Eulerian observer [c]
      !
      !  FT 24.10.2021
      !
      !************************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      !> C pointer pointing to a |lorene| |etdiffrot| object
      TYPE(C_PTR),    INTENT(IN), VALUE :: optr
      !> \(x\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: x
      !> \(y\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: y
      !> \(z\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: z
      !> \(g_{xx}=g_{yy}=g_{zz}\) at \(x,y,z\)
      REAL(C_DOUBLE), INTENT(OUT)       :: g_rr, g_tt, g_pp
      !> Baryon mass density at \(x,y,z\)
      REAL(C_DOUBLE), INTENT(OUT)       :: baryon_density
      !& Relative Lorentz factor between the 4-velocity of the fluid
      !  wrt the Eulerian observer and the 4-velocity of the Eulerian observer
      !  at \(x,y,z\)
      REAL(C_DOUBLE), INTENT(OUT)       :: gamma_euler

    END SUBROUTINE get_diffstar_mass_b


    SUBROUTINE get_diffstar_hydro( optr, &
                                   x, y, z, &
                                   baryon_density, &
                                   energy_density, &
                                   specific_energy, &
                                   pressure, &
                                   v_euler_x, v_euler_y, v_euler_z ) &
      BIND(C, NAME= "get_rotdiff_hydro")

      !***********************************************
      !
      !# Interface to the |lorene| method of class
      !  |etdiffrot| with the same name, that reads the
      !  hydro fields from |lorene|, at the
      !  specified point
      !
      !  - baryon mass density [kg m^{-3}]
      !  - energy density [kg c^2 m^{-3}]
      !  - pressure [kg c^2 m^{-3}]
      !  - specific internal energy [c^2]
      !  - fluid 3-velocity with respect to the
      !    Eulerian observer [c]
      !
      !  FT 24.10.2021
      !
      !***********************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      !> C pointer pointing to a |lorene| |etdiffrot| object
      TYPE(C_PTR),    INTENT(IN), VALUE :: optr
      !> \(x\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: x
      !> \(y\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: y
      !> \(z\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: z
      REAL(C_DOUBLE), INTENT(OUT)       :: baryon_density
      REAL(C_DOUBLE), INTENT(OUT)       :: energy_density
      REAL(C_DOUBLE), INTENT(OUT)       :: specific_energy
      REAL(C_DOUBLE), INTENT(OUT)       :: pressure
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_x
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_y
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_z

    END SUBROUTINE get_diffstar_hydro


    FUNCTION get_diffstar_mass_density( optr, x, y, z ) RESULT( res ) &
      BIND(C, NAME= "get_rotdiff_mass_density")

      !********************************************
      !
      !# Interface to the |lorene| method of class
      !  |etdiffrot| with the same name, that returns
      !  the baryon mass density \([\mathrm{kg}\,
      !  \mathrm{m}^{-3}]\) from |lorene|,
      !  at the specified point
      !
      !  FT 24.10.2021
      !
      !********************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      !> C pointer pointing to a |lorene| |etdiffrot| object
      TYPE(C_PTR),    INTENT(IN),  VALUE :: optr
      !> \(x\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: x
      !> \(y\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: y
      !> \(z\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: z
      !& Baryon mass density \([\mathrm{kg}\, \mathrm{m}^{-3}]\) at the desired
      !  point \((x,y,z)\)
      REAL(C_DOUBLE) :: res

    END FUNCTION get_diffstar_mass_density


    FUNCTION get_diffstar_spatial_metric( optr, x, y, z ) RESULT( res ) &
      BIND(C, NAME= "get_g_rr")

      !************************************************
      !
      !# Interface to the |lorene| method of class
      !  |etdiffrot| with the same name, that returns the
      !  diagonal components of the metric,
      !  all equal to the |lorene| conformal factor to
      !  the 4th power.
      !
      !  FT 24.10.2021
      !
      !************************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      !> C pointer pointing to a |lorene| |etdiffrot| object
      TYPE(C_PTR),    INTENT(IN),  VALUE :: optr
      !> \(x\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: x
      !> \(y\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: y
      !> \(z\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: z
      !& Spatial metric component
      !  \(g_{xx}=g_{yy}=g_{zz}\) at the point \((x,y,z)\)
      REAL(C_DOUBLE) :: res

    END FUNCTION get_diffstar_spatial_metric


    FUNCTION positive_hydro( optr, x, y, z ) RESULT( res ) &
      BIND(C, NAME= "is_hydro_positive_diffrot")

      !************************************************
      !
      !# Interface to the |lorene| method of class
      !  |etdiffrot| with the same name, that returns 1
      !  if the energy density is positive,
      !  and if the specific energy is positive,
      !  and if the pressure is nonpositive,
      !  at the specified point; it returns 0 otherwise
      !
      !  FT 24.10.2021
      !
      !************************************************

      IMPORT :: C_INT, C_DOUBLE, C_PTR

      IMPLICIT NONE

      !> C pointer pointing to a |lorene| |etdiffrot| object
      TYPE(C_PTR),    INTENT(IN),  VALUE :: optr
      !> \(x\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: x
      !> \(y\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: y
      !> \(z\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: z
      !& 1 if the energy density or the specific energy or the pressure
      !  are positive, 0 otherwise
      INTEGER(C_INT) :: res

    END FUNCTION positive_hydro


    SUBROUTINE get_diffstar_params( optr,                           &
                                    omega_c,                        &
                                    mass,                           &
                                    mass_grav,                      &
                                    angular_momentum,               &
                                    tsw,                            &
                                    grv2,                           &
                                    grv3,                           &
                                    r_circ,                         &
                                    surface_area,                   &
                                    r_mean,                         &
                                    r_eq,                           &
                                    r_eq_pi2,                       &
                                    r_eq_pi,                        &
                                    r_eq_3pi2,                      &
                                    r_eq_pole,                      &
                                    r_ratio,                        &
                                    r_isco,                         &
                                    f_isco,                         &
                                    specific_energy_isco,           &
                                    specific_angular_momentum_isco, &
                                    area_radius,                    &
                                    ent_center,                     &
                                    nbar_center,                    &
                                    rho_center,                     &
                                    energy_density_center,          &
                                    specific_energy_center,         &
                                    pressure_center,                &
                                    redshift_eqf,                   &
                                    redshift_eqb,                   &
                                    redshift_pole,                  &
                                    eos,                            &
                                    eos_id,                         &
                                    gamma,                          &
                                    kappa,                          &
                                    npeos,                          &
                                    gamma0,                         &
                                    gamma1,                         &
                                    gamma2,                         &
                                    gamma3,                         &
                                    kappa0,                         &
                                    kappa1,                         &
                                    kappa2,                         &
                                    kappa3,                         &
                                    logP1,                          &
                                    logRho0,                        &
                                    logRho1,                        &
                                    logRho2 )                      &
      BIND(C, NAME= "get_rotdiff_params")

      !**********************************************
      !
      !# Interface to the |lorene| method of class
      !  |etdiffrot| with the same name, that stores
      !  the physical parameters of the binary
      !  system from |lorene| in the desired variables
      !
      !  FT 24.10.2021
      !
      !**********************************************

      IMPORT :: C_INT, C_DOUBLE, C_PTR, C_CHAR

      IMPLICIT NONE

      !> C pointer pointing to a |lorene| |etdiffrot| object
      TYPE(C_PTR),    INTENT(IN), VALUE :: optr
      REAL(C_DOUBLE), INTENT(OUT)       :: omega_c
      REAL(C_DOUBLE), INTENT(OUT)       :: mass
      REAL(C_DOUBLE), INTENT(OUT)       :: mass_grav
      REAL(C_DOUBLE), INTENT(OUT)       :: angular_momentum
      REAL(C_DOUBLE), INTENT(OUT)       :: tsw
      REAL(C_DOUBLE), INTENT(OUT)       :: grv2
      REAL(C_DOUBLE), INTENT(OUT)       :: grv3
      REAL(C_DOUBLE), INTENT(OUT)       :: r_circ
      REAL(C_DOUBLE), INTENT(OUT)       :: surface_area
      REAL(C_DOUBLE), INTENT(OUT)       :: r_mean
      REAL(C_DOUBLE), INTENT(OUT)       :: r_eq
      REAL(C_DOUBLE), INTENT(OUT)       :: r_eq_pi2
      REAL(C_DOUBLE), INTENT(OUT)       :: r_eq_pi
      REAL(C_DOUBLE), INTENT(OUT)       :: r_eq_3pi2
      REAL(C_DOUBLE), INTENT(OUT)       :: r_eq_pole
      REAL(C_DOUBLE), INTENT(OUT)       :: r_ratio
      REAL(C_DOUBLE), INTENT(OUT)       :: r_isco
      REAL(C_DOUBLE), INTENT(OUT)       :: f_isco
      REAL(C_DOUBLE), INTENT(OUT)       :: specific_energy_isco
      REAL(C_DOUBLE), INTENT(OUT)       :: specific_angular_momentum_isco
      REAL(C_DOUBLE), INTENT(OUT)       :: area_radius
      REAL(C_DOUBLE), INTENT(OUT)       :: ent_center
      REAL(C_DOUBLE), INTENT(OUT)       :: nbar_center
      REAL(C_DOUBLE), INTENT(OUT)       :: rho_center
      REAL(C_DOUBLE), INTENT(OUT)       :: energy_density_center
      REAL(C_DOUBLE), INTENT(OUT)       :: specific_energy_center
      REAL(C_DOUBLE), INTENT(OUT)       :: pressure_center
      REAL(C_DOUBLE), INTENT(OUT)       :: redshift_eqf
      REAL(C_DOUBLE), INTENT(OUT)       :: redshift_eqb
      REAL(C_DOUBLE), INTENT(OUT)       :: redshift_pole
      CHARACTER(KIND=C_CHAR), DIMENSION(100), INTENT(OUT):: eos
      INTEGER(C_INT)                    :: eos_id
      REAL(C_DOUBLE), INTENT(OUT)       :: gamma
      REAL(C_DOUBLE), INTENT(OUT)       :: kappa
      INTEGER(C_INT)                    :: npeos
      REAL(C_DOUBLE), INTENT(OUT)       :: gamma0
      REAL(C_DOUBLE), INTENT(OUT)       :: gamma1
      REAL(C_DOUBLE), INTENT(OUT)       :: gamma2
      REAL(C_DOUBLE), INTENT(OUT)       :: gamma3
      REAL(C_DOUBLE), INTENT(OUT)       :: kappa0
      REAL(C_DOUBLE), INTENT(OUT)       :: kappa1
      REAL(C_DOUBLE), INTENT(OUT)       :: kappa2
      REAL(C_DOUBLE), INTENT(OUT)       :: kappa3
      REAL(C_DOUBLE), INTENT(OUT)       :: logP1
      REAL(C_DOUBLE), INTENT(OUT)       :: logRho0
      REAL(C_DOUBLE), INTENT(OUT)       :: logRho1
      REAL(C_DOUBLE), INTENT(OUT)       :: logRho2

    END SUBROUTINE get_diffstar_params


    SUBROUTINE destruct_etdiffrot( optr ) &
      BIND(C, NAME= "destruct_et_diffrot")

      !**********************************************
      !
      !# Interface to the |lorene| method of class
      !  |etdiffrot| with the same name, that destructs
      !  the |lorene| |etdiffrot| object
      !
      ! FT 24.10.2021
      !
      !**********************************************

      IMPORT :: C_PTR

      IMPLICIT NONE

      TYPE(C_PTR), INTENT(IN), VALUE :: optr
      !! C pointer pointing to the |lorene| |etdiffrot| object to destruct

    END SUBROUTINE destruct_etdiffrot


  END INTERFACE


END MODULE diffstar_lorene
