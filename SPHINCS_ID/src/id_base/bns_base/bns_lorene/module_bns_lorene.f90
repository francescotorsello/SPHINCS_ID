! File:         module_bns_lorene.f90
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

MODULE bns_lorene

  !***********************************************************
  !
  !#  This module contains the definition of TYPE [[bnslorene]],
  !   and the SUBROUTINES that bind to the methods
  !   of |lorene|'s class |binns|, defined in
  !   \(\mathrm{\$HOME\_LORENE/Export/BinNS}\)
  !
  !   [|lorene| official repository](https://lorene.obspm.fr/index.html){:target="_blank"}
  !
  !***********************************************************


  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_CHAR, C_NULL_CHAR, &
                                         C_PTR, C_NULL_PTR, C_ASSOCIATED
  USE bns_base,                    ONLY: bnsbase
  USE id_base,                     ONLY: idbase
  USE utility,                     ONLY: itr, ios, err_msg, &
                                         perc, creturn, compute_g4, &
                                         determinant_sym4x4, show_progress
  USE timing,                      ONLY: timer


  IMPLICIT NONE


  !*******************************************************
  !                                                      *
  !            Definition of TYPE bnslorene              *
  !                                                      *
  !   This class reads and stores the |lorene| BNS ID    *
  !                                                      *
  !*******************************************************

  TYPE, EXTENDS(bnsbase):: bnslorene
  !# TYPE representing a binary system of neutron stars (|bns|) produced with
  !  |lorene|


    PRIVATE


    INTEGER:: bns_identifier= 0
    !! Identifier of the bnslorene object
    INTEGER:: eos1_loreneid
    !! |lorene| identifier for the |eos| of star 1
    INTEGER:: eos2_loreneid
    !! |lorene| identifier for the |eos| of star 2

    !
    !-- Hydro fields
    !

    !> 1-D array storing the baryon mass density in the fluid frame [kg m^{-3}]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: baryon_density
    !> 1-D array storing the energy density [kg c^2 m^{-3}]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: energy_density
    !> 1-D array storing the specific internal energy [c^2]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: specific_energy
    !> 1-D array storing the pressure
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: pressure
    !> 1-D array storing the x component of the fluid 3-velocity with respect to
    !  the Eulerian observer [c]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: v_euler_x
    !> 1-D array storing the y component of the fluid 3-velocity with respect to
    !  the Eulerian observer [c]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: v_euler_y
    !> 1-D array storing the z component of the fluid 3-velocity with respect to
    !  the Eulerian observer [c]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: v_euler_z

    !& C pointer to the |lorene|'s |binns| object
    ! N.B. This variable is global. The pointer to the second |lorene| |binns|
    !      object will overwrite the first one, and so on.
    !      This variable stores the pointer to the last defined |lorene| |binns|
    !      object. That's why it is not freed in the destructor of a bns object.
    !      Presently, it has to be freed by the user at the end of the PROGRAM.
    !      See the last part of the PROGRAM in setup_lorene_id.f90, for example.
    TYPE(C_PTR):: bns_ptr


    CONTAINS


    !-------------------!
    !--  SUBROUTINES  --!
    !-------------------!

    PROCEDURE:: derived_type_constructor => construct_bnslorene

    PROCEDURE:: construct_binary
    !! Constructs the |lorene| |binns| object

    PROCEDURE:: destruct_binary
    !! Destructs the |lorene| |binns| object

    PROCEDURE:: allocate_bnslorene_memory
    !! Allocates memory for the [[bnslorene]] member arrays

    PROCEDURE:: deallocate_bnslorene_memory
    !! Deallocates memory for the [[bnslorene]] member arrays

    PROCEDURE:: read_id_params
    !! Imports the parameters of the BNS from |lorene|

    !PROCEDURE:: integrate_field_on_star => integrate_baryon_mass_density
    !# Integrates the |lorene| baryon mass density and computes the
    !  radial mass profile

    PROCEDURE, PUBLIC:: print_id_params
    !! Prints the parameters of the BNS to the standard output

    PROCEDURE:: read_id_int
    !! Stores the ID in the [[bnslorene]] member arrays

    PROCEDURE:: read_id_full      => read_id_full
    PROCEDURE:: read_id_spacetime => read_id_spacetime
    PROCEDURE:: read_id_particles => read_id_particles
    PROCEDURE:: read_id_hydro     => read_id_hydro
    PROCEDURE:: read_id_mass_b    => read_id_mass_b
    PROCEDURE:: read_id_k         => read_id_k

    PROCEDURE:: print_summary_derived => print_summary_bnslorene

    PROCEDURE:: nothing
    !# Procedure that does nothing. It is used to instantiate a deferred
    !  idbase procedure which is not needed in TYPE [[bnslorene]].
    !  It also serves as a placeholder in case the idbase procedure
    !  will be needed in the future.

    PROCEDURE:: initialize_id => nothing

    PROCEDURE, NOPASS:: correct_adm_linear_momentum
    !# Corrects the |sph| |id| so that the linear \(\mathrm{ADM}\) momentum
    !  is zero

    !-----------------!
    !--  FUNCTIONS  --!
    !-----------------!

    !> Returns the |lorene|'s mass density at the desired point
    PROCEDURE:: read_mass_density => read_bnslorene_mass_density

    !> Returns the |lorene|'s pressure at the desired point
    PROCEDURE:: read_pressure => read_bnslorene_pressure

    !> Returns the |lorene|'s conformally flat spatial ADM metric
    PROCEDURE:: read_spatial_metric

    !& Returns 1 if the energy density or the specific energy or the pressure
    !  are negative
    PROCEDURE:: test_position => is_hydro_positive

    !PROCEDURE, NOPASS:: derived_type_constructor => construct_bnslorene2

    !
    !-- Overloaded FUNCTION to access the fields as arrays and as values
    !

    GENERIC, PUBLIC:: get_field => get_fa, get_fv
    !# GENERIC PROCEDURE, overloded to access the [[bnslorene]]-member variables
    !  as arrays and as values
    PROCEDURE::       get_fa    => get_field_array
    !! Access the [[bnslorene]]-member arrays
    PROCEDURE::       get_fv    => get_field_value
    !! Access the components of the [[bnslorene]]-member arrays

    !
    !-- FUNCTIONS that access member variables
    !
    PROCEDURE:: get_eos1_id => get_eos1_loreneid
    !! Returns the |lorene| identifier for the EOS of star 1
    PROCEDURE:: get_eos2_id => get_eos2_loreneid
    !! Returns the |lorene| identifier for the EOS of star 2

    PROCEDURE:: return_eos_parameters => get_eos_parameters

    PROCEDURE, PUBLIC:: get_eos1_loreneid
    !! Returns [[bnslorene:eos1_loreneid]]
    PROCEDURE, PUBLIC:: get_eos2_loreneid
    !! Returns [[bnslorene:eos2_loreneid]]

    PROCEDURE, PUBLIC:: get_bns_identifier
    !! Returns [[bnslorene:bns_identifier]]

    !PROCEDURE, PUBLIC:: get_bns_ptr

    !PROCEDURE:: derived_type_destructor => destruct_bnslorene

    !PROCEDURE:: derived_type_destructor => destruct_bnslorene
    FINAL:: destruct_bnslorene
    !! Finalizer (Destructor) of a [[bnslorene]] object

  END TYPE bnslorene

  !
  !-- Interface of the TYPE bnslorene (i.e., declaration of the constructor)
  !-- (see https://dannyvanpoucke.be/oop-fortran-tut4-en/)
  !
  !INTERFACE bnslorene
  !!! Interface of TYPE [[bnslorene]]
  !
  !  !MODULE PROCEDURE:: construct_bnslorene
  !  MODULE PROCEDURE:: construct_bnslorene
  !  !! Constructs a [[bnslorene]] object
  !
  !END INTERFACE bnslorene

  !
  !-- Interfaces of the constructor and destructor of the TYPE bnslorene
  !
  INTERFACE


   ! MODULE FUNCTION construct_bnslorene( &!derived_type,
   ! filename ) RESULT( foo )
   ! !! Constructs a [[bnslorene]] object
   !
   !   CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: filename
   !   !! |lorene| binary file containing the spectral BNS ID
   !   !CLASS(bnslorene):: derived_type
   !   !! Constructed [[bnslorene]] object
   !   CLASS(idbase), ALLOCATABLE:: foo
   !
   ! END FUNCTION construct_bnslorene


    MODULE SUBROUTINE construct_bnslorene( derived_type, filename )
    !# Prints a summary of the physical properties the system
    !  to the standard output and, optionally, to a formatted file whose name
    !  is given as the optional argument `filename`

      CHARACTER(LEN=*), INTENT( IN ), OPTIONAL :: filename
      !! |lorene| binary file containing the spectral BNS ID
      CLASS(bnslorene), INTENT( OUT ):: derived_type
      !! Constructed [[bnslorene]] object

    END SUBROUTINE construct_bnslorene


    MODULE SUBROUTINE destruct_bnslorene( this )
    !! Destruct a [[bnslorene]] object

      TYPE(bnslorene), INTENT( IN OUT ):: this
      !! [[bnslorene]] object to be destructed

    END SUBROUTINE destruct_bnslorene

  END INTERFACE

  !
  !-- Interfaces of the methods of the TYPE bnslorene
  !-- Their implementations are in submodule_bnslorene_methods.f90
  !
  INTERFACE


    !------------------------------!
    !--  OVERRIDING SUBROUTINES  --!
    !------------------------------!


    MODULE SUBROUTINE print_summary_bnslorene( this, filename )
    !# Prints a summary of the physical properties of the |bns| produced by
    !  |lorene| to the standard output and, optionally, to a formatted file
    !  whose name is given as the optional argument `filename`


      CLASS(bnslorene), INTENT( IN ):: this
      CHARACTER( LEN= * ), INTENT( INOUT ), OPTIONAL:: filename
      !! Name of the formatted file to print the summary to

    END SUBROUTINE print_summary_bnslorene

    !
    !-- SUBROUTINES
    !
    MODULE SUBROUTINE construct_binary( this, resu_file )
    !! Interface of the subroutine that constructs the |lorene| |binns| object

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),                     INTENT( IN OUT )      :: this
      !> |lorene| binary file containing the spectral BNS ID
      CHARACTER(KIND= C_CHAR, LEN=*), INTENT( IN ), OPTIONAL:: resu_file

    END SUBROUTINE construct_binary


    MODULE SUBROUTINE destruct_binary( this )
    !! Destructs a |lorene| |binns| object

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene), INTENT( IN OUT ):: this

    END SUBROUTINE destruct_binary


    MODULE SUBROUTINE allocate_bnslorene_memory( this, d )
    !! Allocates allocatable arrays member of a [[bnslorene]] object

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene), INTENT( IN OUT ):: this
      !> Dimension of the arrays
      INTEGER,    INTENT( IN )    :: d

    END SUBROUTINE allocate_bnslorene_memory


    MODULE SUBROUTINE deallocate_bnslorene_memory( this )
    !! Deallocates allocatable arrays member of a [[bnslorene]] object

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene), INTENT( IN OUT ):: this

    END SUBROUTINE deallocate_bnslorene_memory


    MODULE SUBROUTINE read_id_params( this )
    !! Imports the BNS parameters from |lorene|

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene), INTENT( IN OUT ):: this

    END SUBROUTINE read_id_params


    MODULE SUBROUTINE print_id_params( this )
    !! Prints the BNS parameters to the standard output

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene), INTENT( IN OUT ):: this

    END SUBROUTINE print_id_params


  !  MODULE SUBROUTINE integrate_baryon_mass_density( this, center, radius, &
  !                                                   central_density, &
  !                                                   dr, dth, dphi, &
  !                                                   mass, mass_profile, &
  !                                                   mass_profile_idx )
  !  !# Integrates the |lorene| baryon mass density to compute the radial mass
  !  !  profile. TODO: Improve integration algorithm.
  !
  !    !> [[bnslorene]] object which this PROCEDURE is a member of
  !    CLASS(bnslorene), INTENT( IN OUT )      :: this
  !    !& Array to store the indices for array mass_profile, sorted so that
  !    !  mass_profile[mass_profile_idx] is in increasing order
  !    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT( IN OUT ):: mass_profile_idx
  !    !> Center of the star
  !    DOUBLE PRECISION, INTENT( IN )    :: center
  !    !> Central density of the star
  !    DOUBLE PRECISION, INTENT( IN )    :: central_density
  !    !> Radius of the star
  !    DOUBLE PRECISION, INTENT( IN )    :: radius
  !    !> Integration steps
  !    DOUBLE PRECISION, INTENT( IN )    :: dr, dth, dphi
  !    !> Integrated mass of the star
  !    DOUBLE PRECISION, INTENT( IN OUT ):: mass
  !    !> Array storing the radial mass profile of the star
  !    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT( IN OUT ):: &
  !                                     mass_profile
  !
  !  END SUBROUTINE integrate_baryon_mass_density


    MODULE SUBROUTINE read_id_int( this, n, x, y, z )
    !! Stores the ID in the [[bnslorene]] member arrays

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),                     INTENT( IN OUT ):: this
      INTEGER, INTENT( IN ):: n
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: x
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: y
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: z

    END SUBROUTINE read_id_int

    ! BE CAREFUL! Look at the following page:
    !
    ! https://www.ibm.com/support/knowledgecenter/SSAT4T_15.1.5/com.ibm.xlf1515.lelinux.doc/language_ref/allocobj.html

    ! where you can find the following statement,
    !
    ! "On procedure entry, the allocation status of an allocatable dummy
    !  argument becomes that of the associated actual argument. If the
    !  dummy argument is INTENT(OUT) and the associated actual argument is
    !  allocated, the actual argument is deallocated on procedure invocation
    !  so that the dummy argument has an allocation status of not allocated.
    !  If the dummy argument is not INTENT(OUT) and the actual argument is
    !  allocated, the value of the dummy argument is that of the associated
    !  actual argument."

    ! Hence, the intent of allocatable array arguments  has to be IN OUT,
    ! not OUT. The array arguments are not allocatable anymore


    MODULE SUBROUTINE read_id_full( this, n, x, y, z, &
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
    !# Stores the ID in non [[bnslorene]]-member arrays with the same shape as the
    !  [[bnslorene]] member arrays

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),                     INTENT( IN OUT ):: this
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

    END SUBROUTINE read_id_full


    MODULE SUBROUTINE read_id_spacetime( this, nx, ny, nz, &
                                              pos, &
                                              lapse, &
                                              shift, &
                                              g, &
                                              ek )
    !# Stores the spacetime ID in multi-dimensional arrays needed to compute
    !  the BSSN variables and constraints

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),                           INTENT( IN OUT ):: this
      INTEGER,                              INTENT( IN )    :: nx
      INTEGER,                              INTENT( IN )    :: ny
      INTEGER,                              INTENT( IN )    :: nz
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN )    :: pos
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: lapse
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: shift
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: g
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: ek

    END SUBROUTINE read_id_spacetime


    MODULE SUBROUTINE read_id_hydro( this, nx, ny, nz, &
                                             pos, &
                                             baryon_density, &
                                             energy_density, &
                                             specific_energy, &
                                             pressure, &
                                             u_euler )
    !# Stores the hydro ID in the arrays needed to compute the constraints
    !  on the refined mesh

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),                           INTENT( IN OUT ):: this
      INTEGER,                              INTENT( IN )    :: nx
      INTEGER,                              INTENT( IN )    :: ny
      INTEGER,                              INTENT( IN )    :: nz
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN )    :: pos
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: baryon_density
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: energy_density
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: specific_energy
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: pressure
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: u_euler

    END SUBROUTINE read_id_hydro


    MODULE SUBROUTINE read_id_particles( this, n, x, y, z, &
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

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),                     INTENT( IN OUT ):: this
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

    END SUBROUTINE read_id_particles


    MODULE SUBROUTINE read_id_mass_b( this, x, y, z, &
                                        g, &
                                        baryon_density, &
                                        gamma_euler )
    !! Stores the hydro ID in the arrays needed to compute the baryon mass

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),       INTENT( IN OUT ):: this
      DOUBLE PRECISION, INTENT( IN )    :: x
      DOUBLE PRECISION, INTENT( IN )    :: y
      DOUBLE PRECISION, INTENT( IN )    :: z
      DOUBLE PRECISION, DIMENSION(6), INTENT( OUT ):: g
      DOUBLE PRECISION, INTENT( OUT ):: baryon_density
      DOUBLE PRECISION, INTENT( OUT ):: gamma_euler

    END SUBROUTINE read_id_mass_b


    MODULE SUBROUTINE read_id_k( this, n, x, y, z,&
                                         k_xx, k_xy, k_xz, &
                                         k_yy, k_yz, k_zz )
   !! Stores the components of the extrinsic curvature in arrays

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),                     INTENT( IN OUT ):: this
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

    END SUBROUTINE read_id_k


    !
    !-- FUNCTIONS
    !
    MODULE FUNCTION read_bnslorene_mass_density( this, x, y, z ) RESULT( res )
    !! Returns the |lorene| baryon mass density at a point \((x,y,z)\)

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),     INTENT( IN )         :: this
      !> \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: x
      !> \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: y
      !> \(z\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: z
      !> Baryon mass density at \((x,y,z)\)
      DOUBLE PRECISION:: res

    END FUNCTION read_bnslorene_mass_density


    MODULE FUNCTION read_bnslorene_pressure( this, x, y, z ) RESULT( res )
    !! Returns the |lorene| pressure at a point \((x,y,z)\)

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene), INTENT( IN )       :: this
      !> \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: x
      !> \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: y
      !> \(z\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: z
      !> Pressure at \((x,y,z)\)
      DOUBLE PRECISION:: res

    END FUNCTION read_bnslorene_pressure


    MODULE FUNCTION read_spatial_metric( this, x, y, z ) RESULT( res )
    !# Returns the |lorene| conformally flat spatial metric component
    !  \(g_{xx}=g_{yy}=g_{zz}\) at a point \((x,y,z)\)

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),     INTENT( IN )       :: this
      !> \(x\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT( IN ), VALUE:: x
      !> \(y\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT( IN ), VALUE:: y
      !> \(z\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT( IN ), VALUE:: z
      !> \(g_{xx}=g_{yy}=g_{zz}\) at \((x,y,z)\)
      REAL(C_DOUBLE):: res

    END FUNCTION read_spatial_metric


    MODULE FUNCTION is_hydro_positive( this, x, y, z ) RESULT( res )
    !# Returns 1 if the energy density or the specific energy or the pressure
    !  are negative, 0 otherwise

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),     INTENT( IN )       :: this
      !> \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: x
      !> \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: y
      !> \(z\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: z
      !& `.TRUE.` if the energy density or the specific energy or the pressure
      !  are negative, `.FALSE.` otherwise
      LOGICAL:: res

    END FUNCTION is_hydro_positive


    MODULE FUNCTION get_field_array( this, field ) RESULT( field_array )
    !! Returns the [[bnslorene]] member arrays named field

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),          INTENT( IN )             :: this
      !> Name of the desired [[bnslorene]] member array
      CHARACTER( LEN= : ), INTENT( IN ), ALLOCATABLE:: field
      !> Desired [[bnslorene]] member array
      DOUBLE PRECISION, DIMENSION(:),    ALLOCATABLE:: field_array

    END FUNCTION get_field_array


    MODULE FUNCTION get_field_value( this, field, n ) RESULT( field_value )
    !! Returns the component n of the [[bnslorene]] member arrays named field

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),    INTENT( IN )             :: this
      !> Name of the desired [[bnslorene]] member array
      CHARACTER( LEN= : ), INTENT( IN ), ALLOCATABLE:: field
      !> Component of the desired [[bnslorene]] member array
      INTEGER,             INTENT( IN )             :: n
      !> Component n of the desired [[bnslorene]] member array
      DOUBLE PRECISION                              :: field_value

    END FUNCTION get_field_value


    MODULE FUNCTION get_bns_identifier( this )

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene), INTENT( IN ):: this
      ! Result
      DOUBLE PRECISION:: get_bns_identifier

    END FUNCTION get_bns_identifier


    MODULE FUNCTION get_eos1_loreneid( this )

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene), INTENT( IN ):: this
      ! Result
      INTEGER:: get_eos1_loreneid

    END FUNCTION get_eos1_loreneid


    MODULE FUNCTION get_eos2_loreneid( this )

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene), INTENT( IN ):: this
      ! Result
      INTEGER:: get_eos2_loreneid

    END FUNCTION get_eos2_loreneid


    MODULE SUBROUTINE get_eos_parameters( this, i_matter, eos_params )

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene), INTENT( IN ):: this
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT):: eos_params
      !# Array containing the parameters of the |eos| for the `i_matter`-th
      !  matter object

    END SUBROUTINE get_eos_parameters


    MODULE SUBROUTINE correct_adm_linear_momentum &
    ( npart, pos, nlrf, u, pr, vel_u, theta, nstar, nu, g_xx, g_xy, g_xz,  &
      g_yy, g_yz, g_zz, lapse, shift_x, shift_y, shift_z, adm_mom_error, &
      adm_mass )
    !# Post-process the |sph| |id|; for example, correct for the residual
    !  ADM linear momentum.

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
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN)   :: g_xx
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN)   :: g_xy
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN)   :: g_xz
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN)   :: g_yy
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN)   :: g_yz
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN)   :: g_zz
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN)   :: lapse
      !! Lapse function on the particles
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN)   :: shift_x
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN)   :: shift_y
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN)   :: shift_z
      DOUBLE PRECISION, DIMENSION(3),       INTENT(IN)   :: adm_mom_error
      DOUBLE PRECISION,                     INTENT(IN)   :: adm_mass

    END SUBROUTINE correct_adm_linear_momentum


    MODULE SUBROUTINE nothing( this, flag, switch )
    !# Procedure that does nothing. It is used to instantiate a deferred
    !  idbase procedure which is not needed in TYPE [[bnslorene]].
    !  It also serves as a placeholder in case the idbase procedure
    !  will be needed in the future.

      CLASS(bnslorene), INTENT(INOUT)      :: this
      INTEGER,          INTENT(IN)         :: flag
      !! Identifies what kind of initialization has to be done
      LOGICAL,          INTENT(IN), OPTIONAL:: switch
      !! If `.TRUE.`, switch to a different initialization

    END SUBROUTINE nothing


    !MODULE FUNCTION get_bns_ptr( this )
    !
    !  ! Argument
    !  CLASS(bnslorene), INTENT( IN ):: this
    !  ! Result
    !  TYPE(C_PTR):: get_bns_ptr
    !
    !END FUNCTION get_bns_ptr


  END INTERFACE


  !---------------------------------------------------------------------!
  !--  PRIVATE interfaces to the methods of |lorene|'s class |binns|  --!
  !---------------------------------------------------------------------!


  PRIVATE:: construct_bin_ns, get_lorene_id, get_lorene_id_spacetime, &
            get_lorene_id_particles, get_lorene_id_mass_b, &
            get_lorene_id_hydro, get_lorene_id_k, get_lorene_mass_density, &
            get_lorene_pressure, get_lorene_spatial_metric, &
            positive_hydro, get_lorene_id_params, destruct_bin_ns


  INTERFACE


    FUNCTION construct_bin_ns( c_resu_file ) RESULT( optr ) &
      BIND(C, NAME= "construct_bin_ns")

      !***********************************************
      !
      !# Interface to the |lorene| method of class
      !  |binns| with the same name, that constructs
      !  the |lorene| |binns| object
      !
      !  FT
      !
      !***********************************************

      IMPORT :: C_PTR, C_CHAR

      IMPLICIT NONE

      !& C string of the name of the |lorene| binary file storing the spectral
      !  BNS ID
      CHARACTER(KIND= C_CHAR), DIMENSION(*), INTENT(IN), OPTIONAL :: &
                                                              c_resu_file
      !> C pointer pointing to the constructed |lorene| |binns| object
      TYPE(C_PTR) :: optr

    END FUNCTION construct_bin_ns


    SUBROUTINE get_lorene_id( optr, &
                              x, y, z, &
                              lapse, &
                              shift_x, shift_y, shift_z, &
                              g_diag, &
                              k_xx, k_xy, k_xz, &
                              k_yy, k_yz, k_zz, &
                              baryon_density, &
                              energy_density, &
                              specific_energy, &
                              pressure, &
                              v_euler_x, v_euler_y, v_euler_z ) &
      BIND(C, NAME= "get_lorene_id")

      !*************************************************
      !
      !# Interface to the |lorene| method of class
      !  |binns| with the same name, that reads the full
      !  |lorene| ID at the specified point.
      !  That is, reads the metric fields, the
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
      !  FT
      !
      !*************************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      !> C pointer pointing to a |lorene| |binns| object
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
      REAL(C_DOUBLE), INTENT(OUT)       :: g_diag
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

    END SUBROUTINE get_lorene_id


    SUBROUTINE get_lorene_id_spacetime( optr, &
                                        x, y, z, &
                                        lapse, &
                                        shift_x, shift_y, shift_z, &
                                        g_diag, &
                                        k_xx, k_xy, k_xz, &
                                        k_yy, k_yz, k_zz ) &
      BIND(C, NAME= "get_lorene_id_spacetime")

      !*************************************************
      !
      !# Interface to the |lorene| method of class
      !  |binns| with the same name, that reads the
      !  metric fields and the components
      !  of the extrinsic curvature [c/km] from |lorene|,
      !  at the specified point
      !
      !  FT
      !
      !*************************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      !> C pointer pointing to a |lorene| |binns| object
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
      REAL(C_DOUBLE), INTENT(OUT)       :: g_diag
      REAL(C_DOUBLE), INTENT(OUT)       :: k_xx
      REAL(C_DOUBLE), INTENT(OUT)       :: k_xy
      REAL(C_DOUBLE), INTENT(OUT)       :: k_xz
      REAL(C_DOUBLE), INTENT(OUT)       :: k_yy
      REAL(C_DOUBLE), INTENT(OUT)       :: k_yz
      REAL(C_DOUBLE), INTENT(OUT)       :: k_zz

    END SUBROUTINE get_lorene_id_spacetime


    SUBROUTINE get_lorene_id_particles( optr, &
                                        x, y, z, &
                                        lapse, &
                                        shift_x, shift_y, shift_z, &
                                        g_diag, &
                                        baryon_density, &
                                        energy_density, &
                                        specific_energy, &
                                        pressure, &
                                        v_euler_x, v_euler_y, v_euler_z ) &
      BIND(C, NAME= "get_lorene_id_particles")

      !**********************************************
      !
      !# Interface to the |lorene| method of class
      !  |binns| with the same name, that reads the
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
      !  FT
      !
      !**********************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      !> C pointer pointing to a |lorene| |binns| object
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
      REAL(C_DOUBLE), INTENT(OUT)       :: g_diag
      REAL(C_DOUBLE), INTENT(OUT)       :: baryon_density
      REAL(C_DOUBLE), INTENT(OUT)       :: energy_density
      REAL(C_DOUBLE), INTENT(OUT)       :: specific_energy
      REAL(C_DOUBLE), INTENT(OUT)       :: pressure
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_x
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_y
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_z

    END SUBROUTINE get_lorene_id_particles


    SUBROUTINE get_lorene_id_mass_b( optr, &
                                     x, y, z, &
                                     g_diag, &
                                     baryon_density, &
                                     gamma_euler ) &
      BIND(C, NAME= "get_lorene_id_mass_b")

      !************************************************
      !
      !# Interface to the |lorene| method of class
      !  |binns| with the same name, that reads the
      !  hydro fields and the metric fields
      !  from |lorene|, at the specified point,
      !  needed to compute the baryon mass.
      !
      !  - shift vector [c]
      !  - baryon mass density [kg m^{-3}]
      !  - fluid 3-velocity with respect to the
      !    Eulerian observer [c]
      !
      !  FT
      !
      !************************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      !> C pointer pointing to a |lorene| |binns| object
      TYPE(C_PTR),    INTENT(IN), VALUE :: optr
      !> \(x\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: x
      !> \(y\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: y
      !> \(z\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: z
      !> \(g_{xx}=g_{yy}=g_{zz}\) at \(x,y,z\)
      REAL(C_DOUBLE), INTENT(OUT)       :: g_diag
      !> Baryon mass density at \(x,y,z\)
      REAL(C_DOUBLE), INTENT(OUT)       :: baryon_density
      !& Relative Lorentz factor between the 4-velocity of the fluid
      !  wrt the Eulerian observer and the 4-velocity of the Eulerian observer
      !  at \(x,y,z\)
      REAL(C_DOUBLE), INTENT(OUT)       :: gamma_euler

    END SUBROUTINE get_lorene_id_mass_b


    SUBROUTINE get_lorene_id_hydro( optr, &
                                    x, y, z, &
                                    baryon_density, &
                                    energy_density, &
                                    specific_energy, &
                                    pressure, &
                                    v_euler_x, v_euler_y, v_euler_z ) &
      BIND(C, NAME= "get_lorene_id_hydro")

      !***********************************************
      !
      !# Interface to the |lorene| method of class
      !  |binns| with the same name, that reads the
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
      !  FT
      !
      !***********************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      !> C pointer pointing to a |lorene| |binns| object
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

    END SUBROUTINE get_lorene_id_hydro


    SUBROUTINE get_lorene_id_k( optr, &
                                x, y, z, &
                                k_xx, k_xy, k_xz, &
                                k_yy, k_yz, k_zz ) &
      BIND(C, NAME= "get_lorene_id_k")

      !***********************************************
      !
      !# Interface to the |lorene| method of class
      !  |binns| with the same name, that reads the
      !  components of the extrinsic
      !  curvature [c/km] from |lorene|, at the
      !  specified point
      !
      !  FT
      !
      !***********************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      !> C pointer pointing to a |lorene| |binns| object
      TYPE(C_PTR),    INTENT(IN), VALUE :: optr
      !> \(x\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: x
      !> \(y\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: y
      !> \(z\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: z
      REAL(C_DOUBLE), INTENT(OUT)       :: k_xx
      REAL(C_DOUBLE), INTENT(OUT)       :: k_xy
      REAL(C_DOUBLE), INTENT(OUT)       :: k_xz
      REAL(C_DOUBLE), INTENT(OUT)       :: k_yy
      REAL(C_DOUBLE), INTENT(OUT)       :: k_yz
      REAL(C_DOUBLE), INTENT(OUT)       :: k_zz

    END SUBROUTINE get_lorene_id_k


    FUNCTION get_lorene_mass_density( optr, x, y, z ) RESULT( res ) &
      BIND(C, NAME= "get_mass_density")

      !********************************************
      !
      !# Interface to the |lorene| method of class
      !  |binns| with the same name, that returns
      !  the baryon mass density \([\mathrm{kg}\,
      !  \mathrm{m}^{-3}]\) from |lorene|,
      !  at the specified point
      !
      !  FT
      !
      !********************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      !> C pointer pointing to a |lorene| |binns| object
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

    END FUNCTION get_lorene_mass_density


    FUNCTION get_lorene_pressure( optr, x, y, z ) RESULT( res ) &
      BIND(C, NAME= "get_bns_pressure")

      !********************************************
      !
      !# Interface to the |lorene| method of class
      !  |binns| with the same name, that returns
      !  the pressure \([\mathrm{kg}\,
      !  c^2 \mathrm{m}^{-3}]\) from |lorene|,
      !  at the specified point
      !
      !  FT 11.02.2022
      !
      !********************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      !> C pointer pointing to a |lorene| |binns| object
      TYPE(C_PTR),    INTENT(IN),  VALUE :: optr
      !> \(x\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: x
      !> \(y\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: y
      !> \(z\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: z
      !& Pressure \([\mathrm{kg}\,c^2\, \mathrm{m}^{-3}]\) at the desired
      !  point \((x,y,z)\)
      REAL(C_DOUBLE) :: res

    END FUNCTION get_lorene_pressure


    FUNCTION get_lorene_spatial_metric( optr, x, y, z ) RESULT( res ) &
      BIND(C, NAME= "get_lorene_id_g")

      !************************************************
      !
      !# Interface to the |lorene| method of class
      !  |binns| with the same name, that returns the
      !  diagonal components of the metric,
      !  all equal to the |lorene| conformal factor to
      !  the 4th power.
      !
      !  FT
      !
      !************************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      !> C pointer pointing to a |lorene| |binns| object
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

    END FUNCTION get_lorene_spatial_metric


    FUNCTION positive_hydro( optr, x, y, z ) RESULT( res ) &
      BIND(C, NAME= "is_hydro_positive")

      !************************************************
      !
      !# Interface to the |lorene| method of class
      !  |binns| with the same name, that returns 1
      !  if the energy density is positive,
      !  and if the specific energy is positive,
      !  and if the pressure is positive,
      !  at the specified point; it returns 0 otherwise
      !
      !  FT 12.03.2021
      !
      !************************************************

      IMPORT :: C_INT, C_DOUBLE, C_PTR

      IMPLICIT NONE

      !> C pointer pointing to a |lorene| |binns| object
      TYPE(C_PTR),    INTENT(IN),  VALUE :: optr
      !> \(x\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: x
      !> \(y\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: y
      !> \(z\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: z
      !& 1 if the energy density or the specific energy or the pressure
      !  are positve, 0 otherwise
      INTEGER(C_INT) :: res

    END FUNCTION positive_hydro


    SUBROUTINE get_lorene_id_params( optr, &
                                     angular_vel, &
                                     distance, &
                                     distance_com, &
                                     mass1, &
                                     mass2, &
                                     mass_grav1, &
                                     mass_grav2, &
                                     adm_mass, &
                                     linear_momentum_x, &
                                     linear_momentum_y, &
                                     linear_momentum_z, &
                                     angular_momentum_x, &
                                     angular_momentum_y, &
                                     angular_momentum_z, &
                                     area_radius1, &
                                     radius1_x_comp, &
                                     radius1_y, &
                                     radius1_z, &
                                     radius1_x_opp, &
                                     center1_x, &
                                     barycenter1_x, &
                                     area_radius2, &
                                     radius2_x_comp, &
                                     radius2_y, &
                                     radius2_z, &
                                     radius2_x_opp, &
                                     center2_x, &
                                     barycenter2_x, &
                                     ent_center1, &
                                     nbar_center1, &
                                     rho_center1, &
                                     energy_density_center1, &
                                     specific_energy_center1, &
                                     pressure_center1, &
                                     ent_center2, &
                                     nbar_center2, &
                                     rho_center2, &
                                     energy_density_center2, &
                                     specific_energy_center2, &
                                     pressure_center2, &
                                     eos1, &
                                     eos2, &
                                     eos1_id, &
                                     eos2_id, &
                                     gamma_1, &
                                     kappa_1, &
                                     gamma_2, &
                                     kappa_2, &
                                     npeos_1, &
                                     gamma0_1, &
                                     gamma1_1, &
                                     gamma2_1, &
                                     gamma3_1, &
                                     kappa0_1, &
                                     kappa1_1, &
                                     kappa2_1, &
                                     kappa3_1, &
                                     logP1_1,  &
                                     logRho0_1, &
                                     logRho1_1, &
                                     logRho2_1, &
                                     npeos_2,  &
                                     gamma0_2, &
                                     gamma1_2, &
                                     gamma2_2, &
                                     gamma3_2, &
                                     kappa0_2, &
                                     kappa1_2, &
                                     kappa2_2, &
                                     kappa3_2, &
                                     logP1_2,  &
                                     logRho0_2, &
                                     logRho1_2, &
                                     logRho2_2 ) &
      BIND(C, NAME= "get_lorene_id_params")

      !**********************************************
      !
      !# Interface to the |lorene| method of class
      !  |binns| with the same name, that stores
      !  the physical parameters of the binary
      !  system from |lorene| in the desired variables
      !
      !  FT
      !
      !**********************************************

      IMPORT :: C_INT, C_DOUBLE, C_PTR, C_CHAR

      IMPLICIT NONE

      !> C pointer pointing to a |lorene| |binns| object
      TYPE(C_PTR),    INTENT(IN), VALUE :: optr
      REAL(C_DOUBLE), INTENT(OUT)       :: angular_vel
      REAL(C_DOUBLE), INTENT(OUT)       :: distance
      REAL(C_DOUBLE), INTENT(OUT)       :: distance_com
      REAL(C_DOUBLE), INTENT(OUT)       :: mass1
      REAL(C_DOUBLE), INTENT(OUT)       :: mass2
      REAL(C_DOUBLE), INTENT(OUT)       :: mass_grav1
      REAL(C_DOUBLE), INTENT(OUT)       :: mass_grav2
      REAL(C_DOUBLE), INTENT(OUT)       :: adm_mass
      REAL(C_DOUBLE), INTENT(OUT)       :: linear_momentum_x
      REAL(C_DOUBLE), INTENT(OUT)       :: linear_momentum_y
      REAL(C_DOUBLE), INTENT(OUT)       :: linear_momentum_z
      REAL(C_DOUBLE), INTENT(OUT)       :: angular_momentum_x
      REAL(C_DOUBLE), INTENT(OUT)       :: angular_momentum_y
      REAL(C_DOUBLE), INTENT(OUT)       :: angular_momentum_z
      REAL(C_DOUBLE), INTENT(OUT)       :: area_radius1
      REAL(C_DOUBLE), INTENT(OUT)       :: radius1_x_comp
      REAL(C_DOUBLE), INTENT(OUT)       :: radius1_y
      REAL(C_DOUBLE), INTENT(OUT)       :: radius1_z
      REAL(C_DOUBLE), INTENT(OUT)       :: radius1_x_opp
      REAL(C_DOUBLE), INTENT(OUT)       :: center1_x
      REAL(C_DOUBLE), INTENT(OUT)       :: barycenter1_x
      REAL(C_DOUBLE), INTENT(OUT)       :: area_radius2
      REAL(C_DOUBLE), INTENT(OUT)       :: radius2_x_comp
      REAL(C_DOUBLE), INTENT(OUT)       :: radius2_y
      REAL(C_DOUBLE), INTENT(OUT)       :: radius2_z
      REAL(C_DOUBLE), INTENT(OUT)       :: radius2_x_opp
      REAL(C_DOUBLE), INTENT(OUT)       :: center2_x
      REAL(C_DOUBLE), INTENT(OUT)       :: barycenter2_x
      REAL(C_DOUBLE), INTENT(OUT)       :: ent_center1
      REAL(C_DOUBLE), INTENT(OUT)       :: nbar_center1
      REAL(C_DOUBLE), INTENT(OUT)       :: rho_center1
      REAL(C_DOUBLE), INTENT(OUT)       :: energy_density_center1
      REAL(C_DOUBLE), INTENT(OUT)       :: specific_energy_center1
      REAL(C_DOUBLE), INTENT(OUT)       :: pressure_center1
      REAL(C_DOUBLE), INTENT(OUT)       :: ent_center2
      REAL(C_DOUBLE), INTENT(OUT)       :: nbar_center2
      REAL(C_DOUBLE), INTENT(OUT)       :: rho_center2
      REAL(C_DOUBLE), INTENT(OUT)       :: energy_density_center2
      REAL(C_DOUBLE), INTENT(OUT)       :: specific_energy_center2
      REAL(C_DOUBLE), INTENT(OUT)       :: pressure_center2
      CHARACTER(KIND=C_CHAR), DIMENSION(100), INTENT(OUT):: eos1
      CHARACTER(KIND=C_CHAR), DIMENSION(100), INTENT(OUT):: eos2
      INTEGER(C_INT)                    :: eos1_id
      INTEGER(C_INT)                    :: eos2_id
      REAL(C_DOUBLE), INTENT(OUT)       :: gamma_1
      REAL(C_DOUBLE), INTENT(OUT)       :: kappa_1
      REAL(C_DOUBLE), INTENT(OUT)       :: gamma_2
      REAL(C_DOUBLE), INTENT(OUT)       :: kappa_2
      INTEGER(C_INT)                    :: npeos_1
      REAL(C_DOUBLE), INTENT(OUT)       :: gamma0_1
      REAL(C_DOUBLE), INTENT(OUT)       :: gamma1_1
      REAL(C_DOUBLE), INTENT(OUT)       :: gamma2_1
      REAL(C_DOUBLE), INTENT(OUT)       :: gamma3_1
      REAL(C_DOUBLE), INTENT(OUT)       :: kappa0_1
      REAL(C_DOUBLE), INTENT(OUT)       :: kappa1_1
      REAL(C_DOUBLE), INTENT(OUT)       :: kappa2_1
      REAL(C_DOUBLE), INTENT(OUT)       :: kappa3_1
      REAL(C_DOUBLE), INTENT(OUT)       :: logP1_1
      REAL(C_DOUBLE), INTENT(OUT)       :: logRho0_1
      REAL(C_DOUBLE), INTENT(OUT)       :: logRho1_1
      REAL(C_DOUBLE), INTENT(OUT)       :: logRho2_1
      INTEGER(C_INT)                    :: npeos_2
      REAL(C_DOUBLE), INTENT(OUT)       :: gamma0_2
      REAL(C_DOUBLE), INTENT(OUT)       :: gamma1_2
      REAL(C_DOUBLE), INTENT(OUT)       :: gamma2_2
      REAL(C_DOUBLE), INTENT(OUT)       :: gamma3_2
      REAL(C_DOUBLE), INTENT(OUT)       :: kappa0_2
      REAL(C_DOUBLE), INTENT(OUT)       :: kappa1_2
      REAL(C_DOUBLE), INTENT(OUT)       :: kappa2_2
      REAL(C_DOUBLE), INTENT(OUT)       :: kappa3_2
      REAL(C_DOUBLE), INTENT(OUT)       :: logP1_2
      REAL(C_DOUBLE), INTENT(OUT)       :: logRho0_2
      REAL(C_DOUBLE), INTENT(OUT)       :: logRho1_2
      REAL(C_DOUBLE), INTENT(OUT)       :: logRho2_2

    END SUBROUTINE get_lorene_id_params


    SUBROUTINE destruct_bin_ns( optr ) &
      BIND(C, NAME= "destruct_bin_ns")

      !**********************************************
      !
      !# Interface to the |lorene| method of class
      !  |binns| with the same name, that destructs
      !  the |lorene| |binns| object
      !
      ! FT
      !
      !**********************************************

      IMPORT :: C_PTR

      IMPLICIT NONE

      !> C pointer pointing to the |lorene| |binns| object to destruct
      TYPE(C_PTR), INTENT(IN), VALUE :: optr

    END SUBROUTINE destruct_bin_ns


  END INTERFACE


END MODULE bns_lorene
