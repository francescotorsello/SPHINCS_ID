! File:         module_bns_fuka.f90
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

MODULE bns_fuka

  !***********************************************************
  !
  !#  This module contains the definition of TYPE bnsfuka,
  !   and the SUBROUTINES that bind to the methods
  !   of |fuka|'s class |binns|
  !
  !   [|fuka| official site](https://kadath.obspm.fr/fuka/#){:target="_blank"}
  !
  !***********************************************************


  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_CHAR, C_NULL_CHAR, &
                                         C_PTR, C_NULL_PTR, C_ASSOCIATED
  USE bns_base,                    ONLY: bnsbase
  USE id_base,                     ONLY: idbase
  USE utility,                     ONLY: itr, ios, err_msg, test_status, &
                                         perc, creturn, compute_g4, &
                                         determinant_sym4x4, show_progress
  USE timing,                      ONLY: timer


  IMPLICIT NONE


  !********************************************************
  !                                                       *
  !            Definition of TYPE bnsfuka                 *
  !                                                       *
  !   This class reads and stores the |fuka| |bns| |id|   *
  !                                                       *
  !********************************************************

  TYPE, EXTENDS(bnsbase):: bnsfuka
  !# TYPE representing a binary system of neutron stars (|bns|) produced with
  !  |fuka|


    PRIVATE


    !> Identifier of the bnsfuka object
    INTEGER:: bns_identifier= 0
    !> |fuka| identifiers for the EoS
    INTEGER:: eos1_fukaid, eos2_fukaid

    !
    !-- Spacetime fields
    !

    !> 1-D array storing the lapse function
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: lapse
    !> 1-D array storing the x component of the shift vector [c]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: shift_x
    !> 1-D array storing the y component of the shift vector [c]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: shift_y
    !> 1-D array storing the z component of the shift vector [c]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: shift_z
    !> 1-D array storing the xx component of the spatial metric [pure number]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: g_xx
    !> 1-D array storing the xy component of the spatial metric [pure number]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: g_xy
    !> 1-D array storing the xz component of the spatial metric [pure number]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: g_xz
    !> 1-D array storing the yy component of the spatial metric [pure number]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: g_yy
    !> 1-D array storing the yz component of the spatial metric [pure number]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: g_yz
    !> 1-D array storing the zz component of the spatial metric [pure number]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: g_zz
    !& 1-D array storing the xx component of the extrinsic curvature
    !  [c/MSun_geo]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: k_xx
    !& 1-D array storing the xy component of the extrinsic curvature
    !  [c/MSun_geo]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: k_xy
    !& 1-D array storing the xz component of the extrinsic curvature
    !  [c/MSun_geo]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: k_xz
    !& 1-D array storing the yy component of the extrinsic curvature
    !  [c/MSun_geo]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: k_yy
    !& 1-D array storing the yz component of the extrinsic curvature
    !  [c/MSun_geo]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: k_yz
    !& 1-D array storing the zz component of the extrinsic curvature
    !  [c/MSun_geo]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: k_zz

    !
    !-- Hydro fields
    !

    !> 1-D array storing the baryon mass density in the fluid frame [kg m^{-3}]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: mass_density
    !> 1-D array storing the energy density [kg c^2 m^{-3}]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: energy_density
    !> 1-D array storing the specific internal energy [c^2]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: specific_energy
    !> 1-D array storing the x component of the fluid 3-velocity with respect to
    !  the Eulerian observer [c]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: v_euler_x
    !> 1-D array storing the y component of the fluid 3-velocity with respect to
    !  the Eulerian observer [c]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: v_euler_y
    !> 1-D array storing the z component of the fluid 3-velocity with respect to
    !  the Eulerian observer [c]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: v_euler_z

    !& C pointer to the |fuka|'s |binns| object
    ! N.B. This variable is global. The pointer to the second |fuka| |binns|
    !      object will overwrite the first one, and so on.
    !      This variable stores the pointer to the last defined |fuka| |binns|
    !      object. That's why it is not freed in the destructor of a bns object.
    !      Presently, it has to be freed by the user at the end of the PROGRAM.
    !      See the last part of the PROGRAM in sphincs_id.f90, for example.
    TYPE(C_PTR):: bns_ptr


    CONTAINS


    !-------------------!
    !--  SUBROUTINES  --!
    !-------------------!

    PROCEDURE:: derived_type_constructor => construct_bnsfuka

    PROCEDURE:: construct_binary
    !! Constructs the |fuka| |binns| object

    PROCEDURE:: destruct_binary
    !! Destructs the |fuka| |binns| object

    PROCEDURE:: allocate_bnsfuka_memory
    !! Allocates memory for the [[bnsfuka]] member arrays

    PROCEDURE:: deallocate_bnsfuka_memory
    !! Deallocates memory for the [[bnsfuka]] member arrays

    PROCEDURE:: read_fuka_id_params
    !! Imports the parameters of the |bns| from |fuka|

    !PROCEDURE:: integrate_field_on_star => integrate_baryon_mass_density
    !# Integrates the |fuka| baryon mass density and computes the
    !  radial mass profile

    PROCEDURE, PUBLIC:: print_id_params
    !! Prints the parameters of the |bns| to the standard output

    PROCEDURE:: read_fuka_id_member
    !! Stores the |id| in the [[bnsfuka]] member arrays

    PROCEDURE:: read_id_full      => read_fuka_id_full
    PROCEDURE:: read_id_spacetime => read_fuka_id_spacetime
    PROCEDURE:: read_id_particles => read_fuka_id_particles
    PROCEDURE:: read_id_hydro     => read_fuka_id_hydro
    PROCEDURE:: read_id_mass_b    => read_fuka_id_mass_b
    PROCEDURE:: read_id_k         => read_fuka_id_k

    PROCEDURE:: print_summary_derived => print_summary_bnsfuka

    !-----------------!
    !--  FUNCTIONS  --!
    !-----------------!

    !> Returns the |fuka|'s mass density at the given point
    PROCEDURE:: read_mass_density => read_fuka_mass_density

    !> Returns the |fuka|'s conformally flat spatial ADM metric
    PROCEDURE:: read_fuka_spatial_metric

    !& Returns 1 if the energy density or the specific energy or the pressure
    !  are negative
    PROCEDURE:: test_position => is_hydro_positive

    !PROCEDURE, NOPASS:: derived_type_constructor => construct_bnsfuka2

    !
    !-- Overloaded FUNCTION to access the fields as arrays and as values
    !

    GENERIC, PUBLIC:: get_field => get_fa, get_fv
    !# GENERIC PROCEDURE, overloded to access the [[bnsfuka]]-member variables
    !  as arrays and as values
    PROCEDURE::       get_fa    => get_field_array
    !! Access the [[bnsfuka]]-member arrays
    PROCEDURE::       get_fv    => get_field_value
    !! Access the components of the [[bnsfuka]]-member arrays

    !
    !-- FUNCTIONS that access member variables
    !
    PROCEDURE:: get_eos1_id => get_eos1_fukaid
    !! Returns the |fuka| identifier for the EOS of star 1
    PROCEDURE:: get_eos2_id => get_eos2_fukaid
    !! Returns the |fuka| identifier for the EOS of star 2

    PROCEDURE:: return_eos_parameters => get_eos_parameters

    PROCEDURE, PUBLIC:: get_eos1_fukaid
    !! Returns [[bnsfuka:eos1_fukaid]]
    PROCEDURE, PUBLIC:: get_eos2_fukaid
    !! Returns [[bnsfuka:eos2_fukaid]]

    PROCEDURE, PUBLIC:: get_bns_identifier
    !! Returns [[bnsfuka:bns_identifier]]

    !PROCEDURE, PUBLIC:: get_bns_ptr

    !PROCEDURE:: derived_type_destructor => destruct_bnsfuka

    !PROCEDURE:: derived_type_destructor => destruct_bnsfuka
    FINAL:: destruct_bnsfuka
    !! Finalizer (Destructor) of a [[bnsfuka]] object

  END TYPE bnsfuka

  !
  !-- Interface of the TYPE bnsfuka (i.e., declaration of the constructor)
  !-- (see https://dannyvanpoucke.be/oop-fortran-tut4-en/)
  !
  !INTERFACE bnsfuka
  !!! Interface of TYPE [[bnsfuka]]
  !
  !  !MODULE PROCEDURE:: construct_bnsfuka
  !  MODULE PROCEDURE:: construct_bnsfuka
  !  !! Constructs a [[bnsfuka]] object
  !
  !END INTERFACE bnsfuka

  !
  !-- Interfaces of the constructor and destructor of the TYPE bnsfuka
  !
  INTERFACE


   ! MODULE FUNCTION construct_bnsfuka( &!derived_type,
   ! filename ) RESULT( foo )
   ! !! Constructs a [[bnsfuka]] object
   !
   !   CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: filename
   !   !! |fuka| binary file containing the spectral |bns| |id|
   !   !CLASS(bnsfuka):: derived_type
   !   !! Constructed [[bnsfuka]] object
   !   CLASS(idbase), ALLOCATABLE:: foo
   !
   ! END FUNCTION construct_bnsfuka


    MODULE SUBROUTINE construct_bnsfuka( derived_type, filename )
    !# Prints a summary of the physical properties the system
    !  to the standard output and, optionally, to a formatted file whose name
    !  is given as the optional argument `filename`

      CHARACTER(LEN=*), INTENT( IN ), OPTIONAL :: filename
      !! |fuka| binary file containing the spectral |bns| |id|
      CLASS(bnsfuka), INTENT( OUT ):: derived_type
      !! Constructed [[bnsfuka]] object

    END SUBROUTINE construct_bnsfuka


    MODULE SUBROUTINE destruct_bnsfuka( this )
    !! Destruct a [[bnsfuka]] object

      TYPE(bnsfuka), INTENT( IN OUT ):: this
      !! [[bnsfuka]] object to be destructed

    END SUBROUTINE destruct_bnsfuka

  END INTERFACE

  !
  !-- Interfaces of the methods of the TYPE bnsfuka
  !-- Their implementations are in submodule_bnsfuka_methods.f90
  !
  INTERFACE


    !------------------------------!
    !--  OVERRIDING SUBROUTINES  --!
    !------------------------------!


    MODULE SUBROUTINE print_summary_bnsfuka( this, filename )
    !# Prints a summary of the physical properties of the |bns| produced by
    !  |fuka| to the standard output and, optionally, to a formatted file
    !  whose name is given as the optional argument `filename`


      CLASS(bnsfuka), INTENT( IN ):: this
      CHARACTER( LEN= * ), INTENT( INOUT ), OPTIONAL:: filename
      !! Name of the formatted file to print the summary to

    END SUBROUTINE print_summary_bnsfuka

    !
    !-- SUBROUTINES
    !
    MODULE SUBROUTINE construct_binary( this, resu_file )
    !! Interface of the subroutine that constructs the |fuka| |binns| object

      !> [[bnsfuka]] object which this PROCEDURE is a member of
      CLASS(bnsfuka),                     INTENT( IN OUT )      :: this
      !> |fuka| binary file containing the spectral |bns| |id|
      CHARACTER(KIND= C_CHAR, LEN=*), INTENT( IN ), OPTIONAL:: resu_file

    END SUBROUTINE construct_binary


    MODULE SUBROUTINE destruct_binary( this )
    !! Destructs a |fuka| |binns| object

      !> [[bnsfuka]] object which this PROCEDURE is a member of
      CLASS(bnsfuka), INTENT( IN OUT ):: this

    END SUBROUTINE destruct_binary


    MODULE SUBROUTINE allocate_bnsfuka_memory( this, d )
    !! Allocates allocatable arrays member of a [[bnsfuka]] object

      !> [[bnsfuka]] object which this PROCEDURE is a member of
      CLASS(bnsfuka), INTENT( IN OUT ):: this
      !> Dimension of the arrays
      INTEGER,    INTENT( IN )    :: d

    END SUBROUTINE allocate_bnsfuka_memory


    MODULE SUBROUTINE deallocate_bnsfuka_memory( this )
    !! Deallocates allocatable arrays member of a [[bnsfuka]] object

      !> [[bnsfuka]] object which this PROCEDURE is a member of
      CLASS(bnsfuka), INTENT( IN OUT ):: this

    END SUBROUTINE deallocate_bnsfuka_memory


    MODULE SUBROUTINE read_fuka_id_params( this )
    !! Imports the |bns| parameters from |fuka|

      !> [[bnsfuka]] object which this PROCEDURE is a member of
      CLASS(bnsfuka), INTENT( IN OUT ):: this

    END SUBROUTINE read_fuka_id_params


    MODULE SUBROUTINE print_id_params( this )
    !! Prints the |bns| parameters to the standard output

      !> [[bnsfuka]] object which this PROCEDURE is a member of
      CLASS(bnsfuka), INTENT( IN OUT ):: this

    END SUBROUTINE print_id_params


  !  MODULE SUBROUTINE integrate_baryon_mass_density( this, center, radius, &
  !                                                   central_density, &
  !                                                   dr, dth, dphi, &
  !                                                   mass, mass_profile, &
  !                                                   mass_profile_idx )
  !  !# Integrates the |fuka| baryon mass density to compute the radial mass
  !  !  profile. TODO: Improve integration algorithm.
  !
  !    !> [[bnsfuka]] object which this PROCEDURE is a member of
  !    CLASS(bnsfuka), INTENT( IN OUT )      :: this
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


    MODULE SUBROUTINE read_fuka_id_member( this, n, x, y, z )
    !! Stores the |id| in the [[bnsfuka]] member arrays

      !> [[bnsfuka]] object which this PROCEDURE is a member of
      CLASS(bnsfuka),                 INTENT( IN OUT ):: this
      INTEGER, INTENT( IN ):: n
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: x
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: y
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: z

    END SUBROUTINE read_fuka_id_member

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


    MODULE SUBROUTINE read_fuka_id_full( this, n, x, y, z,&
                                         lapse, &
                                         shift_x, shift_y, shift_z, &
                                         g_xx, g_xy, g_xz, &
                                         g_yy, g_yz, g_zz, &
                                         k_xx, k_xy, k_xz, &
                                         k_yy, k_yz, k_zz, &
                                         baryon_density, &
                                         energy_density, &
                                         specific_energy, &
                                         u_euler_x, u_euler_y, u_euler_z )
    !# Stores the |id| in non [[bnsfuka]]-member arrays with the same shape as the
    !  [[bnsfuka]] member arrays

      !> [[bnsfuka]] object which this PROCEDURE is a member of
      CLASS(bnsfuka),                 INTENT( IN OUT ):: this
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
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: u_euler_x
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: u_euler_y
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: u_euler_z

    END SUBROUTINE read_fuka_id_full


    MODULE SUBROUTINE read_fuka_id_spacetime( this, nx, ny, nz, &
                                              pos, &
                                              lapse, &
                                              shift, &
                                              g, &
                                              ek )
    !# Stores the spacetime |id| in multi-dimensional arrays needed to compute
    !  the BSSN variables and constraints

      !> [[bnsfuka]] object which this PROCEDURE is a member of
      CLASS(bnsfuka),                       INTENT( IN OUT ):: this
      INTEGER,                              INTENT( IN )    :: nx
      INTEGER,                              INTENT( IN )    :: ny
      INTEGER,                              INTENT( IN )    :: nz
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN )    :: pos
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: lapse
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: shift
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: g
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: ek

    END SUBROUTINE read_fuka_id_spacetime


    MODULE SUBROUTINE read_fuka_id_hydro( this, nx, ny, nz, &
                                          pos, &
                                          baryon_density, &
                                          energy_density, &
                                          specific_energy, &
                                          pressure, &
                                          u_euler )
    !# Stores the hydro |id| in the arrays needed to compute the constraints
    !  on the refined mesh

      !> [[bnsfuka]] object which this PROCEDURE is a member of
      CLASS(bnsfuka),                       INTENT( IN OUT ):: this
      INTEGER,                              INTENT( IN )    :: nx
      INTEGER,                              INTENT( IN )    :: ny
      INTEGER,                              INTENT( IN )    :: nz
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN )    :: pos
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: baryon_density
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: energy_density
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: specific_energy
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: pressure
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: u_euler

    END SUBROUTINE read_fuka_id_hydro


    MODULE SUBROUTINE read_fuka_id_particles( this, n, x, y, z, &
                                              lapse, &
                                              shift_x, shift_y, shift_z, &
                                              g_xx, g_xy, g_xz, &
                                              g_yy, g_yz, g_zz, &
                                              baryon_density, &
                                              energy_density, &
                                              specific_energy, &
                                              pressure, &
                                              u_euler_x, u_euler_y, u_euler_z )
    !! Stores the hydro |id| in the arrays needed to compute the |sph| |id|

      !> [[bnsfuka]] object which this PROCEDURE is a member of
      CLASS(bnsfuka),                 INTENT( IN OUT ):: this
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

    END SUBROUTINE read_fuka_id_particles


    MODULE SUBROUTINE read_fuka_id_mass_b( this, x, y, z, &
                                           g, &
                                           baryon_density, &
                                           gamma_euler )
    !! Stores the hydro |id| in the arrays needed to compute the baryon mass

      !> [[bnsfuka]] object which this PROCEDURE is a member of
      CLASS(bnsfuka),   INTENT( IN OUT ):: this
      DOUBLE PRECISION, INTENT( IN )    :: x
      DOUBLE PRECISION, INTENT( IN )    :: y
      DOUBLE PRECISION, INTENT( IN )    :: z
      DOUBLE PRECISION, DIMENSION(6), INTENT( OUT ):: g
      DOUBLE PRECISION, INTENT( OUT ):: baryon_density
      DOUBLE PRECISION, INTENT( OUT ):: gamma_euler

    END SUBROUTINE read_fuka_id_mass_b


    MODULE SUBROUTINE read_fuka_id_k( this, n, x, y, z,&
                                      k_xx, k_xy, k_xz, &
                                      k_yy, k_yz, k_zz )
   !! Stores the components of the extrinsic curvature in arrays

      !> [[bnsfuka]] object which this PROCEDURE is a member of
      CLASS(bnsfuka),                 INTENT( IN OUT ):: this
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

    END SUBROUTINE read_fuka_id_k


    !
    !-- FUNCTIONS
    !
    MODULE FUNCTION read_fuka_mass_density( this, x, y, z ) RESULT( res )
    !! Returns the |fuka| baryon mass density at a point \((x,y,z)\)

      !> [[bnsfuka]] object which this PROCEDURE is a member of
      CLASS(bnsfuka),     INTENT( IN )         :: this
      !> \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: x
      !> \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: y
      !> \(z\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: z
      !> Baryon mass density at \((x,y,z)\)
      DOUBLE PRECISION:: res

    END FUNCTION read_fuka_mass_density


    MODULE FUNCTION read_fuka_spatial_metric( this, x, y, z ) RESULT( res )
    !# Returns the |fuka| conformally flat spatial metric component
    !  \(g_{xx}=g_{yy}=g_{zz}\) at a point \((x,y,z)\)

      !> [[bnsfuka]] object which this PROCEDURE is a member of
      CLASS(bnsfuka),     INTENT( IN )       :: this
      !> \(x\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT( IN ), VALUE:: x
      !> \(y\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT( IN ), VALUE:: y
      !> \(z\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT( IN ), VALUE:: z
      !> \(g_{xx}=g_{yy}=g_{zz}\) at \((x,y,z)\)
      REAL(C_DOUBLE):: res

    END FUNCTION read_fuka_spatial_metric


    MODULE FUNCTION is_hydro_positive( this, x, y, z ) RESULT( res )
    !# Returns 1 if the energy density or the specific energy or the pressure
    !  are negative, 0 otherwise

      !> [[bnsfuka]] object which this PROCEDURE is a member of
      CLASS(bnsfuka),     INTENT( IN )       :: this
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
    !! Returns the [[bnsfuka]] member arrays named field

      !> [[bnsfuka]] object which this PROCEDURE is a member of
      CLASS(bnsfuka),          INTENT( IN )             :: this
      !> Name of the desired [[bnsfuka]] member array
      CHARACTER( LEN= : ), INTENT( IN ), ALLOCATABLE:: field
      !> Desired [[bnsfuka]] member array
      DOUBLE PRECISION, DIMENSION(:),    ALLOCATABLE:: field_array

    END FUNCTION get_field_array


    MODULE FUNCTION get_field_value( this, field, n ) RESULT( field_value )
    !! Returns the component n of the [[bnsfuka]] member arrays named field

      !> [[bnsfuka]] object which this PROCEDURE is a member of
      CLASS(bnsfuka),          INTENT( IN )             :: this
      !> Name of the desired [[bnsfuka]] member array
      CHARACTER( LEN= : ), INTENT( IN ), ALLOCATABLE:: field
      !> Component of the desired [[bnsfuka]] member array
      INTEGER,             INTENT( IN )             :: n
      !> Component n of the desired [[bnsfuka]] member array
      DOUBLE PRECISION                              :: field_value

    END FUNCTION get_field_value


    MODULE FUNCTION get_bns_identifier( this )

      !> [[bnsfuka]] object which this PROCEDURE is a member of
      CLASS(bnsfuka), INTENT( IN ):: this
      ! Result
      DOUBLE PRECISION:: get_bns_identifier

    END FUNCTION get_bns_identifier


    MODULE FUNCTION get_eos1_fukaid( this )

      !> [[bnsfuka]] object which this PROCEDURE is a member of
      CLASS(bnsfuka), INTENT( IN ):: this
      ! Result
      INTEGER:: get_eos1_fukaid

    END FUNCTION get_eos1_fukaid


    MODULE FUNCTION get_eos2_fukaid( this )

      !> [[bnsfuka]] object which this PROCEDURE is a member of
      CLASS(bnsfuka), INTENT( IN ):: this
      ! Result
      INTEGER:: get_eos2_fukaid

    END FUNCTION get_eos2_fukaid


    MODULE SUBROUTINE get_eos_parameters( this, i_matter, eos_params )

      !> [[bnsfuka]] object which this PROCEDURE is a member of
      CLASS(bnsfuka), INTENT( IN ):: this
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT):: eos_params
      !# Array containing the parameters of the |eos| for the `i_matter`-th
      !  matter object

    END SUBROUTINE get_eos_parameters


    !MODULE FUNCTION get_bns_ptr( this )
    !
    !  ! Argument
    !  CLASS(bnsfuka), INTENT( IN ):: this
    !  ! Result
    !  TYPE(C_PTR):: get_bns_ptr
    !
    !END FUNCTION get_bns_ptr


  END INTERFACE


  !------------------------------------------------------------------!
  !--  PRIVATE interfaces to the methods of |fuka|'s class |binns|  --!
  !------------------------------------------------------------------!


  PRIVATE:: construct_bns_fuka, get_fuka_id, get_fuka_id_spacetime, &
            get_fuka_id_particles, get_fuka_id_mass_b, &
            get_fuka_id_hydro, get_fuka_id_k, get_fuka_mass_density, &
            get_fuka_spatial_metric, positive_hydro, get_fuka_id_params, &
            destruct_bns_fuka


  INTERFACE


    FUNCTION construct_bns_fuka( fuka_file ) RESULT( optr ) &
      BIND(C, NAME= "construct_bns_fuka")

      !***********************************************
      !
      !# Interface to the |fuka| method of class
      !  |binns| with the same name, that constructs
      !  the |fuka| |binns| object
      !
      !  FT
      !
      !***********************************************

      IMPORT :: C_PTR, C_CHAR

      IMPLICIT NONE

      CHARACTER(KIND= C_CHAR), DIMENSION(*), INTENT(IN) :: fuka_file
      !# C string of the name of the |fuka| binary file storing the spectral
      !  |bns| |id|
      TYPE(C_PTR) :: optr
      !! C pointer pointing to the constructed |fuka| BNS_export object

    END FUNCTION construct_bns_fuka


    SUBROUTINE get_fuka_id( optr, &
                              x, y, z, &
                              lapse, &
                              shift_x, shift_y, shift_z, &
                              psi4, &
                              k_xx, k_xy, k_xz, &
                              k_yy, k_yz, k_zz, &
                              mass_density, &
                              energy_density, &
                              pressure, &
                              v_euler_x, v_euler_y, v_euler_z ) &
      BIND(C, NAME= "get_fuka_id")

      !*************************************************
      !
      !# Interface to the |fuka| method of class
      !  |binns| with the same name, that reads the full
      !  |fuka| |id| at the specified point.
      !  That is, read_fukas the metric fields, the
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

      !> C pointer pointing to a |fuka| |binns| object
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
      REAL(C_DOUBLE), INTENT(OUT)       :: psi4
      REAL(C_DOUBLE), INTENT(OUT)       :: k_xx
      REAL(C_DOUBLE), INTENT(OUT)       :: k_xy
      REAL(C_DOUBLE), INTENT(OUT)       :: k_xz
      REAL(C_DOUBLE), INTENT(OUT)       :: k_yy
      REAL(C_DOUBLE), INTENT(OUT)       :: k_yz
      REAL(C_DOUBLE), INTENT(OUT)       :: k_zz
      REAL(C_DOUBLE), INTENT(OUT)       :: mass_density
      REAL(C_DOUBLE), INTENT(OUT)       :: energy_density
      REAL(C_DOUBLE), INTENT(OUT)       :: pressure
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_x
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_y
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_z

    END SUBROUTINE get_fuka_id


    SUBROUTINE get_fuka_id_spacetime( optr, &
                                      x, y, z, &
                                      lapse, &
                                      shift_x, shift_y, shift_z, &
                                      psi4, &
                                      k_xx, k_xy, k_xz, &
                                      k_yy, k_yz, k_zz ) &
      BIND(C, NAME= "get_fuka_id_spacetime")

      !*************************************************
      !
      !# Interface to the |fuka| method of class
      !  |binns| with the same name, that reads the
      !  metric fields and the components
      !  of the extrinsic curvature [c/km] from |fuka|,
      !  at the specified point
      !
      !  FT
      !
      !*************************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      !> C pointer pointing to a |fuka| |binns| object
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
      REAL(C_DOUBLE), INTENT(OUT)       :: psi4
      REAL(C_DOUBLE), INTENT(OUT)       :: k_xx
      REAL(C_DOUBLE), INTENT(OUT)       :: k_xy
      REAL(C_DOUBLE), INTENT(OUT)       :: k_xz
      REAL(C_DOUBLE), INTENT(OUT)       :: k_yy
      REAL(C_DOUBLE), INTENT(OUT)       :: k_yz
      REAL(C_DOUBLE), INTENT(OUT)       :: k_zz

    END SUBROUTINE get_fuka_id_spacetime


    SUBROUTINE get_fuka_id_particles( optr, &
                                      x, y, z, &
                                      lapse, &
                                      shift_x, shift_y, shift_z, &
                                      psi4, &
                                      mass_density, &
                                      energy_density, &
                                      pressure, &
                                      v_euler_x, v_euler_y, v_euler_z ) &
      BIND(C, NAME= "get_fuka_id_particles")

      !**********************************************
      !
      !# Interface to the |fuka| method of class
      !  |binns| with the same name, that reads the
      !  hydro fields and the metric fields *
      !  from |fuka|, at the specified point
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

      !> C pointer pointing to a |fuka| |binns| object
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
      REAL(C_DOUBLE), INTENT(OUT)       :: psi4
      REAL(C_DOUBLE), INTENT(OUT)       :: mass_density
      REAL(C_DOUBLE), INTENT(OUT)       :: energy_density
      REAL(C_DOUBLE), INTENT(OUT)       :: pressure
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_x
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_y
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_z

    END SUBROUTINE get_fuka_id_particles


    SUBROUTINE get_fuka_id_mass_b( optr, &
                                   x, y, z, &
                                   psi4, &
                                   mass_density, &
                                   gamma_euler ) &
      BIND(C, NAME= "get_fuka_id_massb")

      !************************************************
      !
      !# Interface to the |fuka| method of class
      !  |binns| with the same name, that reads the
      !  hydro fields and the metric fields
      !  from |fuka|, at the specified point,
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

      !> C pointer pointing to a |fuka| |binns| object
      TYPE(C_PTR),    INTENT(IN), VALUE :: optr
      !> \(x\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: x
      !> \(y\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: y
      !> \(z\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: z
      !> \(g_{xx}=g_{yy}=g_{zz}\) at \(x,y,z\)
      REAL(C_DOUBLE), INTENT(OUT)       :: psi4
      !> Baryon mass density at \(x,y,z\)
      REAL(C_DOUBLE), INTENT(OUT)       :: mass_density
      !& Relative Lorentz factor between the 4-velocity of the fluid
      !  wrt the Eulerian observer and the 4-velocity of the Eulerian observer
      !  at \(x,y,z\)
      REAL(C_DOUBLE), INTENT(OUT)       :: gamma_euler

    END SUBROUTINE get_fuka_id_mass_b


    SUBROUTINE get_fuka_id_hydro( optr, &
                                  x, y, z, &
                                  mass_density, &
                                  energy_density, &
                                  pressure, &
                                  v_euler_x, v_euler_y, v_euler_z ) &
      BIND(C, NAME= "get_fuka_id_hydro")

      !***********************************************
      !
      !# Interface to the |fuka| method of class
      !  |binns| with the same name, that reads the
      !  hydro fields from |fuka|, at the
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

      !> C pointer pointing to a |fuka| |binns| object
      TYPE(C_PTR),    INTENT(IN), VALUE :: optr
      !> \(x\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: x
      !> \(y\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: y
      !> \(z\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT(IN), VALUE :: z
      REAL(C_DOUBLE), INTENT(OUT)       :: mass_density
      REAL(C_DOUBLE), INTENT(OUT)       :: energy_density
      REAL(C_DOUBLE), INTENT(OUT)       :: pressure
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_x
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_y
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_z

    END SUBROUTINE get_fuka_id_hydro


    SUBROUTINE get_fuka_id_k( optr, &
                                x, y, z, &
                                k_xx, k_xy, k_xz, &
                                k_yy, k_yz, k_zz ) &
      BIND(C, NAME= "get_fuka_id_k")

      !***********************************************
      !
      !# Interface to the |fuka| method of class
      !  |binns| with the same name, that reads the
      !  components of the extrinsic
      !  curvature [c/km] from |fuka|, at the
      !  specified point
      !
      !  FT
      !
      !***********************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      !> C pointer pointing to a |fuka| |binns| object
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

    END SUBROUTINE get_fuka_id_k


    FUNCTION get_fuka_mass_density( optr, x, y, z ) RESULT( res ) &
      BIND(C, NAME= "get_fuka_mass_density")

      !********************************************
      !
      !#
      !
      !  FT
      !
      !********************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      !> C pointer pointing to a |fuka| |binns| object
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

    END FUNCTION get_fuka_mass_density


    FUNCTION get_fuka_spatial_metric( optr, x, y, z ) RESULT( res ) &
      BIND(C, NAME= "get_fuka_g")

      !************************************************
      !
      !#
      !
      !  FT
      !
      !************************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      !> C pointer pointing to a |fuka| |binns| object
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

    END FUNCTION get_fuka_spatial_metric


    FUNCTION positive_hydro( optr, x, y, z ) RESULT( res ) &
      BIND(C, NAME= "is_fuka_hydro_positive")

      !************************************************
      !
      !#
      !
      !  FT
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


    SUBROUTINE get_fuka_id_params( optr,                   &
                                   angular_vel,            &
                                   distance,               &
                                   mass1,                  &
                                   mass2,                  &
                                   massg1,                 &
                                   massg2,                 &
                                   adm_mass,               &
                                   komar_mass,             &
                                   adm_linear_momentum_x,  &
                                   adm_linear_momentum_y,  &
                                   adm_linear_momentum_z,  &
                                   adm_angular_momentum_z, &
                                   COMx,                   &
                                   COMy,                   &
                                   COMz,                   &
                                   area_radius1,           &
                                   center1_x,              &
                                   area_radius2,           &
                                   center2_x,              &
                                   ent_center1,            &
                                   rho_center1,            &
                                   energy_density_center1, &
                                   ent_center2,            &
                                   rho_center2,            &
                                   energy_density_center2, &
                                   eos_type )              &
      BIND(C, NAME= "get_fuka_id_params")

      !**********************************************
      !
      !#
      !
      !  FT
      !
      !**********************************************

      IMPORT :: C_INT, C_DOUBLE, C_PTR, C_CHAR

      IMPLICIT NONE

      TYPE(C_PTR),    INTENT(IN), VALUE :: optr
      !! C pointer pointing to a |fuka| BNS_export object
      REAL(C_DOUBLE), INTENT(OUT)       :: angular_vel
      REAL(C_DOUBLE), INTENT(OUT)       :: distance
      REAL(C_DOUBLE), INTENT(OUT)       :: mass1
      REAL(C_DOUBLE), INTENT(OUT)       :: mass2
      REAL(C_DOUBLE), INTENT(OUT)       :: massg1
      REAL(C_DOUBLE), INTENT(OUT)       :: massg2
      REAL(C_DOUBLE), INTENT(OUT)       :: adm_mass
      REAL(C_DOUBLE), INTENT(OUT)       :: komar_mass
      REAL(C_DOUBLE), INTENT(OUT)       :: adm_linear_momentum_x
      REAL(C_DOUBLE), INTENT(OUT)       :: adm_linear_momentum_y
      REAL(C_DOUBLE), INTENT(OUT)       :: adm_linear_momentum_z
      REAL(C_DOUBLE), INTENT(OUT)       :: adm_angular_momentum_z
      REAL(C_DOUBLE), INTENT(OUT)       :: COMx
      REAL(C_DOUBLE), INTENT(OUT)       :: COMy
      REAL(C_DOUBLE), INTENT(OUT)       :: COMz
      REAL(C_DOUBLE), INTENT(OUT)       :: area_radius1
      REAL(C_DOUBLE), INTENT(OUT)       :: center1_x
      REAL(C_DOUBLE), INTENT(OUT)       :: area_radius2
      REAL(C_DOUBLE), INTENT(OUT)       :: center2_x
      REAL(C_DOUBLE), INTENT(OUT)       :: ent_center1
      REAL(C_DOUBLE), INTENT(OUT)       :: rho_center1
      REAL(C_DOUBLE), INTENT(OUT)       :: energy_density_center1
      REAL(C_DOUBLE), INTENT(OUT)       :: ent_center2
      REAL(C_DOUBLE), INTENT(OUT)       :: rho_center2
      REAL(C_DOUBLE), INTENT(OUT)       :: energy_density_center2
      CHARACTER(KIND=C_CHAR), DIMENSION(100), INTENT(OUT):: eos_type

    END SUBROUTINE get_fuka_id_params


    SUBROUTINE destruct_bns_fuka( optr ) &
      BIND(C, NAME= "destruct_bns_fuka")

      !**********************************************
      !
      !# Interface to the |fuka| method of class
      !  |binns| with the same name, that destructs
      !  the |fuka| |binns| object
      !
      ! FT
      !
      !**********************************************

      IMPORT :: C_PTR

      IMPLICIT NONE

      !> C pointer pointing to the |fuka| |binns| object to destruct
      TYPE(C_PTR), INTENT(IN), VALUE :: optr

    END SUBROUTINE destruct_bns_fuka


  END INTERFACE


END MODULE bns_fuka
