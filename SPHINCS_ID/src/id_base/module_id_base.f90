! File:         module_id_base.f90
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

MODULE id_base

  !***********************************************************
  !
  !# This MODULE contains the definition of TYPE idbase,
  !  which is an ABSTRACT TYPE representing any possible
  !  type of initial data (|id|) to be set up for |sphincsbssn|.
  !  That is, a binary neutron star system, a rotating
  !  star, a binary black hole system, etc.
  !
  !  PROCEDURES and variables shared by all these types
  !  of |id| should belong to TYPE idbase, as
  !  they are inherited by its EXTENDED TYPES that
  !  represent more specific types of |id|.
  !
  !  FT 24.09.2021
  !
  !***********************************************************


  USE timing, ONLY: timer


  IMPLICIT NONE


  !*************************************************************
  !                                                            *
  !              Definition of TYPE idbase                     *
  !                                                            *
  !   This ABSTRACT TYPE represents a generic |id| for           *
  !   |sphincsbssn| (binary neutron star, rotating star...).   *
  !                                                            *
  !*************************************************************

  TYPE, ABSTRACT:: idbase
  !# Represents a generic |id| for |sphincsbssn| (binary neutron star, rotating
  !  star, etc.)


    PRIVATE


    INTEGER:: n_matter= 0
    !# Number of matter objects belonging the physical system.
    !  For example, n_objects= 2 for a binary system of stars, and n_objects= 1
    !  for a single star or for a binary system of a black hole and a star.


    LOGICAL:: one_lapse
    !# Logical variable that determines if the lapse function \(\alpha=1\),
    !  i.e., if the geodesic gauge is to be used
    LOGICAL:: zero_shift
    !! Logical variable that determines if the shift \(\beta^i=0\)

    LOGICAL:: cold_system
    !# `.TRUE.` if the system is at zero temperature (no thermal component);
    !  `.FALSE.` otherwise


    TYPE(timer), PUBLIC:: construction_timer
    !! Timer that times the construction of the appropriate object


    CONTAINS


    !---------------------------!
    !--  DEFERRED PROCEDURES  --!
    !---------------------------!

    !
    !-- PROCEDURES to read the value of a field at a point
    !

    PROCEDURE(read_double_at_pos),        DEFERRED:: read_mass_density
    !# Returns the baryon mass density at the given point

    PROCEDURE(read_integer_at_pos),       DEFERRED:: test_position
    !# Returns 1 if the position has physically acceptable properties,
    !  0 otherwise

    !
    !-- PROCEDURES to read the value of several fields at several points
    !

    PROCEDURE(read_id_full_int),          DEFERRED:: read_id_full
    !# Reads the full |id|

    PROCEDURE(read_id_particles_int),     DEFERRED:: read_id_particles
    !! Reads the hydro |id| needed to compute the SPH |id|

    PROCEDURE(read_id_mass_b_int),        DEFERRED:: read_id_mass_b
    !! Reads the hydro |id| needed to compute the baryon mass

    PROCEDURE(read_id_spacetime_int),     DEFERRED:: read_id_spacetime
    !# Reads the spacetime |id| needed to compute
    !  the BSSN variables and constraints

    PROCEDURE(read_id_hydro_int),         DEFERRED:: read_id_hydro
    !# Reads the hydro |id| needed to compute the constraints on the refined
    !  mesh

    PROCEDURE(read_id_k_int),             DEFERRED:: read_id_k
    !! Reads the components of the extrinsic curvature

    !
    !-- PROCEDURES returning the values of some parameters of a matter object
    !

    PROCEDURE(return_spatial_extent_int), DEFERRED:: return_spatial_extent
    !# Returns the spatial extent of the matter objects,
    !  returning the array of 6 numbers
    !\(x_{\rm min},x_{\rm max},y_{\rm min},y_{\rm max},z_{\rm min},z_{\rm max}\)

    PROCEDURE(return_double_parameter),   DEFERRED:: return_mass
    !! Returns the masses of the matter objects.

    PROCEDURE(return_position),           DEFERRED:: return_center
    !! Returns the centers of the matter objects.

    PROCEDURE(return_position),           DEFERRED:: return_barycenter
    !! Returns the barycenters (centers of mass) of the matter objects.

    PROCEDURE(return_eos_parameters_int), DEFERRED:: return_eos_parameters
    !# Returns the identification number of the |eos| of the matter objects.
    !  @todo Set up a convention for the identification number

    PROCEDURE(return_string_parameter),   DEFERRED:: return_eos_name
    !! Returns the name of the |eos| of the matter objects.

    !
    !-- PROCEDURE that prints a summary of the physical properties the system
    !-- to the standard output and, optionally, to a formatted file
    !

    PROCEDURE(print_summary_int),         DEFERRED:: print_summary
    !# Prints a summary of the physical properties the system
    !  to the standard output and, optionally, to a formatted file whose name
    !  is given as optional argument

    !
    !-- Constructors and destructors of derived types
    !

    PROCEDURE(derived_type_constructor_int), DEFERRED:: derived_type_constructor
    !# Constructs a TYPE that extends [[idbase]]


    !-------------------------------!
    !--  NON-DEFERRED PROCEDURES  --!
    !-------------------------------!


    PROCEDURE, NON_OVERRIDABLE:: sanity_check
    !# Checks that [[idbase:n_matter]] and the sizes returned by
    ![[idbase:return_spatial_extent]] and [[idbase:get_total_spatial_extent]]
    !  are acceptable. It is called by initialize, after the constructor of the
    !  derived type.


    PROCEDURE, NON_OVERRIDABLE:: initialize
    !# This PROCEDURE calls the constructor of the [[idbase]]-extended type
    !  and the SUBROUTINE [[idbase:sanity_check]] afterwards. It is recommended
    !  to use this SUBROUTINE to construct objects of [[idbase]]-extended type
    !  since the sanity check is performed automatically.


    PROCEDURE, NON_OVERRIDABLE:: get_total_spatial_extent
    !# Returns the spatial extent of the physical system considered,
    !  as the array of 6 numbers
    !\(x_{\rm min},x_{\rm max},y_{\rm min},y_{\rm max},z_{\rm min},z_{\rm max}\)


    PROCEDURE, NON_OVERRIDABLE:: set_n_matter
    !# Sets [[idbase:n_matter]], the number of matter objects in the
    !  physical system, to a value


    PROCEDURE, NON_OVERRIDABLE:: get_n_matter
    !# Returns [[idbase:n_matter]], the number of matter objects in the
    !  physical system


    PROCEDURE, NON_OVERRIDABLE:: set_one_lapse
    !# Sets [[idbase:one_lapse]], the `LOGICAL` variable that determines if
    !  the lapse \(\alpha=1\), i.e., if the geodesic gauge is to be used


    PROCEDURE, NON_OVERRIDABLE:: get_one_lapse
    !# Returns [[idbase:one_lapse]], the `LOGICAL` variable that determines if
    ! the lapse function \(\alpha=1\), i.e., if the geodesic gauge is to be used


    PROCEDURE, NON_OVERRIDABLE:: set_zero_shift
    !# Sets [[idbase:zero_shift]], the `LOGICAL` variable that determines if
    !  the shift \(\beta^i=0\)


    PROCEDURE, NON_OVERRIDABLE:: get_zero_shift
    !# Returns [[idbase:zero_shift]], the `LOGICAL` variable that determines if
    !  the shift \(\beta^i=0\)


    PROCEDURE, NON_OVERRIDABLE:: set_cold_system
    !# Sets [[idbase:cold_system]], the `LOGICAL` variable at specifies if
    !  the system is cold (no thermal component)


    PROCEDURE, NON_OVERRIDABLE:: get_cold_system
    !# Returns [[idbase:cold_system]], the `LOGICAL` variable at specifies if
    !  the system is cold (no thermal component)


    PROCEDURE, NON_OVERRIDABLE:: check_i_matter
    !# Checks that the given index is between 1 and [[idbase:n_matter]],
    !  included. If not, it stops the execution of the program.


    PROCEDURE:: integrate_baryon_mass_density
    !# Integrates the baryon mass density over a matter object, using spherical
    !  coordinates, and computes its radial profile inside the star


  !  PROCEDURE:: estimate_lengthscale_to_resolve_spacetime
    !# Estimate typical length scales, one per each matter object, to be
    !  resolved by the mesh. maybe this one makes more sense if it's deferred,
    !  because it's based on a approximation of the metric? well, this
    !  approximation would be the same for all systems...or would it?


  !  PROCEDURE:: estimate_lengthscale_to_resolve_sph
    !# Estimate typical length scales, one per each matter object, to be
    !  resolved by the particles


  END TYPE idbase


  ABSTRACT INTERFACE


   ! FUNCTION derived_type_constructor_int( &!derived_type,
   ! filename ) RESULT( foo )
   ! !#
   !
   !   IMPORT:: idbase
   !   CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: filename
   !   !CLASS(idbase):: derived_type
   !   CLASS(idbase), ALLOCATABLE:: foo
   !
   ! END FUNCTION derived_type_constructor_int


    SUBROUTINE derived_type_constructor_int( derived_type, filename )
    !# Prints a summary of the physical properties the system
    !  to the standard output and, optionally, to a formatted file whose name
    !  is given as the optional argument `filename`

      IMPORT:: idbase
      CHARACTER(LEN=*), INTENT( IN ), OPTIONAL :: filename
      !! |lorene| binary file containing the spectral DRS |id|
      CLASS(idbase), INTENT( OUT ):: derived_type
      !! Constructed [[diffstarlorene]] object

    END SUBROUTINE derived_type_constructor_int


   ! SUBROUTINE derived_type_destructor_int( derived_type )
   ! !# Prints a summary of the physical properties the system
   ! !  to the standard output and, optionally, to a formatted file whose name
   ! !  is given as the optional argument `filename`
   !
   !   IMPORT:: idbase
   !   CLASS(idbase), INTENT( IN OUT ):: derived_type
   !   !! Constructed [[diffstarlorene]] object
   !
   ! END SUBROUTINE derived_type_destructor_int


    FUNCTION read_double_at_pos( THIS, x, y, z ) RESULT( res )
    !# INTERFACE for a PROCEDURE that returns a DOUBLE PRECISION at a given
    !  position

      IMPORT:: idbase
      CLASS(idbase), INTENT( IN )          :: THIS
      !! Object of class [[idbase]] which this PROCEDURE is a member of
      DOUBLE PRECISION, INTENT( IN ), VALUE:: x
      !! \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: y
      !! \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: z
      !! \(z\) coordinate of the desired point
      DOUBLE PRECISION:: res
      !! Real number at \((x,y,z)\)

    END FUNCTION read_double_at_pos


    FUNCTION read_integer_at_pos( THIS, x, y, z ) RESULT( res )
    !# INTERFACE for a PROCEDURE that returns an INTEGER at a given position

      IMPORT:: idbase
      CLASS(idbase), INTENT( IN )          :: THIS
      !! Object of class [[idbase]] which this PROCEDURE is a member of
      DOUBLE PRECISION, INTENT( IN ), VALUE:: x
      !! \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: y
      !! \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: z
      !! \(z\) coordinate of the desired point
      INTEGER:: res
      !! Integer at \((x,y,z)\)

    END FUNCTION read_integer_at_pos


    FUNCTION return_double_parameter( THIS, i_matter ) RESULT( res )
    !# INTERFACE for a PROCEDURE that returns a DOUBLE PRECISION

      IMPORT:: idbase
      CLASS(idbase), INTENT( IN ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION:: res
      !! Real number. Parameter of the `i_matter`-th matter object

    END FUNCTION return_double_parameter


    FUNCTION return_position( THIS, i_matter ) RESULT( res )
    !# INTERFACE for a PROCEDURE that returns a DOUBLE PRECISION

      IMPORT:: idbase
      CLASS(idbase), INTENT( IN ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION, DIMENSION(3):: res
      !# Centers of the matter objects. The first index runs over the matter
      !  objects, the second index over \((x,y,z)\).

    END FUNCTION return_position


    FUNCTION return_integer_parameter( THIS, i_matter ) RESULT( res )
    !# INTERFACE for a PROCEDURE that returns an INTEGER

      IMPORT:: idbase
      CLASS(idbase), INTENT( IN ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose parameter is to return
      INTEGER:: res
      !! Real number. Parameter of the `i_matter`-th matter object

    END FUNCTION return_integer_parameter


    SUBROUTINE return_eos_parameters_int( THIS, i_matter, eos_params )
    !# INTERFACE for a PROCEDURE that returns an array containing the
    !  parametersf the |eos| for the matter objects

      IMPORT:: idbase
      CLASS(idbase), INTENT( IN ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT):: eos_params
      !# Array containing the parameters of the |eos| for the `i_matter`-th
      !  matter object

    END SUBROUTINE return_eos_parameters_int


    FUNCTION return_string_parameter( THIS, i_matter ) RESULT( string )
    !# INTERFACE for a PROCEDURE that returns a CHARACTER( LEN= : )

      IMPORT:: idbase
      CLASS(idbase), INTENT( IN ):: THIS
      !! [[idbase]] object which this PROCEDURE is a member of
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose string is to return
      CHARACTER( LEN= : ), ALLOCATABLE:: string

    END FUNCTION return_string_parameter


    SUBROUTINE read_id_mass_b_int( THIS, x, y, z, &
                                   g, &
                                   baryon_density, &
                                   gamma_euler )
      !# INTERFACE or the SUBROUTINE reading the hydro |id| needed to compute
      !  the baryon mass

      IMPORT:: idbase
      CLASS(idbase),    INTENT( IN OUT ):: THIS
      !! Object of class [[idbase]] which this PROCEDURE is a member of
      DOUBLE PRECISION, INTENT( IN )    :: x
      DOUBLE PRECISION, INTENT( IN )    :: y
      DOUBLE PRECISION, INTENT( IN )    :: z
      DOUBLE PRECISION, DIMENSION(6), INTENT( OUT ):: g
      DOUBLE PRECISION, INTENT( OUT ):: baryon_density
      DOUBLE PRECISION, INTENT( OUT ):: gamma_euler

    END SUBROUTINE read_id_mass_b_int


    SUBROUTINE read_id_full_int( THIS, n, x, y, z, &
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
     !# INTERFACE or the SUBROUTINE reading the full |id|

      IMPORT:: idbase
      !> [[idbase]] object which this PROCEDURE is a member of
      CLASS(idbase),                  INTENT( IN OUT ):: THIS
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

    END SUBROUTINE read_id_full_int


    SUBROUTINE read_id_spacetime_int( THIS, nx, ny, nz, &
                                              pos, &
                                              lapse, &
                                              shift, &
                                              g, &
                                              ek )
     !# INTERFACE or the SUBROUTINE reading the spacetime |id|

      IMPORT:: idbase
      !> [[idbase]] object which this PROCEDURE is a member of
      CLASS(idbase),                        INTENT( IN OUT ):: THIS
      INTEGER,                              INTENT( IN )    :: nx
      INTEGER,                              INTENT( IN )    :: ny
      INTEGER,                              INTENT( IN )    :: nz
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN )    :: pos
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: lapse
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: shift
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: g
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: ek

    END SUBROUTINE read_id_spacetime_int


    SUBROUTINE read_id_hydro_int( THIS, nx, ny, nz, &
                                             pos, &
                                             baryon_density, &
                                             energy_density, &
                                             specific_energy, &
                                             pressure, &
                                             u_euler )
    !# INTERFACE or the SUBROUTINE reading the the hydro |id| needed to compute
    !  the constraints on the refined mesh

      IMPORT:: idbase
      !> [[idbase]] object which this PROCEDURE is a member of
      CLASS(idbase),                        INTENT( IN OUT ):: THIS
      INTEGER,                              INTENT( IN )    :: nx
      INTEGER,                              INTENT( IN )    :: ny
      INTEGER,                              INTENT( IN )    :: nz
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN )    :: pos
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: baryon_density
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: energy_density
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: specific_energy
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: pressure
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: u_euler

    END SUBROUTINE read_id_hydro_int


    SUBROUTINE read_id_particles_int( THIS, n, x, y, z, &
                                      lapse, &
                                      shift_x, shift_y, shift_z, &
                                      g_xx, g_xy, g_xz, &
                                      g_yy, g_yz, g_zz, &
                                      baryon_density, &
                                      energy_density, &
                                      specific_energy, &
                                      pressure, &
                                      u_euler_x, u_euler_y, u_euler_z )
    !# INTERFACE or the SUBROUTINE reading the hydro |id| needed to compute the
    !  SPH |id|

      IMPORT:: idbase
      !> [[idbase]] object which this PROCEDURE is a member of
      CLASS(idbase),                     INTENT( IN OUT ):: THIS
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
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: baryon_density
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: energy_density
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: specific_energy
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: pressure
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: u_euler_x
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: u_euler_y
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN OUT ):: u_euler_z

    END SUBROUTINE read_id_particles_int


    SUBROUTINE read_id_k_int( THIS, n, x, y, z,&
                              k_xx, k_xy, k_xz, &
                              k_yy, k_yz, k_zz )
    !# INTERFACE or the SUBROUTINE reading the components of the extrinsic
    !  curvature

      IMPORT:: idbase
      !> [[idbase]] object which this PROCEDURE is a member of
      CLASS(idbase),                     INTENT( IN OUT ):: THIS
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

    END SUBROUTINE read_id_k_int


  !  SUBROUTINE integrate_field_int( THIS, center, radius, &
  !                                  central_density, &
  !                                  dr, dth, dphi, &
  !                                  mass, mass_profile, &
  !                                  mass_profile_idx )
  !  !# INTERFACE to the SUBROUTINE integrating the baryon mass density to
  !  !  compute the radial mass profile of a single star.
  !
  !    IMPORT:: idbase
  !    !> Object of class [[idbase]] which this PROCEDURE is a member of
  !    CLASS(idbase), INTENT( IN OUT )   :: THIS
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
  !    !DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT( INOUT ):: &
  !    !                                 mass_profile
  !    DOUBLE PRECISION, DIMENSION(3,0:NINT(radius/dr)), INTENT( OUT ):: &
  !                                         mass_profile
  !    !& Array to store the indices for array mass_profile, sorted so that
  !    !  mass_profile[mass_profile_idx] is in increasing order
  !    !INTEGER, DIMENSION(:), ALLOCATABLE, INTENT( INOUT ):: mass_profile_idx
  !    INTEGER, DIMENSION(0:NINT(radius/dr)), INTENT( OUT ):: mass_profile_idx
  !
  !  END SUBROUTINE integrate_field_int


    FUNCTION return_spatial_extent_int( THIS, i_matter ) RESULT( box )
    !# INTERFACE to the SUBROUTINE that detects the spatial extent of the
    !  matter objects, and returns a 6-dimensional array
    !  containing the coordinates
    !\(x_{\rm min},x_{\rm max},y_{\rm min},y_{\rm max},z_{\rm min},z_{\rm max}\)
    !  of a box **centered at the center of the object** and containing the
    !  system.

      IMPORT:: idbase
      CLASS(idbase), INTENT( IN )   :: THIS
      !! Object of class [[idbase]] which this PROCEDURE is a member of
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose string is to return
      DOUBLE PRECISION, DIMENSION(6):: box
      !# 6-dimensional array containing the coordinates
      !  \(x_{\rm min},x_{\rm max},y_{\rm min},y_{\rm max},
      !    z_{\rm min},z_{\rm max}\)
      !  of a box containing the physical system.

    END FUNCTION return_spatial_extent_int


    SUBROUTINE print_summary_int( THIS, filename )
    !# Prints a summary of the physical properties the system
    !  to the standard output and, optionally, to a formatted file whose name
    !  is given as the optional argument `filename`

      IMPORT:: idbase
      CLASS(idbase), INTENT( IN ):: THIS
      CHARACTER( LEN= * ), INTENT( INOUT ), OPTIONAL:: filename
      !! Name of the formatted file to print the summary to

    END SUBROUTINE print_summary_int


  END INTERFACE


  INTERFACE


    MODULE SUBROUTINE sanity_check( derived_type )
    !# Prints a summary of the physical properties the system
    !  to the standard output and, optionally, to a formatted file whose name
    !  is given as the optional argument `filename`

      !IMPORT:: idbase
      CLASS(idbase), INTENT( IN ):: derived_type

    END SUBROUTINE sanity_check


    MODULE SUBROUTINE initialize( derived_type, filename )
    !# Prints a summary of the physical properties the system
    !  to the standard output and, optionally, to a formatted file whose name
    !  is given as the optional argument `filename`

      CHARACTER(LEN=*), INTENT( IN ), OPTIONAL :: filename
      !! |lorene| binary file containing the spectral DRS ID
      CLASS(idbase), INTENT( OUT ):: derived_type
      !! Constructed [[diffstarlorene]] object

    END SUBROUTINE initialize


    MODULE SUBROUTINE integrate_baryon_mass_density( THIS, center, radius, &
                                                     central_density, &
                                                     dr, dth, dphi, &
                                                     mass, mass_profile, &
                                                     mass_profile_idx )
    !# INTERFACE to the SUBROUTINE integrating the baryon mass density to
    !  compute the radial mass profile of a single star.

      !> Object of class [[idbase]] which this PROCEDURE is a member of
      CLASS(idbase), INTENT( IN OUT )   :: THIS
      !> Center of the star
      DOUBLE PRECISION, INTENT( IN )    :: center
      !> Central density of the star
      DOUBLE PRECISION, INTENT( IN )    :: central_density
      !> Radius of the star
      DOUBLE PRECISION, INTENT( IN )    :: radius
      !> Integration steps
      DOUBLE PRECISION, INTENT( IN )    :: dr, dth, dphi
      !> Integrated mass of the star
      DOUBLE PRECISION, INTENT( IN OUT ):: mass
      !> Array storing the radial mass profile of the star
      !DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT( INOUT ):: &
      !                                 mass_profile
      DOUBLE PRECISION, DIMENSION(3,0:NINT(radius/dr)), INTENT( OUT ):: &
                                           mass_profile
      !& Array to store the indices for array mass_profile, sorted so that
      !  mass_profile[mass_profile_idx] is in increasing order
      !INTEGER, DIMENSION(:), ALLOCATABLE, INTENT( INOUT ):: mass_profile_idx
      INTEGER, DIMENSION(0:NINT(radius/dr)), INTENT( OUT ):: mass_profile_idx

    END SUBROUTINE integrate_baryon_mass_density


    PURE MODULE FUNCTION get_n_matter( THIS )
    !# Returns [[idbase:n_matter]], the number of matter objects in the
    !  physical system

      CLASS(idbase), INTENT( IN ):: THIS
      INTEGER:: get_n_matter
      !! [[idbase:n_matter]], the number of matter objects in the
      !  physical system

    END FUNCTION get_n_matter


    PURE MODULE SUBROUTINE set_n_matter( THIS, value )
    !# Sets [[idbase:n_matter]], the number of matter objects in the
    !  physical system, to the given value

      CLASS(idbase), INTENT( IN OUT ):: THIS
      INTEGER, INTENT( IN ):: value
      !! Value to set [[idbase:n_matter]] to

    END SUBROUTINE set_n_matter


    PURE MODULE FUNCTION get_cold_system( THIS )
    !# Returns [[idbase:cold_system]], the `LOGICAL` variable at specifies if
    !  the system is cold (no thermal component)

      CLASS(idbase), INTENT( IN ):: THIS
      LOGICAL:: get_cold_system
      !! [[idbase:cold_system]], the number of matter objects in the
      !  physical system

    END FUNCTION get_cold_system


    MODULE SUBROUTINE set_cold_system( THIS, value )
    !# Sets [[idbase:cold_system]], the `LOGICAL` variable at specifies if
    !  the system is cold (no thermal component)

      CLASS(idbase), INTENT( IN OUT ):: THIS
      LOGICAL, INTENT( IN ):: value
      !! Value to set [[idbase:cold_system]] to

    END SUBROUTINE set_cold_system


    MODULE SUBROUTINE check_i_matter( THIS, i_matter )
    !# Checks that the given index `i_matter` is between 1 and
    !  [[idbase:n_matter]], included. If not, it stops the execution of the
    !  program.

      CLASS(idbase), INTENT( IN ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      !! Value to be checked

    END SUBROUTINE check_i_matter


    MODULE FUNCTION get_total_spatial_extent( THIS ) RESULT( box )
    !# INTERFACE to the SUBROUTINE that detects the spatial extent of the
    !  physical system considered, and returns a 6-dimensional array
    !  containing the coordinates
    !\(x_{\rm min},x_{\rm max},y_{\rm min},y_{\rm max},z_{\rm min},z_{\rm max}\)
    !  of a box **centered at the center of the object** and containing the
    !  system.

      CLASS(idbase), INTENT( IN )   :: THIS
      !! Object of class [[idbase]] which this PROCEDURE is a member of
      DOUBLE PRECISION, DIMENSION(6):: box
      !# 6-dimensional array containing the coordinates
      !  \(x_{\rm min},x_{\rm max},y_{\rm min},y_{\rm max},
      !    z_{\rm min},z_{\rm max}\)
      !  of a box containing the physical system.

    END FUNCTION get_total_spatial_extent


    PURE MODULE FUNCTION get_one_lapse( THIS )
    !# Returns [[idbase:n_matter]], the number of matter objects in the
    !  physical system

      CLASS(idbase), INTENT( IN ):: THIS
      LOGICAL:: get_one_lapse
      !! [[idbase:n_matter]], the number of matter objects in the
      !  physical system

    END FUNCTION get_one_lapse


    PURE MODULE SUBROUTINE set_one_lapse( THIS, logic )
    !# Sets [[idbase:n_matter]], the number of matter objects in the
    !  physical system, to the given value

      CLASS(idbase), INTENT( IN OUT ):: THIS
      LOGICAL, INTENT( IN ):: logic
      !! Value to set [[idbase:n_matter]] to

    END SUBROUTINE set_one_lapse


    MODULE PURE FUNCTION get_zero_shift( THIS )
    !# Returns [[idbase:n_matter]], the number of matter objects in the
    !  physical system

      CLASS(idbase), INTENT( IN ):: THIS
      LOGICAL:: get_zero_shift
      !! [[idbase:n_matter]], the number of matter objects in the
      !  physical system

    END FUNCTION get_zero_shift


    PURE MODULE SUBROUTINE set_zero_shift( THIS, logic )
    !# Sets [[idbase:n_matter]], the number of matter objects in the
    !  physical system, to the given value

      CLASS(idbase), INTENT( IN OUT ):: THIS
      LOGICAL, INTENT( IN ):: logic
      !! Value to set [[idbase:n_matter]] to

    END SUBROUTINE set_zero_shift


  END INTERFACE


END MODULE id_base
