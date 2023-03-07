! File:         module_id_base.f90
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


  TYPE surface
    LOGICAL:: is_known= .FALSE.
    DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE:: points
    !# Array containing the coordinates of the stars' surfaces
    !  The first index runs over the stars; the second and third run over
    !  the surface points (number of points for \(\theta\) and \(\varphi\));
    !  the fourth index runs overthe Cartesian coordinates of the points
  END TYPE surface


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
    !  For example, n_matter= 2 for a binary system of stars, and n_matter= 1
    !  for a single star or for a binary system of a black hole and a star.

    TYPE(surface), PUBLIC, DIMENSION(:), ALLOCATABLE:: surfaces

    LOGICAL:: one_lapse
    !# Logical variable that determines if the lapse function \(\alpha=1\),
    !  i.e., if the geodesic gauge is to be used
    LOGICAL:: zero_shift
    !! Logical variable that determines if the shift \(\beta^i=0\)

    LOGICAL:: cold_system
    !# `.TRUE.` if the system is at zero temperature (no thermal component);
    !  `.FALSE.` otherwise

    LOGICAL:: estimate_length_scale
    !# `.TRUE.` if a typical length scale equal to the ratio of a field over
    !  its gradient should be computed (usually, the field is the pressure);
    !  `.FALSE.` otherwise


    TYPE(timer), PUBLIC:: construction_timer
    !! Timer that times the construction of the appropriate object


    PROCEDURE(), POINTER, NOPASS, PUBLIC:: finalize_sph_id_ptr
    !# Pointer to a procedure that finalize the |sph| |id|; for example,
    !  correct for the residual ADM linear momentum.


    CONTAINS


    !---------------------------!
    !--  DEFERRED PROCEDURES  --!
    !---------------------------!

    !
    !-- PROCEDURES to read the value of a field at a point
    !

    PROCEDURE(read_double_at_pos),        DEFERRED:: read_mass_density
    !# Returns the baryon mass density from the |id| at the given point

    PROCEDURE(read_double_at_pos),        DEFERRED:: read_pressure
    !# Returns the pressure from te |id| at the given point

    PROCEDURE(read_logical_at_pos),       DEFERRED:: test_position
    !# Returns `.TRUE.` if the position has physically acceptable properties,
    !  `.FALSE.` otherwise

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

    PROCEDURE(return_double_at_object),   DEFERRED:: return_mass
    !! Returns the masses of the matter objects.

    PROCEDURE(return_double),             DEFERRED:: return_adm_mass
    !! Returns the ADM mass of the system

    PROCEDURE(return_position),           DEFERRED:: return_center
    !! Returns the centers of the matter objects.

    PROCEDURE(return_position),           DEFERRED:: return_barycenter
    !! Returns the barycenters (centers of mass) of the matter objects.

    PROCEDURE(return_eos_parameters_int), DEFERRED:: return_eos_parameters
    !# Returns the identification number of the |eos| of the matter objects.
    !  @todo Set up a convention for the identification number

    PROCEDURE(return_string_parameter),   DEFERRED:: return_eos_name
    !! Returns the name of the |eos| of the matter objects.


    PROCEDURE(initialize_id_int), DEFERRED:: initialize_id
    !# Initialize the |id|; for example, set up the lattices around the
    !  stars for the |bns| produced with |fuka|.

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
    !  [[idbase:return_spatial_extent]] and [[idbase:get_total_spatial_extent]]
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
    !# Sets [[idbase:cold_system]], the `LOGICAL` variable that specifies if
    !  the system is cold (no thermal component)


    PROCEDURE, NON_OVERRIDABLE:: get_cold_system
    !# Returns [[idbase:cold_system]], the `LOGICAL` variable that specifies if
    !  the system is cold (no thermal component)


    PROCEDURE, NON_OVERRIDABLE:: set_estimate_length_scale
    !# Sets [[idbase:estimate_length_scale]], the `LOGICAL` variable that
    !  specifies if a typical length scale, equal to the ratio of a field over
    !  its gradient, should be computed


    PROCEDURE, NON_OVERRIDABLE:: get_estimate_length_scale
    !# Returns [[idbase:estimate_length_scale]], the `LOGICAL` variable that
    !  specifies if a typical length scale, equal to the ratio of a field over
    !  its gradient, should be computed


    PROCEDURE, NON_OVERRIDABLE:: check_i_matter
    !# Checks that the given index is between 1 and [[idbase:n_matter]],
    !  included. If not, it stops the execution of the program.


    PROCEDURE:: integrate_baryon_mass_density
    !# Integrates the baryon mass density over a matter object, using spherical
    !  coordinates, and computes its radial profile inside the star


    PROCEDURE:: estimate_lengthscale_field
    !# Estimate typical length scales, one per each matter object, by
    !  computing \(\dfrac{f}{\partial f}\), where \(f\) is a field given
    !  as input, and \(\partial\) represent a derivative of it.
    !  Presently, the derivatives are computed separately along each spatial
    !  dimension, as 1D derivatives.


  END TYPE idbase


  ABSTRACT INTERFACE


    SUBROUTINE derived_type_constructor_int( derived_type, filename, &
                                             eos_filenames )
    !# Construct the DERIVED TYPE that extends [[idbase]]

      IMPORT:: idbase
      CLASS(idbase),    INTENT(OUT)          :: derived_type
      !! Object of DERIVED TYPE that extends [[idbase], to be constructed
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: filename
      !! Name of the file containing the |id|
      CHARACTER(LEN=*), DIMENSION(:), INTENT(IN), OPTIONAL :: eos_filenames
      !# Array of strings containing the names of the files containing the |eos|
      !  to be used for each matter object. If not PRESENT, information from
      !  the file `filename` is used

    END SUBROUTINE derived_type_constructor_int


    FUNCTION read_double_at_pos( this, x, y, z ) RESULT( res )
    !# INTERFACE for a PROCEDURE that returns a DOUBLE PRECISION at a given
    !  position

      IMPORT:: idbase
      CLASS(idbase),    INTENT(IN)       :: this
      !! Object of class [[idbase]] which this PROCEDURE is a member of
      DOUBLE PRECISION, INTENT(IN), VALUE:: x
      !! \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN), VALUE:: y
      !! \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN), VALUE:: z
      !! \(z\) coordinate of the desired point
      DOUBLE PRECISION:: res
      !! Real number at \((x,y,z)\)

    END FUNCTION read_double_at_pos


    FUNCTION read_integer_at_pos( this, x, y, z ) RESULT( res )
    !# INTERFACE for a PROCEDURE that returns an INTEGER at a given position

      IMPORT:: idbase
      CLASS(idbase),    INTENT(IN)       :: this
      !! Object of class [[idbase]] which this PROCEDURE is a member of
      DOUBLE PRECISION, INTENT(IN), VALUE:: x
      !! \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN), VALUE:: y
      !! \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN), VALUE:: z
      !! \(z\) coordinate of the desired point
      INTEGER:: res
      !! Integer at \((x,y,z)\)

    END FUNCTION read_integer_at_pos


    FUNCTION read_logical_at_pos( this, x, y, z ) RESULT( res )
    !# INTERFACE for a PROCEDURE that returns a LOGICAL at a given position

      IMPORT:: idbase
      CLASS(idbase),    INTENT(IN)       :: this
      !! Object of class [[idbase]] which this PROCEDURE is a member of
      DOUBLE PRECISION, INTENT(IN), VALUE:: x
      !! \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN), VALUE:: y
      !! \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN), VALUE:: z
      !! \(z\) coordinate of the desired point
      LOGICAL:: res
      !! Logical at \((x,y,z)\)

    END FUNCTION read_logical_at_pos


    FUNCTION return_double_at_object( this, i_matter ) RESULT( res )
    !# INTERFACE for a PROCEDURE that returns a DOUBLE PRECISION

      IMPORT:: idbase
      CLASS(idbase), INTENT(IN):: this
      INTEGER,       INTENT(IN):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION:: res
      !! Real number. Parameter of the `i_matter`-th matter object

    END FUNCTION return_double_at_object


    FUNCTION return_double( this ) RESULT( res )
    !# INTERFACE for a PROCEDURE that returns a DOUBLE PRECISION

      IMPORT:: idbase
      CLASS(idbase), INTENT(IN):: this
      DOUBLE PRECISION:: res
      !! Real number

    END FUNCTION return_double


    FUNCTION return_position( this, i_matter ) RESULT( res )
    !# INTERFACE for a PROCEDURE that returns a DOUBLE PRECISION

      IMPORT:: idbase
      CLASS(idbase), INTENT(IN):: this
      INTEGER,       INTENT(IN):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION, DIMENSION(3):: res
      !# Centers of the matter objects. The first index runs over the matter
      !  objects, the second index over \((x,y,z)\).

    END FUNCTION return_position


    FUNCTION return_integer_parameter( this, i_matter ) RESULT( res )
    !# INTERFACE for a PROCEDURE that returns an INTEGER

      IMPORT:: idbase
      CLASS(idbase), INTENT(IN):: this
      INTEGER,       INTENT(IN):: i_matter
      !! Index of the matter object whose parameter is to return
      INTEGER:: res
      !! Real number. Parameter of the `i_matter`-th matter object

    END FUNCTION return_integer_parameter


    SUBROUTINE return_eos_parameters_int( this, i_matter, eos_params )
    !# INTERFACE for a PROCEDURE that returns an array containing the
    !  parametersf the |eos| for the matter objects

      IMPORT:: idbase
      CLASS(idbase), INTENT(IN):: this
      INTEGER,       INTENT(IN):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT):: eos_params
      !# Array containing the parameters of the |eos| for the `i_matter`-th
      !  matter object

    END SUBROUTINE return_eos_parameters_int


    FUNCTION return_string_parameter( this, i_matter ) RESULT( string )
    !# INTERFACE for a PROCEDURE that returns a CHARACTER(LEN=:)

      IMPORT:: idbase
      CLASS(idbase), INTENT(IN):: this
      !! [[idbase]] object which this PROCEDURE is a member of
      INTEGER, INTENT(IN):: i_matter
      !! Index of the matter object whose string is to return
      CHARACTER(LEN=:), ALLOCATABLE:: string

    END FUNCTION return_string_parameter


    SUBROUTINE read_id_mass_b_int( this, x, y, z, &
                                   g, &
                                   baryon_density, &
                                   gamma_euler )
      !# INTERFACE or the SUBROUTINE reading the hydro |id| needed to compute
      !  the baryon mass

      IMPORT:: idbase
      CLASS(idbase),    INTENT(INOUT):: this
      !! Object of class [[idbase]] which this PROCEDURE is a member of
      DOUBLE PRECISION, INTENT(IN)    :: x
      DOUBLE PRECISION, INTENT(IN)    :: y
      DOUBLE PRECISION, INTENT(IN)    :: z
      DOUBLE PRECISION, DIMENSION(6), INTENT(OUT):: g
      DOUBLE PRECISION, INTENT(OUT):: baryon_density
      DOUBLE PRECISION, INTENT(OUT):: gamma_euler

    END SUBROUTINE read_id_mass_b_int


    SUBROUTINE read_id_full_int( this, n, x, y, z, &
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
     !# INTERFACE or the SUBROUTINE reading the full |id|

      IMPORT:: idbase
      !> [[idbase]] object which this PROCEDURE is a member of
      CLASS(idbase),                  INTENT(INOUT):: this
      INTEGER,                        INTENT(IN)    :: n
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN)    :: x
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN)    :: y
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN)    :: z
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: lapse
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: shift_x
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: shift_y
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: shift_z
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: g_xx
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: g_xy
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: g_xz
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: g_yy
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: g_yz
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: g_zz
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: k_xx
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: k_xy
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: k_xz
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: k_yy
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: k_yz
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: k_zz
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: baryon_density
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: energy_density
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: specific_energy
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: pressure
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: u_euler_x
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: u_euler_y
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: u_euler_z

    END SUBROUTINE read_id_full_int


    SUBROUTINE read_id_spacetime_int( this, nx, ny, nz, &
                                              pos, &
                                              lapse, &
                                              shift, &
                                              g, &
                                              ek )
     !# INTERFACE or the SUBROUTINE reading the spacetime |id|

      IMPORT:: idbase
      !> [[idbase]] object which this PROCEDURE is a member of
      CLASS(idbase),                        INTENT(INOUT):: this
      INTEGER,                              INTENT(IN)    :: nx
      INTEGER,                              INTENT(IN)    :: ny
      INTEGER,                              INTENT(IN)    :: nz
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(IN)    :: pos
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT(INOUT):: lapse
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(INOUT):: shift
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(INOUT):: g
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(INOUT):: ek

    END SUBROUTINE read_id_spacetime_int


    SUBROUTINE read_id_hydro_int( this, nx, ny, nz, &
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
      CLASS(idbase),                        INTENT(INOUT):: this
      INTEGER,                              INTENT(IN)    :: nx
      INTEGER,                              INTENT(IN)    :: ny
      INTEGER,                              INTENT(IN)    :: nz
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(IN)    :: pos
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT(INOUT):: baryon_density
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT(INOUT):: energy_density
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT(INOUT):: specific_energy
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT(INOUT):: pressure
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(INOUT):: u_euler

    END SUBROUTINE read_id_hydro_int


    SUBROUTINE read_id_particles_int( this, n, x, y, z, &
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
      CLASS(idbase),                  INTENT(INOUT):: this
      INTEGER,                        INTENT(IN)    :: n
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN)    :: x
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN)    :: y
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN)    :: z
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: lapse
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: shift_x
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: shift_y
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: shift_z
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: g_xx
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: g_xy
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: g_xz
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: g_yy
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: g_yz
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: g_zz
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: baryon_density
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: energy_density
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: specific_energy
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: pressure
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: u_euler_x
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: u_euler_y
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: u_euler_z

    END SUBROUTINE read_id_particles_int


    SUBROUTINE read_id_k_int( this, n, x, y, z,&
                              k_xx, k_xy, k_xz, &
                              k_yy, k_yz, k_zz )
    !# INTERFACE or the SUBROUTINE reading the components of the extrinsic
    !  curvature

      IMPORT:: idbase
      !> [[idbase]] object which this PROCEDURE is a member of
      CLASS(idbase),                  INTENT(INOUT):: this
      INTEGER,                        INTENT(IN)    :: n
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN)    :: x
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN)    :: y
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN)    :: z
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: k_xx
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: k_xy
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: k_xz
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: k_yy
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: k_yz
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: k_zz

    END SUBROUTINE read_id_k_int


    FUNCTION return_spatial_extent_int( this, i_matter ) RESULT( box )
    !# INTERFACE to the SUBROUTINE that detects the spatial extent of the
    !  matter objects, and returns a 6-dimensional array
    !  containing the coordinates
    !\(x_{\rm min},x_{\rm max},y_{\rm min},y_{\rm max},z_{\rm min},z_{\rm max}\)
    !  of a box **centered at the center of the object** and containing the
    !  system.

      IMPORT:: idbase
      CLASS(idbase), INTENT(IN):: this
      !! Object of class [[idbase]] which this PROCEDURE is a member of
      INTEGER,       INTENT(IN):: i_matter
      !! Index of the matter object whose string is to return
      DOUBLE PRECISION, DIMENSION(6):: box
      !# 6-dimensional array containing the coordinates
      !  \(x_{\rm min},x_{\rm max},y_{\rm min},y_{\rm max},
      !    z_{\rm min},z_{\rm max}\)
      !  of a box containing the physical system.

    END FUNCTION return_spatial_extent_int


    SUBROUTINE print_summary_int( this, filename )
    !# Prints a summary of the physical properties the system
    !  to the standard output and, optionally, to a formatted file whose name
    !  is given as the optional argument `filename`

      IMPORT:: idbase
      CLASS(idbase),    INTENT(IN):: this
      CHARACTER(LEN=*), INTENT(INOUT), OPTIONAL:: filename
      !! Name of the formatted file to print the summary to

    END SUBROUTINE print_summary_int


    SUBROUTINE initialize_id_int( this, flag, switch )
    !# Initialize the |id|; for example, set up the lattices around the
    !  stars for the |bns| produced with |fuka|.

      IMPORT:: idbase
      CLASS(idbase), INTENT(INOUT)       :: this
      INTEGER,       INTENT(IN)          :: flag
      !! Identifies what kind of initialization has to be done
      LOGICAL,       INTENT(IN), OPTIONAL:: switch
      !! If `.TRUE.`, switch to a different initialization

    END SUBROUTINE initialize_id_int


    SUBROUTINE finalize_sph_id_int &
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

    END SUBROUTINE finalize_sph_id_int


  END INTERFACE


  INTERFACE


    MODULE SUBROUTINE sanity_check( derived_type )
    !# Check that the DERIVED TYPE that extends [[idbase]] is constructed
    !  consistently

      CLASS(idbase), INTENT(IN):: derived_type
      !! Object of DERIVED TYPE that extends [[idbase]]

    END SUBROUTINE sanity_check


    MODULE SUBROUTINE initialize( derived_type, filename, eos_filenames )
    !# Initialize the DERIVED TYPE that extends [[idbase]]

      CLASS(idbase),    INTENT(OUT)          :: derived_type
      !! Object of DERIVED TYPE that extends [[idbase], to be constructed
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: filename
      !! Name of the file containing the |id|
      CHARACTER(LEN=*), DIMENSION(:), INTENT(IN), OPTIONAL :: eos_filenames
      !# Array of strings containing the names of the files containing the |eos|
      !  to be used for each matter object. If not PRESENT, information from
      !  the file `filename` is used

    END SUBROUTINE initialize


    MODULE SUBROUTINE integrate_baryon_mass_density &
      ( this, center, radius, central_density, dr, dth, dphi, &
        mass, mass_profile, mass_profile_idx, radii, surf )
    !# INTERFACE to the SUBROUTINE integrating the baryon mass density to
    !  compute the radial mass profile of a single star.

      !> Object of class [[idbase]] which this PROCEDURE is a member of
      CLASS(idbase),                  INTENT(INOUT):: this
      !> Center of the star
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN)   :: center
      !> Central density of the star
      DOUBLE PRECISION,               INTENT(IN)   :: radius
      !> Integration steps
      DOUBLE PRECISION,               INTENT(IN)   :: central_density
      !> Radius of the star
      DOUBLE PRECISION,               INTENT(IN)   :: dr, dth, dphi
      !> Integrated mass of the star
      DOUBLE PRECISION,               INTENT(INOUT):: mass
      !> Array storing the radial mass profile of the star
      DOUBLE PRECISION, DIMENSION(3,0:NINT(radius/dr)), INTENT(OUT):: &
                                           mass_profile
      !& Array to store the indices for array mass_profile, sorted so that
      !  `mass_profile[mass_profile_idx]` is in increasing order
      INTEGER, DIMENSION(0:NINT(radius/dr)), INTENT(OUT):: mass_profile_idx
      DOUBLE PRECISION, DIMENSION(2), INTENT(IN), OPTIONAL:: radii
      !> Surface of the matter object
      TYPE(surface),                  INTENT(IN), OPTIONAL:: surf

    END SUBROUTINE integrate_baryon_mass_density


    PURE MODULE FUNCTION get_n_matter( this )
    !# Returns [[idbase:n_matter]], the number of matter objects in the
    !  physical system

      CLASS(idbase), INTENT(IN):: this
      INTEGER:: get_n_matter
      !! [[idbase:n_matter]], the number of matter objects in the
      !  physical system

    END FUNCTION get_n_matter


    PURE MODULE SUBROUTINE set_n_matter( this, value )
    !# Sets [[idbase:n_matter]], the number of matter objects in the
    !  physical system, to the given value

      CLASS(idbase), INTENT(INOUT):: this
      INTEGER, INTENT(IN):: value
      !! Value to set [[idbase:n_matter]] to

    END SUBROUTINE set_n_matter


    PURE MODULE FUNCTION get_cold_system( this )
    !# Returns [[idbase:cold_system]], the `LOGICAL` variable at specifies if
    !  the system is cold (no thermal component)

      CLASS(idbase), INTENT(IN):: this
      LOGICAL:: get_cold_system
      !! [[idbase:cold_system]]

    END FUNCTION get_cold_system


    MODULE SUBROUTINE set_cold_system( this, value )
    !# Sets [[idbase:cold_system]], the `LOGICAL` variable at specifies if
    !  the system is cold (no thermal component)

      CLASS(idbase), INTENT(INOUT):: this
      LOGICAL, INTENT(IN):: value
      !! Value to set [[idbase:cold_system]] to

    END SUBROUTINE set_cold_system


    PURE MODULE FUNCTION get_estimate_length_scale( this )
    !# Returns [[idbase:estimate_length_scale]], the `LOGICAL` variable
    !  that specifies if a typical length scale, equal to the ratio
    !  of a field over its gradient, should be computed

      CLASS(idbase), INTENT(IN):: this
      LOGICAL:: get_estimate_length_scale
      !! [[idbase:estimate_length_scale]]

    END FUNCTION get_estimate_length_scale


    MODULE SUBROUTINE set_estimate_length_scale( this, value )
    !# Sets [[idbase:estimate_length_scale]], the `LOGICAL` variable
    !  that specifies if a typical length scale, equal to the ratio
    !  of a field over its gradient, should be computed

      CLASS(idbase), INTENT(INOUT):: this
      LOGICAL, INTENT(IN):: value
      !! Value to set [[idbase:cold_system]] to

    END SUBROUTINE set_estimate_length_scale


    MODULE SUBROUTINE check_i_matter( this, i_matter )
    !# Checks that the given index `i_matter` is between 1 and
    !  [[idbase:n_matter]], included. If not, it stops the execution of the
    !  program.

      CLASS(idbase), INTENT(IN):: this
      INTEGER, INTENT(IN):: i_matter
      !! Value to be checked

    END SUBROUTINE check_i_matter


    MODULE FUNCTION get_total_spatial_extent( this ) RESULT( box )
    !# INTERFACE to the SUBROUTINE that detects the spatial extent of the
    !  physical system considered, and returns a 6-dimensional array
    !  containing the coordinates
    !\(x_{\rm min},x_{\rm max},y_{\rm min},y_{\rm max},z_{\rm min},z_{\rm max}\)
    !  of a box **centered at the center of the object** and containing the
    !  system.

      CLASS(idbase), INTENT(IN)   :: this
      !! Object of class [[idbase]] which this PROCEDURE is a member of
      DOUBLE PRECISION, DIMENSION(6):: box
      !# 6-dimensional array containing the coordinates
      !  \(x_{\rm min},x_{\rm max},y_{\rm min},y_{\rm max},
      !    z_{\rm min},z_{\rm max}\)
      !  of a box containing the physical system.

    END FUNCTION get_total_spatial_extent


    PURE MODULE FUNCTION get_one_lapse( this )
    !# Returns [[idbase:n_matter]], the number of matter objects in the
    !  physical system

      CLASS(idbase), INTENT(IN):: this
      LOGICAL:: get_one_lapse
      !! [[idbase:n_matter]], the number of matter objects in the
      !  physical system

    END FUNCTION get_one_lapse


    PURE MODULE SUBROUTINE set_one_lapse( this, logic )
    !# Sets [[idbase:n_matter]], the number of matter objects in the
    !  physical system, to the given value

      CLASS(idbase), INTENT(INOUT):: this
      LOGICAL, INTENT(IN):: logic
      !! Value to set [[idbase:n_matter]] to

    END SUBROUTINE set_one_lapse


    MODULE PURE FUNCTION get_zero_shift( this )
    !# Returns [[idbase:n_matter]], the number of matter objects in the
    !  physical system

      CLASS(idbase), INTENT(IN):: this
      LOGICAL:: get_zero_shift
      !! [[idbase:n_matter]], the number of matter objects in the
      !  physical system

    END FUNCTION get_zero_shift


    PURE MODULE SUBROUTINE set_zero_shift( this, logic )
    !# Sets [[idbase:n_matter]], the number of matter objects in the
    !  physical system, to the given value

      CLASS(idbase), INTENT(INOUT):: this
      LOGICAL, INTENT(IN):: logic
      !! Value to set [[idbase:n_matter]] to

    END SUBROUTINE set_zero_shift


    MODULE FUNCTION estimate_lengthscale_field( this, get_field, n_mat ) &
      RESULT( scales )
    !# Estimate typical length scales, one per each matter object, by
    !  computing \(\dfrac{f}{\partial f}\), where \(f\) is a field given
    !  as input, and \(\partial\) represent a derivative of it.
    !  Presently, the derivatives are computed separately along each spatial
    !  dimension, as 1D derivatives.

      CLASS(idbase), INTENT(INOUT):: this
      INTERFACE
        FUNCTION get_field( x, y, z ) RESULT( val )
          !! Returns the value of a field at the desired point
          DOUBLE PRECISION, INTENT(IN):: x
          !! \(x\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: y
          !! \(y\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: z
          !! \(z\) coordinate of the desired point
          DOUBLE PRECISION:: val
          !! Value of the field at \((x,y,z)\)
        END FUNCTION get_field
      END INTERFACE

      INTEGER, INTENT(IN):: n_mat
      ! Number of matter objects in the physical ystem
      DOUBLE PRECISION, DIMENSION(n_mat):: scales
      !# Array of the minimum \(\dfrac{f}{\partial f}\) over the lattices that
      !  surround each matter object

    END FUNCTION estimate_lengthscale_field


  END INTERFACE


END MODULE id_base
