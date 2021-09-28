! File:         module_id_base.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

MODULE id_base

  !***********************************************************
  !
  !# This MODULE contains the definition of TYPE idbase,
  !  which is an ABSTRACT TYPE representing any possible
  !  type of initial data (ID) to be set up for SPHINCS_BSSN.
  !  That is, a binary neutron star system, a rotating
  !  star, a binary black hole system, etc.
  !
  !  PROCEDURES and variables shared by all these types
  !  of ID should belong to TYPE idbase, as
  !  they are inherited by its EXTENDED TYPES that
  !  represent more specific types of ID.
  !
  !  FT 24.09.2021
  !
  !***********************************************************


  USE utility,  ONLY: itr, ios, err_msg, test_status, &
                      perc, creturn, compute_g4, &
                      determinant_sym4x4_grid, show_progress


  IMPLICIT NONE


  !**********************************************************
  !                                                         *
  !              Definition of TYPE idbase                  *
  !                                                         *
  !   This ABSTRACT TYPE represents a generic ID for        *
  !   SPHINCS_BSSN (binary neutron star, rotating star...). *
  !                                                         *
  !**********************************************************

  TYPE, ABSTRACT:: idbase
  !# Represents a generic ID for SPHINCS_BSSN (binary neutron star, rotating
  !  star...)


    CONTAINS


    PROCEDURE(read_double_at_pos),    DEFERRED:: read_mass_density
    !# Returns the baryon mass density at the given point

    PROCEDURE(read_integer_at_pos),   DEFERRED:: test_position
    !# Returns 1 if the position has physically acceptable properties,
    !  0 otherwise

    PROCEDURE(read_id_ext_int),       DEFERRED:: read_id_ext
    !# Stores the full ID into arrays
    PROCEDURE(read_id_particles_int), DEFERRED:: read_id_particles
    !! Stores the hydro ID in the arrays needed to compute the SPH ID
    PROCEDURE(read_id_mass_b_int),    DEFERRED:: read_id_mass_b
    !! Stores the hydro ID in the arrays needed to compute the baryon mass
    PROCEDURE(read_id_spacetime_int), DEFERRED:: read_id_spacetime
    !# Stores the spacetime ID in multi-dimensional arrays needed to compute
    !  the BSSN variables and constraints
    PROCEDURE(read_id_hydro_int),     DEFERRED:: read_id_hydro
    !# Stores the hydro ID in the arrays needed to compute the constraints
    !  on the refined mesh
    PROCEDURE(read_id_k_int),         DEFERRED:: read_id_k
    !! Stores the components of the extrinsic curvature in arrays

    PROCEDURE(integrate_field_int),   DEFERRED:: integrate_field_on_star
    !# Integrates a field over the interior of a star in spherical coordinates,
    !  and computes its radial profile inside the star


  END TYPE idbase


  ABSTRACT INTERFACE

    FUNCTION read_double_at_pos( THIS, x, y, z ) RESULT( res )
    !#

      IMPORT:: idbase
      CLASS(idbase), INTENT( IN )          :: THIS
      !! Object of class [[idbase]] which this PROCEDURE is a member of
      DOUBLE PRECISION, INTENT( IN ), VALUE:: x
      !> \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: y
      !> \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: z
      !> \(z\) coordinate of the desired point
      DOUBLE PRECISION:: res
      !> Real number at \((x,y,z)\)

    END FUNCTION read_double_at_pos


    FUNCTION read_integer_at_pos( THIS, x, y, z ) RESULT( res )
    !#

      IMPORT:: idbase
      CLASS(idbase), INTENT( IN )          :: THIS
      !! Object of class [[idbase]] which this PROCEDURE is a member of
      DOUBLE PRECISION, INTENT( IN ), VALUE:: x
      !> \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: y
      !> \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: z
      !> \(z\) coordinate of the desired point
      INTEGER:: res
      !> Integer at \((x,y,z)\)

    END FUNCTION read_integer_at_pos


    SUBROUTINE read_id_mass_b_int( THIS, x, y, z, &
                                   g_xx, &
                                   baryon_density, &
                                   gamma_euler )
      !#

      IMPORT:: idbase
      !> [[idbns]] object which this PROCEDURE is a member of
      CLASS(idbase),    INTENT( IN OUT ):: THIS
      DOUBLE PRECISION, INTENT( IN )    :: x
      DOUBLE PRECISION, INTENT( IN )    :: y
      DOUBLE PRECISION, INTENT( IN)     :: z
      DOUBLE PRECISION, INTENT( IN OUT ):: g_xx
      DOUBLE PRECISION, INTENT( IN OUT ):: baryon_density
      DOUBLE PRECISION, INTENT( IN OUT ):: gamma_euler

    END SUBROUTINE read_id_mass_b_int


    SUBROUTINE read_id_ext_int( THIS, n, x, y, z,&
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
    !# Stores the ID in non [[bns]]-member arrays with the same shape as the
    !  [[bns]] member arrays
      IMPORT:: idbase
      !> [[bns]] object which this PROCEDURE is a member of
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

    END SUBROUTINE read_id_ext_int


    SUBROUTINE read_id_spacetime_int( THIS, nx, ny, nz, &
                                              pos, &
                                              lapse, &
                                              shift, &
                                              g, &
                                              ek )
    !# Stores the spacetime ID in multi-dimensional arrays needed to compute
    !  the BSSN variables and constraints
      IMPORT:: idbase
      !> [[bns]] object which this PROCEDURE is a member of
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
    !# Stores the hydro ID in the arrays needed to compute the constraints
    !  on the refined mesh
      IMPORT:: idbase
      !> [[bns]] object which this PROCEDURE is a member of
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
    !! Stores the hydro ID in the arrays needed to compute the SPH ID
      IMPORT:: idbase
      !> [[bns]] object which this PROCEDURE is a member of
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
    !! Stores the components of the extrinsic curvature in arrays
      IMPORT:: idbase
      !> [[bns]] object which this PROCEDURE is a member of
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


    SUBROUTINE integrate_field_int( THIS, center, radius, &
                                                  central_density, &
                                                  dr, dth, dphi, &
                                                  mass, mass_profile, &
                                                  mass_profile_idx )
    !# Integrates the baryon mass density to compute the radial mass
    !  profile of a single star. TODO: Improve integration algorithm.

      IMPORT:: idbase
      !> Object of class [[idbase]] which this PROCEDURE is a member of
      CLASS(idbase), INTENT( IN OUT )      :: THIS
      !& Array to store the indices for array mass_profile, sorted so that
      !  mass_profile[mass_profile_idx] is in increasing order
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT( IN OUT ):: mass_profile_idx
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
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT( IN OUT ):: &
                                       mass_profile

    END SUBROUTINE integrate_field_int


  END INTERFACE


END MODULE id_base
