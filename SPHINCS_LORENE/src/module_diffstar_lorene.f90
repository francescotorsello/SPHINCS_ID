! File:         module_diffstar_lorene.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

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
  USE utility,                     ONLY: itr, ios, err_msg, test_status, &
                                         perc, creturn, compute_g4, &
                                         determinant_sym4x4_grid, show_progress
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
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: baryon_density
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

    !& C pointer to the |lorene|'s Etdiffrot object
    ! N.B. This variable is global. The pointer to the second |lorene| Etdiffrot
    !      object will overwrite the first one, and so on.
    !      This variable stores the pointer to the last defined |lorene| Etdiffrot
    !      object. That's why it is not freed in the destructor of a bns object.
    !      Presently, it has to be freed by the user at the end of the PROGRAM.
    !      See the last part of the PROGRAM in setup_diffstar.f90, for example.
    TYPE(C_PTR):: bns_ptr

    !> Logical variables to set the geodesic gauge (lapse=1, shift=0)
    LOGICAL, PUBLIC:: one_lapse, zero_shift

    !> Timer that times the construction of the |lorene| Etdiffrot object
    TYPE(timer), PUBLIC:: diffstar_construction_timer


    CONTAINS


    !-------------------!
    !--  SUBROUTINES  --!
    !-------------------!

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

    !PROCEDURE:: integrate_field_on_star => integrate_baryon_mass_density
    !# Integrates the |lorene| baryon mass density and computes the
    !  radial mass profile

    PROCEDURE, PUBLIC:: print_diffstar_params
    !! Prints the parameters of the DRS to the standard output

    PROCEDURE:: import_id_int
    !! Stores the ID in the [[diffstarlorene]] member arrays

    PROCEDURE:: read_id_ext       => import_id_full
    PROCEDURE:: read_id_spacetime => import_id_spacetime
    PROCEDURE:: read_id_particles => import_id_particles
    PROCEDURE:: read_id_hydro     => import_id_hydro
    PROCEDURE:: read_id_mass_b    => import_id_mass_b
    PROCEDURE:: read_id_k         => import_id_k

    !-----------------!
    !--  FUNCTIONS  --!
    !-----------------!

    PROCEDURE:: read_mass_density => import_mass_density
    !! Returns the |lorene|'s mass density at the given point

    PROCEDURE:: import_spatial_metric
    !! Returns the |lorene|'s conformally flat spatial ADM metric

    PROCEDURE:: test_position => is_hydro_negative
    !# Returns 1 if the energy density or the specific energy or the pressure
    !  are negative

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

    PROCEDURE, PUBLIC:: get_eos_loreneid
    !! Returns [[diffstarlorene:eos_loreneid]]

    PROCEDURE, PUBLIC:: get_bns_identifier
    !! Returns [[diffstarlorene:diffstar_identifier]]]

    !PROCEDURE, PUBLIC:: get_bns_ptr

    FINAL:: destruct_diffstarlorene
    !! Finalizer (Destructor) of a [[diffstarlorene]] object

  END TYPE diffstarlorene

  !
  !-- Interface of the TYPE diffstarlorene
  !-- (i.e., declaration of the constructor)
  !
  INTERFACE diffstarlorene
  !! Interface of TYPE [[diffstarlorene]]

    MODULE PROCEDURE:: construct_diffstarlorene
    !! Constructs a [[diffstarlorene]] object

  END INTERFACE diffstarlorene

  !
  !-- Interfaces of the constructor and destructor of the TYPE diffstarlorene
  !
  INTERFACE

    MODULE FUNCTION construct_diffstarlorene( resu_file ) RESULT( drs )
    !! Constructs a [[diffstarlorene]] object

      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: resu_file
      !! |lorene| binary file containing the spectral DRS ID
      TYPE(diffstarlorene):: drs
      !! Constructed [[diffstarlorene]] object

    END FUNCTION construct_diffstarlorene

    MODULE SUBROUTINE destruct_diffstarlorene( THIS )
    !! Destruct a [[diffstarlorene]] object

      TYPE(diffstarlorene), INTENT( IN OUT ):: THIS
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
    MODULE SUBROUTINE construct_drs( THIS, resu_file )
    !! Interface of the subroutine that constructs the |lorene| Etdiffrot object

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),                     INTENT( IN OUT )      :: THIS
      !> |lorene| binary file containing the spectral DRS ID
      CHARACTER(KIND= C_CHAR, LEN=*), INTENT( IN ), OPTIONAL:: resu_file

    END SUBROUTINE construct_drs


    MODULE SUBROUTINE destruct_drs( THIS )
    !! Destructs a |lorene| Etdiffrot object

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene), INTENT( IN OUT ):: THIS

    END SUBROUTINE destruct_drs


    MODULE SUBROUTINE allocate_diffstar_memory( THIS, d )
    !! Allocates allocatable arrays member of a [[diffstarlorene]] object

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene), INTENT( IN OUT ):: THIS
      !> Dimension of the arrays
      INTEGER,    INTENT( IN )    :: d

    END SUBROUTINE allocate_diffstar_memory


    MODULE SUBROUTINE deallocate_diffstar_memory( THIS )
    !! Deallocates allocatable arrays member of a [[diffstarlorene]] object

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene), INTENT( IN OUT ):: THIS

    END SUBROUTINE deallocate_diffstar_memory


    MODULE SUBROUTINE import_diffstar_params( THIS )
    !! Imports the DRS parameters from |lorene|

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene), INTENT( IN OUT ):: THIS

    END SUBROUTINE import_diffstar_params


    MODULE SUBROUTINE print_diffstar_params( THIS )
    !! Prints the DRS parameters to the standard output

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene), INTENT( IN OUT ):: THIS

    END SUBROUTINE print_diffstar_params


    MODULE SUBROUTINE import_id_int( THIS, n, x, y, z )
    !! Stores the ID in the [[diffstarlorene]] member arrays

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),                     INTENT( IN OUT ):: THIS
      INTEGER, INTENT( IN ):: n
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: x
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: y
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: z

    END SUBROUTINE import_id_int


    MODULE SUBROUTINE import_id_full( THIS, n, x, y, z,&
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
    !# Stores the ID in non [[diffstarlorene]]-member arrays with the same
    !  shape as the [[diffstarlorene]] member arrays

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),                     INTENT( IN OUT ):: THIS
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

    END SUBROUTINE import_id_full


    MODULE SUBROUTINE import_id_spacetime( THIS, nx, ny, nz, &
                                              pos, &
                                              lapse, &
                                              shift, &
                                              g, &
                                              ek )
    !# Stores the spacetime ID in multi-dimensional arrays needed to compute
    !  the BSSN variables and constraints

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),                           INTENT( IN OUT ):: THIS
      INTEGER,                              INTENT( IN )    :: nx
      INTEGER,                              INTENT( IN )    :: ny
      INTEGER,                              INTENT( IN )    :: nz
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN )    :: pos
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: lapse
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: shift
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: g
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: ek

    END SUBROUTINE import_id_spacetime


    MODULE SUBROUTINE import_id_hydro( THIS, nx, ny, nz, &
                                             pos, &
                                             baryon_density, &
                                             energy_density, &
                                             specific_energy, &
                                             pressure, &
                                             u_euler )
    !# Stores the hydro ID in the arrays needed to compute the constraints
    !  on the refined mesh

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),                           INTENT( IN OUT ):: THIS
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


    MODULE SUBROUTINE import_id_particles( THIS, n, x, y, z, &
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
      CLASS(diffstarlorene),                     INTENT( IN OUT ):: THIS
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


    MODULE SUBROUTINE import_id_mass_b( THIS, x, y, z, &
                                        g_xx, &
                                        baryon_density, &
                                        gamma_euler )
    !! Stores the hydro ID in the arrays needed to compute the baryon mass

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),       INTENT( IN OUT ):: THIS
      DOUBLE PRECISION, INTENT( IN )    :: x
      DOUBLE PRECISION, INTENT( IN )    :: y
      DOUBLE PRECISION, INTENT( IN)     :: z
      DOUBLE PRECISION, INTENT( IN OUT ):: g_xx
      DOUBLE PRECISION, INTENT( IN OUT ):: baryon_density
      DOUBLE PRECISION, INTENT( IN OUT ):: gamma_euler

    END SUBROUTINE import_id_mass_b


    MODULE SUBROUTINE import_id_k( THIS, n, x, y, z,&
                                         k_xx, k_xy, k_xz, &
                                         k_yy, k_yz, k_zz )
   !! Stores the components of the extrinsic curvature in arrays

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),                     INTENT( IN OUT ):: THIS
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
    MODULE FUNCTION import_mass_density( THIS, x, y, z ) RESULT( res )
    !! Returns the |lorene| baryon mass density at a point \((x,y,z)\)

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),     INTENT( IN )         :: THIS
      !> \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: x
      !> \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: y
      !> \(z\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: z
      !> Baryon mass density at \((x,y,z)\)
      DOUBLE PRECISION:: res

    END FUNCTION import_mass_density


    MODULE FUNCTION import_spatial_metric( THIS, x, y, z ) RESULT( res )
    !# Returns the |lorene| conformally flat spatial metric component
    !  \(g_{xx}=g_{yy}=g_{zz}\) at a point \((x,y,z)\)

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),     INTENT( IN )       :: THIS
      !> \(x\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT( IN ), VALUE:: x
      !> \(y\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT( IN ), VALUE:: y
      !> \(z\) coordinate of the desired point
      REAL(C_DOUBLE), INTENT( IN ), VALUE:: z
      !> \(g_{xx}=g_{yy}=g_{zz}\) at \((x,y,z)\)
      REAL(C_DOUBLE):: res

    END FUNCTION import_spatial_metric


    MODULE FUNCTION is_hydro_negative( THIS, x, y, z ) RESULT( res )
    !# Returns 1 if the energy density or the specific energy or the pressure
    !  are negative, 0 otherwise

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),     INTENT( IN )       :: THIS
      !> \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: x
      !> \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: y
      !> \(z\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: z
      !& 1 if the energy density or the specific energy or the pressure
      !  are negative, 0 otherwise
      INTEGER:: res

    END FUNCTION is_hydro_negative


    MODULE FUNCTION get_field_array( THIS, field ) RESULT( field_array )
    !! Returns the [[diffstarlorene]] member arrays named field

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),          INTENT( IN )             :: THIS
      !> Name of the desired [[diffstarlorene]] member array
      CHARACTER( LEN= : ), INTENT( IN ), ALLOCATABLE:: field
      !> Desired [[diffstarlorene]] member array
      DOUBLE PRECISION, DIMENSION(:),    ALLOCATABLE:: field_array

    END FUNCTION get_field_array


    MODULE FUNCTION get_field_value( THIS, field, n ) RESULT( field_value )
    !! Returns the component n of the [[diffstarlorene]] member arrays named field

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene),          INTENT( IN )             :: THIS
      !> Name of the desired [[diffstarlorene]] member array
      CHARACTER( LEN= : ), INTENT( IN ), ALLOCATABLE:: field
      !> Component of the desired [[diffstarlorene]] member array
      INTEGER,             INTENT( IN )             :: n
      !> Component n of the desired [[diffstarlorene]] member array
      DOUBLE PRECISION                              :: field_value

    END FUNCTION get_field_value


    MODULE FUNCTION get_bns_identifier( THIS )

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_bns_identifier

    END FUNCTION get_bns_identifier


    MODULE FUNCTION get_eos_loreneid( THIS )

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene), INTENT( IN ):: THIS
      ! Result
      INTEGER:: get_eos_loreneid

    END FUNCTION get_eos_loreneid


    !MODULE FUNCTION get_bns_ptr( THIS )
    !
    !  ! Argument
    !  CLASS(diffstarlorene), INTENT( IN ):: THIS
    !  ! Result
    !  TYPE(C_PTR):: get_bns_ptr
    !
    !END FUNCTION get_bns_ptr


  END INTERFACE


  !PRIVATE:: construct_etdiffrot, get_diffstar_full, get_diffstar_spacetime, &
  !          get_diffstar_particles, get_diffstar_mass_b, &
  !          get_diffstar_hydro, get_lorene_mass_density, &
  !          get_lorene_spatial_metric, negative_hydro, get_diffstar_params, &
  !          destruct_etdiffrot



END MODULE diffstar_lorene
