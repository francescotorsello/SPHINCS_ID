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
    TYPE(C_PTR):: diffstar_ptr

    !> Logical variables to set the geodesic gauge (lapse=1, shift=0)
    LOGICAL, PUBLIC:: one_lapse, zero_shift

    !> Timer that times the construction of the |lorene| Etdiffrot object
    TYPE(timer), PUBLIC:: drs_construction_timer


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

    PROCEDURE:: return_eos_parameters => get_eos_parameters

    PROCEDURE, PUBLIC:: get_eos_loreneid
    !! Returns [[diffstarlorene:eos_loreneid]]

    PROCEDURE, PUBLIC:: get_diffstar_identifier
    !! Returns [[diffstarlorene:diffstar_identifier]]]

    !PROCEDURE, PUBLIC:: get_diffstar_ptr

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


    MODULE FUNCTION get_diffstar_identifier( THIS )

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_diffstar_identifier

    END FUNCTION get_diffstar_identifier


    MODULE FUNCTION get_eos_loreneid( THIS )

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene), INTENT( IN ):: THIS
      ! Result
      INTEGER:: get_eos_loreneid

    END FUNCTION get_eos_loreneid


    MODULE SUBROUTINE get_eos_parameters( THIS, i_matter, eos_params )

      !> [[diffstarlorene]] object which this PROCEDURE is a member of
      CLASS(diffstarlorene), INTENT( IN OUT ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT):: eos_params
      !# Array containing the parameters of the |eos| for the DRS

    END SUBROUTINE get_eos_parameters


    !MODULE FUNCTION get_diffstar_ptr( THIS )
    !
    !  ! Argument
    !  CLASS(diffstarlorene), INTENT( IN ):: THIS
    !  ! Result
    !  TYPE(C_PTR):: get_diffstar_ptr
    !
    !END FUNCTION get_diffstar_ptr


  END INTERFACE


  PRIVATE:: construct_etdiffrot, get_diffstar_full, get_diffstar_spacetime, &
            get_diffstar_particles, get_diffstar_mass_b, &
            get_diffstar_hydro, get_diffstar_mass_density, &
            get_diffstar_spatial_metric, negative_hydro, get_diffstar_params, &
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
                                  g_diag, &
                                  k_xx, k_xy, k_xz, &
                                  k_yy, k_yz, k_zz, &
                                  baryon_density, &
                                  energy_density, &
                                  specific_energy, &
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
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_x
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_y
      REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_z

    END SUBROUTINE get_diffstar_full


    SUBROUTINE get_diffstar_spacetime( optr, &
                                       x, y, z, &
                                       lapse, &
                                       shift_x, shift_y, shift_z, &
                                       g_diag, &
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
      REAL(C_DOUBLE), INTENT(OUT)       :: g_diag
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
                                       g_diag, &
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
      REAL(C_DOUBLE), INTENT(OUT)       :: g_diag
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
                                    g_diag, &
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
      REAL(C_DOUBLE), INTENT(OUT)       :: g_diag
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
      BIND(C, NAME= "get_g_diag")

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


    FUNCTION negative_hydro( optr, x, y, z ) RESULT( res ) &
      BIND(C, NAME= "is_hydro_negative")

      !************************************************
      !
      !# Interface to the |lorene| method of class
      !  |etdiffrot| with the same name, that returns 1
      !  if the energy density is nonpositive,
      !  or if the specific energy is nonpositive,
      !  or if the pressure is nonpositive,
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
      !  are negative, 0 otherwise
      INTEGER(C_INT) :: res

    END FUNCTION negative_hydro


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
