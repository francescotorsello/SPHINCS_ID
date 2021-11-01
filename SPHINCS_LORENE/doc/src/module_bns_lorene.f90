! File:         module_bns_lorene.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

MODULE bns_lorene

  !***********************************************************
  !
  !#  This module contains the definition of TYPE bnslorene,
  !   and the SUBROUTINES that bind to the methods
  !   of |lorene|'s class |binns|, defined in
  !   Lorene/Export/BinNS
  !
  !   [|lorene| official repository](https://lorene.obspm.fr/index.html){:target="_blank"}
  !
  !***********************************************************


  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_CHAR, C_NULL_CHAR, &
                                         C_PTR, C_NULL_PTR, C_ASSOCIATED
  USE bns_base,                    ONLY: bnsbase
  USE id_base,                     ONLY: idbase
  USE utility,                     ONLY: itr, ios, err_msg, test_status, &
                                         perc, creturn, compute_g4, &
                                         determinant_sym4x4_grid, show_progress
  USE timing,                      ONLY: timer


  IMPLICIT NONE


  !*******************************************************
  !                                                      *
  !            Definition of TYPE bnslorene              *
  !                                                      *
  !   This class imports and stores the |lorene| BNS ID  *
  !                                                      *
  !*******************************************************

  TYPE, EXTENDS(bnsbase):: bnslorene
  !! TYPE representing a binary system of neutron stars (bns)


    PRIVATE


    !> Identifier of the bnslorene object
    INTEGER:: bns_identifier= 0
    !> |lorene| identifiers for the EoS
    INTEGER:: eos1_loreneid, eos2_loreneid

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

    !& C pointer to the |lorene|'s |binns| object
    ! N.B. This variable is global. The pointer to the second |lorene| |binns|
    !      object will overwrite the first one, and so on.
    !      This variable stores the pointer to the last defined |lorene| |binns|
    !      object. That's why it is not freed in the destructor of a bns object.
    !      Presently, it has to be freed by the user at the end of the PROGRAM.
    !      See the last part of the PROGRAM in setup_lorene_id.f90, for example.
    TYPE(C_PTR):: bns_ptr

    !> Logical variables to set the geodesic gauge (lapse=1, shift=0)
    LOGICAL, PUBLIC:: one_lapse, zero_shift

    !> Timer that times the construction of the |lorene| |binns| object
    TYPE(timer), PUBLIC:: binary_construction_timer


    CONTAINS


    !-------------------!
    !--  SUBROUTINES  --!
    !-------------------!

    PROCEDURE:: construct_binary
    !! Constructs the |lorene| |binns| object

    PROCEDURE:: destruct_binary
    !! Destructs the |lorene| |binns| object

    PROCEDURE:: allocate_lorene_id_memory
    !! Allocates memory for the [[bnslorene]] member arrays

    PROCEDURE:: deallocate_lorene_id_memory
    !! Deallocates memory for the [[bnslorene]] member arrays

    PROCEDURE:: import_id_params
    !! Imports the parameters of the BNS from |lorene|

    !PROCEDURE:: integrate_field_on_star => integrate_baryon_mass_density
    !# Integrates the |lorene| baryon mass density and computes the
    !  radial mass profile

    PROCEDURE, PUBLIC:: print_id_params
    !! Prints the parameters of the BNS to the standard output

    PROCEDURE:: import_id_int
    !! Stores the ID in the [[bnslorene]] member arrays

    PROCEDURE:: read_id_full      => import_id_full
    PROCEDURE:: read_id_spacetime => import_id_spacetime
    PROCEDURE:: read_id_particles => import_id_particles
    PROCEDURE:: read_id_hydro     => import_id_hydro
    PROCEDURE:: read_id_mass_b    => import_id_mass_b
    PROCEDURE:: read_id_k         => import_id_k

    !-----------------!
    !--  FUNCTIONS  --!
    !-----------------!

    !> Returns the |lorene|'s mass density at the given point
    PROCEDURE:: read_mass_density => import_mass_density

    !> Returns the |lorene|'s conformally flat spatial ADM metric
    PROCEDURE:: import_spatial_metric

    !& Returns 1 if the energy density or the specific energy or the pressure
    !  are negative
    PROCEDURE:: test_position => is_hydro_negative

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

    FINAL:: destruct_bnslorene
    !! Finalizer (Destructor) of a [[bnslorene]] object

  END TYPE bnslorene

  !
  !-- Interface of the TYPE bnslorene (i.e., declaration of the constructor)
  !-- (see https://dannyvanpoucke.be/oop-fortran-tut4-en/)
  !
  INTERFACE bnslorene
  !! Interface of TYPE [[bnslorene]]

    MODULE PROCEDURE:: construct_bnslorene
    !! Constructs a [[bnslorene]] object

  END INTERFACE bnslorene

  !
  !-- Interfaces of the constructor and destructor of the TYPE bnslorene
  !
  INTERFACE

   ! MODULE SUBROUTINE construct_bnslorene2( this )
   ! !! Constructs a [[bnslorene]] object
   !
   !   TYPE(bnslorene), INTENT( IN OUT ):: this
   !   !! Constructed [[bnslorene]] object
   !
   ! END SUBROUTINE construct_bnslorene2

    MODULE FUNCTION construct_bnslorene( resu_file ) RESULT( bns_obj )
    !! Constructs a [[bnslorene]] object

      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: resu_file
      !! |lorene| binary file containing the spectral BNS ID
      TYPE(bnslorene):: bns_obj
      !! Constructed [[bnslorene]] object

    END FUNCTION construct_bnslorene

    MODULE SUBROUTINE destruct_bnslorene( THIS )
    !! Destruct a [[bnslorene]] object

      TYPE(bnslorene), INTENT( IN OUT ):: THIS
      !! [[bnslorene]] object to be destructed

    END SUBROUTINE destruct_bnslorene

  END INTERFACE

  !
  !-- Interfaces of the methods of the TYPE bnslorene
  !-- Their implementations are in submodule_bnslorene_methods.f90
  !
  INTERFACE


    !
    !-- SUBROUTINES
    !
    MODULE SUBROUTINE construct_binary( THIS, resu_file )
    !! Interface of the subroutine that constructs the |lorene| |binns| object

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),                     INTENT( IN OUT )      :: THIS
      !> |lorene| binary file containing the spectral BNS ID
      CHARACTER(KIND= C_CHAR, LEN=*), INTENT( IN ), OPTIONAL:: resu_file

    END SUBROUTINE construct_binary


    MODULE SUBROUTINE destruct_binary( THIS )
    !! Destructs a |lorene| |binns| object

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene), INTENT( IN OUT ):: THIS

    END SUBROUTINE destruct_binary


    MODULE SUBROUTINE allocate_lorene_id_memory( THIS, d )
    !! Allocates allocatable arrays member of a [[bnslorene]] object

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene), INTENT( IN OUT ):: THIS
      !> Dimension of the arrays
      INTEGER,    INTENT( IN )    :: d

    END SUBROUTINE allocate_lorene_id_memory


    MODULE SUBROUTINE deallocate_lorene_id_memory( THIS )
    !! Deallocates allocatable arrays member of a [[bnslorene]] object

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene), INTENT( IN OUT ):: THIS

    END SUBROUTINE deallocate_lorene_id_memory


    MODULE SUBROUTINE import_id_params( THIS )
    !! Imports the BNS parameters from |lorene|

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene), INTENT( IN OUT ):: THIS

    END SUBROUTINE import_id_params


    MODULE SUBROUTINE print_id_params( THIS )
    !! Prints the BNS parameters to the standard output

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene), INTENT( IN OUT ):: THIS

    END SUBROUTINE print_id_params


  !  MODULE SUBROUTINE integrate_baryon_mass_density( THIS, center, radius, &
  !                                                   central_density, &
  !                                                   dr, dth, dphi, &
  !                                                   mass, mass_profile, &
  !                                                   mass_profile_idx )
  !  !# Integrates the |lorene| baryon mass density to compute the radial mass
  !  !  profile. TODO: Improve integration algorithm.
  !
  !    !> [[bnslorene]] object which this PROCEDURE is a member of
  !    CLASS(bnslorene), INTENT( IN OUT )      :: THIS
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


    MODULE SUBROUTINE import_id_int( THIS, n, x, y, z )
    !! Stores the ID in the [[bnslorene]] member arrays

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),                     INTENT( IN OUT ):: THIS
      INTEGER, INTENT( IN ):: n
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: x
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: y
      DOUBLE PRECISION, DIMENSION(:), INTENT( IN )    :: z

    END SUBROUTINE import_id_int

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
    !# Stores the ID in non [[bnslorene]]-member arrays with the same shape as the
    !  [[bnslorene]] member arrays

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),                     INTENT( IN OUT ):: THIS
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

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),                           INTENT( IN OUT ):: THIS
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

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),                           INTENT( IN OUT ):: THIS
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

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),                     INTENT( IN OUT ):: THIS
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

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),       INTENT( IN OUT ):: THIS
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

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),                     INTENT( IN OUT ):: THIS
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

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),     INTENT( IN )         :: THIS
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

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),     INTENT( IN )       :: THIS
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

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),     INTENT( IN )       :: THIS
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
    !! Returns the [[bnslorene]] member arrays named field

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),          INTENT( IN )             :: THIS
      !> Name of the desired [[bnslorene]] member array
      CHARACTER( LEN= : ), INTENT( IN ), ALLOCATABLE:: field
      !> Desired [[bnslorene]] member array
      DOUBLE PRECISION, DIMENSION(:),    ALLOCATABLE:: field_array

    END FUNCTION get_field_array


    MODULE FUNCTION get_field_value( THIS, field, n ) RESULT( field_value )
    !! Returns the component n of the [[bnslorene]] member arrays named field

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene),          INTENT( IN )             :: THIS
      !> Name of the desired [[bnslorene]] member array
      CHARACTER( LEN= : ), INTENT( IN ), ALLOCATABLE:: field
      !> Component of the desired [[bnslorene]] member array
      INTEGER,             INTENT( IN )             :: n
      !> Component n of the desired [[bnslorene]] member array
      DOUBLE PRECISION                              :: field_value

    END FUNCTION get_field_value


    MODULE FUNCTION get_bns_identifier( THIS )

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_bns_identifier

    END FUNCTION get_bns_identifier


    MODULE FUNCTION get_eos1_loreneid( THIS )

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene), INTENT( IN ):: THIS
      ! Result
      INTEGER:: get_eos1_loreneid

    END FUNCTION get_eos1_loreneid


    MODULE FUNCTION get_eos2_loreneid( THIS )

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene), INTENT( IN ):: THIS
      ! Result
      INTEGER:: get_eos2_loreneid

    END FUNCTION get_eos2_loreneid


    MODULE SUBROUTINE get_eos_parameters( THIS, i_matter, eos_params )

      !> [[bnslorene]] object which this PROCEDURE is a member of
      CLASS(bnslorene), INTENT( IN ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT):: eos_params
      !# Array containingthe parameters of the |eos| for the `i_matter`-th
      !  matter object

    END SUBROUTINE get_eos_parameters


    !MODULE FUNCTION get_bns_ptr( THIS )
    !
    !  ! Argument
    !  CLASS(bnslorene), INTENT( IN ):: THIS
    !  ! Result
    !  TYPE(C_PTR):: get_bns_ptr
    !
    !END FUNCTION get_bns_ptr


  END INTERFACE


  !------------------------------------------------------------------!
  !--  PRIVATE interfaces to the methods of |lorene|'s class |binns|  --!
  !------------------------------------------------------------------!


  PRIVATE:: construct_bin_ns, get_lorene_id, get_lorene_id_spacetime, &
            get_lorene_id_particles, get_lorene_id_mass_b, &
            get_lorene_id_hydro, get_lorene_id_k, get_lorene_mass_density, &
            get_lorene_spatial_metric, negative_hydro, get_lorene_id_params, &
            destruct_bin_ns


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
                              v_euler_x, v_euler_y, v_euler_z ) &
      BIND(C, NAME= "get_lorene_id")

      !*************************************************
      !
      !# Interface to the |lorene| method of class
      !  |binns| with the same name, that reads the full
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


    FUNCTION negative_hydro( optr, x, y, z ) RESULT( res ) &
      BIND(C, NAME= "negative_hydro")

      !************************************************
      !
      !# Interface to the |lorene| method of class
      !  |binns| with the same name, that returns 1
      !  if the energy density is nonpositive,
      !  or if the specific energy is nonpositive,
      !  or if the pressure is nonpositive,
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
      !  are negative, 0 otherwise
      INTEGER(C_INT) :: res

    END FUNCTION negative_hydro


    SUBROUTINE get_lorene_id_params( optr, &
                                     angular_vel, &
                                     distance, &
                                     distance_com, &
                                     mass1, &
                                     mass2, &
                                     mass_grav1, &
                                     mass_grav2, &
                                     adm_mass, &
                                     angular_momentum, &
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
      REAL(C_DOUBLE), INTENT(OUT)       :: angular_momentum
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
