! File:         module_bns_id.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

MODULE bns_id
!! This module contains the definition of TYPE bns, and the SUBROUTINES that
!! bind to the methods

  !********************************************************
  !                                                       *
  !   This module contains the definition of TYPE bns,    *
  !   and the SUBROUTINES that bind to the methods        *
  !   of LORENE's class Bin_NS, defined in                *
  !   Lorene/Export/BinNS                                 *
  !                                                       *
  !   LORENE official repository:                         *
  !   https://lorene.obspm.fr/index.html                  *
  !                                                       *
  !   Repository to store the LORENE binary neutron       *
  !   stars' initial data:                                *
  !   https://github.com/francescotorsello/LORENE-BNS-ID  *
  !                                                       *
  !********************************************************


  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, &
                                         C_CHAR, C_NULL_CHAR, &
                                         C_PTR, C_NULL_PTR, C_ASSOCIATED
  USE utility,                     ONLY: itr, ios, err_msg, test_status, &
                                         perc, creturn, compute_g4, &
                                         determinant_sym4x4_grid, show_progress                                        
  USE timing,                      ONLY: timer


  IMPLICIT NONE


  !*******************************************************
  !                                                      *
  !     Definition of TYPE bns (binary neutron star)     *
  !                                                      *
  !   This class imports and stores the LORENE BNS ID    *
  !                                                      *
  !*******************************************************

  TYPE:: bns
  !! TYPE representing a binary sstem of neutron stars (bns)


    PRIVATE


    ! Assign a unique identifies to each bns object created
    INTEGER:: bns_identifier= 0
    ! Identifiers for the EoS
    INTEGER:: eos1_id, eos2_id

    !
    !-- Parameters of the binary system
    !
    DOUBLE PRECISION:: angular_vel                  ! [rad/s]   
    ! Distance between the points of maximum baryon density
    DOUBLE PRECISION:: distance                     ! [km]
    ! Distance between the centers of mass
    DOUBLE PRECISION:: distance_com                 ! [Msun_geo]
    DOUBLE PRECISION:: mass1                        ! [Msun]
    DOUBLE PRECISION:: mass2                        ! [Msun]
    DOUBLE PRECISION:: mass_grav1                   ! [Msun]
    DOUBLE PRECISION:: mass_grav2                   ! [Msun]
    DOUBLE PRECISION:: adm_mass                     ! [Msun]
    ! For the following variable,
    ! mOmega= ( angular_vel[km^{-1}] )*( mass_grav1[km] + mass_grav2[km] )
    DOUBLE PRECISION:: mOmega                       ! [pure number]
    ! Estimated time of the merger
    DOUBLE PRECISION:: t_merger                     ! [Msun]
    DOUBLE PRECISION:: angular_momentum= 0.0D0      ! [G Msun^2 /c]
    ! Areal (or circumferential) radius of star 1
    ! Note that these are the areal radii of the stars in the binary system,
    ! which are different than those for an isolated star. The latter are used
    ! in the mass-radius diagrams
    DOUBLE PRECISION:: area_radius1                 ! [Msun_geo]
    ! Radius of star 1, in the x direction, towards the companion
    DOUBLE PRECISION:: radius1_x_comp               ! [Msun_geo]
    ! Radius of star 1, in the y direction
    DOUBLE PRECISION:: radius1_y                    ! [Msun_geo]
    ! Radius of star 1, in the z direction
    DOUBLE PRECISION:: radius1_z                    ! [Msun_geo]
    ! Radius of star 1, in the x direction, opposite to the companion
    DOUBLE PRECISION:: radius1_x_opp                ! [Msun_geo]
    ! Stellar center of star 1 (origin of the LORENE chart centered on star 1)
    DOUBLE PRECISION:: center1_x                    ! [Msun_geo]
    ! Barycenter of star 1
    DOUBLE PRECISION:: barycenter1_x                ! [Msun_geo]
    ! Areal (or circumferential) radius of star 2
    ! Note that these are the areal radii of the stars in the binary system,
    ! which are different than those for an isolated star. The latter are used
    ! in the mass-radius diagrams
    DOUBLE PRECISION:: area_radius2                 ! [Msun_geo]
    ! Radius of star 2, in the x direction, towards the companion
    DOUBLE PRECISION:: radius2_x_comp               ! [Msun_geo]
    ! Radius of star 2, in the y direction
    DOUBLE PRECISION:: radius2_y                    ! [Msun_geo]
    ! Radius of star 2, in the z direction
    DOUBLE PRECISION:: radius2_z                    ! [Msun_geo]
    ! Radius of star, in the x direction, opposite to the companion
    DOUBLE PRECISION:: radius2_x_opp                ! [Msun_geo]
    ! Stellar center of star 2 (origin of the LORENE chart centered on star 2)
    DOUBLE PRECISION:: center2_x                    ! [Msun_geo]
    ! Barycenter of star 2
    DOUBLE PRECISION:: barycenter2_x                ! [Msun_geo]
    ! Central enthalpy for star 1 [c^2]
    DOUBLE PRECISION:: ent_center1 ;
    ! Central baryon number density for star 1 [Msun_geo^-3]
    DOUBLE PRECISION:: nbar_center1 ;
    ! Central baryon mass density for star 1 [Msun Msun_geo^-3]
    DOUBLE PRECISION:: rho_center1 ;
    ! Central energy density for star 1 [Msun c^2 Msun_geo^-3]
    DOUBLE PRECISION:: energy_density_center1 ;
    ! Central specific energy for star 1 [c^2]
    DOUBLE PRECISION:: specific_energy_center1 ;
    ! Central pressure for star 1 [Msun c^2 Msun_geo^-3]
    DOUBLE PRECISION:: pressure_center1 ;
    ! Central enthalpy for star 2 [c^2]
    DOUBLE PRECISION:: ent_center2 ;
    ! Central baryon number density for star 2 [Msun_geo^-3]
    DOUBLE PRECISION:: nbar_center2 ;
    ! Central baryon mass density for star 2 [Msun Msun_geo^-3]
    DOUBLE PRECISION:: rho_center2 ;
    ! Central energy density for star 2 [Msun c^2 Msun_geo^-3]
    DOUBLE PRECISION:: energy_density_center2 ;
    ! Central specific energy for star 2 [c^2]
    DOUBLE PRECISION:: specific_energy_center2 ;
    ! Central pressure for star 2 [Msun c^2 Msun_geo^-3]
    DOUBLE PRECISION:: pressure_center2 ;
    ! Names of the equations of state (EoS) of the two neutron stars
    CHARACTER( LEN=: ), ALLOCATABLE:: eos1
    CHARACTER( LEN=: ), ALLOCATABLE:: eos2

    !
    !-- Parameters of polytropic equations of state for the two NSs
    !
    DOUBLE PRECISION:: gamma_1
    DOUBLE PRECISION:: gamma_2
    DOUBLE PRECISION:: kappa_1
    DOUBLE PRECISION:: kappa_2
    !
    !-- Parameters of the piecewise polytropic equation of state for NS 1
    !
    INTEGER:: npeos_1
    DOUBLE PRECISION:: gamma0_1
    DOUBLE PRECISION:: gamma1_1
    DOUBLE PRECISION:: gamma2_1
    DOUBLE PRECISION:: gamma3_1
    DOUBLE PRECISION:: kappa0_1
    DOUBLE PRECISION:: kappa1_1
    DOUBLE PRECISION:: kappa2_1
    DOUBLE PRECISION:: kappa3_1
    DOUBLE PRECISION:: logP1_1
    DOUBLE PRECISION:: logRho0_1
    DOUBLE PRECISION:: logRho1_1
    DOUBLE PRECISION:: logRho2_1
    !
    !-- Parameters of the piecewise polytropic equation of state for NS 2
    !
    INTEGER:: npeos_2
    DOUBLE PRECISION:: gamma0_2
    DOUBLE PRECISION:: gamma1_2
    DOUBLE PRECISION:: gamma2_2
    DOUBLE PRECISION:: gamma3_2
    DOUBLE PRECISION:: kappa0_2
    DOUBLE PRECISION:: kappa1_2
    DOUBLE PRECISION:: kappa2_2
    DOUBLE PRECISION:: kappa3_2
    DOUBLE PRECISION:: logP1_2
    DOUBLE PRECISION:: logRho0_2
    DOUBLE PRECISION:: logRho1_2
    DOUBLE PRECISION:: logRho2_2

    !
    !-- Spacetime fields
    !

    ! 1-D array storing the values of the lapse function
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: lapse
    ! 1-D arrays storing the values of the shift vector [c]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: shift_x
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: shift_y
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: shift_z
    ! 1-D arrays storing the values of the 3-metric components
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: g_xx
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: g_xy
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: g_xz
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: g_yy
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: g_yz
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: g_zz
    ! 1-D arrays storing the values of the extrinsic curvature components
    ! [c/MSun_geo] (see MODULE constants for the definition of MSun_geo)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: k_xx
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: k_xy
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: k_xz
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: k_yy
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: k_yz
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: k_zz

    !
    !-- Hydro fields
    !

    ! 1-D array storing the baryon mass density in the fluid frame
    ! [kg m^{-3}]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: baryon_density
    ! 1-D array storing the energy density [kg c^2 m^{-3}]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: energy_density
    ! 1-D array storing the specific internal energy [c^2]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: specific_energy
    ! 1-D arrays storing the fluid 3-velocity with respect to the Eulerian
    ! observer [c]
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: v_euler_x
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: v_euler_y
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: v_euler_z

    ! C pointer to the LORENE's Bin_NS object
    ! N.B. This variable is global. The pointer to the second LORENE Bin_NS
    !      object will overwrite the first one, and so on.
    !      This variable stores the pointer to the last defined LORENE Bin_NS
    !      object. That's why it is not freed in the destructor of a bns object.
    !      Presently, it has to be freed by the user at the end of the PROGRAM.
    !      See the last part of the PROGRAM in setup_lorene_id.f90, for example.
    TYPE(C_PTR):: bns_ptr

    ! Variables to set the geodesic gauge (lapse=1, shift=0)
    LOGICAL, PUBLIC:: one_lapse, zero_shift

    TYPE(timer), PUBLIC:: binary_construction_timer


    CONTAINS


    !-------------------!
    !--  SUBROUTINES  --!
    !-------------------!

    PROCEDURE:: construct_binary

    PROCEDURE:: destruct_binary

    PROCEDURE:: allocate_lorene_id_memory

    PROCEDURE:: deallocate_lorene_id_memory

    PROCEDURE:: import_id_params

    PROCEDURE:: integrate_baryon_mass_density

    PROCEDURE, PUBLIC:: print_id_params

    !
    !-- Overloaded SUBROUTINE that imports the LORENE BNS ID
    !-- in different ways: on the gravity grid, on the particles, etc...
    !
    GENERIC, PUBLIC:: import_id => &
                                   ! Store the ID in the bns arrays
                                   import_id_int_ptr, &
                                   ! Store the ID in non-member arrays with
                                   ! the same shape as the bns arrays
                                   import_id_ext_ptr, &
                                   ! Store the hydro ID in the arrays
                                   ! needed for the SPH part of SPHINCS
                                   import_id_particles_ptr, &
                                   ! Store the hydro ID needed to compute
                                   ! the baryon mass in variables
                                   import_id_mass_b_ptr, &
                                   ! Store the spacetime ID in multi-
                                   ! dimensional arrays needed for the
                                   ! spacetime part of SPHINCS
                                   import_id_multid_array_ptr, &
                                   ! Store the hydro ID in the arrays needed
                                   ! for the spacetime part of SPHINCS
                                   import_id_hydro_ptr, &
                                   ! Store the extrinsic curvature in the
                                   ! arrays needed for the SPH part of
                                   ! SPHINCS
                                   import_id_k_ptr
    PROCEDURE::       import_id_int_ptr          => import_id_int
    PROCEDURE::       import_id_ext_ptr          => import_id_ext
    PROCEDURE::       import_id_particles_ptr    => import_id_particles
    PROCEDURE::       import_id_mass_b_ptr       => import_id_mass_b
    PROCEDURE::       import_id_multid_array_ptr => import_id_multid_array
    PROCEDURE::       import_id_hydro_ptr        => import_id_hydro
    PROCEDURE::       import_id_k_ptr            => import_id_k

    !-----------------!
    !--  FUNCTIONS  --!
    !-----------------!

    ! Returns the LORENE's mass density at the given point
    PROCEDURE:: import_mass_density

    ! Returns the LORENE's conformally flat spatial ADM metric
    PROCEDURE:: import_spatial_metric

    ! Returns 1 if the energy density or the specific energy or the pressure
    ! are negative
    PROCEDURE:: is_hydro_negative

    !
    !-- Overloaded FUNCTION to access the fields as arrays and as values
    !
    GENERIC, PUBLIC:: get_field => get_fa, get_fv
    PROCEDURE::       get_fa    => get_field_array
    PROCEDURE::       get_fv    => get_field_value

    !
    !-- FUNCTIONS that access member variables
    !
    PROCEDURE, PUBLIC:: get_bns_identifier
    !PROCEDURE, PUBLIC:: get_bns_ptr

    PROCEDURE, PUBLIC:: get_angular_vel
    PROCEDURE, PUBLIC:: get_distance
    PROCEDURE, PUBLIC:: get_distance_com
    PROCEDURE, PUBLIC:: get_mass1
    PROCEDURE, PUBLIC:: get_mass2
    PROCEDURE, PUBLIC:: get_grav_mass1
    PROCEDURE, PUBLIC:: get_grav_mass2
    PROCEDURE, PUBLIC:: get_adm_mass
    PROCEDURE, PUBLIC:: get_angular_momentum
    PROCEDURE, PUBLIC:: get_radius1_x_comp
    PROCEDURE, PUBLIC:: get_radius1_y
    PROCEDURE, PUBLIC:: get_radius1_z
    PROCEDURE, PUBLIC:: get_radius1_x_opp
    PROCEDURE, PUBLIC:: get_center1_x
    PROCEDURE, PUBLIC:: get_barycenter1_x
    PROCEDURE, PUBLIC:: get_radius2_x_comp
    PROCEDURE, PUBLIC:: get_radius2_y
    PROCEDURE, PUBLIC:: get_radius2_z
    PROCEDURE, PUBLIC:: get_radius2_x_opp
    PROCEDURE, PUBLIC:: get_center2_x
    PROCEDURE, PUBLIC:: get_barycenter2_x
    PROCEDURE, PUBLIC:: get_ent_center1
    PROCEDURE, PUBLIC:: get_nbar_center1
    PROCEDURE, PUBLIC:: get_rho_center1
    PROCEDURE, PUBLIC:: get_energy_density_center1
    PROCEDURE, PUBLIC:: get_specific_energy_center1
    PROCEDURE, PUBLIC:: get_pressure_center1
    PROCEDURE, PUBLIC:: get_ent_center2
    PROCEDURE, PUBLIC:: get_nbar_center2
    PROCEDURE, PUBLIC:: get_rho_center2
    PROCEDURE, PUBLIC:: get_energy_density_center2
    PROCEDURE, PUBLIC:: get_specific_energy_center2
    PROCEDURE, PUBLIC:: get_pressure_center2
    PROCEDURE, PUBLIC:: get_eos1
    PROCEDURE, PUBLIC:: get_eos2
    PROCEDURE, PUBLIC:: get_eos1_id
    PROCEDURE, PUBLIC:: get_eos2_id

    ! PROCEDURES to be used for single polytropic EOS
    PROCEDURE, PUBLIC:: get_gamma_1
    PROCEDURE, PUBLIC:: get_gamma_2
    PROCEDURE, PUBLIC:: get_kappa_1
    PROCEDURE, PUBLIC:: get_kappa_2

    ! PROCEDURES to be used for piecewise polytropic EOS
    PROCEDURE, PUBLIC:: get_npeos_1
    PROCEDURE, PUBLIC:: get_gamma0_1
    PROCEDURE, PUBLIC:: get_gamma1_1
    PROCEDURE, PUBLIC:: get_gamma2_1
    PROCEDURE, PUBLIC:: get_gamma3_1
    PROCEDURE, PUBLIC:: get_kappa0_1
    PROCEDURE, PUBLIC:: get_kappa1_1
    PROCEDURE, PUBLIC:: get_kappa2_1
    PROCEDURE, PUBLIC:: get_kappa3_1
    PROCEDURE, PUBLIC:: get_logP1_1
    PROCEDURE, PUBLIC:: get_logRho0_1
    PROCEDURE, PUBLIC:: get_logRho1_1
    PROCEDURE, PUBLIC:: get_logRho2_1
    PROCEDURE, PUBLIC:: get_npeos_2
    PROCEDURE, PUBLIC:: get_gamma0_2
    PROCEDURE, PUBLIC:: get_gamma1_2
    PROCEDURE, PUBLIC:: get_gamma2_2
    PROCEDURE, PUBLIC:: get_gamma3_2
    PROCEDURE, PUBLIC:: get_kappa0_2
    PROCEDURE, PUBLIC:: get_kappa1_2
    PROCEDURE, PUBLIC:: get_kappa2_2
    PROCEDURE, PUBLIC:: get_kappa3_2
    PROCEDURE, PUBLIC:: get_logP1_2
    PROCEDURE, PUBLIC:: get_logRho0_2
    PROCEDURE, PUBLIC:: get_logRho1_2
    PROCEDURE, PUBLIC:: get_logRho2_2

    ! Destructor
    FINAL:: destruct_bns

  END TYPE bns

  !
  !-- Interface of the TYPE bns (i.e., declaration of the constructor)
  !-- (see https://dannyvanpoucke.be/oop-fortran-tut4-en/)
  !
  INTERFACE bns

    MODULE PROCEDURE:: construct_bns

  END INTERFACE bns

  !
  !-- Interfaces of the constructor and destructor of the TYPE bns
  !
  INTERFACE

    MODULE FUNCTION construct_bns( resu_file ) RESULT( bns_obj )

      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: resu_file
      TYPE(bns):: bns_obj

    END FUNCTION construct_bns

    MODULE SUBROUTINE destruct_bns( THIS )

      TYPE(bns), INTENT( IN OUT ):: THIS

    END SUBROUTINE destruct_bns

  END INTERFACE

  !
  !-- Interfaces of the methods of TYPE bns called by its constructor
  !-- Their implementations are in submodule_bns_constructor.f90
  !
  INTERFACE

    MODULE SUBROUTINE import_id_params( THIS )

      CLASS(bns), INTENT( IN OUT ):: THIS

    END SUBROUTINE import_id_params

  END INTERFACE

  !
  !-- Interfaces of the methods of the TYPE bns
  !-- Their implementations are in submodule_bns_methods.f90
  !
  INTERFACE

    !
    !-- SUBROUTINES
    !
    MODULE SUBROUTINE construct_binary( THIS, resu_file )

      CLASS(bns),                     INTENT( IN OUT )      :: THIS
      CHARACTER(KIND= C_CHAR, LEN=*), INTENT( IN ), OPTIONAL:: resu_file

    END SUBROUTINE construct_binary

    MODULE SUBROUTINE destruct_binary( THIS )

      CLASS(bns), INTENT( IN OUT ):: THIS

    END SUBROUTINE destruct_binary

    MODULE SUBROUTINE allocate_lorene_id_memory( THIS, d )

      CLASS(bns), INTENT( IN OUT ):: THIS
      INTEGER,    INTENT( IN )    :: d  ! dimension of the arrays

    END SUBROUTINE allocate_lorene_id_memory

    MODULE SUBROUTINE deallocate_lorene_id_memory( THIS )

      CLASS(bns), INTENT( IN OUT ):: THIS

    END SUBROUTINE deallocate_lorene_id_memory

    MODULE SUBROUTINE print_id_params( THIS )

      CLASS(bns), INTENT( IN OUT ):: THIS

    END SUBROUTINE print_id_params

    MODULE SUBROUTINE integrate_baryon_mass_density( THIS, center, radius, &
                                                     central_density, &
                                                     dr, dth, dphi, &
                                                     mass, mass_profile, &
                                                     mass_profile_idx )

      CLASS(bns), INTENT( IN OUT )      :: THIS
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT( IN OUT ):: mass_profile_idx
      DOUBLE PRECISION, INTENT( IN )    :: center, central_density, radius, &
                                           dr, dth, dphi
      DOUBLE PRECISION, INTENT( IN OUT ):: mass
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT( IN OUT ):: &
                                       mass_profile

    END SUBROUTINE integrate_baryon_mass_density

    MODULE SUBROUTINE import_id_int( THIS, n, x, y, z )

      CLASS(bns),                     INTENT( IN OUT ):: THIS
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

    MODULE SUBROUTINE import_id_ext( THIS, n, x, y, z,&
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

      CLASS(bns),                     INTENT( IN OUT ):: THIS
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

    END SUBROUTINE import_id_ext

    MODULE SUBROUTINE import_id_multid_array( THIS, nx, ny, nz, &
                                              pos, &
                                              lapse, &
                                              shift, &
                                              g, &
                                              ek )

      CLASS(bns),                           INTENT( IN OUT ):: THIS
      INTEGER,                              INTENT( IN )    :: nx
      INTEGER,                              INTENT( IN )    :: ny
      INTEGER,                              INTENT( IN )    :: nz
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN )    :: pos
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: lapse
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: shift
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: g
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: ek

    END SUBROUTINE import_id_multid_array

    MODULE SUBROUTINE import_id_hydro( THIS, nx, ny, nz, &
                                             pos, &
                                             baryon_density, &
                                             energy_density, &
                                             specific_energy, &
                                             pressure, &
                                             u_euler )

      CLASS(bns),                           INTENT( IN OUT ):: THIS
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

      CLASS(bns),                     INTENT( IN OUT ):: THIS
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

      CLASS(bns),       INTENT( IN OUT ):: THIS
      REAL(C_DOUBLE),   INTENT( IN )    :: x
      REAL(C_DOUBLE),   INTENT( IN )    :: y
      REAL(C_DOUBLE),   INTENT( IN)     :: z
      DOUBLE PRECISION, INTENT( IN OUT ):: g_xx
      DOUBLE PRECISION, INTENT( IN OUT ):: baryon_density
      DOUBLE PRECISION, INTENT( IN OUT ):: gamma_euler

    END SUBROUTINE import_id_mass_b

    MODULE SUBROUTINE import_id_k( THIS, n, x, y, z,&
                                         k_xx, k_xy, k_xz, &
                                         k_yy, k_yz, k_zz )

      CLASS(bns),                     INTENT( IN OUT ):: THIS
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

      ! Arguments
      CLASS(bns),     INTENT( IN )       :: THIS
      REAL(C_DOUBLE), INTENT( IN ), VALUE:: x
      REAL(C_DOUBLE), INTENT( IN ), VALUE:: y
      REAL(C_DOUBLE), INTENT( IN ), VALUE:: z
      ! Result
      REAL(C_DOUBLE):: res

    END FUNCTION import_mass_density

    MODULE FUNCTION import_spatial_metric( THIS, x, y, z ) RESULT( res )

      ! Arguments
      CLASS(bns),     INTENT( IN )       :: THIS
      REAL(C_DOUBLE), INTENT( IN ), VALUE:: x
      REAL(C_DOUBLE), INTENT( IN ), VALUE:: y
      REAL(C_DOUBLE), INTENT( IN ), VALUE:: z
      ! Result
      REAL(C_DOUBLE):: res

    END FUNCTION import_spatial_metric


    MODULE FUNCTION is_hydro_negative( THIS, x, y, z ) RESULT( res )

      ! Arguments
      CLASS(bns),     INTENT( IN )       :: THIS
      REAL(C_DOUBLE), INTENT( IN ), VALUE:: x
      REAL(C_DOUBLE), INTENT( IN ), VALUE:: y
      REAL(C_DOUBLE), INTENT( IN ), VALUE:: z
      ! Result
      INTEGER(C_INT):: res

    END FUNCTION is_hydro_negative

    MODULE FUNCTION get_field_array( THIS, field ) RESULT( field_array )

      ! Arguments
      CLASS(bns),          INTENT( IN )             :: THIS
      CHARACTER( LEN= : ), INTENT( IN ), ALLOCATABLE:: field
      ! Result
      DOUBLE PRECISION, DIMENSION(:),    ALLOCATABLE:: field_array

    END FUNCTION get_field_array

    MODULE FUNCTION get_field_value( THIS, field, n ) RESULT( field_value )

      ! Arguments
      CLASS(bns),          INTENT( IN )             :: THIS
      CHARACTER( LEN= : ), INTENT( IN ), ALLOCATABLE:: field
      INTEGER,             INTENT( IN )             :: n
      ! Result
      DOUBLE PRECISION                              :: field_value

    END FUNCTION get_field_value

    MODULE FUNCTION get_bns_identifier( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_bns_identifier

    END FUNCTION get_bns_identifier

    !MODULE FUNCTION get_bns_ptr( THIS )
    !
    !  ! Argument
    !  CLASS(bns), INTENT( IN ):: THIS
    !  ! Result
    !  TYPE(C_PTR):: get_bns_ptr
    !
    !END FUNCTION get_bns_ptr

    MODULE FUNCTION get_gamma_1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma_1

    END FUNCTION get_gamma_1

    MODULE FUNCTION get_gamma_2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma_2

    END FUNCTION get_gamma_2

    MODULE FUNCTION get_kappa_1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa_1

    END FUNCTION get_kappa_1

    MODULE FUNCTION get_kappa_2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa_2

    END FUNCTION get_kappa_2

    MODULE FUNCTION get_angular_vel( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_angular_vel

    END FUNCTION get_angular_vel

    MODULE FUNCTION get_distance( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_distance

    END FUNCTION get_distance

    MODULE FUNCTION get_distance_com( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_distance_com

    END FUNCTION get_distance_com

    MODULE FUNCTION get_mass1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_mass1

    END FUNCTION get_mass1

    MODULE FUNCTION get_mass2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_mass2

    END FUNCTION get_mass2

    MODULE FUNCTION get_grav_mass1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_grav_mass1

    END FUNCTION get_grav_mass1

    MODULE FUNCTION get_grav_mass2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_grav_mass2

    END FUNCTION get_grav_mass2

    MODULE FUNCTION get_adm_mass( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_adm_mass

    END FUNCTION get_adm_mass

    MODULE FUNCTION get_angular_momentum( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_angular_momentum

    END FUNCTION get_angular_momentum

    MODULE FUNCTION get_radius1_x_comp( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_radius1_x_comp

    END FUNCTION get_radius1_x_comp

    MODULE FUNCTION get_radius1_y( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_radius1_y

    END FUNCTION get_radius1_y

    MODULE FUNCTION get_radius1_z( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_radius1_z

    END FUNCTION get_radius1_z

    MODULE FUNCTION get_radius1_x_opp( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_radius1_x_opp

    END FUNCTION get_radius1_x_opp

    MODULE FUNCTION get_center1_x( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_center1_x

    END FUNCTION get_center1_x

    MODULE FUNCTION get_barycenter1_x( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_barycenter1_x

    END FUNCTION get_barycenter1_x

    MODULE FUNCTION get_radius2_x_comp( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_radius2_x_comp

    END FUNCTION get_radius2_x_comp

    MODULE FUNCTION get_radius2_y( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_radius2_y

    END FUNCTION get_radius2_y

    MODULE FUNCTION get_radius2_z( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_radius2_z

    END FUNCTION get_radius2_z

    MODULE FUNCTION get_radius2_x_opp( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_radius2_x_opp

    END FUNCTION get_radius2_x_opp

    MODULE FUNCTION get_center2_x( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_center2_x

    END FUNCTION get_center2_x

    MODULE FUNCTION get_barycenter2_x( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_barycenter2_x

    END FUNCTION get_barycenter2_x

    MODULE FUNCTION get_ent_center1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_ent_center1

    END FUNCTION get_ent_center1

    MODULE FUNCTION get_nbar_center1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_nbar_center1

    END FUNCTION get_nbar_center1

    MODULE FUNCTION get_rho_center1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_rho_center1

    END FUNCTION get_rho_center1

    MODULE FUNCTION get_energy_density_center1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_energy_density_center1

    END FUNCTION get_energy_density_center1

    MODULE FUNCTION get_specific_energy_center1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_specific_energy_center1

    END FUNCTION get_specific_energy_center1

    MODULE FUNCTION get_pressure_center1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_pressure_center1

    END FUNCTION get_pressure_center1

    MODULE FUNCTION get_ent_center2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_ent_center2

    END FUNCTION get_ent_center2

    MODULE FUNCTION get_nbar_center2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_nbar_center2

    END FUNCTION get_nbar_center2

    MODULE FUNCTION get_rho_center2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_rho_center2

    END FUNCTION get_rho_center2

    MODULE FUNCTION get_energy_density_center2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_energy_density_center2

    END FUNCTION get_energy_density_center2

    MODULE FUNCTION get_specific_energy_center2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_specific_energy_center2

    END FUNCTION get_specific_energy_center2

    MODULE FUNCTION get_pressure_center2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_pressure_center2

    END FUNCTION get_pressure_center2

    MODULE FUNCTION get_eos1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      CHARACTER( LEN= : ), ALLOCATABLE:: get_eos1

    END FUNCTION get_eos1

    MODULE FUNCTION get_eos2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      CHARACTER( LEN= : ), ALLOCATABLE:: get_eos2

    END FUNCTION get_eos2

    MODULE FUNCTION get_eos1_id( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      INTEGER:: get_eos1_id

    END FUNCTION get_eos1_id

    MODULE FUNCTION get_eos2_id( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      INTEGER:: get_eos2_id

    END FUNCTION get_eos2_id

    MODULE FUNCTION get_npeos_1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      INTEGER:: get_npeos_1

    END FUNCTION get_npeos_1

    MODULE FUNCTION get_npeos_2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      INTEGER:: get_npeos_2

    END FUNCTION get_npeos_2

    MODULE FUNCTION get_gamma0_1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma0_1

    END FUNCTION get_gamma0_1

    MODULE FUNCTION get_gamma1_1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma1_1

    END FUNCTION get_gamma1_1

    MODULE FUNCTION get_gamma2_1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma2_1

    END FUNCTION get_gamma2_1

    MODULE FUNCTION get_gamma3_1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma3_1

    END FUNCTION get_gamma3_1

    MODULE FUNCTION get_kappa0_1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa0_1

    END FUNCTION get_kappa0_1

    MODULE FUNCTION get_kappa1_1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa1_1

    END FUNCTION get_kappa1_1

    MODULE FUNCTION get_kappa2_1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa2_1

    END FUNCTION get_kappa2_1

    MODULE FUNCTION get_kappa3_1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa3_1

    END FUNCTION get_kappa3_1

    MODULE FUNCTION get_logP1_1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logP1_1

    END FUNCTION get_logP1_1

    MODULE FUNCTION get_logRho0_1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logRho0_1

    END FUNCTION get_logRho0_1

    MODULE FUNCTION get_logRho1_1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logRho1_1

    END FUNCTION get_logRho1_1

    MODULE FUNCTION get_logRho2_1( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logRho2_1

    END FUNCTION get_logRho2_1

    MODULE FUNCTION get_gamma0_2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma0_2

    END FUNCTION get_gamma0_2

    MODULE FUNCTION get_gamma1_2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma1_2

    END FUNCTION get_gamma1_2

    MODULE FUNCTION get_gamma2_2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma2_2

    END FUNCTION get_gamma2_2

    MODULE FUNCTION get_gamma3_2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma3_2

    END FUNCTION get_gamma3_2

    MODULE FUNCTION get_kappa0_2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa0_2

    END FUNCTION get_kappa0_2

    MODULE FUNCTION get_kappa1_2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa1_2

    END FUNCTION get_kappa1_2

    MODULE FUNCTION get_kappa2_2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa2_2

    END FUNCTION get_kappa2_2

    MODULE FUNCTION get_kappa3_2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa3_2

    END FUNCTION get_kappa3_2

    MODULE FUNCTION get_logP1_2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logP1_2

    END FUNCTION get_logP1_2

    MODULE FUNCTION get_logRho0_2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logRho0_2

    END FUNCTION get_logRho0_2

    MODULE FUNCTION get_logRho1_2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logRho1_2

    END FUNCTION get_logRho1_2

    MODULE FUNCTION get_logRho2_2( THIS )

      ! Argument
      CLASS(bns), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logRho2_2

    END FUNCTION get_logRho2_2

  END INTERFACE

  !----------------------------------------------------------!
  !--  Interfaces to the methods of LORENE's class Bin_NS  --!
  !----------------------------------------------------------!


  INTERFACE


    FUNCTION construct_bin_ns( c_resu_file ) RESULT( optr ) &
      BIND(C, NAME= "construct_bin_ns")

      !************************************************
      !                                               *
      ! Constructs the LORENE Bin_NS object           *
      !                                               *
      ! FT                                            *
      !                                               *
      !************************************************

      IMPORT :: C_PTR, C_CHAR

      IMPLICIT NONE

      ! Argument list
      CHARACTER(KIND= C_CHAR), DIMENSION(*), INTENT(IN), OPTIONAL :: &
                                                              c_resu_file
      ! Function result
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

      !************************************************
      !                                               *
      ! Import the full LORENE ID at the specified    *
      ! point. That is, import the metric fields, the *
      ! components of the extrinsic curvature [c/km], *
      ! and the hydro fields.                         *
      !                                               *
      ! - shift vector [c]                            *
      ! - baryon mass density [kg m^{-3}]             *
      ! - energy density [kg c^2 m^{-3}]              *
      ! - pressure [kg c^2 m^{-3}]                    *
      ! - specific internal energy [c^2]              *
      ! - fluid 3-velocity with respect to the        *
      !   Eulerian observer [c]                       *
      !                                               *
      ! FT                                            *
      !                                               *
      !************************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      TYPE(C_PTR),    INTENT(IN), VALUE :: optr
      REAL(C_DOUBLE), INTENT(IN), VALUE :: x
      REAL(C_DOUBLE), INTENT(IN), VALUE :: y
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

      !************************************************
      !                                               *
      ! Import the metric fields and the components   *
      ! of the extrinsic curvature [c/km] from LORENE,*
      ! at the specified point                        *
      !                                               *
      ! FT                                            *
      !                                               *
      !************************************************


      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      TYPE(C_PTR),    INTENT(IN), VALUE :: optr
      REAL(C_DOUBLE), INTENT(IN), VALUE :: x
      REAL(C_DOUBLE), INTENT(IN), VALUE :: y
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

      !************************************************
      !                                               *
      ! Import the hydro fields and the metric fields *
      ! from LORENE, at the specified point           *
      !                                               *
      ! - shift vector [c]                            *
      ! - baryon mass density [kg m^{-3}]             *
      ! - energy density [kg c^2 m^{-3}]              *
      ! - pressure [kg c^2 m^{-3}]                    *
      ! - specific internal energy [c^2]              *
      ! - fluid 3-velocity with respect to the        *
      !   Eulerian observer [c]                       *
      !                                               *
      ! FT                                            *
      !                                               *
      !************************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      TYPE(C_PTR),    INTENT(IN), VALUE :: optr
      REAL(C_DOUBLE), INTENT(IN), VALUE :: x
      REAL(C_DOUBLE), INTENT(IN), VALUE :: y
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
      !                                               *
      ! Import the hydro fields and the metric fields *
      ! from LORENE, at the specified point,          *
      ! needed to compute the baryon mass.             *
      !                                               *
      ! - shift vector [c]                            *
      ! - baryon mass density [kg m^{-3}]             *
      ! - fluid 3-velocity with respect to the        *
      !   Eulerian observer [c]                       *
      !                                               *
      ! FT                                            *
      !                                               *
      !************************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      TYPE(C_PTR),    INTENT(IN), VALUE :: optr
      REAL(C_DOUBLE), INTENT(IN), VALUE :: x
      REAL(C_DOUBLE), INTENT(IN), VALUE :: y
      REAL(C_DOUBLE), INTENT(IN), VALUE :: z
      REAL(C_DOUBLE), INTENT(OUT)       :: g_diag
      REAL(C_DOUBLE), INTENT(OUT)       :: baryon_density
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

      !************************************************
      !                                               *
      ! Import the hydro fields from LORENE, at the   *
      ! specified point                               *
      !                                               *
      ! - baryon mass density [kg m^{-3}]             *
      ! - energy density [kg c^2 m^{-3}]              *
      ! - pressure [kg c^2 m^{-3}]                    *
      ! - specific internal energy [c^2]              *
      ! - fluid 3-velocity with respect to the        *
      !   Eulerian observer [c]                       *
      !                                               *
      ! FT                                            *
      !                                               *
      !************************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      ! Argument list
      TYPE(C_PTR),    INTENT(IN), VALUE :: optr
      REAL(C_DOUBLE), INTENT(IN), VALUE :: x
      REAL(C_DOUBLE), INTENT(IN), VALUE :: y
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

      !************************************************
      !                                               *
      ! Import the components of the extrinsic        *
      ! curvature [c/km] from LORENE, at the          *
      ! specified point                               *
      !                                               *
      ! FT                                            *
      !                                               *
      !************************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      ! Argument list
      TYPE(C_PTR),    INTENT(IN), VALUE :: optr
      REAL(C_DOUBLE), INTENT(IN), VALUE :: x
      REAL(C_DOUBLE), INTENT(IN), VALUE :: y
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

      !************************************************
      !                                               *
      ! Return the baryon mass density [kg m^{-3}]    *
      ! from LORENE, at the specified point           *
      !                                               *
      ! FT                                            *
      !                                               *
      !************************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      ! Argument list
      TYPE(C_PTR),    INTENT(IN),  VALUE :: optr
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: x
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: y
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: z
      ! Function result
      REAL(C_DOUBLE) :: res

    END FUNCTION get_lorene_mass_density


    FUNCTION get_lorene_spatial_metric( optr, x, y, z ) RESULT( res ) &
      BIND(C, NAME= "get_lorene_id_g")

      !************************************************
      !                                               *
      ! Return the diagonal components of the metric, *
      ! all equal to the LORENE conformal factor to   *
      ! the 4th power.
      !                                               *
      ! FT                                            *
      !                                               *
      !************************************************

      IMPORT :: C_DOUBLE, C_PTR

      IMPLICIT NONE

      ! Argument list
      TYPE(C_PTR),    INTENT(IN),  VALUE :: optr
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: x
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: y
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: z
      ! Function result
      REAL(C_DOUBLE) :: res

    END FUNCTION get_lorene_spatial_metric


    FUNCTION negative_hydro( optr, x, y, z ) RESULT( res ) &
      BIND(C, NAME= "negative_hydro")

      !************************************************
      !                                               *
      ! Return 1 if the energy density is nonpositive,*
      ! or if the specific energy is nonpositive,     *
      ! or if the pressure i nonpositive              *
      ! at the specified point, and 0 otherwise       *
      !                                               *
      ! FT 12.03.2021                                 *
      !                                               *
      !************************************************

      IMPORT :: C_INT, C_DOUBLE, C_PTR

      IMPLICIT NONE

      ! Argument list
      TYPE(C_PTR),    INTENT(IN),  VALUE :: optr
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: x
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: y
      REAL(C_DOUBLE), INTENT(IN),  VALUE :: z
      ! Function result
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

      !************************************************
      !                                               *
      ! Import the physical parameters of the binary  *
      ! system from LORENE                            *
      !                                               *
      ! FT                                            *
      !                                               *
      !************************************************

      IMPORT :: C_INT, C_DOUBLE, C_PTR, C_CHAR

      IMPLICIT NONE

      ! Argument list
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

      !************************************************
      !                                               *
      ! Destructs the LORENE Bin_NS object           *
      !                                               *
      ! FT                                            *
      !                                               *
      !************************************************

      IMPORT :: C_PTR

      IMPLICIT NONE

      ! Argument list
      TYPE(C_PTR), INTENT(IN), VALUE :: optr

    END SUBROUTINE destruct_bin_ns


  END INTERFACE


END MODULE bns_id
