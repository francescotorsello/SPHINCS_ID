! File:         module_ejecta_generic.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

MODULE ejecta_generic

  !********************************************************
  !
  !# This MODULE contains the definition of TYPE ejecta,
  !  which is an ABSTRACT TYPE representing any possible
  !  type of initial data (|id|) for a differentially rotating
  !  star (DRS) to be set up for |sphincsbssn|. That is, DRS |id|
  !  produced with |lorene|, with |fuka|, etc.
  !
  !  PROCEDURES and variables shared by all the types
  !  of DRS |id| should belong to TYPE ejecta, as
  !  they are inherited by its EXTENDED TYPES that
  !  represent more specific types of DRS |id|.
  !
  !  FT 22.10.2021
  !
  !********************************************************


  USE id_base, ONLY: idbase
  USE utility, ONLY: ios, err_msg


  IMPLICIT NONE


  !********************************************************************
  !                                                                   *
  !  Definition of TYPE ejecta  *
  !                                                                   *
  !********************************************************************

  TYPE, EXTENDS(idbase):: ejecta
  !# TYPE for ejecta |id| for |sphincsbssn| prepared on a grid

    INTEGER:: nx_grid
    !! Number of grid points in the \(x\) direction for the grid containing the |id|

    INTEGER:: ny_grid
    !! Number of grid points in the \(y\) direction for the grid containing the |id|

    INTEGER:: nz_grid
    !! Number of grid points in the \(z\) direction for the grid containing the |id|

    INTEGER:: n_gridpoints
    !! Total number of grid points for the grid containing the |id|

    DOUBLE PRECISION:: xL_grid
    !! Minimum \(x\) coordinate on the grid containing the |id|

    DOUBLE PRECISION:: yL_grid
    !! Minimum \(y\) coordinate on the grid containing the |id|

    DOUBLE PRECISION:: zL_grid
    !! Minimum \(z\) coordinate on the grid containing the |id|

    DOUBLE PRECISION:: dx_grid
    !! Spacing on the \(x\)-axis for the grid containing the |id|

    DOUBLE PRECISION:: dy_grid
    !! Spacing on the \(y\)-axis for the grid containing the |id|

    DOUBLE PRECISION:: dz_grid
    !! Spacing on the \(z\)-axis for the grid containing the |id|

    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: grid
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER:: baryon_mass_density

    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: masses
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: sizes
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: centers
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: barycenters

    !--------------------------------!
    !--  Parameters of the ejecta  --!
    !--------------------------------!

    CHARACTER( LEN=: ), ALLOCATABLE:: eos
    !! Name of the equation of state (EoS) of star 1


    !
    !-- Parameters of single polytropic equations of state for the two NSs
    !

    DOUBLE PRECISION:: gamma
    !! Single polytrope: polytropic index

    DOUBLE PRECISION:: kappa
    !! Single polytrope: polytropic constant [pure number]

    !
    !-- Parameters of the piecewise polytropic equation of state for NS 1
    !

    INTEGER:: npeos
    !! Piecewise polytrope: Number of polytropic pieces

    DOUBLE PRECISION:: gamma0
    !! Piecewise polytrope: polytropic index \(\gamma_0\)

    DOUBLE PRECISION:: gamma1
    !! Piecewise polytrope: polytropic index \(\gamma_1\)

    DOUBLE PRECISION:: gamma2
    !! Piecewise polytrope: polytropic index \(\gamma_2\)

    DOUBLE PRECISION:: gamma3
    !! Piecewise polytrope: polytropic index \(\gamma_3\)

    DOUBLE PRECISION:: kappa0
    !# Piecewise polytrope: polytropic constant \(\kappa_0\)
    !  [pure number]

    DOUBLE PRECISION:: kappa1
    !# Piecewise polytrope: polytropic constant \(\kappa_1\)
    !  [pure number]

    DOUBLE PRECISION:: kappa2
    !# Piecewise polytrope: polytropic constant \(\kappa_2\)
    !  [pure number]

    DOUBLE PRECISION:: kappa3
    !# Piecewise polytrope: polytropic constant \(\kappa_3\)
    !  [pure number]

    DOUBLE PRECISION:: logP1
    !# Piecewise polytrope: Base 10 exponent of the pressure at the first
    !  fiducial density (between \(\gamma_0\) and \(\gamma_1\))
    !  \([{\rm dyne/cm^2}]\)

    DOUBLE PRECISION:: logRho0
    !# Piecewise polytrope: Base 10 exponent of the first fiducial density
    !  (between \(\gamma_0\) and \(\gamma_1\)) \([{\rm g/cm^3}]\)

    DOUBLE PRECISION:: logRho1
    !# Piecewise polytrope: Base 10 exponent of the second fiducial density
    !  (between \(\gamma_1\) and \(\gamma_2\)) \([{\rm g/cm^3}]\)

    DOUBLE PRECISION:: logRho2
    !# Piecewise polytrope: Base 10 exponent of the third fiducial density
    !  (between \(\gamma_2\) and \(\gamma_3\)) \([{\rm g/cm^3}]\)



    CONTAINS


    !-------------------!
    !--  SUBROUTINES  --!
    !-------------------!

    PROCEDURE:: derived_type_constructor => construct_ejecta

    !PROCEDURE(get_eos_id_int), DEFERRED:: get_eos_id
    !! Returns an integer that identifies the equation of state

    PROCEDURE:: read_id_full      => interpolate_id_full
    PROCEDURE:: read_id_spacetime => interpolate_id_spacetime
    PROCEDURE:: read_id_particles => interpolate_id_particles
    PROCEDURE:: read_id_hydro     => interpolate_id_hydro
    PROCEDURE:: read_id_mass_b    => interpolate_id_mass_b
    PROCEDURE:: read_id_k         => interpolate_id_k


    !-----------------!
    !--  FUNCTIONS  --!
    !-----------------!

    PROCEDURE:: read_mass_density => interpolate_mass_density
    !! Returns the |lorene|'s mass density at the given point

    !PROCEDURE:: interpolate_spatial_metric
    !! Returns the |lorene|'s conformally flat spatial ADM metric

    PROCEDURE:: test_position => is_hydro_negative
    !# Returns 1 if the energy density or the specific energy or the pressure
    !  are negative


    !
    !-- FUNCTIONS that access PRIVATE member variables
    !

    PROCEDURE:: return_mass                 => get_mass
    PROCEDURE:: return_center               => get_center
    PROCEDURE:: return_barycenter           => get_barycenter
    PROCEDURE:: return_eos_name             => get_eos
    PROCEDURE:: return_spatial_extent       => get_radii
    PROCEDURE:: print_summary               => print_summary_ejecta

    PROCEDURE:: get_eos_id => get_eos_ejectaid
    !! Returns the identifier for the EOS

    PROCEDURE:: return_eos_parameters => get_eos_parameters

    !
    !-- PROCEDURES to be used for single polytropic EOS
    !
    PROCEDURE, PUBLIC:: get_gamma
    !! Returns [[ejecta:gamma]]
    PROCEDURE, PUBLIC:: get_kappa
    !! Returns [[ejecta:kappa]]

    !
    !-- PROCEDURES to be used for piecewise polytropic EOS
    !
    PROCEDURE, PUBLIC:: get_npeos
    !! Returns [[ejecta:npeos]]
    PROCEDURE, PUBLIC:: get_gamma0
    !! Returns [[ejecta:gamma0]]
    PROCEDURE, PUBLIC:: get_gamma1
    !! Returns [[ejecta:gamma1]]
    PROCEDURE, PUBLIC:: get_gamma2
    !! Returns [[ejecta:gamma2]]
    PROCEDURE, PUBLIC:: get_gamma3
    !! Returns [[ejecta:gamma3]]
    PROCEDURE, PUBLIC:: get_kappa0
    !! Returns [[ejecta:kappa0]]
    PROCEDURE, PUBLIC:: get_kappa1
    !! Returns [[ejecta:kappa1]]
    PROCEDURE, PUBLIC:: get_kappa2
    !! Returns [[ejecta:kappa2]]
    PROCEDURE, PUBLIC:: get_kappa3
    !! Returns [[ejecta:kappa3]]
    PROCEDURE, PUBLIC:: get_logP1
    !! Returns [[ejecta:logP1]]
    PROCEDURE, PUBLIC:: get_logRho0
    !! Returns [[ejecta:logRho0]]
    PROCEDURE, PUBLIC:: get_logRho1
    !! Returns [[ejecta:logRho1]]
    PROCEDURE, PUBLIC:: get_logRho2
    !! Returns [[ejecta:logRho2]]

    FINAL:: destruct_ejecta
    !! Finalizer (Destructor) of a [[ejecta]] object

  END TYPE ejecta


  INTERFACE


    !------------------------------!
    !--  OVERRIDING SUBROUTINES  --!
    !------------------------------!


    MODULE SUBROUTINE print_summary_ejecta( THIS, filename )
    !# Prints a summary of the physical properties the system
    !  to the standard output and, optionally, to a formatted file whose name
    !  is given as the optional argument `filename`


      CLASS(ejecta), INTENT( IN ):: THIS
      CHARACTER( LEN= * ), INTENT( INOUT ), OPTIONAL:: filename
      !! Name of the formatted file to print the summary to

    END SUBROUTINE print_summary_ejecta


    !----------------------------!
    !--  OVERRIDING FUNCTIONS  --!
    !----------------------------!


    MODULE FUNCTION get_mass( THIS, i_matter )
    !! Returns [[ejecta:mass]]

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT( IN ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      ! Result
      DOUBLE PRECISION:: get_mass

    END FUNCTION get_mass


    MODULE FUNCTION get_center( THIS, i_matter )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT( IN ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION, DIMENSION(3):: get_center

    END FUNCTION get_center


    MODULE FUNCTION get_barycenter( THIS, i_matter )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT( IN ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION, DIMENSION(3):: get_barycenter

    END FUNCTION get_barycenter


    MODULE FUNCTION get_radii( THIS, i_matter )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT( IN ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose string is to return
      DOUBLE PRECISION, DIMENSION(6):: get_radii

    END FUNCTION get_radii


    !-------------------!
    !--  SUBROUTINES  --!
    !-------------------!

    MODULE SUBROUTINE construct_ejecta( derived_type, filename )
    !! Constructs a [[ejecta]] object
    !# Prints a summary of the physical properties the system
    !  to the standard output and, optionally, to a formatted file whose name
    !  is given as the optional argument `filename`

      CHARACTER(LEN=*), INTENT( IN ), OPTIONAL :: filename
      !! |lorene| binary file containing the spectral DRS ID
      CLASS(ejecta), INTENT( OUT ):: derived_type
      !! Constructed [[ejecta]] object

    END SUBROUTINE construct_ejecta


    MODULE SUBROUTINE destruct_ejecta( THIS )
    !! Destruct a [[ejecta]] object

      TYPE(ejecta), INTENT( IN OUT ):: THIS
      !! [[ejecta]] object to be destructed

    END SUBROUTINE destruct_ejecta


    MODULE SUBROUTINE interpolate_id_full( THIS, n, x, y, z,&
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
    !# Stores the ID in non [[ejecta]]-member arrays with the same
    !  shape as the [[ejecta]] member arrays

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta),                     INTENT( IN OUT ):: THIS
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

    END SUBROUTINE interpolate_id_full


    MODULE SUBROUTINE interpolate_id_spacetime( THIS, nx, ny, nz, &
                                              pos, &
                                              lapse, &
                                              shift, &
                                              g, &
                                              ek )
    !# Stores the spacetime ID in multi-dimensional arrays needed to compute
    !  the BSSN variables and constraints

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta),                           INTENT( IN OUT ):: THIS
      INTEGER,                              INTENT( IN )    :: nx
      INTEGER,                              INTENT( IN )    :: ny
      INTEGER,                              INTENT( IN )    :: nz
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN )    :: pos
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: lapse
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: shift
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: g
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: ek

    END SUBROUTINE interpolate_id_spacetime


    MODULE SUBROUTINE interpolate_id_hydro( THIS, nx, ny, nz, &
                                             pos, &
                                             baryon_density, &
                                             energy_density, &
                                             specific_energy, &
                                             pressure, &
                                             u_euler )
    !# Stores the hydro ID in the arrays needed to compute the constraints
    !  on the refined mesh

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta),                           INTENT( IN OUT ):: THIS
      INTEGER,                              INTENT( IN )    :: nx
      INTEGER,                              INTENT( IN )    :: ny
      INTEGER,                              INTENT( IN )    :: nz
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN )    :: pos
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: baryon_density
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: energy_density
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: specific_energy
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN OUT ):: pressure
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: u_euler

    END SUBROUTINE interpolate_id_hydro


    MODULE SUBROUTINE interpolate_id_particles( THIS, n, x, y, z, &
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

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta),                     INTENT( IN OUT ):: THIS
      INTEGER,                        INTENT( IN )    :: n
      DOUBLE PRECISION,   DIMENSION(:), INTENT( IN )    :: x
      DOUBLE PRECISION,   DIMENSION(:), INTENT( IN )    :: y
      DOUBLE PRECISION,   DIMENSION(:), INTENT( IN )    :: z
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

    END SUBROUTINE interpolate_id_particles


    MODULE SUBROUTINE interpolate_id_mass_b( THIS, x, y, z, &
                                        g, &
                                        baryon_density, &
                                        gamma_euler )
    !! Stores the hydro ID in the arrays needed to compute the baryon mass

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta),       INTENT( IN OUT ):: THIS
      DOUBLE PRECISION, INTENT( IN )    :: x
      DOUBLE PRECISION, INTENT( IN )    :: y
      DOUBLE PRECISION, INTENT( IN)     :: z
      DOUBLE PRECISION, DIMENSION(6), INTENT( OUT ):: g
      DOUBLE PRECISION, INTENT( OUT ):: baryon_density
      DOUBLE PRECISION, INTENT( OUT ):: gamma_euler

    END SUBROUTINE interpolate_id_mass_b


    MODULE SUBROUTINE interpolate_id_k( THIS, n, x, y, z,&
                                         k_xx, k_xy, k_xz, &
                                         k_yy, k_yz, k_zz )
    !! Stores the components of the extrinsic curvature in arrays

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta),                     INTENT( IN OUT ):: THIS
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

    END SUBROUTINE interpolate_id_k


    !
    !-- FUNCTIONS
    !
    MODULE FUNCTION interpolate_mass_density( THIS, x, y, z ) RESULT( res )
    !! Returns the |lorene| baryon mass density at a point \((x,y,z)\)

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta),     INTENT( IN )         :: THIS
      !> \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: x
      !> \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: y
      !> \(z\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: z
      !> Baryon mass density at \((x,y,z)\)
      DOUBLE PRECISION:: res

    END FUNCTION interpolate_mass_density


    MODULE FUNCTION interpolate_spatial_metric( THIS, x, y, z ) RESULT( res )
    !# Returns the |lorene| conformally flat spatial metric component
    !  \(g_{xx}=g_{yy}=g_{zz}\) at a point \((x,y,z)\)

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta),     INTENT( IN )       :: THIS
      !> \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: x
      !> \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: y
      !> \(z\) coordinate of the desired point
      DOUBLE PRECISION, INTENT( IN ), VALUE:: z
      !> \(g_{xx}=g_{yy}=g_{zz}\) at \((x,y,z)\)
      DOUBLE PRECISION:: res

    END FUNCTION interpolate_spatial_metric


    MODULE FUNCTION is_hydro_negative( THIS, x, y, z ) RESULT( res )
    !# Returns 1 if the energy density or the specific energy or the pressure
    !  are negative, 0 otherwise

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta),     INTENT( IN )       :: THIS
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


    !-----------------!
    !--  FUNCTIONS  --!
    !-----------------!


    MODULE PURE FUNCTION get_gamma( THIS )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma

    END FUNCTION get_gamma


    MODULE PURE FUNCTION get_kappa( THIS )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa

    END FUNCTION get_kappa


    MODULE FUNCTION get_eos( THIS, i_matter )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT( IN ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose string is to return
      CHARACTER( LEN= : ), ALLOCATABLE:: get_eos

    END FUNCTION get_eos


    MODULE PURE FUNCTION get_npeos( THIS )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT( IN ):: THIS
      ! Result
      INTEGER:: get_npeos

    END FUNCTION get_npeos


    MODULE PURE FUNCTION get_gamma0( THIS )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma0

    END FUNCTION get_gamma0


    MODULE PURE FUNCTION get_gamma1( THIS )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma1

    END FUNCTION get_gamma1


    MODULE PURE FUNCTION get_gamma2( THIS )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma2

    END FUNCTION get_gamma2


    MODULE PURE FUNCTION get_gamma3( THIS )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma3

    END FUNCTION get_gamma3


    MODULE PURE FUNCTION get_kappa0( THIS )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa0

    END FUNCTION get_kappa0


    MODULE PURE FUNCTION get_kappa1( THIS )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa1

    END FUNCTION get_kappa1


    MODULE PURE FUNCTION get_kappa2( THIS )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa2

    END FUNCTION get_kappa2


    MODULE PURE FUNCTION get_kappa3( THIS )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa3

    END FUNCTION get_kappa3


    MODULE PURE FUNCTION get_logP1( THIS )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logP1

    END FUNCTION get_logP1


    MODULE PURE FUNCTION get_logRho0( THIS )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logRho0

    END FUNCTION get_logRho0


    MODULE PURE FUNCTION get_logRho1( THIS )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logRho1

    END FUNCTION get_logRho1


    MODULE PURE FUNCTION get_logRho2( THIS )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logRho2

    END FUNCTION get_logRho2


    MODULE FUNCTION get_eos_ejectaid( THIS )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT( IN ):: THIS
      ! Result
      INTEGER:: get_eos_ejectaid

    END FUNCTION get_eos_ejectaid


    MODULE SUBROUTINE get_eos_parameters( THIS, i_matter, eos_params )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT( IN ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT):: eos_params
      !# Array containing the parameters of the |eos| for the DRS

    END SUBROUTINE get_eos_parameters


  END INTERFACE


END MODULE ejecta_generic

