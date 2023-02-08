! File:         module_ejecta_generic.f90
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

MODULE ejecta_generic

  !********************************************************
  !
  !# This MODULE contains the definition of TYPE ejecta,
  !  which is an ABSTRACT TYPE representing any possible
  !  type of initial data (|id|) on a Cartesian, uniform grid
  !  to be set up for |sphincsbssn|.
  !
  !  PROCEDURES and variables shared by all the types
  !  of |id| on a Cartesian uniform grid should belong to
  !  TYPE ejecta.
  !
  !  FT xx.11.2021
  !
  !********************************************************


  USE id_base, ONLY: idbase
  USE utility, ONLY: ios, err_msg


  IMPLICIT NONE


  !******************************
  !                             *
  !  Definition of TYPE ejecta  *
  !                             *
  !******************************

  TYPE, EXTENDS(idbase):: ejecta
  !# TYPE for |id| for |sphincsbssn| prepared on a Cartesian, uniform grid

    INTEGER:: nx_grid
    !# Number of grid points in the \(x\) direction for the grid
    !  containing the |id|

    INTEGER:: ny_grid
    !# Number of grid points in the \(y\) direction for the grid
    !  containing the |id|

    INTEGER:: nz_grid
    !# Number of grid points in the \(z\) direction for the grid
    !  containing the |id|

    INTEGER:: n_gridpoints
    !! Total number of grid points for the grid containing the |id|

    DOUBLE PRECISION:: xL_grid
    !! Minimum \(x\) coordinate on the grid containing the |id|

    DOUBLE PRECISION:: yL_grid
    !! Minimum \(y\) coordinate on the grid containing the |id|

    DOUBLE PRECISION:: zL_grid
    !! Minimum \(z\) coordinate on the grid containing the |id|

    DOUBLE PRECISION:: xR_grid
    !! Maximum \(x\) coordinate on the grid containing the |id|

    DOUBLE PRECISION:: yR_grid
    !! Maximum \(y\) coordinate on the grid containing the |id|

    DOUBLE PRECISION:: zR_grid
    !! Maximum \(z\) coordinate on the grid containing the |id|

    DOUBLE PRECISION, PUBLIC:: dx_grid
    !! Spacing on the \(x\)-axis for the grid containing the |id|

    DOUBLE PRECISION:: dy_grid
    !! Spacing on the \(y\)-axis for the grid containing the |id|

    DOUBLE PRECISION:: dz_grid
    !! Spacing on the \(z\)-axis for the grid containing the |id|

    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: grid
    !# Array storing the Cartesian coordinates of the grid points.
    !  The first three indices specify the grid point; the last index
    !  specifies the \(x,y,z\) coordinates

    DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE:: baryon_mass_density
    !# Array storing the baryon mass density at the grid points.
    !  The indices specify the grid point.

    DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE:: specific_energy
    !# Array storing the specific energy at the grid points.
    !  The indices specify the grid point.

    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: vel
    !# Array storing the fluid velocity with respect to the Eulerian observer
    !  at the grid points.
    !  The first three indices specify the grid point; the last index
    !  specifies the \(x,y,z\) components of the velocity.

    DOUBLE PRECISION:: adm_mass
    !! ADM mass of the system \([M_\odot]\)

    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: masses
    !! Masses of the matter objects

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: sizes
    !! Sizes of the matter objects

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: centers
    !! Centers of the matter objects

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: barycenters
    !! Barycenters of the matter objects

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


    INTEGER:: eos_id
    !! |sphincsid| identifier for the |eos|


    CONTAINS


    !-------------------!
    !--  SUBROUTINES  --!
    !-------------------!

    PROCEDURE:: derived_type_constructor => construct_ejecta

    PROCEDURE:: allocate_gridid_memory
    !! Allocates memory for the [[ejecta]] member arrays

    PROCEDURE:: deallocate_gridid_memory
    !! Deallocates memory for the [[ejecta]] member arrays

    PROCEDURE:: read_id_full      => interpolate_id_full
    PROCEDURE:: read_id_spacetime => interpolate_id_spacetime
    PROCEDURE:: read_id_particles => interpolate_id_particles
    PROCEDURE:: read_id_hydro     => interpolate_id_hydro
    PROCEDURE:: read_id_mass_b    => interpolate_id_mass_b
    PROCEDURE:: read_id_k         => interpolate_id_k

    PROCEDURE, NOPASS:: finalize
    !# Corrects the |sph| |id| so that the linear \(\mathrm{ADM}\) momentum
    !  is zero

    PROCEDURE:: nothing
    !# Procedure that does nothing. It is used to instantiate a deferred
    !  idbase procedure which is not needed in TYPE [[ejecta]].
    !  It also serves as a placeholder in case the idbase procedure
    !  will be needed in the future.

    PROCEDURE:: initialize_id => nothing

    !-----------------!
    !--  FUNCTIONS  --!
    !-----------------!

    PROCEDURE:: read_mass_density => interpolate_mass_density
    !! Returns the |lorene|'s mass density at the given point

    PROCEDURE:: read_pressure => interpolate_pressure
    !! Returns the |lorene|'s mass density at the given point

    !PROCEDURE:: interpolate_spatial_metric
    !! Returns the |lorene|'s conformally flat spatial ADM metric

    PROCEDURE:: test_position => is_hydro_positive
    !# Returns 1 if the energy density or the specific energy or the pressure
    !  are negative


    !
    !-- FUNCTIONS that access PRIVATE member variables
    !

    PROCEDURE:: return_mass                 => get_mass
    PROCEDURE:: return_adm_mass             => get_adm_mass
    PROCEDURE:: return_center               => get_center
    PROCEDURE:: return_barycenter           => get_barycenter
    PROCEDURE:: return_eos_name             => get_eos
    PROCEDURE:: return_spatial_extent       => get_radii
    PROCEDURE:: print_summary               => print_summary_ejecta

    !PROCEDURE:: get_eos_id => get_eos_ejectaid
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


    MODULE SUBROUTINE print_summary_ejecta( this, filename )
    !# Prints a summary of the physical properties the system
    !  to the standard output and, optionally, to a formatted file whose name
    !  is given as the optional argument `filename`


      CLASS(ejecta),    INTENT(IN):: this
      CHARACTER(LEN=*), INTENT(INOUT), OPTIONAL:: filename
      !! Name of the formatted file to print the summary to

    END SUBROUTINE print_summary_ejecta


    !----------------------------!
    !--  OVERRIDING FUNCTIONS  --!
    !----------------------------!


    MODULE FUNCTION get_mass( this, i_matter )
    !! Returns [[ejecta:masses]]

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      INTEGER,       INTENT(IN):: i_matter
      ! Result
      DOUBLE PRECISION:: get_mass

    END FUNCTION get_mass


    MODULE FUNCTION get_adm_mass( this )
    !! Returns 0 (the ADM mass is not necessarily known for this TYPE)

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      ! Result
      DOUBLE PRECISION:: get_adm_mass

    END FUNCTION get_adm_mass


    MODULE FUNCTION get_center( this, i_matter )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      INTEGER,       INTENT(IN):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION, DIMENSION(3):: get_center

    END FUNCTION get_center


    MODULE FUNCTION get_barycenter( this, i_matter )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      INTEGER,       INTENT(IN):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION, DIMENSION(3):: get_barycenter

    END FUNCTION get_barycenter


    MODULE FUNCTION get_radii( this, i_matter )

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      INTEGER,       INTENT(IN):: i_matter
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

      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: filename
      !! |lorene| binary file containing the spectral DRS ID
      CLASS(ejecta),    INTENT(OUT):: derived_type
      !! Constructed [[ejecta]] object

    END SUBROUTINE construct_ejecta


    MODULE SUBROUTINE destruct_ejecta( this )
    !! Destruct a [[ejecta]] object

      TYPE(ejecta), INTENT(INOUT):: this
      !! [[ejecta]] object to be destructed

    END SUBROUTINE destruct_ejecta


    MODULE SUBROUTINE allocate_gridid_memory( this, n_matter )
    !! Allocates allocatable arrays member of a [[ejecta]] object

      CLASS(ejecta), INTENT(INOUT):: this
      !! [[ejecta]] object which this PROCEDURE is a member of
      INTEGER,       INTENT(IN):: n_matter
      !! Number of matter objects

    END SUBROUTINE allocate_gridid_memory


    MODULE SUBROUTINE deallocate_gridid_memory( this )
    !! Deallocates allocatable arrays member of a [[ejecta]] object

      CLASS(ejecta), INTENT(INOUT):: this
      !! [[ejecta]] object which this PROCEDURE is a member of

    END SUBROUTINE deallocate_gridid_memory


    MODULE SUBROUTINE interpolate_id_full( this, n, x, y, z,&
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
    !# Stores the ID in non [[ejecta]]-member arrays with the same
    !  shape as the [[ejecta]] member arrays

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta),                  INTENT(INOUT):: this
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

    END SUBROUTINE interpolate_id_full


    MODULE SUBROUTINE interpolate_id_spacetime &
    ( this, nx, ny, nz, pos, lapse, shift, g, ek )
    !# Stores the spacetime ID in multi-dimensional arrays needed to compute
    !  the BSSN variables and constraints

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta),                        INTENT(INOUT):: this
      INTEGER,                              INTENT(IN)    :: nx
      INTEGER,                              INTENT(IN)    :: ny
      INTEGER,                              INTENT(IN)    :: nz
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(IN)    :: pos
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT(INOUT):: lapse
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(INOUT):: shift
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(INOUT):: g
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(INOUT):: ek

    END SUBROUTINE interpolate_id_spacetime


    MODULE SUBROUTINE interpolate_id_hydro( this, nx, ny, nz, &
                                            pos, &
                                            baryon_density, &
                                            energy_density, &
                                            specific_energy, &
                                            pressure, &
                                            u_euler )
    !# Stores the hydro ID in the arrays needed to compute the constraints
    !  on the refined mesh

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta),                        INTENT(INOUT):: this
      INTEGER,                              INTENT(IN)    :: nx
      INTEGER,                              INTENT(IN)    :: ny
      INTEGER,                              INTENT(IN)    :: nz
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(IN)    :: pos
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT(INOUT):: baryon_density
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT(INOUT):: energy_density
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT(INOUT):: specific_energy
      DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT(INOUT):: pressure
      DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(INOUT):: u_euler

    END SUBROUTINE interpolate_id_hydro


    MODULE SUBROUTINE interpolate_id_particles &
    ( this, n, x, y, z, lapse, shift_x, shift_y, shift_z, &
      g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, baryon_density, energy_density, &
      specific_energy, pressure, u_euler_x, u_euler_y, u_euler_z )
    !! Stores the hydro ID in the arrays needed to compute the SPH ID

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta),                  INTENT(INOUT):: this
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

    END SUBROUTINE interpolate_id_particles


    MODULE SUBROUTINE interpolate_id_mass_b( this, x, y, z, &
                                             g, &
                                             baryon_density, &
                                             gamma_euler )
    !! Stores the hydro ID in the arrays needed to compute the baryon mass

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta),    INTENT(INOUT):: this
      DOUBLE PRECISION, INTENT(IN)    :: x
      DOUBLE PRECISION, INTENT(IN)    :: y
      DOUBLE PRECISION, INTENT(IN)    :: z
      DOUBLE PRECISION, DIMENSION(6), INTENT(OUT):: g
      DOUBLE PRECISION, INTENT(OUT):: baryon_density
      DOUBLE PRECISION, INTENT(OUT):: gamma_euler

    END SUBROUTINE interpolate_id_mass_b


    MODULE SUBROUTINE interpolate_id_k( this, n, x, y, z,&
                                        k_xx, k_xy, k_xz, &
                                        k_yy, k_yz, k_zz )
    !! Stores the components of the extrinsic curvature in arrays

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta),                  INTENT(INOUT):: this
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

    END SUBROUTINE interpolate_id_k


    !
    !-- FUNCTIONS
    !
    MODULE FUNCTION interpolate_mass_density( this, x, y, z ) RESULT( res )
    !! Returns the |lorene| baryon mass density at a point \((x,y,z)\)

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta),    INTENT(IN):: this
      !> \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN), VALUE:: x
      !> \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN), VALUE:: y
      !> \(z\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN), VALUE:: z
      !> Baryon mass density at \((x,y,z)\)
      DOUBLE PRECISION:: res

    END FUNCTION interpolate_mass_density


    MODULE FUNCTION interpolate_pressure( this, x, y, z ) RESULT( res )
    !! Returns the |lorene| baryon mass density at a point \((x,y,z)\)

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta),    INTENT(IN):: this
      !> \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN), VALUE:: x
      !> \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN), VALUE:: y
      !> \(z\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN), VALUE:: z
      !> Pressure at \((x,y,z)\)
      DOUBLE PRECISION:: res

    END FUNCTION interpolate_pressure


    MODULE FUNCTION interpolate_spatial_metric( this, x, y, z ) RESULT( res )
    !# Returns the |lorene| conformally flat spatial metric component
    !  \(g_{xx}=g_{yy}=g_{zz}\) at a point \((x,y,z)\)

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta),    INTENT(IN)       :: this
      !> \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN), VALUE:: x
      !> \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN), VALUE:: y
      !> \(z\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN), VALUE:: z
      !> \(g_{xx}=g_{yy}=g_{zz}\) at \((x,y,z)\)
      DOUBLE PRECISION:: res

    END FUNCTION interpolate_spatial_metric


    MODULE FUNCTION is_hydro_positive( this, x, y, z ) RESULT( res )
    !# Returns `.TRUE.` if the energy density or the specific energy or the
    !  pressure are positivee, `.FALSE.` otherwise

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta),    INTENT(IN)       :: this
      !> \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN), VALUE:: x
      !> \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN), VALUE:: y
      !> \(z\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN), VALUE:: z
      !& `.FALSE.` if the energy density or the specific energy or the pressure
      !  are negative, `.TRUE.` otherwise
      LOGICAL:: res

    END FUNCTION is_hydro_positive


    !-----------------!
    !--  FUNCTIONS  --!
    !-----------------!


    MODULE PURE FUNCTION get_gamma( this )
    !! Interface to [[ejecta_generic:get_gamma]]

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      ! Result
      DOUBLE PRECISION:: get_gamma

    END FUNCTION get_gamma


    MODULE PURE FUNCTION get_kappa( this )
    !! Interface to [[ejecta_generic:get_kappa]]

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      ! Result
      DOUBLE PRECISION:: get_kappa

    END FUNCTION get_kappa


    MODULE FUNCTION get_eos( this, i_matter )
    !! Interface to [[ejecta_generic:get_eos]]

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      INTEGER, INTENT(IN):: i_matter
      !! Index of the matter object whose string is to return
      CHARACTER( LEN= : ), ALLOCATABLE:: get_eos

    END FUNCTION get_eos


    MODULE PURE FUNCTION get_npeos( this )
    !! Interface to [[ejecta_generic:get_npeos]]

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      ! Result
      INTEGER:: get_npeos

    END FUNCTION get_npeos


    MODULE PURE FUNCTION get_gamma0( this )
    !! Interface to [[ejecta_generic:get_gamma0]]

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      ! Result
      DOUBLE PRECISION:: get_gamma0

    END FUNCTION get_gamma0


    MODULE PURE FUNCTION get_gamma1( this )
    !! Interface to [[ejecta_generic:get_gamma1]]

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      ! Result
      DOUBLE PRECISION:: get_gamma1

    END FUNCTION get_gamma1


    MODULE PURE FUNCTION get_gamma2( this )
    !! Interface to [[ejecta_generic:get_gamma2]]

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      ! Result
      DOUBLE PRECISION:: get_gamma2

    END FUNCTION get_gamma2


    MODULE PURE FUNCTION get_gamma3( this )
    !! Interface to [[ejecta_generic:get_gamma3]]

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      ! Result
      DOUBLE PRECISION:: get_gamma3

    END FUNCTION get_gamma3


    MODULE PURE FUNCTION get_kappa0( this )
    !! Interface to [[ejecta_generic:get_kappa0]]

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      ! Result
      DOUBLE PRECISION:: get_kappa0

    END FUNCTION get_kappa0


    MODULE PURE FUNCTION get_kappa1( this )
    !! Interface to [[ejecta_generic:get_kappa1]]

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      ! Result
      DOUBLE PRECISION:: get_kappa1

    END FUNCTION get_kappa1


    MODULE PURE FUNCTION get_kappa2( this )
    !! Interface to [[ejecta_generic:get_kappa2]]

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      ! Result
      DOUBLE PRECISION:: get_kappa2

    END FUNCTION get_kappa2


    MODULE PURE FUNCTION get_kappa3( this )
    !! Interface to [[ejecta_generic:get_kappa3]]

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      ! Result
      DOUBLE PRECISION:: get_kappa3

    END FUNCTION get_kappa3


    MODULE PURE FUNCTION get_logP1( this )
    !! Interface to [[ejecta_generic:get_logP1]]

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      ! Result
      DOUBLE PRECISION:: get_logP1

    END FUNCTION get_logP1


    MODULE PURE FUNCTION get_logRho0( this )
    !! Interface to [[ejecta_generic:get_logRho0]]

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      ! Result
      DOUBLE PRECISION:: get_logRho0

    END FUNCTION get_logRho0


    MODULE PURE FUNCTION get_logRho1( this )
    !! Interface to [[ejecta_generic:get_logRho1]]

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      ! Result
      DOUBLE PRECISION:: get_logRho1

    END FUNCTION get_logRho1


    MODULE PURE FUNCTION get_logRho2( this )
    !! Interface to [[ejecta_generic:get_logRho2]]

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      ! Result
      DOUBLE PRECISION:: get_logRho2

    END FUNCTION get_logRho2


    MODULE FUNCTION get_eos_id( this )
    !! Interface to [[ejecta_generic:get_eos_id]]

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      ! Result
      INTEGER:: get_eos_id

    END FUNCTION get_eos_id


    MODULE SUBROUTINE get_eos_parameters( this, i_matter, eos_params )
    !! Interface to [[ejecta_generic:get_eos_parameters]]

      !> [[ejecta]] object which this PROCEDURE is a member of
      CLASS(ejecta), INTENT(IN):: this
      INTEGER, INTENT(IN):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT):: eos_params
      !# Array containing the parameters of the |eos| for the DRS

    END SUBROUTINE get_eos_parameters


    MODULE SUBROUTINE finalize &
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

    END SUBROUTINE finalize


    MODULE SUBROUTINE nothing( this, flag, switch )
    !# Procedure that does nothing. It is used to instantiate a deferred
    !  idbase procedure which is not needed in TYPE [[ejecta]].
    !  It also serves as a placeholder in case the idbase procedure
    !  will be needed in the future.

      CLASS(ejecta), INTENT(INOUT):: this
      INTEGER,       INTENT(IN)   :: flag
      !! Identifies what kind of initialization has to be done
      LOGICAL,       INTENT(IN), OPTIONAL:: switch
      !! If `.TRUE.`, switch to a different initialization

    END SUBROUTINE nothing


  END INTERFACE


END MODULE ejecta_generic

