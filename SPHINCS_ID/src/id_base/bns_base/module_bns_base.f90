! File:         module_bns_base.f90
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

MODULE bns_base

  !********************************************************
  !
  !# This MODULE contains the definition of TYPE bnsbase,
  !  which is an ABSTRACT TYPE representing any possible
  !  type of binary neutron star (|bns|) initial data (|id|)
  !  to be set up for |sphincsbssn|. That is, |bns| |id|
  !  produced with |lorene|, with |fuka|, etc.
  !
  !  PROCEDURES and variables shared by all the types
  !  of |bns| |id| should belong to TYPE bnsbase, as
  !  they are inherited by its EXTENDED TYPES that
  !  represent more specific types of |bns| |id|.
  !
  !  FT 24.09.2021
  !
  !********************************************************


  USE id_base, ONLY: idbase
  USE utility, ONLY: ios, err_msg


  IMPLICIT NONE


  TYPE surface
    DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE:: points
    !# Array containing the coordinates of the stars' surfaces
    !  The first index runs over the stars; the second and third run over
    !  the surface points (number of points for \(\theta\) and \(\varphi\));
    !  the fourth index runs overthe Cartesian coordinates of the points
  END TYPE surface

  !*******************************************************
  !                                                      *
  !   Definition of TYPE bnsbase (binary neutron star)   *
  !                                                      *
  !*******************************************************

  TYPE, ABSTRACT, EXTENDS(idbase):: bnsbase
  !# ABSTRACT TYPE for |bns| |id| (produced with |lorene|, or
  !  with |fuka|, etc.; or produced with the same tool, but read in different
  !  ways, for example by linking to the |lorene| library, or reading the |id|
  !  from a lattice, etc.)


    !-----------------------------!
    !--  Parameters of the BNS  --!
    !-----------------------------!

    !> Angular velocity \([{\rm rad/s}]\)
    DOUBLE PRECISION:: angular_vel

    !> Distance \(d\) between the points of maximum baryon density \([{\rm km}]\)
    DOUBLE PRECISION:: distance

    !> Distance between the centers of mass \([L_\odot(=1.47662503825040{\rm km})]\)
    DOUBLE PRECISION:: distance_com

    DOUBLE PRECISION, DIMENSION(2):: mass
    !! Array containing the baryonic masses \([M_\odot]\)

    DOUBLE PRECISION:: adm_mass
    !! ADM mass of the system \([M_\odot]\)

    !> Baryonic mass of star 1 \([M_\odot]\)
    DOUBLE PRECISION:: mass1

    !> Baryonic mass of star 2 \([M_\odot]\)
    DOUBLE PRECISION:: mass2

    DOUBLE PRECISION, DIMENSION(2):: mass_grav
    !! Array containing the gravitational masses \([M_\odot]\)

    !> Gravitational mass of star 1 \([M_\odot]\)
    DOUBLE PRECISION:: mass_grav1

    !> Gravitational mass of star 2 \([M_\odot]\)
    DOUBLE PRECISION:: mass_grav2

    !& mOmega= ( [[bnsbase:angular_vel]]\([{\rm km^{-1}}]\) )
    !  \(\times\) ( [[bnsbase:mass_grav1]]\([{\rm km}]\)
    !      + [[bnsbase:mass_grav2]]\([{\rm km}]\) ) [pure number]
    !
    !  Constant used in [K. Hotokezaka et al, Phys. Rev. D 87, 024001](https://arxiv.org/abs/1212.0905){:target="_blank"} (see Sec. IIB) to determine when
    !  the BNS is at approximately 3-4 quasicircular orbits from merger. For the
    ! EOS APR4 and ALF2, this requirement is approximately satisfied for
    ! mOmega \(=0.026\); for the EOS H4 and MS1, for mOmega \(=0.025\).
    DOUBLE PRECISION:: mOmega

    !& Estimated time of the merger \([M_\odot]\)
    !  $$
    !  t_\mathrm{merger}=\dfrac{5}{256}
    !  \dfrac{d^4}{M^1_\mathrm{g}M^2_\mathrm{g}(M^1_\mathrm{g}+M^2_\mathrm{g})}
    !  $$
    !  [P. C. Peters, "Gravitational Radiation and the Motion of Two Point
    !  Masses", Phys. Rev. 136, B1224 (1964)](http://gravity.psu.edu/numrel/jclub/jc/Peters_PR_136_B1224_1964.pdf){:target="_blank"}
    DOUBLE PRECISION:: t_merger

    DOUBLE PRECISION:: linear_momentum_x= 0.0D0
    !! \(x\) component of the ADM linear momentum of the system
    !  \([G M_\odot^2/c]\)
    DOUBLE PRECISION:: linear_momentum_y= 0.0D0
    !! \(y\) component of the ADM linear momentum of the system
    !  \([G M_\odot^2/c]\)
    DOUBLE PRECISION:: linear_momentum_z= 0.0D0
    !! \(z\) component of the ADM linear momentum of the system
    !  \([G M_\odot^2/c]\)

    DOUBLE PRECISION:: angular_momentum_x= 0.0D0
    !! \(x\) component of the angular momentum of the BNS system
    !  \([G M_\odot^2/c]\)
    DOUBLE PRECISION:: angular_momentum_y= 0.0D0
    !! \(y\) component of the angular momentum of the BNS system
    !  \([G M_\odot^2/c]\)
    DOUBLE PRECISION:: angular_momentum_z= 0.0D0
    !! \(z\) component of the angular momentum of the BNS system
    !  \([G M_\odot^2/c]\)

    !& Areal (or circumferential) radius of star 1 \([L_\odot]\)
    ! Note that these is the areal radius of the star in the binary system,
    ! which is different than that of an isolated star. The latter is used
    ! in the mass-radius diagrams, together with the gravitatonal mass
    DOUBLE PRECISION:: area_radius1

    TYPE(surface), DIMENSION(2):: surfaces

    DOUBLE PRECISION, DIMENSION(2,6):: radii
    !# Array containing the **signed** radii of the stars
    !  @todo add details

    !> Radius of star 1, in the x direction, towards the companion \([L_\odot]\)
    DOUBLE PRECISION:: radius1_x_comp

    !> Radius of star 1, in the y direction \([L_\odot]\)
    DOUBLE PRECISION:: radius1_y

    !> Radius of star 1, in the z direction \([L_\odot]\)
    DOUBLE PRECISION:: radius1_z

    !> Radius of star 1, in the x direction, opposite to companion \([L_\odot]\)
    DOUBLE PRECISION:: radius1_x_opp

    !& Stellar center of star 1 (origin of the LORENE chart centered on star 1)
    !  \([L_\odot]\)
    DOUBLE PRECISION:: center1_x

    DOUBLE PRECISION, DIMENSION(2,3):: center
    !# Array containing the centers of the stars
    !  @todo add details

    DOUBLE PRECISION, DIMENSION(2,3):: barycenter
    !# Array containing the barycenters of the stars
    !  @todo add details

    DOUBLE PRECISION, DIMENSION(3):: barycenter_system
    !# Array containing the barycenters of the stars
    !  @todo add details

    !> Barycenter of star 1 \([L_\odot]\)
    DOUBLE PRECISION:: barycenter1_x

    !& Areal (or circumferential) radius of star 2 \([L_\odot]\)
    ! Note that these is the areal radius of the star in the binary system,
    ! which is different than that of an isolated star. The latter is used
    ! in the mass-radius diagrams, together with the gravitatonal mass
    DOUBLE PRECISION:: area_radius2

    !> Radius of star 2, in the x direction, towards the companion \([L_\odot]\)
    DOUBLE PRECISION:: radius2_x_comp

    !> Radius of star 2, in the y direction \([L_\odot]\)
    DOUBLE PRECISION:: radius2_y

    !> Radius of star 2, in the z direction \([L_\odot]\)
    DOUBLE PRECISION:: radius2_z

    !> Radius of star 2, in the x direction, opposite to companion \([L_\odot]\)
    DOUBLE PRECISION:: radius2_x_opp

    !& Stellar center of star 2 (origin of the LORENE chart centered on star 2)
    !  \([L_\odot]\)
    DOUBLE PRECISION:: center2_x

    !> Barycenter of star 2 \([L_\odot]\)
    DOUBLE PRECISION:: barycenter2_x

    !> Central enthalpy for star 1 \([c^2]\)
    DOUBLE PRECISION:: ent_center1

    !> Central baryon number density for star 1 \([L_\odot^{-3}]\)
    DOUBLE PRECISION:: nbar_center1

    !> Central baryon mass density for star 1 \([M_\odot L_\odot^{-3}]\)
    DOUBLE PRECISION:: rho_center1

    !> Central energy density for star 1 \([M_\odot c^2 L_\odot^{-3}]\)
    DOUBLE PRECISION:: energy_density_center1

    !> Central specific energy for star 1 \([c^2]\)
    DOUBLE PRECISION:: specific_energy_center1

    !> Central pressure for star 1 \([M_\odot c^2 L_\odot^{-3}]\)
    DOUBLE PRECISION:: pressure_center1

    !> Central enthalpy for star 2 \([c^2]\)
    DOUBLE PRECISION:: ent_center2

    !> Central baryon number density for star 2 \([L_\odot^{-3}]\)
    DOUBLE PRECISION:: nbar_center2

    !> Central baryon mass density for star 2 \([M_\odot L_\odot^{-3}]\)
    DOUBLE PRECISION:: rho_center2

    !> Central energy density for star 2 \([M_\odot c^2 L_\odot^{-3}]\)
    DOUBLE PRECISION:: energy_density_center2

    !> Central specific energy for star 2 \([c^2]\)
    DOUBLE PRECISION:: specific_energy_center2

    !> Central pressure for star 2 \([M_\odot c^2 L_\odot^{-3}]\)
    DOUBLE PRECISION:: pressure_center2

    !> Name of the equation of state (|eos|) of star 1
    CHARACTER(LEN=:), ALLOCATABLE:: eos1

    !> Name of the equation of state (|eos|) of star 2
    CHARACTER(LEN=:), ALLOCATABLE:: eos2

    INTEGER:: eos1_id
    !! |sphincsid| identifier for the |eos| of star 1

    INTEGER:: eos2_id
    !! |sphincsid| identifier for the |eos| of star 2

    !
    !-- Parameters of single polytropic equations of state for the two NSs
    !

    !> Single polytrope: polytropic index for star 1
    DOUBLE PRECISION:: gamma_1

    !> Single polytrope: polytropic index for star 2
    DOUBLE PRECISION:: gamma_2

    !> Single polytrope: polytropic constant for star 1 [pure number]
    DOUBLE PRECISION:: kappa_1

    !> Single polytrope: polytropic constant for star 2 [pure number]
    DOUBLE PRECISION:: kappa_2

    !
    !-- Parameters of the piecewise polytropic equation of state for NS 1
    !

    !> Piecewise polytrope: Number of polytropic pieces for star 1
    INTEGER:: npeos_1

    !> Piecewise polytrope: polytropic index \(\gamma_0\) for star 1
    DOUBLE PRECISION:: gamma0_1

    !> Piecewise polytrope: polytropic index \(\gamma_1\) for star 1
    DOUBLE PRECISION:: gamma1_1

    !> Piecewise polytrope: polytropic index \(\gamma_2\) for star 1
    DOUBLE PRECISION:: gamma2_1

    !> Piecewise polytrope: polytropic index \(\gamma_3\) for star 1
    DOUBLE PRECISION:: gamma3_1

    !& Piecewise polytrope: polytropic constant \(\kappa_0\) for star 1
    !  [pure number]
    DOUBLE PRECISION:: kappa0_1

    !& Piecewise polytrope: polytropic constant \(\kappa_1\) for star 1
    !  [pure number]
    DOUBLE PRECISION:: kappa1_1

    !& Piecewise polytrope: polytropic constant \(\kappa_2\) for star 1
    !  [pure number]
    DOUBLE PRECISION:: kappa2_1

    !& Piecewise polytrope: polytropic constant \(\kappa_3\) for star 1
    !  [pure number]
    DOUBLE PRECISION:: kappa3_1

    !& Piecewise polytrope: Base 10 exponent of the pressure at the first
    !  fiducial density (between \(\gamma_0\) and \(\gamma_1\)) \([{\rm dyne/cm^2}]\)
    !  for star 1
    DOUBLE PRECISION:: logP1_1

    !& Piecewise polytrope: Base 10 exponent of the first fiducial density
    !  (between \(\gamma_0\) and \(\gamma_1\)) \([{\rm g/cm^3}]\) for star 1
    DOUBLE PRECISION:: logRho0_1

    !& Piecewise polytrope: Base 10 exponent of the second fiducial density
    !  (between \(\gamma_1\) and \(\gamma_2\)) \([{\rm g/cm^3}]\) for star 1
    DOUBLE PRECISION:: logRho1_1

    !& Piecewise polytrope: Base 10 exponent of the third fiducial density
    !  (between \(\gamma_2\) and \(\gamma_3\)) \([{\rm g/cm^3}]\) for star 1
    DOUBLE PRECISION:: logRho2_1

    !
    !-- Parameters of the piecewise polytropic equation of state for NS 2
    !

    !> Piecewise polytrope: Number of polytropic pieces for star 2
    INTEGER:: npeos_2

    !> Piecewise polytrope: polytropic index \(\gamma_0\) for star 2
    DOUBLE PRECISION:: gamma0_2

    !> Piecewise polytrope: polytropic index \(\gamma_1\) for star 2
    DOUBLE PRECISION:: gamma1_2

    !> Piecewise polytrope: polytropic index \(\gamma_2\) for star 2
    DOUBLE PRECISION:: gamma2_2

    !> Piecewise polytrope: polytropic index \(\gamma_3\) for star 2
    DOUBLE PRECISION:: gamma3_2

    !& Piecewise polytrope: polytropic constant \(\kappa_0\) for star 2
    !  [pure number]
    DOUBLE PRECISION:: kappa0_2

    !& Piecewise polytrope: polytropic constant \(\kappa_1\) for star 2
    !  [pure number]
    DOUBLE PRECISION:: kappa1_2

    !& Piecewise polytrope: polytropic constant \(\kappa_2\) for star 2
    !  [pure number]
    DOUBLE PRECISION:: kappa2_2

    !& Piecewise polytrope: polytropic constant \(\kappa_3\) for star 2
    !  [pure number]
    DOUBLE PRECISION:: kappa3_2

    !& Piecewise polytrope: Base 10 exponent of the pressure at the first
    !  fiducial density (between \(\gamma_0\) and \(\gamma_1\)) \([{\rm dyne/cm^2}]\)
    !  for star 2
    DOUBLE PRECISION:: logP1_2

    !& Piecewise polytrope: Base 10 exponent of the second fiducial density
    !  (between \(\gamma_1\) and \(\gamma_2\)) \([{\rm g/cm^3}]\) for star 2
    DOUBLE PRECISION:: logRho0_2

    !& Piecewise polytrope: Base 10 exponent of the second fiducial density
    !  (between \(\gamma_1\) and \(\gamma_2\)) \([{\rm g/cm^3}]\) for star 2
    DOUBLE PRECISION:: logRho1_2

    !& Piecewise polytrope: Base 10 exponent of the third fiducial density
    !  (between \(\gamma_2\) and \(\gamma_3\)) \([{\rm g/cm^3}]\) for star 2
    DOUBLE PRECISION:: logRho2_2


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


    CONTAINS


    !-------------------!
    !--  SUBROUTINES  --!
    !-------------------!

    PROCEDURE(get_eos_id_int), DEFERRED:: get_eos1_id
    !! Returns an integer that identifies the equation of state of star 1

    PROCEDURE(get_eos_id_int), DEFERRED:: get_eos2_id
    !! Returns an integer that identifies the equation of state of star 2

    PROCEDURE:: print_summary               => print_summary_bnsbase

    PROCEDURE(print_summary_derived_int),  DEFERRED:: print_summary_derived
    !# Prints a summary of the physical properties the system
    !  to the standard output and, optionally, to a formatted file whose name
    !  is given as optional argument. Printse information relative to
    !  the derived type oly


    !-----------------!
    !--  FUNCTIONS  --!
    !-----------------!

    PROCEDURE, NON_OVERRIDABLE:: find_center
    !# Finds the center of a matter object, as the point where the
    !  density is maximal.

    PROCEDURE, NON_OVERRIDABLE:: find_radius
    !# Finds the radius of a star, relative to a center and along
    !  a direction. The radius is determined as the first point where the
    !  density is zero.

    PROCEDURE, NON_OVERRIDABLE:: find_surface
    !# Finds the surface of a star, using [[bnsbase::find_radius]] along many
    !  directions.

    PROCEDURE, NON_OVERRIDABLE:: find_print_surfaces
    !# Finds the surfaces of the stars, and prints them to a formatted file.

    !
    !-- FUNCTIONS that access PRIVATE member variables
    !

    !PROCEDURE, PUBLIC:: get_bns_identifier
    !PROCEDURE, PUBLIC:: get_bns_ptr
    PROCEDURE:: return_mass                 => get_mass
    PROCEDURE:: return_adm_mass             => get_adm_mass
    PROCEDURE:: return_center               => get_center
    PROCEDURE:: return_barycenter           => get_barycenter
    PROCEDURE:: return_eos_name             => get_eos
    PROCEDURE:: return_spatial_extent       => get_radii

    PROCEDURE, PUBLIC:: get_angular_vel
    !! Returns [[bnsbase:angular_vel]]
    PROCEDURE, PUBLIC:: get_distance
    !! Returns [[bnsbase:distance]]
    PROCEDURE, PUBLIC:: get_distance_com
    !! Returns [[bnsbase:distance_com]]
    PROCEDURE, PUBLIC:: get_mass1
    !! Returns [[bnsbase:mass1]]
    PROCEDURE, PUBLIC:: get_mass2
    !! Returns [[bnsbase:mass2]]
    PROCEDURE, PUBLIC:: get_grav_mass1
    !! Returns [[bnsbase:mass_grav1]]
    PROCEDURE, PUBLIC:: get_grav_mass2
    !! Returns [[bnsbase:mass_grav2]]
    PROCEDURE, PUBLIC:: get_adm_mass
    !! Returns [[bnsbase:adm_mass]]
    PROCEDURE, PUBLIC:: get_linear_momentum
    !# Returns the linear momentum vector
    !  \((\)[[bnsbase:linear_momentum_x]], [[bnsbase:linear_momentum_y]],
    !  [[bnsbase:linear_momentum_z]]\()\)
    PROCEDURE, PUBLIC:: get_angular_momentum
    !# Returns the angular momentum vector
    !  \((\)[[bnsbase:angular_momentum_x]], [[bnsbase:angular_momentum_y]],
    !  [[bnsbase:angular_momentum_z]]\()\)
    PROCEDURE, PUBLIC:: get_radius1_x_comp
    !! Returns [[bnsbase:radius1_x_comp]]
    PROCEDURE, PUBLIC:: get_radius1_y
    !! Returns [[bnsbase:radius1_y]]
    PROCEDURE, PUBLIC:: get_radius1_z
    !! Returns [[bnsbase:radius1_z]]
    PROCEDURE, PUBLIC:: get_radius1_x_opp
    !! Returns [[bnsbase:radius1_x_opp]]
    PROCEDURE, PUBLIC:: get_center1_x
    !! Returns [[bnsbase:center1_x]]
    PROCEDURE, PUBLIC:: get_barycenter1_x
    !! Returns [[bnsbase:barycenter1_x]]
    PROCEDURE, PUBLIC:: get_radius2_x_comp
    !! Returns [[bnsbase:radius2_x_comp]]
    PROCEDURE, PUBLIC:: get_radius2_y
    !! Returns [[bnsbase:radius2_y]]
    PROCEDURE, PUBLIC:: get_radius2_z
    !! Returns [[bnsbase:radius2_y]]
    PROCEDURE, PUBLIC:: get_radius2_x_opp
    !! Returns [[bnsbase:radius2_x_opp]]
    PROCEDURE, PUBLIC:: get_center2_x
    !! Returns [[bnsbase:center2_x]]
    PROCEDURE, PUBLIC:: get_barycenter2_x
    !! Returns [[bnsbase:barycenter2_x]]
    PROCEDURE, PUBLIC:: get_ent_center1
    !! Returns [[bnsbase:ent_center1]]
    PROCEDURE, PUBLIC:: get_nbar_center1
    !! Returns [[bnsbase:nbar_center1]]
    PROCEDURE, PUBLIC:: get_rho_center1
    !! Returns [[bnsbase:rho_center1]]
    PROCEDURE, PUBLIC:: get_energy_density_center1
    !! Returns [[bnsbase:energy_density_center1]]
    PROCEDURE, PUBLIC:: get_specific_energy_center1
    !! Returns [[bnsbase:specific_energy_center1]]
    PROCEDURE, PUBLIC:: get_pressure_center1
    !! Returns [[bnsbase:pressure_center1]]
    PROCEDURE, PUBLIC:: get_ent_center2
    !! Returns [[bnsbase:ent_center2]]
    PROCEDURE, PUBLIC:: get_nbar_center2
    !! Returns [[bnsbase:nbar_center2]]
    PROCEDURE, PUBLIC:: get_rho_center2
    !! Returns [[bnsbase:rho_center2]]
    PROCEDURE, PUBLIC:: get_energy_density_center2
    !! Returns [[bnsbase:energy_density_center2]]
    PROCEDURE, PUBLIC:: get_specific_energy_center2
    !! Returns [[bnsbase:specific_energy_center2]]
    PROCEDURE, PUBLIC:: get_pressure_center2
    !! Returns [[bnsbase:pressure_center2]]
    PROCEDURE, PUBLIC:: get_eos1
    !! Returns [[bnsbase:eos1]]
    PROCEDURE, PUBLIC:: get_eos2
    !! Returns [[bnsbase:eos2]]
    PROCEDURE, PUBLIC:: get_eos_id
    !! Returns [[bnsbase:eos1_id]] or [[bnsbase:eos2_id]]

    !
    !-- PROCEDURES to be used for single polytropic EOS
    !
    PROCEDURE, PUBLIC:: get_gamma_1
    !! Returns [[bnsbase:gamma_1]]
    PROCEDURE, PUBLIC:: get_gamma_2
    !! Returns [[bnsbase:gamma_2]]
    PROCEDURE, PUBLIC:: get_kappa_1
    !! Returns [[bnsbase:kappa_1]]
    PROCEDURE, PUBLIC:: get_kappa_2
    !! Returns [[bnsbase:kappa_2]]

    !
    !-- PROCEDURES to be used for piecewise polytropic EOS
    !
    PROCEDURE, PUBLIC:: get_npeos_1
    !! Returns [[bnsbase:npeos_1]]
    PROCEDURE, PUBLIC:: get_gamma0_1
    !! Returns [[bnsbase:gamma0_1]]
    PROCEDURE, PUBLIC:: get_gamma1_1
    !! Returns [[bnsbase:gamma1_1]]
    PROCEDURE, PUBLIC:: get_gamma2_1
    !! Returns [[bnsbase:gamma2_1]]
    PROCEDURE, PUBLIC:: get_gamma3_1
    !! Returns [[bnsbase:gamma3_1]]
    PROCEDURE, PUBLIC:: get_kappa0_1
    !! Returns [[bnsbase:kappa0_1]]
    PROCEDURE, PUBLIC:: get_kappa1_1
    !! Returns [[bnsbase:kappa1_1]]
    PROCEDURE, PUBLIC:: get_kappa2_1
    !! Returns [[bnsbase:kappa2_1]]
    PROCEDURE, PUBLIC:: get_kappa3_1
    !! Returns [[bnsbase:kappa3_1]]
    PROCEDURE, PUBLIC:: get_logP1_1
    !! Returns [[bnsbase:logP1_1]]
    PROCEDURE, PUBLIC:: get_logRho0_1
    !! Returns [[bnsbase:logRho0_1]]
    PROCEDURE, PUBLIC:: get_logRho1_1
    !! Returns [[bnsbase:logRho1_1]]
    PROCEDURE, PUBLIC:: get_logRho2_1
    !! Returns [[bnsbase:logRho2_1]]
    PROCEDURE, PUBLIC:: get_npeos_2
    !! Returns [[bnsbase:npeos_2]]
    PROCEDURE, PUBLIC:: get_gamma0_2
    !! Returns [[bnsbase:gamma0_2]]
    PROCEDURE, PUBLIC:: get_gamma1_2
    !! Returns [[bnsbase:gamma1_2]]
    PROCEDURE, PUBLIC:: get_gamma2_2
    !! Returns [[bnsbase:gamma2_2]]
    PROCEDURE, PUBLIC:: get_gamma3_2
    !! Returns [[bnsbase:gamma3_2]]
    PROCEDURE, PUBLIC:: get_kappa0_2
    !! Returns [[bnsbase:kappa0_2]]
    PROCEDURE, PUBLIC:: get_kappa1_2
    !! Returns [[bnsbase:kappa1_2]]
    PROCEDURE, PUBLIC:: get_kappa2_2
    !! Returns [[bnsbase:kappa2_2]]
    PROCEDURE, PUBLIC:: get_kappa3_2
    !! Returns [[bnsbase:kappa3_2]]
    PROCEDURE, PUBLIC:: get_logP1_2
    !! Returns [[bnsbase:logP1_2]]
    PROCEDURE, PUBLIC:: get_logRho0_2
    !! Returns [[bnsbase:logRho0_2]]
    PROCEDURE, PUBLIC:: get_logRho1_2
    !! Returns [[bnsbase:logRho1_2]]
    PROCEDURE, PUBLIC:: get_logRho2_2
    !! Returns [[bnsbase:logRho2_2]]


  END TYPE bnsbase


  ABSTRACT INTERFACE

    SUBROUTINE print_summary_derived_int( this, filename )
    !# Prints a summary of the physical properties the system
    !  to the standard output and, optionally, to a formatted file whose name
    !  is given as the optional argument `filename`

      IMPORT:: bnsbase
      CLASS(bnsbase),   INTENT(IN):: this
      !! [[bnsbase]] object owning this PROCEDURE
      CHARACTER(LEN=*), INTENT(INOUT), OPTIONAL:: filename
      !! Name of the formatted file to print the summary to

    END SUBROUTINE print_summary_derived_int

    FUNCTION get_eos_id_int( this )

      IMPORT:: bnsbase
      CLASS(bnsbase),   INTENT(IN):: this
      !! [[bnsbase]] object owning this PROCEDURE
      INTEGER:: get_eos_id_int
      !! Result

    END FUNCTION get_eos_id_int

  END INTERFACE


  INTERFACE

    !------------------------------!
    !--  OVERRIDING SUBROUTINES  --!
    !------------------------------!


    MODULE SUBROUTINE print_summary_bnsbase( this, filename )
    !# Prints a summary of the physical properties the |bns| system
    !  to the standard output and, optionally, to a formatted file whose name
    !  is given as the optional argument `filename`


      CLASS(bnsbase),   INTENT(IN):: this
      !! [[bnsbase]] object owning this PROCEDURE
      CHARACTER(LEN=*), INTENT(INOUT), OPTIONAL:: filename
      !! Name of the formatted file to print the summary to

    END SUBROUTINE print_summary_bnsbase


    !----------------------------!
    !--  OVERRIDING FUNCTIONS  --!
    !----------------------------!


    MODULE FUNCTION get_mass( this, i_matter )

      CLASS(bnsbase),   INTENT(IN):: this
      !! [[bnsbase]] object owning this PROCEDURE
      INTEGER, INTENT(IN):: i_matter
      !! Index identifying the matter object
      DOUBLE PRECISION:: get_mass
      !! Result

    END FUNCTION get_mass


    MODULE FUNCTION get_radii( this, i_matter )

      CLASS(bnsbase),   INTENT(IN):: this
      !! [[bnsbase]] object owning this PROCEDURE
      INTEGER, INTENT(IN):: i_matter
      !! Index of the matter object whose string is to return
      DOUBLE PRECISION, DIMENSION(6):: get_radii
      !! Result

    END FUNCTION get_radii

    MODULE FUNCTION get_center( this, i_matter )

      CLASS(bnsbase),   INTENT(IN):: this
      !! [[bnsbase]] object owning this PROCEDURE
      INTEGER, INTENT(IN):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION, DIMENSION(3):: get_center
      !! Result

    END FUNCTION get_center


    MODULE FUNCTION get_barycenter( this, i_matter )

      CLASS(bnsbase),   INTENT(IN):: this
      !! [[bnsbase]] object owning this PROCEDURE
      INTEGER, INTENT(IN):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION, DIMENSION(3):: get_barycenter
      !! Result

    END FUNCTION get_barycenter


    MODULE FUNCTION get_eos( this, i_matter )

      CLASS(bnsbase),   INTENT(IN):: this
      !! [[bnsbase]] object owning this PROCEDURE
      INTEGER, INTENT(IN):: i_matter
      !! Index of the matter object whose string is to return
      CHARACTER(LEN=:), ALLOCATABLE:: get_eos
      !! Result

    END FUNCTION get_eos


    !-----------------!
    !--  FUNCTIONS  --!
    !-----------------!


    MODULE FUNCTION find_center( this, separation, x_sign, get_density ) &
       RESULT( center )
    !# Finds the center of a star, as the point where the
    !  density is maximal.

      CLASS(bnsbase),   INTENT(IN):: this
      !! [[bnsbase]] object owning this PROCEDURE
      DOUBLE PRECISION, INTENT(IN):: separation
      !! Separation between the stars
      DOUBLE PRECISION, INTENT(IN):: x_sign
      !! Sign of the x coordinates of the point inside the star
      INTERFACE
        FUNCTION get_density_at_pos(x, y, z) RESULT(rho)
          DOUBLE PRECISION, INTENT(IN):: x
          DOUBLE PRECISION, INTENT(IN):: y
          DOUBLE PRECISION, INTENT(IN):: z
          DOUBLE PRECISION:: rho
        END FUNCTION get_density_at_pos
      END INTERFACE
      PROCEDURE(get_density_at_pos), OPTIONAL:: get_density
      DOUBLE PRECISION                       :: center
      !# Center of a star, as the point where the density is maximal.

    END FUNCTION find_center


    MODULE FUNCTION find_radius( this, center, vector, get_density ) &
      RESULT( radius )
    !# Finds the radius of a matter object, relative to a center and along
    !  a direction. The radius is determined as the first point where the
    !  density is zero.

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object owning this PROCEDURE
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: center
      !! Center point relative to which the radius is measured
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: vector
      !# Vector that defines the direction along which to measure the radius.
      !  If not normalized, it will be normalized.
      INTERFACE
        FUNCTION get_density_at_pos(x, y, z) RESULT(rho)
          DOUBLE PRECISION, INTENT(IN):: x
          DOUBLE PRECISION, INTENT(IN):: y
          DOUBLE PRECISION, INTENT(IN):: z
          DOUBLE PRECISION:: rho
        END FUNCTION get_density_at_pos
      END INTERFACE
      PROCEDURE(get_density_at_pos), OPTIONAL:: get_density
      DOUBLE PRECISION                       :: radius
      !# Radius of the star relative to `center`, along the direction
      !  specified by `vector`

    END FUNCTION find_radius


    MODULE SUBROUTINE find_surface( this, center, n_theta, n_phi, surface )
    !# Finds the surface of a star, using [[bnsbase::find_radius]] along many
    !  directions.

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object owning this PROCEDURE
      DOUBLE PRECISION, DIMENSION(3),                  INTENT(IN) :: center
      !! Center point relative to which the radius is measured
      INTEGER,                                         INTENT(IN) :: n_theta
      !! Number of points in \([0,\pi]\) for the colatitude \(\theta\)
      INTEGER,                                         INTENT(IN) :: n_phi
      !! Number of points in \([0,2\pi]\) for the azimuth \(\varphi\)
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT):: surface
      !# Array storing the coordinates of the points on the surface

    END SUBROUTINE find_surface


    MODULE SUBROUTINE find_print_surfaces( this )
    !# Finds the surfaces of the stars, and prints them to a formatted file.

      CLASS(bnsbase), INTENT(INOUT):: this
      !! [[bnsbase]] object owning this PROCEDURE

    END SUBROUTINE find_print_surfaces


    MODULE FUNCTION get_eos_id( this, i_matter )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object owning this PROCEDURE
      INTEGER,        INTENT(IN):: i_matter
      !! Index of the matter object whose string is to return
      INTEGER:: get_eos_id
      !! Result

    END FUNCTION get_eos_id


    MODULE PURE FUNCTION get_gamma_1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_gamma_1
      !! Result

    END FUNCTION get_gamma_1


    MODULE PURE FUNCTION get_gamma_2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_gamma_2
      !! Result

    END FUNCTION get_gamma_2


    MODULE PURE FUNCTION get_kappa_1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_kappa_1
      !! Result

    END FUNCTION get_kappa_1


    MODULE PURE FUNCTION get_kappa_2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_kappa_2
      !! Result

    END FUNCTION get_kappa_2


    MODULE PURE FUNCTION get_angular_vel( this )
    !! Returns angular_vel

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_angular_vel
      !! Result

    END FUNCTION get_angular_vel


    MODULE PURE FUNCTION get_distance( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_distance
      !! Result

    END FUNCTION get_distance


    MODULE PURE FUNCTION get_distance_com( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_distance_com
      !! Result

    END FUNCTION get_distance_com


    MODULE PURE FUNCTION get_mass1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_mass1
      !! Result

    END FUNCTION get_mass1


    MODULE PURE FUNCTION get_mass2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_mass2
      !! Result

    END FUNCTION get_mass2


    MODULE PURE FUNCTION get_grav_mass1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_grav_mass1
      !! Result

    END FUNCTION get_grav_mass1


    MODULE PURE FUNCTION get_grav_mass2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_grav_mass2
      !! Result

    END FUNCTION get_grav_mass2


    MODULE PURE FUNCTION get_adm_mass( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_adm_mass
      !! Result

    END FUNCTION get_adm_mass


    MODULE PURE FUNCTION get_linear_momentum( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_linear_momentum(3)
      !! Result

    END FUNCTION get_linear_momentum


    MODULE PURE FUNCTION get_angular_momentum( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_angular_momentum(3)
      !! Result

    END FUNCTION get_angular_momentum


    MODULE PURE FUNCTION get_radius1_x_comp( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_radius1_x_comp
      !! Result

    END FUNCTION get_radius1_x_comp


    MODULE PURE FUNCTION get_radius1_y( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_radius1_y
      !! Result

    END FUNCTION get_radius1_y


    MODULE PURE FUNCTION get_radius1_z( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_radius1_z
      !! Result

    END FUNCTION get_radius1_z


    MODULE PURE FUNCTION get_radius1_x_opp( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_radius1_x_opp
      !! Result

    END FUNCTION get_radius1_x_opp


    MODULE PURE FUNCTION get_center1_x( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_center1_x
      !! Result

    END FUNCTION get_center1_x


    MODULE PURE FUNCTION get_barycenter1_x( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_barycenter1_x
      !! Result

    END FUNCTION get_barycenter1_x


    MODULE PURE FUNCTION get_radius2_x_comp( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_radius2_x_comp
      !! Result

    END FUNCTION get_radius2_x_comp


    MODULE PURE FUNCTION get_radius2_y( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_radius2_y
      !! Result

    END FUNCTION get_radius2_y


    MODULE PURE FUNCTION get_radius2_z( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_radius2_z
      !! Result

    END FUNCTION get_radius2_z


    MODULE PURE FUNCTION get_radius2_x_opp( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_radius2_x_opp
      !! Result

    END FUNCTION get_radius2_x_opp


    MODULE PURE FUNCTION get_center2_x( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_center2_x
      !! Result

    END FUNCTION get_center2_x


    MODULE PURE FUNCTION get_barycenter2_x( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_barycenter2_x
      !! Result

    END FUNCTION get_barycenter2_x


    MODULE PURE FUNCTION get_ent_center1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_ent_center1
      !! Result

    END FUNCTION get_ent_center1


    MODULE PURE FUNCTION get_nbar_center1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_nbar_center1
      !! Result

    END FUNCTION get_nbar_center1


    MODULE PURE FUNCTION get_rho_center1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_rho_center1
      !! Result

    END FUNCTION get_rho_center1


    MODULE PURE FUNCTION get_energy_density_center1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_energy_density_center1
      !! Result

    END FUNCTION get_energy_density_center1


    MODULE PURE FUNCTION get_specific_energy_center1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_specific_energy_center1
      !! Result

    END FUNCTION get_specific_energy_center1


    MODULE PURE FUNCTION get_pressure_center1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_pressure_center1
      !! Result

    END FUNCTION get_pressure_center1


    MODULE PURE FUNCTION get_ent_center2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_ent_center2
      !! Result

    END FUNCTION get_ent_center2


    MODULE PURE FUNCTION get_nbar_center2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_nbar_center2
      !! Result

    END FUNCTION get_nbar_center2


    MODULE PURE FUNCTION get_rho_center2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_rho_center2
      !! Result

    END FUNCTION get_rho_center2


    MODULE PURE FUNCTION get_energy_density_center2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_energy_density_center2
      !! Result

    END FUNCTION get_energy_density_center2


    MODULE PURE FUNCTION get_specific_energy_center2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_specific_energy_center2
      !! Result

    END FUNCTION get_specific_energy_center2


    MODULE PURE FUNCTION get_pressure_center2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_pressure_center2
      !! Result

    END FUNCTION get_pressure_center2


    MODULE PURE FUNCTION get_eos1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      CHARACTER(LEN=:), ALLOCATABLE:: get_eos1
      !! Result

    END FUNCTION get_eos1


    MODULE PURE FUNCTION get_eos2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      CHARACTER(LEN=:), ALLOCATABLE:: get_eos2
      !! Result

    END FUNCTION get_eos2


    MODULE PURE FUNCTION get_npeos_1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      INTEGER:: get_npeos_1
      !! Result

    END FUNCTION get_npeos_1


    MODULE PURE FUNCTION get_npeos_2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      INTEGER:: get_npeos_2
      !! Result

    END FUNCTION get_npeos_2


    MODULE PURE FUNCTION get_gamma0_1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_gamma0_1
      !! Result

    END FUNCTION get_gamma0_1


    MODULE PURE FUNCTION get_gamma1_1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_gamma1_1
      !! Result

    END FUNCTION get_gamma1_1


    MODULE PURE FUNCTION get_gamma2_1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_gamma2_1
      !! Result

    END FUNCTION get_gamma2_1


    MODULE PURE FUNCTION get_gamma3_1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_gamma3_1
      !! Result

    END FUNCTION get_gamma3_1


    MODULE PURE FUNCTION get_kappa0_1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_kappa0_1
      !! Result

    END FUNCTION get_kappa0_1


    MODULE PURE FUNCTION get_kappa1_1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_kappa1_1
      !! Result

    END FUNCTION get_kappa1_1


    MODULE PURE FUNCTION get_kappa2_1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_kappa2_1
      !! Result

    END FUNCTION get_kappa2_1


    MODULE PURE FUNCTION get_kappa3_1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_kappa3_1
      !! Result

    END FUNCTION get_kappa3_1


    MODULE PURE FUNCTION get_logP1_1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_logP1_1
      !! Result

    END FUNCTION get_logP1_1


    MODULE PURE FUNCTION get_logRho0_1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_logRho0_1
      !! Result

    END FUNCTION get_logRho0_1


    MODULE PURE FUNCTION get_logRho1_1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_logRho1_1
      !! Result

    END FUNCTION get_logRho1_1


    MODULE PURE FUNCTION get_logRho2_1( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_logRho2_1
      !! Result

    END FUNCTION get_logRho2_1


    MODULE PURE FUNCTION get_gamma0_2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_gamma0_2
      !! Result

    END FUNCTION get_gamma0_2


    MODULE PURE FUNCTION get_gamma1_2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_gamma1_2
      !! Result

    END FUNCTION get_gamma1_2


    MODULE PURE FUNCTION get_gamma2_2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_gamma2_2
      !! Result

    END FUNCTION get_gamma2_2


    MODULE PURE FUNCTION get_gamma3_2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_gamma3_2
      !! Result

    END FUNCTION get_gamma3_2


    MODULE PURE FUNCTION get_kappa0_2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_kappa0_2
      !! Result

    END FUNCTION get_kappa0_2


    MODULE PURE FUNCTION get_kappa1_2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_kappa1_2
      !! Result

    END FUNCTION get_kappa1_2


    MODULE PURE FUNCTION get_kappa2_2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_kappa2_2
      !! Result

    END FUNCTION get_kappa2_2


    MODULE PURE FUNCTION get_kappa3_2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_kappa3_2
      !! Result

    END FUNCTION get_kappa3_2


    MODULE PURE FUNCTION get_logP1_2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_logP1_2
      !! Result

    END FUNCTION get_logP1_2


    MODULE PURE FUNCTION get_logRho0_2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_logRho0_2
      !! Result

    END FUNCTION get_logRho0_2


    MODULE PURE FUNCTION get_logRho1_2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_logRho1_2
      !! Result

    END FUNCTION get_logRho1_2


    MODULE PURE FUNCTION get_logRho2_2( this )

      CLASS(bnsbase), INTENT(IN):: this
      !! [[bnsbase]] object which this PROCEDURE is a member of
      DOUBLE PRECISION:: get_logRho2_2
      !! Result

    END FUNCTION get_logRho2_2

  END INTERFACE

END MODULE bns_base

