! File:         module_diffstar_base.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

MODULE diffstar_base

  !********************************************************
  !
  !# This MODULE contains the definition of TYPE diffstarbase,
  !  which is an ABSTRACT TYPE representing any possible
  !  type of initial data (ID) for a differentially rotating
  !  star (DRS) to be set up for |sphincsbssn|. That is, DRS ID
  !  produced with |lorene|, with |fuka|, etc.
  !
  !  PROCEDURES and variables shared by all the types
  !  of DRS ID should belong to TYPE diffstarbase, as
  !  they are inherited by its EXTENDED TYPES that
  !  represent more specific typesof DRS ID.
  !
  !  FT 22.10.2021
  !
  !********************************************************


  USE id_base, ONLY: idbase
  USE utility, ONLY: ios, err_msg


  IMPLICIT NONE


  !********************************************************************
  !                                                                   *
  !  Definition of TYPE diffstarbase  (differentially rotating star)  *
  !                                                                   *
  !********************************************************************

  TYPE, ABSTRACT, EXTENDS(idbase):: diffstarbase
  !# Base TYPE for DRS ID for |sphincsbssn| (produced with |lorene|, or with
  !  |fuka|, etc.; or produced with the same tool, but read in different ways,
  !  for example by linking to the |lorene| library, or reading the ID from
  !  a lattice, etc.)


    !-----------------------------!
    !--  Parameters of the DRS  --!
    !-----------------------------!

    DOUBLE PRECISION:: omega_c
    !! Central angular velocity \([{\rm rad/s}]\)

    DOUBLE PRECISION:: mass
    !! Baryonic mass of DRS \([M_\odot]\)

    DOUBLE PRECISION:: mass_grav
    !! Gravitational mass of DRS \([M_\odot]\)

    DOUBLE PRECISION:: angular_momentum= 0.0D0
    !! Angular momentum of the DRS \([G M_\odot^2/c]\)

    DOUBLE PRECISION:: tsw
    !# Ratio between the kinetic and gravitatial potential energy \(T/W\).
    !
    !  See Section 6 in
    !  [Gourgoulhon et al, Astron.Astrophys.349:851,1999](https://arxiv.org/abs/astro-ph/9907225v1){:target="_blank"}

    DOUBLE PRECISION:: grv2
    !# Error on the virial identity \({\rm GRV2}\).
    !
    !  See Section 3.5 in
    !  [Gourgoulhon et al, Astron.Astrophys.349:851,1999](https://arxiv.org/abs/astro-ph/9907225v1){:target="_blank"}.

    DOUBLE PRECISION:: grv3
    !# Error on the virial identity \({\rm GRV3}\).
    !
    !  See Section 3.5 in
    !  [Gourgoulhon et al, Astron.Astrophys.349:851,1999](https://arxiv.org/abs/astro-ph/9907225v1){:target="_blank"} .
    !
    !  The error is computed as the integral defined
    !  by Eq. (43) of [Gourgoulhon and Bonazzola, Class. Quantum Grav. 11, 443 (1994)](https://iopscience.iop.org/article/10.1088/0264-9381/11/2/015?pageTitle=IOPscience){:target="_blank"}
    !  divided by the integral of the matter terms.

    DOUBLE PRECISION, DIMENSION(3):: center
    !# Array containing the centers of the stars
    !  @todo add details

    DOUBLE PRECISION, DIMENSION(3):: barycenter
    !# Array containing the barycenters of the stars
    !  @todo add details

    DOUBLE PRECISION, DIMENSION(6):: radii

    DOUBLE PRECISION:: r_circ
    !# Circumferential radius

    DOUBLE PRECISION:: r_mean
    !# Mean radius

    DOUBLE PRECISION:: r_eq
    !# Equatorial radius at \(\phi=0\)

    DOUBLE PRECISION:: r_eq_pi2
    !# Equatorial radius at \(\phi=\dfrac{\pi}{2}\)

    DOUBLE PRECISION:: r_eq_pi
    !# Equatorial radius at \(\phi=\pi\)

    DOUBLE PRECISION:: r_eq_3pi2
    !# Equatorial radius at \(\phi=\dfrac{3\pi}{2}\)

    DOUBLE PRECISION:: r_pole
    !# Polar radius

    DOUBLE PRECISION:: r_ratio
    !# Ratio [[diffstarbase:r_pole]]/[[diffstarbase:r_eq]]

    DOUBLE PRECISION:: r_isco
    !# Radius of the Innermost Stable Circular Orbit (ISCO)

    DOUBLE PRECISION:: f_isco
    !# Orbital frequency of the Innermost Stable Circular Orbit (ISCO)

    DOUBLE PRECISION:: specific_energy_isco
    !# Specific energy of a test particle at the Innermost Stable Circular
    !  Orbit (ISCO)

    DOUBLE PRECISION:: specific_angular_momentum_isco
    !# Specific angular momentum of a test particle at the Innermost Stable
    !  Circular Orbit (ISCO)

    DOUBLE PRECISION:: surface_area
    !# Surface area

    DOUBLE PRECISION:: area_radius
    !# Areal (or circumferential) radius of DRS [Msun_geo]
    !  Note that these is the areal radius of the star in the binary system,
    !  which is different than that of an isolated star. The latter is used
    !  in the mass-radius diagrams, together with the gravitatonal mass

    DOUBLE PRECISION:: ent_center
    !! Central enthalpy \([c^2]\)

    DOUBLE PRECISION:: nbar_center
    !! Central baryon number density \([L_\odot^{-3}]\)

    DOUBLE PRECISION:: rho_center
    !! Central baryon mass density \([M_\odot L_\odot^{-3}]\)

    DOUBLE PRECISION:: energy_density_center
    !! Central energy density \([M_\odot c^2 L_\odot^{-3}]\)

    DOUBLE PRECISION:: specific_energy_center
    !! Central specific energy \([c^2]\)

    DOUBLE PRECISION:: pressure_center
    !! Central pressure \([M_\odot c^2 L_\odot^{-3}]\)

    DOUBLE PRECISION:: redshift_eqf
    !! Forward redshift factor at equator

    DOUBLE PRECISION:: redshift_eqb
    !! Backward redshift factor at equator

    DOUBLE PRECISION:: redshift_pole
    !! Redshift factor at North pole

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

    PROCEDURE(get_eos_id_int), DEFERRED:: get_eos_id
    !! Returns an integer that identifies the equation of state

    !PROCEDURE:: integrate_field_on_star => integrate_baryon_mass_density
    !# Integrates the baryon mass density over the starand computes the radial
    !  mass profile


    !-----------------!
    !--  FUNCTIONS  --!
    !-----------------!

    PROCEDURE:: return_mass                 => get_mass
    PROCEDURE:: return_center               => get_center
    PROCEDURE:: return_barycenter           => get_barycenter
    PROCEDURE:: return_eos_name             => get_eos
    PROCEDURE:: return_spatial_extent       => get_radii


    PROCEDURE, PUBLIC:: get_omega_c
    !! Returns [[diffstarbase:omega_c]]
    !PROCEDURE, PUBLIC:: get_mass
    !! Returns [[diffstarbase:mass]]
    PROCEDURE, PUBLIC:: get_mass_grav
    !! Returns [[diffstarbase:mass_grav]]
    PROCEDURE, PUBLIC:: get_angular_momentum
    !! Returns [[diffstarbase:angular_momentum]]
    PROCEDURE, PUBLIC:: get_tsw
    !! Returns [[diffstarbase:tsw]]
    PROCEDURE, PUBLIC:: get_grv2
    !! Returns [[diffstarbase:grv2]]
    PROCEDURE, PUBLIC:: get_grv3
    !! Returns [[diffstarbase:grv3]]
    PROCEDURE, PUBLIC:: get_r_circ
    !! Returns [[diffstarbase:r_circ]]
    PROCEDURE, PUBLIC:: get_r_mean
    !! Returns [[diffstarbase:r_mean]]
    PROCEDURE, PUBLIC:: get_r_eq
    !! Returns [[diffstarbase:r_eq]]
    PROCEDURE, PUBLIC:: get_r_eq_pi2
    !! Returns [[diffstarbase:r_eq_pi2]]
    PROCEDURE, PUBLIC:: get_r_eq_pi
    !! Returns [[diffstarbase:r_eq_pi]]
    PROCEDURE, PUBLIC:: get_r_eq_3pi2
    !! Returns [[diffstarbase:r_eq_3pi2]]
    PROCEDURE, PUBLIC:: get_r_pole
    !! Returns [[diffstarbase:r_pole]]
    PROCEDURE, PUBLIC:: get_r_ratio
    !! Returns [[diffstarbase:r_ratio]]
    PROCEDURE, PUBLIC:: get_r_isco
    !! Returns [[diffstarbase:r_isco]]
    PROCEDURE, PUBLIC:: get_f_isco
    !! Returns [[diffstarbase:f_isco]]
    PROCEDURE, PUBLIC:: get_specific_energy_isco
    !! Returns [[diffstarbase:specific_energy_isco]]
    PROCEDURE, PUBLIC:: get_specific_angular_momentum_isco
    !! Returns [[diffstarbase:specific_angular_momentum_isco]]
    PROCEDURE, PUBLIC:: get_surface_area
    !! Returns [[diffstarbase:surface_area]]
    PROCEDURE, PUBLIC:: get_area_radius
    !! Returns [[diffstarbase:area_radius]]
    PROCEDURE, PUBLIC:: get_ent_center
    !! Returns [[diffstarbase:ent_center]]
    PROCEDURE, PUBLIC:: get_nbar_center
    !! Returns [[diffstarbase:nbar_center]]
    PROCEDURE, PUBLIC:: get_rho_center
    !! Returns [[diffstarbase:rho_center]]
    PROCEDURE, PUBLIC:: get_energy_density_center
    !! Returns [[diffstarbase:energy_density_center]]
    PROCEDURE, PUBLIC:: get_specific_energy_center
    !! Returns [[diffstarbase:specific_energy_center]]
    PROCEDURE, PUBLIC:: get_pressure_center
    !! Returns [[diffstarbase:pressure_center]]
    !PROCEDURE, PUBLIC:: get_eos
    !! Returns [[diffstarbase:eos]]

    !
    !-- PROCEDURES to be used for single polytropic EOS
    !
    PROCEDURE, PUBLIC:: get_gamma
    !! Returns [[diffstarbase:gamma]]
    PROCEDURE, PUBLIC:: get_kappa
    !! Returns [[diffstarbase:kappa]]

    !
    !-- PROCEDURES to be used for piecewise polytropic EOS
    !
    PROCEDURE, PUBLIC:: get_npeos
    !! Returns [[diffstarbase:npeos]]
    PROCEDURE, PUBLIC:: get_gamma0
    !! Returns [[diffstarbase:gamma0]]
    PROCEDURE, PUBLIC:: get_gamma1
    !! Returns [[diffstarbase:gamma1]]
    PROCEDURE, PUBLIC:: get_gamma2
    !! Returns [[diffstarbase:gamma2]]
    PROCEDURE, PUBLIC:: get_gamma3
    !! Returns [[diffstarbase:gamma3]]
    PROCEDURE, PUBLIC:: get_kappa0
    !! Returns [[diffstarbase:kappa0]]
    PROCEDURE, PUBLIC:: get_kappa1
    !! Returns [[diffstarbase:kappa1]]
    PROCEDURE, PUBLIC:: get_kappa2
    !! Returns [[diffstarbase:kappa2]]
    PROCEDURE, PUBLIC:: get_kappa3
    !! Returns [[diffstarbase:kappa3]]
    PROCEDURE, PUBLIC:: get_logP1
    !! Returns [[diffstarbase:logP1]]
    PROCEDURE, PUBLIC:: get_logRho0
    !! Returns [[diffstarbase:logRho0]]
    PROCEDURE, PUBLIC:: get_logRho1
    !! Returns [[diffstarbase:logRho1]]
    PROCEDURE, PUBLIC:: get_logRho2
    !! Returns [[diffstarbase:logRho2]]


  END TYPE diffstarbase


  ABSTRACT INTERFACE

    FUNCTION get_eos_id_int( THIS )

      IMPORT:: diffstarbase
      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      INTEGER:: get_eos_id_int

    END FUNCTION get_eos_id_int

  END INTERFACE


  INTERFACE


 !   MODULE SUBROUTINE integrate_baryon_mass_density( THIS, center, radius, &
 !                                                    central_density, &
 !                                                    dr, dth, dphi, &
 !                                                    mass, mass_profile, &
 !                                                    mass_profile_idx )
 !   !# Integrates the |lorene| baryon mass density to compute the radial mass
 !   !  profile. TODO: Improve integration algorithm.
 !
 !     !> [[diffstarbase]] object which this PROCEDURE is a member of
 !     CLASS(diffstarbase), INTENT( IN OUT )      :: THIS
 !     !& Array to store the indices for array mass_profile, sorted so that
 !     !  mass_profile[mass_profile_idx] is in increasing order
 !     INTEGER, DIMENSION(:), ALLOCATABLE, INTENT( IN OUT ):: mass_profile_idx
 !     !> Center of the star
 !     DOUBLE PRECISION, INTENT( IN )    :: center
 !     !> Central density of the star
 !     DOUBLE PRECISION, INTENT( IN )    :: central_density
 !     !> Radius of the star
 !     DOUBLE PRECISION, INTENT( IN )    :: radius
 !     !> Integration steps
 !     DOUBLE PRECISION, INTENT( IN )    :: dr, dth, dphi
 !     !> Integrated mass of the star
 !     DOUBLE PRECISION, INTENT( IN OUT ):: mass
 !     !> Array storing the radial mass profile of the star
 !     DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT( IN OUT ):: &
 !                                      mass_profile
 !
 !   END SUBROUTINE integrate_baryon_mass_density


    !----------------------------!
    !--  OVERRIDING FUNCTIONS  --!
    !----------------------------!


    MODULE FUNCTION get_mass( THIS, i_matter )
    !! Returns [[diffstarbase:mass]]

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN OUT ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      ! Result
      DOUBLE PRECISION:: get_mass

    END FUNCTION get_mass


    MODULE FUNCTION get_center( THIS, i_matter )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN OUT ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION, DIMENSION(3):: get_center

    END FUNCTION get_center


    MODULE FUNCTION get_barycenter( THIS, i_matter )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN OUT ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose parameter is to return
      DOUBLE PRECISION, DIMENSION(3):: get_barycenter

    END FUNCTION get_barycenter


    MODULE FUNCTION get_radii( THIS, i_matter )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN OUT ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose string is to return
      DOUBLE PRECISION, DIMENSION(6):: get_radii

    END FUNCTION get_radii


    MODULE FUNCTION get_eos( THIS, i_matter )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN OUT ):: THIS
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object whose string is to return
      CHARACTER( LEN= : ), ALLOCATABLE:: get_eos

    END FUNCTION get_eos


    !-----------------!
    !--  FUNCTIONS  --!
    !-----------------!


    MODULE FUNCTION get_gamma( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma

    END FUNCTION get_gamma


    MODULE FUNCTION get_kappa( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa

    END FUNCTION get_kappa


    MODULE FUNCTION get_omega_c( THIS )
    !! Returns [[diffstarbase:omega_c]]

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_omega_c

    END FUNCTION get_omega_c


    MODULE FUNCTION get_mass_grav( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_mass_grav

    END FUNCTION get_mass_grav


    MODULE FUNCTION get_angular_momentum( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_angular_momentum

    END FUNCTION get_angular_momentum


    MODULE FUNCTION get_tsw( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_tsw

    END FUNCTION get_tsw


    MODULE FUNCTION get_grv2( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_grv2

    END FUNCTION get_grv2


    MODULE FUNCTION get_grv3( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_grv3

    END FUNCTION get_grv3


    MODULE FUNCTION get_r_circ( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_r_circ

    END FUNCTION get_r_circ


    MODULE FUNCTION get_r_mean( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_r_mean

    END FUNCTION get_r_mean


    MODULE FUNCTION get_r_eq( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_r_eq

    END FUNCTION get_r_eq


    MODULE FUNCTION get_r_eq_pi2( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_r_eq_pi2

    END FUNCTION get_r_eq_pi2


    MODULE FUNCTION get_r_eq_pi( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_r_eq_pi

    END FUNCTION get_r_eq_pi


    MODULE FUNCTION get_r_eq_3pi2( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_r_eq_3pi2

    END FUNCTION get_r_eq_3pi2


    MODULE FUNCTION get_r_pole( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_r_pole

    END FUNCTION get_r_pole


    MODULE FUNCTION get_r_ratio( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_r_ratio

    END FUNCTION get_r_ratio


    MODULE FUNCTION get_r_isco( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_r_isco

    END FUNCTION get_r_isco


    MODULE FUNCTION get_f_isco( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_f_isco

    END FUNCTION get_f_isco


    MODULE FUNCTION get_specific_energy_isco( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_specific_energy_isco

    END FUNCTION get_specific_energy_isco


    MODULE FUNCTION get_specific_angular_momentum_isco( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_specific_angular_momentum_isco

    END FUNCTION get_specific_angular_momentum_isco


    MODULE FUNCTION get_surface_area( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_surface_area

    END FUNCTION get_surface_area


    MODULE FUNCTION get_area_radius( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_area_radius

    END FUNCTION get_area_radius


    MODULE FUNCTION get_ent_center( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_ent_center

    END FUNCTION get_ent_center


    MODULE FUNCTION get_nbar_center( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_nbar_center

    END FUNCTION get_nbar_center


    MODULE FUNCTION get_rho_center( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_rho_center

    END FUNCTION get_rho_center


    MODULE FUNCTION get_energy_density_center( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_energy_density_center

    END FUNCTION get_energy_density_center


    MODULE FUNCTION get_specific_energy_center( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_specific_energy_center

    END FUNCTION get_specific_energy_center


    MODULE FUNCTION get_pressure_center( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_pressure_center

    END FUNCTION get_pressure_center


    MODULE FUNCTION get_npeos( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      INTEGER:: get_npeos

    END FUNCTION get_npeos


    MODULE FUNCTION get_gamma0( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma0

    END FUNCTION get_gamma0


    MODULE FUNCTION get_gamma1( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma1

    END FUNCTION get_gamma1


    MODULE FUNCTION get_gamma2( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma2

    END FUNCTION get_gamma2


    MODULE FUNCTION get_gamma3( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma3

    END FUNCTION get_gamma3


    MODULE FUNCTION get_kappa0( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa0

    END FUNCTION get_kappa0


    MODULE FUNCTION get_kappa1( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa1

    END FUNCTION get_kappa1


    MODULE FUNCTION get_kappa2( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa2

    END FUNCTION get_kappa2


    MODULE FUNCTION get_kappa3( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa3

    END FUNCTION get_kappa3


    MODULE FUNCTION get_logP1( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logP1

    END FUNCTION get_logP1


    MODULE FUNCTION get_logRho0( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logRho0

    END FUNCTION get_logRho0


    MODULE FUNCTION get_logRho1( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logRho1

    END FUNCTION get_logRho1


    MODULE FUNCTION get_logRho2( THIS )

      !> [[diffstarbase]] object which this PROCEDURE is a member of
      CLASS(diffstarbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logRho2

    END FUNCTION get_logRho2

  END INTERFACE


END MODULE diffstar_base

