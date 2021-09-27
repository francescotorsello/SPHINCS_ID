! File:         module_bns_base.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

MODULE bns_base

  !*******************************************************
  !
  !#
  !
  !   FT 24.09.2021
  !
  !*******************************************************


  USE id_base, ONLY: idbase


  IMPLICIT NONE


  !*******************************************************
  !                                                      *
  !     Definition of TYPE bns (binary neutron star)     *
  !                                                      *
  !   This class imports and stores the LORENE BNS ID    *
  !                                                      *
  !*******************************************************

  TYPE, ABSTRACT, EXTENDS(idbase):: bnsbase


    !> LORENE identifiers for the EoS
    INTEGER:: eos1_id, eos2_id

    !
    !-- Parameters of the binary system
    !

    !> Angular velocity [rad/s]
    DOUBLE PRECISION:: angular_vel
    !> Distance \(d\) between the points of maximum baryon density [km]
    DOUBLE PRECISION:: distance
    !> Distance between the centers of mass [MSun_geo(=1.47662503825040km)]
    DOUBLE PRECISION:: distance_com
    !> Baryonic mass of star 1 [MSun]
    DOUBLE PRECISION:: mass1
    !> Baryonic mass of star 2 [MSun]
    DOUBLE PRECISION:: mass2
    !> Gravitational mass of star 1 [MSun]
    DOUBLE PRECISION:: mass_grav1
    !> Gravitational mass of star 2 [MSun]
    DOUBLE PRECISION:: mass_grav2
    !> ADM mass of the BNS [Msun]
    DOUBLE PRECISION:: adm_mass
    !& mOmega= ( angular_vel[km^{-1}] )*( mass_grav1[km] + mass_grav2[km] )
    !  [pure number]
    DOUBLE PRECISION:: mOmega
    !& Estimated time of the merger [MSun]
    !  $$
    !  t_\mathrm{merger}=\dfrac{5}{256}
    !  \dfrac{d^4}{M^1_\mathrm{g}M^2_\mathrm{g}(M^1_\mathrm{g}+M^2_\mathrm{g})}
    !  $$
    !  P. C. Peters, "Gravitational Radiation and the Motion of Two Point
    !  Masses", Phys. Rev. 136, B1224 (1964)
    !  http://gravity.psu.edu/numrel/jclub/jc/Peters_PR_136_B1224_1964.pdf
    DOUBLE PRECISION:: t_merger
    !> Angular momentum of the BNS system [G Msun^2/c]
    DOUBLE PRECISION:: angular_momentum= 0.0D0
    !& Areal (or circumferential) radius of star 1 [Msun_geo]
    ! Note that these is the areal radius of the star in the binary system,
    ! which is different than that of an isolated star. The latter is used
    ! in the mass-radius diagrams, together with the gravitatonal mass
    DOUBLE PRECISION:: area_radius1
    !> Radius of star 1, in the x direction, towards the companion [Msun_geo]
    DOUBLE PRECISION:: radius1_x_comp
    !> Radius of star 1, in the y direction [Msun_geo]
    DOUBLE PRECISION:: radius1_y
    !> Radius of star 1, in the z direction [Msun_geo]
    DOUBLE PRECISION:: radius1_z
    !> Radius of star 1, in the x direction, opposite to companion [Msun_geo]
    DOUBLE PRECISION:: radius1_x_opp
    !& Stellar center of star 1 (origin of the LORENE chart centered on star 1)
    !  [Msun_geo]
    DOUBLE PRECISION:: center1_x
    !> Barycenter of star 1 [Msun_geo]
    DOUBLE PRECISION:: barycenter1_x
    !& Areal (or circumferential) radius of star 2 [Msun_geo]
    ! Note that these is the areal radius of the star in the binary system,
    ! which is different than that of an isolated star. The latter is used
    ! in the mass-radius diagrams, together with the gravitatonal mass
    DOUBLE PRECISION:: area_radius2
    !> Radius of star 2, in the x direction, towards the companion [Msun_geo]
    DOUBLE PRECISION:: radius2_x_comp
    !> Radius of star 2, in the y direction [Msun_geo]
    DOUBLE PRECISION:: radius2_y
    !> Radius of star 2, in the z direction [Msun_geo]
    DOUBLE PRECISION:: radius2_z
    !> Radius of star 2, in the x direction, opposite to companion [Msun_geo]
    DOUBLE PRECISION:: radius2_x_opp
    !& Stellar center of star 2 (origin of the LORENE chart centered on star 2)
    !  [Msun_geo]
    DOUBLE PRECISION:: center2_x
    !> Barycenter of star 2 [Msun_geo]
    DOUBLE PRECISION:: barycenter2_x
    !> Central enthalpy for star 1 [c^2]
    DOUBLE PRECISION:: ent_center1 ;
    !> Central baryon number density for star 1 [Msun_geo^-3]
    DOUBLE PRECISION:: nbar_center1 ;
    !> Central baryon mass density for star 1 [Msun Msun_geo^-3]
    DOUBLE PRECISION:: rho_center1 ;
    !> Central energy density for star 1 [Msun c^2 Msun_geo^-3]
    DOUBLE PRECISION:: energy_density_center1 ;
    !> Central specific energy for star 1 [c^2]
    DOUBLE PRECISION:: specific_energy_center1 ;
    !> Central pressure for star 1 [Msun c^2 Msun_geo^-3]
    DOUBLE PRECISION:: pressure_center1 ;
    !> Central enthalpy for star 2 [c^2]
    DOUBLE PRECISION:: ent_center2 ;
    !> Central baryon number density for star 2 [Msun_geo^-3]
    DOUBLE PRECISION:: nbar_center2 ;
    !> Central baryon mass density for star 2 [Msun Msun_geo^-3]
    DOUBLE PRECISION:: rho_center2 ;
    !> Central energy density for star 2 [Msun c^2 Msun_geo^-3]
    DOUBLE PRECISION:: energy_density_center2 ;
    !> Central specific energy for star 2 [c^2]
    DOUBLE PRECISION:: specific_energy_center2 ;
    !> Central pressure for star 2 [Msun c^2 Msun_geo^-3]
    DOUBLE PRECISION:: pressure_center2 ;
    !> LORENE name of the equation of state (EoS) of star 1
    CHARACTER( LEN=: ), ALLOCATABLE:: eos1
    !> LORENE name of the equation of state (EoS) of star 2
    CHARACTER( LEN=: ), ALLOCATABLE:: eos2

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
    !  fiducial density (between \(\gamma_0\) and \(\gamma_1\)) [dyne/cm^2]
    !  for star 1
    DOUBLE PRECISION:: logP1_1
    !& Piecewise polytrope: Base 10 exponent of the first fiducial density
    !  (between \(\gamma_0\) and \(\gamma_1\)) [g/cm^3] for star 1
    DOUBLE PRECISION:: logRho0_1
    !& Piecewise polytrope: Base 10 exponent of the second fiducial density
    !  (between \(\gamma_1\) and \(\gamma_2\)) [g/cm^3] for star 1
    DOUBLE PRECISION:: logRho1_1
    !& Piecewise polytrope: Base 10 exponent of the third fiducial density
    !  (between \(\gamma_2\) and \(\gamma_3\)) [g/cm^3] for star 1
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
    !  fiducial density (between \(\gamma_0\) and \(\gamma_1\)) [dyne/cm^2]
    !  for star 2
    DOUBLE PRECISION:: logP1_2
    !& Piecewise polytrope: Base 10 exponent of the second fiducial density
    !  (between \(\gamma_1\) and \(\gamma_2\)) [g/cm^3] for star 2
    DOUBLE PRECISION:: logRho0_2
    !& Piecewise polytrope: Base 10 exponent of the second fiducial density
    !  (between \(\gamma_1\) and \(\gamma_2\)) [g/cm^3] for star 2
    DOUBLE PRECISION:: logRho1_2
    !& Piecewise polytrope: Base 10 exponent of the third fiducial density
    !  (between \(\gamma_2\) and \(\gamma_3\)) [g/cm^3] for star 2
    DOUBLE PRECISION:: logRho2_2



    CONTAINS



      !
      !-- Overloaded FUNCTION to access the fields as arrays and as values
      !

 !     GENERIC, PUBLIC:: get_field => get_fa, get_fv
 !     !# GENERIC PROCEDURE, overloded to access the bns member variables as arrays
 !     !  and as values
 !     PROCEDURE::       get_fa    => get_field_array
 !     !! Access the bns member arrays
 !     PROCEDURE::       get_fv    => get_field_value
      !! Access the components of the bns member arrays

      !
      !-- FUNCTIONS that access member variables
      !

      !PROCEDURE, PUBLIC:: get_bns_identifier
      !PROCEDURE, PUBLIC:: get_bns_ptr

      PROCEDURE, PUBLIC:: get_angular_vel
      !! Returns [[bns:angular_vel]]
      PROCEDURE, PUBLIC:: get_distance
      !! Returns [[bns:distance]]
      PROCEDURE, PUBLIC:: get_distance_com
      !! Returns [[bns:distance_com]]
      PROCEDURE, PUBLIC:: get_mass1
      !! Returns [[bns:mass1]]
      PROCEDURE, PUBLIC:: get_mass2
      !! Returns [[bns:mass2]]
      PROCEDURE, PUBLIC:: get_grav_mass1
      !! Returns [[bns:mass_grav1]]
      PROCEDURE, PUBLIC:: get_grav_mass2
      !! Returns [[bns:mass_grav2]]
      PROCEDURE, PUBLIC:: get_adm_mass
      !! Returns [[bns:adm_mass]]
      PROCEDURE, PUBLIC:: get_angular_momentum
      !! Returns [[bns:angular_momentum]]
      PROCEDURE, PUBLIC:: get_radius1_x_comp
      !! Returns [[bns:radius1_x_comp]]
      PROCEDURE, PUBLIC:: get_radius1_y
      !! Returns [[bns:radius1_y]]
      PROCEDURE, PUBLIC:: get_radius1_z
      !! Returns [[bns:radius1_z]]
      PROCEDURE, PUBLIC:: get_radius1_x_opp
      !! Returns [[bns:radius1_x_opp]]
      PROCEDURE, PUBLIC:: get_center1_x
      !! Returns [[bns:center1_x]]
      PROCEDURE, PUBLIC:: get_barycenter1_x
      !! Returns [[bns:barycenter1_x]]
      PROCEDURE, PUBLIC:: get_radius2_x_comp
      !! Returns [[bns:radius2_x_comp]]
      PROCEDURE, PUBLIC:: get_radius2_y
      !! Returns [[bns:radius2_y]]
      PROCEDURE, PUBLIC:: get_radius2_z
      !! Returns [[bns:radius2_y]]
      PROCEDURE, PUBLIC:: get_radius2_x_opp
      !! Returns [[bns:radius2_x_opp]]
      PROCEDURE, PUBLIC:: get_center2_x
      !! Returns [[bns:center2_x]]
      PROCEDURE, PUBLIC:: get_barycenter2_x
      !! Returns [[bns:barycenter2_x]]
      PROCEDURE, PUBLIC:: get_ent_center1
      !! Returns [[bns:ent_center1]]
      PROCEDURE, PUBLIC:: get_nbar_center1
      !! Returns [[bns:nbar_center1]]
      PROCEDURE, PUBLIC:: get_rho_center1
      !! Returns [[bns:rho_center1]]
      PROCEDURE, PUBLIC:: get_energy_density_center1
      !! Returns [[bns:energy_density_center1]]
      PROCEDURE, PUBLIC:: get_specific_energy_center1
      !! Returns [[bns:specific_energy_center1]]
      PROCEDURE, PUBLIC:: get_pressure_center1
      !! Returns [[bns:pressure_center1]]
      PROCEDURE, PUBLIC:: get_ent_center2
      !! Returns [[bns:ent_center2]]
      PROCEDURE, PUBLIC:: get_nbar_center2
      !! Returns [[bns:nbar_center2]]
      PROCEDURE, PUBLIC:: get_rho_center2
      !! Returns [[bns:rho_center2]]
      PROCEDURE, PUBLIC:: get_energy_density_center2
      !! Returns [[bns:energy_density_center2]]
      PROCEDURE, PUBLIC:: get_specific_energy_center2
      !! Returns [[bns:specific_energy_center2]]
      PROCEDURE, PUBLIC:: get_pressure_center2
      !! Returns [[bns:pressure_center2]]
      PROCEDURE, PUBLIC:: get_eos1
      !! Returns [[bns:eos1]]
      PROCEDURE, PUBLIC:: get_eos2
      !! Returns [[bns:eos2]]
      PROCEDURE, PUBLIC:: get_eos1_id
      !! Returns [[bns:eos1_id]]
      PROCEDURE, PUBLIC:: get_eos2_id
      !! Returns [[bns:eos2_id]]

      !
      !-- PROCEDURES to be used for single polytropic EOS
      !
      PROCEDURE, PUBLIC:: get_gamma_1
      !! Returns [[bns:gamma_1]]
      PROCEDURE, PUBLIC:: get_gamma_2
      !! Returns [[bns:gamma_2]]
      PROCEDURE, PUBLIC:: get_kappa_1
      !! Returns [[bns:kappa_1]]
      PROCEDURE, PUBLIC:: get_kappa_2
      !! Returns [[bns:kappa_2]]

      !
      !-- PROCEDURES to be used for piecewise polytropic EOS
      !
      PROCEDURE, PUBLIC:: get_npeos_1
      !! Returns [[bns:npeos_1]]
      PROCEDURE, PUBLIC:: get_gamma0_1
      !! Returns [[bns:gamma0_1]]
      PROCEDURE, PUBLIC:: get_gamma1_1
      !! Returns [[bns:gamma1_1]]
      PROCEDURE, PUBLIC:: get_gamma2_1
      !! Returns [[bns:gamma2_1]]
      PROCEDURE, PUBLIC:: get_gamma3_1
      !! Returns [[bns:gamma3_1]]
      PROCEDURE, PUBLIC:: get_kappa0_1
      !! Returns [[bns:kappa0_1]]
      PROCEDURE, PUBLIC:: get_kappa1_1
      !! Returns [[bns:kappa1_1]]
      PROCEDURE, PUBLIC:: get_kappa2_1
      !! Returns [[bns:kappa2_1]]
      PROCEDURE, PUBLIC:: get_kappa3_1
      !! Returns [[bns:kappa3_1]]
      PROCEDURE, PUBLIC:: get_logP1_1
      !! Returns [[bns:logP1_1]]
      PROCEDURE, PUBLIC:: get_logRho0_1
      !! Returns [[bns:logRho0_1]]
      PROCEDURE, PUBLIC:: get_logRho1_1
      !! Returns [[bns:logRho1_1]]
      PROCEDURE, PUBLIC:: get_logRho2_1
      !! Returns [[bns:logRho2_1]]
      PROCEDURE, PUBLIC:: get_npeos_2
      !! Returns [[bns:npeos_2]]
      PROCEDURE, PUBLIC:: get_gamma0_2
      !! Returns [[bns:gamma0_2]]
      PROCEDURE, PUBLIC:: get_gamma1_2
      !! Returns [[bns:gamma1_2]]
      PROCEDURE, PUBLIC:: get_gamma2_2
      !! Returns [[bns:gamma2_2]]
      PROCEDURE, PUBLIC:: get_gamma3_2
      !! Returns [[bns:gamma3_2]]
      PROCEDURE, PUBLIC:: get_kappa0_2
      !! Returns [[bns:kappa0_2]]
      PROCEDURE, PUBLIC:: get_kappa1_2
      !! Returns [[bns:kappa1_2]]
      PROCEDURE, PUBLIC:: get_kappa2_2
      !! Returns [[bns:kappa2_2]]
      PROCEDURE, PUBLIC:: get_kappa3_2
      !! Returns [[bns:kappa3_2]]
      PROCEDURE, PUBLIC:: get_logP1_2
      !! Returns [[bns:logP1_2]]
      PROCEDURE, PUBLIC:: get_logRho0_2
      !! Returns [[bns:logRho0_2]]
      PROCEDURE, PUBLIC:: get_logRho1_2
      !! Returns [[bns:logRho1_2]]
      PROCEDURE, PUBLIC:: get_logRho2_2
      !! Returns [[bns:logRho2_2]]


  END TYPE bnsbase


  INTERFACE

  !  MODULE FUNCTION get_field_array( THIS, field ) RESULT( field_array )
  !  !! Returns the [[bns]] member arrays named field
  !
  !    !> [[bns]] object which this PROCEDURE is a member of
  !    CLASS(bnsbase),          INTENT( IN )             :: THIS
  !    !> Name of the desired [[bns]] member array
  !    CHARACTER( LEN= : ), INTENT( IN ), ALLOCATABLE:: field
  !    !> Desired [[bns]] member array
  !    DOUBLE PRECISION, DIMENSION(:),    ALLOCATABLE:: field_array
  !
  !  END FUNCTION get_field_array
  !
  !
  !  MODULE FUNCTION get_field_value( THIS, field, n ) RESULT( field_value )
  !  !! Returns the component n of the [[bns]] member arrays named field
  !
  !    !> [[bns]] object which this PROCEDURE is a member of
  !    CLASS(bnsbase),          INTENT( IN )             :: THIS
  !    !> Name of the desired [[bns]] member array
  !    CHARACTER( LEN= : ), INTENT( IN ), ALLOCATABLE:: field
  !    !> Component of the desired [[bns]] member array
  !    INTEGER,             INTENT( IN )             :: n
  !    !> Component n of the desired [[bns]] member array
  !    DOUBLE PRECISION                              :: field_value
  !
  !  END FUNCTION get_field_value
  !
  !
  !  MODULE FUNCTION get_bns_identifier( THIS )
  !
  !    !> [[bns]] object which this PROCEDURE is a member of
  !    CLASS(bnsbase), INTENT( IN ):: THIS
  !    ! Result
  !    DOUBLE PRECISION:: get_bns_identifier
  !
  !  END FUNCTION get_bns_identifier


    MODULE FUNCTION get_gamma_1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma_1

    END FUNCTION get_gamma_1


    MODULE FUNCTION get_gamma_2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma_2

    END FUNCTION get_gamma_2


    MODULE FUNCTION get_kappa_1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa_1

    END FUNCTION get_kappa_1


    MODULE FUNCTION get_kappa_2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa_2

    END FUNCTION get_kappa_2


    MODULE FUNCTION get_angular_vel( THIS )
    !! Returns angular_vel

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_angular_vel

    END FUNCTION get_angular_vel


    MODULE FUNCTION get_distance( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_distance

    END FUNCTION get_distance


    MODULE FUNCTION get_distance_com( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_distance_com

    END FUNCTION get_distance_com


    MODULE FUNCTION get_mass1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_mass1

    END FUNCTION get_mass1


    MODULE FUNCTION get_mass2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_mass2

    END FUNCTION get_mass2


    MODULE FUNCTION get_grav_mass1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_grav_mass1

    END FUNCTION get_grav_mass1


    MODULE FUNCTION get_grav_mass2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_grav_mass2

    END FUNCTION get_grav_mass2


    MODULE FUNCTION get_adm_mass( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_adm_mass

    END FUNCTION get_adm_mass


    MODULE FUNCTION get_angular_momentum( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_angular_momentum

    END FUNCTION get_angular_momentum


    MODULE FUNCTION get_radius1_x_comp( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_radius1_x_comp

    END FUNCTION get_radius1_x_comp


    MODULE FUNCTION get_radius1_y( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_radius1_y

    END FUNCTION get_radius1_y


    MODULE FUNCTION get_radius1_z( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_radius1_z

    END FUNCTION get_radius1_z


    MODULE FUNCTION get_radius1_x_opp( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_radius1_x_opp

    END FUNCTION get_radius1_x_opp


    MODULE FUNCTION get_center1_x( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_center1_x

    END FUNCTION get_center1_x


    MODULE FUNCTION get_barycenter1_x( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_barycenter1_x

    END FUNCTION get_barycenter1_x


    MODULE FUNCTION get_radius2_x_comp( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_radius2_x_comp

    END FUNCTION get_radius2_x_comp


    MODULE FUNCTION get_radius2_y( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_radius2_y

    END FUNCTION get_radius2_y


    MODULE FUNCTION get_radius2_z( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_radius2_z

    END FUNCTION get_radius2_z


    MODULE FUNCTION get_radius2_x_opp( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_radius2_x_opp

    END FUNCTION get_radius2_x_opp


    MODULE FUNCTION get_center2_x( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_center2_x

    END FUNCTION get_center2_x


    MODULE FUNCTION get_barycenter2_x( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_barycenter2_x

    END FUNCTION get_barycenter2_x


    MODULE FUNCTION get_ent_center1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_ent_center1

    END FUNCTION get_ent_center1


    MODULE FUNCTION get_nbar_center1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_nbar_center1

    END FUNCTION get_nbar_center1


    MODULE FUNCTION get_rho_center1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_rho_center1

    END FUNCTION get_rho_center1


    MODULE FUNCTION get_energy_density_center1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_energy_density_center1

    END FUNCTION get_energy_density_center1


    MODULE FUNCTION get_specific_energy_center1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_specific_energy_center1

    END FUNCTION get_specific_energy_center1


    MODULE FUNCTION get_pressure_center1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_pressure_center1

    END FUNCTION get_pressure_center1


    MODULE FUNCTION get_ent_center2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_ent_center2

    END FUNCTION get_ent_center2


    MODULE FUNCTION get_nbar_center2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_nbar_center2

    END FUNCTION get_nbar_center2


    MODULE FUNCTION get_rho_center2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_rho_center2

    END FUNCTION get_rho_center2


    MODULE FUNCTION get_energy_density_center2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_energy_density_center2

    END FUNCTION get_energy_density_center2


    MODULE FUNCTION get_specific_energy_center2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_specific_energy_center2

    END FUNCTION get_specific_energy_center2


    MODULE FUNCTION get_pressure_center2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_pressure_center2

    END FUNCTION get_pressure_center2


    MODULE FUNCTION get_eos1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      CHARACTER( LEN= : ), ALLOCATABLE:: get_eos1

    END FUNCTION get_eos1


    MODULE FUNCTION get_eos2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      CHARACTER( LEN= : ), ALLOCATABLE:: get_eos2

    END FUNCTION get_eos2


    MODULE FUNCTION get_eos1_id( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      INTEGER:: get_eos1_id

    END FUNCTION get_eos1_id


    MODULE FUNCTION get_eos2_id( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      INTEGER:: get_eos2_id

    END FUNCTION get_eos2_id


    MODULE FUNCTION get_npeos_1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      INTEGER:: get_npeos_1

    END FUNCTION get_npeos_1


    MODULE FUNCTION get_npeos_2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      INTEGER:: get_npeos_2

    END FUNCTION get_npeos_2


    MODULE FUNCTION get_gamma0_1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma0_1

    END FUNCTION get_gamma0_1


    MODULE FUNCTION get_gamma1_1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma1_1

    END FUNCTION get_gamma1_1


    MODULE FUNCTION get_gamma2_1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma2_1

    END FUNCTION get_gamma2_1


    MODULE FUNCTION get_gamma3_1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma3_1

    END FUNCTION get_gamma3_1


    MODULE FUNCTION get_kappa0_1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa0_1

    END FUNCTION get_kappa0_1


    MODULE FUNCTION get_kappa1_1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa1_1

    END FUNCTION get_kappa1_1


    MODULE FUNCTION get_kappa2_1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa2_1

    END FUNCTION get_kappa2_1


    MODULE FUNCTION get_kappa3_1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa3_1

    END FUNCTION get_kappa3_1


    MODULE FUNCTION get_logP1_1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logP1_1

    END FUNCTION get_logP1_1


    MODULE FUNCTION get_logRho0_1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logRho0_1

    END FUNCTION get_logRho0_1


    MODULE FUNCTION get_logRho1_1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logRho1_1

    END FUNCTION get_logRho1_1


    MODULE FUNCTION get_logRho2_1( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logRho2_1

    END FUNCTION get_logRho2_1


    MODULE FUNCTION get_gamma0_2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma0_2

    END FUNCTION get_gamma0_2


    MODULE FUNCTION get_gamma1_2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma1_2

    END FUNCTION get_gamma1_2


    MODULE FUNCTION get_gamma2_2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma2_2

    END FUNCTION get_gamma2_2


    MODULE FUNCTION get_gamma3_2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_gamma3_2

    END FUNCTION get_gamma3_2


    MODULE FUNCTION get_kappa0_2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa0_2

    END FUNCTION get_kappa0_2


    MODULE FUNCTION get_kappa1_2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa1_2

    END FUNCTION get_kappa1_2


    MODULE FUNCTION get_kappa2_2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa2_2

    END FUNCTION get_kappa2_2


    MODULE FUNCTION get_kappa3_2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_kappa3_2

    END FUNCTION get_kappa3_2


    MODULE FUNCTION get_logP1_2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logP1_2

    END FUNCTION get_logP1_2


    MODULE FUNCTION get_logRho0_2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logRho0_2

    END FUNCTION get_logRho0_2


    MODULE FUNCTION get_logRho1_2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logRho1_2

    END FUNCTION get_logRho1_2


    MODULE FUNCTION get_logRho2_2( THIS )

      !> [[bns]] object which this PROCEDURE is a member of
      CLASS(bnsbase), INTENT( IN ):: THIS
      ! Result
      DOUBLE PRECISION:: get_logRho2_2

    END FUNCTION get_logRho2_2

  END INTERFACE

END MODULE bns_base

