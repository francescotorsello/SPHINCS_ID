! File:         module_id_base.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

MODULE id_base

  !*******************************************************
  !
  !#
  !
  !   FT 24.09.2021
  !
  !*******************************************************


  IMPLICIT NONE


  !**********************************************************
  !                                                         *
  !              Definition of TYPE idbase                  *
  !                                                         *
  !   This abstract class represents a generic ID for       *
  !   SPHINCS_BSSN (binary neutron star, rotating star...). *
  !                                                         *
  !**********************************************************

  TYPE, ABSTRACT:: idbase
  !# Represents a generic ID for SPHINCS_BSSN (binary neutron star, rotating
  !  star...)


    CONTAINS


    PROCEDURE:: integrate_baryon_mass_density
    !# Integrates the LORENE baryon mass density and computes the
    !  radial mass profile

  END TYPE idbase


  INTERFACE

    MODULE SUBROUTINE integrate_baryon_mass_density( THIS, center, radius, &
                                                     central_density, &
                                                     dr, dth, dphi, &
                                                     mass, mass_profile, &
                                                     mass_profile_idx )
    !# Integrates the baryon mass density to compute the radial mass
    !  profile of a single star. TODO: Improve integration algorithm.

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

    END SUBROUTINE integrate_baryon_mass_density

  END INTERFACE


END MODULE id_base
