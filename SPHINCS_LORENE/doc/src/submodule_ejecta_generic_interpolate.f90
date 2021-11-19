! File:         submodule_ejecta_generic_interpolate.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (ejecta_generic) ejecta_generic_interpolate

  !****************************************************
  !
  !# Implementation of the methods of TYPE ejecta that
  !  interpolate the data on a grid, to a generic point
  !
  !  FT 19.11.2021
  !
  !****************************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE interpolate_id_full

    !**************************************************
    !
    !# Stores the ID in non-[[diffstarlorene]]-member arrays
    !  with the same shape as the [[diffstarlorene]] member arrays
    !
    !  FT 19.11.2021
    !
    !**************************************************

    IMPLICIT NONE

  END PROCEDURE interpolate_id_full


  MODULE PROCEDURE interpolate_id_spacetime

    !*******************************************************
    !
    !# Stores the spacetime ID in multi-dimensional arrays
    !  needed to compute the BSSN variables and constraints
    !
    !  FT 19.11.2021
    !
    !*******************************************************

    IMPLICIT NONE

  END PROCEDURE interpolate_id_spacetime


  MODULE PROCEDURE interpolate_id_hydro

    !*******************************************************
    !
    !# Stores the hydro ID in the arrays needed to compute
    !  the constraints on the refined mesh
    !
    !  FT 19.11.2021
    !
    !*******************************************************

    IMPLICIT NONE

  END PROCEDURE interpolate_id_hydro


  MODULE PROCEDURE interpolate_id_particles

    !****************************************************
    !
    !# Stores the hydro ID in the arrays needed to
    !  compute the SPH ID
    !
    !  FT 19.11.2020
    !
    !****************************************************

    IMPLICIT NONE

  END PROCEDURE interpolate_id_particles


  MODULE PROCEDURE interpolate_id_mass_b

    !****************************************************
    !
    !# Stores the hydro ID in the arrays needed to
    !  compute the baryon mass, storing it to variables
    !  (not arrays as the others SUBROUTINES in
    !  the [[bns_interpolate]] SUBMODULE).
    !
    !  FT 19.11.2021
    !
    !****************************************************

    IMPLICIT NONE

  END PROCEDURE interpolate_id_mass_b


  MODULE PROCEDURE interpolate_id_k

    !****************************************************
    !
    !# Stores the components of the extrinsic curvature
    !  in arrays
    !
    !  @warning DEPRECATED
    !
    !  FT 19.11.2021
    !
    !****************************************************

    IMPLICIT NONE

  END PROCEDURE interpolate_id_k


  !-----------------!
  !--  FUNCTIONS  --!
  !-----------------!


  MODULE PROCEDURE interpolate_mass_density

    !***********************************************
    !
    !# Returns the |lorene| mass density at the point
    !  given as argument, in units of
    !  \(M_\odot/L_\odot^3\).
    !
    !  FT 19.11.2021
    !
    !***********************************************

    USE Hermite_refine, ONLY: H5_interpolate

    IMPLICIT NONE

  END PROCEDURE interpolate_mass_density


  MODULE PROCEDURE interpolate_spatial_metric

    !***********************************************
    !
    !# Returns the |lorene| conformal factor to the
    !  4th power, equal to the diagonal components
    !  of the conformally flat spatial ADM metric.
    !
    !  FT 19.11.2021
    !
    !***********************************************

    IMPLICIT NONE

  END PROCEDURE interpolate_spatial_metric


  MODULE PROCEDURE is_hydro_negative

    !************************************************
    !
    !# Return 1 if the energy density is nonpositive
    !  or if the specific energy is nonpositive,
    !  or if the pressure is nonpositive
    !  at the specified point
    !
    !  FT 19.11.2021
    !
    !************************************************

    IMPLICIT NONE

  END PROCEDURE is_hydro_negative


END SUBMODULE ejecta_generic_interpolate
