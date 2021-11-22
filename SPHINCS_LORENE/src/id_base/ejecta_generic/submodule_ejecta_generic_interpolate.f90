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

    USE Hermite_refine, ONLY: H3_interpolate
    USE numerics,       ONLY: gf_pointer, pa_pointer


    IMPLICIT NONE


    INTEGER, PARAMETER:: np= 1
    INTEGER, PARAMETER:: nvars= 1
    INTEGER, PARAMETER:: nghost= 0

    DOUBLE PRECISION, DIMENSION(np):: xp, yp, zp
    DOUBLE PRECISION, DIMENSION(1), TARGET:: dens

    TYPE(gf_pointer), DIMENSION(nvars):: density_grid
    TYPE(pa_pointer), DIMENSION(nvars):: density_point

    xp(1)= x
    yp(1)= y
    zp(1)= z

    density_point(1)% p => dens
    density_grid(1)%  p => THIS% baryon_mass_density

    CALL H3_interpolate( &
          np,            &  ! The number of particle positions
          xp, yp, zp,    &  ! The x,y,z arrays of the interpolation points
          THIS% nx_grid, &  ! The dimensions of the grid function
          THIS% xL_grid, &  ! The lower coordinate values of the grid function
          THIS% dx_grid, &  ! The grid spacings
          nghost,        &  ! Number of ghost cells
          THIS% ny_grid, &
          THIS% yL_grid, &
          THIS% dy_grid, &
          nghost,        &
          THIS% nz_grid, &
          THIS% zL_grid, &
          THIS% dz_grid, &
          nghost,        &
          nvars,         &  ! Number of functions to be interpolated
          density_grid,  &  ! Array of pointers to grid functions of the data
          density_point  &  ! Array of pointers to the point array with
                            ! interpolated data
         )

    res= dens(1)

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
