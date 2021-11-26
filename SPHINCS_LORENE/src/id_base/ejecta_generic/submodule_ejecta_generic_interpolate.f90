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

    INTEGER:: i

    DO i= 1, n, 1
      baryon_density(i) = THIS% read_mass_density( x(i), y(i), z(i) )
    ENDDO

    energy_density = 0.0D0
    specific_energy= 0.0D0
    pressure       = 0.0D0
    u_euler_x      = 0.0D0
    u_euler_y      = 0.0D0
    u_euler_z      = 0.0D0

    g_xx= 1.0D0
    g_yy= 1.0D0
    g_zz= 1.0D0
    g_xy= 0.0D0
    g_xz= 0.0D0
    g_yz= 0.0D0

    lapse= 1.0D0
    shift_x= 0.0D0
    shift_y= 0.0D0
    shift_z= 0.0D0

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

    USE tensor,    ONLY: jxx, jxy, jxz, jyy, jyz, jzz

    IMPLICIT NONE

    baryon_density= THIS% read_mass_density( x, y, z )

    g(jxx)= 1.0D0
    g(jyy)= 1.0D0
    g(jzz)= 1.0D0
    g(jxy)= 0.0D0
    g(jxz)= 0.0D0
    g(jyz)= 0.0D0

    gamma_euler= 1.0D0

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

    USE Hermite_refine, ONLY: find_indices
    !USE numerics,       ONLY: gf_pointer, pa_pointer


    IMPLICIT NONE


    INTEGER, PARAMETER:: nghost= 1
    INTEGER:: i, j, k, ierr, sgn_z

    DOUBLE PRECISION:: x0, y0, z0, x1, y1, z1, xd, yd, zd, &
                       c000, c001, c010, c100, c011, c110, c101, c111, &
                       c00, c01, c10, c11, c0, c1, zp

    sgn_z= SIGN(1.0D0,z)
    zp= ABS(z)

    CALL find_indices( x, y, zp, &
                       THIS% nx_grid, THIS% xL_grid, THIS% dx_grid, nghost, &
                       THIS% ny_grid, THIS% yL_grid, THIS% dy_grid, nghost, &
                       THIS% nz_grid, THIS% zL_grid, THIS% dz_grid, nghost, &
                       nghost, i, j, k, ierr )

   ! PRINT *, THIS% xL_grid, THIS% yL_grid, THIS% zL_grid
   ! PRINT *, x, y, z
   ! PRINT *, i, j, k
   ! PRINT *
   !
   ! PRINT *, THIS% grid(i,j,k,1), THIS% grid(i+1,j,k,1)

    IF( i >= THIS% nx_grid )THEN
      i= THIS% nx_grid - 1
    ENDIF
    IF( j >= THIS% ny_grid )THEN
      j= THIS% ny_grid - 1
    ENDIF
    IF( k >= THIS% nz_grid )THEN
      k= THIS% nz_grid - 1
    ENDIF

    !PRINT *, i, j, k
    !PRINT *, THIS% nx_grid, THIS% ny_grid, THIS% nz_grid

    IF( k <= 0 )THEN
      k= 1
    ENDIF
    IF( i <= 0 )THEN
      i= 1
    ENDIF
    IF( j <= 0 )THEN
      j= 1
    ENDIF

    x0= THIS% grid(i,j,k,1)
    x1= THIS% grid(i+1,j,k,1) !+ THIS% dx_grid
    y0= THIS% grid(i,j,k,2)
    y1= THIS% grid(i,j+1,k,2) !+ THIS% dy_grid
    z0= THIS% grid(i,j,k,3)
    z1= THIS% grid(i,j,k+1,3) !+ THIS% dz_grid

    xd= ( x - x0 )/( x1 - x0 )
    yd= ( y - y0 )/( y1 - y0 )
    zd= ( z - z0 )/( z1 - z0 )

  !  PRINT *, x, y, z
  !  PRINT *
  !
  !  PRINT *, xd, yd, zd
  !  PRINT *

    IF( k == 0 )THEN
      c001= 0.0D0 +1.0D0*THIS% baryon_mass_density(i,j,k)
      c101= 0.0D0 +1.0D0*THIS% baryon_mass_density(i+1,j,k)
      c011= 0.0D0 +1.0D0*THIS% baryon_mass_density(i,j+1,k)
      c111= 0.0D0 +1.0D0*THIS% baryon_mass_density(i+1,j+1,k)
    ENDIF
    IF( i >= THIS% nx_grid .OR. i == 0 )THEN
      c001= 0.0D0
      c101= 0.0D0
      c011= 0.0D0
      c111= 0.0D0
    ENDIF
    IF( j >= THIS% ny_grid .OR. j == 0 )THEN
      c001= 0.0D0
      c101= 0.0D0
      c011= 0.0D0
      c111= 0.0D0
    ENDIF
    IF( k >= THIS% nz_grid )THEN
      c001= 0.0D0
      c101= 0.0D0
      c011= 0.0D0
      c111= 0.0D0
    ENDIF

    c000= 0.0D0 +1.0D0*THIS% baryon_mass_density(i,j,k)
    c100= 0.0D0 +1.0D0*THIS% baryon_mass_density(i+1,j,k)
    c001= 0.0D0 +1.0D0*THIS% baryon_mass_density(i,j,k+1)
    c101= 0.0D0 +1.0D0*THIS% baryon_mass_density(i+1,j,k+1)
    c010= 0.0D0 +1.0D0*THIS% baryon_mass_density(i,j+1,k)
    c110= 0.0D0 +1.0D0*THIS% baryon_mass_density(i+1,j+1,k)
    c011= 0.0D0 +1.0D0*THIS% baryon_mass_density(i,j+1,k+1)
    c111= 0.0D0 +1.0D0*THIS% baryon_mass_density(i+1,j+1,k+1)

    c00= c000*( 1.0D0 - xd ) + c100*xd

    c01= c001*( 1.0D0 - xd ) + c101*xd

    c10= c010*( 1.0D0 - xd ) + c110*xd

    c11= c011*( 1.0D0 - xd ) + c111*xd

    c0= c00*( 1.0D0 - yd ) + c10*yd
    c1= c01*( 1.0D0 - yd ) + c11*yd

    res= c0*( 1.0D0 - zd ) + c1*zd

    !IF( res < 1.0D-13 ) res= 0.0D0
    IF( SQRT( ( x - THIS% centers(1,1) )**2.0D0 &
              + ( y - THIS% centers(1,2) )**2.0D0 &
              + ( zp - THIS% centers(1,3) )**2.0D0 ) > 475.0D0 ) res= 0.0D0

   ! PRINT *, c000, &
   !          c100, &
   !          c001, &
   !          c101, &
   !          c010, &
   !          c110, &
   !          c011, &
   !          c111, &
   !          res

  !  INTEGER, PARAMETER:: np= 1
  !  INTEGER, PARAMETER:: nvars= 1
  !  INTEGER, PARAMETER:: nghost= 0
  !
  !  DOUBLE PRECISION, DIMENSION(np):: xp, yp, zp
  !  DOUBLE PRECISION, DIMENSION(1), TARGET:: dens
  !
  !  TYPE(gf_pointer), DIMENSION(nvars):: density_grid
  !  TYPE(pa_pointer), DIMENSION(nvars):: density_point
  !
  !  xp(1)= x
  !  yp(1)= y
  !  zp(1)= z
  !
  !  density_point(1)% p => dens
  !  density_grid(1)%  p => THIS% baryon_mass_density
  !
  !  CALL H3_interpolate( &
  !        np,            &  ! The number of particle positions
  !        xp, yp, zp,    &  ! The x,y,z arrays of the interpolation points
  !        THIS% nx_grid, &  ! The dimensions of the grid function
  !        THIS% xL_grid, &  ! The lower coordinate values of the grid function
  !        THIS% dx_grid, &  ! The grid spacings
  !        nghost,        &  ! Number of ghost cells
  !        THIS% ny_grid, &
  !        THIS% yL_grid, &
  !        THIS% dy_grid, &
  !        nghost,        &
  !        THIS% nz_grid, &
  !        THIS% zL_grid, &
  !        THIS% dz_grid, &
  !        nghost,        &
  !        nvars,         &  ! Number of functions to be interpolated
  !        density_grid,  &  ! Array of pointers to grid functions of the data
  !        density_point  &  ! Array of pointers to the point array with
  !                          ! interpolated data
  !       )
  !
  !  res= dens(1)

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

    DOUBLE PRECISION, DIMENSION(3):: center

    center= THIS% return_center(1)

    IF( THIS% read_mass_density( x, y, z ) <= 2.0D-13 &
        .OR. &
        SQRT( ( x - center(1) )**2 + ( y - center(2) )**2 &
            + ( z - center(3) )**2  ) > 500.0D0 )THEN
      res= 1
    ELSE
      res= 0
    ENDIF

  END PROCEDURE is_hydro_negative


END SUBMODULE ejecta_generic_interpolate
