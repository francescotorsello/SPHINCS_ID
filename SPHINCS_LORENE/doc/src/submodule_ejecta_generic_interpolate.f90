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


  USE constants, ONLY: one


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

    USE constants, ONLY: MSun, amu
    USE numerics,  ONLY: trilinear_interpolation

    IMPLICIT NONE

    INTEGER:: i, j, k, i_eps, i_vel
    DOUBLE PRECISION:: zp, min_eps, min_vel, xtmp, ytmp, ztmp, xsgn, ysgn, zsgn

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile
    LOGICAL:: exist

    DOUBLE PRECISION:: foo(n), foo_exact(n), &
                       foo_grid(THIS% nx_grid, THIS% ny_grid, THIS% nz_grid), &
                       grid_coords(THIS% nx_grid, THIS% ny_grid, THIS% nz_grid, 3), &
                       coords(n,3)

!    DO i= 1, THIS% nx_grid, 1
!      DO j= 1, THIS% ny_grid, 1
!        DO k= 1, THIS% nz_grid, 1
!
!          grid_coords(i,j,k,1)= DBLE(i) - DBLE(THIS% nx_grid)/2.D0
!          grid_coords(i,j,k,2)= DBLE(j) - DBLE(THIS% ny_grid)/2.D0
!          grid_coords(i,j,k,3)= DBLE(k)/2.D0! - DBLE(THIS% nz_grid)/2.D0
!          foo_grid(i,j,k)= (grid_coords(i,j,k,3))**3.D0
!
!        ENDDO
!      ENDDO
!    ENDDO
!
!    foo= 0.D0
!    DO i= 1, n, 1
!
!      !CALL RANDOM_NUMBER( xsgn )
!      !CALL RANDOM_NUMBER( ysgn )
!      !CALL RANDOM_NUMBER( zsgn )
!      CALL RANDOM_NUMBER( xtmp )
!      CALL RANDOM_NUMBER( ytmp )
!      CALL RANDOM_NUMBER( ztmp )
!
!      coords(i,1)= xtmp*DBLE(THIS% nx_grid - 2) - DBLE(THIS% nx_grid)/2.D0 + 2.D0
!      coords(i,2)= ytmp*DBLE(THIS% ny_grid - 2) - DBLE(THIS% ny_grid)/2.D0 + 2.D0
!      coords(i,3)= (- DBLE(THIS% nz_grid)/2.D0 + one)*(one-ztmp) + (DBLE(THIS% nz_grid)/2.D0 - one)*ztmp
!
!      foo(i)= trilinear_interpolation( coords(i,1), coords(i,2), coords(i,3), &
!                    THIS% nx_grid, THIS% ny_grid, THIS% nz_grid, &
!                    grid_coords, foo_grid, &
!                    equator_symmetry= .FALSE., parity= -one, debug= .FALSE. )
!      foo_exact(i)= (coords(i,3))**3.D0
!    ENDDO

 !   min_eps= HUGE(one)
 !   min_vel= HUGE(one)
    ! The density has to be converted in units of the atomic mass unit
    ! TODO: CHECK THAT EVERYTHING ELSE IS CONSISTENT WITH THIS!!
    DO i= 1, n, 1

      zp= z(i)

      baryon_density(i) = THIS% read_mass_density( x(i), y(i), zp )*MSun/amu

      u_euler_x(i)      = trilinear_interpolation( x(i), y(i), zp, &
                                THIS% nx_grid, THIS% ny_grid, THIS% nz_grid, &
                                THIS% grid, THIS% vel(:,:,:,1), &
                                equator_symmetry= .TRUE., parity= one, &
                                debug= .FALSE. )
      u_euler_y(i)      = trilinear_interpolation( x(i), y(i), zp, &
                                THIS% nx_grid, THIS% ny_grid, THIS% nz_grid, &
                                THIS% grid, THIS% vel(:,:,:,2), &
                                equator_symmetry= .TRUE., parity= one, &
                                debug= .FALSE. )
      u_euler_z(i)      = trilinear_interpolation( x(i), y(i), zp, &
                                THIS% nx_grid, THIS% ny_grid, THIS% nz_grid, &
                                THIS% grid, THIS% vel(:,:,:,3), &
                                equator_symmetry= .TRUE., parity= -one, &
                                debug= .FALSE. )

    !  IF( u_euler_x(i) == 0 .AND. u_euler_y(i) == 0 &
    !      .AND. u_euler_z(i) == 0 )THEN
    !    PRINT *, u_euler_x(i), u_euler_y(i), u_euler_z(i)
    !    PRINT *, x(i), y(i), zp
    !    STOP
    !  ENDIF

      specific_energy(i)= trilinear_interpolation( x(i), y(i), zp, &
                                THIS% nx_grid, THIS% ny_grid, THIS% nz_grid, &
                                THIS% grid, THIS% specific_energy, &
                                equator_symmetry= .TRUE., parity= one, &
                                debug= .FALSE. )

      IF( baryon_density(i) == 0.D0 )THEN
        specific_energy(i)= 0.D0
        u_euler_x(i)      = 0.D0
        u_euler_y(i)      = 0.D0
        u_euler_z(i)      = 0.D0
      ENDIF

    !  IF( specific_energy(i) <= 0 )THEN
    !    PRINT *, specific_energy(i)
    !    PRINT *, x(i), y(i), zp
    !    STOP
    !  ENDIF

  !    IF( SQRT( u_euler_x(i)**2.D0 + u_euler_y(i)**2.D0 &
  !            + u_euler_z(i)**2.D0 ) < min_vel &
  !        .AND. baryon_density(i) > 0.D0 )THEN
  !      min_vel= SQRT( u_euler_x(i)**2.D0 + u_euler_y(i)**2.D0 &
  !                   + u_euler_z(i)**2.D0 )
  !      i_vel= i
  !    ENDIF
  !    IF( specific_energy(i) < min_eps .AND. baryon_density(i) > 0.D0 )THEN
  !      min_eps= specific_energy(i)
  !      i_eps= i
  !    ENDIF

    ENDDO

  !  PRINT *
  !  PRINT *, MINVAL( specific_energy, DIM= 1, MASK= baryon_density > 0.D0 ), &
  !           min_eps
  !  PRINT *, x(i_eps), y(i_eps), z(i_eps)
  !
  !  PRINT *, MINVAL( SQRT( (u_euler_x)**2.D0 &
  !                 + (u_euler_y)**2.D0 &
  !                 + (u_euler_z)**2.D0 ), &
  !                   DIM= 1, MASK= baryon_density > 0.D0 ), &
  !           min_vel
  !  PRINT *, x(i_vel), y(i_vel), z(i_vel)
  !
  !  PRINT *

!    finalnamefile= "dbg_interpolation2.dat"
!
!    INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )
!
!    IF( exist )THEN
!        OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
!              FORM= "FORMATTED", &
!              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
!              IOMSG= err_msg )
!    ELSE
!        OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "NEW", &
!              FORM= "FORMATTED", &
!              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
!    ENDIF
!    IF( ios > 0 )THEN
!      PRINT *, "...error when opening " // TRIM(finalnamefile), &
!               ". The error message is", err_msg
!      STOP
!    ENDIF
!
!    DO i= 1, THIS% nx_grid - 1, 1
!      DO j= 1, THIS% ny_grid - 1, 1
!        DO k= 1, THIS% nz_grid - 1, 1
!
!          WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
!            THIS% grid( i, j, k, 1 ), &
!            THIS% grid( i, j, k, 2 ), &
!            THIS% grid( i, j, k, 3 ), &
!            THIS% baryon_mass_density( i, j, k )*Msun/amu, &
!            THIS% read_mass_density( &
!              THIS% grid( i, j, k, 1 ) + THIS% dx_grid/2.D0, &
!              THIS% grid( i, j, k, 2 ), &
!              THIS% grid( i, j, k, 3 ) ), &
!            THIS% grid( i, j, k, 1 ) + THIS% dx_grid/2.D0, &
!            THIS% specific_energy( i, j, k )
!        ENDDO
!      ENDDO
!    ENDDO
!
!    CLOSE( UNIT= 2 )
!
!
!    finalnamefile= "dbg_interpolation.dat"
!
!    INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )
!
!    IF( exist )THEN
!      OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
!            FORM= "FORMATTED", &
!            POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
!            IOMSG= err_msg )
!    ELSE
!      OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "NEW", &
!            FORM= "FORMATTED", &
!            ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
!    ENDIF
!    IF( ios > 0 )THEN
!      PRINT *, "...error when opening " // TRIM(finalnamefile), &
!               ". The error message is", err_msg
!      STOP
!    ENDIF
!
!    DO i= 1, n, 1
!
!      ! IF( coords(i,3) < 0 )THEN
!      !   PRINT *, coords(i,3)
!      !   STOP
!      ! ENDIF
!
!      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
!        i, x(i), y(i), z(i), &
!        baryon_density(i), &
!        u_euler_x(i), &
!        u_euler_y(i), &
!        u_euler_z(i), &
!        specific_energy(i), &
!        coords(i,1), coords(i,2), coords(i,3), &
!        foo(i), foo_exact(i)
!
!    ENDDO
!
!    CLOSE( UNIT= 2 )

    !STOP

    energy_density = 0.D0
    pressure       = 0.D0

    g_xx= one
    g_yy= one
    g_zz= one
    g_xy= 0.D0
    g_xz= 0.D0
    g_yz= 0.D0

    lapse= one
    shift_x= 0.D0
    shift_y= 0.D0
    shift_z= 0.D0

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
    USE constants, ONLY: MSun, amu

    IMPLICIT NONE

    baryon_density= THIS% read_mass_density( x, y, z )

    g(jxx)= one
    g(jyy)= one
    g(jzz)= one
    g(jxy)= 0.D0
    g(jxz)= 0.D0
    g(jyz)= 0.D0

    gamma_euler= one

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
    !# Returns the mass density at the point
    !  given as argument, in units of
    !  \(M_\odot/L_\odot^3\).
    !
    !  FT 19.11.2021
    !
    !***********************************************

    !USE Hermite_refine, ONLY: find_indices
    USE constants, ONLY: pi
    USE numerics,  ONLY: trilinear_interpolation


    IMPLICIT NONE


    INTEGER, PARAMETER:: nghost= 1
    INTEGER:: i, j, k, ierr, sgn_z

    DOUBLE PRECISION:: x0, y0, z0, x1, y1, z1, xd, yd, zd, &
                       c000, c001, c010, c100, c011, c110, c101, c111, &
                       c00, c01, c10, c11, c0, c1, zp, x_ell, y_ell, z_ell, &
                       theta, phi

    zp= z
    res= trilinear_interpolation( x, y, zp, &
                                  THIS% nx_grid, THIS% ny_grid, THIS% nz_grid, &
                                  THIS% grid, THIS% baryon_mass_density, &
                                  equator_symmetry= .TRUE., parity= one, &
                                  debug= .FALSE. )

    IF( x > 0.D0 )THEN

      phi= ATAN( ( y - THIS% centers(1,2) )/( x - THIS% centers(1,1) ) )

    ELSEIF( x < 0.D0 )THEN

      phi= ATAN( ( y - THIS% centers(1,2) )/( x - THIS% centers(1,1) ) ) + pi

    ELSE

      phi= pi/2.D0

    ENDIF

    theta= ACOS( ( z - THIS% centers(1,3) ) &
          /SQRT( ( x - THIS% centers(1,1) )**2.D0 &
               + ( y - THIS% centers(1,2) )**2.D0 &
               + ( z - THIS% centers(1,3) )**2.D0 ) )

    x_ell= THIS% centers(1,1) &
           + MAX(THIS% sizes(1,1),THIS% sizes(1,2))*COS(phi)*SIN(theta)

    y_ell= THIS% centers(1,2) &
           + MAX(THIS% sizes(1,3),THIS% sizes(1,4))*SIN(phi)*SIN(theta)

    z_ell= THIS% centers(1,3) &
           + MAX(THIS% sizes(1,5),THIS% sizes(1,6))*COS(theta)

    IF( SQRT( ( x - THIS% centers(1,1) )**2.D0 &
            + ( y - THIS% centers(1,2) )**2.D0 &
            + ( z - THIS% centers(1,3) )**2.D0 ) >= &
        SQRT( ( x_ell - THIS% centers(1,1) )**2.D0 &
            + ( y_ell - THIS% centers(1,2) )**2.D0 &
            + ( z_ell - THIS% centers(1,3) )**2.D0 ) ) res= 0.D0

    IF( res < 0.D0 ) res= 0.D0

  !  IF(      x > THIS% centers(1,1) + THIS% sizes(1,2) &
  !      .OR. x < THIS% centers(1,1) - THIS% sizes(1,1) &
  !      .OR. y > THIS% centers(1,2) + THIS% sizes(1,4) &
  !      .OR. y < THIS% centers(1,2) - THIS% sizes(1,3) &
  !      .OR. zp > THIS% centers(1,3) + THIS% sizes(1,6) ) res= 0.D0

  !   IF(      x > THIS% xR_grid &
  !       .OR. x < THIS% xL_grid &
  !       .OR. y > THIS% yR_grid &
  !       .OR. y < THIS% yL_grid &
  !       .OR. zp > THIS% zR_grid ) res= 0.D0


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
    !  at the specified point; return 0 otherwise
    !
    !  FT 19.11.2021
    !
    !************************************************

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(3):: center
    DOUBLE PRECISION, DIMENSION(6):: sizes

    center= THIS% return_center(1)
    sizes = THIS% return_spatial_extent(1)

    IF( THIS% read_mass_density( x, y, z ) <= 0.D0 &
        .OR. &
        !SQRT( ( x - center(1) )**2 + ( y - center(2) )**2 &
        !    + ( z - center(3) )**2  ) > 500.D0
             x > center(1) + sizes(1) &
        .OR. x < center(1) - sizes(2) &
        .OR. y > center(2) + sizes(3) &
        .OR. y < center(2) - sizes(4) &
        .OR. ABS(z) > center(3) + sizes(5) &
    )THEN
      res= 1
    ELSE
      res= 0
    ENDIF

  END PROCEDURE is_hydro_negative


END SUBMODULE ejecta_generic_interpolate
