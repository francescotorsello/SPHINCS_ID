! File:         submodule_ejecta_generic_interpolate.f90
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

SUBMODULE (ejecta_generic) interpolate

  !****************************************************
  !
  !# Implementation of the methods of TYPE ejecta that
  !  interpolate the data on a grid, to a generic point
  !
  !  FT 19.11.2021
  !
  !****************************************************


  USE utility,  ONLY: zero, one


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE interpolate_id_full

    !**************************************************
    !
    !# Stores the |id| in non-[[ejecta]]-member arrays
    !  with the same shape as the [[ejecta]] member arrays
    !
    !  FT 19.11.2021
    !
    !**************************************************

    IMPLICIT NONE

  END PROCEDURE interpolate_id_full


  MODULE PROCEDURE interpolate_id_spacetime

    !*******************************************************
    !
    !# Stores the spacetime |id| in multi-dimensional arrays
    !  needed to compute the BSSN variables and constraints
    !
    !  FT 19.11.2021
    !
    !*******************************************************

    USE tensor, ONLY: jxx, jxy, jxz, jyy, jyz, jzz, jx, jy, jz, n_sym4x4

    IMPLICIT NONE

    INTEGER:: i, j, k

    PRINT *, "** N.B. : This SUBROUTINE is not yet implemented!"
    PRINT *, " * Stopping..."
    PRINT *
    STOP

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( nx, ny, nz, this, pos, &
    !$OMP                     lapse, shift, g, ek ) &
    !$OMP             PRIVATE( i, j, k )
    DO k= 1, nz, 1
      DO j= 1, ny, 1
        DO i= 1, nx, 1

          g( i, j, k, jyy )= one
          g( i, j, k, jyy )= one
          g( i, j, k, jzz )= one
          g( i, j, k, jxy )= zero
          g( i, j, k, jxz )= zero
          g( i, j, k, jyz )= zero

          !
          !- Set/unset the geodesic gauge
          !
          IF( this% get_one_lapse() )THEN
            lapse( i, j, k )= one
          ENDIF
          IF( this% get_zero_shift() )THEN
            shift( i, j, k, jx )= zero
            shift( i, j, k, jy )= zero
            shift( i, j, k, jz )= zero
          ENDIF

          !
          !-- Convert the extrinsic curvature from |lorene| units to
          !-- |sphincs| units
          !
          ek( i, j, k, jxx )= zero
          ek( i, j, k, jxy )= zero
          ek( i, j, k, jxz )= zero
          ek( i, j, k, jyy )= zero
          ek( i, j, k, jyz )= zero
          ek( i, j, k, jzz )= zero

        ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END PROCEDURE interpolate_id_spacetime


  MODULE PROCEDURE interpolate_id_hydro

    !*******************************************************
    !
    !# Stores the hydro |id| in the arrays needed to compute
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
    !# Stores the hydro |id| in the arrays needed to
    !  compute the |sph| |id|
    !
    !  FT 19.11.2020
    !
    !****************************************************

    USE constants,  ONLY: MSun, amu
    USE utility,    ONLY: two
    USE numerics,   ONLY: trilinear_interpolation

    IMPLICIT NONE

    LOGICAL, PARAMETER:: debug= .FALSE.
    !

    INTEGER:: i, j, k
    DOUBLE PRECISION:: zp, xtmp, ytmp, ztmp

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile
    LOGICAL:: exist

    DOUBLE PRECISION:: foo(n), foo_exact(n), &
                       foo_grid(this% nx_grid, this% ny_grid, this% nz_grid), &
                       grid_coords(this%nx_grid,this%ny_grid,this%nz_grid,3), &
                       coords(n,3)

    IF( debug )THEN

      DO i= 1, this% nx_grid, 1
        DO j= 1, this% ny_grid, 1
          DO k= 1, this% nz_grid, 1

            grid_coords(i,j,k,1)= DBLE(i) - DBLE(this% nx_grid)/two
            grid_coords(i,j,k,2)= DBLE(j) - DBLE(this% ny_grid)/two
            grid_coords(i,j,k,3)= DBLE(k)/two! - DBLE(this% nz_grid)/two
            foo_grid(i,j,k)= (grid_coords(i,j,k,3))**3.D0

          ENDDO
        ENDDO
      ENDDO

      foo= zero
      DO i= 1, n, 1

        !CALL RANDOM_NUMBER( xsgn )
        !CALL RANDOM_NUMBER( ysgn )
        !CALL RANDOM_NUMBER( zsgn )
        CALL RANDOM_NUMBER( xtmp )
        CALL RANDOM_NUMBER( ytmp )
        CALL RANDOM_NUMBER( ztmp )

        coords(i,1)= xtmp*DBLE(this% nx_grid - 2) &
                     - DBLE(this% nx_grid)/two + two
        coords(i,2)= ytmp*DBLE(this% ny_grid - 2) &
                     - DBLE(this% ny_grid)/two + two
        coords(i,3)= (- DBLE(this% nz_grid)/two + one)*(one-ztmp) &
                      + (DBLE(this% nz_grid)/two - one)*ztmp

        foo(i)= trilinear_interpolation( coords(i,1), coords(i,2), coords(i,3),&
                      this% nx_grid, this% ny_grid, this% nz_grid, &
                      grid_coords, foo_grid, &
                      equator_symmetry= .FALSE., parity= -one, debug= .FALSE. )
        foo_exact(i)= (coords(i,3))**3.D0

      ENDDO

    ENDIF

    ! The density has to be converted in units of the atomic mass unit
    ! TODO: CHECK THAT EVERYTHING ELSE IS CONSISTENT WITH this!!
    DO i= 1, n, 1

      zp= z(i)

      baryon_density(i) = this% read_mass_density( x(i), y(i), zp )*MSun/amu

      u_euler_x(i)      = trilinear_interpolation( x(i), y(i), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% grid, this% vel(:,:,:,1), &
                                equator_symmetry= .TRUE., parity= one, &
                                debug= .FALSE. )
      u_euler_y(i)      = trilinear_interpolation( x(i), y(i), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% grid, this% vel(:,:,:,2), &
                                equator_symmetry= .TRUE., parity= one, &
                                debug= .FALSE. )
      u_euler_z(i)      = trilinear_interpolation( x(i), y(i), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% grid, this% vel(:,:,:,3), &
                                equator_symmetry= .TRUE., parity= -one, &
                                debug= .FALSE. )

      specific_energy(i)= trilinear_interpolation( x(i), y(i), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% grid, this% specific_energy, &
                                equator_symmetry= .TRUE., parity= one, &
                                debug= .FALSE. )

      IF( baryon_density(i) == zero )THEN
        specific_energy(i)= zero
        u_euler_x(i)      = zero
        u_euler_y(i)      = zero
        u_euler_z(i)      = zero
      ENDIF

    ENDDO

    IF( debug )THEN

      finalnamefile= "dbg_interpolation2.dat"

      INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

      IF( exist )THEN
          OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
                FORM= "FORMATTED", &
                POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
                IOMSG= err_msg )
      ELSE
          OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "NEW", &
                FORM= "FORMATTED", &
                ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
      ENDIF
      IF( ios > 0 )THEN
        PRINT *, "...error when opening " // TRIM(finalnamefile), &
                 ". The error message is", err_msg
        STOP
      ENDIF

      DO i= 1, this% nx_grid - 1, 1
        DO j= 1, this% ny_grid - 1, 1
          DO k= 1, this% nz_grid - 1, 1

            WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
              this% grid( i, j, k, 1 ), &
              this% grid( i, j, k, 2 ), &
              this% grid( i, j, k, 3 ), &
              this% baryon_mass_density( i, j, k )*Msun/amu, &
              this% read_mass_density( &
                this% grid( i, j, k, 1 ) + this% dx_grid/two, &
                this% grid( i, j, k, 2 ), &
                this% grid( i, j, k, 3 ) ), &
              this% grid( i, j, k, 1 ) + this% dx_grid/two, &
              this% specific_energy( i, j, k )
          ENDDO
        ENDDO
      ENDDO

      CLOSE( UNIT= 2 )


      finalnamefile= "dbg_interpolation.dat"

      INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

      IF( exist )THEN
        OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
              FORM= "FORMATTED", &
              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
              IOMSG= err_msg )
      ELSE
        OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "NEW", &
              FORM= "FORMATTED", &
              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
      ENDIF
      IF( ios > 0 )THEN
        PRINT *, "...error when opening " // TRIM(finalnamefile), &
                 ". The error message is", err_msg
        STOP
      ENDIF

      DO i= 1, n, 1

        ! IF( coords(i,3) < 0 )THEN
        !   PRINT *, coords(i,3)
        !   STOP
        ! ENDIF

        WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
          i, x(i), y(i), z(i), &
          baryon_density(i), &
          u_euler_x(i), &
          u_euler_y(i), &
          u_euler_z(i), &
          specific_energy(i), &
          coords(i,1), coords(i,2), coords(i,3), &
          foo(i), foo_exact(i)

      ENDDO

      CLOSE( UNIT= 2 )

      !STOP

    ENDIF

    energy_density = zero
    pressure       = zero

    g_xx= one
    g_yy= one
    g_zz= one
    g_xy= zero
    g_xz= zero
    g_yz= zero

    lapse= one
    shift_x= zero
    shift_y= zero
    shift_z= zero

  END PROCEDURE interpolate_id_particles


  MODULE PROCEDURE interpolate_id_mass_b

    !****************************************************
    !
    !# Stores the hydro |id| in the arrays needed to
    !  compute the baryon mass, storing it to variables
    !  (not arrays as the others SUBROUTINES in
    !  the [[ejecta_generic_interpolate]] SUBMODULE).
    !
    !  FT 19.11.2021
    !
    !****************************************************

    USE tensor,    ONLY: jxx, jxy, jxz, jyy, jyz, jzz
    USE constants, ONLY: MSun, amu

    IMPLICIT NONE

    baryon_density= this% read_mass_density( x, y, z )

    g(jxx)= one
    g(jyy)= one
    g(jzz)= one
    g(jxy)= zero
    g(jxz)= zero
    g(jyz)= zero

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

    !USE timing,    ONLY: timer
    USE constants, ONLY: pi
    USE numerics,  ONLY: trilinear_interpolation
    USE utility,   ONLY: spherical_from_cartesian, two


    IMPLICIT NONE

    DOUBLE PRECISION:: zp, x_ell, y_ell, z_ell, theta, phi, r

    !TYPE(timer):: dbg_timer

    !dbg_timer= timer("dbg_timer")
    !CALL dbg_timer% start_timer()

    zp= z
    res= trilinear_interpolation( x, y, zp, &
                                  this% nx_grid, this% ny_grid, this% nz_grid, &
                                  this% grid, this% baryon_mass_density, &
                                  equator_symmetry= .TRUE., parity= one, &
                                  debug= .FALSE. )

    CALL spherical_from_cartesian( x, y, z, &
                  this% centers(1,1), this% centers(1,2), this% centers(1,3), &
                                   r, theta, phi )

    x_ell= this% centers(1,1) &
           + MAX(this% sizes(1,1),this% sizes(1,2))*COS(phi)*SIN(theta)

    y_ell= this% centers(1,2) &
           + MAX(this% sizes(1,3),this% sizes(1,4))*SIN(phi)*SIN(theta)

    z_ell= this% centers(1,3) &
           + MAX(this% sizes(1,5),this% sizes(1,6))*COS(theta)

    IF( r >= SQRT( ( x_ell - this% centers(1,1) )**two &
                 + ( y_ell - this% centers(1,2) )**two &
                 + ( z_ell - this% centers(1,3) )**two ) ) res= zero

    IF( res < zero ) res= zero

  !  IF(      x > this% centers(1,1) + this% sizes(1,2) &
  !      .OR. x < this% centers(1,1) - this% sizes(1,1) &
  !      .OR. y > this% centers(1,2) + this% sizes(1,4) &
  !      .OR. y < this% centers(1,2) - this% sizes(1,3) &
  !      .OR. zp > this% centers(1,3) + this% sizes(1,6) ) res= zero

  !   IF(      x > this% xR_grid &
  !       .OR. x < this% xL_grid &
  !       .OR. y > this% yR_grid &
  !       .OR. y < this% yL_grid &
  !       .OR. zp > this% zR_grid ) res= zero

    !CALL dbg_timer% stop_timer()
    !CALL dbg_timer% print_timer( 2 )
    !STOP

  END PROCEDURE interpolate_mass_density


  MODULE PROCEDURE interpolate_spatial_metric

    !***********************************************
    !
    !# Returns the spatial metric.
    !
    !  FT 19.11.2021
    !
    !***********************************************

    IMPLICIT NONE

  END PROCEDURE interpolate_spatial_metric


  MODULE PROCEDURE is_hydro_positive

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

    center= this% return_center(1)
    sizes = this% return_spatial_extent(1)

    IF( this% read_mass_density( x, y, z ) <= zero &
        .OR. &
        !SQRT( ( x - center(1) )**2 + ( y - center(2) )**2 &
        !    + ( z - center(3) )**2  ) > 50zero
             x > center(1) + sizes(1) &
        .OR. x < center(1) - sizes(2) &
        .OR. y > center(2) + sizes(3) &
        .OR. y < center(2) - sizes(4) &
        .OR. ABS(z) > center(3) + sizes(5) &
    )THEN
      res= .FALSE.
    ELSE
      res= .TRUE.
    ENDIF

  END PROCEDURE is_hydro_positive


END SUBMODULE interpolate
