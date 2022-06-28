! File:         submodule_bns_fuka_interpolate.f90
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

SUBMODULE (bns_fuka) interpolate

  !****************************************************
  !
  !# Implementation of the methods of TYPE bnsfuka that
  !  interpolate the data on a lattice, at a particle
  !  position
  !
  !  FT 28.06.2022
  !
  !****************************************************


  USE utility,  ONLY: zero, one


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE interpolate_fuka_id_particles

    !****************************************************
    !
    !# Stores the hydro |id| in the arrays needed to
    !  compute the |sph| |id|
    !
    !  FT 28.06.2022
    !
    !****************************************************

    USE constants,  ONLY: MSun, amu
    USE utility,    ONLY: two
    USE numerics,   ONLY: trilinear_interpolation

    IMPLICIT NONE

    LOGICAL, PARAMETER:: debug= .FALSE.
    !

    INTEGER:: a, i_star, star!, j, k
    DOUBLE PRECISION:: zp!, xtmp, ytmp, ztmp

    !CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile
    !LOGICAL:: exist

    !DOUBLE PRECISION:: foo(n), foo_exact(n), &
    !                   foo_grid(this% nx_grid, this% ny_grid, this% nz_grid), &
    !                   grid_coords(this%nx_grid,this%ny_grid,this%nz_grid,3), &
    !                   coords(n,3)

    DOUBLE PRECISION, DIMENSION(3):: center
    DOUBLE PRECISION, DIMENSION(6):: sizes

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( n, this, x, y, z, lapse, &
    !$OMP                     shift_x, shift_y, shift_z, &
    !$OMP                     g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, &
    !$OMP                     baryon_density, energy_density, &
    !$OMP                     specific_energy, pressure, &
    !$OMP                     u_euler_x, u_euler_y, u_euler_z ) &
    !$OMP             PRIVATE( a, i_star, star, sizes, center, zp )
    DO a= 1, n, 1

      loop_over_stars: DO i_star= 1, 2, 1

        sizes = this% return_spatial_extent(i_star)
        center= this% return_center(i_star)

        IF( center(1) - sizes(1) < x(a) .AND. x(a) < center(1) + sizes(2) ) &
          star= i_star

      ENDDO loop_over_stars

      zp= z(a)

      baryon_density(a) = trilinear_interpolation( x(a), y(a), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% id_fields(:,:,:,id$x:id$z,star), &
                                this% id_fields(:,:,:,id$massdensity,star), &
                                equator_symmetry= .FALSE., debug= .FALSE. ) &
                                *MSun/amu

      specific_energy(a)= trilinear_interpolation( x(a), y(a), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% id_fields(:,:,:,id$x:id$z,star), &
                                this% id_fields(:,:,:,id$specificenergy,star), &
                                equator_symmetry= .FALSE., debug= .FALSE. )

      pressure(a)       = trilinear_interpolation( x(a), y(a), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% id_fields(:,:,:,id$x:id$z,star), &
                                this% id_fields(:,:,:,id$pressure,star), &
                                equator_symmetry= .FALSE., debug= .FALSE. ) &
                                *MSun/amu

      u_euler_x(a)      = trilinear_interpolation( x(a), y(a), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% id_fields(:,:,:,id$x:id$z,star), &
                                this% id_fields(:,:,:,id$eulvelx,star), &
                                equator_symmetry= .FALSE., debug= .FALSE. )
      u_euler_y(a)      = trilinear_interpolation( x(a), y(a), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% id_fields(:,:,:,id$x:id$z,star), &
                                this% id_fields(:,:,:,id$eulvely,star), &
                                equator_symmetry= .FALSE., debug= .FALSE. )
      u_euler_z(a)      = trilinear_interpolation( x(a), y(a), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% id_fields(:,:,:,id$x:id$z,star), &
                                this% id_fields(:,:,:,id$eulvelz,star), &
                                equator_symmetry= .FALSE., debug= .FALSE. )

      IF( baryon_density(a) == zero )THEN
        specific_energy(a)= zero
        u_euler_x(a)      = zero
        u_euler_y(a)      = zero
        u_euler_z(a)      = zero
      ENDIF

      energy_density(a) = baryon_density(a)*(one + specific_energy(a))

      g_xx(a)           = trilinear_interpolation( x(a), y(a), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% id_fields(:,:,:,id$x:id$z,star), &
                                this% id_fields(:,:,:,id$gxx,star), &
                                equator_symmetry= .FALSE., debug= .FALSE. )

      g_yy(a)= g_xx(a)
      g_zz(a)= g_xx(a)
      g_xy(a)= zero
      g_xz(a)= zero
      g_yz(a)= zero

      lapse(a)          = trilinear_interpolation( x(a), y(a), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% id_fields(:,:,:,id$x:id$z,star), &
                                this% id_fields(:,:,:,id$lapse,star), &
                                equator_symmetry= .FALSE., debug= .FALSE. )
      shift_x(a)        = trilinear_interpolation( x(a), y(a), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% id_fields(:,:,:,id$x:id$z,star), &
                                this% id_fields(:,:,:,id$shiftx,star), &
                                equator_symmetry= .FALSE., debug= .FALSE. )
      shift_y(a)        = trilinear_interpolation( x(a), y(a), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% id_fields(:,:,:,id$x:id$z,star), &
                                this% id_fields(:,:,:,id$shifty,star), &
                                equator_symmetry= .FALSE., debug= .FALSE. )
      shift_z(a)        = trilinear_interpolation( x(a), y(a), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% id_fields(:,:,:,id$x:id$z,star), &
                                this% id_fields(:,:,:,id$shiftz,star), &
                                equator_symmetry= .FALSE., debug= .FALSE. )

    ENDDO
    !$OMP END PARALLEL DO

  END PROCEDURE interpolate_fuka_id_particles


  MODULE PROCEDURE interpolate_fuka_id_mass_b

    !****************************************************
    !
    !# Stores the hydro |id| in the arrays needed to
    !  compute the baryon mass, storing it to variables
    !  (not arrays as the others SUBROUTINES in
    !  the [[ejecta_generic_interpolate]] SUBMODULE).
    !
    !  FT 28.06.2022
    !
    !****************************************************

    USE tensor,    ONLY: jxx, jxy, jxz, jyy, jyz, jzz
    USE constants, ONLY: MSun, amu
    USE numerics,   ONLY: trilinear_interpolation

    IMPLICIT NONE

    INTEGER:: i_star, star

    DOUBLE PRECISION:: zp, veuler_x, veuler_y, veuler_z

    DOUBLE PRECISION, DIMENSION(6):: sizes
    DOUBLE PRECISION, DIMENSION(3):: center

    loop_over_stars: DO i_star= 1, 2, 1

      sizes = this% return_spatial_extent(i_star)
      center= this% return_center(i_star)

      IF( (center(1) - sizes(1) < x) .AND. (x < center(1) + sizes(2)) ) &
        star= i_star

    ENDDO loop_over_stars

    zp= z

    baryon_density= trilinear_interpolation( x, y, zp, &
                          this% nx_grid, this% ny_grid, this% nz_grid, &
                          this% id_fields(:,:,:,id$x:id$z,star), &
                          this% id_fields(:,:,:,id$massdensity,star), &
                          equator_symmetry= .FALSE., debug= .FALSE. )

    g(jxx)= trilinear_interpolation( x, y, zp, &
                  this% nx_grid, this% ny_grid, this% nz_grid, &
                  this% id_fields(:,:,:,id$x:id$z,star), &
                  this% id_fields(:,:,:,id$gxx,star), &
                  equator_symmetry= .FALSE., debug= .FALSE. )
    g(jyy)= g(jxx)
    g(jzz)= g(jxx)
    g(jxy)= zero
    g(jxz)= zero
    g(jyz)= zero

    veuler_x= trilinear_interpolation( x, y, zp, &
                     this% nx_grid, this% ny_grid, this% nz_grid, &
                     this% id_fields(:,:,:,id$x:id$z,star), &
                     this% id_fields(:,:,:,id$eulvelx,star), &
                     equator_symmetry= .FALSE., debug= .FALSE. )
    veuler_y= trilinear_interpolation( x, y, zp, &
                     this% nx_grid, this% ny_grid, this% nz_grid, &
                     this% id_fields(:,:,:,id$x:id$z,star), &
                     this% id_fields(:,:,:,id$eulvely,star), &
                     equator_symmetry= .FALSE., debug= .FALSE. )
    veuler_z= trilinear_interpolation( x, y, zp, &
                     this% nx_grid, this% ny_grid, this% nz_grid, &
                     this% id_fields(:,:,:,id$x:id$z,star), &
                     this% id_fields(:,:,:,id$eulvelz,star), &
                     equator_symmetry= .FALSE., debug= .FALSE. )

    ! See eq.(7.3.13) in Alcubierre, "Introduction to 3+1 Numerical Relativity"
    ! The following formula assumes a conformally flat metric in Cartesian
    ! coordinates
    gamma_euler= one/SQRT( one - g(jxx)*( veuler_x*veuler_x &
                                        + veuler_y*veuler_y &
                                        + veuler_z*veuler_z ) );

  END PROCEDURE interpolate_fuka_id_mass_b


  !-----------------!
  !--  FUNCTIONS  --!
  !-----------------!


  MODULE PROCEDURE interpolate_fuka_mass_density

    !***********************************************
    !
    !# Returns the mass density at the point
    !  given as argument, in units of
    !  \(M_\odot/L_\odot^3\).
    !
    !  FT 28.06.2022
    !
    !***********************************************

    USE constants, ONLY: pi
    USE numerics,  ONLY: trilinear_interpolation
    USE utility,   ONLY: spherical_from_cartesian, two


    IMPLICIT NONE

    INTEGER:: i_star, star

    DOUBLE PRECISION:: zp!, x_ell, y_ell, z_ell, theta, phi, r

    DOUBLE PRECISION, DIMENSION(6):: sizes
    DOUBLE PRECISION, DIMENSION(3):: center

    loop_over_stars: DO i_star= 1, 2, 1

      sizes = this% return_spatial_extent(i_star)
      center= this% return_center(i_star)

      !PRINT *, center(1) - sizes(1)
      !PRINT *, x
      !PRINT *, center(1) + sizes(2)

      IF( (center(1) - sizes(1) < x) .AND. (x < center(1) + sizes(2)) )THEN

        star= i_star
        !PRINT *, star

      ELSE

        star= -1

      ENDIF

    ENDDO loop_over_stars

    IF( star == -1 )THEN
      res= zero
      RETURN
    ENDIF

    !PRINT *, star
    !STOP

    zp= z
    res= trilinear_interpolation( x, y, zp, &
                                  this% nx_grid, this% ny_grid, this% nz_grid, &
                                  this% id_fields(:,:,:,id$x:id$z,star), &
                                  this% id_fields(:,:,:,id$massdensity,star), &
                                  equator_symmetry= .TRUE., parity= one, &
                                  debug= .FALSE. )

    !PRINT *, this% id_fields(:,:,:,id$massdensity,star)
    !STOP

  !  CALL spherical_from_cartesian( x, y, z, &
  !                this% centers(1,1), this% centers(1,2), this% centers(1,3), &
  !                                 r, theta, phi )
  !
  !  x_ell= this% centers(1,1) &
  !         + MAX(this% sizes(1,1),this% sizes(1,2))*COS(phi)*SIN(theta)
  !
  !  y_ell= this% centers(1,2) &
  !         + MAX(this% sizes(1,3),this% sizes(1,4))*SIN(phi)*SIN(theta)
  !
  !  z_ell= this% centers(1,3) &
  !         + MAX(this% sizes(1,5),this% sizes(1,6))*COS(theta)
  !
  !  IF( r >= SQRT( ( x_ell - this% centers(1,1) )**two &
  !               + ( y_ell - this% centers(1,2) )**two &
  !               + ( z_ell - this% centers(1,3) )**two ) ) res= zero
  !
  !  IF( res < zero ) res= zero


  END PROCEDURE interpolate_fuka_mass_density


  MODULE PROCEDURE interpolate_fuka_spatial_metric

    !***********************************************
    !
    !# Returns the spatial metric.
    !
    !  FT 28.06.2022
    !
    !***********************************************

    USE constants, ONLY: pi
    USE numerics,  ONLY: trilinear_interpolation
    USE utility,   ONLY: spherical_from_cartesian, two

    IMPLICIT NONE

    INTEGER:: i_star, star

    DOUBLE PRECISION:: zp!, x_ell, y_ell, z_ell, theta, phi, r

    DOUBLE PRECISION, DIMENSION(6):: sizes
    DOUBLE PRECISION, DIMENSION(3):: center

    loop_over_stars: DO i_star= 1, 2, 1

      sizes = this% return_spatial_extent(i_star)
      center= this% return_center(i_star)

      IF( (center(1) - sizes(1) < x) .AND. (x < center(1) + sizes(2)) ) &
        star= i_star

    ENDDO loop_over_stars

    zp= z
    res= trilinear_interpolation( x, y, zp, &
                                  this% nx_grid, this% ny_grid, this% nz_grid, &
                                  this% id_fields(:,:,:,id$x:id$z,star), &
                                  this% id_fields(:,:,:,id$gxx,star), &
                                  equator_symmetry= .FALSE., debug= .FALSE. )

  END PROCEDURE interpolate_fuka_spatial_metric


  MODULE PROCEDURE interpolate_fuka_pressure

    !***********************************************
    !
    !# Returns the spatial metric.
    !
    !  FT 28.06.2022
    !
    !***********************************************

    USE constants, ONLY: pi
    USE numerics,  ONLY: trilinear_interpolation
    USE utility,   ONLY: spherical_from_cartesian, two

    IMPLICIT NONE

    INTEGER:: i_star, star

    DOUBLE PRECISION:: zp!, x_ell, y_ell, z_ell, theta, phi, r

    DOUBLE PRECISION, DIMENSION(6):: sizes
    DOUBLE PRECISION, DIMENSION(3):: center

    loop_over_stars: DO i_star= 1, 2, 1

      sizes = this% return_spatial_extent(i_star)
      center= this% return_center(i_star)

      IF( (center(1) - sizes(1) < x) .AND. (x < center(1) + sizes(2)) ) &
        star= i_star

    ENDDO loop_over_stars

    zp= z
    res= trilinear_interpolation( x, y, zp, &
                                  this% nx_grid, this% ny_grid, this% nz_grid, &
                                  this% id_fields(:,:,:,id$x:id$z,star), &
                                  this% id_fields(:,:,:,id$pressure,star), &
                                  equator_symmetry= .FALSE., debug= .FALSE. )

  END PROCEDURE interpolate_fuka_pressure


  MODULE PROCEDURE is_hydro_positive_interpolation

    !************************************************
    !
    !# Return 1 if the energy density is nonpositive
    !  or if the specific energy is nonpositive,
    !  or if the pressure is nonpositive
    !  at the specified point; return 0 otherwise
    !
    !  FT 28.06.2022
    !
    !************************************************

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(3):: center
    DOUBLE PRECISION, DIMENSION(6):: sizes

    center= this% return_center(1)
    sizes = this% return_spatial_extent(1)

    IF( this% read_mass_density( x, y, z ) <= zero &
        .OR. x > center(1) + sizes(1) &
        .OR. x < center(1) - sizes(2) &
        .OR. y > center(2) + sizes(3) &
        .OR. y < center(2) - sizes(4) &
        .OR. ABS(z) > center(3) + sizes(5) &
    )THEN
      res= .FALSE.
    ELSE
      res= .TRUE.
    ENDIF

  END PROCEDURE is_hydro_positive_interpolation


END SUBMODULE interpolate
