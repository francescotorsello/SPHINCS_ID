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
  !# Implementation of the methods of TYPE [[bnsfuka]] that
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
    INTEGER:: a, star
    DOUBLE PRECISION:: zp
    !DOUBLE PRECISION, &
    !  DIMENSION(this% nx_grid, this% ny_grid, this% nz_grid, 3):: coords

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( n, this, x, y, z, lapse, &
    !$OMP                     shift_x, shift_y, shift_z, &
    !$OMP                     g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, &
    !$OMP                     baryon_density, energy_density, &
    !$OMP                     specific_energy, pressure, &
    !$OMP                     u_euler_x, u_euler_y, u_euler_z ) &
    !$OMP             PRIVATE( a, star, zp )
    DO a= 1, n, 1

      IF( (this% center(1,1) - this% radii(1,1) <= x(a)) &
          .AND. &
          (x(a) <= this% center(1,1) + this% radii(1,2)) )THEN

        star= 1

      ELSEIF( (this% center(2,1) - this% radii(2,1) <= x(a)) &
              .AND. &
              (x(a) <= this% center(2,1) + this% radii(2,2)) )THEN

        star= 2

      ELSE

        star= -1

      ENDIF

      IF( star == -1 )THEN
        baryon_density(a) = zero
        specific_energy(a)= zero
        pressure(a)       = zero
        u_euler_x(a)      = zero
        u_euler_y(a)      = zero
        u_euler_z(a)      = zero
        energy_density(a) = zero
        g_xx(a)           = zero
        g_yy(a)           = zero
        g_zz(a)           = zero
        g_xy(a)           = zero
        g_xz(a)           = zero
        g_yz(a)           = zero
        lapse(a)          = zero
        shift_x(a)        = zero
        shift_y(a)        = zero
        shift_z(a)        = zero
        CYCLE
      ENDIF

      !coords= this% star_lattice(star)% coords

      zp= z(a)

      baryon_density(a) = trilinear_interpolation( x(a), y(a), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% star_lattice(star)% coords, &
                                this% star_lattice(star)% mass_density, &
                                debug= .TRUE. ) &
                                *MSun/amu

      specific_energy(a)= trilinear_interpolation( x(a), y(a), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% star_lattice(star)% coords, &
                                this% star_lattice(star)% specific_energy )

      pressure(a)       = trilinear_interpolation( x(a), y(a), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% star_lattice(star)% coords, &
                                this% star_lattice(star)% pressure ) &
                                *MSun/amu

      u_euler_x(a)      = trilinear_interpolation( x(a), y(a), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% star_lattice(star)% coords, &
                                this% star_lattice(star)% v_eul_x )
      u_euler_y(a)      = trilinear_interpolation( x(a), y(a), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% star_lattice(star)% coords, &
                                this% star_lattice(star)% v_eul_y )
      u_euler_z(a)      = trilinear_interpolation( x(a), y(a), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% star_lattice(star)% coords, &
                                this% star_lattice(star)% v_eul_z )

      IF( baryon_density(a) == zero )THEN
        specific_energy(a)= zero
        u_euler_x(a)      = zero
        u_euler_y(a)      = zero
        u_euler_z(a)      = zero
      ENDIF

      energy_density(a) = baryon_density(a)*(one + specific_energy(a))

      g_xx(a)           = trilinear_interpolation( x(a), y(a), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% star_lattice(star)% coords, &
                                this% star_lattice(star)% g_xx )

      g_yy(a)= g_xx(a)
      g_zz(a)= g_xx(a)
      g_xy(a)= zero
      g_xz(a)= zero
      g_yz(a)= zero

      lapse(a)          = trilinear_interpolation( x(a), y(a), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% star_lattice(star)% coords, &
                                this% star_lattice(star)% lapse )
      shift_x(a)        = trilinear_interpolation( x(a), y(a), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% star_lattice(star)% coords, &
                                this% star_lattice(star)% shift_x )
      shift_y(a)        = trilinear_interpolation( x(a), y(a), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% star_lattice(star)% coords, &
                                this% star_lattice(star)% shift_y )
      shift_z(a)        = trilinear_interpolation( x(a), y(a), zp, &
                                this% nx_grid, this% ny_grid, this% nz_grid, &
                                this% star_lattice(star)% coords, &
                                this% star_lattice(star)% shift_z )

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

    INTEGER:: star
    DOUBLE PRECISION:: zp, veuler_x, veuler_y, veuler_z
    !DOUBLE PRECISION, &
    !  DIMENSION(this% nx_grid, this% ny_grid, this% nz_grid, 3):: coords

    IF( (this% center(1,1) - this% radii(1,1) <= x) &
        .AND. &
        (x <= this% center(1,1) + this% radii(1,2)) )THEN

      star= 1

    ELSEIF( (this% center(2,1) - this% radii(2,1) <= x) &
            .AND. &
            (x <= this% center(2,1) + this% radii(2,2)) )THEN

      star= 2

    ELSE

      star= -1

    ENDIF

    IF( star == -1 )THEN
      baryon_density= zero
      g(jxx)        = zero
      g(jyy)        = zero
      g(jzz)        = zero
      g(jxy)        = zero
      g(jxz)        = zero
      g(jyz)        = zero
      gamma_euler   = zero
      RETURN
    ENDIF

    !coords= this% star_lattice(star)% coords

    zp= z

    baryon_density= trilinear_interpolation( x, y, zp, &
                          this% nx_grid, this% ny_grid, this% nz_grid, &
                          this% star_lattice(star)% coords, &
                          this% star_lattice(star)% mass_density )

    g(jxx)= trilinear_interpolation( x, y, zp, &
                  this% nx_grid, this% ny_grid, this% nz_grid, &
                  this% star_lattice(star)% coords, &
                  this% star_lattice(star)% g_xx )
    g(jyy)= g(jxx)
    g(jzz)= g(jxx)
    g(jxy)= zero
    g(jxz)= zero
    g(jyz)= zero

    veuler_x= trilinear_interpolation( x, y, zp, &
                     this% nx_grid, this% ny_grid, this% nz_grid, &
                     this% star_lattice(star)% coords, &
                     this% star_lattice(star)% v_eul_x )
    veuler_y= trilinear_interpolation( x, y, zp, &
                     this% nx_grid, this% ny_grid, this% nz_grid, &
                     this% star_lattice(star)% coords, &
                     this% star_lattice(star)% v_eul_y )
    veuler_z= trilinear_interpolation( x, y, zp, &
                     this% nx_grid, this% ny_grid, this% nz_grid, &
                     this% star_lattice(star)% coords, &
                     this% star_lattice(star)% v_eul_z )

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

    USE timing,    ONLY: timer
    USE numerics,  ONLY: trilinear_interpolation
    USE utility,   ONLY: spherical_from_cartesian, two


    IMPLICIT NONE

    INTEGER:: star
    DOUBLE PRECISION:: zp

    !TYPE(timer):: dbg_timer

    !dbg_timer= timer("dbg_timer")
    !CALL dbg_timer% start_timer()

    IF( (this% center(1,1) - this% radii(1,1) <= x) &
        .AND. &
        (x <= this% center(1,1) + this% radii(1,2)) )THEN

      star= 1

    ELSEIF( (this% center(2,1) - this% radii(2,1) <= x) &
            .AND. &
            (x <= this% center(2,1) + this% radii(2,2)) )THEN

      star= 2

    ELSE

      star= -1

    ENDIF

    IF( star == -1 )THEN
      res= zero
      RETURN
    ENDIF

    zp= z

    res= trilinear_interpolation( x, y, zp, &
                                  this% nx_grid, this% ny_grid, this% nz_grid, &
                                  this% star_lattice(star)% coords, &
                                  this% star_lattice(star)% mass_density )
    !CALL dbg_timer% stop_timer()
    !CALL dbg_timer% print_timer( 2 )
    !STOP
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

    INTEGER:: star
    DOUBLE PRECISION:: zp

    IF( (this% center(1,1) - this% radii(1,1) <= x) &
        .AND. &
        (x <= this% center(1,1) + this% radii(1,2)) )THEN

      star= 1

    ELSEIF( (this% center(2,1) - this% radii(2,1) <= x) &
            .AND. &
            (x <= this% center(2,1) + this% radii(2,2)) )THEN

      star= 2

    ELSE

      star= -1

    ENDIF

    IF( star == -1 )THEN
      res= zero
      RETURN
    ENDIF

    zp= z
    res= trilinear_interpolation( x, y, zp, &
                                  this% nx_grid, this% ny_grid, this% nz_grid, &
                                  this% star_lattice(star)% coords, &
                                  this% star_lattice(star)% g_xx )

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

    INTEGER:: star
    DOUBLE PRECISION:: zp

    IF( (this% center(1,1) - this% radii(1,1) <= x) &
        .AND. &
        (x <= this% center(1,1) + this% radii(1,2)) )THEN

      star= 1

    ELSEIF( (this% center(2,1) - this% radii(2,1) <= x) &
            .AND. &
            (x <= this% center(2,1) + this% radii(2,2)) )THEN

      star= 2

    ELSE

      star= -1

    ENDIF

    IF( star == -1 )THEN
      res= zero
      RETURN
    ENDIF

    zp= z
    res= trilinear_interpolation( x, y, zp, &
                                  this% nx_grid, this% ny_grid, this% nz_grid, &
                                  this% star_lattice(star)% coords, &
                                  this% star_lattice(star)% pressure )

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

    LOGICAL:: outside_y, outside_z!, outside_x

    !outside_x=.NOT.( &
    !            (
    !              this% center(1,1) - this% radii(1,1) <= x &
    !              .AND. &
    !              x <= this% center(1,1) + this% radii(1,2) &
    !            ) &
    !            .OR. &
    !            ( &
    !              this% center(2,1) - this% radii(2,1) <= x &
    !              .AND. &
    !              x <= this% center(2,1) + this% radii(2,2) &
    !            ) &
    !          )

    outside_y=.NOT.( &
                ( &
                  this% center(1,2) - this% radii(1,3) <= y &
                  .AND. &
                  y <= this% center(1,2) + this% radii(1,4) &
                ) &
                .OR. &
                ( &
                  this% center(2,2) - this% radii(2,3) <= y &
                  .AND. &
                  y <= this% center(2,2) + this% radii(2,4) &
                ) &
              )

    outside_z=.NOT.( &
                ( &
                  this% center(1,3) - this% radii(1,5) <= z &
                  .AND. &
                  z <= this% center(1,3) + this% radii(1,6) &
                ) &
                .OR. &
                ( &
                  this% center(2,3) - this% radii(2,5) <= z &
                  .AND. &
                  z <= this% center(2,3) + this% radii(2,6) &
                ) &
              )

    IF( this% read_mass_density( x, y, z ) <= zero &
        !.OR. outsize_x &
        .OR. outside_y &
        .OR. outside_z &
    )THEN
      res= .FALSE.
    ELSE
      res= .TRUE.
    ENDIF

  END PROCEDURE is_hydro_positive_interpolation


END SUBMODULE interpolate
