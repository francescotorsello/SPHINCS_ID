! File:         submodule_id_base_mass_profile.f90
! Authors:      Francesco Torsello (FT)
!************************************************************************
! Copyright (C) 2020-2023 Francesco Torsello                            *
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

SUBMODULE (id_base) mass_profile

  !********************************************
  !
  !# Implementation of the method of TYPE idbase
  !  that integrates the baryon mass density to
  !  extract the radial baryon mass profile.
  !
  !  FT 12.07.2021
  !
  !********************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE integrate_baryon_mass_density

    !************************************************
    !
    !# Perform 3D integration over a spherical grid
    !  of the baryon mass density. Output baryon
    !  mass and radial mass profile.
    !
    !  FT 19.02.2021
    !
    !  Upgraded to ellipsoidal grid
    !
    !  FT 15.11.2022
    !
    !************************************************

    USE constants, ONLY: pi
    USE utility,   ONLY: zero, one, two, three, four
    USE NR,        ONLY: indexx
    USE tensor,    ONLY: jxx, jxy, jxz, jyy, jyz, jzz
    USe utility,   ONLY: determinant_sym3x3, cartesian_from_spherical
    USE numerics,  ONLY: bilinear_interpolation

    IMPLICIT NONE

    INTEGER:: r, th, phi
    DOUBLE PRECISION:: rad, rad_coord, colat, long, mass_element, max_radius
    DOUBLE PRECISION:: sq_g, baryon_density, gamma_euler
    DOUBLE PRECISION:: a_x, a_y, a_z, xtemp, ytemp, ztemp
    DOUBLE PRECISION, DIMENSION(6):: g

    CHARACTER(LEN=:), ALLOCATABLE:: surface_type

    LOGICAL, PARAMETER:: debug= .FALSE.

    a_x= one
    a_y= one
    a_z= one
    max_radius= radius
    surface_type= "spherical"
    IF(PRESENT(surf))THEN

      IF(surf% is_known)THEN

        surface_type= "oval"

      ENDIF

    ENDIF

    IF(PRESENT(radii))THEN

      IF(.NOT.PRESENT(surf)) surface_type= "ellipsoidal"

      IF(PRESENT(surf))THEN

        IF(.NOT.(surf% is_known))THEN

          surface_type= "ellipsoidal"

        ENDIF

      ENDIF

      max_radius= MAXVAL([radius,radii(1),radii(2)])

      IF(max_radius == radius)THEN

        a_x= one
        a_y= radii(1)/max_radius
        a_z= radii(2)/max_radius

      ELSEIF(max_radius == radii(1))THEN

        a_x= radius/max_radius
        a_y= one
        a_z= radii(2)/max_radius

      ELSEIF(max_radius == radii(2))THEN

        a_x= radius/max_radius
        a_y= radii(1)/max_radius
        a_z= one

      ENDIF

    ENDIF


    IF(debug) PRINT *, "inside integrate_baryon_mass_density, ", &
                       "surface_type is:", surface_type
   ! STOP

    mass_profile( 1, 0 )= zero
    mass_profile( 2, 0 )= four/three*pi*dr**three*central_density
    mass_profile( 3, 0 )= four/three*pi*dr**three*central_density

    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP             SHARED(dr, dphi, dth, center, max_radius, &
    !$OMP                    mass_profile, this, a_x, a_y, a_z, surf) &
    !$OMP             PRIVATE(r, th, phi, rad, rad_coord, long, colat, sq_g, &
    !$OMP                     gamma_euler, g, baryon_density, mass_element, &
    !$OMP                     mass, xtemp, ytemp, ztemp)
    radius_loop: DO r= 1, NINT(max_radius/dr), 1

      mass= zero
      rad_coord= r*dr

      longitude_loop: DO phi= 1, NINT(two*pi/dphi), 1

        long= phi*dphi

        colatitude_loop: DO th= 1, NINT(pi/two/dth), 1

          colat= th*dth

          ! The definition of the baryon mass for the LORENE ID is in eq.(69)
          ! of Gourgoulhon et al., PRD 63 064029 (2001)
          rad= rad_coord
          IF(PRESENT(surf))THEN

            IF(surf% is_known)THEN

              rad= bilinear_interpolation( colat, long, &
                    SIZE(surf% points(:,1,5)), &
                    SIZE(surf% points(1,:,6)), &
                    surf% points(:,:,5:6), surf% points(:,:,4) )
              rad= rad*rad_coord/max_radius

            ENDIF

          ENDIF

          CALL cartesian_from_spherical( &
            a_x*(rad + dr), colat, long, &
            center(1), center(2), center(3), &
            xtemp, ytemp, ztemp, a_y/a_x, a_z/a_x )

          !CALL this% read_id_mass_b( &
          !             center(1) + a_x*(rad_coord + dr)*SIN(colat)*COS(long), &
          !             center(2) + a_y*(rad_coord + dr)*SIN(colat)*SIN(long), &
          !             center(3) + a_z*(rad_coord + dr)*COS(colat),           &
          !             g, baryon_density, gamma_euler )

          CALL this% read_id_mass_b( xtemp, ytemp, ztemp, &
                                     g, baryon_density, gamma_euler )

          IF(      ISNAN( g(jxx) ) .OR. ISNAN( g(jxy) ) .OR. ISNAN( g(jxz) ) &
              .OR. ISNAN( g(jyy) ) .OR. ISNAN( g(jyz) ) .OR. ISNAN( g(jzz) ) &
              .OR. ISNAN( baryon_density ) .OR. ISNAN( gamma_euler ) ) &
            CYCLE

          ! Compute square root of the determinant of the spatial metric
          CALL determinant_sym3x3(g, sq_g)
          sq_g= SQRT(sq_g)

          !mass_element= a_x*a_y*a_z*(rad_coord**two)*SIN(colat)*dr*dth*dphi &
          !              *sq_g*gamma_euler*baryon_density

          mass_element= a_x*a_y*a_z*(rad**two)*SIN(colat)*dr*dth*dphi &
                        *sq_g*gamma_euler*baryon_density

          mass= mass + two*mass_element

        ENDDO colatitude_loop

      ENDDO longitude_loop

      mass_profile(1, r)= rad_coord
      mass_profile(2, r)= mass

    ENDDO radius_loop
    !$OMP END PARALLEL DO

    DO r= 1, NINT(radius/dr), 1
      mass_profile(3, r)= mass_profile(3, r - 1) + mass_profile(2, r)
    ENDDO

    mass= mass_profile(3, NINT(radius/dr))

    IF( ISNAN(mass) )THEN
      PRINT *, "** ERROR! The integrated mass is a NaN!"
      PRINT *
      STOP
    ELSEIF( mass <= 0 )THEN
      PRINT *, "** ERROR! The integrated mass is mass=", mass, "<= 0!"
      PRINT *
      STOP
    ENDIF

    PRINT *, " * Radius covered by the integration of baryon mass density=", &
             MAXVAL( mass_profile( 1, : ), DIM= 1 )
    PRINT *, " * Integrated baryon mass of the star=", mass
    PRINT *

    CALL indexx( NINT(radius/dr) + 1, mass_profile( 1, : ), mass_profile_idx )


  END PROCEDURE integrate_baryon_mass_density


END SUBMODULE mass_profile


! TODO:  deprecated? It computes the relative Lorentz factor between the fluid
!        and the Eulerian observer from the velocity
!
!        CALL bns_obj% import_id( &
!                 center1 + rad_coord*SIN(lat)*COS(long), &
!                 rad_coord*SIN(lat)*SIN(long), &
!                 rad_coord*COS(lat), &
!                 g_xx, baryon_density, &
!                 gamma_euler )
!
!        ! Compute covariant spatial fluid velocity (metric is diagonal and
!        ! conformally flat)
!        !v_euler_x_l= g_xx*v_euler_x
!        !v_euler_y_l= g_xx*v_euler_y
!        !v_euler_z_l= g_xx*v_euler_z
!        !
!        !! Compute the corresponding Lorentz factor
!        !lorentz_factor= 1.0D0/SQRT( 1.0D0 - ( v_euler_x_l*v_euler_x &
!        !                                    + v_euler_y_l*v_euler_y &
!        !                                    + v_euler_z_l*v_euler_z ) )
!        !
!        !! Compute covariant fluid 4-velocity
!        !u_euler_t_l= lorentz_factor *( - lapse + v_euler_x_l*shift_x &
!        !                                       + v_euler_y_l*shift_y &
!        !                                       + v_euler_z_l*shift_z )
!        !u_euler_x_l= lorentz_factor*v_euler_x_l
!        !u_euler_y_l= lorentz_factor*v_euler_y_l
!        !u_euler_z_l= lorentz_factor*v_euler_z_l
!        !
!        !! Compute vector normal to spacelike hypersurface
!        !! (4-velocity of the Eulerian observer)
!        !n_t= 1.0D0/lapse
!        !n_x= - shift_x/lapse
!        !n_y= - shift_y/lapse
!        !n_z= - shift_z/lapse
!        !
!        !! Compute relative Lorentz factor between 4-velocity of the fluid
!        !! wrt the Eulerian observer and the 4-velocity of the Eulerian observer
!        !lorentz_factor_rel= - ( n_t*u_euler_t_l + n_x*u_euler_x_l &
!        !                      + n_y*u_euler_y_l + n_z*u_euler_z_l )
