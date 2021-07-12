! File:         submodule_bns_mass_profile.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (bns_id) bns_mass_profile

  !*********************************************
  !                                            *
  ! Implementation of the method of TYPE bns   *
  ! that integrates the baryon mass density to *
  ! extract the radial baryon mass profile.    *
  !                                            *
  ! FT 12.07.2021                              *
  !                                            *
  !*********************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE integrate_baryon_mass_density

    !************************************************
    !                                               *
    ! Perform 3D integration over a spherical grid  *
    ! of the baryon mass density. Output baryon     *
    ! mass and radial mass profile.                 *
    !                                               *
    ! FT 19.02.2021                                 *
    !                                               *
    !************************************************

    !$ USE OMP_LIB
    USE constants, ONLY: pi
    USE NR,        ONLY: indexx

    IMPLICIT NONE

    INTEGER:: r, th, phi
    DOUBLE PRECISION:: rad_coord, colat, long, mass_element
    DOUBLE PRECISION:: g_xx, sq_g, baryon_density, gamma_euler
    !DOUBLE PRECISION:: rad

    !rad= 0.0D0

    IF(.NOT.ALLOCATED( mass_profile ))THEN
      ALLOCATE( mass_profile( 3, 0:NINT(radius/dr) ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array mass_profile in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array pos in SUBROUTINE" &
      !                // "place_particles_3D_lattice." )
    ENDIF
    IF(.NOT.ALLOCATED( mass_profile_idx ))THEN
      ALLOCATE( mass_profile_idx( 0:NINT(radius/dr) ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array mass_profile in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array pos in SUBROUTINE" &
      !                // "place_particles_3D_lattice." )
    ENDIF

    mass_profile( 1, 0 )= 0.0D0
    mass_profile( 2, 0 )= 4.0D0/3.0D0*pi*dr**3.0D0*central_density
    mass_profile( 3, 0 )= 4.0D0/3.0D0*pi*dr**3.0D0*central_density

    !$OMP PARALLEL DO SHARED(dr,dphi,dth,center)&
    !$OMP             PRIVATE(r,th,phi,rad_coord,long,colat,sq_g,gamma_euler, &
    !$OMP                     baryon_density,mass_element,mass) !&
  !    !$OMP             SCHEDULE(STATIC,1)
    radius_loop: DO r= 1, NINT(radius/dr), 1

      mass= 0.0D0
      rad_coord= r*dr

  !      !$OMP PARALLEL SECTIONS REDUCTION(+:mass,rad)
      !rad= rad + dr

      longitude_loop: DO phi= 1, NINT(2.0D0*pi/dphi), 1

        long= phi*dphi

        colatitude_loop: DO th= 1, NINT(pi/2.0D0/dth), 1

          colat= th*dth

          ! The definition of the baryon mass for the LORENE ID is in eq.(69)
          ! of Gourgoulhon et al., PRD 63 064029 (2001)

          CALL THIS% import_id( &
                   center + (rad_coord + dr)*SIN(colat)*COS(long), &
                   (rad_coord + dr)*SIN(colat)*SIN(long), &
                   (rad_coord + dr)*COS(colat), &
                   g_xx, baryon_density, gamma_euler )

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
  !        !! Compute vector normal to spacelie hypersurface
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

          ! Compute square root of the determinant of the spatial metric
          sq_g= g_xx*SQRT( g_xx )

          mass_element= (rad_coord**2.0D0)*SIN(colat)*dr*dth*dphi &
                        *sq_g*gamma_euler*baryon_density

          mass= mass + 2.0D0*mass_element

        ENDDO colatitude_loop

      ENDDO longitude_loop
  !      !$OMP END PARALLEL SECTIONS

      mass_profile( 1, r )= rad_coord
      mass_profile( 2, r )= mass

      !PRINT *, rad_coord, mass_profile( 1, r ), mass, mass_profile( 2, r )
      !PRINT *

    ENDDO radius_loop
    !$OMP END PARALLEL DO

    DO r= 1, NINT(radius/dr), 1
      mass_profile( 3, r )= mass_profile( 3, r - 1 ) + mass_profile( 2, r )
    ENDDO

    mass= mass_profile( 3, NINT(radius/dr) )

    PRINT *, "radius covered by the integration=", MAXVAL( mass_profile( 1, : ), DIM= 1 )
    PRINT *, "mass=", mass
    PRINT *

    CALL indexx( NINT(radius/dr) + 1, mass_profile( 1, : ), mass_profile_idx )


  END PROCEDURE integrate_baryon_mass_density


END SUBMODULE bns_mass_profile
