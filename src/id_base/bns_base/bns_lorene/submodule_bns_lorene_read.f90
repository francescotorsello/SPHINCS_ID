! File:         submodule_bnslorene_read.f90
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

SUBMODULE (bns_lorene) read

  !****************************************************
  !
  !# Implementation of the methods of TYPE [[bnslorene]] that
  !  read |bns| data using |lorene|
  !
  !  FT 23.10.2020
  !
  !  Renamed from methods to read upon
  !  improving modularity
  !
  !  OMP parallelized loops that call |lorene|
  !  in all MODULE PROCEDURE
  !
  !  FT 12.07.2021
  !
  !****************************************************


  USE utility,  ONLY: zero, one, two, Msun_geo


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE read_id_int

    !**************************************************
    !
    !# Stores the ID in the [[bnslorene]] member arrays
    !
    !  FT 5.10.2020
    !
    !**************************************************

    IMPLICIT NONE

    IF ( C_ASSOCIATED( this% bns_ptr ) ) THEN

      IF( SIZE( x ) /= SIZE( y ) .OR. SIZE( x ) /= SIZE( z ) &
            .OR. SIZE( y ) /= SIZE( z ) )THEN
        PRINT *, "** ERROR: The sizes of the arrays of positions" &
                 // "passed to read_id are not the same."
        PRINT *
        STOP
      ENDIF

      IF( ALLOCATED( this% lapse )   .AND. &
          ALLOCATED( this% shift_x ) .AND. &
          ALLOCATED( this% shift_y ) .AND. &
          ALLOCATED( this% shift_z ) .AND. &
          ALLOCATED( this% g_xx ) .AND. ALLOCATED( this% g_xy ) .AND. &
          ALLOCATED( this% g_xz ) .AND. ALLOCATED( this% g_yy ) .AND. &
          ALLOCATED( this% g_yz ) .AND. ALLOCATED( this% g_zz ) .AND. &
          ALLOCATED( this% k_xx ) .AND. ALLOCATED( this% k_xy ) .AND. &
          ALLOCATED( this% k_xz ) .AND. ALLOCATED( this% k_yy ) .AND. &
          ALLOCATED( this% k_yz ) .AND. ALLOCATED( this% k_zz ) .AND. &
          ALLOCATED( this% baryon_density )  .AND. &
          ALLOCATED( this% energy_density )  .AND. &
          ALLOCATED( this% specific_energy ) .AND. &
          ALLOCATED( this% v_euler_x )       .AND. &
          ALLOCATED( this% v_euler_y )       .AND. &
          ALLOCATED( this% v_euler_z ) &
      )THEN

        !$OMP PARALLEL DO DEFAULT( NONE ) &
        !$OMP             SHARED( n, this, x, y, z ) &
        !$OMP             PRIVATE( itr )
        read_id_loop: DO itr= 1, n, 1

          ! The coordinates need to be converted from |sphincs| units (Msun_geo)
          ! to |lorene| units (\(\mathrm{km}\)). See MODULE constants for the
          ! definition of Msun_geo
          CALL get_lorene_id( this% bns_ptr, &
                              x( itr )*Msun_geo, &
                              y( itr )*Msun_geo, &
                              z( itr )*Msun_geo, &
                              this% lapse( itr ), &
                              this% shift_x( itr ), &
                              this% shift_y( itr ), &
                              this% shift_z( itr ), &
                              this% g_xx( itr ), &
                              this% k_xx( itr ), &
                              this% k_xy( itr ), &
                              this% k_xz( itr ), &
                              this% k_yy( itr ), &
                              this% k_yz( itr ), &
                              this% k_zz( itr ), &
                              this% baryon_density( itr ), &
                              this% energy_density( itr ), &
                              this% specific_energy( itr ), &
                              this% pressure( itr ), &
                              this% v_euler_x( itr ), &
                              this% v_euler_y( itr ), &
                              this% v_euler_z( itr ) )

        ENDDO read_id_loop
        !$OMP END PARALLEL DO

        DO itr= 1, n, 1

          !
          !-- The following follows from the assumption of conformal
          !-- flatness in |lorene|
          !
          this% g_yy( itr )= this% g_xx( itr )
          this% g_zz( itr )= this% g_xx( itr )
          this% g_xy( itr )= zero
          this% g_xz( itr )= zero
          this% g_yz( itr )= zero

          !
          !- Set/unset the geodesic gauge
          !
          IF( this% get_one_lapse() )THEN
            this% lapse( itr )= one
          ENDIF
          IF( this% get_zero_shift() )THEN
            this% shift_x( itr )= zero
            this% shift_y( itr )= zero
            this% shift_z( itr )= zero
          ENDIF

          !
          !-- Convert the extrinsic curvature from |lorene| units to
          !-- |sphincs| units
          !
          this% k_xx( itr )= this% k_xx( itr )*Msun_geo
          this% k_xy( itr )= this% k_xy( itr )*Msun_geo
          this% k_xz( itr )= this% k_xz( itr )*Msun_geo
          this% k_yy( itr )= this% k_yy( itr )*Msun_geo
          this% k_yz( itr )= this% k_yz( itr )*Msun_geo
          this% k_zz( itr )= this% k_zz( itr )*Msun_geo

          ! Print progress on screen
          perc= 100*itr/n
          IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
            WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
                    creturn//" ", perc, "%"
          ENDIF

        ENDDO
        IF( show_progress ) WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

      ELSE

        PRINT *, "** ERROR: Memory was not allocated before calling " &
                 // "read_id in read_id (TYPE particles)."
        PRINT *
        STOP

      ENDIF

      PRINT *, "** Subroutine read_id executed."
      PRINT *

    ENDIF

  END PROCEDURE read_id_int


  MODULE PROCEDURE read_id_full

    !**************************************************
    !
    !# Stores the ID in non-[[bnslorene]]-member arrays
    !  with the same shape as the [[bnslorene]] member arrays
    !
    !  FT 5.10.2020
    !
    !**************************************************

    IMPLICIT NONE

    IF ( C_ASSOCIATED( this% bns_ptr ) ) THEN

      IF( SIZE( x ) /= SIZE( y ) .OR. SIZE( x ) /= SIZE( z ) &
            .OR. SIZE( y ) /= SIZE( z ) )THEN
        PRINT *, "** ERROR: The sizes of the arrays of positions" &
                 // "passed to read_id are not the same."
        PRINT *
        STOP
      ENDIF

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( n, this, x, y, z, lapse, &
      !$OMP                     shift_x, shift_y, shift_z, &
      !$OMP                     g_xx, k_xx, k_xy, k_xz, k_yy, k_yz, k_zz, &
      !$OMP                     baryon_density, energy_density, &
      !$OMP                     specific_energy, pressure, &
      !$OMP                     u_euler_x, u_euler_y, u_euler_z ) &
      !$OMP             PRIVATE( itr )
      read_id_loop: DO itr= 1, n, 1

        ! The coordinates need to be converted from |sphincs| units (Msun_geo)
        ! to |lorene| units (\(\mathrm{km}\)). See MODULE constants for the definition of
        ! Msun_geo
        CALL get_lorene_id( this% bns_ptr, &
                            x( itr )*Msun_geo, &
                            y( itr )*Msun_geo, &
                            z( itr )*Msun_geo, &
                            lapse( itr ), &
                            shift_x( itr ), &
                            shift_y( itr ), &
                            shift_z( itr ), &
                            g_xx( itr ), &
                            k_xx( itr ), &
                            k_xy( itr ), &
                            k_xz( itr ), &
                            k_yy( itr ), &
                            k_yz( itr ), &
                            k_zz( itr ), &
                            baryon_density( itr ), &
                            energy_density( itr ), &
                            specific_energy( itr ), &
                            pressure( itr ), &
                            u_euler_x( itr ), &
                            u_euler_y( itr ), &
                            u_euler_z( itr ) )

      ENDDO read_id_loop
      !$OMP END PARALLEL DO

      DO itr= 1, n, 1

        !
        !-- The following follows from the assumption of conformal
        !-- flatness in |lorene|
        !
        g_yy( itr )= g_xx( itr )
        g_zz( itr )= g_xx( itr )
        g_xy( itr )= zero
        g_xz( itr )= zero
        g_yz( itr )= zero

        !
        !- Set/unset the geodesic gauge
        !
        IF( this% get_one_lapse() )THEN
          lapse( itr )= one
        ENDIF
        IF( this% get_zero_shift() )THEN
          shift_x( itr )= zero
          shift_y( itr )= zero
          shift_z( itr )= zero
        ENDIF

        !
        !-- Convert the extrinsic curvature from |lorene| units to
        !-- |sphincs| units
        !
        k_xx( itr )= k_xx( itr )*Msun_geo
        k_xy( itr )= k_xy( itr )*Msun_geo
        k_xz( itr )= k_xz( itr )*Msun_geo
        k_yy( itr )= k_yy( itr )*Msun_geo
        k_yz( itr )= k_yz( itr )*Msun_geo
        k_zz( itr )= k_zz( itr )*Msun_geo

        ! Print progress on screen
        perc= 100*itr/n
        IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
          WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
                  creturn//" ", perc, "%"
        ENDIF

      ENDDO
      IF( show_progress ) WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

      PRINT *, "** Subroutine read_id executed."
      PRINT *

    ENDIF

  END PROCEDURE read_id_full


  MODULE PROCEDURE read_id_spacetime

    !*******************************************************
    !
    !# Stores the spacetime ID in multi-dimensional arrays
    !  needed to compute the BSSN variables and constraints
    !
    !  FT 22.11.2020
    !
    !*******************************************************

    USE tensor,    ONLY: jxx, jxy, jxz, &
                         jyy, jyz, jzz, jx, jy, jz, n_sym4x4

    IMPLICIT NONE

    INTEGER:: i, j, k

    DOUBLE PRECISION:: detg
    DOUBLE PRECISION:: detg4
    DOUBLE PRECISION, DIMENSION( :, :, :, : ), ALLOCATABLE:: g4

    ! g4 is allocatable to allocate it on the heap
    ! Allocating it on the stack might exceed stack memory,
    ! causing a segmentation fault
    ALLOCATE( g4( nx, ny, nz, n_sym4x4 ) )

    IF ( C_ASSOCIATED( this% bns_ptr ) ) THEN

      IF( .FALSE. &!SHAPE( pos(:,:,:,1) ) /= SHAPE( lapse ) .OR. &
          !SHAPE( pos(:,:,:,1) ) /= SHAPE( shift(:,:,:,jx) ) & ! .OR. &
        ! SHAPE( pos(:,:,:,1) ) /= SHAPE( g(:,:,:,1) ) .OR. &
        ! SHAPE( pos(:,:,:,1) ) /= SHAPE( k(:,:,:,1) ) &
        )THEN
        PRINT *, "** ERROR: Mismatch in array dimensions" &
                 // "in read_id_spacetime."
        PRINT *
        STOP
      ENDIF

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( nx, ny, nz, this, pos, &
      !$OMP                     lapse, shift, g, ek ) &
      !$OMP             PRIVATE( i, j, k )
      coords_z: DO k= 1, nz, 1
        coords_y: DO j= 1, ny, 1
          coords_x: DO i= 1, nx, 1

            ! The coordinates need to be converted from |sphincs| units (Msun_geo)
            ! to |lorene| units (\(\mathrm{km}\)). See MODULE constants for the definition of
            ! Msun_geo
            CALL get_lorene_id_spacetime( this% bns_ptr, &
                                pos( i, j, k, jx )*Msun_geo, &
                                pos( i, j, k, jy )*Msun_geo, &
                                pos( i, j, k, jz )*Msun_geo, &
                                lapse( i, j, k ), &
                                shift( i, j, k, jx ), &
                                shift( i, j, k, jy ), &
                                shift( i, j, k, jz ), &
                                g( i, j, k, jxx ), &
                                ek( i, j, k, jxx ), &
                                ek( i, j, k, jxy ), &
                                ek( i, j, k, jxz ), &
                                ek( i, j, k, jyy ), &
                                ek( i, j, k, jyz ), &
                                ek( i, j, k, jzz ) )

          ENDDO coords_x
        ENDDO coords_y
      ENDDO coords_z
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( nx, ny, nz, this, pos, &
      !$OMP                     lapse, shift, g, ek, g4 ) &
      !$OMP             PRIVATE( i, j, k, detg, detg4 )
      DO k= 1, nz, 1
        DO j= 1, ny, 1
          DO i= 1, nx, 1

            !
            !-- The following follows from the assumption of
            !-- conformal flatness in |lorene|
            !
            g( i, j, k, jyy )= g( i, j, k, jxx )
            g( i, j, k, jzz )= g( i, j, k, jxx )
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
            ek( i, j, k, jxx )= ek( i, j, k, jxx )*Msun_geo
            ek( i, j, k, jxy )= ek( i, j, k, jxy )*Msun_geo
            ek( i, j, k, jxz )= ek( i, j, k, jxz )*Msun_geo
            ek( i, j, k, jyy )= ek( i, j, k, jyy )*Msun_geo
            ek( i, j, k, jyz )= ek( i, j, k, jyz )*Msun_geo
            ek( i, j, k, jzz )= ek( i, j, k, jzz )*Msun_geo

            detg= 2.0D0*g(i,j,k,jxy)*g(i,j,k,jxz)*g(i,j,k,jyz) &
                  - g(i,j,k,jzz)*g(i,j,k,jxy)**2 + g(i,j,k,jyy) &
                   *( g(i,j,k,jxx)*g(i,j,k,jzz) - g(i,j,k,jxz)**2 ) &
                  - g(i,j,k,jxx)*g(i,j,k,jyz)**2

            IF( ABS( detg ) < 1D-10 )THEN
              PRINT *, "The determinant of the spatial metric " &
                       // "is effectively 0 at the grid point " &
                       // "(ix,iy,iz)= (", i, ",", j,",", k, ")."
              PRINT *, "detg=", detg
              PRINT *
              STOP
            ELSEIF( detg < 0 )THEN
              PRINT *, "The determinant of the spatial metric " &
                       // "is negative at the grid point " &
                       // "(ix,iy,iz)= (", i, ",", j,",", k, ")."
              PRINT *, "detg=", detg
              PRINT *
              STOP
            ENDIF

            CALL compute_g4( lapse(i,j,k), shift(i,j,k,:), &
                             g(i,j,k,:), g4(i,j,k,:) )

            CALL determinant_sym4x4( g4(i,j,k,:), detg4 )

            IF( ABS( detg4 ) < 1D-10 )THEN
              PRINT *, "The determinant of the spacetime metric "&
                       // "is effectively 0 at the grid point " &
                       // "(ix,iy,iz)= (", i, ",", j,",", k, ")."
              PRINT *, "detg4=", detg4
              PRINT *
              STOP
            ELSEIF( detg4 > 0 )THEN
              PRINT *, "The determinant of the spacetime metric "&
                       // "is positive at the grid point " &
                       // "(ix,iy,iz)= (", i, ",", j,",", k, ")."
              PRINT *, "detg4=", detg4
              PRINT *
              STOP
            ENDIF

          ENDDO
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

      PRINT *, "** Subroutine read_id executed."
      PRINT *

    ENDIF

  END PROCEDURE read_id_spacetime


  MODULE PROCEDURE read_id_hydro

    !*******************************************************
    !
    !# Stores the hydro ID in the arrays needed to compute
    !  the constraints on the refined mesh
    !
    !  FT 25.11.2020
    !
    !*******************************************************

    USE tensor,     ONLY: jx, jy, jz

    IMPLICIT NONE

    INTEGER:: ix, iy, iz

    IF ( C_ASSOCIATED( this% bns_ptr ) ) THEN

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( nx, ny, nz, this, pos, &
      !$OMP                     baryon_density, energy_density, &
      !$OMP                     specific_energy, pressure, u_euler ) &
      !$OMP             PRIVATE( ix, iy, iz )
      coords_z: DO iz= 1, nz, 1
        coords_y: DO iy= 1, ny, 1
          coords_x: DO ix= 1, nx, 1

            ! The coordinates need to be converted from |sphincs| units (Msun_geo)
            ! to |lorene| units (\(\mathrm{km}\)). See MODULE constants for the definition of
            ! Msun_geo
            CALL get_lorene_id_hydro( this% bns_ptr, &
                                      pos( ix, iy, iz, jx )*Msun_geo, &
                                      pos( ix, iy, iz, jy )*Msun_geo, &
                                      pos( ix, iy, iz, jz )*Msun_geo, &
                                      baryon_density( ix, iy, iz ), &
                                      energy_density( ix, iy, iz ), &
                                      specific_energy( ix, iy, iz ), &
                                      pressure( ix, iy, iz ), &
                                      u_euler( ix, iy, iz, jx ), &
                                      u_euler( ix, iy, iz, jy ), &
                                      u_euler( ix, iy, iz, jz ) )

          ENDDO coords_x
        ENDDO coords_y
      ENDDO coords_z
      !$OMP END PARALLEL DO

      !      ! Print progress on screen
      !      perc= 100*(nx*ny*(iz - 1) &
      !            + nx*(iy - 1) + ix)/( nx*ny*nz )
      !      IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
      !        WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
      !                creturn//" ", perc, "%"
      !      ENDIF
      !
      !    ENDDO coords_x
      !  ENDDO coords_y
      !ENDDO coords_z
      !IF( show_progress ) WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

      PRINT *, "** Subroutine read_id_hydro executed."
      PRINT *

    ENDIF

  END PROCEDURE read_id_hydro


  MODULE PROCEDURE read_id_particles

    !****************************************************
    !
    !# Stores the hydro ID in the arrays needed to
    !  compute the SPH ID
    !
    !  FT 19.11.2020
    !
    !****************************************************

    USE constants, ONLY: amu
    USE utility,   ONLY: determinant_sym3x3, km2m, g2kg

    IMPLICIT NONE

    INTEGER:: a
    DOUBLE PRECISION:: detg

    IF ( C_ASSOCIATED( this% bns_ptr ) ) THEN

      IF( SIZE( x ) /= SIZE( y ) .OR. SIZE( x ) /= SIZE( z ) &
              .OR. SIZE( y ) /= SIZE( z ) )THEN
        PRINT *, "** ERROR: The sizes of the arrays of positions" &
                 // "passed to read_id are not the same."
        PRINT *
        STOP
      ENDIF

      PRINT *, "** Importing ID on particles..."

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( n, this, x, y, z, lapse, &
      !$OMP                     shift_x, shift_y, shift_z, &
      !$OMP                     g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, &
      !$OMP                     baryon_density, energy_density, &
      !$OMP                     specific_energy, pressure, &
      !$OMP                     u_euler_x, u_euler_y, u_euler_z ) &
      !$OMP             PRIVATE( a, detg )
      read_id_loop: DO a= 1, n, 1

        ! The coordinates need to be converted from |sphincs| units (Msun_geo)
        ! to |lorene| units (\(\mathrm{km}\)). See MODULE constants for the
        ! definition of Msun_geo
        CALL get_lorene_id_particles( this% bns_ptr, &
                                      x(a)*Msun_geo, &
                                      y(a)*Msun_geo, &
                                      z(a)*Msun_geo, &
                                      lapse(a), &
                                      shift_x(a), shift_y(a), shift_z(a), &
                                      g_xx(a), &
                                      baryon_density(a), &
                                      energy_density(a), &
                                      specific_energy(a), &
                                      pressure(a), &
                                      u_euler_x(a), &
                                      u_euler_y(a), &
                                      u_euler_z(a) )

    !  ENDDO read_id_loop
    !  !$OMP END PARALLEL DO
    !
    !  DO a= 1, n, 1

        !
        !-- The following follows from the assumption of conformal
        !-- flatness in |lorene|
        !
        g_yy(a)= g_xx(a)
        g_zz(a)= g_xx(a)
        g_xy(a)= zero
        g_xz(a)= zero
        g_yz(a)= zero

        !
        !- Set/unset the geodesic gauge
        !
        IF( this% get_one_lapse() )THEN
          lapse(a)= one
        ENDIF
        IF( this% get_zero_shift() )THEN
          shift_x(a)= zero
          shift_y(a)= zero
          shift_z(a)= zero
        ENDIF

        CALL determinant_sym3x3( [g_xx(a),g_xy(a),g_xz(a), &
                                  g_yy(a),g_yz(a),g_zz(a)], detg )

        IF( ABS( detg ) < 1D-10 )THEN
          PRINT *, "The determinant of the spatial metric is " &
                   // "effectively 0 at the particle ", a
          PRINT *, "detg=", detg
          PRINT *
          STOP
        ELSEIF( detg < 0 )THEN
          PRINT *, "The determinant of the spatial metric is " &
                   // "negative at the particle ", a
          PRINT *, "detg=", detg
          PRINT *
          STOP
        ENDIF

        ! Print progress on screen
       ! perc= 100*a/n
       ! IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
       !   WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
       !           creturn//" ", perc, "%"
       ! ENDIF

      ENDDO read_id_loop
      !$OMP END PARALLEL DO
      IF( show_progress ) WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

      ! Convert the baryon density and pressure to units of amu (SPH code units)
      baryon_density= baryon_density*((Msun_geo*km2m)**3)/(amu*g2kg)
      pressure      = pressure*((Msun_geo*km2m)**3)/(amu*g2kg)

      PRINT *, "** Subroutine read_id_particles executed."
      PRINT *

    ENDIF

  END PROCEDURE read_id_particles


  MODULE PROCEDURE read_id_mass_b

    !****************************************************
    !
    !# Stores the hydro ID in the arrays needed to
    !  compute the baryon mass, storing it to variables
    !  (not arrays as the others SUBROUTINES in
    !  this SUBMODULE).
    !
    !  FT 15.04.2021
    !
    !****************************************************

    USE utility,  ONLY: density_si2cu
    USE tensor,   ONLY: jxx, jxy, jxz, jyy, jyz, jzz

    IMPLICIT NONE

    IF ( C_ASSOCIATED( this% bns_ptr ) ) THEN

      ! The coordinates need to be converted from |sphincs| units (Msun_geo)
      ! to |lorene| units (\(\mathrm{km}\)).
      ! See MODULE constants for the definition of Msun_geo
      CALL get_lorene_id_mass_b( this% bns_ptr, &
                                 x*Msun_geo, &
                                 y*Msun_geo, &
                                 z*Msun_geo, &
                                 g(jxx), &
                                 baryon_density, &
                                 gamma_euler )

      g(jxy)= zero
      g(jxz)= zero
      g(jyy)= g(jxx)
      g(jyz)= zero
      g(jzz)= g(jxx)

      baryon_density= baryon_density*density_si2cu

    ENDIF

  END PROCEDURE read_id_mass_b


  MODULE PROCEDURE read_id_k

    !****************************************************
    !
    !# Stores the components of the extrinsic curvature
    !  in arrays
    !
    !  FT 25.11.2020
    !
    !****************************************************

    IMPLICIT NONE

    IF ( C_ASSOCIATED( this% bns_ptr ) ) THEN

      IF( SIZE( x ) /= SIZE( y ) .OR. SIZE( x ) /= SIZE( z ) &
              .OR. SIZE( y ) /= SIZE( z ) )THEN
        PRINT *, "** ERROR: The sizes of the arrays of positions" &
                 // "passed to read_id are not the same."
        PRINT *
        STOP
      ENDIF

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( n, this, x, y, z, &
      !$OMP                     k_xx, k_xy, k_xz, k_yy, k_yz, k_zz ) &
      !$OMP             PRIVATE( itr )
      read_id_loop: DO itr= 1, n, 1

        ! The coordinates need to be converted from |sphincs| units (Msun_geo)
        ! to |lorene| units (\(\mathrm{km}\)). See MODULE constants for the definition of
        ! Msun_geo
        CALL get_lorene_id_k( this% bns_ptr, &
                              x( itr )*Msun_geo, &
                              y( itr )*Msun_geo, &
                              z( itr )*Msun_geo, &
                              k_xx( itr ), &
                              k_xy( itr ), &
                              k_xz( itr ), &
                              k_yy( itr ), &
                              k_yz( itr ), &
                              k_zz( itr ) )

      ENDDO read_id_loop
      !$OMP END PARALLEL DO

      DO itr= 1, n, 1

        !
        !-- Convert the extrinsic curvature from |lorene| units to
        !-- |sphincs| units
        !
        k_xx( itr )= k_xx( itr )*Msun_geo
        k_xy( itr )= k_xy( itr )*Msun_geo
        k_xz( itr )= k_xz( itr )*Msun_geo
        k_yy( itr )= k_yy( itr )*Msun_geo
        k_yz( itr )= k_yz( itr )*Msun_geo
        k_zz( itr )= k_zz( itr )*Msun_geo

        ! Print progress on screen
        perc= 100*itr/n
        IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
          WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
                  creturn//" ", perc, "%"
        ENDIF

      ENDDO
      IF( show_progress ) WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

      PRINT *, "** Subroutine read_id_k executed."
      PRINT *

    ENDIF

  END PROCEDURE read_id_k


  !-----------------!
  !--  FUNCTIONS  --!
  !-----------------!


  MODULE PROCEDURE read_bnslorene_mass_density

    !***********************************************
    !
    !# Returns the |lorene| mass density at the point
    !  given as argument, in units of
    !  \(M_\odot/L_\odot^3\).
    !
    !  FT
    !
    !***********************************************

    USE, INTRINSIC:: ISO_C_BINDING, ONLY: C_ASSOCIATED
    USE utility,                    ONLY: density_si2cu

    IMPLICIT NONE

    IF ( C_ASSOCIATED( this% bns_ptr ) )THEN

      ! The coordinates need to be converted from |sphincs| units (Msun_geo)
      ! to |lorene| units (\(\mathrm{km}\)). See MODULE constants for the
      ! definition of Msun_geo
      res= get_lorene_mass_density( this% bns_ptr, &
                                    x*Msun_geo, &
                                    y*Msun_geo, &
                                    z*Msun_geo )*density_si2cu

    ENDIF

  END PROCEDURE read_bnslorene_mass_density


  MODULE PROCEDURE read_bnslorene_pressure

    !***********************************************
    !
    !# Returns the |lorene| pressure at the point
    !  given as argument, in units of
    !  \(M_\odot c^2/L_\odot^3\).
    !
    !  FT 11.02.2022
    !
    !***********************************************

    USE, INTRINSIC:: ISO_C_BINDING, ONLY: C_ASSOCIATED
    USE utility,                    ONLY: density_si2cu

    IMPLICIT NONE

    IF ( C_ASSOCIATED( this% bns_ptr ) )THEN

      ! The coordinates need to be converted from |sphincs| units (Msun_geo)
      ! to |lorene| units (\(\mathrm{km}\)). See MODULE constants for the
      ! definition of Msun_geo
      res= get_lorene_pressure( this% bns_ptr, &
                                x*Msun_geo, &
                                y*Msun_geo, &
                                z*Msun_geo )*density_si2cu

    ENDIF

  END PROCEDURE read_bnslorene_pressure


  MODULE PROCEDURE read_spatial_metric

    !***********************************************
    !
    !# Returns the |lorene| conformal factor to the
    !  4th power, equal to the diagonal components
    !  of the conformally flat spatial ADM metric.
    !
    !  FT 15.04.2021
    !
    !***********************************************

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_ASSOCIATED

    IMPLICIT NONE

    IF ( C_ASSOCIATED( this% bns_ptr ) )THEN

      ! The coordinates need to be converted from |sphincs| units (Msun_geo)
      ! to |lorene| units (\(\mathrm{km}\)). See MODULE constants for the definition of
      ! Msun_geo
      res= get_lorene_spatial_metric( this% bns_ptr, &
                                      x*Msun_geo, &
                                      y*Msun_geo, &
                                      z*Msun_geo )

    ENDIF

  END PROCEDURE read_spatial_metric


  MODULE PROCEDURE is_hydro_positive

    !************************************************
    !
    !# Return .FALSE. if the energy density is nonpositive
    !  or if the specific energy is nonpositive,
    !  or if the pressure is nonpositive
    !  at the specified point; .TRUE. otherwise
    !
    !  FT 12.03.2021
    !
    !************************************************

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_ASSOCIATED

    IMPLICIT NONE

    INTEGER:: tmp

    IF ( C_ASSOCIATED( this% bns_ptr ) )THEN

      ! The coordinates need to be converted from |sphincs| units (Msun_geo)
      ! to |lorene| units (\(\mathrm{km}\)). See MODULE constants for the
      ! definition  of Msun_geo
      tmp= positive_hydro( this% bns_ptr, x*Msun_geo, &
                                          y*Msun_geo, &
                                          z*Msun_geo )

      IF( tmp == 1 )THEN
        res= .TRUE.
      ELSE
        res= .FALSE.
      ENDIF

    ENDIF

  END PROCEDURE is_hydro_positive


END SUBMODULE read
