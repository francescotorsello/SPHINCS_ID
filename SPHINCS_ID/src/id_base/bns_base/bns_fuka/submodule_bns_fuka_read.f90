! File:         submodule_bnsfuka_read.f90
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

SUBMODULE (bns_fuka) read

  !****************************************************
  !
  !# Implementation of the methods of TYPE bnsfuka that
  !  read |bns| data using |fuka|
  !
  !  FT 09.02.2022
  !
  !****************************************************


  USE utility, ONLY: zero, one, two, Msun_geo, is_finite_number


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE read_fuka_id_member

    !**************************************************
    !
    !# Stores the |id| in the [[bnsfuka]] member arrays
    !
    !  FT 09.02.2022
    !
    !**************************************************

    IMPLICIT NONE

    !IF ( C_ASSOCIATED( this% bns_ptr ) ) THEN
    !
    !  IF( SIZE( x ) /= SIZE( y ) .OR. SIZE( x ) /= SIZE( z ) &
    !        .OR. SIZE( y ) /= SIZE( z ) )THEN
    !    PRINT *, "** ERROR: The sizes of the arrays of positions" &
    !             // "passed to read_lorene_id are not the same."
    !    PRINT *
    !    STOP
    !  ENDIF
    !
    !  IF( ALLOCATED( this% lapse )   .AND. &
    !      ALLOCATED( this% shift_x ) .AND. &
    !      ALLOCATED( this% shift_y ) .AND. &
    !      ALLOCATED( this% shift_z ) .AND. &
    !      ALLOCATED( this% g_xx ) .AND. ALLOCATED( this% g_xy ) .AND. &
    !      ALLOCATED( this% g_xz ) .AND. ALLOCATED( this% g_yy ) .AND. &
    !      ALLOCATED( this% g_yz ) .AND. ALLOCATED( this% g_zz ) .AND. &
    !      ALLOCATED( this% k_xx ) .AND. ALLOCATED( this% k_xy ) .AND. &
    !      ALLOCATED( this% k_xz ) .AND. ALLOCATED( this% k_yy ) .AND. &
    !      ALLOCATED( this% k_yz ) .AND. ALLOCATED( this% k_zz ) .AND. &
    !      ALLOCATED( this% baryon_density )  .AND. &
    !      ALLOCATED( this% energy_density )  .AND. &
    !      ALLOCATED( this% specific_energy ) .AND. &
    !      ALLOCATED( this% v_euler_x )       .AND. &
    !      ALLOCATED( this% v_euler_y )       .AND. &
    !      ALLOCATED( this% v_euler_z ) &
    !  )THEN
    !
    !    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !    !$OMP             SHARED( n, this, x, y, z ) &
    !    !$OMP             PRIVATE(itr)
    !    read_fuka_id_loop: DO itr= 1, n, 1
    !
    !      ! The coordinates need to be converted from |sphincs| units (Msun_geo)
    !      ! to |fuka| units (\(\mathrm{km}\)). See MODULE constants for the
    !      ! definition of Msun_geo
    !      CALL get_lorene_id( this% bns_ptr, &
    !                          x(itr), &
    !                          y(itr), &
    !                          z(itr), &
    !                          this% lapse(itr), &
    !                          this% shift_x(itr), &
    !                          this% shift_y(itr), &
    !                          this% shift_z(itr), &
    !                          this% g_xx(itr), &
    !                          this% k_xx(itr), &
    !                          this% k_xy(itr), &
    !                          this% k_xz(itr), &
    !                          this% k_yy(itr), &
    !                          this% k_yz(itr), &
    !                          this% k_zz(itr), &
    !                          this% baryon_density(itr), &
    !                          this% energy_density(itr), &
    !                          this% specific_energy(itr), &
    !                          this% v_euler_x(itr), &
    !                          this% v_euler_y(itr), &
    !                          this% v_euler_z(itr) )
    !
    !    ENDDO read_fuka_id_loop
    !    !$OMP END PARALLEL DO
    !
    !    DO itr= 1, n, 1
    !
    !      !
    !      !-- The following follows from the assumption of conformal
    !      !-- flatness in |fuka|
    !      !
    !      this% g_yy(itr)= this% g_xx(itr)
    !      this% g_zz(itr)= this% g_xx(itr)
    !      this% g_xy(itr)= 0.0D0
    !      this% g_xz(itr)= 0.0D0
    !      this% g_yz(itr)= 0.0D0
    !
    !      !
    !      !- Set/unset the geodesic gauge
    !      !
    !      IF( this% get_one_lapse() )THEN
    !        this% lapse(itr)= 1.0D0
    !      ENDIF
    !      IF( this% get_zero_shift() )THEN
    !        this% shift_x(itr)= 0.0D0
    !        this% shift_y(itr)= 0.0D0
    !        this% shift_z(itr)= 0.0D0
    !      ENDIF
    !
    !      !
    !      !-- Convert the extrinsic curvature from |fuka| units to
    !      !-- |sphincs| units
    !      !
    !      this% k_xx(itr)= this% k_xx(itr)
    !      this% k_xy(itr)= this% k_xy(itr)
    !      this% k_xz(itr)= this% k_xz(itr)
    !      this% k_yy(itr)= this% k_yy(itr)
    !      this% k_yz(itr)= this% k_yz(itr)
    !      this% k_zz(itr)= this% k_zz(itr)
    !
    !      ! Print progress on screen
    !      perc= 100*itr/n
    !      IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
    !        WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
    !                creturn//" ", perc, "%"
    !      ENDIF
    !
    !    ENDDO
    !    IF( show_progress ) WRITE( *, "(A1)", ADVANCE= "NO" ) creturn
    !
    !  ELSE
    !
    !    PRINT *, "** ERROR: Memory was not allocated before calling " &
    !             // "read_fuka_id in read_lorene_id (TYPE particles)."
    !    PRINT *
    !    STOP
    !
    !  ENDIF
    !
    !  PRINT *, "** Subroutine read_lorene_id executed."
    !  PRINT *
    !
    !ENDIF

  END PROCEDURE read_fuka_id_member


  MODULE PROCEDURE read_fuka_id_full

    !**************************************************
    !
    !# Stores the |id| in non-[[bnsfuka]]-member arrays
    !  with the same shape as the [[bnsfuka]] member arrays
    !
    !  Created:     FT 09.02.2022
    !  Last update: FT 27.05.2022
    !
    !**************************************************

    IMPLICIT NONE

    IF ( C_ASSOCIATED( this% bns_ptr ) ) THEN

      IF( SIZE( x ) /= SIZE( y ) .OR. SIZE( x ) /= SIZE( z ) &
            .OR. SIZE( y ) /= SIZE( z ) )THEN
        PRINT *, "** ERROR: The sizes of the arrays of positions" &
                 // "passed to read_fuka_id are not the same."
        PRINT *
        STOP
      ENDIF

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( n, this, x, y, z, lapse, &
      !$OMP                     shift_x, shift_y, shift_z, &
      !$OMP                     g_xx, k_xx, k_xy, k_xz, k_yy, k_yz, k_zz, &
      !$OMP                     baryon_density, energy_density, &
      !$OMP                     specific_energy, &
      !$OMP                     u_euler_x, u_euler_y, u_euler_z ) &
      !$OMP             PRIVATE(itr)
      read_fuka_id_loop: DO itr= 1, n, 1

        ! The coordinates need to be converted from |sphincs| units (Msun_geo)
        ! to |fuka| units (\(\mathrm{km}\)). See MODULE constants for the
        ! definition of Msun_geo
        CALL get_fuka_id( this% bns_ptr, &
                          x(itr), &
                          y(itr), &
                          z(itr), &
                          lapse(itr), &
                          shift_x(itr), &
                          shift_y(itr), &
                          shift_z(itr), &
                          g_xx(itr), &
                          k_xx(itr), &
                          k_xy(itr), &
                          k_xz(itr), &
                          k_yy(itr), &
                          k_yz(itr), &
                          k_zz(itr), &
                          baryon_density(itr), &
                          energy_density(itr), &
                          specific_energy(itr), &
                          u_euler_x(itr), &
                          u_euler_y(itr), &
                          u_euler_z(itr) )

      ENDDO read_fuka_id_loop
      !$OMP END PARALLEL DO

      DO itr= 1, n, 1

        !
        !-- The following follows from the assumption of conformal
        !-- flatness in |fuka|
        !
        g_yy(itr)= g_xx(itr)
        g_zz(itr)= g_xx(itr)
        g_xy(itr)= zero
        g_xz(itr)= zero
        g_yz(itr)= zero

        !
        !- Set/unset the geodesic gauge
        !
        IF( this% get_one_lapse() )THEN
          lapse(itr)= one
        ENDIF
        IF( this% get_zero_shift() )THEN
          shift_x(itr)= zero
          shift_y(itr)= zero
          shift_z(itr)= zero
        ENDIF

        !
        !-- Convert the extrinsic curvature from |fuka| units to
        !-- |sphincs| units
        !
        k_xx(itr)= k_xx(itr)
        k_xy(itr)= k_xy(itr)
        k_xz(itr)= k_xz(itr)
        k_yy(itr)= k_yy(itr)
        k_yz(itr)= k_yz(itr)
        k_zz(itr)= k_zz(itr)

        ! Print progress on screen
       ! perc= 100*itr/n
       ! IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
       !   WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
       !           creturn//" ", perc, "%"
       ! ENDIF

      ENDDO
     ! IF( show_progress ) WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

      PRINT *, "** Subroutine read_fuka_id executed."
      PRINT *

    ENDIF

  END PROCEDURE read_fuka_id_full


  MODULE PROCEDURE read_fuka_id_spacetime

    !*******************************************************
    !
    !# Stores the spacetime |id| in multi-dimensional arrays
    !  needed to compute the BSSN variables and constraints
    !
    !  Created:     FT 09.02.2022
    !  Last update: FT 27.05.2022
    !
    !*******************************************************

    USE tensor,   ONLY: jxx, jxy, jxz, &
                        jyy, jyz, jzz, jx, jy, jz, n_sym4x4
    USE utility,  ONLY: determinant_sym3x3, determinant_sym4x4

    IMPLICIT NONE

    INTEGER:: i, j, k

    DOUBLE PRECISION:: detg
    DOUBLE PRECISION:: detg4
    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: g4

    ! g4 is allocatable to allocate it on the heap
    ! Allocating it on the stack might exceed stack memory,
    ! causing a segmentation fault
    ALLOCATE( g4( nx, ny, nz, n_sym4x4 ) )

    IF ( C_ASSOCIATED( this% bns_ptr ) ) THEN

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( nx, ny, nz, this, pos, &
      !$OMP                     lapse, shift, g, ek, g4 ) &
      !$OMP             PRIVATE( i, j, k, detg, detg4 )
      coords_z: DO k= 1, nz, 1
        coords_y: DO j= 1, ny, 1
          coords_x: DO i= 1, nx, 1

            CALL get_fuka_id_spacetime( this% bns_ptr, &
                                        pos( i, j, k, jx ), &
                                        pos( i, j, k, jy ), &
                                        pos( i, j, k, jz ), &
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

            !
            !-- The following follows from the assumption of
            !-- conformal flatness in |fuka|
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

            CALL determinant_sym3x3( g(i,j,k,:), detg )

            IF( ABS( detg ) < 1D-10 )THEN
              PRINT *, "** ERROR! The determinant of the spatial metric " &
                       // "is effectively 0 at the grid point " &
                       // "(i,j,k)= (", i, ",", j,",", k, ")."
              PRINT *, " * detg=", detg
              PRINT *
              STOP
            ELSEIF( detg < zero )THEN
              PRINT *, "** ERROR! The determinant of the spatial metric " &
                       // "is negative at the grid point " &
                       // "(i,j,k)= (", i, ",", j,",", k, ")."
              PRINT *, " * detg=", detg
              PRINT *
              STOP
            ELSEIF( .NOT.is_finite_number(detg) )THEN
              PRINT *, "** ERROR! The determinant of the spatial metric "&
                       // "is not a finite number at " &
                       // "(i,j,k)= (", i, ",", j,",", k, ")."
              PRINT *, " * detg=", detg
              PRINT *
              STOP
            ENDIF

            CALL compute_g4( lapse(i,j,k), shift(i,j,k,:), &
                             g(i,j,k,:), g4(i,j,k,:) )

            CALL determinant_sym4x4( g4(i,j,k,:), detg4 )

            IF( ABS( detg4 ) < 1D-10 )THEN
              PRINT *, "** ERROR! The determinant of the spacetime metric "&
                       // "is effectively 0 at the grid point " &
                       // "(i,j,k)= (", i, ",", j,",", k, ")."
              PRINT *, " * detg4=", detg4
              PRINT *
              STOP
            ELSEIF( detg4 > 0 )THEN
              PRINT *, "** ERROR! The determinant of the spacetime metric "&
                       // "is positive at the grid point " &
                       // "(i,j,k)= (", i, ",", j,",", k, ")."
              PRINT *, " * detg4=", detg4
              PRINT *
              STOP
            ELSEIF( .NOT.is_finite_number(detg4) )THEN
              PRINT *, "** ERROR! The determinant of the spacetime metric "&
                       // "is not a finite number at " &
                       // "(i,j,k)= (", i, ",", j,",", k, ")."
              PRINT *, " * detg4=", detg4
              PRINT *
              STOP
            ENDIF

          ENDDO coords_x
        ENDDO coords_y
      ENDDO coords_z
      !$OMP END PARALLEL DO

      PRINT *, "** Subroutine read_lorene_id executed."
      PRINT *

    ENDIF

  END PROCEDURE read_fuka_id_spacetime


  MODULE PROCEDURE read_fuka_id_hydro

    !*******************************************************
    !
    !# Stores the hydro |id| in the arrays needed to compute
    !  the constraints on the refined mesh
    !
    !  Created:     FT 09.02.2022
    !  Last update: FT 27.05.2022
    !
    !*******************************************************

    USE tensor,     ONLY: jx, jy, jz

    IMPLICIT NONE

    INTEGER:: i, j, k

    IF ( C_ASSOCIATED( this% bns_ptr ) ) THEN

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( nx, ny, nz, this, pos, &
      !$OMP                     baryon_density, energy_density, &
      !$OMP                     pressure, u_euler ) &
      !$OMP             PRIVATE( i, j, k )
      coords_z: DO i= 1, nz, 1
        coords_y: DO j= 1, ny, 1
          coords_x: DO k= 1, nx, 1

            CALL get_fuka_id_hydro( this% bns_ptr, &
                                    pos( i, j, k, jx ), &
                                    pos( i, j, k, jy ), &
                                    pos( i, j, k, jz ), &
                                    baryon_density( i, j, k ), &
                                    energy_density( i, j, k ), &
                                    pressure( i, j, k ), &
                                    u_euler( i, j, k, jx ), &
                                    u_euler( i, j, k, jy ), &
                                    u_euler( i, j, k, jz ) )

          ENDDO coords_x
        ENDDO coords_y
      ENDDO coords_z
      !$OMP END PARALLEL DO

      PRINT *, "** Subroutine read_fuka_id_hydro executed."
      PRINT *

    ENDIF

  END PROCEDURE read_fuka_id_hydro


  MODULE PROCEDURE read_fuka_id_particles

    !****************************************************
    !
    !# Stores the hydro |id| in the arrays needed to
    !  compute the |sph| |id|
    !
    !  Created:     FT 09.02.2022
    !  Last update: FT 27.05.2022
    !
    !****************************************************

    USE constants,  ONLY: Msun, amu
    USE utility,    ONLY: km2m, g2kg, determinant_sym3x3

    IMPLICIT NONE

    INTEGER:: a
    DOUBLE PRECISION:: detg

    IF ( C_ASSOCIATED( this% bns_ptr ) ) THEN

      IF( SIZE( x ) /= SIZE( y ) .OR. SIZE( x ) /= SIZE( z ) &
              .OR. SIZE( y ) /= SIZE( z ) )THEN
        PRINT *, "** ERROR: The sizes of the arrays of positions" &
                 // "passed to read_lorene_id are not the same."
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
      read_fuka_id_loop: DO a= 1, n, 1

        CALL get_fuka_id_particles( this% bns_ptr, &
                                    x(a), y(a), z(a), &
                                    lapse(a), &
                                    shift_x(a), shift_y(a), shift_z(a), &
                                    g_xx(a), &
                                    baryon_density(a), &
                                    energy_density(a), &
                                    pressure(a), &
                                    u_euler_x(a), &
                                    u_euler_y(a), &
                                    u_euler_z(a) )

        !
        !-- The following follows from the assumption of conformal
        !-- flatness in |fuka|
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
          PRINT *, "** ERROR! The determinant of the spatial metric " &
                   // "is effectively 0 at particle ", a, "."
          PRINT *, " * detg=", detg
          PRINT *
          STOP
        ELSEIF( detg < zero )THEN
          PRINT *, "** ERROR! The determinant of the spatial metric " &
                   // "is negative at particle ", a, "."
          PRINT *, " * detg=", detg
          PRINT *
          STOP
        ELSEIF( .NOT.is_finite_number(detg) )THEN
          PRINT *, "** ERROR! The determinant of the spatial metric "&
                   // "is not a finite number at particle ", a, "."
          PRINT *
          STOP
        ENDIF

      ENDDO read_fuka_id_loop
      !$OMP END PARALLEL DO

      ! Convert the baryon density and pressure to units of amu
      ! (|sph| code units)
      baryon_density= baryon_density*Msun/amu
      pressure      = pressure*MSun/amu

      PRINT *, "** Subroutine read_fuka_id_particles executed."
      PRINT *

    ENDIF

  END PROCEDURE read_fuka_id_particles


  MODULE PROCEDURE read_fuka_id_mass_b

    !****************************************************
    !
    !# Stores the hydro |id| in the arrays needed to
    !  compute the baryon mass, storing it to variables
    !  (not arrays as the others SUBROUTINES in
    !  the [[bns_read]] SUBMODULE).
    !
    !  Created:     FT 09.02.2022
    !  Last update: FT 27.05.2022
    !
    !****************************************************

    USE tensor,   ONLY: jxx, jxy, jxz, jyy, jyz, jzz

    IMPLICIT NONE

    IF ( C_ASSOCIATED( this% bns_ptr ) ) THEN

      ! The coordinates need to be converted from |sphincs| units (Msun_geo)
      ! to |fuka| units (\(\mathrm{km}\)).
      ! See MODULE constants for the definition of Msun_geo
      CALL get_fuka_id_mass_b( this% bns_ptr, x, y, z, &
                               g(jxx), &
                               baryon_density, &
                               gamma_euler )

      g(jxy)= zero
      g(jxz)= zero
      g(jyy)= g(jxx)
      g(jyz)= zero
      g(jzz)= g(jxx)

    ENDIF

  END PROCEDURE read_fuka_id_mass_b


  MODULE PROCEDURE read_fuka_id_k

    !****************************************************
    !
    !# Stores the components of the extrinsic curvature
    !  in arrays
    !
    !  Created:     FT 09.02.2022
    !  Last update: FT 27.05.2022
    !
    !  @warning deprecated?
    !
    !****************************************************

    IMPLICIT NONE

   ! INTEGER:: a

   ! IF ( C_ASSOCIATED( this% bns_ptr ) ) THEN
   !
   !   IF( SIZE( x ) /= SIZE( y ) .OR. SIZE( x ) /= SIZE( z ) &
   !           .OR. SIZE( y ) /= SIZE( z ) )THEN
   !     PRINT *, "** ERROR: The sizes of the arrays of positions" &
   !              // "passed to read_lorene_id are not the same."
   !     PRINT *
   !     STOP
   !   ENDIF
   !
   !   !$OMP PARALLEL DO DEFAULT( NONE ) &
   !   !$OMP             SHARED( n, this, x, y, z, &
   !   !$OMP                     k_xx, k_xy, k_xz, k_yy, k_yz, k_zz ) &
   !   !$OMP             PRIVATE(a)
   !   read_fuka_id_loop: DO a= 1, n, 1
   !
   !     ! The coordinates need to be converted from |sphincs| units (Msun_geo)
   !     ! to |fuka| units (\(\mathrm{km}\)). See MODULE constants for the
   !     ! definition of Msun_geo
   !     CALL get_lorene_id_k( this% bns_ptr, &
   !                           x(a), &
   !                           y(a), &
   !                           z(a), &
   !                           k_xx(a), &
   !                           k_xy(a), &
   !                           k_xz(a), &
   !                           k_yy(a), &
   !                           k_yz(a), &
   !                           k_zz(a) )
   !
   !   ENDDO read_fuka_id_loop
   !   !$OMP END PARALLEL DO
   !
   !   DO a= 1, n, 1
   !
   !     !
   !     !-- Convert the extrinsic curvature from |fuka| units to
   !     !-- |sphincs| units
   !     !
   !     k_xx(a)= k_xx(a)
   !     k_xy(a)= k_xy(a)
   !     k_xz(a)= k_xz(a)
   !     k_yy(a)= k_yy(a)
   !     k_yz(a)= k_yz(a)
   !     k_zz(a)= k_zz(a)
   !
   !     ! Print progress on screen
   !     perc= 100*a/n
   !     IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
   !       WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
   !               creturn//" ", perc, "%"
   !     ENDIF
   !
   !   ENDDO
   !   IF( show_progress ) WRITE( *, "(A1)", ADVANCE= "NO" ) creturn
   !
   !   PRINT *, "** Subroutine read_lorene_id_k executed."
   !   PRINT *
   !
   ! ENDIF

  END PROCEDURE read_fuka_id_k


  !-----------------!
  !--  FUNCTIONS  --!
  !-----------------!


  MODULE PROCEDURE read_fuka_mass_density

    !***********************************************
    !
    !# Returns the |fuka| mass density at the point
    !  given as argument, in units of
    !  \(M_\odot/\ell_\odot^3\).
    !
    !  Created:     FT 09.02.2022
    !  Last update: FT 27.05.2022
    !
    !***********************************************

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_ASSOCIATED
    USE utility,                     ONLY: lorene2hydrobase

    IMPLICIT NONE

    IF ( C_ASSOCIATED( this% bns_ptr ) )THEN

      ! The coordinates need to be converted from |sphincs| units (Msun_geo)
      ! to |fuka| units (\(\mathrm{km}\)). See MODULE constants for the
      ! definition of Msun_geo
      res= get_fuka_mass_density( this% bns_ptr, x, y, z )

    ENDIF

  END PROCEDURE read_fuka_mass_density


  MODULE PROCEDURE read_fuka_pressure

    !***********************************************
    !
    !# Returns the |fuka| pressure at the point
    !  given as argument, in units of
    !  \([\mathrm{kg}\,c^2\, \mathrm{m}^{-3}]\).
    !
    !  Created:     FT 27.05.2022
    !  Last update: FT 27.05.2022
    !
    !***********************************************

    USE, INTRINSIC:: ISO_C_BINDING, ONLY: C_ASSOCIATED

    IMPLICIT NONE

    IF ( C_ASSOCIATED( this% bns_ptr ) )THEN

      ! The coordinates need to be converted from |sphincs| units (Msun_geo)
      ! to |fuka| units (\(\mathrm{km}\)). See MODULE constants for the
      ! definition of Msun_geo
      res= get_fuka_pressure( this% bns_ptr, x, y, z )

    ENDIF

  END PROCEDURE read_fuka_pressure


  MODULE PROCEDURE read_fuka_spatial_metric

    !***********************************************
    !
    !# Returns the |fuka| conformal factor to the
    !  4th power, equal to the diagonal components
    !  of the conformally flat spatial ADM metric.
    !
    !  Created:     FT 09.02.2022
    !  Last update: FT 27.05.2022
    !
    !***********************************************

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_ASSOCIATED

    IMPLICIT NONE

    IF ( C_ASSOCIATED( this% bns_ptr ) )THEN

      ! The coordinates need to be converted from |sphincs| units (Msun_geo)
      ! to |fuka| units (\(\mathrm{km}\)). See MODULE constants for the
      ! definition of Msun_geo
      res= get_fuka_spatial_metric( this% bns_ptr, x, y, z )

    ENDIF

  END PROCEDURE read_fuka_spatial_metric


  MODULE PROCEDURE is_hydro_positive

    !************************************************
    !
    !# Return 1 if the energy density is nonpositive
    !  or if the specific energy is nonpositive,
    !  or if the pressure is nonpositive
    !  at the specified point
    !
    !  Created:     FT 09.02.2022
    !  Last update: FT 27.05.2022
    !
    !************************************************

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_ASSOCIATED

    IMPLICIT NONE

    INTEGER:: tmp

    IF ( C_ASSOCIATED( this% bns_ptr ) )THEN

      ! The coordinates need to be converted from |sphincs| units (Msun_geo)
      ! to |fuka| units (\(\mathrm{km}\)). See MODULE constants for the
      ! definition of Msun_geo
      tmp= positive_hydro( this% bns_ptr, x, y, z )

      IF( tmp == 1 )THEN
        res= .TRUE.
      ELSE
        res= .FALSE.
      ENDIF

    ENDIF

  END PROCEDURE is_hydro_positive


END SUBMODULE read
