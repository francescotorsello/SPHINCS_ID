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
      !$OMP                     specific_energy, pressure, &
      !$OMP                     u_euler_x, u_euler_y, u_euler_z, &
      !$OMP                     g_xy, g_xz, g_yy, g_yz, g_zz ) &
      !$OMP             PRIVATE( itr )
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
                          specific_energy(itr), &
                          pressure(itr), &
                          u_euler_x(itr), &
                          u_euler_y(itr), &
                          u_euler_z(itr) )

        energy_density(itr)= baryon_density(itr)*(one + specific_energy(itr))

    !  ENDDO read_fuka_id_loop
    !  !$OMP END PARALLEL DO
    !
    !  DO itr= 1, n, 1

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
       ! k_xx(itr)= k_xx(itr)
       ! k_xy(itr)= k_xy(itr)
       ! k_xz(itr)= k_xz(itr)
       ! k_yy(itr)= k_yy(itr)
       ! k_yz(itr)= k_yz(itr)
       ! k_zz(itr)= k_zz(itr)

        ! Print progress on screen
       ! perc= 100*itr/n
       ! IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
       !   WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
       !           creturn//" ", perc, "%"
       ! ENDIF

      ENDDO read_fuka_id_loop
      !$OMP END PARALLEL DO
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

    INTEGER, PARAMETER:: mpi_ranks= 40
    LOGICAL, PARAMETER:: debug= .FALSE.

    INTEGER:: i, j, k

    DOUBLE PRECISION:: xmin, xmax, ymin, ymax, zmin, zmax

    DOUBLE PRECISION:: detg
    DOUBLE PRECISION:: detg4
    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: g4
    DOUBLE PRECISION, DIMENSION(nx, ny, nz, 9):: id_tmp

    LOGICAL:: exist

    CHARACTER(LEN=:), ALLOCATABLE:: filename

    ! TODO: parallelize these ones
    xmin= MINVAL( pos(:,1,1,jx), DIM= 1 )
    xmax= MAXVAL( pos(:,1,1,jx), DIM= 1 )
    ymin= MINVAL( pos(1,:,1,jy), DIM= 1 )
    ymax= MAXVAL( pos(1,:,1,jy), DIM= 1 )
    zmin= MINVAL( pos(1,1,:,jz), DIM= 1 )
    zmax= MAXVAL( pos(1,1,:,jz), DIM= 1 )

    CALL this% run_kadath_reader( mpi_ranks, nx, ny, nz, &
                                  xmin, xmax, ymin, ymax, zmin, zmax, &
                                  id_tmp(:,:,:,id$x:id$z), &
                                  lapse, &
                                  shift(:,:,:,jx), &
                                  shift(:,:,:,jy), &
                                  shift(:,:,:,jz), &
                                  g(:,:,:,jxx), &
                                  g(:,:,:,jyy), &
                                  g(:,:,:,jzz), &
                                  g(:,:,:,jxy), &
                                  g(:,:,:,jxz), &
                                  g(:,:,:,jyz), &
                                  ek(:,:,:,jxx), &
                                  ek(:,:,:,jyy), &
                                  ek(:,:,:,jzz), &
                                  ek(:,:,:,jxy), &
                                  ek(:,:,:,jxz), &
                                  ek(:,:,:,jyz), &
                                  id_tmp(:,:,:,4), &
                                  id_tmp(:,:,:,5), &
                                  id_tmp(:,:,:,6), &
                                  id_tmp(:,:,:,7), &
                                  id_tmp(:,:,:,8), &
                                  id_tmp(:,:,:,9), &
                                  this% filename )

    ALLOCATE( g4( nx, ny, nz, n_sym4x4 ) )

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( nx, ny, nz, this, pos, &
    !$OMP                     lapse, shift, g, ek, g4, id_tmp ) &
    !$OMP             PRIVATE( i, j, k, detg, detg4 )
    DO k= 1, nz, 1
      DO j= 1, ny, 1
        DO i= 1, nx, 1

          !lapse          = id_tmp(i,j,k,id$lapse)
          !
          !shift(i,j,k,jx)= id_tmp(i,j,k,id$shiftx)
          !shift(i,j,k,jy)= id_tmp(i,j,k,id$shifty)
          !shift(i,j,k,jz)= id_tmp(i,j,k,id$shiftz)
          !
          !g(i,j,k,jxx)   = id_tmp(i,j,k,id$gxx)
          !
          !ek(i,j,k,jxx)  = id_tmp(i,j,k,id$kxx)
          !ek(i,j,k,jxy)  = id_tmp(i,j,k,id$kxy)
          !ek(i,j,k,jxz)  = id_tmp(i,j,k,id$kxz)
          !ek(i,j,k,jyy)  = id_tmp(i,j,k,id$kyy)
          !ek(i,j,k,jyz)  = id_tmp(i,j,k,id$kyz)
          !ek(i,j,k,jzz)  = id_tmp(i,j,k,id$kzz)
          !
          !!
          !!-- The following follows from the assumption of
          !!-- conformal flatness in |fuka|
          !!
          !g( i, j, k, jyy )= g( i, j, k, jxx )
          !g( i, j, k, jzz )= g( i, j, k, jxx )
          !g( i, j, k, jxy )= zero
          !g( i, j, k, jxz )= zero
          !g( i, j, k, jyz )= zero

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

        ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    PRINT *, "** Subroutine read_fuka_id_spacetime executed."
    PRINT *

    IF( debug )THEN

      filename= "dbg-spaid.dat"

      INQUIRE( FILE= TRIM(filename), EXIST= exist )

      IF( exist )THEN
        OPEN( UNIT= 2, FILE= TRIM(filename), STATUS= "REPLACE", &
              FORM= "FORMATTED", &
              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
              IOMSG= err_msg )
      ELSE
        OPEN( UNIT= 2, FILE= TRIM(filename), STATUS= "NEW", &
              FORM= "FORMATTED", &
              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
      ENDIF
      IF( ios > 0 )THEN
      PRINT *, "...error when opening " // TRIM(filename), &
               ". The error message is", err_msg
      STOP
      ENDIF

      DO k= 1, nz, 1
        DO j= 1, ny, 1
          DO i= 1, nx, 1

            WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
              pos(i,j,k,jx), id_tmp(i,j,k,id$x), &
              pos(i,j,k,jy), id_tmp(i,j,k,id$y), &
              pos(i,j,k,jz), id_tmp(i,j,k,id$z), &
              ABS((pos(i,j,k,jx) - id_tmp(i,j,k,id$x))/id_tmp(i,j,k,id$x)), &
              ABS((pos(i,j,k,jy) - id_tmp(i,j,k,id$y))/id_tmp(i,j,k,id$y)), &
              ABS((pos(i,j,k,jz) - id_tmp(i,j,k,id$z))/id_tmp(i,j,k,id$z))

          ENDDO
        ENDDO
      ENDDO

      CLOSE( UNIT= 2 )

    ENDIF

    RETURN

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

      PRINT *, "** Subroutine read_fuka_id_spacetime executed."
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

    INTEGER, PARAMETER:: mpi_ranks= 40
    !LOGICAL, PARAMETER:: debug= .FALSE.

    INTEGER:: i, j, k

    DOUBLE PRECISION:: xmin, xmax, ymin, ymax, zmin, zmax

    DOUBLE PRECISION, DIMENSION(nx, ny, nz, 19):: id_tmp

    ! TODO: parallelize these ones
    xmin= MINVAL( pos(:,1,1,jx), DIM= 1 )
    xmax= MAXVAL( pos(:,1,1,jx), DIM= 1 )
    ymin= MINVAL( pos(1,:,1,jy), DIM= 1 )
    ymax= MAXVAL( pos(1,:,1,jy), DIM= 1 )
    zmin= MINVAL( pos(1,1,:,jz), DIM= 1 )
    zmax= MAXVAL( pos(1,1,:,jz), DIM= 1 )

    CALL this% run_kadath_reader( mpi_ranks, nx, ny, nz, &
                                  xmin, xmax, ymin, ymax, zmin, zmax, &
                                  id_tmp(:,:,:,id$x:id$z), &
                                  id_tmp(:,:,:,id$lapse), &
                                  id_tmp(:,:,:,id$shiftx), &
                                  id_tmp(:,:,:,id$shifty), &
                                  id_tmp(:,:,:,id$shiftz), &
                                  id_tmp(:,:,:,id$gxx), &
                                  id_tmp(:,:,:,id$gxy), &
                                  id_tmp(:,:,:,id$gxz), &
                                  id_tmp(:,:,:,id$gyy), &
                                  id_tmp(:,:,:,id$gyz), &
                                  id_tmp(:,:,:,id$gzz), &
                                  id_tmp(:,:,:,id$kxx), &
                                  id_tmp(:,:,:,id$kxy), &
                                  id_tmp(:,:,:,id$kxz), &
                                  id_tmp(:,:,:,id$kyy), &
                                  id_tmp(:,:,:,id$kyz), &
                                  id_tmp(:,:,:,id$kzz), &
                                  baryon_density, &
                                  specific_energy, &
                                  pressure, &
                                  u_euler(:,:,:,jx), &
                                  u_euler(:,:,:,jy), &
                                  u_euler(:,:,:,jz), &
                                  this% filename )

  !  !$OMP PARALLEL DO DEFAULT( NONE ) &
  !  !$OMP             SHARED( nx, ny, nz, this, pos, &
  !  !$OMP                     baryon_density, specific_energy, pressure, &
  !  !$OMP                     energy_density, u_euler, id_tmp ) &
  !  !$OMP             PRIVATE( i, j, k )
  !  DO k= 1, nz, 1
  !    DO j= 1, ny, 1
  !      DO i= 1, nx, 1
  !
  !        baryon_density(i,j,k) = id_tmp(i,j,k,id$massdensity)
  !        specific_energy(i,j,k)= id_tmp(i,j,k,id$specificenergy)
  !        pressure(i,j,k)       = id_tmp(i,j,k,id$pressure)
  !
  !        energy_density(i,j,k) = baryon_density(i,j,k) &
  !                                 *(one + specific_energy(i,j,k))
  !
  !        u_euler(i,j,k,jx)     = id_tmp(i,j,k,id$eulvelx)
  !        u_euler(i,j,k,jy)     = id_tmp(i,j,k,id$eulvely)
  !        u_euler(i,j,k,jz)     = id_tmp(i,j,k,id$eulvelz)
  !
  !      ENDDO
  !    ENDDO
  !  ENDDO
  !  !$OMP END PARALLEL DO

    RETURN

    IF ( C_ASSOCIATED( this% bns_ptr ) ) THEN

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( nx, ny, nz, this, pos, &
      !$OMP                     baryon_density, energy_density, &
      !$OMP                     specific_energy, pressure, u_euler ) &
      !$OMP             PRIVATE( i, j, k )
      coords_z: DO i= 1, nz, 1
        coords_y: DO j= 1, ny, 1
          coords_x: DO k= 1, nx, 1

            CALL get_fuka_id_hydro( this% bns_ptr, &
                                    pos( i, j, k, jx ), &
                                    pos( i, j, k, jy ), &
                                    pos( i, j, k, jz ), &
                                    baryon_density( i, j, k ), &
                                    specific_energy( i, j, k ), &
                                    pressure( i, j, k ), &
                                    u_euler( i, j, k, jx ), &
                                    u_euler( i, j, k, jy ), &
                                    u_euler( i, j, k, jz ) )

            energy_density(i, j, k)= baryon_density(i, j, k) &
                                     *(one + specific_energy(i, j, k))

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
                                    specific_energy(a), &
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

        ! Convert the baryon density and pressure to units of amu
        ! (|sph| code units)
        baryon_density(a)= baryon_density(a)*Msun/amu
        pressure      (a)= pressure(a)*MSun/amu
        energy_density(a)= baryon_density(a)*(one + specific_energy(a))

      ENDDO read_fuka_id_loop
      !$OMP END PARALLEL DO

     ! PRINT *, "energy_density=", energy_density(100)
     ! PRINT *, "baryon_density=", baryon_density(100)
     ! PRINT *, "pressure=", pressure(100)
     ! PRINT *

      ! Convert the baryon density and pressure to units of amu
      ! (|sph| code units)
     ! baryon_density= baryon_density*Msun/amu
     ! pressure      = pressure*MSun/amu
     ! energy_density= baryon_density*(one + specific_energy)

     ! PRINT *, "energy_density/baryon_density=", &
     !           specific_energy(100)/baryon_density(100)
     ! PRINT *, "specific_energy=", energy_density(100)
     ! PRINT *, "baryon_density=", baryon_density(100)
     ! PRINT *, "pressure=", pressure(100)
     ! PRINT *
     ! STOP

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


  MODULE PROCEDURE run_kadath_reader

    !************************************************
    !
    !# Calls the MPI-parallelized version of
    !  the function KadathExportBNS within Kadath
    !
    !  Created:     FT 28.06.2022
    !  Last update: FT 28.06.2022
    !
    !************************************************

#ifdef __INTEL_COMPILER

  USE IFPORT, ONLY: CHANGEDIRQQ

#endif

    IMPLICIT NONE

    INTEGER, PARAMETER:: unit_par = 3480

    INTEGER:: ios, i_char, i_file, i_field, i, j, k
    INTEGER:: nlines, nlines_prev
    INTEGER, DIMENSION(mpi_ranks):: unit_rank
    INTEGER:: unit_rank_prev

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: grid_tmp

    LOGICAL:: exist
    LOGICAL(4):: status

    CHARACTER( LEN= : ), ALLOCATABLE:: filename_par
    CHARACTER( LEN= : ), ALLOCATABLE:: filename_id
    CHARACTER( LEN= : ), ALLOCATABLE:: filename_rank
    CHARACTER( LEN= : ), ALLOCATABLE:: dir_id
    CHARACTER( LEN= 3 ):: mpi_ranks_str

    find_name_loop: DO i_char= LEN(filename), 1, -1

      IF( filename(i_char:i_char) == "/" )THEN
        filename_id= TRIM(filename(i_char+1:LEN(filename)))
        dir_id     = TRIM(filename(1:i_char))
        EXIT
      ENDIF

    ENDDO find_name_loop

    filename_par= "lattice_par.dat"

    INQUIRE( FILE= TRIM(dir_id//filename_par), EXIST= exist )

    IF( exist )THEN
        OPEN( UNIT= unit_par, FILE= TRIM(dir_id//filename_par), &
              STATUS= "REPLACE", FORM= "FORMATTED", &
              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
              IOMSG= err_msg )
    ELSE
        OPEN( UNIT= unit_par, FILE= TRIM(dir_id//filename_par), &
              STATUS= "NEW", FORM= "FORMATTED", &
              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening " // TRIM(dir_id//filename_par), &
               ". The error message is", err_msg
      STOP
    ENDIF

    WRITE( UNIT = unit_par, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      TRIM(filename_id)
    IF( ios > 0 )THEN
      PRINT *, "...error when writing the arrays in ", &
               TRIM(dir_id//filename_par), ". The error message is", err_msg
      STOP
    ENDIF
    WRITE( UNIT = unit_par, IOSTAT = ios, IOMSG = err_msg, FMT = * ) xmin
    IF( ios > 0 )THEN
      PRINT *, "...error when writing the arrays in ", &
               TRIM(dir_id//filename_par), ". The error message is", err_msg
      STOP
    ENDIF
    WRITE( UNIT = unit_par, IOSTAT = ios, IOMSG = err_msg, FMT = * ) xmax
    IF( ios > 0 )THEN
      PRINT *, "...error when writing the arrays in ", &
               TRIM(dir_id//filename_par), ". The error message is", err_msg
      STOP
    ENDIF
    WRITE( UNIT = unit_par, IOSTAT = ios, IOMSG = err_msg, FMT = * ) ymin
    IF( ios > 0 )THEN
      PRINT *, "...error when writing the arrays in ", &
               TRIM(dir_id//filename_par), ". The error message is", err_msg
      STOP
    ENDIF
    WRITE( UNIT = unit_par, IOSTAT = ios, IOMSG = err_msg, FMT = * ) ymax
    IF( ios > 0 )THEN
      PRINT *, "...error when writing the arrays in ", &
               TRIM(dir_id//filename_par), ". The error message is", err_msg
      STOP
    ENDIF
    WRITE( UNIT = unit_par, IOSTAT = ios, IOMSG = err_msg, FMT = * ) zmin
    IF( ios > 0 )THEN
      PRINT *, "...error when writing the arrays in ", &
               TRIM(dir_id//filename_par), ". The error message is", err_msg
      STOP
    ENDIF
    WRITE( UNIT = unit_par, IOSTAT = ios, IOMSG = err_msg, FMT = * ) zmax
    IF( ios > 0 )THEN
      PRINT *, "...error when writing the arrays in ", &
               TRIM(dir_id//filename_par), ". The error message is", err_msg
      STOP
    ENDIF
    WRITE( UNIT = unit_par, IOSTAT = ios, IOMSG = err_msg, FMT = * ) nx
    IF( ios > 0 )THEN
      PRINT *, "...error when writing the arrays in ", &
               TRIM(dir_id//filename_par), ". The error message is", err_msg
      STOP
    ENDIF
    WRITE( UNIT = unit_par, IOSTAT = ios, IOMSG = err_msg, FMT = * ) ny
    IF( ios > 0 )THEN
      PRINT *, "...error when writing the arrays in ", &
               TRIM(dir_id//filename_par), ". The error message is", err_msg
      STOP
    ENDIF
    WRITE( UNIT = unit_par, IOSTAT = ios, IOMSG = err_msg, FMT = * ) nz
    IF( ios > 0 )THEN
      PRINT *, "...error when writing the arrays in ", &
               TRIM(dir_id//filename_par), ". The error message is", err_msg
      STOP
    ENDIF

    CLOSE( UNIT= unit_par )

    ! Change working directory to where the FUKA ID files and the
    ! Kadath reader are stored (they must be in the same directory)

#ifdef __INTEL_COMPILER

  status= CHANGEDIRQQ("/disk/stero-1/ftors/SPHINCS/sphincs_repository/SPHINCS_ID/"//dir_id)
  IF( status == .FALSE. )THEN
  PRINT *, "** ERROR! Unable to change directory in SUBROUTINE ", &
           "set_up_lattices_around_stars!"
  PRINT *, " * Stopping..."
  PRINT *
  STOP
  ENDIF

#endif

#ifdef __GFORTRAN__

  CALL CHDIR(dir_id)

#endif

    ! Run the MPI parallelized Kadath reader
    IF( mpi_ranks <= 9   ) WRITE( mpi_ranks_str, '(I1)' ) mpi_ranks
    IF( mpi_ranks >= 10  ) WRITE( mpi_ranks_str, '(I2)' ) mpi_ranks
    IF( mpi_ranks >= 100 ) WRITE( mpi_ranks_str, '(I3)' ) mpi_ranks
    CALL EXECUTE_COMMAND_LINE("mpirun -np "//TRIM(mpi_ranks_str)// &
                              " export_bns_test")

    ! Delete the parameter file that specifies the lattice
    CALL EXECUTE_COMMAND_LINE("rm -f "//TRIM(filename_par))

    ! Allocate memory
    IF(.NOT.ALLOCATED(grid_tmp)) ALLOCATE(grid_tmp(nx*ny*nz, n_fields_fuka))

    ! Read the ID from the ASCII files printed by the reader
    ! (one per MPI rank)
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( grid_tmp, mpi_ranks ) &
    !$OMP             PRIVATE( i_file, i, filename_rank, mpi_ranks_str, &
    !$OMP                      unit_rank, nlines, ios, nlines_prev, exist, &
    !$OMP                      err_msg, unit_rank_prev )
    loop_over_id_files: DO i_file= 0, mpi_ranks - 1, 1

      IF( i_file <= 9   ) WRITE( mpi_ranks_str, '(I1)' ) i_file
      IF( i_file >= 10  ) WRITE( mpi_ranks_str, '(I2)' ) i_file
      IF( i_file >= 100 ) WRITE( mpi_ranks_str, '(I3)' ) i_file
      filename_rank= "id-"//TRIM(mpi_ranks_str)//".dat"

      unit_rank(i_file + 1)= 8346 + i_file

      INQUIRE( FILE= TRIM(filename_rank), EXIST= exist )
      IF( exist )THEN
        OPEN( unit_rank(i_file + 1), FILE= TRIM(filename_rank), &
              FORM= "FORMATTED", ACTION= "READ" )
      ELSE
        PRINT *
        PRINT *, "** ERROR: ", TRIM(filename_rank), &
                 " file not found!"
        PRINT *
        STOP
      ENDIF

      ! Get total number of lines in the file
      nlines = 0
      DO
        READ( unit_rank(i_file + 1), * , IOSTAT= ios )
        IF ( ios /= 0 ) EXIT
        nlines = nlines + 1
      ENDDO
      CLOSE( UNIT= unit_rank(i_file + 1) )
      OPEN( unit_rank(i_file + 1), FILE= TRIM(filename_rank), &
            FORM= "FORMATTED", ACTION= "READ" )

      ! Neglect header
      nlines= nlines - 1

      !
      IF( i_file == mpi_ranks - 1 )THEN

        unit_rank_prev= unit_rank(i_file + 1) + 1

        OPEN( unit_rank_prev, FILE= TRIM("id-0.dat"), &
              FORM= "FORMATTED", ACTION= "READ" )

        nlines_prev = 0
        DO
          READ( unit_rank_prev, * , IOSTAT= ios )
          IF ( ios /= 0 ) EXIT
          nlines_prev = nlines_prev + 1
        ENDDO

        CLOSE(unit_rank_prev)

        ! Neglect header
        nlines_prev= nlines_prev - 1

      ELSE

        nlines_prev= nlines

      ENDIF

      ! Skip header
      READ( unit_rank(i_file + 1), * )
      DO i= 1, nlines, 1
        READ( UNIT= unit_rank(i_file + 1), FMT= *, IOSTAT = ios, &
              IOMSG= err_msg ) grid_tmp( nlines_prev*i_file + i, : )
      ENDDO

      ! Close file and delete it
      CLOSE( unit_rank(i_file + 1) )
      CALL EXECUTE_COMMAND_LINE("rm -f "//TRIM(filename_rank))

    ENDDO loop_over_id_files
    !$OMP END PARALLEL DO

    ! Store fields in desired format (needed by trilinear_interpolation
    ! in MODULE numerics)

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( this, grid_tmp, nx, ny, nz, coords, lapse, &
    !$OMP                     shift_x, shift_y, shift_z, g_xx, g_xy, g_xz, &
    !$OMP                     g_yy, g_yz, g_zz, k_xx, k_xy, k_xz, k_yy, k_yz, &
    !$OMP                     k_zz, mass_density, specific_energy, pressure, &
    !$OMP                     v_eul_x, v_eul_y, v_eul_z ) &
    !$OMP             PRIVATE( i, j, k )
    DO k= 1, nz, 1
      DO j= 1, ny, 1
        DO i= 1, nx, 1
          !DO i_field= 1, n_fields_fuka
          !
          !  id_fields( i, j, k, i_field )= &
          !        grid_tmp( (k-1)*ny*nx + (j-1)*nx + i, i_field )
          !
          !ENDDO
          coords    (i,j,k,id$x)= grid_tmp( (k-1)*ny*nx + (j-1)*nx + i, id$x )
          coords    (i,j,k,id$y)= grid_tmp( (k-1)*ny*nx + (j-1)*nx + i, id$y )
          coords    (i,j,k,id$z)= grid_tmp( (k-1)*ny*nx + (j-1)*nx + i, id$z )
          lapse          (i,j,k)= grid_tmp( (k-1)*ny*nx + (j-1)*nx + i, &
            id$lapse )
          shift_x        (i,j,k)= grid_tmp( (k-1)*ny*nx + (j-1)*nx + i, &
            id$shiftx )
          shift_y        (i,j,k)= grid_tmp( (k-1)*ny*nx + (j-1)*nx + i, &
            id$shifty )
          shift_z        (i,j,k)= grid_tmp( (k-1)*ny*nx + (j-1)*nx + i, &
            id$shiftz )
          g_xx           (i,j,k)= grid_tmp( (k-1)*ny*nx + (j-1)*nx + i, &
            id$gxx )
          g_xy           (i,j,k)= grid_tmp( (k-1)*ny*nx + (j-1)*nx + i, &
            id$gxy )
          g_xz           (i,j,k)= grid_tmp( (k-1)*ny*nx + (j-1)*nx + i, &
            id$gxz )
          g_yy           (i,j,k)= grid_tmp( (k-1)*ny*nx + (j-1)*nx + i, &
            id$gyy )
          g_yz           (i,j,k)= grid_tmp( (k-1)*ny*nx + (j-1)*nx + i, &
            id$gyz )
          g_zz           (i,j,k)= grid_tmp( (k-1)*ny*nx + (j-1)*nx + i, &
            id$gzz )
          mass_density   (i,j,k)= grid_tmp( (k-1)*ny*nx + (j-1)*nx + i, &
            id$massdensity )
          specific_energy(i,j,k)= grid_tmp( (k-1)*ny*nx + (j-1)*nx + i, &
            id$specificenergy )
          pressure       (i,j,k)= grid_tmp( (k-1)*ny*nx + (j-1)*nx + i, &
            id$pressure  )
          v_eul_x        (i,j,k)= grid_tmp( (k-1)*ny*nx + (j-1)*nx + i, &
            id$eulvelx )
          v_eul_y        (i,j,k)= grid_tmp( (k-1)*ny*nx + (j-1)*nx + i, &
            id$eulvely )
          v_eul_z        (i,j,k)= grid_tmp( (k-1)*ny*nx + (j-1)*nx + i, &
            id$eulvelz )
        ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Change working directory back to HOME_SPHINCS_ID

#ifdef __INTEL_COMPILER

  status= CHANGEDIRQQ("/disk/stero-1/ftors/SPHINCS/sphincs_repository/SPHINCS_ID/")
  IF( status == .FALSE. )THEN
  PRINT *, "** ERROR! Unable to change directory in SUBROUTINE ", &
       "set_up_lattices_around_stars!"
  PRINT *, " * Stopping..."
  PRINT *
  STOP
  ENDIF

#endif

#ifdef __GFORTRAN__

  CALL CHDIR("/disk/stero-1/ftors/SPHINCS/sphincs_repository/SPHINCS_ID/")

#endif

  END PROCEDURE run_kadath_reader


  MODULE PROCEDURE set_up_lattices_around_stars

    !***********************************************
    !
    !#
    !
    !  FT 27.06.2022
    !
    !***********************************************


    IMPLICIT NONE


    DOUBLE PRECISION, PARAMETER:: stretch= 1.05D0
    !! The lattices' sizes will be 5% larger than the radii of the stars

    INTEGER:: nx
    INTEGER:: ny
    INTEGER:: nz
    INTEGER:: i_star, i, j, k
    INTEGER:: mpi_ranks
    !

    DOUBLE PRECISION:: xmin, xmax, ymin, ymax, zmin, zmax
    DOUBLE PRECISION, DIMENSION(6):: sizes
    DOUBLE PRECISION, DIMENSION(3):: center

    LOGICAL:: exist

    CHARACTER(LEN=:), ALLOCATABLE:: filename_id

#ifdef MPI_ranks

  mpi_ranks= MPI_ranks

#else

  PRINT *, "** ERROR! No value assigned to the variable MPI_ranks in the ", &
           "SConstruct file! Please assign a value to it!"
  PRINT *, " * Stopping..."
  PRINT *
  STOP

#endif

    nx= this% nx_grid
    ny= this% ny_grid
    nz= this% nz_grid

    ! Allocate and initialize member arrays
    !ALLOCATE( this% id_fields( nx, ny, nz, n_fields_fuka, 2 ) )

    loop_over_stars: DO i_star= 1, 2, 1

      CALL this% star_lattice(i_star)% allocate_lattice_memory(nx,ny,nz)

      sizes = this% return_spatial_extent(i_star)
      center= this% return_center(i_star)

      ! Determine boundaries of the lattices
      xmin= center(1) - stretch*sizes(1)
      xmax= center(1) + stretch*sizes(2)
      ymin= center(2) - stretch*sizes(3)
      ymax= center(2) + stretch*sizes(4)
      zmin= center(3) - stretch*sizes(5)
      zmax= center(3) + stretch*sizes(6)

      CALL this% run_kadath_reader( mpi_ranks, nx, ny, nz                    , &
                                  xmin, xmax, ymin, ymax, zmin, zmax         , &
                                  this% star_lattice(i_star)% coords         , &
                                  this% star_lattice(i_star)% lapse          , &
                                  this% star_lattice(i_star)% shift_x        , &
                                  this% star_lattice(i_star)% shift_y        , &
                                  this% star_lattice(i_star)% shift_z        , &
                                  this% star_lattice(i_star)% g_xx           , &
                                  this% star_lattice(i_star)% g_xy           , &
                                  this% star_lattice(i_star)% g_xz           , &
                                  this% star_lattice(i_star)% g_yy           , &
                                  this% star_lattice(i_star)% g_yz           , &
                                  this% star_lattice(i_star)% g_zz           , &
                                  this% star_lattice(i_star)% k_xx           , &
                                  this% star_lattice(i_star)% k_xy           , &
                                  this% star_lattice(i_star)% k_xz           , &
                                  this% star_lattice(i_star)% k_yy           , &
                                  this% star_lattice(i_star)% k_yz           , &
                                  this% star_lattice(i_star)% k_zz           , &
                                  this% star_lattice(i_star)% mass_density   , &
                                  this% star_lattice(i_star)% specific_energy, &
                                  this% star_lattice(i_star)% pressure       , &
                                  this% star_lattice(i_star)% v_eul_x        , &
                                  this% star_lattice(i_star)% v_eul_y        , &
                                  this% star_lattice(i_star)% v_eul_z        , &
                                  this% filename )

      PRINT *, "** ID stored on a fine lattice around star ", i_star
      PRINT *

    ENDDO loop_over_stars

  !  filename_id= "dbg-id.dat"
  !
  !  INQUIRE( FILE= TRIM(filename_id), EXIST= exist )
  !
  !  IF( exist )THEN
  !    OPEN( UNIT= 2, FILE= TRIM(filename_id), STATUS= "REPLACE", &
  !          FORM= "FORMATTED", &
  !          POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
  !          IOMSG= err_msg )
  !  ELSE
  !    OPEN( UNIT= 2, FILE= TRIM(filename_id), STATUS= "NEW", &
  !          FORM= "FORMATTED", &
  !          ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
  !  ENDIF
  !  IF( ios > 0 )THEN
  !  PRINT *, "...error when opening " // TRIM(filename_id), &
  !           ". The error message is", err_msg
  !  STOP
  !  ENDIF
  !
  !  DO i_star= 1, 2, 1
  !    DO k= 1, nz, 1
  !      DO j= 1, ny, 1
  !        DO i= 1, nx, 1
  !
  !          WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
  !            this% id_fields( i, j, k, :, i_star )
  !
  !        ENDDO
  !      ENDDO
  !    ENDDO
  !  ENDDO
  !
  !  CLOSE( UNIT= 2 )

    !STOP

  END PROCEDURE set_up_lattices_around_stars


END SUBMODULE read
