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


  !USE constants, ONLY: Msun_geo


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

    !IF ( C_ASSOCIATED( THIS% bns_ptr ) ) THEN
    !
    !  IF( SIZE( x ) /= SIZE( y ) .OR. SIZE( x ) /= SIZE( z ) &
    !        .OR. SIZE( y ) /= SIZE( z ) )THEN
    !    PRINT *, "** ERROR: The sizes of the arrays of positions" &
    !             // "passed to read_lorene_id are not the same."
    !    PRINT *
    !    STOP
    !  ENDIF
    !
    !  IF( ALLOCATED( THIS% lapse )   .AND. &
    !      ALLOCATED( THIS% shift_x ) .AND. &
    !      ALLOCATED( THIS% shift_y ) .AND. &
    !      ALLOCATED( THIS% shift_z ) .AND. &
    !      ALLOCATED( THIS% g_xx ) .AND. ALLOCATED( THIS% g_xy ) .AND. &
    !      ALLOCATED( THIS% g_xz ) .AND. ALLOCATED( THIS% g_yy ) .AND. &
    !      ALLOCATED( THIS% g_yz ) .AND. ALLOCATED( THIS% g_zz ) .AND. &
    !      ALLOCATED( THIS% k_xx ) .AND. ALLOCATED( THIS% k_xy ) .AND. &
    !      ALLOCATED( THIS% k_xz ) .AND. ALLOCATED( THIS% k_yy ) .AND. &
    !      ALLOCATED( THIS% k_yz ) .AND. ALLOCATED( THIS% k_zz ) .AND. &
    !      ALLOCATED( THIS% baryon_density )  .AND. &
    !      ALLOCATED( THIS% energy_density )  .AND. &
    !      ALLOCATED( THIS% specific_energy ) .AND. &
    !      ALLOCATED( THIS% v_euler_x )       .AND. &
    !      ALLOCATED( THIS% v_euler_y )       .AND. &
    !      ALLOCATED( THIS% v_euler_z ) &
    !  )THEN
    !
    !    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !    !$OMP             SHARED( n, THIS, x, y, z ) &
    !    !$OMP             PRIVATE( itr )
    !    read_fuka_id_loop: DO itr= 1, n, 1
    !
    !      ! The coordinates need to be converted from |sphincs| units (Msun_geo)
    !      ! to |fuka| units (\(\mathrm{km}\)). See MODULE constants for the
    !      ! definition of Msun_geo
    !      CALL get_lorene_id( THIS% bns_ptr, &
    !                          x( itr )*Msun_geo, &
    !                          y( itr )*Msun_geo, &
    !                          z( itr )*Msun_geo, &
    !                          THIS% lapse( itr ), &
    !                          THIS% shift_x( itr ), &
    !                          THIS% shift_y( itr ), &
    !                          THIS% shift_z( itr ), &
    !                          THIS% g_xx( itr ), &
    !                          THIS% k_xx( itr ), &
    !                          THIS% k_xy( itr ), &
    !                          THIS% k_xz( itr ), &
    !                          THIS% k_yy( itr ), &
    !                          THIS% k_yz( itr ), &
    !                          THIS% k_zz( itr ), &
    !                          THIS% baryon_density( itr ), &
    !                          THIS% energy_density( itr ), &
    !                          THIS% specific_energy( itr ), &
    !                          THIS% v_euler_x( itr ), &
    !                          THIS% v_euler_y( itr ), &
    !                          THIS% v_euler_z( itr ) )
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
    !      THIS% g_yy( itr )= THIS% g_xx( itr )
    !      THIS% g_zz( itr )= THIS% g_xx( itr )
    !      THIS% g_xy( itr )= 0.0D0
    !      THIS% g_xz( itr )= 0.0D0
    !      THIS% g_yz( itr )= 0.0D0
    !
    !      !
    !      !- Set/unset the geodesic gauge
    !      !
    !      IF( THIS% get_one_lapse() )THEN
    !        THIS% lapse( itr )= 1.0D0
    !      ENDIF
    !      IF( THIS% get_zero_shift() )THEN
    !        THIS% shift_x( itr )= 0.0D0
    !        THIS% shift_y( itr )= 0.0D0
    !        THIS% shift_z( itr )= 0.0D0
    !      ENDIF
    !
    !      !
    !      !-- Convert the extrinsic curvature from |fuka| units to
    !      !-- |sphincs| units
    !      !
    !      THIS% k_xx( itr )= THIS% k_xx( itr )*Msun_geo
    !      THIS% k_xy( itr )= THIS% k_xy( itr )*Msun_geo
    !      THIS% k_xz( itr )= THIS% k_xz( itr )*Msun_geo
    !      THIS% k_yy( itr )= THIS% k_yy( itr )*Msun_geo
    !      THIS% k_yz( itr )= THIS% k_yz( itr )*Msun_geo
    !      THIS% k_zz( itr )= THIS% k_zz( itr )*Msun_geo
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
    !  FT 09.02.2022
    !
    !**************************************************

    IMPLICIT NONE

   ! IF ( C_ASSOCIATED( THIS% bns_ptr ) ) THEN
   !
   !   IF( SIZE( x ) /= SIZE( y ) .OR. SIZE( x ) /= SIZE( z ) &
   !         .OR. SIZE( y ) /= SIZE( z ) )THEN
   !     PRINT *, "** ERROR: The sizes of the arrays of positions" &
   !              // "passed to read_lorene_id are not the same."
   !     PRINT *
   !     STOP
   !   ENDIF
   !
   !   !$OMP PARALLEL DO DEFAULT( NONE ) &
   !   !$OMP             SHARED( n, THIS, x, y, z, lapse, &
   !   !$OMP                     shift_x, shift_y, shift_z, &
   !   !$OMP                     g_xx, k_xx, k_xy, k_xz, k_yy, k_yz, k_zz, &
   !   !$OMP                     baryon_density, energy_density, &
   !   !$OMP                     specific_energy, &
   !   !$OMP                     u_euler_x, u_euler_y, u_euler_z ) &
   !   !$OMP             PRIVATE( itr )
   !   read_fuka_id_loop: DO itr= 1, n, 1
   !
   !     ! The coordinates need to be converted from |sphincs| units (Msun_geo)
   !     ! to |fuka| units (\(\mathrm{km}\)). See MODULE constants for the definition of
   !     ! Msun_geo
   !     CALL get_lorene_id( THIS% bns_ptr, &
   !                         x( itr )*Msun_geo, &
   !                         y( itr )*Msun_geo, &
   !                         z( itr )*Msun_geo, &
   !                         lapse( itr ), &
   !                         shift_x( itr ), &
   !                         shift_y( itr ), &
   !                         shift_z( itr ), &
   !                         g_xx( itr ), &
   !                         k_xx( itr ), &
   !                         k_xy( itr ), &
   !                         k_xz( itr ), &
   !                         k_yy( itr ), &
   !                         k_yz( itr ), &
   !                         k_zz( itr ), &
   !                         baryon_density( itr ), &
   !                         energy_density( itr ), &
   !                         specific_energy( itr ), &
   !                         u_euler_x( itr ), &
   !                         u_euler_y( itr ), &
   !                         u_euler_z( itr ) )
   !
   !   ENDDO read_fuka_id_loop
   !   !$OMP END PARALLEL DO
   !
   !   DO itr= 1, n, 1
   !
   !     !
   !     !-- The following follows from the assumption of conformal
   !     !-- flatness in |fuka|
   !     !
   !     g_yy( itr )= g_xx( itr )
   !     g_zz( itr )= g_xx( itr )
   !     g_xy( itr )= 0.0D0
   !     g_xz( itr )= 0.0D0
   !     g_yz( itr )= 0.0D0
   !
   !     !
   !     !- Set/unset the geodesic gauge
   !     !
   !     IF( THIS% get_one_lapse() )THEN
   !       lapse( itr )= 1.0D0
   !     ENDIF
   !     IF( THIS% get_zero_shift() )THEN
   !       shift_x( itr )= 0.0D0
   !       shift_y( itr )= 0.0D0
   !       shift_z( itr )= 0.0D0
   !     ENDIF
   !
   !     !
   !     !-- Convert the extrinsic curvature from |fuka| units to
   !     !-- |sphincs| units
   !     !
   !     k_xx( itr )= k_xx( itr )*Msun_geo
   !     k_xy( itr )= k_xy( itr )*Msun_geo
   !     k_xz( itr )= k_xz( itr )*Msun_geo
   !     k_yy( itr )= k_yy( itr )*Msun_geo
   !     k_yz( itr )= k_yz( itr )*Msun_geo
   !     k_zz( itr )= k_zz( itr )*Msun_geo
   !
   !     ! Print progress on screen
   !     perc= 100*itr/n
   !     IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
   !       WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
   !               creturn//" ", perc, "%"
   !     ENDIF
   !
   !   ENDDO
   !   IF( show_progress ) WRITE( *, "(A1)", ADVANCE= "NO" ) creturn
   !
   !   PRINT *, "** Subroutine read_lorene_id executed."
   !   PRINT *
   !
   ! ENDIF

  END PROCEDURE read_fuka_id_full


  MODULE PROCEDURE read_fuka_id_spacetime

    !*******************************************************
    !
    !# Stores the spacetime |id| in multi-dimensional arrays
    !  needed to compute the BSSN variables and constraints
    !
    !  FT 09.02.2022
    !
    !*******************************************************

    USE tensor,    ONLY: jxx, jxy, jxz, &
                         jyy, jyz, jzz, jx, jy, jz, n_sym4x4

    IMPLICIT NONE

 !   INTEGER:: i, j, k
 !
 !   DOUBLE PRECISION:: detg
 !   DOUBLE PRECISION:: detg4
 !   DOUBLE PRECISION, DIMENSION( :, :, :, : ), ALLOCATABLE:: g4
 !
 !   ! g4 is allocatable to allocate it on the heap
 !   ! Allocating it on the stack might exceed stack memory,
 !   ! causing a segmentation fault
 !   ALLOCATE( g4( nx, ny, nz, n_sym4x4 ) )
 !
 !   IF ( C_ASSOCIATED( THIS% bns_ptr ) ) THEN
 !
 !     IF( .FALSE. &!SHAPE( pos(:,:,:,1) ) /= SHAPE( lapse ) .OR. &
 !         !SHAPE( pos(:,:,:,1) ) /= SHAPE( shift(:,:,:,jx) ) & ! .OR. &
 !       ! SHAPE( pos(:,:,:,1) ) /= SHAPE( g(:,:,:,1) ) .OR. &
 !       ! SHAPE( pos(:,:,:,1) ) /= SHAPE( k(:,:,:,1) ) &
 !       )THEN
 !       PRINT *, "** ERROR: Mismatch in array dimensions" &
 !                // "in read_fuka_id_spacetime."
 !       PRINT *
 !       STOP
 !     ENDIF
 !
 !     !$OMP PARALLEL DO DEFAULT( NONE ) &
 !     !$OMP             SHARED( nx, ny, nz, THIS, pos, &
 !     !$OMP                     lapse, shift, g, ek ) &
 !     !$OMP             PRIVATE( i, j, k )
 !     coords_z: DO k= 1, nz, 1
 !       coords_y: DO j= 1, ny, 1
 !         coords_x: DO i= 1, nx, 1
 !
 !           ! The coordinates need to be converted from |sphincs| units (Msun_geo)
 !           ! to |fuka| units (\(\mathrm{km}\)). See MODULE constants for the definition of
 !           ! Msun_geo
 !           CALL get_lorene_id_spacetime( THIS% bns_ptr, &
 !                               pos( i, j, k, jx )*Msun_geo, &
 !                               pos( i, j, k, jy )*Msun_geo, &
 !                               pos( i, j, k, jz )*Msun_geo, &
 !                               lapse( i, j, k ), &
 !                               shift( i, j, k, jx ), &
 !                               shift( i, j, k, jy ), &
 !                               shift( i, j, k, jz ), &
 !                               g( i, j, k, jxx ), &
 !                               ek( i, j, k, jxx ), &
 !                               ek( i, j, k, jxy ), &
 !                               ek( i, j, k, jxz ), &
 !                               ek( i, j, k, jyy ), &
 !                               ek( i, j, k, jyz ), &
 !                               ek( i, j, k, jzz ) )
 !
 !         ENDDO coords_x
 !       ENDDO coords_y
 !     ENDDO coords_z
 !     !$OMP END PARALLEL DO
 !
 !     DO k= 1, nz, 1
 !       DO j= 1, ny, 1
 !         DO i= 1, nx, 1
 !
 !           !
 !           !-- The following follows from the assumption of
 !           !-- conformal flatness in |fuka|
 !           !
 !           g( i, j, k, jyy )= g( i, j, k, jxx )
 !           g( i, j, k, jzz )= g( i, j, k, jxx )
 !           g( i, j, k, jxy )= 0.0D0
 !           g( i, j, k, jxz )= 0.0D0
 !           g( i, j, k, jyz )= 0.0D0
 !
 !           !
 !           !- Set/unset the geodesic gauge
 !           !
 !           IF( THIS% get_one_lapse() )THEN
 !             lapse( i, j, k )= 1.0D0
 !           ENDIF
 !           IF( THIS% get_zero_shift() )THEN
 !             shift( i, j, k, jx )= 0.0D0
 !             shift( i, j, k, jy )= 0.0D0
 !             shift( i, j, k, jz )= 0.0D0
 !           ENDIF
 !
 !           !
 !           !-- Convert the extrinsic curvature from |fuka| units to
 !           !-- |sphincs| units
 !           !
 !           ek( i, j, k, jxx )= ek( i, j, k, jxx )*Msun_geo
 !           ek( i, j, k, jxy )= ek( i, j, k, jxy )*Msun_geo
 !           ek( i, j, k, jxz )= ek( i, j, k, jxz )*Msun_geo
 !           ek( i, j, k, jyy )= ek( i, j, k, jyy )*Msun_geo
 !           ek( i, j, k, jyz )= ek( i, j, k, jyz )*Msun_geo
 !           ek( i, j, k, jzz )= ek( i, j, k, jzz )*Msun_geo
 !
 !           detg= 2.0D0*g(i,j,k,jxy)*g(i,j,k,jxz)*g(i,j,k,jyz) &
 !                 - g(i,j,k,jzz)*g(i,j,k,jxy)**2 + g(i,j,k,jyy) &
 !                  *( g(i,j,k,jxx)*g(i,j,k,jzz) - g(i,j,k,jxz)**2 ) &
 !                 - g(i,j,k,jxx)*g(i,j,k,jyz)**2
 !
 !           IF( ABS( detg ) < 1D-10 )THEN
 !             PRINT *, "The determinant of the spatial metric " &
 !                      // "is effectively 0 at the grid point " &
 !                      // "(ix,iy,iz)= (", i, ",", j,",", k, ")."
 !             PRINT *, "detg=", detg
 !             PRINT *
 !             STOP
 !           ELSEIF( detg < 0 )THEN
 !             PRINT *, "The determinant of the spatial metric " &
 !                      // "is negative at the grid point " &
 !                      // "(ix,iy,iz)= (", i, ",", j,",", k, ")."
 !             PRINT *, "detg=", detg
 !             PRINT *
 !             STOP
 !           ENDIF
 !
 !           CALL compute_g4( lapse(i,j,k), shift(i,j,k,:), &
 !                            g(i,j,k,:), g4(i,j,k,:) )
 !
 !           CALL determinant_sym4x4( g4(i,j,k,:), detg4 )
 !
 !           IF( ABS( detg4 ) < 1D-10 )THEN
 !             PRINT *, "The determinant of the spacetime metric "&
 !                      // "is effectively 0 at the grid point " &
 !                      // "(ix,iy,iz)= (", i, ",", j,",", k, ")."
 !             PRINT *, "detg4=", detg4
 !             PRINT *
 !             STOP
 !           ELSEIF( detg4 > 0 )THEN
 !             PRINT *, "The determinant of the spacetime metric "&
 !                      // "is positive at the grid point " &
 !                      // "(ix,iy,iz)= (", i, ",", j,",", k, ")."
 !             PRINT *, "detg4=", detg4
 !             PRINT *
 !             STOP
 !           ENDIF
 !
 !           ! Print progress on screen
 !           perc= 100*( nx*ny*(k - 1) + nx*(j - 1) + i )/( nx*ny*nz )
 !           !perc2= 100.0*DBLE(nx*ny*(iz - 1) + nx*(iy - 1) + ix) &
 !           !       /DBLE( nx*ny*nz )
 !           !perc= 100*cnt/( nx*ny*nz )
 !           IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
 !             WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
 !                     creturn//" ", perc, "%"
 !             !WRITE( *, "(A2,F5.2,A1)", ADVANCE= "NO" ) &
 !             !        creturn//" ", perc2, "%"
 !           ENDIF
 !
 !         ENDDO
 !       ENDDO
 !     ENDDO
 !     IF( show_progress ) WRITE( *, "(A1)", ADVANCE= "NO" ) creturn
 !
 !     PRINT *, "** Subroutine read_lorene_id executed."
 !     PRINT *
 !
 !   ENDIF

  END PROCEDURE read_fuka_id_spacetime


  MODULE PROCEDURE read_fuka_id_hydro

    !*******************************************************
    !
    !# Stores the hydro |id| in the arrays needed to compute
    !  the constraints on the refined mesh
    !
    !  FT 09.02.2022
    !
    !*******************************************************

    USE tensor,     ONLY: jx, jy, jz

    IMPLICIT NONE

   ! INTEGER:: ix, iy, iz
   !
   ! IF ( C_ASSOCIATED( THIS% bns_ptr ) ) THEN
   !
   !   !$OMP PARALLEL DO DEFAULT( NONE ) &
   !   !$OMP             SHARED( nx, ny, nz, THIS, pos, &
   !   !$OMP                     baryon_density, energy_density, &
   !   !$OMP                     specific_energy, pressure, u_euler ) &
   !   !$OMP             PRIVATE( ix, iy, iz )
   !   coords_z: DO iz= 1, nz, 1
   !     coords_y: DO iy= 1, ny, 1
   !       coords_x: DO ix= 1, nx, 1
   !
   !         ! The coordinates need to be converted from |sphincs| units (Msun_geo)
   !         ! to |fuka| units (\(\mathrm{km}\)). See MODULE constants for the definition of
   !         ! Msun_geo
   !         CALL get_lorene_id_hydro( THIS% bns_ptr, &
   !                           pos( ix, iy, iz, jx )*Msun_geo, &
   !                           pos( ix, iy, iz, jy )*Msun_geo, &
   !                           pos( ix, iy, iz, jz )*Msun_geo, &
   !                           baryon_density( ix, iy, iz ), &
   !                           energy_density( ix, iy, iz ), &
   !                           specific_energy( ix, iy, iz ), &
   !                           pressure( ix, iy, iz ), &
   !                           u_euler( ix, iy, iz, jx ), &
   !                           u_euler( ix, iy, iz, jy ), &
   !                           u_euler( ix, iy, iz, jz ) )
   !
   !       ENDDO coords_x
   !     ENDDO coords_y
   !   ENDDO coords_z
   !   !$OMP END PARALLEL DO
   !
   !   !      ! Print progress on screen
   !   !      perc= 100*(nx*ny*(iz - 1) &
   !   !            + nx*(iy - 1) + ix)/( nx*ny*nz )
   !   !      IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
   !   !        WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
   !   !                creturn//" ", perc, "%"
   !   !      ENDIF
   !   !
   !   !    ENDDO coords_x
   !   !  ENDDO coords_y
   !   !ENDDO coords_z
   !   !IF( show_progress ) WRITE( *, "(A1)", ADVANCE= "NO" ) creturn
   !
   !   PRINT *, "** Subroutine read_lorene_id_hydro executed."
   !   PRINT *
   !
   ! ENDIF

  END PROCEDURE read_fuka_id_hydro


  MODULE PROCEDURE read_fuka_id_particles

    !****************************************************
    !
    !# Stores the hydro |id| in the arrays needed to
    !  compute the |sph| |id|
    !
    !  FT 09.02.2022
    !
    !****************************************************

    USE constants,  ONLY: amu
    USE utility,    ONLY: Msun_geo, km2m, g2kg

    IMPLICIT NONE

!    DOUBLE PRECISION:: detg
!
!    IF ( C_ASSOCIATED( THIS% bns_ptr ) ) THEN
!
!      IF( SIZE( x ) /= SIZE( y ) .OR. SIZE( x ) /= SIZE( z ) &
!              .OR. SIZE( y ) /= SIZE( z ) )THEN
!        PRINT *, "** ERROR: The sizes of the arrays of positions" &
!                 // "passed to read_lorene_id are not the same."
!        PRINT *
!        STOP
!      ENDIF
!
!      PRINT *, "** Importing |id| on particles..."
!
!      !$OMP PARALLEL DO DEFAULT( NONE ) &
!      !$OMP             SHARED( n, THIS, x, y, z, lapse, &
!      !$OMP                     shift_x, shift_y, shift_z, &
!      !$OMP                     g_xx, &
!      !$OMP                     baryon_density, energy_density, &
!      !$OMP                     specific_energy, pressure, &
!      !$OMP                     u_euler_x, u_euler_y, u_euler_z ) &
!      !$OMP             PRIVATE( itr )
!      read_fuka_id_loop: DO itr= 1, n, 1
!
!        ! The coordinates need to be converted from |sphincs| units (Msun_geo)
!        ! to |fuka| units (\(\mathrm{km}\)). See MODULE constants for the definition of
!        ! Msun_geo
!        CALL get_lorene_id_particles( THIS% bns_ptr, &
!                                      x( itr )*Msun_geo, &
!                                      y( itr )*Msun_geo, &
!                                      z( itr )*Msun_geo, &
!                                      lapse( itr ), &
!                                      shift_x( itr ), &
!                                      shift_y( itr ), &
!                                      shift_z( itr ), &
!                                      g_xx( itr ), &
!                                      baryon_density( itr ), &
!                                      energy_density( itr ), &
!                                      specific_energy( itr ), &
!                                      pressure( itr ), &
!                                      u_euler_x( itr ), &
!                                      u_euler_y( itr ), &
!                                      u_euler_z( itr ) )
!
!      ENDDO read_fuka_id_loop
!      !$OMP END PARALLEL DO
!
!      DO itr= 1, n, 1
!
!        !
!        !-- The following follows from the assumption of conformal
!        !-- flatness in |fuka|
!        !
!        g_yy( itr )= g_xx( itr )
!        g_zz( itr )= g_xx( itr )
!        g_xy( itr )= 0.0D0
!        g_xz( itr )= 0.0D0
!        g_yz( itr )= 0.0D0
!
!        !
!        !- Set/unset the geodesic gauge
!        !
!        IF( THIS% get_one_lapse() )THEN
!          lapse( itr )= 1.0D0
!        ENDIF
!        IF( THIS% get_zero_shift() )THEN
!          shift_x( itr )= 0.0D0
!          shift_y( itr )= 0.0D0
!          shift_z( itr )= 0.0D0
!        ENDIF
!
!        detg= 2*g_xy(itr)*g_xz(itr)*g_yz(itr) &
!              - g_zz(itr)*g_xy(itr)**2 &
!              + g_yy(itr)*( g_xx(itr)*g_zz(itr) - g_xz(itr)**2 ) &
!              - g_xx(itr)*g_yz(itr)**2
!
!        IF( ABS( detg ) < 1D-10 )THEN
!          PRINT *, "The determinant of the spatial metric is " &
!                   // "effectively 0 at the particle ", itr
!          PRINT *, "detg=", detg
!          PRINT *
!          STOP
!        ELSEIF( detg < 0 )THEN
!          PRINT *, "The determinant of the spatial metric is " &
!                   // "negative at the particle ", itr
!          PRINT *, "detg=", detg
!          PRINT *
!          STOP
!        ENDIF
!
!        ! Print progress on screen
!        perc= 100*itr/n
!        IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
!          WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
!                  creturn//" ", perc, "%"
!        ENDIF
!
!      ENDDO
!      IF( show_progress ) WRITE( *, "(A1)", ADVANCE= "NO" ) creturn
!
!      ! Convert the baryon density and pressure to units of amu (|sph| code units)
!      baryon_density= baryon_density*((Msun_geo*km2m)**3)/(amu*g2kg)
!      pressure      = pressure*((Msun_geo*km2m)**3)/(amu*g2kg)
!
!      PRINT *, "** Subroutine read_fuka_id_particles executed."
!      PRINT *
!
!    ENDIF

  END PROCEDURE read_fuka_id_particles


  MODULE PROCEDURE read_fuka_id_mass_b

    !****************************************************
    !
    !# Stores the hydro |id| in the arrays needed to
    !  compute the baryon mass, storing it to variables
    !  (not arrays as the others SUBROUTINES in
    !  the [[bns_read]] SUBMODULE).
    !
    !  FT 09.02.2022
    !
    !****************************************************

    USE utility,  ONLY: lorene2hydrobase
    USE tensor,   ONLY: jxx, jxy, jxz, jyy, jyz, jzz

    IMPLICIT NONE

   ! IF ( C_ASSOCIATED( THIS% bns_ptr ) ) THEN
   !
   !   ! The coordinates need to be converted from |sphincs| units (Msun_geo)
   !   ! to |fuka| units (\(\mathrm{km}\)).
   !   ! See MODULE constants for the definition of Msun_geo
   !   CALL get_lorene_id_mass_b( THIS% bns_ptr, &
   !                                 x*Msun_geo, &
   !                                 y*Msun_geo, &
   !                                 z*Msun_geo, &
   !                                 g(jxx), &
   !                                 baryon_density, &
   !                                 gamma_euler )
   !
   !   g(jxy)= 0.0D0
   !   g(jxz)= 0.0D0
   !   g(jyy)= g(jxx)
   !   g(jyz)= 0.0D0
   !   g(jzz)= g(jxx)
   !
   !   baryon_density= baryon_density*lorene2hydrobase
   !
   ! ENDIF

  END PROCEDURE read_fuka_id_mass_b


  MODULE PROCEDURE read_fuka_id_k

    !****************************************************
    !
    !# Stores the components of the extrinsic curvature
    !  in arrays
    !
    !  FT 09.02.2022
    !
    !****************************************************

    IMPLICIT NONE

   ! IF ( C_ASSOCIATED( THIS% bns_ptr ) ) THEN
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
   !   !$OMP             SHARED( n, THIS, x, y, z, &
   !   !$OMP                     k_xx, k_xy, k_xz, k_yy, k_yz, k_zz ) &
   !   !$OMP             PRIVATE( itr )
   !   read_fuka_id_loop: DO itr= 1, n, 1
   !
   !     ! The coordinates need to be converted from |sphincs| units (Msun_geo)
   !     ! to |fuka| units (\(\mathrm{km}\)). See MODULE constants for the definition of
   !     ! Msun_geo
   !     CALL get_lorene_id_k( THIS% bns_ptr, &
   !                           x( itr )*Msun_geo, &
   !                           y( itr )*Msun_geo, &
   !                           z( itr )*Msun_geo, &
   !                           k_xx( itr ), &
   !                           k_xy( itr ), &
   !                           k_xz( itr ), &
   !                           k_yy( itr ), &
   !                           k_yz( itr ), &
   !                           k_zz( itr ) )
   !
   !   ENDDO read_fuka_id_loop
   !   !$OMP END PARALLEL DO
   !
   !   DO itr= 1, n, 1
   !
   !     !
   !     !-- Convert the extrinsic curvature from |fuka| units to
   !     !-- |sphincs| units
   !     !
   !     k_xx( itr )= k_xx( itr )*Msun_geo
   !     k_xy( itr )= k_xy( itr )*Msun_geo
   !     k_xz( itr )= k_xz( itr )*Msun_geo
   !     k_yy( itr )= k_yy( itr )*Msun_geo
   !     k_yz( itr )= k_yz( itr )*Msun_geo
   !     k_zz( itr )= k_zz( itr )*Msun_geo
   !
   !     ! Print progress on screen
   !     perc= 100*itr/n
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
    !  \(M_\odot/L_\odot^3\).
    !
    !  FT 09.02.2022
    !
    !***********************************************

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_ASSOCIATED
    USE utility,                     ONLY: lorene2hydrobase

    IMPLICIT NONE

  !  IF ( C_ASSOCIATED( THIS% bns_ptr ) )THEN
  !
  !    ! The coordinates need to be converted from |sphincs| units (Msun_geo)
  !    ! to |fuka| units (\(\mathrm{km}\)). See MODULE constants for the definition of
  !    ! Msun_geo
  !    res= get_lorene_mass_density( THIS% bns_ptr, &
  !                                  x*Msun_geo, &
  !                                  y*Msun_geo, &
  !                                  z*Msun_geo )*lorene2hydrobase
  !
  !  ENDIF

  END PROCEDURE read_fuka_mass_density


  MODULE PROCEDURE read_fuka_spatial_metric

    !***********************************************
    !
    !# Returns the |fuka| conformal factor to the
    !  4th power, equal to the diagonal components
    !  of the conformally flat spatial ADM metric.
    !
    !  FT 09.02.2022
    !
    !***********************************************

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_ASSOCIATED

    IMPLICIT NONE

  !  IF ( C_ASSOCIATED( THIS% bns_ptr ) )THEN
  !
  !    ! The coordinates need to be converted from |sphincs| units (Msun_geo)
  !    ! to |fuka| units (\(\mathrm{km}\)). See MODULE constants for the definition of
  !    ! Msun_geo
  !    res= get_lorene_spatial_metric( THIS% bns_ptr, &
  !                                    x*Msun_geo, &
  !                                    y*Msun_geo, &
  !                                    z*Msun_geo )
  !
  !  ENDIF

  END PROCEDURE read_fuka_spatial_metric


  MODULE PROCEDURE is_hydro_positive

    !************************************************
    !
    !# Return 1 if the energy density is nonpositive
    !  or if the specific energy is nonpositive,
    !  or if the pressure is nonpositive
    !  at the specified point
    !
    !  FT 09.02.2022
    !
    !************************************************

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_ASSOCIATED

    IMPLICIT NONE

  !  IF ( C_ASSOCIATED( THIS% bns_ptr ) )THEN
  !
  !    ! The coordinates need to be converted from |sphincs| units (Msun_geo)
  !    ! to |fuka| units (\(\mathrm{km}\)). See MODULE constants for the definition of
  !    ! Msun_geo
  !    res= negative_hydro( THIS% bns_ptr, &
  !                                  x*Msun_geo, &
  !                                  y*Msun_geo, &
  !                                  z*Msun_geo )
  !
  !  ENDIF

  END PROCEDURE is_hydro_positive


END SUBMODULE read
