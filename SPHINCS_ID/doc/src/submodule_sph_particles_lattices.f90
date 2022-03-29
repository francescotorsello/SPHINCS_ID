! File:         submodule_sph_particles_lattices.f90
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

SUBMODULE (sph_particles) lattices

  !***************************************************
  !
  !# This SUBMODULE contains the implementation of
  !  the methods of TYPE sph_particles
  !  that place particles on a lattice within a
  !  matter object.
  !
  !  FT 12.07.2021
  !
  !***************************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE place_particles_lattice

    !*********************************************************
    !
    !# Places paricles on a 3D lattice containing both stars
    !
    !  FT 5.10.2020
    !
    !*********************************************************

    USE constants,    ONLY: pi, third

    IMPLICIT NONE

    INTEGER:: i, j, k, sgn, nx, ny, nz, npart_half
    INTEGER:: npart_tmp

    DOUBLE PRECISION:: dx, dy, dz, vol, vol_a
    DOUBLE PRECISION:: xtemp, ytemp, ztemp, zlim
    DOUBLE PRECISION:: thres_baryon_density
    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: pos_tmp

    PRINT *, "** Executing the place_particles_lattice " &
             // "subroutine..."
    PRINT *

    !
    !-- Set the boundary in z
    !
    IF( ABS(zmax) > ABS(zmin) )THEN
      zlim= zmax
    ELSE
      zlim= zmin
    ENDIF

    !
    !-- Compute number of lattice points (for now, equal in each direction)
    !
    nx= FLOOR(stretch*(6.0D0*DBLE(npart_des)/pi)**third)
    IF( MOD( nx, 2 ) /= 0 ) nx= nx + 1
    ny= nx
    nz= nx

    !
    !-- Consistency checks
    !
    IF( MOD( nz, 2 ) /= 0 )THEN
      PRINT *
      PRINT *, "** ERROR! nz should be even!"
      PRINT *
      STOP
    ENDIF

    IF( nx == 0 .OR. ny == 0 .OR. nz == 0 )THEN
      PRINT *
      PRINT *, "** ERROR! nx, ny, nz are 0!"
      PRINT *
      STOP
    ENDIF

    PRINT *, " * nx= ny= nz=", nx
    PRINT *

    !
    !-- Compute lattice steps
    !
    dx= ABS(xmax - xmin)/DBLE( nx )
    dy= ABS(ymax - ymin)/DBLE( ny )
    dz= ABS(zlim)/DBLE( nz/2 )

    PRINT *, " * dx=", dx,  ", dy=", dx,  ", dz=", dz
    PRINT *

    npart_tmp = nx*ny*nz

    PRINT *, " * Number of lattice points= nx*ny*nz=", npart_tmp
    PRINT *

    !
    !-- Set the threshold above which a lattice point is
    !-- promoted to a particle
    !
    IF( THIS% use_thres )THEN
      thres_baryon_density= central_density/thres
    ELSE
      thres_baryon_density= 0.0D0
    ENDIF

    IF(.NOT.ALLOCATED( pos_tmp ))THEN
      ALLOCATE( pos_tmp( 3, nx, ny, nz ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array pos_tmp in SUBROUTINE" &
                  // "place_particles_3D_lattice. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array pos in SUBROUTINE" &
      !                // "place_particles_3D_lattice." )
    ENDIF
    ! Initializing the array pos to 0
    pos_tmp= HUGE(0.0D0)

    !---------------------------------------------------------!
    !--  Storing the particle positions into the array pos  --!
    !--  symmetrically w.r.t. the xy plane                  --!
    !---------------------------------------------------------!

    PRINT *, " * Placing particles on the lattice..."
    PRINT *

    !
    !-- Choose the larger value for the boundary in z
    !
    IF( zlim == zmin )THEN
      sgn= - 1
    ELSE
      sgn= 1
    ENDIF

    !
    !-- Place the first half of the particle (above or below the xy plane,
    !-- depending on the variable sgn)
    !
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( nx, ny, nz, dx, dy, dz, sgn, &
    !$OMP                     pos_tmp, thres_baryon_density, xmin, ymin ) &
    !$OMP             PRIVATE( i, j, k, xtemp, ytemp, ztemp )
    particle_pos_z: DO k= 1, nz/2, 1

      ztemp= sgn*( dz/2 + ( k - 1 )*dz )

      particle_pos_y: DO j= 1, ny, 1

        ytemp= ymin + ( j - 1 )*dy

        particle_pos_x: DO i= 1, nx, 1

          xtemp= xmin + dx/2 + ( i - 1 )*dx

          !
          !-- Promote a lattice point to a particle,
          !-- if the mass density is higher than the threshold
          !
          IF( get_density( xtemp, ytemp, ztemp ) &
                                  > thres_baryon_density &
              .AND. &
              validate_position( xtemp, ytemp, ztemp ) == 0 )THEN

            pos_tmp( 1, i, j, k )= xtemp
            pos_tmp( 2, i, j, k )= ytemp
            pos_tmp( 3, i, j, k )= ztemp

          ENDIF

         ENDDO particle_pos_x
      ENDDO particle_pos_y
    ENDDO particle_pos_z
    !$OMP END PARALLEL DO

    npart_out= 0
    DO k= 1, nz, 1

      DO j= 1, ny, 1

        DO i= 1, nx, 1

          IF( pos_tmp( 1, i, j, k ) < HUGE(0.0D0) )THEN

            npart_out= npart_out + 1

            pos( 1, npart_out )= pos_tmp( 1, i, j, k )
            pos( 2, npart_out )= pos_tmp( 2, i, j, k )
            pos( 3, npart_out )= pos_tmp( 3, i, j, k )

          ENDIF

         ENDDO
      ENDDO
    ENDDO
    npart_half= npart_out
    IF( npart_half == 0 )THEN
      PRINT *, "** There are no particles! Execution stopped..."
      PRINT *
      STOP
    ENDIF

    DEALLOCATE( pos_tmp )

    !
    !-- Place the second half of the particles, mirroring the first half
    !-- w.r.t the xy plane
    !
    particle_pos_z_mirror: DO k= 1, npart_half, 1

      xtemp=   pos( 1, k )
      ytemp=   pos( 2, k )
      ztemp= - pos( 3, k )

      npart_out= npart_out + 1
      pos( 1, npart_out )= xtemp
      pos( 2, npart_out )= ytemp
      pos( 3, npart_out )= ztemp

      !ENDIF

      ! Print progress on screen, every 10%
      !perc= 50 + 50*k/( npart_half )
      !IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
      !   WRITE( *, "(A2,I3,A1)", ADVANCE= "NO" ) &
      !           creturn//" ", perc, "%"
      !ENDIF
    ENDDO particle_pos_z_mirror
    !WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

    !
    !-- Consistency checks
    !
    IF( npart_out /= 2*npart_half )THEN
      PRINT *
      PRINT *, "** ERROR: The number of particles ", npart_out, &
               " is not the expected value ", 2*npart_half
      PRINT *
      STOP
    ENDIF

    DO k= 1, npart_half, 1
      IF( pos( 3, k ) /= - pos( 3, npart_half + k ) )THEN
        PRINT *
        PRINT *, "** ERROR: The lattice is not mirrored " &
                 // "by the xy plane."
        PRINT *
        STOP
      ENDIF
    ENDDO

    PRINT *, " * Particles placed. Number of particles=", &
             npart_out, "=", DBLE(npart_out)/DBLE(npart_tmp), &
             " of the points in lattice."
    PRINT *

    !
    !-- Computing total volume and volume per particle
    !
    IF(.NOT.ALLOCATED( pvol ))THEN
      ALLOCATE( pvol( npart_out ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array pvol ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for array v_euler_parts_z" )
    ENDIF

    vol  = (xmax - xmin)*(ymax - ymin)*2*ABS(zlim)
    vol_a= vol/npart_tmp

    pvol= vol_a

    ! Consistency check for the particle volume
    IF( ABS( vol_a - dx*dy*dz ) > 1.0D-9 )THEN
      PRINT *, " * The particle volume vol_a=", vol_a, "Msun_geo^3"
      PRINT *, " is not equal to dx*dy*dz=", dx*dy*dz, "Msun_geo^3."
      PRINT *
      STOP
    ENDIF

    PRINT *, " * Total volume of the lattices=", vol, "Msun_geo^3"
    PRINT *, " * Particle volume=", vol_a, "Msun_geo^3"
    PRINT *

    PRINT *, "** Subroutine place_particles_3D_lattice executed."
    PRINT *

  END PROCEDURE place_particles_lattice


END SUBMODULE lattices
