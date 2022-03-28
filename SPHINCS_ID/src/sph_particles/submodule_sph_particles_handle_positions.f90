! File:         submodule_sph_particles_handle_positions.f90
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

SUBMODULE (sph_particles) handle_positions

  !**************************************************
  !
  !# This SUBMODULE contains the implementation of
  !  the PROCEDURES to handle particle positions.
  !
  !  FT 24.03.2022
  !
  !**************************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE find_particles_above_xy_plane

    !*************************************************************
    !
    !# Find the particles above the \(xy\) plane
    !
    !  FT 25.03.2022
    !
    !*************************************************************

    USE constants,  ONLY: zero

    IMPLICIT NONE

    INTEGER:: a

    INTEGER, DIMENSION(npart):: above_xy_plane

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile
    LOGICAL:: exist

    above_xy_plane= zero
    npart_above_xy= zero
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( pos, above_xy_plane, npart ) &
    !$OMP             PRIVATE( a ) &
    !$OMP             REDUCTION(+:npart_above_xy)
    DO a= 1, npart, 1

      IF( pos(3,a) > zero )THEN

        npart_above_xy= npart_above_xy + 1
        above_xy_plane(a)= a

      ENDIF

    ENDDO
    !$OMP END PARALLEL DO

    ALLOCATE(above_xy_plane_a(npart_above_xy))
    above_xy_plane_a= PACK( above_xy_plane, above_xy_plane /= zero )

  END PROCEDURE find_particles_above_xy_plane


  MODULE PROCEDURE reflect_particles_xy_plane

    !*************************************************************
    !
    !# Reflect the particle with z>0 with respect to the xy plane
    !
    !  FT 25.03.2022
    !
    !*************************************************************

    IMPLICIT NONE

    INTEGER:: a

    DOUBLE PRECISION, DIMENSION(3,npart):: pos_tmp
    DOUBLE PRECISION, DIMENSION(npart)  :: nu_tmp

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile
    LOGICAL:: exist

    IF( PRESENT(nu) .NEQV. PRESENT(nu_below) )THEN
      PRINT *, "** ERROR! In SUBROUTINE reflect_particles_xy_plane, the ", &
               "arguments 'nu' and 'nu_below must be either both present or ", &
               "both absent."
      PRINT *, " * Stopping..."
      PRINT *
      STOP
    ENDIF

    IF( npart/2 /= npart_above_xy )THEN

      PRINT *, "** WARNING! Mismatch in the number of particles above the xy ",&
               "plane in SUBROUTINE reflect_particles_xy_plane!"
      PRINT *, " * npart/2= ", npart/2
      PRINT *, " * npart_above_xy= ", npart_above_xy
      PRINT *, " * If you are inside the APM, you're safe since this is ", &
               "taken care of. Otherwise, you may want to double check that ", &
               "you know what's going on."
      PRINT *

    ENDIF

    pos_tmp= pos
    IF( PRESENT(nu) ) nu_tmp = nu

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( pos, pos_tmp, above_xy_plane_a, npart_above_xy, &
    !$OMP                     nu, nu_tmp ) &
    !$OMP             PRIVATE( a )
    DO a= 1, npart_above_xy, 1

      pos( :, a )= pos_tmp( :, above_xy_plane_a(a) )
      IF( PRESENT(nu) ) nu( a )= nu_tmp( above_xy_plane_a(a) )

    ENDDO
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( pos, pos_below, npart_above_xy, nu, nu_below ) &
    !$OMP             PRIVATE( a )
    DO a= 1, npart_above_xy, 1
      pos_below( 1, a )=   pos( 1, a )
      pos_below( 2, a )=   pos( 2, a )
      pos_below( 3, a )= - pos( 3, a )
      IF( PRESENT(nu) ) nu_below( a )= nu( a )
    ENDDO
    !$OMP END PARALLEL DO

  END PROCEDURE reflect_particles_xy_plane


  MODULE PROCEDURE impose_equatorial_plane_symmetry

    !*************************************************************
    !
    !# Mirror the particle with z>0 with respect to the xy plane,
    !  to impose the equatorial-plane symmetry
    !
    !  FT 1.09.2021
    !
    !*************************************************************

    USE analyze,    ONLY: COM

    IMPLICIT NONE

    INTEGER:: a, npart_above_xy
    DOUBLE PRECISION:: com_x, com_y, com_z, com_d

    INTEGER, DIMENSION(:), ALLOCATABLE:: above_xy_plane_a

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos_below
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nu_below

  !  DOUBLE PRECISION, DIMENSION(3,npart_real+npart_ghost):: pos_tmp
  !  DOUBLE PRECISION, DIMENSION(npart_real+npart_ghost)  :: nu_tmp

    IF( MOD(npart,2) /= 0 )THEN
      PRINT *, "** ERROR! If the equatorial symmetry has to be imposed, ", &
               "the particle number must be even."
      PRINT *, " * npart= ", npart
      PRINT *, " * Stopping..."
      PRINT *
    ENDIF

  !  IF( PRESENT(pos_prev) )THEN
  !
  !    ! If some of the particles crossed the xy plane top-down in the
  !    ! last step, replace their z coordinate with their previous
  !    ! z coordinate
  !    !$OMP PARALLEL DO DEFAULT( NONE ) &
  !    !$OMP             SHARED( pos, pos_prev, npart ) &
  !    !$OMP             PRIVATE( a )
  !    DO a= 1, npart, 1
  !
  !      IF( (pos_prev(3,a) > zero) .AND. (pos(3,a) <= zero) )THEN
  !
  !        pos(3,a)= pos_prev(3,a)
  !
  !      ENDIF
  !
  !    ENDDO
  !    !$OMP END PARALLEL DO
  !
  !  ENDIF

  !  pos_tmp= pos
  !  nu_tmp = nu

   ! itr= 0
   ! DO a= 1, npart, 1
   !
   !   IF( pos_tmp( 3, a ) > zero &
   !       .AND. &
   !       itr <= npart/2 )THEN
   !
   !     itr= itr + 1
   !     pos( 1, itr )= pos_tmp( 1, a )
   !     pos( 2, itr )= pos_tmp( 2, a )
   !     pos( 3, itr )= pos_tmp( 3, a )
   !     IF( PRESENT(nu) ) nu( itr )= nu_tmp( a )
   !
   !   ENDIF
   !
   ! ENDDO
   ! npart_above_xy= itr

    CALL find_particles_above_xy_plane( npart, pos, npart_above_xy, &
                                        above_xy_plane_a )

   ! IF(npart/2 /= npart_above_xy )THEN
   !
   !
   !
   ! ENDIF
    ALLOCATE( pos_below(3,npart_above_xy) )

    IF( PRESENT(nu) )THEN

      ALLOCATE( nu_below(npart_above_xy) )
      CALL reflect_particles_xy_plane( npart, pos, pos_below, npart_above_xy, &
                                       above_xy_plane_a, nu, nu_below )

    ELSE

      CALL reflect_particles_xy_plane( npart, pos, pos_below, npart_above_xy, &
                                       above_xy_plane_a )

    ENDIF

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( pos, npart_above_xy, nu, &
    !$OMP                     pos_below, nu_below ) &
    !$OMP             PRIVATE( a )
    DO a= 1, npart_above_xy, 1
      pos( :, npart_above_xy + a )= pos_below( :, a )
      IF( PRESENT(nu) )THEN
        nu( npart_above_xy + a )= nu_below( a )
      ENDIF
    ENDDO
    !$OMP END PARALLEL DO

    IF( PRESENT(verbose) .AND. verbose .EQV. .TRUE. )THEN

      CALL COM( npart, pos, nu, & ! input
                com_x, com_y, com_z, com_d ) ! output

      PRINT *, "** After mirroring particles:"
      IF( PRESENT(com_star) ) PRINT *, &
               " * x coordinate of the center of mass of the star, ", &
               "from LORENE: com_star= ", com_star, "Msun_geo"
      PRINT *, " * x coordinate of the center of mass of the particle ", &
               "distribution: com_x= ", com_x, "Msun_geo"
      PRINT *, " * y coordinate of the center of mass of the particle ", &
               "distribution: com_y= ", com_y, "Msun_geo"
      PRINT *, " * z coordinate of the center of mass of the particle ", &
               "distribution: com_z= ", com_z, "Msun_geo"
      PRINT *, " * Distance of the center of mass of the particle ", &
               "distribution from the  origin: com_d= ", com_d
      IF( PRESENT(com_star) ) PRINT *, " * |com_x-com_star/com_star|=", &
               ABS( com_x-com_star )/ABS( com_star + 1 )
      PRINT *

    ENDIF

  END PROCEDURE impose_equatorial_plane_symmetry


END SUBMODULE handle_positions
