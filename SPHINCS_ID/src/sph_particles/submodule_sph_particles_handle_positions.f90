! File:         submodule_sph_particles_handle_positions.f90
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

    USE utility,  ONLY: zero

    IMPLICIT NONE

    INTEGER:: a, npart_diff

    INTEGER, DIMENSION(npart):: above_xy_plane, below_xy_plane

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: &
      above_xy_plane_a_tmp, below_xy_plane_a_tmp

    above_xy_plane= zero
    below_xy_plane= zero
    npart_above_xy= zero
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( pos, above_xy_plane, below_xy_plane, npart ) &
    !$OMP             PRIVATE( a ) &
    !$OMP             REDUCTION(+:npart_above_xy)
    DO a= 1, npart, 1

      IF( pos(3,a) > zero )THEN

        npart_above_xy= npart_above_xy + 1
        above_xy_plane(a)= a

      ELSEIF( pos(3,a) < zero )THEN

        below_xy_plane(a)= a

      ENDIF

    ENDDO
    !$OMP END PARALLEL DO

    npart_diff= 0
    IF( npart/2 /= npart_above_xy )THEN

      PRINT *, "** WARNING! Mismatch in the number of particles above the xy ",&
               "plane in SUBROUTINE reflect_particles_xy_plane!"
      PRINT *, " * npart/2= ", npart/2
      PRINT *, " * npart_above_xy= ", npart_above_xy
      PRINT *, " * If you are inside the APM iteration, you're safe since ", &
               "this is taken care of."
      PRINT *, "   Otherwise, you may want to double check that ", &
               "you know what's going on."
      PRINT *

      IF( npart/2 > npart_above_xy )THEN
        npart_diff= npart/2 - npart_above_xy
      ENDIF

      npart_above_xy= npart/2

    ENDIF

    ALLOCATE(above_xy_plane_a(npart_above_xy))

    above_xy_plane_a_tmp= PACK( above_xy_plane, above_xy_plane /= zero )
    below_xy_plane_a_tmp= PACK( below_xy_plane, below_xy_plane /= zero )

    IF( npart_diff > 0 )THEN
      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( pos, above_xy_plane_a, above_xy_plane_a_tmp, npart ) &
      !$OMP             PRIVATE( a )
      DO a= 1, SIZE(above_xy_plane_a_tmp), 1

         above_xy_plane_a(a)= above_xy_plane_a_tmp(a)

      ENDDO
      !$OMP END PARALLEL DO
      IF( npart_above_xy /= SIZE(above_xy_plane_a_tmp) + npart_diff )THEN
         PRINT *, "** ERROR! Mismatch in the number of particles above the xy ",&
                  "plane in SUBROUTINE reflect_particles_xy_plane!"
         PRINT *, " * SIZE(above_xy_plane_a_tmp)= ", SIZE(above_xy_plane_a_tmp)
         PRINT *, " * npart_diff= ", npart_diff
         PRINT *, " * npart_above_xy= ", npart_above_xy
         PRINT *, " * Stopping..."
         PRINT *
         STOP
      ENDIF
      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( pos, above_xy_plane_a, above_xy_plane_a_tmp, &
      !$OMP                     npart_above_xy, below_xy_plane_a_tmp ) &
      !$OMP             PRIVATE( a )
      DO a= SIZE(above_xy_plane_a_tmp) + 1, npart_above_xy, 1

         above_xy_plane_a(a)= below_xy_plane_a_tmp(a)

      ENDDO
      !$OMP END PARALLEL DO
    ELSE
      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( pos, above_xy_plane_a, above_xy_plane_a_tmp, npart_above_xy ) &
      !$OMP             PRIVATE( a )
      DO a= 1, npart_above_xy, 1

         above_xy_plane_a(a)= above_xy_plane_a_tmp(a)

      ENDDO
      !$OMP END PARALLEL DO
    ENDIF

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

    !CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile
    !LOGICAL:: exist

    IF( PRESENT(nu) .NEQV. PRESENT(nu_below) )THEN
      PRINT *, "** ERROR! In SUBROUTINE reflect_particles_xy_plane, the ", &
               "arguments 'nu' and 'nu_below must be either both present or ", &
               "both absent."
      PRINT *, " * Stopping..."
      PRINT *
      STOP
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

  !  DOUBLE PRECISION, DIMENSION(3,npart+npart_ghost):: pos_tmp
  !  DOUBLE PRECISION, DIMENSION(npart+npart_ghost)  :: nu_tmp

    IF( MOD(npart,2) /= 0 )THEN
      PRINT *, "** ERROR! If the equatorial symmetry has to be imposed, ", &
               "the particle number must be even."
      PRINT *, " * npart= ", npart
      PRINT *, " * Stopping..."
      PRINT *
      STOP
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

    IF( PRESENT(verbose) )THEN
      IF( verbose .EQV. .TRUE. )THEN

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
    ENDIF

  END PROCEDURE impose_equatorial_plane_symmetry


  MODULE PROCEDURE check_particle_positions

    !*************************************************
    !
    !# Check that the particles are not at the same
    !  positions
    !
    !  FT 1.9.2021
    !
    !*************************************************

    USE NR, ONLY: indexx

    IMPLICIT NONE

    INTEGER:: itr
    !! Iterator
    INTEGER:: itr2
    !! Iterator
    INTEGER:: x_idx
    !# Index at which a new value of the \(x\) coordinate appears,
    !  in the array `pos` sorted so that the \(x\) coordinate does not decrease
    INTEGER, DIMENSION(:), ALLOCATABLE:: x_sort
    !# Array storing the sorted indices of array `pos`, so that the \(x\)
    !  coordinate of the particles is in nondecreasing order
    INTEGER, DIMENSION(:), ALLOCATABLE:: x_number
    !# Array storing, for each \(x\) coordinate, the number of particles
    !  having that \(x\) coordinate

    PRINT *, "** Checking that there are not multiple particles", &
             " at the same position..."
    PRINT *

    ALLOCATE( x_sort( npart ) )
    ALLOCATE( x_number( npart ) )

    ! Sort x coordinates of the particles
    CALL indexx( npart, pos( 1, : ), x_sort )

    x_number= 1
    itr2= 1
    ! Find the number of times each x appears
    ! TODO: parallelize this
    DO itr= 1, npart - 1, 1

      IF( pos( 1, x_sort(itr) ) == &
          pos( 1, x_sort(itr+1) ) )THEN

        x_number(itr2)= x_number(itr2) + 1

      ELSE

        itr2= itr2 + 1

      ENDIF

    ENDDO
    x_number= x_number(1:itr2)

    IF( SUM( x_number ) /= npart )THEN

      PRINT *, "** ERROR! The sum of the numbers of particles with the same", &
               " x is not equal to the particle number."
      PRINT *, " * SUM( x_number )=", SUM( x_number ), ", ", &
               "npart=", npart
      PRINT *, " * Stopping..."
      PRINT *
      STOP

    ENDIF

    IF( PRESENT(debug) )THEN

      IF( debug .EQV. .TRUE. )THEN
      ! These two IF statements are nested because gfortran doesn't like them
      ! together when debug is not PRESENT

        !$OMP PARALLEL DO DEFAULT( NONE ) &
        !$OMP             SHARED( pos, x_sort, x_number ) &
        !$OMP             PRIVATE( itr, itr2, x_idx )
        DO itr= 1, SIZE(x_number), 1

          IF( itr == 1 )THEN
            x_idx= 1
          ELSE
            x_idx= SUM(x_number(1:itr-1)) + 1
          ENDIF

          DO itr2= x_idx, x_idx + x_number(itr) - 2, 1

            ! If they do not have the same x
            IF( pos( 1, x_sort(itr2) ) /= &
                pos( 1, x_sort(itr2+1) ) )THEN

              PRINT *, "** ERROR! ", "The two particles ", x_sort(itr2), &
                       " and", x_sort(itr2+1), &
                       " do not have the same x, but should!"
              PRINT *, pos( :, x_sort(itr2) )
              PRINT *, pos( :, x_sort(itr2+1) )
              PRINT *, " * Stopping..."
              PRINT *
              STOP

            ENDIF

          ENDDO
        ENDDO
        !$OMP END PARALLEL DO

      ENDIF

    ENDIF

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( pos, x_sort, x_number ) &
    !$OMP             PRIVATE( itr, itr2, x_idx )
    DO itr= 1, SIZE(x_number), 1

      IF( itr == 1 )THEN
        x_idx= 1
      ELSE
        x_idx= SUM(x_number(1:itr-1)) + 1
      ENDIF

      DO itr2= x_idx, x_idx + x_number(itr) - 2, 1

        ! If they have the same y
        IF( pos( 2, x_sort(itr2) ) == &
            pos( 2, x_sort(itr2+1) ) )THEN

          ! If they have the same z
          IF( pos( 3, x_sort(itr2) ) == &
              pos( 3, x_sort(itr2+1) ) )THEN

            ! They are the same
            PRINT *, "** ERROR! ", "The two particles ", x_sort(itr2), &
                     " and", x_sort(itr2+1), " have the same coordinates!"
            PRINT *, pos( :, x_sort(itr2) )
            PRINT *, pos( :, x_sort(itr2+1) )
            PRINT *, " * Stopping..."
            PRINT *
            STOP

          ENDIF
        ENDIF

      ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    DEALLOCATE( x_sort )
    DEALLOCATE( x_number )


  END PROCEDURE check_particle_positions


  MODULE PROCEDURE check_particle_position

    !*****************************************************
    !
    !# Return the number of times that pos_a appears
    !  in the array pos
    !  @todo This algorithm scales as O(npart**2)
    !        if used in a loop over the particles...
    !        To be documented, after it's fixed
    !
    !  FT 13.10.2021
    !
    !*****************************************************

    !USE NR,             ONLY: indexx

    IMPLICIT NONE


    INTEGER:: itr, itr2, size_x!, cnt
    INTEGER, DIMENSION(npart):: x_sort, cnts
    INTEGER, DIMENSION(npart):: x_number
    INTEGER, DIMENSION(:), ALLOCATABLE:: x_number_filt

    ! Sort x coordinates of the particles
    !CALL indexx( npart, pos( 1, : ), x_sort )

    x_number= 0
    itr2= 0
    ! Find the number of times that the x coordinate of pos_a appears in pos
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( pos, pos_a, x_sort, x_number, npart ) &
    !$OMP             PRIVATE( itr )
    DO itr= 1, npart, 1

      IF( pos( 1, itr ) == pos_a( 1 ) )THEN

        !itr2= itr2 + 1
        x_number(itr)= itr

      !ELSEIF( pos( 1, x_sort(itr) ) > pos_a( 1 ) )THEN
      !
      !  EXIT

      ENDIF

    ENDDO
    !$OMP END PARALLEL DO
    x_number_filt= PACK( x_number, x_number /= 0 )
    size_x= SIZE(x_number_filt)

    cnts= 0
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( pos, pos_a, x_sort, x_number_filt, size_x, cnts )&
    !$OMP             PRIVATE( itr )
    DO itr= 1, size_x, 1

      ! If they have the same y
      IF( pos( 2, x_number_filt(itr) ) == pos_a( 2 ) )THEN

        ! If they have the same z
        IF( pos( 3, x_number_filt(itr) ) == pos_a( 3 ) )THEN

          cnts(itr)= cnts(itr) + 1

        ENDIF
      ENDIF

    ENDDO
    !$OMP END PARALLEL DO

    cnt= SUM( cnts )

  END PROCEDURE check_particle_position


  MODULE PROCEDURE correct_center_of_mass

    !***********************************************************
    !
    !# Translate the particles so that their center of mass
    !  coincides with the center of mass of the star, given by
    !  |id|
    !
    !  FT 1.09.2021
    !
    !***********************************************************

    USE analyze,    ONLY: COM
    USE utility,    ONLY: is_finite_number, zero

    IMPLICIT NONE

    INTEGER:: a, i
    DOUBLE PRECISION:: com_x, com_y, com_z, com_d
    DOUBLE PRECISION, DIMENSION(3):: pos_corr_tmp

    CALL COM( npart, pos, nu, &            ! input
              com_x, com_y, com_z, com_d ) ! output

    IF( PRESENT(verbose) )THEN
      IF( verbose .EQV. .TRUE. )THEN
        PRINT *, "** Before center of mass correction:"
        PRINT *, " * x coordinate of the center of mass of the star, ", &
                 "from the ID: com_star= ", com_star(1), "Msun_geo"
        PRINT *, " * y coordinate of the center of mass of the star, ", &
                 "from the ID: com_star= ", com_star(2), "Msun_geo"
        PRINT *, " * z coordinate of the center of mass of the star, ", &
                 "from the ID: com_star= ", com_star(3), "Msun_geo"
        PRINT *, " * x coordinate of the center of mass of the particle ", &
                 "distribution: com_x= ", com_x, "Msun_geo"
        PRINT *, " * y coordinate of the center of mass of the particle ", &
                 "distribution: com_y= ", com_y, "Msun_geo"
        PRINT *, " * z coordinate of the center of mass of the particle ", &
                 "distribution: com_z= ", com_z, "Msun_geo"
        PRINT *, " * Distance of the center of mass of the particle ", &
                 "distribution from the  origin: com_d= ", com_d
        PRINT *, " * |com_x-com_star_x/com_star_x|=", &
                 ABS( com_x-com_star(1) )/ABS( com_star(1) + 1 )
        PRINT *
      ENDIF
    ENDIF

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( pos, com_star, &
    !$OMP                     com_x, com_y, com_z, npart ) &
    !$OMP             PRIVATE( pos_corr_tmp, a )
    DO a= 1, npart, 1

      pos_corr_tmp(1)= pos(1,a) - ( com_x - com_star(1) )
      pos_corr_tmp(2)= pos(2,a) - ( com_y - com_star(2) )
      pos_corr_tmp(3)= pos(3,a) - ( com_z - com_star(3) )

      DO i= 1, 3, 1
        IF( .NOT.is_finite_number(pos_corr_tmp(i)) )THEN
          PRINT *, pos_corr_tmp
          STOP
        ENDIF
      ENDDO

      IF( get_density( &
                  pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) > zero &
          .AND. &
          validate_pos( &
                  pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) &
      )THEN

        pos(:,a)= pos_corr_tmp

      ENDIF

    ENDDO
    !$OMP END PARALLEL DO

    CALL COM( npart, pos, nu, & ! input
              com_x, com_y, com_z, com_d ) ! output

    IF( PRESENT(verbose) )THEN
      IF( verbose .EQV. .TRUE. )THEN
        PRINT *, "** After center of mass correction:"
        PRINT *, " * x coordinate of the center of mass of the star, ", &
                 "from the ID: com_star= ", com_star(1), "Msun_geo"
        PRINT *, " * y coordinate of the center of mass of the star, ", &
                 "from the ID: com_star= ", com_star(2), "Msun_geo"
        PRINT *, " * z coordinate of the center of mass of the star, ", &
                 "from the ID: com_star= ", com_star(3), "Msun_geo"
        PRINT *, " * x coordinate of the center of mass of the particle ", &
                 "distribution: com_x= ", com_x, "Msun_geo"
        PRINT *, " * y coordinate of the center of mass of the particle ", &
                 "distribution: com_y= ", com_y, "Msun_geo"
        PRINT *, " * z coordinate of the center of mass of the particle ", &
                 "distribution: com_z= ", com_z, "Msun_geo"
        PRINT *, " * Distance of the center of mass of the particle ", &
                 "distribution from the  origin: com_d= ", com_d
        PRINT *, " * |com_x-com_star_x/com_star_x|=", &
                 ABS( com_x-com_star(1) )/ABS( com_star(1) + 1 )
        PRINT *
      ENDIF
    ENDIF

  END PROCEDURE correct_center_of_mass


  MODULE PROCEDURE get_neighbours_bf

    !**************************************************************
    !
    !# just for test purposes: get neighbours of particle ipart in
    !  a "brute force" way; ipart is ALSO on the neighbour list;
    !  SKR 8.2.2010
    !
    !  Removed ipart from its own neighbors' list
    !  FT 04.06.2021
    !
    !**************************************************************

    IMPLICIT NONE

    INTEGER:: a
    DOUBLE PRECISION:: diff(dimensions), d2, r_int2

    ! square of interaction radius
    r_int2= (2.D0*h(ipart))**2

    nnei= 0
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP             SHARED(pos,dimensions,ipart,npart,r_int2,nnei,neilist)&
    !$OMP             PRIVATE(a,diff,d2)
    DO a= 1, npart, 1

      IF( a /= ipart )THEN

        diff= pos(1:dimensions,a)-pos(1:dimensions,ipart)
        d2= DOT_PRODUCT(diff,diff)

        ! neighbour?
        IF(d2 < r_int2)THEN
          nnei= nnei + 1
          neilist(nnei)= a
        ENDIF

      ENDIF

    ENDDO
    !$OMP END PARALLEL DO

  END PROCEDURE get_neighbours_bf


  MODULE PROCEDURE get_object_of_particle

    !**************************************************************
    !
    !# Returns the number of the matter object asociated with the particle
    !  number given as input. Example: give number \(n\) as input; this
    !  particle number corresponds to a particle on matter object \(m\).
    !  This functions returns \(m\).
    !
    !  FT 05.12.2022
    !
    !**************************************************************

    IMPLICIT NONE

    INTEGER:: i_matter

    DO i_matter= 1, this% n_matter, 1

      IF(a > this% npart_fin(i_matter - 1) &
         .AND. a <= this% npart_fin(i_matter))THEN

        get_object_of_particle= i_matter
        RETURN

      ENDIF

    ENDDO

  END PROCEDURE get_object_of_particle


END SUBMODULE handle_positions
