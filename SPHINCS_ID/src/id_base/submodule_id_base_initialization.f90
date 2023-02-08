! File:         submodule_id_base_initialization.f90
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

SUBMODULE (id_base) initialization

  !********************************************
  !
  !# Implementation of the methods of TYPE
  !  [[idbase]] that initialize objects of
  !  [[idbase]]-extended TYPE
  !
  !  FT 8.11.2021
  !
  !********************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE sanity_check

    !************************************************
    !
    !# Checks that [[idbase:n_matter]] and the sizes
    !  returned by [[idbase:return_spatial_extent]] and
    !  [[idbase:get_total_spatial_extent]]
    !  are acceptable. It is called by initialize,
    !  after the constructor of the derived type.
    !
    !  FT 8.11.2021
    !
    !************************************************

    IMPLICIT NONE

    INTEGER:: i_matter, itr
    DOUBLE PRECISION, DIMENSION(derived_type% n_matter,6):: sizes
    DOUBLE PRECISION, DIMENSION(6):: total_sizes
    DOUBLE PRECISION, DIMENSION(derived_type% n_matter,3):: centers


    ! Check that n_matter is strictly positive
    IF( derived_type% n_matter <= 0 )THEN

      PRINT *, "** ERROR! n_matter, the number of matter objects in the", &
               " physical system, is nonpositive: n_matter= ", &
               derived_type% n_matter
      PRINT *, "   Please assign the appropriate strictly positive value", &
               " in the constructor of the TYPE that extends idbase."
      PRINT *, " * Stopping..."
      PRINT *
      STOP

    ENDIF


    ! Check that the sizes of each matter object are strictly positive
    DO i_matter= 1, derived_type% n_matter, 1

      sizes(i_matter,:)  = derived_type% return_spatial_extent(i_matter)
      centers(i_matter,:)= derived_type% return_center(i_matter)

      check_for_negative_size: DO itr= 1, 6, 1

        IF( sizes(i_matter,itr) <= 0 )THEN

          PRINT *, "** ERROR! The size n., ", itr, " of matter object ", &
                   i_matter, " is nonpositive: sizes(", itr, ")=", &
                   sizes(i_matter,itr)
          PRINT *, "   Please assign the appropriate strictly positive value", &
                   " in the constructor of the TYPE that extends idbase."
          PRINT *, " * Stopping..."
          PRINT *
          STOP

        ENDIF

      ENDDO check_for_negative_size

    !  check_for_consistent_centers: DO itr= 1, 3, 1
    !
    !    IF( ABS(centers(i_matter,itr)) <= sizes(i_matter,itr) &
    !        .OR. &
    !        ABS(centers(i_matter,itr)) >= sizes(i_matter,itr+1) )THEN
    !
    !      PRINT *, "** ERROR! The coordinate n., ", itr, " of the center", &
    !               " of matter object ", i_matter, &
    !               " is not bracketed by the sizes of the object! "
    !      PRINT *, " * Absolute value of the ", itr, &
    !               "coordinate of the center: ", ABS(centers(i_matter,itr))
    !      PRINT *, " * Sizes of the object in the direction ", itr, ": ", &
    !               sizes(i_matter,itr), sizes(i_matter,itr+1)
    !      PRINT *, "   Please assign the appropriate coordinate to the", &
    !               " center and the sizes", &
    !               " in the constructor of the TYPE that extends idbase."
    !      PRINT *, " * Stopping..."
    !      PRINT *
    !      STOP
    !
    !    ENDIF
    !
    !  ENDDO check_for_consistent_centers

    ENDDO

    ! Check that the sizes of the physical system are strictly positive
    total_sizes= derived_type% get_total_spatial_extent()

    DO itr= 1, 6, 2

      IF( total_sizes(itr+1) <= total_sizes(itr) )THEN

        PRINT *, "** ERROR! The size n. ", itr, " of the physical system ", &
                 " is larger than the size n.", itr+1
        PRINT *, " * total_sizes(", itr, ")=", total_sizes(itr)
        PRINT *, " * total_sizes(", itr+1, ")=", total_sizes(itr+1)
        PRINT *, "   Please assign the appropriate strictly positive value", &
                 " in the constructor of the TYPE that extends idbase."
        PRINT *, " * Stopping..."
        PRINT *
        STOP

      ENDIF

    ENDDO


    ! Check that the sizes of the matter objects are within the sizes of the
    ! physical system

    DO itr= 1, 6, 2

      IF( MINVAL( centers(:,CEILING(DBLE(itr)/DBLE(2))) - sizes(:,itr) ) &
          < total_sizes(itr) &
          .OR. &
          MAXVAL( centers(:,CEILING(DBLE(itr+1)/DBLE(2))) + sizes(:,itr+1) ) &
                    > total_sizes(itr+1))THEN


        PRINT *, "** ERROR! A matter object", &
                 " is not contained within the given size of the", &
                 " physical system."
        PRINT *, " * 'Left' size of the matter object= ", &
              MINVAL( centers(:,CEILING(DBLE(itr)/DBLE(2))) - sizes(:,itr) )
        PRINT *, " * 'Left' size of the physical system= ", total_sizes(itr)
        PRINT *, " * 'Right' size of the matter object= ", &
              MAXVAL( centers(:,CEILING(DBLE(itr+1)/DBLE(2))) + sizes(:,itr+1) )
        PRINT *, " * 'Right' size of the physical system=", total_sizes(itr+1)
        PRINT *
        PRINT *, "   Please assign the appropriate sizes", &
                 " in the constructor of the TYPE that extends idbase."
        PRINT *, " * Stopping..."
        PRINT *
        STOP

      ENDIF

    ENDDO


  END PROCEDURE sanity_check


  MODULE PROCEDURE initialize

    !************************************************
    !
    !# This PROCEDURE calls the constructor of the
    !  [[idbase]]-extended type and the SUBROUTINE
    !  [[idbase:sanity_check]] afterwards. It is recommended
    !  to use this SUBROUTINE to construct objects of
    !  [[idbase]]-extended type since the sanity check is
    !  performed automatically.
    !
    !  FT 8.11.2021
    !
    !************************************************


    IMPLICIT NONE


    CALL derived_type% derived_type_constructor( filename )

    !derived_type% finalize_sph_id_ptr => derived_type% finalize_sph_id

    CALL derived_type% sanity_check()


  END PROCEDURE initialize


  !MODULE PROCEDURE finalize_sph
  !
  !  !************************************************
  !  !
  !  !# This PROCEDURE calls the constructor of the
  !  !  [[idbase]]-extended type and the SUBROUTINE
  !  !  [[idbase:sanity_check]] afterwards. It is recommended
  !  !  to use this SUBROUTINE to construct objects of
  !  !  [[idbase]]-extended type since the sanity check is
  !  !  performed automatically.
  !  !
  !  !  FT 8.11.2021
  !  !
  !  !************************************************
  !
  !
  !  IMPLICIT NONE
  !
  !
  !  CALL derived_type% finalize_sph_id()
  !
  !
  !END PROCEDURE finalize_sph


END SUBMODULE initialization
