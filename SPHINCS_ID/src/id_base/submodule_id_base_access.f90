! File:         submodule_id_base_access.f90
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

SUBMODULE (id_base) access

  !********************************************
  !
  !# Implementation of the methods of TYPE
  !  idbase to access PRIVATE data
  !
  !  FT 28.10.2021
  !
  !********************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE set_n_matter

    !************************************************
    !
    !# Sets [[idbase:n_matter]] to the given value
    !
    !  FT 28.10.2021
    !
    !************************************************

    IMPLICIT NONE

    this% n_matter= value

  END PROCEDURE set_n_matter


  MODULE PROCEDURE get_n_matter

    !************************************************
    !
    !# Returns [[idbase:n_matter]], the number of
    !  matter objects in the physical system
    !
    !  FT 28.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_n_matter= this% n_matter

  END PROCEDURE get_n_matter


  MODULE PROCEDURE set_one_lapse

    !************************************************
    !
    !# Sets [[idbase:one_lapse]] to the given value
    !
    !  FT 3.11.2021
    !
    !************************************************

    IMPLICIT NONE

    this% one_lapse= logic

  END PROCEDURE set_one_lapse


  MODULE PROCEDURE get_one_lapse

    !************************************************
    !
    !# Returns [[idbase:one_lapse]]
    !
    !  FT 3.11.2021
    !
    !************************************************

    IMPLICIT NONE

    get_one_lapse= this% one_lapse

  END PROCEDURE get_one_lapse


  MODULE PROCEDURE set_zero_shift

    !************************************************
    !
    !# Sets [[idbase:zero_shift]] to the given value
    !
    !  FT 3.11.2021
    !
    !************************************************

    IMPLICIT NONE

    this% zero_shift= logic

  END PROCEDURE set_zero_shift


  MODULE PROCEDURE get_zero_shift

    !************************************************
    !
    !# Returns [[idbase:zero_shift]]
    !
    !  FT 3.11.2021
    !
    !************************************************

    IMPLICIT NONE

    get_zero_shift= this% zero_shift

  END PROCEDURE get_zero_shift


  MODULE PROCEDURE set_cold_system

    !************************************************
    !
    !# Sets [[idbase:cold_system]] to the given value
    !
    !  FT 19.01.2021
    !
    !************************************************

    IMPLICIT NONE

    IF( .NOT.(value .EQV. .TRUE.) .AND. .NOT.(value .EQV. .FALSE.) )THEN

      PRINT *, "** ERROR in SUBROUTINE set_cold_system! The LOGICAL ", &
               "variable cold_system should be .TRUE. or .FALSE. ."
      PRINT *, " * Stopping..."
      PRINT *
      STOP

    ELSE

      this% cold_system= value

    ENDIF

  END PROCEDURE set_cold_system


  MODULE PROCEDURE get_cold_system

    !************************************************
    !
    !# Returns [[idbase:cold_system]], the `LOGICAL`
    !  variable at specifies if the system is cold
    !  (no thermal component)
    !
    !  FT 19.01.2021
    !
    !************************************************

    IMPLICIT NONE

    get_cold_system= this% cold_system

  END PROCEDURE get_cold_system


  MODULE PROCEDURE set_estimate_length_scale

    !************************************************
    !
    !# Sets [[idbase:estimate_length_scale]] to the
    !  given value
    !
    !  FT 18.02.2021
    !
    !************************************************

    IMPLICIT NONE

    IF( .NOT.(value .EQV. .TRUE.) .AND. .NOT.(value .EQV. .FALSE.) )THEN

      PRINT *, "** ERROR in SUBROUTINE set_estimate_length_scale! The ", &
               "LOGICAL variable cold_system should be .TRUE. or .FALSE. ."
      PRINT *, " * Stopping..."
      PRINT *
      STOP

    ELSE

      this% estimate_length_scale= value

    ENDIF

  END PROCEDURE set_estimate_length_scale


  MODULE PROCEDURE get_estimate_length_scale

    !************************************************
    !
    !# Returns [[idbase:estimate_length_scale]], the
    !  `LOGICAL` variable that specifies if a typical
    !  length scale, equal to the ratio of a field over
    !  its gradient, should be computed
    !
    !  FT 18.02.2021
    !
    !************************************************

    IMPLICIT NONE

    get_estimate_length_scale= this% estimate_length_scale

  END PROCEDURE get_estimate_length_scale


  MODULE PROCEDURE check_i_matter

    !************************************************
    !
    !# Checks that the given index `i_matter` is
    !  between 1 and [[idbase:n_matter]], included.
    !  If not, it stops the execution of the program.
    !
    !  FT 2.11.2021
    !
    !************************************************

    IMPLICIT NONE

    IF( i_matter < 1 )THEN

      PRINT *, "** ERROR! The index of the matter object is lower than 1. "
      PRINT *, "   It should be between 1 and n_matter= ", this% n_matter
      PRINT *, "   The index is ", i_matter
      PRINT *, "   Stopping..."
      PRINT *
      STOP

    ELSEIF( i_matter > this% n_matter )THEN

      PRINT *, "** ERROR! The index of the matter object is larger than ", &
               this% n_matter
      PRINT *, "   It should be between 1 and n_matter= ", this% n_matter
      PRINT *, "   The index is ", i_matter
      PRINT *, "   Stopping..."
      PRINT *
      STOP

    ENDIF

  END PROCEDURE check_i_matter


  MODULE PROCEDURE get_total_spatial_extent

    !************************************************
    !
    !# Return the coordinates
    !\(x_{\rm min},x_{\rm max},y_{\rm min},y_{\rm max},z_{\rm min},z_{\rm max}\)
    !  of a box containing the entire physical system
    !
    !  FT 28.10.2021
    !
    !************************************************

    USE utility, ONLY: one, five, ten

    IMPLICIT NONE

    DOUBLE PRECISION, PARAMETER:: stretch= five/(ten*ten*ten)

    INTEGER:: i_matter
    DOUBLE PRECISION, DIMENSION(3):: center_matter
    DOUBLE PRECISION, DIMENSION(6):: size_matter

    box(1)=   HUGE(one)
    box(3)=   HUGE(one)
    box(5)=   HUGE(one)
    box(2)= - HUGE(one)
    box(4)= - HUGE(one)
    box(6)= - HUGE(one)

    DO i_matter= 1, this% get_n_matter(), 1

      size_matter  = this% return_spatial_extent( i_matter )*(one + stretch)
      center_matter= this% return_center( i_matter )

      IF( center_matter(1) - size_matter(1) < box(1) ) &
                          box(1) = center_matter(1) - size_matter(1)
      IF( center_matter(1) + size_matter(2) > box(2) ) &
                          box(2) = center_matter(1) + size_matter(2)
      IF( center_matter(2) - size_matter(3) < box(3) ) &
                          box(3) = center_matter(2) - size_matter(3)
      IF( center_matter(2) + size_matter(4) > box(4) ) &
                          box(4) = center_matter(2) + size_matter(4)
      IF( center_matter(3) - size_matter(5) < box(5) ) &
                          box(5) = center_matter(3) - size_matter(5)
      IF( center_matter(3) + size_matter(6) > box(6) ) &
                          box(6) = center_matter(3) + size_matter(6)

    ENDDO

  END PROCEDURE get_total_spatial_extent


END SUBMODULE access
