! File:         submodule_standard_tpo_formulation_io.f90
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

SUBMODULE (standard_tpo_formulation) io

  !***************************************************
  !
  !# This submodule contains the implementation of the
  !  methods of TYPE tpo_formulation that handle I/O
  !  (input/output)
  !
  !  FT 5.11.2021
  !
  !***************************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE print_summary

    !************************************************
    !
    !# Prints a summary of the properties of the refined mesh,
    !  and optionally, to a formatted file whose name
    !  is given as the optional argument `filename`
    !
    !  FT 5.11.2021
    !
    !************************************************

    USE tensor, ONLY: jx, jy, jz

    IMPLICIT NONE

    INTEGER:: l, last_level, i_matter

    last_level= THIS% get_nlevels()

    PRINT *, " * Spacetime:"
    PRINT *
    PRINT *, "   Number of refinement levels= ", last_level
    PRINT *
    PRINT *, "   Number of grid points on each level= ", &
             THIS% get_ngrid_x( 1 ), "**3"
    PRINT *
    DO l= 1, last_level, 1
      PRINT *, "   Resolution on level ", l, "= ", THIS% get_dx(l)
    ENDDO
    PRINT *
    DO l= 1, last_level, 1
      PRINT *, "   x boundary of level ", l, "= ", THIS% get_xR(l)
      PRINT *, "   y boundary of level ", l, "= ", THIS% get_yR(l)
      PRINT *, "   z boundary of level ", l, "= ", THIS% get_zR(l)
    ENDDO
    PRINT *
    DO i_matter= 1, THIS% n_matter, 1
      PRINT *, "   Number of grid points across the x-axis-diameter of ", &
               "matter object ", i_matter, "=", THIS% npoints_xaxis(i_matter)
    ENDDO
    PRINT *
    DO l= 1, last_level, 1
      PRINT *, "   x component of the ADM momentum ", &
               "on level    ", &
               l, "= ", THIS% MC_int(l,jx)
      PRINT *, "   y component of the ADM momentum ", &
               "on level    ", &
               l, "= ", THIS% MC_int(l,jy)
      PRINT *, "   z component of the ADM momentum ", &
               "on level    ", &
               l, "= ", THIS% MC_int(l,jz)
      PRINT *, "   Euclidean norm of the ADM momentum ", &
               "on level ", l, "= ", NORM2( THIS% MC_int(l,:), DIM=1 )
    ENDDO
    PRINT *

  END PROCEDURE print_summary


END SUBMODULE io