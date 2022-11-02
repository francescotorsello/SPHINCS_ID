! File:         module_cauchy_convergence_test.f90
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

MODULE cauchy_convergence_test

  !***********************************************************
  !
  !# This MODULE collects PROCEDURES used in the PROGRAM
  !  convergence_test
  !
  !  FT 22.09.2022
  !
  !  @note This MODULE was created on 22.09.2022 to collect
  !        already existing PROCEDURES
  !
  !***********************************************************


  USE standard_tpo_formulation, ONLY: tpo
  USE utility,                  ONLY: spacetime_path, ios, err_msg, &
                                      export_constraints_xy, &
                                      export_constraints_x, run_id


  IMPLICIT NONE


  INTEGER, PARAMETER:: use_constraints_on_mesh          = 1
  !# Parameter that identifies the case where the Cauchy convergence test
  !  is performed using the constraints computed using ID read on the
  !  refined mesh. Neither the particle ID nor the mapping are used.
  INTEGER, PARAMETER:: use_constraints_with_mapped_hydro= 2
  !# Parameter that identifies the case where the Cauchy convergence test
  !  is performed using the constraints computed using the hydro ID mapped from
  !  the particles to the refined mesh.

  DOUBLE PRECISION, PARAMETER:: tol= 1.D-10
  !# Tolerance used as an upper bound for the relative difference between
  !  the coordinates of 2 grid points, to consider them the same point


  INTERFACE find_shared_grid
  !# Generic PROCEDURE to compute the grid points shared by all the grids
  !  used to perform a Cauchy convergence test


    MODULE SUBROUTINE find_shared_grid_unknown_sol &
      ( tpo_coarse, tpo_medium, tpo_fine, num, den, ref_lev, shared_grid )
    !# Find the grid points shared by the 3 grids used in the
    !  Cauchy convergence test, when the exact solution is not
    !  known. The ratio between the grid spacings is `num/den`

      CLASS(tpo),       INTENT(INOUT):: tpo_coarse
      CLASS(tpo),       INTENT(INOUT):: tpo_medium
      CLASS(tpo),       INTENT(INOUT):: tpo_fine
      DOUBLE PRECISION, INTENT(IN)   :: num
      DOUBLE PRECISION, INTENT(IN)   :: den
      INTEGER,          INTENT(IN)   :: ref_lev

      DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(OUT):: &
      shared_grid

    END SUBROUTINE find_shared_grid_unknown_sol


    MODULE SUBROUTINE find_shared_grid_known_sol &
      ( tpo_coarse, tpo_fine, num, den, ref_lev, shared_grid )
    !# Find the grid points shared by the 2 grids used in the
    !  Cauchy convergence test, when the exact solution is
    !  known. The ratio between the grid spacings is `num/den`

      CLASS(tpo),       INTENT(INOUT):: tpo_coarse
      CLASS(tpo),       INTENT(INOUT):: tpo_fine
      DOUBLE PRECISION, INTENT(IN)   :: num
      DOUBLE PRECISION, INTENT(IN)   :: den
      INTEGER,          INTENT(IN)   :: ref_lev

      DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(OUT):: &
      shared_grid

    END SUBROUTINE find_shared_grid_known_sol


  END INTERFACE find_shared_grid


  INTERFACE perform_cauchy_convergence_test
  !# Generic PROCEDURE to perform the Cauchy convergence test


    MODULE SUBROUTINE perform_cauchy_convergence_test_unknown_sol &
      ( tpo_coarse, tpo_medium, tpo_fine, use_constraints, num, den, ref_lev )
    !# Perform the Cauchy convergence test when the exact solution is not
    !  known. The ratio between the grid spacings is `num/den`.

      CLASS(tpo),       INTENT(INOUT):: tpo_coarse
      CLASS(tpo),       INTENT(INOUT):: tpo_medium
      CLASS(tpo),       INTENT(INOUT):: tpo_fine
      INTEGER,          INTENT(IN)   :: use_constraints
      DOUBLE PRECISION, INTENT(IN)   :: num
      DOUBLE PRECISION, INTENT(IN)   :: den
      INTEGER,          INTENT(IN)   :: ref_lev

    END SUBROUTINE perform_cauchy_convergence_test_unknown_sol


    MODULE SUBROUTINE perform_cauchy_convergence_test_known_sol &
      ( tpo_coarse, tpo_fine, use_constraints, num, den, ref_lev )
    !# Perform the Cauchy convergence test when the exact solution is
    !  known. The ratio between the grid spacings is `num/den`.

      IMPLICIT NONE

      CLASS(tpo),       INTENT(INOUT):: tpo_coarse
      CLASS(tpo),       INTENT(INOUT):: tpo_fine
      INTEGER,          INTENT(IN)   :: use_constraints
      DOUBLE PRECISION, INTENT(IN)   :: num
      DOUBLE PRECISION, INTENT(IN)   :: den
      INTEGER,          INTENT(IN)   :: ref_lev

    END SUBROUTINE perform_cauchy_convergence_test_known_sol


  END INTERFACE perform_cauchy_convergence_test


END MODULE cauchy_convergence_test
