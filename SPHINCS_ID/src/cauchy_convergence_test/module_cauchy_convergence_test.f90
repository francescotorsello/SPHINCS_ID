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


  IMPLICIT NONE


  DOUBLE PRECISION, PARAMETER:: tol= 1.D-10


  INTERFACE find_shared_grid


    MODULE SUBROUTINE find_shared_grid_unknown_sol &
    ( tpo_coarse, tpo_medium, tpo_fine, num, den, ref_lev, shared_grid )
    !#

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
    !#

      CLASS(tpo),       INTENT(INOUT):: tpo_coarse
      CLASS(tpo),       INTENT(INOUT):: tpo_fine
      DOUBLE PRECISION, INTENT(INOUT):: num
      DOUBLE PRECISION, INTENT(IN)   :: den
      INTEGER,          INTENT(IN)   :: ref_lev

      DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(OUT):: &
      shared_grid

    END SUBROUTINE find_shared_grid_known_sol


  END INTERFACE find_shared_grid


  CONTAINS


  MODULE PROCEDURE find_shared_grid_unknown_sol

    !***********************************************************
    !
    !# Find the grid points shared by the 3 grids used in
    !  Cauchy convergence test, when the exact solution is not
    !  known.
    !  The ratio between the grid spacings is
    !  `num/den`
    !
    !  @todo add computation that justifies the algorithm
    !
    !***********************************************************

    IMPLICIT NONE

    INTEGER:: nx
    INTEGER:: ny
    INTEGER:: nz

    INTEGER:: i, j, k

    DOUBLE PRECISION, DIMENSION(3):: point_medium
    DOUBLE PRECISION, DIMENSION(3):: point_fine

    nx= tpo_coarse% get_ngrid_x(ref_lev)
    ny= tpo_coarse% get_ngrid_y(ref_lev)
    nz= tpo_coarse% get_ngrid_z(ref_lev)

    nx= FLOOR( DBLE( nx - 1 )/den**2 ) + 1
    ny= FLOOR( DBLE( ny - 1 )/den**2 ) + 1
    nz= FLOOR( DBLE( nz - 1 )/den**2 ) + 1

    ALLOCATE( shared_grid( nx, ny, nz, 3 ) )

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( tpo_coarse, tpo_medium, tpo_fine, ref_lev, &
    !$OMP                     shared_grid, den, num, nx, ny, nz ) &
    !$OMP             PRIVATE( i, j, k, point_medium, point_fine )
    DO k= 0, nz - 1, 1
      DO j= 0, ny - 1, 1
        DO i= 0, nx - 1, 1

          shared_grid( 1 + i, 1 + j, 1 + k, : ) = &
                     tpo_coarse% get_grid_point( 1 + INT(den**2)*i, &
                                                 1 + INT(den**2)*j, &
                                                 1 + INT(den**2)*k, ref_lev )

          point_medium= tpo_medium% get_grid_point( 1 + INT(num*den)*i, &
                                                    1 + INT(num*den)*j, &
                                                    1 + INT(num*den)*k, ref_lev)

          point_fine  = tpo_fine% get_grid_point( 1 + INT(num**2)*i, &
                                                  1 + INT(num**2)*j, &
                                                  1 + INT(num**2)*k, ref_lev )

          IF(  ABS(shared_grid(1+i,1+j,1+k, 1)-point_medium(1)) &
              /ABS(point_medium(1)) > tol &
          .OR. ABS(shared_grid(1+i,1+j,1+k, 1)-point_fine(1)) &
              /ABS(point_fine(1)) > tol &
          .OR. ABS(shared_grid(1+i,1+j,1+k, 2)-point_medium(2)) &
              /ABS(point_medium(2)) > tol &
          .OR. ABS(shared_grid(1+i,1+j,1+k, 2)-point_fine(2)) &
              /ABS(point_fine(2)) > tol &
          .OR. ABS(shared_grid(1+i,1+j,1+k, 3)-point_medium(3)) &
              /ABS(point_medium(3)) > tol &
          .OR. ABS(shared_grid(1+i,1+j,1+k, 3)-point_fine(3)) &
              /ABS(point_fine(3)) > tol &

          )THEN

            PRINT *
            PRINT *, "**ERROR in SUBROUTINE find_shared_grid_unknown_sol! ", &
                     "The grid functions in the Cauchy ", &
                     "convergence test are not evaluated at the ", &
                     "same grid point at (i,j,k)=(", i, j, k, ")."
            PRINT *, shared_grid(1+i,1+j,1+k, 1), point_medium(1), point_fine(1)
            PRINT *, shared_grid(1+i,1+j,1+k, 2), point_medium(2), point_fine(2)
            PRINT *, shared_grid(1+i,1+j,1+k, 3), point_medium(3), point_fine(3)
            PRINT *
            STOP

          ENDIF

        ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END PROCEDURE find_shared_grid_unknown_sol


  MODULE PROCEDURE find_shared_grid_known_sol

    !***********************************************************
    !
    !# This MODULE collects PROCEDURES used in the PROGRAM
    !  convergence_test
    !
    !***********************************************************

  END PROCEDURE find_shared_grid_known_sol


END MODULE cauchy_convergence_test
