! File:         submodule_cauchy_convergence_test_shared_grid.f90
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

SUBMODULE (cauchy_convergence_test) shared_grid

  !********************************************
  !
  !# This submodule contains the implementation
  !  of the PROCEDURES in MODULE
  !  cauchy_convergence_test that find the
  !  shared grids
  !
  !  FT 22.09.2022
  !
  !********************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE find_shared_grid_unknown_sol

    !***********************************************************
    !
    !# Find the grid points shared by the 3 grids used in the
    !  Cauchy convergence test, when the exact solution is not
    !  known. The ratio between the grid spacings is `num/den`
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

    IF( ALLOCATED(shared_grid) ) DEALLOCATE(shared_grid)
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

          IF(  ABS(shared_grid(1+i,1+j,1+k, 1) - point_medium(1)) &
              /ABS(point_medium(1)) > tol &
          .OR. ABS(shared_grid(1+i,1+j,1+k, 1) - point_fine(1)) &
              /ABS(point_fine(1)) > tol &
          .OR. ABS(shared_grid(1+i,1+j,1+k, 2) - point_medium(2)) &
              /ABS(point_medium(2)) > tol &
          .OR. ABS(shared_grid(1+i,1+j,1+k, 2) - point_fine(2)) &
              /ABS(point_fine(2)) > tol &
          .OR. ABS(shared_grid(1+i,1+j,1+k, 3) - point_medium(3)) &
              /ABS(point_medium(3)) > tol &
          .OR. ABS(shared_grid(1+i,1+j,1+k, 3) - point_fine(3)) &
              /ABS(point_fine(3)) > tol &

          )THEN

            PRINT *
            PRINT *, "** ERROR in SUBROUTINE find_shared_grid_unknown_sol! ", &
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
    !# Find the grid points shared by the 2 grids used in the
    !  Cauchy convergence test, when the exact solution is
    !  known. The ratio between the grid spacings is `num/den`
    !
    !  @todo add computation that justifies the algorithm
    !
    !***********************************************************

    IMPLICIT NONE

    INTEGER:: nx
    INTEGER:: ny
    INTEGER:: nz

    INTEGER:: i, j, k

    DOUBLE PRECISION, DIMENSION(3):: point_fine

    nx= tpo_coarse% get_ngrid_x(ref_lev)
    ny= tpo_coarse% get_ngrid_y(ref_lev)
    nz= tpo_coarse% get_ngrid_z(ref_lev)

    nx= FLOOR( DBLE( nx - 1 )/den ) + 1
    ny= FLOOR( DBLE( ny - 1 )/den ) + 1
    nz= FLOOR( DBLE( nz - 1 )/den ) + 1

    IF( ALLOCATED(shared_grid) ) DEALLOCATE(shared_grid)
    ALLOCATE( shared_grid( nx, ny, nz, 3 ) )

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( tpo_coarse, tpo_fine, ref_lev, &
    !$OMP                     shared_grid, den, num, nx, ny, nz ) &
    !$OMP             PRIVATE( i, j, k, point_fine )
    DO k= 0, nz - 1, 1
      DO j= 0, ny - 1, 1
        DO i= 0, nx - 1, 1

          shared_grid( 1 + i, 1 + j, 1 + k, : ) = &
                      tpo_coarse% get_grid_point( 1 + INT(den)*i, &
                                                  1 + INT(den)*j, &
                                                  1 + INT(den)*k, ref_lev )
          point_fine= tpo_fine% get_grid_point( 1 + INT(num)*i, &
                                                1 + INT(num)*j, &
                                                1 + INT(num)*k, ref_lev )

          IF(  ABS(shared_grid(1+i,1+j,1+k, 1) - point_fine(1)) &
              /ABS(point_fine(1)) > tol &
          .OR. ABS(shared_grid(1+i,1+j,1+k, 2) - point_fine(2)) &
              /ABS(point_fine(2)) > tol &
          .OR. ABS(shared_grid(1+i,1+j,1+k, 3) - point_fine(3)) &
              /ABS(point_fine(3)) > tol &

          )THEN

            PRINT *
            PRINT *, "** ERROR in SUBROUTINE find_shared_grid_known_sol! ", &
                     "The grid functions in the Cauchy ", &
                     "convergence test are not evaluated at the ", &
                     "same grid point at (i,j,k)=(", i, j, k, ")."
            PRINT *, shared_grid(1+i,1+j,1+k, 1), point_fine(1)
            PRINT *, shared_grid(1+i,1+j,1+k, 2), point_fine(2)
            PRINT *, shared_grid(1+i,1+j,1+k, 3), point_fine(3)
            PRINT *
            STOP

          ENDIF

        ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END PROCEDURE find_shared_grid_known_sol


END SUBMODULE shared_grid
