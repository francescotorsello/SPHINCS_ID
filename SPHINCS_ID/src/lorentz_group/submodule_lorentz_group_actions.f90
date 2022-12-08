! File:         submodule_lorentz_group_actions.f90
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

SUBMODULE(lorentz_group) actions


  !***********************************************************
  !
  !# This submodule contains the implementations of the
  !  actions of the TYPES representing the members of
  !  the Lorentz group defined in MODULE lorentz_group,
  !  presently boosts and spatial rotations, to geometrical
  !  objects
  !
  !  FT 08.12.2022
  !
  !***********************************************************


  IMPLICIT NONE



  CONTAINS



  MODULE PROCEDURE apply_boost_to_vector

    !***********************************************************
    !
    !# Implements the action of a boost on a \(4-\)vector
    !
    !  FT 08.12.2022
    !
    !  @todo implement matrix multiplications to make all this
    !        computations simpler and safer?
    !
    !***********************************************************


    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(3):: u_spatial

    u_spatial= [u(1), u(2), u(3)]

    !PRINT *, u
    !PRINT *, u_spatial

    boosted_u(0)= this% lambda*u(0) &
                + euclidean_inner_product( this% p, u_spatial )

    boosted_u(1)= this% p(1)*u(0) &
                + euclidean_inner_product( &
                    [this% lambda_s(1), this% lambda_s(2), this% lambda_s(3)], &
                    u_spatial )

    boosted_u(2)= this% p(2)*u(0) &
                + euclidean_inner_product( &
                    [this% lambda_s(2), this% lambda_s(4), this% lambda_s(5)], &
                    u_spatial )

    boosted_u(3)= this% p(3)*u(0) &
                + euclidean_inner_product( &
                    [this% lambda_s(3), this% lambda_s(5), this% lambda_s(6)], &
                    u_spatial )

  END PROCEDURE apply_boost_to_vector


  MODULE PROCEDURE apply_boost_similarity

    !***********************************************************
    !
    !# Implements the action of a boost as a similarity
    !  on a linear operator
    !
    !  FT 08.12.2022
    !
    !***********************************************************


    IMPLICIT NONE

    INTEGER:: i, j
    DOUBLE PRECISION, DIMENSION(4):: row
    DOUBLE PRECISION, DIMENSION(4):: column
    DOUBLE PRECISION, DIMENSION(4,4):: tmp(0:3,0:3)

    DO i= 0, 3, 1
      DO j= 0, 3, 1

        row   = this% inv_matrix(i,:)
        column= t(:,j)

        tmp(i,j)= row_by_column(row,column)

      ENDDO
    ENDDO

    DO i= 0, 3, 1
      DO j= 0, 3, 1

        row   = tmp(i,:)
        column= this% matrix(:,j)

        boosted_t(i,j)= row_by_column(row,column)

      ENDDO
    ENDDO

  END PROCEDURE apply_boost_similarity


  MODULE PROCEDURE apply_boost_congruence

    !***********************************************************
    !
    !# Implements the action of a boost as a congruence
    !  on a metric
    !
    !  FT 08.12.2022
    !
    !***********************************************************


    IMPLICIT NONE

    INTEGER:: i, j
    DOUBLE PRECISION, DIMENSION(4):: row
    DOUBLE PRECISION, DIMENSION(4):: column
    DOUBLE PRECISION, DIMENSION(4,4):: tmp(0:3,0:3)

    DO i= 0, 3, 1
      DO j= 0, 3, 1

        row   = this% matrix(:,i) ! Transpose
        column= t(:,j)

        tmp(i,j)= row_by_column(row,column)

      ENDDO
    ENDDO

    DO i= 0, 3, 1
      DO j= 0, 3, 1

        row   = tmp(i,:)
        column= this% matrix(:,j)

        boosted_t(i,j)= row_by_column(row,column)

      ENDDO
    ENDDO

  END PROCEDURE apply_boost_congruence


END SUBMODULE actions
