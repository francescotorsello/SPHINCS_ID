! File:         submodule_lorentz_group_boost_actions.f90
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



  MODULE PROCEDURE apply_to_vector

    !***********************************************************
    !
    !# Implements the action of a boost on a \(4-\)vector
    !
    !  FT 08.12.2022
    !
    !***********************************************************


    IMPLICIT NONE

    INTEGER:: i
    DOUBLE PRECISION, DIMENSION(4):: row(0:3)
    DOUBLE PRECISION, DIMENSION(4):: column(0:3)

    !u_spatial= [u(1), u(2), u(3)]
    !
    !transformed_u(0)= this% lambda*u(0) &
    !            + euclidean_inner_product( this% p, u_spatial )
    !
    !transformed_u(1)= this% p(1)*u(0) &
    !            + euclidean_inner_product( &
    !                [this% lambda_s(1), this% lambda_s(2), this% lambda_s(3)], &
    !                u_spatial )
    !
    !transformed_u(2)= this% p(2)*u(0) &
    !            + euclidean_inner_product( &
    !                [this% lambda_s(2), this% lambda_s(4), this% lambda_s(5)], &
    !                u_spatial )
    !
    !transformed_u(3)= this% p(3)*u(0) &
    !            + euclidean_inner_product( &
    !                [this% lambda_s(3), this% lambda_s(5), this% lambda_s(6)], &
    !                u_spatial )

    column= u

    DO i= 0, 3, 1
      row         = this% matrix(i,:)
      transformed_u(i)= row_by_column(row,column)
    ENDDO

  END PROCEDURE apply_to_vector


  MODULE PROCEDURE apply_as_similarity_to_tensor

    !***********************************************************
    !
    !# Implements the action of a boost as a similarity
    !  on a linear operator
    !
    !  FT 08.12.2022
    !
    !***********************************************************

    IMPLICIT NONE

    transformed_t= square_matrix_multiplication( &
                square_matrix_multiplication(this% inv_matrix,t), this% matrix)

  END PROCEDURE apply_as_similarity_to_tensor


  MODULE PROCEDURE apply_as_similarity_to_symrank2_tensor

    !***********************************************************
    !
    !# Implements the action of a boost as a similarity
    !  on a \(10\)-vector storing the components of a symmetric
    !  \(4\times 4\) tensor
    !
    !  FT 09.12.2022
    !
    !***********************************************************

    USE metric_on_particles,  ONLY: gvec2mat, mat2gvec

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(4,4):: t_mat(0:3,0:3)

    CALL gvec2mat(t,t_mat)

    CALL mat2gvec(transformed_t, this% apply_as_similarity_to_tensor(t_mat))

  END PROCEDURE apply_as_similarity_to_symrank2_tensor


  MODULE PROCEDURE apply_as_congruence_to_tensor

    !***********************************************************
    !
    !# Implements the action of a boost as a congruence
    !  on a metric
    !
    !  FT 08.12.2022
    !
    !***********************************************************

    IMPLICIT NONE

    transformed_t= square_matrix_multiplication( &
                square_matrix_multiplication(this% tr_matrix,t), this% matrix)

  END PROCEDURE apply_as_congruence_to_tensor


  MODULE PROCEDURE apply_as_congruence_to_symrank2_tensor

    !***********************************************************
    !
    !# Implements the action of a boost as a similarity
    !  on a \(10\)-vector storing the components of a symmetric
    !  \(4\times 4\) tensor
    !
    !  FT 09.12.2022
    !
    !***********************************************************

    USE metric_on_particles,  ONLY: gvec2mat, mat2gvec

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(4,4):: t_mat(0:3,0:3)

    CALL gvec2mat(t, t_mat)

    CALL mat2gvec(transformed_t, this% apply_as_congruence_to_tensor(t_mat))

  END PROCEDURE apply_as_congruence_to_symrank2_tensor


END SUBMODULE actions
