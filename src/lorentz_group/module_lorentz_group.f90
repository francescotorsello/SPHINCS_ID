! File:         module_lorentz_group.f90
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

MODULE lorentz_group


  !***********************************************************
  !
  !# This module contains the definitions of the TYPES
  !  representing objects in the Lorentz group, that is,
  !  the isometry group of the Lorentz metric
  !
  !  FT 23.02.2022
  !
  !***********************************************************


  USE utility,  ONLY: zero, one
  USE tensor,   ONLY: n_sym3x3, n_sym4x4


  IMPLICIT NONE


  DOUBLE PRECISION, DIMENSION(n_sym4x4), PARAMETER:: &
    eta= [-one, zero, zero, zero, &
                one,  zero, zero, &
                      one,  zero, &
                            one]
  !! Minkowski metric


  TYPE, ABSTRACT:: lorentz_transformation
  !! TYPE representing a 4D, proper, orthochronous Lorentz transformation

    PRIVATE

    DOUBLE PRECISION, DIMENSION(4,4):: matrix(0:3,0:3)
    !! \(4\times 4\) matrix representing the Lorentz transformation
    DOUBLE PRECISION, DIMENSION(4,4):: inv_matrix(0:3,0:3)
    !! \(4\times 4\) matrix representing the inverse Lorentz transformation
    DOUBLE PRECISION, DIMENSION(4,4):: tr_matrix(0:3,0:3)
    !# Transpose of the \(4\times 4\) matrix representing the Lorentz
    !  transformation


    CONTAINS


    PROCEDURE, PUBLIC, NON_OVERRIDABLE:: apply_to_vector
    !! Action of the [[lorentz_transformation]] on a \(4\)-vector

    GENERIC, PUBLIC:: apply_as_similarity => apply_as_similarity_to_tensor, &
                                          apply_as_similarity_to_symrank2_tensor
    !# Generic procedure to apply the [[lorentz_transformation]] as a similarity

    PROCEDURE, NON_OVERRIDABLE:: apply_as_similarity_to_tensor
    !# Action of the [[lorentz_transformation]] as a similarity on a generic
    !  tensor

    PROCEDURE, NON_OVERRIDABLE:: apply_as_similarity_to_symrank2_tensor
    !# Action of the [[lorentz_transformation]] as a similarity on a
    !  \(10\)-vector storing the components of a symmetric \(4\times 4\) tensor

    GENERIC, PUBLIC:: apply_as_congruence => apply_as_congruence_to_tensor, &
                                          apply_as_congruence_to_symrank2_tensor
    !# Generic procedure to apply the [[lorentz_transformation]] as a congruence

    PROCEDURE, NON_OVERRIDABLE:: apply_as_congruence_to_tensor
    !# Action of the [[lorentz_transformation]] as a congruence on a generic
    !  purely covariant tensor

    PROCEDURE, NON_OVERRIDABLE:: apply_as_congruence_to_symrank2_tensor
    !# Action of the [[lorentz_transformation]] as a congruence on a
    !  \(10\)-vector storing the components of a symmetric, purely covariant,
    !  \(4\times 4\) tensor

  END TYPE lorentz_transformation


  INTERFACE


    MODULE FUNCTION apply_to_vector(this, u) RESULT(transformed_u)
    !! Action of the [[lorentz_transformation]] on a \(4\)-vector

      CLASS(lorentz_transformation),           INTENT(IN):: this
      !! [[lorentz_transformation]] object to apply
      DOUBLE PRECISION, DIMENSION(4), INTENT(IN):: u(0:3)
      !! \(4\)-vector to be boosted
      DOUBLE PRECISION, DIMENSION(4):: transformed_u(0:3)
      !! Boosted \(4\)-vector

    END FUNCTION apply_to_vector


    MODULE FUNCTION apply_as_similarity_to_tensor(this, t) RESULT(transformed_t)
    !! Action of the [[lorentz_transformation]] on a \(4\)-vector

      CLASS(lorentz_transformation),             INTENT(IN):: this
      !! [[lorentz_transformation]] object to apply
      DOUBLE PRECISION, DIMENSION(4,4), INTENT(IN):: t(0:3,0:3)
      !! \(4\times 4\) tensor to be boosted
      DOUBLE PRECISION, DIMENSION(4,4):: transformed_t(0:3,0:3)
      !! Boosted \(4\times 4\) tensor

    END FUNCTION apply_as_similarity_to_tensor


    MODULE FUNCTION apply_as_congruence_to_tensor(this, t) RESULT(transformed_t)
    !! Action of the [[lorentz_transformation]] on a \(4\)-vector

      CLASS(lorentz_transformation),             INTENT(IN):: this
      !! [[lorentz_transformation]] object to apply
      DOUBLE PRECISION, DIMENSION(4,4), INTENT(IN):: t(0:3,0:3)
      !! \(4\times 4\) tensor to be boosted
      DOUBLE PRECISION, DIMENSION(4,4):: transformed_t(0:3,0:3)
      !! Boosted \(4\times 4\) tensor

    END FUNCTION apply_as_congruence_to_tensor


    MODULE FUNCTION apply_as_similarity_to_symrank2_tensor(this, t) &
      RESULT(transformed_t)
    !! Action of the [[lorentz_transformation]] on a \(4\)-vector

      CLASS(lorentz_transformation),                  INTENT(IN):: this
      !! [[lorentz_transformation]] object to apply
      DOUBLE PRECISION, DIMENSION(n_sym4x4), INTENT(IN):: t
      !# \(10\)-vector storing the components of the symmetric \(4\times 4\)
      !  tensor to be boosted
      DOUBLE PRECISION, DIMENSION(n_sym4x4):: transformed_t
      !# \(10\)-vector storing the components of the boosted symmetric
      !  \(4\times 4\) tensor

    END FUNCTION apply_as_similarity_to_symrank2_tensor


    MODULE FUNCTION apply_as_congruence_to_symrank2_tensor(this, t) &
      RESULT(transformed_t)
    !! Action of the [[lorentz_transformation]] on a \(4\)-vector

      CLASS(lorentz_transformation),                  INTENT(IN):: this
      !! [[lorentz_transformation]] object to apply
      DOUBLE PRECISION, DIMENSION(n_sym4x4), INTENT(IN):: t
      !# \(10\)-vector storing the components of the symmetric \(4\times 4\)
      !  tensor to be boosted
      DOUBLE PRECISION, DIMENSION(n_sym4x4):: transformed_t
      !# \(10\)-vector storing the components of the boosted symmetric
      !  \(4\times 4\) tensor

    END FUNCTION apply_as_congruence_to_symrank2_tensor


  END INTERFACE


  !--------------!
  !--  BOOSTS  --!
  !--------------!

  TYPE, EXTENDS(lorentz_transformation):: lorentz_boost
  !! TYPE representing a Lorentz boost

    PRIVATE

    DOUBLE PRECISION, DIMENSION(3):: v
    !! Spatial velocity that determines the boost

    DOUBLE PRECISION:: v_speed
    !! Euclidean norm of [[lorentz_boost:v]] (its speed)

    DOUBLE PRECISION:: lambda
    !! Lorentz factor

    DOUBLE PRECISION, DIMENSION(3):: p
    !! Spatial vector equal to \(\lambda \,v\)

    DOUBLE PRECISION, DIMENSION(n_sym3x3):: lambda_s
    !! Spatial part of the Lorentz boost
    DOUBLE PRECISION, DIMENSION(n_sym3x3):: inv_lambda_s
    !! Spatial part of the inverse Lorentz boost


    CONTAINS


    PROCEDURE:: compute_boost_matrices
    !# Computes the spatial part of the matrix of the Lorentz
    !  boost, and its whole matrix, starting from the vector
    !  \(p\)

    PROCEDURE, PUBLIC:: get_lambda
    !! Returns the Lorentz factor [[lorentz_boost:lambda]]


  END TYPE lorentz_boost


  INTERFACE lorentz_boost

    MODULE FUNCTION construct_boost(v) RESULT(boost)

      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: v
      !! Spatial velocity that determines the boost
      TYPE(lorentz_boost):: boost
      !! [[lorentz_boost]] object to be constructed

    END FUNCTION construct_boost

    MODULE FUNCTION construct_boost_components(vx, vy, vz) RESULT(boost)

      DOUBLE PRECISION, INTENT(IN):: vx
      !! \(x\) component of the spatial velocity that determines the boost
      DOUBLE PRECISION, INTENT(IN):: vy
      !! \(y\) component of the spatial velocity that determines the boost
      DOUBLE PRECISION, INTENT(IN):: vz
      !! \(z\) component of the spatial velocity that determines the boost
      TYPE(lorentz_boost):: boost
      !! [[lorentz_boost]] object to be constructed

    END FUNCTION construct_boost_components

  END INTERFACE lorentz_boost


  INTERFACE

    MODULE PURE SUBROUTINE compute_boost_matrices(this, p, lambda_s, matrix)
    !! Compute the matrices for the [[lorentz_boost]]

      CLASS(lorentz_boost),                  INTENT(INOUT):: this
      DOUBLE PRECISION, DIMENSION(3),        INTENT(IN)   :: p
      !! [[lorentz_boost:lambda]]\(\times\)[[lorentz_boost:v]]
      DOUBLE PRECISION, DIMENSION(n_sym3x3), INTENT(OUT)  :: lambda_s
      !! Spatial part of the Lorentz boost
      DOUBLE PRECISION, DIMENSION(4,4),      INTENT(OUT)  :: matrix(0:3,0:3)
      !! \(4\times 4\) matrix representing the Lorentz boost

    END SUBROUTINE compute_boost_matrices

    PURE MODULE FUNCTION get_lambda(this) RESULT(lambda)
    !! Returns the Lorentz factor [[lorentz_boost:lambda]]

      CLASS(lorentz_boost), INTENT(IN):: this
      !! [[lorentz_boost]] object owning this FUNCTION
      DOUBLE PRECISION:: lambda
      !! Lorentz factor [[lorentz_boost:lambda]]

    END FUNCTION get_lambda

  END INTERFACE


  !------------------------!
  !--  SPATIAL ROTATIONS --!
  !------------------------!

  TYPE, EXTENDS(lorentz_transformation):: spatial_rotation
  !! TYPE representing a spatial rotation

    PRIVATE

    DOUBLE PRECISION, DIMENSION(3):: euler_angles
    !# Euler angles that define the rotation around the \(x,y,z\) axes,
    !  in this order

    DOUBLE PRECISION, DIMENSION(3,3):: r_x
    !! Rotation operator around the \(x\) axis
    DOUBLE PRECISION, DIMENSION(3,3):: r_y
    !! Rotation operator around the \(y\) axis
    DOUBLE PRECISION, DIMENSION(3,3):: r_z
    !! Rotation operator around the \(z\) axis
    DOUBLE PRECISION, DIMENSION(3,3):: r
    !! Full rotation operator

    DOUBLE PRECISION, DIMENSION(3,3):: tr_r
    !! Transpose of the full inverse rotation operator

    DOUBLE PRECISION, DIMENSION(3,3):: inv_r_x
    !! Inverse rotation operator around the \(x\) axis
    DOUBLE PRECISION, DIMENSION(3,3):: inv_r_y
    !! Inverse rotation operator around the \(y\) axis
    DOUBLE PRECISION, DIMENSION(3,3):: inv_r_z
    !! Inverse rotation operator around the \(z\) axis
    DOUBLE PRECISION, DIMENSION(3,3):: inv_r
    !! Full inverse rotation operator


    CONTAINS


    PROCEDURE:: compute_rotation_matrices
    !! Computes the spatial part of the matrix of the Lorentz
    !  boost, and its whole matrix, starting from the vector
    !  \(p\)


  END TYPE spatial_rotation


  INTERFACE spatial_rotation

    MODULE FUNCTION construct_rotation(euler_angles) RESULT(rotation)

      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: euler_angles
      !# Euler angles that define the rotation around the \(x,y,z\) axes,
      !  in this order
      TYPE(spatial_rotation):: rotation
      !! [[spatial_rotation]] object to be constructed

    END FUNCTION construct_rotation

    MODULE FUNCTION construct_rotation_angles(alpha, beta, gamma) &
      RESULT(rotation)

      DOUBLE PRECISION, INTENT(IN):: alpha
      !# Euler angle that defines the rotation around the \(x\) axis
      DOUBLE PRECISION, INTENT(IN):: beta
      !# Euler angle that defines the rotation around the \(y\) axis
      DOUBLE PRECISION, INTENT(IN):: gamma
      !# Euler angle that defines the rotation around the \(z\) axis
      TYPE(spatial_rotation):: rotation
      !! [[spatial_rotation]] object to be constructed

    END FUNCTION construct_rotation_angles

  END INTERFACE spatial_rotation


  INTERFACE

    MODULE SUBROUTINE compute_rotation_matrices &
      (this, euler_angles, r_x, r_y, r_z, r, matrix, &
       inv_r_x, inv_r_y, inv_r_z, inv_r, inv_matrix)
    !! Compute the matrices for the [[spatial_rotation]]

      CLASS(spatial_rotation),               INTENT(INOUT):: this
      !! [[spatial_rotation]] object to compuet the matrices for
      DOUBLE PRECISION, DIMENSION(3),        INTENT(IN)   :: euler_angles
      !! [[lorentz_boost:lambda]]\(\times\)[[lorentz_boost:v]]
      DOUBLE PRECISION, DIMENSION(3,3),      INTENT(OUT)  :: r_x
      !! Rotation operator around the \(x\) axis
      DOUBLE PRECISION, DIMENSION(3,3),      INTENT(OUT)  :: r_y
      !! Rotation operator around the \(y\) axis
      DOUBLE PRECISION, DIMENSION(3,3),      INTENT(OUT)  :: r_z
      !! Rotation operator around the \(z\) axis
      DOUBLE PRECISION, DIMENSION(3,3),      INTENT(OUT)  :: r
      !! Full \(3\times 3\) rotation operator
      DOUBLE PRECISION, DIMENSION(4,4),      INTENT(OUT)  :: matrix(0:3,0:3)
      !! \(4\times 4\) matrix representing the Lorentz boost
      DOUBLE PRECISION, DIMENSION(3,3),      INTENT(OUT)  :: inv_r_x
      !! Inverse rotation operator around the \(x\) axis
      DOUBLE PRECISION, DIMENSION(3,3),      INTENT(OUT)  :: inv_r_y
      !! Inverse rotation operator around the \(y\) axis
      DOUBLE PRECISION, DIMENSION(3,3),      INTENT(OUT)  :: inv_r_z
      !! Inverse rotation operator around the \(z\) axis
      DOUBLE PRECISION, DIMENSION(3,3),      INTENT(OUT)  :: inv_r
      !! Inverse of the full \(3\times 3\) rotation operator
      DOUBLE PRECISION, DIMENSION(4,4),      INTENT(OUT)  :: inv_matrix(0:3,0:3)
      !! Inverse of the \(4\times 4\) matrix representing the Lorentz boost

    END SUBROUTINE compute_rotation_matrices

  END INTERFACE



  CONTAINS



  PURE FUNCTION euclidean_inner_product( u, v ) RESULT( inner_product )

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: u
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: v
    DOUBLE PRECISION:: inner_product

    inner_product= u(1)*v(1) + u(2)*v(2) + u(3)*v(3)

  END FUNCTION euclidean_inner_product


  PURE FUNCTION minkowski_inner_product( u, v ) RESULT( inner_product )

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(4), INTENT(IN):: u(0:3)
    DOUBLE PRECISION, DIMENSION(4), INTENT(IN):: v(0:3)
    DOUBLE PRECISION:: inner_product

    inner_product= - u(0)*v(0) + u(1)*v(1) + u(2)*v(2) + u(3)*v(3)

  END FUNCTION minkowski_inner_product


  PURE FUNCTION minkowski_sqnorm( u ) RESULT( sqnorm )

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(4), INTENT(IN):: u(0:3)
    DOUBLE PRECISION:: sqnorm

    sqnorm= minkowski_inner_product(u,u)

  END FUNCTION minkowski_sqnorm


  PURE FUNCTION row_by_column( u, v ) RESULT( res )

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(:), INTENT(IN):: u
    !! Row
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN):: v
    !! Column
    DOUBLE PRECISION:: res

    INTEGER:: i, n

    n= SIZE(u)

    res= 0.D0
    DO i= 1, n, 1
      res= res + u(i)*v(i)
    ENDDO

  END FUNCTION row_by_column


  FUNCTION square_matrix_multiplication( a, b ) RESULT( c )

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN):: a
    !! First matrix
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN):: b
    !! Second matrix

    DOUBLE PRECISION, DIMENSION(SIZE(a(1,:))):: row
    DOUBLE PRECISION, DIMENSION(SIZE(a(1,:))):: column
    DOUBLE PRECISION, DIMENSION(SIZE(a(1,:)),SIZE(a(1,:))):: c

    INTEGER:: i, j, n

    n= SIZE(a(1,:))

    IF(n /= SIZE(a(:,1)))THEN
      PRINT *, "** ERROR! Matrix a is not square!"
      PRINT *, " * Number of rows=", SIZE(a(:,1))
      PRINT *, "   Number of columns=", n
      PRINT *, " * Stopping..."
      PRINT *
      STOP
    ENDIF
    IF(SIZE(b(1,:)) /= SIZE(b(:,1)))THEN
      PRINT *, "** ERROR! Matrix b is not square!"
      PRINT *, " * Number of rows=", SIZE(b(:,1))
      PRINT *, "   Number of columns=", SIZE(b(1,:))
      PRINT *, " * Stopping..."
      PRINT *
      STOP
    ENDIF
    IF(n /= SIZE(b(:,1)))THEN
      PRINT *, "** ERROR! Matrix a and b do not have compatible dimensions!"
      PRINT *, "   Number of columns of a=", n
      PRINT *, " * Number of rows of b=", SIZE(b(:,1))
      PRINT *, " * Stopping..."
      PRINT *
      STOP
    ENDIF

    DO i= 1, n, 1
      DO j= 1, n, 1

        row   = a(i,:)
        column= b(:,j)

        c(i,j)= row_by_column(row,column)

      ENDDO
    ENDDO


  END FUNCTION square_matrix_multiplication


END MODULE lorentz_group
