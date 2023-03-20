! File:         submodule_lorentz_group_constructors.f90
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

SUBMODULE(lorentz_group) constructors


  !***********************************************************
  !
  !# This submodule contains the implementations of the
  !  constructors of the TYPES representing the members of
  !  the Lorentz group defined in MODULE lorentz_group,
  !  presently boosts and spatial rotations
  !
  !  FT 08.12.2022
  !
  !***********************************************************


  USE utility,  ONLY: zero, one


  IMPLICIT NONE



  CONTAINS



  !------------!
  !-- BOOSTS --!
  !------------!


  MODULE PROCEDURE construct_boost

    !***********************************************************
    !
    !# Construct a boost transformation from the spatial velocity
    !
    !  FT 08.12.2022
    !
    !***********************************************************

    IMPLICIT NONE

    INTEGER:: i, j

    DOUBLE PRECISION:: v_sqnorm, delta
    DOUBLE PRECISION, DIMENSION(4,4):: identity(0:3,0:3)
    LOGICAL:: is_identity

    v_sqnorm= euclidean_inner_product(v,v)
    boost% v_speed= SQRT(v_sqnorm)

    IF( boost% v_speed >= one )THEN
      PRINT *, "** ERROR! The speed given to construct_boost is larger", &
               " than, or equal to, 1!"
      PRINT *, "|v|^2=", boost% v_speed
      PRINT *, "* Stopping..."
      PRINT *
      STOP
    ENDIF

    boost% v= v

    boost% lambda= one/SQRT( one - v_sqnorm )

    boost% p= boost% lambda*boost% v

    CALL boost% compute_boost_matrices(boost% p, boost% lambda_s, boost% matrix)
    CALL boost% compute_boost_matrices(- boost% p, boost% inv_lambda_s, &
                                         boost% inv_matrix)

    boost% tr_matrix= TRANSPOSE(boost% matrix)

    identity= square_matrix_multiplication(boost% inv_matrix, boost% matrix)
    is_identity= .TRUE.
    DO i= 0, 3, 1
      DO j= 0, 3, 1
        IF( i == j )THEN
          delta= one
        ELSE
          delta= zero
        ENDIF
        is_identity= is_identity .AND. (identity(i,j) - delta) < 1.D-12
      ENDDO
    ENDDO
    IF( .NOT.is_identity )THEN

      PRINT *, "** ERROR! boost% inv_matrix is not the inverse of ", &
               "boost% matrix! Their product is:"
      PRINT *, "   id(0,:)=", identity(0,:)
      PRINT *, "   id(1,:)=", identity(1,:)
      PRINT *, "   id(2,:)=", identity(2,:)
      PRINT *, "   id(3,:)=", identity(3,:)
      PRINT *
      STOP

    ENDIF

  END PROCEDURE construct_boost


  MODULE PROCEDURE construct_boost_components

    !***********************************************************
    !
    !# Construct a boost transformation from the components
    !  of the spatial velocity
    !
    !  FT 10.12.2022
    !
    !***********************************************************

    IMPLICIT NONE

    boost= lorentz_boost([vx,vy,vz])

  END PROCEDURE construct_boost_components


  MODULE PROCEDURE compute_boost_matrices

    !***********************************************************
    !
    !# Computes the spatial part of the matrix of the Lorentz
    !  boost, and its whole matrix, starting from the vector
    !  \(p\)
    !
    !  FT 09.12.2022
    !
    !***********************************************************

    USE utility,  ONLY: one

    IMPLICIT NONE

    DOUBLE PRECISION:: lambda, lambda_plus_one

    lambda= SQRT(one + euclidean_inner_product(p,p))
    lambda_plus_one= lambda + one

    lambda_s(1)= one + (p(1)**2)  /lambda_plus_one
    lambda_s(2)=       (p(1)*p(2))/lambda_plus_one
    lambda_s(3)=       (p(1)*p(3))/lambda_plus_one
    lambda_s(4)= one + (p(2)**2)  /lambda_plus_one
    lambda_s(5)=       (p(2)*p(3))/lambda_plus_one
    lambda_s(6)= one + (p(3)**2)  /lambda_plus_one

    matrix(0,0)  = lambda
    matrix(0,1:3)= p
    matrix(1:3,0)= p
    matrix(1,1:3)= lambda_s(1:3)
    matrix(2,1:3)= [lambda_s(2), lambda_s(4), lambda_s(5)]
    matrix(3,1:3)= [lambda_s(3), lambda_s(5), lambda_s(6)]

  END PROCEDURE compute_boost_matrices


  !------------------------!
  !--  SPATIAL ROTATIONS --!
  !------------------------!


  MODULE PROCEDURE construct_rotation

    !***********************************************************
    !
    !# Construct a spatial rotation from a vector containing
    !  the Euler angles
    !
    !  FT 09.12.2022
    !
    !***********************************************************

    IMPLICIT NONE

    INTEGER:: i, j

    DOUBLE PRECISION:: v_sqnorm, delta
    DOUBLE PRECISION, DIMENSION(4,4):: identity(0:3,0:3)
    LOGICAL:: is_identity

    rotation% euler_angles= euler_angles

    CALL rotation% compute_rotation_matrices(rotation% euler_angles, &
                    rotation% r_x, rotation% r_y, rotation% r_z, &
                    rotation% r, rotation% matrix, &
                    rotation% inv_r_x, rotation% inv_r_y, rotation% inv_r_z, &
                    rotation% inv_r, rotation% inv_matrix)

    rotation% tr_r= TRANSPOSE(rotation% r)
    rotation% tr_matrix= TRANSPOSE(rotation% matrix)

    identity= square_matrix_multiplication(rotation% inv_matrix, &
                                           rotation% matrix)

    is_identity= .TRUE.
    DO i= 0, 3, 1
      DO j= 0, 3, 1
        IF( i == j )THEN
          delta= one
        ELSE
          delta= zero
        ENDIF
        is_identity= is_identity .AND. (identity(i,j) - delta) < 1.D-12
      ENDDO
    ENDDO
    IF( .NOT.is_identity )THEN

      PRINT *, "** ERROR! rotation% inv_matrix is not the inverse of ", &
               "rotation% matrix! Their product is:"
      PRINT *, "   id(0,:)=", identity(0,:)
      PRINT *, "   id(1,:)=", identity(1,:)
      PRINT *, "   id(2,:)=", identity(2,:)
      PRINT *, "   id(3,:)=", identity(3,:)
      PRINT *
      STOP

    ENDIF

    is_identity= .TRUE.
    DO i= 1, 3, 1
      DO j= 1, 3, 1
        is_identity= is_identity .AND. &
                     (rotation% inv_r(i,j) - rotation% tr_r(i,j)) < 1.D-12
      ENDDO
    ENDDO
    IF( .NOT.is_identity )THEN

      PRINT *, "** ERROR! rotation% inv_r is not the same as ", &
               "rotation% tr_r! "
      PRINT *
      STOP

    ENDIF

  END PROCEDURE construct_rotation


  MODULE PROCEDURE construct_rotation_angles

    !***********************************************************
    !
    !# Construct a spatial rotation from the Euler angles
    !
    !  FT 10.12.2022
    !
    !***********************************************************

    IMPLICIT NONE

    rotation= spatial_rotation([alpha,beta,gamma])

  END PROCEDURE construct_rotation_angles


  MODULE PROCEDURE compute_rotation_matrices

    !***********************************************************
    !
    !# Computes the spatial part of the matrix of the Lorentz
    !  rotation, and its whole matrix, starting from the vector
    !  \(p\)
    !
    !  FT 09.12.2022
    !
    !***********************************************************

    IMPLICIT NONE

    r_x(1,:)= [one, zero, zero]
    r_x(2,:)= [zero, COS(euler_angles(1)), - SIN(euler_angles(1))]
    r_x(3,:)= [zero, SIN(euler_angles(1)),   COS(euler_angles(1))]

    r_y(1,:)= [  COS(euler_angles(2)), zero, SIN(euler_angles(2))]
    r_y(2,:)= [  zero, one, zero]
    r_y(3,:)= [- SIN(euler_angles(2)), zero, COS(euler_angles(2))]

    r_z(1,:)= [COS(euler_angles(3)), - SIN(euler_angles(3)), zero]
    r_z(2,:)= [SIN(euler_angles(3)),   COS(euler_angles(3)), zero]
    r_z(3,:)= [zero, zero, one]

    r= square_matrix_multiplication(r_z,square_matrix_multiplication(r_y,r_x))

    matrix(0,0)= one
    matrix(0,1:3)= [zero,zero,zero]
    matrix(1:3,0)= [zero,zero,zero]
    matrix(1,1:3)= r(1,1:3)
    matrix(2,1:3)= r(2,1:3)
    matrix(3,1:3)= r(3,1:3)

    !
    !-- Computation of the inverse matrices
    !
    inv_r_x(1,:)= [one, zero, zero]
    inv_r_x(2,:)= [zero, COS(-euler_angles(1)), - SIN(-euler_angles(1))]
    inv_r_x(3,:)= [zero, SIN(-euler_angles(1)),   COS(-euler_angles(1))]

    inv_r_y(1,:)= [  COS(-euler_angles(2)), zero, SIN(-euler_angles(2))]
    inv_r_y(2,:)= [  zero, one, zero]
    inv_r_y(3,:)= [- SIN(-euler_angles(2)), zero, COS(-euler_angles(2))]

    inv_r_z(1,:)= [COS(-euler_angles(3)), - SIN(-euler_angles(3)), zero]
    inv_r_z(2,:)= [SIN(-euler_angles(3)),   COS(-euler_angles(3)), zero]
    inv_r_z(3,:)= [zero, zero, one]

    inv_r= square_matrix_multiplication(inv_r_x, &
                                square_matrix_multiplication(inv_r_y,inv_r_z))

    inv_matrix(0,0)= one
    inv_matrix(0,1:3)= [zero,zero,zero]
    inv_matrix(1:3,0)= [zero,zero,zero]
    inv_matrix(1,1:3)= inv_r(1,1:3)
    inv_matrix(2,1:3)= inv_r(2,1:3)
    inv_matrix(3,1:3)= inv_r(3,1:3)

  END PROCEDURE compute_rotation_matrices


  MODULE PROCEDURE get_lambda

    !***********************************************************
    !
    !# Returns the Lorentz factor [[lorentz_boost:lambda]]
    !
    !  FT 21.02.2023
    !
    !***********************************************************

    IMPLICIT NONE

    lambda= this% lambda

  END PROCEDURE get_lambda


END SUBMODULE constructors
