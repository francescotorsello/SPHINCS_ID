! File:         module_lorentz_group.f90
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


  TYPE lorentz_boost
  !!

    PRIVATE

    DOUBLE PRECISION, DIMENSION(3):: v
    !! Spatial velocity that determines the boost

    !DOUBLE PRECISION, DIMENSION(n_sym4x4):: g
    !! Spacetime Lorentzian metric

    DOUBLE PRECISION:: lambda
    !! Lorentz factor

    DOUBLE PRECISION, DIMENSION(3):: p
    !! Spatial vector equal to \(\lambda \,v\)

    DOUBLE PRECISION, DIMENSION(n_sym3x3):: lambda_s
    !! Spatial part of the Lorentz boost

    DOUBLE PRECISION, DIMENSION(4,4):: matrix(0:3,0:3)


    CONTAINS


    GENERIC, PUBLIC:: apply => apply_boost_to_vector, &
                               apply_boost_to_rank2_tensor!, &
                               !apply_boost_to_symrank2_tensor
    !# Generic procedure to apply the Lorentz boost to different geometric
    !  objects

    PROCEDURE:: apply_boost_to_vector
    !! Action of the [[lorentz_boost]] on a \(4\)-vector

    PROCEDURE:: apply_boost_to_rank2_tensor
    !! Action of the [[lorentz_boost]] on a \(4\)-vector

    !PROCEDURE:: apply_boost_to_symrank2_tensor
    !! Action of the [[lorentz_boost]] on a \(4\)-vector


  END TYPE lorentz_boost


  TYPE spatial_rotation
  !!

    PRIVATE

  END TYPE spatial_rotation


  TYPE lorentz_transformation
  !!

    PRIVATE

  END TYPE lorentz_transformation


  INTERFACE lorentz_boost

    MODULE FUNCTION construct_boost(v) RESULT(boost)

      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: v
      !! Spatial velocity that determines the boost
      TYPE(lorentz_boost):: boost
      !! [[lorentz_boost]] object to be constructed

    END FUNCTION construct_boost

  END INTERFACE lorentz_boost


  INTERFACE

    MODULE FUNCTION apply_boost_to_vector(this, u) RESULT(boosted_u)
    !! Action of the [[lorentz_boost]] on a \(4\)-vector

      CLASS(lorentz_boost), INTENT(IN):: this
      !! [[lorentz_boost]] object to apply
      DOUBLE PRECISION, DIMENSION(4), INTENT(IN):: u(0:3)
      !! \(4\)-vector to be boosted
      DOUBLE PRECISION, DIMENSION(4):: boosted_u(0:3)
      !! Boosted \(4\)-vector

    END FUNCTION apply_boost_to_vector

    MODULE FUNCTION apply_boost_to_rank2_tensor(this, t) RESULT(boosted_t)
    !! Action of the [[lorentz_boost]] on a \(4\)-vector

      CLASS(lorentz_boost), INTENT(IN):: this
      !! [[lorentz_boost]] object to apply
      DOUBLE PRECISION, DIMENSION(4,4), INTENT(IN):: t(0:3,0:3)
      !! \(4\)-vector to be boosted
      DOUBLE PRECISION, DIMENSION(4,4):: boosted_t(0:3,0:3)
      !! Boosted \(4\)-vector

    END FUNCTION apply_boost_to_rank2_tensor

  !  MODULE FUNCTION apply_boost_to_symrank2_tensor(this, t) RESULT(boosted_t)
  !  !! Action of the [[lorentz_boost]] on a \(4\)-vector
  !
  !    CLASS(lorentz_boost), INTENT(IN):: this
  !    !! [[lorentz_boost]] object to apply
  !    DOUBLE PRECISION, DIMENSION(10), INTENT(IN):: t
  !    !! \(4\)-vector to be boosted
  !    DOUBLE PRECISION, DIMENSION(4,4):: boosted_t(0:3,0:3)
  !    !! Boosted \(4\)-vector
  !
  !  END FUNCTION apply_boost_to_symrank2_tensor

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


END MODULE lorentz_group
