! File:         submodule_lorentz_group_constructors.f90
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


  IMPLICIT NONE



  CONTAINS



  MODULE PROCEDURE construct_boost

    !***********************************************************
    !
    !# Construct a boost transformation
    !
    !  FT 08.12.2022
    !
    !***********************************************************

    USE utility,  ONLY: one

    IMPLICIT NONE

    DOUBLE PRECISION:: v_sqnorm, lambda_plus_one

    v_sqnorm= euclidean_inner_product(v,v)

    IF( v_sqnorm >= one )THEN
      PRINT *, "** ERROR! The speed given to construct_boost is larger", &
               " than, or equal to, 1!"
      PRINT *, "|v|^2=", v_sqnorm
      PRINT *, "* Stopping..."
      PRINT *
    ENDIF

    boost% v= v

    boost% lambda= one/SQRT( one - v_sqnorm )
    lambda_plus_one= boost% lambda + one

    boost% p= boost% lambda*boost% v

    boost% lambda_s(1)= one + (boost% p(1)**2)         /lambda_plus_one
    boost% lambda_s(2)=       (boost% p(1)*boost% p(2))/lambda_plus_one
    boost% lambda_s(3)=       (boost% p(1)*boost% p(3))/lambda_plus_one
    boost% lambda_s(4)= one + (boost% p(2)**2)         /lambda_plus_one
    boost% lambda_s(5)=       (boost% p(2)*boost% p(3))/lambda_plus_one
    boost% lambda_s(6)= one + (boost% p(3)**2)         /lambda_plus_one

    boost% matrix(0,0)  = boost% lambda
    boost% matrix(0,1:3)= boost% p
    boost% matrix(1:3,0)= boost% p
    boost% matrix(1,1:3)= boost% lambda_s(1:3)
    boost% matrix(2,1:3)= &
      [boost% lambda_s(2),boost% lambda_s(4),boost% lambda_s(5)]
    boost% matrix(3,1:3)= &
      [boost% lambda_s(3),boost% lambda_s(5),boost% lambda_s(6)]

  END PROCEDURE construct_boost


END SUBMODULE constructors
