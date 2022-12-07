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

    DOUBLE PRECISION:: p
    !! Spatial vector equal to \(\lambda \,v\)

    DOUBLE PRECISION, DIMENSION(n_sym3x3):: lambda_s
    !! Spatial part of the Lorentz boost


    CONTAINS


    !GENERIC, PUBLIC:: apply => apply_to_vector, apply_to_rank2_tensor


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

    MODULE FUNCTION construct_boost( v ) RESULT( boost )

      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: v
      !! Spatial velocity that determines the boost
      TYPE(lorentz_boost):: boost
      !! [[lorentz_boost]] object to be constructed

    END FUNCTION construct_boost

  END INTERFACE lorentz_boost


END MODULE lorentz_group
