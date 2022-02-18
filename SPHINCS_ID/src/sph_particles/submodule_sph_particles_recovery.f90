! File:         submodule_sph_particles_recovery.f90
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

SUBMODULE (sph_particles) recovery

  !************************************************
  !
  !# This SUBMODULE contains the implementation
  !  of the method test_recovery of TYPE particles.
  !
  !  FT 18.02.2020
  !
  !************************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE test_recovery

    !************************************************
    !
    !# Tests the recovery. Computes the conserved
    !  variables from the physical ones, and then the
    !  physical ones from the conserved ones. It then
    !  compares the variables computed with the
    !  recovery PROCEDURES, with those computed with
    !  |sphincsid|. @todo add reference for recovery
    !
    !  FT 18.02.2020
    !
    !************************************************

    IMPLICIT NONE



  END PROCEDURE test_recovery


END SUBMODULE recovery
