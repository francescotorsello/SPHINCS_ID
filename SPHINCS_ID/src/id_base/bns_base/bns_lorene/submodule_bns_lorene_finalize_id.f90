! File:         submodule_bnslorene_finalize_id.f90
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

SUBMODULE (bns_lorene) finalize_id

  !*********************************************************
  !
  !# Implementation of the PROCEDURE that act on the |id|
  !  after it is set up on the mesh and/or on the particles,
  !  to finalize its preparation
  !
  !  FT 14.04.2022
  !
  !*********************************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE correct_adm_linear_momentum

    !***********************************************
    !
    !# Correct the velocity and the generalized
    !  Lorentz factor, so that the \(\mathrm{ADM}\)
    !  linear momentum of the |bns| is 0
    !
    !  FT 14.04.2022
    !
    !***********************************************

    IMPLICIT NONE



  END PROCEDURE correct_adm_linear_momentum


END SUBMODULE finalize_id
