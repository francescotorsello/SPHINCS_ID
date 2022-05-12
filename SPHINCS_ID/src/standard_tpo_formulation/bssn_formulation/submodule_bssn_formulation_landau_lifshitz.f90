! File:         submodule_bssn_formulation_landau_lifshitz.f90
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

SUBMODULE (bssn_formulation) landau_lifshitz

  !************************************************
  !
  !# Implementation of the method of TYPE bssn
  !  that computes the Ricci tensor and scalar
  !
  !  FT 10.02.2022
  !
  !************************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  SUBROUTINE compute_landau_lifshitz_pseudotensor()

    !************************************************
    !
    !# Computes the Landau-Lifshitz pseudotensor
    !  on the mesh, using the definition in eq.(96.9)
    !  in Landau and Lifshitz, The Classical Theory
    !  of Fields, Pergamon Press(1975)
    !
    !  FT 11.05.2022
    !
    !************************************************



  END SUBROUTINE compute_landau_lifshitz_pseudotensor


END SUBMODULE landau_lifshitz
