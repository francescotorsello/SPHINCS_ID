! File:         submodule_bns_fuka_io.f90
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

SUBMODULE (bns_fuka) io

  !********************************************
  !
  !# This submodule contains the implementation of the
  !  methods of TYPE [[bnsfuka]] that handle I/O (input/output)
  !
  !  FT 05.11.2021
  !
  !********************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE print_summary_bnsfuka

    !************************************************
    !
    !# Prints a summary of the physical properties the |bns| system
    !  produced by |fuka| to the standard output and, optionally,
    !  to a formatted file whose name is given as the optional
    !  argument `filename`
    !
    !  FT 4.02.2022
    !
    !************************************************

    IMPLICIT NONE

  !  PRINT *, "   * Binary system of neutron stars produced by FUKA:"
  !  PRINT *
  !  PRINT *, "   ADM linear momentum of the system=(", this% linear_momentum_x,&
  !           ", "
  !  PRINT *, "                                      ", this% linear_momentum_y,&
  !           ", "
  !  PRINT *, "                                      ", this% linear_momentum_z,&
  !           ") Msun*c"
  !  PRINT *
  !  PRINT *, "   Bowen-York angular momentum of the system= (", &
  !           this% angular_momentum_x, &
  !           ", "
  !  PRINT *, "                                               ", &
  !           this% angular_momentum_y, &
  !           ", "
  !  PRINT *, "                                               ", &
  !           this% angular_momentum_z, ") G*Msun^2/c"
  !  PRINT *


  END PROCEDURE print_summary_bnsfuka


  MODULE PROCEDURE print_bns_properties

    !****************************************************
    !
    !# Print the parameters of the binary neutron
    !  stars' initial data computed by |fuka|
    !
    !  FT 4.02.2022
    !
    !****************************************************


  END PROCEDURE print_bns_properties


END SUBMODULE io
