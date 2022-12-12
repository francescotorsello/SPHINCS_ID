! File:         construct_eccentric_binary.f90
! Author:       Francesco Torsello (FT)
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

PROGRAM construct_eccentric_binary

  !*****************************************************
  !
  !# Read two |tov| star BSSN and SPH ID, and construct
  !  an approximate eccentric binary system with a
  !  velocity determined with Newtonian mechanics and
  !  Newtonian gravitational law
  !
  !  FT 12.12.2022
  !
  !*****************************************************

  USE sph_variables,       ONLY: npart, &  ! particle number
                                 pos_u, &  ! particle positions
                                 vel_u, &  ! particle velocities in
                                           ! coordinate frame
                                 nlrf,  &  ! baryon number density in
                                           ! local rest frame
                                 !ehat,  &  ! canonical energy per baryon
                                 nu,    &  ! canonical baryon number per
                                           ! particle
                                 Theta, &  ! Generalized Lorentz factor
                                 h,     &  ! Smoothing length
                                 Pr,    &  ! Pressure
                                 u,     &  ! Internal energy in local rest
                                           ! frame (no kinetic energy)
                                 temp,  &  ! Temperature
                                 av,    &  ! Dissipation
                                 ye,    &  ! Electron fraction
                                 divv,  &  ! Divergence of velocity vel_u
                                 allocate_SPH_memory, &
                                 deallocate_SPH_memory, &
                                 n1, n2
  USE input_output, ONLY: read_sphincs_dump, write_sphincs_dump

  IMPLICIT NONE

  INTEGER, PARAMETER:: max_npart= 5.D+7

  CHARACTER(LEN=:), ALLOCATABLE:: filename

  npart= max_npart

  CALL set_units('NSM')

  CALL allocate_SPH_memory

  filename= 'tov-id-files/NSxx.00000'

  CALL read_sphincs_dump(filename)

END PROGRAM construct_eccentric_binary
