! File:         submodule_sph_particles_access.f90
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

SUBMODULE (sph_particles) access

  !**************************************************
  !
  !# This SUBMODULE contains the implementation of
  !  the methods of TYPE sph_particles
  !  that allow to access PRIVATE members.
  !
  !  FT 12.07.2021
  !
  !**************************************************


  IMPLICIT NONE


  CONTAINS


  !-----------------!
  !--  FUNCTIONS  --!
  !-----------------!


  MODULE PROCEDURE get_npart

    !************************************************
    !
    !# Returns the total number of particles
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    n_part= this% npart

  END PROCEDURE get_npart


  MODULE PROCEDURE get_npart_i

    !************************************************
    !
    !# Returns the number of particles on the
    !  \(i^{th}\) matter object
    !
    !  FT 10.11.2021
    !
    !************************************************

    IMPLICIT NONE

    n_part= this% npart_i(i_matter)

  END PROCEDURE get_npart_i


  MODULE PROCEDURE get_nuratio

    !************************************************
    !
    !# Returns the baryon number ratio on the stars
    !
    !  FT 27.04.2021
    !
    !************************************************

    IMPLICIT NONE

    nuratio= this% nuratio

  END PROCEDURE get_nuratio


  MODULE PROCEDURE get_nuratio_i

    !************************************************
    !
    !# Returns the baryon number ratio on the
    !  \(i^{th}\) matter object
    !
    !  FT 10.11.2021
    !
    !************************************************

    IMPLICIT NONE

    nuratio= this% nuratio_i(i_matter)

  END PROCEDURE get_nuratio_i


  MODULE PROCEDURE get_n_matter

    !************************************************
    !
    !# Returns the number of matter objects
    !
    !  FT 25.02.2021
    !
    !************************************************

    IMPLICIT NONE

    n_matter= this% n_matter

  END PROCEDURE get_n_matter


  MODULE PROCEDURE get_pos

    !************************************************
    !
    !# Returns the array of particle positions
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    pos_u= this% pos

  END PROCEDURE get_pos


  MODULE PROCEDURE get_vel

    !************************************************
    !
    !# Returns the array of coordinate 3-velocity of
    !  particles
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    vel= this% v(1:3,:)

  END PROCEDURE get_vel


  MODULE PROCEDURE get_nstar

    !************************************************
    !
    !# Returns the array of proper baryon density
    !  in the local rest frame
    !
    !  FT 25.02.2022
    !
    !************************************************

    IMPLICIT NONE

    nstar= this% nstar

  END PROCEDURE get_nstar


  MODULE PROCEDURE get_nstar_sph

    !************************************************
    !
    !# Returns the array of SPH density estimate of
    !  proper baryon density in the local rest frame
    !
    !  FT 25.02.2022
    !
    !************************************************

    IMPLICIT NONE

    nstar_sph= this% nstar_int

END PROCEDURE get_nstar_sph


  MODULE PROCEDURE get_nlrf

    !************************************************
    !
    !# Returns the array of baryon density in the
    !  local rest frame
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    nlrf= this% nlrf

  END PROCEDURE get_nlrf


  MODULE PROCEDURE get_nlrf_sph

    !************************************************
    !
    !# Returns the array of SPH estimate of the
    !  baryon density in the local rest frame
    !
    !  FT 25.02.2022
    !
    !************************************************

    IMPLICIT NONE

    nlrf_sph= this% nlrf_int

  END PROCEDURE get_nlrf_sph


  MODULE PROCEDURE get_nu

    !************************************************
    !
    !# Returns the array of baryon per particle
    ! [baryon (Msun_geo)^{-3}]
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    nu= this% nu

  END PROCEDURE get_nu


  MODULE PROCEDURE get_u

    !************************************************
    !
    !# Returns the array of specific internal
    !  energy [c^2]
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    u= this% specific_energy

  END PROCEDURE get_u


  MODULE PROCEDURE get_u_sph

    !************************************************
    !
    !# Returns the array of specific internal
    !  energy computed from the SPH estimate
    !  of the baryon number density in the local
    !  rest frame [c^2]
    !
    !  FT 25.02.2022
    !
    !************************************************

    IMPLICIT NONE

    u_sph= this% u_pwp

  END PROCEDURE get_u_sph


  MODULE PROCEDURE get_pressure

    !************************************************
    !
    !# Returns the array of pressure
    !  \([kg c^2 m^{-3}]\)
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    pressure= this% pressure

  END PROCEDURE get_pressure


  MODULE PROCEDURE get_pressure_cu

    !************************************************
    !
    !# Returns the array of pressure in code units
    ! [amu*c**2/(Msun_geo**3)]
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    pressure_cu= this% pressure_cu

  END PROCEDURE get_pressure_cu


  MODULE PROCEDURE get_theta

    !************************************************
    !
    !# Returns the array of generalized Lorentz
    ! factor
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    theta= this% Theta

  END PROCEDURE get_theta


  MODULE PROCEDURE get_h

    !************************************************
    !
    !# Returns the array of initial guess for the
    ! smoothing length [Msun_geo]
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    h= this% h

  END PROCEDURE get_h


  MODULE PROCEDURE get_lapse

    !************************************************
    !
    !# Returns the lapse function on the particles
    !
    !  FT 25.02.2022
    !
    !************************************************

    IMPLICIT NONE

    lapse= this% lapse

  END PROCEDURE get_lapse


  MODULE PROCEDURE get_shift

    !************************************************
    !
    !# Returns the shift vector on the particles
    !
    !  FT 25.02.2022
    !
    !************************************************

    IMPLICIT NONE

    shift(1,:)= this% shift_x
    shift(2,:)= this% shift_y
    shift(3,:)= this% shift_z

  END PROCEDURE get_shift


  MODULE PROCEDURE get_g3

    !************************************************
    !
    !# Returns the spatial metric on the particles
    !
    !  FT 25.02.2022
    !
    !************************************************

    IMPLICIT NONE

    g3(1,:)= this% g_xx
    g3(2,:)= this% g_xy
    g3(3,:)= this% g_xz
    g3(4,:)= this% g_yy
    g3(5,:)= this% g_yz
    g3(6,:)= this% g_zz

  END PROCEDURE get_g3


  MODULE PROCEDURE is_empty

    !************************************************
    !
    !# Returns the variable empty_object
    !  @warning experimental, not actively used in
    !          the code yet
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    answer= this% empty_object

  END PROCEDURE is_empty


END SUBMODULE access
