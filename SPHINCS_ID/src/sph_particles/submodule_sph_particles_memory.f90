! File:         submodule_sph_particles_memory.f90
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

SUBMODULE (sph_particles) memory

  !***************************************************
  !
  !# This SUBMODULE contains the implementation of
  !  the methods of TYPE sph_particles
  !  that place particles on 1 or 2 lattices around
  !  the stars.
  !
  !  FT 12.07.2021
  !
  !***************************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE allocate_particles_memory

    !************************************************
    !
    !# Allocate memory for the LORENE ID on the
    !  particles
    !
    !  FT 10.11.2020
    !
    !************************************************

    IMPLICIT NONE

    PRINT *, "** Executing allocate_particles_memory."

    IF(.NOT.ALLOCATED( this% pos ))THEN
      ALLOCATE( this% pos( 3, this% npart ), STAT= ios, &
            ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array pos ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array pos" )
    ENDIF
    IF(.NOT.ALLOCATED( this% lapse ))THEN
      ALLOCATE( this% lapse( this% npart ), STAT= ios, &
            ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array lapse ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array lapse" )
    ENDIF
    IF(.NOT.ALLOCATED( this% shift_x ))THEN
      ALLOCATE( this% shift_x( this% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array shift_x ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for shift_x" )
    ENDIF
    IF(.NOT.ALLOCATED( this% shift_y ))THEN
      ALLOCATE( this% shift_y( this% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array shift_y ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for shift_y" )
    ENDIF
    IF(.NOT.ALLOCATED( this% shift_z ))THEN
      ALLOCATE( this% shift_z( this% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array shift_z ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for shift_z" )
    ENDIF
    IF(.NOT.ALLOCATED( this% g_xx ))THEN
      ALLOCATE( this% g_xx( this% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_xx ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array g_xx" )
    ENDIF
    IF(.NOT.ALLOCATED( this% g_xy ))THEN
      ALLOCATE( this% g_xy( this% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_xy ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for array g_xy" )
    ENDIF
    IF(.NOT.ALLOCATED( this% g_xz ))THEN
      ALLOCATE( this% g_xz( this% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_xz ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array g_xz" )
    ENDIF
    IF(.NOT.ALLOCATED( this% g_yy ))THEN
      ALLOCATE( this% g_yy( this% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_yy ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for array g_yy" )
    ENDIF
    IF(.NOT.ALLOCATED( this% g_yz ))THEN
      ALLOCATE( this% g_yz( this% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_yz ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array g_yz" )
    ENDIF
    IF(.NOT.ALLOCATED( this% g_zz ))THEN
      ALLOCATE( this% g_zz( this% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_zz ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array g_zz" )
    ENDIF
    IF(.NOT.ALLOCATED( this% baryon_density ))THEN
      ALLOCATE( this% baryon_density( this% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array baryon_density ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !     "...allocation error for array baryon_density" )
    ENDIF
    IF(.NOT.ALLOCATED( this% energy_density ))THEN
      ALLOCATE( this% energy_density( this% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array energy_density ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !             "...allocation error for array energy_density" )
    ENDIF
    IF(.NOT.ALLOCATED( this% specific_energy ))THEN
      ALLOCATE( this% specific_energy( this% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array specific_energy ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !    "...allocation error for array specific_energy" )
    ENDIF
    IF(.NOT.ALLOCATED( this% pressure ))THEN
      ALLOCATE( this% pressure( this% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array pressure ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array pressure" )
    ENDIF
    IF(.NOT.ALLOCATED( this% pressure_sph ))THEN
      ALLOCATE( this% pressure_sph( this% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array pressure_sph ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array pressure_sph" )
    ENDIF
    IF(.NOT.ALLOCATED( this% v_euler_x ))THEN
      ALLOCATE( this% v_euler_x( this% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array v_euler_x ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array v_euler_x" )
    ENDIF
    IF(.NOT.ALLOCATED( this% v_euler_y ))THEN
      ALLOCATE( this% v_euler_y( this% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array v_euler_y ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for array v_euler_y" )
    ENDIF
    IF(.NOT.ALLOCATED( this% v_euler_z ))THEN
      ALLOCATE( this% v_euler_z( this% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array v_euler_z ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for array v_euler_z" )
    ENDIF
    IF(.NOT.ALLOCATED( this% Theta ))THEN
      ALLOCATE( this% Theta( this% npart ), STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
      PRINT *, "...allocation error for array Theta in SUBROUTINE" &
               // "allocate_particles_memory. ", &
               "The error message is", err_msg
        STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( this% nstar ))THEN
        ALLOCATE( this% nstar( this% npart ), STAT= ios )
    ENDIF
    IF( ios > 0 )THEN
       PRINT *, '...allocation error for nstar'
       STOP
    ENDIF
    IF(.NOT.ALLOCATED( this% nstar_sph ))THEN
        ALLOCATE( this% nstar_sph( this% npart ), STAT= ios )
    ENDIF
    IF( ios > 0 )THEN
       PRINT *, '...allocation error for nstar_sph'
       STOP
    ENDIF
    IF(.NOT.ALLOCATED( this% particle_density ))THEN
        ALLOCATE( this% particle_density( this% npart ), STAT= ios )
    ENDIF
    IF( ios > 0 )THEN
       PRINT *, '...allocation error for particle_density'
       STOP
    ENDIF
    IF(.NOT.ALLOCATED( this% particle_density_sph ))THEN
        ALLOCATE( this% particle_density_sph( this% npart ), STAT= ios )
    ENDIF
    IF( ios > 0 )THEN
       PRINT *, '...allocation error for particle_density_sph'
       STOP
    ENDIF
    IF(.NOT.ALLOCATED( this% u_sph ))THEN
      ALLOCATE( this% u_sph( this% npart ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array u_sph in SUBROUTINE" &
                  // "allocate_particles_memory. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( this% nlrf ))THEN
      ALLOCATE( this% nlrf( this% npart ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array nlrf in SUBROUTINE" &
                  // "allocate_particles_memory. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( this% nlrf_sph ))THEN
      ALLOCATE( this% nlrf_sph( this% npart ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array nlrf_sph in SUBROUTINE" &
                  // "allocate_particles_memory. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( this% enthalpy ))THEN
      ALLOCATE( this% enthalpy( this% npart ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array enthalpy in SUBROUTINE" &
                  // "allocate_particles_memory. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( this% v ))THEN
      ALLOCATE( this% v( 0:3, this% npart ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array v in SUBROUTINE" &
                  // "allocate_particles_memory. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( this% Ye ))THEN
      ALLOCATE( this% Ye( this% npart ), STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
      PRINT *, "...allocation error for array Ye in SUBROUTINE" &
               // "allocate_particles_memory. ", &
               "The error message is", err_msg
        STOP
      ENDIF
    ENDIF

    PRINT *, "** Subroutine allocate_particles_memory executed."
    PRINT *

  END PROCEDURE allocate_particles_memory


  MODULE PROCEDURE deallocate_particles_memory

    !*************************************************
    !
    !# Deallocate memory for the LORENE ID on the
    !  particles
    !
    !  FT 12.07.2021 (this was part of the destructor
    !                 of TYPE [[particles]]
    !                 before this date)
    !
    !*************************************************

    IMPLICIT NONE

    IF( ALLOCATED( this% pos ))THEN
      DEALLOCATE( this% pos, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array pos. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% baryon_density ))THEN
      DEALLOCATE( this% baryon_density, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array baryon_density. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% energy_density ))THEN
      DEALLOCATE( this% energy_density, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array energy_density. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% specific_energy ))THEN
      DEALLOCATE( this% specific_energy, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array specific_energy. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% pressure ))THEN
      DEALLOCATE( this% pressure, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array pressure. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% pressure_sph ))THEN
      DEALLOCATE( this% pressure_sph, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array pressure_sph. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% v_euler_x ))THEN
      DEALLOCATE( this% v_euler_x, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array v_euler_x. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% v_euler_y ))THEN
      DEALLOCATE( this% v_euler_y, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array v_euler_y. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% v_euler_z ))THEN
      DEALLOCATE( this% v_euler_z, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array v_euler_z. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% lapse ))THEN
      DEALLOCATE( this% lapse, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array lapse. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% shift_x ))THEN
      DEALLOCATE( this% shift_x, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array shift_x. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% shift_y ))THEN
      DEALLOCATE( this% shift_y, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array shift_y. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% shift_z ))THEN
      DEALLOCATE( this% shift_z, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array shift_z. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% g_xx ))THEN
      DEALLOCATE( this% g_xx, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_xx. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% g_xy ))THEN
      DEALLOCATE( this% g_xy, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_xy. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% g_xz ))THEN
      DEALLOCATE( this% g_xz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_xz. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% g_yy ))THEN
      DEALLOCATE( this% g_yy, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_yy. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% g_yz ))THEN
      DEALLOCATE( this% g_yz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_yz. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% g_zz ))THEN
      DEALLOCATE( this% g_zz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_zz. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% nlrf ))THEN
      DEALLOCATE( this% nlrf, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array nlrf. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% nu ))THEN
      DEALLOCATE( this% nu, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array nu. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% Theta ))THEN
      DEALLOCATE( this% Theta, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array Theta. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( this% v ))THEN
      ALLOCATE( this% v( 0:3, this% npart ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array v ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% h ))THEN
      DEALLOCATE( this% h, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array h. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% v ))THEN
      DEALLOCATE( this% v, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array v. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% Ye ))THEN
      DEALLOCATE( this% Ye, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array Ye. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% nstar ))THEN
      DEALLOCATE( this% nstar, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array nstar. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% nstar_sph ))THEN
      DEALLOCATE( this% nstar_sph, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array nstar_sph. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% particle_density ))THEN
      DEALLOCATE( this% particle_density, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array particle_density. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% particle_density_sph ))THEN
      DEALLOCATE( this% particle_density_sph, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array particle_density_sph. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% u_sph ))THEN
      DEALLOCATE( this% u_sph, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array u_sph. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% nlrf_sph ))THEN
      DEALLOCATE( this% nlrf_sph, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array nlrf_sph. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% enthalpy ))THEN
      DEALLOCATE( this% enthalpy, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array enthalpy. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF( ALLOCATED( this% v ))THEN
      DEALLOCATE( this% v, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array v. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF

  END PROCEDURE deallocate_particles_memory


END SUBMODULE memory
