! File:         submodule_sph_particles_memory.f90
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


  MODULE PROCEDURE allocate_lorene_id_parts_memory

    !************************************************
    !
    !# Allocate memory for the LORENE ID on the
    !  particles
    !
    !  FT 10.11.2020
    !
    !************************************************

    IMPLICIT NONE

    PRINT *, "** Executing allocate_lorene_id_parts_memory."

    IF(.NOT.ALLOCATED( THIS% pos ))THEN
      ALLOCATE( THIS% pos( 3, THIS% npart ), STAT= ios, &
            ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array pos ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array pos" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% lapse ))THEN
      ALLOCATE( THIS% lapse( THIS% npart ), STAT= ios, &
            ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array lapse ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array lapse" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% shift_x ))THEN
      ALLOCATE( THIS% shift_x( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array shift_x ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for shift_x" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% shift_y ))THEN
      ALLOCATE( THIS% shift_y( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array shift_y ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for shift_y" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% shift_z ))THEN
      ALLOCATE( THIS% shift_z( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array shift_z ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for shift_z" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_xx ))THEN
      ALLOCATE( THIS% g_xx( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_xx ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array g_xx" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_xy ))THEN
      ALLOCATE( THIS% g_xy( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_xy ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for array g_xy" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_xz ))THEN
      ALLOCATE( THIS% g_xz( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_xz ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array g_xz" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_yy ))THEN
      ALLOCATE( THIS% g_yy( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_yy ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for array g_yy" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_yz ))THEN
      ALLOCATE( THIS% g_yz( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_yz ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array g_yz" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_zz ))THEN
      ALLOCATE( THIS% g_zz( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_zz ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array g_zz" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% baryon_density ))THEN
      ALLOCATE( THIS% baryon_density( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array baryon_density ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !     "...allocation error for array baryon_density" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% energy_density ))THEN
      ALLOCATE( THIS% energy_density( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array energy_density ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !             "...allocation error for array energy_density" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% specific_energy ))THEN
      ALLOCATE( THIS% specific_energy( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array specific_energy ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !    "...allocation error for array specific_energy" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% pressure ))THEN
      ALLOCATE( THIS% pressure( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array pressure ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array pressure" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% pressure_cu ))THEN
      ALLOCATE( THIS% pressure_cu( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array pressure_cu ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array pressure_cu" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% v_euler_x ))THEN
      ALLOCATE( THIS% v_euler_x( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array v_euler_x ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array v_euler_x" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% v_euler_y ))THEN
      ALLOCATE( THIS% v_euler_y( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array v_euler_y ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for array v_euler_y" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% v_euler_z ))THEN
      ALLOCATE( THIS% v_euler_z( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array v_euler_z ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for array v_euler_z" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% nstar ))THEN
        ALLOCATE( THIS% nstar( THIS% npart ), STAT= ios )
    ENDIF
    IF( ios > 0 )THEN
       PRINT *, '...allocation error for nstar'
       STOP
    ENDIF
    IF(.NOT.ALLOCATED( THIS% nstar_int ))THEN
        ALLOCATE( THIS% nstar_int( THIS% npart ), STAT= ios )
    ENDIF
    IF( ios > 0 )THEN
       PRINT *, '...allocation error for nstar_int'
       STOP
    ENDIF
    IF(.NOT.ALLOCATED( THIS% particle_density ))THEN
        ALLOCATE( THIS% particle_density( THIS% npart ), STAT= ios )
    ENDIF
    IF( ios > 0 )THEN
       PRINT *, '...allocation error for particle_density'
       STOP
    ENDIF
    IF(.NOT.ALLOCATED( THIS% particle_density_int ))THEN
        ALLOCATE( THIS% particle_density_int( THIS% npart ), STAT= ios )
    ENDIF
    IF( ios > 0 )THEN
       PRINT *, '...allocation error for particle_density_int'
       STOP
    ENDIF
    IF(.NOT.ALLOCATED( THIS% pmass ))THEN
      ALLOCATE( THIS% pmass( THIS% npart ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array pmass in SUBROUTINE" &
                  // " allocate_lorene_id_memory. ", &
                  "The STAT variable is", ios, ". ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( THIS% u_pwp ))THEN
      ALLOCATE( THIS% u_pwp( THIS% npart ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array u_pwp in SUBROUTINE" &
                  // "allocate_lorene_id_memory. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( THIS% nlrf_int ))THEN
      ALLOCATE( THIS% nlrf_int( THIS% npart ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array nlrf_int in SUBROUTINE" &
                  // "allocate_lorene_id_memory. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF

    PRINT *, "** Subroutine allocate_lorene_id_memory executed."
    PRINT *

  END PROCEDURE allocate_lorene_id_parts_memory


  MODULE PROCEDURE deallocate_lorene_id_parts_memory

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

    IF( ALLOCATED( THIS% pos ))THEN
      DEALLOCATE( THIS% pos, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array pos. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array pos in SUBROUTINE"&
      !                // "destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% baryon_density ))THEN
      DEALLOCATE( THIS% baryon_density, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array baryon_density. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array " &
      !                // "baryon_density in SUBROUTINE " &
      !                // "destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% energy_density ))THEN
      DEALLOCATE( THIS% energy_density, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array energy_density. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array " &
      !                // "energy_density in SUBROUTINE " &
      !                // "destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% specific_energy ))THEN
      DEALLOCATE( THIS% specific_energy, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array specific_energy. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array " &
      !                // "specific_energy in SUBROUTINE " &
      !                // "destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% pressure ))THEN
      DEALLOCATE( THIS% pressure, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array pressure. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array " &
      !                // "pressure in SUBROUTINE " &
      !                // "destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% pressure_cu ))THEN
      DEALLOCATE( THIS% pressure_cu, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array pressure_cu. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array " &
      !                // "pressure_cu in SUBROUTINE " &
      !                // "destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% v_euler_x ))THEN
      DEALLOCATE( THIS% v_euler_x, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array v_euler_x. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array " &
      !                // "v_euler_x in SUBROUTINE " &
      !                // "destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% v_euler_y ))THEN
      DEALLOCATE( THIS% v_euler_y, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array v_euler_y. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array " &
      !                // "v_euler_y in SUBROUTINE " &
      !                // "destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% v_euler_z ))THEN
      DEALLOCATE( THIS% v_euler_z, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array v_euler_z. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array " &
      !                // "v_euler_z in SUBROUTINE " &
      !                // "destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% lapse ))THEN
      DEALLOCATE( THIS% lapse, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array lapse. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array lapse in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% shift_x ))THEN
      DEALLOCATE( THIS% shift_x, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array shift_x. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array shift_x in "&
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% shift_y ))THEN
      DEALLOCATE( THIS% shift_y, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array shift_y. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array shift_y in "&
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% shift_z ))THEN
      DEALLOCATE( THIS% shift_z, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array shift_z. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array shift_z in "&
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% g_xx ))THEN
      DEALLOCATE( THIS% g_xx, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_xx. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_xx in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% g_xy ))THEN
      DEALLOCATE( THIS% g_xy, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_xy. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_xy in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% g_xz ))THEN
      DEALLOCATE( THIS% g_xz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_xz. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_xz in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% g_yy ))THEN
      DEALLOCATE( THIS% g_yy, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_yy. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_yy in " &
      !                // "SUBROUTINE estruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% g_yz ))THEN
      DEALLOCATE( THIS% g_yz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_yz. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_yz in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% g_zz ))THEN
      DEALLOCATE( THIS% g_zz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_zz. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_zz in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% nlrf ))THEN
      DEALLOCATE( THIS% nlrf, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array nlrf. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array nlrf in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% nu ))THEN
      DEALLOCATE( THIS% nu, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array nu. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array nu in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% Theta ))THEN
      DEALLOCATE( THIS% Theta, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array Theta. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array Theta in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% h ))THEN
      DEALLOCATE( THIS% h, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array h. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array h in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% v ))THEN
      DEALLOCATE( THIS% v, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array v. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% Ye ))THEN
      DEALLOCATE( THIS% Ye, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array Ye. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% nstar ))THEN
      DEALLOCATE( THIS% nstar, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array nstar. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% nstar_int ))THEN
      DEALLOCATE( THIS% nstar_int, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array nstar_int. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% particle_density ))THEN
      DEALLOCATE( THIS% particle_density, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array particle_density. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% particle_density_int ))THEN
      DEALLOCATE( THIS% particle_density_int, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array particle_density_int. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% pmass ))THEN
      DEALLOCATE( THIS% pmass, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array pmass. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% u_pwp ))THEN
      DEALLOCATE( THIS% u_pwp, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array u_pwp. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% nlrf_int ))THEN
      DEALLOCATE( THIS% nlrf_int, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array u_pwp. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF

  END PROCEDURE deallocate_lorene_id_parts_memory


END SUBMODULE memory
