! File:         submodule_diffstar_lorene_memory.f90
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

SUBMODULE (diffstar_lorene) memory

  !***********************************************
  !
  !# Implementation of the methods of TYPE diffstar
  !  that (de)allocate memory
  !
  ! FT 25.10.2021
  !
  !***********************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE allocate_diffstar_memory

    !***********************************************
    !
    !# Allocate the memory to store the LORENE ID
    !  in the member arrays
    !
    !  FT 25.10.2021
    !
    !***********************************************

    IMPLICIT NONE

    !PRINT *, "** Executing the allocate_diffstar_memory subroutine..."

    IF(.NOT.ALLOCATED( THIS% lapse ))THEN
      ALLOCATE( THIS% lapse( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array lapse. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    !  CALL test_status( ios, err_msg, &
    !              "...allocation error for array lapse" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% shift_x ))THEN
      ALLOCATE( THIS% shift_x( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array shift_x. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array shift_x" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% shift_y ))THEN
      ALLOCATE( THIS% shift_y( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array shift_y. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array shift_y" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% shift_z ))THEN
      ALLOCATE( THIS% shift_z( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array shift_z. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array shift_z" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_xx ))THEN
      ALLOCATE( THIS% g_xx( d ), STAT= ios, &
          ERRMSG = err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array g_xx. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array g_xx" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_xy ))THEN
      ALLOCATE( THIS% g_xy( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array g_xy. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array g_xy" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_xz ))THEN
      ALLOCATE( THIS% g_xz( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array g_xz. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array g_xz" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_yy ))THEN
      ALLOCATE( THIS% g_yy( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array g_yy. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array g_yy" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_yz ))THEN
      ALLOCATE( THIS% g_yz( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array g_yz. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array g_yz" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_zz ))THEN
      ALLOCATE( THIS% g_zz( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array g_zz. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array g_zz" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% k_xx ))THEN
      ALLOCATE( THIS% k_xx( d ), STAT= ios, &

          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array k_xx. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array k_xx" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% k_xy ))THEN
      ALLOCATE( THIS% k_xy( d ), STAT= ios, &

          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array k_xy. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array k_xy" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% k_xz ))THEN
      ALLOCATE( THIS% k_xz( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array k_xz. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array k_xz" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% k_yy ))THEN
      ALLOCATE( THIS%  k_yy( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array k_yy ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array k_yy" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% k_yz ))THEN
      ALLOCATE( THIS% k_yz( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array k_yz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array k_yz" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% k_zz ))THEN
      ALLOCATE( THIS% k_zz( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array k_zz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array k_zz" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% baryon_density ))THEN
      ALLOCATE( THIS% baryon_density( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array baryon_density ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array baryon_density" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% energy_density ))THEN
      ALLOCATE( THIS% energy_density( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array energy_density ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array energy_density" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% specific_energy ))THEN
      ALLOCATE( THIS% specific_energy( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array specific_energy ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array specific_energy" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% v_euler_x ))THEN
      ALLOCATE( THIS% v_euler_x( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array v_euler_x ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array v_euler_x" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% v_euler_y ))THEN
      ALLOCATE( THIS% v_euler_y( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array v_euler_y ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array v_euler_y" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% v_euler_z ))THEN
      ALLOCATE( THIS% v_euler_z( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array v_euler_z ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array v_euler_z" )
    ENDIF

    IF( SIZE( THIS% lapse ) /= d )THEN
      PRINT *, "** ERROR in memory allocation in allocate_diffstar_memory"
    ENDIF

    !PRINT *, "** Subroutine allocate_diffstar_memory executed."
    !PRINT *

  END PROCEDURE allocate_diffstar_memory


  MODULE PROCEDURE deallocate_diffstar_memory

    !***********************************************
    !
    !# Deallocate the memory for the member arrays
    !
    !  FT 25.10.2021
    !
    !***********************************************

    IMPLICIT NONE

    !PRINT *, "** Executing the deallocate_diffstar_memory subroutine..."

    IF(ALLOCATED( THIS% lapse ))THEN
      DEALLOCATE( THIS% lapse, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array lapse ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                  "...deallocation error for array lapse" )
    ENDIF
    IF(ALLOCATED( THIS% shift_x ))THEN
      DEALLOCATE( THIS% shift_x, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array shift_x ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array shift_x" )
    ENDIF
    IF(ALLOCATED( THIS% shift_y ))THEN
      DEALLOCATE( THIS% shift_y, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array shift_y ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array shift_y" )
    ENDIF
    IF(ALLOCATED( THIS% shift_z ))THEN
      DEALLOCATE( THIS% shift_z, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array shift_z ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array shift_z" )
    ENDIF
    IF(ALLOCATED( THIS% g_xx ))THEN
      DEALLOCATE( THIS% g_xx, STAT= ios, ERRMSG = err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_xx ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_xx" )
    ENDIF
    IF(ALLOCATED( THIS% g_xy ))THEN
      DEALLOCATE( THIS% g_xy, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_xy ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !               "...deallocation error for array g_xy" )
    ENDIF
    IF(ALLOCATED( THIS% g_xz ))THEN
      DEALLOCATE( THIS% g_xz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_xz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_xz" )
    ENDIF
    IF(ALLOCATED( THIS% g_yy ))THEN
      DEALLOCATE( THIS% g_yy, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_yy ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_yy" )
    ENDIF
    IF(ALLOCATED( THIS% g_yz ))THEN
      DEALLOCATE( THIS% g_yz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_yz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_yz" )
    ENDIF
    IF(ALLOCATED( THIS% g_zz ))THEN
      DEALLOCATE( THIS% g_zz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_zz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_zz" )
    ENDIF
    IF(ALLOCATED( THIS% k_xx ))THEN
      DEALLOCATE( THIS% k_xx, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array k_xx ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array k_xx" )
    ENDIF
    IF(ALLOCATED( THIS% k_xy ))THEN
      DEALLOCATE( THIS% k_xy, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array k_xy ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array k_xy" )
    ENDIF
    IF(ALLOCATED( THIS% k_xz ))THEN
      DEALLOCATE( THIS% k_xz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array k_xz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array k_xz" )
    ENDIF
    IF(ALLOCATED( THIS% k_yy ))THEN
      DEALLOCATE( THIS% k_yy, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array k_yy ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array k_yy" )
    ENDIF
    IF(ALLOCATED( THIS% k_yz ))THEN
      DEALLOCATE( THIS% k_yz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array k_yz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array k_yz" )
    ENDIF
    IF(ALLOCATED( THIS% k_zz ))THEN
      DEALLOCATE( THIS% k_zz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array k_zz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array k_zz" )
    ENDIF
    IF(ALLOCATED( THIS% baryon_density ))THEN
      DEALLOCATE( THIS% baryon_density, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array baryon_density ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...deallocation error for array baryon_density" )
    ENDIF
    IF(ALLOCATED( THIS% energy_density ))THEN
      DEALLOCATE( THIS% energy_density, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array energy_density ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...deallocation error for array energy_density" )
    ENDIF
    IF(ALLOCATED( THIS% specific_energy ))THEN
      DEALLOCATE( THIS% specific_energy, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array specific_energy ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...deallocation error for array specific_energy" )
    ENDIF
    IF(ALLOCATED( THIS% v_euler_x ))THEN
      DEALLOCATE( THIS% v_euler_x, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array v_euler_x ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...deallocation error for array v_euler_x" )
    ENDIF
    IF(ALLOCATED( THIS% v_euler_y ))THEN
      DEALLOCATE( THIS% v_euler_y, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array v_euler_y ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v_euler_y" )
    ENDIF
    IF(ALLOCATED( THIS% v_euler_z ))THEN
      DEALLOCATE( THIS% v_euler_z, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array v_euler_z ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v_euler_z" )
    ENDIF

    !PRINT *, "** Subroutine deallocate_diffstar_memory executed."
    !PRINT *

  END PROCEDURE deallocate_diffstar_memory


END SUBMODULE memory
