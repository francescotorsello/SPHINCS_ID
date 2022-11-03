! File:         submodule_bns_lorene_memory.f90
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

SUBMODULE (bns_lorene) memory

  !***********************************************
  !
  !# Implementation of the methods of TYPE [[bnslorene]]
  !  that (de)allocate memory
  !
  ! FT 9.07.2021
  !
  !***********************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE allocate_bnslorene_memory

    !***********************************************
    !
    !# Allocate the memory to store the \lorene| |id|
    !  in the member arrays
    !
    !  FT 17.09.2020
    !
    !***********************************************

    IMPLICIT NONE

    !PRINT *, "** Executing the allocate_bnslorene_memory subroutine..."

    IF(.NOT.ALLOCATED( this% lapse ))THEN
      ALLOCATE( this% lapse( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array lapse. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    !  CALL test_status( ios, err_msg, &
    !              "...allocation error for array lapse" )
    ENDIF
    IF(.NOT.ALLOCATED( this% shift_x ))THEN
      ALLOCATE( this% shift_x( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array shift_x. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array shift_x" )
    ENDIF
    IF(.NOT.ALLOCATED( this% shift_y ))THEN
      ALLOCATE( this% shift_y( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array shift_y. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array shift_y" )
    ENDIF
    IF(.NOT.ALLOCATED( this% shift_z ))THEN
      ALLOCATE( this% shift_z( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array shift_z. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array shift_z" )
    ENDIF
    IF(.NOT.ALLOCATED( this% g_xx ))THEN
      ALLOCATE( this% g_xx( d ), STAT= ios, &
          ERRMSG = err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array g_xx. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array g_xx" )
    ENDIF
    IF(.NOT.ALLOCATED( this% g_xy ))THEN
      ALLOCATE( this% g_xy( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array g_xy. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array g_xy" )
    ENDIF
    IF(.NOT.ALLOCATED( this% g_xz ))THEN
      ALLOCATE( this% g_xz( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array g_xz. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array g_xz" )
    ENDIF
    IF(.NOT.ALLOCATED( this% g_yy ))THEN
      ALLOCATE( this% g_yy( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array g_yy. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array g_yy" )
    ENDIF
    IF(.NOT.ALLOCATED( this% g_yz ))THEN
      ALLOCATE( this% g_yz( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array g_yz. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array g_yz" )
    ENDIF
    IF(.NOT.ALLOCATED( this% g_zz ))THEN
      ALLOCATE( this% g_zz( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array g_zz. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array g_zz" )
    ENDIF
    IF(.NOT.ALLOCATED( this% k_xx ))THEN
      ALLOCATE( this% k_xx( d ), STAT= ios, &

          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array k_xx. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array k_xx" )
    ENDIF
    IF(.NOT.ALLOCATED( this% k_xy ))THEN
      ALLOCATE( this% k_xy( d ), STAT= ios, &

          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array k_xy. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array k_xy" )
    ENDIF
    IF(.NOT.ALLOCATED( this% k_xz ))THEN
      ALLOCATE( this% k_xz( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array k_xz. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array k_xz" )
    ENDIF
    IF(.NOT.ALLOCATED( this% k_yy ))THEN
      ALLOCATE( this%  k_yy( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array k_yy ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array k_yy" )
    ENDIF
    IF(.NOT.ALLOCATED( this% k_yz ))THEN
      ALLOCATE( this% k_yz( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array k_yz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array k_yz" )
    ENDIF
    IF(.NOT.ALLOCATED( this% k_zz ))THEN
      ALLOCATE( this% k_zz( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array k_zz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array k_zz" )
    ENDIF
    IF(.NOT.ALLOCATED( this% baryon_density ))THEN
      ALLOCATE( this% baryon_density( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array baryon_density ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array baryon_density" )
    ENDIF
    IF(.NOT.ALLOCATED( this% energy_density ))THEN
      ALLOCATE( this% energy_density( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array energy_density ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array energy_density" )
    ENDIF
    IF(.NOT.ALLOCATED( this% specific_energy ))THEN
      ALLOCATE( this% specific_energy( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array specific_energy ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array specific_energy" )
    ENDIF
    IF(.NOT.ALLOCATED( this% v_euler_x ))THEN
      ALLOCATE( this% v_euler_x( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array v_euler_x ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array v_euler_x" )
    ENDIF
    IF(.NOT.ALLOCATED( this% v_euler_y ))THEN
      ALLOCATE( this% v_euler_y( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array v_euler_y ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array v_euler_y" )
    ENDIF
    IF(.NOT.ALLOCATED( this% v_euler_z ))THEN
      ALLOCATE( this% v_euler_z( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array v_euler_z ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array v_euler_z" )
    ENDIF

    IF( SIZE( this% lapse ) /= d )THEN
      PRINT *, "** ERROR in memory allocation in allocate_bnslorene_memory"
    ENDIF

    !PRINT *, "** Subroutine allocate_bnslorene_memory executed."
    !PRINT *

  END PROCEDURE allocate_bnslorene_memory


  MODULE PROCEDURE deallocate_bnslorene_memory

    !***********************************************
    !
    !# Deallocate the memory for the member arrays
    !
    !  FT 17.09.2020
    !
    !***********************************************

    IMPLICIT NONE

    !PRINT *, "** Executing the deallocate_bnslorene_memory subroutine..."

    IF(ALLOCATED( this% lapse ))THEN
      DEALLOCATE( this% lapse, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array lapse ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                  "...deallocation error for array lapse" )
    ENDIF
    IF(ALLOCATED( this% shift_x ))THEN
      DEALLOCATE( this% shift_x, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array shift_x ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array shift_x" )
    ENDIF
    IF(ALLOCATED( this% shift_y ))THEN
      DEALLOCATE( this% shift_y, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array shift_y ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array shift_y" )
    ENDIF
    IF(ALLOCATED( this% shift_z ))THEN
      DEALLOCATE( this% shift_z, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array shift_z ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array shift_z" )
    ENDIF
    IF(ALLOCATED( this% g_xx ))THEN
      DEALLOCATE( this% g_xx, STAT= ios, ERRMSG = err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_xx ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_xx" )
    ENDIF
    IF(ALLOCATED( this% g_xy ))THEN
      DEALLOCATE( this% g_xy, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_xy ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !               "...deallocation error for array g_xy" )
    ENDIF
    IF(ALLOCATED( this% g_xz ))THEN
      DEALLOCATE( this% g_xz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_xz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_xz" )
    ENDIF
    IF(ALLOCATED( this% g_yy ))THEN
      DEALLOCATE( this% g_yy, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_yy ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_yy" )
    ENDIF
    IF(ALLOCATED( this% g_yz ))THEN
      DEALLOCATE( this% g_yz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_yz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_yz" )
    ENDIF
    IF(ALLOCATED( this% g_zz ))THEN
      DEALLOCATE( this% g_zz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_zz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_zz" )
    ENDIF
    IF(ALLOCATED( this% k_xx ))THEN
      DEALLOCATE( this% k_xx, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array k_xx ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array k_xx" )
    ENDIF
    IF(ALLOCATED( this% k_xy ))THEN
      DEALLOCATE( this% k_xy, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array k_xy ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array k_xy" )
    ENDIF
    IF(ALLOCATED( this% k_xz ))THEN
      DEALLOCATE( this% k_xz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array k_xz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array k_xz" )
    ENDIF
    IF(ALLOCATED( this% k_yy ))THEN
      DEALLOCATE( this% k_yy, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array k_yy ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array k_yy" )
    ENDIF
    IF(ALLOCATED( this% k_yz ))THEN
      DEALLOCATE( this% k_yz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array k_yz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array k_yz" )
    ENDIF
    IF(ALLOCATED( this% k_zz ))THEN
      DEALLOCATE( this% k_zz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array k_zz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array k_zz" )
    ENDIF
    IF(ALLOCATED( this% baryon_density ))THEN
      DEALLOCATE( this% baryon_density, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array baryon_density ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...deallocation error for array baryon_density" )
    ENDIF
    IF(ALLOCATED( this% energy_density ))THEN
      DEALLOCATE( this% energy_density, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array energy_density ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...deallocation error for array energy_density" )
    ENDIF
    IF(ALLOCATED( this% specific_energy ))THEN
      DEALLOCATE( this% specific_energy, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array specific_energy ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...deallocation error for array specific_energy" )
    ENDIF
    IF(ALLOCATED( this% v_euler_x ))THEN
      DEALLOCATE( this% v_euler_x, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array v_euler_x ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...deallocation error for array v_euler_x" )
    ENDIF
    IF(ALLOCATED( this% v_euler_y ))THEN
      DEALLOCATE( this% v_euler_y, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array v_euler_y ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v_euler_y" )
    ENDIF
    IF(ALLOCATED( this% v_euler_z ))THEN
      DEALLOCATE( this% v_euler_z, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array v_euler_z ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v_euler_z" )
    ENDIF

    !PRINT *, "** Subroutine deallocate_bnslorene_memory executed."
    !PRINT *

  END PROCEDURE deallocate_bnslorene_memory


END SUBMODULE memory
