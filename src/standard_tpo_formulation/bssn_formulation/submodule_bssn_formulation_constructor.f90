! File:         submodule_bssn_formulation_constructor.f90
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

SUBMODULE (bssn_formulation) constructor

  !************************************************
  !
  !# Implementation of the constructor and
  !  destructor of TYPE [[bssn]]
  !
  !  FT 23.10.2020
  !
  !  Updated to support mesh refinement
  !
  !  FT 26.03.2021
  !
  !************************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE construct_bssn

    !****************************************************
    !
    !# This constructor of TYPE bssn calls the
    !  SUBROUTINES that rely on an bns object, and
    !  allocates memory. It constructs the grid
    !  using the number of grid points along each axis.
    !
    !  FT 23.10.2020
    !
    !****************************************************

    USE McLachlan_refine, ONLY: initialize_bssn
    USE BSSN_refine,      ONLY: deallocate_bssn
    USE mesh_refinement,  ONLY: levels, allocate_grid_function
    USE Extract_Mass,     ONLY: radius2
    USE utility,          ONLY: zero, one
    USE options,          ONLY: integ_meth

    IMPLICIT NONE

    ! integ_meth is a parameter set in SPHINCS_fm_input.dat,
    ! needed during the time integration. It is not needed at the level
    ! of the ID, but it must be set, otherwise the allocation of Ztmp
    ! does not work (see allocate_Ztmp). Ztmp s needed during the integration,
    ! so it should not be needed at the level of the ID. Nevertheless,
    ! some SUBROUTINES complain if it is not allocated (TODO: double check).
    ! When SPHINCS_ID is called with the parameter run_sph=.FALSE., the
    ! file SPHINCS_fm_input.dat is not read, hence integ_meth is not set.
    ! That is why it is set here.
    IF(integ_meth /= 1 .AND. integ_meth /= 2) integ_meth= 1

    ! Initialize the timer
    bssnid% bssn_computer_timer= timer( "bssn_computer_timer" )

    ! Construct the gravity grid and read the ID on it,
    ! in standard 3+1 formulation
    IF( PRESENT(dx) .AND. PRESENT(dy) .AND. PRESENT(dz) )THEN

      CALL bssnid% setup_standard_tpo_variables(id, dx, dy, dz)

    ELSE

      CALL bssnid% setup_standard_tpo_variables(id)

    ENDIF

    ! Read and store the BSSN parameters
    CALL initialize_bssn()
    CALL deallocate_bssn()
    DEALLOCATE(levels)

    IF( .NOT.ALLOCATED( bssnid% GC_int ))THEN
      ALLOCATE( bssnid% GC_int( bssnid% nlevels, 3 ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array GC_loo. ", &
                 "The error message is", err_msg
        STOP
      ENDIF
    ENDIF
    bssnid% GC_int= HUGE(one)

    IF( .NOT.ALLOCATED( bssnid% GC_parts_int ))THEN
      ALLOCATE( bssnid% GC_parts_int( bssnid% nlevels, 3 ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array GC_loo. ", &
                 "The error message is", err_msg
        STOP
      ENDIF
    ENDIF
    bssnid% GC_parts_int= HUGE(one)

    ! radius2 is the extraction radius. If not set here, then it is 0 by default
    ! and the metric is not interpolated on the particle in
    ! get_metric_on_particles
    radius2= HUGE(DBLE(1.0D0))

    PRINT *
    PRINT *, " * Ready to compute BSSN variables."
    PRINT *

  END PROCEDURE construct_bssn

  !
  !-- Keeping the following two SUBROUTINES separate in case it is needed
  !-- to add other PROCEDURES to the destructor (probably superfluous...)
  !
  MODULE PROCEDURE destruct_bssn

    !**************************************************
    !
    !# Destructor of the EXTENDED TYPE bssn
    !
    !  FT
    !
    !**************************************************

    IMPLICIT NONE

    CALL deallocate_bssn_fields(this)

  END PROCEDURE destruct_bssn


  MODULE PROCEDURE destructor

    !**************************************************
    !
    !# Destructor of EXTENDED TYPE bssn
    !
    !  FT
    !
    !**************************************************

    IMPLICIT NONE

    CALL this% destruct_bssn()

    CALL this% deallocate_standard_tpo_variables()

  END PROCEDURE destructor


END SUBMODULE constructor
