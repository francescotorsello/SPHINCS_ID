! File:         submodule_standard_tpo_formulation_standard_tpo_variables.f90
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

SUBMODULE (standard_tpo_formulation) standard_tpo_variables

  !****************************************************
  !                                                   *
  !# Implementation of the methods of TYPE formul_tpo *
  !  that are called from the constructors and        *
  !  destructors of its EXTENDED TYPES                *
  !                                                   *
  !  FT 22.10.2020                                    *
  !                                                   *
  !****************************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE setup_standard_tpo_variables

    !*************************************************
    !                                                *
    !# Read the gravity grid parameters, computes    *
    !  gravity grid coordinates, imports the         *
    !  spacetime ID on the gravity grid, and         *
    !  performs some checks on it.                   *
    !  Its input includes the numbers of grid points *
    !  per axis, contrary to                         *
    !  construct_formul_tpo_bns_grid                 *
    !  where those numbers are replaced by the grid  *
    !  spacings.                                     *
    !                                                *
    !  FT 22.10.2020                                 *
    !                                                *
    !*************************************************

    USE mesh_refinement,  ONLY: levels, nlevels, initialize_grid, &
                                allocate_grid_function, &
                                deallocate_grid_function, &
                                coords, rad_coord
    !USE NaNChecker, ONLY: Check_Grid_Function_for_NAN
    USE tensor,           ONLY: jxx, jxy, jxz, &
                                jyy, jyz, jzz, n_sym3x3
    USE utility,          ONLY: determinant_sym3x3, one

    IMPLICIT NONE

    ! Index running over the refinement levels
    INTEGER:: l
    ! Indices running over the grids
    INTEGER:: i, j, k, i_matter

    ! Determinant of the standard 3+1 spatial metric
    DOUBLE PRECISION:: detg

    DOUBLE PRECISION, DIMENSION(6):: system_size
    DOUBLE PRECISION, DIMENSION(id% get_n_matter(),6):: sizes

    ! Get the number of matter objects in the physical system
    ftpo% n_matter= id% get_n_matter()

    !
    !-- Initialize timers
    !
    ftpo% grid_timer    = timer( "grid_timer" )
    ftpo% importer_timer= timer( "importer_timer" )

    CALL ftpo% grid_timer% start_timer()

    IF( PRESENT(dx) .AND. PRESENT(dy) .AND. PRESENT(dz) )THEN

      CALL initialize_grid( dx, dy, dz )

    ELSE

      CALL initialize_grid()

    ENDIF

    !PRINT *, ABS(id% get_center1_x()) + id% get_radius1_x_opp()
    !PRINT *, ABS(id% get_center2_x()) + id% get_radius2_x_opp()
    !PRINT *, ABS(levels(nlevels)% xR)

    !
    !-- Check that the stars are inside the finest refinement lvel
    !

    system_size= id% get_total_spatial_extent()

    IF( MAXVAL( ABS(system_size) ) > ABS(levels(nlevels)% xR) )THEN

      PRINT *
      PRINT *, "** The innermost, finest refinement level does not contain ", &
               "the entire system."
      PRINT *, "   Boundary of the innermost, finest level: ", &
               ABS(levels(nlevels)% xR), " Msun_geo"
      PRINT *, "   Size of the system: ", MAXVAL( ABS(system_size) ), &
               " Msun_geo"
      PRINT *, "   Please make the boundary of the innermost, finest level, ", &
               "larger than ", MAXVAL( ABS(system_size) ), &
               " Msun_geo"
      PRINT *, "   Stopping..."
      PRINT *
      STOP

    ENDIF

    CALL allocate_grid_function( ftpo% coords,    "coords_id", 3 )
    CALL allocate_grid_function( ftpo% rad_coord, 'rad_coord_id', 1 )

    ftpo% nlevels= nlevels
    ftpo% levels = levels

    ALLOCATE( ftpo% npoints_xaxis( ftpo% n_matter ) )

    DO i_matter= 1, ftpo% n_matter, 1

      sizes(i_matter,:)= id% return_spatial_extent(i_matter)

      ftpo% npoints_xaxis(i_matter)= FLOOR( ( sizes(i_matter,1) &
                                            + sizes(i_matter,2) ) &
                                              /ftpo% get_dx( ftpo% nlevels ) )

    ENDDO

    ref_levels: DO l= 1, ftpo% nlevels

      ftpo% coords%    levels(l)% var= coords%    levels(l)% var
      ftpo% rad_coord% levels(l)% var= rad_coord% levels(l)% var

    ENDDO ref_levels
    CALL deallocate_grid_function ( coords, 'coords' )
    CALL deallocate_grid_function ( rad_coord, 'rad_coord' )

    !
    !-- Allocating the memory for the grid functions
    !-- storing the spacetime ID at the grid points
    !
    CALL allocate_grid_function( ftpo% lapse,      "lapse_id",      1 )
    CALL allocate_grid_function( ftpo% shift_u,    "shift_u_id",    3 )
    CALL allocate_grid_function( ftpo% g_phys3_ll, "g_phys3_ll_id", 6 )
    CALL allocate_grid_function( ftpo% K_phys3_ll, "K_phys3_ll_id", 6 )

    CALL ftpo% grid_timer% stop_timer()

    !
    !-- Import the spacetime ID on the refined mesh,
    !-- and time the process
    !
    PRINT *
    PRINT *, "** Importing the spacetime ID on the refined mesh..."
    PRINT *
    CALL ftpo% importer_timer% start_timer()

    ref_levels2: DO l= 1, ftpo% nlevels, 1

      PRINT *, " * Importing on refinement level l=", l, "..."

      CALL id% read_id_spacetime( ftpo% get_ngrid_x(l), &
                                  ftpo% get_ngrid_y(l), &
                                  ftpo% get_ngrid_z(l), &
                                  ftpo% coords%     levels(l)% var, &
                                  ftpo% lapse%      levels(l)% var, &
                                  ftpo% shift_u%    levels(l)% var, &
                                  ftpo% g_phys3_ll% levels(l)% var, &
                                  ftpo% K_phys3_ll% levels(l)% var )

    ENDDO ref_levels2

    CALL ftpo% importer_timer% stop_timer()

    PRINT *, " * Spacetime ID imported on the gravity grid."

    !
    !-- Check that the imported ID does not contain NaNs
    !
    !CALL Check_Grid_Function_for_NAN( ftpo% lapse, "lapse" )
    !CALL Check_Grid_Function_for_NAN( ftpo% shift_u(:,:,:,jx), &
    !                                                    "shift_u_x" )
    !CALL Check_Grid_Function_for_NAN( ftpo% shift_u(:,:,:,jy), &
    !                                                    "shift_u_y" )
    !CALL Check_Grid_Function_for_NAN( ftpo% shift_u(:,:,:,jz), &
    !                                                    "shift_u_z" )
    !CALL Check_Grid_Function_for_NAN( ftpo% g_phys3_ll(:,:,:,jxx), &
    !                                                    "g_phys3_ll_jxx" )
    !CALL Check_Grid_Function_for_NAN( ftpo% g_phys3_ll(:,:,:,jxy), &
    !                                                    "g_phys3_ll_jxy" )
    !CALL Check_Grid_Function_for_NAN( ftpo% g_phys3_ll(:,:,:,jxz), &
    !                                                    "g_phys3_ll_jxz" )
    !CALL Check_Grid_Function_for_NAN( ftpo% g_phys3_ll(:,:,:,jyy), &
    !                                                    "g_phys3_ll_jyy" )
    !CALL Check_Grid_Function_for_NAN( ftpo% g_phys3_ll(:,:,:,jyz), &
    !                                                    "g_phys3_ll_jyz" )
    !CALL Check_Grid_Function_for_NAN( ftpo% g_phys3_ll(:,:,:,jzz), &
    !                                                    "g_phys3_ll_jzz" )
    !CALL Check_Grid_Function_for_NAN( ftpo% K_phys3_ll(:,:,:,jxx), &
    !                                                    "K_phys3_ll_jxx" )
    !CALL Check_Grid_Function_for_NAN( ftpo% K_phys3_ll(:,:,:,jxy), &
    !                                                    "K_phys3_ll_jxy" )
    !CALL Check_Grid_Function_for_NAN( ftpo% K_phys3_ll(:,:,:,jxz), &
    !                                                    "K_phys3_ll_jxz" )
    !CALL Check_Grid_Function_for_NAN( ftpo% K_phys3_ll(:,:,:,jyy), &
    !                                                    "K_phys3_ll_jyy" )
    !CALL Check_Grid_Function_for_NAN( ftpo% K_phys3_ll(:,:,:,jyz), &
    !                                                    "K_phys3_ll_jyz" )
    !CALL Check_Grid_Function_for_NAN( ftpo% K_phys3_ll(:,:,:,jzz), &
    !                                                    "K_phys3_ll_jzz" )

    !
    !-- Check that the determinant of the spatial metric is
    !-- strictly positive
    !
    DO l= 1, ftpo% nlevels, 1
      DO k= 1, ftpo% get_ngrid_z(l), 1
        DO j= 1, ftpo% get_ngrid_y(l), 1
          DO i= 1, ftpo% get_ngrid_x(l), 1

            CALL determinant_sym3x3( &
                              ftpo% g_phys3_ll% levels(l)% var(i,j,k,:), detg )

            IF( detg < 1D-10 )THEN

              PRINT *, "** ERROR! construct_formul_tpo_bns: The " &
                       // "determinant of the spatial metric is " &
                       // "effectively 0 at the grid point " &
                       // "(i,j,k)= (", i, ",", j,",",k, "), " &
                       // "(x,y,z)= ", "(", &
                       ftpo% coords% levels(l)% var( i, j, k, 1 ), ",", &
                       ftpo% coords% levels(l)% var( i, j, k, 2 ), ",", &
                       ftpo% coords% levels(l)% var( i, j, k, 3 ), ")."
              PRINT *
              PRINT *, ftpo% get_ngrid_x(l), ftpo% get_ngrid_y(l), &
                       ftpo% get_ngrid_z(l)
              PRINT *
              PRINT *, "detg=", detg
              PRINT *
              PRINT *, "g_xx=", ftpo% g_phys3_ll% levels(l)% var(i,j,k,jxx)
              PRINT *, "g_xy=", ftpo% g_phys3_ll% levels(l)% var(i,j,k,jxy)
              PRINT *, "g_xz=", ftpo% g_phys3_ll% levels(l)% var(i,j,k,jxz)
              PRINT *, "g_yy=", ftpo% g_phys3_ll% levels(l)% var(i,j,k,jyy)
              PRINT *, "g_yz=", ftpo% g_phys3_ll% levels(l)% var(i,j,k,jyz)
              PRINT *, "g_zz=", ftpo% g_phys3_ll% levels(l)% var(i,j,k,jzz)
              STOP

            ELSEIF( detg < 0 )THEN

              PRINT *, "** ERROR! construct_formul_tpo_bns: The " &
                       // "determinant of the spatial metric is " &
                       // "negative at the grid point " &
                       // "(i,j,k)= (", i, ",", j,",",k, "), " &
                       // "(x,y,z)= ", "(", &
                       ftpo% coords% levels(l)% var( i, j, k, 1 ), ",", &
                       ftpo% coords% levels(l)% var( i, j, k, 2 ), ",", &
                       ftpo% coords% levels(l)% var( i, j, k, 3 ), ")."
              PRINT *
              PRINT *, ftpo% get_ngrid_x(l), ftpo% get_ngrid_y(l), &
                       ftpo% get_ngrid_z(l)
              PRINT *
              PRINT *, "detg=", detg
              PRINT *
              PRINT *, "g_xx=", ftpo% g_phys3_ll% levels(l)% var(i,j,k,jxx)
              PRINT *, "g_xy=", ftpo% g_phys3_ll% levels(l)% var(i,j,k,jxy)
              PRINT *, "g_xz=", ftpo% g_phys3_ll% levels(l)% var(i,j,k,jxz)
              PRINT *, "g_yy=", ftpo% g_phys3_ll% levels(l)% var(i,j,k,jyy)
              PRINT *, "g_yz=", ftpo% g_phys3_ll% levels(l)% var(i,j,k,jyz)
              PRINT *, "g_zz=", ftpo% g_phys3_ll% levels(l)% var(i,j,k,jzz)
              STOP

            ENDIF

          ENDDO
        ENDDO
      ENDDO
    ENDDO

    IF( .NOT.ALLOCATED( ftpo% HC_int ))THEN
      ALLOCATE( ftpo% HC_int( ftpo% nlevels ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array HC_loo. ", &
                 "The error message is", err_msg
        STOP
      ENDIF
    ENDIF
    IF( .NOT.ALLOCATED( ftpo% MC_int ))THEN
      ALLOCATE( ftpo% MC_int( ftpo% nlevels, 3 ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array MC_loo. ", &
                 "The error message is", err_msg
        STOP
      ENDIF
    ENDIF
    ftpo% HC_int= HUGE(one)
    ftpo% MC_int= HUGE(one)

    IF( .NOT.ALLOCATED( ftpo% HC_parts_int ))THEN
      ALLOCATE( ftpo% HC_parts_int( ftpo% nlevels ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array MC_loo. ", &
                 "The error message is", err_msg
        STOP
      ENDIF
    ENDIF
    IF( .NOT.ALLOCATED( ftpo% MC_parts_int ))THEN
      ALLOCATE( ftpo% MC_parts_int( ftpo% nlevels, 3 ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array MC_loo. ", &
                 "The error message is", err_msg
        STOP
      ENDIF
    ENDIF
    ftpo% HC_parts_int= HUGE(one)
    ftpo% MC_parts_int= HUGE(one)

    PRINT *, " * Checked that the determinant of the spatial metric is", &
             " strictly positive."
    PRINT *

  END PROCEDURE setup_standard_tpo_variables


  MODULE PROCEDURE deallocate_standard_tpo_variables

    !***************************************************
    !                                                  *
    !# Core of the destructors of TYPES derived from   *
    !  formul_tpo. Their destructors should call this  *
    !  SUBROUTINE. It deallocates memory.              *
    !                                                  *
    !  FT                                              *
    !                                                  *
    !***************************************************

    USE mesh_refinement, ONLY: deallocate_grid_function

    IMPLICIT NONE

    IF( ALLOCATED( ftpo% coords% levels ) )THEN
      CALL deallocate_grid_function( ftpo% coords, "coords_id" )
    ENDIF

    IF( ALLOCATED( ftpo% rad_coord% levels ) )THEN
      CALL deallocate_grid_function( ftpo% rad_coord, "rad_coord_id" )
    ENDIF

    IF( ALLOCATED( ftpo% lapse% levels ) )THEN
      CALL deallocate_grid_function( ftpo% lapse, "lapse_id" )
    ENDIF

    IF( ALLOCATED( ftpo% shift_u% levels ) )THEN
      CALL deallocate_grid_function( ftpo% shift_u, "shift_u_id" )
    ENDIF

    IF( ALLOCATED( ftpo% g_phys3_ll% levels ) )THEN
      CALL deallocate_grid_function( ftpo% g_phys3_ll, "g_phys3_ll_id" )
    ENDIF

    IF( ALLOCATED( ftpo% K_phys3_ll% levels ) )THEN
      CALL deallocate_grid_function( ftpo% K_phys3_ll, "K_phys3_ll_id" )
    ENDIF

    IF( ALLOCATED( ftpo% HC% levels ) )THEN
      CALL deallocate_grid_function( ftpo% HC, "HC_id" )
    ENDIF

    IF( ALLOCATED( ftpo% HC_parts% levels ) )THEN
      CALL deallocate_grid_function( ftpo% HC_parts, "HC_parts_id" )
    ENDIF

    IF( ALLOCATED( ftpo% MC% levels ) )THEN
      CALL deallocate_grid_function( ftpo% MC, "MC_id" )
    ENDIF

    IF( ALLOCATED( ftpo% MC_parts% levels ) )THEN
      CALL deallocate_grid_function( ftpo% MC_parts, "MC_parts_id" )
    ENDIF

  END PROCEDURE deallocate_standard_tpo_variables


END SUBMODULE standard_tpo_variables
