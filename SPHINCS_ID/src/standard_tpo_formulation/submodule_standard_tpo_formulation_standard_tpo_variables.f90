! File:         submodule_standard_tpo_formulation_standard_tpo_variables.f90
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

SUBMODULE (standard_tpo_formulation) standard_tpo_variables

  !****************************************************
  !
  !# Implementation of the methods of TYPE formul_tpo
  !  that are called from the constructors and
  !  destructors of its EXTENDED TYPES
  !
  !  FT 22.10.2020
  !
  !****************************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE setup_standard_tpo_variables

    !*************************************************
    !
    !# Read the gravity grid parameters, computes
    !  gravity grid coordinates, imports the
    !  spacetime ID on the gravity grid, and
    !  performs some checks on it.
    !  Its input includes the numbers of grid points
    !  per axis, contrary to
    !  construct_formul_tpo_bns_grid
    !  where those numbers are replaced by the grid
    !  spacings.
    !
    !  FT 22.10.2020
    !  Last updated: FT 05.07.2022
    !
    !*************************************************

    USE mesh_refinement,  ONLY: levels, nlevels, initialize_grid, &
                                allocate_grid_function, &
                                deallocate_grid_function, &
                                coords, rad_coord
    USE tensor,           ONLY: jx, jy, jz, jxx, jxy, jxz, &
                                jyy, jyz, jzz, n_sym3x3
    USE utility,          ONLY: zero, determinant_sym3x3, one, flag$tpo, &
                                scan_3d_array_for_nans


    IMPLICIT NONE


    INTEGER, SAVE:: tpo_id_number= -1
    !! Negative integer that identifies the [[tpo]] object

    ! Index running over the refinement levels
    INTEGER:: l
    ! Indices running over the grids
    INTEGER:: i, j, k, i_matter

    ! Determinant of the standard 3+1 spatial metric
    DOUBLE PRECISION:: detg

    DOUBLE PRECISION, DIMENSION(6):: system_size
    DOUBLE PRECISION, DIMENSION(id% get_n_matter(),6):: sizes


    tpof% tpo_id_number= tpo_id_number
    tpo_id_number      = tpo_id_number - 1

    ! Get the number of matter objects in the physical system
    tpof% n_matter= id% get_n_matter()

    !
    !-- Initialize timers
    !
    tpof% grid_timer    = timer( "grid_timer" )
    tpof% importer_timer= timer( "importer_timer" )

    CALL tpof% grid_timer% start_timer()

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

    CALL allocate_grid_function( tpof% coords,    "coords_id", 3 )
    CALL allocate_grid_function( tpof% rad_coord, 'rad_coord_id', 1 )

    tpof% nlevels= nlevels
    tpof% levels = levels

    ALLOCATE( tpof% npoints_xaxis( tpof% n_matter ) )

    DO i_matter= 1, tpof% n_matter, 1

      sizes(i_matter,:)= id% return_spatial_extent(i_matter)

      tpof% npoints_xaxis(i_matter)= FLOOR( ( sizes(i_matter,1) &
                                            + sizes(i_matter,2) ) &
                                              /tpof% get_dx( tpof% nlevels ) )

    ENDDO

    ref_levels: DO l= 1, tpof% nlevels

      tpof% coords%    levels(l)% var= coords%    levels(l)% var
      tpof% rad_coord% levels(l)% var= rad_coord% levels(l)% var

    ENDDO ref_levels
    CALL deallocate_grid_function ( coords, 'coords' )
    CALL deallocate_grid_function ( rad_coord, 'rad_coord' )

    !
    !-- Allocating the memory for the grid functions
    !-- storing the spacetime ID at the grid points
    !
    CALL allocate_grid_function( tpof% lapse,      "lapse_id",      1 )
    CALL allocate_grid_function( tpof% shift_u,    "shift_u_id",    3 )
    CALL allocate_grid_function( tpof% g_phys3_ll, "g_phys3_ll_id", 6 )
    CALL allocate_grid_function( tpof% K_phys3_ll, "K_phys3_ll_id", 6 )

    CALL tpof% grid_timer% stop_timer()

    CALL id% initialize_id(tpof% tpo_id_number, switch=.TRUE.)
    !CALL id% initialize_id(flag$tpo)

    !
    !-- Import the spacetime ID on the refined mesh,
    !-- and time the process
    !
    PRINT *
    PRINT *, "** Importing the spacetime ID on the refined mesh..."
    PRINT *
    CALL tpof% importer_timer% start_timer()

    ref_levels2: DO l= 1, tpof% nlevels, 1

      PRINT *, " * Importing on refinement level l=", l, "..."

      CALL id% initialize_id(l)

      CALL id% read_id_spacetime( tpof% get_ngrid_x(l), &
                                  tpof% get_ngrid_y(l), &
                                  tpof% get_ngrid_z(l), &
                                  tpof% coords%     levels(l)% var, &
                                  tpof% lapse%      levels(l)% var, &
                                  tpof% shift_u%    levels(l)% var, &
                                  tpof% g_phys3_ll% levels(l)% var, &
                                  tpof% K_phys3_ll% levels(l)% var )

    ENDDO ref_levels2

    CALL tpof% importer_timer% stop_timer()

    PRINT *, " * Spacetime ID imported on the gravity grid."

    !
    !-- Ensure that the standard 3+1 ID does not contain NaNs,
    !-- and that the determinant of the spatial metric is
    !-- strictly positive
    !
    PRINT *, "** Ensuring that the ID does not have any NaNs or infinities, ", &
             "and that the determinant of the spatial metric is strictly ", &
             "positive..."

    DO l= 1, tpof% nlevels, 1

      ASSOCIATE( nx     => tpof% get_ngrid_x(l), &
                 ny     => tpof% get_ngrid_y(l), &
                 nz     => tpof% get_ngrid_z(l), &
                 coords => tpof% coords%     levels(l)% var, &
                 lapse  => tpof% lapse%      levels(l)% var, &
                 shift  => tpof% shift_u%    levels(l)% var, &
                 g      => tpof% g_phys3_ll% levels(l)% var, &
                 eK     => tpof% K_phys3_ll% levels(l)% var )

      CALL scan_3d_array_for_nans( nx, ny, nz, lapse, "lapse" )

      CALL scan_3d_array_for_nans( nx, ny, nz, shift(:,:,:,jx), &
                                   "shift(:,:,:,jx)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, shift(:,:,:,jy), &
                                   "shift(:,:,:,jy)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, shift(:,:,:,jz), &
                                   "shift(:,:,:,jz)" )

      CALL scan_3d_array_for_nans( nx, ny, nz, g(:,:,:,jxx), &
                                   "g_phys3_ll(:,:,:,jxx)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, g(:,:,:,jxy), &
                                   "g_phys3_ll(:,:,:,jxy)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, g(:,:,:,jxz), &
                                   "g_phys3_ll(:,:,:,jxz)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, g(:,:,:,jyy), &
                                   "g_phys3_ll(:,:,:,jyy)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, g(:,:,:,jyz), &
                                   "g_phys3_ll(:,:,:,jyz)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, g(:,:,:,jzz), &
                                   "g_phys3_ll(:,:,:,jzz)" )

      CALL scan_3d_array_for_nans( nx, ny, nz, eK(:,:,:,jxx), &
                                   "K_phys3_ll(:,:,:,jxx)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, eK(:,:,:,jxy), &
                                   "K_phys3_ll(:,:,:,jxy)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, eK(:,:,:,jxz), &
                                   "K_phys3_ll(:,:,:,jxz)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, eK(:,:,:,jyy), &
                                   "K_phys3_ll(:,:,:,jyy)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, eK(:,:,:,jyz), &
                                   "K_phys3_ll(:,:,:,jyz)" )
      CALL scan_3d_array_for_nans( nx, ny, nz, eK(:,:,:,jzz), &
                                   "K_phys3_ll(:,:,:,jzz)" )

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( tpof, l ) &
      !$OMP             PRIVATE( i, j, k, detg )
      DO k= 1, nz, 1
        DO j= 1, ny, 1
          DO i= 1, nx, 1

            CALL determinant_sym3x3(g(i,j,k,:), detg)

            IF( detg < 1.D-10 )THEN

              PRINT *, "** ERROR in setup_standard_tpo_variables: The " &
                       // "determinant of the spatial metric is " &
                       // "effectively 0 at the grid point " &
                       // "(i,j,k)= (", i, ",", j,",",k, "), " &
                       // "(x,y,z)= ", "(", &
                       coords(i, j, k, 1), ",", &
                       coords(i, j, k, 2), ",", &
                       coords(i, j, k, 3), ")."
              PRINT *
              PRINT *, "   nx, ny, nz =", nx, ny, nz
              PRINT *
              PRINT *, "   detg=", detg
              PRINT *
              PRINT *, "   g_xx=", g(i,j,k,jxx)
              PRINT *, "   g_xy=", g(i,j,k,jxy)
              PRINT *, "   g_xz=", g(i,j,k,jxz)
              PRINT *, "   g_yy=", g(i,j,k,jyy)
              PRINT *, "   g_yz=", g(i,j,k,jyz)
              PRINT *, "   g_zz=", g(i,j,k,jzz)
              PRINT *
              STOP

            ELSEIF( detg < zero )THEN

              PRINT *, "** ERROR in setup_standard_tpo_variables: The " &
                       // "determinant of the spatial metric is " &
                       // "negative at the grid point " &
                       // "(i,j,k)= (", i, ",", j,",",k, "), " &
                       // "(x,y,z)= ", "(", &
                       coords(i, j, k, 1), ",", &
                       coords(i, j, k, 2), ",", &
                       coords(i, j, k, 3), ")."
              PRINT *
              PRINT *, "   nx, ny, nz =", nx, ny, nz
              PRINT *
              PRINT *, "   detg=", detg
              PRINT *
              PRINT *, "   g_xx=", g(i,j,k,jxx)
              PRINT *, "   g_xy=", g(i,j,k,jxy)
              PRINT *, "   g_xz=", g(i,j,k,jxz)
              PRINT *, "   g_yy=", g(i,j,k,jyy)
              PRINT *, "   g_yz=", g(i,j,k,jyz)
              PRINT *, "   g_zz=", g(i,j,k,jzz)
              PRINT *
              STOP

            ENDIF

          ENDDO
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

      END ASSOCIATE

    ENDDO

    PRINT *, "...the standard 3+1 ID does not contain NaNs or infinites, ", &
             "and the determinant of the spatial metric is strictly positive."
    PRINT *

    !
    !-- Initialize the arrays containing the integrals of the constraints over
    !-- the refined mesh
    !

    IF( .NOT.ALLOCATED( tpof% HC_int ))THEN
      ALLOCATE( tpof% HC_int( tpof% nlevels ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array HC_int. ", &
                 "The error message is", err_msg
        STOP
      ENDIF
    ENDIF
    IF( .NOT.ALLOCATED( tpof% MC_int ))THEN
      ALLOCATE( tpof% MC_int( tpof% nlevels, 3 ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array MC_int. ", &
                 "The error message is", err_msg
        STOP
      ENDIF
    ENDIF
    tpof% HC_int= HUGE(one)
    tpof% MC_int= HUGE(one)

    IF( .NOT.ALLOCATED( tpof% HC_parts_int ))THEN
      ALLOCATE( tpof% HC_parts_int( tpof% nlevels ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array MC_int. ", &
                 "The error message is", err_msg
        STOP
      ENDIF
    ENDIF
    IF( .NOT.ALLOCATED( tpof% MC_parts_int ))THEN
      ALLOCATE( tpof% MC_parts_int( tpof% nlevels, 3 ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array MC_int. ", &
                 "The error message is", err_msg
        STOP
      ENDIF
    ENDIF
    tpof% HC_parts_int= HUGE(one)
    tpof% MC_parts_int= HUGE(one)

  END PROCEDURE setup_standard_tpo_variables


  MODULE PROCEDURE deallocate_standard_tpo_variables

    !***************************************************
    !
    !# Core of the destructors of TYPES derived from
    !  formul_tpo. Their destructors should call this
    !  SUBROUTINE. It deallocates memory.
    !
    !  FT
    !
    !***************************************************

    USE mesh_refinement, ONLY: deallocate_grid_function

    IMPLICIT NONE

    IF( ALLOCATED( tpof% coords% levels ) )THEN
      CALL deallocate_grid_function( tpof% coords, "coords_id" )
    ENDIF

    IF( ALLOCATED( tpof% rad_coord% levels ) )THEN
      CALL deallocate_grid_function( tpof% rad_coord, "rad_coord_id" )
    ENDIF

    IF( ALLOCATED( tpof% lapse% levels ) )THEN
      CALL deallocate_grid_function( tpof% lapse, "lapse_id" )
    ENDIF

    IF( ALLOCATED( tpof% shift_u% levels ) )THEN
      CALL deallocate_grid_function( tpof% shift_u, "shift_u_id" )
    ENDIF

    IF( ALLOCATED( tpof% g_phys3_ll% levels ) )THEN
      CALL deallocate_grid_function( tpof% g_phys3_ll, "g_phys3_ll_id" )
    ENDIF

    IF( ALLOCATED( tpof% K_phys3_ll% levels ) )THEN
      CALL deallocate_grid_function( tpof% K_phys3_ll, "K_phys3_ll_id" )
    ENDIF

    IF( ALLOCATED( tpof% HC% levels ) )THEN
      CALL deallocate_grid_function( tpof% HC, "HC_id" )
    ENDIF

    IF( ALLOCATED( tpof% HC_parts% levels ) )THEN
      CALL deallocate_grid_function( tpof% HC_parts, "HC_parts_id" )
    ENDIF

    IF( ALLOCATED( tpof% MC% levels ) )THEN
      CALL deallocate_grid_function( tpof% MC, "MC_id" )
    ENDIF

    IF( ALLOCATED( tpof% MC_parts% levels ) )THEN
      CALL deallocate_grid_function( tpof% MC_parts, "MC_parts_id" )
    ENDIF

  END PROCEDURE deallocate_standard_tpo_variables


END SUBMODULE standard_tpo_variables
