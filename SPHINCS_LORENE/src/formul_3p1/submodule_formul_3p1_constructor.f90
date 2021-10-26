! File:         submodule_formul_3p1_methods.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (formul_3p1_id) formul_3p1_constructor

  !****************************************************
  !                                                   *
  !# Implementation of the methods of TYPE formul_3p1 *
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


  MODULE PROCEDURE construct_formul_3p1_bns

    !*************************************************
    !                                                *
    !# Read the gravity grid parameters, computes    *
    !  gravity grid coordinates, imports the LORENE  *
    !  spacetime ID on the gravity grid, and         *
    !  performs some checks on it.                   *
    !  Its input includes the numbers of grid points *
    !  per axis, contrary to                         *
    !  construct_formul_3p1_bns_grid                 *
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
    USE tensor,           ONLY: itt, itx, ity, itz, ixx, ixy, &
                                ixz, iyy, iyz, izz, jxx, jxy, jxz, &
                                jyy, jyz, jzz, jx, jy, jz, n_sym3x3
    USE utility,          ONLY: determinant_sym3x3_grid

    IMPLICIT NONE

    ! Index running over the refinement levels
    INTEGER:: l
    ! Indices running over the grids
    INTEGER:: i, j, k

    ! Determinant of the standard 3+1 spatial metric
    DOUBLE PRECISION:: detg

    !
    !-- Initialize timers
    !
    f3p1% grid_timer    = timer( "grid_timer" )
    f3p1% importer_timer= timer( "importer_timer" )

    CALL f3p1% grid_timer% start_timer()

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
  !  IF( MAX( ABS(id% get_center1_x()) + id% get_radius1_x_opp(), &
  !           ABS(id% get_center2_x()) + id% get_radius2_x_opp() ) &
  !           > ABS(levels(nlevels)% xR) )THEN

    IF( MAXVAL( ABS(id% total_spatial_extent) ) > ABS(levels(nlevels)% xR) )THEN

      PRINT *
      PRINT *, "** The innermost, finest refinement level does not contain ", &
               "the entire system."
      PRINT *, "   Boundary of the innermost, finest level: ", &
               ABS(levels(nlevels)% xR), " Msun_geo"
      PRINT *, "   Size of the system: ", &
            MAXVAL( ABS(id% total_spatial_extent) ) > ABS(levels(nlevels)% xR), &
               " Msun_geo"
      PRINT *, "   Please make the boundary of the innermost, finest level, ", &
               "larger than ", &
            MAXVAL( ABS(id% total_spatial_extent) ) > ABS(levels(nlevels)% xR), &
               " Msun_geo"
      PRINT *, "   Stopping..."
      PRINT *
      STOP

    ENDIF

    CALL allocate_grid_function( f3p1% coords,    "coords_id", 3 )
    CALL allocate_grid_function( f3p1% rad_coord, 'rad_coord_id', 1 )

    f3p1% nlevels= nlevels
    f3p1% levels = levels

    ref_levels: DO l= 1, f3p1% nlevels

      f3p1% coords%    levels(l)% var= coords%    levels(l)% var
      f3p1% rad_coord% levels(l)% var= rad_coord% levels(l)% var

    ENDDO ref_levels
    CALL deallocate_grid_function ( coords, 'coords' )
    CALL deallocate_grid_function ( rad_coord, 'rad_coord' )

    !
    !-- Allocating the memory for the grid functions
    !-- storing the LORENE spacetime ID at the grid points
    !
    CALL allocate_grid_function( f3p1% lapse,      "lapse_id",      1 )
    CALL allocate_grid_function( f3p1% shift_u,    "shift_u_id",    3 )
    CALL allocate_grid_function( f3p1% g_phys3_ll, "g_phys3_ll_id", 6 )
    CALL allocate_grid_function( f3p1% K_phys3_ll, "K_phys3_ll_id", 6 )

    CALL f3p1% grid_timer% stop_timer()

    !
    !-- Import the LORENE spacetime ID on the refined mesh,
    !-- and time the process
    !
    PRINT *
    PRINT *, "** Importing the LORENE spacetime ID on the refined mesh..."
    PRINT *
    CALL f3p1% importer_timer% start_timer()

    !PRINT *, "nlevels=", f3p1% nlevels

    ref_levels2: DO l= 1, f3p1% nlevels, 1

      PRINT *, " * Importing on refinement level l=", l, "..."

      !PRINT *, "get_ngrid_x=", f3p1% get_ngrid_x(l)
      !PRINT *, "get_ngrid_y=", f3p1% get_ngrid_y(l)
      !PRINT *, "get_ngrid_z=", f3p1% get_ngrid_z(l)
      !PRINT *

      !PRINT *, f3p1% coords% levels(l)% var

      CALL id% read_id_spacetime( f3p1% get_ngrid_x(l), &
                             f3p1% get_ngrid_y(l), &
                             f3p1% get_ngrid_z(l), &
                             f3p1% coords%     levels(l)% var, &
                             f3p1% lapse%      levels(l)% var, &
                             f3p1% shift_u%    levels(l)% var, &
                             f3p1% g_phys3_ll% levels(l)% var, &
                             f3p1% K_phys3_ll% levels(l)% var )

    ENDDO ref_levels2

    CALL f3p1% importer_timer% stop_timer()

    PRINT *, " * LORENE spacetime ID imported on the gravity grid."

    !
    !-- Check that the imported ID does not contain NaNs
    !
    !CALL Check_Grid_Function_for_NAN( f3p1% lapse, "lapse" )
    !CALL Check_Grid_Function_for_NAN( f3p1% shift_u(:,:,:,jx), &
    !                                                    "shift_u_x" )
    !CALL Check_Grid_Function_for_NAN( f3p1% shift_u(:,:,:,jy), &
    !                                                    "shift_u_y" )
    !CALL Check_Grid_Function_for_NAN( f3p1% shift_u(:,:,:,jz), &
    !                                                    "shift_u_z" )
    !CALL Check_Grid_Function_for_NAN( f3p1% g_phys3_ll(:,:,:,jxx), &
    !                                                    "g_phys3_ll_jxx" )
    !CALL Check_Grid_Function_for_NAN( f3p1% g_phys3_ll(:,:,:,jxy), &
    !                                                    "g_phys3_ll_jxy" )
    !CALL Check_Grid_Function_for_NAN( f3p1% g_phys3_ll(:,:,:,jxz), &
    !                                                    "g_phys3_ll_jxz" )
    !CALL Check_Grid_Function_for_NAN( f3p1% g_phys3_ll(:,:,:,jyy), &
    !                                                    "g_phys3_ll_jyy" )
    !CALL Check_Grid_Function_for_NAN( f3p1% g_phys3_ll(:,:,:,jyz), &
    !                                                    "g_phys3_ll_jyz" )
    !CALL Check_Grid_Function_for_NAN( f3p1% g_phys3_ll(:,:,:,jzz), &
    !                                                    "g_phys3_ll_jzz" )
    !CALL Check_Grid_Function_for_NAN( f3p1% K_phys3_ll(:,:,:,jxx), &
    !                                                    "K_phys3_ll_jxx" )
    !CALL Check_Grid_Function_for_NAN( f3p1% K_phys3_ll(:,:,:,jxy), &
    !                                                    "K_phys3_ll_jxy" )
    !CALL Check_Grid_Function_for_NAN( f3p1% K_phys3_ll(:,:,:,jxz), &
    !                                                    "K_phys3_ll_jxz" )
    !CALL Check_Grid_Function_for_NAN( f3p1% K_phys3_ll(:,:,:,jyy), &
    !                                                    "K_phys3_ll_jyy" )
    !CALL Check_Grid_Function_for_NAN( f3p1% K_phys3_ll(:,:,:,jyz), &
    !                                                    "K_phys3_ll_jyz" )
    !CALL Check_Grid_Function_for_NAN( f3p1% K_phys3_ll(:,:,:,jzz), &
    !                                                    "K_phys3_ll_jzz" )

    !
    !-- Check that the determinant of the spatial metric is
    !-- strictly positive
    !
    DO l= 1, f3p1% nlevels, 1
      DO k= 1, f3p1% get_ngrid_z(l), 1
        DO j= 1, f3p1% get_ngrid_y(l), 1
          DO i= 1, f3p1% get_ngrid_x(l), 1

            CALL determinant_sym3x3_grid( i, j, k, &
                                     f3p1% g_phys3_ll% levels(l)% var, &
                                     detg )

            IF( detg < 1D-10 )THEN

              PRINT *, "** ERROR! construct_formul_3p1_bns: The " &
                       // "determinant of the spatial metric is " &
                       // "effectively 0 at the grid point " &
                       // "(ix,iy,iz)= (", i, ",", j,",",k, "), " &
                       // "(x,y,z)= ", "(", &
                       f3p1% coords% levels(l)% var( i, j, k, 1 ), ",", &
                       f3p1% coords% levels(l)% var( i, j, k, 2 ), ",", &
                       f3p1% coords% levels(l)% var( i, j, k, 3 ), ")."
              PRINT *
              PRINT *, f3p1% get_ngrid_x(l), f3p1% get_ngrid_y(l), &
                       f3p1% get_ngrid_z(l)
              PRINT *
              PRINT *, "detg=", detg
              PRINT *
              PRINT *, "g_xx=", f3p1% g_phys3_ll% levels(l)% var(i,j,k,jxx)
              PRINT *, "g_xy=", f3p1% g_phys3_ll% levels(l)% var(i,j,k,jxy)
              PRINT *, "g_xz=", f3p1% g_phys3_ll% levels(l)% var(i,j,k,jxz)
              PRINT *, "g_yy=", f3p1% g_phys3_ll% levels(l)% var(i,j,k,jyy)
              PRINT *, "g_yz=", f3p1% g_phys3_ll% levels(l)% var(i,j,k,jyz)
              PRINT *, "g_zz=", f3p1% g_phys3_ll% levels(l)% var(i,j,k,jzz)
              STOP

            ELSEIF( detg < 0 )THEN

              PRINT *, "** ERROR! construct_formul_3p1_bns: The " &
                       // "determinant of the spatial metric is " &
                       // "negative at the grid point " &
                       // "(ix,iy,iz)= (", i, ",", j,",",k, "), " &
                       // "(x,y,z)= ", "(", &
                       f3p1% coords% levels(l)% var( i, j, k, 1 ), ",", &
                       f3p1% coords% levels(l)% var( i, j, k, 2 ), ",", &
                       f3p1% coords% levels(l)% var( i, j, k, 3 ), ")."
              PRINT *
              PRINT *, f3p1% get_ngrid_x(l), f3p1% get_ngrid_y(l), &
                       f3p1% get_ngrid_z(l)
              PRINT *
              PRINT *, "detg=", detg
              PRINT *
              PRINT *, "g_xx=", f3p1% g_phys3_ll% levels(l)% var(i,j,k,jxx)
              PRINT *, "g_xy=", f3p1% g_phys3_ll% levels(l)% var(i,j,k,jxy)
              PRINT *, "g_xz=", f3p1% g_phys3_ll% levels(l)% var(i,j,k,jxz)
              PRINT *, "g_yy=", f3p1% g_phys3_ll% levels(l)% var(i,j,k,jyy)
              PRINT *, "g_yz=", f3p1% g_phys3_ll% levels(l)% var(i,j,k,jyz)
              PRINT *, "g_zz=", f3p1% g_phys3_ll% levels(l)% var(i,j,k,jzz)
              STOP

            ENDIF

          ENDDO
        ENDDO
      ENDDO
    ENDDO

    PRINT *, " * Checked that the determinant of the spatial metric is", &
             " strictly positive."
    PRINT *

  END PROCEDURE construct_formul_3p1_bns


  MODULE PROCEDURE destruct_formul_3p1

    !***************************************************
    !                                                  *
    !# Core of the destructors of TYPES derived from   *
    !  formul_3p1. Their destructors should call this  *
    !  SUBROUTINE. It deallocates memory.              *
    !                                                  *
    !  FT                                              *
    !                                                  *
    !***************************************************

    USE mesh_refinement, ONLY: deallocate_grid_function

    IMPLICIT NONE

    IF( ALLOCATED( f3p1% coords% levels ) )THEN
      CALL deallocate_grid_function( f3p1% coords, "coords_id" )
    ENDIF

    IF( ALLOCATED( f3p1% rad_coord% levels ) )THEN
      CALL deallocate_grid_function( f3p1% rad_coord, "rad_coord_id" )
    ENDIF

    IF( ALLOCATED( f3p1% lapse% levels ) )THEN
      CALL deallocate_grid_function( f3p1% lapse, "lapse_id" )
    ENDIF

    IF( ALLOCATED( f3p1% shift_u% levels ) )THEN
      CALL deallocate_grid_function( f3p1% shift_u, "shift_u_id" )
    ENDIF

    IF( ALLOCATED( f3p1% g_phys3_ll% levels ) )THEN
      CALL deallocate_grid_function( f3p1% g_phys3_ll, "g_phys3_ll_id" )
    ENDIF

    IF( ALLOCATED( f3p1% K_phys3_ll% levels ) )THEN
      CALL deallocate_grid_function( f3p1% K_phys3_ll, "K_phys3_ll_id" )
    ENDIF

    IF( ALLOCATED( f3p1% HC% levels ) )THEN
      CALL deallocate_grid_function( f3p1% HC, "HC_id" )
    ENDIF

    IF( ALLOCATED( f3p1% HC_parts% levels ) )THEN
      CALL deallocate_grid_function( f3p1% HC_parts, "HC_parts_id" )
    ENDIF

    IF( ALLOCATED( f3p1% MC% levels ) )THEN
      CALL deallocate_grid_function( f3p1% MC, "MC_id" )
    ENDIF

    IF( ALLOCATED( f3p1% MC_parts% levels ) )THEN
      CALL deallocate_grid_function( f3p1% MC_parts, "MC_parts_id" )
    ENDIF

  END PROCEDURE destruct_formul_3p1


END SUBMODULE formul_3p1_constructor
