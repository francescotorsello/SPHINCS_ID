! File:         submodule_formul_3p1_methods.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (formul_3p1_id) formul_3p1_methods

  !***************************************************
  !                                                  *
  ! Implementation of the methods of TYPE formul_3p1 *
  !                                                  *
  ! FT 22.10.2020                                    *
  !                                                  *
  !***************************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE construct_formul_3p1_bns

    !************************************************
    !                                               *
    ! Read the gravity grid parameters, computes    *
    ! gravity grid coordinates, imports the LORENE  *
    ! spacetime ID on the gravity grid, and         *
    ! performs some checks on it.                   *
    ! Its input includes the numbers of grid points *
    ! per axis, contrary to                         *
    ! construct_formul_3p1_bns_grid                 *
    ! where those numbers are replaced by the grid  *
    ! spacings.                                     *
    !                                               *
    ! FT 22.10.2020                                 *
    !                                               *
    !************************************************

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
    f3p1_obj% grid_timer    = timer( "grid_timer" )
    f3p1_obj% importer_timer= timer( "importer_timer" )

    !
    !-- Read the grid parameters
    !
    CALL initialize_grid()

    CALL allocate_grid_function( f3p1_obj% coords,    "coords_id", 3 )
    CALL allocate_grid_function( f3p1_obj% rad_coord, 'rad_coord_id', 1 )

    f3p1_obj% nlevels= nlevels
    f3p1_obj% levels = levels

    ref_levels: DO l= 1, f3p1_obj% nlevels

      f3p1_obj% coords%    levels(l)% var= coords%    levels(l)% var
      f3p1_obj% rad_coord% levels(l)% var= rad_coord% levels(l)% var

    ENDDO ref_levels
    CALL deallocate_grid_function ( coords, 'coords' )
    CALL deallocate_grid_function ( rad_coord, 'rad_coord' )

    !
    !-- Allocating the memory for the grid functions
    !-- storing the LORENE spacetime ID at the grid points
    !
    CALL allocate_grid_function( f3p1_obj% lapse,      "lapse_id",      1 )
    CALL allocate_grid_function( f3p1_obj% shift_u,    "shift_u_id",    3 )
    CALL allocate_grid_function( f3p1_obj% g_phys3_ll, "g_phys3_ll_id", 6 )
    CALL allocate_grid_function( f3p1_obj% K_phys3_ll, "K_phys3_ll_id", 6 )

    !
    !-- Import the LORENE spacetime ID on the refined mesh,
    !-- and time the process
    !
    PRINT *
    PRINT *, "** Importing the LORENE spacetime ID on the refined mesh..."
    PRINT *
    CALL f3p1_obj% importer_timer% start_timer()

    !PRINT *, "nlevels=", f3p1_obj% nlevels

    ref_levels2: DO l= 1, f3p1_obj% nlevels, 1

      PRINT *, " * Importing on refinement level l=", l, "..."

      !PRINT *, "get_ngrid_x=", f3p1_obj% get_ngrid_x(l)
      !PRINT *, "get_ngrid_y=", f3p1_obj% get_ngrid_y(l)
      !PRINT *, "get_ngrid_z=", f3p1_obj% get_ngrid_z(l)
      !PRINT *

      !PRINT *, f3p1_obj% coords% levels(l)% var

      CALL bns_obj% import_id( f3p1_obj% get_ngrid_x(l), &
                               f3p1_obj% get_ngrid_y(l), &
                               f3p1_obj% get_ngrid_z(l), &
                               f3p1_obj% coords%     levels(l)% var, &
                               f3p1_obj% lapse%      levels(l)% var, &
                               f3p1_obj% shift_u%    levels(l)% var, &
                               f3p1_obj% g_phys3_ll% levels(l)% var, &
                               f3p1_obj% K_phys3_ll% levels(l)% var )

    ENDDO ref_levels2

    CALL f3p1_obj% importer_timer% stop_timer()

    PRINT *, " * LORENE spacetime ID imported on the gravity grid."

    !
    !-- Check that the imported ID does not contain NaNs
    !
    !CALL Check_Grid_Function_for_NAN( f3p1_obj% lapse, "lapse" )
    !CALL Check_Grid_Function_for_NAN( f3p1_obj% shift_u(:,:,:,jx), &
    !                                                    "shift_u_x" )
    !CALL Check_Grid_Function_for_NAN( f3p1_obj% shift_u(:,:,:,jy), &
    !                                                    "shift_u_y" )
    !CALL Check_Grid_Function_for_NAN( f3p1_obj% shift_u(:,:,:,jz), &
    !                                                    "shift_u_z" )
    !CALL Check_Grid_Function_for_NAN( f3p1_obj% g_phys3_ll(:,:,:,jxx), &
    !                                                    "g_phys3_ll_jxx" )
    !CALL Check_Grid_Function_for_NAN( f3p1_obj% g_phys3_ll(:,:,:,jxy), &
    !                                                    "g_phys3_ll_jxy" )
    !CALL Check_Grid_Function_for_NAN( f3p1_obj% g_phys3_ll(:,:,:,jxz), &
    !                                                    "g_phys3_ll_jxz" )
    !CALL Check_Grid_Function_for_NAN( f3p1_obj% g_phys3_ll(:,:,:,jyy), &
    !                                                    "g_phys3_ll_jyy" )
    !CALL Check_Grid_Function_for_NAN( f3p1_obj% g_phys3_ll(:,:,:,jyz), &
    !                                                    "g_phys3_ll_jyz" )
    !CALL Check_Grid_Function_for_NAN( f3p1_obj% g_phys3_ll(:,:,:,jzz), &
    !                                                    "g_phys3_ll_jzz" )
    !CALL Check_Grid_Function_for_NAN( f3p1_obj% K_phys3_ll(:,:,:,jxx), &
    !                                                    "K_phys3_ll_jxx" )
    !CALL Check_Grid_Function_for_NAN( f3p1_obj% K_phys3_ll(:,:,:,jxy), &
    !                                                    "K_phys3_ll_jxy" )
    !CALL Check_Grid_Function_for_NAN( f3p1_obj% K_phys3_ll(:,:,:,jxz), &
    !                                                    "K_phys3_ll_jxz" )
    !CALL Check_Grid_Function_for_NAN( f3p1_obj% K_phys3_ll(:,:,:,jyy), &
    !                                                    "K_phys3_ll_jyy" )
    !CALL Check_Grid_Function_for_NAN( f3p1_obj% K_phys3_ll(:,:,:,jyz), &
    !                                                    "K_phys3_ll_jyz" )
    !CALL Check_Grid_Function_for_NAN( f3p1_obj% K_phys3_ll(:,:,:,jzz), &
    !                                                    "K_phys3_ll_jzz" )

    !
    !-- Check that the determinant of the spatial metric is
    !-- strictly positive
    !
    DO l= 1, f3p1_obj% nlevels, 1
      DO k= 1, f3p1_obj% get_ngrid_z(l), 1
        DO j= 1, f3p1_obj% get_ngrid_y(l), 1
          DO i= 1, f3p1_obj% get_ngrid_x(l), 1

            CALL determinant_sym3x3_grid( i, j, k, &
                                     f3p1_obj% g_phys3_ll% levels(l)% var, &
                                     detg )

            IF( detg < 1D-10 )THEN

              PRINT *, "** ERROR! construct_formul_3p1_bns: The " &
                       // "determinant of the spatial metric is " &
                       // "effectively 0 at the grid point " &
                       // "(ix,iy,iz)= (", i, ",", j,",",k, "), " &
                       // "(x,y,z)= ", "(", &
                       f3p1_obj% coords% levels(l)% var( i, j, k, 1 ), ",", &
                       f3p1_obj% coords% levels(l)% var( i, j, k, 2 ), ",", &
                       f3p1_obj% coords% levels(l)% var( i, j, k, 3 ), ")."
              PRINT *
              PRINT *, f3p1_obj% get_ngrid_x(l), f3p1_obj% get_ngrid_y(l), &
                       f3p1_obj% get_ngrid_z(l)
              PRINT *
              PRINT *, "detg=", detg
              PRINT *
              PRINT *, "g_xx=", f3p1_obj% g_phys3_ll% levels(l)% var(i,j,k,jxx)
              PRINT *, "g_xy=", f3p1_obj% g_phys3_ll% levels(l)% var(i,j,k,jxy)
              PRINT *, "g_xz=", f3p1_obj% g_phys3_ll% levels(l)% var(i,j,k,jxz)
              PRINT *, "g_yy=", f3p1_obj% g_phys3_ll% levels(l)% var(i,j,k,jyy)
              PRINT *, "g_yz=", f3p1_obj% g_phys3_ll% levels(l)% var(i,j,k,jyz)
              PRINT *, "g_zz=", f3p1_obj% g_phys3_ll% levels(l)% var(i,j,k,jzz)
              STOP

            ELSEIF( detg < 0 )THEN

              PRINT *, "** ERROR! construct_formul_3p1_bns: The " &
                       // "determinant of the spatial metric is " &
                       // "negative at the grid point " &
                       // "(ix,iy,iz)= (", i, ",", j,",",k, "), " &
                       // "(x,y,z)= ", "(", &
                       f3p1_obj% coords% levels(l)% var( i, j, k, 1 ), ",", &
                       f3p1_obj% coords% levels(l)% var( i, j, k, 2 ), ",", &
                       f3p1_obj% coords% levels(l)% var( i, j, k, 3 ), ")."
              PRINT *
              PRINT *, f3p1_obj% get_ngrid_x(l), f3p1_obj% get_ngrid_y(l), &
                       f3p1_obj% get_ngrid_z(l)
              PRINT *
              PRINT *, "detg=", detg
              PRINT *
              PRINT *, "g_xx=", f3p1_obj% g_phys3_ll% levels(l)% var(i,j,k,jxx)
              PRINT *, "g_xy=", f3p1_obj% g_phys3_ll% levels(l)% var(i,j,k,jxy)
              PRINT *, "g_xz=", f3p1_obj% g_phys3_ll% levels(l)% var(i,j,k,jxz)
              PRINT *, "g_yy=", f3p1_obj% g_phys3_ll% levels(l)% var(i,j,k,jyy)
              PRINT *, "g_yz=", f3p1_obj% g_phys3_ll% levels(l)% var(i,j,k,jyz)
              PRINT *, "g_zz=", f3p1_obj% g_phys3_ll% levels(l)% var(i,j,k,jzz)
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


 ! MODULE PROCEDURE construct_formul_3p1_bns_spacings
 !
 !   !************************************************
 !   !                                               *
 !   ! Read the gravity grid parameters, computes    *
 !   ! gravity grid coordinates, imports the LORENE  *
 !   ! spacetime ID on the gravity grid, and         *
 !   ! performs some checks on it.                   *
 !   ! Its input includes the grid spacings,         *
 !   ! contrary to                                   *
 !   ! construct_formul_3p1_bns                      *
 !   ! where the grid spacings are replaced by the   *
 !   ! numbers of grid points per axis.              *
 !   !                                               *
 !   ! FT 8.12.2020                                  *
 !   !                                               *
 !   !************************************************
 !
 !   USE grav_grid,  ONLY: ngrid_x, ngrid_y, ngrid_z, &
 !                         xR, xL, yR, yL, zR, zL, &
 !                         grid_file, gravity_grid, &
 !                         check_parameter, &
 !                         dx_1, dy_1, dz_1, dgmin, dgmin_1, tol
 !   !USE NaNChecker, ONLY: Check_Grid_Function_for_NAN
 !   USE tensor,     ONLY: itt, itx, ity, itz, ixx, ixy, &
 !                         ixz, iyy, iyz, izz, jxx, jxy, jxz, &
 !                         jyy, jyz, jzz, jx, jy, jz, n_sym3x3
 !
 !   IMPLICIT NONE
 !
 !   INTEGER:: ix, iy, iz, nx, ny, nz
 !
 !   DOUBLE PRECISION:: xtemp, ytemp, ztemp, detg
 !
 !   CHARACTER( LEN= : ), ALLOCATABLE :: field
 !
 !   LOGICAL:: file_exists
 !
 !   !
 !   !-- Initialize timers
 !   !
 !   f3p1_obj% grid_timer    = timer( "grid_timer" )
 !   f3p1_obj% importer_timer= timer( "importer_timer" )
 !
 !   !
 !   !-- Read the grid parameters
 !   !
 !   grid_file= 'gravity_grid_parameters.dat'
 !
 !   INQUIRE( FILE= grid_file, EXIST= file_exists)
 !   IF( file_exists )THEN
 !      OPEN( 10, FILE= grid_file, STATUS= 'OLD' )
 !   ELSE
 !     grid_file= '../../BSSN/gravity_grid_parameters.dat'
 !     INQUIRE( FILE= grid_file, EXIST= file_exists )
 !     IF( file_exists )THEN
 !       OPEN( 10, FILE= grid_file, STATUS= 'OLD' )
 !     ELSE
 !      PRINT *
 !      PRINT *, '...grid parameter file not found...'
 !      PRINT *
 !      STOP
 !     ENDIF
 !   ENDIF
 !
 !   READ( UNIT= 10, NML= gravity_grid )
 !   CLOSE( UNIT= 10 )
 !
 !   !
 !   !-- Check the grid parameters
 !   !
 !   CALL check_parameter(ngrid_x,'ngrid_x')
 !   CALL check_parameter(ngrid_y,'ngrid_y')
 !   CALL check_parameter(ngrid_z,'ngrid_z')
 !   CALL check_parameter(xL,'xL')
 !   CALL check_parameter(xR,'xR')
 !   CALL check_parameter(yL,'yL')
 !   CALL check_parameter(yR,'yR')
 !   CALL check_parameter(zL,'zL')
 !   CALL check_parameter(zR,'zR')
 !
 !   ! Check whether they are the same
 !   IF( ABS(dy/dx - 1.0D0) > tol .OR. &
 !       ABS(dz/dx - 1.0D0) > tol )THEN
 !     PRINT*
 !     PRINT*,'...we currently assume that '
 !     PRINT*,'...dx=dy=dz '
 !     PRINT*,'...stopping...'
 !     STOP
 !   ENDIF
 !
 !   nx= CEILING( 1 + ( xR - xL )/dx )
 !   ny= CEILING( 1 + ( yR - yL )/dy )
 !   nz= CEILING( 1 + ( zR - zL )/dz )
 !
 !   f3p1_obj% ngrid_x= nx
 !   f3p1_obj% ngrid_y= ny
 !   f3p1_obj% ngrid_z= nz
 !
 !   f3p1_obj% dx= dx
 !   f3p1_obj% dy= dy
 !   f3p1_obj% dz= dz
 !
 !   ! Inverse (frequently needed in mapping)
 !   dx_1= 1.0D0/dx
 !   dy_1= 1.0D0/dy
 !   dz_1= 1.0D0/dz
 !
 !   f3p1_obj% dx_1= dx_1
 !   f3p1_obj% dy_1= dy_1
 !   f3p1_obj% dz_1= dz_1
 !
 !   ! Also store minimum grid spacing
 !   dgmin= MIN(dx,dy,dz)
 !
 !   ! ...and its inverse
 !   dgmin_1= 1.0D0/dgmin
 !
 !   !
 !   !-- Allocating the memory for the arrays x, y, z
 !   !-- storing the Cartesian coordinates of the grid points
 !   !
 !   IF( .NOT.ALLOCATED( f3p1_obj% grid ) )THEN
 !       ALLOCATE( f3p1_obj% grid( 3, nx, ny, nz ), &
 !                                           STAT= ios, ERRMSG= err_msg )
 !       IF( ios > 0 )THEN
 !         PRINT *, "...allocation error for array grid", &
 !                  ". The error message is", err_msg
 !         STOP
 !       ENDIF
 !       !CALL test_status( ios, err_msg, &
 !       !                "...allocation error for array grid in SUBROUTINE" &
 !       !                // "construct_formul_3p1." )
 !   ENDIF
 !
 !   !
 !   !-- Storing the coordinates of the grid points
 !   !
 !   PRINT *, "** Setting up the gravity grid..."
 !   CALL f3p1_obj% grid_timer% start_timer()
 !   coords_z: DO iz= 1, nz, 1
 !
 !     ztemp= zL + DBLE( iz - 1 )*dz
 !
 !     coords_y: DO iy= 1, ny, 1
 !
 !     ytemp= yL + DBLE( iy - 1 )*dy
 !
 !       coords_x: DO ix= 1, nx, 1
 !
 !         f3p1_obj% grid( 1, ix, iy, iz )= xL + DBLE( ix - 1 )*dx
 !         f3p1_obj% grid( 2, ix, iy, iz )= ytemp
 !         f3p1_obj% grid( 3, ix, iy, iz )= ztemp
 !
 !       ENDDO coords_x
 !     ENDDO coords_y
 !   ENDDO coords_z
 !   CALL f3p1_obj% grid_timer% stop_timer()
 !   PRINT *, " * Gravity grid set up."
 !   PRINT *
 !
 !   !
 !   !-- Allocating the memory for the arrays
 !   !-- storing the LORENE spacetime ID at the grid points
 !   !
 !   IF(.NOT.ALLOCATED( f3p1_obj% lapse ))THEN
 !     ALLOCATE( f3p1_obj% lapse( nx, ny, nz ), STAT= ios, ERRMSG= err_msg )
 !     IF( ios > 0 )THEN
 !       PRINT *, "...allocation error for array lapse", &
 !                ". The error message is", err_msg
 !       STOP
 !     ENDIF
 !     !CALL test_status( ios, err_msg, &
 !     !            "...allocation error for array lapse in SUBROUTINE" &
 !     !            // "construct_formul_3p1." )
 !   ENDIF
 !   IF(.NOT.ALLOCATED( f3p1_obj% shift_u ))THEN
 !     ALLOCATE( f3p1_obj% shift_u( nx, ny, nz, 3 ), STAT= ios, ERRMSG= err_msg )
 !     IF( ios > 0 )THEN
 !       PRINT *, "...allocation error for array shift_u", &
 !                ". The error message is", err_msg
 !       STOP
 !     ENDIF
 !     !CALL test_status( ios, err_msg, &
 !     !            "...allocation error for array shift_u in SUBROUTINE" &
 !     !            // "construct_formul_3p1." )
 !   ENDIF
 !   IF(.NOT.ALLOCATED( f3p1_obj% g_phys3_ll ))THEN
 !     ALLOCATE( f3p1_obj% g_phys3_ll( nx, ny, nz, 10 ), &
 !                                             STAT= ios, ERRMSG= err_msg )
 !     IF( ios > 0 )THEN
 !       PRINT *, "...allocation error for array g_phys3_ll", &
 !                ". The error message is", err_msg
 !       STOP
 !     ENDIF
 !     !CALL test_status( ios, err_msg, &
 !     !        "...allocation error for array g_phys3_ll in SUBROUTINE" &
 !     !        // "construct_formul_3p1." )
 !   ENDIF
 !   IF(.NOT.ALLOCATED( f3p1_obj% K_phys3_ll ))THEN
 !     ALLOCATE( f3p1_obj% K_phys3_ll( nx, ny, nz, 10 ), &
 !                                             STAT= ios, ERRMSG= err_msg )
 !     IF( ios > 0 )THEN
 !       PRINT *, "...allocation error for array K_phys3_ll", &
 !                ". The error message is", err_msg
 !       STOP
 !     ENDIF
 !     !CALL test_status( ios, err_msg, &
 !     !        "...allocation error for array K_phys3_ll in SUBROUTINE" &
 !     !        // "construct_formul_3p1." )
 !   ENDIF
 !
 !   !
 !   !-- Import the LORENE spacetime ID on the gravity grid
 !   !
 !   PRINT *, "** Importing the LORENE spacetime ID on the gravity grid..."
 !   PRINT *
 !   CALL f3p1_obj% importer_timer% start_timer()
 !   CALL bns_obj% import_id( nx, ny, nz, &
 !                            f3p1_obj% grid, &
 !                            f3p1_obj% lapse, &
 !                            f3p1_obj% shift_u, &
 !                            f3p1_obj% g_phys3_ll, &
 !                            f3p1_obj% K_phys3_ll )
 !   CALL f3p1_obj% importer_timer% stop_timer()
 !   PRINT *, " * LORENE spacetime ID imported on the gravity grid."
 !   PRINT *
 !
 !   !
 !   !-- Check that the imported ID does not contain NaNs
 !   !
 !   !CALL Check_Grid_Function_for_NAN( f3p1_obj% lapse, "lapse" )
 !   !
 !   !CALL Check_Grid_Function_for_NAN( f3p1_obj% shift_u(:,:,:,jx), &
 !   !                                                    "shift_u_x" )
 !   !CALL Check_Grid_Function_for_NAN( f3p1_obj% shift_u(:,:,:,jy), &
 !   !                                                    "shift_u_y" )
 !   !CALL Check_Grid_Function_for_NAN( f3p1_obj% shift_u(:,:,:,jz), &
 !   !                                                    "shift_u_z" )
 !   !CALL Check_Grid_Function_for_NAN( f3p1_obj% g_phys3_ll(:,:,:,jxx), &
 !   !                                                    "g_phys3_ll_jxx" )
 !   !CALL Check_Grid_Function_for_NAN( f3p1_obj% g_phys3_ll(:,:,:,jxy), &
 !   !                                                    "g_phys3_ll_jxy" )
 !   !CALL Check_Grid_Function_for_NAN( f3p1_obj% g_phys3_ll(:,:,:,jxz), &
 !   !                                                    "g_phys3_ll_jxz" )
 !   !CALL Check_Grid_Function_for_NAN( f3p1_obj% g_phys3_ll(:,:,:,jyy), &
 !   !                                                    "g_phys3_ll_jyy" )
 !   !CALL Check_Grid_Function_for_NAN( f3p1_obj% g_phys3_ll(:,:,:,jyz), &
 !   !                                                    "g_phys3_ll_jyz" )
 !   !CALL Check_Grid_Function_for_NAN( f3p1_obj% g_phys3_ll(:,:,:,jzz), &
 !   !                                                    "g_phys3_ll_jzz" )
 !   !CALL Check_Grid_Function_for_NAN( f3p1_obj% K_phys3_ll(:,:,:,jxx), &
 !   !                                                    "K_phys3_ll_jxx" )
 !   !CALL Check_Grid_Function_for_NAN( f3p1_obj% K_phys3_ll(:,:,:,jxy), &
 !   !                                                    "K_phys3_ll_jxy" )
 !   !CALL Check_Grid_Function_for_NAN( f3p1_obj% K_phys3_ll(:,:,:,jxz), &
 !   !                                                    "K_phys3_ll_jxz" )
 !   !CALL Check_Grid_Function_for_NAN( f3p1_obj% K_phys3_ll(:,:,:,jyy), &
 !   !                                                    "K_phys3_ll_jyy" )
 !   !CALL Check_Grid_Function_for_NAN( f3p1_obj% K_phys3_ll(:,:,:,jyz), &
 !   !                                                    "K_phys3_ll_jyz" )
 !   !CALL Check_Grid_Function_for_NAN( f3p1_obj% K_phys3_ll(:,:,:,jzz), &
 !   !                                                    "K_phys3_ll_jzz" )
 !
 !   !
 !   !-- Check that the determinant of the spatial metric is
 !   !-- strictly positive
 !   !
 !   DO iz= 1, f3p1_obj% ngrid_z, 1
 !     DO iy= 1, f3p1_obj% ngrid_y, 1
 !       DO ix= 1, f3p1_obj% ngrid_x, 1
 !
 !         detg= 2.0D0*f3p1_obj% g_phys3_ll(ix,iy,iz,jxy) &
 !                *f3p1_obj% g_phys3_ll(ix,iy,iz,jxz) &
 !                *f3p1_obj% g_phys3_ll(ix,iy,iz,jyz) &
 !               - f3p1_obj% g_phys3_ll(ix,iy,iz,jzz) &
 !                *f3p1_obj% g_phys3_ll(ix,iy,iz,jxy)**2 &
 !               + f3p1_obj% g_phys3_ll(ix,iy,iz,jyy) &
 !                *( f3p1_obj% g_phys3_ll(ix,iy,iz,jxx) &
 !                  *f3p1_obj% g_phys3_ll(ix,iy,iz,jzz) &
 !                 - f3p1_obj% g_phys3_ll(ix,iy,iz,jxz)**2 ) &
 !               - f3p1_obj% g_phys3_ll(ix,iy,iz,jxx) &
 !                *f3p1_obj% g_phys3_ll(ix,iy,iz,jyz)**2
 !
 !         IF( detg < 1D-10 )THEN
 !           PRINT *, "** ERROR! construct_formul_3p1_bns: The " &
 !                    // "determinant of the spatial metric is " &
 !                    // "effectively 0 at the grid point " &
 !                    // "(ix,iy,iz)= (", ix, ",", iy,",",iz, "), " &
 !                    // "(x,y,z)= ", "(", &
 !                    f3p1_obj% grid( 1, ix, iy, iz ), ",", &
 !                    f3p1_obj% grid( 2, ix, iy, iz ), ",", &
 !                    f3p1_obj% grid( 3, ix, iy, iz ), ")."
 !           PRINT *
 !           PRINT *, f3p1_obj% ngrid_x, f3p1_obj% ngrid_y, &
 !                    f3p1_obj% ngrid_z
 !           PRINT *
 !           PRINT *, "detg=", detg
 !           PRINT *
 !           PRINT *, "g_xx=", f3p1_obj% g_phys3_ll(ix,iy,iz,jxx)
 !           PRINT *, "g_xy=", f3p1_obj% g_phys3_ll(ix,iy,iz,jxy)
 !           PRINT *, "g_xz=", f3p1_obj% g_phys3_ll(ix,iy,iz,jxz)
 !           PRINT *, "g_yy=", f3p1_obj% g_phys3_ll(ix,iy,iz,jyy)
 !           PRINT *, "g_yz=", f3p1_obj% g_phys3_ll(ix,iy,iz,jyz)
 !           PRINT *, "g_zz=", f3p1_obj% g_phys3_ll(ix,iy,iz,jzz)
 !           STOP
 !         ELSEIF( detg < 0 )THEN
 !           PRINT *, "** ERROR! construct_formul_3p1_bns: The " &
 !                    // "determinant of the spatial metric is " &
 !                    // "negative at the grid point " &
 !                    // "(ix,iy,iz)= (", ix, ",", iy,",",iz, "), " &
 !                    // "(x,y,z)= ", "(", &
 !                    f3p1_obj% grid( 1, ix, iy, iz ), ",", &
 !                    f3p1_obj% grid( 2, ix, iy, iz ), ",", &
 !                    f3p1_obj% grid( 3, ix, iy, iz ), ")."
 !           PRINT *, "detg=", detg
 !           PRINT *
 !           PRINT *, "g_xx=", f3p1_obj% g_phys3_ll(ix,iy,iz,jxx)
 !           PRINT *, "g_xy=", f3p1_obj% g_phys3_ll(ix,iy,iz,jxy)
 !           PRINT *, "g_xz=", f3p1_obj% g_phys3_ll(ix,iy,iz,jxz)
 !           PRINT *, "g_yy=", f3p1_obj% g_phys3_ll(ix,iy,iz,jyy)
 !           PRINT *, "g_yz=", f3p1_obj% g_phys3_ll(ix,iy,iz,jyz)
 !           PRINT *, "g_zz=", f3p1_obj% g_phys3_ll(ix,iy,iz,jzz)
 !           STOP
 !         ENDIF
 !
 !       ENDDO
 !     ENDDO
 !   ENDDO
 !
 ! END PROCEDURE construct_formul_3p1_bns_spacings


  MODULE PROCEDURE analyze_constraint

    IMPLICIT NONE

    INTEGER:: cnt_m7, cnt_m6, cnt_m5, cnt_m4, cnt_m3, cnt_m2, cnt_m1, cnt_0, &
              cnt_p1, cnt_p2, cnt_p3, cnt_oo, grid_points, i, j, k, &
              unit_analysis, nx, ny, nz

    DOUBLE PRECISION:: tmp, total

    LOGICAL:: exist

    IF( THIS% export_constraints_details )THEN
      !
      !-- Export the constraint analysis to a formatted file
      !
      unit_analysis= 20120

      INQUIRE( FILE= TRIM(name_analysis), EXIST= exist )

      IF( exist )THEN
          OPEN( UNIT= unit_analysis, FILE= TRIM(name_analysis), &
                STATUS= "REPLACE", &
                FORM= "FORMATTED", &
                POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
                IOMSG= err_msg )
      ELSE
          OPEN( UNIT= unit_analysis, FILE= TRIM(name_analysis), STATUS= "NEW", &
                FORM= "FORMATTED", &
                ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
      ENDIF
      IF( ios > 0 )THEN
        PRINT *, "... error when opening " // TRIM(name_analysis), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when opening " &
      !         // TRIM(name_analysis) )

      WRITE( UNIT = unit_analysis, IOSTAT = ios, IOMSG = err_msg, &
      FMT = * ) &
      "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id
      WRITE( UNIT= unit_analysis, IOSTAT = ios, &
             IOMSG = err_msg, FMT = * ) &
      "# The rows contain the points (x,y,z) at which ", name_constraint, &
      " has values in: (-oo,1D-7], [1D-7,1D-6], [1D-6,1D-5], [1D-5,1D-4]" &
      // ", [1D-4,1D-3], [1D-3,1D-2], [1D-2,1D-1], [1D-1,1], [1,1D+1]" &
      // ", [1D+1,1D+2], [1D+2,1D+3], [1D+3,+oo)"
    ENDIF

    CALL THIS% abs_values_in( 0.0D0, 1.0D-7, constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_m7 )

    CALL THIS% abs_values_in( 1.0D-7, 1.0D-6, constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_m6 )

    CALL THIS% abs_values_in( 1.0D-6, 1.0D-5, constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_m5 )

    CALL THIS% abs_values_in( 1.0D-5, 1.0D-4, constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_m4 )

    CALL THIS% abs_values_in( 1.0D-4, 1.0D-3, constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_m3 )

    CALL THIS% abs_values_in( 1.0D-3, 1.0D-2, constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_m2 )

    CALL THIS% abs_values_in( 1.0D-2, 1.0D-1, constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_m1 )

    CALL THIS% abs_values_in( 1.0D-1, 1.0D0, constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_0 )

    CALL THIS% abs_values_in( 1.0D0, 1D+1, constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_p1 )

    CALL THIS% abs_values_in( 1.0D+1, 1.0D+2, constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_p2 )

    CALL THIS% abs_values_in( 1.0D+2, 1.0D+3, constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_p3 )

    CALL THIS% abs_values_in( 1.0D+3, HUGE(DBLE(1)), constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_oo )

    CLOSE( UNIT= unit_analysis )

    nx= THIS% get_ngrid_x(l)
    ny= THIS% get_ngrid_y(l)
    nz= THIS% get_ngrid_z(l)
    grid_points= nx*ny*nz

    !
    !-- Compute the l2 norm of the constraints
    !
    l2_norm= 0
    DO k= 1, nz, 1
      DO j= 1, ny, 1
        DO i= 1, nx, 1
          l2_norm= l2_norm + constraint(i,j,k)*constraint(i,j,k)
        ENDDO
      ENDDO
    ENDDO
    l2_norm= SQRT( l2_norm/grid_points )

    !
    !-- Compute the loo norm (supremum norm) of the constraints
    !
    loo_norm= 0
    DO k= 1, nz, 1
      DO j= 1, ny, 1
        DO i= 1, nx, 1
          tmp= ABS( constraint(i,j,k) )
          IF( tmp > loo_norm )THEN
            loo_norm= tmp
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    !
    !-- Write a summary of the results to the logfile
    !
    WRITE( UNIT= unit_logfile, FMT = * ) "# The absolute values of ", &
         name_constraint, &
         " on the gravity grid are in the following intervals, on the ", &
         "given percentage of grid points:"
    WRITE( UNIT= unit_logfile, FMT = * ) ""
    WRITE( UNIT= unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (-oo,1D-7]: ", &
        100.0D0*DBLE(cnt_m7)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D-7,1D-6]: ", &
        100.0D0*DBLE(cnt_m6)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D-6,1D-5]: ", &
        100.0D0*DBLE(cnt_m5)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D-5,1D-4]: ", &
        100.0D0*DBLE(cnt_m4)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D-4,1D-3]: ", &
        100.0D0*DBLE(cnt_m3)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D-3,1D-2]: ", &
        100.0D0*DBLE(cnt_m2)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D-2,1D-1]: ", &
        100.0D0*DBLE(cnt_m1)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D-1,1]: ",   &
        100.0D0*DBLE(cnt_0)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1,10]: ",      &
        100.0D0*DBLE(cnt_p1)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (10,1D+2]: ",   &
        100.0D0*DBLE(cnt_p2)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D+2,1D+3]: ", &
        100.0D0*DBLE(cnt_p3)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D+3,+oo]: ",  &
        100.0D0*DBLE(cnt_oo)/DBLE(grid_points), "%"
    WRITE( UNIT = unit_logfile, FMT = * )
    WRITE( UNIT = unit_logfile, FMT = * ) "# l2-norm of ", name_constraint,&
                               " over the gravity grid= ", l2_norm
    WRITE( UNIT = unit_logfile, FMT = * )
    WRITE( UNIT = unit_logfile, FMT = * ) &
                       "# loo-norm (supremum of the absolute values) of ", &
                       name_constraint,&
                       " over the gravity grid= ", loo_norm
    WRITE( UNIT = unit_logfile, FMT = * )

    total= 100.0D0*DBLE(cnt_m7)/DBLE(grid_points) + &
           100.0D0*DBLE(cnt_m6)/DBLE(grid_points) + &
           100.0D0*DBLE(cnt_m5)/DBLE(grid_points) + &
           100.0D0*DBLE(cnt_m4)/DBLE(grid_points) + &
           100.0D0*DBLE(cnt_m3)/DBLE(grid_points) + &
           100.0D0*DBLE(cnt_m2)/DBLE(grid_points) + &
           100.0D0*DBLE(cnt_m1)/DBLE(grid_points) + &
           100.0D0*DBLE(cnt_0) /DBLE(grid_points) + &
           100.0D0*DBLE(cnt_p1)/DBLE(grid_points) + &
           100.0D0*DBLE(cnt_p2)/DBLE(grid_points) + &
           100.0D0*DBLE(cnt_p3)/DBLE(grid_points) + &
           100.0D0*DBLE(cnt_oo)/DBLE(grid_points)

    IF( total - DBLE(100) > 1.0D-4 )THEN
      PRINT *, " * WARNING! The percentages of the absolute values of ", &
               name_constraint, &
               " in the given intervals do not sum up to 100%. ", &
               " They sum up to", total, "%"
      PRINT *, "Check what happens in the SUBROUTINES analyze_constraint ", &
               "and abs_values_in."
    ENDIF

  END PROCEDURE analyze_constraint


  MODULE PROCEDURE abs_values_in

    !**************************************************
    !                                                 *
    ! Set "cnt" equal to the number of times that the *
    ! absolute value of "constraint" is in            *
    ! (lower_bound,upper_bound].                      *
    ! Depending on "export", it prints to file the    *
    ! grid points at which this happens.              *
    !                                                 *
    ! FT 24.03.2021                                   *
    !                                                 *
    !**************************************************

    USE tensor, ONLY: jx, jy, jz

    IMPLICIT NONE

    INTEGER:: i, j, k, nx, ny, nz

    nx= THIS% get_ngrid_x(l)
    ny= THIS% get_ngrid_y(l)
    nz= THIS% get_ngrid_z(l)

    cnt= 0
    DO k= 1, nz, 1
      DO j= 1, ny, 1
        DO i= 1, nx, 1
          IF( ABS( constraint(i,j,k) ) >= lower_bound .AND. &
              ABS( constraint(i,j,k) ) < upper_bound )THEN
            cnt= cnt + 1
            IF( export )THEN
              WRITE( UNIT= unit_analysis, IOSTAT = ios, &
                     IOMSG = err_msg, FMT = "(F17.13)", &
                     ADVANCE= "NO" ) &!"(E15.6)" &
                     THIS% coords% levels(l)% var( i, j, k, jx )
              WRITE( UNIT= unit_analysis, FMT= "(A2)", &
                     ADVANCE= "NO" ) "  "
              WRITE( UNIT= unit_analysis, IOSTAT = ios, &
                     IOMSG = err_msg, FMT = "(F17.13)", &
                     ADVANCE= "NO" ) &
                     THIS% coords% levels(l)% var( i, j, k, jy )
              WRITE( UNIT= unit_analysis, FMT= "(A2)", &
                     ADVANCE= "NO" ) "  "
              WRITE( UNIT= unit_analysis, IOSTAT = ios, &
                     IOMSG = err_msg, FMT = "(F17.13)", &
                     ADVANCE= "NO" ) &
                     THIS% coords% levels(l)% var( i, j, k, jz )
              WRITE( UNIT= unit_analysis, FMT= "(A2)", &
                     ADVANCE= "NO" ) "  "
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    IF( export .AND. cnt == 0 )THEN
      WRITE( UNIT= unit_analysis, FMT= "(I1)", ADVANCE= "NO" ) 0
      WRITE( UNIT= unit_analysis, FMT= "(I1)", ADVANCE= "NO" ) 0
      WRITE( UNIT= unit_analysis, FMT= "(I1)", ADVANCE= "NO" ) 0
    ENDIF
    WRITE( UNIT= unit_analysis, FMT= * ) ""

  END PROCEDURE abs_values_in


  MODULE PROCEDURE destruct_formul_3p1

    !**************************************************
    !                                                 *
    ! Core of the destructors of TYPES derived from   *
    ! formul_3p1. Their destructors should call this  *
    ! SUBROUTINE. It deallocates memory.              *
    !                                                 *
    ! FT                                              *
    !                                                 *
    !**************************************************

    USE mesh_refinement, ONLY: deallocate_grid_function

    IMPLICIT NONE

    IF( ALLOCATED( f3p1_obj% coords% levels ) )THEN
      CALL deallocate_grid_function( f3p1_obj% coords, "coords_id" )
    ENDIF

    IF( ALLOCATED( f3p1_obj% rad_coord% levels ) )THEN
      CALL deallocate_grid_function( f3p1_obj% rad_coord, "rad_coord_id" )
    ENDIF

    IF( ALLOCATED( f3p1_obj% lapse% levels ) )THEN
      CALL deallocate_grid_function( f3p1_obj% lapse, "lapse_id" )
    ENDIF

    IF( ALLOCATED( f3p1_obj% shift_u% levels ) )THEN
      CALL deallocate_grid_function( f3p1_obj% shift_u, "shift_u_id" )
    ENDIF

    IF( ALLOCATED( f3p1_obj% g_phys3_ll% levels ) )THEN
      CALL deallocate_grid_function( f3p1_obj% g_phys3_ll, "g_phys3_ll_id" )
    ENDIF

    IF( ALLOCATED( f3p1_obj% K_phys3_ll% levels ) )THEN
      CALL deallocate_grid_function( f3p1_obj% K_phys3_ll, "K_phys3_ll_id" )
    ENDIF

    IF( ALLOCATED( f3p1_obj% HC% levels ) )THEN
      CALL deallocate_grid_function( f3p1_obj% HC, "HC_id" )
    ENDIF

    IF( ALLOCATED( f3p1_obj% HC_parts% levels ) )THEN
      CALL deallocate_grid_function( f3p1_obj% HC_parts, "HC_parts_id" )
    ENDIF

    IF( ALLOCATED( f3p1_obj% MC% levels ) )THEN
      CALL deallocate_grid_function( f3p1_obj% MC, "MC_id" )
    ENDIF

    IF( ALLOCATED( f3p1_obj% MC_parts% levels ) )THEN
      CALL deallocate_grid_function( f3p1_obj% MC_parts, "MC_parts_id" )
    ENDIF

  END PROCEDURE destruct_formul_3p1


  !-----------------!
  !--  FUNCTIONS  --!
  !-----------------!


  MODULE PROCEDURE get_grid_point

    !*************************************************
    !                                                *
    ! Returns the array with the coordinates of the  *
    ! grid point                                     *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    USE tensor, ONLY: jx, jy, jz

    IMPLICIT NONE

    IF( i > THIS% levels(l)% ngrid_x )THEN
      PRINT *, "** ERROR in get_grid_point: i=", i, "> ngrid_x=", &
               THIS% levels(l)% ngrid_x, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( j > THIS% levels(l)% ngrid_y )THEN
      PRINT *, "** ERROR in get_grid_point j=", j, "> ngrid_y=", &
               THIS% levels(l)% ngrid_y, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( k > THIS% levels(l)% ngrid_z )THEN
      PRINT *, "** ERROR in get_grid_point k=", k, "> ngrid_z=", &
               THIS% levels(l)% ngrid_z, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF

    grid_point(1)= THIS% coords% levels(l)% var( i, j, k, jx )
    grid_point(2)= THIS% coords% levels(l)% var( i, j, k, jy )
    grid_point(3)= THIS% coords% levels(l)% var( i, j, k, jz )

  END PROCEDURE get_grid_point


  MODULE PROCEDURE get_nlevels

    !**************************************************
    !                                                 *
    ! Returns the number of refinement levels nlevels *
    !                                                 *
    ! FT 26.03.2021                                   *
    !                                                 *
    !**************************************************

    IMPLICIT NONE

    nlevels= THIS% nlevels

  END PROCEDURE get_nlevels


  MODULE PROCEDURE get_levels

    !**************************************************
    !                                                 *
    ! Returns the data structure levels               *
    !                                                 *
    ! FT 26.03.2021                                   *
    !                                                 *
    !**************************************************

    IMPLICIT NONE

    levels= THIS% levels

  END PROCEDURE get_levels


  MODULE PROCEDURE get_dx

    !*************************************************
    !                                                *
    ! Returns the grid spacing on the x axis         *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    dx= THIS% levels(l)% dx

  END PROCEDURE get_dx


  MODULE PROCEDURE get_dy

    !*************************************************
    !                                                *
    ! Returns the grid spacing on the y axis         *
    !                                                *
    ! FT 26.03.2021                                  *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    dy= THIS% levels(l)% dy

  END PROCEDURE get_dy


  MODULE PROCEDURE get_dz

    !*************************************************
    !                                                *
    ! Returns the grid spacing on the z axis         *
    !                                                *
    ! FT 26.03.2021                                  *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    dz= THIS% levels(l)% dz

  END PROCEDURE get_dz


  MODULE PROCEDURE get_ngrid_x

    !*************************************************
    !                                                *
    ! Returns the number of grid points on the x     *
    ! axis                                           *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    ngrid_x= THIS% levels(l)% ngrid_x

  END PROCEDURE get_ngrid_x


  MODULE PROCEDURE get_ngrid_y

    !*************************************************
    !                                                *
    ! Returns the number of grid points on the y     *
    ! axis                                           *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    ngrid_y= THIS% levels(l)% ngrid_y

  END PROCEDURE get_ngrid_y


  MODULE PROCEDURE get_ngrid_z

    !*************************************************
    !                                                *
    ! Returns the number of grid points on the z     *
    ! axis                                           *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    ngrid_z= THIS% levels(l)% ngrid_z

  END PROCEDURE get_ngrid_z


  MODULE PROCEDURE get_HC

    !**************************************************
    !                                                 *
    ! Returns the value of the Hamiltonian constraint *
    ! at the specified grid point                     *
    !                                                 *
    ! FT                                              *
    !                                                 *
    !**************************************************

    IMPLICIT NONE

    IF( i > THIS% levels(l)% ngrid_x )THEN
      PRINT *, "** ERROR in get_HC: i=", i, "> ngrid_x=", &
               THIS% levels(l)% ngrid_x, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( j > THIS% levels(l)% ngrid_y )THEN
      PRINT *, "** ERROR in get_HC: j=", j, "> ngrid_y=", &
               THIS% levels(l)% ngrid_y, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( k > THIS% levels(l)% ngrid_z )THEN
      PRINT *, "** ERROR in get_HC: k=", k, "> ngrid_z=", &
               THIS% levels(l)% ngrid_z, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF

    HC_value= THIS% HC% levels(l)% var( i, j, k )

  END PROCEDURE get_HC


  MODULE PROCEDURE get_MC

    !**************************************************
    !                                                 *
    ! Returns the array of values of the momentum     *
    ! constraint at the specified grid point          *
    !                                                 *
    ! FT                                              *
    !                                                 *
    !**************************************************

    USE tensor, ONLY: jx, jy, jz

    IMPLICIT NONE

    IF( i > THIS% levels(l)% ngrid_x )THEN
      PRINT *, "** ERROR in get_HC: i=", i, "> ngrid_x=", &
               THIS% levels(l)% ngrid_x, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( j > THIS% levels(l)% ngrid_y )THEN
      PRINT *, "** ERROR in get_HC: j=", j, "> ngrid_y=", &
               THIS% levels(l)% ngrid_y, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( k > THIS% levels(l)% ngrid_z )THEN
      PRINT *, "** ERROR in get_HC: k=", k, "> ngrid_z=", &
               THIS% levels(l)% ngrid_z, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF

    MC_value(1)= THIS% MC% levels(l)% var( i, j, k, jx )
    MC_value(2)= THIS% MC% levels(l)% var( i, j, k, jy )
    MC_value(3)= THIS% MC% levels(l)% var( i, j, k, jz )

  END PROCEDURE get_MC


  MODULE PROCEDURE get_HC_parts

    !**************************************************
    !                                                 *
    ! Returns the value of the Hamiltonian constraint *
    ! computed with particle data, at the specified   *
    ! grid point                                      *
    !                                                 *
    ! FT                                              *
    !                                                 *
    !**************************************************

    IMPLICIT NONE

    IF( i > THIS% levels(l)% ngrid_x )THEN
      PRINT *, "** ERROR in get_HC: i=", i, "> ngrid_x=", &
               THIS% levels(l)% ngrid_x, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( j > THIS% levels(l)% ngrid_y )THEN
      PRINT *, "** ERROR in get_HC: j=", j, "> ngrid_y=", &
               THIS% levels(l)% ngrid_y, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( k > THIS% levels(l)% ngrid_z )THEN
      PRINT *, "** ERROR in get_HC: k=", k, "> ngrid_z=", &
               THIS% levels(l)% ngrid_z, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF

    HC_value= THIS% HC_parts% levels(l)% var( i, j, k )

  END PROCEDURE get_HC_parts


  MODULE PROCEDURE get_MC_parts

    !**************************************************
    !                                                 *
    ! Returns the value of the momentum constraint    *
    ! computed with particle data, at the specified   *
    ! grid point                                      *
    !                                                 *
    ! FT                                              *
    !                                                 *
    !**************************************************

    USE tensor, ONLY: jx, jy, jz

    IMPLICIT NONE

    IF( i > THIS% levels(l)% ngrid_x )THEN
      PRINT *, "** ERROR in get_HC: i=", i, "> ngrid_x=", &
               THIS% levels(l)% ngrid_x, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( j > THIS% levels(l)% ngrid_y )THEN
      PRINT *, "** ERROR in get_HC: j=", j, "> ngrid_y=", &
               THIS% levels(l)% ngrid_y, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( k > THIS% levels(l)% ngrid_z )THEN
      PRINT *, "** ERROR in get_HC: k=", k, "> ngrid_z=", &
               THIS% levels(l)% ngrid_z, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF

    MC_value(1)= THIS% MC_parts% levels(l)% var( i, j, k, jx )
    MC_value(2)= THIS% MC_parts% levels(l)% var( i, j, k, jy )
    MC_value(3)= THIS% MC_parts% levels(l)% var( i, j, k, jz )

  END PROCEDURE get_MC_parts


END SUBMODULE formul_3p1_methods
