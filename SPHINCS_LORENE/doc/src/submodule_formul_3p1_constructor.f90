! File:         submodule_formul_3p1_methods.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (formul_3p1_id) formul_3p1_constructor

  !***************************************************
  !                                                  *
  ! Implementation of the methods of TYPE formul_3p1 *
  ! that are called from the constructors and        *
  ! destructors of its EXTENDED TYPES                *
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

    IF( PRESENT(dx) .AND. PRESENT(dy) .AND. PRESENT(dz) )THEN

      CALL initialize_grid( dx, dy, dz )

    ELSE

      CALL initialize_grid()

    ENDIF

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


!  MODULE PROCEDURE construct_formul_3p1_bns_spacings
!
!    !************************************************
!    !                                               *
!    ! Read the gravity grid parameters, computes    *
!    ! gravity grid coordinates, imports the LORENE  *
!    ! spacetime ID on the gravity grid, and         *
!    ! performs some checks on it.                   *
!    ! Its input includes the grid spacings,         *
!    ! contrary to                                   *
!    ! construct_formul_3p1_bns                      *
!    ! where the grid spacings are replaced by the   *
!    ! numbers of grid points per axis.              *
!    !                                               *
!    ! FT 8.12.2020                                  *
!    !                                               *
!    ! Made compatible with mesh refinement.         *
!    !                                               *
!    ! FT 21.06.2021                                 *
!    !                                               *
!    !************************************************
!
!    USE grav_grid,  ONLY: ngrid_x, ngrid_y, ngrid_z, &
!                          xR, xL, yR, yL, zR, zL, &
!                          grid_file, gravity_grid, &
!                          check_parameter, &
!                          dx_1, dy_1, dz_1, dgmin, dgmin_1, tol
!    !USE NaNChecker, ONLY: Check_Grid_Function_for_NAN
!    USE tensor,     ONLY: itt, itx, ity, itz, ixx, ixy, &
!                          ixz, iyy, iyz, izz, jxx, jxy, jxz, &
!                          jyy, jyz, jzz, jx, jy, jz, n_sym3x3
!
!    IMPLICIT NONE
!
!    INTEGER:: ix, iy, iz, nx, ny, nz
!
!    DOUBLE PRECISION:: xtemp, ytemp, ztemp, detg
!
!    CHARACTER( LEN= : ), ALLOCATABLE :: field
!
!    LOGICAL:: file_exists
!
!    !
!    !-- Initialize timers
!    !
!    f3p1_obj% grid_timer    = timer( "grid_timer" )
!    f3p1_obj% importer_timer= timer( "importer_timer" )
!
!    CALL initialize_grid( dx, dy, dz )
!
!    CALL allocate_grid_function( f3p1_obj% coords,    "coords_id", 3 )
!    CALL allocate_grid_function( f3p1_obj% rad_coord, 'rad_coord_id', 1 )
!
!    f3p1_obj% nlevels= nlevels
!    f3p1_obj% levels = levels
!
!    ref_levels: DO l= 1, f3p1_obj% nlevels
!
!      f3p1_obj% coords%    levels(l)% var= coords%    levels(l)% var
!      f3p1_obj% rad_coord% levels(l)% var= rad_coord% levels(l)% var
!
!    ENDDO ref_levels
!    CALL deallocate_grid_function ( coords, 'coords' )
!    CALL deallocate_grid_function ( rad_coord, 'rad_coord' )
!
!
!
!    !
!    !-- Read the grid parameters
!    !
!!    grid_file= 'gravity_grid_parameters.dat'
!!
!!    INQUIRE( FILE= grid_file, EXIST= file_exists)
!!    IF( file_exists )THEN
!!       OPEN( 10, FILE= grid_file, STATUS= 'OLD' )
!!    ELSE
!!      grid_file= '../../BSSN/gravity_grid_parameters.dat'
!!      INQUIRE( FILE= grid_file, EXIST= file_exists )
!!      IF( file_exists )THEN
!!        OPEN( 10, FILE= grid_file, STATUS= 'OLD' )
!!      ELSE
!!       PRINT *
!!       PRINT *, '...grid parameter file not found...'
!!       PRINT *
!!       STOP
!!      ENDIF
!!    ENDIF
!!
!!    READ( UNIT= 10, NML= gravity_grid )
!!    CLOSE( UNIT= 10 )
!!
!!    !
!!    !-- Check the grid parameters
!!    !
!!    CALL check_parameter(ngrid_x,'ngrid_x')
!!    CALL check_parameter(ngrid_y,'ngrid_y')
!!    CALL check_parameter(ngrid_z,'ngrid_z')
!!    CALL check_parameter(xL,'xL')
!!    CALL check_parameter(xR,'xR')
!!    CALL check_parameter(yL,'yL')
!!    CALL check_parameter(yR,'yR')
!!    CALL check_parameter(zL,'zL')
!!    CALL check_parameter(zR,'zR')
!!
!!    ! Check whether they are the same
!!    IF( ABS(dy/dx - 1.0D0) > tol .OR. &
!!        ABS(dz/dx - 1.0D0) > tol )THEN
!!      PRINT*
!!      PRINT*,'...we currently assume that '
!!      PRINT*,'...dx=dy=dz '
!!      PRINT*,'...stopping...'
!!      STOP
!!    ENDIF
!!
!!    nx= CEILING( 1 + ( xR - xL )/dx )
!!    ny= CEILING( 1 + ( yR - yL )/dy )
!!    nz= CEILING( 1 + ( zR - zL )/dz )
!!
!!    f3p1_obj% ngrid_x= nx
!!    f3p1_obj% ngrid_y= ny
!!    f3p1_obj% ngrid_z= nz
!!
!!    f3p1_obj% dx= dx
!!    f3p1_obj% dy= dy
!!    f3p1_obj% dz= dz
!!
!!    ! Inverse (frequently needed in mapping)
!!    dx_1= 1.0D0/dx
!!    dy_1= 1.0D0/dy
!!    dz_1= 1.0D0/dz
!!
!!    f3p1_obj% dx_1= dx_1
!!    f3p1_obj% dy_1= dy_1
!!    f3p1_obj% dz_1= dz_1
!!
!!    ! Also store minimum grid spacing
!!    dgmin= MIN(dx,dy,dz)
!!
!!    ! ...and its inverse
!!    dgmin_1= 1.0D0/dgmin
!!
!!    !
!!    !-- Allocating the memory for the arrays x, y, z
!!    !-- storing the Cartesian coordinates of the grid points
!!    !
!!    IF( .NOT.ALLOCATED( f3p1_obj% grid ) )THEN
!!        ALLOCATE( f3p1_obj% grid( 3, nx, ny, nz ), &
!!                                            STAT= ios, ERRMSG= err_msg )
!!        IF( ios > 0 )THEN
!!          PRINT *, "...allocation error for array grid", &
!!                   ". The error message is", err_msg
!!          STOP
!!        ENDIF
!!        !CALL test_status( ios, err_msg, &
!!        !                "...allocation error for array grid in SUBROUTINE" &
!!        !                // "construct_formul_3p1." )
!!    ENDIF
!!
!!    !
!!    !-- Storing the coordinates of the grid points
!!    !
!!    PRINT *, "** Setting up the gravity grid..."
!!    CALL f3p1_obj% grid_timer% start_timer()
!!    coords_z: DO iz= 1, nz, 1
!!
!!      ztemp= zL + DBLE( iz - 1 )*dz
!!
!!      coords_y: DO iy= 1, ny, 1
!!
!!      ytemp= yL + DBLE( iy - 1 )*dy
!!
!!        coords_x: DO ix= 1, nx, 1
!!
!!          f3p1_obj% grid( 1, ix, iy, iz )= xL + DBLE( ix - 1 )*dx
!!          f3p1_obj% grid( 2, ix, iy, iz )= ytemp
!!          f3p1_obj% grid( 3, ix, iy, iz )= ztemp
!!
!!        ENDDO coords_x
!!      ENDDO coords_y
!!    ENDDO coords_z
!!    CALL f3p1_obj% grid_timer% stop_timer()
!!    PRINT *, " * Gravity grid set up."
!!    PRINT *
!!
!!    !
!!    !-- Allocating the memory for the arrays
!!    !-- storing the LORENE spacetime ID at the grid points
!!    !
!!    IF(.NOT.ALLOCATED( f3p1_obj% lapse ))THEN
!!      ALLOCATE( f3p1_obj% lapse( nx, ny, nz ), STAT= ios, ERRMSG= err_msg )
!!      IF( ios > 0 )THEN
!!        PRINT *, "...allocation error for array lapse", &
!!                 ". The error message is", err_msg
!!        STOP
!!      ENDIF
!!      !CALL test_status( ios, err_msg, &
!!      !            "...allocation error for array lapse in SUBROUTINE" &
!!      !            // "construct_formul_3p1." )
!!    ENDIF
!!    IF(.NOT.ALLOCATED( f3p1_obj% shift_u ))THEN
!!      ALLOCATE( f3p1_obj% shift_u( nx, ny, nz, 3 ), STAT= ios, ERRMSG= err_msg )
!!      IF( ios > 0 )THEN
!!        PRINT *, "...allocation error for array shift_u", &
!!                 ". The error message is", err_msg
!!        STOP
!!      ENDIF
!!      !CALL test_status( ios, err_msg, &
!!      !            "...allocation error for array shift_u in SUBROUTINE" &
!!      !            // "construct_formul_3p1." )
!!    ENDIF
!!    IF(.NOT.ALLOCATED( f3p1_obj% g_phys3_ll ))THEN
!!      ALLOCATE( f3p1_obj% g_phys3_ll( nx, ny, nz, 10 ), &
!!                                              STAT= ios, ERRMSG= err_msg )
!!      IF( ios > 0 )THEN
!!        PRINT *, "...allocation error for array g_phys3_ll", &
!!                 ". The error message is", err_msg
!!        STOP
!!      ENDIF
!!      !CALL test_status( ios, err_msg, &
!!      !        "...allocation error for array g_phys3_ll in SUBROUTINE" &
!!      !        // "construct_formul_3p1." )
!!    ENDIF
!!    IF(.NOT.ALLOCATED( f3p1_obj% K_phys3_ll ))THEN
!!      ALLOCATE( f3p1_obj% K_phys3_ll( nx, ny, nz, 10 ), &
!!                                              STAT= ios, ERRMSG= err_msg )
!!      IF( ios > 0 )THEN
!!        PRINT *, "...allocation error for array K_phys3_ll", &
!!                 ". The error message is", err_msg
!!        STOP
!!      ENDIF
!!      !CALL test_status( ios, err_msg, &
!!      !        "...allocation error for array K_phys3_ll in SUBROUTINE" &
!!      !        // "construct_formul_3p1." )
!!    ENDIF
!!
!!    !
!!    !-- Import the LORENE spacetime ID on the gravity grid
!!    !
!!    PRINT *, "** Importing the LORENE spacetime ID on the gravity grid..."
!!    PRINT *
!!    CALL f3p1_obj% importer_timer% start_timer()
!!    CALL bns_obj% import_id( nx, ny, nz, &
!!                             f3p1_obj% grid, &
!!                             f3p1_obj% lapse, &
!!                             f3p1_obj% shift_u, &
!!                             f3p1_obj% g_phys3_ll, &
!!                             f3p1_obj% K_phys3_ll )
!!    CALL f3p1_obj% importer_timer% stop_timer()
!!    PRINT *, " * LORENE spacetime ID imported on the gravity grid."
!!    PRINT *
!!
!!    !
!!    !-- Check that the imported ID does not contain NaNs
!!    !
!!    !CALL Check_Grid_Function_for_NAN( f3p1_obj% lapse, "lapse" )
!!    !
!!    !CALL Check_Grid_Function_for_NAN( f3p1_obj% shift_u(:,:,:,jx), &
!!    !                                                    "shift_u_x" )
!!    !CALL Check_Grid_Function_for_NAN( f3p1_obj% shift_u(:,:,:,jy), &
!!    !                                                    "shift_u_y" )
!!    !CALL Check_Grid_Function_for_NAN( f3p1_obj% shift_u(:,:,:,jz), &
!!    !                                                    "shift_u_z" )
!!    !CALL Check_Grid_Function_for_NAN( f3p1_obj% g_phys3_ll(:,:,:,jxx), &
!!    !                                                    "g_phys3_ll_jxx" )
!!    !CALL Check_Grid_Function_for_NAN( f3p1_obj% g_phys3_ll(:,:,:,jxy), &
!!    !                                                    "g_phys3_ll_jxy" )
!!    !CALL Check_Grid_Function_for_NAN( f3p1_obj% g_phys3_ll(:,:,:,jxz), &
!!    !                                                    "g_phys3_ll_jxz" )
!!    !CALL Check_Grid_Function_for_NAN( f3p1_obj% g_phys3_ll(:,:,:,jyy), &
!!    !                                                    "g_phys3_ll_jyy" )
!!    !CALL Check_Grid_Function_for_NAN( f3p1_obj% g_phys3_ll(:,:,:,jyz), &
!!    !                                                    "g_phys3_ll_jyz" )
!!    !CALL Check_Grid_Function_for_NAN( f3p1_obj% g_phys3_ll(:,:,:,jzz), &
!!    !                                                    "g_phys3_ll_jzz" )
!!    !CALL Check_Grid_Function_for_NAN( f3p1_obj% K_phys3_ll(:,:,:,jxx), &
!!    !                                                    "K_phys3_ll_jxx" )
!!    !CALL Check_Grid_Function_for_NAN( f3p1_obj% K_phys3_ll(:,:,:,jxy), &
!!    !                                                    "K_phys3_ll_jxy" )
!!    !CALL Check_Grid_Function_for_NAN( f3p1_obj% K_phys3_ll(:,:,:,jxz), &
!!    !                                                    "K_phys3_ll_jxz" )
!!    !CALL Check_Grid_Function_for_NAN( f3p1_obj% K_phys3_ll(:,:,:,jyy), &
!!    !                                                    "K_phys3_ll_jyy" )
!!    !CALL Check_Grid_Function_for_NAN( f3p1_obj% K_phys3_ll(:,:,:,jyz), &
!!    !                                                    "K_phys3_ll_jyz" )
!!    !CALL Check_Grid_Function_for_NAN( f3p1_obj% K_phys3_ll(:,:,:,jzz), &
!!    !                                                    "K_phys3_ll_jzz" )
!!
!!    !
!!    !-- Check that the determinant of the spatial metric is
!!    !-- strictly positive
!!    !
!!    DO iz= 1, f3p1_obj% ngrid_z, 1
!!      DO iy= 1, f3p1_obj% ngrid_y, 1
!!        DO ix= 1, f3p1_obj% ngrid_x, 1
!!
!!          detg= 2.0D0*f3p1_obj% g_phys3_ll(ix,iy,iz,jxy) &
!!                 *f3p1_obj% g_phys3_ll(ix,iy,iz,jxz) &
!!                 *f3p1_obj% g_phys3_ll(ix,iy,iz,jyz) &
!!                - f3p1_obj% g_phys3_ll(ix,iy,iz,jzz) &
!!                 *f3p1_obj% g_phys3_ll(ix,iy,iz,jxy)**2 &
!!                + f3p1_obj% g_phys3_ll(ix,iy,iz,jyy) &
!!                 *( f3p1_obj% g_phys3_ll(ix,iy,iz,jxx) &
!!                   *f3p1_obj% g_phys3_ll(ix,iy,iz,jzz) &
!!                  - f3p1_obj% g_phys3_ll(ix,iy,iz,jxz)**2 ) &
!!                - f3p1_obj% g_phys3_ll(ix,iy,iz,jxx) &
!!                 *f3p1_obj% g_phys3_ll(ix,iy,iz,jyz)**2
!!
!!          IF( detg < 1D-10 )THEN
!!            PRINT *, "** ERROR! construct_formul_3p1_bns: The " &
!!                     // "determinant of the spatial metric is " &
!!                     // "effectively 0 at the grid point " &
!!                     // "(ix,iy,iz)= (", ix, ",", iy,",",iz, "), " &
!!                     // "(x,y,z)= ", "(", &
!!                     f3p1_obj% grid( 1, ix, iy, iz ), ",", &
!!                     f3p1_obj% grid( 2, ix, iy, iz ), ",", &
!!                     f3p1_obj% grid( 3, ix, iy, iz ), ")."
!!            PRINT *
!!            PRINT *, f3p1_obj% ngrid_x, f3p1_obj% ngrid_y, &
!!                     f3p1_obj% ngrid_z
!!            PRINT *
!!            PRINT *, "detg=", detg
!!            PRINT *
!!            PRINT *, "g_xx=", f3p1_obj% g_phys3_ll(ix,iy,iz,jxx)
!!            PRINT *, "g_xy=", f3p1_obj% g_phys3_ll(ix,iy,iz,jxy)
!!            PRINT *, "g_xz=", f3p1_obj% g_phys3_ll(ix,iy,iz,jxz)
!!            PRINT *, "g_yy=", f3p1_obj% g_phys3_ll(ix,iy,iz,jyy)
!!            PRINT *, "g_yz=", f3p1_obj% g_phys3_ll(ix,iy,iz,jyz)
!!            PRINT *, "g_zz=", f3p1_obj% g_phys3_ll(ix,iy,iz,jzz)
!!            STOP
!!          ELSEIF( detg < 0 )THEN
!!            PRINT *, "** ERROR! construct_formul_3p1_bns: The " &
!!                     // "determinant of the spatial metric is " &
!!                     // "negative at the grid point " &
!!                     // "(ix,iy,iz)= (", ix, ",", iy,",",iz, "), " &
!!                     // "(x,y,z)= ", "(", &
!!                     f3p1_obj% grid( 1, ix, iy, iz ), ",", &
!!                     f3p1_obj% grid( 2, ix, iy, iz ), ",", &
!!                     f3p1_obj% grid( 3, ix, iy, iz ), ")."
!!            PRINT *, "detg=", detg
!!            PRINT *
!!            PRINT *, "g_xx=", f3p1_obj% g_phys3_ll(ix,iy,iz,jxx)
!!            PRINT *, "g_xy=", f3p1_obj% g_phys3_ll(ix,iy,iz,jxy)
!!            PRINT *, "g_xz=", f3p1_obj% g_phys3_ll(ix,iy,iz,jxz)
!!            PRINT *, "g_yy=", f3p1_obj% g_phys3_ll(ix,iy,iz,jyy)
!!            PRINT *, "g_yz=", f3p1_obj% g_phys3_ll(ix,iy,iz,jyz)
!!            PRINT *, "g_zz=", f3p1_obj% g_phys3_ll(ix,iy,iz,jzz)
!!            STOP
!!          ENDIF
!!
!!        ENDDO
!!      ENDDO
!!    ENDDO
!
!  END PROCEDURE construct_formul_3p1_bns_spacings


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


END SUBMODULE formul_3p1_constructor
