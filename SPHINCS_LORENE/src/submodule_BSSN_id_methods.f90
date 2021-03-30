! File:         submodule_BSSN_id_methods.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (formul_bssn_id) bssn_id_methods

  !************************************************
  !                                               *
  ! Implementation of the methods of TYPE bssn_id *
  !                                               *
  ! FT 23.10.2020                                 *
  !                                              *
  ! Updated to mesh refinement                   *
  !                                              *
  ! FT 26.03.2021                                *
  !                                               *
  !************************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE compute_and_export_bssn_variables

    !************************************************
    !                                               *
    ! Compute, stores and export the BSSN variables *
    ! to a binary file to be read by the evolution  *
    ! code in SPHINCS                               *
    !                                               *
    ! FT 23.10.2020                                 *
    !                                               *
    !************************************************

    !USE NaNChecker,          ONLY: Check_Grid_Function_for_NAN
    USE mesh_refinement,            ONLY: nlevels, rad_coord, &
                                          allocate_grid_function, &
                                          deallocate_grid_function
    USE tensor,                     ONLY: itt, itx, ity, itz, ixx, ixy, &
                                          ixz, iyy, iyz, izz, jxx, jxy, jxz, &
                                          jyy, jyz, jzz, jx, jy, jz, n_sym3x3
    USE ADM_refine,                 ONLY: lapse, dt_lapse, shift_u, dt_shift_u, &
                                          K_phys3_ll, g_phys3_ll, &
                                          allocate_ADM, deallocate_ADM
    USE McLachlan_refine,           ONLY: allocate_Ztmp, deallocate_Ztmp, &
                                          ADM_to_BSSN
    USE Tmunu_refine,               ONLY: allocate_Tmunu, deallocate_Tmunu, &
                                          Tmunu_ll
    USE GravityAcceleration_refine, ONLY: dt_ehat_grav, dt_S_grav_l, &
                                          d_g_phys4_lll, &
                                          allocate_GravityAcceleration, &
                                          deallocate_GravityAcceleration
    USE options,                    ONLY: basename
    USE constants,                  ONLY: Msun_geo
    !
    !-- Use the arrays from the MODULE BSSN to store the BSSN variables
    !-- for the LORENE ID on the grid, and the SUBROUTINE write_BSSN_dump
    !-- to export them to the binary file needed by the evolution code
    !-- in SPHINCS
    !
    USE BSSN_refine, ONLY: allocate_BSSN, deallocate_BSSN, &
                           Gamma_u,          & ! Conformal connection
                           phi,              & ! Conformal factor
                           trK,              & ! Trace of extrinsic curvature
                           A_BSSN3_ll,       & ! Conformal traceless
                                               ! extrinsic curvature
                           g_BSSN3_ll,       & ! Conformal metric
                           Theta_Z4,         & ! Vector in the CCZ4 formulation
                                               ! Used because ADM_TO_BSSN
                                               ! calls SUBROUTINES that need it
                                               ! as input; however, it is not
                                               ! evolved in BSSN
                           lapse_A_BSSN,     & ! Time derivative of lapse
                           shift_B_BSSN_u,   & ! Time derivativeof shift
                           write_BSSN_dump

    IMPLICIT NONE

    ! The flag call_flag is set different than 0 if the SUBROUTINE
    ! compute_and_export_SPH_variables is called
    INTEGER, SAVE:: call_flag= 0

    INTEGER:: ix, iy, iz, i, j, k, allocation_status, l
    DOUBLE PRECISION:: detg


    PRINT *, "** Computing and exporting BSSN ID..."

    ! Allocate memory for the ADM MODULE variables (this has to be done since
    ! the MODULE SUBROUTINES need them; not allocating it results in a
    ! segmentation fault)
    PRINT *
    PRINT *, " * Allocating needed memory..."
    PRINT *

    !ngrid_x= THIS% ngrid_x
    !ngrid_y= THIS% ngrid_y
    !ngrid_z= THIS% ngrid_z
    !dx= THIS% dx
    !dy= THIS% dy
    !dz= THIS% dz
    !dx_1= THIS% dx_1
    !dy_1= THIS% dy_1
    !dz_1= THIS% dz_1

    CALL allocate_ADM()
    CALL allocate_BSSN()

    ! Allocate temporary memory for time integration
    CALL allocate_Ztmp()

    ! Allocate memory for the stress-energy tensor (used in
    ! write_BSSN_dump)
    CALL allocate_Tmunu()

    ! Allocate memory for the derivatives of the ADM variables
    CALL allocate_GravityAcceleration()

    CALL allocate_grid_function ( rad_coord, 'rad_coord' )

    ! Set the stress-energy tensor to 0 (here dummy)
    !Tmunu_ll(:,:,:,:)= 0.0D0

    !
    !-- Allocate memory for the coordinate radius
    !-- (this is needed by ADM_to_BSSN_args)
    !
    !IF( .NOT. ALLOCATED( rad_coord ) )THEN
    !    ALLOCATE( rad_coord( THIS% ngrid_x, THIS% ngrid_y, &
    !                         THIS% ngrid_z ), STAT= allocation_status )
    !ENDIF
    !IF( allocation_status > 0 )THEN
    !   PRINT *, '...allocation error for rad_coord'
    !   STOP
    !ENDIF
    !
    !DO k= 1, THIS% ngrid_z
    !  DO j= 1, THIS% ngrid_y
    !    DO i= 1, THIS% ngrid_x
    !      rad_coord( i, j, k )= SQRT( (xL+(i-1)*dx)**2 &
    !                                + (yL+(j-1)*dy)**2 &
    !                                + (zL+(k-1)*dz)**2 )
    !    ENDDO
    !  ENDDO
    !ENDDO

    ! Assign values to the MODULE variables, in order to call ADM_to_BSSN
    ref_levels: DO l= 1, THIS% nlevels

      Tmunu_ll% levels(l)% var= 0.0D0

      dt_lapse% levels(l)% var  = 0.0D0
      dt_shift_u% levels(l)% var= 0.0D0

      rad_coord% levels(l)% var = THIS% rad_coord% levels(l)% var
      lapse% levels(l)% var     = THIS% lapse% levels(l)% var
      shift_u% levels(l)% var   = THIS% shift_u% levels(l)% var
      g_phys3_ll% levels(l)% var= THIS% g_phys3_ll% levels(l)% var
      K_phys3_ll% levels(l)% var= THIS% K_phys3_ll% levels(l)% var

    ENDDO ref_levels

    !
    !-- Compute BSSN variables, and time the process
    !-- The BSSN variables are stored in the MODULE variables since
    !-- write_BSSN_dump need them
    !
    PRINT *, " * Computing BSSN variables..."
    PRINT *
    CALL THIS% bssn_computer_timer% start_timer()
    CALL ADM_to_BSSN()
    !CALL ADM_to_BSSN_args( &
    !  THIS% dx, THIS% dy, THIS% dz, &
    !  ! ADM variables (input)
    !  THIS% g_phys3_ll(:,:,:,jxx), THIS% g_phys3_ll(:,:,:,jxy), &
    !  THIS% g_phys3_ll(:,:,:,jxz), THIS% g_phys3_ll(:,:,:,jyy), &
    !  THIS% g_phys3_ll(:,:,:,jyz), THIS% g_phys3_ll(:,:,:,jzz), &
    !  THIS% K_phys3_ll(:,:,:,jxx), THIS% K_phys3_ll(:,:,:,jxy), &
    !  THIS% K_phys3_ll(:,:,:,jxz), THIS% K_phys3_ll(:,:,:,jyy), &
    !  THIS% K_phys3_ll(:,:,:,jyz), THIS% K_phys3_ll(:,:,:,jzz), &
    !  THIS% lapse(:,:,:), &
    !  THIS% shift_u(:,:,:,jx), &
    !  THIS% shift_u(:,:,:,jy), &
    !  THIS% shift_u(:,:,:,jz), &
    !  dt_lapse(:,:,:), &
    !  dt_shift_u(:,:,:,jx), dt_shift_u(:,:,:,jy), dt_shift_u(:,:,:,jz), &
    !  ! BSSN variables (output)
    !  g_BSSN3_ll(:,:,:,jxx), g_BSSN3_ll(:,:,:,jxy), &
    !  g_BSSN3_ll(:,:,:,jxz), g_BSSN3_ll(:,:,:,jyy), &
    !  g_BSSN3_ll(:,:,:,jyz), g_BSSN3_ll(:,:,:,jzz), &
    !  A_BSSN3_ll(:,:,:,jxx), A_BSSN3_ll(:,:,:,jxy), &
    !  A_BSSN3_ll(:,:,:,jxz), A_BSSN3_ll(:,:,:,jyy), &
    !  A_BSSN3_ll(:,:,:,jyz), A_BSSN3_ll(:,:,:,jzz), &
    !  phi(:,:,:), trK(:,:,:), Theta_Z4(:,:,:), &
    !  lapse_A_BSSN(:,:,:), &
    !  shift_B_BSSN_u(:,:,:,jx), shift_B_BSSN_u(:,:,:,jy), &
    !  shift_B_BSSN_u(:,:,:,jz), &
    !  Gamma_u(:,:,:,jx), Gamma_u(:,:,:,jy), &
    !  Gamma_u(:,:,:,jz) &
    !)
    CALL THIS% bssn_computer_timer% stop_timer()

    ! Set the MODULE variables equal to the TYPE variables
    !lapse= THIS% lapse
    !shift_u= THIS% shift_u

    !
    !-- Check the BSSN MODULE variables for NaNs
    !
    !CALL Check_Grid_Function_for_NAN( lapse, "lapse" )
    !CALL Check_Grid_Function_for_NAN( shift_u(:,:,:,jx), "shift_u_x" )
    !CALL Check_Grid_Function_for_NAN( shift_u(:,:,:,jy), "shift_u_y" )
    !CALL Check_Grid_Function_for_NAN( shift_u(:,:,:,jz), "shift_u_z" )
    !CALL Check_Grid_Function_for_NAN( g_BSSN3_ll(:,:,:,jxx), &
    !                                                    "g_BSSN3_ll_jxx" )
    !CALL Check_Grid_Function_for_NAN( g_BSSN3_ll(:,:,:,jxy), &
    !                                                    "g_BSSN3_ll_jxy" )
    !CALL Check_Grid_Function_for_NAN( g_BSSN3_ll(:,:,:,jxz), &
    !                                                    "g_BSSN3_ll_jxz" )
    !CALL Check_Grid_Function_for_NAN( g_BSSN3_ll(:,:,:,jyy), &
    !                                                    "g_BSSN3_ll_jyy" )
    !CALL Check_Grid_Function_for_NAN( g_BSSN3_ll(:,:,:,jyz), &
    !                                                    "g_BSSN3_ll_jyz" )
    !CALL Check_Grid_Function_for_NAN( g_BSSN3_ll(:,:,:,jzz), &
    !                                                    "g_BSSN3_ll_jzz" )
    !CALL Check_Grid_Function_for_NAN( A_BSSN3_ll(:,:,:,jxx), &
    !                                                    "A_BSSN3_ll_jxx" )
    !CALL Check_Grid_Function_for_NAN( A_BSSN3_ll(:,:,:,jxy), &
    !                                                    "A_BSSN3_ll_jxy" )
    !CALL Check_Grid_Function_for_NAN( A_BSSN3_ll(:,:,:,jxz), &
    !                                                    "A_BSSN3_ll_jxz" )
    !CALL Check_Grid_Function_for_NAN( A_BSSN3_ll(:,:,:,jyy), &
    !                                                    "A_BSSN3_ll_jyy" )
    !CALL Check_Grid_Function_for_NAN( A_BSSN3_ll(:,:,:,jyz), &
    !                                                    "A_BSSN3_ll_jyz" )
    !CALL Check_Grid_Function_for_NAN( A_BSSN3_ll(:,:,:,jzz), &
    !                                                    "A_BSSN3_ll_jzz" )
    !CALL Check_Grid_Function_for_NAN( phi, "phi" )
    !CALL Check_Grid_Function_for_NAN( trK, "trK" )
    !CALL Check_Grid_Function_for_NAN( Gamma_u(:,:,:,jx), "Gamma_u_x" )
    !CALL Check_Grid_Function_for_NAN( Gamma_u(:,:,:,jy), "Gamma_u_y" )
    !CALL Check_Grid_Function_for_NAN( Gamma_u(:,:,:,jz), "Gamma_u_z" )

    !
    !-- Setting the TYPE variables equal to the MODULE variables
    !
    ref_levels2: DO l= 1, THIS% nlevels

      THIS% Gamma_u% levels(l)% var   = Gamma_u% levels(l)% var
      THIS% phi% levels(l)% var       = phi% levels(l)% var
      THIS% trK% levels(l)% var       = trK% levels(l)% var
      THIS% g_BSSN3_ll% levels(l)% var= g_BSSN3_ll% levels(l)% var
      THIS% A_BSSN3_ll% levels(l)% var= A_BSSN3_ll% levels(l)% var

    ENDDO ref_levels2

    ! Write BSSN ID to a binary file to be read by the evolution code
    ! in SPHINCS
    IF( THIS% export_bin )THEN
      !IF( PRESENT(namefile) )THEN
        !CALL write_BSSN_dump( namefile )
        CALL write_BSSN_dump()
      !ELSE
      !  CALL write_BSSN_dump()
      !ENDIF
    ENDIF

    !
    !-- Deallocate MODULE variables
    !
    CALL deallocate_ADM()
    CALL deallocate_Ztmp()
    CALL deallocate_Tmunu()
    CALL deallocate_GravityAcceleration()
    CALL deallocate_BSSN()
    CALL deallocate_grid_function( rad_coord, 'rad_coord' )
    !CALL deallocate_gravity_grid()

    call_flag= call_flag + 1
    THIS% call_flag= call_flag

    PRINT *, "** BSSN ID computed."
    PRINT *

  END PROCEDURE compute_and_export_bssn_variables


!  MODULE PROCEDURE read_bssn_dump_print_formatted
!
!    !************************************************
!    !                                               *
!    ! Read the BSSN ID from the binary file output  *
!    ! by write_BSSN_dump, and print it to a         *
!    ! formatted file                                *
!    !                                               *
!    ! FT 08.02.2021                                 *
!    !                                               *
!    !************************************************
!
!    !USE grav_grid,           ONLY: ngrid_x, ngrid_y, ngrid_z, &
!    !                               !dx, dy, dz, dx_1, dy_1, dz_1, &
!    !                               !xR, xL, yR, yL, zR, zL, &
!    !                               !rad_coord, &
!    !                               deallocate_gravity_grid
!    USE mesh_refinement,     ONLY: levels
!    USE tensor,              ONLY: jxx, jxy, jxz, jyy, jyz, jzz, jx, jy, jz
!    USE ADM,                 ONLY: lapse, shift_u, &
!                                   allocate_ADM, deallocate_ADM
!    !USE McLachlan,           ONLY: initialize_BSSN, allocate_Ztmp, &
!    !                               deallocate_Ztmp, ADM_to_BSSN, &
!    !                               ADM_to_BSSN_args
!    !USE Tmunu,               ONLY: allocate_Tmunu, deallocate_Tmunu
!    !USE GravityAcceleration, ONLY: dt_ehat_grav, dt_S_grav_l, &
!    !                               d_g_phys4_lll, &
!    !                               allocate_GravityAcceleration, &
!    !                               deallocate_GravityAcceleration
!    USE BSSN,       ONLY: allocate_BSSN, deallocate_BSSN, &
!                          Gamma_u,          & ! Conformal connection
!                          phi,              & ! Conformal factor
!                          trK,              & ! Trace of extrinsic curvature
!                          A_BSSN3_ll,       & ! Conformal traceless
!                                              ! extrinsic curvature
!                          g_BSSN3_ll,       & ! Conformal metric
!                          Theta_Z4,         & ! Vector in the CCZ4 formulation.
!                                              ! Loaded here because ADM_TO_BSSN
!                                              ! calls SUBROUTINES that need it
!                                              ! as input; however, it is not
!                                              ! evolved in BSSN
!                          lapse_A_BSSN,     & ! Time derivative of lapse
!                          shift_B_BSSN_u,   & ! Time derivativeof shift
!                          read_BSSN_dump
!
!    IMPLICIT NONE
!
!    INTEGER:: ix, iy, iz, min_ix_y, min_iy_y, min_iz_y, &
!              min_ix_z, min_iy_z, min_iz_z
!
!    DOUBLE PRECISION:: min_abs_y, min_abs_z
!    DOUBLE PRECISION, DIMENSION( :, :, :, : ), ALLOCATABLE:: abs_grid
!
!    LOGICAL:: exist
!
!    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile
!
!    PRINT *, "** Executing the read_bssn_dump_print_formatted subroutine..."
!
!    !ngrid_x= THIS% ngrid_x
!    !ngrid_y= THIS% ngrid_y
!    !ngrid_z= THIS% ngrid_z
!    !dx= THIS% dx
!    !dy= THIS% dy
!    !dz= THIS% dz
!    !dx_1= THIS% dx_1
!    !dy_1= THIS% dy_1
!    !dz_1= THIS% dz_1
!
!    levels= THIS% levels
!
!    CALL allocate_ADM()
!    CALL allocate_BSSN()
!
!    ! Allocate temporary memory for time integration
!    !CALL allocate_Ztmp()
!
!    ! Allocate memory for the derivatives of the ADM variables
!    !CALL allocate_GravityAcceleration()
!
!    CALL read_BSSN_dump( ngrid_x, ngrid_y, ngrid_z, 00000, namefile_bin )
!
!    ! Being abs_grid a local array, it is good practice to allocate it on the
!    ! heap, otherwise it will be stored on the stack which has a very limited
!    ! size. This results in a segmentation fault.
!    ALLOCATE( abs_grid( 3, THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z ) )
!
!    IF( THIS% call_flag == 0 )THEN
!      PRINT *, "** The SUBROUTINE print_formatted_lorene_id_bssn_variables ", &
!        " must be called after compute_and_export_bssn_variables, otherwise", &
!        " there are no bssn fields to export to the formatted file."
!      PRINT *, "   Aborting."
!      PRINT *
!      STOP
!    ENDIF
!
!    IF( PRESENT(namefile) )THEN
!      finalnamefile= namefile
!    ELSE
!      finalnamefile= "bssn_vars.dat"
!    ENDIF
!
!    INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )
!
!    IF( exist )THEN
!      OPEN( UNIT= 20, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
!            FORM= "FORMATTED", &
!            POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
!            IOMSG= err_msg )
!    ELSE
!      OPEN( UNIT= 20, FILE= TRIM(finalnamefile), STATUS= "NEW", &
!      FORM= "FORMATTED", &
!            ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
!    ENDIF
!    IF( ios > 0 )THEN
!      PRINT *, "...error when opening ", TRIM(finalnamefile), &
!               ". The error message is", err_msg
!      STOP
!    ENDIF
!    !CALL test_status( ios, err_msg, "...error when opening " &
!    !         // TRIM(finalnamefile) )
!
!    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
!    "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id
!    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
!    "# Values of the fields (including coordinates) exported by LORENE "&
!    // "on each grid point"
!    IF( ios > 0 )THEN
!      PRINT *, "...error when writing line 1 in ", TRIM(finalnamefile), &
!               ". The error message is", err_msg
!      STOP
!    ENDIF
!    !CALL test_status( ios, err_msg, "...error when writing line 1 in "&
!    !         // TRIM(finalnamefile) )
!    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
!    "# column:      1        2       3       4       5", &
!    "       6       7       8", &
!    "       9       10      11", &
!    "       12      13      14", &
!    "       15      16      17      18      19", &
!    "       20      21      22", &
!    "       23      24"
!    IF( ios > 0 )THEN
!      PRINT *, "...error when writing line 2 in ", TRIM(finalnamefile), &
!               ". The error message is", err_msg
!      STOP
!    ENDIF
!    !CALL test_status( ios, err_msg, "...error when writing line 2 in "&
!    !        // TRIM(finalnamefile) )
!    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
!    "#      x [km]       y [km]       z [km]       lapse", &
!    "       shift_x [c]    shift_y [c]    shift_z [c]", &
!    "       conformal factor phi        trace of extr. curv. trK", &
!    "       g_BSSN_xx       g_BSSN_xy      g_BSSN_xz", &
!    "       g_BSSN_yy       g_BSSN_yz      g_BSSN_zz", &
!    "       A_BSSN_xx       A_BSSN_xy      A_BSSN_xz    ", &
!    "       A_BSSN_yy       A_BSSN_yz      A_BSSN_zz", &
!    "       Gamma_u_x       Gamma_u_y      Gamma_u_z"
!    IF( ios > 0 )THEN
!      PRINT *, "...error when writing line 3 in ", TRIM(finalnamefile), &
!               ". The error message is", err_msg
!      STOP
!    ENDIF
!    !CALL test_status( ios, err_msg, "...error when writing line 3 in "&
!    !        // TRIM(finalnamefile) )
!
!    DO iz= 1, THIS% ngrid_z, 1
!      DO iy= 1, THIS% ngrid_y, 1
!        DO ix= 1, THIS% ngrid_x, 1
!          abs_grid( 1, ix, iy, iz )= ABS( THIS% grid( 1, ix, iy, iz ) )
!          abs_grid( 2, ix, iy, iz )= ABS( THIS% grid( 2, ix, iy, iz ) )
!          abs_grid( 3, ix, iy, iz )= ABS( THIS% grid( 3, ix, iy, iz ) )
!        ENDDO
!      ENDDO
!    ENDDO
!
!    min_abs_y= 1D+20
!    min_abs_z= 1D+20
!    DO iz= 1, THIS% ngrid_z, 1
!      DO iy= 1, THIS% ngrid_y, 1
!        DO ix= 1, THIS% ngrid_x, 1
!          IF( ABS( THIS% grid( 2, ix, iy, iz ) ) < min_abs_y )THEN
!            min_abs_y= ABS( THIS% grid( 2, ix, iy, iz ) )
!            min_ix_y= ix
!            min_iy_y= iy
!            min_iz_y= iz
!          ENDIF
!          IF( ABS( THIS% grid( 3, ix, iy, iz ) ) < min_abs_z )THEN
!            min_abs_z= ABS( THIS% grid( 3, ix, iy, iz ) )
!            min_ix_z= ix
!            min_iy_z= iy
!            min_iz_z= iz
!          ENDIF
!        ENDDO
!      ENDDO
!    ENDDO
!
!    coords_z: DO iz= 1, THIS% ngrid_z, 1
!      coords_y: DO iy= 1, THIS% ngrid_y, 1
!        coords_x: DO ix= 1, THIS% ngrid_x, 1
!
!          IF( THIS% export_form_xy .AND. &
!              THIS% grid( 3, ix, iy, iz ) /= &
!              THIS% grid( 3, min_ix_z, min_iy_z, min_iz_z ) )THEN
!            CYCLE
!          ENDIF
!          IF( THIS% export_form_x .AND. &
!              ( THIS% grid( 3, ix, iy, iz ) /= &
!                THIS% grid( 3, min_ix_z, min_iy_z, min_iz_z ) &
!                .OR. &
!                THIS% grid( 2, ix, iy, iz ) /= &
!                THIS% grid( 2, min_ix_y, min_iy_y, min_iz_y ) ) )THEN
!            CYCLE
!          ENDIF
!
!          WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * )&
!              THIS% grid( 1, ix, iy, iz ), &
!              THIS% grid( 2, ix, iy, iz ), &
!              THIS% grid( 3, ix, iy, iz ), &
!              lapse( ix, iy, iz ), &
!              shift_u( ix, iy, iz, jx ), &
!              shift_u( ix, iy, iz, jy ), &
!              shift_u( ix, iy, iz, jz ), &
!              phi(ix,iy,iz), &
!              trK(ix,iy,iz), &
!              g_BSSN3_ll( ix, iy, iz, jxx ), &
!              g_BSSN3_ll( ix, iy, iz, jxy ), &
!              g_BSSN3_ll( ix, iy, iz, jxz ), &
!              g_BSSN3_ll( ix, iy, iz, jyy ), &
!              g_BSSN3_ll( ix, iy, iz, jyz ), &
!              g_BSSN3_ll( ix, iy, iz, jzz ), &
!              A_BSSN3_ll( ix, iy, iz, jxx ), &
!              A_BSSN3_ll( ix, iy, iz, jxy ), &
!              A_BSSN3_ll( ix, iy, iz, jxz ), &
!              A_BSSN3_ll( ix, iy, iz, jyy ), &
!              A_BSSN3_ll( ix, iy, iz, jyz ), &
!              A_BSSN3_ll( ix, iy, iz, jzz ), &
!              Gamma_u( ix, iy, iz, jx ), &
!              Gamma_u( ix, iy, iz, jy ), &
!              Gamma_u( ix, iy, iz, jz )
!
!          IF( ios > 0 )THEN
!            PRINT *, "...error when writing the arrays in ", &
!                     TRIM(finalnamefile), ". The error message is", err_msg
!            STOP
!          ENDIF
!          !CALL test_status( ios, err_msg, "...error when writing " &
!          !                  // "the arrays in " // TRIM(namefile) )
!
!        ENDDO coords_x
!      ENDDO coords_y
!    ENDDO coords_z
!
!    CLOSE( UNIT= 20 )
!
!    !
!    !-- Deallocate MODULE variables
!    !
!    CALL deallocate_ADM()
!    !CALL deallocate_Ztmp()
!    !CALL deallocate_GravityAcceleration()
!    CALL deallocate_BSSN()
!    !CALL deallocate_gravity_grid()
!
!    PRINT *, " * LORENE BSSN ID on the gravity grid saved to formatted " &
!             // "file ", TRIM(namefile)
!
!    PRINT *, "** Subroutine read_bssn_dump_print_formatted " &
!             // "executed."
!    PRINT *
!
!  END PROCEDURE read_bssn_dump_print_formatted


  MODULE PROCEDURE print_formatted_lorene_id_bssn_variables

    !************************************************
    !                                               *
    ! Print the BSSN ID, computed on the gravity    *
    ! grid from the LORENE ID, in a formatted file  *
    !                                               *
    ! FT 26.10.2020                                 *
    !                                               *
    !************************************************

    USE tensor,              ONLY: itt, itx, ity, itz, ixx, ixy, &
                                   ixz, iyy, iyz, izz, jxx, jxy, jxz, &
                                   jyy, jyz, jzz, jx, jy, jz
    USE mesh_refinement,     ONLY: output_1D, output_2D

    IMPLICIT NONE

    INTEGER:: i, j, k, min_ix_y, min_iy_y, min_iz_y, &
              min_ix_z, min_iy_z, min_iz_z

    DOUBLE PRECISION:: min_abs_y, min_abs_z
    DOUBLE PRECISION, DIMENSION( :, :, :, : ), ALLOCATABLE:: abs_grid

    LOGICAL:: exist

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    ! Being abs_grid a local array, it is good practice to allocate it on the
    ! heap, otherwise it will be stored on the stack which has a very limited
    ! size. This results in a segmentation fault.
    !ALLOCATE( abs_grid( 3, THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z ) )

    PRINT *, "** Executing the print_formatted_lorene_id_BSSN_variables " &
             // "subroutine..."

    IF( THIS% call_flag == 0 )THEN
      PRINT *, "** The SUBROUTINE print_formatted_lorene_id_bssn_variables ", &
        " must be called after compute_and_export_bssn_variables, otherwise", &
        " there are no bssn fields to export to the formatted file."
      PRINT *, "   Aborting."
      PRINT *
      STOP
    ENDIF

    IF( PRESENT(namefile) )THEN
      finalnamefile= namefile
    ELSE
      finalnamefile= "lorene-bns-id-bssn-form.dat"
    ENDIF

    INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

    IF( exist )THEN
      OPEN( UNIT= 20, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
            FORM= "FORMATTED", &
            POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
            IOMSG= err_msg )
    ELSE
      OPEN( UNIT= 20, FILE= TRIM(finalnamefile), STATUS= "NEW", &
      FORM= "FORMATTED", &
            ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening ", TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when opening " &
    !         // TRIM(finalnamefile) )

    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id
    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Values of the fields (including coordinates) exported by LORENE "&
    // "on each grid point"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 1 in ", TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 1 in "&
    !         // TRIM(finalnamefile) )
    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# column:      1        2       3       4       5", &
    "       6       7       8", &
    "       9       10      11", &
    "       12      13      14", &
    "       15      16      17      18      19", &
    "       20      21      22", &
    "       23      24"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 2 in ", TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 2 in "&
    !        // TRIM(finalnamefile) )
    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "#      x [km]       y [km]       z [km]       lapse", &
    "       shift_x [c]    shift_y [c]    shift_z [c]", &
    "       conformal factor phi        trace of extr. curv. trK", &
    "       g_BSSN_xx       g_BSSN_xy      g_BSSN_xz", &
    "       g_BSSN_yy       g_BSSN_yz      g_BSSN_zz", &
    "       A_BSSN_xx       A_BSSN_xy      A_BSSN_xz    ", &
    "       A_BSSN_yy       A_BSSN_yz      A_BSSN_zz", &
    "       Gamma_u_x       Gamma_u_y      Gamma_u_z"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 3 in ", TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 3 in "&
    !        // TRIM(finalnamefile) )

 !   DO iz= 1, THIS% ngrid_z, 1
 !     DO iy= 1, THIS% ngrid_y, 1
 !       DO ix= 1, THIS% ngrid_x, 1
 !         abs_grid( 1, ix, iy, iz )= ABS( THIS% grid( 1, ix, iy, iz ) )
 !         abs_grid( 2, ix, iy, iz )= ABS( THIS% grid( 2, ix, iy, iz ) )
 !         abs_grid( 3, ix, iy, iz )= ABS( THIS% grid( 3, ix, iy, iz ) )
 !       ENDDO
 !     ENDDO
 !   ENDDO
 !
 !   min_abs_y= 1D+20
 !   min_abs_z= 1D+20
 !   DO iz= 1, THIS% ngrid_z, 1
 !     DO iy= 1, THIS% ngrid_y, 1
 !       DO ix= 1, THIS% ngrid_x, 1
 !         IF( ABS( THIS% grid( 2, ix, iy, iz ) ) < min_abs_y )THEN
 !           min_abs_y= ABS( THIS% grid( 2, ix, iy, iz ) )
 !           min_ix_y= ix
 !           min_iy_y= iy
 !           min_iz_y= iz
 !         ENDIF
 !         IF( ABS( THIS% grid( 3, ix, iy, iz ) ) < min_abs_z )THEN
 !           min_abs_z= ABS( THIS% grid( 3, ix, iy, iz ) )
 !           min_ix_z= ix
 !           min_iy_z= iy
 !           min_iz_z= iz
 !         ENDIF
 !       ENDDO
 !     ENDDO
 !   ENDDO
 !
 !   coords_z: DO iz= 1, THIS% ngrid_z, 1
 !     coords_y: DO iy= 1, THIS% ngrid_y, 1
 !       coords_x: DO ix= 1, THIS% ngrid_x, 1
 !
 !         IF( THIS% export_form_xy .AND. &
 !             THIS% grid( 3, ix, iy, iz ) /= &
 !             THIS% grid( 3, min_ix_z, min_iy_z, min_iz_z ) )THEN
 !           CYCLE
 !         ENDIF
 !         IF( THIS% export_form_x .AND. &
 !             ( THIS% grid( 3, ix, iy, iz ) /= &
 !               THIS% grid( 3, min_ix_z, min_iy_z, min_iz_z ) &
 !               .OR. &
 !               THIS% grid( 2, ix, iy, iz ) /= &
 !               THIS% grid( 2, min_ix_y, min_iy_y, min_iz_y ) ) )THEN
 !           CYCLE
 !         ENDIF
 !
 !         WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * )&
 !             THIS% grid( 1, ix, iy, iz ), &
 !             THIS% grid( 2, ix, iy, iz ), &
 !             THIS% grid( 3, ix, iy, iz ), &
 !             THIS% lapse( ix, iy, iz ), &
 !             THIS% shift_u( ix, iy, iz, jx ), &
 !             THIS% shift_u( ix, iy, iz, jy ), &
 !             THIS% shift_u( ix, iy, iz, jz ), &
 !             THIS% phi(ix,iy,iz), &
 !             THIS% trK(ix,iy,iz), &
 !             THIS% g_BSSN3_ll( ix, iy, iz, jxx ), &
 !             THIS% g_BSSN3_ll( ix, iy, iz, jxy ), &
 !             THIS% g_BSSN3_ll( ix, iy, iz, jxz ), &
 !             THIS% g_BSSN3_ll( ix, iy, iz, jyy ), &
 !             THIS% g_BSSN3_ll( ix, iy, iz, jyz ), &
 !             THIS% g_BSSN3_ll( ix, iy, iz, jzz ), &
 !             THIS% A_BSSN3_ll( ix, iy, iz, jxx ), &
 !             THIS% A_BSSN3_ll( ix, iy, iz, jxy ), &
 !             THIS% A_BSSN3_ll( ix, iy, iz, jxz ), &
 !             THIS% A_BSSN3_ll( ix, iy, iz, jyy ), &
 !             THIS% A_BSSN3_ll( ix, iy, iz, jyz ), &
 !             THIS% A_BSSN3_ll( ix, iy, iz, jzz ), &
 !             THIS% Gamma_u( ix, iy, iz, jx ), &
 !             THIS% Gamma_u( ix, iy, iz, jy ), &
 !             THIS% Gamma_u( ix, iy, iz, jz )
 !
 !         IF( ios > 0 )THEN
 !           PRINT *, "...error when writing the arrays in ", &
 !                    TRIM(finalnamefile), ". The error message is", err_msg
 !           STOP
 !         ENDIF
 !         !CALL test_status( ios, err_msg, "...error when writing " &
 !         !                  // "the arrays in " // TRIM(finalnamefile) )
 !
 !       ENDDO coords_x
 !     ENDDO coords_y
 !   ENDDO coords_z

    IF( THIS% export_form_xy )THEN
      CALL output_2D( THIS% coords, 1, 3, 1, .TRUE. )
    ENDIF

    CLOSE( UNIT= 20 )

    PRINT *, " * LORENE BSSN ID on the gravity grid saved to formatted " &
             // "file ", TRIM(finalnamefile)

    PRINT *, "** Subroutine print_formatted_lorene_id_BSSN_variables " &
             // "executed."
    PRINT *

  END PROCEDURE print_formatted_lorene_id_bssn_variables


  MODULE PROCEDURE compute_and_export_bssn_constraints_grid

    !**************************************************
    !                                                 *
    ! Compute, store, analyze and export the BSSN     *
    ! constraints to a formatted file. The computaton *
    ! is done by importing the LORENE hydro ID on the *
    ! gravity grid, without any information on the    *
    ! particles.                                      *
    !                                                 *
    ! FT 1.02.2021                                    *
    !                                                 *
    !**************************************************

    USE constants,         ONLY: c_light2, cm2m, MSun, g2kg, m2cm, &
                                 lorene2hydrobase, MSun_geo, pi
    USE matrix,            ONLY: invert_4x4_matrix
    USE tensor,            ONLY: itt, itx, ity, itz, ixx, ixy, &
                                 ixz, iyy, iyz, izz, jxx, jxy, jxz, &
                                 jyy, jyz, jzz, jx, jy, jz, &
                                 it, ix, iy, iz, n_sym3x3, n_sym4x4
    USE mesh_refinement,   ONLY: allocate_grid_function, &
                                 levels, nlevels
    USE McLachlan_refine,  ONLY: BSSN_CONSTRAINTS_INTERIOR

    IMPLICIT NONE

    INTEGER:: i, j, k, allocation_status, fd_lim, l
    INTEGER, DIMENSION(3) :: imin, imax
    INTEGER:: unit_logfile, &
              min_ix_y, min_iy_y, min_iz_y, &
              min_ix_z, min_iy_z, min_iz_z

    DOUBLE PRECISION:: min_abs_y, min_abs_z
    DOUBLE PRECISION, DIMENSION( :, :, :, : ), ALLOCATABLE:: abs_grid

    TYPE(grid_function_scalar):: baryon_density
    TYPE(grid_function_scalar):: energy_density
    TYPE(grid_function_scalar):: specific_energy
    TYPE(grid_function_scalar):: pressure
    TYPE(grid_function):: v_euler
    TYPE(grid_function):: v_euler_l
    TYPE(grid_function):: u_euler_l
    TYPE(grid_function_scalar):: lorentz_factor
    DOUBLE PRECISION:: u_euler_norm= 0.0D0
    DOUBLE PRECISION:: detg4
    ! Spacetime metric
    TYPE(grid_function):: g4
    ! Stress-energy tensor
    TYPE(grid_function):: Tmunu_ll
    ! Spacetime metric as a 4x4 matrix
    DOUBLE PRECISION, DIMENSION( 4, 4 ):: g4temp
    ! Inverse spacetime metric as a 4x4 matrix
    DOUBLE PRECISION, DIMENSION( 4, 4 ):: ig4

    ! Declaration of debug variables needed to compute the Hamiltonian
    ! constraint directly, without calling the Cactus-bound SUBROUTINE
    ! BSSN_CONSTRAINTS_INTERIOR
    TYPE(grid_function_scalar):: HC_hand
    TYPE(grid_function_scalar):: HC_rho
    TYPE(grid_function_scalar):: HC_trK
    TYPE(grid_function_scalar):: HC_A
    TYPE(grid_function_scalar):: HC_derphi

    CHARACTER( LEN= : ), ALLOCATABLE:: name_constraint
    CHARACTER( LEN= : ), ALLOCATABLE:: name_analysis
    CHARACTER( LEN= : ), ALLOCATABLE:: finalname_logfile
    CHARACTER( LEN= 2 ):: n_reflev

    LOGICAL:: exist
    LOGICAL, PARAMETER:: debug= .FALSE.

    levels = THIS% levels
    nlevels= THIS% nlevels

    CALL allocate_grid_function( baryon_density, "baryon_density", 1 )
    CALL allocate_grid_function( energy_density, "energy_density", 1 )
    CALL allocate_grid_function( specific_energy, "specific_energy", 1 )
    CALL allocate_grid_function( pressure, "pressure", 1 )

    CALL allocate_grid_function( v_euler, "v_euler", 3 )
    CALL allocate_grid_function( v_euler_l, "v_euler_l", 3 )
    CALL allocate_grid_function( u_euler_l, "u_euler_l", 4 )
    CALL allocate_grid_function( lorentz_factor, "lorentz_factor", 1 )

    CALL allocate_grid_function( g4, "g4", n_sym4x4 )
    CALL allocate_grid_function( Tmunu_ll, "Tmunu_ll", n_sym4x4 )

    CALL allocate_grid_function( HC_hand, "HC_hand", 1 )
    CALL allocate_grid_function( HC_rho, "HC_rho", 1 )
    CALL allocate_grid_function( HC_trK, "HC_trK", 1 )
    CALL allocate_grid_function( HC_A, "HC_A", 1 )
    CALL allocate_grid_function( HC_derphi, "HC_derphi", 1 )

    CALL allocate_grid_function( THIS% HC, "HC_id", 1 )
    CALL allocate_grid_function( THIS% MC, "MC_id", 3 )
    CALL allocate_grid_function( THIS% GC, "GC_id", 3 )

    !
    !-- Import the hydro LORENE ID on the gravity grid
    !
    PRINT *, "** Importing LORENE hydro ID on the gravity grid..."
    PRINT *
    ref_levels: DO l= 1, THIS% nlevels, 1

      PRINT *, " * Importing on refinement level l=", l, "..."

      CALL bns_obj% import_id( THIS% get_ngrid_x(l), &
                               THIS% get_ngrid_y(l), &
                               THIS% get_ngrid_z(l), &
                               THIS% coords% levels(l)% var, &
                               baryon_density% levels(l)% var, &
                               energy_density% levels(l)% var, &
                               specific_energy% levels(l)% var, &
                               pressure% levels(l)% var, &
                               v_euler% levels(l)% var )

    ENDDO ref_levels
    PRINT *, " * LORENE hydro ID imported."
    PRINT *

    !---------------------------!
    !--  Compute constraints  --!
    !---------------------------!

    !
    !-- Compute the fluid 4-velocity in the coordinate frame
    !
    PRINT *, "** Computing fluid 4-velocity wrt Eulerian observer..."

    ref_levels2: DO l= 1, THIS% nlevels

      ASSOCIATE( v_euler_l      => v_euler_l% levels(l)% var, &
                 u_euler_l      => u_euler_l% levels(l)% var, &
                 v_euler        => v_euler% levels(l)% var, &
                 lorentz_factor => lorentz_factor% levels(l)% var, &
                 lapse          => THIS% lapse% levels(l)% var, &
                 shift_u        => THIS% shift_u% levels(l)% var, &
                 g_phys3_ll     => THIS% g_phys3_ll% levels(l)% var, &
                 g4             => g4% levels(l)% var, &
                 Tmunu_ll       => Tmunu_ll% levels(l)% var, &
                 energy_density => energy_density% levels(l)% var, &
                 pressure       => pressure% levels(l)% var &
      )

        DO k= 1, THIS% get_ngrid_z(l), 1
          DO j= 1, THIS% get_ngrid_y(l), 1
            DO i= 1, THIS% get_ngrid_x(l), 1

              !energy_density( i, j, k )= baryon_density( i, j, k ) &
              !                            + ( specific_energy(i,j,k) + 1.0 ) &
              !                                 *baryon_density( i, j, k )

              v_euler_l(i,j,k,jx)= g_phys3_ll(i,j,k,jxx)*v_euler(i,j,k,jx) &
                                 + g_phys3_ll(i,j,k,jxy)*v_euler(i,j,k,jy) &
                                 + g_phys3_ll(i,j,k,jxz)*v_euler(i,j,k,jz)
              v_euler_l(i,j,k,jy)= g_phys3_ll(i,j,k,jxy)*v_euler(i,j,k,jx) &
                                 + g_phys3_ll(i,j,k,jyy)*v_euler(i,j,k,jy) &
                                 + g_phys3_ll(i,j,k,jyz)*v_euler(i,j,k,jz)
              v_euler_l(i,j,k,jz)= g_phys3_ll(i,j,k,jxz)*v_euler(i,j,k,jx) &
                                 + g_phys3_ll(i,j,k,jyz)*v_euler(i,j,k,jy) &
                                 + g_phys3_ll(i,j,k,jzz)*v_euler(i,j,k,jz)

              lorentz_factor( i, j, k )= 1.0D0/SQRT( 1.0D0 &
                              - ( v_euler_l(i,j,k,jx)*v_euler(i,j,k,jx) &
                                + v_euler_l(i,j,k,jy)*v_euler(i,j,k,jy) &
                                + v_euler_l(i,j,k,jz)*v_euler(i,j,k,jz) ) )


              u_euler_l(i,j,k,it)= lorentz_factor( i, j, k ) &
                 *( - lapse( i, j, k ) &
                    + v_euler_l( i, j, k, jx )*shift_u( i, j, k, jx ) &
                    + v_euler_l( i, j, k, jy )*shift_u( i, j, k, jy ) &
                    + v_euler_l( i, j, k, jz )*shift_u( i, j, k, jz ) )
              u_euler_l(i,j,k,ix)= lorentz_factor( i, j, k ) &
                                     *v_euler_l( i, j, k, jx )
              u_euler_l(i,j,k,iy)= lorentz_factor( i, j, k ) &
                                     *v_euler_l( i, j, k, jy )
              u_euler_l(i,j,k,iz)= lorentz_factor( i, j, k ) &
                                     *v_euler_l( i, j, k, jz )

              CALL compute_g4( i, j, k, lapse, shift_u, g_phys3_ll, g4 )

              CALL determinant_sym4x4_grid( i, j, k, g4, detg4 )

              IF( ABS( detg4 ) < 1.0D-10 )THEN
                  PRINT *, "The determinant of the spacetime metric "&
                           // "is effectively 0 at the grid point " &
                           // "(i,j,k)= (", i, ",", j, ",", k, &
                              ")."
                  PRINT *, "detg4=", detg4
                  PRINT *
                  STOP
              ELSEIF( detg4 > 0.0D0 )THEN
                  PRINT *, "The determinant of the spacetime metric "&
                           // "is positive at the grid point " &
                           // "(i,j,k)= (", i, ",", j, ",", k, &
                              ")."
                  PRINT *, "detg4=", detg4
                  PRINT *
                  STOP
              ENDIF

              g4temp(1,1)= g4(i,j,k,itt)
              g4temp(1,2)= g4(i,j,k,itx)
              g4temp(1,3)= g4(i,j,k,ity)
              g4temp(1,4)= g4(i,j,k,itz)

              g4temp(2,1)= g4(i,j,k,itx)
              g4temp(2,2)= g4(i,j,k,ixx)
              g4temp(2,3)= g4(i,j,k,ixy)
              g4temp(2,4)= g4(i,j,k,ixz)

              g4temp(3,1)= g4(i,j,k,ity)
              g4temp(3,2)= g4(i,j,k,ixy)
              g4temp(3,3)= g4(i,j,k,iyy)
              g4temp(3,4)= g4(i,j,k,iyz)

              g4temp(4,1)= g4(i,j,k,itz)
              g4temp(4,2)= g4(i,j,k,ixz)
              g4temp(4,3)= g4(i,j,k,iyz)
              g4temp(4,4)= g4(i,j,k,izz)

              CALL invert_4x4_matrix( g4temp, ig4 )

              u_euler_norm= ig4(it,it)* &
                            u_euler_l(i,j,k,it)*u_euler_l(i,j,k,it) &
                          + 2.0D0*ig4(it,ix)* &
                            u_euler_l(i,j,k,it)*u_euler_l(i,j,k,ix) &
                          + 2.0D0*ig4(it,iy)* &
                            u_euler_l(i,j,k,it)*u_euler_l(i,j,k,iy) &
                          + 2.0D0*ig4(it,iz)* &
                            u_euler_l(i,j,k,it)*u_euler_l(i,j,k,iz) &
                          + ig4(ix,ix)* &
                            u_euler_l(i,j,k,ix)*u_euler_l(i,j,k,ix) &
                          + 2.0D0*ig4(ix,iy)* &
                            u_euler_l(i,j,k,ix)*u_euler_l(i,j,k,iy) &
                          + 2.0D0*ig4(ix,iz)* &
                            u_euler_l(i,j,k,ix)*u_euler_l(i,j,k,iz) &
                          + ig4(iy,iy)* &
                            u_euler_l(i,j,k,iy)*u_euler_l(i,j,k,iy) &
                          + 2.0D0*ig4(iy,iz)* &
                            u_euler_l(i,j,k,iy)*u_euler_l(i,j,k,iz) &
                          + 2.0D0*ig4(iz,iz)* &
                            u_euler_l(i,j,k,iz)*u_euler_l(i,j,k,iz)

              IF( ABS( u_euler_norm + 1.0D0 ) > 1.0D-4 )THEN
                  PRINT *, "** ERROR! The fluid 4-velocity in the " &
                           // "coordinate frame does not have norm -1. " &
                           // "The norm is", u_euler_norm
                  STOP
              ENDIF

            ! Print progress on screen
            perc= 100*(THIS% get_ngrid_x(l)*THIS% get_ngrid_y(l)*(k - 1) &
                  + THIS% get_ngrid_x(l)*(j - 1) + i) &
                  /( THIS% get_ngrid_x(l)*THIS% get_ngrid_y(l)* &
                     THIS% get_ngrid_z(l) )
            IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
              WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) creturn//" ", perc, "%"
            ENDIF

            ENDDO
          ENDDO
        ENDDO
        WRITE( *, "(A1)", ADVANCE= "NO" ) creturn
        PRINT *, " * Fluid 4-velocity wrt Eulerian observer computed."
        PRINT *

        ! Note that the units used in the spacetime part of SPHINCS are the
        ! same units as in the HydroBase thorn in the Einstein Toolkit.
        ! Such units can be found here, https://einsteintoolkit.org/thornguide/EinsteinBase/HydroBase/documentation.html
        ! The order of magnitude of the energy density can be found in
        ! https://www.ias.ac.in/article/fulltext/pram/084/05/0927-0941,
        ! and it is 150 MeV fm^{-3} ~ (2.4*10^{-11}J) / (10^{-45}m^3)
        !                           = 2.4*10^34 J m^{-3}

        !
        !-- Compute the stress-energy tensor
        !
        PRINT *, "** Computing stress-energy tensor..."
        Tmunu_ll= 0.0
        DO k= 1, THIS% get_ngrid_z(l), 1
          DO j= 1, THIS% get_ngrid_y(l), 1
            DO i= 1, THIS% get_ngrid_x(l), 1

              Tmunu_ll(i,j,k,itt)= lorene2hydrobase*( &
                      ( energy_density(i,j,k) + pressure(i,j,k) ) &
                      *u_euler_l(i,j,k,it)*u_euler_l(i,j,k,it) &
                      + pressure(i,j,k)*g4(i,j,k,itt) &
                       )

              Tmunu_ll(i,j,k,itx)= lorene2hydrobase*( &
                      ( energy_density(i,j,k) + pressure(i,j,k) ) &
                      *u_euler_l(i,j,k,it)*u_euler_l(i,j,k,ix) &
                      + pressure(i,j,k)*g4(i,j,k,itx) &
                       )

              Tmunu_ll(i,j,k,ity)= lorene2hydrobase*( &
                      ( energy_density(i,j,k) + pressure(i,j,k) ) &
                      *u_euler_l(i,j,k,it)*u_euler_l(i,j,k,iy) &
                      + pressure(i,j,k)*g4(i,j,k,ity) &
                       )

              Tmunu_ll(i,j,k,itz)= lorene2hydrobase*( &
                      ( energy_density(i,j,k) + pressure(i,j,k) ) &
                      *u_euler_l(i,j,k,it)*u_euler_l(i,j,k,iz) &
                      + pressure(i,j,k)*g4(i,j,k,itz) &
                       )

              Tmunu_ll(i,j,k,ixx)= lorene2hydrobase*( &
                      ( energy_density(i,j,k) + pressure(i,j,k) ) &
                      *u_euler_l(i,j,k,ix)*u_euler_l(i,j,k,ix) &
                      + pressure(i,j,k)*g4(i,j,k,ixx) &
                       )

              Tmunu_ll(i,j,k,ixy)= lorene2hydrobase*( &
                      ( energy_density(i,j,k) + pressure(i,j,k) ) &
                      *u_euler_l(i,j,k,ix)*u_euler_l(i,j,k,iy) &
                      + pressure(i,j,k)*g4(i,j,k,ixy) &
                       )

              Tmunu_ll(i,j,k,ixz)= lorene2hydrobase*( &
                      ( energy_density(i,j,k) + pressure(i,j,k) ) &
                      *u_euler_l(i,j,k,ix)*u_euler_l(i,j,k,iz) &
                      + pressure(i,j,k)*g4(i,j,k,ixz) &
                       )

              Tmunu_ll(i,j,k,iyy)= lorene2hydrobase*( &
                      ( energy_density(i,j,k) + pressure(i,j,k) ) &
                      *u_euler_l(i,j,k,iy)*u_euler_l(i,j,k,iy) &
                      + pressure(i,j,k)*g4(i,j,k,iyy)  &
                       )

              Tmunu_ll(i,j,k,iyz)= lorene2hydrobase*( &
                      ( energy_density(i,j,k) + pressure(i,j,k) ) &
                      *u_euler_l(i,j,k,iy)*u_euler_l(i,j,k,iz) &
                      + pressure(i,j,k)*g4(i,j,k,iyz) &
                       )

              Tmunu_ll(i,j,k,izz)= lorene2hydrobase*( &
                      ( energy_density(i,j,k) + pressure(i,j,k) ) &
                      *u_euler_l(i,j,k,iz)*u_euler_l(i,j,k,iz) &
                      + pressure(i,j,k)*g4(i,j,k,izz) &
                       )

            ! Print progress on screen
            perc= 100*(THIS% get_ngrid_x(l)*THIS% get_ngrid_y(l)*(k - 1) &
                  + THIS% get_ngrid_x(l)*(j - 1) + i) &
                  /( THIS% get_ngrid_x(l)* THIS% get_ngrid_y(l)* &
                     THIS% get_ngrid_z(l) )
            IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
              WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
                      creturn//" ", perc, "%"
            ENDIF

            ENDDO
          ENDDO
        ENDDO
      END ASSOCIATE
    ENDDO ref_levels2
    WRITE( *, "(A1)", ADVANCE= "NO" ) creturn
    PRINT *, " * Stress-energy tensor computed."
    PRINT *

    ! In debug mode, compute the Hamiltonian constraint by hand
    IF( debug )THEN

      DO l= 1, THIS% nlevels, 1

        ASSOCIATE( HC_rho         => HC_rho% levels(l)%var, &
                   HC_trK         => HC_trK% levels(l)%var, &
                   HC_A           => HC_A% levels(l)%var, &
                   HC_derphi      => HC_derphi% levels(l)%var, &
                   HC_hand        => HC_hand% levels(l)%var, &
                   phi            => THIS% phi% levels(l)%var, &
                   trK            => THIS% trK% levels(l)%var, &
                   A_BSSN3_ll     => THIS% A_BSSN3_ll% levels(l)%var, &
                   energy_density => energy_density% levels(l)% var, &
                   pressure       => pressure% levels(l)% var &
        )

          HC_rho= 0.0D0
          HC_trK= 0.0D0
          HC_A= 0.0D0
          HC_derphi= 0.0D0
          HC_hand= 0.0D0
          fd_lim= 5

          DO k= fd_lim, THIS% get_ngrid_z(l) - fd_lim, 1
            DO j= fd_lim, THIS% get_ngrid_y(l) - fd_lim, 1
              DO i= fd_lim, THIS% get_ngrid_x(l) - fd_lim, 1

                HC_rho( i, j, k )= 2.0D0*pi*EXP(5.0D0*phi( i, j, k )) &
                                      *lorene2hydrobase*energy_density( i, j, k )

                HC_trK( i, j, k )= - EXP(5.0D0*phi( i, j, k ))/12.0D0 &
                                        *trK( i, j, k )**2

                HC_A( i, j, k )= EXP(5.0D0*phi( i, j, k ))/8.0D0 &
                 *( A_BSSN3_ll(i, j, k,jxx)*A_BSSN3_ll(i, j, k,jxx) &
                  + A_BSSN3_ll(i, j, k,jxy)*A_BSSN3_ll(i, j, k,jxy) &
                  + A_BSSN3_ll(i, j, k,jxz)*A_BSSN3_ll(i, j, k,jxz) &
                  + A_BSSN3_ll(i, j, k,jxy)*A_BSSN3_ll(i, j, k,jxy) &
                  + A_BSSN3_ll(i, j, k,jyy)*A_BSSN3_ll(i, j, k,jyy) &
                  + A_BSSN3_ll(i, j, k,jyz)*A_BSSN3_ll(i, j, k,jyz) &
                  + A_BSSN3_ll(i, j, k,jxz)*A_BSSN3_ll(i, j, k,jxz) &
                  + A_BSSN3_ll(i, j, k,jyz)*A_BSSN3_ll(i, j, k,jyz) &
                  + A_BSSN3_ll(i, j, k,jzz)*A_BSSN3_ll(i, j, k,jzz) &
                  )

                ! Second derivative of conformal factor with fourth-order FD
                !HC_derphi( ix, iy, iz )= &
                !                   ( -      EXP(THIS% phi( ix + 2, iy, iz )) &
                !                     + 16.0*EXP(THIS% phi( ix + 1, iy, iz )) &
                !                     - 30.0*EXP(THIS% phi( ix    , iy, iz )) &
                !                     + 16.0*EXP(THIS% phi( ix - 1, iy, iz )) &
                !                     -      EXP(THIS% phi( ix - 2, iy, iz )) &
                !                     -      EXP(THIS% phi( ix, iy + 2, iz )) &
                !                     + 16.0*EXP(THIS% phi( ix, iy + 1, iz )) &
                !                     - 30.0*EXP(THIS% phi( ix, iy, iz )) &
                !                     + 16.0*EXP(THIS% phi( ix, iy - 1, iz )) &
                !                     -      EXP(THIS% phi( ix, iy - 2, iz )) &
                !                     -      EXP(THIS% phi( ix, iy, iz + 2 )) &
                !                     + 16.0*EXP(THIS% phi( ix, iy, iz + 1 )) &
                !                     - 30.0*EXP(THIS% phi( ix, iy, iz )) &
                !                     + 16.0*EXP(THIS% phi( ix, iy, iz - 1 )) &
                !                     -      EXP(THIS% phi( ix, iy, iz - 2 )) )&
                !                     /(12.0*THIS% dx**2)

                ! Second derivative of conformal factor with eighth-order FD
                HC_derphi( i, j, k )= ( &
                                - DBLE(1.0/560.0)*EXP(phi(i + 4, j, k)) &
                                + DBLE(8.0/315.0)*EXP(phi(i + 3, j, k)) &
                                - DBLE(1.0/5.0  )*EXP(phi(i + 2, j, k)) &
                                + DBLE(8.0/5.0  )*EXP(phi(i + 1, j, k)) &
                                - DBLE(205.0/72.0)*EXP(phi(i, j, k)) &
                                + DBLE(8.0/5.0  )*EXP(phi(i - 1, j, k)) &
                                - DBLE(1.0/5.0  )*EXP(phi(i - 2, j, k)) &
                                + DBLE(8.0/315.0)*EXP(phi(i - 3, j, k)) &
                                - DBLE(1.0/560.0)*EXP(phi(i - 4, j, k)) &
                                - DBLE(1.0/560.0)*EXP(phi(i, j + 4, k)) &
                                + DBLE(8.0/315.0)*EXP(phi(i, j + 3, k)) &
                                - DBLE(1.0/5.0  )*EXP(phi(i, j + 2, k)) &
                                + DBLE(8.0/5.0  )*EXP(phi(i, j + 1, k)) &
                                - DBLE(205.0/72.0)*EXP(phi(i, j, k)) &
                                + DBLE(8.0/5.0  )*EXP(phi(i, j - 1, k)) &
                                - DBLE(1.0/5.0  )*EXP(phi(i, j - 2, k)) &
                                + DBLE(8.0/315.0)*EXP(phi(i, j - 3, k)) &
                                - DBLE(1.0/560.0)*EXP(phi(i, j - 4, k)) &
                                - DBLE(1.0/560.0)*EXP(phi(i, j, k + 4)) &
                                + DBLE(8.0/315.0)*EXP(phi(i, j, k + 3)) &
                                - DBLE(1.0/5.0  )*EXP(phi(i, j, k + 2)) &
                                + DBLE(8.0/5.0  )*EXP(phi(i, j, k + 1)) &
                                - DBLE(205.0/72.0)*EXP(phi(i, j, k)) &
                                + DBLE(8.0/5.0  )*EXP(phi(i, j, k - 1)) &
                                - DBLE(1.0/5.0  )*EXP(phi(i, j, k - 2)) &
                                + DBLE(8.0/315.0)*EXP(phi(i, j, k - 3)) &
                                - DBLE(1.0/560.0)*EXP(phi(i, j, k - 4)) )&
                                /(THIS% levels(l)% dx**2)


                HC_hand( i, j, k )= HC_rho( i, j, k ) + &
                                    HC_trK( i, j, k ) + &
                                    HC_A( i, j, k )   + &
                                    HC_derphi( i, j, k )

              ENDDO
            ENDDO
          ENDDO
        END ASSOCIATE
      ENDDO

    ENDIF

    !
    !-- Compute the BSSN constraints by calling the Cactus-bound procedure
    !-- BSSN_CONSTRAINTS_INTERIOR
    !
    PRINT *, "** Computing contraints..."
    DO l= 1, THIS% nlevels, 1

      ASSOCIATE( lapse      => THIS% lapse% levels(l)% var, &
                 shift_u    => THIS% shift_u% levels(l)% var, &
                 phi        => THIS% phi% levels(l)% var, &
                 trK        => THIS% trK% levels(l)% var, &
                 g_BSSN3_ll => THIS% g_BSSN3_ll% levels(l)% var, &
                 A_BSSN3_ll => THIS% A_BSSN3_ll% levels(l)% var, &
                 Gamma_u    => THIS% Gamma_u% levels(l)% var, &
                 Tmunu_ll   => Tmunu_ll% levels(l)% var, &
                 HC         => THIS% HC% levels(l)% var, &
                 MC         => THIS% MC% levels(l)% var, &
                 GC         => THIS% GC% levels(l)% var &
      )

        imin(1) = THIS% levels(l)% nghost_x
        imin(2) = THIS% levels(l)% nghost_y
        imin(3) = THIS% levels(l)% nghost_z
        imax(1) = THIS% get_ngrid_x(l) - THIS% levels(l)% nghost_x - 1
        imax(2) = THIS% get_ngrid_y(l) - THIS% levels(l)% nghost_y - 1
        imax(3) = THIS% get_ngrid_z(l) - THIS% levels(l)% nghost_z - 1

        HC= 0.0D0
        MC= 0.0D0
        GC= 0.0D0
        CALL BSSN_CONSTRAINTS_INTERIOR( &
          !
          !-- Input
          !
          THIS% get_ngrid_x(l), THIS% get_ngrid_y(l), THIS% get_ngrid_z(l), &
          imin, imax, &
          THIS% get_dx(l), THIS% get_dy(l), THIS% get_dz(l), &
          g_BSSN3_ll(:,:,:,jxx), g_BSSN3_ll(:,:,:,jxy), &
          g_BSSN3_ll(:,:,:,jxz), g_BSSN3_ll(:,:,:,jyy), &
          g_BSSN3_ll(:,:,:,jyz), g_BSSN3_ll(:,:,:,jzz), &
          A_BSSN3_ll(:,:,:,jxx), A_BSSN3_ll(:,:,:,jxy), &
          A_BSSN3_ll(:,:,:,jxz), A_BSSN3_ll(:,:,:,jyy), &
          A_BSSN3_ll(:,:,:,jyz), A_BSSN3_ll(:,:,:,jzz), &
          trK(:,:,:), phi(:,:,:), &
          Gamma_u(:,:,:,jx), &
          Gamma_u(:,:,:,jy), &
          Gamma_u(:,:,:,jz), &
          Tmunu_ll(:,:,:,itt), &
          Tmunu_ll(:,:,:,itx), &
          Tmunu_ll(:,:,:,ity), &
          Tmunu_ll(:,:,:,itz), &
          Tmunu_ll(:,:,:,ixx), &
          Tmunu_ll(:,:,:,ixy), &
          Tmunu_ll(:,:,:,ixz), &
          Tmunu_ll(:,:,:,iyy), &
          Tmunu_ll(:,:,:,iyz), &
          Tmunu_ll(:,:,:,izz), &
          lapse(:,:,:), &
          shift_u(:,:,:,jx), &
          shift_u(:,:,:,jy), &
          shift_u(:,:,:,jz), &
          !
          !-- Output
          !
          ! Connection constraints
          GC(:,:,:,jx), GC(:,:,:,jy), &
          GC(:,:,:,jz), &
          ! Hamiltonian and momentum constraints
          HC(:,:,:), &
          MC(:,:,:,jx), &
          MC(:,:,:,jy), &
          MC(:,:,:,jz) &
        )
      END ASSOCIATE
    ENDDO
    PRINT *, " * Constraints computed."
    PRINT *

    !---------------------------------------------------------!
    !--  Analyze constraints, and print to formatted files  --!
    !---------------------------------------------------------!

    DO l= 1, THIS% nlevels, 1

      ASSOCIATE( HC         => THIS% HC% levels(l)% var, &
                 MC         => THIS% MC% levels(l)% var, &
                 GC         => THIS% GC% levels(l)% var &
      )
        !
        !-- Export the constraint statistics to a formatted file
        !
        unit_logfile= 2891

        IF( l > 9 )THEN
          WRITE( n_reflev, "(I2)" ) l
        ELSE
          WRITE( n_reflev, "(I1)" ) l
        ENDIF

        finalname_logfile= TRIM(name_logfile)//"-reflev"//TRIM(n_reflev)//".log"

        INQUIRE( FILE= TRIM(finalname_logfile), EXIST= exist )

        IF( exist )THEN
            OPEN( UNIT= unit_logfile, FILE= TRIM(finalname_logfile), &
                  STATUS= "REPLACE", &
                  FORM= "FORMATTED", &
                  POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
                  IOMSG= err_msg )
        ELSE
            OPEN( UNIT= unit_logfile, FILE= TRIM(finalname_logfile), &
                  STATUS= "NEW", &
                  FORM= "FORMATTED", &
                  ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
        ENDIF
        IF( ios > 0 )THEN
          PRINT *, "...error when opening ", TRIM(finalname_logfile), &
                   ". The error message is", err_msg
          STOP
        ENDIF
        !CALL test_status( ios, err_msg, "...error when opening " &
        !         // TRIM(name_logfile) )

        IF( .NOT.ALLOCATED( THIS% HC_l2 ))THEN
          ALLOCATE( THIS% HC_l2( THIS% nlevels ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array HC_l2. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF
        IF( .NOT.ALLOCATED( THIS% MC_l2 ))THEN
          ALLOCATE( THIS% MC_l2( THIS% nlevels, 3 ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array MC_l2. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF
        IF( .NOT.ALLOCATED( THIS% GC_l2 ))THEN
          ALLOCATE( THIS% GC_l2( THIS% nlevels, 3 ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array GC_l2. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF
        IF( .NOT.ALLOCATED( THIS% HC_loo ))THEN
          ALLOCATE( THIS% HC_loo( THIS% nlevels ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array HC_loo. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF
        IF( .NOT.ALLOCATED( THIS% MC_loo ))THEN
          ALLOCATE( THIS% MC_loo( THIS% nlevels, 3 ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array MC_loo. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF
        IF( .NOT.ALLOCATED( THIS% GC_loo ))THEN
          ALLOCATE( THIS% GC_loo( THIS% nlevels, 3 ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array GC_loo. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF

        WRITE( UNIT = unit_logfile, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id

        PRINT *, "** Analyzing constraints on refinement level ", l, "..."

        name_analysis= "bssn-hc-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the Hamiltonian constraint"
        CALL THIS% analyze_constraint( &
             l, &
             HC, name_constraint, unit_logfile, name_analysis, &
             THIS% HC_l2(l), THIS% HC_loo(l) )

        name_analysis= "bssn-mc1-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the first component of the momentum constraint"
        CALL THIS% analyze_constraint( &
             l, &
             MC(:,:,:,jx), name_constraint, unit_logfile, name_analysis, &
             THIS% MC_l2(l,jx), THIS% MC_loo(l,jx) )

        name_analysis= "bssn-mc2-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the second component of the momentum constraint"
        CALL THIS% analyze_constraint( &
             l, &
             MC(:,:,:,jy), name_constraint, unit_logfile, name_analysis, &
             THIS% MC_l2(l,jy), THIS% MC_loo(l,jy) )

        name_analysis= "bssn-mc3-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the third component of the momentum constraint"
        CALL THIS% analyze_constraint( &
             l, &
             MC(:,:,:,jz), name_constraint, unit_logfile, name_analysis, &
             THIS% MC_l2(l,jz), THIS% MC_loo(l,jz) )

        name_analysis= "bssn-gc1-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the first component of the connection constraint"
        CALL THIS% analyze_constraint( &
             l, &
             GC(:,:,:,jx), name_constraint, unit_logfile, name_analysis, &
             THIS% GC_l2(l,jx), THIS% GC_loo(l,jx) )

        name_analysis= "bssn-gc2-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the second component of the connection constraint"
        CALL THIS% analyze_constraint( &
             l, &
             GC(:,:,:,jy), name_constraint, unit_logfile, name_analysis, &
             THIS% GC_l2(l,jy), THIS% GC_loo(l,jy) )

        name_analysis= "bssn-gc3-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the third component of the connection constraint"
        CALL THIS% analyze_constraint( &
             l, &
             GC(:,:,:,jz), name_constraint, unit_logfile, name_analysis, &
             THIS% GC_l2(l,jz), THIS% GC_loo(l,jz) )

        CLOSE( UNIT= unit_logfile )

      PRINT *, " * Constraints analyzed. Summary of results saved to ", &
               finalname_logfile
      PRINT *
      END ASSOCIATE
    ENDDO

    IF( THIS% export_constraints )THEN

      PRINT *, "** Printing constraints to file ", TRIM(namefile), "..."

      !
      !-- Export the constraints to a formatted file
      !
      INQUIRE( FILE= TRIM(namefile), EXIST= exist )

      IF( exist )THEN
        OPEN( UNIT= 20, FILE= TRIM(namefile), STATUS= "REPLACE", &
              FORM= "FORMATTED", &
              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
              IOMSG= err_msg )
      ELSE
        OPEN( UNIT= 20, FILE= TRIM(namefile), STATUS= "NEW", &
        FORM= "FORMATTED", &
              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
      ENDIF
      IF( ios > 0 )THEN
        PRINT *, "...error when opening ", TRIM(namefile), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when opening " &
      !         // TRIM(namefile) )

      WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id
      WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "# Values of the stress-energy tensor and the BSSN constraints" &
      // " for the LORENE ID " &
      // "on selected grid points"
      IF( ios > 0 )THEN
        PRINT *, "...error when writing line 1 in ", TRIM(namefile), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when writing line 1 in "&
      !         // TRIM(namefile) )
      WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "# column:      1        2       3       4       5", &
      "       6       7       8       9       10", &
      "       11       12       13       14       15", &
      "       16       17       18       19       20"
      IF( ios > 0 )THEN
        PRINT *, "...error when writing line 2 in ", TRIM(namefile), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when writing line 2 in "&
      !        // TRIM(namefile) )
      WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "#      refinement level    x   y   z   Stress-energy (10 components)   "&
      // "Hamiltonian constraint       " &
      // "Momentum constraint (three components)       " &
      // "Connection constraint (three components)"
      IF( ios > 0 )THEN
        PRINT *, "...error when writing line 3 in ", TRIM(namefile), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when writing line 3 in "&
      !        // TRIM(namefile) )

      DO l= 1, THIS% nlevels, 1

        ASSOCIATE( lapse           => THIS% lapse% levels(l)% var, &
                   shift_u         => THIS% shift_u% levels(l)% var, &
                   phi             => THIS% phi% levels(l)% var, &
                   trK             => THIS% trK% levels(l)% var, &
                   g_BSSN3_ll      => THIS% g_BSSN3_ll% levels(l)% var, &
                   A_BSSN3_ll      => THIS% A_BSSN3_ll% levels(l)% var, &
                   g_phys3_ll      => THIS% g_phys3_ll% levels(l)% var, &
                   k_phys3_ll      => THIS% k_phys3_ll% levels(l)% var, &
                   Gamma_u         => THIS% Gamma_u% levels(l)% var, &
                   Tmunu_ll        => Tmunu_ll% levels(l)% var, &
                   v_euler_l       => v_euler_l% levels(l)% var, &
                   u_euler_l       => u_euler_l% levels(l)% var, &
                   v_euler         => v_euler% levels(l)% var, &
                   lorentz_factor  => lorentz_factor% levels(l)% var, &
                   HC              => THIS% HC% levels(l)% var, &
                   MC              => THIS% MC% levels(l)% var, &
                   GC              => THIS% GC% levels(l)% var, &
                   HC_rho          => HC_rho% levels(l)%var, &
                   HC_trK          => HC_trK% levels(l)%var, &
                   HC_A            => HC_A% levels(l)%var, &
                   HC_derphi       => HC_derphi% levels(l)%var, &
                   HC_hand         => HC_hand% levels(l)%var, &
                   g4              => g4% levels(l)% var, &
                   baryon_density  => baryon_density% levels(l)% var, &
                   specific_energy => specific_energy% levels(l)% var, &
                   energy_density  => energy_density% levels(l)% var, &
                   pressure        => pressure% levels(l)% var &
        )

          ! Being abs_grid a local array, it is good practice to allocate it on
          ! the heap, otherwise it will be stored on the stack which has a very
          ! limited size. This results in a segmentation fault.
          IF( ALLOCATED( abs_grid ) )THEN
            DEALLOCATE( abs_grid )
          ENDIF
          ALLOCATE( abs_grid( THIS% get_ngrid_x(l), THIS% get_ngrid_y(l), &
                              THIS% get_ngrid_z(l), 3 ) )

          DO k= 1, THIS% get_ngrid_z(l), 1
            DO j= 1, THIS% get_ngrid_y(l), 1
              DO i= 1, THIS% get_ngrid_x(l), 1

                abs_grid( i, j, k, jx )= &
                            ABS( THIS% coords% levels(l)% var( i, j, k, jx ) )
                abs_grid( i, j, k, jy )= &
                            ABS( THIS% coords% levels(l)% var( i, j, k, jy ) )
                abs_grid( i, j, k, jz )= &
                            ABS( THIS% coords% levels(l)% var( i, j, k, jz ) )

              ENDDO
            ENDDO
          ENDDO

          min_abs_y= 1D+20
          min_abs_z= 1D+20
          DO k= 1, THIS% get_ngrid_z(l), 1
            DO j= 1, THIS% get_ngrid_y(l), 1
              DO i= 1, THIS% get_ngrid_x(l), 1

                IF( ABS( THIS% coords% levels(l)% var( i, j, k, jy ) ) &
                    < min_abs_y )THEN
                  min_abs_y= ABS( THIS% coords% levels(l)% var( i, j, k, jy ) )
                  min_ix_y= i
                  min_iy_y= j
                  min_iz_y= k
                ENDIF

                IF( ABS( THIS% coords% levels(l)% var( i, j, k, jz ) ) &
                    < min_abs_z )THEN
                  min_abs_z= ABS( THIS% coords% levels(l)% var( i, j, k, jz ) )
                  min_ix_z= i
                  min_iy_z= j
                  min_iz_z= k
                ENDIF

              ENDDO
            ENDDO
          ENDDO

          DO k= 1, THIS% get_ngrid_z(l), 1

            IF( MOD( k, THIS% cons_step ) /= 0 ) CYCLE

            DO j= 1, THIS% get_ngrid_y(l), 1

              IF( MOD( j, THIS% cons_step ) /= 0 ) CYCLE

              DO i= 1, THIS% get_ngrid_x(l), 1

                IF( MOD( i, THIS% cons_step ) /= 0 ) CYCLE

                IF( THIS% export_constraints_xy .AND. &
                    ( THIS% coords% levels(l)% var( i, j, k, jz ) /= &
                      THIS% coords% levels(l)% var( min_ix_z, min_iy_z, &
                                                    min_iz_z, jz ) ) )THEN
                  CYCLE
                ENDIF
                IF( THIS% export_constraints_x .AND. &
                    ( THIS% coords% levels(l)% var( i, j, k, jz ) /= &
                      THIS% coords% levels(l)% var( min_ix_z, min_iy_z, &
                                                    min_iz_z, jz ) &
                      .OR. &
                      THIS% coords% levels(l)% var( i, j, k, jy ) /= &
                      THIS% coords% levels(l)% var( min_ix_y, min_iy_y, &
                                                    min_iz_y, jy ) ) )THEN
                  CYCLE
                ENDIF

                IF( debug )THEN
                  WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, &
                           FMT = * )&
                    l, &
                    THIS% coords% levels(l)% var( i, j, k, jx ), &
                    THIS% coords% levels(l)% var( i, j, k, jy ), &
                    THIS% coords% levels(l)% var( i, j, k, jz ), &
                    baryon_density( i, j, k ), &
                    energy_density( i, j, k ), &
                    specific_energy( i, j, k ), &
                    pressure( ix, iy, iz ), &
                    u_euler_l( i, j, k, it ), &
                    u_euler_l( i, j, k, ix ), &
                    u_euler_l( i, j, k, iy ), &
                    u_euler_l( i, j, k, iz ), &
                    u_euler_l( i, j, k, it ), &
                    u_euler_l( i, j, k, ix ), &
                    u_euler_l( i, j, k, iy ), &
                    u_euler_l( i, j, k, iz ), &
                    v_euler( i, j, k, jx ), &
                    v_euler( i, j, k, jy ), &
                    v_euler( i, j, k, jz ), &
                    Tmunu_ll( i, j, k, itt ), &
                    Tmunu_ll( i, j, k, itx ), &
                    Tmunu_ll( i, j, k, ity ), &
                    Tmunu_ll( i, j, k, itz ), &
                    Tmunu_ll( i, j, k, ixx ), &
                    Tmunu_ll( i, j, k, ixy ), &
                    Tmunu_ll( i, j, k, ixz ), &
                    Tmunu_ll( i, j, k, iyy ), &
                    Tmunu_ll( i, j, k, iyz ), &
                    Tmunu_ll( i, j, k, izz ), &
                    HC( i, j, k ), &
                    HC_hand( i, j, k ), &
                    HC_rho( i, j, k ), &
                    HC_trK( i, j, k ), &
                    HC_A( i, j, k ), &
                    HC_derphi( i, j, k ), &
                    lorentz_factor( i, j, k ), &
                    lapse( ix, iy, iz ), &
                    shift_u( i, j, k, jx ), &
                    shift_u( i, j, k, jy ), &
                    shift_u( i, j, k, jz ), &
                    g4( i, j, k, ixx ), &
                    g4( i, j, k, ixy ), &
                    g4( i, j, k, ixz ), &
                    g4( i, j, k, iyy ), &
                    g4( i, j, k, iyz ), &
                    g4( i, j, k, izz ), &
                    !g_BSSN3_ll( i, j, k, jxx ), &
                    !g_BSSN3_ll( i, j, k, jxy ), &
                    !g_BSSN3_ll( i, j, k, jxz ), &
                    !g_BSSN3_ll( i, j, k, jyy ), &
                    !g_BSSN3_ll( i, j, k, jyz ), &
                    !g_BSSN3_ll( i, j, k, jzz ), &
                    k_phys3_ll( i, j, k, jxx ), &
                    k_phys3_ll( i, j, k, jxy ), &
                    k_phys3_ll( i, j, k, jxz ), &
                    k_phys3_ll( i, j, k, jyy ), &
                    k_phys3_ll( i, j, k, jyz ), &
                    k_phys3_ll( i, j, k, jzz ), &
                    A_BSSN3_ll( i, j, k, jxx ), &
                    A_BSSN3_ll( i, j, k, jxy ), &
                    A_BSSN3_ll( i, j, k, jxz ), &
                    A_BSSN3_ll( i, j, k, jyy ), &
                    A_BSSN3_ll( i, j, k, jyz ), &
                    A_BSSN3_ll( i, j, k, jzz ), &
                    trK( i, j, k ), &
                    phi( i, j, k ), &
                    Gamma_u( i, j, k, 1 ), &
                    Gamma_u( i, j, k, 2 ), &
                    Gamma_u( i, j, k, 3 )
                ELSE
                  WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * )&
                    l, &
                    THIS% coords% levels(l)% var( i, j, k, jx ), &
                    THIS% coords% levels(l)% var( i, j, k, jy ), &
                    THIS% coords% levels(l)% var( i, j, k, jz ), &
                    Tmunu_ll( i, j, k, itt ), &
                    Tmunu_ll( i, j, k, itx ), &
                    Tmunu_ll( i, j, k, ity ), &
                    Tmunu_ll( i, j, k, itz ), &
                    Tmunu_ll( i, j, k, ixx ), &
                    Tmunu_ll( i, j, k, ixy ), &
                    Tmunu_ll( i, j, k, ixz ), &
                    Tmunu_ll( i, j, k, iyy ), &
                    Tmunu_ll( i, j, k, iyz ), &
                    Tmunu_ll( i, j, k, izz ), &
                    HC( i, j, k ), &
                    MC( i, j, k, jx ), &
                    MC( i, j, k, jy ), &
                    MC( i, j, k, jz ), &
                    GC( i, j, k, jx ), &
                    GC( i, j, k, jy ), &
                    GC( i, j, k, jz )
                ENDIF

                IF( ios > 0 )THEN
                  PRINT *, "...error when writing the arrays in ", &
                           TRIM(namefile), &
                           ". The error message is", err_msg
                  STOP
                ENDIF
                !CALL test_status( ios, err_msg, &
                !          "...error in writing " &
                !          // "the arrays in " // TRIM(namefile) )
              ENDDO
            ENDDO
          ENDDO
        END ASSOCIATE
      ENDDO

      CLOSE( UNIT= 20 )

      PRINT *, " * Printed."
      PRINT *

    ENDIF

    !DEALLOCATE( baryon_density )
    !DEALLOCATE( energy_density )
    !DEALLOCATE( specific_energy )
    !DEALLOCATE( pressure )
    !DEALLOCATE( v_euler )
    !!DEALLOCATE( u_coord )
    !!DEALLOCATE( u_coord_l )
    !DEALLOCATE( g4 )
    !DEALLOCATE( g4temp )
    !DEALLOCATE( ig4 )
    !DEALLOCATE( Tmunu_ll )

  END PROCEDURE compute_and_export_bssn_constraints_grid


  MODULE PROCEDURE compute_and_export_bssn_constraints_particles

    !**************************************************
    !                                                 *
    ! Compute, store and export the BSSN constraints  *
    ! to a formatted file. The computaton is done     *
    ! mapping the physical metric from the gravity    *
    ! to the particles, computing e stress-energy     *
    ! tensor on the particles, and mapping it to the  *
    ! gravity grid.                                   *
    ! TODO: use the SPH density to compute the        *
    !       stress-energy tensor, rather than the     *
    !       LORENE density                            *
    !                                                 *
    ! FT 1.02.2021                                    *
    !                                                 *
    !**************************************************

    USE constants,            ONLY: c_light2, cm2m, MSun, g2kg, m2cm, Msun_geo
    USE units,                ONLY: set_units, m0c2_cu
    USE tensor,               ONLY: itt, itx, ity, itz, ixx, ixy, &
                                    ixz, iyy, iyz, izz, jxx, jxy, jxz, &
                                    jyy, jyz, jzz, jx, jy, jz, &
                                    n_sym3x3, n_sym4x4

    USE mesh_refinement,             ONLY: allocate_grid_function, levels, &
                                           rad_coord, deallocate_grid_function
    USE ADM_refine,                  ONLY: lapse, dt_lapse, shift_u, &
                                           dt_shift_u, &
                                           K_phys3_ll, g_phys3_ll, &
                                           allocate_ADM, deallocate_ADM
    USE BSSN_refine,                 ONLY: allocate_BSSN, deallocate_BSSN
    USE Tmunu_refine,                ONLY: Tmunu_ll, allocate_Tmunu, &
                                           deallocate_Tmunu
    USE McLachlan_refine,            ONLY: BSSN_CONSTRAINTS_INTERIOR, &
                                           allocate_Ztmp, deallocate_Ztmp
    USE GravityAcceleration_refine,  ONLY: allocate_GravityAcceleration, &
                                           deallocate_GravityAcceleration
    USE Extract_Mass,                ONLY: radius2


    USE input_output,         ONLY: read_options
    USE options,              ONLY: ikernel, ndes, metric_type
    USE sph_variables,        ONLY: npart, &  ! particle number
                                    pos_u, &  ! particle positions
                                    vel_u, &  ! particle velocities in
                                              ! coordinate frame
                                    nlrf,  &  ! baryon number density in
                                              ! local rest frame
                                    ehat,  &  ! canonical energy per baryon
                                    nu,    &  ! canonical baryon number per
                                              ! particle
                                    Theta, &  ! Generalized Lorentz factor
                                    h,     &  ! Smoothing length
                                    Pr,    &  ! Pressure
                                    u,     &  ! Internal energy in local rest
                                              ! frame (no kinetic energy)
                                    temp,  &  ! Temperature
                                    av,    &  ! Dissipation
                                    Ye,    &  ! Electron fraction
                                    divv,  &  ! Divergence of velocity vel_u
                                    Nstar, &  ! Comput.frame baryon number
                                              ! density
                                    allocate_SPH_memory, &
                                    deallocate_SPH_memory
    USE RCB_tree_3D,          ONLY: allocate_RCB_tree_memory_3D,&
                                    deallocate_RCB_tree_memory_3D, iorig
    USE kernel_table,         ONLY: ktable
    USE gradient,             ONLY: allocate_gradient, deallocate_gradient
    USE sphincs_sph,          ONLY: density, flag_dead_ll_cells
    USE set_h,                ONLY: exact_nei_tree_update
    USE alive_flag,           ONLY: alive

    USE map_particles_2_grid, ONLY: map_2_grid
    USE metric_on_particles,  ONLY: allocate_metric_on_particles, &
                                    deallocate_metric_on_particles, &
                                    get_metric_on_particles
    USE particle_mesh,        ONLY: deallocate_all_lists, &
                                    deallocate_flag_nei_cell, &
                                    deallocate_pp_g

    IMPLICIT NONE

    INTEGER:: i, j, k, l, a, allocation_status
    INTEGER, DIMENSION(3) :: imin, imax
    INTEGER:: unit_logfile, min_ix_y, min_iy_y, min_iz_y, &
              min_ix_z, min_iy_z, min_iz_z
    INTEGER:: itr, min_y_index
    INTEGER, SAVE:: counter= 1


    DOUBLE PRECISION:: min_abs_y, min_abs_z
    DOUBLE PRECISION, DIMENSION( :, :, :, : ), ALLOCATABLE:: abs_grid
    DOUBLE PRECISION, DIMENSION( :, : ), ALLOCATABLE:: abs_pos

    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nlrf_loc
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nu_loc
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: u_loc
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_loc
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos_loc
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: vel_loc
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: theta_loc
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: sph_density

    CHARACTER( LEN= : ), ALLOCATABLE:: name_constraint
    CHARACTER( LEN= : ), ALLOCATABLE:: name_analysis
    CHARACTER( LEN= : ), ALLOCATABLE:: finalname_logfile
    CHARACTER( 2 ):: n_reflev

    LOGICAL:: exist
    LOGICAL, PARAMETER:: debug= .FALSE.

    levels = THIS% levels
    ! radius2 is the extraction radius. If not set here, then it is 0 by default
    ! and the metric is not interpolate on the particle in
    ! get_metric_on_particles
    radius2= HUGE(DBLE(1.0D0))

    CALL allocate_grid_function( THIS% HC_parts, "HC_parts_ID", 1 )
    CALL allocate_grid_function( THIS% MC_parts, "MC_parts_ID", 3 )
    CALL allocate_grid_function( THIS% GC_parts, "GC_parts_ID", 3 )

    PRINT *, "Mapping hydro fields from particles to grid..."

    CALL allocate_ADM()
    CALL allocate_BSSN()

    ! Allocate temporary memory for time integration
    CALL allocate_Ztmp()

    ! Allocate memory for the stress-energy tensor (used in
    ! write_BSSN_dump)
    CALL allocate_Tmunu()

    ! Allocate memory for the derivatives of the ADM variables
    CALL allocate_GravityAcceleration()

    CALL allocate_grid_function( rad_coord, 'rad_coord', 1 )

    ! Initialize the stress-energy tensor to 0
    DO l= 1, THIS% nlevels, 1
      Tmunu_ll%   levels(l)% var= 0.0D0
      rad_coord%  levels(l)% var= THIS% rad_coord%  levels(l)% var
      g_phys3_ll% levels(l)% var= THIS% g_phys3_ll% levels(l)% var
      shift_u%    levels(l)% var= THIS% shift_u%    levels(l)% var
      lapse%      levels(l)% var= THIS% lapse%      levels(l)% var
    ENDDO

    !
    !-- Allocate memory for the coordinate radius
    !-- (this is needed by ADM_to_BSSN)
    !
    !IF( .NOT. ALLOCATED( rad_coord ) )THEN
    !    ALLOCATE( rad_coord( THIS% ngrid_x, THIS% ngrid_y, &
    !                         THIS% ngrid_z ), STAT= allocation_status )
    !ENDIF
    !IF( allocation_status > 0 )THEN
    !   PRINT *, '...allocation error for rad_coord'
    !   STOP
    !ENDIF
    !
    !DO iz= 1, THIS% ngrid_z
    !  DO iy= 1, THIS% ngrid_y
    !    DO ix= 1, THIS% ngrid_x
    !      rad_coord( ix, iy, iz )= SQRT( (xL+(ix-1)*dx)**2 &
    !                                   + (yL+(iy-1)*dy)**2 &
    !                                   + (zL+(iz-1)*dz)**2 )
    !    ENDDO
    !  ENDDO
    !ENDDO

    !dt_lapse  = 0.0D0
    !dt_shift_u= 0.0D0
    !lapse= THIS% lapse
    !shift_u= THIS% shift_u
    !g_phys3_ll= THIS% g_phys3_ll

    npart= parts_obj% get_npart()

    IF( .NOT. ALLOCATED( nlrf_loc ) )THEN
        ALLOCATE( nlrf_loc( npart ), STAT= allocation_status )
    ENDIF
    IF( allocation_status > 0 )THEN
       PRINT *, '...allocation error for nlrf_loc'
       STOP
    ENDIF
    IF( .NOT. ALLOCATED( nu_loc ) )THEN
        ALLOCATE( nu_loc( npart ), STAT= allocation_status )
    ENDIF
    IF( allocation_status > 0 )THEN
       PRINT *, '...allocation error for nu_loc'
       STOP
    ENDIF
    IF( .NOT. ALLOCATED( u_loc ) )THEN
        ALLOCATE( u_loc( npart ), STAT= allocation_status )
    ENDIF
    IF( allocation_status > 0 )THEN
       PRINT *, '...allocation error for u_loc'
       STOP
    ENDIF
    IF( .NOT. ALLOCATED( pressure_loc ) )THEN
        ALLOCATE( pressure_loc( npart ), STAT= allocation_status )
    ENDIF
    IF( allocation_status > 0 )THEN
       PRINT *, '...allocation error for pressure_loc'
       STOP
    ENDIF
    IF( .NOT. ALLOCATED( theta_loc ) )THEN
        ALLOCATE( theta_loc( npart ), STAT= allocation_status )
    ENDIF
    IF( allocation_status > 0 )THEN
       PRINT *, '...allocation error for theta_loc'
       STOP
    ENDIF
    IF( .NOT. ALLOCATED( pos_loc ) )THEN
        ALLOCATE( pos_loc( 3, npart ), STAT= allocation_status )
    ENDIF
    IF( allocation_status > 0 )THEN
       PRINT *, '...allocation error for pos_loc'
       STOP
    ENDIF
    IF( .NOT. ALLOCATED( vel_loc ) )THEN
        ALLOCATE( vel_loc( 3, npart ), STAT= allocation_status )
    ENDIF
    IF( allocation_status > 0 )THEN
       PRINT *, '...allocation error for vel_loc'
       STOP
    ENDIF
    IF( .NOT. ALLOCATED( sph_density ) )THEN
        ALLOCATE( sph_density( npart ), STAT= allocation_status )
    ENDIF
    IF( allocation_status > 0 )THEN
       PRINT *, '...allocation error for sph_density'
       STOP
    ENDIF

    !IF( .NOT. ALLOCATED( tmp ) )THEN
    !    ALLOCATE( tmp( npart ), STAT= allocation_status )
    !ENDIF
    !IF( allocation_status > 0 )THEN
    !   PRINT *, '...allocation error for tmp'
    !   STOP
    !ENDIF
    !IF( .NOT. ALLOCATED( tmp2 ) )THEN
    !    ALLOCATE( tmp2( 3, npart ), STAT= allocation_status )
    !ENDIF
    !IF( allocation_status > 0 )THEN
    !   PRINT *, '...allocation error for tmp2'
    !   STOP
    !ENDIF

    ! Set the SPH density to 0 by default
    sph_density= 0.0D0

    CALL set_units('NSM')
    CALL read_options

    CALL allocate_SPH_memory

    IF( debug ) PRINT *, "-2"

    CALL allocate_RCB_tree_memory_3D(npart)
    iorig(1:npart)= (/ (a,a=1,npart) /)

    IF( debug ) PRINT *, "-1"

    h= parts_obj% get_h()

    IF( counter == 1 )THEN
      ! tabulate kernel, get ndes
      CALL ktable(ikernel,ndes)
    ENDIF

    IF( debug ) PRINT *, "0"

    nu_loc      = parts_obj% get_nu()
    pos_loc     = parts_obj% get_pos()
    vel_loc     = parts_obj% get_vel()
    u_loc       = parts_obj% get_u()
    nlrf_loc    = parts_obj% get_nlrf()
    theta_loc   = parts_obj% get_theta()
    pressure_loc= parts_obj% get_pressure_cu()

    IF( debug ) PRINT *, "1"

    PRINT *, " * Allocating needed memory..."
    PRINT *

    ! flag that particles are 'alive'
    ALLOCATE( alive( npart ) )
    alive( 1:npart )= 1

    CALL allocate_gradient( npart )

    IF( debug ) PRINT *, "2"

    CALL allocate_metric_on_particles( npart )

    IF( debug ) PRINT *, "3"

    !---------------------------!
    !--  Compute constraints  --!
    !---------------------------!

    PRINT *, " * Mapping metric from the grid to the particles..."
    PRINT *
    CALL get_metric_on_particles( npart, &
                                  pos_loc )

    IF( debug ) PRINT *, "4"

    !
    !-- Seems like computing neighbors and SPH density is not needed to map
    !-- the stress-energy tensor from the particles to the grid
    !

    PRINT *, " * Computing neighbours..."
    PRINT *
    CALL exact_nei_tree_update( ndes,    &
                                npart,   &
                                pos_loc, &
                                nu_loc )

    IF( debug ) PRINT *, "5"

    PRINT *, " * Computing SPH density..."
    PRINT *
    nu = nu_loc
    pos_u = pos_loc
    vel_u = vel_loc
    u = u_loc
    nlrf = nlrf_loc
    Theta = theta_loc
    Pr = pressure_loc
    CALL density( npart,   &
                  pos_loc, &
                  sph_density )

    IF( debug ) PRINT *, "6"

    IF( debug .AND. .TRUE. ) PRINT *, "npart= ", npart
    IF( debug .AND. .TRUE. ) PRINT *, "nu_loc= ", nu_loc(npart/2)
    IF( debug .AND. .TRUE. ) PRINT *, "pos_loc= ", pos_loc(2,npart/2)
    IF( debug .AND. .TRUE. ) PRINT *, "vel_loc= ", vel_loc(2,npart/2)
    IF( debug .AND. .TRUE. ) PRINT *, "u_loc= ", u_loc(npart/2)
    IF( debug .AND. .TRUE. ) PRINT *, "nlrf_loc= ", nlrf_loc(npart/2)
    IF( debug .AND. .TRUE. ) PRINT *, "theta_loc= ", theta_loc(npart/2)
    IF( debug .AND. .TRUE. ) PRINT *, "pressure_loc= ", pressure_loc(npart/2)
    IF( debug .AND. .TRUE. ) PRINT *

    PRINT *, " * Mapping stress-energy tensor from the particles to the grid..."
    PRINT *
    CALL map_2_grid( npart        , &
                     nu_loc       , &
                     pos_loc      , &
                     vel_loc      , &
                     u_loc        , &
                     nlrf_loc     , &
                     theta_loc    , &
                     pressure_loc )

    IF( debug ) PRINT *, "7"

    !
    !-- Deallocate SPH MODULE variables
    !
    CALL deallocate_grid_function ( rad_coord, 'rad_coord' )
    CALL deallocate_flag_nei_cell
    CALL deallocate_pp_g
    CALL deallocate_all_lists
    CALL deallocate_metric_on_particles
    CALL deallocate_gradient
    DEALLOCATE( alive )
    !DEALLOCATE(W_no_norm)
    !DEALLOCATE(dWdv_no_norm)
    !DEALLOCATE(fmass)
    !DEALLOCATE(fpoten)
    !DEALLOCATE(dphidh)
    CALL deallocate_RCB_tree_memory_3D
    CALL deallocate_SPH_memory

    !
    !-- Compute the BSSN constraints by calling the Cactus-bound procedure
    !-- BSSN_CONSTRAINTS_INTERIOR
    !
    PRINT *, " * Computing contraints using particle data..."
    DO l= 1, THIS% nlevels, 1

      ASSOCIATE( lapse      => THIS% lapse% levels(l)% var, &
                 shift_u    => THIS% shift_u% levels(l)% var, &
                 phi        => THIS% phi% levels(l)% var, &
                 trK        => THIS% trK% levels(l)% var, &
                 g_BSSN3_ll => THIS% g_BSSN3_ll% levels(l)% var, &
                 A_BSSN3_ll => THIS% A_BSSN3_ll% levels(l)% var, &
                 Gamma_u    => THIS% Gamma_u% levels(l)% var, &
                 Tmunu_ll   => Tmunu_ll% levels(l)% var, &
                 HC_parts   => THIS% HC_parts % levels(l)% var, &
                 MC_parts   => THIS% MC_parts % levels(l)% var, &
                 GC_parts   => THIS% GC_parts % levels(l)% var &
      )

        imin(1) = THIS% levels(l)% nghost_x
        imin(2) = THIS% levels(l)% nghost_y
        imin(3) = THIS% levels(l)% nghost_z
        imax(1) = THIS% get_ngrid_x(l) - THIS% levels(l)% nghost_x - 1
        imax(2) = THIS% get_ngrid_y(l) - THIS% levels(l)% nghost_y - 1
        imax(3) = THIS% get_ngrid_z(l) - THIS% levels(l)% nghost_z - 1

        HC_parts= 0.0D0
        MC_parts= 0.0D0
        GC_parts= 0.0D0
        CALL BSSN_CONSTRAINTS_INTERIOR( &
          !
          !-- Input
          !
          THIS% get_ngrid_x(l), THIS% get_ngrid_y(l), THIS% get_ngrid_z(l), &
          imin, imax, &
          THIS% get_dx(l), THIS% get_dy(l), THIS% get_dz(l), &
          g_BSSN3_ll(:,:,:,jxx), g_BSSN3_ll(:,:,:,jxy), &
          g_BSSN3_ll(:,:,:,jxz), g_BSSN3_ll(:,:,:,jyy), &
          g_BSSN3_ll(:,:,:,jyz), g_BSSN3_ll(:,:,:,jzz), &
          A_BSSN3_ll(:,:,:,jxx), A_BSSN3_ll(:,:,:,jxy), &
          A_BSSN3_ll(:,:,:,jxz), A_BSSN3_ll(:,:,:,jyy), &
          A_BSSN3_ll(:,:,:,jyz), A_BSSN3_ll(:,:,:,jzz), &
          trK(:,:,:), phi(:,:,:), &
          Gamma_u(:,:,:,jx), &
          Gamma_u(:,:,:,jy), Gamma_u(:,:,:,jz), &
          Tmunu_ll(:,:,:,itt), &
          Tmunu_ll(:,:,:,itx), &
          Tmunu_ll(:,:,:,ity), &
          Tmunu_ll(:,:,:,itz), &
          Tmunu_ll(:,:,:,ixx), &
          Tmunu_ll(:,:,:,ixy), &
          Tmunu_ll(:,:,:,ixz), &
          Tmunu_ll(:,:,:,iyy), &
          Tmunu_ll(:,:,:,iyz), &
          Tmunu_ll(:,:,:,izz), &
          lapse(:,:,:), shift_u(:,:,:,jx), &
          shift_u(:,:,:,jy), shift_u(:,:,:,jz), &
          !
          !-- Output
          !
          ! Connection constraints
          GC_parts(:,:,:,jx), &
          GC_parts(:,:,:,jy), &
          GC_parts(:,:,:,jz), &
          ! Hamiltonian and momentum constraints
          HC_parts(:,:,:), &
          MC_parts(:,:,:,jx), &
          MC_parts(:,:,:,jy), &
          MC_parts(:,:,:,jz) &
        )
      END ASSOCIATE
    ENDDO
    PRINT *, " * Constraints computed."
    PRINT *

    IF( debug ) PRINT *, "0"

    !---------------------------------------------------------!
    !--  Analyze constraints, and print to formatted files  --!
    !---------------------------------------------------------!

    DO l= 1, THIS% nlevels, 1

      ASSOCIATE( HC_parts => THIS% HC_parts% levels(l)% var, &
                 MC_parts => THIS% MC_parts% levels(l)% var, &
                 GC_parts => THIS% GC_parts% levels(l)% var &
      )

        unit_logfile= 2791

        IF( l > 9 )THEN
          WRITE( n_reflev, "(I2)" ) l
        ELSE
          WRITE( n_reflev, "(I1)" ) l
        ENDIF

        finalname_logfile= TRIM(name_logfile)//"-reflev"//TRIM(n_reflev)//".log"

        INQUIRE( FILE= TRIM(finalname_logfile), EXIST= exist )

        IF( debug ) PRINT *, "1"

        IF( exist )THEN
            OPEN( UNIT= unit_logfile, FILE= TRIM(finalname_logfile), &
                  STATUS= "REPLACE", &
                  FORM= "FORMATTED", &
                  POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
                  IOMSG= err_msg )
        ELSE
            OPEN( UNIT= unit_logfile, FILE= TRIM(finalname_logfile), &
                  STATUS= "NEW", &
                  FORM= "FORMATTED", &
                  ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
        ENDIF
        IF( ios > 0 )THEN
          PRINT *, "...error when opening ", TRIM(finalname_logfile), &
                   ". The error message is", err_msg
          STOP
        ENDIF
        !CALL test_status( ios, err_msg, "...error when opening " &
        !                  // TRIM(name_logfile) )

        IF( debug ) PRINT *, "2"

        IF( .NOT.ALLOCATED( THIS% HC_parts_l2 ))THEN
          ALLOCATE( THIS% HC_parts_l2( THIS% nlevels ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array HC_parts_l2. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF
        IF( .NOT.ALLOCATED( THIS% MC_parts_l2 ))THEN
          ALLOCATE( THIS% MC_parts_l2( THIS% nlevels, 3 ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array MC_parts_l2. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF
        IF( .NOT.ALLOCATED( THIS% GC_parts_l2 ))THEN
          ALLOCATE( THIS% GC_parts_l2( THIS% nlevels, 3 ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array GC_parts_l2. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF
        IF( .NOT.ALLOCATED( THIS% HC_parts_loo ))THEN
          ALLOCATE( THIS% HC_parts_loo( THIS% nlevels ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array HC_parts_loo. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF
        IF( .NOT.ALLOCATED( THIS% MC_parts_loo ))THEN
          ALLOCATE( THIS% MC_parts_loo( THIS% nlevels, 3 ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array MC_parts_loo. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF
        IF( .NOT.ALLOCATED( THIS% GC_parts_loo ))THEN
          ALLOCATE( THIS% GC_parts_loo( THIS% nlevels, 3 ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array GC_parts_loo. ", &
                     "The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !                "...deallocation error for array HC" )
        ENDIF

        WRITE( UNIT = unit_logfile, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id

        IF( debug ) PRINT *, "3"

        PRINT *, "** Analyzing constraints on refinement level ", l, "..."

        name_analysis= "bssn-hc-parts-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the Hamiltonian constraint"
        CALL THIS% analyze_constraint( &
             l, &
             HC_parts, name_constraint, unit_logfile, name_analysis, &
             THIS% HC_parts_l2(l), THIS% HC_parts_loo(l) )

        IF( debug ) PRINT *, "3.1"

        name_analysis= "bssn-mc1-parts-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the first component of the momentum constraint"
        CALL THIS% analyze_constraint( &
             l, &
             MC_parts(:,:,:,jx), name_constraint, unit_logfile, name_analysis, &
             THIS% MC_parts_l2(l,jx), THIS% MC_parts_loo(l,jx) )

        name_analysis= "bssn-mc2-parts-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the second component of the momentum constraint"
        CALL THIS% analyze_constraint( &
             l, &
             MC_parts(:,:,:,jy), name_constraint, unit_logfile, name_analysis, &
             THIS% MC_parts_l2(l,jy), THIS% MC_parts_loo(l,jy) )

        name_analysis= "bssn-mc3-parts-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the third component of the momentum constraint"
        CALL THIS% analyze_constraint( &
             l, &
             MC_parts(:,:,:,jz), name_constraint, unit_logfile, name_analysis, &
             THIS% MC_parts_l2(l,jz), THIS% MC_parts_loo(l,jz) )

        name_analysis= "bssn-gc1-parts-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the first component of the connection constraint"
        CALL THIS% analyze_constraint( &
             l, &
             GC_parts(:,:,:,jx), name_constraint, unit_logfile, name_analysis, &
             THIS% GC_parts_l2(l,jx), THIS% GC_parts_loo(l,jx) )

        IF( debug ) PRINT *, "3.7"

        name_analysis= "bssn-gc2-parts-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the second component of the connection constraint"
        CALL THIS% analyze_constraint( &
             l, &
             GC_parts(:,:,:,jy), name_constraint, unit_logfile, name_analysis, &
             THIS% GC_parts_l2(l,jy), THIS% GC_parts_loo(l,jy) )

        name_analysis= "bssn-gc3-parts-analysis-reflev"//TRIM(n_reflev)//".dat"
        name_constraint= "the third component of the connection constraint"
        CALL THIS% analyze_constraint( &
             l, &
             GC_parts(:,:,:,jz), name_constraint, unit_logfile, name_analysis, &
             THIS% GC_parts_l2(l,jz), THIS% GC_parts_loo(l,jz) )

        IF( debug ) PRINT *, "4"

        CLOSE( UNIT= unit_logfile )

        PRINT *, " * Constraints analyzed. Summary of results saved to ", &
                 finalname_logfile
        PRINT *
      END ASSOCIATE
    ENDDO

    IF( THIS% export_constraints )THEN

      PRINT *, " * Printing constraints to file ", TRIM(namefile), "..."

      !
      !-- Export the constraints to a formatted file
      !

      INQUIRE( FILE= TRIM(namefile), EXIST= exist )

      IF( debug ) PRINT *, "1"

      IF( exist )THEN
          OPEN( UNIT= 21, FILE= TRIM(namefile), STATUS= "REPLACE", &
                FORM= "FORMATTED", &
                POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
                IOMSG= err_msg )
      ELSE
          OPEN( UNIT= 21, FILE= TRIM(namefile), STATUS= "NEW", &
          FORM= "FORMATTED", &
                ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
      ENDIF
      IF( ios > 0 )THEN
        PRINT *, "...error when opening ", TRIM(namefile), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when opening " &
      !         // TRIM(namefile) )

      IF( debug ) PRINT *, "2"

      WRITE( UNIT = 21, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id
      WRITE( UNIT = 21, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "# Values of the BSSN constraints computed with the mapping routines ", &
      "for the LORENE ID on selected grid points"
      IF( ios > 0 )THEN
        PRINT *, "...error when writing line 1 in ", TRIM(namefile), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when writing line 1 in "&
      !         // TRIM(namefile) )
      WRITE( UNIT = 21, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "# column:      1        2       3       4       5", &
      "       6       7       8       9       10", &
      "       11       12       13       14       15", &
      "       16       17       18       19       20"
      IF( ios > 0 )THEN
        PRINT *, "...error when writing line 2 in ", TRIM(namefile), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when writing line 2 in "&
      !        // TRIM(namefile) )
      WRITE( UNIT = 21, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "#      refinement level    x   y   z   Stress-energy (10 components)   "&
      // "Hamiltonian constraint       " &
      // "Momentum constraint (three components)       " &
      // "Connection constraint (three components)"
      IF( ios > 0 )THEN
        PRINT *, "...error when writing line 3 in ", TRIM(namefile), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when writing line 3 in "&
      !        // TRIM(namefile) )

      IF( debug ) PRINT *, "3"

      DO l= 1, THIS% nlevels, 1

        ASSOCIATE( lapse           => THIS% lapse% levels(l)% var, &
                   shift_u         => THIS% shift_u% levels(l)% var, &
                   phi             => THIS% phi% levels(l)% var, &
                   trK             => THIS% trK% levels(l)% var, &
                   g_BSSN3_ll      => THIS% g_BSSN3_ll% levels(l)% var, &
                   A_BSSN3_ll      => THIS% A_BSSN3_ll% levels(l)% var, &
                   g_phys3_ll      => THIS% g_phys3_ll% levels(l)% var, &
                   k_phys3_ll      => THIS% k_phys3_ll% levels(l)% var, &
                   Gamma_u         => THIS% Gamma_u% levels(l)% var, &
                   Tmunu_ll        => Tmunu_ll% levels(l)% var, &
                   HC_parts        => THIS% HC_parts% levels(l)% var, &
                   MC_parts        => THIS% MC_parts% levels(l)% var, &
                   GC_parts        => THIS% GC_parts% levels(l)% var &
        )

          ! Being abs_grid a local array, it is good practice to allocate it on
          ! the heap, otherwise it will be stored on the stack which has a very
          ! limited size. This results in a segmentation fault.
          IF( ALLOCATED( abs_grid ) )THEN
            DEALLOCATE( abs_grid )
          ENDIF
          ALLOCATE( abs_grid( THIS% get_ngrid_x(l), THIS% get_ngrid_y(l), &
                              THIS% get_ngrid_z(l), 3 ) )

          DO k= 1, THIS% get_ngrid_z(l), 1
            DO j= 1, THIS% get_ngrid_y(l), 1
              DO i= 1, THIS% get_ngrid_x(l), 1

                abs_grid( i, j, k, jx )= &
                            ABS( THIS% coords% levels(l)% var( i, j, k, jx ) )
                abs_grid( i, j, k, jy )= &
                            ABS( THIS% coords% levels(l)% var( i, j, k, jy ) )
                abs_grid( i, j, k, jz )= &
                            ABS( THIS% coords% levels(l)% var( i, j, k, jz ) )

              ENDDO
            ENDDO
          ENDDO

          min_abs_y= 1D+20
          min_abs_z= 1D+20
          DO k= 1, THIS% get_ngrid_z(l), 1
            DO j= 1, THIS% get_ngrid_y(l), 1
              DO i= 1, THIS% get_ngrid_x(l), 1

                IF( ABS( THIS% coords% levels(l)% var( i, j, k, jy ) ) &
                    < min_abs_y )THEN
                  min_abs_y= ABS( THIS% coords% levels(l)% var( i, j, k, jy ) )
                  min_ix_y= i
                  min_iy_y= j
                  min_iz_y= k
                ENDIF

                IF( ABS( THIS% coords% levels(l)% var( i, j, k, jz ) ) &
                    < min_abs_z )THEN
                  min_abs_z= ABS( THIS% coords% levels(l)% var( i, j, k, jz ) )
                  min_ix_z= i
                  min_iy_z= j
                  min_iz_z= k
                ENDIF

              ENDDO
            ENDDO
          ENDDO

          DO k= 1, THIS% get_ngrid_z(l), 1

            IF( MOD( k, THIS% cons_step ) /= 0 ) CYCLE

            DO j= 1, THIS% get_ngrid_y(l), 1

              IF( MOD( j, THIS% cons_step ) /= 0 ) CYCLE

              DO i= 1, THIS% get_ngrid_x(l), 1

                IF( MOD( i, THIS% cons_step ) /= 0 ) CYCLE

                IF( THIS% export_constraints_xy .AND. &
                    ( THIS% coords% levels(l)% var( i, j, k, jz ) /= &
                      THIS% coords% levels(l)% var( min_ix_z, min_iy_z, &
                                                    min_iz_z, jz ) ) )THEN
                  CYCLE
                ENDIF
                IF( THIS% export_constraints_x .AND. &
                    ( THIS% coords% levels(l)% var( i, j, k, jz ) /= &
                      THIS% coords% levels(l)% var( min_ix_z, min_iy_z, &
                                                    min_iz_z, jz ) &
                      .OR. &
                      THIS% coords% levels(l)% var( i, j, k, jy ) /= &
                      THIS% coords% levels(l)% var( min_ix_y, min_iy_y, &
                                                    min_iz_y, jy ) ) )THEN
                  CYCLE
                ENDIF

                IF( debug )THEN
                  WRITE( UNIT = 21, IOSTAT = ios, IOMSG = err_msg, FMT = * )&
                    l, &
                    THIS% coords% levels(l)% var( i, j, k, jx ), &
                    THIS% coords% levels(l)% var( i, j, k, jy ), &
                    THIS% coords% levels(l)% var( i, j, k, jz ), &
                    THIS% coords% levels(l)% var( i, j, k, jx ), &
                    THIS% coords% levels(l)% var( i, j, k, jy ), &
                    THIS% coords% levels(l)% var( i, j, k, jz ), &
                    THIS% coords% levels(l)% var( i, j, k, jx ), &
                    THIS% coords% levels(l)% var( i, j, k, jy ), &
                    THIS% coords% levels(l)% var( i, j, k, jz ), &
                    THIS% coords% levels(l)% var( i, j, k, jx ), &
                    THIS% coords% levels(l)% var( i, j, k, jy ), &
                    THIS% coords% levels(l)% var( i, j, k, jz ), &
                    THIS% coords% levels(l)% var( i, j, k, jx ), &
                    THIS% coords% levels(l)% var( i, j, k, jy ), &
                    THIS% coords% levels(l)% var( i, j, k, jz ), &
                    THIS% coords% levels(l)% var( i, j, k, jx ), &
                    THIS% coords% levels(l)% var( i, j, k, jy ), &
                    THIS% coords% levels(l)% var( i, j, k, jz ), & ! columns 18
                    !pos_loc( 1, ix, iy, iz ), &
                    !pos_loc( 2, ix, iy, iz ), &
                    !pos_loc( 3, ix, iy, iz ), &
                    !vel_loc( 1, ix, iy, iz ), &
                    !vel_loc( 2, ix, iy, iz ), &
                    !vel_loc( 3, ix, iy, iz ), &
                    !nu_loc( ix, iy, iz ), &
                    !u_loc( ix, iy, iz ), &
                    !nlrf_loc( ix, iy, iz ), &
                    !theta_loc( ix, iy, iz ), &
                    !pressure_loc( ix, iy, iz ), &
                    Tmunu_ll( i, j, k, itt ), &
                    Tmunu_ll( i, j, k, itx ), &
                    Tmunu_ll( i, j, k, ity ), &
                    Tmunu_ll( i, j, k, itz ), &
                    Tmunu_ll( i, j, k, ixx ), &
                    Tmunu_ll( i, j, k, ixy ), &
                    Tmunu_ll( i, j, k, ixz ), &
                    Tmunu_ll( i, j, k, iyy ), &
                    Tmunu_ll( i, j, k, iyz ), &
                    Tmunu_ll( i, j, k, izz ), &
                    HC_parts( i, j, k ), &
                    MC_parts( i, j, k, jx ), &
                    MC_parts( i, j, k, jy ), &
                    MC_parts( i, j, k, jz ), &
                    GC_parts( i, j, k, jx ), &
                    GC_parts( i, j, k, jy ), &
                    GC_parts( i, j, k, jz ), &
                    lapse( i, j, k ), &
                    shift_u( i, j, k, jx ), &
                    shift_u( i, j, k, jy ), &
                    shift_u( i, j, k, jz ), &
                    g_BSSN3_ll( i, j, k, jxx ), &
                    g_BSSN3_ll( i, j, k, jxy ), &
                    g_BSSN3_ll( i, j, k, jxz ), &
                    g_BSSN3_ll( i, j, k, jyy ), &
                    g_BSSN3_ll( i, j, k, jyz ), &
                    g_BSSN3_ll( i, j, k, jzz ), &
                    k_phys3_ll( i, j, k, jxx ), &
                    k_phys3_ll( i, j, k, jxy ), &
                    k_phys3_ll( i, j, k, jxz ), &
                    k_phys3_ll( i, j, k, jyy ), &
                    k_phys3_ll( i, j, k, jyz ), &
                    k_phys3_ll( i, j, k, jzz ), &
                    A_BSSN3_ll( i, j, k, jxx ), &
                    A_BSSN3_ll( i, j, k, jxy ), &
                    A_BSSN3_ll( i, j, k, jxz ), &
                    A_BSSN3_ll( i, j, k, jyy ), &
                    A_BSSN3_ll( i, j, k, jyz ), &
                    A_BSSN3_ll( i, j, k, jzz ), &
                    trK( i, j, k ), &
                    phi( i, j, k ), &
                    Gamma_u( i, j, k, 1 ), &
                    Gamma_u( i, j, k, 2 ), &
                    Gamma_u( i, j, k, 3 )
                ELSE
                  WRITE( UNIT = 21, IOSTAT = ios, IOMSG = err_msg, FMT = * )&
                    l, &
                    THIS% coords% levels(l)% var( i, j, k, jx ), &
                    THIS% coords% levels(l)% var( i, j, k, jy ), &
                    THIS% coords% levels(l)% var( i, j, k, jz ), &
                    Tmunu_ll( i, j, k, itt ), &
                    Tmunu_ll( i, j, k, itx ), &
                    Tmunu_ll( i, j, k, ity ), &
                    Tmunu_ll( i, j, k, itz ), &
                    Tmunu_ll( i, j, k, ixx ), &
                    Tmunu_ll( i, j, k, ixy ), &
                    Tmunu_ll( i, j, k, ixz ), &
                    Tmunu_ll( i, j, k, iyy ), &
                    Tmunu_ll( i, j, k, iyz ), &
                    Tmunu_ll( i, j, k, izz ), &
                    HC_parts( i, j, k ), &
                    MC_parts( i, j, k, jx ), &
                    MC_parts( i, j, k, jy ), &
                    MC_parts( i, j, k, jz ), &
                    GC_parts( i, j, k, jx ), &
                    GC_parts( i, j, k, jy ), &
                    GC_parts( i, j, k, jz )
                ENDIF

                IF( ios > 0 )THEN
                  PRINT *, "...error when writing the arrays in ", &
                           TRIM(namefile), &
                           ". The error message is", err_msg
                  STOP
                ENDIF
                !CALL test_status( ios, err_msg, &
                !                  "...error in writing " &
                !                  // "the arrays in " // TRIM(namefile) )
              ENDDO
            ENDDO
          ENDDO
        END ASSOCIATE
      ENDDO

      IF( debug ) PRINT *, "4"

      CLOSE( UNIT= 21 )

      PRINT *, " * Printed."
      PRINT *

      PRINT *, " * Printing sph density to file ", TRIM(namefile_sph), "..."

      INQUIRE( FILE= TRIM(namefile_sph), EXIST= exist )

      IF( exist )THEN
          OPEN( UNIT= 2, FILE= TRIM(namefile_sph), STATUS= "REPLACE", &
                FORM= "FORMATTED", &
                POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
                IOMSG= err_msg )
      ELSE
          OPEN( UNIT= 2, FILE= TRIM(namefile_sph), STATUS= "NEW", &
                FORM= "FORMATTED", &
                ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
      ENDIF
      IF( ios > 0 )THEN
        PRINT *, "...error when opening ", TRIM(namefile_sph), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when opening " &
      !                  // TRIM(namefile_sph) )

      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id

      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "# Values of the SPH density"
      IF( ios > 0 )THEN
        PRINT *, "...error when writing line 1 in ", TRIM(namefile_sph), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when writing line 1 in "&
      !        // TRIM(namefile_sph) )

      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "# column:      1        2       3       4"

      IF( ios > 0 )THEN
        PRINT *, "...error when writing line 2 in ", TRIM(namefile_sph), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when writing line 2 in "&
      !        // TRIM(namefile_sph) )

      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "#      particle      x [km]       y [km]       z [km]       ", &
      "SPH density"

      IF( ios > 0 )THEN
        PRINT *, "...error when writing line 3 in ", TRIM(namefile_sph), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when writing line 3 in "&
      !        // TRIM(namefile_sph) )

      !IF( ALLOCATED( abs_grid ) )THEN
      !  DEALLOCATE( abs_grid )
      !ENDIF
      !ALLOCATE( abs_grid( npart, 3 ) )

      !DO itr = 1, npart, 1
      !  abs_pos( itr, jx )= ABS( pos_loc( 1, itr ) )
      !  abs_pos( itr, jy )= ABS( pos_loc( 2, itr ) )
      !  abs_pos( itr, jz )= ABS( pos_loc( 3, itr ) )
      !ENDDO

      min_y_index= 0
      min_abs_y= 1D+20
      DO itr = 1, npart, 1
        IF( ABS( pos_loc( 2, itr ) ) < min_abs_y )THEN
          min_abs_y= ABS( pos_loc( 2, itr ) )
          min_y_index= itr
        ENDIF
      ENDDO

      min_abs_z= MINVAL( abs_pos( 3, : ) )

      write_data_loop: DO itr = 1, npart, 1

        IF( THIS% export_form_xy .AND. pos_loc( 3, itr ) /= min_abs_z )THEN
          CYCLE
        ENDIF
        IF( THIS% export_form_x .AND. ( pos_loc( 3, itr ) /= min_abs_z &
            .OR. pos_loc( 2, itr ) /= pos_loc( 2, min_y_index ) ) )THEN
          CYCLE
        ENDIF
        WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
          itr, &
          pos_loc( 1, itr ), &
          pos_loc( 2, itr ), &
          pos_loc( 3, itr ), &
          sph_density( itr )

        IF( ios > 0 )THEN
          PRINT *, "...error when writing the arrays in ", TRIM(namefile_sph), &
                   ". The error message is", err_msg
          STOP
        ENDIF
        !CALL test_status( ios, err_msg, "...error when writing " &
        !         // "the arrays in " // TRIM(namefile_sph) )
      ENDDO write_data_loop

      CLOSE( UNIT= 2 )

      PRINT *, " * Printed."
      PRINT *

    ENDIF

    !
    !-- Deallocate spacetime MODULE variables
    !
    CALL deallocate_ADM()
    CALL deallocate_Ztmp()
    CALL deallocate_Tmunu()
    CALL deallocate_GravityAcceleration()
    CALL deallocate_BSSN()
    !CALL deallocate_gravity_grid()

    ! Count the number of times that this SUBROUTINE is called, since the
    ! kernel has to be tabulated only once in the present implementation
    counter= counter+ 1

  END PROCEDURE compute_and_export_bssn_constraints_particles


  MODULE PROCEDURE deallocate_bssn_fields

    !**************************************************
    !                                                 *
    ! Deallocate BSSN memory                          *
    !                                                 *
    ! FT                                              *
    !                                                 *
    !**************************************************

    USE mesh_refinement, ONLY: deallocate_grid_function

    IMPLICIT NONE

    IF( ALLOCATED( THIS% Gamma_u% levels ) )THEN
      CALL deallocate_grid_function( THIS% Gamma_u, "Gamma_u_id" )
    ENDIF

    IF( ALLOCATED( THIS% phi% levels ) )THEN
      CALL deallocate_grid_function( THIS% phi, "phi_id" )
    ENDIF

    IF( ALLOCATED( THIS% trK% levels ) )THEN
      CALL deallocate_grid_function( THIS% trK, "trK_id" )
    ENDIF

    IF( ALLOCATED( THIS% A_BSSN3_ll% levels ) )THEN
      CALL deallocate_grid_function( THIS% A_BSSN3_ll, "A_BSSN3_ll_id" )
    ENDIF

    IF( ALLOCATED( THIS% g_BSSN3_ll% levels ) )THEN
      CALL deallocate_grid_function( THIS% g_BSSN3_ll, "g_BSSN3_ll_id" )
    ENDIF

    IF( ALLOCATED( THIS% GC% levels ) )THEN
      CALL deallocate_grid_function( THIS% GC, "GC_id" )
    ENDIF

  END PROCEDURE deallocate_bssn_fields


  !
  !-- Keeping the following two SUBROUTINES separate in case it is needed
  !-- to add other PROCEDURES to the destructor (probably superfluous...)
  !
  MODULE PROCEDURE destruct_bssn_id

    !**************************************************
    !                                                 *
    ! Finalizer for members of the extended class     *
    ! bssn_id, not the primitive class formul_3p1     *
    !                                                 *
    ! FT                                              *
    !                                                 *
    !**************************************************

    IMPLICIT NONE

    !PRINT *, "Inside destruct_BSSN_id"
    CALL deallocate_bssn_fields( THIS )

  END PROCEDURE destruct_bssn_id


  MODULE PROCEDURE destructor

    !**************************************************
    !                                                 *
    ! Destructor of TYPE bssn_id                      *
    !                                                 *
    ! FT                                              *
    !                                                 *
    !**************************************************

    IMPLICIT NONE

    !PRINT *, "Inside destructor"
    CALL destruct_bssn_id( THIS )
    CALL destruct_formul_3p1( THIS )

  END PROCEDURE destructor


END SUBMODULE bssn_id_methods
