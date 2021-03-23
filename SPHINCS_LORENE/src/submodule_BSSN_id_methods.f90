! File:         submodule_BSSN_id_methods.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (formul_bssn_id) bssn_id_methods

  !************************************************
  !                                               *
  ! Implementation of the methods of TYPE bssn_id *
  !                                               *
  ! FT 23.10.2020                                 *
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


    USE grav_grid,           ONLY: ngrid_x, ngrid_y, ngrid_z, &
                                   dx, dy, dz, dx_1, dy_1, dz_1, &
                                   xR, xL, yR, yL, zR, zL, &
                                   rad_coord, &
                                   deallocate_gravity_grid
    !USE NaNChecker,          ONLY: Check_Grid_Function_for_NAN
    USE tensor,              ONLY: itt, itx, ity, itz, ixx, ixy, &
                                   ixz, iyy, iyz, izz, jxx, jxy, jxz, &
                                   jyy, jyz, jzz, jx, jy, jz, n_sym3x3
    USE ADM,                 ONLY: lapse, dt_lapse, shift_u, dt_shift_u, &
                                   K_phys3_ll, g_phys3_ll, &
                                   allocate_ADM, deallocate_ADM
    USE McLachlan,           ONLY: initialize_BSSN, allocate_Ztmp, &
                                   deallocate_Ztmp, ADM_to_BSSN, &
                                   ADM_to_BSSN_args
    USE Tmunu,               ONLY: allocate_Tmunu, deallocate_Tmunu, Tmunu_ll
    USE GravityAcceleration, ONLY: dt_ehat_grav, dt_S_grav_l, &
                                   d_g_phys4_lll, &
                                   allocate_GravityAcceleration, &
                                   deallocate_GravityAcceleration
    USE options,             ONLY: basename
    USE constants,           ONLY: Msun_geo
    !
    !-- Use the arrays from the MODULE BSSN to store the BSSN variables
    !-- for the LORENE ID on the grid, and the SUBROUTINE write_BSSN_dump
    !-- to export them to the binary file needed by the evolution code
    !-- in SPHINCS
    !
    USE BSSN,   ONLY: allocate_BSSN, deallocate_BSSN, &
                      Gamma_u,          & ! Conformal connection
                      phi,              & ! Conformal factor
                      trK,              & ! Trace of extrinsic curvature
                      A_BSSN3_ll,       & ! Conformal traceless
                                          ! extrinsic curvature
                      g_BSSN3_ll,       & ! Conformal metric
                      Theta_Z4,         & ! Vector in the CCZ4 formulation.
                                          ! Loaded here because ADM_TO_BSSN
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

    INTEGER:: ix, iy, iz, i, j, k, allocation_status
    DOUBLE PRECISION:: detg


    PRINT *, "** Computing and exporting BSSN ID..."

    ! Allocate memory for the ADM MODULE variables (this has to be done since
    ! the MODULE SUBROUTINES need them; not allocating it results in a
    ! segmentation fault)
    PRINT *
    PRINT *, " * Allocating needed memory..."
    PRINT *

    ngrid_x= THIS% ngrid_x
    ngrid_y= THIS% ngrid_y
    ngrid_z= THIS% ngrid_z
    dx= THIS% dx
    dy= THIS% dy
    dz= THIS% dz
    dx_1= THIS% dx_1
    dy_1= THIS% dy_1
    dz_1= THIS% dz_1

    CALL allocate_ADM()
    CALL allocate_BSSN()

    ! Allocate temporary memory for time integration
    CALL allocate_Ztmp()

    ! Allocate memory for the stress-energy tensor (used in
    ! write_BSSN_dump)
    CALL allocate_Tmunu()

    ! Allocate memory for the derivatives of the ADM variables
    CALL allocate_GravityAcceleration()

    ! Set the stress-energy tensor to 0 (here dummy)
    Tmunu_ll(:,:,:,:)= 0.0D0

    !
    !-- Allocate memory for the coordinate radius
    !-- (this is needed by ADM_to_BSSN_args)
    !
    IF( .NOT. ALLOCATED( rad_coord ) )THEN
        ALLOCATE( rad_coord( THIS% ngrid_x, THIS% ngrid_y, &
                             THIS% ngrid_z ), STAT= allocation_status )
    ENDIF
    IF( allocation_status > 0 )THEN
       PRINT *, '...allocation error for rad_coord'
       STOP
    ENDIF

    DO k= 1, THIS% ngrid_z
      DO j= 1, THIS% ngrid_y
        DO i= 1, THIS% ngrid_x
          rad_coord( i, j, k )= SQRT( (xL+(i-1)*dx)**2 &
                                    + (yL+(j-1)*dy)**2 &
                                    + (zL+(k-1)*dz)**2 )
        ENDDO
      ENDDO
    ENDDO

    dt_lapse  = 0.0D0
    dt_shift_u= 0.0D0

    !
    !-- Compute BSSN variables, and time the process
    !-- The BSSN variables are stored in the MODULE variables since
    !-- write_BSSN_dump need them
    !
    PRINT *, " * Computing BSSN variables..."
    PRINT *
    CALL THIS% bssn_computer_timer% start_timer()
    CALL ADM_to_BSSN_args( &
      THIS% dx, THIS% dy, THIS% dz, &
      ! ADM variables (input)
      THIS% g_phys3_ll(:,:,:,jxx), THIS% g_phys3_ll(:,:,:,jxy), &
      THIS% g_phys3_ll(:,:,:,jxz), THIS% g_phys3_ll(:,:,:,jyy), &
      THIS% g_phys3_ll(:,:,:,jyz), THIS% g_phys3_ll(:,:,:,jzz), &
      THIS% K_phys3_ll(:,:,:,jxx), THIS% K_phys3_ll(:,:,:,jxy), &
      THIS% K_phys3_ll(:,:,:,jxz), THIS% K_phys3_ll(:,:,:,jyy), &
      THIS% K_phys3_ll(:,:,:,jyz), THIS% K_phys3_ll(:,:,:,jzz), &
      THIS% lapse(:,:,:), &
      THIS% shift_u(:,:,:,jx), &
      THIS% shift_u(:,:,:,jy), &
      THIS% shift_u(:,:,:,jz), &
      dt_lapse(:,:,:), &
      dt_shift_u(:,:,:,jx), dt_shift_u(:,:,:,jy), dt_shift_u(:,:,:,jz), &
      ! BSSN variables (output)
      g_BSSN3_ll(:,:,:,jxx), g_BSSN3_ll(:,:,:,jxy), &
      g_BSSN3_ll(:,:,:,jxz), g_BSSN3_ll(:,:,:,jyy), &
      g_BSSN3_ll(:,:,:,jyz), g_BSSN3_ll(:,:,:,jzz), &
      A_BSSN3_ll(:,:,:,jxx), A_BSSN3_ll(:,:,:,jxy), &
      A_BSSN3_ll(:,:,:,jxz), A_BSSN3_ll(:,:,:,jyy), &
      A_BSSN3_ll(:,:,:,jyz), A_BSSN3_ll(:,:,:,jzz), &
      phi(:,:,:), trK(:,:,:), Theta_Z4(:,:,:), &
      lapse_A_BSSN(:,:,:), &
      shift_B_BSSN_u(:,:,:,jx), shift_B_BSSN_u(:,:,:,jy), &
      shift_B_BSSN_u(:,:,:,jz), &
      Gamma_u(:,:,:,jx), Gamma_u(:,:,:,jy), &
      Gamma_u(:,:,:,jz) &
    )
    CALL THIS% bssn_computer_timer% stop_timer()

    ! Set the MODULE variables equal to the TYPE variables
    lapse= THIS% lapse
    shift_u= THIS% shift_u

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
    !-- Setting the local variables equal to the MODULE variables
    !
    THIS% g_BSSN3_ll= g_BSSN3_ll
    THIS% A_BSSN3_ll= A_BSSN3_ll
    THIS% phi       = phi
    THIS% trK       = trK
    THIS% Gamma_u   = Gamma_u

    ! Write BSSN ID to a binary file to be read by the evolution code
    ! in SPHINCS
    IF( THIS% export_bin )THEN
      IF( PRESENT(namefile) )THEN
        CALL write_BSSN_dump( namefile )
      ELSE
        CALL write_BSSN_dump()
      ENDIF
    ENDIF

    !
    !-- Deallocate MODULE variables
    !
    CALL deallocate_ADM()
    CALL deallocate_Ztmp()
    CALL deallocate_Tmunu()
    CALL deallocate_GravityAcceleration()
    CALL deallocate_BSSN()
    CALL deallocate_gravity_grid()

    call_flag= call_flag + 1
    THIS% call_flag= call_flag

    PRINT *, "** BSSN ID computed."
    PRINT *

  END PROCEDURE compute_and_export_bssn_variables


  MODULE PROCEDURE read_bssn_dump_print_formatted

    !************************************************
    !                                               *
    ! Read the BSSN ID from the binary file output  *
    ! by write_BSSN_dump, and print it to a         *
    ! formatted file                                *
    !                                               *
    ! FT 08.02.2021                                 *
    !                                               *
    !************************************************

    USE grav_grid,           ONLY: ngrid_x, ngrid_y, ngrid_z, &
                                   !dx, dy, dz, dx_1, dy_1, dz_1, &
                                   !xR, xL, yR, yL, zR, zL, &
                                   !rad_coord, &
                                   deallocate_gravity_grid
    USE tensor,              ONLY: jxx, jxy, jxz, jyy, jyz, jzz, jx, jy, jz
    USE ADM,                 ONLY: lapse, shift_u, &
                                   allocate_ADM, deallocate_ADM
    !USE McLachlan,           ONLY: initialize_BSSN, allocate_Ztmp, &
    !                               deallocate_Ztmp, ADM_to_BSSN, &
    !                               ADM_to_BSSN_args
    !USE Tmunu,               ONLY: allocate_Tmunu, deallocate_Tmunu
    !USE GravityAcceleration, ONLY: dt_ehat_grav, dt_S_grav_l, &
    !                               d_g_phys4_lll, &
    !                               allocate_GravityAcceleration, &
    !                               deallocate_GravityAcceleration
    USE BSSN,       ONLY: allocate_BSSN, deallocate_BSSN, &
                          Gamma_u,          & ! Conformal connection
                          phi,              & ! Conformal factor
                          trK,              & ! Trace of extrinsic curvature
                          A_BSSN3_ll,       & ! Conformal traceless
                                              ! extrinsic curvature
                          g_BSSN3_ll,       & ! Conformal metric
                          Theta_Z4,         & ! Vector in the CCZ4 formulation.
                                              ! Loaded here because ADM_TO_BSSN
                                              ! calls SUBROUTINES that need it
                                              ! as input; however, it is not
                                              ! evolved in BSSN
                          lapse_A_BSSN,     & ! Time derivative of lapse
                          shift_B_BSSN_u,   & ! Time derivativeof shift
                          read_BSSN_dump

    IMPLICIT NONE

    INTEGER:: ix, iy, iz, min_ix_y, min_iy_y, min_iz_y, &
              min_ix_z, min_iy_z, min_iz_z

    DOUBLE PRECISION:: min_abs_y, min_abs_z
    DOUBLE PRECISION, DIMENSION( :, :, :, : ), ALLOCATABLE:: abs_grid

    LOGICAL:: exist

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    PRINT *, "** Executing the read_bssn_dump_print_formatted subroutine..."

    ngrid_x= THIS% ngrid_x
    ngrid_y= THIS% ngrid_y
    ngrid_z= THIS% ngrid_z
    !dx= THIS% dx
    !dy= THIS% dy
    !dz= THIS% dz
    !dx_1= THIS% dx_1
    !dy_1= THIS% dy_1
    !dz_1= THIS% dz_1

    CALL allocate_ADM()
    CALL allocate_BSSN()

    ! Allocate temporary memory for time integration
    !CALL allocate_Ztmp()

    ! Allocate memory for the derivatives of the ADM variables
    !CALL allocate_GravityAcceleration()

    CALL read_BSSN_dump( ngrid_x, ngrid_y, ngrid_z, 00000, namefile_bin )

    ! Being abs_grid a local array, it is good practice to allocate it on the
    ! heap, otherwise it will be stored on the stack which has a very limited
    ! size. This results in a segmentation fault.
    ALLOCATE( abs_grid( 3, THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z ) )

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
      finalnamefile= "bssn_vars.dat"
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

    DO iz= 1, THIS% ngrid_z, 1
      DO iy= 1, THIS% ngrid_y, 1
        DO ix= 1, THIS% ngrid_x, 1
          abs_grid( 1, ix, iy, iz )= ABS( THIS% grid( 1, ix, iy, iz ) )
          abs_grid( 2, ix, iy, iz )= ABS( THIS% grid( 2, ix, iy, iz ) )
          abs_grid( 3, ix, iy, iz )= ABS( THIS% grid( 3, ix, iy, iz ) )
        ENDDO
      ENDDO
    ENDDO

    min_abs_y= 1D+20
    min_abs_z= 1D+20
    DO iz= 1, THIS% ngrid_z, 1
      DO iy= 1, THIS% ngrid_y, 1
        DO ix= 1, THIS% ngrid_x, 1
          IF( ABS( THIS% grid( 2, ix, iy, iz ) ) < min_abs_y )THEN
            min_abs_y= ABS( THIS% grid( 2, ix, iy, iz ) )
            min_ix_y= ix
            min_iy_y= iy
            min_iz_y= iz
          ENDIF
          IF( ABS( THIS% grid( 3, ix, iy, iz ) ) < min_abs_z )THEN
            min_abs_z= ABS( THIS% grid( 3, ix, iy, iz ) )
            min_ix_z= ix
            min_iy_z= iy
            min_iz_z= iz
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    coords_z: DO iz= 1, THIS% ngrid_z, 1
      coords_y: DO iy= 1, THIS% ngrid_y, 1
        coords_x: DO ix= 1, THIS% ngrid_x, 1

          IF( THIS% export_form_xy .AND. &
              THIS% grid( 3, ix, iy, iz ) /= &
              THIS% grid( 3, min_ix_z, min_iy_z, min_iz_z ) )THEN
            CYCLE
          ENDIF
          IF( THIS% export_form_x .AND. &
              ( THIS% grid( 3, ix, iy, iz ) /= &
                THIS% grid( 3, min_ix_z, min_iy_z, min_iz_z ) &
                .OR. &
                THIS% grid( 2, ix, iy, iz ) /= &
                THIS% grid( 2, min_ix_y, min_iy_y, min_iz_y ) ) )THEN
            CYCLE
          ENDIF

          WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * )&
              THIS% grid( 1, ix, iy, iz ), &
              THIS% grid( 2, ix, iy, iz ), &
              THIS% grid( 3, ix, iy, iz ), &
              lapse( ix, iy, iz ), &
              shift_u( ix, iy, iz, jx ), &
              shift_u( ix, iy, iz, jy ), &
              shift_u( ix, iy, iz, jz ), &
              phi(ix,iy,iz), &
              trK(ix,iy,iz), &
              g_BSSN3_ll( ix, iy, iz, jxx ), &
              g_BSSN3_ll( ix, iy, iz, jxy ), &
              g_BSSN3_ll( ix, iy, iz, jxz ), &
              g_BSSN3_ll( ix, iy, iz, jyy ), &
              g_BSSN3_ll( ix, iy, iz, jyz ), &
              g_BSSN3_ll( ix, iy, iz, jzz ), &
              A_BSSN3_ll( ix, iy, iz, jxx ), &
              A_BSSN3_ll( ix, iy, iz, jxy ), &
              A_BSSN3_ll( ix, iy, iz, jxz ), &
              A_BSSN3_ll( ix, iy, iz, jyy ), &
              A_BSSN3_ll( ix, iy, iz, jyz ), &
              A_BSSN3_ll( ix, iy, iz, jzz ), &
              Gamma_u( ix, iy, iz, jx ), &
              Gamma_u( ix, iy, iz, jy ), &
              Gamma_u( ix, iy, iz, jz )

          IF( ios > 0 )THEN
            PRINT *, "...error when writing the arrays in ", &
                     TRIM(finalnamefile), ". The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, "...error when writing " &
          !                  // "the arrays in " // TRIM(namefile) )

        ENDDO coords_x
      ENDDO coords_y
    ENDDO coords_z

    CLOSE( UNIT= 20 )

    !
    !-- Deallocate MODULE variables
    !
    CALL deallocate_ADM()
    !CALL deallocate_Ztmp()
    !CALL deallocate_GravityAcceleration()
    CALL deallocate_BSSN()
    !CALL deallocate_gravity_grid()

    PRINT *, " * LORENE BSSN ID on the gravity grid saved to formatted " &
             // "file ", TRIM(namefile)

    PRINT *, "** Subroutine read_bssn_dump_print_formatted " &
             // "executed."
    PRINT *

  END PROCEDURE read_bssn_dump_print_formatted


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

    IMPLICIT NONE

    INTEGER:: ix, iy, iz, min_ix_y, min_iy_y, min_iz_y, &
              min_ix_z, min_iy_z, min_iz_z

    DOUBLE PRECISION:: min_abs_y, min_abs_z
    DOUBLE PRECISION, DIMENSION( :, :, :, : ), ALLOCATABLE:: abs_grid

    LOGICAL:: exist

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    ! Being abs_grid a local array, it is good practice to allocate it on the
    ! heap, otherwise it will be stored on the stack which has a very limited
    ! size. This results in a segmentation fault.
    ALLOCATE( abs_grid( 3, THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z ) )

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

    DO iz= 1, THIS% ngrid_z, 1
      DO iy= 1, THIS% ngrid_y, 1
        DO ix= 1, THIS% ngrid_x, 1
          abs_grid( 1, ix, iy, iz )= ABS( THIS% grid( 1, ix, iy, iz ) )
          abs_grid( 2, ix, iy, iz )= ABS( THIS% grid( 2, ix, iy, iz ) )
          abs_grid( 3, ix, iy, iz )= ABS( THIS% grid( 3, ix, iy, iz ) )
        ENDDO
      ENDDO
    ENDDO

    min_abs_y= 1D+20
    min_abs_z= 1D+20
    DO iz= 1, THIS% ngrid_z, 1
      DO iy= 1, THIS% ngrid_y, 1
        DO ix= 1, THIS% ngrid_x, 1
          IF( ABS( THIS% grid( 2, ix, iy, iz ) ) < min_abs_y )THEN
            min_abs_y= ABS( THIS% grid( 2, ix, iy, iz ) )
            min_ix_y= ix
            min_iy_y= iy
            min_iz_y= iz
          ENDIF
          IF( ABS( THIS% grid( 3, ix, iy, iz ) ) < min_abs_z )THEN
            min_abs_z= ABS( THIS% grid( 3, ix, iy, iz ) )
            min_ix_z= ix
            min_iy_z= iy
            min_iz_z= iz
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    coords_z: DO iz= 1, THIS% ngrid_z, 1
      coords_y: DO iy= 1, THIS% ngrid_y, 1
        coords_x: DO ix= 1, THIS% ngrid_x, 1

          IF( THIS% export_form_xy .AND. &
              THIS% grid( 3, ix, iy, iz ) /= &
              THIS% grid( 3, min_ix_z, min_iy_z, min_iz_z ) )THEN
            CYCLE
          ENDIF
          IF( THIS% export_form_x .AND. &
              ( THIS% grid( 3, ix, iy, iz ) /= &
                THIS% grid( 3, min_ix_z, min_iy_z, min_iz_z ) &
                .OR. &
                THIS% grid( 2, ix, iy, iz ) /= &
                THIS% grid( 2, min_ix_y, min_iy_y, min_iz_y ) ) )THEN
            CYCLE
          ENDIF

          WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * )&
              THIS% grid( 1, ix, iy, iz ), &
              THIS% grid( 2, ix, iy, iz ), &
              THIS% grid( 3, ix, iy, iz ), &
              THIS% lapse( ix, iy, iz ), &
              THIS% shift_u( ix, iy, iz, jx ), &
              THIS% shift_u( ix, iy, iz, jy ), &
              THIS% shift_u( ix, iy, iz, jz ), &
              THIS% phi(ix,iy,iz), &
              THIS% trK(ix,iy,iz), &
              THIS% g_BSSN3_ll( ix, iy, iz, jxx ), &
              THIS% g_BSSN3_ll( ix, iy, iz, jxy ), &
              THIS% g_BSSN3_ll( ix, iy, iz, jxz ), &
              THIS% g_BSSN3_ll( ix, iy, iz, jyy ), &
              THIS% g_BSSN3_ll( ix, iy, iz, jyz ), &
              THIS% g_BSSN3_ll( ix, iy, iz, jzz ), &
              THIS% A_BSSN3_ll( ix, iy, iz, jxx ), &
              THIS% A_BSSN3_ll( ix, iy, iz, jxy ), &
              THIS% A_BSSN3_ll( ix, iy, iz, jxz ), &
              THIS% A_BSSN3_ll( ix, iy, iz, jyy ), &
              THIS% A_BSSN3_ll( ix, iy, iz, jyz ), &
              THIS% A_BSSN3_ll( ix, iy, iz, jzz ), &
              THIS% Gamma_u( ix, iy, iz, jx ), &
              THIS% Gamma_u( ix, iy, iz, jy ), &
              THIS% Gamma_u( ix, iy, iz, jz )

          IF( ios > 0 )THEN
            PRINT *, "...error when writing the arrays in ", &
                     TRIM(finalnamefile), ". The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, "...error when writing " &
          !                  // "the arrays in " // TRIM(finalnamefile) )

        ENDDO coords_x
      ENDDO coords_y
    ENDDO coords_z

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

    USE constants,  ONLY: c_light2, cm2m, MSun, g2kg, m2cm, &
                          lorene2hydrobase, MSun_geo, pi
    USE matrix,     ONLY: invert_4x4_matrix
    USE grav_grid,  ONLY: ngrid_x, ngrid_y, ngrid_z, dx, dy, dz, &
                          dx_1, dy_1, dz_1
    USE tensor,     ONLY: itt, itx, ity, itz, ixx, ixy, &
                          ixz, iyy, iyz, izz, jxx, jxy, jxz, &
                          jyy, jyz, jzz, jx, jy, jz, n_sym3x3, n_sym4x4
    USE McLachlan,  ONLY: ghost_size, BSSN_CONSTRAINTS_INTERIOR

    IMPLICIT NONE

    INTEGER:: ix, iy, iz, allocation_status, fd_lim
    INTEGER, DIMENSION(3) :: imin, imax
    INTEGER:: unit_logfile, min_ix_y, min_iy_y, min_iz_y, &
              min_ix_z, min_iy_z, min_iz_z

    DOUBLE PRECISION:: min_abs_y, min_abs_z
    DOUBLE PRECISION, DIMENSION( :, :, :, : ), ALLOCATABLE:: abs_grid

    DOUBLE PRECISION, DIMENSION( THIS% ngrid_x, &
                                 THIS% ngrid_y, &
                                 THIS% ngrid_z ):: baryon_density
    DOUBLE PRECISION, DIMENSION( THIS% ngrid_x, &
                                 THIS% ngrid_y, &
                                 THIS% ngrid_z ):: energy_density
    DOUBLE PRECISION, DIMENSION( THIS% ngrid_x, &
                                 THIS% ngrid_y, &
                                 THIS% ngrid_z ):: specific_energy
    DOUBLE PRECISION, DIMENSION( THIS% ngrid_x, &
                                 THIS% ngrid_y, &
                                 THIS% ngrid_z ):: pressure
    DOUBLE PRECISION, DIMENSION( THIS% ngrid_x, &
                                 THIS% ngrid_y, &
                                 THIS% ngrid_z, 3 ):: v_euler
    DOUBLE PRECISION, DIMENSION( THIS% ngrid_x, &
                                 THIS% ngrid_y, &
                                 THIS% ngrid_z, 3 ):: v_euler_l
    DOUBLE PRECISION, DIMENSION( THIS% ngrid_x, &
                                 THIS% ngrid_y, &
                                 THIS% ngrid_z, 0:3 ):: u_euler_l
    DOUBLE PRECISION, DIMENSION( THIS% ngrid_x, &
                                 THIS% ngrid_y, &
                                 THIS% ngrid_z ):: lorentz_factor
    DOUBLE PRECISION:: u_euler_norm= 0.0D0
    DOUBLE PRECISION:: detg4
    ! Spacetime metric
    DOUBLE PRECISION, DIMENSION( THIS% ngrid_x, &
                                 THIS% ngrid_y, &
                                 THIS% ngrid_z, &
                                 n_sym4x4):: g4
    ! Stress-energy tensor
    DOUBLE PRECISION, DIMENSION( THIS% ngrid_x, &
                                 THIS% ngrid_y, &
                                 THIS% ngrid_z, &
                                 n_sym4x4):: Tmunu_ll
    ! Spacetime metric as a 4x4 matrix
    DOUBLE PRECISION, DIMENSION( 4, 4 ):: g4temp
    ! Inverse spacetime metric as a 4x4 matrix
    DOUBLE PRECISION, DIMENSION( 4, 4 ):: ig4

    ! Declaration of debug variables needed to compute the Hamiltonian
    ! constraint directly, without calling the Cactus-bound SUBROUTINE
    ! BSSN_CONSTRAINTS_INTERIOR
    DOUBLE PRECISION, DIMENSION( THIS% ngrid_x, &
                                 THIS% ngrid_y, &
                                 THIS% ngrid_z ):: HC_hand
    DOUBLE PRECISION, DIMENSION( THIS% ngrid_x, &
                                 THIS% ngrid_y, &
                                 THIS% ngrid_z ):: HC_rho
    DOUBLE PRECISION, DIMENSION( THIS% ngrid_x, &
                                 THIS% ngrid_y, &
                                 THIS% ngrid_z ):: HC_trK
    DOUBLE PRECISION, DIMENSION( THIS% ngrid_x, &
                                 THIS% ngrid_y, &
                                 THIS% ngrid_z ):: HC_A
    DOUBLE PRECISION, DIMENSION( THIS% ngrid_x, &
                                 THIS% ngrid_y, &
                                 THIS% ngrid_z ):: HC_derphi

    CHARACTER( LEN= : ), ALLOCATABLE:: name_constraint
    CHARACTER( LEN= : ), ALLOCATABLE:: name_analysis

    LOGICAL:: exist
    LOGICAL, PARAMETER:: debug= .FALSE.

    ! Being abs_grid a local array, it is good practice to allocate it on the
    ! heap, otherwise it will be stored on the stack which has a very limited
    ! size. This results in a segmentation fault.
    ALLOCATE( abs_grid( 3, THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z ) )

    ngrid_x= THIS% ngrid_x
    ngrid_y= THIS% ngrid_y
    ngrid_z= THIS% ngrid_z
    dx= THIS% dx
    dy= THIS% dy
    dz= THIS% dz
    dx_1= THIS% dx_1
    dy_1= THIS% dy_1
    dz_1= THIS% dz_1

    !
    !-- Keeping the following lines commented, in case the arrays have to be
    !-- allocated on the heap with ALLOCATE
    !

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

    ALLOCATE( THIS% HC( THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z ) )
    ALLOCATE( THIS% MC( THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, 3 ) )
    ALLOCATE( THIS% GC( THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, 3 ) )
    !ALLOCATE( Tmunu_ll( &
    !             THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, n_sym4x4 ) )

    !
    !-- Import the hydro on the grid from the LORENE ID directly
    !

    !ALLOCATE( baryon_density( &
    !            THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z ) )
    !ALLOCATE( energy_density( &
    !            THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z ) )
    !ALLOCATE( specific_energy( &
    !            THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z ) )
    !ALLOCATE( pressure( &
    !            THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z ) )
    !ALLOCATE( v_euler( &
    !            THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, 3 ) )
    !ALLOCATE( v_euler_l( &
    !            THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, 3 ) )
    !ALLOCATE( lorentz_factor( &
    !            THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z ) )
    !!ALLOCATE( u_euler0( &
    !!            THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z ) )
    !!ALLOCATE( v_coord( &
    !!            THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, 0:3 ) )
    !!ALLOCATE( u_coord( &
    !!            THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, 0:3 ) )
    !!ALLOCATE( u_coord_l( &
    !!            THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, 0:3 ) )
    !ALLOCATE( u_euler_l( &
    !            THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, 0:3 ) )
    !ALLOCATE( g4( &
    !            THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, n_sym4x4 ) )
    !ALLOCATE( g4temp( 4, 4 ) )
    !ALLOCATE( ig4( 4, 4 ) )

    !
    !-- Import the hydro LORENE ID on the gravity grid
    !
    PRINT *, "** Importing LORENE hydro ID on the gravity grid..."
    CALL bns_obj% import_id( THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, &
                             THIS% grid, &
                             baryon_density, &
                             energy_density, &
                             specific_energy, &
                             pressure, &
                             v_euler )
    PRINT *, " * LORENE hydro ID imported."
    PRINT *

    !
    !-- Replace the points with negative hydro fields near the surface
    !-- with vacuum
    !
    PRINT *, "** Cleaning LORENE hydro ID around the surfaces of the stars..."
    DO iz= 1, THIS% ngrid_z, 1
      DO iy= 1, THIS% ngrid_y, 1
        DO ix= 1, THIS% ngrid_x, 1

          IF(      baryon_density ( ix, iy, iz ) < 0.0D0 &
              .OR. energy_density ( ix, iy, iz ) < 0.0D0 &
              .OR. specific_energy( ix, iy, iz ) < 0.0D0 &
              .OR. pressure       ( ix, iy, iz ) < 0.0D0 )THEN
              baryon_density ( ix, iy, iz )= 0.0D0
              energy_density ( ix, iy, iz )= 0.0D0
              specific_energy( ix, iy, iz )= 0.0D0
              pressure       ( ix, iy, iz )= 0.0D0
              v_euler        ( ix, iy, iz, : )= 0.0D0
          ENDIF

        ! Print progress on screen
        perc= 100*(THIS% ngrid_x*THIS% ngrid_y*(iz - 1) &
              + THIS% ngrid_x*(iy - 1) + ix) &
              /( THIS% ngrid_x* THIS% ngrid_y*THIS% ngrid_z )
        IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
          WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
                  creturn//" ", perc, "%"
        ENDIF

        ENDDO
      ENDDO
    ENDDO
    WRITE( *, "(A1)", ADVANCE= "NO" ) creturn
    PRINT *, " * LORENE hydro ID cleaned."
    PRINT *

    !THIS% shift_u= 0.0D0*2.15D-1*THIS% shift_u

    !---------------------------!
    !--  Compute constraints  --!
    !---------------------------!

    !
    !-- Compute the fluid 4-velocity in the coordinate frame
    !
    PRINT *, "** Computing fluid 4-velocity wrt Eulerian observer..."
    DO iz= 1, THIS% ngrid_z, 1
      DO iy= 1, THIS% ngrid_y, 1
        DO ix= 1, THIS% ngrid_x, 1

          !energy_density( ix, iy, iz )= baryon_density( ix, iy, iz ) &
          !                            + ( specific_energy(ix,iy,iz) + 1.0 ) &
          !                                 *baryon_density( ix, iy, iz )

          v_euler_l(ix,iy,iz,1)= THIS% g_phys3_ll(ix,iy,iz,jxx)&
                                      *v_euler(ix,iy,iz,1)&
                               + THIS% g_phys3_ll(ix,iy,iz,jxy)&
                                      *v_euler(ix,iy,iz,2)&
                               + THIS% g_phys3_ll(ix,iy,iz,jxz)&
                                      *v_euler(ix,iy,iz,3)
          v_euler_l(ix,iy,iz,2)= THIS% g_phys3_ll(ix,iy,iz,jxy)&
                                      *v_euler(ix,iy,iz,1)&
                               + THIS% g_phys3_ll(ix,iy,iz,jyy)&
                                      *v_euler(ix,iy,iz,2)&
                               + THIS% g_phys3_ll(ix,iy,iz,jyz)&
                                      *v_euler(ix,iy,iz,3)
          v_euler_l(ix,iy,iz,3)= THIS% g_phys3_ll(ix,iy,iz,jxz)&
                                      *v_euler(ix,iy,iz,1)&
                               + THIS% g_phys3_ll(ix,iy,iz,jyz)&
                                      *v_euler(ix,iy,iz,2)&
                               + THIS% g_phys3_ll(ix,iy,iz,jzz)&
                                      *v_euler(ix,iy,iz,3)

          lorentz_factor( ix, iy, iz )= 1.0D0/SQRT( 1.0D0 &
                              - ( v_euler_l(ix,iy,iz,1)*v_euler(ix,iy,iz,1) &
                                + v_euler_l(ix,iy,iz,2)*v_euler(ix,iy,iz,2) &
                                + v_euler_l(ix,iy,iz,3)*v_euler(ix,iy,iz,3) ) )

          u_euler_l(ix,iy,iz,0)= lorentz_factor( ix, iy, iz ) &
             *( - THIS% lapse( ix, iy, iz ) &
                + v_euler_l( ix, iy, iz, 1 )*THIS% shift_u( ix, iy, iz, 1 ) &
                + v_euler_l( ix, iy, iz, 2 )*THIS% shift_u( ix, iy, iz, 2 ) &
                + v_euler_l( ix, iy, iz, 3 )*THIS% shift_u( ix, iy, iz, 3 ) )
          u_euler_l(ix,iy,iz,1)= lorentz_factor( ix, iy, iz ) &
                                 *v_euler_l( ix, iy, iz, 1 )
          u_euler_l(ix,iy,iz,2)= lorentz_factor( ix, iy, iz ) &
                                 *v_euler_l( ix, iy, iz, 2 )
          u_euler_l(ix,iy,iz,3)= lorentz_factor( ix, iy, iz ) &
                                 *v_euler_l( ix, iy, iz, 3 )

          CALL compute_g4( ix, iy, iz, THIS% lapse, THIS% shift_u, &
                           THIS% g_phys3_ll, g4 )

          CALL determinant_sym4x4_grid( ix, iy, iz, g4, detg4 )

          IF( ABS( detg4 ) < 1.0D-10 )THEN
              PRINT *, "The determinant of the spacetime metric "&
                       // "is effectively 0 at the grid point " &
                       // "(ix,iy,iz)= (", ix, ",", iy,",",iz, &
                          ")."
              PRINT *, "detg4=", detg4
              PRINT *
              STOP
          ELSEIF( detg4 > 0.0D0 )THEN
              PRINT *, "The determinant of the spacetime metric "&
                       // "is positive at the grid point " &
                       // "(ix,iy,iz)= (", ix, ",", iy,",",iz, &
                          ")."
              PRINT *, "detg4=", detg4
              PRINT *
              STOP
          ENDIF

          g4temp(1,1)= g4(ix,iy,iz,itt)
          g4temp(1,2)= g4(ix,iy,iz,itx)
          g4temp(1,3)= g4(ix,iy,iz,ity)
          g4temp(1,4)= g4(ix,iy,iz,itz)

          g4temp(2,1)= g4(ix,iy,iz,itx)
          g4temp(2,2)= g4(ix,iy,iz,ixx)
          g4temp(2,3)= g4(ix,iy,iz,ixy)
          g4temp(2,4)= g4(ix,iy,iz,ixz)

          g4temp(3,1)= g4(ix,iy,iz,ity)
          g4temp(3,2)= g4(ix,iy,iz,ixy)
          g4temp(3,3)= g4(ix,iy,iz,iyy)
          g4temp(3,4)= g4(ix,iy,iz,iyz)

          g4temp(4,1)= g4(ix,iy,iz,itz)
          g4temp(4,2)= g4(ix,iy,iz,ixz)
          g4temp(4,3)= g4(ix,iy,iz,iyz)
          g4temp(4,4)= g4(ix,iy,iz,izz)

          CALL invert_4x4_matrix( g4temp, ig4 )

          u_euler_norm= ig4(1,1)* &
                        u_euler_l(ix,iy,iz,0)*u_euler_l(ix,iy,iz,0) &
                      + 2.0D0*ig4(1,2)* &
                        u_euler_l(ix,iy,iz,0)*u_euler_l(ix,iy,iz,1) &
                      + 2.0D0*ig4(1,3)* &
                        u_euler_l(ix,iy,iz,0)*u_euler_l(ix,iy,iz,2) &
                      + 2.0D0*ig4(1,4)* &
                        u_euler_l(ix,iy,iz,0)*u_euler_l(ix,iy,iz,3) &
                      + ig4(2,2)* &
                        u_euler_l(ix,iy,iz,1)*u_euler_l(ix,iy,iz,1) &
                      + 2.0D0*ig4(2,3)* &
                        u_euler_l(ix,iy,iz,1)*u_euler_l(ix,iy,iz,2) &
                      + 2.0D0*ig4(2,4)* &
                        u_euler_l(ix,iy,iz,1)*u_euler_l(ix,iy,iz,3) &
                      + ig4(3,3)* &
                        u_euler_l(ix,iy,iz,2)*u_euler_l(ix,iy,iz,2) &
                      + 2.0D0*ig4(3,4)* &
                        u_euler_l(ix,iy,iz,2)*u_euler_l(ix,iy,iz,3) &
                      + 2.0D0*ig4(4,4)* &
                        u_euler_l(ix,iy,iz,3)*u_euler_l(ix,iy,iz,3)

          IF( ABS( u_euler_norm + 1.0D0 ) > 1.0D-4 )THEN
              PRINT *, "** ERROR! The fluid 4-velocity in the " &
                       // "coordinate frame does not have norm -1. " &
                       // "The norm is", u_euler_norm
              STOP
          ENDIF

        ! Print progress on screen
        perc= 100*(THIS% ngrid_x*THIS% ngrid_y*(iz - 1) &
              + THIS% ngrid_x*(iy - 1) + ix) &
              /( THIS% ngrid_x*THIS% ngrid_y*THIS% ngrid_z )
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
    DO iz= 1, THIS% ngrid_z, 1
      DO iy= 1, THIS% ngrid_y, 1
        DO ix= 1, THIS% ngrid_x, 1

          Tmunu_ll(ix,iy,iz,itt)= lorene2hydrobase*( &
                  ( energy_density(ix,iy,iz) + pressure(ix,iy,iz) ) &
                  *u_euler_l(ix,iy,iz,0)*u_euler_l(ix,iy,iz,0) &
                  + pressure(ix,iy,iz)*g4(ix,iy,iz,itt) &
                   )

          Tmunu_ll(ix,iy,iz,itx)= lorene2hydrobase*( &
                  ( energy_density(ix,iy,iz) + pressure(ix,iy,iz) ) &
                  *u_euler_l(ix,iy,iz,0)*u_euler_l(ix,iy,iz,1) &
                  + pressure(ix,iy,iz)*g4(ix,iy,iz,itx) &
                   )

          Tmunu_ll(ix,iy,iz,ity)= lorene2hydrobase*( &
                  ( energy_density(ix,iy,iz) + pressure(ix,iy,iz) ) &
                  *u_euler_l(ix,iy,iz,0)*u_euler_l(ix,iy,iz,2) &
                  + pressure(ix,iy,iz)*g4(ix,iy,iz,ity) &
                   )

          Tmunu_ll(ix,iy,iz,itz)= lorene2hydrobase*( &
                  ( energy_density(ix,iy,iz) + pressure(ix,iy,iz) ) &
                  *u_euler_l(ix,iy,iz,0)*u_euler_l(ix,iy,iz,3) &
                  + pressure(ix,iy,iz)*g4(ix,iy,iz,itz) &
                   )

          Tmunu_ll(ix,iy,iz,ixx)= lorene2hydrobase*( &
                  ( energy_density(ix,iy,iz) + pressure(ix,iy,iz) ) &
                  *u_euler_l(ix,iy,iz,1)*u_euler_l(ix,iy,iz,1) &
                  + pressure(ix,iy,iz)*g4(ix,iy,iz,ixx) &
                   )

          Tmunu_ll(ix,iy,iz,ixy)= lorene2hydrobase*( &
                  ( energy_density(ix,iy,iz) + pressure(ix,iy,iz) ) &
                  *u_euler_l(ix,iy,iz,1)*u_euler_l(ix,iy,iz,2) &
                  + pressure(ix,iy,iz)*g4(ix,iy,iz,ixy) &
                   )

          Tmunu_ll(ix,iy,iz,ixz)= lorene2hydrobase*( &
                  ( energy_density(ix,iy,iz) + pressure(ix,iy,iz) ) &
                  *u_euler_l(ix,iy,iz,1)*u_euler_l(ix,iy,iz,3) &
                  + pressure(ix,iy,iz)*g4(ix,iy,iz,ixz) &
                   )

          Tmunu_ll(ix,iy,iz,iyy)= lorene2hydrobase*( &
                  ( energy_density(ix,iy,iz) + pressure(ix,iy,iz) ) &
                  *u_euler_l(ix,iy,iz,2)*u_euler_l(ix,iy,iz,2) &
                  + pressure(ix,iy,iz)*g4(ix,iy,iz,iyy)  &
                   )

          Tmunu_ll(ix,iy,iz,iyz)= lorene2hydrobase*( &
                  ( energy_density(ix,iy,iz) + pressure(ix,iy,iz) ) &
                  *u_euler_l(ix,iy,iz,2)*u_euler_l(ix,iy,iz,3) &
                  + pressure(ix,iy,iz)*g4(ix,iy,iz,iyz) &
                   )

          Tmunu_ll(ix,iy,iz,izz)= lorene2hydrobase*( &
                  ( energy_density(ix,iy,iz) + pressure(ix,iy,iz) ) &
                  *u_euler_l(ix,iy,iz,3)*u_euler_l(ix,iy,iz,3) &
                  + pressure(ix,iy,iz)*g4(ix,iy,iz,izz) &
                   )

        ! Print progress on screen
        perc= 100*(THIS% ngrid_x*THIS% ngrid_y*(iz - 1) &
              + THIS% ngrid_x*(iy - 1) + ix) &
              /( THIS% ngrid_x* THIS% ngrid_y*THIS% ngrid_z )
        IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
          WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
                  creturn//" ", perc, "%"
        ENDIF

        ENDDO
      ENDDO
    ENDDO
    WRITE( *, "(A1)", ADVANCE= "NO" ) creturn
    PRINT *, " * Stress-energy tensor computed."
    PRINT *

    ! In debug mode, compute the Hamiltonian constraint by hand
    IF( debug )THEN

      HC_rho= 0.0D0
      HC_trK= 0.0D0
      HC_A= 0.0D0
      HC_derphi= 0.0D0
      HC_hand= 0.0D0
      fd_lim= 5
      DO iz= fd_lim, THIS% ngrid_z - fd_lim, 1
        DO iy= fd_lim, THIS% ngrid_y - fd_lim, 1
          DO ix= fd_lim, THIS% ngrid_x - fd_lim, 1

            HC_rho( ix, iy, iz )= 2.0D0*pi*EXP(5.0D0*THIS% phi( ix, iy, iz )) &
                                  *lorene2hydrobase*energy_density( ix, iy, iz )

            HC_trK( ix, iy, iz )= - EXP(5.0D0*THIS% phi( ix, iy, iz ))/12.0D0 &
                                    *THIS% trK( ix, iy, iz )**2

            HC_A( ix, iy, iz )= EXP(5.0D0*THIS% phi( ix, iy, iz ))/8.0D0 &
             *( THIS% A_BSSN3_ll(ix,iy,iz,jxx)*THIS% A_BSSN3_ll(ix,iy,iz,jxx) &
              + THIS% A_BSSN3_ll(ix,iy,iz,jxy)*THIS% A_BSSN3_ll(ix,iy,iz,jxy) &
              + THIS% A_BSSN3_ll(ix,iy,iz,jxz)*THIS% A_BSSN3_ll(ix,iy,iz,jxz) &
              + THIS% A_BSSN3_ll(ix,iy,iz,jxy)*THIS% A_BSSN3_ll(ix,iy,iz,jxy) &
              + THIS% A_BSSN3_ll(ix,iy,iz,jyy)*THIS% A_BSSN3_ll(ix,iy,iz,jyy) &
              + THIS% A_BSSN3_ll(ix,iy,iz,jyz)*THIS% A_BSSN3_ll(ix,iy,iz,jyz) &
              + THIS% A_BSSN3_ll(ix,iy,iz,jxz)*THIS% A_BSSN3_ll(ix,iy,iz,jxz) &
              + THIS% A_BSSN3_ll(ix,iy,iz,jyz)*THIS% A_BSSN3_ll(ix,iy,iz,jyz) &
              + THIS% A_BSSN3_ll(ix,iy,iz,jzz)*THIS% A_BSSN3_ll(ix,iy,iz,jzz) &
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
            HC_derphi( ix, iy, iz )= ( &
                            - DBLE(1.0/560.0)*EXP(THIS% phi(ix + 4,iy,iz)) &
                            + DBLE(8.0/315.0)*EXP(THIS% phi(ix + 3,iy,iz)) &
                            - DBLE(1.0/5.0  )*EXP(THIS% phi(ix + 2,iy,iz)) &
                            + DBLE(8.0/5.0  )*EXP(THIS% phi(ix + 1,iy,iz)) &
                            - DBLE(205.0/72.0)*EXP(THIS% phi(ix,iy,iz)) &
                            + DBLE(8.0/5.0  )*EXP(THIS% phi(ix - 1,iy,iz)) &
                            - DBLE(1.0/5.0  )*EXP(THIS% phi(ix - 2,iy,iz)) &
                            + DBLE(8.0/315.0)*EXP(THIS% phi(ix - 3,iy,iz)) &
                            - DBLE(1.0/560.0)*EXP(THIS% phi(ix - 4,iy,iz)) &
                            - DBLE(1.0/560.0)*EXP(THIS% phi(ix,iy + 4,iz)) &
                            + DBLE(8.0/315.0)*EXP(THIS% phi(ix,iy + 3,iz)) &
                            - DBLE(1.0/5.0  )*EXP(THIS% phi(ix,iy + 2,iz)) &
                            + DBLE(8.0/5.0  )*EXP(THIS% phi(ix,iy + 1,iz)) &
                            - DBLE(205.0/72.0)*EXP(THIS% phi(ix,iy,iz)) &
                            + DBLE(8.0/5.0  )*EXP(THIS% phi(ix,iy - 1,iz)) &
                            - DBLE(1.0/5.0  )*EXP(THIS% phi(ix,iy - 2,iz)) &
                            + DBLE(8.0/315.0)*EXP(THIS% phi(ix,iy - 3,iz)) &
                            - DBLE(1.0/560.0)*EXP(THIS% phi(ix,iy - 4,iz)) &
                            - DBLE(1.0/560.0)*EXP(THIS% phi(ix,iy,iz + 4)) &
                            + DBLE(8.0/315.0)*EXP(THIS% phi(ix,iy,iz + 3)) &
                            - DBLE(1.0/5.0  )*EXP(THIS% phi(ix,iy,iz + 2)) &
                            + DBLE(8.0/5.0  )*EXP(THIS% phi(ix,iy,iz + 1)) &
                            - DBLE(205.0/72.0)*EXP(THIS% phi(ix,iy,iz)) &
                            + DBLE(8.0/5.0  )*EXP(THIS% phi(ix,iy,iz - 1)) &
                            - DBLE(1.0/5.0  )*EXP(THIS% phi(ix,iy,iz - 2)) &
                            + DBLE(8.0/315.0)*EXP(THIS% phi(ix,iy,iz - 3)) &
                            - DBLE(1.0/560.0)*EXP(THIS% phi(ix,iy,iz - 4)) )&
                            /(THIS% dx**2)


            HC_hand( ix, iy, iz )= HC_rho( ix, iy, iz ) + &
                                   HC_trK( ix, iy, iz ) + &
                                   HC_A( ix, iy, iz )   + &
                                   HC_derphi( ix, iy, iz )

          ENDDO
        ENDDO
      ENDDO

    ENDIF

    !
    !-- Compute the BSSN constraints by calling the Cactus-bound procedure
    !-- BSSN_CONSTRAINTS_INTERIOR
    !
    imin = ghost_size
    imax(1) = THIS% ngrid_x - ghost_size - 1
    imax(2) = THIS% ngrid_y - ghost_size - 1
    imax(3) = THIS% ngrid_z - ghost_size - 1

    THIS% HC= 0.0D0
    THIS% MC= 0.0D0
    THIS% GC= 0.0D0
    PRINT *, "** Computing contraints..."
    CALL BSSN_CONSTRAINTS_INTERIOR( &
      !
      !-- Input
      !
      THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, &
      imin, imax, &
      THIS% dx, THIS% dy, THIS% dz, &
      THIS% g_BSSN3_ll(:,:,:,jxx), THIS% g_BSSN3_ll(:,:,:,jxy), &
      THIS% g_BSSN3_ll(:,:,:,jxz), THIS% g_BSSN3_ll(:,:,:,jyy), &
      THIS% g_BSSN3_ll(:,:,:,jyz), THIS% g_BSSN3_ll(:,:,:,jzz), &
      THIS% A_BSSN3_ll(:,:,:,jxx), THIS% A_BSSN3_ll(:,:,:,jxy), &
      THIS% A_BSSN3_ll(:,:,:,jxz), THIS% A_BSSN3_ll(:,:,:,jyy), &
      THIS% A_BSSN3_ll(:,:,:,jyz), THIS% A_BSSN3_ll(:,:,:,jzz), &
      THIS% trK(:,:,:), THIS% phi(:,:,:), &
      THIS% Gamma_u(:,:,:,jx), &
      THIS% Gamma_u(:,:,:,jy), &
      THIS% Gamma_u(:,:,:,jz), &
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
      THIS% lapse(:,:,:), &
      THIS% shift_u(:,:,:,jx), &
      THIS% shift_u(:,:,:,jy), &
      THIS% shift_u(:,:,:,jz), &
      !
      !-- Output
      !
      ! Connection constraints
      THIS% GC(:,:,:,jx), THIS% GC(:,:,:,jy), &
      THIS% GC(:,:,:,jz), &
      ! Hamiltonian and momentum constraints
      THIS% HC(:,:,:), &
      THIS% MC(:,:,:,jx), &
      THIS% MC(:,:,:,jy), &
      THIS% MC(:,:,:,jz) &
    )
    PRINT *, " * Constraints computed."
    PRINT *

    !---------------------------------------------------------!
    !--  Analyze constraints, and print to formatted files  --!
    !---------------------------------------------------------!

    !
    !-- Export the constraint statistics to a formatted file
    !
    unit_logfile= 2891

    INQUIRE( FILE= TRIM(name_logfile), EXIST= exist )

    IF( exist )THEN
        OPEN( UNIT= unit_logfile, FILE= TRIM(name_logfile), &
              STATUS= "REPLACE", &
              FORM= "FORMATTED", &
              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
              IOMSG= err_msg )
    ELSE
        OPEN( UNIT= unit_logfile, FILE= TRIM(name_logfile), &
              STATUS= "NEW", &
              FORM= "FORMATTED", &
              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening ", TRIM(name_logfile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when opening " &
    !         // TRIM(name_logfile) )

    WRITE( UNIT = unit_logfile, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id

    PRINT *, "** Analyzing constraints..."

    name_analysis= "bssn-hc-analysis.dat"
    name_constraint= "the Hamiltonian constraint"
    CALL THIS% analyze_constraint( &
         THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, &
         THIS% HC(:,:,:), name_constraint, unit_logfile, name_analysis )

    name_analysis= "bssn-mc1-analysis.dat"
    name_constraint= "the first component of the momentum constraint"
    CALL THIS% analyze_constraint( &
         THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, &
         THIS% MC(:,:,:,1), name_constraint, unit_logfile, name_analysis )

    name_analysis= "bssn-mc2-analysis.dat"
    name_constraint= "the second component of the momentum constraint"
    CALL THIS% analyze_constraint( &
         THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, &
         THIS% MC(:,:,:,2), name_constraint, unit_logfile, name_analysis )

    name_analysis= "bssn-mc3-analysis.dat"
    name_constraint= "the third component of the momentum constraint"
    CALL THIS% analyze_constraint( &
         THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, &
         THIS% MC(:,:,:,3), name_constraint, unit_logfile, name_analysis )

    name_analysis= "bssn-gc1-analysis.dat"
    name_constraint= "the first component of the connection constraint"
    CALL THIS% analyze_constraint( &
         THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, &
         THIS% GC(:,:,:,1), name_constraint, unit_logfile, name_analysis )

    name_analysis= "bssn-gc2-analysis.dat"
    name_constraint= "the second component of the connection constraint"
    CALL THIS% analyze_constraint( &
         THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, &
         THIS% GC(:,:,:,2), name_constraint, unit_logfile, name_analysis )

    name_analysis= "bssn-gc3-analysis.dat"
    name_constraint= "the third component of the connection constraint"
    CALL THIS% analyze_constraint( &
         THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, &
         THIS% GC(:,:,:,3), name_constraint, unit_logfile, name_analysis )

    CLOSE( UNIT= unit_logfile )

    PRINT *, " * Constraints analyzed. Results saved in " &
             // "bssn-constraints-statistics-*.log files."
    PRINT *

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
      "#      x   y   z   Stress-energy (10 components)   " &
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

      DO iz= 1, THIS% ngrid_z, 1
        DO iy= 1, THIS% ngrid_y, 1
          DO ix= 1, THIS% ngrid_x, 1
            abs_grid( 1, ix, iy, iz )= ABS( THIS% grid( 1, ix, iy, iz ) )
            abs_grid( 2, ix, iy, iz )= ABS( THIS% grid( 2, ix, iy, iz ) )
            abs_grid( 3, ix, iy, iz )= ABS( THIS% grid( 3, ix, iy, iz ) )
          ENDDO
        ENDDO
      ENDDO

      min_abs_y= 1D+20
      min_abs_z= 1D+20
      DO iz= 1, THIS% ngrid_z, 1
        DO iy= 1, THIS% ngrid_y, 1
          DO ix= 1, THIS% ngrid_x, 1
            IF( ABS( THIS% grid( 2, ix, iy, iz ) ) < min_abs_y )THEN
              min_abs_y= ABS( THIS% grid( 2, ix, iy, iz ) )
              min_ix_y= ix
              min_iy_y= iy
              min_iz_y= iz
            ENDIF
            IF( ABS( THIS% grid( 3, ix, iy, iz ) ) < min_abs_z )THEN
              min_abs_z= ABS( THIS% grid( 3, ix, iy, iz ) )
              min_ix_z= ix
              min_iy_z= iy
              min_iz_z= iz
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DO iz= 1, THIS% ngrid_z, 1

        IF( MOD( iz, THIS% cons_step ) /= 0 ) CYCLE

        DO iy= 1, THIS% ngrid_y, 1

          IF( MOD( iy, THIS% cons_step ) /= 0 ) CYCLE

          DO ix= 1, THIS% ngrid_x, 1

            IF( MOD( ix, THIS% cons_step ) /= 0 ) CYCLE

            IF( THIS% export_constraints_xy .AND. &
                THIS% grid( 3, ix, iy, iz ) /= &
                THIS% grid( 3, min_ix_z, min_iy_z, min_iz_z ) )THEN
              CYCLE
            ENDIF
            IF( THIS% export_constraints_x .AND. &
                ( THIS% grid( 3, ix, iy, iz ) /= &
                  THIS% grid( 3, min_ix_z, min_iy_z, min_iz_z ) &
                  .OR. &
                  THIS% grid( 2, ix, iy, iz ) /= &
                  THIS% grid( 2, min_ix_y, min_iy_y, min_iz_y ) ) )THEN
              CYCLE
            ENDIF

            IF( debug )THEN
              WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, &
                       FMT = * )&
                THIS% grid( 1, ix, iy, iz ), &
                THIS% grid( 2, ix, iy, iz ), &
                THIS% grid( 3, ix, iy, iz ), &
                baryon_density( ix, iy, iz ), &
                energy_density( ix, iy, iz ), &
                specific_energy( ix, iy, iz ), &
                pressure( ix, iy, iz ), &
                u_euler_l(ix,iy,iz,0), &
                u_euler_l(ix,iy,iz,1), &
                u_euler_l(ix,iy,iz,2), &
                u_euler_l(ix,iy,iz,3), &
                u_euler_l(ix,iy,iz,0), &
                u_euler_l(ix,iy,iz,1), &
                u_euler_l(ix,iy,iz,2), &
                u_euler_l(ix,iy,iz,3), &
                v_euler(ix,iy,iz,1), &
                v_euler(ix,iy,iz,2), &
                v_euler(ix,iy,iz,3), &
                Tmunu_ll( ix, iy, iz, itt ), &
                Tmunu_ll( ix, iy, iz, itx ), &
                Tmunu_ll( ix, iy, iz, ity ), &
                Tmunu_ll( ix, iy, iz, itz ), &
                Tmunu_ll( ix, iy, iz, ixx ), &
                Tmunu_ll( ix, iy, iz, ixy ), &
                Tmunu_ll( ix, iy, iz, ixz ), &
                Tmunu_ll( ix, iy, iz, iyy ), &
                Tmunu_ll( ix, iy, iz, iyz ), &
                Tmunu_ll( ix, iy, iz, izz ), &
                THIS% HC( ix, iy, iz ), &
                HC_hand( ix, iy, iz ), &
                HC_rho( ix, iy, iz ), &
                HC_trK( ix, iy, iz ), &
                HC_A( ix, iy, iz ), &
                HC_derphi( ix, iy, iz ), &
                lorentz_factor( ix, iy, iz ), &
                THIS% lapse( ix, iy, iz ), &
                THIS% shift_u( ix, iy, iz, jx ), &
                THIS% shift_u( ix, iy, iz, jy ), &
                THIS% shift_u( ix, iy, iz, jz ), &
                g4( ix, iy, iz, ixx ), &
                g4( ix, iy, iz, ixy ), &
                g4( ix, iy, iz, ixz ), &
                g4( ix, iy, iz, iyy ), &
                g4( ix, iy, iz, iyz ), &
                g4( ix, iy, iz, izz ), &
                !THIS% g_BSSN3_ll( ix, iy, iz, jxx ), &
                !THIS% g_BSSN3_ll( ix, iy, iz, jxy ), &
                !THIS% g_BSSN3_ll( ix, iy, iz, jxz ), &
                !THIS% g_BSSN3_ll( ix, iy, iz, jyy ), &
                !THIS% g_BSSN3_ll( ix, iy, iz, jyz ), &
                !THIS% g_BSSN3_ll( ix, iy, iz, jzz ), &
                THIS% k_phys3_ll( ix, iy, iz, jxx ), &
                THIS% k_phys3_ll( ix, iy, iz, jxy ), &
                THIS% k_phys3_ll( ix, iy, iz, jxz ), &
                THIS% k_phys3_ll( ix, iy, iz, jyy ), &
                THIS% k_phys3_ll( ix, iy, iz, jyz ), &
                THIS% k_phys3_ll( ix, iy, iz, jzz ), &
                THIS% A_BSSN3_ll( ix, iy, iz, jxx ), &
                THIS% A_BSSN3_ll( ix, iy, iz, jxy ), &
                THIS% A_BSSN3_ll( ix, iy, iz, jxz ), &
                THIS% A_BSSN3_ll( ix, iy, iz, jyy ), &
                THIS% A_BSSN3_ll( ix, iy, iz, jyz ), &
                THIS% A_BSSN3_ll( ix, iy, iz, jzz ), &
                THIS% trK( ix, iy, iz ), &
                THIS% phi( ix, iy, iz ), &
                THIS% Gamma_u( ix, iy, iz, 1 ), &
                THIS% Gamma_u( ix, iy, iz, 2 ), &
                THIS% Gamma_u( ix, iy, iz, 3 )
            ELSE
              WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * )&
                THIS% grid( 1, ix, iy, iz ), &
                THIS% grid( 2, ix, iy, iz ), &
                THIS% grid( 3, ix, iy, iz ), &
                Tmunu_ll( ix, iy, iz, itt ), &
                Tmunu_ll( ix, iy, iz, itx ), &
                Tmunu_ll( ix, iy, iz, ity ), &
                Tmunu_ll( ix, iy, iz, itz ), &
                Tmunu_ll( ix, iy, iz, ixx ), &
                Tmunu_ll( ix, iy, iz, ixy ), &
                Tmunu_ll( ix, iy, iz, ixz ), &
                Tmunu_ll( ix, iy, iz, iyy ), &
                Tmunu_ll( ix, iy, iz, iyz ), &
                Tmunu_ll( ix, iy, iz, izz ), &
                THIS% HC( ix, iy, iz ), &
                THIS% MC( ix, iy, iz, 1 ), &
                THIS% MC( ix, iy, iz, 2 ), &
                THIS% MC( ix, iy, iz, 3 ), &
                THIS% Gamma_u( ix, iy, iz, 1 ), &
                THIS% Gamma_u( ix, iy, iz, 2 ), &
                THIS% Gamma_u( ix, iy, iz, 3 )
            ENDIF

            IF( ios > 0 )THEN
              PRINT *, "...error when writing the arrays in ", TRIM(namefile), &
                       ". The error message is", err_msg
              STOP
            ENDIF
            !CALL test_status( ios, err_msg, &
            !          "...error in writing " &
            !          // "the arrays in " // TRIM(namefile) )
          ENDDO
        ENDDO
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

    USE constants,  ONLY: c_light2, cm2m, MSun, g2kg, m2cm, Msun_geo
    USE units,      ONLY: set_units, m0c2_cu
    USE grav_grid,  ONLY: ngrid_x, ngrid_y, ngrid_z, dx, dy, dz, &
                          dx_1, dy_1, dz_1, rad_coord, xR, xL, yR, yL, zR, zL, &
                          deallocate_gravity_grid
    USE tensor,     ONLY: itt, itx, ity, itz, ixx, ixy, &
                          ixz, iyy, iyz, izz, jxx, jxy, jxz, &
                          jyy, jyz, jzz, jx, jy, jz, n_sym3x3, n_sym4x4
    USE McLachlan,  ONLY: ghost_size, BSSN_CONSTRAINTS_INTERIOR, &
                          allocate_Ztmp, deallocate_Ztmp
    USE ADM,                  ONLY: lapse, dt_lapse, shift_u, dt_shift_u, &
                                    K_phys3_ll, g_phys3_ll, &
                                    allocate_ADM, deallocate_ADM
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
    USE map_particles_2_grid, ONLY: map_2_grid
    USE metric_on_particles,  ONLY: allocate_metric_on_particles, &
                                    deallocate_metric_on_particles, &
                                    get_metric_on_particles
    USE particle_mesh,        ONLY: deallocate_all_lists, &
                                    deallocate_flag_nei_cell, &
                                    deallocate_pp_g
    USE Tmunu,                ONLY: Tmunu_ll, allocate_Tmunu, deallocate_Tmunu
    USE GravityAcceleration,  ONLY: allocate_GravityAcceleration, &
                                    deallocate_GravityAcceleration
    USE alive_flag,           ONLY: alive

    IMPLICIT NONE

    INTEGER:: ix, iy, iz, a, allocation_status
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

    LOGICAL:: exist
    LOGICAL, PARAMETER:: debug= .FALSE.

    ! Being abs_grid a local array, it is good practice to allocate it on the
    ! heap, otherwise it will be stored on the stack which has a very limited
    ! size. This results in a segmentation fault.
    ALLOCATE( abs_grid( 3, THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z ) )

    ngrid_x= THIS% ngrid_x
    ngrid_y= THIS% ngrid_y
    ngrid_z= THIS% ngrid_z
    dx= THIS% dx
    dy= THIS% dy
    dz= THIS% dz
    dx_1= THIS% dx_1
    dy_1= THIS% dy_1
    dz_1= THIS% dz_1

    ALLOCATE( THIS% HC_parts( THIS% ngrid_x, THIS% ngrid_y, &
                                   THIS% ngrid_z ) )
    ALLOCATE( THIS% MC_parts( THIS% ngrid_x, THIS% ngrid_y, &
                                   THIS% ngrid_z, 3 ) )
    ALLOCATE( THIS% GC_parts( THIS% ngrid_x, THIS% ngrid_y, &
                                   THIS% ngrid_z, 3 ) )

    PRINT *, "Mapping hydro fields from particles to grid..."

    CALL allocate_ADM()
    !CALL allocate_BSSN()

    ! Allocate temporary memory for time integration
    CALL allocate_Ztmp()

    ! Allocate memory for the stress-energy tensor (used in
    ! write_BSSN_dump)
    CALL allocate_Tmunu()

    ! Allocate memory for the derivatives of the ADM variables
    CALL allocate_GravityAcceleration()

    ! Initialize the stress-energy tensor to 0
    Tmunu_ll(:,:,:,:)= 0.0D0

    !
    !-- Allocate memory for the coordinate radius
    !-- (this is needed by ADM_to_BSSN)
    !
    IF( .NOT. ALLOCATED( rad_coord ) )THEN
        ALLOCATE( rad_coord( THIS% ngrid_x, THIS% ngrid_y, &
                             THIS% ngrid_z ), STAT= allocation_status )
    ENDIF
    IF( allocation_status > 0 )THEN
       PRINT *, '...allocation error for rad_coord'
       STOP
    ENDIF

    DO iz= 1, THIS% ngrid_z
      DO iy= 1, THIS% ngrid_y
        DO ix= 1, THIS% ngrid_x
          rad_coord( ix, iy, iz )= SQRT( (xL+(ix-1)*dx)**2 &
                                       + (yL+(iy-1)*dy)**2 &
                                       + (zL+(iz-1)*dz)**2 )
        ENDDO
      ENDDO
    ENDDO

    dt_lapse  = 0.0D0
    dt_shift_u= 0.0D0
    lapse= THIS% lapse
    shift_u= THIS% shift_u
    g_phys3_ll= THIS% g_phys3_ll

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
    ALLOCATE( abs_pos( 3, npart ) )

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

    IF( debug .AND. .FALSE. ) PRINT *, "npart= ", npart
    IF( debug .AND. .FALSE. ) PRINT *, "nu_loc= ", nu_loc(npart/2)
    IF( debug .AND. .FALSE. ) PRINT *, "pos_loc= ", pos_loc(2,npart/2)
    IF( debug .AND. .FALSE. ) PRINT *, "vel_loc= ", vel_loc(2,npart/2)
    IF( debug .AND. .FALSE. ) PRINT *, "u_loc= ", u_loc(npart/2)
    IF( debug .AND. .FALSE. ) PRINT *, "nlrf_loc= ", nlrf_loc(npart/2)
    IF( debug .AND. .FALSE. ) PRINT *, "theta_loc= ", theta_loc(npart/2)
    IF( debug .AND. .FALSE. ) PRINT *, "pressure_loc= ", pressure_loc(npart/2)
    IF( debug .AND. .FALSE. ) PRINT *

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
    imin = ghost_size
    imax(1) = THIS% ngrid_x - ghost_size - 1
    imax(2) = THIS% ngrid_y - ghost_size - 1
    imax(3) = THIS% ngrid_z - ghost_size - 1

    THIS% HC_parts= 0.0
    THIS% MC_parts= 0.0
    THIS% GC_parts= 0.0

    PRINT *, " * Computing contraints using particle data..."
    CALL BSSN_CONSTRAINTS_INTERIOR( &
      !
      !-- Input
      !
      THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, &
      imin, imax, &
      THIS% dx, THIS% dy, THIS% dz, &
      THIS% g_BSSN3_ll(:,:,:,jxx), THIS% g_BSSN3_ll(:,:,:,jxy), &
      THIS% g_BSSN3_ll(:,:,:,jxz), THIS% g_BSSN3_ll(:,:,:,jyy), &
      THIS% g_BSSN3_ll(:,:,:,jyz), THIS% g_BSSN3_ll(:,:,:,jzz), &
      THIS% A_BSSN3_ll(:,:,:,jxx), THIS% A_BSSN3_ll(:,:,:,jxy), &
      THIS% A_BSSN3_ll(:,:,:,jxz), THIS% A_BSSN3_ll(:,:,:,jyy), &
      THIS% A_BSSN3_ll(:,:,:,jyz), THIS% A_BSSN3_ll(:,:,:,jzz), &
      THIS% trK(:,:,:), THIS% phi(:,:,:), &
      THIS% Gamma_u(:,:,:,jx), &
      THIS% Gamma_u(:,:,:,jy), THIS% Gamma_u(:,:,:,jz), &
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
      THIS% lapse(:,:,:), THIS% shift_u(:,:,:,jx), &
      THIS% shift_u(:,:,:,jy), THIS% shift_u(:,:,:,jz), &
      !
      !-- Output
      !
      ! Connection constraints
      THIS% GC_parts(:,:,:,jx), &
      THIS% GC_parts(:,:,:,jy), &
      THIS% GC_parts(:,:,:,jz), &
      ! Hamiltonian and momentum constraints
      THIS% HC_parts(:,:,:), &
      THIS% MC_parts(:,:,:,jx), &
      THIS% MC_parts(:,:,:,jy), &
      THIS% MC_parts(:,:,:,jz) &
    )
    PRINT *, " * Constraints computed."
    PRINT *

    IF( debug ) PRINT *, "0"

    !---------------------------------------------------------!
    !--  Analyze constraints, and print to formatted files  --!
    !---------------------------------------------------------!

    unit_logfile= 2891

    INQUIRE( FILE= TRIM(name_logfile), EXIST= exist )

    IF( debug ) PRINT *, "1"

    IF( exist )THEN
        OPEN( UNIT= unit_logfile, FILE= TRIM(name_logfile), &
              STATUS= "REPLACE", &
              FORM= "FORMATTED", &
              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
              IOMSG= err_msg )
    ELSE
        OPEN( UNIT= unit_logfile, FILE= TRIM(name_logfile), &
              STATUS= "NEW", &
              FORM= "FORMATTED", &
              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening ", TRIM(name_logfile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when opening " &
    !                  // TRIM(name_logfile) )

    IF( debug ) PRINT *, "2"

    WRITE( UNIT = unit_logfile, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id

    IF( debug ) PRINT *, "3"

    name_analysis= "bssn-hc-parts-analysis.dat"
    name_constraint= "the Hamiltonian constraint"
    CALL THIS% analyze_constraint( &
         THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, &
         THIS% HC_parts(:,:,:), name_constraint, unit_logfile, &
         name_analysis )

    name_analysis= "bssn-mc1-parts-analysis.dat"
    name_constraint= "the first component of the momentum constraint"
    CALL THIS% analyze_constraint( &
         THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, &
         THIS% MC_parts(:,:,:,1), name_constraint, unit_logfile, &
         name_analysis )

    name_analysis= "bssn-mc2-parts-analysis.dat"
    name_constraint= "the second component of the momentum constraint"
    CALL THIS% analyze_constraint( &
         THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, &
         THIS% MC_parts(:,:,:,2), name_constraint, unit_logfile, &
         name_analysis )

    name_analysis= "bssn-mc3-parts-analysis.dat"
    name_constraint= "the third component of the momentum constraint"
    CALL THIS% analyze_constraint( &
         THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, &
         THIS% MC_parts(:,:,:,3), name_constraint, unit_logfile, &
         name_analysis )

    name_analysis= "bssn-gc1-parts-analysis.dat"
    name_constraint= "the first component of the connection constraint"
    CALL THIS% analyze_constraint( &
         THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, &
         THIS% GC_parts(:,:,:,1), name_constraint, unit_logfile, &
         name_analysis )

    name_analysis= "bssn-gc2-parts-analysis.dat"
    name_constraint= "the second component of the connection constraint"
    CALL THIS% analyze_constraint( &
         THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, &
         THIS% GC_parts(:,:,:,2), name_constraint, unit_logfile, &
         name_analysis )

    name_analysis= "bssn-gc3-parts-analysis.dat"
    name_constraint= "the third component of the connection constraint"
    CALL THIS% analyze_constraint( &
         THIS% ngrid_x, THIS% ngrid_y, THIS% ngrid_z, &
         THIS% GC_parts(:,:,:,3), name_constraint, unit_logfile, &
         name_analysis )

    IF( debug ) PRINT *, "4"

    CLOSE( UNIT= unit_logfile )

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
      "# Values of the BSSN constraints for the LORENE ID "&
      // "on selected grid points"
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
      "#      x   y   z   Stress-energy (10 components)   " &
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

      DO iz= 1, THIS% ngrid_z, 1
        DO iy= 1, THIS% ngrid_y, 1
          DO ix= 1, THIS% ngrid_x, 1
            abs_grid( 1, ix, iy, iz )= ABS( THIS% grid( 1, ix, iy, iz ) )
            abs_grid( 2, ix, iy, iz )= ABS( THIS% grid( 2, ix, iy, iz ) )
            abs_grid( 3, ix, iy, iz )= ABS( THIS% grid( 3, ix, iy, iz ) )
          ENDDO
        ENDDO
      ENDDO

      min_abs_y= 1D+20
      min_abs_z= 1D+20
      DO iz= 1, THIS% ngrid_z, 1
        DO iy= 1, THIS% ngrid_y, 1
          DO ix= 1, THIS% ngrid_x, 1
            IF( ABS( THIS% grid( 2, ix, iy, iz ) ) < min_abs_y )THEN
              min_abs_y= ABS( THIS% grid( 2, ix, iy, iz ) )
              min_ix_y= ix
              min_iy_y= iy
              min_iz_y= iz
            ENDIF
            IF( ABS( THIS% grid( 3, ix, iy, iz ) ) < min_abs_z )THEN
              min_abs_z= ABS( THIS% grid( 3, ix, iy, iz ) )
              min_ix_z= ix
              min_iy_z= iy
              min_iz_z= iz
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DO iz= 1, THIS% ngrid_z, 1

        IF( MOD( iz, THIS% cons_step ) /= 0 ) CYCLE

        DO iy= 1, THIS% ngrid_y, 1

          IF( MOD( iy, THIS% cons_step ) /= 0 ) CYCLE

          DO ix= 1, THIS% ngrid_x, 1

            IF( MOD( ix, THIS% cons_step ) /= 0 ) CYCLE

            IF( THIS% export_constraints_xy .AND. &
                THIS% grid( 3, ix, iy, iz ) /= &
                THIS% grid( 3, min_ix_z, min_iy_z, min_iz_z ) )THEN
              CYCLE
            ENDIF
            IF( THIS% export_constraints_x .AND. &
                ( THIS% grid( 3, ix, iy, iz ) /= &
                  THIS% grid( 3, min_ix_z, min_iy_z, min_iz_z ) &
                  .OR. &
                  THIS% grid( 2, ix, iy, iz ) /= &
                  THIS% grid( 2, min_ix_y, min_iy_y, min_iz_y ) ) )THEN
              CYCLE
            ENDIF

            IF( debug )THEN
            WRITE( UNIT = 21, IOSTAT = ios, IOMSG = err_msg, FMT = * )&
              THIS% grid( 1, ix, iy, iz ), &
              THIS% grid( 2, ix, iy, iz ), &
              THIS% grid( 3, ix, iy, iz ), &
              THIS% grid( 1, ix, iy, iz ), &
              THIS% grid( 2, ix, iy, iz ), &
              THIS% grid( 3, ix, iy, iz ), &
              THIS% grid( 1, ix, iy, iz ), &
              THIS% grid( 2, ix, iy, iz ), &
              THIS% grid( 3, ix, iy, iz ), &
              THIS% grid( 1, ix, iy, iz ), &
              THIS% grid( 2, ix, iy, iz ), &
              THIS% grid( 3, ix, iy, iz ), &
              THIS% grid( 1, ix, iy, iz ), &
              THIS% grid( 2, ix, iy, iz ), &
              THIS% grid( 3, ix, iy, iz ), &
              THIS% grid( 1, ix, iy, iz ), &
              THIS% grid( 2, ix, iy, iz ), &
              THIS% grid( 3, ix, iy, iz ), & ! column 18
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
              Tmunu_ll( ix, iy, iz, itt ), &
              Tmunu_ll( ix, iy, iz, itx ), &
              Tmunu_ll( ix, iy, iz, ity ), &
              Tmunu_ll( ix, iy, iz, itz ), &
              Tmunu_ll( ix, iy, iz, ixx ), &
              Tmunu_ll( ix, iy, iz, ixy ), &
              Tmunu_ll( ix, iy, iz, ixz ), &
              Tmunu_ll( ix, iy, iz, iyy ), &
              Tmunu_ll( ix, iy, iz, iyz ), &
              Tmunu_ll( ix, iy, iz, izz ), &
              THIS% HC_parts( ix, iy, iz ), &
              THIS% MC_parts( ix, iy, iz, jx ), &
              THIS% MC_parts( ix, iy, iz, jy ), &
              THIS% MC_parts( ix, iy, iz, jz ), &
              THIS% GC_parts( ix, iy, iz, jx ), &
              THIS% GC_parts( ix, iy, iz, jy ), &
              THIS% GC_parts( ix, iy, iz, jz ), &
              THIS% lapse( ix, iy, iz ), &
              THIS% shift_u( ix, iy, iz, jx ), &
              THIS% shift_u( ix, iy, iz, jy ), &
              THIS% shift_u( ix, iy, iz, jz ), &
              THIS% g_BSSN3_ll( ix, iy, iz, jxx ), &
              THIS% g_BSSN3_ll( ix, iy, iz, jxy ), &
              THIS% g_BSSN3_ll( ix, iy, iz, jxz ), &
              THIS% g_BSSN3_ll( ix, iy, iz, jyy ), &
              THIS% g_BSSN3_ll( ix, iy, iz, jyz ), &
              THIS% g_BSSN3_ll( ix, iy, iz, jzz ), &
              THIS% k_phys3_ll( ix, iy, iz, jxx ), &
              THIS% k_phys3_ll( ix, iy, iz, jxy ), &
              THIS% k_phys3_ll( ix, iy, iz, jxz ), &
              THIS% k_phys3_ll( ix, iy, iz, jyy ), &
              THIS% k_phys3_ll( ix, iy, iz, jyz ), &
              THIS% k_phys3_ll( ix, iy, iz, jzz ), &
              THIS% A_BSSN3_ll( ix, iy, iz, jxx ), &
              THIS% A_BSSN3_ll( ix, iy, iz, jxy ), &
              THIS% A_BSSN3_ll( ix, iy, iz, jxz ), &
              THIS% A_BSSN3_ll( ix, iy, iz, jyy ), &
              THIS% A_BSSN3_ll( ix, iy, iz, jyz ), &
              THIS% A_BSSN3_ll( ix, iy, iz, jzz ), &
              THIS% trK( ix, iy, iz ), &
              THIS% phi( ix, iy, iz ), &
              THIS% Gamma_u(ix,iy,iz,1), &
              THIS% Gamma_u(ix,iy,iz,2), &
              THIS% Gamma_u(ix,iy,iz,3)
            ELSE
              WRITE( UNIT = 21, IOSTAT = ios, IOMSG = err_msg, FMT = * )&
                THIS% grid( 1, ix, iy, iz ), &
                THIS% grid( 2, ix, iy, iz ), &
                THIS% grid( 3, ix, iy, iz ), &
                Tmunu_ll( ix, iy, iz, itt ), &
                Tmunu_ll( ix, iy, iz, itx ), &
                Tmunu_ll( ix, iy, iz, ity ), &
                Tmunu_ll( ix, iy, iz, itz ), &
                Tmunu_ll( ix, iy, iz, ixx ), &
                Tmunu_ll( ix, iy, iz, ixy ), &
                Tmunu_ll( ix, iy, iz, ixz ), &
                Tmunu_ll( ix, iy, iz, iyy ), &
                Tmunu_ll( ix, iy, iz, iyz ), &
                Tmunu_ll( ix, iy, iz, izz ), &
                THIS% HC_parts( ix, iy, iz ), &
                THIS% MC_parts( ix, iy, iz, jx ), &
                THIS% MC_parts( ix, iy, iz, jy ), &
                THIS% MC_parts( ix, iy, iz, jz ), &
                THIS% GC_parts( ix, iy, iz, jx ), &
                THIS% GC_parts( ix, iy, iz, jy ), &
                THIS% GC_parts( ix, iy, iz, jz )
            ENDIF

            IF( ios > 0 )THEN
              PRINT *, "...error when writing the arrays in ", TRIM(namefile), &
                       ". The error message is", err_msg
              STOP
            ENDIF
            !CALL test_status( ios, err_msg, &
            !                  "...error in writing " &
            !                  // "the arrays in " // TRIM(namefile) )
          ENDDO
        ENDDO
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

      DO itr = 1, npart, 1
        abs_pos( 1, itr )= ABS( pos_loc( 1, itr ) )
        abs_pos( 2, itr )= ABS( pos_loc( 2, itr ) )
        abs_pos( 3, itr )= ABS( pos_loc( 3, itr ) )
      ENDDO

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
    !CALL deallocate_BSSN()
    CALL deallocate_gravity_grid()

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

    IMPLICIT NONE

    IF(ALLOCATED( THIS% Gamma_u ))THEN
      DEALLOCATE( THIS% Gamma_u, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...deallocation error for array Gamma_u. ", &
                 "The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array Gamma_u" )
    ENDIF
    IF(ALLOCATED( THIS% phi ))THEN
      DEALLOCATE( THIS% phi, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...deallocation error for array phi. ", &
                 "The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array phi" )
    ENDIF
    IF(ALLOCATED( THIS% trK ))THEN
      DEALLOCATE( THIS% trK, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...deallocation error for array trK. ", &
                 "The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array trK" )
    ENDIF
    IF(ALLOCATED( THIS% A_BSSN3_ll ))THEN
      DEALLOCATE( THIS% A_BSSN3_ll, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...deallocation error for array A_BSSN3_ll. ", &
                 "The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array A_BSSN3_ll" )
    ENDIF
    IF(ALLOCATED( THIS% g_BSSN3_ll ))THEN
      DEALLOCATE( THIS% g_BSSN3_ll, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...deallocation error for array g_BSSN3_ll. ", &
                 "The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_BSSN3_ll" )
    ENDIF
    IF(ALLOCATED( THIS% GC ))THEN
      DEALLOCATE( THIS% GC, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...deallocation error for array GC. ", &
                 "The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array GC" )
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
