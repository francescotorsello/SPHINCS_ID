! File:         submodule_bssn_formulation_bssn_variables.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (bssn_formulation) bssn_variables

  !************************************************
  !
  !# Implementation of the methods of TYPE bssn
  !  that compute the |bssn| variables
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


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE compute_and_export_bssn_variables

    !************************************************
    !
    !# Compute, stores and prints the BSSN variables
    !  to a binary file to be read by the evolution
    !  code SPHINCS_BSSN
    !
    !  FT 23.10.2020
    !
    !************************************************

    !USE NaNChecker,          ONLY: Check_Grid_Function_for_NAN
    USE mesh_refinement,            ONLY: nlevels, levels, rad_coord, &
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
    USE GravityAcceleration_refine, ONLY: allocate_GravityAcceleration, &
                                          deallocate_GravityAcceleration
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
                           !Theta_Z4,         & ! Vector in the CCZ4 formulation
                                               ! Used because ADM_TO_BSSN
                                               ! calls SUBROUTINES that need it
                                               ! as input; however, it is not
                                               ! evolved in BSSN
                           !lapse_A_BSSN,     & ! Time derivative of lapse
                           !shift_B_BSSN_u,   & ! Time derivativeof shift
                           write_BSSN_dump

    IMPLICIT NONE

    ! The flag call_flag is set different than 0 if the SUBROUTINE
    ! compute_and_export_SPH_variables is called
    INTEGER, SAVE:: call_flag= 0

    INTEGER:: l

    PRINT *, "** Computing and exporting BSSN ID..."

    ! Allocate memory for the ADM MODULE variables (this has to be done since
    ! the MODULE SUBROUTINES need them; not allocating it results in a
    ! segmentation fault)
    PRINT *
    PRINT *, " * Allocating needed memory..."
    PRINT *

    ALLOCATE ( levels( THIS% nlevels ), STAT=ios )
    IF( ios > 0 )THEN
     PRINT*,'...allocation error for levels'
     STOP
    ENDIF
    levels = THIS% levels
    nlevels= THIS% nlevels

    !DO l= 1, THIS% nlevels, 1
    !  levels(l)% ngrid_x= THIS% levels(l)% ngrid_x
    !  levels(l)% ngrid_x= THIS% levels(l)% ngrid_x
    !  levels(l)% ngrid_x= THIS% levels(l)% ngrid_x
    !ENDDO

    CALL allocate_ADM()
    CALL allocate_BSSN()

    ! Allocate temporary memory for time integration
    CALL allocate_Ztmp()

    ! Allocate memory for the stress-energy tensor (used in write_BSSN_dump)
    CALL allocate_Tmunu()

    ! Allocate memory for the derivatives of the ADM variables
    ! Here there is not enough memory for 7rl with 291**3 points
    CALL allocate_GravityAcceleration()

    CALL allocate_grid_function ( rad_coord, 'rad_coord' )

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

    !CALL deallocate_ADM()
    CALL deallocate_Ztmp()
    CALL deallocate_Tmunu()
    CALL deallocate_GravityAcceleration()

    !
    !-- Setting the TYPE variables equal to the MODULE variables
    !
    CALL allocate_bssn_fields( THIS )
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
      IF( PRESENT(namefile) )THEN
        CALL write_BSSN_dump( namefile )
        !CALL write_BSSN_dump()
      ELSE
        CALL write_BSSN_dump()
      ENDIF
    ENDIF

    !
    !-- Deallocate MODULE variables
    !
    CALL deallocate_ADM()
    !CALL deallocate_Ztmp()
    !CALL deallocate_Tmunu()
    !CALL deallocate_GravityAcceleration()
    CALL deallocate_BSSN()
    CALL deallocate_grid_function( rad_coord, 'rad_coord' )
    !CALL deallocate_gravity_grid()
    DEALLOCATE( levels )

    call_flag= call_flag + 1
    THIS% call_flag= call_flag

    PRINT *, "** BSSN ID computed."
    PRINT *

  END PROCEDURE compute_and_export_bssn_variables


  MODULE PROCEDURE read_bssn_dump_print_formatted

    !************************************************
    !
    !# Read the BSSN ID from the binary file output
    !  by write_BSSN_dump, and print it to a
    !  formatted file
    !
    !  FT 08.02.2021
    !
    !************************************************

    USE mesh_refinement,  ONLY: levels, nlevels
    USE tensor,           ONLY: jxx, jxy, jxz, jyy, jyz, jzz, jx, jy, jz
    USE ADM_refine,       ONLY: lapse, shift_u, &
                                   allocate_ADM, deallocate_ADM
    USE BSSN_refine,  ONLY: allocate_BSSN, deallocate_BSSN, &
                            Gamma_u,        & ! Conformal connection
                            phi,            & ! Conformal factor
                            trK,            & ! Trace of extrinsic curvature
                            A_BSSN3_ll,     & ! Conformal traceless
                                              ! extrinsic curvature
                            g_BSSN3_ll,     & ! Conformal metric
                            !Theta_Z4,       & ! Vector in the CCZ4 formulation.
                                              ! Loaded here because ADM_TO_BSSN
                                              ! calls SUBROUTINES that need it
                                              ! as input; however, it is not
                                              ! evolved in BSSN
                            !lapse_A_BSSN,   & ! Time derivative of lapse
                            !shift_B_BSSN_u, & ! Time derivativeof shift
                            read_BSSN_dump

    IMPLICIT NONE

    INTEGER:: i, j, k, l, min_ix_y, min_iy_y, min_iz_y, &
              min_ix_z, min_iy_z, min_iz_z

    DOUBLE PRECISION:: min_abs_y, min_abs_z

    LOGICAL:: exist

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    PRINT *, "** Executing the read_bssn_dump_print_formatted subroutine..."

    levels = THIS% levels
    nlevels= THIS% nlevels

    CALL allocate_ADM()
    CALL allocate_BSSN()

    CALL read_BSSN_dump( 00000, namefile_bin )

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
    "#      refinement level    x [km]       y [km]       z [km]       lapse", &
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

    DO l= 1, THIS% nlevels, 1

      ASSOCIATE( lapse      => lapse% levels(l)% var, &
                 shift_u    => shift_u% levels(l)% var, &
                 g_BSSN3_ll => g_BSSN3_ll% levels(l)% var, &
                 A_BSSN3_ll => A_BSSN3_ll% levels(l)% var, &
                 phi        => phi% levels(l)% var, &
                 trK        => trK% levels(l)% var, &
                 Gamma_u    => Gamma_u% levels(l)% var &
      )

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

          DO j= 1, THIS% get_ngrid_y(l), 1

            DO i= 1, THIS% get_ngrid_x(l), 1

              IF( THIS% export_form_xy .AND. &
                  ( THIS% coords% levels(l)% var( i, j, k, jz ) /= &
                    THIS% coords% levels(l)% var( min_ix_z, min_iy_z, &
                                                  min_iz_z, jz ) ) )THEN
                CYCLE
              ENDIF
              IF( THIS% export_form_x .AND. &
                  ( THIS% coords% levels(l)% var( i, j, k, jz ) /= &
                    THIS% coords% levels(l)% var( min_ix_z, min_iy_z, &
                                                  min_iz_z, jz ) &
                    .OR. &
                    THIS% coords% levels(l)% var( i, j, k, jy ) /= &
                    THIS% coords% levels(l)% var( min_ix_y, min_iy_y, &
                                                  min_iz_y, jy ) ) )THEN
                CYCLE
              ENDIF

              WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * )&
                  l, &
                  THIS% coords% levels(l)% var( i, j, k, jx ), &
                  THIS% coords% levels(l)% var( i, j, k, jy ), &
                  THIS% coords% levels(l)% var( i, j, k, jz ), &
                  lapse( i, j, k ), &
                  shift_u( i, j, k, jx ), &
                  shift_u( i, j, k, jy ), &
                  shift_u( i, j, k, jz ), &
                  phi( i, j, k ), &
                  trK( i, j, k ), &
                  g_BSSN3_ll( i, j, k, jxx ), &
                  g_BSSN3_ll( i, j, k, jxy ), &
                  g_BSSN3_ll( i, j, k, jxz ), &
                  g_BSSN3_ll( i, j, k, jyy ), &
                  g_BSSN3_ll( i, j, k, jyz ), &
                  g_BSSN3_ll( i, j, k, jzz ), &
                  A_BSSN3_ll( i, j, k, jxx ), &
                  A_BSSN3_ll( i, j, k, jxy ), &
                  A_BSSN3_ll( i, j, k, jxz ), &
                  A_BSSN3_ll( i, j, k, jyy ), &
                  A_BSSN3_ll( i, j, k, jyz ), &
                  A_BSSN3_ll( i, j, k, jzz ), &
                  Gamma_u( i, j, k, jx ), &
                  Gamma_u( i, j, k, jy ), &
                  Gamma_u( i, j, k, jz )

              IF( ios > 0 )THEN
                PRINT *, "...error when writing the arrays in ", &
                         TRIM(finalnamefile), ". The error message is", err_msg
                STOP
              ENDIF
              !CALL test_status( ios, err_msg, "...error when writing " &
              !                  // "the arrays in " // TRIM(namefile) )

            ENDDO
          ENDDO
        ENDDO
      END ASSOCIATE
    ENDDO

    CLOSE( UNIT= 20 )

    !
    !-- Deallocate MODULE variables
    !
    CALL deallocate_ADM()
    CALL deallocate_BSSN()

    PRINT *, " * LORENE BSSN ID on the refined mesh, to be supplied to ", &
             "SPHINCS_BSSN, printed to formatted file ", TRIM(namefile)

    PRINT *, "** Subroutine read_bssn_dump_print_formatted " &
             // "executed."
    PRINT *

  END PROCEDURE read_bssn_dump_print_formatted


  MODULE PROCEDURE print_formatted_lorene_id_bssn_variables

    !************************************************
    !
    !# Print the BSSN ID, computed on the gravity
    !  grid, to a formatted file
    !
    !  FT 26.10.2020
    !
    !************************************************

    USE tensor,              ONLY: itt, itx, ity, itz, ixx, ixy, &
                                   ixz, iyy, iyz, izz, jxx, jxy, jxz, &
                                   jyy, jyz, jzz, jx, jy, jz

    IMPLICIT NONE

    INTEGER:: i, j, k, l, min_ix_y, min_iy_y, min_iz_y, &
              min_ix_z, min_iy_z, min_iz_z

    DOUBLE PRECISION:: min_abs_y, min_abs_z

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
    "       23      24    25"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 2 in ", TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 2 in "&
    !        // TRIM(finalnamefile) )
    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "#      refinement level    x [km]       y [km]       z [km]       lapse", &
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

    DO l= 1, THIS% nlevels, 1

      ASSOCIATE( lapse           => THIS% lapse% levels(l)% var, &
                 shift_u         => THIS% shift_u% levels(l)% var, &
                 phi             => THIS% phi% levels(l)% var, &
                 trK             => THIS% trK% levels(l)% var, &
                 g_BSSN3_ll      => THIS% g_BSSN3_ll% levels(l)% var, &
                 A_BSSN3_ll      => THIS% A_BSSN3_ll% levels(l)% var, &
                 Gamma_u         => THIS% Gamma_u% levels(l)% var &
      )


        !DO iz= 1, THIS% ngrid_z, 1
        !  DO iy= 1, THIS% ngrid_y, 1
        !    DO ix= 1, THIS% ngrid_x, 1
        !      abs_grid( 1, ix, iy, iz )= ABS( THIS% grid( 1, ix, iy, iz ) )
        !      abs_grid( 2, ix, iy, iz )= ABS( THIS% grid( 2, ix, iy, iz ) )
        !      abs_grid( 3, ix, iy, iz )= ABS( THIS% grid( 3, ix, iy, iz ) )
        !    ENDDO
        !  ENDDO
        !ENDDO

        min_abs_y= HUGE(1.0D0)
        min_abs_z= HUGE(1.0D0)
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

          DO j= 1, THIS% get_ngrid_y(l), 1

            DO i= 1, THIS% get_ngrid_x(l), 1

              IF( THIS% export_form_xy .AND. &
                  ( THIS% coords% levels(l)% var( i, j, k, jz ) /= &
                    THIS% coords% levels(l)% var( min_ix_z, min_iy_z, &
                                                  min_iz_z, jz ) ) )THEN
                CYCLE
              ENDIF
              IF( THIS% export_form_x .AND. &
                  ( THIS% coords% levels(l)% var( i, j, k, jz ) /= &
                    THIS% coords% levels(l)% var( min_ix_z, min_iy_z, &
                                                  min_iz_z, jz ) &
                    .OR. &
                    THIS% coords% levels(l)% var( i, j, k, jy ) /= &
                    THIS% coords% levels(l)% var( min_ix_y, min_iy_y, &
                                                  min_iz_y, jy ) ) )THEN
                CYCLE
              ENDIF

              WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * )&
                l, &
                THIS% coords% levels(l)% var( i, j, k, jx ), &
                THIS% coords% levels(l)% var( i, j, k, jy ), &
                THIS% coords% levels(l)% var( i, j, k, jz ), &
                lapse( i, j, k ), &
                shift_u( i, j, k, jx ), &
                shift_u( i, j, k, jy ), &
                shift_u( i, j, k, jz ), &
                phi( i, j, k ), &
                trK( i, j, k ), &
                g_BSSN3_ll( i, j, k, jxx ), &
                g_BSSN3_ll( i, j, k, jxy ), &
                g_BSSN3_ll( i, j, k, jxz ), &
                g_BSSN3_ll( i, j, k, jyy ), &
                g_BSSN3_ll( i, j, k, jyz ), &
                g_BSSN3_ll( i, j, k, jzz ), &
                A_BSSN3_ll( i, j, k, jxx ), &
                A_BSSN3_ll( i, j, k, jxy ), &
                A_BSSN3_ll( i, j, k, jxz ), &
                A_BSSN3_ll( i, j, k, jyy ), &
                A_BSSN3_ll( i, j, k, jyz ), &
                A_BSSN3_ll( i, j, k, jzz ), &
                Gamma_u( i, j, k, jx ), &
                Gamma_u( i, j, k, jy ), &
                Gamma_u( i, j, k, jz )

              IF( ios > 0 )THEN
                PRINT *, "...error when writing the arrays in ", &
                         TRIM(finalnamefile), ". The error message is", err_msg
                STOP
              ENDIF
              !CALL test_status( ios, err_msg, "...error when writing " &
              !                  // "the arrays in " // TRIM(finalnamefile) )

            ENDDO
          ENDDO
        ENDDO
      END ASSOCIATE
    ENDDO

    CLOSE( UNIT= 20 )

    PRINT *, " * LORENE BSSN ID on the gravity grid saved to formatted " &
             // "file ", TRIM(finalnamefile)

    PRINT *, "** Subroutine print_formatted_lorene_id_BSSN_variables " &
             // "executed."
    PRINT *

  END PROCEDURE print_formatted_lorene_id_bssn_variables


END SUBMODULE bssn_variables
