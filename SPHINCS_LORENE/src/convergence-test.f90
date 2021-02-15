! File:         convergence_test.f90
! Author:       Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

PROGRAM convergence_test

  !*****************************************************
  !                                                    *
  ! Make a convergence test to check the validity of   *
  ! the code SPHINCS_LORENE.                           *
  !                                                    *
  ! FT 8.12.2020                                       *
  !                                                    *
  !*****************************************************

  USE sphincs_lorene

  IMPLICIT NONE

  ! Maximum length for strings, and for the number of imported binaries
  INTEGER, PARAMETER:: max_length= 50
  ! Maximum number of binary systems
  INTEGER, PARAMETER:: max_n_bns= 50
  ! Maximum number of particle distributions
  INTEGER, PARAMETER:: max_n_parts= 250
  ! Loop limits for BSSN objects (for debugging; 3 is for production)
  INTEGER, PARAMETER:: min_bssn= 1
  INTEGER, PARAMETER:: max_bssn= 3
  ! Number of binary systems of neutron stars (BNS) to import
  INTEGER:: n_bns
  ! Export the constraints every constraints_step-th step
  INTEGER:: constraints_step

  ! Matrix storing the information on how to place particles for each bns
  ! object. Row i contains information about the i^th bns object.
  INTEGER, PARAMETER:: test_int= - 112
  INTEGER, DIMENSION( max_n_bns, max_n_parts ):: placer= test_int

  ! Grid spacing for the first BSSN object; the other two will have
  ! original_dx/2 and original_dx/4 as grid spacings
  DOUBLE PRECISION:: original_dx
  ! Rational ratio between the large grid spacing and the medium one,
  ! equal to the ratio between the medium grid spacing nd the small one
  DOUBLE PRECISION:: numerator_ratio_dx
  DOUBLE PRECISION:: denominator_ratio_dx
  DOUBLE PRECISION:: ratio_dx

  ! Strings storing different names for output files
  CHARACTER( LEN= 500 ):: namefile_parts, namefile_sph
  CHARACTER( LEN= 500 ):: namefile_bssn, name_logfile
  ! Array of strings storing the names of the LORENE BNS ID binary files
  CHARACTER( LEN= max_length ), DIMENSION( max_length ):: filenames= "0"
  ! String storing the local path to the directory where the
  ! LORENE BNS ID files are stored
  CHARACTER( LEN= max_length ):: common_path

  LOGICAL, PARAMETER:: debug= .FALSE.
  ! Logical variable to check if files exist
  LOGICAL:: exist
  ! Logical variables to steer the execution
  LOGICAL:: export_bin, export_form, export_form_xy, export_form_x, &
            compute_constraints, export_constraints_xy, &
            export_constraints_x, export_constraints, &
            export_constraints_details, compute_parts_constraints, &
            one_lapse, zero_shift

  TYPE( timer ):: execution_timer

  ! Declaration of the bns (binary neutron star) object containing the
  ! LORENE ID
  TYPE( bns ):: binary
  ! Declaration of the particles object containing the particle distribution
  TYPE( particles ):: particles_dist
  ! Declaration of the 3-component array storing the 3 bssn_id objects,
  ! containing the BSSN variables on the gravity grid
  TYPE( bssn_id ), DIMENSION(3):: bssn_forms

  ! Namelist containing parameters read from lorene_bns_id_parameters.par
  ! by the SUBROUTINE read_bns_id_parameters of this PROGRAM
  NAMELIST /bns_parameters/ n_bns, common_path, filenames, placer, &
                            export_bin, export_form, export_form_xy, &
                            export_form_x, export_constraints_xy, &
                            export_constraints_x, compute_constraints, &
                            export_constraints, export_constraints_details, &
                            constraints_step, compute_parts_constraints, &
                            numerator_ratio_dx, denominator_ratio_dx, &
                            one_lapse, zero_shift, show_progress

  !---------------------------!
  !--  End of declarations  --!
  !---------------------------!

  ! Conversions of some polytropic constants from LORENE units to
  ! SPHINCS units, and vice versa
  !gamma= 2
  !PRINT *, 0.0332278*k_lorene2hydrobase( gamma )
  !PRINT *
  !gamma= 2.75
  !PRINT *, 30000.0D0/k_lorene2hydrobase( gamma )
  !gamma= 2.75
  !PRINT *, 0.0332278*k_lorene2hydrobase( gamma )
  !STOP

  CALL DATE_AND_TIME( date, time, zone, values )
  run_id= date // "-" // time
  !PRINT *, run_id
  !STOP

  execution_timer= timer( "execution_timer" )
  CALL execution_timer% start_timer()

  CALL read_bns_id_parameters()

  ratio_dx= numerator_ratio_dx/denominator_ratio_dx

  ! Check that ratio_dx > 1
  IF( ratio_dx <= 1.0D0 )THEN
    PRINT *, "** ERROR! numerator_ratio_dx has to be larger than ", &
             "denominator_ratio_dx. The current values are ", &
             numerator_ratio_dx, " and ", denominator_ratio_dx, &
             ", respectively."
    PRINT *
    STOP
  ENDIF

  !
  !-- Construct the bns object from the LORENE binary file
  !
  binary= bns( TRIM(common_path)//"/"//TRIM(filenames( 1 )) )
  ! Set the variables to decide on using the geodesic gauge or not
  ! (lapse=1, shift=0)
  binary% one_lapse = one_lapse
  binary% zero_shift= zero_shift

  !
  !-- Construct the particles object from the bns object
  !
  IF( compute_parts_constraints )THEN
    PRINT *, "===================================================" &
             // "==============="
    PRINT *, " Placing particles"
    PRINT *, "===================================================" &
             // "==============="
    PRINT *
    particles_dist= particles( binary, placer( 1, 1 ) )

    !
    !-- Compute the SPH variables
    !
    PRINT *, "===================================================" &
             // "====================="
    PRINT *, " Computing SPH variables "
    PRINT *, "===================================================" &
             // "====================="
    PRINT *
    WRITE( namefile_parts, "(A1,I1,A1,I1,A1)" ) &
                                "l", &
                                1, "-", 1, "."
    particles_dist% export_bin= export_bin
    CALL particles_dist% compute_and_export_SPH_variables( namefile_parts )

    !
    !-- Print the particle initial data to a formatted file
    !
    IF( export_form )THEN
      WRITE( namefile_parts, "(A34)" ) &
                             "lorene-bns-id-particles-form_1.dat"
      particles_dist% export_form_xy= export_form_xy
      particles_dist% export_form_x = export_form_x
      CALL particles_dist% print_formatted_lorene_id_particles( namefile_parts )
    ENDIF
  ENDIF

  !
  !-- Construct the bssn_id objects from the bns object
  !
  construct_bssn_loop: DO itr3 = min_bssn, max_bssn, 1

    PRINT *, "===================================================" &
             // "==============="
    PRINT *, " Setting up BSSN object ", itr3
    PRINT *, "===================================================" &
             // "==============="
    PRINT *

    IF( itr3 == 1 )THEN

      bssn_forms( itr3 )= BSSN_id( binary )
      original_dx= bssn_forms( itr3 )% get_x_spacing()

    ELSE
      IF( itr3 == min_bssn )THEN

        bssn_forms( 1 )= BSSN_id( binary )
        original_dx= bssn_forms( 1 )% get_x_spacing()

      ENDIF
      bssn_forms( itr3 )= BSSN_id( binary, &
                                    original_dx/( ratio_dx**( itr3 - 1 ) ), &
                                    original_dx/( ratio_dx**( itr3 - 1 ) ), &
                                    original_dx/( ratio_dx**( itr3 - 1 ) ) )

      IF( bssn_forms( itr3 )% get_x_spacing() /= &
          original_dx/( ratio_dx**( itr3 - 1 ) ) )THEN
        PRINT *, " ** ERROR! The grid spacing #", itr3, ",", &
                 bssn_forms( itr3 )% get_x_spacing(), &
                 " is not equal to dx/4, with dx=", &
                 original_dx
        STOP
      ENDIF
    ENDIF

    PRINT *, "** The grid spacing is dx=", bssn_forms( itr3 )% get_x_spacing()
    PRINT *, "** The number of grid points for dx is:", &
                              bssn_forms( itr3 )% get_ngrid_x(), "**3"
    PRINT *

  ENDDO construct_bssn_loop

  !
  !-- Compute the BSSN variables
  !
  compute_export_bssn_loop: DO itr3 = min_bssn, max_bssn, 1
    PRINT *, "===================================================" &
             // "==============="
    PRINT *, " Computing BSSN variables for BSSN formulation", itr3
    PRINT *, "===================================================" &
             // "==============="
    PRINT *
    WRITE( namefile_bssn, "(A6,I1,A4)" ) "BSSN_l", itr3, ".bin"
    bssn_forms( itr3 )% export_bin= export_bin
    bssn_forms( itr3 )% export_form_xy= export_form_xy
    bssn_forms( itr3 )% export_form_x = export_form_x
    CALL bssn_forms( itr3 )% &
                        compute_and_export_3p1_variables( namefile_bssn )
  ENDDO compute_export_bssn_loop

  !
  !-- Print the BSSN initial data to a formatted file
  !
  IF( export_form )THEN
    export_bssn_loop: DO itr3 = min_bssn, max_bssn, 1
      WRITE( namefile_bssn, "(A24,I1,A4)" ) &
                            "lorene-bns-id-bssn-form_", itr3, ".dat"
      bssn_forms( itr3 )% export_form_xy= export_form_xy
      bssn_forms( itr3 )% export_form_x = export_form_x
      CALL bssn_forms( itr3 )% &
                  print_formatted_lorene_id_3p1_variables( namefile_bssn )
    ENDDO export_bssn_loop
  ENDIF

  !
  !-- Compute the BSSN constraints
  !
  compute_export_bssn_constraints_loop: DO itr3 = min_bssn, max_bssn, 1

    bssn_forms( itr3 )% cons_step= constraints_step
    bssn_forms( itr3 )% export_constraints= export_constraints
    bssn_forms( itr3 )% export_constraints_details= &
                        export_constraints_details
    bssn_forms( itr3 )% export_constraints_xy= export_constraints_xy
    bssn_forms( itr3 )% export_constraints_x = export_constraints_x

    IF( compute_constraints )THEN

      PRINT *, "===================================================" &
               // "==============="
      PRINT *, " Computing BSSN constraints for BSSN formulation", itr3
      PRINT *, "===================================================" &
               // "==============="
      PRINT *

      WRITE( namefile_bssn, "(A17,I1,A4)" ) "bssn-constraints-", itr3, ".dat"
      WRITE( name_logfile, "(A28,I1,A4)" ) &
                          "bssn-constraints-statistics-", itr3, ".log"

      CALL bssn_forms( itr3 )% &
                  compute_and_export_3p1_constraints( binary, &
                                                      namefile_bssn, &
                                                      name_logfile )

    ENDIF
    IF( compute_parts_constraints )THEN

      PRINT *, "===================================================" &
               // "==============="
      PRINT *, " Computing BSSN constraints with particle data for BSSN", &
               " formulation", itr3
      PRINT *, "===================================================" &
               // "==============="
      PRINT *

      WRITE( namefile_bssn, "(A23,I1,A4)" ) "bssn-constraints-parts-", &
                                           itr3, ".dat"
      WRITE( namefile_sph, "(A12,I1,A4)" ) "sph-density-", itr3, ".dat"
      WRITE( name_logfile, "(A34,I1,A4)" ) &
                           "bssn-constraints-parts-statistics-", itr3, ".log"

      CALL bssn_forms( itr3 )% &
                  compute_and_export_3p1_constraints( particles_dist, &
                                                      namefile_bssn, &
                                                      namefile_sph, &
                                                      name_logfile )

    ENDIF

  ENDDO compute_export_bssn_constraints_loop

  !
  !-- The BSSN formulations on grids with different resolutions are ready
  !-- to be used in a convergence test
  !
  IF( debug )THEN
    PRINT *, bssn_forms( 1 )% get_grid_point( 1 + 1, 1 + 2, 1 + 1 )
    PRINT *, bssn_forms( 2 )% get_grid_point( 1 + 2, 1 + 4, 1 + 2 )
    PRINT *, bssn_forms( 3 )% get_grid_point( 1 + 4, 1 + 8, 1 + 4 )
    PRINT *

    PRINT *, bssn_forms( 1 )% get_grid_point( 1 + 2, 1 + 1, 1 + 2 )
    PRINT *, bssn_forms( 2 )% get_grid_point( 1 + 4, 1 + 2, 1 + 4 )
    PRINT *, bssn_forms( 3 )% get_grid_point( 1 + 8, 1 + 4, 1 + 8 )
    PRINT *

    PRINT *, ABS(bssn_forms( 1 )% get_HC( 1 + 30,   1 + 25,   1 + 17 ))
    PRINT *, ABS(bssn_forms( 2 )% get_HC( 1 + 2*30, 1 + 2*25, 1 + 2*17 ))
    PRINT *, ABS(bssn_forms( 3 )% get_HC( 1 + 4*30, 1 + 4*25, 1 + 4*17 ))
    PRINT *
    PRINT *, ABS(bssn_forms( 1 )% get_HC( 1 + 30,   1 + 25,   1 + 17 )) &
             - ABS(bssn_forms( 2 )% get_HC( 1 + 2*30, 1 + 2*25, 1 + 2*17 ))
    PRINT *, ABS(bssn_forms( 2 )% get_HC( 1 + 2*30, 1 + 2*25, 1 + 2*17 )) &
             - ABS(bssn_forms( 3 )% get_HC( 1 + 4*30, 1 + 4*25, 1 + 4*17 ))
    STOP
  ENDIF

  !
  !-- Perform the convergence test with the appropriate constraints
  !
  IF( compute_constraints )THEN
    PRINT *, "** Performing convergence test with constraints computed ", &
             "without particle data."
    PRINT *
    CALL cauchy_convergence_test_unknown( bssn_forms(1), bssn_forms(2), &
                                          bssn_forms(3), 1 )
    CALL cauchy_convergence_test_known( bssn_forms(2), bssn_forms(3), 1 )
  ENDIF
  IF( compute_parts_constraints )THEN
    PRINT *, "** Performing convergence test with constraints computed ", &
             "with particle data."
    PRINT *
    CALL cauchy_convergence_test_unknown( bssn_forms(1), bssn_forms(2), &
                                          bssn_forms(3), 2 )
    CALL cauchy_convergence_test_known( bssn_forms(2), bssn_forms(3), 2 )
  ENDIF

  CALL execution_timer% stop_timer()

  !
  !-- Print the timers
  !
  PRINT *, "** Timing."
  PRINT *
  PRINT *
  PRINT *, " * BSSN formulation with spacing:", bssn_forms( 1 )% get_x_spacing()
  !original_dx
  CALL bssn_forms( 1 )% grid_timer% print_timer( 2 )
  CALL bssn_forms( 1 )% importer_timer% print_timer( 2 )
  CALL bssn_forms( 1 )% bssn_computer_timer% print_timer( 2 )
  PRINT *
  PRINT *, " * BSSN formulation with spacing:", bssn_forms( 2 )% get_x_spacing()
  !original_dx/2
  CALL bssn_forms( 2 )% grid_timer% print_timer( 2 )
  CALL bssn_forms( 2 )% importer_timer% print_timer( 2 )
  CALL bssn_forms( 2 )% bssn_computer_timer% print_timer( 2 )
  PRINT *
  PRINT *, " * BSSN formulation with spacing:", bssn_forms( 3 )% get_x_spacing()
  !original_dx/4
  CALL bssn_forms( 3 )% grid_timer% print_timer( 2 )
  CALL bssn_forms( 3 )% importer_timer% print_timer( 2 )
  CALL bssn_forms( 3 )% bssn_computer_timer% print_timer( 2 )
  PRINT *
  PRINT *, " * Total:"
  CALL execution_timer% print_timer( 2 )
  PRINT *

  !
  !-- Destruct the LORENE Bin_NS object by hand, since the pointer to it is
  !-- global (because it is bound to C++) and cannot be nullified by the
  !-- destructor of bns. In case of multiple bns objects, this would lead
  !-- to problems...
  !-- TODO: fix this
  !
  CALL binary% destruct_binary()


  CONTAINS


  SUBROUTINE cauchy_convergence_test_known( formul_dx, formul_dx2, &
                                            use_constraints )

    IMPLICIT NONE

    CLASS(formul_3p1), INTENT( IN OUT ):: formul_dx, formul_dx2
    INTEGER:: use_constraints

    INTEGER:: ix, iy, iz, nx, ny, nz, unit_cauchy_ct, unit_cauchy_parts_ct, &
              min_ix_y, min_iy_y, min_iz_y, &
              min_ix_z, min_iy_z, min_iz_z

    DOUBLE PRECISION:: min_abs_y, min_abs_z
    DOUBLE PRECISION, DIMENSION( :, :, :, : ), ALLOCATABLE:: abs_grid

    DOUBLE PRECISION, PARAMETER:: tiny_real= 1D-30
    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: grid_dx
    DOUBLE PRECISION, DIMENSION(3):: point_dx2
    DOUBLE PRECISION, DIMENSION(formul_dx% ngrid_x, &
                                formul_dx% ngrid_y, &
                                formul_dx% ngrid_z):: convergence_factor

    CHARACTER( LEN=: ), ALLOCATABLE:: name_cauchy_ct, name_cauchy_parts_ct

    LOGICAL:: exist

    nx= formul_dx% get_ngrid_x()
    ny= formul_dx% get_ngrid_y()
    nz= formul_dx% get_ngrid_z()

    nx= FLOOR( DBLE( nx - 1 )/denominator_ratio_dx ) + 1
    ny= FLOOR( DBLE( ny - 1 )/denominator_ratio_dx ) + 1
    nz= FLOOR( DBLE( nz - 1 )/denominator_ratio_dx ) + 1

    ALLOCATE( grid_dx( 3, nx, ny, nz ) )
    ALLOCATE( abs_grid( 3, nx, ny, nz ) )

    PRINT *, "** Computing convergence factor..."

    choose_constraints: SELECT CASE( use_constraints )

    CASE(1)

      shared_grid_loops1: DO iz= 0, nz - 1, 1
        DO iy= 0, ny - 1, 1
          DO ix= 0, nx - 1, 1

            grid_dx( :, 1 + ix, 1 + iy, 1 + iz ) = &
                        formul_dx%  get_grid_point(  &
                                      1 + INT(denominator_ratio_dx)*ix, &
                                      1 + INT(denominator_ratio_dx)*iy, &
                                      1 + INT(denominator_ratio_dx)*iz   )
                        point_dx2= formul_dx2% get_grid_point( &
                                      1 + INT(numerator_ratio_dx)*ix, &
                                      1 + INT(numerator_ratio_dx)*iy, &
                                      1 + INT(numerator_ratio_dx)*iz )

            IF(ABS(grid_dx( 1, 1 + ix, 1 + iy, 1 + iz )-point_dx2(1)) > 1D-10 &
          .OR. ABS(grid_dx( 2, 1 + ix, 1 + iy, 1 + iz )-point_dx2(2)) > 1D-10 &
          .OR. ABS(grid_dx( 3, 1 + ix, 1 + iy, 1 + iz )-point_dx2(3)) > 1D-10 &

           )THEN

              PRINT *, "**ERROR! The grid functions in the Cauchy ", &
                       "convergence test are not evaluated at the ", &
                       "same grid point at (ix,iy,iz)=(", &
                       ix, iy, iz, ")."
              PRINT *, grid_dx( 1, 1 + ix, 1 + iy, 1 + iz ), point_dx2(1)
              PRINT *, grid_dx( 2, 1 + ix, 1 + iy, 1 + iz ), point_dx2(2)
              PRINT *, grid_dx( 3, 1 + ix, 1 + iy, 1 + iz ), point_dx2(3)
              PRINT *
              STOP

            ENDIF

            convergence_factor( 1 + ix, 1 + iy, 1 + iz )= &
             LOG( &
             ABS( &
             ( formul_dx%  get_HC( 1 + INT(denominator_ratio_dx)*ix, &
                                   1 + INT(denominator_ratio_dx)*iy, &
                                   1 + INT(denominator_ratio_dx)*iz ))&
            /( formul_dx2% get_HC( 1 + INT(numerator_ratio_dx)*ix, &
                                   1 + INT(numerator_ratio_dx)*iy, &
                                   1 + INT(numerator_ratio_dx)*iz ) &
             + 0*tiny_real ) &
             ) )/LOG(ratio_dx)

          ENDDO
        ENDDO
      ENDDO shared_grid_loops1
      PRINT *, " * Convergence factor computed."
      PRINT *

      unit_cauchy_ct= 3109
      name_cauchy_ct= "cauchy_convergence_test_known.dat"

      INQUIRE( FILE= TRIM(name_cauchy_ct), EXIST= exist )

      IF( exist )THEN
        OPEN( UNIT= unit_cauchy_ct, FILE= TRIM(name_cauchy_ct), &
              STATUS= "REPLACE", FORM= "FORMATTED", &
              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
      ELSE
        OPEN( UNIT= unit_cauchy_ct, FILE= TRIM(name_cauchy_ct), &
              STATUS= "NEW", FORM= "FORMATTED", &
              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
      ENDIF
      CALL test_status( ios, err_msg, "...error when opening " &
               // TRIM(name_cauchy_ct) )

      WRITE( UNIT = unit_cauchy_ct, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id
      WRITE( UNIT= unit_cauchy_ct, IOSTAT = ios, &
             IOMSG = err_msg, FMT = * ) &
      "# Cauchy convergence test. "
      WRITE( UNIT = unit_cauchy_ct, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "# column:      1        2       3       4"
      WRITE( UNIT = unit_cauchy_ct, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "#      x [km]       y [km]       z [km]       " &
      //"convergence factor [pure number]"

      DO iz= 1, nx, 1
        DO iy= 1, ny, 1
          DO ix= 1, nz, 1
            abs_grid( 1, ix, iy, iz )= ABS( grid_dx( 1, ix, iy, iz ) )
            abs_grid( 2, ix, iy, iz )= ABS( grid_dx( 2, ix, iy, iz ) )
            abs_grid( 3, ix, iy, iz )= ABS( grid_dx( 3, ix, iy, iz ) )
          ENDDO
        ENDDO
      ENDDO

      min_abs_y= 1D+20
      min_abs_z= 1D+20
      DO iz= 1, nx, 1
        DO iy= 1, ny, 1
          DO ix= 1, nz, 1
            IF( ABS( grid_dx( 2, ix, iy, iz ) ) < min_abs_y )THEN
              min_abs_y= ABS( grid_dx( 2, ix, iy, iz ) )
              min_ix_y= ix
              min_iy_y= iy
              min_iz_y= iz
            ENDIF
            IF( ABS( grid_dx( 3, ix, iy, iz ) ) < min_abs_z )THEN
              min_abs_z= ABS( grid_dx( 3, ix, iy, iz ) )
              min_ix_z= ix
              min_iy_z= iy
              min_iz_z= iz
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DO iz= 1, nz, 1
        DO iy= 1, ny, 1
          DO ix= 1, nx, 1
            IF( .FALSE. .AND. export_constraints_xy .AND. &
                grid_dx( 3, ix, iy, iz ) /= &
                grid_dx( 3, min_ix_z, min_iy_z, min_iz_z ) )THEN
              CYCLE
            ENDIF
            IF( .FALSE. .AND. export_constraints_x .AND. &
                ( grid_dx( 3, ix, iy, iz ) /= &
                  grid_dx( 3, min_ix_z, min_iy_z, min_iz_z ) &
                  .OR. &
                  grid_dx( 2, ix, iy, iz ) /= &
                  grid_dx( 2, min_ix_y, min_iy_y, min_iz_y ) ) )THEN
              CYCLE
            ENDIF
            WRITE( UNIT = unit_cauchy_ct, IOSTAT = ios, &
                   IOMSG = err_msg, FMT = * )&
                grid_dx( 1, ix, iy, iz ), &
                grid_dx( 2, ix, iy, iz ), &
                grid_dx( 3, ix, iy, iz ), &
                convergence_factor( ix, iy, iz )
            CALL test_status( ios, err_msg, "...error in writing " &
                        // "the arrays in " // TRIM(name_cauchy_ct) )
          ENDDO
        ENDDO
      ENDDO

      CLOSE( UNIT= unit_cauchy_ct )

      PRINT *, " * Convergence factor exported to formatted file ", &
                                                      TRIM(name_cauchy_ct)
      PRINT *

    CASE(2)

      shared_grid_loops2: DO iz= 0, nz - 1, 1
        DO iy= 0, ny - 1, 1
          DO ix= 0, nx - 1, 1

            grid_dx( :, 1 + ix, 1 + iy, 1 + iz ) = &
                        formul_dx%  get_grid_point(  &
                                      1 + INT(denominator_ratio_dx)*ix, &
                                      1 + INT(denominator_ratio_dx)*iy, &
                                      1 + INT(denominator_ratio_dx)*iz   )
                        point_dx2= formul_dx2% get_grid_point( &
                                      1 + INT(numerator_ratio_dx)*ix, &
                                      1 + INT(numerator_ratio_dx)*iy, &
                                      1 + INT(numerator_ratio_dx)*iz )

            IF(ABS(grid_dx( 1, 1 + ix, 1 + iy, 1 + iz )-point_dx2(1)) > 1D-10 &
          .OR. ABS(grid_dx( 2, 1 + ix, 1 + iy, 1 + iz )-point_dx2(2)) > 1D-10 &
          .OR. ABS(grid_dx( 3, 1 + ix, 1 + iy, 1 + iz )-point_dx2(3)) > 1D-10 &

           )THEN

              PRINT *, "**ERROR! The grid functions in the Cauchy ", &
                       "convergence test are not evaluated at the ", &
                       "same grid point at (ix,iy,iz)=(", &
                       ix, iy, iz, ")."
              PRINT *, grid_dx( 1, 1 + ix, 1 + iy, 1 + iz ), point_dx2(1)
              PRINT *, grid_dx( 2, 1 + ix, 1 + iy, 1 + iz ), point_dx2(2)
              PRINT *, grid_dx( 3, 1 + ix, 1 + iy, 1 + iz ), point_dx2(3)
              PRINT *
              STOP

            ENDIF

            convergence_factor( 1 + ix, 1 + iy, 1 + iz )= &
             LOG( ABS( &
               formul_dx%  get_HC_parts( 1 + INT(denominator_ratio_dx)*ix, &
                                         1 + INT(denominator_ratio_dx)*iy, &
                                         1 + INT(denominator_ratio_dx)*iz )&
            /( formul_dx2% get_HC_parts( 1 + INT(numerator_ratio_dx)*ix, &
                                         1 + INT(numerator_ratio_dx)*iy, &
                                         1 + INT(numerator_ratio_dx)*iz )  &
             + 0*tiny_real ) &
             ) )/LOG(ratio_dx)

          ENDDO
        ENDDO
      ENDDO shared_grid_loops2
      PRINT *, " * Convergence factor computed."
      PRINT *

      unit_cauchy_parts_ct= 3111
      name_cauchy_parts_ct= "cauchy_convergence_test_known_parts.dat"

      INQUIRE( FILE= TRIM(name_cauchy_parts_ct), EXIST= exist )

      IF( exist )THEN
        OPEN( UNIT= unit_cauchy_parts_ct, FILE= TRIM(name_cauchy_parts_ct), &
              STATUS= "REPLACE", FORM= "FORMATTED", &
              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
      ELSE
        OPEN( UNIT= unit_cauchy_parts_ct, FILE= TRIM(name_cauchy_parts_ct), &
              STATUS= "NEW", FORM= "FORMATTED", &
              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
      ENDIF
      CALL test_status( ios, err_msg, "...error when opening " &
               // TRIM(name_cauchy_parts_ct) )

      WRITE( UNIT = unit_cauchy_parts_ct, IOSTAT = ios, IOMSG = err_msg, &
             FMT = * ) &
      "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id
      WRITE( UNIT= unit_cauchy_parts_ct, IOSTAT = ios, &
             IOMSG = err_msg, FMT = * ) &
      "# Cauchy convergence test. "
      WRITE( UNIT = unit_cauchy_parts_ct, IOSTAT = ios, IOMSG = err_msg, &
             FMT = * ) &
      "# column:      1        2       3       4"
      WRITE( UNIT = unit_cauchy_parts_ct, IOSTAT = ios, IOMSG = err_msg, &
             FMT = * ) &
      "#      x [km]       y [km]       z [km]       " &
      //"convergence factor [pure number]"

      DO iz= 1, nx, 1
        DO iy= 1, ny, 1
          DO ix= 1, nz, 1
            abs_grid( 1, ix, iy, iz )= ABS( grid_dx( 1, ix, iy, iz ) )
            abs_grid( 2, ix, iy, iz )= ABS( grid_dx( 2, ix, iy, iz ) )
            abs_grid( 3, ix, iy, iz )= ABS( grid_dx( 3, ix, iy, iz ) )
          ENDDO
        ENDDO
      ENDDO

      min_abs_y= 1D+20
      min_abs_z= 1D+20
      DO iz= 1, nx, 1
        DO iy= 1, ny, 1
          DO ix= 1, nz, 1
            IF( ABS( grid_dx( 2, ix, iy, iz ) ) < min_abs_y )THEN
              min_abs_y= ABS( grid_dx( 2, ix, iy, iz ) )
              min_ix_y= ix
              min_iy_y= iy
              min_iz_y= iz
            ENDIF
            IF( ABS( grid_dx( 3, ix, iy, iz ) ) < min_abs_z )THEN
              min_abs_z= ABS( grid_dx( 3, ix, iy, iz ) )
              min_ix_z= ix
              min_iy_z= iy
              min_iz_z= iz
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DO iz= 1, nz, 1
        DO iy= 1, ny, 1
          DO ix= 1, nx, 1
            IF( .FALSE. .AND. export_constraints_xy .AND. &
                grid_dx( 3, ix, iy, iz ) /= &
                grid_dx( 3, min_ix_z, min_iy_z, min_iz_z ) )THEN
              CYCLE
            ENDIF
            IF( .FALSE. .AND. export_constraints_x .AND. &
                ( grid_dx( 3, ix, iy, iz ) /= &
                  grid_dx( 3, min_ix_z, min_iy_z, min_iz_z ) &
                  .OR. &
                  grid_dx( 2, ix, iy, iz ) /= &
                  grid_dx( 2, min_ix_y, min_iy_y, min_iz_y ) ) )THEN
              CYCLE
            ENDIF
            WRITE( UNIT = unit_cauchy_parts_ct, IOSTAT = ios, &
                   IOMSG = err_msg, FMT = * )&
                grid_dx( 1, ix, iy, iz ), &
                grid_dx( 2, ix, iy, iz ), &
                grid_dx( 3, ix, iy, iz ), &
                convergence_factor( ix, iy, iz )
            CALL test_status( ios, err_msg, "...error in writing " &
                        // "the arrays in " // TRIM(name_cauchy_parts_ct) )
          ENDDO
        ENDDO
      ENDDO

      CLOSE( UNIT= unit_cauchy_parts_ct )

      PRINT *, " * Convergence factor exported to formatted file ", &
                                                  TRIM(name_cauchy_parts_ct)
      PRINT *

    CASE DEFAULT

      PRINT *, "** There is no well defined algorithm " &
               // "corresponding to the number", use_constraints
      PRINT *, " * Please set use_constraints to 1 or 2."
      STOP

    END SELECT choose_constraints

    DEALLOCATE( grid_dx )

  END SUBROUTINE cauchy_convergence_test_known


  SUBROUTINE cauchy_convergence_test_unknown( formul_dx, formul_dx2, &
                                              formul_dx4, use_constraints )

    IMPLICIT NONE

    CLASS(formul_3p1), INTENT( IN OUT ):: formul_dx, formul_dx2, formul_dx4
    INTEGER:: use_constraints

    INTEGER:: ix, iy, iz, nx, ny, nz, unit_cauchy_ct, unit_cauchy_parts_ct, &
              min_ix_y, min_iy_y, min_iz_y, &
              min_ix_z, min_iy_z, min_iz_z

    DOUBLE PRECISION:: min_abs_y, min_abs_z
    DOUBLE PRECISION, DIMENSION( :, :, :, : ), ALLOCATABLE:: abs_grid

    DOUBLE PRECISION, PARAMETER:: tiny_real= 1D-30
    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: grid_dx
    DOUBLE PRECISION, DIMENSION(3):: point_dx2
    DOUBLE PRECISION, DIMENSION(3):: point_dx4
    DOUBLE PRECISION, DIMENSION(formul_dx% ngrid_x, &
                                formul_dx% ngrid_y, &
                                formul_dx% ngrid_z):: convergence_factor

    CHARACTER( LEN=: ), ALLOCATABLE:: name_cauchy_ct, name_cauchy_parts_ct

    LOGICAL:: exist

    nx= formul_dx% get_ngrid_x()
    ny= formul_dx% get_ngrid_y()
    nz= formul_dx% get_ngrid_z()

    nx= FLOOR( DBLE( nx - 1 )/denominator_ratio_dx**2 ) + 1
    ny= FLOOR( DBLE( ny - 1 )/denominator_ratio_dx**2 ) + 1
    nz= FLOOR( DBLE( nz - 1 )/denominator_ratio_dx**2 ) + 1

    ALLOCATE( grid_dx( 3, nx, ny, nz ) )
    ALLOCATE( abs_grid( 3, nx, ny, nz ) )

    PRINT *, "** Computing convergence factor..."

    choose_constraints: SELECT CASE( use_constraints )

    CASE(1)

      shared_grid_loops3: DO iz= 0, nz - 1, 1
        DO iy= 0, ny - 1, 1
          DO ix= 0, nx - 1, 1

           ! grid_dx( :, 1 + ix, 1 + iy, 1 + iz ) = &
           !            formul_dx%  get_grid_point(  &
           !                             1 + ix,   1 + iy,   1 + iz   )
           ! point_dx2= formul_dx2% get_grid_point( &
           !                             1 + 2*ix, 1 + 2*iy, 1 + 2*iz )
           ! point_dx4= formul_dx4% get_grid_point( &
           !                             1 + 4*ix, 1 + 4*iy, 1 + 4*iz )
           !
           ! IF( grid_dx( 1, 1 + ix, 1 + iy, 1 + iz ) /= point_dx2(1) &
           !.OR. grid_dx( 1, 1 + ix, 1 + iy, 1 + iz ) /= point_dx4(1) &
           !.OR. grid_dx( 2, 1 + ix, 1 + iy, 1 + iz ) /= point_dx2(2) &
           !.OR. grid_dx( 2, 1 + ix, 1 + iy, 1 + iz ) /= point_dx4(2) &
           !.OR. grid_dx( 3, 1 + ix, 1 + iy, 1 + iz ) /= point_dx2(3) &
           !.OR. grid_dx( 3, 1 + ix, 1 + iy, 1 + iz ) /= point_dx4(3) &

            grid_dx( :, 1 + ix, 1 + iy, 1 + iz ) = &
                       formul_dx%  get_grid_point(  &
                                      1 + INT(denominator_ratio_dx**2)*ix, &
                                      1 + INT(denominator_ratio_dx**2)*iy, &
                                      1 + INT(denominator_ratio_dx**2)*iz   )
            point_dx2= formul_dx2% get_grid_point( &
                                      1 + INT(numerator_ratio_dx* &
                                          denominator_ratio_dx)*ix, &
                                      1 + INT(numerator_ratio_dx* &
                                          denominator_ratio_dx)*iy, &
                                      1 + INT(numerator_ratio_dx* &
                                          denominator_ratio_dx)*iz )
            point_dx4= formul_dx4% get_grid_point( &
                                      1 + INT(numerator_ratio_dx**2)*ix, &
                                      1 + INT(numerator_ratio_dx**2)*iy, &
                                      1 + INT(numerator_ratio_dx**2)*iz )

            IF(ABS(grid_dx( 1, 1 + ix, 1 + iy, 1 + iz )-point_dx2(1)) > 1D-10 &
          .OR. ABS(grid_dx( 1, 1 + ix, 1 + iy, 1 + iz )-point_dx4(1)) > 1D-10 &
          .OR. ABS(grid_dx( 2, 1 + ix, 1 + iy, 1 + iz )-point_dx2(2)) > 1D-10 &
          .OR. ABS(grid_dx( 2, 1 + ix, 1 + iy, 1 + iz )-point_dx4(2)) > 1D-10 &
          .OR. ABS(grid_dx( 3, 1 + ix, 1 + iy, 1 + iz )-point_dx2(3)) > 1D-10 &
          .OR. ABS(grid_dx( 3, 1 + ix, 1 + iy, 1 + iz )-point_dx4(3)) > 1D-10 &

           )THEN

              PRINT *, "**ERROR! The grid functions in the Cauchy ", &
                       "convergence test are not evaluated at the ", &
                       "same grid point at (ix,iy,iz)=(", &
                       ix, iy, iz, ")."
              PRINT *, grid_dx( 1, 1 + ix, 1 + iy, 1 + iz ), point_dx2(1), &
                       point_dx4(1)
              PRINT *, grid_dx( 2, 1 + ix, 1 + iy, 1 + iz ), point_dx2(2), &
                       point_dx4(2)
              PRINT *, grid_dx( 3, 1 + ix, 1 + iy, 1 + iz ), point_dx2(3), &
                       point_dx4(3)
              PRINT *
              STOP

            ENDIF

            !convergence_factor( 1 + ix, 1 + iy, 1 + iz )= &
            ! LOG( &
            ! ABS( &
            ! ( ABS(formul_dx%  get_HC( 1 + ix,   1 + iy,   1 + iz )) &
            ! - ABS(formul_dx2% get_HC( 1 + 2*ix, 1 + 2*iy, 1 + 2*iz )))&
            !/( ABS(formul_dx2% get_HC( 1 + 2*ix, 1 + 2*iy, 1 + 2*iz )) &
            ! - ABS(formul_dx4% get_HC( 1 + 4*ix, 1 + 4*iy, 1 + 4*iz ) )&
            ! + 0*tiny_real ) &
            ! ) &
            ! )/ln2!/LOG(2.0)

            convergence_factor( 1 + ix, 1 + iy, 1 + iz )= &
             LOG( &
             ABS( &
             ( ABS(formul_dx%  get_HC( 1 + INT(denominator_ratio_dx**2)*ix, &
                                       1 + INT(denominator_ratio_dx**2)*iy, &
                                       1 + INT(denominator_ratio_dx**2)*iz )) &
             - ABS(formul_dx2% get_HC( 1 + INT(numerator_ratio_dx* &
                                           denominator_ratio_dx)*ix, &
                                       1 + INT(numerator_ratio_dx* &
                                           denominator_ratio_dx)*iy, &
                                       1 + INT(numerator_ratio_dx* &
                                           denominator_ratio_dx)*iz )))&
            /( ABS(formul_dx2% get_HC( 1 + INT(numerator_ratio_dx* &
                                           denominator_ratio_dx)*ix, &
                                       1 + INT(numerator_ratio_dx* &
                                           denominator_ratio_dx)*iy, &
                                       1 + INT(numerator_ratio_dx* &
                                           denominator_ratio_dx)*iz )) &
             - ABS(formul_dx4% get_HC( 1 + INT(numerator_ratio_dx**2)*ix, &
                                       1 + INT(numerator_ratio_dx**2)*iy, &
                                       1 + INT(numerator_ratio_dx**2)*iz ) )&
             + 0*tiny_real ) &
             ) &
             )/LOG(ratio_dx)

          ENDDO
        ENDDO
      ENDDO shared_grid_loops3
      PRINT *, " * Convergence factor computed."
      PRINT *

      unit_cauchy_ct= 3108
      name_cauchy_ct= "cauchy_convergence_test_unknown.dat"

      INQUIRE( FILE= TRIM(name_cauchy_ct), EXIST= exist )

      IF( exist )THEN
        OPEN( UNIT= unit_cauchy_ct, FILE= TRIM(name_cauchy_ct), &
              STATUS= "REPLACE", FORM= "FORMATTED", &
              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
      ELSE
        OPEN( UNIT= unit_cauchy_ct, FILE= TRIM(name_cauchy_ct), &
              STATUS= "NEW", FORM= "FORMATTED", &
              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
      ENDIF
      CALL test_status( ios, err_msg, "...error when opening " &
               // TRIM(name_cauchy_ct) )

      WRITE( UNIT = unit_cauchy_ct, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id
      WRITE( UNIT= unit_cauchy_ct, IOSTAT = ios, &
             IOMSG = err_msg, FMT = * ) &
      "# Cauchy convergence test. "
      WRITE( UNIT = unit_cauchy_ct, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "# column:      1        2       3       4"
      WRITE( UNIT = unit_cauchy_ct, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      "#      x [km]       y [km]       z [km]       " &
      //"convergence factor [pure number]"

      DO iz= 1, nx, 1
        DO iy= 1, ny, 1
          DO ix= 1, nz, 1
            abs_grid( 1, ix, iy, iz )= ABS( grid_dx( 1, ix, iy, iz ) )
            abs_grid( 2, ix, iy, iz )= ABS( grid_dx( 2, ix, iy, iz ) )
            abs_grid( 3, ix, iy, iz )= ABS( grid_dx( 3, ix, iy, iz ) )
          ENDDO
        ENDDO
      ENDDO

      min_abs_y= 1D+20
      min_abs_z= 1D+20
      DO iz= 1, nx, 1
        DO iy= 1, ny, 1
          DO ix= 1, nz, 1
            IF( ABS( grid_dx( 2, ix, iy, iz ) ) < min_abs_y )THEN
              min_abs_y= ABS( grid_dx( 2, ix, iy, iz ) )
              min_ix_y= ix
              min_iy_y= iy
              min_iz_y= iz
            ENDIF
            IF( ABS( grid_dx( 3, ix, iy, iz ) ) < min_abs_z )THEN
              min_abs_z= ABS( grid_dx( 3, ix, iy, iz ) )
              min_ix_z= ix
              min_iy_z= iy
              min_iz_z= iz
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DO iz= 1, nz, 1
        DO iy= 1, ny, 1
          DO ix= 1, nx, 1
            IF( .FALSE. .AND. export_constraints_xy .AND. &
                grid_dx( 3, ix, iy, iz ) /= &
                grid_dx( 3, min_ix_z, min_iy_z, min_iz_z ) )THEN
              CYCLE
            ENDIF
            IF( .FALSE. .AND. export_constraints_x .AND. &
                ( grid_dx( 3, ix, iy, iz ) /= &
                  grid_dx( 3, min_ix_z, min_iy_z, min_iz_z ) &
                  .OR. &
                  grid_dx( 2, ix, iy, iz ) /= &
                  grid_dx( 2, min_ix_y, min_iy_y, min_iz_y ) ) )THEN
              CYCLE
            ENDIF
            WRITE( UNIT = unit_cauchy_ct, IOSTAT = ios, &
                   IOMSG = err_msg, FMT = * )&
                grid_dx( 1, ix, iy, iz ), &
                grid_dx( 2, ix, iy, iz ), &
                grid_dx( 3, ix, iy, iz ), &
                convergence_factor( ix, iy, iz )
            CALL test_status( ios, err_msg, "...error in writing " &
                        // "the arrays in " // TRIM(name_cauchy_ct) )
          ENDDO
        ENDDO
      ENDDO

      CLOSE( UNIT= unit_cauchy_ct )

      PRINT *, " * Convergence factor exported to formatted file ", &
                                                      TRIM(name_cauchy_ct)
      PRINT *

    CASE(2)

      shared_grid_loops4: DO iz= 0, nz - 1, 1
        DO iy= 0, ny - 1, 1
          DO ix= 0, nx - 1, 1

         ! grid_dx( :, 1 + ix, 1 + iy, 1 + iz ) = &
         !            formul_dx%  get_grid_point(  &
         !                             1 + ix,   1 + iy,   1 + iz   )
         ! point_dx2= formul_dx2% get_grid_point( &
         !                             1 + 2*ix, 1 + 2*iy, 1 + 2*iz )
         ! point_dx4= formul_dx4% get_grid_point( &
         !                             1 + 4*ix, 1 + 4*iy, 1 + 4*iz )
         !
         ! IF( grid_dx( 1, 1 + ix, 1 + iy, 1 + iz ) /= point_dx2(1) &
         !.OR. grid_dx( 1, 1 + ix, 1 + iy, 1 + iz ) /= point_dx4(1) &
         !.OR. grid_dx( 2, 1 + ix, 1 + iy, 1 + iz ) /= point_dx2(2) &
         !.OR. grid_dx( 2, 1 + ix, 1 + iy, 1 + iz ) /= point_dx4(2) &
         !.OR. grid_dx( 3, 1 + ix, 1 + iy, 1 + iz ) /= point_dx2(3) &
         !.OR. grid_dx( 3, 1 + ix, 1 + iy, 1 + iz ) /= point_dx4(3) &

          grid_dx( :, 1 + ix, 1 + iy, 1 + iz ) = &
                     formul_dx%  get_grid_point(  &
                                    1 + INT(denominator_ratio_dx**2)*ix, &
                                    1 + INT(denominator_ratio_dx**2)*iy, &
                                    1 + INT(denominator_ratio_dx**2)*iz   )
          point_dx2= formul_dx2% get_grid_point( &
                                    1 + INT(numerator_ratio_dx* &
                                        denominator_ratio_dx)*ix, &
                                    1 + INT(numerator_ratio_dx* &
                                        denominator_ratio_dx)*iy, &
                                    1 + INT(numerator_ratio_dx* &
                                        denominator_ratio_dx)*iz )
          point_dx4= formul_dx4% get_grid_point( &
                                    1 + INT(numerator_ratio_dx**2)*ix, &
                                    1 + INT(numerator_ratio_dx**2)*iy, &
                                    1 + INT(numerator_ratio_dx**2)*iz )

          IF(ABS(grid_dx( 1, 1 + ix, 1 + iy, 1 + iz )-point_dx2(1)) > 1D-10 &
        .OR. ABS(grid_dx( 1, 1 + ix, 1 + iy, 1 + iz )-point_dx4(1)) > 1D-10 &
        .OR. ABS(grid_dx( 2, 1 + ix, 1 + iy, 1 + iz )-point_dx2(2)) > 1D-10 &
        .OR. ABS(grid_dx( 2, 1 + ix, 1 + iy, 1 + iz )-point_dx4(2)) > 1D-10 &
        .OR. ABS(grid_dx( 3, 1 + ix, 1 + iy, 1 + iz )-point_dx2(3)) > 1D-10 &
        .OR. ABS(grid_dx( 3, 1 + ix, 1 + iy, 1 + iz )-point_dx4(3)) > 1D-10 &

         )THEN

            PRINT *, "**ERROR! The grid functions in the Cauchy ", &
                     "convergence test are not evaluated at the ", &
                     "same grid point at (ix,iy,iz)=(", &
                     ix, iy, iz, ")."
            PRINT *, grid_dx( 1, 1 + ix, 1 + iy, 1 + iz ), point_dx2(1), &
                     point_dx4(1)
            PRINT *, grid_dx( 2, 1 + ix, 1 + iy, 1 + iz ), point_dx2(2), &
                     point_dx4(2)
            PRINT *, grid_dx( 3, 1 + ix, 1 + iy, 1 + iz ), point_dx2(3), &
                     point_dx4(3)
            PRINT *
            STOP

          ENDIF

          convergence_factor( 1 + ix, 1 + iy, 1 + iz )= &
           LOG( &
           ABS( &
           ( ABS(formul_dx%  get_HC_parts( &
                                1 + INT(denominator_ratio_dx**2)*ix, &
                                1 + INT(denominator_ratio_dx**2)*iy, &
                                1 + INT(denominator_ratio_dx**2)*iz )) &
           - ABS(formul_dx2% get_HC_parts( &
                                1 + INT(numerator_ratio_dx* &
                                        denominator_ratio_dx)*ix, &
                                1 + INT(numerator_ratio_dx* &
                                        denominator_ratio_dx)*iy, &
                                1 + INT(numerator_ratio_dx* &
                                        denominator_ratio_dx)*iz )))&
          /( ABS(formul_dx2% get_HC_parts( &
                                1 + INT(numerator_ratio_dx* &
                                        denominator_ratio_dx)*ix, &
                                1 + INT(numerator_ratio_dx* &
                                        denominator_ratio_dx)*iy, &
                                1 + INT(numerator_ratio_dx* &
                                        denominator_ratio_dx)*iz )) &
           - ABS(formul_dx4% get_HC_parts( &
                                1 + INT(numerator_ratio_dx**2)*ix, &
                                1 + INT(numerator_ratio_dx**2)*iy, &
                                1 + INT(numerator_ratio_dx**2)*iz ) )&
           + 0*tiny_real ) &
           ) &
           )/LOG(ratio_dx)

          ENDDO
        ENDDO
      ENDDO shared_grid_loops4
      PRINT *, " * Convergence factor computed."
      PRINT *

      unit_cauchy_parts_ct= 3110
      name_cauchy_parts_ct= "cauchy_convergence_test_unknown_parts.dat"

      INQUIRE( FILE= TRIM(name_cauchy_parts_ct), EXIST= exist )

      IF( exist )THEN
        OPEN( UNIT= unit_cauchy_parts_ct, FILE= TRIM(name_cauchy_parts_ct), &
              STATUS= "REPLACE", FORM= "FORMATTED", &
              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
      ELSE
        OPEN( UNIT= unit_cauchy_parts_ct, FILE= TRIM(name_cauchy_parts_ct), &
              STATUS= "NEW", FORM= "FORMATTED", &
              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
      ENDIF
      CALL test_status( ios, err_msg, "...error when opening " &
               // TRIM(name_cauchy_parts_ct) )

      WRITE( UNIT = unit_cauchy_parts_ct, IOSTAT = ios, IOMSG = err_msg, &
             FMT = * ) &
      "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id
      WRITE( UNIT= unit_cauchy_parts_ct, IOSTAT = ios, &
             IOMSG = err_msg, FMT = * ) &
      "# Cauchy convergence test. "
      WRITE( UNIT = unit_cauchy_parts_ct, IOSTAT = ios, IOMSG = err_msg, &
             FMT = * ) &
      "# column:      1        2       3       4"
      WRITE( UNIT = unit_cauchy_parts_ct, IOSTAT = ios, IOMSG = err_msg, &
             FMT = * ) &
      "#      x [km]       y [km]       z [km]       " &
      //"convergence factor [pure number]"

      DO iz= 1, nx, 1
        DO iy= 1, ny, 1
          DO ix= 1, nz, 1
            abs_grid( 1, ix, iy, iz )= ABS( grid_dx( 1, ix, iy, iz ) )
            abs_grid( 2, ix, iy, iz )= ABS( grid_dx( 2, ix, iy, iz ) )
            abs_grid( 3, ix, iy, iz )= ABS( grid_dx( 3, ix, iy, iz ) )
          ENDDO
        ENDDO
      ENDDO

      min_abs_y= 1D+20
      min_abs_z= 1D+20
      DO iz= 1, nx, 1
        DO iy= 1, ny, 1
          DO ix= 1, nz, 1
            IF( ABS( grid_dx( 2, ix, iy, iz ) ) < min_abs_y )THEN
              min_abs_y= ABS( grid_dx( 2, ix, iy, iz ) )
              min_ix_y= ix
              min_iy_y= iy
              min_iz_y= iz
            ENDIF
            IF( ABS( grid_dx( 3, ix, iy, iz ) ) < min_abs_z )THEN
              min_abs_z= ABS( grid_dx( 3, ix, iy, iz ) )
              min_ix_z= ix
              min_iy_z= iy
              min_iz_z= iz
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DO iz= 1, nz, 1
        DO iy= 1, ny, 1
          DO ix= 1, nx, 1
            IF( .FALSE. .AND. export_constraints_xy .AND. &
                grid_dx( 3, ix, iy, iz ) /= &
                grid_dx( 3, min_ix_z, min_iy_z, min_iz_z ) )THEN
              CYCLE
            ENDIF
            IF( .FALSE. .AND. export_constraints_x .AND. &
                ( grid_dx( 3, ix, iy, iz ) /= &
                  grid_dx( 3, min_ix_z, min_iy_z, min_iz_z ) &
                  .OR. &
                  grid_dx( 2, ix, iy, iz ) /= &
                  grid_dx( 2, min_ix_y, min_iy_y, min_iz_y ) ) )THEN
              CYCLE
            ENDIF
            WRITE( UNIT = unit_cauchy_parts_ct, IOSTAT = ios, &
                   IOMSG = err_msg, FMT = * )&
                grid_dx( 1, ix, iy, iz ), &
                grid_dx( 2, ix, iy, iz ), &
                grid_dx( 3, ix, iy, iz ), &
                convergence_factor( ix, iy, iz )
            CALL test_status( ios, err_msg, "...error in writing " &
                        // "the arrays in " // TRIM(name_cauchy_parts_ct) )
          ENDDO
        ENDDO
      ENDDO

      CLOSE( UNIT= unit_cauchy_parts_ct )

      PRINT *, " * Convergence factor exported to formatted file ", &
                                                  TRIM(name_cauchy_parts_ct)
      PRINT *

    CASE DEFAULT

      PRINT *, "** There is no well defined algorithm " &
               // "corresponding to the number", use_constraints
      PRINT *, " * Please set use_constraints to 1 or 2."
      STOP

    END SELECT choose_constraints

    DEALLOCATE( grid_dx )

  END SUBROUTINE cauchy_convergence_test_unknown


  SUBROUTINE read_bns_id_parameters()

    IMPLICIT NONE

    INTEGER:: stat
    INTEGER, PARAMETER:: unit_parameters= 17

    CHARACTER( LEN= : ), ALLOCATABLE:: lorene_bns_id_parameters
    CHARACTER( LEN= : ), ALLOCATABLE:: msg

    lorene_bns_id_parameters= 'lorene_bns_id_parameters.par'

    INQUIRE( FILE= lorene_bns_id_parameters, EXIST= file_exists )

    IF( file_exists )THEN

     OPEN( UNIT= unit_parameters, FILE= lorene_bns_id_parameters, &
           STATUS= 'OLD' )

    ELSE

     PRINT*
     PRINT*,'** ERROR: ', lorene_bns_id_parameters, " file not found!"
     PRINT*
     STOP

    ENDIF

    READ( UNIT= unit_parameters, NML= bns_parameters, IOSTAT= stat, &
          IOMSG= msg )

    IF( stat /= 0 )THEN
      PRINT *, "** ERROR: Error in reading ",lorene_bns_id_parameters,&
               ". The IOSTAT variable is ", stat, &
               "The error message is", msg
      STOP
    ENDIF

    CLOSE( UNIT= unit_parameters )

    DO itr= 1, max_length, 1
      IF( TRIM(filenames(itr)).NE."0" )THEN
          cnt= cnt + 1
      ENDIF
    ENDDO
    IF( cnt.NE.n_bns )THEN
      PRINT *, "** ERROR! The number of file names is", cnt, &
               "and n_bns=", n_bns, ". The two should be the same."
      PRINT *
      STOP
    ENDIF

   !DO itr= 1, n_bns, 1
   !  DO itr2= 1, max_n_parts, 1
   !    IF( placer( itr, itr2 ) == test_int )THEN
   !      PRINT *
   !      PRINT *, "** ERROR! The array placer does not have ", &
   !               "enough components to specify all the desired ", &
   !               "particle distributions. Specify the ", &
   !               "components in file lorene_bns_id_particles.par"
   !      PRINT *
   !      STOP
   !    ENDIF
   !  ENDDO
   !ENDDO

  END SUBROUTINE read_bns_id_parameters

END PROGRAM convergence_test
