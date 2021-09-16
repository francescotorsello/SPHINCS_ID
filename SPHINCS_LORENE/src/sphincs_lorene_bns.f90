! File:         setup_lorene_id.f90
! Author:       Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

PROGRAM sphincs_lorene_bns

  !*****************************************************
  !                                                    *
  ! Use the MODULE sphincs_lorene to export binary     *
  ! files containing the initial data (ID) required    *
  ! by the evolution code in SPHINCS, and built using  *
  ! the binary files produced by LORENE and containing *
  ! the binary neutron stars (BNS) ID.                 *
  !                                                    *
  ! FT 28.10.2020                                      *
  !                                                    *
  !*****************************************************

  USE IFPORT,         ONLY: MAKEDIRQQ
  USE sphincs_lorene
  USE constants,      ONLY: lorene2hydrobase, c_light2, k_lorene2hydrobase, &
                            k_lorene2hydrobase_piecewisepolytrope, MSun_geo, &
                            kg2g, m2cm, m0c2

  IMPLICIT NONE

  ! Maximum length for strings, and for the number of imported binaries
  INTEGER, PARAMETER:: max_length= 50
  ! Maximum number of binary systems
  INTEGER, PARAMETER:: max_n_bns= 50
  ! Maximum number of particle distributions
  INTEGER, PARAMETER:: max_n_parts= 250
  ! Number of binary systems of neutron stars (BNS) to import
  INTEGER:: n_bns, ref_lev
  ! Export the constraints every constraints_step-th step
  INTEGER:: constraints_step, last_level

  ! Matrix storing the information on how to place particles for each bns
  ! object. Row i contains information about the i^th bns object.
  INTEGER, PARAMETER:: test_int= - 112
  INTEGER, DIMENSION( max_n_bns, max_n_parts ):: placer= test_int

  ! Rational ratio between the large grid spacing and the medium one,
  ! equal to the ratio between the medium grid spacing nd the small one
  ! Not used in this PROGRAM, but needed since the PROGRAM reads the same
  ! parameter ile as the convergence_test PROGRAM
  DOUBLE PRECISION:: numerator_ratio_dx
  DOUBLE PRECISION:: denominator_ratio_dx

  ! Strings storing different names for output files
  CHARACTER( LEN= 500 ):: namefile_parts, namefile_parts_bin, namefile_sph
  CHARACTER( LEN= 500 ):: namefile_bssn, namefile_bssn_bin, name_logfile
  ! Array of strings storing the names of the LORENE BNS ID binary files
  CHARACTER( LEN= max_length ), DIMENSION( max_length ):: filenames= "0"
  ! String storing the local path to the directory where the
  ! LORENE BNS ID files are stored
  CHARACTER( LEN= max_length ):: common_path
  ! String storing the local path to the directory where the
  ! SPH output is to be saved
  CHARACTER( LEN= max_length ):: sph_path
  ! String storing the local path to the directory where the
  ! spacetime output is to be saved
  CHARACTER( LEN= max_length ):: spacetime_path

  LOGICAL:: exist
  LOGICAL(4):: dir_out
  ! Logical variables to steer the execution
  LOGICAL:: export_bin, export_form, export_form_xy, export_form_x, &
            compute_constraints, export_constraints_xy, &
            export_constraints_x, export_constraints, &
            export_constraints_details, compute_parts_constraints, &
            one_lapse, zero_shift, run_sph, run_spacetime

  TYPE( timer ):: execution_timer

  ! Declaration of the allocatable array storing the bns objects,
  ! containing the LORENE ID for different BNS
  TYPE( bns ),       DIMENSION(:),   ALLOCATABLE:: binaries
  ! Declaration of the allocatable array storing the particles objects,
  ! containing the particle distributions for each bns object.
  ! Multiple particle objects can contain different particle distributions
  ! for the same bns object.
  TYPE( particles ), DIMENSION(:,:), ALLOCATABLE:: particles_dist
  ! Declaration of the allocatable array storing the bssn_id objects,
  ! containing the BSSN variables on the gravity grid ofr each bns object
  TYPE( bssn_id ),   DIMENSION(:),   ALLOCATABLE:: bssn_forms

  ! Namelist containing parameters read from lorene_bns_id_parameters.par
  ! by the SUBROUTINE read_bns_id_parameters of this PROGRAM
  NAMELIST /bns_parameters/ n_bns, common_path, filenames, placer, &
                            export_bin, export_form, export_form_xy, &
                            export_form_x, export_constraints_xy, &
                            export_constraints_x, compute_constraints, &
                            export_constraints, export_constraints_details, &
                            constraints_step, compute_parts_constraints, &
                            numerator_ratio_dx, denominator_ratio_dx, ref_lev, &
                            one_lapse, zero_shift, show_progress, &
                            run_sph, run_spacetime, sph_path, spacetime_path

  !---------------------------!
  !--  End of declarations  --!
  !---------------------------!

  !PRINT *, lorene2hydrobase
  !PRINT *, 2.45191D-4/lorene2hydrobase/1000
  !PRINT *, LOG10(2.45191D-4/lorene2hydrobase/1000)
  !PRINT *
  !PRINT *, LOG10(10**(34.616)/c_light2)
  !STOP

  !PRINT *, "** Polytropic constant used for gamma= 2.75 single polytrope:"
  !PRINT *, "   k used in LORENE= ", 0.01691726009823966
  !PRINT *, "   k converted in SPHINCS units= ", &
  !                               0.01691726009823966*k_lorene2hydrobase(2.75D0)
  !PRINT *
  !PRINT *, "** Polytropic constant used for gamma= 2 single polytrope:"
  !PRINT *, "   k used in LORENE= ", 0.02686965902663748
  !PRINT *, "   k converted in SPHINCS units= ", &
  !                               0.02686965902663748*k_lorene2hydrobase(2.0D0)
  !PRINT *
  !PRINT *, "** Polytropic constant used for the crust in PWP:"
  !PRINT *, "   k used in LORENE= ", 3.99874D-8
  !PRINT *, "   k converted in SPHINCS units= ", &
  !                3.99874D-8*k_lorene2hydrobase_piecewisepolytrope(1.35692395D0)
  !PRINT *
  !PRINT *, "** Polytropic constant used for the crust in PWP:"
  !PRINT *, "   k used in LORENE= ", 8.948185D-2
  !PRINT *, "   k converted in SPHINCS units= ", &
  !                               8.948185D-2*k_lorene2hydrobase(1.35692395D0)
  !PRINT *
  !PRINT *, "   k used in LORENE, corresponding to k-100 in SPHINCS units= ", &
  !         100/k_lorene2hydrobase(2.0D0)
  ! Our testbed cases are gamma= 2.75, k= 30000; and gamma=2, k= 100
  ! in SPHINCS units
  ! 7.901e+14 density for 1.4 GRAVITATIONAL mass, poly 2
  ! 1.4-1.4 systems for both ; 1.6-1.6 ; 1.2-1.8 GRAVIATIONAL masses
  !STOP

  CALL DATE_AND_TIME( date, time, zone, values )
  run_id= date // "-" // time

  execution_timer= timer( "execution_timer" )
  CALL execution_timer% start_timer()

  CALL read_bns_id_parameters()

  !
  !-- Check that the specified subdirectories exist. If not, create them
  !-- TODO: this compils with ifort, but not with gfortran
  !
  !INQUIRE( FILE= TRIM(sph_path)//"/.", EXIST= exist )
  INQUIRE( DIRECTORY= TRIM(sph_path), EXIST= exist )
  IF( .NOT.exist ) dir_out= MAKEDIRQQ( TRIM(sph_path) )
  IF( .NOT.dir_out )THEN
    PRINT *, "** ERROR! Failed to create subdirectory ", TRIM(sph_path)
    PRINT *, "Stopping..."
    PRINT *
    STOP
  ENDIF

  !INQUIRE( FILE= TRIM(spacetime_path)//"/.", EXIST= exist )
  INQUIRE( DIRECTORY= TRIM(spacetime_path), EXIST= exist )
  IF( .NOT.exist ) dir_out= MAKEDIRQQ( TRIM(spacetime_path) )
  IF( .NOT.dir_out )THEN
    PRINT *, "** ERROR! Failed to create subdirectory ", TRIM(sph_path)
    PRINT *, "Stopping..."
    PRINT *
    STOP
  ENDIF

  ! Allocate needed memory
  ALLOCATE( binaries      ( n_bns ) )
  ALLOCATE( particles_dist( n_bns, max_n_parts ) )
  ALLOCATE( bssn_forms    ( n_bns ) )

  !
  !-- Construct the LORENE ID from the LORENE binary files
  !
  build_bns_loop: DO itr= 1, n_bns, 1
    binaries( itr )= bns( TRIM(common_path)//TRIM(filenames( itr )) )
    ! Set the variables to decide on using the geodesic gauge or not
    ! (lapse=1, shift=0)
    binaries( itr )% one_lapse = one_lapse
    binaries( itr )% zero_shift= zero_shift
  ENDDO build_bns_loop

  IF( run_sph )THEN

    !
    !-- Construct the particles objects from the bns objects
    !
    place_hydro_id_loops: DO itr3= 1, n_bns, 1
      part_distribution_loop: DO itr4= 1, max_n_parts, 1
        IF( placer( itr3, itr4 ) == test_int )THEN
          EXIT part_distribution_loop
        ELSE

          PRINT *, "===================================================" &
                   // "==============="
          PRINT *, " Placing particles for BNS", itr3, &
                   ", distribution", itr4
          PRINT *, "===================================================" &
                   // "==============="
          PRINT *
          particles_dist( itr3, itr4 )= particles( binaries( itr3 ), &
                                                   placer( itr3, itr4 ) )

        ENDIF
      ENDDO part_distribution_loop
    ENDDO place_hydro_id_loops

    compute_export_sph_loops: DO itr3= 1, n_bns, 1
      part_distribution_loop2: DO itr4= 1, max_n_parts, 1
        IF( placer( itr3, itr4 ) == test_int )THEN
          EXIT part_distribution_loop2
          ! Experimental: empty particles object
          !particles_dist( itr, itr2 )= particles()
        ELSE

          PRINT *, "===================================================" &
                   // "====================="
          PRINT *, " Computing SPH variables for BNS", itr3, &
                   ", distribution", itr4
          PRINT *, "===================================================" &
                   // "====================="
          PRINT *
          !WRITE( namefile_parts_bin, "(A1,I1,A1,I1,A1)" ) &
          !                            "l", &
          !                            itr3, "-", itr4, "."
          WRITE( namefile_parts_bin, "(A5)" ) "NSNS."
          namefile_parts_bin= TRIM( sph_path ) // TRIM( namefile_parts_bin )

          particles_dist( itr3, itr4 )% export_bin= export_bin
          particles_dist( itr3, itr4 )% export_form_xy= export_form_xy
          particles_dist( itr3, itr4 )% export_form_x = export_form_x

          CALL particles_dist( itr3, itr4 )% &
               compute_and_export_SPH_variables( namefile_parts_bin )
          !IF( particles_dist( itr3, itr4 )% export_bin )THEN
          !  WRITE( namefile_parts, "(A10,I1,A1,I1,A4)" ) &
          !                  "sph_vars-", itr3, "-", itr4, ".dat"
          !  CALL particles_dist( itr3, itr4 )% &
          !                  read_sphincs_dump_print_formatted( &
          !                                namefile_parts_bin, namefile_parts )
          !ENDIF

        ENDIF
      ENDDO part_distribution_loop2
    ENDDO compute_export_sph_loops

    !
    !-- Print the particle initial data to a formatted file
    !
    IF( export_form )THEN
      export_sph_loops: DO itr3= 1, n_bns, 1
        DO itr4= 1, max_n_parts, 1
          IF( placer( itr3, itr4 ) == test_int )THEN
            EXIT
            ! Experimental: empty particles object
            !particles_dist( itr, itr2 )= particles()
          ELSE
            WRITE( namefile_parts, "(A29,I1,A1,I1,A4)" ) &
                                   "lorene-bns-id-particles-form_", &
                                   itr3, "-", itr4, ".dat"
            namefile_parts= TRIM( sph_path ) // TRIM( namefile_parts )
            CALL particles_dist( itr3, itr4 )% &
                 print_formatted_lorene_id_particles( namefile_parts )
          ENDIF
        ENDDO
      ENDDO export_sph_loops
    ENDIF

    !PRINT *, "===================================================" &
    !         // "================================================"
    !PRINT *, " Timing "
    !PRINT *, "===================================================" &
    !         // "================================================"
    !PRINT *
    !PRINT *
    !PRINT *, " * SPH:"
    !CALL particles_dist( 1, 1 )% placer_timer% print_timer( 2 )
    !CALL particles_dist( 1, 1 )% apm1_timer% print_timer( 2 )
    !CALL particles_dist( 1, 1 )% apm2_timer% print_timer( 2 )
    !CALL particles_dist( 1, 1 )% importer_timer% print_timer( 2 )
    !CALL particles_dist( 1, 1 )% sph_computer_timer% print_timer( 2 )
    !PRINT *

  ENDIF

  IF( run_spacetime )THEN

    !
    !-- Construct the bssn_id objects from the bns objects
    !
    place_spacetime_id_loop: DO itr3 = 1, n_bns, 1
      PRINT *, "===================================================" &
               // "==============="
      PRINT *, " Setting up BSSN object for BNS", itr3
      PRINT *, "===================================================" &
               // "==============="
      PRINT *
      bssn_forms( itr3 )= bssn_id( binaries( itr3 ) )
    ENDDO place_spacetime_id_loop

    !
    !-- Compute the BSSN initial data, optionally export it to a binary file
    !-- readable by SPHINCS_BSSN, and optionally read the content of such binary
    !-- file and print it to a formatted file (the latter for debugging)
    !
    compute_export_bssn_loop: DO itr3 = 1, n_bns, 1
      PRINT *, "===================================================" &
               // "==============="
      PRINT *, " Computing BSSN variables for BNS", itr3
      PRINT *, "===================================================" &
               // "==============="
      PRINT *
      WRITE( namefile_bssn_bin, "(A15)" ) "BSSN_vars.00000"!"BSSN_l", itr3, ".bin""(A6,I1,A4)"
      namefile_bssn_bin= TRIM( spacetime_path ) // TRIM( namefile_bssn_bin )

      bssn_forms( itr3 )% export_form_xy= export_form_xy
      bssn_forms( itr3 )% export_form_x = export_form_x
      bssn_forms( itr3 )% export_bin= export_bin

      CALL bssn_forms( itr3 )% &
                          compute_and_export_3p1_variables( namefile_bssn_bin )
      !IF( bssn_forms( itr3 )% export_bin )THEN
      !  WRITE( namefile_bssn, "(A10,I1,A4)" ) "bssn_vars-", itr3, ".dat"
      !  CALL bssn_forms( itr3 )% &
      !        read_bssn_dump_print_formatted( namefile_bssn_bin, namefile_bssn )
      !ENDIF
    ENDDO compute_export_bssn_loop

    !
    !-- Print the BSSN initial data to a formatted file
    !
    IF( export_form )THEN
      export_bssn_loop: DO itr3 = 1, n_bns, 1
        WRITE( namefile_bssn, "(A24,I1,A4)" ) &
                              "lorene-bns-id-bssn-form_", itr3, ".dat"

        namefile_bssn= TRIM( spacetime_path ) // TRIM( namefile_bssn )

        CALL bssn_forms( itr3 )% &
                    print_formatted_lorene_id_3p1_variables( namefile_bssn )
      ENDDO export_bssn_loop
    ENDIF

    !
    !-- Compute the BSSN constraints
    !
    compute_export_bssn_constraints_loop: DO itr3 = 1, n_bns, 1

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

          WRITE( namefile_bssn, "(A17,I1,A4)" ) "bssn-constraints-", itr3, &
                                                ".dat"
          WRITE( name_logfile, "(A28,I1)" ) &
                              "bssn-constraints-statistics-", itr3

          namefile_bssn= TRIM( spacetime_path ) // TRIM( namefile_bssn )
          name_logfile = TRIM( spacetime_path ) // TRIM( name_logfile )

          CALL bssn_forms( itr3 )% &
                      compute_and_export_3p1_constraints( binaries( itr3 ), &
                                                          namefile_bssn, &
                                                          name_logfile )

        ENDIF

        part_distribution_loop3: DO itr4= 1, max_n_parts, 1

          IF( placer( itr3, itr4 ) == test_int )THEN
            EXIT
            ! Experimental: empty particles object
            !particles_dist( itr, itr2 )= particles()
          ELSE

          IF( compute_parts_constraints .AND. run_sph )THEN

            PRINT *, "===================================================" &
                     // "================================================"
            PRINT *, " Computing BSSN constraints for BSSN", &
                     " formulation", itr3, "with particle distribution", itr4
            PRINT *, "===================================================" &
                     // "================================================"
            PRINT *

            WRITE( namefile_bssn, "(A23,I1,A1,I1,A4)" ) &
                                                  "bssn-constraints-parts-", &
                                                  itr3, "-", itr4, ".dat"
            WRITE( namefile_sph, "(A12,I1,A1,I1,A4)" ) "sph-density-", itr3, &
                                                  "-", itr4, ".dat"
            WRITE( name_logfile, "(A34,I1,A1,I1,A4)" ) &
                                 "bssn-constraints-parts-statistics-", itr3, &
                                 "-", itr4

            namefile_bssn= TRIM( spacetime_path ) // TRIM( namefile_bssn )
            namefile_sph = TRIM( sph_path ) // TRIM( namefile_sph )
            name_logfile = TRIM( spacetime_path ) // TRIM( name_logfile )

            CALL bssn_forms( itr3 )% &
                        compute_and_export_3p1_constraints( &
                                                particles_dist( itr3, itr4 ), &
                                                namefile_bssn, &
                                                name_logfile )

          ENDIF
        ENDIF

      ENDDO part_distribution_loop3
    ENDDO compute_export_bssn_constraints_loop

  ENDIF

  CALL execution_timer% stop_timer()

  CALL DATE_AND_TIME( date, time, zone, values )
  end_time= date // "-" // time

  !
  !-- Print the timers
  !
  PRINT *, "===================================================" &
           // "================================================"
  PRINT *, " Timing "
  PRINT *, "===================================================" &
           // "================================================"
  PRINT *
  PRINT *, " * LORENE:"
  CALL binaries( 1 )% binary_construction_timer% print_timer( 2 )
  PRINT *
  IF( run_sph )THEN
    PRINT *, " * SPH:"
    CALL particles_dist( 1, 1 )% placer_timer% print_timer( 2 )
    CALL particles_dist( 1, 1 )% same_particle_timer% print_timer( 2 )
    CALL particles_dist( 1, 1 )% apm1_timer% print_timer( 2 )
    CALL particles_dist( 1, 1 )% apm2_timer% print_timer( 2 )
    CALL particles_dist( 1, 1 )% importer_timer% print_timer( 2 )
    CALL particles_dist( 1, 1 )% sph_computer_timer% print_timer( 2 )
    PRINT *
  ENDIF
  IF( run_spacetime )THEN
    PRINT *, " * Spacetime:"
    CALL bssn_forms( 1 )% grid_timer% print_timer( 2 )
    CALL bssn_forms( 1 )% importer_timer% print_timer( 2 )
    CALL bssn_forms( 1 )% bssn_computer_timer% print_timer( 2 )
    PRINT *
  ENDIF
  PRINT *, " * Total:"
  CALL execution_timer% print_timer( 2 )
  PRINT *

  !
  !-- Print a summary
  !
  PRINT *, "===================================================" &
           // "================================================"
  PRINT *, " Summary "
  PRINT *, "===================================================" &
           // "================================================"
  PRINT *
  PRINT *, " * Binary system of neutron stars:"
  PRINT *
  PRINT *, "   Used binary file produced by LORENE: " &
           // TRIM(common_path)//TRIM(filenames( 1 ))
  PRINT *
  PRINT *, "   Baryon mass of neutron star 1= ", binaries( 1 )% get_mass1(), &
           "Msun"
  PRINT *, "   Baryon mass of neutron star 2= ", binaries( 1 )% get_mass1(), &
           "Msun"
  PRINT *, "   Gravitational mass of neutron star 1= ", &
           binaries( 1 )% get_grav_mass1(), "Msun"
  PRINT *, "   Gravitational mass of neutron star 2= ", &
           binaries( 1 )% get_grav_mass1(), "Msun"
  PRINT *
  PRINT *, "   Equatorial (not areal) radius of neutron star 1 towards " &
           // "companion= ", &
               binaries( 1 )% get_radius1_x_comp(), "Msun_geo"
  PRINT *, "   Equatorial (not areal) radius of neutron star 2 towards " &
           // "companion= ", &
               binaries( 1 )% get_radius2_x_comp(), "Msun_geo"
  PRINT *, "   Equatorial (not areal) radius of neutron star 1 opposite to " &
           // "companion= ", &
               binaries( 1 )% get_radius1_x_opp(), "Msun_geo"
  PRINT *, "   Equatorial (not areal) radius of neutron star 2 opposite to " &
           // "companion= ", &
               binaries( 1 )% get_radius2_x_opp(), "Msun_geo"
  PRINT *, "   Radius (not areal) along y of neutron star 1= ", &
               binaries( 1 )% get_radius1_y(), "Msun_geo"
  PRINT *, "   Radius (not areal) along y of neutron star 2= ", &
               binaries( 1 )% get_radius2_y(), "Msun_geo"
  PRINT *, "   Radius (not areal) along y of neutron star 1= ", &
               binaries( 1 )% get_radius1_z(), "Msun_geo"
  PRINT *, "   Radius (not areal) along y of neutron star 2= ", &
               binaries( 1 )% get_radius2_z(), "Msun_geo"
  PRINT *
  PRINT *, "   EOS for neutron star 1= ", &
               binaries( 1 )% get_eos1()
  PRINT *, "   EOS for neutron star 2= ", &
               binaries( 1 )% get_eos2()
  PRINT *
  PRINT *, "   Central baryon mass density for star 1= ", &
               binaries( 1 )% get_rho_center1(), "Msun/Msun_geo**3= ", &
               binaries( 1 )% get_rho_center1() &
               /lorene2hydrobase*kg2g/(m2cm**3), "g cm^{-3}"
  PRINT *, "   Central baryon mass density for star 2= ", &
               binaries( 1 )% get_rho_center2(), "Msun/Msun_geo**3= ", &
               binaries( 1 )% get_rho_center2() &
               /lorene2hydrobase*kg2g/(m2cm**3), "g cm^{-3}"
  PRINT *
  IF( run_sph )THEN
    PRINT *, " * SPH:"
    PRINT *
    PRINT *, "   Total particle number= ", particles_dist( 1, 1 )% get_npart()
    PRINT *, "   Particle number on star 1: npart1=", &
                                          particles_dist( 1, 1 )% get_npart1()
    PRINT *, "   Particle number on star 2: npart2=", &
                                          particles_dist( 1, 1 )% get_npart2()
    PRINT *, "   Particle number ratio: npart1/npart2= ", &
                            DBLE(particles_dist( 1, 1 )% get_npart1()) &
                           /DBLE(particles_dist( 1, 1 )% get_npart2())
    PRINT *, "   Star mass ratio: mass1/mass2= ", &
                           binaries( 1 )% get_mass1()/binaries( 1 )% get_mass2()
    PRINT *
    PRINT *, "   Baryon number ratio over both stars=", &
             particles_dist( 1, 1 )% get_nuratio()
    PRINT *, "   Baryon number ratio on star 1=", &
             particles_dist( 1, 1 )% get_nuratio1()
    PRINT *, "   Baryon number ratio on star 2=", &
             particles_dist( 1, 1 )% get_nuratio2()
    PRINT *
  ENDIF
  IF( run_spacetime )THEN
    PRINT *, " * Spacetime:"
    PRINT *
    PRINT *, "   Number of refinement levels= ", bssn_forms( 1 )% get_nlevels()
    PRINT *
    PRINT *, "   Number of grid points on each level= ", &
             bssn_forms( 1 )% get_ngrid_x( 1 ), "**3"
    PRINT *
    DO itr= 1, bssn_forms( 1 )% get_nlevels(), 1
      PRINT *, "   Resolution on level ", itr, "= ", &
               bssn_forms( 1 )% get_dx(itr)
    ENDDO
    PRINT *
    DO itr= 1, bssn_forms( 1 )% get_nlevels(), 1
      PRINT *, "   x boundary of level ", itr, "= ", &
               bssn_forms( 1 )% get_xR(itr)
      PRINT *, "   y boundary of level ", itr, "= ", &
               bssn_forms( 1 )% get_yR(itr)
      PRINT *, "   z boundary of level ", itr, "= ", &
               bssn_forms( 1 )% get_zR(itr)
               !bssn_forms( 1 )% get_dx(itr)* &
               !(bssn_forms( 1 )% get_ngrid_x(itr)-2.0D0)/2.0D0, "Msun_geo= ", &
               !bssn_forms( 1 )% get_dx(itr)* &
               !(bssn_forms( 1 )% get_ngrid_x(itr)-2.0D0)/2.0D0*Msun_geo, "km "
    ENDDO
    PRINT *
    last_level= bssn_forms( 1 )% get_nlevels()
    PRINT *, "   Number of grid points across the x-axis-diameter of star 1=", &
             FLOOR( ( binaries( 1 )% get_radius1_x_comp() + &
             binaries( 1 )% get_radius1_x_opp() ) &
             /bssn_forms( 1 )% get_dx( last_level ) )
    PRINT *, "   Number of grid points across the x-axis-diameter of star 2=", &
             FLOOR( ( binaries( 1 )% get_radius2_x_comp() + &
            binaries( 1 )% get_radius2_x_opp() ) &
            /bssn_forms( 1 )% get_dx( last_level ) )
    PRINT *
  ENDIF
  PRINT *, "** Run started on ", run_id, " and ended on ", end_time
  PRINT *

  !
  !-- Deallocate memory
  !
  DO itr= 1, n_bns, 1
    !
    !-- Destruct the LORENE Bin_NS object by hand, since the pointer to it is
    !-- global (because it is bound to C++) and cannot be nullified by the
    !-- destructor of bns. In case of multiple bns objects, this would lead
    !-- to problems...
    !-- TODO: fix this
    !
    CALL binaries( itr )% destruct_binary()
  ENDDO
  IF( ALLOCATED( binaries ) )THEN
    DEALLOCATE( binaries )
  ENDIF
  IF( ALLOCATED( particles_dist ) )THEN
    DEALLOCATE( particles_dist )
  ENDIF
  IF( ALLOCATED( bssn_forms ) )THEN
    DEALLOCATE( bssn_forms )
  ENDIF


  CONTAINS


  SUBROUTINE read_bns_id_parameters()

    IMPLICIT NONE

    INTEGER:: stat

    CHARACTER( LEN= : ), ALLOCATABLE:: lorene_bns_id_parameters
    CHARACTER( LEN= 100 ):: msg

    lorene_bns_id_parameters= 'sphincs_lorene_bns_parameters.par'

    INQUIRE( FILE= lorene_bns_id_parameters, EXIST= file_exists )
    IF( file_exists )THEN
     OPEN( 17, FILE= lorene_bns_id_parameters, STATUS= 'OLD' )
    ELSE
     PRINT*
     PRINT*,'** ERROR: ', lorene_bns_id_parameters, " file not found!"
     PRINT*
     STOP
    ENDIF

    READ( 17, NML= bns_parameters, IOSTAT= stat, IOMSG= msg )
      IF( stat /= 0 )THEN
        PRINT *, "** ERROR: Error in reading ",lorene_bns_id_parameters,&
                 ". The IOSTAT variable is ", stat, &
                 "The error message is", msg
        STOP
      ENDIF
    CLOSE( 17 )

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

  END SUBROUTINE

END PROGRAM sphincs_lorene_bns
