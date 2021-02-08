! File:         setup_lorene_id.f90
! Author:       Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

PROGRAM setup_lorene_id

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

  USE sphincs_lorene

  IMPLICIT NONE

  ! Maximum length for strings, and for the number of imported binaries
  INTEGER, PARAMETER:: max_length= 50
  ! Maximum number of binary systems
  INTEGER, PARAMETER:: max_n_bns= 50
  ! Maximum number of particle distributions
  INTEGER, PARAMETER:: max_n_parts= 250
  ! Number of binary systems of neutron stars (BNS) to import
  INTEGER:: n_bns
  ! Export the constraints every constraints_step-th step
  INTEGER:: constraints_step

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
  CHARACTER( LEN= 500 ):: namefile_parts
  CHARACTER( LEN= 500 ):: namefile_bssn, name_logfile
  ! Array of strings storing the names of the LORENE BNS ID binary files
  CHARACTER( LEN= max_length ), DIMENSION( max_length ):: filenames= "0"
  ! String storing the local path to the directory where the
  ! LORENE BNS ID files are stored
  CHARACTER( LEN= max_length ):: common_path

  ! Logical variable to check if files exist
  LOGICAL:: exist
  ! Logical variables to steer the execution
  LOGICAL:: export_bin, export_form, export_form_xy, export_form_x, &
            compute_constraints, export_constraints_xy, &
            export_constraints_x, export_constraints, &
            export_constraints_details, compute_parts_constraints

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
                            numerator_ratio_dx, denominator_ratio_dx, &
                            show_progress

  !---------------------------!
  !--  End of declarations  --!
  !---------------------------!

  CALL DATE_AND_TIME( date, time, zone, values )
  run_id= date // "-" // time

  execution_timer= timer( "execution_timer" )
  CALL execution_timer% start_timer()

  CALL read_bns_id_parameters()

  ! Allocate needed memory
  ALLOCATE( binaries      ( n_bns ) )
  ALLOCATE( particles_dist( n_bns, max_n_parts ) )
  ALLOCATE( bssn_forms    ( n_bns ) )

  !
  !-- Construct the LORENE ID from the LORENE binary files
  !
  build_bns_loop: DO itr= 1, n_bns, 1
    binaries( itr )= bns( TRIM(common_path)//"/"//TRIM(filenames( itr )) )
  ENDDO build_bns_loop

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
        WRITE( namefile_parts, "(A1,I1,A1,I1,A1)" ) &
                                    "l", &
                                    itr3, "-", itr4, "."
        particles_dist( itr3, itr4 )% export_bin= export_bin
        CALL particles_dist( itr3, itr4 )% &
             compute_and_export_SPH_variables( namefile_parts )

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
          particles_dist( itr3, itr4 )% export_form_xy= export_form_xy
          particles_dist( itr3, itr4 )% export_form_x = export_form_x
          CALL particles_dist( itr3, itr4 )% &
               print_formatted_lorene_id_particles( namefile_parts )
        ENDIF
      ENDDO
    ENDDO export_sph_loops
  ENDIF

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
    bssn_forms( itr3 )= BSSN_id( binaries( itr3 ) )
  ENDDO place_spacetime_id_loop

  compute_export_bssn_loop: DO itr3 = 1, n_bns, 1
    PRINT *, "===================================================" &
             // "==============="
    PRINT *, " Computing BSSN variables for BNS", itr3
    PRINT *, "===================================================" &
             // "==============="
    PRINT *
    WRITE( namefile_bssn, "(A6,I1,A4)" ) "BSSN_l", itr3, ".bin"
    bssn_forms( itr3 )% export_bin= export_bin
    CALL bssn_forms( itr3 )% &
                        compute_and_export_3p1_variables( namefile_bssn )
    IF( bssn_forms( itr3 )% export_bin )THEN
      WRITE( namefile_bssn, "(A10,I1,A4)" ) "bssn_vars-", itr3, ".dat"
      bssn_forms( itr3 )% export_form_xy= export_form_xy
      bssn_forms( itr3 )% export_form_x = export_form_x
      CALL bssn_forms( itr3 )% read_bssn_dump_print_formatted( namefile_bssn )
    ENDIF
  ENDDO compute_export_bssn_loop

  !
  !-- Print the BSSN initial data to a formatted file
  !
  IF( export_form )THEN
    export_bssn_loop: DO itr3 = 1, n_bns, 1
      WRITE( namefile_bssn, "(A24,I1,A4)" ) &
                            "lorene-bns-id-bssn-form_", itr3, ".dat"
      CALL bssn_forms( itr3 )% &
                  print_formatted_lorene_id_3p1_variables( namefile_bssn )
    ENDDO export_bssn_loop
  ENDIF

  !
  !-- Compute the BSSN constraints
  !
  IF( compute_constraints )THEN
    compute_export_bssn_constraints_loop: DO itr3 = 1, n_bns, 1
      PRINT *, "===================================================" &
               // "==============="
      PRINT *, " Computing BSSN constraints for BNS", itr3
      PRINT *, "===================================================" &
               // "==============="
      PRINT *
      WRITE( namefile_bssn, "(A17,I1,A4)" ) "bssn-constraints-", &
                                           itr3, ".dat"
      WRITE( name_logfile, "(A28,I1,A4)" ) &
                          "bssn-constraints-statistics-", itr3, ".log"
      bssn_forms( itr3 )% cons_step= constraints_step
      bssn_forms( itr3 )% export_constraints= export_constraints
      bssn_forms( itr3 )% export_constraints_details= &
                          export_constraints_details
      bssn_forms( itr3 )% export_constraints_xy= export_constraints_xy
      bssn_forms( itr3 )% export_constraints_x = export_constraints_x
      CALL bssn_forms( itr3 )% &
                  compute_and_export_3p1_constraints( binaries( itr3 ), &
                                                      namefile_bssn, &
                                                      name_logfile )
    ENDDO compute_export_bssn_constraints_loop
  ENDIF

  CALL execution_timer% stop_timer()

  !
  !-- Print the timers
  !
  PRINT *, "** Timing."
  PRINT *
  PRINT *, " * SPH:"
  CALL particles_dist( 1, 1 )% placer_timer% print_timer( 2 )
  CALL particles_dist( 1, 1 )% importer_timer% print_timer( 2 )
  CALL particles_dist( 1, 1 )% sph_computer_timer% print_timer( 2 )
  PRINT *
  PRINT *, " * Gravity:"
  CALL bssn_forms( 1 )% grid_timer% print_timer( 2 )
  CALL bssn_forms( 1 )% importer_timer% print_timer( 2 )
  CALL bssn_forms( 1 )% bssn_computer_timer% print_timer( 2 )
  PRINT *
  PRINT *, " * Total:"
  CALL execution_timer% print_timer( 2 )
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
    CHARACTER( LEN= : ), ALLOCATABLE:: msg

    lorene_bns_id_parameters= 'lorene_bns_id_parameters.par'

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

END PROGRAM setup_lorene_id
