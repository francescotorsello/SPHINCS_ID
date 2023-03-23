! File:         sphincs_id.f90
! Author:       Francesco Torsello (FT)
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

PROGRAM sphincs_id

  !*****************************************************
  !
  !# Set up the |sph| and |bssn| |id| to be read
  !  by |sphincsbssn|.
  !
  !  FT 28.10.2020
  !
  !*****************************************************


  USE, INTRINSIC:: ISO_FORTRAN_ENV,  ONLY: COMPILER_VERSION, COMPILER_OPTIONS

#ifdef __INTEL_COMPILER

  USE IFPORT,                  ONLY: MAKEDIRQQ

#endif

#if flavour == 1

  USE sphincs_id_full,         ONLY: allocate_idbase

#elif flavour == 2

  USE sphincs_id_lorene,       ONLY: allocate_idbase

#elif flavour == 3

  USE sphincs_id_fuka,         ONLY: allocate_idbase

#elif flavour == 4

  USE sphincs_id_interpolate,  ONLY: allocate_idbase

#endif

  USE id_base,          ONLY: idbase, initialize
  USE sph_particles,    ONLY: particles
  USE bssn_formulation, ONLY: bssn
  USE timing,           ONLY: timer
  USE utility,          ONLY: date, time, zone, values, run_id, itr, itr3, &
                              itr4, hostname, version, &
                              test_status, show_progress, end_time, &
                              read_sphincs_id_parameters, &
                              !----------
                              n_id, common_path, filenames, &
                              eos_filenames, placer, &
                              export_bin, export_form, export_form_xy, &
                              export_form_x, export_constraints_xy, &
                              export_constraints_x, compute_constraints, &
                              export_constraints, export_constraints_details, &
                              constraints_step, compute_parts_constraints, &
                              numerator_ratio_dx, denominator_ratio_dx, &
                              one_lapse, zero_shift, show_progress, ref_lev, &
                              run_sph, run_spacetime, sph_path, &
                              spacetime_path, estimate_length_scale, &
                              test_int, max_n_parts


  IMPLICIT NONE


  INTEGER:: i_matter
  !! Index running over the number of physical systems

  DOUBLE PRECISION, DIMENSION(3):: adm_mom_m2p
  !# ADM linear momentum of the fluid computed using the metric mapped
  !  with the mesh-to-particle mapping

  CHARACTER(LEN=:), DIMENSION(:), ALLOCATABLE:: systems, systems_name
  !! String storing the name of the phyical systems
  CHARACTER(LEN=500):: namefile_parts
  !# String storing the name for the formatted file containing the |sph|
  !  particle |id|
  CHARACTER(LEN=500):: namefile_parts_bin
  !# String storing the name for the binary file containing the |sph|
  !  particle |id|
  CHARACTER(LEN=500):: namefile_sph
  !# String storing the name for ??
  !
  CHARACTER(LEN=500):: namefile_recovery
  !# String storing the name for the formatted file containing the data
  !  from the recovery test
  CHARACTER(LEN=500):: namefile_bssn
  !# String storing the name for the formatted file containing the |bssn| |id|
  CHARACTER(LEN=500):: namefile_bssn_bin
  !# String storing the name for the binary file containing the |bssn| |id|
  CHARACTER(LEN=500):: name_logfile
  !# String storing the name for the formatted file containing a summary about
  !  the |bssn| constraints violations

  LOGICAL:: exist
  LOGICAL(4):: dir_out

  TYPE id
    CLASS(idbase), ALLOCATABLE:: idata
  END TYPE id
  TYPE(id), DIMENSION(:), ALLOCATABLE:: ids

  TYPE(particles), DIMENSION(:,:), ALLOCATABLE:: particles_dist
  !# Array storing the particles objects,
  !  containing the particle distributions for each idbase object.
  !  Multiple particle objects can contain different particle distributions
  !  for the same idbase object.

  TYPE(bssn), DIMENSION(:), ALLOCATABLE:: bssn_forms
  !# Array storing the bssn objects,
  !  containing the BSSN variables on the gravity grid for each idbase object

  TYPE(timer):: execution_timer

  !---------------------------!
  !--  End of declarations  --!
  !---------------------------!

  CALL DATE_AND_TIME( date, time, zone, values )
  run_id= date//"-"//time

  !
  !-- Convert preprocessor strings to Fortran strings
  !

  !CALL HOSTNM( hostname )
#ifdef host

#ifdef __GFORTRAN__

# define stringize_start(x) "&
# define stringize_end(x) &x"

  hostname= stringize_start(host)
stringize_end(host)

#else

#define stringize(x) tostring(x)
#define tostring(x) #x

  hostname= stringize(host)

#endif

#else

  hostname= "unspecified host."

#endif

#ifdef vers

#ifdef __GFORTRAN__

# define stringize_start(x) "&
# define stringize_end(x) &x"

  version= stringize_start(vers)
stringize_end(vers)

#else

#define stringize(x) tostring(x)
#define tostring(x) #x

  version= stringize(vers)

#endif

#else

  hostname= "unspecified version."

#endif


  PRINT *, "  ________________________________________________________________ "
  PRINT *, "             ____________  ________  __________    __ ___          "
  PRINT *, "            / ___/ _  / /_/ / / __ \/ ___/ ___/   / / __ \         "
  PRINT *, "           (__  ) ___/ __  / / / / / /__(__  )___/ / /_/ /         "
  PRINT *, "          /____/_/  /_/ /_/_/_/ /_/____/____/___/_/_____/          "
  PRINT *
  PRINT *, "  Smoothed Particle Hydrodynamics IN Curved Spacetime              "
  PRINT *, "  Initial Data builder, ", TRIM(version)
  PRINT *
  PRINT *, "  SPHINCS_ID  Copyright (C) 2020, 2021, 2022  Francesco Torsello   "
  PRINT *
  PRINT *, "  SPHINCS_ID is free software: you can redistribute it and/or      "
  PRINT *, "  modify it under the terms of the GNU General Public License      "
  PRINT *, "  as published by the Free Software Foundation, either version     "
  PRINT *, "  of the License, or (at your option) any later version.           "
  PRINT *
  PRINT *, "  SPHINCS_ID is distributed in the hope that it will be useful,    "
  PRINT *, "  but WITHOUT ANY WARRANTY; without even the implied warranty of   "
  PRINT *, "  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU "
  PRINT *, "  General Public License for more details.                         "
  PRINT *
  PRINT *, "  You should have received a copy of the GNU General Public License"
  PRINT *, "  along with SPHINCS_ID. If not, see https://www.gnu.org/licenses/."
  PRINT *, "  The copy of the GNU General Public License should be in the file "
  PRINT *, "  'COPYING'.                                                       "
  PRINT *, "  ________________________________________________________________ "
  PRINT *
  PRINT *, "  SPHINCS_ID was compiled with: "
  PRINT *, COMPILER_VERSION()
  PRINT *
  PRINT *, "  using the options: "
  PRINT *, COMPILER_OPTIONS()
  PRINT *
  PRINT *, "  SPHINCS_ID was run on: ", TRIM(hostname)
  PRINT *, "  ________________________________________________________________ "
  PRINT *
  PRINT *, "  Run id: ", run_id
  PRINT *, "  ________________________________________________________________ "
  PRINT *

  execution_timer= timer( "execution_timer" )
  CALL execution_timer% start_timer()

  CALL read_sphincs_id_parameters()

  !
  !-- Check that the specified subdirectories exist. If not, create them.
  !-- The reason for having different versions for different compilers
  !-- is that different compilers require different versions of the
  !-- INQUIRE command. please see below.
  !

#ifdef __INTEL_COMPILER

  INQUIRE( DIRECTORY= TRIM(sph_path), EXIST= exist )
  IF( .NOT.exist )THEN
    dir_out= MAKEDIRQQ( TRIM(sph_path) )
  ELSE
    dir_out= .TRUE.
  ENDIF
  IF( .NOT.dir_out )THEN
    PRINT *, "** ERROR! Failed to create subdirectory ", TRIM(sph_path)
    PRINT *, "Stopping..."
    PRINT *
    STOP
  ENDIF

  INQUIRE( DIRECTORY= TRIM(spacetime_path), EXIST= exist )
  IF( .NOT.exist )THEN
    dir_out= MAKEDIRQQ( TRIM(spacetime_path) )
  ELSE
    dir_out= .TRUE.
  ENDIF
  IF( .NOT.dir_out )THEN
    PRINT *, "** ERROR! Failed to create subdirectory ", TRIM(sph_path)
    PRINT *, "Stopping..."
    PRINT *
    STOP
  ENDIF

#endif

#ifdef __GFORTRAN__

  INQUIRE( FILE= TRIM(sph_path)//"/.", EXIST= exist )
  IF( .NOT.exist )THEN
    CALL EXECUTE_COMMAND_LINE("mkdir "//TRIM(sph_path))
  ENDIF
  INQUIRE( FILE= TRIM(sph_path)//"/.", EXIST= exist )
  IF( .NOT.exist )THEN
    PRINT *, "** ERROR! Failed to create subdirectory ", TRIM(sph_path)
    PRINT *, " * Please create it and re-run the executable."
    PRINT *, " * Stopping..."
    PRINT *
    STOP
  ENDIF

  INQUIRE( FILE= TRIM(spacetime_path)//"/.", EXIST= exist )
  IF( .NOT.exist )THEN
    CALL EXECUTE_COMMAND_LINE("mkdir "//TRIM(spacetime_path))
  ENDIF
  INQUIRE( FILE= TRIM(spacetime_path)//"/.", EXIST= exist )
  IF( .NOT.exist )THEN
    PRINT *, "** ERROR! Failed to create subdirectory ", TRIM(spacetime_path)
    PRINT *, " * Please create it and re-run the executable."
    PRINT *, " * Stopping..."
    PRINT *
    STOP
  ENDIF

#endif


  !
  !-- Check that some parameters are consistent
  !
  IF( ( (compute_parts_constraints .EQV. .TRUE.) .AND. (run_sph .EQV. .FALSE.) &
        .AND. (run_spacetime .EQV. .TRUE.) ) )THEN

    PRINT *
    PRINT *, "** WARNING: The variable `compute_parts_constraints` is ", &
             ".TRUE., but `run_sph` is .FALSE. . "
    PRINT *, "   Please set both variables to .TRUE. to compute the ", &
             "constraints using the particle data mapped to the mesh."
    PRINT *
    STOP

  ENDIF


  !
  !-- Allocate needed arrays
  !
  ALLOCATE( CHARACTER(5):: systems(n_id) )
  ALLOCATE( CHARACTER(5):: systems_name(n_id) )

  ALLOCATE( ids( n_id ) )
  ALLOCATE( particles_dist( n_id, max_n_parts ) )
  ALLOCATE( bssn_forms    ( n_id ) )


  !
  !-- Construct the idbase objects
  !
  build_idbase_loop: DO itr= 1, n_id, 1

    CALL allocate_idbase &
      ( ids(itr)% idata, TRIM(filenames(itr)), systems(itr), systems_name(itr) )

    PRINT *, "===================================================" &
             //"==============="
    PRINT *, " Constructing idbase object for "//systems(itr), itr
    PRINT *, "===================================================" &
             //"==============="
    PRINT *

    CALL ids(itr)% idata% &
      initialize( TRIM(common_path)//TRIM(filenames(itr)), eos_filenames )

    CALL ids(itr)% idata% set_one_lapse( one_lapse )
    CALL ids(itr)% idata% set_zero_shift( zero_shift )
    CALL ids(itr)% idata% set_estimate_length_scale( estimate_length_scale )

  ENDDO build_idbase_loop


  IF( run_spacetime )THEN

    !
    !-- Construct the bssn objects from the bns objects
    !
    construct_spacetime_id_loop: DO itr3 = 1, n_id, 1
      PRINT *, "===================================================" &
               //"==============="
      PRINT *, " Setting up BSSN object for "//systems(itr3), itr3
      PRINT *, "===================================================" &
               //"==============="
      PRINT *
      bssn_forms(itr3)= bssn( ids(itr3)% idata )
    ENDDO construct_spacetime_id_loop

    !
    !-- Compute the BSSN initial data, optionally print it to a binary file
    !-- readable by SPHINCS_BSSN, and optionally read the content of such binary
    !-- file and print it to a formatted file (the latter for debugging)
    !
    compute_print_bssn_loop: DO itr3 = 1, n_id, 1
      PRINT *, "===================================================" &
               //"==============="
      PRINT *, " Computing BSSN variables for "//systems(itr3), itr3
      PRINT *, "===================================================" &
               //"==============="
      PRINT *
      WRITE( namefile_bssn_bin, "(A15)" ) "BSSN_vars.00000"
      !"BSSN_l", itr3, ".bin""(A6,I1,A4)"
      namefile_bssn_bin= TRIM( spacetime_path )//TRIM( namefile_bssn_bin )

      bssn_forms(itr3)% export_form_xy= export_form_xy
      bssn_forms(itr3)% export_form_x = export_form_x
      bssn_forms(itr3)% export_bin    = export_bin

      CALL bssn_forms(itr3)% &
        compute_and_print_tpo_variables( namefile_bssn_bin )
      !IF( bssn_forms(itr3)% export_bin )THEN
      !  WRITE( namefile_bssn, "(A10,I1,A4)" ) "bssn_vars-", itr3, ".dat"
      !  CALL bssn_forms(itr3)% read_bssn_dump_print_formatted &
      !    ( namefile_bssn_bin, namefile_bssn )
      !ENDIF
    ENDDO compute_print_bssn_loop

    !
    !-- Print the BSSN initial data to a formatted file
    !
    IF( export_form )THEN
      export_bssn_loop: DO itr3 = 1, n_id, 1
        WRITE( namefile_bssn, "(A8,I1,A4)" ) &
                              "bssn-id_", itr3, ".dat"
        namefile_bssn= TRIM( spacetime_path )//TRIM( namefile_bssn )

        CALL bssn_forms(itr3)% &
          print_formatted_id_tpo_variables( namefile= namefile_bssn )

      ENDDO export_bssn_loop
    ENDIF

  ENDIF


  IF( run_sph )THEN

    !
    !-- Construct the particles objects
    !
    place_hydro_id_loops: DO itr3= 1, n_id, 1
      part_distribution_loop: DO itr4= 1, max_n_parts, 1

        IF( placer( itr3, itr4 ) == test_int )THEN

          EXIT part_distribution_loop

        ELSE

          PRINT *, "===================================================" &
                   //"==============="
          PRINT *, " Placing particles for "//systems(itr3), itr3, &
                   ", distribution", itr4
          PRINT *, "===================================================" &
                   //"==============="
          PRINT *

          particles_dist( itr3, itr4 )= &
            particles( ids(itr3)% idata, placer( itr3, itr4 ) )

          !namefile_parts_bin= "sph-output/NSNS-2M-closeghost.00000"
          !particles_dist( itr3, itr4 )= &
          !  particles( ids(itr3)% idata, namefile_parts_bin )

        ENDIF

      ENDDO part_distribution_loop
    ENDDO place_hydro_id_loops

    !namefile_parts_bin= "sph-output/NSNS-2M-closeghost.00000"
    !namefile_parts= "try.dat"
    !CALL particles_dist(1,1)% read_sphincs_dump_print_formatted &
    !  ( namefile_parts_bin, namefile_parts )
    !STOP

  ENDIF

  IF( run_spacetime .AND. estimate_length_scale )THEN

    !
    !-- Compute the Ricci scalar on the mesh, to estimate typical length scale
    !-- of the system
    !
    compute_ricci_loop: DO itr3 = 1, n_id, 1

      PRINT *, "===================================================" &
               //"==============="
      PRINT *, " Computing Ricci tensor and scalar for "//systems(itr3), itr3
      PRINT *, "===================================================" &
               //"==============="
      PRINT *

      CALL bssn_forms(itr3)% compute_ricci()

    ENDDO compute_ricci_loop

    STOP

  ENDIF

  IF( .NOT.estimate_length_scale )THEN

    IF( run_sph )THEN

      compute_print_sph_loops: DO itr3= 1, n_id, 1
        part_distribution_loop2: DO itr4= 1, max_n_parts, 1
          IF( placer( itr3, itr4 ) == test_int )THEN
            EXIT part_distribution_loop2
            ! Experimental: empty particles object
            !particles_dist( itr, itr2 )= particles()
          ELSE

            PRINT *, "===================================================" &
                     //"====================="
            PRINT *, " Computing SPH variables for "//systems(itr3), itr3, &
                     ", distribution", itr4
            PRINT *, "===================================================" &
                     //"====================="
            PRINT *
            !WRITE( namefile_parts_bin, "(A1,I1,A1,I1,A1)" ) &
            !                            "l", &
            !                            itr3, "-", itr4, "."
            WRITE( namefile_parts_bin, "(A5)" ) systems_name(itr3)
            namefile_parts_bin= TRIM( sph_path )//TRIM( namefile_parts_bin )

            particles_dist( itr3, itr4 )% export_bin    = export_bin
            particles_dist( itr3, itr4 )% export_form_xy= export_form_xy
            particles_dist( itr3, itr4 )% export_form_x = export_form_x

            CALL particles_dist( itr3, itr4 )% &
              compute_and_print_sph_variables( namefile_parts_bin )
            !IF( particles_dist( itr3, itr4 )% export_bin )THEN
            !  WRITE( namefile_parts, "(A10,I1,A1,I1,A4)" ) &
            !                  "sph_vars-", itr3, "-", itr4, ".dat"
            !  CALL particles_dist( itr3, itr4 )% &
            !    read_sphincs_dump_print_formatted &
            !      ( namefile_parts_bin, namefile_parts )
            !ENDIF

          ENDIF
        ENDDO part_distribution_loop2
      ENDDO compute_print_sph_loops

      !
      !-- Print the particle initial data to a formatted file
      !
      IF( export_form )THEN
        print_sph_loops: DO itr3= 1, n_id, 1
          DO itr4= 1, max_n_parts, 1
            IF( placer( itr3, itr4 ) == test_int )THEN
              EXIT
              ! Experimental: empty particles object
              !particles_dist( itr, itr2 )= particles()
            ELSE
              WRITE( namefile_parts, "(A7,I1,A1,I1,A4)" ) &
                                     "sph-id_", &
                                     itr3, "-", itr4, ".dat"
              namefile_parts= TRIM( sph_path )//TRIM( namefile_parts )
              CALL particles_dist( itr3, itr4 )% &
                print_formatted_id_particles( namefile_parts )
            ENDIF
          ENDDO
        ENDDO print_sph_loops
      ENDIF

    ENDIF


    IF( run_spacetime )THEN

      !
      !-- Compute the BSSN constraints
      !
      compute_print_bssn_constraints_loop: DO itr3 = 1, n_id, 1

        bssn_forms(itr3)% cons_step         = constraints_step
        bssn_forms(itr3)% export_constraints= export_constraints
        bssn_forms(itr3)% export_constraints_details= &
                          export_constraints_details
        bssn_forms(itr3)% export_constraints_xy= export_constraints_xy
        bssn_forms(itr3)% export_constraints_x = export_constraints_x

        IF( compute_constraints )THEN

          PRINT *, "===================================================" &
                   //"==============="
          PRINT *, " Computing BSSN constraints for BSSN formulation", itr3
          PRINT *, "===================================================" &
                   //"==============="
          PRINT *

          WRITE( namefile_bssn, "(A17,I1,A4)" ) "bssn-constraints-", itr3, &
                                                ".dat"
          WRITE( name_logfile, "(A28,I1)" ) &
                              "bssn-constraints-statistics-", itr3

          namefile_bssn= TRIM(spacetime_path)//TRIM(namefile_bssn)
          name_logfile = TRIM(spacetime_path)//TRIM(name_logfile)

          CALL bssn_forms(itr3)% compute_and_print_tpo_constraints &
            ( ids(itr3)% idata, namefile_bssn, name_logfile )

        ENDIF

        part_distribution_loop3: DO itr4= 1, max_n_parts, 1

          IF( placer( itr3, itr4 ) == test_int )THEN
            EXIT
            ! Experimental: empty particles object
            !particles_dist( itr, itr2 )= particles()
          ELSE

            IF( compute_parts_constraints .AND. run_sph )THEN

              PRINT *, "===================================================" &
                       //"================================================"
              PRINT *, " Computing BSSN constraints for BSSN", &
                       " formulation", itr3, "with particle distribution", itr4
              PRINT *, "===================================================" &
                       //"================================================"
              PRINT *

              WRITE( namefile_bssn, "(A23,I1,A1,I1,A4)" ) &
                                                    "bssn-constraints-parts-", &
                                                    itr3, "-", itr4, ".dat"
              WRITE( namefile_sph, "(A12,I1,A1,I1,A4)" ) "sph-density-", itr3, &
                                                    "-", itr4, ".dat"
              WRITE( name_logfile, "(A34,I1,A1,I1,A4)" ) &
                                   "bssn-constraints-parts-statistics-", itr3, &
                                   "-", itr4

              namefile_bssn= TRIM( spacetime_path )//TRIM( namefile_bssn )
              namefile_sph = TRIM( sph_path )//TRIM( namefile_sph )
              name_logfile = TRIM( spacetime_path )//TRIM( name_logfile )

              CALL bssn_forms(itr3)% compute_and_print_tpo_constraints &
                ( particles_dist( itr3, itr4 ), namefile_bssn, name_logfile )

            ENDIF
          ENDIF

        ENDDO part_distribution_loop3

      ENDDO compute_print_bssn_constraints_loop

      !
      !-- Test recovery using the mesh-2-particle mapping
      !
      IF( run_sph )THEN

        test_recovery_m2p: DO itr3 = 1, n_id, 1

          part_distribution_loop4: DO itr4= 1, max_n_parts, 1
            IF( placer( itr3, itr4 ) == test_int )THEN
              EXIT part_distribution_loop4
              ! Experimental: empty particles object
              !particles_dist( itr, itr2 )= particles()
            ELSE

              IF(particles_dist( itr3, itr4 )% get_compose_eos()) CYCLE
              ! TODO: as of 03.03.2023, SPHINCS_BSSN does not support tabulated
              !       EOS, hence the recovery should not be called when using
              !       tabulated EOS

              PRINT *, "===================================================" &
                       //"================================================"
              PRINT *, " Testing recovery using mesh-to-particle mapping, for",&
                       " BSSN formulation", itr3, &
                       "with particle distribution", itr4
              PRINT *, "===================================================" &
                       //"================================================"
              PRINT *

              WRITE( namefile_recovery, "(A18,I1,A4)" ) &
                              "recovery-test-m2p_", itr3
              namefile_recovery= TRIM( spacetime_path ) &
                               //TRIM( namefile_recovery )
              CALL bssn_forms(itr3)% test_recovery_m2p &
                ( particles_dist( itr3, itr4 ), namefile_recovery )

              PRINT *, "===================================================" &
                       //"================================================"
              PRINT *, " Estimating the ADM momentum of the fluid using ", &
                       " the metric mapped with mesh-to-particle mapping, for",&
                       " BSSN formulation", itr3, &
                       "with particle distribution", itr4
              PRINT *, "===================================================" &
                       //"================================================"
              PRINT *

              CALL bssn_forms(itr3)% compute_adm_momentum_fluid_m2p &
                ( particles_dist( itr3, itr4 ), adm_mom_m2p )

            ENDIF

          ENDDO part_distribution_loop4

        ENDDO test_recovery_m2p

      ENDIF

    ENDIF

  ENDIF

  CALL execution_timer% stop_timer()

  CALL DATE_AND_TIME( date, time, zone, values )
  end_time= date//"-"//time

  !
  !-- Print the timers
  !

  DO itr= 1, n_id, 1

    PRINT *, "===================================================" &
             //"================================================"
    PRINT *, " Timing for physical system ", itr
    PRINT *, "===================================================" &
             //"================================================"
    PRINT *
    PRINT *, " * ID:"
    CALL ids(itr)% idata% construction_timer% print_timer( 2 )
    PRINT *
    IF( run_sph )THEN
      PRINT *, " * SPH:"
      CALL particles_dist(itr,1)% placer_timer% print_timer( 2 )
      CALL particles_dist(itr,1)% same_particle_timer% print_timer( 2 )
      DO i_matter= 1, ids(itr)% idata% get_n_matter(), 1
        CALL particles_dist(itr,1)% apm_timers(i_matter)% print_timer( 2 )
      ENDDO
      CALL particles_dist(itr,1)% importer_timer% print_timer( 2 )
      CALL particles_dist(itr,1)% sph_computer_timer% print_timer( 2 )
      PRINT *
    ENDIF
    IF( run_spacetime )THEN
      PRINT *, " * Spacetime:"
      CALL bssn_forms(itr)% grid_timer% print_timer( 2 )
      CALL bssn_forms(itr)% importer_timer% print_timer( 2 )
      CALL bssn_forms(itr)% bssn_computer_timer% print_timer( 2 )
      PRINT *
    ENDIF
    PRINT *, " * Total:"
    CALL execution_timer% print_timer( 2 )
    PRINT *

  ENDDO

  !
  !-- Print a summary
  !
  DO itr= 1, n_id, 1

    PRINT *, "===================================================" &
             //"================================================"
    PRINT *, " Summary for physical system ", itr
    PRINT *, "===================================================" &
             //"================================================"
    PRINT *
    PRINT *, "   Used ID data file: " &
             //TRIM(common_path)//TRIM(filenames(itr))
    PRINT *

    CALL ids(itr)% idata% print_summary()

    IF( run_sph )THEN

      CALL particles_dist(itr,1)% print_summary()

    ENDIF

    IF( run_spacetime )THEN

      CALL bssn_forms(itr)% print_summary()
      IF( run_sph )THEN
        PRINT *, "   Estimate of the ADM momentum of the fluid computed ", &
                 "using the SPH hydro fields and the metric mapped with ", &
                 "mesh-to-particle mapping= "
        PRINT *, "   (", adm_mom_m2p(1), ","
        PRINT *, "    ", adm_mom_m2p(2), ","
        PRINT *, "    ", adm_mom_m2p(3), ") Msun*c"
        PRINT *
      ENDIF

    ENDIF

  ENDDO
  PRINT *, "** Run started on ", run_id, " and ended on ", end_time
  PRINT *

  !
  !-- Deallocate memory
  !
  DO itr= 1, n_id, 1
    !
    !-- Destruct the LORENE Bin_NS object by hand, since the pointer to it is
    !-- global (because it is bound to C++) and cannot be nullified by the
    !-- destructor of bns. In case of multiple bns objects, this would lead
    !-- to problems...
    !-- TODO: fix this
    !
    !CALL binaries(itr)% destruct_binary()
  ENDDO
  IF( ALLOCATED( ids ) )THEN
    DEALLOCATE( ids )
  ENDIF
  IF( ALLOCATED( particles_dist ) )THEN
    DEALLOCATE( particles_dist )
  ENDIF
  IF( ALLOCATED( bssn_forms ) )THEN
    DEALLOCATE( bssn_forms )
  ENDIF


END PROGRAM sphincs_id
