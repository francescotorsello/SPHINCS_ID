! File:         convergence_test.f90
! Author:       Francesco Torsello (FT)
!************************************************************************
! Copyright (C) 2020, 2021, 2022 Francesco Torsello                     *
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

PROGRAM convergence_test

  !*****************************************************
  !
  !# Make a convergence test to check the validity of
  !  the code |sphincsid|.
  !
  !  FT 8.12.2020
  !
  !*****************************************************

#ifdef __INTEL_COMPILER
  USE IFPORT,          ONLY: MAKEDIRQQ
#endif

#if flavour == 1
  USE sphincs_id_lorene,        ONLY: allocate_idbase
#elif flavour == 2
  USE sphincs_id_interpolate,   ONLY: allocate_idbase
#endif

  USE id_base,                  ONLY: idbase, initialize
  USE sph_particles,            ONLY: particles
  USE bssn_formulation,         ONLY: bssn
  USE standard_tpo_formulation, ONLY: tpo
  USE timing,                   ONLY: timer
  USE utility,                  ONLY: date, time, zone, values, run_id, itr, &
                                      itr3, file_exists, cnt, ios, err_msg, &
                                      test_status, show_progress, end_time

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
  ! Refinement level over which to compute contraints (make input parameter)
  INTEGER:: ref_lev
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

  CHARACTER( LEN= : ), DIMENSION(:), ALLOCATABLE:: systems, systems_name
  !! String storing the name of the phyical systems
  CHARACTER( LEN= 500 ):: namefile_parts
  !# String storing the name for the formatted file containing the |sph|
  !  particle |id|
  CHARACTER( LEN= 500 ):: namefile_parts_bin
  !# String storing the name for the binary file containing the |sph|
  !  particle |id|
  CHARACTER( LEN= 500 ):: namefile_sph
  !# String storing the name for ??
  !
  CHARACTER( LEN= 500 ):: namefile_bssn
  !# String storing the name for the formatted file containing the |bssn| |id|
  CHARACTER( LEN= 500 ):: namefile_bssn_bin
  !# String storing the name for the binary file containing the |bssn| |id|
  CHARACTER( LEN= 500 ):: name_logfile
  !# String storing the name for the formatted file containing a summary about
  !  the |bssn| constraints violations
  CHARACTER( LEN= max_length ), DIMENSION( max_length ):: filenames= "0"
  ! Array of strings storing the names of the |id| files
  CHARACTER( LEN= max_length ):: common_path
  !# String storing the local path to the directory where the |id| files
  !  are stored
  CHARACTER( LEN= max_length ):: sph_path
  !# String storing the local path to the directory where the
  !  SPH output is to be saved
  CHARACTER( LEN= max_length ):: spacetime_path
  !# String storing the local path to the directory where the
  !  spacetime output is to be saved

  LOGICAL, PARAMETER:: debug= .FALSE.
  LOGICAL:: exist
  LOGICAL(4):: dir_out

  ! Logical variables to steer the execution
  LOGICAL:: export_bin, export_form, export_form_xy, export_form_x, &
            compute_constraints, export_constraints_xy, &
            export_constraints_x, export_constraints, &
            export_constraints_details, compute_parts_constraints, &
            one_lapse, zero_shift, run_sph, run_spacetime

  TYPE( timer ):: execution_timer

  ! Declaration of the bns (binary neutron star) object containing the
  ! formatted ID
  !TYPE( bnslorene ):: binary
  CLASS( idbase ), ALLOCATABLE:: idata
  ! Declaration of the particles object containing the particle distribution
  TYPE( particles ):: particles_dist
  ! Declaration of the 3-component array storing the 3 bssn objects,
  ! containing the BSSN variables on the gravity grid
  TYPE( bssn ), DIMENSION(3):: bssn_forms

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

  ! Conversions of some polytropic constants from formatted units to
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

  PRINT *, "  ________________________________________________________________ "
  PRINT *, "             ____________  ________  __________    __ ___          "
  PRINT *, "            / ___/ _  / /_/ / / __ \/ ___/ ___/   / / __ \         "
  PRINT *, "           (__  ) ___/ __  / / / / / /__(__  )___/ / /_/ /         "
  PRINT *, "          /____/_/  /_/ /_/_/_/ /_/____/____/___/_/_____/          "
  PRINT *
  PRINT *, "  Smoothed Particle Hydrodynamics IN Curved Spacetime              "
  PRINT *, "  Initial Data builder, v1.0 - Cauchy convergence test             "
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
  PRINT *, "  You should have received a copy of the GNU General Public        "
  PRINT *, "  along with SPHINCS_ID. If not, see https://www.gnu.org/licenses/."
  PRINT *, "  The copy of the GNU General Public License should be in the file "
  PRINT *, "  'COPYING'.                                                       "
  PRINT *, "  ________________________________________________________________ "
  PRINT *
  PRINT *, "  Run id: ", run_id
  PRINT *, "  ________________________________________________________________ "
  PRINT *

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
    PRINT *, "** ERROR! Directory ", TRIM(sph_path), " does not exist!"
    PRINT *, "   Please create it and re-run the executable. Stopping..."
    STOP
  ENDIF

  INQUIRE( FILE= TRIM(spacetime_path)//"/.", EXIST= exist )
  IF( .NOT.exist )THEN
    PRINT *, "** ERROR! Directory ", TRIM(spacetime_path), " does not exist!"
    PRINT *, "   Please create it and re-run the executable. Stopping..."
    STOP
  ENDIF

#endif

  ALLOCATE( CHARACTER(5):: systems(1) )
  ALLOCATE( CHARACTER(5):: systems_name(1) )

  !
  !-- Construct the idbase objects
  !
  CALL allocate_idbase( idata, TRIM(filenames(1)), systems(1), systems_name(1) )
  CALL idata% initialize( TRIM(common_path)//TRIM(filenames(1)) )
  ! Set the variables to decide on using the geodesic gauge or not
  ! (lapse=1, shift=0)
  CALL idata% set_one_lapse ( one_lapse )
  CALL idata% set_zero_shift( zero_shift )

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
    particles_dist= particles( idata, placer( 1, 1 ) )

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
    WRITE( namefile_parts_bin, "(A5)" ) systems_name(1)
    namefile_parts_bin= TRIM( sph_path ) // TRIM( namefile_parts_bin )

    particles_dist% export_bin    = export_bin
    particles_dist% export_form_xy= export_form_xy
    particles_dist% export_form_x = export_form_x
    CALL particles_dist% compute_and_export_SPH_variables( namefile_parts )

    !
    !-- Print the particle initial data to a formatted file
    !
    IF( export_form )THEN
      WRITE( namefile_parts, "(A34)" ) &
                             "lorene-bns-id-particles-form_1.dat"
      namefile_parts= TRIM( sph_path ) // TRIM( namefile_parts )
      CALL particles_dist% print_formatted_id_particles( namefile_parts )
    ENDIF

  ENDIF

  !
  !-- Construct the bssn objects from the bns object
  !
  construct_bssn_loop: DO itr3 = min_bssn, max_bssn, 1

    PRINT *, "===================================================" &
             // "==============="
    PRINT *, " Setting up BSSN object ", itr3
    PRINT *, "===================================================" &
             // "==============="
    PRINT *

    IF( itr3 == 1 )THEN

      bssn_forms( itr3 )= bssn( idata )
      original_dx= bssn_forms( itr3 )% get_dx(ref_lev)

    ELSE
      IF( itr3 == min_bssn )THEN

        bssn_forms( 1 )= bssn( idata )
        original_dx= bssn_forms( 1 )% get_dx(ref_lev)

      ENDIF
      bssn_forms( itr3 )= bssn( idata, &
                                   original_dx/( ratio_dx**( itr3 - 1 ) ), &
                                   original_dx/( ratio_dx**( itr3 - 1 ) ), &
                                   original_dx/( ratio_dx**( itr3 - 1 ) ) )

      IF( bssn_forms( itr3 )% get_dx(ref_lev) /= &
          original_dx/( ratio_dx**( itr3 - 1 ) ) )THEN
        PRINT *, " ** ERROR! The grid spacing #", itr3, ",", &
                 bssn_forms( itr3 )% get_dx(ref_lev), &
                 " is not equal to dx/", ratio_dx**( itr3 - 1 ), "= ", &
                 original_dx/( ratio_dx**( itr3 - 1 ) )
        STOP
      ENDIF
    ENDIF

    PRINT *, "** The grid spacing is dx=", bssn_forms( itr3 )% get_dx(ref_lev)
    PRINT *, "** The number of grid points for dx is:", &
                            bssn_forms( itr3 )% get_ngrid_x(ref_lev), "**3"
    PRINT *

  ENDDO construct_bssn_loop

  IF( debug )THEN
    PRINT *, "bssn_forms( 1 )% get_ngrid_x=", bssn_forms( 1 )% get_ngrid_x(ref_lev)
    PRINT *, "bssn_forms( 1 )% get_ngrid_y=", bssn_forms( 1 )% get_ngrid_y(ref_lev)
    PRINT *, "bssn_forms( 1 )% get_ngrid_z=", bssn_forms( 1 )% get_ngrid_z(ref_lev)
    PRINT *, "bssn_forms( 2 )% get_ngrid_x=", bssn_forms( 2 )% get_ngrid_x(ref_lev)
    PRINT *, "bssn_forms( 2 )% get_ngrid_y=", bssn_forms( 2 )% get_ngrid_y(ref_lev)
    PRINT *, "bssn_forms( 2 )% get_ngrid_z=", bssn_forms( 2 )% get_ngrid_z(ref_lev)
    PRINT *, "bssn_forms( 3 )% get_ngrid_x=", bssn_forms( 3 )% get_ngrid_x(ref_lev)
    PRINT *, "bssn_forms( 3 )% get_ngrid_y=", bssn_forms( 3 )% get_ngrid_y(ref_lev)
    PRINT *, "bssn_forms( 3 )% get_ngrid_z=", bssn_forms( 3 )% get_ngrid_z(ref_lev)
    PRINT *
    PRINT *, "bssn_forms( 1 )% get_dx ", bssn_forms( 1 )% get_dx(ref_lev)
    PRINT *, "bssn_forms( 2 )% get_dy ", bssn_forms( 2 )% get_dx(ref_lev)
    PRINT *, "bssn_forms( 3 )% get_dz ", bssn_forms( 3 )% get_dx(ref_lev)
    PRINT *
    !STOP
  ENDIF

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
    WRITE( namefile_bssn_bin, "(A6,I1,A4)" ) "BSSN_l", itr3, ".bin"
    namefile_bssn_bin= TRIM( spacetime_path ) // TRIM( namefile_bssn_bin )

    bssn_forms( itr3 )% export_bin= export_bin
    bssn_forms( itr3 )% export_form_xy= export_form_xy
    bssn_forms( itr3 )% export_form_x = export_form_x
    CALL bssn_forms( itr3 )% &
                        compute_and_export_tpo_variables( namefile_bssn_bin )
  ENDDO compute_export_bssn_loop

  !
  !-- Print the BSSN initial data to a formatted file
  !
  IF( export_form )THEN
    export_bssn_loop: DO itr3 = min_bssn, max_bssn, 1
      WRITE( namefile_bssn, "(A24,I1,A4)" ) &
                            "lorene-bns-id-bssn-form_", itr3, ".dat"
      namefile_bssn= TRIM( spacetime_path ) // TRIM( namefile_bssn )

      CALL bssn_forms( itr3 )% &
                  print_formatted_id_tpo_variables( namefile_bssn )
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
                          "bssn-constraints-statistics-", itr3

      namefile_bssn= TRIM( spacetime_path ) // TRIM( namefile_bssn )
      name_logfile = TRIM( spacetime_path ) // TRIM( name_logfile )

      CALL bssn_forms( itr3 )% &
                  compute_and_export_tpo_constraints( idata, &
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

      namefile_bssn= TRIM( spacetime_path ) // TRIM( namefile_bssn )
      namefile_sph = TRIM( sph_path ) // TRIM( namefile_sph )
      name_logfile = TRIM( spacetime_path ) // TRIM( name_logfile )

      CALL bssn_forms( itr3 )% &
                  compute_and_export_tpo_constraints( particles_dist, &
                                                      namefile_bssn, &
                                                      name_logfile )

    ENDIF

  ENDDO compute_export_bssn_constraints_loop

! Here the bug appears. Since I am only computing the constraints with particles
! it must be in there
!STOP

  !
  !-- The BSSN formulations on grids with different resolutions are ready
  !-- to be used in a convergence test
  !
  IF( debug .AND. .FALSE. )THEN
    PRINT *, bssn_forms( 1 )% get_grid_point( 1 + 1, 1 + 2, 1 + 1, ref_lev )
    PRINT *, bssn_forms( 2 )% get_grid_point( 1 + 2, 1 + 4, 1 + 2, ref_lev )
    PRINT *, bssn_forms( 3 )% get_grid_point( 1 + 4, 1 + 8, 1 + 4, ref_lev )
    PRINT *

    PRINT *, bssn_forms( 1 )% get_grid_point( 1 + 2, 1 + 1, 1 + 2, ref_lev )
    PRINT *, bssn_forms( 2 )% get_grid_point( 1 + 4, 1 + 2, 1 + 4, ref_lev )
    PRINT *, bssn_forms( 3 )% get_grid_point( 1 + 8, 1 + 4, 1 + 8, ref_lev )
    PRINT *

    PRINT *, ABS(bssn_forms( 1 )% get_HC( 1 + 30,   1 + 25,   1 + 17, ref_lev ))
    PRINT *, ABS(bssn_forms( 2 )% get_HC( 1 + 2*30, 1 + 2*25, 1 + 2*17, ref_lev ))
    PRINT *, ABS(bssn_forms( 3 )% get_HC( 1 + 4*30, 1 + 4*25, 1 + 4*17, ref_lev ))
    PRINT *
    PRINT *, ABS(bssn_forms( 1 )% get_HC( 1 + 30,   1 + 25,   1 + 17, ref_lev )) &
           - ABS(bssn_forms( 2 )% get_HC( 1 + 2*30, 1 + 2*25, 1 + 2*17, ref_lev ))
    PRINT *, ABS(bssn_forms( 2 )% get_HC( 1 + 2*30, 1 + 2*25, 1 + 2*17, ref_lev )) &
           - ABS(bssn_forms( 3 )% get_HC( 1 + 4*30, 1 + 4*25, 1 + 4*17, ref_lev ))
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
                                          bssn_forms(3), 1, ref_lev )
    CALL cauchy_convergence_test_known( bssn_forms(2), bssn_forms(3), 1, ref_lev )

  ENDIF

  IF( compute_parts_constraints )THEN

    PRINT *, "** Performing convergence test with constraints computed ", &
             "with particle data."
    PRINT *
    CALL cauchy_convergence_test_unknown( bssn_forms(1), bssn_forms(2), &
                                          bssn_forms(3), 2, ref_lev )
    CALL cauchy_convergence_test_known( bssn_forms(2), bssn_forms(3), 2, ref_lev )

  ENDIF

  CALL execution_timer% stop_timer()

  CALL DATE_AND_TIME( date, time, zone, values )
  end_time= date // "-" // time

!STOP

  !
  !-- Print the timers
  !
  PRINT *, "===================================================" &
           // "================================================"
  PRINT *, " Timing and summaries"
  PRINT *, "===================================================" &
           // "================================================"
  PRINT *
  PRINT *
  CALL idata% print_summary()
  !PRINT *, " * BSSN formulation with uniform resolution:", &
  !         bssn_forms( 1 )% get_dx(ref_lev)
  !PRINT *, "    and number of points:", bssn_forms( 1 )% get_ngrid_x(ref_lev), &
  !         "**3"
  !original_dx
  CALL bssn_forms( 1 )% print_summary()
  CALL bssn_forms( 1 )% grid_timer% print_timer( 2 )
  CALL bssn_forms( 1 )% importer_timer% print_timer( 2 )
  CALL bssn_forms( 1 )% bssn_computer_timer% print_timer( 2 )
  PRINT *
  !PRINT *, " * BSSN formulation with uniform resolution:", &
  !         bssn_forms( 2 )% get_dx(ref_lev)
  !PRINT *, "    and number of points:", bssn_forms( 2 )% get_ngrid_x(ref_lev), &
  !         "**3"
  !original_dx/2
  CALL bssn_forms( 2 )% print_summary()
  CALL bssn_forms( 2 )% grid_timer% print_timer( 2 )
  CALL bssn_forms( 2 )% importer_timer% print_timer( 2 )
  CALL bssn_forms( 2 )% bssn_computer_timer% print_timer( 2 )
  PRINT *
  !PRINT *, " * BSSN formulation with uniform resolution:", &
  !         bssn_forms( 3 )% get_dx(ref_lev)
  !PRINT *, "    and number of points:", bssn_forms( 3 )% get_ngrid_x(ref_lev), &
  !         "**3"
  !original_dx/4
  CALL bssn_forms( 3 )% print_summary()
  CALL bssn_forms( 3 )% grid_timer% print_timer( 2 )
  CALL bssn_forms( 3 )% importer_timer% print_timer( 2 )
  CALL bssn_forms( 3 )% bssn_computer_timer% print_timer( 2 )
  PRINT *
  PRINT *, " * Total:"
  CALL execution_timer% print_timer( 2 )
  PRINT *

  PRINT *
  PRINT *, "** Run started on ", run_id, " and ended on ", end_time
  PRINT *

  !
  !-- Destruct the formatted Bin_NS object by hand, since the pointer to it is
  !-- global (because it is bound to C++) and cannot be nullified by the
  !-- destructor of bns. In case of multiple bns objects, this would lead
  !-- to problems...
  !-- TODO: fix this
  !
  !CALL binary% destruct_binary()


  CONTAINS


  SUBROUTINE cauchy_convergence_test_known( formul_dx, formul_dx2, &
                                            use_constraints, ref_lev )

    IMPLICIT NONE

    CLASS(tpo), INTENT( IN OUT ):: formul_dx, formul_dx2
    INTEGER, INTENT( IN ):: use_constraints, ref_lev

    INTEGER:: ix, iy, iz, nx, ny, nz, unit_cauchy_ct, unit_cauchy_parts_ct, &
              min_ix_y, min_iy_y, min_iz_y, &
              min_ix_z, min_iy_z, min_iz_z

    DOUBLE PRECISION:: min_abs_y, min_abs_z
    DOUBLE PRECISION, DIMENSION( :, :, :, : ), ALLOCATABLE:: abs_grid

    DOUBLE PRECISION, PARAMETER:: tiny_real= 1D-30
    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: grid_dx
    DOUBLE PRECISION, DIMENSION(3):: point_dx2
    DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE:: convergence_factor

    CHARACTER( LEN=: ), ALLOCATABLE:: name_cauchy_ct, name_cauchy_parts_ct

    LOGICAL:: exist

    nx= formul_dx% get_ngrid_x(ref_lev)
    ny= formul_dx% get_ngrid_y(ref_lev)
    nz= formul_dx% get_ngrid_z(ref_lev)

    nx= FLOOR( DBLE( nx - 1 )/denominator_ratio_dx ) + 1
    ny= FLOOR( DBLE( ny - 1 )/denominator_ratio_dx ) + 1
    nz= FLOOR( DBLE( nz - 1 )/denominator_ratio_dx ) + 1

    ALLOCATE( grid_dx( 3, nx, ny, nz ) )
    ALLOCATE( abs_grid( 3, nx, ny, nz ) )
    ALLOCATE( convergence_factor( nx, ny, nz ) )

    PRINT *, "** Computing convergence factor..."

    choose_constraints: SELECT CASE( use_constraints )

    CASE(1)

      shared_grid_loops1: DO iz= 0, nz - 1, 1
        DO iy= 0, ny - 1, 1
          DO ix= 0, nx - 1, 1

            grid_dx( :, 1 + ix, 1 + iy, 1 + iz ) = &
                        formul_dx% get_grid_point(  &
                                    1 + INT(denominator_ratio_dx)*ix, &
                                    1 + INT(denominator_ratio_dx)*iy, &
                                    1 + INT(denominator_ratio_dx)*iz, ref_lev )
                        point_dx2= formul_dx2% get_grid_point( &
                                    1 + INT(numerator_ratio_dx)*ix, &
                                    1 + INT(numerator_ratio_dx)*iy, &
                                    1 + INT(numerator_ratio_dx)*iz, ref_lev )

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
                                   1 + INT(denominator_ratio_dx)*iz, ref_lev ))&
            /( formul_dx2% get_HC( 1 + INT(numerator_ratio_dx)*ix, &
                                   1 + INT(numerator_ratio_dx)*iy, &
                                   1 + INT(numerator_ratio_dx)*iz, ref_lev ) &
             + tiny_real ) &
             ) + tiny_real)/LOG(ratio_dx)

          ENDDO
        ENDDO
      ENDDO shared_grid_loops1
      PRINT *, " * Convergence factor computed."
      PRINT *

      unit_cauchy_ct= 3109
      name_cauchy_ct= TRIM(spacetime_path)//"cauchy_convergence_test_known.dat"

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
      IF( ios > 0 )THEN
        PRINT *, "...error when opening ", TRIM(name_cauchy_ct), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when opening " &
      !         // TRIM(name_cauchy_ct) )

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
            IF( ios > 0 )THEN
              PRINT *, "...error when writing he arrays in ", &
                       TRIM(name_cauchy_ct), &
                       ". The error message is", err_msg
              STOP
            ENDIF
            !CALL test_status( ios, err_msg, "...error in writing " &
            !            // "the arrays in " // TRIM(name_cauchy_ct) )
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
                                      1 + INT(denominator_ratio_dx)*iz, ref_lev )
                        point_dx2= formul_dx2% get_grid_point( &
                                      1 + INT(numerator_ratio_dx)*ix, &
                                      1 + INT(numerator_ratio_dx)*iy, &
                                      1 + INT(numerator_ratio_dx)*iz, ref_lev )

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
                                         1 + INT(denominator_ratio_dx)*iz, ref_lev )&
            /( formul_dx2% get_HC_parts( 1 + INT(numerator_ratio_dx)*ix, &
                                         1 + INT(numerator_ratio_dx)*iy, &
                                         1 + INT(numerator_ratio_dx)*iz, ref_lev )  &
             + tiny_real ) &
             ) + tiny_real )/LOG(ratio_dx)

          ENDDO
        ENDDO
      ENDDO shared_grid_loops2
      PRINT *, " * Convergence factor computed."
      PRINT *

      unit_cauchy_parts_ct= 3111
      name_cauchy_parts_ct= TRIM(spacetime_path) &
                            // "cauchy_convergence_test_known_parts.dat"

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
      IF( ios > 0 )THEN
        PRINT *, "...error when opening ", TRIM(name_cauchy_parts_ct), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when opening " &
      !         // TRIM(name_cauchy_parts_ct) )

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
            IF( ios > 0 )THEN
              PRINT *, "...error when writing e arrays in ", &
                       TRIM(name_cauchy_parts_ct), &
                       ". The error message is", err_msg
              STOP
            ENDIF
            !CALL test_status( ios, err_msg, "...error in writing " &
            !            // "the arrays in " // TRIM(name_cauchy_parts_ct) )
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
                                              formul_dx4, use_constraints, &
                                              ref_lev )

    IMPLICIT NONE

    CLASS(tpo), INTENT( IN OUT ):: formul_dx, formul_dx2, formul_dx4
    INTEGER, INTENT( IN ):: use_constraints, ref_lev

    INTEGER:: ix, iy, iz, nx, ny, nz, unit_cauchy_ct, unit_cauchy_parts_ct, &
              min_ix_y, min_iy_y, min_iz_y, &
              min_ix_z, min_iy_z, min_iz_z

    DOUBLE PRECISION:: min_abs_y, min_abs_z
    DOUBLE PRECISION, DIMENSION( :, :, :, : ), ALLOCATABLE:: abs_grid

    DOUBLE PRECISION, PARAMETER:: tiny_real= 1D-30
    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: grid_dx
    DOUBLE PRECISION, DIMENSION(3):: point_dx2
    DOUBLE PRECISION, DIMENSION(3):: point_dx4
    DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE:: convergence_factor

    CHARACTER( LEN=: ), ALLOCATABLE:: name_cauchy_ct, name_cauchy_parts_ct

    LOGICAL:: exist

    nx= formul_dx% get_ngrid_x(ref_lev)
    ny= formul_dx% get_ngrid_y(ref_lev)
    nz= formul_dx% get_ngrid_z(ref_lev)

    nx= FLOOR( DBLE( nx - 1 )/denominator_ratio_dx**2 ) + 1
    ny= FLOOR( DBLE( ny - 1 )/denominator_ratio_dx**2 ) + 1
    nz= FLOOR( DBLE( nz - 1 )/denominator_ratio_dx**2 ) + 1

    ALLOCATE( grid_dx( 3, nx, ny, nz ) )
    ALLOCATE( abs_grid( 3, nx, ny, nz ) )
    ALLOCATE( convergence_factor( nx, ny, nz ) )

    PRINT *, "** Computing convergence factor..."

    IF( debug )THEN
      PRINT *, "formul_dx%  get_ngrid_x=", formul_dx%  get_ngrid_x(ref_lev)
      PRINT *, "formul_dx%  get_ngrid_y=", formul_dx%  get_ngrid_y(ref_lev)
      PRINT *, "formul_dx%  get_ngrid_z=", formul_dx%  get_ngrid_z(ref_lev)
      PRINT *, "formul_dx2% get_ngrid_x=", formul_dx2% get_ngrid_x(ref_lev)
      PRINT *, "formul_dx2% get_ngrid_y=", formul_dx2% get_ngrid_y(ref_lev)
      PRINT *, "formul_dx2% get_ngrid_z=", formul_dx2% get_ngrid_z(ref_lev)
      PRINT *, "formul_dx4% get_ngrid_x=", formul_dx4% get_ngrid_x(ref_lev)
      PRINT *, "formul_dx4% get_ngrid_y=", formul_dx4% get_ngrid_y(ref_lev)
      PRINT *, "formul_dx4% get_ngrid_z=", formul_dx4% get_ngrid_z(ref_lev)
      PRINT *
      PRINT *, "formul_dx%  get_dx=", formul_dx %  get_dx(ref_lev)
      PRINT *, "formul_dx2% get_dy=", formul_dx2%  get_dx(ref_lev)
      PRINT *, "formul_dx4% get_dz=", formul_dx4%  get_dx(ref_lev)
      PRINT *
      !STOP
    ENDIF

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
                                1 + INT(denominator_ratio_dx**2)*iz, ref_lev )
            point_dx2= formul_dx2% get_grid_point( &
                                1 + INT(numerator_ratio_dx* &
                                    denominator_ratio_dx)*ix, &
                                1 + INT(numerator_ratio_dx* &
                                    denominator_ratio_dx)*iy, &
                                1 + INT(numerator_ratio_dx* &
                                    denominator_ratio_dx)*iz, ref_lev )
            point_dx4= formul_dx4% get_grid_point( &
                                1 + INT(numerator_ratio_dx**2)*ix, &
                                1 + INT(numerator_ratio_dx**2)*iy, &
                                1 + INT(numerator_ratio_dx**2)*iz, ref_lev )

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
                               1 + INT(denominator_ratio_dx**2)*iz, ref_lev )) &
             - ABS(formul_dx2% get_HC( 1 + INT(numerator_ratio_dx* &
                                           denominator_ratio_dx)*ix, &
                                   1 + INT(numerator_ratio_dx* &
                                       denominator_ratio_dx)*iy, &
                                   1 + INT(numerator_ratio_dx* &
                                       denominator_ratio_dx)*iz, ref_lev ))) &
            /( ABS(formul_dx2% get_HC( 1 + INT(numerator_ratio_dx* &
                                       denominator_ratio_dx)*ix, &
                                   1 + INT(numerator_ratio_dx* &
                                       denominator_ratio_dx)*iy, &
                                   1 + INT(numerator_ratio_dx* &
                                       denominator_ratio_dx)*iz, ref_lev )) &
             - ABS(formul_dx4% get_HC( 1 + INT(numerator_ratio_dx**2)*ix, &
                               1 + INT(numerator_ratio_dx**2)*iy, &
                               1 + INT(numerator_ratio_dx**2)*iz, ref_lev ) ) &
             + tiny_real ) &
             ) + tiny_real&
             )/LOG(ratio_dx)

          ENDDO
        ENDDO
      ENDDO shared_grid_loops3
      PRINT *, " * Convergence factor computed."
      PRINT *

      unit_cauchy_ct= 3108
      name_cauchy_ct= TRIM(spacetime_path) &
                      //"cauchy_convergence_test_unknown.dat"

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
      IF( ios > 0 )THEN
        PRINT *, "...error when opening ", TRIM(name_cauchy_ct), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when opening " &
      !         // TRIM(name_cauchy_ct) )

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
            IF( ios > 0 )THEN
              PRINT *, "...error when writing the arrays in ", &
                       TRIM(name_cauchy_ct), &
                       ". The error message is", err_msg
              STOP
            ENDIF
            !CALL test_status( ios, err_msg, "...error in writing " &
            !            // "the arrays in " // TRIM(name_cauchy_ct) )
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
                                1 + INT(denominator_ratio_dx**2)*iz, ref_lev )
          point_dx2= formul_dx2% get_grid_point( &
                                1 + INT(numerator_ratio_dx* &
                                    denominator_ratio_dx)*ix, &
                                1 + INT(numerator_ratio_dx* &
                                    denominator_ratio_dx)*iy, &
                                1 + INT(numerator_ratio_dx* &
                                    denominator_ratio_dx)*iz, ref_lev )
          point_dx4= formul_dx4% get_grid_point( &
                                1 + INT(numerator_ratio_dx**2)*ix, &
                                1 + INT(numerator_ratio_dx**2)*iy, &
                                1 + INT(numerator_ratio_dx**2)*iz, ref_lev )

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
                                1 + INT(denominator_ratio_dx**2)*iz, ref_lev )) &
           - ABS(formul_dx2% get_HC_parts( &
                                1 + INT(numerator_ratio_dx* &
                                        denominator_ratio_dx)*ix, &
                                1 + INT(numerator_ratio_dx* &
                                        denominator_ratio_dx)*iy, &
                                1 + INT(numerator_ratio_dx* &
                                        denominator_ratio_dx)*iz, ref_lev )))&
          /( ABS(formul_dx2% get_HC_parts( &
                                1 + INT(numerator_ratio_dx* &
                                        denominator_ratio_dx)*ix, &
                                1 + INT(numerator_ratio_dx* &
                                        denominator_ratio_dx)*iy, &
                                1 + INT(numerator_ratio_dx* &
                                        denominator_ratio_dx)*iz, ref_lev )) &
           - ABS(formul_dx4% get_HC_parts( &
                                1 + INT(numerator_ratio_dx**2)*ix, &
                                1 + INT(numerator_ratio_dx**2)*iy, &
                                1 + INT(numerator_ratio_dx**2)*iz, ref_lev ) )&
           + tiny_real ) &
           ) + tiny_real &
           )/LOG(ratio_dx)

          ENDDO
        ENDDO
      ENDDO shared_grid_loops4
      PRINT *, " * Convergence factor computed."
      PRINT *

      unit_cauchy_parts_ct= 3110
      name_cauchy_parts_ct= TRIM(spacetime_path) &
                            //"cauchy_convergence_test_unknown_parts.dat"

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
      IF( ios > 0 )THEN
        PRINT *, "...error when opening ", TRIM(name_cauchy_parts_ct), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when opening " &
      !         // TRIM(name_cauchy_parts_ct) )

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
            IF( ios > 0 )THEN
              PRINT *, "...error when writing the arrays in ", &
                       TRIM(name_cauchy_parts_ct), &
                       ". The error message is", err_msg
              STOP
            ENDIF
            !CALL test_status( ios, err_msg, "...error in writing " &
            !            // "the arrays in " // TRIM(name_cauchy_parts_ct) )
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

    lorene_bns_id_parameters= 'sphincs_id_parameters.dat'

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
