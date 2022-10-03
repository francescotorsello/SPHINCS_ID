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
  !# Make a Cauchy convergence test to check the
  !  validity of the |id| set up by |sphincsid|.
  !
  !  FT 8.12.2020
  !
  !*****************************************************

#ifdef __INTEL_COMPILER

  USE IFPORT,          ONLY: MAKEDIRQQ

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

  USE id_base,                  ONLY: idbase, initialize
  USE sph_particles,            ONLY: particles
  USE bssn_formulation,         ONLY: bssn
  USE standard_tpo_formulation, ONLY: tpo
  USE tensor,                   ONLY: jy, jz
  USE timing,                   ONLY: timer
  USE utility,                  ONLY: date, time, zone, values, run_id, &
                                      itr3, &
                                      test_status, show_progress, end_time, &
                                      read_sphincs_id_parameters, one, &
                                      !----------
                                      n_id, common_path, filenames, placer, &
                                      export_bin, export_form, export_form_xy, &
                                      export_form_x, export_constraints_xy, &
                                      export_constraints_x, &
                                      compute_constraints, &
                                      export_constraints, &
                                      export_constraints_details, &
                                      constraints_step, &
                                      compute_parts_constraints, &
                                      numerator_ratio_dx, &
                                      denominator_ratio_dx, &
                                      one_lapse, zero_shift, show_progress, &
                                      run_sph, run_spacetime, sph_path, &
                                      spacetime_path, estimate_length_scale, &
                                      test_int, max_n_parts, ref_lev
  USE cauchy_convergence_test,  ONLY: find_shared_grid, &
                                      perform_cauchy_convergence_test, &
                                      use_constraints_on_mesh, &
                                      use_constraints_with_mapped_hydro
  USE ISO_FORTRAN_ENV,          ONLY: COMPILER_VERSION, COMPILER_OPTIONS

  IMPLICIT NONE

  ! Loop limits for BSSN objects (for debugging; 3 is for production)
  INTEGER, PARAMETER:: min_bssn= 1
  INTEGER, PARAMETER:: max_bssn= 3

  DOUBLE PRECISION, PARAMETER:: tol_dx= 1.D-10
  DOUBLE PRECISION, PARAMETER:: tol_coord= 1.D-10

  ! Grid spacing for the first BSSN object; the other two will have
  ! original_dx/2 and original_dx/4 as grid spacings
  DOUBLE PRECISION:: original_dx
  DOUBLE PRECISION:: ratio_dx

  DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: shared_grid

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
  CHARACTER( LEN= 500 ):: hostname
  !# String storing the name of the host machine

  LOGICAL, PARAMETER:: debug= .FALSE.
  LOGICAL:: exist
  LOGICAL(4):: dir_out

  CLASS(idbase), ALLOCATABLE:: idata

  TYPE(particles):: particles_dist
  !# Array storing the particles objects,
  !  containing the particle distributions for each idbase object.
  !  Multiple particle objects can contain different particle distributions
  !  for the same idbase object.

  TYPE(bssn), DIMENSION(3):: bssn_forms
  !# Array storing the bssn objects,
  !  containing the BSSN variables on the gravity grid for each idbase object

  TYPE(timer):: execution_timer

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

  !CALL HOSTNM( hostname )

  PRINT *, "  ________________________________________________________________ "
  PRINT *, "             ____________  ________  __________    __ ___          "
  PRINT *, "            / ___/ _  / /_/ / / __ \/ ___/ ___/   / / __ \         "
  PRINT *, "           (__  ) ___/ __  / / / / / /__(__  )___/ / /_/ /         "
  PRINT *, "          /____/_/  /_/ /_/_/_/ /_/____/____/___/_/_____/          "
  PRINT *
  PRINT *, "  Smoothed Particle Hydrodynamics IN Curved Spacetime              "
  PRINT *, "  Initial Data builder, v1.6 - Cauchy convergence test             "
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
  !PRINT *
  !PRINT *, "  This run was done on: "
  !PRINT *, TRIM(hostname)
  PRINT *, "  ________________________________________________________________ "
  PRINT *
  PRINT *, "  Run id: ", run_id
  PRINT *, "  ________________________________________________________________ "
  PRINT *

  execution_timer= timer( "execution_timer" )
  CALL execution_timer% start_timer()

  CALL read_sphincs_id_parameters()

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
  PRINT *, "===================================================" &
           // "==============="
  PRINT *, " Constructing idbase object for "//systems(1)
  PRINT *, "===================================================" &
           // "==============="
  PRINT *
  CALL idata% initialize( TRIM(common_path)//TRIM(filenames(1)) )
  ! Set the variables to decide on using the geodesic gauge or not
  ! (lapse=1, shift=0)
  CALL idata% set_one_lapse ( one_lapse )
  CALL idata% set_zero_shift( zero_shift )


  !
  !-- Construct the bssn objects from the idbase object
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

      IF( ABS( bssn_forms( itr3 )% get_dx(ref_lev) - &
          original_dx/( ratio_dx**( itr3 - 1 ) ) ) &
          /(original_dx/( ratio_dx**( itr3 - 1 ) )) > tol_dx &
      )THEN
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

  ! Find the grid points shared by the grids
  CALL find_shared_grid( bssn_forms(min_bssn), bssn_forms(min_bssn + 1), &
                         bssn_forms(min_bssn + 2), &
                         numerator_ratio_dx, denominator_ratio_dx, ref_lev, &
                         shared_grid )

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
    WRITE( namefile_bssn_bin, "(A10,I1,A4)" ) "BSSN_vars-", itr3, ".bin"
    namefile_bssn_bin= TRIM( spacetime_path ) // TRIM( namefile_bssn_bin )

    bssn_forms( itr3 )% export_bin= export_bin
    bssn_forms( itr3 )% export_form_xy= export_form_xy
    bssn_forms( itr3 )% export_form_x = export_form_x
    CALL bssn_forms( itr3 )% &
                        compute_and_print_tpo_variables( namefile_bssn_bin )
  ENDDO compute_export_bssn_loop

  !
  !-- Print the BSSN initial data to a formatted file
  !
  IF( export_form )THEN
    export_bssn_loop: DO itr3 = min_bssn, max_bssn, 1
      WRITE( namefile_bssn, "(A10,I1,A4)" ) &
                            "bssn-form_", itr3, ".dat"
      namefile_bssn= TRIM( spacetime_path ) // TRIM( namefile_bssn )

      CALL bssn_forms( itr3 )% &
                  print_formatted_id_tpo_variables( namefile= namefile_bssn )
    ENDDO export_bssn_loop
  ENDIF

  !
  !-- Construct the particles object from the idbase object
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
    CALL particles_dist% compute_and_print_sph_variables( namefile_parts )

    !
    !-- Print the particle initial data to a formatted file
    !
    IF( export_form )THEN
      WRITE( namefile_parts, "(A34)" ) &
                             "sph-form_1.dat"
      namefile_parts= TRIM( sph_path ) // TRIM( namefile_parts )
      CALL particles_dist% print_formatted_id_particles( namefile_parts )
    ENDIF

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
                  compute_and_print_tpo_constraints( idata, &
                                                     namefile_bssn, &
                                                     name_logfile, &
                                                     shared_grid )

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
                  compute_and_print_tpo_constraints( particles_dist, &
                                                     namefile_bssn, &
                                                     name_logfile, &
                                                     shared_grid )

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

    CALL perform_cauchy_convergence_test &
      ( bssn_forms(1), bssn_forms(2), bssn_forms(3), use_constraints_on_mesh, &
        numerator_ratio_dx, denominator_ratio_dx, ref_lev )

    CALL perform_cauchy_convergence_test &
      ( bssn_forms(2), bssn_forms(3), use_constraints_on_mesh, &
        numerator_ratio_dx, denominator_ratio_dx, ref_lev )

  ENDIF

  IF( compute_parts_constraints )THEN

    PRINT *, "** Performing convergence test with constraints computed ", &
             "with particle data."
    PRINT *

    CALL perform_cauchy_convergence_test &
      ( bssn_forms(1), bssn_forms(2), bssn_forms(3), &
        use_constraints_with_mapped_hydro, &
        numerator_ratio_dx, denominator_ratio_dx, ref_lev )

    CALL perform_cauchy_convergence_test &
      ( bssn_forms(2), bssn_forms(3), &
        use_constraints_with_mapped_hydro, &
        numerator_ratio_dx, denominator_ratio_dx, ref_lev )

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
  !-- destructor of bns. In case of multiple idbase objects, this would lead
  !-- to problems...
  !-- TODO: fix this
  !
  !CALL binary% destruct_binary()


END PROGRAM convergence_test
