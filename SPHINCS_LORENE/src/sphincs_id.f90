! File:         sphincs_id.f90
! Author:       Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

!************************************************************************
! Copyright (C) 2022 Francesco Torsello                                 *
!                                                                       *
! This program is free software: you can redistribute it and/or modify  *
! it under the terms of the GNU General Public License as published by  *
! the Free Software Foundation, either version 3 of the License, or     *
! (at your option) any later version.                                   *
!                                                                       *
! This program is distributed in the hope that it will be useful,       *
! but WITHOUT ANY WARRANTY; without even the implied warranty of        *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          *
! GNU General Public License for more details.                          *
!                                                                       *
! You should have received a copy of the GNU General Public License     *
! along with this program. If not, see <https://www.gnu.org/licenses/>. *
!************************************************************************

PROGRAM sphincs_id

  !*****************************************************
  !
  !# Use the MODULE sphincs_lorene to export binary
  !  files containing the initial data (|id|) required
  !  by the evolution code in SPHINCS, and built using
  !  the binary files produced by LORENE and containing
  !  the binary neutron stars (BNS) |id|.
  !
  !  FT 28.10.2020
  !
  !*****************************************************

#ifdef __INTEL_COMPILER
  USE IFPORT,          ONLY: MAKEDIRQQ
#endif

#if flavour == 1
  USE sphincs_id_lorene,  ONLY: allocate_idbase
#elif flavour == 2
  USE sphincs_id_interpolate,  ONLY: allocate_idbase
#endif

  USE utility,         ONLY: date, time, zone, values, run_id, itr, itr3, &
                             itr4, file_exists, cnt, &
                             test_status, show_progress, end_time
  USE timing,          ONLY: timer
  USE id_base,         ONLY: idbase, initialize
  USE particles_id,    ONLY: particles
  USE formul_bssn_id,  ONLY: bssn_id
  USE constants,       ONLY: lorene2hydrobase, c_light2, k_lorene2hydrobase, &
                             k_lorene2hydrobase_piecewisepolytrope, MSun_geo, &
                             kg2g, m2cm, m0c2

  IMPLICIT NONE


  INTEGER, PARAMETER:: max_length= 50
  !! Maximum length for strings
  INTEGER, PARAMETER:: max_n_bns= 50
  ! Maximum number of physical systems
  INTEGER, PARAMETER:: max_n_parts= 250
  !! Maximum number of particle distributions

  INTEGER, PARAMETER:: test_int= - 112
  INTEGER, DIMENSION( max_n_bns, max_n_parts ):: placer= test_int
  !# Matrix storing the information on how to place particles for each bns
  !  object. Row i contains information about the i^th bns object.

  INTEGER:: n_bns
  !! Number of physical systems to set up
  INTEGER:: i_matter
  !! Index running over the number of physical systems
  INTEGER:: ref_lev
  !! Number of refinement levels
  INTEGER:: constraints_step
  !! Export the constraints every constraints_step-th step


  DOUBLE PRECISION:: numerator_ratio_dx
  !# Numerator of the rational ratio between the large grid spacing and the
  !  medium one,equal to the ratio between the medium grid spacing nd the small
  !  one. Not used in this PROGRAM, but needed since the PROGRAM reads the same
  !  parameter file as the convergence_test PROGRAM
  DOUBLE PRECISION:: denominator_ratio_dx
  !# Denominator of the rational ratio between the large grid spacing and the
  !  medium one,equal to the ratio between the medium grid spacing nd the small
  !  one. Not used in this PROGRAM, but needed since the PROGRAM reads the same
  !  parameter file as the convergence_test PROGRAM


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
  ! containing the LORENE |id| for different BNS
  !TYPE( bnslorene ), DIMENSION(:), ALLOCATABLE:: binaries

  TYPE id
    CLASS( idbase ), ALLOCATABLE:: idata
  END TYPE id
  TYPE( id ), DIMENSION(:), ALLOCATABLE:: ids
  !CLASS(idbase), POINTER:: foo
  !CLASS( idbase ), DIMENSION(:), ALLOCATABLE:: ids
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

  PRINT *, "  ____________________________________________________ "
  PRINT *, "       ____________  ________  __________    __ ___    "
  PRINT *, "      / ___/ _  / /_/ / / __ \/ ___/ ___/   / / __ \   "
  PRINT *, "     (__  ) ___/ __  / / / / / /__(__  )___/ / /_/ /   "
  PRINT *, "    /____/_/  /_/ /_/_/_/ /_/____/____/___/_/_____/    "
  PRINT *
  PRINT *, "  Smoothed Particle Hydrodynamics IN Curved Spacetime  "
  PRINT *, "  Initial Data builder  "
  PRINT *
  PRINT *, "  Author:       Francesco Torsello  "
  PRINT *, "  Copyright:    GNU General Public License (GPLv3)     "
  PRINT *
  PRINT *, "  SPHINCS_ID  Copyright (C) 2022  Francesco Torsello   "
  PRINT *
  PRINT *, "  This program is free software: you can redistribute  "
  PRINT *, "  it and/or modify it under the terms of the GNU       "
  PRINT *, "  General Public License as published by the Free      "
  PRINT *, "  Software Foundation, either version 3 of the License,"
  PRINT *, "  or (at your option) any later version.               "
  PRINT *
  PRINT *, "  This program is distributed in the hope that it will "
  PRINT *, "  be useful, but WITHOUT ANY WARRANTY; without even   "
  PRINT *, "  the implied warranty of MERCHANTABILITY or FITNESS   "
  PRINT *, "  FOR A PARTICULAR PURPOSE. See the GNU General Public "
  PRINT *, "  License for more details.                            "
  PRINT *
  PRINT *, "  You should have received a copy of the GNU General   "
  PRINT *, "  Public License along with this program. If not, see  "
  PRINT *, "  https://www.gnu.org/licenses/.  "
  PRINT *, "  ____________________________________________________ "
  PRINT *
  PRINT *, "  Run id: ", run_id
  PRINT *, "  ____________________________________________________ "
  PRINT *

  execution_timer= timer( "execution_timer" )
  CALL execution_timer% start_timer()

  CALL read_bns_id_parameters()

  !
  !-- Check that the specified subdirectories exist. If not, create them
  !-- TODO: this compils with ifort, but not with gfortran
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

  ALLOCATE( CHARACTER(5):: systems(n_bns) )
  ALLOCATE( CHARACTER(5):: systems_name(n_bns) )

  !DO itr= 1, n_bns, 1
  !  systems(itr)= filenames(itr)(1:5)
  !  IF( systems(itr) /= bnslo .AND. systems(itr) /= drslo )THEN
  !    PRINT *, "** ERROR! Unrecognized physical system ", systems(itr), ",",&
  !             " system number", itr, "."
  !    PRINT *
  !    PRINT *, "   Please specify the type of physical system in the first 5",&
  !             " characters of the name of the file containing the initial", &
  !             " data."
  !    PRINT *
  !    PRINT *, "   The 5-character names, and associated physical systems,", &
  !             " supported by this version of SPHINCS_ID are:"
  !    PRINT *
  !    PRINT *, "   1. BNSLO: Binary Neutron Stars produced with LORENE"
  !    PRINT *, "   2. DRSLO: Differentially Rotating Star produced with LORENE"
  !    PRINT *
  !    STOP
  !  ENDIF
  !ENDDO

  ! Allocate needed memory
  !ALLOCATE( binaries      ( n_bns ) )
  !ALLOCATE( diffrotstars  ( n_bns ) )
  !IF( TRIM(systems(1)) == bnslo )THEN
  !  ALLOCATE( bnslorene:: idata( n_bns ) )
  !ELSEIF( TRIM(systems(1)) == drslo )THEN
  !  ALLOCATE( diffstarlorene:: idata( n_bns ) )
  !ELSE
  !  PRINT *, "** ERROR! Unknown name for the physical system: ", TRIM(systems(1))
  !  PRINT *, "   Set the variable 'system' in the parameter file ", &
  !           "sphincs_lorene_parameters.par to one of the values listed there."
  !  PRINT *, "   Stopping..."
  !  PRINT *
  !  STOP
  !ENDIF
  ALLOCATE( ids( n_bns ) )
  ALLOCATE( particles_dist( n_bns, max_n_parts ) )
  ALLOCATE( bssn_forms    ( n_bns ) )

  !
  !-- Construct the LORENE |id| from the LORENE binary files
  !
 ! build_bns_loop: DO itr= 1, n_bns, 1
 !   binaries( itr )= bnslorene( TRIM(common_path)//TRIM(filenames( itr )) )
 !   ! Set the variables to decide on using the geodesic gauge or not
 !   ! (lapse=1, shift=0)
 !   binaries( itr )% one_lapse = one_lapse
 !   binaries( itr )% zero_shift= zero_shift
 ! ENDDO build_bns_loop

  !ALLOCATE(foo, source = bnslorene( TRIM(common_path)//TRIM(filenames( 1 ))) )

  build_drs_loop: DO itr= 1, n_bns, 1
   ! IF( systems(itr) == bnslo )THEN
   !   !idata( itr )=
   !   !ALLOCATE(foo, source = bnslorene( TRIM(common_path)//TRIM(filenames( 1 ))) )
   !   ALLOCATE( bnslorene:: ids(itr)% idata )
   !   !CALL ids(itr)% idata% initialize( TRIM(common_path)//TRIM(filenames(itr)))
   ! ELSEIF( systems(itr) == drslo )THEN
   !   !idata( itr )=
   !   ALLOCATE( diffstarlorene:: ids(itr)% idata )
   !   !CALL ids(itr)% idata% initialize( TRIM(common_path)//TRIM(filenames(itr)))
   !   !ids(itr)% idata = diffstarlorene( TRIM(common_path)//TRIM(filenames( itr )) )
   ! ELSE
   !   PRINT *, "** ERROR! Unknown name for the physical system: ", systems(itr)
   !   PRINT *, "   Set the variable 'system' in the parameter file ", &
   !            "sphincs_lorene_parameters.par to one of the values listed there."
   !   PRINT *, "   Stopping..."
   !   PRINT *
   !   STOP
   ! ENDIF
    CALL allocate_idbase( ids(itr)% idata, TRIM(filenames(itr)), &
                          systems(itr), systems_name(itr) )
    CALL ids(itr)% idata% initialize( TRIM(common_path)//TRIM(filenames(itr)) )

    !ids(itr)% idata= initialize( ids(itr)% idata, TRIM(common_path)//TRIM(filenames(itr)) )
    !ids(itr)% idata= ids(1)% idata_constructor( file )
    ! Set the variables to decide on using the geodesic gauge or not
    ! (lapse=1, shift=0)
    !CALL idata( itr )% set_one_lapse( one_lapse )
    !CALL idata( itr )% set_zero_shift( zero_shift )
    CALL ids(itr)% idata% set_one_lapse( one_lapse )
    CALL ids(itr)% idata% set_zero_shift( zero_shift )
  ENDDO build_drs_loop

  !PRINT *, ids(1)% idata% read_mass_density( 2.0D0, 3.0D0, 5.0D0 )

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
          PRINT *, " Placing particles for "//systems(itr3), itr3, &
                   ", distribution", itr4
          PRINT *, "===================================================" &
                   // "==============="
          PRINT *
          particles_dist( itr3, itr4 )= particles( ids(itr3)% idata, &
                                                   placer( itr3, itr4 ) )

        ENDIF
      ENDDO part_distribution_loop
    ENDDO place_hydro_id_loops

    !STOP

    compute_export_sph_loops: DO itr3= 1, n_bns, 1
      part_distribution_loop2: DO itr4= 1, max_n_parts, 1
        IF( placer( itr3, itr4 ) == test_int )THEN
          EXIT part_distribution_loop2
          ! Experimental: empty particles object
          !particles_dist( itr, itr2 )= particles()
        ELSE

          PRINT *, "===================================================" &
                   // "====================="
          PRINT *, " Computing SPH variables for "//systems(itr3), itr3, &
                   ", distribution", itr4
          PRINT *, "===================================================" &
                   // "====================="
          PRINT *
          !WRITE( namefile_parts_bin, "(A1,I1,A1,I1,A1)" ) &
          !                            "l", &
          !                            itr3, "-", itr4, "."
          WRITE( namefile_parts_bin, "(A5)" ) systems_name(itr3)
          namefile_parts_bin= TRIM( sph_path ) // TRIM( namefile_parts_bin )

          particles_dist( itr3, itr4 )% export_bin    = export_bin
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

  ENDIF

  IF( run_spacetime )THEN

    !
    !-- Construct the bssn_id objects from the bns objects
    !
    place_spacetime_id_loop: DO itr3 = 1, n_bns, 1
      PRINT *, "===================================================" &
               // "==============="
      PRINT *, " Setting up BSSN object for "//systems(itr3), itr3
      PRINT *, "===================================================" &
               // "==============="
      PRINT *
      bssn_forms( itr3 )= bssn_id( ids(itr3)% idata )
    ENDDO place_spacetime_id_loop

    !
    !-- Compute the BSSN initial data, optionally export it to a binary file
    !-- readable by SPHINCS_BSSN, and optionally read the content of such binary
    !-- file and print it to a formatted file (the latter for debugging)
    !
    compute_export_bssn_loop: DO itr3 = 1, n_bns, 1
      PRINT *, "===================================================" &
               // "==============="
      PRINT *, " Computing BSSN variables for "//systems(itr3), itr3
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
                      compute_and_export_3p1_constraints( ids(itr3)% idata, &
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

  DO itr= 1, n_bns, 1

    PRINT *, "===================================================" &
             // "================================================"
    PRINT *, " Timing for physical system ", itr
    PRINT *, "===================================================" &
             // "================================================"
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
  DO itr= 1, n_bns, 1

    PRINT *, "===================================================" &
             // "================================================"
    PRINT *, " Summary for physical system ", itr
    PRINT *, "===================================================" &
             // "================================================"
    PRINT *
    PRINT *, "   Used ID data file: " &
             // TRIM(common_path)//TRIM(filenames(itr))
    PRINT *

    CALL ids(itr)% idata% print_summary()

    IF( run_sph )THEN

      CALL particles_dist(itr,1)% print_summary()

    ENDIF

    IF( run_spacetime )THEN

      CALL bssn_forms(itr)% print_summary()

    ENDIF

  ENDDO
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
    !CALL binaries( itr )% destruct_binary()
  ENDDO
  !IF( ALLOCATED( binaries ) )THEN
  !  DEALLOCATE( binaries )
  !ENDIF
  !IF( ALLOCATED( idata ) )THEN
  !  DEALLOCATE( idata )
  !ENDIF
  IF( ALLOCATED( ids ) )THEN
    DEALLOCATE( ids )
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

END PROGRAM sphincs_id
