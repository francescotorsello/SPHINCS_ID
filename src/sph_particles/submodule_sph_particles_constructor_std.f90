! File:         submodule_sph_particles_constructor_std.f90
! Authors:      Francesco Torsello (FT)
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

SUBMODULE (sph_particles) constructor_std

  !************************************************
  !
  !# This SUBMODULE contains the implementation
  !  of the constructor and the
  !  destructor of TYPE sph_particles.
  !
  !  FT 16.10.2020
  !
  !************************************************


  IMPLICIT NONE


  TYPE parts_i

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos_i
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pvol_i
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: h_i
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nu_i
!    TYPE(eos):: eos_i


!    CONTAINS
!
!
!    PROCEDURE:: compute_pressure

  END TYPE parts_i

!  INTERFACE
!    MODULE SUBROUTINE compute_pressure( this, parts, npart, x, y, z, nlrf, &
!      pressure )
!      CLASS(parts_i),   INTENT(IN)   :: this
!      CLASS(particles), INTENT(INOUT):: parts
!      INTEGER,          INTENT(IN)   :: npart
!      !! Returns the baryon mass density at the desired point
!      DOUBLE PRECISION, INTENT(IN)   :: x(npart)
!      !! \(x\) coordinate of the desired point
!      DOUBLE PRECISION, INTENT(IN)   :: y(npart)
!      !! \(y\) coordinate of the desired point
!      DOUBLE PRECISION, INTENT(IN)   :: z(npart)
!      !! \(z\) coordinate of the desired point
!      DOUBLE PRECISION, INTENT(IN)   :: nlrf(npart)
!      DOUBLE PRECISION, INTENT(INOUT):: pressure(npart)
!      !! Pressure at \((x,y,z)\)
!    END SUBROUTINE compute_pressure
!  END INTERFACE
!
!
!  PROCEDURE(), POINTER:: get_pressure_i


  CONTAINS


 ! MODULE PROCEDURE compute_pressure
 ! !! Compute the pressure from the given input
 !
 !   IMPLICIT NONE
 !
 !   DOUBLE PRECISION, DIMENSION(npart):: tmp, tmp2, tmp3
 !
 !   CALL parts% compute_sph_hydro(1, npart, &
 !     this% eos_i, nlrf, tmp, pressure, tmp2, tmp3 )
 !
 ! END PROCEDURE compute_pressure


  !MODULE PROCEDURE construct_particles_idase_empty
  !
  !    !************************************************
  !    !
  !    !# The constructor of an empty particle object.
  !    !
  !    !  FT 02.11.2020
  !    !
  !    !************************************************
  !
  !
  !    IMPLICIT NONE
  !
  !
  !    parts% empty_object= .TRUE.
  !
  !    parts% npart_temp= 0
  !
  !END PROCEDURE construct_particles_idase_empt


  MODULE PROCEDURE construct_particles_std

    !**************************************************
    !
    !# The constructor of TYPE particles is supposed
    !  to set up aparticle distribution by assigning
    !  the particle positions, their baryon numbers
    !  nu and first guesses for their smoothing lengths h.
    !  It also sets up the unit system and the kernel.
    !
    !  After the particle distribution is set up,
    !  it assigns the |id| to the particles.
    !  It does NOT compute the |sph| variables and it
    !  does NOT set up the neighbors' tree. The latter
    !  two things are delegated to the specific methods
    !  of TYPE particles that need them.
    !
    !  FT 17.10.2020
    !
    !  @note Last updated: FT 27.10.2022
    !
    !**************************************************

    USE constants,          ONLY: amu, Msun, pi, half, third
    USE utility,            ONLY: zero, one, two, ten, Msun_geo, flag$sph
    USE kernel_table,       ONLY: ktable
    USE input_output,       ONLY: read_options
    USE units,              ONLY: set_units
    USE options,            ONLY: ikernel, ndes
    USE alive_flag,         ONLY: alive
    USE analyze,            ONLY: COM
    USE utility,            ONLY: spherical_from_cartesian, &
                                  spatial_vector_norm_sym3x3, sph_path, &
                                  scan_1d_array_for_nans, eos$tabu$compose


    IMPLICIT NONE


    INTEGER, PARAMETER:: unit_pos= 2289
    INTEGER, PARAMETER:: unit_pos_out= 8754
    DOUBLE PRECISION, PARAMETER:: tol_equal_mass= 5.0D-3
    ! Tolerance for the difference between the masse of the stars
    ! for a BNS system, to determine if a BNS is equal-mass or not.
    ! If the relative difference between the masses of the stars is lower
    ! than tol_equal_mass, then the BNS is considered equal_mass.

    ! The variable counter counts how many times the PROCEDURE
    ! construct_particles_idase is called
    INTEGER, SAVE:: counter= 1
    INTEGER:: npart_des, a, max_steps, nx_gh, ny_gh, nz_gh, i_matter, &
              nlines, npart_tmp, n_matter_tmp

    ! Maximum length for strings, and for the number of imported binaries
    INTEGER, PARAMETER:: max_length= 50
    ! APM parameters
    INTEGER:: apm_max_it, max_inc, print_step
    ! Array storing the columns of the file parts_pos (defined below) that
    ! contain the particle positions
    INTEGER, DIMENSION(3):: columns
    INTEGER:: column_nu, header_lines, n_cols
    INTEGER:: tmp
    INTEGER, DIMENSION(id% get_n_matter()):: npart_des_i
    ! Temporary array storing the number of particles on each matter object
    INTEGER, DIMENSION(:), ALLOCATABLE:: npart_i_tmp

    DOUBLE PRECISION:: thres, nu_ratio_des, ghost_dist
    DOUBLE PRECISION:: xmin, xmax, ymin, ymax, zmin, zmax, stretch
    DOUBLE PRECISION:: upper_bound, lower_bound, upper_factor, lower_factor, &
                       last_r
    !DOUBLE PRECISION:: pvol_tmp
    DOUBLE PRECISION:: max_mass, total_mass
    !DOUBLE PRECISION:: ratio_npart_des_real, pmass_des
    DOUBLE PRECISION:: min_eps, min_vel, theta_a, phi_a, r_a!, rad_part

    DOUBLE PRECISION, DIMENSION(id% get_n_matter())  :: central_density
    DOUBLE PRECISION, DIMENSION(id% get_n_matter(),3):: center
    DOUBLE PRECISION, DIMENSION(id% get_n_matter(),3):: barycenter
    DOUBLE PRECISION, DIMENSION(id% get_n_matter(),6):: sizes
    DOUBLE PRECISION, DIMENSION(id% get_n_matter())  :: ghost_dists
    DOUBLE PRECISION, DIMENSION(id% get_n_matter())  :: lapse_lengthscales
    DOUBLE PRECISION, DIMENSION(id% get_n_matter())  :: g00_lengthscales

    DOUBLE PRECISION:: nuratio_thres, nuratio_des
    DOUBLE PRECISION:: min_lapse, min_g00_abs, shift_norm
    !DOUBLE PRECISION:: com_x, com_y, com_z, com_d

    TYPE(parts_i), DIMENSION(id% get_n_matter()):: parts_all

    ! String storing the name of the directory storing the files containing
    ! the particle distributions
    CHARACTER(LEN= max_length):: parts_pos_path
    ! String storing the name of the file containing the particle positions
    CHARACTER(LEN= max_length):: parts_pos
    ! Final name for the file containing the particle positions
    CHARACTER(LEN=:), ALLOCATABLE:: parts_pos_namefile

    CHARACTER(LEN=:), ALLOCATABLE:: parts_out_namefile
    !! Name for the file to print the final particle distribution and nu
    CHARACTER(LEN=3):: str_i
    ! String storing the local path to the directory where the
    ! |lorene| BNS ID files are stored
    CHARACTER(LEN= max_length):: compose_path
    ! String storing the names of the |lorene| BNS ID binary files
    CHARACTER(LEN= max_length):: compose_filename

    CHARACTER(LEN= max_length):: filename_apm_pos_id, filename_apm_pos, &
                                 filename_apm_results

    CHARACTER(LEN= max_length):: filename_mass_profile, &
                                 filename_shells_radii, filename_shells_pos

    LOGICAL:: file_exists, use_thres, redistribute_nu, correct_nu, &
              compose_eos, exist, randomize_phi, randomize_theta, &
              randomize_r, mass_it, adapt_ghosts, move_away_ghosts, &
              read_nu, reflect_particles_x, use_pressure

    LOGICAL, PARAMETER:: debug= .FALSE.

    LOGICAL, DIMENSION(max_length):: apm_iterate, use_atmosphere, &
                                     remove_atmosphere

    !TYPE procedure_pointer
    !  PROCEDURE, POINTER:: proc
    !END TYPE procedure_pointer

    !TYPE(procedure_pointer), DIMENSION(:), ALLOCATABLE:: compute_pressure_i

    ! Get the number of matter objects in the physical system
    parts% n_matter= id% get_n_matter()

    ! Get the the logical variable at specifies if the system is cold
    ! (no thermal component)
    parts% cold_system= id% get_cold_system()

    !
    !-- Initialize the timers
    !
    ALLOCATE( parts% apm_timers(parts% n_matter) )
    parts% placer_timer       = timer( "placer_timer" )
    parts% importer_timer     = timer( "importer_timer" )
    parts% sph_computer_timer = timer( "sph_computer_timer" )
    parts% same_particle_timer= timer( "same_particle_timer" )
    DO i_matter= 1, parts% n_matter, 1
      IF( parts% n_matter <= 9 ) WRITE( str_i, "(I1)" ) i_matter
      IF( parts% n_matter >= 10 .AND. parts% n_matter <= 99 ) &
                                                WRITE( str_i, "(I2)" ) i_matter
      IF( parts% n_matter >= 100 .AND. parts% n_matter <= 999 ) &
                                                WRITE( str_i, "(I3)" ) i_matter
      parts% apm_timers(i_matter)= timer( "apm_timer"//TRIM(str_i) )
    ENDDO

    ! Declare this object as non-empty (experimental)
    parts% empty_object= .FALSE.

    !
    !-- Read the options and parameters for the particle distributions
    !
    CALL read_particles_options()

    !
    !-- Read needed data from the idbase object
    !
    parts% nbar_tot       = zero
    parts% npart          = 0
    parts% distribution_id= dist

    ALLOCATE( parts% masses (parts% n_matter) )
    ALLOCATE( parts% all_eos(parts% n_matter) )
    ALLOCATE( parts% npart_i(0:parts% n_matter) )
    ALLOCATE( npart_i_tmp(0:parts% n_matter) )
    ALLOCATE( parts% nbar_i(parts% n_matter) )
    ALLOCATE( parts% nuratio_i(parts% n_matter) )
    ALLOCATE( parts% mass_ratios(parts% n_matter) )
    ALLOCATE( parts% mass_fractions(parts% n_matter) )
    ALLOCATE( parts% barycenter(parts% n_matter,3) )
    ALLOCATE( parts% surfaces (parts% n_matter) )

    parts% npart_i(0)= 0
    npart_i_tmp(0)   = 0
    parts% nbar_i    = zero
    parts% nuratio_i = zero

    loop_over_matter_objects: DO i_matter= 1, parts% n_matter, 1

      parts% adm_mass          = id% return_adm_mass()
      parts% masses(i_matter)  = id% return_mass(i_matter)
      center(i_matter,:)       = id% return_center(i_matter)
      !central_density(i_matter)= id% read_mass_density( center(i_matter,1), &
      !                                                  center(i_matter,2), &
      !                                                  center(i_matter,3) )
      barycenter(i_matter,:)       = id% return_barycenter(i_matter)
      parts% barycenter(i_matter,:)= barycenter(i_matter,:)
      sizes(i_matter, :)           = id% return_spatial_extent(i_matter)

      parts% all_eos(i_matter)% eos_name= id% return_eos_name(i_matter)
      CALL id% return_eos_parameters( i_matter, &
                                      parts% all_eos(i_matter)% eos_parameters )                          

      IF(parts% all_eos(i_matter)% eos_parameters(1) == eos$tabu$compose)THEN

        IF(ALLOCATED(id% tab_eos))THEN

          IF(ALLOCATED(id% tab_eos(i_matter)% table_eos))THEN

            parts% all_eos(i_matter)% table_eos= &
              id% tab_eos(i_matter)% table_eos

          ELSE

            PRINT *, "** ERROR! The EOS for matter object ", i_matter, &
                     " is supposed to be tabulated, since its EOS ", &
                     "identification number is ", &
                     parts% all_eos(i_matter)% eos_parameters(1), ", ", &
                     "but the table has not been read."
            PRINT *, " * Please read the EOS table within the constructor ", &
                     " of the appropriate TYPE that EXTENDS idbase."
            PRINT *, " * Stopping..."
            PRINT *
            STOP

          ENDIF

        ELSE

          PRINT *, "** ERROR! The EOS for matter object ", i_matter, &
                   " is supposed to be tabulated, since its EOS ", &
                   "identification number is ", &
                   parts% all_eos(i_matter)% eos_parameters(1), ", ", &
                   "but the table has not been allocated."
          PRINT *, " * Please allocate the EOS table within the constructor ", &
                   " of the appropriate TYPE that EXTENDS idbase."
          PRINT *, " * Stopping..."
          PRINT *
          STOP

        ENDIF

      ENDIF

    ENDDO loop_over_matter_objects

    ! Compute desired particle numbers based on mass ratios
    max_mass  = MAXVAL(parts% masses)
    total_mass= SUM(parts% masses)
    DO i_matter= 1, parts% n_matter, 1

      parts% mass_ratios(i_matter)   = parts% masses(i_matter)/max_mass
      parts% mass_fractions(i_matter)= parts% masses(i_matter)/total_mass
      npart_des_i(i_matter)          = &
                          NINT(parts% mass_fractions(i_matter)*DBLE(npart_des))
      ghost_dists(i_matter)          = &
                          ghost_dist/parts% mass_fractions(i_matter)
      tmp= 2*npart_des_i(i_matter)
      ALLOCATE( parts_all(i_matter)% pos_i  ( 3, tmp ) )
      ALLOCATE( parts_all(i_matter)% pvol_i ( tmp ) )
      ALLOCATE( parts_all(i_matter)% h_i    ( tmp ) )
      ALLOCATE( parts_all(i_matter)% nu_i   ( tmp ) )
      !parts_all(i_matter)% eos_i= parts% all_eos(i_matter)

    ENDDO
    ghost_dists(1)= ghost_dist
    !   IF( parts% redistribute_nu )THEN
    !     thres= 100.0D0*parts% nu_ratio
    !   ENDIF


    !
    !-- Check that the EOS specified in the ID file is the same as the one
    !-- specified in the parameter file SPHINCS_fm_input.dat
    !
    CALL check_eos()

    CALL id% initialize_id( flag$sph )

    DO i_matter= 1, parts% n_matter, 1

      central_density(i_matter)= id% read_mass_density( center(i_matter,1), &
                                                        center(i_matter,2), &
                                                        center(i_matter,3) )

    ENDDO

    !
    !-- Copy the surfaces of the matter objects, if they are known
    !
    IF(ALLOCATED(id% surfaces))THEN

      DO i_matter= 1, parts% n_matter, 1

        IF(ALLOCATED(id% surfaces(i_matter)% points))THEN

          parts% surfaces(i_matter)= id% surfaces(i_matter)

        ELSE

          parts% surfaces(i_matter)% is_known= .FALSE.

        ENDIF

      ENDDO

    ELSE

      DO i_matter= 1, parts% n_matter, 1

        parts% surfaces(i_matter)% is_known= .FALSE.

      ENDDO

    ENDIF


    parts% post_process_sph_id => id% finalize_sph_id_ptr


    ! TODO: Add check that the number of rows in placer is the same as the
    !       number of bns objects, and that all bns have a value for placer

    !
    !-- Choose particle placer
    !
    choose_particle_placer: SELECT CASE( dist )

    CASE( id_particles_from_formatted_file )


      CALL read_particles_from_formatted_file()                                        


    CASE( id_particles_on_lattice )


      CALL place_particles_on_lattices()


    CASE( id_particles_on_ellipsoidal_surfaces )


      CALL place_particles_on_ellipsoidal_surfaces()


    CASE DEFAULT

      PRINT *, "** There is no implemented particle placer " &
               // "corresponding to the number", dist
      PRINT *, " * Stopping..."
      PRINT *
      STOP

    END SELECT choose_particle_placer

    PRINT *, " * Particles placed. Number of particles=", parts% npart
    DO i_matter= 1, parts% n_matter, 1
      PRINT *, " * Number of particles on object ", i_matter, "=", &
               parts% npart_i(i_matter)
      PRINT *
    ENDDO
    PRINT *

    !---------------------------------------------------------------!
    !--  At this point, the particles are placed without the APM  --!
    !---------------------------------------------------------------!

    ! Check that there aren't particles with the same coordinates
    CALL parts% same_particle_timer% start_timer()
    check_particles_loop: DO i_matter= 1, parts% n_matter, 1

      CALL check_particle_positions( parts% npart_i(i_matter), &
                                     parts_all(i_matter)% pos_i )

    ENDDO check_particles_loop
    CALL parts% same_particle_timer% stop_timer()


    !
    !-- APM iteration
    !
    matter_objects_apm_loop: DO i_matter= 1, parts% n_matter, 1

      run_apm: IF( apm_iterate(i_matter) )THEN

        PRINT *
        PRINT *, "** Placing particles on matter object", i_matter, &
                 "using the APM..."
        PRINT *

        PRINT *, " * Particle number on object ", i_matter," before the APM=", &
                 parts% npart_i(i_matter)
        PRINT *

        IF( i_matter <= 9 ) WRITE( str_i, '(I1)' ) i_matter
        IF( i_matter >= 10 .AND. parts% n_matter <= 99 ) &
                                              WRITE( str_i, '(I2)' ) i_matter
        IF( i_matter >= 100 .AND. parts% n_matter <= 999 ) &
                                              WRITE( str_i, '(I3)' ) i_matter

        filename_apm_pos_id = TRIM(sph_path)//"apm_pos_id"//TRIM(str_i)//".dat"
        filename_apm_pos    = TRIM(sph_path)//"apm_pos"//TRIM(str_i)//".dat"
        filename_apm_results= TRIM(sph_path)//"apm_results"//TRIM(str_i)//".dat"

        !get_pressure_i => parts_all(i_matter)% compute_pressure

        CALL parts% apm_timers(i_matter)% start_timer()
        CALL perform_apm( &
                    ! PROCEDURES to get the density and pressure at a point
                    import_density, get_nstar_id, &
                    import_pressure_id, compute_pressure, &
                    ! Arguments pertaining to the matter object
                    parts% npart_i(i_matter), &
                    parts_all(i_matter)% pos_i, &
                    parts_all(i_matter)% pvol_i, &
                    parts_all(i_matter)% h_i, &
                    parts_all(i_matter)% nu_i, &
                    center(i_matter,:), barycenter(i_matter,:), &
                    parts% masses(i_matter), &
                    sizes(i_matter, :), &
                    parts% all_eos(i_matter), &
                    ! Steering parameters for the APM iteration
                    apm_max_it, max_inc, mass_it, parts% correct_nu, &
                    nuratio_thres, nuratio_des, use_pressure, &
                    ! Arguments pertaining to the ghost particles
                    adapt_ghosts, move_away_ghosts, &
                    nx_gh, ny_gh, nz_gh, ghost_dists(i_matter), &
                    ! Arguments pertaining to the atmosphere
                    use_atmosphere(i_matter), &
                    remove_atmosphere(i_matter), &
                    ! Arguments pertaining to input/output
                    print_step, filename_apm_pos_id, &
                    filename_apm_pos, filename_apm_results, &
                    ! Optional argument
                    validate_position, &
                    parts% surfaces(i_matter) )
        CALL parts% apm_timers(i_matter)% stop_timer()

        IF( debug ) PRINT *, "average nu= ", &
          SUM(parts_all(i_matter)% nu_i, DIM= 1)/SIZE(parts_all(i_matter)% nu_i)

        PRINT *, "** Particles placed on matter object", i_matter, &
                 "according to the APM."
        PRINT *

        ! If there are 2 matter objects...
        two_matter_objects: IF( i_matter == 1 .AND. parts% n_matter == 2 )THEN

          ! with practically the same mass, and the physical system
          ! is symmetric wrt the yz plane (in which case the user should set
          ! the reflect_particles_x to .TRUE. in the parameter file)
          equal_masses_apm: &
          IF( ABS(parts% masses(1) - parts% masses(2)) &
             /parts% masses(2) <= tol_equal_mass .AND. reflect_particles_x )THEN

            ! ...reflect particles

            DEALLOCATE(parts_all(2)% pos_i)
            DEALLOCATE(parts_all(2)% pvol_i)
            DEALLOCATE(parts_all(2)% h_i)
            DEALLOCATE(parts_all(2)% nu_i)

            CALL reflect_particles_yz_plane( parts_all(1)% pos_i,   &
                                             parts_all(1)% pvol_i,  &
                                             parts_all(1)% nu_i,    &
                                             parts_all(1)% h_i,     &
                                             parts% npart_i(1),     &
                                             parts_all(2)% pos_i,   &
                                             parts_all(2)% pvol_i,  &
                                             parts_all(2)% nu_i,    &
                                             parts_all(2)% h_i,     &
                                             parts% npart_i(2) )

            PRINT *, "** Particles placed on star 1 according to the APM,", &
                     " and reflected about the yz plane onto star 2."
            PRINT *

            EXIT

          ENDIF equal_masses_apm

        ENDIF two_matter_objects

      ENDIF run_apm

    ENDDO matter_objects_apm_loop

    PRINT *, " * Particle numbers after the APM=", &
             parts% npart_i(1:parts% n_matter)
    PRINT *

    parts% npart= SUM( parts% npart_i )

    ALLOCATE(parts% npart_fin(0:parts% n_matter))
    parts% npart_fin(0)= 0
    DO i_matter= 1, parts% n_matter, 1

      parts% npart_fin(i_matter)= SUM(parts% npart_i(1:i_matter))

    ENDDO

    !
    !-- Assign particle properties to the TYPE-bound variables
    !

    IF( ALLOCATED(parts% h) ) DEALLOCATE( parts% h )
    ALLOCATE( parts% h( parts% npart ), &
              STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...allocation error for array h in SUBROUTINE" &
                // "place_particles_. ", &
                "The error message is", err_msg
       STOP
    ENDIF
    IF( ALLOCATED(parts% pvol) ) DEALLOCATE( parts% pvol )
    ALLOCATE( parts% pvol( parts% npart ), &
              STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...allocation error for array pvol in SUBROUTINE" &
                // "place_particles_. ", &
                "The error message is", err_msg
       STOP
    ENDIF
    IF( ALLOCATED(parts% nu) ) DEALLOCATE( parts% nu )
      ALLOCATE( parts% nu( parts% npart ), &
                STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...allocation error for array nu in SUBROUTINE" &
                // "place_particles_. ", &
                "The error message is", err_msg
       STOP
    ENDIF
    IF( ALLOCATED(parts% pos) ) DEALLOCATE( parts% pos )
    ALLOCATE( parts% pos( 3, parts% npart ), &
              STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...allocation error for array pos in SUBROUTINE" &
                // "place_particles_. ", &
                "The error message is", err_msg
       STOP
    ENDIF

    parts% nbar_tot= zero
    DO i_matter= 1, parts% n_matter, 1

      ASSOCIATE( npart_in  => parts% npart_fin(i_matter-1) + 1, &
                 npart_fin => parts% npart_fin(i_matter) )

        parts% pos( :, npart_in : npart_fin )= parts_all(i_matter)% pos_i

        parts% nu( npart_in : npart_fin )= parts_all(i_matter)% nu_i

        parts% h( npart_in : npart_fin )= &
                  parts_all(i_matter)% h_i(1:parts% npart_i(i_matter))

        parts% pvol( npart_in : npart_fin )= &
                (parts_all(i_matter)% pvol_i(1:parts% npart_i(i_matter)))

        parts% nbar_i(i_matter)= SUM( parts% nu( npart_in : npart_fin ), DIM=1 )

        parts% nbar_tot= parts% nbar_tot + parts% nbar_i(i_matter)

      END ASSOCIATE

    ENDDO

    PRINT *, " * Final particle distribution determined. Number of particles=",&
             parts% npart
    DO i_matter= 1, parts% n_matter, 1
      PRINT *, " * Number of particles on object ", i_matter, "=", &
               parts% npart_i(i_matter)
      PRINT *
    ENDDO
    PRINT *

    CALL COM( parts% npart, parts% pos, parts% nu, &
              parts% barycenter_system(1), &
              parts% barycenter_system(2), &
              parts% barycenter_system(3), &
              parts% barycenter_system(4) )

    parts_out_namefile= "final_pos_nu.dat"

    PRINT *, "** Printing final particle positions and nu to file ", &
             TRIM(sph_path)//TRIM(parts_out_namefile), "..."

    INQUIRE( FILE= TRIM(sph_path)//TRIM(parts_out_namefile), EXIST= exist )

    IF( exist )THEN
      OPEN( UNIT= unit_pos_out, &
            FILE= TRIM(sph_path)//TRIM(parts_out_namefile), &
            STATUS= "REPLACE", FORM= "FORMATTED", &
            POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
            IOMSG= err_msg )
    ELSE
      OPEN( UNIT= unit_pos_out, &
            FILE= TRIM(sph_path)//TRIM(parts_out_namefile), &
            STATUS= "NEW", FORM= "FORMATTED", &
            ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening " // &
               TRIM(sph_path)//TRIM(parts_out_namefile), &
               ". The error message is", err_msg
      STOP
    ENDIF

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Particle distribution. The first line contains the number of matter ", &
    "objects (e.g. 2 for a binary neutron star, 1 for a single star or a ", &
    "neutron star-black hole system), and the number of particles on each", &
    "matter object. The next lines contain the positions and the ", &
    "baryon numbers of the particles."
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 1 in " // TRIM(parts_out_namefile), &
               ". The error message is", err_msg
      STOP
    ENDIF

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# column from the second row:      1        2       3       4"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 2 in " // TRIM(parts_out_namefile), &
               ". The error message is", err_msg
      STOP
    ENDIF

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "#      x [Msun_geo]       y [Msun_geo]       z [Msun_geo]      nu"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 3 in " // TRIM(parts_out_namefile), &
               ". The error message is", err_msg
      STOP
    ENDIF

    WRITE( UNIT = unit_pos_out, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      parts% n_matter, parts% npart_i(1:parts% n_matter)

    DO a= 1, parts% npart, 1
      WRITE( UNIT = unit_pos_out, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        parts% pos(1,a), parts% pos(2,a), parts% pos(3,a), parts% nu(a)
    ENDDO

    CLOSE( UNIT= unit_pos_out )

    PRINT *, " * final particle positions printed to file ", &
             TRIM(sph_path)//TRIM(parts_out_namefile)
    PRINT *

    !
    !-- Printouts about the baryon number ratios
    !
    DO i_matter= 1, parts% n_matter, 1

      ASSOCIATE( npart_in  => parts% npart_fin(i_matter-1) + 1, &
                 npart_fin => parts% npart_fin(i_matter) )

        parts% nuratio_i(i_matter)= &
          MAXVAL( parts% nu(npart_in:npart_fin), DIM= 1 ) &
         /MINVAL( parts% nu(npart_in:npart_fin), DIM= 1 )
        PRINT *, " * Maximum n. baryon per particle (nu) on object", i_matter, &
                            "=", MAXVAL( parts% nu(npart_in:npart_fin), DIM= 1 )
        PRINT *, " * Minimum n. baryon per particle (nu) on object", i_matter, &
                            "=", MINVAL( parts% nu(npart_in:npart_fin), DIM= 1 )
        PRINT *, " * Ratio between the two=", parts% nuratio_i(i_matter)
        PRINT *

        PRINT *, " * Number of baryons on object", i_matter, "=", &
                 parts% nbar_i(i_matter)
        PRINT *, " * Total mass of the baryons on object", i_matter, "=", &
                 parts% nbar_i(i_matter)*amu/Msun, "Msun =", &
                 parts% nbar_i(i_matter)*amu/Msun/parts% masses(i_matter), &
                 "of the baryon mass of object", i_matter, "."
        PRINT *

      END ASSOCIATE

    ENDDO

    parts% nuratio= MAXVAL( parts% nu, DIM= 1 )/MINVAL( parts% nu, DIM= 1 )
    PRINT *, " * Baryon number ratio across the stars=", parts% nuratio
    PRINT *
    PRINT *, " * Total mass of the baryons=", &
             parts% nbar_tot*amu/Msun, "Msun =", &
             parts% nbar_tot*amu/Msun/(SUM(parts% masses, DIM=1)), &
             "of the total baryon mass."
    PRINT *

    !
    !-- Adjusting the baryon number per particle uniformly so that
    !-- the baryon mass is correct.
    !
    IF( parts% correct_nu )THEN

      parts% nbar_tot= zero
      DO i_matter= 1, parts% n_matter, 1

        ASSOCIATE( npart_in  => parts% npart_fin(i_matter-1) + 1, &
                   npart_fin => parts% npart_fin(i_matter) )

        parts% nu( npart_in:npart_fin )= parts% nu( npart_in:npart_fin ) &
                    /(parts% nbar_i(i_matter)*amu/Msun/parts% masses(i_matter))

        parts% nbar_i(i_matter)= parts% nbar_i(i_matter) &
                      /(parts% nbar_i(i_matter)*amu/Msun/parts% masses(i_matter))

        parts% nbar_tot= parts% nbar_tot + parts% nbar_i(i_matter)

        PRINT *, " * Number of corrected baryons on object", i_matter, "=", &
                 parts% nbar_i(i_matter)
        PRINT *, " * Total mass of the corrected baryons object", i_matter, &
                 "=", parts% nbar_i(i_matter)*amu/Msun, "Msun =", &
                 parts% nbar_i(i_matter)*amu/Msun/parts% masses(i_matter), &
                 "of the baryon mass of object", i_matter, "."

        END ASSOCIATE

      ENDDO

      PRINT *, " * Total number of corrected baryons=", parts% nbar_tot
      PRINT *, " * Total mass of the corrected baryons=", &
               parts% nbar_tot*amu/Msun, "Msun =", &
               parts% nbar_tot*amu/Msun/(SUM(parts% masses, DIM=1)), &
               "of the total baryon mass."
      PRINT *

    ENDIF

    !--------------------------------------------------------------------!
    !--  At this point, the final particle distribution is determined, --!
    !--  and nu and the first guess for h are assigned.                --!
    !--  Now the ID can be read on the particle positions.             --!
    !--------------------------------------------------------------------!

    ! Allocate needed memory
    CALL parts% allocate_particles_memory()

    ! flag that particles are 'alive'
    ALLOCATE( alive( parts% npart ) )
    alive( 1:parts% npart )= 1

    IF( debug ) PRINT *, "33"

    !
    !-- Read the needed ID on the particles, and time the process
    !
    PRINT *, "** Assigning the ID to the particles..."
    PRINT *

    CALL parts% importer_timer% start_timer()
    CALL id% read_id_particles( parts% npart, &
                                parts% pos(1,:), &
                                parts% pos(2,:), &
                                parts% pos(3,:), &
                                parts% lapse, &
                                parts% shift_x, &
                                parts% shift_y, &
                                parts% shift_z, &
                                parts% g_xx, &
                                parts% g_xy, &
                                parts% g_xz, &
                                parts% g_yy, &
                                parts% g_yz, &
                                parts% g_zz, &
                                parts% baryon_density, &
                                parts% energy_density, &
                                parts% specific_energy, &
                                parts% pressure, &
                                parts% v_euler_x, &
                                parts% v_euler_y, &
                                parts% v_euler_z )
    CALL parts% importer_timer% stop_timer()

    IF( debug ) PRINT *, "34"

    !-----------------------------------------------------------------------!
    ! If an atmosphere was used during the APM iteration, and kept, assign  !
    ! the minimum specific internal energy and the minimum velocity, to it. !
    ! N.B. The velocity has a hard-wired direction to reproduce counter-    !
    !      clockwise rotation.                                              !
    !-----------------------------------------------------------------------!

    matter_objects_atmo_loop: DO i_matter= 1, parts% n_matter, 1

    ASSOCIATE( npart_in  => parts% npart_fin(i_matter-1) + 1, &
               npart_fin => parts% npart_fin(i_matter) )

      IF( use_atmosphere(i_matter) .AND. .NOT.remove_atmosphere(i_matter) )THEN

        min_eps= MINVAL( parts% specific_energy(npart_in:npart_fin), &
                         DIM= 1, &
                MASK= parts% baryon_density(npart_in:npart_fin) > zero )
        min_vel= MINVAL( SQRT( &
                       (parts% v_euler_x(npart_in:npart_fin))**2 &
                     + (parts% v_euler_y(npart_in:npart_fin))**2 &
                     + (parts% v_euler_z(npart_in:npart_fin))**2 ), &
                         DIM= 1, &
                MASK= parts% baryon_density(npart_in:npart_fin) > zero )

        particle_loop2: DO a= npart_in, npart_fin, 1

          IF( parts% baryon_density(a) <= zero )THEN

            CALL spherical_from_cartesian( &
                  parts% pos(1,a), parts% pos(2,a), parts% pos(3,a), &
                  center(i_matter,1), center(i_matter,2), center(i_matter,3), &
                  r_a, theta_a, phi_a )

            parts% specific_energy(a)= min_eps

            parts% v_euler_x(a)      = &
      ( min_vel*SIN(theta_a - pi*half)*COS(phi_a) + parts% shift_x(a) ) &
              /parts% lapse(a)
            parts% v_euler_y(a)      = &
      ( min_vel*SIN(theta_a - pi*half)*SIN(phi_a) + parts% shift_y(a) ) &
              /parts% lapse(a)
            parts% v_euler_z(a)      = &
              ( min_vel*COS(theta_a - pi*half) + parts% shift_z(a) ) &
              /parts% lapse(a)

          ENDIF

        ENDDO particle_loop2

      ENDIF

    END ASSOCIATE

    ENDDO matter_objects_atmo_loop


    !
    !-- Ensure that the ID does not contain NaNs or infinities
    !
    PRINT *, "** Ensuring that the ID does not have any NaNs or infinities..."

    CALL scan_1d_array_for_nans( parts% npart, parts% lapse, "lapse" )

    CALL scan_1d_array_for_nans( parts% npart, parts% shift_x, "shift_x" )
    CALL scan_1d_array_for_nans( parts% npart, parts% shift_y, "shift_y" )
    CALL scan_1d_array_for_nans( parts% npart, parts% shift_z, "shift_z" )

    CALL scan_1d_array_for_nans( parts% npart, parts% g_xx, "g_xx" )
    CALL scan_1d_array_for_nans( parts% npart, parts% g_xy, "g_xy" )
    CALL scan_1d_array_for_nans( parts% npart, parts% g_xz, "g_xz" )
    CALL scan_1d_array_for_nans( parts% npart, parts% g_yy, "g_yy" )
    CALL scan_1d_array_for_nans( parts% npart, parts% g_yz, "g_yz" )
    CALL scan_1d_array_for_nans( parts% npart, parts% g_zz, "g_zz" )

    CALL scan_1d_array_for_nans( parts% npart, &
                                 parts% baryon_density, "baryon_density" )
    CALL scan_1d_array_for_nans( parts% npart, &
                                 parts% energy_density, "energy_density" )
    CALL scan_1d_array_for_nans( parts% npart, &
                                 parts% specific_energy, "specific_energy" )
    CALL scan_1d_array_for_nans( parts% npart, &
                                 parts% pressure, "pressure" )

    CALL scan_1d_array_for_nans( parts% npart, &
                                 parts% v_euler_x, "v_euler_x" )
    CALL scan_1d_array_for_nans( parts% npart, &
                                 parts% v_euler_y, "v_euler_y" )
    CALL scan_1d_array_for_nans( parts% npart, &
                                 parts% v_euler_z, "v_euler_z" )

    PRINT *, " * the ID does not have NaNs or infinities."
    PRINT *


    !
    !-- Compute typical length-scale approximating g_00 with the Newtonian
    !-- potential
    !
    DO i_matter= 1, parts% n_matter, 1

      !ASSOCIATE( npart_in   => parts% npart_i(i_matter-1) + 1, &
      !           npart_fin  => parts% npart_i(i_matter-1) +    &
      !                         parts% npart_i(i_matter) )
      ASSOCIATE( npart_in  => parts% npart_fin(i_matter-1) + 1, &
                 npart_fin => parts% npart_fin(i_matter) )

        min_g00_abs= HUGE(one)
        DO itr= npart_in, npart_fin, 1

          CALL spatial_vector_norm_sym3x3( &
                  [parts% g_xx(itr), parts% g_xy(itr), parts% g_xz(itr), &
                   parts% g_yy(itr), parts% g_yz(itr), parts% g_zz(itr)], &
              [parts% shift_x(itr), parts% shift_y(itr), parts% shift_z(itr)], &
                   shift_norm )

          IF( min_g00_abs > parts% lapse(itr)**2 - shift_norm )THEN
            min_g00_abs= parts% lapse(itr)**2 - shift_norm
          ENDIF

        ENDDO
        min_lapse= MINVAL( parts% lapse, DIM= 1 )
        IF( one == min_lapse )THEN
          lapse_lengthscales= HUGE(one)/(ten*ten*ten)
        ELSE
          lapse_lengthscales= two*parts% masses(i_matter)/( one - min_lapse )
        ENDIF
        IF( one == min_g00_abs )THEN
          g00_lengthscales= HUGE(one)/(ten*ten*ten)
        ELSE
          g00_lengthscales= two*parts% masses(i_matter)/( one - min_g00_abs )
        ENDIF

      END ASSOCIATE

    ENDDO
    PRINT *, "** Approximating the g_00 component of the metric as a ", &
             "Newtonian potential (!) and neglecting the shift (!), ", &
             "the minimum lengthscales given by ", &
             "the lapse on each matter object are: "
    DO i_matter= 1, parts% n_matter, 1
      PRINT *, " * Matter object ", i_matter, "=", &
               lapse_lengthscales(i_matter), "Msun_geo=", &
               lapse_lengthscales(i_matter)*Msun_geo, "km"
    ENDDO
    PRINT *
    PRINT *, "** Approximating the g_00 component of the metric as a ", &
             "Newtonian potential (!), ", &
             "the minimum lengthscales given by ", &
             "g_00 on each matter object are: "
    DO i_matter= 1, parts% n_matter, 1
      PRINT *, " * Matter object ", i_matter, "=", &
               g00_lengthscales(i_matter), "Msun_geo=", &
               g00_lengthscales(i_matter)*Msun_geo, "km"
    ENDDO
    PRINT *


    ! Increase the counter that identifies the particle distribution
    counter= counter + 1



    CONTAINS



    FUNCTION import_density( x, y, z ) RESULT( density )
    !! Wrapper function to read the baryon mass density from the |id|

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: x
      DOUBLE PRECISION, INTENT(IN):: y
      DOUBLE PRECISION, INTENT(IN):: z
      DOUBLE PRECISION:: density

      density= id% read_mass_density( x, y, z )

    END FUNCTION import_density


    FUNCTION import_pressure_id( x, y, z ) RESULT( pressure )
    !! Wrapper function to read the pressure from the |id|

      USE constants,  ONLY: Msun, amu

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: x
      DOUBLE PRECISION, INTENT(IN):: y
      DOUBLE PRECISION, INTENT(IN):: z
      DOUBLE PRECISION:: pressure

      pressure= id% read_pressure( x, y, z )*Msun/amu

    END FUNCTION import_pressure_id


    SUBROUTINE import_id( x, y, z, &
                          sqdetg, &
                          baryon_density, &
                          gamma_euler )
    !# Wrapper function to read the ID necessary to compute the relativistic
    !  baryonic mass

      USE utility,  ONLY: determinant_sym3x3

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: x
      DOUBLE PRECISION, INTENT(IN) :: y
      DOUBLE PRECISION, INTENT(IN) :: z
      DOUBLE PRECISION, INTENT(OUT):: sqdetg
      DOUBLE PRECISION, INTENT(OUT):: baryon_density
      DOUBLE PRECISION, INTENT(OUT):: gamma_euler

      DOUBLE PRECISION, DIMENSION(6) :: g

      CALL id% read_id_mass_b( x, y, z, &
                               g, &
                               baryon_density, &
                               gamma_euler )

      CALL determinant_sym3x3(g,sqdetg)
      sqdetg= SQRT(sqdetg)

    END SUBROUTINE import_id


    SUBROUTINE integrate_mass_density &
      ( center, radius, central_density, dr, dth, dphi, mass, mass_profile, &
        mass_profile_idx, radii, surf )
    !# Wrapper function to integrate the relativistic baryonic mass density

      IMPLICIT NONE

      !> Center of the star
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN)    :: center
      !> Central density of the star
      DOUBLE PRECISION, INTENT(IN)    :: central_density
      !> Radius of the star
      DOUBLE PRECISION, INTENT(IN)    :: radius
      !> Integration steps
      DOUBLE PRECISION, INTENT(IN)    :: dr, dth, dphi
      !> Integrated mass of the star
      DOUBLE PRECISION, INTENT(INOUT):: mass
      !> Array storing the radial mass profile of the star
      !DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT):: &
      !                                 mass_profile
      DOUBLE PRECISION, DIMENSION(3,0:NINT(radius/dr)), INTENT(OUT):: &
                                           mass_profile
      !& Array to store the indices for array mass_profile, sorted so that
      !  mass_profile[mass_profile_idx] is in increasing order
      !INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT):: mass_profile_idx
      INTEGER, DIMENSION(0:NINT(radius/dr)), INTENT(OUT):: mass_profile_idx
      DOUBLE PRECISION, DIMENSION(2), INTENT(IN), OPTIONAL:: radii
      !> Surface of the matter object
      TYPE(surface),                  INTENT(IN), OPTIONAL:: surf

      CALL id% integrate_baryon_mass_density &
        ( center, radius, central_density, dr, dth, dphi, &
          mass, mass_profile, mass_profile_idx, radii, surf )

    END SUBROUTINE integrate_mass_density


    FUNCTION validate_position( x, y, z ) RESULT( answer )
    !! Wrapper function to validate a position

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: x
      DOUBLE PRECISION, INTENT(IN):: y
      DOUBLE PRECISION, INTENT(IN):: z
      LOGICAL:: answer

      answer= id% test_position( x, y, z )

    END FUNCTION validate_position


    SUBROUTINE correct_center_of_mass_of_system( npart, pos, nu, &
                                                 com_system )
    !! Set the \(COM\) of the system to `com_system`

      IMPLICIT NONE

      INTEGER, INTENT(IN):: npart
      DOUBLE PRECISION, INTENT(IN)   :: com_system(3)
      DOUBLE PRECISION, INTENT(INOUT):: nu(npart)
      DOUBLE PRECISION, INTENT(INOUT):: pos(3,npart)

      DOUBLE PRECISION:: nstar_id(npart)
      DOUBLE PRECISION:: nstar_eul_id(npart)
      DOUBLE PRECISION:: nu_eul(npart)

      INTEGER:: a, itr2

      PRINT *, "1"

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( npart, pos ) &
      !$OMP             PRIVATE( a, itr2 )
      find_nan_in_pos: DO a= 1, npart, 1

        DO itr2= 1, 3, 1
          IF( .NOT.is_finite_number(pos(itr2,a)) )THEN
            PRINT *, "** ERROR! pos(", itr2, a, ")= ", pos(itr2,a), &
                     " is not a finite number!"
            PRINT *, " * Stopping.."
            PRINT *
            STOP
          ENDIF
        ENDDO

      ENDDO find_nan_in_pos
      !$OMP END PARALLEL DO

      CALL correct_center_of_mass( npart, pos, nu, import_density, &
                                   validate_position, com_system, &
                                   verbose= .TRUE. )

      PRINT *, "4"

    END SUBROUTINE correct_center_of_mass_of_system


    SUBROUTINE get_nstar_id(npart, x, y, z, nstar_sph, nstar_id, nlrf_sph, sqg)
    !! Wrapper function to compute the relativistic baryon mass density

      IMPLICIT NONE

      INTEGER,          INTENT(IN):: npart
      DOUBLE PRECISION, INTENT(IN):: x(npart)
      DOUBLE PRECISION, INTENT(IN):: y(npart)
      DOUBLE PRECISION, INTENT(IN):: z(npart)
      DOUBLE PRECISION, INTENT(IN):: nstar_sph(npart)
      DOUBLE PRECISION, INTENT(OUT):: nstar_id(npart)
      DOUBLE PRECISION, INTENT(OUT):: nlrf_sph(npart)
      DOUBLE PRECISION, INTENT(OUT):: sqg(npart)

      DOUBLE PRECISION, DIMENSION(npart):: lapse, &
                                           shift_x, shift_y, shift_z, &
                                           g_xx, g_xy, g_xz, &
                                           g_yy, g_yz, g_zz, &
                                           baryon_density, &
                                           energy_density, &
                                           !specific_energy, &
                                           pressure, &
                                           v_euler_x, v_euler_y, v_euler_z

      ! compute_sph_hydro in SUBMODULE sph_particles@sph_variables requires the
      ! knowledge of parts% specific_energy, to compute the SPH pressure for a
      ! hot system in the APM, when using the real pressure to compute the
      ! artifical pressure.
      ! That is why we allocate (if necessary) and assign values to
      ! parts% specific_energy
      ! TODO: Another strategy would be adding the specific energy as an
      !       optional argument to compute_sph_hydro.
      ! TODO: Note that, since the pressure from the ID is not known for the
      !       ejecta, the SPH pressure computed in the APM cannot be compared
      !       with the pressure from the ID. One need to compute the pressure
      !       also using the same internalenergy, but the density from the ID
      !       (not the SPH density)
      IF(ALLOCATED(parts% specific_energy))THEN

        IF(SIZE(parts% specific_energy) /= npart)THEN

          DEALLOCATE(parts% specific_energy)
          ALLOCATE(parts% specific_energy(npart))

        ENDIF

      ELSE

        ALLOCATE(parts% specific_energy(npart))

      ENDIF

      CALL id% read_id_particles( npart, x, y, z, &
                                  lapse, shift_x, shift_y, shift_z, &
                                  g_xx, g_xy, g_xz, &
                                  g_yy, g_yz, g_zz, &
                                  baryon_density, &
                                  energy_density, &
                                  parts% specific_energy, &
                                  pressure, &
                                  v_euler_x, v_euler_y, v_euler_z )

 !     !$OMP PARALLEL DO DEFAULT( NONE ) &
 !     !$OMP             SHARED( npart, nlrf_id, baryon_density ) &
 !     !$OMP             PRIVATE( a )
 !     DO a= 1, npart, 1
 !
 !       nlrf_id(a)= baryon_density(a)
 !
 !     ENDDO
 !     !$OMP END PARALLEL DO

      CALL compute_nstar_id( npart, lapse, shift_x, shift_y, &
                             shift_z, v_euler_x, v_euler_y, v_euler_z, &
                             g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, &
                             baryon_density, nstar_sph, nstar_id, nlrf_sph, &
                             sqg )

    END SUBROUTINE get_nstar_id


    SUBROUTINE compute_nstar_id( npart, lapse, shift_x, shift_y, &
                                 shift_z, v_euler_x, v_euler_y, v_euler_z, &
                                 g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, &
                                 baryon_density, nstar_sph, nstar_id, nlrf_sph,&
                                 sqg )

      !**************************************************************
      !
      !# Compute nstar_id, the relativistic baryon mass density,
      !  given the required |id| as input
      !
      !  FT 31.08.2021
      !
      !**************************************************************

      USE tensor,   ONLY: jx, jy, jz, n_sym4x4
      USE utility,  ONLY: compute_g4, determinant_sym4x4, &
                          spacetime_vector_norm_sym4x4, zero, one, two

      IMPLICIT NONE

      INTEGER, INTENT(IN):: npart
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: lapse
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: shift_x
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: shift_y
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: shift_z
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: v_euler_x
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: v_euler_y
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: v_euler_z
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: g_xx
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: g_xy
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: g_xz
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: g_yy
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: g_yz
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: g_zz
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: baryon_density
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: nstar_sph
      DOUBLE PRECISION, DIMENSION(npart), INTENT(OUT):: nstar_id
      DOUBLE PRECISION, DIMENSION(npart), INTENT(OUT):: nlrf_sph
      DOUBLE PRECISION, DIMENSION(npart), INTENT(OUT):: sqg

      INTEGER:: a, i!mus, nus
      DOUBLE PRECISION:: det, sq_g, Theta_a
      DOUBLE PRECISION, DIMENSION(0:3,npart):: vel
      !DOUBLE PRECISION:: g4(0:3,0:3)
      DOUBLE PRECISION:: g4(n_sym4x4)

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( npart, lapse, shift_x, shift_y, shift_z, &
      !$OMP                     v_euler_x, v_euler_y, v_euler_z, &
      !$OMP                     g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, &
      !$OMP                     baryon_density, vel, nstar_id, nstar_sph, &
      !$OMP                     nlrf_sph, sqg ) &
      !$OMP             PRIVATE( a, det, sq_g, Theta_a, g4 )
      DO a= 1, npart, 1

        ! Coordinate velocity of the fluid [c]
        vel(0,a) = one
        vel(jx,a)= lapse(a)*v_euler_x(a) - shift_x(a)
        vel(jy,a)= lapse(a)*v_euler_y(a) - shift_y(a)
        vel(jz,a)= lapse(a)*v_euler_z(a) - shift_z(a)

        CALL compute_g4( lapse(a), [shift_x(a),shift_y(a),shift_z(a)], &
                         [g_xx(a),g_xy(a),g_xz(a),g_yy(a),g_yz(a),g_zz(a)], g4 )

        CALL determinant_sym4x4( g4, det )

        IF( ABS(det) < 1D-10 )THEN
          PRINT *, "ERROR! The determinant of the spacetime metric is " &
                   // "effectively 0 at particle ", a
          PRINT *
          STOP
        ELSEIF( det > 0 )THEN
          PRINT *, "ERROR! The determinant of the spacetime metric is " &
                   // "positive at particle ", a
          PRINT *
          STOP
        ELSEIF( .NOT.is_finite_number(det) )THEN
          PRINT *, "ERROR! The determinant is ", det, "at particle ", a
          PRINT *
          STOP
        ENDIF
        sq_g= SQRT(-det)

        !
        !-- Generalized Lorentz factor
        !
        Theta_a= zero
        CALL spacetime_vector_norm_sym4x4( g4, vel(:,a), Theta_a )
        IF( .NOT.is_finite_number(Theta_a) )THEN
          PRINT *, "** ERROR! The spacetime norm of vel is ", Theta_a, &
                   "at particle ", a, &
                   "in SUBROUTINE compute_nstar_id"
          PRINT *, " * Stopping..."
          PRINT *
          STOP
        ENDIF

        Theta_a= one/SQRT(-Theta_a)
        IF( .NOT.is_finite_number(Theta_a) )THEN
          PRINT *, "** ERROR! The generalized Lorentz factor is ", Theta_a, &
                   "at particle ", a, &
                   "in SUBROUTINE compute_nstar_id"
          PRINT *, " * Stopping..."
          PRINT *
          STOP
        ENDIF
        IF( Theta_a < one )THEN
          PRINT *, "** ERROR! The generalized Lorentz factor is ", Theta_a, &
                   "< 1 at particle ", a, &
                   "in SUBROUTINE compute_nstar_id"
          PRINT *, " * Stopping..."
          PRINT *
          STOP
        ENDIF

        nstar_id(a)= sq_g*Theta_a*baryon_density(a)
        nlrf_sph(a)= nstar_sph(a)/(sq_g*Theta_a)
        sqg(a)     = sq_g

      ENDDO
      !$OMP END PARALLEL DO

    END SUBROUTINE compute_nstar_id


    SUBROUTINE read_particles_options

      !**************************************************************
      !
      !# Read the parameters in the file sphincs_id_particles.dat
      !
      !  FT 2022
      !
      !**************************************************************

      IMPLICIT NONE

      INTEGER, PARAMETER:: unit_particles= 6534

      NAMELIST /sphincs_id_particles/ &
                parts_pos_path, parts_pos, columns, header_lines, n_cols, &
                read_nu, column_nu, stretch, &
                use_thres, thres, nu_ratio_des, redistribute_nu, correct_nu, &
                compose_eos, compose_path, compose_filename, &
                npart_des, last_r, upper_bound, lower_bound, &
                upper_factor, lower_factor, max_steps, &
                randomize_phi, randomize_theta, randomize_r, &
                apm_iterate, apm_max_it, max_inc, mass_it, &
                nuratio_thres, reflect_particles_x, nx_gh, ny_gh, nz_gh, &
                use_atmosphere, remove_atmosphere, nuratio_des, print_step, &
                ghost_dist, adapt_ghosts, move_away_ghosts, use_pressure

      parts% sphincs_id_particles= 'sphincs_id_particles.dat'

      INQUIRE( FILE= parts% sphincs_id_particles, EXIST= file_exists )
      IF( file_exists )THEN
        OPEN( unit_particles, FILE= parts% sphincs_id_particles, STATUS= 'OLD' )
      ELSE
        PRINT *
        PRINT *, "** ERROR: ", parts% sphincs_id_particles, " file not found!"
        PRINT *
        STOP
      ENDIF

      apm_iterate      = .FALSE.
      use_atmosphere   = .FALSE.
      remove_atmosphere= .FALSE.

      compose_path    = "compose_path is a deprecated variable"
      compose_filename= "compose_filename is a deprecated variable"

      READ(unit_particles, NML= sphincs_id_particles)
      CLOSE(unit_particles)

      parts% use_thres          = use_thres
      parts% correct_nu         = correct_nu
      parts% compose_eos        = compose_eos
      parts% compose_path       = compose_path
      parts% compose_filename   = compose_filename
      parts% redistribute_nu    = redistribute_nu
      parts% nu_ratio_des       = nu_ratio_des
      parts% reflect_particles_x= reflect_particles_x
      parts% randomize_phi      = randomize_phi
      parts% randomize_theta    = randomize_theta
      parts% randomize_r        = randomize_r
      ! APM parameters
      ALLOCATE( parts% apm_iterate(parts% n_matter) )
      parts% apm_iterate   = apm_iterate(1:parts% n_matter)
      parts% use_atmosphere= use_atmosphere
      parts% read_nu       = read_nu

      parts_pos_namefile= TRIM(parts_pos_path)//TRIM(parts_pos)

      !
      !-- Check that the parameters are acceptable
      !
      IF( upper_bound <= lower_bound )THEN
        PRINT *
        PRINT *, "** ERROR in ", parts% sphincs_id_particles, &
                 "upper_bound should be greater than lower_bound!"
        PRINT *
        STOP
      ENDIF
      IF( upper_factor < 1.0D0 )THEN
        PRINT *
        PRINT *, "** ERROR in ", parts% sphincs_id_particles, &
                 "upper_factor should be greater than or equal to 1!"
        PRINT *
        STOP
      ENDIF
      IF( lower_factor > 1 )THEN
        PRINT *
        PRINT *, "** ERROR in ", parts% sphincs_id_particles, &
                 "lower_factor should be smaller than or equal to 1!"
        PRINT *
        STOP
      ENDIF
      IF( max_steps < 10 )THEN
        PRINT *
        PRINT *, "** ERROR in ", parts% sphincs_id_particles, &
                 "max_steps should be an integer greater than or equal to 10!"
        PRINT *
        STOP
      ENDIF
      IF( last_r < 0.95D0 .OR. last_r > 1.0D0 )THEN
        PRINT *
        PRINT *, "** ERROR in ", parts% sphincs_id_particles, &
                 "last_r should be greater than or equal to 0.95, ", &
                 "and lower than or equal to 1!"
        PRINT *
        STOP
      ENDIF
      IF( apm_max_it < 0 .OR. max_inc < 0 .OR. nuratio_thres < 0 &
          .OR. nuratio_des < 0 .OR. nx_gh < 0 .OR. ny_gh < 0 .OR. nz_gh < 0 )THEN
        PRINT *
        PRINT *, "** ERROR in ", parts% sphincs_id_particles, &
                 "the numeric parameters for the APM method should be positive!"
        PRINT *
        STOP
      ENDIF
      IF( nuratio_des >= nuratio_thres )THEN
        PRINT *
        PRINT *, "** ERROR in ", parts% sphincs_id_particles, &
                 "nuratio_des has to be stricly lower than nuratio_thres!"
        PRINT *
        STOP
      ENDIF
      IF( print_step < 0 )THEN
        PRINT *
        PRINT *, "** ERROR in ", parts% sphincs_id_particles, &
                 "print_step has to be a positive integer or zero!"
        PRINT *
        STOP
      ENDIF
      IF( ghost_dist < zero )THEN
        PRINT *
        PRINT *, "** ERROR in ", parts% sphincs_id_particles, &
                 "ghost_dist has to be a positive double precision or zero!"
        PRINT *
        STOP
      ENDIF

      ! setup unit system
      CALL set_units('NSM')
      CALL read_options

      ! tabulate kernel, get ndes
      CALL ktable( ikernel, ndes )

    END SUBROUTINE read_particles_options


    SUBROUTINE check_eos

      !**************************************************************
      !
      !# Check that the supplied |eos| parameters are consistent
      !  with the |eos| used to compute the |id|
      !
      !  FT xx.09.2022
      !
      !**************************************************************

      USE utility,  ONLY: eos$poly, eos$pwpoly, eos$tabu$compose
      USE options,  ONLY: eos_str, eos_type

      IMPLICIT NONE

      IF( (eos_type /= 'Poly') .AND. (eos_type /= 'pwp') )THEN
        PRINT *, "** ERROR! Unkown EOS specified in parameter file ", &
                 "SPHINCS_fm_input.dat."
        PRINT *, " * The currently supported EOS types are 'Poly' for a ", &
                 "polytropic EOS, and 'pwp' for a piecewise polytropic EOS."
        PRINT *
        PRINT *, " * EOS from the parameter file SPHINCS_fm_input.dat: ", &
                 eos_type
        PRINT *, " * Stopping..."
        PRINT *
        STOP
      ENDIF

      DO i_matter= 1, parts% n_matter, 1

        IF( parts% all_eos(i_matter)% eos_parameters(1) == eos$poly )THEN

          IF( compose_eos )THEN
            PRINT *, "** ERROR! On matter object ", i_matter, &
                     ", the EOS taken from the ID is a single polytrope, ", &
                     "so the parameter compose_eos should be set to .FALSE. ", &
                     "in sphincs_id_particles.dat."
            PRINT *
            PRINT *, " * EOS from the ID: ", &
                     parts% all_eos(i_matter)% eos_name
            PRINT *, " * Stopping..."
            PRINT *
            STOP
          ENDIF

          IF( eos_type == 'pwp' )THEN
            PRINT *, "** ERROR! On matter object ", i_matter, &
                     ", the EOS taken from the ID is not the same as the ",&
                     "one specified in parameter file SPHINCS_fm_input.dat."
            PRINT *
            PRINT *, " * EOS from the ID: ", &
                     parts% all_eos(i_matter)% eos_name
            PRINT *, " * EOS from the parameter file SPHINCS_fm_input.dat: ", &
                     eos_type
            PRINT *, " * Stopping..."
            PRINT *
            STOP
          ENDIF

        ENDIF

        IF( parts% all_eos(i_matter)% eos_parameters(1) == eos$pwpoly )THEN

          IF( compose_eos )THEN
            PRINT *, "** ERROR! On matter object ", i_matter, &
                     ", the EOS taken from the ID is a piecewise polytrope, ", &
                     "so the parameter compose_eos should be set to .FALSE. ", &
                     "in sphincs_id_particles.dat."
            PRINT *
            PRINT *, " * EOS from the ID: ", &
                     parts% all_eos(i_matter)% eos_name
            PRINT *, " * Stopping..."
            PRINT *
            STOP
          ENDIF

          IF( eos_type == 'Poly' )THEN
            PRINT *, "** ERROR! On matter object ", i_matter, &
                     ", the EOS taken from the ID is not the same as the ",&
                     "one specified in parameter file SPHINCS_fm_input.dat."
            PRINT *
            PRINT *, " * EOS from the ID: ", &
                     parts% all_eos(i_matter)% eos_name
            PRINT *, " * EOS from the parameter file SPHINCS_fm_input.dat: ", &
                     eos_type
            PRINT *, " * Stopping..."
            PRINT *
            STOP
          ENDIF

          IF( (parts% all_eos(i_matter)% eos_name .LT. eos_str)&
              .OR. &
              (parts% all_eos(i_matter)% eos_name .GT. eos_str)&
          )THEN

            PRINT *, "** ERROR! On matter object ", i_matter, &
                     ", the EOS taken from the ID is not the same as the ",&
                     "one specified in parameter file SPHINCS_fm_input.dat."
            PRINT *
            PRINT *, " * EOS from the ID: ", parts% all_eos(i_matter)% eos_name
            PRINT *, " * EOS from the parameter file SPHINCS_fm_input.dat: ", &
                     eos_str
            PRINT *, " * Stopping..."
            PRINT *
            STOP

          ENDIF

        ENDIF

      ENDDO

    END SUBROUTINE check_eos


    SUBROUTINE read_particles_from_formatted_file

      !**************************************************************
      !
      !# Read particles from formatted file, and
      !  reflect  particles with respect to the yz plane in the case
      !  of equal-mass binaries
      !
      !  FT 21.10.2022
      !
      !**************************************************************

      IMPLICIT NONE

      PRINT *, " * Reading particle positions from formatted file " &
               // TRIM(parts_pos_namefile)
      PRINT *

      INQUIRE( FILE= TRIM(parts_pos_namefile), EXIST= exist )

      IF( exist )THEN
        OPEN( UNIT= unit_pos, FILE= TRIM(parts_pos_namefile), &
              FORM= "FORMATTED", ACTION= "READ", IOSTAT= ios, &
              IOMSG= err_msg )
        IF( ios > 0 )THEN
          PRINT *, "...error when opening " // TRIM(parts_pos_namefile), &
                  ". The error message is", err_msg
          STOP
        ENDIF
      ELSE
        PRINT *, "** ERROR! Unable to find file " // TRIM(parts_pos_namefile)
        STOP
      ENDIF

      ! Get total number of lines in the file
      nlines = 0
      DO
        READ( unit_pos, * , IOSTAT= ios )
        IF ( ios /= 0 ) EXIT
        nlines = nlines + 1
      ENDDO

      IF( debug ) PRINT *, "nlines=", nlines

      CLOSE( UNIT= unit_pos )

      ! Set the total number of particles to the number of lines in the file,
      ! minus the number of header lines, minus the line containing the number
      ! of particles on each matter object
      npart_tmp= nlines - header_lines - 1

      IF( debug ) PRINT *, "npart_tmp=", npart_tmp

      ! Read all particle positions, and nu, if present
      OPEN( UNIT= unit_pos, FILE= TRIM(parts_pos_namefile), &
            FORM= "FORMATTED", ACTION= "READ" )

      ! Skip header
      DO itr= 1, header_lines, 1
        READ( unit_pos, * )
      ENDDO

      ! Read the number of matter objects and the particle numbers on each
      ! matter object
      READ( UNIT= unit_pos, FMT= *, IOSTAT = ios, IOMSG= err_msg ) &
              n_matter_tmp, npart_i_tmp(1:parts% n_matter)

      IF( ios > 0 )THEN
        PRINT *, "...error when reading " // TRIM(parts_pos_namefile), &
                " at particle ", itr,". The status variable is ", ios, &
                ". The error message is", err_msg
        STOP
      ENDIF

      ! Check that the numbers of matter objects is consistent
      IF( n_matter_tmp /= parts% n_matter )THEN
        PRINT *, "** ERROR! The numbers of matter objects", &
                 " in file ", TRIM(parts_pos_namefile), ", equal to ", &
                 n_matter_tmp, ", is not consistent", &
                 " with the one corresponding to ID file, equal to", &
                 parts% n_matter
        PRINT *, "   Stopping..."
        PRINT *
        STOP
      ENDIF

      ! Check that the numbers of particles are consistent
      IF( npart_tmp /= SUM(npart_i_tmp) )THEN
        PRINT *, "** ERROR! The numbers of particles on each matter object", &
                 " do not add up to the total number of particles, in file ", &
                 TRIM(parts_pos_namefile)
        PRINT *, " * npart_tmp= ", npart_tmp
        PRINT *, " * npart_i_tmp=", npart_i_tmp
        PRINT *, " * SUM(npart_i_tmp)=", SUM(npart_i_tmp)
        PRINT *, " * Stopping..."
        PRINT *
        STOP
      ENDIF

      parts% npart_i(1:parts% n_matter)= npart_i_tmp(1:parts% n_matter)
      parts% npart = SUM(parts% npart_i)

      CALL parts% placer_timer% start_timer()
      matter_objects_formatted_file_loop: DO i_matter= 1, parts% n_matter, 1

        ASSOCIATE( nline_in   => header_lines + 1 + &
                                 npart_i_tmp(i_matter)*(i_matter-1) + 1, &
                   nline_fin  => header_lines + 1 + &
                                 npart_i_tmp(i_matter) + &
                                 npart_i_tmp(i_matter)*(i_matter-1) )

          ! Determine boundaries of the lattices
          xmin= ABS(center(i_matter, 1)) + sizes(i_matter, 1)
          xmax= ABS(center(i_matter, 1)) + sizes(i_matter, 2)
          ymin= ABS(center(i_matter, 2)) + sizes(i_matter, 3)
          ymax= ABS(center(i_matter, 2)) + sizes(i_matter, 4)
          zmin= ABS(center(i_matter, 3)) + sizes(i_matter, 5)
          zmax= ABS(center(i_matter, 3)) + sizes(i_matter, 6)

          CALL parts% read_particles_formatted_file &
            ( unit_pos, nline_in, nline_fin, &
              xmin, xmax, ymin, ymax, zmin, zmax, &
              parts_all(i_matter)% pos_i, &
              parts_all(i_matter)% pvol_i, &
              parts_all(i_matter)% nu_i, &
              parts_all(i_matter)% h_i )

          IF(.NOT.parts% read_nu) parts_all(i_matter)% nu_i= &
              parts% masses(i_matter)/npart_i_tmp(i_matter)

          ! Now that the real particle numbers are known, reallocate the arrays
          ! to the appropriate sizes. Note that, if the APM is performed,
          ! this step will be done after it as well
          ! TODO: maybe it is not necessary for the arrays pvol_i, nu_i and h_i?
          parts_all(i_matter)% pos_i = &
                    parts_all(i_matter)% pos_i( :, 1:parts% npart_i(i_matter) )
          parts_all(i_matter)% pvol_i = &
                      parts_all(i_matter)% pvol_i( 1:parts% npart_i(i_matter) )
          parts_all(i_matter)% nu_i = &
                      parts_all(i_matter)% nu_i( 1:parts% npart_i(i_matter) )
          parts_all(i_matter)% h_i = &
                      parts_all(i_matter)% h_i( 1:parts% npart_i(i_matter) )

          CALL impose_equatorial_plane_symmetry &
            ( npart_i_tmp(i_matter), parts_all(i_matter)% pos_i, &
                                     parts_all(i_matter)% nu_i )

          PRINT *, " * Maximum n. baryon per particle (nu) on object", &
                   i_matter, "=", MAXVAL( parts_all(i_matter)% nu_i, DIM= 1 )
          PRINT *, " * Minimum n. baryon per particle (nu) on object", &
                   i_matter, "=", MINVAL( parts_all(i_matter)% nu_i, DIM= 1 )
          PRINT *, " * Ratio between the two=", &
                   MAXVAL( parts_all(i_matter)% nu_i, DIM= 1 )&
                  /MINVAL( parts_all(i_matter)% nu_i, DIM= 1 )
          PRINT *

          two_matter_objects_read: &
          IF( i_matter == 1 .AND. parts% n_matter == 2 )THEN

            ! with practically the same mass, and the physical system
            ! is symmetric wrt the yz plane (in which case the user should set
            ! the reflect_particles_x to .TRUE. in the parameter file)
            equal_masses_read: &
            IF( ABS(parts% masses(1) - parts% masses(2)) &
               /parts% masses(2) <= tol_equal_mass .AND. reflect_particles_x )THEN

              ! ...reflect particles

              DEALLOCATE(parts_all(2)% pos_i)
              DEALLOCATE(parts_all(2)% pvol_i)
              DEALLOCATE(parts_all(2)% h_i)
              DEALLOCATE(parts_all(2)% nu_i)

              CALL reflect_particles_yz_plane( parts_all(1)% pos_i,   &
                                               parts_all(1)% pvol_i,  &
                                               parts_all(1)% nu_i,    &
                                               parts_all(1)% h_i,     &
                                               parts% npart_i(1),     &
                                               parts_all(2)% pos_i,   &
                                               parts_all(2)% pvol_i,  &
                                               parts_all(2)% nu_i,    &
                                               parts_all(2)% h_i,     &
                                               parts% npart_i(2) )



              PRINT *, "** Particles placed on star 1, read from formatted ", &
                       " file and reflected about the yz plane onto star 2."
              PRINT *

              EXIT

            ENDIF equal_masses_read

          ENDIF two_matter_objects_read

        END ASSOCIATE

      ENDDO matter_objects_formatted_file_loop
      CALL parts% placer_timer% stop_timer()

      CLOSE( unit= unit_pos )

      IF( debug )THEN
        PRINT *, "parts% npart_i_tmp=", npart_i_tmp
        PRINT *, "parts% npart_i=", parts% npart_i
        PRINT *, "parts% npart=", parts% npart
        PRINT *
      ENDIF

    END SUBROUTINE read_particles_from_formatted_file


    SUBROUTINE place_particles_on_lattices

      !**************************************************************
      !
      !# Place particles on lattices, one per matter object, and
      !  reflect  particles with respect to the yz plane in the case
      !  of equal-mass binaries
      !
      !  FT 24.10.2022
      !
      !**************************************************************

      IMPLICIT NONE

      PRINT *, " * Placing particles on lattices, one around each ", &
               "matter object."
      PRINT *

      ! Place particles, and time the process

      CALL parts% placer_timer% start_timer()
      matter_objects_lattices_loop: DO i_matter= 1, parts% n_matter, 1

        ! Determine boundaries of the lattices
        xmin= center(i_matter, 1) - stretch*sizes(i_matter, 1)
        xmax= center(i_matter, 1) + stretch*sizes(i_matter, 2)
        ymin= center(i_matter, 2) - stretch*sizes(i_matter, 3)
        ymax= center(i_matter, 2) + stretch*sizes(i_matter, 4)
        zmin= center(i_matter, 3) - stretch*sizes(i_matter, 5)
        zmax= center(i_matter, 3) + stretch*sizes(i_matter, 6)

        central_density(i_matter)= id% read_mass_density &
          ( center(i_matter, 1), center(i_matter, 2), center(i_matter, 3) )

        CALL parts% place_particles_lattice( central_density(i_matter), &
                                             xmin, xmax, ymin, &
                                             ymax, zmin, zmax, &
                                             npart_des_i(i_matter), &
                                             parts% npart_i(i_matter), &
                                             stretch, thres, &
                                             parts_all(i_matter)% pos_i, &
                                             parts_all(i_matter)% pvol_i, &
                                             parts_all(i_matter)% nu_i, &
                                             parts_all(i_matter)% h_i, &
                                             import_density, &
                                             import_id, &
                                             validate_position )

        ! Now that the real particle numbers are known, reallocate the arrays
        ! to the appropriate sizes. Note that, if the APM is performed,
        ! this step will be done after it as well
        ! TODO: maybe this is not necessary for the arrays pvol_i, nu_i and h_i?
        parts_all(i_matter)% pos_i = &
                    parts_all(i_matter)% pos_i( :, 1:parts% npart_i(i_matter) )
        parts_all(i_matter)% pvol_i = &
                    parts_all(i_matter)% pvol_i( 1:parts% npart_i(i_matter) )
        parts_all(i_matter)% nu_i = &
                    parts_all(i_matter)% nu_i( 1:parts% npart_i(i_matter) )
        parts_all(i_matter)% h_i = &
                    parts_all(i_matter)% h_i( 1:parts% npart_i(i_matter) )

        ! If there are 2 matter objects...
        equal_masses_lattices: &
        IF( i_matter == 1 .AND. parts% n_matter == 2 )THEN

          ! ...with practically the same mass, and the physical system is
          ! symmetric wrt the yz plane (in which case the user should
          ! set reflect_particles_x in the parameter file)...
          IF( ABS(parts% masses(1) - parts% masses(2)) &
            /parts% masses(2) <= tol_equal_mass .AND. reflect_particles_x )THEN

            CALL reflect_particles_yz_plane( parts_all(1)% pos_i,   &
                                             parts_all(1)% pvol_i,  &
                                             parts_all(1)% nu_i,    &
                                             parts_all(1)% h_i,     &
                                             parts% npart_i(1),     &
                                             parts_all(2)% pos_i,   &
                                             parts_all(2)% pvol_i,  &
                                             parts_all(2)% nu_i,    &
                                             parts_all(2)% h_i,     &
                                             parts% npart_i(2) )

            EXIT

          ENDIF

        ENDIF equal_masses_lattices

      ENDDO matter_objects_lattices_loop
      CALL parts% placer_timer% stop_timer()

      parts% npart= SUM( parts% npart_i )

      IF( debug ) PRINT *, "10"


    END SUBROUTINE place_particles_on_lattices


    SUBROUTINE place_particles_on_ellipsoidal_surfaces

      !**************************************************************
      !
      !# Place particles on ellipsoidal surfaces, and
      !  reflect  particles with respect to the yz plane in the case
      !  of equal-mass binaries
      !
      !  FT 24.10.2022
      !
      !**************************************************************

      IMPLICIT NONE

      CHARACTER(LEN=:), ALLOCATABLE:: surface_geometry

      PRINT *, "** Placing equal-mass particles on surfaces, " &
               // "taking into account the mass profile of the stars."
      PRINT *

      ! Here the particle mass is computed using the radial mass profile
      ! of the star, so nu should not be redistributed to achieve a given
      ! particle mass ratio
      !    IF( parts% redistribute_nu .EQV. .TRUE. )THEN
      !        parts% redistribute_nu= .FALSE.
      !    ENDIF

      ! Place particles, and time the process
      CALL parts% placer_timer% start_timer()

      matter_objects_sphersurfaces_loop: DO i_matter= 1, parts% n_matter, 1

        IF( i_matter <= 9 ) WRITE( str_i, '(I1)' ) i_matter
        IF( i_matter >= 10 .AND. parts% n_matter <= 99 ) &
                                              WRITE( str_i, '(I2)' ) i_matter
        IF( i_matter >= 100 .AND. parts% n_matter <= 999 ) &
                                              WRITE( str_i, '(I3)' ) i_matter

        IF( parts% surfaces(i_matter)% is_known )THEN

          surface_geometry="oval_surfaces"

        ELSE

          surface_geometry="ellipsoidal_surfaces"

        ENDIF

        filename_mass_profile= TRIM(sph_path)//TRIM(surface_geometry) &
          //"_mass_profile"//TRIM(str_i)//".dat"
        filename_shells_radii= TRIM(sph_path)//TRIM(surface_geometry) &
          //"_radii"//TRIM(str_i)//".dat"
        filename_shells_pos  = TRIM(sph_path)//TRIM(surface_geometry) &
          //"_pos"//TRIM(str_i)//".dat"

        CALL parts% place_particles_ellipsoidal_surfaces( &
                                              parts% masses(i_matter), &
                                              MAXVAL(sizes(i_matter,1:2)), &
                                              center(i_matter,:), &
                                              central_density(i_matter), &
                                              npart_des_i(i_matter), &
                                              parts% npart_i(i_matter), &
                                              parts_all(i_matter)% pos_i, &
                                              parts_all(i_matter)% pvol_i, &
                                              parts_all(i_matter)% nu_i, &
                                              parts_all(i_matter)% h_i, &
                                              last_r, &
                                              upper_bound, lower_bound, &
                                              upper_factor, lower_factor,&
                                              max_steps, &
                                              filename_mass_profile, &
                                              filename_shells_radii, &
                                              filename_shells_pos, &
                                              import_density, &
                                              integrate_mass_density, &
                                              import_id, &
                                        validate_position= validate_position, &
            radii= [MAXVAL(sizes(i_matter,3:4)),MAXVAL(sizes(i_matter,5:6))], &
            surf= parts% surfaces(i_matter) )

        ! Now that the real particle numbers are known, reallocate the arrays
        ! to the appropriate sizes. Note that, if the APM is performed,
        ! this step will be done after it as well
        ! TODO: maybe this is not necessary for the arrays pvol_i, nu_i and h_i?
      !   parts_all(itr)% pos_i = &
      !                     parts_all(itr)% pos_i( :, 1:parts% npart_i(itr) )
      !   parts_all(itr)% pvol_i = &
      !                     parts_all(itr)% pvol_i( 1:parts% npart_i(itr) )
      !   parts_all(itr)% nu_i = &
      !                     parts_all(itr)% nu_i( 1:parts% npart_i(itr) )
      !   parts_all(itr)% h_i = &
      !                     parts_all(itr)% h_i( 1:parts% npart_i(itr) )

        ! If there are 2 matter objects...
        equal_masses: IF( i_matter == 1 .AND. parts% n_matter == 2 )THEN

          ! ...with practically the same mass, and the physical system is
          ! symmetric wrt the yz plane (in which case the user should
          ! set reflect_particles_x in the parameter file)...
          IF( ABS(parts% masses(1) - parts% masses(2)) &
              /parts% masses(2) <= tol_equal_mass .AND. reflect_particles_x )THEN

            CALL reflect_particles_yz_plane( parts_all(1)% pos_i,   &
                                             parts_all(1)% pvol_i,  &
                                             parts_all(1)% nu_i,    &
                                             parts_all(1)% h_i,     &
                                             parts% npart_i(1),     &
                                             parts_all(2)% pos_i,   &
                                             parts_all(2)% pvol_i,  &
                                             parts_all(2)% nu_i,    &
                                             parts_all(2)% h_i,     &
                                             parts% npart_i(2) )

            EXIT

          ENDIF

        ENDIF equal_masses

        !STOP

      ENDDO matter_objects_sphersurfaces_loop
      CALL parts% placer_timer% stop_timer()

      DO i_matter= 1, parts% n_matter, 1

        parts_all(i_matter)% pos_i = &
                    parts_all(i_matter)% pos_i( :, 1:parts% npart_i(i_matter) )
        parts_all(i_matter)% pvol_i = &
                    parts_all(i_matter)% pvol_i( 1:parts% npart_i(i_matter) )
        parts_all(i_matter)% nu_i = &
                    parts_all(i_matter)% nu_i( 1:parts% npart_i(i_matter) )
      ENDDO

      parts% npart= SUM( parts% npart_i )


    END SUBROUTINE place_particles_on_ellipsoidal_surfaces


    SUBROUTINE compute_nstar_eul_id( npart, &
                                     v_euler_x, v_euler_y, v_euler_z, &
                                     g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, &
                                     baryon_density, nstar_eul_id )

      !**************************************************************
      !
      !# Compute nstar_eul_id, the relativistic baryon mass density
      !  seen by the Eulerian observer, given the |id|
      !
      !  FT 31.08.2021
      !
      !**************************************************************

      USE tensor,   ONLY: jx, jy, jz, n_sym4x4
      USE utility,  ONLY: compute_g4, determinant_sym3x3, &
                          spatial_vector_norm_sym3x3, zero, one, two

      IMPLICIT NONE

      INTEGER, INTENT(IN):: npart
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: v_euler_x
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: v_euler_y
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: v_euler_z
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: g_xx
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: g_xy
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: g_xz
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: g_yy
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: g_yz
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: g_zz
      DOUBLE PRECISION, DIMENSION(npart), INTENT(IN):: baryon_density
      DOUBLE PRECISION, DIMENSION(npart), INTENT(OUT):: nstar_eul_id

      INTEGER:: a, i!mus, nus
      DOUBLE PRECISION:: det, sq_g, v_euler_norm2, gamma_eul_a
      DOUBLE PRECISION, DIMENSION(0:3,npart):: vel
      !DOUBLE PRECISION:: g4(0:3,0:3)
      DOUBLE PRECISION:: g4(n_sym4x4)

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( npart, &
      !$OMP                     v_euler_x, v_euler_y, v_euler_z, &
      !$OMP                     g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, &
      !$OMP                     baryon_density, nstar_eul_id ) &
      !$OMP             PRIVATE( a, det, sq_g, v_euler_norm2, gamma_eul_a )
      DO a= 1, npart, 1

        CALL determinant_sym3x3( &
              [g_xx(a),g_xy(a),g_xz(a),g_yy(a),g_yz(a),g_zz(a)], det )

        IF( ABS(det) < 1D-10 )THEN
          PRINT *, "ERROR! The determinant of the spatial metric is " &
                   // "effectively 0 at particle ", a
          PRINT *
          STOP
        ELSEIF( det < 0 )THEN
          PRINT *, "ERROR! The determinant of the spatial metric is " &
                   // "negative at particle ", a
          PRINT *
          STOP
        ELSEIF( .NOT.is_finite_number(det) )THEN
          PRINT *, "ERROR! The determinant is ", det, "at particle ", a
          PRINT *
          STOP
        ENDIF
        sq_g= SQRT(det)

        !
        !-- Generalized Lorentz factor
        !
        v_euler_norm2= zero
        CALL spatial_vector_norm_sym3x3( &
             [g_xx(a),g_xy(a),g_xz(a),g_yy(a),g_yz(a),g_zz(a)], &
             [v_euler_x(a),v_euler_y(a),v_euler_z(a)], v_euler_norm2 )
        IF( .NOT.is_finite_number(v_euler_norm2) )THEN
          PRINT *, "** ERROR! The spatial norm of v_euler is ", v_euler_norm2, &
                   "at particle ", a, &
                   "in SUBROUTINE compute_nstar_eul_id"
          PRINT *, " * Stopping..."
          PRINT *
          STOP
        ENDIF

        gamma_eul_a= one/SQRT(one - v_euler_norm2)
        IF( .NOT.is_finite_number(gamma_eul_a) )THEN
          PRINT *, "** ERROR! The Lorentz factor is ", gamma_eul_a, &
                   "at particle ", a, &
                   "in SUBROUTINE compute_nstar_eul_id"
          PRINT *, " * Stopping..."
          PRINT *
          STOP
        ENDIF
        IF( gamma_eul_a < one )THEN
          PRINT *, "** ERROR! The Lorentz factor is ", gamma_eul_a, &
                   "< 1 at particle ", a, &
                   "in SUBROUTINE compute_nstar_eul_id"
          PRINT *, " * Stopping..."
          PRINT *
          STOP
        ENDIF

        nstar_eul_id(a)= sq_g*gamma_eul_a*baryon_density(a)

      ENDDO
      !$OMP END PARALLEL DO

    END SUBROUTINE compute_nstar_eul_id


    SUBROUTINE compute_pressure( npart, x, y, z, nlrf, eqos, pressure, verbose )
    !! Wrapper function to compute the pressure from the given input

      IMPLICIT NONE

      INTEGER,          INTENT(IN)   :: npart
      !! Returns the baryon mass density at the desired point
      DOUBLE PRECISION, INTENT(IN)   :: x(npart)
      !! \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN)   :: y(npart)
      !! \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN)   :: z(npart)
      !! \(z\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN)   :: nlrf(npart)
      !! Baryon mass density in the local rest frame
      TYPE(eos),        INTENT(IN)   :: eqos
      !! |eos| to use
      DOUBLE PRECISION, INTENT(INOUT):: pressure(npart)
      !! Pressure at \((x,y,z)\)
      LOGICAL, INTENT(IN), OPTIONAL:: verbose

      DOUBLE PRECISION, DIMENSION(npart):: tmp, tmp2, tmp3
      LOGICAL:: verb

      IF(PRESENT(verbose))THEN
        verb= verbose
      ELSE
        verb=.TRUE.
      ENDIF

      CALL parts% compute_sph_hydro( 1, npart, &
        eqos, nlrf, tmp, pressure, tmp2, tmp3, verb )

    END SUBROUTINE compute_pressure


    SUBROUTINE reflect_particles_yz_plane( pos_star1, pvol_star1, &
                                           nu_star1, h_star1, npart_star1, &
                                           pos_star2, pvol_star2, &
                                           nu_star2, h_star2, npart_star2 )

      !**************************************************************
      !
      !# Reflect particles of star 1
      !  with respect to the \(yz\) plane and place them on star 2
      !
      !  FT 07.02.2022
      !
      !**************************************************************

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN):: pos_star1
      !! Array where to store the particle positions for star 1
      DOUBLE PRECISION, DIMENSION(:),   INTENT(IN):: pvol_star1
      !! Array where to store the particle volumes for star 1
      DOUBLE PRECISION, DIMENSION(:),   INTENT(IN):: nu_star1
      !! Array where to store the particle baryon number for star 1
      DOUBLE PRECISION, DIMENSION(:),   INTENT(IN):: h_star1
      !! Array where to store the particle smoothing lengths for star 1
      INTEGER,                          INTENT(IN):: npart_star1
      !! Variable where to store the particle number for star 1
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT):: pos_star2
      !! Array where to store the particle positions for star 2
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, INTENT(INOUT):: pvol_star2
      !! Array where to store the particle volumes for star 2
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, INTENT(INOUT):: nu_star2
      !! Array where to store the particle baryon number for star 2
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, INTENT(INOUT):: h_star2
      !! Array where to store the particle smoothing lengths for star 2
      INTEGER,                                       INTENT(INOUT):: npart_star2
      !! Variable where to store the particle number for star 2


      IF(ALLOCATED(pos_star2)) DEALLOCATE(pos_star2)
      ALLOCATE( pos_star2( 3, npart_star1 ), STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array pos_star2 in SUBROUTINE" &
                 // "reflect_particles_yz_plane. ", &
                 "The error message is", err_msg
        STOP
      ENDIF

      IF(ALLOCATED(pvol_star2)) DEALLOCATE(pvol_star2)
      ALLOCATE( pvol_star2( npart_star1 ), STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array pvol_star2 in SUBROUTINE" &
                 // "reflect_particles_yz_plane. ", &
                 "The error message is", err_msg
        STOP
      ENDIF

      IF(ALLOCATED(nu_star2)) DEALLOCATE(nu_star2)
      ALLOCATE( nu_star2( npart_star1 ), STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array nu_star2 in SUBROUTINE" &
                 // "reflect_particles_yz_plane. ", &
                 "The error message is", err_msg
        STOP
      ENDIF

      IF(ALLOCATED(h_star2)) DEALLOCATE(h_star2)
      ALLOCATE( h_star2( npart_star1 ), STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array h_star2 in SUBROUTINE" &
                 // "reflect_particles_yz_plane. ", &
                 "The error message is", err_msg
        STOP
      ENDIF

      PRINT *, " * Reflecting particles with respect to the yz plane..."
      PRINT *

      ! Reflect the particles on matter object 1, and their properties,
      ! to matter object 2
      npart_star2= npart_star1
      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( pos_star1, pos_star2, pvol_star1, pvol_star2, &
      !$OMP                     nu_star1, nu_star2, h_star1, h_star2, &
      !$OMP                     npart_star2 ) &
      !$OMP             PRIVATE( a )
      DO a= 1, npart_star2, 1

        pos_star2(1,a)= - pos_star1(1,a)
        pos_star2(2,a)=   pos_star1(2,a)
        pos_star2(3,a)=   pos_star1(3,a)
        pvol_star2 (a)=   pvol_star1(a)
        nu_star2   (a)=   nu_star1  (a)
        h_star2    (a)=   h_star1   (a)

      ENDDO
      !$OMP END PARALLEL DO


    END SUBROUTINE reflect_particles_yz_plane


  END PROCEDURE construct_particles_std


  MODULE PROCEDURE destruct_particles

    !*********************************************
    !
    !# Destructor of a particles object
    !
    !  FT
    !
    !*********************************************

    IMPLICIT NONE

    CALL this% deallocate_particles_memory()


  END PROCEDURE destruct_particles


END SUBMODULE constructor_std




!  DEPRECATED? This is a relic from when nu was re-assigned to the particles
!
!    IF( this% redistribute_nu )THEN
!
!      !---------------------------------------------------------------------!
!      !--  Assignment of nu on the stars, with the purpose                --!
!      !--  of having a more uniform nu over the particles without losing  --!
!      !--  baryon mass. This is used only on the lattice, optionally.     --!
!      !---------------------------------------------------------------------!
!
!      IF( this% distribution_id == id_particles_on_ellipsoidal_surfaces )THEN
!        PRINT *, "** ERROR! Particle placer ", this% distribution_id, &
!                 " is not compatible with redistribute_nu= .TRUE."
!        PRINT *, " * Check the parameter file lorene_bns_id_particles.par. ", &
!                 "Stopping..."
!        PRINT *
!        STOP
!      ENDIF
!
!      nu_max1= nlrf( this% baryon_density_index( this% npart1 ) )&
!              *this% pvol( this% npart1 ) &
!              *Theta( this% baryon_density_index( this% npart1 ) )&
!              *sq_det_g4( this% baryon_density_index( this% npart1 ) )
!      nu_max2= nlrf( this% baryon_density_index( this% npart ) )&
!              *this% pvol( this% npart ) &
!              *Theta( this% baryon_density_index( this% npart ) )&
!              *sq_det_g4( this% baryon_density_index( this% npart ) )
!
!      nu_thres1= nu_max1/this% nu_ratio
!      nu_thres2= nu_max2/this% nu_ratio
!
!      ! Reset the total baryon number to 0 (necessary), and nu to an arbitrary
!      ! value (to make debugging easier)
!
!      nu= one
!      this% nu= one
!      this% nbar_tot= zero
!      this% nbar1= zero
!      this% nbar2= zero
!
!      cnt1= 0
!      compute_nu_on_particles_star1: DO itr= this% npart1, 1, -1
!
!        cnt1= cnt1 + 1
!
!        nu_tmp= nlrf( this% baryon_density_index( itr ) ) &
!                *this% pvol(itr) &
!                *Theta( this% baryon_density_index( itr ) )&
!                *sq_det_g4( this% baryon_density_index( itr ) )
!
!        !IF( itr == this% npart1 ) nu_max= nu_tmp ! move this out of the loop
!
!        IF( nu_tmp > nu_thres1 )THEN
!          nu( this% baryon_density_index( itr ) )      = nu_tmp
!          this% nu( this% baryon_density_index( itr ) )= nu_tmp
!        ELSE
!          nu( this% baryon_density_index( itr ) )      = nu_thres1
!          this% nu( this% baryon_density_index( itr ) )= nu_thres1
!        ENDIF
!
!        this% nbar1= this% nbar1 + &
!                     this% nu( this% baryon_density_index( itr ) )
!
!        IF( this% nbar1*amu/MSun > this% masses(1) )THEN
!          EXIT
!        ENDIF
!
!      ENDDO compute_nu_on_particles_star1
!
!      cnt2= 0
!      compute_nu_on_particles_star2: DO itr= this% npart, this% npart1 + 1, -1
!
!        cnt2= cnt2 + 1
!
!        nu_tmp= nlrf( this% baryon_density_index( itr ) ) &
!                *this% pvol(itr) &
!                *Theta( this% baryon_density_index( itr ) ) &
!                *sq_det_g4( this% baryon_density_index( itr ) )
!
!        !IF( itr == this% npart ) nu_max= nu_tmp
!
!        IF( nu_tmp > nu_thres2 )THEN
!          nu( this% baryon_density_index( itr ) )      = nu_tmp
!          this% nu( this% baryon_density_index( itr ) )= nu_tmp
!        ELSE
!          nu( this% baryon_density_index( itr ) )      = nu_thres2
!          this% nu( this% baryon_density_index( itr ) )= nu_thres2
!        ENDIF
!
!        this% nbar2= this% nbar2 + &
!                     this% nu( this% baryon_density_index( itr ) )
!
!        IF( this% nbar2*amu/MSun > this% masses(2) )THEN
!          EXIT
!        ENDIF
!
!      ENDDO compute_nu_on_particles_star2
!      this% nbar_tot= this% nbar1 + this% nbar2
!
!      !
!      !-- Reshape MODULE variables
!      !
!
!      CALL this% reshape_sph_field( pos_u, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( vel_u, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( Theta, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( h, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( nlrf, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( u, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( Pr, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( nu, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( temp, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( av, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( divv, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      !
!      !-- Reshape TYPE member SPH variables
!      !
!
!      CALL this% reshape_sph_field( this% pos, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% v, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% v_euler_x, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% v_euler_y, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% v_euler_z, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% Theta, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% h, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% baryon_density, cnt1, &
!                                    cnt2, this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% nlrf, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% energy_density, cnt1, &
!                                    cnt2, this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% specific_energy, cnt1, &
!                                    cnt2, this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% pressure, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% pressure_sph, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% nu, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% pvol, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      !
!      !-- Reshape TYPE member spacetime variables
!      !
!
!      CALL this% reshape_sph_field( this% lapse, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% shift_x, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% shift_y, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% shift_z, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% g_xx, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% g_xy, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% g_xz, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% g_yy, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% g_yz, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      CALL this% reshape_sph_field( this% g_zz, cnt1, cnt2, &
!                                    this% baryon_density_index )
!
!      !
!      !-- Reassign particle numbers
!      !
!
!      npart= cnt1 + cnt2
!      this% npart= npart
!      this% npart1= cnt1
!      this% npart2= cnt2
!      n1= this% npart1
!      n2= this% npart2
!
!      PRINT *, " * Particles replaced after reassigning nu."
!      PRINT *, " * New number of particles=", this% npart
!      PRINT *
!      PRINT *, " * Number of particles on NS 1=", this% npart1
!      PRINT *, " * Number of particles on NS 2=", this% npart2
!      PRINT *

!  DEPRECATED? This is a relic from when nu was re-assigned to the particles
!              It was executed at the end of the constructor
!
!  IF( parts% redistribute_nu )THEN
!
!    ! Index particles on star 1 in increasing order of nu
!
!    CALL indexx( parts% npart1, &
!                 parts% baryon_density( 1 : parts% npart1 ), &
!                 parts% baryon_density_index( 1 : parts% npart1 ) )
!
!    ! Index particles on star 2 in increasing order of nu
!
!    CALL indexx( parts% npart2, &
!                 parts% baryon_density( parts% npart1 + 1 : &
!                                                  parts% npart ), &
!                 parts% baryon_density_index( parts% npart1 + 1 : &
!                                                  parts% npart ) )
!
!    ! Shift indices on star 2 by npart1 since all the arrays store
!    ! the quantities on star 1 first, and then on star 2
!
!    parts% baryon_density_index( parts% npart1 + 1 : &
!                                     parts% npart )= &
!                   parts% npart1 + &
!                   parts% baryon_density_index( parts% npart1 + 1 : &
!                                                    parts% npart )
!
!  ENDIF
