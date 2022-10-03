! File:         submodule_sph_particles_constructor_std.f90
! Authors:      Francesco Torsello (FT)
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


  CONTAINS


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
  !END PROCEDURE construct_particles_idase_empty


  MODULE PROCEDURE construct_particles_std

    !**************************************************
    !
    !# The constructor performs all the tasks needed
    !  to set up the particle distribution with the
    !  |id| evaluated on it. It calls all the PROCEDURES
    !  that rely on an object of STATIC TYPE idbase.
    !
    !  @todo assign sub-tasks to separate SUBROUTINES
    !        contained in this SUBMODULE
    !
    !  FT 17.10.2020
    !
    !  @note Last updated: FT 25.04.2022
    !
    !**************************************************

    !USE NaNChecker,    ONLY: Check_Array_for_NAN
    USE constants,          ONLY: amu, pi, half, third
    USE utility,            ONLY: zero, one, two, ten, Msun_geo, flag$sph
    !USE NR,                 ONLY: indexx
    USE kernel_table,       ONLY: ktable
    USE input_output,       ONLY: read_options
    USE units,              ONLY: set_units
    USE options,            ONLY: ikernel, ndes, eos_str, eos_type
    USE alive_flag,         ONLY: alive
    USE analyze,            ONLY: COM
    USE utility,            ONLY: spherical_from_cartesian, &
                                  spatial_vector_norm_sym3x3, sph_path

    IMPLICIT NONE

    INTEGER, PARAMETER:: unit_pos= 2289
    ! Variable storing the number of column where nu is written
    DOUBLE PRECISION, PARAMETER:: tol_equal_mass= 5.0D-3
    ! Tolerance for the difference between the masse of the stars
    ! for a BNS system, to determine if a BNS is equal-mass or not.
    ! If the relative difference between the masses of the stars is lower
    ! than tol_equal_mass, then the BNS is considered equal_mass.

    ! The variable counter counts how many times the PROCEDURE
    ! construct_particles_idase is called
    INTEGER, SAVE:: counter= 1
    INTEGER:: npart_des, a, max_steps, nlines, header_lines, n_cols, &
              npart_tmp, nx_gh, ny_gh, nz_gh, i_matter

    ! Maximum length for strings, and for the number of imported binaries
    INTEGER, PARAMETER:: max_length= 50
    ! APM parameters
    INTEGER:: apm_max_it, max_inc, print_step
    INTEGER:: column_nu
    ! Temporary number of matter objects
    INTEGER:: n_matter_tmp, tmp
    ! Array storing the columns of the file parts_pos (defined below) that
    ! contain the particle positions
    INTEGER, DIMENSION(3):: columns
    INTEGER, DIMENSION(id% get_n_matter()):: npart_des_i
    ! Temporary array storing the number of particles on each matter object
    INTEGER, DIMENSION(:), ALLOCATABLE:: npart_i_tmp

    DOUBLE PRECISION:: thres, nu_ratio_des, ghost_dist
    DOUBLE PRECISION:: xmin, xmax, ymin, ymax, zmin, zmax, stretch
    DOUBLE PRECISION:: upper_bound, lower_bound, upper_factor, lower_factor, &
                       last_r
    DOUBLE PRECISION:: pvol_tmp
    DOUBLE PRECISION:: max_mass, total_mass
    !DOUBLE PRECISION:: ratio_npart_des_real, pmass_des
    DOUBLE PRECISION:: min_eps, min_vel, theta_a, phi_a, r_a!, rad_part

    DOUBLE PRECISION, DIMENSION(id% get_n_matter())  :: central_density
    DOUBLE PRECISION, DIMENSION(id% get_n_matter(),3):: center
    DOUBLE PRECISION, DIMENSION(id% get_n_matter(),3):: barycenter
    DOUBLE PRECISION, DIMENSION(id% get_n_matter(),6):: sizes
    DOUBLE PRECISION, DIMENSION(id% get_n_matter())  :: lapse_lengthscales
    DOUBLE PRECISION, DIMENSION(id% get_n_matter())  :: g00_lengthscales

    DOUBLE PRECISION, DIMENSION( :, : ), ALLOCATABLE:: tmp_pos
    DOUBLE PRECISION:: nuratio_thres, nuratio_des
    DOUBLE PRECISION:: min_lapse, min_g00_abs, shift_norm
    !DOUBLE PRECISION:: com_x, com_y, com_z, com_d

    TYPE parts_i
      DOUBLE PRECISION, DIMENSION( :, : ), ALLOCATABLE:: pos_i
      DOUBLE PRECISION, DIMENSION( : ),    ALLOCATABLE:: pvol_i
      DOUBLE PRECISION, DIMENSION( : ),    ALLOCATABLE:: pmass_i
      DOUBLE PRECISION, DIMENSION( : ),    ALLOCATABLE:: h_i
      DOUBLE PRECISION, DIMENSION( : ),    ALLOCATABLE:: nu_i
    END TYPE

    TYPE(parts_i), DIMENSION(id% get_n_matter()):: parts_all

    ! String storing the name of the directory storing the files containing
    ! the particle distributions
    CHARACTER( LEN= max_length ):: parts_pos_path
    ! String storing the name of the file containing the particle positions
    CHARACTER( LEN= max_length ):: parts_pos
    ! Final name for the file containing the particle positions
    CHARACTER( LEN= : ), ALLOCATABLE:: parts_pos_namefile
    CHARACTER( LEN= 3 ):: str_i
    ! String storing the local path to the directory where the
    ! |lorene| BNS ID files are stored
    CHARACTER( LEN= max_length ):: compose_path
    ! String storing the names of the |lorene| BNS ID binary files
    CHARACTER( LEN= max_length ):: compose_filename

    CHARACTER( LEN= max_length ):: filename_apm_pos_id, filename_apm_pos, &
                                   filename_apm_results

    CHARACTER( LEN= max_length ):: filename_mass_profile, &
                                   filename_shells_radii, filename_shells_pos

    LOGICAL:: file_exists, use_thres, redistribute_nu, correct_nu, &
              compose_eos, exist, randomize_phi, randomize_theta, &
              randomize_r, mass_it, &
              read_nu, reflect_particles_x

    LOGICAL, PARAMETER:: debug= .FALSE.

    LOGICAL, DIMENSION(id% get_n_matter()):: apm_iterate, use_atmosphere, &
                                             remove_atmosphere

    NAMELIST /sphincs_id_particles/ &
              parts_pos_path, parts_pos, columns, header_lines, n_cols, &
              read_nu, column_nu, &
              stretch, &
              use_thres, thres, nu_ratio_des, redistribute_nu, correct_nu, &
              compose_eos, compose_path, compose_filename, &
              npart_des, last_r, upper_bound, lower_bound, &
              upper_factor, lower_factor, max_steps, &
              randomize_phi, randomize_theta, randomize_r, &
              apm_iterate, apm_max_it, max_inc, mass_it, &
              nuratio_thres, reflect_particles_x, nx_gh, ny_gh, nz_gh, &
              use_atmosphere, remove_atmosphere, nuratio_des, print_step, &
              ghost_dist

    ! Get the number of matter objects in the physical system
    parts% n_matter= id% get_n_matter()

    ! Get the the logical variable at specifies if the system is cold
    ! (no thermal component)
    parts% cold_system= id% get_cold_system()

    ALLOCATE( parts% apm_timers(parts% n_matter) )

    !
    !-- Initialize the timers
    !
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
    !-- Read the parameters of the particle distributions
    !
    parts% sphincs_id_particles= 'sphincs_id_particles.dat'

    INQUIRE( FILE= parts% sphincs_id_particles, EXIST= file_exists )
    IF( file_exists )THEN
      OPEN( 10, FILE= parts% sphincs_id_particles, STATUS= 'OLD' )
    ELSE
      PRINT *
      PRINT *, "** ERROR: ", parts% sphincs_id_particles, &
              " file not found!"
      PRINT *
      STOP
    ENDIF

    READ( 10, NML= sphincs_id_particles )
    CLOSE( 10 )

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
    ALLOCATE( parts% apm_iterate( parts% n_matter ) )
    parts% apm_iterate   = apm_iterate
    parts% use_atmosphere= use_atmosphere
    parts% read_nu       = read_nu

    parts_pos_namefile= TRIM(parts_pos_path)//TRIM(parts_pos)

    !
    !-- Check that the parameters are acceptable
    !

    IF( upper_bound <= lower_bound )THEN
      PRINT *
      PRINT *, "** ERROR in lorene_bns_id_particles.par: ", &
               "upper_bound should be greater than lower_bound!"
      PRINT *
      STOP
    ENDIF
    IF( upper_factor < 1.0D0 )THEN
      PRINT *
      PRINT *, "** ERROR in lorene_bns_id_particles.par: ", &
               "upper_factor should be greater than or equal to 1!"
      PRINT *
      STOP
    ENDIF
    IF( lower_factor > 1 )THEN
      PRINT *
      PRINT *, "** ERROR in lorene_bns_id_particles.par: ", &
               "lower_factor should be smaller than or equal to 1!"
      PRINT *
      STOP
    ENDIF
    IF( max_steps < 10 )THEN
      PRINT *
      PRINT *, "** ERROR in lorene_bns_id_particles.par: ", &
               "max_steps should be an integer greater than or equal to 10!"
      PRINT *
      STOP
    ENDIF
    IF( last_r < 0.95D0 .OR. last_r > 1.0D0 )THEN
      PRINT *
      PRINT *, "** ERROR in lorene_bns_id_particles.par: ", &
               "last_r should be greater than or equal to 0.95, ", &
               "and lower than or equal to 1!"
      PRINT *
      STOP
    ENDIF
    IF( apm_max_it < 0 .OR. max_inc < 0 .OR. nuratio_thres < 0 &
        .OR. nuratio_des < 0 .OR. nx_gh < 0 .OR. ny_gh < 0 .OR. nz_gh < 0 )THEN
      PRINT *
      PRINT *, "** ERROR in lorene_bns_id_particles.par: ", &
               "the numeric parameters for the APM method should be positive!"
      PRINT *
      STOP
    ENDIF
    IF( nuratio_des >= nuratio_thres )THEN
      PRINT *
      PRINT *, "** ERROR in lorene_bns_id_particles.par: ", &
               "nuratio_des has to be stricly lower than nuratio_thres!"
      PRINT *
      STOP
    ENDIF
    IF( print_step < 0 )THEN
      PRINT *
      PRINT *, "** ERROR in sphincs_id_particles.dat: ", &
               "print_step has to be a positive integer or zero!"
      PRINT *
      STOP
    ENDIF
    IF( ghost_dist < zero )THEN
      PRINT *
      PRINT *, "** ERROR in sphincs_id_particles.dat: ", &
               "ghost_dist has to be a positive double precision or zero!"
      PRINT *
      STOP
    ENDIF

    ! setup unit system
    CALL set_units('NSM')
    CALL read_options

    ! tabulate kernel, get ndes
    CALL ktable( ikernel, ndes )

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

    CALL id% initialize_id( flag$sph )

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

    parts% npart_i(0)= 0
    npart_i_tmp(0)   = 0
    parts% nbar_i    = zero
    parts% nuratio_i = zero

    DO i_matter= 1, parts% n_matter, 1

      parts% adm_mass          = id% adm_mass
      parts% masses(i_matter)  = id% return_mass(i_matter)
      center(i_matter,:)       = id% return_center(i_matter)
      central_density(i_matter)= id% read_mass_density( center(i_matter,1), &
                                                        center(i_matter,2), &
                                                        center(i_matter,3) )
      barycenter(i_matter,:)   = id% return_barycenter(i_matter)
      parts% barycenter(i_matter,:)= barycenter(i_matter,:)
      sizes(i_matter, :)       = id% return_spatial_extent(i_matter)

      parts% all_eos(i_matter)% eos_name= id% return_eos_name(i_matter)
      CALL id% return_eos_parameters( i_matter, &
                                      parts% all_eos(i_matter)% eos_parameters )

    ENDDO

    ! Compute desired particle numbers based on mass ratios
    max_mass  = MAXVAL( parts% masses )
    total_mass= SUM( parts% masses )
    DO i_matter= 1, parts% n_matter, 1
      parts% mass_ratios(i_matter)   = parts% masses(i_matter)/max_mass
      parts% mass_fractions(i_matter)= parts% masses(i_matter)/total_mass
      npart_des_i(i_matter)          = &
                          NINT(parts% mass_fractions(i_matter)*DBLE(npart_des))
      tmp= 2*npart_des_i(i_matter)
      ALLOCATE( parts_all(i_matter)% pos_i  ( 3, tmp ) )
      ALLOCATE( parts_all(i_matter)% pvol_i ( tmp ) )
      ALLOCATE( parts_all(i_matter)% pmass_i( tmp ) )
      ALLOCATE( parts_all(i_matter)% h_i    ( tmp ) )
      ALLOCATE( parts_all(i_matter)% nu_i   ( tmp ) )
    ENDDO
    !   IF( parts% redistribute_nu )THEN
    !     thres= 100.0D0*parts% nu_ratio
    !   ENDIF

    parts% post_process_sph_id => id% finalize_sph_id_ptr

    !PRINT *, parts% all_eos(1)% eos_parameters
    !STOP

    DO i_matter= 1, parts% n_matter, 1

      IF( parts% all_eos(i_matter)% eos_parameters(1) == DBLE(1) )THEN

        IF( eos_type == 'pwp' )THEN
          PRINT *, "** ERROR! On matter object ", i_matter, &
                   ", the EOS taken from the ID is not the same as the ",&
                   "one specified in parameter file SPHINCS_fm_input.dat."
          PRINT *
          PRINT *, " * EOS from the ID: ", &
                   parts% all_eos(i_matter)% eos_name
          PRINT *, " * EOS from the parameter file SPHINCS_fm_input.dat: ", &
                   eos_type
          PRINT *, "Stopping..."
          PRINT *
          STOP
        ENDIF

      ENDIF

      IF( parts% all_eos(i_matter)% eos_parameters(1) == DBLE(110) )THEN

        IF( eos_type == 'Poly' )THEN
          PRINT *, "** ERROR! On matter object ", i_matter, &
                   ", the EOS taken from the ID is not the same as the ",&
                   "one specified in parameter file SPHINCS_fm_input.dat."
          PRINT *
          PRINT *, " * EOS from the ID: ", &
                   parts% all_eos(i_matter)% eos_name
          PRINT *, " * EOS from the parameter file SPHINCS_fm_input.dat: ", &
                   eos_type
          PRINT *, "Stopping..."
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
          PRINT *, "Stopping..."
          PRINT *
          STOP

        ENDIF

      ENDIF

    ENDDO


    ! TODO: Add check that the number of rows in placer is the same as the
    !       number of bns objects, and that all bns have a value for placer

    !
    !-- Choose particle placer
    !

    choose_particle_placer: SELECT CASE( dist )

    CASE( id_particles_from_file )
    ! Read particles from formatted file

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

      ! Allocate the temporary array, with fixed size, to store data
      ALLOCATE( tmp_pos( n_cols, 2*npart_tmp ) )
      tmp_pos= 0.0D0

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
                 " do not add up to the total number of particles, in file", &
                 TRIM(parts_pos_namefile)
        PRINT *, "   Stopping..."
        PRINT *
        STOP
      ENDIF

      ! Read the data into the temporary array
      DO itr= 1, npart_tmp, 1

        READ( UNIT= unit_pos, FMT= *, IOSTAT = ios, IOMSG= err_msg ) &
          tmp_pos( :, itr )

        IF( ios > 0 )THEN
          PRINT *, "...error when reading " // TRIM(parts_pos_namefile), &
                  " at particle ", itr,". The status variable is ", ios, &
                  ". The error message is", err_msg
          STOP
        ENDIF

      ENDDO

      CLOSE( UNIT= unit_pos )


      ! Check that the positions are within the sizes of the matter objects.
      ! This checks that the positions read from the formatted
      ! file are compatible with the sizes given by the idbase object

      DO itr= 1, parts% n_matter, 1

        IF( .NOT.use_atmosphere(itr) )THEN

          ASSOCIATE( npart_in   => npart_i_tmp(itr)*(itr-1) + 1, &
                     npart_fin  => npart_i_tmp(itr) + npart_i_tmp(itr)*(itr-1) )

          !  PRINT *, ABS( MINVAL( tmp_pos(1,npart_in:npart_fin) ) ) > &
          !          ABS(center(itr,1)) + sizes(itr, 1)
          !  PRINT *, ABS( MAXVAL( tmp_pos(1,npart_in:npart_fin) ) ) > &
          !          ABS(center(itr,1)) + sizes(itr, 2)
          !  PRINT *, ABS( MINVAL( tmp_pos(2,npart_in:npart_fin) ) ) > &
          !          ABS(center(itr,2)) + sizes(itr, 3)
          !  PRINT *, ABS( MAXVAL( tmp_pos(2,npart_in:npart_fin) ) ) > &
          !          ABS(center(itr,2)) + sizes(itr, 4)
          !  PRINT *, ABS( MINVAL( tmp_pos(3,npart_in:npart_fin) ) ) > &
          !          ABS(center(itr,3)) + sizes(itr, 5)
          !  PRINT *, ABS( MAXVAL( tmp_pos(3,npart_in:npart_fin) ) ) > &
          !          ABS(center(itr,3)) + sizes(itr, 6)
          !
          !  PRINT *, ABS( MINVAL( tmp_pos(1,npart_in:npart_fin) ) )
          !  PRINT *, ABS(center(itr,1)) + sizes(itr, 1)
          !  PRINT *, ABS( MAXVAL( tmp_pos(1,npart_in:npart_fin) ) )
          !  PRINT *, ABS(center(itr,1)) + sizes(itr, 2)
          !  PRINT *, ABS( MINVAL( tmp_pos(2,npart_in:npart_fin) ) )
          !  PRINT *, ABS(center(itr,2)) + sizes(itr, 3)
          !  PRINT *, ABS( MAXVAL( tmp_pos(2,npart_in:npart_fin) ) )
          !  PRINT *, ABS(center(itr,2)) + sizes(itr, 4)
          !  PRINT *, ABS( MINVAL( tmp_pos(3,npart_in:npart_fin) ) )
          !  PRINT *, ABS(center(itr,3)) + sizes(itr, 5)
          !  PRINT *, ABS( MAXVAL( tmp_pos(3,npart_in:npart_fin) ) )
          !  PRINT *, ABS(center(itr,3)) + sizes(itr, 6)

            IF( ABS( MINVAL( tmp_pos(1,npart_in:npart_fin) ) ) > &
                        ABS(center(itr,1)) + sizes(itr, 1) &
                .OR. &
                ABS( MAXVAL( tmp_pos(1,npart_in:npart_fin) ) ) > &
                        ABS(center(itr,1)) + sizes(itr, 2) &
                .OR. &
                ABS( MINVAL( tmp_pos(2,npart_in:npart_fin) ) ) > &
                        ABS(center(itr,2)) + sizes(itr, 3) &
                .OR. &
                ABS( MAXVAL( tmp_pos(2,npart_in:npart_fin) ) ) > &
                        ABS(center(itr,2)) + sizes(itr, 4) &
                .OR. &
                ABS( MINVAL( tmp_pos(3,npart_in:npart_fin) ) ) > &
                        ABS(center(itr,3)) + sizes(itr, 5) &
                .OR. &
                ABS( MAXVAL( tmp_pos(3,npart_in:npart_fin) ) ) > &
                        ABS(center(itr,3)) + sizes(itr, 6) &

            )THEN

              PRINT *, "** ERROR! The positions of the particles on object ", &
                       itr, ", read from file "// TRIM(parts_pos_namefile), &
                       " are not compatible with the ", &
                       "physical system read from file. Stopping..."
              PRINT *
              STOP

            ENDIF

          END ASSOCIATE

        ENDIF
      ENDDO

      IF( debug ) PRINT *, "parts% npart_i_tmp=", npart_i_tmp

      ! Impose equatorial plane symmetry on each object
      ! @TODO: make this optional
      DO itr= 1, parts% n_matter, 1

        ASSOCIATE( npart_in   => npart_i_tmp(itr-1) + 1, &
                   npart_fin  => npart_i_tmp(itr-1) + npart_i_tmp(itr) )

          IF( debug )THEN
            PRINT *, "npart_in=", npart_in
            PRINT *, "npart_fin=", npart_fin
            PRINT *
          ENDIF

          IF( read_nu )THEN

            CALL impose_equatorial_plane_symmetry( npart_i_tmp(itr), &
                                            tmp_pos(1:3,npart_in:npart_fin), &
                                            tmp_pos(4,npart_in:npart_fin) )

          ELSE

            CALL impose_equatorial_plane_symmetry( npart_i_tmp(itr), &
                                            tmp_pos(1:3,npart_in:npart_fin) )

          ENDIF

          !parts% npart_i(itr)= 2*parts% npart_i(itr)
          parts_all(itr)% pos_i= &
                              tmp_pos(1:3,npart_in:npart_in+npart_i_tmp(itr)-1)
          IF( read_nu ) parts_all(itr)% nu_i= &
                              tmp_pos(4,npart_in:npart_in+npart_i_tmp(itr)-1)

        END ASSOCIATE

      ENDDO
      parts% npart_i(1:parts% n_matter)= npart_i_tmp(1:parts% n_matter)
      parts% npart = SUM(parts% npart_i)

      IF( debug )THEN
        PRINT *, "parts% npart_i_tmp=", npart_i_tmp
        PRINT *, "parts% npart_i=", parts% npart_i
        PRINT *, "parts% npart=", parts% npart
        PRINT *
      ENDIF

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
          DEALLOCATE(parts_all(2)% pmass_i)
          DEALLOCATE(parts_all(2)% nu_i)

          CALL reflect_particles_yz_plane( parts_all(1)% pos_i,   &
                                           parts_all(1)% pvol_i,  &
                                           parts_all(1)% pmass_i, &
                                           parts_all(1)% nu_i,    &
                                           parts_all(1)% h_i,     &
                                           parts% npart_i(1),     &
                                           parts_all(2)% pos_i,   &
                                           parts_all(2)% pvol_i,  &
                                           parts_all(2)% pmass_i, &
                                           parts_all(2)% nu_i,    &
                                           parts_all(2)% h_i,     &
                                           parts% npart_i(2) )

          PRINT *, "** Particles placed on star 1 according to the APM,", &
                   " and reflected about the yz plane onto star 2."
          PRINT *

        ENDIF equal_masses_read

      ENDIF two_matter_objects_read

      PRINT *, " * Particle positions read. Number of particles=", &
               parts% npart
      PRINT *
      DO itr= 1, parts% n_matter, 1
        PRINT *, " * Number of particles on matter object ", itr, "=", &
                 parts% npart_i(itr)
      ENDDO
      PRINT *

      !
      !-- Computing volume per particle
      !
   !   IF(.NOT.ALLOCATED( parts% pvol ))THEN
   !     ALLOCATE( parts% pvol( parts% npart ), STAT= ios, &
   !             ERRMSG= err_msg )
   !     IF( ios > 0 )THEN
   !       PRINT *, "...allocation error for array pvol ", &
   !                ". The error message is", err_msg
   !       STOP
   !     ENDIF
   !     !CALL test_status( ios, err_msg, &
   !     !        "...allocation error for array v_euler_z" )
   !   ENDIF                                                  

      ! First guess of the particle volume and mass (the first will be computed
      ! exactly later, as the cube of the exact smoothing length). The particle
      ! volume guess determines the first guess for the smoothing length
      ! The particle mass is computed if nu is read from the file
      DO itr= 1, parts% n_matter, 1

        ASSOCIATE( npart_in   => npart_i_tmp(itr)*(itr-1) + 1, &
                   npart_fin  => npart_i_tmp(itr) + npart_i_tmp(itr)*(itr-1) )

          DEALLOCATE( parts_all(itr)% h_i )
          IF(.NOT.ALLOCATED( parts_all(itr)% h_i ))THEN
            ALLOCATE( parts_all(itr)% h_i( npart_i_tmp(itr) ), &
                      STAT= ios, ERRMSG= err_msg )
            IF( ios > 0 )THEN
              PRINT *, "...allocation error for array h_i ", &
                       ". The error message is", err_msg
              STOP
            ENDIF
            !CALL test_status( ios, err_msg, &
            !        "...allocation error for array v_euler_z" )
          ENDIF
          DEALLOCATE( parts_all(itr)% pvol_i )
          IF(.NOT.ALLOCATED( parts_all(itr)% pvol_i ))THEN
            ALLOCATE( parts_all(itr)% pvol_i( npart_i_tmp(itr) ), &
                      STAT= ios, ERRMSG= err_msg )
            IF( ios > 0 )THEN
              PRINT *, "...allocation error for array pvol_i ", &
                       ". The error message is", err_msg
              STOP
            ENDIF
            !CALL test_status( ios, err_msg, &
            !        "...allocation error for array v_euler_z" )
          ENDIF

          pvol_tmp= 0.D0
         !DO a= npart_in, npart_fin - 1, 1
         !
         !  pvol_tmp= pvol_tmp + ABS( parts% pos(3,a + 1) &
         !                          - parts% pos(3,a) )
         !
         !ENDDO
          DO a= 1, npart_i_tmp(itr) - 1, 1

            pvol_tmp= pvol_tmp + ABS( parts_all(itr)% pos_i(3,a + 1) &
                                    - parts_all(itr)% pos_i(3,a) )

          ENDDO
          pvol_tmp= pvol_tmp/( npart_i_tmp(itr) - 1 )

          !parts% pvol(npart_in:npart_fin)= 2.D0*pvol_tmp**3.D0
          parts_all(itr)% pvol_i= 2.D0*pvol_tmp**3.D0
          DO a= 1, npart_i_tmp(itr), 1
            parts_all(itr)% h_i(a)   = (parts_all(itr)% pvol_i(a))**third
          ENDDO

          !IF( read_nu ) parts% pmass(npart_in:npart_fin)= &
          !              parts% nu(npart_in:npart_fin)*amu
          IF( read_nu ) parts_all(itr)% pmass_i= parts_all(itr)% nu_i*amu

        END ASSOCIATE

        PRINT *, " * Maximum n. baryon per particle (nu) on object", itr, &
                            "=", MAXVAL( parts_all(itr)% nu_i, DIM= 1 )
        PRINT *, " * Minimum n. baryon per particle (nu) on object", itr, &
                            "=", MINVAL( parts_all(itr)% nu_i, DIM= 1 )
        PRINT *, " * Ratio between the two=", &
                 MAXVAL( parts_all(itr)% nu_i, DIM= 1 )&
                /MINVAL( parts_all(itr)% nu_i, DIM= 1 )
        PRINT *

      ENDDO


    CASE( id_particles_on_lattice )


      PRINT *, " * Placing particles on lattices, one around each ", &
               "matter object."
      PRINT *

      ! Place particles, and time the proces

      CALL parts% placer_timer% start_timer()
      matter_objects_lattices_loop: DO itr= 1, parts% n_matter, 1

        ! Determine boundaries of the lattices
        xmin= center(itr, 1) - stretch*sizes(itr, 1)
        xmax= center(itr, 1) + stretch*sizes(itr, 2)
        ymin= center(itr, 2) - stretch*sizes(itr, 3)
        ymax= center(itr, 2) + stretch*sizes(itr, 4)
        zmin= center(itr, 3) - stretch*sizes(itr, 5)
        zmax= center(itr, 3) + stretch*sizes(itr, 6)

        central_density(itr)= id% read_mass_density( center(itr, 1), &
                                                     center(itr, 2), &
                                                     center(itr, 3) )

        CALL parts% place_particles_lattice( central_density(itr), &
                                             xmin, xmax, ymin, &
                                             ymax, zmin, zmax, &
                                             npart_des_i(itr), &
                                             parts% npart_i(itr), &
                                             stretch, thres, &
                                             parts_all(itr)% pos_i, &
                                             parts_all(itr)% pvol_i, &
                                             import_density, &
                                             validate_position )

        ! Now that the real particle numbers are known, reallocate the arrays
        ! to the appropriate sizes. Note that, if the APM is performed,
        ! this step will be done after it as well
        parts_all(itr)% pos_i = &
                          parts_all(itr)% pos_i( :, 1:parts% npart_i(itr) )
        parts_all(itr)% pvol_i = &
                          parts_all(itr)% pvol_i( 1:parts% npart_i(itr) )

        IF(.NOT.ALLOCATED( parts_all(itr)% h_i ))THEN
          ALLOCATE( parts_all(itr)% h_i( parts% npart_i(itr) ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array h_i ", &
                     ". The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !        "...allocation error for array v_euler_z" )
        ENDIF
        DO a= 1, parts% npart_i(itr), 1
          parts_all(itr)% h_i(a)   = (parts_all(itr)% pvol_i(a))**third
        ENDDO

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
                                             parts_all(1)% pmass_i, &
                                             parts_all(1)% nu_i,    &
                                             parts_all(1)% h_i,     &
                                             parts% npart_i(1),     &
                                             parts_all(2)% pos_i,   &
                                             parts_all(2)% pvol_i,  &
                                             parts_all(2)% pmass_i, &
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


    CASE( id_particles_on_spherical_surfaces )


      PRINT *, "** Placing equal-mass particles on spherical surfaces, " &
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

        filename_mass_profile= &
          TRIM(sph_path)//"spherical_surfaces_mass_profile"//TRIM(str_i)//".dat"
        filename_shells_radii= &
          TRIM(sph_path)//"spherical_surfaces_radii"//TRIM(str_i)//".dat"
        filename_shells_pos  = &
          TRIM(sph_path)//"spherical_surfaces_pos"//TRIM(str_i)//".dat"

        CALL parts% place_particles_spherical_surfaces( &
                                              parts% masses(i_matter), &
                                              MAXVAL(sizes(i_matter,1:2)), &
                                              center(i_matter,1), &
                                              central_density(i_matter), &
                                              npart_des_i(i_matter), &
                                              parts% npart_i(i_matter), &
                                              parts_all(i_matter)% pos_i, &
                                              parts_all(i_matter)% pvol_i, &
                                              parts_all(i_matter)% pmass_i, &
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
                                              validate_position )

      !
      !-- Experimental code to get the desired number of particles for
      !-- matter objects which have an irregular geometry
      !-- Right now, the desired and real particle numbers are almost the same
      !-- only if the geometry of the system is elliptical
      !

      !  DO
      !
      !    ratio_npart_des_real= &
      !          ABS( DBLE(parts% npart_i(i_matter) - npart_des_i(i_matter)) ) &
      !          /DBLE(npart_des_i(i_matter))
      !
      !    IF( ratio_npart_des_real <= 0.15D0 )THEN
      !      EXIT
      !    ELSE
      !      pmass_des= SUM(parts_all(i_matter)% pmass_i) &
      !                 /SIZE(parts_all(i_matter)% pmass_i) &
      !                 *parts% npart_i(i_matter)/npart_des_i(i_matter)
      !    ENDIF
      !
      !    CALL parts% place_particles_spherical_surfaces( &
      !                                          parts% masses(i_matter), &
      !                                          MAXVAL(sizes(i_matter,1:2)), &
      !                                          center(i_matter,1), &
      !                                          central_density(i_matter), &
      !                                          npart_des_i(i_matter), &
      !                                          parts% npart_i(i_matter), &
      !                                          parts_all(i_matter)% pos_i, &
      !                                          parts_all(i_matter)% pvol_i, &
      !                                          parts_all(i_matter)% pmass_i, &
      !                                          last_r, &
      !                                          upper_bound, lower_bound, &
      !                                          upper_factor, lower_factor,&
      !                                          max_steps, &
      !                                          filename_mass_profile, &
      !                                          filename_shells_radii, &
      !                                          filename_shells_pos, &
      !                                          import_density, &
      !                                          integrate_mass_density, &
      !                                          import_id, &
      !                                          validate_position, &
      !                                          pmass_des )
      !
      !  ENDDO

        IF(.NOT.ALLOCATED( parts_all(i_matter)% h_i ))THEN
          ALLOCATE( parts_all(i_matter)% h_i( parts% npart_i(i_matter) ), &
                    STAT= ios, ERRMSG= err_msg )
          IF( ios > 0 )THEN
            PRINT *, "...allocation error for array h_i ", &
                     ". The error message is", err_msg
            STOP
          ENDIF
          !CALL test_status( ios, err_msg, &
          !        "...allocation error for array v_euler_z" )
        ENDIF
        DO a= 1, parts% npart_i(i_matter), 1
          parts_all(i_matter)% h_i(a)= (parts_all(i_matter)% pvol_i(a))**third
        ENDDO

        ! If there are 2 matter objects...
        equal_masses: IF( i_matter == 1 .AND. parts% n_matter == 2 )THEN

          ! ...with practically the same mass, and the physical system is
          ! symmetric wrt the yz plane (in which case the user should
          ! set reflect_particles_x in the parameter file)...
          IF( ABS(parts% masses(1) - parts% masses(2)) &
              /parts% masses(2) <= tol_equal_mass .AND. reflect_particles_x )THEN

            CALL reflect_particles_yz_plane( parts_all(1)% pos_i,   &
                                             parts_all(1)% pvol_i,  &
                                             parts_all(1)% pmass_i, &
                                             parts_all(1)% nu_i,    &
                                             parts_all(1)% h_i,     &
                                             parts% npart_i(1),     &
                                             parts_all(2)% pos_i,   &
                                             parts_all(2)% pvol_i,  &
                                             parts_all(2)% pmass_i, &
                                             parts_all(2)% nu_i,    &
                                             parts_all(2)% h_i,     &
                                             parts% npart_i(2) )

            EXIT

          ENDIF

        ENDIF equal_masses

        ! Now that the real particle numbers are known, reallocate the arrays
        ! with the appropriate sizes. If the APM is performed, this step is
        ! done after it too
        !PRINT *, "parts% npart_i(i_matter)=", parts% npart_i(i_matter)
        !STOP

      ENDDO matter_objects_sphersurfaces_loop
      CALL parts% placer_timer% stop_timer()

      DO i_matter= 1, parts% n_matter, 1

        parts_all(i_matter)% pos_i = &
                    parts_all(i_matter)% pos_i( :, 1:parts% npart_i(i_matter) )
        parts_all(i_matter)% pvol_i = &
                    parts_all(i_matter)% pvol_i( 1:parts% npart_i(i_matter) )
        parts_all(i_matter)% pmass_i = &
                    parts_all(i_matter)% pmass_i( 1:parts% npart_i(i_matter) )
        parts_all(i_matter)% nu_i = &
                    parts_all(i_matter)% pmass_i( 1:parts% npart_i(i_matter) )
      ENDDO

      parts% npart= SUM( parts% npart_i )

      PRINT *, " * Particles placed. Number of particles=", parts% npart
      DO itr= 1, parts% n_matter, 1
        PRINT *, " * Number of particles on object ", itr, "=", &
                 parts% npart_i(itr)
        PRINT *
      ENDDO
      PRINT *
      !STOP

    CASE DEFAULT

      PRINT *, "** There is no implemented particle placer " &
               // "corresponding to the number", dist
      PRINT *, " * Stopping..."
      PRINT *
      STOP

    END SELECT choose_particle_placer

    !--------------------------------------------------------------!
    !--  At this point,the particles are placed without the APM  --!
    !--------------------------------------------------------------!

    ! Reshape the arrays pos and pvol by deleting the unnecessary elements
  !  parts% pos = parts% pos( :, 1:parts% npart )
  !  parts% pvol= parts% pvol( 1:parts% npart )

    ! Check that there aren't particles with the same coordinates
    CALL parts% same_particle_timer% start_timer()
    check_particles_loop: DO i_matter= 1, parts% n_matter, 1
      CALL check_particle_positions( parts% npart_i(i_matter), &
                                     parts_all(i_matter)% pos_i )
    ENDDO check_particles_loop
    CALL parts% same_particle_timer% stop_timer()

    ! TODO: make the following an idbase type-bound procedure


    !
    !-- APM iteration
    !
    matter_objects_apm_loop: DO i_matter= 1, parts% n_matter, 1

      run_apm: IF( apm_iterate(i_matter) )THEN

        PRINT *
        PRINT *, "** Placing particles on matter object", i_matter, &
                 "using the APM..."
        PRINT *

        IF( i_matter <= 9 ) WRITE( str_i, '(I1)' ) i_matter
        IF( i_matter >= 10 .AND. parts% n_matter <= 99 ) &
                                              WRITE( str_i, '(I2)' ) i_matter
        IF( i_matter >= 100 .AND. parts% n_matter <= 999 ) &
                                              WRITE( str_i, '(I3)' ) i_matter

        filename_apm_pos_id = TRIM(sph_path)//"apm_pos_id"//TRIM(str_i)//".dat"
        filename_apm_pos    = TRIM(sph_path)//"apm_pos"//TRIM(str_i)//".dat"
        filename_apm_results= TRIM(sph_path)//"apm_results"//TRIM(str_i)//".dat"

        ! Matter object 1
        CALL parts% apm_timers(i_matter)% start_timer()
        CALL parts% perform_apm( &
                    import_density, get_nstar_id, &
                    parts% npart_i(i_matter), &
                    parts_all(i_matter)% pos_i, &
                    parts_all(i_matter)% pvol_i, &
                    parts_all(i_matter)% h_i, &
                    parts_all(i_matter)% nu_i, &
                    center(i_matter,:), barycenter(i_matter,:), &
                    parts% masses(i_matter), &
                    sizes(i_matter, :), &
                    apm_max_it, max_inc, mass_it, parts% correct_nu, &
                    nuratio_thres, nuratio_des, nx_gh, ny_gh, nz_gh, &
                    ghost_dist, &
                    use_atmosphere(i_matter), &
                    remove_atmosphere(i_matter), &
                    print_step, &
                    filename_apm_pos_id, filename_apm_pos, &
                    filename_apm_results, validate_position )
        CALL parts% apm_timers(i_matter)% stop_timer()

        parts_all(i_matter)% pmass_i = &
                    parts_all(i_matter)% nu_i( 1:parts% npart_i(i_matter) )*amu

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
            DEALLOCATE(parts_all(2)% pmass_i)
            DEALLOCATE(parts_all(2)% nu_i)

            CALL reflect_particles_yz_plane( parts_all(1)% pos_i,   &
                                             parts_all(1)% pvol_i,  &
                                             parts_all(1)% pmass_i, &
                                             parts_all(1)% nu_i,    &
                                             parts_all(1)% h_i,     &
                                             parts% npart_i(1),     &
                                             parts_all(2)% pos_i,   &
                                             parts_all(2)% pvol_i,  &
                                             parts_all(2)% pmass_i, &
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

    parts% npart= SUM( parts% npart_i )

    !
    !-- Reassign TYPE member variables after the APM iteration
    !

    IF( ALLOCATED(parts% h) ) DEALLOCATE( parts% h )
    ALLOCATE( parts% h( parts% npart ), &
              STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...allocation error for array pmass in SUBROUTINE" &
                // "place_particles_. ", &
                "The error message is", err_msg
       STOP
    ENDIF

    DO i_matter= 1, parts% n_matter, 1
      !IF( apm_iterate(i_matter) )THEN
        parts% h( parts% npart_i(i_matter-1) + 1: &
                  parts% npart_i(i_matter-1) + parts% npart_i(i_matter) )= &
                  parts_all(i_matter)% h_i(1:parts% npart_i(i_matter))
      !ENDIF
    ENDDO

    ALLOCATE( parts% pvol( parts% npart ), &
              STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...allocation error for array pvol in SUBROUTINE" &
                // "place_particles_. ", &
                "The error message is", err_msg
       STOP
    ENDIF

    DO i_matter= 1, parts% n_matter, 1
      parts% pvol( parts% npart_i(i_matter-1) + 1: &
                   parts% npart_i(i_matter-1) + parts% npart_i(i_matter) )= &
                   (parts_all(i_matter)% h_i(1:parts% npart_i(i_matter)))**3.0D0
    ENDDO

    IF( ALLOCATED(parts% nu) ) DEALLOCATE( parts% nu )
      ALLOCATE( parts% nu( parts% npart ), &
                STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...allocation error for array pmass in SUBROUTINE" &
                // "place_particles_. ", &
                "The error message is", err_msg
       STOP
    ENDIF

    DO i_matter= 1, parts% n_matter, 1

      ! If nu is known already at this stage, either because computed by the APM
      ! or because it is read from file, assign it to the TYPE member array
      IF( apm_iterate(i_matter) &
          .OR. &
          ( parts% distribution_id == id_particles_from_file .AND. read_nu ) &
          .OR. &
          ( parts% distribution_id == id_particles_on_spherical_surfaces ) &
      )THEN
        parts% nu( parts% npart_i(i_matter-1) + 1: &
                   parts% npart_i(i_matter-1) + parts% npart_i(i_matter) )= &
                   parts_all(i_matter)% nu_i
      ENDIF

    ENDDO

    !DEALLOCATE( parts% pos )
    ALLOCATE( parts% pos( 3, parts% npart ), &
              STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...allocation error for array pos in SUBROUTINE" &
                // "place_particles_. ", &
                "The error message is", err_msg
       STOP
    ENDIF

    DO i_matter= 1, parts% n_matter, 1
      parts% pos( :, parts% npart_i(i_matter-1) + 1: &
                     parts% npart_i(i_matter-1) + parts% npart_i(i_matter) )= &
                     parts_all(i_matter)% pos_i
    ENDDO

    ALLOCATE( parts% pmass( parts% npart ), &
              STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
       PRINT *, "...allocation error for array pmass in SUBROUTINE" &
                // "place_particles_. ", &
                "The error message is", err_msg
       STOP
    ENDIF
    parts% pmass= zero

    DO i_matter= 1, parts% n_matter, 1

      particles_not_on_lattice_or_apm: &
      IF( parts% distribution_id == id_particles_on_spherical_surfaces &
          .OR. ( parts% distribution_id == id_particles_from_file &
                 .AND. read_nu ) &
          .OR. apm_iterate(i_matter) )THEN
      ! If the particles are not placed on lattices, or the APM was performed...

            ! ...assign pmass to the TYPE member array
            parts% pmass( parts% npart_i(i_matter-1) + 1: &
                      parts% npart_i(i_matter-1) + parts% npart_i(i_matter) )= &
                      parts_all(i_matter)% pmass_i

      ENDIF particles_not_on_lattice_or_apm

    ENDDO

  !  !$OMP PARALLEL DO DEFAULT( NONE ) &
  !  !$OMP             SHARED( parts ) &
  !  !$OMP             PRIVATE( a, itr2 )
  !  find_nan_in_pos: DO a= 1, parts% npart, 1
  !
  !    DO itr2= 1, 3, 1
  !      IF( .NOT.is_finite_number(parts% pos(itr2,a)) )THEN
  !        PRINT *, "** ERROR! pos(", itr2, a, ")= ", parts% pos(itr2,a), &
  !                 " is not a finite number!"
  !        PRINT *, " * Stopping.."
  !        PRINT *
  !        STOP
  !      ENDIF
  !    ENDDO
  !
  !  ENDDO find_nan_in_pos
  !  !$OMP END PARALLEL DO

    !
    !-- Set the total center of mass of the system at the
    !-- Cartesian origin
    !
    ! TODO: The idbase object should tell the location of the total
    !       computing frame center of mass to the particle object

   ! DO i_matter= 1, parts% n_matter, 1
   !   parts% barycenter_system(1:3)= &
   !           parts% barycenter_system(1:3) + barycenter(i_matter,1:3)
   ! ENDDO
   ! PRINT *, parts% barycenter_system(1:3)
   ! STOP

    IF( dist /= id_particles_on_lattice )THEN

      CALL correct_center_of_mass( parts% npart, parts% pos, parts% nu, &
                                   import_density, &
                                   validate_position, [zero,zero,zero], &
                                   verbose= .FALSE. )

      CALL correct_center_of_mass( parts% npart, parts% pos, parts% nu, &
                                   import_density, &
                                   validate_position, [zero,zero,zero], &
                                   verbose= .FALSE. )

      CALL correct_center_of_mass( parts% npart, parts% pos, parts% nu, &
                                   import_density, &
                                   validate_position, [zero,zero,zero], &
                                   verbose= .TRUE. )

      CALL COM( & ! input
                parts% npart, parts% pos, parts% nu, &
                ! output
                parts% barycenter_system(1), parts% barycenter_system(2), &
                parts% barycenter_system(3), parts% barycenter_system(4) )

    ENDIF

    !PRINT *, com_x, com_y, com_z, com_d

    !STOP

    PRINT *, " * Final particle distribution prepared. Number of particles=", &
             parts% npart
    DO i_matter= 1, parts% n_matter, 1
      PRINT *, " * Number of particles on object ", i_matter, "=", &
               parts% npart_i(i_matter)
      PRINT *
    ENDDO
    PRINT *

    !-----------------------------------------------------------------------!
    !--  At this point, the particles are placed with or without the APM  --!
    !-----------------------------------------------------------------------!

    ! Allocate needed memory
    CALL allocate_particles_memory( parts )

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
                                parts% pos( 1, : ), &
                                parts% pos( 2, : ), &
                                parts% pos( 3, : ), &
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

    !PRINT *, MAXVAL( parts% baryon_density, DIM= 1 )
    !PRINT *, MAXVAL( parts% pressure, DIM= 1 )
    !STOP

    IF( debug ) PRINT *, "34"

    !-----------------------------------------------------------------------!
    ! If an atmosphere was used during the APM iteration, and kept, assign  !
    ! the minimum specific internal energy and the minimum velocity, to it. !
    ! N.B. The velocity has an hard-wired direction to reproduce counter-   !
    !      clockwise rotation.                                              !
    !-----------------------------------------------------------------------!

    matter_objects_atmo_loop: DO i_matter= 1, parts% n_matter, 1

    ASSOCIATE( npart_in   => parts% npart_i(i_matter-1) + 1, &
               npart_fin  => parts% npart_i(i_matter-1) +    &
                             parts% npart_i(i_matter) )

      IF( use_atmosphere(i_matter) .AND. .NOT.remove_atmosphere(i_matter) )THEN

        min_eps= MINVAL( parts% specific_energy(npart_in:npart_fin), &
                         DIM= 1, &
                MASK= parts% baryon_density(npart_in:npart_fin) > zero )
        min_vel= MINVAL( SQRT( &
                       (parts% v_euler_x(npart_in:npart_fin))**2.0D0 &
                     + (parts% v_euler_y(npart_in:npart_fin))**2.0D0 &
                     + (parts% v_euler_z(npart_in:npart_fin))**2.0D0 ), &
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
    !-- Check that the imported ID does not contain NaNs
    !
    !CALL Check_Array_for_NAN( parts% npart, parts% lapse, &
    !                                         "lapse" )
    !CALL Check_Array_for_NAN( parts% npart, parts% shift_x, &
    !                                         "shift_x" )
    !CALL Check_Array_for_NAN( parts% npart, parts% shift_y, &
    !                                         "shift_y" )
    !CALL Check_Array_for_NAN( parts% npart, parts% shift_z, &
    !                                         "shift_z" )
    !CALL Check_Array_for_NAN( parts% npart, parts% g_xx, &
    !                                         "g_xx" )
    !CALL Check_Array_for_NAN( parts% npart, parts% g_xy, &
    !                                         "g_xy" )
    !CALL Check_Array_for_NAN( parts% npart, parts% g_xz, &
    !                                         "g_xz" )
    !CALL Check_Array_for_NAN( parts% npart, parts% g_yy, &
    !                                         "g_yy" )
    !CALL Check_Array_for_NAN( parts% npart, parts% g_yz, &
    !                                         "g_yz" )
    !CALL Check_Array_for_NAN( parts% npart, parts% g_zz, &
    !                                         "g_zz" )
    !CALL Check_Array_for_NAN( parts% npart, &
    !        parts% baryon_density, "baryon_density" )
    !CALL Check_Array_for_NAN( parts% npart, &
    !        parts% energy_density, "energy_density" )
    !CALL Check_Array_for_NAN( parts% npart, &
    !        parts% specific_energy, "specific_energy" )
    !CALL Check_Array_for_NAN( parts% npart, &
    !               parts% pressure, "pressure" )
    !CALL Check_Array_for_NAN( parts% npart, &
    !              parts% v_euler_x, "v_euler_x" )
    !CALL Check_Array_for_NAN( parts% npart, &
    !              parts% v_euler_y, "v_euler_y" )
    !CALL Check_Array_for_NAN( parts% npart, &
    !              parts% v_euler_z, "v_euler_z" )

   ! IF(.NOT.ALLOCATED( parts% baryon_density_index ))THEN
   !   ALLOCATE( parts% baryon_density_index( parts% npart ), &
   !             STAT= ios, ERRMSG= err_msg )
   !   IF( ios > 0 )THEN
   !      PRINT *, "...allocation error for array baryon_density_index in " &
   !               // "SUBROUTINE construct_particles_idase. ", &
   !               "The error message is", err_msg
   !      STOP
   !   ENDIF
   !   !CALL test_status( ios, err_msg, &
   !   !                "...allocation error for array pos in SUBROUTINE" &
   !   !                // "place_particles_3D_lattice." )
   ! ENDIF

    !
    !-- Compute typical length-scale approximating g_00 with the Newtonian
    !-- potential
    !
    DO i_matter= 1, parts% n_matter, 1

      ASSOCIATE( npart_in   => parts% npart_i(i_matter-1) + 1, &
                 npart_fin  => parts% npart_i(i_matter-1) +    &
                               parts% npart_i(i_matter) )

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

    !PRINT *, "End of particle constructor"

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



    CONTAINS



    FUNCTION import_density( x, y, z ) RESULT( density )

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: x
      DOUBLE PRECISION, INTENT(IN):: y
      DOUBLE PRECISION, INTENT(IN):: z
      DOUBLE PRECISION:: density

      density= id% read_mass_density( x, y, z )

    END FUNCTION import_density


    SUBROUTINE import_id( x, y, z, &
                          sqdetg, &
                          baryon_density, &
                          gamma_euler )

      USE utility,  ONLY: determinant_sym3x3

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT( IN ) :: x
      DOUBLE PRECISION, INTENT( IN ) :: y
      DOUBLE PRECISION, INTENT( IN ) :: z
      DOUBLE PRECISION, INTENT( OUT ):: sqdetg
      DOUBLE PRECISION, INTENT( OUT ):: baryon_density
      DOUBLE PRECISION, INTENT( OUT ):: gamma_euler

      DOUBLE PRECISION, DIMENSION(6) :: g

      CALL id% read_id_mass_b( x, y, z, &
                               g, &
                               baryon_density, &
                               gamma_euler )

      CALL determinant_sym3x3(g,sqdetg)
      sqdetg= SQRT(sqdetg)

    END SUBROUTINE import_id


    SUBROUTINE integrate_mass_density( center, radius, &
                                       central_density, &
                                       dr, dth, dphi, &
                                       mass, mass_profile, &
                                       mass_profile_idx )

      IMPLICIT NONE

      !> Center of the star
      DOUBLE PRECISION, INTENT( IN )    :: center
      !> Central density of the star
      DOUBLE PRECISION, INTENT( IN )    :: central_density
      !> Radius of the star
      DOUBLE PRECISION, INTENT( IN )    :: radius
      !> Integration steps
      DOUBLE PRECISION, INTENT( IN )    :: dr, dth, dphi
      !> Integrated mass of the star
      DOUBLE PRECISION, INTENT( IN OUT ):: mass
      !> Array storing the radial mass profile of the star
      !DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT( INOUT ):: &
      !                                 mass_profile
      DOUBLE PRECISION, DIMENSION(3,0:NINT(radius/dr)), INTENT( OUT ):: &
                                           mass_profile
      !& Array to store the indices for array mass_profile, sorted so that
      !  mass_profile[mass_profile_idx] is in increasing order
      !INTEGER, DIMENSION(:), ALLOCATABLE, INTENT( INOUT ):: mass_profile_idx
      INTEGER, DIMENSION(0:NINT(radius/dr)), INTENT( OUT ):: mass_profile_idx

      CALL id% integrate_baryon_mass_density( center, radius, &
                                              central_density, &
                                              dr, dth, dphi, &
                                              mass, mass_profile, &
                                              mass_profile_idx )

    END SUBROUTINE integrate_mass_density


    FUNCTION validate_position( x, y, z ) RESULT( answer )

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: x
      DOUBLE PRECISION, INTENT(IN):: y
      DOUBLE PRECISION, INTENT(IN):: z
      LOGICAL:: answer

      answer= id% test_position( x, y, z )

    END FUNCTION validate_position


    SUBROUTINE correct_center_of_mass_of_system( npart, pos, nu, &
                                                 com_system )

      IMPLICIT NONE

      INTEGER, INTENT(IN):: npart
      DOUBLE PRECISION, INTENT(IN):: com_system(3)
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


    SUBROUTINE get_nstar_id( npart, x, y, z, nstar_sph, nstar_id, nlrf_sph,sqg )

      IMPLICIT NONE

      INTEGER, INTENT(IN):: npart
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
                                           specific_energy, &
                                           pressure, &
                                           v_euler_x, v_euler_y, v_euler_z

      CALL id% read_id_particles( npart, x, y, z, &
                                  lapse, shift_x, shift_y, shift_z, &
                                  g_xx, g_xy, g_xz, &
                                  g_yy, g_yz, g_zz, &
                                  baryon_density, &
                                  energy_density, &
                                  specific_energy, &
                                  pressure, &
                                  v_euler_x, v_euler_y, v_euler_z )

      CALL compute_nstar_id( npart, lapse, shift_x, shift_y, &
                             shift_z, v_euler_x, v_euler_y, v_euler_z, &
                             g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, &
                             baryon_density, nstar_sph, nstar_id, nlrf_sph, &
                             sqg )

      !CALL compute_nstar_eul_id( npart, &
      !                           v_euler_x, v_euler_y, v_euler_z, &
      !                           g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, &
      !                           baryon_density, nstar_eul_id )

    END SUBROUTINE get_nstar_id


    SUBROUTINE compute_nstar_id( npart, lapse, shift_x, shift_y, &
                                 shift_z, v_euler_x, v_euler_y, v_euler_z, &
                                 g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, &
                                 baryon_density, nstar_sph, nstar_id, nlrf_sph,&
                                 sqg )

      !**************************************************************
      !
      !# Compute nstar_id, the proper baryon mass density, given the
      !  |id|
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

        gamma_eul_a= one/SQRT(one-v_euler_norm2)
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


    SUBROUTINE reflect_particles_yz_plane( pos_star1, pvol_star1, pmass_star1, &
                                           nu_star1, h_star1, npart_star1,     &
                                           pos_star2, pvol_star2, pmass_star2, &
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
      DOUBLE PRECISION, DIMENSION(:),   INTENT(IN):: pmass_star1
      !! Array where to store the particle masses for star 1
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
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, INTENT(INOUT):: pmass_star2
      !! Array where to store the particle masses for star 2
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, INTENT(INOUT):: nu_star2
      !! Array where to store the particle baryon number for star 2
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, INTENT(INOUT):: h_star2
      !! Array where to store the particle smoothing lengths for star 2
      INTEGER,                                       INTENT(INOUT):: npart_star2
      !! Variable where to store the particle number for star 2


      IF(.NOT.ALLOCATED( pos_star2 ))THEN
        ALLOCATE( pos_star2( 3, npart_star1 ), STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array pos_star2 in SUBROUTINE" &
                    // "reflect_particles_yz_plane. ", &
                    "The error message is", err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( pvol_star2 ))THEN
        ALLOCATE( pvol_star2( npart_star1 ), STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array pvol_star2 in SUBROUTINE" &
                    // "reflect_particles_yz_plane. ", &
                    "The error message is", err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( pmass_star2 ))THEN
        ALLOCATE( pmass_star2( npart_star1 ), STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array pmass_star2 in SUBROUTINE" &
                    // "reflect_particles_yz_plane. ", &
                    "The error message is", err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( nu_star2 ))THEN
        ALLOCATE( nu_star2( npart_star1 ), STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array nu_star2 in SUBROUTINE" &
                    // "reflect_particles_yz_plane. ", &
                    "The error message is", err_msg
           STOP
        ENDIF
      ENDIF
      IF(.NOT.ALLOCATED( h_star2 ))THEN
        ALLOCATE( h_star2( npart_star1 ), STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array h_star2 in SUBROUTINE" &
                    // "reflect_particles_yz_plane. ", &
                    "The error message is", err_msg
           STOP
        ENDIF
      ENDIF

      PRINT *, " * Reflecting particles about the yz plane..."
      PRINT *

      ! Reflect the particles on matter object 1, and their properties,
      ! to matter object 2
      pos_star2(1,:)= - pos_star1(1,:) !+ ( ABS(center(2,:)) - ABS(center(1,:)) )
      pos_star2(2,:)=   pos_star1(2,:)
      pos_star2(3,:)=   pos_star1(3,:)
      pvol_star2    =   pvol_star1
      pmass_star2   =   pmass_star1
      nu_star2      =   nu_star1
      h_star2       =   h_star1
      npart_star2   =   npart_star1


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
