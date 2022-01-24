! File:         submodule_sph_particles_constructor.f90
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

SUBMODULE (sph_particles) constructor

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


  MODULE PROCEDURE construct_particles

    !**************************************************
    !
    !# The constructor performs all the tasks needed
    !  to set up the particle distribution with the
    !  |id| evaluated on it. It calls all the PROCEDURES
    !  that rely on an object of STATIC TYPE idbase.
    !
    !  @todo assign sub-tasks to separate SUBROUTINES
    !        CONTAINED in this SUBMODULE
    !
    !  FT 17.10.2020
    !
    !  @note Last updated: FT 17.01.2022
    !
    !**************************************************

    !USE NaNChecker,    ONLY: Check_Array_for_NAN
    USE constants,      ONLY: amu, pi, zero, half, third
    USE NR,             ONLY: indexx
    USE kernel_table,   ONLY: ktable
    USE input_output,   ONLY: read_options
    USE units,          ONLY: set_units
    USE options,        ONLY: ikernel, ndes
    USE alive_flag,     ONLY: alive
    USE utility,        ONLY: spherical_from_cartesian

    IMPLICIT NONE

    ! The variable counter counts how many times the PROCEDURE
    ! construct_particles_idase is called
    INTEGER, SAVE:: counter= 1
    INTEGER:: npart_des, a, max_steps, nlines, header_lines, n_cols, &
              npart_tmp, nx_gh, ny_gh, nz_gh, i_matter

    ! Maximum length for strings, and for the number of imported binaries
    INTEGER, PARAMETER:: max_length= 50
    ! APM parameters
    INTEGER:: apm_max_it, max_inc
    INTEGER, PARAMETER:: unit_pos= 2289
    ! Variable storing the number of column where nu is written
    INTEGER:: column_nu
    ! Temporary number of matter objects
    INTEGER:: n_matter_tmp, tmp
    ! Array storing the columns of the file parts_pos (defined below) that
    ! contain the particle positions
    INTEGER, DIMENSION(3):: columns
    INTEGER, DIMENSION(id% get_n_matter()):: npart_des_i
    ! Temporary array storing the number of particles on each matter object
    INTEGER, DIMENSION(:), ALLOCATABLE:: npart_i_tmp

    DOUBLE PRECISION:: thres, nu_ratio_des
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

    DOUBLE PRECISION, DIMENSION( :, : ), ALLOCATABLE:: tmp_pos
    DOUBLE PRECISION:: nuratio_thres, nuratio_des

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

    NAMELIST /bns_particles/ &
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
              use_atmosphere, remove_atmosphere, nuratio_des

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
    !-- Read needed data from the idbase object
    !

    parts% nbar_tot       = 0.0D0
    parts% npart          = 0
    parts% distribution_id= dist

    ALLOCATE( parts% masses (parts% n_matter) )
    ALLOCATE( parts% all_eos(parts% n_matter) )
    ALLOCATE( parts% npart_i(0:parts% n_matter) )
    ALLOCATE( npart_i_tmp(0:parts% n_matter) )
    ALLOCATE( parts% nbar_i(parts% n_matter) )
    ALLOCATE( parts% nuratio_i(parts% n_matter) )
    ALLOCATE( parts% mass_ratios (parts% n_matter) )
    ALLOCATE( parts% mass_fractions (parts% n_matter) )
    parts% npart_i(0)= 0
    npart_i_tmp(0)   = 0
    parts% nbar_i    = 0.0D0
    parts% nuratio_i = 0.0D0

    DO i_matter= 1, parts% n_matter, 1

      parts% masses(i_matter)  = id% return_mass(i_matter)
      center(i_matter,:)       = id% return_center(i_matter)
      central_density(i_matter)= id% read_mass_density( center(i_matter,1), &
                                                        center(i_matter,2), &
                                                        center(i_matter,3) )
      barycenter(i_matter,:)   = id% return_barycenter(i_matter)
      sizes(i_matter, :)       = id% return_spatial_extent(i_matter)

      parts% all_eos(i_matter)% eos_name= id% return_eos_name(i_matter)
      CALL id% return_eos_parameters( i_matter, &
                                      parts% all_eos(i_matter)% eos_parameters )

    ENDDO

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

    READ( 10, NML= bns_particles )
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

    ! setup unit system
    CALL set_units('NSM')
    CALL read_options

    ! tabulate kernel, get ndes
    CALL ktable( ikernel, ndes )

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

           ! PRINT *, ABS( MINVAL( tmp_pos(1,npart_in:npart_fin) ) ) > &
           !         ABS(center(itr,1)) + sizes(itr, 1)
           ! PRINT *, ABS( MAXVAL( tmp_pos(1,npart_in:npart_fin) ) ) > &
           !         ABS(center(itr,1)) + sizes(itr, 2)
           ! PRINT *, ABS( MINVAL( tmp_pos(2,npart_in:npart_fin) ) ) > &
           !         ABS(center(itr,2)) + sizes(itr, 3)
           ! PRINT *, ABS( MAXVAL( tmp_pos(2,npart_in:npart_fin) ) ) > &
           !         ABS(center(itr,2)) + sizes(itr, 4)
           ! PRINT *, ABS( MINVAL( tmp_pos(3,npart_in:npart_fin) ) ) > &
           !         ABS(center(itr,3)) + sizes(itr, 5)
           ! PRINT *, ABS( MAXVAL( tmp_pos(3,npart_in:npart_fin) ) ) > &
           !         ABS(center(itr,3)) + sizes(itr, 6)
           !
           ! PRINT *, ABS( MINVAL( tmp_pos(1,npart_in:npart_fin) ) )
           ! PRINT *, ABS(center(itr,1)) + sizes(itr, 1)
           ! PRINT *, ABS( MAXVAL( tmp_pos(1,npart_in:npart_fin) ) )
           ! PRINT *, ABS(center(itr,1)) + sizes(itr, 2)
           ! PRINT *, ABS( MINVAL( tmp_pos(2,npart_in:npart_fin) ) )
           ! PRINT *, ABS(center(itr,2)) + sizes(itr, 3)
           ! PRINT *, ABS( MAXVAL( tmp_pos(2,npart_in:npart_fin) ) )
           ! PRINT *, ABS(center(itr,2)) + sizes(itr, 4)
           ! PRINT *, ABS( MINVAL( tmp_pos(3,npart_in:npart_fin) ) )
           ! PRINT *, ABS(center(itr,3)) + sizes(itr, 5)
           ! PRINT *, ABS( MAXVAL( tmp_pos(3,npart_in:npart_fin) ) )
           ! PRINT *, ABS(center(itr,3)) + sizes(itr, 6)

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

      ! Allocating the memory for the array pos( 3, npart )
  !    IF(.NOT.ALLOCATED( parts% pos ))THEN
  !      ALLOCATE( parts% pos( 3, parts% npart ), STAT= ios, &
  !                ERRMSG= err_msg )
  !      IF( ios > 0 )THEN
  !         PRINT *, "...allocation error for array pos in SUBROUTINE" &
  !                  // ". ", &
  !                  "The error message is", err_msg
  !         STOP
  !      ENDIF
  !      !CALL test_status( ios, err_msg, &
  !      !                "...allocation error for array pos in SUBROUTINE" &
  !      !                // "place_particles_3D_lattice." )
  !    ENDIF
  !    IF( read_nu .AND. .NOT.ALLOCATED( parts% nu ))THEN
  !      ALLOCATE( parts% nu( parts% npart ), STAT= ios, &
  !                ERRMSG= err_msg )
  !      IF( ios > 0 )THEN
  !         PRINT *, "...allocation error for array nu in SUBROUTINE" &
  !                  // ". ", &
  !                  "The error message is", err_msg
  !         STOP
  !      ENDIF
  !      !CALL test_status( ios, err_msg, &
  !      !                "...allocation error for array pos in SUBROUTINE" &
  !      !                // "place_particles_3D_lattice." )
  !    ENDIF
  !    IF( read_nu .AND. .NOT.ALLOCATED( parts% pmass ))THEN
  !      ALLOCATE( parts% pmass( parts% npart ), STAT= ios, &
  !                ERRMSG= err_msg )
  !      IF( ios > 0 )THEN
  !         PRINT *, "...allocation error for array pmass in SUBROUTINE" &
  !                  // ". ", &
  !                  "The error message is", err_msg
  !         STOP
  !      ENDIF
  !      !CALL test_status( ios, err_msg, &
  !      !                "...allocation error for array pos in SUBROUTINE" &
  !      !                // "place_particles_3D_lattice." )
  !    ENDIF

      !parts% pos= tmp_pos(1:3,1:parts% npart)
      !IF( read_nu ) parts% nu= tmp_pos(4,1:parts% npart)

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
   !     !        "...allocation error for array v_euler_parts_z" )
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
            !        "...allocation error for array v_euler_parts_z" )
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
            !        "...allocation error for array v_euler_parts_z" )
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

      ! Assign values to the arrays storing the data for each matter object
      ! @todo: this is repetition. You should use directly the arrays associated
      !        to each matter object in the previous code, not the final arrays
      !        parts% pos, parts% pvol, parts% pmass and parts% nu
   !   DO itr= 1, parts% n_matter, 1
   !
   !     ASSOCIATE( npart_in   => npart_i_tmp(itr)*(itr-1) + 1, &
   !                npart_fin  => npart_i_tmp(itr) + npart_i_tmp(itr)*(itr-1) )
   !
   !       parts_all(itr)% pos_i = parts% pos( :, npart_in:npart_fin )
   !       parts_all(itr)% pvol_i= parts% pvol( npart_in:npart_fin )
   !       parts_all(itr)% h_i   = (parts_all(itr)% pvol_i)**third
   !
   !       IF( read_nu )THEN
   !
   !         parts_all(itr)% nu_i   = parts% nu( npart_in:npart_fin )
   !         parts_all(itr)% pmass_i= parts% pmass( npart_in:npart_fin )
   !
   !         PRINT *, " * Maximum n. baryon per particle (nu) on object", itr, &
   !                             "=", MAXVAL( parts_all(itr)% nu_i, DIM= 1 )
   !         PRINT *, " * Minimum n. baryon per particle (nu) on object", itr, &
   !                             "=", MINVAL( parts_all(itr)% nu_i, DIM= 1 )
   !         PRINT *, " * Ratio between the two=", &
   !                  MAXVAL( parts_all(itr)% nu_i, DIM= 1 )&
   !                 /MINVAL( parts_all(itr)% nu_i, DIM= 1 )
   !         PRINT *
   !
   !       ENDIF
   !
   !     END ASSOCIATE
   !
   !   ENDDO


    CASE( id_particles_on_lattice )


      PRINT *, " * Placing particles on lattices, one around each ", &
               "matter object."
      PRINT *

      ! Place particles, and time the proces

      CALL parts% placer_timer% start_timer()
      DO itr= 1, parts% n_matter, 1

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
                                             check_negative_hydro )

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
          !        "...allocation error for array v_euler_parts_z" )
        ENDIF
        DO a= 1, parts% npart_i(itr), 1
          parts_all(itr)% h_i(a)   = (parts_all(itr)% pvol_i(a))**third
        ENDDO

      ENDDO
      CALL parts% placer_timer% stop_timer()

      parts% npart= SUM( parts% npart_i )

      IF( debug ) PRINT *, "10"

      !
      !-- Assign TYPE member variables
      !-- @todo: this shouldn´t be needed here, since these array should be
      !--        allocated and populated after the APM code
      !

   !   IF(.NOT.ALLOCATED( parts% pos ))THEN
   !     ALLOCATE( parts% pos( 3, parts% npart ), &
   !               STAT= ios, ERRMSG= err_msg )
   !     IF( ios > 0 )THEN
   !        PRINT *, "...allocation error for array pos in SUBROUTINE" &
   !                 // "place_particles_. ", &
   !                 "The error message is", err_msg
   !        STOP
   !     ENDIF
   !     !CALL test_status( ios, err_msg, &
   !     !                "...allocation error for array pos in SUBROUTINE" &
   !     !                // "place_particles_3D_lattice." )
   !   ENDIF
   !   !parts% pos( :, 1:parts% npart1 )= pos1
   !   !parts% pos( :, parts% npart1 + 1:parts% npart )= pos2
   !   DO itr= 1, parts% n_matter, 1
   !     parts% pos( :, parts% npart_i(itr-1) + 1: &
   !                    parts% npart_i(itr-1) + parts% npart_i(itr) )= &
   !                    parts_all(itr)% pos_i
   !   ENDDO
   !
   !   IF(.NOT.ALLOCATED( parts% pvol ))THEN
   !     ALLOCATE( parts% pvol( parts% npart ), &
   !               STAT= ios, ERRMSG= err_msg )
   !     IF( ios > 0 )THEN
   !        PRINT *, "...allocation error for array pvol in SUBROUTINE" &
   !                 // "place_particles_. ", &
   !                 "The error message is", err_msg
   !        STOP
   !     ENDIF
   !     !CALL test_status( ios, err_msg, &
   !     !                "...allocation error for array pos in SUBROUTINE" &
   !     !                // "place_particles_3D_lattice." )
   !   ENDIF
   !   !parts% pvol( 1:parts% npart1 )= pvol1
   !   !parts% pvol( parts% npart1 + 1:parts% npart )= pvol2
   !   DO itr= 1, parts% n_matter, 1
   !     parts% pvol( parts% npart_i(itr-1) + 1: &
   !                  parts% npart_i(itr-1) + parts% npart_i(itr) )= &
   !                  parts_all(itr)% pvol_i
   !   ENDDO


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

        filename_mass_profile= "spherical_surfaces_mass_profile"//TRIM(str_i)//&
                               ".dat"
        filename_shells_radii= "spherical_surfaces_radii"//TRIM(str_i)//".dat"
        filename_shells_pos  = "spherical_surfaces_pos"//TRIM(str_i)//".dat"

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
                                              check_negative_hydro )

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
      !                                          check_negative_hydro, &
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
          !        "...allocation error for array v_euler_parts_z" )
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
              /parts% masses(2) <= 0.005 .AND. reflect_particles_x )THEN

            IF(.NOT.ALLOCATED( parts_all(2)% pos_i ))THEN
              ALLOCATE( parts_all(2)% pos_i( 3, parts% npart_i(1) ), &
                        STAT= ios, ERRMSG= err_msg )
              IF( ios > 0 )THEN
                 PRINT *, "...allocation error for array pos in SUBROUTINE" &
                          // "place_particles_. ", &
                          "The error message is", err_msg
                 STOP
              ENDIF
            ENDIF
            IF(.NOT.ALLOCATED( parts_all(2)% pvol_i ))THEN
              ALLOCATE( parts_all(2)% pvol_i( parts% npart_i(1) ), &
                        STAT= ios, ERRMSG= err_msg )
              IF( ios > 0 )THEN
                 PRINT *, "...allocation error for array pvol in SUBROUTINE" &
                          // "place_particles_. ", &
                          "The error message is", err_msg
                 STOP
              ENDIF
            ENDIF
            IF(.NOT.ALLOCATED( parts_all(2)% pmass_i ))THEN
              ALLOCATE( parts_all(2)% pmass_i( parts% npart_i(1) ), &
                        STAT= ios, ERRMSG= err_msg )
              IF( ios > 0 )THEN
                 PRINT *, "...allocation error for array pmass in SUBROUTINE" &
                          // "place_particles_. ", &
                          "The error message is", err_msg
                 STOP
              ENDIF
            ENDIF

            ! Reflect the particles on matter object 1, and their properties,
            ! to matter object 2
            parts_all(2)% pos_i(1,:)= - parts_all(1)% pos_i(1,:)
            parts_all(2)% pos_i(2,:)=   parts_all(1)% pos_i(2,:)
            parts_all(2)% pos_i(3,:)=   parts_all(1)% pos_i(3,:)
            parts_all(2)% pvol_i    =   parts_all(1)% pvol_i
            parts_all(2)% h_i       =   parts_all(1)% h_i
            parts_all(2)% pmass_i   =   parts_all(1)% pmass_i
            parts% npart_i(2)       =   parts% npart_i(1)

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

    !  ASSOCIATE( npart_in   => parts% npart_i(i_matter-1) + 1, &
    !             npart_fin  => parts% npart_i(i_matter-1) +    &
    !                           parts% npart_i(i_matter) )

      IF( apm_iterate(i_matter) )THEN

      !  IF(.NOT.ALLOCATED( parts% h ))THEN
      !    ALLOCATE( parts% h( parts% npart ), STAT= ios, &
      !              ERRMSG= err_msg )
      !    IF( ios > 0 )THEN
      !       PRINT *, "...allocation error for array h in SUBROUTINE ", &
      !                "construct_particles_idase. The error message is",&
      !                err_msg
      !       STOP
      !    ENDIF
      !  ENDIF
      !
      !  IF(.NOT.ALLOCATED( parts% nu ))THEN
      !    ALLOCATE( parts% nu( parts% npart ), STAT= ios, &
      !              ERRMSG= err_msg )
      !    IF( ios > 0 )THEN
      !       PRINT *, "...allocation error for array nu in SUBROUTINE ", &
      !                "construct_particles_idase. The error message is",&
      !                err_msg
      !       STOP
      !    ENDIF
      !  ENDIF

        PRINT *
        PRINT *, "** Placing particles on matter object", i_matter, &
                 "using the APM..."
        PRINT *

        IF( i_matter <= 9 ) WRITE( str_i, '(I1)' ) i_matter
        IF( i_matter >= 10 .AND. parts% n_matter <= 99 ) &
                                              WRITE( str_i, '(I2)' ) i_matter
        IF( i_matter >= 100 .AND. parts% n_matter <= 999 ) &
                                              WRITE( str_i, '(I3)' ) i_matter

        filename_apm_pos_id = "apm_pos_id"//TRIM(str_i)//".dat"
        filename_apm_pos    = "apm_pos"//TRIM(str_i)//".dat"
        filename_apm_results= "apm_results"//TRIM(str_i)//".dat"

        ! Matter object 1
        CALL parts% apm_timers(i_matter)% start_timer()
        CALL parts% perform_apm( &
                    import_density, get_nstar_p, &
                    parts% npart_i(i_matter), &
                    !npart_i_tmp(i_matter), &
                    !parts% pos(:,npart_in:npart_fin), &
                    !parts% pvol(npart_in:npart_fin), &
                    !parts% h(npart_in:npart_fin), &
                    !parts% nu(npart_in:npart_fin), &
                    parts_all(i_matter)% pos_i, &
                    parts_all(i_matter)% pvol_i, &
                    parts_all(i_matter)% h_i, &
                    parts_all(i_matter)% nu_i, &
                    center(i_matter,:), barycenter(i_matter,:), &
                    parts% masses(i_matter), &
                    sizes(i_matter, :), &
                    apm_max_it, max_inc, mass_it, parts% correct_nu, &
                    nuratio_thres, nuratio_des, nx_gh, ny_gh, nz_gh, &
                    use_atmosphere(i_matter), &
                    remove_atmosphere(i_matter), &
                    filename_apm_pos_id, filename_apm_pos, &
                    filename_apm_results, check_negative_hydro )
        CALL parts% apm_timers(i_matter)% stop_timer()

        parts_all(i_matter)% pmass_i = &
                    parts_all(i_matter)% nu_i( 1:parts% npart_i(i_matter) )*amu

        IF( debug ) PRINT *, "average nu= ", &
          SUM(parts_all(i_matter)% nu_i, DIM= 1)/SIZE(parts_all(i_matter)% nu_i)

        PRINT *, "** Particles placed on matter object", i_matter, &
                 "according to the APM."
        PRINT *

        !DEALLOCATE

      !  parts_all(i_matter)% pos_i = &
      !            parts% pos(:,npart_in:npart_in + parts% npart_i_tmp(i_matter))
      !  parts_all(i_matter)% h_i = &
      !            parts% h(npart_in:npart_in + parts% npart_i_tmp(i_matter))
      !  parts_all(i_matter)% nu_i = &
      !            parts% nu(npart_in:npart_in + parts% npart_i_tmp(i_matter))

        ! If there are 2 matter objects...
        equal_masses_apm: IF( i_matter == 1 .AND. parts% n_matter == 2 )THEN

          ! with practically the same mass, and the physical system
          ! is symmetric wrt the yz plane (in which case the user should set
          ! the reflect_particles_x to .TRUE. in the parameter file)
          IF( ABS(parts% mass_ratios(1) - parts% mass_ratios(2) ) &
                          /parts% mass_ratios(2) <= 0.005 .AND. &
                          reflect_particles_x )THEN

            !parts% pos(1,parts% npart_i(1)+1:parts% npart)= &
            !                          - parts% pos(1,1:parts% npart_i(1))
            !parts% pos(2,parts% npart_i(1)+1:parts% npart)= &
            !                            parts% pos(2,1:parts% npart_i(1))
            !parts% pos(3,parts% npart_i(1)+1:parts% npart)= &
            !                            parts% pos(3,1:parts% npart_i(1))
            !
            !parts% nu(parts% npart_i(1)+1:parts% npart)= &
            !                            parts% nu(1:parts% npart_i(1))
            !
            !parts% h(parts% npart_i(1)+1:parts% npart)= &
            !                            parts% h(1:parts% npart_i(1))

       !     parts% pos(1,parts% npart_i(1)+1:parts% npart(1)+parts% npart_i_tmp(1))= &
       !                               - parts% pos(1,1:parts% npart_i_tmp(1))
       !     parts% pos(2,parts% npart_i(1)+1:parts% npart(1)+parts% npart_i_tmp(1))= &
       !                                 parts% pos(2,1:parts% npart_i_tmp(1))
       !     parts% pos(3,parts% npart_i(1)+1:parts% npart(1)+parts% npart_i_tmp(1))= &
       !                                 parts% pos(3,1:parts% npart_i_tmp(1))
       !
       !     parts% nu(parts% npart_i(1)+1:parts% npart(1)+parts% npart_i_tmp(1))= &
       !                                 parts% nu(1:parts% npart_i_tmp(1))
       !
       !     parts% h(parts% npart_i(1)+1:parts% npart(1)+parts% npart_i_tmp(1))= &
       !                                 parts% h(1:parts% npart_i_tmp(1))

            !parts% npart_i(2)= parts% npart_i(1)
            !parts% npart= parts% npart_i(1) + parts% npart_i(1)

        !    parts% npart_i_tmp(2)= parts% npart_i_tmp(1)
            ! ...reflect particles
            parts_all(2)% pos_i(1,:)= - parts_all(1)% pos_i(1,:)
            parts_all(2)% pos_i(2,:)=   parts_all(1)% pos_i(2,:)
            parts_all(2)% pos_i(3,:)=   parts_all(1)% pos_i(3,:)
            parts_all(2)% pvol_i    =   parts_all(1)% pvol_i
            parts_all(2)% h_i       =   parts_all(1)% h_i
            parts_all(2)% pmass_i   =   parts_all(1)% pmass_i
            parts_all(2)% nu_i      =   parts_all(1)% nu_i
            parts% npart_i(2)       =   parts% npart_i(1)

            PRINT *, "** Particles placed on star 1 according to the APM", &
                     " reflected about the yz plane onto star 2."
            PRINT *

            EXIT

          ENDIF

        ENDIF equal_masses_apm

      ENDIF

   !   END ASSOCIATE

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
                      parts_all(i_matter)% h_i
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
                   (parts_all(i_matter)% h_i)**3.0D0
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
          .OR. ( parts% distribution_id == 0 .AND. read_nu ) )THEN
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

    PRINT *, " * Particles placed according to the APM. Number of particles=", &
             parts% npart
    DO i_matter= 1, parts% n_matter, 1
      PRINT *, " * Number of particles on object ", i_matter, "=", &
               parts% npart_i(i_matter)
      PRINT *
    ENDDO
    PRINT *

    !STOP

    !-----------------------------------------------------------------------!
    !--  At this point, the particles are placed with or without the APM  --!
    !-----------------------------------------------------------------------!

    ! Allocate needed memory
    CALL allocate_lorene_id_parts_memory( parts )

    ! flag that particles are 'alive'
    ALLOCATE( alive( parts% npart ) )
    alive( 1:parts% npart )= 1

    IF( debug ) PRINT *, "33"

    !
    !-- Import the needed ID on the particles, and time the process
    !
    PRINT *, "** Assigning the ID to the particles..."
    PRINT *

    CALL parts% importer_timer% start_timer()
    CALL id% read_id_particles( parts% npart, &
                                parts% pos( 1, : ), &
                                parts% pos( 2, : ), &
                                parts% pos( 3, : ), &
                                parts% lapse_parts, &
                                parts% shift_parts_x, &
                                parts% shift_parts_y, &
                                parts% shift_parts_z, &
                                parts% g_xx_parts, &
                                parts% g_xy_parts, &
                                parts% g_xz_parts, &
                                parts% g_yy_parts, &
                                parts% g_yz_parts, &
                                parts% g_zz_parts, &
                                parts% baryon_density_parts, &
                                parts% energy_density_parts, &
                                parts% specific_energy_parts, &
                                parts% pressure_parts, &
                                parts% v_euler_parts_x, &
                                parts% v_euler_parts_y, &
                                parts% v_euler_parts_z )
    CALL parts% importer_timer% stop_timer()

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

        min_eps= MINVAL( parts% specific_energy_parts(npart_in:npart_fin), &
                         DIM= 1, &
                MASK= parts% baryon_density_parts(npart_in:npart_fin) > 0.0D0 )
        min_vel= MINVAL( SQRT( &
                       (parts% v_euler_parts_x(npart_in:npart_fin))**2.0D0 &
                     + (parts% v_euler_parts_y(npart_in:npart_fin))**2.0D0 &
                     + (parts% v_euler_parts_z(npart_in:npart_fin))**2.0D0 ), &
                         DIM= 1, &
                MASK= parts% baryon_density_parts(npart_in:npart_fin) > 0.0D0 )

        particle_loop2: DO a= npart_in, npart_fin, 1

          IF( parts% baryon_density_parts(a) <= 0.0D0 )THEN

          !  IF( parts% pos(1,a) > 0.0D0 )THEN
          !
          !    phi_a= ATAN( &
          !                ( parts% pos(2,a) - center(i_matter,2) ) &
          !               /( parts% pos(1,a) - center(i_matter,1) ) &
          !              )
          !
          !  ELSEIF( parts% pos(1,a) < 0.0D0 )THEN
          !
          !    phi_a= ATAN( &
          !                ( parts% pos(2,a) - center(i_matter,2) ) &
          !               /( parts% pos(1,a) - center(i_matter,1) ) &
          !              ) + pi
          !
          !  ELSE
          !
          !    phi_a= pi/2.0D0
          !
          !  ENDIF
          !
          !  theta_a= ACOS( &
          !              ( parts% pos(3,a) - center(i_matter,3) ) &
          !              /SQRT( &
          !                ( parts% pos(1,a) - center(i_matter,1) )**2.0D0 &
          !              + ( parts% pos(2,a) - center(i_matter,2) )**2.0D0 &
          !              + ( parts% pos(3,a) - center(i_matter,3) )**2.0D0 &
          !              ) &
          !            )

            CALL spherical_from_cartesian( &
                  parts% pos(1,a), parts% pos(2,a), parts% pos(3,a), &
                  center(i_matter,1), center(i_matter,2), center(i_matter,3), &
                  r_a, theta_a, phi_a )

            parts% specific_energy_parts(a)= min_eps

            parts% v_euler_parts_x(a)      = &
      ( min_vel*SIN(theta_a - pi*half)*COS(phi_a) + parts% shift_parts_x(a) ) &
              /parts% lapse_parts(a)
            parts% v_euler_parts_y(a)      = &
      ( min_vel*SIN(theta_a - pi*half)*SIN(phi_a) + parts% shift_parts_y(a) ) &
              /parts% lapse_parts(a)
            parts% v_euler_parts_z(a)      = &
              ( min_vel*COS(theta_a - pi*half) + parts% shift_parts_z(a) ) &
              /parts% lapse_parts(a)

          ENDIF

        ENDDO particle_loop2

      ENDIF

    END ASSOCIATE

    ENDDO matter_objects_atmo_loop

    !
    !-- Check that the imported ID does not contain NaNs
    !
    !CALL Check_Array_for_NAN( parts% npart, parts% lapse_parts, &
    !                                         "lapse_parts" )
    !CALL Check_Array_for_NAN( parts% npart, parts% shift_parts_x, &
    !                                         "shift_parts_x" )
    !CALL Check_Array_for_NAN( parts% npart, parts% shift_parts_y, &
    !                                         "shift_parts_y" )
    !CALL Check_Array_for_NAN( parts% npart, parts% shift_parts_z, &
    !                                         "shift_parts_z" )
    !CALL Check_Array_for_NAN( parts% npart, parts% g_xx_parts, &
    !                                         "g_xx_parts" )
    !CALL Check_Array_for_NAN( parts% npart, parts% g_xy_parts, &
    !                                         "g_xy_parts" )
    !CALL Check_Array_for_NAN( parts% npart, parts% g_xz_parts, &
    !                                         "g_xz_parts" )
    !CALL Check_Array_for_NAN( parts% npart, parts% g_yy_parts, &
    !                                         "g_yy_parts" )
    !CALL Check_Array_for_NAN( parts% npart, parts% g_yz_parts, &
    !                                         "g_yz_parts" )
    !CALL Check_Array_for_NAN( parts% npart, parts% g_zz_parts, &
    !                                         "g_zz_parts" )
    !CALL Check_Array_for_NAN( parts% npart, &
    !        parts% baryon_density_parts, "baryon_density_parts" )
    !CALL Check_Array_for_NAN( parts% npart, &
    !        parts% energy_density_parts, "energy_density_parts" )
    !CALL Check_Array_for_NAN( parts% npart, &
    !        parts% specific_energy_parts, "specific_energy_parts" )
    !CALL Check_Array_for_NAN( parts% npart, &
    !               parts% pressure_parts, "pressure_parts" )
    !CALL Check_Array_for_NAN( parts% npart, &
    !              parts% v_euler_parts_x, "v_euler_parts_x" )
    !CALL Check_Array_for_NAN( parts% npart, &
    !              parts% v_euler_parts_y, "v_euler_parts_y" )
    !CALL Check_Array_for_NAN( parts% npart, &
    !              parts% v_euler_parts_z, "v_euler_parts_z" )

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

    !PRINT *, "baryon_density_index"
    !DO itr= 1, parts% npart1, 1
    !  PRINT *, parts% baryon_density_index( itr ), &
    !           parts% baryon_density_parts( itr )
    !ENDDO
    !    PRINT *, "baryon_density_parts in ascending order"
    !DO itr= 1, parts% npart1, 1
    !  PRINT *, parts% baryon_density_parts( &
    !                                parts% baryon_density_index( itr ) )
    !ENDDO
    !PRINT *, "baryon_density_parts in descending order"
    !DO itr= parts% npart1, 1, -1
    !  PRINT *, parts% baryon_density_parts( &
    !                                parts% baryon_density_index( itr ) )
    !ENDDO
    ! Ok it seems working

!    PRINT *, "** Computing typical length scale for the change in pressure", &
!             " on the x axis."
!    PRINT *
!
!    ALLOCATE( abs_pos( 3, parts% npart ) )
!
!    DO itr = 1, parts% npart, 1
!      abs_pos( 1, itr )= ABS( parts% pos( 1, itr ) )
!      abs_pos( 2, itr )= ABS( parts% pos( 2, itr ) )
!      abs_pos( 3, itr )= ABS( parts% pos( 3, itr ) )
!    ENDDO
!
!    min_y_index= 0
!    min_abs_y= 1D+20
!    DO itr = 1, parts% npart, 1
!      IF( ABS( parts% pos( 2, itr ) ) < min_abs_y )THEN
!        min_abs_y= ABS( parts% pos( 2, itr ) )
!        min_y_index= itr
!      ENDIF
!    ENDDO
!
!    min_z_index= 0
!    min_abs_z= 1D+20
!    DO itr = 1, parts% npart, 1
!      IF( ABS( parts% pos( 3, itr ) ) < min_abs_z )THEN
!        min_abs_z= ABS( parts% pos( 3, itr ) )
!        min_z_index= itr
!      ENDIF
!    ENDDO
!
!    min_abs_z= MINVAL( abs_pos( 3, : ) )
!
!    !PRINT *, "1"
!
!    cntr1= 0
!    cntr2= 0
!    DO itr = 1, parts% npart, 1
!      IF( parts% pos( 3, itr ) == min_abs_z &
!          .AND. &
!          ABS( ( parts% pos( 2, itr ) - &
!                 parts% pos( 2, min_y_index ) )/ &
!                 parts% pos( 2, min_y_index ) ) < 1.0D-5 &
!      )THEN
!
!        IF( parts% pos( 1, itr ) < 0 )THEN
!          cntr1= cntr1 + 1
!        ELSEIF( parts% pos( 1, itr ) > 0 )THEN
!          cntr2= cntr2 + 1
!        ENDIF
!
!      ENDIF
!    ENDDO
!    !PRINT *, "cntr1= ", cntr1
!    !PRINT *, "cntr2= ", cntr2
!
!    ALLOCATE( parts% pos_x1( cntr1 ) )
!    ALLOCATE( parts% pos_x2( cntr2 ) )
!    ALLOCATE( parts% pressure_parts_x1( cntr1 ) )
!    ALLOCATE( parts% pressure_parts_x2( cntr2 ) )
!    ALLOCATE( parts% pressure_parts_x_der1( cntr1 - 5 ) )
!    ALLOCATE( parts% pressure_parts_x_der2( cntr2 - 5 ) )
!    ALLOCATE( parts% pressure_length_scale_x1( cntr1 - 5 ) )
!    ALLOCATE( parts% pressure_length_scale_x2( cntr2 - 5 ) )
!
!    !PRINT *, "2"
!
!    itr_1= 0
!    itr_2= 0
!    DO itr = 1, parts% npart, 1
!      IF( parts% pos( 3, itr ) == min_abs_z &
!          .AND. &
!          ABS( ( parts% pos( 2, itr ) - &
!                 parts% pos( 2, min_y_index ) )/ &
!                 parts% pos( 2, min_y_index ) ) < 1.0D-5 &
!        )THEN
!
!        IF( parts% pos( 1, itr ) < 0 )THEN
!          itr_1= itr_1 + 1
!          parts% pos_x1( itr_1 )= parts% pos( 1, itr )
!          parts% pressure_parts_x1( itr_1 )= &
!                                              parts% pressure_parts( itr )
!        ELSEIF( parts% pos( 1, itr ) > 0 )THEN
!          itr_2= itr_2 + 1
!          parts% pos_x2( itr_2 )= parts% pos( 1, itr )
!          parts% pressure_parts_x2( itr_2 )= &
!                                              parts% pressure_parts( itr )
!        ENDIF
!
!      ENDIF
!    ENDDO
!
!    !PRINT *, "3"
!
!    DO itr= 3, cntr1 - 3, 1
!      parts% pressure_parts_x_der1( itr - 2 )=&
!                     ( + parts% pressure_parts_x1( itr - 2 )/12.0D0 &
!                       - 2.0*parts% pressure_parts_x1( itr - 1 )/3.0D0 &
!                       + 2.0*parts% pressure_parts_x1( itr + 1 )/3.0D0 &
!                       - parts% pressure_parts_x1( itr + 2 )/12.0D0 )&
!                       /( Msun_geo*km2m*ABS( parts% pos_x1( itr ) - &
!                                             parts% pos_x1( itr - 1 ) ) )
!
!      parts% pressure_length_scale_x1( itr - 2 )= &
!                          ABS( parts% pressure_parts_x1( itr - 2 )/ &
!                               parts% pressure_parts_x_der1( itr - 2 ) )
!
!      !PRINT *, "p1=", parts% pressure_parts_x1( itr - 2 )
!      !PRINT *, "p_r1=", parts% pressure_parts_x_der1( itr - 2 )
!      !PRINT *, "p/p_r1=", parts% pressure_length_scale_x1( itr - 2 )
!      !PRINT *
!
!    ENDDO
!    DO itr= 3, cntr2 - 3, 1
!      parts% pressure_parts_x_der2( itr - 2 )=&
!                     ( + parts% pressure_parts_x2( itr - 2 )/12.0D0 &
!                       - 2.0*parts% pressure_parts_x2( itr - 1 )/3.0D0 &
!                       + 2.0*parts% pressure_parts_x2( itr + 1 )/3.0D0 &
!                       - parts% pressure_parts_x2( itr + 2 )/12.0D0 )&
!                       /( Msun_geo*km2m*ABS( parts% pos_x2( itr ) - &
!                                             parts% pos_x2( itr - 1 ) ) )
!
!      parts% pressure_length_scale_x2( itr - 2 )= &
!                          ABS( parts% pressure_parts_x2( itr - 2 )/ &
!                               parts% pressure_parts_x_der2( itr - 2 ) )
!
!      !PRINT *, "p2=", parts% pressure_parts_x2( itr - 2 )
!      !PRINT *, "p_r2=", parts% pressure_parts_x_der2( itr - 2 )
!      !PRINT *, "p/p_r2=", parts% pressure_length_scale_x2( itr - 2 )
!      !PRINT *
!
!    ENDDO
!
!    PRINT *, " * Maximum typical length scale for change in pressure", &
!             " along the x axis for NS 1= ", &
!             MAXVAL( parts% pressure_length_scale_x1, DIM= 1 )/km2m, " km"
!    PRINT *, " * Minimum typical length scale for change in pressure", &
!             " along the x axis for NS 1= ", &
!             MINVAL( parts% pressure_length_scale_x1, DIM= 1 )/km2m, " km"
!    PRINT *
!    PRINT *, " * Maximum typical length scale for change in pressure", &
!             " along the x axis for NS 2= ", &
!             MAXVAL( parts% pressure_length_scale_x2, DIM= 1 )/km2m, " km"
!    PRINT *, " * Minimum typical length scale for change in pressure", &
!             " along the x axis for NS 2= ", &
!             MINVAL( parts% pressure_length_scale_x2, DIM= 1 )/km2m, " km"
!    PRINT *

    ! Increase the counter that identifies the particle distribution
    counter= counter + 1

    !PRINT *, "End of particle constructor"

  !  IF( parts% redistribute_nu )THEN
  !
  !    ! Index particles on star 1 in increasing order of nu
  !
  !    CALL indexx( parts% npart1, &
  !                 parts% baryon_density_parts( 1 : parts% npart1 ), &
  !                 parts% baryon_density_index( 1 : parts% npart1 ) )
  !
  !    ! Index particles on star 2 in increasing order of nu
  !
  !    CALL indexx( parts% npart2, &
  !                 parts% baryon_density_parts( parts% npart1 + 1 : &
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

    ! TODO: fix this by removing the abs_pos array
 !   IF( debug )THEN
 !
 !     namefile= "dbg-hydro.dat"
 !
 !     INQUIRE( FILE= TRIM(namefile), EXIST= exist )
 !
 !     IF( exist )THEN
 !         OPEN( UNIT= 2, FILE= TRIM(namefile), STATUS= "REPLACE", &
 !               FORM= "FORMATTED", &
 !               POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
 !               IOMSG= err_msg )
 !     ELSE
 !         OPEN( UNIT= 2, FILE= TRIM(namefile), STATUS= "NEW", &
 !               FORM= "FORMATTED", &
 !               ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
 !     ENDIF
 !     IF( ios > 0 )THEN
 !       PRINT *, "...error when opening " // TRIM(namefile), &
 !                ". The error message is", err_msg
 !       STOP
 !     ENDIF
 !     !CALL test_status( ios, err_msg, "...error when opening " &
 !     !                  // TRIM(namefile) )
 !
 !     WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
 !     "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id
 !
 !     WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
 !     "# Values of the fields (including coordinates) exported by |lorene| "&
 !     // "on each grid point"
 !     IF( ios > 0 )THEN
 !       PRINT *, "...error when writing line 1 in " // TRIM(namefile), &
 !                ". The error message is", err_msg
 !       STOP
 !     ENDIF
 !     !CALL test_status( ios, err_msg, "...error when writing line 1 in "&
 !     !        // TRIM(namefile) )
 !
 !     WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
 !     "# column:      1        2       3       4       5", &
 !     "       6       7       8", &
 !     "       9       10      11", &
 !     "       12      13      14", &
 !     "       15      16      17      18"
 !
 !     IF( ios > 0 )THEN
 !       PRINT *, "...error when writing line 2 in " // TRIM(namefile), &
 !                ". The error message is", err_msg
 !       STOP
 !     ENDIF
 !     !CALL test_status( ios, err_msg, "...error when writing line 2 in "&
 !     !            // TRIM(namefile) )
 !
 !     WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
 !     "#      grid point      x [km]       y [km]       z [km]       lapse", &
 !     "       shift_x [c]    shift_y [c]    shift_z [c]", &
 !     "       baryon density in the local rest frame [kg m^{-3}$]", &
 !     "       energy density [c^2]", &
 !     "       specific energy [c^2]", &
 !     "       pressure [Pa]", &
 !     "       fluid 3-velocity wrt the Eulerian observer (3 columns) [c]", &
 !     "       fluid coordinate 3-velocity vel_u (3 columns) [c]", &
 !     "       baryon number per particle nu", &
 !     "       baryon density in the local rest frame nlrf [baryon/Msun_geo^3]", &
 !     "       electron fraction", &
 !     "       generalized Lorentz factor Theta"
 !     IF( ios > 0 )THEN
 !       PRINT *, "...error when writing line 3 in " // TRIM(namefile), &
 !                ". The error message is", err_msg
 !       STOP
 !     ENDIF
 !     !CALL test_status( ios, err_msg, "...error when writing line 3 in "&
 !     !          // TRIM(namefile) )
 !
 !     DO itr = 1, parts% npart, 1
 !       abs_pos( 1, itr )= ABS( parts% pos( 1, itr ) )
 !       abs_pos( 2, itr )= ABS( parts% pos( 2, itr ) )
 !       abs_pos( 3, itr )= ABS( parts% pos( 3, itr ) )
 !     ENDDO
 !
 !     min_y_index= 0
 !     min_abs_y= 1D+20
 !     DO itr = 1, parts% npart, 1
 !       IF( ABS( parts% pos( 2, itr ) ) < min_abs_y )THEN
 !         min_abs_y= ABS( parts% pos( 2, itr ) )
 !         min_y_index= itr
 !       ENDIF
 !     ENDDO
 !
 !     min_abs_z= MINVAL( abs_pos( 3, : ) )
 !
 !     write_data_loop: DO itr = 1, parts% npart, 1
 !
 !       IF( parts% export_form_xy .AND. &
 !           parts% pos( 3, itr ) /= min_abs_z )THEN
 !         CYCLE
 !       ENDIF
 !       IF( parts% export_form_x .AND. &
 !           ( parts% pos( 3, itr ) /= min_abs_z &
 !             .OR. &
 !             parts% pos( 2, itr ) /= parts% pos( 2, min_y_index ) ) &
 !       )THEN
 !         CYCLE
 !       ENDIF
 !       WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
 !         itr, &
 !         parts% pos( 1, itr ), &
 !         parts% pos( 2, itr ), &
 !         parts% pos( 3, itr ), &
 !         parts% lapse_parts( itr ), &
 !         parts% shift_parts_x( itr ), &
 !         parts% shift_parts_y( itr ), &
 !         parts% shift_parts_z( itr ), &
 !         parts% baryon_density_parts( itr ), &
 !         parts% energy_density_parts( itr ), &
 !         parts% specific_energy_parts( itr ), &
 !         parts% pressure_parts( itr ), &
 !         parts% v_euler_parts_x( itr ), &
 !         parts% v_euler_parts_y( itr ), &
 !         parts% v_euler_parts_z( itr )
 !
 !     IF( ios > 0 )THEN
 !       PRINT *, "...error when writing the arrays in " // TRIM(namefile), &
 !                ". The error message is", err_msg
 !       STOP
 !     ENDIF
 !     !CALL test_status( ios, err_msg, "...error when writing " &
 !     !         // "the arrays in " // TRIM(finalnamefile) )
 !     ENDDO write_data_loop
 !
 !     CLOSE( UNIT= 2 )
 !
 !   ENDIF



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

      USE matrix, ONLY: determinant_3x3_sym_matrix

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

      CALL determinant_3x3_sym_matrix(g,sqdetg)
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


    FUNCTION check_negative_hydro( x, y, z ) RESULT( answer )

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: x
      DOUBLE PRECISION, INTENT(IN):: y
      DOUBLE PRECISION, INTENT(IN):: z
      INTEGER:: answer

      answer= id% test_position( x, y, z )

    END FUNCTION check_negative_hydro


    SUBROUTINE get_nstar_p( npart_real, x, y, z, nstar_p )

      IMPLICIT NONE

      INTEGER, INTENT(IN):: npart_real
      DOUBLE PRECISION, INTENT(IN):: x(npart_real)
      DOUBLE PRECISION, INTENT(IN):: y(npart_real)
      DOUBLE PRECISION, INTENT(IN):: z(npart_real)
      DOUBLE PRECISION, INTENT(OUT):: nstar_p(npart_real)

      DOUBLE PRECISION, DIMENSION(npart_real):: lapse, &
                                                shift_x, shift_y, shift_z, &
                                                g_xx, g_xy, g_xz, &
                                                g_yy, g_yz, g_zz, &
                                                baryon_density, &
                                                energy_density, &
                                                specific_energy, &
                                                pressure, &
                                                v_euler_x, v_euler_y, v_euler_z

      CALL id% read_id_particles( npart_real, x, y, z, &
                                  lapse, shift_x, shift_y, shift_z, &
                                  g_xx, g_xy, g_xz, &
                                  g_yy, g_yz, g_zz, &
                                  baryon_density, &
                                  energy_density, &
                                  specific_energy, &
                                  pressure, &
                                  v_euler_x, v_euler_y, v_euler_z )

      CALL compute_nstar_p( npart_real, lapse, shift_x, shift_y, &
                            shift_z, v_euler_x, v_euler_y, v_euler_z, &
                            g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, &
                            baryon_density, nstar_p )

    END SUBROUTINE get_nstar_p


    SUBROUTINE compute_nstar_p( npart_real, lapse, shift_x, shift_y, &
                                shift_z, v_euler_x, v_euler_y, v_euler_z, &
                                g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, &
                                baryon_density, nstar_p )

      !**************************************************************
      !
      !# Compute nstar_p, the proper baryon mass density, given the
      !  |lorene| ID
      !
      !  FT 31.08.2021
      !
      !**************************************************************

      !USE constants, ONLY: MSun, Msun_geo, km2m, amu, g2kg
      USE matrix,    ONLY: determinant_4x4_matrix

      IMPLICIT NONE

      INTEGER, INTENT(IN):: npart_real
      DOUBLE PRECISION, DIMENSION(npart_real), INTENT(IN):: lapse
      DOUBLE PRECISION, DIMENSION(npart_real), INTENT(IN):: shift_x
      DOUBLE PRECISION, DIMENSION(npart_real), INTENT(IN):: shift_y
      DOUBLE PRECISION, DIMENSION(npart_real), INTENT(IN):: shift_z
      DOUBLE PRECISION, DIMENSION(npart_real), INTENT(IN):: v_euler_x
      DOUBLE PRECISION, DIMENSION(npart_real), INTENT(IN):: v_euler_y
      DOUBLE PRECISION, DIMENSION(npart_real), INTENT(IN):: v_euler_z
      DOUBLE PRECISION, DIMENSION(npart_real), INTENT(IN):: g_xx
      DOUBLE PRECISION, DIMENSION(npart_real), INTENT(IN):: g_xy
      DOUBLE PRECISION, DIMENSION(npart_real), INTENT(IN):: g_xz
      DOUBLE PRECISION, DIMENSION(npart_real), INTENT(IN):: g_yy
      DOUBLE PRECISION, DIMENSION(npart_real), INTENT(IN):: g_yz
      DOUBLE PRECISION, DIMENSION(npart_real), INTENT(IN):: g_zz
      DOUBLE PRECISION, DIMENSION(npart_real), INTENT(IN):: baryon_density
      DOUBLE PRECISION, DIMENSION(npart_real), INTENT(OUT):: nstar_p

      INTEGER:: a, mus, nus
      DOUBLE PRECISION:: det, sq_g, Theta_a
      DOUBLE PRECISION, DIMENSION(0:3,npart_real):: vel
      DOUBLE PRECISION:: g4(0:3,0:3)

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( npart_real, lapse, shift_x, shift_y, shift_z, &
      !$OMP                     v_euler_x, v_euler_y, v_euler_z, &
      !$OMP                     g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, &
      !$OMP                     baryon_density, vel, nstar_p ) &
      !$OMP             PRIVATE( a, det, sq_g, Theta_a, g4 )
      DO a= 1, npart_real, 1

        ! Coordinate velocity of the fluid [c]
        vel(0,a)= 1.0D0
        vel(1,a)= lapse(a)*v_euler_x(a)- shift_x(a)
        vel(2,a)= lapse(a)*v_euler_y(a)- shift_y(a)
        vel(3,a)= lapse(a)*v_euler_z(a)- shift_z(a)

        !
        !-- Metric as matrix for easy manipulation
        !
        g4(0,0)= - lapse(a)**2 + g_xx(a)*shift_x(a)*shift_x(a)&
               + 2*g_xy(a)*shift_x(a)*shift_y(a) &
               + 2*g_xz(a)*shift_x(a)*shift_z(a) &
               + g_yy(a)*shift_y(a)*shift_y(a) &
               + 2*g_yz(a)*shift_y(a)*shift_z(a) &
               + g_zz(a)*shift_z(a)*shift_z(a)
        g4(0,1)= g_xx(a)*shift_x(a) + g_xy(a)*shift_y(a) + g_xz(a)*shift_z(a)
        g4(0,2)= g_xy(a)*shift_x(a) + g_yy(a)*shift_y(a) + g_yz(a)*shift_z(a)
        g4(0,3)= g_xz(a)*shift_x(a) + g_yz(a)*shift_y(a) + g_zz(a)*shift_z(a)

        g4(1,0)= g_xx(a)*shift_x(a) + g_xy(a)*shift_y(a) + g_xz(a)*shift_z(a)
        g4(1,1)= g_xx(a)
        g4(1,2)= g_xy(a)
        g4(1,3)= g_xz(a)

        g4(2,0)= g_xy(a)*shift_x(a) + g_yy(a)*shift_y(a) + g_yz(a)*shift_z(a)
        g4(2,1)= g_xy(a)
        g4(2,2)= g_yy(a)
        g4(2,3)= g_yz(a)

        g4(3,0)= g_xz(a)*shift_x(a) + g_yz(a)*shift_y(a) + g_zz(a)*shift_z(a)
        g4(3,1)= g_xz(a)
        g4(3,2)= g_yz(a)
        g4(3,3)= g_zz(a)

        ! sqrt(-det(g4))
        CALL determinant_4x4_matrix(g4,det)
        IF( ABS(det) < 1D-10 )THEN
            PRINT *, "The determinant of the spacetime metric is " &
                     // "effectively 0 at particle ", a
            STOP
        ELSEIF( det > 0 )THEN
            PRINT *, "The determinant of the spacetime metric is " &
                     // "positive at particle ", a
            STOP
        ENDIF
        sq_g= SQRT(-det)

        !
        !-- Generalized Lorentz factor
        !
        Theta_a= 0.D0
        DO nus=0,3
          DO mus=0,3
            Theta_a= Theta_a &
                     + g4(mus,nus)*vel(mus,a)*vel(nus,a)
          ENDDO
        ENDDO
        Theta_a= 1.0D0/SQRT(-Theta_a)
        !Theta(a)= Theta_a

        nstar_p(a)= sq_g*Theta_a*baryon_density(a)
                    !*((Msun_geo*km2m)**3)/(amu*g2kg)

      ENDDO
      !$OMP END PARALLEL DO

    END SUBROUTINE compute_nstar_p


  END PROCEDURE construct_particles


  MODULE PROCEDURE destruct_particles

    !*********************************************
    !
    !# Destructor of a particles object
    !
    !  FT
    !
    !*********************************************

    IMPLICIT NONE

    CALL THIS% deallocate_lorene_id_parts_memory()


  END PROCEDURE destruct_particles


END SUBMODULE constructor
