! File:         submodule_particles_constructor.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (particles_id) particles_constructor

  !************************************************
  !
  !# This SUBMODULE contains the implementation
  !  of the constructor and the
  !  destructor of TYPE particles.
  !
  !  FT 16.10.2020
  !
  !************************************************


  IMPLICIT NONE


  CONTAINS


  !MODULE PROCEDURE construct_particles_empty
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
  !    parts_obj% empty_object= .TRUE.
  !
  !    parts_obj% npart_temp= 0
  !
  !END PROCEDURE construct_particles_empty


  MODULE PROCEDURE construct_particles

    !**************************************************
    !
    !# The constructor performs all the tasks needed
    !  to set up the particle distribution with the
    !  LORENE ID on it. It calls all the PROCEDURES
    !  that rely on an object of TYPE bns.
    !
    !  @todo assign sub-tasks to separate SUBROUTINES
    !        CONTAINED in this SUBMODULE
    !
    !  FT 17.10.2020
    !
    !**************************************************

    !USE NaNChecker, ONLY: Check_Array_for_NAN
    USE constants,      ONLY: Msun_geo, km2m, amu
    USE NR,             ONLY: indexx
    USE kernel_table,   ONLY: ktable
    USE input_output,   ONLY: read_options
    USE units,          ONLY: set_units
    USE options,        ONLY: ikernel, ndes
    USE alive_flag,     ONLY: alive

    IMPLICIT NONE

    ! The variable counter counts how many times the PROCEDURE
    ! construct_particles is called
    INTEGER, SAVE:: counter= 1
    INTEGER:: nx, ny, nz, &
              npart_approx, npart2_approx, max_steps, &
              nlines, header_lines, n_cols, npart_tmp, npart1_tmp, npart2_tmp, &
              nx_gh, ny_gh, nz_gh
    ! Maximum length for strings, and for the number of imported binaries
    INTEGER, PARAMETER:: max_length= 50
    ! APM parameters
    INTEGER:: apm_max_it, max_inc, n_particles_first_shell
    INTEGER, PARAMETER:: unit_pos= 2289
    ! Variable storing the number of column where nu is written
    INTEGER:: column_nu
    ! Array storing the columns of the file parts_pos (defined below) that
    ! contain the particle positions
    INTEGER, DIMENSION(3):: columns

    DOUBLE PRECISION:: thres, nu_ratio
    DOUBLE PRECISION:: xmin, xmax, ymin, ymax, zmin, zmax, stretch
    DOUBLE PRECISION:: xmin1, xmax1, ymin1, ymax1, zmin1, zmax1
    DOUBLE PRECISION:: xmin2, xmax2, ymin2, ymax2, zmin2, zmax2
    DOUBLE PRECISION:: center1, center2, radius1, radius2, com1, com2
    DOUBLE PRECISION:: upper_bound, lower_bound, upper_factor, lower_factor, &
                       last_r
    DOUBLE PRECISION:: pvol_tmp

    DOUBLE PRECISION, DIMENSION( :, : ), ALLOCATABLE:: tmp_pos
    DOUBLE PRECISION, DIMENSION( :, : ), ALLOCATABLE:: tmp_pos2
    DOUBLE PRECISION, DIMENSION( :, : ), ALLOCATABLE:: pos1, pos2
    DOUBLE PRECISION, DIMENSION( : ),    ALLOCATABLE:: pvol1, pvol2, &
                                                       pmass1, pmass2
    DOUBLE PRECISION:: nuratio_thres, nuratio_des

    ! String storing the name of the directory storing the files containing
    ! the particle distributions
    CHARACTER( LEN= max_length ):: parts_pos_path
    ! String storing the name of the file containing the particle positions
    CHARACTER( LEN= max_length ):: parts_pos
    ! Final name for the file containing the particle positions
    CHARACTER( LEN= : ), ALLOCATABLE:: parts_pos_namefile
    ! String storing the local path to the directory where the
    ! LORENE BNS ID files are stored
    CHARACTER( LEN= max_length ):: compose_path
    ! String storing the names of the LORENE BNS ID binary files
    CHARACTER( LEN= max_length ):: compose_filename

    CHARACTER( LEN= max_length ):: filename_apm_pos_id, filename_apm_pos, &
                                   filename_apm_results

    CHARACTER( LEN= max_length ):: filename_mass_profile, &
                                   filename_shells_radii, filename_shells_pos

    LOGICAL:: file_exists, use_thres, redistribute_nu, correct_nu, &
              compose_eos, exist, randomize_phi, randomize_theta, &
              randomize_r, apm_iterate1, apm_iterate2, mass_it, &
              read_nu, reflect_particles_x

    LOGICAL, PARAMETER:: debug= .FALSE.

    NAMELIST /bns_particles/ &
              parts_pos_path, parts_pos, columns, header_lines, n_cols, &
              read_nu, column_nu, &
              stretch, &
              nx, ny, nz, &
              use_thres, thres, nu_ratio, redistribute_nu, correct_nu, &
              compose_eos, compose_path, compose_filename, &
              npart_approx, last_r, upper_bound, lower_bound, &
              upper_factor, lower_factor, max_steps, &
              randomize_phi, randomize_theta, randomize_r, &
              apm_iterate1, apm_iterate2, apm_max_it, max_inc, mass_it, &
              nuratio_thres, reflect_particles_x, nx_gh, ny_gh, nz_gh, &
              nuratio_des

    !
    !-- Initialize the timers
    !
    parts_obj% placer_timer       = timer( "placer_timer" )
    parts_obj% apm1_timer         = timer( "apm_star1_timer" )
    parts_obj% apm2_timer         = timer( "apm_star2_timer" )
    parts_obj% importer_timer     = timer( "importer_timer" )
    parts_obj% sph_computer_timer = timer( "sph_computer_timer" )
    parts_obj% same_particle_timer= timer( "same_particle_timer" )

    ! Declare this object as non-empty (experimental)
    parts_obj% empty_object= .FALSE.

    parts_obj% mass1          = bns_obj% get_mass1()
    parts_obj% mass2          = bns_obj% get_mass2()
    center1                   = bns_obj% get_center1_x()
    center2                   = bns_obj% get_center2_x()
    com1                      = bns_obj% get_barycenter1_x()
    com2                      = bns_obj% get_barycenter2_x()
    radius1                   = bns_obj% get_radius1_x_comp()
    radius2                   = bns_obj% get_radius2_x_comp()
    parts_obj% nbar_tot       = 0.0D0
    parts_obj% nbar1          = 0.0D0
    parts_obj% nbar2          = 0.0D0
    parts_obj% npart          = 0.0D0
    parts_obj% distribution_id= dist

    parts_obj% eos1= bns_obj% get_eos1()
    parts_obj% eos2= bns_obj% get_eos2()

    parts_obj% eos1_id= bns_obj% get_eos1_id()
    parts_obj% eos2_id= bns_obj% get_eos2_id()

    parts_obj% gamma_sp1= bns_obj% get_gamma_1()
    parts_obj% kappa_sp1= bns_obj% get_kappa_1()
    parts_obj% gamma_sp2= bns_obj% get_gamma_2()
    parts_obj% kappa_sp2= bns_obj% get_kappa_2()

    !
    !-- Read the parameters of the particle distributions
    !
    parts_obj% lorene_bns_id_parfile= 'sphincs_lorene_bns_particles.par'

    INQUIRE( FILE= parts_obj% lorene_bns_id_parfile, EXIST= file_exists )
    IF( file_exists )THEN
     OPEN( 10, FILE= parts_obj% lorene_bns_id_parfile, STATUS= 'OLD' )
    ELSE
     PRINT *
     PRINT *, "** ERROR: ", parts_obj% lorene_bns_id_parfile, &
              " file not found!"
     PRINT *
     STOP
    ENDIF

    READ( 10, NML= bns_particles )
    CLOSE( 10 )

    parts_obj% use_thres          = use_thres
    parts_obj% correct_nu         = correct_nu
    parts_obj% compose_eos        = compose_eos
    parts_obj% compose_path       = compose_path
    parts_obj% compose_filename   = compose_filename
    parts_obj% redistribute_nu    = redistribute_nu
    parts_obj% nu_ratio           = nu_ratio
    parts_obj% reflect_particles_x= reflect_particles_x
    parts_obj% randomize_phi      = randomize_phi
    parts_obj% randomize_theta    = randomize_theta
    parts_obj% randomize_r        = randomize_r
    ! APM parameters
    parts_obj% apm_iterate1   = apm_iterate1
    parts_obj% apm_iterate2   = apm_iterate2
    !parts_obj% apm_max_it   = apm_max_it
    !parts_obj% max_inc      = max_inc
    !parts_obj% mass_it      = mass_it
    !parts_obj% nuratio_thres= nuratio_thres
    parts_obj% read_nu       = read_nu

    parts_pos_namefile= TRIM(parts_pos_path)//TRIM(parts_pos)

    IF( parts_obj% redistribute_nu )THEN
      thres= 100.0D0*parts_obj% nu_ratio
    ENDIF

    IF( MOD( nz, 2 ) /= 0 )THEN
      PRINT *
      PRINT *, "** ERROR in lorene_bns_id_particles.par: nz should be even!"
      PRINT *
      STOP
    ENDIF

    IF( nx == 0 .OR. ny == 0 .OR. nz == 0 )THEN
      PRINT *
      PRINT *, "** ERROR in lorene_bns_id_particles.par: ", &
               "nx, ny, nz cannot be 0!"
      PRINT *
      STOP
    ENDIF

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
    CALL ktable(ikernel,ndes)

    ! TODO: Add check that the number of rows in placer is the same as the
    !       number of bns objects, and that all bns have a value for placer

    !
    !-- Choose particle placer
    !
    choose_particle_placer: SELECT CASE( dist )

    CASE(0)

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

      nlines = 0
      DO
        READ( unit_pos, * , IOSTAT= ios )
        IF ( ios /= 0 ) EXIT
        nlines = nlines + 1
      ENDDO

      CLOSE( UNIT= unit_pos )

      npart_tmp= nlines - header_lines

      !PRINT *, "nlines=", nlines
      !PRINT *, "header_lines=", header_lines
      !PRINT *, "npart_tmp=", npart_tmp
      !PRINT *

      OPEN( UNIT= unit_pos, FILE= TRIM(parts_pos_namefile), &
            FORM= "FORMATTED", ACTION= "READ" )

      ! Skip header
      DO itr= 1, header_lines, 1
        READ( unit_pos, * )
      ENDDO

      ! Allocate the temporary array to store data
      ALLOCATE( tmp_pos( n_cols, 2*npart_tmp ) )
      tmp_pos= 0.0D0

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

      ! Allocate the temporary array to store data
      IF( read_nu )THEN
        ALLOCATE( tmp_pos2( 4, 2*npart_tmp ) )
      ELSE
        ALLOCATE( tmp_pos2( 3, 2*npart_tmp ) )
      ENDIF
      tmp_pos2= 0.0D0

      ! Separate particle positions on star 1 and star 2,
      ! and compute the temporary npart1 and npart2 (before mirroring)
      npart1_tmp= 0
      DO itr= 1, npart_tmp, 1

        IF( tmp_pos(columns(1),itr) < 0.0D0 )THEN

          npart1_tmp= npart1_tmp + 1
          tmp_pos2(1,npart1_tmp)= tmp_pos(columns(1),itr)
          tmp_pos2(2,npart1_tmp)= tmp_pos(columns(2),itr)
          tmp_pos2(3,npart1_tmp)= tmp_pos(columns(3),itr)
          IF( read_nu ) tmp_pos2(4,npart1_tmp)= tmp_pos(column_nu,itr)

        ENDIF

      ENDDO

      npart2_tmp= 0
      DO itr= 1, npart_tmp, 1

        IF( tmp_pos(columns(1),itr) > 0.0D0 )THEN

          npart2_tmp= npart2_tmp + 1
          tmp_pos2(1,npart1_tmp+npart2_tmp)= tmp_pos(columns(1),itr)
          tmp_pos2(2,npart1_tmp+npart2_tmp)= tmp_pos(columns(2),itr)
          tmp_pos2(3,npart1_tmp+npart2_tmp)= tmp_pos(columns(3),itr)
          IF( read_nu ) tmp_pos2(4,npart1_tmp+npart2_tmp)=tmp_pos(column_nu,itr)

        ENDIF

      ENDDO

      IF( npart1_tmp + npart2_tmp /= npart_tmp )THEN
        PRINT *, "** ERROR! parts_obj% npart1 + parts_obj% npart2 /= npart_tmp"
        PRINT *
        PRINT *, "   npart1_tmp= ", npart1_tmp
        PRINT *, "   npart2_tmp= ", npart2_tmp
        PRINT *, "   npart1_tmp + npart2_tmp= ", npart1_tmp + npart2_tmp
        PRINT *, "   npart_tmp= ", npart_tmp
        PRINT *
        STOP
      ENDIF

      ! Check that the positions are within the stars read from the LORENE
      ! binary file. This checks that the positions read from the formatted
      ! file are compatible with the binary file read

      ! Star 1
      IF( MINVAL( ABS( tmp_pos2(1,1:npart1_tmp) ) ) < ABS(center1) - &
                                           bns_obj% get_radius1_x_comp() &
          .OR. &
          MAXVAL( ABS( tmp_pos2(1,1:npart1_tmp) ) ) > ABS(center1) + &
                                           bns_obj% get_radius1_x_opp() &
          .OR. &
          ABS( MINVAL( tmp_pos2(2,1:npart1_tmp) ) ) > bns_obj% get_radius1_y() &
          .OR. &
          ABS( MAXVAL( tmp_pos2(2,1:npart1_tmp) ) ) > bns_obj% get_radius1_y() &
          .OR. &
          ABS( MINVAL( tmp_pos2(3,1:npart1_tmp) ) ) > bns_obj% get_radius1_z() &
          .OR. &
          ABS( MAXVAL( tmp_pos2(3,1:npart1_tmp) ) ) > bns_obj% get_radius1_z() &
          .OR. &
          MINVAL( ABS( tmp_pos2(1,1:npart1_tmp) ) ) > ABS(center1) - &
                                           0.95*bns_obj% get_radius1_x_comp() &
          .OR. &
          MAXVAL( ABS( tmp_pos2(1,1:npart1_tmp) ) ) < ABS(center1) + &
                                           0.95*bns_obj% get_radius1_x_opp() &
          .OR. &
          ABS( MINVAL( tmp_pos2(2,1:npart1_tmp) ) ) < &
                      0.95*bns_obj% get_radius1_y() &
          .OR. &
          ABS( MAXVAL( tmp_pos2(2,1:npart1_tmp) ) ) < &
                      0.95*bns_obj% get_radius1_y() &
          .OR. &
          ABS( MINVAL( tmp_pos2(3,1:npart1_tmp) ) ) < &
                      0.95*bns_obj% get_radius1_z() &
          .OR. &
          ABS( MAXVAL( tmp_pos2(3,1:npart1_tmp) ) ) < &
                      0.95*bns_obj% get_radius1_z() &

      )THEN

        PRINT *, "** ERROR! The positions of the particles on star 1, ", &
                 "read from file " &
                 // TRIM(parts_pos_namefile), " are not compatible with the ", &
                 "binary system read from the LORENE binary file. Stopping..."
        PRINT *
        STOP

      ENDIF

      ! Star 2
      IF( MINVAL( ABS( tmp_pos2(1,npart1_tmp:npart_tmp) ) ) < ABS(center2) - &
                                           bns_obj% get_radius2_x_comp() &
          .OR. &
          MAXVAL( ABS( tmp_pos2(1,npart1_tmp:npart_tmp) ) ) > ABS(center2) + &
                                           bns_obj% get_radius2_x_opp() &
          .OR. &
          ABS( MINVAL( tmp_pos2(2,npart1_tmp:npart_tmp) ) ) > &
                      bns_obj% get_radius2_y() &
          .OR. &
          ABS( MAXVAL( tmp_pos2(2,npart1_tmp:npart_tmp) ) ) > &
                      bns_obj% get_radius2_y() &
          .OR. &
          ABS( MINVAL( tmp_pos2(3,npart1_tmp:npart_tmp) ) ) > &
                      bns_obj% get_radius2_z() &
          .OR. &
          ABS( MAXVAL( tmp_pos2(3,npart1_tmp:npart_tmp) ) ) > &
                      bns_obj% get_radius2_z() &
          .OR. &
          MINVAL( ABS( tmp_pos2(1,npart1_tmp:npart_tmp) ) ) > ABS(center2) - &
                                           0.95*bns_obj% get_radius2_x_comp() &
          .OR. &
          MAXVAL( ABS( tmp_pos2(1,npart1_tmp:npart_tmp) ) ) < ABS(center2) + &
                                           0.95*bns_obj% get_radius2_x_opp() &
          .OR. &
          ABS( MINVAL( tmp_pos2(2,npart1_tmp:npart_tmp) ) ) < &
                      0.95*bns_obj% get_radius2_y() &
          .OR. &
          ABS( MAXVAL( tmp_pos2(2,npart1_tmp:npart_tmp) ) ) < &
                      0.95*bns_obj% get_radius2_y() &
          .OR. &
          ABS( MINVAL( tmp_pos2(3,npart1_tmp:npart_tmp) ) ) < &
                      0.95*bns_obj% get_radius2_z() &
          .OR. &
          ABS( MAXVAL( tmp_pos2(3,npart1_tmp:npart_tmp) ) ) < &
                      0.95*bns_obj% get_radius2_z() &

      )THEN

        PRINT *, "** ERROR! The positions of the particles on star 2, ", &
                 "read from file " &
                 // TRIM(parts_pos_namefile), " are not compatible with the ", &
                 "binary system read from the LORENE binary file. Stopping..."
        PRINT *
        STOP

      ENDIF

      !DO itr= 1, npart1_tmp, 1
      !  IF( tmp_pos2(1,itr) <  )
      !ENDDO

      ! Mirror the particles on star 1

      tmp_pos(columns(1),:)= tmp_pos2(1,:)
      tmp_pos(columns(2),:)= tmp_pos2(2,:)
      tmp_pos(columns(3),:)= tmp_pos2(3,:)
      IF( read_nu ) tmp_pos(column_nu,:) = tmp_pos2(4,:)

      parts_obj% npart1= 0
      DO itr= 1, npart1_tmp, 1

        IF( tmp_pos(columns(3),itr) > 0 )THEN

          parts_obj% npart1= parts_obj% npart1 + 1
          tmp_pos2(1,parts_obj% npart1)= tmp_pos(columns(1),itr)
          tmp_pos2(2,parts_obj% npart1)= tmp_pos(columns(2),itr)
          tmp_pos2(3,parts_obj% npart1)= tmp_pos(columns(3),itr)
          IF( read_nu ) tmp_pos2(4,parts_obj% npart1)= tmp_pos(column_nu,itr)

        ENDIF

      ENDDO

      DO itr= 1, parts_obj% npart1, 1

        tmp_pos2(1,parts_obj% npart1+itr)=   tmp_pos2(1,itr)
        tmp_pos2(2,parts_obj% npart1+itr)=   tmp_pos2(2,itr)
        tmp_pos2(3,parts_obj% npart1+itr)= - tmp_pos2(3,itr)
        IF( read_nu ) tmp_pos2(4,parts_obj% npart1+itr)= tmp_pos2(4,itr)

      ENDDO

      parts_obj% npart1= 2*parts_obj% npart1

      parts_obj% npart2= 0
      DO itr= npart1_tmp + 1, npart_tmp, 1

        IF( tmp_pos(columns(3),itr) > 0 )THEN

          parts_obj% npart2= parts_obj% npart2 + 1
          tmp_pos2(1,parts_obj% npart1+parts_obj% npart2)= tmp_pos(columns(1),itr)
          tmp_pos2(2,parts_obj% npart1+parts_obj% npart2)= tmp_pos(columns(2),itr)
          tmp_pos2(3,parts_obj% npart1+parts_obj% npart2)= tmp_pos(columns(3),itr)
          IF( read_nu ) tmp_pos2(4,parts_obj% npart1+parts_obj% npart2)= tmp_pos(column_nu,itr)

        ENDIF

      ENDDO

      DO itr= 1, parts_obj% npart2, 1

        tmp_pos2(1,parts_obj% npart1+parts_obj% npart2+itr)=   &
                                            tmp_pos2(1,parts_obj% npart1+itr)
        tmp_pos2(2,parts_obj% npart1+parts_obj% npart2+itr)=   &
                                            tmp_pos2(2,parts_obj% npart1+itr)
        tmp_pos2(3,parts_obj% npart1+parts_obj% npart2+itr)= &
                                          - tmp_pos2(3,parts_obj% npart1+itr)
        IF( read_nu ) tmp_pos2(4,parts_obj% npart1+parts_obj% npart2+itr)= &
                                            tmp_pos2(4,parts_obj% npart1+itr)

      ENDDO

      parts_obj% npart2= 2*parts_obj% npart2
      parts_obj% npart = parts_obj% npart1 + parts_obj% npart2

      !PRINT *, tmp_pos(:,1)
      ! Allocating the memory for the array pos( 3, npart )
      IF(.NOT.ALLOCATED( parts_obj% pos ))THEN
        ALLOCATE( parts_obj% pos( 3, parts_obj% npart ), STAT= ios, &
                  ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array pos in SUBROUTINE" &
                    // ". ", &
                    "The error message is", err_msg
           STOP
        ENDIF
        !CALL test_status( ios, err_msg, &
        !                "...allocation error for array pos in SUBROUTINE" &
        !                // "place_particles_3D_lattice." )
      ENDIF
      IF( read_nu .AND. .NOT.ALLOCATED( parts_obj% nu ))THEN
        ALLOCATE( parts_obj% nu( parts_obj% npart ), STAT= ios, &
                  ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array nu in SUBROUTINE" &
                    // ". ", &
                    "The error message is", err_msg
           STOP
        ENDIF
        !CALL test_status( ios, err_msg, &
        !                "...allocation error for array pos in SUBROUTINE" &
        !                // "place_particles_3D_lattice." )
      ENDIF
      IF( read_nu .AND. .NOT.ALLOCATED( parts_obj% pmass ))THEN
        ALLOCATE( parts_obj% pmass( parts_obj% npart ), STAT= ios, &
                  ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array pmass in SUBROUTINE" &
                    // ". ", &
                    "The error message is", err_msg
           STOP
        ENDIF
        !CALL test_status( ios, err_msg, &
        !                "...allocation error for array pos in SUBROUTINE" &
        !                // "place_particles_3D_lattice." )
      ENDIF

      !---------------------------------------------------------!
      !--  Storing the particle positions into the array pos  --!
      !--  symmetrically w.r.t. the xy plane                  --!
      !---------------------------------------------------------!

  !    ! Particles with z > 0 for star 1
  !    parts_obj% pos(:,1:parts_obj% npart1/2)= &
  !                                  tmp_pos2(:,1:parts_obj% npart1/2)
  !
  !    ! Particles with z < 0 for star 1
  !    parts_obj% pos(1:2,parts_obj% npart1/2+1:parts_obj% npart1)= &
  !                                  tmp_pos2(1:2,1:parts_obj% npart1/2)
  !
  !    parts_obj% pos(3,parts_obj% npart1/2+1:parts_obj% npart1)= &
  !                                - tmp_pos2(3,1:parts_obj% npart1/2)
  !
  !    ! Particles with z > 0 for star 2
  !    parts_obj% pos(:,parts_obj% npart1+1: &
  !                     parts_obj% npart1+parts_obj% npart2/2)= &
  !    tmp_pos2(:,parts_obj% npart1/2+1:parts_obj% npart1/2+parts_obj% npart2/2)
  !
  !    ! Particles with z < 0 for star 2
  !    parts_obj% pos(1:2,parts_obj% npart1+parts_obj% npart2/2+1: &
  !                     parts_obj% npart)= &
  !    tmp_pos2(1:2,parts_obj% npart1/2+1:parts_obj% npart1/2+parts_obj% npart2/2)
  !
  !    parts_obj% pos(3,parts_obj% npart1+parts_obj% npart2/2+1: &
  !                     parts_obj% npart)= &
  !    tmp_pos2(3,parts_obj% npart1/2+1:parts_obj% npart1/2+parts_obj% npart2/2)

      parts_obj% pos= tmp_pos2(1:3,1:parts_obj% npart)
      IF( read_nu ) parts_obj% nu= tmp_pos2(4,1:parts_obj% npart)

      PRINT *, " * Particle positions read. Number of particles=", &
               parts_obj% npart
      PRINT *
      PRINT *, " * Number of particles on NS 1=", parts_obj% npart1
      PRINT *, " * Number of particles on NS 2=", parts_obj% npart2
      PRINT *

      !
      !-- Computing volume per particle
      !
      IF(.NOT.ALLOCATED( parts_obj% pvol ))THEN
        ALLOCATE( parts_obj% pvol( parts_obj% npart ), STAT= ios, &
                ERRMSG= err_msg )
        IF( ios > 0 )THEN
          PRINT *, "...allocation error for array pvol ", &
                   ". The error message is", err_msg
          STOP
        ENDIF
        !CALL test_status( ios, err_msg, &
        !        "...allocation error for array v_euler_parts_z" )
      ENDIF

      ! First guess of the particle volume (it will be computed exactly later)

      pvol_tmp= 0.0D0
      DO itr= 1, parts_obj% npart - 1, 1

        pvol_tmp= pvol_tmp + ABS( parts_obj% pos(3,itr + 1) &
                                - parts_obj% pos(3,itr) )

      ENDDO
      pvol_tmp= pvol_tmp/( parts_obj% npart - 1 )

      parts_obj% pvol= 2.0D0*pvol_tmp**3.0D0

      IF( parts_obj% mass1 > parts_obj% mass2 )THEN

        ! mass_ratio < 1
        parts_obj% mass_ratio= parts_obj% mass2/parts_obj% mass1

      ELSE

        ! mass_ratio < 1
        parts_obj% mass_ratio= parts_obj% mass1/parts_obj% mass2

      ENDIF

      parts_obj% pmass= parts_obj% nu * amu

      !STOP

    CASE(1)

      PRINT *, " * Placing particles on one lattice around the stars."
      PRINT *

      !parts_obj% nx= nx
      !parts_obj% ny= ny
      !parts_obj% nz= nz

      !
      !-- Determine boundaries of the single lattice around the stars (Msun_geo)
      !
      xmin=   bns_obj% get_center1_x() - &
                                stretch*MAX( bns_obj% get_radius1_x_comp(), &
                                             bns_obj% get_radius1_x_opp() )
      xmax=   bns_obj% get_center2_x() + &
                                stretch*MAX( bns_obj% get_radius2_x_comp(), &
                                             bns_obj% get_radius2_x_opp() )
      ymin= - stretch*bns_obj% get_radius1_y()
      ymax=   stretch*bns_obj% get_radius2_y()
      zmin= - stretch*bns_obj% get_radius1_z()
      zmax=   stretch*bns_obj% get_radius2_z()

      ! Place particles, and time the process
      CALL parts_obj% placer_timer% start_timer()
      CALL parts_obj% place_particles_lattice( xmin, xmax, ymin, &
                                               ymax, zmin, zmax, &
                                               nx, ny, nz, &
                                               thres, bns_obj )
      CALL parts_obj% placer_timer% stop_timer()

    CASE(2)

      PRINT *, " * Placing particles on two lattices, " &
               // "one around each star."
      PRINT *

      !parts_obj% nx= nx
      !parts_obj% ny= ny
      !parts_obj% nz= nz

      !
      !-- Determine boundaries of the two lattices around the stars (Msun_geo)
      !
      xmin1=   bns_obj% get_center1_x() - &
                                stretch*MAX( bns_obj% get_radius1_x_comp(), &
                                             bns_obj% get_radius1_x_opp() )
      xmax1=   bns_obj% get_center1_x() + &
                                stretch*MAX( bns_obj% get_radius1_x_comp(), &
                                             bns_obj% get_radius1_x_opp() )
      ymin1= - stretch*bns_obj% get_radius1_y()
      ymax1=   stretch*bns_obj% get_radius1_y()
      zmin1= - stretch*bns_obj% get_radius1_z()
      zmax1=   stretch*bns_obj% get_radius1_z()

      xmin2=   bns_obj% get_center2_x() - &
                                stretch*MAX( bns_obj% get_radius2_x_comp(), &
                                             bns_obj% get_radius2_x_opp() )
      xmax2=   bns_obj% get_center2_x() + &
                                stretch*MAX( bns_obj% get_radius2_x_comp(), &
                                             bns_obj% get_radius2_x_opp() )
      ymin2= - stretch*bns_obj% get_radius2_y()
      ymax2=   stretch*bns_obj% get_radius2_y()
      zmin2= - stretch*bns_obj% get_radius2_z()
      zmax2=   stretch*bns_obj% get_radius2_z()

      ! Place particles, and time the process
      CALL parts_obj% placer_timer% start_timer()
      CALL parts_obj% place_particles_lattices( xmin1, xmax1, ymin1, &
                                                  ymax1, zmin1, zmax1, &
                                                  xmin2, xmax2, ymin2, &
                                                  ymax2, zmin2, zmax2, &
                                                  nx, ny, nz, &
                                                  thres, bns_obj )
      CALL parts_obj% placer_timer% stop_timer()

    CASE(3)

      PRINT *, "** Placing equal-mass particles on spherical surfaces, " &
               // "taking into account the mass profile of the stars."
      PRINT *

      ! Here the particle mass is computed using the radial mass profile
      ! of the star, so nu should not be redistributed to achieve a given
      ! particle mass ratio
      IF( parts_obj% redistribute_nu == .TRUE. )THEN
          parts_obj% redistribute_nu= .FALSE.
      ENDIF

      first_star_more_massive: IF( parts_obj% mass1 > parts_obj% mass2 )THEN

        filename_mass_profile= "spherical_surfaces_mass_profile2.dat"
        filename_shells_radii= "spherical_surfaces_radii2.dat"
        filename_shells_pos  = "spherical_surfaces_pos2.dat"

        n_particles_first_shell= 4

        ! Place particles, and time the process
        CALL parts_obj% placer_timer% start_timer()
        CALL parts_obj% place_particles_spherical_shells( parts_obj% mass2, &
                                                    radius2, center2, &
                                                    npart_approx, &
                                                    parts_obj% npart2, &
                                                    pos2, pvol2, pmass2, &
                                                    bns_obj, &
                                                    last_r, &
                                                    upper_bound, lower_bound, &
                                                    upper_factor, lower_factor,&
                                                    max_steps, &
                                                    filename_mass_profile, &
                                                    filename_shells_radii, &
                                                    filename_shells_pos )

        ! mass_ratio < 1
        parts_obj% mass_ratio= parts_obj% mass2/parts_obj% mass1

        equal_masses: IF( parts_obj% mass_ratio >= 0.995 .AND. &
            parts_obj% mass_ratio <= 1.005 .AND. reflect_particles_x )THEN

          IF(.NOT.ALLOCATED( pos1 ))THEN
            ALLOCATE( pos1( 3, parts_obj% npart2 ), &
                      STAT= ios, ERRMSG= err_msg )
            IF( ios > 0 )THEN
               PRINT *, "...allocation error for array pos in SUBROUTINE" &
                        // "place_particles_. ", &
                        "The error message is", err_msg
               STOP
            ENDIF
          ENDIF
          IF(.NOT.ALLOCATED( pvol1 ))THEN
            ALLOCATE( pvol1( parts_obj% npart2 ), &
                      STAT= ios, ERRMSG= err_msg )
            IF( ios > 0 )THEN
               PRINT *, "...allocation error for array pvol in SUBROUTINE" &
                        // "place_particles_. ", &
                        "The error message is", err_msg
               STOP
            ENDIF
          ENDIF
          IF(.NOT.ALLOCATED( pmass1 ))THEN
            ALLOCATE( pmass1( parts_obj% npart2 ), &
                      STAT= ios, ERRMSG= err_msg )
            IF( ios > 0 )THEN
               PRINT *, "...allocation error for array pmass in SUBROUTINE" &
                        // "place_particles_. ", &
                        "The error message is", err_msg
               STOP
            ENDIF
          ENDIF
          pos1(1,:)= - pos2(1,:)
          pos1(2,:)=   pos2(2,:)
          pos1(3,:)=   pos2(3,:)
          pvol1 = pvol2
          pmass1= pmass2
          parts_obj% npart1= parts_obj% npart2

        ELSE

          IF( parts_obj% mass_ratio >= 0.95 .AND. &
              parts_obj% mass_ratio <= 1.05 )THEN
            npart2_approx= npart_approx/parts_obj% mass_ratio
          ELSE
            npart2_approx= parts_obj% npart1/parts_obj% mass_ratio
          ENDIF

          filename_mass_profile= "spherical_surfaces_mass_profile1.dat"
          filename_shells_radii= "spherical_surfaces_radii1.dat"
          filename_shells_pos  = "spherical_surfaces_pos1.dat"

          n_particles_first_shell= n_particles_first_shell/parts_obj% mass_ratio

          CALL parts_obj% place_particles_spherical_shells( parts_obj% mass1, &
                                                radius1, center1, &
                                                npart2_approx, &
                                                parts_obj% npart1, &
                                                pos1, pvol1, pmass1, &
                                                bns_obj, &
                                                last_r, &
                                                upper_bound, lower_bound, &
                                                upper_factor, lower_factor,&
                                                max_steps, &
                                                filename_mass_profile, &
                                                filename_shells_radii, &
                                                filename_shells_pos )

        ENDIF equal_masses

        CALL parts_obj% placer_timer% stop_timer()

        parts_obj% npart= parts_obj% npart1 + parts_obj% npart2

      ELSE

        filename_mass_profile= "spherical_surfaces_mass_profile1.dat"
        filename_shells_radii= "spherical_surfaces_radii1.dat"
        filename_shells_pos  = "spherical_surfaces_pos1.dat"

        n_particles_first_shell= 4

        ! Place particles, and time the process
        CALL parts_obj% placer_timer% start_timer()

        !DO

        CALL parts_obj% place_particles_spherical_shells( parts_obj% mass1, &
                                              radius1, center1, &
                                              npart_approx, &
                                              parts_obj% npart1, &
                                              pos1, pvol1, pmass1, &
                                              bns_obj, &
                                              last_r, &
                                              upper_bound, lower_bound, &
                                              upper_factor, lower_factor,&
                                              max_steps, &
                                              filename_mass_profile, &
                                              filename_shells_radii, &
                                              filename_shells_pos )

        IF( debug ) PRINT *, "30"

        ! mass_ratio < 1
        parts_obj% mass_ratio= parts_obj% mass1/parts_obj% mass2

        IF( debug ) PRINT *, "31"

        equal_masses2: IF( parts_obj% mass_ratio >= 0.995 .AND. &
            parts_obj% mass_ratio <= 1.005 .AND. reflect_particles_x )THEN

          IF(.NOT.ALLOCATED( pos2 ))THEN
            ALLOCATE( pos2( 3, parts_obj% npart1 ), &
                      STAT= ios, ERRMSG= err_msg )
            IF( ios > 0 )THEN
               PRINT *, "...allocation error for array pos in SUBROUTINE" &
                        // "place_particles_. ", &
                        "The error message is", err_msg
               STOP
            ENDIF
          ENDIF
          IF(.NOT.ALLOCATED( pvol2 ))THEN
            ALLOCATE( pvol2( parts_obj% npart1 ), &
                      STAT= ios, ERRMSG= err_msg )
            IF( ios > 0 )THEN
               PRINT *, "...allocation error for array pvol in SUBROUTINE" &
                        // "place_particles_. ", &
                        "The error message is", err_msg
               STOP
            ENDIF
          ENDIF
          IF(.NOT.ALLOCATED( pmass2 ))THEN
            ALLOCATE( pmass2( parts_obj% npart1 ), &
                      STAT= ios, ERRMSG= err_msg )
            IF( ios > 0 )THEN
               PRINT *, "...allocation error for array pmass in SUBROUTINE" &
                        // "place_particles_. ", &
                        "The error message is", err_msg
               STOP
            ENDIF
          ENDIF
          pos2(1,:)= - pos1(1,:)
          pos2(2,:)=   pos1(2,:)
          pos2(3,:)=   pos1(3,:)
          pvol2 = pvol1
          pmass2= pmass1
          parts_obj% npart2= parts_obj% npart1

        ELSE

          IF( parts_obj% mass_ratio >= 0.95 .AND. &
              parts_obj% mass_ratio <= 1.05 )THEN
            npart2_approx= npart_approx/parts_obj% mass_ratio
          ELSE
            npart2_approx= parts_obj% npart1/parts_obj% mass_ratio
          ENDIF

          filename_mass_profile= "spherical_surfaces_mass_profile2.dat"
          filename_shells_radii= "spherical_surfaces_radii2.dat"
          filename_shells_pos  = "spherical_surfaces_pos2.dat"

          n_particles_first_shell= 4!n_particles_first_shell/parts_obj% mass_ratio

          IF( debug ) PRINT *, "32"

          CALL parts_obj% place_particles_spherical_shells( parts_obj% mass2, &
                                                radius2, center2, &
                                                npart2_approx, &
                                                parts_obj% npart2, &
                                                pos2, pvol2, pmass2, &
                                                bns_obj, &
                                                last_r, &
                                                upper_bound, lower_bound, &
                                                upper_factor, lower_factor,&
                                                max_steps, &
                                                filename_mass_profile, &
                                                filename_shells_radii, &
                                                filename_shells_pos )

        ENDIF equal_masses2

        CALL parts_obj% placer_timer% stop_timer()

        parts_obj% npart= parts_obj% npart1 + parts_obj% npart2

      ENDIF first_star_more_massive

      !
      !-- Assign TYPE member variables
      !

      IF(.NOT.ALLOCATED( parts_obj% pos ))THEN
        ALLOCATE( parts_obj% pos( 3, parts_obj% npart ), &
                  STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array pos in SUBROUTINE" &
                    // "place_particles_. ", &
                    "The error message is", err_msg
           STOP
        ENDIF
        !CALL test_status( ios, err_msg, &
        !                "...allocation error for array pos in SUBROUTINE" &
        !                // "place_particles_3D_lattice." )
      ENDIF
      parts_obj% pos( :, 1:parts_obj% npart1 )= pos1
      parts_obj% pos( :, parts_obj% npart1 + 1:parts_obj% npart )= pos2

      IF(.NOT.ALLOCATED( parts_obj% pvol ))THEN
        ALLOCATE( parts_obj% pvol( parts_obj% npart ), &
                  STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array pvol in SUBROUTINE" &
                    // "place_particles_. ", &
                    "The error message is", err_msg
           STOP
        ENDIF
        !CALL test_status( ios, err_msg, &
        !                "...allocation error for array pos in SUBROUTINE" &
        !                // "place_particles_3D_lattice." )
      ENDIF
      parts_obj% pvol( 1:parts_obj% npart1 )= pvol1
      parts_obj% pvol( parts_obj% npart1 + 1:parts_obj% npart )= pvol2

      IF(.NOT.ALLOCATED( parts_obj% pmass ))THEN
        ALLOCATE( parts_obj% pmass( parts_obj% npart ), &
                  STAT= ios, ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array pmass in SUBROUTINE" &
                    // "place_particles_. ", &
                    "The error message is", err_msg
           STOP
        ENDIF
        !CALL test_status( ios, err_msg, &
        !                "...allocation error for array pos in SUBROUTINE" &
        !                // "place_particles_3D_lattice." )
      ENDIF
      parts_obj% pmass( 1:parts_obj% npart1 )= pmass1
      parts_obj% pmass( parts_obj% npart1 + 1:parts_obj% npart )= pmass2

      PRINT *, " * Particles placed. Number of particles=", parts_obj% npart
      PRINT *, " * Number of particles on NS 1=", parts_obj% npart1
      PRINT *, " * Number of particles on NS 2=", parts_obj% npart2
      PRINT *
      !STOP

    CASE DEFAULT

      PRINT *, "** There is no implemented particle placer " &
               // "corresponding to the number", dist
      PRINT *, " * Stopping..."
      STOP

    END SELECT choose_particle_placer

    ! Reshape the arrays pos and pvol by deleting the unnecessary elements
    parts_obj% pos = parts_obj% pos( :, 1:parts_obj% npart )
    parts_obj% pvol= parts_obj% pvol( 1:parts_obj% npart )

    ! Check that there aren't particles with the same coordinates
    CALL parts_obj% same_particle_timer% start_timer()
    CALL check_particle_positions( parts_obj% npart, parts_obj% pos )
    CALL parts_obj% same_particle_timer% stop_timer()

    IF( apm_iterate1 )THEN

      PRINT *
      PRINT *, "** Placing particles on star 1 using the APM..."
      PRINT *

      IF(.NOT.ALLOCATED( parts_obj% h ))THEN
        ALLOCATE( parts_obj% h( parts_obj% npart ), STAT= ios, &
                  ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array h in SUBROUTINE ", &
                    "construct_particles. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF

      IF(.NOT.ALLOCATED( parts_obj% nu ))THEN
        ALLOCATE( parts_obj% nu( parts_obj% npart ), STAT= ios, &
                  ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array nu in SUBROUTINE ", &
                    "construct_particles. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF

      filename_apm_pos_id = "apm_pos_id1.dat"
      filename_apm_pos    = "apm_pos1.dat"
      filename_apm_results= "apm_results1.dat"

      ! Star 1
      CALL parts_obj% apm1_timer% start_timer()
      CALL parts_obj% perform_apm( &
                  bns_obj, &
                  parts_obj% pos(:,1:parts_obj% npart1), &
                  parts_obj% pvol(1:parts_obj% npart1), &
                  parts_obj% h(1:parts_obj% npart1), &
                  parts_obj% nu(1:parts_obj% npart1), &
                  center1, com1, parts_obj% mass1, &
                  apm_max_it, max_inc, mass_it, parts_obj% correct_nu, &
                  nuratio_thres, nuratio_des, nx_gh, ny_gh, nz_gh, &
                  filename_apm_pos_id, filename_apm_pos, filename_apm_results )
      CALL parts_obj% apm1_timer% stop_timer()

      PRINT *, "** Particles placed on star 1 according to the APM."
      PRINT *

      IF( parts_obj% mass_ratio >= 0.995 .AND. &
          parts_obj% mass_ratio <= 1.005 .AND. reflect_particles_x )THEN

        parts_obj% pos(1,parts_obj% npart1+1:parts_obj% npart)= &
                                  - parts_obj% pos(1,1:parts_obj% npart1)
        parts_obj% pos(2,parts_obj% npart1+1:parts_obj% npart)= &
                                    parts_obj% pos(2,1:parts_obj% npart1)
        parts_obj% pos(3,parts_obj% npart1+1:parts_obj% npart)= &
                                    parts_obj% pos(3,1:parts_obj% npart1)

        parts_obj% nu(parts_obj% npart1+1:parts_obj% npart)= &
                                    parts_obj% nu(1:parts_obj% npart1)

        parts_obj% h(parts_obj% npart1+1:parts_obj% npart)= &
                                    parts_obj% h(1:parts_obj% npart1)

        parts_obj% npart2= parts_obj% npart1
        parts_obj% npart= parts_obj% npart1 + parts_obj% npart1

      !ELSEIF( ( parts_obj% mass_ratio <= 0.995 .OR. &
      !        parts_obj% mass_ratio >= 1.005 ) .AND. reflect_particles_x )THEN
      !
      !  PRINT *, "** ERROR! The two stars are not the same. The particles", &
      !           " on star 1 cannot be reflected with respect to the yz ", &
      !           " plane to become the particles star 2."
      !  PRINT *, "   Please, choose an equal-mass system, or set the ", &
      !           "   variable reflect_particles_x to .FALSE. in the file", &
      !           "   lorene_bns_id_particles ."
      !  PRINT *
      !  STOP
      PRINT *, "** Particles placed on star 1 according to the APM", &
               " reflected about the yz plane onto star 2."
      PRINT *

      ENDIF
    ENDIF
    IF( apm_iterate2 .AND. .NOT.(parts_obj% mass_ratio >= 0.995 .AND. &
        parts_obj% mass_ratio <= 1.005 .AND. reflect_particles_x) )THEN

      PRINT *
      PRINT *, "** Placing particles on star 2 using the APM..."
      PRINT *

      IF(.NOT.ALLOCATED( parts_obj% h ))THEN
        ALLOCATE( parts_obj% h( parts_obj% npart ), STAT= ios, &
                  ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array h in SUBROUTINE ", &
                    "construct_particles. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF

      IF(.NOT.ALLOCATED( parts_obj% nu ))THEN
        ALLOCATE( parts_obj% nu( parts_obj% npart ), STAT= ios, &
                  ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array nu in SUBROUTINE ", &
                    "construct_particles. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF

      filename_apm_pos_id = "apm_pos_id2.dat"
      filename_apm_pos    = "apm_pos2.dat"
      filename_apm_results= "apm_results2.dat"

      ! Star 2
      CALL parts_obj% apm2_timer% start_timer()
      CALL parts_obj% perform_apm( &
                bns_obj, &
                parts_obj% pos(:,parts_obj% npart1+1:parts_obj% npart), &
                parts_obj% pvol(parts_obj% npart1+1:parts_obj% npart), &
                parts_obj% h(parts_obj% npart1+1:parts_obj% npart), &
                parts_obj% nu(parts_obj% npart1+1:parts_obj% npart), &
                center2, com2, parts_obj% mass2, &
                apm_max_it, max_inc, mass_it, parts_obj% correct_nu, &
                nuratio_thres, nuratio_des, nx_gh, ny_gh, nz_gh, &
                filename_apm_pos_id, filename_apm_pos, filename_apm_results )
      CALL parts_obj% apm2_timer% stop_timer()

      PRINT *, "** Particles placed on star 2 according to the APM."
      PRINT *

    ENDIF

    ! Allocate needed memory
    CALL allocate_lorene_id_parts_memory( parts_obj )

    ! flag that particles are 'alive'
    ALLOCATE( alive( parts_obj% npart ) )
    alive( 1:parts_obj% npart )= 1

    IF( debug ) PRINT *, "33"

    !
    !-- Import the needed LORENE ID on the particles, and time the process
    !
    PRINT *, "** Importing the LORENE ID on the particles..."

    CALL parts_obj% importer_timer% start_timer()
    CALL bns_obj% import_id( parts_obj% npart, &
                             parts_obj% pos( 1, : ), &
                             parts_obj% pos( 2, : ), &
                             parts_obj% pos( 3, : ), &
                             parts_obj% lapse_parts, &
                             parts_obj% shift_parts_x, &
                             parts_obj% shift_parts_y, &
                             parts_obj% shift_parts_z, &
                             parts_obj% g_xx_parts, &
                             parts_obj% g_xy_parts, &
                             parts_obj% g_xz_parts, &
                             parts_obj% g_yy_parts, &
                             parts_obj% g_yz_parts, &
                             parts_obj% g_zz_parts, &
                             parts_obj% baryon_density_parts, &
                             parts_obj% energy_density_parts, &
                             parts_obj% specific_energy_parts, &
                             parts_obj% pressure_parts, &
                             parts_obj% v_euler_parts_x, &
                             parts_obj% v_euler_parts_y, &
                             parts_obj% v_euler_parts_z )
    CALL parts_obj% importer_timer% stop_timer()

    IF( debug ) PRINT *, "34"

    !
    !-- Check that the imported ID does not contain NaNs
    !
    !CALL Check_Array_for_NAN( parts_obj% npart, parts_obj% lapse_parts, &
    !                                         "lapse_parts" )
    !CALL Check_Array_for_NAN( parts_obj% npart, parts_obj% shift_parts_x, &
    !                                         "shift_parts_x" )
    !CALL Check_Array_for_NAN( parts_obj% npart, parts_obj% shift_parts_y, &
    !                                         "shift_parts_y" )
    !CALL Check_Array_for_NAN( parts_obj% npart, parts_obj% shift_parts_z, &
    !                                         "shift_parts_z" )
    !CALL Check_Array_for_NAN( parts_obj% npart, parts_obj% g_xx_parts, &
    !                                         "g_xx_parts" )
    !CALL Check_Array_for_NAN( parts_obj% npart, parts_obj% g_xy_parts, &
    !                                         "g_xy_parts" )
    !CALL Check_Array_for_NAN( parts_obj% npart, parts_obj% g_xz_parts, &
    !                                         "g_xz_parts" )
    !CALL Check_Array_for_NAN( parts_obj% npart, parts_obj% g_yy_parts, &
    !                                         "g_yy_parts" )
    !CALL Check_Array_for_NAN( parts_obj% npart, parts_obj% g_yz_parts, &
    !                                         "g_yz_parts" )
    !CALL Check_Array_for_NAN( parts_obj% npart, parts_obj% g_zz_parts, &
    !                                         "g_zz_parts" )
    !CALL Check_Array_for_NAN( parts_obj% npart, &
    !        parts_obj% baryon_density_parts, "baryon_density_parts" )
    !CALL Check_Array_for_NAN( parts_obj% npart, &
    !        parts_obj% energy_density_parts, "energy_density_parts" )
    !CALL Check_Array_for_NAN( parts_obj% npart, &
    !        parts_obj% specific_energy_parts, "specific_energy_parts" )
    !CALL Check_Array_for_NAN( parts_obj% npart, &
    !               parts_obj% pressure_parts, "pressure_parts" )
    !CALL Check_Array_for_NAN( parts_obj% npart, &
    !              parts_obj% v_euler_parts_x, "v_euler_parts_x" )
    !CALL Check_Array_for_NAN( parts_obj% npart, &
    !              parts_obj% v_euler_parts_y, "v_euler_parts_y" )
    !CALL Check_Array_for_NAN( parts_obj% npart, &
    !              parts_obj% v_euler_parts_z, "v_euler_parts_z" )

    IF(.NOT.ALLOCATED( parts_obj% baryon_density_index ))THEN
      ALLOCATE( parts_obj% baryon_density_index( parts_obj% npart ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array baryon_density_index in " &
                  // "SUBROUTINE construct_particles. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array pos in SUBROUTINE" &
      !                // "place_particles_3D_lattice." )
    ENDIF

    !PRINT *, "baryon_density_index"
    !DO itr= 1, parts_obj% npart1, 1
    !  PRINT *, parts_obj% baryon_density_index( itr ), &
    !           parts_obj% baryon_density_parts( itr )
    !ENDDO
    !    PRINT *, "baryon_density_parts in ascending order"
    !DO itr= 1, parts_obj% npart1, 1
    !  PRINT *, parts_obj% baryon_density_parts( &
    !                                parts_obj% baryon_density_index( itr ) )
    !ENDDO
    !PRINT *, "baryon_density_parts in descending order"
    !DO itr= parts_obj% npart1, 1, -1
    !  PRINT *, parts_obj% baryon_density_parts( &
    !                                parts_obj% baryon_density_index( itr ) )
    !ENDDO
    ! Ok it seems working

!    PRINT *, "** Computing typical length scale for the change in pressure", &
!             " on the x axis."
!    PRINT *
!
!    ALLOCATE( abs_pos( 3, parts_obj% npart ) )
!
!    DO itr = 1, parts_obj% npart, 1
!      abs_pos( 1, itr )= ABS( parts_obj% pos( 1, itr ) )
!      abs_pos( 2, itr )= ABS( parts_obj% pos( 2, itr ) )
!      abs_pos( 3, itr )= ABS( parts_obj% pos( 3, itr ) )
!    ENDDO
!
!    min_y_index= 0
!    min_abs_y= 1D+20
!    DO itr = 1, parts_obj% npart, 1
!      IF( ABS( parts_obj% pos( 2, itr ) ) < min_abs_y )THEN
!        min_abs_y= ABS( parts_obj% pos( 2, itr ) )
!        min_y_index= itr
!      ENDIF
!    ENDDO
!
!    min_z_index= 0
!    min_abs_z= 1D+20
!    DO itr = 1, parts_obj% npart, 1
!      IF( ABS( parts_obj% pos( 3, itr ) ) < min_abs_z )THEN
!        min_abs_z= ABS( parts_obj% pos( 3, itr ) )
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
!    DO itr = 1, parts_obj% npart, 1
!      IF( parts_obj% pos( 3, itr ) == min_abs_z &
!          .AND. &
!          ABS( ( parts_obj% pos( 2, itr ) - &
!                 parts_obj% pos( 2, min_y_index ) )/ &
!                 parts_obj% pos( 2, min_y_index ) ) < 1.0D-5 &
!      )THEN
!
!        IF( parts_obj% pos( 1, itr ) < 0 )THEN
!          cntr1= cntr1 + 1
!        ELSEIF( parts_obj% pos( 1, itr ) > 0 )THEN
!          cntr2= cntr2 + 1
!        ENDIF
!
!      ENDIF
!    ENDDO
!    !PRINT *, "cntr1= ", cntr1
!    !PRINT *, "cntr2= ", cntr2
!
!    ALLOCATE( parts_obj% pos_x1( cntr1 ) )
!    ALLOCATE( parts_obj% pos_x2( cntr2 ) )
!    ALLOCATE( parts_obj% pressure_parts_x1( cntr1 ) )
!    ALLOCATE( parts_obj% pressure_parts_x2( cntr2 ) )
!    ALLOCATE( parts_obj% pressure_parts_x_der1( cntr1 - 5 ) )
!    ALLOCATE( parts_obj% pressure_parts_x_der2( cntr2 - 5 ) )
!    ALLOCATE( parts_obj% pressure_length_scale_x1( cntr1 - 5 ) )
!    ALLOCATE( parts_obj% pressure_length_scale_x2( cntr2 - 5 ) )
!
!    !PRINT *, "2"
!
!    itr_1= 0
!    itr_2= 0
!    DO itr = 1, parts_obj% npart, 1
!      IF( parts_obj% pos( 3, itr ) == min_abs_z &
!          .AND. &
!          ABS( ( parts_obj% pos( 2, itr ) - &
!                 parts_obj% pos( 2, min_y_index ) )/ &
!                 parts_obj% pos( 2, min_y_index ) ) < 1.0D-5 &
!        )THEN
!
!        IF( parts_obj% pos( 1, itr ) < 0 )THEN
!          itr_1= itr_1 + 1
!          parts_obj% pos_x1( itr_1 )= parts_obj% pos( 1, itr )
!          parts_obj% pressure_parts_x1( itr_1 )= &
!                                              parts_obj% pressure_parts( itr )
!        ELSEIF( parts_obj% pos( 1, itr ) > 0 )THEN
!          itr_2= itr_2 + 1
!          parts_obj% pos_x2( itr_2 )= parts_obj% pos( 1, itr )
!          parts_obj% pressure_parts_x2( itr_2 )= &
!                                              parts_obj% pressure_parts( itr )
!        ENDIF
!
!      ENDIF
!    ENDDO
!
!    !PRINT *, "3"
!
!    DO itr= 3, cntr1 - 3, 1
!      parts_obj% pressure_parts_x_der1( itr - 2 )=&
!                     ( + parts_obj% pressure_parts_x1( itr - 2 )/12.0D0 &
!                       - 2.0*parts_obj% pressure_parts_x1( itr - 1 )/3.0D0 &
!                       + 2.0*parts_obj% pressure_parts_x1( itr + 1 )/3.0D0 &
!                       - parts_obj% pressure_parts_x1( itr + 2 )/12.0D0 )&
!                       /( Msun_geo*km2m*ABS( parts_obj% pos_x1( itr ) - &
!                                             parts_obj% pos_x1( itr - 1 ) ) )
!
!      parts_obj% pressure_length_scale_x1( itr - 2 )= &
!                          ABS( parts_obj% pressure_parts_x1( itr - 2 )/ &
!                               parts_obj% pressure_parts_x_der1( itr - 2 ) )
!
!      !PRINT *, "p1=", parts_obj% pressure_parts_x1( itr - 2 )
!      !PRINT *, "p_r1=", parts_obj% pressure_parts_x_der1( itr - 2 )
!      !PRINT *, "p/p_r1=", parts_obj% pressure_length_scale_x1( itr - 2 )
!      !PRINT *
!
!    ENDDO
!    DO itr= 3, cntr2 - 3, 1
!      parts_obj% pressure_parts_x_der2( itr - 2 )=&
!                     ( + parts_obj% pressure_parts_x2( itr - 2 )/12.0D0 &
!                       - 2.0*parts_obj% pressure_parts_x2( itr - 1 )/3.0D0 &
!                       + 2.0*parts_obj% pressure_parts_x2( itr + 1 )/3.0D0 &
!                       - parts_obj% pressure_parts_x2( itr + 2 )/12.0D0 )&
!                       /( Msun_geo*km2m*ABS( parts_obj% pos_x2( itr ) - &
!                                             parts_obj% pos_x2( itr - 1 ) ) )
!
!      parts_obj% pressure_length_scale_x2( itr - 2 )= &
!                          ABS( parts_obj% pressure_parts_x2( itr - 2 )/ &
!                               parts_obj% pressure_parts_x_der2( itr - 2 ) )
!
!      !PRINT *, "p2=", parts_obj% pressure_parts_x2( itr - 2 )
!      !PRINT *, "p_r2=", parts_obj% pressure_parts_x_der2( itr - 2 )
!      !PRINT *, "p/p_r2=", parts_obj% pressure_length_scale_x2( itr - 2 )
!      !PRINT *
!
!    ENDDO
!
!    PRINT *, " * Maximum typical length scale for change in pressure", &
!             " along the x axis for NS 1= ", &
!             MAXVAL( parts_obj% pressure_length_scale_x1, DIM= 1 )/km2m, " km"
!    PRINT *, " * Minimum typical length scale for change in pressure", &
!             " along the x axis for NS 1= ", &
!             MINVAL( parts_obj% pressure_length_scale_x1, DIM= 1 )/km2m, " km"
!    PRINT *
!    PRINT *, " * Maximum typical length scale for change in pressure", &
!             " along the x axis for NS 2= ", &
!             MAXVAL( parts_obj% pressure_length_scale_x2, DIM= 1 )/km2m, " km"
!    PRINT *, " * Minimum typical length scale for change in pressure", &
!             " along the x axis for NS 2= ", &
!             MINVAL( parts_obj% pressure_length_scale_x2, DIM= 1 )/km2m, " km"
!    PRINT *

    ! Increase the counter that identifies the particle distribution
    counter= counter + 1

    !PRINT *, "End of particle constructor"

    IF( parts_obj% redistribute_nu )THEN

      ! Index particles on star 1 in increasing order of nu

      CALL indexx( parts_obj% npart1, &
                   parts_obj% baryon_density_parts( 1 : parts_obj% npart1 ), &
                   parts_obj% baryon_density_index( 1 : parts_obj% npart1 ) )

      ! Index particles on star 2 in increasing order of nu

      CALL indexx( parts_obj% npart2, &
                   parts_obj% baryon_density_parts( parts_obj% npart1 + 1 : &
                                                    parts_obj% npart ), &
                   parts_obj% baryon_density_index( parts_obj% npart1 + 1 : &
                                                    parts_obj% npart ) )

      ! Shift indices on star 2 by npart1 since all the arrays store
      ! the quantities on star 1 first, and then on star 2

      parts_obj% baryon_density_index( parts_obj% npart1 + 1 : &
                                       parts_obj% npart )= &
                     parts_obj% npart1 + &
                     parts_obj% baryon_density_index( parts_obj% npart1 + 1 : &
                                                      parts_obj% npart )

    ENDIF

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
 !     "# Values of the fields (including coordinates) exported by LORENE "&
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
 !     DO itr = 1, parts_obj% npart, 1
 !       abs_pos( 1, itr )= ABS( parts_obj% pos( 1, itr ) )
 !       abs_pos( 2, itr )= ABS( parts_obj% pos( 2, itr ) )
 !       abs_pos( 3, itr )= ABS( parts_obj% pos( 3, itr ) )
 !     ENDDO
 !
 !     min_y_index= 0
 !     min_abs_y= 1D+20
 !     DO itr = 1, parts_obj% npart, 1
 !       IF( ABS( parts_obj% pos( 2, itr ) ) < min_abs_y )THEN
 !         min_abs_y= ABS( parts_obj% pos( 2, itr ) )
 !         min_y_index= itr
 !       ENDIF
 !     ENDDO
 !
 !     min_abs_z= MINVAL( abs_pos( 3, : ) )
 !
 !     write_data_loop: DO itr = 1, parts_obj% npart, 1
 !
 !       IF( parts_obj% export_form_xy .AND. &
 !           parts_obj% pos( 3, itr ) /= min_abs_z )THEN
 !         CYCLE
 !       ENDIF
 !       IF( parts_obj% export_form_x .AND. &
 !           ( parts_obj% pos( 3, itr ) /= min_abs_z &
 !             .OR. &
 !             parts_obj% pos( 2, itr ) /= parts_obj% pos( 2, min_y_index ) ) &
 !       )THEN
 !         CYCLE
 !       ENDIF
 !       WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
 !         itr, &
 !         parts_obj% pos( 1, itr ), &
 !         parts_obj% pos( 2, itr ), &
 !         parts_obj% pos( 3, itr ), &
 !         parts_obj% lapse_parts( itr ), &
 !         parts_obj% shift_parts_x( itr ), &
 !         parts_obj% shift_parts_y( itr ), &
 !         parts_obj% shift_parts_z( itr ), &
 !         parts_obj% baryon_density_parts( itr ), &
 !         parts_obj% energy_density_parts( itr ), &
 !         parts_obj% specific_energy_parts( itr ), &
 !         parts_obj% pressure_parts( itr ), &
 !         parts_obj% v_euler_parts_x( itr ), &
 !         parts_obj% v_euler_parts_y( itr ), &
 !         parts_obj% v_euler_parts_z( itr )
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


END SUBMODULE particles_constructor