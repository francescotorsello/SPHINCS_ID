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


  !MODULE PROCEDURE construct_particles_bnsbase_empty
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
  !END PROCEDURE construct_particles_bnsbase_empty


  MODULE PROCEDURE construct_particles_bnsbase

    !**************************************************
    !
    !# The constructor performs all the tasks needed
    !  to set up the particle distribution with the
    !  |lorene| ID on it. It calls all the PROCEDURES
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
    USE tensor,         ONLY: jx, jy, jz

    IMPLICIT NONE

    ! The variable counter counts how many times the PROCEDURE
    ! construct_particles_bnsbase is called
    INTEGER, SAVE:: counter= 1
    INTEGER:: nx, ny, nz, &
              npart_approx, npart2_approx, max_steps, &
              nlines, header_lines, n_cols, npart_tmp, npart1_tmp, npart2_tmp, &
              nx_gh, ny_gh, nz_gh
    ! Maximum length for strings, and for the number of imported binaries
    INTEGER, PARAMETER:: max_length= 50
    ! APM parameters
    INTEGER:: apm_max_it, max_inc!, n_particles_first_shell
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
    DOUBLE PRECISION:: central_density1, central_density2
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
    parts% placer_timer       = timer( "placer_timer" )
    parts% apm1_timer         = timer( "apm_star1_timer" )
    parts% apm2_timer         = timer( "apm_star2_timer" )
    parts% importer_timer     = timer( "importer_timer" )
    parts% sph_computer_timer = timer( "sph_computer_timer" )
    parts% same_particle_timer= timer( "same_particle_timer" )

    ! Declare this object as non-empty (experimental)
    parts% empty_object= .FALSE.

    parts% mass1          = bnsb% get_mass1()
    parts% mass2          = bnsb% get_mass2()
    center1                   = bnsb% get_center1_x()
    center2                   = bnsb% get_center2_x()
    central_density1          = bnsb% get_rho_center1()
    central_density2          = bnsb% get_rho_center2()
    com1                      = bnsb% get_barycenter1_x()
    com2                      = bnsb% get_barycenter2_x()
    radius1                   = bnsb% get_radius1_x_comp()
    radius2                   = bnsb% get_radius2_x_comp()
    parts% nbar_tot       = 0.0D0
    parts% nbar1          = 0.0D0
    parts% nbar2          = 0.0D0
    parts% npart          = 0.0D0
    parts% distribution_id= dist

    parts% eos1= bnsb% get_eos1()
    parts% eos2= bnsb% get_eos2()

    parts% eos1_id= bnsb% get_eos1_id()
    parts% eos2_id= bnsb% get_eos2_id()

    parts% gamma_sp1= bnsb% get_gamma_1()
    parts% kappa_sp1= bnsb% get_kappa_1()
    parts% gamma_sp2= bnsb% get_gamma_2()
    parts% kappa_sp2= bnsb% get_kappa_2()

    !
    !-- Read the parameters of the particle distributions
    !
    parts% lorene_bns_id_parfile= 'sphincs_lorene_bns_particles.par'

    INQUIRE( FILE= parts% lorene_bns_id_parfile, EXIST= file_exists )
    IF( file_exists )THEN
     OPEN( 10, FILE= parts% lorene_bns_id_parfile, STATUS= 'OLD' )
    ELSE
     PRINT *
     PRINT *, "** ERROR: ", parts% lorene_bns_id_parfile, &
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
    parts% nu_ratio           = nu_ratio
    parts% reflect_particles_x= reflect_particles_x
    parts% randomize_phi      = randomize_phi
    parts% randomize_theta    = randomize_theta
    parts% randomize_r        = randomize_r
    ! APM parameters
    parts% apm_iterate1   = apm_iterate1
    parts% apm_iterate2   = apm_iterate2
    !parts% apm_max_it   = apm_max_it
    !parts% max_inc      = max_inc
    !parts% mass_it      = mass_it
    !parts% nuratio_thres= nuratio_thres
    parts% read_nu       = read_nu

    parts_pos_namefile= TRIM(parts_pos_path)//TRIM(parts_pos)

    IF( parts% redistribute_nu )THEN
      thres= 100.0D0*parts% nu_ratio
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
        PRINT *, "** ERROR! parts% npart1 + parts% npart2 /= npart_tmp"
        PRINT *
        PRINT *, "   npart1_tmp= ", npart1_tmp
        PRINT *, "   npart2_tmp= ", npart2_tmp
        PRINT *, "   npart1_tmp + npart2_tmp= ", npart1_tmp + npart2_tmp
        PRINT *, "   npart_tmp= ", npart_tmp
        PRINT *
        STOP
      ENDIF

      ! Check that the positions are within the stars read from the |lorene|
      ! binary file. This checks that the positions read from the formatted
      ! file are compatible with the binary file read

      ! Star 1
      IF( MINVAL( ABS( tmp_pos2(1,1:npart1_tmp) ) ) < ABS(center1) - &
                                           bnsb% get_radius1_x_comp() &
          .OR. &
          MAXVAL( ABS( tmp_pos2(1,1:npart1_tmp) ) ) > ABS(center1) + &
                                           bnsb% get_radius1_x_opp() &
          .OR. &
          ABS( MINVAL( tmp_pos2(2,1:npart1_tmp) ) ) > bnsb% get_radius1_y() &
          .OR. &
          ABS( MAXVAL( tmp_pos2(2,1:npart1_tmp) ) ) > bnsb% get_radius1_y() &
          .OR. &
          ABS( MINVAL( tmp_pos2(3,1:npart1_tmp) ) ) > bnsb% get_radius1_z() &
          .OR. &
          ABS( MAXVAL( tmp_pos2(3,1:npart1_tmp) ) ) > bnsb% get_radius1_z() &
          .OR. &
          MINVAL( ABS( tmp_pos2(1,1:npart1_tmp) ) ) > ABS(center1) - &
                                           0.95*bnsb% get_radius1_x_comp() &
          .OR. &
          MAXVAL( ABS( tmp_pos2(1,1:npart1_tmp) ) ) < ABS(center1) + &
                                           0.95*bnsb% get_radius1_x_opp() &
          .OR. &
          ABS( MINVAL( tmp_pos2(2,1:npart1_tmp) ) ) < &
                      0.95*bnsb% get_radius1_y() &
          .OR. &
          ABS( MAXVAL( tmp_pos2(2,1:npart1_tmp) ) ) < &
                      0.95*bnsb% get_radius1_y() &
          .OR. &
          ABS( MINVAL( tmp_pos2(3,1:npart1_tmp) ) ) < &
                      0.95*bnsb% get_radius1_z() &
          .OR. &
          ABS( MAXVAL( tmp_pos2(3,1:npart1_tmp) ) ) < &
                      0.95*bnsb% get_radius1_z() &

      )THEN

        PRINT *, "** ERROR! The positions of the particles on star 1, ", &
                 "read from file " &
                 // TRIM(parts_pos_namefile), " are not compatible with the ", &
                 "binary system read from the |lorene| binary file. Stopping..."
        PRINT *
        STOP

      ENDIF

      ! Star 2
      IF( MINVAL( ABS( tmp_pos2(1,npart1_tmp+1:npart_tmp) ) ) < ABS(center2) - &
                                           bnsb% get_radius2_x_comp() &
          .OR. &
          MAXVAL( ABS( tmp_pos2(1,npart1_tmp+1:npart_tmp) ) ) > ABS(center2) + &
                                           bnsb% get_radius2_x_opp() &
          .OR. &
          ABS( MINVAL( tmp_pos2(2,npart1_tmp+1:npart_tmp) ) ) > &
                      bnsb% get_radius2_y() &
          .OR. &
          ABS( MAXVAL( tmp_pos2(2,npart1_tmp+1:npart_tmp) ) ) > &
                      bnsb% get_radius2_y() &
          .OR. &
          ABS( MINVAL( tmp_pos2(3,npart1_tmp+1:npart_tmp) ) ) > &
                      bnsb% get_radius2_z() &
          .OR. &
          ABS( MAXVAL( tmp_pos2(3,npart1_tmp+1:npart_tmp) ) ) > &
                      bnsb% get_radius2_z() &
          .OR. &
          MINVAL( ABS( tmp_pos2(1,npart1_tmp+1:npart_tmp) ) ) > ABS(center2) - &
                                           0.95*bnsb% get_radius2_x_comp() &
          .OR. &
          MAXVAL( ABS( tmp_pos2(1,npart1_tmp+1:npart_tmp) ) ) < ABS(center2) + &
                                           0.95*bnsb% get_radius2_x_opp() &
          .OR. &
          ABS( MINVAL( tmp_pos2(2,npart1_tmp+1:npart_tmp) ) ) < &
                      0.95*bnsb% get_radius2_y() &
          .OR. &
          ABS( MAXVAL( tmp_pos2(2,npart1_tmp+1:npart_tmp) ) ) < &
                      0.95*bnsb% get_radius2_y() &
          .OR. &
          ABS( MINVAL( tmp_pos2(3,npart1_tmp+1:npart_tmp) ) ) < &
                      0.95*bnsb% get_radius2_z() &
          .OR. &
          ABS( MAXVAL( tmp_pos2(3,npart1_tmp+1:npart_tmp) ) ) < &
                      0.95*bnsb% get_radius2_z() &

      )THEN

        PRINT *, "** ERROR! The positions of the particles on star 2, ", &
                 "read from file " &
                 // TRIM(parts_pos_namefile), " are not compatible with the ", &
                 "binary system read from the |lorene| binary file. Stopping..."
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

      parts% npart1= 0
      DO itr= 1, npart1_tmp, 1

        IF( tmp_pos(columns(3),itr) > 0 )THEN

          parts% npart1= parts% npart1 + 1
          tmp_pos2(1,parts% npart1)= tmp_pos(columns(1),itr)
          tmp_pos2(2,parts% npart1)= tmp_pos(columns(2),itr)
          tmp_pos2(3,parts% npart1)= tmp_pos(columns(3),itr)
          IF( read_nu ) tmp_pos2(4,parts% npart1)= tmp_pos(column_nu,itr)

        ENDIF

      ENDDO

      DO itr= 1, parts% npart1, 1

        tmp_pos2(1,parts% npart1+itr)=   tmp_pos2(1,itr)
        tmp_pos2(2,parts% npart1+itr)=   tmp_pos2(2,itr)
        tmp_pos2(3,parts% npart1+itr)= - tmp_pos2(3,itr)
        IF( read_nu ) tmp_pos2(4,parts% npart1+itr)= tmp_pos2(4,itr)

      ENDDO

      parts% npart1= 2*parts% npart1

      parts% npart2= 0
      DO itr= npart1_tmp + 1, npart_tmp, 1

        IF( tmp_pos(columns(3),itr) > 0 )THEN

          parts% npart2= parts% npart2 + 1
          tmp_pos2(1,parts% npart1+parts% npart2)= tmp_pos(columns(1),itr)
          tmp_pos2(2,parts% npart1+parts% npart2)= tmp_pos(columns(2),itr)
          tmp_pos2(3,parts% npart1+parts% npart2)= tmp_pos(columns(3),itr)
          IF( read_nu ) tmp_pos2(4,parts% npart1+parts% npart2)= tmp_pos(column_nu,itr)

        ENDIF

      ENDDO

      DO itr= 1, parts% npart2, 1

        tmp_pos2(1,parts% npart1+parts% npart2+itr)=   &
                                            tmp_pos2(1,parts% npart1+itr)
        tmp_pos2(2,parts% npart1+parts% npart2+itr)=   &
                                            tmp_pos2(2,parts% npart1+itr)
        tmp_pos2(3,parts% npart1+parts% npart2+itr)= &
                                          - tmp_pos2(3,parts% npart1+itr)
        IF( read_nu ) tmp_pos2(4,parts% npart1+parts% npart2+itr)= &
                                            tmp_pos2(4,parts% npart1+itr)

      ENDDO

      parts% npart2= 2*parts% npart2
      parts% npart = parts% npart1 + parts% npart2

      !PRINT *, tmp_pos(:,1)
      ! Allocating the memory for the array pos( 3, npart )
      IF(.NOT.ALLOCATED( parts% pos ))THEN
        ALLOCATE( parts% pos( 3, parts% npart ), STAT= ios, &
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
      IF( read_nu .AND. .NOT.ALLOCATED( parts% nu ))THEN
        ALLOCATE( parts% nu( parts% npart ), STAT= ios, &
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
      IF( read_nu .AND. .NOT.ALLOCATED( parts% pmass ))THEN
        ALLOCATE( parts% pmass( parts% npart ), STAT= ios, &
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
  !    parts% pos(:,1:parts% npart1/2)= &
  !                                  tmp_pos2(:,1:parts% npart1/2)
  !
  !    ! Particles with z < 0 for star 1
  !    parts% pos(1:2,parts% npart1/2+1:parts% npart1)= &
  !                                  tmp_pos2(1:2,1:parts% npart1/2)
  !
  !    parts% pos(3,parts% npart1/2+1:parts% npart1)= &
  !                                - tmp_pos2(3,1:parts% npart1/2)
  !
  !    ! Particles with z > 0 for star 2
  !    parts% pos(:,parts% npart1+1: &
  !                     parts% npart1+parts% npart2/2)= &
  !    tmp_pos2(:,parts% npart1/2+1:parts% npart1/2+parts% npart2/2)
  !
  !    ! Particles with z < 0 for star 2
  !    parts% pos(1:2,parts% npart1+parts% npart2/2+1: &
  !                     parts% npart)= &
  !    tmp_pos2(1:2,parts% npart1/2+1:parts% npart1/2+parts% npart2/2)
  !
  !    parts% pos(3,parts% npart1+parts% npart2/2+1: &
  !                     parts% npart)= &
  !    tmp_pos2(3,parts% npart1/2+1:parts% npart1/2+parts% npart2/2)

      parts% pos= tmp_pos2(1:3,1:parts% npart)
      IF( read_nu ) parts% nu= tmp_pos2(4,1:parts% npart)

      PRINT *, " * Particle positions read. Number of particles=", &
               parts% npart
      PRINT *
      PRINT *, " * Number of particles on NS 1=", parts% npart1
      PRINT *, " * Number of particles on NS 2=", parts% npart2
      PRINT *

      !
      !-- Computing volume per particle
      !
      IF(.NOT.ALLOCATED( parts% pvol ))THEN
        ALLOCATE( parts% pvol( parts% npart ), STAT= ios, &
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
      DO itr= 1, parts% npart - 1, 1

        pvol_tmp= pvol_tmp + ABS( parts% pos(3,itr + 1) &
                                - parts% pos(3,itr) )

      ENDDO
      pvol_tmp= pvol_tmp/( parts% npart - 1 )

      parts% pvol= 2.0D0*pvol_tmp**3.0D0

      IF( parts% mass1 > parts% mass2 )THEN

        ! mass_ratio < 1
        parts% mass_ratio= parts% mass2/parts% mass1

      ELSE

        ! mass_ratio < 1
        parts% mass_ratio= parts% mass1/parts% mass2

      ENDIF

      parts% pmass= parts% nu * amu

      !STOP

    CASE(1)

      PRINT *, " * Placing particles on one lattice around the stars."
      PRINT *

      !
      !-- Determine boundaries of the single lattice around the stars (Msun_geo)
      !
      xmin=   bnsb% get_center1_x() - &
                                stretch*MAX( bnsb% get_radius1_x_comp(), &
                                             bnsb% get_radius1_x_opp() )
      xmax=   bnsb% get_center2_x() + &
                                stretch*MAX( bnsb% get_radius2_x_comp(), &
                                             bnsb% get_radius2_x_opp() )
      ymin= - stretch*bnsb% get_radius1_y()
      ymax=   stretch*bnsb% get_radius2_y()
      zmin= - stretch*bnsb% get_radius1_z()
      zmax=   stretch*bnsb% get_radius2_z()

      ! Place particles, and time the process
    !  CALL parts% placer_timer% start_timer()
    !  CALL parts% place_particles_lattice( xmin, xmax, ymin, &
    !                                           ymax, zmin, zmax, &
    !                                           nx, ny, nz, &
    !                                           thres, bnsb )
    !  CALL parts% placer_timer% stop_timer()

    DO itr= 1, bnsb% n_matter, 1

      !xmin= bnsb% centers(1,jx) - stretch*bnsb% larger_radii_x(1)
      !xmax= bnsb% centers(1,jx) - stretch*bnsb% larger_radii_x(1)
      !ymin= bnsb% centers(1,jy) - stretch*bnsb% larger_radii_y(1)
      !ymax= bnsb% centers(1,jy) - stretch*bnsb% larger_radii_y(1)
      !zmin= bnsb% centers(1,jz) - stretch*bnsb% larger_radii_z(1)
      !zmax= bnsb% centers(1,jz) - stretch*bnsb% larger_radii_z(1)

    ENDDO

    CALL parts% placer_timer% start_timer()
    !DO itr= 1, bnsb% n_matter, 1
    !
    !  CALL parts% place_particles_lattice(
    !      bnsb% get_center1_x() - stretch*MAX( bnsb% get_radius1_x_comp(), &
    !                                           bnsb% get_radius1_x_opp() )
    !      bnsb% get_center1_x() - stretch*MAX( bnsb% get_radius1_x_comp(), &
    !                                           bnsb% get_radius1_x_opp() )
    !  )
    !
    !ENDDO
    CALL parts% placer_timer% stop_timer()

    CASE(2)

      PRINT *, " * Placing particles on two lattices, " &
               // "one around each star."
      PRINT *

      !parts% nx= nx
      !parts% ny= ny
      !parts% nz= nz

      !
      !-- Determine boundaries of the two lattices around the stars (Msun_geo)
      !
      xmin1=   bnsb% get_center1_x() - &
                                stretch*MAX( bnsb% get_radius1_x_comp(), &
                                             bnsb% get_radius1_x_opp() )
      xmax1=   bnsb% get_center1_x() + &
                                stretch*MAX( bnsb% get_radius1_x_comp(), &
                                             bnsb% get_radius1_x_opp() )
      ymin1= - stretch*bnsb% get_radius1_y()
      ymax1=   stretch*bnsb% get_radius1_y()
      zmin1= - stretch*bnsb% get_radius1_z()
      zmax1=   stretch*bnsb% get_radius1_z()

      xmin2=   bnsb% get_center2_x() - &
                                stretch*MAX( bnsb% get_radius2_x_comp(), &
                                             bnsb% get_radius2_x_opp() )
      xmax2=   bnsb% get_center2_x() + &
                                stretch*MAX( bnsb% get_radius2_x_comp(), &
                                             bnsb% get_radius2_x_opp() )
      ymin2= - stretch*bnsb% get_radius2_y()
      ymax2=   stretch*bnsb% get_radius2_y()
      zmin2= - stretch*bnsb% get_radius2_z()
      zmax2=   stretch*bnsb% get_radius2_z()

      ! Place particles, and time the process
      CALL parts% placer_timer% start_timer()
      CALL parts% place_particles_lattices( xmin1, xmax1, ymin1, &
                                                  ymax1, zmin1, zmax1, &
                                                  xmin2, xmax2, ymin2, &
                                                  ymax2, zmin2, zmax2, &
                                                  nx, ny, nz, &
                                                  thres, bnsb )
      CALL parts% placer_timer% stop_timer()

    CASE(3)

      PRINT *, "** Placing equal-mass particles on spherical surfaces, " &
               // "taking into account the mass profile of the stars."
      PRINT *

      ! Here the particle mass is computed using the radial mass profile
      ! of the star, so nu should not be redistributed to achieve a given
      ! particle mass ratio
      IF( parts% redistribute_nu .EQV. .TRUE. )THEN
          parts% redistribute_nu= .FALSE.
      ENDIF

      ! TODO: Change back the inequality from < to > if the IF statement!
      !       Changed for debugging purposes
      first_star_more_massive: IF( parts% mass1 > parts% mass2 )THEN

        filename_mass_profile= "spherical_surfaces_mass_profile2.dat"
        filename_shells_radii= "spherical_surfaces_radii2.dat"
        filename_shells_pos  = "spherical_surfaces_pos2.dat"

        ! Place particles, and time the process
        CALL parts% placer_timer% start_timer()
        CALL parts% place_particles_spherical_surfaces( parts% mass2, &
                                                    radius2, center2, &
                                                    central_density2, &
                                                    npart_approx, &
                                                    parts% npart2, &
                                                    pos2, pvol2, pmass2, &
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

        ! mass_ratio < 1
        parts% mass_ratio= parts% mass2/parts% mass1

        equal_masses: IF( parts% mass_ratio >= 0.995 .AND. &
            parts% mass_ratio <= 1.005 .AND. reflect_particles_x )THEN

          IF(.NOT.ALLOCATED( pos1 ))THEN
            ALLOCATE( pos1( 3, parts% npart2 ), &
                      STAT= ios, ERRMSG= err_msg )
            IF( ios > 0 )THEN
               PRINT *, "...allocation error for array pos in SUBROUTINE" &
                        // "place_particles_. ", &
                        "The error message is", err_msg
               STOP
            ENDIF
          ENDIF
          IF(.NOT.ALLOCATED( pvol1 ))THEN
            ALLOCATE( pvol1( parts% npart2 ), &
                      STAT= ios, ERRMSG= err_msg )
            IF( ios > 0 )THEN
               PRINT *, "...allocation error for array pvol in SUBROUTINE" &
                        // "place_particles_. ", &
                        "The error message is", err_msg
               STOP
            ENDIF
          ENDIF
          IF(.NOT.ALLOCATED( pmass1 ))THEN
            ALLOCATE( pmass1( parts% npart2 ), &
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
          parts% npart1= parts% npart2

        ELSE

          IF( parts% mass_ratio >= 0.95 .AND. &
              parts% mass_ratio <= 1.05 )THEN
            npart2_approx= npart_approx/parts% mass_ratio
          ELSE
            npart2_approx= parts% npart1/parts% mass_ratio
          ENDIF

          filename_mass_profile= "spherical_surfaces_mass_profile1.dat"
          filename_shells_radii= "spherical_surfaces_radii1.dat"
          filename_shells_pos  = "spherical_surfaces_pos1.dat"

          CALL parts% place_particles_spherical_surfaces( parts% mass1,&
                                                radius1, center1, &
                                                central_density1, &
                                                npart2_approx, &
                                                parts% npart1, &
                                                pos1, pvol1, pmass1, &
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

        ENDIF equal_masses

        CALL parts% placer_timer% stop_timer()

        parts% npart= parts% npart1 + parts% npart2

      ELSE

        filename_mass_profile= "spherical_surfaces_mass_profile1.dat"
        filename_shells_radii= "spherical_surfaces_radii1.dat"
        filename_shells_pos  = "spherical_surfaces_pos1.dat"

        ! Place particles, and time the process
        CALL parts% placer_timer% start_timer()

        !DO

        CALL parts% place_particles_spherical_surfaces( parts% mass1, &
                                              radius1, center1, &
                                              central_density1, &
                                              npart_approx, &
                                              parts% npart1, &
                                              pos1, pvol1, pmass1, &
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

        IF( debug ) PRINT *, "30"

        ! mass_ratio < 1
        parts% mass_ratio= parts% mass1/parts% mass2

        IF( debug ) PRINT *, "31"

        equal_masses2: IF( parts% mass_ratio >= 0.995 .AND. &
            parts% mass_ratio <= 1.005 .AND. reflect_particles_x )THEN

          IF(.NOT.ALLOCATED( pos2 ))THEN
            ALLOCATE( pos2( 3, parts% npart1 ), &
                      STAT= ios, ERRMSG= err_msg )
            IF( ios > 0 )THEN
               PRINT *, "...allocation error for array pos in SUBROUTINE" &
                        // "place_particles_. ", &
                        "The error message is", err_msg
               STOP
            ENDIF
          ENDIF
          IF(.NOT.ALLOCATED( pvol2 ))THEN
            ALLOCATE( pvol2( parts% npart1 ), &
                      STAT= ios, ERRMSG= err_msg )
            IF( ios > 0 )THEN
               PRINT *, "...allocation error for array pvol in SUBROUTINE" &
                        // "place_particles_. ", &
                        "The error message is", err_msg
               STOP
            ENDIF
          ENDIF
          IF(.NOT.ALLOCATED( pmass2 ))THEN
            ALLOCATE( pmass2( parts% npart1 ), &
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
          parts% npart2= parts% npart1

        ELSE

          IF( parts% mass_ratio >= 0.95 .AND. &
              parts% mass_ratio <= 1.05 )THEN
            npart2_approx= npart_approx/parts% mass_ratio
          ELSE
            npart2_approx= parts% npart1/parts% mass_ratio
          ENDIF

          filename_mass_profile= "spherical_surfaces_mass_profile2.dat"
          filename_shells_radii= "spherical_surfaces_radii2.dat"
          filename_shells_pos  = "spherical_surfaces_pos2.dat"

          IF( debug ) PRINT *, "32"

          CALL parts% place_particles_spherical_surfaces( parts% mass2,&
                                                radius2, center2, &
                                                central_density2, &
                                                npart2_approx, &
                                                parts% npart2, &
                                                pos2, pvol2, pmass2, &
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

        ENDIF equal_masses2

        CALL parts% placer_timer% stop_timer()

        parts% npart= parts% npart1 + parts% npart2

      ENDIF first_star_more_massive

      !
      !-- Assign TYPE member variables
      !

      IF(.NOT.ALLOCATED( parts% pos ))THEN
        ALLOCATE( parts% pos( 3, parts% npart ), &
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
      parts% pos( :, 1:parts% npart1 )= pos1
      parts% pos( :, parts% npart1 + 1:parts% npart )= pos2

      IF(.NOT.ALLOCATED( parts% pvol ))THEN
        ALLOCATE( parts% pvol( parts% npart ), &
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
      parts% pvol( 1:parts% npart1 )= pvol1
      parts% pvol( parts% npart1 + 1:parts% npart )= pvol2

      IF(.NOT.ALLOCATED( parts% pmass ))THEN
        ALLOCATE( parts% pmass( parts% npart ), &
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
      parts% pmass( 1:parts% npart1 )= pmass1
      parts% pmass( parts% npart1 + 1:parts% npart )= pmass2

      PRINT *, " * Particles placed. Number of particles=", parts% npart
      PRINT *, " * Number of particles on NS 1=", parts% npart1
      PRINT *, " * Number of particles on NS 2=", parts% npart2
      PRINT *
      !STOP

    CASE DEFAULT

      PRINT *, "** There is no implemented particle placer " &
               // "corresponding to the number", dist
      PRINT *, " * Stopping..."
      STOP

    END SELECT choose_particle_placer

    !----------------------------------------------!
    !--  At this point,the particles are placed  --!
    !----------------------------------------------!

    ! Reshape the arrays pos and pvol by deleting the unnecessary elements
    parts% pos = parts% pos( :, 1:parts% npart )
    parts% pvol= parts% pvol( 1:parts% npart )

    ! Check that there aren't particles with the same coordinates
    CALL parts% same_particle_timer% start_timer()
    CALL check_particle_positions( parts% npart, parts% pos )
    CALL parts% same_particle_timer% stop_timer()

    IF( apm_iterate1 )THEN

      PRINT *
      PRINT *, "** Placing particles on star 1 using the APM..."
      PRINT *

      IF(.NOT.ALLOCATED( parts% h ))THEN
        ALLOCATE( parts% h( parts% npart ), STAT= ios, &
                  ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array h in SUBROUTINE ", &
                    "construct_particles_bnsbase. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF

      IF(.NOT.ALLOCATED( parts% nu ))THEN
        ALLOCATE( parts% nu( parts% npart ), STAT= ios, &
                  ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array nu in SUBROUTINE ", &
                    "construct_particles_bnsbase. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF

      filename_apm_pos_id = "apm_pos_id1.dat"
      filename_apm_pos    = "apm_pos1.dat"
      filename_apm_results= "apm_results1.dat"

      ! Star 1
      CALL parts% apm1_timer% start_timer()
      CALL parts% perform_apm( &
                  import_density, get_nstar_p, &
                  parts% pos(:,1:parts% npart1), &
                  parts% pvol(1:parts% npart1), &
                  parts% h(1:parts% npart1), &
                  parts% nu(1:parts% npart1), &
                  center1, com1, parts% mass1, &
                  bnsb% get_radius1_x_comp(), &
                  bnsb% get_radius1_x_opp(), &
                  bnsb% get_radius1_y(), &
                  bnsb% get_radius1_z(), &
                  apm_max_it, max_inc, mass_it, parts% correct_nu, &
                  nuratio_thres, nuratio_des, nx_gh, ny_gh, nz_gh, &
                  filename_apm_pos_id, filename_apm_pos, filename_apm_results, &
                  check_negative_hydro )
      CALL parts% apm1_timer% stop_timer()

      PRINT *, "** Particles placed on star 1 according to the APM."
      PRINT *

      IF( parts% mass_ratio >= 0.995 .AND. &
          parts% mass_ratio <= 1.005 .AND. reflect_particles_x )THEN

        parts% pos(1,parts% npart1+1:parts% npart)= &
                                  - parts% pos(1,1:parts% npart1)
        parts% pos(2,parts% npart1+1:parts% npart)= &
                                    parts% pos(2,1:parts% npart1)
        parts% pos(3,parts% npart1+1:parts% npart)= &
                                    parts% pos(3,1:parts% npart1)

        parts% nu(parts% npart1+1:parts% npart)= &
                                    parts% nu(1:parts% npart1)

        parts% h(parts% npart1+1:parts% npart)= &
                                    parts% h(1:parts% npart1)

        parts% npart2= parts% npart1
        parts% npart= parts% npart1 + parts% npart1

      !ELSEIF( ( parts% mass_ratio <= 0.995 .OR. &
      !        parts% mass_ratio >= 1.005 ) .AND. reflect_particles_x )THEN
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
    IF( apm_iterate2 .AND. .NOT.(parts% mass_ratio >= 0.995 .AND. &
        parts% mass_ratio <= 1.005 .AND. reflect_particles_x) )THEN

      PRINT *
      PRINT *, "** Placing particles on star 2 using the APM..."
      PRINT *

      IF(.NOT.ALLOCATED( parts% h ))THEN
        ALLOCATE( parts% h( parts% npart ), STAT= ios, &
                  ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array h in SUBROUTINE ", &
                    "construct_particles_bnsbase. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF

      IF(.NOT.ALLOCATED( parts% nu ))THEN
        ALLOCATE( parts% nu( parts% npart ), STAT= ios, &
                  ERRMSG= err_msg )
        IF( ios > 0 )THEN
           PRINT *, "...allocation error for array nu in SUBROUTINE ", &
                    "construct_particles_bnsbase. The error message is",&
                    err_msg
           STOP
        ENDIF
      ENDIF

      filename_apm_pos_id = "apm_pos_id2.dat"
      filename_apm_pos    = "apm_pos2.dat"
      filename_apm_results= "apm_results2.dat"

      ! Star 2
      CALL parts% apm2_timer% start_timer()
      CALL parts% perform_apm( &
                import_density, get_nstar_p, &
                parts% pos(:,parts% npart1+1:parts% npart), &
                parts% pvol(parts% npart1+1:parts% npart), &
                parts% h(parts% npart1+1:parts% npart), &
                parts% nu(parts% npart1+1:parts% npart), &
                center2, com2, parts% mass2, &
                bnsb% get_radius2_x_comp(), &
                bnsb% get_radius2_x_opp(), &
                bnsb% get_radius2_y(), &
                bnsb% get_radius2_z(), &
                apm_max_it, max_inc, mass_it, parts% correct_nu, &
                nuratio_thres, nuratio_des, nx_gh, ny_gh, nz_gh, &
                filename_apm_pos_id, filename_apm_pos, filename_apm_results, &
                check_negative_hydro )
      CALL parts% apm2_timer% stop_timer()

      PRINT *, "** Particles placed on star 2 according to the APM."
      PRINT *

    ENDIF

    ! Allocate needed memory
    CALL allocate_lorene_id_parts_memory( parts )

    ! flag that particles are 'alive'
    ALLOCATE( alive( parts% npart ) )
    alive( 1:parts% npart )= 1

    IF( debug ) PRINT *, "33"

    !
    !-- Import the needed |lorene| ID on the particles, and time the process
    !
    PRINT *, "** Importing the |lorene| ID on the particles..."

    CALL parts% importer_timer% start_timer()
    CALL bnsb% read_id_particles( parts% npart, &
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

    IF(.NOT.ALLOCATED( parts% baryon_density_index ))THEN
      ALLOCATE( parts% baryon_density_index( parts% npart ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array baryon_density_index in " &
                  // "SUBROUTINE construct_particles_bnsbase. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array pos in SUBROUTINE" &
      !                // "place_particles_3D_lattice." )
    ENDIF

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

    IF( parts% redistribute_nu )THEN

      ! Index particles on star 1 in increasing order of nu

      CALL indexx( parts% npart1, &
                   parts% baryon_density_parts( 1 : parts% npart1 ), &
                   parts% baryon_density_index( 1 : parts% npart1 ) )

      ! Index particles on star 2 in increasing order of nu

      CALL indexx( parts% npart2, &
                   parts% baryon_density_parts( parts% npart1 + 1 : &
                                                    parts% npart ), &
                   parts% baryon_density_index( parts% npart1 + 1 : &
                                                    parts% npart ) )

      ! Shift indices on star 2 by npart1 since all the arrays store
      ! the quantities on star 1 first, and then on star 2

      parts% baryon_density_index( parts% npart1 + 1 : &
                                       parts% npart )= &
                     parts% npart1 + &
                     parts% baryon_density_index( parts% npart1 + 1 : &
                                                      parts% npart )

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

      density= bnsb% read_mass_density( x, y, z )

    END FUNCTION import_density


    SUBROUTINE import_id( x, y, z, &
                          g_xx, &
                          baryon_density, &
                          gamma_euler )

      IMPLICIT NONE

      DOUBLE PRECISION,   INTENT( IN )    :: x
      DOUBLE PRECISION,   INTENT( IN )    :: y
      DOUBLE PRECISION,   INTENT( IN)     :: z
      DOUBLE PRECISION, INTENT( IN OUT ):: g_xx
      DOUBLE PRECISION, INTENT( IN OUT ):: baryon_density
      DOUBLE PRECISION, INTENT( IN OUT ):: gamma_euler

      CALL bnsb% read_id_mass_b( x, y, z, &
                               g_xx, &
                               baryon_density, &
                               gamma_euler  )

    END SUBROUTINE import_id


    SUBROUTINE integrate_mass_density( center, radius, &
                                  central_density, &
                                  dr, dth, dphi, &
                                  mass, mass_profile, &
                                  mass_profile_idx )

      IMPLICIT NONE

      !& Array to store the indices for array mass_profile, sorted so that
      !  mass_profile[mass_profile_idx] is in increasing order
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT( IN OUT ):: mass_profile_idx
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
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT( IN OUT ):: &
                                       mass_profile

      CALL bnsb% integrate_baryon_mass_density( center, radius, &
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

      answer= bnsb% test_position( x, y, z )

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

      CALL bnsb% read_id_particles( npart_real, x, y, z, &
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

      USE constants, ONLY: Msun_geo, km2m, amu, g2kg
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

        nstar_p(a)= sq_g*Theta_a*baryon_density(a)*((Msun_geo*km2m)**3) &
                    /(amu*g2kg)

        IF( ISNAN( nstar_p( a ) ) )THEN
          PRINT *, "** ERROR! nstar_p(", a, ") is a NaN!", &
                   " Stopping.."
          PRINT *
          STOP
        ENDIF
        IF( nstar_p( a ) == 0 )THEN
          PRINT *, "** ERROR! nstar_p(", a, ")= 0 on a real particle!"
          !PRINT *, " * Particle position: x=", all_pos(1,a), &
          !         ", y=", all_pos(2,a), ", z=", all_pos(3,a)
          PRINT *, "   sq_g=", sq_g
          PRINT *, "   Theta_a=", Theta_a
          PRINT *, "   baryon_density(", a, ")=", baryon_density(a)
          PRINT *, " * Stopping.."
          PRINT *
          STOP
        ENDIF

      ENDDO
      !$OMP END PARALLEL DO

    END SUBROUTINE compute_nstar_p


  END PROCEDURE construct_particles_bnsbase


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
