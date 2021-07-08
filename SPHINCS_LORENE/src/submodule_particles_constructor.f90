! File:         submodule_particles_constructor.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (particles_id) particles_constructor

  !************************************************
  !                                               *
  ! Implementation of the constructor of TYPE     *
  ! particles, and all the PROCEDURES it calls.   *
  !                                               *
  ! FT 16.10.2020                                 *
  !                                               *
  !************************************************


  IMPLICIT NONE


  CONTAINS


  !MODULE PROCEDURE construct_particles_empty
  !
  !    !************************************************
  !    !                                               *
  !    ! The constructor of an empty particle object.  *
  !    !                                               *
  !    ! FT 02.11.2020                                 *
  !    !                                               *
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

    !************************************************
    !                                               *
    ! The constructor performs all the tasks needed *
    ! to set up the particle distribution with the  *
    ! LORENE ID on it. It calls all the PROCEDURES  *
    ! that rely on an object of TYPE bns.           *
    !                                               *
    ! FT 17.10.2020                                 *
    !                                               *
    !************************************************

    !USE NaNChecker, ONLY: Check_Array_for_NAN
    USE constants,      ONLY: Msun_geo, km2m
    USE NR,             ONLY: indexx
    USE kernel_table,   ONLY: ktable
    USE input_output,   ONLY: read_options
    USE units,          ONLY: umass, set_units
    USE options,        ONLY: ikernel, ndes
    USE alive_flag,     ONLY: alive

    IMPLICIT NONE

    ! The variable counter counts how many times the PROCEDURE
    ! construct_particles is called
    INTEGER, SAVE:: counter= 1
    INTEGER:: nx, ny, nz, min_y_index, min_z_index, cntr1, cntr2, &
              npart_approx, npart2_approx, itr_1, itr_2, max_steps, itr2, &
              nlines, header_lines, n_cols, npart_tmp, npart1_tmp, npart2_tmp, &
              nx_gh, ny_gh, nz_gh
    ! Maximum length for strings, and for the number of imported binaries
    INTEGER, PARAMETER:: max_length= 50
    INTEGER, DIMENSION( : ), ALLOCATABLE:: x_sort, y_sort, z_sort
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
    DOUBLE PRECISION:: min_abs_y, min_abs_z
    DOUBLE PRECISION:: upper_bound, lower_bound, upper_factor, lower_factor, &
                       last_r
    DOUBLE PRECISION:: pvol_tmp
    DOUBLE PRECISION, DIMENSION( :, : ), ALLOCATABLE:: abs_pos
    DOUBLE PRECISION, DIMENSION( :, : ), ALLOCATABLE:: tmp_pos
    DOUBLE PRECISION, DIMENSION( :, : ), ALLOCATABLE:: tmp_pos2
    DOUBLE PRECISION, DIMENSION( :, : ), ALLOCATABLE:: pos1, pos2
    DOUBLE PRECISION, DIMENSION( : ),    ALLOCATABLE:: pvol1, pvol2, &
                                                       pmass1, pmass2
    DOUBLE PRECISION:: nuratio_thres

    CHARACTER( LEN= : ), ALLOCATABLE:: namefile
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
              randomize_r, apm_iterate, mass_it, find_npart, read_nu, &
              reflect_particles_x
    LOGICAL, DIMENSION( : ), ALLOCATABLE:: negative_hydro

    NAMELIST /bns_particles/ &
              parts_pos_path, parts_pos, columns, header_lines, n_cols, &
              read_nu, column_nu, &
              stretch, &
              nx, ny, nz, &
              use_thres, thres, nu_ratio, redistribute_nu, correct_nu, &
              compose_eos, compose_path, compose_filename, &
              npart_approx, last_r, upper_bound, lower_bound, &
              upper_factor, lower_factor, max_steps, &
              randomize_phi, randomize_theta, randomize_r, find_npart, &
              apm_iterate, apm_max_it, max_inc, mass_it, nuratio_thres, &
              reflect_particles_x, nx_gh, ny_gh, nz_gh

    !
    !-- Initialize the timers
    !
    parts_obj% placer_timer      = timer( "placer_timer" )
    parts_obj% apm1_timer        = timer( "apm_star1_timer" )
    parts_obj% apm2_timer        = timer( "apm_star2_timer" )
    parts_obj% importer_timer    = timer( "importer_timer" )
    parts_obj% sph_computer_timer= timer( "sph_computer_timer" )

    ! Declare this object as non-empty (experimental)
    parts_obj% empty_object= .FALSE.
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
    parts_obj% randomize_phi  = randomize_phi
    parts_obj% randomize_theta= randomize_theta
    parts_obj% randomize_r    = randomize_r

    !
    !-- Read the parameters of the particle distributions
    !
    parts_obj% lorene_bns_id_parfile= 'lorene_bns_id_particles.par'

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

    parts_obj% use_thres= use_thres
    parts_obj% correct_nu= correct_nu
    parts_obj% compose_eos= compose_eos
    parts_obj% compose_path= compose_path
    parts_obj% compose_filename= compose_filename
    parts_obj% redistribute_nu= redistribute_nu
    parts_obj% nu_ratio= nu_ratio
    ! APM parameters
    parts_obj% apm_iterate   = apm_iterate
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

    ! setup unit system
    CALL set_units('NSM')
    CALL read_options       ! TODO: set units and read options only once in the constructor

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

        IF( tmp_pos(columns(1),itr) < 0 )THEN

          npart1_tmp= npart1_tmp + 1
          tmp_pos2(1,npart1_tmp)= tmp_pos(columns(1),itr)
          tmp_pos2(2,npart1_tmp)= tmp_pos(columns(2),itr)
          tmp_pos2(3,npart1_tmp)= tmp_pos(columns(3),itr)
          IF( read_nu ) tmp_pos2(4,npart1_tmp)= tmp_pos(column_nu,itr)

        ENDIF

      ENDDO

      npart2_tmp= 0
      DO itr= 1, npart_tmp, 1

        IF( tmp_pos(columns(1),itr) > 0 )THEN

          npart2_tmp= npart2_tmp + 1
          tmp_pos2(1,npart1_tmp+npart2_tmp)= tmp_pos(columns(1),itr)
          tmp_pos2(2,npart1_tmp+npart2_tmp)= tmp_pos(columns(2),itr)
          tmp_pos2(3,npart1_tmp+npart2_tmp)= tmp_pos(columns(3),itr)
          IF( read_nu ) tmp_pos2(4,npart1_tmp+npart2_tmp)= tmp_pos(column_nu,itr)

        ENDIF

      ENDDO

      IF( npart1_tmp + npart2_tmp /= npart_tmp )THEN
        PRINT *, "** ERROR! parts_obj% npart1 + parts_obj% npart2 /= npart_tmp"
        PRINT *
        PRINT *, "   parts_obj% npart1= ", npart1_tmp
        PRINT *, "   parts_obj% npart2= ", npart2_tmp
        PRINT *, "   parts_obj% npart1 + parts_obj% npart2= ", &
                 npart1_tmp + npart2_tmp
        PRINT *, "   npart_tmp= ", npart_tmp
        PRINT *
        STOP
      ENDIF

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

      pvol_tmp= 0
      DO itr= 1, parts_obj% npart - 1, 1

        pvol_tmp= pvol_tmp + ABS( parts_obj% pos(3,itr + 1) &
                                - parts_obj% pos(3,itr) )

      ENDDO
      pvol_tmp= pvol_tmp/( parts_obj% npart - 1)

      parts_obj% pvol= pvol_tmp

      IF( parts_obj% mass1 > parts_obj% mass2 )THEN

        ! mass_ratio < 1
        parts_obj% mass_ratio= parts_obj% mass2/parts_obj% mass1

      ELSE

        ! mass_ratio < 1
        parts_obj% mass_ratio= parts_obj% mass1/parts_obj% mass2

      ENDIF

      !STOP

    CASE(1)

      PRINT *, " * Placing particles on one lattice around the stars."
      PRINT *

      parts_obj% nx= nx
      parts_obj% ny= ny
      parts_obj% nz= nz

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
      CALL parts_obj% place_particles_3dlattice( xmin, xmax, ymin, &
                                                 ymax, zmin, zmax, thres, &
                                                 bns_obj )
      CALL parts_obj% placer_timer% stop_timer()

    CASE(2)

      PRINT *, " * Placing particles on two lattices, " &
               // "one around each star."
      PRINT *

      parts_obj% nx= nx
      parts_obj% ny= ny
      parts_obj% nz= nz

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
      CALL parts_obj% place_particles_3dlattices( xmin1, xmax1, ymin1, &
                                                  ymax1, zmin1, zmax1, &
                                                  xmin2, xmax2, ymin2, &
                                                  ymax2, zmin2, zmax2, thres, &
                                                  bns_obj )
      CALL parts_obj% placer_timer% stop_timer()

    CASE(3)

      PRINT *, "** Placing equal-mass particles on spherical shells, " &
               // "taking into account the mass profile of the stars."
      PRINT *

      IF( parts_obj% redistribute_nu == .TRUE. )THEN
        parts_obj% redistribute_nu= .FALSE.
      ENDIF

      IF( parts_obj% mass1 > parts_obj% mass2 )THEN

        filename_mass_profile= "shells_mass_profile2.dat"
        filename_shells_radii= "shells_radii2.dat"
        filename_shells_pos  = "shells_pos2.dat"

        n_particles_first_shell= 4

        ! Place particles, and time the process
        CALL parts_obj% placer_timer% start_timer()
        CALL parts_obj% place_particles_spherical_shells( parts_obj% mass2, &
                                                    radius2, center2, &
                                                    npart_approx, &
                                                    parts_obj% npart2, &
                                                    pos2, pvol2, pmass2, &
                                                    thres, &
                                                    bns_obj, &
                                                    last_r, &
                                                    upper_bound, lower_bound, &
                                                    upper_factor, lower_factor,&
                                                    max_steps, &
                                                    n_particles_first_shell, &
                                                    find_npart, &
                                                    filename_mass_profile, &
                                                    filename_shells_radii, &
                                                    filename_shells_pos )

        ! mass_ratio < 1
        parts_obj% mass_ratio= parts_obj% mass2/parts_obj% mass1

        IF( parts_obj% mass_ratio >= 0.995 .AND. &
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

          filename_mass_profile= "shells_mass_profile1.dat"
          filename_shells_radii= "shells_radii1.dat"
          filename_shells_pos  = "shells_pos1.dat"

          n_particles_first_shell= n_particles_first_shell/parts_obj% mass_ratio

          CALL parts_obj% place_particles_spherical_shells( parts_obj% mass1, &
                                                radius1, center1, &
                                                npart2_approx, &
                                                parts_obj% npart1, &
                                                pos1, pvol1, pmass1, &
                                                thres, &
                                                bns_obj, &
                                                last_r, &
                                                upper_bound, lower_bound, &
                                                upper_factor, lower_factor,&
                                                max_steps, &
                                                n_particles_first_shell, &
                                                find_npart, &
                                                filename_mass_profile, &
                                                filename_shells_radii, &
                                                filename_shells_pos )

        ENDIF

        CALL parts_obj% placer_timer% stop_timer()

        parts_obj% npart= parts_obj% npart1 + parts_obj% npart2

      ELSE

        filename_mass_profile= "shells_mass_profile1.dat"
        filename_shells_radii= "shells_radii1.dat"
        filename_shells_pos  = "shells_pos1.dat"

        n_particles_first_shell= 4

        ! Place particles, and time the process
        CALL parts_obj% placer_timer% start_timer()

        !DO

        CALL parts_obj% place_particles_spherical_shells( parts_obj% mass1, &
                                              radius1, center1, &
                                              npart_approx, &
                                              parts_obj% npart1, &
                                              pos1, pvol1, pmass1, &
                                              thres, &
                                              bns_obj, &
                                              last_r, &
                                              upper_bound, lower_bound, &
                                              upper_factor, lower_factor,&
                                              max_steps, &
                                              n_particles_first_shell, &
                                              find_npart, &
                                              filename_mass_profile, &
                                              filename_shells_radii, &
                                              filename_shells_pos )

        ! mass_ratio < 1
        parts_obj% mass_ratio= parts_obj% mass1/parts_obj% mass2
        !npart2_approx= npart_approx/parts_obj% mass_ratio
        !npart2_approx= MIN(npart_approx,parts_obj% npart1)/(parts_obj% mass_ratio)

        IF( parts_obj% mass_ratio >= 0.995 .AND. &
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

        !ELSEIF( ( parts_obj% mass_ratio <= 0.995 .OR. &
        !        parts_obj% mass_ratio >= 1.005 ) .AND. reflect_particles_x )THEN
        !
        !  PRINT *, "** ERROR! The two stars are not the same. The particles", &
        !           " on star 1 cannot be reflected with respect to the yz", &
        !           " plane to become the particles on star 2."
        !  PRINT *, "   Please, choose an equal-mass system, or set the", &
        !           " variable reflect_particles_x to .FALSE. in the file", &
        !           " lorene_bns_id_particles ."
        !  PRINT *
        !  STOP

        ELSE

          IF( parts_obj% mass_ratio >= 0.95 .AND. &
              parts_obj% mass_ratio <= 1.05 )THEN
            npart2_approx= npart_approx/parts_obj% mass_ratio
          ELSE
            npart2_approx= parts_obj% npart1/parts_obj% mass_ratio
          ENDIF

          filename_mass_profile= "shells_mass_profile2.dat"
          filename_shells_radii= "shells_radii2.dat"
          filename_shells_pos  = "shells_pos2.dat"

          n_particles_first_shell= 4!n_particles_first_shell/parts_obj% mass_ratio

          CALL parts_obj% place_particles_spherical_shells( parts_obj% mass2, &
                                                radius2, center2, &
                                                npart2_approx, &
                                                parts_obj% npart2, &
                                                pos2, pvol2, pmass2, &
                                                thres, &
                                                bns_obj, &
                                                last_r, &
                                                upper_bound, lower_bound, &
                                                upper_factor, lower_factor,&
                                                max_steps, &
                                                n_particles_first_shell, &
                                                find_npart, &
                                                filename_mass_profile, &
                                                filename_shells_radii, &
                                                filename_shells_pos )

        ENDIF

          !IF( parts_obj% npart1/parts_obj% npart2 >= &
          !    0.9D0*parts_obj% mass_ratio .AND. &
          !    parts_obj% npart1/parts_obj% npart2 <= &
          !    1.1D0*parts_obj% mass_ratio &
          !)THEN
          !  EXIT
          !ENDIF

          !IF( parts_obj% npart1/parts_obj% npart2 )THEN
          !  CALL parts_obj% place_particles_spherical_shells( parts_obj% mass2, &
          !                                          radius2, center2, &
          !                                          npart2_approx, &
          !                                          parts_obj% npart2, &
          !                                          pos2, pvol2, pmass2, &
          !                                          thres, &
          !                                          bns_obj, &
          !                                          last_r, &
          !                                          upper_bound, lower_bound, &
          !                                          upper_factor, lower_factor,&
          !                                          max_steps )
        !ENDIF
        !ENDDO

        CALL parts_obj% placer_timer% stop_timer()

        parts_obj% npart= parts_obj% npart1 + parts_obj% npart2

        !PRINT *
        !PRINT *, "npart1= ", parts_obj% npart1
        !PRINT *, "npart_approx= ", npart_approx
        !PRINT *, "npart2= ", parts_obj% npart2
        !PRINT *, "npart2_approx= ", npart2_approx
        !PRINT *, "npart= ", parts_obj% npart
        !PRINT *
        !STOP

      ENDIF

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
      STOP

    END SELECT choose_particle_placer

    ! Reshape the arrays pos and pvol by deleting the unnecessary elements
    parts_obj% pos = parts_obj% pos( :, 1:parts_obj% npart )
    parts_obj% pvol= parts_obj% pvol( 1:parts_obj% npart )

    IF(.NOT.ALLOCATED( x_sort ))THEN
      ALLOCATE( x_sort( parts_obj% npart ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array x_sort in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( y_sort ))THEN
      ALLOCATE( y_sort( parts_obj% npart ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array y_sort in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( z_sort ))THEN
      ALLOCATE( z_sort( parts_obj% npart ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array z_sort in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    ! Sort x, y, z coordinates of the particles
    CALL indexx( parts_obj% npart, parts_obj% pos( 1, : ), x_sort )
    CALL indexx( parts_obj% npart, parts_obj% pos( 2, : ), y_sort )
    CALL indexx( parts_obj% npart, parts_obj% pos( 3, : ), z_sort )

    ! Check that there aren't particles with the same coordinates
    PRINT *, "** Checking that there are no multiple particles", &
             " at the same position..."
    DO itr= 1, parts_obj% npart - 1, 1
      IF( parts_obj% pos( 1, x_sort(itr) ) == &
                parts_obj% pos( 1, x_sort(itr + 1) ) .AND. &
          parts_obj% pos( 2, x_sort(itr) ) == &
                parts_obj% pos( 2, x_sort(itr + 1) ) .AND. &
          parts_obj% pos( 3, x_sort(itr) ) == &
                parts_obj% pos( 3, x_sort(itr + 1) ) )THEN
        PRINT *, "** ERROR in SUBROUTINE place_particles_3dlattices! ", &
                 "The two particles ", itr, " and", itr - 1, " have the same ",&
                 "coordinates!"
        PRINT *, parts_obj% pos( 1, x_sort(itr) ), parts_obj% pos( 1, x_sort(itr) ), parts_obj% pos( 1, x_sort(itr) )
        PRINT *, parts_obj% pos( 1, x_sort(itr + 1) ), parts_obj% pos( 1, x_sort(itr + 1) ), parts_obj% pos( 1, x_sort(itr + 1) )
        STOP
      ENDIF
    ENDDO
    DO itr= 1, parts_obj% npart - 1, 1
      IF( parts_obj% pos( 1, y_sort(itr) ) == &
                parts_obj% pos( 1, y_sort(itr + 1) ) .AND. &
          parts_obj% pos( 2, y_sort(itr) ) == &
                parts_obj% pos( 2, y_sort(itr + 1) ) .AND. &
          parts_obj% pos( 3, y_sort(itr) ) == &
                parts_obj% pos( 3, y_sort(itr + 1) ) )THEN
        PRINT *, "** ERROR in SUBROUTINE place_particles_3dlattices! ", &
                 "The two particles ", itr, " and", itr - 1, " have the same ",&
                 "coordinates!"
        PRINT *, parts_obj% pos( 1, y_sort(itr) ), parts_obj% pos( 1, y_sort(itr) ), parts_obj% pos( 1, y_sort(itr) )
        PRINT *, parts_obj% pos( 1, y_sort(itr + 1) ), parts_obj% pos( 1, y_sort(itr + 1) ), parts_obj% pos( 1, y_sort(itr + 1) )
        STOP
      ENDIF
    ENDDO
    DO itr= 1, parts_obj% npart - 1, 1
      IF( parts_obj% pos( 1, z_sort(itr) ) == &
                parts_obj% pos( 1, z_sort(itr + 1) ) .AND. &
          parts_obj% pos( 2, z_sort(itr) ) == &
                parts_obj% pos( 2, z_sort(itr + 1) ) .AND. &
          parts_obj% pos( 3, z_sort(itr) ) == &
                parts_obj% pos( 3, z_sort(itr + 1) ) )THEN
        PRINT *, "** ERROR in SUBROUTINE place_particles_3dlattices! ", &
                 "The two particles ", itr, " and", itr - 1, " have the same ",&
                 "coordinates!"
        PRINT *, parts_obj% pos( 1, z_sort(itr) ), parts_obj% pos( 1, z_sort(itr) ), parts_obj% pos( 1, z_sort(itr) )
        PRINT *, parts_obj% pos( 1, z_sort(itr + 1) ), parts_obj% pos( 1, z_sort(itr + 1) ), parts_obj% pos( 1, z_sort(itr + 1) )
        STOP
      ENDIF
    ENDDO

    !DO itr= 1, parts_obj% npart, 1
    !  DO itr2= itr + 1, parts_obj% npart, 1
    !    IF( parts_obj% pos( 1, itr ) == parts_obj% pos( 1, itr2 ) .AND. &
    !        parts_obj% pos( 2, itr ) == parts_obj% pos( 2, itr2 ) .AND. &
    !        parts_obj% pos( 3, itr ) == parts_obj% pos( 3, itr2 ) )THEN
    !      PRINT *, "** ERROR in SUBROUTINE place_particles_3dlattices! ", &
    !               "The two particles ", itr, " and", itr2, " have the same ", &
    !               "coordinates!"
    !      STOP
    !    ENDIF
    !  ENDDO
    !ENDDO

    IF( apm_iterate )THEN

      PRINT *
      PRINT *, " ** Placing particles using the APM..."
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
                  nuratio_thres, nx_gh, ny_gh, nz_gh, &
                  filename_apm_pos_id, filename_apm_pos, filename_apm_results )
      CALL parts_obj% apm1_timer% stop_timer()

      PRINT *, "APM done for star 1"
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

      ELSE

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
                  nuratio_thres, nx_gh, ny_gh, nz_gh, &
                  filename_apm_pos_id, filename_apm_pos, filename_apm_results )
        CALL parts_obj% apm2_timer% stop_timer()

        PRINT *, "APM done for star 2"
        PRINT *

      ENDIF

      PRINT *, " ** Particles placed according to the APM."
      PRINT *

    ENDIF

    ! Allocate needed memory
    CALL allocate_lorene_id_parts_memory( parts_obj )

    ! flag that particles are 'alive'
    ALLOCATE( alive( parts_obj% npart ) )
    alive( 1:parts_obj% npart )= 1

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

    PRINT *, "** Computing typical length scale for the change in pressure", &
             " on the x axis."
    PRINT *

    ALLOCATE( abs_pos( 3, parts_obj% npart ) )

    DO itr = 1, parts_obj% npart, 1
      abs_pos( 1, itr )= ABS( parts_obj% pos( 1, itr ) )
      abs_pos( 2, itr )= ABS( parts_obj% pos( 2, itr ) )
      abs_pos( 3, itr )= ABS( parts_obj% pos( 3, itr ) )
    ENDDO

    min_y_index= 0
    min_abs_y= 1D+20
    DO itr = 1, parts_obj% npart, 1
      IF( ABS( parts_obj% pos( 2, itr ) ) < min_abs_y )THEN
        min_abs_y= ABS( parts_obj% pos( 2, itr ) )
        min_y_index= itr
      ENDIF
    ENDDO

    min_z_index= 0
    min_abs_z= 1D+20
    DO itr = 1, parts_obj% npart, 1
      IF( ABS( parts_obj% pos( 3, itr ) ) < min_abs_z )THEN
        min_abs_z= ABS( parts_obj% pos( 3, itr ) )
        min_z_index= itr
      ENDIF
    ENDDO

    min_abs_z= MINVAL( abs_pos( 3, : ) )

    !PRINT *, "1"

    cntr1= 0
    cntr2= 0
    DO itr = 1, parts_obj% npart, 1
      IF( parts_obj% pos( 3, itr ) == min_abs_z &
          .AND. &
          ABS( ( parts_obj% pos( 2, itr ) - &
                 parts_obj% pos( 2, min_y_index ) )/ &
                 parts_obj% pos( 2, min_y_index ) ) < 1.0D-5 &
      )THEN

        IF( parts_obj% pos( 1, itr ) < 0 )THEN
          cntr1= cntr1 + 1
        ELSEIF( parts_obj% pos( 1, itr ) > 0 )THEN
          cntr2= cntr2 + 1
        ENDIF

      ENDIF
    ENDDO
    !PRINT *, "cntr1= ", cntr1
    !PRINT *, "cntr2= ", cntr2

    ALLOCATE( parts_obj% pos_x1( cntr1 ) )
    ALLOCATE( parts_obj% pos_x2( cntr2 ) )
    ALLOCATE( parts_obj% pressure_parts_x1( cntr1 ) )
    ALLOCATE( parts_obj% pressure_parts_x2( cntr2 ) )
    ALLOCATE( parts_obj% pressure_parts_x_der1( cntr1 - 5 ) )
    ALLOCATE( parts_obj% pressure_parts_x_der2( cntr2 - 5 ) )
    ALLOCATE( parts_obj% pressure_length_scale_x1( cntr1 - 5 ) )
    ALLOCATE( parts_obj% pressure_length_scale_x2( cntr2 - 5 ) )

    !PRINT *, "2"

    itr_1= 0
    itr_2= 0
    DO itr = 1, parts_obj% npart, 1
      IF( parts_obj% pos( 3, itr ) == min_abs_z &
          .AND. &
          ABS( ( parts_obj% pos( 2, itr ) - &
                 parts_obj% pos( 2, min_y_index ) )/ &
                 parts_obj% pos( 2, min_y_index ) ) < 1.0D-5 &
        )THEN

        IF( parts_obj% pos( 1, itr ) < 0 )THEN
          itr_1= itr_1 + 1
          parts_obj% pos_x1( itr_1 )= parts_obj% pos( 1, itr )
          parts_obj% pressure_parts_x1( itr_1 )= &
                                              parts_obj% pressure_parts( itr )
        ELSEIF( parts_obj% pos( 1, itr ) > 0 )THEN
          itr_2= itr_2 + 1
          parts_obj% pos_x2( itr_2 )= parts_obj% pos( 1, itr )
          parts_obj% pressure_parts_x2( itr_2 )= &
                                              parts_obj% pressure_parts( itr )
        ENDIF

      ENDIF
    ENDDO

    !PRINT *, "3"

    DO itr= 3, cntr1 - 3, 1
      parts_obj% pressure_parts_x_der1( itr - 2 )=&
                     ( + parts_obj% pressure_parts_x1( itr - 2 )/12.0D0 &
                       - 2.0*parts_obj% pressure_parts_x1( itr - 1 )/3.0D0 &
                       + 2.0*parts_obj% pressure_parts_x1( itr + 1 )/3.0D0 &
                       - parts_obj% pressure_parts_x1( itr + 2 )/12.0D0 )&
                       /( Msun_geo*km2m*ABS( parts_obj% pos_x1( itr ) - &
                                             parts_obj% pos_x1( itr - 1 ) ) )

      parts_obj% pressure_length_scale_x1( itr - 2 )= &
                          ABS( parts_obj% pressure_parts_x1( itr - 2 )/ &
                               parts_obj% pressure_parts_x_der1( itr - 2 ) )

      !PRINT *, "p1=", parts_obj% pressure_parts_x1( itr - 2 )
      !PRINT *, "p_r1=", parts_obj% pressure_parts_x_der1( itr - 2 )
      !PRINT *, "p/p_r1=", parts_obj% pressure_length_scale_x1( itr - 2 )
      !PRINT *

    ENDDO
    DO itr= 3, cntr2 - 3, 1
      parts_obj% pressure_parts_x_der2( itr - 2 )=&
                     ( + parts_obj% pressure_parts_x2( itr - 2 )/12.0D0 &
                       - 2.0*parts_obj% pressure_parts_x2( itr - 1 )/3.0D0 &
                       + 2.0*parts_obj% pressure_parts_x2( itr + 1 )/3.0D0 &
                       - parts_obj% pressure_parts_x2( itr + 2 )/12.0D0 )&
                       /( Msun_geo*km2m*ABS( parts_obj% pos_x2( itr ) - &
                                             parts_obj% pos_x2( itr - 1 ) ) )

      parts_obj% pressure_length_scale_x2( itr - 2 )= &
                          ABS( parts_obj% pressure_parts_x2( itr - 2 )/ &
                               parts_obj% pressure_parts_x_der2( itr - 2 ) )

      !PRINT *, "p2=", parts_obj% pressure_parts_x2( itr - 2 )
      !PRINT *, "p_r2=", parts_obj% pressure_parts_x_der2( itr - 2 )
      !PRINT *, "p/p_r2=", parts_obj% pressure_length_scale_x2( itr - 2 )
      !PRINT *

    ENDDO

    PRINT *, " * Maximum typical length scale for change in pressure", &
             " along the x axis for NS 1= ", &
             MAXVAL( parts_obj% pressure_length_scale_x1, DIM= 1 )/km2m, " km"
    PRINT *, " * Minimum typical length scale for change in pressure", &
             " along the x axis for NS 1= ", &
             MINVAL( parts_obj% pressure_length_scale_x1, DIM= 1 )/km2m, " km"
    PRINT *
    PRINT *, " * Maximum typical length scale for change in pressure", &
             " along the x axis for NS 2= ", &
             MAXVAL( parts_obj% pressure_length_scale_x2, DIM= 1 )/km2m, " km"
    PRINT *, " * Minimum typical length scale for change in pressure", &
             " along the x axis for NS 2= ", &
             MINVAL( parts_obj% pressure_length_scale_x2, DIM= 1 )/km2m, " km"
    PRINT *

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

!===============================================================================
!===============================DEBUGGING=======================================
!===============================================================================

    namefile= "dbg-hydro.dat"

    INQUIRE( FILE= TRIM(namefile), EXIST= exist )

    IF( exist )THEN
        OPEN( UNIT= 2, FILE= TRIM(namefile), STATUS= "REPLACE", &
              FORM= "FORMATTED", &
              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
              IOMSG= err_msg )
    ELSE
        OPEN( UNIT= 2, FILE= TRIM(namefile), STATUS= "NEW", &
              FORM= "FORMATTED", &
              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening " // TRIM(namefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when opening " &
    !                  // TRIM(namefile) )

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Values of the fields (including coordinates) exported by LORENE "&
    // "on each grid point"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 1 in " // TRIM(namefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 1 in "&
    !        // TRIM(namefile) )

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# column:      1        2       3       4       5", &
    "       6       7       8", &
    "       9       10      11", &
    "       12      13      14", &
    "       15      16      17      18"

    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 2 in " // TRIM(namefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 2 in "&
    !            // TRIM(namefile) )

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "#      grid point      x [km]       y [km]       z [km]       lapse", &
    "       shift_x [c]    shift_y [c]    shift_z [c]", &
    "       baryon density in the local rest frame [kg m^{-3}$]", &
    "       energy density [c^2]", &
    "       specific energy [c^2]", &
    "       pressure [Pa]", &
    "       fluid 3-velocity wrt the Eulerian observer (3 columns) [c]", &
    "       fluid coordinate 3-velocity vel_u (3 columns) [c]", &
    "       baryon number per particle nu", &
    "       baryon density in the local rest frame nlrf [baryon/Msun_geo^3]", &
    "       electron fraction", &
    "       generalized Lorentz factor Theta"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 3 in " // TRIM(namefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 3 in "&
    !          // TRIM(namefile) )

    DO itr = 1, parts_obj% npart, 1
      abs_pos( 1, itr )= ABS( parts_obj% pos( 1, itr ) )
      abs_pos( 2, itr )= ABS( parts_obj% pos( 2, itr ) )
      abs_pos( 3, itr )= ABS( parts_obj% pos( 3, itr ) )
    ENDDO

    min_y_index= 0
    min_abs_y= 1D+20
    DO itr = 1, parts_obj% npart, 1
      IF( ABS( parts_obj% pos( 2, itr ) ) < min_abs_y )THEN
        min_abs_y= ABS( parts_obj% pos( 2, itr ) )
        min_y_index= itr
      ENDIF
    ENDDO

    min_abs_z= MINVAL( abs_pos( 3, : ) )

    write_data_loop: DO itr = 1, parts_obj% npart, 1

      IF( parts_obj% export_form_xy .AND. &
          parts_obj% pos( 3, itr ) /= min_abs_z )THEN
        CYCLE
      ENDIF
      IF( parts_obj% export_form_x .AND. &
          ( parts_obj% pos( 3, itr ) /= min_abs_z &
          .OR. parts_obj% pos( 2, itr ) /= parts_obj% pos( 2, min_y_index ) ) )THEN
        CYCLE
      ENDIF
      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        itr, &
        parts_obj% pos( 1, itr ), &
        parts_obj% pos( 2, itr ), &
        parts_obj% pos( 3, itr ), &
        parts_obj% lapse_parts( itr ), &
        parts_obj% shift_parts_x( itr ), &
        parts_obj% shift_parts_y( itr ), &
        parts_obj% shift_parts_z( itr ), &
        parts_obj% baryon_density_parts( itr ), &
        parts_obj% energy_density_parts( itr ), &
        parts_obj% specific_energy_parts( itr ), &
        parts_obj% pressure_parts( itr ), &
        parts_obj% v_euler_parts_x( itr ), &
        parts_obj% v_euler_parts_y( itr ), &
        parts_obj% v_euler_parts_z( itr )

    IF( ios > 0 )THEN
      PRINT *, "...error when writing the arrays in " // TRIM(namefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing " &
    !         // "the arrays in " // TRIM(finalnamefile) )
    ENDDO write_data_loop

    CLOSE( UNIT= 2 )

  END PROCEDURE construct_particles


  MODULE PROCEDURE place_particles_3dlattice

    !*****************************************************
    !                                                    *
    ! Compute positions in a 3D lattice around the stars *
    !                                                    *
    ! FT 5.10.2020                                       *
    !                                                    *
    !*****************************************************

    USE constants,    ONLY: Msun_geo

    IMPLICIT NONE

    INTEGER:: ix, iy, iz, sgn, npart_half

    DOUBLE PRECISION:: dx, dy, dz
    DOUBLE PRECISION:: xtemp, ytemp, ztemp, zlim
    DOUBLE PRECISION:: max_baryon_density, thres_baryon_density

    PRINT *, "** Executing the place_particles_3dlattice " &
             // "subroutine..."
    PRINT *

    !
    !-- Set the boundary in z
    !
    IF( ABS(zmax) > ABS(zmin) )THEN
      zlim= zmax
    ELSE
      zlim= zmin
    ENDIF

    !
    !-- Compute lattice steps
    !
    dx= ABS(xmax - xmin)/DBLE( THIS% nx )
    dy= ABS(ymax - ymin)/DBLE( THIS% ny )
    dz= ABS(zlim)/DBLE( THIS% nz/2 )

    PRINT *, " * dx=", dx,  ", dy=", dx,  ", dz=", dz
    PRINT *, " * nx=", THIS% nx,  ", ny=", THIS% ny,  ", nz=", THIS% nz
    PRINT *

    THIS% npart_temp = THIS% nx*THIS% ny*THIS% nz

    PRINT *, " * Number of lattice points= nx*ny*nz=", THIS% npart_temp
    PRINT *

    !
    !-- Compute the mass density at the center of the stars
    !
    max_baryon_density= MAX( bns_obj% get_rho_center1(), &
                             bns_obj% get_rho_center2() )

    !
    !-- Set the threshold above which a lattice point is
    !-- promoted to a particle
    !
    IF( THIS% use_thres )THEN
      thres_baryon_density= max_baryon_density/thres
    ELSE
      thres_baryon_density= 0.0D0
    ENDIF

    ! Allocating the memory for the array pos( 3, npart_temp )
    ! Note that after determining npart, the array pos is reshaped into
    ! pos( 3, npart )
    IF(.NOT.ALLOCATED( THIS% pos ))THEN
      ALLOCATE( THIS% pos( 3, THIS% npart_temp ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array pos in SUBROUTINE" &
                  // "place_particles_3D_lattice. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array pos in SUBROUTINE" &
      !                // "place_particles_3D_lattice." )
    ENDIF
    ! Initializing the array pos to 0
    THIS% pos= 0.0D0

    !---------------------------------------------------------!
    !--  Storing the particle positions into the array pos  --!
    !--  symmetrically w.r.t. the xy plane                  --!
    !---------------------------------------------------------!

    PRINT *, " * Placing particles around NSs..."
    PRINT *

    THIS% npart= 0
    THIS% npart1= 0
    THIS% npart2= 0
    !
    !-- Choose the larger value for the boundary in z
    !
    IF( zlim == zmin )THEN
      sgn= - 1
    ELSE
      sgn= 1
    ENDIF

    !
    !-- Place the first half of the particle (above or below the xy plane)
    !
    particle_pos_z: DO iz= 1, THIS% nz/2, 1

      ztemp= sgn*( dz/2 + ( iz - 1 )*dz )

      particle_pos_y: DO iy= 1, THIS% ny, 1

        ytemp= ymin + ( iy - 1 )*dy

        particle_pos_x: DO ix= 1, THIS% nx, 1

          xtemp= xmin + dx/2 + ( ix - 1 )*dx

          !
          !-- Promote a lattice point to a particle,
          !-- if the mass density is higher than the threshold
          !
          IF( bns_obj% import_mass_density( xtemp, ytemp, ztemp ) &
                                  > thres_baryon_density &
              .AND. &
              bns_obj% is_hydro_negative( xtemp, ytemp, ztemp ) == 0 )THEN

            THIS% npart= THIS% npart + 1
            IF( xtemp < 0 )THEN
              THIS% npart1= THIS% npart1 + 1
            ELSEIF( xtemp > 0 )THEN
              THIS% npart2= THIS% npart2 + 1
            ENDIF
            THIS% pos( 1, THIS% npart )= xtemp
            THIS% pos( 2, THIS% npart )= ytemp
            THIS% pos( 3, THIS% npart )= ztemp

          ENDIF

          ! Print progress on screen, every 10%
          perc= 50*( THIS% nx*THIS% ny*iz + THIS% nx*iy + ix )/ &
                  ( THIS% nx*THIS% ny*THIS% nz/2 )
          IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
            WRITE( *, "(A2,I3,A1)", ADVANCE= "NO" ) &
                   creturn//" ", perc, "%"
          ENDIF

         ENDDO particle_pos_x
      ENDDO particle_pos_y
    ENDDO particle_pos_z
    WRITE( *, "(A1)", ADVANCE= "NO" ) creturn
    npart_half= THIS% npart
    IF( npart_half == 0 )THEN
      PRINT *, "** There are no particles! Execution stopped..."
      PRINT *
      STOP
    ENDIF
    !
    !-- Place the second half of the particles, mirroring the first half
    !-- w.r.t the xy plane
    !
    particle_pos_z_mirror: DO iz= 1, npart_half, 1

      xtemp=   THIS% pos( 1, iz )
      ytemp=   THIS% pos( 2, iz )
      ztemp= - THIS% pos( 3, iz )

      ! TODO: is this check needed?
      !IF( import_mass_density( xtemp, ytemp, ztemp ) &
      !                               > thres_baryon_density )THEN

      THIS% npart= THIS% npart + 1
      IF( xtemp < 0 )THEN
        THIS% npart1= THIS% npart1 + 1
      ELSEIF( xtemp > 0 )THEN
        THIS% npart2= THIS% npart2 + 1
      ENDIF
      THIS% pos( 1, THIS% npart )= xtemp
      THIS% pos( 2, THIS% npart )= ytemp
      THIS% pos( 3, THIS% npart )= ztemp

      !ENDIF

      ! Print progress on screen, every 10%
      perc= 50 + 50*iz/( npart_half )
      IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
         WRITE( *, "(A2,I3,A1)", ADVANCE= "NO" ) &
                 creturn//" ", perc, "%"
      ENDIF
    ENDDO particle_pos_z_mirror
    WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

    !
    !-- Consistency checks
    !
    IF( THIS% npart /= 2*npart_half )THEN
      PRINT *
      PRINT *, "** ERROR: The number of particles ", THIS% npart, &
               " is not the expected value ", 2*npart_half
      PRINT *
      STOP
    ENDIF

    DO iz= 1, npart_half, 1
      IF( THIS% pos( 3, iz ) /= - THIS% pos( 3, npart_half + iz ) )THEN
        PRINT *
        PRINT *, "** ERROR: The lattice is not mirrored " &
                 // "by the xy plane."
        PRINT *
        STOP
      ENDIF
    ENDDO

    IF( THIS% npart1 + THIS% npart2 /= THIS% npart )THEN
      PRINT *, "** ERROR: npart1 + npart2 /= npart"
      PRINT *, " * npart1=", THIS% npart1
      PRINT *, " * npart2=", THIS% npart2
      PRINT *, " * npart1 + npart2=", THIS% npart1 + THIS% npart2
      PRINT *, " * npart=", THIS% npart
      STOP
    ENDIF

    PRINT *, " * Particles placed. Number of particles=", &
             THIS% npart, "=", DBLE(THIS% npart)/DBLE(THIS% npart_temp), &
             " of the points in lattice."
    PRINT *
    PRINT *, " * Number of particles on NS 1=", THIS% npart1
    PRINT *, " * Number of particles on NS 2=", THIS% npart2
    PRINT *

    !
    !-- Computing total volume and volume per particle
    !
    IF(.NOT.ALLOCATED( THIS% pvol ))THEN
      ALLOCATE( THIS% pvol( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array pvol ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for array v_euler_parts_z" )
    ENDIF

    THIS% vol  = (xmax - xmin)*(ymax - ymin)*2*ABS(zlim)
    THIS% vol_a= THIS% vol/THIS% npart_temp

    THIS% pvol= THIS% vol_a

    ! Consistency check for the particle volume
    IF( ABS( THIS% vol_a - dx*dy*dz ) > 1D-9 )THEN
      PRINT *, " * The particle volume vol_a=", THIS% vol_a, "Msun_geo^3"
      PRINT *, " is not equal to dx*dy*dz=", dx*dy*dz, "Msun_geo^3."
      PRINT *
      STOP
    ENDIF

    PRINT *, " * Total volume of the lattices=", THIS% vol, "Msun_geo^3"
    PRINT *, " * Particle volume=", THIS% vol_a, "Msun_geo^3"
    PRINT *

    PRINT *, "** Subroutine place_particles_3D_lattice executed."
    PRINT *

  END PROCEDURE place_particles_3dlattice


  MODULE PROCEDURE place_particles_3dlattices

    !*****************************************************
    !                                                    *
    ! Compute positions in two 3D lattices each one      *
    ! around each star                                   *
    !                                                    *
    ! FT 19.10.2020                                      *
    !                                                    *
    !*****************************************************

    USE constants, ONLY: Msun_geo

    IMPLICIT NONE

    INTEGER:: ix, iy, iz, nx2, ny2, nz2, sgn, npart_half, npart_half2, itr2

    DOUBLE PRECISION:: dx1, dy1, dz1, dx2, dy2, dz2
    DOUBLE PRECISION:: xtemp, ytemp, ztemp, zlim, zlim2
    DOUBLE PRECISION:: max_baryon_density1, thres_baryon_density1
    DOUBLE PRECISION:: max_baryon_density2, thres_baryon_density2
    ! Variable used to compute the volume of a particle in an alternative way
    ! to perform a consistency check
    DOUBLE PRECISION:: vol_a_alt1, vol_a_alt2

    CHARACTER( LEN= 50 ):: lorene_bns_id_parfile

    PRINT *, "** Executing the place_particles_3dlattices " &
             // "subroutine..."
    PRINT *

    !
    !-- Set the boundaries in z, for lattice 1 and lattice 2
    !
    IF( ABS(zmax1) > ABS(zmin1) )THEN
      zlim= zmax1
    ELSE
      zlim= zmin1
    ENDIF
    IF( ABS(zmax2) > ABS(zmin2) )THEN
      zlim2= zmax2
    ELSE
      zlim2= zmin2
    ENDIF

    THIS% mass1= bns_obj% get_mass1()
    THIS% mass2= bns_obj% get_mass2()

    IF( THIS% mass1 > THIS% mass2 )THEN

      ! mass_ratio < 1
      THIS% mass_ratio= THIS% mass2/THIS% mass1
      !
      !-- Compute lattices' steps
      !
      THIS% nx2= THIS% nx
      THIS% ny2= THIS% ny
      THIS% nz2= THIS% nz
      dx2= ABS(xmax2 - xmin2)/DBLE( THIS% nx2 )
      dy2= ABS(ymax2 - ymin2)/DBLE( THIS% ny2 )
      dz2= ABS(zlim2)/DBLE( THIS% nz2/2 )

      dx1= dx2*(THIS% mass_ratio**(1.0D0/3.0D0))
      dy1= dy2*(THIS% mass_ratio**(1.0D0/3.0D0))
      dz1= dz2*(THIS% mass_ratio**(1.0D0/3.0D0))
      THIS% nx1= NINT( ABS(xmax1 - xmin1)/dx1 ) + 1
      THIS% ny1= NINT( ABS(ymax1 - ymin1)/dy1 ) + 1
      THIS% nz1= NINT( 2*ABS(zlim)/dz1 ) + 1

    ELSE

      ! mass_ratio < 1
      THIS% mass_ratio= THIS% mass1/THIS% mass2
      !
      !-- Compute lattices' steps
      !
      THIS% nx1= THIS% nx
      THIS% ny1= THIS% ny
      THIS% nz1= THIS% nz
      dx1= ABS(xmax1 - xmin1)/DBLE( THIS% nx1 )
      dy1= ABS(ymax1 - ymin1)/DBLE( THIS% ny1 )
      dz1= ABS(zlim)/DBLE( THIS% nz1/2 )

      dx2= dx1*(THIS% mass_ratio**(1.0D0/3.0D0))
      dy2= dy1*(THIS% mass_ratio**(1.0D0/3.0D0))
      dz2= dz1*(THIS% mass_ratio**(1.0D0/3.0D0))
      THIS% nx2= NINT( ABS(xmax2 - xmin2)/dx2 ) + 1
      THIS% ny2= NINT( ABS(ymax2 - ymin2)/dy2 ) + 1
      THIS% nz2= NINT( 2*ABS(zlim2)/dz2 ) + 1

    ENDIF

    ! Set the number of particles in the z direction to an even number
    ! since half of the particles are above the xy plane, and half below it
    IF( MOD( THIS% nz2, 2 ) /= 0 )THEN
      THIS% nz2= THIS% nz2 - 1
    ENDIF

    PRINT *, " * dx1=", dx1,  ", dy1=", dx1,  ", dz1=", dz1
    PRINT *, " * dx2=", dx2,  ", dy2=", dx2,  ", dz2=", dz2
    PRINT *, " * nx1=", THIS% nx1, ", ny1=", THIS% ny1, ", nz1=", THIS% nz1
    PRINT *, " * nx2=", THIS% nx2, ", ny2=", THIS% ny2, ", nz2=", THIS% nz2
    PRINT *

    !PRINT *, " * xmin1=", xmin1, ", xmax1=", xmax1
    !PRINT *, " * xmin2=", xmin2, ", xmax2=", xmax2
    !PRINT *
    !STOP

    ! Compute number of lattice points (temporary particle number)
    THIS% npart1_temp = THIS% nx*THIS% ny*THIS% nz !+ THIS% nx*THIS% ny
    THIS% npart2_temp = THIS% nx2*THIS% ny2*THIS% nz2 !+ nx2*ny2
    THIS% npart_temp  = THIS% npart1_temp + THIS% npart2_temp

    PRINT *, " * Number of points for lattice 1= nx1*ny1*nz1=", &
             THIS% npart1_temp
    PRINT *, " * Number of points for lattice 2= nx2*ny2*nz2=", &
             THIS% npart2_temp
    PRINT *

    !
    !-- Compute the mass density at the center of the stars
    !

    ! N.B. The following two densities are in LORENE units [kg m^{-3}]
    !max_baryon_density1= bns_obj% import_mass_density( &
    !                                 bns_obj% get_center1_x(), 0.0D0, 0.0D0 )
    !max_baryon_density2= bns_obj% import_mass_density( &
    !                                 bns_obj% get_center2_x(), 0.0D0, 0.0D0 )

    ! The following two density ar in SPHINCS units [Msun Msun_geo^{-3}]
    max_baryon_density1= bns_obj% get_rho_center1()
    max_baryon_density2= bns_obj% get_rho_center2()

    !
    !-- Set the thresholds above which a lattice point is
    !-- promoted to a particle
    !
    IF( THIS% use_thres )THEN
      thres_baryon_density1= max_baryon_density1/thres
      thres_baryon_density2= max_baryon_density2/thres
    ELSE
      thres_baryon_density1= 0.0D0
      thres_baryon_density2= 0.0D0
    ENDIF

    ! Allocating the memory for the array pos( 3, npart_temp )
    ! Note that after determining npart, the array pos is reshaped into
    ! pos( 3, npart )
    IF(.NOT.ALLOCATED( THIS% pos ))THEN
      ALLOCATE( THIS% pos( 3, THIS% npart_temp ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array pos in SUBROUTINE" &
                  // "place_particles_3D_lattices. ", &
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

    !
    !-- Placing particles on NS 1
    !
    PRINT *, "Placing particles on NS 1..."
    PRINT *
    THIS% npart= 0
    THIS% npart1= 0
    !
    !-- Choose the larger value for the boundary in z
    !
    IF( zlim == zmin1 )THEN
      sgn= - 1
    ELSE
      sgn= 1
    ENDIF
    !
    !-- Place the first half of the particle (above or below the xy plane)
    !
    particle_pos_z1: DO iz= 1, THIS% nz/2, 1

      ztemp= sgn*( dz1/2 + ( iz - 1 )*dz1 )

      particle_pos_y1: DO iy= 1, THIS% ny, 1

        ytemp= ymin1 + dy1/2 + ( iy - 1 )*dy1

        particle_pos_x1: DO ix= 1, THIS% nx, 1

          xtemp= xmin1 + dx1/2 + ( ix - 1 )*dx1

          !
          !-- Promote a lattice point to a particle,
          !-- if the mass density is higher than the threshold
          !
          IF( bns_obj% import_mass_density( xtemp, ytemp, ztemp ) &
                                > thres_baryon_density1 &
              .AND. &
              bns_obj% is_hydro_negative( xtemp, ytemp, ztemp ) == 0 )THEN

            THIS% npart = THIS% npart + 1
            THIS% npart1= THIS% npart1 + 1
            THIS% pos( 1, THIS% npart )= xtemp
            THIS% pos( 2, THIS% npart )= ytemp
            THIS% pos( 3, THIS% npart )= ztemp

          ENDIF

          ! Print progress on screen, every 10%
          perc= 50*( THIS% nx*THIS% ny*iz + THIS% nx*iy + ix )/ &
                  ( THIS% nx*THIS% ny*THIS% nz/2 )
          IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
            WRITE( *, "(A2,I3,A1)", ADVANCE= "NO" ) &
                   creturn//" ", perc, "%"
           ENDIF
         ENDDO particle_pos_x1
      ENDDO particle_pos_y1
    ENDDO particle_pos_z1
    WRITE( *, "(A1)", ADVANCE= "NO" ) creturn
    npart_half= THIS% npart
    IF( npart_half == 0 )THEN
      PRINT *, "** There are no particles on star 1! Execution stopped..."
      PRINT *
      STOP
    ENDIF
    !
    !-- Place the second half of the particles, mirroring the first half
    !-- w.r.t the xy plane
    !
    particle_pos_z1_mirror: DO iz= 1, npart_half, 1

      xtemp=   THIS% pos( 1, iz )
      ytemp=   THIS% pos( 2, iz )
      ztemp= - THIS% pos( 3, iz )

      ! TODO: is this check needed?
      !IF( import_mass_density( xtemp, ytemp, ztemp ) &
      !                               > thres_baryon_density1 )THEN

      THIS% npart = THIS% npart + 1
      THIS% npart1= THIS% npart1 + 1
      THIS% pos( 1, THIS% npart )= xtemp
      THIS% pos( 2, THIS% npart )= ytemp
      THIS% pos( 3, THIS% npart )= ztemp

      !ENDIF

      perc= 50 + 50*iz/( npart_half )
      IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
        WRITE( *, "(A2,I3,A1)", ADVANCE= "NO" ) &
               creturn//" ", perc, "%"
      ENDIF

    ENDDO particle_pos_z1_mirror
    WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

    !
    !-- Place the particles on the xy plane
    !
    !ztemp= 0.0D0
    !
    !particle_pos_y1_xy: DO iy= 1, THIS% ny, 1
    !
    !  ytemp= ymin1 + dy/2 + ( iy - 1 )*dy
    !
    !  particle_pos_x1_xy: DO ix= 1, THIS% nx, 1
    !
    !    xtemp= xmin1 + dx/2 + ( ix - 1 )*dx
    !
    !    !
    !    !-- Promote a lattice point to a particle,
    !    !-- if the mass density is higher than the threshold
    !    !
    !    IF( bns_obj% import_mass_density( xtemp, ytemp, ztemp ) &
    !                          > thres_baryon_density1 )THEN
    !
    !      THIS% npart = THIS% npart + 1
    !      THIS% npart1= THIS% npart1 + 1
    !      THIS% pos( 1, THIS% npart )= xtemp
    !      THIS% pos( 2, THIS% npart )= ytemp
    !      THIS% pos( 3, THIS% npart )= ztemp
    !
    !    ENDIF
    !
    !    ! Print progress on screen, every 10%
    !    perc= 50*( THIS% nx*THIS% ny*iz + THIS% nx*iy + ix )/ &
    !            ( THIS% nx*THIS% ny*THIS% nz/2 )
    !    IF( MOD( perc, 10 ) == 0 )THEN
    !      WRITE( *, "(A2,I3,A1)", ADVANCE= "NO" ) &
    !             creturn//" ", perc, "%"
    !     ENDIF
    !   ENDDO particle_pos_x1_xy
    !ENDDO particle_pos_y1_xy
    !WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

    !
    !-- Placing particles on NS 2 with the same algorithm as for NS 1
    !
    PRINT *, "Placing particles on NS 2..."
    PRINT *
    THIS% npart2= 0
    IF( zlim2 == zmin2 )THEN
      sgn= - 1
    ELSE
      sgn= 1
    ENDIF
    particle_pos_z2: DO iz= 1, THIS% nz2/2, 1

      ztemp= sgn*( dz2/2 + ( iz - 1 )*dz2 )

      particle_pos_y2: DO iy= 1, THIS% ny2, 1

        ytemp= ymin2 + dy2/2 + ( iy - 1 )*dy2

        particle_pos_x2: DO ix= 1, THIS% nx2, 1

          xtemp= xmin2 + dx2/2 + ( ix - 1 )*dx2

          IF( bns_obj% import_mass_density( xtemp, ytemp, ztemp ) &
                                  > thres_baryon_density2 &
              .AND. &
              bns_obj% is_hydro_negative( xtemp, ytemp, ztemp ) == 0 )THEN

            THIS% npart = THIS% npart + 1
            THIS% npart2= THIS% npart2 + 1
            THIS% pos( 1, THIS% npart )= xtemp
            THIS% pos( 2, THIS% npart )= ytemp
            THIS% pos( 3, THIS% npart )= ztemp

          ENDIF

          ! Print progress on screen, every 10%
          perc= 50*( THIS% nx2*THIS% ny2*( iz - 1 ) + THIS% nx2*( iy - 1 ) &
                + ix )/( THIS% nx2*THIS% ny2*THIS% nz2/2 )
          IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
            WRITE( *, "(A2,I3,A1)", ADVANCE= "NO" ) &
                   creturn//" ", perc, "%"
          ENDIF
        ENDDO particle_pos_x2
      ENDDO particle_pos_y2
    ENDDO particle_pos_z2
    WRITE( *, "(A1)", ADVANCE= "NO" ) creturn
    npart_half2= THIS% npart
    IF( npart_half2 == 2*npart_half )THEN
      PRINT *, "** There are no particles on star 2! Execution stopped..."
      PRINT *
      STOP
    ENDIF
    particle_pos_z2_mirror: DO iz= 2*npart_half + 1, npart_half2, 1

      xtemp=   THIS% pos( 1, iz )
      ytemp=   THIS% pos( 2, iz )
      ztemp= - THIS% pos( 3, iz )

      !IF( import_mass_density( xtemp, ytemp, ztemp ) &
      !                               > thres_baryon_density2 )THEN

      THIS% npart = THIS% npart + 1
      THIS% npart2= THIS% npart2 + 1
      THIS% pos( 1, THIS% npart )= xtemp
      THIS% pos( 2, THIS% npart )= ytemp
      THIS% pos( 3, THIS% npart )= ztemp

      !ENDIF

      ! Print progress on screen, every 10%
      perc= 50 + 50*( iz - 2*npart_half + 1 ) &
                    /( npart_half2 - 2*npart_half )
      IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
        WRITE( *, "(A2,I3,A1)", ADVANCE= "NO" ) &
                creturn//" ", perc, "%"
      ENDIF

    ENDDO particle_pos_z2_mirror
    WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

    !
    !-- Place the particles on the xy plane
    !
    !ztemp= 0.0D0
    !
    !particle_pos_y2_xy: DO iy= 1, ny2, 1
    !
    !  ytemp= ymin2 + dy/2 + ( iy - 1 )*dy
    !
    !  particle_pos_x2_xy: DO ix= 1, nx2, 1
    !
    !    xtemp= xmin2 + dx/2 + ( ix - 1 )*dx
    !
    !    IF( bns_obj% import_mass_density( xtemp, ytemp, ztemp ) &
    !                            > thres_baryon_density2 )THEN
    !
    !      THIS% npart = THIS% npart + 1
    !      THIS% npart2= THIS% npart2 + 1
    !      THIS% pos( 1, THIS% npart )= xtemp
    !      THIS% pos( 2, THIS% npart )= ytemp
    !      THIS% pos( 3, THIS% npart )= ztemp
    !
    !    ENDIF
    !
    !    ! Print progress on screen, every 10%
    !    perc= 50*( nx2*ny2*( iz - 1 ) + nx2*( iy - 1 ) + ix )&
    !          /( nx2*ny2*nz2/2 )
    !    IF( MOD( perc, 10 ) == 0 )THEN
    !      WRITE( *, "(A2,I3,A1)", ADVANCE= "NO" ) &
    !             creturn//" ", perc, "%"
    !    ENDIF
    !  ENDDO particle_pos_x2_xy
    !ENDDO particle_pos_y2_xy
    !WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

    !
    !-- Consistency checks
    !
    IF( THIS% npart /= ( 2*( npart_half2 - npart_half ) ) )THEN
      PRINT *
      PRINT *, "** ERROR: The number of particles ", THIS% npart, &
               " is not the expected value ", &
               2*( npart_half2 - npart_half )
      PRINT *
      STOP
    ENDIF

    DO iz= 1, npart_half, 1
      IF( THIS% pos( 3, iz ) /= - THIS% pos( 3, npart_half + iz ) )THEN
        PRINT *
        PRINT *, "** ERROR: The lattice around NS 1 are not mirrored " &
                 // "by the xy plane."
        PRINT *
        STOP
      ENDIF
    ENDDO
    DO iz= 2*npart_half + 1, npart_half2, 1
      IF( THIS% pos( 3, iz ) /= &
          - THIS% pos( 3, ( npart_half2 - 2*npart_half ) + iz ) )THEN
        PRINT *
        PRINT *, "** ERROR: The lattice around NS 2 are not mirrored " &
                 // "by the xy plane."
        PRINT *
        STOP
      ENDIF
    ENDDO

    IF( THIS% npart1 + THIS% npart2 /= THIS% npart )THEN
      PRINT *, "** ERROR: npart1 + npart2 /= npart"
      PRINT *, " * npart1=", THIS% npart1
      PRINT *, " * npart2=", THIS% npart2
      PRINT *, " * npart1 + npart2=", THIS% npart1 + THIS% npart2
      PRINT *, " * npart=", THIS% npart
      STOP
    ENDIF

    !
    !-- Printouts
    !
    PRINT *, " * Particles placed. Number of particles=", &
             THIS% npart, "=", DBLE(THIS% npart)/DBLE(THIS% npart_temp), &
             " of the points in lattices."
    PRINT *
    PRINT *, " * Number of particles on NS 1=", THIS% npart1, "=", &
             DBLE(THIS% npart1)/DBLE(THIS% npart1_temp), &
             " of the points in the first lattice."
    PRINT *, " * Number of particles on NS 2=", THIS% npart2, "=", &
             DBLE(THIS% npart2)/DBLE(THIS% npart2_temp), &
             " of the points in the second lattice."
    PRINT *

    !
    !-- Computing total volume and volume per particle
    !
    IF(.NOT.ALLOCATED( THIS% pvol ))THEN
      ALLOCATE( THIS% pvol( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array pvol ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for array v_euler_parts_z" )
    ENDIF

    THIS% vol1_a= dx1*dy1*dz1
    THIS% vol1 = (xmax1 - xmin1)*(ymax1 - ymin1)*2*ABS(zlim)
    !THIS% vol2 = THIS% npart2_temp * THIS% vol_a
    !THIS% vol  = THIS% vol1 + THIS% vol2
    vol_a_alt1  = THIS% vol1/THIS% npart1_temp

    THIS% vol2_a= dx2*dy2*dz2
    THIS% vol2 =  dx2*THIS% nx2*dy2*THIS% ny2*dz2*THIS% nz2
    !THIS% vol2 = (xmax2 - xmin2)*(ymax2 - ymin2)*2*ABS(zlim2)
    !THIS% vol2 = THIS% npart2_temp * THIS% vol_a2
    !THIS% vol  = THIS% vol1 + THIS% vol2
    vol_a_alt2  = THIS% vol2/THIS% npart2_temp

    THIS% pvol( 1:THIS% npart1 )              = THIS% vol1_a
    THIS% pvol( THIS% npart1 + 1:THIS% npart )= THIS% vol2_a

    THIS% vol= THIS% vol1 + THIS% vol2

    ! Consistency check for the particle volume
    IF( ABS( THIS% vol1_a - vol_a_alt1 ) > 1D-7 )THEN
      PRINT *, " * The particle volume vol_a_alt1=", vol_a_alt1, "Msun_geo^3"
      PRINT *, " is not equal to dx1*dy1*dz1=", THIS% vol1_a, "Msun_geo^3."
      PRINT *
      STOP
    ENDIF
    ! Consistency check for the particle volume
    IF( ABS( THIS% vol2_a - vol_a_alt2 ) > 1D-7 )THEN
      PRINT *, " * The particle volume vol_a_alt2=", vol_a_alt2, "Msun_geo^3"
      PRINT *, " is not equal to dx2*dy2*dz2=", THIS% vol2_a, "Msun_geo^3."
      PRINT *
      STOP
    ENDIF

    PRINT *, " * Total volume of the lattices=", THIS% vol, "Msun_geo^3"
    PRINT *, " * Particle volume on NS 1=", THIS% vol1_a, "Msun_geo^3"
    PRINT *, " * Particle volume on NS 2=", THIS% vol2_a, "Msun_geo^3"
    PRINT *

    PRINT *, "** Subroutine place_particles_3dlattices " &
             // "executed."
    PRINT *

  END PROCEDURE place_particles_3dlattices

  MODULE PROCEDURE analyze_hydro

    !************************************************
    !                                               *
    ! Export the points where some of the hydro     *
    ! fields are negative to a formatted file       *
    !                                               *
    ! FT 5.12.2020                                  *
    !                                               *
    !************************************************

    IMPLICIT NONE

    LOGICAL:: exist, negative_hydro

    INQUIRE( FILE= TRIM(namefile), EXIST= exist )

    IF( exist )THEN
      OPEN( UNIT= 20, FILE= TRIM(namefile), STATUS= "REPLACE", &
            FORM= "FORMATTED", &
            POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
            IOMSG= err_msg )
    ELSE
      OPEN( UNIT= 20, FILE= TRIM(namefile), STATUS= "NEW", &
      FORM= "FORMATTED", &
            ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
       PRINT *, "...error when opening ",  TRIM(namefile), &
                " The error message is", err_msg
       STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when opening " &
    !         // TRIM(namefile) )

    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Points where some of the hydro fields are negative. "
    IF( ios > 0 )THEN
       PRINT *, "...error when writing line 1 in ",  TRIM(namefile), &
                " The error message is", err_msg
       STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 1 in "&
    !         // TRIM(namefile) )
    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# column:      1        2       3"
    IF( ios > 0 )THEN
       PRINT *, "...error when writing line 2 in ",  TRIM(namefile), &
                " The error message is", err_msg
       STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 2 in "&
    !        // TRIM(namefile) )
    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "#      x   y   z"
    IF( ios > 0 )THEN
       PRINT *, "...error when writing line 3 in ",  TRIM(namefile), &
                " The error message is", err_msg
       STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 3 in "&
    !        // TRIM(namefile) )

    DO itr= 1, THIS% npart, 1
      IF( THIS% baryon_density_parts ( itr ) < 0 .OR. &
          THIS% energy_density_parts ( itr ) < 0 .OR. &
          THIS% specific_energy_parts( itr ) < 0 .OR. &
          THIS% pressure_parts       ( itr ) < 0 )THEN

        negative_hydro= .TRUE.

        WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, &
               FMT = * )&
            THIS% pos( 1, itr ), &
            THIS% pos( 2, itr ), &
            THIS% pos( 3, itr )

        IF( ios > 0 )THEN
          PRINT *, "...error when writing the arrays in ", TRIM(namefile), &
                   " The error message is", err_msg
          STOP
        ENDIF
        !CALL test_status( ios, err_msg, "...error in writing "&
        !                // "the arrays in " // TRIM(namefile) )
      ENDIF
    ENDDO

    CLOSE( UNIT= 20 )

    IF( negative_hydro )THEN
      PRINT *, "** WARNING! Some of the hydro fields are negative on", &
               " some of the particles! See the file ", namefile, &
               " for the positions of such particles."
      PRINT *
    ELSE
      PRINT *, " * The hydro fields are positive on the particles."
      PRINT *
    ENDIF

  END PROCEDURE analyze_hydro


  MODULE PROCEDURE allocate_lorene_id_parts_memory

    !************************************************
    !                                               *
    ! Allocate memory for the LORENE ID on the      *
    ! particles                                     *
    !                                               *
    ! FT 10.11.2020                                 *
    !                                               *
    !************************************************

    IMPLICIT NONE

    PRINT *, "** Executing allocate_lorene_id_parts_memory."

    IF(.NOT.ALLOCATED( THIS% pos ))THEN
      ALLOCATE( THIS% pos( 3, THIS% npart ), STAT= ios, &
            ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array pos ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array pos" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% lapse_parts ))THEN
      ALLOCATE( THIS% lapse_parts( THIS% npart ), STAT= ios, &
            ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array lapse_parts ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array lapse_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% shift_parts_x ))THEN
      ALLOCATE( THIS% shift_parts_x( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array shift_parts_x ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for shift_parts_x" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% shift_parts_y ))THEN
      ALLOCATE( THIS% shift_parts_y( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array shift_parts_y ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for shift_parts_y" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% shift_parts_z ))THEN
      ALLOCATE( THIS% shift_parts_z( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array shift_parts_z ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for shift_parts_z" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_xx_parts ))THEN
      ALLOCATE( THIS% g_xx_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_xx_parts ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array g_xx_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_xy_parts ))THEN
      ALLOCATE( THIS% g_xy_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_xy_parts ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for array g_xy_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_xz_parts ))THEN
      ALLOCATE( THIS% g_xz_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_xz_parts ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array g_xz_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_yy_parts ))THEN
      ALLOCATE( THIS% g_yy_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_yy_parts ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for array g_yy_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_yz_parts ))THEN
      ALLOCATE( THIS% g_yz_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_yz_parts ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array g_yz_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_zz_parts ))THEN
      ALLOCATE( THIS% g_zz_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_zz_parts ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array g_zz_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% baryon_density_parts ))THEN
      ALLOCATE( THIS% baryon_density_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array baryon_density_parts ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !     "...allocation error for array baryon_density_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% energy_density_parts ))THEN
      ALLOCATE( THIS% energy_density_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array energy_density_parts ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !             "...allocation error for array energy_density_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% specific_energy_parts ))THEN
      ALLOCATE( THIS% specific_energy_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array specific_energy_parts ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !    "...allocation error for array specific_energy_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% pressure_parts ))THEN
      ALLOCATE( THIS% pressure_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array pressure_parts ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array pressure_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% pressure_parts_cu ))THEN
      ALLOCATE( THIS% pressure_parts_cu( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array pressure_parts_cu ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array pressure_parts_cu" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% v_euler_parts_x ))THEN
      ALLOCATE( THIS% v_euler_parts_x( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array v_euler_parts_x ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array v_euler_parts_x" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% v_euler_parts_y ))THEN
      ALLOCATE( THIS% v_euler_parts_y( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array v_euler_parts_y ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for array v_euler_parts_y" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% v_euler_parts_z ))THEN
      ALLOCATE( THIS% v_euler_parts_z( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array v_euler_parts_z ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for array v_euler_parts_z" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% nstar ))THEN
        ALLOCATE( THIS% nstar( THIS% npart ), STAT= ios )
    ENDIF
    IF( ios > 0 )THEN
       PRINT *, '...allocation error for nstar'
       STOP
    ENDIF
    IF(.NOT.ALLOCATED( THIS% nstar_int ))THEN
        ALLOCATE( THIS% nstar_int( THIS% npart ), STAT= ios )
    ENDIF
    IF( ios > 0 )THEN
       PRINT *, '...allocation error for nstar_int'
       STOP
    ENDIF
    IF(.NOT.ALLOCATED( THIS% particle_density ))THEN
        ALLOCATE( THIS% particle_density( THIS% npart ), STAT= ios )
    ENDIF
    IF( ios > 0 )THEN
       PRINT *, '...allocation error for particle_density'
       STOP
    ENDIF
    IF(.NOT.ALLOCATED( THIS% particle_density_int ))THEN
        ALLOCATE( THIS% particle_density_int( THIS% npart ), STAT= ios )
    ENDIF
    IF( ios > 0 )THEN
       PRINT *, '...allocation error for particle_density_int'
       STOP
    ENDIF
    IF(.NOT.ALLOCATED( THIS% pmass ))THEN
      ALLOCATE( THIS% pmass( THIS% npart ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array pmass in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF

    PRINT *, "** Subroutine allocate_lorene_id_parts_memory executed."
    PRINT *

  END PROCEDURE allocate_lorene_id_parts_memory


END SUBMODULE particles_constructor
