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

    USE NaNChecker, ONLY: Check_Array_for_NAN
    USE constants,  ONLY: Msun_geo

    IMPLICIT NONE

    ! The variable counter counts how many times the PROCEDURE
    ! construct_particles is called
    INTEGER, SAVE:: counter= 1
    INTEGER:: nx, ny, nz

    DOUBLE PRECISION:: thres
    DOUBLE PRECISION:: xmin, xmax, ymin, ymax, zmin, zmax, stretch
    DOUBLE PRECISION:: xmin1, xmax1, ymin1, ymax1, zmin1, zmax1
    DOUBLE PRECISION:: xmin2, xmax2, ymin2, ymax2, zmin2, zmax2

    CHARACTER( LEN= : ), ALLOCATABLE:: namefile

    LOGICAL:: file_exists

    NAMELIST /bns_particles/ &
              stretch, &
              nx, ny, nz, &
              thres

    !
    !-- Initialize the timers
    !
    parts_obj% placer_timer      = timer( "placer_timer" )
    parts_obj% importer_timer    = timer( "importer_timer" )
    parts_obj% sph_computer_timer= timer( "sph_computer_timer" )

    ! Declare this object as non-empty (experimental)
    parts_obj% empty_object= .FALSE.

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

    IF( MOD( nz, 2 ) /= 0 )THEN
     PRINT *
     PRINT *, "** ERROR: nz should be even!"
     PRINT *
     STOP
    ENDIF

    IF( nx == 0 .OR. ny == 0 .OR. nz == 0 )THEN
     PRINT *
     PRINT *, "** ERROR: nx, ny, nz cannot be 0!"
     PRINT *
     STOP
    ENDIF

    ! TODO: Add check that the number of rows in placer is the same as the
    !       number of bns objects, and that all bns have a value for placer

    !
    !-- Choose particle placer
    !
    choose_particle_placer: SELECT CASE( dist )

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

    CASE DEFAULT

      PRINT *, "** There is no implemented particle placer " &
               // "corresponding to the number", dist
      STOP

    END SELECT choose_particle_placer

    ! Reshape the array pos by deleting the unnecessary elements
    parts_obj% pos= parts_obj% pos( :, 1:parts_obj% npart )

    ! Allocate needed memory
    CALL allocate_lorene_id_parts_memory( parts_obj )

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
    CALL Check_Array_for_NAN( parts_obj% npart, parts_obj% lapse_parts, &
                                             "lapse_parts" )
    CALL Check_Array_for_NAN( parts_obj% npart, parts_obj% shift_parts_x, &
                                             "shift_parts_x" )
    CALL Check_Array_for_NAN( parts_obj% npart, parts_obj% shift_parts_y, &
                                             "shift_parts_y" )
    CALL Check_Array_for_NAN( parts_obj% npart, parts_obj% shift_parts_z, &
                                             "shift_parts_z" )
    CALL Check_Array_for_NAN( parts_obj% npart, parts_obj% g_xx_parts, &
                                             "g_xx_parts" )
    CALL Check_Array_for_NAN( parts_obj% npart, parts_obj% g_xy_parts, &
                                             "g_xy_parts" )
    CALL Check_Array_for_NAN( parts_obj% npart, parts_obj% g_xz_parts, &
                                             "g_xz_parts" )
    CALL Check_Array_for_NAN( parts_obj% npart, parts_obj% g_yy_parts, &
                                             "g_yy_parts" )
    CALL Check_Array_for_NAN( parts_obj% npart, parts_obj% g_yz_parts, &
                                             "g_yz_parts" )
    CALL Check_Array_for_NAN( parts_obj% npart, parts_obj% g_zz_parts, &
                                             "g_zz_parts" )
    CALL Check_Array_for_NAN( parts_obj% npart, &
            parts_obj% baryon_density_parts, "baryon_density_parts" )
    CALL Check_Array_for_NAN( parts_obj% npart, &
            parts_obj% energy_density_parts, "energy_density_parts" )
    CALL Check_Array_for_NAN( parts_obj% npart, &
            parts_obj% specific_energy_parts, "specific_energy_parts" )
    CALL Check_Array_for_NAN( parts_obj% npart, &
                   parts_obj% pressure_parts, "pressure_parts" )
    CALL Check_Array_for_NAN( parts_obj% npart, &
                  parts_obj% v_euler_parts_x, "v_euler_parts_x" )
    CALL Check_Array_for_NAN( parts_obj% npart, &
                  parts_obj% v_euler_parts_y, "v_euler_parts_y" )
    CALL Check_Array_for_NAN( parts_obj% npart, &
                  parts_obj% v_euler_parts_z, "v_euler_parts_z" )

    !
    !-- Replace the points with negative hydro fields near the surface
    !-- with vacuum
    !
    PRINT *, "** Cleaning LORENE hydro ID around the surfaces of the stars..."
    DO itr= 1, parts_obj% npart, 1

          IF(      parts_obj% baryon_density_parts ( itr ) < 0.0D0 &
              .OR. parts_obj% energy_density_parts ( itr ) < 0.0D0 &
              .OR. parts_obj% specific_energy_parts( itr ) < 0.0D0 &
              .OR. parts_obj% pressure_parts       ( itr ) < 0.0D0 )THEN
              parts_obj% baryon_density_parts ( itr )= 0.0D0
              parts_obj% energy_density_parts ( itr )= 0.0D0
              parts_obj% specific_energy_parts( itr )= 0.0D0
              parts_obj% pressure_parts       ( itr )= 0.0D0
              parts_obj% v_euler_parts_x      ( itr )= 0.0D0
              parts_obj% v_euler_parts_y      ( itr )= 0.0D0
              parts_obj% v_euler_parts_z      ( itr )= 0.0D0
          ENDIF

        ! Print progress on screen
        perc= 100*( itr/parts_obj% npart )
        IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
          WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
                  creturn//" ", perc, "%"
        ENDIF

    ENDDO
    WRITE( *, "(A1)", ADVANCE= "NO" ) creturn
    PRINT *, " * LORENE hydro ID cleaned."
    PRINT *

    !
    !-- Check that the baryon density, the energy density, the specific
    !-- energy, and the pressure are positive. If not, print the coordinates
    !-- of the points to a formatted file
    !
    !namefile= "negative-hydro-particles.dat"
    !CALL parts_obj% analyze_hydro( namefile )

    ! Increase the counter that identifies the particle distribution
    counter= counter + 1

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

    !
    !-- Compute the mass density almost at the center of the stars
    !
    max_baryon_density= MAX( &
                bns_obj%import_mass_density( - bns_obj% get_distance()/2.0D0, &
                                               DBLE(0), DBLE(0) ), &
                bns_obj%import_mass_density(   bns_obj% get_distance()/2.0D0, &
                                               DBLE(0), DBLE(0) ) )

    !
    !-- Set the threshold above which a lattice point is
    !-- promoted to a particle
    !
    thres_baryon_density= max_baryon_density/thres

    ! Allocating the memory for the array pos( 3, npart_temp )
    ! Note that after determining npart, the array pos is reshaped into
    ! pos( 3, npart )
    IF(.NOT.ALLOCATED( THIS% pos ))THEN
      ALLOCATE( THIS% pos( 3, THIS% npart_temp ), STAT= ios, &
                ERRMSG= err_msg )
      CALL test_status( ios, err_msg, &
                      "...allocation error for array pos in SUBROUTINE" &
                      // "place_particles_3D_lattice." )
    ENDIF
    ! Initializing the array pos to 0
    THIS% pos= 0.0D0

    !---------------------------------------------------------!
    !--  Storing the particle positions into the array pos  --!
    !--  symmetrically w.r.t. the xy plane                  --!
    !---------------------------------------------------------!

    PRINT *, "Placing particles around NSs..."
    PRINT *

    THIS% npart= 0
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
                                  > thres_baryon_density )THEN

            THIS% npart= THIS% npart + 1
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

    PRINT *, " * Particles placed. Number of particles=", &
             THIS% npart, "=", DBLE(THIS% npart)/DBLE(THIS% npart_temp), &
             " of the points in lattice."
    PRINT *

    !
    !-- Computing total volume and volume per particle
    !
    THIS% vol  = (xmax - xmin)*(ymax - ymin)*2*ABS(zlim)
    THIS% vol_a= THIS% vol/THIS% npart_temp

    ! Consistency check for the volume
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

    INTEGER:: ix, iy, iz, nx2, ny2, nz2, sgn, npart_half, npart_half2

    DOUBLE PRECISION:: dx, dy, dz
    DOUBLE PRECISION:: xtemp, ytemp, ztemp, zlim, zlim2
    DOUBLE PRECISION:: max_baryon_density1, thres_baryon_density1
    DOUBLE PRECISION:: max_baryon_density2, thres_baryon_density2
    ! Variable used to compute the volume of a particle in an alternative way
    ! to perform a consistency check
    DOUBLE PRECISION:: vol_a_alt

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

    !
    !-- Compute lattices' steps
    !
    dx= ABS(xmax1 - xmin1)/DBLE( THIS% nx )
    dy= ABS(ymax1 - ymin1)/DBLE( THIS% ny )
    dz= ABS(zlim)/DBLE( THIS% nz/2 )

    !
    !-- Compute the number of particles for the second lattice,
    !-- using the same lattice step
    !
    nx2= NINT( ABS(xmax2 - xmin2)/dx )
    ny2= NINT( ABS(ymax2 - ymin2)/dy )
    nz2= NINT( 2*ABS(zlim2)/dz )

    ! Set the number of particles in the z direction to an even number
    ! since half of the particles are above the xy plane, and half below it
    IF( MOD( nz2, 2 ) /= 0 )THEN
      nz2= nz2 - 1
    ENDIF

    PRINT *, " * dx=", dx,  ", dy=", dx,  ", dz=", dz
    PRINT *, " * nx=", THIS% nx,  ", ny=", THIS% ny,  ", nz=", THIS% nz
    PRINT *, " * nx2=", nx2, ", ny2=", ny2, ", nz2=", nz2
    PRINT *

    ! Compute number of lattice points (temporary particle number)
    THIS% npart1_temp = THIS% nx*THIS% ny*THIS% nz !+ THIS% nx*THIS% ny
    THIS% npart2_temp = nx2*ny2*nz2 !+ nx2*ny2
    THIS% npart_temp  = THIS% npart1_temp + THIS% npart2_temp

    PRINT *, " * Number of points for lattice 1= nx1*ny1*nz1=", &
             THIS% npart1_temp
    PRINT *, " * Number of points for lattice 2= nx2*ny2*nz2=", &
             THIS% npart2_temp
    PRINT *

    !
    !-- Compute the mass density almost at the center of the stars
    !
    max_baryon_density1= bns_obj% import_mass_density( &
                                              ( xmin1 + xmax1 )/2.0D0, &
                                              DBLE(0), DBLE(0) )
    max_baryon_density2= bns_obj% import_mass_density( &
                                              ( xmin2 + xmax2 )/2.0D0, &
                                              DBLE(0), DBLE(0) )

    !
    !-- Set the thresholds above which a lattice point is
    !-- promoted to a particle
    !
    thres_baryon_density1= max_baryon_density1/thres
    thres_baryon_density2= max_baryon_density2/thres

    ! Allocating the memory for the array pos( 3, npart_temp )
    ! Note that after determining npart, the array pos is reshaped into
    ! pos( 3, npart )
    IF(.NOT.ALLOCATED( THIS% pos ))THEN
      ALLOCATE( THIS% pos( 3, THIS% npart_temp ), STAT= ios, &
                ERRMSG= err_msg )
      CALL test_status( ios, err_msg, &
                      "...allocation error for array pos in SUBROUTINE" &
                      // "place_particles_3D_lattice." )
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

      ztemp= sgn*( dz/2 + ( iz - 1 )*dz )

      particle_pos_y1: DO iy= 1, THIS% ny, 1

        ytemp= ymin1 + dy/2 + ( iy - 1 )*dy

        particle_pos_x1: DO ix= 1, THIS% nx, 1

          xtemp= xmin1 + dx/2 + ( ix - 1 )*dx

          !
          !-- Promote a lattice point to a particle,
          !-- if the mass density is higher than the threshold
          !
          IF( bns_obj% import_mass_density( xtemp, ytemp, ztemp ) &
                                > thres_baryon_density1 )THEN

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
    particle_pos_z2: DO iz= 1, nz2/2, 1

      ztemp= sgn*( dz/2 + ( iz - 1 )*dz )

      particle_pos_y2: DO iy= 1, ny2, 1

        ytemp= ymin2 + dy/2 + ( iy - 1 )*dy

        particle_pos_x2: DO ix= 1, nx2, 1

          xtemp= xmin2 + dx/2 + ( ix - 1 )*dx

          IF( bns_obj% import_mass_density( xtemp, ytemp, ztemp ) &
                                  > thres_baryon_density2 )THEN

            THIS% npart = THIS% npart + 1
            THIS% npart2= THIS% npart2 + 1
            THIS% pos( 1, THIS% npart )= xtemp
            THIS% pos( 2, THIS% npart )= ytemp
            THIS% pos( 3, THIS% npart )= ztemp

          ENDIF

          ! Print progress on screen, every 10%
          perc= 50*( nx2*ny2*( iz - 1 ) + nx2*( iy - 1 ) + ix )&
                /( nx2*ny2*nz2/2 )
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

    PRINT *, " * Particles placed. Number of particles=", &
             THIS% npart, "=", DBLE(THIS% npart)/DBLE(THIS% npart_temp), &
             " of the points in lattices."
    PRINT *

    !
    !-- Computing total volume and volume per particle
    !
    THIS% vol_a= dx*dy*dz
    THIS% vol1 = (xmax1 - xmin1)*(ymax1 - ymin1)*2*ABS(zlim)
    THIS% vol2 = THIS% npart2_temp * THIS% vol_a
    THIS% vol  = THIS% vol1 + THIS% vol2
    vol_a_alt  = THIS% vol/THIS% npart_temp

    ! Consistency check
    IF( ABS( THIS% vol_a - vol_a_alt ) > 1D-15 )THEN
      PRINT *, " * The particle volume vol_a_alt=", vol_a_alt, "Msun_geo^3"
      PRINT *, " is not equal to dx*dy*dz=", THIS% vol_a, "Msun_geo^3."
      PRINT *
      STOP
    ENDIF

    PRINT *, " * Total volume of the lattices=", THIS% vol, "Msun_geo^3"
    PRINT *, " * Particle volume=", THIS% vol_a, "Msun_geo^3"
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
    CALL test_status( ios, err_msg, "...error when opening " &
             // TRIM(namefile) )

    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Points where some of the hydro fields are negative. "
    CALL test_status( ios, err_msg, "...error when writing line 1 in "&
             // TRIM(namefile) )
    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# column:      1        2       3"
    CALL test_status( ios, err_msg, "...error when writing line 2 in "&
            // TRIM(namefile) )
    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "#      x   y   z"
    CALL test_status( ios, err_msg, "...error when writing line 3 in "&
            // TRIM(namefile) )

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

        CALL test_status( ios, err_msg, "...error in writing "&
                        // "the arrays in " // TRIM(namefile) )
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
      CALL test_status( ios, err_msg, &
                  "...allocation error for array pos" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% lapse_parts ))THEN
      ALLOCATE( THIS% lapse_parts( THIS% npart ), STAT= ios, &
            ERRMSG= err_msg )
      CALL test_status( ios, err_msg, &
                  "...allocation error for array lapse_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% shift_parts_x ))THEN
      ALLOCATE( THIS% shift_parts_x( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      CALL test_status( ios, err_msg, &
                      "...allocation error for shift_parts_x" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% shift_parts_y ))THEN
      ALLOCATE( THIS% shift_parts_y( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      CALL test_status( ios, err_msg, &
                      "...allocation error for shift_parts_y" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% shift_parts_z ))THEN
      ALLOCATE( THIS% shift_parts_z( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      CALL test_status( ios, err_msg, &
                      "...allocation error for shift_parts_z" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_xx_parts ))THEN
      ALLOCATE( THIS% g_xx_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      CALL test_status( ios, err_msg, &
                      "...allocation error for array g_xx_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_xy_parts ))THEN
      ALLOCATE( THIS% g_xy_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      CALL test_status( ios, err_msg, &
                      "...allocation error for array g_xy_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_xz_parts ))THEN
      ALLOCATE( THIS% g_xz_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      CALL test_status( ios, err_msg, &
                      "...allocation error for array g_xz_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_yy_parts ))THEN
      ALLOCATE( THIS% g_yy_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      CALL test_status( ios, err_msg, &
                      "...allocation error for array g_yy_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_yz_parts ))THEN
      ALLOCATE( THIS% g_yz_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      CALL test_status( ios, err_msg, &
                      "...allocation error for array g_yz_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_zz_parts ))THEN
      ALLOCATE( THIS% g_zz_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      CALL test_status( ios, err_msg, &
                      "...allocation error for array g_zz_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% baryon_density_parts ))THEN
      ALLOCATE( THIS% baryon_density_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      CALL test_status( ios, err_msg, &
                   "...allocation error for array baryon_density_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% energy_density_parts ))THEN
      ALLOCATE( THIS% energy_density_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      CALL test_status( ios, err_msg, &
                   "...allocation error for array energy_density_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% specific_energy_parts ))THEN
      ALLOCATE( THIS% specific_energy_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      CALL test_status( ios, err_msg, &
                  "...allocation error for array specific_energy_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% pressure_parts ))THEN
      ALLOCATE( THIS% pressure_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      CALL test_status( ios, err_msg, &
                  "...allocation error for array pressure_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% pressure_parts_cu ))THEN
      ALLOCATE( THIS% pressure_parts_cu( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      CALL test_status( ios, err_msg, &
                  "...allocation error for array pressure_parts_cu" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% v_euler_parts_x ))THEN
      ALLOCATE( THIS% v_euler_parts_x( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      CALL test_status( ios, err_msg, &
                      "...allocation error for array v_euler_parts_x" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% v_euler_parts_y ))THEN
      ALLOCATE( THIS% v_euler_parts_y( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      CALL test_status( ios, err_msg, &
                      "...allocation error for array v_euler_parts_y" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% v_euler_parts_z ))THEN
      ALLOCATE( THIS% v_euler_parts_z( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      CALL test_status( ios, err_msg, &
                      "...allocation error for array v_euler_parts_z" )
    ENDIF

    PRINT *, "** Subroutine allocate_lorene_id_parts_memory executed."
    PRINT *

  END PROCEDURE allocate_lorene_id_parts_memory

END SUBMODULE particles_constructor
