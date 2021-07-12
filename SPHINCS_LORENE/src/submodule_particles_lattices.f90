! File:         submodule_particles_lattices.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (particles_id) particles_lattices

  !***************************************************
  !                                                  *
  ! Implementation of the methods of TYPE particles  *
  ! that place particles on 1 or 2 lattices around   *
  ! the stars.                                       *
  !                                                  *
  ! FT 12.07.2021                                    *
  !                                                  *
  !***************************************************


  IMPLICIT NONE


  CONTAINS


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


END SUBMODULE particles_lattices
