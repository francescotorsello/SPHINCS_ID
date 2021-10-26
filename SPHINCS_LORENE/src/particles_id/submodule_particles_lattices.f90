! File:         submodule_particles_lattices.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (particles_id) particles_lattices

  !***************************************************
  !
  !# This SUBMODULE contains the implementation of
  !  the methods of TYPE particles
  !  that place particles on 1 or 2 lattices around
  !  the stars.
  !
  !  FT 12.07.2021
  !
  !***************************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE place_particles_lattice

    !*********************************************************
    !
    !# Places paricles on a 3D lattice containing both stars
    !
    !  FT 5.10.2020
    !
    !*********************************************************

    USE constants,    ONLY: Msun_geo, pi, third

    IMPLICIT NONE

    INTEGER:: i, j, k, sgn, nx, ny, nz, npart_half
    INTEGER:: npart_tmp

    DOUBLE PRECISION:: dx, dy, dz
    DOUBLE PRECISION:: xtemp, ytemp, ztemp, zlim
    DOUBLE PRECISION:: thres_baryon_density
    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: pos_tmp

    PRINT *, "** Executing the place_particles_lattice " &
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
    !-- Compute number of lattice points (for now, equal in each direction)
    !
    nx= CEILING(stretch*(6.0D0*npart_des/pi)**third)
    ny= nx
    nz= nx
    PRINT *, " * nx= ny= nz=", nx
    PRINT *

    !
    !-- Compute lattice steps
    !
    dx= ABS(xmax - xmin)/DBLE( nx )
    dy= ABS(ymax - ymin)/DBLE( ny )
    dz= ABS(zlim)/DBLE( nz/2 )

    PRINT *, " * dx=", dx,  ", dy=", dx,  ", dz=", dz
    PRINT *

    npart_tmp = nx*ny*nz

    PRINT *, " * Number of lattice points= nx*ny*nz=", npart_tmp
    PRINT *

    !
    !-- Set the threshold above which a lattice point is
    !-- promoted to a particle
    !
    IF( THIS% use_thres )THEN
      thres_baryon_density= central_density/thres
    ELSE
      thres_baryon_density= 0.0D0
    ENDIF

    ! Allocating the memory for the array pos( 3, npart_tmp )
    ! Note that after determining npart, the array pos is reshaped into
    ! pos( 3, npart )
    IF(.NOT.ALLOCATED( THIS% pos ))THEN
      ALLOCATE( THIS% pos( 3, npart_tmp ), STAT= ios, &
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

    IF(.NOT.ALLOCATED( pos_tmp ))THEN
      ALLOCATE( pos_tmp( 3, nx, ny, nz ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array pos_tmp in SUBROUTINE" &
                  // "place_particles_3D_lattice. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array pos in SUBROUTINE" &
      !                // "place_particles_3D_lattice." )
    ENDIF
    ! Initializing the array pos to 0
    pos_tmp= HUGE(0.0D0)

    !---------------------------------------------------------!
    !--  Storing the particle positions into the array pos  --!
    !--  symmetrically w.r.t. the xy plane                  --!
    !---------------------------------------------------------!

    PRINT *, " * Placing particles on the lattice..."
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
    !-- Place the first half of the particle (above or below the xy plane,
    !-- depending on the variable sgn)
    !
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( nx, ny, nz, dx, dy, dz, sgn, &
    !$OMP                     pos_tmp, thres_baryon_density, xmin, ymin ) &
    !$OMP             PRIVATE( i, j, k, xtemp, ytemp, ztemp )
    particle_pos_z: DO k= 1, nz/2, 1

      ztemp= sgn*( dz/2 + ( k - 1 )*dz )

      particle_pos_y: DO j= 1, ny, 1

        ytemp= ymin + ( j - 1 )*dy

        particle_pos_x: DO i= 1, nx, 1

          xtemp= xmin + dx/2 + ( i - 1 )*dx

          !
          !-- Promote a lattice point to a particle,
          !-- if the mass density is higher than the threshold
          !
          IF( get_density( xtemp, ytemp, ztemp ) &
                                  > thres_baryon_density &
              .AND. &
              validate_position( xtemp, ytemp, ztemp ) == 0 )THEN

            pos_tmp( 1, i, j, k )= xtemp
            pos_tmp( 2, i, j, k )= ytemp
            pos_tmp( 3, i, j, k )= ztemp

          ENDIF

         ENDDO particle_pos_x
      ENDDO particle_pos_y
    ENDDO particle_pos_z
    !$OMP END PARALLEL DO

    npart_out= 0
    DO k= 1, nz, 1

      DO j= 1, ny, 1

        DO i= 1, nx, 1

          IF( pos_tmp( 1, i, j, k ) < HUGE(0.0D0) )THEN

            npart_out= npart_out + 1

            THIS% pos( 1, THIS% npart )= pos_tmp( 1, i, j, k )
            THIS% pos( 2, THIS% npart )= pos_tmp( 2, i, j, k )
            THIS% pos( 3, THIS% npart )= pos_tmp( 3, i, j, k )

          ENDIF

         ENDDO
      ENDDO
    ENDDO
    npart_half= npart_out
    IF( npart_half == 0 )THEN
      PRINT *, "** There are no particles! Execution stopped..."
      PRINT *
      STOP
    ENDIF

    DEALLOCATE( pos_tmp )

    !
    !-- Place the second half of the particles, mirroring the first half
    !-- w.r.t the xy plane
    !
    particle_pos_z_mirror: DO k= 1, npart_half, 1

      xtemp=   THIS% pos( 1, k )
      ytemp=   THIS% pos( 2, k )
      ztemp= - THIS% pos( 3, k )

      npart_out= npart_out + 1
      THIS% pos( 1, THIS% npart )= xtemp
      THIS% pos( 2, THIS% npart )= ytemp
      THIS% pos( 3, THIS% npart )= ztemp

      !ENDIF

      ! Print progress on screen, every 10%
      !perc= 50 + 50*k/( npart_half )
      !IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
      !   WRITE( *, "(A2,I3,A1)", ADVANCE= "NO" ) &
      !           creturn//" ", perc, "%"
      !ENDIF
    ENDDO particle_pos_z_mirror
    !WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

    !
    !-- Consistency checks
    !
    IF( npart_out /= 2*npart_half )THEN
      PRINT *
      PRINT *, "** ERROR: The number of particles ", npart_out, &
               " is not the expected value ", 2*npart_half
      PRINT *
      STOP
    ENDIF

    DO k= 1, npart_half, 1
      IF( THIS% pos( 3, k ) /= - THIS% pos( 3, npart_half + k ) )THEN
        PRINT *
        PRINT *, "** ERROR: The lattice is not mirrored " &
                 // "by the xy plane."
        PRINT *
        STOP
      ENDIF
    ENDDO

    PRINT *, " * Particles placed. Number of particles=", &
             npart_out, "=", DBLE(npart_out)/DBLE(npart_tmp), &
             " of the points in lattice."
    PRINT *

    !
    !-- Computing total volume and volume per particle
    !
    IF(.NOT.ALLOCATED( pvol ))THEN
      ALLOCATE( pvol( npart_out ), STAT= ios, &
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
    THIS% vol_a= THIS% vol/npart_tmp

    pvol= THIS% vol_a

    ! Consistency check for the particle volume
    IF( ABS( THIS% vol_a - dx*dy*dz ) > 1.0D-9 )THEN
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

  END PROCEDURE place_particles_lattice


  MODULE PROCEDURE place_particles_lattices

    !****************************************************
    !
    !# Places paricles on two 3D lattices, each one
    !  containing a star
    !
    !  FT 19.10.2020
    !
    !****************************************************

    USE constants, ONLY: Msun_geo

    IMPLICIT NONE

    INTEGER:: i, j, k, sgn, npart_half, npart_half2
    INTEGER:: npart_tmp, npart1_temp, npart2_temp
    INTEGER:: nx1, ny1, nz1, nx2, ny2, nz2

    DOUBLE PRECISION:: dx1, dy1, dz1, dx2, dy2, dz2
    DOUBLE PRECISION:: xtemp, ytemp, ztemp, zlim, zlim2
    DOUBLE PRECISION:: max_baryon_density1, thres_baryon_density1
    DOUBLE PRECISION:: max_baryon_density2, thres_baryon_density2
    ! Variable used to compute the volume of a particle in an alternative way
    ! to perform a consistency check
    DOUBLE PRECISION:: vol_a_alt1, vol_a_alt2

    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: pos_tmp

    PRINT *, "** Executing the place_particles_lattices " &
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

    THIS% mass1= bnsb% get_mass1()
    THIS% mass2= bnsb% get_mass2()

    IF( THIS% mass1 > THIS% mass2 )THEN

      ! mass_ratio < 1
      THIS% mass_ratio= THIS% mass2/THIS% mass1
      !
      !-- Compute lattices' steps
      !
      nx2= nx
      ny2= ny
      nz2= nz
      dx2= ABS(xmax2 - xmin2)/DBLE( nx2 )
      dy2= ABS(ymax2 - ymin2)/DBLE( ny2 )
      dz2= ABS(zlim2)/DBLE( nz2/2 )

      dx1= dx2*(THIS% mass_ratio**(1.0D0/3.0D0))
      dy1= dy2*(THIS% mass_ratio**(1.0D0/3.0D0))
      dz1= dz2*(THIS% mass_ratio**(1.0D0/3.0D0))
      nx1= NINT( ABS(xmax1 - xmin1)/dx1 ) + 1
      ny1= NINT( ABS(ymax1 - ymin1)/dy1 ) + 1
      nz1= NINT( 2*ABS(zlim)/dz1 ) + 1

    ELSE

      ! mass_ratio < 1
      THIS% mass_ratio= THIS% mass1/THIS% mass2
      !
      !-- Compute lattices' steps
      !
      nx1= nx
      ny1= ny
      nz1= nz
      dx1= ABS(xmax1 - xmin1)/DBLE( nx1 )
      dy1= ABS(ymax1 - ymin1)/DBLE( ny1 )
      dz1= ABS(zlim)/DBLE( nz1/2 )

      dx2= dx1*(THIS% mass_ratio**(1.0D0/3.0D0))
      dy2= dy1*(THIS% mass_ratio**(1.0D0/3.0D0))
      dz2= dz1*(THIS% mass_ratio**(1.0D0/3.0D0))
      nx2= NINT( ABS(xmax2 - xmin2)/dx2 ) + 1
      ny2= NINT( ABS(ymax2 - ymin2)/dy2 ) + 1
      nz2= NINT( 2*ABS(zlim2)/dz2 ) + 1

    ENDIF

    ! Set the number of particles in the z direction to an even number
    ! since half of the particles are above the xy plane, and half below it
    IF( MOD( nz2, 2 ) /= 0 )THEN
      nz2= nz2 - 1
    ENDIF

    PRINT *, " * dx1=", dx1,  ", dy1=", dx1,  ", dz1=", dz1
    PRINT *, " * dx2=", dx2,  ", dy2=", dx2,  ", dz2=", dz2
    PRINT *, " * nx1=", nx1, ", ny1=", ny1, ", nz1=", nz1
    PRINT *, " * nx2=", nx2, ", ny2=", ny2, ", nz2=", nz2
    PRINT *

    !PRINT *, " * xmin1=", xmin1, ", xmax1=", xmax1
    !PRINT *, " * xmin2=", xmin2, ", xmax2=", xmax2
    !PRINT *
    !STOP

    ! Compute number of lattice points (temporary particle number)
    npart1_temp = nx*ny*nz !+ nx*ny
    npart2_temp = nx2*ny2*nz2 !+ nx2*ny2
    npart_tmp  = npart1_temp + npart2_temp

    PRINT *, " * Number of points for lattice 1= nx1*ny1*nz1=", &
             npart1_temp
    PRINT *, " * Number of points for lattice 2= nx2*ny2*nz2=", &
             npart2_temp
    PRINT *

    !
    !-- Compute the mass density at the center of the stars
    !

    ! The following two density ar in SPHINCS units [Msun Msun_geo^{-3}]
    max_baryon_density1= bnsb% get_rho_center1()
    max_baryon_density2= bnsb% get_rho_center2()

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

    ! Allocating the memory for the array pos( 3, npart_tmp )
    ! Note that after determining npart, the array pos is reshaped into
    ! pos( 3, npart )
    IF(.NOT.ALLOCATED( THIS% pos ))THEN
      ALLOCATE( THIS% pos( 3, npart_tmp ), STAT= ios, &
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

    IF(.NOT.ALLOCATED( pos_tmp ))THEN
      ALLOCATE( pos_tmp( 3, nx, ny, nz ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array pos_tmp in SUBROUTINE" &
                  // "place_particles_3D_lattices. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array pos in SUBROUTINE" &
      !                // "place_particles_3D_lattice." )
    ENDIF
    ! Initializing the array pos to 0
    pos_tmp= HUGE(0.0D0)

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
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( nx, ny, nz, bnsb, dx1, dy1, dz1, sgn, &
    !$OMP                     pos_tmp, thres_baryon_density1, xmin1, ymin1 ) &
    !$OMP             PRIVATE( i, j, k, xtemp, ytemp, ztemp )
    particle_pos_z1: DO k= 1, nz/2, 1

      ztemp= sgn*( dz1/2 + ( k - 1 )*dz1 )

      particle_pos_y1: DO j= 1, ny, 1

        ytemp= ymin1 + dy1/2 + ( j - 1 )*dy1

        particle_pos_x1: DO i= 1, nx, 1

          xtemp= xmin1 + dx1/2 + ( i - 1 )*dx1

          !
          !-- Promote a lattice point to a particle,
          !-- if the mass density is higher than the threshold
          !
          IF( bnsb% read_mass_density( xtemp, ytemp, ztemp ) &
                                > thres_baryon_density1 &
              .AND. &
              bnsb% test_position( xtemp, ytemp, ztemp ) == 0 )THEN

            !THIS% npart = THIS% npart + 1
            !THIS% npart1= THIS% npart1 + 1
            !THIS% pos( 1, THIS% npart )= xtemp
            !THIS% pos( 2, THIS% npart )= ytemp
            !THIS% pos( 3, THIS% npart )= ztemp
            pos_tmp( 1, i, j, k )= xtemp
            pos_tmp( 2, i, j, k )= ytemp
            pos_tmp( 3, i, j, k )= ztemp

          ENDIF

          ! Print progress on screen, every 10%
          !perc= 50*( nx*ny*k + nx*j + i )/ &
          !        ( nx*ny*nz/2 )
          !IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
          !  WRITE( *, "(A2,I3,A1)", ADVANCE= "NO" ) &
          !         creturn//" ", perc, "%"
          ! ENDIF
        ENDDO particle_pos_x1
      ENDDO particle_pos_y1
    ENDDO particle_pos_z1
    !$OMP END PARALLEL DO
    !WRITE( *, "(A1)", ADVANCE= "NO" ) creturn


    DO k= 1, nz, 1

      DO j= 1, ny, 1

        DO i= 1, nx, 1

          IF( pos_tmp( 1, i, j, k ) < HUGE(0.0D0) )THEN

            THIS% npart= THIS% npart + 1
            THIS% npart1= THIS% npart1 + 1
            THIS% pos( 1, THIS% npart )= pos_tmp( 1, i, j, k )
            THIS% pos( 2, THIS% npart )= pos_tmp( 2, i, j, k )
            THIS% pos( 3, THIS% npart )= pos_tmp( 3, i, j, k )

          ENDIF

         ENDDO
      ENDDO
    ENDDO
    npart_half= THIS% npart
    IF( npart_half == 0 )THEN
      PRINT *, "** There are no particles on star 1! Execution stopped..."
      PRINT *
      STOP
    ENDIF

    DEALLOCATE( pos_tmp )

    !
    !-- Place the second half of the particles, mirroring the first half
    !-- w.r.t the xy plane
    !
    particle_pos_z1_mirror: DO k= 1, npart_half, 1

      xtemp=   THIS% pos( 1, k )
      ytemp=   THIS% pos( 2, k )
      ztemp= - THIS% pos( 3, k )

      ! TODO: is this check needed?
      !IF( import_mass_density( xtemp, ytemp, ztemp ) &
      !                               > thres_baryon_density1 )THEN

      THIS% npart = THIS% npart + 1
      THIS% npart1= THIS% npart1 + 1
      THIS% pos( 1, THIS% npart )= xtemp
      THIS% pos( 2, THIS% npart )= ytemp
      THIS% pos( 3, THIS% npart )= ztemp

      !ENDIF

      perc= 50 + 50*k/( npart_half )
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
    !particle_pos_y1_xy: DO j= 1, ny, 1
    !
    !  ytemp= ymin1 + dy/2 + ( j - 1 )*dy
    !
    !  particle_pos_x1_xy: DO i= 1, nx, 1
    !
    !    xtemp= xmin1 + dx/2 + ( i - 1 )*dx
    !
    !    !
    !    !-- Promote a lattice point to a particle,
    !    !-- if the mass density is higher than the threshold
    !    !
    !    IF( bnsb% import_mass_density( xtemp, ytemp, ztemp ) &
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
    !    perc= 50*( nx*ny*k + nx*j + i )/ &
    !            ( nx*ny*nz/2 )
    !    IF( MOD( perc, 10 ) == 0 )THEN
    !      WRITE( *, "(A2,I3,A1)", ADVANCE= "NO" ) &
    !             creturn//" ", perc, "%"
    !     ENDIF
    !   ENDDO particle_pos_x1_xy
    !ENDDO particle_pos_y1_xy
    !WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

    IF(.NOT.ALLOCATED( pos_tmp ))THEN
      ALLOCATE( pos_tmp( 3, nx2, ny2, nz2 ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array pos_tmp in SUBROUTINE" &
                  // "place_particles_3D_lattices. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array pos in SUBROUTINE" &
      !                // "place_particles_3D_lattice." )
    ENDIF
    ! Initializing the array pos to 0
    pos_tmp= HUGE(0.0D0)

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
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( nx2, ny2, nz2, bnsb, dx2, dy2, dz2, sgn, &
    !$OMP                     pos_tmp, thres_baryon_density2, xmin2, ymin2 ) &
    !$OMP             PRIVATE( i, j, k, xtemp, ytemp, ztemp )
    particle_pos_z2: DO k= 1, nz2/2, 1

      ztemp= sgn*( dz2/2 + ( k - 1 )*dz2 )

      particle_pos_y2: DO j= 1, ny2, 1

        ytemp= ymin2 + dy2/2 + ( j - 1 )*dy2

        particle_pos_x2: DO i= 1, nx2, 1

          xtemp= xmin2 + dx2/2 + ( i - 1 )*dx2

          IF( bnsb% read_mass_density( xtemp, ytemp, ztemp ) &
                                  > thres_baryon_density2 &
              .AND. &
              bnsb% test_position( xtemp, ytemp, ztemp ) == 0 )THEN

            !THIS% npart = THIS% npart + 1
            !THIS% npart2= THIS% npart2 + 1
            !THIS% pos( 1, THIS% npart )= xtemp
            !THIS% pos( 2, THIS% npart )= ytemp
            !THIS% pos( 3, THIS% npart )= ztemp
            pos_tmp( 1, i, j, k )= xtemp
            pos_tmp( 2, i, j, k )= ytemp
            pos_tmp( 3, i, j, k )= ztemp

          ENDIF

          ! Print progress on screen, every 10%
          !perc= 50*( nx2*ny2*( k - 1 ) + nx2*( j - 1 ) &
          !      + i )/( nx2*ny2*nz2/2 )
          !IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
          !  WRITE( *, "(A2,I3,A1)", ADVANCE= "NO" ) &
          !         creturn//" ", perc, "%"
          !ENDIF
        ENDDO particle_pos_x2
      ENDDO particle_pos_y2
    ENDDO particle_pos_z2
    !$OMP END PARALLEL DO
    !WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

    DO k= 1, nz2, 1

      DO j= 1, ny2, 1

        DO i= 1, nx2, 1

          IF( pos_tmp( 1, i, j, k ) < HUGE(0.0D0) )THEN

            THIS% npart= THIS% npart + 1
            THIS% npart2= THIS% npart2 + 1
            THIS% pos( 1, THIS% npart )= pos_tmp( 1, i, j, k )
            THIS% pos( 2, THIS% npart )= pos_tmp( 2, i, j, k )
            THIS% pos( 3, THIS% npart )= pos_tmp( 3, i, j, k )

          ENDIF

         ENDDO
      ENDDO
    ENDDO
    npart_half2= THIS% npart
    IF( npart_half2 == 2*npart_half )THEN
      PRINT *, "** There are no particles on star 2! Execution stopped..."
      PRINT *
      STOP
    ENDIF

    DEALLOCATE( pos_tmp )

    particle_pos_z2_mirror: DO k= 2*npart_half + 1, npart_half2, 1

      xtemp=   THIS% pos( 1, k )
      ytemp=   THIS% pos( 2, k )
      ztemp= - THIS% pos( 3, k )

      !IF( import_mass_density( xtemp, ytemp, ztemp ) &
      !                               > thres_baryon_density2 )THEN

      THIS% npart = THIS% npart + 1
      THIS% npart2= THIS% npart2 + 1
      THIS% pos( 1, THIS% npart )= xtemp
      THIS% pos( 2, THIS% npart )= ytemp
      THIS% pos( 3, THIS% npart )= ztemp

      !ENDIF

      ! Print progress on screen, every 10%
      perc= 50 + 50*( k - 2*npart_half + 1 ) &
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
    !particle_pos_y2_xy: DO j= 1, ny2, 1
    !
    !  ytemp= ymin2 + dy/2 + ( j - 1 )*dy
    !
    !  particle_pos_x2_xy: DO i= 1, nx2, 1
    !
    !    xtemp= xmin2 + dx/2 + ( i - 1 )*dx
    !
    !    IF( bnsb% import_mass_density( xtemp, ytemp, ztemp ) &
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
    !    perc= 50*( nx2*ny2*( k - 1 ) + nx2*( j - 1 ) + i )&
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

    DO k= 1, npart_half, 1
      IF( THIS% pos( 3, k ) /= - THIS% pos( 3, npart_half + k ) )THEN
        PRINT *
        PRINT *, "** ERROR: The lattice around NS 1 are not mirrored " &
                 // "by the xy plane."
        PRINT *
        STOP
      ENDIF
    ENDDO
    DO k= 2*npart_half + 1, npart_half2, 1
      IF( THIS% pos( 3, k ) /= &
          - THIS% pos( 3, ( npart_half2 - 2*npart_half ) + k ) )THEN
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
             THIS% npart, "=", DBLE(THIS% npart)/DBLE(npart_tmp), &
             " of the points in lattices."
    PRINT *
    PRINT *, " * Number of particles on NS 1=", THIS% npart1, "=", &
             DBLE(THIS% npart1)/DBLE(npart1_temp), &
             " of the points in the first lattice."
    PRINT *, " * Number of particles on NS 2=", THIS% npart2, "=", &
             DBLE(THIS% npart2)/DBLE(npart2_temp), &
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
    !THIS% vol2 = npart2_temp * THIS% vol_a
    !THIS% vol  = THIS% vol1 + THIS% vol2
    vol_a_alt1  = THIS% vol1/npart1_temp

    THIS% vol2_a= dx2*dy2*dz2
    THIS% vol2 =  dx2*nx2*dy2*ny2*dz2*nz2
    !THIS% vol2 = (xmax2 - xmin2)*(ymax2 - ymin2)*2*ABS(zlim2)
    !THIS% vol2 = npart2_temp * THIS% vol_a2
    !THIS% vol  = THIS% vol1 + THIS% vol2
    vol_a_alt2  = THIS% vol2/npart2_temp

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

    PRINT *, "** Subroutine place_particles_lattices " &
             // "executed."
    PRINT *

  END PROCEDURE place_particles_lattices


END SUBMODULE particles_lattices
