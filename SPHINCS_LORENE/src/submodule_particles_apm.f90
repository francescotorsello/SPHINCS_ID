! File:         submodule_particles_apm.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (particles_id) particles_apm

  !****************************************************
  !                                                   *
  ! Implementation of the method                      *
  ! compute_and_export_SPH_variables_apm              *
  ! of TYPE particles.                                *
  !                                                   *
  ! FT 04.06.2021                                     *
  !                                                   *
  !****************************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE compute_and_export_SPH_variables_apm

    !****************************************************
    !                                                   *
    ! Compute  the particle positions as follows:       *
    !                                                   *
    !   1. Take initial particle distribution as input  *
    !   2. Assume that the particles have the same mass *
    !   3. Do the APM iteration so that the final       *
    !      particle number density matches the baryon   *
    !      density in the star                          *
    !   4. Correct the particle masses ONCE in order    *
    !      to match the density even better. Since we   *
    !      don't want a large mass ratio, we impose a   *
    !      maximum mass ratio when performing this      *
    !      correction.                                  *
    !                                                   *
    ! After this procedure, the resulting particle      *
    ! distribution has positions and baryon numbers     *
    ! that kernel-estimate very well the mass density   *
    ! of the star, and has a low mass ratio.            *
    !                                                   *
    ! This procedure assigns positions and nu. After it *
    ! is performed, the other SPH quantities are        *
    ! computed, and then they are printed to a binary   *
    ! file ready to be used by SPHINCS_BSSN.            *
    !                                                   *
    ! FT 04.06.2021                                     *
    !                                                   *
    !****************************************************

    USE constants, ONLY: third

    IMPLICIT NONE

    INTEGER, PARAMETER:: max_npart   = 5D+6
    LOGICAL, PARAMETER:: debug= .TRUE.

    INTEGER:: a, itr            ! iterators
    INTEGER:: npart_real, npart_real_half, npart_ghost
    INTEGER:: nx, ny, nz, i, j, k

    DOUBLE PRECISION:: smaller_radius, larger_radius, radius_y, radius_z
    DOUBLE PRECISION:: h_max, h_av, eps, delta
    DOUBLE PRECISION:: xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, &
                       rad_x, rad_y, rad_z
    DOUBLE PRECISION:: xtemp, ytemp, ztemp, x_ell, y_ell, z_ell

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos_tmp
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: ghost_pos

    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: h_guess
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: h_tmp

    LOGICAL:: exist

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    IF( debug ) PRINT *, "0"

    npart_real= SIZE( pos_input(1,:) )

    IF(.NOT.ALLOCATED( h_guess ))THEN
      ALLOCATE( h_guess( npart_real ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array h_guess in SUBROUTINE ", &
                  "compute_and_export_SPH_variables_apm. The error message is",&
                  err_msg
         STOP
      ENDIF
    ENDIF

    DO a= 1, npart_real, 1
      h_guess(a)= pvol(a)**third
    ENDDO

    !--------------------------------------------------------------------!
    !-- Store particles above xy plane as the first half of the array, --!
    !-- and mirror them to the second half                             --!
    !--------------------------------------------------------------------!

    pos_tmp= pos_input
    h_tmp= h_guess
    itr= 0

    DO a= 1, npart_real, 1
      IF( pos_tmp( 3, a ) > 0.0D0 )THEN
        itr= itr + 1
        pos_input( 1, itr )= pos_tmp( 1, a )
        pos_input( 2, itr )= pos_tmp( 2, a )
        pos_input( 3, itr )= pos_tmp( 3, a )
        h_guess( itr )     = h_tmp( itr )
      ENDIF
    ENDDO
    npart_real_half= itr

    DO a= 1, npart_real_half, 1
      pos_input( 1, npart_real_half + a )=   pos_input( 1, a )
      pos_input( 2, npart_real_half + a )=   pos_input( 2, a )
      pos_input( 3, npart_real_half + a )= - pos_input( 3, a )
      h_guess( npart_real_half + a )     =   h_guess( a )
    ENDDO
    npart_real= 2*npart_real_half

    IF( debug ) PRINT *, "1"

    !--------------------------------------------------------------------!
    !-- Find the maximum and the average smoothing length of the       --!
    !-- particles whose distance from the center is higher than        --!
    !-- radius_z, and use them to place ghost particles a little more  --!
    !-- outside than the surface of the particles.                     --!
    !--------------------------------------------------------------------!

    smaller_radius= ABS( MINVAL( pos_input( 1, : ), DIM= 1 ) - center )
    larger_radius = ABS( center - MAXVAL( pos_input( 1, : ), DIM= 1 ) )
    radius_y= ABS( MAXVAL( pos_input( 2, : ), DIM= 1 ) )
    radius_z= ABS( MAXVAL( pos_input( 3, : ), DIM= 1 ) )

    h_max= 0.0D0
    h_av = 0.0D0
    itr  = 0
    DO a= 1, npart_real, 1

      IF( SQRT( ( pos_input( 1, a ) - center )**2.0D0 &
                + pos_input( 2, a )**2.0D0 &
                + pos_input( 3, a )**2.0D0 ) > 0.99D0*radius_z )THEN

        itr= itr + 1
        IF( h_guess(a) > h_max )THEN
          h_max= h_guess(a)
        ENDIF
        h_av= h_av + h_guess(a)

      ENDIF

    ENDDO
    h_av= h_av/itr

    IF( debug ) PRINT *, "2"

    !-------------------------------!
    !--  Placing ghost particles  --!
    !-------------------------------!

    IF(.NOT.ALLOCATED( ghost_pos ))THEN
      ALLOCATE( ghost_pos( 3, max_npart ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array ghost_pos in SUBROUTINE ", &
                  "compute_and_export_SPH_variables_apm. The error message is",&
                  err_msg
         STOP
      ENDIF
    ENDIF

    ghost_pos= 0.0D0

    PRINT *, " * Placing ghost particles on a lattice between ellipsodial ", &
             "surfaces..."
    PRINT *

    nx= 150
    ny= 150
    nz= 150
    eps= 5.0D-1
    xmin= center - larger_radius*( 1.0D0 + eps )
    xmax= center + larger_radius*( 1.0D0 + eps )
    ymin= - radius_y*( 1.0D0 + eps )
    ymax=   radius_y*( 1.0D0 + eps )
    zmin= - radius_z*( 1.0D0 + eps )
    zmax=   radius_z*( 1.0D0 + eps )
    dx= ABS( xmax - xmin )/DBLE( nx )
    dy= ABS( ymax - ymin )/DBLE( ny )
    dz= ABS( zmax - zmin )/DBLE( nz )
    delta= 0.25D0

    rad_x= larger_radius + h_av/1.0D0
    rad_y= radius_y + h_av/1.0D0
    rad_z= radius_z + h_av/1.0D0

    itr= 0
    DO k= 1, nz, 1

      ztemp= zmin + dz/2 + ( k - 1 )*dz

      DO j= 1, ny, 1

        ytemp= ymin + dy/2 + ( j - 1 )*dy

        DO i= 1, nx, 1

          xtemp= xmin + dx/2 + ( i - 1 )*dx

          x_ell= center + rad_x*COS(ATAN( ytemp/xtemp )) &
                 *SIN(ACOS(ztemp/SQRT( ( xtemp - center )**2.0D0 &
                                       + ytemp**2.0D0 + ztemp**2.0D0 )))

          y_ell= rad_y*SIN(ATAN( ytemp/xtemp )) &
                 *SIN(ACOS(ztemp/SQRT( ( xtemp - center )**2.0D0 &
                                       + ytemp**2.0D0 + ztemp**2.0D0 )))

          z_ell= rad_z*( ztemp/SQRT( ( xtemp - center )**2.0D0 &
                                     + ytemp**2.0D0 + ztemp**2.0D0 ))

          IF( binary% import_mass_density( xtemp, ytemp, ztemp ) <= 0.0D0 &
              .AND. &
              ! TODO: understand why delta is needed...
              SQRT( ( xtemp - center )**2.0D0 + ytemp**2.0D0 &
                    + ztemp**2.0D0 ) <= &
                    1.1D0*SQRT( ( x_ell - center )**2.0D0 &
                                + delta*y_ell**2.0D0 + z_ell**2.0D0 ) &
              .AND. &
              SQRT( ( xtemp - center )**2.0D0 + ytemp**2.0D0 &
                    + ztemp**2.0D0 ) >= &
              SQRT( ( x_ell - center )**2.0D0 + delta*y_ell**2.0D0 &
                    + z_ell**2.0D0 ) &
          )THEN

            itr= itr + 1
            ghost_pos( 1, itr )= xtemp
            ghost_pos( 2, itr )= ytemp
            ghost_pos( 3, itr )= ztemp

          ENDIF

         ENDDO
      ENDDO
    ENDDO
    npart_ghost= itr
    IF( npart_ghost == 0 )THEN
      PRINT *, "** ERROR: No ghost particles were placed. Stopping.."
      PRINT *
      STOP
    ENDIF
    ghost_pos = ghost_pos( :, 1:npart_ghost )

    PRINT *, " * ", npart_ghost, " ghost particles placed around ", &
             npart_real, "real particles."
    PRINT *

    PRINT *, " * Printing ghost particles to file..."

    finalnamefile= "ghost_pos.dat"

    INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

    IF( exist )THEN
        OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
              FORM= "FORMATTED", &
              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
              IOMSG= err_msg )
    ELSE
        OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "NEW", &
              FORM= "FORMATTED", &
              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening " // TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF

    DO a= 1, npart_real, 1
      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        1, a, &
        pos_input( 1, a ), &
        pos_input( 2, a ), &
        pos_input( 3, a )
    ENDDO

    DO a= 1, npart_ghost, 1
      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        2, a, &
        ghost_pos( 1, a ), &
        ghost_pos( 2, a ), &
        ghost_pos( 3, a )
    ENDDO

    CLOSE( UNIT= 2 )

    PRINT *, " * Printed."

    STOP

  END PROCEDURE compute_and_export_SPH_variables_apm


END SUBMODULE particles_apm
