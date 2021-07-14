! File:         submodule_particles_apm.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (particles_id) particles_apm

  !****************************************************
  !                                                   *
  ! Implementation of the method                      *
  ! perform_apm of TYPE particles.                    *
  !                                                   *
  ! FT 04.06.2021                                     *
  !                                                   *
  !****************************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE perform_apm

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

    USE constants,           ONLY: half, third, Msun, Msun_geo, km2m, g2kg, &
                                   amu, pi

    USE sph_variables,       ONLY: allocate_sph_memory, deallocate_sph_memory, &
                                   npart, h, nu, Theta, &
                                   divv,av,Pr,ye,temp,nlrf,&
                                   u,tterm,tgrav,tkin,&
                                   escap,t,n1,n2,pos_u,&
                                   npm,Nstar, S_l, ehat, cs, Kent, vel_u
    USE metric_on_particles, ONLY: allocate_metric_on_particles, &
                                   deallocate_metric_on_particles
    USE gradient,            ONLY: allocate_gradient, deallocate_gradient
    !USE alive_flag,          ONLY: alive
    USE set_h,               ONLY: exact_nei_tree_update, posmash
    USE RCB_tree_3D,         ONLY: allocate_RCB_tree_memory_3D, iorig, &
                                   deallocate_RCB_tree_memory_3D
    !USE kernel_table,        ONLY: ktable
    !USE input_output,        ONLY: read_options
    USE units,               ONLY: umass!, set_units
    !USE options,             ONLY: ikernel, ndes

    USE APM,                 ONLY: density_loop, position_correction, assign_h
    USE analyze,             ONLY: COM
    USE matrix,              ONLY: determinant_4x4_matrix

    USE sphincs_sph,         ONLY: density

    IMPLICIT NONE

    INTEGER,          PARAMETER:: max_npart   = 5D+6
    INTEGER,          PARAMETER:: nn_des      = 301
    INTEGER,          PARAMETER:: m_max_it    = 50
    DOUBLE PRECISION, PARAMETER:: tol= 1.0D-3
    DOUBLE PRECISION, PARAMETER:: iter_tol= 2.0D-2
    LOGICAL,          PARAMETER:: debug= .TRUE.

    INTEGER:: a, a2, itr, itr2, n_inc            ! iterators
    INTEGER:: npart_real, npart_real_half, npart_ghost, npart_all, npart_missing
    INTEGER:: nx, ny, nz, i, j, k, nus, mus
    INTEGER:: a_numin, a_numin2, a_numax, a_numax2

    DOUBLE PRECISION:: smaller_radius, larger_radius, radius_y, radius_z
    DOUBLE PRECISION:: h_max, h_av, eps!, delta
    DOUBLE PRECISION:: xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, &
                       rad_x, rad_y, rad_z, com_x, com_y, com_z, com_d
    DOUBLE PRECISION:: max_r_real, min_r_ghost, r_real, r_ghost, max_z_real
    DOUBLE PRECISION:: xtemp, ytemp, ztemp, x_ell, y_ell, z_ell
    DOUBLE PRECISION:: min_nu, max_nu, min_nu2, max_nu2
    DOUBLE PRECISION:: det, sq_g, Theta_a
    ! The value of nu equal for all the particles, used during the APM iteration
    DOUBLE PRECISION:: nu_all
    DOUBLE PRECISION:: err_N_mean_min, err_N_mean_min_old, err_N_mean, &
                       err_mean_old, err_n_min, err_N_max, dN, dNstar, &
                       nstar_p_err, nstar_real_err, r_tmp, dN_max, dN_av
    DOUBLE PRECISION:: art_pr_max
    DOUBLE PRECISION:: nu_tot, nu_ratio, nu_tmp2
    DOUBLE PRECISION:: variance_nu, stddev_nu, mean_nu
    DOUBLE PRECISION:: dist
    DOUBLE PRECISION:: rand_num, rand_num2

    INTEGER, DIMENSION(:), ALLOCATABLE:: neighbors_lists
    INTEGER, DIMENSION(:), ALLOCATABLE:: n_neighbors

    DOUBLE PRECISION, DIMENSION(3):: pos_corr_tmp
    DOUBLE PRECISION, DIMENSION(3):: pos_maxerr
    DOUBLE PRECISION:: g4(0:3,0:3)
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos_tmp
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: ghost_pos
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: all_pos
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: all_pos_tmp
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: all_pos_best
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: all_pos_tmp2
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: correction_pos
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: vel

    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: h_guess
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: h_tmp

    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar_p
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar_real
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: art_pr
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: freeze

    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nu_tmp
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nu_one

    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar_int
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: particle_density_final

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: nearest_neighbors

    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: lapse, &
                                               shift_x, shift_y, shift_z, &
                                               g_xx, g_xy, g_xz, &
                                               g_yy, g_yz, g_zz, &
                                               baryon_density, &
                                               energy_density, &
                                               specific_energy, &
                                               pressure, &
                                               v_euler_x, v_euler_y, v_euler_z

    LOGICAL:: exist

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    IF( debug ) PRINT *, "0"

    npart_real= SIZE( pos_input(1,:) )

    IF( mass == THIS% mass2 )THEN
      PRINT *, "npart_real=", npart_real
    ENDIF

    !---------------------------------------!
    !-- Allocate, assign and test h_guess --!
    !---------------------------------------!

    IF(.NOT.ALLOCATED( h_guess ))THEN
      ALLOCATE( h_guess( max_npart ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array h_guess in SUBROUTINE ", &
                  "perform_apm. The error message is",&
                  err_msg
         STOP
      ENDIF
    ENDIF

    h_guess= 0.0D0
    DO a= 1, npart_real, 1
      h_guess(a)= 3.0D0*(pvol(a)**third)
      IF( ISNAN( h_guess(a) ) )THEN
        PRINT *, " ** ERROR! h_guess(", a, &
                 ") is a NaN in SUBROUTINE perform_apm!"
        PRINT *, "pvol(", a, ")=", pvol(a)
        PRINT *, "    Stopping..."
        PRINT *
        STOP
      ENDIF
      IF( h_guess( a ) <= 0.0D0 )THEN
        PRINT *, "** ERROR! h_guess(", a, ") is zero or negative!"
        PRINT *, "   pvol(", a, ")=", pvol(a)
        PRINT *, "   Stopping..."
        PRINT *
        STOP
      ENDIF
    ENDDO

    IF( debug ) PRINT *, "0.5"

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

    !smaller_radius= ABS( MINVAL( pos_input( 1, : ), DIM= 1 ) - center )
    !larger_radius = ABS( center - MAXVAL( pos_input( 1, : ), DIM= 1 ) )
    !radius_y= ABS( MAXVAL( pos_input( 2, : ), DIM= 1 ) )
    !radius_z= ABS( MAXVAL( pos_input( 3, : ), DIM= 1 ) )

    IF( pos_input( 1, 10 ) < 0 )THEN

      smaller_radius= MIN( binary% get_radius1_x_comp(), &
                           binary% get_radius1_x_opp() )
      larger_radius = MAX( binary% get_radius1_x_comp(), &
                           binary% get_radius1_x_opp() )
      radius_y= binary% get_radius1_y()
      radius_z= binary% get_radius1_z()

    ELSE

      smaller_radius= MIN( binary% get_radius2_x_comp(), &
                           binary% get_radius2_x_opp() )
      larger_radius = MAX( binary% get_radius2_x_comp(), &
                           binary% get_radius2_x_opp() )
      radius_y= binary% get_radius2_y()
      radius_z= binary% get_radius2_z()

    ENDIF

    h_max= 0.0D0
    h_av = 0.0D0
    itr  = 0
    max_z_real= ABS( MAXVAL( pos_input( 3, : ), DIM= 1 ) )
    DO a= 1, npart_real, 1

      IF( SQRT( ( pos_input( 1, a ) - center )**2.0D0 &
                + pos_input( 2, a )**2.0D0 &
                + pos_input( 3, a )**2.0D0 ) > 0.99D0*max_z_real )THEN

        itr= itr + 1
        IF( h_guess(a) > h_max )THEN
          h_max= h_guess(a)
        ENDIF
        h_av= h_av + h_guess(a)

      ENDIF

    ENDDO
    h_av= h_av/itr
    IF( debug ) PRINT *, "h_av=", h_av
    IF( debug ) PRINT *

    IF( debug ) PRINT *, "2"

    !-------------------------------!
    !--  Placing ghost particles  --!
    !-------------------------------!

    IF(.NOT.ALLOCATED( ghost_pos ))THEN
      ALLOCATE( ghost_pos( 3, max_npart ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array ghost_pos in SUBROUTINE ", &
                  "perform_apm. The error message is",&
                  err_msg
         STOP
      ENDIF
    ENDIF

    ghost_pos= 0.0D0

    PRINT *, " * Placing ghost particles on a lattice between ellipsodial ", &
             "surfaces..."
    PRINT *

    max_r_real= 0.0D0
    DO itr= 1, npart_real, 1

      r_real= SQRT( ( pos_input( 1, itr ) - center )**2.0D0 &
                  + pos_input( 2, itr )**2.0D0 + pos_input( 3, itr )**2.0D0 )
      IF( r_real > max_r_real ) max_r_real= r_real

    ENDDO

    nx= nx_gh
    ny= ny_gh
    nz= nz_gh
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
    !delta= 1.0D0

    rad_x= larger_radius + h_av/1.0D0
    rad_y= radius_y + h_av/1.0D0
    rad_z= radius_z + h_av/1.0D0

    IF( debug ) PRINT *, "larger_radius= ", larger_radius
    IF( debug ) PRINT *, "radius_y= ", radius_y
    IF( debug ) PRINT *, "radius_z= ", radius_z
    IF( debug ) PRINT *, "rad_x= ", rad_x
    IF( debug ) PRINT *, "rad_y= ", rad_y
    IF( debug ) PRINT *, "rad_z= ", rad_z
    IF( debug ) PRINT *

   ! DO

      itr= 0
      DO k= 1, nz, 1

        ztemp= zmin + dz/2.0D0 + DBLE( k - 1 )*dz

        DO j= 1, ny, 1

          ytemp= ymin + dy/2.0D0 + DBLE( j - 1 )*dy

          DO i= 1, nx, 1

            xtemp= xmin + dx/2.0D0 + DBLE( i - 1 )*dx

            x_ell= center + rad_x*COS(ATAN( ytemp/( xtemp - center ) )) &
                   *SIN(ACOS(ztemp/SQRT( ( xtemp - center )**2.0D0 &
                                         + ytemp**2.0D0 + ztemp**2.0D0 )))

            y_ell= rad_y*SIN(ATAN( ytemp/( xtemp - center ) )) &
                   *SIN(ACOS(ztemp/SQRT( ( xtemp - center )**2.0D0 &
                                         + ytemp**2.0D0 + ztemp**2.0D0 )))

            z_ell= rad_z*( ztemp/SQRT( ( xtemp - center )**2.0D0 &
                                     + ytemp**2.0D0 + ztemp**2.0D0 ) )

            IF( SQRT( ( xtemp - center )**2.0D0 + ytemp**2.0D0 &
                      + ztemp**2.0D0 ) <= &
                      1.1D0*SQRT( ( x_ell - center )**2.0D0 &
                                  + y_ell**2.0D0 + z_ell**2.0D0 ) &
                .AND. &
                SQRT( ( xtemp - center )**2.0D0 + ytemp**2.0D0 &
                      + ztemp**2.0D0 ) >= &
                SQRT( ( x_ell - center )**2.0D0 + y_ell**2.0D0 &
                      + z_ell**2.0D0 ) &
                .AND. &
                binary% import_mass_density( xtemp, ytemp, ztemp ) <= 0.0D0 &
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
      !ghost_pos = ghost_pos( :, 1:npart_ghost )

      ! Test that there are no real particles too close to the ghost particles
    !  max_r_real= 0.0D0
    !  DO itr= 1, npart_real, 1
    !
    !    r_real= SQRT( ( pos_input( 1, itr ) - center )**2.0D0 &
    !                + pos_input( 2, itr )**2.0D0 + pos_input( 3, itr )**2.0D0 )
    !    IF( r_real > max_r_real ) max_r_real= r_real
    !
    !  ENDDO
    !
    !  min_r_ghost= HUGE(0.0D0)
    !  DO itr= 1, npart_ghost, 1
    !
    !    r_ghost= SQRT( ( ghost_pos( 1, itr ) - center )**2.0D0 &
    !                 + ghost_pos( 2, itr )**2.0D0 + ghost_pos( 3, itr )**2.0D0 )
    !    IF( r_ghost < min_r_ghost ) min_r_ghost= r_ghost
    !
    !  ENDDO
    !
    !  PRINT *, max_r_real, min_r_ghost
    !
    !  IF( min_r_ghost - max_r_real < 0 .OR. &
    !      ABS( min_r_ghost - max_r_real ) < larger_radius*0.005D0 )THEN
    !
    !    rad_x= rad_x*( 1.0D0 + 0.005D0 )
    !    rad_y= rad_y*( 1.0D0 + 0.005D0 )
    !    rad_z= rad_z*( 1.0D0 + 0.005D0 )
    !
    !  ELSE
    !
    !    EXIT
    !
    !  ENDIF
    !
    !ENDDO
    ghost_pos = ghost_pos( :, 1:npart_ghost )

    PRINT *, " * ", npart_ghost, " ghost particles placed around ", &
             npart_real, "real particles."
    PRINT *

    PRINT *, " * Printing ghost particles to file..."

    IF( PRESENT(namefile_pos_id) )THEN
      finalnamefile= namefile_pos_id
    ELSE
      finalnamefile= "apm_pos_id.dat"
    ENDIF

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

    PRINT *, " * Positions of ghost and real particles printed to ", &
             finalnamefile, " ."

    !STOP

    npart_all= npart_real + npart_ghost

    IF(.NOT.ALLOCATED( all_pos ))THEN
      ALLOCATE( all_pos( 3, npart_all ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array ghost_pos in SUBROUTINE ", &
                  "perform_apm. The error message is",&
                  err_msg
         STOP
      ENDIF
    ENDIF

    all_pos( :, 1:npart_real )          = pos_input
    all_pos( :, npart_real+1:npart_all )= ghost_pos

    h_guess= h_guess(1:npart_all)
    h_guess(npart_real+1:npart_all)= ( dx*dy*dz )**third

    PRINT *, "SIZE(npart_all)= ", SIZE(h_guess)
    PRINT *

    !----------------------------!
    !-- Allocate needed memory --!
    !----------------------------!

    ! setup unit system
    !CALL set_units('NSM')
    !CALL read_options       ! TODO: set units and read options only once in the constructor

    npart= npart_all
    CALL allocate_SPH_memory

    CALL allocate_RCB_tree_memory_3D(npart)
    iorig(1:npart)= (/ (a,a=1,npart) /)

    IF( debug ) PRINT *, "10"

    ! tabulate kernel, get ndes
    !CALL ktable(ikernel,ndes)

    ! flag that particles are 'alive'
    !IF( .NOT.ALLOCATED( alive ) ) ALLOCATE( alive( npart ) )
    !alive= 1

    CALL allocate_gradient( npart )
    CALL allocate_metric_on_particles( npart )

    !------------------------------------------!
    !-- Apply the artificial pressure method --!
    !------------------------------------------!

    PRINT *, "** Setting up ID for APM iteration..."
    PRINT *

    PRINT *, " * Assign h..."
    PRINT *
    CALL assign_h( nn_des, &
                   npart_all, &
                   all_pos, h_guess, &
                   h )

    DO a= 1, npart_all, 1

      IF( ISNAN( h( a ) ) )THEN
        PRINT *, "** ERROR! h(", a, ") is a NaN!"
        PRINT *, " * h_guess(", a, ")= ", h_guess(a)
        PRINT *, " * all_pos(:,", a, ")= ", all_pos(:,a)
        PRINT *, " Stopping..."
        PRINT *
        STOP
      ENDIF
      IF( h( a ) <= 0.0D0 )THEN
        PRINT *, "** ERROR! h(", a, ") is zero or negative!"
        PRINT *, " * h_guess(", a, ")= ", h_guess(a)
        PRINT *, " * all_pos(:,", a, ")= ", all_pos(:,a)
        PRINT *, " * h(", a, ")= ", h(a)
        PRINT *, " Stopping..."
        PRINT *
        STOP
      ENDIF

    ENDDO

    PRINT *, " * Measure SPH particle number density..."
    PRINT *

    IF(.NOT.ALLOCATED( nstar_real ))THEN
      ALLOCATE( nstar_real( npart_all ), STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array nstar_real in SUBROUTINE ", &
                  "perform_apm. The error message is",&
                  err_msg
         STOP
      ENDIF
    ENDIF

    nu= 1.0D0
    CALL density_loop( npart_all, all_pos, &    ! input
                       nu, h, nstar_real )      ! output

    DO a= 1, npart_all, 1

      IF( ISNAN( nstar_real( a ) ) )THEN

        PRINT *, "** WARNING! nstar_real(", a, ") is a NaN!", &
                 "   Changing the values of nstar_real and h to one taken", &
                 "   from a neighboring particle."
        IF( debug ) PRINT *, " * h(", a, ")=", h(a)
        IF( debug ) PRINT *, " * nu(", a, ")=", nu(a)
        IF( debug ) PRINT *, " * all_pos(", a, ")=", all_pos(:,a)
        IF( debug ) PRINT *, " * r(", a, ")=", &
                              SQRT( ( all_pos(1,a) - center )**2.0D0 &
                              + all_pos(2,a)**2.0D0 + all_pos(3,a)**2.0D0 )
        PRINT *, " * Check if the smoothing length is 0 for some particles,", &
                 "   and if so, make its initial guess, h_guess, a bit larger."
        PRINT *

      ENDIF

    ENDDO

    IF( debug ) PRINT *, "4"

    IF(.NOT.ALLOCATED( nstar_p ))THEN
      ALLOCATE( nstar_p( npart_all ), STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array nstar_p in SUBROUTINE ", &
                  "perform_apm. The error message is",&
                  err_msg
         STOP
      ENDIF
    ENDIF

    IF( debug ) PRINT *, "5"

    ALLOCATE( lapse          (npart_real) )
    ALLOCATE( shift_x        (npart_real) )
    ALLOCATE( shift_y        (npart_real) )
    ALLOCATE( shift_z        (npart_real) )
    ALLOCATE( g_xx           (npart_real) )
    ALLOCATE( g_xy           (npart_real) )
    ALLOCATE( g_xz           (npart_real) )
    ALLOCATE( g_yy           (npart_real) )
    ALLOCATE( g_yz           (npart_real) )
    ALLOCATE( g_zz           (npart_real) )
    ALLOCATE( baryon_density (npart_real) )
    ALLOCATE( energy_density (npart_real) )
    ALLOCATE( specific_energy(npart_real) )
    ALLOCATE( pressure       (npart_real) )
    ALLOCATE( v_euler_x      (npart_real) )
    ALLOCATE( v_euler_y      (npart_real) )
    ALLOCATE( v_euler_z      (npart_real) )

    IF(.NOT.ALLOCATED( vel ))THEN
      ALLOCATE( vel( 0:3, npart_real ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array vel in SUBROUTINE ", &
                  "perform_apm. The error message is",&
                  err_msg
         STOP
      ENDIF
    ENDIF

    IF( debug ) PRINT *, "6"

    CALL binary% import_id( npart_real, all_pos(1,1:npart_real), &
                            all_pos(2,1:npart_real), &
                            all_pos(3,1:npart_real), &
                            lapse, shift_x, shift_y, shift_z, &
                            g_xx, g_xy, g_xz, &
                            g_yy, g_yz, g_zz, &
                            baryon_density, &
                            energy_density, &
                            specific_energy, &
                            pressure, &
                            v_euler_x, v_euler_y, v_euler_z )

    IF( debug ) PRINT *, "7"

    max_nu= 0.0D0
    min_nu= 1.0D60

    IF(.NOT.ALLOCATED( art_pr ))THEN
      ALLOCATE( art_pr( npart_all ), STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array art_pr in SUBROUTINE ", &
                  "perform_apm. The error message is",&
                  err_msg
         STOP
      ENDIF
    ENDIF

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
      Theta(a)= Theta_a

      nstar_p(a)= sq_g*Theta_a*baryon_density(a)*((Msun_geo*km2m)**3)/(amu*g2kg)

      IF( ISNAN( nstar_p( a ) ) )THEN
        PRINT *, "** ERROR! nstar_p(", a, ") is a NaN!", &
                 " Stopping.."
        PRINT *
        STOP
      ENDIF
      IF( nstar_p( a ) == 0 )THEN
        PRINT *, "** ERROR! When setting up ID: ", &
                 "nstar_p(", a, ")= 0 on a real particle!"
        PRINT *, " * Particle position: x=", all_pos(1,a), &
                 ", y=", all_pos(2,a), ", z=", all_pos(3,a)
        PRINT *, "   sq_g=", sq_g
        PRINT *, "   Theta_a=", Theta_a
        PRINT *, "   baryon_density(a)=", baryon_density(a)
        PRINT *, " * Stopping.."
        PRINT *
        STOP
      ENDIF

    ENDDO

    IF( debug ) PRINT *, "8"

    !----------------------------------------------------!
    !-- enforce centre of mass after having changed nu --!
    !----------------------------------------------------!

    CALL COM( npart_real, all_pos(:,1:npart_real), nu(1:npart_real), & ! input
              com_x, com_y, com_z, com_d ) ! output

    IF( debug ) PRINT *, "9"

    PRINT *, "** Before center of mass correction:"
    PRINT *, " * x coordinate of the center of mass of the star, ", &
             "from LORENE: com_star= ", com_star, "Msun_geo"
    PRINT *, " * x coordinate of the center of mass of the particle ", &
             "distribution: com_x= ", com_x, "Msun_geo"
    PRINT *, " * y coordinate of the center of mass of the particle ", &
             "distribution: com_y= ", com_y, "Msun_geo"
    PRINT *, " * z coordinate of the center of mass of the particle ", &
             "distribution: com_z= ", com_z, "Msun_geo"
    PRINT *, " * Distance of the center of mass of the particle ", &
             "distribution from the  origin: com_d= ", com_d
    PRINT *, " * |com_x-com_star/com_star|=", &
             ABS( com_x-com_star )/ABS( com_star )
    PRINT *

    DO a= 1, npart_real, 1

      pos_corr_tmp(1)= all_pos(1,a) - ( com_x - com_star )
      pos_corr_tmp(2)= all_pos(2,a) - com_y
      pos_corr_tmp(3)= all_pos(3,a) - com_z

      IF( binary% import_mass_density( &
                  pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) > 0.0D0 &
          .AND. &
          binary% is_hydro_negative( &
                  pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) == 0 &
      )THEN

        all_pos(:,a)= pos_corr_tmp

      ENDIF

    ENDDO

    IF( debug ) PRINT *, "10"

    CALL COM( npart_real, all_pos(:,1:npart_real), nu(1:npart_real), & ! input
              com_x, com_y, com_z, com_d ) ! output

    IF( debug ) PRINT *, "9"

    PRINT *, "** After center of mass correction:"
    PRINT *, " * x coordinate of the center of mass of the star, ", &
             "from LORENE: com_star= ", com_star, "Msun_geo"
    PRINT *, " * x coordinate of the center of mass of the particle ", &
             "distribution: com_x= ", com_x, "Msun_geo"
    PRINT *, " * y coordinate of the center of mass of the particle ", &
             "distribution: com_y= ", com_y, "Msun_geo"
    PRINT *, " * z coordinate of the center of mass of the particle ", &
             "distribution: com_z= ", com_z, "Msun_geo"
    PRINT *, " * Distance of the center of mass of the particle ", &
             "distribution from the  origin: com_d= ", com_d
    PRINT *, " * |com_x-com_star/com_star|=", &
             ABS( com_x-com_star )/ABS( com_star )
    PRINT *

    !-----------------------------------------------------------------------!
    !-- Mirror the positions after having repositioned the center of mass --!
    !-----------------------------------------------------------------------!

    pos_tmp= all_pos(:,1:npart_real)
    nu_tmp= nu
    itr= 0
    DO a= 1, npart_real, 1
      IF( pos_tmp( 3, a ) > 0.0D0 )THEN
        itr= itr + 1
        all_pos( 1, itr )= pos_tmp( 1, a )
        all_pos( 2, itr )= pos_tmp( 2, a )
        all_pos( 3, itr )= pos_tmp( 3, a )
        nu( itr )        = nu_tmp( a )
      ENDIF
    ENDDO
    npart_real_half= itr

    DO a= 1, npart_real_half, 1
      all_pos( 1, npart_real_half + a )=   all_pos( 1, a )
      all_pos( 2, npart_real_half + a )=   all_pos( 2, a )
      all_pos( 3, npart_real_half + a )= - all_pos( 3, a )
      nu( npart_real_half + a )        =   nu( a )
    ENDDO

    PRINT *, "** After mirroring particles:"
    PRINT *, " * x coordinate of the center of mass of the star, ", &
             "from LORENE: com_star= ", com_star, "Msun_geo"
    PRINT *, " * x coordinate of the center of mass of the particle ", &
             "distribution: com_x= ", com_x, "Msun_geo"
    PRINT *, " * y coordinate of the center of mass of the particle ", &
             "distribution: com_y= ", com_y, "Msun_geo"
    PRINT *, " * z coordinate of the center of mass of the particle ", &
             "distribution: com_z= ", com_z, "Msun_geo"
    PRINT *, " * Distance of the center of mass of the particle ", &
             "distribution from the  origin: com_d= ", com_d
    PRINT *, " * |com_x-com_star/com_star|=", &
             ABS( com_x-com_star )/ABS( com_star )
    PRINT *

    PRINT *, " * ID set up for the APM iteration."
    PRINT *

    !-------------------------------------------------!
    !--               APM iteration                 --!
    !-- Assume equal mass particles, and move them  --!
    !-- so that the SPH kernel estimate of the mass --!
    !-- density matches the star mass density as    --!
    !-- well as reasonably possible.                --!
    !-------------------------------------------------!

    PRINT *, " * Performing APM iteration..."
    PRINT *

    ALLOCATE( freeze( npart_all ) )
    ALLOCATE( correction_pos( 3, npart_all ) )
    ALLOCATE( all_pos_tmp( 3, npart_all ) )
    ALLOCATE( all_pos_tmp2( 3, npart_all ) )

    ! Set the particles to be equal-mass
    nu_all= (mass/DBLE(npart_real))*umass/amu
    nu= nu_all
    DO a= 1, npart_all
      IF( a < npart_real )THEN
        freeze(a)= 0
      ELSE
        freeze(a)= 1
      ENDIF
    ENDDO

    CALL COM( npart_real, all_pos(:,1:npart_real), nu(1:npart_real), & ! input
              com_x, com_y, com_z, com_d ) ! output

    IF( debug ) PRINT *, "9"

    PRINT *, "** Before center of mass correction:"
    PRINT *, " * x coordinate of the center of mass of the star, ", &
             "from LORENE: com_star= ", com_star, "Msun_geo"
    PRINT *, " * x coordinate of the center of mass of the particle ", &
             "distribution: com_x= ", com_x, "Msun_geo"
    PRINT *, " * y coordinate of the center of mass of the particle ", &
             "distribution: com_y= ", com_y, "Msun_geo"
    PRINT *, " * z coordinate of the center of mass of the particle ", &
             "distribution: com_z= ", com_z, "Msun_geo"
    PRINT *, " * Distance of the center of mass of the particle ", &
             "distribution from the  origin: com_d= ", com_d
    PRINT *, " * |com_x-com_star/com_star|=", &
             ABS( com_x-com_star )/ABS( com_star )
    PRINT *

    DO a= 1, npart_real, 1

      pos_corr_tmp(1)= all_pos(1,a) - ( com_x - com_star )
      pos_corr_tmp(2)= all_pos(2,a) - com_y
      pos_corr_tmp(3)= all_pos(3,a) - com_z

      IF( binary% import_mass_density( &
                  pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) > 0.0D0 &
          .AND. &
          binary% is_hydro_negative( &
                  pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) == 0 &
      )THEN

        all_pos(:,a)= pos_corr_tmp

      ENDIF

    ENDDO

    IF( debug ) PRINT *, "10"

    CALL COM( npart_real, all_pos(:,1:npart_real), nu(1:npart_real), & ! input
              com_x, com_y, com_z, com_d ) ! output

    IF( debug ) PRINT *, "9"

    PRINT *, "** After center of mass correction:"
    PRINT *, " * x coordinate of the center of mass of the star, ", &
             "from LORENE: com_star= ", com_star, "Msun_geo"
    PRINT *, " * x coordinate of the center of mass of the particle ", &
             "distribution: com_x= ", com_x, "Msun_geo"
    PRINT *, " * y coordinate of the center of mass of the particle ", &
             "distribution: com_y= ", com_y, "Msun_geo"
    PRINT *, " * z coordinate of the center of mass of the particle ", &
             "distribution: com_z= ", com_z, "Msun_geo"
    PRINT *, " * Distance of the center of mass of the particle ", &
             "distribution from the  origin: com_d= ", com_d
    PRINT *, " * |com_x-com_star/com_star|=", &
             ABS( com_x-com_star )/ABS( com_star )
    PRINT *

    all_pos_tmp2= -1.0D0
    PRINT *, " * The APM iteration starts here."
    PRINT *

    n_inc= 0
    err_N_mean_min= HUGE(1.0D0)
    apm_iteration: DO itr= 1, apm_max_it, 1

      IF( MOD( itr, 15 ) == 0 )THEN

        IF( debug ) PRINT *, "printing positions to file..."

        IF( PRESENT(namefile_pos) )THEN
          finalnamefile= namefile_pos
        ELSE
          finalnamefile= "apm_pos.dat"
        ENDIF

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
            all_pos( 1, a ), &
            all_pos( 2, a ), &
            all_pos( 3, a )
        ENDDO

        DO a= npart_real + 1, npart_all, 1
          WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
            2, a, &
            all_pos( 1, a ), &
            all_pos( 2, a ), &
            all_pos( 3, a )
        ENDDO

        CLOSE( UNIT= 2 )

      ENDIF

      IF( debug ) PRINT *, "mirroring particles..."

      ! Find particles above the xy plane
      all_pos_tmp= all_pos(:,1:npart_real)
      itr2= 0
      DO a= 1, npart_real, 1

        IF( all_pos_tmp( 3, a ) > 0.0D0 &
            .AND. &
            itr2 < npart_real/2 &
        )THEN
          itr2= itr2 + 1
          all_pos( 1, itr2 )= all_pos_tmp( 1, a )
          all_pos( 2, itr2 )= all_pos_tmp( 2, a )
          all_pos( 3, itr2 )= all_pos_tmp( 3, a )
        ENDIF

      ENDDO
      npart_real_half= itr2

      ! If some of the particles crossed the xy plane top-down in the
      ! last step, replace them with their previous position
      ! above the xy plane
   !   IF( npart_real_half < npart_real/2 )THEN
   !
   !     npart_missing= npart_real/2 - npart_real_half
   !
   !     DO a= npart_real_half + 1, npart_real/2, 1
   !
   !       all_pos( :, a )= all_pos_tmp2( :, a )
   !
   !     ENDDO
   !
   !   ENDIF

      ! Mirror the particles above the xy plane, to below the xy plane
      DO a= 1, npart_real/2, 1
        all_pos( 1, npart_real/2 + a )=   all_pos( 1, a )
        all_pos( 2, npart_real/2 + a )=   all_pos( 2, a )
        all_pos( 3, npart_real/2 + a )= - all_pos( 3, a )
      ENDDO

      IF( debug )THEN

        CALL COM( npart_real, all_pos(:,1:npart_real), nu(1:npart_real), & ! input
                  com_x, com_y, com_z, com_d ) ! output

        PRINT *, "** After center of mass correction:"
        PRINT *, " * x coordinate of the center of mass of the star, ", &
                 "from LORENE: com_star= ", com_star, "Msun_geo"
        PRINT *, " * x coordinate of the center of mass of the particle ", &
                 "distribution: com_x= ", com_x, "Msun_geo"
        PRINT *, " * y coordinate of the center of mass of the particle ", &
                 "distribution: com_y= ", com_y, "Msun_geo"
        PRINT *, " * z coordinate of the center of mass of the particle ", &
                 "distribution: com_z= ", com_z, "Msun_geo"
        PRINT *, " * Distance of the center of mass of the particle ", &
                 "distribution from the  origin: com_d= ", com_d
        PRINT *, " * |com_x-com_star/com_star|=", &
                 ABS( com_x-com_star )/ABS( com_star )
        PRINT *

        finalnamefile= "dbg-pos.dat"

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

        DO a= 1, npart_real/2, 1
          WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
            1, a, &
            all_pos( 1, a ), &
            all_pos( 2, a ), &
            all_pos( 3, a )
        ENDDO

        DO a= npart_real/2+1, npart_real, 1
          WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
            2, a, &
            all_pos( 1, a ), &
            all_pos( 2, a ), &
            all_pos( 3, a )
        ENDDO

        DO a= npart_real+1, npart_all, 1
          WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
            3, a, &
            all_pos( 1, a ), &
            all_pos( 2, a ), &
            all_pos( 3, a )
        ENDDO

        CLOSE( UNIT= 2 )

      ENDIF

      IF( debug ) PRINT *, "assign h..."

      h_guess(1:npart_real)= h
      !h_guess(npart_real+1:npart_all)= dx*dy*dz
      CALL assign_h( nn_des, &
                     npart_all, &
                     all_pos, h_guess, &
                     h )

      find_problem_in_h: DO a= 1, npart_all, 1

        IF( ISNAN( h( a ) ) )THEN
          PRINT *, "** ERROR! h(", a, ") is a NaN!"
          PRINT *, " * h_guess(", a, ")= ", h_guess(a)
          PRINT *, " * all_pos(:,", a, ")= ", all_pos(:,a)
          PRINT *, " Stopping..."
          PRINT *
          STOP
        ENDIF
        IF( h( a ) <= 0.0D0 )THEN
          PRINT *, "** ERROR! h(", a, ") is zero or negative!"
          PRINT *, " * h_guess(", a, ")= ", h_guess(a)
          PRINT *, " * all_pos(:,", a, ")= ", all_pos(:,a)
          PRINT *, " * h(", a, ")= ", h(a)
          PRINT *, " Stopping..."
          PRINT *
          STOP
        ENDIF

      ENDDO find_problem_in_h

      IF( debug ) PRINT *, "density_loop..."

      CALL density_loop( npart_all, all_pos, &    ! input
                         nu, h, nstar_real )      ! output

      PRINT *, "npart_real= ", npart_real
      PRINT *, "npart_all= ", npart_all
      PRINT *

      find_nan_in_nstar_real: DO a= 1, npart_all, 1

        IF( ISNAN( nstar_real( a ) ) )THEN

          PRINT *, "** WARNING! nstar_real(", a, ") is a NaN!", &
                   "   Changing the values of nstar_real and h to one taken", &
                   "   from a neighboring particle."
          IF( debug ) PRINT *, " * h(", a, ")=", h(a)
          IF( debug ) PRINT *, " * nu(", a, ")=", nu(a)
          IF( debug ) PRINT *, " * all_pos(", a, ")=", all_pos(:,a)
          IF( debug ) PRINT *, " * r(", a, ")=", &
                                SQRT( ( all_pos(1,a) - center )**2.0D0 &
                                + all_pos(2,a)**2.0D0 + all_pos(3,a)**2.0D0 )
          PRINT *
          STOP

          IF( a == 1 )THEN
            DO a2= 2, npart_all, 1
              IF( .NOT.ISNAN( nstar_real( a2 ) ) )THEN
                nstar_real( a )= nstar_real( a2 )
                h( a )= h( a2 )
                EXIT
              ENDIF
            ENDDO
          ELSEIF( npart_real == a )THEN
            nstar_real( a )= nstar_real( a - 1 )
            h( a )= h( a - 1 )
          ELSEIF( npart_real + 1 == a )THEN
            nstar_real( a )= nstar_real( a + 1 )
            h( a )= h( a + 1 )
          !ELSEIF( npart_all == a )THEN
          !  nstar_real( a )= nstar_real( a - 1 )
          ELSE
            nstar_real( a )= nstar_real( a - 1 )
            h( a )= h( a - 1 )
          ENDIF
          ! TODO: here you need a recursive SUBROUTINE
          !CALL RANDOM_NUMBER( rand_num )
          !CALL RANDOM_NUMBER( rand_num2 )
          !IF( rand_num2 < half )  rel_sign= - 1
          !IF( rand_num2 >= half ) rel_sign=   1
          !pos_corr_tmp= all_pos(:,a)*( 1.0D0 + rel_sign*rand_num/100.0D0 )
          !IF( binary% import_mass_density( &
          !        pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) > 0.0D0 &
          !    .AND. &
          !    binary% is_hydro_negative( &
          !        pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) == 0 &
          !)THEN
          !
          !  all_pos(:,a)= pos_corr_tmp
          !
          !ENDIF
          !CALL density_loop( npart_all, all_pos, &    ! input
          !                   nu, h, nstar_real )      ! output

          !PRINT *, "** ERROR! nstar_real(", a, ") is a NaN!"
          !PRINT *, " * h(", a, ")=", h(a)
          !PRINT *, " * nu(", a, ")=", nu(a)
          !PRINT *, " * all_pos(", a, ")=", all_pos(:,a)
          !PRINT *, " * r(", a, ")=", SQRT( ( all_pos(1,a) - center )**2.0D0 &
          !                      + all_pos(2,a)**2.0D0 + all_pos(3,a)**2.0D0 )
          !PRINT *, " Stopping.."
          !PRINT *
          !STOP

        ENDIF

      ENDDO find_nan_in_nstar_real

      CALL binary% import_id( npart_real, all_pos(1,1:npart_real), &
                              all_pos(2,1:npart_real), &
                              all_pos(3,1:npart_real), &
                              lapse, shift_x, shift_y, shift_z, &
                              g_xx, g_xy, g_xz, &
                              g_yy, g_yz, g_zz, &
                              baryon_density, &
                              energy_density, &
                              specific_energy, &
                              pressure, &
                              v_euler_x, v_euler_y, v_euler_z )

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
        Theta(a)= Theta_a

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
          PRINT *, " * Particle position: x=", all_pos(1,a), &
                   ", y=", all_pos(2,a), ", z=", all_pos(3,a)
          PRINT *, "   sq_g=", sq_g
          PRINT *, "   Theta_a=", Theta_a
          PRINT *, "   baryon_density(a)=", baryon_density(a)
          PRINT *, " * Stopping.."
          PRINT *
          STOP
        ENDIF

      ENDDO

      art_pr_max= 0.0D0
      err_N_max=  0.0D0
      err_N_min=  1.D30
      err_N_mean= 0.0D0
      DO a= 1, npart_real, 1

        dNstar= ( nstar_real(a) - nstar_p(a) )/nstar_p(a)
        art_pr(a) = MAX( 1.0D0 + dNstar, 0.1D0 )
        art_pr_max= MAX( art_pr_max, art_pr(a) )
        IF( ABS(dNstar) > err_N_max )THEN
          err_N_max     = ABS(dNstar)
          pos_maxerr    = all_pos(:,a)
          nstar_real_err= nstar_real(a)
          nstar_p_err   = nstar_p(a)
        ENDIF
        !err_N_max = MAX( err_N_max, ABS(dNstar) )
        err_N_min = MIN( err_N_min, ABS(dNstar) )
        err_N_mean= err_N_mean + ABS(dNstar)

        IF( ISNAN(dNstar) )THEN
          PRINT *, "dNstar is a NaN at particle ", a
          PRINT *, "nstar_real is ", nstar_real(a)
          PRINT *, "nstar_p is ", nstar_p(a)
          STOP
        ENDIF

      ENDDO

      IF( ISNAN( art_pr_max ) )THEN
        PRINT *, "** ERROR! art_pr_max is a NaN!", &
                 " Stopping.."
        PRINT *
        STOP
      ENDIF
      IF( ISNAN( nu_all ) )THEN
        PRINT *, "** ERROR! nu_all is a NaN!", &
                 " Stopping.."
        PRINT *
        STOP
      ENDIF

      nstar_p( npart_real+1:npart_all )= 0.0D0
      art_pr ( npart_real+1:npart_all )= 3.0D0*art_pr_max

      PRINT *, "Before calling position_correction"

      !PRINT *, npart_all
      !PRINT *, SIZE(all_pos(1,:))
      !PRINT *, SIZE(h)
      !PRINT *, SIZE(art_pr)
      !PRINT *, SIZE(nstar_real)
      !PRINT *, SIZE(correction_pos(1,:))

      find_nan_in_art_pr: DO a= 1, npart_all, 1

        IF( ISNAN( art_pr( a ) ) )THEN
          PRINT *, "** ERROR! art_pr(", a, ") is a NaN!", &
                   " Stopping.."
          PRINT *
          STOP
        ENDIF
        IF( art_pr( a ) > HUGE(1.0D0) )THEN
          PRINT *, "** ERROR! art_pr(", a, ")= ", art_pr( a ), " is infinite!",&
                   " Stopping.."
          PRINT *
          STOP
        ENDIF

      ENDDO find_nan_in_art_pr

      find_nan_in_all_pos: DO a= 1, npart_all, 1

        DO itr2= 1, 3, 1
          IF( ISNAN( all_pos( itr2, a ) ) )THEN
            PRINT *, "** ERROR! all_pos(", itr2, ",", a, ") is a NaN!", &
                     " Stopping.."
            PRINT *
            STOP
          ENDIF
        ENDDO

      ENDDO find_nan_in_all_pos

      all_pos_tmp2= all_pos

      CALL position_correction( npart_all, &
                                all_pos, h, nu_all, art_pr, nstar_real, &
                                correction_pos )

      find_nan_in_correction_pos: DO a= 1, npart_all, 1

        DO itr2= 1, 3, 1
          IF( ISNAN( correction_pos( itr2, a ) ) )THEN
            PRINT *, "** ERROR! correction_pos(", itr2, ",", a, ") is a NaN!"
            PRINT *, " *        Particle position: x=", all_pos(1,a), &
                     ", y=", all_pos(2,a), ", z=", all_pos(3,a)
            r_tmp= SQRT( ( all_pos(1,a) - center )**2.0D0 + &
                           all_pos(2,a)**2.0D0 + all_pos(3,a)**2.0D0 )
            PRINT *, " *        Particle radius: r=", r_tmp, &
                     "=", r_tmp/larger_radius*100.0D0, &
                     "% of the larger radius of the star."
            PRINT *, " *        Particle colatitude: theta=", &
                     ACOS( all_pos(3,a)/r_tmp )/pi, " pi"
            PRINT *, " *        Particle longitude: phi=", &
                     ATAN( all_pos(2,a)/all_pos(1,a) )/pi, " pi"
            PRINT *, " * Stopping.."
            PRINT *
            STOP
          ENDIF
        ENDDO

      ENDDO find_nan_in_correction_pos

      IF( debug ) PRINT *, "After calling position_correction"

      itr2= 0
      DO a= 1, npart_real, 1
        pos_corr_tmp= all_pos(:,a) + correction_pos(:,a)
        IF( &!binary% import_mass_density( &
            !                    all_pos(1,a), all_pos(2,a), all_pos(3,a) ) > 0 &
            !.AND. &
            binary% import_mass_density( &
                  pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) > 0.0D0 &
            .AND. &
            nstar_p(a) > 0.0D0 &
            .AND. &
            binary% is_hydro_negative( &
                  pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) == 0 &
        )THEN
          itr2= itr2 + 1
          !IF( binary% is_hydro_negative( &
          !                        all_pos(1,a), all_pos(2,a), all_pos(3,a) ) > 0 &
          !    .OR. &
          !    binary% is_hydro_negative( &
          !                        pos_tmp(1), pos_tmp(2), pos_tmp(3) ) > 0 &
          !)THEN
            all_pos(:,a)= pos_corr_tmp
          !ELSE
          !  all_pos(:,a)= all_pos(:,a) + correction_pos(:,a)/2.0D0
          !ENDIF
        ENDIF
      ENDDO

      ! If some of the particles crossed the xy plane top-down in the
      ! last step, reflect them back above the xy plane
      DO a= 1, npart_real, 1

        IF( all_pos_tmp2( 3, a ) > 0 .AND. &
            all_pos( 3, a ) <= 0 )THEN

          all_pos( 3, a )= all_pos_tmp2( 3, a )

        ENDIF

      ENDDO

      IF( debug ) PRINT *, "After correcting positions"

      err_N_mean= err_N_mean/DBLE(npart_real)

      err_N_mean_min= MIN( err_N_mean, err_N_mean_min )

      IF( err_N_mean_min == err_N_mean )THEN
        all_pos_best= all_pos
      ENDIF

      !PRINT *, "itr2= ", itr2
      PRINT *, ' * Done with position update #: ', itr
      PRINT *
      PRINT *, '   err_N_max=  ', err_N_max
      PRINT *, "   at pos=", pos_maxerr
      PRINT *, "   with r/r_x_opp= ", SQRT( &
                                    ( ABS(pos_maxerr(1)) - ABS(center) )**2.0D0 &
                                          + pos_maxerr(2)**2.0D0 &
                                          + pos_maxerr(3)**2.0D0 ) &
                                    /smaller_radius
      PRINT *
      PRINT *, "   nstar_real_err= ", nstar_real_err
      PRINT *, "   nstar_p_err   = ", nstar_p_err
      PRINT *
      PRINT *, '   err_N_min=  ', err_N_min
      PRINT *, '   err_N_mean= ', err_N_mean
      PRINT *, '   err_N_mean_min= ', err_N_mean_min
      PRINT *

      PRINT *, "   max_nu=", MAXVAL( nstar_p(1:npart_real)/nstar_real(1:npart_real), DIM= 1 )
      PRINT *, "   min_nu=", MINVAL( nstar_p(1:npart_real)/nstar_real(1:npart_real), DIM= 1 )
      PRINT *, "   max_nu/min_nu=", MAXVAL( nstar_p(1:npart_real)/nstar_real(1:npart_real), DIM= 1 )/MINVAL( nstar_p(1:npart_real)/nstar_real(1:npart_real), DIM= 1 )
      PRINT *

      ! Exit condition
      IF( err_N_mean > err_mean_old )THEN
        n_inc= n_inc + 1
      ENDIF
      PRINT *, "n_inc= ", n_inc
      PRINT *
      !IF( ABS( err_N_mean - err_mean_old )/ABS( err_mean_old ) < iter_tol &
      !    .AND. &
      !    err_N_max < 10.0D0 &
      !)THEN
      !  n_inc= n_inc + 1
      !  PRINT *, "n_inc/max_inc= ", n_inc, "/", max_inc
      !  PRINT *, "ABS( err_N_mean - err_mean_old )/ABS(err_mean_old)= ", &
      !           ABS( err_N_mean - err_mean_old )/ABS(err_mean_old)
      !ENDIF
      !IF( ABS(err_N_mean_min - err_N_mean_min_old)/ABS(err_N_mean_min_old) &
      !      < 10.0D0*iter_tol &
      !    .AND. &
      !    ABS(err_N_mean - err_N_mean_min)/ABS(err_N_mean_min) < iter_tol &
      !)THEN
      !  n_inc= n_inc + 1
      !  PRINT *, "n_inc/max_inc= ", n_inc, "/", max_inc
      !  PRINT *, err_N_mean, "err_N_mean"
      !  PRINT *, err_N_mean_min, "err_N_mean_min"
      !  PRINT *, "ABS(err_N_mean - err_N_mean_min)/ABS(err_N_mean_min)= ", &
      !           ABS(err_N_mean - err_N_mean_min)/ABS(err_N_mean_min), " < ", &
      !           iter_tol
      !ELSE
      !  n_inc= 0
      !ENDIF
      IF( n_inc == max_inc ) EXIT
      err_mean_old      = err_N_mean
      err_N_mean_min_old= err_N_mean_min

    ENDDO apm_iteration

    PRINT *, "Iteration completed."
    PRINT *

    ! Now get rid of the ghost particles
    IF(.NOT.ALLOCATED( pos ))THEN
      ALLOCATE( pos( 3, npart_real ), STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array pos in SUBROUTINE ", &
                  "perform_apm. The error message is",&
                  err_msg
         STOP
      ENDIF
    ENDIF

    pos= all_pos( :, 1:npart_real )
    npart= npart_real
    PRINT *, npart

    h_guess= h

    !----------------------------!
    !-- enforce centre of mass --!
    !----------------------------!

    CALL COM( npart_real, pos, nu(1:npart_real), & ! input
              com_x, com_y, com_z, com_d ) ! output

    IF( debug ) PRINT *, "9"

    PRINT *, "** Before center of mass correction:"
    PRINT *, " * x coordinate of the center of mass of the star, ", &
             "from LORENE: com_star= ", com_star, "Msun_geo"
    PRINT *, " * x coordinate of the center of mass of the particle ", &
             "distribution: com_x= ", com_x, "Msun_geo"
    PRINT *, " * y coordinate of the center of mass of the particle ", &
             "distribution: com_y= ", com_y, "Msun_geo"
    PRINT *, " * z coordinate of the center of mass of the particle ", &
             "distribution: com_z= ", com_z, "Msun_geo"
    PRINT *, " * Distance of the center of mass of the particle ", &
             "distribution from the  origin: com_d= ", com_d
    PRINT *, " * |com_x-com_star/com_star|=", &
             ABS( com_x-com_star )/ABS( com_star )
    PRINT *

    DO a= 1, npart_real, 1

      pos_corr_tmp(1)= pos(1,a) - ( com_x - com_star )
      pos_corr_tmp(2)= pos(2,a) - com_y
      pos_corr_tmp(3)= pos(3,a) - com_z

      IF( binary% import_mass_density( &
                  pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) > 0.0D0 &
          .AND. &
          binary% is_hydro_negative( &
                  pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) == 0 &
      )THEN

        pos(:,a)= pos_corr_tmp

      ENDIF

    ENDDO

    IF( debug ) PRINT *, "10"

    CALL COM( npart_real, pos, nu(1:npart_real), & ! input
              com_x, com_y, com_z, com_d ) ! output

    IF( debug ) PRINT *, "9"

    PRINT *, "** After center of mass correction:"
    PRINT *, " * x coordinate of the center of mass of the star, ", &
             "from LORENE: com_star= ", com_star, "Msun_geo"
    PRINT *, " * x coordinate of the center of mass of the particle ", &
             "distribution: com_x= ", com_x, "Msun_geo"
    PRINT *, " * y coordinate of the center of mass of the particle ", &
             "distribution: com_y= ", com_y, "Msun_geo"
    PRINT *, " * z coordinate of the center of mass of the particle ", &
             "distribution: com_z= ", com_z, "Msun_geo"
    PRINT *, " * Distance of the center of mass of the particle ", &
             "distribution from the  origin: com_d= ", com_d
    PRINT *, " * |com_x-com_star/com_star|=", &
             ABS( com_x-com_star )/ABS( com_star )
    PRINT *

    !-----------------------------------------------------------------------!
    !-- Mirror the positions after having repositioned the center of mass --!
    !-----------------------------------------------------------------------!

    pos_tmp= pos
    nu_tmp= nu
    itr= 0
    DO a= 1, npart_real, 1
      IF( pos_tmp( 3, a ) > 0.0D0 )THEN
        itr= itr + 1
        pos( 1, itr )= pos_tmp( 1, a )
        pos( 2, itr )= pos_tmp( 2, a )
        pos( 3, itr )= pos_tmp( 3, a )
        nu( itr )    = nu_tmp( a )
      ENDIF
    ENDDO
    npart_real_half= itr

    DO a= 1, npart_real_half, 1
      pos( 1, npart_real_half + a )=   pos( 1, a )
      pos( 2, npart_real_half + a )=   pos( 2, a )
      pos( 3, npart_real_half + a )= - pos( 3, a )
      nu( npart_real_half + a )    =   nu( a )
    ENDDO

    PRINT *, "** After mirroring particles:"
    PRINT *, " * x coordinate of the center of mass of the star, ", &
             "from LORENE: com_star= ", com_star, "Msun_geo"
    PRINT *, " * x coordinate of the center of mass of the particle ", &
             "distribution: com_x= ", com_x, "Msun_geo"
    PRINT *, " * y coordinate of the center of mass of the particle ", &
             "distribution: com_y= ", com_y, "Msun_geo"
    PRINT *, " * z coordinate of the center of mass of the particle ", &
             "distribution: com_z= ", com_z, "Msun_geo"
    PRINT *, " * Distance of the center of mass of the particle ", &
             "distribution from the  origin: com_d= ", com_d
    PRINT *, " * |com_x-com_star/com_star|=", &
             ABS( com_x - com_star )/ABS( com_star )
    PRINT *

    IF( PRESENT(namefile_pos) )THEN
      finalnamefile= namefile_pos
    ELSE
      finalnamefile= "apm_pos.dat"
    ENDIF

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
        pos( 1, a ), &
        pos( 2, a ), &
        pos( 3, a )
    ENDDO

    DO a= npart_real + 1, npart_all, 1
      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        2, a, &
        all_pos( 1, a ), &
        all_pos( 2, a ), &
        all_pos( 3, a )
    ENDDO

    CLOSE( UNIT= 2 )

    !-------------------------------------------------------------------!
    !-- now assign baryon number to match profile as good as possible --!
    !-------------------------------------------------------------------!

    PRINT *, " * Assign baryon number..."
    PRINT *

    IF( debug ) PRINT *, "1"

    h      = h(1:npart_real)
    h_guess= h_guess(1:npart_real)
    nu     = nu(1:npart_real)

    CALL assign_h( nn_des, &           !
                   npart_real, &        !
                   pos, h_guess, & ! Input
                   h )                 ! Output

    IF( debug ) PRINT *, "2"

    ! Measure SPH particle number density
    nu= 1.0D0
    CALL density_loop( npart_real, pos, &    ! input
                       nu, h, nstar_real )      ! output

    IF( debug ) PRINT *, "3"

    CALL binary% import_id( npart_real, pos(1,:), &
                            pos(2,:), &
                            pos(3,:), &
                            lapse, shift_x, shift_y, shift_z, &
                            g_xx, g_xy, g_xz, &
                            g_yy, g_yz, g_zz, &
                            baryon_density, &
                            energy_density, &
                            specific_energy, &
                            pressure, &
                            v_euler_x, v_euler_y, v_euler_z )

    IF( debug ) PRINT *, "3.5"

    DO a= 1, npart_real, 1

      IF( baryon_density(a) <= 0 )THEN
        PRINT *, " * The baryon density is 0 at particle ", a
        PRINT *
        STOP
      ENDIF

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
      Theta(a)= Theta_a

      nstar_p(a)= sq_g*Theta_a*baryon_density(a)*((Msun_geo*km2m)**3)/(amu*g2kg)

    ENDDO

    nu= nu_all
    PRINT *, "nu_all= ", nu_all

    nu_tot= 0.0D0
    DO a= 1, npart_real, 1
      nu_tot= nu_tot + nu(a)
    ENDDO

    PRINT *, "nu_tot=", nu_tot
    PRINT *, "mass estimate= ", nu_tot*amu/MSun, "=", &
             100.0D0*nu_tot*amu/MSun/mass, "% of the LORENE baryon mass"
    PRINT *

    IF( debug ) PRINT *, "4"

    PRINT *, "npart_real= ", npart_real
    PRINT *, "SIZE(nu)= ", SIZE(nu)
    PRINT *

    nu_ratio= MAXVAL( nu, DIM= 1 )/MINVAL( nu, DIM= 1 )
    PRINT *, "nu_ratio before correction = ", nu_ratio
    PRINT *

    !----------------------!
    !-- Correcting nu... --!
    !----------------------!

    PRINT *, "Correcting nu"
    PRINT *

    DO a= 1, npart_real, 1

      !IF( a == 1 ) PRINT *, "1"

      nu_tmp2= nu(a)
      nu(a)= nstar_p(a)/nstar_real(a)
      !IF( MAXVAL( nu, DIM= 1 )/MINVAL( nu, DIM= 1 ) > nuratio_thres*nu_ratio )THEN
        IF( nu(a) > nu_tmp2*SQRT(nuratio_thres) ) nu(a)= &
                                          nu_tmp2*SQRT(nuratio_thres)
        IF( nu(a) < nu_tmp2/SQRT(nuratio_thres) ) nu(a)= &
                                          nu_tmp2/SQRT(nuratio_thres)
      !ENDIF

    ENDDO

    nu_ratio= MAXVAL( nu, DIM= 1 )/MINVAL( nu, DIM= 1 )
    PRINT *, "nu_ratio after correction = ", nu_ratio

    max_nu= 0.0D0
    min_nu= HUGE(1.0D0)
    DO a= 1, npart_real, 1
       IF( nu(a) > max_nu )THEN
         max_nu= nu(a)
         a_numax= a
       ENDIF
       IF( nu(a) < min_nu )THEN
         min_nu= nu(a)
         a_numin= a
       ENDIF
    ENDDO

    PRINT *, "Baryon number assigned."
    PRINT *

    PRINT *, "max_nu=", max_nu
    PRINT *, "        at ", pos(:, a_numax), " r= ", &
             NORM2( pos(:, a_numax) )/larger_radius
    PRINT *, "min_nu=", min_nu
    PRINT *, "        at ", pos(:, a_numin), " r= ", &
             NORM2( pos(:, a_numin) )/larger_radius
    PRINT *, "max_nu/min_nu=", max_nu/min_nu
    PRINT *

    IF( mass_it )THEN

      ! just a few iterations to NOT get the nu-ratio too large
      mass_iteration: DO itr= 1, m_max_it, 1

         ! measure density
         CALL density_loop( npart_real, pos, &    ! input
                            nu, h, nstar_real )      ! output


         CALL binary% import_id( npart_real, pos(1,:), &
                                 pos(2,:), &
                                 pos(3,:), &
                                 lapse, shift_x, shift_y, shift_z, &
                                 g_xx, g_xy, g_xz, &
                                 g_yy, g_yz, g_zz, &
                                 baryon_density, &
                                 energy_density, &
                                 specific_energy, &
                                 pressure, &
                                 v_euler_x, v_euler_y, v_euler_z )

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
           Theta(a)= Theta_a

           nstar_p(a)= sq_g*Theta_a*baryon_density(a)*((Msun_geo*km2m)**3)/(amu*g2kg)

         ENDDO

         !nstar_p( npart_real+1:npart_all )= 0.0D0

         ! get RELATIVE nu's right
         dN_av= 0.0D0
         max_nu= 0.0D0
         min_nu= HUGE(1.0D0)
         DO a= 1, npart_real, 1
            dN=    (nstar_real(a)-nstar_p(a))/nstar_p(a)
            nu(a)= nu(a)*(1.0D0 - dN)
            dN_av= dN_av + dN
            IF( nu(a) > max_nu )THEN
              max_nu= nu(a)
              a_numax= a
            ENDIF
            IF( nu(a) < min_nu )THEN
              min_nu= nu(a)
              a_numin= a
            ENDIF
         ENDDO
         dN_av= dN_av/DBLE(npart_real)

         ! exit condition
         IF( dN_av < tol )EXIT

      ENDDO mass_iteration

    ENDIF

    PRINT *, "max_nu=", max_nu
    PRINT *, "        at ", pos(:, a_numax), " r= ", &
             NORM2( pos(:, a_numax) )/larger_radius
    PRINT *, "min_nu=", min_nu
    PRINT *, "        at ", pos(:, a_numin), " r= ", &
             NORM2( pos(:, a_numin) )/larger_radius
    PRINT *, "max_nu/min_nu=", max_nu/min_nu
    PRINT *

    max_nu2= 0.0D0
    min_nu2= HUGE(1.0D0)
    DO a= 1, npart_real, 1
       IF( nu(a) > max_nu2 .AND. a /= a_numax )THEN
         max_nu2= nu(a)
         a_numax2= a
       ENDIF
       IF( nu(a) < min_nu2 .AND. a /= a_numin )THEN
         min_nu2= nu(a)
         a_numin2= a
       ENDIF
    ENDDO

    PRINT *, "Excluding the absolute max and min of nu:"
    PRINT *
    PRINT *, "max_nu=", max_nu2
    PRINT *, "        at ", pos(:, a_numax2), " r= ", &
             NORM2( pos(:, a_numax2) )/larger_radius
    PRINT *, "min_nu=", min_nu2
    PRINT *, "        at ", pos(:, a_numin2), " r= ", &
             NORM2( pos(:, a_numin2) )/larger_radius
    PRINT *, "max_nu/min_nu=", max_nu2/min_nu2
    PRINT *

    nu_tot= 0.0D0
    DO a= 1, npart_real, 1
      nu_tot= nu_tot + nu(a)
    ENDDO
    mean_nu= nu_tot/npart_real

    variance_nu = 0.0                       ! compute variance
    DO a = 1, npart_real, 1
      variance_nu = variance_nu + (nu(a) - mean_nu)**2.0D0
    END DO
    variance_nu = variance_nu / DBLE(npart_real - 1)
    stddev_nu   = SQRT(variance_nu)            ! compute standard deviation

    PRINT *, "mean_nu=", mean_nu
    PRINT *, "variance_nu=", variance_nu
    PRINT *, "stddev_nu=", stddev_nu
    PRINT *, "stddev_nu/mean_nu=", stddev_nu/mean_nu
    PRINT *

    PRINT *, "Before correcting nu to match the mass of the star..."
    PRINT *
    PRINT *, "nu_tot=", nu_tot
    PRINT *, "mass estimate= ", nu_tot*amu/MSun, "=", &
             100.0D0*nu_tot*amu/MSun/mass, "% of the LORENE baryon mass"
    PRINT *

    !----------------------------!
    !-- enforce centre of mass --!
    !----------------------------!

    CALL COM( npart_real, pos, nu(1:npart_real), & ! input
              com_x, com_y, com_z, com_d ) ! output

    IF( debug ) PRINT *, "9"

    PRINT *, "** Before center of mass correction:"
    PRINT *, " * x coordinate of the center of mass of the star, ", &
             "from LORENE: com_star= ", com_star, "Msun_geo"
    PRINT *, " * x coordinate of the center of mass of the particle ", &
             "distribution: com_x= ", com_x, "Msun_geo"
    PRINT *, " * y coordinate of the center of mass of the particle ", &
             "distribution: com_y= ", com_y, "Msun_geo"
    PRINT *, " * z coordinate of the center of mass of the particle ", &
             "distribution: com_z= ", com_z, "Msun_geo"
    PRINT *, " * Distance of the center of mass of the particle ", &
             "distribution from the  origin: com_d= ", com_d
    PRINT *, " * |com_x-com_star/com_star|=", &
             ABS( com_x-com_star )/ABS( com_star )
    PRINT *

    DO a= 1, npart_real, 1

      pos_corr_tmp(1)= pos(1,a) - ( com_x - com_star )
      pos_corr_tmp(2)= pos(2,a) - com_y
      pos_corr_tmp(3)= pos(3,a) - com_z

      IF( binary% import_mass_density( &
                  pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) > 0.0D0 &
          .AND. &
          binary% is_hydro_negative( &
                  pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) == 0 &
      )THEN

        pos(:,a)= pos_corr_tmp

      ENDIF

    ENDDO

    IF( debug ) PRINT *, "10"

    CALL COM( npart_real, pos, nu(1:npart_real), & ! input
              com_x, com_y, com_z, com_d ) ! output

    IF( debug ) PRINT *, "9"

    PRINT *, "** After center of mass correction:"
    PRINT *, " * x coordinate of the center of mass of the star, ", &
             "from LORENE: com_star= ", com_star, "Msun_geo"
    PRINT *, " * x coordinate of the center of mass of the particle ", &
             "distribution: com_x= ", com_x, "Msun_geo"
    PRINT *, " * y coordinate of the center of mass of the particle ", &
             "distribution: com_y= ", com_y, "Msun_geo"
    PRINT *, " * z coordinate of the center of mass of the particle ", &
             "distribution: com_z= ", com_z, "Msun_geo"
    PRINT *, " * Distance of the center of mass of the particle ", &
             "distribution from the  origin: com_d= ", com_d
    PRINT *, " * |com_x-com_star/com_star|=", &
             ABS( com_x-com_star )/ABS( com_star )
    PRINT *

    !-----------------------------------------------------------------------!
    !-- Mirror the positions after having repositioned the center of mass --!
    !-----------------------------------------------------------------------!

    pos_tmp= pos
    nu_tmp= nu
    itr= 0
    DO a= 1, npart_real, 1
      IF( pos_tmp( 3, a ) > 0.0D0 )THEN
        itr= itr + 1
        pos( 1, itr )= pos_tmp( 1, a )
        pos( 2, itr )= pos_tmp( 2, a )
        pos( 3, itr )= pos_tmp( 3, a )
        nu( itr )    = nu_tmp( a )
      ENDIF
    ENDDO
    npart_real_half= itr

    DO a= 1, npart_real_half, 1
      pos( 1, npart_real_half + a )=   pos( 1, a )
      pos( 2, npart_real_half + a )=   pos( 2, a )
      pos( 3, npart_real_half + a )= - pos( 3, a )
      nu( npart_real_half + a )    =   nu( a )
    ENDDO

    PRINT *, "** After mirroring particles:"
    PRINT *, " * x coordinate of the center of mass of the star, ", &
             "from LORENE: com_star= ", com_star, "Msun_geo"
    PRINT *, " * x coordinate of the center of mass of the particle ", &
             "distribution: com_x= ", com_x, "Msun_geo"
    PRINT *, " * y coordinate of the center of mass of the particle ", &
             "distribution: com_y= ", com_y, "Msun_geo"
    PRINT *, " * z coordinate of the center of mass of the particle ", &
             "distribution: com_z= ", com_z, "Msun_geo"
    PRINT *, " * Distance of the center of mass of the particle ", &
             "distribution from the  origin: com_d= ", com_d
    PRINT *, " * |com_x-com_star/com_star|=", &
             ABS( com_x - com_star )/ABS( com_star )
    PRINT *

    IF( correct_nu )THEN
      !THIS% nlrf=
      !        THIS% nlrf/(THIS% nbar_tot*amu/Msun/(THIS% mass1 + THIS% mass2))
      !nlrf= nlrf/(THIS% nbar_tot*amu/Msun/(THIS% mass1 + THIS% mass2))
      !THIS% nu= THIS% nu/(THIS% nbar_tot*amu/Msun/(THIS% mass1 + THIS% mass2))
      !nu= nu/(THIS% nbar_tot*amu/Msun/(THIS% mass1 + THIS% mass2))
      !THIS% nbar_tot= &
      !    THIS% nbar_tot/(THIS% nbar_tot*amu/Msun/(THIS% mass1 + THIS% mass2))
      nu= nu/(nu_tot*amu/Msun/mass)
      nu_tot= 0.0D0
      DO a= 1, npart_real, 1
        nu_tot= nu_tot + nu(a)
      ENDDO

      PRINT *, "After correcting nu to match the mass of the star..."
      PRINT *
      PRINT *, "nu_tot=", nu_tot
      PRINT *, "mass estimate= ", nu_tot*amu/MSun, "=", &
               100.0D0*nu_tot*amu/MSun/mass, "% of the LORENE baryon mass"
      PRINT *

      !----------------------------!
      !-- enforce centre of mass --!
      !----------------------------!

      CALL COM( npart_real, pos, nu(1:npart_real), & ! input
                com_x, com_y, com_z, com_d ) ! output

      IF( debug ) PRINT *, "9"

      PRINT *, "** Before center of mass correction:"
      PRINT *, " * x coordinate of the center of mass of the star, ", &
               "from LORENE: com_star= ", com_star, "Msun_geo"
      PRINT *, " * x coordinate of the center of mass of the particle ", &
               "distribution: com_x= ", com_x, "Msun_geo"
      PRINT *, " * y coordinate of the center of mass of the particle ", &
               "distribution: com_y= ", com_y, "Msun_geo"
      PRINT *, " * z coordinate of the center of mass of the particle ", &
               "distribution: com_z= ", com_z, "Msun_geo"
      PRINT *, " * Distance of the center of mass of the particle ", &
               "distribution from the  origin: com_d= ", com_d
      PRINT *, " * |com_x-com_star/com_star|=", &
               ABS( com_x-com_star )/ABS( com_star )
      PRINT *

      DO a= 1, npart_real, 1

        pos_corr_tmp(1)= pos(1,a) - ( com_x - com_star )
        pos_corr_tmp(2)= pos(2,a) - com_y
        pos_corr_tmp(3)= pos(3,a) - com_z

        IF( binary% import_mass_density( &
                    pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) > 0.0D0 &
            .AND. &
            binary% is_hydro_negative( &
                    pos_corr_tmp(1), pos_corr_tmp(2), pos_corr_tmp(3) ) == 0 &
        )THEN

          pos(:,a)= pos_corr_tmp

        ENDIF

      ENDDO

      IF( debug ) PRINT *, "10"

      CALL COM( npart_real, pos, nu(1:npart_real), & ! input
                com_x, com_y, com_z, com_d ) ! output

      IF( debug ) PRINT *, "9"

      PRINT *, "** After center of mass correction:"
      PRINT *, " * x coordinate of the center of mass of the star, ", &
               "from LORENE: com_star= ", com_star, "Msun_geo"
      PRINT *, " * x coordinate of the center of mass of the particle ", &
               "distribution: com_x= ", com_x, "Msun_geo"
      PRINT *, " * y coordinate of the center of mass of the particle ", &
               "distribution: com_y= ", com_y, "Msun_geo"
      PRINT *, " * z coordinate of the center of mass of the particle ", &
               "distribution: com_z= ", com_z, "Msun_geo"
      PRINT *, " * Distance of the center of mass of the particle ", &
               "distribution from the  origin: com_d= ", com_d
      PRINT *, " * |com_x-com_star/com_star|=", &
               ABS( com_x-com_star )/ABS( com_star )
      PRINT *

      !-----------------------------------------------------------------------!
      !-- Mirror the positions after having repositioned the center of mass --!
      !-----------------------------------------------------------------------!

      pos_tmp= pos
      nu_tmp= nu
      itr= 0
      DO a= 1, npart_real, 1
        IF( pos_tmp( 3, a ) > 0.0D0 )THEN
          itr= itr + 1
          pos( 1, itr )= pos_tmp( 1, a )
          pos( 2, itr )= pos_tmp( 2, a )
          pos( 3, itr )= pos_tmp( 3, a )
          nu( itr )    = nu_tmp( a )
        ENDIF
      ENDDO
      npart_real_half= itr

      DO a= 1, npart_real_half, 1
        pos( 1, npart_real_half + a )=   pos( 1, a )
        pos( 2, npart_real_half + a )=   pos( 2, a )
        pos( 3, npart_real_half + a )= - pos( 3, a )
        nu( npart_real_half + a )    =   nu( a )
      ENDDO

    ENDIF

    !-------------------!
    !-- monitoring... --!
    !-------------------!

    ! measure density
    CALL density_loop( npart_real, pos, &    ! input
                       nu, h, nstar_real )      ! output

    CALL binary% import_id( npart_real, pos(1,:), &
                            pos(2,:), &
                            pos(3,:), &
                            lapse, shift_x, shift_y, shift_z, &
                            g_xx, g_xy, g_xz, &
                            g_yy, g_yz, g_zz, &
                            baryon_density, &
                            energy_density, &
                            specific_energy, &
                            pressure, &
                            v_euler_x, v_euler_y, v_euler_z )

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
      Theta(a)= Theta_a

      nstar_p(a)= sq_g*Theta_a*baryon_density(a)*((Msun_geo*km2m)**3)/(amu*g2kg)

      IF( ISNAN( nstar_p( a ) ) )THEN
        PRINT *, "** ERROR! nstar_p(", a, ") is a NaN!", &
                 " Stopping.."
        PRINT *
        STOP
      ENDIF
      IF( nstar_p( a ) == 0 )THEN
        PRINT *, "** ERROR! When setting up ID: ", &
                 "nstar_p(", a, ")= 0 on a real particle!"
        PRINT *, " * Particle position: x=", all_pos(1,a), &
                 ", y=", all_pos(2,a), ", z=", all_pos(3,a)
        PRINT *, "   sq_g=", sq_g
        PRINT *, "   Theta_a=", Theta_a
        PRINT *, "   baryon_density(a)=", baryon_density(a)
        PRINT *, " * Stopping.."
        PRINT *
        STOP
      ENDIF

    ENDDO

    !nstar_p( npart_real+1:npart_all )= 0.0D0

    ! get RELATIVE nu's right
    dN_av= 0.0D0
        dN_max= 0.0D0
    DO a= 1, npart_real, 1
      dN=     ABS(nstar_real(a)-nstar_p(a))/nstar_p(a)
      dN_max= MAX(dN_max,dN)
      dN_av=  dN_av + dN
    ENDDO
    dN_av= dN_av/DBLE(npart_real)
    PRINT *,'...dN_max ', dN_max
    PRINT *,'...dN_av  ', dN_av
    PRINT *

    IF( debug ) PRINT *, "100"

    IF( .NOT.ALLOCATED( nstar_int ) ) ALLOCATE( nstar_int( npart_real ) )

    IF( debug ) PRINT *, "101"

    CALL exact_nei_tree_update( nn_des, &
                                npart_real, &
                                pos, nu )

    !IF( mass == THIS% mass2 ) STOP

    IF( debug ) PRINT *, "102"

    CALL density( npart_real, pos, nstar_int )

    IF( debug ) PRINT *, "103"

  !  PRINT *, "** Finding nearest neighbors..."

    ALLOCATE( neighbors_lists( npart_real ) )
    ALLOCATE( n_neighbors( npart_real ) )
    ALLOCATE( nearest_neighbors( 2, npart_real ) )

    neighbors_lists= 0
    n_neighbors= 0
    nearest_neighbors(1,:)= 0
    nearest_neighbors(2,:)= HUGE(1.0D0)

 !   find_neighbors: DO a= 1, npart_real, 1
 !
 !     CALL get_neighbours_bf( a, npart_real, pos, h, 3, &           ! Input
 !                             n_neighbors(a), neighbors_lists(:) )  ! Output
 !
 !     DO itr= 1, n_neighbors(a), 1
 !
 !       dist= NORM2( pos(:,a) - pos(:,neighbors_lists(itr)) )
 !
 !       !PRINT *, "dist= ", dist
 !       !PRINT *, "nearest_neighbors(2,a)= ", nearest_neighbors(2,a)
 !       !PRINT *, "dist < nearest_neighbors(2,a)= ", dist < nearest_neighbors(2,a)
 !
 !       IF( dist < nearest_neighbors(2,a) )THEN
 !
 !         nearest_neighbors(1,a)= neighbors_lists(itr)
 !         nearest_neighbors(2,a)= dist
 !
 !       ENDIF
 !
 !     ENDDO
 !
 !     !PRINT *, "dist= ", dist
 !     !PRINT *, "nearest_neighbors(2,a)= ", nearest_neighbors(2,a)
 !     !PRINT *
 !
 !   ENDDO find_neighbors
 !
 !   PRINT *, " * Nearest neighbors found. "
 !   PRINT *, " * Average number of neighbors= ", DBLE(SUM(n_neighbors))/DBLE(npart_real)
 !   PRINT *

    IF( debug ) PRINT *, "0"

    IF( PRESENT(namefile_results) )THEN
      finalnamefile= namefile_results
    ELSE
      finalnamefile= "apm_results.dat"
    ENDIF

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

    IF( debug ) PRINT *, "1"

    IF( .NOT.ALLOCATED( nu_one ) ) ALLOCATE( nu_one( npart_real ) )
    IF( .NOT.ALLOCATED( particle_density_final ) ) &
      ALLOCATE( particle_density_final( npart_real ) )
    nu_one= 1.0D0
    CALL density_loop( npart_real, pos, &    ! input
                       nu_one, h, particle_density_final )      ! output

    DO a= 1, npart_real, 1
      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        a, &
        pos( 1, a ), pos( 2, a ), pos( 3, a ), &
        nstar_p( a ), &
        nstar_int( a ), &
        particle_density_final( a ), &
        particle_density_final( a )*nstar_p( 1 )/particle_density_final( 1 ), &
        ABS(nstar_real(a)-nstar_p(a))/nstar_p(a), &
        nu(a), &
        nearest_neighbors(2,a)
    ENDDO

    CLOSE( UNIT= 2 )

    IF( debug ) PRINT *, "2"

    pos_input= pos
    h_output = h
    nu_output= nu

    IF( debug ) PRINT *, "3"

  !  PRINT *, ALLOCATED(pos_u)
  !  PRINT *, ALLOCATED(vel_u)
  !  PRINT *, ALLOCATED(S_l)
  !  PRINT *, ALLOCATED(u)
  !  PRINT *, ALLOCATED(h)
  !  PRINT *, ALLOCATED(Nstar)
  !  PRINT *, ALLOCATED(nlrf)
  !  PRINT *, ALLOCATED(temp)
  !  PRINT *, ALLOCATED(Ye)
  !  PRINT *, ALLOCATED(nu)
  !  PRINT *, ALLOCATED(ehat)
  !  PRINT *, ALLOCATED(theta)
  !  PRINT *, ALLOCATED(Pr)
  !  PRINT *, ALLOCATED(cs)
  !  PRINT *, ALLOCATED(av)
  !  PRINT *, ALLOCATED(Kent)
  !  PRINT *, ALLOCATED(divv)
  !  PRINT *

    DEALLOCATE( posmash )
    CALL deallocate_metric_on_particles()
    IF( debug ) PRINT *, "4"
    CALL deallocate_gradient()
    IF( debug ) PRINT *, "5"
    CALL deallocate_RCB_tree_memory_3D()
    IF( debug ) PRINT *, "6"
    CALL deallocate_SPH_memory()

    !STOP

  END PROCEDURE perform_apm


  SUBROUTINE get_neighbours_bf(ipart,npart,pos,h,dimensions,nnei,neilist)

    !**************************************************************
    !                                                             *
    ! just for test purposes: get neighbours of particle ipart in *
    ! a "brute force" way; ipart is ALSO on the neighbour list;   *
    ! SKR 8.2.2010                                                *
    !                                                             *
    ! Removed ipart from its own neighbors' list                  *
    ! FT 04.06.2021                                               *
    !                                                             *
    !**************************************************************

    IMPLICIT NONE

    INTEGER,INTENT(IN)::          ipart,npart,dimensions
    DOUBLE PRECISION,INTENT(IN):: pos(dimensions,npart),h(npart)
    INTEGER,INTENT(OUT)::         nnei,neilist(npart)
    INTEGER a
    DOUBLE PRECISION diff(dimensions),d2,r_int2

    ! square of interaction radius
    r_int2= (2.0D0*h(ipart))**2

    nnei= 0
  !    !$OMP PARALLEL DO SHARED(pos,dimensions,ipart,npart,r_int2,nnei,neilist)&
  !    !$OMP             PRIVATE(a,diff,d2)
    DO a= 1, npart, 1

      IF( a /= ipart )THEN

        diff= pos(1:dimensions,a)-pos(1:dimensions,ipart)
        d2= DOT_PRODUCT(diff,diff)

        ! neighbour?
        IF(d2 < r_int2)THEN
          nnei= nnei + 1
          neilist(nnei)= a
        ENDIF

      ENDIF

    ENDDO
  !    !$OMP END PARALLEL DO

  END SUBROUTINE get_neighbours_bf


END SUBMODULE particles_apm
