! File:         proto_artificial-pressure.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

PROGRAM proto_apm

  !*****************************************************
  !                                                    *
  ! Prototype for the development of the artificial    *
  ! pressure method for the LORENE ID                  *
  !                                                    *
  ! FT 10.05.2021                                      *
  !                                                    *
  !*****************************************************

  USE bns_id,              ONLY: bns
  !USE particles_id,        ONLY: particles
  USE RCB_tree_3D,         ONLY: allocate_RCB_tree_memory_3D,&
                                 iorig
  USE kernel_table,        ONLY: ktable
  USE sph_variables,       ONLY: Rstar,divv,av,Pr,ye,temp,nlrf,&
                                 u,Theta,vel_u,tterm,tgrav,tkin,&
                                 escap,t,n1,n2,pos_u,h,nu,npart,&
                                 npm,Nstar,allocate_sph_memory
  USE options,             ONLY: basename,ikernel,ndes,av_max,icow
  USE input_output,        ONLY: read_options
  USE units,               ONLY: umass,m0c2_CU,set_units

  USE APM,                 ONLY: density_loop, position_correction, assign_h
  USE set_h,               ONLY: exact_nei_tree_update
  USE analyze,             ONLY: COM
  USE matrix,              ONLY: determinant_4x4_matrix
  USE constants,           ONLY: MSun, Msun_geo, km2m, g2kg, amu, pi, half, &
                                 m0c2, lorene2hydrobase, c_light
  USE NR,                  ONLY: indexx
  USE metric_on_particles, ONLY: allocate_metric_on_particles, &
                                 deallocate_metric_on_particles
  USE gradient,            ONLY: allocate_gradient, deallocate_gradient
  USE alive_flag,          ONLY: alive
  USE set_h,               ONLY: exact_nei_tree_update

  USE pwp_EOS,             ONLY: select_EOS_parameters, gen_pwp_eos
  USE sphincs_sph,         ONLY: density

  IMPLICIT NONE

  INTEGER, PARAMETER:: unit_id     = 23
  INTEGER, PARAMETER:: max_npart   = 5D+6
  INTEGER, PARAMETER:: apm_max_it  = 1500
  INTEGER, PARAMETER:: m_max_it    = 50
  INTEGER, PARAMETER:: max_inc     = 50
  INTEGER, PARAMETER:: nn_des      = 301
  DOUBLE PRECISION, PARAMETER:: tol= 1.0D-3
  DOUBLE PRECISION, PARAMETER:: iter_tol= 1.0D-2
  LOGICAL, PARAMETER:: run_apm     = .FALSE.
  LOGICAL, PARAMETER:: post_correction= .FALSE.

  INTEGER:: npart_real, tmp, ios, itr, itr2, a, nout, nus, mus, npart_eq, &
            npart_ghost_shell, npart_ghost, npart_all, r, th, phi, npart1, &
            n_inc, nx, ny, nz, i, j, k, tmp2, a_numax, a_numin, &
            a_numax2, a_numin2
  INTEGER, DIMENSION(:), ALLOCATABLE:: x_sort, xy_sort, xyz_sort, lim

  DOUBLE PRECISION:: com_x, com_y, com_z, com_d, det, sq_g, Theta_a, &
                     nu_corr, max_nu, min_nu, smaller_radius, larger_radius, &
                     dr, rad, col, col_tmp, alpha, phase_th, rand_num, &
                     center, long, xtemp, ytemp, ztemp, radius_z, rel_sign, &
                     phase_phi, h_max, h_av, dNstar, art_pr_max, err_N_max , &
                     err_N_min, err_N_mean, nu_all, mass_star, err_mean_old, &
                     rad_x, rad_y, rad_z, radius_y, nstar_real_err, &
                     nstar_p_err, dx, dy, dz, xmin, xmax, ymin, ymax, &
                     zmin, zmax, eps, x_ell, y_ell, z_ell, delta, dN, dN_av, &
                     dN_max, spec_ener, press, mass_dens, tmp3, tmp4, nu_tot, &
                     mass, mean_nu, variance_nu, stddev_nu, max_nu2, min_nu2, &
                     nu_tmp, nu_ratio, err_N_mean_min
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: lapse, &
                     shift_x, shift_y, shift_z, &
                     g_xx, g_xy, g_xz, &
                     g_yy, g_yz, g_zz, &
                     baryon_density, &
                     energy_density, &
                     specific_energy, &
                     pressure, &
                     v_euler_x, v_euler_y, v_euler_z
  DOUBLE PRECISION, DIMENSION(3):: pos_tmp
  DOUBLE PRECISION, DIMENSION(3):: pos_maxerr
  DOUBLE PRECISION:: g4(0:3,0:3)
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: lapse_parts
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: shift_parts_x
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: shift_parts_y
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: shift_parts_z
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: baryon_density_parts
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: energy_density_parts
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: specific_energy_parts
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_parts
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: v_euler_parts_x
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: v_euler_parts_y
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: v_euler_parts_z
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: v
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nu0
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nlrf0
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: Ye0
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: Theta0
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: sph_density
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar0
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: h0

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos_star1
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: sorted_pos
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: ghost_pos
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: all_pos
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: correction_pos
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar_real
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar_p
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: art_pr
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: freeze
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: h_guess
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar_int
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nu_one
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: particle_density_final

  LOGICAL:: exist

  CHARACTER( LEN= : ), ALLOCATABLE:: namefile
  CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile
  CHARACTER( LEN= : ), ALLOCATABLE:: err_msg

  TYPE( bns ):: binary

  !---------------------------!
  !--  End of declarations  --!
  !---------------------------!

  !CALL DATE_AND_TIME( date, time, zone, values )
  !run_id= date // "-" // time

!  CALL select_EOS_parameters( 'APR4' )
!
!  finalnamefile= "apr4_sphincs.dat"
!
!  INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )
!
!  IF( exist )THEN
!      OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
!            FORM= "FORMATTED", &
!            POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
!            IOMSG= err_msg )
!  ELSE
!      OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "NEW", &
!            FORM= "FORMATTED", &
!            ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
!  ENDIF
!  IF( ios > 0 )THEN
!    PRINT *, "...error when opening " // TRIM(finalnamefile), &
!             ". The error message is", err_msg
!    STOP
!  ENDIF
!
!  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
!    "   baryon mass density [g cm^{-3}]    pressure [dyne cm^{-2}]", &
!    "   specific energy (or internal energy) [c^2]"
!
!  DO a= 1, 50000, 1
!
!    mass_dens= ( 1.0D10 + (DBLE(a)/5.0D5)*( 1.0D20 - 1.0D10 ) )*lorene2hydrobase
!
!    CALL gen_pwp_eos( mass_dens, spec_ener, press, 0.0D0, tmp3, tmp4 )
!
!    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
!      (mass_dens/lorene2hydrobase)/1000.0D0, &
!      (press/lorene2hydrobase*(c_light/100.0D0)**2.0D0)*10, &
!      spec_ener
!
!  ENDDO
!
!  CLOSE( UNIT= 2 )
!
!  STOP

  !PRINT *, "Hello world!"
  PRINT *, "** Beginning of PROGRAM proto_apm."
  PRINT *

  namefile= "poly2-75_1.2-1.8_45km.bin"

  binary= bns( namefile )

  ALLOCATE( pos(3,max_npart) )
  ALLOCATE( lapse_parts(max_npart) )
  ALLOCATE( shift_parts_x(max_npart) )
  ALLOCATE( shift_parts_y(max_npart) )
  ALLOCATE( shift_parts_z(max_npart) )
  ALLOCATE( baryon_density_parts(max_npart) )
  ALLOCATE( energy_density_parts(max_npart) )
  ALLOCATE( specific_energy_parts(max_npart) )
  ALLOCATE( pressure_parts(max_npart) )
  ALLOCATE( v_euler_parts_x(max_npart) )
  ALLOCATE( v_euler_parts_y(max_npart) )
  ALLOCATE( v_euler_parts_z(max_npart) )
  ALLOCATE( v(3,max_npart) )
  ALLOCATE( nu0(max_npart) )
  ALLOCATE( nlrf0(max_npart) )
  ALLOCATE( Ye0(max_npart) )
  ALLOCATE( Theta0(max_npart) )
  ALLOCATE( sph_density(max_npart) )
  ALLOCATE( nstar0(max_npart) )
  ALLOCATE( h0(max_npart) )

  finalnamefile= "lorene-bns-id-particles-shells.dat"

  INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

  IF( exist )THEN
    OPEN( UNIT= unit_id, FILE= TRIM(finalnamefile), &
          FORM= "FORMATTED", ACTION= "READ", IOSTAT= ios, &
          IOMSG= err_msg )
    IF( ios > 0 )THEN
      PRINT *, "...error when opening " // TRIM(finalnamefile), &
              ". The error message is", err_msg
      STOP
    ENDIF
  ELSE
    PRINT *, "** ERROR! Unable to find file " // TRIM(finalnamefile)
    STOP
  ENDIF

  PRINT *, " * Reading file " // TRIM(finalnamefile) // "..."
  npart_real= 0
  ! Skipping the first 4 header rows
  READ( unit_id, * )
  READ( unit_id, * )
  READ( unit_id, * )
  READ( unit_id, * )
  ! Reading particle data
  DO itr= 1, max_npart, 1
    READ( UNIT= unit_id, FMT= *, IOSTAT = ios, IOMSG= err_msg ) &
      tmp, &
      pos( 1, itr ), &
      pos( 2, itr ), &
      pos( 3, itr ), &
      lapse_parts( itr ), &
      shift_parts_x( itr ), &
      shift_parts_y( itr ), &
      shift_parts_z( itr ), &
      baryon_density_parts( itr ), &
      energy_density_parts( itr ), &
      specific_energy_parts( itr ), &
      pressure_parts( itr ), &
      v_euler_parts_x( itr ), &
      v_euler_parts_y( itr ), &
      v_euler_parts_z( itr ), &
      v( 1, itr ), &
      v( 2, itr ), &
      v( 3, itr ), &
      nu0( itr ), &
      nlrf0( itr ), &
      Ye0( itr ), &
      Theta0( itr ), &
      sph_density( itr ), &
      nstar0( itr ), &
      h0( itr )
    IF( ios > 0 )THEN
      PRINT *, "...error when reading " // TRIM(finalnamefile), &
              " at step ", itr,". The status variable is ", ios, &
              ". The error message is", err_msg
      STOP
    ENDIF
    IF( ios < 0 )THEN
      PRINT *, " * Reached end of file " // TRIM(finalnamefile)
      PRINT *
      EXIT
    ENDIF
    npart_real= npart_real + 1
  ENDDO
  IF( tmp /= npart_real )THEN
    PRINT *, "** ERROR! Mismatch in the values of the particle number read", &
             " from " // TRIM(finalnamefile), " in two ways. "
    PRINT *, "tmp= ", tmp, ", npart_real= ", npart_real
    PRINT *, "Stopping..."
    STOP
  ENDIF

  CLOSE( unit_id )

  pos                  = pos                  ( :, 1:npart_real )
  lapse_parts          = lapse_parts          ( 1:npart_real )
  shift_parts_x        = shift_parts_x        ( 1:npart_real )
  shift_parts_y        = shift_parts_y        ( 1:npart_real )
  shift_parts_z        = shift_parts_z        ( 1:npart_real )
  baryon_density_parts = baryon_density_parts ( 1:npart_real )
  energy_density_parts = energy_density_parts ( 1:npart_real )
  specific_energy_parts= specific_energy_parts( 1:npart_real )
  pressure_parts       = pressure_parts       ( 1:npart_real )
  v_euler_parts_x      = v_euler_parts_x      ( 1:npart_real )
  v_euler_parts_y      = v_euler_parts_y      ( 1:npart_real )
  v_euler_parts_z      = v_euler_parts_z      ( 1:npart_real )
  v                    = v                    ( :, 1:npart_real )
  nu0                  = nu0                  ( 1:npart_real )
  nlrf0                = nlrf0                ( 1:npart_real )
  Ye0                  = Ye0                  ( 1:npart_real )
  Theta0               = Theta0               ( 1:npart_real )
  sph_density          = sph_density          ( 1:npart_real )
  nstar0               = nstar0               ( 1:npart_real )
  h0                   = h0                   ( 1:npart_real )

  !---------------------------!
  !-- Place ghost particles --!
  !---------------------------!

  !PRINT *, "Sorting positions"
  !PRINT *
  !
  !finalnamefile= "sorted_pos.dat"
  !
  !INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )
  !
  !IF( exist )THEN
  !    OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
  !          FORM= "FORMATTED", &
  !          POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
  !          IOMSG= err_msg )
  !ELSE
  !    OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "NEW", &
  !          FORM= "FORMATTED", &
  !          ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
  !ENDIF
  !IF( ios > 0 )THEN
  !  PRINT *, "...error when opening " // TRIM(finalnamefile), &
  !           ". The error message is", err_msg
  !  STOP
  !ENDIF
  !
  !
  !ALLOCATE( sorted_pos(3,npart_real) )
  !ALLOCATE( x_sort(npart_real) )
  !ALLOCATE( xy_sort(npart_real) )
  !ALLOCATE( xyz_sort(npart_real) )
  !sorted_pos= 0.0D0
  !x_sort    = 0.0D0
  !xy_sort   = 0.0D0
  !xyz_sort  = 0.0D0
  !
  !PRINT *, sorted_pos(:,1)
  !PRINT *, sorted_pos(:,npart_real)
  !PRINT *
  !
  !! Sort particle positions along x
  !CALL indexx( npart_real, pos( 1, : ), x_sort )
  !
  !DO a= 1, npart_real, 1
  !  sorted_pos(1,a)= pos(1,x_sort(a))
  !  sorted_pos(2,a)= pos(2,x_sort(a))
  !  sorted_pos(3,a)= pos(3,x_sort(a))
  !ENDDO
  !PRINT *, "After sorting x"
  !
  !! Why do you need to subsort the posiions over y and z? x is perhaps enough
  !
  !!PRINT *, sorted_pos(:,1)
  !!PRINT *, sorted_pos(:,npart_real)
  !!PRINT *
  !DO a= 1, npart_real, 1
  !  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
  !    a, &
  !    sorted_pos( 1, a ), &
  !    sorted_pos( 2, a ), &
  !    sorted_pos( 3, a )
  !ENDDO

  ! Sort along y
  !ALLOCATE( lim(npart_real) )
  !lim(0)= 0
  !itr= 1
  !DO a= 1, npart_real - 1, 1
  !  IF( sorted_pos(1,a) /= sorted_pos(1,a + 1) )THEN
  !    lim(itr)= a
  !    CALL indexx( lim(itr)-lim(itr-1), pos( 2, lim(itr-1)+1:lim(itr) ), &
  !                 xy_sort(lim(itr-1)+1:lim(itr)) )
  !    xy_sort(lim(itr-1)+1:lim(itr))= xy_sort(lim(itr-1)+1:lim(itr)) + lim(itr-1)
  !    itr= itr + 1
  !  ENDIF
  !ENDDO
  !!PRINT *, x_sort
  !DO a= 1, npart_real, 1
  !  sorted_pos(1,a)= sorted_pos(1,xy_sort(a))
  !  sorted_pos(2,a)= sorted_pos(2,xy_sort(a))
  !  sorted_pos(3,a)= sorted_pos(3,xy_sort(a))
  !ENDDO
  !
  !DO a= 1, npart_real, 1
  !  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
  !    a, &
  !    sorted_pos( 1, a ), &
  !    sorted_pos( 2, a ), &
  !    sorted_pos( 3, a )
  !ENDDO
  !
  !!DO a= 1, itr - 1, 1
  !!  DO itr2= 1, lim(itr)-lim(itr-1), 1
  !!    sorted_pos(1,itr2)= pos(1,xy_sort(itr2))
  !!  ENDDO
  !!ENDDO
  !PRINT *, "After sorting y"
  !PRINT *, sorted_pos(:,1)
  !PRINT *, sorted_pos(:,npart_real)
  !PRINT *

  ! Sort along z
  !lim(0)= 0
  !itr= 1
  !DO a= 1, npart_real - 1, 1
  !  IF( sorted_pos(2,a) /= sorted_pos(2,a + 1) )THEN
  !    lim(itr)= a
  !    CALL indexx( lim(itr)-lim(itr-1), pos( 3, lim(itr-1)+1:lim(itr) ), &
  !                 xyz_sort(lim(itr-1)+1:lim(itr)) )
  !    xyz_sort(lim(itr-1)+1:lim(itr))= xyz_sort(lim(itr-1)+1:lim(itr)) + lim(itr-1)
  !    itr= itr + 1
  !  ENDIF
  !ENDDO
  !DO a= 1, npart_real, 1
  !  sorted_pos(1,a)= pos(1,xyz_sort(a))
  !  sorted_pos(2,a)= pos(2,xyz_sort(a))
  !  sorted_pos(3,a)= pos(3,xyz_sort(a))
  !ENDDO
  !PRINT *, "After sorting z"
  !PRINT *, sorted_pos(:,1)
  !PRINT *, sorted_pos(:,npart_real)
  !PRINT *
  !
  !PRINT *, "End of sorting algorithm"
  !PRINT *

  !CLOSE( UNIT= 2 )

  ALLOCATE( ghost_pos( 3, max_npart ) )
  ALLOCATE( pos_star1( 3, npart_real ) )
  ghost_pos= 0.0D0
  pos_star1= 0.0D0
  itr= 0
  DO a= 1, npart_real, 1
    IF( pos( 1, a ) < 0.0D0 )THEN
      itr= itr + 1
      pos_star1( 1, itr )= pos( 1 ,a )
      pos_star1( 2, itr )= pos( 2 ,a )
      pos_star1( 3, itr )= pos( 3 ,a )
    ENDIF
  ENDDO
  npart_real= itr
  pos_star1 = pos_star1( :, 1:npart_real )
  !smaller_radius= binary% get_radius1_x_opp()
  !larger_radius = binary% get_radius1_x_comp()
  mass          = binary% get_mass1()
  mass_star     = binary% get_mass1()
  center        = binary% get_center1_x()
  smaller_radius= ABS( MINVAL( pos_star1( 1, : ), DIM= 1 ) - center )
  larger_radius = ABS( center - MAXVAL( pos_star1( 1, : ), DIM= 1 ) )
  radius_y= ABS( MAXVAL( pos_star1( 2, : ), DIM= 1 ) )
  radius_z= ABS( MAXVAL( pos_star1( 3, : ), DIM= 1 ) )
  dr= ( larger_radius - smaller_radius )/15.0D0
  npart_eq= 110
  alpha= 2*pi/DBLE(npart_eq)
  !npart_ghost_shell= ( npart_eq**2 )/2
  itr= 1

  PRINT *, "center=", center
  PRINT *, "smaller_radius=", smaller_radius
  PRINT *, "larger_radius=", larger_radius
  PRINT *, "radius_y=", radius_y
  PRINT *, "radius_z=", radius_z
  PRINT *

  !PRINT *, MAXVAL( pos_star1( 1, : ), DIM= 1 ), &
  !         smaller_radius, larger_radius, center, dr, alpha

  ! You need to place ghost particles one smoothing length outside of the
  ! surface

  ! Find the maximum of the smoothing lencth of the particles whose distance form thecenter is higher than radius_z, and place particles outside of that

  h_max= 0.0D0
  h_av = 0.0D0
  itr  = 0
  DO a= 1, npart_real, 1
    IF( SQRT( pos_star1( 1, a )**2.0D0 + pos_star1( 2, a )**2.0D0 + pos_star1( 3, a )**2.0D0 ) &
        > 0.99D0*radius_z )THEN
      itr= itr + 1
      IF( h0(a) > h_max )THEN
        h_max= h0(a)
      ENDIF
      h_av= h_av + h0(a)
    ENDIF
  ENDDO
  h_av= h_av/itr

!  PRINT *, " * Placing ghost particles..."
!
!  itr= 1
!  DO r= 1, 30, 1
!
!    rad_x= smaller_radius + 2*DBLE( r - 1 )*dr !+ h_av
!    rad_y= radius_y + 2*DBLE( r - 1 )*dr !+ h_av
!    rad_z= radius_z + 2*DBLE( r - 1 )*dr !+ h_av
!
!    DO th= 1, npart_eq/2, 1
!
!      CALL RANDOM_NUMBER( phase_phi )
!      phase_phi= phase_phi*alpha
!
!      col= ( th - 1 )*alpha + alpha/2.0D0
!      CALL RANDOM_NUMBER( phase_th )
!      CALL RANDOM_NUMBER( rand_num )
!      IF( rand_num >= half ) rel_sign=  1
!      IF( rand_num < half )  rel_sign= -1
!
!      col_tmp= col*( 1.0D0 + rel_sign*0.05D0*phase_th )
!
!      IF( col_tmp < pi .AND. col_tmp > 0 )THEN
!
!        col= col_tmp
!
!      ENDIF
!
!
!      DO phi= 1, npart_eq, 1
!
!        long= phase_phi + phi*alpha
!
!        xtemp= center + rad_x*COS(long)*SIN(col)
!        ytemp= rad_y*SIN(long)*SIN(col)
!        ztemp= rad_z*COS(col)
!
!        IF( binary% import_mass_density( xtemp, ytemp, ztemp ) <= 0.0D0 &
!        )THEN
!
!          ghost_pos( 1, itr )= xtemp
!          ghost_pos( 2, itr )= ytemp
!          ghost_pos( 3, itr )= ztemp
!
!          itr= itr + 1
!
!        ENDIF
!
!      ENDDO
!    ENDDO
!  ENDDO
!  npart_ghost= itr - 1
!  IF( npart_ghost == 0 )THEN
!    PRINT *, "No ghost particles were placed. Stopping.."
!    PRINT *
!    STOP
!  ENDIF
!  ghost_pos = ghost_pos( :, 1:npart_ghost )

  !IF( run_apm )THEN

    PRINT *, "Placing ghost particles on a lattice..."

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

    PRINT *, "rad_x= ", rad_x
    PRINT *, "rad_y= ", rad_y
    PRINT *, "rad_z= ", rad_z

    itr= 0
    DO k= 1, nz, 1

      ztemp= zmin + dz/2 + ( k - 1 )*dz

      DO j= 1, ny, 1

        ytemp= ymin + dy/2 + ( j - 1 )*dy

        DO i= 1, nx, 1

          xtemp= xmin + dx/2 + ( i - 1 )*dx

          x_ell= center + rad_x*COS(ATAN( ytemp/xtemp )) &
                 *SIN(ACOS(ztemp/SQRT( ( xtemp - center )**2.0D0 + ytemp**2.0D0 + ztemp**2.0D0 )))

          y_ell= rad_y*SIN(ATAN( ytemp/xtemp )) &
                 *SIN(ACOS(ztemp/SQRT( ( xtemp - center )**2.0D0 + ytemp**2.0D0 + ztemp**2.0D0 )))

          z_ell= rad_z*( ztemp/SQRT( ( xtemp - center )**2.0D0 + ytemp**2.0D0 + ztemp**2.0D0 ))

          IF( binary% import_mass_density( xtemp, ytemp, ztemp ) <= 0.0D0 &
              .AND. &
              !SQRT( ( xtemp - center )**2.0D0 + ytemp**2.0D0 + ztemp**2.0D0 ) <= &
              !      larger_radius*( 1.0D0 + eps ) &
              SQRT( ( xtemp - center )**2.0D0 + ytemp**2.0D0 + ztemp**2.0D0 ) <= &
              1.1D0*SQRT( ( x_ell - center )**2.0D0 + delta*y_ell**2.0D0 + z_ell**2.0D0 ) &
              .AND. &
              !SQRT( ( xtemp - center )**2.0D0 + ytemp**2.0D0 + ztemp**2.0D0 ) >= &
              !      larger_radius + h_av/10.0D0 &
              SQRT( ( xtemp - center )**2.0D0 + ytemp**2.0D0 + ztemp**2.0D0 ) >= &
              SQRT( ( x_ell - center )**2.0D0 + delta*y_ell**2.0D0 + z_ell**2.0D0 ) &
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
      PRINT *, "No ghost particles were placed. Stopping.."
      PRINT *
      STOP
    ENDIF
    ghost_pos = ghost_pos( :, 1:npart_ghost )

    PRINT *, " * ", npart_ghost, " ghost particles placed."
    PRINT *, "npart_real=", npart_real
    PRINT *, "SIZE(ghost_pos)=", SIZE(ghost_pos)
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
        pos( 1, a ), &
        pos( 2, a ), &
        pos( 3, a )
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

    npart_all= npart_real + npart_ghost
    ALLOCATE( all_pos( 3, npart_all ) )

    all_pos( :, 1:npart_real )          = pos_star1
    all_pos( :, npart_real+1:npart_all )= ghost_pos

    PRINT *, "npart_real=", npart_real
    PRINT *, "npart_ghost=", npart_ghost
    PRINT *, "npart_all=", npart_all

  !ENDIF

  !-------------------------------!
  !-- Needed calls from SPHINCS --!
  !-------------------------------!

  ! setup unit system
  CALL set_units('NSM')
  CALL read_options

  npart= npart_all
  CALL allocate_SPH_memory

  !pos_u( :, 1:npart_real )= all_pos( :, 1:npart_real )
  !vel_u( :, 1:npart_real )= v  ( :, 1:npart_real )
  !h    ( 1:npart_real )   = h0 ( 1:npart_real )

  CALL allocate_RCB_tree_memory_3D(npart)
  iorig(1:npart)= (/ (a,a=1,npart) /)

  !CALL allocate_metric_on_particles(npart)
  ! for now: just constant
  !av(1:npart)= av_max

  ! tabulate kernel, get ndes
  CALL ktable(ikernel,ndes)

  ! flag that particles are 'alive'
  ALLOCATE( alive( npart ) )
  alive( 1:npart )= 1

  CALL allocate_gradient( npart )
  CALL allocate_metric_on_particles( npart )

  !PRINT *, " * Computing neighbours..."
  !PRINT *
  !
  !CALL exact_nei_tree_update( ndes,    &
  !                            npart,   &
  !                            all_pos, &
  !                            nu )

  !-----------------------------------------------------------------------!
  !-- At this point I can apply the artificial pressure method (can I?) --!
  !-----------------------------------------------------------------------!

  IF( run_apm )THEN

    PRINT *, " * Setting up ID for APM iteration..."
    PRINT *

    ! In setup_tov_star distribute_and_stretch is called
    ! In distribute_and_stretch, setup_uniform_sphere is called
    ! setup_uniform_sphere first place positions in a cubic lattice or a
    ! close-packed lattice, then assigns the smoothing length
    ! The ghost particles are included already here

    PRINT *, "assign h..."
    PRINT *
    CALL assign_h( nn_des, &
                   npart_all, &
                   all_pos, h0, &
                   h )
    PRINT *, "npart_real=", npart_real
    PRINT *, "npart_all=", npart_all
    PRINT *, "npart=", npart
    PRINT *

    !STOP
    PRINT *, "density loop..."
    PRINT *
    ! measure SPH-particle number density
    ALLOCATE( nstar_real( npart_all ) )
    nu= 1.0D0
    CALL density_loop( npart_all, all_pos, &    ! input
                       nu, h, nstar_real )  ! output

    ALLOCATE( nstar_p( npart_all ) )
    ! The following stands for get_profile_density
    nstar_p( 1:npart_real )      = nstar0( 1:npart_real )
    nstar_p( npart_real+1:npart_all )= 0.0D0

    ! In setup_uniform_sphere, get_profile_density is called
    ! this computed nstar, but we have it from the ID file

    ! assign nu's
    nu= nstar_p/nstar_real

    ! In setup_uniform_sphere, reset_COM is called
    !CALL reset_COM(npart,nu,pos_u)
  !  CALL COM( npart, all_pos, nu, &       ! input
  !            com_x, com_y, com_z, com_d) ! output
  !
  !  DO a= 1, npart, 1
  !     all_pos(1,a)= all_pos(1,a) - com_x
  !     all_pos(2,a)= all_pos(2,a) - com_y
  !     all_pos(3,a)= all_pos(3,a) - com_z
  !  ENDDO

    ! at this point, the particle should be mirrored (?)

    PRINT *, "density loop..."
    PRINT *
    ! Re-estimate nu
    CALL density_loop( npart_all, all_pos, &    ! input
                       nu, h, nstar_real )  ! output

    PRINT *, "import ID..."
    PRINT *
    ! Here you should get the LORENE density on the new positions
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

    max_nu= 0.0D0
    min_nu= 1.0D60

    ALLOCATE( art_pr( npart_all ) )

    DO a= 1, npart_real, 1

      ! Coordinate velocity of the fluid [c]
      vel_u(0,a)= 1.0D0
      vel_u(1,a)= lapse(a)*v_euler_x(a)- shift_x(a)
      vel_u(2,a)= lapse(a)*v_euler_y(a)- shift_y(a)
      vel_u(3,a)= lapse(a)*v_euler_z(a)- shift_z(a)

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
                   + g4(mus,nus)*vel_u(mus,a)*vel_u(nus,a)
        ENDDO
      ENDDO
      Theta_a= 1.0D0/SQRT(-Theta_a)
      Theta(a)= Theta_a

      nstar_p(a)= sq_g*Theta_a*baryon_density(a)*((Msun_geo*km2m)**3)/(amu*g2kg)

      ! Recompute nu taking into account the nstar just computed from th star
      ! Maybe this step s not needed, and the artificial pressure can be
      ! computed at this point
  !    nu_corr= 1.0D0 + (nstar_real(a) - nstar(a))/nstar_real(a)
  !    nu_corr= MAX(nu_corr,0.2D0)
  !    nu_corr= MIN(nu_corr,5.0D0)
  !    nu(a)  = nu_corr*nu(a)
  !
  !    ! baryon numbers
  !    max_nu= MAX(nu(a),max_nu)
  !    min_nu= MIN(nu(a),min_nu)

      ! inside
      !IF( NORM2(pos_u(:,a)) < Rstar )THEN
      !   dNstar=     (Nstar_real(a)-Nstar_P(a))/Nstar_P(a)
      !   aPr(a)=     MAX(1.0D0 + dNstar,0.1D0)
      !   aPr_max=    MAX(aPr_max,aPr(a))
      !   err_N_max=  MAX(err_N_max,ABS(dNstar))
      !   err_N_min=  MIN(err_N_min,ABS(dNstar))
      !   err_N_mean= err_N_mean + ABS(dNstar)
      !   n_inside=   n_inside + 1
      !   ! outside
      !ELSE
      !   IF( it == 1 )THEN
      !      aPr(a)= aPr_0
      !   ELSE
      !      aPr(a)= aPr_outside
      !   ENDIF
      !ENDIF

      !dNstar= ( nstar_real(a) - nstar_p(a) )/nstar_p(a)
      !art_pr(a) = MAX( 1.0D0 + dNstar, 0.1D0 )
      !art_pr_max= MAX( art_pr_max, art_pr(a) )
      !err_N_max = MAX( err_N_max, ABS(dNstar) )
      !err_N_min = MIN( err_N_min, ABS(dNstar) )
      !err_N_mean= err_N_mean + ABS(dNstar)

    ENDDO

    ! Here finishes setup_uniform_sphere

    PRINT *, " * ID set up for the APM iteration."
    PRINT *

    !-------------------------------------------------!
    !-- iterate to (close to) uniform particle mass --!
    !-------------------------------------------------!

    PRINT *, " * Performing APM iteration..."
    PRINT *

    ALLOCATE( freeze( npart_all ) )
    ALLOCATE( correction_pos( 3, npart_all ) )
    ALLOCATE( h_guess( npart_all ) )

    nu_all= (mass_star/DBLE(npart_real))*umass/amu
    nu= nu_all
    DO a= 1, npart_all
      IF( a < npart_real )THEN
        freeze(a)= 0
      ELSE
        freeze(a)= 1
      ENDIF
    ENDDO

    PRINT *, "iterating..."

    n_inc= 0
    err_N_mean_min= HUGE(1.0D0)
    apm_iteration: DO itr= 1, apm_max_it, 1

      PRINT *, "assign h..."

      h_guess= h
      CALL assign_h( nn_des, &
                     npart_all, &
                     all_pos, h_guess, &
                     h )

      PRINT *, "density_loop..."

      CALL density_loop( npart_all, all_pos, &    ! input
                         nu, h, nstar_real )      ! output

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
        vel_u(0,a)= 1.0D0
        vel_u(1,a)= lapse(a)*v_euler_x(a)- shift_x(a)
        vel_u(2,a)= lapse(a)*v_euler_y(a)- shift_y(a)
        vel_u(3,a)= lapse(a)*v_euler_z(a)- shift_z(a)

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
                     + g4(mus,nus)*vel_u(mus,a)*vel_u(nus,a)
          ENDDO
        ENDDO
        Theta_a= 1.0D0/SQRT(-Theta_a)
        Theta(a)= Theta_a

        nstar_p(a)= sq_g*Theta_a*baryon_density(a)*((Msun_geo*km2m)**3)/(amu*g2kg)

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

        IF( ISNAN(nstar_real(a)) )THEN
          PRINT *, "nstar_real is a NaN at particle ", a
          PRINT *, "nstar_p is ", nstar_p(a)
          STOP
        ENDIF
        IF( ISNAN(nstar_p(a)) )THEN
          PRINT *, "nstar_p is a NaN at particle ", a
          PRINT *, "nstar_real is ", nstar_real(a)
          STOP
        ENDIF
        IF( ISNAN(dNstar) )THEN
          PRINT *, "dNstar is a NaN at particle ", a
          PRINT *, "nstar_real is ", nstar_real(a)
          PRINT *, "nstar_p is ", nstar_p(a)
          STOP
        ENDIF

      ENDDO

      nstar_p( npart_real+1:npart_all )= 0.0D0
      art_pr ( npart_real+1:npart_all )= art_pr_max

      PRINT *, "Before calling position_correction"

      !PRINT *, npart_all
      !PRINT *, SIZE(all_pos(1,:))
      !PRINT *, SIZE(h)
      !PRINT *, SIZE(art_pr)
      !PRINT *, SIZE(nstar_real)
      !PRINT *, SIZE(correction_pos(1,:))

      CALL position_correction( npart_all, &
                                all_pos, h, nu_all, art_pr, nstar_real, &
                                correction_pos )

      PRINT *, "After calling position_correction"

      ! TODO: mirror the particles here (or only afterwards?)
      itr2= 0
      DO a= 1, npart_real, 1
        pos_tmp= all_pos(:,a) + correction_pos(:,a)
        IF( &!binary% import_mass_density( &
            !                    all_pos(1,a), all_pos(2,a), all_pos(3,a) ) > 0 &
            !.AND. &
            binary% import_mass_density( &
                                pos_tmp(1), pos_tmp(2), pos_tmp(3) ) > 0.0D0 &
            .AND. &
            nstar_p(a) > 0.0D0 &
            .AND. &
            binary% is_hydro_negative( &
                                pos_tmp(1), pos_tmp(2), pos_tmp(3) ) == 0 &
        )THEN
          itr2= itr2 + 1
          !IF( binary% is_hydro_negative( &
          !                        all_pos(1,a), all_pos(2,a), all_pos(3,a) ) > 0 &
          !    .OR. &
          !    binary% is_hydro_negative( &
          !                        pos_tmp(1), pos_tmp(2), pos_tmp(3) ) > 0 &
          !)THEN
            all_pos(:,a)= pos_tmp
          !ELSE
          !  all_pos(:,a)= all_pos(:,a) + correction_pos(:,a)/2.0D0
          !ENDIF
        ENDIF
      ENDDO

      !PRINT *, "err_N_mean= ", err_N_mean
      !npart_real= itr2
      err_N_mean= err_N_mean/DBLE(npart_real)

      err_N_mean_min= MIN( err_N_mean, err_N_mean_min )

      !PRINT *, "itr2= ", itr2
      PRINT *, "npart_real= ", npart_real
      PRINT *
      PRINT *, '...done with position update #: ', itr
      PRINT *
      PRINT *, '...err_N_max=  ', err_N_max
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
      PRINT *, '...err_N_min=  ', err_N_min
      PRINT *, '...err_N_mean= ', err_N_mean
      PRINT *, '...err_N_mean_min= ', err_N_mean_min
      PRINT *

      PRINT *, "max_nu=", MAXVAL( nstar_p(1:npart_real)/nstar_real(1:npart_real), DIM= 1 )
      PRINT *, "min_nu=", MINVAL( nstar_p(1:npart_real)/nstar_real(1:npart_real), DIM= 1 )
      PRINT *, "max_nu/min_nu=", MAXVAL( nstar_p(1:npart_real)/nstar_real(1:npart_real), DIM= 1 )/MINVAL( nstar_p(1:npart_real)/nstar_real(1:npart_real), DIM= 1 )
      PRINT *

      ! exit condition
      !IF( err_N_mean > err_mean_old ) n_inc= n_inc + 1
      !IF( ABS( err_N_mean - err_mean_old )/ABS( err_mean_old ) < iter_tol &
      !    .AND. &
      !    err_N_max < 10.0D0 &
      !)THEN
      !  n_inc= n_inc + 1
      !  PRINT *, "n_inc/max_inc= ", n_inc, "/", max_inc
      !  PRINT *, "ABS( err_N_mean - err_mean_old )/ABS(err_mean_old)= ", &
      !           ABS( err_N_mean - err_mean_old )/ABS(err_mean_old)
      !ENDIF
      IF( ABS(err_N_mean - err_N_mean_min)/ABS(err_N_mean_min) < iter_tol &
      )THEN
        n_inc= n_inc + 1
        PRINT *, "n_inc/max_inc= ", n_inc, "/", max_inc
        PRINT *, "ABS(err_N_mean - err_N_mean_min)/ABS(err_N_mean_min)= ", &
                 ABS(err_N_mean - err_N_mean_min)/ABS(err_N_mean_min), " < ", &
                 iter_tol
      ENDIF
      IF( n_inc == max_inc ) EXIT
      err_mean_old= err_N_mean

    ENDDO apm_iteration

    PRINT *, "Iteration completed."
    PRINT *

    finalnamefile= "apm_pos.dat"

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

    ! Now get rid of the ghost particles
    pos= all_pos( :, 1:npart_real )
    npart= npart_real
    PRINT *, npart

    h_guess= h

  ELSE

    nu_all= (mass_star/DBLE(npart_real))*umass/amu

    finalnamefile= "apm_pos.dat"

    INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

    IF( exist )THEN
      OPEN( UNIT= unit_id, FILE= TRIM(finalnamefile), &
            FORM= "FORMATTED", ACTION= "READ", IOSTAT= ios, &
            IOMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...error when opening " // TRIM(finalnamefile), &
                ". The error message is", err_msg
        STOP
      ENDIF
    ELSE
      PRINT *, "** ERROR! Unable to find file " // TRIM(finalnamefile)
      STOP
    ENDIF

    PRINT *, " * Reading file " // TRIM(finalnamefile) // "..."
    ! Reading particle data
    DO itr= 1, npart_real, 1
      READ( UNIT= unit_id, FMT= *, IOSTAT = ios, IOMSG= err_msg ) &
        tmp, tmp2, &
        pos( 1, itr ), &
        pos( 2, itr ), &
        pos( 3, itr )

      IF( ios > 0 )THEN
        PRINT *, "...error when reading " // TRIM(finalnamefile), &
                " at step ", itr,". The status variable is ", ios, &
                ". The error message is", err_msg
        STOP
      ENDIF
    ENDDO
    IF( tmp2 /= npart_real )THEN
      PRINT *, "** ERROR! Mismatch in the values of the particle number read", &
               " from " // TRIM(finalnamefile), " in two ways. "
      PRINT *, "tmp= ", tmp, ", npart_real= ", npart_real
      PRINT *, "Stopping..."
      STOP
    ENDIF

    PRINT *, " * Positions of real particles read."

    DO itr= npart_real + 1, npart_all, 1
      READ( UNIT= unit_id, FMT= *, IOSTAT = ios, IOMSG= err_msg ) &
        tmp, tmp2, &
        ghost_pos( 1, itr - npart_real ), &
        ghost_pos( 2, itr - npart_real ), &
        ghost_pos( 3, itr - npart_real )

      IF( ios > 0 )THEN
        PRINT *, "...error when reading " // TRIM(finalnamefile), &
                " at step ", itr,". The status variable is ", ios, &
                ". The error message is", err_msg
        STOP
      ENDIF
    ENDDO

    CLOSE( unit_id )

    PRINT *, " * Positions of ghost particles read."

    ALLOCATE( h_guess( npart_real ) )

    h_guess(1:npart_real)= h0
    !h_guess(npart_real+1:npart_all)= h0(1)

    ALLOCATE( nstar_real( npart_real ) )
    ALLOCATE( nstar_p( npart_real ) )

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

  ENDIF

  !----------------------------------------------------!
  !-- at this point, pos contains the real particles --!
  !----------------------------------------------------!

  !----------------------------!
  !-- enforce centre of mass --!
  !----------------------------!

  ! Rescale to the center of mass
  ! TODO: this rescales the COM to (0,0,0), but this is wrong for the binary
  !       commented for now
!  CALL COM( npart, pos, nu, &       ! input
!            com_x, com_y, com_z, com_d) ! output
!
!  DO a= 1, npart, 1
!     all_pos(1,a)= all_pos(1,a) - com_x
!     all_pos(2,a)= all_pos(2,a) - com_y
!     all_pos(3,a)= all_pos(3,a) - com_z
!  ENDDO

  !-------------------------------------------------------------------!
  !-- now assign baryon number to match profile as good as possible --!
  !-------------------------------------------------------------------!

  PRINT *, " * Assign baryon number..."
  PRINT *

  PRINT *, "1"

  h      = h(1:npart_real)
  h_guess= h_guess(1:npart_real)
  nu     = nu(1:npart_real)

  CALL assign_h( nn_des, &           !
                 npart_real, &        !
                 pos, h_guess, & ! Input
                 h )                 ! Output

  PRINT *, "2"

  ! Measure SPH particle number density
  nu= 1.0D0
  CALL density_loop( npart_real, pos, &    ! input
                     nu, h, nstar_real )      ! output

  PRINT *, "3"

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

  PRINT *, "3.5"

  DO a= 1, npart_real, 1

    IF( baryon_density(a) <= 0 )THEN
      PRINT *, " * The baryon density is 0 at particle ", a
      PRINT *
      STOP
    ENDIF

    ! Coordinate velocity of the fluid [c]
    vel_u(0,a)= 1.0D0
    vel_u(1,a)= lapse(a)*v_euler_x(a)- shift_x(a)
    vel_u(2,a)= lapse(a)*v_euler_y(a)- shift_y(a)
    vel_u(3,a)= lapse(a)*v_euler_z(a)- shift_z(a)

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
                 + g4(mus,nus)*vel_u(mus,a)*vel_u(nus,a)
      ENDDO
    ENDDO
    Theta_a= 1.0D0/SQRT(-Theta_a)
    Theta(a)= Theta_a

    nstar_p(a)= sq_g*Theta_a*baryon_density(a)*((Msun_geo*km2m)**3)/(amu*g2kg)

  ENDDO

  !nstar_p( npart_real+1:npart_all )= 0.0D0

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

  PRINT *, "4"

  PRINT *, "npart_real= ", npart_real
  PRINT *, "SIZE(nu)= ", SIZE(nu)
  PRINT *

  nu_ratio= MAXVAL( nu, DIM= 1 )/MINVAL( nu, DIM= 1 )
  PRINT *, "nu_ratio before correction = ", nu_ratio
  PRINT *

  !-----------------------------------------!
  !-- monitoring before correction  nu... --!
  !-----------------------------------------!

 ! ! measure density
 ! CALL density_loop( npart_real, pos, &    ! input
 !                    nu, h, nstar_real )      ! output
 !
 ! CALL binary% import_id( npart_real, pos(1,:), &
 !                         pos(2,:), &
 !                         pos(3,:), &
 !                         lapse, shift_x, shift_y, shift_z, &
 !                         g_xx, g_xy, g_xz, &
 !                         g_yy, g_yz, g_zz, &
 !                         baryon_density, &
 !                         energy_density, &
 !                         specific_energy, &
 !                         pressure, &
 !                         v_euler_x, v_euler_y, v_euler_z )
 !
 ! DO a= 1, npart_real, 1
 !
 !   ! Coordinate velocity of the fluid [c]
 !   vel_u(0,a)= 1.0D0
 !   vel_u(1,a)= lapse(a)*v_euler_x(a)- shift_x(a)
 !   vel_u(2,a)= lapse(a)*v_euler_y(a)- shift_y(a)
 !   vel_u(3,a)= lapse(a)*v_euler_z(a)- shift_z(a)
 !
 !   !
 !   !-- Metric as matrix for easy manipulation
 !   !
 !   g4(0,0)= - lapse(a)**2 + g_xx(a)*shift_x(a)*shift_x(a)&
 !          + 2*g_xy(a)*shift_x(a)*shift_y(a) &
 !          + 2*g_xz(a)*shift_x(a)*shift_z(a) &
 !          + g_yy(a)*shift_y(a)*shift_y(a) &
 !          + 2*g_yz(a)*shift_y(a)*shift_z(a) &
 !          + g_zz(a)*shift_z(a)*shift_z(a)
 !   g4(0,1)= g_xx(a)*shift_x(a) + g_xy(a)*shift_y(a) + g_xz(a)*shift_z(a)
 !   g4(0,2)= g_xy(a)*shift_x(a) + g_yy(a)*shift_y(a) + g_yz(a)*shift_z(a)
 !   g4(0,3)= g_xz(a)*shift_x(a) + g_yz(a)*shift_y(a) + g_zz(a)*shift_z(a)
 !
 !   g4(1,0)= g_xx(a)*shift_x(a) + g_xy(a)*shift_y(a) + g_xz(a)*shift_z(a)
 !   g4(1,1)= g_xx(a)
 !   g4(1,2)= g_xy(a)
 !   g4(1,3)= g_xz(a)
 !
 !   g4(2,0)= g_xy(a)*shift_x(a) + g_yy(a)*shift_y(a) + g_yz(a)*shift_z(a)
 !   g4(2,1)= g_xy(a)
 !   g4(2,2)= g_yy(a)
 !   g4(2,3)= g_yz(a)
 !
 !   g4(3,0)= g_xz(a)*shift_x(a) + g_yz(a)*shift_y(a) + g_zz(a)*shift_z(a)
 !   g4(3,1)= g_xz(a)
 !   g4(3,2)= g_yz(a)
 !   g4(3,3)= g_zz(a)
 !
 !   ! sqrt(-det(g4))
 !   CALL determinant_4x4_matrix(g4,det)
 !   IF( ABS(det) < 1D-10 )THEN
 !       PRINT *, "The determinant of the spacetime metric is " &
 !                // "effectively 0 at particle ", a
 !       STOP
 !   ELSEIF( det > 0 )THEN
 !       PRINT *, "The determinant of the spacetime metric is " &
 !                // "positive at particle ", a
 !       STOP
 !   ENDIF
 !   sq_g= SQRT(-det)
 !
 !   !
 !   !-- Generalized Lorentz factor
 !   !
 !   Theta_a= 0.D0
 !   DO nus=0,3
 !     DO mus=0,3
 !       Theta_a= Theta_a &
 !                + g4(mus,nus)*vel_u(mus,a)*vel_u(nus,a)
 !     ENDDO
 !   ENDDO
 !   Theta_a= 1.0D0/SQRT(-Theta_a)
 !   Theta(a)= Theta_a
 !
 !   nstar_p(a)= sq_g*Theta_a*baryon_density(a)*((Msun_geo*km2m)**3)/(amu*g2kg)
 !
 ! ENDDO
 !
 ! !nstar_p( npart_real+1:npart_all )= 0.0D0
 !
 ! ! get RELATIVE nu's right
 ! dN_av= 0.0D0
 !     dN_max= 0.0D0
 ! DO a= 1, npart_real, 1
 !   dN=     ABS(nstar_real(a)-nstar_p(a))/nstar_p(a)
 !   dN_max= MAX(dN_max,dN)
 !   dN_av=  dN_av + dN
 ! ENDDO
 ! dN_av= dN_av/DBLE(npart_real)
 ! PRINT*,'...dN_max ',dN_max
 ! PRINT*,'...dN_av  ',dN_av
 !
 ! IF( .NOT.ALLOCATED( nstar_int ) ) ALLOCATE( nstar_int( npart_real ) )
 !
 ! CALL exact_nei_tree_update( nn_des, &           !
 !                             npart_real, &        !
 !                             pos, nu )
 !
 ! CALL density( npart_real, pos, nstar_int )
 !
 ! PRINT *, "0"
 !
 ! finalnamefile= "densities-before-nu-correction.dat"
 !
 ! INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )
 !
 ! IF( exist )THEN
 !     OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
 !           FORM= "FORMATTED", &
 !           POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
 !           IOMSG= err_msg )
 ! ELSE
 !     OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "NEW", &
 !           FORM= "FORMATTED", &
 !           ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
 ! ENDIF
 ! IF( ios > 0 )THEN
 !   PRINT *, "...error when opening " // TRIM(finalnamefile), &
 !            ". The error message is", err_msg
 !   STOP
 ! ENDIF
 !
 ! PRINT *, "1"
 !
 ! IF( .NOT.ALLOCATED( nu_one ) ) ALLOCATE( nu_one( npart_real ) )
 ! IF( .NOT.ALLOCATED( particle_density_final ) ) &
 !   ALLOCATE( particle_density_final( npart_real ) )
 ! nu_one= 1.0D0
 ! CALL density_loop( npart_real, pos, &    ! input
 !                    nu_one, h, particle_density_final )      ! output
 !
 ! DO a= 1, npart_real, 1
 !   WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
 !     a, &
 !     pos( 1, a ), pos( 2, a ), pos( 3, a ), &
 !     nstar_p( a ), &
 !     nstar_int( a ), &
 !     particle_density_final( a ), &
 !     particle_density_final( a )*nstar_p( 1 )/particle_density_final( 1 ), &
 !     ABS(nstar_real(a)-nstar_p(a))/nstar_p(a), &
 !     nu(a)
 ! ENDDO
 !
 ! CLOSE( UNIT= 2 )

  nu= nu_all

  !----------------------!
  !-- Correcting nu... --!
  !----------------------!

  PRINT *, "Correcting nu"
  PRINT *

  !nu(1:npart_real)= nstar_p(1:npart_real)/nstar_real(1:npart_real)
  !nu= nu(1:npart_real)

  !nu(1:npart_all)= nstar_p(1:npart_all)/nstar_real(1:npart_all)
  !nu= nstar_p/nstar_real
  !
  !DO a= 1, npart_real, 1
  !
  !  nu_tmp= nu(a)
  !  nu(a)= nstar_p(a)/nstar_real(a)
  !  IF( MAXVAL( nu, DIM= 1 )/MINVAL( nu, DIM= 1 ) > 3.0D0*nu_ratio )THEN
  !    nu(a)= nu_tmp
  !  ENDIF
  !
  !ENDDO

  DO a= 1, npart_real, 1

    nu_tmp= nu(a)
    nu(a)= nstar_p(a)/nstar_real(a)
    IF( MAXVAL( nu, DIM= 1 )/MINVAL( nu, DIM= 1 ) > 3.0D0*nu_ratio )THEN
      nu(a)= nu_tmp
    ENDIF

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

  !STOP

  IF( post_correction )THEN

    ! just a few iterations to NOT get the nu-ratio too large
    mass_it: DO itr= 1, m_max_it, 1

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
         vel_u(0,a)= 1.0D0
         vel_u(1,a)= lapse(a)*v_euler_x(a)- shift_x(a)
         vel_u(2,a)= lapse(a)*v_euler_y(a)- shift_y(a)
         vel_u(3,a)= lapse(a)*v_euler_z(a)- shift_z(a)

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
                      + g4(mus,nus)*vel_u(mus,a)*vel_u(nus,a)
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

    ENDDO mass_it

  ENDIF


  ! TODO: Enforce center of mass
  !       Mirror particles

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

  PRINT *, "nu_tot=", nu_tot
  PRINT *, "mass estimate= ", nu_tot*amu/MSun, "=", &
           100.0D0*nu_tot*amu/MSun/mass, "% of the LORENE baryon mass"
  PRINT *
  PRINT *, "mean_nu=", mean_nu
  PRINT *, "variance_nu=", variance_nu
  PRINT *, "stddev_nu=", stddev_nu
  PRINT *, "stddev_nu/mean_nu=", stddev_nu/mean_nu
  PRINT *

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
    vel_u(0,a)= 1.0D0
    vel_u(1,a)= lapse(a)*v_euler_x(a)- shift_x(a)
    vel_u(2,a)= lapse(a)*v_euler_y(a)- shift_y(a)
    vel_u(3,a)= lapse(a)*v_euler_z(a)- shift_z(a)

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
                 + g4(mus,nus)*vel_u(mus,a)*vel_u(nus,a)
      ENDDO
    ENDDO
    Theta_a= 1.0D0/SQRT(-Theta_a)
    Theta(a)= Theta_a

    nstar_p(a)= sq_g*Theta_a*baryon_density(a)*((Msun_geo*km2m)**3)/(amu*g2kg)

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
  PRINT*,'...dN_max ',dN_max
  PRINT*,'...dN_av  ',dN_av

  IF( .NOT.ALLOCATED( nstar_int ) ) ALLOCATE( nstar_int( npart_real ) )

  CALL exact_nei_tree_update( nn_des, &           !
                              npart_real, &        !
                              pos, nu )

  CALL density( npart_real, pos, nstar_int )

PRINT *, "0"

  finalnamefile= "densities.dat"

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

PRINT *, "1"

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
      nu(a)
  ENDDO

  CLOSE( UNIT= 2 )

  ! Here setup_uniform_sphere ends

  PRINT *, "** End of PROGRAM proto_apm."
  PRINT *

  CONTAINS


END PROGRAM proto_apm
