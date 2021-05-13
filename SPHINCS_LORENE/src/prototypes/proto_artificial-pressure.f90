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

  USE bns_id,         ONLY: bns
  !USE particles_id,   ONLY: particles
  USE RCB_tree_3D,    ONLY: allocate_RCB_tree_memory_3D,&
                            iorig
  USE kernel_table,   ONLY: ktable
  USE sph_variables,  ONLY: Rstar,divv,av,Pr,ye,temp,nlrf,&
                            u,Theta,vel_u,tterm,tgrav,tkin,&
                            escap,t,n1,n2,pos_u,h,nu,npart,&
                            npm,Nstar,allocate_sph_memory
  USE options,        ONLY: basename,ikernel,ndes,av_max,icow
  USE input_output,   ONLY: read_options
  USE units,          ONLY: umass,m0c2_CU,set_units

  USE APM,            ONLY: density_loop, position_correction, assign_h
  USE set_h,          ONLY: exact_nei_tree_update
  USE analyze,        ONLY: COM
  USE matrix,         ONLY: determinant_4x4_matrix
  USE constants,      ONLY: MSun, Msun_geo, km2m, g2kg, amu, pi, half
  USE NR,             ONLY: indexx

  IMPLICIT NONE

  INTEGER, PARAMETER:: unit_id  = 23
  INTEGER, PARAMETER:: max_npart= 5D+6
  INTEGER:: npart_tmp, tmp, ios, itr, a, nout, nus, mus, npart_eq, &
            npart_ghost_shell, npart_ghost, npart_all, r, th, phi, npart1
  INTEGER, DIMENSION(:), ALLOCATABLE:: x_sort, xy_sort, xyz_sort, lim

  DOUBLE PRECISION:: com_x, com_y, com_z, com_d, det, sq_g, Theta_a, &
                     nu_corr, max_nu, min_nu, smaller_radius, larger_radius, &
                     dr, rad, col, col_tmp, alpha, phase_th, rand_num, &
                     center, long, xtemp, ytemp, ztemp, radius_z, rel_sign, &
                     phase_phi, h_max, h_av, dNstar, art_pr_max, err_N_max , &
                     err_N_min, err_N_mean, nu_all, mass_star
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

  finalnamefile= "lorene-bns-id-particles-shells-2.dat"

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
  npart_tmp= 0
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
    npart_tmp= npart_tmp + 1
  ENDDO
  IF( tmp /= npart_tmp )THEN
    PRINT *, "** ERROR! Mismatch in the values of the particle number read", &
             " from " // TRIM(finalnamefile), " in two ways. "
    PRINT *, "tmp= ", tmp, ", npart_tmp= ", npart_tmp
    PRINT *, "Stopping..."
    STOP
  ENDIF

  CLOSE( unit_id )

  pos                  = pos                  ( :, 1:npart_tmp )
  lapse_parts          = lapse_parts          ( 1:npart_tmp )
  shift_parts_x        = shift_parts_x        ( 1:npart_tmp )
  shift_parts_y        = shift_parts_y        ( 1:npart_tmp )
  shift_parts_z        = shift_parts_z        ( 1:npart_tmp )
  baryon_density_parts = baryon_density_parts ( 1:npart_tmp )
  energy_density_parts = energy_density_parts ( 1:npart_tmp )
  specific_energy_parts= specific_energy_parts( 1:npart_tmp )
  pressure_parts       = pressure_parts       ( 1:npart_tmp )
  v_euler_parts_x      = v_euler_parts_x      ( 1:npart_tmp )
  v_euler_parts_y      = v_euler_parts_y      ( 1:npart_tmp )
  v_euler_parts_z      = v_euler_parts_z      ( 1:npart_tmp )
  v                    = v                    ( :, 1:npart_tmp )
  nu0                  = nu0                  ( 1:npart_tmp )
  nlrf0                = nlrf0                ( 1:npart_tmp )
  Ye0                  = Ye0                  ( 1:npart_tmp )
  Theta0               = Theta0               ( 1:npart_tmp )
  sph_density          = sph_density          ( 1:npart_tmp )
  nstar0               = nstar0               ( 1:npart_tmp )
  h0                   = h0                   ( 1:npart_tmp )

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
  !ALLOCATE( sorted_pos(3,npart_tmp) )
  !ALLOCATE( x_sort(npart_tmp) )
  !ALLOCATE( xy_sort(npart_tmp) )
  !ALLOCATE( xyz_sort(npart_tmp) )
  !sorted_pos= 0.0D0
  !x_sort    = 0.0D0
  !xy_sort   = 0.0D0
  !xyz_sort  = 0.0D0
  !
  !PRINT *, sorted_pos(:,1)
  !PRINT *, sorted_pos(:,npart_tmp)
  !PRINT *
  !
  !! Sort particle positions along x
  !CALL indexx( npart_tmp, pos( 1, : ), x_sort )
  !
  !DO a= 1, npart_tmp, 1
  !  sorted_pos(1,a)= pos(1,x_sort(a))
  !  sorted_pos(2,a)= pos(2,x_sort(a))
  !  sorted_pos(3,a)= pos(3,x_sort(a))
  !ENDDO
  !PRINT *, "After sorting x"
  !
  !! Why do you need to subsort the posiions over y and z? x is perhaps enough
  !
  !!PRINT *, sorted_pos(:,1)
  !!PRINT *, sorted_pos(:,npart_tmp)
  !!PRINT *
  !DO a= 1, npart_tmp, 1
  !  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
  !    a, &
  !    sorted_pos( 1, a ), &
  !    sorted_pos( 2, a ), &
  !    sorted_pos( 3, a )
  !ENDDO

  ! Sort along y
  !ALLOCATE( lim(npart_tmp) )
  !lim(0)= 0
  !itr= 1
  !DO a= 1, npart_tmp - 1, 1
  !  IF( sorted_pos(1,a) /= sorted_pos(1,a + 1) )THEN
  !    lim(itr)= a
  !    CALL indexx( lim(itr)-lim(itr-1), pos( 2, lim(itr-1)+1:lim(itr) ), &
  !                 xy_sort(lim(itr-1)+1:lim(itr)) )
  !    xy_sort(lim(itr-1)+1:lim(itr))= xy_sort(lim(itr-1)+1:lim(itr)) + lim(itr-1)
  !    itr= itr + 1
  !  ENDIF
  !ENDDO
  !!PRINT *, x_sort
  !DO a= 1, npart_tmp, 1
  !  sorted_pos(1,a)= sorted_pos(1,xy_sort(a))
  !  sorted_pos(2,a)= sorted_pos(2,xy_sort(a))
  !  sorted_pos(3,a)= sorted_pos(3,xy_sort(a))
  !ENDDO
  !
  !DO a= 1, npart_tmp, 1
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
  !PRINT *, sorted_pos(:,npart_tmp)
  !PRINT *

  ! Sort along z
  !lim(0)= 0
  !itr= 1
  !DO a= 1, npart_tmp - 1, 1
  !  IF( sorted_pos(2,a) /= sorted_pos(2,a + 1) )THEN
  !    lim(itr)= a
  !    CALL indexx( lim(itr)-lim(itr-1), pos( 3, lim(itr-1)+1:lim(itr) ), &
  !                 xyz_sort(lim(itr-1)+1:lim(itr)) )
  !    xyz_sort(lim(itr-1)+1:lim(itr))= xyz_sort(lim(itr-1)+1:lim(itr)) + lim(itr-1)
  !    itr= itr + 1
  !  ENDIF
  !ENDDO
  !DO a= 1, npart_tmp, 1
  !  sorted_pos(1,a)= pos(1,xyz_sort(a))
  !  sorted_pos(2,a)= pos(2,xyz_sort(a))
  !  sorted_pos(3,a)= pos(3,xyz_sort(a))
  !ENDDO
  !PRINT *, "After sorting z"
  !PRINT *, sorted_pos(:,1)
  !PRINT *, sorted_pos(:,npart_tmp)
  !PRINT *
  !
  !PRINT *, "End of sorting algorithm"
  !PRINT *

  !CLOSE( UNIT= 2 )

  ALLOCATE( ghost_pos( 3, max_npart ) )
  ALLOCATE( pos_star1( 3, npart_tmp ) )
  ghost_pos= 0.0D0
  pos_star1= 0.0D0
  itr= 1
  DO a= 1, npart_tmp, 1
    IF( pos( 1, a ) < 0.0D0 )THEN
      pos_star1( 1, itr )= pos( 1 ,a )
      pos_star1( 2, itr )= pos( 2 ,a )
      pos_star1( 3, itr )= pos( 3 ,a )
      itr= itr + 1
    ENDIF
  ENDDO
  npart1= itr
  pos_star1 = pos_star1( :, 1:npart1 - 1 )
  !smaller_radius= binary% get_radius1_x_opp()
  !larger_radius = binary% get_radius1_x_comp()
  center        = binary% get_center1_x()
  smaller_radius= ABS( MINVAL( pos_star1( 1, : ), DIM= 1 ) - center )
  larger_radius = ABS( center - MAXVAL( pos_star1( 1, : ), DIM= 1 ) )
  radius_z= ABS( MAXVAL( pos_star1( 3, : ), DIM= 1 ) )
  dr= ( larger_radius - smaller_radius )/5.0D0
  npart_eq= 50
  alpha= 2*pi/DBLE(npart_eq)
  !npart_ghost_shell= ( npart_eq**2 )/2
  itr= 1

  !PRINT *, MAXVAL( pos_star1( 1, : ), DIM= 1 ), &
  !         smaller_radius, larger_radius, center, dr, alpha

  ! You need to place ghost particles one smoothing length outside of the
  ! surface

  ! Find the maximum of the smoothing lencth of the particles whose distance form thecenter is higher than radius_z, and place particles outside of that

  h_max= 0.0D0
  h_av = 0.0D0
  itr  = 0
  DO a= 1, npart_tmp, 1
    IF( SQRT( pos( 1, a )**2.0D0 + pos( 2, a )**2.0D0 + pos( 3, a )**2.0D0 ) &
        > 0.99D0*radius_z )THEN
      itr= itr + 1
      IF( h0(a) > h_max )THEN
        h_max= h0(a)
      ENDIF
      h_av= h_av + h0(a)
    ENDIF
  ENDDO
  h_av= h_av/itr

  PRINT *, " * Placing ghost particles..."

  DO r= 1, 10, 1

    rad= radius_z + h_av + DBLE( r )*dr

    DO th= 1, npart_eq/2, 1

      CALL RANDOM_NUMBER( phase_phi )
      phase_phi= phase_phi*alpha

      col= ( th - 1 )*alpha + alpha/2.0D0
      CALL RANDOM_NUMBER( phase_th )
      CALL RANDOM_NUMBER( rand_num )
      IF( rand_num >= half ) rel_sign=  1
      IF( rand_num < half )  rel_sign= -1

      col_tmp= col*( 1.0D0 + rel_sign*0.05D0*phase_th )

      IF( col_tmp < pi .AND. col_tmp > 0 )THEN

        col= col_tmp

      ENDIF


      DO phi= 1, npart_eq, 1

        long= phase_phi + phi*alpha

        xtemp= center + rad*COS(long)*SIN(col)
        ytemp= rad*SIN(long)*SIN(col)
        ztemp= rad*COS(col)

        IF( binary% import_mass_density( xtemp, ytemp, ztemp ) <= 0.0D0 &
        )THEN

          ghost_pos( 1, itr )= xtemp
          ghost_pos( 2, itr )= ytemp
          ghost_pos( 3, itr )= ztemp

          itr= itr + 1

        ENDIF

      ENDDO
    ENDDO
  ENDDO
  npart_ghost= itr
  ghost_pos = ghost_pos( :, 1:npart_ghost )

  PRINT *, " * ", npart_ghost, " ghost particles placed."
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

  DO a= 1, npart_ghost, 1
    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      a, &
      ghost_pos( 1, a ), &
      ghost_pos( 2, a ), &
      ghost_pos( 3, a )
  ENDDO

  CLOSE( UNIT= 2 )

  PRINT *, " * Printed."

  npart_all= npart + npart_ghost
  ALLOCATE( all_pos( 3, npart_all ) )

  all_pos( :, 1:npart_tmp )          = pos
  all_pos( :, npart_tmp+1:npart_all )= ghost_pos

  !-------------------------------!
  !-- Needed calls from SPHINCS --!
  !-------------------------------!

  ! setup unit system
  CALL set_units('NSM')

  !-------------------------------!
  !-- read main options for run --!
  !-------------------------------!
  CALL read_options

  npart= npart_all
  CALL allocate_SPH_memory

  pos_u( :, 1:npart_tmp )= all_pos( :, 1:npart_tmp )
  vel_u( :, 1:npart_tmp )= v  ( :, 1:npart_tmp )
  h    ( 1:npart_tmp )   = h0 ( 1:npart_tmp )

  CALL allocate_RCB_tree_memory_3D(npart)
  iorig(1:npart)= (/ (a,a=1,npart) /)

  !CALL allocate_metric_on_particles(npart)
  ! for now: just constant
  !av(1:npart)= av_max

  ! tabulate kernel, get ndes
  CALL ktable(ikernel,ndes)

  !-----------------------------------------------------------------------!
  !-- At this point I can apply the artificial pressure method (can I?) --!
  !-----------------------------------------------------------------------!

  PRINT *, " * Setting up ID for APM iteration..."
  PRINT *

  ! In setup_tov_star distribute_and_stretch is called
  ! In distribute_and_stretch, setup_uniform_sphere is called
  ! setup_uniform_sphere first place positions in a cubic lattice or a
  ! close-packed lattice, then assigns the smoothing length
  ! The ghost particles are included already here

  PRINT *, "assign h..."
  PRINT *
  CALL assign_h( npart_all, &
                 npart, &
                 all_pos, h0, &
                 h )
  PRINT *, "npart_tmp=", npart_tmp
  PRINT *, "npart_all=", npart_all
  PRINT *, "npart=", npart
  PRINT *

  !STOP
  PRINT *, "density loop..."
  PRINT *
  ! measure SPH-particle number density
  ALLOCATE( Nstar_real(npart) )
  nu= 1.0D0
  CALL density_loop( npart, all_pos, &    ! input
                     nu, h, nstar_real )  ! output

  ALLOCATE( nstar_p( npart ) )
  ! The following stands for get_profile_density
  nstar_p( 1:npart_tmp )      = nstar0( 1:npart_tmp )
  nstar_p( npart_tmp+1:npart )= 0.0D0

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
  CALL density_loop( npart, all_pos, &    ! input
                     nu, h, nstar_real )  ! output

  PRINT *, "import ID..."
  PRINT *
  ! Here you should get the LORENE density on the new positions
  ALLOCATE( lapse          (npart) )
  ALLOCATE( shift_x        (npart) )
  ALLOCATE( shift_y        (npart) )
  ALLOCATE( shift_z        (npart) )
  ALLOCATE( g_xx           (npart) )
  ALLOCATE( g_xy           (npart) )
  ALLOCATE( g_xz           (npart) )
  ALLOCATE( g_yy           (npart) )
  ALLOCATE( g_yz           (npart) )
  ALLOCATE( g_zz           (npart) )
  ALLOCATE( baryon_density (npart) )
  ALLOCATE( energy_density (npart) )
  ALLOCATE( specific_energy(npart) )
  ALLOCATE( pressure       (npart) )
  ALLOCATE( v_euler_x      (npart) )
  ALLOCATE( v_euler_y      (npart) )
  ALLOCATE( v_euler_z      (npart) )

  CALL binary% import_id( npart, all_pos(1,1:npart_tmp), &
                          all_pos(2,1:npart_tmp), &
                          all_pos(3,1:npart_tmp), &
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

  DO a= 1, npart_tmp, 1

    ! Coordinate velocity of the fluid [c]
    vel_u(0,itr)= 1.0D0
    vel_u(1,itr)= lapse(a)*v_euler_x(a)- shift_x(a)
    vel_u(2,itr)= lapse(a)*v_euler_y(a)- shift_y(a)
    vel_u(3,itr)= lapse(a)*v_euler_z(a)- shift_z(a)

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

  mass_star= binary% get_mass1()
  nu_all= (mass_star/DBLE(npart_tmp))*MSun/amu
  nu= nu_all
  DO a= 1, npart_all
    IF( a < npart_tmp )THEN
      freeze(a)= 0
    ELSE
      freeze(a)= 1
    ENDIF
  ENDDO

  PRINT *, "iterating..."

  apm_iteration: DO itr= 1, 1, 1

    CALL assign_h( npart_all, &
                   npart, &
                   all_pos, h0, &
                   h )

    CALL density_loop( npart, all_pos, &    ! input
                       nu, h, nstar_real )  ! output

    CALL binary% import_id( npart, all_pos(1,1:npart_tmp), &
                            all_pos(2,1:npart_tmp), &
                            all_pos(3,1:npart_tmp), &
                            lapse, shift_x, shift_y, shift_z, &
                            g_xx, g_xy, g_xz, &
                            g_yy, g_yz, g_zz, &
                            baryon_density, &
                            energy_density, &
                            specific_energy, &
                            pressure, &
                            v_euler_x, v_euler_y, v_euler_z )

    DO a= 1, npart_tmp, 1

      ! Coordinate velocity of the fluid [c]
      vel_u(0,itr)= 1.0D0
      vel_u(1,itr)= lapse(a)*v_euler_x(a)- shift_x(a)
      vel_u(2,itr)= lapse(a)*v_euler_y(a)- shift_y(a)
      vel_u(3,itr)= lapse(a)*v_euler_z(a)- shift_z(a)

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
    DO a= 1, npart_tmp, 1

      dNstar= ( nstar_real(a) - nstar_p(a) )/nstar_p(a)
      art_pr(a) = MAX( 1.0D0 + dNstar, 0.1D0 )
      art_pr_max= MAX( art_pr_max, art_pr(a) )
      err_N_max = MAX( err_N_max, ABS(dNstar) )
      err_N_min = MIN( err_N_min, ABS(dNstar) )
      err_N_mean= err_N_mean + ABS(dNstar)

    ENDDO

    nstar_p( npart_tmp+1:npart_all )= 0.0D0
    art_pr ( npart_tmp+1:npart_all )= art_pr_max

    CALL position_correction( npart_all, &
                              all_pos, h, nu_all, art_pr, nstar_real, &
                              correction_pos )

    DO a=1,npart
      pos_tmp= all_pos(:,a) + correction_pos(:,a)
      IF( binary% import_mass_density( &
                              all_pos(1,a), all_pos(2,a), all_pos(3,a) ) > 0 &
          .AND. &
          binary% import_mass_density( &
                              pos_tmp(1), pos_tmp(2), pos_tmp(3) ) > 0 &
      )THEN
        all_pos(:,a)= pos_tmp
      ENDIF
    ENDDO

  ENDDO apm_iteration

  PRINT *, "Iteration completed."
  PRINT *

  ! Now get rid of the ghost particles
  pos= all_pos( :, 1:npart_tmp )
  npart= npart_tmp
  PRINT *, npart

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

  CALL assign_h( npart_all, &
                 npart, &
                 pos, h0, &
                 h )

  CALL density_loop( npart, pos, &    ! input
                     nu, h, nstar_real )  ! output

  CALL binary% import_id( npart, pos(1,1:npart), &
                          pos(2,1:npart), &
                          pos(3,1:npart), &
                          lapse, shift_x, shift_y, shift_z, &
                          g_xx, g_xy, g_xz, &
                          g_yy, g_yz, g_zz, &
                          baryon_density, &
                          energy_density, &
                          specific_energy, &
                          pressure, &
                          v_euler_x, v_euler_y, v_euler_z )

  DO a= 1, npart, 1

    ! Coordinate velocity of the fluid [c]
    vel_u(0,itr)= 1.0D0
    vel_u(1,itr)= lapse(a)*v_euler_x(a)- shift_x(a)
    vel_u(2,itr)= lapse(a)*v_euler_y(a)- shift_y(a)
    vel_u(3,itr)= lapse(a)*v_euler_z(a)- shift_z(a)

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

  nu= nstar_p/nstar_real

  PRINT *, "Baryon number assigned."
  PRINT *

!  PRINT *, "max_nu=", max_nu
!  PRINT *, "min_nu=", min_nu
!  PRINT *, "max_nu/min_nu=", max_nu/min_nu
!  PRINT *

  ! Here setup_uniform_sphere ends

  PRINT *, "** End of PROGRAM proto_apm."
  PRINT *

  CONTAINS


END PROGRAM proto_apm
