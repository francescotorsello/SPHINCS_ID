! File:         submodule_particles_constructor.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (particles_id) stretched_lattice

  !************************************************
  !                                               *
  ! Implementation of the constructor of TYPE     *
  ! particles, and all the PROCEDURES it calls.   *
  !                                               *
  ! FT 19.04.2020                                 *
  !                                               *
  !************************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE place_particles_stretched_lattice

    !************************************************
    !                                               *
    !     *
    !     *
    !                                               *
    ! FT 5.12.2020                                  *
    !                                               *
    !************************************************

    !$ USE OMP_LIB
    USE constants, ONLY: pi, MSun, MSun_geo, km2m, kg2g, lorene2hydrobase, &
                         golden_ratio
    USE matrix,    ONLY: determinant_4x4_matrix
    USE NR,        ONLY: indexx

    IMPLICIT NONE

    INTEGER:: npart1_tmp, npart2_tmp, npart1_radius, npart2_radius, nradii1, &
              nradii2, tmp, cnt, cnt2, npart1_eqplane, npart2_eqplane, &
              nradii1_plane, nradii2_plane, itr3, mass_index, r, th, phi, &
              npart1_half, shell_index
    INTEGER, DIMENSION(:), ALLOCATABLE:: mass_profile_idx
    INTEGER, DIMENSION(:), ALLOCATABLE:: npart_shell, npart_shelleq
    DOUBLE PRECISION:: radius1, radius2, alpha1, alpha2, itr, itr2, &
                       mass_step1, mass_step2, mass_tmp, rad_step, vol_tmp, &
                       mass_shell, rad_coord, lat, long, &
                       lapse, shift_x, shift_y, shift_z, &
                       g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, &
                       det, sq_g, baryon_density, &
                       v_euler_x, v_euler_y, v_euler_z, &
                       v_euler_x_l, v_euler_y_l, v_euler_z_l, &
                       u_euler_t_l, u_euler_x_l, u_euler_y_l, u_euler_z_l, &
                       lorentz_factor, lorentz_factor_rel, &
                       n_t, n_x, n_y, n_z, gamma_euler, &
                       lat_step, long_step, m_p, m_shell, dr, rad, mass, &
                       rad_coord2, center1, mass_element, colat, phase, &
                       central_density
    DOUBLE PRECISION:: g4(0:3,0:3)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: mass_fractions1, &
                                                  mass_fractions2
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: mass_profile
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: shell_radii, shell_masses, &
                                                  alpha

    LOGICAL:: adapt_rad_step, exist

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    TYPE:: colatitude_pos_shell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: colatitudes
    END TYPE

    TYPE(colatitude_pos_shell), DIMENSION(:), ALLOCATABLE:: colatitude_pos

    npart1_tmp= 1D+6
    !npart2_tmp= 1D+5*1.8D0/1.2D0

    THIS% mass1= bns_obj% get_mass1()
    !THIS% mass2= bns_obj% get_mass2()
    radius1    = bns_obj% get_radius1_x_comp()
    !radius2    = bns_obj% get_radius2_x_comp()

    m_p= THIS% mass1/npart1_tmp

    npart1_radius= NINT( radius1* &
                   (npart1_tmp/(4.0D0/3.0D0*pi*radius1**3.0D0))**(1.0D0/3.0D0) )
    !npart2_radius= NINT( radius2* &
    !               (npart2_tmp/(4.0D0/3.0D0*pi*radius2**3.0D0))**(1.0D0/3.0D0) )

    !npart1_eqplane= NINT( pi*radius1**2* &
    !               (npart1_tmp/(4.0D0/3.0D0*pi*radius1**3.0D0))**(2.0D0/3.0D0) )
    !npart2_eqplane= NINT( pi*radius2**2* &
    !               (npart2_tmp/(4.0D0/3.0D0*pi*radius2**3.0D0))**(2.0D0/3.0D0) )

!    nradii1= CEILING( DBLE(npart1_tmp)/DBLE(npart1_radius) )
!    !nradii2= CEILING( DBLE(npart2_tmp)/DBLE(npart2_radius) )
!
!    nradii1_plane= CEILING( SQRT(DBLE(nradii1)) )
!    !nradii2_plane= CEILING( SQRT(DBLE(nradii2)) )
!
!    IF( MOD( nradii1_plane, 2 ) /= 0 ) nradii1_plane= nradii1_plane + 1
!    !IF( MOD( nradii2_plane, 2 ) /= 0 ) nradii2_plane= nradii2_plane + 1
!
!    alpha1= 2.0D0*pi/DBLE(nradii1_plane)
!    !alpha2= 2.0D0*pi/DBLE(nradii2_plane)
!    !alpha1= pi/DBLE(nradii1)*(1 + SQRT( DBLE(1 + 4*nradii1) ))
!    !alpha2= pi/DBLE(nradii2)*(1 + SQRT( DBLE(1 + 4*nradii2) ))
!
!    ALLOCATE( colatitude_pos( nradii1_plane ) )
!    DO itr= 1, nradii1_plane, 1
!      colatitude_pos( itr )= ACOS( 2.0D0*itr/nradii1_plane - 1.0D0 )
!    ENDDO
!
!    PRINT *, "npart1_radius=", npart1_radius
!    !PRINT *, "npart2_radius=", npart2_radius
!    PRINT *, "nradii1_plane=", nradii1_plane
!    !PRINT *, "nradii2_plane=", nradii2_plane
!    PRINT *, "npart1_tmp=", npart1_tmp
!    !PRINT *, "npart2_tmp=", npart2_tmp
!
!    PRINT *, "nradii1=", nradii1
!    !PRINT *, "nradii2=", nradii2
!
!    PRINT *, "nradii1*npart1_radius=", nradii1*npart1_radius
!    !PRINT *, "nradii2*npart2_radius=", nradii2*npart2_radius
!
!    PRINT *, "alpha1=", alpha1
!    !PRINT *, "alpha2=", alpha2
!    PRINT *
!
!    PRINT *, "nradii1_plane*alpha1/(2*pi)=", nradii1_plane*alpha1/(2*pi)
!    !PRINT *, "nradii2_plane*alpha2/(2*pi)=", nradii2_plane*alpha2/(2*pi)
!    PRINT *
!
!    PRINT *, "nradii1_eq=", 2*pi/alpha1
!    PRINT *, "nradii1_mer=", pi/alpha1
!    PRINT *, "nradii1_plane*nradii1_mer=", 2*2*pi/alpha1*pi/alpha1
!    PRINT *, "nradii1_plane*nradii1_mer*npart1_radius=", &
!             2*2*pi/alpha1*pi/alpha1*npart1_radius
!    !PRINT *, "nradii1_eq=", 2*pi/alpha2
!    !PRINT *, "nradii1_mer=", pi/alpha2
!    !PRINT *, "nradii1_plane*nradii1_mer=", 2*2*pi/alpha2*pi/alpha2
!    !PRINT *, "nradii1_plane*nradii1_mer*npart1_radius=", &
!    !         2*2*pi/alpha2*pi/alpha2*n    part2_radius
!    PRINT *
!
!    tmp= 0
!    cnt= 0
!    cnt2=0
!    DO itr= 0, 2*pi - alpha1, alpha1
!      cnt= cnt + 1
!      DO itr2= alpha1/2, pi-alpha1/2, alpha1
!        IF(itr==0)THEN
!          cnt2=cnt2+1
!        ENDIF
!        tmp= tmp + npart1_radius
!      ENDDO
!    ENDDO
!    PRINT *, "itr=", itr/(2*pi-alpha1)
!    PRINT *, "itr2=", itr2/(pi-alpha1/2)
!    PRINT *, "tmp=", tmp*2
!    PRINT *, "cnt=", cnt
!    PRINT *, "cnt2=", cnt2
!    !PRINT *, "cnt*2=", cnt*2
!    !PRINT *, (nradii1 - cnt*2)
!    !PRINT *, (nradii1 - cnt*2)*npart1_radius + tmp*2
!    PRINT *
!
!    ALLOCATE( mass_fractions1(npart1_radius) )
!    ALLOCATE( mass_fractions2(npart2_radius) )
!
!    mass_step1= THIS% mass1/npart1_radius
!    mass_step2= THIS% mass2/npart2_radius
!
!    DO itr= 1, npart1_radius, 1
!      mass_fractions1(itr)= itr*mass_step1
!    ENDDO
!    !DO itr= 1, npart2_radius, 1
!    !  mass_fractions2(itr)= itr*mass_step2
!    !ENDDO
!
!    IF( mass_fractions1(npart1_radius) /= THIS% mass1 )THEN
!      PRINT *, "** ERROR in ! The mass partition for star 1 is incorrect."
!      STOP
!    ENDIF
!    !IF( mass_fractions2(npart2_radius) /= THIS% mass2 )THEN
!    !  PRINT *, "** ERROR in ! The mass partition for star 2 is incorrect."
!    !  STOP
!    !ENDIF

    IF(.NOT.ALLOCATED( shell_radii ))THEN
      ALLOCATE( shell_radii( npart1_radius ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array shell_radii in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array pos in SUBROUTINE" &
      !                // "place_particles_3D_lattice." )
    ENDIF
    IF(.NOT.ALLOCATED( shell_masses ))THEN
      ALLOCATE( shell_masses( npart1_radius ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array shell_masses in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array pos in SUBROUTINE" &
      !                // "place_particles_3D_lattice." )
    ENDIF

    PRINT *, "Integrating the baryon mass density to get the mass profile..."
    PRINT *

    rad_step= radius1/500.0D0
    lat_step= pi/2.0D0/250.0D0
    long_step= 2.0D0*pi/500.0D0
    center1= bns_obj% get_center1_x()
    central_density= bns_obj% get_rho_center1()
    CALL bns_obj% integrate_baryon_mass_density( center1, radius1, &
                                                 central_density, &
                                                 rad_step, lat_step, long_step, &
                                                 mass, mass_profile, &
                                                 mass_profile_idx )

    DO itr= 1, npart1_radius, 1

      shell_radii( itr )= radius1*itr/npart1_radius

    ENDDO
    shell_index= 1
    itr2= 0
    DO itr= 0, NINT(radius1/rad_step), 1

      IF( mass_profile( 1, mass_profile_idx(itr) ) &
          >= shell_radii( shell_index ) )THEN

        shell_masses( shell_index )= SUM( mass_profile( 2, &
                      mass_profile_idx(itr2):mass_profile_idx(itr) ), DIM= 1 )

        IF( shell_index == npart1_radius )THEN
          EXIT
        ELSE
         itr2= itr + 1
         shell_index= shell_index + 1
        ENDIF

      ENDIF

    ENDDO

    ALLOCATE( npart_shell(npart1_radius) )
    ALLOCATE( npart_shelleq(npart1_radius) )
    ALLOCATE( alpha(npart1_radius) )
    ALLOCATE( colatitude_pos( npart1_radius ) )

    DO itr= 1, npart1_radius, 1

      npart_shelleq( itr )= CEILING( SQRT(DBLE(shell_masses( itr )/m_p)) )
      IF( MOD( npart_shelleq( itr ), 2 ) /= 0 )THEN
        npart_shelleq( itr )= npart_shelleq( itr ) + 1
      ENDIF
      alpha( itr )= 2.0D0*pi/DBLE(npart_shelleq( itr ))
      npart_shell( itr )= npart_shelleq( itr )**2.0D0

      ALLOCATE( colatitude_pos( itr )% colatitudes( npart_shelleq( itr ) ) )

      DO itr2= 1, npart_shelleq( itr ), 1

        colatitude_pos( itr )% colatitudes( itr2 )= &
                      ACOS( 2.0D0*itr2/(npart_shelleq( itr ) + 1.0D0 ) - 1.0D0 )

      ENDDO

    ENDDO
    THIS% npart= SUM( npart_shell, DIM= 1 )

    !PRINT *, colatitude_pos( 1 )% colatitudes
    !STOP

    ! Place the particles for one star only, since the subroutine will place
    ! particles for one star

    ! Allocating the memory for the array pos( 3, npart_temp )
    ! Note that after determining npart, the array pos is reshaped into
    ! pos( 3, npart )
    IF(.NOT.ALLOCATED( THIS% pos ))THEN
      ALLOCATE( THIS% pos( 3, THIS% npart ), &
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

    ! Latitude first, longitude second
    !mass_index= 1
    !shell_radii= 1.0D0
    !mass_tmp= 0.0D0
    !vol_tmp= 0.0D0

    ! Assume same mass for particles
    !m_p             = THIS% mass1/npart1_tmp
    !m_shell         = m_p*nradii1
    !shell_radii( npart1_radius )= radius1*0.99D0!(3.0D0*(THIS% mass1/1000)/(4.0D0*pi&
                      !*bns_obj% import_mass_density( bns_obj% get_center1_x(), &
                      !0.0D0, 0.0D0)))**(1.0D0/3.0D0)
    !rad_coord       = shell_radii( npart1_radius )
    !PRINT *, rad_coord
    !PRINT *, bns_obj% get_center1_x() + rad_coord
    !PRINT *, m_shell*npart1_radius
    !PRINT *

  !   DO itr= npart1_radius - 1, 1, -1
  !
  !     CALL bns_obj% import_id( bns_obj% get_center1_x() &
  !                              + shell_radii( itr + 1 ), &
  !                              0.0D0, &
  !                              0.0D0, &
  !                              lapse, shift_x, shift_y, shift_z, &
  !                              g_xx, baryon_density, &
  !                              v_euler_x, v_euler_y, v_euler_z, &
  !                              gamma_euler )
  !
  !     ! Compute covariant spatial fluid velocity (metric is diagonal and
  !     ! conformally flat)
  !     !v_euler_x_l= g_xx*v_euler_x
  !     !v_euler_y_l= g_xx*v_euler_y
  !     !v_euler_z_l= g_xx*v_euler_z
  !     !
  !     !! Compute the corresponding Lorentz factor
  !     !lorentz_factor= 1.0D0/SQRT( 1.0D0 - ( v_euler_x_l*v_euler_x &
  !     !                                    + v_euler_y_l*v_euler_y &
  !     !                                    + v_euler_z_l*v_euler_z ) )
  !     !
  !     !! Compute covariant fluid 4-velocity
  !     !u_euler_t_l= lorentz_factor *( - lapse + v_euler_x_l*shift_x &
  !     !                                       + v_euler_y_l*shift_y &
  !     !                                       + v_euler_z_l*shift_z )
  !     !u_euler_x_l= lorentz_factor*v_euler_x_l
  !     !u_euler_y_l= lorentz_factor*v_euler_y_l
  !     !u_euler_z_l= lorentz_factor*v_euler_z_l
  !     !
  !     !! Compute vector normal to spacelike hypersurface
  !     !! (4-velocity of the Eulerian observer)
  !     !n_t= 1.0D0/lapse
  !     !n_x= - shift_x/lapse
  !     !n_y= - shift_y/lapse
  !     !n_z= - shift_z/lapse
  !     !
  !     !! Compute relative Lorentz factor between 4-velocity of the fluid
  !     !! wrt the Eulerian observer and the 4-velocity of the Eulerian observer
  !     !lorentz_factor_rel= - ( n_t*u_euler_t_l + n_x*u_euler_x_l &
  !     !                      + n_y*u_euler_y_l + n_z*u_euler_z_l )
  !     !
  !     !PRINT *, lorentz_factor_rel
  !     !PRINT *, gamma_euler
  !     !STOP
  !
  !     ! Compute square root of the determinant of the spatial metric
  !     sq_g= g_xx*SQRT( g_xx )
  !
  !     dr= m_shell/( baryon_density*sq_g*gamma_euler &
  !                   *4.0D0*pi*shell_radii( itr + 1 )**2.0D0 )
  !
  !     PRINT *, m_shell, baryon_density, sq_g, gamma_euler, rad_coord
  !     PRINT *, "dr=", dr
  !     shell_radii( itr )= shell_radii( itr + 1 ) - dr
  !     PRINT *, "shell_radii=", shell_radii( itr )
  !
  !     PRINT *
  !   ENDDO
  !   PRINT *, shell_radii
  !   STOP

    !radius_loop: DO rad_coord= rad_step, radius1, rad_step
    !
    !  ! The definition of the baryon mass for the LORENE ID is in eq.(69)
    !  ! of Gourgoulhon et al., PRD 63 064029 (2001)
    !
    !  CALL bns_obj% import_id( &
    !                      bns_obj% get_center1_x() + rad_coord - rad_step, &
    !                      0.0D0, &
    !                      0.0D0, &
    !                      lapse, shift_x, shift_y, shift_z, &
    !                      g_xx, baryon_density, &
    !                      v_euler_x, v_euler_y, v_euler_z, &
    !                      gamma_euler )
    !
    !  ! Compute covariant spatial fluid velocity (metric is diagonal and
    !  ! conformally flat)
    !  v_euler_x_l= g_xx*v_euler_x
    !  v_euler_y_l= g_xx*v_euler_y
    !  v_euler_z_l= g_xx*v_euler_z
    !
    !  ! Compute the corresponding Lorentz factor
    !  lorentz_factor= 1.0D0/SQRT( 1.0D0 - ( v_euler_x_l*v_euler_x &
    !                                      + v_euler_y_l*v_euler_y &
    !                                      + v_euler_z_l*v_euler_z ) )
    !
    !  ! Compute covariant fluid 4-velocity
    !  u_euler_t_l= lorentz_factor *( - lapse + v_euler_x_l*shift_x &
    !                                         + v_euler_y_l*shift_y &
    !                                         + v_euler_z_l*shift_z )
    !  u_euler_x_l= lorentz_factor*v_euler_x_l
    !  u_euler_y_l= lorentz_factor*v_euler_y_l
    !  u_euler_z_l= lorentz_factor*v_euler_z_l
    !
    !  ! Compute vector normal to spacelie hypersurface
    !  ! (4-velocity of the Eulerian observer)
    !  n_t= 1.0D0/lapse
    !  n_x= - shift_x/lapse
    !  n_y= - shift_y/lapse
    !  n_z= - shift_z/lapse
    !
    !  ! Compute relative Lorentz factor between 4-velocity of the fluid
    !  ! wrt the Eulerian observer and the 4-velocity of the Eulerian observer
    !  lorentz_factor_rel= - ( n_t*u_euler_t_l + n_x*u_euler_x_l &
    !                        + n_y*u_euler_y_l + n_z*u_euler_z_l )
    !
    !  !IF( gamma_euler /= lorentz_factor_rel )THEN
    !  !  PRINT *, "gamma_euler=", gamma_euler
    !  !  PRINT *, "lorentz_factor_rel=", lorentz_factor_rel
    !  !  STOP
    !  !ENDIF
    !  !STOP
    !
    !  ! Compute square root of the determinant of the spatial metric
    !  sq_g= g_xx*SQRT( g_xx )
    !
    !  ! Compute mass of a spherical shell the the assumption of
    !  ! spherical symmetry
    !  mass_shell= 4.0D0*pi*(rad_coord**2.0D0)*rad_step &
    !              *sq_g*gamma_euler*baryon_density
    !              !*(4.0D0/3.0D0)*pi*sq_g &
    !              !*(rad_coord**3.0D0 - (rad_coord - rad_step)**3.0D0)
    !
    !  mass_tmp= mass_tmp + mass_shell
    !
    !  vol_tmp= vol_tmp &
    !           + 4.0D0*pi*rad_coord**2.0D0*rad_step*sq_g
    !
    !  PRINT *, bns_obj% get_center1_x() + rad_coord*COS(0.0D0)*COS(0.0D0), &
    !           rad_coord, &
    !           bns_obj% import_mass_density( &
    !                 bns_obj% get_center1_x() + rad_coord, 0.0D0, 0.0D0 ), &
    !           + 4.0D0*pi*rad_coord**2.0D0*rad_step*sq_g, &
    !           mass_shell, mass_tmp
    !
    !  IF( mass_tmp >= mass_fractions1( mass_index ) )THEN
    !    shell_radii( mass_index )= rad_coord
    !    IF( mass_index == npart1_radius )THEN
    !      EXIT
    !    ELSE
    !      mass_index= mass_index + 1
    !    ENDIF
    !  ENDIF
    !
    !ENDDO radius_loop
  !rad_step= radius1/1000.0D0
  !lat_step= pi/2.0D0/250.0D0
  !long_step= 2.0D0*pi/1000.0D0



  PRINT *, "Print mass profile to file..."

  finalnamefile= "mass_profile.dat"

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

  write_data_loop: DO itr = 1, NINT(radius1/rad_step), 1

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      mass_profile( 1, mass_profile_idx(itr) ), &
      mass_profile( 2, mass_profile_idx(itr) ), &
      mass_profile( 3, mass_profile_idx(itr) )

    IF( ios > 0 )THEN
      PRINT *, "...error when writing the arrays in " &
               // TRIM(finalnamefile), ". The error message is", err_msg
      STOP
    ENDIF

  ENDDO write_data_loop

  CLOSE( UNIT= 2 )

  PRINT *, "Print shell radii to file..."

  !mass_index= 1
  !DO itr = 1, NINT(radius1/rad_step), 1
  !
  !  IF( mass_profile( 3, mass_profile_idx(itr) ) &
  !      >= mass_fractions1( mass_index ) )THEN
  !   shell_radii( mass_index )= mass_profile( 1, mass_profile_idx(itr) )
  !   IF( mass_index == npart1_radius )THEN
  !     EXIT
  !   ELSE
  !    mass_index= mass_index + 1
  !   ENDIF
  !  ENDIF
  !
  !ENDDO

  finalnamefile= "shell_radii.dat"

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

  DO itr = 1, npart1_radius, 1

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
      shell_radii( itr )

    IF( ios > 0 )THEN
      PRINT *, "...error when writing the arrays in " &
               // TRIM(finalnamefile), ". The error message is", err_msg
      STOP
    ENDIF

  ENDDO

  CLOSE( UNIT= 2 )

  PRINT *, "Assigning first half of particle positions..."

    PRINT *, SUM( npart_shell, DIM= 1 )

    THIS% npart= 0
    THIS% pos= 0
    phase= 0
    DO r= 1, npart1_radius, 1
      !phase= phase + (r - 1.0D0)*2*pi/(SQRT(2.0D0)/(7.0D0/5.0D0)*nradii1_plane)
      phase= phase + r*alpha(r)/golden_ratio
      DO th= 1, npart_shelleq( r )/2, 1
        DO phi= 1, npart_shelleq( r ), 1 !part_shelleq( r ) is even, see above

          THIS% npart= THIS% npart + 1
          THIS% pos( 1, THIS% npart )= &
                shell_radii(r)*COS(phase + phi*alpha(r)) &
                              *SIN(colatitude_pos(r)% colatitudes(th))
          THIS% pos( 2, THIS% npart )= &
                shell_radii(r)*SIN(phase + phi*alpha(r)) &
                              *SIN(colatitude_pos(r)% colatitudes(th))
          THIS% pos( 3, THIS% npart )= &
                shell_radii(r)*COS(colatitude_pos(r)% colatitudes(th))

        ENDDO
      ENDDO
    ENDDO
    PRINT *, "npart=", THIS% npart

    PRINT *, "Mirroring particles..."

    PRINT *, "npart/2=", THIS% npart
    npart1_half= THIS% npart
    DO itr= 1, npart1_half, 1
      THIS% npart= THIS% npart + 1
      THIS% pos( 1, THIS% npart )=   THIS% pos( 1, itr )
      THIS% pos( 2, THIS% npart )=   THIS% pos( 2, itr )
      THIS% pos( 3, THIS% npart )= - THIS% pos( 3, itr )
    ENDDO
    PRINT *, "npart=", THIS% npart

  PRINT *, "Printing particle positions to file..."

    finalnamefile= "tmp_pos.dat"

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

    DO itr = 1, THIS% npart, 1

      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        THIS% pos( 1, itr ), THIS% pos( 2, itr ), THIS% pos( 3, itr )

      IF( ios > 0 )THEN
        PRINT *, "...error when writing the arrays in " &
                 // TRIM(finalnamefile), ". The error message is", err_msg
        STOP
      ENDIF

    ENDDO

    CLOSE( UNIT= 2 )

  PRINT *, "STOP..."

    STOP

    IF( THIS% mass1 > THIS% mass2 )THEN

      ! mass_ratio < 1
      THIS% mass_ratio= THIS% mass2/THIS% mass1


    ELSE

      ! mass_ratio < 1
      THIS% mass_ratio= THIS% mass1/THIS% mass2


    ENDIF

  END PROCEDURE place_particles_stretched_lattice


END SUBMODULE stretched_lattice
