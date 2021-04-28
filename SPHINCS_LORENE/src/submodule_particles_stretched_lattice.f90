! File:         submodule_particles_constructor.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (particles_id) stretched_lattice

  !************************************************
  !                                               *
  !
  !*
  !                                               *
  ! FT 19.04.2021                                 *
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
    ! FT 19.04.2021                                 *
    !                                               *
    !************************************************

    !$ USE OMP_LIB
    USE constants, ONLY: pi, MSun, MSun_geo, km2m, kg2g, lorene2hydrobase, &
                         golden_ratio, third
    USE matrix,    ONLY: determinant_4x4_matrix
    USE NR,        ONLY: indexx

    IMPLICIT NONE

    INTEGER:: n_shells, itr2, itr3, mass_index, npart_half, npart_tmp, cnt, &
              shell_index, r, th, phi, i_shell, npart_test, npart_shell_tmp, &
              cnt2, rel_sign, cnt3, dim_seed
    INTEGER, PARAMETER:: max_length= 5D+6
    INTEGER, DIMENSION(:), ALLOCATABLE:: mass_profile_idx, seed
    INTEGER, DIMENSION(:), ALLOCATABLE:: npart_shell, npart_shelleq

    DOUBLE PRECISION:: xtemp, ytemp, ztemp, m_p, central_density, &
                       dr, dth, dphi, phase, phase_th, mass, baryon_density, &
                       dr_shells, dth_shells, dphi_shells, col, rad, &
                       g_xx, gamma_euler, proper_volume, mass_test, mass_test2,&
                       proper_volume_test, npart_shell_kept, &
                       upper_bound, lower_bound, rand_num, rand_num2, &
                       upper_bound_tmp, lower_bound_tmp
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: mass_profile
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: shell_radii, shell_masses, &
                                                  alpha, m_parts, vol_shell, &
                                                  vol_shell2, mass_shell, &
                                                  mass_shell2

    LOGICAL:: exist, high_mass, low_mass, kept_all

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    TYPE:: colatitude_pos_shell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: colatitudes
    END TYPE

    TYPE(colatitude_pos_shell), DIMENSION(:), ALLOCATABLE:: colatitude_pos

    TYPE:: pos_on_shells
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos_shell
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pos_th
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pos_phi
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pvol_shell
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pvol_shell2
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: g_xx
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: baryon_density
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: gamma_euler
    END TYPE

    TYPE(pos_on_shells), DIMENSION(:), ALLOCATABLE:: pos_shells

    CALL RANDOM_SEED( SIZE= dim_seed )
    ALLOCATE( seed( dim_seed ) )
    seed( 1 )= 0
    seed( 2 )= 1
    DO itr= 3, dim_seed
      seed( itr )= seed( itr - 1 ) + seed( itr - 2 )
    ENDDO
    CALL RANDOM_SEED( PUT= seed )

    !npart_approx= 1D+6
    !npart2_tmp= 1D+5*1.8D0/1.2D0

    !THIS% mass1= bns_obj% get_mass1()
    !THIS% mass2= bns_obj% get_mass2()
    !radius    = bns_obj% get_radius1_x_comp()
    !radius2    = bns_obj% get_radius2_x_comp()

    m_p= mass_star/npart_approx

    n_shells= NINT( radius* &
                  (npart_approx/(4.0D0/3.0D0*pi*radius**3.0D0))**third )
    !npart2_radius= NINT( radius2* &
    !               (npart2_tmp/(4.0D0/3.0D0*pi*radius2**3.0D0))**(1.0D0/3.0D0) )

    !npart1_eqplane= NINT( pi*radius**2* &
    !               (npart1_tmp/(4.0D0/3.0D0*pi*radius1**3.0D0))**(2.0D0/3.0D0) )
    !npart2_eqplane= NINT( pi*radius2**2* &
    !               (npart2_tmp/(4.0D0/3.0D0*pi*radius2**3.0D0))**(2.0D0/3.0D0) )

!    nradii1= CEILING( DBLE(npart1_tmp)/DBLE(n_shells) )
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
!    PRINT *, "n_shells=", n_shells
!    !PRINT *, "npart2_radius=", npart2_radius
!    PRINT *, "nradii1_plane=", nradii1_plane
!    !PRINT *, "nradii2_plane=", nradii2_plane
!    PRINT *, "npart1_tmp=", npart1_tmp
!    !PRINT *, "npart2_tmp=", npart2_tmp
!
!    PRINT *, "nradii1=", nradii1
!    !PRINT *, "nradii2=", nradii2
!
!    PRINT *, "nradii1*n_shells=", nradii1*n_shells
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
!    PRINT *, "nradii1_plane*nradii1_mer*n_shells=", &
!             2*2*pi/alpha1*pi/alpha1*n_shells
!    !PRINT *, "nradii1_eq=", 2*pi/alpha2
!    !PRINT *, "nradii1_mer=", pi/alpha2
!    !PRINT *, "nradii1_plane*nradii1_mer=", 2*2*pi/alpha2*pi/alpha2
!    !PRINT *, "nradii1_plane*nradii1_mer*n_shells=", &
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
!        tmp= tmp + n_shells
!      ENDDO
!    ENDDO
!    PRINT *, "itr=", itr/(2*pi-alpha1)
!    PRINT *, "itr2=", itr2/(pi-alpha1/2)
!    PRINT *, "tmp=", tmp*2
!    PRINT *, "cnt=", cnt
!    PRINT *, "cnt2=", cnt2
!    !PRINT *, "cnt*2=", cnt*2
!    !PRINT *, (nradii1 - cnt*2)
!    !PRINT *, (nradii1 - cnt*2)*n_shells + tmp*2
!    PRINT *
!
!    ALLOCATE( mass_fractions1(n_shells) )
!    ALLOCATE( mass_fractions2(npart2_radius) )
!
!    mass_step1= THIS% mass1/n_shells
!    mass_step2= THIS% mass2/npart2_radius
!
!    DO itr= 1, n_shells, 1
!      mass_fractions1(itr)= itr*mass_step1
!    ENDDO
!    !DO itr= 1, npart2_radius, 1
!    !  mass_fractions2(itr)= itr*mass_step2
!    !ENDDO
!
!    IF( mass_fractions1(n_shells) /= THIS% mass1 )THEN
!      PRINT *, "** ERROR in ! The mass partition for star 1 is incorrect."
!      STOP
!    ENDIF
!    !IF( mass_fractions2(npart2_radius) /= THIS% mass2 )THEN
!    !  PRINT *, "** ERROR in ! The mass partition for star 2 is incorrect."
!    !  STOP
!    !ENDIF

    IF(.NOT.ALLOCATED( shell_radii ))THEN
      ALLOCATE( shell_radii( n_shells ), STAT= ios, &
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
      ALLOCATE( shell_masses( n_shells ), STAT= ios, &
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
    IF(.NOT.ALLOCATED( vol_shell ))THEN
      ALLOCATE( vol_shell( n_shells ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array vol_shell in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array pos in SUBROUTINE" &
      !                // "place_particles_3D_lattice." )
    ENDIF
    IF(.NOT.ALLOCATED( vol_shell2 ))THEN
      ALLOCATE( vol_shell2( n_shells ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array vol_shell2 in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array pos in SUBROUTINE" &
      !                // "place_particles_3D_lattice." )
    ENDIF
    IF(.NOT.ALLOCATED( mass_shell ))THEN
      ALLOCATE( mass_shell( n_shells ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array vol_shell in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array pos in SUBROUTINE" &
      !                // "place_particles_3D_lattice." )
    ENDIF
    IF(.NOT.ALLOCATED( mass_shell2 ))THEN
      ALLOCATE( mass_shell2( n_shells ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array vol_shell2 in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array pos in SUBROUTINE" &
      !                // "place_particles_3D_lattice." )
    ENDIF
    IF(.NOT.ALLOCATED( m_parts ))THEN
      ALLOCATE( m_parts( n_shells ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array m_parts in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array pos in SUBROUTINE" &
      !                // "place_particles_3D_lattice." )
    ENDIF

    PRINT *, " ** Integrating the baryon mass density to get the mass profile..."
    PRINT *

    CALL OMP_SET_NUM_THREADS(80)

    dr             = radius/500.0D0
    dth            = pi/2.0D0/250.0D0
    dphi           = 2.0D0*pi/500.0D0
    central_density= bns_obj% get_rho_center1()
    CALL bns_obj% integrate_baryon_mass_density( center, radius, &
                                                 central_density, &
                                                 dr, dth, dphi, &
                                                 mass, mass_profile, &
                                                 mass_profile_idx )

    mass_profile( 2:3, : )= mass_profile( 2:3, : )*mass_star/mass

    PRINT *, "mass_profile( 2, mass_profile_idx(NINT(radius/dr)) )= ", &
             mass_profile( 2, mass_profile_idx(NINT(radius/dr)) )

    shell_radii= 0.0D0
    DO itr= 1, n_shells, 1

      shell_radii( itr )= (radius/DBLE(n_shells))*DBLE( itr - 1/2 )

    ENDDO
    shell_index= 1
    itr2= 0
    shell_masses= 0.0D0
    DO itr= 0, NINT(radius/dr), 1

      IF( mass_profile( 1, mass_profile_idx(itr) ) &
          >= shell_radii( shell_index ) + radius/DBLE(2*n_shells) )THEN

        IF( shell_index == n_shells )THEN

          shell_masses( shell_index )= SUM( mass_profile( 2, &
           mass_profile_idx(itr2):mass_profile_idx(NINT(radius/dr)) ), DIM= 1 )

          EXIT

        ELSE

         shell_masses( shell_index )= SUM( mass_profile( 2, &
                       mass_profile_idx(itr2):mass_profile_idx(itr) ), DIM= 1 )

         itr2= itr + 1
         shell_index= shell_index + 1

        ENDIF

      ENDIF

    ENDDO
    IF( ABS( SUM( shell_masses, DIM= 1 ) - mass_star )/mass_star > 5.0D-3 )THEN
      PRINT *, " ** The sum of the masses of the shells do not add up to the ", &
               "mass of the star. Stopping..."
      PRINT *, "SUM( shell_masses )= ", SUM( shell_masses, DIM=1 )
      PRINT *, "mass_star= ", mass_star
      PRINT *
      STOP
    ENDIF

    ALLOCATE( npart_shell( n_shells ) )
    ALLOCATE( npart_shelleq( n_shells ) )
    ALLOCATE( alpha( n_shells ) )
    ALLOCATE( colatitude_pos( n_shells ) )
    ALLOCATE( pos_shells( n_shells ) )

    npart_shelleq= 0
    DO itr= 1, n_shells, 1

      npart_shelleq( itr )= CEILING( SQRT(DBLE(shell_masses( itr )/m_p)) )
      !IF( itr == n_shells - 1 ) npart_shelleq( itr )= &
      !                      npart_shelleq( n_shells - 2 )
      !                      !MAXVAL( npart_shelleq )
      !                      !NINT( 1.5D0*npart_shelleq( itr ) )
      !IF( itr == n_shells ) npart_shelleq( itr )= &
      !                      !npart_shelleq( n_shells - 2 )
      !                      MAXVAL( npart_shelleq )
      !                      !NINT( 4.0D0*npart_shelleq( itr ) )
      IF( MOD( npart_shelleq( itr ), 2 ) /= 0 )THEN
        npart_shelleq( itr )= npart_shelleq( itr ) + 1
      ENDIF
      IF( npart_shelleq( itr ) == 0 )THEN
        alpha( itr )= 0.0D0
      ELSE
        alpha( itr )= 2.0D0*pi/DBLE(npart_shelleq( itr ))
      ENDIF
      npart_shell( itr )= npart_shelleq( itr )**2.0D0

      ALLOCATE( colatitude_pos( itr )% colatitudes( npart_shelleq( itr ) ) )

      DO itr2= 1, npart_shelleq( itr ), 1

        colatitude_pos( itr )% colatitudes( itr2 )= &
                ACOS( 2.0D0*itr2/(npart_shelleq( itr ) + 1.0D0 ) - 1.0D0 )

      ENDDO

      ALLOCATE( pos_shells( itr )% pos_shell( 3, npart_shell( itr ) ) )
      ALLOCATE( pos_shells( itr )% pvol_shell( npart_shell( itr ) ) )
      ALLOCATE( pos_shells( itr )% g_xx( npart_shell( itr ) ) )
      ALLOCATE( pos_shells( itr )% baryon_density( npart_shell( itr ) ) )
      ALLOCATE( pos_shells( itr )% gamma_euler( npart_shell( itr ) ) )

    ENDDO
    npart_tmp= SUM( npart_shell, DIM= 1 )

    PRINT *, "npart before importing ID=", SUM( npart_shell, DIM= 1 )

    DO r= 1, n_shells, 1
      IF( npart_shell( r ) == 0 )THEN
        m_parts( r )= m_parts( r - 1 )
      ELSE
        m_parts( r )= shell_masses( r )/npart_shell( r )
      ENDIF
    ENDDO
    PRINT *, npart_shell
    PRINT *
    PRINT *, m_parts
    PRINT *
    PRINT *, "nu_ratio= ", MAXVAL(m_parts)/MINVAL(m_parts)
    PRINT *

    !PRINT *, colatitude_pos( 1 )% colatitudes
    !STOP

    ! Place the particles for one star only, since the subroutine will place
    ! particles for one star

    ! Allocating the memory for the array pos( 3, npart_temp )
    ! Note that after determining npart, the array pos is reshaped into
    ! pos( 3, npart )
    !IF(.NOT.ALLOCATED( pvol ))THEN
    !  ALLOCATE( pvol( npart_tmp ), &
    !            STAT= ios, ERRMSG= err_msg )
    !  IF( ios > 0 )THEN
    !     PRINT *, "...allocation error for array pvol in SUBROUTINE" &
    !              // "place_particles_. ", &
    !              "The error message is", err_msg
    !     STOP
    !  ENDIF
    !  !CALL test_status( ios, err_msg, &
    !  !                "...allocation error for array pos in SUBROUTINE" &
    !  !                // "place_particles_3D_lattice." )
    !ENDIF

    ! Latitude first, longitude second
    !mass_index= 1
    !shell_radii= 1.0D0
    !mass_tmp= 0.0D0
    !vol_tmp= 0.0D0

    ! Assume same mass for particles
    !m_p             = THIS% mass1/npart1_tmp
    !m_shell         = m_p*nradii1
    !shell_radii( n_shells )= radius1*0.99D0!(3.0D0*(THIS% mass1/1000)/(4.0D0*pi&
                      !*bns_obj% import_mass_density( bns_obj% get_center1_x(), &
                      !0.0D0, 0.0D0)))**(1.0D0/3.0D0)
    !rad_coord       = shell_radii( n_shells )
    !PRINT *, rad_coord
    !PRINT *, bns_obj% get_center1_x() + rad_coord
    !PRINT *, m_shell*n_shells
    !PRINT *

  !   DO itr= n_shells - 1, 1, -1
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

  write_data_loop: DO itr = 1, NINT(radius/dr), 1

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
  !   IF( mass_index == n_shells )THEN
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

  DO itr = 1, n_shells, 1

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

    !DO r= 1, n_shells, 1
    !  pos_shells(r)% pos_shell= 0.0D0
    !  pos_shells(r)% pos_phi= -1.0D0
    !  pos_shells(r)% pos_th= -1.0D0
    !  pos_shells(r)% pvol_shell= 0.0D0
    !  pos_shells(r)% pvol_shell2= 0.0D0
    !  pos_shells(r)% g_xx= 0.0D0
    !  pos_shells(r)% baryon_density= 0.0D0
    !  pos_shells(r)% gamma_euler= 0.0D0
    !  m_parts( r )= m_p
    !  npart_shelleq( r )= CEILING( SQRT(DBLE(shell_masses( r )/m_parts( r ))) )
    !ENDDO

    DO r= 1, n_shells, 1
      IF( ALLOCATED( pos_shells( r )% pos_shell ) ) &
        DEALLOCATE( pos_shells( r )% pos_shell )
        PRINT *, "0.1"
      IF( ALLOCATED( pos_shells( r )% pvol_shell ) ) &
        DEALLOCATE( pos_shells( r )% pvol_shell )
                PRINT *, "0.2"
      IF( ALLOCATED( pos_shells( r )% pvol_shell2 ) ) &
        DEALLOCATE( pos_shells( r )% pvol_shell2 )
                PRINT *, "0.3"
      IF( ALLOCATED( pos_shells( r )% g_xx ) )&
        DEALLOCATE( pos_shells( r )% g_xx )
                PRINT *, "0.4"
      IF( ALLOCATED( pos_shells( r )% baryon_density ) ) &
        DEALLOCATE( pos_shells( r )% baryon_density )
                PRINT *, "0.5"
      IF( ALLOCATED( pos_shells( r )% gamma_euler ) ) &
        DEALLOCATE( pos_shells( r )% gamma_euler )
                      PRINT *, "0.6"
      IF( ALLOCATED( pos_shells( r )% pos_th ) ) &
        DEALLOCATE( pos_shells( r )% pos_th )
                      PRINT *, "0.7"
      IF( ALLOCATED( pos_shells( r )% pos_phi ) ) &
        DEALLOCATE( pos_shells( r )% pos_phi )

      ALLOCATE( pos_shells( r )% pos_shell     ( 3, max_length ) )
      ALLOCATE( pos_shells( r )% pvol_shell    ( max_length ) )
      ALLOCATE( pos_shells( r )% pvol_shell2   ( max_length ) )
      ALLOCATE( pos_shells( r )% g_xx          ( max_length ) )
      ALLOCATE( pos_shells( r )% baryon_density( max_length ) )
      ALLOCATE( pos_shells( r )% gamma_euler   ( max_length ) )
      ALLOCATE( pos_shells( r )% pos_th        ( max_length ) )
      ALLOCATE( pos_shells( r )% pos_phi       ( max_length ) )

      pos_shells(r)% pos_shell= 0.0D0
      pos_shells(r)% pos_phi= -1.0D0
      pos_shells(r)% pos_th= -1.0D0
      pos_shells(r)% pvol_shell= 0.0D0
      pos_shells(r)% pvol_shell2= 0.0D0
      pos_shells(r)% g_xx= 0.0D0
      pos_shells(r)% baryon_density= 0.0D0
      pos_shells(r)% gamma_euler= 0.0D0
      m_parts( r )= m_p
      npart_shelleq( r )= CEILING( SQRT(DBLE(shell_masses( r )/m_parts( r ))) )

    ENDDO

    pos  = 0.0D0
    pmass = 0.0D0
    phase= 0.0D0
    proper_volume= 0.0D0
    vol_shell= 0.0D0
    vol_shell2= 0.0D0
    dr_shells= radius/n_shells
    npart_out= 0
    upper_bound= 1.025D0
    lower_bound= 0.975D0
    upper_bound_tmp= upper_bound
    lower_bound_tmp= lower_bound
    r= 1
    cnt2= 0

    !----------------------!
    !--  Main iteration  --!
    !----------------------!

    CALL OMP_SET_NUM_THREADS(80)

    DO

PRINT *, "Start of iteration, shell ", r, "iteration ", cnt2 + 1
PRINT *
      !npart_shelleq( r )= NINT( correction*npart_shelleq( r ) )
      IF( MOD( npart_shelleq( r ), 2 ) /= 0 )THEN
        npart_shelleq( r )= npart_shelleq( r ) + 1
      ENDIF
      IF( npart_shelleq( r ) == 0 )THEN
        alpha( r )= 0.0D0
      ELSE
        alpha( r )= 2.0D0*pi/DBLE(npart_shelleq( r ))
      ENDIF
      npart_shell( r )= npart_shelleq( r )**2.0D0

      IF( ALLOCATED( colatitude_pos( r )% colatitudes ) ) &
        DEALLOCATE( colatitude_pos( r )% colatitudes )
      ALLOCATE( colatitude_pos( r )% colatitudes( npart_shelleq( r ) ) )

      DO itr2= 1, npart_shelleq( r ), 1

        colatitude_pos( r )% colatitudes( itr2 )= &
                      ACOS( 2.0D0*itr2/(npart_shelleq( r ) + 1.0D0 ) - 1.0D0 )
        !colatitude_pos( r )% colatitudes( itr2 )= &
        !              colatitude_pos( r )% colatitudes( itr2 ) &
        !              *( 1 + rel_sign*0.05D0*phase_th )

      ENDDO

      PRINT *, "Before deallocation, shell ", r, "iteration ", cnt2 + 1
      PRINT *

      !IF( ALLOCATED( pos_shells( r )% pos_shell ) ) &
      !   DEALLOCATE( pos_shells( r )% pos_shell )
      !  PRINT *, "0.1"
      !IF( ALLOCATED( pos_shells( r )% pvol_shell ) ) &
      !   DEALLOCATE( pos_shells( r )% pvol_shell )
      !          PRINT *, "0.2"
      !IF( ALLOCATED( pos_shells( r )% pvol_shell2 ) ) &
      !   DEALLOCATE( pos_shells( r )% pvol_shell2 )
      !          PRINT *, "0.3"
      !IF( ALLOCATED( pos_shells( r )% g_xx ) )&
      !   DEALLOCATE( pos_shells( r )% g_xx )
      !          PRINT *, "0.4"
      !IF( ALLOCATED( pos_shells( r )% baryon_density ) ) &
      !   DEALLOCATE( pos_shells( r )% baryon_density )
      !          PRINT *, "0.5"
      !IF( ALLOCATED( pos_shells( r )% gamma_euler ) ) &
      !   DEALLOCATE( pos_shells( r )% gamma_euler )
      !                  PRINT *, "0.6"
      !IF( ALLOCATED( pos_shells( r )% pos_th ) ) &
      !   DEALLOCATE( pos_shells( r )% pos_th )
      !                  PRINT *, "0.7"
      !IF( ALLOCATED( pos_shells( r )% pos_phi ) ) &
      !   DEALLOCATE( pos_shells( r )% pos_phi )

PRINT *, "Before allocation, shell ", r, "iteration ", cnt2 + 1
PRINT *

      !ALLOCATE( pos_shells( r )% pos_shell( 3, npart_shell( r ) ) )
      !ALLOCATE( pos_shells( r )% pvol_shell( npart_shell( r ) ) )
      !ALLOCATE( pos_shells( r )% pvol_shell2( npart_shell( r ) ) )
      !ALLOCATE( pos_shells( r )% g_xx( npart_shell( r ) ) )
      !ALLOCATE( pos_shells( r )% baryon_density( npart_shell( r ) ) )
      !ALLOCATE( pos_shells( r )% gamma_euler( npart_shell( r ) ) )
      !ALLOCATE( pos_shells( r )% pos_th( npart_shelleq( r ) ) )
      !ALLOCATE( pos_shells( r )% pos_phi( npart_shelleq( r ) ) )

      !phase= phase + r*alpha(r)/golden_ratio
      !CALL RANDOM_NUMBER( phase )
      !phase= phase*alpha(r)
      !dphi_shells= alpha(r)
      rad= shell_radii(r)
      itr= 1

      npart_shell_tmp= npart_shell( r )
      !PRINT *, npart_shell_tmp

      !CALL RANDOM_NUMBER( phase_th )
      !phase_th= phase_th*ABS( MINVAL( colatitude_pos(r)% colatitudes, DIM= 1 ) &
      !                     - pi/2.0D0 )*4.0D0/5.0D0

      PRINT *, "Right before OMP, shell ", r, "iteration ", cnt2 + 1
      PRINT *

      !$OMP PARALLEL DO DEFAULT(NONE), &
      !$OMP             PRIVATE(phase,col,xtemp,ytemp,ztemp, &
      !$OMP                     dphi_shells,dth_shells, &
      !$OMP                     th,phi,rand_num2,phase_th,rel_sign), &
      !$OMP             SHARED(r,npart_shelleq,center,rad,alpha,pos_shells, &
      !$OMP                    npart_shell,colatitude_pos,bns_obj, &
      !$OMP                    dr_shells), &
      !$OMP             REDUCTION(+:npart_out,itr)
      DO th= 1, npart_shelleq( r )/2, 1 !npart_shelleq( r ) is even, see above

        CALL RANDOM_NUMBER( phase )
        phase= phase*alpha(r)
        dphi_shells= alpha(r)

        col= colatitude_pos(r)% colatitudes(th)
        !IF( th == 1 )THEN
        !
        !  !dth_shells= pi - ( col + colatitude_pos(r)% colatitudes(th+1) )/2.0D0
        !  dth_shells= 2.0D0*( col - colatitude_pos(r)% colatitudes(th+1) )
        !
        !ELSEIF( th == npart_shelleq( r )/2 )THEN
        !
        !  !dth_shells= ( colatitude_pos(r)% colatitudes(th-1) + col - pi )/2.0D0
        !  dth_shells= 2.0D0*( colatitude_pos(r)% colatitudes(th-1) - col )
        !
        !ELSE
        !
        !  dth_shells= ( colatitude_pos(r)% colatitudes(th-1) + col )/2.0D0 &
        !            - ( col + colatitude_pos(r)% colatitudes(th+1) )/2.0D0
        !
        !ENDIF
        PRINT *, "Right before longitude loop"

        DO phi= 1, npart_shelleq( r ), 1

          CALL RANDOM_NUMBER( phase_th )
          CALL RANDOM_NUMBER( rand_num2 )
          IF( rand_num2 >= 0.5D0 ) rel_sign=  1
          IF( rand_num2 < 0.5D0 )  rel_sign= -1

          IF( col*( 1 + rel_sign*0.05D0*phase_th ) < pi .AND. &
              col*( 1 + rel_sign*0.05D0*phase_th ) > pi/2.0D0 )THEN

            col= col*( 1 + rel_sign*0.05D0*phase_th )

          ENDIF

        PRINT *, "2.1"

          xtemp= center + rad*COS(phase + phi*alpha(r))*SIN(col)
          ytemp= rad*SIN(phase + phi*alpha(r))*SIN(col)
          ztemp= rad*COS(col)

PRINT *, "2.2"

          CALL bns_obj% import_id( &
                   xtemp, ytemp, ztemp, &
                   pos_shells(r)% g_xx( itr ), &
                   pos_shells(r)% baryon_density( itr ), &
                   pos_shells(r)% gamma_euler( itr ) )

          IF( pos_shells(r)% baryon_density( itr ) > 0.0D0 &
              .AND. &
              bns_obj% is_hydro_negative( xtemp, ytemp, ztemp ) == 0 )THEN

            IF( pos_shells(r)% baryon_density( itr ) == 0 )THEN
              PRINT *, "When placing first half of particles"
              PRINT *, r, pos_shells(r)% pos_shell( 1, itr ), &
                       pos_shells(r)% pos_shell( 2, itr ), &
                       pos_shells(r)% pos_shell( 3, itr ), &
                       pos_shells(r)% baryon_density( itr )
            ENDIF

PRINT *, "2.3"

            npart_out= npart_out + 1
            pos_shells(r)% pos_shell( 1, itr )= xtemp
            pos_shells(r)% pos_shell( 2, itr )= ytemp
            pos_shells(r)% pos_shell( 3, itr )= ztemp
            !pos_shells(r)% pos_th( th )= col
            !pos_shells(r)% pos_phi( phi )= phase + phi*alpha(r)

         !   IF( th == 1 )THEN
         !
         !     !dth_shells= pi - ( col + colatitude_pos(r)% colatitudes(th+1) )/2.0D0
         !     dth_shells= 2.0D0*ABS( col - &
         !               ( col + colatitude_pos(r)% colatitudes(th + 1) )/2.0D0 )
         !
         !     !PRINT *, "case 1. ", "dt_shells= ", dth_shells, &
         !     !                     "col= ", col, &
         !     !                     "colatitude_pos(r)% colatitudes(th + 1)=", &
         !     !                     colatitude_pos(r)% colatitudes(th + 1)
         !     !PRINT *
         !
         !   ELSEIF( th == npart_shelleq( r )/2 )THEN
         !
         !     !dth_shells= ( colatitude_pos(r)% colatitudes(th-1) + col - pi )/2.0D0
         !     dth_shells= 2.0D0*ABS( ( colatitude_pos(r)% colatitudes(th - 1) &
         !               + col )/2.0D0 - col )
         !
         !     !PRINT *, "case 2. ", "dt_shells= ", dth_shells, &
         !     !                     "col= ", col, &
         !     !                     "colatitude_pos(r)% colatitudes(th - 1)=", &
         !     !                     colatitude_pos(r)% colatitudes(th - 1)
         !     !PRINT *
         !
         !   ELSE
         !
         !     dth_shells= ABS( &
         !             ( colatitude_pos(r)% colatitudes(th + 1) + col )/2.0D0 &
         !           - ( col + colatitude_pos(r)% colatitudes(th - 1) )/2.0D0 )
         !
         !
         !     !PRINT *, "case 3. ", "dt_shells= ", dth_shells, &
         !     !                     "col= ", col, &
         !     !                     "colatitude_pos(r)% colatitudes(th + 1)=", &
         !     !                     colatitude_pos(r)% colatitudes(th + 1), &
         !     !                     "colatitude_pos(r)% colatitudes(th - 1)=", &
         !     !                     colatitude_pos(r)% colatitudes(th - 1)
         !     !PRINT *
         !
         !   ENDIF
         !
         !   !pvol( npart_out )= rad**2.0D0*SIN(col) &
         !   !                   *dr_shells*dth_shells*dphi_shells
         !
         !   pos_shells(r)% pvol_shell2( itr )= rad**2.0D0*SIN(col) &
         !                      *dr_shells*dth_shells*dphi_shells! &
         !                      !*pos_shells(r)% g_xx( itr ) &
         !                      !*SQRT(pos_shells(r)% g_xx( itr ))
         !
         !   IF( pos_shells(r)% pvol_shell2( itr ) <= 0 )THEN
         !     PRINT *, "When placing first half of particles"
         !     PRINT *, "pos_shells(", r, ")% pvol_shell2( ", itr, " ) =", &
         !              pos_shells(r)% pvol_shell2( itr )
         !     STOP
         !   ENDIF

PRINT *, "2.4"

            itr= itr + 1

          ELSE

            npart_shell( r )= npart_shell( r ) - 2

          ENDIF

          ! Print progress on screen, every 10%
          !perc= 50*( THIS% nx*THIS% ny*iz + THIS% nx*iy + ix )/ &
          !        ( THIS% nx*THIS% ny*THIS% nz/2 )
          !IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
          !  WRITE( *, "(A2,I3,A1)", ADVANCE= "NO" ) &
          !         creturn//" ", perc, "%"
          !ENDIF

        ENDDO
        PRINT *, "Right after longitude loop"
      ENDDO
      !$OMP END PARALLEL DO
      PRINT *, "Right after OMP"
      IF( itr /= npart_shell( r )/2 + 1 )THEN
        PRINT *, "** ERROR! Mismatch in the particle counters on shell ", r
        PRINT *, "itr=", itr, "npart_shell( r )/2", npart_shell( r )/2
        PRINT *, "itr should be equal to npart_shell( r )/2 + 1. Stopping..."
        PRINT *
        STOP
      ENDIF

      !WRITE( *, "(A2,I4,I4,F3.2)", ADVANCE= "NO" ) &
      !          creturn//" ", npart_shell_tmp, npart_shell( r ), &
      !                        DBLE(npart_shell( r ))/DBLE(npart_shell_tmp)

      PRINT *, "Right before safety check"

      IF( npart_shell( r ) < 0 ) npart_shell( r )= 0
      IF( npart_shell( r ) == 0 )THEN
        m_parts( r )= m_parts( r - 1 )
        PRINT *, " * Placed", npart_shell( r ), &
                 " particles on spherical shell ", r, " of ", n_shells
        IF( r == n_shells )THEN
          EXIT
        ELSE
          PRINT *, "r=", r
          m_parts( r )= m_parts( r - 1 )
          r= r + 1
          cnt2 = 0
          upper_bound_tmp= upper_bound
          lower_bound_tmp= lower_bound
          CYCLE
        ENDIF
      ELSE
        m_parts( r )= shell_masses( r )/DBLE(npart_shell( r ))
      ENDIF

      PRINT *, "Right before correction of particle number"

      IF( r > 1 )THEN

        high_mass= m_parts( r )/m_parts( r - 1 ) > upper_bound_tmp
        low_mass = m_parts( r )/m_parts( r - 1 ) < lower_bound_tmp
        kept_all = npart_shell_kept == 1.0D0
        npart_shell_kept= DBLE(npart_shell( r ))/DBLE(npart_shell_tmp)

        PRINT *, "cnt2=", cnt2
        PRINT *, "upper_bound_tmp=", upper_bound_tmp
        PRINT *, "lower_bound_tmp=", lower_bound_tmp
        PRINT *, "n_shells=", n_shells
        PRINT *, "r=", r
        PRINT *, "npart_shell( r )=", npart_shell( r )
        PRINT *, "npart_shell_tmp=", npart_shell_tmp
        PRINT *, "npart_shell_kept=", npart_shell_kept
        PRINT *, "high_mass=", high_mass
        PRINT *, "low_mass=", low_mass
        PRINT *, "kept_all=", kept_all
        PRINT *, "m_parts( r )=", m_parts( r )
        PRINT *, "m_parts( r - 1 )=", m_parts( r - 1 )
        PRINT *, " m_parts( r )/m_parts( r - 1 )= ",  &
                                   m_parts( r )/m_parts( r - 1 )
        PRINT *

        IF( high_mass .AND. kept_all )THEN
PRINT *, "case 1"

          cnt2= cnt2 + 1
PRINT *, "1.2"
          IF( cnt2 > 100 )THEN
            upper_bound_tmp= upper_bound_tmp*1.01D0
            lower_bound_tmp= lower_bound_tmp*0.99D0
            cnt2= 1
          ENDIF
PRINT *, "1.4"
          npart_out= npart_out - ( itr - 1 )
PRINT *, "1.6"
          CALL RANDOM_NUMBER( rand_num )
          CALL RANDOM_NUMBER( rand_num2 )
          npart_shelleq( r )= npart_shelleq( r ) + 1*NINT( 1 + 1.0*rand_num ) &
                                                 + 1*NINT( 1 + 1.0*rand_num2 )
PRINT *, "1.8"
          IF( npart_shelleq( r ) == 0 .OR. npart_shell( r ) == 0 )THEN
            CALL RANDOM_NUMBER( rand_num )
            CALL RANDOM_NUMBER( rand_num2 )
            IF( rand_num2 < 0.5D0 )  rel_sign= - 1
            IF( rand_num2 >= 0.5D0 ) rel_sign=   1
            npart_shelleq( r )= npart_shelleq( r - 1 ) &
                              + rel_sign*NINT( 1 + rand_num )
          ENDIF
PRINT *, "1.9"
          CYCLE

        ELSEIF( low_mass .AND. kept_all )THEN
        PRINT *, "case 2"

          cnt2= cnt2 + 1
          IF( cnt2 > 100 )THEN
            upper_bound_tmp= upper_bound_tmp*1.01D0
            lower_bound_tmp= lower_bound_tmp*0.99D0
            cnt2= 1
          ENDIF

          npart_out= npart_out - ( itr - 1 )

          CALL RANDOM_NUMBER( rand_num )
          CALL RANDOM_NUMBER( rand_num2 )
          npart_shelleq( r )= npart_shelleq( r ) - 1*NINT( 1 + 1.0*rand_num ) &
                                                 - 1*NINT( 1 + 1.0*rand_num2 )

          IF( npart_shelleq( r ) == 0 .OR. npart_shell( r ) == 0 )THEN
            CALL RANDOM_NUMBER( rand_num )
            CALL RANDOM_NUMBER( rand_num2 )
            IF( rand_num2 < 0.5D0 )  rel_sign= - 1
            IF( rand_num2 >= 0.5D0 ) rel_sign=   1
            npart_shelleq( r )= npart_shelleq( r - 1 ) &
                              + rel_sign*NINT( 1 + rand_num )
          ENDIF

          CYCLE

        ELSEIF( high_mass .AND. .NOT.kept_all ) THEN
        PRINT *, "case 3"

          cnt2= cnt2 + 1
          IF( cnt2 > 100 )THEN
            upper_bound_tmp= upper_bound_tmp*1.01D0
            lower_bound_tmp= lower_bound_tmp*0.99D0
            cnt2= 1
          ENDIF

          npart_out= npart_out - ( itr - 1 )

          CALL RANDOM_NUMBER( rand_num )
          CALL RANDOM_NUMBER( rand_num2 )
          IF( rand_num2 < 0.5D0 )  rel_sign= - 1
          IF( rand_num2 >= 0.5D0 ) rel_sign=   1
          npart_shelleq( r )= CEILING( SQRT( &
                                (shell_masses( r )/m_parts( r - 1 )) &
                                /npart_shell_kept &
                              ) ) + rel_sign*NINT( 1 + rand_num )

          IF( npart_shelleq( r ) == 0 .OR. npart_shell( r ) == 0 )THEN
            CALL RANDOM_NUMBER( rand_num )
            CALL RANDOM_NUMBER( rand_num2 )
            IF( rand_num2 < 0.5D0 )  rel_sign= - 1
            IF( rand_num2 >= 0.5D0 ) rel_sign=   1
            npart_shelleq( r )= npart_shelleq( r - 1 ) &
                              + rel_sign*NINT( 1 + rand_num )
          ENDIF

          CYCLE

        ELSEIF( low_mass .AND. .NOT.kept_all ) THEN
        PRINT *, "case 4"

          cnt2= cnt2 + 1
          IF( cnt2 > 100 )THEN
            upper_bound_tmp= upper_bound_tmp*1.01D0
            lower_bound_tmp= lower_bound_tmp*0.99D0
            cnt2= 1
          ENDIF

          npart_out= npart_out - ( itr - 1 )

          CALL RANDOM_NUMBER( rand_num )
          CALL RANDOM_NUMBER( rand_num2 )
          IF( rand_num2 < 0.5D0 )  rel_sign= - 1
          IF( rand_num2 >= 0.5D0 ) rel_sign=   1
          npart_shelleq( r )= CEILING( SQRT( &
                                (shell_masses( r )/m_parts( r - 1 )) &
                                /npart_shell_kept &
                              ) ) + rel_sign*NINT( 1 + rand_num )

          IF( npart_shelleq( r ) == 0 .OR. npart_shell( r ) == 0 )THEN
            CALL RANDOM_NUMBER( rand_num )
            CALL RANDOM_NUMBER( rand_num2 )
            IF( rand_num2 < 0.5D0 )  rel_sign= - 1
            IF( rand_num2 >= 0.5D0 ) rel_sign=   1
            npart_shelleq( r )= npart_shelleq( r - 1 ) &
                              + rel_sign*NINT( 1 + rand_num )
          ENDIF

          PRINT *, "Right after correction of particle number"

          CYCLE

        ENDIF

      ENDIF

      PRINT *, " * Placed", npart_shell( r ), &
               " particles on spherical shell ", r, " of ", n_shells
      IF( r > 1 ) PRINT *, "   Ratio of particle masses on last 2 shells: ", &
                           "   m_parts(", r, ")/m_parts(", r - 1, ")= ",  &
                           m_parts( r )/m_parts( r - 1 )
      PRINT *

      IF( r == n_shells ) EXIT

      !IF( r == 2 ) STOP

      r= r + 1
      cnt2 = 0
      upper_bound_tmp= upper_bound
      lower_bound_tmp= lower_bound

      !pos_shells( r )% pos_shell= &
      !                pos_shells( r )% pos_shell( :, 1:npart_shell( r ) )
      !pos_shells( r )% pvol_shell= &
      !                pos_shells( r )% pvol_shell( 1:npart_shell( r ) )
      !pos_shells( r )% pvol_shell2= &
      !                pos_shells( r )% pvol_shell2( 1:npart_shell( r ) )
      !pos_shells( r )% g_xx= &
      !                pos_shells( r )% g_xx( 1:npart_shell( r ) )
      !pos_shells( r )% baryon_density= &
      !                pos_shells( r )% baryon_density( 1:npart_shell( r ) )
      !pos_shells( r )% gamma_euler= &
      !                pos_shells( r )% gamma_euler( 1:npart_shell( r ) )
      !pos_shells( r )% pos_th= &
      !                pos_shells( r )% pos_th( 1:npart_shell( r ) )
      !pos_shells( r )% pos_phi= &
      !                pos_shells( r )% pos_phi( 1:npart_shell( r ) )

    ENDDO
    PRINT *, "npart=", npart_out
    PRINT *, npart_shell
    PRINT *
    PRINT *, m_parts
    PRINT *
    PRINT *, "particle mass ratio= ", MAXVAL(m_parts)/MINVAL(m_parts)
    npart_test= 0
    DO r= 1, n_shells, 1
      npart_test= npart_test + npart_shell( r )
    ENDDO
    IF( npart_test/2 /= npart_out )THEN
      PRINT *, "** ERROR! The sum of the particles on the shells is not ", &
               "equal to the total number of particles. Stopping.."
      PRINT *, "npart_test", npart_test/2, ", npart_out=", npart_out
      PRINT *
      STOP
    ENDIF
    STOP

    PRINT *, "Mirroring particles..."

    !PRINT *, "npart/2=", npart_out
    !npart_half= npart_out
    !DO itr= 1, npart_half, 1
    !  npart_out= npart_out + 1
    !  pos_shells(r)% pos_shell( 1, npart_out )=   pos_shells(r)% pos_shell( 1, npart_out )
    !  pos_shells(r)% pos_shell( 2, npart_out )=   pos_shells(r)% pos_shell( 2, npart_out )
    !  pos_shells(r)% pos_shell( 3, npart_out )= - pos_shells(r)% pos_shell( 3, npart_out )
    !  !pvol( npart_out )  =   pvol( itr )
    !ENDDO
    !PRINT *, "npart=", npart_out
    !IF( SUM( npart_shell, DIM=1 ) /= npart_out )THEN
    !  PRINT *, "** ERROR! The sum of the particles on the shells is not ", &
    !           "equal to the total number of particles. Stopping.."
    !  PRINT *, "SUM( npart_shell )", SUM( npart_shell, DIM=1 ), &
    !           ", npart_out=", npart_out
    !  PRINT *
    !  STOP
    !ENDIF

    PRINT *, "npart/2=", npart_out
    !npart_half= npart_out
    DO r= 1, n_shells, 1
      DO itr= 1, npart_shell( r )/2, 1
        npart_out= npart_out + 1
        pos_shells(r)% pos_shell( 1, npart_shell( r )/2 + itr )= &
                                          pos_shells(r)% pos_shell( 1, itr )
        pos_shells(r)% pos_shell( 2, npart_shell( r )/2 + itr )= &
                                          pos_shells(r)% pos_shell( 2, itr )
        pos_shells(r)% pos_shell( 3, npart_shell( r )/2 + itr )= &
                                        - pos_shells(r)% pos_shell( 3, itr )
        pos_shells(r)% g_xx( npart_shell( r )/2 + itr )= &
                                                  pos_shells(r)% g_xx( itr )
        pos_shells(r)% baryon_density( npart_shell( r )/2 + itr )= &
                                        pos_shells(r)% baryon_density( itr )
        pos_shells(r)% gamma_euler( npart_shell( r )/2 + itr )= &
                                          pos_shells(r)% gamma_euler( itr )
        !pos_shells(r)% pos_th( npart_shell( r )/2 + itr )= &
        !                                  pi - pos_shells(r)% pos_th( itr )
        !pos_shells(r)% pos_phi( npart_shell( r )/2 + itr )= &
        !                                  pos_shells(r)% pos_phi( itr )
        pos_shells(r)% pvol_shell2( npart_shell( r )/2 + itr )= &
                                           pos_shells(r)% pvol_shell2( itr )
        !pvol( npart_out )  =   pvol( itr )
        IF( pos_shells(r)% baryon_density( itr ) == 0 )THEN
          PRINT *, "When mirroring particles"
          PRINT *, r, pos_shells(r)% pos_shell( 1, itr ), &
                   pos_shells(r)% pos_shell( 2, itr ), &
                   pos_shells(r)% pos_shell( 3, itr ), &
                   pos_shells(r)% baryon_density( itr )
          STOP
        ENDIF
        !IF( pos_shells(r)% pvol_shell2( itr ) < 0 )THEN
        !  PRINT *, "When mirroring particles"
        !  PRINT *, "pos_shells(", r, ")% pvol_shell2( ", itr, " ) =", &
        !           pos_shells(r)% pvol_shell2( itr )
        !  STOP
        !ENDIF
      ENDDO
    ENDDO
    PRINT *, "npart=", npart_out
    npart_test= 0
    DO r= 1, n_shells, 1
      npart_test= npart_test + npart_shell( r )
    ENDDO
    IF( npart_test /= npart_out )THEN
      PRINT *, "** ERROR! The sum of the particles on the shells is not ", &
               "equal to the total number of particles. Stopping.."
      PRINT *, "npart_test", npart_test, ", npart_out=", npart_out
      PRINT *
      STOP
    ENDIF
    IF( SUM( npart_shell, DIM=1 ) /= npart_out )THEN
      PRINT *, "** ERROR! The sum of the particles on the shells is not ", &
               "equal to the total number of particles. Stopping.."
      PRINT *, "SUM( npart_shell )", SUM( npart_shell, DIM=1 ), &
               ", npart_out=", npart_out
      PRINT *
      STOP
    ENDIF

  ! DO r= 1, n_shells, 1
  !
  !   dphi_shells= alpha(r)
  !   rad= shell_radii(r)
  !
  !   DO th= 1, npart_shelleq( r )/2, 1
  !
  !     col= pos_shells(r)% pos_th( th )
  !     IF( th == 1 )THEN
  !
  !       !dth_shells= pi - ( col + colatitude_pos(r)% colatitudes(th+1) )/2.0D0
  !       dth_shells= 2.0D0*( col - pos_shells(r)% pos_th( th + 1 ) )
  !
  !     ELSEIF( th == npart_shelleq( r )/2 )THEN
  !
  !       !dth_shells= ( colatitude_pos(r)% colatitudes(th-1) + col - pi )/2.0D0
  !       dth_shells= 2.0D0*( pos_shells(r)% pos_th( th - 1 ) - col )
  !
  !     ELSE
  !
  !       dth_shells= ( pos_shells(r)% pos_th( th + 1 ) + col )/2.0D0 &
  !                 - ( col + pos_shells(r)% pos_th( th - 1 ) )/2.0D0
  !
  !     ENDIF
  !
  !     DO phi= 1, npart_shelleq( r ), 1
  !
  !       PRINT *, "rad= ", rad, "col= ", col
  !       PRINT *, "dr_shells= ", dr_shells, "dth_shells= ", dth_shells, &
  !                "dphi_shells= ", dphi_shells
  !
  !       pos_shells(r)% pvol_shell2( itr )= rad**2.0D0*SIN(col) &
  !                          *dr_shells*dth_shells*dphi_shells! &
  !                          !*pos_shells(r)% g_xx( itr ) &
  !                          !*SQRT(pos_shells(r)% g_xx( itr ))
  !       itr= itr + 1
  !
  !     ENDDO
  !   ENDDO
  !   !DO itr= 1, npart_shell( r ), 1
  !   !
  !   !  !PRINT *, pos_shells(r)% pos_th
  !   !
  !   !  col= pos_shells(r)% pos_th( itr )
  !   !  IF( col == MAXVAL( pos_shells(r)% pos_th ) )THEN
  !   !
  !   !    !dth_shells= pi - ( col + colatitude_pos(r)% colatitudes(th+1) )/2.0D0
  !   !    dth_shells= 2.0D0*( col - pos_shells(r)% pos_th( itr + 1 ) )
  !   !
  !   !  ELSEIF( col == MINVAL( pos_shells(r)% pos_th ) )THEN
  !   !
  !   !    !dth_shells= ( colatitude_pos(r)% colatitudes(th-1) + col - pi )/2.0D0
  !   !    dth_shells= 2.0D0*( pos_shells(r)% pos_th( npart_shelleq( r ) - 1 ) - col )
  !   !
  !   !  ELSE
  !   !
  !   !    dth_shells= ( pos_shells(r)% pos_th( itr + 1 ) + col )/2.0D0 &
  !   !              - ( col + pos_shells(r)% pos_th( itr - 1 ) )/2.0D0
  !   !
  !   !  ENDIF
  !   !
  !   !  IF( .TRUE. )THEN
  !   !    PRINT *, "rad= ", rad, "col= ", col, &
  !   !  !"pos_shells(r)% pos_th( itr - 1 )=", pos_shells(r)% pos_th( itr - 1 ), &
  !   !  "pos_shells(r)% pos_th( itr + 1 )=", pos_shells(r)% pos_th( npart_shelleq( r ) + 1 )
  !   !    PRINT *, "dr_shells= ", dr_shells, "dth_shells= ", dth_shells, &
  !   !             "dphi_shells= ", dphi_shells
  !   !  ENDIF
  !   !
  !   !  pos_shells(r)% pvol_shell2( itr )= rad**2.0D0*SIN(col) &
  !   !                     *dr_shells*dth_shells*dphi_shells! &
  !   !                     !*pos_shells(r)% g_xx( itr ) &
  !   !                     !*SQRT(pos_shells(r)% g_xx( itr ))
  !   !  IF( pos_shells(r)% pvol_shell2( itr ) <= 0.0D0 )THEN
  !   !    PRINT *, "When computing particle volume"
  !   !    PRINT *, "pos_shells(", r, ")% pvol_shell2( ", itr, " ) =", &
  !   !             pos_shells(r)% pvol_shell2( itr )
  !   !    STOP
  !   !  ENDIF
  !   !ENDDO
  !   !PRINT *
  ! ENDDO
    !DO r= 1, n_shells, 1
    !  DO itr= 1, npart_shell( r )/2, 1
    !    pos_shells(r)% pvol_shell2( npart_shell( r )/2 + itr )= &
    !                                       pos_shells(r)% pvol_shell2( itr )
    !  ENDDO
    !ENDDO
  !  PRINT *
  !  STOP

    mass_test= 0.0D0
    mass_test2= 0.0D0
    proper_volume_test= 0.0D0
    proper_volume= 0.0D0
    i_shell= 1
    DO r= 1, n_shells, 1
      !DO itr= i_shell, (i_shell - 1) + (npart_shelleq(r)**2.0D0)/2.0D0, 1
      !  CALL bns_obj% import_id( &
      !           pos( 1, itr ), pos( 2, itr ), pos( 3, itr ), &
      !           g_xx, baryon_density, gamma_euler )
      !
      !  pvol( itr )= m_parts( r )/( baryon_density*g_xx*SQRT(g_xx)*gamma_euler )
      !
      !  proper_volume_test= proper_volume_test + 2.0D0*pvol( itr )*g_xx*SQRT(g_xx)
      !  mass_test= mass_test + &
      !              2.0D0*baryon_density*pvol( itr )*g_xx*SQRT(g_xx)*gamma_euler
      !ENDDO
      !i_shell= i_shell + (npart_shelleq(r)**2.0D0)/2.0D0
      DO itr= 1, npart_shell( r ), 1

        IF( pos_shells(r)% baryon_density( itr ) == 0 )THEN
          PRINT *, "When computing particle volume"
          PRINT *, r, itr, pos_shells(r)% pos_shell( 1, itr ), &
                   pos_shells(r)% pos_shell( 2, itr ), &
                   pos_shells(r)% pos_shell( 3, itr ), &
                   pos_shells(r)% baryon_density( itr )
        ENDIF
        IF( pos_shells(r)% pvol_shell2( itr ) <= 0.0D0 )THEN
          PRINT *, "When computing particle volume"
          PRINT *, "pos_shells(", r, ")% pvol_shell2( ", itr, " ) =", &
                   pos_shells(r)% pvol_shell2( itr )
          STOP
        ENDIF

        pos_shells(r)% pvol_shell( itr )= m_parts( r ) &
                          /( pos_shells(r)% baryon_density( itr ) &
                            *pos_shells(r)% g_xx( itr ) &
                            *SQRT(pos_shells(r)% g_xx( itr )) &
                            *pos_shells(r)% gamma_euler( itr ) )

        proper_volume_test= proper_volume_test + &
                            pos_shells(r)% pvol_shell( itr )! &
                            !*pos_shells(r)% g_xx( itr ) &
                            !*SQRT(pos_shells(r)% g_xx( itr ))

        proper_volume= proper_volume + pos_shells(r)% pvol_shell2( itr )
        !PRINT *, proper_volume

        mass_test= mass_test + pos_shells(r)% baryon_density( itr ) &
                  *pos_shells(r)% pvol_shell( itr ) &
                  *pos_shells(r)% g_xx( itr ) &
                  *SQRT(pos_shells(r)% g_xx( itr )) &
                  *pos_shells(r)% gamma_euler( itr )
        mass_test2= mass_test2 + pos_shells(r)% baryon_density( itr ) &
                  *pos_shells(r)% pvol_shell2( itr ) &
                  *pos_shells(r)% g_xx( itr ) &
                  *SQRT(pos_shells(r)% g_xx( itr )) &
                  *pos_shells(r)% gamma_euler( itr )

      ENDDO
    ENDDO

    PRINT *, mass_test, mass_test2, mass_star, proper_volume_test, proper_volume
    PRINT *

    mass_test = 0.0D0
    mass_test2= 0.0D0
    proper_volume_test = 0.0D0
    proper_volume= 0.0D0
    vol_shell( r )  = 0.0D0
    vol_shell2( r ) = 0.0D0
    mass_shell( r ) = 0.0D0
    mass_shell2( r )= 0.0D0
    DO r= 1, n_shells, 1
      DO itr= 1, npart_shell( r ), 1
        IF( pos_shells(r)% pvol_shell2( itr ) <= 0 )THEN
          PRINT *, "When computing shell volumes and masses"
          PRINT *, "pos_shells(", r, ")% pvol_shell2( ", itr, " ) =", &
                   pos_shells(r)% pvol_shell2( itr )
          STOP
        ENDIF
        vol_shell( r )  = vol_shell( r )  + pos_shells(r)% pvol_shell( itr )
        vol_shell2( r ) = vol_shell2( r ) + pos_shells(r)% pvol_shell2( itr )
        mass_shell( r ) = mass_shell( r ) + &
                          pos_shells(r)% baryon_density( itr ) &
                         *pos_shells(r)% pvol_shell( itr ) &
                         *pos_shells(r)% g_xx( itr ) &
                         *SQRT(pos_shells(r)% g_xx( itr )) &
                         *pos_shells(r)% gamma_euler( itr )
        mass_shell2( r )= mass_shell2( r ) + &
                          pos_shells(r)% baryon_density( itr ) &
                         *pos_shells(r)% pvol_shell2( itr ) &
                         *pos_shells(r)% g_xx( itr ) &
                         *SQRT(pos_shells(r)% g_xx( itr )) &
                         *pos_shells(r)% gamma_euler( itr )
      ENDDO
      mass_test= mass_test + mass_shell( r )
      mass_test2= mass_test2 + mass_shell2( r )
      proper_volume= proper_volume + vol_shell( r )
      proper_volume_test= proper_volume_test + vol_shell2( r )
      IF( r > 1 )THEN
        PRINT *, "shell", r
        PRINT *, "  shell volumes:", vol_shell( r ), vol_shell2( r ), &
                 4.0D0/3.0D0*pi* &
                 ( shell_radii( r )**3.0D0 - shell_radii( r - 1 )**3.0D0 )
        PRINT *, "  shell masses:", mass_shell( r ), mass_shell2( r ), &
                 shell_masses( r )
        PRINT *
      ENDIF
    ENDDO
    PRINT *
    PRINT *, "masses of the star:", mass_test, mass_test2, mass_star
    PRINT *, "volumes of the star:", proper_volume, proper_volume_test, &
                                     4.0D0/3.0D0*pi*radius**3.0D0
    PRINT *
    !STOP

    IF(.NOT.ALLOCATED( pos ))THEN
      ALLOCATE( pos( 3, npart_out ), &
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
    IF(.NOT.ALLOCATED( pvol ))THEN
      ALLOCATE( pvol( npart_out ), &
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
    IF(.NOT.ALLOCATED( pmass ))THEN
      ALLOCATE( pmass( npart_out ), &
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

    cnt= 0
    npart_shell_tmp= 0
    DO r= 1, n_shells, 1
      DO itr= 1, npart_shell( r ), 1
        pos( 1, itr + npart_shell_tmp )= pos_shells(r)% pos_shell( 1, itr )
        pos( 2, itr + npart_shell_tmp )= pos_shells(r)% pos_shell( 2, itr )
        pos( 3, itr + npart_shell_tmp )= pos_shells(r)% pos_shell( 3, itr )
        pvol( itr + npart_shell_tmp )  = pos_shells(r)% pvol_shell( itr )
        pmass( itr + npart_shell_tmp ) = m_parts( r )
        cnt= cnt + 1
      ENDDO
      npart_shell_tmp= cnt
    ENDDO
    IF( cnt /= npart_out )THEN
      PRINT *, "** ERROR! The sum of the particles on the shells is not ", &
               "equal to the total number of particles. Stopping.."
      PRINT *, "cnt", cnt, ", npart_out=", npart_out
      PRINT *
      STOP
    ENDIF



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

    !DO itr = 1, npart_out, 1
    !
    !  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    !    pos( 1, itr ), pos( 2, itr ), pos( 3, itr )
    !
    !  IF( ios > 0 )THEN
    !    PRINT *, "...error when writing the arrays in " &
    !             // TRIM(finalnamefile), ". The error message is", err_msg
    !    STOP
    !  ENDIF
    !
    !ENDDO

    DO r= 1, n_shells, 1
      DO itr= 1, npart_shell( r ), 1

        WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
          r, pos_shells(r)% pos_shell( 1, itr ), &
          pos_shells(r)% pos_shell( 2, itr ), &
          pos_shells(r)% pos_shell( 3, itr )

        IF( ios > 0 )THEN
          PRINT *, "...error when writing the arrays in " &
                   // TRIM(finalnamefile), ". The error message is", err_msg
          STOP
        ENDIF

      ENDDO
    ENDDO

    CLOSE( UNIT= 2 )

    PRINT *, " * Particles placed."

    STOP

  END PROCEDURE place_particles_stretched_lattice


END SUBMODULE stretched_lattice
