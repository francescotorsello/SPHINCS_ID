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
                         golden_ratio
    USE matrix,    ONLY: determinant_4x4_matrix
    USE NR,        ONLY: indexx

    IMPLICIT NONE

    INTEGER:: n_shells, itr2, itr3, mass_index, npart_half, npart, &
              shell_index, r, th, phi
    INTEGER, DIMENSION(:), ALLOCATABLE:: mass_profile_idx
    INTEGER, DIMENSION(:), ALLOCATABLE:: npart_shell, npart_shelleq

    DOUBLE PRECISION:: m_p, central_density, &
                       dr, dth, dphi, phase, mass
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: mass_profile
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: shell_radii, shell_masses, &
                                                  alpha

    LOGICAL:: exist

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    TYPE:: colatitude_pos_shell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: colatitudes
    END TYPE

    TYPE(colatitude_pos_shell), DIMENSION(:), ALLOCATABLE:: colatitude_pos

    !npart_approx= 1D+6
    !npart2_tmp= 1D+5*1.8D0/1.2D0

    !THIS% mass1= bns_obj% get_mass1()
    !THIS% mass2= bns_obj% get_mass2()
    !radius    = bns_obj% get_radius1_x_comp()
    !radius2    = bns_obj% get_radius2_x_comp()

    m_p= THIS% mass1/npart_approx

    n_shells= NINT( radius* &
                  (npart_approx/(4.0D0/3.0D0*pi*radius**3.0D0))**(1.0D0/3.0D0) )
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

    PRINT *, " ** Integrating the baryon mass density to get the mass profile..."
    PRINT *

    dr             = radius/500.0D0
    dth            = pi/2.0D0/250.0D0
    dphi           = 2.0D0*pi/500.0D0
    central_density= bns_obj% get_rho_center1()
    CALL bns_obj% integrate_baryon_mass_density( center, radius, &
                                                 central_density, &
                                                 dr, dth, dphi, &
                                                 mass, mass_profile, &
                                                 mass_profile_idx )

    DO itr= 1, n_shells, 1

      shell_radii( itr )= radius*itr/n_shells

    ENDDO
    shell_index= 1
    itr2= 0
    DO itr= 0, NINT(radius/dr), 1

      IF( mass_profile( 1, mass_profile_idx(itr) ) &
          >= shell_radii( shell_index ) )THEN

        shell_masses( shell_index )= SUM( mass_profile( 2, &
                      mass_profile_idx(itr2):mass_profile_idx(itr) ), DIM= 1 )

        IF( shell_index == n_shells )THEN
          EXIT
        ELSE
         itr2= itr + 1
         shell_index= shell_index + 1
        ENDIF

      ENDIF

    ENDDO

    ALLOCATE( npart_shell(n_shells) )
    ALLOCATE( npart_shelleq(n_shells) )
    ALLOCATE( alpha(n_shells) )
    ALLOCATE( colatitude_pos( n_shells ) )

    DO itr= 1, n_shells, 1

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
    npart= SUM( npart_shell, DIM= 1 )

    !PRINT *, colatitude_pos( 1 )% colatitudes
    !STOP

    ! Place the particles for one star only, since the subroutine will place
    ! particles for one star

    ! Allocating the memory for the array pos( 3, npart_temp )
    ! Note that after determining npart, the array pos is reshaped into
    ! pos( 3, npart )
    IF(.NOT.ALLOCATED( THIS% pos ))THEN
      ALLOCATE( THIS% pos( 3, npart ), &
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

    PRINT *, SUM( npart_shell, DIM= 1 )

    npart_out= 0
    THIS% pos= 0
    phase= 0
    DO r= 1, n_shells, 1
      !phase= phase + (r - 1.0D0)*2*pi/(SQRT(2.0D0)/(7.0D0/5.0D0)*nradii1_plane)
      phase= phase + r*alpha(r)/golden_ratio
      DO th= 1, npart_shelleq( r )/2, 1
        DO phi= 1, npart_shelleq( r ), 1 !part_shelleq( r ) is even, see above

          npart_out= npart_out + 1
          THIS% pos( 1, npart_out )= &
                center + shell_radii(r)*COS(phase + phi*alpha(r)) &
                              *SIN(colatitude_pos(r)% colatitudes(th))
          THIS% pos( 2, npart_out )= &
                shell_radii(r)*SIN(phase + phi*alpha(r)) &
                              *SIN(colatitude_pos(r)% colatitudes(th))
          THIS% pos( 3, npart_out )= &
                shell_radii(r)*COS(colatitude_pos(r)% colatitudes(th))

        ENDDO
      ENDDO
    ENDDO
    PRINT *, "npart=", npart_out

    PRINT *, "Mirroring particles..."

    PRINT *, "npart/2=", npart_out
    npart_half= npart_out
    DO itr= 1, npart_half, 1
      npart_out= npart_out + 1
      THIS% pos( 1, npart_out )=   THIS% pos( 1, itr )
      THIS% pos( 2, npart_out )=   THIS% pos( 2, itr )
      THIS% pos( 3, npart_out )= - THIS% pos( 3, itr )
    ENDDO
    PRINT *, "npart=", npart_out

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

    DO itr = 1, npart_out, 1

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

  END PROCEDURE place_particles_stretched_lattice


END SUBMODULE stretched_lattice
