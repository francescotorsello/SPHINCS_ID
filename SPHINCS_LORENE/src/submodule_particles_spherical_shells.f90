! File:         submodule_particles_constructor.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (particles_id) spherical_shells

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


  MODULE PROCEDURE place_particles_spherical_shells

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
                         golden_ratio, third, half, amu, g2kg, sixth
    USE matrix,    ONLY: determinant_4x4_matrix
    USE NR,        ONLY: indexx

    IMPLICIT NONE

    INTEGER:: n_shells, itr2, itr3, mass_index, npart_half, npart_tmp, cnt, &
              shell_index, r, th, phi, i_shell, npart_test, npart_shell_tmp, &
              cnt2, rel_sign, cnt3, dim_seed, r_cnt, first_shell, prev_shell, &
              npart_discard
    !INTEGER, PARAMETER:: max_length= 5D+6
    INTEGER, DIMENSION(:), ALLOCATABLE:: mass_profile_idx, seed
    INTEGER, DIMENSION(:), ALLOCATABLE:: npart_shell, npart_shelleq

    DOUBLE PRECISION:: xtemp, ytemp, ztemp, m_p, central_density, &
                       dr, dth, dphi, phase, phase_th, mass, baryon_density, &
                       dr_shells, dth_shells, dphi_shells, col, rad, &
                       g_xx, gamma_euler, proper_volume, mass_test, mass_test2,&
                       proper_volume_test, npart_shell_kept, &
                       rand_num, rand_num2, delta_r, shell_thickness, &
                       upper_bound_tmp, lower_bound_tmp, col_tmp, &
                       surface_density, density_step, n_shells_tmp, &
                       gxx_tmp, baryon_density_tmp, gamma_euler_tmp, rho_tmp
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: mass_profile, &
                                                    particle_profile
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: shell_radii, shell_masses, &
                                                  alpha, m_parts, vol_shell, &
                                                  vol_shell2, mass_shell, &
                                                  mass_shell2, shell_scales

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

    DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE:: pos_shell_tmp
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: g_xx_tmp
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: bar_density_tmp
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: gam_euler_tmp
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pvol_tmp

    LOGICAL, PARAMETER:: debug= .FALSE.

    PRINT *, "** Executing the place_particles_shells..."
    PRINT *

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

    !n_shells= NINT( radius* &
    !              (npart_approx/(4.0D0/3.0D0*pi*radius**3.0D0))**third )
    !PRINT *, n_shells

    n_shells_tmp= 0.0D0
    particle_profile= 0.0D0
    shell_index= 1
    IF(.NOT.ALLOCATED( particle_profile ))THEN
      ALLOCATE( particle_profile( 2, 500 ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array particle_profile in" &
                  // "SUBROUTINE place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array pos in SUBROUTINE" &
      !                // "place_particles_3D_lattice." )
    ENDIF
    DO r= 1, 500, 1

      n_shells_tmp= n_shells_tmp + &
                      radius/500*( ( bns_obj% import_mass_density( &
                                     center + r*radius/500, 0.0D0, 0.0D0 ) &
                                     )/m_p )**third
      particle_profile( 1, r )= r*radius/500
      particle_profile( 2, r )= n_shells_tmp

    ENDDO
    n_shells= NINT( n_shells_tmp )
    !PRINT *, n_shells_tmp, n_shells, shell_index
    !STOP

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
    IF(.NOT.ALLOCATED( shell_scales ))THEN
      ALLOCATE( shell_scales( n_shells ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array shell_scales in SUBROUTINE" &
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
      ALLOCATE( m_parts( 0:n_shells ), STAT= ios, &
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

    central_density= bns_obj% get_rho_center1()

    surface_density= bns_obj% import_mass_density( center + radius, &
                                                   0.0D0, 0.0D0 )
    density_step= ( central_density - surface_density )/(n_shells)
    shell_radii= 0.0D0

    !--------------------------------------------------------------!
    !-- Place shells based on constant intervals of mass density --!
    !-- this does the opposite of what I want....
    !--------------------------------------------------------------!

    !shell_index= 2
    !DO itr= 0, 200, 1
    !
    !  IF( bns_obj% import_mass_density( center + itr*radius/200, 0.0D0, 0.0D0 )&
    !    <= central_density - ( shell_index - 1 )*density_step &
    !  )THEN
    !
    !    shell_radii( shell_index )= itr*radius/200
    !    shell_index= shell_index + 1
    !
    !    !shell_radii( itr + 1 )= &
    !    !  ( ( central_density - itr*density_step)*Msun/amu )**(-third)
    !
    !  ENDIF
    !
    !ENDDO
    !shell_radii( 1 )= shell_radii( 2 )/2.0D0
    !PRINT *, n_shells, central_density, surface_density, density_step
    !PRINT *, radius
    !PRINT *, shell_radii
    !STOP

    !-----------------------------------------------------!
    !-- Place shells based on mass density a that point --!
    !-----------------------------------------------------!

!    DO itr= 0, n_shells - 1, 1
!
!      shell_scales( itr + 1 )= &
!          ( ( central_density - itr*density_step )/m_p )**(-third)
!
!    ENDDO
!    shell_scales= shell_scales - shell_scales( 1 )
!    shell_scales= shell_scales/shell_scales( n_shells )
!    PRINT *, shell_scales
!    PRINT *, SUM( shell_scales , DIM= 1 )
!    !STOP
!    DO itr= 0, n_shells - 1, 1
!
!      shell_radii( itr + 1 )= 1.0D0 &!last_r*radius &
!          *( ( central_density - itr*density_step )/m_p )**(-third)
!
!    ENDDO
!    !shell_radii= shell_radii*(last_r*radius)/shell_radii( n_shells )
!    PRINT *, n_shells, central_density, surface_density, density_step
!    PRINT *, radius
!    PRINT *, shell_radii
!    !STOP
!    shell_index= 1
!    DO r= 1, 500, 1
!
!      IF( particle_profile( 2, r ) > shell_index )THEN
!        shell_radii( shell_index ) = particle_profile( 1, r )
!        shell_index= shell_index + 1
!      ENDIF
!
!    ENDDO
!    IF( shell_radii( n_shells ) > 0.98D0*radius )THEN
!      shell_radii( n_shells )= shell_radii( n_shells - 1 ) &
!                               + ( shell_radii( n_shells ) &
!                                 - shell_radii( n_shells - 1 ) )/2.0D0
!    ENDIF
!    !PRINT *, radius
!    !PRINT *, shell_radii
!    !STOP

    !---------------------------------------------------------------------!
    !-- Place shells based on mass density a that point, second version --!
    !---------------------------------------------------------------------!

    shell_radii= 0
    shell_radii(1)= ( central_density/m_p )**(-third)
    DO itr= 2, n_shells, 1

      rho_tmp= bns_obj% import_mass_density( center + shell_radii( itr - 1 ), &
                                             0.0D0, 0.0D0 )

      IF( rho_tmp == 0 )THEN
        shell_radii= shell_radii*itr/n_shells
      ENDIF

      shell_radii( itr )= shell_radii( itr - 1 ) + ( rho_tmp/m_p )**(-third)

    ENDDO
    shell_radii= shell_radii*(radius*last_r/shell_radii(n_shells))
    PRINT *, "n_shells= ", n_shells
    PRINT *, "shell_radii= ", shell_radii
    PRINT *
    !STOP
!    shell_scales= shell_scales - shell_scales( 1 )
!    shell_scales= shell_scales/shell_scales( n_shells )
!    PRINT *, shell_scales
!    PRINT *, SUM( shell_scales , DIM= 1 )
!    !STOP
!    DO itr= 0, n_shells - 1, 1
!
!      shell_radii( itr + 1 )= 1.0D0 &!last_r*radius &
!          *( ( central_density - itr*density_step )/m_p )**(-third)
!
!    ENDDO
!    !shell_radii= shell_radii*(last_r*radius)/shell_radii( n_shells )
!    PRINT *, n_shells, central_density, surface_density, density_step
!    PRINT *, radius
!    PRINT *, shell_radii
!    !STOP
!    shell_index= 1
!    DO r= 1, 500, 1
!
!      IF( particle_profile( 2, r ) > shell_index )THEN
!        shell_radii( shell_index ) = particle_profile( 1, r )
!        shell_index= shell_index + 1
!      ENDIF
!
!    ENDDO
!    IF( shell_radii( n_shells ) > 0.98D0*radius )THEN
!      shell_radii( n_shells )= shell_radii( n_shells - 1 ) &
!                               + ( shell_radii( n_shells ) &
!                                 - shell_radii( n_shells - 1 ) )/2.0D0
!    ENDIF
    !PRINT *, radius
    !PRINT *, shell_radii
    !STOP

    PRINT *, " * Integrating the baryon mass density to get the mass profile..."
    PRINT *

    !CALL OMP_SET_NUM_THREADS(80)

    dr             = radius/500.0D0
    dth            = pi/2.0D0/250.0D0
    dphi           = 2.0D0*pi/500.0D0
    CALL bns_obj% integrate_baryon_mass_density( center, radius, &
                                                 central_density, &
                                                 dr, dth, dphi, &
                                                 mass, mass_profile, &
                                                 mass_profile_idx )

    mass_profile( 2:3, : )= mass_profile( 2:3, : )*mass_star/mass

    !PRINT *, "mass_profile( 2, mass_profile_idx(NINT(radius/dr)) )= ", &
    !         mass_profile( 2, mass_profile_idx(NINT(radius/dr)) )
    !DO itr= 1, NINT(radius/dr), 1
    !  PRINT *, mass_profile( 3, mass_profile_idx(itr) )
    !ENDDO
    !STOP

    !-------------------------------------------------!
    !-- Place shells based on constant step on mass --!
    !-------------------------------------------------!

 !   shell_radii= 0.0D0
 !   DO itr= 1, n_shells, 1
 !
 !     shell_radii( itr )= (( radius*last_r )/DBLE(n_shells))*DBLE( itr - 1/2 )
 !
 !   ENDDO
 !   shell_thickness= shell_radii( 2 ) - shell_radii( 1 )

    !surface_density= bns_obj% import_mass_density( center + radius, &
    !                                               0.0D0, 0.0D0 )
    !density_step= ( central_density - surface_density )/DBLE( n_shells - 1 )

    shell_index= 1
    itr2= 0
    shell_masses= 0.0D0
    DO itr= 0, NINT(radius/dr), 1

      IF( shell_index == n_shells )THEN

        shell_masses( shell_index )= SUM( mass_profile( 2, &
         mass_profile_idx(itr2):mass_profile_idx(NINT(radius/dr)-1) ), DIM= 1 )

        EXIT

      ENDIF

      IF( mass_profile( 1, mass_profile_idx(itr) ) &
          >= shell_radii( shell_index ) &!+ radius/DBLE(2*n_shells)
      )THEN

       shell_masses( shell_index )= SUM( mass_profile( 2, &
                     mass_profile_idx(itr2):mass_profile_idx(itr) ), DIM= 1 )

       itr2= itr + 1
       shell_index= shell_index + 1

      ENDIF

    ENDDO
    IF( ABS( SUM( shell_masses, DIM= 1 ) - mass_star )/mass_star > 5.0D-3 )THEN
      PRINT *, " ** The masses of the shells do not add up to the ", &
               "mass of the star. Stopping..."
      PRINT *, "SUM( shell_masses )= ", SUM( shell_masses, DIM=1 )
      PRINT *, "mass_star= ", mass_star
      PRINT *, "shell_masses=", shell_masses
      PRINT *
      STOP
    ENDIF
    !STOP

    ALLOCATE( npart_shell( n_shells ) )
    ALLOCATE( npart_shelleq( n_shells ) )
    ALLOCATE( alpha( n_shells ) )
    ALLOCATE( colatitude_pos( n_shells ) )
    ALLOCATE( pos_shells( n_shells ) )

  CALL bns_obj% import_id( &
           center, 0.0D0, 0.0D0, &
           gxx_tmp, &
           baryon_density_tmp, &
           gamma_euler_tmp )

    npart_shelleq(1)= CEILING( SQRT(DBLE(2*shell_masses(1)/m_p)) )
    PRINT *, npart_shelleq(1), ( npart_shelleq(1)**2 )/2
    IF( MOD( npart_shelleq(1), 2 ) /= 0 )THEN
      npart_shelleq(1)= npart_shelleq(1) + 1
    ENDIF
    IF( MOD( npart_shelleq(1)/2, 2 ) /= 0 )THEN
      npart_shelleq(1)= 2*( npart_shelleq(1)/2 + 1 )
    ENDIF
    PRINT *, ( npart_shelleq(1)**2 )/2
    PRINT *, ( 4.0D0/3.0D0*pi*shell_radii(1)**3.0D0 )* &
             gxx_tmp*SQRT(gxx_tmp)* &
             baryon_density_tmp*gamma_euler_tmp/(shell_masses(1)/(( npart_shelleq(1)**2 )/2))
    PRINT *
    PRINT *, shell_masses(1), ( 4.0D0/3.0D0*pi*shell_radii(1)**3.0D0 )* &
             gxx_tmp*SQRT(gxx_tmp)* &
             baryon_density_tmp*gamma_euler_tmp, (( 4.0D0/3.0D0*pi*shell_radii(1)**3.0D0 )* &
             gxx_tmp*SQRT(gxx_tmp)* &
             baryon_density_tmp*gamma_euler_tmp)/shell_masses(1)
    PRINT *
    PRINT *, gxx_tmp*SQRT(gxx_tmp)*baryon_density_tmp*gamma_euler_tmp, &
             shell_masses(1)/( 4.0D0/3.0D0*pi*shell_radii(1)**3.0D0 ), &
             gxx_tmp*SQRT(gxx_tmp)*baryon_density_tmp*gamma_euler_tmp/ &
             (shell_masses(1)/( 4.0D0/3.0D0*pi*shell_radii(1)**3.0D0 ))
    !STOP

    npart_shelleq= 0
    DO r= 1, n_shells, 1

      npart_shelleq( r )= CEILING( SQRT(DBLE(2*shell_masses( r )/m_p)) )
      !IF( itr == n_shells - 1 ) npart_shelleq( itr )= &
      !                      npart_shelleq( n_shells - 2 )
      !                      !MAXVAL( npart_shelleq )
      !                      !NINT( 1.5D0*npart_shelleq( itr ) )
      !IF( itr == n_shells ) npart_shelleq( itr )= &
      !                      !npart_shelleq( n_shells - 2 )
      !                      MAXVAL( npart_shelleq )
      !                      !NINT( 4.0D0*npart_shelleq( itr ) )
      IF( MOD( npart_shelleq( r ), 2 ) /= 0 )THEN
        npart_shelleq( r )= npart_shelleq( r ) + 1
      ENDIF
      IF( MOD( npart_shelleq( r )/2, 2 ) /= 0 )THEN
        npart_shelleq( r )= 2*( npart_shelleq( r )/2 + 1 )
      ENDIF
      IF( MOD( npart_shelleq( r ), 4 ) /= 0 )THEN
        PRINT *, "npart_shelleq(", r, ")=", npart_shelleq( r ), &
                 " is not divisible by 4 in the initialization loop. ", &
                 " Check the algorithm. Stopping..."
        STOP
      ENDIF
      IF( npart_shelleq( r ) == 0 )THEN
        alpha( r )= 0.0D0
      ELSE
        alpha( r )= 2.0D0*pi/DBLE(npart_shelleq( r ))
      ENDIF
      npart_shell( r )= ( npart_shelleq( r )**2.0D0 )/2.0D0

      ALLOCATE( colatitude_pos( r )% colatitudes( npart_shelleq( r )/4 ) )

      DO itr2= 1, npart_shelleq( r )/4, 1

        colatitude_pos( r )% colatitudes( itr2 )= &
                      ACOS( 2.0D0*itr2/(npart_shelleq( r )/2 + 1.0D0 )&
                          - 1.0D0 )
      ENDDO

      ALLOCATE( pos_shells( r )% pos_shell( 3, npart_shell( r ) ) )
      ALLOCATE( pos_shells( r )% pvol_shell( npart_shell( r ) ) )
      ALLOCATE( pos_shells( r )% g_xx( npart_shell( r ) ) )
      ALLOCATE( pos_shells( r )% baryon_density( npart_shell( r ) ) )
      ALLOCATE( pos_shells( r )% gamma_euler( npart_shell( r ) ) )

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
    m_parts( 0 )= m_parts( NINT(DBLE(n_shells/2)) )!first_shell= NINT(DBLE(n_shells/2))
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

  IF( PRESENT(filename_mass_profile) )THEN
    finalnamefile= filename_mass_profile
  ELSE
    finalnamefile= "mass_profile.dat"
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

  IF( PRESENT(filename_shells_radii) )THEN
    finalnamefile= filename_shells_radii
  ELSE
    finalnamefile= "shell_radii.dat"
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

  !STOP

  PRINT *, " * Assigning first half of particle positions..."

    DO r= 1, n_shells, 1
      IF( ALLOCATED( pos_shells( r )% pos_shell ) ) &
        DEALLOCATE( pos_shells( r )% pos_shell )
        !PRINT *, "0.1"
      IF( ALLOCATED( pos_shells( r )% pvol_shell ) ) &
        DEALLOCATE( pos_shells( r )% pvol_shell )
                !PRINT *, "0.2"
      IF( ALLOCATED( pos_shells( r )% pvol_shell2 ) ) &
        DEALLOCATE( pos_shells( r )% pvol_shell2 )
                !PRINT *, "0.3"
      IF( ALLOCATED( pos_shells( r )% g_xx ) )&
        DEALLOCATE( pos_shells( r )% g_xx )
                !PRINT *, "0.4"
      IF( ALLOCATED( pos_shells( r )% baryon_density ) ) &
        DEALLOCATE( pos_shells( r )% baryon_density )
                !PRINT *, "0.5"
      IF( ALLOCATED( pos_shells( r )% gamma_euler ) ) &
        DEALLOCATE( pos_shells( r )% gamma_euler )
                      !PRINT *, "0.6"
      IF( ALLOCATED( pos_shells( r )% pos_th ) ) &
        DEALLOCATE( pos_shells( r )% pos_th )
                      !PRINT *, "0.7"
      IF( ALLOCATED( pos_shells( r )% pos_phi ) ) &
        DEALLOCATE( pos_shells( r )% pos_phi )

      ALLOCATE( pos_shells( r )% pos_shell     ( 3, npart_approx ) )
      ALLOCATE( pos_shells( r )% pvol_shell    ( npart_approx ) )
      ALLOCATE( pos_shells( r )% pvol_shell2   ( npart_approx ) )
      ALLOCATE( pos_shells( r )% g_xx          ( npart_approx ) )
      ALLOCATE( pos_shells( r )% baryon_density( npart_approx ) )
      ALLOCATE( pos_shells( r )% gamma_euler   ( npart_approx ) )
      ALLOCATE( pos_shells( r )% pos_th        ( npart_approx ) )
      ALLOCATE( pos_shells( r )% pos_phi       ( npart_approx ) )

      pos_shells(r)% pos_shell= 0.0D0
      pos_shells(r)% pos_phi= -1.0D0
      pos_shells(r)% pos_th= -1.0D0
      pos_shells(r)% pvol_shell= 0.0D0
      pos_shells(r)% pvol_shell2= 0.0D0
      pos_shells(r)% g_xx= 0.0D0
      pos_shells(r)% baryon_density= 0.0D0
      pos_shells(r)% gamma_euler= 0.0D0
      m_parts( r )= m_p
      npart_shelleq( r )= CEILING( SQRT(DBLE(2*shell_masses( r )/m_parts( r ))))

    ENDDO
    npart_shelleq( 1 )= n_particles_first_shell !4

    pos  = 0.0D0
    pmass = 0.0D0
    phase= 0.0D0
    proper_volume= 0.0D0
    vol_shell= 0.0D0
    vol_shell2= 0.0D0
    dr_shells= radius/n_shells
    npart_out= 0
    !upper_bound= 1.025D0
    !lower_bound= 0.975D0
    !upper_factor= 1.01D0
    !lower_factor= 0.99D0
    upper_bound_tmp= upper_bound
    lower_bound_tmp= lower_bound
    r= CEILING(DBLE(n_shells)/2.0D0)!first_shell
    cnt2= 0
    r_cnt= 1

    ALLOCATE( pos_shell_tmp  ( 3, CEILING(SQRT(DBLE(2*npart_approx))), CEILING(SQRT(DBLE(2*npart_approx))) ) )
    ALLOCATE( g_xx_tmp       ( CEILING(SQRT(DBLE(2*npart_approx))), CEILING(SQRT(DBLE(2*npart_approx))) ) )
    ALLOCATE( bar_density_tmp( CEILING(SQRT(DBLE(2*npart_approx))), CEILING(SQRT(DBLE(2*npart_approx))) ) )
    ALLOCATE( gam_euler_tmp  ( CEILING(SQRT(DBLE(2*npart_approx))), CEILING(SQRT(DBLE(2*npart_approx))) ) )
    ALLOCATE( pvol_tmp       ( CEILING(SQRT(DBLE(2*npart_approx))), CEILING(SQRT(DBLE(2*npart_approx))) ) )

    !----------------------!
    !--  Main iteration  --!
    !----------------------!

    !CALL OMP_SET_NUM_THREADS(80)

    IF( debug ) PRINT *, THIS% randomize_phi
    IF( debug ) PRINT *, THIS% randomize_theta
    IF( debug ) PRINT *, THIS% randomize_r

    find_desired_npart: DO

      !DO r= 1, n_shells, 1
      !  m_parts( r )= m_p
      !  npart_shelleq( r )= &
      !            CEILING( SQRT(DBLE(2*shell_masses( r )/m_parts( r ))))
      !ENDDO
      !npart_shelleq( 1 )= n_particles_first_shell !4
      !
      !r= 1

    place_particles_on_shells: DO

      pos_shell_tmp  = HUGE(0.0D0)
      g_xx_tmp       = 0.0D0
      bar_density_tmp= 0.0D0
      gam_euler_tmp  = 0.0D0
      pvol_tmp       = 0.0D0

      !IF( r_cnt == 2 ) r= 1

!PRINT *, "Start of iteration, shell ", r, "iteration ", cnt2 + 1
!PRINT *
      !npart_shelleq( r )= NINT( correction*npart_shelleq( r ) )
      !IF( r == 1 ) npart_shelleq( r )= npart_shelleq( r )/1.5
      IF( MOD( npart_shelleq( r ), 2 ) /= 0 )THEN
        CALL RANDOM_NUMBER( rand_num2 )
        IF( rand_num2 >= half ) rel_sign=  1
        IF( rand_num2 < half )  rel_sign= -1
        npart_shelleq( r )= npart_shelleq( r ) + rel_sign
      ENDIF
      IF( MOD( npart_shelleq( r )/2, 2 ) /= 0 )THEN
        CALL RANDOM_NUMBER( rand_num2 )
        IF( rand_num2 >= half ) rel_sign=  1
        IF( rand_num2 < half )  rel_sign= -1
        npart_shelleq( r )= 2*( npart_shelleq( r )/2 + rel_sign )
      ENDIF
      IF( MOD( npart_shelleq( r ), 4 ) /= 0 )THEN
        PRINT *, "npart_shelleq(", r, ")=", npart_shelleq( r ), &
                 " is not divisible by 4 in the main iteration. ", &
                 " Check the algorithm. Stopping..."
        STOP
      ENDIF
      IF( npart_shelleq( r ) == 0 )THEN
        alpha( r )= 0.0D0
      ELSE
        alpha( r )= 2.0D0*pi/DBLE(npart_shelleq( r ))
      ENDIF
      npart_shell( r )= ( npart_shelleq( r )**2.0D0 )/2.0D0

      IF( ALLOCATED( colatitude_pos( r )% colatitudes ) ) &
        DEALLOCATE( colatitude_pos( r )% colatitudes )
      ALLOCATE( colatitude_pos( r )% colatitudes( npart_shelleq( r )/4 ) )

      DO itr2= 1, npart_shelleq( r )/4, 1

        colatitude_pos( r )% colatitudes( itr2 )= &
                      ACOS( 2.0D0*DBLE(itr2)/ &
                            (DBLE(npart_shelleq( r )/2) + 0.625D0 ) - 1.0D0 )
          !            alpha( r )*1.0D0/2.0D0 + ( itr2 - 1 )*alpha( r )
        !  ACOS( 2.0D0*( 1.0D0 - COS( pi/3.0D0*( 2.0D0/3.0D0 + DBLE(itr2 - 1)*DBLE(npart_shelleq( r )/4 + 1.0D0 -(1.0D0/2.0D0)-(2.0D0/3.0D0) )/DBLE(npart_shelleq( r )/4 - 1.0D0 ) ) &
        !                   /DBLE(npart_shelleq( r )/4 + 1.0D0 ) ) ) &
        !      - 1.0D0 )
              !5.0D0/12.0D0
        !colatitude_pos( r )% colatitudes( itr2 )= &
        !              colatitude_pos( r )% colatitudes( itr2 ) &
        !              *( 1 + rel_sign*0.05D0*phase_th )
        IF( colatitude_pos( r )% colatitudes( itr2 ) <= pi/2 .OR. &
            colatitude_pos( r )% colatitudes( itr2 ) >= pi &
        )THEN
          PRINT *, "The colatitudes are not in the interval (pi/2,pi). ", &
                   "Stopping..."
          STOP
        ENDIF
      ENDDO
      !PRINT *, npart_shell( r ), npart_shelleq( r )
      !PRINT *

      !phase= phase + r*alpha(r)/golden_ratio
      !CALL RANDOM_NUMBER( phase )
      !phase= phase*alpha(r)
      !dphi_shells= alpha(r)
      !rad= shell_radii(r)
      itr= 0
      npart_discard= 0

      npart_shell_tmp= npart_shell( r )
      !PRINT *, npart_shell_tmp

      !CALL RANDOM_NUMBER( phase_th )
      !phase_th= phase_th*ABS( MINVAL( colatitude_pos(r)% colatitudes, DIM= 1 ) &
      !                     - pi/2.0D0 )*4.0D0/5.0D0

      !PRINT *, "Right before OMP, shell ", r, "iteration ", cnt2 + 1
      !PRINT *

      !$OMP PARALLEL DO DEFAULT(NONE), &
      !$OMP             PRIVATE(phase,col,col_tmp,xtemp,ytemp,ztemp, &
      !$OMP                     dphi_shells,dth_shells, delta_r, &
      !$OMP                     th,phi,rand_num2,phase_th,rel_sign), &
      !$OMP             SHARED(r,npart_shelleq,center,rad,alpha,pos_shells, &
      !$OMP                    colatitude_pos,bns_obj,n_shells, &
      !$OMP                    dr_shells, shell_radii, shell_thickness, THIS, &
      !$OMP                    g_xx_tmp, bar_density_tmp, gam_euler_tmp, &
      !$OMP                    pos_shell_tmp, pvol_tmp ), &
      !$OMP             REDUCTION(+:npart_discard)
      DO th= 1, npart_shelleq( r )/4, 1 !npart_shelleq( r ) is even, see above

        !PRINT *, "itr= ", itr
        !STOP

        IF( THIS% randomize_phi )THEN

          CALL RANDOM_NUMBER( phase )
          phase= phase*alpha(r)

        ENDIF

        dphi_shells= alpha(r)

        !col= colatitude_pos(r)% colatitudes(th)
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
        !PRINT *, "Right before longitude loop"

        DO phi= 1, npart_shelleq( r ), 1

          col= colatitude_pos(r)% colatitudes(th)
          IF( THIS% randomize_theta )THEN

            CALL RANDOM_NUMBER( phase_th )
            CALL RANDOM_NUMBER( rand_num2 )
            IF( rand_num2 >= half ) rel_sign=  1
            IF( rand_num2 < half )  rel_sign= -1

            col_tmp= col*( 1.0D0 + rel_sign*0.05D0*phase_th )

            IF( col_tmp < pi .AND. col_tmp > pi/2.0D0 )THEN

              col= col_tmp

            ENDIF

          ENDIF

          rad= shell_radii(r)
          IF( THIS% randomize_r )THEN

            CALL RANDOM_NUMBER( delta_r )
            CALL RANDOM_NUMBER( rand_num2 )
            IF( rand_num2 >= half ) rel_sign=  1
            IF( rand_num2 < half )  rel_sign= -1

            !PRINT *, rad
            IF( r/n_shells < 0.75D0 )THEN
              rad= rad + rel_sign*delta_r*0.35D0*dr_shells
            ENDIF
            !PRINT *, rad, shell_thickness, rel_sign
            !PRINT *

          ENDIF

          IF( rad < 0 )THEN
            PRINT *, " * ERROR! rad < 0. Check the computation of the radial", &
                     " coordinates of the particles. Stopping.."
            STOP
          ENDIF

        !PRINT *, "2.1"

          xtemp= center + rad*COS(phase + phi*alpha(r))*SIN(col)
          ytemp= rad*SIN(phase + phi*alpha(r))*SIN(col)
          ztemp= rad*COS(col)

          IF( ISNAN( xtemp ) )THEN
            PRINT *, "** ERROR when placing first half of the particles! ", &
                     "xtemp is a NaN. Stopping.."
            STOP
          ENDIF
          IF( ISNAN( ytemp ) )THEN
            PRINT *, "** ERROR when placing first half of the particles! ", &
                     "ytemp is a NaN. Stopping.."
            STOP
          ENDIF
          IF( ISNAN( ztemp ) )THEN
            PRINT *, "** ERROR when placing first half of the particles! ", &
                     "ztemp is a NaN. Stopping.."
            STOP
          ENDIF

!PRINT *, "2.2"

          !CALL bns_obj% import_id( &
          !         xtemp, ytemp, ztemp, &
          !         pos_shells(r)% g_xx( itr + 1 ), &
          !         pos_shells(r)% baryon_density( itr + 1 ), &
          !         pos_shells(r)% gamma_euler( itr + 1 ) )

          CALL bns_obj% import_id( &
                   xtemp, ytemp, ztemp, &
                   g_xx_tmp( th, phi ), &
                   bar_density_tmp( th, phi ), &
                   gam_euler_tmp( th, phi ) )

          IF( bar_density_tmp( th, phi ) > 0.0D0 &
              !pos_shells(r)% baryon_density( itr + 1 ) > 0.0D0 &
              .AND. &
              bns_obj% is_hydro_negative( xtemp, ytemp, ztemp ) == 0 )THEN

            IF( bar_density_tmp( th, phi ) == 0 )THEN
                !pos_shells(r)% baryon_density( itr + 1 ) == 0 )THEN
              PRINT *, "When placing first half of particles"
              PRINT *, r, th, phi, xtemp, ytemp, ztemp, &
                       bar_density_tmp( th, phi )
            ENDIF

!PRINT *, "2.3"

            !npart_out= npart_out + 1
            !pos_shells(r)% pos_shell( 1, itr + 1 )= xtemp
            !pos_shells(r)% pos_shell( 2, itr + 1 )= ytemp
            !pos_shells(r)% pos_shell( 3, itr + 1 )= ztemp

            pos_shell_tmp( 1, th, phi )= xtemp
            pos_shell_tmp( 2, th, phi )= ytemp
            pos_shell_tmp( 3, th, phi )= ztemp


            IF( th == 1 )THEN

              !dth_shells= pi - ( col + colatitude_pos(r)% colatitudes(th+1) )/2.0D0
              dth_shells= 2.0D0*ABS( col - &
                        ( col + colatitude_pos(r)% colatitudes(th + 1) )/2.0D0 )

              !PRINT *, "case 1. ", "dt_shells= ", dth_shells, &
              !                     "col= ", col, &
              !                     "colatitude_pos(r)% colatitudes(th + 1)=", &
              !                     colatitude_pos(r)% colatitudes(th + 1)
              !PRINT *

            ELSEIF( th == npart_shelleq( r )/4 )THEN

              !dth_shells= ( colatitude_pos(r)% colatitudes(th-1) + col - pi )/2.0D0
              dth_shells= 2.0D0*ABS( ( colatitude_pos(r)% colatitudes(th - 1) &
                        + col )/2.0D0 - col )

              !PRINT *, "case 2. ", "dt_shells= ", dth_shells, &
              !                     "col= ", col, &
              !                     "colatitude_pos(r)% colatitudes(th - 1)=", &
              !                     colatitude_pos(r)% colatitudes(th - 1)
              !PRINT *

            ELSE

              dth_shells= ABS( &
                      ( colatitude_pos(r)% colatitudes(th + 1) + col )/2.0D0 &
                    - ( col + colatitude_pos(r)% colatitudes(th - 1) )/2.0D0 )


              !PRINT *, "case 3. ", "dt_shells= ", dth_shells, &
              !                     "col= ", col, &
              !                     "colatitude_pos(r)% colatitudes(th + 1)=", &
              !                     colatitude_pos(r)% colatitudes(th + 1), &
              !                     "colatitude_pos(r)% colatitudes(th - 1)=", &
              !                     colatitude_pos(r)% colatitudes(th - 1)
              !PRINT *

            ENDIF

            !pvol( npart_out )= rad**2.0D0*SIN(col) &
            !                   *dr_shells*dth_shells*dphi_shells

         ! This is commented out to parallelize the loop
         !   pos_shells(r)% pvol_shell2( itr + 1 )= rad**2.0D0*SIN(col) &
         !                      *dr_shells*dth_shells*dphi_shells! &
         !                      !*pos_shells(r)% g_xx( itr ) &
         !                      !*SQRT(pos_shells(r)% g_xx( itr ))

             pvol_tmp( th, phi )= rad**2.0D0*SIN(col) &
                                     *dr_shells*dth_shells*dphi_shells! &

            IF( pvol_tmp( th, phi ) <= 0 )THEN
                ! pos_shells(r)% pvol_shell2( itr + 1 ) <= 0 )THEN
              PRINT *, "When placing first half of particles"
              PRINT *, "pvol_tmp( ", r, th, phi, " ) =", &
                       pvol_tmp( th, phi )
              STOP
            ENDIF

!PRINT *, "2.4"
            !itr= itr + 1
            !PRINT *, itr, npart_out, col/pi, phi*alpha(r)/pi

          ELSE

            npart_discard= npart_discard + 2
            !PRINT *, "removed"

          ENDIF

          ! Print progress on screen, every 10%
          !perc= 50*( THIS% nx*THIS% ny*iz + THIS% nx*iy + ix )/ &
          !        ( THIS% nx*THIS% ny*THIS% nz/2 )
          !IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
          !  WRITE( *, "(A2,I3,A1)", ADVANCE= "NO" ) &
          !         creturn//" ", perc, "%"
          !ENDIF

        ENDDO
        !PRINT *, "Right after longitude loop"
      ENDDO
      !$OMP END PARALLEL DO

      npart_shell( r )= npart_shell( r ) - npart_discard

      !PRINT *, "Right after OMP"

    !  IF( itr /= npart_shell( r )/2 )THEN
    !    PRINT *, "** ERROR! Mismatch in the particle counters on shell ", r
    !    PRINT *, "itr=", itr, "npart_shell( r )/2=", npart_shell( r )/2
    !    PRINT *, "itr should be equal to npart_shell( r )/2. Stopping..."
    !    PRINT *
    !    STOP
    !  ENDIF
      !PRINT *, npart_out

      !WRITE( *, "(A2,I4,I4,F3.2)", ADVANCE= "NO" ) &
      !          creturn//" ", npart_shell_tmp, npart_shell( r ), &
      !                        DBLE(npart_shell( r ))/DBLE(npart_shell_tmp)

      !PRINT *, "Right before safety check"

      IF( npart_shell( r ) < 0 ) npart_shell( r )= 0
      IF( npart_shell( r ) == 0 )THEN
        m_parts( r )= m_parts( prev_shell )
        PRINT *, " * Placed", npart_shell( r )/2, &
                 " particles on one emisphere of spherical shell ", r, &
                 " out of ", n_shells
        IF( r == 1 )THEN
          EXIT
        ELSEIF( r < CEILING(DBLE(n_shells)/2.0D0) )THEN
          !PRINT *, "r=", r
          r= r - 1
          cnt2 = 0
          upper_bound_tmp= upper_bound
          lower_bound_tmp= lower_bound
          CYCLE
        ELSEIF( r == n_shells )THEN
          !PRINT *, "r=", r
          r= CEILING(DBLE(n_shells)/2.0D0) - 1
          r_cnt= r_cnt + 1
          cnt2 = 0
          upper_bound_tmp= upper_bound
          lower_bound_tmp= lower_bound
          CYCLE
        ELSEIF( r >= CEILING(DBLE(n_shells)/2.0D0) )THEN
          !PRINT *, "r=", r
          r= r + 1
          cnt2 = 0
          upper_bound_tmp= upper_bound
          lower_bound_tmp= lower_bound
          CYCLE
        ENDIF
      ELSE
        m_parts( r )= shell_masses( r )/DBLE(npart_shell( r ))
      ENDIF

      !IF( r == 2 .AND. m_parts( r ) > m_parts( r - 1 ) )THEN
      !  npart_shelleq( r )= &
      !              NINT( m_parts( r )/m_parts( r - 1 )*npart_shelleq( r ) )
      !  npart_out= npart_out - npart_shell( r )/2
      !  CYCLE
      !ENDIF

      IF( debug ) PRINT *, "10.9"

      itr= 0
      DO th= 1, npart_shelleq( r )/4, 1
        DO phi= 1, npart_shelleq( r ), 1

          IF( pos_shell_tmp( 1, th, phi ) < HUGE(0.0D0) )THEN

            itr= itr + 1
            pos_shells(r)% pos_shell( 1, itr )= &
                                              pos_shell_tmp( 1, th, phi )
            pos_shells(r)% pos_shell( 2, itr )= &
                                              pos_shell_tmp( 2, th, phi )
            pos_shells(r)% pos_shell( 3, itr )= &
                                              pos_shell_tmp( 3, th, phi )

            pos_shells(r)% g_xx( itr )          = &
                                              g_xx_tmp( th, phi )
            pos_shells(r)% baryon_density( itr )= &
                                              bar_density_tmp( th, phi )
            pos_shells(r)% gamma_euler( itr )   = &
                                              gam_euler_tmp( th, phi )
            pos_shells(r)% pvol_shell2( itr )   = &
                                              pvol_tmp( th, phi )


          ENDIF

        ENDDO
      ENDDO
      npart_out= npart_out + itr

      IF( debug ) PRINT *, "11"

      !PRINT *, "Right before correction of particle number"
      !PRINT *, npart_out

      !IF( r > 1 )THEN
      IF( r /= CEILING(DBLE(n_shells)/2.0D0) )THEN

        !IF( r_cnt == 2 )THEN
        !  high_mass= m_parts( r )/m_parts( first_shell ) > upper_bound_tmp
        !  low_mass = m_parts( r )/m_parts( first_shell ) < lower_bound_tmp
        !ELSE

        IF( r < CEILING(DBLE(n_shells)/2.0D0) )THEN
          prev_shell= r + 1
        ELSEIF( r > CEILING(DBLE(n_shells)/2.0D0) )THEN
          prev_shell= r - 1
        ELSEIF( r == 1 )THEN
          EXIT
        ENDIF

          high_mass= m_parts( r )/m_parts( prev_shell ) > upper_bound_tmp
          low_mass = m_parts( r )/m_parts( prev_shell ) < lower_bound_tmp
        !ENDIF
        npart_shell_kept= DBLE(npart_shell( r ))/DBLE(npart_shell_tmp)
        kept_all = npart_shell_kept == 1.0D0

        !PRINT *, "cnt2=", cnt2
        !PRINT *, "upper_bound_tmp=", upper_bound_tmp
        !PRINT *, "lower_bound_tmp=", lower_bound_tmp
        !PRINT *, "n_shells=", n_shells
        !PRINT *, "r=", r
        !PRINT *, "npart_shell( r )=", npart_shell( r )
        !PRINT *, "npart_shell_tmp=", npart_shell_tmp
        !PRINT *, "npart_shell_kept=", npart_shell_kept
        !PRINT *, "high_mass=", high_mass
        !PRINT *, "low_mass=", low_mass
        !PRINT *, "kept_all=", kept_all
        !PRINT *, "m_parts( r )=", m_parts( r )
        !PRINT *, "m_parts( r - 1 )=", m_parts( r - 1 )
        !PRINT *, " m_parts( r )/m_parts( r - 1 )= ",  &
        !                           m_parts( r )/m_parts( r - 1 )
        !PRINT *

        IF( high_mass .AND. kept_all )THEN
!PRINT *, "case 1"

          cnt2= cnt2 + 1
!PRINT *, "1.2"
          IF( cnt2 > max_steps )THEN
            upper_bound_tmp= upper_bound_tmp*upper_factor
            lower_bound_tmp= lower_bound_tmp*lower_factor
            IF( r > 0.8D0*n_shells .AND. &
                m_parts( r )/m_parts( prev_shell ) > 1.1D0*upper_bound &
                !upper_bound_tmp > 1.1D0*upper_bound &
            )THEN
              CALL RANDOM_NUMBER( rand_num2 )
              rand_num= NINT( 5.0D0*( rand_num2 + 1.0D0 ) )
              npart_shelleq( r )= rand_num*npart_shelleq( r )
                                !  npart_shelleq( r - 1 ) &
                                !+ rel_sign*NINT( 1 + rand_num )
              upper_bound_tmp= upper_bound
              lower_bound_tmp= lower_bound
              npart_out= npart_out - npart_shell( r )/2
              cnt2= 1
              CYCLE
            ENDIF
            cnt2= 1
          ENDIF
!PRINT *, "1.4"
          npart_out= npart_out - npart_shell( r )/2
!PRINT *, "1.6"
          CALL RANDOM_NUMBER( rand_num )
          CALL RANDOM_NUMBER( rand_num2 )
          npart_shelleq( r )= npart_shelleq( r ) + 1*NINT( 1 + 1.0*rand_num ) &
                                                 + 1*NINT( 1 + 1.0*rand_num2 )
!PRINT *, "1.8"
          IF( npart_shelleq( r ) == 0 .OR. npart_shell( r ) == 0 )THEN
            CALL RANDOM_NUMBER( rand_num )
            CALL RANDOM_NUMBER( rand_num2 )
            IF( rand_num2 < half )  rel_sign= - 1
            IF( rand_num2 >= half ) rel_sign=   1
            npart_shelleq( r )= npart_shelleq( r - 1 ) &
                              + rel_sign*NINT( 1 + rand_num )
          ENDIF
!PRINT *, "1.9"
          CYCLE

        ELSEIF( low_mass .AND. kept_all )THEN
        !PRINT *, "case 2"

          cnt2= cnt2 + 1
          IF( cnt2 > max_steps )THEN
            upper_bound_tmp= upper_bound_tmp*upper_factor
            lower_bound_tmp= lower_bound_tmp*lower_factor
            IF( r > 0.8D0*n_shells .AND. &
                m_parts( r )/m_parts( prev_shell ) < 0.9D0*lower_bound &
                !lower_bound_tmp < 0.9D0*lower_bound &
            )THEN
              CALL RANDOM_NUMBER( rand_num2 )
              rand_num= NINT( 5.0D0*( rand_num2 + 1.0D0 ) )
              npart_shelleq( r )= npart_shelleq( r )/rand_num
                                !  npart_shelleq( r - 1 ) &
                                !+ rel_sign*NINT( 1 + rand_num )
              upper_bound_tmp= upper_bound
              lower_bound_tmp= lower_bound
              npart_out= npart_out - npart_shell( r )/2
              cnt2= 1
              CYCLE
            ENDIF
            cnt2= 1
          ENDIF

          npart_out= npart_out - npart_shell( r )/2

          CALL RANDOM_NUMBER( rand_num )
          CALL RANDOM_NUMBER( rand_num2 )
          npart_shelleq( r )= npart_shelleq( r ) - 1*NINT( 1 + 1.0*rand_num ) &
                                                 - 1*NINT( 1 + 1.0*rand_num2 )

          IF( npart_shelleq( r ) == 0 .OR. npart_shell( r ) == 0 )THEN
            CALL RANDOM_NUMBER( rand_num )
            CALL RANDOM_NUMBER( rand_num2 )
            IF( rand_num2 < half )  rel_sign= - 1
            IF( rand_num2 >= half ) rel_sign=   1
            npart_shelleq( r )= npart_shelleq( r - 1 ) &
                              + rel_sign*NINT( 1 + rand_num )
          ENDIF

          CYCLE

        ELSEIF( high_mass .AND. .NOT.kept_all ) THEN
        !PRINT *, "case 3"

          cnt2= cnt2 + 1
          IF( cnt2 > max_steps )THEN
            upper_bound_tmp= upper_bound_tmp*upper_factor
            lower_bound_tmp= lower_bound_tmp*lower_factor
            IF( r > 0.8D0*n_shells .AND. &
                m_parts( r )/m_parts( prev_shell ) > 1.1D0*upper_bound &
                !upper_bound_tmp > 1.1D0*upper_bound &
            )THEN
              CALL RANDOM_NUMBER( rand_num2 )
              rand_num= NINT( 5.0D0*( rand_num2 + 1.0D0 ) )
              npart_shelleq( r )= rand_num*npart_shelleq( r )
                                !  npart_shelleq( r - 1 ) &
                                !+ rel_sign*NINT( 1 + rand_num )
              upper_bound_tmp= upper_bound
              lower_bound_tmp= lower_bound
              npart_out= npart_out - npart_shell( r )/2
              cnt2= 1
              CYCLE
            ENDIF
            cnt2= 1
          ENDIF

          npart_out= npart_out - npart_shell( r )/2

          CALL RANDOM_NUMBER( rand_num )
          CALL RANDOM_NUMBER( rand_num2 )
          IF( rand_num2 < half )  rel_sign= - 1
          IF( rand_num2 >= half ) rel_sign=   1
          npart_shelleq( r )= CEILING( SQRT( &
                                2*(shell_masses( r )/m_parts( r - 1 )) &
                                /npart_shell_kept &
                              ) ) + rel_sign*NINT( 1 + rand_num )

          IF( npart_shelleq( r ) == 0 .OR. npart_shell( r ) == 0 )THEN
            CALL RANDOM_NUMBER( rand_num )
            CALL RANDOM_NUMBER( rand_num2 )
            IF( rand_num2 < half )  rel_sign= - 1
            IF( rand_num2 >= half ) rel_sign=   1
            npart_shelleq( r )= npart_shelleq( r - 1 ) &
                              + rel_sign*NINT( 1 + rand_num )
          ENDIF

          CYCLE

        ELSEIF( low_mass .AND. .NOT.kept_all ) THEN
        !PRINT *, "case 4"

          cnt2= cnt2 + 1
          IF( cnt2 > max_steps )THEN
            upper_bound_tmp= upper_bound_tmp*upper_factor
            lower_bound_tmp= lower_bound_tmp*lower_factor
            IF( r > 0.8D0*n_shells .AND. &
                m_parts( r )/m_parts( prev_shell ) < 0.9D0*lower_bound &
                !lower_bound_tmp < 0.9D0*lower_bound &
            )THEN
              CALL RANDOM_NUMBER( rand_num2 )
              rand_num= NINT( 5.0D0*( rand_num2 + 1.0D0 ) )
              npart_shelleq( r )= npart_shelleq( r )/rand_num
                                !  npart_shelleq( r - 1 ) &
                                !+ rel_sign*NINT( 1 + rand_num )
              upper_bound_tmp= upper_bound
              lower_bound_tmp= lower_bound
              npart_out= npart_out - npart_shell( r )/2
              cnt2= 1
              CYCLE
            ENDIF
            cnt2= 1
          ENDIF

          npart_out= npart_out - npart_shell( r )/2

          CALL RANDOM_NUMBER( rand_num )
          CALL RANDOM_NUMBER( rand_num2 )
          IF( rand_num2 < half )  rel_sign= - 1
          IF( rand_num2 >= half ) rel_sign=   1
          npart_shelleq( r )= CEILING( SQRT( &
                                2*(shell_masses( r )/m_parts( r - 1 )) &
                                /npart_shell_kept &
                              ) ) + rel_sign*NINT( 1 + rand_num )

          IF( npart_shelleq( r ) == 0 .OR. npart_shell( r ) == 0 )THEN
            CALL RANDOM_NUMBER( rand_num )
            CALL RANDOM_NUMBER( rand_num2 )
            IF( rand_num2 < half )  rel_sign= - 1
            IF( rand_num2 >= half ) rel_sign=   1
            npart_shelleq( r )= npart_shelleq( r - 1 ) &
                              + rel_sign*NINT( 1 + rand_num )
          ENDIF

          !PRINT *, "Right after correction of particle number"

          CYCLE

        ENDIF

        !IF( shell_radii( r )/radius > 0.85D0 )THEN
        !  npart_shelleq( r )= npart_shelleq( r )*5
        !ENDIF

      ENDIF

      !IF( r_cnt > 1 )THEN
      !  npart_test= 0
      !  DO itr= 1, r, 1
      !    npart_test= npart_test + npart_shell( itr )
      !  ENDDO
      !  IF( npart_test/2 /= npart_out )THEN
      !    PRINT *, "** ERROR! The sum of the particles on the shells is not ", &
      !             "equal to the total number of particles. Stopping.."
      !    PRINT *, "npart_test=", npart_test/2, ", npart_out=", npart_out
      !    PRINT *
      !    STOP
      !  ENDIF
      !ENDIF

      IF( debug ) PRINT *, "10"

      PRINT *, " * Placed", npart_shell( r )/2, &
               " particles on one emisphere of spherical shell ", r, &
               " out of ", n_shells
      IF( r /= CEILING(DBLE(n_shells)/2.0D0) ) PRINT *, "   Ratio of particle masses on last 2 shells: ", &
                           "   m_parts(", r, ")/m_parts(", prev_shell, ")= ",  &
                           m_parts( r )/m_parts( prev_shell )
      PRINT *, "   Shell radius= ", shell_radii( r )/radius*100.0D0, &
               "% of the radius of the star"
      PRINT *, "   Placed", npart_out, " particles overall, so far."
      PRINT *

      !IF( r == n_shells ) EXIT
      IF( r == n_shells )THEN
        r= CEILING(DBLE(n_shells)/2.0D0) - 1
        r_cnt= r_cnt + 1
        cnt2 = 0
        upper_bound_tmp= upper_bound
        lower_bound_tmp= lower_bound
        !PRINT *, "last shell"
      ELSEIF( r == 1 )THEN
        !PRINT *, "exit"
        EXIT
      ELSEIF( r < CEILING(DBLE(n_shells)/2.0D0) )THEN
        r= r - 1
        r_cnt= r_cnt + 1
        cnt2 = 0
        upper_bound_tmp= upper_bound
        lower_bound_tmp= lower_bound
        !PRINT *, "inner layers"
      ELSEIF( r >= CEILING(DBLE(n_shells)/2.0D0) )THEN
        r= r + 1
        r_cnt= r_cnt + 1
        cnt2 = 0
        upper_bound_tmp= upper_bound
        lower_bound_tmp= lower_bound
        !PRINT *, "outer layers"
      ENDIF

      IF( debug ) PRINT *, "12"

      !PRINT *, npart_shelleq( r ), (npart_shelleq( r )**2)/2, npart_shell( r )
      !PRINT *, colatitude_pos( r )% colatitudes

      !IF( r == 2 ) STOP

   !   r= r + 1
   !   r_cnt= r_cnt + 1
   !   cnt2 = 0
   !   upper_bound_tmp= upper_bound
   !   lower_bound_tmp= lower_bound

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

    ENDDO place_particles_on_shells
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
      PRINT *, "npart_test=", npart_test/2, ", npart_out=", npart_out
      PRINT *
      STOP
    ENDIF

    IF( ( 2*npart_out >= 1.225D0*npart_approx .OR. &
          2*npart_out <= 0.775D0*npart_approx ) .AND. find_npart )THEN
      r= 1
      r_cnt= 1
      cnt2 = 0
      upper_bound_tmp= upper_bound
      lower_bound_tmp= lower_bound
      npart_out= 0
    ELSE
      EXIT
    ENDIF

    IF( debug ) PRINT *, "13"

    !IF( find_npart )THEN
    !  IF( 2*npart_out >= 1.2D0*npart_approx )THEN
    !    !m_p= 1.1D0*m_p
    !    !upper_bound    = 1.1D0*upper_bound
    !    !lower_bound    = 1.1D0*lower_bound
    !    upper_bound_tmp= upper_bound
    !    lower_bound_tmp= lower_bound
    !  ELSEIF( 2*npart_out <= 0.8D0*npart_approx )THEN
    !    !m_p= 0.9D0*m_p
    !    !upper_bound    = 0.9D0*upper_bound
    !    !lower_bound    = 0.9D0*lower_bound
    !    upper_bound_tmp= upper_bound
    !    lower_bound_tmp= lower_bound
    !  ELSE
    !    EXIT
    !  ENDIF
    !  r= 1
    !  r_cnt= 1
    !  cnt2 = 0
    !  !upper_bound_tmp= upper_bound
    !  !lower_bound_tmp= lower_bound
    !  npart_out= 0
    !ELSE
    !  EXIT
    !ENDIF

    ENDDO find_desired_npart
    !STOP

    DEALLOCATE( pos_shell_tmp   )
    DEALLOCATE( g_xx_tmp        )
    DEALLOCATE( bar_density_tmp )
    DEALLOCATE( gam_euler_tmp   )
    DEALLOCATE( pvol_tmp        )

    IF( debug ) PRINT *, "14"

    PRINT *, "Mirroring particles..."

    PRINT *, "npart/2=", npart_out
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
        !PRINT *, pos_shells(r)% pos_shell( 1, 1:npart_shell( r ) )
        !PRINT *
        !PRINT *, pos_shells(r)% pos_shell( 2, 1:npart_shell( r ) )
        !PRINT *
        !PRINT *, pos_shells(r)% pos_shell( 3, 1:npart_shell( r ) )
        !PRINT *
        IF( pos_shells(r)% baryon_density( itr ) == 0 )THEN
          PRINT *, "When mirroring particles"
          PRINT *, r, itr, pos_shells(r)% pos_shell( 1, itr ), &
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
        !IF( pos_shells(r)% pvol_shell2( itr ) <= 0.0D0 )THEN
        !  PRINT *, "When computing particle volume"
        !  PRINT *, "pos_shells(", r, ")% pvol_shell2( ", itr, " ) =", &
        !           pos_shells(r)% pvol_shell2( itr )
        !  STOP
        !ENDIF

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
        !IF( pos_shells(r)% pvol_shell2( itr ) <= 0 )THEN
        !  PRINT *, "When computing shell volumes and masses"
        !  PRINT *, "pos_shells(", r, ")% pvol_shell2( ", itr, " ) =", &
        !           pos_shells(r)% pvol_shell2( itr )
        !  STOP
        !ENDIF
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

    !DO r= 1, n_shells, 1
    !  DO itr= 1, npart_shell( r ), 1
    !
    !    PRINT*, (m_parts( r )*MSun/amu)/pos_shells(r)% pvol_shell( itr ) &
    !            /(pos_shells(r)% g_xx( itr ) &
    !              *SQRT(pos_shells(r)% g_xx( itr )) &
    !              *pos_shells(r)% gamma_euler( itr )), &
    !            pos_shells(r)% baryon_density( itr )*MSun/amu
    !  ENDDO
    !ENDDO
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
        pvol( itr + npart_shell_tmp )  = pos_shells(r)% pvol_shell2( itr )
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

    IF( PRESENT(filename_shells_pos) )THEN
      finalnamefile= filename_shells_pos
    ELSE
      finalnamefile= "shells_pos.dat"
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

    !CALL OMP_SET_NUM_THREADS(80)

    PRINT *, " * SUBROUTINE place_particles_spherical_shells executed."
    PRINT *

    !STOP

  END PROCEDURE place_particles_spherical_shells


END SUBMODULE spherical_shells
