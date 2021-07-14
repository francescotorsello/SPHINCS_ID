! File:         submodule_particles_sph_variables.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (particles_id) particles_sph_variables

  !****************************************************
  !                                                   *
  ! Implementation of the methods of TYPE particles   *
  ! that compute, print and read the SPH variables.   *
  !                                                   *
  ! FT 16.10.2020                                     *
  !                                                   *
  ! Renamed from particles_methods to                 *
  ! particles_sph_variables upon improving modularity *
  !                                                   *
  ! FT 12.07.2021                                     *
  !                                                   *
  !****************************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE compute_and_export_SPH_variables

    !************************************************
    !                                               *
    ! Compute the SPH quantities from the LORENE    *
    ! ID, and export it to a binary file with       *
    ! write_SPHINCS_dump, and to a formatted file   *
    !                                               *
    ! FT 18.09.2020                                 *
    !                                               *
    !************************************************

    USE constants,           ONLY: km2cm, km2m, m2cm, g2kg, amu, MSun_geo, &
                                   third, kg2g, Msun
    USE units,               ONLY: set_units, umass
    USE matrix,              ONLY: determinant_4x4_matrix
    USE sph_variables,       ONLY: npart, &  ! particle number
                                   n1,    &  ! particle number for star 1
                                   n2,    &  ! particle number for star 2
                                   pos_u, &  ! particle positions
                                   vel_u, &  ! particle velocities in
                                             ! coordinate frame
                                   nlrf,  &  ! baryon number density in
                                             ! local rest frame
                                   ehat,  &  ! canonical energy per baryon
                                   nu,    &  ! canonical baryon number per
                                             ! particle
                                   Theta, &  ! Generalized Lorentz factor
                                   h,     &  ! Smoothing length
                                   Pr,    &  ! Pressure
                                   u,     &  ! Internal energy in local rest
                                             ! frame (no kinetic energy)
                                   temp,  &  ! Temperature
                                   av,    &  ! Dissipation
                                   ye,    &  ! Electron fraction
                                   divv,  &  ! Divergence of velocity vel_u
                                   allocate_SPH_memory, &
                                   deallocate_SPH_memory
    USE metric_on_particles, ONLY: allocate_metric_on_particles, &
                                   deallocate_metric_on_particles, &
                                   sq_det_g4
    USE options,             ONLY: basename
    USE input_output,        ONLY: dcount, write_SPHINCS_dump, read_options
    USE NR,                  ONLY: indexx

    USE RCB_tree_3D,         ONLY: allocate_RCB_tree_memory_3D,&
                                   deallocate_RCB_tree_memory_3D, iorig
    USE APM,                 ONLY: density_loop
    USE kernel_table,        ONLY: ktable
    USE options,             ONLY: ikernel, ndes
    USE set_h,               ONLY: exact_nei_tree_update
    USE gradient,            ONLY: allocate_gradient, deallocate_gradient
    USE sphincs_sph,         ONLY: density!, flag_dead_ll_cells
    USE alive_flag,          ONLY: alive
    USE APM,                 ONLY: assign_h

    IMPLICIT NONE

    ! The flag call_flag is set different than 0 if the SUBROUTINE
    ! compute_and_export_SPH_variables is called
    INTEGER, SAVE:: call_flag= 0

    ! Spacetime indices \mu and \nu
    INTEGER:: nus, mus, cnt1, cnt2, a

    DOUBLE PRECISION:: g4(0:3,0:3)
    DOUBLE PRECISION:: det,sq_g, Theta_a, temp1, temp2, nu_max1, nu_max2, &
                       nu_tmp, nu_thres1, nu_thres2
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: tmp
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: h_guess

    LOGICAL, PARAMETER:: debug= .FALSE.

    CHARACTER( LEN= : ), ALLOCATABLE:: compose_namefile
    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile


    PRINT *, "** Executing the compute_and_export_SPH_variables " &
             // "subroutine..."
    PRINT *

    !
    !-- Set up the MODULE variables in MODULE sph_variables
    !-- (used by write_SPHINCS_dump)
    !
    npart= THIS% npart
    n1= THIS% npart1
    n2= THIS% npart2

    CALL set_units('NSM')
    CALL read_options

    CALL allocate_SPH_memory
    CALL allocate_metric_on_particles( THIS% npart )

    IF( debug ) PRINT *, "1"

    IF(.NOT.ALLOCATED( THIS% nu ))THEN
      ALLOCATE( THIS% nu( THIS% npart ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array nu ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array nu" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% nlrf ))THEN
      ALLOCATE( THIS% nlrf( THIS% npart ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array nlrf ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array nlrf" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% Theta ))THEN
      ALLOCATE( THIS% Theta( THIS% npart ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array Theta ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array Theta" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% v ))THEN
      ALLOCATE( THIS% v( 0:3, THIS% npart ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array v ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array v" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% h ))THEN
      ALLOCATE( THIS% h( THIS% npart ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array h ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array h" )
    ENDIF

    IF( debug ) PRINT *, "2"

    !
    !-- Compute SPH quantities
    !
    CALL THIS% sph_computer_timer% start_timer()
    compute_SPH_variables_on_particles: DO itr= 1, THIS% npart, 1

      ! Particle positions [Msun_geo]
      pos_u(1,itr)= THIS% pos(1,itr)
      pos_u(2,itr)= THIS% pos(2,itr)
      pos_u(3,itr)= THIS% pos(3,itr)

      ! Coordinate velocity of the fluid [c]
      THIS% v(0,itr)= 1.0D0
      THIS% v(1,itr)= THIS% lapse_parts(itr)*THIS% v_euler_parts_x(itr) &
                    - THIS% shift_parts_x(itr)
      THIS% v(2,itr)= THIS% lapse_parts(itr)*THIS% v_euler_parts_y(itr) &
                    - THIS% shift_parts_y(itr)
      THIS% v(3,itr)= THIS% lapse_parts(itr)*THIS% v_euler_parts_z(itr) &
                    - THIS% shift_parts_z(itr)
      vel_u(1,itr)  = THIS% v(1,itr)
      vel_u(2,itr)  = THIS% v(2,itr)
      vel_u(3,itr)  = THIS% v(3,itr)

      !
      !-- Metric as matrix for easy manipulation
      !
      g4(0,0)= - THIS% lapse_parts(itr)**2 &
             + THIS% g_xx_parts(itr)*THIS% shift_parts_x(itr) &
              *THIS% shift_parts_x(itr)&
             + 2 * THIS% g_xy_parts(itr)*THIS% shift_parts_x(itr) &
              *THIS% shift_parts_y(itr)&
             + 2 * THIS% g_xz_parts(itr)*THIS% shift_parts_x(itr) &
              *THIS% shift_parts_z(itr)&
             + THIS% g_yy_parts(itr)*THIS% shift_parts_y(itr) &
              *THIS% shift_parts_y(itr)&
             + 2 * THIS% g_yz_parts(itr)*THIS% shift_parts_y(itr) &
              *THIS% shift_parts_z(itr)&
             + THIS% g_zz_parts(itr)*THIS% shift_parts_z(itr) &
              *THIS% shift_parts_z(itr)
      g4(0,1)= THIS% g_xx_parts(itr)*THIS% shift_parts_x(itr) &
             + THIS% g_xy_parts(itr)*THIS% shift_parts_y(itr) &
             + THIS% g_xz_parts(itr)*THIS% shift_parts_z(itr)
      g4(0,2)= THIS% g_xy_parts(itr)*THIS% shift_parts_x(itr) &
             + THIS% g_yy_parts(itr)*THIS% shift_parts_y(itr) &
             + THIS% g_yz_parts(itr)*THIS% shift_parts_z(itr)
      g4(0,3)= THIS% g_xz_parts(itr)*THIS% shift_parts_x(itr) &
             + THIS% g_yz_parts(itr)*THIS% shift_parts_y(itr) &
             + THIS% g_zz_parts(itr)*THIS% shift_parts_z(itr)

      g4(1,0)= THIS% g_xx_parts(itr)*THIS% shift_parts_x(itr) &
             + THIS% g_xy_parts(itr)*THIS% shift_parts_y(itr) &
             + THIS% g_xz_parts(itr)*THIS% shift_parts_z(itr)
      g4(1,1)= THIS% g_xx_parts(itr)
      g4(1,2)= THIS% g_xy_parts(itr)
      g4(1,3)= THIS% g_xz_parts(itr)

      g4(2,0)= THIS% g_xy_parts(itr)*THIS% shift_parts_x(itr) &
             + THIS% g_yy_parts(itr)*THIS% shift_parts_y(itr) &
             + THIS% g_yz_parts(itr)*THIS% shift_parts_z(itr)
      g4(2,1)= THIS% g_xy_parts(itr)
      g4(2,2)= THIS% g_yy_parts(itr)
      g4(2,3)= THIS% g_yz_parts(itr)

      g4(3,0)= THIS% g_xz_parts(itr)*THIS% shift_parts_x(itr) &
             + THIS% g_yz_parts(itr)*THIS% shift_parts_y(itr) &
             + THIS% g_zz_parts(itr)*THIS% shift_parts_z(itr)
      g4(3,1)= THIS% g_xz_parts(itr)
      g4(3,2)= THIS% g_yz_parts(itr)
      g4(3,3)= THIS% g_zz_parts(itr)

      ! sqrt(-det(g4))
      CALL determinant_4x4_matrix(g4,det)
      IF( ABS(det) < 1D-10 )THEN
          PRINT *, "The determinant of the spacetime metric is " &
                   // "effectively 0 at particle ", itr
          STOP
      ELSEIF( det > 0 )THEN
          PRINT *, "The determinant of the spacetime metric is " &
                   // "positive at particle ", itr
          STOP
      ENDIF
      sq_g= SQRT(-det)
      sq_det_g4(itr)= sq_g

      !
      !-- Generalized Lorentz factor
      !
      Theta_a= 0.D0
      DO nus=0,3
        DO mus=0,3
          Theta_a= Theta_a &
                   + g4(mus,nus)*THIS% v(mus,itr)*THIS% v(nus,itr)
        ENDDO
      ENDDO
      Theta_a= 1.0D0/SQRT(-Theta_a)
      Theta(itr)= Theta_a
      THIS% Theta(itr)= Theta_a

      ! This is a first guess for the smoothing length
      ! The real smoothing length is such that the kernel sees
      ! only 300 neighbors, and it is computed in exact_nei_tree_update,
      ! MODULE set_h.
      ! N.B. This first guess is very important, as it affects the algorithm in
      !      exact_tree_nei_update. Here it is set to 3 times the size of a
      !      particle
      !IF( itr <= THIS% npart1 )THEN
      !  THIS% h(itr)= 3.0*(THIS% vol1_a)**third
      !ELSE
      !  THIS% h(itr)= 3.0*(THIS% vol2_a)**third
      !ENDIF

      IF( debug ) PRINT *, "3"

      IF( .NOT.ALLOCATED( THIS% pvol ) )THEN
        PRINT *, "** ERROR! The array pvol is not allocated. ", &
                 "Stopping..."
        PRINT *
        STOP
      ENDIF

      IF( .NOT.THIS% apm_iterate )THEN

        THIS% h(itr)= 3.0D0*(THIS% pvol(itr))**third
        h(itr)= THIS% h(itr)
        ! /(Msun_geo**3)
      ENDIF

      ! Baryon density in the local rest frame [baryon (Msun_geo)^{-3}]
      ! Computed from the LORENE baryon mass density in [kg/m^3]
      nlrf(itr)= THIS% baryon_density_parts(itr)*((Msun_geo*km2m)**3)/(amu*g2kg)
      THIS% nlrf(itr)= nlrf(itr)

      THIS% nstar( itr )= THIS% nlrf(itr)*Theta_a*sq_g
      THIS% pmass( itr )= THIS% nstar( itr )*THIS% pvol( itr )
      THIS% particle_density( itr )= (THIS% nstar( itr ))/( THIS% pmass(itr) )

      ! Internal energy per baryon (specific internal energy)
      ! In module_TOV.f90, this quantity is computed as follows,
      ! eps= pressure/(Gammap-1.0D0)/density
      !    = energy_density/mass_density
      ! In LORENE, the specific energy is computed as
      ! e_spec= (energy_density - mass_density)/mass_density

      ! Specific internal energy [c^2]
      u(itr)= THIS% specific_energy_parts(itr)

      ! Pressure [amu*c**2/(Msun_geo**3)]
      !          dimensions: [(M/L**3)*L**2/T**2]= [M/(L*T**2)], same as
      !                      energy density
      Pr(itr)= THIS% pressure_parts(itr)*((Msun_geo*km2m)**3)/(amu*g2kg)
      THIS% pressure_parts_cu(itr)= Pr(itr)

      ! Baryon per particle
      !IF( itr <= THIS% npart1 )THEN
      !  nu(itr)= nlrf(itr)*THIS% vol1_a*Theta_a*sq_g
      !  THIS% nbar1= THIS% nbar1 + nu(itr)
      !ELSE
      !  nu(itr)= nlrf(itr)*THIS% vol2_a*Theta_a*sq_g
      !  THIS% nbar2= THIS% nbar2 + nu(itr)
      !ENDIF
      !THIS% nu(itr)= nu(itr)
      !THIS% nbar_tot= THIS% nbar_tot + THIS% nu(itr)

      ! THIS% mass1*umass converts mass1 in g; umass is Msun in g
      ! amu is the atomic mass unit in g
      !temp1= THIS% mass1*umass/amu/THIS% npart1
      !temp2= THIS% mass2*umass/amu/THIS% npart2
      !IF( nu(itr) /= temp1 )THEN
      !  PRINT *, "The two formulas to compute the baryons per particle ", &
      !           "do not give consistent results at particle ", itr, "."
      !  PRINT *, "nu(itr)*THIS% vol_a=", nu(itr)
      !  PRINT *, "THIS% mass1*umass/THIS% npart1/amu=", temp1
      !  PRINT *, "THIS% mass2*umass/THIS% npart2/amu=", temp2
      !  PRINT *
      !  STOP
      !ENDIF

      IF( debug )THEN
        IF( itr >= THIS% npart/10 - 200 .AND. itr <= THIS% npart/10 )THEN
          PRINT "(A15,E15.4)", "det=", det
          PRINT "(A15,E15.4)", "sq_g=", sq_g
          PRINT "(A15,E15.4)", "sq_det_g4(itr)=", sq_det_g4(itr)
          PRINT "(A15,E15.4)", "amu=", amu
          PRINT "(A15,E15.4)", "g2kg=", g2kg
          PRINT "(A15,E15.4)", "(1477**3)=", (1477.0)**3
          PRINT *
          PRINT "(A15,E15.4)", "nbar(a)=", THIS% baryon_density_parts(itr)
          PRINT "(A25,E15.4)", "THIS% nlrf(a)=", THIS% nlrf(itr)
          PRINT "(A25,E15.4)", "THIS% nu(a)=", THIS% nu(itr)
          PRINT "(A25,E15.4)", "THIS% pressure_parts_cu(a)=", &
                                        THIS% pressure_parts_cu(itr)
          PRINT *
          PRINT "(A15,E15.4)", "theta(a)=", THIS% Theta(itr)
          PRINT "(A15,E15.4)", "sq_det_g4(a)=", sq_det_g4(itr)
          PRINT "(A15,E15.4)", "vol_a=", THIS% vol_a
          PRINT *
          PRINT *
        ENDIF
      ENDIF

      ! Temperature: here dummy
      temp(itr)=  1.0D0

      ! Dissipation parameter
      av(itr)=    1.0D0

      ! Velocity divergence
      divv(itr)=  0.D0

    ENDDO compute_SPH_variables_on_particles
    !THIS% nbar_tot= THIS% nbar1 + THIS% nbar2

    IF( debug ) PRINT *, "1"

    !---------------------------------------------------------------------!
    !--  Assignment of nu on the stars, with the purpose                --!
    !--  of having a more uniform nu over the particles without losing  --!
    !--  baryon mass.                                                   --!
    !---------------------------------------------------------------------!

    IF( THIS% redistribute_nu )THEN

      IF( THIS% distribution_id == 3 )THEN
        PRINT *, "** ERROR! Particle placer ", THIS% distribution_id, &
                 " is not compatible with redistribute_nu= .TRUE."
        PRINT *, " * Check the parameter file lorene_bns_id_particles.par. ", &
                 "Stopping..."
        PRINT *
        STOP
      ENDIF

      ! Index particles on star 1 in increasing order of nu

      !CALL indexx( THIS% npart1, &
      !             THIS% baryon_density_parts( 1 : THIS% npart1 ), &
      !             !THIS% nu( 1 : THIS% npart1 ), &
      !             THIS% baryon_density_index( 1 : THIS% npart1 ) )

      ! Index particles on star 2 in increasing order of nu

      !CALL indexx( THIS% npart2, &
      !             THIS% baryon_density_parts( THIS% npart1 + 1 : &
      !                                         THIS% npart ), &
      !             !THIS% nu( THIS% npart1 + 1 : THIS% npart ), &
      !             THIS% baryon_density_index( THIS% npart1 + 1 : &
      !                                         THIS% npart ) )

      ! Shift indices on star 2 by npart1 since all the arrays store
      ! the quantities on star 1 first, and then on star 2

      !THIS% baryon_density_index( THIS% npart1 + 1 : THIS% npart )= &
      !            THIS% npart1 + &
      !            THIS% baryon_density_index( THIS% npart1 + 1 : THIS% npart )

      ! Store the maximum values of nu on the two stars, and the thresholds

      !nu_max1= MAXVAL( THIS% nu( 1 : THIS% npart1 ) )
      !nu_max2= MAXVAL( THIS% nu( THIS% npart1 + 1 : THIS% npart ) )

      nu_max1= nlrf( THIS% baryon_density_index( THIS% npart1 ) )&
              *THIS% pvol( THIS% npart1 ) &
              *Theta( THIS% baryon_density_index( THIS% npart1 ) )&
              *sq_det_g4( THIS% baryon_density_index( THIS% npart1 ) )
      nu_max2= nlrf( THIS% baryon_density_index( THIS% npart ) )&
              *THIS% pvol( THIS% npart ) &
              *Theta( THIS% baryon_density_index( THIS% npart ) )&
              *sq_det_g4( THIS% baryon_density_index( THIS% npart ) )

      nu_thres1= nu_max1/THIS% nu_ratio
      nu_thres2= nu_max2/THIS% nu_ratio

      ! Reset the total baryon number to 0 (necessary), and nu to an arbitrary
      ! value (to make debugging easier)

      nu= 1.0D0
      THIS% nu= 1.0D0
      THIS% nbar_tot= 0.0D0
      THIS% nbar1= 0.0D0
      THIS% nbar2= 0.0D0

      cnt1= 0
      compute_nu_on_particles_star1: DO itr= THIS% npart1, 1, -1

        cnt1= cnt1 + 1

        nu_tmp= nlrf( THIS% baryon_density_index( itr ) ) &
                *THIS% pvol(itr) &
                *Theta( THIS% baryon_density_index( itr ) )&
                *sq_det_g4( THIS% baryon_density_index( itr ) )

        !IF( itr == THIS% npart1 ) nu_max= nu_tmp ! move this out of the loop

        IF( nu_tmp > nu_thres1 )THEN
          nu( THIS% baryon_density_index( itr ) )      = nu_tmp
          THIS% nu( THIS% baryon_density_index( itr ) )= nu_tmp
        ELSE
          nu( THIS% baryon_density_index( itr ) )      = nu_thres1
          THIS% nu( THIS% baryon_density_index( itr ) )= nu_thres1
        ENDIF

        THIS% nbar1= THIS% nbar1 + &
                     THIS% nu( THIS% baryon_density_index( itr ) )

        IF( THIS% nbar1*amu/MSun > THIS% mass1 )THEN
          EXIT
        ENDIF

      ENDDO compute_nu_on_particles_star1

      cnt2= 0
      compute_nu_on_particles_star2: DO itr= THIS% npart, THIS% npart1 + 1, -1

        cnt2= cnt2 + 1

        nu_tmp= nlrf( THIS% baryon_density_index( itr ) ) &
                *THIS% pvol(itr) &
                *Theta( THIS% baryon_density_index( itr ) ) &
                *sq_det_g4( THIS% baryon_density_index( itr ) )

        !IF( itr == THIS% npart ) nu_max= nu_tmp

        IF( nu_tmp > nu_thres2 )THEN
          nu( THIS% baryon_density_index( itr ) )      = nu_tmp
          THIS% nu( THIS% baryon_density_index( itr ) )= nu_tmp
        ELSE
          nu( THIS% baryon_density_index( itr ) )      = nu_thres2
          THIS% nu( THIS% baryon_density_index( itr ) )= nu_thres2
        ENDIF

        THIS% nbar2= THIS% nbar2 + &
                     THIS% nu( THIS% baryon_density_index( itr ) )

        IF( THIS% nbar2*amu/MSun > THIS% mass2 )THEN
          EXIT
        ENDIF

      ENDDO compute_nu_on_particles_star2
      THIS% nbar_tot= THIS% nbar1 + THIS% nbar2
      !PRINT *, cnt1, cnt2
      !PRINT *, nu_thres1, nu_thres2, nu_tmp
      !STOP

      !
      !-- Reshape MODULE variables
      !

      CALL THIS% reshape_sph_field( pos_u, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( vel_u, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( Theta, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( h, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( nlrf, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( u, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( Pr, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( nu, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( temp, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( av, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( divv, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      !
      !-- Reshape TYPE member SPH variables
      !

      CALL THIS% reshape_sph_field( THIS% pos, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% v, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% v_euler_parts_x, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% v_euler_parts_y, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% v_euler_parts_z, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% Theta, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% h, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% baryon_density_parts, cnt1, &
                                    cnt2, THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% nlrf, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% energy_density_parts, cnt1, &
                                    cnt2, THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% specific_energy_parts, cnt1, &
                                    cnt2, THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% pressure_parts, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% pressure_parts_cu, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% nu, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% pvol, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      !
      !-- Reshape TYPE member spacetime variables
      !

      CALL THIS% reshape_sph_field( THIS% lapse_parts, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% shift_parts_x, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% shift_parts_y, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% shift_parts_z, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% g_xx_parts, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% g_xy_parts, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% g_xz_parts, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% g_yy_parts, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% g_yz_parts, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      CALL THIS% reshape_sph_field( THIS% g_zz_parts, cnt1, cnt2, &
                                    THIS% baryon_density_index )

      !
      !-- Reassign particle numbers
      !

      npart= cnt1 + cnt2
      THIS% npart= npart
      THIS% npart1= cnt1
      THIS% npart2= cnt2
      n1= THIS% npart1
      n2= THIS% npart2

      PRINT *, " * Particles replaced after reassigning nu."
      PRINT *, " * New number of particles=", THIS% npart
      PRINT *
      PRINT *, " * Number of particles on NS 1=", THIS% npart1
      PRINT *, " * Number of particles on NS 2=", THIS% npart2
      PRINT *

    ELSEIF( THIS% distribution_id == 0 .AND. THIS% read_nu )THEN

      ! Do nothing, nu is already read fom file
      nu= THIS% nu
      THIS% nbar1= SUM( THIS% nu(1:THIS% npart1), DIM= 1 )
      THIS% nbar2= SUM( THIS% nu(THIS% npart1+1:THIS% npart), DIM= 1 )
      THIS% nbar_tot= THIS% nbar1 + THIS% nbar2

    ELSEIF( THIS% apm_iterate )THEN

      ! Do nothing, nu and h are already computed in the APM iteration
      nu= THIS% nu
      h = THIS% h
      THIS% nbar1= SUM( THIS% nu(1:THIS% npart1), DIM= 1 )
      THIS% nbar2= SUM( THIS% nu(THIS% npart1+1:THIS% npart), DIM= 1 )
      THIS% nbar_tot= THIS% nbar1 + THIS% nbar2

    ELSE

      IF( .FALSE. .AND. THIS% distribution_id == 3 )THEN

        DO itr= 1, THIS% npart1, 1
          nu(itr)= THIS% pmass( itr )*MSun/amu
          THIS% nu(itr)= nu(itr)
          THIS% nbar1= THIS% nbar1 + nu(itr)
        ENDDO
        DO itr= THIS% npart1 + 1, THIS% npart, 1
          nu(itr)= THIS% pmass( itr )*MSun/amu
          THIS% nu(itr)= nu(itr)
          THIS% nbar2= THIS% nbar2 + nu(itr)
        ENDDO
        THIS% nbar_tot= THIS% nbar1 + THIS% nbar2

      ELSE

        DO itr= 1, THIS% npart1, 1
          nu(itr)= nlrf(itr)*THIS% pvol(itr)*Theta( itr )*sq_det_g4( itr )
          THIS% nu(itr)= nu(itr)
          THIS% nbar1= THIS% nbar1 + nu(itr)
        ENDDO
        DO itr= THIS% npart1 + 1, THIS% npart, 1
          nu(itr)= nlrf(itr)*THIS% pvol(itr)*Theta( itr )*sq_det_g4( itr )
          THIS% nu(itr)= nu(itr)
          THIS% nbar2= THIS% nbar2 + nu(itr)
        ENDDO
        THIS% nbar_tot= THIS% nbar1 + THIS% nbar2

      ENDIF

    ENDIF

    IF(.NOT.ALLOCATED( THIS% Ye ))THEN
      ALLOCATE( THIS% Ye( THIS% npart ), STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array Ye ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array lapse_parts" )
    ENDIF

    IF( THIS% compose_eos )THEN
      compose_namefile= TRIM(THIS% compose_path)//TRIM(THIS% compose_filename)
      CALL THIS% read_compose_composition( compose_namefile )
      CALL THIS% compute_Ye()
    ENDIF

    assign_ye_on_particles: DO itr= 1, THIS% npart, 1

      ! Electron fraction
      IF( THIS% compose_eos )THEN
        Ye(itr)= THIS% Ye(itr)
      ELSE
        Ye(itr)= 0.0D0
        THIS% Ye(itr)= 0.0D0
      ENDIF

    ENDDO assign_ye_on_particles
    CALL THIS% sph_computer_timer% stop_timer()

    !
    !-- Printouts
    !
    !"(A28,E15.8,A10)"
    PRINT *, " * Maximum baryon density on star 1= ", &
              MAXVAL( THIS% baryon_density_parts(1:THIS% npart1), DIM= 1 ) &
              *kg2g/(m2cm**3), " g cm^{-3}"
    PRINT *, " * Minimum baryon density= on star 1 ", &
              MINVAL( THIS% baryon_density_parts(1:THIS% npart1), DIM= 1 ) &
              *kg2g/(m2cm**3), " g cm^{-3}"
    PRINT *, " * Ratio between the two=", &
             MAXVAL( THIS% baryon_density_parts(1:THIS% npart1), DIM= 1 )/ &
             MINVAL( THIS% baryon_density_parts(1:THIS% npart1), DIM= 1 )
    PRINT *
    PRINT *, " * Maximum baryon density on star 2= ", &
     MAXVAL( THIS% baryon_density_parts(THIS% npart1+1:THIS% npart), DIM= 1 ) &
     *kg2g/(m2cm**3), " g cm^{-3}"
    PRINT *, " * Minimum baryon density on star 2= ", &
     MINVAL( THIS% baryon_density_parts(THIS% npart1+1:THIS% npart), DIM= 1 ) &
     *kg2g/(m2cm**3), " g cm^{-3}"
    PRINT *, " * Ratio between the two=", &
     MAXVAL( THIS% baryon_density_parts(THIS% npart1+1:THIS% npart), DIM= 1 )/ &
     MINVAL( THIS% baryon_density_parts(THIS% npart1+1:THIS% npart), DIM= 1 )
    PRINT *

    PRINT *, " * Maximum nlrf on star 1=", &
             MAXVAL( THIS% nlrf(1:THIS% npart1), DIM= 1 ), &
             "baryon Msun_geo^{-3}"
    PRINT *, " * Minimum nlrf on star 1=", &
             MINVAL( THIS% nlrf(1:THIS% npart1), DIM= 1 ), &
             "baryon Msun_geo^{-3}"
    PRINT *, " * Ratio between the two=", &
             MAXVAL( THIS% nlrf(1:THIS% npart1), DIM= 1 )/ &
             MINVAL( THIS% nlrf(1:THIS% npart1), DIM= 1 )
    PRINT *
    PRINT *, " * Maximum nlrf on star 2=", &
             MAXVAL( THIS% nlrf(THIS% npart1+1:THIS% npart), DIM= 1 ), &
             "baryon Msun_geo^{-3}"
    PRINT *, " * Minimum nlrf on star 2=", &
             MINVAL( THIS% nlrf(THIS% npart1+1:THIS% npart), DIM= 1 ), &
             "baryon Msun_geo^{-3}"
    PRINT *, " * Ratio between the two=", &
             MAXVAL( THIS% nlrf(THIS% npart1+1:THIS% npart), DIM= 1 )/ &
             MINVAL( THIS% nlrf(THIS% npart1+1:THIS% npart), DIM= 1 )
    PRINT *

    THIS% nuratio1= MAXVAL( THIS% nu(1:THIS% npart1), DIM= 1 )/ &
                    MINVAL( THIS% nu(1:THIS% npart1), DIM= 1 )
    PRINT *, " * Maximum n. baryon per particle on star 1 (nu)=", &
             MAXVAL( THIS% nu(1:THIS% npart1), DIM= 1 )
    PRINT *, " * Minimum n. baryon per particle on star 1 (nu)=", &
             MINVAL( THIS% nu(1:THIS% npart1), DIM= 1 )
    PRINT *, " * Ratio between the two=", THIS% nuratio1
    PRINT *
    THIS% nuratio2= MAXVAL( THIS% nu(THIS% npart1+1:THIS% npart), DIM= 1 ) &
                    /MINVAL( THIS% nu(THIS% npart1+1:THIS% npart), DIM= 1 )
    PRINT *, " * Maximum n. baryon per particle on star 2 (nu)=", &
             MAXVAL( THIS% nu(THIS% npart1+1:THIS% npart), DIM= 1 )
    PRINT *, " * Minimum n. baryon per particle on star 2 (nu)=", &
             MINVAL( THIS% nu(THIS% npart1+1:THIS% npart), DIM= 1 )
    PRINT *, " * Ratio between the two=", THIS% nuratio2
    PRINT *
    THIS% nuratio= MAXVAL( THIS% nu, DIM= 1 )/MINVAL( THIS% nu, DIM= 1 )
    PRINT *, " * Baryon number ratio across the stars=", THIS% nuratio
    PRINT *
    PRINT *, " * Number of baryons on star 1=", THIS% nbar1
    PRINT *, " * Total mass of the baryons on star 1=", &
             THIS% nbar1*amu/Msun, "Msun =", &
             THIS% nbar1*amu/Msun/THIS% mass1, &
             "of the LORENE baryon mass of star 1."
    PRINT *, " * Number of baryons on star 2=", THIS% nbar2
    PRINT *, " * Total mass of the baryons on star 2=", &
             THIS% nbar2*amu/Msun, "Msun =", &
             THIS% nbar2*amu/Msun/THIS% mass2, &
             "of the LORENE baryon mass of star 2."
    PRINT *, " * Total number of baryons=", THIS% nbar_tot
    PRINT *, " * Total mass of the baryons=", &
             THIS% nbar_tot*amu/Msun, "Msun =", &
             THIS% nbar_tot*amu/Msun/(THIS% mass1 + THIS% mass2), &
             "of the total LORENE baryon mass."
    PRINT *

    !
    !-- Adjusting the baryon number per particle uniformly so that
    !-- the baryon mass is correct, but the ratio between nu_max and nu_min
    !-- does not change.
    !-- nlrf is not to be rescaled, according to Stephan, since:
    !--   (i)  it is directly computed from the LORENE ID and should therefore
    !--        be consistent with it
    !--   (ii) it is anyway immediately recomputed in SPHINCS_BSSN
    !
    IF( THIS% correct_nu )THEN
      !THIS% nlrf=
      !        THIS% nlrf/(THIS% nbar_tot*amu/Msun/(THIS% mass1 + THIS% mass2))
      !nlrf= nlrf/(THIS% nbar_tot*amu/Msun/(THIS% mass1 + THIS% mass2))
      !THIS% nu= THIS% nu/(THIS% nbar_tot*amu/Msun/(THIS% mass1 + THIS% mass2))
      !nu= nu/(THIS% nbar_tot*amu/Msun/(THIS% mass1 + THIS% mass2))
      !THIS% nbar_tot= &
      !    THIS% nbar_tot/(THIS% nbar_tot*amu/Msun/(THIS% mass1 + THIS% mass2))
      THIS% nu( 1:THIS% npart1 )= &
                THIS% nu( 1:THIS% npart1 )/(THIS% nbar1*amu/Msun/THIS% mass1)
      nu( 1:THIS% npart1 )= &
                nu( 1:THIS% npart1 )/(THIS% nbar1*amu/Msun/THIS% mass1)
      THIS% nbar1= THIS% nbar1/(THIS% nbar1*amu/Msun/THIS% mass1)

      THIS% nu( THIS% npart1 + 1:THIS% npart )= &
                THIS% nu( THIS% npart1 + 1:THIS% npart ) &
                /(THIS% nbar2*amu/Msun/THIS% mass2)
      nu( THIS% npart1 + 1:THIS% npart )= &
                nu( THIS% npart1 + 1:THIS% npart ) &
                /(THIS% nbar2*amu/Msun/THIS% mass2)
      THIS% nbar2= THIS% nbar2/(THIS% nbar2*amu/Msun/THIS% mass2)
      THIS% nbar_tot= THIS% nbar1 + THIS% nbar2

      PRINT *, " * Number of corrected baryons on star 1=", THIS% nbar1
      PRINT *, " * Total mass of the corrected baryons on star 1=", &
               THIS% nbar1*amu/Msun, "Msun =", &
               THIS% nbar1*amu/Msun/THIS% mass1, &
               "of the LORENE baryon mass of star 1."
      PRINT *, " * Number of corrected baryons on star 2=", THIS% nbar2
      PRINT *, " * Total mass of the corrected baryons on star 2=", &
               THIS% nbar2*amu/Msun, "Msun =", &
               THIS% nbar2*amu/Msun/THIS% mass2, &
               "of the LORENE baryon mass of star 2."
      PRINT *, " * Total number of corrected baryons=", THIS% nbar_tot
      PRINT *, " * Total mass of the corrected baryons=", &
               THIS% nbar_tot*amu/Msun, "Msun =", &
               THIS% nbar_tot*amu/Msun/(THIS% mass1 + THIS% mass2), &
               "of the total LORENE baryon mass."
      PRINT *
    ENDIF

    !
    !-- Exporting the SPH ID to a binary file, for evolution
    !
    IF( THIS% export_bin )THEN

      IF( PRESENT(namefile) )THEN

        finalnamefile= TRIM( namefile ) // "00000"
        dcount= -1 ! since it is increased before writing
        CALL write_SPHINCS_dump( finalnamefile )

      ELSE

        basename= "NSNS."
        dcount= -1 ! since it is increased before writing
        CALL write_SPHINCS_dump()

      ENDIF

    ENDIF

    CALL allocate_RCB_tree_memory_3D(npart)
    iorig(1:npart)= (/ (a,a=1,npart) /)

    !IF( call_flag == 0 )THEN
    !  ! tabulate kernel, get ndes
    !  CALL ktable(ikernel,ndes)
    !ENDIF

    ! flag that particles are 'alive'
    !ALLOCATE( alive( npart ) )
    !alive( 1:npart )= 1

    CALL allocate_gradient( npart )

    IF( debug ) PRINT *, "-1"

    PRINT *, " * Computing neighbours..."
    PRINT *
    CALL exact_nei_tree_update( ndes,        &
                                THIS% npart, &
                                THIS% pos,   &
                                THIS% nu )

    THIS% h= h

    !STOP

    PRINT *, " * Computing SPH density by kernel interpolation..."
    PRINT *
    ! density calls dens_ll_cell, which computes nstar on particle a as
    ! Na=     Na + nu(b)*Wab_ha, so this is nstar= nlrf*sq_g*Theta
    ! It has to be compared with particle_density_int= nlrf*sq_g*Theta
    CALL density( THIS% npart, &
                  THIS% pos,   &
                  THIS% nstar_int )

    PRINT *, " * Computing particle number density by kernel interpolation..."
    PRINT *
    nu= 1.0D0
    CALL density_loop( THIS% npart, THIS% pos, nu, h, &
                       THIS% particle_density_int )

    PRINT *, " * Deallocating MODULE variables..."
    PRINT *
    CALL deallocate_metric_on_particles
    CALL deallocate_gradient
    DEALLOCATE( alive )
    CALL deallocate_RCB_tree_memory_3D
    CALL deallocate_SPH_memory

    call_flag= call_flag + 1
    THIS% call_flag= call_flag

    PRINT *, "** Subroutine compute_and_export_SPH_variables executed."
    PRINT *

  END PROCEDURE compute_and_export_SPH_variables


  MODULE PROCEDURE read_sphincs_dump_print_formatted

    !************************************************
    !                                               *
    ! Read the SPH ID from the binary file output   *
    ! by write_SPHINCS_dump, and print it to a      *
    ! formatted file                                *
    !                                               *
    ! FT 12.02.2021                                 *
    !                                               *
    !************************************************

    USE sph_variables,       ONLY: npart, &  ! particle number
                                   pos_u, &  ! particle positions
                                   vel_u, &  ! particle velocities in
                                             ! coordinate frame
                                   nlrf,  &  ! baryon number density in
                                             ! local rest frame
                                   ehat,  &  ! canonical energy per baryon
                                   nu,    &  ! canonical baryon number per
                                             ! particle
                                   Theta, &  ! Generalized Lorentz factor
                                   h,     &  ! Smoothing length
                                   Pr,    &  ! Pressure
                                   u,     &  ! Internal energy in local rest
                                             ! frame (no kinetic energy)
                                   temp,  &  ! Temperature
                                   av,    &  ! Dissipation
                                   ye,    &  ! Electron fraction
                                   divv,  &  ! Divergence of velocity vel_u
                                   allocate_SPH_memory, &
                                   deallocate_SPH_memory
    USE metric_on_particles, ONLY: allocate_metric_on_particles, &
                                   deallocate_metric_on_particles, &
                                   sq_det_g4
    USE options,             ONLY: basename
    USE input_output,        ONLY: set_units, dcount, read_SPHINCS_dump

    IMPLICIT NONE

    INTEGER:: itr, min_y_index

    DOUBLE PRECISION:: min_abs_y, min_abs_z1, min_abs_z2
    DOUBLE PRECISION, DIMENSION( :, : ), ALLOCATABLE:: abs_pos

    LOGICAL:: exist

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    PRINT *, "** Executing the read_sphincs_dump_print_formatted subroutine..."

    !
    !-- Set up the MODULE variables in MODULE sph_variables
    !-- (used by write_SPHINCS_dump)
    !
    npart= THIS% npart

    CALL set_units('NSM')

    CALL allocate_SPH_memory
    CALL allocate_metric_on_particles( THIS% npart )

    finalnamefile= TRIM(namefile_bin)//"00000"
    CALL read_SPHINCS_dump( finalnamefile )

    ! Being abs_grid a local array, it is good practice to allocate it on the
    ! heap, otherwise it will be stored on the stack which has a very limited
    ! size. This results in a segmentation fault.
    ALLOCATE( abs_pos( 3, THIS% npart ) )

    IF( THIS% call_flag == 0 )THEN
      PRINT *, "** The SUBROUTINE print_formatted_lorene_id_particles must", &
               " be called after compute_and_export_SPH_variables, otherwise", &
               " there are no SPH fields to export to the formatted file."
      PRINT *, "   Aborting."
      PRINT *
      STOP
    ENDIF

    IF( PRESENT(namefile) )THEN
      finalnamefile= namefile
    ELSE
      finalnamefile= "sph_vars.dat"
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
    !CALL test_status( ios, err_msg, "...error when opening " &
    !                  // TRIM(finalnamefile) )

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Values of the fields (including coordinates) exported by LORENE "&
    // "on each grid point"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 1 in " // TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 1 in "&
    !        // TRIM(finalnamefile) )

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# column:      1        2       3       4       5", &
    "       6       7       8", &
    "       9       10      11", &
    "       12      13      14", &
    "       15      16      17"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 2 in " // TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 2 in "&
    !        // TRIM(finalnamefile) )

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "#      grid point      x [km]       y [km]       z [km]       lapse", &
    "       shift_x [c]    shift_y [c]    shift_z [c]", &
    "       baryon density in the local rest frame [kg m^{-3}$]", &
    "       energy density [c^2]", &
    "       specific energy [c^2]", &
    "       pressure [Pa]", &
    "       fluid 3-velocity wrt the Eulerian observer (3 columns) [c]", &
    "       fluid coordinate 3-velocity vel_u (3 columns) [c]", &
    "       baryon number per particle nu", &
    "       baryon density in the local rest frame nlrf [baryon/cm^3]", &
    "       generalized Lorentz factor Theta"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 3 in " // TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 3 in "&
    !        // TRIM(finalnamefile) )

    DO itr = 1, THIS% npart, 1
      abs_pos( 1, itr )= ABS( THIS% pos( 1, itr ) )
      abs_pos( 2, itr )= ABS( THIS% pos( 2, itr ) )
      abs_pos( 3, itr )= ABS( THIS% pos( 3, itr ) )
    ENDDO

    min_y_index= 0
    min_abs_y= 1D+20
    DO itr = 1, THIS% npart, 1
      IF( ABS( THIS% pos( 2, itr ) ) < min_abs_y )THEN
        min_abs_y= ABS( THIS% pos( 2, itr ) )
        min_y_index= itr
      ENDIF
    ENDDO

    min_abs_z1= MINVAL( abs_pos( 3, 1:THIS% npart1 ) )
    min_abs_z2= MINVAL( abs_pos( 3, THIS% npart1+1:THIS% npart ) )

    write_data_loop: DO itr = 1, THIS% npart, 1

      IF( THIS% export_form_xy .AND. &
          !( THIS% pos( 3, itr ) /= min_abs_z1 .AND. &
          !  THIS% pos( 3, itr ) /= min_abs_z2 ) &
          ( THIS% pos( 3, itr ) >=  0.5D0 .OR. &
            THIS% pos( 3, itr ) <= -0.5D0 ) &
      )THEN
        CYCLE
      ENDIF
      IF( THIS% export_form_x .AND. &
          !( THIS% pos( 3, itr ) /= min_abs_z1 &
          !.OR. THIS% pos( 3, itr ) /= min_abs_z2 &
          !.OR. THIS% pos( 2, itr ) /= THIS% pos( 2, min_y_index ) ) )THEN
          ( THIS% pos( 3, itr ) >=  0.5D0 .OR. &
            THIS% pos( 3, itr ) <= -0.5D0 .OR. &
            THIS% pos( 2, itr ) >=  0.5D0 .OR. &
            THIS% pos( 2, itr ) <= -0.5D0 ) &
      )THEN
        CYCLE
      ENDIF
      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        itr, &
        pos_u( 1, itr ), &
        pos_u( 2, itr ), &
        pos_u( 3, itr ), &
        vel_u( 1, itr ), &
        vel_u( 2, itr ), &
        vel_u( 3, itr ), &
        h( itr ), &
        u( itr ), &
        nu( itr ), &
        nlrf( itr ), &
        temp( itr ), &
        av( itr ), &
        ye( itr ), &
        divv( itr ), &
        Theta( itr ), &
        Pr( itr )

      IF( ios > 0 )THEN
        PRINT *, "...error when writing the arrays in " &
                 // TRIM(finalnamefile), ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when writing " &
      !       // "the arrays in " // TRIM(finalnamefile) )
    ENDDO write_data_loop

    CLOSE( UNIT= 2 )

    !
    !-- Deallocate MODULE variables
    !
    CALL deallocate_metric_on_particles
    CALL deallocate_SPH_memory

    PRINT *, " * LORENE SPH ID on the particles saved to formatted " &
             // "file", TRIM(namefile)

    PRINT *, "** Subroutine read_sphincs_dump_print_formatted " &
             // "executed."
    PRINT *

  END PROCEDURE read_sphincs_dump_print_formatted


  MODULE PROCEDURE print_formatted_lorene_id_particles

    !************************************************
    !                                               *
    ! Print the LORENE ID on the particles in a     *
    ! formatted file                                *
    !                                               *
    ! FT 18.09.2020                                 *
    !                                               *
    !************************************************

    USE constants, ONLY: c_light2, cm2m

    IMPLICIT NONE

    INTEGER:: itr, min_y_index

    DOUBLE PRECISION:: min_abs_y, min_abs_z1, min_abs_z2
    DOUBLE PRECISION, DIMENSION( :, : ), ALLOCATABLE:: abs_pos

    LOGICAL:: exist

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    ! Being abs_pos a local array, it is good practice to allocate it on the
    ! heap, otherwise it will be stored on the stack which has a very limited
    ! size. This results in a segmentation fault.
    ALLOCATE( abs_pos( 3, THIS% npart ) )

    PRINT *, "** Executing the print_formatted_lorene_id_particles " &
             // "subroutine..."
    PRINT *

    IF( THIS% call_flag == 0 )THEN
      PRINT *, "** The SUBROUTINE print_formatted_lorene_id_particles must", &
               " be called after compute_and_export_SPH_variables, otherwise", &
               " there are no SPH fields to export to the formatted file."
      PRINT *, "   Aborting."
      PRINT *
      STOP
    ENDIF

    IF( PRESENT(namefile) )THEN
      finalnamefile= namefile
    ELSE
      finalnamefile= "lorene-bns-id-particles-form.dat"
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
    !CALL test_status( ios, err_msg, "...error when opening " &
    !                  // TRIM(finalnamefile) )

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Values of the fields (including coordinates) exported by LORENE "&
    // "on each grid point"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 1 in " // TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 1 in "&
    !        // TRIM(finalnamefile) )

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# column:      1        2       3       4       5", &
    "       6       7       8", &
    "       9       10      11", &
    "       12      13      14", &
    "       15      16      17      18     19     20      21", &
    "       22      23      24      25     26     27      28      29"

    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 2 in " // TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 2 in "&
    !            // TRIM(finalnamefile) )

    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
  "#      grid point      x [km]       y [km]       z [km]       lapse", &
  "       shift_x [c]    shift_y [c]    shift_z [c]", &
  "       baryon density in the local rest frame [kg m^{-3}$]", &
  "       energy density [c^2]", &
  "       specific energy [c^2]", &
  "       pressure [Pa]", &
  "       fluid 3-velocity wrt the Eulerian observer (3 columns) [c]", &
  "       fluid coordinate 3-velocity vel_u (3 columns) [c]", &
  "       baryon number per particle nu", &
  "       baryon density in the local rest frame nlrf [baryon/Msun_geo^3]", &
  "       electron fraction", &
  "       generalized Lorentz factor Theta", &
  "       computing frame baryon number density", &
  "       computing frame baryon number density from kernel interpolation", &
  "       smoothing length", &
  "       particle density [particle/Msun_geo^3]", &
  "       particle volume [1/Msun_geo^3]", &
  "       particle mass [Msun]"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 3 in " // TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 3 in "&
    !          // TRIM(finalnamefile) )

    !DO itr = 1, THIS% npart, 1
    !  abs_pos( 1, itr )= ABS( THIS% pos( 1, itr ) )
    !  abs_pos( 2, itr )= ABS( THIS% pos( 2, itr ) )
    !  abs_pos( 3, itr )= ABS( THIS% pos( 3, itr ) )
    !ENDDO
    !
    !min_y_index= 0
    !min_abs_y= 1D+20
    !DO itr = 1, THIS% npart, 1
    !  IF( ABS( THIS% pos( 2, itr ) ) < min_abs_y )THEN
    !    min_abs_y= ABS( THIS% pos( 2, itr ) )
    !    min_y_index= itr
    !  ENDIF
    !ENDDO
    !
    !min_abs_z1= MINVAL( abs_pos( 3, 1:THIS% npart1 ) )
    !min_abs_z2= MINVAL( abs_pos( 3, THIS% npart1+1:THIS% npart ) )

    write_data_loop: DO itr = 1, THIS% npart, 1

      IF( THIS% export_form_xy .AND. &
          !( THIS% pos( 3, itr ) /= min_abs_z1 .AND. &
          !  THIS% pos( 3, itr ) /= min_abs_z2 ) &
          ( THIS% pos( 3, itr ) >=  0.5D0 .OR. &
            THIS% pos( 3, itr ) <= -0.5D0 ) &
      )THEN
        CYCLE
      ENDIF
      IF( THIS% export_form_x .AND. &
          !( THIS% pos( 3, itr ) /= min_abs_z1 &
          !.OR. THIS% pos( 3, itr ) /= min_abs_z2 &
          !.OR. THIS% pos( 2, itr ) /= THIS% pos( 2, min_y_index ) ) )THEN
          ( THIS% pos( 3, itr ) >=  0.5D0 .OR. &
            THIS% pos( 3, itr ) <= -0.5D0 .OR. &
            THIS% pos( 2, itr ) >=  0.5D0 .OR. &
            THIS% pos( 2, itr ) <= -0.5D0 ) &
      )THEN
        CYCLE
      ENDIF
      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        itr, &
        THIS% pos( 1, itr ), &
        THIS% pos( 2, itr ), &
        THIS% pos( 3, itr ), &
        THIS% lapse_parts( itr ), &
        THIS% shift_parts_x( itr ), &
        THIS% shift_parts_y( itr ), &
        THIS% shift_parts_z( itr ), &
        THIS% baryon_density_parts( itr ), &
        THIS% energy_density_parts( itr ), &
        THIS% specific_energy_parts( itr ), &
        THIS% pressure_parts( itr ), &
        THIS% v_euler_parts_x( itr ), &
        THIS% v_euler_parts_y( itr ), &
        THIS% v_euler_parts_z( itr ), &
        THIS% v( 1, itr ), &
        THIS% v( 2, itr ), &
        THIS% v( 3, itr ), &
        THIS% nu( itr ), &
        THIS% nlrf( itr ), &
        THIS% Ye( itr ), &
        THIS% Theta( itr ), &
        THIS% nstar( itr ), &
        THIS% nstar_int( itr ), &
        THIS% h( itr ), &
        THIS% particle_density( itr ), &
        THIS% particle_density_int( itr ), &
        THIS% pvol( itr ), &
        THIS% pmass( itr )

    IF( ios > 0 )THEN
      PRINT *, "...error when writing the arrays in " // TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing " &
    !         // "the arrays in " // TRIM(finalnamefile) )
    ENDDO write_data_loop

    CLOSE( UNIT= 2 )

    PRINT *, " * LORENE ID on particles saved to formatted file ", &
             TRIM(finalnamefile)
    PRINT *

    PRINT *, "** Subroutine print_formatted_lorene_id_particles executed."
    PRINT *

  END PROCEDURE print_formatted_lorene_id_particles


  MODULE PROCEDURE analyze_hydro

    !************************************************
    !                                               *
    ! Export the points where some of the hydro     *
    ! fields are negative to a formatted file       *
    !                                               *
    ! FT 5.12.2020                                  *
    !                                               *
    !************************************************

    IMPLICIT NONE

    LOGICAL:: exist, negative_hydro

    INQUIRE( FILE= TRIM(namefile), EXIST= exist )

    IF( exist )THEN
      OPEN( UNIT= 20, FILE= TRIM(namefile), STATUS= "REPLACE", &
            FORM= "FORMATTED", &
            POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
            IOMSG= err_msg )
    ELSE
      OPEN( UNIT= 20, FILE= TRIM(namefile), STATUS= "NEW", &
      FORM= "FORMATTED", &
            ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
       PRINT *, "...error when opening ",  TRIM(namefile), &
                " The error message is", err_msg
       STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when opening " &
    !         // TRIM(namefile) )

    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Points where some of the hydro fields are negative. "
    IF( ios > 0 )THEN
       PRINT *, "...error when writing line 1 in ",  TRIM(namefile), &
                " The error message is", err_msg
       STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 1 in "&
    !         // TRIM(namefile) )
    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# column:      1        2       3"
    IF( ios > 0 )THEN
       PRINT *, "...error when writing line 2 in ",  TRIM(namefile), &
                " The error message is", err_msg
       STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 2 in "&
    !        // TRIM(namefile) )
    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "#      x   y   z"
    IF( ios > 0 )THEN
       PRINT *, "...error when writing line 3 in ",  TRIM(namefile), &
                " The error message is", err_msg
       STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 3 in "&
    !        // TRIM(namefile) )

    DO itr= 1, THIS% npart, 1
      IF( THIS% baryon_density_parts ( itr ) < 0 .OR. &
          THIS% energy_density_parts ( itr ) < 0 .OR. &
          THIS% specific_energy_parts( itr ) < 0 .OR. &
          THIS% pressure_parts       ( itr ) < 0 )THEN

        negative_hydro= .TRUE.

        WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, &
               FMT = * )&
            THIS% pos( 1, itr ), &
            THIS% pos( 2, itr ), &
            THIS% pos( 3, itr )

        IF( ios > 0 )THEN
          PRINT *, "...error when writing the arrays in ", TRIM(namefile), &
                   " The error message is", err_msg
          STOP
        ENDIF
        !CALL test_status( ios, err_msg, "...error in writing "&
        !                // "the arrays in " // TRIM(namefile) )
      ENDIF
    ENDDO

    CLOSE( UNIT= 20 )

    IF( negative_hydro )THEN
      PRINT *, "** WARNING! Some of the hydro fields are negative on", &
               " some of the particles! See the file ", namefile, &
               " for the positions of such particles."
      PRINT *
    ELSE
      PRINT *, " * The hydro fields are positive on the particles."
      PRINT *
    ENDIF

  END PROCEDURE analyze_hydro


  !MODULE PROCEDURE write_lorene_bns_id_dump
  !
  !    !*************************************************
  !    !                                                *
  !    ! Returns the array of initial guess for the     *
  !    ! smoothing length                               *
  !    !                                                *
  !    ! FT                                             *
  !    !                                                *
  !    !*************************************************
  !
  !    USE input_output
  !    USE options, ONLY: basename
  !
  !    INTEGER:: a
  !
  !    LOGICAL:: exist
  !
  !    ! TODO: THIS OPTIONAL ARGUMENT DOES NOT WORK...
  !    IF( .NOT.PRESENT(TRIM(namefile)) )THEN
  !            TRIM(namefile)= "lorene-bns-id-particles-form.dat"
  !    ENDIF
  !
  !    INQUIRE( FILE= TRIM(namefile), EXIST= exist )
  !
  !    !PRINT *, TRIM(namefile)
  !    !PRINT *
  !
  !    IF( exist )THEN
  !        OPEN( UNIT= 3, FILE= TRIM(namefile), STATUS= "REPLACE", &
  !              FORM= "UNFORMATTED", &
  !              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
  !              IOMSG= err_msg )
  !    ELSE
  !        OPEN( UNIT= 3, FILE= TRIM(namefile), STATUS= "NEW", &
  !              FORM= "UNFORMATTED", &
  !              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
  !    ENDIF
  !    IF( ios > 0 )THEN
  !      PRINT *, "..error when opening " // TRIM(namefile)
  !               ". The error message is", err_msg
  !      STOP
  !    ENDIF
  !
  !    ! update dump counter
  !    dcount= dcount + 1
  !    ! construct file name
  !    basename= 'lbns.'
  !    CALL construct_filename(dcount,filename)
  !
  !    ! nlrf & nu are LARGE numbers --> scale for SPLASH
  !    nlrf= nlrf*m0c2_CU
  !    nu=   nu*amu/umass
  !
  !    ! write in MAGMA-type format
  !    WRITE( UNIT= 3, IOSTAT = ios, IOMSG = err_msg ) &
  !            npart,       & ! number of particles
  !            rstar,mstar, & ! radius and mass of the star
  !                           ! obsolete (see module_sph_variables)
  !            n1,n2,       & ! obsolete (see module_sph_variables)
  !            npm,         & ! obsolete (see module_sph_variables)
  !            t,           & ! time
  !            ( h(a), a=1, npart ),      & ! smoothing length
  !            escap,tkin,tgrav,tterm, & ! obsolete (see module_sph_variables)
  !            ( pos_u(1,a), a=1, npart), & ! particle positions
  !            ( pos_u(2,a), a=1, npart ),&
  !            ( pos_u(3,a), a=1, npart), &
  !            ( vel_u(1,a), a=1, npart ),& ! spatial coordinate velocity
  !            ( vel_u(2,a), a=1, npart), & ! of particles
  !            ( vel_u(3,a), a=1, npart ),&
  !            ( u(a), a=1, npart),       &
  !            ( nu(a), a=1, npart ),     &
  !            ( nlrf(a), a=1, npart),    &
  !            ( temp(a), a=1, npart ),   &
  !            ( Ye(a), a=1, npart),      &
  !            ( av(a), a=1, npart ),     & ! = 1
  !            ( divv(a), a=1, npart ),   & ! = 0
  !            ( Theta(a), a=1, npart ),  &
  !            ( Pr(a), a=1, npart )
  !            !
  !            !-- leave here for potential later use
  !            !
  !            !(pmasspm(a),a=1,npm),&
  !            !(pmpos(1,a),a=1,npm),(pmpos(2,a),a=1,npm),&
  !            !(pmpos(3,a),a=1,npm),(pmvel(1,a),a=1,npm),&
  !            !(pmvel(2,a),a=1,npm),(pmvel(3,a),a=1,npm),&
  !            !(pmdvdt(1,a),a=1,npm),(pmdvdt(2,a),a=1,npm),&
  !            !(pmdvdt(3,a),a=1,npm)
  !          IF( ios > 0 )THEN
  !            PRINT *, "..error when writing in " // TRIM(namefile)
  !                     ". The error message is", err_msg
  !            STOP
  !          ENDIF
  !          !CALL test_status( ios, err_msg, "...error when writing in " &
  !         !            // TRIM(namefile) )
  !
  !    CLOSE( UNIT= 3 )
  !
  !END PROCEDURE write_lorene_bns_id_dump


END SUBMODULE particles_sph_variables

