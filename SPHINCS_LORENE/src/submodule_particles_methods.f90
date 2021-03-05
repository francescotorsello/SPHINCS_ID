! File:         submodule_particles_methods.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (particles_id) particles_methods

  !***************************************************
  !                                                  *
  ! Implementation of the methods of TYPE particles. *
  !                                                  *
  ! FT 16.10.2020                                    *
  !                                                  *
  !***************************************************


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
    USE input_output,        ONLY: dcount, write_SPHINCS_dump

    IMPLICIT NONE

    ! The flag call_flag is set different than 0 if the SUBROUTINE
    ! compute_and_export_SPH_variables is called
    INTEGER, SAVE:: call_flag= 0

    ! Spacetime indices \mu and \nu
    INTEGER:: nus, mus

    DOUBLE PRECISION:: g4(0:3,0:3)
    DOUBLE PRECISION:: det,sq_g, Theta_a, temp1, temp2

    LOGICAL, PARAMETER:: debug= .FALSE.

    CHARACTER( LEN= : ), ALLOCATABLE:: compose_namefile


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

    CALL allocate_SPH_memory
    CALL allocate_metric_on_particles( THIS% npart )

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
    IF(.NOT.ALLOCATED( THIS% Ye ))THEN
      ALLOCATE( THIS% Ye( THIS% npart ), STAT= ios, &
            ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array Ye ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array lapse_parts" )
    ENDIF

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
      THIS% h(itr)= 3.0*(THIS% vol_a/(Msun_geo**3))**third
      h(itr)= THIS% h(itr)

      ! Baryon density in the local rest frame [baryon (Msun_geo)^{-3}]
      ! Computed from the LORENE baryon mass density in [kg/m^3]
      nlrf(itr)= THIS% baryon_density_parts(itr)*((Msun_geo*km2m)**3)/(amu*g2kg)
      THIS% nlrf(itr)= nlrf(itr)

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
      nu(itr)= nlrf(itr)*THIS% vol_a*Theta_a*sq_g
      THIS% nu(itr)= nu(itr)
      THIS% nbar_tot= THIS% nbar_tot + THIS% nu(itr)

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

    IF( THIS% compose_eos )THEN
      compose_namefile= TRIM(THIS% compose_path)//TRIM(THIS% compose_filename)
      CALL THIS% read_compose_composition( compose_namefile )
      CALL THIS% compute_Ye()
    ENDIF

    assign_ye_on_particles: DO itr= 1, THIS% npart, 1

      ! Electron fraction
      IF( THIS% compose_eos )THEN
        ye(itr)= THIS% Ye(itr)
      ELSE
        ye(itr)= 0.0D0
        THIS% Ye(itr)= 0.0D0
      ENDIF

    ENDDO assign_ye_on_particles
    CALL THIS% sph_computer_timer% stop_timer()

    PRINT "(A28,E15.8,A10)", "  * Maximum baryon density= ", &
             MAXVAL( THIS% baryon_density_parts, DIM= 1 )*kg2g/(m2cm**3), &
             " g cm^{-3}"
    PRINT "(A28,E15.8,A10)", "  * Minimum baryon density= ", &
             MINVAL( THIS% baryon_density_parts, DIM= 1 )*kg2g/(m2cm**3), &
             " g cm^{-3}"
    PRINT *, " * Ratio between the two=", &
             MAXVAL( THIS% baryon_density_parts, DIM= 1 )/ &
             MINVAL( THIS% baryon_density_parts, DIM= 1 )
    PRINT *

    PRINT *, " * Maximum nlrf=", MAXVAL( THIS% nlrf, DIM= 1 ), &
             "baryon Msun_geo^{-3}"
    PRINT *, " * Minimum nlrf=", MINVAL( THIS% nlrf, DIM= 1 ), &
             "baryon Msun_geo^{-3}"
    PRINT *, " * Ratio between the two=", &
                      MAXVAL( THIS% nlrf, DIM= 1 )/MINVAL( THIS% nlrf, DIM= 1 )
    PRINT *

    PRINT *, " * Maximum n. baryon per particle (nu)=", &
              MAXVAL( THIS% nu, DIM= 1 )
    PRINT *, " * Minimum n. baryon per particle (nu)=", &
              MINVAL( THIS% nu, DIM= 1 )
    PRINT *, " * Ratio between the two=", &
                      MAXVAL( THIS% nu, DIM= 1 )/MINVAL( THIS% nu, DIM= 1 )
    PRINT *, " * Total number of baryons=", THIS% nbar_tot
    PRINT *
    PRINT *, " * Total mass of the considered baryons=", &
             THIS% nbar_tot*amu/Msun, "Msun =", &
             THIS% nbar_tot*amu/Msun/(THIS% mass1 + THIS% mass2), &
             "of the total baryon mass."
    PRINT *

    !
    !-- Adding the missing baryons into the particles uniformly so that
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
      THIS% nu= THIS% nu/(THIS% nbar_tot*amu/Msun/(THIS% mass1 + THIS% mass2))
      nu= nu/(THIS% nbar_tot*amu/Msun/(THIS% mass1 + THIS% mass2))
      THIS% nbar_tot= &
          THIS% nbar_tot/(THIS% nbar_tot*amu/Msun/(THIS% mass1 + THIS% mass2))
      PRINT *, " * Total mass of the considered corrected baryons=", &
               THIS% nbar_tot*amu/Msun, "Msun =", &
               THIS% nbar_tot*amu/Msun/(THIS% mass1 + THIS% mass2), &
               "of the total baryon mass."
      PRINT *
    ENDIF

    !
    !-- Exporting the SPH ID to a binary file, for evolution
    !
    IF( THIS% export_bin )THEN
      IF( PRESENT(namefile) )THEN
        basename= namefile
      ELSE
        basename= "NSNS."
      ENDIF
      dcount= -1 ! since it is increased before writing
      CALL write_SPHINCS_dump()
    ENDIF

    CALL deallocate_metric_on_particles
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

    DOUBLE PRECISION:: min_abs_y, min_abs_z
    DOUBLE PRECISION, DIMENSION( :, : ), ALLOCATABLE:: abs_pos

    LOGICAL:: exist

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    PRINT *, "** Executing the read_bssn_dump_print_formatted subroutine..."

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

    min_abs_z= MINVAL( abs_pos( 3, : ) )

    write_data_loop: DO itr = 1, THIS% npart, 1

      IF( THIS% export_form_xy .AND. THIS% pos( 3, itr ) /= min_abs_z )THEN
        CYCLE
      ENDIF
      IF( THIS% export_form_x .AND. ( THIS% pos( 3, itr ) /= min_abs_z &
          .OR. THIS% pos( 2, itr ) /= THIS% pos( 2, min_y_index ) ) )THEN
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


  MODULE PROCEDURE read_compose_composition

    !************************************************
    !                                               *
    ! Read the electron fraction Y_e = n_e/n_b,     *
    ! with n_e electron number density and n_b      *
    ! baryon number density, from the .compo file   *
    ! taken from the CompOSE database of EoS.       *
    ! Y_e is given as a function of T, n_b, Y_q on  *
    ! a grid; the comutation of Y_e on the stars is *
    ! done by the SUBROUTINE compute_Ye_on_stars.   *
    !                                               *
    ! FT 1.03.2021                                  *
    !                                               *
    !************************************************

    USE constants, ONLY: fm2cm, cm2km, km2Msun_geo

    IMPLICIT NONE

    INTEGER:: itr, min_y_index, i_T, i_nb, i_yq, i_phase, n_pairs, i_e, cntr, &
              i_n, Y_n, &
              n_quad, &
              i_i, A_i, Z_i, Y_i, i_leptons, i_ns
    INTEGER, PARAMETER:: unit_compose= 56
    INTEGER, PARAMETER:: max_length_eos= 10000

    DOUBLE PRECISION:: min_abs_y, min_abs_z, m_n, m_p, &
                       p_nb, s_nb, mub_mn, muq_mn, mul_mn, f_nbmn, e_nbmn, h
    DOUBLE PRECISION, DIMENSION( :, : ), ALLOCATABLE:: abs_pos

    !DOUBLE PRECISION, DIMENSION( : ), ALLOCATABLE:: n_b, Y_e

    LOGICAL:: exist

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    PRINT *, "** Executing the read_compose_composition subroutine..."

    ALLOCATE( THIS% nb_table( max_length_eos ) )
    ALLOCATE( THIS% Ye_table( max_length_eos ) )
    THIS% nb_table= 0.0D0
    THIS% Ye_table= 0.0D0


    IF( PRESENT(namefile) )THEN
      finalnamefile= TRIM(namefile)//".beta"
    ELSE
      finalnamefile= "../../CompOSE_EOS/SFHO_with_electrons/eos.beta"
    ENDIF

    INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

    IF( exist )THEN
      OPEN( UNIT= unit_compose, FILE= TRIM(finalnamefile), &
            FORM= "FORMATTED", ACTION= "READ", IOSTAT= ios, &
            IOMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...error when opening " // TRIM(finalnamefile), &
                ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when opening " &
      !                  // TRIM(finalnamefile) )
    ELSE
      PRINT *, "** ERROR! Unable to find file " // TRIM(finalnamefile)
      STOP
    ENDIF

    !READ( UNIT= unit_compose, FMT= *, IOSTAT = ios, IOMSG= err_msg ) &
    !                                                    m_n, m_p, i_leptons

    PRINT *, " * Reading file " // TRIM(finalnamefile) // "..."
    cntr= 0
    read_compose_beta: DO itr= 1, max_length_eos, 1
      READ( UNIT= unit_compose, FMT= *, IOSTAT = ios, IOMSG= err_msg ) &
                        ! Variables for .thermo.ns file
                        !i_T, i_nb, i_yq, &
                        !p_nb, s_nb, mub_mn, muq_mn, mul_mn, f_nbmn, e_nbmn, &
                        !i_ns, THIS% Ye_table(itr), h
                        ! Variables for .beta file
                        THIS% nb_table(itr), &
                        THIS% Ye_table(itr)
                        ! Variables for .compo file
                        !i_T, i_nb, i_yq, &
                        !i_phase, n_pairs, &
                        !i_e, Y_e(itr)!, &
                        !i_n, Y_n, &
                        !n_quad, &
                        !i_i, A_i, Z_i, Y_i
      IF( ios > 0 )THEN
        PRINT *, "...error when opening " // TRIM(finalnamefile), &
                ". The error message is", err_msg
        STOP
      ENDIF
      IF( ios < 0 )THEN
        PRINT *, " * Reached end of file " // TRIM(finalnamefile)
        EXIT
      ENDIF
      cntr= cntr + 1
    ENDDO read_compose_beta
    !PRINT *, "cntr= ", cntr
    ! Reallocate the arrays to delete the empty elements
    THIS% nb_table= THIS% nb_table( 1:cntr )/(fm2cm**3*cm2km**3*km2Msun_geo**3)
    THIS% Ye_table= THIS% Ye_table( 1:cntr )
    !PRINT *, i_T, i_nb, i_yq, i_phase, n_pairs, i_e
    !PRINT *, "SIZE(n_b)= ", SIZE(THIS% n_b), "SIZE(Y_e)= ", SIZE(THIS% Y_e)
    !PRINT *, "n_b(1)= ", THIS% n_b(1), "Y_e(1)= ", THIS% Y_e(1)
    !PRINT *, "n_b(cntr)= ", THIS% n_b(cntr), "Y_e(cntr)= ", THIS% Y_e(cntr)
    !STOP

    PRINT *, "** Subroutine read_compose_composition executed."
    PRINT *

  END PROCEDURE read_compose_composition


  MODULE PROCEDURE compute_Ye

    !************************************************
    !                                               *
    ! Interpolate the electron fraction             *
    ! Y_e = n_e/n_b                                 *
    ! at the particle positions, using the data     *
    ! read by read_compose_composition.             *
    !                                               *
    ! FT 3.03.2021                                  *
    !                                               *
    !************************************************

    IMPLICIT NONE

    INTEGER:: d, itr2

    DOUBLE PRECISION:: min_nb_table, max_nb_table

    d= SIZE(THIS% nb_table)
    min_nb_table= MINVAL( THIS% nb_table, 1 )
    max_nb_table= MAXVAL( THIS% nb_table, 1 )

    particle_loop: DO itr= 1, THIS% npart, 1

      ! There may be particles with initially negative hydro fields,
      ! close to the surface of the stars. The present
      ! version of the code sets their hydro fields to 0.
      ! In turn, this is a problem if we read Ye from the .beta file
      ! from the CompOSe database since the value of 0
      ! for the baryon number density (used in the interpolation
      ! of Ye) is not included in the .beta file.
      ! In other words, Ye cannot be interpolated on those particles,
      ! and the present version of the code sets Ye to 0 as well.
      IF( THIS% nlrf(itr) == 0.0D0 )THEN
        THIS% Ye(itr)= 0.0D0
        CYCLE
      ELSEIF( THIS% nlrf(itr) < min_nb_table )THEN
        PRINT *, "** ERROR! The value of nlrf(", itr, ")=", THIS% nlrf(itr), &
                 "is lower than the minimum value in the table =", min_nb_table
        PRINT *, " * Is nlrf computed when you call this SUBROUTINE? " // &
                 "If yes, please select a table with a wider range."
        STOP
      ELSEIF( THIS% nlrf(itr) > max_nb_table )THEN
        PRINT *, "** ERROR! The value of nlrf(", itr, ")=", THIS% nlrf(itr), &
                 "is larger than the maximum value in the table =", max_nb_table
        PRINT *, " * Is nlrf computed when you call this SUBROUTINE? " // &
                 "If yes, please select a table with a wider range."
        STOP
      ENDIF

      Ye_linear_interpolation_loop: DO itr2= 1, d - 1, 1

        IF( THIS% nb_table(itr2) < THIS% nlrf(itr) .AND. &
            THIS% nlrf(itr) < THIS% nb_table(itr2 + 1) )THEN

          THIS% Ye(itr)= THIS% Ye_table(itr2) &
                         + (THIS% Ye_table(itr2 + 1) - THIS% Ye_table(itr2))/ &
                         (THIS% nb_table(itr2 + 1) - THIS% nb_table(itr2)) &
                         *(THIS% nlrf(itr) - THIS% nb_table(itr2))
          EXIT

        ENDIF

      ENDDO Ye_linear_interpolation_loop

      !PRINT *, "Ye(", itr, ")=", THIS% Ye(itr)
      !PRINT *

    ENDDO particle_loop

  END PROCEDURE compute_Ye


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

    DOUBLE PRECISION:: min_abs_y, min_abs_z
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
    "       15      16      17      18"

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
  "       baryon density in the local rest frame nlrf [baryon/cm^3]", &
  "       electron fraction", &
  "       generalized Lorentz factor Theta"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 3 in " // TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 3 in "&
    !          // TRIM(finalnamefile) )

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

    min_abs_z= MINVAL( abs_pos( 3, : ) )

    write_data_loop: DO itr = 1, THIS% npart, 1

      IF( THIS% export_form_xy .AND. THIS% pos( 3, itr ) /= min_abs_z )THEN
        CYCLE
      ENDIF
      IF( THIS% export_form_x .AND. ( THIS% pos( 3, itr ) /= min_abs_z &
          .OR. THIS% pos( 2, itr ) /= THIS% pos( 2, min_y_index ) ) )THEN
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
        THIS% Theta( itr )

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


  MODULE PROCEDURE destruct_particles

    !*************************************************
    !                                                *
    ! Destructor of a particles object               *
    ! It deallocates memory                          *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    IF( ALLOCATED( THIS% pos ))THEN
      DEALLOCATE( THIS% pos, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array pos. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array pos in SUBROUTINE"&
      !                // "destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% baryon_density_parts ))THEN
      DEALLOCATE( THIS% baryon_density_parts, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array baryon_density_parts. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array " &
      !                // "baryon_density_parts in SUBROUTINE " &
      !                // "destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% energy_density_parts ))THEN
      DEALLOCATE( THIS% energy_density_parts, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array energy_density_parts. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array " &
      !                // "energy_density_parts in SUBROUTINE " &
      !                // "destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% specific_energy_parts ))THEN
      DEALLOCATE( THIS% specific_energy_parts, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array specific_energy_parts. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array " &
      !                // "specific_energy_parts in SUBROUTINE " &
      !                // "destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% pressure_parts ))THEN
      DEALLOCATE( THIS% pressure_parts, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array pressure_parts. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array " &
      !                // "pressure_parts in SUBROUTINE " &
      !                // "destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% pressure_parts_cu ))THEN
      DEALLOCATE( THIS% pressure_parts_cu, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array pressure_parts_cu. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array " &
      !                // "pressure_parts_cu in SUBROUTINE " &
      !                // "destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% v_euler_parts_x ))THEN
      DEALLOCATE( THIS% v_euler_parts_x, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array v_euler_parts_x. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array " &
      !                // "v_euler_parts_x in SUBROUTINE " &
      !                // "destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% v_euler_parts_y ))THEN
      DEALLOCATE( THIS% v_euler_parts_y, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array v_euler_parts_y. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array " &
      !                // "v_euler_parts_y in SUBROUTINE " &
      !                // "destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% v_euler_parts_z ))THEN
      DEALLOCATE( THIS% v_euler_parts_z, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array v_euler_parts_z. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array " &
      !                // "v_euler_parts_z in SUBROUTINE " &
      !                // "destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% lapse_parts ))THEN
      DEALLOCATE( THIS% lapse_parts, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array lapse_parts. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array lapse_parts in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% shift_parts_x ))THEN
      DEALLOCATE( THIS% shift_parts_x, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array shift_parts_x. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array shift_parts_x in "&
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% shift_parts_y ))THEN
      DEALLOCATE( THIS% shift_parts_y, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array shift_parts_y. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array shift_parts_y in "&
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% shift_parts_z ))THEN
      DEALLOCATE( THIS% shift_parts_z, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array shift_parts_z. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array shift_parts_z in "&
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% g_xx_parts ))THEN
      DEALLOCATE( THIS% g_xx_parts, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_xx_parts. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_xx_parts in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% g_xy_parts ))THEN
      DEALLOCATE( THIS% g_xy_parts, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_xy_parts. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_xy_parts in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% g_xz_parts ))THEN
      DEALLOCATE( THIS% g_xz_parts, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_xz_parts. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_xz_parts in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% g_yy_parts ))THEN
      DEALLOCATE( THIS% g_yy_parts, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_yy_parts. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_yy_parts in " &
      !                // "SUBROUTINE estruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% g_yz_parts ))THEN
      DEALLOCATE( THIS% g_yz_parts, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_yz_parts. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_yz_parts in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% g_zz_parts ))THEN
      DEALLOCATE( THIS% g_zz_parts, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_zz_parts. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_zz_parts in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% nlrf ))THEN
      DEALLOCATE( THIS% nlrf, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array nlrf. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array nlrf in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% nu ))THEN
      DEALLOCATE( THIS% nu, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array nu. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array nu in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% Theta ))THEN
      DEALLOCATE( THIS% Theta, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array Theta. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array Theta in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% h ))THEN
      DEALLOCATE( THIS% h, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array h. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array h in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% v ))THEN
      DEALLOCATE( THIS% v, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array v. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF

  END PROCEDURE destruct_particles


  MODULE PROCEDURE is_empty

    !*************************************************
    !                                                *
    ! Returns the variable empty_object              *
    ! (experimental)                                 *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    answer= THIS% empty_object

  END PROCEDURE is_empty


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


  !-----------------!
  !--  FUNCTIONS  --!
  !-----------------!


  MODULE PROCEDURE get_npart

    !*************************************************
    !                                                *
    ! Returns the total number of particles          *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    n_part= THIS% npart

  END PROCEDURE get_npart


  MODULE PROCEDURE get_pos

    !*************************************************
    !                                                *
    ! Returns the array of particle positions        *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    pos_u= THIS% pos

  END PROCEDURE get_pos


  MODULE PROCEDURE get_vel

    !*************************************************
    !                                                *
    ! Returns the array of coordinate 3-velocity of  *
    ! particles                                      *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    vel= THIS% v(1:3,:)

  END PROCEDURE get_vel


  MODULE PROCEDURE get_nlrf

    !*************************************************
    !                                                *
    ! Returns the array of baryon density in the     *
    ! local rest frame                               *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    nlrf= THIS% nlrf

  END PROCEDURE get_nlrf


  MODULE PROCEDURE get_nu

    !*************************************************
    !                                                *
    ! Returns the array of baryon per particle       *
    ! [baryon (Msun_geo)^{-3}]                       *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    nu= THIS% nu

  END PROCEDURE get_nu


  MODULE PROCEDURE get_u

    !*************************************************
    !                                                *
    ! Returns the array of specific internal         *
    ! energy [c^2]                                   *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    u= THIS% specific_energy_parts

  END PROCEDURE get_u


  MODULE PROCEDURE get_pressure

    !*************************************************
    !                                                *
    ! Returns the array of pressure [kg c^2 m^{-3}]  *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    pressure= THIS% pressure_parts

  END PROCEDURE get_pressure


  MODULE PROCEDURE get_pressure_cu

    !*************************************************
    !                                                *
    ! Returns the array of pressure in code units    *
    ! [amu*c**2/(Msun_geo**3)]                       *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    pressure_cu= THIS% pressure_parts_cu

  END PROCEDURE get_pressure_cu


  MODULE PROCEDURE get_theta

    !*************************************************
    !                                                *
    ! Returns the array of generalized Lorentz       *
    ! factor                                         *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    theta= THIS% Theta

  END PROCEDURE get_theta


  MODULE PROCEDURE get_h

    !*************************************************
    !                                                *
    ! Returns the array of initial guess for the     *
    ! smoothing length [Msun_geo]                    *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    h= THIS% h

  END PROCEDURE get_h


END SUBMODULE particles_methods
