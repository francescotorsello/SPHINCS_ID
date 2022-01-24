! File:         submodule_sph_particles_io.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (sph_particles) io

  !***************************************************
  !
  !# This submodule contains the implementation of the
  !  methods of TYPE sph_particles that handle I/O
  !  (input/output)
  !
  !  FT 5.11.2021
  !
  !***************************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE print_summary

    !************************************************
    !
    !# Prints a summary of the properties of the
    !  |sph| particle distribution, optionally, to
    !  a formatted file whose name
    !  is given as the optional argument `filename`
    !
    !  FT 5.11.2021
    !
    !************************************************

    IMPLICIT NONE

    INTEGER:: i_matter

    PRINT *, " * SPH:"
    PRINT *
    PRINT *, "   Total particle number= ", THIS% npart
    DO i_matter= 1, THIS% n_matter, 1
      PRINT *, "   Particle number on matter object    ", i_matter, "=", &
                                            THIS% npart_i(i_matter)
    ENDDO
    PRINT *
    DO i_matter= 1, THIS% n_matter, 1
      PRINT *, "   Mass fraction of matter object      ", i_matter, "=", &
               THIS% mass_fractions(i_matter)
      PRINT *, "   Particle fraction of matter object  ", i_matter, "=", &
               DBLE(THIS% npart_i(i_matter))/DBLE(THIS% npart)
      PRINT *, "   Baryon number ratio on matter object", i_matter, "=", &
               THIS% nuratio_i(i_matter)
    ENDDO
    PRINT *

    PRINT *, "   Baryon number ratio across all matter objects    =", &
             THIS% nuratio
    PRINT *

  END PROCEDURE print_summary


  MODULE PROCEDURE read_sphincs_dump_print_formatted

    !************************************************
    !
    !# Read the SPH ID from the binary file output
    !  by write_SPHINCS_dump, and print it to a
    !  formatted file
    !
    !  FT 12.02.2021
    !
    !************************************************

    USE sph_variables,       ONLY: npart, &  ! particle number
                                   pos_u, &  ! particle positions
                                   vel_u, &  ! particle velocities in
                                             ! coordinate frame
                                   nlrf,  &  ! baryon number density in
                                             ! local rest frame
                                   !ehat,  &  ! canonical energy per baryon
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
                                   deallocate_metric_on_particles
    USE input_output,        ONLY: set_units, read_SPHINCS_dump

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
      PRINT *, "** The SUBROUTINE print_formatted_id_particles must", &
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

    min_abs_z1= MINVAL( abs_pos( 3, 1:THIS% npart_i(1) ) )
    min_abs_z2= MINVAL( abs_pos( 3, THIS% npart_i(1)+1:THIS% npart ) )

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


  MODULE PROCEDURE print_formatted_id_particles

    !************************************************
    !
    !# Print the LORENE ID on the particles in a
    !  formatted file
    !
    !  FT 18.09.2020
    !
    !************************************************

    IMPLICIT NONE

    INTEGER:: itr

    DOUBLE PRECISION, DIMENSION( :, : ), ALLOCATABLE:: abs_pos

    LOGICAL:: exist

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    ! Being abs_pos a local array, it is good practice to allocate it on the
    ! heap, otherwise it will be stored on the stack which has a very limited
    ! size. This results in a segmentation fault.
    ALLOCATE( abs_pos( 3, THIS% npart ) )

    PRINT *, "** Executing the print_formatted_id_particles " &
             // "subroutine..."
    PRINT *

    IF( THIS% call_flag == 0 )THEN
      PRINT *, "** The SUBROUTINE print_formatted_id_particles must", &
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
        itr, &                                                             ! 1
        THIS% pos( 1, itr ), &                                             ! 2
        THIS% pos( 2, itr ), &                                             ! 3
        THIS% pos( 3, itr ), &                                             ! 4
        THIS% lapse_parts( itr ), &                                        ! 5
        THIS% shift_parts_x( itr ), &                                      ! 6
        THIS% shift_parts_y( itr ), &                                      ! 7
        THIS% shift_parts_z( itr ), &                                      ! 8
        THIS% baryon_density_parts( itr ), &                               ! 9
        THIS% energy_density_parts( itr ), &                               ! 10
        THIS% specific_energy_parts( itr ), &                              ! 11
        THIS% pressure_parts( itr ), &                                     ! 12
        THIS% pressure_parts_cu( itr ), &                                  ! 13
        THIS% v_euler_parts_x( itr ), &                                    ! 14
        THIS% v_euler_parts_y( itr ), &                                    ! 15
        THIS% v_euler_parts_z( itr ), &                                    ! 16
        THIS% v( 1, itr ), &                                               ! 17
        THIS% v( 2, itr ), &                                               ! 18
        THIS% v( 3, itr ), &                                               ! 19
        THIS% nu( itr ), &                                                 ! 20
        THIS% nlrf( itr ), &                                               ! 21
        THIS% Ye( itr ), &                                                 ! 22
        THIS% Theta( itr ), &                                              ! 23
        THIS% nstar( itr ), &                                              ! 24
        THIS% nstar_int( itr ), &                                          ! 25
        THIS% h( itr ), &                                                  ! 26
        THIS% particle_density( itr ), &                                   ! 27
        THIS% particle_density_int( itr ), &                               ! 28
        THIS% pvol( itr ), &                                               ! 29
        THIS% pmass( itr ), &                                              ! 30
        THIS% u_pwp( itr ), &                                              ! 31
        THIS% nlrf_int( itr )                                              ! 32

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

    PRINT *, "** Subroutine print_formatted_id_particles executed."
    PRINT *

  END PROCEDURE print_formatted_id_particles


  MODULE PROCEDURE analyze_hydro

    !************************************************
    !
    !# Export the points where some of the hydro
    !  fields are negative to a formatted file
    !  (deprecated?)
    !
    !  FT 5.12.2020
    !
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
  !    !
  !    !# Returns the array of initial guess for the
  !    !  smoothing length
  !    !
  !    !  FT
  !    !
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


END SUBMODULE io
