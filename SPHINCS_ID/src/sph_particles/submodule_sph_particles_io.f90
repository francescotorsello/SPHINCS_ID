! File:         submodule_sph_particles_io.f90
! Authors:      Francesco Torsello (FT)
!************************************************************************
! Copyright (C) 2020, 2021, 2022 Francesco Torsello                     *
!                                                                       *
! This file is part of SPHINCS_ID                                       *
!                                                                       *
! SPHINCS_ID is free software: you can redistribute it and/or modify    *
! it under the terms of the GNU General Public License as published by  *
! the Free Software Foundation, either version 3 of the License, or     *
! (at your option) any later version.                                   *
!                                                                       *
! SPHINCS_ID is distributed in the hope that it will be useful,         *
! but WITHOUT ANY WARRANTY; without even the implied warranty of        *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          *
! GNU General Public License for more details.                          *
!                                                                       *
! You should have received a copy of the GNU General Public License     *
! along with SPHINCS_ID. If not, see <https://www.gnu.org/licenses/>.   *
! The copy of the GNU General Public License should be in the file      *
! 'COPYING'.                                                            *
!************************************************************************

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


  USE constants,  ONLY: amu, c_light2
  USE utility,    ONLY: km2m, m2cm, Msun_geo, one


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

    DOUBLE PRECISION:: max_nlrf_id, max_nlrf_sph, min_nlrf_sph
    DOUBLE PRECISION:: max_pr_id,   max_pr_sph,   min_pr_sph
    DOUBLE PRECISION:: max_u_id,    max_u_sph,    min_u_sph

    PRINT *, " * SPH:"
    PRINT *
    PRINT *, "   Total particle number= ", this% npart
    DO i_matter= 1, this% n_matter, 1
      PRINT *, "   Particle number on matter object    ", i_matter, "=", &
                                            this% npart_i(i_matter)
    ENDDO
    PRINT *
    DO i_matter= 1, this% n_matter, 1
      PRINT *, "   Mass fraction of matter object      ", i_matter, "=", &
               this% mass_fractions(i_matter)
      PRINT *, "   Particle fraction of matter object  ", i_matter, "=", &
               DBLE(this% npart_i(i_matter))/DBLE(this% npart)
      PRINT *, "   Baryon number ratio on matter object", i_matter, "=", &
               this% nuratio_i(i_matter)
    ENDDO
    PRINT *

    PRINT *, "   Baryon number ratio across all matter objects    =", &
             this% nuratio
    PRINT *

    PRINT *, "   Center of mass of the entire particle distribution    ="
    PRINT *, "   (", this% barycenter_system(1), ","
    PRINT *, "    ", this% barycenter_system(2), ","
    PRINT *, "    ", this% barycenter_system(3), ") Msun"
    PRINT *, "   Its distance from the origin is: ", &
             this% barycenter_system(4), " Msun"
    PRINT *

    DO i_matter= 1, this% n_matter, 1

      ASSOCIATE( npart_in   => this% npart_i(i_matter-1) + 1, &
                 npart_fin  => this% npart_i(i_matter-1) +    &
                               this% npart_i(i_matter) )

        max_nlrf_id = MAXVAL( this% nlrf(npart_in:npart_fin), DIM= 1 )
        max_nlrf_sph= MAXVAL( this% nlrf_int(npart_in:npart_fin), DIM= 1 )
        min_nlrf_sph= MINVAL( this% nlrf_int(npart_in:npart_fin), DIM= 1 )

        max_pr_id = MAXVAL( this% pressure(npart_in:npart_fin), DIM= 1 )
        max_pr_sph= MAXVAL( this% pressure_cu(npart_in:npart_fin), DIM= 1 )
        min_pr_sph= MINVAL( this% pressure_cu(npart_in:npart_fin), DIM= 1 )

        max_u_id = MAXVAL( this% specific_energy(npart_in:npart_fin), DIM= 1 )
        max_u_sph= MAXVAL( this% u_pwp(npart_in:npart_fin), DIM= 1 )
        min_u_sph= MINVAL( this% u_pwp(npart_in:npart_fin), DIM= 1 )

        PRINT *, "   * On matter object ", i_matter, ":"
        PRINT *

        PRINT *, "   Maximum nlrf from the ID = ", max_nlrf_id &
                 /((Msun_geo*km2m*m2cm)**3)*amu, " g cm^{-3}"
        PRINT *, "   Maximum SPH interpolated nlrf =", max_nlrf_sph &
                 /((Msun_geo*km2m*m2cm)**3)*amu, " g cm^{-3}"
        PRINT *, "   Minimum SPH interpolated nlrf =", min_nlrf_sph &
                 /((Msun_geo*km2m*m2cm)**3)*amu, " g cm^{-3}"
        PRINT *, "   Relative difference between the two maximum nlrf=", &
                 ABS( max_nlrf_sph - max_nlrf_id )/max_nlrf_id
        PRINT *, "   Ratio between maximum and minimum nlrf=", &
                 max_nlrf_sph/min_nlrf_sph
        PRINT *

        PRINT *, "   Maximum pressure from the ID = ", max_pr_id &
                 *amu*c_light2/((Msun_geo*km2m*m2cm)**3), " Ba"
        PRINT *, "   Maximum pressure from SPH interpolated density and EOS =",&
                 max_pr_sph*amu*c_light2/((Msun_geo*km2m*m2cm)**3), " Ba"
        PRINT *, "   Minimum pressure from SPH interpolated density and EOS =",&
                 min_pr_sph*amu*c_light2/((Msun_geo*km2m*m2cm)**3), " Ba"
        PRINT *, "   Relative difference between the maximum pressures=", &
                 ABS( max_pr_sph - max_pr_id )/max_pr_id
        PRINT *, "   Ratio between maximum and minimum pressure=", &
                 max_pr_sph/min_pr_sph
        PRINT *

        PRINT *, "   Maximum specific internal energy from the ID =", &
                 max_u_id, " c^2"
        PRINT *, "   Maximum specific internal energy from SPH interpolated ", &
                 "density and EOS =", max_u_sph, " c^2"
        PRINT *, "   Minimum specific internal energy from SPH interpolated ", &
                 "density and EOS =", min_u_sph, " c^2"
        PRINT *, "   Relative difference between the maximum specific ", &
                 "internal energies=", ABS( max_u_sph - max_u_id )/max_u_id
        PRINT *, "   Ratio between maximum and minimum specific internal ", &
                 "energy=", max_u_sph/min_u_sph
        PRINT *

      END ASSOCIATE

    ENDDO

    DO i_matter= 1, this% n_matter, 1

      PRINT *, "   SPH estimate of the ADM linear momentum computed using ", &
               "the canonical momentum per baryon, on matter object", &
               i_matter,"= "
      PRINT *, "   (", this% adm_linear_momentum_i(i_matter, 1), ","
      PRINT *, "    ", this% adm_linear_momentum_i(i_matter, 2), ","
      PRINT *, "    ", this% adm_linear_momentum_i(i_matter, 3), ") Msun*c"
      PRINT *

    ENDDO
    PRINT *, "   SPH estimate of the ADM momentum of the fluid ", &
             "computed using the canonical momentum per baryon= "
    PRINT *, "   (", this% adm_linear_momentum_fluid(1), ","
    PRINT *, "    ", this% adm_linear_momentum_fluid(2), ","
    PRINT *, "    ", this% adm_linear_momentum_fluid(3), ") Msun*c"
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
                                   deallocate_SPH_memory, &
                                   n1, n2
    !USE metric_on_particles, ONLY: allocate_metric_on_particles, &
    !                               deallocate_metric_on_particles
    USE input_output,        ONLY: set_units, read_SPHINCS_dump
    USE alive_flag,          ONLY: alive

    IMPLICIT NONE

    INTEGER, PARAMETER:: max_npart= 2.D+7

    INTEGER:: itr

    LOGICAL:: exist, final_save_data

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    IF( PRESENT(save_data) )THEN
      final_save_data= save_data
    ELSE
      final_save_data= .FALSE.
    ENDIF

    PRINT *, "** Executing the read_sphincs_dump_print_formatted subroutine..."

    !
    !-- Set up the MODULE variables in MODULE sph_variables
    !-- (used by write_SPHINCS_dump)
    !
    !npart= 10*this% npart
    npart= max_npart

    CALL set_units('NSM')

    CALL allocate_SPH_memory
    !CALL allocate_metric_on_particles( npart )

    finalnamefile= TRIM(namefile_bin)!//"00000"
    CALL read_SPHINCS_dump( finalnamefile )

    PRINT *, " * Maximum interpolated nlrf =", &
             MAXVAL( nlrf(1:npart), DIM= 1 ) &
             /((Msun_geo*km2m*m2cm)**3)*amu, " g cm^{-3}"
    PRINT *, " * Minimum interpolated nlrf =", &
             MINVAL( nlrf(1:npart), DIM= 1 ) &
             /((Msun_geo*km2m*m2cm)**3)*amu, " g cm^{-3}"
    PRINT *, " * Ratio between the two=", &
             MAXVAL( nlrf(1:npart), DIM= 1 )/ &
             MINVAL( nlrf(1:npart), DIM= 1 )
    PRINT *

    PRINT *, " * Maximum pressure =", &
             MAXVAL( Pr(1:npart), DIM= 1 ) &
             *amu*c_light2/((Msun_geo*km2m*m2cm)**3), " Ba"
    PRINT *, " * Minimum pressure =", &
             MINVAL( Pr(1:npart), DIM= 1 ) &
             *amu*c_light2/((Msun_geo*km2m*m2cm)**3), " Ba"
    PRINT *, " * Ratio between the two=", &
             MAXVAL( Pr(1:npart), DIM= 1 )/ &
             MINVAL( Pr(1:npart), DIM= 1 )
    PRINT *

    PRINT *, " * Maximum specific internal energy =", &
             MAXVAL( u(1:npart), DIM= 1 ), " c^2"
    PRINT *, " * Minimum specific internal energy =", &
             MINVAL( u(1:npart), DIM= 1 ), " c^2"
    PRINT *, " * Ratio between the two=", &
             MAXVAL( u(1:npart), DIM= 1 )/ &
             MINVAL( u(1:npart), DIM= 1 )
    PRINT *

    PRINT *, " * Maximum nu =", &
             MAXVAL( nu(1:npart), DIM= 1 )
    PRINT *, " * Minimum nu =", &
             MINVAL( nu(1:npart), DIM= 1 )
    PRINT *, " * Ratio between the two=", &
             MAXVAL( nu(1:npart), DIM= 1 )/ &
             MINVAL( nu(1:npart), DIM= 1 )
    PRINT *

    !STOP
    IF( final_save_data )THEN

      ALLOCATE( this% npart_i(2) )
      ALLOCATE( alive( npart ) )
      alive( 1:npart )= 1

      this% npart      = npart
      CALL this% allocate_particles_memory

      this% npart_i(1) = n1
      this% npart_i(2) = n2
      this% pos        = pos_u(:,1:npart)
      this% v(0,:)     = one
      this% v(1:3,:)   = vel_u(:,1:npart)
      this% u_pwp      = u(1:npart)
      this% nu         = nu(1:npart)
      this% h          = h(1:npart)
      this% nlrf_int   = nlrf(1:npart)
      this% pressure_cu= Pr(1:npart)
      this% Ye         = Ye(1:npart)
      this% theta      = Theta(1:npart)

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
    "# Values of the fields (including coordinates) on the particles "
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
    "#      particle      x [Msun_geo]       y [Msun_geo]       z [Msun_geo]", &
    "       lapse     shift_x [c]    shift_y [c]    shift_z [c]", &
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

    write_data_loop: DO itr = 1, npart, 1

      IF( this% export_form_xy .AND. &
          !( pos_u( 3, itr ) /= min_abs_z1 .AND. &
          !  pos_u( 3, itr ) /= min_abs_z2 ) &
          ( pos_u( 3, itr ) >=  0.5D0 .OR. &
            pos_u( 3, itr ) <= -0.5D0 ) &
      )THEN
        CYCLE
      ENDIF
      IF( this% export_form_x .AND. &
          !( pos_u( 3, itr ) /= min_abs_z1 &
          !.OR. pos_u( 3, itr ) /= min_abs_z2 &
          !.OR. pos_u( 2, itr ) /= this% pos( 2, min_y_index ) ) )THEN
          ( pos_u( 3, itr ) >=  0.5D0 .OR. &
            pos_u( 3, itr ) <= -0.5D0 .OR. &
            pos_u( 2, itr ) >=  0.5D0 .OR. &
            pos_u( 2, itr ) <= -0.5D0 ) &
      )THEN
        CYCLE
      ENDIF
      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        itr, &                                  ! 1     33
        pos_u( 1, itr ), &                      ! 2     34
        pos_u( 2, itr ), &                      ! 3     35
        pos_u( 3, itr ), &                      ! 4     36
        vel_u( 1, itr ), &                      ! 5     37
        vel_u( 2, itr ), &                      ! 6     38
        vel_u( 3, itr ), &                      ! 7     39
        h( itr ), &                             ! 8     40
        u( itr ), &                             ! 9     41
        nu( itr ), &                            ! 10    42
        nlrf( itr ), &                          ! 11    43
        temp( itr ), &                          ! 12    44
        av( itr ), &                            ! 13    45
        ye( itr ), &                            ! 14    46
        divv( itr ), &                          ! 15    47
        Theta( itr ), &                         ! 16    48
        Pr( itr )                               ! 17    49
        !this% pos( 1, itr ), &                  ! 2     34
        !this% pos( 2, itr ), &                  ! 3     35
        !this% pos( 3, itr ), &                  ! 4     36
        !this% v( 1, itr ), &                    ! 5     37
        !this% v( 2, itr ), &                    ! 6     38
        !this% v( 3, itr ), &                    ! 7     39
        !this% h( itr ), &                       ! 8     40
        !this% u_pwp( itr ), &                   ! 9     41
        !this% nu( itr ), &                      ! 10    42
        !this% nlrf_int( itr ), &                ! 11    43
        !temp( itr ), &                          ! 12    44
        !av( itr ), &                            ! 13    45
        !ye( itr ), &                            ! 14    46
        !divv( itr ), &                          ! 15    47
        !this% theta( itr ), &                   ! 16    48
        !this% pressure_cu( itr )                ! 17    49

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
    !CALL deallocate_metric_on_particles
    CALL deallocate_SPH_memory

    PRINT *, " * LORENE SPH ID on the particles saved to formatted " &
             // "file ", TRIM(namefile)

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
    ALLOCATE( abs_pos( 3, this% npart ) )

    PRINT *, "** Executing the print_formatted_id_particles " &
             // "subroutine..."
    PRINT *

    IF( this% call_flag == 0 )THEN
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

    !DO itr = 1, this% npart, 1
    !  abs_pos( 1, itr )= ABS( this% pos( 1, itr ) )
    !  abs_pos( 2, itr )= ABS( this% pos( 2, itr ) )
    !  abs_pos( 3, itr )= ABS( this% pos( 3, itr ) )
    !ENDDO
    !
    !min_y_index= 0
    !min_abs_y= 1D+20
    !DO itr = 1, this% npart, 1
    !  IF( ABS( this% pos( 2, itr ) ) < min_abs_y )THEN
    !    min_abs_y= ABS( this% pos( 2, itr ) )
    !    min_y_index= itr
    !  ENDIF
    !ENDDO
    !
    !min_abs_z1= MINVAL( abs_pos( 3, 1:this% npart1 ) )
    !min_abs_z2= MINVAL( abs_pos( 3, this% npart1+1:this% npart ) )

    write_data_loop: DO itr = 1, this% npart, 1

      IF( this% export_form_xy .AND. &
          !( this% pos( 3, itr ) /= min_abs_z1 .AND. &
          !  this% pos( 3, itr ) /= min_abs_z2 ) &
          ( this% pos( 3, itr ) >=  0.5D0 .OR. &
            this% pos( 3, itr ) <= -0.5D0 ) &
      )THEN
        CYCLE
      ENDIF
      IF( this% export_form_x .AND. &
          !( this% pos( 3, itr ) /= min_abs_z1 &
          !.OR. this% pos( 3, itr ) /= min_abs_z2 &
          !.OR. this% pos( 2, itr ) /= this% pos( 2, min_y_index ) ) )THEN
          ( this% pos( 3, itr ) >=  0.5D0 .OR. &
            this% pos( 3, itr ) <= -0.5D0 .OR. &
            this% pos( 2, itr ) >=  0.5D0 .OR. &
            this% pos( 2, itr ) <= -0.5D0 ) &
      )THEN
        CYCLE
      ENDIF
      WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        itr, &                                                       ! 1
        this% pos( 1, itr ), &                                       ! 2
        this% pos( 2, itr ), &                                       ! 3
        this% pos( 3, itr ), &                                       ! 4
        this% lapse( itr ), &                                        ! 5
        this% shift_x( itr ), &                                      ! 6
        this% shift_y( itr ), &                                      ! 7
        this% shift_z( itr ), &                                      ! 8
        this% baryon_density( itr ), &                               ! 9
        this% energy_density( itr ), &                               ! 10
        this% specific_energy( itr ), &                              ! 11
        this% pressure( itr ), &                                     ! 12
        this% pressure_cu( itr ), &                                  ! 13
        this% v_euler_x( itr ), &                                    ! 14
        this% v_euler_y( itr ), &                                    ! 15
        this% v_euler_z( itr ), &                                    ! 16
        this% v( 1, itr ), &                                         ! 17
        this% v( 2, itr ), &                                         ! 18
        this% v( 3, itr ), &                                         ! 19
        this% nu( itr ), &                                           ! 20
        this% nlrf( itr ), &                                         ! 21
        this% Ye( itr ), &                                           ! 22
        this% Theta( itr ), &                                        ! 23
        this% nstar( itr ), &                                        ! 24
        this% nstar_int( itr ), &                                    ! 25
        this% h( itr ), &                                            ! 26
        this% particle_density( itr ), &                             ! 27
        this% particle_density_int( itr ), &                         ! 28
        this% pvol( itr ), &                                         ! 29
        this% pmass( itr ), &                                        ! 30
        this% u_pwp( itr ), &                                        ! 31
        this% nlrf_int( itr )                                        ! 32

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

    DO itr= 1, this% npart, 1

      IF( this% baryon_density ( itr ) < 0 .OR. &
          this% energy_density ( itr ) < 0 .OR. &
          this% specific_energy( itr ) < 0 .OR. &
          this% pressure       ( itr ) < 0 )THEN

        negative_hydro= .TRUE.

        WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, &
               FMT = * )&

            this% pos( 1, itr ), &
            this% pos( 2, itr ), &
            this% pos( 3, itr )

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
  !    ! TODO: this OPTIONAL ARGUMENT DOES NOT WORK...
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
