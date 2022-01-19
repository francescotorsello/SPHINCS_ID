! File:         submodule_sph_particles_memory.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (sph_particles) memory

  !***************************************************
  !
  !# This SUBMODULE contains the implementation of
  !  the methods of TYPE sph_particles
  !  that place particles on 1 or 2 lattices around
  !  the stars.
  !
  !  FT 12.07.2021
  !
  !***************************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE allocate_lorene_id_parts_memory

    !************************************************
    !
    !# Allocate memory for the LORENE ID on the
    !  particles
    !
    !  FT 10.11.2020
    !
    !************************************************

    IMPLICIT NONE

    PRINT *, "** Executing allocate_lorene_id_parts_memory."

    IF(.NOT.ALLOCATED( THIS% pos ))THEN
      ALLOCATE( THIS% pos( 3, THIS% npart ), STAT= ios, &
            ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array pos ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array pos" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% lapse_parts ))THEN
      ALLOCATE( THIS% lapse_parts( THIS% npart ), STAT= ios, &
            ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array lapse_parts ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array lapse_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% shift_parts_x ))THEN
      ALLOCATE( THIS% shift_parts_x( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array shift_parts_x ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for shift_parts_x" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% shift_parts_y ))THEN
      ALLOCATE( THIS% shift_parts_y( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array shift_parts_y ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for shift_parts_y" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% shift_parts_z ))THEN
      ALLOCATE( THIS% shift_parts_z( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array shift_parts_z ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for shift_parts_z" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_xx_parts ))THEN
      ALLOCATE( THIS% g_xx_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_xx_parts ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array g_xx_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_xy_parts ))THEN
      ALLOCATE( THIS% g_xy_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_xy_parts ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for array g_xy_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_xz_parts ))THEN
      ALLOCATE( THIS% g_xz_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_xz_parts ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array g_xz_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_yy_parts ))THEN
      ALLOCATE( THIS% g_yy_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_yy_parts ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for array g_yy_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_yz_parts ))THEN
      ALLOCATE( THIS% g_yz_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_yz_parts ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array g_yz_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_zz_parts ))THEN
      ALLOCATE( THIS% g_zz_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array g_zz_parts ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array g_zz_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% baryon_density_parts ))THEN
      ALLOCATE( THIS% baryon_density_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array baryon_density_parts ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !     "...allocation error for array baryon_density_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% energy_density_parts ))THEN
      ALLOCATE( THIS% energy_density_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array energy_density_parts ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !             "...allocation error for array energy_density_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% specific_energy_parts ))THEN
      ALLOCATE( THIS% specific_energy_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array specific_energy_parts ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !    "...allocation error for array specific_energy_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% pressure_parts ))THEN
      ALLOCATE( THIS% pressure_parts( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array pressure_parts ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array pressure_parts" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% pressure_parts_cu ))THEN
      ALLOCATE( THIS% pressure_parts_cu( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array pressure_parts_cu ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array pressure_parts_cu" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% v_euler_parts_x ))THEN
      ALLOCATE( THIS% v_euler_parts_x( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array v_euler_parts_x ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array v_euler_parts_x" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% v_euler_parts_y ))THEN
      ALLOCATE( THIS% v_euler_parts_y( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array v_euler_parts_y ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for array v_euler_parts_y" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% v_euler_parts_z ))THEN
      ALLOCATE( THIS% v_euler_parts_z( THIS% npart ), STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...allocation error for array v_euler_parts_z ", &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...allocation error for array v_euler_parts_z" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% nstar ))THEN
        ALLOCATE( THIS% nstar( THIS% npart ), STAT= ios )
    ENDIF
    IF( ios > 0 )THEN
       PRINT *, '...allocation error for nstar'
       STOP
    ENDIF
    IF(.NOT.ALLOCATED( THIS% nstar_int ))THEN
        ALLOCATE( THIS% nstar_int( THIS% npart ), STAT= ios )
    ENDIF
    IF( ios > 0 )THEN
       PRINT *, '...allocation error for nstar_int'
       STOP
    ENDIF
    IF(.NOT.ALLOCATED( THIS% particle_density ))THEN
        ALLOCATE( THIS% particle_density( THIS% npart ), STAT= ios )
    ENDIF
    IF( ios > 0 )THEN
       PRINT *, '...allocation error for particle_density'
       STOP
    ENDIF
    IF(.NOT.ALLOCATED( THIS% particle_density_int ))THEN
        ALLOCATE( THIS% particle_density_int( THIS% npart ), STAT= ios )
    ENDIF
    IF( ios > 0 )THEN
       PRINT *, '...allocation error for particle_density_int'
       STOP
    ENDIF
    IF(.NOT.ALLOCATED( THIS% pmass ))THEN
      ALLOCATE( THIS% pmass( THIS% npart ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array pmass in SUBROUTINE" &
                  // " allocate_lorene_id_parts_memory. ", &
                  "The STAT variable is", ios, ". ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( THIS% u_pwp ))THEN
      ALLOCATE( THIS% u_pwp( THIS% npart ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array u_pwp in SUBROUTINE" &
                  // "allocate_lorene_id_parts_memory. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF
    IF(.NOT.ALLOCATED( THIS% nlrf_int ))THEN
      ALLOCATE( THIS% nlrf_int( THIS% npart ), &
                STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array nlrf_int in SUBROUTINE" &
                  // "allocate_lorene_id_parts_memory. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    ENDIF

    PRINT *, "** Subroutine allocate_lorene_id_parts_memory executed."
    PRINT *

  END PROCEDURE allocate_lorene_id_parts_memory


  MODULE PROCEDURE deallocate_lorene_id_parts_memory

    !*************************************************
    !
    !# Deallocate memory for the LORENE ID on the
    !  particles
    !
    !  FT 12.07.2021 (this was part of the destructor
    !                 of TYPE [[particles]]
    !                 before this date)
    !
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
    IF( ALLOCATED( THIS% Ye ))THEN
      DEALLOCATE( THIS% Ye, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array Ye. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% nstar ))THEN
      DEALLOCATE( THIS% nstar, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array nstar. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% nstar_int ))THEN
      DEALLOCATE( THIS% nstar_int, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array nstar_int. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% particle_density ))THEN
      DEALLOCATE( THIS% particle_density, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array particle_density. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% particle_density_int ))THEN
      DEALLOCATE( THIS% particle_density_int, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array particle_density_int. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% pmass ))THEN
      DEALLOCATE( THIS% pmass, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array pmass. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% u_pwp ))THEN
      DEALLOCATE( THIS% u_pwp, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array u_pwp. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF
    IF( ALLOCATED( THIS% nlrf_int ))THEN
      DEALLOCATE( THIS% nlrf_int, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array u_pwp. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v in " &
      !                // "SUBROUTINE destruct_particles." )
    ENDIF

  END PROCEDURE deallocate_lorene_id_parts_memory


END SUBMODULE memory
