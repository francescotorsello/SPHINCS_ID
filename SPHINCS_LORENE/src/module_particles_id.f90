! File:         module_particles_id.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

MODULE particles_id


  !************************************************************
  !                                                           *
  !   This module contains the definition of TYPE particles   *
  !                                                           *
  !************************************************************


  USE utility, ONLY: itr, ios, err_msg, test_status, &
                     perc, creturn, run_id, show_progress
  USE bns_id,  ONLY: bns
  USE timing,  ONLY: timer


  IMPLICIT NONE


  !INTEGER:: max_particles


  !**********************************************************
  !                                                         *
  !              Definition of TYPE particles               *
  !                                                         *
  ! This class places the SPH particles, imports            *
  ! the LORENE BNS ID on the particle positions, stores     *
  ! it, computes the relevant SPH fields and exports it to  *
  ! both a formatted, and a binary file for evolution       *
  !                                                         *
  !**********************************************************

  TYPE:: particles


    PRIVATE


    INTEGER:: npart, npart1, npart2, npart_temp, npart1_temp, npart2_temp
    INTEGER:: nx, ny, nz
    ! The flag call_flag is set different than 0 if the SUBROUTINE
    ! compute_and_export_SPH_variables is called
    INTEGER:: call_flag= 0

    INTEGER, DIMENSION(:), ALLOCATABLE:: filt_pos

    !
    !-- Hydro variables on the particles
    !
    ! 2-D array storing the position of the particles
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos
    ! 1-D array storing the position of the particles on the x axis for S 1
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pos_x1
    ! 1-D array storing the position of the particles on the x axis for NS2
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pos_x2
    ! 1-D array storing the baryon mass density in the fluid frame
    ! [kg m^{-3}]
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: baryon_density_parts
    ! 1-D array storing the energy density [kg c^2 m^{-3}]
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: energy_density_parts
    ! 1-D array storing the specific internal energy [c^2]
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: specific_energy_parts
    ! 1-D array storing the pressure [kg c^2 m^{-3}]
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_parts
    ! 1-D array storing the pressure on the x axis [kg c^2 m^{-3}] for NS 1
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_parts_x1
    ! 1-D array storing the pressure on the x axis [kg c^2 m^{-3}] for NS 2
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_parts_x2
    ! 1-D array storing the first derivative of the pressure
    ! along the x axis [kg c^2 m^{-3}] for NS 1
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_parts_x_der1
    ! 1-D array storing the first derivative of the pressure
    ! along the x axis [kg c^2 m^{-3}] for NS2
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_parts_x_der2
    ! 1-D array storing the typical length scale for the pressure change
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_length_scale_x1
    ! 1-D array storing the typical length scale for the pressure change
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_length_scale_x2
    ! 1-D array storing the pressure in code units [amu c^2/(Msun_geo**3)]
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_parts_cu
    ! 1-D arrays storing the components of the fluid 3-velocity wrt
    ! the Eulerian observer [c]
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: v_euler_parts_x
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: v_euler_parts_y
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: v_euler_parts_z

    !
    !-- Spacetime fields
    !
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: lapse_parts
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: shift_parts_x
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: shift_parts_y
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: shift_parts_z
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: g_xx_parts
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: g_xy_parts
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: g_yy_parts
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: g_xz_parts
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: g_yz_parts
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: g_zz_parts

    !
    !-- SPH fields
    !
    ! 1-D array storing baryon density in the local rest frame
    ! baryon (Msun_geo)^{-3}
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nlrf
    ! 1-D array storing the baryon number per particle
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nu
    ! 2-D array storing the coordinate fluid 4-velocity [c]
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: v
    ! 1-D array storing the generalized Lorentz factor
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: Theta
    ! 1-D array storing the smoothing length
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: h
    ! Baryonic masses of the neutron stars [Msun]
    DOUBLE PRECISION:: mass1, mass2
    ! Total grid volume
    DOUBLE PRECISION:: vol, vol1, vol2
    ! Volume per particle
    DOUBLE PRECISION:: vol_a, vol1_a, vol2_a
    ! Total baryon number
    DOUBLE PRECISION:: nbar_tot

    CHARACTER( LEN= 50 ):: lorene_bns_id_parfile

    LOGICAL:: empty_object
    LOGICAL, PUBLIC:: export_bin
    LOGICAL, PUBLIC:: export_form_xy, export_form_x

    TYPE(timer), PUBLIC:: placer_timer
    TYPE(timer), PUBLIC:: importer_timer
    TYPE(timer), PUBLIC:: sph_computer_timer


    CONTAINS


    !-------------------!
    !--  SUBROUTINES  --!
    !-------------------!

    PROCEDURE:: place_particles_3dlattice

    PROCEDURE:: place_particles_3dlattices

    PROCEDURE:: allocate_lorene_id_parts_memory

    PROCEDURE, PUBLIC:: analyze_hydro

    PROCEDURE, PUBLIC:: compute_and_export_SPH_variables

    PROCEDURE, PUBLIC:: read_sphincs_dump_print_formatted

    PROCEDURE, PUBLIC:: print_formatted_lorene_id_particles

    PROCEDURE, PUBLIC:: is_empty

    !PROCEDURE, PUBLIC:: write_lorene_bns_id_dump

    !-----------------!
    !--  FUNCTIONS  --!
    !-----------------!

    PROCEDURE, PUBLIC:: get_npart
    PROCEDURE, PUBLIC:: get_pos
    PROCEDURE, PUBLIC:: get_vel
    PROCEDURE, PUBLIC:: get_nlrf
    PROCEDURE, PUBLIC:: get_nu
    PROCEDURE, PUBLIC:: get_u
    PROCEDURE, PUBLIC:: get_pressure
    PROCEDURE, PUBLIC:: get_pressure_cu
    PROCEDURE, PUBLIC:: get_theta
    PROCEDURE, PUBLIC:: get_h

    ! Destructor
    FINAL:: destruct_particles

  END TYPE particles

  !
  !-- Interface of the TYPE particles (i.e., declaration of the constructor)
  !-- Multiple procedures in this interface would overload the constructor.
  !-- Such procedures must have distingushable interfaces, in particular
  !-- distinguishable arguments)
  !
  INTERFACE particles

    MODULE PROCEDURE construct_particles

  END INTERFACE particles

  !
  !-- Interface of the constructor of TYPE particles
  !-- Its implementation is in submodule_particles_constructor.f90
  !
  INTERFACE

    MODULE FUNCTION construct_particles( bns_obj, dist ) RESULT ( parts_obj )
        CLASS(bns), INTENT( IN OUT ):: bns_obj
        INTEGER,    INTENT( IN )    :: dist
        TYPE(particles)             :: parts_obj

    END FUNCTION construct_particles

   !MODULE FUNCTION construct_particles_empty() &
   !                    RESULT ( parts_sl_obj )
   !    TYPE(particles)             :: parts_sl_obj
   !
   !END FUNCTION construct_particles_empty

  END INTERFACE

  !
  !-- Interfaces of the methods of TYPE particles called by its constructor
  !-- Their implementations are in submodule_particles_constructor.f90
  !
  INTERFACE

    MODULE SUBROUTINE place_particles_3dlattice( THIS, &
                                  xmin, xmax, ymin, ymax, zmin, zmax, &
                                  thres, bns_obj )

      CLASS(particles), INTENT( IN OUT ):: THIS
      CLASS(bns),       INTENT( IN OUT ):: bns_obj
      DOUBLE PRECISION, INTENT( IN )    :: xmin, xmax, ymin, &
                                           ymax, zmin, zmax, thres

    END SUBROUTINE place_particles_3dlattice

    MODULE SUBROUTINE place_particles_3dlattices( THIS, &
                                  xmin1, xmax1, ymin1, ymax1, zmin1, zmax1, &
                                  xmin2, xmax2, ymin2, ymax2, zmin2, zmax2, &
                                  thres, bns_obj )

      CLASS(particles), INTENT( IN OUT ):: THIS
      CLASS(bns),       INTENT( IN OUT ):: bns_obj
      DOUBLE PRECISION, INTENT( IN )    :: xmin1, xmax1, ymin1, &
                                           ymax1, zmin1, zmax1
      DOUBLE PRECISION, INTENT( IN )    :: xmin2, xmax2, ymin2, &
                                           ymax2, zmin2, zmax2, thres

    END SUBROUTINE place_particles_3dlattices

    MODULE SUBROUTINE allocate_lorene_id_parts_memory( THIS )

      CLASS(particles), INTENT( IN OUT ):: THIS

    END SUBROUTINE allocate_lorene_id_parts_memory

  END INTERFACE

  !
  !-- Interfaces of the methods of TYPE particles
  !-- Their implementations are in module_particles_methods.f90
  !
  INTERFACE

    MODULE SUBROUTINE analyze_hydro( THIS, namefile )

      CLASS(particles),    INTENT( IN OUT ):: THIS
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE analyze_hydro

    MODULE SUBROUTINE compute_and_export_SPH_variables( THIS, namefile )

      CLASS(particles),    INTENT( IN OUT ):: THIS
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE compute_and_export_SPH_variables

    MODULE SUBROUTINE read_sphincs_dump_print_formatted( THIS, namefile_bin, &
                                                               namefile )

      CLASS(particles),    INTENT( IN OUT )           :: THIS
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile_bin
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE read_sphincs_dump_print_formatted

    MODULE SUBROUTINE print_formatted_lorene_id_particles( THIS, namefile )

      CLASS(particles),    INTENT( IN OUT )           :: THIS
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE print_formatted_lorene_id_particles

    MODULE SUBROUTINE destruct_particles( THIS )

      TYPE(particles), INTENT( IN OUT ):: THIS

    END SUBROUTINE destruct_particles

    MODULE FUNCTION is_empty( THIS ) RESULT( answer )

      CLASS(particles), INTENT( IN ):: THIS
      LOGICAL:: answer

    END FUNCTION is_empty

   !MODULE SUBROUTINE write_lorene_bns_id_dump( THIS, namefile )
   !
   !    CLASS(particles),    INTENT( IN )               :: THIS
   !    CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile
   !
   !END SUBROUTINE write_lorene_bns_id_dump

    MODULE FUNCTION get_npart( THIS ) RESULT( n_part )

      ! Arguments
      CLASS(particles), INTENT( IN OUT ):: THIS
      ! Result
      INTEGER:: n_part

    END FUNCTION get_npart

    MODULE FUNCTION get_pos( THIS ) RESULT( pos_u )

      ! Arguments
      CLASS(particles), INTENT( IN OUT ):: THIS
      ! Result
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos_u

    END FUNCTION get_pos

    MODULE FUNCTION get_vel( THIS ) RESULT( vel )

      ! Arguments
      CLASS(particles), INTENT( IN OUT ):: THIS
      ! Result
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: vel

    END FUNCTION get_vel

    MODULE FUNCTION get_nlrf( THIS ) RESULT( nlrf )

      ! Arguments
      CLASS(particles), INTENT( IN OUT ):: THIS
      ! Result
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: nlrf

    END FUNCTION get_nlrf

    MODULE FUNCTION get_nu( THIS ) RESULT( nu )

      ! Arguments
      CLASS(particles), INTENT( IN OUT ):: THIS
      ! Result
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: nu

    END FUNCTION get_nu

    MODULE FUNCTION get_u( THIS ) RESULT( u )

      ! Arguments
      CLASS(particles), INTENT( IN OUT ):: THIS
      ! Result
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: u

    END FUNCTION get_u

    MODULE FUNCTION get_pressure( THIS ) RESULT( pressure )

      ! Arguments
      CLASS(particles), INTENT( IN OUT ):: THIS
      ! Result
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: pressure

    END FUNCTION get_pressure

    MODULE FUNCTION get_pressure_cu( THIS ) RESULT( pressure_cu )

      ! Arguments
      CLASS(particles), INTENT( IN OUT ):: THIS
      ! Result
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: pressure_cu

    END FUNCTION get_pressure_cu

    MODULE FUNCTION get_theta( THIS ) RESULT( theta )

      ! Arguments
      CLASS(particles), INTENT( IN OUT ):: THIS
      ! Result
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: theta

    END FUNCTION get_theta

    MODULE FUNCTION get_h( THIS ) RESULT( h )

      ! Arguments
      CLASS(particles), INTENT( IN OUT ):: THIS
      ! Result
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: h

    END FUNCTION get_h

  END INTERFACE

END MODULE particles_id
