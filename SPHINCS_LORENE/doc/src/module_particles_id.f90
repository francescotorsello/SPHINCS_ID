! File:         module_particles_id.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

MODULE particles_id


  !***********************************************************
  !
  !# This module contains the definition of TYPE particles
  !
  !***********************************************************


  USE utility, ONLY: itr, ios, err_msg, test_status, &
                     perc, creturn, run_id, show_progress
  USE bns_id,  ONLY: bns
  USE timing,  ONLY: timer


  IMPLICIT NONE


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
  !! TYPE representing a particle distribution


    PRIVATE

    INTEGER:: npart
    !! Total particle number
    INTEGER:: npart1
    !! Particle number for star 1
    INTEGER:: npart2
    !! Particle number for star 2
    INTEGER:: npart_temp, npart1_temp, npart2_temp
    !! Various particle numbers used internally by the TYPE
    !INTEGER:: nx, ny, nz, nx1, ny1, nz1, nx2, ny2, nz2
    INTEGER:: distribution_id
    !! Identification number for the particle distribution
    INTEGER:: eos1_id
    !! LORENE identification number for the EOS of star 1
    INTEGER:: eos2_id
    !! LORENE identification number for the EOS of star 1
    INTEGER:: call_flag= 0
    ! Flag that is set different than 0 if the SUBROUTINE
    ! compute_and_export_SPH_variables is called

    INTEGER, DIMENSION(:), ALLOCATABLE:: baryon_density_index
    !# Array storing the indices to use with [[particles:baryon_density_parts]]
    !  to sort the elements of [[particles:baryon_density_parts]] in increasing
    !  order

    !INTEGER, DIMENSION(:), ALLOCATABLE:: filt_pos

    !
    !-- Hydro variables on the particles
    !
    !> 2-D array storing the particle positions
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos
    !> 1-D array storing the position of the particles on the x axis for S 1
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pos_x1
    !> 1-D array storing the position of the particles on the x axis for NS2
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pos_x2
    !& 1-D array storing the baryon mass density in the fluid frame
    !  \([\mathrm{kg}\,\mathrm{m}^{-3}]\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: baryon_density_parts
    !> 1-D array storing the energy density
    !  \([\mathrm{kg}\,c^2\,\mathrm{m}^{-3}]\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: energy_density_parts
    !> 1-D array storing the specific internal energy \([c^2]\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: specific_energy_parts
    !& 1-D array storing the specific internal energy \([c^2]\) computed using
    !  formula (9) in Read et al., Phys.Rev.D79:124032,2009,
    !  [arXiv:0812.2163][https://arxiv.org/abs/0812.2163]{:target="_blank"}
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: u_pwp
    !> 1-D array storing the pressure \([\mathrm{kg}\,c^2\,\mathrm{m}^{-3}]\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_parts
    !& 1-D array storing the pressure on the x axis
    !  \([\mathrm{kg}\,c^2\,\mathrm{m}^{-3}]\) for NS 1
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_parts_x1
    !& 1-D array storing the pressure on the x axis
    !  \([\mathrm{kg}\,c^2\,\mathrm{m}^{-3}]\) for NS 2
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_parts_x2
    !& 1-D array storing the first derivative of the pressure
    !  along the x axis \([\mathrm{kg}\,c^2\,\mathrm{m}^{-3}]\) for NS 1
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_parts_x_der1
    !& 1-D array storing the first derivative of the pressure
    !  along the x axis \([\mathrm{kg}\,c^2\,\mathrm{m}^{-3}]\) for NS2
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_parts_x_der2
    !> 1-D array storing the typical length scale for the pressure change
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_length_scale_x1
    !> 1-D array storing the typical length scale for the pressure change
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_length_scale_x2
    !& 1-D array storing the pressure in code units
    !  \([\mathrm{amu}\,c^2\,\mathrm{L_\odot}^{-3}]\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_parts_cu
    !& 1-D array storing the x component of the fluid 3-velocity wrt
    !  the Eulerian observer \([c]\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: v_euler_parts_x
    !& 1-D array storing the y component of the fluid 3-velocity wrt
    !  the Eulerian observer \([c]\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: v_euler_parts_y
    !& 1-D array storing the z component of the fluid 3-velocity wrt
    !  the Eulerian observer \([c]\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: v_euler_parts_z

    !
    !-- Arrays to store the electron fraction Ye as a function of the
    !-- baryon number density for beta-equilibrated EoS at T~0,
    !-- imported from the CompOSE database's and software's files
    !

    !& Array storing the values of the baryon number density in the CompOSE
    !  table. @todo ADD UNITS
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nb_table
    !> Array storing the values of the electron fraction in the CompOSE table
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: Ye_table

    !
    !-- Spacetime fields
    !

    !> Array storing the values of the lapse function on the particles
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: lapse_parts
    !& Array storing the values of the x component of the shift vector
    !  on the particles
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: shift_parts_x
    !& Array storing the values of the y component of the shift vector
    !  on the particles
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: shift_parts_y
    !& Array storing the values of the z component of the shift vector
    !  on the particles
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: shift_parts_z
    !& Array storing the values of the xx component of the spatial metric
    !  on the particles
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: g_xx_parts
    !& Array storing the values of the xy component of the spatial metric
    !  on the particles
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: g_xy_parts
    !& Array storing the values of the xz component of the spatial metric
    !  on the particles
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: g_xz_parts
    !& Array storing the values of the xz component of the spatial metric
    !  on the particles
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: g_yy_parts
    !& Array storing the values of the yz component of the spatial metric
    !  on the particles
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: g_yz_parts
    !& Array storing the values of the zz component of the spatial metric
    !  on the particles
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: g_zz_parts

    !
    !-- SPH fields
    !

    !& 1-D array storing baryon density in the local rest frame
    !  \([\mathrm{baryon}\, (L_\odot)^{-3}]\), computed directly from
    !  the LORENE density
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nlrf
    !& 1-D array storing baryon density in the local rest frame
    !  \([\mathrm{baryon}\, (L_\odot)^{-3}]\), computed from the kernel
    !  interpolated proper baryon number density [[particles:nstar_int]]
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nlrf_int
    !> 1-D array storing the baryon number per particle
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nu
    !& 1-D array storing the SPH estimate of the proper baryon number density
    !  \([\mathrm{baryon}\, (L_\odot)^{-3}]\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar
    !& 1-D array storing the particle number density
    !  \([\mathrm{particle}\, (L_\odot)^{-3}]\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: particle_density
    !& 1-D array storing the SPH estimate of the proper baryon number density,
    !  from kernel interpolation \([\mathrm{baryon}\, (L_\odot)^{-3}]\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar_int
    !& 1-D array storing the SPH estimate of the particle number density, from
    !  kernel interpolation \([\mathrm{particle}\, (L_\odot)^{-3}]\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: particle_density_int
    !> 2-D array storing the coordinate fluid 4-velocity \([c]\)
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: v
    !> 1-D array storing the generalized Lorentz factor
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: Theta
    !> 1-D array storing the electron fraction
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: Ye
    !> 1-D array storing the smoothing length \(L_\odot\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: h
    !> 1-D array storing the particle volumes \(L_\odot^3\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pvol
    !> 1-D array storing the particle masses \(M_\odot\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pmass
    !> Baryonic mass of of star 1 \(M_\odot\)
    DOUBLE PRECISION:: mass1
    !> Baryonic mass of of star 2 \(M_\odot\)
    DOUBLE PRECISION:: mass2
    !& Ratio of baryonic masses of the stars \(M_\odot\)
    !  @warning always \(< 1\)
    DOUBLE PRECISION:: mass_ratio
    !> Total grid volume
    !DOUBLE PRECISION:: vol, vol1, vol2
    !> Volume per particle
    !DOUBLE PRECISION:: vol_a, vol1_a, vol2_a
    !> Ratio between the max and min of the baryon number per particle
    DOUBLE PRECISION:: nu_ratio
    !> Total baryon number
    DOUBLE PRECISION:: nbar_tot
    !> Baryon number on star 1
    DOUBLE PRECISION:: nbar1
    !> Baryon number on star 2
    DOUBLE PRECISION:: nbar2
    !> Baryon number ratio on both stars
    DOUBLE PRECISION:: nuratio
    !> Baryon number ratio on star 1
    DOUBLE PRECISION:: nuratio1
    !> Baryon number ratio on star 2
    DOUBLE PRECISION:: nuratio2

    !CHARACTER( LEN= 50 ):: lorene_bns_id_parfile

    !> String storing the local path to the directory containing the CompOSE EOS
    CHARACTER( LEN= : ), ALLOCATABLE:: compose_path
    !> String storing the subpath of compose_path to the CompOSE file with
    !  .beta extension
    CHARACTER( LEN= : ), ALLOCATABLE:: compose_filename

    CHARACTER( LEN= : ), ALLOCATABLE:: eos1
    CHARACTER( LEN= : ), ALLOCATABLE:: eos2

    LOGICAL:: empty_object
    LOGICAL, PUBLIC:: export_bin
    LOGICAL, PUBLIC:: export_form_xy, export_form_x
    LOGICAL:: use_thres
    LOGICAL:: redistribute_nu
    LOGICAL:: correct_nu
    LOGICAL:: compose_eos
    LOGICAL:: randomize_phi, randomize_theta, randomize_r
    LOGICAL:: apm_iterate1, apm_iterate2
    LOGICAL:: read_nu
    LOGICAL:: reflect_particles_x

    TYPE(timer), PUBLIC:: placer_timer
    TYPE(timer), PUBLIC:: same_particle_timer
    TYPE(timer), PUBLIC:: apm1_timer
    TYPE(timer), PUBLIC:: apm2_timer
    TYPE(timer), PUBLIC:: importer_timer
    TYPE(timer), PUBLIC:: sph_computer_timer


    CONTAINS


    !-------------------!
    !--  SUBROUTINES  --!
    !-------------------!

    PROCEDURE:: place_particles_3dlattice

    PROCEDURE:: place_particles_3dlattices

    PROCEDURE:: place_particles_spherical_shells

    PROCEDURE:: perform_apm

    GENERIC, PUBLIC:: reshape_sph_field => reshape_sph_field_1d_ptr, &
                                           reshape_sph_field_2d_ptr
    PROCEDURE:: reshape_sph_field_1d_ptr => reshape_sph_field_1d
    PROCEDURE:: reshape_sph_field_2d_ptr => reshape_sph_field_2d


    PROCEDURE:: allocate_lorene_id_parts_memory

    PROCEDURE:: deallocate_lorene_id_parts_memory

    PROCEDURE:: read_compose_composition

    PROCEDURE:: compute_Ye

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
    PROCEDURE, PUBLIC:: get_npart1
    PROCEDURE, PUBLIC:: get_npart2
    PROCEDURE, PUBLIC:: get_nuratio
    PROCEDURE, PUBLIC:: get_nuratio1
    PROCEDURE, PUBLIC:: get_nuratio2
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


    !-------------------!
    !--  SUBROUTINES  --!
    !-------------------!


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


    MODULE SUBROUTINE place_particles_spherical_shells( THIS, &
                                  mass_star, radius, center, npart_approx, &
                                  npart_out, pos, pvol, pmass, bns_obj, &
                                  last_r, upper_bound, lower_bound, &
                                  upper_factor, lower_factor, max_steps, &
                                  filename_mass_profile, filename_shells_radii,&
                                  filename_shells_pos )

      CLASS(particles), INTENT( IN OUT ):: THIS
      CLASS(bns),       INTENT( IN OUT ):: bns_obj
      INTEGER,          INTENT( IN )    :: npart_approx, max_steps
      INTEGER,          INTENT( OUT )   :: npart_out
      DOUBLE PRECISION, INTENT( IN )    :: mass_star, radius, center, &
                                           last_r
      DOUBLE PRECISION, INTENT( INOUT ) :: upper_bound, lower_bound
      DOUBLE PRECISION, INTENT( IN )    :: upper_factor, lower_factor
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT( OUT ):: pos
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, INTENT( OUT ):: pvol
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, INTENT( OUT ):: pmass
      CHARACTER( LEN= * ), INTENT( INOUT ), OPTIONAL :: filename_mass_profile
      CHARACTER( LEN= * ), INTENT( INOUT ), OPTIONAL :: filename_shells_radii
      CHARACTER( LEN= * ), INTENT( INOUT ), OPTIONAL :: filename_shells_pos

    END SUBROUTINE place_particles_spherical_shells


    MODULE SUBROUTINE reshape_sph_field_1d( THIS, field, new_size1, new_size2, &
                                            index_array )

      CLASS(particles), INTENT( IN OUT ):: THIS
      INTEGER,                        INTENT( IN ):: new_size1
      INTEGER,                        INTENT( IN ):: new_size2
      INTEGER,          DIMENSION(:), INTENT( IN ):: index_array
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT( IN OUT ):: field

    END SUBROUTINE reshape_sph_field_1d


    MODULE SUBROUTINE reshape_sph_field_2d( THIS, field, new_size1, new_size2, &
                                            index_array )

      CLASS(particles), INTENT( IN OUT ):: THIS
      INTEGER,                        INTENT( IN ):: new_size1
      INTEGER,                        INTENT( IN ):: new_size2
      INTEGER,          DIMENSION(:), INTENT( IN ):: index_array
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT( IN OUT ):: field

    END SUBROUTINE reshape_sph_field_2d


    MODULE SUBROUTINE allocate_lorene_id_parts_memory( THIS )

      CLASS(particles), INTENT( IN OUT ):: THIS

    END SUBROUTINE allocate_lorene_id_parts_memory


    MODULE SUBROUTINE deallocate_lorene_id_parts_memory( THIS )

      CLASS(particles), INTENT( IN OUT ):: THIS

    END SUBROUTINE deallocate_lorene_id_parts_memory


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

    MODULE SUBROUTINE perform_apm( THIS, &
                                   binary, &
                                   pos_input, &
                                   pvol, h_output, nu_output, &
                                   center, &
                                   com_star, &
                                   mass, &
                                   apm_max_it, max_inc, &
                                   mass_it, correct_nu, nuratio_thres, &
                                   nuratio_des, &
                                   nx_gh, ny_gh, nz_gh, &
                                   namefile_pos_id, namefile_pos, &
                                   namefile_results )

      CLASS(particles),                 INTENT( INOUT ):: THIS
      CLASS(bns),                       INTENT( INOUT ):: binary
      DOUBLE PRECISION, DIMENSION(:,:), INTENT( INOUT ):: pos_input
      DOUBLE PRECISION, DIMENSION(:),   INTENT( INOUT ):: pvol
      DOUBLE PRECISION, DIMENSION(:),   INTENT( OUT )  :: h_output
      DOUBLE PRECISION, DIMENSION(:),   INTENT( OUT )  :: nu_output
      DOUBLE PRECISION,                 INTENT( IN )   :: center
      DOUBLE PRECISION,                 INTENT( IN )   :: com_star
      DOUBLE PRECISION,                 INTENT( IN )   :: mass
      INTEGER,                          INTENT( IN )   :: apm_max_it
      INTEGER,                          INTENT( IN )   :: max_inc
      LOGICAL,                          INTENT( IN )   :: mass_it
      LOGICAL,                          INTENT( IN )   :: correct_nu
      DOUBLE PRECISION,                 INTENT( IN )   :: nuratio_thres
      DOUBLE PRECISION,                 INTENT( IN )   :: nuratio_des
      INTEGER,                          INTENT( IN )   :: nx_gh, ny_gh, nz_gh
      CHARACTER( LEN= * ),              INTENT( INOUT ), OPTIONAL :: &
                                                            namefile_pos_id
      CHARACTER( LEN= * ),              INTENT( INOUT ), OPTIONAL :: &
                                                            namefile_pos
      CHARACTER( LEN= * ),              INTENT( INOUT ), OPTIONAL :: &
                                                            namefile_results

    END SUBROUTINE perform_apm

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

    MODULE SUBROUTINE read_compose_composition( THIS, namefile )

      CLASS(particles),    INTENT( IN OUT )           :: THIS
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE read_compose_composition

    MODULE SUBROUTINE compute_Ye( THIS )!, nlrf, Ye )

      CLASS(particles),    INTENT( IN OUT )           :: THIS
      !DOUBLE PRECISION, DIMENSION( : ), INTENT( IN ):: nlrf
      !DOUBLE PRECISION, DIMENSION( : ), INTENT( OUT ):: Ye

    END SUBROUTINE compute_Ye

    MODULE SUBROUTINE destruct_particles( THIS )

      TYPE(particles), INTENT( IN OUT ):: THIS

    END SUBROUTINE destruct_particles

    !-----------------!
    !--  FUNCTIONS  --!
    !-----------------!

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

    MODULE FUNCTION get_npart1( THIS ) RESULT( n_part )

      ! Arguments
      CLASS(particles), INTENT( IN OUT ):: THIS
      ! Result
      INTEGER:: n_part

    END FUNCTION get_npart1

    MODULE FUNCTION get_npart2( THIS ) RESULT( n_part )

      ! Arguments
      CLASS(particles), INTENT( IN OUT ):: THIS
      ! Result
      INTEGER:: n_part

    END FUNCTION get_npart2

    MODULE FUNCTION get_nuratio( THIS ) RESULT( nuratio )

      ! Arguments
      CLASS(particles), INTENT( IN OUT ):: THIS
      ! Result
      DOUBLE PRECISION:: nuratio

    END FUNCTION get_nuratio

    MODULE FUNCTION get_nuratio1( THIS ) RESULT( nuratio1 )

      ! Arguments
      CLASS(particles), INTENT( IN OUT ):: THIS
      ! Result
      DOUBLE PRECISION:: nuratio1

    END FUNCTION get_nuratio1

    MODULE FUNCTION get_nuratio2( THIS ) RESULT( nuratio2 )

      ! Arguments
      CLASS(particles), INTENT( IN OUT ):: THIS
      ! Result
      DOUBLE PRECISION:: nuratio2

    END FUNCTION get_nuratio2

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



  CONTAINS



  SUBROUTINE check_particle_positions( npart, pos, debug )

    !*************************************************
    !                                                *
    ! Check that the particles do not have the same  *
    ! positions                                      *
    !                                                *
    ! FT 1.9.2021                                    *
    !                                                *
    !*************************************************

    USE NR,             ONLY: indexx

    IMPLICIT NONE

    INTEGER, INTENT(IN):: npart
    LOGICAL, INTENT(IN), OPTIONAL:: debug
    DOUBLE PRECISION, DIMENSION(3,npart), INTENT(IN):: pos

    INTEGER:: itr, itr2, x_idx
    INTEGER, DIMENSION(npart):: x_sort
    INTEGER, DIMENSION(:), ALLOCATABLE:: x_number

    PRINT *, "** Checking that there are not multiple particles", &
             " at the same position..."
    PRINT *

    ALLOCATE( x_number( npart ) )

    ! Sort x coordinates of the particles
    CALL indexx( npart, pos( 1, : ), x_sort )

    x_number= 1
    itr2= 1
    ! Find the number of times each x appears
    DO itr= 1, npart - 1, 1

      IF( pos( 1, x_sort(itr) ) == &
          pos( 1, x_sort(itr+1) ) )THEN

        x_number(itr2)= x_number(itr2) + 1

      ELSE

        itr2= itr2 + 1

      ENDIF

    ENDDO
    x_number= x_number(1:itr2)

    IF( SUM( x_number ) /= npart )THEN

      PRINT *, "** ERROR! The sum of the numbers of particles with the same", &
               " x is not equal to the particle number."
      PRINT *, " * SUM( x_number )=", SUM( x_number ), ", ", &
               "npart=", npart
      PRINT *, " * Stopping..."
      PRINT *
      STOP

    ENDIF

    IF( PRESENT(debug) .AND. debug == .TRUE. )THEN

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( pos, x_sort, x_number ) &
      !$OMP             PRIVATE( itr, itr2, x_idx )
      DO itr= 1, SIZE(x_number), 1

        IF( itr == 1 )THEN
          x_idx= 1
        ELSE
          x_idx= SUM(x_number(1:itr-1)) + 1
        ENDIF

        DO itr2= x_idx, x_idx + x_number(itr) - 2, 1

          ! If they do not have the same x
          IF( pos( 1, x_sort(itr2) ) /= &
              pos( 1, x_sort(itr2+1) ) )THEN

            PRINT *, "** ERROR! ", "The two particles ", x_sort(itr2), &
                     " and", x_sort(itr2+1), &
                     " do not have the same x, but should!"
            PRINT *, pos( :, x_sort(itr2) )
            PRINT *, pos( :, x_sort(itr2+1) )
            PRINT *, " * Stopping..."
            PRINT *
            STOP

          ENDIF

        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

    ENDIF

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( pos, x_sort, x_number ) &
    !$OMP             PRIVATE( itr, itr2, x_idx )
    DO itr= 1, SIZE(x_number), 1

      IF( itr == 1 )THEN
        x_idx= 1
      ELSE
        x_idx= SUM(x_number(1:itr-1)) + 1
      ENDIF

      DO itr2= x_idx, x_idx + x_number(itr) - 2, 1

        ! If they have the same y
        IF( pos( 2, x_sort(itr2) ) == &
            pos( 2, x_sort(itr2+1) ) )THEN

          ! If they have the same z
          IF( pos( 3, x_sort(itr2) ) == &
              pos( 3, x_sort(itr2+1) ) )THEN

            ! They are the same
            PRINT *, "** ERROR! ", "The two particles ", x_sort(itr2), &
                     " and", x_sort(itr2+1), " have the same coordinates!"
            PRINT *, pos( :, x_sort(itr2) )
            PRINT *, pos( :, x_sort(itr2+1) )
            PRINT *, " * Stopping..."
            PRINT *
            STOP

          ENDIF
        ENDIF

      ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    DEALLOCATE( x_number )

  END SUBROUTINE check_particle_positions


END MODULE particles_id
