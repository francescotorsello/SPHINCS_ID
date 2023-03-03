! File:         module_sph_particles.f90
! Authors:      Francesco Torsello (FT)
!************************************************************************
! Copyright (C) 2020-2023 Francesco Torsello                            *
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

MODULE sph_particles


  !***********************************************************
  !
  !# This module contains the definition of TYPES,
  !  PROCEDURES and variables needed to set up the |sph| |id|
  !
  !***********************************************************


  USE utility,                  ONLY: itr, ios, err_msg, test_status, &
                                      is_finite_number, perc, creturn, &
                                      run_id, show_progress
  USE id_base,                  ONLY: idbase, surface
  USE timing,                   ONLY: timer


  IMPLICIT NONE


  INTEGER, PARAMETER:: id_particles_from_formatted_file  = 0
  !! Identifier for a particle distribution read from formatted file
  INTEGER, PARAMETER:: id_particles_on_lattice           = 1
  !! Identifier for a particle distribution on a lattice
  INTEGER, PARAMETER:: id_particles_on_ellipsoidal_surfaces= 2
  !! Identifier for particle distribution on ellipsoidal surfaces

  INTEGER, PARAMETER:: max_it_tree= 1
  !# When computing the neighbours' tree with the SUBROUTINE
  !  exact_nei_tree_update, it can happen that some particles do not have
  !  exactly 300 neighours. If this happens, the smoothing lenghts of the
  !  particles in the tree-cell are increased by a factor of 3, and used
  !  as a new guess to recompute the entire tree.
  !  [[sph_particles:max_it_tree]] specifies how many times the process should
  !  be iterated, if needed.
  !
  !  Note that, at the end of the iteration, the smoothing lengths are checked
  !  for them being NaNs or 0. For the particles at which this happens,
  !  the smoothing lengths are computed brute-force using the SUBROUTINE
  !  [[sph_particles:find_h_backup]], so that such particles have exactly
  !  300 neighbours.
  !
  !  @warning This is DEPRECATED because the problem in the neighbor's tree
  !           was solved.


  !
  !-- TODO: Make the following parameters member of TYPE eos?
  !--       They belong there, but putting them there decreases readability
  !--       of the code
  !--       Rewrite this as an enumeration?
  !--       How to write enumerations in Fortran?
  !
  INTEGER,               PARAMETER:: eos$id     = 1
  !! First component of array [[eos:eos_parameters]] for a polytropic |eos|
  INTEGER,               PARAMETER:: poly$gamma = 2
  !! Second component of array [[eos:eos_parameters]] for a polytropic |eos|
  INTEGER,               PARAMETER:: poly$kappa = 3
  !! Third component of array [[eos:eos_parameters]] for a polytropic |eos|
  INTEGER,               PARAMETER:: pwp$npoly  = 2
  !# Second component of array [[eos:eos_parameters]] for a piecewise
  !  polytropic |eos|
  INTEGER, DIMENSION(4), PARAMETER:: pwp$gamma  = [3,4,5,6]
  !# Third to sixth component of array [[eos:eos_parameters]] for a piecewise
  !  polytropic |eos|
  INTEGER, DIMENSION(4), PARAMETER:: pwp$kappa  = [7,8,9,10]
  !# Seventh to tenth component of array [[eos:eos_parameters]] for a piecewise
  !  polytropic |eos|
  INTEGER,               PARAMETER:: pwp$log10p1  = 11
  !# Eleventh component of array [[eos:eos_parameters]] for a piecewise
  !  polytropic |eos|
  INTEGER,               PARAMETER:: pwp$log10rho0= 12
  !# Twelfth component of array [[eos:eos_parameters]] for a piecewise
  !  polytropic |eos|
  INTEGER,               PARAMETER:: pwp$log10rho1= 13
  !# Thirteenth component of array [[eos:eos_parameters]] for a piecewise
  !  polytropic |eos|
  INTEGER,               PARAMETER:: pwp$log10rho2= 14
  !# Fourteenth component of array [[eos:eos_parameters]] for a piecewise
  !  polytropic |eos|

  TYPE:: eos
  !! Data structure representing an |eos|

    CHARACTER(LEN=:), ALLOCATABLE:: eos_name
    !! The |eos| name
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: eos_parameters
    !# Array containing the |eos| parameters, in the following order:
    !
    !  - For a polytropic |eos|:
    !    [ eos identification number, \(\gamma\),\(\kappa\) ]
    !  - For a piecewise polytropic |eos|:
    !    [ eos identification number, number of polytropic pieces, \(\gamma_0\),
    !   \(\gamma_1\), \(\gamma_2\), \(\gamma_3\), \(\kappa_0\), \(\kappa_1\),
    !    \(\kappa_2\), \(\kappa_3\), \(\log_{10}(p_1)\), \(\log_{10}(\rho_0)\),
    !    \(\log_{10}(\rho_1)\), \(\log_{10}(\rho_2)\) ]
    !  - For a tabulated |eos|: [ eos identification number ]
    !
    ! \(\gamma\) is the adimensional polytropic exponent, \(\kappa\) is the
    ! polytropic constant in units of
    ! \(\left(M_\odot L_\odot^{-3}\right)^{1-\gamma}\). Pressure and baryon
    ! mass density have the same units \(M_\odot L_\odot^{-3}\) since \(c^2=1\).

  END TYPE eos


  !***********************************************************
  !                                                          *
  !              Definition of TYPE particles                *
  !                                                          *
  ! This class places the |sph| particles, imports           *
  ! the |id| on the particle positions, stores               *
  ! it, computes the relevant |sph| fields and exports it to *
  ! both a formatted, and a binary file for evolution        *
  !                                                          *
  !***********************************************************

  TYPE:: particles
  !! TYPE representing a |sph| particle distribution


    PRIVATE


    INTEGER:: npart
    !! Total particle number
    INTEGER:: n_matter
    !! Number of matter objects in the physical system
    INTEGER, DIMENSION(:), ALLOCATABLE:: npart_i
    !! Array storing the particle numbers for the matter objects
    INTEGER, DIMENSION(:), ALLOCATABLE:: npart_fin
    !! Array storing the index of the last particle for each matter objects
    INTEGER:: distribution_id
    !! Identification number for the particle distribution
    INTEGER:: call_flag= 0
    !# Flag that is set different than 0 if the SUBROUTINE
    !  compute_and_print_sph_variables is called
    LOGICAL:: cold_system
    !# `.TRUE.` if the system is at zero temperature (no thermal component);
    !  `.FALSE.` otherwise
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: barycenter
    !# Array storing the centers of mass of the matter objects
    DOUBLE PRECISION, DIMENSION(4):: barycenter_system
    !# Array storing the center of mass of the **entire particle distribution**

    !DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: surfaces
    !# Array storing the surfaces of the matter objects

    INTEGER, DIMENSION(:), ALLOCATABLE:: baryon_density_index
    !# Array storing the indices to use with [[particles:baryon_density]]
    !  to sort the elements of [[particles:baryon_density]] in increasing
    !  order

    !
    !-- Hydro variables on the particles
    !

    !> 2-D array storing the particle positions
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos
    !& 1-D array storing the baryon mass density in the fluid frame
    !  \([\mathrm{kg}\,\mathrm{m}^{-3}]\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: baryon_density
    !> 1-D array storing the energy density
    !  \([\mathrm{kg}\,c^2\,\mathrm{m}^{-3}]\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: energy_density
    !> 1-D array storing the specific internal energy \([c^2]\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: specific_energy
    !& 1-D array storing the specific internal energy \([c^2]\) computed using
    !  formula (9) in Read et al., Phys.Rev.D79:124032,2009,
    !  [arXiv:0812.2163][https://arxiv.org/abs/0812.2163]{:target="_blank"}
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: u_sph
    !> 1-D array storing the pressure \([\mathrm{kg}\,c^2\,\mathrm{m}^{-3}]\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure
    !& 1-D array storing the pressure in code units
    !  \([\mathrm{amu}\,c^2\,\mathrm{L_\odot}^{-3}]\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_sph
    !& 1-D array storing the x component of the fluid 3-velocity wrt
    !  the Eulerian observer \([c]\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: v_euler_x
    !& 1-D array storing the y component of the fluid 3-velocity wrt
    !  the Eulerian observer \([c]\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: v_euler_y
    !& 1-D array storing the z component of the fluid 3-velocity wrt
    !  the Eulerian observer \([c]\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: v_euler_z

    !
    !-- Arrays to store the electron fraction Ye as a function of the
    !-- baryon number density for beta-equilibrated EoS at T~0,
    !-- imported from the CompOSE database's and software's files
    !

    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nb_table
    !# Array storing the values of the baryon number density in the |compose|
    !  table. @todo ADD UNITS
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: Ye_table
    !! Array storing the values of the electron fraction in the |compose| table

    !
    !-- Spacetime fields
    !

    !> Array storing the values of the lapse function on the particles
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: lapse
    !& Array storing the values of the x component of the shift vector
    !  on the particles
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: shift_x
    !& Array storing the values of the y component of the shift vector
    !  on the particles
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: shift_y
    !& Array storing the values of the z component of the shift vector
    !  on the particles
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: shift_z
    !& Array storing the values of the xx component of the spatial metric
    !  on the particles
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: g_xx
    !& Array storing the values of the xy component of the spatial metric
    !  on the particles
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: g_xy
    !& Array storing the values of the xz component of the spatial metric
    !  on the particles
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: g_xz
    !& Array storing the values of the xz component of the spatial metric
    !  on the particles
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: g_yy
    !& Array storing the values of the yz component of the spatial metric
    !  on the particles
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: g_yz
    !& Array storing the values of the zz component of the spatial metric
    !  on the particles
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: g_zz

    !
    !-- SPH fields
    !

    !& 1-D array storing baryon density in the local rest frame
    !  \([\mathrm{baryon}\, (L_\odot)^{-3}]\), computed directly from
    !  the |lorene| density
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nlrf
    !& 1-D array storing baryon density in the local rest frame
    !  \([\mathrm{baryon}\, (L_\odot)^{-3}]\), computed from the kernel
    !  interpolated proper baryon number density [[particles:nstar_sph]]
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nlrf_sph
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
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar_sph
    !& 1-D array storing the SPH estimate of the particle number density, from
    !  kernel interpolation \([\mathrm{particle}\, (L_\odot)^{-3}]\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: particle_density_sph
    !> 1-D array storing the enthalpy \([??]\) @todo add units
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: enthalpy
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
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: masses
    !> 1-D array storing the surfacesof the matter objects
    TYPE(surface), DIMENSION(:), ALLOCATABLE:: surfaces
    !& Ratio between baryonic masses of the matter objects and the maximum
    !  baryonic mass among them @warning always \(< 1\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: mass_ratios
    !& Ratio between baryonic masses of the matter objects and the total
    !  baryonic mass of all the matter objects @warning always \(< 1\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: mass_fractions
    !> Total baryon number
    DOUBLE PRECISION:: nbar_tot
    !> Baryon numbers of the matter objects
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nbar_i
    !& Ratio between the max and min of the baryon number per particle, over
    !  all the matter objects
    DOUBLE PRECISION:: nuratio
    !& Desired ratio between the max and min of the baryon number per particle,
    !  over all the matter objects. **Only used when redistribute_nu is .TRUE.**
    !  @warning Almost deprecated
    DOUBLE PRECISION:: nu_ratio_des
    !> Baryon number ratios on the matter objects
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nuratio_i

    DOUBLE PRECISION:: adm_mass
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: adm_linear_momentum_i
    !# Estimate of the \(\mathrm{ADM}\) linear momentum of each matter object,
    !  computed from the canonical momentum per baryon on the particles
    DOUBLE PRECISION, DIMENSION(3)               :: adm_linear_momentum_fluid
    !# Estimate of the \(\mathrm{ADM}\) linear momentum of the fluid computed
    !  from the canonical momentum per baryon on the particles

    !
    !-- Strings
    !

    CHARACTER( LEN= 50 ):: sphincs_id_particles
    !! String containing the name of the particles parameter file

    CHARACTER(LEN=:), ALLOCATABLE:: compose_path
    !# String storing the local path to the directory containing the
    !  CompOSE |eos|
    CHARACTER(LEN=:), ALLOCATABLE:: compose_filename
    !# String storing the subpath of compose_path to the CompOSE file with
    !  .beta extension

    !
    !-- Equations of state
    !

    TYPE(eos), DIMENSION(:), ALLOCATABLE:: all_eos
    !# Array of TYPE [[eos]] containing the |eos| information for all the matter
    !  objects

    !
    !-- Procedure pointers
    !

    !PROCEDURE(post_process_sph_id_int), POINTER, NOPASS:: post_process_sph_id
    PROCEDURE(), POINTER, NOPASS:: post_process_sph_id
    !# Pointer to a procedure that post_process the |sph| |id|; for example,
    !  correct for the residual ADM linear momentum.

    !
    !-- Steering variables
    !

    !> `.TRUE.` if the object is empty, `.FALSE.` if it's not empty
    LOGICAL:: empty_object
    !& `.TRUE.` if the binary files for SPHINCS_BSSN are to be exported,
    !  `.FALSE.` otherwise
    LOGICAL, PUBLIC:: export_bin
    !& `.TRUE.` if the ID in the formatted files is to be on the xy plane only,
    !  `.FALSE.` otherwise
    LOGICAL, PUBLIC:: export_form_xy
    !& `.TRUE.` if the ID in the formatted files is to be on the x axis only,
    !  `.FALSE.` otherwise
    LOGICAL, PUBLIC:: export_form_x
    !& `.TRUE.` if the threshold on the baryon mass density should e applied
    !  when placing particles on lattices, `.FALSE.` otherwise
    LOGICAL:: use_thres
    !& `.TRUE.` if the baryon number per particle should be reassigned, trying
    !  to obtain a baryon number ratio no larger than nu_ratio,
    !  when placing particles on lattices; `.FALSE.` otherwise
    LOGICAL:: redistribute_nu
    !& `.TRUE.` if the baryon number per particle should be corrected to account
    !  for the total baryon masses of the stars, `.FALSE.` otherwise
    LOGICAL:: correct_nu
    !& `.TRUE.` if the electron fraction \(Y_e\) should be read from the CompOSE
    !  table with extension.beta, `.FALSE.` otherwise
    !  @todo Chamge name of this variable to assign_Ye_compose. Check that
    !        the used EOS is indeed the one used to read \(Y_e\)
    LOGICAL:: compose_eos
    !& `.TRUE.` if the particle positions on ellipsoidal surfaces are randomized
    !  in the \(\phi\) direction, `.FALSE.` otherwise
    LOGICAL:: randomize_phi
    !& `.TRUE.` if the particle positions on ellipsoidal surfaces are randomized
    !  in the \(\theta\) direction, `.FALSE.` otherwise
    LOGICAL:: randomize_theta
    !& `.TRUE.` if the particle positions on ellipsoidal surfaces are randomized
    !  in the \(r\) direction, `.FALSE.` otherwise
    LOGICAL:: randomize_r
    !& `.TRUE.` if the Artificial Pressure Method (APM) has to be applied to the
    !  particles, `.FALSE.` otherwise
    LOGICAL, DIMENSION(:), ALLOCATABLE:: apm_iterate
    !& `.TRUE.` to allow the particles to move where the density is 0 during the
    !  Artificial Pressure Method (APM) iteration. This can be useful when the
    !  system has an irregular geometry, as, for example, an ejecta
    !  `.FALSE.` otherwise
    LOGICAL, DIMENSION(:), ALLOCATABLE:: use_atmosphere
    !& `.TRUE.` if the baryon number per particle \(\nu\) has to be read from
    !  the formatted file containing the particle positions, `.FALSE.` otherwise
    LOGICAL:: read_nu
    !& `.TRUE.` if the particles on star 2 should be the reflection of the
    !  particles on star 1 with respect to the \(yz\) plane, only if the baryon
    !  masses of the stars differe less than \(0.2\%\); `.FALSE.` otherwise
    LOGICAL:: reflect_particles_x

    !
    !-- Timers
    !

    !> Timer that times how long it takes to place particles on the stars
    TYPE(timer), PUBLIC:: placer_timer
    !& Timer that times how long it takes to check if there are multiple
    !  particles at the same positions
    TYPE(timer), PUBLIC:: same_particle_timer
    !& Timer that times how long it takes to perform the APM on the matter
    !  objects
    TYPE(timer), DIMENSION(:), ALLOCATABLE, PUBLIC:: apm_timers
    !& Timer that times how long it takes to import the \(\texttt{|lorene|}\) ID
    !  at the particle positions
    TYPE(timer), PUBLIC:: importer_timer
    !& Timer that times how long it takes to compute the SPH variables at the
    !  particle positions
    TYPE(timer), PUBLIC:: sph_computer_timer


    CONTAINS


    !-------------------!
    !--  SUBROUTINES  --!
    !-------------------!

    PROCEDURE:: read_particles_formatted_file
    !# Read particle positions and nu from a formatted file with the following
    !  format:
    !
    !   1. The first line must contain the integer number of objects for the
    !      system, and the particle numbers on each object;
    !      each column separated by one tab.
    !   2. The other lines must contain at least 3 columns with the x, y, z
    !      coordinates of the particle positions. An optional 4th column can
    !      contain the values of nu on the particles. If the 4th column is
    !      not present, nu will be roughly estimated using the density read
    !      from the ID and the average volume per particle; the latter is
    !      computed by computing the average particle separation along the
    !      z direction.

    PROCEDURE:: place_particles_lattice
    !! Places particles on a single lattice that surrounds both stars

    PROCEDURE:: place_particles_ellipsoidal_surfaces
    !! Places particles on ellipsoidal surfaces on one star

    PROCEDURE, NOPASS:: impose_equatorial_plane_symmetry
    !# Reflects the positions of the particles on a matter object with respect
    !  to the \(xy\) plane

    PROCEDURE:: get_object_of_particle
    !# Returns the number of the matter object asociated with the particle
    !  number given as input. Example: give number \(n\) as input; this
    !  particle number corresponds to a particle on matter object \(m\).
    !  This functions returns \(m\).

    PROCEDURE, NOPASS:: perform_apm
    !! Performs the Artificial Pressure Method (APM) on one star's particles

  !  GENERIC:: reshape_sph_field => reshape_sph_field_1d_ptr, &
  !                                 reshape_sph_field_2d_ptr
  !  !# GENERIC PROCEDURE, overloded to reallocate 1d and 2d arrays
  !  PROCEDURE:: reshape_sph_field_1d_ptr => reshape_sph_field_1d
  !  !! Reallocates a 1d array
  !  PROCEDURE:: reshape_sph_field_2d_ptr => reshape_sph_field_2d
  !  !! Reallocates a 2d array

    PROCEDURE:: allocate_particles_memory
    !! Allocates memory for the [[particles]] member arrays

    PROCEDURE:: deallocate_particles_memory
    !! Deallocates memory for the [[particles]] member arrays

    PROCEDURE:: compute_sph_hydro
    !# Computes the hydro fields on a section of the particles specified as
    !  input.
    !  First, computes the |sph| pressure starting from the |sph| baryon mass
    !  density, and the specific internal
    !  energy. The pressure is computed differently for different |eos|, and
    !  for cold and hot systems.
    !  Then computes the enthalpy and the sound speed accordingly.

    PROCEDURE:: read_compose_composition
    !! Reads the \(Y_e(n_b)\) table in the CompOSE file with extension .beta

    PROCEDURE:: compute_Ye
    !# Interpates linearly the electron fraction \(Y_e\) at the particle
    !  densities; that is, assigns \(Y_e\) at the particle positions

    PROCEDURE:: test_recovery
    !# Computes the conserved variables from the physical ones, and vice versa,
    !  to test that the recovered physical variables are the same to those
    !  computed from the |id|. @todo add reference for recovery

    !PROCEDURE:: compute_adm_momentum
    !# Computes an estimate of the \(\mathrm{ADM}\) linear momentum using
    !  the canonical momentum per baryon on the particles


    !
    !-- PUBLIC SUBROUTINES
    !

    PROCEDURE, PUBLIC:: analyze_hydro
    !# Scans the hydro fields taken from \(\texttt{|lorene|}\) to look
    !  for negative or zero values

    PROCEDURE, PUBLIC:: compute_and_print_sph_variables
    !# Computes the SPH variables at the particle positions, and optionally
    !  prints them to a binary file to be read by \(\texttt{SPHINCS_BSSN}\)
    !  and \(\texttt{splash}\), and to a formatted file to be read by
    !  \(\texttt{gnuplot}\), by calling
    !  [[particles:print_formatted_id_particles]]

    PROCEDURE, PUBLIC:: read_sphincs_dump_print_formatted
    !# Reads the binary ID file printed by
    !  [[particles:compute_and_print_sph_variables]], and prints it
    !  to a formatted file. @todo use this procedure in a second constructor

    PROCEDURE, PUBLIC:: print_formatted_id_particles
    !! Prints the |sph| |id| to a formatted file

    PROCEDURE, PUBLIC:: print_summary
    !! Prints the SPH ID to a formatted file

    PROCEDURE, PUBLIC:: is_empty
    !# Returns `.TRUE` if the [[particles]] object is empty, `.FALSE` otherwise
    !  @warning experimental, not actively used in the code yet

    !PROCEDURE, PUBLIC:: write_lorene_bns_id_dump

    !-----------------!
    !--  FUNCTIONS  --!
    !-----------------!

    PROCEDURE, PUBLIC:: get_npart
    !! Returns [[particles:npart]]
    PROCEDURE, PUBLIC:: get_npart_i
    !! Returns the number of particles on the object `i_matter`
    PROCEDURE, PUBLIC:: get_n_matter
    !! Returns [[particles:n_matter]]
    PROCEDURE, PUBLIC:: get_nuratio
    !! Returns [[particles:nuratio]]
    PROCEDURE, PUBLIC:: get_nuratio_i
    !! Returns the baryon number ratio on the object `i_matter`
    PROCEDURE, PUBLIC:: get_pos
    !! Returns [[particles:pos]]
    PROCEDURE, PUBLIC:: get_vel
    !! Returns [[particles:v]]
    PROCEDURE, PUBLIC:: get_nstar
    !! Returns [[particles:nstar]]
    PROCEDURE, PUBLIC:: get_nstar_sph
    !! Returns [[particles:nstar_sph]]
    PROCEDURE, PUBLIC:: get_nlrf
    !! Returns [[particles:nlrf]]
    PROCEDURE, PUBLIC:: get_nlrf_sph
    !! Returns [[particles:nlrf_sph]]
    PROCEDURE, PUBLIC:: get_nu
    !! Returns [[particles:nu]]
    PROCEDURE, PUBLIC:: get_u
    !! Returns [[particles:specific_energy]]
    PROCEDURE, PUBLIC:: get_u_sph
    !! Returns [[particles:u_sph]]
    PROCEDURE, PUBLIC:: get_pressure
    !! Returns [[particles:pressure]]
    PROCEDURE, PUBLIC:: get_pressure_sph
    !! Returns [[particles:pressure_sph]]
    PROCEDURE, PUBLIC:: get_theta
    !! Returns [[particles:theta]]
    PROCEDURE, PUBLIC:: get_h
    !! Returns [[particles:h]]
    PROCEDURE, PUBLIC:: get_lapse
    !! Returns [[particles:lapse]]
    PROCEDURE, PUBLIC:: get_shift
    !! Returns ([[particles:shift_x]],[[particles:shift_y]],[[particles:shift_z]])
    PROCEDURE, PUBLIC:: get_g3
    !! Returns ([[particles:g_xx]],[[particles:g_xy]],[[particles:g_xz]],[[particles:g_yy]],[[particles:g_yz]],[[particles:g_zz]])
    PROCEDURE, PUBLIC:: get_compose_eos
    !! Returns ([[particles:compose_eos]])

    FINAL:: destruct_particles
    !! Finalizer (Destructor) of [[particles]] object

  END TYPE particles

  !
  !-- Interface of the TYPE particles (i.e., declaration of the constructor)
  !-- Multiple procedures in this interface would overload the constructor.
  !-- Such procedures must have distingushable interfaces, in particular
  !-- distinguishable arguments)
  !
  INTERFACE particles
  !! Interface of TYPE [[particles]]

    MODULE PROCEDURE construct_particles_std
    !! Standard constructor for a [[particles]] object
    MODULE PROCEDURE construct_particles_bin
    !! Constructs a [[particles]] object from an |id| binary file

  END INTERFACE particles


  ABSTRACT INTERFACE

    SUBROUTINE post_process_sph_id_int &
      ( npart, pos, nlrf, u, pr, vel_u, theta, nstar, nu )
    !# Post-process the |sph| |id|; for example, correct for the residual
    !  ADM linear momentum.

    !IMPORT:: particles
    !CLASS(particles),                     INTENT(IN)   :: this
    INTEGER,                              INTENT(IN)   :: npart
    !! Particle number
    DOUBLE PRECISION, DIMENSION(3,npart), INTENT(INOUT):: pos
    !! Particle positions
    DOUBLE PRECISION, DIMENSION(npart),   INTENT(INOUT):: nlrf
    !! Baryon density in the local rest frame on the particles
    DOUBLE PRECISION, DIMENSION(npart),   INTENT(INOUT):: u
    !! Specific internal energy on the particles
    DOUBLE PRECISION, DIMENSION(npart),   INTENT(INOUT):: pr
    !! Pressure on the particles
    DOUBLE PRECISION, DIMENSION(3,npart), INTENT(INOUT):: vel_u
    !! Spatial velocity in the computing frame on the particles
    DOUBLE PRECISION, DIMENSION(npart),   INTENT(INOUT):: theta
    !! Generalized Lorentz factor on the particles
    DOUBLE PRECISION, DIMENSION(npart),   INTENT(INOUT):: nstar
    !! Proper baryon density in the local rest frame on the particles
    DOUBLE PRECISION, DIMENSION(npart),   INTENT(INOUT):: nu
    !! Baryon number per particle

    END SUBROUTINE post_process_sph_id_int

  END INTERFACE


  !
  !-- Interface of the constructor of TYPE particles
  !-- Its implementation is in submodule_particles_constructor.f90
  !
  INTERFACE


    MODULE FUNCTION construct_particles_std( id, dist ) RESULT ( parts )
    !! Standard constructor for a [[particles]] object

        CLASS(idbase), INTENT(INOUT):: id
        !# [[idbase]] object representing the BNS for which we want to place
        !  particles
        INTEGER,       INTENT(IN)   :: dist
        !# Identifier of the desired particle distribution:
        !
        !  - 0: Read particle positions (and optionally the baryon number per
        !       particle \(\nu\)) from a formatted file
        !
        !  - 1: Place particles on several lattices, each one surrounding a
        !       matter object
        !
        !  - 3: Place particles on ellipsoidal surfaces inside each matter object
        !
        TYPE(particles):: parts
        !! Constructed [[particles]] object

    END FUNCTION construct_particles_std


    MODULE FUNCTION construct_particles_bin( id, namefile ) RESULT ( parts )
    !! Constructs a [[particles]] object from an |id| binary file

        CLASS(idbase), INTENT(INOUT):: id
        !# [[idbase]] object representing the BNS for which we want to place
        !  particles
        CHARACTER(LEN=*), INTENT(INOUT) :: namefile
        !! Name of the |id| binary file
        TYPE(particles):: parts
        !! Constructed [[particles]] object

    END FUNCTION construct_particles_bin


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


    MODULE SUBROUTINE read_particles_formatted_file &
      ( this, parts_pos_unit, nline_in, nline_fin, &
        xmin, xmax, ymin, ymax, zmin, zmax, &
        pos, pvol, nu, h )
    !# Read particle positions and nu from a formatted file with the following
    !  format:
    !
    !   1. The first line must contain the integer number of objects for the
    !      system, and the particle numbers on each object;
    !      each column separated by one tab.
    !   2. The other lines must contain at least 3 columns with the x, y, z
    !      coordinates of the particle positions. An optional 4th column can
    !      contain the values of nu on the particles. If the 4th column is
    !      not present, nu will be roughly estimated using the density read
    !      from the ID and the average volume per particle; the latter is
    !      computed by computing the average particle separation along the
    !      z direction.


      CLASS(particles), INTENT(INOUT):: this
      !! [[particles]] object which this PROCEDURE is a member of
      !CHARACTER(LEN=*), INTENT(IN)   :: parts_pos_namefile
      !# String containing the name of the formatted file containing the
      !  particle positions and, optionally, nu
      INTEGER:: parts_pos_unit
      !# Unit number of the formatted file containing the particle positions
      !  and, optionally, nu
      INTEGER:: nline_in
      !! First line containing the relevant data
      INTEGER:: nline_fin
      !! Last line containing the relevant data
      DOUBLE PRECISION, INTENT(IN)    :: xmin
      !! Left \(x\) boundary of the lattice
      DOUBLE PRECISION, INTENT(IN)    :: xmax
      !! Right \(x\) boundary of the lattice
      DOUBLE PRECISION, INTENT(IN)    :: ymin
      !! Left \(y\) boundary of the lattice
      DOUBLE PRECISION, INTENT(IN)    :: ymax
      !! Right \(y\) boundary of the lattice
      DOUBLE PRECISION, INTENT(IN)    :: zmin
      !! Left \(z\) boundary of the lattice
      DOUBLE PRECISION, INTENT(IN)    :: zmax
      !! Right \(z\) boundary of the lattice
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT):: pos
      !> Array storing the particle positions
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, INTENT(INOUT):: pvol
      !! Array storing the particle volumes
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, INTENT(INOUT):: nu
      !! Array storing the particle baryon masses
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, INTENT(INOUT):: h
      !! Array storing the initial guess for the particle smoothing lengths

    END SUBROUTINE read_particles_formatted_file


    MODULE SUBROUTINE place_particles_lattice( this, &
                                  central_density, &
                                  xmin, xmax, ymin, ymax, zmin, zmax, &
                                  npart_des, npart_out, stretch, &
                                  thres, pos, pvol, nu, h, &
                                  get_density, get_id, validate_position )
    !! Places particles on a lattice containing a physical object

      CLASS(particles), INTENT(INOUT):: this
      !! [[particles]] object which this PROCEDURE is a member of
      DOUBLE PRECISION, INTENT(IN)    :: central_density
      !! Maximum baryon mass density of the system
      DOUBLE PRECISION, INTENT(IN)    :: xmin
      !! Left \(x\) boundary of the lattice
      DOUBLE PRECISION, INTENT(IN)    :: xmax
      !! Right \(x\) boundary of the lattice
      DOUBLE PRECISION, INTENT(IN)    :: ymin
      !! Left \(y\) boundary of the lattice
      DOUBLE PRECISION, INTENT(IN)    :: ymax
      !! Right \(y\) boundary of the lattice
      DOUBLE PRECISION, INTENT(IN)    :: zmin
      !! Left \(z\) boundary of the lattice
      DOUBLE PRECISION, INTENT(IN)    :: zmax
      !! Right \(z\) boundary of the lattice
      INTEGER,          INTENT(IN)    :: npart_des
      !! Desired particle number
      INTEGER,          INTENT(OUT)   :: npart_out
      !! Real, output particle number
      DOUBLE PRECISION, INTENT(IN)    :: stretch
      !! Stretching factor fo the lattice. `xmin` to `zmax` are multiplied by it
      DOUBLE PRECISION, INTENT(IN)    :: thres
      !# (~rho_max)/thres is the minimum mass density considered
      ! when placing particles. Used only when redistribute_nu is
      ! .FALSE. . When redistribute_nu is .TRUE. thres= 100*nu_ratio
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT):: pos
      !> Array storing the particle positions
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, INTENT(INOUT):: pvol
      !! Array storing the particle volumes
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, INTENT(INOUT):: nu
      !! Array storing the particle baryon masses
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, INTENT(INOUT):: h
      !! Array storing the initial guess for the particle smoothing lengths
      INTERFACE
        FUNCTION get_density( x, y, z ) RESULT( density )
          !! Returns the baryon mass density at the desired point
          DOUBLE PRECISION, INTENT(IN):: x
          !! \(x\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: y
          !! \(y\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: z
          !! \(z\) coordinate of the desired point
          DOUBLE PRECISION:: density
          !! Baryon mass density at \((x,y,z)\)
        END FUNCTION get_density
      END INTERFACE
      INTERFACE
        SUBROUTINE get_id( x, y, z , sqdetg, baryon_density, gamma_euler )
          !! Returns the baryon mass density at the desired point
          DOUBLE PRECISION, INTENT(IN):: x
          !! \(x\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: y
          !! \(y\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: z
          !! \(z\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(OUT):: sqdetg
          DOUBLE PRECISION, INTENT(OUT):: baryon_density
          DOUBLE PRECISION, INTENT(OUT):: gamma_euler
        END SUBROUTINE get_id
      END INTERFACE
      INTERFACE
        FUNCTION validate_position_int( x, y, z ) RESULT( answer )
        !! Returns 1 if the position is not valid, 0 otherwise
          DOUBLE PRECISION, INTENT(IN):: x
          !! \(x\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: y
          !! \(y\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: z
          !! \(z\) coordinate of the desired point
          LOGICAL:: answer
          !! `.TRUE.` if the position is valid, `.FALSE.` otherwise
        END FUNCTION validate_position_int
      END INTERFACE
      !> Returns 1 if the position is not valid, 0 otherwise
      PROCEDURE(validate_position_int), OPTIONAL:: validate_position

    END SUBROUTINE place_particles_lattice


    MODULE SUBROUTINE place_particles_ellipsoidal_surfaces( this, &
                                  mass_star, radius, center, &
                                  central_density, npart_des, &
                                  npart_out, pos, pvol, nu, h, &
                                  last_r, upper_bound, lower_bound, &
                                  upper_factor, lower_factor, max_steps, &
                                  filename_mass_profile, &
                                  filename_shells_radii, &
                                  filename_shells_pos, &
                                  get_density, integrate_density, &
                                  get_id, validate_position, &
                                  radii, surf )
    !! Places particles on ellipsoidal surfaces on one star

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(INOUT):: this
      !& [[idbase]] object needed to access the BNS data
      !  @TODO Remove the [[idbase]] argument as done in SUBROUTINE perform_apm
      !CLASS(idbase),       INTENT(INOUT):: id
      !> Approximate particle number on the star
      INTEGER,          INTENT(IN)    :: npart_des
      !> Final number of particles on the star
      INTEGER,          INTENT(OUT)   :: npart_out
      !& If, after max_steps, the iteration did not converge,
      !  multiply upper_bound by upper_factor, and lower_bound
      !  by lower_factor. max_steps >= 10. 100 is a nice value
      INTEGER,          INTENT(IN)    :: max_steps
      !> Baryonic mass of the star
      DOUBLE PRECISION, INTENT(IN)    :: mass_star
      !> Radius of the star in the x direction towards the companion
      DOUBLE PRECISION, INTENT(IN)    :: radius
      !& \(x|) coordinate of the center of the star, i.e.,
      !  of the point with highest density
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: center
      !> Central density of the star, i.e., highest density
      DOUBLE PRECISION, INTENT(IN)    :: central_density
      !> Radius of the last ellipsoidal surface
      DOUBLE PRECISION, INTENT(IN)    :: last_r
      !& If, after max_steps, the iteration did not converge,
      !  multiply upper_bound by upper_factor, and lower_bound
      !  by lower_factor. upper_factor >= 1, usually an increase of 1% works
      DOUBLE PRECISION, INTENT(IN)    :: upper_factor
      !& If, after max_steps, the iteration did not converge,
      !  multiply upper_bound by upper_factor, and lower_bound
      !  by lower_factor. lower_factor <= 1, usually a decrease of 1% works
      DOUBLE PRECISION, INTENT(IN)    :: lower_factor
      !& Desired upper bound for the differences between particle
      !  masses on neighbouring ellipsoidal surfaces
      DOUBLE PRECISION, INTENT(INOUT) :: upper_bound
      !& Desired lower bound for the differences between particle
      !  masses on neighbouring ellipsoidal surfaces
      DOUBLE PRECISION, INTENT(INOUT) :: lower_bound
      !> Array string the final positions
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT):: pos
      !> Array soring the inal particle volumes
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, INTENT(INOUT):: pvol
      !> Array storing the final particle masses
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, INTENT(INOUT):: nu
      !! Array storing the particle baryon masses
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, INTENT(INOUT):: h
      !! Array storing the initial guess for the particle smoothing lengths
      CHARACTER(LEN=*), INTENT(INOUT), OPTIONAL :: filename_mass_profile
      !& Name of the file to store the surface radii
      !  @TODO change name of variable to filename_surfaces_radii
      CHARACTER(LEN=*), INTENT(INOUT), OPTIONAL :: filename_shells_radii
      !& Name of the file to store the final particle positions
      !  @TODO change name of variable to filename_surfaces_pos
      CHARACTER(LEN=*), INTENT(INOUT), OPTIONAL :: filename_shells_pos
      INTERFACE
        FUNCTION get_density( x, y, z ) RESULT( density )
          !! Returns the baryon mass density at the desired point
          DOUBLE PRECISION, INTENT(IN):: x
          !! \(x\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: y
          !! \(y\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: z
          !! \(z\) coordinate of the desired point
          DOUBLE PRECISION:: density
          !! Baryon mass density at \((x,y,z)\)
        END FUNCTION get_density
      END INTERFACE
      INTERFACE
        SUBROUTINE get_id( x, y, z , sqdetg, baryon_density, gamma_euler )
          !! Returns the baryon mass density at the desired point
          DOUBLE PRECISION, INTENT(IN):: x
          !! \(x\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: y
          !! \(y\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: z
          !! \(z\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(OUT):: sqdetg
          DOUBLE PRECISION, INTENT(OUT):: baryon_density
          DOUBLE PRECISION, INTENT(OUT):: gamma_euler
        END SUBROUTINE get_id
      END INTERFACE
      INTERFACE
        SUBROUTINE integrate_density &
          ( center, radius, central_density, dr, dth, dphi, &
            mass, mass_profile, mass_profile_idx, radii, surf )
          IMPORT:: surface
          !> Center of the star
          DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: center
          !> Central density of the star
          DOUBLE PRECISION, INTENT(IN)   :: central_density
          !> Radius of the star
          DOUBLE PRECISION, INTENT(IN)   :: radius
          !> Integration steps
          DOUBLE PRECISION, INTENT(IN)   :: dr, dth, dphi
          !> Integrated mass of the star
          DOUBLE PRECISION, INTENT(INOUT):: mass
          !> Array storing the radial mass profile of the star
          !DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT):: &
          !                                 mass_profile
          DOUBLE PRECISION, DIMENSION(3,0:NINT(radius/dr)), INTENT(OUT):: &
                                               mass_profile
          !& Array to store the indices for array mass_profile, sorted so that
          !  mass_profile[mass_profile_idx] is in increasing order
          !INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT)::mass_profile_idx
          INTEGER, DIMENSION(0:NINT(radius/dr)), INTENT(OUT):: &
                                               mass_profile_idx

          DOUBLE PRECISION, DIMENSION(2), INTENT(IN), OPTIONAL:: radii
          !> Surface of the matter object
          TYPE(surface),                  INTENT(IN), OPTIONAL:: surf
        END SUBROUTINE integrate_density
      END INTERFACE
      INTERFACE
        FUNCTION validate_position_int( x, y, z ) RESULT( answer )
        !! Returns 1 if the position is not valid, 0 otherwise
          DOUBLE PRECISION, INTENT(IN):: x
          !! \(x\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: y
          !! \(y\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: z
          !! \(z\) coordinate of the desired point
          LOGICAL:: answer
          !! `.TRUE.` if the position is valid, `.FALSE.` otherwise
        END FUNCTION validate_position_int
      END INTERFACE
      !> Returns 1 if the position is not valid, 0 otherwise
      PROCEDURE(validate_position_int), OPTIONAL:: validate_position
      !DOUBLE PRECISION, INTENT(IN),   OPTIONAL:: pmass_des
      DOUBLE PRECISION, DIMENSION(2), INTENT(IN), OPTIONAL:: radii
      !> Surface of the matter object
      TYPE(surface),                  INTENT(IN), OPTIONAL :: surf

    END SUBROUTINE place_particles_ellipsoidal_surfaces


    MODULE SUBROUTINE impose_equatorial_plane_symmetry( npart, pos, &
                                                 nu, com_star, verbose )
    !# Mirror the particle with z>0 with respect to the xy plane,
    !  to impose the equatorial-plane symmetry

      !CLASS(particles), INTENT(INOUT):: this
      INTEGER, INTENT(INOUT):: npart
      DOUBLE PRECISION, INTENT(IN), OPTIONAL:: com_star
      LOGICAL, INTENT(IN), OPTIONAL:: verbose

      DOUBLE PRECISION, DIMENSION(3,npart), INTENT(INOUT):: pos
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(INOUT), OPTIONAL:: nu

    END SUBROUTINE impose_equatorial_plane_symmetry


    MODULE SUBROUTINE find_particles_above_xy_plane( npart, pos, &
                                            npart_above_xy, above_xy_plane_a )
    !# Mirror the particle with z>0 with respect to the xy plane,
    !  to impose the equatorial-plane symmetry

      INTEGER, INTENT(IN) :: npart
      INTEGER, INTENT(OUT):: npart_above_xy

      DOUBLE PRECISION, DIMENSION(3,npart), INTENT(IN):: pos

      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT):: above_xy_plane_a

    END SUBROUTINE find_particles_above_xy_plane


    MODULE SUBROUTINE reflect_particles_xy_plane( npart, pos, pos_below,  &
                                        npart_above_xy, above_xy_plane_a, &
                                        nu, nu_below )
    !# Mirror the particle with z>0 with respect to the xy plane,
    !  to impose the equatorial-plane symmetry

      INTEGER, INTENT(IN):: npart
      INTEGER, INTENT(INOUT):: npart_above_xy
      INTEGER, DIMENSION(npart_above_xy), INTENT(IN):: above_xy_plane_a

      DOUBLE PRECISION, DIMENSION(3,npart),          INTENT(INOUT) :: pos
      DOUBLE PRECISION, DIMENSION(3,npart_above_xy), INTENT(OUT)   :: pos_below
      DOUBLE PRECISION, DIMENSION(npart),        INTENT(INOUT), OPTIONAL:: nu
      DOUBLE PRECISION, DIMENSION(npart_above_xy), INTENT(OUT), OPTIONAL:: &
      nu_below

    END SUBROUTINE reflect_particles_xy_plane


    MODULE PURE FUNCTION get_object_of_particle( this, a )
    !# Returns the number of the matter object asociated with the particle
    !  number given as input. Example: give number \(n\) as input; this
    !  particle number corresponds to a particle on matter object \(m\).
    !  This functions returns \(m\).

      CLASS(particles), INTENT(IN):: this
      !! [[particles]] object which this PROCEDURE is a member of
      INTEGER,          INTENT(IN):: a
      !! Particle number
      INTEGER:: get_object_of_particle

    END FUNCTION get_object_of_particle


!    MODULE SUBROUTINE reshape_sph_field_1d( this, field, new_size1, new_size2, &
!                                            index_array )
!    !! Reallocates a 1d array
!
!      !> [[particles]] object which this PROCEDURE is a member of
!      CLASS(particles), INTENT(INOUT):: this
!      !> New particle number on star 1
!      INTEGER,                        INTENT(IN):: new_size1
!      !> New particle number on star 2
!      INTEGER,                        INTENT(IN):: new_size2
!      !> Array to select elements to keep in the reshaped array
!      INTEGER,          DIMENSION(:), INTENT(IN):: index_array
!      !> 1D array to reshape
!      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(INOUT):: field
!
!    END SUBROUTINE reshape_sph_field_1d
!
!
!    MODULE SUBROUTINE reshape_sph_field_2d( this, field, new_size1, new_size2, &
!                                            index_array )
!    !! Reallocates a 2d array
!
!      !> [[particles]] object which this PROCEDURE is a member of
!      CLASS(particles), INTENT(INOUT):: this
!      !> New particle number on star 1
!      INTEGER,                        INTENT(IN):: new_size1
!      !> New particle number on star 2
!      INTEGER,                        INTENT(IN):: new_size2
!      !> Array to select elements to keep in the reshaped array
!      INTEGER,          DIMENSION(:), INTENT(IN):: index_array
!      !> 2D array to reshape
!      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT):: field
!
!    END SUBROUTINE reshape_sph_field_2d


    MODULE SUBROUTINE allocate_particles_memory( this )
    !! Allocates allocatable arrays member of a [[particles]] object

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(INOUT):: this

    END SUBROUTINE allocate_particles_memory


    MODULE SUBROUTINE deallocate_particles_memory( this )
    !! Deallocates allocatable arrays member of a [[particles]] object

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(INOUT):: this

    END SUBROUTINE deallocate_particles_memory


  END INTERFACE

  !
  !-- Interfaces of the methods of TYPE particles
  !-- Their implementations are in module_particles_methods.f90
  !
  INTERFACE

    MODULE SUBROUTINE analyze_hydro( this, namefile )
    !# Scans the hydro fields taken from \(\texttt{|lorene|}\) to look
    !  for negative or zero values

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles),    INTENT(INOUT):: this
      !& Name of the formatted file where the particle positions at which
      !  some of the hydro fields are negative or zero are printed to
      CHARACTER(LEN=*), INTENT(INOUT), OPTIONAL :: namefile

    END SUBROUTINE analyze_hydro

    MODULE SUBROUTINE compute_and_print_sph_variables( this, namefile )
    !# Computes the SPH variables at the particle positions, and optionally
    !  prints them to a binary file to be read by \(\texttt{SPHINCS_BSSN}\)
    !  and \(\texttt{splash}\), and to a formatted file to be read by
    !  \(\texttt{gnuplot}\), by calling
    !  [[particles:print_formatted_id_particles]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(INOUT):: this
      !> Name of the formatted file where the SPH ID is printed to
      CHARACTER(LEN=*), INTENT(INOUT), OPTIONAL :: namefile

    END SUBROUTINE compute_and_print_sph_variables

    MODULE SUBROUTINE compute_sph_hydro( this, npart_in, npart_fin, &
      eqos, nlrf, u, pr, enthalpy, cs, verbose )
    !# Computes the hydro fields on a section of the particles specified as
    !  input.
    !  First, computes the |sph| pressure starting from the |sph| baryon mass
    !  density, and the specific internal
    !  energy. The pressure is computed differently for different |eos|, and
    !  for cold and hot systems.
    !  Then computes the enthalpy and the sound speed accordingly.

      CLASS(particles), INTENT(INOUT):: this
      !! [[particles]] object which this PROCEDURE is a member of
      INTEGER,          INTENT(IN):: npart_in
      !! First index of the desired section of the particles
      INTEGER,          INTENT(IN):: npart_fin
      !! Last index of the desired section of the particles
      CLASS(eos),       INTENT(IN):: eqos
      !! |eos| to be used
      DOUBLE PRECISION, DIMENSION(npart_fin - npart_in + 1), INTENT(IN)   :: &
      nlrf
      !! Baryon mass density in the local rest frame
      DOUBLE PRECISION, DIMENSION(npart_fin - npart_in + 1), INTENT(INOUT):: u
      !! Specific internal energy
      DOUBLE PRECISION, DIMENSION(npart_fin - npart_in + 1), INTENT(INOUT):: Pr
      !! Pressure
      DOUBLE PRECISION, DIMENSION(npart_fin - npart_in + 1), INTENT(INOUT):: &
      enthalpy
      !! Enthalpy
      DOUBLE PRECISION, DIMENSION(npart_fin - npart_in + 1), INTENT(INOUT):: cs
      !! Speed of sound
      LOGICAL, INTENT(IN), OPTIONAL:: verbose

    END SUBROUTINE compute_sph_hydro

    MODULE SUBROUTINE perform_apm( &
      ! PROCEDURES to get the density and pressure at a point
      get_density, get_nstar_id, &
      get_pressure_id, compute_pressure, &
      ! Arguments pertaining to the matter object
      npart_output, pos_input, pvol, h_output, nu_output, &
      center, com_star, mass, sizes, eqos, &
      ! Steering parameters for the APM iteration
      apm_max_it, max_inc, mass_it, correct_nu, &
      nuratio_thres, nuratio_des, use_pressure, &
      ! Arguments pertaining to the ghost particles
      adapt_ghosts, move_away_ghosts, &
      nx_gh, ny_gh, nz_gh, ghost_dist, &
      ! Arguments pertaining to the atmosphere
      use_atmosphere, remove_atmosphere, &
      ! Arguments pertaining to input/output
      print_step, namefile_pos_id, namefile_pos, namefile_results, &
      ! Optional argument
      validate_position, surf )
    !! Performs the Artificial Pressure Method (APM)

      !CLASS(particles),     INTENT(INOUT):: this
      !! [[particles]] object which this PROCEDURE is a member of
      INTERFACE
        FUNCTION get_density( x, y, z ) RESULT( density )
          !! Returns the baryon mass density at the desired point
          DOUBLE PRECISION, INTENT(IN):: x
          !! \(x\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: y
          !! \(y\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: z
          !! \(z\) coordinate of the desired point
          DOUBLE PRECISION:: density
          !! Baryon mass density at \((x,y,z)\)
        END FUNCTION get_density
      END INTERFACE
      INTERFACE
        SUBROUTINE get_nstar_id( npart, x, y, z, nstar_sph, nstar_id, &
                                 nlrf_sph, sqg )
        !! Computes the proper baryon number density at the particle positions
          INTEGER, INTENT(IN):: npart
          !! Number of real particles (i.e., no ghost particles included here)
          DOUBLE PRECISION, INTENT(IN):: x(npart)
          !! Array of \(x\) coordinates
          DOUBLE PRECISION, INTENT(IN):: y(npart)
          !! Array of \(y\) coordinates
          DOUBLE PRECISION, INTENT(IN):: z(npart)
          !! Array of \(z\) coordinates
          DOUBLE PRECISION, INTENT(IN):: nstar_sph(npart)
          !! |sph| proper baryon density
          DOUBLE PRECISION, INTENT(OUT):: nstar_id(npart)
          !! Array to store the computed proper baryon number density
          DOUBLE PRECISION, INTENT(OUT):: nlrf_sph(npart)
          !# Array to store the local rest frame baryon density computed from
          !  the |sph| proper baryon density
          DOUBLE PRECISION, INTENT(OUT):: sqg(npart)
          !# Square root of minus the determinant of the spacetime metric
        END SUBROUTINE get_nstar_id
      END INTERFACE
      INTERFACE
        SUBROUTINE compute_pressure &
          ( npart, x, y, z, nlrf, eqos, pressure, verbose )
          IMPORT:: eos
          INTEGER,          INTENT(IN):: npart
          !! Returns the baryon mass density at the desired point
          DOUBLE PRECISION, INTENT(IN):: x(npart)
          !! \(x\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: y(npart)
          !! \(y\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: z(npart)
          !! \(z\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN)   :: nlrf(npart)
          !! Baryon mass density in the local rest frame
          TYPE(eos),        INTENT(IN)   :: eqos
          !! |eos| to use
          DOUBLE PRECISION, INTENT(INOUT):: pressure(npart)
          !! Baryon mass density at \((x,y,z)\)
          LOGICAL, INTENT(IN), OPTIONAL:: verbose
          !# If .TRUE., print informative standard output about how the
          !  pressure is computed. Default is .TRUE.
        END SUBROUTINE compute_pressure
      END INTERFACE
      INTERFACE
        FUNCTION get_pressure_id( x, y, z ) RESULT( pressure )
          !! Returns the baryon mass density at the desired point
          DOUBLE PRECISION, INTENT(IN):: x
          !! \(x\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: y
          !! \(y\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: z
          !! \(z\) coordinate of the desired point
          DOUBLE PRECISION:: pressure
          !! Baryon mass density at \((x,y,z)\)
        END FUNCTION get_pressure_id
      END INTERFACE
      INTERFACE
        FUNCTION validate_position_int( x, y, z ) RESULT( answer )
        !! Returns 1 if the position is not valid, 0 otherwise
          DOUBLE PRECISION, INTENT(IN):: x
          !! \(x\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: y
          !! \(y\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: z
          !! \(z\) coordinate of the desired point
          LOGICAL:: answer
          !! `.TRUE.` if the position is valid, `.FALSE.` otherwise
        END FUNCTION validate_position_int
      END INTERFACE
      !> Returns 1 if the position is not valid, 0 otherwise
      PROCEDURE(validate_position_int), OPTIONAL:: validate_position
      !> Initial particle number
      INTEGER,                          INTENT(INOUT):: npart_output
      !& Prints the particle positions to a formatted file every
      !  print_step steps
      INTEGER,                          INTENT(INOUT):: print_step
      !> Initial particle positions
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT):: pos_input
      !> Initial particle volume
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, INTENT(INOUT):: pvol
      !& Array to store the smoothing lengths computed at the end of the
      !  APM iteration
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, INTENT(INOUT):: h_output
      !& Array to store the baryon number per particle computed at the end of
      !  the APM iteration
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, INTENT(INOUT):: nu_output
      !> Center of the matter object, from the |id|
      DOUBLE PRECISION, DIMENSION(3),   INTENT(IN)   :: center
      !> Center of mass of the matter object, from the |id|
      DOUBLE PRECISION, DIMENSION(3),   INTENT(INOUT):: com_star
      !> Mass of the matter object
      DOUBLE PRECISION,                 INTENT(IN)   :: mass
      !> Sizes of the matter object
      DOUBLE PRECISION, DIMENSION(6),   INTENT(IN)   :: sizes
      !> |eos| to use when computing the pressure
      TYPE(eos),                        INTENT(IN)   :: eqos
      !> Maximum number of APM iterations, irrespective of the EXIT condition
      INTEGER,                          INTENT(IN)   :: apm_max_it
      !& Sets the EXIT condition: If the average over all the
      !  particles of the relative error in the density estimate
      !  grows max_inc times, exit the iteration.
      INTEGER,                          INTENT(IN)   :: max_inc
      !& If .TRUE. performs a second iteration after the APM one, without moving
      !  the particles, changing their mass in order to better match
      !  the star density. The mass ratio grows very fast in all the tried
      !  experiments, hence the suggested value is .FALSE.
      LOGICAL,                          INTENT(IN)   :: mass_it
      !& If .TRUE., the baryon number per particle nu is corrected
      !  to include the total baryonic masses of the
      !  stars.
      LOGICAL,                          INTENT(IN)   :: correct_nu
      !& Maximum mass ratio (equivalently baryon number ratio)
      !  to be used in the one-time-only final correction
      !  of the particle masses to match the star density even
      !  better (without moving the particles)
      DOUBLE PRECISION,                 INTENT(IN)   :: nuratio_thres
      !& Sets the EXIT condition: If the baryon number ratio
      !  is within 2.5% of nuratio_des, exit the iteration
      !  Set nuratio_des to 0 to deactivate and exit the APM
      !  iteration using max_inc
      DOUBLE PRECISION,                 INTENT(IN)   :: nuratio_des
      !& If .TRUE., uses the physical pressure computed with
      !  the |eos| using the SPH estimate of the density [[particles:nlrf_sph]],
      !  to compute the artificial pressure. Otherwise, the
      !  density variable [[particles:nstar_sph]] is used to
      !  compute the artificial pressure
      LOGICAL,                          INTENT(IN)   :: use_pressure
      !& If .TRUE., the ghost particles will be placed and have
      !  a baryon number such to reproduce the density of the
      !  outermost layers (r > 99% of the minimum radius) of
      !  the object. If .TRUE., the arguments nx_gh, ny_gh,
      !  nz_gh, ghost_dist are ignored; if .FALSE., they are
      !  instead used to place the ghost particles
      LOGICAL,                          INTENT(IN)   :: adapt_ghosts
      !& If .TRUE., the ghost particles will slowly move away
      !  from the surface of the star during the iteration
      !  (depending on the EXIT condition chosen, this happens
      !  in slightly different ways), to allow for the real
      !  particles to get closer to the surface
      LOGICAL,                          INTENT(IN)   :: move_away_ghosts
      !> Number of lattice points in the x direction for ghosts
      INTEGER,                          INTENT(IN)   :: nx_gh
      !> Number of lattice points in the y direction for ghosts
      INTEGER,                          INTENT(IN)   :: ny_gh
      !> Number of lattice points in the z direction for ghosts
      INTEGER,                          INTENT(IN)   :: nz_gh
      !& Distance between the ghost particles and the surface
      !  of the matter object considered (star, ejecta, etc...)
      DOUBLE PRECISION,                 INTENT(INOUT):: ghost_dist
      !& If .TRUE., allows particles to move where the density
      !  is 0, and displace them using only the smoothing term.
      !  This can be useful when the system has an irregular
      !  geometry, as, for example, an ejecta; `.FALSE.` otherwise
      LOGICAL,                          INTENT(INOUT):: use_atmosphere
      !& If .TRUE., removes the particles placed where the density is 0,
      !  at the end of the APM iteration; `.FALSE.` otherwise
      LOGICAL,                          INTENT(INOUT):: remove_atmosphere
      !& Name for the formatted file where the initial particle positions
      !  and the ghost positions will be printed
      CHARACTER(LEN=*),              INTENT(INOUT), OPTIONAL :: &
                                                            namefile_pos_id
      !& Name for the formatted file where the particle positions
      !  and the ghost positions will be printed every 15 iterations
      CHARACTER(LEN=*),              INTENT(INOUT), OPTIONAL :: &
                                                            namefile_pos
      !& Name for the formatted file where various quantities related
      !  to the particle distribution, the baryon number particle and the
      !  kernel estimate of the density will be printed at the end of
      !  the APM iteration
      CHARACTER(LEN=*),              INTENT(INOUT), OPTIONAL :: &
                                                            namefile_results
      !> Surface of the matter object
      TYPE(surface),                 INTENT(IN)   , OPTIONAL :: surf

    END SUBROUTINE perform_apm


    MODULE SUBROUTINE test_recovery( this, npart, pos, nlrf, u, pr, vel_u, &
                                     theta, nstar )!, namefile )
    !# Tests the recovery. Computes the conserved variables from the physical
    !  ones, and then the physical ones from the conserved ones. It then
    !  compares the variables computed with the recovery PROCEDURES, with
    !  those computed with |sphincsid|. @todo add reference for recovery

      CLASS(particles),                     INTENT(INOUT):: this
      !! [[particles]] object which this PROCEDURE is a member of
      INTEGER,                              INTENT(IN)   :: npart
      !! Particle number
      DOUBLE PRECISION, DIMENSION(3,npart), INTENT(IN)   :: pos
      !! Particle positions
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN)   :: nlrf
      !! Baryon density in the local rest frame on the particles
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN)   :: u
      !! Specific internal energy on the particles
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN)   :: pr
      !! Pressure on the particles
      DOUBLE PRECISION, DIMENSION(3,npart), INTENT(IN)   :: vel_u
      !! Spatial velocity in the computing frame on the particles
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN)   :: theta
      !! Generalized Lorentz factor on the particles
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN)   :: nstar
      !! Proper baryon density in the local rest frame on the particles
      !DOUBLE PRECISION, DIMENSION(3,npart), INTENT(IN)   :: s_l
      !! Canonical momentum on the particles
      !DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN)   :: e_hat
      !! Canonical energy on the particles
      !CHARACTER(LEN=*),                  INTENT(INOUT), OPTIONAL :: namefile
      !! Name of the formatted file where the data is printed

    END SUBROUTINE test_recovery


    MODULE SUBROUTINE read_sphincs_dump_print_formatted( this, namefile_bin, &
                                                         namefile, save_data )
    !# Reads the binary ID file printed by
    !  [[particles:compute_and_print_sph_variables]]
    !   and prints the data stored in it to a formatted file

      CLASS(particles),    INTENT(INOUT)           :: this
      !! [[particles]] object which this PROCEDURE is a member of
      CHARACTER(LEN=*), INTENT(INOUT), OPTIONAL :: namefile_bin
      !! Name of the binary file to be read
      CHARACTER(LEN=*), INTENT(INOUT), OPTIONAL :: namefile
      !! Name of the formatted file to be printed
      LOGICAL,             INTENT(IN),    OPTIONAL :: save_data
      !! If `.TRUE.`, saves the read data into the TYPE member variables

    END SUBROUTINE read_sphincs_dump_print_formatted


    MODULE SUBROUTINE print_formatted_id_particles( this, namefile )
    !! Prints the SPH ID to a formatted file

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles),    INTENT(INOUT)           :: this
      !> Name of the formatted output file
      CHARACTER(LEN=*), INTENT(INOUT), OPTIONAL :: namefile

    END SUBROUTINE print_formatted_id_particles


    MODULE SUBROUTINE read_compose_composition( this, namefile )
    !! Reads the \(Y_e(n_b)\) table in the CompOSE file with extension .beta

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles),    INTENT(INOUT)           :: this
      !& To read the file great_eos.beta in directory compose_path/GREAT_EoS,
      !  namefile="GREAT_EoS/great_eos"
      CHARACTER(LEN=*), INTENT(INOUT), OPTIONAL :: namefile

    END SUBROUTINE read_compose_composition


    MODULE SUBROUTINE compute_Ye( this )!, nlrf, Ye )
    !# Interpolates linearly the electron fraction \(Y_e\) at the particle
    !  densities; that is, assigns \(Y_e\) at the particle positions

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles),    INTENT(INOUT)           :: this
      !DOUBLE PRECISION, DIMENSION( : ), INTENT(IN):: nlrf
      !DOUBLE PRECISION, DIMENSION( : ), INTENT(OUT):: Ye

    END SUBROUTINE compute_Ye


    MODULE SUBROUTINE print_summary( this )!, filename )
    !# Prints a summary of the properties of the |sph| particle
    !  distribution, optionally, to a formatted file whose name
    !  is given as the optional argument `filename`


      CLASS(particles), INTENT(INOUT):: this
      !CHARACTER(LEN=*), INTENT(INOUT), OPTIONAL:: filename
      !! Name of the formatted file to print the summary to

    END SUBROUTINE print_summary


    MODULE SUBROUTINE destruct_particles( this )
    !> Finalizer (Destructor) of [[particles]] object

      !> [[particles]] object which this PROCEDURE is a member of
      TYPE(particles), INTENT(INOUT):: this

    END SUBROUTINE destruct_particles


    !-----------------!
    !--  FUNCTIONS  --!
    !-----------------!


    MODULE FUNCTION is_empty( this ) RESULT( answer )
    !# Returns `.TRUE` if the [[particles]] object is empty, `.FALSE` otherwise
    !  @warning experimental, not actively used in the code yet

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(IN):: this
      !> `.TRUE` if the [[particles]] object is empty, `.FALSE` otherwise
      LOGICAL:: answer

    END FUNCTION is_empty


   !MODULE SUBROUTINE write_lorene_bns_id_dump( this, namefile )
   !
   !    CLASS(particles),    INTENT(IN)               :: this
   !    CHARACTER(LEN=*), INTENT(INOUT), OPTIONAL :: namefile
   !
   !END SUBROUTINE write_lorene_bns_id_dump


    MODULE PURE FUNCTION get_n_matter( this ) RESULT( n_matter )
    !! Returns [[particles:n_matter]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(IN):: this
      !> [[particles:n_matter]]
      INTEGER:: n_matter

    END FUNCTION get_n_matter


    MODULE PURE FUNCTION get_npart( this ) RESULT( n_part )
    !! Returns [[particles:npart]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(IN):: this
      !> [[particles:npart]]
      INTEGER:: n_part

    END FUNCTION get_npart


    MODULE PURE FUNCTION get_npart_i( this, i_matter ) RESULT( n_part )
    !! Returns the number of particles on the object `i_matter`

      CLASS(particles), INTENT(IN):: this
      !! [[particles]] object which this PROCEDURE is a member of
      INTEGER, INTENT(IN):: i_matter
      !! Index of the matter object
      INTEGER:: n_part
      !! Number of particles on the object `i_matter`

    END FUNCTION get_npart_i


    MODULE PURE FUNCTION get_nuratio( this ) RESULT( nuratio )
    !! Returns [[particles:nuratio]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(IN):: this
      !> [[particles:nuratio]]
      DOUBLE PRECISION:: nuratio

    END FUNCTION get_nuratio


    MODULE PURE FUNCTION get_nuratio_i( this, i_matter ) RESULT( nuratio )
    !! Returns the baryon number ratio on the object `i_matter`

      CLASS(particles), INTENT(IN):: this
      !! [[particles]] object which this PROCEDURE is a member of
      INTEGER, INTENT(IN):: i_matter
      !! Index of the matter object
      DOUBLE PRECISION:: nuratio
      !! Baryon number ratio on the object `i_matter`

    END FUNCTION get_nuratio_i


    MODULE PURE FUNCTION get_pos( this ) RESULT( pos_u )
    !! Returns [[particles:pos]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(IN):: this
      !> [[particles:pos]]
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos_u

    END FUNCTION get_pos


    MODULE PURE FUNCTION get_vel( this ) RESULT( vel )
    !! Returns [[particles:v]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(IN):: this
      !> [[particles:v]]
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: vel

    END FUNCTION get_vel


    MODULE PURE FUNCTION get_nstar( this ) RESULT( nstar )
    !! Returns [[particles:nstar]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(IN):: this
      !> [[particles:nstar]]
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: nstar

    END FUNCTION get_nstar


    MODULE PURE FUNCTION get_nstar_sph( this ) RESULT( nstar_sph )
    !! Returns [[particles:nstar_sph]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(IN):: this
      !> [[particles:nstar_sph]]
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: nstar_sph

    END FUNCTION get_nstar_sph


    MODULE PURE FUNCTION get_nlrf( this ) RESULT( nlrf )
    !! Returns [[particles:nlrf]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(IN):: this
      !> [[particles:nlrf]]
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: nlrf

    END FUNCTION get_nlrf


    MODULE PURE FUNCTION get_nlrf_sph( this ) RESULT( nlrf_sph )
    !! Returns [[particles:nlrf_sph]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(IN):: this
      !> [[particles:nlrf_sph]]
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: nlrf_sph

    END FUNCTION get_nlrf_sph


    MODULE PURE FUNCTION get_nu( this ) RESULT( nu )
    !! Returns [[particles:nu]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(IN):: this
      !> [[particles:nu]]
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: nu

    END FUNCTION get_nu


    MODULE PURE FUNCTION get_u( this ) RESULT( u )
    !! Returns [[particles:specific_energy]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(IN):: this
      !> [[particles:specific_energy]]
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: u

    END FUNCTION get_u


    MODULE PURE FUNCTION get_u_sph( this ) RESULT( u_sph )
    !! Returns [[particles:u_sph]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(IN):: this
      !> [[particles:u_sph]]
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: u_sph

    END FUNCTION get_u_sph


    MODULE PURE FUNCTION get_pressure( this ) RESULT( pressure )
    !! Returns [[particles:pressure]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(IN):: this
      !> [[particles:pressure]]
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: pressure

    END FUNCTION get_pressure


    MODULE PURE FUNCTION get_pressure_sph( this ) RESULT( pressure_sph )
    !! Returns [[particles:pressure_sph]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(IN):: this
      !> [[particles:pressure_sph]]
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: pressure_sph

    END FUNCTION get_pressure_sph


    MODULE PURE FUNCTION get_theta( this ) RESULT( theta )
    !! Returns [[particles:theta]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(IN):: this
      !> [[particles:theta]]
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: theta

    END FUNCTION get_theta


    MODULE PURE FUNCTION get_h( this ) RESULT( h )
    !! Returns [[particles:h]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(IN):: this
      !> [[particles:h]]
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: h

    END FUNCTION get_h


    MODULE PURE FUNCTION get_lapse( this ) RESULT( lapse )
    !! Returns [[particles:lapse]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(IN):: this
      !> [[particles:lapse]]
      DOUBLE PRECISION, DIMENSION(this% npart):: lapse

    END FUNCTION get_lapse


    MODULE PURE FUNCTION get_shift( this ) RESULT( shift )
    !! Returns \([[particles:shift_x]],[[particles:shift_y]],[[particles:shift_z]]\)

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(IN):: this
      !> ([[particles:shift_x]],[[particles:shift_y]],[[particles:shift_z]])
      DOUBLE PRECISION, DIMENSION(3,this% npart):: shift

    END FUNCTION get_shift


    MODULE PURE FUNCTION get_g3( this ) RESULT( g3 )
    !! Returns the spatial metric on the particles

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(IN):: this
      !> ([[particles:g_xx]],[[particles:g_xy]],[[particles:g_xz]],[[particles:g_yy]],[[particles:g_yz]],[[particles:g_zz]])
      DOUBLE PRECISION, DIMENSION(6,this% npart):: g3

    END FUNCTION get_g3


    MODULE PURE FUNCTION get_compose_eos( this ) RESULT( compose_eos )
    !! Returns [[particles:compose_eos]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT(IN):: this
      !> ([[particles:compose_eos]])
      LOGICAL:: compose_eos

    END FUNCTION get_compose_eos


  END INTERFACE


  !---------------------------------!
  !--  NON-TYPE-BOUND PROCEDURES  --!
  !---------------------------------!


  INTERFACE


    MODULE SUBROUTINE check_particle_positions( npart, pos, debug )
    !! Check that the particles are not at the same positions

      INTEGER, INTENT(IN):: npart
      !! Number of particles
      LOGICAL, INTENT(IN), OPTIONAL:: debug
      !! `TRUE` to debug the SUBROUTINE, `FALSE` otherwise
      DOUBLE PRECISION, DIMENSION(3,npart), INTENT(IN):: pos
      !! Array of particle positions

    END SUBROUTINE check_particle_positions


    MODULE FUNCTION check_particle_position( npart, pos, pos_a ) RESULT( cnt )
    !# Return the number of times that pos_a appears in the array pos
    !  @todo This algorithm scales as O(npart**2)
    !        if used in a loop over the particles...
    !        To be documented, after it's fixed

      INTEGER, INTENT(IN):: npart
      DOUBLE PRECISION, DIMENSION(3,npart), INTENT(IN):: pos
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: pos_a
      INTEGER:: cnt

    END FUNCTION check_particle_position


    MODULE SUBROUTINE correct_center_of_mass &
    ( npart, pos, nu, get_density, validate_pos, com_star, verbose )
    !# Translate the particles so that their center of mass
    !  coincides with the center of mass of the star, given by |id|

      IMPLICIT NONE

      INTEGER, INTENT(IN):: npart
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: com_star
      LOGICAL, INTENT(IN), OPTIONAL:: verbose

      INTERFACE
        FUNCTION get_density( x, y, z ) RESULT( density )
          DOUBLE PRECISION, INTENT(IN):: x
          DOUBLE PRECISION, INTENT(IN):: y
          DOUBLE PRECISION, INTENT(IN):: z
          DOUBLE PRECISION:: density
        END FUNCTION
      END INTERFACE
      INTERFACE
        FUNCTION validate_pos( x, y, z ) RESULT( answer )
          DOUBLE PRECISION, INTENT(IN):: x
          DOUBLE PRECISION, INTENT(IN):: y
          DOUBLE PRECISION, INTENT(IN):: z
          LOGICAL:: answer
        END FUNCTION
      END INTERFACE

      DOUBLE PRECISION, DIMENSION(3,npart), INTENT(INOUT):: pos
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(INOUT):: nu

    END SUBROUTINE correct_center_of_mass


    MODULE SUBROUTINE get_neighbours_bf &
    (ipart,npart,pos,h,dimensions,nnei,neilist)
    !# Get neighbours of particle ipart in a "brute force" way

      IMPLICIT NONE

      INTEGER,          INTENT(IN) :: ipart, npart, dimensions
      DOUBLE PRECISION, INTENT(IN) :: pos(dimensions,npart), h(npart)
      INTEGER,          INTENT(OUT):: nnei, neilist(npart)

    END SUBROUTINE get_neighbours_bf


    MODULE SUBROUTINE compute_and_print_quality_indicators &
    (npart, pos, h, nu, nstar, path)
    !# Compute the quality indicators, referring to
    !
    !  Daniel J. Price, Smoothed Particle Hydrodynamics and
    !  Magnetohydrodynamics.
    !  Journal of Computational Physics, 231, 3, 759-794 (2012)
    !  DOI: 10.1016/j.jcp.2010.12.011
    !  http://arxiv.org/abs/1012.1885
    !  eqs.(64), (67) and (74-75)
    !
    !  Rosswog, S. SPH Methods in the Modelling of Compact Objects.
    !  Living Rev Comput Astrophys 1, 1 (2015).
    !  https://doi.org/10.1007/lrca-2015-1.
    !  https://arxiv.org/abs/1406.4224.
    !  eqs.(6) and (9)

      IMPLICIT NONE

      INTEGER,                                INTENT(IN):: npart
      DOUBLE PRECISION, DIMENSION(3,npart),   INTENT(IN):: pos
      DOUBLE PRECISION, DIMENSION(npart),     INTENT(IN):: h
      DOUBLE PRECISION, DIMENSION(npart),     INTENT(IN):: nu
      DOUBLE PRECISION, DIMENSION(npart),     INTENT(IN):: nstar
      CHARACTER(LEN=*), OPTIONAL,             INTENT(IN):: path
      !# Path to which saving the output file containing the quality
      !  indicators.

    END SUBROUTINE compute_and_print_quality_indicators


  END INTERFACE


  INTERFACE compute_adm_momentum_fluid
  !# Computes an estimate of the \(\mathrm{ADM}\) linear momentum of the fluid

    MODULE PROCEDURE compute_adm_momentum_fluid_fields
    !# Computes an estimate of the \(\mathrm{ADM}\) linear momentum of the
    !  fluid using the canonical momentum per baryon on the particles
    MODULE PROCEDURE compute_adm_momentum_fluid_canmom
    !# Computes an estimate of the \(\mathrm{ADM}\) linear momentum of the
    !  fluid using the |sph| fields

  END INTERFACE compute_adm_momentum_fluid


  INTERFACE


    MODULE SUBROUTINE compute_adm_momentum_fluid_fields &
    ( npart, g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, lapse, &
      shift_x, shift_y, shift_z, nu, theta, nlrf, pr, u, vel_u, adm_mom )
    !# Computes an estimate of the \(\mathrm{ADM}\) linear momentum of the
    !  fluid using the |sph| fields
    !  @todo add reference

      INTEGER,                              INTENT(IN) :: npart
      !! Particle number
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN) :: g_xx
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN) :: g_xy
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN) :: g_xz
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN) :: g_yy
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN) :: g_yz
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN) :: g_zz
      !! Spatial metric on the particles
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN) :: lapse
      !! Lapse function on the particles
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN) :: shift_x
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN) :: shift_y
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN) :: shift_z
      !! Shift vector on the particles
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN) :: nu
      !! Baryon number per particle
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN) :: theta
      !! Generalized Lorentz factor on the particles
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN) :: nlrf
      !! Baryon density in the local rest frame on the particles
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN) :: pr
      !! Pressure on the particles
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN) :: u
      !! Specific internal energy on the particles
      DOUBLE PRECISION, DIMENSION(3,npart), INTENT(IN) :: vel_u
      !! Spatial velocity in the computing frame on the particles
      DOUBLE PRECISION, DIMENSION(3),       INTENT(OUT):: adm_mom
      !! ADM linear momentum

    END SUBROUTINE compute_adm_momentum_fluid_fields


    MODULE SUBROUTINE compute_adm_momentum_fluid_canmom &
    ( npart, g3, lapse, shift, nu, s_l, adm_mom )
    !# Computes an estimate of the \(\mathrm{ADM}\) linear momentum of the
    !  fluid using the canonical momentum per baryon on the particles
    !  @todo add reference

      INTEGER,                              INTENT(IN) :: npart
      !! Particle number
      DOUBLE PRECISION, DIMENSION(6,npart), INTENT(IN) :: g3
      !! Spatial metric on the particles
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN) :: lapse
      !! Lapse function on the particles
      DOUBLE PRECISION, DIMENSION(3,npart), INTENT(IN) :: shift
      !! Shift vector on the particles
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN) :: nu
      !! Baryon number per particle
      DOUBLE PRECISION, DIMENSION(3,npart), INTENT(IN) :: s_l
      !! Canonical momentum per baryon on the particles
      DOUBLE PRECISION, DIMENSION(3),       INTENT(OUT):: adm_mom
      !! ADM linear momentum

    END SUBROUTINE compute_adm_momentum_fluid_canmom


  END INTERFACE



  CONTAINS



  FUNCTION find_h_backup( a, npart, pos, ndes ) RESULT( h )

    !**************************************************************
    !
    !# Backup method to find the smoothing length via brute force
    !  if the optimized method gives 0.
    !  It sets the smoothing lengths to the distance between the
    !  particle and the ndes-th closest neighbour.
    !
    !  FT 24.11.2021
    !
    !**************************************************************

    USE NR,        ONLY: select
    USE constants, ONLY: half

    IMPLICIT NONE

    INTEGER,          INTENT(IN):: a
    !! Index of the particle whose smoothing length is to be computed
    INTEGER,          INTENT(IN):: npart
    !! Number of particles
    INTEGER,          INTENT(IN):: ndes
    !! Desired number of neighbours
    DOUBLE PRECISION, DIMENSION(3,npart), INTENT(IN):: pos
    !! Array containing particle positions

    DOUBLE PRECISION:: h
    !! Smoothing length

    INTEGER:: b
    !! Particle index running over all particles, except particle `a`
    DOUBLE PRECISION, DIMENSION(npart):: dist2
    !! Square norm of the distance vector between the particles `a` and `b`
    DOUBLE PRECISION, DIMENSION(3):: dist
    !! Distance vector between the particles `a` and `b`

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( pos, a, npart, dist2 ) &
    !$OMP             PRIVATE( b, dist )
    DO b= 1, npart, 1

      !IF( a /= b )THEN

        dist(:)= pos(:,b) - pos(:,a)
        dist2(b)= DOT_PRODUCT(dist,dist)

      !ENDIF

    ENDDO
    !$OMP END PARALLEL DO

    ! ndes+1 is used, rather tan ndes, because the particle itself is included
    ! in the distance array dist2
    h= half*SQRT( select(ndes+1, npart, dist2) )

  END FUNCTION find_h_backup


!  SUBROUTINE COM_1PN( npart, pos, vel, nu, nlrf, u, nstar, &
!                      sq_det_g4, g4, &
!                      pn_com_x, pn_com_y, pn_com_z, pn_com_d, mass )
!
!    !************************************************************
!    !                                                           *
!    ! monitor PN baryonic center of mass; PD 30.09.2021         *
!    !                                                           *
!    ! See eqs. (7.4), (8.136), (8.146) and (8.4.7) in           *
!    ! Eric Poisson, Clifford M. Will                            *
!    ! "Gravity: Newtonian, Post-Newtonian, Relativistic",       *
!    ! Cambridge University Press, 2014. FT 26.01.2022           *
!    !                                                           *
!    !************************************************************
!
!    USE tensor,    ONLY: n_sym4x4, ixx, ixy, ixz, iyy, iyz, izz, jx, jy, jz
!    USE constants, ONLY: zero, third, half, one, two, three
!
!    IMPLICIT NONE
!
!    INTEGER,          INTENT(IN) :: npart
!    ! Number of particles
!    DOUBLE PRECISION, INTENT(IN) :: pos(3,npart)
!    ! Particles' positions
!    DOUBLE PRECISION, INTENT(IN) :: vel(3,npart)
!    ! Particles' 3-velocities in the SPH computing frame
!    DOUBLE PRECISION, INTENT(IN) :: nu(npart)
!    ! Particles' baryon numbers
!    DOUBLE PRECISION, INTENT(IN) :: nlrf(npart)
!    ! Particles densities in the local rest frame
!    DOUBLE PRECISION, INTENT(IN) :: u(npart)
!    ! Particles' specific internal energy
!    DOUBLE PRECISION, INTENT(IN) :: nstar(npart)
!    ! Particles' baryon density
!    DOUBLE PRECISION, INTENT(IN) :: sq_det_g4(npart)
!    ! Square root of minus the determinant of the spacetime metric, at the
!    ! particle positions
!    DOUBLE PRECISION, INTENT(IN) :: g4(n_sym4x4,npart)
!    ! Spacetime metric at the particle positions
!    DOUBLE PRECISION, INTENT(OUT):: pn_com_x
!    ! x component of the 1PN center of mass
!    DOUBLE PRECISION, INTENT(OUT):: pn_com_y
!    ! y component of the 1PN center of mass
!    DOUBLE PRECISION, INTENT(OUT):: pn_com_z
!    ! z component of the 1PN center of mass
!    DOUBLE PRECISION, INTENT(OUT):: pn_com_d
!    ! Distance of the 1PN center of mass from the origin
!
!    INTEGER         :: a
!    ! Index running over the particles
!    DOUBLE PRECISION, INTENT(OUT):: mass
!    ! Total mass of the particles
!    DOUBLE PRECISION:: v_sqnorm
!    ! Squared norm of the particles' 3-velocities in the SPH computing frame
!    DOUBLE PRECISION:: u_pot
!    ! Potential U in the formulas. See Exercise 8.1 on p.410 in the reference
!    ! cited in this SUBROUTINE's description
!    DOUBLE PRECISION:: pi_pot
!    ! Potential Pi in the formulas. See line right below eq.(8.7c) on p.373
!    ! in the reference cited in this SUBROUTINE's description
!    DOUBLE PRECISION:: nu_pot
!    ! Potential nu in the formulas. See lines right above eq.(8.148) on p.409
!    ! in the reference cited in this SUBROUTINE's description
!
!    mass    = zero
!    pn_com_x= zero
!    pn_com_y= zero
!    pn_com_z= zero
!
!    DO a= 1, npart, 1
!
!      v_sqnorm= vel(jx,a)**2 + vel(jy,a)**2 + vel(jz,a)**2
!      !v_sqnorm= g4(ixx,a)*vel(jx,a)**2 + two*g4(ixy,a)*vel(jx,a)*vel(jy,a) &
!      !        + two*g4(ixz,a)*vel(jx,a)*vel(jz,a) + g4(iyy,a)*vel(jy,a)**2 &
!      !        + two*g4(iyz,a)*vel(jy,a)*vel(jz,a) + g4(izz,a)*vel(jz,a)**2
!
!      u_pot   = half*( sq_det_g4(a) - one )
!
!      pi_pot  = u(a)*nlrf(a)/nstar(a)
!      !pi_pot  = u(a)*( one - half*v_sqnorm - three*u_pot )
!
!      nu_pot  = one + half*v_sqnorm - half*u_pot + pi_pot
!
!      !mass    = mass + nu(a)*nu_pot
!      mass    = mass + nlrf(a)*( nu(a)/nstar(a) ) &
!                /( one - half*v_sqnorm - three*u_pot )*nu_pot
!
!!IF( one - half*v_sqnorm - three*u_pot < zero )THEN
!!  PRINT *, one - half*v_sqnorm - three*u_pot
!!  PRINT *, v_sqnorm
!!  PRINT *, u_pot
!!  STOP
!!ENDIF
!
!      !pn_com_x= pn_com_x + nu(a)*pos(jx,a)*nu_pot
!      !pn_com_y= pn_com_y + nu(a)*pos(jy,a)*nu_pot
!      !pn_com_z= pn_com_z + nu(a)*pos(jz,a)*nu_pot
!      pn_com_x= pn_com_x + nlrf(a)*( nu(a)/nstar(a) ) &
!                /( one - half*v_sqnorm - three*u_pot )*pos(jx,a)*nu_pot
!      pn_com_y= pn_com_y + nlrf(a)*( nu(a)/nstar(a) ) &
!                /( one - half*v_sqnorm - three*u_pot )*pos(jy,a)*nu_pot
!      pn_com_z= pn_com_z + nlrf(a)*( nu(a)/nstar(a) ) &
!                /( one - half*v_sqnorm - three*u_pot )*pos(jz,a)*nu_pot
!
!    ENDDO
!
!    pn_com_x  = pn_com_x/mass
!    pn_com_y  = pn_com_y/mass
!    pn_com_z  = pn_com_z/mass
!    pn_com_d  = SQRT(pn_com_x**2 + pn_com_y**2 + pn_com_z**2)
!
!  END SUBROUTINE COM_1PN
!
!
!  SUBROUTINE momentum_1pn( npart, pos, vel, nu, rho, u, pressure, nstar, h, &
!                           sq_det_g4, g4, &
!                           p_x, p_y, p_z )
!
!    !************************************************************
!    !                                                           *
!    ! Compute the first-order Post-Newtonian three-momentum     *
!    ! of the spacetime.                                         *
!    !                                                           *
!    ! See eqs. (8.145) in                                       *
!    ! Eric Poisson, Clifford M. Will                            *
!    ! "Gravity: Newtonian, Post-Newtonian, Relativistic",       *
!    ! Cambridge University Press, 2014                          *
!    !                                                           *
!    ! @todo: adapt comments to be read by FORD                  *
!    !        (FORtran Documenter)? This requires minimal        *
!    !        modifications                                      *
!    !                                                           *
!    ! FT 26.01.2022                                             *
!    !                                                           *
!    !************************************************************
!
!    USE tensor,       ONLY: n_sym4x4, ixx, ixy, ixz, iyy, iyz, izz, jx, jy, jz
!    USE constants,    ONLY: zero, third, half, one, two, three, amu, m0c2
!    USE RCB_tree_3D,  ONLY: nfinal, nprev, nic, lpart, rpart, iorig
!    USE SPHINCS_SPH,  ONLY: all_clists, ncand
!    USE kernel_table, ONLY: dv_table, dv_table_1, &
!                            W_no_norm, n_tab_entry!, dWdv_no_norm
!
!    IMPLICIT NONE
!
!    INTEGER,          INTENT(IN) :: npart
!    ! Number of particles
!    DOUBLE PRECISION, INTENT(IN) :: pos(3,npart)
!    ! Particles' positions
!    DOUBLE PRECISION, INTENT(IN) :: vel(3,npart)
!    ! Particles' 3-velocities in the SPH computing frame
!    DOUBLE PRECISION, INTENT(IN) :: nu(npart)
!    ! Particles' baryon numbers
!    DOUBLE PRECISION, INTENT(IN) :: rho(npart)
!    ! Particles densities
!    DOUBLE PRECISION, INTENT(IN) :: u(npart)
!    ! Particles' specific internal energy
!    DOUBLE PRECISION, INTENT(IN) :: pressure(npart)
!    ! Particles' pressure
!    DOUBLE PRECISION, INTENT(IN) :: nstar(npart)
!    ! Particles' baryon density
!    DOUBLE PRECISION, INTENT(IN) :: h(npart)
!    ! Particles' smoothing lengths
!    DOUBLE PRECISION, INTENT(IN) :: sq_det_g4(npart)
!    ! Square root of minus the determinant of the spacetime metric, at the
!    ! particle positions
!    DOUBLE PRECISION, INTENT(IN) :: g4(n_sym4x4,npart)
!    ! Spacetime metric at the particle positions
!    DOUBLE PRECISION, INTENT(OUT):: p_x
!    ! x component of the 1PN spacetime momentum
!    DOUBLE PRECISION, INTENT(OUT):: p_y
!    ! y component of the 1PN spacetime momentum
!    DOUBLE PRECISION, INTENT(OUT):: p_z
!    ! z component of the 1PN spacetime momentum
!
!    INTEGER         :: a
!    ! Index running over the particles
!    INTEGER ill, itot, indexx, index1, l, k, b
!    DOUBLE PRECISION:: mass
!    ! Total mass of the particles
!    DOUBLE PRECISION:: v_sqnorm
!    ! Squared norm of the particles' 3-velocities in the SPH computing frame
!    DOUBLE PRECISION:: u_pot
!    ! Potential U in the formulas. See Exercise 8.1 on p.410 in the reference
!    ! cited in this SUBROUTINE's description
!    DOUBLE PRECISION:: pi_pot
!    ! Potential Pi in the formulas. See line right below eq.(8.7c) on p.373
!    ! in the reference cited in this SUBROUTINE's description
!    DOUBLE PRECISION:: nu_pot
!    ! Potential nu in the formulas. See lines right above eq.(8.148) on p.409
!    ! in the reference cited in this SUBROUTINE's description
!    DOUBLE PRECISION:: phi_pot(3,npart)
!    ! Nonlocal potential phi in the formulas. See eq.(8.8) on p.374 in the
!    ! reference cited in this SUBROUTINE's description
!    !DOUBLE PRECISION:: phi_pot_integrand(3,npart)
!    ! Integrand of the nonlocal potential phi in the formulas. See eq.(8.8)
!    ! on p.374 in the reference cited in this SUBROUTINE's description
!    DOUBLE PRECISION:: vel_cov(0:3,npart)
!    ! Covariant particles' 4-velocities in the SPH computing frame
!    DOUBLE PRECISION:: dx, dy, dz, va, ha, ha_1, ha_3, phi_pot_integ(3), &
!                       Wi, Wi1, dvv, Wab_ha
!
!    PRINT *, "nu      =", nu      (1)
!    PRINT *, "rho     =", rho     (1)
!    PRINT *, "u       =", u       (1)
!    PRINT *, "pressure=", pressure(1)
!    PRINT *, "nstar   =", nstar   (1)
!    PRINT *
!    !STOP
!
!    !
!    !-- Computation of the nonlocal potential phi_pot
!    !
!
!    phi_pot= zero
!
!   ! DO a= 1, npart, 1
!   !
!   !   ! For each a, we should sum over the other particles (neighbours?)
!   !   phi_pot= phi_pot + phi_pot_integrand(a,b)
!   !
!   ! ENDDO
!
!    PRINT *, "** Computing SPH integral estimate of the vector potential ", &
!             "Phi^i..."
!
!    !$OMP PARALLEL DO DEFAULT(SHARED) &
!    !$OMP PRIVATE(ill,itot) &
!    !$OMP SCHEDULE(STATIC)
!    ll_cell_loop: DO ill= 1, nfinal, 1
!
!     ! if empty: skip
!     itot= ill + nprev
!     IF( nic(itot) == 0 ) CYCLE
!
!     ! particle content in this cell
!     particle_loop: DO l= lpart(itot), rpart(itot)
!
!      a   = iorig(l)
!
!      ha  = h(a)
!      ha_1= one/ha
!      ha_3= ha_1*ha_1*ha_1
!
!      cand_loop: DO k= 1, ncand(ill), 1
!
!        b= all_clists(ill)% list(k)
!
!        ! Distances (ATTENTION: flatspace version !!!)
!        dx= pos(1,a) - pos(1,b)
!        dy= pos(2,a) - pos(2,b)
!        dz= pos(3,a) - pos(3,b)
!        va= SQRT(dx*dx + dy*dy + dz*dz)*ha_1
!
!        ! get interpolation indices
!        indexx= MIN( INT(va*dv_table_1), n_tab_entry )
!        index1= MIN( indexx + 1, n_tab_entry )
!
!        ! get tabulated values
!        Wi = W_no_norm(indexx)
!        Wi1= W_no_norm(index1)
!
!        ! interpolate
!        dvv   = ( va - DBLE(indexx)*dv_table )*dv_table_1
!        Wab_ha= Wi + ( Wi1 - Wi )*dvv
!
!        ! sum up for number density
!        phi_pot(1,a)= phi_pot(1,a) + phi_pot_integrand(1,a,b)*Wab_ha &
!                                    *nu(a)/nstar(a)
!        phi_pot(2,a)= phi_pot(2,a) + phi_pot_integrand(2,a,b)*Wab_ha &
!                                    *nu(a)/nstar(a)
!        phi_pot(3,a)= phi_pot(3,a) + phi_pot_integrand(3,a,b)*Wab_ha &
!                                    *nu(a)/nstar(a)
!
!        IF( ISNAN( phi_pot(1,a) ) .OR. ISNAN( phi_pot(2,a) ) &
!          .OR. ISNAN( phi_pot(3,a) ) )THEN
!
!          PRINT *
!          PRINT *, "phi_pot_integrand(1,a,b)= ", phi_pot_integrand(1,a,b)
!          PRINT *, "phi_pot_integrand(2,a,b)= ", phi_pot_integrand(2,a,b)
!          PRINT *, "phi_pot_integrand(3,a,b)= ", phi_pot_integrand(3,a,b)
!          PRINT *, "Wab_ha= ", Wab_ha
!          PRINT *, "dvv= ", dvv
!          PRINT *, "dx= ", dx
!          PRINT *, "dy= ", dy
!          PRINT *, "dz= ", dz
!          PRINT *, "phi_pot(:,a)= ", phi_pot(:,a)
!          PRINT *
!          STOP
!
!        ENDIF
!
!      ENDDO cand_loop
!
!      ! normalize with ha's
!      phi_pot(:,a)= phi_pot(:,a)!*ha_3
!
!     ENDDO particle_loop
!
!    ENDDO ll_cell_loop
!    !$OMP END PARALLEL DO
!    PRINT *, "...done."
!    PRINT *
!
!    !
!    !-- Computation of the local potentials, and the spacetime momentum
!    !
!
!    p_x = zero
!    p_y = zero
!    p_z = zero
!    mass= zero
!
!    DO a= 1, npart, 1
!
!      v_sqnorm= vel(jx,a)**two + vel(jy,a)**two + vel(jz,a)**two
!      !v_sqnorm= g4(ixx,a)*vel(jx,a)**two + two*g4(ixy,a)*vel(jx,a)*vel(jy,a) &
!      !        + two*g4(ixz,a)*vel(jx,a)*vel(jz,a) + g4(iyy,a)*vel(jy,a)**two &
!      !        + two*g4(iyz,a)*vel(jy,a)*vel(jz,a) + g4(izz,a)*vel(jz,a)**two
!
!      u_pot= half*( sq_det_g4(a) - one )
!
!      pi_pot= u(a)*( one - half*v_sqnorm - three*u_pot )!*rho(a)/nstar(a)
!
!      nu_pot= one + half*v_sqnorm - half*u_pot + pi_pot
!
!      mass= mass + nu(a)*nu_pot
!
!      p_x= p_x + nu(a)*vel(jx,a)*( nu_pot + pressure(a)/nstar(a) ) &
!               - half*phi_pot(jx,a)*nu(a)/nstar(a)
!      p_y= p_y + nu(a)*vel(jy,a)*( nu_pot + pressure(a)/nstar(a) ) &
!               - half*phi_pot(jy,a)*nu(a)/nstar(a)
!      p_z= p_z + nu(a)*vel(jz,a)*( nu_pot + pressure(a)/nstar(a) ) &
!               - half*phi_pot(jz,a)*nu(a)/nstar(a)
!
!     ! PRINT *, "u_pot =", u_pot
!     ! PRINT *, "pi_pot=", pi_pot
!     ! PRINT *, "nu_pot=", nu_pot
!     ! PRINT *, "mass  =", mass
!     ! PRINT *, "pressure(a)/nstar(a)= ", pressure(a)/nstar(a)
!     ! PRINT *, "phi_pot(:,a)=", phi_pot(:,a)
!     ! PRINT *
!     ! STOP
!
!    ENDDO
!
!    PRINT *, "u_pot =", u_pot
!    PRINT *, "pi_pot=", pi_pot
!    PRINT *, "nu_pot=", nu_pot
!    PRINT *, "mass  =", mass
!    PRINT *, "phi_pot(:,100)=", phi_pot(:,100)
!    PRINT *, "pressure(100)/nstar(100)= ", pressure(100)/nstar(100)
!    PRINT *, "p/M=", p_x/mass, p_y/mass, p_z/mass
!    PRINT *
!    !STOP
!
!
!    CONTAINS
!
!
!    FUNCTION phi_pot_integrand(i,a,b) RESULT(res)
!
!      !************************************************************
!      !                                                           *
!      ! Integrand of the nonlocal potential phi in the formulas.  *
!      ! See eq.(8.8) on p.374 in the reference cited in this      *
!      ! SUBROUTINE's description                                  *
!      !                                                           *
!      ! FT 26.01.2022                                             *
!      !                                                           *
!      !************************************************************
!
!      USE tensor,    ONLY: lower_index_4vector
!      USE constants, ONLY: zero
!
!      IMPLICIT NONE
!
!      INTEGER, INTENT(IN):: a
!      ! Index of the particle with non-dummy position
!      INTEGER, INTENT(IN):: b
!      ! Index of the particle with dummy position (the dummy position is
!      ! the integratin variable)
!      INTEGER, INTENT(IN):: i
!      ! Spatial index
!      DOUBLE PRECISION:: res
!      ! Vector integrand
!
!      IF( a == b )THEN
!        res= zero
!        RETURN
!      ENDIF
!
!      CALL lower_index_4vector( [one,vel(:,a)], g4(:,a), vel_cov(:,a) )
!
!      res= nu(b)*( vel_cov(jx,b)*(pos(jx,a)-pos(jx,b)) + &
!      !res= ( vel_cov(jx,b)*(pos(jx,a)-pos(jx,b)) + &
!                   vel_cov(jy,b)*(pos(jy,a)-pos(jy,b)) + &
!                   vel_cov(jz,b)*(pos(jz,a)-pos(jz,b)) &
!                 )* &
!                 (pos(i,a)-pos(i,b))/(NORM2(pos(:,a)-pos(:,b))**three)
!
!    END FUNCTION
!
!
!  END SUBROUTINE momentum_1pn


END MODULE sph_particles
