! File:         module_particles_id.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

MODULE particles_id


  !***********************************************************
  !
  !# This module contains the definition of TYPE particles
  !
  !***********************************************************


  USE utility,        ONLY: itr, ios, err_msg, test_status, &
                            perc, creturn, run_id, show_progress
  USE id_base,        ONLY: idbase
  USE timing,         ONLY: timer


  IMPLICIT NONE


  INTEGER, PARAMETER:: id_particles_from_file            = 0
  !! Identifier for a particle distribution read from formatted file
  INTEGER, PARAMETER:: id_particles_on_lattice           = 1
  !! Identifier for a prticle distribution on a lattice
  INTEGER, PARAMETER:: id_particles_on_spherical_surfaces= 2
  !! Identifier for particle distribution on spherical surfaces


  TYPE eos
  !! Data structure representing an |eos|
    CHARACTER( LEN= : ), ALLOCATABLE:: eos_name
    !! The |eos| name
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: eos_parameters
    !# Array containing the |eos| parameters, in the following order:
    !
    !  - For a polytropic |eos|:
    !    [ eos identification number, \(\gamma\),\(\kappa\) ]
    !  - For a piecewise polytropic |eos|:
    !    [ eos identification number, number of polytropic pieces, \(\gamma_0\),
    !   \(\gamma_1\), \(\gamma_2\), \(\gamma_3\), \(\kappa_0\), \(\kappa_1\),
    !    \(\kappa_2\), \(\kappa_3\), \(\log(p_1)\), \(\log(\rho_0)\),
    !    \(\log(\rho_1)\), \(\log(\rho_2)\) ]
    !  - For a tabulated |eos|: [ eos identification number ]
    !
    ! \(\gamma\) is the adimensional polytropic exponent, \(\kappa\) is the
    ! polytropic constant in units of
    ! \(\left(M_\odot L_\odot^{-3}\right)^{1-\gamma}\). Pressure and baryon
    ! mass density have the same units \(M_\odot L_\odot^{-3}\) since \(c^2=1\).
  END TYPE


  !**********************************************************
  !                                                         *
  !              Definition of TYPE particles               *
  !                                                         *
  ! This class places the SPH particles, imports            *
  ! the |id| on the particle positions, stores              *
  ! it, computes the relevant SPH fields and exports it to  *
  ! both a formatted, and a binary file for evolution       *
  !                                                         *
  !**********************************************************

  TYPE:: particles
  !! TYPE representing a |sph| particle distribution


    PRIVATE


    INTEGER:: npart
    !! Total particle number
    INTEGER:: n_matter
    !! Number of matter objects in the physical system
    INTEGER, DIMENSION(:), ALLOCATABLE:: npart_i
    !! Array storing the particle numbers for the matter objects
    INTEGER:: distribution_id
    !! Identification number for the particle distribution
    INTEGER:: call_flag= 0
    !# Flag that is set different than 0 if the SUBROUTINE
    !  compute_and_export_SPH_variables is called

    INTEGER, DIMENSION(:), ALLOCATABLE:: baryon_density_index
    !# Array storing the indices to use with [[particles:baryon_density_parts]]
    !  to sort the elements of [[particles:baryon_density_parts]] in increasing
    !  order

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
    !DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_parts_x1
    !& 1-D array storing the pressure on the x axis
    !  \([\mathrm{kg}\,c^2\,\mathrm{m}^{-3}]\) for NS 2
    !DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_parts_x2
    !& 1-D array storing the first derivative of the pressure
    !  along the x axis \([\mathrm{kg}\,c^2\,\mathrm{m}^{-3}]\) for NS 1
    !DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_parts_x_der1
    !& 1-D array storing the first derivative of the pressure
    !  along the x axis \([\mathrm{kg}\,c^2\,\mathrm{m}^{-3}]\) for NS2
    !DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_parts_x_der2
    !> 1-D array storing the typical length scale for the pressure change
    !DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_length_scale_x1
    !> 1-D array storing the typical length scale for the pressure change
    !DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_length_scale_x2
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
    !  the |lorene| density
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
    !> Baryonic masses of the matter objects \(M_\odot\)
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: masses
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

    !
    !-- Strings
    !

    !> String containing the name of the particles parameter file
    CHARACTER( LEN= 50 ):: lorene_bns_id_parfile

    !& String storing the local path to the directory containing the
    !  CompOSE |eos|
    CHARACTER( LEN= : ), ALLOCATABLE:: compose_path
    !> String storing the subpath of compose_path to the CompOSE file with
    !  .beta extension
    CHARACTER( LEN= : ), ALLOCATABLE:: compose_filename

    !
    !-- Equations of state
    !

    TYPE(eos), DIMENSION(:), ALLOCATABLE:: all_eos

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
    !& `.TRUE.` if the particle positions on spherical surfaces are randomized
    !  in the \(\phi\) direction, `.FALSE.` otherwise
    LOGICAL:: randomize_phi
    !& `.TRUE.` if the particle positions on spherical surfaces are randomized
    !  in the \(\theta\) direction, `.FALSE.` otherwise
    LOGICAL:: randomize_theta
    !& `.TRUE.` if the particle positions on spherical surfaces are randomized
    !  in the \(r\) direction, `.FALSE.` otherwise
    LOGICAL:: randomize_r
    !& `.TRUE.` if the Artificial Pressure Method (APM) has to be applied to the
    !  particles on star 1, `.FALSE.` otherwise
    LOGICAL, DIMENSION(:), ALLOCATABLE:: apm_iterate
    !& `.TRUE.` if the baryon number per particle \(\nu\) has to be read from the
    !  formatted file containing the particle positions, `.FALSE.` otherwise
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

    PROCEDURE:: place_particles_lattice
    !! Places particles on a single lattice that surrounds both stars

    PROCEDURE:: place_particles_spherical_surfaces
    !! Places particles on spherical surfaces on one star

    PROCEDURE, NOPASS:: impose_equatorial_plane_symmetry
    !# Reflects the positions of the particles on a matter object with respect
    !  to the \(xy\) plane

    PROCEDURE, NOPASS:: perform_apm
    !! Performs the Artificial Pressure Method (APM) on one star's particles

  !  GENERIC:: reshape_sph_field => reshape_sph_field_1d_ptr, &
  !                                 reshape_sph_field_2d_ptr
  !  !# GENERIC PROCEDURE, overloded to reallocate 1d and 2d arrays
  !  PROCEDURE:: reshape_sph_field_1d_ptr => reshape_sph_field_1d
  !  !! Reallocates a 1d array
  !  PROCEDURE:: reshape_sph_field_2d_ptr => reshape_sph_field_2d
  !  !! Reallocates a 2d array

    PROCEDURE:: allocate_lorene_id_parts_memory
    !! Allocates memory for the [[particles]] member arrays

    PROCEDURE:: deallocate_lorene_id_parts_memory
    !! Deallocates memory for the [[particles]] member arrays

    PROCEDURE:: read_compose_composition
    !! Reads the \(Y_e(n_b)\) table in the CompOSE file with extension .beta

    PROCEDURE:: compute_Ye
    !# Interpates linearly the electron fraction \(Y_e\) at the particle
    !  densities; that is, assigns \(Y_e\) at the particle positions

    PROCEDURE, PUBLIC:: analyze_hydro
    !# Scans the hydro fields taken from \(\texttt{|lorene|}\) to look
    !  for negative or zero values

    PROCEDURE, PUBLIC:: compute_and_export_SPH_variables
    !# Computes the SPH variables at the particle positions, and optionally
    !  prints them to a binary file to be read by \(\texttt{SPHINCS_BSSN}\)
    !  and \(\texttt{splash}\), and to a formatted file to be read by
    !  \(\texttt{gnuplot}\), by calling
    !  [[particles:print_formatted_lorene_id_particles]]

    PROCEDURE, PUBLIC:: read_sphincs_dump_print_formatted
    !# Reads the binary ID file printed by
    !  [[particles:compute_and_export_SPH_variables]]

    PROCEDURE, PUBLIC:: print_formatted_lorene_id_particles
    !! Prints the SPH ID to a formatted file

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
    PROCEDURE, PUBLIC:: get_nuratio
    !! Returns [[particles:nuratio]]
    PROCEDURE, PUBLIC:: get_nuratio_i
    !! Returns the baryon number ratio on the object `i_matter`
    PROCEDURE, PUBLIC:: get_pos
    !! Returns [[particles:pos]]
    PROCEDURE, PUBLIC:: get_vel
    !! Returns [[particles:v]]
    PROCEDURE, PUBLIC:: get_nlrf
    !! Returns [[particles:nlrf]]
    PROCEDURE, PUBLIC:: get_nu
    !! Returns [[particles:nu]]
    PROCEDURE, PUBLIC:: get_u
    !! Returns [[particles:specific_energy_parts]]
    PROCEDURE, PUBLIC:: get_pressure
    !! Returns [[particles:pressure_parts]]
    PROCEDURE, PUBLIC:: get_pressure_cu
    !! Returns [[particles:pressure_parts_cu]]
    PROCEDURE, PUBLIC:: get_theta
    !! Returns [[particles:theta]]
    PROCEDURE, PUBLIC:: get_h
    !! Returns [[particles:h]]

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

    MODULE PROCEDURE construct_particles
    !! Constructs a [[particles]] object

  END INTERFACE particles

  !
  !-- Interface of the constructor of TYPE particles
  !-- Its implementation is in submodule_particles_constructor.f90
  !
  INTERFACE

    MODULE FUNCTION construct_particles( id, dist ) RESULT ( parts )
    !! Constructs a [[particles]] object

        CLASS(idbase), INTENT( IN OUT ):: id
        !# [[idbase]] object representing the BNS for which we want to place
        !  particles
        INTEGER,    INTENT( IN )    :: dist
        !# Identifier of the desired particle distribution:
        !
        !  - 0: Read particle positions (and optionally the baryon number per
        !     particle \(\nu\)) from a formatted file
        !
        !  - 1: Place particles on a single lattice that surrounds both stars
        !
        !  - 2: Place particles on two lattices, each one surrounding a star
        !
        !  - 3: Place particles on spherical surfaces inside the stars
        !
        !  @warning Method 1 is almost deprecated, since method 2 is effectively
        !           an improvement of method 1
        TYPE(particles)             :: parts
        !! Constructed [[particles]] object

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


    MODULE SUBROUTINE place_particles_lattice( THIS, &
                                  central_density, &
                                  xmin, xmax, ymin, ymax, zmin, zmax, &
                                  npart_des, npart_out, stretch, &
                                  thres, pos, pvol, &
                                  get_density, validate_position )
    !! Places particles on a lattice containing a physical object

      CLASS(particles), INTENT( IN OUT ):: THIS
      !! [[particles]] object which this PROCEDURE is a member of
      DOUBLE PRECISION, INTENT( IN )    :: central_density
      !! Maximum baryon mass density of the system
      INTEGER,          INTENT( IN )    :: npart_des
      !! Desired particle number
      DOUBLE PRECISION, INTENT( IN )    :: xmin
      !! Left \(x\) boundary of the lattice
      DOUBLE PRECISION, INTENT( IN )    :: xmax
      !! Right \(x\) boundary of the lattice
      DOUBLE PRECISION, INTENT( IN )    :: ymin
      !! Left \(y\) boundary of the lattice
      DOUBLE PRECISION, INTENT( IN )    :: ymax
      !! Right \(y\) boundary of the lattice
      DOUBLE PRECISION, INTENT( IN )    :: zmin
      !! Left \(z\) boundary of the lattice
      DOUBLE PRECISION, INTENT( IN )    :: zmax
      !! Right \(z\) boundary of the lattice
      DOUBLE PRECISION, INTENT( IN )    :: stretch
      !! Stretching factor fo the lattice. `xmin` to `zmax` are multiplied by it
      DOUBLE PRECISION, INTENT( IN )    :: thres
      !# (~rho_max)/thres is the minimum mass density considered
      ! when placing particles. Used only when redistribute_nu is
      ! .FALSE. . When redistribute_nu is .TRUE. thres= 100*nu_ratio
      INTEGER,          INTENT( OUT )    :: npart_out
      !! Real, output particle number
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT( OUT ):: pos
      !> Array soring the inal particle volumes
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT( OUT ):: pvol
      !! Array storing the final particle volumes
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
        FUNCTION validate_position_int( x, y, z ) RESULT( answer )
        !! Returns 1 if the position is not valid, 0 otherwise
          DOUBLE PRECISION, INTENT(IN):: x
          !! \(x\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: y
          !! \(y\) coordinate of the desired point
          DOUBLE PRECISION, INTENT(IN):: z
          !! \(z\) coordinate of the desired point
          INTEGER:: answer
          !! 1 if the position is not valid, 0 otherwise
        END FUNCTION validate_position_int
      END INTERFACE
      !> Returns 1 if the position is not valid, 0 otherwise
      PROCEDURE(validate_position_int), OPTIONAL:: validate_position

    END SUBROUTINE place_particles_lattice


    MODULE SUBROUTINE place_particles_spherical_surfaces( THIS, &
                                  mass_star, radius, center, &
                                  central_density, npart_des, &
                                  npart_out, pos, pvol, pmass, &
                                  last_r, upper_bound, lower_bound, &
                                  upper_factor, lower_factor, max_steps, &
                                  filename_mass_profile, &
                                  filename_shells_radii, &
                                  filename_shells_pos, &
                                  get_density, integrate_density, &
                                  get_id, validate_position, pmass_des )
    !! Places particles on spherical surfaces on one star

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT( IN OUT ):: THIS
      !& [[idbase]] object needed to access the BNS data
      !  @TODO Remove the [[idbase]] argument as done in SUBROUTINE perform_apm
      !CLASS(idbase),       INTENT( IN OUT ):: id
      !> Approximate particle number on the star
      INTEGER,          INTENT( IN )    :: npart_des
      !> Final number of particles on the star
      INTEGER,          INTENT( OUT )   :: npart_out
      !& If, after max_steps, the iteration did not converge,
      !  multiply upper_bound by upper_factor, and lower_bound
      !  by lower_factor. max_steps >= 10. 100 is a nice value
      INTEGER,          INTENT( IN )    :: max_steps
      !> Baryonic mass of the star
      DOUBLE PRECISION, INTENT( IN )    :: mass_star
      !> Radius of the star in the x direction towards the companion
      DOUBLE PRECISION, INTENT( IN )    :: radius
      !& \(x|) coordinate of the center of the star, i.e.,
      !  of the point with highest density
      DOUBLE PRECISION, INTENT( IN )    :: center
      !> Central density of the star, i.e., highest density
      DOUBLE PRECISION, INTENT( IN )    :: central_density
      !> Radius of the last spherical surface
      DOUBLE PRECISION, INTENT( IN )    :: last_r
      !& If, after max_steps, the iteration did not converge,
      !  multiply upper_bound by upper_factor, and lower_bound
      !  by lower_factor. upper_factor >= 1, usually an increase of 1% works
      DOUBLE PRECISION, INTENT( IN )    :: upper_factor
      !& If, after max_steps, the iteration did not converge,
      !  multiply upper_bound by upper_factor, and lower_bound
      !  by lower_factor. lower_factor <= 1, usually a decrease of 1% works
      DOUBLE PRECISION, INTENT( IN )    :: lower_factor
      !& Desired upper bound for the differences between particle
      !  masses on neighbouring spherical surfaces
      DOUBLE PRECISION, INTENT( INOUT ) :: upper_bound
      !& Desired lower bound for the differences between particle
      !  masses on neighbouring spherical surfaces
      DOUBLE PRECISION, INTENT( INOUT ) :: lower_bound
      !> Array string the final positions
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT( INOUT ):: pos
      !> Array soring the inal particle volumes
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, INTENT( INOUT ):: pvol
      !> Array storing the final particle masses
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, INTENT( INOUT ):: pmass
      !> Name of the file to store the radial mass profile
      CHARACTER( LEN= * ), INTENT( INOUT ), OPTIONAL :: filename_mass_profile
      !& Name of the file to store the surface radii
      !  @TODO change name of variable to filename_surfaces_radii
      CHARACTER( LEN= * ), INTENT( INOUT ), OPTIONAL :: filename_shells_radii
      !& Name of the file to store the final particle positions
      !  @TODO change name of variable to filename_surfaces_pos
      CHARACTER( LEN= * ), INTENT( INOUT ), OPTIONAL :: filename_shells_pos
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
          DOUBLE PRECISION, INTENT( OUT ):: sqdetg
          DOUBLE PRECISION, INTENT( OUT ):: baryon_density
          DOUBLE PRECISION, INTENT( OUT ):: gamma_euler
        END SUBROUTINE get_id
      END INTERFACE
      INTERFACE
        SUBROUTINE integrate_density( center, radius, &
                                      central_density, &
                                      dr, dth, dphi, &
                                      mass, mass_profile, &
                                      mass_profile_idx )
          !& Array to store the indices for array mass_profile, sorted so that
          !  mass_profile[mass_profile_idx] is in increasing order
          INTEGER, DIMENSION(:), ALLOCATABLE, INTENT( IN OUT )::mass_profile_idx
          !> Center of the star
          DOUBLE PRECISION, INTENT( IN )    :: center
          !> Central density of the star
          DOUBLE PRECISION, INTENT( IN )    :: central_density
          !> Radius of the star
          DOUBLE PRECISION, INTENT( IN )    :: radius
          !> Integration steps
          DOUBLE PRECISION, INTENT( IN )    :: dr, dth, dphi
          !> Integrated mass of the star
          DOUBLE PRECISION, INTENT( IN OUT ):: mass
          !> Array storing the radial mass profile of the star
          DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT( IN OUT ):: &
                                           mass_profile
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
          INTEGER:: answer
          !! 1 if the position is not valid, 0 otherwise
        END FUNCTION validate_position_int
      END INTERFACE
      !> Returns 1 if the position is not valid, 0 otherwise
      PROCEDURE(validate_position_int), OPTIONAL:: validate_position
      DOUBLE PRECISION, INTENT( IN ),   OPTIONAL:: pmass_des


    END SUBROUTINE place_particles_spherical_surfaces


    MODULE SUBROUTINE impose_equatorial_plane_symmetry( npart, pos, &
                                                 nu, com_star, verbose )
    !# Mirror the particle with z>0 with respect to the xy plane,
    !  to impose the equatorial-plane symmetry

      !CLASS(particles), INTENT( IN OUT ):: THIS
      INTEGER, INTENT(INOUT):: npart
      DOUBLE PRECISION, INTENT(IN), OPTIONAL:: com_star
      LOGICAL, INTENT(IN), OPTIONAL:: verbose

      DOUBLE PRECISION, DIMENSION(3,npart), INTENT(INOUT):: pos
      DOUBLE PRECISION, DIMENSION(npart),   INTENT(INOUT), OPTIONAL:: nu

    END SUBROUTINE impose_equatorial_plane_symmetry


!    MODULE SUBROUTINE reshape_sph_field_1d( THIS, field, new_size1, new_size2, &
!                                            index_array )
!    !! Reallocates a 1d array
!
!      !> [[particles]] object which this PROCEDURE is a member of
!      CLASS(particles), INTENT( IN OUT ):: THIS
!      !> New particle number on star 1
!      INTEGER,                        INTENT( IN ):: new_size1
!      !> New particle number on star 2
!      INTEGER,                        INTENT( IN ):: new_size2
!      !> Array to select elements to keep in the reshaped array
!      INTEGER,          DIMENSION(:), INTENT( IN ):: index_array
!      !> 1D array to reshape
!      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT( IN OUT ):: field
!
!    END SUBROUTINE reshape_sph_field_1d
!
!
!    MODULE SUBROUTINE reshape_sph_field_2d( THIS, field, new_size1, new_size2, &
!                                            index_array )
!    !! Reallocates a 2d array
!
!      !> [[particles]] object which this PROCEDURE is a member of
!      CLASS(particles), INTENT( IN OUT ):: THIS
!      !> New particle number on star 1
!      INTEGER,                        INTENT( IN ):: new_size1
!      !> New particle number on star 2
!      INTEGER,                        INTENT( IN ):: new_size2
!      !> Array to select elements to keep in the reshaped array
!      INTEGER,          DIMENSION(:), INTENT( IN ):: index_array
!      !> 2D array to reshape
!      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT( IN OUT ):: field
!
!    END SUBROUTINE reshape_sph_field_2d


    MODULE SUBROUTINE allocate_lorene_id_parts_memory( THIS )
    !! Allocates allocatable arrays member of a [[particles]] object

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT( IN OUT ):: THIS

    END SUBROUTINE allocate_lorene_id_parts_memory


    MODULE SUBROUTINE deallocate_lorene_id_parts_memory( THIS )
    !! Deallocates allocatable arrays member of a [[particles]] object

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT( IN OUT ):: THIS

    END SUBROUTINE deallocate_lorene_id_parts_memory


  END INTERFACE

  !
  !-- Interfaces of the methods of TYPE particles
  !-- Their implementations are in module_particles_methods.f90
  !
  INTERFACE

    MODULE SUBROUTINE analyze_hydro( THIS, namefile )
    !# Scans the hydro fields taken from \(\texttt{|lorene|}\) to look
    !  for negative or zero values

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles),    INTENT( IN OUT ):: THIS
      !& Name of the formatted file where the particle positions at which
      !  some of the hydro fields are negative or zero are printed to
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE analyze_hydro

    MODULE SUBROUTINE compute_and_export_SPH_variables( THIS, namefile )
    !# Computes the SPH variables at the particle positions, and optionally
    !  prints them to a binary file to be read by \(\texttt{SPHINCS_BSSN}\)
    !  and \(\texttt{splash}\), and to a formatted file to be read by
    !  \(\texttt{gnuplot}\), by calling
    !  [[particles:print_formatted_lorene_id_particles]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles),    INTENT( IN OUT ):: THIS
      !> Name of the formatted file where the SPH ID is printed to
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE compute_and_export_SPH_variables

    MODULE SUBROUTINE perform_apm( get_density, get_nstar_p, &
                                   pos_input, &
                                   pvol, h_output, nu_output, &
                                   center, &
                                   com_star, &
                                   mass, &
                                   sizes, &!radx_opp, &
                                   !rady, radz, &
                                   apm_max_it, max_inc, &
                                   mass_it, correct_nu, nuratio_thres, &
                                   nuratio_des, &
                                   nx_gh, ny_gh, nz_gh, &
                                   use_atmosphere, &
                                   namefile_pos_id, namefile_pos, &
                                   namefile_results, &
                                   validate_position )
    !! Performs the Artificial Pressure Method (APM) on one star's particles

      !> [[particles]] object which this PROCEDURE is a member of
      !CLASS(particles),                 INTENT( INOUT ):: THIS
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
        SUBROUTINE get_nstar_p( npart_real, x, y, z, nstar_p )
        !! Computes the proper baryon number density at the particle positions
          INTEGER, INTENT(IN):: npart_real
          !! Number of real particles (i.e., no ghost particles included here)
          DOUBLE PRECISION, INTENT(IN):: x(npart_real)
          !! Array of \(x\) coordinates
          DOUBLE PRECISION, INTENT(IN):: y(npart_real)
          !! Array of \(y\) coordinates
          DOUBLE PRECISION, INTENT(IN):: z(npart_real)
          !! Array of \(z\) coordinates
          DOUBLE PRECISION, INTENT(OUT):: nstar_p(npart_real)
          !! Array to store the computed proper baryon number density
        END SUBROUTINE get_nstar_p
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
          INTEGER:: answer
          !! 1 if the position is not valid, 0 otherwise
        END FUNCTION validate_position_int
      END INTERFACE
      !> Returns 1 if the position is not valid, 0 otherwise
      PROCEDURE(validate_position_int), OPTIONAL:: validate_position
      !> Initial particle positions
      DOUBLE PRECISION, DIMENSION(:,:), INTENT( INOUT ):: pos_input
      !> Initial particle volume
      DOUBLE PRECISION, DIMENSION(:),   INTENT( INOUT ):: pvol
      !& Array to store the smoothing lengths computed at the end of the
      !  APM iteration
      DOUBLE PRECISION, DIMENSION(:),   INTENT( OUT )  :: h_output
      !& Array to store the baryon number per particle computed at the end of
      !  the APM iteration
      DOUBLE PRECISION, DIMENSION(:),   INTENT( OUT )  :: nu_output
      !> Center of the star (point of highest density), computed by |lorene|
      DOUBLE PRECISION, DIMENSION(3),   INTENT( IN )   :: center
      !> Center of mass of the star, computed by |lorene|
      DOUBLE PRECISION, DIMENSION(3),   INTENT( INOUT ):: com_star
      !> Mass of the star
      DOUBLE PRECISION,                 INTENT( IN )   :: mass
      !> Radius of the star in the x direction, towards the companion
      DOUBLE PRECISION, DIMENSION(6),   INTENT( IN )   :: sizes
      !> Radius of the star in the x direction, opposite to companion
      !DOUBLE PRECISION,                 INTENT( IN )   :: radx_opp
      !> Radius of the star in the y direction
      !DOUBLE PRECISION,                 INTENT( IN )   :: rady
      !> Radius of the star in the z direction
      !DOUBLE PRECISION,                 INTENT( IN )   :: radz
      !> Maximum number of APM iterations, irrespective of the EXIT condition
      INTEGER,                          INTENT( IN )   :: apm_max_it
      !& Sets the EXIT condition: If the average over all the
      !  particles of the relative error in the density estimate
      !  grows max_inc times, exit the iteration.
      INTEGER,                          INTENT( IN )   :: max_inc
      !& If .TRUE. performs a second iteration after the APM one, without moving
      !  the particles, changing their mass in order to better match
      !  the star density. The mass ratio grows very fast in all the tried
      !  experiments, hence the suggested value is .FALSE.
      LOGICAL,                          INTENT( IN )   :: mass_it
      !& If .TRUE., the baryon number per particle nu is corrected
      !  to include the total baryonic masses of the
      !  stars.
      LOGICAL,                          INTENT( IN )   :: correct_nu
      !& Maximum mass ratio (equivalently baryon number ratio)
      !  to be used in the one-time-only final correction
      !  of the particle masses to match the star density even
      !  better (without moving the particles)
      DOUBLE PRECISION,                 INTENT( IN )   :: nuratio_thres
      !& Sets the EXIT condition: If the baryon number ratio
      !  is within 2.5% of nuratio_des, exit the iteration
      !  Set nuratio_des to 0 to deactivate and exit the APM
      !  iteration using max_inc
      DOUBLE PRECISION,                 INTENT( IN )   :: nuratio_des
      !> Number of lattice points in the x direction for ghosts
      INTEGER,                          INTENT( IN )   :: nx_gh
      !> Number of lattice points in the y direction for ghosts
      INTEGER,                          INTENT( IN )   :: ny_gh
      !> Number of lattice points in the z direction for ghosts
      INTEGER,                          INTENT( IN )   :: nz_gh
      !> `.TRUE.` if an atmosphere should be used during the APM, to allow
      !  the real aprticles more freedom to move around and adjust;
      !  `.FALSE.` otherwise
      LOGICAL,                          INTENT( INOUT ):: use_atmosphere
      !> Name for the formatted file where the initial particle positions
      !  and the ghost positions will be printed
      CHARACTER( LEN= * ),              INTENT( INOUT ), OPTIONAL :: &
                                                            namefile_pos_id
      !> Name for the formatted file where the particle positions
      !  and the ghost positions will be printed every 15 iterations
      CHARACTER( LEN= * ),              INTENT( INOUT ), OPTIONAL :: &
                                                            namefile_pos
      !> Name for the formatted file where various quantities related
      !  to the particle distribution, the baryon number particle and the
      !  kernel estimate of the density will be printed at the end of
      !  the APM iteration
      CHARACTER( LEN= * ),              INTENT( INOUT ), OPTIONAL :: &
                                                            namefile_results


    END SUBROUTINE perform_apm


    MODULE SUBROUTINE read_sphincs_dump_print_formatted( THIS, namefile_bin, &
                                                               namefile )
    !# Reads the binary ID file printed by
    !  [[particles:compute_and_export_SPH_variables]]
    !   and prints the data stored in it to a formatted file

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles),    INTENT( IN OUT )           :: THIS
      !> Name of the binary file to be read
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile_bin
      !> Name of the formatted file to be printed
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE read_sphincs_dump_print_formatted


    MODULE SUBROUTINE print_formatted_lorene_id_particles( THIS, namefile )
    !! Prints the SPH ID to a formatted file

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles),    INTENT( IN OUT )           :: THIS
      !> Name of the formatted output file
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE print_formatted_lorene_id_particles


    MODULE SUBROUTINE read_compose_composition( THIS, namefile )
    !! Reads the \(Y_e(n_b)\) table in the CompOSE file with extension .beta

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles),    INTENT( IN OUT )           :: THIS
      !& To read the file great_eos.beta in directory compose_path/GREAT_EoS,
      !  namefile="GREAT_EoS/great_eos"
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE read_compose_composition


    MODULE SUBROUTINE compute_Ye( THIS )!, nlrf, Ye )
    !# Interpolates linearly the electron fraction \(Y_e\) at the particle
    !  densities; that is, assigns \(Y_e\) at the particle positions

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles),    INTENT( IN OUT )           :: THIS
      !DOUBLE PRECISION, DIMENSION( : ), INTENT( IN ):: nlrf
      !DOUBLE PRECISION, DIMENSION( : ), INTENT( OUT ):: Ye

    END SUBROUTINE compute_Ye


    MODULE SUBROUTINE print_summary( THIS, filename )
    !# Prints a summary of the properties of the |sph| particle
    !  distribution, optionally, to a formatted file whose name
    !  is given as the optional argument `filename`


      CLASS(particles), INTENT( IN OUT ):: THIS
      CHARACTER( LEN= * ), INTENT( INOUT ), OPTIONAL:: filename
      !! Name of the formatted file to print the summary to

    END SUBROUTINE print_summary


    MODULE SUBROUTINE destruct_particles( THIS )
    !> Finalizer (Destructor) of [[particles]] object

      !> [[particles]] object which this PROCEDURE is a member of
      TYPE(particles), INTENT( IN OUT ):: THIS

    END SUBROUTINE destruct_particles

    !-----------------!
    !--  FUNCTIONS  --!
    !-----------------!

    MODULE FUNCTION is_empty( THIS ) RESULT( answer )
    !# Returns `.TRUE` if the [[particles]] object is empty, `.FALSE` otherwise
    !  @warning experimental, not actively used in the code yet

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT( IN ):: THIS
      !> `.TRUE` if the [[particles]] object is empty, `.FALSE` otherwise
      LOGICAL:: answer

    END FUNCTION is_empty

   !MODULE SUBROUTINE write_lorene_bns_id_dump( THIS, namefile )
   !
   !    CLASS(particles),    INTENT( IN )               :: THIS
   !    CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile
   !
   !END SUBROUTINE write_lorene_bns_id_dump

    MODULE PURE FUNCTION get_npart( THIS ) RESULT( n_part )
    !! Returns [[particles:npart]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT( IN ):: THIS
      !> [[particles:npart]]
      INTEGER:: n_part

    END FUNCTION get_npart

    MODULE PURE FUNCTION get_npart_i( THIS, i_matter ) RESULT( n_part )
    !! Returns the number of particles on the object `i_matter`

      CLASS(particles), INTENT( IN ):: THIS
      !! [[particles]] object which this PROCEDURE is a member of
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object
      INTEGER:: n_part
      !! Number of particles on the object `i_matter`

    END FUNCTION get_npart_i


    MODULE PURE FUNCTION get_nuratio( THIS ) RESULT( nuratio )
    !! Returns [[particles:nuratio]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT( IN ):: THIS
      !> [[particles:nuratio]]
      DOUBLE PRECISION:: nuratio

    END FUNCTION get_nuratio


    MODULE PURE FUNCTION get_nuratio_i( THIS, i_matter ) RESULT( nuratio )
    !! Returns the baryon number ratio on the object `i_matter`

      CLASS(particles), INTENT( IN ):: THIS
      !! [[particles]] object which this PROCEDURE is a member of
      INTEGER, INTENT( IN ):: i_matter
      !! Index of the matter object
      INTEGER:: nuratio
      !! Baryon number ratio on the object `i_matter`

    END FUNCTION get_nuratio_i


    MODULE PURE FUNCTION get_pos( THIS ) RESULT( pos_u )
    !! Returns [[particles:pos]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT( IN ):: THIS
      !> [[particles:pos]]
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos_u

    END FUNCTION get_pos

    MODULE PURE FUNCTION get_vel( THIS ) RESULT( vel )
    !! Returns [[particles:v]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT( IN ):: THIS
      !> [[particles:v]]
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: vel

    END FUNCTION get_vel

    MODULE PURE FUNCTION get_nlrf( THIS ) RESULT( nlrf )
    !! Returns [[particles:nlrf]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT( IN ):: THIS
      !> [[particles:nlrf]]
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: nlrf

    END FUNCTION get_nlrf

    MODULE PURE FUNCTION get_nu( THIS ) RESULT( nu )
    !! Returns [[particles:nu]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT( IN ):: THIS
      !> [[particles:nu]]
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: nu

    END FUNCTION get_nu

    MODULE PURE FUNCTION get_u( THIS ) RESULT( u )
    !! Returns [[particles:specific_energy_parts]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT( IN ):: THIS
      !> [[particles:specific_energy_parts]]
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: u

    END FUNCTION get_u

    MODULE PURE FUNCTION get_pressure( THIS ) RESULT( pressure )
    !! Returns [[particles:pressure_parts]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT( IN ):: THIS
      !> [[particles:pressure_parts]]
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: pressure

    END FUNCTION get_pressure

    MODULE PURE FUNCTION get_pressure_cu( THIS ) RESULT( pressure_cu )
    !! Returns [[particles:pressure_parts_cu]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT( IN ):: THIS
      !> [[particles:pressure_parts_cu]]
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: pressure_cu

    END FUNCTION get_pressure_cu

    MODULE PURE FUNCTION get_theta( THIS ) RESULT( theta )
    !! Returns [[particles:theta]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT( IN ):: THIS
      !> [[particles:theta]]
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: theta

    END FUNCTION get_theta

    MODULE PURE FUNCTION get_h( THIS ) RESULT( h )
    !! Returns [[particles:h]]

      !> [[particles]] object which this PROCEDURE is a member of
      CLASS(particles), INTENT( IN ):: THIS
      !> [[particles:h]]
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: h

    END FUNCTION get_h

  END INTERFACE



  CONTAINS



  SUBROUTINE check_particle_positions( npart, pos, debug )

    !*************************************************
    !
    !# Check that the particles are not at the same
    !  positions
    !
    !  FT 1.9.2021
    !
    !*************************************************

    USE NR, ONLY: indexx

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

    IF( PRESENT(debug) .AND. debug .EQV. .TRUE. )THEN

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


  FUNCTION check_particle_position( npart, pos, pos_a ) RESULT( cnt )

    !*****************************************************
    !
    !# Return the number of times that pos_a appears
    !  in the array pos
    !
    !  FT 13.10.2021
    !
    !*****************************************************

    !USE NR,             ONLY: indexx

    IMPLICIT NONE

    INTEGER, INTENT(IN):: npart
    DOUBLE PRECISION, DIMENSION(3,npart), INTENT(IN):: pos
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: pos_a
    INTEGER:: cnt

    INTEGER:: itr, itr2, size_x!, cnt
    INTEGER, DIMENSION(npart):: x_sort, cnts
    INTEGER, DIMENSION(npart):: x_number
    INTEGER, DIMENSION(:), ALLOCATABLE:: x_number_filt

    ! Sort x coordinates of the particles
    !CALL indexx( npart, pos( 1, : ), x_sort )

    x_number= 0
    itr2= 0
    ! Find the number of times that the x coordinate of pos_a appears in pos
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( pos, pos_a, x_sort, x_number, npart ) &
    !$OMP             PRIVATE( itr )
    DO itr= 1, npart, 1

      IF( pos( 1, itr ) == pos_a( 1 ) )THEN

        !itr2= itr2 + 1
        x_number(itr)= itr

      !ELSEIF( pos( 1, x_sort(itr) ) > pos_a( 1 ) )THEN
      !
      !  EXIT

      ENDIF

    ENDDO
    !$OMP END PARALLEL DO
    x_number_filt= PACK( x_number, x_number /= 0 )
    size_x= SIZE(x_number_filt)

    cnts= 0
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( pos, pos_a, x_sort, x_number_filt, size_x, cnts )&
    !$OMP             PRIVATE( itr )
    DO itr= 1, size_x, 1

      ! If they have the same y
      IF( pos( 2, x_number_filt(itr) ) == pos_a( 2 ) )THEN

        ! If they have the same z
        IF( pos( 3, x_number_filt(itr) ) == pos_a( 3 ) )THEN

          cnts(itr)= cnts(itr) + 1

        ENDIF
      ENDIF

    ENDDO
    !$OMP END PARALLEL DO

    cnt= SUM( cnts )

  END FUNCTION check_particle_position


  FUNCTION find_h_backup( a, npart, pos, ndes ) RESULT( h )

    !**************************************************************
    !
    !# Backup method to find the smoothing length via brute force
    !  if the optimized method gives 0
    !
    !  FT 24.11.2021
    !
    !**************************************************************

    USE NR, ONLY: select

    IMPLICIT NONE

    INTEGER,          INTENT(IN):: a, npart, ndes
    DOUBLE PRECISION, DIMENSION(3,npart), INTENT(IN):: pos

    DOUBLE PRECISION:: h

    INTEGER:: b
    DOUBLE PRECISION, DIMENSION(npart):: dist2
    DOUBLE PRECISION, DIMENSION(3):: dist

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( pos, a, npart, dist2 ) &
    !$OMP             PRIVATE( b, dist )
    DO b= 1, npart, 1

      IF( a /= b )THEN

        dist(:)= pos(:,b) - pos(:,a)
        dist2(b)= DOT_PRODUCT(dist,dist)

      ENDIF

    ENDDO
    !$OMP END PARALLEL DO

    h= 0.5D0*SQRT( select(ndes, npart, dist2) )

  END FUNCTION find_h_backup


END MODULE particles_id
