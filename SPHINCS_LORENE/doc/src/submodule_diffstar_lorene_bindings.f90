! File:         submodule_bns_bindings.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (diffstar_lorene) diffstar_lorene_bindings

        !***********************************************
        !                                              *
        !# Define SUBROUTINES that create and destroy  *
        !  a Lorene::|etdiffrot| object from LORENE,     *
        !  and call its member methods                 *
        !  that extracts the LORENE initial data       *
        !  from the spectral data saved in the 'resu'  *
        !  binary files produced by LORENE             *
        !                                              *
        !  For more information on how to link C
        !  functions to Fortran,
        !  [see here](https://community.intel.com/t5/Intel-Fortran-Compiler/Calling-C-cpp-bin_nsects-from-a-Fortran-subroutine/td-p/1110556){:target="_blank"}.
        !                                              *
        !  FT 24.10.2021                               *
        !                                              *
        !***********************************************

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, &
                                           C_CHAR, C_NULL_CHAR, &
                                           C_PTR, C_NULL_PTR, C_ASSOCIATED

    IMPLICIT NONE


    !---------------------------------------------------------------!
    !--  Interfaces to the methods of |lorene|'s class |etdiffrot|  --!
    !---------------------------------------------------------------!


    INTERFACE


      FUNCTION construct_etdiffrot( c_resu_file ) RESULT( optr ) &
        BIND(C, NAME= "construct_et_diffrot")

        !***********************************************
        !
        !# Interface to the |lorene| method of class
        !  |etdiffrot| with the same name, that constructs
        !  the |lorene| |etdiffrot| object
        !
        !  FT 24.10.2021
        !
        !***********************************************

        IMPORT :: C_PTR, C_CHAR

        IMPLICIT NONE

        !& C string of the name of the |lorene| binary file storing the spectral
        !  DRS ID
        CHARACTER(KIND= C_CHAR), DIMENSION(*), INTENT(IN), OPTIONAL :: &
                                                                c_resu_file
        !> C pointer pointing to the constructed |lorene| |etdiffrot| object
        TYPE(C_PTR) :: optr

      END FUNCTION construct_etdiffrot


      SUBROUTINE get_diffstar_full( optr, &
                                    x, y, z, &
                                    lapse, &
                                    shift_x, shift_y, shift_z, &
                                    g_diag, &
                                    k_xx, k_xy, k_xz, &
                                    k_yy, k_yz, k_zz, &
                                    baryon_density, &
                                    energy_density, &
                                    specific_energy, &
                                    v_euler_x, v_euler_y, v_euler_z ) &
        BIND(C, NAME= "get_rotdiff_id")

        !*************************************************
        !
        !# Interface to the |lorene| method of class
        !  |etdiffrot| with the same name, that reads the full
        !  |lorene| ID at the specified point.
        !  That is, imports the metric fields, the
        !  components of the extrinsic curvature [c/km],
        !  and the hydro fields.
        !
        !  - shift vector [c]
        !  - baryon mass density [kg m^{-3}]
        !  - energy density [kg c^2 m^{-3}]
        !  - pressure [kg c^2 m^{-3}]
        !  - specific internal energy [c^2]
        !  - fluid 3-velocity with respect to the
        !    Eulerian observer [c]
        !
        !  FT 24.10.2021
        !
        !*************************************************

        IMPORT :: C_DOUBLE, C_PTR

        IMPLICIT NONE

        !> C pointer pointing to a |lorene| |etdiffrot| object
        TYPE(C_PTR),    INTENT(IN), VALUE :: optr
        !> \(x\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN), VALUE :: x
        !> \(y\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN), VALUE :: y
        !> \(z\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN), VALUE :: z
        REAL(C_DOUBLE), INTENT(OUT)       :: lapse
        REAL(C_DOUBLE), INTENT(OUT)       :: shift_x
        REAL(C_DOUBLE), INTENT(OUT)       :: shift_y
        REAL(C_DOUBLE), INTENT(OUT)       :: shift_z
        REAL(C_DOUBLE), INTENT(OUT)       :: g_diag
        REAL(C_DOUBLE), INTENT(OUT)       :: k_xx
        REAL(C_DOUBLE), INTENT(OUT)       :: k_xy
        REAL(C_DOUBLE), INTENT(OUT)       :: k_xz
        REAL(C_DOUBLE), INTENT(OUT)       :: k_yy
        REAL(C_DOUBLE), INTENT(OUT)       :: k_yz
        REAL(C_DOUBLE), INTENT(OUT)       :: k_zz
        REAL(C_DOUBLE), INTENT(OUT)       :: baryon_density
        REAL(C_DOUBLE), INTENT(OUT)       :: energy_density
        REAL(C_DOUBLE), INTENT(OUT)       :: specific_energy
        REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_x
        REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_y
        REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_z

      END SUBROUTINE get_diffstar_full


      SUBROUTINE get_diffstar_spacetime( optr, &
                                         x, y, z, &
                                         lapse, &
                                         shift_x, shift_y, shift_z, &
                                         g_diag, &
                                         k_xx, k_xy, k_xz, &
                                         k_yy, k_yz, k_zz ) &
        BIND(C, NAME= "get_rotdiff_spacetime")

        !*************************************************
        !
        !# Interface to the |lorene| method of class
        !  |etdiffrot| with the same name, that reads the
        !  metric fields and the components
        !  of the extrinsic curvature [c/km] from |lorene|,
        !  at the specified point
        !
        !  FT 24.10.2021
        !
        !*************************************************

        IMPORT :: C_DOUBLE, C_PTR

        IMPLICIT NONE

        !> C pointer pointing to a |lorene| |etdiffrot| object
        TYPE(C_PTR),    INTENT(IN), VALUE :: optr
        !> \(x\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN), VALUE :: x
        !> \(y\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN), VALUE :: y
        !> \(z\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN), VALUE :: z
        REAL(C_DOUBLE), INTENT(OUT)       :: lapse
        REAL(C_DOUBLE), INTENT(OUT)       :: shift_x
        REAL(C_DOUBLE), INTENT(OUT)       :: shift_y
        REAL(C_DOUBLE), INTENT(OUT)       :: shift_z
        REAL(C_DOUBLE), INTENT(OUT)       :: g_diag
        REAL(C_DOUBLE), INTENT(OUT)       :: k_xx
        REAL(C_DOUBLE), INTENT(OUT)       :: k_xy
        REAL(C_DOUBLE), INTENT(OUT)       :: k_xz
        REAL(C_DOUBLE), INTENT(OUT)       :: k_yy
        REAL(C_DOUBLE), INTENT(OUT)       :: k_yz
        REAL(C_DOUBLE), INTENT(OUT)       :: k_zz

      END SUBROUTINE get_diffstar_spacetime


      SUBROUTINE get_diffstar_particles( optr, &
                                         x, y, z, &
                                         lapse, &
                                         shift_x, shift_y, shift_z, &
                                         g_diag, &
                                         baryon_density, &
                                         energy_density, &
                                         specific_energy, &
                                         pressure, &
                                         v_euler_x, v_euler_y, v_euler_z ) &
        BIND(C, NAME= "get_rotdiff_particles")

        !**********************************************
        !
        !# Interface to the |lorene| method of class
        !  |etdiffrot| with the same name, that reads the
        !  hydro fields and the metric fields *
        !  from |lorene|, at the specified point
        !
        !  - shift vector [c]
        !  - baryon mass density [kg m^{-3}]
        !  - energy density [kg c^2 m^{-3}]
        !  - pressure [kg c^2 m^{-3}]
        !  - specific internal energy [c^2]
        !  - fluid 3-velocity with respect to the
        !    Eulerian observer [c]
        !
        !  FT 24.10.2021
        !
        !**********************************************

        IMPORT :: C_DOUBLE, C_PTR

        IMPLICIT NONE

        !> C pointer pointing to a |lorene| |etdiffrot| object
        TYPE(C_PTR),    INTENT(IN), VALUE :: optr
        !> \(x\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN), VALUE :: x
        !> \(y\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN), VALUE :: y
        !> \(z\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN), VALUE :: z
        REAL(C_DOUBLE), INTENT(OUT)       :: lapse
        REAL(C_DOUBLE), INTENT(OUT)       :: shift_x
        REAL(C_DOUBLE), INTENT(OUT)       :: shift_y
        REAL(C_DOUBLE), INTENT(OUT)       :: shift_z
        REAL(C_DOUBLE), INTENT(OUT)       :: g_diag
        REAL(C_DOUBLE), INTENT(OUT)       :: baryon_density
        REAL(C_DOUBLE), INTENT(OUT)       :: energy_density
        REAL(C_DOUBLE), INTENT(OUT)       :: specific_energy
        REAL(C_DOUBLE), INTENT(OUT)       :: pressure
        REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_x
        REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_y
        REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_z

      END SUBROUTINE get_diffstar_particles


      SUBROUTINE get_diffstar_mass_b( optr, &
                                      x, y, z, &
                                      g_diag, &
                                      baryon_density, &
                                      gamma_euler ) &
        BIND(C, NAME= "get_rotdiff_mass_b")

        !************************************************
        !
        !# Interface to the |lorene| method of class
        !  |etdiffrot| with the same name, that reads the
        !  hydro fields and the metric fields
        !  from |lorene|, at the specified point,
        !  needed to compute the baryon mass.
        !
        !  - shift vector [c]
        !  - baryon mass density [kg m^{-3}]
        !  - fluid 3-velocity with respect to the
        !    Eulerian observer [c]
        !
        !  FT 24.10.2021
        !
        !************************************************

        IMPORT :: C_DOUBLE, C_PTR

        IMPLICIT NONE

        !> C pointer pointing to a |lorene| |etdiffrot| object
        TYPE(C_PTR),    INTENT(IN), VALUE :: optr
        !> \(x\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN), VALUE :: x
        !> \(y\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN), VALUE :: y
        !> \(z\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN), VALUE :: z
        !> \(g_{xx}=g_{yy}=g_{zz}\) at \(x,y,z\)
        REAL(C_DOUBLE), INTENT(OUT)       :: g_diag
        !> Baryon mass density at \(x,y,z\)
        REAL(C_DOUBLE), INTENT(OUT)       :: baryon_density
        !& Relative Lorentz factor between the 4-velocity of the fluid
        !  wrt the Eulerian observer and the 4-velocity of the Eulerian observer
        !  at \(x,y,z\)
        REAL(C_DOUBLE), INTENT(OUT)       :: gamma_euler

      END SUBROUTINE get_diffstar_mass_b


      SUBROUTINE get_diffstar_hydro( optr, &
                                     x, y, z, &
                                     baryon_density, &
                                     energy_density, &
                                     specific_energy, &
                                     pressure, &
                                     v_euler_x, v_euler_y, v_euler_z ) &
        BIND(C, NAME= "get_rotdiff_hydro")

        !***********************************************
        !
        !# Interface to the |lorene| method of class
        !  |etdiffrot| with the same name, that reads the
        !  hydro fields from |lorene|, at the
        !  specified point
        !
        !  - baryon mass density [kg m^{-3}]
        !  - energy density [kg c^2 m^{-3}]
        !  - pressure [kg c^2 m^{-3}]
        !  - specific internal energy [c^2]
        !  - fluid 3-velocity with respect to the
        !    Eulerian observer [c]
        !
        !  FT 24.10.2021
        !
        !***********************************************

        IMPORT :: C_DOUBLE, C_PTR

        IMPLICIT NONE

        !> C pointer pointing to a |lorene| |etdiffrot| object
        TYPE(C_PTR),    INTENT(IN), VALUE :: optr
        !> \(x\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN), VALUE :: x
        !> \(y\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN), VALUE :: y
        !> \(z\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN), VALUE :: z
        REAL(C_DOUBLE), INTENT(OUT)       :: baryon_density
        REAL(C_DOUBLE), INTENT(OUT)       :: energy_density
        REAL(C_DOUBLE), INTENT(OUT)       :: specific_energy
        REAL(C_DOUBLE), INTENT(OUT)       :: pressure
        REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_x
        REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_y
        REAL(C_DOUBLE), INTENT(OUT)       :: v_euler_z

      END SUBROUTINE get_diffstar_hydro


      FUNCTION get_lorene_mass_density( optr, x, y, z ) RESULT( res ) &
        BIND(C, NAME= "get_mass_density")

        !********************************************
        !
        !# Interface to the |lorene| method of class
        !  |etdiffrot| with the same name, that returns
        !  the baryon mass density \([\mathrm{kg}\,
        !  \mathrm{m}^{-3}]\) from |lorene|,
        !  at the specified point
        !
        !  FT 24.10.2021
        !
        !********************************************

        IMPORT :: C_DOUBLE, C_PTR

        IMPLICIT NONE

        !> C pointer pointing to a |lorene| |etdiffrot| object
        TYPE(C_PTR),    INTENT(IN),  VALUE :: optr
        !> \(x\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN),  VALUE :: x
        !> \(y\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN),  VALUE :: y
        !> \(z\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN),  VALUE :: z
        !& Baryon mass density \([\mathrm{kg}\, \mathrm{m}^{-3}]\) at the desired
        !  point \((x,y,z)\)
        REAL(C_DOUBLE) :: res

      END FUNCTION get_lorene_mass_density


      FUNCTION get_lorene_spatial_metric( optr, x, y, z ) RESULT( res ) &
        BIND(C, NAME= "get_g_diag")

        !************************************************
        !
        !# Interface to the |lorene| method of class
        !  |etdiffrot| with the same name, that returns the
        !  diagonal components of the metric,
        !  all equal to the |lorene| conformal factor to
        !  the 4th power.
        !
        !  FT 24.10.2021
        !
        !************************************************

        IMPORT :: C_DOUBLE, C_PTR

        IMPLICIT NONE

        !> C pointer pointing to a |lorene| |etdiffrot| object
        TYPE(C_PTR),    INTENT(IN),  VALUE :: optr
        !> \(x\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN),  VALUE :: x
        !> \(y\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN),  VALUE :: y
        !> \(z\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN),  VALUE :: z
        !& Spatial metric component
        !  \(g_{xx}=g_{yy}=g_{zz}\) at the point \((x,y,z)\)
        REAL(C_DOUBLE) :: res

      END FUNCTION get_lorene_spatial_metric


      FUNCTION negative_hydro( optr, x, y, z ) RESULT( res ) &
        BIND(C, NAME= "is_hydro_negative")

        !************************************************
        !
        !# Interface to the |lorene| method of class
        !  |etdiffrot| with the same name, that returns 1
        !  if the energy density is nonpositive,
        !  or if the specific energy is nonpositive,
        !  or if the pressure is nonpositive,
        !  at the specified point; it returns 0 otherwise
        !
        !  FT 24.10.2021
        !
        !************************************************

        IMPORT :: C_INT, C_DOUBLE, C_PTR

        IMPLICIT NONE

        !> C pointer pointing to a |lorene| |etdiffrot| object
        TYPE(C_PTR),    INTENT(IN),  VALUE :: optr
        !> \(x\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN),  VALUE :: x
        !> \(y\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN),  VALUE :: y
        !> \(z\) coordinate of the desired point
        REAL(C_DOUBLE), INTENT(IN),  VALUE :: z
        !& 1 if the energy density or the specific energy or the pressure
        !  are negative, 0 otherwise
        INTEGER(C_INT) :: res

      END FUNCTION negative_hydro


      SUBROUTINE get_diffstar_params( optr,                           &
                                      omega_c,                        &
                                      mass,                           &
                                      mass_grav,                      &
                                      angular_momentum,               &
                                      tsw,                            &
                                      grv2,                           &
                                      grv3,                           &
                                      r_circ,                         &
                                      surface_area,                   &
                                      r_mean,                         &
                                      r_eq,                           &
                                      r_eq_pi2,                       &
                                      r_eq_pi,                        &
                                      r_eq_3pi2,                      &
                                      r_eq_pole,                      &
                                      r_ratio,                        &
                                      r_isco,                         &
                                      f_isco,                         &
                                      specific_energy_isco,           &
                                      specific_angular_momentum_isco, &
                                      area_radius,                    &
                                      ent_center,                     &
                                      nbar_center,                    &
                                      rho_center,                     &
                                      energy_density_center,          &
                                      specific_energy_center,         &
                                      pressure_center,                &
                                      redshift_eqf,                   &
                                      redshift_eqb,                   &
                                      redshift_pole,                  &
                                      eos,                            &
                                      eos_id,                         &
                                      gamma,                          &
                                      kappa,                          &
                                      npeos,                          &
                                      gamma0,                         &
                                      gamma1,                         &
                                      gamma2,                         &
                                      gamma3,                         &
                                      kappa0,                         &
                                      kappa1,                         &
                                      kappa2,                         &
                                      kappa3,                         &
                                      logP1,                          &
                                      logRho0,                        &
                                      logRho1 )                       &
        BIND(C, NAME= "get_rotdiff_params")

        !**********************************************
        !
        !# Interface to the |lorene| method of class
        !  |etdiffrot| with the same name, that stores
        !  the physical parameters of the binary
        !  system from |lorene| in the desired variables
        !
        !  FT 24.10.2021
        !
        !**********************************************

        IMPORT :: C_INT, C_DOUBLE, C_PTR, C_CHAR

        IMPLICIT NONE

        !> C pointer pointing to a |lorene| |etdiffrot| object
        TYPE(C_PTR),    INTENT(IN), VALUE :: optr
        REAL(C_DOUBLE), INTENT(OUT)       :: omega_c
        REAL(C_DOUBLE), INTENT(OUT)       :: mass
        REAL(C_DOUBLE), INTENT(OUT)       :: mass_grav
        REAL(C_DOUBLE), INTENT(OUT)       :: angular_momentum
        REAL(C_DOUBLE), INTENT(OUT)       :: tsw
        REAL(C_DOUBLE), INTENT(OUT)       :: grv2
        REAL(C_DOUBLE), INTENT(OUT)       :: grv3
        REAL(C_DOUBLE), INTENT(OUT)       :: r_circ
        REAL(C_DOUBLE), INTENT(OUT)       :: surface_area
        REAL(C_DOUBLE), INTENT(OUT)       :: r_mean
        REAL(C_DOUBLE), INTENT(OUT)       :: r_eq
        REAL(C_DOUBLE), INTENT(OUT)       :: r_eq_pi2
        REAL(C_DOUBLE), INTENT(OUT)       :: r_eq_pi
        REAL(C_DOUBLE), INTENT(OUT)       :: r_eq_3pi2
        REAL(C_DOUBLE), INTENT(OUT)       :: r_eq_pole
        REAL(C_DOUBLE), INTENT(OUT)       :: r_ratio
        REAL(C_DOUBLE), INTENT(OUT)       :: r_isco
        REAL(C_DOUBLE), INTENT(OUT)       :: f_isco
        REAL(C_DOUBLE), INTENT(OUT)       :: specific_energy_isco
        REAL(C_DOUBLE), INTENT(OUT)       :: specific_angular_momentum_isco
        REAL(C_DOUBLE), INTENT(OUT)       :: area_radius
        REAL(C_DOUBLE), INTENT(OUT)       :: ent_center
        REAL(C_DOUBLE), INTENT(OUT)       :: nbar_center
        REAL(C_DOUBLE), INTENT(OUT)       :: rho_center
        REAL(C_DOUBLE), INTENT(OUT)       :: energy_density_center
        REAL(C_DOUBLE), INTENT(OUT)       :: specific_energy_center
        REAL(C_DOUBLE), INTENT(OUT)       :: pressure_center
        REAL(C_DOUBLE), INTENT(OUT)       :: redshift_eqf
        REAL(C_DOUBLE), INTENT(OUT)       :: redshift_eqb
        REAL(C_DOUBLE), INTENT(OUT)       :: redshift_pole
        CHARACTER(KIND=C_CHAR), DIMENSION(100), INTENT(OUT):: eos
        INTEGER(C_INT)                    :: eos_id
        REAL(C_DOUBLE), INTENT(OUT)       :: gamma
        REAL(C_DOUBLE), INTENT(OUT)       :: kappa
        INTEGER(C_INT)                    :: npeos
        REAL(C_DOUBLE), INTENT(OUT)       :: gamma0
        REAL(C_DOUBLE), INTENT(OUT)       :: gamma1
        REAL(C_DOUBLE), INTENT(OUT)       :: gamma2
        REAL(C_DOUBLE), INTENT(OUT)       :: gamma3
        REAL(C_DOUBLE), INTENT(OUT)       :: kappa0
        REAL(C_DOUBLE), INTENT(OUT)       :: kappa1
        REAL(C_DOUBLE), INTENT(OUT)       :: kappa2
        REAL(C_DOUBLE), INTENT(OUT)       :: kappa3
        REAL(C_DOUBLE), INTENT(OUT)       :: logP1
        REAL(C_DOUBLE), INTENT(OUT)       :: logRho0
        REAL(C_DOUBLE), INTENT(OUT)       :: logRho1

      END SUBROUTINE get_diffstar_params


      SUBROUTINE destruct_etdiffrot( optr ) &
        BIND(C, NAME= "destruct_et_diffrot")

        !**********************************************
        !
        !# Interface to the |lorene| method of class
        !  |etdiffrot| with the same name, that destructs
        !  the |lorene| |etdiffrot| object
        !
        ! FT 24.10.2021
        !
        !**********************************************

        IMPORT :: C_PTR

        IMPLICIT NONE

        TYPE(C_PTR), INTENT(IN), VALUE :: optr
        !! C pointer pointing to the |lorene| |etdiffrot| object to destruct

      END SUBROUTINE destruct_etdiffrot


    END INTERFACE


END SUBMODULE diffstar_lorene_bindings
