! File:         module_bssn_formulation.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

MODULE bssn_formulation

  !********************************************
  !
  !# This module contains the definition of
  !  the EXTENDED TYPE [[bssn]], representing
  !  the |id| for the |bssn| formulation of
  !  the |ee|
  !
  !********************************************


  USE id_base,                  ONLY: idbase
  USE standard_tpo_formulation, ONLY: tpo
  USE sph_particles,            ONLY: particles
  USE mesh_refinement,          ONLY: grid_function_scalar, grid_function
  USE timing,                   ONLY: timer
  USE utility,                  ONLY: ios, err_msg, perc, creturn, run_id, &
                                      test_status, compute_g4, &
                                      determinant_sym4x4_grid, show_progress


  IMPLICIT NONE


  !********************************************************
  !                                                       *
  !              Definition of TYPE bssn                  *
  !                                                       *
  ! This class extends the ABSTRACT TYPE tpo by           *
  ! implementing its deferred methods such that the BSSN  *
  ! variables are computed on the grid for the ID,        *
  ! stored, exported to a binary file for evolution and   *
  ! to a formatted file. The BSSN constraints can also    *
  ! be computed in different ways, analyzed, and exported *
  ! in different ways.                                    *
  !                                                       *
  !********************************************************

  TYPE, EXTENDS(tpo):: bssn
  !# TYPE representing the |id| for the |bssn| formulation
  !  of the Einstein equations


    INTEGER:: call_flag= 0
    !# Flag set to a value different than 0 if the SUBROUTINE
    !  compute_and_export_bssn_variables is called

    !
    !-- Arrays storing the BSSN variables for the LORENE ID on the grid
    !

    TYPE(grid_function):: Gamma_u
    !! Conformal connection \(\bar{\Gamma} ^i_{jk}\)

    TYPE(grid_function_scalar):: phi
    !! Conformal factor \(\phi \)

    TYPE(grid_function_scalar):: trK
    !! Trace of extrinsic curvature \(K \)

    TYPE(grid_function):: A_BSSN3_ll
    !! Conformal traceless extrinsic curvature \(A_{ij} \)

    TYPE(grid_function):: g_BSSN3_ll
    !! Conformal spatial metric \(\gamma_{ij} \)

    !
    !-- Connection constraints and its l2 norm and loo norm
    !

    TYPE(grid_function):: GC
    !! Connection constraint computed with the |id| on the mesh
    TYPE(grid_function):: GC_parts
    !# Connection constraint computed with the |bssn| |id| on the mesh, and
    !  the hydrodynamical |id| mapped from the particles to the mesh

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: GC_l2
    !# \(\ell_2\) norm of the connection constraint computed
    !  with the |id| on the mesh
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: GC_parts_l2
    !# \(\ell_2\) norm of the connection constraint computed with the |bssn|
    !  |id| on the mesh, and the hydrodynamical |id| mapped from the particles
    !  to the mesh
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: GC_loo
    !# \(\ell_\infty\) norm, i.e., supremum of the absolute value, of the
    !  connection constraint computed with the |id| on the mesh
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: GC_parts_loo
    !# \(\ell_\infty\) norm, i.e., supremum of the absolute value, of the
    !  connection constraint computed with the |bssn| |id| on the mesh,
    !  and the hydrodynamical |id| mapped from the particles to the mesh

    LOGICAL, PUBLIC:: export_bin
    !# `.TRUE.` if the binary files for SPHINCS_BSSN are to be exported,
    !  `.FALSE.` otherwise
    LOGICAL, PUBLIC:: export_form_xy
    !# `.TRUE.` if the |id| in the formatted files is to be on the xy plane
    !  only, `.FALSE.` otherwise
    LOGICAL, PUBLIC:: export_form_x
    !# `.TRUE.` if the |id| in the formatted files is to be on the x axis
    !  only, `.FALSE.` otherwise

    TYPE(timer):: bssn_computer_timer
    !# Timer that times how long it takes to compute the |bssn| variables on
    !  the refined mesh


    CONTAINS


    !-------------------!
    !--  SUBROUTINES  --!
    !-------------------!


    PROCEDURE :: define_allocate_fields => allocate_bssn_fields
    !! Allocates memory for the [[bssn]] member grid functions

    PROCEDURE :: deallocate_fields => deallocate_bssn_fields
    !! Deallocates memory for the [[bssn]] member arrays

    PROCEDURE :: compute_and_export_tpo_variables &
                    => compute_and_export_bssn_variables
    !# Computes the |bssn| variables at the particle positions, and optionally
    !  prints them to a binary file to be read by \(\texttt{SPHINCS_BSSN}\)
    !  and \(\texttt{splash}\), and to a formatted file to be read by
    !  \(\texttt{gnuplot}\), by calling
    !  [[bssn:print_formatted_lorene_id_tpo_variables]]

    PROCEDURE, PUBLIC :: read_bssn_dump_print_formatted
    !# Reads the binary |id| file printed by
    !  [[bssn:compute_and_export_tpo_variables]]

    PROCEDURE :: print_formatted_lorene_id_tpo_variables &
                    => print_formatted_lorene_id_bssn_variables
    !! Prints the |bssn| |id| to a formatted file

    PROCEDURE :: compute_and_export_tpo_constraints_grid &
                    => compute_and_export_bssn_constraints_grid
    !# Computes the |bssn| constraints using the full |id| on the refined mesh,
    !  prints a summary with the statistics for the constraints. Optionally,
    !  prints the constraints to a formatted file to be read by
    !  \(\texttt{gnuplot}\), and analyze the constraints by calling
    !  [[tpo:analyze_constraint]]

    PROCEDURE :: compute_and_export_tpo_constraints_particles &
                    => compute_and_export_bssn_constraints_particles
    !# Computes the |bssn| constraints using the |bssn| |id| on the refined
    !  mesh and the hydrodynamical |id| mapped from the particles to the mesh,
    !  prints a summary with the statistics for the constraints. Optionally,
    !  prints the constraints to a formatted file to be read by
    !  \(\texttt{gnuplot}\), and analyze the constraints by calling
    !  [[tpo:analyze_constraint]]

    PROCEDURE :: destruct_bssn
    !# Destructor for the EXTENDED TYPE bssn, not ABSTRACT TYPE tpo

    FINAL     :: destructor
    !# Destructor; finalizes members from both TYPES tpo and bssn,
    !  by calling [[tpo:destruct_tpo]] and
    !  [[bssn:destruct_bssn]]

  END TYPE bssn

  !
  !-- Interface of the TYPE bssn
  !-- (i.e., declaration of the overloaded constructor)
  !
  INTERFACE bssn

    MODULE PROCEDURE:: construct_bssn
    !# Constructs the bssn object from the number of grid points
    !  along each axis

  END INTERFACE bssn

  !
  !-- Interface of the constructor of TYPE bssn
  !-- Its implementation is in submodule_BSSN_id_constructor.f90
  !
  INTERFACE

    MODULE FUNCTION construct_bssn( id, dx, dy, dz ) RESULT ( bssnid )
    !# Constructs the [[bssn]] object from the number of grid points
    !  along each axis

      CLASS(idbase), INTENT( INOUT ):: id
      !! [[idbase]] object to use to construct the [[bssn]] object
      TYPE(bssn)              :: bssnid
      !! [[bssn]] object to be constructed
      DOUBLE PRECISION, OPTIONAL :: dx, dy, dz
      !! Mesh spacings @todo for which refinement level?

    END FUNCTION construct_bssn

  END INTERFACE

  !
  !-- Interfaces of the methods of TYPE [[bssn]]
  !-- Their implementations are in submodule_BSSN_id_methods.f90
  !
  INTERFACE


    MODULE SUBROUTINE allocate_bssn_fields( THIS )
    !! Interface to [[bssn:define_allocate_fields]]

      CLASS(bssn), INTENT( IN OUT ):: THIS
      !! [[bssn]] object to which this PROCEDURE is bound

    END SUBROUTINE allocate_bssn_fields


    MODULE SUBROUTINE compute_and_export_bssn_variables( THIS, namefile )
    !! Interface to [[bssn:compute_and_export_tpo_variables]]

      CLASS(bssn),      INTENT( IN OUT )           :: THIS
      !! [[bssn]] object to which this PROCEDURE is bound
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE compute_and_export_bssn_variables


    MODULE SUBROUTINE read_bssn_dump_print_formatted( THIS, namefile_bin, &
                                                            namefile )
    !! Interface to [[bssn:read_bssn_dump_print_formatted]]

      CLASS(bssn),      INTENT( IN OUT )           :: THIS
      !! [[bssn]] object to which this PROCEDURE is bound
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile_bin
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE read_bssn_dump_print_formatted


    MODULE SUBROUTINE print_formatted_lorene_id_bssn_variables( THIS, &
                                                                namefile )
    !! Interface to [[bssn:print_formatted_lorene_id_tpo_variables]]

      CLASS(bssn),      INTENT( IN OUT )           :: THIS
      !! [[bssn]] object to which this PROCEDURE is bound
      CHARACTER( LEN= * ), INTENT( IN OUT ), OPTIONAL :: namefile

    END SUBROUTINE print_formatted_lorene_id_bssn_variables


    MODULE SUBROUTINE compute_and_export_bssn_constraints_grid( THIS, &
                                                           id, &
                                                           namefile, &
                                                           name_logfile )
    !! Interface to [[bssn:compute_and_export_tpo_constraints_grid]]

      CLASS(bssn),      INTENT( IN OUT ):: THIS
      !! [[bssn]] object to which this PROCEDURE is bound
      CLASS(idbase),      INTENT( IN OUT ):: id
      !! [[idbase]] object used to read the hydrodynamical |id| to the mesh
      CHARACTER( LEN= * ), INTENT( IN OUT ):: namefile
      CHARACTER( LEN= * ), INTENT( IN OUT ):: name_logfile

    END SUBROUTINE compute_and_export_bssn_constraints_grid


    MODULE SUBROUTINE compute_and_export_bssn_constraints_particles( THIS, &
                                                           parts_obj, &
                                                           namefile, &
                                                           name_logfile )
    !! Interface to [[bssn:compute_and_export_tpo_constraints_particles]]

      CLASS(bssn),      INTENT( IN OUT ):: THIS
      !! [[bssn]] object to which this PROCEDURE is bound
      CLASS(particles),    INTENT( IN OUT ):: parts_obj
      !! [[particles]] object used to map the hydrodynamical |id| to the mesh
      CHARACTER( LEN= * ), INTENT( IN OUT ):: namefile
      CHARACTER( LEN= * ), INTENT( IN OUT ):: name_logfile

    END SUBROUTINE compute_and_export_bssn_constraints_particles


    MODULE SUBROUTINE deallocate_bssn_fields( THIS )
    !! Interface to [[bssn:deallocate_fields]]

      CLASS(bssn), INTENT( IN OUT ):: THIS
      !! [[bssn]] object to which this PROCEDURE is bound

    END SUBROUTINE deallocate_bssn_fields


    MODULE SUBROUTINE destruct_bssn( THIS )
    !! Interface to [[bssn:destruct_bssn]]

      CLASS(bssn), INTENT( IN OUT ):: THIS
      !! [[bssn]] object to which this PROCEDURE is bound

    END SUBROUTINE destruct_bssn


    MODULE SUBROUTINE destructor( THIS )
    !! Interface to [[bssn:destructor]]

      TYPE(bssn), INTENT( IN OUT ):: THIS
      !! [[bssn]] object to which this PROCEDURE is bound, to be destructed

    END SUBROUTINE destructor


  END INTERFACE


END MODULE bssn_formulation
