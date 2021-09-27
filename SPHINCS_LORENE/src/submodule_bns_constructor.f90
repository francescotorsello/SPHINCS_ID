! File:         submodule_bns_constructor.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (bns_id) bns_constructor

  !*********************************************************
  !
  !# Implementation of the constructor and
  !  destructor of TYPE [[bns]], and of the [[bns]]-member
  !  PROCEDURES that call the C-bound PROCEDURES
  !  constructig and destructing the LORENE
  !  Bin_NS object
  !
  !  FT 23.10.2020
  !
  !*********************************************************


  IMPLICIT NONE


  CONTAINS


  !
  !-- Implementation of the constructor of the bns object
  !
  MODULE PROCEDURE construct_bns

    !****************************************************
    !
    !# Constructs an object of TYPE [[bns]]
    !
    !  FT
    !
    !****************************************************

    IMPLICIT NONE

    INTEGER, SAVE:: bns_counter= 1

    !DOUBLE PRECISION:: tmp

    bns_obj% binary_construction_timer= timer( "binary_construction_timer" )

    ! Construct LORENE Bin_NS object
    IF( PRESENT( resu_file ) )THEN
        CALL bns_obj% construct_binary( resu_file )
    ELSE
        CALL bns_obj% construct_binary()
    ENDIF

    ! Import the parameters of the binary system
    CALL import_id_params( bns_obj )

    ! Assign a unique identifier to the bns object
    bns_obj% bns_identifier= bns_counter
    bns_counter= bns_counter + 1

    ! Do not use the geodesic gauge by default
    bns_obj% one_lapse = .FALSE.
    bns_obj% zero_shift= .FALSE.

    !PRINT *, "End of bns constructor."
    !PRINT *

  END PROCEDURE construct_bns


  !
  !-- Implementation of the destructor of the bns object
  !
  MODULE PROCEDURE destruct_bns

    !***********************************************
    !
    !# Destructs an object of TYPE [[bns]]
    !
    !  FT
    !
    !***********************************************

    IMPLICIT NONE

    !PRINT *, "Inside destructor of bns."
    !PRINT *

    ! Deallocate memory
    CALL THIS% deallocate_lorene_id_memory()

  END PROCEDURE destruct_bns


  MODULE PROCEDURE construct_binary

    !***********************************************
    !
    !# Construct the LORENE Bin_NS object
    !
    !  FT
    !
    !***********************************************

    IMPLICIT NONE

    CHARACTER(KIND= C_CHAR, LEN= 7):: default_case
    LOGICAL:: exist

    !PRINT *, "** Executing the construct_binary subroutine..."

#ifdef __INTEL_COMPILER

    IF ( C_ASSOCIATED( THIS% bns_ptr ) ) THEN

      CALL destruct_bin_ns( THIS% bns_ptr )

    ENDIF

#endif

    !
    !-- If the name of the LORENE binary file resu_file is given as argument to
    !-- construct_binary, use it. Otherwise, give the string "read_it"
    !-- to construct_bin_ns as argument, which makes LORENE read the name of
    !-- the file from the parameter file read_bin_ns.par
    !
    IF( PRESENT( resu_file ) )THEN

      INQUIRE( FILE= resu_file, EXIST= exist )

      IF( exist )THEN

        CALL THIS% binary_construction_timer% start_timer()
        THIS% bns_ptr = construct_bin_ns( resu_file//C_NULL_CHAR )
        CALL THIS% binary_construction_timer% stop_timer()

      ELSE

        PRINT *, "** ERROR in SUBROUTINE construct_binary: file ", &
                 resu_file, " cannot be found!"
        PRINT *
        STOP

      ENDIF

    ELSE

      default_case= "read_it"
      CALL THIS% binary_construction_timer% start_timer()
      THIS% bns_ptr = construct_bin_ns( default_case//C_NULL_CHAR )
      CALL THIS% binary_construction_timer% stop_timer()

    ENDIF

    !PRINT *, "** Subroutine construct_binary executed."
    !PRINT *

  END PROCEDURE construct_binary


  MODULE PROCEDURE destruct_binary

    !************************************************
    !
    !# Destructs the LORENE Bin_NS object and frees
    !  the pointer [[bns:bns_ptr]] pointing to it
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    !PRINT *, "** Executing the destruct_binary subroutine."

    IF ( C_ASSOCIATED( THIS% bns_ptr ) ) THEN

      CALL destruct_bin_ns( THIS% bns_ptr )
      THIS% bns_ptr = C_NULL_PTR

    ENDIF

    !PRINT *, "** Subroutine destruct_binary executed."
    !PRINT *

  END PROCEDURE destruct_binary


END SUBMODULE bns_constructor
