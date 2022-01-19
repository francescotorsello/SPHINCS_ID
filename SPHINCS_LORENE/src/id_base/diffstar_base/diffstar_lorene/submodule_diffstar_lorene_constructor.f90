! File:         submodule_diffstar_lorene_constructor.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (diffstar_lorene) constructor

  !*********************************************************
  !
  !# Implementation of the constructor and
  !  destructor of TYPE [[diffstarlorene]], and of the
  !  [[diffstarlorene]]-member
  !  PROCEDURES that call the C-bound PROCEDURES
  !  constructig and destructing the |lorene|
  !  |etdiffrot| object
  !
  !  FT 25.10.2021
  !
  !*********************************************************


  USE, INTRINSIC:: ISO_C_BINDING, ONLY: C_ASSOCIATED, C_NULL_CHAR, C_NULL_PTR


  IMPLICIT NONE


  CONTAINS


  !
  !-- Implementation of the constructor of the bns object
  !
  MODULE PROCEDURE construct_diffstarlorene

    !****************************************************
    !
    !# Constructs an object of TYPE [[diffstarlorene]]
    !
    !  FT 25.10.2021
    !
    !****************************************************

    IMPLICIT NONE

    INTEGER, SAVE:: diffstar_counter= 1

    CALL derived_type% set_n_matter(1)

    derived_type% construction_timer= timer( "drs_construction_timer" )

    ! Construct |lorene| |etdiffrot| object
    IF( PRESENT( filename ) )THEN
        CALL derived_type% construct_drs( filename )
    ELSE
        CALL derived_type% construct_drs()
    ENDIF

    ! Import the parameters of the binary system
    CALL import_diffstar_params( derived_type )

    ! Assign a unique identifier to the bns object
    derived_type% diffstar_identifier= diffstar_counter
    diffstar_counter= diffstar_counter + 1

    ! Do not use the geodesic gauge by default
    CALL derived_type% set_one_lapse ( .FALSE. )
    CALL derived_type% set_zero_shift( .FALSE. )

  END PROCEDURE construct_diffstarlorene


  !
  !-- Implementation of the destructor of the bns object
  !
  MODULE PROCEDURE destruct_diffstarlorene

    !***********************************************
    !
    !# Destructs an object of TYPE [[diffstarlorene]]
    !
    !  FT 25.10.2021
    !
    !***********************************************

    IMPLICIT NONE

    ! Deallocate memory
    CALL THIS% deallocate_diffstar_memory()

  END PROCEDURE destruct_diffstarlorene


  MODULE PROCEDURE construct_drs

    !***********************************************
    !
    !# Construct the |lorene| |etdiffrot| object
    !
    !  FT 25.10.2021
    !
    !***********************************************

    IMPLICIT NONE

    CHARACTER(KIND= C_CHAR, LEN= 7):: default_case
    LOGICAL:: exist

    !PRINT *, "** Executing the construct_binary subroutine..."

#ifdef __INTEL_COMPILER

    IF ( C_ASSOCIATED( THIS% diffstar_ptr ) ) THEN

      CALL destruct_etdiffrot( THIS% diffstar_ptr )

    ENDIF

#endif

    !
    !-- If the name of the |lorene| binary file resu_file is given as argument to
    !-- construct_binary, use it. Otherwise, give the string "read_it"
    !-- to construct_drs as argument, which makes |lorene| read the name of
    !-- the file from the parameter file read_bin_ns.par
    !
    IF( PRESENT( resu_file ) )THEN

      INQUIRE( FILE= resu_file, EXIST= exist )

      IF( exist )THEN

        CALL THIS% construction_timer% start_timer()
        THIS% diffstar_ptr = construct_etdiffrot( resu_file//C_NULL_CHAR )
        CALL THIS% construction_timer% stop_timer()

      ELSE

        PRINT *, "** ERROR in SUBROUTINE construct_binary: file ", &
                 resu_file, " cannot be found!"
        PRINT *
        STOP

      ENDIF

    ELSE

      default_case= "read_it"
      CALL THIS% construction_timer% start_timer()
      THIS% diffstar_ptr = construct_etdiffrot( default_case//C_NULL_CHAR )
      CALL THIS% construction_timer% stop_timer()

    ENDIF

  END PROCEDURE construct_drs


  MODULE PROCEDURE destruct_drs

    !************************************************
    !
    !# Destructs the |lorene| |etdiffrot| object and frees
    !  the pointer [[diffstar:diffstar_ptr]] pointing to it
    !
    !  FT 25.10.2021
    !
    !************************************************

    IMPLICIT NONE

    IF ( C_ASSOCIATED( THIS% diffstar_ptr ) ) THEN

      CALL destruct_etdiffrot( THIS% diffstar_ptr )
      THIS% diffstar_ptr = C_NULL_PTR

    ENDIF

  END PROCEDURE destruct_drs


END SUBMODULE constructor
