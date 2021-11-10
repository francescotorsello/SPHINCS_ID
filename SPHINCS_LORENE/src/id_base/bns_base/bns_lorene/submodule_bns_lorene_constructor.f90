! File:         submodule_bns_constructor.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (bns_lorene) bns_lorene_constructor

  !*********************************************************
  !
  !# Implementation of the constructor and
  !  destructor of TYPE [[bnslorene]], and of the
  !  [[bnslorene]]-member
  !  PROCEDURES that call the C-bound PROCEDURES
  !  constructig and destructing the |lorene|
  !  |binns| object
  !
  !  FT 23.10.2020
  !
  !*********************************************************


  IMPLICIT NONE


  CONTAINS


 ! MODULE PROCEDURE construct_bnslorene2
 !
 !   !****************************************************
 !   !
 !   !# Constructs an object of TYPE [[bnslorene]]
 !   !
 !   !  FT
 !   !
 !   !****************************************************
 !
 !   IMPLICIT NONE
 !
 !   CHARACTER(LEN=10) :: resu_file
 !
 !   derived_type => construct_bnslorene( resu_file )
 !
 ! END PROCEDURE construct_bnslorene2


  !
  !-- Implementation of the constructor of the bns object
  !
  MODULE PROCEDURE construct_bnslorene

    !****************************************************
    !
    !# Constructs an object of TYPE [[bnslorene]]
    !
    !  FT
    !
    !****************************************************

    IMPLICIT NONE

    INTEGER, SAVE:: bns_counter= 1

    !DOUBLE PRECISION:: tmp

    CALL derived_type% set_n_matter(2)

    derived_type% construction_timer= timer( "binary_construction_timer" )

    ! Construct |lorene| |binns| object
    IF( PRESENT( filename ) )THEN
        CALL derived_type% construct_binary( filename )
    ELSE
        CALL derived_type% construct_binary()
    ENDIF

    ! Import the parameters of the binary system
    CALL import_id_params( derived_type )

    ! Assign a unique identifier to the bns object
    derived_type% bns_identifier= bns_counter
    bns_counter= bns_counter + 1

    ! Do not use the geodesic gauge by default
    CALL derived_type% set_one_lapse ( .FALSE. )
    CALL derived_type% set_zero_shift( .FALSE. )

    !foo= derived_type

  END PROCEDURE construct_bnslorene


  !
  !-- Implementation of the destructor of the bns object
  !
  MODULE PROCEDURE destruct_bnslorene

    !***********************************************
    !
    !# Destructs an object of TYPE [[bnslorene]]
    !
    !  FT
    !
    !***********************************************

    IMPLICIT NONE

    !PRINT *, "Inside destructor of bns."
    !PRINT *

    ! Deallocate memory
    CALL THIS% deallocate_lorene_id_memory()

  END PROCEDURE destruct_bnslorene


  MODULE PROCEDURE construct_binary

    !***********************************************
    !
    !# Construct the |lorene| |binns| object
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
    !-- If the name of the |lorene| binary file resu_file is given as argument to
    !-- construct_binary, use it. Otherwise, give the string "read_it"
    !-- to construct_bin_ns as argument, which makes |lorene| read the name of
    !-- the file from the parameter file read_bin_ns.par
    !
    IF( PRESENT( resu_file ) )THEN

      INQUIRE( FILE= resu_file, EXIST= exist )

      IF( exist )THEN

        CALL THIS% construction_timer% start_timer()
        THIS% bns_ptr = construct_bin_ns( resu_file//C_NULL_CHAR )
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
      THIS% bns_ptr = construct_bin_ns( default_case//C_NULL_CHAR )
      CALL THIS% construction_timer% stop_timer()

    ENDIF

    !PRINT *, "** Subroutine construct_binary executed."
    !PRINT *

  END PROCEDURE construct_binary


  MODULE PROCEDURE destruct_binary

    !************************************************
    !
    !# Destructs the |lorene| |binns| object and frees
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


END SUBMODULE bns_lorene_constructor
