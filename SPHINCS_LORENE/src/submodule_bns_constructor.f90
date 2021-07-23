! File:         submodule_bns_constructor.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (bns_id) bns_constructor

  !*****************************************************
  !                                                    *
  ! Implementation the constructor of TYPE bns and the *
  ! PROCEDURES it calls                                *
  !                                                    *
  ! FT 23.10.2020                                      *
  !                                                    *
  !*****************************************************


  IMPLICIT NONE


  CONTAINS


  !
  !-- Implementation of the constructor of the bns object
  !
  MODULE PROCEDURE construct_bns

    IMPLICIT NONE

    INTEGER, SAVE:: bns_counter= 1

    !DOUBLE PRECISION:: tmp

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

    IMPLICIT NONE

    !PRINT *, "Inside destructor of bns."
    !PRINT *

    ! Deallocate memory, if allocated
    CALL THIS% deallocate_lorene_id_memory()

  END PROCEDURE destruct_bns


  MODULE PROCEDURE construct_binary

    IMPLICIT NONE

    CHARACTER(KIND= C_CHAR, LEN= 7):: default_case

    LOGICAL:: exist

    !PRINT *, "** Executing the construct_binary subroutine..."

    IF ( C_ASSOCIATED( THIS% bns_ptr ) ) THEN

      CALL destruct_bin_ns( THIS% bns_ptr )

    ENDIF

    !
    !-- If the name of the LORENE binary file resu_file is given as argument to
    !-- construct_binary, use it. Otherwise, give the string "read_it"
    !-- to construct_bin_ns as argument, which makes LORENE read the name of
    !-- the file from the parameter file read_bin_ns.par
    !
    IF( PRESENT( resu_file ) )THEN

      INQUIRE( FILE= resu_file, EXIST= exist )

      IF( exist )THEN

        THIS% bns_ptr = construct_bin_ns( resu_file//C_NULL_CHAR )

      ELSE

        PRINT *, "** ERROR in bns SUBROUTINE construct_binary: File ", &
                 resu_file, "cannot be found!"
        PRINT *
        STOP

      ENDIF

    ELSE

      default_case= "read_it"
      THIS% bns_ptr = construct_bin_ns( default_case//C_NULL_CHAR )

    ENDIF

    PRINT *, "17"


    !PRINT *, "** Subroutine construct_binary executed."
    !PRINT *

  END PROCEDURE construct_binary


  MODULE PROCEDURE destruct_binary

    !************************************************
    !                                               *
    ! Destruct the LORENE Bin_NS object and free    *
    ! the pointeri pointing to it                   *
    !                                               *
    ! FT                                            *
    !                                               *
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
