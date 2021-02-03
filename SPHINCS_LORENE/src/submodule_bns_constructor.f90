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

    DOUBLE PRECISION:: tmp

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

    INTEGER:: itr

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

        PRINT *, "** ERROR: File ", resu_file, "cannot be found!"
        PRINT *
        STOP

      ENDIF

    ELSE

      default_case= "read_it"
      THIS% bns_ptr = construct_bin_ns( default_case//C_NULL_CHAR )

    ENDIF

    !PRINT *, "** Subroutine construct_binary executed."
    !PRINT *

  END PROCEDURE construct_binary

  MODULE PROCEDURE import_id_params

    !*****************************************************
    !                                                    *
    ! Store the parameters of the binary neutron         *
    ! stars' LORENE ID into member variables             *
    !                                                    *
    ! FT 5.10.2020                                       *
    !                                                    *
    !*****************************************************

    USE constants, ONLY: Msun_geo

    IMPLICIT NONE

    PRINT *, "** Executing the import_lorene_id_params subroutine..."

    CALL get_lorene_id_params( THIS% bns_ptr, &
                               THIS% angular_vel, &
                               THIS% distance, &
                               THIS% distance_com, &
                               THIS% mass1, &
                               THIS% mass2, &
                               THIS% adm_mass, &
                               THIS% angular_momentum, &
                               THIS% radius1_x_comp, &
                               THIS% radius1_y, &
                               THIS% radius1_z, &
                               THIS% radius1_x_opp, &
                               THIS% radius2_x_comp, &
                               THIS% radius2_y, &
                               THIS% radius2_z, &
                               THIS% radius2_x_opp, &
                               THIS% eos1, &
                               THIS% eos2, &
                               THIS% npeos_1, &
                               THIS% gamma0_1, &
                               THIS% gamma1_1, &
                               THIS% gamma2_1, &
                               THIS% gamma3_1, &
                               THIS% kappa0_1, &
                               THIS% kappa1_1, &
                               THIS% kappa2_1, &
                               THIS% kappa3_1, &
                               THIS% logP1_1,  &
                               THIS% logRho0_1,&
                               THIS% logRho1_1,&
                               THIS% logRho2_1,&
                               THIS% npeos_2,  &
                               THIS% gamma0_2, &
                               THIS% gamma1_2, &
                               THIS% gamma2_2, &
                               THIS% gamma3_2, &
                               THIS% kappa0_2, &
                               THIS% kappa1_2, &
                               THIS% kappa2_2, &
                               THIS% kappa3_2, &
                               THIS% logP1_2,  &
                               THIS% logRho0_2,&
                               THIS% logRho1_2,&
                               THIS% logRho2_2)

    ! Convert distances from LORENE units (km) to SPHINCS units (Msun_geo)
    ! See MODULE constants for the definition of Msun_geo
    THIS% distance      = THIS% distance/Msun_geo
    THIS% distance_com  = THIS% distance_com/Msun_geo
    THIS% radius1_x_comp= THIS% radius1_x_comp/Msun_geo
    THIS% radius1_y     = THIS% radius1_y/Msun_geo
    THIS% radius1_z     = THIS% radius1_z/Msun_geo
    THIS% radius1_x_opp = THIS% radius1_x_opp/Msun_geo
    THIS% radius2_x_comp= THIS% radius2_x_comp/Msun_geo
    THIS% radius2_y     = THIS% radius2_y/Msun_geo
    THIS% radius2_z     = THIS% radius2_z/Msun_geo
    THIS% radius2_x_opp = THIS% radius2_x_opp/Msun_geo

    CALL print_id_params( THIS )

    PRINT *, "** Subroutine import_lorene_id_params executed."
    PRINT *

  END PROCEDURE import_id_params

END SUBMODULE bns_constructor
