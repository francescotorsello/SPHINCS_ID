! File:         submodule_bns_bindings.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (bns_id) bns_bindings

        !**********************************************
        !                                             *
        ! Define SUBROUTINES that create and destroy  *
        ! a Lorene::Bin_NS bin_ns object from LORENE, *
        ! and call its get_fields member method       *
        ! that extracts the LORENE initial data       *
        ! from the spectral data saved in the 'resu'  *
        ! binary files produced by LORENE             *
        !                                             *
        ! https://community.intel.com/t5/Intel-Fortran-Compiler/Calling-C-cpp-bin_nsects-from-a-Fortran-subroutine/td-p/1110556
        !                                             *
        ! FT 5.10.2020                                *
        !                                             *
        !                                             *
        !                                             *
        ! FT 9.11.2020                                *
        !                                             *
        !**********************************************

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, &
                                           C_CHAR, C_NULL_CHAR, &
                                           C_PTR, C_NULL_PTR, C_ASSOCIATED

    IMPLICIT NONE


    INTERFACE

        FUNCTION construct_bin_ns( c_resu_file ) RESULT( optr ) &
            BIND(C, NAME= "construct_bin_ns")

            IMPORT :: C_PTR, C_CHAR

            IMPLICIT NONE

            ! Argument list
            CHARACTER(KIND= C_CHAR), DIMENSION(*), INTENT(IN), OPTIONAL :: &
                                                                    c_resu_file
            ! Function result
            TYPE(C_PTR) :: optr

        END FUNCTION construct_bin_ns

    END INTERFACE

    INTERFACE

        SUBROUTINE get_lorene_id( optr, &
                                  x, y, z, &
                                  lapse, &
                                  shift_x, shift_y, shift_z, &
                                  g_diag, &
                                  k_xx, k_xy, k_xz, &
                                  k_yy, k_yz, k_zz, &
                                  baryon_density, &
                                  energy_density, &
                                  specific_energy, &
                                  u_euler_x, u_euler_y, u_euler_z ) &
            BIND(C, NAME= "get_lorene_id")

            IMPORT :: C_DOUBLE, C_PTR

            IMPLICIT NONE

            ! Argument list
            TYPE(C_PTR),    INTENT(IN),  VALUE :: optr
            REAL(C_DOUBLE), INTENT(IN),  VALUE :: x
            REAL(C_DOUBLE), INTENT(IN),  VALUE :: y
            REAL(C_DOUBLE), INTENT(IN),  VALUE :: z
            REAL(C_DOUBLE), INTENT(OUT) :: lapse
            REAL(C_DOUBLE), INTENT(OUT) :: shift_x
            REAL(C_DOUBLE), INTENT(OUT) :: shift_y
            REAL(C_DOUBLE), INTENT(OUT) :: shift_z
            REAL(C_DOUBLE), INTENT(OUT) :: g_diag
            REAL(C_DOUBLE), INTENT(OUT) :: k_xx
            REAL(C_DOUBLE), INTENT(OUT) :: k_xy
            REAL(C_DOUBLE), INTENT(OUT) :: k_xz
            REAL(C_DOUBLE), INTENT(OUT) :: k_yy
            REAL(C_DOUBLE), INTENT(OUT) :: k_yz
            REAL(C_DOUBLE), INTENT(OUT) :: k_zz
            REAL(C_DOUBLE), INTENT(OUT) :: baryon_density
            REAL(C_DOUBLE), INTENT(OUT) :: energy_density
            REAL(C_DOUBLE), INTENT(OUT) :: specific_energy
            REAL(C_DOUBLE), INTENT(OUT) :: u_euler_x
            REAL(C_DOUBLE), INTENT(OUT) :: u_euler_y
            REAL(C_DOUBLE), INTENT(OUT) :: u_euler_z

        END SUBROUTINE get_lorene_id

    END INTERFACE

    INTERFACE

        SUBROUTINE get_lorene_id_no_k( optr, &
                                       x, y, z, &
                                       lapse, &
                                       shift_x, shift_y, shift_z, &
                                       g_diag, &
                                       baryon_density, &
                                       energy_density, &
                                       specific_energy, &
                                       u_euler_x, u_euler_y, u_euler_z ) &
            BIND(C, NAME= "get_lorene_id_particles")

            IMPORT :: C_DOUBLE, C_PTR

            IMPLICIT NONE

            ! Argument list
            TYPE(C_PTR),    INTENT(IN),  VALUE :: optr
            REAL(C_DOUBLE), INTENT(IN),  VALUE :: x
            REAL(C_DOUBLE), INTENT(IN),  VALUE :: y
            REAL(C_DOUBLE), INTENT(IN),  VALUE :: z
            REAL(C_DOUBLE), INTENT(OUT) :: lapse
            REAL(C_DOUBLE), INTENT(OUT) :: shift_x
            REAL(C_DOUBLE), INTENT(OUT) :: shift_y
            REAL(C_DOUBLE), INTENT(OUT) :: shift_z
            REAL(C_DOUBLE), INTENT(OUT) :: g_diag
            REAL(C_DOUBLE), INTENT(OUT) :: baryon_density
            REAL(C_DOUBLE), INTENT(OUT) :: energy_density
            REAL(C_DOUBLE), INTENT(OUT) :: specific_energy
            REAL(C_DOUBLE), INTENT(OUT) :: u_euler_x
            REAL(C_DOUBLE), INTENT(OUT) :: u_euler_y
            REAL(C_DOUBLE), INTENT(OUT) :: u_euler_z

        END SUBROUTINE get_lorene_id_no_k

    END INTERFACE

    INTERFACE

        FUNCTION get_lorene_mass_density( optr, x, y, z ) RESULT( res ) &
            BIND(C, NAME= "get_mass_density")

            IMPORT :: C_DOUBLE, C_PTR

            IMPLICIT NONE

            ! Argument list
            TYPE(C_PTR),    INTENT(IN),  VALUE :: optr
            REAL(C_DOUBLE), INTENT(IN),  VALUE :: x
            REAL(C_DOUBLE), INTENT(IN),  VALUE :: y
            REAL(C_DOUBLE), INTENT(IN),  VALUE :: z
            ! Function result
            REAL(C_DOUBLE) :: res

        END FUNCTION get_lorene_mass_density

    END INTERFACE

    INTERFACE

        SUBROUTINE get_lorene_id_params( optr, &
                                         angular_vel, &
                                         distance, &
                                         distance_com, &
                                         mass1, &
                                         mass2, &
                                         adm_mass, &
                                         angular_momentum, &
                                         radius1_x_comp, &
                                         radius1_y, &
                                         radius1_z, &
                                         radius1_x_opp, &
                                         radius2_x_comp, &
                                         radius2_y, &
                                         radius2_z, &
                                         radius2_x_opp ) &
            BIND(C, NAME= "get_lorene_id_params")

            IMPORT :: C_DOUBLE, C_PTR

            IMPLICIT NONE

            ! Argument list
            TYPE(C_PTR),    INTENT(IN),  VALUE :: optr
            REAL(C_DOUBLE), INTENT(OUT) :: angular_vel
            REAL(C_DOUBLE), INTENT(OUT) :: distance
            REAL(C_DOUBLE), INTENT(OUT) :: distance_com
            REAL(C_DOUBLE), INTENT(OUT) :: mass1
            REAL(C_DOUBLE), INTENT(OUT) :: mass2
            REAL(C_DOUBLE), INTENT(OUT) :: adm_mass
            REAL(C_DOUBLE), INTENT(OUT) :: angular_momentum
            REAL(C_DOUBLE), INTENT(OUT) :: radius1_x_comp
            REAL(C_DOUBLE), INTENT(OUT) :: radius1_y
            REAL(C_DOUBLE), INTENT(OUT) :: radius1_z
            REAL(C_DOUBLE), INTENT(OUT) :: radius1_x_opp
            REAL(C_DOUBLE), INTENT(OUT) :: radius2_x_comp
            REAL(C_DOUBLE), INTENT(OUT) :: radius2_y
            REAL(C_DOUBLE), INTENT(OUT) :: radius2_z
            REAL(C_DOUBLE), INTENT(OUT) :: radius2_x_opp

        END SUBROUTINE get_lorene_id_params

    END INTERFACE

    INTERFACE

        SUBROUTINE destruct_bin_ns( optr ) &
            BIND(C, NAME= "destruct_bin_ns")

            IMPORT :: C_PTR

            IMPLICIT NONE

            ! Argument list
            TYPE(C_PTR), INTENT(IN), VALUE :: optr

        END SUBROUTINE destruct_bin_ns

    END INTERFACE

END SUBMODULE bns_bindings

!TYPE(C_PTR), SAVE :: bin_ns = C_NULL_PTR
!
!
!PUBLIC :: construct_binary
!!PUBLIC :: import_id
!!PUBLIC :: import_id_particles
!!PUBLIC :: import_id_params
!!PUBLIC :: import_mass_density
!PUBLIC :: destruct_binary
!
!
!CONTAINS
!
!
!SUBROUTINE construct_binary( resu_file )
!
!    IMPLICIT NONE
!
!    INTEGER:: itr
!    CHARACTER(KIND= C_CHAR, LEN= 7):: default_case
!    LOGICAL:: exist
!    ! Argument list
!    CHARACTER(KIND= C_CHAR, LEN=*), INTENT(IN), OPTIONAL:: resu_file
!
!    PRINT *, "** Executing the construct_binary subroutine..."
!
!    IF ( C_ASSOCIATED( bns_ptr ) ) THEN
!
!        CALL destruct_bin_ns( bns_ptr )
!
!    ENDIF
!
!    IF( PRESENT( resu_file ) )THEN
!
!        INQUIRE( FILE= resu_file, EXIST= exist )
!
!        IF( exist )THEN
!
!            bns_ptr = construct_bin_ns( resu_file//C_NULL_CHAR )
!
!        ELSE
!
!            PRINT *, "** ERROR: File ", resu_file, "cannot be found!"
!            PRINT *
!            STOP
!
!        ENDIF
!
!    ELSE
!
!        default_case= "read_it"
!        bns_ptr = construct_bin_ns( default_case//C_NULL_CHAR )
!
!    ENDIF
!
!    PRINT *, "** Subroutine construct_binary executed."
!    PRINT *
!
!    RETURN
!
!END SUBROUTINE construct_binary
!
!SUBROUTINE import_id( x, y, z, &
!                      lapse, &
!                      shift_x, shift_y, shift_z, &
!                      g_diag, &
!                      k_xx, k_xy, k_xz, &
!                      k_yy, k_yz, k_zz, &
!                      baryon_density, &
!                      energy_density, &
!                      specific_energy, &
!                      u_euler_x, u_euler_y, u_euler_z )
!
!    USE, INTRINSIC :: ISO_C_BINDING
!
!    ! Argument list
!    REAL(C_DOUBLE), INTENT(IN),  VALUE :: x
!    REAL(C_DOUBLE), INTENT(IN),  VALUE :: y
!    REAL(C_DOUBLE), INTENT(IN),  VALUE :: z
!    REAL(C_DOUBLE), INTENT(OUT) :: lapse
!    REAL(C_DOUBLE), INTENT(OUT) :: shift_x
!    REAL(C_DOUBLE), INTENT(OUT) :: shift_y
!    REAL(C_DOUBLE), INTENT(OUT) :: shift_z
!    REAL(C_DOUBLE), INTENT(OUT) :: g_diag
!    REAL(C_DOUBLE), INTENT(OUT) :: k_xx
!    REAL(C_DOUBLE), INTENT(OUT) :: k_xy
!    REAL(C_DOUBLE), INTENT(OUT) :: k_xz
!    REAL(C_DOUBLE), INTENT(OUT) :: k_yy
!    REAL(C_DOUBLE), INTENT(OUT) :: k_yz
!    REAL(C_DOUBLE), INTENT(OUT) :: k_zz
!    REAL(C_DOUBLE), INTENT(OUT) :: baryon_density
!    REAL(C_DOUBLE), INTENT(OUT) :: energy_density
!    REAL(C_DOUBLE), INTENT(OUT) :: specific_energy
!    REAL(C_DOUBLE), INTENT(OUT) :: u_euler_x
!    REAL(C_DOUBLE), INTENT(OUT) :: u_euler_y
!    REAL(C_DOUBLE), INTENT(OUT) :: u_euler_z
!
!    IF ( C_ASSOCIATED( bin_ns ) ) THEN
!
!        CALL get_lorene_id( bin_ns, &
!                            x, y, z, &
!                            lapse, &
!                            shift_x, shift_y, shift_z, &
!                            g_diag, &
!                            k_xx, k_xy, k_xz, &
!                            k_yy, k_yz, k_zz, &
!                            baryon_density, &
!                            energy_density, &
!                            specific_energy, &
!                            u_euler_x, u_euler_y, u_euler_z )
!
!    ENDIF
!
!    RETURN
!
!END SUBROUTINE import_id
!
!SUBROUTINE import_id_particles( x, y, z, &
!                      lapse, &
!                      shift_x, shift_y, shift_z, &
!                      g_diag, &
!                      baryon_density, &
!                      energy_density, &
!                      specific_energy, &
!                      u_euler_x, u_euler_y, u_euler_z )
!
!    USE, INTRINSIC :: ISO_C_BINDING
!
!    ! Argument list
!    REAL(C_DOUBLE), INTENT(IN),  VALUE :: x
!    REAL(C_DOUBLE), INTENT(IN),  VALUE :: y
!    REAL(C_DOUBLE), INTENT(IN),  VALUE :: z
!    REAL(C_DOUBLE), INTENT(OUT) :: lapse
!    REAL(C_DOUBLE), INTENT(OUT) :: shift_x
!    REAL(C_DOUBLE), INTENT(OUT) :: shift_y
!    REAL(C_DOUBLE), INTENT(OUT) :: shift_z
!    REAL(C_DOUBLE), INTENT(OUT) :: g_diag
!    REAL(C_DOUBLE), INTENT(OUT) :: baryon_density
!    REAL(C_DOUBLE), INTENT(OUT) :: energy_density
!    REAL(C_DOUBLE), INTENT(OUT) :: specific_energy
!    REAL(C_DOUBLE), INTENT(OUT) :: u_euler_x
!    REAL(C_DOUBLE), INTENT(OUT) :: u_euler_y
!    REAL(C_DOUBLE), INTENT(OUT) :: u_euler_z
!
!    IF ( C_ASSOCIATED( bin_ns ) ) THEN
!
!        CALL get_lorene_id_no_k( bin_ns, &
!                            x, y, z, &
!                            lapse, &
!                            shift_x, shift_y, shift_z, &
!                            g_diag, &
!                            baryon_density, &
!                            energy_density, &
!                            specific_energy, &
!                            u_euler_x, u_euler_y, u_euler_z )
!
!    ENDIF
!
!    RETURN
!
!END SUBROUTINE import_id_particles
!
!FUNCTION import_mass_density( x, y, z ) RESULT( res )
!
!    USE, INTRINSIC :: ISO_C_BINDING
!
!    ! Argument list
!    REAL(C_DOUBLE), INTENT(IN),  VALUE :: x
!    REAL(C_DOUBLE), INTENT(IN),  VALUE :: y
!    REAL(C_DOUBLE), INTENT(IN),  VALUE :: z
!    ! Function result
!    REAL(C_DOUBLE):: res
!
!    IF ( C_ASSOCIATED( bin_ns ) ) THEN
!
!        res= get_lorene_mass_density( bin_ns, x, y, z )
!
!    ENDIF
!
!    RETURN
!
!END FUNCTION import_mass_density
!
!SUBROUTINE import_id_params( angular_vel, &
!                             distance, &
!                             distance_com, &
!                             mass1, &
!                             mass2, &
!                             adm_mass, &
!                             angular_momentum, &
!                             radius1_x_comp, &
!                             radius1_y, &
!                             radius1_z, &
!                             radius1_x_opp, &
!                             radius2_x_comp, &
!                             radius2_y, &
!                             radius2_z, &
!                             radius2_x_opp )
!
!    USE, INTRINSIC :: ISO_C_BINDING
!
!    ! Argument list
!    REAL(C_DOUBLE), INTENT(OUT) :: angular_vel
!    REAL(C_DOUBLE), INTENT(OUT) :: distance
!    REAL(C_DOUBLE), INTENT(OUT) :: distance_com
!    REAL(C_DOUBLE), INTENT(OUT) :: mass1
!    REAL(C_DOUBLE), INTENT(OUT) :: mass2
!    REAL(C_DOUBLE), INTENT(OUT) :: adm_mass
!    REAL(C_DOUBLE), INTENT(OUT) :: angular_momentum
!    REAL(C_DOUBLE), INTENT(OUT) :: radius1_x_comp
!    REAL(C_DOUBLE), INTENT(OUT) :: radius1_y
!    REAL(C_DOUBLE), INTENT(OUT) :: radius1_z
!    REAL(C_DOUBLE), INTENT(OUT) :: radius1_x_opp
!    REAL(C_DOUBLE), INTENT(OUT) :: radius2_x_comp
!    REAL(C_DOUBLE), INTENT(OUT) :: radius2_y
!    REAL(C_DOUBLE), INTENT(OUT) :: radius2_z
!    REAL(C_DOUBLE), INTENT(OUT) :: radius2_x_opp
!
!    IF ( C_ASSOCIATED( bin_ns ) ) THEN
!
!        CALL get_lorene_id_params( bin_ns, &
!                                   angular_vel, &
!                                   distance, &
!                                   distance_com, &
!                                   mass1, &
!                                   mass2, &
!                                   adm_mass, &
!                                   angular_momentum, &
!                                   radius1_x_comp, &
!                                   radius1_y, &
!                                   radius1_z, &
!                                   radius1_x_opp, &
!                                   radius2_x_comp, &
!                                   radius2_y, &
!                                   radius2_z, &
!                                   radius2_x_opp )
!
!    ENDIF
!
!    RETURN
!
!END SUBROUTINE import_id_params
!
!SUBROUTINE destruct_binary( bns_ptr )
!
!    TYPE(C_PTR), INTENT(IN OUT) :: bns_ptr
!
!    PRINT *, "** Executing the destruct_binary subroutine."
!
!    IF ( C_ASSOCIATED( bin_ns ) ) THEN
!
!        CALL destruct_bin_ns( bin_ns )
!        bin_ns = C_NULL_PTR
!
!    ENDIF
!
!    PRINT *, "** Subroutine destruct_binary executed."
!    PRINT *
!
!    RETURN
!
!END SUBROUTINE destruct_binary
