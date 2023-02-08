! File:         submodule_sph_particles_constructor_bin.f90
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

SUBMODULE (sph_particles) constructor_bin

  !************************************************
  !
  !# This SUBMODULE contains the implementation
  !  of the constructor of TYPE sph_particles
  !  from an |id| binary file.
  !
  !  FT 01.03.2022
  !
  !************************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE construct_particles_bin

    !**************************************************
    !
    !# Constructs a [[particles]] object from an |id|
    !  binary file
    !
    !  FT 01.03.2022
    !
    !**************************************************


    USE tensor,         ONLY: n_sym4x4
    USE utility,        ONLY: compute_g4, determinant_sym4x4, &
                              spacetime_vector_norm_sym4x4,zero, one
    USE input_output,   ONLY: read_options
    USE options,        ONLY: eos_str
    USE units,          ONLY: set_units
    USE pwp_EOS,        ONLY: select_EOS_parameters
    USE sph_variables,  ONLY: allocate_SPH_memory, deallocate_SPH_memory


    IMPLICIT NONE


    INTEGER:: a

    DOUBLE PRECISION:: g4(n_sym4x4), det, v_norm
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: tmp1
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: tmp2
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: tmp3
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: tmp4
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: tmp5
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: tmp6
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: tmp7
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: tmp8
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: tmp9
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: tmp10

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: sq_detg4


    PRINT *, "** Assigning the ID to the particles..."
    PRINT *

    CALL set_units('NSM')
    CALL read_options

    CALL parts% read_sphincs_dump_print_formatted( namefile_bin= namefile, &
                                                   save_data   = .TRUE. )

    ALLOCATE( tmp1 ( parts% npart ) )
    ALLOCATE( tmp2 ( parts% npart ) )
    ALLOCATE( tmp3 ( parts% npart ) )
    ALLOCATE( tmp4 ( parts% npart ) )
    ALLOCATE( tmp5 ( parts% npart ) )
    ALLOCATE( tmp6 ( parts% npart ) )
    ALLOCATE( tmp7 ( parts% npart ) )
    ALLOCATE( tmp8 ( parts% npart ) )
    ALLOCATE( tmp9 ( parts% npart ) )
    ALLOCATE( tmp10( parts% npart ) )

    CALL id% read_id_particles( parts% npart, &
                                parts% pos( 1, : ), &
                                parts% pos( 2, : ), &
                                parts% pos( 3, : ), &
                                parts% lapse, &
                                parts% shift_x, &
                                parts% shift_y, &
                                parts% shift_z, &
                                parts% g_xx, &
                                parts% g_xy, &
                                parts% g_xz, &
                                parts% g_yy, &
                                parts% g_yz, &
                                parts% g_zz, &
                                tmp4, &
                                tmp5, &
                                tmp6, &
                                tmp7, &
                                tmp8, &
                                tmp9, &
                                tmp10 )

    DEALLOCATE( tmp1  )
    DEALLOCATE( tmp2  )
    DEALLOCATE( tmp3  )
    DEALLOCATE( tmp4  )
    DEALLOCATE( tmp5  )
    DEALLOCATE( tmp6  )
    DEALLOCATE( tmp7  )
    DEALLOCATE( tmp8  )
    DEALLOCATE( tmp9  )
    DEALLOCATE( tmp10 )

    ALLOCATE( sq_detg4( parts% npart ) )

    DO a= 1, parts% npart, 1

      CALL compute_g4( parts% lapse(a), &
            [parts% shift_x(a), parts% shift_y(a), parts% shift_z(a)], &
            [parts% g_xx(a), parts% g_xy(a), parts% g_xz(a), &
             parts% g_yy(a), parts% g_yz(a), parts% g_zz(a)], &
             g4(1:n_sym4x4) )

      CALL determinant_sym4x4( g4(1:n_sym4x4), det )
      IF( ABS(det) < 1D-10 )THEN
          PRINT *, "** ERROR! The determinant of the spacetime metric is " &
                   // "effectively 0 at particle ", a
          PRINT *
          STOP
      ELSEIF( det > 0 )THEN
          PRINT *, "** ERROR! The determinant of the spacetime metric is " &
                   // "positive at particle ", a
          PRINT *
          STOP
      ENDIF
      sq_detg4(a)= SQRT(-det)

      IF( parts% theta(a) < zero )THEN
        PRINT *, "** ERROR! The computing frame generalized Lorentz factor ", &
                 "is negative at particle  ", a
        PRINT *, " * Its value is ", parts% theta(a)
        PRINT *, " * Stopping.."
        PRINT *
        STOP
      ELSEIF( parts% theta(a) == zero )THEN
        PRINT *, "** ERROR! The computing frame generalized Lorentz factor ", &
                 "is zero at particle ", a
        PRINT *, " * Its value is ", parts% theta(a)
        PRINT *, " * Stopping.."
        PRINT *
        STOP
      ENDIF

      v_norm= zero
      CALL spacetime_vector_norm_sym4x4( g4, parts% v(:,a), v_norm )
      IF( v_norm > zero )THEN
        PRINT *, "** ERROR! The computing frame particle 4-velocity is ", &
                 "spacelike at particle ", a
        PRINT *, " * Its norm is ", v_norm
        PRINT *, " * Stopping.."
        PRINT *
        STOP
      ELSEIF( v_norm == zero )THEN
        PRINT *, "** ERROR! The computing frame particle 4-velocity is ", &
                 "null at particle ", a
        PRINT *, " * Its norm is ", v_norm
        PRINT *, " * Stopping.."
        PRINT *
        STOP
      ENDIF

    ENDDO

    parts% nstar_sph= parts% nlrf_sph*sq_detg4*parts% theta

    CALL allocate_sph_memory

    !CALL OMP_SET_NUM_THREADS(1)

    CALL select_EOS_parameters( eos_str )

    CALL parts% test_recovery( parts% npart,       &
                               parts% pos,         &
                               parts% nlrf_sph,    &
                               parts% u_sph,       &
                               parts% pressure_sph, &
                               parts% v(1:3,:),    &
                               parts% theta,       &
                               parts% nstar_sph )

    DEALLOCATE( sq_detg4 )
    CALL deallocate_sph_memory

  END PROCEDURE construct_particles_bin


END SUBMODULE constructor_bin
