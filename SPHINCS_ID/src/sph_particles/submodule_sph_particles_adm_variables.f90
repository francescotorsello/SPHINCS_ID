! File:         submodule_sph_particles_adm_variables.f90
! Authors:      Francesco Torsello (FT)
!************************************************************************
! Copyright (C) 2020, 2021, 2022 Francesco Torsello                     *
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

SUBMODULE (sph_particles) adm_variables

  !************************************************
  !
  !# This SUBMODULE contains the implementation
  !  of the method compute_adm_variables of TYPE particles.
  !
  !  FT 08.04.2022
  !
  !************************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE compute_adm_momentum_fluid_canmom

    !************************************************
    !
    !# Computes an estimate of the \(\mathrm{ADM}\)
    !  linear momentum using the canonical momentum
    !  per baryon on the particles
    !  @todo add reference
    !
    !  FT 08.04.2022
    !
    !************************************************

    USE tensor,               ONLY: jx, jy, jz
    USE constants,            ONLY: amu, MSun
    USE utility,              ONLY: spatial_vector_norm_sym3x3, &
                                    zero, one, two

    IMPLICIT NONE

    INTEGER, PARAMETER:: unit_recovery= 34956

    INTEGER:: a, j

    DOUBLE PRECISION:: det, shift_norm2

    LOGICAL, PARAMETER:: debug= .FALSE.

    PRINT *, " * Estimating the ADM linear momentum of the fluid using the ", &
             "canonical SPH momentum per baryon on the particles..."
    PRINT *

    adm_mom= zero
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP             SHARED( npart, nu, lapse, shift, g3, s_l ) &
    !$OMP             PRIVATE( a, det, shift_norm2, j ) &
    !$OMP             REDUCTION( +: adm_mom )
    DO a= 1, npart, 1

      CALL spatial_vector_norm_sym3x3( g3(:,a), shift(:,a), shift_norm2 )

      DO j= jx, jz, 1

        adm_mom(j)= adm_mom(j) &
            - ( nu(a)*amu/Msun )*( shift_norm2/(lapse(a)**two) - one )*s_l(j,a)

      ENDDO

    ENDDO
    !$OMP END PARALLEL DO

    IF( debug ) PRINT *, "4"

  END PROCEDURE compute_adm_momentum_fluid_canmom


  MODULE PROCEDURE compute_adm_momentum_fluid_fields

    !************************************************
    !
    !# Computes an estimate of the \(\mathrm{ADM}\)
    !  linear momentum using the |sph| fields
    !  on the particles
    !  @todo add reference
    !
    !  FT 12.04.2022
    !
    !************************************************

    USE tensor,               ONLY: jx, jy, jz, n_sym4x4, lower_index_4vector
    USE constants,            ONLY: amu, MSun
    USE utility,              ONLY: compute_g4, determinant_sym4x4, &
                                    spatial_vector_norm_sym3x3, zero, one, two

    IMPLICIT NONE

    INTEGER, PARAMETER:: unit_recovery= 34956

    INTEGER:: a, j

    DOUBLE PRECISION:: det, shift_norm2

    DOUBLE PRECISION, DIMENSION(n_sym4x4,npart)  :: g4
    DOUBLE PRECISION, DIMENSION(0:3)             :: v_u
    DOUBLE PRECISION, DIMENSION(0:3,npart)       :: v_l

    LOGICAL, PARAMETER:: debug= .FALSE.

    adm_mom= zero
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP             SHARED( npart, nu, lapse, shift_x, shift_y, shift_z, &
    !$OMP                     theta, u, pr, nlrf, vel_u, &
    !$OMP                     v_l, g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, g4 ) &
    !$OMP             PRIVATE( a, det, v_u, shift_norm2, j ) &
    !$OMP             REDUCTION( +: adm_mom )
    DO a= 1, npart, 1

      CALL compute_g4( lapse(a), [shift_x(a),shift_y(a),shift_z(a)], &
                       [g_xx(a),g_xy(a),g_xz(a),g_yy(a),g_yz(a),g_zz(a)], &
                       g4(:,a) )

      CALL determinant_sym4x4( g4(:,a), det )
      IF( ABS(det) < 1D-10 )THEN
        PRINT *, "** ERROR! The determinant of the spacetime metric is " &
                 // "effectively 0 at particle ", a
        PRINT *, " * det= ", det
        PRINT *, " * Stopping..."
        STOP
      ELSEIF( det > 0 )THEN
        PRINT *, "** ERROR! The determinant of the spacetime metric is " &
                 // "positive at particle ", a
        PRINT *, " * det= ", det
        PRINT *, " * Stopping..."
        STOP
      ENDIF

      v_u= [one, vel_u(:,a)]
      CALL lower_index_4vector( v_u, g4(:,a), v_l(:,a) )

      CALL spatial_vector_norm_sym3x3( &
                    [g_xx(a),g_xy(a),g_xz(a),g_yy(a),g_yz(a),g_zz(a)], &
                    [shift_x(a),shift_y(a),shift_z(a)], shift_norm2 )

      DO j= jx, jz, 1

        adm_mom(j)= adm_mom(j) &
            - ( nu(a)*amu/Msun )*( shift_norm2/(lapse(a)**two) - one ) &
              *theta(a)*( one + u(a) + pr(a)/nlrf(a) )*v_l(j,a)

      ENDDO

    ENDDO
    !$OMP END PARALLEL DO

  END PROCEDURE compute_adm_momentum_fluid_fields


END SUBMODULE adm_variables
