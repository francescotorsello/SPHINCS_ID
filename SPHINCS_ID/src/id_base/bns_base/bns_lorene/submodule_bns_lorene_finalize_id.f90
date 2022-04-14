! File:         submodule_bnslorene_finalize_id.f90
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

SUBMODULE (bns_lorene) finalize_id

  !*********************************************************
  !
  !# Implementation of the PROCEDURE that act on the |id|
  !  after it is set up on the mesh and/or on the particles,
  !  to finalize its preparation for |sphincsbssn|
  !
  !  FT 14.04.2022
  !
  !*********************************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE correct_adm_linear_momentum

    !***********************************************
    !
    !# Correct the velocity and the generalized
    !  Lorentz factor, so that the \(\mathrm{ADM}\)
    !  linear momentum of the |bns| is 0
    !  @todo: add reference
    !
    !  FT 14.04.2022
    !
    !***********************************************

    USE constants,  ONLY: zero, one, two
    USE tensor,     ONLY: jx, jy, jz, n_sym4x4
    USE matrix,     ONLY: invert_3x3_matrix
    USE utility,    ONLY: compute_g4, spacetime_vector_norm_sym4x4, &
                          spatial_vector_norm_sym3x3, is_finite_number

    IMPLICIT NONE

    INTEGER:: a, j

    DOUBLE PRECISION:: shift_norm2
    DOUBLE PRECISION, DIMENSION(3,npart):: vel_u_corr
    DOUBLE PRECISION:: g4(n_sym4x4)
    DOUBLE PRECISION, DIMENSION(3,3):: g3mat
    DOUBLE PRECISION, DIMENSION(3,3):: g3mat_inv

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( npart, lapse, shift_x, shift_y, shift_z, vel_u, &
    !$OMP                     vel_u_corr, g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, &
    !$OMP                     theta, adm_mom_error ) &
    !$OMP             PRIVATE( a, j, shift_norm2, g4, g3mat, g3mat_inv &
    !$OMP                       )
    DO a= 1, npart, 1

      CALL spatial_vector_norm_sym3x3( &
                    [g_xx(a),g_xy(a),g_xz(a),g_yy(a),g_yz(a),g_zz(a)], &
                    [shift_x(a),shift_y(a),shift_z(a)], shift_norm2 )

      g3mat(1,1)= g_xx(a)
      g3mat(1,2)= g_xy(a)
      g3mat(1,3)= g_xz(a)
      g3mat(2,1)= g_xy(a)
      g3mat(2,2)= g_yy(a)
      g3mat(2,3)= g_yz(a)
      g3mat(3,1)= g_xz(a)
      g3mat(3,2)= g_yz(a)
      g3mat(3,3)= g_zz(a)

      CALL invert_3x3_matrix( g3mat, g3mat_inv )

      ! That's the correction to the COVARIANT velocity, isn't it?

      vel_u_corr(1,a)= adm_mom_error(1) &
                     /( two*g3mat_inv(1,1)*NORM2(adm_mom_error)**two ) &
                     *( - one + SQRT( one - &
                        two*two*g3mat_inv(1,1)*NORM2(adm_mom_error)**two &
                        *( -lapse(a)**two + shift_norm2 ) ) )

      vel_u_corr(2,a)= vel_u_corr(1,a)*adm_mom_error(2)/adm_mom_error(1)
      vel_u_corr(3,a)= vel_u_corr(1,a)*adm_mom_error(3)/adm_mom_error(1)

      IF( a == 1 )THEN
        PRINT *, "adm_mom_error=", adm_mom_error
        PRINT *, "vel_u_corr=", vel_u_corr
      ENDIF

      DO j= jx, jz, 1

        vel_u(j,a)= vel_u(j,a) + vel_u_corr(j,a)

      ENDDO

      CALL compute_g4( lapse(a), [shift_x(a),shift_y(a),shift_z(a)], &
                       [g_xx(a),g_xy(a),g_xz(a),g_yy(a),g_yz(a),g_zz(a)], g4 )

      CALL spacetime_vector_norm_sym4x4( g4, [one, vel_u(:,a)], theta(a) )
      IF( .NOT.is_finite_number(theta(a)) )THEN
        PRINT *, "** ERROR! The spacetime norm of vel_u is ", theta(a), &
                 "at particle ", a, &
                 "in SUBROUTINE correct_adm_linear_momentum"
        PRINT *, " * Stopping..."
        PRINT *
        STOP
      ENDIF
      IF( theta(a) > zero )THEN
        PRINT *, "** ERROR! The spacetime norm of vel_u is ", theta(a), &
                 "at particle ", a, "(that is, vel_u is spacelike)", &
                 "in SUBROUTINE correct_adm_linear_momentum"
        PRINT *, " * Stopping..."
        PRINT *
        STOP
      ENDIF
      IF( theta(a) == zero )THEN
        PRINT *, "** ERROR! The spacetime norm of vel_u is ", theta(a), &
                 "at particle ", a, "(that is, vel_u is null)", &
                 "in SUBROUTINE correct_adm_linear_momentum"
        PRINT *, " * Stopping..."
        PRINT *
        STOP
      ENDIF

      theta(a)= one/SQRT(-theta(a))
      IF( .NOT.is_finite_number(theta(a)) )THEN
        PRINT *, "** ERROR! The generalized Lorentz factor is ", theta(a), &
                 "at particle ", a, &
                 "in SUBROUTINE correct_adm_linear_momentum"
        PRINT *, " * Stopping..."
        PRINT *
        STOP
      ENDIF
      IF( theta(a) < one )THEN
        PRINT *, "** ERROR! The generalized Lorentz factor is ", theta(a), &
                 "< 1 at particle ", a, &
                 "in SUBROUTINE correct_adm_linear_momentum"
        PRINT *, " * Stopping..."
        PRINT *
        STOP
      ENDIF

    ENDDO
    !$OMP END PARALLEL DO

    STOP


  END PROCEDURE correct_adm_linear_momentum


END SUBMODULE finalize_id
