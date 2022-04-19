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

    USE constants,  ONLY: zero, one, two, amu, Msun
    USE tensor,     ONLY: jx, jy, jz, n_sym4x4, itt, itx, ity, itz, ixx, ixy, &
                          ixz, iyy, iyz, izz, raise_index_4vector,  &
                          lower_index_4vector
    USE matrix,     ONLY: invert_3x3_matrix, invert_4x4_matrix
    USE utility,    ONLY: compute_g4, spacetime_vector_norm_sym4x4, &
                          spatial_vector_norm_sym3x3, is_finite_number, &
                          determinant_sym4x4

    IMPLICIT NONE

    INTEGER:: a, j

    DOUBLE PRECISION:: shift_norm2, shift_delta
    DOUBLE PRECISION, DIMENSION(3,npart):: vel_l_corr
    DOUBLE PRECISION, DIMENSION(3):: delta!, com_vel!, adm_mom
    DOUBLE PRECISION:: g4(n_sym4x4), det, den
    DOUBLE PRECISION, DIMENSION(3,3):: g3mat
    DOUBLE PRECISION, DIMENSION(3,3):: g3mat_inv
    DOUBLE PRECISION, DIMENSION(4,4):: g4mat
    DOUBLE PRECISION, DIMENSION(4,4):: g4mat_inv
    DOUBLE PRECISION, DIMENSION(0:3):: v_l
    DOUBLE PRECISION, DIMENSION(0:3):: v_u

    !DOUBLE PRECISION, DIMENSION(3):: com_vel


    !
    !-- Simply subtract the velocity of the center of mass from the velocity
    !-- of the particles, and recompute the generalized Lorentz factor theta
    !-- It doesn't work
    !

!    com_vel= adm_mom_error/adm_mass
!    PRINT *, "com_vel= ", com_vel
!
!    !$OMP PARALLEL DO DEFAULT( NONE ) &
!    !$OMP             SHARED( npart, lapse, shift_x, shift_y, shift_z, vel_u, &
!    !$OMP                     g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, &
!    !$OMP                     theta, adm_mom_error, com_vel ) &
!    !$OMP             PRIVATE( a, g4, det )
!    DO a= 1, npart, 1
!
!      vel_u(:,a)= vel_u(:,a) - com_vel
!
!      CALL compute_g4( lapse(a), [shift_x(a),shift_y(a),shift_z(a)], &
!                       [g_xx(a),g_xy(a),g_xz(a),g_yy(a),g_yz(a),g_zz(a)], g4 )
!
!      CALL determinant_sym4x4( g4, det )
!      IF( ABS(det) < 1D-10 )THEN
!        PRINT *, "** ERROR! The determinant of the spacetime metric is " &
!                 // "effectively 0 at particle ", a
!        PRINT *, " * det= ", det, &
!                 "in SUBROUTINE correct_adm_linear_momentum"
!        PRINT *, " * Stopping..."
!        STOP
!      ELSEIF( det > 0 )THEN
!        PRINT *, "** ERROR! The determinant of the spacetime metric is " &
!                 // "positive at particle ", a
!        PRINT *, " * det= ", det, &
!                 "in SUBROUTINE correct_adm_linear_momentum"
!        PRINT *, " * Stopping..."
!        STOP
!      ENDIF
!
!      CALL spacetime_vector_norm_sym4x4( g4, [one, vel_u(:,a)], theta(a) )
!      IF( .NOT.is_finite_number(theta(a)) )THEN
!        PRINT *, "** ERROR! The spacetime norm of vel_u is ", theta(a), &
!                 "at particle ", a, &
!                 "in SUBROUTINE correct_adm_linear_momentum"
!        PRINT *, " * Stopping..."
!        PRINT *
!        STOP
!      ENDIF
!      IF( theta(a) > zero )THEN
!        PRINT *, "** ERROR! The spacetime norm of vel_u is ", theta(a), &
!                 "at particle ", a, "(that is, vel_u is spacelike)", &
!                 "in SUBROUTINE correct_adm_linear_momentum"
!        PRINT *, " * Stopping..."
!        PRINT *
!        STOP
!      ENDIF
!      IF( theta(a) == zero )THEN
!        PRINT *, "** ERROR! The spacetime norm of vel_u is ", theta(a), &
!                 "at particle ", a, "(that is, vel_u is null)", &
!                 "in SUBROUTINE correct_adm_linear_momentum"
!        PRINT *, " * Stopping..."
!        PRINT *
!        STOP
!      ENDIF
!
!      theta(a)= one/SQRT(-theta(a))
!      IF( .NOT.is_finite_number(theta(a)) )THEN
!        PRINT *, "** ERROR! The generalized Lorentz factor is ", theta(a), &
!                 "at particle ", a, &
!                 "in SUBROUTINE correct_adm_linear_momentum"
!        PRINT *, " * Stopping..."
!        PRINT *
!        STOP
!      ENDIF
!      IF( theta(a) < one )THEN
!        PRINT *, "** ERROR! The generalized Lorentz factor is ", theta(a), &
!                 "< 1 at particle ", a, &
!                 "in SUBROUTINE correct_adm_linear_momentum"
!        PRINT *, " * Stopping..."
!        PRINT *
!        STOP
!      ENDIF
!
!    ENDDO
!    !$OMP END PARALLEL DO


 !   adm_mom= zero
 !   !$OMP PARALLEL DO DEFAULT(NONE) &
 !   !$OMP             SHARED( npart, nu, lapse, shift_x, shift_y, shift_z, &
 !   !$OMP                     theta, u, pr, nlrf, vel_u, &
 !   !$OMP                     g_xx, g_xy, g_xz, g_yy, g_yz, g_zz ) &
 !   !$OMP             PRIVATE( a, det, v_u, shift_norm2, j, g4, v_l ) &
 !   !$OMP             REDUCTION( +: adm_mom )
 !   DO a= 1, npart, 1
 !
 !     CALL compute_g4( lapse(a), [shift_x(a),shift_y(a),shift_z(a)], &
 !                      [g_xx(a),g_xy(a),g_xz(a),g_yy(a),g_yz(a),g_zz(a)], &
 !                      g4 )
 !
 !     CALL determinant_sym4x4( g4, det )
 !     IF( ABS(det) < 1D-10 )THEN
 !       PRINT *, "** ERROR! The determinant of the spacetime metric is " &
 !                // "effectively 0 at particle ", a
 !       PRINT *, " * det= ", det
 !       PRINT *, " * Stopping..."
 !       STOP
 !     ELSEIF( det > 0 )THEN
 !       PRINT *, "** ERROR! The determinant of the spacetime metric is " &
 !                // "positive at particle ", a
 !       PRINT *, " * det= ", det
 !       PRINT *, " * Stopping..."
 !       STOP
 !     ENDIF
 !
 !     v_u= [one, vel_u(:,a)]
 !     CALL lower_index_4vector( v_u, g4, v_l )
 !
 !     CALL spatial_vector_norm_sym3x3( &
 !                   [g_xx(a),g_xy(a),g_xz(a),g_yy(a),g_yz(a),g_zz(a)], &
 !                   [shift_x(a),shift_y(a),shift_z(a)], shift_norm2 )
 !
 !     DO j= jx, jz, 1
 !
 !       adm_mom(j)= adm_mom(j) &
 !           - ( nu(a)*amu/Msun )*( shift_norm2/(lapse(a)**two) - one ) &
 !             *theta(a)*( one + u(a) + pr(a)/nlrf(a) )*v_l(j)
 !
 !     ENDDO
 !
 !     !IF( a == 1 )THEN
 !     !  PRINT *, "v_l inside compute_adm= ", v_l
 !     !  PRINT *
 !     !ENDIF
 !
 !   ENDDO
 !   !$OMP END PARALLEL DO
 !
 !   PRINT *, "adm_mom before correction=", adm_mom
 !   PRINT *

!===============================================================================

    den= zero
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP             SHARED( npart, nu, lapse, shift_x, shift_y, shift_z, &
    !$OMP                     theta, u, pr, nlrf, vel_u, &
    !$OMP                     v_l, g_xx, g_xy, g_xz, g_yy, g_yz, g_zz ) &
    !$OMP             PRIVATE( a, det, v_u, shift_norm2, j, g4 ) &
    !$OMP             REDUCTION( +: den )
    DO a= 1, npart, 1

      CALL spatial_vector_norm_sym3x3( &
                    [g_xx(a),g_xy(a),g_xz(a),g_yy(a),g_yz(a),g_zz(a)], &
                    [shift_x(a),shift_y(a),shift_z(a)], shift_norm2 )

      den= den - ( nu(a)*amu/Msun )*theta(a) &
                *( shift_norm2/(lapse(a)**two) - one ) &
                *( one + u(a) + pr(a)/nlrf(a) )

    ENDDO
    !$OMP END PARALLEL DO

    delta(1)= - adm_mom_error(1)/den
    delta(2)= - adm_mom_error(2)/den
    delta(3)= - adm_mom_error(3)/den

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( npart, lapse, shift_x, shift_y, shift_z, vel_u, &
    !$OMP                     vel_l_corr, g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, &
    !$OMP                     theta, adm_mom_error, nu, nstar, nlrf, delta ) &
    !$OMP             PRIVATE( a, j, shift_norm2, g4, g3mat, g3mat_inv, det, &
    !$OMP                      v_l, v_u, g4mat, g4mat_inv, shift_delta )
    DO a= 1, npart, 1

      CALL spatial_vector_norm_sym3x3( &
                    [g_xx(a),g_xy(a),g_xz(a),g_yy(a),g_yz(a),g_zz(a)], &
                    [shift_x(a),shift_y(a),shift_z(a)], shift_norm2 )

      CALL compute_g4( lapse(a), [shift_x(a),shift_y(a),shift_z(a)], &
                       [g_xx(a),g_xy(a),g_xz(a),g_yy(a),g_yz(a),g_zz(a)], g4 )

      CALL determinant_sym4x4( g4, det )
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

      !nu(a)= nu(a)/theta(a)

      g3mat(1,1)= g_xx(a)
      g3mat(1,2)= g_xy(a)
      g3mat(1,3)= g_xz(a)
      g3mat(2,1)= g_xy(a)
      g3mat(2,2)= g_yy(a)
      g3mat(2,3)= g_yz(a)
      g3mat(3,1)= g_xz(a)
      g3mat(3,2)= g_yz(a)
      g3mat(3,3)= g_zz(a)

      g4mat(1,1)= g4(itt)
      g4mat(1,2)= g4(itx)
      g4mat(1,3)= g4(ity)
      g4mat(1,4)= g4(itz)
      g4mat(2,1)= g4(itx)
      g4mat(2,2)= g4(ixx)
      g4mat(2,3)= g4(ixy)
      g4mat(2,4)= g4(ixz)
      g4mat(3,1)= g4(ity)
      g4mat(3,2)= g4(ixy)
      g4mat(3,3)= g4(iyy)
      g4mat(3,4)= g4(iyz)
      g4mat(4,1)= g4(itz)
      g4mat(4,2)= g4(ixz)
      g4mat(4,3)= g4(iyz)
      g4mat(4,4)= g4(izz)

      CALL invert_3x3_matrix( g3mat, g3mat_inv )
      CALL invert_4x4_matrix( g4mat, g4mat_inv )

      CALL lower_index_4vector( [one, vel_u(:,a)], g4, v_l )

      !shift_delta= shift_x(a)*delta(jx) + shift_y(a)*delta(jy) &
      !             + shift_z(a)*delta(jz)

      !vel_l_corr(jx,a)= delta(jx) &
      !               /( two*g3mat_inv(1,1)*NORM2(delta)**two ) &
      !               *( - ( one + shift_delta ) &
      !                  + SQRT( ( one + shift_delta )**two &
      !                  + two*two*g3mat_inv(1,1)*NORM2(delta)**two &
      !                    *lapse(a)**two ) )

      !vel_l_corr(jx,a)= delta(jx) &
      !               /( two*( one + g3mat_inv(1,1)*NORM2(delta)**two ) ) &
      !               *( - shift_delta &
      !                  + SQRT( shift_delta**two &
      !                  + two*two*( one + g3mat_inv(1,1)*NORM2(delta)**two ) &
      !                    *lapse(a)**two ) )

      vel_l_corr(jx,a)= delta(jx)

      vel_l_corr(jy,a)= vel_l_corr(jx,a)*delta(jy)/delta(jx)
      vel_l_corr(jz,a)= vel_l_corr(jx,a)*delta(jz)/delta(jx)

      DO j= jx, jz, 1

        v_l(j)= v_l(j) + vel_l_corr(j,a)

      ENDDO

      v_l(0)= ( one - g4mat_inv(1,2)*v_l(jx) - g4mat_inv(1,3)*v_l(jy) &
                - g4mat_inv(1,4)*v_l(jz) )/g4mat_inv(1,1)

      !IF( a == 1 )THEN
      !  PRINT *, "v_l inside correct_adm= ", v_l
      !  PRINT *
      !ENDIF

      CALL raise_index_4vector( v_l, &
              [g4mat_inv(1,1),g4mat_inv(1,2),g4mat_inv(1,3), &
               g4mat_inv(1,4),g4mat_inv(2,2),g4mat_inv(2,3), &
               g4mat_inv(2,4),g4mat_inv(3,3),g4mat_inv(3,4),g4mat_inv(4,4)], &
               v_u )

      IF( ABS( v_u(0) - one ) > 1.D-10 )THEN
        PRINT *, "** ERROR! The 0 component of the corrected computing frame ",&
                 "velocity at particle ", a, "is not 1 ", &
                 "in SUBROUTINE correct_adm_linear_momentum."
        PRINT *, " * v(0)= ", v_u(0)
        PRINT *, " * Stopping..."
        PRINT *
        STOP
      ENDIF

      vel_u(:,a)= v_u(1:3)

   !   CALL spacetime_vector_norm_sym4x4( g4, v_u, theta(a) )
   !   IF( .NOT.is_finite_number(theta(a)) )THEN
   !     PRINT *, "** ERROR! The spacetime norm of vel_u is ", theta(a), &
   !              "at particle ", a, &
   !              "in SUBROUTINE correct_adm_linear_momentum"
   !     PRINT *, " * Stopping..."
   !     PRINT *
   !     STOP
   !   ENDIF
   !   IF( theta(a) > zero )THEN
   !     PRINT *, "** ERROR! The spacetime norm of vel_u is ", theta(a), &
   !              "at particle ", a, "(that is, vel_u is spacelike) ", &
   !              "in SUBROUTINE correct_adm_linear_momentum"
   !     PRINT *, " * Stopping..."
   !     PRINT *
   !     STOP
   !   ENDIF
   !   IF( theta(a) == zero )THEN
   !     PRINT *, "** ERROR! The spacetime norm of vel_u is ", theta(a), &
   !              "at particle ", a, "(that is, vel_u is null) ", &
   !              "in SUBROUTINE correct_adm_linear_momentum"
   !     PRINT *, " * Stopping..."
   !     PRINT *
   !     STOP
   !   ENDIF
   !
   !   theta(a)= one/SQRT(-theta(a))
   !   IF( .NOT.is_finite_number(theta(a)) )THEN
   !     PRINT *, "** ERROR! The generalized Lorentz factor is ", theta(a), &
   !              "at particle ", a, &
   !              "in SUBROUTINE correct_adm_linear_momentum"
   !     PRINT *, " * Stopping..."
   !     PRINT *
   !     STOP
   !   ENDIF
   !   IF( theta(a) < one )THEN
   !     PRINT *, "** ERROR! The generalized Lorentz factor is ", theta(a), &
   !              "< 1 at particle ", a, &
   !              "in SUBROUTINE correct_adm_linear_momentum"
   !     PRINT *, " * Stopping..."
   !     PRINT *
   !     STOP
   !   ENDIF

      !nstar(a)= nlrf(a)*SQRT(-det)*theta(a)
      !nu(a)= nu(a)*theta(a)

    ENDDO
    !$OMP END PARALLEL DO

!===============================================================================

 !   adm_mom= zero
 !   !$OMP PARALLEL DO DEFAULT(NONE) &
 !   !$OMP             SHARED( npart, nu, lapse, shift_x, shift_y, shift_z, &
 !   !$OMP                     theta, u, pr, nlrf, vel_u, &
 !   !$OMP                     g_xx, g_xy, g_xz, g_yy, g_yz, g_zz ) &
 !   !$OMP             PRIVATE( a, det, v_u, shift_norm2, j, g4, v_l ) &
 !   !$OMP             REDUCTION( +: adm_mom )
 !   DO a= 1, npart, 1
 !
 !     CALL compute_g4( lapse(a), [shift_x(a),shift_y(a),shift_z(a)], &
 !                      [g_xx(a),g_xy(a),g_xz(a),g_yy(a),g_yz(a),g_zz(a)], &
 !                      g4 )
 !
 !     CALL determinant_sym4x4( g4, det )
 !     IF( ABS(det) < 1D-10 )THEN
 !       PRINT *, "** ERROR! The determinant of the spacetime metric is " &
 !                // "effectively 0 at particle ", a
 !       PRINT *, " * det= ", det
 !       PRINT *, " * Stopping..."
 !       STOP
 !     ELSEIF( det > 0 )THEN
 !       PRINT *, "** ERROR! The determinant of the spacetime metric is " &
 !                // "positive at particle ", a
 !       PRINT *, " * det= ", det
 !       PRINT *, " * Stopping..."
 !       STOP
 !     ENDIF
 !
 !     v_u= [one, vel_u(:,a)]
 !     CALL lower_index_4vector( v_u, g4, v_l )
 !
 !     CALL spatial_vector_norm_sym3x3( &
 !                   [g_xx(a),g_xy(a),g_xz(a),g_yy(a),g_yz(a),g_zz(a)], &
 !                   [shift_x(a),shift_y(a),shift_z(a)], shift_norm2 )
 !
 !     DO j= jx, jz, 1
 !
 !       adm_mom(j)= adm_mom(j) &
 !           - ( nu(a)*amu/Msun )*( shift_norm2/(lapse(a)**two) - one ) &
 !             *theta(a)*( one + u(a) + pr(a)/nlrf(a) )*v_l(j)
 !
 !     ENDDO
 !
 !     !IF( a == 1 )THEN
 !     !  PRINT *, "v_l inside compute_adm= ", v_l
 !     !  PRINT *
 !     !ENDIF
 !
 !   ENDDO
 !   !$OMP END PARALLEL DO
 !
 !   PRINT *, "adm_mom after correction=", adm_mom
 !   PRINT *
 !   !STOP


  END PROCEDURE correct_adm_linear_momentum


END SUBMODULE finalize_id
