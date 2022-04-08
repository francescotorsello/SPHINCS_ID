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


  MODULE PROCEDURE compute_adm_momentum

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

    USE units,                ONLY: m0c2_cu
    USE recovery,             ONLY: phys_2_cons, cons_2_phys
    USE tensor,               ONLY: jx, jy, jz, n_sym4x4
    USE constants,            ONLY: zero, one, two, amu, MSun
    USE deactivate_particles, ONLY: nlrf_fb, u_fb, pr_fb, vel_u_fb, theta_fb, &
                                    cs_fb
    USE metric_on_particles,  ONLY: allocate_metric_on_particles, &
                                    deallocate_metric_on_particles, &
                                    g4_ll
    USE utility,              ONLY: compute_g4, determinant_sym4x4, &
                                    spatial_vector_norm_sym3x3
    !USE tmp,                  ONLY: fill_arrays

    IMPLICIT NONE

    INTEGER, PARAMETER:: unit_recovery= 34956

    INTEGER:: i_matter, a, a_max, j

    DOUBLE PRECISION:: det, p_max, shift_norm2

    DOUBLE PRECISION, DIMENSION(npart)  :: nlrf_rec
    DOUBLE PRECISION, DIMENSION(npart)  :: u_rec
    DOUBLE PRECISION, DIMENSION(npart)  :: pr_rec
    DOUBLE PRECISION, DIMENSION(3,npart):: vel_u_rec
    DOUBLE PRECISION, DIMENSION(npart)  :: theta_rec
    DOUBLE PRECISION, DIMENSION(npart)  :: nstar_rec
    DOUBLE PRECISION, DIMENSION(3,npart):: s_l_rec
    DOUBLE PRECISION, DIMENSION(npart)  :: e_hat_rec

    LOGICAL:: exist

    CHARACTER( LEN= 2 ):: i_mat
    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    LOGICAL, PARAMETER:: debug= .FALSE.

    PRINT *, " * Estimating ADM linear momentum using the canonical SPH ", &
             "momentum per baryon on the particles..."
    PRINT *

    IF( .NOT.ALLOCATED(g4_ll) )THEN
      CALL allocate_metric_on_particles(npart)
    ENDIF

    DO a= 1, npart, 1

      CALL compute_g4( lapse(a), shift(:,a), &
            [this% g_xx(a), this% g_xy(a), this% g_xz(a), &
             this% g_yy(a), this% g_yz(a), this% g_zz(a)], &
             g4_ll(1:n_sym4x4,a) )

      CALL determinant_sym4x4( g4_ll(1:n_sym4x4,a), det )
      IF( ABS(det) < 1D-10 )THEN
          PRINT *, "** ERROR! The determinant of the spacetime metric is " &
                   // "effectively 0 at particle ", a
          STOP
      ELSEIF( det > 0 )THEN
          PRINT *, "** ERROR! The determinant of the spacetime metric is " &
                   // "positive at particle ", a
          STOP
      ENDIF

    ENDDO

    ALLOCATE( nlrf_fb (npart) )
    ALLOCATE( u_fb    (npart) )
    ALLOCATE( pr_fb   (npart) )
    ALLOCATE( vel_u_fb(3,npart) )
    ALLOCATE( theta_fb(npart) )
    ALLOCATE( cs_fb(npart) )
    nlrf_fb = nlrf
    u_fb    = u
    pr_fb   = pr
    vel_u_fb= vel_u(1:3,:)
    theta_fb= theta
    ! TODO: set the sound speed properly. From pwp_eos MODULE:
    ! enth= 1.0D0 + u + rho_rest/P
    ! cs=   SQRT((Gamma*P_cold + Gamma_th*P_th)/(rho_rest*enth))
    cs_fb   = one


    ! Initialize local arrays
    nlrf_rec = zero
    u_rec    = zero
    pr_rec   = zero
    vel_u_rec= zero
    theta_rec= zero
    nstar_rec= zero
    s_l_rec  = zero
    e_hat_rec= zero

    IF( debug ) PRINT *, "1"

    !
    !-- Compute conserved fields from physical fields
    !
    CALL phys_2_cons( npart, nlrf, u, pr, vel_u, &
                      ! following is output
                      nstar_rec, s_l_rec, e_hat_rec )

    IF( debug ) PRINT *, "2"

    !
    !-- Compute the ADM linear momentum
    !
    adm_mom= zero
    !$OMP PARALLEL DO SHARED( npart, nu, lapse, shift, s_l_rec, this ) &
    !$OMP             PRIVATE( a, shift_norm2 ) &
    !$OMP             REDUCTION( +: adm_mom )
    DO a= 1, npart, 1

      DO j= 1, 3, 1

        CALL spatial_vector_norm_sym3x3( &
                        [this% g_xx(a), this% g_xy(a), this% g_xz(a), &
                         this% g_yy(a), this% g_yz(a), this% g_zz(a)], &
                        shift(:,a), shift_norm2 )

        adm_mom(j)= adm_mom(j) &
            - nu(a)*( shift_norm2/(lapse(a)**two) - one )*s_l_rec(j,a)*amu/Msun

      ENDDO

    ENDDO
    !$OMP END PARALLEL DO

    IF( debug ) PRINT *, "3"

    DEALLOCATE( nlrf_fb )
    DEALLOCATE( u_fb )
    DEALLOCATE( pr_fb )
    DEALLOCATE( vel_u_fb )
    DEALLOCATE( theta_fb )
    DEALLOCATE( cs_fb )

    CALL deallocate_metric_on_particles

    IF( debug ) PRINT *, "4"

  END PROCEDURE compute_adm_momentum


END SUBMODULE adm_variables
