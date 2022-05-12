! File:         submodule_bssn_formulation_ricci.f90
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

SUBMODULE (bssn_formulation) ricci

  !************************************************
  !
  !# Implementation of the method of TYPE bssn
  !  that computes the Ricci tensor and scalar
  !
  !  FT 10.02.2022
  !
  !************************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE compute_ricci

    !************************************************
    !
    !# Computes the Ricci tensor and the Ricci
    !  scalar on the mesh
    !
    !  FT 10.02.2022
    !
    !************************************************

    USE constants,       ONLY: zero, one, two, three, ten, Msun_geo
    USE tensor,          ONLY: jx, jy, jz, jxx, jxy, jxz, jyy, jyz, jzz
    USE mesh_refinement, ONLY: allocate_grid_function, levels, nlevels

    IMPLICIT NONE

    INTEGER:: l
    INTEGER, DIMENSION(3) :: imin, imax

    DOUBLE PRECISION:: max_ricci, min_lapse

    ALLOCATE ( levels( this% nlevels ), STAT=ios )
    IF( ios > 0 )THEN
     PRINT*,'...allocation error for levels'
     STOP
    ENDIF
    levels = this% levels
    nlevels= this% nlevels

    CALL allocate_grid_function( this% Ricci_ll, "Ricci_ll_id", 6 )
    CALL allocate_grid_function( this% Ricci_scalar, "Ricci_scalar_id", 1 )

    DO l= 1, this% nlevels, 1

      ASSOCIATE( lapse        => this% lapse% levels(l)% var, &
                 shift_u      => this% shift_u% levels(l)% var, &
                 phi          => this% phi% levels(l)% var, &
                 trK          => this% trK% levels(l)% var, &
                 g_BSSN3_ll   => this% g_BSSN3_ll% levels(l)% var, &
                 A_BSSN3_ll   => this% A_BSSN3_ll% levels(l)% var, &
                 Gamma_u      => this% Gamma_u% levels(l)% var, &
                 Ricci_ll     => this% Ricci_ll% levels(l)% var, &
                 Ricci_scalar => this% Ricci_scalar% levels(l)% var &
      )

        imin(1)= this% levels(l)% nghost_x
        imin(2)= this% levels(l)% nghost_y
        imin(3)= this% levels(l)% nghost_z
        imax(1)= this% get_ngrid_x(l) - this% levels(l)% nghost_x - 1
        imax(2)= this% get_ngrid_y(l) - this% levels(l)% nghost_y - 1
        imax(3)= this% get_ngrid_z(l) - this% levels(l)% nghost_z - 1

        Ricci_ll    = zero
        Ricci_scalar= zero
        CALL bssn_ricci_interior( &
            this% get_ngrid_x(l), this% get_ngrid_y(l), this% get_ngrid_z(l), &
            imin, imax, &
            this% get_dx(l), this% get_dy(l), this% get_dz(l), &
            g_BSSN3_ll(:,:,:,jxx), g_BSSN3_ll(:,:,:,jxy), &
            g_BSSN3_ll(:,:,:,jxz), g_BSSN3_ll(:,:,:,jyy), &
            g_BSSN3_ll(:,:,:,jyz), g_BSSN3_ll(:,:,:,jzz), &
            A_BSSN3_ll(:,:,:,jxx), A_BSSN3_ll(:,:,:,jxy), &
            A_BSSN3_ll(:,:,:,jxz), A_BSSN3_ll(:,:,:,jyy), &
            A_BSSN3_ll(:,:,:,jyz), A_BSSN3_ll(:,:,:,jzz), &
            trK(:,:,:), phi(:,:,:), &
            Gamma_u(:,:,:,jx), &
            Gamma_u(:,:,:,jy), &
            Gamma_u(:,:,:,jz), &
            Ricci_ll(:,:,:,jxx), Ricci_ll(:,:,:,jxy), Ricci_ll(:,:,:,jxz), &
            Ricci_ll(:,:,:,jyy), Ricci_ll(:,:,:,jyz), Ricci_ll(:,:,:,jzz), &
            Ricci_scalar(:,:,:) )

      END ASSOCIATE

    ENDDO
    PRINT *, " * Ricci tensor and scalar computed."
    PRINT *

    ! From Wikipedia:
    ! The sectional curvature of an n-sphere of radius r is K = 1/r^2.
    ! Hence the scalar curvature is R = n(n − 1)/r^2.
    ! Find real reference for this.
    max_ricci= zero
    min_lapse= HUGE(one)
    ASSOCIATE( Ricci_scalar => this% Ricci_scalar% levels(this% nlevels)% var, &
               lapse        => this% lapse% levels(this% nlevels)% var )

      !IF( MAXVAL( ABS(Ricci_scalar) ) > max_ricci )THEN

        max_ricci= MAXVAL( ABS(Ricci_scalar) )

      !ENDIF

      !IF( MINVAL( ABS(lapse) ) < min_lapse )THEN
      !
      !  min_lapse= MINVAL( ABS(lapse) )
      !
      !ENDIF

    END ASSOCIATE
    PRINT *, "* Maximum value of the Ricci scalar over the finest refinement ", &
             "level= ", max_ricci, "Msun_geo^{-2}"
    PRINT *, "* Radius of a 3-sphere with the same Ricci scalar=", &
             SQRT(three*two/max_ricci), "Msun_geo=", &
             SQRT(three*two/max_ricci)*Msun_geo*ten*ten*ten, "m"
    !PRINT *, min_lapse

  END PROCEDURE compute_ricci


END SUBMODULE ricci
