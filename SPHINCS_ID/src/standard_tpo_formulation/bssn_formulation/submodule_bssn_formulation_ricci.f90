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

    USE constants,       ONLY: zero, one
    USE tensor,          ONLY: jx, jy, jz, jxx, jxy, jxz, jyy, jyz, jzz
    USE mesh_refinement, ONLY: allocate_grid_function, levels, nlevels

    IMPLICIT NONE

    INTEGER:: l
    INTEGER, DIMENSION(3) :: imin, imax

    DOUBLE PRECISION:: max_ricci

    ALLOCATE ( levels( THIS% nlevels ), STAT=ios )
    IF( ios > 0 )THEN
     PRINT*,'...allocation error for levels'
     STOP
    ENDIF
    levels = THIS% levels
    nlevels= THIS% nlevels

    CALL allocate_grid_function( THIS% Ricci_ll, "Ricci_ll_id", 6 )
    CALL allocate_grid_function( THIS% Ricci_scalar, "Ricci_scalar_id", 1 )

    DO l= 1, THIS% nlevels, 1

      ASSOCIATE( lapse        => THIS% lapse% levels(l)% var, &
                 shift_u      => THIS% shift_u% levels(l)% var, &
                 phi          => THIS% phi% levels(l)% var, &
                 trK          => THIS% trK% levels(l)% var, &
                 g_BSSN3_ll   => THIS% g_BSSN3_ll% levels(l)% var, &
                 A_BSSN3_ll   => THIS% A_BSSN3_ll% levels(l)% var, &
                 Gamma_u      => THIS% Gamma_u% levels(l)% var, &
                 Ricci_ll     => THIS% Ricci_ll% levels(l)% var, &
                 Ricci_scalar => THIS% Ricci_scalar% levels(l)% var &
      )

        imin(1)= THIS% levels(l)% nghost_x
        imin(2)= THIS% levels(l)% nghost_y
        imin(3)= THIS% levels(l)% nghost_z
        imax(1)= THIS% get_ngrid_x(l) - THIS% levels(l)% nghost_x - 1
        imax(2)= THIS% get_ngrid_y(l) - THIS% levels(l)% nghost_y - 1
        imax(3)= THIS% get_ngrid_z(l) - THIS% levels(l)% nghost_z - 1

        Ricci_ll    = zero
        Ricci_scalar= zero
        CALL bssn_ricci_interior( &
            THIS% get_ngrid_x(l), THIS% get_ngrid_y(l), THIS% get_ngrid_z(l), &
            imin, imax, &
            THIS% get_dx(l), THIS% get_dy(l), THIS% get_dz(l), &
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

    DO l= 1, THIS% nlevels, 1

      ASSOCIATE( Ricci_scalar => THIS% Ricci_scalar% levels(l)% var )

        IF( MAXVAL( ABS(Ricci_scalar) ) > max_ricci )THEN

          max_ricci= MAXVAL( ABS(Ricci_scalar) )

        ENDIF

      END ASSOCIATE

    ENDDO
    PRINT *, max_ricci, SQRT(one/max_ricci)

  END PROCEDURE compute_ricci

END SUBMODULE ricci
