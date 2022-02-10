! File:         submodule_id_base_length_scale.f90
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

SUBMODULE (id_base) length_scale

  !********************************************
  !
  !# Implementation of the method of TYPE idbase
  !  that estimates typical length scales, one
  !  per each matter object, by computing
  !  \(\dfrac{f}{\partial f}\), where \(f\) is a
  !  field given as input, and \(\partial\)
  !  represent a derivative of it.
  !
  !  FT 10.02.2022
  !
  !********************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE estimate_lengthscale_field

    !************************************************
    !
    !# Estimate typical length scales, one per each
    !  matter object, by computing \(\dfrac{f}{\partial f}\),
    !  where \(f\) is a field given as input, and \(\partial\)
    !  represent a derivative of it.
    !  Presently, the derivatives are computed separately
    !  along each spatial dimension, as 1D derivatives.
    !
    !  FT 10.02.2022
    !
    !************************************************


    IMPLICIT NONE

    INTEGER, PARAMETER:: n= 200
    !! Number of grid points along the shortest size of the matter object

    !INTEGER:: n_mat
    !! Number of matter objects in the physical system
    INTEGER:: i_mat
    !! Index running over the matter objects
    INTEGER:: i
    INTEGER:: j
    INTEGER:: k

    INTEGER:: nx(n_mat)
    !!
    INTEGER:: ny(n_mat)
    !!
    INTEGER:: nz(n_mat)
    !!

    DOUBLE PRECISION:: xL(n_mat)
    !! Left boundaries of the lattices in the \(x\) direction
    DOUBLE PRECISION:: xR(n_mat)
    !! Right boundaries of the lattices in the \(x\) direction
    DOUBLE PRECISION:: yL(n_mat)
    !! Left boundaries of the lattices in the \(y\) direction
    DOUBLE PRECISION:: yR(n_mat)
    !! Right boundaries of the lattices in the \(y\) direction
    DOUBLE PRECISION:: zL(n_mat)
    !! Left boundaries of the lattices in the \(z\) direction
    DOUBLE PRECISION:: zR(n_mat)
    !! Right boundaries of the lattices in the \(z\) direction
    DOUBLE PRECISION:: sizes(6)
    !! Temporary array to store the sizes of the matter objects
    DOUBLE PRECISION:: dx(n_mat)
    !! Uniform spacings of the lattices

    DOUBLE PRECISION:: x
    DOUBLE PRECISION:: y
    DOUBLE PRECISION:: z

    TYPE field
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE:: val
    END TYPE field
    TYPE(field), DIMENSION(n_mat):: field_mat

    matter_objects_loop: DO i_mat= 1, n_mat, 1

      sizes= THIS% return_spatial_extent(i_mat)
      xL(i_mat)= sizes(1)
      xR(i_mat)= sizes(2)
      yL(i_mat)= sizes(3)
      yR(i_mat)= sizes(4)
      zL(i_mat)= sizes(5)
      zR(i_mat)= sizes(6)

      dx(i_mat)= MINVAL( [(xR(i_mat)-xL(i_mat))/DBLE(n), &
                          (yR(i_mat)-yL(i_mat))/DBLE(n), &
                          (zR(i_mat)-zL(i_mat))/DBLE(n)] )

      nx(i_mat)= xR(i_mat)-xL(i_mat)/dx(i_mat)
      ny(i_mat)= yR(i_mat)-yL(i_mat)/dx(i_mat)
      nz(i_mat)= zR(i_mat)-zL(i_mat)/dx(i_mat)

      ALLOCATE( field_mat(i_mat)% val(nx(i_mat), ny(i_mat), nz(i_mat)) )

      lattice_loop: DO i= 1, n, 1

        x= xL(i_mat) + DBLE(i)*dx(i_mat)

        DO j= 1, n, 1

          y= yL(i_mat) + DBLE(j)*dx(i_mat)

          DO k= 1, n, 1

            z= zL(i_mat) + DBLE(z)*dx(i_mat)

            field_mat% val(i,j,k)= get_field( x, y, z )

          ENDDO
        ENDDO
      END DO lattice_loop

      DEALLOCATE( field_mat(i_mat)% val )

    ENDDO matter_objects_loop


  END PROCEDURE estimate_lengthscale_field


END SUBMODULE length_scale
