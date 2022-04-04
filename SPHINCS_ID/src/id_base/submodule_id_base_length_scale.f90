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

    USE constants, ONLY: zero, half, two, three, ten, Msun_geo

    IMPLICIT NONE

    INTEGER, PARAMETER:: n= 350
    !! Number of grid points along the shortest size of the matter object

    INTEGER, PARAMETER:: nghost= 4

    !INTEGER:: n_mat
    !! Number of matter objects in the physical system
    INTEGER:: i_mat
    !! Index running over the matter objects
    INTEGER:: i
    INTEGER:: j
    INTEGER:: k
    INTEGER:: ig
    INTEGER:: jg
    INTEGER:: kg

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
    DOUBLE PRECISION:: center(3)
    !! Temporary array to store the centers of the matter objects
    DOUBLE PRECISION:: dx(n_mat)
    !! Uniform spacings of the lattices

    DOUBLE PRECISION:: x
    DOUBLE PRECISION:: y
    DOUBLE PRECISION:: z

    TYPE field
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE:: val
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE:: der
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE:: ratio
    END TYPE field
    TYPE(field), DIMENSION(n_mat):: field_mat

    matter_objects_loop: DO i_mat= 1, n_mat, 1

      sizes= THIS% return_spatial_extent(i_mat)
      center= THIS% return_center(i_mat)
      xL(i_mat)= center(1) - sizes(1)
      xR(i_mat)= center(1) + sizes(2)
      yL(i_mat)= center(2) - sizes(3)
      yR(i_mat)= center(2) + sizes(4)
      zL(i_mat)= center(3) - sizes(5)
      zR(i_mat)= center(3) + sizes(6)

      dx(i_mat)= MINVAL( [(xR(i_mat)-xL(i_mat))/DBLE(n), &
                          (yR(i_mat)-yL(i_mat))/DBLE(n), &
                          (zR(i_mat)-zL(i_mat))/DBLE(n)] )

      PRINT *, " * Lattice step on matter object ", i_mat, "=", &
               dx(i_mat)*Msun_geo*ten*ten*ten, "m"

      nx(i_mat)= NINT( (xR(i_mat)-xL(i_mat))/dx(i_mat) )
      ny(i_mat)= NINT( (yR(i_mat)-yL(i_mat))/dx(i_mat) )
      nz(i_mat)= NINT( (zR(i_mat)-zL(i_mat))/dx(i_mat) )

      ALLOCATE( field_mat(i_mat)% val( nx(i_mat) + 2*nghost, &
                                       ny(i_mat) + 2*nghost, &
                                       nz(i_mat) + 2*nghost ) )
      field_mat(i_mat)% val= zero

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( nx, ny, nz, i_mat, xL, xR, yL, yR, zL, zR, &
      !$OMP                     dx, field_mat ) &
      !$OMP             PRIVATE( i, j, k, x, y, z )
      lattice_loop_x: DO i= 1, nx(i_mat) + 2*nghost, 1

        x= xL(i_mat) + DBLE(i-1)*dx(i_mat)

        lattice_loop_y: DO j= 1, ny(i_mat) + 2*nghost, 1

          y= yL(i_mat) + DBLE(j-1)*dx(i_mat)

          lattice_loop_z: DO k= 1, nz(i_mat) + 2*nghost, 1

            z= zL(i_mat) + DBLE(k-1)*dx(i_mat)

            IF( get_field( x, y, z ) > zero )THEN

              field_mat(i_mat)% val(i,j,k)= get_field( x, y, z )

            ENDIF

          ENDDO lattice_loop_z
        ENDDO lattice_loop_y
      ENDDO lattice_loop_x
      !$OMP END PARALLEL DO

      ! Derivatives
      ALLOCATE( field_mat(i_mat)% der( nx(i_mat), ny(i_mat), nz(i_mat) ) )
      ALLOCATE( field_mat(i_mat)% ratio( nx(i_mat), ny(i_mat), nz(i_mat) ) )

      field_mat(i_mat)% der  = zero
      field_mat(i_mat)% ratio= zero

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( nx, ny, nz, i_mat, xL, xR, yL, yR, zL, zR, &
      !$OMP                     dx, field_mat ) &
      !$OMP             PRIVATE( i, j, k, ig, jg, kg, x, y, z )
      lattice_loop_x_der: DO i= 1, nx(i_mat), 1

        x= xL(i_mat) + DBLE(i-1)*dx(i_mat)
        ig= nghost + i

        lattice_loop_y_der: DO j= 1, ny(i_mat), 1

          y= yL(i_mat) + DBLE(j-1)*dx(i_mat)
          jg= nghost + j

          lattice_loop_z_der: DO k= 1, nz(i_mat), 1

            z= zL(i_mat) + DBLE(k-1)*dx(i_mat)
            kg= nghost + k

            field_mat(i_mat)% der(i,j,k)= &
  !( field_mat(i_mat)% val(ig+1,jg,kg) - field_mat(i_mat)% val(ig-1,jg,kg) )*half
   !         - field_mat(i_mat)% val(ig+2,jg,kg)/(two*two*three) &
   !         + field_mat(i_mat)% val(ig+1,jg,kg)*two/three &
   !         - field_mat(i_mat)% val(ig-1,jg,kg)*two/three &
   !         + field_mat(i_mat)% val(ig-2,jg,kg)/(two*two*three) &
            ( field_mat(i_mat)% val(ig+3,jg,kg)/(two*three*ten) &
            - field_mat(i_mat)% val(ig+2,jg,kg)*three/(two*ten) &
            + field_mat(i_mat)% val(ig+1,jg,kg)*three/(two*two) &
            - field_mat(i_mat)% val(ig-1,jg,kg)*three/(two*two) &
            + field_mat(i_mat)% val(ig-2,jg,kg)*three/(two*ten) &
            - field_mat(i_mat)% val(ig-3,jg,kg)/(two*three*ten) )/dx(i_mat)

            IF( field_mat(i_mat)% der(i,j,k) /= zero )THEN

              field_mat(i_mat)% ratio(i,j,k)= &
                field_mat(i_mat)% val(i,j,k)/ABS(field_mat(i_mat)% der(i,j,k))

            ENDIF

          ENDDO lattice_loop_z_der
        ENDDO lattice_loop_y_der
      ENDDO lattice_loop_x_der
      !$OMP END PARALLEL DO

      scales(i_mat)= MINVAL( field_mat(i_mat)% ratio, &
                             field_mat(i_mat)% ratio > 0 )

      !PRINT *, SUM( field_mat(i_mat)% ratio )/(nx(i_mat)*ny(i_mat)*nz(i_mat))
      !STOP

      DEALLOCATE( field_mat(i_mat)% val )
      DEALLOCATE( field_mat(i_mat)% der )
      DEALLOCATE( field_mat(i_mat)% ratio )

    ENDDO matter_objects_loop
    PRINT *

  END PROCEDURE estimate_lengthscale_field


END SUBMODULE length_scale
