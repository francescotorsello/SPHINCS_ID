! File:         submodule_bns_base_find_radii.f90
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

SUBMODULE(bns_base) find_radii

  !********************************************************
  !
  !# This MODULE contains the implementation of the 
  !  methods of TYPE bns_base that deal with the
  !  computation of the radii of the stars.
  !
  !  FT 27.09.2022
  !
  !********************************************************


  USE utility, ONLY: zero, one, two, three, four, five, ten


  IMPLICIT NONE


  DOUBLE PRECISION, PARAMETER:: tol= 1.D-5



  CONTAINS



  MODULE PROCEDURE find_radius

    !********************************************************
    !
    !# 
    !
    !  FT 27.09.2022
    !
    !********************************************************

    INTEGER,          PARAMETER:: n_pts   = 1.5D+3
    DOUBLE PRECISION, PARAMETER:: min_dist= 1.D-2
    DOUBLE PRECISION, PARAMETER:: x_ini   = 130.D0
    DOUBLE PRECISION, PARAMETER:: min_dens= 1.D-8
    LOGICAL,          PARAMETER:: debug= .FALSE.

    INTEGER:: i
    DOUBLE PRECISION:: vector_norm, x
    DOUBLE PRECISION, DIMENSION(3):: point_left, point_right, point, versor, &
                                     point_mean
    DOUBLE PRECISION:: rho, rho_left, rho_right, x_left, x_right, x_mean, &
                       rho_mean
  

    !INTERFACE
    !  FUNCTION get_density_at_pos(x, y, z) RESULT(rho)
    !    DOUBLE PRECISION, INTENT(IN):: x
    !    DOUBLE PRECISION, INTENT(IN):: y
    !    DOUBLE PRECISION, INTENT(IN):: z
    !    DOUBLE PRECISION:: rho
    !  END FUNCTION get_density_at_pos
    !END INTERFACE
    PROCEDURE(), POINTER:: return_density

    IF(PRESENT(get_density))THEN

       return_density => read_density_opt
      
    ELSE

       return_density => read_density

    ENDIF

    radius= 1.45D0

    vector_norm= NORM2(vector - center)
    versor= vector/vector_norm

    x_left = one
    x_right= x_ini
    DO

      point_left = line(x_left)
      point_right= line(x_right)
      CALL return_density(point_left(1), point_left(2), point_left(3), rho_left)
      CALL return_density(point_right(1), point_right(2), point_right(3), &
                          rho_right)

      IF( rho_left > min_dens .AND. rho_right <= min_dens )THEN

        IF( NORM2(point_left - point_right) < min_dist )THEN
          EXIT
        ENDIF
        x_mean= (x_left + x_right)/two
        point_mean= line(x_mean)
        CALL return_density(point_mean(1), point_mean(2), point_mean(3), &
                            rho_mean)
        IF(rho_mean > min_dens)THEN
          x_left= x_mean
        ELSE
          x_right= x_mean
        ENDIF

      ELSE

        PRINT *
        PRINT *, "** ERROR in SUBROUTINE find_radius in SUBMODULE ", &
                  "bns_base@find_radii!"
        PRINT *, "x_left=", x_left
        PRINT *, "x_right=", x_right
        PRINT *, "point_left=", point_left
        PRINT *, "point_right=", point_right
        PRINT *, "rho_left=", rho_left
        PRINT *, "rho_right=", rho_right
        PRINT *
        STOP

      ENDIF

    ENDDO

    IF(debug) PRINT *, x_left
    IF(debug) PRINT *, x_right
    IF(debug) PRINT *, point_left - center
    IF(debug) PRINT *, point_right - center

    DO i= 0, n_pts, 1

       x= x_left + DBLE(i/n_pts)*(x_right - x_left)
       point= line(x)
       CALL return_density(point(1), point(2), point(3), rho)

       IF( rho <= min_dens )THEN
         radius= NORM2(point - center)
         EXIT
       ENDIF

    ENDDO


    CONTAINS

    
    PURE FUNCTION line(x) RESULT(point)
    !# The line \(c+xv\) parametrized by x, with \(c\) being `center`
    !  and \(v\) being `vector`

      DOUBLE PRECISION, INTENT(IN):: x
      !! Parameter along the line
      DOUBLE PRECISION, DIMENSION(3):: point
      !! Coordinates of the point on the line at the given \(x\)

      point= center + x*versor

    END FUNCTION line


    SUBROUTINE read_density(x, y, z, rho)
    !# Reads the density using the method [[id_base:read_mass_density]].
    !  This is a wrapper needed as the TARGET of a PROCEDURE POINTER
    !  (PROCEDURE POINTERS cannot point to TYPE-BOUND PROCEDURES unfortunately)

      DOUBLE PRECISION, INTENT(IN) :: x
      DOUBLE PRECISION, INTENT(IN) :: y
      DOUBLE PRECISION, INTENT(IN) :: z
      DOUBLE PRECISION, INTENT(OUT):: rho

      rho= this% read_mass_density(x, y, z)

    END SUBROUTINE read_density


    SUBROUTINE read_density_opt(x, y, z, rho)
    !# Reads the density using the method passed optionally as an argument.

      DOUBLE PRECISION, INTENT(IN) :: x
      DOUBLE PRECISION, INTENT(IN) :: y
      DOUBLE PRECISION, INTENT(IN) :: z
      DOUBLE PRECISION, INTENT(OUT):: rho

      rho= get_density(x, y, z)

    END SUBROUTINE read_density_opt


  END PROCEDURE find_radius


END SUBMODULE find_radii

