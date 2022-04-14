! File:         submodule_standard_tpo_formulation_access.f90
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

SUBMODULE (standard_tpo_formulation) access

  !****************************************************
  !                                                   *
  ! Implementation of the methods of TYPE tpo_formulation  *
  ! that allow to access PRIVATE members.             *
  !                                                   *
  ! FT 12.07.2021                                     *
  !                                                   *
  !****************************************************


  IMPLICIT NONE


  CONTAINS


  !-----------------!
  !--  FUNCTIONS  --!
  !-----------------!


  MODULE PROCEDURE get_grid_point

    !*************************************************
    !                                                *
    ! Returns the array with the coordinates of the  *
    ! grid point                                     *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    USE tensor, ONLY: jx, jy, jz

    IMPLICIT NONE

    IF( i > this% levels(l)% ngrid_x )THEN
      PRINT *, "** ERROR in get_grid_point: i=", i, "> ngrid_x=", &
               this% levels(l)% ngrid_x, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( j > this% levels(l)% ngrid_y )THEN
      PRINT *, "** ERROR in get_grid_point j=", j, "> ngrid_y=", &
               this% levels(l)% ngrid_y, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( k > this% levels(l)% ngrid_z )THEN
      PRINT *, "** ERROR in get_grid_point k=", k, "> ngrid_z=", &
               this% levels(l)% ngrid_z, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF

    grid_point(1)= this% coords% levels(l)% var( i, j, k, jx )
    grid_point(2)= this% coords% levels(l)% var( i, j, k, jy )
    grid_point(3)= this% coords% levels(l)% var( i, j, k, jz )

  END PROCEDURE get_grid_point


  MODULE PROCEDURE get_nlevels

    !**************************************************
    !                                                 *
    ! Returns the number of refinement levels nlevels *
    !                                                 *
    ! FT 26.03.2021                                   *
    !                                                 *
    !**************************************************

    IMPLICIT NONE

    nlevels= this% nlevels

  END PROCEDURE get_nlevels


  MODULE PROCEDURE get_levels

    !**************************************************
    !                                                 *
    ! Returns the data structure levels               *
    !                                                 *
    ! FT 26.03.2021                                   *
    !                                                 *
    !**************************************************

    IMPLICIT NONE

    levels= this% levels

  END PROCEDURE get_levels


  MODULE PROCEDURE get_dx

    !*************************************************
    !                                                *
    ! Returns the grid spacing on the x axis         *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    dx= this% levels(l)% dx

  END PROCEDURE get_dx


  MODULE PROCEDURE get_dy

    !*************************************************
    !                                                *
    ! Returns the grid spacing on the y axis         *
    !                                                *
    ! FT 26.03.2021                                  *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    dy= this% levels(l)% dy

  END PROCEDURE get_dy


  MODULE PROCEDURE get_dz

    !*************************************************
    !                                                *
    ! Returns the grid spacing on the z axis         *
    !                                                *
    ! FT 26.03.2021                                  *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    dz= this% levels(l)% dz

  END PROCEDURE get_dz


  MODULE PROCEDURE get_ngrid_x

    !*************************************************
    !                                                *
    ! Returns the number of grid points on the x     *
    ! axis                                           *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    ngrid_x= this% levels(l)% ngrid_x

  END PROCEDURE get_ngrid_x


  MODULE PROCEDURE get_ngrid_y

    !*************************************************
    !                                                *
    ! Returns the number of grid points on the y     *
    ! axis                                           *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    ngrid_y= this% levels(l)% ngrid_y

  END PROCEDURE get_ngrid_y


  MODULE PROCEDURE get_ngrid_z

    !*************************************************
    !                                                *
    ! Returns the number of grid points on the z     *
    ! axis                                           *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    ngrid_z= this% levels(l)% ngrid_z

  END PROCEDURE get_ngrid_z


  MODULE PROCEDURE get_xR

    !***********************************************
    !                                              *
    ! Returns the x boundary of refinement level l *
    !                                              *
    ! FT                                           *
    !                                              *
    !***********************************************

    IMPLICIT NONE

    xR= this% levels(l)% xR

  END PROCEDURE get_xR


  MODULE PROCEDURE get_yR

    !***********************************************
    !                                              *
    ! Returns the y boundary of refinement level l *
    !                                              *
    ! FT                                           *
    !                                              *
    !***********************************************

    IMPLICIT NONE

    yR= this% levels(l)% yR

  END PROCEDURE get_yR


  MODULE PROCEDURE get_zR

    !***********************************************
    !                                              *
    ! Returns the z boundary of refinement level l *
    !                                              *
    ! FT                                           *
    !                                              *
    !***********************************************

    IMPLICIT NONE

    zR= this% levels(l)% zR

  END PROCEDURE get_zR


  MODULE PROCEDURE get_HC

    !**************************************************
    !                                                 *
    ! Returns the value of the Hamiltonian constraint *
    ! at the specified grid point                     *
    !                                                 *
    ! FT                                              *
    !                                                 *
    !**************************************************

    IMPLICIT NONE

    IF( i > this% levels(l)% ngrid_x )THEN
      PRINT *, "** ERROR in get_HC: i=", i, "> ngrid_x=", &
               this% levels(l)% ngrid_x, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( j > this% levels(l)% ngrid_y )THEN
      PRINT *, "** ERROR in get_HC: j=", j, "> ngrid_y=", &
               this% levels(l)% ngrid_y, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( k > this% levels(l)% ngrid_z )THEN
      PRINT *, "** ERROR in get_HC: k=", k, "> ngrid_z=", &
               this% levels(l)% ngrid_z, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF

    HC_value= this% HC% levels(l)% var( i, j, k )

  END PROCEDURE get_HC


  MODULE PROCEDURE get_MC

    !**************************************************
    !                                                 *
    ! Returns the array of values of the momentum     *
    ! constraint at the specified grid point          *
    !                                                 *
    ! FT                                              *
    !                                                 *
    !**************************************************

    USE tensor, ONLY: jx, jy, jz

    IMPLICIT NONE

    IF( i > this% levels(l)% ngrid_x )THEN
      PRINT *, "** ERROR in get_HC: i=", i, "> ngrid_x=", &
               this% levels(l)% ngrid_x, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( j > this% levels(l)% ngrid_y )THEN
      PRINT *, "** ERROR in get_HC: j=", j, "> ngrid_y=", &
               this% levels(l)% ngrid_y, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( k > this% levels(l)% ngrid_z )THEN
      PRINT *, "** ERROR in get_HC: k=", k, "> ngrid_z=", &
               this% levels(l)% ngrid_z, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF

    MC_value(1)= this% MC% levels(l)% var( i, j, k, jx )
    MC_value(2)= this% MC% levels(l)% var( i, j, k, jy )
    MC_value(3)= this% MC% levels(l)% var( i, j, k, jz )

  END PROCEDURE get_MC


  MODULE PROCEDURE get_HC_parts

    !**************************************************
    !                                                 *
    ! Returns the value of the Hamiltonian constraint *
    ! computed with particle data, at the specified   *
    ! grid point                                      *
    !                                                 *
    ! FT                                              *
    !                                                 *
    !**************************************************

    IMPLICIT NONE

    IF( i > this% levels(l)% ngrid_x )THEN
      PRINT *, "** ERROR in get_HC: i=", i, "> ngrid_x=", &
               this% levels(l)% ngrid_x, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( j > this% levels(l)% ngrid_y )THEN
      PRINT *, "** ERROR in get_HC: j=", j, "> ngrid_y=", &
               this% levels(l)% ngrid_y, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( k > this% levels(l)% ngrid_z )THEN
      PRINT *, "** ERROR in get_HC: k=", k, "> ngrid_z=", &
               this% levels(l)% ngrid_z, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF

    HC_value= this% HC_parts% levels(l)% var( i, j, k )

  END PROCEDURE get_HC_parts


  MODULE PROCEDURE get_MC_parts

    !**************************************************
    !                                                 *
    ! Returns the value of the momentum constraint    *
    ! computed with particle data, at the specified   *
    ! grid point                                      *
    !                                                 *
    ! FT                                              *
    !                                                 *
    !**************************************************

    USE tensor, ONLY: jx, jy, jz

    IMPLICIT NONE

    IF( i > this% levels(l)% ngrid_x )THEN
      PRINT *, "** ERROR in get_HC: i=", i, "> ngrid_x=", &
               this% levels(l)% ngrid_x, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( j > this% levels(l)% ngrid_y )THEN
      PRINT *, "** ERROR in get_HC: j=", j, "> ngrid_y=", &
               this% levels(l)% ngrid_y, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( k > this% levels(l)% ngrid_z )THEN
      PRINT *, "** ERROR in get_HC: k=", k, "> ngrid_z=", &
               this% levels(l)% ngrid_z, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF

    MC_value(1)= this% MC_parts% levels(l)% var( i, j, k, jx )
    MC_value(2)= this% MC_parts% levels(l)% var( i, j, k, jy )
    MC_value(3)= this% MC_parts% levels(l)% var( i, j, k, jz )

  END PROCEDURE get_MC_parts


END SUBMODULE access
