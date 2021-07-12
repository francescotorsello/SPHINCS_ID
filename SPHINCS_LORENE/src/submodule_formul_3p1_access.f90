! File:         submodule_formul_3p1_access.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (formul_3p1_id) formul_3p1_access

  !****************************************************
  !                                                   *
  ! Implementation of the methods of TYPE formul_3p1  *
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

    IF( i > THIS% levels(l)% ngrid_x )THEN
      PRINT *, "** ERROR in get_grid_point: i=", i, "> ngrid_x=", &
               THIS% levels(l)% ngrid_x, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( j > THIS% levels(l)% ngrid_y )THEN
      PRINT *, "** ERROR in get_grid_point j=", j, "> ngrid_y=", &
               THIS% levels(l)% ngrid_y, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( k > THIS% levels(l)% ngrid_z )THEN
      PRINT *, "** ERROR in get_grid_point k=", k, "> ngrid_z=", &
               THIS% levels(l)% ngrid_z, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF

    grid_point(1)= THIS% coords% levels(l)% var( i, j, k, jx )
    grid_point(2)= THIS% coords% levels(l)% var( i, j, k, jy )
    grid_point(3)= THIS% coords% levels(l)% var( i, j, k, jz )

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

    nlevels= THIS% nlevels

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

    levels= THIS% levels

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

    dx= THIS% levels(l)% dx

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

    dy= THIS% levels(l)% dy

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

    dz= THIS% levels(l)% dz

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

    ngrid_x= THIS% levels(l)% ngrid_x

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

    ngrid_y= THIS% levels(l)% ngrid_y

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

    ngrid_z= THIS% levels(l)% ngrid_z

  END PROCEDURE get_ngrid_z


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

    IF( i > THIS% levels(l)% ngrid_x )THEN
      PRINT *, "** ERROR in get_HC: i=", i, "> ngrid_x=", &
               THIS% levels(l)% ngrid_x, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( j > THIS% levels(l)% ngrid_y )THEN
      PRINT *, "** ERROR in get_HC: j=", j, "> ngrid_y=", &
               THIS% levels(l)% ngrid_y, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( k > THIS% levels(l)% ngrid_z )THEN
      PRINT *, "** ERROR in get_HC: k=", k, "> ngrid_z=", &
               THIS% levels(l)% ngrid_z, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF

    HC_value= THIS% HC% levels(l)% var( i, j, k )

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

    IF( i > THIS% levels(l)% ngrid_x )THEN
      PRINT *, "** ERROR in get_HC: i=", i, "> ngrid_x=", &
               THIS% levels(l)% ngrid_x, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( j > THIS% levels(l)% ngrid_y )THEN
      PRINT *, "** ERROR in get_HC: j=", j, "> ngrid_y=", &
               THIS% levels(l)% ngrid_y, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( k > THIS% levels(l)% ngrid_z )THEN
      PRINT *, "** ERROR in get_HC: k=", k, "> ngrid_z=", &
               THIS% levels(l)% ngrid_z, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF

    MC_value(1)= THIS% MC% levels(l)% var( i, j, k, jx )
    MC_value(2)= THIS% MC% levels(l)% var( i, j, k, jy )
    MC_value(3)= THIS% MC% levels(l)% var( i, j, k, jz )

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

    IF( i > THIS% levels(l)% ngrid_x )THEN
      PRINT *, "** ERROR in get_HC: i=", i, "> ngrid_x=", &
               THIS% levels(l)% ngrid_x, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( j > THIS% levels(l)% ngrid_y )THEN
      PRINT *, "** ERROR in get_HC: j=", j, "> ngrid_y=", &
               THIS% levels(l)% ngrid_y, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( k > THIS% levels(l)% ngrid_z )THEN
      PRINT *, "** ERROR in get_HC: k=", k, "> ngrid_z=", &
               THIS% levels(l)% ngrid_z, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF

    HC_value= THIS% HC_parts% levels(l)% var( i, j, k )

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

    IF( i > THIS% levels(l)% ngrid_x )THEN
      PRINT *, "** ERROR in get_HC: i=", i, "> ngrid_x=", &
               THIS% levels(l)% ngrid_x, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( j > THIS% levels(l)% ngrid_y )THEN
      PRINT *, "** ERROR in get_HC: j=", j, "> ngrid_y=", &
               THIS% levels(l)% ngrid_y, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF
    IF( k > THIS% levels(l)% ngrid_z )THEN
      PRINT *, "** ERROR in get_HC: k=", k, "> ngrid_z=", &
               THIS% levels(l)% ngrid_z, "on refinement level l=", l
      PRINT *
      STOP
    ENDIF

    MC_value(1)= THIS% MC_parts% levels(l)% var( i, j, k, jx )
    MC_value(2)= THIS% MC_parts% levels(l)% var( i, j, k, jy )
    MC_value(3)= THIS% MC_parts% levels(l)% var( i, j, k, jz )

  END PROCEDURE get_MC_parts


END SUBMODULE formul_3p1_access
