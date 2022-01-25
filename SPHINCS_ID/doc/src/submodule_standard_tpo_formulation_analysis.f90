! File:         submodule_standard_tpo_formulation_analysis.f90
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

SUBMODULE (standard_tpo_formulation) analysis

  !****************************************************
  !
  !# Implementation of the methods of TYPE
  !  tpo_formulation that analyze a grid function.
  !
  !  FT 12.07.2021
  !
  !****************************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE analyze_constraint

    !****************************************************
    !
    !#
    !
    !  FT
    !
    !****************************************************

    IMPLICIT NONE

    INTEGER:: cnt_m7, cnt_m6, cnt_m5, cnt_m4, cnt_m3, cnt_m2, cnt_m1, cnt_0, &
              cnt_p1, cnt_p2, cnt_p3, cnt_oo, grid_points, i, j, k, &
              unit_analysis, nx, ny, nz

    DOUBLE PRECISION:: tmp, total

    LOGICAL:: exist
    !LOGICAL, PARAMETER:: DEBUG= .FALSE.

    IF( THIS% export_constraints_details )THEN
      !
      !-- Export the constraint analysis to a formatted file
      !
      unit_analysis= 20120

      INQUIRE( FILE= TRIM(name_analysis), EXIST= exist )

      IF( exist )THEN
          OPEN( UNIT= unit_analysis, FILE= TRIM(name_analysis), &
                STATUS= "REPLACE", &
                FORM= "FORMATTED", &
                POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
                IOMSG= err_msg )
      ELSE
          OPEN( UNIT= unit_analysis, FILE= TRIM(name_analysis), STATUS= "NEW", &
                FORM= "FORMATTED", &
                ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
      ENDIF
      IF( ios > 0 )THEN
        PRINT *, "... error when opening " // TRIM(name_analysis), &
                 ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when opening " &
      !         // TRIM(name_analysis) )

      WRITE( UNIT = unit_analysis, IOSTAT = ios, IOMSG = err_msg, &
      FMT = * ) &
      "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id
      WRITE( UNIT= unit_analysis, IOSTAT = ios, &
             IOMSG = err_msg, FMT = * ) &
      "# The rows contain the points (l,x,y,z), with l refinement level, ", &
      "at which ", name_constraint, &
      " has values in: (-oo,1D-7], [1D-7,1D-6], [1D-6,1D-5], [1D-5,1D-4]" &
      // ", [1D-4,1D-3], [1D-3,1D-2], [1D-2,1D-1], [1D-1,1], [1,1D+1]" &
      // ", [1D+1,1D+2], [1D+2,1D+3], [1D+3,+oo)"
    ENDIF

    CALL THIS% abs_values_in( 0.0D0, 1.0D-7, constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_m7 )

    CALL THIS% abs_values_in( 1.0D-7, 1.0D-6, constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_m6 )

    CALL THIS% abs_values_in( 1.0D-6, 1.0D-5, constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_m5 )

    CALL THIS% abs_values_in( 1.0D-5, 1.0D-4, constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_m4 )

    CALL THIS% abs_values_in( 1.0D-4, 1.0D-3, constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_m3 )

    CALL THIS% abs_values_in( 1.0D-3, 1.0D-2, constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_m2 )

    CALL THIS% abs_values_in( 1.0D-2, 1.0D-1, constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_m1 )

    CALL THIS% abs_values_in( 1.0D-1, 1.0D0, constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_0 )

    CALL THIS% abs_values_in( 1.0D0, 1D+1, constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_p1 )

    CALL THIS% abs_values_in( 1.0D+1, 1.0D+2, constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_p2 )

    CALL THIS% abs_values_in( 1.0D+2, 1.0D+3, constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_p3 )

    CALL THIS% abs_values_in( 1.0D+3, HUGE(DBLE(1.0D0)), constraint, l, &
                              THIS% export_constraints_details, &
                              unit_analysis, cnt_oo )

    CLOSE( UNIT= unit_analysis )

    IF( THIS% export_constraints_details )THEN
      PRINT *, " * The details about the absolute values of ", &
               name_constraint, " are printed to ", name_analysis
    ENDIF

    nx= THIS% get_ngrid_x(l)
    ny= THIS% get_ngrid_y(l)
    nz= THIS% get_ngrid_z(l)
    grid_points= nx*ny*nz

    !
    !-- Compute the l2 norm of the constraints
    !
    l2_norm= 0.0D0
    DO k= 1, nz, 1
      DO j= 1, ny, 1
        DO i= 1, nx, 1
          l2_norm= l2_norm + constraint(i,j,k)*constraint(i,j,k)
        ENDDO
      ENDDO
    ENDDO
    l2_norm= SQRT( l2_norm/grid_points )

    !
    !-- Compute the loo norm (supremum norm) of the constraints
    !
    loo_norm= 0.0D0
    DO k= 1, nz, 1
      DO j= 1, ny, 1
        DO i= 1, nx, 1
          tmp= ABS( constraint(i,j,k) )
          IF( tmp > loo_norm )THEN
            loo_norm= tmp
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    !
    !-- Write a summary of the results to the logfile
    !
    WRITE( UNIT= unit_logfile, FMT = * ) "# The absolute values of ", &
         name_constraint, &
         " on the gravity grid are in the following intervals, on the ", &
         "given percentage of grid points:"
    WRITE( UNIT= unit_logfile, FMT = * ) ""
    WRITE( UNIT= unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (-oo,1D-7]: ", &
        100.0D0*DBLE(cnt_m7)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D-7,1D-6]: ", &
        100.0D0*DBLE(cnt_m6)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D-6,1D-5]: ", &
        100.0D0*DBLE(cnt_m5)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D-5,1D-4]: ", &
        100.0D0*DBLE(cnt_m4)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D-4,1D-3]: ", &
        100.0D0*DBLE(cnt_m3)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D-3,1D-2]: ", &
        100.0D0*DBLE(cnt_m2)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D-2,1D-1]: ", &
        100.0D0*DBLE(cnt_m1)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D-1,1]: ",   &
        100.0D0*DBLE(cnt_0)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1,10]: ",      &
        100.0D0*DBLE(cnt_p1)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (10,1D+2]: ",   &
        100.0D0*DBLE(cnt_p2)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D+2,1D+3]: ", &
        100.0D0*DBLE(cnt_p3)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D+3,+oo]: ",  &
        100.0D0*DBLE(cnt_oo)/DBLE(grid_points), "%"
    WRITE( UNIT = unit_logfile, FMT = * )
    WRITE( UNIT = unit_logfile, FMT = * ) "# l2-norm of ", name_constraint,&
                               " over the gravity grid= ", l2_norm
    WRITE( UNIT = unit_logfile, FMT = * )
    WRITE( UNIT = unit_logfile, FMT = * ) &
                       "# loo-norm (supremum of the absolute values) of ", &
                       name_constraint,&
                       " over the gravity grid= ", loo_norm
    WRITE( UNIT = unit_logfile, FMT = * )

    total= 100.0D0*DBLE(cnt_m7)/DBLE(grid_points) + &
           100.0D0*DBLE(cnt_m6)/DBLE(grid_points) + &
           100.0D0*DBLE(cnt_m5)/DBLE(grid_points) + &
           100.0D0*DBLE(cnt_m4)/DBLE(grid_points) + &
           100.0D0*DBLE(cnt_m3)/DBLE(grid_points) + &
           100.0D0*DBLE(cnt_m2)/DBLE(grid_points) + &
           100.0D0*DBLE(cnt_m1)/DBLE(grid_points) + &
           100.0D0*DBLE(cnt_0) /DBLE(grid_points) + &
           100.0D0*DBLE(cnt_p1)/DBLE(grid_points) + &
           100.0D0*DBLE(cnt_p2)/DBLE(grid_points) + &
           100.0D0*DBLE(cnt_p3)/DBLE(grid_points) + &
           100.0D0*DBLE(cnt_oo)/DBLE(grid_points)

    IF( total - DBLE(100) > 1.0D-4 )THEN
      PRINT *, " * WARNING! The percentages of the absolute values of ", &
               name_constraint, &
               " in the given intervals do not sum up to 100%. ", &
               " They sum up to", total, "%"
      PRINT *, "Check what happens in the SUBROUTINES analyze_constraint ", &
               "and abs_values_in."
    ENDIF

  END PROCEDURE analyze_constraint


  MODULE PROCEDURE abs_values_in

    !**************************************************
    !                                                 *
    ! Set "cnt" equal to the number of times that the *
    ! absolute value of "constraint" is in            *
    ! (lower_bound,upper_bound].                      *
    ! Depending on "export", it prints to file the    *
    ! grid points at which this happens.              *
    !                                                 *
    ! FT 24.03.2021                                   *
    !                                                 *
    !**************************************************

    USE tensor, ONLY: jx, jy, jz

    IMPLICIT NONE

    INTEGER:: i, j, k, nx, ny, nz

    nx= THIS% get_ngrid_x(l)
    ny= THIS% get_ngrid_y(l)
    nz= THIS% get_ngrid_z(l)

    cnt= 0
    DO k= 1, nz, 1
      DO j= 1, ny, 1
        DO i= 1, nx, 1
          IF( ABS( constraint(i,j,k) ) >= lower_bound .AND. &
              ABS( constraint(i,j,k) ) < upper_bound )THEN
            cnt= cnt + 1
            IF( export )THEN
              WRITE( UNIT= unit_analysis, IOSTAT = ios, &
                     IOMSG = err_msg, FMT = "(F17.13)", &
                     ADVANCE= "NO" ) &!"(E15.6)" &
                     l
              WRITE( UNIT= unit_analysis, IOSTAT = ios, &
                     IOMSG = err_msg, FMT = "(F17.13)", &
                     ADVANCE= "NO" ) &!"(E15.6)" &
                     THIS% coords% levels(l)% var( i, j, k, jx )
              WRITE( UNIT= unit_analysis, FMT= "(A2)", &
                     ADVANCE= "NO" ) "  "
              WRITE( UNIT= unit_analysis, IOSTAT = ios, &
                     IOMSG = err_msg, FMT = "(F17.13)", &
                     ADVANCE= "NO" ) &
                     THIS% coords% levels(l)% var( i, j, k, jy )
              WRITE( UNIT= unit_analysis, FMT= "(A2)", &
                     ADVANCE= "NO" ) "  "
              WRITE( UNIT= unit_analysis, IOSTAT = ios, &
                     IOMSG = err_msg, FMT = "(F17.13)", &
                     ADVANCE= "NO" ) &
                     THIS% coords% levels(l)% var( i, j, k, jz )
              WRITE( UNIT= unit_analysis, FMT= "(A2)", &
                     ADVANCE= "NO" ) "  "
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    IF( export .AND. cnt == 0 )THEN
      WRITE( UNIT= unit_analysis, FMT= "(I1)", ADVANCE= "NO" ) 0
      WRITE( UNIT= unit_analysis, FMT= "(I1)", ADVANCE= "NO" ) 0
      WRITE( UNIT= unit_analysis, FMT= "(I1)", ADVANCE= "NO" ) 0
    ENDIF
    WRITE( UNIT= unit_analysis, FMT= * ) ""

  END PROCEDURE abs_values_in


END SUBMODULE analysis
