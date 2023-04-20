! File:         submodule_standard_tpo_formulation_analysis.f90
! Authors:      Francesco Torsello (FT)
!************************************************************************
! Copyright (C) 2020-2023 Francesco Torsello                            *
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
    !# Count on how many grid points the argument 'constraint'
    !  has values lying in the intervals \( (-oo,10^{-7}],
    !  [10^{-7},10^{-6}],
    !  [10^{-6},10^{-5}], [10^{-5},10^{-4}], [10^{-4},10^{-3}],
    !  [10^{-3},10^{-2}], [10^{-2},10^{-1}], [10^{-1},1],
    !  [1,10^1], [10^1,10^2], [10^2,10^3],
    !  [10^3,+oo) \)
    !
    !  FT
    !
    !****************************************************

    USE constants, ONLY: pi
    USE utility,   ONLY: determinant_sym3x3, zero, one, two, four, ten

    IMPLICIT NONE

    INTEGER:: cnt_m7, cnt_m6, cnt_m5, cnt_m4, cnt_m3, cnt_m2, cnt_m1, cnt_0, &
              cnt_p1, cnt_p2, cnt_p3, cnt_oo, grid_points, i, j, k, &
              unit_analysis, nx, ny, nz

    DOUBLE PRECISION:: total, dx, dy, dz, detg3

    LOGICAL:: exist
    !LOGICAL, PARAMETER:: DEBUG= .FALSE.

    IF( this% export_constraints_details )THEN
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

    CALL this% abs_values_in( zero, 1.0D-7, constraint, l, &
                              this% export_constraints_details, &
                              unit_analysis, cnt_m7 )

    CALL this% abs_values_in( 1.0D-7, 1.0D-6, constraint, l, &
                              this% export_constraints_details, &
                              unit_analysis, cnt_m6 )

    CALL this% abs_values_in( 1.0D-6, 1.0D-5, constraint, l, &
                              this% export_constraints_details, &
                              unit_analysis, cnt_m5 )

    CALL this% abs_values_in( 1.0D-5, 1.0D-4, constraint, l, &
                              this% export_constraints_details, &
                              unit_analysis, cnt_m4 )

    CALL this% abs_values_in( 1.0D-4, 1.0D-3, constraint, l, &
                              this% export_constraints_details, &
                              unit_analysis, cnt_m3 )

    CALL this% abs_values_in( 1.0D-3, 1.0D-2, constraint, l, &
                              this% export_constraints_details, &
                              unit_analysis, cnt_m2 )

    CALL this% abs_values_in( 1.0D-2, 1.0D-1, constraint, l, &
                              this% export_constraints_details, &
                              unit_analysis, cnt_m1 )

    CALL this% abs_values_in( 1.0D-1, 1.0D0, constraint, l, &
                              this% export_constraints_details, &
                              unit_analysis, cnt_0 )

    CALL this% abs_values_in( 1.0D0, 1D+1, constraint, l, &
                              this% export_constraints_details, &
                              unit_analysis, cnt_p1 )

    CALL this% abs_values_in( 1.0D+1, 1.0D+2, constraint, l, &
                              this% export_constraints_details, &
                              unit_analysis, cnt_p2 )

    CALL this% abs_values_in( 1.0D+2, 1.0D+3, constraint, l, &
                              this% export_constraints_details, &
                              unit_analysis, cnt_p3 )

    CALL this% abs_values_in( 1.0D+3, HUGE(DBLE(one)), constraint, l, &
                              this% export_constraints_details, &
                              unit_analysis, cnt_oo )

    CLOSE( UNIT= unit_analysis )

    IF( this% export_constraints_details )THEN
      PRINT *, " * The details about the absolute values of ", &
               name_constraint, " are printed to ", name_analysis
    ENDIF

    nx= this% get_ngrid_x(l)
    ny= this% get_ngrid_y(l)
    nz= this% get_ngrid_z(l)
    dx= this% get_dx(l)
    dy= this% get_dy(l)
    dz= this% get_dz(l)
    grid_points= nx*ny*nz

    !
    !-- Compute the l2 norm of the constraints
    !
    l2_norm= zero
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( nx, ny, nz, constraint ) &
    !$OMP             PRIVATE( i, j, k ) &
    !$OMP             REDUCTION( +: l2_norm )
    DO k= 1, nz, 1
      DO j= 1, ny, 1
        DO i= 1, nx, 1
          l2_norm= l2_norm + constraint(i,j,k)*constraint(i,j,k)
        ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    l2_norm= SQRT( l2_norm/grid_points )

    !
    !-- Compute a rough estimate of the integral of the constraints
    !
    integral= zero
    IF( PRESENT(source) )THEN

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( nx, ny, nz, dx, dy, dz, this, &
      !$OMP                     constraint, source, l ) &
      !$OMP             PRIVATE( i, j, k, detg3 ) &
      !$OMP             REDUCTION( +: integral )
      DO k= 1, nz, 1
        DO j= 1, ny, 1
          DO i= 1, nx, 1

            CALL determinant_sym3x3( &
                              this% g_phys3_ll% levels(l)% var(i,j,k,:), detg3 )

            integral= integral &
                    + dx*dy*dz*SQRT(detg3)*( constraint(i,j,k) - source(i,j,k) )

          ENDDO
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

    ELSE

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( nx, ny, nz, dx, dy, dz, this, &
      !$OMP                     constraint, l ) &
      !$OMP             PRIVATE( i, j, k, detg3 ) &
      !$OMP             REDUCTION( +: integral )
      DO k= 1, nz, 1
        DO j= 1, ny, 1
          DO i= 1, nx, 1

            CALL determinant_sym3x3( &
                              this% g_phys3_ll% levels(l)% var(i,j,k,:), detg3 )

            integral= integral &
                      + dx*dy*dz*SQRT(detg3)*constraint(i,j,k)

          ENDDO
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

    ENDIF
    integral= integral/(two*four*pi)

    !
    !-- Compute the loo norm (supremum norm) of the constraints
    !
    loo_norm= zero
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( nx, ny, nz, constraint ) &
    !$OMP             PRIVATE( i, j, k ) &
    !$OMP             REDUCTION( MAX: loo_norm )
    DO k= 1, nz, 1
      DO j= 1, ny, 1
        DO i= 1, nx, 1
          loo_norm= MAX( loo_norm, ABS( constraint(i,j,k) ) )
        ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    !
    !-- Write a summary of the results to the logfile
    !
    WRITE( UNIT= unit_logfile, FMT = * ) "# The absolute values of ", &
         name_constraint, &
         " on the gravity grid are in the following intervals, on the ", &
         "given percentage of grid points:"
    WRITE( UNIT= unit_logfile, FMT = * ) ""
    WRITE( UNIT= unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (-oo,1D-7]: ", &
        ten*ten*DBLE(cnt_m7)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D-7,1D-6]: ", &
        ten*ten*DBLE(cnt_m6)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D-6,1D-5]: ", &
        ten*ten*DBLE(cnt_m5)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D-5,1D-4]: ", &
        ten*ten*DBLE(cnt_m4)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D-4,1D-3]: ", &
        ten*ten*DBLE(cnt_m3)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D-3,1D-2]: ", &
        ten*ten*DBLE(cnt_m2)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D-2,1D-1]: ", &
        ten*ten*DBLE(cnt_m1)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D-1,1]: ",   &
        ten*ten*DBLE(cnt_0)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1,10]: ",      &
        ten*ten*DBLE(cnt_p1)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (10,1D+2]: ",   &
        ten*ten*DBLE(cnt_p2)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D+2,1D+3]: ", &
        ten*ten*DBLE(cnt_p3)/DBLE(grid_points) , "%"
    WRITE( UNIT = unit_logfile, FMT = "(A15,F5.1,A1)" ) "  (1D+3,+oo]: ",  &
        ten*ten*DBLE(cnt_oo)/DBLE(grid_points), "%"
    WRITE( UNIT = unit_logfile, FMT = * )
    WRITE( UNIT = unit_logfile, FMT = * ) "# l2-norm of ", name_constraint,&
                               " over the gravity grid= ", l2_norm
    WRITE( UNIT = unit_logfile, FMT = * )
    WRITE( UNIT = unit_logfile, FMT = * ) "# Integral of ", name_constraint,&
                               " over the gravity grid= ", integral
    WRITE( UNIT = unit_logfile, FMT = * )
    WRITE( UNIT = unit_logfile, FMT = * ) &
                       "# loo-norm (supremum of the absolute values) of ", &
                       name_constraint,&
                       " over the gravity grid= ", loo_norm
    WRITE( UNIT = unit_logfile, FMT = * )

    total= ten*ten*DBLE(cnt_m7)/DBLE(grid_points) + &
           ten*ten*DBLE(cnt_m6)/DBLE(grid_points) + &
           ten*ten*DBLE(cnt_m5)/DBLE(grid_points) + &
           ten*ten*DBLE(cnt_m4)/DBLE(grid_points) + &
           ten*ten*DBLE(cnt_m3)/DBLE(grid_points) + &
           ten*ten*DBLE(cnt_m2)/DBLE(grid_points) + &
           ten*ten*DBLE(cnt_m1)/DBLE(grid_points) + &
           ten*ten*DBLE(cnt_0) /DBLE(grid_points) + &
           ten*ten*DBLE(cnt_p1)/DBLE(grid_points) + &
           ten*ten*DBLE(cnt_p2)/DBLE(grid_points) + &
           ten*ten*DBLE(cnt_p3)/DBLE(grid_points) + &
           ten*ten*DBLE(cnt_oo)/DBLE(grid_points)

    IF( total - ten*ten > 1.0D-4 )THEN
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

    nx= this% get_ngrid_x(l)
    ny= this% get_ngrid_y(l)
    nz= this% get_ngrid_z(l)

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
                     this% coords% levels(l)% var( i, j, k, jx )
              WRITE( UNIT= unit_analysis, FMT= "(A2)", &
                     ADVANCE= "NO" ) "  "
              WRITE( UNIT= unit_analysis, IOSTAT = ios, &
                     IOMSG = err_msg, FMT = "(F17.13)", &
                     ADVANCE= "NO" ) &
                     this% coords% levels(l)% var( i, j, k, jy )
              WRITE( UNIT= unit_analysis, FMT= "(A2)", &
                     ADVANCE= "NO" ) "  "
              WRITE( UNIT= unit_analysis, IOSTAT = ios, &
                     IOMSG = err_msg, FMT = "(F17.13)", &
                     ADVANCE= "NO" ) &
                     this% coords% levels(l)% var( i, j, k, jz )
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
    !WRITE( UNIT= unit_analysis, FMT= * ) ""

  END PROCEDURE abs_values_in


END SUBMODULE analysis
