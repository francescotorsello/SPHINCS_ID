! File:         module_utility.f90
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

MODULE utility

  !***********************************************************************
  !
  !# This MODULE contains useful miscellaneous PROCEDURES and variables
  !
  !***********************************************************************


  USE matrix, ONLY: determinant_4x4_matrix


  IMPLICIT NONE


  INTEGER:: itr
  !! Iterator for loops
  INTEGER:: itr3
  !! Iterator for loops
  INTEGER:: itr4
  !! Iterator for loops
  INTEGER:: ios
  !! Variable to store the state of I/O
  INTEGER:: cnt= 0
  !! Counter

  !
  !-- Variables to print progress on screen
  !

  INTEGER:: perc
  !! Percentage of execution time (integer)
  DOUBLE PRECISION:: perc2
  !! Percentage of execution time (double)
  CHARACTER, PARAMETER:: creturn= ACHAR(13)
  !! Carriage return

  LOGICAL:: file_exists
  !! `.TRUE.` if a given file exists; `.FALSE.` otherwise
  LOGICAL:: show_progress
  !# `.TRUE.` if loop progress is to be printed to standard output;
  !  `.FALSE.` otherwise

  CHARACTER( LEN= : ), ALLOCATABLE:: err_msg
  !! String storing error messages

  !
  !-- Variables used to set the run_id
  !

  CHARACTER(8)  :: date
  !! Date when the run starts
  CHARACTER(10) :: time
  !! Time when the run starts
  CHARACTER(5)  :: zone
  !! Place where the run runs
  INTEGER, DIMENSION(8) :: values
  !# An integer array of 8 elements described below:
  !
  !   1. The year as a 4-digit integer
  !   2. The month as an integer from 1 to 12
  !   3. The day of the moneth as an integer from 1 to 31
  !   4. The time difference, in minutes, with respect to UTC
  !      (Coordinated Universal Time)
  !   5. The hour of the day as an integer from 1 to 23
  !   6. The minutes of the hour as an integer from 1 to 59
  !   7. The second of the minute as an integer from 0 to 60
  !   8. The millisecond of the second as an integer from 0 to 999
  CHARACTER( LEN= 19 ):: run_id
  !! Identification string for the run
  CHARACTER( LEN= 19 ):: end_time
  !! Time when the run ends


  CONTAINS


  SUBROUTINE test_status( io_stat, io_msg, opt_msg )

    !***********************************************
    !
    !# Test if a status variable is 0 or not
    !
    !  FT 17.09.2020
    !
    !***********************************************

    IMPLICIT NONE

    INTEGER,               INTENT(IN)           :: io_stat
    !! Status variable
    CHARACTER( LEN= 100 ), INTENT(IN)           :: io_msg
    !! Status message
    CHARACTER( LEN= * ),   INTENT(IN), OPTIONAL :: opt_msg
    !! Optional status message

    IF( io_stat > 0 )THEN

      PRINT *
      PRINT *, "***** ERROR! IOSTAT > 0. ", &
               "The error message is: ", io_msg
      IF( PRESENT( opt_msg ) )THEN
        PRINT *, opt_msg
      ENDIF
      PRINT *
      STOP

    ENDIF

  END SUBROUTINE test_status


  PURE FUNCTION is_finite_number( x ) RESULT( res )

    !***********************************************
    !
    !# Test if a double precision is a finite number
    !
    !  FT 11.02.2022
    !
    !***********************************************

    USE, INTRINSIC:: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE

    DOUBLE PRECISION, INTENT(IN):: x
    LOGICAL:: res

    res= (.NOT.ISNAN(x)) .AND. IEEE_IS_FINITE(x)

  END FUNCTION is_finite_number


  PURE SUBROUTINE compute_g4( lapse, shift, g3, g4 )

    !***********************************************
    !
    !# Computes the spacetime metric from lapse,
    !  shift and spatial metric
    !
    !  FT 27.11.2020
    !
    !  Generalized to not be bound to the mesh
    !
    !  FT 07.02.2022
    !
    !***********************************************

    USE constants,  ONLY: two
    USE tensor,     ONLY: itt, itx, ity, itz, ixx, ixy, &
                          ixz, iyy, iyz, izz, jxx, jxy, jxz, &
                          jyy, jyz, jzz, jx, jy, jz

    IMPLICIT NONE

    DOUBLE PRECISION,                INTENT(IN)   :: lapse
    !! Lapse function
    DOUBLE PRECISION, DIMENSION(3),  INTENT(IN)   :: shift
    !! Contravariant shift vector
    DOUBLE PRECISION, DIMENSION(6),  INTENT(IN)   :: g3
    !! Covariant spatial metric
    DOUBLE PRECISION, DIMENSION(10), INTENT(INOUT):: g4
    !! Covariant spacetime metric

    g4(itt)= - lapse*lapse + g3(jxx)*shift(jx)*shift(jx)     &
                           + g3(jxy)*shift(jx)*shift(jy)*two &
                           + g3(jxz)*shift(jx)*shift(jz)*two &
                           + g3(jyy)*shift(jy)*shift(jy)     &
                           + g3(jyz)*shift(jy)*shift(jz)*two &
                           + g3(jzz)*shift(jz)*shift(jz)

    g4(itx)= g3(jxx)*shift(jx) + g3(jxy)*shift(jy) + g3(jxz)*shift(jz)
    g4(ity)= g3(jxy)*shift(jx) + g3(jyy)*shift(jy) + g3(jyz)*shift(jz)
    g4(itz)= g3(jxz)*shift(jx) + g3(jyz)*shift(jy) + g3(jzz)*shift(jz)

    g4(ixx)= g3(jxx)
    g4(ixy)= g3(jxy)
    g4(ixz)= g3(jxz)
    g4(iyy)= g3(jyy)
    g4(iyz)= g3(jyz)
    g4(izz)= g3(jzz)

  END SUBROUTINE compute_g4


  SUBROUTINE determinant_sym4x4( A, det )

    !****************************************************************
    !
    !# Compute the determinant of a \(4\times 4\) symmetric matrix
    !  field, given as a 10-vector
    !
    !  FT 26.01.2022
    !
    !  Generalized to not be bound to the mesh
    !
    !  FT 07.02.2022
    !
    !****************************************************************

    USE tensor, ONLY: itt, itx, ity, itz, ixx, ixy, ixz, iyy, iyz, izz, n_sym4x4

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: A(:)
    !# The \(4\times 4\) symmetric matrix, given as a 10-vector.
    !  The first 3 components run over the numbers of grid points
    !  along each axis. The fourth index runs over the number of
    !  independent components of the \(4\times 4\) symmetric matrix.
    DOUBLE PRECISION, INTENT(OUT):: det
    !! Determinant of the \(4\times 4\) symmetric matrix

    IF( SIZE(A) /= n_sym4x4 )THEN
      PRINT *, "** ERROR in determinant_sym4x4_grid in MODULE utility.", &
               " This subroutine needs a symmetric matrix with 10 components,",&
               " and a ", SIZE(A), "component matrix was given instead."
      STOP
    ENDIF

    det=   A(itt)*(A(ixx)*(A(iyy)*A(izz) &
         - A(iyz)* A(iyz)) &
         + A(ixy)*(A(iyz)* A(ixz) &
         - A(ixy)* A(izz)) &
         + A(ixz)*(A(ixy)* A(iyz) &
         - A(iyy)* A(ixz))) &
         - A(itx)*(A(itx)*(A(iyy)*A(izz) &
         - A(iyz)* A(iyz)) &
         + A(ixy)*(A(iyz)* A(itz) &
         - A(ity)* A(izz)) &
         + A(ixz)*(A(ity)* A(iyz) &
         - A(iyy)* A(itz))) &
         + A(ity)*(A(itx)*(A(ixy)*A(izz) &
         - A(iyz)* A(ixz)) &
         + A(ixx)*(A(iyz)* A(itz) &
         - A(ity)* A(izz)) &
         + A(ixz)*(A(ity)* A(ixz) &
         - A(ixy)* A(itz))) &
         - A(itz)*(A(itx)*(A(ixy)*A(iyz) &
         - A(iyy)* A(ixz)) &
         + A(ixx)*(A(iyy)* A(itz) &
         - A(ity)* A(iyz)) &
         + A(ixy)*(A(ity)* A(ixz) &
         - A(ixy)* A(itz)))

  END SUBROUTINE determinant_sym4x4


  SUBROUTINE spacetime_vector_norm_sym4x4( g4, v, norm )

    !****************************************************************
    !
    !# Compute the spacetime norm of a vector, using the metric
    !  given as an array of 10 components
    !
    !  FT 07.02.2022
    !
    !****************************************************************

    USE tensor,    ONLY: itt, itx, ity, itz, ixx, ixy, ixz, iyy, iyz, izz, &
                         it, ix, iy, iz, n_sym4x4
    USE constants, ONLY: two

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: g4(itt:izz)
    !# The \(4\times 4\) spacetime metric, given as a 10-vector.
    DOUBLE PRECISION, INTENT(IN):: v(it:iz)
    !# The \(4\)-vector whose norm has to be computed.
    DOUBLE PRECISION, INTENT(OUT):: norm
    !! Spacetime norm of the vector v.

    IF( SIZE(g4) /= n_sym4x4 )THEN
      PRINT *, "** ERROR in determinant_sym4x4_grid in MODULE utility.", &
               " This subroutine needs a symmetric matrix with 10 components,",&
               " and a ", SIZE(g4), "component matrix was given instead."
      STOP
    ENDIF

    norm= g4(itt)*v(it)*v(it)     + two*g4(itx)*v(it)*v(ix) &
        + two*g4(ity)*v(it)*v(iy) + two*g4(itz)*v(it)*v(iz) &
        + g4(ixx)*v(ix)*v(ix)     + two*g4(ixy)*v(ix)*v(iy) &
        + two*g4(ixz)*v(ix)*v(iz) + g4(iyy)*v(iy)*v(iy)     &
        + two*g4(iyz)*v(iy)*v(iz) + g4(izz)*v(iz)*v(iz)

  END SUBROUTINE spacetime_vector_norm_sym4x4


  SUBROUTINE determinant_sym3x3( A, det )

    !****************************************************************
    !
    !# Compute the determinant of a \(3\times 3\) symmetric matrix
    !  field, given as a 6-vector, at a given grid point
    !
    !  FT 26.03.2021
    !
    !  Generalized to not be bound to the mesh
    !
    !  FT 07.02.2022
    !
    !****************************************************************

    USE tensor, ONLY: jxx, jxy, jxz, jyy, jyz, jzz, n_sym3x3

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: A(jxx:jzz)
    !# The \(3\times 3\) symmetric matrix, given as a 6-vector.
    DOUBLE PRECISION, INTENT(OUT):: det
    !! Determinant of the \(3\times 3\) symmetric matrix

    IF( SIZE(A) /= n_sym3x3 )THEN
      PRINT *, "** ERROR in determinant_sym3x3_grid in MODULE utility.", &
               " This subroutine needs an array with 6 components,",&
               " and a ", SIZE(A), "component array was given instead."
      STOP
    ENDIF

    det=   A(jxx)*A(jyy)*A(jzz) &
         + A(jxy)*A(jyz)*A(jxz) &
         + A(jxz)*A(jxy)*A(jyz) &
         - A(jxy)*A(jxy)*A(jzz) &
         - A(jxz)*A(jyy)*A(jxz) &
         - A(jxx)*A(jyz)*A(jyz)

  END SUBROUTINE determinant_sym3x3


  PURE SUBROUTINE spherical_from_cartesian( x, y, z, xo, yo, zo, r, theta, phi )

    !****************************************************************
    !
    !# Compute the spherical polar coordinates of a point \(p\)
    !  relative to a point \(O\), starting from the Cartesian
    !  coordinates of the points \(p\) and \(O\)
    !
    !  FT 14.01.2021
    !
    !****************************************************************

    USE constants, ONLY: zero, pi

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: x
    !! \(x\) coordinate of the point \(p\)
    DOUBLE PRECISION, INTENT(IN):: y
    !! \(y\) coordinate of the point \(p\)
    DOUBLE PRECISION, INTENT(IN):: z
    !! \(z\) coordinate of the point \(p\)

    DOUBLE PRECISION, INTENT(IN):: xo
    !! \(x\) coordinate of the point \(O\)
    DOUBLE PRECISION, INTENT(IN):: yo
    !! \(y\) coordinate of the point \(O\)
    DOUBLE PRECISION, INTENT(IN):: zo
    !! \(z\) coordinate of the point \(O\)

    DOUBLE PRECISION, INTENT(OUT):: r
    !! \(r\) coordinate of the point \(p\), relative to \(O\)
    DOUBLE PRECISION, INTENT(OUT):: theta
    !! \(\theta\) coordinate (colatitude) of the point \(p\), relative to \(O\)
    DOUBLE PRECISION, INTENT(OUT):: phi
    !! \(\phi\) coordinate (azimuth) of the point \(p\), relative to \(O\)

    DOUBLE PRECISION:: xd, yd, zd

    xd= x - xo
    yd= y - yo
    zd= z - zo

    IF( x > zero )THEN

      phi= ATAN( yd/xd )

    ELSEIF( x < zero )THEN

      phi= ATAN( yd/xd ) + pi

    ELSE

      phi= pi/2.D0

    ENDIF

    r= SQRT( xd**2.D0 + yd**2.D0 + zd**2.D0 )

    theta= ACOS( zd/r )

  END SUBROUTINE spherical_from_cartesian


END MODULE utility
