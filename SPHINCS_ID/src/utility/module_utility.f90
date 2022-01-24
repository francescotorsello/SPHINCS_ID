! File:         module_utility.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

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


  PURE SUBROUTINE compute_g4( i, j, k, lapse, shift_u, g_phys3_ll, g4 )

    !***********************************************
    !
    !# Computes the spacetime metric from lapse,
    !  shift and spatial metric
    !
    !  FT 27.11.2020
    !
    !***********************************************

    USE constants,  ONLY: two
    USE tensor,     ONLY: itt, itx, ity, itz, ixx, ixy, &
                          ixz, iyy, iyz, izz, jxx, jxy, jxz, &
                          jyy, jyz, jzz, jx, jy, jz

    IMPLICIT NONE

    INTEGER, INTENT( IN ):: i, j, k
    !! Indices of the grid point
    DOUBLE PRECISION, DIMENSION(:,:,:),   INTENT( IN ):: lapse
    !! Lapse function
    DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN ):: shift_u
    !! Shift vector
    DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN ):: g_phys3_ll
    !! Spatial metric
    DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: g4
    !! Spacetime metric

    g4(i,j,k,itt)= - lapse(i,j,k)*lapse(i,j,k) &
                      + g_phys3_ll(i,j,k,jxx) &
                       *shift_u(i,j,k,jx) &
                       *shift_u(i,j,k,jx) &
                      + two*g_phys3_ll(i,j,k,jxy) &
                       *shift_u(i,j,k,jx) &
                       *shift_u(i,j,k,jy) &
                      + two*g_phys3_ll(i,j,k,jxz) &
                       *shift_u(i,j,k,jx) &
                       *shift_u(i,j,k,jz) &
                      + g_phys3_ll(i,j,k,jyy) &
                       *shift_u(i,j,k,jy) &
                       *shift_u(i,j,k,jy) &
                      + two*g_phys3_ll(i,j,k,jyz) &
                       *shift_u(i,j,k,jy) &
                       *shift_u(i,j,k,jz) &
                      + g_phys3_ll(i,j,k,jzz) &
                       *shift_u(i,j,k,jz) &
                       *shift_u(i,j,k,jz)

    g4(i,j,k,itx)=   g_phys3_ll(i,j,k,jxx)*shift_u(i,j,k,jx) &
                      + g_phys3_ll(i,j,k,jxy)*shift_u(i,j,k,jy) &
                      + g_phys3_ll(i,j,k,jxz)*shift_u(i,j,k,jz)

    g4(i,j,k,ity)=   g_phys3_ll(i,j,k,jxy)*shift_u(i,j,k,jx) &
                      + g_phys3_ll(i,j,k,jyy)*shift_u(i,j,k,jy) &
                      + g_phys3_ll(i,j,k,jyz)*shift_u(i,j,k,jz)

    g4(i,j,k,itz)=   g_phys3_ll(i,j,k,jxz)*shift_u(i,j,k,jx) &
                      + g_phys3_ll(i,j,k,jyz)*shift_u(i,j,k,jy) &
                      + g_phys3_ll(i,j,k,jzz)*shift_u(i,j,k,jz)

    g4(i,j,k,ixx)= g_phys3_ll(i,j,k,jxx)
    g4(i,j,k,ixy)= g_phys3_ll(i,j,k,jxy)
    g4(i,j,k,ixz)= g_phys3_ll(i,j,k,jxz)
    g4(i,j,k,iyy)= g_phys3_ll(i,j,k,jyy)
    g4(i,j,k,iyz)= g_phys3_ll(i,j,k,jyz)
    g4(i,j,k,izz)= g_phys3_ll(i,j,k,jzz)

  END SUBROUTINE compute_g4


  SUBROUTINE determinant_sym4x4_grid( i, j, k, A, det )

    !****************************************************************
    !
    !# Compute the determinant of a \(4\times 4\) symmetric matrix
    !  field, given as a 10-vector, at a given grid point
    !
    !  FT
    !
    !****************************************************************

    USE tensor, ONLY: itt, itx, ity, itz, ixx, ixy, ixz, iyy, iyz, izz, n_sym4x4

    IMPLICIT NONE

    INTEGER:: i, j, k
    !! Indices of the grid point
    INTEGER, DIMENSION(4):: components
    !# Array containing the shape of the \(4\times 4\) symmetric matrix.
    !  The first 3 components are equal to the numbers of grid points
    !  along each axis. The fourth component is equal to the number of
    !  independent components of the 4x4 symmetric matrix.
    DOUBLE PRECISION, INTENT(IN):: A(:,:,:,:)
    !# The \(4\times 4\) symmetric matrix, given as a 10-vector.
    !  The first 3 components run over the numbers of grid points
    !  along each axis. The fourth index runs over the number of
    !  independent components of the \(4\times 4\) symmetric matrix.
    DOUBLE PRECISION, INTENT(OUT):: det
    !! Determinant of the \(4\times 4\) symmetric matrix

    components= SHAPE( A )

    IF( components(4) /= n_sym4x4 )THEN
      PRINT *, "** ERROR in determinant_sym4x4_grid in MODULE utility.", &
               " This subroutine needs a symmetric matrix with 10 components,",&
               " and a ", components, "component matrix was given instead."
      STOP
    ENDIF

    det=   A(i,j,k,itt)*(A(i,j,k,ixx)*(A(i,j,k,iyy)*A(i,j,k,izz) &
         - A(i,j,k,iyz)*A(i,j,k,iyz)) &
         + A(i,j,k,ixy)*(A(i,j,k,iyz)*A(i,j,k,ixz) &
         - A(i,j,k,ixy)*A(i,j,k,izz)) &
         + A(i,j,k,ixz)*(A(i,j,k,ixy)*A(i,j,k,iyz) &
         - A(i,j,k,iyy)*A(i,j,k,ixz))) &
         - A(i,j,k,itx)*(A(i,j,k,itx)*(A(i,j,k,iyy)*A(i,j,k,izz) &
         - A(i,j,k,iyz)*A(i,j,k,iyz)) &
         + A(i,j,k,ixy)*(A(i,j,k,iyz)*A(i,j,k,itz) &
         - A(i,j,k,ity)*A(i,j,k,izz)) &
         + A(i,j,k,ixz)*(A(i,j,k,ity)*A(i,j,k,iyz) &
         - A(i,j,k,iyy)*A(i,j,k,itz))) &
         + A(i,j,k,ity)*(A(i,j,k,itx)*(A(i,j,k,ixy)*A(i,j,k,izz) &
         - A(i,j,k,iyz)*A(i,j,k,ixz)) &
         + A(i,j,k,ixx)*(A(i,j,k,iyz)*A(i,j,k,itz) &
         - A(i,j,k,ity)*A(i,j,k,izz)) &
         + A(i,j,k,ixz)*(A(i,j,k,ity)*A(i,j,k,ixz) &
         - A(i,j,k,ixy)*A(i,j,k,itz))) &
         - A(i,j,k,itz)*(A(i,j,k,itx)*(A(i,j,k,ixy)*A(i,j,k,iyz) &
         - A(i,j,k,iyy)*A(i,j,k,ixz)) &
         + A(i,j,k,ixx)*(A(i,j,k,iyy)*A(i,j,k,itz) &
         - A(i,j,k,ity)*A(i,j,k,iyz)) &
         + A(i,j,k,ixy)*(A(i,j,k,ity)*A(i,j,k,ixz) &
         - A(i,j,k,ixy)*A(i,j,k,itz)))

  END SUBROUTINE determinant_sym4x4_grid


  SUBROUTINE determinant_sym3x3_grid( i, j, k, A, det )

    !****************************************************************
    !
    !# Compute the determinant of a \(3\times 3\) symmetric matrix
    !  field, given as a 6-vector, at a given grid point
    !
    !  FT 26.03.2021
    !
    !****************************************************************

    USE tensor, ONLY: jxx, jxy, jxz, jyy, jyz, jzz, n_sym3x3

    IMPLICIT NONE

    INTEGER:: i, j, k
    !! Indices of the grid point
    INTEGER, DIMENSION(4):: components
    !# Array containing the shape of the \(3\times 3\) symmetric matrix.
    !  The first 3 components are equal to the numbers of grid points
    !  along each axis. The fourth component is equal to the number of
    !  independent components of the \(3\times 3\) symmetric matrix.
    DOUBLE PRECISION, INTENT(IN):: A(:,:,:,:)
    !# The \(3\times 3\) symmetric matrix, given as a 6-vector.
    !  The first 3 components run over the numbers of grid points
    !  along each axis. The fourth index runs over the number of
    !  independent components of the \(3\times 3\) symmetric matrix.
    DOUBLE PRECISION, INTENT(OUT):: det
    !! Determinant of the \(3\times 3\) symmetric matrix

    components= SHAPE( A )

    IF( components(4) /= n_sym3x3 )THEN
      PRINT *, "** ERROR in determinant_sym3x3_grid in MODULE utility.", &
               " This subroutine needs a symmetric matrix with 6 components,",&
               " and a ", components, "component matrix was given instead."
      STOP
    ENDIF

    det=   A(i,j,k,jxx)*A(i,j,k,jyy)*A(i,j,k,jzz) &
         + A(i,j,k,jxy)*A(i,j,k,jyz)*A(i,j,k,jxz) &
         + A(i,j,k,jxz)*A(i,j,k,jxy)*A(i,j,k,jyz) &
         - A(i,j,k,jxy)*A(i,j,k,jxy)*A(i,j,k,jzz) &
         - A(i,j,k,jxz)*A(i,j,k,jyy)*A(i,j,k,jxz) &
         - A(i,j,k,jxx)*A(i,j,k,jyz)*A(i,j,k,jyz)

  END SUBROUTINE determinant_sym3x3_grid


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
