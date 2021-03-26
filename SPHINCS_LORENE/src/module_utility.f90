! File:         module_utility.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

MODULE utility

  !***********************************************************************
  !                                                                      *
  !  This module contains useful miscellaneous PROCEDURES and variables  *
  !                                                                      *
  !***********************************************************************


  USE matrix, ONLY: determinant_4x4_matrix


  IMPLICIT NONE


  INTEGER:: itr, itr3, itr4     ! iterators for loops
  INTEGER:: ios                 ! variable to store the state of I/O
  INTEGER:: cnt= 0              ! counter

  ! Variables to print progress on screen
  INTEGER:: perc
  DOUBLE PRECISION:: perc2
  CHARACTER, PARAMETER:: creturn= ACHAR(13)   ! Carriage return

  LOGICAL:: file_exists, show_progress

  ! String storing error messages
  CHARACTER( LEN= : ), ALLOCATABLE:: err_msg

  ! Variables used to set the run_id
  CHARACTER(8)  :: date
  CHARACTER(10) :: time
  CHARACTER(5)  :: zone
  INTEGER, DIMENSION(8) :: values
  CHARACTER( LEN= 19 ):: run_id


  CONTAINS


  SUBROUTINE test_status( io_stat, io_msg, opt_msg )

    !************************************************
    !                                               *
    ! Test if a status variable is 0 or not         *
    !                                               *
    ! FT 17.09.2020                                 *
    !                                               *
    !************************************************

    IMPLICIT NONE

    INTEGER,               INTENT(IN)           :: io_stat
    CHARACTER( LEN= 100 ), INTENT(IN)           :: io_msg
    CHARACTER( LEN= * ),   INTENT(IN), OPTIONAL :: opt_msg

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

  SUBROUTINE compute_g4( ix, iy, iz, lapse, shift_u, g_phys3_ll, g4 )

    !************************************************
    !                                               *
    ! Computes the spacetime metric from lapse,     *
    ! shift and spatial metric
    !                                               *
    ! FT 27.11.2020                                 *
    !                                               *
    !************************************************

    USE tensor,     ONLY: itt, itx, ity, itz, ixx, ixy, &
                          ixz, iyy, iyz, izz, jxx, jxy, jxz, &
                          jyy, jyz, jzz, jx, jy, jz, n_sym3x3, n_sym4x4

    IMPLICIT NONE

    INTEGER:: ix, iy, iz

    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT( IN OUT ):: &
    lapse
    DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: &
    shift_u
    DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: &
    g_phys3_ll
    DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT( IN OUT ):: g4

    g4(ix,iy,iz,itt)= - lapse(ix,iy,iz)*lapse(ix,iy,iz) &
                      + g_phys3_ll(ix,iy,iz,jxx) &
                       *shift_u(ix,iy,iz,jx) &
                       *shift_u(ix,iy,iz,jx) &
                      + 2.0D0*g_phys3_ll(ix,iy,iz,jxy) &
                       *shift_u(ix,iy,iz,jx) &
                       *shift_u(ix,iy,iz,jy) &
                      + 2.0D0*g_phys3_ll(ix,iy,iz,jxz) &
                       *shift_u(ix,iy,iz,jx) &
                       *shift_u(ix,iy,iz,jz) &
                      + g_phys3_ll(ix,iy,iz,jyy) &
                       *shift_u(ix,iy,iz,jy) &
                       *shift_u(ix,iy,iz,jy) &
                      + 2.0D0*g_phys3_ll(ix,iy,iz,jyz) &
                       *shift_u(ix,iy,iz,jy) &
                       *shift_u(ix,iy,iz,jz) &
                      + g_phys3_ll(ix,iy,iz,jzz) &
                       *shift_u(ix,iy,iz,jz) &
                       *shift_u(ix,iy,iz,jz)

    g4(ix,iy,iz,itx)=   g_phys3_ll(ix,iy,iz,jxx)*shift_u(ix,iy,iz,jx) &
                      + g_phys3_ll(ix,iy,iz,jxy)*shift_u(ix,iy,iz,jy) &
                      + g_phys3_ll(ix,iy,iz,jxz)*shift_u(ix,iy,iz,jz)

    g4(ix,iy,iz,ity)=   g_phys3_ll(ix,iy,iz,jxy)*shift_u(ix,iy,iz,jx) &
                      + g_phys3_ll(ix,iy,iz,jyy)*shift_u(ix,iy,iz,jy) &
                      + g_phys3_ll(ix,iy,iz,jyz)*shift_u(ix,iy,iz,jz)

    g4(ix,iy,iz,itz)=   g_phys3_ll(ix,iy,iz,jxz)*shift_u(ix,iy,iz,jx) &
                      + g_phys3_ll(ix,iy,iz,jyz)*shift_u(ix,iy,iz,jy) &
                      + g_phys3_ll(ix,iy,iz,jzz)*shift_u(ix,iy,iz,jz)

    g4(ix,iy,iz,ixx)= g_phys3_ll(ix,iy,iz,jxx)
    g4(ix,iy,iz,ixy)= g_phys3_ll(ix,iy,iz,jxy)
    g4(ix,iy,iz,ixz)= g_phys3_ll(ix,iy,iz,jxz)
    g4(ix,iy,iz,iyy)= g_phys3_ll(ix,iy,iz,jyy)
    g4(ix,iy,iz,iyz)= g_phys3_ll(ix,iy,iz,jyz)
    g4(ix,iy,iz,izz)= g_phys3_ll(ix,iy,iz,jzz)

  END SUBROUTINE compute_g4

  SUBROUTINE determinant_sym4x4_grid( ix, iy, iz, A, det )

    !*****************************************************************
    !                                                                *
    ! Compute the determinant of a 4x4 symmetric matrix field at a   *
    ! given grid point                                               *
    !                                                                *
    !*****************************************************************

    USE tensor, ONLY: itt, itx, ity, itz, ixx, ixy, ixz, iyy, iyz, izz, n_sym4x4

    IMPLICIT NONE

    INTEGER:: ix, iy, iz
    INTEGER, DIMENSION(4):: components
    DOUBLE PRECISION, INTENT(IN):: A(:,:,:,:)
    DOUBLE PRECISION, INTENT(OUT):: det

    components= SHAPE( A )

    IF( components(4) /= n_sym4x4 )THEN
      PRINT *, "** ERROR in determinant_sym4x4_grid in MODULE utility.", &
               " This subroutine needs a symmetric matrix with 10 components,",&
               " and a ", components, "component matrix was given instead."
      STOP
    ENDIF

    det=   A(ix,iy,iz,itt)*(A(ix,iy,iz,ixx)*(A(ix,iy,iz,iyy)*A(ix,iy,iz,izz) &
         - A(ix,iy,iz,iyz)*A(ix,iy,iz,iyz)) &
         + A(ix,iy,iz,ixy)*(A(ix,iy,iz,iyz)*A(ix,iy,iz,ixz) &
         - A(ix,iy,iz,ixy)*A(ix,iy,iz,izz)) &
         + A(ix,iy,iz,ixz)*(A(ix,iy,iz,ixy)*A(ix,iy,iz,iyz) &
         - A(ix,iy,iz,iyy)*A(ix,iy,iz,ixz))) &
         - A(ix,iy,iz,itx)*(A(ix,iy,iz,itx)*(A(ix,iy,iz,iyy)*A(ix,iy,iz,izz) &
         - A(ix,iy,iz,iyz)*A(ix,iy,iz,iyz)) &
         + A(ix,iy,iz,ixy)*(A(ix,iy,iz,iyz)*A(ix,iy,iz,itz) &
         - A(ix,iy,iz,ity)*A(ix,iy,iz,izz)) &
         + A(ix,iy,iz,ixz)*(A(ix,iy,iz,ity)*A(ix,iy,iz,iyz) &
         - A(ix,iy,iz,iyy)*A(ix,iy,iz,itz))) &
         + A(ix,iy,iz,ity)*(A(ix,iy,iz,itx)*(A(ix,iy,iz,ixy)*A(ix,iy,iz,izz) &
         - A(ix,iy,iz,iyz)*A(ix,iy,iz,ixz)) &
         + A(ix,iy,iz,ixx)*(A(ix,iy,iz,iyz)*A(ix,iy,iz,itz) &
         - A(ix,iy,iz,ity)*A(ix,iy,iz,izz)) &
         + A(ix,iy,iz,ixz)*(A(ix,iy,iz,ity)*A(ix,iy,iz,ixz) &
         - A(ix,iy,iz,ixy)*A(ix,iy,iz,itz))) &
         - A(ix,iy,iz,itz)*(A(ix,iy,iz,itx)*(A(ix,iy,iz,ixy)*A(ix,iy,iz,iyz) &
         - A(ix,iy,iz,iyy)*A(ix,iy,iz,ixz)) &
         + A(ix,iy,iz,ixx)*(A(ix,iy,iz,iyy)*A(ix,iy,iz,itz) &
         - A(ix,iy,iz,ity)*A(ix,iy,iz,iyz)) &
         + A(ix,iy,iz,ixy)*(A(ix,iy,iz,ity)*A(ix,iy,iz,ixz) &
         - A(ix,iy,iz,ixy)*A(ix,iy,iz,itz)))

  END SUBROUTINE determinant_sym4x4_grid

  SUBROUTINE determinant_sym3x3_grid( i, j, k, A, det )

    !*****************************************************************
    !                                                                *
    ! Compute the determinant of a 3x3 symmetric matrix field at a   *
    ! given grid point                                               *
    !                                                                *
    ! FT 26.03.2021                                                  *
    !                                                                *
    !*****************************************************************


    USE tensor, ONLY: jxx, jxy, jxz, jyy, jyz, jzz, n_sym3x3

    IMPLICIT NONE

    INTEGER:: i, j, k
    INTEGER, DIMENSION(4):: components
    DOUBLE PRECISION, INTENT(IN):: A(:,:,:,:)
    DOUBLE PRECISION, INTENT(OUT):: det

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


END MODULE utility
