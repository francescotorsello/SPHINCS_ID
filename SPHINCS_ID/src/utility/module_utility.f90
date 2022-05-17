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


  USE matrix,     ONLY: determinant_4x4_matrix
  USE constants,  ONLY: G_Msun, c_light2, MSun


  IMPLICIT NONE

  ! This should go in MODULE constants
  DOUBLE PRECISION, PARAMETER :: zero            = 0.D0
  DOUBLE PRECISION, PARAMETER :: one             = 1.D0
  DOUBLE PRECISION, PARAMETER :: two             = 2.D0
  DOUBLE PRECISION, PARAMETER :: three           = 3.D0
  DOUBLE PRECISION, PARAMETER :: four            = 4.D0
  DOUBLE PRECISION, PARAMETER :: five            = 5.D0
  DOUBLE PRECISION, PARAMETER :: ten             = 10.D0
  DOUBLE PRECISION, PARAMETER :: golden_ratio    = 1.618033988749894D0
  DOUBLE PRECISION, PARAMETER :: km2m            = 1.D+3
  DOUBLE PRECISION, PARAMETER :: m2cm            = 1D+2
  DOUBLE PRECISION, PARAMETER :: g2kg            = 1D-3
  DOUBLE PRECISION, PARAMETER :: kg2g            = 1D+3
  DOUBLE PRECISION, PARAMETER :: MSun_geo        = G_Msun/c_light2/1.0D5! in km;
  !# Msun_geo = 1.47662503825040 km
  !  see https://einsteintoolkit.org/thornguide/EinsteinBase/HydroBase/documentation.html
  DOUBLE PRECISION, PARAMETER :: km2Msun_geo     = 1.0D0/MSun_geo
  DOUBLE PRECISION, PARAMETER :: lorene2hydrobase= &
                                              (MSun_geo*km2m)**3/(MSun*g2kg)
  !# lorene2hydrobase= (1477m)^3 / (2*10^30kg) ! 1.6110591665D-21
  !  new value 1.6186541582311746851140226630074e-21, different by 0.5%
  !  lorene2hydrobase is the conversion factor for the baryon mass density

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


  INTEGER, PARAMETER:: max_length= 50
  !! Maximum length for strings
  INTEGER, PARAMETER:: max_n_id= 50
  ! Maximum number of physical systems
  INTEGER, PARAMETER:: max_n_parts= 250
  !! Maximum number of particle distributions

  INTEGER:: n_id
  !! Number of physical systems to set up
  INTEGER:: ref_lev
  !! Number of refinement levels
  INTEGER:: constraints_step
  !! Export the constraints every constraints_step-th step

  INTEGER, PARAMETER:: test_int= - 112
  INTEGER, DIMENSION( max_n_id, max_n_parts ):: placer= test_int
  !# Matrix storing the information on how to place particles for each bns
  !  object. Row i contains information about the i^th bns object.


  DOUBLE PRECISION:: numerator_ratio_dx
  !# Numerator of the rational ratio between the large grid spacing and the
  !  medium one,equal to the ratio between the medium grid spacing nd the small
  !  one. Not used in this PROGRAM, but needed since the PROGRAM reads the same
  !  parameter file as the convergence_test PROGRAM
  DOUBLE PRECISION:: denominator_ratio_dx
  !# Denominator of the rational ratio between the large grid spacing and the
  !  medium one,equal to the ratio between the medium grid spacing nd the small
  !  one. Not used in this PROGRAM, but needed since the PROGRAM reads the same
  !  parameter file as the convergence_test PROGRAM

  ! Logical variables to steer the execution
  LOGICAL:: export_bin, export_form, export_form_xy, export_form_x, &
            compute_constraints, export_constraints_xy, &
            export_constraints_x, export_constraints, &
            export_constraints_details, compute_parts_constraints, &
            one_lapse, zero_shift, run_sph, run_spacetime, estimate_length_scale

  CHARACTER( LEN= max_length ), DIMENSION( max_length ):: filenames= "0"
  !! Array of strings storing the names of the |id| files
  CHARACTER( LEN= max_length ):: common_path
  !# String storing the local path to the directory where the |id| files
  !  are stored
  CHARACTER( LEN= max_length ):: sph_path
  !# String storing the local path to the directory where the
  !  SPH output is to be saved
  CHARACTER( LEN= max_length ):: spacetime_path
  !# String storing the local path to the directory where the
  !  spacetime output is to be saved


  CONTAINS


  SUBROUTINE read_sphincs_id_parameters()

    !***********************************************
    !
    !# Read the parameters to steer SPHINCS_ID
    !
    !  FT
    !
    !***********************************************

    IMPLICIT NONE

    INTEGER:: stat
    INTEGER, PARAMETER:: unit_parameters= 17

    CHARACTER( LEN= : ), ALLOCATABLE:: sphincs_id_parameters_namefile
    CHARACTER( LEN= 100 ):: msg

    ! Namelist containing parameters read from sphincs_id_parameters.par
    ! by the SUBROUTINE read_sphincs_id_parameters of this PROGRAM
    NAMELIST /sphincs_id_parameters/ &
              n_id, common_path, filenames, placer, &
              export_bin, export_form, export_form_xy, &
              export_form_x, export_constraints_xy, &
              export_constraints_x, compute_constraints, &
              export_constraints, export_constraints_details, &
              constraints_step, compute_parts_constraints, &
              numerator_ratio_dx, denominator_ratio_dx, ref_lev, &
              one_lapse, zero_shift, show_progress, &
              run_sph, run_spacetime, sph_path, spacetime_path, &
              estimate_length_scale

    sphincs_id_parameters_namefile= 'sphincs_id_parameters.dat'

    INQUIRE( FILE= sphincs_id_parameters_namefile, EXIST= file_exists )
    IF( file_exists )THEN

     OPEN( 17, FILE= sphincs_id_parameters_namefile, STATUS= 'OLD' )

    ELSE

     PRINT*
     PRINT*,'** ERROR: ', sphincs_id_parameters_namefile, " file not found!"
     PRINT*
     STOP

    ENDIF

    READ( UNIT= unit_parameters, NML= sphincs_id_parameters, IOSTAT= stat, &
          IOMSG= msg )

      IF( stat /= 0 )THEN
        PRINT *, "** ERROR: Error in reading ", sphincs_id_parameters_namefile,&
                 ". The IOSTAT variable is ", stat, &
                 "The error message is", msg
        STOP
      ENDIF

    CLOSE( UNIT= unit_parameters )

    DO itr= 1, max_length, 1
      IF( TRIM(filenames(itr)).NE."0" )THEN
        cnt= cnt + 1
      ENDIF
    ENDDO
    IF( cnt.NE.n_id )THEN
      PRINT *, "** ERROR! The number of file names is", cnt, &
               "and n_id=", n_id, ". The two should be the same."
      PRINT *
      STOP
    ENDIF

   !DO itr= 1, n_id, 1
   !  DO itr2= 1, max_n_parts, 1
   !    IF( placer( itr, itr2 ) == test_int )THEN
   !      PRINT *
   !      PRINT *, "** ERROR! The array placer does not have ", &
   !               "enough components to specify all the desired ", &
   !               "particle distributions. Specify the ", &
   !               "components in file lorene_bns_id_particles.par"
   !      PRINT *
   !      STOP
   !    ENDIF
   !  ENDDO
   !ENDDO

  END SUBROUTINE read_sphincs_id_parameters


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


  SUBROUTINE compute_tpo_metric( g4, lapse, shift, g3 )

    !***********************************************
    !
    !# Computes the lapse,shift and spatial metric
    !  from the covariant spacetime metric
    !
    !  FT 12.04.2022
    !
    !***********************************************

    USE matrix,     ONLY: invert_3x3_matrix
    USE tensor,     ONLY: itt, itx, ity, itz, ixx, ixy, &
                          ixz, iyy, iyz, izz, jxx, jxy, jxz, &
                          jyy, jyz, jzz, jx, jy, jz, n_sym4x4, n_sym3x3

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(n_sym4x4),  INTENT(IN) :: g4
    !! Covariant spacetime metric
    DOUBLE PRECISION,                       INTENT(OUT):: lapse
    !! Lapse function
    DOUBLE PRECISION, DIMENSION(3),         INTENT(OUT):: shift
    !! Contravariant shift vector
    DOUBLE PRECISION, DIMENSION(n_sym3x3),  INTENT(OUT):: g3
    !! Covariant spatial metric

    DOUBLE PRECISION, DIMENSION(3,3):: gmat
    DOUBLE PRECISION, DIMENSION(3,3):: gmat_inv

    g3(jxx)= g4(ixx)
    g3(jxy)= g4(ixy)
    g3(jxz)= g4(ixz)
    g3(jyy)= g4(iyy)
    g3(jyz)= g4(iyz)
    g3(jzz)= g4(izz)

    gmat(1,1)= g3(jxx)
    gmat(1,2)= g3(jxy)
    gmat(1,3)= g3(jxz)
    gmat(2,1)= g3(jxy)
    gmat(2,2)= g3(jyy)
    gmat(2,3)= g3(jyz)
    gmat(3,1)= g3(jxz)
    gmat(3,2)= g3(jyz)
    gmat(3,3)= g3(jzz)

    CALL invert_3x3_matrix( gmat, gmat_inv )

    shift(jx)= gmat_inv(jx,jx)*g4(itx) + gmat_inv(jx,jy)*g4(ity) &
             + gmat_inv(jx,jz)*g4(itz)
    shift(jy)= gmat_inv(jy,jx)*g4(itx) + gmat_inv(jy,jy)*g4(ity) &
             + gmat_inv(jy,jz)*g4(itz)
    shift(jz)= gmat_inv(jz,jx)*g4(itx) + gmat_inv(jz,jy)*g4(ity) &
             + gmat_inv(jz,jz)*g4(itz)

    lapse= SQRT( shift(jx)*g4(itx) + shift(jy)*g4(ity) + shift(jz)*g4(itz) &
                 - g4(itt) )

  END SUBROUTINE compute_tpo_metric


  SUBROUTINE invert_sym4x4( A, iA )

    !****************************************************************
    !
    !# Invert a \(4\times 4\) symemtric matrix stored as a \(10\)-vector
    !  @note The inverse (and the adjugate) of a symmetric matrix
    !        is symmetric
    !
    !  FT 25.04.2022
    !
    !****************************************************************

    USE tensor,               ONLY: n_sym4x4
    USE metric_on_particles,  ONLY: gvec2mat, mat2gvec
    USE matrix,               ONLY: invert_4x4_matrix

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: A(:)
    !# The \(4\times 4\) symmetric matrix, given as a 10-vector.
    !  The first 3 components run over the numbers of grid points
    !  along each axis. The fourth index runs over the number of
    !  independent components of the \(4\times 4\) symmetric matrix.
    DOUBLE PRECISION, INTENT(OUT):: iA(n_sym4x4)
    !# Inverse of the \(4\times 4\) symmetric matrix, given as input.

    DOUBLE PRECISION:: Amat(4,4)
    !! The \(4\times 4\) symmetric matrix as a matrix
    DOUBLE PRECISION:: iAmat(4,4)
    !! The inverse of the \(4\times 4\) symmetric matrix, as a matrix

    IF( SIZE(A) /= n_sym4x4 )THEN
      PRINT *, "** ERROR in determinant_sym4x4_grid in MODULE utility.", &
               " This subroutine needs a symmetric matrix with 10 components,",&
               " and a ", SIZE(A), "component matrix was given instead."
      STOP
    ENDIF

    CALL gvec2mat(A,Amat)

    CALL invert_4x4_matrix(Amat,iAmat)

    CALL mat2gvec(iA,iAmat)

  END SUBROUTINE invert_sym4x4


  SUBROUTINE invert_sym3x3( A, iA )

    !****************************************************************
    !
    !# Invert a \(3\times 3\) symemtric matrix stored as a \(6\)-vector
    !  @note The inverse (and the adjugate) of a symmetric matrix
    !        is symmetric
    !
    !  FT 25.04.2022
    !
    !****************************************************************

    USE tensor, ONLY: n_sym3x3
    USE matrix, ONLY: invert_3x3_matrix

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: A(:)
    !# The \(3\times 3\) symmetric matrix, given as a \(6\)-vector.
    DOUBLE PRECISION, INTENT(OUT):: iA(n_sym3x3)
    !# Inverse of the \(3\times 3\) symmetric matrix, given as input.

    DOUBLE PRECISION:: Amat(3,3)
    !! The \(3\times 3\) symmetric matrix as a matrix
    DOUBLE PRECISION:: iAmat(3,3)
    !! The inverse of the \(3\times 3\) symmetric matrix, as a matrix

    IF( SIZE(A) /= n_sym3x3 )THEN
      PRINT *, "** ERROR in determinant_sym3x3_grid in MODULE utility.", &
               " This subroutine needs a symmetric matrix with 6 components,",&
               " and a ", SIZE(A), "component matrix was given instead."
      STOP
    ENDIF

    CALL vec2mat_sym3x3(A,Amat)

    CALL invert_3x3_matrix(Amat,iAmat)

    CALL mat2vec_sym3x3(iAmat,iA)

  END SUBROUTINE invert_sym3x3


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


  SUBROUTINE spacetime_vector_norm_sym4x4( g4, v, norm )

    !****************************************************************
    !
    !# Compute the spacetime squared norm of a vector, using the
    !  metric given as an array of 10 components
    !
    !  FT 07.02.2022
    !
    !****************************************************************

    USE tensor,    ONLY: itt, itx, ity, itz, ixx, ixy, ixz, iyy, iyz, izz, &
                         it, ix, iy, iz, n_sym4x4

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: g4(itt:izz)
    !# The \(4\times 4\) spacetime metric, given as a 10-vector.
    DOUBLE PRECISION, INTENT(IN):: v(it:iz)
    !# The \(4\)-vector whose norm has to be computed.
    DOUBLE PRECISION, INTENT(OUT):: norm
    !! Spacetime norm of the vector v.

    IF( SIZE(g4) /= n_sym4x4 )THEN
      PRINT *, "** ERROR in spacetime_vector_norm_sym4x4 in MODULE utility.", &
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


  SUBROUTINE spatial_vector_norm_sym3x3( g3, v, norm )

    !****************************************************************
    !
    !# Compute the spatial squared norm of a vector, using the
    !  spatial metric given as an array of 6 components
    !
    !  FT 14.02.2022
    !
    !****************************************************************

    USE tensor,    ONLY: jxx, jxy, jxz, jyy, jyz, jzz, jx, jy, jz, n_sym3x3

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: g3(jxx:jzz)
    !# The \(3\times 3\) spacetime metric, given as a 6-vector.
    DOUBLE PRECISION, INTENT(IN):: v(jx:jz)
    !# The \(3\)-vector whose norm has to be computed.
    DOUBLE PRECISION, INTENT(OUT):: norm
    !! Spatial norm of the vector v.

    IF( SIZE(g3) /= n_sym3x3 )THEN
      PRINT *, "** ERROR in spatial_vector_norm_sym3x3 in MODULE utility.", &
               " This subroutine needs a symmetric matrix with 6 components,",&
               " and a ", SIZE(g3), "component matrix was given instead."
      STOP
    ENDIF

    norm= g3(jxx)*v(jx)*v(jx)     + two*g3(jxy)*v(jx)*v(jy) &
        + two*g3(jxz)*v(jx)*v(jz) + g3(jyy)*v(jy)*v(jy)     &
        + two*g3(jyz)*v(jy)*v(jz) + g3(jzz)*v(jz)*v(jz)

  END SUBROUTINE spatial_vector_norm_sym3x3


  PURE SUBROUTINE vec2mat_sym3x3( vec, mat )

    !********************************************
    !
    !# Write the components of symmetric \(3\times 3\)
    !  matrix given as a \(6\)-vector, into a
    !  \(3\times 3\) matrix
    !
    ! FT 25.04.2022
    !
    !*********************************************

    USE tensor,    ONLY: jxx, jxy, jxz, jyy, jyz, jzz, jx, jy, jz, n_sym3x3

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN)  :: vec(n_sym3x3)
    DOUBLE PRECISION, INTENT(OUT) :: mat(3,3)

    mat(1,1)= vec(jxx)
    mat(1,2)= vec(jxy)
    mat(1,3)= vec(jxz)

    mat(2,1)= mat(1,2)
    mat(2,2)= vec(jyy)
    mat(2,3)= vec(jyz)

    mat(3,1)= mat(1,3)
    mat(3,2)= mat(2,3)
    mat(3,3)= vec(jzz)

  END SUBROUTINE vec2mat_sym3x3


  PURE SUBROUTINE mat2vec_sym3x3( mat, vec )

    !************************************
    !                                   *
    ! transform symmetric 4x4-matrix    *
    ! into vector; SKR 30.11.2017       *
    !                                   *
    !************************************

    USE tensor,    ONLY: jxx, jxy, jxz, jyy, jyz, jzz, jx, jy, jz, n_sym3x3

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN)  :: mat(3,3)
    DOUBLE PRECISION, INTENT(OUT) :: vec(n_sym3x3)

    vec(jxx)= mat(1,1)
    vec(jxy)= mat(1,2)
    vec(jxz)= mat(1,3)
    vec(jyy)= mat(2,2)
    vec(jyz)= mat(2,3)
    vec(jzz)= mat(3,3)

  END SUBROUTINE mat2vec_sym3x3


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

    USE constants, ONLY: pi

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


  FUNCTION k_lorene2hydrobase( gamma )

    DOUBLE PRECISION :: gamma
    DOUBLE PRECISION :: k_lorene2hydrobase

    ! LORENE's EOS is in terms on number density n = rho/m_nucleon:
    ! P = K n^Gamma
    ! to convert to SI units:
    ! K_SI(n) = K_LORENE rho_nuc c^2 / n_nuc^gamma
    ! Converting this to be in terms of the mass density rho = n m_nucleon gets
    ! changes n_nuc to rho_nuc:
    ! K_SI(rho) = K_LORENE c^2 / rho_nuc^(gamma-1)
    ! In SI units P has units of M / (L T^2) and rho has units of M/L^3 thus
    ! K_SI has units of (L^3/M)^Gamma M/(L T^2).
    ! In Cactus units P and rho have the same units thus K_Cactus is unitless.
    ! Conversion between K_SI and K_Cactus thus amounts to dividing out the
    ! units of the SI quantity.

    !k_lorene2hydrobase= ((c_light*cm2m)**6/ &
    !              ( (G_grav*(cm2m**3)/(g2kg))**3*(MSun*g2kg)**2*(1.66E+17) )) &
    !              **( gamma - 1.0D0 )

    ! Our testbed cases are gamma= 2.75, k= 30000; and gamma=2, k=100
    ! in SPHINCS units

    k_lorene2hydrobase= &
                        ( (MSun*g2kg)/((MSun_geo*km2m)**3*(1.66D+17)) ) &
                        **( gamma - 1.0D0 )


  END FUNCTION


  FUNCTION k_lorene2hydrobase_piecewisepolytrope( gamma0 )

    DOUBLE PRECISION :: gamma0
    DOUBLE PRECISION :: k_lorene2hydrobase_piecewisepolytrope

    ! LORENE has K0 in units of (g cm^{-3})^{1-gamma0} for the piecewise
    ! polytropes. This factor writes it in SPHINCS units

    k_lorene2hydrobase_piecewisepolytrope= &
                        ( MSun/((MSun_geo*km2m*m2cm)**3) ) &
                        **( gamma0 - 1.0D0 )


  END FUNCTION


END MODULE utility
