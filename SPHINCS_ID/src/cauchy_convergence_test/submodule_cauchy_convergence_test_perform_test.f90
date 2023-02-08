! File:         submodule_cauchy_convergence_test_perform_test.f90
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

SUBMODULE (cauchy_convergence_test) perform_test

  !********************************************
  !
  !# This submodule contains the implementation
  !  of the PROCEDURES in MODULE
  !  cauchy_convergence_test that perform
  !  the Cauchy convergence test
  !
  !  FT 22.09.2022
  !
  !********************************************


  USE tensor,  ONLY: jx, jy, jz
  USE utility, ONLY: one


  IMPLICIT NONE


  DOUBLE PRECISION, PARAMETER:: tiny_real= 1D-30
  LOGICAL,          PARAMETER:: debug= .FALSE.


  INTERFACE
    MODULE SUBROUTINE get_scalar_at_grid_point( tpof, i, j, k, l, scalar )
    !! Returns the value of a scalar field at the desired point
      CLASS(tpo), INTENT(INOUT):: tpof
      INTEGER, INTENT(IN):: i
      !! \(x\) index of the desired point
      INTEGER, INTENT(IN):: j
      !! \(y\) index of the desired point
      INTEGER, INTENT(IN):: k
      !! \(z\) index of the desired point
      INTEGER, INTENT(IN):: l
      !! Index of the refinement level
      DOUBLE PRECISION, INTENT(OUT):: scalar
      !! Value of the scalar field at \((i,j,k)\)
    END SUBROUTINE get_scalar_at_grid_point
  END INTERFACE
  PROCEDURE(get_scalar_at_grid_point), POINTER:: get_ham_ptr

  INTERFACE compute_convergence_factor
  !# Generic PROCEDURE to compute the Cauchy convergence factor


    MODULE SUBROUTINE compute_convergence_factor_unknown_sol &
      ( nx, ny, nz, num, den, ref_lev, tpo_coarse, tpo_medium, tpo_fine, &
        get_hc, convergence_factor )
    !# Compute the Cauchy convergence factor when the exact solution is not
    !  known.

      INTEGER,                               INTENT(IN)   :: nx, ny, nz
      DOUBLE PRECISION,                      INTENT(IN)   :: num
      DOUBLE PRECISION,                      INTENT(IN)   :: den
      INTEGER,                               INTENT(IN)   :: ref_lev
      CLASS(tpo),                            INTENT(INOUT):: tpo_coarse
      CLASS(tpo),                            INTENT(INOUT):: tpo_medium
      CLASS(tpo),                            INTENT(INOUT):: tpo_fine
      PROCEDURE(get_scalar_at_grid_point),   POINTER, INTENT(IN):: get_hc
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT):: &
      convergence_factor

    END SUBROUTINE compute_convergence_factor_unknown_sol


    MODULE SUBROUTINE compute_convergence_factor_known_sol &
      ( nx, ny, nz, num, den, ref_lev, tpo_coarse, tpo_fine, &
        get_hc, convergence_factor )
    !# Compute the Cauchy convergence factor when the exact solution is known.

      INTEGER,                               INTENT(IN)   :: nx, ny, nz
      DOUBLE PRECISION,                      INTENT(IN)   :: num
      DOUBLE PRECISION,                      INTENT(IN)   :: den
      INTEGER,                               INTENT(IN)   :: ref_lev
      CLASS(tpo),                            INTENT(INOUT):: tpo_coarse
      CLASS(tpo),                            INTENT(INOUT):: tpo_fine
      PROCEDURE(get_scalar_at_grid_point),   POINTER, INTENT(IN):: get_hc
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT):: &
      convergence_factor

    END SUBROUTINE compute_convergence_factor_known_sol


  END INTERFACE compute_convergence_factor



  CONTAINS



  SUBROUTINE get_ham( tpof, i, j, k, l, hc )
  !# Wrapper SUBROUTINE to get the value of the Hamiltonian constraint
  !  computed using ID read on the refined mesh,
  !  at a given mesh point `(i,j,k)`, on the refinement level `l`, for the
  !  [[tpo]] object `tpof`

    IMPLICIT NONE

    CLASS(tpo), INTENT(INOUT):: tpof
    !! [[tpo]] object to use
    INTEGER, INTENT(IN):: i
    !! \(x\) index of the desired point
    INTEGER, INTENT(IN):: j
    !! \(y\) index of the desired point
    INTEGER, INTENT(IN):: k
    !! \(z\) index of the desired point
    INTEGER, INTENT(IN):: l
    !! Index of the refinement level
    DOUBLE PRECISION, INTENT(OUT):: hc
    !! Value of the Hamiltonian constraint at \((i,j,k)\)

    hc= tpof% get_HC( i, j, k, l )

  END SUBROUTINE get_ham


  SUBROUTINE get_ham_parts( tpof, i, j, k, l, hc )
  !# Wrapper SUBROUTINE to get the value of the Hamiltonian constraint
  !  computed using using the hydro ID mapped from the particles to the refined
  !  mesh, at a given mesh point `(i,j,k)`, on the refinement level `l`, for the
  !  [[tpo]] object `tpof`

    IMPLICIT NONE

    CLASS(tpo), INTENT(INOUT):: tpof
    !! [[tpo]] object to use
    INTEGER, INTENT(IN):: i
    !! \(x\) index of the desired point
    INTEGER, INTENT(IN):: j
    !! \(y\) index of the desired point
    INTEGER, INTENT(IN):: k
    !! \(z\) index of the desired point
    INTEGER, INTENT(IN):: l
    !! Index of the refinement level
    DOUBLE PRECISION, INTENT(OUT):: hc
    !! Value of the Hamiltonian constraint at \((i,j,k)\)

    hc= tpof% get_HC_parts( i, j, k, l )

  END SUBROUTINE get_ham_parts


  MODULE PROCEDURE compute_convergence_factor_unknown_sol
  !# Compute the Cauchy convergence factor when the exact solution is not known.

    IMPLICIT NONE

    INTEGER:: i, j, k
    DOUBLE PRECISION:: ratio_dx, hc_coarse, hc_medium, hc_fine

    PRINT *, "** Computing convergence factor..."

    IF( ALLOCATED(convergence_factor) ) DEALLOCATE(convergence_factor)
    ALLOCATE( convergence_factor( nx, ny, nz ) )

    ratio_dx= num/den

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( tpo_coarse, tpo_medium, tpo_fine, ref_lev, &
    !$OMP                     convergence_factor, den, num, nx, ny, nz, &
    !$OMP                     ratio_dx, get_ham_ptr ) &
    !$OMP             PRIVATE( i, j, k, hc_coarse, hc_medium, hc_fine )
    shared_grid_loop: DO k= 0, nz - 1, 1
      DO j= 0, ny - 1, 1
        DO i= 0, nx - 1, 1

          CALL get_ham_ptr( tpo_coarse, 1 + INT(den)*i, &
                            1 + INT(den)*j, &
                            1 + INT(den)*k, ref_lev, hc_coarse )

          CALL get_ham_ptr( tpo_medium, 1 + INT(den)*i, &
                            1 + INT(den)*j, &
                            1 + INT(den)*k, ref_lev, hc_medium )

          CALL get_ham_ptr( tpo_coarse, 1 + INT(den)*i, &
                            1 + INT(den)*j, &
                            1 + INT(den)*k, ref_lev, hc_fine )

          convergence_factor( 1 + i, 1 + j, 1 + k )= &
           LOG( ABS( ( hc_coarse - hc_medium ) &
                    /( hc_medium - hc_fine + tiny_real ) ) &
                + tiny_real &
           )/LOG(ratio_dx)
           !LOG( ABS( ( ABS(hc_coarse) - ABS(hc_medium) ) &
           !         /( ABS(hc_medium) - ABS(hc_fine) + tiny_real ) ) &
           !     + tiny_real &
           !)/LOG(ratio_dx)

        ENDDO
      ENDDO
    ENDDO shared_grid_loop
    !$OMP END PARALLEL DO
    PRINT *, " * Convergence factor computed."
    PRINT *

  END PROCEDURE compute_convergence_factor_unknown_sol


  MODULE PROCEDURE compute_convergence_factor_known_sol
  !# Compute the Cauchy convergence factor when the exact solution is known.

    IMPLICIT NONE

    INTEGER:: i, j, k
    DOUBLE PRECISION:: ratio_dx, hc_coarse, hc_fine

    PRINT *, "** Computing convergence factor..."

    IF( ALLOCATED(convergence_factor) ) DEALLOCATE(convergence_factor)
    ALLOCATE( convergence_factor( nx, ny, nz ) )

    ratio_dx= num/den

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( tpo_coarse, tpo_fine, ref_lev, &
    !$OMP                     convergence_factor, den, num, nx, ny, nz, &
    !$OMP                     ratio_dx, get_ham_ptr ) &
    !$OMP             PRIVATE( i, j, k, hc_coarse, hc_fine )
    shared_grid_loop: DO k= 0, nz - 1, 1
      DO j= 0, ny - 1, 1
        DO i= 0, nx - 1, 1

          CALL get_ham_ptr( tpo_coarse, 1 + INT(den)*i, &
                            1 + INT(den)*j, &
                            1 + INT(den)*k, ref_lev, hc_coarse )

          CALL get_ham_ptr( tpo_coarse, 1 + INT(den)*i, &
                            1 + INT(den)*j, &
                            1 + INT(den)*k, ref_lev, hc_fine )

          convergence_factor( 1 + i, 1 + j, 1 + k )= &
           LOG( ABS( hc_coarse/(hc_fine + tiny_real) ) + tiny_real ) &
           /LOG(ratio_dx)

        ENDDO
      ENDDO
    ENDDO shared_grid_loop
    !$OMP END PARALLEL DO
    PRINT *, " * Convergence factor computed."
    PRINT *

  END PROCEDURE compute_convergence_factor_known_sol


  MODULE PROCEDURE perform_cauchy_convergence_test_unknown_sol

    !***********************************************************
    !
    !# Perform the Cauchy convergence test when the exact solution is not
    !  known. The ratio between the grid spacings is `num/den`.
    !
    !***********************************************************

    IMPLICIT NONE

    INTEGER:: nx, ny, nz, unit_cauchy_ct

    DOUBLE PRECISION:: ratio_dx

    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: shared_grid
    DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE:: convergence_factor

    CHARACTER( LEN=: ), ALLOCATABLE:: name_cauchy_ct

    ratio_dx= num/den

    CALL find_shared_grid( tpo_coarse, tpo_medium, tpo_fine, num, den, &
                           ref_lev, shared_grid )

    nx= SIZE(shared_grid(:,1,1,jx))
    ny= SIZE(shared_grid(1,:,1,jy))
    nz= SIZE(shared_grid(1,1,:,jz))

    choose_constraints: SELECT CASE( use_constraints )

    CASE(use_constraints_on_mesh)

      get_ham_ptr => get_ham

      unit_cauchy_ct= 3108
      name_cauchy_ct= TRIM(spacetime_path) &
                      //"cauchy_convergence_test_unknown.dat"

    CASE(use_constraints_with_mapped_hydro)

      get_ham_ptr => get_ham_parts

      unit_cauchy_ct= 3110
      name_cauchy_ct= TRIM(spacetime_path) &
                      //"cauchy_convergence_test_unknown_parts.dat"

    CASE DEFAULT

      PRINT *, "** There is no well defined algorithm " &
               // "corresponding to the number", use_constraints
      PRINT *, " * Please set use_constraints to 1 or 2."
      PRINT *, "   1: Use the constraints computed entirely on the mesh"
      PRINT *, "   2: Use the constraints computed on the mesh, using ", &
               "the hydro data mapped from the particles."
      PRINT *
      STOP

    END SELECT choose_constraints

    CALL compute_convergence_factor_unknown_sol( nx, ny, nz, num, den, &
                            ref_lev, tpo_coarse, tpo_medium, tpo_fine, &
                            get_ham_ptr, convergence_factor )

    !
    !-- Print results to formatted file
    !
    CALL print_convergence_factor( nx, ny, nz, shared_grid, convergence_factor,&
                                   unit_cauchy_ct, name_cauchy_ct )

    IF( ALLOCATED(shared_grid) ) DEALLOCATE( shared_grid )
    IF( ALLOCATED(convergence_factor) ) DEALLOCATE( convergence_factor )


  END PROCEDURE perform_cauchy_convergence_test_unknown_sol


  MODULE PROCEDURE perform_cauchy_convergence_test_known_sol

    !***********************************************************
    !
    !# Perform the Cauchy convergence test when the exact solution is
    !  known. The ratio between the grid spacings is `num/den`.
    !
    !***********************************************************

    IMPLICIT NONE

    INTEGER:: nx, ny, nz, unit_cauchy_ct

    DOUBLE PRECISION:: ratio_dx

    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: shared_grid
    DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE:: convergence_factor

    CHARACTER( LEN=: ), ALLOCATABLE:: name_cauchy_ct

    ratio_dx= num/den

    CALL find_shared_grid( tpo_coarse, tpo_fine, num, den, &
                           ref_lev, shared_grid )

    nx= SIZE(shared_grid(:,1,1,jx))
    ny= SIZE(shared_grid(1,:,1,jy))
    nz= SIZE(shared_grid(1,1,:,jz))

    PRINT *, "** Computing convergence factor..."

    choose_constraints: SELECT CASE( use_constraints )

    CASE(use_constraints_on_mesh)

      get_ham_ptr => get_ham

      unit_cauchy_ct= 3109
      name_cauchy_ct= TRIM(spacetime_path) &
                      //"cauchy_convergence_test_unknown.dat"

    CASE(use_constraints_with_mapped_hydro)

      get_ham_ptr => get_ham_parts

      unit_cauchy_ct= 3111
      name_cauchy_ct= TRIM(spacetime_path) &
                      //"cauchy_convergence_test_unknown_parts.dat"

    CASE DEFAULT

      PRINT *, "** There is no well defined algorithm " &
               // "corresponding to the number", use_constraints
      PRINT *, " * Please set use_constraints to 1 or 2."
      PRINT *, "   1: Use the constraints computed entirely on the mesh"
      PRINT *, "   2: Use the constraints computed on the mesh, using ", &
               "the hydro data mapped from the particles."
      PRINT *
      STOP

    END SELECT choose_constraints

    ALLOCATE( convergence_factor( nx, ny, nz ) )

    CALL compute_convergence_factor_known_sol( nx, ny, nz, num, den, &
                                               ref_lev, tpo_coarse, tpo_fine, &
                                               get_ham_ptr, convergence_factor )

    !
    !-- Print results to formatted file
    !
    CALL print_convergence_factor( nx, ny, nz, shared_grid, convergence_factor,&
                                   unit_cauchy_ct, name_cauchy_ct )

    IF( ALLOCATED(shared_grid) ) DEALLOCATE( shared_grid )
    IF( ALLOCATED(convergence_factor) ) DEALLOCATE( convergence_factor )

  END PROCEDURE perform_cauchy_convergence_test_known_sol


  SUBROUTINE print_convergence_factor &
    ( nx, ny, nz, shared_grid, convergence_factor, unit, filename )

    !***********************************************************
    !
    !# Print the Cauchy convergence factor to a formatted file
    !
    !***********************************************************

    IMPLICIT NONE

    INTEGER, INTENT(IN):: nx, ny, nz
    DOUBLE PRECISION, DIMENSION(nx,ny,nz,3), INTENT(IN):: shared_grid
    DOUBLE PRECISION, DIMENSION(nx,ny,nz),   INTENT(IN):: convergence_factor
    INTEGER, INTENT(IN):: unit
    LOGICAL:: exist
    CHARACTER( LEN=: ), ALLOCATABLE:: filename

    INTEGER:: i, j, k
    DOUBLE PRECISION:: min_abs_y, min_abs_z

    INQUIRE( FILE= TRIM(filename), EXIST= exist )

    IF( exist )THEN
      OPEN( UNIT= unit, FILE= TRIM(filename), &
            STATUS= "REPLACE", FORM= "FORMATTED", &
            POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ELSE
      OPEN( UNIT= unit, FILE= TRIM(filename), &
            STATUS= "NEW", FORM= "FORMATTED", &
            ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening ", TRIM(filename), &
               ". The error message is", err_msg
      STOP
    ENDIF

    WRITE( UNIT = unit, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id
    WRITE( UNIT= unit, IOSTAT = ios, &
           IOMSG = err_msg, FMT = * ) &
    "# Cauchy convergence test. "
    WRITE( UNIT = unit, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# column:      1        2       3       4"
    WRITE( UNIT = unit, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "#      x [km]       y [km]       z [km]       " &
    //"convergence factor [pure number]"

    min_abs_y= HUGE(one)
    min_abs_z= HUGE(one)
    DO j= 1, ny, 1
      IF( ABS( shared_grid(1,j,1,jy) ) < ABS( min_abs_y ) )THEN
        min_abs_y= shared_grid(1,j,1,jy)
      ENDIF
    ENDDO
    DO k= 1, nz, 1
      IF( ABS( shared_grid(1,1,k,jz) ) < ABS( min_abs_z ) )THEN
        min_abs_z= shared_grid(1,1,k,jz)
      ENDIF
    ENDDO

    DO k= 1, nz, 1
      DO j= 1, ny, 1
        DO i= 1, nx, 1

          IF( .FALSE. .AND. export_constraints_xy &
              .AND. &
              ABS(shared_grid(i,j,k,jz) - min_abs_z)/ABS(min_abs_z) > tol &
          )THEN

            CYCLE

          ENDIF
          IF( .FALSE. .AND. export_constraints_x &
              .AND. &
              ( ABS(shared_grid(i,j,k,jz) - min_abs_z)/ABS(min_abs_z) > tol &
              .OR. &
              ABS(shared_grid(i,j,k,jy) - min_abs_y)/ABS(min_abs_y) > tol ) &
          )THEN

            CYCLE

          ENDIF

          WRITE( UNIT = unit, IOSTAT = ios, &
                 IOMSG = err_msg, FMT = * )&
              shared_grid( i, j, k, jx ), &
              shared_grid( i, j, k, jy ), &
              shared_grid( i, j, k, jz ), &
              convergence_factor( i, j, k )

          IF( ios > 0 )THEN
            PRINT *, "...error when writing the arrays in ", &
                     TRIM(filename), &
                     ". The error message is", err_msg
            STOP
          ENDIF

        ENDDO
      ENDDO
    ENDDO

    CLOSE( UNIT= unit )

    PRINT *, " * Convergence factor printed to formatted file ", &
             TRIM(filename)
    PRINT *

  END SUBROUTINE print_convergence_factor


END SUBMODULE perform_test
