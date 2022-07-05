! File:         submodule_bssn_formulation_bssn_variables.f90
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

SUBMODULE (bssn_formulation) bssn_variables

  !************************************************
  !
  !# Implementation of the methods of TYPE bssn
  !  that compute the |bssn| variables
  !
  !  FT 23.10.2020
  !
  !  Updated to support mesh refinement
  !
  !  FT 26.03.2021
  !
  !************************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE compute_and_print_bssn_variables

    !************************************************
    !
    !# Compute, stores and prints the BSSN variables
    !  to a binary file to be read by the evolution
    !  code SPHINCS_BSSN
    !
    !  Created:      FT 23.10.2020
    !  Last updated: FT 05.07.2022
    !
    !************************************************

    !USE NaNChecker,          ONLY: Check_Grid_Function_for_NAN
    USE tensor,                     ONLY: jx, jy, jz, &
                                          jxx, jxy, jxz, jyy, jyz, jzz
    USE utility,                    ONLY: zero, is_finite_number
    USE mesh_refinement,            ONLY: nlevels, levels, rad_coord, &
                                          allocate_grid_function, &
                                          deallocate_grid_function
    USE ADM_refine,                 ONLY: lapse, shift_u!, &
                                          !dt_lapse, dt_shift_u, &
                                          !K_phys3_ll, g_phys3_ll, &
                                          !allocate_ADM, deallocate_ADM
    USE McLachlan_refine,           ONLY: allocate_Ztmp, deallocate_Ztmp, &
                                          ADM_to_BSSN
    USE Tmunu_refine,               ONLY: allocate_Tmunu, deallocate_Tmunu, &
                                          Tmunu_ll
    USE GravityAcceleration_refine, ONLY: allocate_GravityAcceleration, &
                                          deallocate_GravityAcceleration
    !
    !-- Use the arrays from the MODULE BSSN to store the BSSN variables
    !-- for the LORENE ID on the grid, and the SUBROUTINE write_BSSN_dump
    !-- to export them to the binary file needed by the evolution code
    !-- in SPHINCS
    !
    USE BSSN_refine, ONLY: allocate_BSSN, deallocate_BSSN, &
                           Gamma_u,          & ! Conformal connection
                           phi,              & ! Conformal factor
                           trK,              & ! Trace of extrinsic curvature
                           A_BSSN3_ll,       & ! Conformal traceless
                                               ! extrinsic curvature
                           g_BSSN3_ll,       & ! Conformal metric
                           !Theta_Z4,         & ! Vector in the CCZ4 formulation
                                               ! Used because ADM_TO_BSSN
                                               ! calls SUBROUTINES that need it
                                               ! as input; however, it is not
                                               ! evolved in BSSN
                           !lapse_A_BSSN,     & ! Time derivative of lapse
                           !shift_B_BSSN_u,   & ! Time derivativeof shift
                           write_BSSN_dump

    IMPLICIT NONE

    ! The flag call_flag is set different than 0 if the SUBROUTINE
    ! compute_and_print_SPH_variables is called
    INTEGER, SAVE:: call_flag= 0

    INTEGER:: i, j, k, l

    TYPE(grid_function_scalar):: dt_lapse
    TYPE(grid_function)       :: dt_shift_u

    PRINT *, "** Computing and exporting BSSN ID..."

    ! Allocate memory for the ADM MODULE variables (this has to be done since
    ! the MODULE SUBROUTINES need them; not allocating it results in a
    ! segmentation fault)
    PRINT *
    PRINT *, " * Allocating needed memory..."
    PRINT *

    ALLOCATE ( levels( this% nlevels ), STAT=ios )
    IF( ios > 0 )THEN
     PRINT*,'...allocation error for levels'
     STOP
    ENDIF
    levels = this% levels
    nlevels= this% nlevels

    !DO l= 1, this% nlevels, 1
    !  levels(l)% ngrid_x= this% levels(l)% ngrid_x
    !  levels(l)% ngrid_x= this% levels(l)% ngrid_x
    !  levels(l)% ngrid_x= this% levels(l)% ngrid_x
    !ENDDO

    !CALL allocate_ADM()
    CALL allocate_BSSN()

    ! Allocate temporary memory for time integration
    CALL allocate_Ztmp()

    ! Allocate memory for the stress-energy tensor (used in write_BSSN_dump)
    CALL allocate_Tmunu()

    ! Allocate memory for the derivatives of the ADM variables
    CALL allocate_GravityAcceleration()

    !CALL allocate_grid_function( rad_coord, 'rad_coord' )
    CALL allocate_grid_function( lapse,      'lapse_tmp' )
    CALL allocate_grid_function( shift_u,    'shift_u_tmp', 3 )
    CALL allocate_grid_function( dt_lapse,   'dt_lapse_tmp' )
    CALL allocate_grid_function( dt_shift_u, 'dt_shift_u_tmp', 3 )

    ! Assign values to the MODULE variables, in order to call ADM_to_BSSN
    ref_levels: DO l= 1, this% nlevels

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( this, Tmunu_ll, dt_lapse, dt_shift_u, &
      !$OMP                     rad_coord, l ) &
      !$OMP             PRIVATE( i, j, k )
      DO k= 1, this% levels(l)% ngrid_z, 1
        DO j= 1, this% levels(l)% ngrid_y, 1
          DO i= 1, this% levels(l)% ngrid_x, 1

            Tmunu_ll% levels(l)%   var(i,j,k,:)= zero

            dt_lapse% levels(l)%   var(i,j,k)  = zero
            dt_shift_u% levels(l)% var(i,j,k,:)= zero

            !rad_coord% levels(l)%  var(i,j,k)  = &
            !  this% rad_coord% levels(l)% var(i,j,k)
            !lapse% levels(l)%      var(i,j,k)  = &
            !  this% lapse% levels(l)% var(i,j,k)
            !shift_u% levels(l)%    var(i,j,k,:)= &
            !  this% shift_u% levels(l)% var(i,j,k,:)
            !g_phys3_ll% levels(l)% var(i,j,k,:)= &
            !  this% g_phys3_ll% levels(l)% var(i,j,k,:)
            !K_phys3_ll% levels(l)% var(i,j,k,:)= &
            !  this% K_phys3_ll% levels(l)% var(i,j,k,:)

          ENDDO
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

    ENDDO ref_levels

    !
    !-- Compute BSSN variables, and time the process
    !-- The BSSN variables are stored in the MODULE variables since
    !-- write_BSSN_dump need them
    !
    PRINT *, " * Computing BSSN variables..."
    PRINT *
    CALL this% bssn_computer_timer% start_timer()
    !CALL ADM_to_BSSN()
    ref_levels2: DO l= 1, this% nlevels, 1

      CALL standard_tpo_to_bssn( l, &
        this% levels(l)% ngrid_x, this% levels(l)% ngrid_y, &
        this% levels(l)% ngrid_z, &
        this% levels(l)% dx, this% levels(l)% dy, this% levels(l)% dz, &
        this% levels(l)% nghost_x, this% levels(l)% nghost_y, &
        this% levels(l)% nghost_z, &
        ! Standard 3+1 variables (input)
        this% g_phys3_ll% levels(l)% var(:,:,:,jxx), &
        this% g_phys3_ll% levels(l)% var(:,:,:,jxy), &
        this% g_phys3_ll% levels(l)% var(:,:,:,jxz), &
        this% g_phys3_ll% levels(l)% var(:,:,:,jyy), &
        this% g_phys3_ll% levels(l)% var(:,:,:,jyz), &
        this% g_phys3_ll% levels(l)% var(:,:,:,jzz), &
        this% K_phys3_ll% levels(l)% var(:,:,:,jxx), &
        this% K_phys3_ll% levels(l)% var(:,:,:,jxy), &
        this% K_phys3_ll% levels(l)% var(:,:,:,jxz), &
        this% K_phys3_ll% levels(l)% var(:,:,:,jyy), &
        this% K_phys3_ll% levels(l)% var(:,:,:,jyz), &
        this% K_phys3_ll% levels(l)% var(:,:,:,jzz), &
        this% lapse     % levels(l)% var(:,:,:), &
        this% shift_u   % levels(l)% var(:,:,:,jx), &
        this% shift_u   % levels(l)% var(:,:,:,jy), &
        this% shift_u   % levels(l)% var(:,:,:,jz), &
        dt_lapse        % levels(l)% var, &
        dt_shift_u      % levels(l)% var(:,:,:,jx), &
        dt_shift_u      % levels(l)% var(:,:,:,jy), &
        dt_shift_u      % levels(l)% var(:,:,:,jz), &
        this% rad_coord % levels(l)% var &
      )

    ENDDO ref_levels2
    CALL this% bssn_computer_timer% stop_timer()

    !CALL deallocate_ADM()
    CALL deallocate_Ztmp()
    CALL deallocate_Tmunu()
    CALL deallocate_GravityAcceleration()

    !
    !-- Check the BSSN MODULE variables for NaNs
    !
    CALL check_bssn_id_for_NaNs()

    !
    !-- Setting the TYPE variables equal to the MODULE variables
    !
    CALL allocate_bssn_fields( this )
    ref_levels4: DO l= 1, this% nlevels

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( this, Gamma_u, phi, trK, g_BSSN3_ll, &
      !$OMP                     A_BSSN3_ll, lapse, shift_u, l ) &
      !$OMP             PRIVATE( i, j, k )
      DO k= 1, this% levels(l)% ngrid_z, 1
        DO j= 1, this% levels(l)% ngrid_y, 1
          DO i= 1, this% levels(l)% ngrid_x, 1

            this% Gamma_u%    levels(l)% var(i,j,k,:)= &
                  Gamma_u%    levels(l)% var(i,j,k,:)
            this% phi%        levels(l)% var(i,j,k)  = &
                  phi%        levels(l)% var(i,j,k)
            this% trK%        levels(l)% var(i,j,k)  = &
                  trK%        levels(l)% var(i,j,k)
            this% g_BSSN3_ll% levels(l)% var(i,j,k,:)= &
                  g_BSSN3_ll% levels(l)% var(i,j,k,:)
            this% A_BSSN3_ll% levels(l)% var(i,j,k,:)= &
                  A_BSSN3_ll% levels(l)% var(i,j,k,:)

                  lapse%      levels(l)% var(i,j,k)  = &
            this% lapse%      levels(l)% var(i,j,k)
                  shift_u%    levels(l)% var(i,j,k,:)= &
            this% shift_u%    levels(l)% var(i,j,k,:)

          ENDDO
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

    ENDDO ref_levels4

    ! Write BSSN ID to a binary file to be read by the evolution code
    ! in SPHINCS
    IF( this% export_bin )THEN
      IF( PRESENT(namefile) )THEN
        CALL write_BSSN_dump( namefile )
        !CALL write_BSSN_dump()
      ELSE
        CALL write_BSSN_dump()
      ENDIF
    ENDIF

    !
    !-- Deallocate MODULE variables
    !
    !CALL deallocate_ADM()
    !CALL deallocate_Ztmp()
    !CALL deallocate_Tmunu()
    !CALL deallocate_GravityAcceleration()
    CALL deallocate_BSSN()
    !CALL deallocate_grid_function( rad_coord, 'rad_coord' )
    !CALL deallocate_gravity_grid()
    DEALLOCATE( levels )

    call_flag= call_flag + 1
    this% call_flag= call_flag

    PRINT *, "** BSSN ID computed."
    PRINT *


    CONTAINS


    SUBROUTINE check_bssn_id_for_NaNs

      IMPLICIT NONE

      ref_levels: DO l= 1, this% nlevels

        !$OMP PARALLEL DO DEFAULT( NONE ) &
        !$OMP             SHARED( this, Gamma_u, phi, trK, g_BSSN3_ll, &
        !$OMP                     A_BSSN3_ll, l ) &
        !$OMP             PRIVATE( i, j, k )
        DO k= 1, this% levels(l)% ngrid_z, 1
          DO j= 1, this% levels(l)% ngrid_y, 1
            DO i= 1, this% levels(l)% ngrid_x, 1

              IF( .NOT.is_finite_number(Gamma_u% levels(l)% var(i,j,k,jx)) )THEN
                PRINT *, "** ERROR! Gamma_u is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * Gamma_u(l,i,j,k)=", &
                         Gamma_u% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF( .NOT.is_finite_number(Gamma_u% levels(l)% var(i,j,k,jy)) )THEN
                PRINT *, "** ERROR! Gamma_u is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * Gamma_u(l,i,j,k)=", &
                         Gamma_u% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF( .NOT.is_finite_number(Gamma_u% levels(l)% var(i,j,k,jz)) )THEN
                PRINT *, "** ERROR! Gamma_u is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * Gamma_u(l,i,j,k)=", &
                         Gamma_u% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF( .NOT.is_finite_number(phi% levels(l)% var(i,j,k)) )THEN
                PRINT *, "** ERROR! phi is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * phi(l,i,j,k)=", phi% levels(l)% var(i,j,k)
                PRINT *
                STOP
              ENDIF
              IF( .NOT.is_finite_number(trK% levels(l)% var(i,j,k)) )THEN
                PRINT *, "** ERROR! trK is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * trK(l,i,j,k)=", trK% levels(l)% var(i,j,k)
                PRINT *
                STOP
              ENDIF
              IF( .NOT.is_finite_number(g_BSSN3_ll% levels(l)% var(i,j,k,jxx)) )&
              THEN
                PRINT *, "** ERROR! g_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * g_BSSN3_ll(l,i,j,k)=", &
                         g_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF( .NOT.is_finite_number(g_BSSN3_ll% levels(l)% var(i,j,k,jxy)) )&
              THEN
                PRINT *, "** ERROR! g_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * g_BSSN3_ll(l,i,j,k)=", &
                         g_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF( .NOT.is_finite_number(g_BSSN3_ll% levels(l)% var(i,j,k,jxz)) )&
              THEN
                PRINT *, "** ERROR! g_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * g_BSSN3_ll(l,i,j,k)=", &
                         g_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF( .NOT.is_finite_number(g_BSSN3_ll% levels(l)% var(i,j,k,jyy)) )&
              THEN
                PRINT *, "** ERROR! g_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * g_BSSN3_ll(l,i,j,k)=", &
                         g_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF( .NOT.is_finite_number(g_BSSN3_ll% levels(l)% var(i,j,k,jyz)) )&
              THEN
                PRINT *, "** ERROR! g_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * g_BSSN3_ll(l,i,j,k)=", &
                         g_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF( .NOT.is_finite_number(g_BSSN3_ll% levels(l)% var(i,j,k,jzz)) )&
              THEN
                PRINT *, "** ERROR! g_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * g_BSSN3_ll(l,i,j,k)=", &
                         g_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF( .NOT.is_finite_number(A_BSSN3_ll% levels(l)% var(i,j,k,jxx)) )&
              THEN
                PRINT *, "** ERROR! A_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * A_BSSN3_ll(l,i,j,k)=", &
                         A_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF( .NOT.is_finite_number(A_BSSN3_ll% levels(l)% var(i,j,k,jxy)) )&
              THEN
                PRINT *, "** ERROR! A_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * A_BSSN3_ll(l,i,j,k)=", &
                         A_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF( .NOT.is_finite_number(A_BSSN3_ll% levels(l)% var(i,j,k,jxz)) )&
              THEN
                PRINT *, "** ERROR! A_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * A_BSSN3_ll(l,i,j,k)=", &
                         A_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF( .NOT.is_finite_number(A_BSSN3_ll% levels(l)% var(i,j,k,jyy)) )&
              THEN
                PRINT *, "** ERROR! A_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * A_BSSN3_ll(l,i,j,k)=", &
                         A_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF( .NOT.is_finite_number(A_BSSN3_ll% levels(l)% var(i,j,k,jyz)) )&
              THEN
                PRINT *, "** ERROR! A_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * A_BSSN3_ll(l,i,j,k)=", &
                         A_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF( .NOT.is_finite_number(A_BSSN3_ll% levels(l)% var(i,j,k,jzz)) )&
              THEN
                PRINT *, "** ERROR! A_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * A_BSSN3_ll(l,i,j,k)=", &
                         A_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF

            ENDDO
          ENDDO
        ENDDO
        !$OMP END PARALLEL DO

      ENDDO ref_levels

    END SUBROUTINE check_bssn_id_for_NaNs


  END PROCEDURE compute_and_print_bssn_variables


  SUBROUTINE standard_tpo_to_bssn( l, nx, ny, nz, dx, dy, dz, ngx, ngy, ngz, &
    gxx, gxy, gxz, gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz, kzz, alp, betax, &
    betay, betaz, dtalp, dtbetax, dtbetay, dtbetaz, r )

    !************************************************
    !
    !# Compute the BSSN variables starting from the
    !  standard 3+1 (aka ADM) variables
    !  This is basically a version of ADM_to_BSSN
    !  from MODULE McLachlan_refine that allows for
    !  array arguments.
    !
    !  FT 05.07.2022
    !  FT 05.07.2022
    !
    !************************************************

    USE tensor,           ONLY: jx, jy, jz, jxx, jxy, jxz, jyy, jyz, jzz
    USE BSSN_refine,      ONLY: g_BSSN3_ll, A_BSSN3_ll, lapse_A_BSSN, &
                                shift_B_BSSN_u, Theta_Z4, Gamma_u, phi, trK

    USE McLachlan_refine, ONLY: ADM_TO_BSSN_EVERYWHERE, ADM_TO_BSSN_INTERIOR

    IMPLICIT NONE

    !TYPE(bssn), INTENT(INOUT):: bssnf
    INTEGER,          INTENT(IN):: l, nx ,ny, nz, ngx, ngy, ngz
    DOUBLE PRECISION, INTENT(IN):: dx, dy, dz
    DOUBLE PRECISION, INTENT(IN):: gxx(nx,ny,nz), gxy(nx,ny,nz), &
                                   gxz(nx,ny,nz), gyy(nx,ny,nz), &
                                   gyz(nx,ny,nz), gzz(nx,ny,nz)
    DOUBLE PRECISION, INTENT(IN):: kxx(nx,ny,nz), kxy(nx,ny,nz), &
                                   kxz(nx,ny,nz), kyy(nx,ny,nz), &
                                   kyz(nx,ny,nz), kzz(nx,ny,nz)
    DOUBLE PRECISION, INTENT(IN):: alp(nx,ny,nz)
    DOUBLE PRECISION, INTENT(IN):: betax(nx,ny,nz), betay(nx,ny,nz), &
                                   betaz(nx,ny,nz)
    DOUBLE PRECISION, INTENT(IN):: dtalp(nx,ny,nz)
    DOUBLE PRECISION, INTENT(IN):: dtbetax(nx,ny,nz), dtbetay(nx,ny,nz), &
                                   dtbetaz(nx,ny,nz)
    DOUBLE PRECISION, INTENT(IN):: r(nx,ny,nz)

    !----------------------------------------------------------------!
    !-- It is assumed that the outer boundary width is the same as --!
    !-- the ghostsize                                              --!
    !----------------------------------------------------------------!
    INTEGER, DIMENSION(3) :: imin, imax
    DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: tmp1, tmp2, tmp3, tmp4

    !--------------------------------------------------------------!
    !-- These are to be used in C, hence the zero based indexing --!
    !--------------------------------------------------------------!
    imin= 0
    imax(1)= nx - 1
    imax(2)= ny - 1
    imax(3)= nz - 1

    ALLOCATE( tmp1(nx,ny,nz), tmp2(nx,ny,nz), tmp3(nx,ny,nz), tmp4(nx,ny,nz) )

    CALL ADM_TO_BSSN_EVERYWHERE( nx, ny, nz, imin, imax, dx, dy, dz, &
                                 gxx, gxy, gxz, gyy, gyz, gzz, &
                                 kxx, kxy, kxz, kyy, kyz, kzz, &
                                 alp, betax, betay, betaz, &
                                 g_BSSN3_ll%levels(l)%var(:,:,:,jxx), &
                                 g_BSSN3_ll%levels(l)%var(:,:,:,jxy), &
                                 g_BSSN3_ll%levels(l)%var(:,:,:,jxz), &
                                 g_BSSN3_ll%levels(l)%var(:,:,:,jyy), &
                                 g_BSSN3_ll%levels(l)%var(:,:,:,jyz), &
                                 g_BSSN3_ll%levels(l)%var(:,:,:,jzz), &
                                 A_BSSN3_ll%levels(l)%var(:,:,:,jxx), &
                                 A_BSSN3_ll%levels(l)%var(:,:,:,jxy), &
                                 A_BSSN3_ll%levels(l)%var(:,:,:,jxz), &
                                 A_BSSN3_ll%levels(l)%var(:,:,:,jyy), &
                                 A_BSSN3_ll%levels(l)%var(:,:,:,jyz), &
                                 A_BSSN3_ll%levels(l)%var(:,:,:,jzz), &
                                 tmp1(:,:,:), tmp2(:,:,:), &
                                 tmp3(:,:,:), tmp4(:,:,:), &
                                 phi%levels(l)%var(:,:,:), &
                                 trK%levels(l)%var(:,:,:), &
                                 Theta_Z4%levels(l)%var(:,:,:) )

    DEALLOCATE( tmp1, tmp2, tmp3, tmp4 )

    imin = [ ngx, ngy, ngz ]
    imax(1) = nx - ngx - 1
    imax(2) = ny - ngy - 1
    imax(3) = nz - ngz - 1

    CALL ADM_TO_BSSN_INTERIOR( nx, ny, nz, imin, imax, dx, dy, dz, &
                               gxx, gxy, gxz, gyy, gyz, gzz, &
                               kxx, kxy, kxz, kyy, kyz, kzz, &
                               alp, betax, betay, betaz, &
                               dtalp, dtbetax, dtbetay, dtbetaz, &
                               g_BSSN3_ll% levels(l)% var(:,:,:,jxx), &
                               g_BSSN3_ll% levels(l)% var(:,:,:,jxy), &
                               g_BSSN3_ll% levels(l)% var(:,:,:,jxz), &
                               g_BSSN3_ll% levels(l)% var(:,:,:,jyy), &
                               g_BSSN3_ll% levels(l)% var(:,:,:,jyz), &
                               g_BSSN3_ll% levels(l)% var(:,:,:,jzz), &
                               phi% levels(l)% var(:,:,:), &
                               r, &
                               lapse_A_BSSN% levels(l)% var(:,:,:), &
                               shift_B_BSSN_u% levels(l)% var(:,:,:,jx), &
                               shift_B_BSSN_u% levels(l)% var(:,:,:,jy), &
                               shift_B_BSSN_u% levels(l)% var(:,:,:,jz), &
                               Gamma_u% levels(l)% var(:,:,:,jx), &
                               Gamma_u% levels(l)% var(:,:,:,jy), &
                               Gamma_u% levels(l)% var(:,:,:,jz) )

  END SUBROUTINE standard_tpo_to_bssn


  MODULE PROCEDURE read_bssn_dump_print_formatted

    !************************************************
    !
    !# Read the BSSN ID from the binary file output
    !  by write_BSSN_dump, and print it to a
    !  formatted file
    !
    !  FT 08.02.2021
    !
    !************************************************

    USE mesh_refinement,  ONLY: levels, nlevels
    USE tensor,           ONLY: jxx, jxy, jxz, jyy, jyz, jzz, jx, jy, jz
    USE ADM_refine,       ONLY: lapse, shift_u, &
                                allocate_ADM, deallocate_ADM
    USE BSSN_refine,  ONLY: allocate_BSSN, deallocate_BSSN, &
                            Gamma_u,        & ! Conformal connection
                            phi,            & ! Conformal factor
                            trK,            & ! Trace of extrinsic curvature
                            A_BSSN3_ll,     & ! Conformal traceless
                                              ! extrinsic curvature
                            g_BSSN3_ll,     & ! Conformal metric
                            !Theta_Z4,       & ! Vector in the CCZ4 formulation.
                                              ! Loaded here because ADM_TO_BSSN
                                              ! calls SUBROUTINES that need it
                                              ! as input; however, it is not
                                              ! evolved in BSSN
                            !lapse_A_BSSN,   & ! Time derivative of lapse
                            !shift_B_BSSN_u, & ! Time derivativeof shift
                            read_BSSN_dump

    IMPLICIT NONE

    INTEGER:: i, j, k, l, min_ix_y, min_iy_y, min_iz_y, &
              min_ix_z, min_iy_z, min_iz_z

    DOUBLE PRECISION:: min_abs_y, min_abs_z

    LOGICAL:: exist

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    PRINT *, "** Executing the read_bssn_dump_print_formatted subroutine..."

    levels = this% levels
    nlevels= this% nlevels

    CALL allocate_ADM()
    CALL allocate_BSSN()

    CALL read_BSSN_dump( 00000, namefile_bin )

    IF( this% call_flag == 0 )THEN
      PRINT *, "** The SUBROUTINE print_formatted_id_bssn_variables ", &
        " must be called after compute_and_print_bssn_variables, otherwise", &
        " there are no bssn fields to export to the formatted file."
      PRINT *, "   Aborting."
      PRINT *
      STOP
    ENDIF

    IF( PRESENT(namefile) )THEN
      finalnamefile= namefile
    ELSE
      finalnamefile= "bssn_vars.dat"
    ENDIF

    INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

    IF( exist )THEN
      OPEN( UNIT= 20, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
            FORM= "FORMATTED", &
            POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
            IOMSG= err_msg )
    ELSE
      OPEN( UNIT= 20, FILE= TRIM(finalnamefile), STATUS= "NEW", &
      FORM= "FORMATTED", &
            ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening ", TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when opening " &
    !         // TRIM(finalnamefile) )

    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id
    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Values of the fields (including coordinates) exported by LORENE "&
    // "on each grid point"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 1 in ", TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 1 in "&
    !         // TRIM(finalnamefile) )
    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# column:      1        2       3       4       5", &
    "       6       7       8", &
    "       9       10      11", &
    "       12      13      14", &
    "       15      16      17      18      19", &
    "       20      21      22", &
    "       23      24"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 2 in ", TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 2 in "&
    !        // TRIM(finalnamefile) )
    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "#      refinement level    x [km]       y [km]       z [km]       lapse", &
    "       shift_x [c]    shift_y [c]    shift_z [c]", &
    "       conformal factor phi        trace of extr. curv. trK", &
    "       g_BSSN_xx       g_BSSN_xy      g_BSSN_xz", &
    "       g_BSSN_yy       g_BSSN_yz      g_BSSN_zz", &
    "       A_BSSN_xx       A_BSSN_xy      A_BSSN_xz    ", &
    "       A_BSSN_yy       A_BSSN_yz      A_BSSN_zz", &
    "       Gamma_u_x       Gamma_u_y      Gamma_u_z"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 3 in ", TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 3 in "&
    !        // TRIM(finalnamefile) )

    DO l= 1, this% nlevels, 1

      ASSOCIATE( lapse      => lapse% levels(l)% var, &
                 shift_u    => shift_u% levels(l)% var, &
                 g_BSSN3_ll => g_BSSN3_ll% levels(l)% var, &
                 A_BSSN3_ll => A_BSSN3_ll% levels(l)% var, &
                 phi        => phi% levels(l)% var, &
                 trK        => trK% levels(l)% var, &
                 Gamma_u    => Gamma_u% levels(l)% var &
      )

        min_abs_y= 1D+20
        min_abs_z= 1D+20
        DO k= 1, this% get_ngrid_z(l), 1
          DO j= 1, this% get_ngrid_y(l), 1
            DO i= 1, this% get_ngrid_x(l), 1

              IF( ABS( this% coords% levels(l)% var( i, j, k, jy ) ) &
                  < min_abs_y )THEN
                min_abs_y= ABS( this% coords% levels(l)% var( i, j, k, jy ) )
                min_ix_y= i
                min_iy_y= j
                min_iz_y= k
              ENDIF

              IF( ABS( this% coords% levels(l)% var( i, j, k, jz ) ) &
                  < min_abs_z )THEN
                min_abs_z= ABS( this% coords% levels(l)% var( i, j, k, jz ) )
                min_ix_z= i
                min_iy_z= j
                min_iz_z= k
              ENDIF

            ENDDO
          ENDDO
        ENDDO

        DO k= 1, this% get_ngrid_z(l), 1

          DO j= 1, this% get_ngrid_y(l), 1

            DO i= 1, this% get_ngrid_x(l), 1

              IF( this% export_form_xy .AND. &
                  ( this% coords% levels(l)% var( i, j, k, jz ) /= &
                    this% coords% levels(l)% var( min_ix_z, min_iy_z, &
                                                  min_iz_z, jz ) ) )THEN
                CYCLE
              ENDIF
              IF( this% export_form_x .AND. &
                  ( this% coords% levels(l)% var( i, j, k, jz ) /= &
                    this% coords% levels(l)% var( min_ix_z, min_iy_z, &
                                                  min_iz_z, jz ) &
                    .OR. &
                    this% coords% levels(l)% var( i, j, k, jy ) /= &
                    this% coords% levels(l)% var( min_ix_y, min_iy_y, &
                                                  min_iz_y, jy ) ) )THEN
                CYCLE
              ENDIF

              WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * )&
                  l, &
                  this% coords% levels(l)% var( i, j, k, jx ), &
                  this% coords% levels(l)% var( i, j, k, jy ), &
                  this% coords% levels(l)% var( i, j, k, jz ), &
                  lapse( i, j, k ), &
                  shift_u( i, j, k, jx ), &
                  shift_u( i, j, k, jy ), &
                  shift_u( i, j, k, jz ), &
                  phi( i, j, k ), &
                  trK( i, j, k ), &
                  g_BSSN3_ll( i, j, k, jxx ), &
                  g_BSSN3_ll( i, j, k, jxy ), &
                  g_BSSN3_ll( i, j, k, jxz ), &
                  g_BSSN3_ll( i, j, k, jyy ), &
                  g_BSSN3_ll( i, j, k, jyz ), &
                  g_BSSN3_ll( i, j, k, jzz ), &
                  A_BSSN3_ll( i, j, k, jxx ), &
                  A_BSSN3_ll( i, j, k, jxy ), &
                  A_BSSN3_ll( i, j, k, jxz ), &
                  A_BSSN3_ll( i, j, k, jyy ), &
                  A_BSSN3_ll( i, j, k, jyz ), &
                  A_BSSN3_ll( i, j, k, jzz ), &
                  Gamma_u( i, j, k, jx ), &
                  Gamma_u( i, j, k, jy ), &
                  Gamma_u( i, j, k, jz )

              IF( ios > 0 )THEN
                PRINT *, "...error when writing the arrays in ", &
                         TRIM(finalnamefile), ". The error message is", err_msg
                STOP
              ENDIF
              !CALL test_status( ios, err_msg, "...error when writing " &
              !                  // "the arrays in " // TRIM(namefile) )

            ENDDO
          ENDDO
        ENDDO
      END ASSOCIATE
    ENDDO

    CLOSE( UNIT= 20 )

    !
    !-- Deallocate MODULE variables
    !
    CALL deallocate_ADM()
    CALL deallocate_BSSN()

    PRINT *, " * LORENE BSSN ID on the refined mesh, to be supplied to ", &
             "SPHINCS_BSSN, printed to formatted file ", TRIM(namefile)

    PRINT *, "** Subroutine read_bssn_dump_print_formatted " &
             // "executed."
    PRINT *

  END PROCEDURE read_bssn_dump_print_formatted


  MODULE PROCEDURE print_formatted_id_bssn_variables

    !************************************************
    !
    !# Print the BSSN ID, computed on the gravity
    !  grid, to a formatted file
    !
    !  FT 26.10.2020
    !
    !************************************************

    USE tensor,              ONLY: jxx, jxy, jxz, jyy, jyz, jzz, jx, jy, jz

    IMPLICIT NONE

    INTEGER:: i, j, k, l, min_ix_y, min_iy_y, min_iz_y, &
              min_ix_z, min_iy_z, min_iz_z

    DOUBLE PRECISION:: min_abs_y, min_abs_z

    LOGICAL:: exist

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    ! Being abs_grid a local array, it is good practice to allocate it on the
    ! heap, otherwise it will be stored on the stack which has a very limited
    ! size. This results in a segmentation fault.
    !ALLOCATE( abs_grid( 3, this% ngrid_x, this% ngrid_y, this% ngrid_z ) )

    PRINT *, "** Executing the print_formatted_id_BSSN_variables " &
             // "subroutine..."

    IF( this% call_flag == 0 )THEN
      PRINT *, "** The SUBROUTINE print_formatted_id_bssn_variables ", &
        " must be called after compute_and_print_bssn_variables, otherwise", &
        " there are no bssn fields to export to the formatted file."
      PRINT *, "   Aborting."
      PRINT *
      STOP
    ENDIF

    IF( PRESENT(namefile) )THEN
      finalnamefile= namefile
    ELSE
      finalnamefile= "lorene-bns-id-bssn-form.dat"
    ENDIF

    INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

    IF( exist )THEN
      OPEN( UNIT= 20, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
            FORM= "FORMATTED", &
            POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
            IOMSG= err_msg )
    ELSE
      OPEN( UNIT= 20, FILE= TRIM(finalnamefile), STATUS= "NEW", &
      FORM= "FORMATTED", &
            ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening ", TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when opening " &
    !         // TRIM(finalnamefile) )

    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Run ID [ccyymmdd-hhmmss.sss]: " // run_id
    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# Values of the fields (including coordinates) exported by LORENE "&
    // "on each grid point"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 1 in ", TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 1 in "&
    !         // TRIM(finalnamefile) )
    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "# column:      1        2       3       4       5", &
    "       6       7       8", &
    "       9       10      11", &
    "       12      13      14", &
    "       15      16      17      18      19", &
    "       20      21      22", &
    "       23      24    25"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 2 in ", TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 2 in "&
    !        // TRIM(finalnamefile) )
    WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "#      refinement level    x [km]       y [km]       z [km]       lapse", &
    "       shift_x [c]    shift_y [c]    shift_z [c]", &
    "       conformal factor phi        trace of extr. curv. trK", &
    "       g_BSSN_xx       g_BSSN_xy      g_BSSN_xz", &
    "       g_BSSN_yy       g_BSSN_yz      g_BSSN_zz", &
    "       A_BSSN_xx       A_BSSN_xy      A_BSSN_xz    ", &
    "       A_BSSN_yy       A_BSSN_yz      A_BSSN_zz", &
    "       Gamma_u_x       Gamma_u_y      Gamma_u_z"
    IF( ios > 0 )THEN
      PRINT *, "...error when writing line 3 in ", TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF
    !CALL test_status( ios, err_msg, "...error when writing line 3 in "&
    !        // TRIM(finalnamefile) )

    DO l= 1, this% nlevels, 1

      ASSOCIATE( lapse           => this% lapse% levels(l)% var, &
                 shift_u         => this% shift_u% levels(l)% var, &
                 phi             => this% phi% levels(l)% var, &
                 trK             => this% trK% levels(l)% var, &
                 g_BSSN3_ll      => this% g_BSSN3_ll% levels(l)% var, &
                 A_BSSN3_ll      => this% A_BSSN3_ll% levels(l)% var, &
                 Gamma_u         => this% Gamma_u% levels(l)% var &
      )


        !DO iz= 1, this% ngrid_z, 1
        !  DO iy= 1, this% ngrid_y, 1
        !    DO ix= 1, this% ngrid_x, 1
        !      abs_grid( 1, ix, iy, iz )= ABS( this% grid( 1, ix, iy, iz ) )
        !      abs_grid( 2, ix, iy, iz )= ABS( this% grid( 2, ix, iy, iz ) )
        !      abs_grid( 3, ix, iy, iz )= ABS( this% grid( 3, ix, iy, iz ) )
        !    ENDDO
        !  ENDDO
        !ENDDO

        min_abs_y= HUGE(1.0D0)
        min_abs_z= HUGE(1.0D0)
        DO k= 1, this% get_ngrid_z(l), 1
          DO j= 1, this% get_ngrid_y(l), 1
            DO i= 1, this% get_ngrid_x(l), 1

              IF( ABS( this% coords% levels(l)% var( i, j, k, jy ) ) &
                  < min_abs_y )THEN
                min_abs_y= ABS( this% coords% levels(l)% var( i, j, k, jy ) )
                min_ix_y= i
                min_iy_y= j
                min_iz_y= k
              ENDIF

              IF( ABS( this% coords% levels(l)% var( i, j, k, jz ) ) &
                  < min_abs_z )THEN
                min_abs_z= ABS( this% coords% levels(l)% var( i, j, k, jz ) )
                min_ix_z= i
                min_iy_z= j
                min_iz_z= k
              ENDIF

            ENDDO
          ENDDO
        ENDDO

        DO k= 1, this% get_ngrid_z(l), 1

          DO j= 1, this% get_ngrid_y(l), 1

            DO i= 1, this% get_ngrid_x(l), 1

              IF( this% export_form_xy .AND. &
                  ( this% coords% levels(l)% var( i, j, k, jz ) /= &
                    this% coords% levels(l)% var( min_ix_z, min_iy_z, &
                                                  min_iz_z, jz ) ) )THEN
                CYCLE
              ENDIF
              IF( this% export_form_x .AND. &
                  ( this% coords% levels(l)% var( i, j, k, jz ) /= &
                    this% coords% levels(l)% var( min_ix_z, min_iy_z, &
                                                  min_iz_z, jz ) &
                    .OR. &
                    this% coords% levels(l)% var( i, j, k, jy ) /= &
                    this% coords% levels(l)% var( min_ix_y, min_iy_y, &
                                                  min_iz_y, jy ) ) )THEN
                CYCLE
              ENDIF

              WRITE( UNIT = 20, IOSTAT = ios, IOMSG = err_msg, FMT = * )&
                l, &
                this% coords% levels(l)% var( i, j, k, jx ), &
                this% coords% levels(l)% var( i, j, k, jy ), &
                this% coords% levels(l)% var( i, j, k, jz ), &
                lapse( i, j, k ), &
                shift_u( i, j, k, jx ), &
                shift_u( i, j, k, jy ), &
                shift_u( i, j, k, jz ), &
                phi( i, j, k ), &
                trK( i, j, k ), &
                g_BSSN3_ll( i, j, k, jxx ), &
                g_BSSN3_ll( i, j, k, jxy ), &
                g_BSSN3_ll( i, j, k, jxz ), &
                g_BSSN3_ll( i, j, k, jyy ), &
                g_BSSN3_ll( i, j, k, jyz ), &
                g_BSSN3_ll( i, j, k, jzz ), &
                A_BSSN3_ll( i, j, k, jxx ), &
                A_BSSN3_ll( i, j, k, jxy ), &
                A_BSSN3_ll( i, j, k, jxz ), &
                A_BSSN3_ll( i, j, k, jyy ), &
                A_BSSN3_ll( i, j, k, jyz ), &
                A_BSSN3_ll( i, j, k, jzz ), &
                Gamma_u( i, j, k, jx ), &
                Gamma_u( i, j, k, jy ), &
                Gamma_u( i, j, k, jz )

              IF( ios > 0 )THEN
                PRINT *, "...error when writing the arrays in ", &
                         TRIM(finalnamefile), ". The error message is", err_msg
                STOP
              ENDIF
              !CALL test_status( ios, err_msg, "...error when writing " &
              !                  // "the arrays in " // TRIM(finalnamefile) )

            ENDDO
          ENDDO
        ENDDO
      END ASSOCIATE
    ENDDO

    CLOSE( UNIT= 20 )

    PRINT *, " * LORENE BSSN ID on the gravity grid saved to formatted " &
             // "file ", TRIM(finalnamefile)

    PRINT *, "** Subroutine print_formatted_id_BSSN_variables " &
             // "executed."
    PRINT *

  END PROCEDURE print_formatted_id_bssn_variables


END SUBMODULE bssn_variables
