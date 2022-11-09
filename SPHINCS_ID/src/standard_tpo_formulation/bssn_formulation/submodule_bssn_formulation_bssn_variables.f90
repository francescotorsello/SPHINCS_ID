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
  !# Implementation of the methods of TYPE [[bssn]]
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
    !  @note Tested with 7 301**3 refinement levels
    !        on 05.07.2022. No memory issues.
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
    CALL deallocate_grid_function( lapse,      'lapse_tmp' )
    CALL deallocate_grid_function( shift_u,    'shift_u_tmp' )
    CALL deallocate_grid_function( dt_lapse,   'dt_lapse_tmp' )
    CALL deallocate_grid_function( dt_shift_u, 'dt_shift_u_tmp' )
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

      !************************************************
      !
      !# Search the |bssn| fields for NaNs, and STOP
      !  if one NaN is found, printing an informative
      !  message to the standard output.
      !
      !  Created:      FT 05.07.2022
      !  Last updated: FT 05.07.2022
      !
      !************************************************

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
              IF(.NOT.is_finite_number(g_BSSN3_ll% levels(l)% var(i,j,k,jxx)))&
              THEN
                PRINT *, "** ERROR! g_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * g_BSSN3_ll(l,i,j,k)=", &
                         g_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF(.NOT.is_finite_number(g_BSSN3_ll% levels(l)% var(i,j,k,jxy)))&
              THEN
                PRINT *, "** ERROR! g_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * g_BSSN3_ll(l,i,j,k)=", &
                         g_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF(.NOT.is_finite_number(g_BSSN3_ll% levels(l)% var(i,j,k,jxz)))&
              THEN
                PRINT *, "** ERROR! g_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * g_BSSN3_ll(l,i,j,k)=", &
                         g_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF(.NOT.is_finite_number(g_BSSN3_ll% levels(l)% var(i,j,k,jyy)))&
              THEN
                PRINT *, "** ERROR! g_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * g_BSSN3_ll(l,i,j,k)=", &
                         g_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF(.NOT.is_finite_number(g_BSSN3_ll% levels(l)% var(i,j,k,jyz)))&
              THEN
                PRINT *, "** ERROR! g_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * g_BSSN3_ll(l,i,j,k)=", &
                         g_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF(.NOT.is_finite_number(g_BSSN3_ll% levels(l)% var(i,j,k,jzz)))&
              THEN
                PRINT *, "** ERROR! g_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * g_BSSN3_ll(l,i,j,k)=", &
                         g_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF(.NOT.is_finite_number(A_BSSN3_ll% levels(l)% var(i,j,k,jxx)))&
              THEN
                PRINT *, "** ERROR! A_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * A_BSSN3_ll(l,i,j,k)=", &
                         A_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF(.NOT.is_finite_number(A_BSSN3_ll% levels(l)% var(i,j,k,jxy)))&
              THEN
                PRINT *, "** ERROR! A_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * A_BSSN3_ll(l,i,j,k)=", &
                         A_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF(.NOT.is_finite_number(A_BSSN3_ll% levels(l)% var(i,j,k,jxz)))&
              THEN
                PRINT *, "** ERROR! A_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * A_BSSN3_ll(l,i,j,k)=", &
                         A_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF(.NOT.is_finite_number(A_BSSN3_ll% levels(l)% var(i,j,k,jyy)))&
              THEN
                PRINT *, "** ERROR! A_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * A_BSSN3_ll(l,i,j,k)=", &
                         A_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF(.NOT.is_finite_number(A_BSSN3_ll% levels(l)% var(i,j,k,jyz)))&
              THEN
                PRINT *, "** ERROR! A_BSSN3_ll is not a finite number at ", &
                         "   (l,i,j,k)= (", l, ", ", i, ", ", j, ", ", k, ")"
                PRINT *, " * A_BSSN3_ll(l,i,j,k)=", &
                         A_BSSN3_ll% levels(l)% var(i,j,k,:)
                PRINT *
                STOP
              ENDIF
              IF(.NOT.is_finite_number(A_BSSN3_ll% levels(l)% var(i,j,k,jzz)))&
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

    USE utility,          ONLY: zero
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
                                 g_BSSN3_ll% levels(l)% var(:,:,:,jxx), &
                                 g_BSSN3_ll% levels(l)% var(:,:,:,jxy), &
                                 g_BSSN3_ll% levels(l)% var(:,:,:,jxz), &
                                 g_BSSN3_ll% levels(l)% var(:,:,:,jyy), &
                                 g_BSSN3_ll% levels(l)% var(:,:,:,jyz), &
                                 g_BSSN3_ll% levels(l)% var(:,:,:,jzz), &
                                 A_BSSN3_ll% levels(l)% var(:,:,:,jxx), &
                                 A_BSSN3_ll% levels(l)% var(:,:,:,jxy), &
                                 A_BSSN3_ll% levels(l)% var(:,:,:,jxz), &
                                 A_BSSN3_ll% levels(l)% var(:,:,:,jyy), &
                                 A_BSSN3_ll% levels(l)% var(:,:,:,jyz), &
                                 A_BSSN3_ll% levels(l)% var(:,:,:,jzz), &
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

   !!$OMP PARALLEL DO DEFAULT( NONE ) &
   !!$OMP             SHARED( Gamma_u, l ) &
   !!$OMP             PRIVATE( i, j, k )
   !DO k= 1, nz, 1
   !  DO j= 1, ny, 1
   !    DO i= 1, nx, 1
   !
   !      Gamma_u% levels(l)% var(i,j,k,jx)= zero
   !      Gamma_u% levels(l)% var(i,j,k,jy)= zero
   !      Gamma_u% levels(l)% var(i,j,k,jz)= zero
   !
   !    ENDDO
   !  ENDDO
   !ENDDO
   !!$OMP END PARALLEL DO
    Gamma_u% levels(l)% var= zero

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


END SUBMODULE bssn_variables
