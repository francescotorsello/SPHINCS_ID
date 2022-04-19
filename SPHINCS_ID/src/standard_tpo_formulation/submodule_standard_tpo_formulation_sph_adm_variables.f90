! File:         submodule_standard_tpo_formulation_sph_adm_variables.f90
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

SUBMODULE (standard_tpo_formulation) sph_adm_variables

  !************************************************
  !
  !# This SUBMODULE contains the implementation
  !  of the methods that compute an estimate of
  !  the ADM variables on the |sph| fluid, using
  !  the metric mapped from the mesh to the particles
  !
  !  FT 12.04.2020
  !
  !************************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE compute_adm_momentum_fluid_m2p

    !************************************************
    !
    !# Computes an estimate of the \(\mathrm{ADM}\)
    !  linear momentum using the |sph| fields
    !  on the particles
    !  @todo add reference
    !
    !  FT 12.04.2022
    !
    !************************************************

    USE tensor,               ONLY: jx, jy, jz, n_sym4x4, lower_index_4vector
    USE constants,            ONLY: zero, one, two, amu, MSun
    USE utility,              ONLY: compute_g4, determinant_sym4x4, &
                                    spatial_vector_norm_sym3x3, &
                                    compute_tpo_metric

    USE metric_on_particles,  ONLY: allocate_metric_on_particles, &
                                    deallocate_metric_on_particles, &
                                    get_metric_on_particles, g4_ll

    USE ADM_refine,           ONLY: lapse, shift_u, g_phys3_ll, &
                                    allocate_ADM, deallocate_ADM
    USE BSSN_refine,          ONLY: allocate_BSSN, deallocate_BSSN
    USE mesh_refinement,      ONLY: nlevels, levels, rad_coord, coords, &
                                    allocate_grid_function, &
                                    deallocate_grid_function
    USE gradient,             ONLY: allocate_gradient, deallocate_gradient
    USE McLachlan_refine,     ONLY: allocate_Ztmp, deallocate_Ztmp
    USE alive_flag,           ONLY: alive
    USE input_output,         ONLY: read_options
    USE units,                ONLY: set_units
    USE sph_variables,        ONLY: npart, allocate_SPH_memory, &
                                    deallocate_SPH_memory

    IMPLICIT NONE

    INTEGER, PARAMETER:: unit_recovery= 34956

    INTEGER:: a, j, l, i_matter

    DOUBLE PRECISION:: det, shift_norm2

    !DOUBLE PRECISION, DIMENSION(n_sym4x4,npart)  :: g4
    DOUBLE PRECISION, DIMENSION(0:3)             :: v_u
    DOUBLE PRECISION, DIMENSION(0:3,parts% get_npart()):: v_l

    DOUBLE PRECISION              :: lapse_loc
    DOUBLE PRECISION, DIMENSION(3):: shift
    DOUBLE PRECISION, DIMENSION(6):: g3

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: vel_loc
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nu_loc
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nlrf_loc
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: u_loc
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pr_loc
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: theta_loc

    DOUBLE PRECISION, DIMENSION(parts% get_n_matter(),3):: adm_mom_i

    LOGICAL, PARAMETER:: debug= .FALSE.

    npart    = parts% get_npart()
    pos      = parts% get_pos()
    vel_loc  = parts% get_vel()
    nu_loc   = parts% get_nu()
    nlrf_loc = parts% get_nlrf()
    u_loc    = parts% get_u()
    pr_loc   = parts% get_pressure_cu()
    theta_loc= parts% get_theta()

    ALLOCATE ( levels( this% nlevels ), STAT=ios )
    IF( ios > 0 )THEN
     PRINT*,'...allocation error for levels'
     STOP
    ENDIF
    nlevels= this% nlevels
    levels = this% levels
    coords = this% coords

    CALL allocate_ADM()
    CALL allocate_BSSN()

    ! Allocate temporary memory for time integration
    CALL allocate_Ztmp()

    ! Allocate memory for the derivatives of the ADM variables
    ! CALL allocate_GravityAcceleration()

    CALL allocate_grid_function( rad_coord, 'rad_coord', 1 )

    ! Initialize the stress-energy tensor to 0
    DO l= 1, this% nlevels, 1
      rad_coord%  levels(l)% var= this% rad_coord%  levels(l)% var
      g_phys3_ll% levels(l)% var= this% g_phys3_ll% levels(l)% var
      shift_u%    levels(l)% var= this% shift_u%    levels(l)% var
      lapse%      levels(l)% var= this% lapse%      levels(l)% var
    ENDDO

    CALL set_units('NSM')
    CALL read_options

    CALL allocate_SPH_memory

    ! flag that particles are 'alive'
    IF( .NOT.ALLOCATED( alive ) ) ALLOCATE( alive( npart ) )
    alive( 1:npart )= 1

    CALL allocate_gradient( npart )

    IF( ALLOCATED(g4_ll) )THEN
      DEALLOCATE(g4_ll)
    ENDIF
    CALL allocate_metric_on_particles(npart)

    PRINT *, " * Mapping metric from the grid to the particles..."
    PRINT *
    CALL get_metric_on_particles( npart, pos )

    adm_mom= zero
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( npart, parts, g4_ll, vel_loc, v_l, nu_loc, &
    !$OMP                     theta_loc, pr_loc, nlrf_loc, u_loc ) &
    !$OMP             PRIVATE( a, det, v_u, shift_norm2, j, g3, lapse_loc, &
    !$OMP                      shift ) &
    !$OMP             REDUCTION( +: adm_mom )
    DO a= 1, npart, 1

      CALL compute_tpo_metric( g4_ll(:,a), lapse_loc, shift, g3 )

      v_u= [one, vel_loc(1:3,a)]
      CALL lower_index_4vector( v_u, g4_ll(:,a), v_l(:,a) )

      CALL spatial_vector_norm_sym3x3( g3, shift, shift_norm2 )

      DO j= jx, jz, 1

        adm_mom(j)= adm_mom(j) &
            - ( nu_loc(a)*amu/Msun )*( shift_norm2/(lapse_loc**two) - one ) &
              *theta_loc(a)*( one + u_loc(a) + pr_loc(a)/nlrf_loc(a) )*v_l(j,a)

      ENDDO

    ENDDO
    !$OMP END PARALLEL DO

    adm_mom= zero
    matter_objects_loop: DO i_matter= 1, parts% get_n_matter(), 1

      adm_mom_i(i_matter,:)= zero
      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( parts, g4_ll, vel_loc, &
      !$OMP                     v_l, nu_loc, theta_loc, pr_loc, nlrf_loc, &
      !$OMP                     u_loc, i_matter ) &
      !$OMP             PRIVATE( a, det, v_u, shift_norm2, j, g3, lapse_loc, &
      !$OMP                      shift ) &
      !$OMP             REDUCTION( +: adm_mom_i )
      DO a= parts% get_npart_i(i_matter-1) + 1, &
            parts% get_npart_i(i_matter-1) + parts% get_npart_i(i_matter), 1

        CALL compute_tpo_metric( g4_ll(:,a), lapse_loc, shift, g3 )

        v_u= [one, vel_loc(1:3,a)]
        CALL lower_index_4vector( v_u, g4_ll(:,a), v_l(:,a) )

        CALL spatial_vector_norm_sym3x3( g3, shift, shift_norm2 )

        DO j= jx, jz, 1

          adm_mom_i(i_matter,j)= adm_mom_i(i_matter,j) &
            - ( nu_loc(a)*amu/Msun )*( shift_norm2/(lapse_loc**two) - one ) &
              *theta_loc(a)*( one + u_loc(a) + pr_loc(a)/nlrf_loc(a) )*v_l(j,a)

        ENDDO

      ENDDO
      !$OMP END PARALLEL DO

      PRINT *, "   Estimate of the ADM momentum of matter object", i_matter, &
               " computed using the SPH hydro fields and", &
               " the metric mapped with mesh-to-particle mapping= "
      PRINT *, "   (", adm_mom_i(i_matter,1), ","
      PRINT *, "    ", adm_mom_i(i_matter,2), ","
      PRINT *, "    ", adm_mom_i(i_matter,3), ") Msun*c"
      PRINT *

      adm_mom= adm_mom + adm_mom_i(i_matter,:)

    ENDDO matter_objects_loop

    PRINT *, "   Estimate of the ADM momentum of the fluid computed using ",&
             "the SPH hydro fields and the metric mapped with ", &
             "mesh-to-particle mapping= "
    PRINT *, "   (", adm_mom(1), ","
    PRINT *, "    ", adm_mom(2), ","
    PRINT *, "    ", adm_mom(3), ") Msun*c"
    PRINT *

    CALL deallocate_metric_on_particles
    CALL deallocate_gradient
    DEALLOCATE(alive)
    CALL deallocate_sph_memory

    CALL deallocate_grid_function ( rad_coord, 'rad_coord' )
    CALL deallocate_ADM()
    CALL deallocate_Ztmp()
    !CALL deallocate_GravityAcceleration()
    CALL deallocate_BSSN()
    DEALLOCATE( levels )
    DEALLOCATE( pos     )
    DEALLOCATE( vel_loc )
    DEALLOCATE( nu_loc    )
    DEALLOCATE( nlrf_loc  )
    DEALLOCATE( u_loc     )
    DEALLOCATE( pr_loc    )
    DEALLOCATE( theta_loc )

  END PROCEDURE compute_adm_momentum_fluid_m2p


END SUBMODULE sph_adm_variables
