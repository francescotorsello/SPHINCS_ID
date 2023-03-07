! File:         submodule_ejecta_generic_constructor.f90
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

SUBMODULE (ejecta_generic) constructor

  !*********************************************************
  !
  !# Implementation of the constructor and
  !  destructor of TYPE [[ejecta]]
  !
  !  FT xx.11.2021
  !
  !*********************************************************


  IMPLICIT NONE


  CONTAINS


  !
  !-- Implementation of the constructor of the ejecta object
  !
  MODULE PROCEDURE construct_ejecta

    !****************************************************
    !
    !# Constructs an object of TYPE [[ejecta]]
    !  @todo to be OMP parallelized
    !
    !  FT xx.11.2021
    !
    !****************************************************

    USE constants, ONLY: pi
    USE utility,   ONLY: zero, one, two, four, ten, eos$poly, eos$pwpoly, &
                         sph_path
    USE NR,        ONLY: indexx
    USE pwp_EOS,   ONLY: get_Gamma0, get_Gamma1, get_Gamma2, get_Gamma3, &
                         get_K0, get_K1, get_K2, get_K3, get_p1, &
                         get_rho_0, get_rho_1, get_rho_2, select_EOS_parameters
    USE timing,    ONLY: timer

    IMPLICIT NONE


    INTEGER, PARAMETER:: unit_pos= 2589
    DOUBLE PRECISION, PARAMETER:: atmosphere_density= 1.0439859633622731D-17

    INTEGER:: header_lines= 2 ! TODO: give this as input
    INTEGER:: nlines, ntmp
    INTEGER:: i_matter, n_matter_loc, itr, i, j, k
   ! INTEGER, DIMENSION(:), ALLOCATABLE:: x_sorted, y_sorted, z_sorted
    INTEGER, DIMENSION(:), ALLOCATABLE:: mass_profile_idx

    DOUBLE PRECISION:: xtmp, ytmp, ztmp, &
                       rhotmp, epstmp, vxtmp, vytmp, vztmp, &
                       dr, dphi, dth
    DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE:: grid_tmp
    !DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE:: rho_tmp
    DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE:: mass_profile

    LOGICAL:: exist

    CHARACTER(LEN=:), ALLOCATABLE:: finalnamefile

    PRINT *, " * Reading ID on Cartesian, uniform grid from formatted file " &
             // TRIM(filename), "..."

    CALL derived_type% set_n_matter(1) ! TODO: give this as argument
    n_matter_loc= derived_type% get_n_matter()

    CALL derived_type% set_cold_system(.FALSE.)

    derived_type% construction_timer= timer( "ejecta_construction_timer" )

    CALL derived_type% construction_timer% start_timer()

    INQUIRE( FILE= TRIM(filename), EXIST= exist )

    IF( exist )THEN
      OPEN( UNIT= unit_pos, FILE= TRIM(filename), &
            FORM= "FORMATTED", ACTION= "READ", IOSTAT= ios, &
            IOMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...error when opening " // TRIM(filename), &
                ". The error message is", err_msg
        STOP
      ENDIF
    ELSE
      PRINT *, "** ERROR! Unable to find file " // TRIM(filename)
      STOP
    ENDIF

    ! Get total number of lines in the file
    nlines = 0
    DO
      READ( unit_pos, * , IOSTAT= ios )
      IF ( ios /= 0 ) EXIT
      nlines = nlines + 1
    ENDDO

    CLOSE( UNIT= unit_pos )

    derived_type% n_gridpoints= nlines - header_lines

    ! Allocate the temporary array to store data
    ALLOCATE( grid_tmp( 2*derived_type% n_gridpoints, 8 ) )
    grid_tmp= zero

    ! Read the ID
    OPEN( UNIT= unit_pos, FILE= TRIM(filename), &
          FORM= "FORMATTED", ACTION= "READ" )

    ! Skip header
    DO itr= 1, header_lines, 1
      READ( unit_pos, * )
    ENDDO

    ! Read the data into the temporary array
    ntmp= 0
    DO itr= 1, derived_type% n_gridpoints, 1

      READ( UNIT= unit_pos, FMT= *, IOSTAT = ios, IOMSG= err_msg ) &
        xtmp, ytmp, ztmp, rhotmp, epstmp, vxtmp, vytmp, vztmp

      IF( ztmp > 0 )THEN
        ntmp= ntmp + 1
        grid_tmp( ntmp, 1 )= xtmp
        grid_tmp( ntmp, 2 )= ytmp
        grid_tmp( ntmp, 3 )= ztmp
        grid_tmp( ntmp, 4 )= rhotmp
        grid_tmp( ntmp, 5 )= epstmp
        grid_tmp( ntmp, 6 )= vxtmp
        grid_tmp( ntmp, 7 )= vytmp
        grid_tmp( ntmp, 8 )= vztmp
      ENDIF

      IF( ios > 0 )THEN
        PRINT *, "...error when reading " // TRIM(filename), &
                " at particle ", itr,". The status variable is ", ios, &
                ". The error message is", err_msg
        STOP
      ENDIF

    ENDDO
    derived_type% n_gridpoints= ntmp
    grid_tmp= grid_tmp( 1:2*derived_type% n_gridpoints, : )

    CLOSE( UNIT= unit_pos )

    DO itr= 1, SIZE(grid_tmp(:,1)), 1
      IF( grid_tmp(itr,1) > grid_tmp(1,1) )THEN
        derived_type% dx_grid= grid_tmp(itr,1) - grid_tmp(1,1)
        EXIT
      ENDIF
    ENDDO
    DO itr= 1, SIZE(grid_tmp(:,2)), 1
      IF( grid_tmp(itr,2) > grid_tmp(1,2) )THEN
        derived_type% dy_grid= grid_tmp(itr,2) - grid_tmp(1,2)
        EXIT
      ENDIF
    ENDDO
    DO itr= 1, SIZE(grid_tmp(:,3)), 1
      IF( grid_tmp(itr,3) > grid_tmp(1,3) )THEN
        derived_type% dz_grid= grid_tmp(itr,3) - grid_tmp(1,3)
        EXIT
      ENDIF
    ENDDO

    derived_type% xL_grid= MINVAL(grid_tmp( :, 1 ))
    derived_type% yL_grid= MINVAL(grid_tmp( :, 2 ))
    derived_type% zL_grid= MINVAL(grid_tmp( :, 3 ), grid_tmp( :, 3 ) /= 0 )
    derived_type% xR_grid= MAXVAL(grid_tmp( :, 1 ))
    derived_type% yR_grid= MAXVAL(grid_tmp( :, 2 ))
    derived_type% zR_grid= MAXVAL(grid_tmp( :, 3 ) )

    derived_type% nx_grid= NINT((MAXVAL(grid_tmp(:,1)) - derived_type% xL_grid)&
                                /derived_type% dx_grid + 1)
    derived_type% ny_grid= NINT((MAXVAL(grid_tmp(:,2)) - derived_type% yL_grid)&
                                /derived_type% dy_grid + 1 )
    derived_type% nz_grid= NINT((MAXVAL(grid_tmp(:,3)) - derived_type% zL_grid)&
                                /derived_type% dz_grid + 1 )

    PRINT *, " * ID on Cartesian, uniform grid read."
    PRINT *
    PRINT *, "** Grid information:"
    PRINT *
    PRINT *, " * Grid size in x direction:", &
             derived_type% xL_grid, derived_type% xR_grid
    PRINT *, " * Grid size in y direction:", &
             derived_type% yL_grid, derived_type% yR_grid
    PRINT *, " * Grid size in z direction:", &
             derived_type% zL_grid, derived_type% zR_grid
    PRINT *
    PRINT *, " * Grid spacing in the x direction:", derived_type% dx_grid
    PRINT *, " * Grid spacing in the y direction:", derived_type% dy_grid
    PRINT *, " * Grid spacing in the z direction:", derived_type% dz_grid
    PRINT *
    PRINT *, " * Number of grid points in the x direction:", &
             derived_type% nx_grid
    PRINT *, " * Number of grid points in the y direction:", &
             derived_type% ny_grid
    PRINT *, " * Number of grid points in the z direction:", &
             derived_type% nz_grid
    PRINT *

    PRINT *, "** Checking that the grid dimensions are consistent with the " &
             //"number of lines in the file..."
    ! Check that the grid dimensions are consistent
    IF( derived_type% nx_grid*derived_type% ny_grid*derived_type% nz_grid &
        /= derived_type% n_gridpoints )THEN

      PRINT *, derived_type% nx_grid
      PRINT *, derived_type% ny_grid
      PRINT *, derived_type% nz_grid
      PRINT *, derived_type% nx_grid*derived_type% ny_grid*derived_type% nz_grid
      PRINT *, derived_type% n_gridpoints
      STOP

    ENDIF

    PRINT *, "** Checking that the number of grid points, the grid spacings " &
             //"and the grid sizes are consistent..."
    ! Check that nx dx and the grid extent are consistent
    ztmp= derived_type% zL_grid
    DO k= 1, derived_type% nz_grid - 1, 1
      ztmp= ztmp + derived_type% dz_grid
    ENDDO
    IF( ABS( ztmp - MAXVAL(grid_tmp( :, 3 )) ) &
        > derived_type% dz_grid/1.0D+6 )THEN
      PRINT *, "** ERROR! ztmp=", ztmp
      PRINT *, "          zR_grid=", MAXVAL(grid_tmp( :, 3 ))
      STOP
    ENDIF
    ytmp= derived_type% yL_grid
    DO j= 1, derived_type% ny_grid - 1, 1
      ytmp= ytmp + derived_type% dy_grid
    ENDDO
    IF( ABS( ytmp - MAXVAL(grid_tmp( :, 2 )) ) &
        > derived_type% dz_grid/1.0D+6 )THEN
      PRINT *, "** ERROR! ytmp=", ytmp
      PRINT *, "          yR_grid=", MAXVAL(grid_tmp( :, 2 ))
      STOP
    ENDIF
    xtmp= derived_type% xL_grid
    DO i= 1, derived_type% nx_grid - 1, 1
      xtmp= xtmp + derived_type% dx_grid
    ENDDO
    IF( ABS( xtmp - MAXVAL(grid_tmp( :, 1 )) ) &
        > derived_type% dx_grid/1.0D+6 )THEN
      PRINT *, "** ERROR! xtmp=", xtmp
      PRINT *, "          xR_grid=", MAXVAL(grid_tmp( :, 1 ))
      STOP
    ENDIF
    PRINT *

    ! Allocate and initialize member arrays
    CALL derived_type% allocate_gridid_memory( n_matter_loc )

    derived_type% grid= zero
    derived_type% baryon_mass_density= zero
    derived_type% specific_energy= zero
    derived_type% vel= zero

    ! Store the ID into the member arrays
    DO i= 1, derived_type% nx_grid, 1
      DO j= 1, derived_type% ny_grid, 1
        DO k= 1, derived_type% nz_grid, 1

          derived_type% grid( i, j, k, 1 )= &
            grid_tmp( (i-1)*(derived_type% ny_grid)*(derived_type% nz_grid) &
                      + (j-1)*(derived_type% nz_grid) + k, 1 )
          derived_type% grid( i, j, k, 2 )= &
            grid_tmp( (i-1)*(derived_type% ny_grid)*(derived_type% nz_grid) &
                      + (j-1)*(derived_type% nz_grid) + k, 2 )
          derived_type% grid( i, j, k, 3 )= &
            grid_tmp( (i-1)*(derived_type% ny_grid)*(derived_type% nz_grid) &
                      + (j-1)*(derived_type% nz_grid) + k, 3 )

          derived_type% baryon_mass_density( i, j, k )= &
            grid_tmp( (i-1)*(derived_type% ny_grid)*(derived_type% nz_grid) &
                      + (j-1)*(derived_type% nz_grid) + k, 4 )

          derived_type% specific_energy( i, j, k )= &
            grid_tmp( (i-1)*(derived_type% ny_grid)*(derived_type% nz_grid) &
                      + (j-1)*(derived_type% nz_grid) + k, 5 )

          derived_type% vel( i, j, k, 1 )= &
            grid_tmp( (i-1)*(derived_type% ny_grid)*(derived_type% nz_grid) &
                      + (j-1)*(derived_type% nz_grid) + k, 6 )
          derived_type% vel( i, j, k, 2 )= &
            grid_tmp( (i-1)*(derived_type% ny_grid)*(derived_type% nz_grid) &
                      + (j-1)*(derived_type% nz_grid) + k, 7 )
          derived_type% vel( i, j, k, 3 )= &
            grid_tmp( (i-1)*(derived_type% ny_grid)*(derived_type% nz_grid) &
                      + (j-1)*(derived_type% nz_grid) + k, 8 )

        ENDDO
      ENDDO
    ENDDO

    ! Get rid of the atmosphere coming from a mesh-based simulation, if present
    DO i= 1, derived_type% nx_grid, 1
      DO j= 1, derived_type% ny_grid, 1
        DO k= 1, derived_type% nz_grid, 1

          IF( derived_type% baryon_mass_density( i, j, k ) &
              <= atmosphere_density )THEN

            derived_type% baryon_mass_density( i, j, k )= zero

            derived_type% specific_energy( i, j, k )= zero

            derived_type% vel( i, j, k, : )= zero

          ENDIF

          !derived_type% baryon_mass_density( i, j, k )= &
          !  MAX( zero, &
          !  derived_type% baryon_mass_density( i, j, k ) - atmosphere_density )
          !
          !IF( derived_type% baryon_mass_density( i, j, k ) == zero )THEN
          !
          !  derived_type% specific_energy( i, j, k )= zero
          !
          !  derived_type% vel( i, j, k, : )= zero
          !
          !ENDIF

        ENDDO
      ENDDO
    ENDDO


    ! Assign ID properties to member arrays
    DO i_matter= 1, n_matter_loc, 1

      !derived_type% masses(i_matter)= zero
      derived_type% centers(i_matter,:)= zero
      derived_type% barycenters(i_matter,:)= zero
      derived_type% sizes(i_matter,:)= [ &
                  !1.3D0*SQRT( ABS(MAXVAL(grid_tmp( :, 1 )))**two &
                  !    + ABS(MAXVAL(grid_tmp( :, 2 )))**two ), &
                  !1.3D0*SQRT( ABS(MAXVAL(grid_tmp( :, 1 )))**two &
                  !    + ABS(MAXVAL(grid_tmp( :, 2 )))**two ), &
                  !1.3D0*SQRT( ABS(MAXVAL(grid_tmp( :, 1 )))**two &
                  !    + ABS(MAXVAL(grid_tmp( :, 2 )))**two ), &
                  !1.3D0*SQRT( ABS(MAXVAL(grid_tmp( :, 1 )))**two &
                  !    + ABS(MAXVAL(grid_tmp( :, 2 )))**two ), &
                  !1.3D0*SQRT( ABS(MAXVAL(grid_tmp( :, 1 )))**two &
                  !    + ABS(MAXVAL(grid_tmp( :, 3 )))**two ), &
                  !1.3D0*SQRT( ABS(MAXVAL(grid_tmp( :, 1 )))**two &
                  !    + ABS(MAXVAL(grid_tmp( :, 3 )))**two ) ]
                                         ABS(derived_type% xL_grid), &
                                         ABS(MAXVAL(grid_tmp( :, 1 ))), &
                                         ABS(derived_type% yL_grid), &
                                         ABS(MAXVAL(grid_tmp( :, 2 ))), &
                                         ABS(MAXVAL(grid_tmp( :, 3 ))), &
                                         ABS(MAXVAL(grid_tmp( :, 3 ))) ]

    ENDDO

    finalnamefile= TRIM(sph_path)//"pos_ejecta.dat"

    INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

    IF( exist )THEN
        OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
              FORM= "FORMATTED", &
              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
              IOMSG= err_msg )
    ELSE
        OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "NEW", &
              FORM= "FORMATTED", &
              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening " // TRIM(finalnamefile), &
               ". The error message is", err_msg
      STOP
    ENDIF

    DO i= 1, derived_type% nx_grid - 1, 1
      DO j= 1, derived_type% ny_grid - 1, 1
        DO k= 1, derived_type% nz_grid - 1, 1

          WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
            derived_type% grid( i, j, k, 1 ), &
            derived_type% grid( i, j, k, 2 ), &
            derived_type% grid( i, j, k, 3 ), &
            derived_type% baryon_mass_density( i, j, k ), &
            derived_type% read_mass_density( &
              derived_type% grid( i, j, k, 1 ) + derived_type% dx_grid/two, &
              derived_type% grid( i, j, k, 2 ), &
              derived_type% grid( i, j, k, 3 ) ), &
            derived_type% grid( i, j, k, 1 ) + derived_type% dx_grid/two, &
            derived_type% specific_energy( i, j, k )
        ENDDO
      ENDDO
    ENDDO

    CLOSE( UNIT= 2 )

    ! Assign total mass to member variable
    dr             = derived_type% dx_grid/four
    dth            = pi/two/ten*ten
    dphi           = two*pi/ten*ten

    ALLOCATE( mass_profile( 3, 0:NINT(ABS(MAXVAL(grid_tmp( :, 1 )))/dr) ), &
              STAT= ios, ERRMSG= err_msg )
    ALLOCATE( mass_profile_idx( 0:NINT(ABS(MAXVAL(grid_tmp( :, 1 )))/dr) ), &
              STAT= ios, ERRMSG= err_msg )

    CALL derived_type% integrate_baryon_mass_density( &
                            derived_type% centers(1,:), &
                            MAXVAL(grid_tmp( :, 1 )), &
                            derived_type% read_mass_density( &
                              derived_type% centers(1,1), &
                              derived_type% centers(1,2), &
                              derived_type% centers(1,3) ), &
                            dr, dth, dphi, &
                            derived_type% masses(1), mass_profile, &
                            mass_profile_idx )!, &
                            !radii= [MAXVAL(derived_type% sizes(1,3:4)), &
                            !        MAXVAL(derived_type% sizes(1,5:6))] )

    ! Set the EOS parameters
    derived_type% eos_id= eos$pwpoly

    CALL select_EOS_parameters("APR4")

    derived_type% npeos  = 3
    derived_type% gamma0 = get_Gamma0()
    derived_type% gamma1 = get_Gamma1()
    derived_type% gamma2 = get_Gamma2()
    derived_type% gamma3 = get_Gamma3()
    derived_type% kappa0 = get_K0()
    derived_type% kappa1 = get_K1()
    derived_type% kappa2 = get_K2()
    derived_type% kappa3 = get_K3()
    derived_type% logP1  = LOG10(get_p1())
    derived_type% logRho0= LOG10(get_rho_0())
    derived_type% logRho1= LOG10(get_rho_1())
    derived_type% logRho2= LOG10(get_rho_2())

    derived_type% finalize_sph_id_ptr => finalize

    CALL derived_type% construction_timer% stop_timer()

  END PROCEDURE construct_ejecta


  MODULE PROCEDURE finalize

    !***********************************************
    !
    !#
    !
    !  FT 14.04.2022
    !
    !***********************************************

    IMPLICIT NONE

    ! Temporary implementation, to avoid warnings about unused variables

    pos  = pos
    nlrf = nlrf
    nu   = nu
    pr   = pr
    vel_u= vel_u
    theta= theta
    nstar= nstar
    u    = u

  END PROCEDURE finalize


  MODULE PROCEDURE nothing

    !***********************************************
    !
    !# Procedure that does nothing. It is used to instantiate a deferred
    !  idbase procedure which s not needed in TYPE [[ejecta_generic]].
    !  It also serves as a placeholder in case the idbase procedure
    !  will be needed in the future.
    !
    !  FT 15.09.2022
    !
    !***********************************************

    IMPLICIT NONE

  END PROCEDURE nothing


  !
  !-- Implementation of the destructor of the bns object
  !
  MODULE PROCEDURE destruct_ejecta

    !****************************************************
    !
    !# Destructs an object of TYPE [[ejecta]]
    !
    !  FT xx.11.2021
    !
    !****************************************************

    IMPLICIT NONE

    CALL this% deallocate_gridid_memory()


  END PROCEDURE destruct_ejecta


END SUBMODULE constructor
