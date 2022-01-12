! File:         submodule_ejecta_generic_constructor.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (ejecta_generic) ejecta_generic_constructor

  !*********************************************************
  !
  !# Implementation of the constructor and
  !  destructor of TYPE [[ejecta]], and of the
  !  [[ejecta]]-member
  !  PROCEDURES that call the C-bound PROCEDURES
  !  constructig and destructing the |lorene|
  !  |etdiffrot| object
  !
  !  FT 19.11.2021
  !
  !*********************************************************


  IMPLICIT NONE


  CONTAINS


  !
  !-- Implementation of the constructor of the bns object
  !
  MODULE PROCEDURE construct_ejecta

    !****************************************************
    !
    !# Constructs an object of TYPE [[ejecta]]
    !
    !  FT 19.11.2021
    !
    !****************************************************

    USE constants, ONLY: pi
    USE NR,        ONLY: indexx
    USE pwp_EOS,   ONLY: get_Gamma0, get_Gamma1, get_Gamma2, get_Gamma3, &
                         get_K0, get_K1, get_K2, get_K3, get_p1, &
                         get_rho_0, get_rho_1, get_rho_2, select_EOS_parameters

    IMPLICIT NONE


    INTEGER, PARAMETER:: unit_pos= 2589
    DOUBLE PRECISION, PARAMETER:: atmosphere_density= 1.0439859633622731D-17

    INTEGER:: header_lines= 2 ! TODO: give this as input
    INTEGER:: nlines, ntmp
    INTEGER:: i_matter, n_matter_loc, itr, i, j, k

    DOUBLE PRECISION:: xtmp, ytmp, ztmp, rhotmp, epstmp, vxtmp, vytmp, vztmp, &
                       dr, dphi, dth
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: grid_tmp
    DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE:: rho_tmp
    INTEGER, DIMENSION(:), ALLOCATABLE:: x_sorted, y_sorted, z_sorted

    LOGICAL:: exist

    CHARACTER(LEN=:), ALLOCATABLE:: finalnamefile

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: mass_profile
    INTEGER, DIMENSION(:), ALLOCATABLE:: mass_profile_idx

    PRINT *, " * Reading ejecta ID from formatted file " &
             // TRIM(filename)
    PRINT *

    CALL derived_type% set_n_matter(1) ! TODO: give this as argument
    n_matter_loc= derived_type% get_n_matter()

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
    grid_tmp= 0.0D0

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

    !DO itr= 1, derived_type% n_gridpoints, 1
    !
    !  grid_tmp( derived_type% n_gridpoints + itr, 1 )=   grid_tmp( itr, 1 )
    !  grid_tmp( derived_type% n_gridpoints + itr, 2 )=   grid_tmp( itr, 2 )
    !  grid_tmp( derived_type% n_gridpoints + itr, 3 )= - grid_tmp( itr, 3 )
    !
    !ENDDO

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

    derived_type% nx_grid= NINT( ( MAXVAL( grid_tmp(:,1) ) - derived_type% xL_grid ) &
                           /derived_type% dx_grid + 1 )
    derived_type% ny_grid= NINT( ( MAXVAL( grid_tmp(:,2) ) - derived_type% yL_grid ) &
                           /derived_type% dy_grid + 1 )
    derived_type% nz_grid= NINT( ( MAXVAL( grid_tmp(:,3) ) - derived_type% zL_grid ) &
                           /derived_type% dz_grid + 1 )

    PRINT *, ( MAXVAL( grid_tmp(:,1) ) - derived_type% xL_grid )
    PRINT *, ( MAXVAL( grid_tmp(:,1) ) - derived_type% xL_grid )/derived_type% dx_grid
    PRINT *
    PRINT *, MAXVAL( grid_tmp(:,1) )
    PRINT *, MAXVAL( grid_tmp(:,2) )
    PRINT *, MAXVAL( grid_tmp(:,3) )
    PRINT *
    PRINT *, MINVAL( grid_tmp(:,1) ), derived_type% xL_grid
    PRINT *, MINVAL( grid_tmp(:,2) ), derived_type% yL_grid
    PRINT *, MINVAL( grid_tmp(:,3), grid_tmp(:,3) /= 0 ), derived_type% zL_grid
    PRINT *, MAXVAL( grid_tmp(:,1) ), derived_type% xR_grid
    PRINT *, MAXVAL( grid_tmp(:,2) ), derived_type% yR_grid
    PRINT *, MAXVAL( grid_tmp(:,3) ), derived_type% zR_grid
    PRINT *
    PRINT *, derived_type% dx_grid
    PRINT *, derived_type% dy_grid
    PRINT *, derived_type% dz_grid
    PRINT *
    PRINT *, derived_type% nx_grid
    PRINT *, derived_type% ny_grid
    PRINT *, derived_type% nz_grid
    PRINT *

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

  !  finalnamefile= "pos_ejecta.dat"
  !
  !  INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )
  !
  !  IF( exist )THEN
  !    OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
  !          FORM= "FORMATTED", &
  !          POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
  !          IOMSG= err_msg )
  !  ELSE
  !    OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "NEW", &
  !          FORM= "FORMATTED", &
  !          ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
  !  ENDIF
  !  IF( ios > 0 )THEN
  !    PRINT *, "...error when opening " // TRIM(finalnamefile), &
  !             ". The error message is", err_msg
  !    STOP
  !  ENDIF
  !
  !  DO i= 1, 2*derived_type% n_gridpoints, 1
  !
  !    WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
  !      grid_tmp(i,1), grid_tmp(i,2), grid_tmp(i,3)
  !
  !  ENDDO
  !
  !  CLOSE( UNIT= 2 )


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

    ! Store the grid on the member array
    ALLOCATE( derived_type% grid( derived_type% nx_grid, &
                                  derived_type% ny_grid, &
                                  derived_type% nz_grid, 3 ) )
    derived_type% grid= 0.0D0

    ALLOCATE( derived_type% baryon_mass_density( derived_type% nx_grid, &
                                                 derived_type% ny_grid, &
                                                 derived_type% nz_grid ) )
    derived_type% baryon_mass_density= 0.0D0

    ALLOCATE( derived_type% specific_energy( derived_type% nx_grid, &
                                             derived_type% ny_grid, &
                                             derived_type% nz_grid ) )
    derived_type% specific_energy= 0.0D0

    ALLOCATE( derived_type% vel( derived_type% nx_grid, &
                                 derived_type% ny_grid, &
                                 derived_type% nz_grid, 3 ) )
    derived_type% vel= 0.0D0

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

    ! Get rid of the athmosphere
    DO i= 1, derived_type% nx_grid, 1
      DO j= 1, derived_type% ny_grid, 1
        DO k= 1, derived_type% nz_grid, 1

          IF( derived_type% baryon_mass_density( i, j, k ) &
              <= atmosphere_density )THEN

            derived_type% baryon_mass_density( i, j, k )= 0.0D0

            derived_type% specific_energy( i, j, k )= 0.0D0

            derived_type% vel( i, j, k, : )= 0.0D0

          ENDIF

          !derived_type% baryon_mass_density( i, j, k )= &
          !  MAX( 0.0D0, &
          !  derived_type% baryon_mass_density( i, j, k ) - atmosphere_density )
          !
          !IF( derived_type% baryon_mass_density( i, j, k ) == 0.0D0 )THEN
          !
          !  derived_type% specific_energy( i, j, k )= 0.0D0
          !
          !  derived_type% vel( i, j, k, : )= 0.0D0
          !
          !ENDIF

        ENDDO
      ENDDO
    ENDDO

    PRINT *, derived_type% grid( 1, 1, 1, : )
    PRINT *, derived_type% grid( 2, 1, 1, : )
    PRINT *, derived_type% grid( 1, 2, 1, : )
    PRINT *, derived_type% grid( 1, 1, 2, : )
    PRINT *

    ALLOCATE( derived_type% masses(n_matter_loc) )
    ALLOCATE( derived_type% sizes(n_matter_loc,6) )
    ALLOCATE( derived_type% centers(n_matter_loc,3) )
    ALLOCATE( derived_type% barycenters(n_matter_loc,3) )

    DO i_matter= 1, n_matter_loc, 1

      !derived_type% masses(i_matter)= 0.0D0
      derived_type% centers(i_matter,:)= 0.0D0
      derived_type% barycenters(i_matter,:)= 0.0D0
      derived_type% sizes(i_matter,:)= [ &
                  !1.3D0*SQRT( ABS(MAXVAL(grid_tmp( :, 1 )))**2.0D0 &
                  !    + ABS(MAXVAL(grid_tmp( :, 2 )))**2.0D0 ), &
                  !1.3D0*SQRT( ABS(MAXVAL(grid_tmp( :, 1 )))**2.0D0 &
                  !    + ABS(MAXVAL(grid_tmp( :, 2 )))**2.0D0 ), &
                  !1.3D0*SQRT( ABS(MAXVAL(grid_tmp( :, 1 )))**2.0D0 &
                  !    + ABS(MAXVAL(grid_tmp( :, 2 )))**2.0D0 ), &
                  !1.3D0*SQRT( ABS(MAXVAL(grid_tmp( :, 1 )))**2.0D0 &
                  !    + ABS(MAXVAL(grid_tmp( :, 2 )))**2.0D0 ), &
                  !1.3D0*SQRT( ABS(MAXVAL(grid_tmp( :, 1 )))**2.0D0 &
                  !    + ABS(MAXVAL(grid_tmp( :, 3 )))**2.0D0 ), &
                  !1.3D0*SQRT( ABS(MAXVAL(grid_tmp( :, 1 )))**2.0D0 &
                  !    + ABS(MAXVAL(grid_tmp( :, 3 )))**2.0D0 ) ]
                                         ABS(derived_type% xL_grid), &
                                         ABS(MAXVAL(grid_tmp( :, 1 ))), &
                                         ABS(derived_type% yL_grid), &
                                         ABS(MAXVAL(grid_tmp( :, 2 ))), &
                                         ABS(MAXVAL(grid_tmp( :, 3 ))), &
                                         ABS(MAXVAL(grid_tmp( :, 3 ))) ]

    ENDDO

    finalnamefile= "pos_ejecta.dat"

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
              derived_type% grid( i, j, k, 1 ) + derived_type% dx_grid/2.0D0, &
              derived_type% grid( i, j, k, 2 ), &
              derived_type% grid( i, j, k, 3 ) ), &
            derived_type% grid( i, j, k, 1 ) + derived_type% dx_grid/2.0D0, &
            derived_type% specific_energy( i, j, k )
        ENDDO
      ENDDO
    ENDDO

    CLOSE( UNIT= 2 )

    dr             = derived_type% dx_grid/4.0D0
    dth            = pi/2.0D0/100.0D0
    dphi           = 2.0D0*pi/100.0D0
    CALL derived_type% integrate_baryon_mass_density( &
                            derived_type% centers(1,1), &
                            ABS(MAXVAL(grid_tmp( :, 1 ))), &
                            derived_type% read_mass_density( &
                              derived_type% centers(1,1), &
                              derived_type% centers(1,2), &
                              derived_type% centers(1,3) ), &
                            dr, dth, dphi, &
                            derived_type% masses(1), mass_profile, &
                            mass_profile_idx )

    derived_type% eos_ejectaid= 110

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

    !PRINT *, "sizes=", derived_type% sizes(1,:)

    ! Determine nx_grid, ny_grid, nz_grid TODO: parallelize everything

!    ALLOCATE( x_sorted( 2*derived_type% n_gridpoints ) )
!    ALLOCATE( y_sorted( 2*derived_type% n_gridpoints ) )
!    ALLOCATE( z_sorted( 2*derived_type% n_gridpoints ) )
!    x_sorted= 0.0D0
!    y_sorted= 0.0D0
!    z_sorted= 0.0D0
!
!    ! Sort the x coordinates of the grid
!    CALL indexx( 2*derived_type% n_gridpoints, grid_tmp( :, 1 ), x_sorted )
!
!    ! Count how many different x coordinates there are
!    derived_type% nx_grid= 1
!    DO itr= 1, 2*derived_type% n_gridpoints - 1, 1
!
!      IF( grid_tmp(x_sorted(itr),1) /= grid_tmp(x_sorted(itr+1),1) )THEN
!        derived_type% nx_grid= derived_type% nx_grid + 1
!      ENDIF
!
!    ENDDO
!
!    ! Sort the y coordinates of the grid
!    CALL indexx( 2*derived_type% n_gridpoints, grid_tmp( :, 2 ), y_sorted )
!
!    ! Count how many different y coordinates there are
!    derived_type% ny_grid= 1
!    DO itr= 1, 2*derived_type% n_gridpoints - 1, 1
!
!      IF( grid_tmp(y_sorted(itr),2) /= grid_tmp(y_sorted(itr+1),2) )THEN
!        derived_type% ny_grid= derived_type% ny_grid + 1
!      ENDIF
!
!    ENDDO
!
!    ! Sort the z coordinates of the grid
!    CALL indexx( 2*derived_type% n_gridpoints, grid_tmp( :, 3 ), z_sorted )
!
!    ! Count how many different z coordinates there are
!    derived_type% nz_grid= 1
!    DO itr= 1, derived_type% n_gridpoints - 1, 1
!
!      IF( grid_tmp(z_sorted(itr),3) /= grid_tmp(z_sorted(itr+1),3) )THEN
!        derived_type% nz_grid= derived_type% nz_grid + 1
!      ENDIF
!
!    ENDDO
!    derived_type% nz_grid= 2*derived_type% nz_grid
!
!    ! Check that the grid dimensions are consistent
!    IF( derived_type% nx_grid*derived_type% ny_grid*derived_type% nz_grid &
!        /= 2*derived_type% n_gridpoints )THEN
!
!      PRINT *, derived_type% nx_grid
!      PRINT *, derived_type% ny_grid
!      PRINT *, derived_type% nz_grid
!      PRINT *, derived_type% nx_grid*derived_type% ny_grid*derived_type% nz_grid
!      PRINT *, 2*derived_type% n_gridpoints
!      PRINT *, derived_type% n_gridpoints
!      STOP
!
!    ENDIF
!
!    derived_type% dx_grid= (MAXVAL(grid_tmp( :, 1 ))-MINVAL(grid_tmp( :, 1 )))/(derived_type% nx_grid)
!    derived_type% dy_grid= (MAXVAL(grid_tmp( :, 2 ))-MINVAL(grid_tmp( :, 2 )))/(derived_type% ny_grid)
!    derived_type% dz_grid= (MAXVAL(grid_tmp( :, 3 ))-MINVAL(grid_tmp( :, 3 )))/(derived_type% nz_grid)
!
!
!
!    derived_type% n_gridpoints=2.0D0*derived_type% n_gridpoints
!
!    ! Check that nx dx and the grid extent are consistent
!    ztmp= derived_type% zL_grid
!    DO k= 1, derived_type% nz_grid, 1
!      ztmp= ztmp + derived_type% dz_grid
!    ENDDO
!    IF( ABS( ztmp - MAXVAL(grid_tmp( :, 3 )) ) &
!        > derived_type% dz_grid/1.0D+6 )THEN
!      PRINT *, "** ERROR! ztmp=", ztmp
!      PRINT *, "          zR_grid=", MAXVAL(grid_tmp( :, 3 ))
!      STOP
!    ENDIF
!    ytmp= derived_type% yL_grid
!    DO j= 1, derived_type% ny_grid, 1
!      ytmp= ytmp + derived_type% dy_grid
!    ENDDO
!    IF( ABS( ytmp - MAXVAL(grid_tmp( :, 2 )) ) &
!        > derived_type% dz_grid/1.0D+6 )THEN
!      PRINT *, "** ERROR! ytmp=", ytmp
!      PRINT *, "          yR_grid=", MAXVAL(grid_tmp( :, 2 ))
!      STOP
!    ENDIF
!    xtmp= derived_type% xL_grid
!    DO i= 1, derived_type% nx_grid, 1
!      xtmp= xtmp + derived_type% dx_grid
!    ENDDO
!    IF( ABS( xtmp - MAXVAL(grid_tmp( :, 1 )) ) &
!        > derived_type% dx_grid/1.0D+6 )THEN
!      PRINT *, "** ERROR! xtmp=", xtmp
!      PRINT *, "          xR_grid=", MAXVAL(grid_tmp( :, 1 ))
!      STOP
!    ENDIF
!
!    PRINT *, derived_type% nx_grid
!    PRINT *, derived_type% ny_grid
!    PRINT *, derived_type% nz_grid
!    PRINT *
!    PRINT *, MAXVAL( grid_tmp(:,1) )
!    PRINT *, MAXVAL( grid_tmp(:,2) )
!    PRINT *, MAXVAL( grid_tmp(:,3) )
!    PRINT *
!    PRINT *, MINVAL( grid_tmp(:,1) )
!    PRINT *, MINVAL( grid_tmp(:,2) )
!    PRINT *, MINVAL( grid_tmp(:,3) )
!    PRINT *
!    PRINT *, derived_type% xL_grid
!    PRINT *, derived_type% yL_grid
!    PRINT *, derived_type% zL_grid
!    PRINT *
!    PRINT *, grid_tmp(1,1)
!    PRINT *, grid_tmp(1,1)
!    PRINT *, grid_tmp(1,1)
!    PRINT *
!    PRINT *, derived_type% dx_grid
!    PRINT *, derived_type% dy_grid
!    PRINT *, derived_type% dz_grid
!    PRINT *
!
!    !STOP
!
!    !derived_type% nz_grid= 2.0D0*derived_type% nz_grid
!
!    ! Store the grid on the member array
!    ALLOCATE( derived_type% grid(derived_type% nx_grid, &
!                                 derived_type% ny_grid, &
!                                 derived_type% nz_grid, 3) )
!  !
!  !  DO k= 1, derived_type% nz_grid, 1
!  !    DO j= 1, derived_type% ny_grid, 1
!  !      DO i= 1, derived_type% nx_grid, 1
!  !
!  !        derived_type% grid(i,j,k,:)= &
!  !              grid_tmp((k-1)*(derived_type% ny_grid)*(derived_type% nx_grid) &
!  !                     + (j-1)*(derived_type% nx_grid) + i, :)
!  !
!  !      ENDDO
!  !    ENDDO
!  !  ENDDO
!
!
!    STOP
!
!    ! Read the ID
!    OPEN( UNIT= unit_pos, FILE= TRIM(filename), &
!          FORM= "FORMATTED", ACTION= "READ" )
!
!    ! Skip header
!    DO itr= 1, header_lines, 1
!      READ( unit_pos, * )
!    ENDDO
!
!    ! Allocate the arrays to store data
!    ALLOCATE( derived_type% baryon_mass_density( derived_type% nx_grid, &
!                                                 derived_type% ny_grid, &
!                                                 derived_type% nz_grid ) )
!    derived_type% baryon_mass_density= 0.0D0
!
!    ! Read the data into the array
!    DO i= 1, derived_type% nx_grid, 1
!      DO j= 1, derived_type% ny_grid, 1
!        DO k= 1, derived_type% nz_grid/2.0D0, 1
!
!        READ( UNIT= unit_pos, FMT= *, IOSTAT = ios, IOMSG= err_msg ) &
!          xtmp, ytmp, ztmp, rhotmp
!
!          IF( ztmp > 0 )THEN
!
!            derived_type% grid( i, j, derived_type% nz_grid/2.0D0 + k, 1 )= xtmp
!            derived_type% grid( i, j, derived_type% nz_grid/2.0D0 + k, 2 )= ytmp
!            derived_type% grid( i, j, derived_type% nz_grid/2.0D0 + k, 3 )= ztmp
!            derived_type% baryon_mass_density( i, j, derived_type% nz_grid/2.0D0 + k )= rhotmp
!
!          ENDIF
!
!          IF( ios > 0 )THEN
!            PRINT *, "...error when reading " // TRIM(filename), &
!                    " at particle ", itr,". The status variable is ", ios, &
!                    ". The error message is", err_msg
!            STOP
!          ENDIF
!
!        ENDDO
!      ENDDO
!    ENDDO
!    DO i= 1, derived_type% nx_grid, 1
!      DO j= 1, derived_type% ny_grid, 1
!        DO k= 1, derived_type% nz_grid/2.0D0, 1
!
!          derived_type% grid( i, j, k, 1 )= &
!                    derived_type% grid( i, j, derived_type% nz_grid - k+1, 1 )
!          derived_type% grid( i, j, k, 2 )= &
!                    derived_type% grid( i, j, derived_type% nz_grid - k+1, 2 )
!          derived_type% grid( i, j, k, 3 )= &
!                  - derived_type% grid( i, j, derived_type% nz_grid - k+1, 3 )
!          derived_type% baryon_mass_density( i, j, k )= &
!        derived_type% baryon_mass_density( i, j, derived_type% nz_grid - k+1 )
!
!        ENDDO
!      ENDDO
!    ENDDO
!
!    CLOSE( UNIT= unit_pos )
!
!    !derived_type% zL_grid= MINVAL(derived_type% grid( 1, 1, :, 3 ))
!    !ztmp= derived_type% zL_grid
!    !DO k= 1, derived_type% nz_grid, 1
!    !  ztmp= ztmp + derived_type% dz_grid
!    !ENDDO
!    !IF( ABS( ztmp - MAXVAL(derived_type% grid( 1, 1, :, 3 )) ) &
!    !    > derived_type% dz_grid/1.0D+6 )THEN
!    !  PRINT *, "** ERROR! ztmp=", ztmp
!    !  PRINT *, "          zR_grid=", MAXVAL(derived_type% grid( 1, 1, :, 3 ))
!    !  STOP
!    !ENDIF
!
!    PRINT *, derived_type% nx_grid
!    PRINT *, derived_type% ny_grid
!    PRINT *, derived_type% nz_grid
!    PRINT *
!    PRINT *, SIZE( derived_type% grid( :, 1, 1, 1 ) )
!    PRINT *, SIZE( derived_type% grid( 1, :, 1, 1 ) )
!    PRINT *, SIZE( derived_type% grid( 1, 1, :, 1 ) )
!    PRINT *
!    PRINT *, MAXVAL( derived_type% grid( :, 1, 1, 1 ) )
!    PRINT *, MAXVAL( derived_type% grid( 1, :, 1, 2 ) )
!    PRINT *, MAXVAL( derived_type% grid( 1, 1, :, 3 ) )
!    PRINT *
!    PRINT *, MINVAL( derived_type% grid( :, 1, 1, 1 ) )
!    PRINT *, MINVAL( derived_type% grid( 1, :, 1, 2 ) )
!    PRINT *, MINVAL( derived_type% grid( 1, 1, :, 3 ) )
!    PRINT *
!    PRINT *, derived_type% xL_grid
!    PRINT *, derived_type% yL_grid
!    PRINT *, derived_type% zL_grid
!    PRINT *
!    PRINT *, derived_type% grid( 1,1,1, 1 )
!    PRINT *, derived_type% grid( 1,1,1, 2 )
!    PRINT *, derived_type% grid( 1,1,1, 3 )
!    PRINT *
!    PRINT *, derived_type% dx_grid
!    PRINT *, derived_type% dy_grid
!    PRINT *, derived_type% dz_grid
!    PRINT *
!
!    derived_type% xL_grid= MINVAL( derived_type% grid( :, 1, 1, 1 ) )
!    derived_type% yL_grid= MINVAL( derived_type% grid( 1, :, 1, 2 ) )
!    derived_type% zL_grid= MINVAL( derived_type% grid( 1, 1, :, 3 ) )
!
!    ALLOCATE( rho_tmp(derived_type% nx_grid, derived_type% ny_grid, derived_type% nz_grid) )
!    rho_tmp= 0.0D0
!
!    PRINT *, derived_type% read_mass_density( &
!       derived_type% grid( 1,1,1, 1 ) + derived_type% dx_grid/2.0D0, &
!       derived_type% grid( 1,1,1, 2 ) + derived_type% dy_grid/2.0D0, &
!       derived_type% grid( 1,1,1, 3 ) + derived_type% dz_grid/2.0D0 &
!     )
!
!    STOP
!
!    finalnamefile= "pos_ejecta.dat"
!
!    INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )
!
!    IF( exist )THEN
!        OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
!              FORM= "FORMATTED", &
!              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
!              IOMSG= err_msg )
!    ELSE
!        OPEN( UNIT= 2, FILE= TRIM(finalnamefile), STATUS= "NEW", &
!              FORM= "FORMATTED", &
!              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
!    ENDIF
!    IF( ios > 0 )THEN
!      PRINT *, "...error when opening " // TRIM(finalnamefile), &
!               ". The error message is", err_msg
!      STOP
!    ENDIF
!
!    DO i= 1, derived_type% nx_grid - 1, 1
!      DO j= 1, derived_type% ny_grid - 1, 1
!        DO k= 1, derived_type% nz_grid - 1, 1
!
!          WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
!            derived_type% grid( i, j, k, 1 ), &
!            derived_type% grid( i, j, k, 2 ), &
!            derived_type% grid( i, j, k, 3 )
!
!        ENDDO
!      ENDDO
!    ENDDO
!
!    CLOSE( UNIT= 2 )
!
!    STOP

  END PROCEDURE construct_ejecta


  !
  !-- Implementation of the destructor of the bns object
  !
  MODULE PROCEDURE destruct_ejecta

    !****************************************************
    !
    !# Destructs an object of TYPE [[ejecta]]
    !
    !  FT 19.11.2021
    !
    !****************************************************

    IMPLICIT NONE


  END PROCEDURE destruct_ejecta


END SUBMODULE ejecta_generic_constructor
