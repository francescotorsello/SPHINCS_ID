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

    USE NR, ONLY: indexx

    IMPLICIT NONE


    INTEGER, PARAMETER:: unit_pos= 2589

    INTEGER:: header_lines= 2 ! TODO: give this as input
    INTEGER:: nlines, ntmp
    INTEGER:: i_matter, n_matter_loc, itr, i, j, k

    DOUBLE PRECISION:: xtmp, ytmp, ztmp
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: grid_tmp
    DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE:: rho_tmp
    INTEGER, DIMENSION(:), ALLOCATABLE:: grid_sorted

    LOGICAL:: exist

    CHARACTER(LEN=:), ALLOCATABLE:: finalnamefile

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
    ALLOCATE( grid_tmp( 2*derived_type% n_gridpoints, 3 ) )
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
        xtmp, ytmp, ztmp

      IF( ztmp > 0 )THEN
        ntmp= ntmp + 1
        grid_tmp( ntmp, 1 )= xtmp
        grid_tmp( ntmp, 2 )= ytmp
        grid_tmp( ntmp, 3 )= ztmp
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

    DO itr= 1, derived_type% n_gridpoints, 1

      grid_tmp( derived_type% n_gridpoints + itr, 1 )=   grid_tmp( itr, 1 )
      grid_tmp( derived_type% n_gridpoints + itr, 2 )=   grid_tmp( itr, 2 )
      grid_tmp( derived_type% n_gridpoints + itr, 3 )= - grid_tmp( itr, 3 )

    ENDDO

    ALLOCATE( derived_type% masses(n_matter_loc) )
    ALLOCATE( derived_type% sizes(n_matter_loc,6) )
    ALLOCATE( derived_type% centers(n_matter_loc,3) )
    ALLOCATE( derived_type% barycenters(n_matter_loc,3) )

    DO i_matter= 1, n_matter_loc, 1

      derived_type% masses(i_matter)= 0.0D0
      derived_type% centers(i_matter,:)= 0.0D0
      derived_type% barycenters(i_matter,:)= 0.0D0
      derived_type% sizes(i_matter,:)= [ ABS(MINVAL(grid_tmp( :, 1 ))), &
                                         ABS(MAXVAL(grid_tmp( :, 1 ))), &
                                         ABS(MINVAL(grid_tmp( :, 2 ))), &
                                         ABS(MAXVAL(grid_tmp( :, 2 ))), &
                                         ABS(MINVAL(grid_tmp( :, 3 ))), &
                                         ABS(MAXVAL(grid_tmp( :, 3 ))) ]

    ENDDO

    derived_type% xL_grid= MINVAL(grid_tmp( :, 1 ))
    derived_type% yL_grid= MINVAL(grid_tmp( :, 2 ))
    derived_type% zL_grid= MINVAL(grid_tmp( :, 3 ))

    ! Determine nx_grid, ny_grid, nz_grid TODO: parallelize everything

    ALLOCATE( grid_sorted( 2*derived_type% n_gridpoints ) )
    grid_sorted= 0.0D0

    ! Sort the x coordinates of the grid
    CALL indexx( 2*derived_type% n_gridpoints, grid_tmp( :, 1 ), grid_sorted )

    ! Count how many different x coordinates there are
    derived_type% nx_grid= 1
    DO itr= 1, 2*derived_type% n_gridpoints - 1, 1

      IF( grid_tmp(grid_sorted(itr),1) /= grid_tmp(grid_sorted(itr+1),1) )THEN
        derived_type% nx_grid= derived_type% nx_grid + 1
      ENDIF

    ENDDO

    ! Sort the y coordinates of the grid
    CALL indexx( 2*derived_type% n_gridpoints, grid_tmp( :, 2 ), grid_sorted )

    ! Count how many different y coordinates there are
    derived_type% ny_grid= 1
    DO itr= 1, 2*derived_type% n_gridpoints - 1, 1

      IF( grid_tmp(grid_sorted(itr),2) /= grid_tmp(grid_sorted(itr+1),2) )THEN
        derived_type% ny_grid= derived_type% ny_grid + 1
      ENDIF

    ENDDO

    ! Sort the z coordinates of the grid
    CALL indexx( 2*derived_type% n_gridpoints, grid_tmp( :, 3 ), grid_sorted )

    ! Count how many different z coordinates there are
    derived_type% nz_grid= 1
    DO itr= 1, derived_type% n_gridpoints - 1, 1

      IF( grid_tmp(grid_sorted(itr),3) /= grid_tmp(grid_sorted(itr+1),3) )THEN
        derived_type% nz_grid= derived_type% nz_grid + 1
      ENDIF

    ENDDO
    derived_type% nz_grid= 2*derived_type% nz_grid

    ! Check that the grid dimensions are consistent
    IF( derived_type% nx_grid*derived_type% ny_grid*derived_type% nz_grid &
        /= 2*derived_type% n_gridpoints )THEN

      PRINT *, derived_type% nx_grid
      PRINT *, derived_type% ny_grid
      PRINT *, derived_type% nz_grid
      PRINT *, derived_type% nx_grid*derived_type% ny_grid*derived_type% nz_grid
      PRINT *, 2*derived_type% n_gridpoints
      PRINT *, derived_type% n_gridpoints
      STOP

    ENDIF

    derived_type% dx_grid= (MAXVAL(grid_tmp( :, 1 ))-MINVAL(grid_tmp( :, 1 )))/(derived_type% nx_grid)
    derived_type% dy_grid= (MAXVAL(grid_tmp( :, 2 ))-MINVAL(grid_tmp( :, 2 )))/(derived_type% ny_grid)
    derived_type% dz_grid= (MAXVAL(grid_tmp( :, 3 ))-MINVAL(grid_tmp( :, 3 )))/(derived_type% nz_grid)






    derived_type% n_gridpoints=2.0D0*derived_type% n_gridpoints

    ! TODO: Add test that nx dx and the grid extent are consistent
    ztmp= derived_type% zL_grid
    DO k= 1, derived_type% nz_grid, 1
      ztmp= ztmp + derived_type% dz_grid
    ENDDO
    IF( ABS( ztmp - MAXVAL(grid_tmp( :, 3 )) ) &
        > derived_type% dz_grid/1.0D+6 )THEN
      PRINT *, "** ERROR! ztmp=", ztmp
      PRINT *, "          zR_grid=", MAXVAL(grid_tmp( :, 3 ))
      STOP
    ENDIF
    ytmp= derived_type% yL_grid
    DO j= 1, derived_type% ny_grid, 1
      ytmp= ytmp + derived_type% dy_grid
    ENDDO
    IF( ABS( ytmp - MAXVAL(grid_tmp( :, 2 )) ) &
        > derived_type% dz_grid/1.0D+6 )THEN
      PRINT *, "** ERROR! ytmp=", ytmp
      PRINT *, "          yR_grid=", MAXVAL(grid_tmp( :, 2 ))
      STOP
    ENDIF
    xtmp= derived_type% xL_grid
    DO i= 1, derived_type% nx_grid, 1
      xtmp= xtmp + derived_type% dx_grid
    ENDDO
    IF( ABS( xtmp - MAXVAL(grid_tmp( :, 1 )) ) &
        > derived_type% dx_grid/1.0D+6 )THEN
      PRINT *, "** ERROR! xtmp=", xtmp
      PRINT *, "          xR_grid=", MAXVAL(grid_tmp( :, 1 ))
      STOP
    ENDIF

    !derived_type% nz_grid= 2.0D0*derived_type% nz_grid

    ! Store the grid on the member array
    ALLOCATE( derived_type% grid(derived_type% nx_grid, &
                                 derived_type% ny_grid, &
                                 derived_type% nz_grid, 3) )
  !
  !  DO k= 1, derived_type% nz_grid, 1
  !    DO j= 1, derived_type% ny_grid, 1
  !      DO i= 1, derived_type% nx_grid, 1
  !
  !        derived_type% grid(i,j,k,:)= &
  !              grid_tmp((k-1)*(derived_type% ny_grid)*(derived_type% nx_grid) &
  !                     + (j-1)*(derived_type% nx_grid) + i, :)
  !
  !      ENDDO
  !    ENDDO
  !  ENDDO

    ! Read the ID
    OPEN( UNIT= unit_pos, FILE= TRIM(filename), &
          FORM= "FORMATTED", ACTION= "READ" )

    ! Skip header
    DO itr= 1, header_lines, 1
      READ( unit_pos, * )
    ENDDO

    ! Allocate the arrays to store data
    ALLOCATE( derived_type% baryon_mass_density( derived_type% nx_grid, &
                                                 derived_type% ny_grid, &
                                                 derived_type% nz_grid ) )
    derived_type% baryon_mass_density= 0.0D0

    ! Read the data into the array
    DO i= 1, derived_type% nx_grid, 1
      DO j= 1, derived_type% ny_grid, 1
        DO k= 1, derived_type% nz_grid/2.0D0, 1

          READ( UNIT= unit_pos, FMT= *, IOSTAT = ios, IOMSG= err_msg ) &
            derived_type% grid( i, j, k, 1 ), &
            derived_type% grid( i, j, k, 2 ), &
            derived_type% grid( i, j, k, 3 ), &
            derived_type% baryon_mass_density( i, j, k )

          IF( ios > 0 )THEN
            PRINT *, "...error when reading " // TRIM(filename), &
                    " at particle ", itr,". The status variable is ", ios, &
                    ". The error message is", err_msg
            STOP
          ENDIF

        ENDDO
      ENDDO
    ENDDO
    DO i= 1, derived_type% nx_grid, 1
      DO j= 1, derived_type% ny_grid, 1
        DO k= 1, derived_type% nz_grid/2.0D0, 1

          derived_type% grid( i, j, derived_type% nz_grid/2.0D0 + k, 1 )= &
                          derived_type% grid( i, j, k, 1 )
          derived_type% grid( i, j, derived_type% nz_grid/2.0D0 + k, 2 )= &
                          derived_type% grid( i, j, k, 2 )
          derived_type% grid( i, j, derived_type% nz_grid/2.0D0 + k, 3 )= &
                        - derived_type% grid( i, j, k, 3 )
          derived_type% baryon_mass_density( i, j, derived_type% nz_grid/2.0D0 + k )= &
                            derived_type% baryon_mass_density( i, j, k )

        ENDDO
      ENDDO
    ENDDO

    CLOSE( UNIT= unit_pos )

    !derived_type% zL_grid= MINVAL(derived_type% grid( 1, 1, :, 3 ))
    !ztmp= derived_type% zL_grid
    !DO k= 1, derived_type% nz_grid, 1
    !  ztmp= ztmp + derived_type% dz_grid
    !ENDDO
    !IF( ABS( ztmp - MAXVAL(derived_type% grid( 1, 1, :, 3 )) ) &
    !    > derived_type% dz_grid/1.0D+6 )THEN
    !  PRINT *, "** ERROR! ztmp=", ztmp
    !  PRINT *, "          zR_grid=", MAXVAL(derived_type% grid( 1, 1, :, 3 ))
    !  STOP
    !ENDIF

    PRINT *, derived_type% nx_grid
    PRINT *, derived_type% ny_grid
    PRINT *, derived_type% nz_grid
    PRINT *
    PRINT *, SIZE( derived_type% grid( :, 1, 1, 1 ) )
    PRINT *, SIZE( derived_type% grid( 1, :, 1, 1 ) )
    PRINT *, SIZE( derived_type% grid( 1, 1, :, 1 ) )
    PRINT *
    PRINT *, MAXVAL( derived_type% grid( :, 1, 1, 1 ) )
    PRINT *, MAXVAL( derived_type% grid( 1, :, 1, 2 ) )
    PRINT *, MAXVAL( derived_type% grid( 1, 1, :, 3 ) )
    PRINT *
    PRINT *, MINVAL( derived_type% grid( :, 1, 1, 1 ) )
    PRINT *, MINVAL( derived_type% grid( 1, :, 1, 2 ) )
    PRINT *, MINVAL( derived_type% grid( 1, 1, :, 3 ) )
    PRINT *
    PRINT *, derived_type% dx_grid
    PRINT *, derived_type% dy_grid
    PRINT *, derived_type% dz_grid
    PRINT *

    ALLOCATE( rho_tmp(derived_type% nx_grid, derived_type% ny_grid, derived_type% nz_grid) )
    rho_tmp= 0.0D0

!$OMP PARALLEL DO SHARED(derived_type,grid_tmp,rho_tmp) &
!$OMP             PRIVATE(i,j,k)
DO k= 1, derived_type% nz_grid, 1
  DO j= 1, derived_type% ny_grid, 1
    DO i= 1, derived_type% nx_grid, 1

IF( grid_tmp((k-1)*(derived_type% ny_grid)*(derived_type% nx_grid) &
    + (j-1)*(derived_type% nx_grid) + i,3) < MAXVAL(grid_tmp(:,3)) - derived_type% dz_grid &
    .AND. &
    grid_tmp((k-1)*(derived_type% ny_grid)*(derived_type% nx_grid) &
            + (j-1)*(derived_type% nx_grid) + i,2) < MAXVAL(grid_tmp(:,2)) - derived_type% dy_grid &
    .AND. &
    grid_tmp((k-1)*(derived_type% ny_grid)*(derived_type% nx_grid) &
    + (j-1)*(derived_type% nx_grid) + i,1) < MAXVAL(grid_tmp(:,1)) - derived_type% dx_grid &
    .AND. &
    grid_tmp((k-1)*(derived_type% ny_grid)*(derived_type% nx_grid) &
            + (j-1)*(derived_type% nx_grid) + i,3) > MINVAL(grid_tmp(:,3)) + 2.0D0*derived_type% dz_grid &
    .AND. &
    grid_tmp((k-1)*(derived_type% ny_grid)*(derived_type% nx_grid) &
            + (j-1)*(derived_type% nx_grid) + i,2) > MINVAL(grid_tmp(:,2)) + 2.0D0*derived_type% dy_grid &
    .AND. &
    grid_tmp((k-1)*(derived_type% ny_grid)*(derived_type% nx_grid) &
            + (j-1)*(derived_type% nx_grid) + i,1) > MINVAL(grid_tmp(:,1)) + 2.0D0*derived_type% dx_grid )THEN

       rho_tmp(i,j,k)= derived_type% read_mass_density( &
          grid_tmp((k-1)*(derived_type% ny_grid)*(derived_type% nx_grid) &
                             + (j-1)*(derived_type% nx_grid) + i,1), &
          grid_tmp((k-1)*(derived_type% ny_grid)*(derived_type% nx_grid) &
                             + (j-1)*(derived_type% nx_grid) + i,2), &
          grid_tmp((k-1)*(derived_type% ny_grid)*(derived_type% nx_grid) &
                             + (j-1)*(derived_type% nx_grid) + i,3) &
        )!, &
      ENDIF

ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO

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

    DO k= 10, derived_type% nz_grid - 10, 1
      DO j= 10, derived_type% ny_grid - 10, 1
        DO i= 10, derived_type% nx_grid - 10, 1

          WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
          !  derived_type% grid( i, j, k, 1 ), &
          !  derived_type% grid( i, j, k, 2 ), &
          !  derived_type% grid( i, j, k, 3 ), &
          !  derived_type% baryon_mass_density( i, j, k ), &
          !  derived_type% read_mass_density( &
          !    derived_type% grid( i, j, k, 1 ), &
          !    derived_type% grid( i, j, k, 2 ), &
          !    derived_type% grid( i, j, k, 3 ) &
          !  )!, &
            grid_tmp((k-1)*(derived_type% ny_grid)*(derived_type% nx_grid) &
                                 + (j-1)*(derived_type% nx_grid) + i,1), &
            grid_tmp((k-1)*(derived_type% ny_grid)*(derived_type% nx_grid) &
                                 + (j-1)*(derived_type% nx_grid) + i,2), &
            grid_tmp((k-1)*(derived_type% ny_grid)*(derived_type% nx_grid) &
                                 + (j-1)*(derived_type% nx_grid) + i,3), &

            derived_type% baryon_mass_density( i, j, k ), &
            rho_tmp(i,j,k)
            !derived_type% read_mass_density( &
            !  grid_tmp((k-1)*(derived_type% ny_grid)*(derived_type% nx_grid) &
            !                     + (j-1)*(derived_type% nx_grid) + i,1), &
            !  grid_tmp((k-1)*(derived_type% ny_grid)*(derived_type% nx_grid) &
            !                     + (j-1)*(derived_type% nx_grid) + i,2), &
            !  grid_tmp((k-1)*(derived_type% ny_grid)*(derived_type% nx_grid) &
            !                     + (j-1)*(derived_type% nx_grid) + i,3) &
            !)!, &
           ! derived_type% read_mass_density( &
           !   derived_type% grid( i, j, k, 1 ), &
           !   derived_type% grid( i, j, k, 2 ), &
           !   derived_type% grid( i, j, k, 3 ) + derived_type% dz_grid/2.0D0 &
           ! ), &
           ! derived_type% read_mass_density( &
           !   derived_type% grid( i, j, k, 1 ), &
           !   derived_type% grid( i, j, k, 2 ) + derived_type% dy_grid/2.0D0, &
           !   derived_type% grid( i, j, k, 3 ) &
           ! ), &
           ! derived_type% read_mass_density( &
           !   derived_type% grid( i, j, k, 1 ) + derived_type% dx_grid/2.0D0, &
           !   derived_type% grid( i, j, k, 2 ) + derived_type% dy_grid/2.0D0, &
           !   derived_type% grid( i, j, k, 3 ) &
           ! ), &
           ! derived_type% read_mass_density( &
           !   derived_type% grid( i, j, k, 1 ) + derived_type% dx_grid/2.0D0, &
           !   derived_type% grid( i, j, k, 2 ), &
           !   derived_type% grid( i, j, k, 3 ) + derived_type% dz_grid/2.0D0 &
           ! ), &
           ! derived_type% read_mass_density( &
           !   derived_type% grid( i, j, k, 1 ), &
           !   derived_type% grid( i, j, k, 2 ) + derived_type% dy_grid/2.0D0, &
           !   derived_type% grid( i, j, k, 3 ) + derived_type% dz_grid/2.0D0 &
           ! ), &
           ! derived_type% read_mass_density( &
           !   derived_type% grid( i, j, k, 1 ) + derived_type% dx_grid/2.0D0, &
           !   derived_type% grid( i, j, k, 2 ) + derived_type% dy_grid/2.0D0, &
           !   derived_type% grid( i, j, k, 3 ) + derived_type% dz_grid/2.0D0 &
           ! )

        ENDDO
      ENDDO
    ENDDO

    CLOSE( UNIT= 2 )

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
