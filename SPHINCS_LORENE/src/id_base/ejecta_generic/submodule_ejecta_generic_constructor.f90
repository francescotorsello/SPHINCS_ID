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
    INTEGER:: nlines
    INTEGER:: i_matter, n_matter_loc, itr, i, j, k

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: grid_tmp
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

    ! Read the ID
    OPEN( UNIT= unit_pos, FILE= TRIM(filename), &
          FORM= "FORMATTED", ACTION= "READ" )

    ! Skip header
    DO itr= 1, header_lines, 1
      READ( unit_pos, * )
    ENDDO

    ! Allocate the temporary array to store data
    ALLOCATE( grid_tmp( derived_type% n_gridpoints, 3 ) )
    grid_tmp= 0.0D0

    ! Read the data into the temporary array
    DO itr= 1, derived_type% n_gridpoints, 1

      READ( UNIT= unit_pos, FMT= *, IOSTAT = ios, IOMSG= err_msg ) &
        grid_tmp( itr, : )

      IF( ios > 0 )THEN
        PRINT *, "...error when reading " // TRIM(filename), &
                " at particle ", itr,". The status variable is ", ios, &
                ". The error message is", err_msg
        STOP
      ENDIF

    ENDDO

    CLOSE( UNIT= unit_pos )

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

    ALLOCATE( grid_sorted( derived_type% n_gridpoints ) )
    grid_sorted= 0.0D0

    ! Sort the x coordinates of the grid
    CALL indexx( derived_type% n_gridpoints, grid_tmp( :, 1 ), grid_sorted )

    ! Count how many different x coordinates there are
    derived_type% nx_grid= 1
    DO itr= 1, derived_type% n_gridpoints - 1, 1

      IF( grid_tmp(grid_sorted(itr),1) /= grid_tmp(grid_sorted(itr+1),1) )THEN
        derived_type% nx_grid= derived_type% nx_grid + 1
      ENDIF

    ENDDO

    ! Sort the y coordinates of the grid
    CALL indexx( derived_type% n_gridpoints, grid_tmp( :, 2 ), grid_sorted )

    ! Count how many different y coordinates there are
    derived_type% ny_grid= 1
    DO itr= 1, derived_type% n_gridpoints - 1, 1

      IF( grid_tmp(grid_sorted(itr),2) /= grid_tmp(grid_sorted(itr+1),2) )THEN
        derived_type% ny_grid= derived_type% ny_grid + 1
      ENDIF

    ENDDO

    ! Sort the z coordinates of the grid
    CALL indexx( derived_type% n_gridpoints, grid_tmp( :, 3 ), grid_sorted )

    ! Count how many different z coordinates there are
    derived_type% nz_grid= 1
    DO itr= 1, derived_type% n_gridpoints - 1, 1

      IF( grid_tmp(grid_sorted(itr),3) /= grid_tmp(grid_sorted(itr+1),3) )THEN
        derived_type% nz_grid= derived_type% nz_grid + 1
      ENDIF

    ENDDO

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

    derived_type% dx_grid= (MAXVAL(grid_tmp( :, 1 ))-MINVAL(grid_tmp( :, 1 )))/(derived_type% nx_grid-1)
    derived_type% dy_grid= (MAXVAL(grid_tmp( :, 2 ))-MINVAL(grid_tmp( :, 2 )))/(derived_type% ny_grid-1)
    derived_type% dz_grid= (MAXVAL(grid_tmp( :, 3 ))-MINVAL(grid_tmp( :, 3 )))/(derived_type% nz_grid-1)

    ! TODO: Add test that nx dx and the grid extent are consistent

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
    ALLOCATE( derived_type% baryon_mass_density(derived_type% nx_grid, &
                                                derived_type% ny_grid, &
                                                derived_type% nz_grid ) )
    grid_tmp= 0.0D0

    ! Read the data into the array
    DO k= 1, derived_type% nz_grid, 1
      DO j= 1, derived_type% ny_grid, 1
        DO i= 1, derived_type% nx_grid, 1

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

    CLOSE( UNIT= unit_pos )

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

    DO k= 1, derived_type% nz_grid, 1
      DO j= 1, derived_type% ny_grid, 1
        DO i= 1, derived_type% nx_grid, 1

          WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
            derived_type% grid( i, j, k, 1 ), &
            derived_type% grid( i, j, k, 2 ), &
            derived_type% grid( i, j, k, 3 ), &
            derived_type% baryon_mass_density( i, j, k )

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
