! File:         submodule_bnsfuka_constructor.f90
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

SUBMODULE (bns_fuka) constructor

  !*********************************************************
  !
  !# Implementation of the constructor and
  !  destructor of TYPE [[bnsfuka]], and of the
  !  [[bnsfuka]]-member
  !  PROCEDURES that call the C-bound PROCEDURES
  !  constructig and destructing the |fuka|
  !  |binns| object
  !
  !  FT 23.10.2020
  !
  !*********************************************************


  IMPLICIT NONE


  CONTAINS


  !
  !-- Implementation of the constructor of the bns object
  !
  MODULE PROCEDURE construct_bnsfuka

    !****************************************************
    !
    !# Constructs an object of TYPE [[bnsfuka]]
    !
    !  Created:     FT 09.02.2022
    !  Last update: FT 27.05.2022
    !
    !****************************************************

    USE utility,  ONLY: ten, Msun_geo

    IMPLICIT NONE

    INTEGER, SAVE:: bns_counter= 1

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: length_scale_pressure

    CALL derived_type% set_n_matter(2)
    CALL derived_type% set_cold_system(.TRUE.)

    derived_type% construction_timer= timer( "binary_construction_timer" )

    ! Construct |fuka| bns_export object
    CALL derived_type% construct_binary( filename )

    ! Read the parameters of the binary system
    CALL read_fuka_id_params( derived_type )

    ! Assign a unique identifier to the bns_fuka object
    derived_type% bns_identifier= bns_counter
    bns_counter= bns_counter + 1

    ! Do not use the geodesic gauge by default
    CALL derived_type% set_one_lapse ( .FALSE. )
    CALL derived_type% set_zero_shift( .FALSE. )

    ! Since Kadath is not thread-safe, we cannot parallelize it using OMP
    ! within SPHINCS_ID. Hence, we chose to make a system call to a program
    ! within Kadath that reads the ID from the FUKA output file and prints it
    ! on a lattice. The ID on the particles will be interplated from ths fine
    ! lattice.
    CALL set_up_lattices_around_stars()

    ! Compute typical length scales of the system using the pressure
    IF( derived_type% get_estimate_length_scale() )THEN

      ALLOCATE( length_scale_pressure(derived_type% get_n_matter()) )
      length_scale_pressure= derived_type% estimate_lengthscale_field( &
                                                get_pressure, &
                                                derived_type% get_n_matter() )

      PRINT *, " * Minimum length scale to resolve on star 1, based on ", &
               "pressure= ", length_scale_pressure(1)*Msun_geo*ten*ten*ten, "m"
      PRINT *, " * Minimum length scale to resolve on star 2, based on ", &
               "pressure= ", length_scale_pressure(1)*Msun_geo*ten*ten*ten, "m"
      PRINT *

    ENDIF

    ! Assign PROCEDURE POINTER to the desired PROCEDURE
    derived_type% finalize_sph_id_ptr => finalize

    CONTAINS

    FUNCTION get_pressure( x, y, z ) RESULT( val )
    !! Returns the value of the pressure at the desired point

      DOUBLE PRECISION, INTENT(IN):: x
      !! \(x\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN):: y
      !! \(y\) coordinate of the desired point
      DOUBLE PRECISION, INTENT(IN):: z
      !! \(z\) coordinate of the desired point
      DOUBLE PRECISION:: val
      !! Pressure at \((x,y,z)\)

      val= derived_type% read_fuka_pressure( x, y, z )

    END FUNCTION get_pressure


    SUBROUTINE set_up_lattices_around_stars()

      !***********************************************
      !
      !#
      !
      !  FT 27.06.2022
      !
      !***********************************************

      !USE OMP_LIB, ONLY: OMP_GET_NUM_PROCS

#ifdef __INTEL_COMPILER

  USE IFPORT, ONLY: CHANGEDIRQQ

#endif

      IMPLICIT NONE

      INTEGER, PARAMETER:: nx       = 25
      INTEGER, PARAMETER:: ny       = 25
      INTEGER, PARAMETER:: nz       = 25
      INTEGER, PARAMETER:: unit_par = 3480
      INTEGER, PARAMETER:: mpi_ranks= 40

      DOUBLE PRECISION, PARAMETER:: stretch= 1.05D0
      !! The lattices' sizes will be 5% larger than the radii of the stars

      INTEGER:: i_star, ios, i_char, i_file, i, j, k
      INTEGER:: nlines, nlines_prev
      INTEGER, DIMENSION(mpi_ranks):: unit_rank
      INTEGER:: unit_rank_prev

      DOUBLE PRECISION:: xmin, xmax, ymin, ymax, zmin, zmax
      DOUBLE PRECISION, DIMENSION(6):: sizes
      DOUBLE PRECISION, DIMENSION(3):: center
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: grid_tmp

      LOGICAL:: exist
      LOGICAL(4):: status

      CHARACTER( LEN= : ), ALLOCATABLE:: filename_par
      CHARACTER( LEN= : ), ALLOCATABLE:: filename_id
      CHARACTER( LEN= : ), ALLOCATABLE:: filename_rank
      CHARACTER( LEN= : ), ALLOCATABLE:: dir_id
      CHARACTER( LEN= : ), ALLOCATABLE:: working_directory
      CHARACTER( LEN= 3 ):: mpi_ranks_str

      !num_processors= OMP_GET_NUM_PROCS()
      !PRINT*, num_processors
      !STOP

      ALLOCATE( derived_type% id_fields( nx, ny, nz, n_fields_fuka, 2 ) )

      loop_over_stars: DO i_star= 1, 2, 1

        sizes = derived_type% return_spatial_extent(i_star)
        center= derived_type% return_center(i_star)

        ! Determine boundaries of the lattices
        xmin= center(1) - stretch*sizes(1)
        xmax= center(1) + stretch*sizes(2)
        ymin= center(2) - stretch*sizes(3)
        ymax= center(2) + stretch*sizes(4)
        zmin= center(3) - stretch*sizes(5)
        zmax= center(3) + stretch*sizes(6)

        filename_par= "lattice_par.dat"

        find_name_loop: DO i_char= LEN(filename), 1, -1

          IF( filename(i_char:i_char) == "/" )THEN
            filename_id= TRIM(filename(i_char+1:LEN(filename)))
            dir_id     = TRIM(filename(1:i_char))
            EXIT
          ENDIF

        ENDDO find_name_loop

        INQUIRE( FILE= TRIM(dir_id//filename_par), EXIST= exist )

        IF( exist )THEN
            OPEN( UNIT= unit_par, FILE= TRIM(dir_id//filename_par), &
                  STATUS= "REPLACE", FORM= "FORMATTED", &
                  POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
                  IOMSG= err_msg )
        ELSE
            OPEN( UNIT= unit_par, FILE= TRIM(dir_id//filename_par), &
                  STATUS= "NEW", FORM= "FORMATTED", &
                  ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
        ENDIF
        IF( ios > 0 )THEN
          PRINT *, "...error when opening " // TRIM(dir_id//filename_par), &
                   ". The error message is", err_msg
          STOP
        ENDIF

        WRITE( UNIT = unit_par, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
          TRIM(filename_id)
        IF( ios > 0 )THEN
          PRINT *, "...error when writing the arrays in ", &
                   TRIM(dir_id//filename_par), ". The error message is", err_msg
          STOP
        ENDIF
        WRITE( UNIT = unit_par, IOSTAT = ios, IOMSG = err_msg, FMT = * ) xmin
        IF( ios > 0 )THEN
          PRINT *, "...error when writing the arrays in ", &
                   TRIM(dir_id//filename_par), ". The error message is", err_msg
          STOP
        ENDIF
        WRITE( UNIT = unit_par, IOSTAT = ios, IOMSG = err_msg, FMT = * ) xmax
        IF( ios > 0 )THEN
          PRINT *, "...error when writing the arrays in ", &
                   TRIM(dir_id//filename_par), ". The error message is", err_msg
          STOP
        ENDIF
        WRITE( UNIT = unit_par, IOSTAT = ios, IOMSG = err_msg, FMT = * ) ymin
        IF( ios > 0 )THEN
          PRINT *, "...error when writing the arrays in ", &
                   TRIM(dir_id//filename_par), ". The error message is", err_msg
          STOP
        ENDIF
        WRITE( UNIT = unit_par, IOSTAT = ios, IOMSG = err_msg, FMT = * ) ymax
        IF( ios > 0 )THEN
          PRINT *, "...error when writing the arrays in ", &
                   TRIM(dir_id//filename_par), ". The error message is", err_msg
          STOP
        ENDIF
        WRITE( UNIT = unit_par, IOSTAT = ios, IOMSG = err_msg, FMT = * ) zmin
        IF( ios > 0 )THEN
          PRINT *, "...error when writing the arrays in ", &
                   TRIM(dir_id//filename_par), ". The error message is", err_msg
          STOP
        ENDIF
        WRITE( UNIT = unit_par, IOSTAT = ios, IOMSG = err_msg, FMT = * ) zmax
        IF( ios > 0 )THEN
          PRINT *, "...error when writing the arrays in ", &
                   TRIM(dir_id//filename_par), ". The error message is", err_msg
          STOP
        ENDIF
        WRITE( UNIT = unit_par, IOSTAT = ios, IOMSG = err_msg, FMT = * ) nx
        IF( ios > 0 )THEN
          PRINT *, "...error when writing the arrays in ", &
                   TRIM(dir_id//filename_par), ". The error message is", err_msg
          STOP
        ENDIF
        WRITE( UNIT = unit_par, IOSTAT = ios, IOMSG = err_msg, FMT = * ) ny
        IF( ios > 0 )THEN
          PRINT *, "...error when writing the arrays in ", &
                   TRIM(dir_id//filename_par), ". The error message is", err_msg
          STOP
        ENDIF
        WRITE( UNIT = unit_par, IOSTAT = ios, IOMSG = err_msg, FMT = * ) nz
        IF( ios > 0 )THEN
          PRINT *, "...error when writing the arrays in ", &
                   TRIM(dir_id//filename_par), ". The error message is", err_msg
          STOP
        ENDIF

        CLOSE( UNIT= unit_par )

        CALL derived_type% run_kadath_reader( mpi_ranks, nx, ny, nz, &
                                    derived_type% id_fields(:,:,:,:,i_star), &
                                    filename_par, dir_id )

        PRINT *, "** ID stored on a fine lattice around star ", i_star
        PRINT *

      ENDDO loop_over_stars

    !  filename_id= "dbg-id.dat"
    !
    !  INQUIRE( FILE= TRIM(filename_id), EXIST= exist )
    !
    !  IF( exist )THEN
    !    OPEN( UNIT= 2, FILE= TRIM(filename_id), STATUS= "REPLACE", &
    !          FORM= "FORMATTED", &
    !          POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
    !          IOMSG= err_msg )
    !  ELSE
    !    OPEN( UNIT= 2, FILE= TRIM(filename_id), STATUS= "NEW", &
    !          FORM= "FORMATTED", &
    !          ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    !  ENDIF
    !  IF( ios > 0 )THEN
    !  PRINT *, "...error when opening " // TRIM(filename_id), &
    !           ". The error message is", err_msg
    !  STOP
    !  ENDIF
    !
    !  DO i_star= 1, 2, 1
    !    DO k= 1, nz, 1
    !      DO j= 1, ny, 1
    !        DO i= 1, nx, 1
    !
    !          WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    !            derived_type% id_fields( i, j, k, :, i_star )
    !
    !        ENDDO
    !      ENDDO
    !    ENDDO
    !  ENDDO
    !
    !  CLOSE( UNIT= 2 )

      STOP

    END SUBROUTINE set_up_lattices_around_stars


  END PROCEDURE construct_bnsfuka




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


  !
  !-- Implementation of the destructor of the bns object
  !
  MODULE PROCEDURE destruct_bnsfuka

    !***********************************************
    !
    !# Destructs an object of TYPE [[bnsfuka]]
    !
    !  Created:     FT 09.02.2022
    !  Last update: FT 27.05.2022
    !
    !***********************************************

    IMPLICIT NONE

    ! Deallocate memory
    CALL this% deallocate_bnsfuka_memory()

  END PROCEDURE destruct_bnsfuka


  MODULE PROCEDURE construct_binary

    !***********************************************
    !
    !# Construct the |fuka| bns_export object
    !
    !  Created:     FT 09.02.2022
    !  Last update: FT 27.05.2022
    !
    !***********************************************

    IMPLICIT NONE

    LOGICAL:: exist

#ifdef __INTEL_COMPILER

    IF ( C_ASSOCIATED( this% bns_ptr ) ) THEN

      CALL destruct_bns_fuka( this% bns_ptr )

    ENDIF

#endif

    INQUIRE( FILE= fukafile, EXIST= exist )

    IF( exist )THEN

      CALL this% construction_timer% start_timer()
      this% bns_ptr = construct_bns_fuka( fukafile//C_NULL_CHAR )
      CALL this% construction_timer% stop_timer()

    ELSE

      PRINT *, "** ERROR in SUBROUTINE construct_binary: file ", &
               fukafile, " cannot be found!"
      PRINT *
      STOP

    ENDIF


  END PROCEDURE construct_binary


  MODULE PROCEDURE destruct_binary

    !************************************************
    !
    !# Destructs the |fuka| bns_export object and frees
    !  the pointer [[bns:bns_ptr]] pointing to it
    !
    !  Created:     FT 09.02.2022
    !  Last update: FT 27.05.2022
    !
    !************************************************

    IMPLICIT NONE


    IF ( C_ASSOCIATED( this% bns_ptr ) ) THEN

      CALL destruct_bns_fuka( this% bns_ptr )
      this% bns_ptr = C_NULL_PTR

    ENDIF

  END PROCEDURE destruct_binary


END SUBMODULE constructor
