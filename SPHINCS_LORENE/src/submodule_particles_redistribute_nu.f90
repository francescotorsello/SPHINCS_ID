! File:         submodule_particles_redistribute_nu.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (particles_id) particles_redistribute_nu

  !***************************************************
  !
  !# This SUBMODULE contains the implementation of
  !  the methods of TYPE particles
  !  that reallocate the sph variables and
  !  redistribute nu (baryon number per particle)
  !  on the particles.
  !
  !  These methods find application when one wants to
  !  decrease the particle mass ratio with particles
  !  on lattices.
  !
  !  They DON'T HAVE ANYTHING to do with the APM.
  !
  !  TODO: Add other SUBROUTINES to improve
  !        modularity in SUBMODULE
  !        particles_sph_variables
  !
  !  FT 12.07.2021
  !
  !***************************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE reshape_sph_field_1d

    !************************************************
    !
    !# Read the SPH ID from the binary file output
    !  by write_SPHINCS_dump, and print it to a
    !  formatted file
    !
    !  FT 31.03.2021
    !
    !************************************************

    IMPLICIT NONE

    INTEGER:: itr, i_tmp
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: tmp

    ALLOCATE( tmp( new_size1 + new_size2 ), STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
      PRINT *, "...allocation error for array tmp in SUBROUTINE ", &
               "reshape_sph_field_1d. The error message is", err_msg
      STOP
    ENDIF
    i_tmp= 0
    DO itr= THIS% npart1, THIS% npart1 - new_size1 + 1, -1

      i_tmp= i_tmp + 1
      tmp( i_tmp )= field( index_array( itr ) )
      !IF( itr == THIS% npart1 - new_size1 + 1 )THEN
      !  PRINT *, i_tmp
      !ENDIF

    ENDDO
    DO itr= THIS% npart, THIS% npart - new_size2 + 1, -1

      i_tmp= i_tmp + 1
      tmp( i_tmp )= field( index_array( itr ) )
      !IF( itr == THIS% npart )THEN
      !  PRINT *, i_tmp
      !ENDIF
      !IF( itr == THIS% npart - new_size2 + 1 )THEN
      !  PRINT *, i_tmp
      !  PRINT *, new_size1 + new_size2
      !ENDIF

    ENDDO

    DEALLOCATE( field, STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
      PRINT *, "...deallocation error for array field in SUBROUTINE ", &
               "reshape_sph_field_1d. The error message is", err_msg
      STOP
    ENDIF
    ALLOCATE( field( new_size1 + new_size2 ), STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
      PRINT *, "...allocation error for array field in SUBROUTINE ", &
               "reshape_sph_field_1d. The error message is", err_msg
      STOP
    ENDIF

    DO itr= 1, new_size1 + new_size2, 1
      field( itr )= tmp( itr )
    ENDDO

  END PROCEDURE reshape_sph_field_1d


  MODULE PROCEDURE reshape_sph_field_2d

    !************************************************
    !
    !# Read the SPH ID from the binary file output
    !  by write_SPHINCS_dump, and print it to a
    !  formatted file
    !
    !  FT 31.03.2021
    !
    !************************************************

    IMPLICIT NONE

    INTEGER:: itr, i_tmp, itr2
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: tmp

    ALLOCATE( tmp( 3, new_size1 + new_size2 ), STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
      PRINT *, "...allocation error for array tmp in SUBROUTINE ", &
               "reshape_sph_field_2d. The error message is", err_msg
      STOP
    ENDIF
    DO itr2= 1, 3, 1
      i_tmp= 0
      DO itr= THIS% npart1, THIS% npart1 - new_size1 + 1, -1

        i_tmp= i_tmp + 1
        tmp( itr2, i_tmp )= field( itr2, index_array( itr ) )
        !PRINT *, field( itr2, index_array( itr ) )
        !PRINT *, index_array( itr )
        !PRINT *, itr
        !EXIT

      ENDDO
      !PRINT *, i_tmp
      !PRINT *, new_size1
      DO itr= THIS% npart, THIS% npart - new_size2 + 1, -1

        i_tmp= i_tmp + 1
        tmp( itr2, i_tmp )= field( itr2, index_array( itr ) )
        !PRINT *, field( itr2, index_array( itr ) )
        !PRINT *, index_array( itr )
        !PRINT *, itr
        !IF( field( itr2, index_array( itr ) ) < 0 )THEN
        !  PRINT *, "The x coordinate of the second star is negative...", &
        !           "something is wrong"
        !  STOP
        !ENDIF
        !STOP

      ENDDO
      !PRINT *, i_tmp - new_size1
      !PRINT *, new_size2
    ENDDO

    DEALLOCATE( field, STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
      PRINT *, "...deallocation error for array field in SUBROUTINE ", &
               "reshape_sph_field_2d. The error message is", err_msg
      STOP
    ENDIF

    ALLOCATE( field( 3, new_size1 + new_size2 ), STAT= ios, ERRMSG= err_msg )
    IF( ios > 0 )THEN
      PRINT *, "...allocation error for array field in SUBROUTINE ", &
               "reshape_sph_field_2d. The error message is", err_msg
      STOP
    ENDIF

    DO itr2= 1, 3, 1
      DO itr= 1, new_size1 + new_size2, 1
        field( itr2, itr )= tmp( itr2, itr )
      ENDDO
    ENDDO

  END PROCEDURE reshape_sph_field_2d


END SUBMODULE particles_redistribute_nu
