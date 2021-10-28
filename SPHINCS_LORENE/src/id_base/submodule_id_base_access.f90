! File:         submodule_id_base_access.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (id_base) id_base_access

  !********************************************
  !
  !# Implementation of the methods of TYPE
  !  idbase to access PRIVATE data
  !
  !  FT 28.10.2021
  !
  !********************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE set_n_matter

    !************************************************
    !
    !# Sets [[idbase:n_matter]] to the given value
    !
    !  FT 28.10.2021
    !
    !************************************************

    IMPLICIT NONE

    THIS% n_matter= value

  END PROCEDURE set_n_matter


  MODULE PROCEDURE get_n_matter

    !************************************************
    !
    !# Returns [[idbase:n_matter]], the number of
    !  matter objects in the physical system
    !
    !  FT 28.10.2021
    !
    !************************************************

    IMPLICIT NONE

    get_n_matter= THIS% n_matter

  END PROCEDURE get_n_matter


  MODULE PROCEDURE get_total_spatial_extent

    !************************************************
    !
    !# Return the coordinates
    !\(x_{\rm min},x_{\rm max},y_{\rm min},y_{\rm max},z_{\rm min},z_{\rm max}\)
    !  of a box containing the entire physical system
    !
    !  FT 28.10.2021
    !
    !************************************************

    IMPLICIT NONE

    INTEGER:: i_matter
    DOUBLE PRECISION, DIMENSION(3):: center_matter
    DOUBLE PRECISION, DIMENSION(6):: size_matter

    box(1)=   HUGE(1.0D0)
    box(3)=   HUGE(1.0D0)
    box(5)=   HUGE(1.0D0)
    box(2)= - HUGE(1.0D0)
    box(4)= - HUGE(1.0D0)
    box(6)= - HUGE(1.0D0)

    DO i_matter= 1, THIS% get_n_matter(), 1

      size_matter  = THIS% return_spatial_extent( i_matter )
      center_matter= THIS% return_center( i_matter )

      IF( size_matter(1) < box(1) ) box(1) = center_matter(1) - size_matter(1)
      IF( size_matter(2) > box(2) ) box(2) = center_matter(1) + size_matter(2)
      IF( size_matter(3) < box(3) ) box(3) = center_matter(2) - size_matter(3)
      IF( size_matter(4) > box(4) ) box(4) = center_matter(2) + size_matter(4)
      IF( size_matter(5) < box(5) ) box(5) = center_matter(3) - size_matter(5)
      IF( size_matter(6) > box(6) ) box(6) = center_matter(3) + size_matter(6)

    ENDDO

  END PROCEDURE get_total_spatial_extent



END SUBMODULE id_base_access
