! File:         submodule_bns_base_geometry.f90
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

SUBMODULE(bns_base) geometry

  !********************************************************
  !
  !# This MODULE contains the implementation of the 
  !  methods of TYPE bns_base that deal with the
  !  computation of the radii of the stars.
  !
  !  FT 27.09.2022
  !
  !********************************************************


  USE utility, ONLY: zero, one, two, three, four, five, ten


  IMPLICIT NONE


  DOUBLE PRECISION, PARAMETER:: tol= 1.D-5



  CONTAINS



  MODULE PROCEDURE find_print_surfaces

    !********************************************************
    !
    !# Finds the surfaces of the stars, and prints them to
    !  a formatted file.
    !
    !  FT 18.02.2022
    !
    !********************************************************

    IMPLICIT NONE

    INTEGER, PARAMETER:: n_theta  = 180
    INTEGER, PARAMETER:: n_phi    = 360
    INTEGER, PARAMETER:: unit_dump= 2764

    INTEGER:: i_matter, n_matter, i, j, ios
    LOGICAL:: exist
    CHARACTER(LEN=3):: str_i
    CHARACTER(LEN=:), ALLOCATABLE:: finalnamefile

    n_matter= this% get_n_matter()

    DO i_matter= 1, n_matter, 1

      PRINT *, " * Finding surface for star ", i_matter, "..."
      CALL this% find_surface(this% return_center(i_matter), &
                              n_theta, n_phi, &
                              this% surfaces(i_matter)% points)
      PRINT *, "   ...done."

      IF( i_matter <= 9 ) WRITE( str_i, '(I1)' ) i_matter
      IF( i_matter >= 10 .AND. n_matter <= 99 ) &
        WRITE( str_i, '(I2)' ) i_matter
      IF( i_matter >= 100 .AND. n_matter <= 999 ) &
        WRITE( str_i, '(I3)' ) i_matter

      finalnamefile= "bns_star_surface-"//TRIM(str_i)//".dat"

      INQUIRE(FILE= TRIM(finalnamefile), EXIST= exist)

      IF(exist)THEN
        OPEN( UNIT= unit_dump, FILE= TRIM(finalnamefile), STATUS= "REPLACE", &
              FORM= "FORMATTED", POSITION= "REWIND", ACTION= "WRITE", &
              IOSTAT= ios, IOMSG= err_msg )
      ELSE
        OPEN( UNIT= unit_dump, FILE= TRIM(finalnamefile), STATUS= "NEW", &
              FORM= "FORMATTED", &
              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
      ENDIF
      IF(ios > 0)THEN
        PRINT *, "...error when opening " // TRIM(finalnamefile), &
                 ". The error message is", err_msg
        STOP
      ENDIF

      DO i= 1, n_theta, 1
        DO j= 1, n_phi, 1
          WRITE( UNIT = unit_dump, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
            i_matter, &
            this% surfaces(i_matter)% points(i,j,1), &
            this% surfaces(i_matter)% points(i,j,2), &
            this% surfaces(i_matter)% points(i,j,3), &
            this% surfaces(i_matter)% points(i,j,4), &
            this% surfaces(i_matter)% points(i,j,5), &
            this% surfaces(i_matter)% points(i,j,6)
        ENDDO
      ENDDO

      CLOSE(UNIT= unit_dump)

    ENDDO
    PRINT *

  END PROCEDURE find_print_surfaces


  MODULE PROCEDURE find_surface

    !********************************************************
    !
    !# Finds the surface of a star, using [[bnsbase::find_radius]]
    !  along many directions.
    !
    !  FT 18.02.2022
    !
    !********************************************************

    USE constants,  ONLY: pi
    USE utility,    ONLY: zero, one, two, cartesian_from_spherical

    IMPLICIT NONE

    INTEGER:: i, j
    DOUBLE PRECISION:: theta, phi
    DOUBLE PRECISION, DIMENSION(3):: direction_vector
    DOUBLE PRECISION:: radius

    ALLOCATE(surface(n_theta, n_phi,6))

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( n_theta, n_phi, surface, center, this ) &
    !$OMP             PRIVATE( i, j, theta, phi, direction_vector, radius )
    colatitude_loop: DO i= 1, n_theta, 1

      theta= DBLE(i)/DBLE(n_theta)*pi

      azimuth_loop: DO j= 1, n_phi, 1

        phi= DBLE(j)/DBLE(n_phi)*two*pi

        CALL cartesian_from_spherical(one, theta, phi, zero, zero, zero, &
          direction_vector(1), direction_vector(2), direction_vector(3) )

        radius= this% find_radius(center, direction_vector)

        CALL cartesian_from_spherical(radius, theta, phi, &
                                center(1), center(2), center(3), &
                                surface(i,j,1), surface(i,j,2), surface(i,j,3))
        surface(i,j,4)= radius
        surface(i,j,5)= theta
        surface(i,j,6)= phi

      ENDDO azimuth_loop

    ENDDO colatitude_loop
    !$OMP END PARALLEL DO

  END PROCEDURE find_surface


  MODULE PROCEDURE find_radius

    !********************************************************
    !
    !# Finds the radius of a matter object, relative to a center and along
    !  a direction. The radius is determined as the first point where the
    !  density is zero.
    !
    !  FT 27.09.2022
    !
    !********************************************************

    IMPLICIT NONE

    !INTEGER,          PARAMETER:: n_pts   = 5D+4
    DOUBLE PRECISION, PARAMETER:: min_dist= 1.D-6
    DOUBLE PRECISION, PARAMETER:: x_ini   = 130.D0
    DOUBLE PRECISION, PARAMETER:: min_dens= 1.D-12
    LOGICAL,          PARAMETER:: debug= .FALSE.

    DOUBLE PRECISION:: vector_norm
    DOUBLE PRECISION, DIMENSION(3):: point_left, point_right, versor, point_mean
    DOUBLE PRECISION:: rho_left, rho_right, x_left, x_right, x_mean, rho_mean
  
    PROCEDURE(), POINTER:: return_density

    IF(PRESENT(get_density))THEN

       return_density => read_density_opt
      
    ELSE

       return_density => read_density

    ENDIF

    radius= 0.D0

    vector_norm= NORM2(vector - center)
    versor= vector/vector_norm

    x_left = one
    x_right= x_ini
    DO

      point_left = line(x_left)
      point_right= line(x_right)
      CALL return_density(point_left(1), point_left(2), point_left(3), rho_left)
      CALL return_density(point_right(1), point_right(2), point_right(3), &
                          rho_right)

      IF( rho_left > min_dens .AND. rho_right <= min_dens )THEN

        IF( NORM2(point_left - point_right) < min_dist )THEN
          EXIT
        ENDIF
        x_mean= (x_left + x_right)/two
        point_mean= line(x_mean)
        CALL return_density(point_mean(1), point_mean(2), point_mean(3), &
                            rho_mean)
        IF(rho_mean > min_dens)THEN
          x_left = x_mean
        ELSE
          x_right= x_mean
        ENDIF

      ELSEIF( rho_left > min_dens .AND. rho_right > min_dens )THEN

        x_right= two*x_right

      ELSE

        PRINT *
        PRINT *, "** ERROR in SUBROUTINE find_radius in SUBMODULE ", &
                  "bns_base@geometry!"
        PRINT *, "x_left=", x_left
        PRINT *, "x_right=", x_right
        PRINT *, "point_left=", point_left
        PRINT *, "point_right=", point_right
        PRINT *, "rho_left=", rho_left
        PRINT *, "rho_right=", rho_right
        PRINT *
        STOP

      ENDIF

    ENDDO

    IF(debug) PRINT *, x_left
    IF(debug) PRINT *, x_right
    IF(debug) PRINT *, point_left - center
    IF(debug) PRINT *, point_right - center

    radius= NORM2(point_left - center)

!    DO i= 0, n_pts, 1
!    ! This loop cannot be OMP-parallelized because FUKA is not thread-safe
!    ! Luckily, it does not take too long, unless n_pts is set too large

!       x= x_left + DBLE(i)/DBLE(n_pts)*(x_right - x_left)
!       point= line(x)
!       CALL return_density(point(1), point(2), point(3), rho)

!       IF( rho <= min_dens )THEN
!         radius= NORM2(point - center)
!         EXIT
!       ENDIF

!    ENDDO


    CONTAINS

    
    PURE FUNCTION line(x) RESULT(point)
    !# The line \(c+xv\) parametrized by x, with \(c\) being `center`
    !  and \(v\) being `vector`

      DOUBLE PRECISION, INTENT(IN):: x
      !! Parameter along the line
      DOUBLE PRECISION, DIMENSION(3):: point
      !! Coordinates of the point on the line at the given \(x\)

      point= center + x*versor

    END FUNCTION line


    SUBROUTINE read_density(x, y, z, rho)
    !# Reads the density using the method [[id_base:read_mass_density]].
    !  This is a wrapper needed as the TARGET of a PROCEDURE POINTER
    !  (PROCEDURE POINTERS cannot point to TYPE-BOUND PROCEDURES unfortunately)

      DOUBLE PRECISION, INTENT(IN) :: x
      DOUBLE PRECISION, INTENT(IN) :: y
      DOUBLE PRECISION, INTENT(IN) :: z
      DOUBLE PRECISION, INTENT(OUT):: rho

      rho= this% read_mass_density(x, y, z)

    END SUBROUTINE read_density


    SUBROUTINE read_density_opt(x, y, z, rho)
    !# Reads the density using the method passed optionally as an argument.

      DOUBLE PRECISION, INTENT(IN) :: x
      DOUBLE PRECISION, INTENT(IN) :: y
      DOUBLE PRECISION, INTENT(IN) :: z
      DOUBLE PRECISION, INTENT(OUT):: rho

      rho= get_density(x, y, z)

    END SUBROUTINE read_density_opt


  END PROCEDURE find_radius


  MODULE PROCEDURE find_center

    !********************************************************
    !
    !# Finds the center of a star, as the point where the
    !  density is maximal.
    !
    !  FT 27.09.2022
    !
    !********************************************************

    IMPLICIT NONE

    INTEGER, PARAMETER:: n_pts= 5000
    LOGICAL, PARAMETER:: debug= .FALSE.

    INTEGER:: i
    DOUBLE PRECISION:: x
    DOUBLE PRECISION:: rho, rho_max, x_left, x_right
  
    PROCEDURE(), POINTER:: return_density

    IF( x_sign /= one .AND. x_sign /= -one )THEN
       PRINT *, "** ERROR in FUNCTION find_center! The argument x_sign ", &
                "should be a DOUBLE PRECISION 1 or -1."
       PRINT *, " * Stopping..."
       PRINT *
    ENDIF

    IF(PRESENT(get_density))THEN

       return_density => read_density_opt
      
    ELSE

       return_density => read_density

    ENDIF

    x_left = (x_sign - one)*separation
    x_right= (x_sign + one)*separation
    rho_max= -HUGE(one)
    DO i= 0, n_pts, 1

      x= x_left + DBLE(i)/DBLE(n_pts)*(x_right - x_left)
      CALL return_density(x, zero, zero, rho)

      IF( rho >= rho_max )THEN
        rho_max= rho
        center = x
      ENDIF

    ENDDO


    CONTAINS


    SUBROUTINE read_density(x, y, z, rho)
    !# Reads the density using the method [[id_base:read_mass_density]].
    !  This is a wrapper needed as the TARGET of a PROCEDURE POINTER
    !  (PROCEDURE POINTERS cannot point to TYPE-BOUND PROCEDURES unfortunately)

      DOUBLE PRECISION, INTENT(IN) :: x
      DOUBLE PRECISION, INTENT(IN) :: y
      DOUBLE PRECISION, INTENT(IN) :: z
      DOUBLE PRECISION, INTENT(OUT):: rho

      rho= this% read_mass_density(x, y, z)

    END SUBROUTINE read_density


    SUBROUTINE read_density_opt(x, y, z, rho)
    !# Reads the density using the method passed optionally as an argument.

      DOUBLE PRECISION, INTENT(IN) :: x
      DOUBLE PRECISION, INTENT(IN) :: y
      DOUBLE PRECISION, INTENT(IN) :: z
      DOUBLE PRECISION, INTENT(OUT):: rho

      rho= get_density(x, y, z)

    END SUBROUTINE read_density_opt


  END PROCEDURE find_center


END SUBMODULE geometry

