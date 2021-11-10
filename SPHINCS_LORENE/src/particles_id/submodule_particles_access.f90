! File:         submodule_particles_access.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (particles_id) particles_access

  !**************************************************
  !
  !# This SUBMODULE contains the implementation of
  !  the methods of TYPE particles
  !  that allow to access PRIVATE members.
  !
  !  FT 12.07.2021
  !
  !**************************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE impose_equatorial_plane_symmetry

    !*************************************************************
    !
    !# Mirror the particle with z>0 with respect to the xy plane,
    !  to impose the equatorial-plane symmetry
    !
    !  FT 1.09.2021
    !
    !*************************************************************

    USE analyze, ONLY: COM

    IMPLICIT NONE

   ! INTEGER, INTENT(IN):: npart
   ! DOUBLE PRECISION, INTENT(IN), OPTIONAL:: com_star
   ! LOGICAL, INTENT(IN), OPTIONAL:: verbose
   !
   ! DOUBLE PRECISION, DIMENSION(3,npart), INTENT(INOUT):: pos
   ! DOUBLE PRECISION, DIMENSION(npart),   INTENT(INOUT):: nu

    INTEGER:: a, itr, npart_half
    DOUBLE PRECISION:: com_x, com_y, com_z, com_d

    DOUBLE PRECISION, DIMENSION(3,npart):: pos_tmp
    DOUBLE PRECISION, DIMENSION(npart)  :: nu_tmp

    pos_tmp= pos
    nu_tmp= nu
    itr= 0
    DO a= 1, npart, 1
      IF( pos_tmp( 3, a ) > 0.0D0 &
          .AND. &
          itr < npart/2 )THEN
        itr= itr + 1
        pos( 1, itr )= pos_tmp( 1, a )
        pos( 2, itr )= pos_tmp( 2, a )
        pos( 3, itr )= pos_tmp( 3, a )
        IF( PRESENT(nu) ) nu( itr )= nu_tmp( a )
      ENDIF
    ENDDO
    npart_half= itr

    ! If some of the particles crossed the xy plane top-down in the
    ! last step, replace them with their previous position
    ! above the xy plane
  !   IF( npart_half < npart/2 )THEN
  !
  !     npart_missing= npart/2 - npart_half
  !
  !     DO a= npart_half + 1, npart/2, 1
  !
  !       pos( :, a )= all_pos_tmp2( :, a )
  !
  !     ENDDO
  !
  !   ENDIF

    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( pos, npart_half, nu ) &
    !$OMP             PRIVATE( a )
    DO a= 1, npart_half, 1
      pos( 1, npart_half + a )=   pos( 1, a )
      pos( 2, npart_half + a )=   pos( 2, a )
      pos( 3, npart_half + a )= - pos( 3, a )
      IF( PRESENT(nu) ) nu( npart_half + a )=   nu( a )
    ENDDO
    !$OMP END PARALLEL DO

    npart= 2*npart_half

    IF( PRESENT(verbose) .AND. verbose .EQV. .TRUE. )THEN

      CALL COM( npart, pos, nu, & ! input
                com_x, com_y, com_z, com_d ) ! output

      PRINT *, "** After mirroring particles:"
      IF( PRESENT(com_star) ) PRINT *, &
               " * x coordinate of the center of mass of the star, ", &
               "from LORENE: com_star= ", com_star, "Msun_geo"
      PRINT *, " * x coordinate of the center of mass of the particle ", &
               "distribution: com_x= ", com_x, "Msun_geo"
      PRINT *, " * y coordinate of the center of mass of the particle ", &
               "distribution: com_y= ", com_y, "Msun_geo"
      PRINT *, " * z coordinate of the center of mass of the particle ", &
               "distribution: com_z= ", com_z, "Msun_geo"
      PRINT *, " * Distance of the center of mass of the particle ", &
               "distribution from the  origin: com_d= ", com_d
      IF( PRESENT(com_star) ) PRINT *, " * |com_x-com_star/com_star|=", &
               ABS( com_x-com_star )/ABS( com_star + 1 )
      PRINT *

    ENDIF

  END PROCEDURE impose_equatorial_plane_symmetry


  !-----------------!
  !--  FUNCTIONS  --!
  !-----------------!


  MODULE PROCEDURE get_npart

    !************************************************
    !
    !# Returns the total number of particles
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    n_part= THIS% npart

  END PROCEDURE get_npart


  MODULE PROCEDURE get_npart1

    !************************************************
    !
    !# Returns the number of particles on star 1
    !
    !  FT 27.04.2021
    !
    !************************************************

    IMPLICIT NONE

    n_part= THIS% npart_i(1)

  END PROCEDURE get_npart1


  MODULE PROCEDURE get_npart2

    !************************************************
    !
    !# Returns the number of particles on star 2
    !
    !  FT 27.04.2021
    !
    !************************************************

    IMPLICIT NONE

    n_part= THIS% npart_i(2)

  END PROCEDURE get_npart2


  MODULE PROCEDURE get_nuratio

    !************************************************
    !
    !# Returns the baryon number ratio on the stars
    !
    !  FT 27.04.2021
    !
    !************************************************

    IMPLICIT NONE

    nuratio= THIS% nuratio

  END PROCEDURE get_nuratio


  MODULE PROCEDURE get_nuratio1

    !************************************************
    !
    !# Returns the baryon number ratio on star 1
    !
    !  FT 27.04.2021
    !
    !************************************************

    IMPLICIT NONE

    nuratio1= THIS% nuratio_i(1)

  END PROCEDURE get_nuratio1


  MODULE PROCEDURE get_nuratio2

    !************************************************
    !
    !# Returns the baryon number ratio on star 2
    !
    !  FT 27.04.2021
    !
    !************************************************

    IMPLICIT NONE

    nuratio2= THIS% nuratio_i(2)

  END PROCEDURE get_nuratio2


  MODULE PROCEDURE get_pos

    !************************************************
    !
    !# Returns the array of particle positions
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    pos_u= THIS% pos

  END PROCEDURE get_pos


  MODULE PROCEDURE get_vel

    !************************************************
    !
    !# Returns the array of coordinate 3-velocity of
    ! particles
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    vel= THIS% v(1:3,:)

  END PROCEDURE get_vel


  MODULE PROCEDURE get_nlrf

    !************************************************
    !
    !# Returns the array of baryon density in the
    ! local rest frame
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    nlrf= THIS% nlrf

  END PROCEDURE get_nlrf


  MODULE PROCEDURE get_nu

    !************************************************
    !
    !# Returns the array of baryon per particle
    ! [baryon (Msun_geo)^{-3}]
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    nu= THIS% nu

  END PROCEDURE get_nu


  MODULE PROCEDURE get_u

    !************************************************
    !
    !# Returns the array of specific internal
    ! energy [c^2]
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    u= THIS% specific_energy_parts

  END PROCEDURE get_u


  MODULE PROCEDURE get_pressure

    !************************************************
    !
    !# Returns the array of pressure [kg c^2 m^{-3}]
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    pressure= THIS% pressure_parts

  END PROCEDURE get_pressure


  MODULE PROCEDURE get_pressure_cu

    !************************************************
    !
    !# Returns the array of pressure in code units
    ! [amu*c**2/(Msun_geo**3)]
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    pressure_cu= THIS% pressure_parts_cu

  END PROCEDURE get_pressure_cu


  MODULE PROCEDURE get_theta

    !************************************************
    !
    !# Returns the array of generalized Lorentz
    ! factor
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    theta= THIS% Theta

  END PROCEDURE get_theta


  MODULE PROCEDURE get_h

    !************************************************
    !
    !# Returns the array of initial guess for the
    ! smoothing length [Msun_geo]
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    h= THIS% h

  END PROCEDURE get_h


  MODULE PROCEDURE is_empty

    !************************************************
    !
    !# Returns the variable empty_object
    !  @warning experimental, not actively used in
    !          the code yet
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    answer= THIS% empty_object

  END PROCEDURE is_empty


END SUBMODULE particles_access
