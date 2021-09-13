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

    n_part= THIS% npart1

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

    n_part= THIS% npart2

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

    nuratio1= THIS% nuratio1

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

    nuratio2= THIS% nuratio2

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
    ! (experimental)
    !
    !  FT
    !
    !************************************************

    IMPLICIT NONE

    answer= THIS% empty_object

  END PROCEDURE is_empty


END SUBMODULE particles_access
