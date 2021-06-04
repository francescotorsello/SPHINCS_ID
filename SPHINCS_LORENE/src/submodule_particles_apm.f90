! File:         submodule_particles_apm.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (particles_id) particles_apm

  !****************************************************
  !                                                   *
  ! Implementation of the method                      *
  ! compute_and_export_SPH_variables_apm              *
  ! of TYPE particles.                                *
  !                                                   *
  ! FT 04.06.2021                                     *
  !                                                   *
  !****************************************************


  IMPLICIT NONE


  CONTAINS


  MODULE PROCEDURE compute_and_export_SPH_variables_apm

    !****************************************************
    !                                                   *
    ! Compute  the particle positions as follows:       *
    !                                                   *
    !   1. Take initial particle distribution as input  *
    !   2. Assume that the particles have the same mass *
    !   3. Do the APM iteration so that the final       *
    !      particle number density matches the baryon   *
    !      density in the star                          *
    !   4. Correct the particle masses ONCE in order    *
    !      to match the density even better. Since we   *
    !      don't want a large mass ratio, we impose a   *
    !      maximum mass ratio when performing this      *
    !      correction.                                  *
    !                                                   *
    ! After this procedure, the resulting particle      *
    ! distribution has positions and baryon numbers     *
    ! that kernel-estimate very well the mass density   *
    ! of the star, and has a low mass ratio.            *
    !                                                   *
    ! This procedure assigns positions and nu. After it *
    ! is performed, the other SPH quantities are        *
    ! computed, and then they are printed to a binary   *
    ! file ready to be used by SPHINCS_BSSN.            *
    !                                                   *
    ! FT 04.06.2021                                     *
    !                                                   *
    !****************************************************

  END PROCEDURE compute_and_export_SPH_variables_apm


END SUBMODULE particles_apm
