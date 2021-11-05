! File:         submodule_particles_io.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (particles_id) particles_io

  !***************************************************
  !
  !# This submodule contains the implementation of the
  !  methods of TYPE particles that handle I/O (input/output)
  !
  !  FT 5.11.2021
  !
  !***************************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE print_summary

    !************************************************
    !
    !# Prints a summary of the properties of the |sph| particle
    !  distribution, optionally, to a formatted file whose name
    !  is given as the optional argument `filename`
    !
    !  FT 5.11.2021
    !
    !************************************************

    IMPLICIT NONE

    INTEGER:: i_matter

    PRINT *, " * SPH:"
    PRINT *
    PRINT *, "   Total particle number= ", THIS% npart
    DO i_matter= 1, THIS% n_matter, 1
      PRINT *, "   Particle number on matter object ", i_matter, "= ", &
                                            THIS% npart_i(i_matter)
    ENDDO
    PRINT *
    DO i_matter= 1, THIS% n_matter, 1
      PRINT *, "   Mass fraction of matter object", i_matter, "=", &
               THIS% mass_fractions(i_matter)
      PRINT *, "   Particle fraction of matter object", i_matter, "=", &
               THIS% npart_i(i_matter)/THIS% npart
      PRINT *, "   Baryon number ratio on matter object", i_matter, "=", &
               THIS% nuratio_i(i_matter)
    ENDDO
    PRINT *

    PRINT *, "   Baryon number ratio across all matter objects=", THIS% nuratio
    PRINT *

  END PROCEDURE print_summary


END SUBMODULE particles_io
