! File:         submodule_ejecta_generic_io.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (ejecta_generic) ejecta_generic_io

  !***************************************************
  !
  !# This submodule contains the implementation of the
  !  methods of TYPE ejecta that handle I/O (input/output)
  !
  !  FT xx.11.2021
  !
  !***************************************************


  IMPLICIT NONE


  CONTAINS


  !------------------------------!
  !--  OVERRIDING SUBROUTINES  --!
  !------------------------------!


  MODULE PROCEDURE print_summary_ejecta

    !************************************************
    !
    !# Prints a summary of the physical properties the system
    !  to the standard output and, optionally, to a formatted
    !  file whose name is given as the optional argument `filename`
    !  @todo to be implemented
    !
    !  FT xx.11.2021
    !
    !************************************************

    IMPLICIT NONE

  END PROCEDURE print_summary_ejecta


END SUBMODULE ejecta_generic_io
