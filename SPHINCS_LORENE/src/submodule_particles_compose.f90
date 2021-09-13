! File:         submodule_particles_compose.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (particles_id) particles_compose

  !***************************************************
  !
  !# This SUBMODULE contains the implementation of
  !  the methods of TYPE particles
  !  that compute Ye on the particles, using the
  !  data from the CompOSE database
  !
  !  https://compose.obspm.fr/
  !
  !  FT 12.07.2021
  !
  !***************************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE read_compose_composition

    !************************************************
    !
    !# Reads the electron fraction Y_e = n_e/n_b,
    !  with n_e electron number density and n_b
    !  baryon number density, from the .compo file
    !  taken from the CompOSE database of EoS.
    !  Y_e is given as a function of T, n_b, Y_q on
    !  a grid; the computation of Ye on the stars is
    !  done by the SUBROUTINE compute_Ye_on_stars.
    !
    !  FT 1.03.2021
    !
    !************************************************

    USE constants, ONLY: fm2cm, cm2km, km2Msun_geo

    IMPLICIT NONE

    ! The commented variables might be useful in the future

    INTEGER:: itr, cntr!, &
              !i_t, i_nb, i_phase, n_pairs, i_e, i_n, Y_n, &
              !n_quad, &
              !i_i, A_i, Z_i, Y_i, i_leptons, i_ns
    INTEGER, PARAMETER:: unit_compose= 56
    INTEGER, PARAMETER:: max_length_eos= 10000

    !DOUBLE PRECISION:: m_n, m_p, &
    !                   p_nb, s_nb, mub_mn, muq_mn, mul_mn, f_nbmn, e_nbmn, h

    !DOUBLE PRECISION, DIMENSION( : ), ALLOCATABLE:: n_b, Y_e

    LOGICAL:: exist

    CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile

    PRINT *, "** Executing the read_compose_composition subroutine..."

    ALLOCATE( THIS% nb_table( max_length_eos ) )
    ALLOCATE( THIS% Ye_table( max_length_eos ) )
    THIS% nb_table= 0.0D0
    THIS% Ye_table= 0.0D0


    IF( PRESENT(namefile) )THEN
      finalnamefile= TRIM(namefile)//".beta"
    ELSE
      finalnamefile= "../../CompOSE_EOS/SFHO_with_electrons/eos.beta"
    ENDIF

    INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

    IF( exist )THEN
      OPEN( UNIT= unit_compose, FILE= TRIM(finalnamefile), &
            FORM= "FORMATTED", ACTION= "READ", IOSTAT= ios, &
            IOMSG= err_msg )
      IF( ios > 0 )THEN
        PRINT *, "...error when opening " // TRIM(finalnamefile), &
                ". The error message is", err_msg
        STOP
      ENDIF
      !CALL test_status( ios, err_msg, "...error when opening " &
      !                  // TRIM(finalnamefile) )
    ELSE
      PRINT *, "** ERROR! Unable to find file " // TRIM(finalnamefile)
      STOP
    ENDIF

    !READ( UNIT= unit_compose, FMT= *, IOSTAT = ios, IOMSG= err_msg ) &
    !                                                    m_n, m_p, i_leptons

    PRINT *, " * Reading file " // TRIM(finalnamefile) // "..."
    cntr= 0
    read_compose_beta: DO itr= 1, max_length_eos, 1
      READ( UNIT= unit_compose, FMT= *, IOSTAT = ios, IOMSG= err_msg ) &
                        ! Variables for .thermo.ns file
                        !i_T, i_nb, i_yq, &
                        !p_nb, s_nb, mub_mn, muq_mn, mul_mn, f_nbmn, e_nbmn, &
                        !i_ns, THIS% Ye_table(itr), h
                        ! Variables for .beta file
                        THIS% nb_table(itr), &
                        THIS% Ye_table(itr)
                        ! Variables for .compo file
                        !i_T, i_nb, i_yq, &
                        !i_phase, n_pairs, &
                        !i_e, Y_e(itr)!, &
                        !i_n, Y_n, &
                        !n_quad, &
                        !i_i, A_i, Z_i, Y_i
      IF( ios > 0 )THEN
        PRINT *, "...error when reading " // TRIM(finalnamefile), &
                ". The error message is", err_msg
        STOP
      ENDIF
      IF( ios < 0 )THEN
        PRINT *, " * Reached end of file " // TRIM(finalnamefile)
        EXIT
      ENDIF
      cntr= cntr + 1
    ENDDO read_compose_beta
    !PRINT *, "cntr= ", cntr
    ! Reallocate the arrays to delete the empty elements
    THIS% nb_table= THIS% nb_table( 1:cntr )/(fm2cm**3*cm2km**3*km2Msun_geo**3)
    THIS% Ye_table= THIS% Ye_table( 1:cntr )
    !PRINT *, i_T, i_nb, i_yq, i_phase, n_pairs, i_e
    !PRINT *, "SIZE(n_b)= ", SIZE(THIS% n_b), "SIZE(Y_e)= ", SIZE(THIS% Y_e)
    !PRINT *, "n_b(1)= ", THIS% n_b(1), "Y_e(1)= ", THIS% Y_e(1)
    !PRINT *, "n_b(cntr)= ", THIS% n_b(cntr), "Y_e(cntr)= ", THIS% Y_e(cntr)
    !STOP
    CLOSE( unit_compose )

    PRINT *, "** Subroutine read_compose_composition executed."
    PRINT *

  END PROCEDURE read_compose_composition


  MODULE PROCEDURE compute_Ye

    !************************************************
    !
    !# Interpolates the electron fraction
    !  Y_e = n_e/n_b
    !  at the particle positions, using the data
    !  read by read_compose_composition.
    !
    !  FT 3.03.2021
    !
    !************************************************

    IMPLICIT NONE

    INTEGER:: d, itr2

    DOUBLE PRECISION:: min_nb_table, max_nb_table

    d= SIZE(THIS% nb_table)
    min_nb_table= MINVAL( THIS% nb_table, 1 )
    max_nb_table= MAXVAL( THIS% nb_table, 1 )

    particle_loop: DO itr= 1, THIS% npart, 1

      ! There may be particles with initially negative hydro fields,
      ! close to the surface of the stars. The present
      ! version of the code sets their hydro fields to 0.
      ! In turn, this is a problem if we read Ye from the .beta file
      ! from the CompOSe database since the value of 0
      ! for the baryon number density (used in the interpolation
      ! of Ye) is not included in the .beta file.
      ! In other words, Ye cannot be interpolated on those particles,
      ! and the present version of the code sets Ye to 0 as well.

      ! The problem above is solved: particles are not placed if the hydro
      ! is negative or zero.

      !IF( THIS% nlrf(itr) == 0.0D0 )THEN
        !THIS% Ye(itr)= 0.0D0
        !CYCLE
      !ELSE
      IF( THIS% nlrf(itr) < min_nb_table )THEN
        PRINT *, "** ERROR! The value of nlrf(", itr, ")=", THIS% nlrf(itr), &
                 "is lower than the minimum value in the table =", min_nb_table
        PRINT *, " * Is nlrf computed when you call this SUBROUTINE? " // &
                 "If yes, please select a table with a wider range."
        STOP
      ELSEIF( THIS% nlrf(itr) > max_nb_table )THEN
        PRINT *, "** ERROR! The value of nlrf(", itr, ")=", THIS% nlrf(itr), &
                 "is larger than the maximum value in the table =", max_nb_table
        PRINT *, " * Is nlrf computed when you call this SUBROUTINE? " // &
                 "If yes, please select a table with a wider range."
        STOP
      ENDIF

      Ye_linear_interpolation_loop: DO itr2= 1, d - 1, 1

        IF( THIS% nb_table(itr2) < THIS% nlrf(itr) .AND. &
            THIS% nlrf(itr) < THIS% nb_table(itr2 + 1) )THEN

          THIS% Ye(itr)= THIS% Ye_table(itr2) &
                         + (THIS% Ye_table(itr2 + 1) - THIS% Ye_table(itr2))/ &
                         (THIS% nb_table(itr2 + 1) - THIS% nb_table(itr2)) &
                         *(THIS% nlrf(itr) - THIS% nb_table(itr2))
          EXIT

        ENDIF

      ENDDO Ye_linear_interpolation_loop

      !PRINT *, "Ye(", itr, ")=", THIS% Ye(itr)
      !PRINT *

    ENDDO particle_loop

  END PROCEDURE compute_Ye


END SUBMODULE particles_compose
