! File:         module_SPHINCS_LORENE.f90
! Author:       Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

MODULE sphincs_lorene


    !**********************************************
    !                                             *
    ! This module uses the other                  *
    ! modules needed to import and set up the     *
    ! LORENE ID                                   *
    !                                             *
    ! FT 23.10.2020                               *
    !                                             *
    !**********************************************


    USE utility,        ONLY: date, time, zone, values, run_id, itr, itr3, &
                              itr4, ios, err_msg, file_exists, cnt, &
                              test_status, show_progress
    USE timing,         ONLY: timer
    USE bns_id,         ONLY: bns
    USE particles_id,   ONLY: particles
    USE formul_3p1_id,  ONLY: formul_3p1
    USE formul_bssn_id, ONLY: bssn_id


    IMPLICIT NONE


END MODULE sphincs_lorene
