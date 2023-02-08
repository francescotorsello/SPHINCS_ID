! File:         submodule_sph_particles_quality_indicators.f90
! Authors:      Francesco Torsello (FT)
!************************************************************************
! Copyright (C) 2020-2023 Francesco Torsello                            *
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

SUBMODULE(sph_particles) quality_indicators

  !****************************************************
  !
  !# This SUBMODULE contains the implementation of
  !  the methods of TYPE particles
  !  that computes the quality indicators, referring to
  !
  !  Daniel J. Price, Smoothed Particle Hydrodynamics and
  !  Magnetohydrodynamics.
  !  Journal of Computational Physics, 231, 3, 759-794 (2012)
  !  DOI: 10.1016/j.jcp.2010.12.011
  !  http://arxiv.org/abs/1012.1885
  !  eqs.(64), (67) and (74-75)
  !
  !  Rosswog, S. SPH Methods in the Modelling of Compact Objects.
  !  Living Rev Comput Astrophys 1, 1 (2015).
  !  https://doi.org/10.1007/lrca-2015-1.
  !  https://arxiv.org/abs/1406.4224.
  !  eqs.(6) and (9)
  !
  !  FT 05.10.2022
  !
  !****************************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE compute_and_print_quality_indicators

    !****************************************************
    !
    !# Compute the quality indicators, referring to
    !
    !  Daniel J. Price, Smoothed Particle Hydrodynamics and
    !  Magnetohydrodynamics.
    !  Journal of Computational Physics, 231, 3, 759-794 (2012)
    !  DOI: 10.1016/j.jcp.2010.12.011
    !  http://arxiv.org/abs/1012.1885
    !  eqs.(64), (67) and (74-75)
    !
    !  Rosswog, S. SPH Methods in the Modelling of Compact Objects.
    !  Living Rev Comput Astrophys 1, 1 (2015).
    !  https://doi.org/10.1007/lrca-2015-1.
    !  https://arxiv.org/abs/1406.4224.
    !  eqs.(6) and (9)
    !
    !  FT 05.10.2022
    !
    !****************************************************

    USE kernel_table, ONLY: dv_table, dv_table_1,&
                            W_no_norm, n_tab_entry, &!dWdv_no_norm, &
                            interp_gradW_table,interp_W_gradW_table
    USE sphincs_sph,  ONLY: ncand, all_clists
    USE RCB_tree_3D,  ONLY: iorig, nic, nfinal, nprev, lpart, rpart
    USE options,      ONLY: ndes
    USE set_h,        ONLY: exact_nei_tree_update
    USE utility,      ONLY: zero, one, two

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(npart)    :: qi_1
    DOUBLE PRECISION, DIMENSION(3,npart)  :: qi_2
    DOUBLE PRECISION, DIMENSION(3,npart)  :: qi_3
    DOUBLE PRECISION, DIMENSION(3,3,npart):: qi_4
    DOUBLE PRECISION, DIMENSION(3,npart)  :: qi_ham

    INTEGER:: a, b, ill, l, itot, inde, index1, k

    DOUBLE PRECISION:: ha, ha_1, ha_3, ha_4, va, xa, ya, za, &
                       hb, hb_1, hb_3, hb_4, xb, yb, zb, rab, rab2, rab_1, &
                       ha2, ha2_4, hb2_4, dx, dy, dz

    DOUBLE PRECISION:: Wab_ha, Wi, Wi1, dvv, &
                       grW_ha_x, grW_ha_y, grW_ha_z, &
                       grW_hb_x, grW_hb_y, grW_hb_z, eab(3), &
                       Wa, grW, grWa, grWb, vb, vol_b

    LOGICAL, PARAMETER:: debug= .TRUE.

    INTEGER                      :: ios, unit_qi
    LOGICAL                      :: exist
    CHARACTER(LEN=:), ALLOCATABLE:: namefile, err_msg

    PRINT *, " * Computing the quality indicators..."

    CALL exact_nei_tree_update( ndes, npart, pos, nu )
    ! exact_nei_tree_update updates h

    qi_1  = zero
    qi_2  = zero
    qi_3  = zero
    qi_4  = zero
    qi_ham= zero
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( nfinal, nprev, iorig, lpart, rpart, nic,nu,nstar, &
    !$OMP                     ncand, all_clists, W_no_norm, pos, h, &
    !$OMP                     qi_1, qi_2, qi_3, qi_4, qi_ham ) &
    !$OMP             PRIVATE( ill, itot, a, b, l, &
    !$OMP                      ha, ha_1, ha_3, ha_4, ha2, ha2_4, &
    !$OMP                      hb, hb_1, hb_3, hb_4, hb2_4, &
    !$OMP                      xa, ya, za, xb, yb, zb, dx, dy, dz, rab2, &
    !$OMP                      inde, va, dv_table_1, index1, Wi, &
    !$OMP                      dvv, dv_table, Wab_ha, Wa, grW, Wi1, &
    !$OMP                      grW_ha_x, grW_ha_y, grW_ha_z, grWa, grWb, &
    !$OMP                      grW_hb_x, grW_hb_y, grW_hb_z, eab, rab, vb, &
    !$OMP                      rab_1, vol_b )
    ll_cell_loop: DO ill= 1, nfinal

      itot= nprev + ill
      IF( nic(itot) == 0 ) CYCLE

      particle_loop: DO l= lpart(itot), rpart(itot)

        a    = iorig(l)

        ha   = h(a)
        ha_1 = one/ha
        ha_3 = ha_1*ha_1*ha_1
        ha_4 = ha_3*ha_1
        ha2  = ha*ha
        ha2_4= two*two*ha2

        xa= pos(1,a)
        ya= pos(2,a)
        za= pos(3,a)

        !prgNa= (pr_sph(a)*sqg(a)/(nstar_sph(a)/m0c2_cu)**2)

        cand_loop: DO k= 1, ncand(ill)

          b= all_clists(ill)% list(k)

          hb   = h(b)
          hb_1 = one/hb
          hb_3 = hb_1*hb_1*hb_1
          hb_4 = hb_3*hb_1
          hb2_4= two*two*hb*hb

          xb= pos(1,b)  ! CONTRA-variant
          yb= pos(2,b)
          zb= pos(3,b)

          ! potentially bail out
          dx= xa - xb
          dy= ya - yb
          dz= za - zb

          rab2 = dx*dx + dy*dy + dz*dz
          rab  = SQRT(rab2) + 1.D-30
          rab_1= 1.D0/rab
          IF( rab2 > ha2_4 .AND. rab2 > hb2_4 ) CYCLE

          !--------!
          !-- ha --!
          !--------!
          va= rab*ha_1

          ! get interpolation indices
          inde  = MIN(INT(va*dv_table_1),n_tab_entry)
          index1= MIN(inde + 1,n_tab_entry)

          ! get tabulated values
          Wi = W_no_norm(inde)
          Wi1= W_no_norm(index1)

          ! interpolate
          dvv   = (va - DBLE(inde)*dv_table)*dv_table_1
          Wab_ha= Wi + (Wi1 - Wi)*dvv

          ! unit vector ("a-b" --> from b to a)
          eab(1)= dx*rab_1
          eab(2)= dy*rab_1
          eab(3)= dz*rab_1

          !--------!
          !-- ha --!
          !--------!

          ! kernel and its gradient
          !DIR$ INLINE
          CALL interp_W_gradW_table( va, Wa, grW )
          Wa = Wa*ha_3
          grW= grW*ha_4

          ! nabla_a Wab(ha)
          grW_ha_x= grW*eab(1)
          grW_ha_y= grW*eab(2)
          grW_ha_z= grW*eab(3)

    !      grWa= grW_ha_x*eab(1) + &
    !            grW_ha_y*eab(2) + &
    !            grW_ha_z*eab(3)

          !--------!
          !-- hb --!
          !--------!
          vb= rab*hb_1

          ! kernel and its gradient
          !DIR$ INLINE
          CALL interp_gradW_table(vb,grW)
          grW= grW*hb_4

          ! nabla_a Wab(hb)
          grW_hb_x= grW*eab(1)
          grW_hb_y= grW*eab(2)
          grW_hb_z= grW*eab(3)

          vol_b= nu(b)/nstar(b)

    !      grWb= grW_hb_x*eab(1) + &
    !            grW_hb_y*eab(2) + &
    !            grW_hb_z*eab(3)

          qi_1(a)= qi_1(a) + Wab_ha*vol_b

          qi_2(1,a)= qi_2(1,a) + (-dx)*Wab_ha*vol_b
          qi_2(2,a)= qi_2(2,a) + (-dy)*Wab_ha*vol_b
          qi_2(3,a)= qi_2(3,a) + (-dz)*Wab_ha*vol_b

          qi_3(1,a)= qi_3(1,a) + grW_ha_x*vol_b
          qi_3(2,a)= qi_3(2,a) + grW_ha_y*vol_b
          qi_3(3,a)= qi_3(3,a) + grW_ha_z*vol_b

          qi_4(1,1,a)= qi_4(1,1,a) + (-dx)*grW_ha_x*vol_b
          qi_4(1,2,a)= qi_4(1,2,a) + (-dx)*grW_ha_y*vol_b
          qi_4(1,3,a)= qi_4(1,3,a) + (-dx)*grW_ha_z*vol_b
          qi_4(2,1,a)= qi_4(2,1,a) + (-dy)*grW_ha_x*vol_b
          qi_4(2,2,a)= qi_4(2,2,a) + (-dy)*grW_ha_y*vol_b
          qi_4(2,3,a)= qi_4(2,3,a) + (-dy)*grW_ha_z*vol_b
          qi_4(3,1,a)= qi_4(3,1,a) + (-dz)*grW_ha_x*vol_b
          qi_4(3,2,a)= qi_4(3,2,a) + (-dz)*grW_ha_y*vol_b
          qi_4(3,3,a)= qi_4(3,3,a) + (-dz)*grW_ha_z*vol_b

          qi_ham(1,a)= qi_ham(1,a) &
                       + nu(b)*grW_ha_x*(one/nstar(a) + nstar(a)/nstar(b)**2)
          qi_ham(2,a)= qi_ham(2,a) &
                       + nu(b)*grW_ha_y*(one/nstar(a) + nstar(a)/nstar(b)**2)
          qi_ham(3,a)= qi_ham(3,a) &
                       + nu(b)*grW_ha_z*(one/nstar(a) + nstar(a)/nstar(b)**2)

        ENDDO cand_loop

        !qi(a)= qi(a)/(nu(a)/nstar(a))

      ENDDO particle_loop

    ENDDO ll_cell_loop
    !$OMP END PARALLEL DO

    PRINT *, " * Printing the quality indicators to file..."

    IF(PRESENT(path))THEN
      namefile= path//"quality_indicators.dat"
    ELSE
      namefile= "quality_indicators.dat"
    ENDIF
    unit_qi= 279465

    INQUIRE( FILE= TRIM(namefile), EXIST= exist )

    IF( exist )THEN
        OPEN( UNIT= unit_qi, FILE= TRIM(namefile), STATUS= "REPLACE", &
              FORM= "FORMATTED", &
              POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
              IOMSG= err_msg )
    ELSE
        OPEN( UNIT= unit_qi, FILE= TRIM(namefile), STATUS= "NEW", &
              FORM= "FORMATTED", &
              ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
    ENDIF
    IF( ios > 0 )THEN
      PRINT *, "...error when opening " // TRIM(namefile), &
               ". The error message is", err_msg
      STOP
    ENDIF

    DO a= 1, npart, 1
      WRITE( UNIT = unit_qi, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
        a, &           ! 1
        pos(1,a), &    ! 2
        pos(2,a), &    ! 3
        pos(3,a), &    ! 4
        qi_1(a), &     ! 5
        qi_2(1,a), &   ! 6
        qi_2(2,a), &   ! 7
        qi_2(3,a), &   ! 8
        qi_3(1,a), &   ! 9
        qi_3(2,a), &   ! 10
        qi_3(3,a), &   ! 11
        qi_4(1,1,a), & ! 12
        qi_4(1,2,a), & ! 13
        qi_4(1,3,a), & ! 14
        qi_4(2,1,a), & ! 15
        qi_4(2,2,a), & ! 16
        qi_4(2,3,a), & ! 17
        qi_4(3,1,a), & ! 18
        qi_4(3,2,a), & ! 19
        qi_4(3,3,a), & ! 20
        qi_ham(1,a), & ! 21
        qi_ham(2,a), & ! 22
        qi_ham(3,a)    ! 23
    ENDDO

    CLOSE( UNIT= unit_qi )

    PRINT *, "   quality indicators printed to file ", namefile
    PRINT *


  END PROCEDURE compute_and_print_quality_indicators


END SUBMODULE quality_indicators

