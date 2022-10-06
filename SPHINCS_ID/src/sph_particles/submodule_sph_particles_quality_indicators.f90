! File:         submodule_sph_particles_quality_indicators.f90
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

MODULE quality_indicators

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
  !  eqs.(64) and (67)
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


  SUBROUTINE compute_interpolation_qi1(npart, pos, h, nu, nlrf, qi)

    !****************************************************
    !
    !#
    !
    !  FT 05.10.2022
    !
    !****************************************************

    USE kernel_table, ONLY: &!dWdv_no_norm, &
                            dv_table, dv_table_1,&
                            W_no_norm, n_tab_entry, &
                            interp_gradW_table,interp_W_gradW_table
    USE sphincs_sph,  ONLY: ncand, all_clists
    USE RCB_tree_3D,  ONLY: iorig, nic, nfinal, nprev, lpart, &
                            rpart
    USE utility,      ONLY: zero, one, two, sph_path

    IMPLICIT NONE

    INTEGER, INTENT(IN):: npart

    DOUBLE PRECISION, DIMENSION(3,npart), INTENT(IN)   :: pos
    DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN)   :: h
    DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN)   :: nu
    DOUBLE PRECISION, DIMENSION(npart),   INTENT(IN)   :: nlrf
    DOUBLE PRECISION, DIMENSION(npart),   INTENT(INOUT):: qi
    DOUBLE PRECISION, DIMENSION(npart):: qi2, qi3, qi4, ham_qi

    INTEGER:: a, b, ill, l, itot, inde, index1, k

    DOUBLE PRECISION:: ha, ha_1, ha_3, ha_4, va, xa, ya, za, &
                       hb, hb_1, hb_3, hb_4, xb, yb, zb, rab, rab2, rab_1, &
                       ha2, ha2_4, hb2_4, dx, dy, dz
    !DOUBLE PRECISION:: mat_xx, mat_xy, mat_xz, mat_yy, mat_yz, mat_zz
    DOUBLE PRECISION:: &!Wdx, Wdy, Wdz, ddx, ddy, ddz, Wab, &
                       Wab_ha, Wi, Wi1, dvv, &
                       grW_ha_x, grW_ha_y, grW_ha_z, &
                       grW_hb_x, grW_hb_y, grW_hb_z, eab(3), &
                       Wa, grW, grWa, grWb, vb

    LOGICAL, PARAMETER:: debug= .TRUE.

    INTEGER                      :: ios, unit_qi
    LOGICAL                      :: exist
    CHARACTER(LEN=:), ALLOCATABLE:: namefile, err_msg

    qi    = zero
    qi2   = zero
    qi3   = zero
    qi4   = zero
    ham_qi= zero
    !$OMP PARALLEL DO DEFAULT( NONE ) &
    !$OMP             SHARED( nfinal, nprev, iorig, lpart, rpart, nic,nu,nlrf, &
    !$OMP                     ncand, all_clists, W_no_norm, pos, h, qi, qi2, &
    !$OMP                     qi3, qi4, ham_qi ) &
    !$OMP             PRIVATE( ill, itot, a, b, l, &
    !$OMP                      ha, ha_1, ha_3, ha_4, ha2, ha2_4, &
    !$OMP                      hb, hb_1, hb_3, hb_4, hb2_4, &
    !$OMP                      xa, ya, za, xb, yb, zb, dx, dy, dz, rab2, &
    !$OMP                      inde, va, dv_table_1, index1, Wi, &
    !$OMP                      dvv, dv_table, Wab_ha, Wa, grW, Wi1, &
    !$OMP                      grW_ha_x, grW_ha_y, grW_ha_z, grWa, grWb, &
    !$OMP                      grW_hb_x, grW_hb_y, grW_hb_z, eab, rab, vb, &
    !$OMP                      rab_1 )
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

          IF(debug .AND. .NOT.ALLOCATED(all_clists))THEN
            PRINT *, "all_clists is not allocated"
            PRINT *, "nfinal"
            STOP
          ENDIF
          IF(debug .AND. .NOT.ALLOCATED(all_clists(ill)% list))THEN
            PRINT *, "list is not allocated"
            PRINT *, "ill=", ill
            PRINT *, "nfinal"
            STOP
          ENDIF

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

    !      grWb= grW_hb_x*eab(1) + &
    !            grW_hb_y*eab(2) + &
    !            grW_hb_z*eab(3)

          qi(a)    = qi(a)  + Wab_ha*(nu(b)/nlrf(b))!/hb_3
          qi2(a)   = qi2(a) + (-dx)*Wab_ha*(nu(b)/nlrf(b))!/hb_3
          qi3(a)   = qi3(a) + grW_ha_x*(nu(b)/nlrf(b))
          qi4(a)   = qi4(a) + (-dx)*grW_ha_x*(nu(b)/nlrf(b))
          ham_qi(a)= ham_qi(a) + nu(b)*grW_ha_x*(one/nlrf(a) + nlrf(a)/nlrf(b)**2)

        ENDDO cand_loop

        !qi(a)= qi(a)/(nu(a)/nlrf(a))

      ENDDO particle_loop

    ENDDO ll_cell_loop
    !$OMP END PARALLEL DO

    PRINT *
    PRINT *, " * Printing the interpolation quality indicator Q_1 to file..."

    namefile= TRIM(sph_path)//"quality_indicators.dat"
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
        a, &
        pos(1, a), &
        pos(2, a), &
        pos(3, a), &
        qi(a), qi2(a), qi3(a), qi4(a), ham_qi(a)
    ENDDO

    CLOSE( UNIT= unit_qi )

    PRINT *, "   interpolation quality indicator Q_1 printed to file ", namefile
    PRINT *


  END SUBROUTINE compute_interpolation_qi1


  SUBROUTINE compute_interpolation_qi2

    !****************************************************
    !
    !#
    !
    !  FT 05.10.2022
    !
    !****************************************************

    IMPLICIT NONE

  END SUBROUTINE compute_interpolation_qi2


  SUBROUTINE compute_interpolation_qi3

    !****************************************************
    !
    !#
    !
    !  FT 05.10.2022
    !
    !****************************************************

    IMPLICIT NONE

  END SUBROUTINE compute_interpolation_qi3


  SUBROUTINE compute_gradient_qi4

    !****************************************************
    !
    !#
    !
    !  FT 05.10.2022
    !
    !****************************************************

    IMPLICIT NONE

  END SUBROUTINE compute_gradient_qi4


END MODULE quality_indicators

