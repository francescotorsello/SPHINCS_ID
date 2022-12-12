! File:         test_lorentz_group.f90
! Author:       Francesco Torsello (FT)
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

PROGRAM test_lorentz_group

  !*****************************************************
  !
  !# Tests the lorentz_group module
  !
  !  FT 12.12.2022
  !
  !*****************************************************

  USE lorentz_group,    ONLY: lorentz_boost, spatial_rotation, &
                              minkowski_sqnorm, eta

  USE constants,  ONLY: pi

  IMPLICIT NONE

  TYPE(lorentz_boost)   :: boost
  TYPE(spatial_rotation):: rotation
  DOUBLE PRECISION, DIMENSION(3):: v
  DOUBLE PRECISION, DIMENSION(4):: u, boosted_u
  DOUBLE PRECISION, DIMENSION(4,4):: eta_mat(0:3,0:3)
  DOUBLE PRECISION, DIMENSION(4,4):: boosted_eta_mat(0:3,0:3)
  DOUBLE PRECISION, DIMENSION(10):: boosted_eta_vec

  !---------------------------!
  !--  End of declarations  --!
  !---------------------------!

  !PRINT *, lorene2hydrobase
  !PRINT *, 2.45191D-4/lorene2hydrobase/1000
  !PRINT *, LOG10(2.45191D-4/lorene2hydrobase/1000)
  !PRINT *
  !PRINT *, LOG10(10**(34.616)/c_light2)
  !STOP

  !PRINT *, "** Polytropic constant used for gamma= 2.75 single polytrope:"
  !PRINT *, "   k used in LORENE= ", 0.01691726009823966
  !PRINT *, "   k converted in SPHINCS units= ", &
  !                               0.01691726009823966*k_lorene2hydrobase(2.75D0)
  !PRINT *
  !PRINT *, "** Polytropic constant used for gamma= 2 single polytrope:"
  !PRINT *, "   k used in LORENE= ", 0.02686965902663748
  !PRINT *, "   k converted in SPHINCS units= ", &
  !                               0.02686965902663748*k_lorene2hydrobase(2.0D0)
  !PRINT *
  !PRINT *, "** Polytropic constant used for the crust in PWP:"
  !PRINT *, "   k used in LORENE= ", 3.99874D-8
  !PRINT *, "   k converted in SPHINCS units= ", &
  !                3.99874D-8*k_lorene2hydrobase_piecewisepolytrope(1.35692395D0)
  !PRINT *
  !PRINT *, "** Polytropic constant used for the crust in PWP:"
  !PRINT *, "   k used in LORENE= ", 8.948185D-2
  !PRINT *, "   k converted in SPHINCS units= ", &
  !                               8.948185D-2*k_lorene2hydrobase(1.35692395D0)
  !PRINT *
  !PRINT *, "   k used in LORENE, corresponding to k-100 in SPHINCS units= ", &
  !         100/k_lorene2hydrobase(2.0D0)
  ! Our testbed cases are gamma= 2.75, k= 30000; and gamma=2, k= 100
  ! in SPHINCS units
  ! 7.901e+14 density for 1.4 GRAVITATIONAL mass, poly 2
  ! 1.4-1.4 systems for both ; 1.6-1.6 ; 1.2-1.8 GRAVIATIONAL masses
  !STOP

  !PRINT *, 1.283004487272563D54*amu/(MSun_geo*km2m*m2cm)**3
  !PRINT *, ( 661708760715581.D0 - 661747751578110.D0 )/661747751578110.D0
  !PRINT *, ( 664672071917413.D0 - 661747751578110.D0 )/661747751578110.D0
  !STOP

  eta_mat(:,0)= [-1.D0,0.D0,0.D0,0.D0]
  eta_mat(:,1)= [0.D0,1.D0,0.D0,0.D0]
  eta_mat(:,2)= [0.D0,0.D0,1.D0,0.D0]
  eta_mat(:,3)= [0.D0,0.D0,0.D0,1.D0]

  v= [0.56D0,0.56D0,0.56D0]
  boost= lorentz_boost(0.56D0,0.56D0,0.56D0)

  u= [1.D0,1.D0,0.D0,0.D0]
  boosted_u= boost% apply_to_vector(u)
  boosted_eta_mat= boost% apply_as_congruence(eta_mat)
  boosted_eta_vec= boost% apply_as_congruence(eta)
  PRINT *, "boosted u=", boosted_u
  PRINT *
  PRINT *, "norm of u=", minkowski_sqnorm(u)
  PRINT *, "norm of boosted u=", minkowski_sqnorm(boosted_u)
  PRINT *
  PRINT *, "boosted eta_mat(0,:)=", boosted_eta_mat(0,:)
  PRINT *, "boosted eta_mat(1,:)=", boosted_eta_mat(1,:)
  PRINT *, "boosted eta_mat(2,:)=", boosted_eta_mat(2,:)
  PRINT *, "boosted eta_mat(3,:)=", boosted_eta_mat(3,:)
  PRINT *
  PRINT *, "boosted eta_vec=", boosted_eta_vec
  PRINT *

  v= [pi/4.D0, pi/3.D0, pi/2.D0]
  rotation= spatial_rotation(pi/4.D0, pi/3.D0, pi/2.D0)

  u= [1.D0,1.D0,0.D0,0.D0]
  boosted_u= rotation% apply_to_vector(u)
  boosted_eta_mat= rotation% apply_as_congruence(eta_mat)
  boosted_eta_vec= rotation% apply_as_congruence(eta)
  PRINT *, "rotated u=", boosted_u
  PRINT *
  PRINT *, "norm of u=", minkowski_sqnorm(u)
  PRINT *, "norm of rotated u=", minkowski_sqnorm(boosted_u)
  PRINT *
  PRINT *, "rotated eta_mat(0,:)=", boosted_eta_mat(0,:)
  PRINT *, "rotated eta_mat(1,:)=", boosted_eta_mat(1,:)
  PRINT *, "rotated eta_mat(2,:)=", boosted_eta_mat(2,:)
  PRINT *, "rotated eta_mat(3,:)=", boosted_eta_mat(3,:)
  PRINT *
  PRINT *, "rotated eta_vec=", boosted_eta_vec
  PRINT *

END PROGRAM test_lorentz_group
