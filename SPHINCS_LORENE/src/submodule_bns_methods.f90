! File:         submodule_bns_methods.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (bns_id) bns_methods

  !*****************************************************
  !                                                    *
  ! Implementation of the methods of TYPE bns          *
  !                                                    *
  ! FT 23.10.2020                                      *
  !                                                    *
  !*****************************************************


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE import_id_int

    !*****************************************************
    !                                                    *
    ! Store the ID in the bns arrays                     *
    !                                                    *
    ! FT 5.10.2020                                       *
    !                                                    *
    !*****************************************************

    USE constants, ONLY: Msun_geo

    IMPLICIT NONE

    IF ( C_ASSOCIATED( THIS% bns_ptr ) ) THEN

      IF( SIZE( x ) /= SIZE( y ) .OR. SIZE( x ) /= SIZE( z ) &
            .OR. SIZE( y ) /= SIZE( z ) )THEN
        PRINT *, "** ERROR: The sizes of the arrays of positions" &
                 // "passed to import_lorene_id are not the same."
        PRINT *
        STOP
      ENDIF

      IF( ALLOCATED( THIS% lapse )   .AND. &
          ALLOCATED( THIS% shift_x ) .AND. &
          ALLOCATED( THIS% shift_y ) .AND. &
          ALLOCATED( THIS% shift_z ) .AND. &
          ALLOCATED( THIS% g_xx ) .AND. ALLOCATED( THIS% g_xy ) .AND. &
          ALLOCATED( THIS% g_xz ) .AND. ALLOCATED( THIS% g_yy ) .AND. &
          ALLOCATED( THIS% g_yz ) .AND. ALLOCATED( THIS% g_zz ) .AND. &
          ALLOCATED( THIS% k_xx ) .AND. ALLOCATED( THIS% k_xy ) .AND. &
          ALLOCATED( THIS% k_xz ) .AND. ALLOCATED( THIS% k_yy ) .AND. &
          ALLOCATED( THIS% k_yz ) .AND. ALLOCATED( THIS% k_zz ) .AND. &
          ALLOCATED( THIS% baryon_density )  .AND. &
          ALLOCATED( THIS% energy_density )  .AND. &
          ALLOCATED( THIS% specific_energy ) .AND. &
          ALLOCATED( THIS% v_euler_x )       .AND. &
          ALLOCATED( THIS% v_euler_y )       .AND. &
          ALLOCATED( THIS% v_euler_z ) &
      )THEN

        import_id_loop: DO itr= 1, n, 1

          ! The coordinates need to be converted from SPHINCS units (Msun_geo)
          ! to LORENE units (km). See MODULE constants for the definition of
          ! Msun_geo
          CALL get_lorene_id( THIS% bns_ptr, &
                              x( itr )*Msun_geo, &
                              y( itr )*Msun_geo, &
                              z( itr )*Msun_geo, &
                              THIS% lapse( itr ), &
                              THIS% shift_x( itr ), &
                              THIS% shift_y( itr ), &
                              THIS% shift_z( itr ), &
                              THIS% g_xx( itr ), &
                              THIS% k_xx( itr ), &
                              THIS% k_xy( itr ), &
                              THIS% k_xz( itr ), &
                              THIS% k_yy( itr ), &
                              THIS% k_yz( itr ), &
                              THIS% k_zz( itr ), &
                              THIS% baryon_density( itr ), &
                              THIS% energy_density( itr ), &
                              THIS% specific_energy( itr ), &
                              THIS% v_euler_x( itr ), &
                              THIS% v_euler_y( itr ), &
                              THIS% v_euler_z( itr ) )

          !
          !-- The following follows from the assumption of conformal
          !-- flatness in LORENE
          !
          THIS% g_yy( itr )= THIS% g_xx( itr )
          THIS% g_zz( itr )= THIS% g_xx( itr )
          THIS% g_xy( itr )= 0.0D0
          THIS% g_xz( itr )= 0.0D0
          THIS% g_yz( itr )= 0.0D0

          !
          !- Set/unset the geodesic gauge
          !
          IF( THIS% one_lapse )THEN
            THIS% lapse( itr )= 1.0D0
          ENDIF
          IF( THIS% zero_shift )THEN
            THIS% shift_x( itr )= 0.0D0
            THIS% shift_y( itr )= 0.0D0
            THIS% shift_z( itr )= 0.0D0
          ENDIF

          !
          !-- Convert the extrinsic curvature from LORENE units to
          !-- SPHINCS units
          !
          THIS% k_xx( itr )= THIS% k_xx( itr )*Msun_geo
          THIS% k_xy( itr )= THIS% k_xy( itr )*Msun_geo
          THIS% k_xz( itr )= THIS% k_xz( itr )*Msun_geo
          THIS% k_yy( itr )= THIS% k_yy( itr )*Msun_geo
          THIS% k_yz( itr )= THIS% k_yz( itr )*Msun_geo
          THIS% k_zz( itr )= THIS% k_zz( itr )*Msun_geo

          ! Print progress on screen
          perc= 100*itr/n
          IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
            WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
                    creturn//" ", perc, "%"
          ENDIF

        ENDDO import_id_loop
        IF( show_progress ) WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

      ELSE

        PRINT *, "** ERROR: Memory was not allocated before calling " &
                 // "import_id in import_lorene_id (TYPE particles)."
        PRINT *
        STOP

      ENDIF

      PRINT *, "** Subroutine import_lorene_id executed."
      PRINT *

    ENDIF

  END PROCEDURE import_id_int


  MODULE PROCEDURE import_id_ext

    !*****************************************************
    !                                                    *
    ! Store the ID in non-member arrays with             *
    ! the same shape as the bns arrays                   *
    !                                                    *
    ! FT 5.10.2020                                       *
    !                                                    *
    !*****************************************************

    USE constants, ONLY: Msun_geo

    IMPLICIT NONE

    IF ( C_ASSOCIATED( THIS% bns_ptr ) ) THEN

      IF( SIZE( x ) /= SIZE( y ) .OR. SIZE( x ) /= SIZE( z ) &
            .OR. SIZE( y ) /= SIZE( z ) )THEN
        PRINT *, "** ERROR: The sizes of the arrays of positions" &
                 // "passed to import_lorene_id are not the same."
        PRINT *
        STOP
      ENDIF

      import_id_loop: DO itr= 1, n, 1

        ! The coordinates need to be converted from SPHINCS units (Msun_geo)
        ! to LORENE units (km). See MODULE constants for the definition of
        ! Msun_geo
        CALL get_lorene_id( THIS% bns_ptr, &
                            x( itr )*Msun_geo, &
                            y( itr )*Msun_geo, &
                            z( itr )*Msun_geo, &
                            lapse( itr ), &
                            shift_x( itr ), &
                            shift_y( itr ), &
                            shift_z( itr ), &
                            g_xx( itr ), &
                            k_xx( itr ), &
                            k_xy( itr ), &
                            k_xz( itr ), &
                            k_yy( itr ), &
                            k_yz( itr ), &
                            k_zz( itr ), &
                            baryon_density( itr ), &
                            energy_density( itr ), &
                            specific_energy( itr ), &
                            u_euler_x( itr ), &
                            u_euler_y( itr ), &
                            u_euler_z( itr ) )

        !
        !-- The following follows from the assumption of conformal
        !-- flatness in LORENE
        !
        g_yy( itr )= g_xx( itr )
        g_zz( itr )= g_xx( itr )
        g_xy( itr )= 0.0D0
        g_xz( itr )= 0.0D0
        g_yz( itr )= 0.0D0

        !
        !- Set/unset the geodesic gauge
        !
        IF( THIS% one_lapse )THEN
          lapse( itr )= 1.0D0
        ENDIF
        IF( THIS% zero_shift )THEN
          shift_x( itr )= 0.0D0
          shift_y( itr )= 0.0D0
          shift_z( itr )= 0.0D0
        ENDIF

        !
        !-- Convert the extrinsic curvature from LORENE units to
        !-- SPHINCS units
        !
        k_xx( itr )= k_xx( itr )*Msun_geo
        k_xy( itr )= k_xy( itr )*Msun_geo
        k_xz( itr )= k_xz( itr )*Msun_geo
        k_yy( itr )= k_yy( itr )*Msun_geo
        k_yz( itr )= k_yz( itr )*Msun_geo
        k_zz( itr )= k_zz( itr )*Msun_geo

        ! Print progress on screen
        perc= 100*itr/n
        IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
          WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
                  creturn//" ", perc, "%"
        ENDIF

      ENDDO import_id_loop
      IF( show_progress ) WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

      PRINT *, "** Subroutine import_lorene_id executed."
      PRINT *

    ENDIF

  END PROCEDURE import_id_ext


  MODULE PROCEDURE import_id_multid_array

    !*****************************************************
    !                                                    *
    ! Store the spacetime ID in multi-dimensional arrays *
    ! needed for the spacetime part of SPHINCS           *
    !                                                    *
    ! FT 22.11.2020                                      *
    !                                                    *
    !*****************************************************

    USE constants, ONLY: Msun_geo

    USE tensor,    ONLY: itt, itx, ity, itz, ixx, ixy, &
                         ixz, iyy, iyz, izz, jxx, jxy, jxz, &
                         jyy, jyz, jzz, jx, jy, jz, n_sym4x4

    IMPLICIT NONE

    INTEGER:: ix, iy, iz

    DOUBLE PRECISION:: detg
    DOUBLE PRECISION:: detg4
    DOUBLE PRECISION, DIMENSION( :, :, :, : ), ALLOCATABLE:: g4

    ! g4 is allocatable to allocate it on the heap
    ! Allocating it on the stack might exceed stack memory,
    ! causing a segmentation fault
    ALLOCATE( g4( nx, ny, nz, n_sym4x4 ) )

    IF ( C_ASSOCIATED( THIS% bns_ptr ) ) THEN

      IF( .FALSE. &!SHAPE( pos(:,:,:,1) ) /= SHAPE( lapse ) .OR. &
          !SHAPE( pos(:,:,:,1) ) /= SHAPE( shift(:,:,:,jx) ) & ! .OR. &
        ! SHAPE( pos(:,:,:,1) ) /= SHAPE( g(:,:,:,1) ) .OR. &
        ! SHAPE( pos(:,:,:,1) ) /= SHAPE( k(:,:,:,1) ) &
        )THEN
        PRINT *, "** ERROR: Mismatch in array dimensions" &
                 // "in import_id_multid_array."
        PRINT *
        STOP
      ENDIF

      coords_z: DO iz= 1, nz, 1
        coords_y: DO iy= 1, ny, 1
          coords_x: DO ix= 1, nx, 1

            ! The coordinates need to be converted from SPHINCS units (Msun_geo)
            ! to LORENE units (km). See MODULE constants for the definition of
            ! Msun_geo
            CALL get_lorene_id_spacetime( THIS% bns_ptr, &
                                pos( ix, iy, iz, jx )*Msun_geo, &
                                pos( ix, iy, iz, jy )*Msun_geo, &
                                pos( ix, iy, iz, jz )*Msun_geo, &
                                lapse( ix, iy, iz ), &
                                shift( ix, iy, iz, jx ), &
                                shift( ix, iy, iz, jy ), &
                                shift( ix, iy, iz, jz ), &
                                g( ix, iy, iz, jxx ), &
                                k( ix, iy, iz, jxx ), &
                                k( ix, iy, iz, jxy ), &
                                k( ix, iy, iz, jxz ), &
                                k( ix, iy, iz, jyy ), &
                                k( ix, iy, iz, jyz ), &
                                k( ix, iy, iz, jzz ) )

            !
            !-- The following follows from the assumption of
            !-- conformal flatness in LORENE
            !
            g( ix, iy, iz, jyy )= g( ix, iy, iz, jxx )
            g( ix, iy, iz, jzz )= g( ix, iy, iz, jxx )
            g( ix, iy, iz, jxy )= 0.0D0
            g( ix, iy, iz, jxz )= 0.0D0
            g( ix, iy, iz, jyz )= 0.0D0

            !
            !- Set/unset the geodesic gauge
            !
            IF( THIS% one_lapse )THEN
              lapse( ix, iy, iz )= 1.0D0
            ENDIF
            IF( THIS% zero_shift )THEN
              shift( ix, iy, iz, jx )= 0.0D0
              shift( ix, iy, iz, jy )= 0.0D0
              shift( ix, iy, iz, jz )= 0.0D0
            ENDIF

            !
            !-- Convert the extrinsic curvature from LORENE units to
            !-- SPHINCS units
            !
            k( ix, iy, iz, jxx )= k( ix, iy, iz, jxx )*Msun_geo
            k( ix, iy, iz, jxy )= k( ix, iy, iz, jxy )*Msun_geo
            k( ix, iy, iz, jxz )= k( ix, iy, iz, jxz )*Msun_geo
            k( ix, iy, iz, jyy )= k( ix, iy, iz, jyy )*Msun_geo
            k( ix, iy, iz, jyz )= k( ix, iy, iz, jyz )*Msun_geo
            k( ix, iy, iz, jzz )= k( ix, iy, iz, jzz )*Msun_geo

!IF( ix < 10 .AND. iy == 1 .AND. iz == 1 )THEN
!                PRINT *, "5"
!              ENDIF

            !IF( ix == 1 .AND. iy == 1 .AND. iz == 54 )THEN
            !  PRINT *, "g_xx=", g(1,1,54,jxx)
            !  PRINT *, "g_xy=", g(1,1,54,jxy)
            !  PRINT *, "g_xz=", g(1,1,54,jxz)
            !  PRINT *, "g_yy=", g(1,1,54,jyy)
            !  PRINT *, "g_yz=", g(1,1,54,jyz)
            !  PRINT *, "g_zz=", g(1,1,54,jzz)
            !  PRINT *
            !ENDIF
            !IF( ix == 1 .AND. iy == 1 .AND. iz == 55 )THEN
            !  PRINT *, "g_xx=", g(1,1,55,jxx)
            !  PRINT *, "g_xy=", g(1,1,55,jxy)
            !  PRINT *, "g_xz=", g(1,1,55,jxz)
            !  PRINT *, "g_yy=", g(1,1,55,jyy)
            !  PRINT *, "g_yz=", g(1,1,55,jyz)
            !  PRINT *, "g_zz=", g(1,1,55,jzz)
            !  PRINT *
            !ENDIF

            detg= 2.0D0*g(ix,iy,iz,jxy)*g(ix,iy,iz,jxz)*g(ix,iy,iz,jyz) &
                  - g(ix,iy,iz,jzz)*g(ix,iy,iz,jxy)**2 + g(ix,iy,iz,jyy) &
                   *( g(ix,iy,iz,jxx)*g(ix,iy,iz,jzz) - g(ix,iy,iz,jxz)**2 ) &
                  - g(ix,iy,iz,jxx)*g(ix,iy,iz,jyz)**2

!IF( ix < 10 .AND. iy == 1 .AND. iz == 1 )THEN
!                PRINT *, "6"
!              ENDIF

            IF( ABS( detg ) < 1D-10 )THEN
              PRINT *, "The determinant of the spatial metric " &
                       // "is effectively 0 at the grid point " &
                       // "(ix,iy,iz)= (", ix, ",", iy,",",iz, ")."
              PRINT *, "detg=", detg
              PRINT *
              STOP
            ELSEIF( detg < 0 )THEN
              PRINT *, "The determinant of the spatial metric " &
                       // "is negative at the grid point " &
                       // "(ix,iy,iz)= (", ix, ",", iy,",",iz, ")."
              PRINT *, "detg=", detg
              PRINT *
              STOP
            ENDIF

!IF( ix < 10 .AND. iy == 1 .AND. iz == 1 )THEN
!  PRINT *, "7"
!ENDIF

            CALL compute_g4( ix, iy, iz, lapse, shift, g, g4 )

!IF( ix < 10 .AND. iy == 1 .AND. iz == 1 )THEN
!  PRINT *, "8"
!ENDIF

            CALL determinant_sym4x4_grid( ix, iy, iz, g4, detg4 )

!IF( ix < 10 .AND. iy == 1 .AND. iz == 1 )THEN
!  PRINT *, "9"
!ENDIF

            IF( ABS( detg4 ) < 1D-10 )THEN
              PRINT *, "The determinant of the spacetime metric "&
                       // "is effectively 0 at the grid point " &
                       // "(ix,iy,iz)= (", ix, ",", iy,",",iz, ")."
              PRINT *, "detg4=", detg4
              PRINT *
              STOP
            ELSEIF( detg4 > 0 )THEN
              PRINT *, "The determinant of the spacetime metric "&
                       // "is positive at the grid point " &
                       // "(ix,iy,iz)= (", ix, ",", iy,",",iz, ")."
              PRINT *, "detg4=", detg4
              PRINT *
              STOP
            ENDIF

!IF( ix < 10 .AND. iy == 1 .AND. iz == 1 )THEN
!  PRINT *, "10"
!ENDIF
!IF( ix < 20 .AND. iy == 1 .AND. iz == 1 )THEN
!  PRINT *, "ix=", ix, ", iy=", iy, ", iz=", iz
!ENDIF

            ! Print progress on screen
            perc= 100*( nx*ny*(iz - 1) + nx*(iy - 1) + ix )/( nx*ny*nz )
            !perc2= 100.0*DBLE(nx*ny*(iz - 1) + nx*(iy - 1) + ix)/DBLE( nx*ny*nz )
            !perc= 100*cnt/( nx*ny*nz )
            IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
              WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
                      creturn//" ", perc, "%"
              !WRITE( *, "(A2,F5.2,A1)", ADVANCE= "NO" ) &
              !        creturn//" ", perc2, "%"
            ENDIF

          ENDDO coords_x
        ENDDO coords_y
      ENDDO coords_z
      IF( show_progress ) WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

      PRINT *, "** Subroutine import_lorene_id executed."
      PRINT *

    ENDIF

  END PROCEDURE import_id_multid_array


  MODULE PROCEDURE import_id_hydro

    !*****************************************************
    !                                                    *
    ! Store the hydro ID in the arrays needed            *
    ! for the spacetime part of SPHINCS                  *
    !                                                    *
    ! FT 25.11.2020                                      *
    !                                                    *
    !*****************************************************

    USE constants,  ONLY: Msun_geo
    USE tensor,     ONLY: jx, jy, jz

    IMPLICIT NONE

    INTEGER:: ix, iy, iz

    IF ( C_ASSOCIATED( THIS% bns_ptr ) ) THEN

      coords_z: DO iz= 1, nz, 1
        coords_y: DO iy= 1, ny, 1
          coords_x: DO ix= 1, nx, 1

            ! The coordinates need to be converted from SPHINCS units (Msun_geo)
            ! to LORENE units (km). See MODULE constants for the definition of
            ! Msun_geo
            CALL get_lorene_id_hydro( THIS% bns_ptr, &
                              pos( ix, iy, iz, jx )*Msun_geo, &
                              pos( ix, iy, iz, jy )*Msun_geo, &
                              pos( ix, iy, iz, jz )*Msun_geo, &
                              baryon_density( ix, iy, iz ), &
                              energy_density( ix, iy, iz ), &
                              specific_energy( ix, iy, iz ), &
                              pressure( ix, iy, iz ), &
                              u_euler( ix, iy, iz, 1 ), &
                              u_euler( ix, iy, iz, 2 ), &
                              u_euler( ix, iy, iz, 3 ) )

            ! Print progress on screen
            perc= 100*(nx*ny*(iz - 1) &
                  + nx*(iy - 1) + ix)/( nx*ny*nz )
            IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
              WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
                      creturn//" ", perc, "%"
            ENDIF

          ENDDO coords_x
        ENDDO coords_y
      ENDDO coords_z
      IF( show_progress ) WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

      PRINT *, "** Subroutine import_lorene_id_hydro executed."
      PRINT *

    ENDIF

  END PROCEDURE import_id_hydro


  MODULE PROCEDURE import_id_particles

    !*****************************************************
    !                                                    *
    ! Import the LORENE ID on the particles, ignoring    *
    ! the extrinsic curvature which is not needed.       *
    !                                                    *
    ! FT 19.11.2020                                      *
    !                                                    *
    !*****************************************************

    USE constants, ONLY: Msun_geo

    IMPLICIT NONE

    DOUBLE PRECISION:: detg

    IF ( C_ASSOCIATED( THIS% bns_ptr ) ) THEN

      IF( SIZE( x ) /= SIZE( y ) .OR. SIZE( x ) /= SIZE( z ) &
              .OR. SIZE( y ) /= SIZE( z ) )THEN
        PRINT *, "** ERROR: The sizes of the arrays of positions" &
                 // "passed to import_lorene_id are not the same."
        PRINT *
        STOP
      ENDIF

      PRINT *, "** Importing ID on particles..."

      import_id_loop: DO itr= 1, n, 1

        ! The coordinates need to be converted from SPHINCS units (Msun_geo)
        ! to LORENE units (km). See MODULE constants for the definition of
        ! Msun_geo
        CALL get_lorene_id_particles( THIS% bns_ptr, &
                                      x( itr )*Msun_geo, &
                                      y( itr )*Msun_geo, &
                                      z( itr )*Msun_geo, &
                                      lapse( itr ), &
                                      shift_x( itr ), &
                                      shift_y( itr ), &
                                      shift_z( itr ), &
                                      g_xx( itr ), &
                                      baryon_density( itr ), &
                                      energy_density( itr ), &
                                      specific_energy( itr ), &
                                      pressure( itr ), &
                                      u_euler_x( itr ), &
                                      u_euler_y( itr ), &
                                      u_euler_z( itr ) )

        !
        !-- The following follows from the assumption of conformal
        !-- flatness in LORENE
        !
        g_yy( itr )= g_xx( itr )
        g_zz( itr )= g_xx( itr )
        g_xy( itr )= 0.0D0
        g_xz( itr )= 0.0D0
        g_yz( itr )= 0.0D0

        !
        !- Set/unset the geodesic gauge
        !
        IF( THIS% one_lapse )THEN
          lapse( itr )= 1.0D0
        ENDIF
        IF( THIS% zero_shift )THEN
          shift_x( itr )= 0.0D0
          shift_y( itr )= 0.0D0
          shift_z( itr )= 0.0D0
        ENDIF

        detg= 2*g_xy(itr)*g_xz(itr)*g_yz(itr) &
              - g_zz(itr)*g_xy(itr)**2 &
              + g_yy(itr)*( g_xx(itr)*g_zz(itr) - g_xz(itr)**2 ) &
              - g_xx(itr)*g_yz(itr)**2

        IF( ABS( detg ) < 1D-10 )THEN
          PRINT *, "The determinant of the spatial metric is " &
                   // "effectively 0 at the particle ", itr
          PRINT *, "detg=", detg
          PRINT *
          STOP
        ELSEIF( detg < 0 )THEN
          PRINT *, "The determinant of the spatial metric is " &
                   // "negative at the particle ", itr
          PRINT *, "detg=", detg
          PRINT *
          STOP
        ENDIF

        ! Print progress on screen
        perc= 100*itr/n
        IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
          WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
                  creturn//" ", perc, "%"
        ENDIF

      ENDDO import_id_loop
      IF( show_progress ) WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

      PRINT *, "** Subroutine import_id_particles executed."
      PRINT *

    ENDIF

  END PROCEDURE import_id_particles


  MODULE PROCEDURE import_id_mass_b

    !*****************************************************
    !                                                    *
    ! Import the LORENE ID needed to compute the baryon  *
    ! mass, storing it to variables (not arrays as the   *
    ! others importing SUBROUTINES).                     *
    !                                                    *
    ! FT 15.04.2021                                      *
    !                                                    *
    !*****************************************************

    USE constants, ONLY: Msun_geo, lorene2hydrobase

    IMPLICIT NONE

    DOUBLE PRECISION:: detg

    IF ( C_ASSOCIATED( THIS% bns_ptr ) ) THEN

      ! The coordinates need to be converted from SPHINCS units (Msun_geo)
      ! to LORENE units (km). See MODULE constants for the definition of
      ! Msun_geo
      CALL get_lorene_id_mass_b( THIS% bns_ptr, &
                                    x*Msun_geo, &
                                    y*Msun_geo, &
                                    z*Msun_geo, &
                                    g_xx, &
                                    baryon_density, &
                                    gamma_euler )

      baryon_density= baryon_density*lorene2hydrobase

    ENDIF

  END PROCEDURE import_id_mass_b


  MODULE PROCEDURE import_id_k

    !*****************************************************
    !                                                    *
    ! Store the extrinsic curvature in the               *
    ! arrays needed for the SPH part of                  *
    ! SPHINCS                                            *
    !                                                    *
    ! FT 25.11.2020                                      *
    !                                                    *
    !*****************************************************

    USE constants, ONLY: Msun_geo

    IMPLICIT NONE

    DOUBLE PRECISION:: detg

    IF ( C_ASSOCIATED( THIS% bns_ptr ) ) THEN

      IF( SIZE( x ) /= SIZE( y ) .OR. SIZE( x ) /= SIZE( z ) &
              .OR. SIZE( y ) /= SIZE( z ) )THEN
        PRINT *, "** ERROR: The sizes of the arrays of positions" &
                 // "passed to import_lorene_id are not the same."
        PRINT *
        STOP
      ENDIF

      import_id_loop: DO itr= 1, n, 1

        ! The coordinates need to be converted from SPHINCS units (Msun_geo)
        ! to LORENE units (km). See MODULE constants for the definition of
        ! Msun_geo
        CALL get_lorene_id_k( THIS% bns_ptr, &
                              x( itr )*Msun_geo, &
                              y( itr )*Msun_geo, &
                              z( itr )*Msun_geo, &
                              k_xx( itr ), &
                              k_xy( itr ), &
                              k_xz( itr ), &
                              k_yy( itr ), &
                              k_yz( itr ), &
                              k_zz( itr ) )

        !
        !-- Convert the extrinsic curvature from LORENE units to
        !-- SPHINCS units
        !
        k_xx( itr )= k_xx( itr )*Msun_geo
        k_xy( itr )= k_xy( itr )*Msun_geo
        k_xz( itr )= k_xz( itr )*Msun_geo
        k_yy( itr )= k_yy( itr )*Msun_geo
        k_yz( itr )= k_yz( itr )*Msun_geo
        k_zz( itr )= k_zz( itr )*Msun_geo

        ! Print progress on screen
        perc= 100*itr/n
        IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
          WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
                  creturn//" ", perc, "%"
        ENDIF

      ENDDO import_id_loop
      IF( show_progress ) WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

      PRINT *, "** Subroutine import_lorene_id_k executed."
      PRINT *

    ENDIF

  END PROCEDURE import_id_k


  MODULE PROCEDURE print_id_params

    !*****************************************************
    !                                                    *
    ! Print the parameters of the binary neutron         *
    ! stars' initial data computed by LORENE             *
    !                                                    *
    ! FT 8.10.2020                                       *
    !                                                    *
    !*****************************************************

    USE constants, ONLY: k_lorene2hydrobase, Msun_geo, km2m, m2cm, kg2g, &
                         lorene2hydrobase

    IMPLICIT NONE

    IF( THIS% angular_momentum == 0.0D0 )THEN

      PRINT *
      PRINT *, " ** The parameters have not ben read yet. ", &
          "Call the SUBROUTINE import_lorene_id_params to read them."
      PRINT *

    ELSE

      PRINT *
      PRINT *, " ** The parameters of the binary system are:"
      PRINT *
      PRINT *, " Distance between the points of highest density = ",&
               THIS% distance, " M_sun^geo"
      PRINT *, " Distance between the barycenters = ", &
               THIS% distance_com, " M_sun^geo"
      PRINT *
      PRINT *, " Mass of NS 1 = ", THIS% mass1, " M_sun"
      PRINT *, " Mass of NS 2 = ", THIS% mass2, " M_sun"
      PRINT *, " ADM mass = ", THIS% adm_mass, " M_sun"
      PRINT *
      PRINT *, " Stellar center of NS 1 = ", THIS% center1_x, " M_sun^geo"
      PRINT *, " Stellar center of NS 2 = ", THIS% center2_x, " M_sun^geo"
      PRINT *, " Barycenter of NS 1 = ", THIS% barycenter1_x, " M_sun^geo"
      PRINT *, " Barycenter of NS 2 = ", THIS% barycenter2_x, " M_sun^geo"
      PRINT *, " Angular velocity = ", THIS% angular_vel, " rad/s"
      PRINT *, " Angular momentum of the system = ", &
               THIS% angular_momentum, " G M_sun^2 /c"
      PRINT *
      PRINT *, " Radii of star 1: "
      PRINT *, "  x direction, towards companion = ", &
               THIS% radius1_x_comp, " M_sun^geo"
      PRINT *, "  x direction, opposite to companion = ", &
               THIS% radius1_x_opp, " M_sun^geo"
      PRINT *, "  y direction = ", THIS% radius1_y, " M_sun^geo"
      PRINT *, "  z direction = ", THIS% radius1_z, " M_sun^geo"
      PRINT *, " Radii of star 2 :"
      PRINT *, "  x direction, towards companion = ", &
               THIS% radius2_x_comp, " M_sun^geo"
      PRINT *, "  x direction, opposite to companion = ", &
               THIS% radius2_x_opp, " M_sun^geo"
      PRINT *, "  y direction = ", THIS% radius2_y, " M_sun^geo"
      PRINT *, "  z direction = ", THIS% radius2_z, " M_sun^geo"
      PRINT *
      PRINT *, " Hydro quantities at the center of star 1: "
      PRINT *, "  Central enthalpy = ", THIS% ent_center1, " c^2"
      PRINT *, "  Central baryon number density = ", THIS% nbar_center1, &
               " (M_sun^geo)^{-3} =", &
               THIS% nbar_center1/(MSun_geo*km2m*m2cm)**3, "cm^{-3}"
      PRINT *, "  Central baryon mass density = ", THIS% rho_center1, &
               " M_sun^geo (M_sun^geo)^{-3} =", &
               THIS% rho_center1/lorene2hydrobase*kg2g/(m2cm**3), "g cm^{-3}"
      PRINT *, "  Central energy density = ", THIS% energy_density_center1, &
               " M_sun^geo c^2 (M_sun^geo)^{-3}", &
               THIS% energy_density_center1/lorene2hydrobase*kg2g/(m2cm**3), &
               "g c^2 cm^{-3}"
      PRINT *, "  Central specific energy = ", THIS% specific_energy_center1, &
               " c^2"
      PRINT *, "  Central pressure = ", THIS% pressure_center1, &
               " M_sun^geo c^2 (M_sun^geo)^{-3}", &
               THIS% pressure_center1/lorene2hydrobase*kg2g/(m2cm**3), &
               "g c^2 cm^{-3}"
      PRINT *, " Hydro quantities at the center of star 2: "
      PRINT *, "  Central enthalpy = ", THIS% ent_center2, " c^2"
      PRINT *, "  Central baryon number density = ", THIS% nbar_center2, &
               " (M_sun^geo)^{-3} =", &
               THIS% nbar_center2/(MSun_geo*km2m*m2cm)**3, "cm^{-3}"
      PRINT *, "  Central baryon mass density = ", THIS% rho_center2, &
               " M_sun^geo (M_sun^geo)^{-3} =", &
               THIS% rho_center2/lorene2hydrobase*kg2g/(m2cm**3), "g cm^{-3}"
      PRINT *, "  Central energy density = ", THIS% energy_density_center2, &
               " M_sun^geo c^2 (M_sun^geo)^{-3}", &
               THIS% energy_density_center2/lorene2hydrobase*kg2g/(m2cm**3), &
               "g c^2 cm^{-3}"
      PRINT *, "  Central specific energy = ", THIS% specific_energy_center2, &
               " c^2"
      PRINT *, "  Central pressure = ", THIS% pressure_center2, &
               " M_sun^geo c^2 (M_sun^geo)^{-3}", &
               THIS% pressure_center2/lorene2hydrobase*kg2g/(m2cm**3), &
               "g c^2 cm^{-3}"
      PRINT *
      IF( show_progress ) &
        PRINT *, " Equations of state for star 1 (EOS1) = ", TRIM(THIS% eos1)
      IF( show_progress ) &
        PRINT *, " Equations of state for star 2 (EOS2) = ", TRIM(THIS% eos2)
      PRINT *
      STOP

      IF( THIS% gamma0_1 == 0 )THEN ! If the EOS is polytropic

        PRINT *, " Parameters for EOS1: "
        PRINT *, "  Polytopic index gamma_1 = ", THIS% gamma_1
        PRINT *, "  Pressure coefficient = ",&
                 THIS% kappa_1/k_lorene2hydrobase( THIS% gamma_1 ), &
                 "rho_nuc c^2 / n_nuc^gamma_1 = ", THIS% kappa_1, &
                 "[pure number]"
        PRINT *, " Parameters for EOS2: "
        PRINT *, "  Polytopic index gamma_2 = ", THIS% gamma_2
        PRINT *, "  Pressure coefficient = ",&
                 THIS% kappa_2/k_lorene2hydrobase( THIS% gamma_2 ), &
                 "rho_nuc c^2 / n_nuc^gamma_2 = ", THIS% kappa_2, &
                 "[pure number]"
        PRINT *

      ELSEIF( THIS% gamma0_1 /= 0 )THEN ! If the EOS is piecewise polytropic

        PRINT *, " Parameters for EOS1: "
        PRINT *, "  Number of polytropic indexes = ", THIS% npeos_1
        PRINT *, "  Polytopic index gamma0_1 = ", THIS% gamma0_1
        PRINT *, "  Polytopic index gamma1_1 = ", THIS% gamma1_1
        PRINT *, "  Polytopic index gamma2_1 = ", THIS% gamma2_1
        PRINT *, "  Polytopic index gamma3_1 = ", THIS% gamma3_1
        PRINT *, "  Pressure coefficient for the crust (here from SLy) = ",&
                 THIS% kappa0_1/k_lorene2hydrobase( THIS% gamma0_1 ), &
                 "rho_nuc c^2 / n_nuc^gamma0_1 = ", THIS% kappa0_1, &
                 "[pure number]"
        PRINT *, "  Pressure coefficient for the first polytrope = ",&
                 THIS% kappa1_1/k_lorene2hydrobase( THIS% gamma1_1 ), &
                 "rho_nuc c^2 / n_nuc^gamma1_1", THIS% kappa1_1, &
                 "[pure number]"
        PRINT *, "  Pressure coefficient for the second polytrope = ",&
                 THIS% kappa2_1/k_lorene2hydrobase( THIS% gamma2_1 ), &
                 "rho_nuc c^2 / n_nuc^gamma2_1", THIS% kappa2_1, &
                 "[pure number]"
        PRINT *, "  Pressure coefficient for the third polytrope = ",&
                 THIS% kappa3_1/k_lorene2hydrobase( THIS% gamma3_1 ), &
                 "rho_nuc c^2 / n_nuc^gamma3_1", THIS% kappa3_1, &
                 "[pure number]"
        PRINT *, "  Base 10 exponent of the pressure at the first fiducial " &
                 // "density (between gamma_0 and gamma_1) (dyne/cm^2)= ", &
                 THIS% logP1_1
        PRINT *, "  Base 10 exponent of first fiducial density (g/cm^3) = ", &
                 THIS% logRho0_1
        PRINT *, "  Base 10 exponent of second fiducial density (g/cm^3) = ",&
                 THIS% logRho1_1
        PRINT *, "  Base 10 exponent of third fiducial density (g/cm^3) = ", &
                 THIS% logRho2_1
        PRINT *
        PRINT *, " Parameters for EOS2: "
        PRINT *, "  Number of polytropic indexes = ", THIS% npeos_2
        PRINT *, "  Polytopic index gamma0_2 = ", THIS% gamma0_2
        PRINT *, "  Polytopic index gamma1_2 = ", THIS% gamma1_2
        PRINT *, "  Polytopic index gamma2_2 = ", THIS% gamma2_2
        PRINT *, "  Polytopic index gamma3_2 = ", THIS% gamma3_2
        PRINT *, "  Pressure coefficient for the crust (here from SLy) = ",&
                 THIS% kappa0_2/k_lorene2hydrobase( THIS% gamma0_2 ), &
                 "rho_nuc c^2 / n_nuc^gamma0_2 = ", THIS% kappa0_2, &
                 "[pure number]"
        PRINT *, "  Pressure coefficient for the first polytrope = ",&
                 THIS% kappa1_2/k_lorene2hydrobase( THIS% gamma1_2 ), &
                 "rho_nuc c^2 / n_nuc^gamma1_2", THIS% kappa1_2, &
                 "[pure number]"
        PRINT *, "  Pressure coefficient for the second polytrope = ",&
                 THIS% kappa2_2/k_lorene2hydrobase( THIS% gamma2_2 ), &
                 "rho_nuc c^2 / n_nuc^gamma2_2", THIS% kappa2_2, &
                 "[pure number]"
        PRINT *, "  Pressure coefficient for the third polytrope = ",&
                 THIS% kappa3_2/k_lorene2hydrobase( THIS% gamma3_2 ), &
                 "rho_nuc c^2 / n_nuc^gamma3_2", THIS% kappa3_2, &
                 "[pure number]"
        PRINT *, "  Base 10 exponent of the pressure at the first fiducial " &
                 // "density (between gamma_0 and gamma_1) (dyne/cm^2)= ", &
                 THIS% logP1_2
        PRINT *, "  Base 10 exponent of first fiducial density (g/cm^3) = ", &
                 THIS% logRho0_2
        PRINT *, "  Base 10 exponent of second fiducial density (g/cm^3) = ",&
                 THIS% logRho1_2
        PRINT *, "  Base 10 exponent of third fiducial density (g/cm^3) = ", &
                 THIS% logRho2_2
        PRINT *

      ELSE

        PRINT *, "** ERROR in SUBROUTINE import_lorene_id_params!", &
                 " The equation of state is unknown!"
        STOP

      ENDIF

    ENDIF

  END PROCEDURE print_id_params


  MODULE PROCEDURE integrate_baryon_mass_density

    !************************************************
    !                                               *
    ! Perform 3D integration over a spherical grid  *
    ! of the baryon mass density. Output baryon     *
    ! mass and radial mass profile.                 *
    !                                               *
    ! FT 19.02.2021                                 *
    !                                               *
    !************************************************

    !$ USE OMP_LIB
    USE constants, ONLY: pi
    USE NR,        ONLY: indexx

    IMPLICIT NONE

    INTEGER:: r, th, phi
    DOUBLE PRECISION:: rad_coord, colat, long, mass_element
    DOUBLE PRECISION:: g_xx, sq_g, baryon_density, gamma_euler
    !DOUBLE PRECISION:: rad

    !rad= 0.0D0

    IF(.NOT.ALLOCATED( mass_profile ))THEN
      ALLOCATE( mass_profile( 3, 0:NINT(radius/dr) ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array mass_profile in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array pos in SUBROUTINE" &
      !                // "place_particles_3D_lattice." )
    ENDIF
    IF(.NOT.ALLOCATED( mass_profile_idx ))THEN
      ALLOCATE( mass_profile_idx( 0:NINT(radius/dr) ), STAT= ios, &
                ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array mass_profile in SUBROUTINE" &
                  // "place_particles_. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...allocation error for array pos in SUBROUTINE" &
      !                // "place_particles_3D_lattice." )
    ENDIF

    mass_profile( 1, 0 )= 0.0D0
    mass_profile( 2, 0 )= 4.0D0/3.0D0*pi*dr**3.0D0*central_density
    mass_profile( 3, 0 )= 4.0D0/3.0D0*pi*dr**3.0D0*central_density

    !$OMP PARALLEL DO SHARED(dr,dphi,dth,center)&
    !$OMP             PRIVATE(r,th,phi,rad_coord,long,colat,sq_g,gamma_euler, &
    !$OMP                     baryon_density,mass_element,mass) !&
!    !$OMP             SCHEDULE(STATIC,1)
    radius_loop: DO r= 1, NINT(radius/dr), 1

      mass= 0.0D0
      rad_coord= r*dr

!      !$OMP PARALLEL SECTIONS REDUCTION(+:mass,rad)
      !rad= rad + dr

      longitude_loop: DO phi= 1, NINT(2.0D0*pi/dphi), 1

        long= phi*dphi

        colatitude_loop: DO th= 1, NINT(pi/2.0D0/dth), 1

          colat= th*dth

          ! The definition of the baryon mass for the LORENE ID is in eq.(69)
          ! of Gourgoulhon et al., PRD 63 064029 (2001)

          CALL THIS% import_id( &
                   center + (rad_coord + dr)*SIN(colat)*COS(long), &
                   (rad_coord + dr)*SIN(colat)*SIN(long), &
                   (rad_coord + dr)*COS(colat), &
                   g_xx, baryon_density, gamma_euler )

!        CALL bns_obj% import_id( &
!                 center1 + rad_coord*SIN(lat)*COS(long), &
!                 rad_coord*SIN(lat)*SIN(long), &
!                 rad_coord*COS(lat), &
!                 g_xx, baryon_density, &
!                 gamma_euler )
!
!        ! Compute covariant spatial fluid velocity (metric is diagonal and
!        ! conformally flat)
!        !v_euler_x_l= g_xx*v_euler_x
!        !v_euler_y_l= g_xx*v_euler_y
!        !v_euler_z_l= g_xx*v_euler_z
!        !
!        !! Compute the corresponding Lorentz factor
!        !lorentz_factor= 1.0D0/SQRT( 1.0D0 - ( v_euler_x_l*v_euler_x &
!        !                                    + v_euler_y_l*v_euler_y &
!        !                                    + v_euler_z_l*v_euler_z ) )
!        !
!        !! Compute covariant fluid 4-velocity
!        !u_euler_t_l= lorentz_factor *( - lapse + v_euler_x_l*shift_x &
!        !                                       + v_euler_y_l*shift_y &
!        !                                       + v_euler_z_l*shift_z )
!        !u_euler_x_l= lorentz_factor*v_euler_x_l
!        !u_euler_y_l= lorentz_factor*v_euler_y_l
!        !u_euler_z_l= lorentz_factor*v_euler_z_l
!        !
!        !! Compute vector normal to spacelie hypersurface
!        !! (4-velocity of the Eulerian observer)
!        !n_t= 1.0D0/lapse
!        !n_x= - shift_x/lapse
!        !n_y= - shift_y/lapse
!        !n_z= - shift_z/lapse
!        !
!        !! Compute relative Lorentz factor between 4-velocity of the fluid
!        !! wrt the Eulerian observer and the 4-velocity of the Eulerian observer
!        !lorentz_factor_rel= - ( n_t*u_euler_t_l + n_x*u_euler_x_l &
!        !                      + n_y*u_euler_y_l + n_z*u_euler_z_l )

          ! Compute square root of the determinant of the spatial metric
          sq_g= g_xx*SQRT( g_xx )

          mass_element= (rad_coord**2.0D0)*SIN(colat)*dr*dth*dphi &
                        *sq_g*gamma_euler*baryon_density

          mass= mass + 2.0D0*mass_element

        ENDDO colatitude_loop

      ENDDO longitude_loop
!      !$OMP END PARALLEL SECTIONS

      mass_profile( 1, r )= rad_coord
      mass_profile( 2, r )= mass

      !PRINT *, rad_coord, mass_profile( 1, r ), mass, mass_profile( 2, r )
      !PRINT *

    ENDDO radius_loop
    !$OMP END PARALLEL DO

    DO r= 1, NINT(radius/dr), 1
      mass_profile( 3, r )= mass_profile( 3, r - 1 ) + mass_profile( 2, r )
    ENDDO

    mass= mass_profile( 3, NINT(radius/dr) )

    PRINT *, "radius covered by the integration=", MAXVAL( mass_profile( 1, : ), DIM= 1 )
    PRINT *, "mass=", mass
    PRINT *

    CALL indexx( NINT(radius/dr) + 1, mass_profile( 1, : ), mass_profile_idx )


  END PROCEDURE integrate_baryon_mass_density


  MODULE PROCEDURE destruct_binary

    !************************************************
    !                                               *
    ! Destruct the LORENE Bin_NS object and free    *
    ! the pointeri pointing to it                   *
    !                                               *
    ! FT                                            *
    !                                               *
    !************************************************

    IMPLICIT NONE

    !PRINT *, "** Executing the destruct_binary subroutine."

    IF ( C_ASSOCIATED( THIS% bns_ptr ) ) THEN

      CALL destruct_bin_ns( THIS% bns_ptr )
      THIS% bns_ptr = C_NULL_PTR

    ENDIF

    !PRINT *, "** Subroutine destruct_binary executed."
    !PRINT *

  END PROCEDURE destruct_binary

  MODULE PROCEDURE allocate_lorene_id_memory

    !************************************************
    !                                               *
    ! Allocate the memory to store the LORENE ID    *
    ! in the member arrays                          *
    !                                               *
    ! FT 17.09.2020                                 *
    !                                               *
    !************************************************

    IMPLICIT NONE

    PRINT *, "** Executing the allocate_lorene_id_memory subroutine..."

    IF(.NOT.ALLOCATED( THIS% lapse ))THEN
      ALLOCATE( THIS% lapse( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array lapse. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
    !  CALL test_status( ios, err_msg, &
    !              "...allocation error for array lapse" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% shift_x ))THEN
      ALLOCATE( THIS% shift_x( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array shift_x. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array shift_x" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% shift_y ))THEN
      ALLOCATE( THIS% shift_y( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array shift_y. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array shift_y" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% shift_z ))THEN
      ALLOCATE( THIS% shift_z( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array shift_z. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array shift_z" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_xx ))THEN
      ALLOCATE( THIS% g_xx( d ), STAT= ios, &
          ERRMSG = err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array g_xx. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array g_xx" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_xy ))THEN
      ALLOCATE( THIS% g_xy( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array g_xy. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array g_xy" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_xz ))THEN
      ALLOCATE( THIS% g_xz( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array g_xz. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array g_xz" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_yy ))THEN
      ALLOCATE( THIS% g_yy( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array g_yy. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array g_yy" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_yz ))THEN
      ALLOCATE( THIS% g_yz( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array g_yz. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array g_yz" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% g_zz ))THEN
      ALLOCATE( THIS% g_zz( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array g_zz. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array g_zz" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% k_xx ))THEN
      ALLOCATE( THIS% k_xx( d ), STAT= ios, &

          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array k_xx. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array k_xx" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% k_xy ))THEN
      ALLOCATE( THIS% k_xy( d ), STAT= ios, &

          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array k_xy. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array k_xy" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% k_xz ))THEN
      ALLOCATE( THIS% k_xz( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array k_xz. ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array k_xz" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% k_yy ))THEN
      ALLOCATE( THIS%  k_yy( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array k_yy ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array k_yy" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% k_yz ))THEN
      ALLOCATE( THIS% k_yz( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array k_yz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array k_yz" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% k_zz ))THEN
      ALLOCATE( THIS% k_zz( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array k_zz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array k_zz" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% baryon_density ))THEN
      ALLOCATE( THIS% baryon_density( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array baryon_density ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array baryon_density" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% energy_density ))THEN
      ALLOCATE( THIS% energy_density( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array energy_density ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array energy_density" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% specific_energy ))THEN
      ALLOCATE( THIS% specific_energy( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array specific_energy ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array specific_energy" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% v_euler_x ))THEN
      ALLOCATE( THIS% v_euler_x( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array v_euler_x ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array v_euler_x" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% v_euler_y ))THEN
      ALLOCATE( THIS% v_euler_y( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array v_euler_y ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array v_euler_y" )
    ENDIF
    IF(.NOT.ALLOCATED( THIS% v_euler_z ))THEN
      ALLOCATE( THIS% v_euler_z( d ), STAT= ios, &
          ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...allocation error for array v_euler_z ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !            "...allocation error for array v_euler_z" )
    ENDIF

    IF( SIZE( THIS% lapse ) /= d )THEN
      PRINT *, "** ERROR in memory allocation in allocate_lorene_id_memory"
    ENDIF

    PRINT *, "** Subroutine allocate_lorene_id_memory executed."
    PRINT *

  END PROCEDURE allocate_lorene_id_memory


  MODULE PROCEDURE deallocate_lorene_id_memory

    !************************************************
    !                                               *
    ! Deallocate the memory for the member arrays   *
    !                                               *
    ! FT 17.09.2020                                 *
    !                                               *
    !************************************************

    IMPLICIT NONE

    PRINT *, "** Executing the deallocate_lorene_id_memory subroutine..."

    IF(ALLOCATED( THIS% lapse ))THEN
      DEALLOCATE( THIS% lapse, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array lapse ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                  "...deallocation error for array lapse" )
    ENDIF
    IF(ALLOCATED( THIS% shift_x ))THEN
      DEALLOCATE( THIS% shift_x, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array shift_x ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array shift_x" )
    ENDIF
    IF(ALLOCATED( THIS% shift_y ))THEN
      DEALLOCATE( THIS% shift_y, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array shift_y ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array shift_y" )
    ENDIF
    IF(ALLOCATED( THIS% shift_z ))THEN
      DEALLOCATE( THIS% shift_z, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array shift_z ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array shift_z" )
    ENDIF
    IF(ALLOCATED( THIS% g_xx ))THEN
      DEALLOCATE( THIS% g_xx, STAT= ios, ERRMSG = err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_xx ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_xx" )
    ENDIF
    IF(ALLOCATED( THIS% g_xy ))THEN
      DEALLOCATE( THIS% g_xy, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_xy ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !               "...deallocation error for array g_xy" )
    ENDIF
    IF(ALLOCATED( THIS% g_xz ))THEN
      DEALLOCATE( THIS% g_xz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_xz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_xz" )
    ENDIF
    IF(ALLOCATED( THIS% g_yy ))THEN
      DEALLOCATE( THIS% g_yy, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_yy ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_yy" )
    ENDIF
    IF(ALLOCATED( THIS% g_yz ))THEN
      DEALLOCATE( THIS% g_yz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_yz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_yz" )
    ENDIF
    IF(ALLOCATED( THIS% g_zz ))THEN
      DEALLOCATE( THIS% g_zz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array g_zz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array g_zz" )
    ENDIF
    IF(ALLOCATED( THIS% k_xx ))THEN
      DEALLOCATE( THIS% k_xx, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array k_xx ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array k_xx" )
    ENDIF
    IF(ALLOCATED( THIS% k_xy ))THEN
      DEALLOCATE( THIS% k_xy, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array k_xy ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array k_xy" )
    ENDIF
    IF(ALLOCATED( THIS% k_xz ))THEN
      DEALLOCATE( THIS% k_xz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array k_xz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array k_xz" )
    ENDIF
    IF(ALLOCATED( THIS% k_yy ))THEN
      DEALLOCATE( THIS% k_yy, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array k_yy ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array k_yy" )
    ENDIF
    IF(ALLOCATED( THIS% k_yz ))THEN
      DEALLOCATE( THIS% k_yz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array k_yz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array k_yz" )
    ENDIF
    IF(ALLOCATED( THIS% k_zz ))THEN
      DEALLOCATE( THIS% k_zz, STAT= ios, ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array k_zz ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array k_zz" )
    ENDIF
    IF(ALLOCATED( THIS% baryon_density ))THEN
      DEALLOCATE( THIS% baryon_density, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array baryon_density ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...deallocation error for array baryon_density" )
    ENDIF
    IF(ALLOCATED( THIS% energy_density ))THEN
      DEALLOCATE( THIS% energy_density, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array energy_density ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...deallocation error for array energy_density" )
    ENDIF
    IF(ALLOCATED( THIS% specific_energy ))THEN
      DEALLOCATE( THIS% specific_energy, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array specific_energy ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...deallocation error for array specific_energy" )
    ENDIF
    IF(ALLOCATED( THIS% v_euler_x ))THEN
      DEALLOCATE( THIS% v_euler_x, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array v_euler_x ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !        "...deallocation error for array v_euler_x" )
    ENDIF
    IF(ALLOCATED( THIS% v_euler_y ))THEN
      DEALLOCATE( THIS% v_euler_y, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array v_euler_y ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v_euler_y" )
    ENDIF
    IF(ALLOCATED( THIS% v_euler_z ))THEN
      DEALLOCATE( THIS% v_euler_z, STAT= ios, &
              ERRMSG= err_msg )
      IF( ios > 0 )THEN
         PRINT *, "...deallocation error for array v_euler_z ", &
                  "The error message is", err_msg
         STOP
      ENDIF
      !CALL test_status( ios, err_msg, &
      !                "...deallocation error for array v_euler_z" )
    ENDIF

    PRINT *, "** Subroutine deallocate_lorene_id_memory executed."
    PRINT *

  END PROCEDURE deallocate_lorene_id_memory


  !-----------------!
  !--  FUNCTIONS  --!
  !-----------------!


  MODULE PROCEDURE import_mass_density

    !************************************************
    !                                               *
    ! Returns the LORENE mass density at the point  *
    ! given as argument, in units of                *
    ! MSun/(MSun_geo^3).                            *
    !                                               *
    ! FT                                            *
    !                                               *
    !************************************************

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_ASSOCIATED
    USE constants,                   ONLY: Msun_geo, lorene2hydrobase

    IMPLICIT NONE

    IF ( C_ASSOCIATED( THIS% bns_ptr ) )THEN

      ! The coordinates need to be converted from SPHINCS units (Msun_geo)
      ! to LORENE units (km). See MODULE constants for the definition of
      ! Msun_geo
      res= get_lorene_mass_density( THIS% bns_ptr, &
                                    x*Msun_geo, &
                                    y*Msun_geo, &
                                    z*Msun_geo )*lorene2hydrobase

    ENDIF

  END PROCEDURE import_mass_density


  MODULE PROCEDURE import_spatial_metric

    !************************************************
    !                                               *
    ! Returns the LORENE conformal factor to the    *
    ! 4th power, equal to the diagonal components   *
    ! of the conformally flat spatial ADM metric.   *
    !                                               *
    ! FT 15.04.2021                                 *
    !                                               *
    !************************************************

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_ASSOCIATED
    USE constants,                   ONLY: Msun_geo

    IMPLICIT NONE

    IF ( C_ASSOCIATED( THIS% bns_ptr ) )THEN

      ! The coordinates need to be converted from SPHINCS units (Msun_geo)
      ! to LORENE units (km). See MODULE constants for the definition of
      ! Msun_geo
      res= get_lorene_spatial_metric( THIS% bns_ptr, &
                                      x*Msun_geo, &
                                      y*Msun_geo, &
                                      z*Msun_geo )

    ENDIF

  END PROCEDURE import_spatial_metric


  MODULE PROCEDURE is_hydro_negative

    !************************************************
    !                                               *
    ! Return 1 if the energy density is nonpositive,*
    ! or if the specific energy is nonpositive,     *
    ! or if the pressure is nonpositive             *
    ! at the specified point                        *
    !                                               *
    ! FT 12.03.2021                                 *
    !                                               *
    !************************************************

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_ASSOCIATED
    USE constants,                   ONLY: Msun_geo

    IMPLICIT NONE

    IF ( C_ASSOCIATED( THIS% bns_ptr ) )THEN

      ! The coordinates need to be converted from SPHINCS units (Msun_geo)
      ! to LORENE units (km). See MODULE constants for the definition of
      ! Msun_geo
      res= negative_hydro( THIS% bns_ptr, &
                                    x*Msun_geo, &
                                    y*Msun_geo, &
                                    z*Msun_geo )

    ENDIF

  END PROCEDURE is_hydro_negative


  MODULE PROCEDURE get_field_array

    !************************************************
    !                                               *
    ! Returns one of the member arrays, selected    *
    ! with the string input.                        *
    !                                               *
    ! FT                                            *
    !                                               *
    !************************************************

    IMPLICIT NONE

    select_field: SELECT CASE( field )

    CASE( "lapse" )

      field_array= THIS% lapse

    CASE( "shift_x" )

      field_array= THIS% shift_x

    CASE( "shift_y" )

      field_array= THIS% shift_y

    CASE( "shift_z" )

      field_array= THIS% shift_z

    CASE( "g_xx" )

      field_array= THIS% g_xx

    CASE( "g_xy" )

      field_array= THIS% g_xy

    CASE( "g_xz" )

      field_array= THIS% g_xz

    CASE( "g_yy" )

      field_array= THIS% g_yy

    CASE( "g_yz" )

      field_array= THIS% g_yz

    CASE( "g_zz" )

      field_array= THIS% g_zz

    CASE( "k_xx" )

      field_array= THIS% k_xx

    CASE( "k_xy" )

      field_array= THIS% k_xy

    CASE( "k_xz" )

      field_array= THIS% k_xz

    CASE( "k_yy" )

      field_array= THIS% k_yy

    CASE( "k_yz" )

      field_array= THIS% k_yz

    CASE( "k_zz" )

      field_array= THIS% k_zz

    CASE( "baryon_density" )

      field_array= THIS% baryon_density

    CASE( "energy_density" )

      field_array= THIS% energy_density

    CASE( "specific_energy" )

      field_array= THIS% specific_energy

    CASE( "v_euler_x" )

      field_array= THIS% v_euler_x

    CASE( "v_euler_y" )

      field_array= THIS% v_euler_y

    CASE( "v_euler_z" )

      field_array= THIS% v_euler_z

    CASE DEFAULT

      PRINT *, "** There is no field named ", field, "in TYPE bns."
      STOP

    END SELECT select_field

  END PROCEDURE get_field_array


  MODULE PROCEDURE get_field_value

    !*************************************************
    !                                                *
    ! Returns the value of one of the member arrays, *
    ! selected with the string input, at the point   *
    ! given as argument.                             *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    select_field: SELECT CASE( field )

    CASE( "lapse" )

      field_value= THIS% lapse( n )

    CASE( "shift_x" )

      field_value= THIS% shift_x( n )

    CASE( "shift_y" )

      field_value= THIS% shift_y( n )

    CASE( "shift_z" )

      field_value= THIS% shift_z( n )

    CASE( "g_xx" )

      field_value= THIS% g_xx( n )

    CASE( "g_xy" )

      field_value= THIS% g_xy( n )

    CASE( "g_xz" )

      field_value= THIS% g_xz( n )

    CASE( "g_yy" )

      field_value= THIS% g_yy( n )

    CASE( "g_yz" )

      field_value= THIS% g_yz( n )

    CASE( "g_zz" )

      field_value= THIS% g_zz( n )

    CASE( "k_xx" )

      field_value= THIS% k_xx( n )

    CASE( "k_xy" )

      field_value= THIS% k_xy( n )

    CASE( "k_xz" )

      field_value= THIS% k_xz( n )

    CASE( "k_yy" )

      field_value= THIS% k_yy( n )

    CASE( "k_yz" )

      field_value= THIS% k_yz( n )

    CASE( "k_zz" )

      field_value= THIS% k_zz( n )

    CASE( "baryon_density" )

      field_value= THIS% baryon_density( n )

    CASE( "energy_density" )

      field_value= THIS% energy_density( n )

    CASE( "specific_energy" )

      field_value= THIS% specific_energy( n )

    CASE( "v_euler_x" )

      field_value= THIS% v_euler_x( n )

    CASE( "v_euler_y" )

      field_value= THIS% v_euler_y( n )

    CASE( "v_euler_z" )

      field_value= THIS% v_euler_z( n )

    CASE DEFAULT

      PRINT *, "** There is no field named ", field, "in TYPE bns."
      STOP

    END SELECT select_field

  END PROCEDURE get_field_value


  MODULE PROCEDURE get_bns_identifier

    !*************************************************
    !                                                *
    ! Returns the value of bns_identifier, the       *
    ! integer identifier of the bns object           *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_bns_identifier= THIS% bns_identifier

  END PROCEDURE get_bns_identifier


  !MODULE PROCEDURE get_bns_ptr
  !
  !  !*************************************************
  !  !                                                *
  !  ! Returns the value of bns_ptr, the C pointer    *
  !  ! to the LORENE's Bin_NS object                  *
  !  ! N.B. This variable is global. The pointer      *
  !  !      to the second LORENE Bin_NS object will   *
  !  !      overwrite the first one, and so on.       *
  !  !      This variable stores the pointer to       *
  !  !      the last defined LORENE Bin_NS object.    *
  !  !      That's why it is not freed in the         *
  !  !      destructor of a bns object. Presently, it *
  !  !      has to be freed by the user at the end of *
  !  !      the PROGRAM. See the last part of the     *
  !  !      PROGRAM in setup_lorene_id.f90, for       *
  !  !      example.                                  *
  !  !                                                *
  !  ! FT                                             *
  !  !                                                *
  !  !*************************************************
  !
  !  IMPLICIT NONE
  !
  !  get_bns_ptr= THIS% bns_ptr
  !
  !END PROCEDURE get_bns_ptr


  MODULE PROCEDURE get_gamma_1

    !*************************************************
    !                                                *
    ! Returns the value of gamma_1, the              *
    ! polytropic index for NS 1 with polytropic EOS, *
    ! not piecewise polytropic EOS                   *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_gamma_1= THIS% gamma_1

  END PROCEDURE get_gamma_1


  MODULE PROCEDURE get_gamma_2

    !*************************************************
    !                                                *
    ! Returns the value of gamma_2, the              *
    ! polytropic index for NS 2 with polytropic EOS, *
    ! not piecewise polytropic EOS                   *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_gamma_2= THIS% gamma_2

  END PROCEDURE get_gamma_2


  MODULE PROCEDURE get_kappa_1

    !*************************************************
    !                                                *
    ! Returns the value of kappa_1, the              *
    ! polytropic constant for NS 1 with polytropic   *
    ! EOS, not piecewise polytropic EOS              *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_kappa_1= THIS% kappa_1

  END PROCEDURE get_kappa_1


  MODULE PROCEDURE get_kappa_2

    !*************************************************
    !                                                *
    ! Returns the value of kappa_2, the              *
    ! polytropic constant for NS 2 with polytropic   *
    ! EOS, not piecewise polytropic EOS              *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_kappa_2= THIS% kappa_2

  END PROCEDURE get_kappa_2


  MODULE PROCEDURE get_angular_vel

    !*************************************************
    !                                                *
    ! Returns the angular velocity of the system     *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_angular_vel= THIS% angular_vel

  END PROCEDURE get_angular_vel


  MODULE PROCEDURE get_distance

    !*************************************************
    !                                                *
    ! Returns the distance between the NSs           *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_distance= THIS% distance

  END PROCEDURE get_distance


  MODULE PROCEDURE get_distance_com

    !*************************************************
    !                                                *
    ! Returns the distance between the centers of    *
    ! mass of the NSs                                *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_distance_com= THIS% distance_com

  END PROCEDURE get_distance_com


  MODULE PROCEDURE get_mass1

    !*************************************************
    !                                                *
    ! Returns the mass of NS 1 [Msun]                *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_mass1= THIS% mass1

  END PROCEDURE get_mass1


  MODULE PROCEDURE get_mass2

    !*************************************************
    !                                                *
    ! Returns the mass of NS 2 [Msun]                *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_mass2= THIS% mass2

  END PROCEDURE get_mass2


  MODULE PROCEDURE get_adm_mass

    !*************************************************
    !                                                *
    ! Returns the ADM mass of the system             *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_adm_mass= THIS% adm_mass

  END PROCEDURE get_adm_mass


  MODULE PROCEDURE get_angular_momentum

    !*************************************************
    !                                                *
    ! Returns the angular momentum of the system     *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_angular_momentum= THIS% angular_momentum

  END PROCEDURE get_angular_momentum


  MODULE PROCEDURE get_radius1_x_comp

    !*************************************************
    !                                                *
    ! Returns the radius of NS 1 along the x axis    *
    ! on the side of the companion                   *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_radius1_x_comp= THIS% radius1_x_comp

  END PROCEDURE get_radius1_x_comp


  MODULE PROCEDURE get_radius1_y

    !*************************************************
    !                                                *
    ! Returns the radius of NS 1 along the y axis    *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_radius1_y= THIS% radius1_y

  END PROCEDURE get_radius1_y


  MODULE PROCEDURE get_radius1_z

    !*************************************************
    !                                                *
    ! Returns the radius of NS 1 along the z axis    *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_radius1_z= THIS% radius1_z

  END PROCEDURE get_radius1_z


  MODULE PROCEDURE get_radius1_x_opp

    !*************************************************
    !                                                *
    ! Returns the radius of NS 1 along the x axis    *
    ! on the side opposite to the companion          *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_radius1_x_opp= THIS% radius1_x_opp

  END PROCEDURE get_radius1_x_opp


  MODULE PROCEDURE get_center1_x

    !*************************************************
    !                                                *
    ! Returns the stellar center of NS 1, i.e., the  *
    ! origin of the LORENE chart centered on NS 1    *
    !                                                *
    ! FT 09.02.2021                                  *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_center1_x= THIS% center1_x

  END PROCEDURE get_center1_x


  MODULE PROCEDURE get_barycenter1_x

    !*************************************************
    !                                                *
    ! Returns the barycenter of NS 1                 *
    !                                                *
    ! FT 09.02.2021                                  *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_barycenter1_x= THIS% barycenter1_x

  END PROCEDURE get_barycenter1_x


  MODULE PROCEDURE get_radius2_x_comp

    !*************************************************
    !                                                *
    ! Returns the radius of NS 2 along the x axis    *
    ! on the side of the companion                   *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_radius2_x_comp= THIS% radius2_x_comp

  END PROCEDURE get_radius2_x_comp


  MODULE PROCEDURE get_radius2_y

    !*************************************************
    !                                                *
    ! Returns the radius of NS 2 along the y axis    *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_radius2_y= THIS% radius2_y

  END PROCEDURE get_radius2_y


  MODULE PROCEDURE get_radius2_z

    !*************************************************
    !                                                *
    ! Returns the radius of NS 2 along the z axis    *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_radius2_z= THIS% radius2_z

  END PROCEDURE get_radius2_z


  MODULE PROCEDURE get_radius2_x_opp

    !*************************************************
    !                                                *
    ! Returns the radius of NS 2 along the x axis    *
    ! on the side opposite to the companion          *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_radius2_x_opp= THIS% radius2_x_opp

  END PROCEDURE get_radius2_x_opp


  MODULE PROCEDURE get_center2_x

    !*************************************************
    !                                                *
    ! Returns the stellar center of NS 2, i.e., the  *
    ! origin of the LORENE chart centered on NS 2    *
    !                                                *
    ! FT 09.02.2021                                  *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_center2_x= THIS% center2_x

  END PROCEDURE get_center2_x


  MODULE PROCEDURE get_barycenter2_x

    !*************************************************
    !                                                *
    ! Returns the barycenter of NS 2                 *
    !                                                *
    ! FT 09.02.2021                                  *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_barycenter2_x= THIS% barycenter2_x

  END PROCEDURE get_barycenter2_x


  MODULE PROCEDURE get_ent_center1

    !*************************************************
    !                                                *
    ! Returns the central enthalpy of NS 1           *
    !                                                *
    ! FT 12.02.2021                                  *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_ent_center1= THIS% ent_center1

  END PROCEDURE get_ent_center1


  MODULE PROCEDURE get_nbar_center1

    !*************************************************
    !                                                *
    ! Returns the central baryon number density      *
    ! of NS 1                                        *
    !                                                *
    ! FT 12.02.2021                                  *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_nbar_center1= THIS% nbar_center1

  END PROCEDURE get_nbar_center1


  MODULE PROCEDURE get_rho_center1

    !*************************************************
    !                                                *
    ! Returns the central baryon mass density        *
    ! of NS 1                                        *
    !                                                *
    ! FT 12.02.2021                                  *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_rho_center1= THIS% rho_center1

  END PROCEDURE get_rho_center1


  MODULE PROCEDURE get_energy_density_center1

    !*************************************************
    !                                                *
    ! Returns the central energy density of NS 1     *
    !                                                *
    ! FT 12.02.2021                                  *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_energy_density_center1= THIS% energy_density_center1

  END PROCEDURE get_energy_density_center1


  MODULE PROCEDURE get_specific_energy_center1

    !*************************************************
    !                                                *
    ! Returns the central specific energy of NS 1    *
    !                                                *
    ! FT 12.02.2021                                  *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_specific_energy_center1= THIS% specific_energy_center1

  END PROCEDURE get_specific_energy_center1


  MODULE PROCEDURE get_pressure_center1

    !*************************************************
    !                                                *
    ! Returns the central pressure of NS 1           *
    !                                                *
    ! FT 12.02.2021                                  *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_pressure_center1= THIS% pressure_center1

  END PROCEDURE get_pressure_center1


  MODULE PROCEDURE get_ent_center2

    !*************************************************
    !                                                *
    ! Returns the central enthalpy of NS 2           *
    !                                                *
    ! FT 12.02.2021                                  *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_ent_center2= THIS% ent_center2

  END PROCEDURE get_ent_center2


  MODULE PROCEDURE get_nbar_center2

    !*************************************************
    !                                                *
    ! Returns the central baryon number density      *
    ! of NS 2                                        *
    !                                                *
    ! FT 12.02.2021                                  *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_nbar_center2= THIS% nbar_center2

  END PROCEDURE get_nbar_center2


  MODULE PROCEDURE get_rho_center2

    !*************************************************
    !                                                *
    ! Returns the central baryon mass density        *
    ! of NS 2                                        *
    !                                                *
    ! FT 12.02.2021                                  *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_rho_center2= THIS% rho_center2

  END PROCEDURE get_rho_center2


  MODULE PROCEDURE get_energy_density_center2

    !*************************************************
    !                                                *
    ! Returns the central energy density of NS 2     *
    !                                                *
    ! FT 12.02.2021                                  *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_energy_density_center2= THIS% energy_density_center2

  END PROCEDURE get_energy_density_center2


  MODULE PROCEDURE get_specific_energy_center2

    !*************************************************
    !                                                *
    ! Returns the central specific energy of NS 2    *
    !                                                *
    ! FT 12.02.2021                                  *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_specific_energy_center2= THIS% specific_energy_center2

  END PROCEDURE get_specific_energy_center2


  MODULE PROCEDURE get_pressure_center2

    !*************************************************
    !                                                *
    ! Returns the central pressure of NS 2           *
    !                                                *
    ! FT 12.02.2021                                  *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_pressure_center2= THIS% pressure_center2

  END PROCEDURE get_pressure_center2


  MODULE PROCEDURE get_eos1

    !*************************************************
    !                                                *
    ! Returns the name of the EOS for NS 1           *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_eos1= THIS% eos1

  END PROCEDURE get_eos1


  MODULE PROCEDURE get_eos2

    !*************************************************
    !                                                *
    ! Returns the name of the EOS for NS 2           *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_eos2= THIS% eos2

  END PROCEDURE get_eos2


  MODULE PROCEDURE get_npeos_1

    !*************************************************
    !                                                *
    ! Returns the identifier of the EOS for NS 1     *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_npeos_1= THIS% npeos_1

  END PROCEDURE get_npeos_1


  MODULE PROCEDURE get_npeos_2

    !*************************************************
    !                                                *
    ! Returns the identifier of the EOS for NS 2     *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_npeos_2= THIS% npeos_2

  END PROCEDURE get_npeos_2


  MODULE PROCEDURE get_gamma0_1

    !*************************************************
    !                                                *
    ! Returns the value of gamma0_1, the crust's     *
    ! polytropic index for NS 1 with piecewise       *
    ! polytropic EOS                                 *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_gamma0_1= THIS% gamma0_1

  END PROCEDURE get_gamma0_1


  MODULE PROCEDURE get_gamma0_2

    !*************************************************
    !                                                *
    ! Returns the value of gamma0_2, the crust's     *
    ! polytropic index for NS 2 with piecewise       *
    ! polytropic EOS                                 *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_gamma0_2= THIS% gamma0_2

  END PROCEDURE get_gamma0_2


  MODULE PROCEDURE get_gamma1_1

    !*************************************************
    !                                                *
    ! Returns the value of gamma1_1, the first       *
    ! polytropic index for NS 1 with piecewise       *
    ! polytropic EOS                                 *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_gamma1_1= THIS% gamma1_1

  END PROCEDURE get_gamma1_1


  MODULE PROCEDURE get_gamma1_2

    !*************************************************
    !                                                *
    ! Returns the value of gamma1_2, the first       *
    ! polytropic index for NS 2 with piecewise       *
    ! polytropic EOS                                 *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_gamma1_2= THIS% gamma1_2

  END PROCEDURE get_gamma1_2


  MODULE PROCEDURE get_gamma2_1

    !*************************************************
    !                                                *
    ! Returns the value of gamma2_1, the second      *
    ! polytropic index for NS 2 with piecewise       *
    ! polytropic EOS                                 *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_gamma2_1= THIS% gamma2_1

  END PROCEDURE get_gamma2_1


  MODULE PROCEDURE get_gamma2_2

    !*************************************************
    !                                                *
    ! Returns the value of gamma2_2, the second      *
    ! polytropic index for NS 2 with piecewise       *
    ! polytropic EOS                                 *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_gamma2_2= THIS% gamma2_2

  END PROCEDURE get_gamma2_2


  MODULE PROCEDURE get_gamma3_1

    !*************************************************
    !                                                *
    ! Returns the value of gamma3_1, the third       *
    ! polytropic index for NS 1 with piecewise       *
    ! polytropic EOS (innermost index)               *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_gamma3_1= THIS% gamma3_1

  END PROCEDURE get_gamma3_1


  MODULE PROCEDURE get_gamma3_2

    !*************************************************
    !                                                *
    ! Returns the value of gamma3_2, the third       *
    ! polytropic index for NS 2 with piecewise       *
    ! polytropic EOS (innermost index)               *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_gamma3_2= THIS% gamma3_2

  END PROCEDURE get_gamma3_2


  MODULE PROCEDURE get_kappa0_1

    !*************************************************
    !                                                *
    ! Returns the value of kappa0_1, the crust's     *
    ! polytropic constant for NS 1 with piecewise    *
    ! polytropic EOS                                 *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_kappa0_1= THIS% kappa0_1

  END PROCEDURE get_kappa0_1


  MODULE PROCEDURE get_kappa1_1

    !*************************************************
    !                                                *
    ! Returns the value of kappa1_1, the first       *
    ! polytropic constant for NS 1 with piecewise    *
    ! polytropic EOS                                 *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_kappa1_1= THIS% kappa1_1

  END PROCEDURE get_kappa1_1


  MODULE PROCEDURE get_kappa2_1

    !*************************************************
    !                                                *
    ! Returns the value of kappa2_1, the second      *
    ! polytropic constant for NS 1 with piecewise    *
    ! polytropic EOS                                 *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_kappa2_1= THIS% kappa2_1

  END PROCEDURE get_kappa2_1


  MODULE PROCEDURE get_kappa3_1

    !*************************************************
    !                                                *
    ! Returns the value of kappa3_1, the third       *
    ! polytropic constant for NS 1 with piecewise    *
    ! polytropic EOS                                 *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_kappa3_1= THIS% kappa3_1

  END PROCEDURE get_kappa3_1


  MODULE PROCEDURE get_kappa0_2

    !*************************************************
    !                                                *
    ! Returns the value of kappa0_2, the crust's     *
    ! polytropic constant for NS 2 with piecewise    *
    ! polytropic EOS                                 *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_kappa0_2= THIS% kappa0_2

  END PROCEDURE get_kappa0_2


  MODULE PROCEDURE get_kappa1_2

    !*************************************************
    !                                                *
    ! Returns the value of kappa1_2, the first       *
    ! polytropic constant for NS 2 with piecewise    *
    ! polytropic EOS                                 *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_kappa1_2= THIS% kappa1_2

  END PROCEDURE get_kappa1_2


  MODULE PROCEDURE get_kappa2_2

    !*************************************************
    !                                                *
    ! Returns the value of kappa2_2, the second      *
    ! polytropic constant for NS 2 with piecewise    *
    ! polytropic EOS                                 *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_kappa2_2= THIS% kappa2_2

  END PROCEDURE get_kappa2_2


  MODULE PROCEDURE get_kappa3_2

    !*************************************************
    !                                                *
    ! Returns the value of kappa3_2, the third       *
    ! polytropic constant for NS 2 with piecewise    *
    ! polytropic EOS                                 *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_kappa3_2= THIS% kappa3_2

  END PROCEDURE get_kappa3_2


  MODULE PROCEDURE get_logp1_1

    !*************************************************
    !                                                *
    ! Returns the value of logp1_1, the base 10      *
    ! logarithm of the pressure where the gamma1_1   *
    ! polytrope starts, for NS 1 with piecewise      *
    ! polytropic EOS                                 *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_logp1_1= THIS% logp1_1

  END PROCEDURE get_logp1_1


  MODULE PROCEDURE get_logp1_2

    !*************************************************
    !                                                *
    ! Returns the value of logp1_2, the base 10      *
    ! logarithm of the pressure where the gamma1_2   *
    ! polytrope starts, for NS 2 with piecewise      *
    ! polytropic EOS                                 *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_logp1_2= THIS% logp1_2

  END PROCEDURE get_logp1_2


  MODULE PROCEDURE get_logRho0_1

    !*************************************************
    !                                                *
    ! Returns the value of logRho0_1, the base 10    *
    ! logarithm of the mass density where the        *
    ! gamma1_1 polytrope starts, for NS 1 with       *
    ! piecewise polytropic EOS                       *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_logRho0_1= THIS% logRho0_1

  END PROCEDURE get_logRho0_1


  MODULE PROCEDURE get_logRho0_2

    !*************************************************
    !                                                *
    ! Returns the value of logRho0_2, the base 10    *
    ! logarithm of the mass density where the        *
    ! gamma1_2 polytrope starts, for NS 2 with       *
    ! piecewise polytropic EOS                       *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_logRho0_2= THIS% logRho0_2

  END PROCEDURE get_logRho0_2


  MODULE PROCEDURE get_logRho1_1

    !*************************************************
    !                                                *
    ! Returns the value of logRho1_1, the base 10    *
    ! logarithm of the mass density where the        *
    ! gamma2_1 polytrope starts, for NS 1 with       *
    ! piecewise polytropic EOS                       *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_logRho1_1= THIS% logRho1_1

  END PROCEDURE get_logRho1_1


  MODULE PROCEDURE get_logRho1_2

    !*************************************************
    !                                                *
    ! Returns the value of logRho1_2, the base 10    *
    ! logarithm of the mass density where the        *
    ! gamma2_2 polytrope starts, for NS 2 with       *
    ! piecewise polytropic EOS                       *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_logRho1_2= THIS% logRho1_2

  END PROCEDURE get_logRho1_2


  MODULE PROCEDURE get_logRho2_1

    !*************************************************
    !                                                *
    ! Returns the value of logRho2_1, the base 10    *
    ! logarithm of the mass density where the        *
    ! gamma3_1 polytrope starts, for NS 1 with       *
    ! piecewise polytropic EOS                       *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_logRho2_1= THIS% logRho2_1

  END PROCEDURE get_logRho2_1


  MODULE PROCEDURE get_logRho2_2

    !*************************************************
    !                                                *
    ! Returns the value of logRho2_2, the base 10    *
    ! logarithm of the mass density where the        *
    ! gamma3_2 polytrope starts, for NS 2 with       *
    ! piecewise polytropic EOS                       *
    !                                                *
    ! FT                                             *
    !                                                *
    !*************************************************

    IMPLICIT NONE

    get_logRho2_2= THIS% logRho2_2

  END PROCEDURE get_logRho2_2


END SUBMODULE bns_methods
