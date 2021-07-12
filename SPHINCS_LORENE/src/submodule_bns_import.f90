! File:         submodule_bns_methods.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (bns_id) bns_import

  !*****************************************************
  !                                                    *
  ! Implementation of the methods of TYPE bns that     *
  ! import BNS data using LORENE                       *
  !                                                    *
  ! FT 23.10.2020                                      *
  !                                                    *
  ! Renamed from bns_methods to bns_import upon        *
  ! improving modularity                               *
  !                                                    *
  ! FT 12.07.2021                                      *
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


END SUBMODULE bns_import
