! File:         submodule_bns_methods.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

SUBMODULE (diffstar_lorene) diffstar_lorene_import

  !****************************************************
  !
  !# Implementation of the methods of TYPE bns that
  !  import BNS data using |lorene|
  !
  !  FT 25.10.2021
  !
  !****************************************************


  USE, INTRINSIC:: ISO_C_BINDING, ONLY: C_ASSOCIATED


  IMPLICIT NONE


  CONTAINS


  !-------------------!
  !--  SUBROUTINES  --!
  !-------------------!


  MODULE PROCEDURE import_id_int

    !**************************************************
    !
    !# Stores the ID in the [[diffstarlorene]] member arrays
    !
    !  FT 25.10.2021
    !
    !**************************************************

    USE constants, ONLY: Msun_geo

    IMPLICIT NONE

    IF ( C_ASSOCIATED( THIS% diffstar_ptr ) ) THEN

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

        !$OMP PARALLEL DO DEFAULT( NONE ) &
        !$OMP             SHARED( n, THIS, x, y, z ) &
        !$OMP             PRIVATE( itr )
        import_id_loop: DO itr= 1, n, 1

          ! The coordinates need to be converted from |sphincs| units (Msun_geo)
          ! to |lorene| units (\(\mathrm{km}\)). See MODULE constants for the definition of
          ! Msun_geo
          CALL get_diffstar_full( THIS% diffstar_ptr, &
                                  x( itr )*Msun_geo, &
                                  y( itr )*Msun_geo, &
                                  z( itr )*Msun_geo, &
                                  THIS% lapse( itr ), &
                                  THIS% shift_x( itr ), &
                                  THIS% shift_y( itr ), &
                                  THIS% shift_z( itr ), &
                                  THIS% g_xx( itr ), &
                                  THIS% g_yy( itr ), &
                                  THIS% g_zz( itr ), &
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

        ENDDO import_id_loop
        !$OMP END PARALLEL DO

        DO itr= 1, n, 1

          !
          !-- The following follows from the assumption of conformal
          !-- flatness in |lorene|
          !
          THIS% g_xy( itr )= 0.0D0
          THIS% g_xz( itr )= 0.0D0
          THIS% g_yz( itr )= 0.0D0

          !
          !- Set/unset the geodesic gauge
          !
          IF( THIS% get_one_lapse() )THEN
            THIS% lapse( itr )= 1.0D0
          ENDIF
          IF( THIS% get_zero_shift() )THEN
            THIS% shift_x( itr )= 0.0D0
            THIS% shift_y( itr )= 0.0D0
            THIS% shift_z( itr )= 0.0D0
          ENDIF

          !
          !-- Convert the extrinsic curvature from |lorene| units to
          !-- |sphincs| units
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

        ENDDO
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


  MODULE PROCEDURE import_id_full

    !**************************************************
    !
    !# Stores the ID in non-[[diffstarlorene]]-member arrays
    !  with the same shape as the [[diffstarlorene]] member arrays
    !
    !  FT 25.10.2021
    !
    !**************************************************

    USE constants, ONLY: Msun_geo

    IMPLICIT NONE

    IF ( C_ASSOCIATED( THIS% diffstar_ptr ) ) THEN

      IF( SIZE( x ) /= SIZE( y ) .OR. SIZE( x ) /= SIZE( z ) &
            .OR. SIZE( y ) /= SIZE( z ) )THEN
        PRINT *, "** ERROR: The sizes of the arrays of positions" &
                 // "passed to import_lorene_id are not the same."
        PRINT *
        STOP
      ENDIF

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( n, THIS, x, y, z, lapse, &
      !$OMP                     shift_x, shift_y, shift_z, &
      !$OMP                     g_xx, g_yy, g_zz, &
      !$OMP                     k_xx, k_xy, k_xz, k_yy, k_yz, k_zz, &
      !$OMP                     baryon_density, energy_density, &
      !$OMP                     specific_energy, &
      !$OMP                     u_euler_x, u_euler_y, u_euler_z ) &
      !$OMP             PRIVATE( itr )
      import_id_loop: DO itr= 1, n, 1

        ! The coordinates need to be converted from |sphincs| units (Msun_geo)
        ! to |lorene| units (\(\mathrm{km}\)). See MODULE constants for the definition of
        ! Msun_geo
        CALL get_diffstar_full( THIS% diffstar_ptr, &
                                x( itr )*Msun_geo, &
                                y( itr )*Msun_geo, &
                                z( itr )*Msun_geo, &
                                lapse( itr ), &
                                shift_x( itr ), &
                                shift_y( itr ), &
                                shift_z( itr ), &
                                g_xx( itr ), &
                                g_yy( itr ), &
                                g_zz( itr ), &
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

      ENDDO import_id_loop
      !$OMP END PARALLEL DO

      DO itr= 1, n, 1

        !
        !-- The following follows from the maximal-slicing-quasi-isotropic
        !-- coordinates used in |lorene|
        !
        g_xy( itr )= 0.0D0
        g_xz( itr )= 0.0D0
        g_yz( itr )= 0.0D0

        !
        !- Set/unset the geodesic gauge
        !
        IF( THIS% get_one_lapse() )THEN
          lapse( itr )= 1.0D0
        ENDIF
        IF( THIS% get_zero_shift() )THEN
          shift_x( itr )= 0.0D0
          shift_y( itr )= 0.0D0
          shift_z( itr )= 0.0D0
        ENDIF

        !
        !-- Convert the extrinsic curvature from |lorene| units to
        !-- |sphincs| units
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

      ENDDO
      IF( show_progress ) WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

      PRINT *, "** Subroutine import_lorene_id executed."
      PRINT *

    ENDIF

  END PROCEDURE import_id_full


  MODULE PROCEDURE import_id_spacetime

    !*******************************************************
    !
    !# Stores the spacetime ID in multi-dimensional arrays
    !  needed to compute the BSSN variables and constraints
    !
    !  FT 22.11.2020
    !
    !*******************************************************

    USE constants, ONLY: Msun_geo

    USE tensor,    ONLY: itt, itx, ity, itz, ixx, ixy, &
                         ixz, iyy, iyz, izz, jxx, jxy, jxz, &
                         jyy, jyz, jzz, jx, jy, jz, n_sym4x4

    IMPLICIT NONE

    INTEGER:: i, j, k

    DOUBLE PRECISION:: detg
    DOUBLE PRECISION:: detg4
    DOUBLE PRECISION, DIMENSION( :, :, :, : ), ALLOCATABLE:: g4

    ! g4 is allocatable to allocate it on the heap
    ! Allocating it on the stack might exceed stack memory,
    ! causing a segmentation fault
    ALLOCATE( g4( nx, ny, nz, n_sym4x4 ) )

    IF ( C_ASSOCIATED( THIS% diffstar_ptr ) ) THEN

      IF( .FALSE. &!SHAPE( pos(:,:,:,1) ) /= SHAPE( lapse ) .OR. &
          !SHAPE( pos(:,:,:,1) ) /= SHAPE( shift(:,:,:,jx) ) & ! .OR. &
        ! SHAPE( pos(:,:,:,1) ) /= SHAPE( g(:,:,:,1) ) .OR. &
        ! SHAPE( pos(:,:,:,1) ) /= SHAPE( k(:,:,:,1) ) &
        )THEN
        PRINT *, "** ERROR: Mismatch in array dimensions" &
                 // "in import_id_spacetime."
        PRINT *
        STOP
      ENDIF

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( nx, ny, nz, THIS, pos, &
      !$OMP                     lapse, shift, g, ek ) &
      !$OMP             PRIVATE( i, j, k )
      coords_z: DO k= 1, nz, 1
        coords_y: DO j= 1, ny, 1
          coords_x: DO i= 1, nx, 1

            ! The coordinates need to be converted from |sphincs| units (Msun_geo)
            ! to |lorene| units (\(\mathrm{km}\)). See MODULE constants for the definition of
            ! Msun_geo
            CALL get_diffstar_spacetime( THIS% diffstar_ptr, &
                                         pos( i, j, k, jx )*Msun_geo, &
                                         pos( i, j, k, jy )*Msun_geo, &
                                         pos( i, j, k, jz )*Msun_geo, &
                                         lapse( i, j, k ), &
                                         shift( i, j, k, jx ), &
                                         shift( i, j, k, jy ), &
                                         shift( i, j, k, jz ), &
                                         g( i, j, k, jxx ), &
                                         g( i, j, k, jyy ), &
                                         g( i, j, k, jzz ), &
                                         ek( i, j, k, jxx ), &
                                         ek( i, j, k, jxy ), &
                                         ek( i, j, k, jxz ), &
                                         ek( i, j, k, jyy ), &
                                         ek( i, j, k, jyz ), &
                                         ek( i, j, k, jzz ) )

          ENDDO coords_x
        ENDDO coords_y
      ENDDO coords_z
      !$OMP END PARALLEL DO

      DO k= 1, nz, 1
        DO j= 1, ny, 1
          DO i= 1, nx, 1

            !
            !-- The following follows from the assumption of
            !-- conformal flatness in |lorene|
            !
            g( i, j, k, jxy )= 0.0D0
            g( i, j, k, jxz )= 0.0D0
            g( i, j, k, jyz )= 0.0D0

            !
            !- Set/unset the geodesic gauge
            !
            IF( THIS% get_one_lapse() )THEN
              lapse( i, j, k )= 1.0D0
            ENDIF
            IF( THIS% get_zero_shift() )THEN
              shift( i, j, k, jx )= 0.0D0
              shift( i, j, k, jy )= 0.0D0
              shift( i, j, k, jz )= 0.0D0
            ENDIF

            !
            !-- Convert the extrinsic curvature from |lorene| units to
            !-- |sphincs| units
            !
            ek( i, j, k, jxx )= ek( i, j, k, jxx )*Msun_geo
            ek( i, j, k, jxy )= ek( i, j, k, jxy )*Msun_geo
            ek( i, j, k, jxz )= ek( i, j, k, jxz )*Msun_geo
            ek( i, j, k, jyy )= ek( i, j, k, jyy )*Msun_geo
            ek( i, j, k, jyz )= ek( i, j, k, jyz )*Msun_geo
            ek( i, j, k, jzz )= ek( i, j, k, jzz )*Msun_geo

            detg= 2.0D0*g(i,j,k,jxy)*g(i,j,k,jxz)*g(i,j,k,jyz) &
                  - g(i,j,k,jzz)*g(i,j,k,jxy)**2 + g(i,j,k,jyy) &
                   *( g(i,j,k,jxx)*g(i,j,k,jzz) - g(i,j,k,jxz)**2 ) &
                  - g(i,j,k,jxx)*g(i,j,k,jyz)**2

            IF( ABS( detg ) < 1D-10 )THEN
              PRINT *, "** ERROR! import_id_spacetime: ", &
                       "The determinant of the spatial metric " &
                       // "is effectively 0 at the grid point " &
                       // "(i,j,k)= (", i, ",", j,",", k, "), " &
                       // "(x,y,z)= ", "(", &
                       pos( i, j, k, jx ),  ",", &
                       pos( i, j, k, jy ),  ",", &
                       pos( i, j, k, jz ), ")."
              PRINT *, "detg=", detg
              PRINT *
              PRINT *, "g_xx=", g(i,j,k,jxx)
              PRINT *, "g_xy=", g(i,j,k,jxy)
              PRINT *, "g_xz=", g(i,j,k,jxz)
              PRINT *, "g_yy=", g(i,j,k,jyy)
              PRINT *, "g_yz=", g(i,j,k,jyz)
              PRINT *, "g_zz=", g(i,j,k,jzz)
              STOP
            ELSEIF( detg < 0 )THEN
              PRINT *, "** ERROR! import_id_spacetime: ", &
                       "The determinant of the spatial metric " &
                       // "is negative at the grid point " &
                       // "(i,j,k)= (", i, ",", j,",", k, "), " &
                       // "(x,y,z)= ", "(", &
                       pos( i, j, k, jx ),  ",", &
                       pos( i, j, k, jy ),  ",", &
                       pos( i, j, k, jz ), ")."
              PRINT *, "detg=", detg
              PRINT *
              PRINT *, "g_xx=", g(i,j,k,jxx)
              PRINT *, "g_xy=", g(i,j,k,jxy)
              PRINT *, "g_xz=", g(i,j,k,jxz)
              PRINT *, "g_yy=", g(i,j,k,jyy)
              PRINT *, "g_yz=", g(i,j,k,jyz)
              PRINT *, "g_zz=", g(i,j,k,jzz)
              STOP
            ENDIF

            CALL compute_g4( i, j, k, lapse, shift, g, g4 )

            CALL determinant_sym4x4_grid( i, j, k, g4, detg4 )

            IF( ABS( detg4 ) < 1D-10 )THEN
              PRINT *, "** ERROR! import_id_spacetime: ", &
                       "The determinant of the spacetime metric "&
                       // "is effectively 0 at the grid point " &
                       // "(i,j,k)= (", i, ",", j,",", k, "), " &
                       // "(x,y,z)= ", "(", &
                       pos( i, j, k, jx ),  ",", &
                       pos( i, j, k, jy ),  ",", &
                       pos( i, j, k, jz ), ")."
              PRINT *, "detg4=", detg4
              PRINT *
              STOP
            ELSEIF( detg4 > 0 )THEN
              PRINT *, "** ERROR! import_id_spacetime: ", &
                       "The determinant of the spacetime metric "&
                       // "is positive at the grid point " &
                       // "(i,j,k)= (", i, ",", j,",", k, "), " &
                       // "(x,y,z)= ", "(", &
                       pos( i, j, k, jx ),  ",", &
                       pos( i, j, k, jy ),  ",", &
                       pos( i, j, k, jz ), ")."
              PRINT *, "detg4=", detg4
              PRINT *
              STOP
            ENDIF

            ! Print progress on screen
            perc= 100*( nx*ny*(k - 1) + nx*(j - 1) + i )/( nx*ny*nz )
            !perc2= 100.0*DBLE(nx*ny*(iz - 1) + nx*(iy - 1) + ix) &
            !       /DBLE( nx*ny*nz )
            !perc= 100*cnt/( nx*ny*nz )
            IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
              WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
                      creturn//" ", perc, "%"
              !WRITE( *, "(A2,F5.2,A1)", ADVANCE= "NO" ) &
              !        creturn//" ", perc2, "%"
            ENDIF

          ENDDO
        ENDDO
      ENDDO
      IF( show_progress ) WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

      PRINT *, "** Subroutine import_lorene_id executed."
      PRINT *

    ENDIF

  END PROCEDURE import_id_spacetime


  MODULE PROCEDURE import_id_hydro

    !*******************************************************
    !
    !# Stores the hydro ID in the arrays needed to compute
    !  the constraints on the refined mesh
    !
    !  FT 25.11.2020
    !
    !*******************************************************

    USE constants,  ONLY: Msun_geo
    USE tensor,     ONLY: jx, jy, jz

    IMPLICIT NONE

    INTEGER:: ix, iy, iz

    IF ( C_ASSOCIATED( THIS% diffstar_ptr ) ) THEN

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( nx, ny, nz, THIS, pos, &
      !$OMP                     baryon_density, energy_density, &
      !$OMP                     specific_energy, pressure, u_euler ) &
      !$OMP             PRIVATE( ix, iy, iz )
      coords_z: DO iz= 1, nz, 1
        coords_y: DO iy= 1, ny, 1
          coords_x: DO ix= 1, nx, 1

            ! The coordinates need to be converted from |sphincs| units (Msun_geo)
            ! to |lorene| units (\(\mathrm{km}\)). See MODULE constants for the definition of
            ! Msun_geo
            CALL get_diffstar_hydro( THIS% diffstar_ptr, &
                                     pos( ix, iy, iz, jx )*Msun_geo, &
                                     pos( ix, iy, iz, jy )*Msun_geo, &
                                     pos( ix, iy, iz, jz )*Msun_geo, &
                                     baryon_density( ix, iy, iz ), &
                                     energy_density( ix, iy, iz ), &
                                     specific_energy( ix, iy, iz ), &
                                     pressure( ix, iy, iz ), &
                                     u_euler( ix, iy, iz, jx ), &
                                     u_euler( ix, iy, iz, jy ), &
                                     u_euler( ix, iy, iz, jz ) )

          ENDDO coords_x
        ENDDO coords_y
      ENDDO coords_z
      !$OMP END PARALLEL DO

      !      ! Print progress on screen
      !      perc= 100*(nx*ny*(iz - 1) &
      !            + nx*(iy - 1) + ix)/( nx*ny*nz )
      !      IF( show_progress .AND. MOD( perc, 10 ) == 0 )THEN
      !        WRITE( *, "(A2,I2,A1)", ADVANCE= "NO" ) &
      !                creturn//" ", perc, "%"
      !      ENDIF
      !
      !    ENDDO coords_x
      !  ENDDO coords_y
      !ENDDO coords_z
      !IF( show_progress ) WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

      PRINT *, "** Subroutine import_lorene_id_hydro executed."
      PRINT *

    ENDIF

  END PROCEDURE import_id_hydro


  MODULE PROCEDURE import_id_particles

    !****************************************************
    !
    !# Stores the hydro ID in the arrays needed to
    !  compute the SPH ID
    !
    !  FT 19.11.2020
    !
    !****************************************************

    USE constants, ONLY: Msun_geo, km2m, g2kg, amu

    IMPLICIT NONE

    DOUBLE PRECISION:: detg

    IF ( C_ASSOCIATED( THIS% diffstar_ptr ) ) THEN

      IF( SIZE( x ) /= SIZE( y ) .OR. SIZE( x ) /= SIZE( z ) &
              .OR. SIZE( y ) /= SIZE( z ) )THEN
        PRINT *, "** ERROR: The sizes of the arrays of positions" &
                 // "passed to import_lorene_id are not the same."
        PRINT *
        STOP
      ENDIF

      PRINT *, "** Importing ID on particles..."

      !$OMP PARALLEL DO DEFAULT( NONE ) &
      !$OMP             SHARED( n, THIS, x, y, z, lapse, &
      !$OMP                     shift_x, shift_y, shift_z, &
      !$OMP                     g_xx, g_yy, g_zz, &
      !$OMP                     baryon_density, energy_density, &
      !$OMP                     specific_energy, pressure, &
      !$OMP                     u_euler_x, u_euler_y, u_euler_z ) &
      !$OMP             PRIVATE( itr )
      import_id_loop: DO itr= 1, n, 1

        ! The coordinates need to be converted from |sphincs| units (Msun_geo)
        ! to |lorene| units (\(\mathrm{km}\)). See MODULE constants for the definition of
        ! Msun_geo
        CALL get_diffstar_particles( THIS% diffstar_ptr, &
                                     x( itr )*Msun_geo, &
                                     y( itr )*Msun_geo, &
                                     z( itr )*Msun_geo, &
                                     lapse( itr ), &
                                     shift_x( itr ), &
                                     shift_y( itr ), &
                                     shift_z( itr ), &
                                     g_xx( itr ), &
                                     g_yy( itr ), &
                                     g_zz( itr ), &
                                     baryon_density( itr ), &
                                     energy_density( itr ), &
                                     specific_energy( itr ), &
                                     pressure( itr ), &
                                     u_euler_x( itr ), &
                                     u_euler_y( itr ), &
                                     u_euler_z( itr ) )

      ENDDO import_id_loop
      !$OMP END PARALLEL DO

      DO itr= 1, n, 1

        !
        !-- The following follows from the assumption of conformal
        !-- flatness in |lorene|
        !
        g_xy( itr )= 0.0D0
        g_xz( itr )= 0.0D0
        g_yz( itr )= 0.0D0

        !
        !- Set/unset the geodesic gauge
        !
        IF( THIS% get_one_lapse() )THEN
          lapse( itr )= 1.0D0
        ENDIF
        IF( THIS% get_zero_shift() )THEN
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

      ENDDO
      IF( show_progress ) WRITE( *, "(A1)", ADVANCE= "NO" ) creturn

      ! Convert the baryon density to units of amu (SPH code units)
      baryon_density= baryon_density*((Msun_geo*km2m)**3)/(amu*g2kg)

      PRINT *, "** Subroutine import_id_particles executed."
      PRINT *

    ENDIF

  END PROCEDURE import_id_particles


  MODULE PROCEDURE import_id_mass_b

    !****************************************************
    !
    !# Stores the hydro ID in the arrays needed to
    !  compute the baryon mass, storing it to variables
    !  (not arrays as the others SUBROUTINES in
    !  the [[bns_import]] SUBMODULE).
    !
    !  FT 15.04.2021
    !
    !****************************************************

    USE constants, ONLY: Msun_geo, lorene2hydrobase
    USE tensor,    ONLY: jxx, jxy, jxz, jyy, jyz, jzz

    IMPLICIT NONE

    IF ( C_ASSOCIATED( THIS% diffstar_ptr ) ) THEN

      ! The coordinates need to be converted from |sphincs| units (Msun_geo)
      ! to |lorene| units (\(\mathrm{km}\)).
      ! See MODULE constants for the definition of Msun_geo
      CALL get_diffstar_mass_b( THIS% diffstar_ptr, &
                                x*Msun_geo, &
                                y*Msun_geo, &
                                z*Msun_geo, &
                                g(jxx), &
                                g(jyy), &
                                g(jzz), &
                                baryon_density, &
                                gamma_euler )

      g(jxy)= 0.0D0
      g(jxz)= 0.0D0
      g(jyz)= 0.0D0

      baryon_density= baryon_density*lorene2hydrobase

    ENDIF

  END PROCEDURE import_id_mass_b


  MODULE PROCEDURE import_id_k

    !****************************************************
    !
    !# Stores the components of the extrinsic curvature
    !  in arrays
    !
    !  @warning DEPRECATED
    !
    !  FT 25.10.2021
    !
    !****************************************************


  END PROCEDURE import_id_k


  !-----------------!
  !--  FUNCTIONS  --!
  !-----------------!


  MODULE PROCEDURE import_mass_density

    !***********************************************
    !
    !# Returns the |lorene| mass density at the point
    !  given as argument, in units of
    !  \(M_\odot/L_\odot^3\).
    !
    !  FT
    !
    !***********************************************

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_ASSOCIATED
    USE constants,                   ONLY: Msun_geo, lorene2hydrobase

    IMPLICIT NONE

    IF ( C_ASSOCIATED( THIS% diffstar_ptr ) )THEN

      ! The coordinates need to be converted from |sphincs| units (Msun_geo)
      ! to |lorene| units (\(\mathrm{km}\)). See MODULE constants for the definition of
      ! Msun_geo
      res= get_diffstar_mass_density( THIS% diffstar_ptr, &
                                      x*Msun_geo, &
                                      y*Msun_geo, &
                                      z*Msun_geo )*lorene2hydrobase

    ENDIF

  END PROCEDURE import_mass_density


  MODULE PROCEDURE import_spatial_metric

    !***********************************************
    !
    !# Returns the |lorene| conformal factor to the
    !  4th power, equal to the diagonal components
    !  of the conformally flat spatial ADM metric.
    !
    !  FT 15.04.2021
    !
    !***********************************************

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_ASSOCIATED
    USE constants,                   ONLY: Msun_geo

    IMPLICIT NONE

    IF ( C_ASSOCIATED( THIS% diffstar_ptr ) )THEN

      ! The coordinates need to be converted from |sphincs| units (Msun_geo)
      ! to |lorene| units (\(\mathrm{km}\)). See MODULE constants for the definition of
      ! Msun_geo
      res= get_diffstar_spatial_metric( THIS% diffstar_ptr, &
                                        x*Msun_geo, &
                                        y*Msun_geo, &
                                        z*Msun_geo )

    ENDIF

  END PROCEDURE import_spatial_metric


  MODULE PROCEDURE is_hydro_negative

    !************************************************
    !
    !# Return 1 if the energy density is nonpositive
    !  or if the specific energy is nonpositive,
    !  or if the pressure is nonpositive
    !  at the specified point
    !
    !  FT 25.10.2021
    !
    !************************************************

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_ASSOCIATED
    USE constants,                   ONLY: Msun_geo

    IMPLICIT NONE

    IF ( C_ASSOCIATED( THIS% diffstar_ptr ) )THEN

      ! The coordinates need to be converted from |sphincs| units (Msun_geo)
      ! to |lorene| units (\(\mathrm{km}\)). See MODULE constants for the definition of
      ! Msun_geo
      res= negative_hydro( THIS% diffstar_ptr, &
                           x*Msun_geo, &
                           y*Msun_geo, &
                           z*Msun_geo )

    ENDIF

  END PROCEDURE is_hydro_negative


END SUBMODULE diffstar_lorene_import
