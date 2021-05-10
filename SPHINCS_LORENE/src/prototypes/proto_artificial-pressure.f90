! File:         proto_artificial-pressure.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

PROGRAM proto_apm

  !*****************************************************
  !                                                    *
  ! Prototype for the development of the artificial    *
  ! pressure method for the LORENE ID                  *
  !                                                    *
  ! FT 10.05.2021                                      *
  !                                                    *
  !*****************************************************

  !USE particles_id,   ONLY: particles

  IMPLICIT NONE

  INTEGER, PARAMETER:: unit_id  = 23
  INTEGER, PARAMETER:: max_npart= 5D+6
  INTEGER:: npart, tmp, ios, itr

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: pos
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: lapse_parts
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: shift_parts_x
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: shift_parts_y
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: shift_parts_z
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: baryon_density_parts
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: energy_density_parts
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: specific_energy_parts
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: pressure_parts
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: v_euler_parts_x
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: v_euler_parts_y
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: v_euler_parts_z
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: v
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nu
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nlrf
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: Ye
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: Theta
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: sph_density
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: h

  LOGICAL:: exist

  CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile
  CHARACTER( LEN= : ), ALLOCATABLE:: err_msg

  !---------------------------!
  !--  End of declarations  --!
  !---------------------------!

  !CALL DATE_AND_TIME( date, time, zone, values )
  !run_id= date // "-" // time


  !PRINT *, "Hello world!"
  PRINT *, "** Beginning of PROGRAM proto_apm."
  PRINT *

  ALLOCATE( pos(3,max_npart) )
  ALLOCATE( lapse_parts(max_npart) )
  ALLOCATE( shift_parts_x(max_npart) )
  ALLOCATE( shift_parts_y(max_npart) )
  ALLOCATE( shift_parts_z(max_npart) )
  ALLOCATE( baryon_density_parts(max_npart) )
  ALLOCATE( energy_density_parts(max_npart) )
  ALLOCATE( specific_energy_parts(max_npart) )
  ALLOCATE( pressure_parts(max_npart) )
  ALLOCATE( v_euler_parts_x(max_npart) )
  ALLOCATE( v_euler_parts_y(max_npart) )
  ALLOCATE( v_euler_parts_z(max_npart) )
  ALLOCATE( v(3,max_npart) )
  ALLOCATE( nu(max_npart) )
  ALLOCATE( nlrf(max_npart) )
  ALLOCATE( Ye(max_npart) )
  ALLOCATE( Theta(max_npart) )
  ALLOCATE( sph_density(max_npart) )
  ALLOCATE( nstar(max_npart) )
  ALLOCATE( h(max_npart) )

  finalnamefile= "lorene-bns-id-particles.dat"

  INQUIRE( FILE= TRIM(finalnamefile), EXIST= exist )

  IF( exist )THEN
    OPEN( UNIT= unit_id, FILE= TRIM(finalnamefile), &
          FORM= "FORMATTED", ACTION= "READ", IOSTAT= ios, &
          IOMSG= err_msg )
    IF( ios > 0 )THEN
      PRINT *, "...error when opening " // TRIM(finalnamefile), &
              ". The error message is", err_msg
      STOP
    ENDIF
  ELSE
    PRINT *, "** ERROR! Unable to find file " // TRIM(finalnamefile)
    STOP
  ENDIF

  PRINT *, " * Reading file " // TRIM(finalnamefile) // "..."
  npart= 0
  ! Skipping the first 4 header rows
  READ( unit_id, * )
  READ( unit_id, * )
  READ( unit_id, * )
  READ( unit_id, * )
  ! Reading particle data
  DO itr= 1, max_npart, 1
    READ( UNIT= unit_id, FMT= *, IOSTAT = ios, IOMSG= err_msg ) &
      tmp, &
      pos( 1, itr ), &
      pos( 2, itr ), &
      pos( 3, itr ), &
      lapse_parts( itr ), &
      shift_parts_x( itr ), &
      shift_parts_y( itr ), &
      shift_parts_z( itr ), &
      baryon_density_parts( itr ), &
      energy_density_parts( itr ), &
      specific_energy_parts( itr ), &
      pressure_parts( itr ), &
      v_euler_parts_x( itr ), &
      v_euler_parts_y( itr ), &
      v_euler_parts_z( itr ), &
      v( 1, itr ), &
      v( 2, itr ), &
      v( 3, itr ), &
      nu( itr ), &
      nlrf( itr ), &
      Ye( itr ), &
      Theta( itr ), &
      sph_density( itr ), &
      nstar( itr ), &
      h( itr )
    IF( ios > 0 )THEN
      PRINT *, "...error when reading " // TRIM(finalnamefile), &
              " at step ", itr,". The status variable is ", ios, &
              ". The error message is", err_msg
      STOP
    ENDIF
    IF( ios < 0 )THEN
      PRINT *, " * Reached end of file " // TRIM(finalnamefile)
      PRINT *
      EXIT
    ENDIF
    npart= npart + 1
  ENDDO

  CLOSE( unit_id )

  pos                  = pos                  ( :, 1:npart )
  lapse_parts          = lapse_parts          ( 1:npart )
  shift_parts_x        = shift_parts_x        ( 1:npart )
  shift_parts_y        = shift_parts_y        ( 1:npart )
  shift_parts_z        = shift_parts_z        ( 1:npart )
  baryon_density_parts = baryon_density_parts ( 1:npart )
  energy_density_parts = energy_density_parts ( 1:npart )
  specific_energy_parts= specific_energy_parts( 1:npart )
  pressure_parts       = pressure_parts       ( 1:npart )
  v_euler_parts_x      = v_euler_parts_x      ( 1:npart )
  v_euler_parts_y      = v_euler_parts_y      ( 1:npart )
  v_euler_parts_z      = v_euler_parts_z      ( 1:npart )
  v                    = v                    ( :, 1:npart )
  nu                   = nu                   ( 1:npart )
  nlrf                 = nlrf                 ( 1:npart )
  Ye                   = Ye                   ( 1:npart )
  Theta                = Theta                ( 1:npart )
  sph_density          = sph_density          ( 1:npart )
  nstar                = nstar                ( 1:npart )
  h                    = h                    ( 1:npart )

  PRINT *, "** End of PROGRAM proto_apm."
  PRINT *

  CONTAINS


END PROGRAM proto_apm
