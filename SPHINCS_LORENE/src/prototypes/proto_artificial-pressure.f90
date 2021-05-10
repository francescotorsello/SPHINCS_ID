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

  USE bns_id,         ONLY: bns
  !USE particles_id,   ONLY: particles
  USE RCB_tree_3D,    ONLY: allocate_RCB_tree_memory_3D,&
                            iorig
  USE kernel_table,   ONLY: ktable
  USE sph_variables,  ONLY: Rstar,divv,av,Pr,ye,temp,nlrf,&
                            u,Theta,vel_u,tterm,tgrav,tkin,&
                            escap,t,n1,n2,pos_u,h,nu,npart,&
                            npm,Nstar,allocate_sph_memory
  USE options,        ONLY: basename,ikernel,ndes,av_max,icow
  USE input_output,   ONLY: read_options
  USE units,          ONLY: umass,m0c2_CU,set_units

  USE APM,            ONLY: density_loop,position_correction,&
                            assign_h
  USE set_h,          ONLY: exact_nei_tree_update
  USE analyze,        ONLY: COM

  IMPLICIT NONE

  INTEGER, PARAMETER:: unit_id  = 23
  INTEGER, PARAMETER:: max_npart= 5D+6
  INTEGER:: npart_tmp, tmp, ios, itr, a, nout

  DOUBLE PRECISION com_x,com_y,com_z,com_d
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
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nu0
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nlrf0
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: Ye0
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: Theta0
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: sph_density
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: nstar0
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: h0

  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: Nstar_real

  LOGICAL:: exist

  CHARACTER( LEN= : ), ALLOCATABLE:: namefile
  CHARACTER( LEN= : ), ALLOCATABLE:: finalnamefile
  CHARACTER( LEN= : ), ALLOCATABLE:: err_msg

  TYPE( bns ):: binary

  !---------------------------!
  !--  End of declarations  --!
  !---------------------------!

  !CALL DATE_AND_TIME( date, time, zone, values )
  !run_id= date // "-" // time


  !PRINT *, "Hello world!"
  PRINT *, "** Beginning of PROGRAM proto_apm."
  PRINT *

  namefile= "poly2-75_1.2-1.8_45km.bin"

  binary= bns( namefile )

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
  ALLOCATE( nu0(max_npart) )
  ALLOCATE( nlrf0(max_npart) )
  ALLOCATE( Ye0(max_npart) )
  ALLOCATE( Theta0(max_npart) )
  ALLOCATE( sph_density(max_npart) )
  ALLOCATE( nstar0(max_npart) )
  ALLOCATE( h0(max_npart) )

  finalnamefile= "lorene-bns-id-particles-shells.dat"

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
  npart_tmp= 0
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
      nu0( itr ), &
      nlrf0( itr ), &
      Ye0( itr ), &
      Theta0( itr ), &
      sph_density( itr ), &
      nstar0( itr ), &
      h0( itr )
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
    npart_tmp= npart_tmp + 1
  ENDDO
  IF( tmp /= npart_tmp )THEN
    PRINT *, "** ERROR! Mismatch in the values of the particle number read", &
             " from " // TRIM(finalnamefile), "in two ways. "
    PRINT *, "tmp= ", tmp, ", npart_tmp= ", npart_tmp
    PRINT *, "Stopping..."
    STOP
  ENDIF

  CLOSE( unit_id )

  pos                  = pos                  ( :, 1:npart_tmp )
  lapse_parts          = lapse_parts          ( 1:npart_tmp )
  shift_parts_x        = shift_parts_x        ( 1:npart_tmp )
  shift_parts_y        = shift_parts_y        ( 1:npart_tmp )
  shift_parts_z        = shift_parts_z        ( 1:npart_tmp )
  baryon_density_parts = baryon_density_parts ( 1:npart_tmp )
  energy_density_parts = energy_density_parts ( 1:npart_tmp )
  specific_energy_parts= specific_energy_parts( 1:npart_tmp )
  pressure_parts       = pressure_parts       ( 1:npart_tmp )
  v_euler_parts_x      = v_euler_parts_x      ( 1:npart_tmp )
  v_euler_parts_y      = v_euler_parts_y      ( 1:npart_tmp )
  v_euler_parts_z      = v_euler_parts_z      ( 1:npart_tmp )
  v                    = v                    ( :, 1:npart_tmp )
  nu0                  = nu0                  ( 1:npart_tmp )
  nlrf0                = nlrf0                ( 1:npart_tmp )
  Ye0                  = Ye0                  ( 1:npart_tmp )
  Theta0               = Theta0               ( 1:npart_tmp )
  sph_density          = sph_density          ( 1:npart_tmp )
  nstar0               = nstar0               ( 1:npart_tmp )
  h0                   = h0                   ( 1:npart_tmp )

  !-------------------------------!
  !-- Needed calls from SPHINCS --!
  !-------------------------------!

  ! setup unit system
  CALL set_units('NSM')

  !-------------------------------!
  !-- read main options for run --!
  !-------------------------------!
  CALL read_options

  npart= npart_tmp
  CALL allocate_SPH_memory

  pos_u                = pos                  ( :, 1:npart_tmp )
  vel_u                = v                    ( :, 1:npart_tmp )
  h                    = h0                   ( 1:npart_tmp )

  CALL allocate_RCB_tree_memory_3D(npart)
  iorig(1:npart)= (/ (a,a=1,npart) /)

  !CALL allocate_metric_on_particles(npart)
  ! for now: just constant
  !av(1:npart)= av_max

  ! tabulate kernel, get ndes
  CALL ktable(ikernel,ndes)

  !-----------------------------------------------------------------------!
  !-- At this point I can apply the artificial pressure method (can I?) --!
  !-----------------------------------------------------------------------!

  ! In setup_tov_star distribute_and_stretch is called
  ! In distribute_and_stretch, setup_uniform_sphere is called
  ! setup_uniform_sphere first place positions in a cubic lattice or a
  ! close-packed lattice, then assigns the smoothing length

  PRINT *, "npart_tmp=", npart_tmp
  PRINT *
  CALL assign_h(npart_tmp,npart,pos_u,h0,h)

  !STOP

  ! measure SPH-particle number density
  ALLOCATE( Nstar_real(max_npart) )
  nu= 1.0D0
  CALL density_loop(npart,pos_u,nu,h,Nstar)

  ! In setup_uniform_sphere, get_profile_density is called
  ! this computed nstar, but we have it from the ID file

  ! assign nu's
  nu= nstar/Nstar

  ! In setup_uniform_sphere, reset_COM is called
  !CALL reset_COM(npart,nu,pos_u)
  CALL COM(npart,pos_u,nu,com_x,com_y,com_z,com_d)
  DO a=1,npart
     pos_u(1,a)= pos_u(1,a) - com_x
     pos_u(2,a)= pos_u(2,a) - com_y
     pos_u(3,a)= pos_u(3,a) - com_z
  ENDDO

  ! at this point, the particle should be mirrored (?)
  ! TODO: distinguish the two stars and write the prototype to handle one star
  !       at the time

  ! Re-estimate nu
  CALL density_loop(npart,pos_u,nu,h,Nstar)

  ! Here you should get the LORENE ensity on the new positions...this means you
  ! need to compile LORENE and SPHINCS_LORENE...however, you might be able to
  ! use the SUBROUTINES bound to LORENE without creating a bns object..
  ! no, you need a bns object

  PRINT *, "** End of PROGRAM proto_apm."
  PRINT *

  CONTAINS


END PROGRAM proto_apm
