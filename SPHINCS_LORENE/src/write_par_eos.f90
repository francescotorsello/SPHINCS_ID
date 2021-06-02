! File:         write_par_eos.f90
! Authors:      Francesco Torsello (FT)
! Copyright:    GNU General Public License (GPLv3)

PROGRAM write_par_eos

  USE constants, ONLY: lorene2hydrobase, kg2g, m2cm, &
                       k_lorene2hydrobase_piecewisepolytrope
  USE pwp_EOS,   ONLY: Gamma0, K0, select_EOS_parameters, gen_pwp_eos, &
                       get_rho_0, get_Gamma1, get_Gamma2, get_Gamma3, &
                       get_K1, get_K2, get_K3

  IMPLICIT NONE

  INTEGER:: ios
  LOGICAL:: exist
  CHARACTER( LEN= : ), ALLOCATABLE:: namefile
  CHARACTER( LEN= : ), ALLOCATABLE:: err_msg
  CHARACTER( LEN= 4 ):: eos

  WRITE(*,'("Please write a 4 character string containing the name of the piecewise polytropic EOS: ",\)')
  READ(*,'(A)') eos

  CALL select_EOS_parameters(eos)

  namefile= "par_eos.d"

  INQUIRE( FILE= TRIM(namefile), EXIST= exist )

  IF( exist )THEN
      OPEN( UNIT= 2, FILE= TRIM(namefile), STATUS= "REPLACE", &
            FORM= "FORMATTED", &
            POSITION= "REWIND", ACTION= "WRITE", IOSTAT= ios, &
            IOMSG= err_msg )
  ELSE
      OPEN( UNIT= 2, FILE= TRIM(namefile), STATUS= "NEW", &
            FORM= "FORMATTED", &
            ACTION= "WRITE", IOSTAT= ios, IOMSG= err_msg )
  ENDIF
  IF( ios > 0 )THEN
    PRINT *, "...error when opening " // TRIM(namefile), &
             ". The error message is", err_msg
    STOP
  ENDIF

  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "110        Type of the EOS (cf. documentation of Eos::eos_multi_poly)"

  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "Multipolytropic ", eos, " EOS"

  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "4    		npoly,	  number of polytropes"

  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    Gamma0, "     gamma_0,  crust (here from SLy)"
  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    get_Gamma1(), "		gamma_1,  array of adiabatic indexes (from crust to core)"
  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    get_Gamma2(), "		gamma_2"
  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    get_Gamma3(), "		gamma_3"
  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    K0/k_lorene2hydrobase_piecewisepolytrope(Gamma0), &
    "     kappa0,     pressure coefficient for the crust (here from SLy)"
  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    LOG10( K0/k_lorene2hydrobase_piecewisepolytrope(Gamma0) &
           *( get_rho_0()/lorene2hydrobase*kg2g/(m2cm**3.0D0) )**(Gamma0) ), &
    "     log10(P0/c^2),	log of pressure between gam_0 and gam_1 (dyne/cm^2/c_cgs^2)"

  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    LOG10( get_rho_0()/lorene2hydrobase*kg2g/(m2cm**3.0D0) ), &
    "     log10(rho_0),	array of the exponent of fiducial densities (g/cm^3)"
  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "14.7       log10(rho_1)"
  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "15.0       log10(rho_2)"

  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "0.         decInc,   array (size npeos - 1) of the percentage which "
  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "0.                   determines the terminating enthalpy for lower "
  WRITE( UNIT = 2, IOSTAT = ios, IOMSG = err_msg, FMT = * ) &
    "0.                   density and the starting enthalpy for higher density."

  CLOSE( UNIT= 2 )

END PROGRAM write_par_eos
