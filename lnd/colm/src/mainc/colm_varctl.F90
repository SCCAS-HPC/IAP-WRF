#include <define.h>

module colm_varctl

   use precision
   implicit none

   integer :: irad              ! solar radiation frequency (iterations)
   logical :: wrtdia            ! true => write global average diagnostics to std out
   logical :: csm_doflxave      ! true => only communicate with flux coupler on albedo calc time steps
   logical :: doalb             ! albedo 
   integer :: nsrest            ! nsrest=0: initial run; nsrest=1: continual run

   integer, public, parameter :: nsrStartup  = 0                          ! Startup from initial conditions
   integer, public, parameter :: nsrContinue = 1                          ! Continue from restart files
   integer, public, parameter :: nsrBranch   = 2                          ! Branch from restart files

   logical, public :: create_glacier_mec_landunit = .false.               ! glacier_mec landunit is not created

   logical :: debug = .false.     ! control debug information output

   logical :: lcycdump  = .false. ! logical flag do enable cycdump function.
   integer :: ncycdump  = 96      ! number of steps to dump
   integer :: lucycdump = 100     ! logical unit of fcycdump
   character(len=255) :: fcycdump = "./CoLM-cycdump"  ! file name of cycdump (restart & forcing)

! history output options
   integer, parameter :: nMaxHist = 6

   logical :: lhist_yearly  = .false. ! yearly history files
   logical :: lhist_monthly = .true.  ! monthly history files
   logical :: lhist_daily   = .false. ! daily history files
   logical :: lhist_3hourly = .false. ! 3hourly history files
   logical :: lhist_hourly  = .false. ! hourly history files
   logical :: lhist_steply  = .false. ! steply history files

#ifdef RTM
   integer :: rtm_dtime = 3600*3      ! time interval(seconds) to call RTM
#endif

#ifdef SPINUP
   integer :: spinup_curryr = 1000
   integer :: spinup_begdat = 1990
   integer :: spinup_enddat = 1990
#endif

#ifdef COUP_CSM
   integer :: lnd_cflux_year = 0.               ! year to turn on carbon flux exchanging 
!!! Only for CPL6
   character(len=31) :: co2_option = "co2prog"  ! 'co2prog','co2diag'
!!! End for CPL6
#endif

   real(r8) :: po2  = 0.209                     ! constant atmospheric partial pressure  O2 (mol/mol)
   real(r8) :: pco2 = 355.0E-06                 ! constant atmospheric partial pressure CO2 (mol/mol)
!  real(r8) :: pco2 = 284.725E-06               ! constant atmospheric partial pressure CO2 (mol/mol)

#ifdef CPL7
!!! Only for CPL7
   real(r8) :: co2_ppmv = 367.0                 ! atmospheric CO2 molar ratio (by volume) (umol/mol)
   character(len=16) :: co2_type = 'constant'   ! values of 'prognostic','diagnostic','constant'
!!! End for CPL7
#endif

   logical :: lcarbon_emission = .false.

end module colm_varctl
