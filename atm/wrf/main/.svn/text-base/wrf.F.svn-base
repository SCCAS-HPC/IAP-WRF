!WRF:DRIVER_LAYER:MAIN
!

PROGRAM wrf

#ifdef CCSMCOUPLED
#ifdef SEQ_MCT
#else
#endif
#else
   USE atm_comp_mct, only : atm_init_mct, atm_dfi, atm_run_mct, atm_final_mct

!<DESCRIPTION>
! Main program of WRF model.  Responsible for starting up the model, reading in (and
! broadcasting for distributed memory) configuration data, defining and initializing
! the top-level domain, either from initial or restart data, setting up time-keeping, and
! then calling the <a href=integrate.html>integrate</a> routine to advance the domain
! to the ending time of the simulation. After the integration is completed, the model
! is properly shut down.
!
!</DESCRIPTION>

   IMPLICIT NONE

#ifdef _OPENMP
  CALL setfeenv()
#endif

  ! Set up WRF model.  
  CALL atm_init_mct

  ! Run digital filter initialization if requested.
  CALL atm_dfi_mct

  ! WRF model time-stepping.  Calls integrate().  
  CALL atm_run_mct

  ! WRF model clean-up.  This calls MPI_FINALIZE() for DM parallel runs.  
  CALL atm_final_mct
#endif
END PROGRAM wrf


