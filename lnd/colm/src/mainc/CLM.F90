#include <define.h>

#ifndef CPL7

  PROGRAM CoLM
! ======================================================================
! The Common Land Model was developed in cooperation with
!     Beijing Normal University                             (Dai)
!     Georgia Institute of Technology                       (Dickinson)
!     National Center for Atmospheric Research              (Bonan, Oleson)
!     University of Arizona                                 (Zeng)
!     University of Texas at Austin                         (Yang)
!     GSFC/NASA                                             (Houser, Bosilovich)
!     COLA                                                  (Dirmeyer, Schlosser)
!     Colorado State University at Fort Collins             (Denning, Baker)
!
! Reference: 
!     [1] Dai et al., 2003: The Common Land Model (CLM). 
!         Bull. of Amer. Meter. Soc., 84: 1013-1023
!     [2] Dai et al., 2004: A two-big-leaf model for canopy temperature,
!         photosynthesis and stomatal conductance. J. Climate, 17:2281-2299
!
!     Author: Yongjiu Dai, January 2004
! ======================================================================

      use precision
      use spmd
      use spmd_decomp
      use timemgr
      use forcedata
      use colm_cplMod
      use colm_csmMod
      use colm_rtmMod
      use colm_ioMod
      use colm_varMod
      use colm_varctl
      use landuse
      use nchistMod
      use shr_msg_mod        ! CSM message passing routines and variables 
      use shr_sys_mod        ! CSM system module

      implicit none

!local variables:

      logical :: lwrite                          ! true: write out frequency
      logical :: doalb_lnd                       ! true if time for surface albedo calculation
      logical :: dolai_lnd                       ! true if time for time-varying vegetation paramter
      logical :: dosst_lnd                       ! true if time for update sst/ice/snow

      integer :: i,j,k,l,m                       ! looping indices

!ROUTINE:

    ! Determine input/output units 
      CALL shr_msg_stdio ('lnd')

#ifdef COUP_CSM
    ! Initialize MPI communication groups for flux coupler
      CALL csm_setup(p_comm)
#endif

#ifdef SPMD 
    ! Initialize MPI spmd environment
      CALL spmd_init
#endif

      if (p_master) then 
         write(6,*) '--------Now SPMD is OK---------'
         call shr_sys_flush(6) 
      endif 

    ! Read namelist
      CALL readnml

      CALL colm_var_alloc1

    ! Read surface data, to construct grid info
      CALL readgridat

    ! MPI decomposition
      CALL task_decomp

#ifdef RTM
    ! Initialize RTM model
      CALL colm_rtm_init
#endif

#ifdef COUP_CSM
    ! land-atmos flux initialize
      CALL colm_cpl_init
#endif

    ! Initialize time-constant and time-varying variables
      CALL colm_var_alloc2

      CALL readinidat

      CALL colm_var_dealloc1

    ! Initialize land use module

      CALL landuse_init

#ifdef COUP_CSM
    ! Initialize flux coupler communication
      CALL csm_initialize(irad, eccen, obliqr, lambm0, mvelpp)

    ! advance one step to continue run
      if(nsrest.gt.0) then
         CALL TICKTIME
         istep = istep + 1
      endif

    ! Send first land model data to flux coupler.
      CALL csm_sendalb

      if (p_master) then 
         write(6,*) ' CoLM Main -------- irad: ', irad
         call shr_sys_flush(6) 
      endif 
#else
    ! Open forcing data for offline run
      CALL open_forcedata

    ! Calendar for NEXT time step (offline)
      CALL TICKTIME
      istep = istep + 1
#endif

    ! Initialize netCDF history output
      CALL nchist_init

! ======================================================================
! begin time stepping loop
! ======================================================================

      DO

       ! doalb_lnd is true when the next time step is a radiation time step
         doalb_lnd = .true.
         dolai_lnd = .true.
         dosst_lnd = .false.

#ifdef COUP_CSM
       ! doalb is true when the next time step is a radiation time step
       ! this allows for the fact that an atmospheric model may not do
       ! the radiative calculations every time step. for example:
       !      nstep dorad doalb
       !        1     F     F
       !        2     F     T
       !        3     T     F
       ! The following expression for doalb is specific to CAM

         doalb = ((irad==1 .and. istep+1/=1) .or. (mod(istep,irad)==0 .and. istep+1/=1))

       ! Determine if information should be sent/received to/from flux coupler
         CALL csm_dosndrcv(doalb)
    
       ! Get atmospheric state and fluxes from flux coupler
         if (dorecv) then
            CALL csm_recv
            if (csmstop_now) EXIT
         endif

         CALL colm_cpl_a2l
#else
         if(istep.gt.mstep) EXIT

       ! Read forcing data for offline run
         CALL read_forcedata(dolai_lnd)
#endif

         oro(:) = 1.

       ! Call clm driver
         CALL CLMDRIVER(dolai_lnd,doalb_lnd,dosst_lnd)
 
       ! Mapping subgrid patch [numpatch] vector of subgrid points to grid average
         CALL FLUXAVE

#ifdef DGVM
         CALL LPJDRIVER
#endif

#ifdef RTM
         CALL colm_rtm_drv
#endif

#ifdef COUP_CSM
       ! Average fluxes over interval if appropriate
       ! Surface states sent to the flux coupler states are not time averaged

         CALL colm_cpl_l2a

         if (csm_doflxave) CALL csm_flxave
       
       ! Send fields to flux coupler
       ! Send states[n] (except for snow[n-1]), time averaged fluxes for [n,n-1,n-2],
       ! albedos[n+1], and ocnrof_vec[n]

         if (dosend) CALL csm_send
#endif

#ifdef MYBUG
         if (p_master) then
            write(6,'(a, a10, i10, a8, i4.4, a8, i3.3, a8, i5.5)' )  &
                    '--------------------CoLM Main Count : Steps & Time : ', &
                    'Nstep=', istep, 'Year=', idate(1), 'Days=', idate(2), 'Secs=', idate(3) 

            call shr_sys_flush(6)
         endif
#endif
 
         CALL writehistdat
 
       ! Calendar for NEXT time step
         CALL TICKTIME

         istep = istep + 1

      END DO

! ======================================================================
! end of time stepping loop
! ======================================================================

#ifdef COUP_CSM
      CALL csm_shutdown

      CALL colm_cpl_exit
#else
      CALL close_forcedata
#endif
     
      CALL colm_var_dealloc2

#ifdef RTM
      CALL colm_rtm_exit
#endif

      CALL landuse_exit

      CALL nchist_exit

#ifdef SPMD
      CALL spmd_var_dealloc

      CALL spmd_exit
#endif

      write(6,*) 'CoLM Execution Completed'

  END PROGRAM CoLM

#endif
