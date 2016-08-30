module seq_comm_mct
!---------------------------------------------------------------------
!
! Purpose: MCT utitlity functions used in sequential CCSM.
!          Note that if no MPI, will call MCTs fake version 
!          (including mpif.h) will be utilized
!       
! Author: R. Jacob
!
!---------------------------------------------------------------------

  use mct_mod, only : mct_world_init, mct_die
  use shr_sys_mod, only : shr_sys_abort, shr_sys_flush
  use shr_mpi_mod, only : shr_mpi_chkerr, shr_mpi_bcast, shr_mpi_max
  use shr_file_mod, only : shr_file_getUnit, shr_file_freeUnit

  implicit none
  save
  private

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public seq_comm_init
  public seq_comm_iamin
  public seq_comm_iamroot
  public seq_comm_mpicom
  public seq_comm_iam
  public seq_comm_gloiam
  public seq_comm_gloroot
  public seq_comm_cplpe
  public seq_comm_cmppe
  public seq_comm_setptrs
  public seq_comm_setnthreads
  public seq_comm_getnthreads
  public seq_comm_printcomms

!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

  integer, public :: logunit  = 6     ! log unit number
  integer, public :: loglevel = 1     ! log level

!  NOTE: the LNDID must be 1 as this is currently hardwired in
!  clm for lnd/rtm mapping.  This constraint can be removed when
!  clm gets the compid from the driver.

  integer, parameter, public :: ncomps = 19   ! by juanxiong he for wrf/cam coupling 
                               ! by huiqun hao for dc
  integer, parameter, public :: LNDID = 1
  integer, parameter, public :: ATMID = 2
  integer, parameter, public :: OCNID = 3
  integer, parameter, public :: ICEID = 4
  integer, parameter, public :: GLCID = 5
  integer, parameter, public :: CPLID = 6
  integer, parameter, public :: GLOID = 7
  integer, parameter, public :: CPLATMID = 8
  integer, parameter, public :: CPLLNDID = 9
  integer, parameter, public :: CPLICEID = 10
  integer, parameter, public :: CPLOCNID = 11
  integer, parameter, public :: CPLGLCID = 12
  integer, parameter, public :: WRFID = 13    ! by juanxiong he for wrf/cam coupling
  integer, parameter, public :: CPLWRFID = 14  ! by juanxiong he for wrf/cam coupling
  integer, parameter, public :: GEAID = 15    ! by juanxiong he for wrf/camcoupling
  integer, parameter, public :: CPLGEAID = 16  ! by juanxiong he for wrf/cam coupling
  integer, parameter, public :: SRDID = 17  ! by huiqun hao for dc 
  integer, parameter, public :: SRDATMID = 18  ! by huiqun hao for dc 
  integer, parameter, public :: SRDWRFID = 19  ! by huiqun hao for dc 

  character(len=8),parameter,private :: IDname(ncomps) = &
    (/ '   LND  ','   ATM  ','   OCN  ','   ICE  ','   GLC  ', &
       '   CPL  ',' GLOBAL ',' CPLATM ',' CPLLND ',' CPLICE ', &
       ' CPLOCN ',' CPLGLC ','   WRF  ',' CPLWRF ','   GEA  ', &
       ' CPLGEA ','   SRD  ',' SRDATM ',' SRDWRF '/)   ! by juanxiong he for wrf/cam coupling
	                                                   ! by huiqun hao for dc
 
  type seq_comm_type
    character(len=8) :: name   ! my name, see IDname above
    integer :: ID              ! my id number, see parameters above
    integer :: mpicom          ! mpicom
    integer :: mpigrp          ! mpigrp
    integer :: npes            ! number of pes in comm
    integer :: nthreads        ! number of omp threads per pe
    integer :: iam             ! my pe number in mpicom
    logical :: iamroot         ! am i the root pe in mpicom
    integer :: gloiam          ! my pe number in mpi_comm_world
    integer :: gloroot         ! the global pe number of each comps root on all pes
    integer :: pethreads       ! max number of threads on my pe
    integer :: cplpe           ! a common pe in mpicom from the cpl group for join mpicoms
    integer :: cmppe           ! a common pe in mpicom from the component group for join mpicoms
    integer :: mpicomcart      ! mpicomcart for wrf, added by juanxiong he
    integer :: mpicomcart_periodic    ! mpicomcart_periodic, added by juanxiong he
    integer :: mpicomx         ! mpicomx for wrf,    added by juanxiong he
    integer :: mpicomy         ! mpicomy for wrf,    added by juanxiong he
    integer :: ntasks_x        ! ntasks_x for wrf,   added by juanxiong he
    integer :: ntasks_y        ! ntasks_y for wrf,   added by juanxiong he
    integer :: gmpicomcart      ! mpicomcart for wrf, added by juanxiong he
    integer :: gmpicomcart_periodic    ! mpicomcart_periodic, added by juanxiong he
    integer :: gmpicomx         ! mpicomx for wrf,    added by juanxiong he
    integer :: gmpicomy         ! mpicomy for wrf,    added by juanxiong he
    integer :: gntasks_x        ! ntasks_x for wrf,   added by juanxiong he
    integer :: gntasks_y        ! ntasks_y for wrf,   added by juanxiong he

    logical :: set             ! has this datatype been set
  end type seq_comm_type

  type(seq_comm_type) :: seq_comms(ncomps)

  character(*),parameter :: F11 = "(a,a,'(',i3,a,')',a,   3i6,' (',a,i6,')',' (',a,i3,')')"
  character(*),parameter :: F12 = "(a,a,'(',i3,a,')',a,2i6,6x,' (',a,i6,')',' (',a,i3,')','(',a,2i6,')')"
  character(*),parameter :: F13 = "(a,a,'(',i3,a,')',a,2i6,6x,' (',a,i6,')',' (',a,i3,')')"
  integer :: Global_Comm
#include <mpif.h>  

!=======================================================================
contains
!=======================================================================
  
  subroutine seq_comm_init(Comm_in, nmlfile, atm_petlist, lnd_petlist, ice_petlist, ocn_petlist, glc_petlist, &
                            wrf_petlist,gea_petlist,srd_petlist) !juanxiong he
    !----------------------------------------------------------
    !
    ! Arguments
    implicit none
    integer, intent(in) :: Comm_in
    character(len=*), intent(IN) :: nmlfile
    integer, pointer, optional             :: atm_petlist(:)
    integer, pointer, optional             :: lnd_petlist(:)
    integer, pointer, optional             :: ice_petlist(:)
    integer, pointer, optional             :: ocn_petlist(:)
    integer, pointer, optional             :: glc_petlist(:)
    integer, pointer, optional             :: wrf_petlist(:) !juanxiong he
    integer, pointer, optional             :: gea_petlist(:) !juanxiong he
	integer, pointer, optional             :: srd_petlist(:) !huiqun hao
    !
    ! Local variables
    !
    integer :: ierr,n
    character(*),parameter :: subName =   '(seq_comm_init) '
    integer :: mpi_group_world   ! GLOBAL_COMM group
    integer :: mype,numpes,myncomps,max_threads,gloroot
    integer :: amin,amax,astr
    integer :: lmin,lmax,lstr
    integer :: imin,imax,istr
    integer :: omin,omax,ostr
    integer :: gmin,gmax,gstr
    integer :: cmin,cmax,cstr
    integer :: wmin,wmax,wstr    ! by juanxiong he for wrf/cam coupling
    integer :: gemin,gemax,gestr    ! by juanxiong he for wrf/cam coupling
	integer :: smin,smax,sstr    ! by huiqun hao for dc
    integer :: pelist(3,1)       ! start, stop, stride for group
    integer, pointer :: comps(:) ! array with component ids
    integer, pointer :: comms(:) ! array with mpicoms
    integer :: onecomm           ! single comm for "old" mct init
    integer :: nu, i
    logical,save :: first_pass = .true.   ! 
    integer :: &
         atm_ntasks, atm_rootpe, atm_pestride, atm_nthreads, &
         lnd_ntasks, lnd_rootpe, lnd_pestride, lnd_nthreads, &
         ice_ntasks, ice_rootpe, ice_pestride, ice_nthreads, &
         glc_ntasks, glc_rootpe, glc_pestride, glc_nthreads, &
         ocn_ntasks, ocn_rootpe, ocn_pestride, ocn_nthreads, &
         cpl_ntasks, cpl_rootpe, cpl_pestride, cpl_nthreads, &
         wrf_ntasks, wrf_rootpe, wrf_pestride, wrf_nthreads, &  ! by juanxiong he for wrf/cam coupling
         gea_ntasks, gea_rootpe, gea_pestride, gea_nthreads, &  ! by juanxiong he for wrf/cam coupling
		 srd_ntasks, srd_rootpe, srd_pestride, srd_nthreads  ! by huiqun hao for dc
    namelist /ccsm_pes/  &
         atm_ntasks, atm_rootpe, atm_pestride, atm_nthreads, &
         lnd_ntasks, lnd_rootpe, lnd_pestride, lnd_nthreads, &
         ice_ntasks, ice_rootpe, ice_pestride, ice_nthreads, &
         glc_ntasks, glc_rootpe, glc_pestride, glc_nthreads, &
         ocn_ntasks, ocn_rootpe, ocn_pestride, ocn_nthreads, &
         cpl_ntasks, cpl_rootpe, cpl_pestride, cpl_nthreads, &
         wrf_ntasks, wrf_rootpe, wrf_pestride, wrf_nthreads, & ! by juanxiong he for wrf/cam coupling 
         gea_ntasks, gea_rootpe, gea_pestride, gea_nthreads, & ! by juanxiong he for wrf/cam coupling 
		 srd_ntasks, srd_rootpe, srd_pestride, srd_nthreads ! by huiqun hao for dc
    !----------------------------------------------------------

    ! make sure this is first pass and set comms unset
    if (.not. first_pass) then
       write(logunit,*) trim(subname),' ERROR seq_comm_init already called '
       call shr_sys_abort()
    endif
    first_pass = .false.
    Global_Comm = Comm_in

    do n = 1,ncomps
       seq_comms(n)%set = .false.
       seq_comms(n)%mpicom = MPI_COMM_NULL    ! do some initialization here 
       seq_comms(n)%iam = -1
       seq_comms(n)%iamroot = .false.
       seq_comms(n)%npes = -1
       seq_comms(n)%nthreads = -1
       seq_comms(n)%gloiam = -1
       seq_comms(n)%gloroot = -1
       seq_comms(n)%pethreads = -1
       seq_comms(n)%cplpe = -1
       seq_comms(n)%cmppe = -1
    enddo

    ! Initialize MPI
    ! Note that if no MPI, will call MCTs fake version

    call mpi_comm_rank(GLOBAL_COMM, mype  , ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank comm_world')
    call mpi_comm_size(GLOBAL_COMM, numpes, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_size comm_world')

    ! Initialize gloiam on all IDs

    do n = 1,ncomps
       seq_comms(n)%gloiam = mype
    enddo

    ! Initialize IDs

    if (mype == 0) then

       ! Set default values
       atm_ntasks=numpes; atm_rootpe=0; atm_pestride=1; atm_nthreads=1
       lnd_ntasks=numpes; lnd_rootpe=0; lnd_pestride=1; lnd_nthreads=1
       ice_ntasks=numpes; ice_rootpe=0; ice_pestride=1; ice_nthreads=1
       glc_ntasks=numpes; glc_rootpe=0; glc_pestride=1; glc_nthreads=1
       ocn_ntasks=numpes; ocn_rootpe=0; ocn_pestride=1; ocn_nthreads=1
       cpl_ntasks=numpes; cpl_rootpe=0; cpl_pestride=1; cpl_nthreads=1
       wrf_ntasks=numpes; wrf_rootpe=0; wrf_pestride=1; wrf_nthreads=1 ! by juanxiong he for wrf/cam coupling
       gea_ntasks=numpes; gea_rootpe=0; gea_pestride=1; gea_nthreads=1 ! by juanxiong he for wrf/cam coupling
	   srd_ntasks=numpes; srd_rootpe=0; srd_pestride=1; srd_nthreads=1 ! by huiqun hao for dc

       ! Read namelist if it exists
       nu = shr_file_getUnit()
       open(nu,file=trim(nmlfile),status='old', iostat=ierr)
       if (ierr == 0) then
          ierr = 1
          do while( ierr > 0 )
             read(nu, nml=ccsm_pes, iostat=ierr)
          end do
          close(nu)
       end if
       call shr_file_freeUnit(nu)

       !--- validation of inputs ---
       ! rootpes >= 0
       if (atm_rootpe < 0 .or. lnd_rootpe < 0 .or. ice_rootpe < 0 .or. &
           ocn_rootpe < 0 .or. glc_rootpe < 0 .or. cpl_rootpe < 0 .or. &
           wrf_rootpe < 0 .or. gea_rootpe < 0 .or. srd_rootpe < 0) then    ! by juanxiong he for wrf/cam coupling
		                                           ! by huiqun hao for dc
          write(logunit,*) trim(subname),' ERROR: rootpes must be >= 0'
          call shr_sys_abort()
       endif

!       ! nthreads = 1, temporary
!       if (atm_nthreads /= 1 .or. lnd_nthreads /= 1 .or. ice_nthreads /= 1 .or. &
!           ocn_nthreads /= 1 .or. cpl_nthreads /= 1) then
!          write(logunit,*) trim(subname),' ERROR: nthreads must be 1'
!          call shr_sys_abort()
!       endif

!       ! nthreads should be 1 or something consistent, compute max nthreads
!       amax = max(atm_nthreads,lnd_nthreads)
!       amax = max(amax        ,ice_nthreads)
!       amax = max(amax        ,ocn_nthreads)
!       amax = max(amax        ,cpl_nthreads)

!       ! check that everything is either 1 or max nthreads
!       if ((atm_nthreads /= 1 .and. atm_nthreads /= amax) .or. &
!           (lnd_nthreads /= 1 .and. lnd_nthreads /= amax) .or. &
!           (ice_nthreads /= 1 .and. ice_nthreads /= amax) .or. &
!           (ocn_nthreads /= 1 .and. ocn_nthreads /= amax) .or. &
!           (cpl_nthreads /= 1 .and. cpl_nthreads /= amax)) then
!          write(logunit,*) trim(subname),' ERROR: nthreads must be consistent'
!          call shr_sys_abort()
!       endif

    endif

    ! NOTE: Only valid on root pe due to namelist read above, bcast below
!    print *, 'srd_rootpe',srd_rootpe,'srd_ntasks',srd_ntasks,'srd_pestride',srd_pestride
    if (mype == 0) then
       amin = atm_rootpe
       amax = atm_rootpe + (atm_ntasks-1)*atm_pestride
       astr = atm_pestride
       
       lmin = lnd_rootpe
       lmax = lnd_rootpe + (lnd_ntasks-1)*lnd_pestride
       lstr = lnd_pestride
       
       imin = ice_rootpe
       imax = ice_rootpe + (ice_ntasks-1)*ice_pestride
       istr = ice_pestride
       
       omin = ocn_rootpe
       omax = ocn_rootpe + (ocn_ntasks-1)*ocn_pestride
       ostr = ocn_pestride

       gmin = glc_rootpe
       gmax = glc_rootpe + (glc_ntasks-1)*glc_pestride
       gstr = glc_pestride

       cmin = cpl_rootpe
       cmax = cpl_rootpe + (cpl_ntasks-1)*cpl_pestride
       cstr = cpl_pestride

       wmin = wrf_rootpe
       wmax = wrf_rootpe + (wrf_ntasks-1)*wrf_pestride  ! by juanxiong he for wrf/cam coupling       
       wstr = wrf_pestride

       gemin = gea_rootpe
       gemax = gea_rootpe + (gea_ntasks-1)*gea_pestride  ! by juanxiong he for wrf/cam coupling       
       gestr = gea_pestride
	   
	   smin = srd_rootpe
       smax = srd_rootpe + (srd_ntasks-1)*srd_pestride  ! by huiqun hao for dc       
       sstr = srd_pestride
    end if

    ! create petlist for ESMF components
    if(present(atm_petlist)) then
        call shr_mpi_bcast(atm_ntasks, GLOBAL_COMM, 'atm_ntasks')
        call shr_mpi_bcast(atm_rootpe, GLOBAL_COMM, 'atm_rootpe')
        call shr_mpi_bcast(atm_pestride, GLOBAL_COMM, 'atm_pestride')
        allocate(atm_petlist(atm_ntasks))
        do i = 1, atm_ntasks
            atm_petlist(i) = atm_rootpe + (i-1)*atm_pestride
        enddo
    endif

    if(present(lnd_petlist)) then
        call shr_mpi_bcast(lnd_ntasks, GLOBAL_COMM, 'lnd_ntasks')
        call shr_mpi_bcast(lnd_rootpe, GLOBAL_COMM, 'lnd_rootpe')
        call shr_mpi_bcast(lnd_pestride, GLOBAL_COMM, 'lnd_pestride')
        allocate(lnd_petlist(lnd_ntasks))
        do i = 1, lnd_ntasks
            lnd_petlist(i) = lnd_rootpe + (i-1)*lnd_pestride
        enddo
    endif

    if(present(ice_petlist)) then
        call shr_mpi_bcast(ice_ntasks, GLOBAL_COMM, 'ice_ntasks')
        call shr_mpi_bcast(ice_rootpe, GLOBAL_COMM, 'ice_rootpe')
        call shr_mpi_bcast(ice_pestride, GLOBAL_COMM, 'ice_pestride')
        allocate(ice_petlist(ice_ntasks))
        do i = 1, ice_ntasks
            ice_petlist(i) = ice_rootpe + (i-1)*ice_pestride
        enddo
    endif

    if(present(ocn_petlist)) then
        call shr_mpi_bcast(ocn_ntasks, GLOBAL_COMM, 'ocn_ntasks')
        call shr_mpi_bcast(ocn_rootpe, GLOBAL_COMM, 'ocn_rootpe')
        call shr_mpi_bcast(ocn_pestride, GLOBAL_COMM, 'ocn_pestride')
        allocate(ocn_petlist(ocn_ntasks))
        do i = 1, ocn_ntasks
            ocn_petlist(i) = ocn_rootpe + (i-1)*ocn_pestride
        enddo
    endif

    if(present(glc_petlist)) then
        call shr_mpi_bcast(glc_ntasks, GLOBAL_COMM, 'glc_ntasks')
        call shr_mpi_bcast(glc_rootpe, GLOBAL_COMM, 'glc_rootpe')
        call shr_mpi_bcast(glc_pestride, GLOBAL_COMM, 'glc_pestride')
        allocate(glc_petlist(glc_ntasks))
        do i = 1, glc_ntasks
            glc_petlist(i) = glc_rootpe + (i-1)*glc_pestride
        enddo
    endif
       
    ! juanxiong he for wrf/ccsm4 coupling
    if(present(wrf_petlist)) then
        call shr_mpi_bcast(wrf_ntasks, GLOBAL_COMM, 'wrf_ntasks')
        call shr_mpi_bcast(wrf_rootpe, GLOBAL_COMM, 'wrf_rootpe')
        call shr_mpi_bcast(wrf_pestride, GLOBAL_COMM, 'wrf_pestride')
        allocate(wrf_petlist(wrf_ntasks))
        do i = 1, wrf_ntasks
            wrf_petlist(i) = wrf_rootpe + (i-1)*wrf_pestride
        enddo
    endif

    ! juanxiong he for wrf/ccsm4 coupling
    if(present(gea_petlist)) then
        call shr_mpi_bcast(gea_ntasks, GLOBAL_COMM, 'gea_ntasks')
        call shr_mpi_bcast(gea_rootpe, GLOBAL_COMM, 'gea_rootpe')
        call shr_mpi_bcast(gea_pestride, GLOBAL_COMM, 'gea_pestride')
        allocate(gea_petlist(gea_ntasks))
        do i = 1, gea_ntasks
            gea_petlist(i) = gea_rootpe + (i-1)*gea_pestride
        enddo
    endif
	
	! huiqun hao for dc
    if(present(srd_petlist)) then
        call shr_mpi_bcast(srd_ntasks, GLOBAL_COMM, 'srd_ntasks')
        call shr_mpi_bcast(srd_rootpe, GLOBAL_COMM, 'srd_rootpe')
        call shr_mpi_bcast(srd_pestride, GLOBAL_COMM, 'srd_pestride')
        allocate(srd_petlist(srd_ntasks))
        do i = 1, srd_ntasks
            srd_petlist(i) = srd_rootpe + (i-1)*srd_pestride
        enddo
    endif
    
    call shr_mpi_bcast(atm_nthreads,GLOBAL_COMM,'atm_nthreads')
    call shr_mpi_bcast(lnd_nthreads,GLOBAL_COMM,'lnd_nthreads')
    call shr_mpi_bcast(ocn_nthreads,GLOBAL_COMM,'ocn_nthreads')
    call shr_mpi_bcast(ice_nthreads,GLOBAL_COMM,'ice_nthreads')
    call shr_mpi_bcast(glc_nthreads,GLOBAL_COMM,'glc_nthreads')
    call shr_mpi_bcast(cpl_nthreads,GLOBAL_COMM,'cpl_nthreads')
    call shr_mpi_bcast(wrf_nthreads,GLOBAL_COMM,'wrf_nthreads')  ! juanxiong he for wrf/ccsm4 coupling
    call shr_mpi_bcast(gea_nthreads,GLOBAL_COMM,'gea_nthreads')  ! juanxiong he for geatm/ccsm4 coupling
	call shr_mpi_bcast(srd_nthreads,GLOBAL_COMM,'srd_nthreads')  ! huiqun hao for dc

    ! Create MPI communicator groups

    if (mype == 0) then
       pelist(1,1) = 0
       pelist(2,1) = numpes-1
       pelist(3,1) = 1
    end if
    call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
    call seq_comm_setcomm(GLOID,pelist)

    if (mype == 0) then
       pelist(1,1) = amin
       pelist(2,1) = amax
       pelist(3,1) = astr
    end if
    call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
    call seq_comm_setcomm(ATMID,pelist,atm_nthreads)

    if (mype == 0) then
       pelist(1,1) = lmin
       pelist(2,1) = lmax
       pelist(3,1) = lstr
    end if
    call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
    call seq_comm_setcomm(LNDID,pelist,lnd_nthreads)

    if (mype == 0) then
       pelist(1,1) = imin
       pelist(2,1) = imax
       pelist(3,1) = istr
    end if
    call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
    call seq_comm_setcomm(ICEID,pelist,ice_nthreads)

    if (mype == 0) then
       pelist(1,1) = gmin
       pelist(2,1) = gmax
       pelist(3,1) = gstr
    end if
    call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
    call seq_comm_setcomm(GLCID,pelist,glc_nthreads)

    if (mype == 0) then
       pelist(1,1) = omin
       pelist(2,1) = omax
       pelist(3,1) = ostr
    end if
    call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
    call seq_comm_setcomm(OCNID,pelist,ocn_nthreads)

    if (mype == 0) then
       pelist(1,1) = cmin
       pelist(2,1) = cmax
       pelist(3,1) = cstr
    end if
    call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
    call seq_comm_setcomm(CPLID,pelist,cpl_nthreads)

    if (mype == 0) then             ! by juanxiong he for wrf/cam coupling
       pelist(1,1) = wmin            ! by juanxiong he for wrf/cam coupling 
       pelist(2,1) = wmax            ! by juanxiong he for wrf/cam coupling
       pelist(3,1) = wstr            ! by juanxiong he for wrf/cam coupling
    end if                          ! by juanxiong he for wrf/cam coupling 
    call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr) ! by juanxiong he for wrf/cam coupling
    call seq_comm_setcomm(WRFID,pelist,wrf_nthreads) ! by juanxiong he for wrf/cam coupling   

    if (mype == 0) then             ! by juanxiong he for geatm/cam coupling
       pelist(1,1) = gemin            ! by juanxiong he for geatm/cam coupling 
       pelist(2,1) = gemax            ! by juanxiong he for geatm/cam coupling
       pelist(3,1) = gestr            ! by juanxiong he for geatm/cam coupling
    end if                          ! by juanxiong he for geatm/cam coupling 
	
!	if (mype == 0) write(221,*) 'pelist',pelist
	
    call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr) ! by juanxiong he for geatm/cam coupling
    call seq_comm_setcomm(GEAID,pelist,gea_nthreads) ! by juanxiong he for geatm/cam coupling   
	
!    if (mype == 0) write(211,*) 'gea_ntasks',gea_ntasks,'seq_comms(GEAID)%mpicom',seq_comms(GEAID)%mpicom
	
	if (mype == 0) then             ! by huiqun hao for dc
       pelist(1,1) = smin            ! by huiqun hao for dc 
       pelist(2,1) = smax            ! by huiqun hao for dc
       pelist(3,1) = sstr            ! by huiqun hao for dc
    end if                          ! by huiqun hao for dc 
	
!	if (mype == 0) write(222,*) 'pelist',pelist
!	print

    call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr) ! by huiqun hao for dc
!    write(logunit,*) 'srd pelist',pelist
    call seq_comm_setcomm(SRDID,pelist,srd_nthreads) ! by huiqun hao for dc   
	
!    if (mype == 0) write(212,*) 'srd_ntasks',srd_ntasks,'seq_comms(SRDID)%mpicom',seq_comms(SRDID)%mpicom
	
    call seq_comm_joincomm(CPLID,ATMID,CPLATMID)
    call seq_comm_joincomm(CPLID,LNDID,CPLLNDID)
    call seq_comm_joincomm(CPLID,ICEID,CPLICEID)
    call seq_comm_joincomm(CPLID,OCNID,CPLOCNID)
    call seq_comm_joincomm(CPLID,GLCID,CPLGLCID)
    call seq_comm_joincomm(CPLID,WRFID,CPLWRFID) ! by juanxiong he for wrf/cam coupling
    call seq_comm_joincomm(CPLID,GEAID,CPLGEAID) ! by juanxiong he for geatm/cam coupling
	call seq_comm_joincomm(SRDID,ATMID,SRDATMID) ! by huiqun hao for dc
	call seq_comm_joincomm(SRDID,WRFID,SRDWRFID) ! by huiqun hao for dc

    max_threads = -1
    do n = 1,ncomps
       max_threads = max(max_threads,seq_comms(n)%nthreads)
    enddo
    do n = 1,ncomps
       seq_comms(n)%pethreads = max_threads
    enddo

    ! compute each components root pe global id and broadcast so all pes have info
    do n = 1,ncomps
       gloroot = -999
       if (seq_comms(n)%iamroot) gloroot = seq_comms(n)%gloiam
       call shr_mpi_max(gloroot,seq_comms(n)%gloroot,GLOBAL_COMM, &
                        trim(subname)//' gloroot',all=.true.)
    enddo

    ! Initialize MCT

    ! add up valid comps on local pe

    myncomps = 0
    do n = 1,ncomps
       if (seq_comms(n)%mpicom /= MPI_COMM_NULL) then
          myncomps = myncomps + 1
       endif
    enddo

    ! set comps and comms

    allocate(comps(myncomps),comms(myncomps),stat=ierr)
    if(ierr/=0) call mct_die(subName,'allocate comps comms',ierr)

    myncomps = 0
    do n = 1,ncomps
       if (seq_comms(n)%mpicom /= MPI_COMM_NULL) then
          myncomps = myncomps + 1
          if (myncomps > size(comps)) then
             write(logunit,*) trim(subname),' ERROR in myncomps ',myncomps,size(comps)
             call shr_sys_abort()
          endif
          comps(myncomps) = seq_comms(n)%ID
          comms(myncomps) = seq_comms(n)%mpicom
          onecomm = seq_comms(n)%mpicom   ! if one unique comm per pe, then pick any
       endif
    enddo

    if (myncomps /= size(comps)) then
       write(logunit,*) trim(subname),' ERROR in myncomps ',myncomps,size(comps)
       call shr_sys_abort()
    endif

    call mct_world_init(ncomps, GLOBAL_COMM, comms, comps)

    deallocate(comps,comms)

    call seq_comm_printcomms()

  end subroutine seq_comm_init

!---------------------------------------------------------
  subroutine seq_comm_setcomm(ID,pelist,nthreads)

    implicit none
    integer,intent(IN) :: ID
    integer,intent(IN) :: pelist(:,:)
    integer,intent(IN),optional :: nthreads

    integer :: mpigrp_world
    integer :: mpigrp
    integer :: mpicom
    integer :: ierr
    character(*),parameter :: subName =   '(seq_comm_setcomm) '

    integer :: ntasks_x,ntasks_y,ntasks ! added by juanxiong he 
    integer :: gntasks_x,gntasks_y,gntasks ! added by juanxiong he 
    INTEGER :: M, N, MINI, mytask, mytask_x, mytask_y,comdup !added by juanxiong he
    INTEGER, DIMENSION(2) :: dims, coords  ! juanxiong he
    LOGICAL, DIMENSION(2) :: isperiodic  !juanxiong he


    if (ID < 1 .or. ID > ncomps) then
       write(logunit,*) subname,' ID out of range, abort ',ID
       call shr_sys_abort()
    endif 

    call mpi_comm_group(GLOBAL_COMM, mpigrp_world, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_group mpigrp_world')
    call mpi_group_range_incl(mpigrp_world, 1, pelist, mpigrp,ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_group_range_incl mpigrp')
    call mpi_comm_create(GLOBAL_COMM, mpigrp, mpicom, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_create mpigrp')

    seq_comms(ID)%set = .true.
    seq_comms(ID)%ID = ID
    seq_comms(ID)%name = IDname(ID)
    seq_comms(ID)%mpicom = mpicom
    seq_comms(ID)%mpigrp = mpigrp
    if (present(nthreads)) then
       seq_comms(ID)%nthreads = nthreads
    else
       seq_comms(ID)%nthreads = 1
    endif

    if (mpicom /= MPI_COMM_NULL) then
       call mpi_comm_size(mpicom,seq_comms(ID)%npes,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_size')
       call mpi_comm_rank(mpicom,seq_comms(ID)%iam,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank')
       if (seq_comms(ID)%iam == 0) then
          seq_comms(ID)%iamroot = .true.
       else
          seq_comms(ID)%iamroot = .false.
       endif
    else
       seq_comms(ID)%npes = -1
       seq_comms(ID)%iam = -1
       seq_comms(ID)%nthreads = 1
       seq_comms(ID)%iamroot = .false.
    endif

! ---------------------------------------------------------------------
! added by juanxiong he, begin
! ---------------------------------------------------------------------
    seq_comms(ID)%mpicomcart_periodic = 0
    seq_comms(ID)%mpicomcart = 0
    seq_comms(ID)%ntasks_x = 0
    seq_comms(ID)%ntasks_y = 0
    seq_comms(ID)%mpicomx = 0
    seq_comms(ID)%mpicomy = 0

    if(ID.eq.WRFID .and. seq_comms(ID)%iam >= 0) then
      ntasks=seq_comms(ID)%npes
      MINI = 2*ntasks
      ntasks_x = 1
      ntasks_y = ntasks
      DO M = 1, ntasks
        IF ( MOD( ntasks, M ) .EQ. 0 ) THEN
          N = ntasks / M
          IF ( ABS(M-N) .LT. MINI                &
               .AND. M .GE. 1            &
               .AND. N .GE. 1            &
             ) THEN
            MINI = ABS(M-N)
            ntasks_x = M
            ntasks_y = N
!            ntasks_x = 32
!            ntasks_y = 16
          ENDIF
        ENDIF
      ENDDO


     if (mpicom /= MPI_COMM_NULL) then
      dims(1) = ntasks_y  ! rows
      dims(2) = ntasks_x  ! columns
      isperiodic(1) = .true.
      isperiodic(2) = .true.
      CALL mpi_cart_create( mpicom, 2, dims, isperiodic, .false., seq_comms(ID)%mpicomcart_periodic, ierr )

      dims(1) = ntasks_y  ! rows
      dims(2) = ntasks_x  ! columns
      isperiodic(1) = .false.
      isperiodic(2) = .false.
      CALL mpi_cart_create( seq_comms(ID)%mpicom, 2, dims, isperiodic, .false., seq_comms(ID)%mpicomcart, ierr )
      CALL mpi_comm_rank( seq_comms(ID)%mpicomcart, mytask, ierr )
      CALL mpi_cart_coords(seq_comms(ID)%mpicomcart, mytask, 2, coords, ierr )

      mytask_x = coords(2)   ! col task (x)
      mytask_y = coords(1)   ! row task (y)

      CALL MPI_Comm_dup( seq_comms(ID)%mpicomcart, comdup, ierr )
      CALL MPI_Comm_split(comdup,mytask_y,mytask,seq_comms(ID)%mpicomx,ierr)
      CALL MPI_Comm_split(comdup,mytask_x,mytask,seq_comms(ID)%mpicomy,ierr)

       seq_comms(ID)%ntasks_x=ntasks_x
       seq_comms(ID)%ntasks_y=ntasks_y
     endif

    endif

! ---------------------------------------------------------------------
! added by juanxiong he, end
! ---------------------------------------------------------------------

    if (seq_comms(ID)%iamroot) then
       write(logunit,F11) trim(subname),'  initialize ID ',ID,seq_comms(ID)%name, &
         ' pelist   =',pelist,' npes =',seq_comms(ID)%npes,' nthreads =',seq_comms(ID)%nthreads
    endif

  end subroutine seq_comm_setcomm

!---------------------------------------------------------
  subroutine seq_comm_joincomm(ID1,ID2,ID)

    implicit none
    integer,intent(IN) :: ID1    ! src id
    integer,intent(IN) :: ID2    ! srd id
    integer,intent(IN) :: ID     ! computed id

    integer :: mpigrp
    integer :: mpicom
    integer :: ierr
    integer :: n,nsize
    integer,allocatable :: pe_t1(:),pe_t2(:)
    character(*),parameter :: subName =   '(seq_comm_joincomm) '

    ! check that IDs are in valid range, that ID1 and ID2 have
    ! been set, and that ID has not been set

    if (ID1 < 1 .or. ID1 > ncomps) then
       write(logunit,*) subname,' ID1 out of range, abort ',ID1
       call shr_sys_abort()
    endif 
    if (ID2 < 1 .or. ID2 > ncomps) then
       write(logunit,*) subname,' ID2 out of range, abort ',ID2
       call shr_sys_abort()
    endif 
    if (ID < 1 .or. ID > ncomps) then
       write(logunit,*) subname,' ID out of range, abort ',ID
       call shr_sys_abort()
    endif
    if (.not. seq_comms(ID1)%set .or. .not. seq_comms(ID2)%set) then
       write(logunit,*) subname,' ID1 or ID2 not set ',ID1,ID2
       call shr_sys_abort()
    endif
    if (seq_comms(ID)%set) then
       write(logunit,*) subname,' ID already set ',ID
       call shr_sys_abort()
    endif

    call mpi_group_union(seq_comms(ID1)%mpigrp,seq_comms(ID2)%mpigrp,mpigrp,ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_union mpigrp')
    call mpi_comm_create(GLOBAL_COMM, mpigrp, mpicom, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_create mpigrp')

    seq_comms(ID)%set = .true.
    seq_comms(ID)%ID = ID
    seq_comms(ID)%name = IDname(ID)
    seq_comms(ID)%mpicom = mpicom
    seq_comms(ID)%mpigrp = mpigrp
    seq_comms(ID)%nthreads = max(seq_comms(ID1)%nthreads,seq_comms(ID2)%nthreads)
    seq_comms(ID)%nthreads = max(seq_comms(ID)%nthreads,1)

    if (mpicom /= MPI_COMM_NULL) then
       call mpi_comm_size(mpicom,seq_comms(ID)%npes,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_size')
       call mpi_comm_rank(mpicom,seq_comms(ID)%iam,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank')
       if (seq_comms(ID)%iam == 0) then
          seq_comms(ID)%iamroot = .true.
       else
          seq_comms(ID)%iamroot = .false.
       endif
    else
       seq_comms(ID)%npes = -1
       seq_comms(ID)%iam = -1
       seq_comms(ID)%iamroot = .false.
    endif

! needs to be excluded until mpi_group_size is added to serial mpi in mct
#if (1 == 0)
    if (loglevel > 3) then
       ! some debug code to prove the join is working ok
       ! when joining mpicomms, the local rank may be quite different
       !   from either the global or local ranks of the joining comms
       call mpi_group_size(seq_comms(ID1)%mpigrp,nsize,ierr)
       allocate(pe_t1(nsize),pe_t2(nsize))
       do n = 1,nsize
          pe_t1(n) = n-1
          pe_t2(n) = -1
       enddo
       call mpi_group_translate_ranks(seq_comms(ID1)%mpigrp, nsize, pe_t1, mpigrp, pe_t2, ierr)
       write(logunit,*) 'ID1      ranks ',pe_t1
       write(logunit,*) 'ID1-JOIN ranks ',pe_t2
       deallocate(pe_t1,pe_t2)

       call mpi_group_size(seq_comms(ID2)%mpigrp,nsize,ierr)
       allocate(pe_t1(nsize),pe_t2(nsize))
       do n = 1,nsize
          pe_t1(n) = n-1
          pe_t2(n) = -1
       enddo
       call mpi_group_translate_ranks(seq_comms(ID2)%mpigrp, nsize, pe_t1, mpigrp, pe_t2, ierr)
       write(logunit,*) 'ID2      ranks ',pe_t1
       write(logunit,*) 'ID2-JOIN ranks ',pe_t2
       deallocate(pe_t1,pe_t2)
    endif
#endif

    allocate(pe_t1(1),pe_t2(1))
    pe_t1(1) = 0
    call mpi_group_translate_ranks(seq_comms(ID1)%mpigrp, 1, pe_t1, mpigrp, pe_t2, ierr)
    seq_comms(ID)%cplpe = pe_t2(1)
    pe_t1(1) = 0
    call mpi_group_translate_ranks(seq_comms(ID2)%mpigrp, 1, pe_t1, mpigrp, pe_t2, ierr)
    seq_comms(ID)%cmppe = pe_t2(1)
    deallocate(pe_t1,pe_t2)

    if (seq_comms(ID)%iamroot) then
       if (loglevel > 1) then
          write(logunit,F12) trim(subname),' initialize ID ',ID,seq_comms(ID)%name, &
          ' join IDs =',ID1,ID2,' npes =',seq_comms(ID)%npes, &
          ' nthreads =',seq_comms(ID)%nthreads, &
          ' cpl/cmp pes =',seq_comms(ID)%cplpe,seq_comms(ID)%cmppe
       else
          write(logunit,F13) trim(subname),' initialize ID ',ID,seq_comms(ID)%name, &
          ' join IDs =',ID1,ID2,' npes =',seq_comms(ID)%npes, &
          ' nthreads =',seq_comms(ID)%nthreads
       endif
    endif

  end subroutine seq_comm_joincomm

!---------------------------------------------------------
  subroutine seq_comm_printcomms()

    implicit none
    character(*),parameter :: subName =   '(seq_comm_printcomms) '
    integer :: m,n,mype,npes,cpes,ierr
    character(len=256) :: iamstring
    character(*),parameter :: F01 = "(4x,a4,4x   ,40(1x,a8))"
    character(*),parameter :: F02 = "(4x,i4,3x,a1,40(2x,i6,1x))"
    character(*),parameter :: F03 = "(4x,i4,3x,a1,a)"

    call mpi_comm_size(GLOBAL_COMM, npes  , ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_size comm_world')
    call mpi_comm_rank(GLOBAL_COMM, mype  , ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank comm_world')

    call shr_sys_flush(logunit)
    call mpi_barrier(GLOBAL_COMM,ierr)
    if (mype == 0) then
!       do n = 1,ncomps
!          call mpi_group_size(seq_comms(n)%mpigrp, cpes, ierr)
!          call shr_mpi_chkerr(ierr,subname//' mpi_group_size')
!          write(logunit,*) trim(subName),' comp ntasks,nthreads ',n,trim(seq_comms(n)%name), &
!             seq_comms(n)%npes,seq_comms(n)%nthreads
!       enddo
       write(logunit,*) ' '
       write(logunit,*) trim(subName),' ID layout : global pes vs local pe for each ID'
       write(logunit,F01) ' gpe',(seq_comms(n)%name,n=1,ncomps),'nthrds'
       write(logunit,F01) ' ---',(' ------ '       ,n=1,ncomps),'------'
       call shr_sys_flush(logunit)
    endif
    iamstring = ' '
    do n = 1,ncomps
       if (seq_comms(n)%iam >= 0) then
          write(iamstring((n-1)*9+1:n*9),"(2x,i6,1x)") seq_comms(n)%iam
       endif
    enddo
    n = ncomps + 1
    write(iamstring((n-1)*9+1:n*9),"(2x,i6,1x)") seq_comms(GLOID)%pethreads

    call shr_sys_flush(logunit)
    call mpi_barrier(GLOBAL_COMM,ierr)
    do m = 0,npes-1
       if (mype == m) then
!          write(logunit,F02) mype,':',(seq_comms(n)%iam,n=1,ncomps)
          write(logunit,F03) mype,':',trim(iamstring)
          if (m == npes-1) then
             write(logunit,*) ' '
          endif
       endif
       call shr_sys_flush(logunit)
       call mpi_barrier(GLOBAL_COMM,ierr)
    enddo

  end subroutine seq_comm_printcomms

!---------------------------------------------------------
  subroutine seq_comm_setptrs(ID,mpicom,mpigrp,npes,nthreads,iam,iamroot,gloiam,gloroot,cplpe,cmppe,pethreads,&
                              mpicomcart,mpicomcart_periodic,mpicomx,mpicomy,ntasks_x,ntasks_y,&
                              gmpicomcart, gmpicomcart_periodic, gmpicomx, gmpicomy, gntasks_x, gntasks_y)  ! juanxiong he 

    implicit none
    integer,intent(in) :: ID
    integer,intent(out),optional :: mpicom
    integer,intent(out),optional :: mpigrp
    integer,intent(out),optional :: npes
    integer,intent(out),optional :: nthreads
    integer,intent(out),optional :: iam
    logical,intent(out),optional :: iamroot
    integer,intent(out),optional :: gloiam
    integer,intent(out),optional :: gloroot
    integer,intent(out),optional :: cplpe
    integer,intent(out),optional :: cmppe
    integer,intent(out),optional :: pethreads
    integer,intent(out),optional :: mpicomcart ! juanxiong he
    integer,intent(out),optional :: mpicomcart_periodic ! juanxiong he
    integer,intent(out),optional :: mpicomx  ! juanxiong he 
    integer,intent(out),optional :: mpicomy  ! juanxiong he 
    integer,intent(out),optional :: ntasks_x  ! juanxiong he 
    integer,intent(out),optional :: ntasks_y  ! juanxiong he 
    integer,intent(out),optional :: gmpicomcart ! juanxiong he
    integer,intent(out),optional :: gmpicomcart_periodic ! juanxiong he
    integer,intent(out),optional :: gmpicomx  ! juanxiong he 
    integer,intent(out),optional :: gmpicomy  ! juanxiong he 
    integer,intent(out),optional :: gntasks_x  ! juanxiong he 
    integer,intent(out),optional :: gntasks_y  ! juanxiong he 

    character(*),parameter :: subName =   '(seq_comm_setptrs) '

    if (ID < 1 .or. ID > ncomps) then
       write(logunit,*) subname,' ID out of range, return ',ID,seq_comms(ID)%name
       return
    endif 

    if (present(mpicom)) then
       mpicom = seq_comms(ID)%mpicom
    endif

    if (present(mpigrp)) then
       mpigrp = seq_comms(ID)%mpigrp
    endif

    if (present(npes)) then
       npes = seq_comms(ID)%npes
    endif

    if (present(nthreads)) then
       nthreads = seq_comms(ID)%nthreads
    endif

    if (present(iam)) then
       iam = seq_comms(ID)%iam
    endif

    if (present(iamroot)) then
       iamroot = seq_comms(ID)%iamroot
    endif

    if (present(gloiam)) then
       gloiam = seq_comms(ID)%gloiam
    endif

    if (present(gloroot)) then
       gloroot = seq_comms(ID)%gloroot
    endif

    if (present(cplpe)) then
       cplpe = seq_comms(ID)%cplpe
    endif

    if (present(cmppe)) then
       cmppe = seq_comms(ID)%cmppe
    endif

    if (present(pethreads)) then
       pethreads = seq_comms(ID)%pethreads
    endif

    if (present(mpicomcart)) then   !juanxiong he
       mpicomcart = seq_comms(ID)%mpicomcart
    endif

    if (present(mpicomcart_periodic)) then !juanxiong he
       mpicomcart_periodic = seq_comms(ID)%mpicomcart_periodic
    endif

    if (present(mpicomx)) then !juanxiong he
       mpicomx = seq_comms(ID)%mpicomx
    endif

    if (present(mpicomy)) then !juanxiong he
       mpicomy = seq_comms(ID)%mpicomy
    endif

    if (present(ntasks_x)) then !juanxiong he
       ntasks_x = seq_comms(ID)%ntasks_x
    endif

    if (present(ntasks_y)) then !juanxiong he
       ntasks_y = seq_comms(ID)%ntasks_y
    endif

    if (present(gmpicomcart)) then   !juanxiong he
       gmpicomcart = seq_comms(ID)%gmpicomcart
    endif

    if (present(gmpicomcart_periodic)) then !juanxiong he
       gmpicomcart_periodic = seq_comms(ID)%gmpicomcart_periodic
    endif

    if (present(gmpicomx)) then !juanxiong he
       gmpicomx = seq_comms(ID)%gmpicomx
    endif

    if (present(gmpicomy)) then !juanxiong he
       gmpicomy = seq_comms(ID)%gmpicomy
    endif

    if (present(gntasks_x)) then !juanxiong he
       gntasks_x = seq_comms(ID)%gntasks_x
    endif

    if (present(gntasks_y)) then !juanxiong he
       gntasks_y = seq_comms(ID)%gntasks_y
    endif

  end subroutine seq_comm_setptrs
!---------------------------------------------------------
  subroutine seq_comm_setnthreads(nthreads)

    implicit none
    integer,intent(in) :: nthreads
    character(*),parameter :: subName =   '(seq_comm_setnthreads) '

#ifdef _OPENMP
    if (nthreads < 1) then
       call shr_sys_abort(subname//' ERROR: nthreads less than one')
    endif
    call omp_set_num_threads(nthreads)
#endif

  end subroutine seq_comm_setnthreads
!---------------------------------------------------------
  integer function seq_comm_getnthreads()

    implicit none
    integer :: omp_get_num_threads
    character(*),parameter :: subName =   '(seq_comm_getnthreads) '

    seq_comm_getnthreads = -1
#ifdef _OPENMP
!$OMP PARALLEL
    seq_comm_getnthreads = omp_get_num_threads()
!$OMP END PARALLEL
#endif

  end function seq_comm_getnthreads
!---------------------------------------------------------
  logical function seq_comm_iamin(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_iamin) '

    if (seq_comms(ID)%iam >= 0) then
       seq_comm_iamin = .true.
    else
       seq_comm_iamin = .false.
    endif

  end function seq_comm_iamin
!---------------------------------------------------------
  logical function seq_comm_iamroot(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_iamroot) '

    seq_comm_iamroot = seq_comms(ID)%iamroot

  end function seq_comm_iamroot
!---------------------------------------------------------
  integer function seq_comm_mpicom(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_mpicom) '

    seq_comm_mpicom = seq_comms(ID)%mpicom

  end function seq_comm_mpicom
!---------------------------------------------------------
  integer function seq_comm_iam(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_iam) '

    seq_comm_iam = seq_comms(ID)%iam

  end function seq_comm_iam
!---------------------------------------------------------
  integer function seq_comm_gloiam(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_gloiam) '

    seq_comm_gloiam = seq_comms(ID)%gloiam

  end function seq_comm_gloiam
!---------------------------------------------------------
  integer function seq_comm_gloroot(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_gloroot) '

    seq_comm_gloroot = seq_comms(ID)%gloroot

  end function seq_comm_gloroot
!---------------------------------------------------------
  integer function seq_comm_cplpe(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_cplpe) '

    seq_comm_cplpe = seq_comms(ID)%cplpe

  end function seq_comm_cplpe
!---------------------------------------------------------
  integer function seq_comm_cmppe(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_cmppe) '

    seq_comm_cmppe = seq_comms(ID)%cmppe

  end function seq_comm_cmppe
!---------------------------------------------------------
end module seq_comm_mct
