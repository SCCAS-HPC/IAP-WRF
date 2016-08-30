!----------------------------------------------------------------------------
!    This program bases on the original GEATM.  It refactorys the original
!    GEATM with the new modulization and stratification. The new top layer
!    program is added to embed GEATM wihtin CESM. The output method is 
!    is replaced with the new method to make sure the efficiency.
!     
!    Juanxiong He, 2013-06-30
!----------------------------------------------------------------------------
module geatm_comp_mct
  use mct_mod
  use esmf_mod
  use seq_flds_mod
  use seq_flds_indices
  use seq_cdata_mod
  use seq_infodata_mod
  use seq_timemgr_mod
  use seq_comm_mct, only:seq_comm_iamroot,seq_comm_iamin
  use shr_kind_mod     , only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use shr_file_mod     , only: shr_file_getunit, shr_file_freeunit, &
                               shr_file_setLogUnit, shr_file_setLogLevel, &
                               shr_file_getLogUnit, shr_file_getLogLevel, &
                               shr_file_setIO
  use shr_sys_mod      , only: shr_sys_flush, shr_sys_abort
  
  use parall
  use inout
  use gaschem
  use aqueous_chem
  use getclouddepth_c
  use scavrat_c
  use wetdepaer
  use wetdepgas
  use hgchem
  use adv_mark
  use chemprod_mark
  use diffusion_mark
  use getsmark
  use inout2
  use inout3
  use tropause
  use termbal_model
  use setwindboundary
  use setboundary
  use setboundtop_model
  use nbalance
  use modis
  use model
  use getdu_model
  use getboundary
  use eddyz_model
  use drydep
  use diffusion
  use convect43c
  use alld4main
  use advection

  public :: geatm_init_mct
  public :: geatm_run_mct
  public :: geatm_final_mct

  private :: geatm_SetgsMap_cam
  private :: geatm_import_mct
  private :: geatm_export_mct
  private :: geatm_domain_cam
 
  integer, parameter  :: nlen = 256     ! Length of character strings

!
! Time averaged counter for flux fields
!
  integer :: avg_count
  integer, parameter :: cam_to_geatm_points = 1000000
  integer,dimension(1:cam_to_geatm_points) :: remap_cami, remap_camj
  real*8,dimension(1:cam_to_geatm_points) :: cam_weight
  integer,dimension(:,:),allocatable ::cam_local_points 

  type wrfgrid_c
  integer :: num_wrfgrid_soil_levels, num_wrfgrid_levels, num_wrfgrid_lat, num_wrfgrid_lon
  real(8),dimension(:),allocatable::sigma
  real(8),dimension(:,:,:),allocatable::u3d,v3d,t3d,qv3d,qi3d,qc3d,rh3d,p3d,z3d,taucldv3d,taucldi3d,&
                                     soilt,soilm,soildepth,soilthick
  real(8),dimension(:,:),allocatable::ps,ht,swdown,raincv,rainncv,clflo,clfmi,clfhi,&
                                   pblh,rmol,ust,u10,v10,t2,q2,rh2,&
                                   tsk,xland, sst, xice, snowh
  end type wrfgrid_c

  type camgrid_c
  integer :: num_camgrid_soil_levels, num_camgrid_levels, num_camgrid_lat, num_camgrid_lon
  real(8),dimension(:,:,:),allocatable::u3d,v3d,t3d,qv3d,qi3d,qc3d,rh3d,p3d,z3d,taucldv3d,taucldi3d,&
                                     soilt,soilm,soildepth,soilthick
  real(8),dimension(:,:),allocatable::xlat,xlon,ps,ht,swdown,raincv,rainncv,clflo,clfmi,clfhi,&
                                   pblh,rmol,ust,u10,v10,t2,q2,rh2,&
                                   tsk,xland, sst, xice, snowh
  end type camgrid_c
  type(camgrid_c):: camgrid
  type(wrfgrid_c):: wrfgrid  

contains

subroutine geatm_init_mct(EClock, cdata_ge, x2c_c1, x2c_c2, c2x_c)
  use geatm_vartype
 
    include 'mpif.h'
 
    type(ESMF_Clock),intent(in)                 :: EClock
    type(seq_cdata), intent(inout)              :: cdata_ge
    type(mct_aVect), intent(inout)              :: x2c_c1
    type(mct_aVect), intent(inout)              :: x2c_c2
    type(mct_aVect), intent(inout)              :: c2x_c   

    type(mct_gsMap), pointer   :: gsMap_ge
    type(mct_gGrid), pointer   :: dom_ge
    type(seq_infodata_type),pointer :: infodata
   
    integer :: GEATMID    
    integer :: mpicom_atm
    integer :: lsize    
    logical :: exists           ! true if file exists    
    integer :: stepno           ! time step                      
    integer :: dtime_sync       ! integer timestep size
    integer :: currentymd       ! current year-month-day
    integer :: dtime            ! time step increment (sec)
    integer :: geatm_cpl_dt       ! driver geatm coupling time step 
    integer :: dtime_geatm        ! Time-step increment (sec)
    integer :: ymd              ! CAM current date (YYYYMMDD)
    integer :: yr               ! CAM current year
    integer :: mon              ! CAM current month
    integer :: day              ! CAM current day
    integer :: tod              ! CAM current time of day (sec)
    integer :: start_ymd        ! Start date (YYYYMMDD)
    integer :: start_tod        ! Start time of day (sec)
    integer :: ref_ymd          ! Reference date (YYYYMMDD)
    integer :: ref_tod          ! Reference time of day (sec)
    integer :: stop_ymd         ! Stop date (YYYYMMDD)
    integer :: stop_tod         ! Stop time of day (sec)
    integer :: shrlogunit,shrloglev ! old values
    logical,save :: first_time = .true.
    character(len=SHR_KIND_CS) :: calendar  ! Calendar type
    character(len=SHR_KIND_CS) :: starttype ! infodata start type
    integer :: lbnum
    integer :: i, j, k, ne, iiaer, ig, ierr, iisize
    integer :: iulog
    integer,dimension(:),allocatable :: cpuid
    integer,dimension(2) ::     buffer
    integer status(MPI_STATUS_SIZE)    
    include 'chm1.inc'
    include 'gas1.inc'
    include 'params1'
    include 'tuv.inc'

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)

       call seq_cdata_setptrs(cdata_ge, ID=GEATMID, mpicom=mpicom_atm,&
                              ntasks=numprocs,gsMap=gsMap_ge, dom=dom_ge, infodata=infodata)

if(first_time) then
	    
       if (seq_comm_iamroot(GEATMID)) then
          inquire(file='gea_modelio.nml',exist=exists)
          if (exists) then
             iulog = shr_file_getUnit()
             call shr_file_setIO('gea_modelio.nml',iulog)
          endif
          write(iulog,*) "GEATM initialization"
       endif
       
       print *,'***mpi_communicator_geatm***',mpicom_atm

!--------------------------------------------------------------------
!   parameter reading
!--------------------------------------------------------------------       
    OPEN(31,FILE='input.dat',form='formatted')
      
    read(31,*)iyear1,imonth1,idate1,ihour1
    read(31,*) camgrid%num_camgrid_lon, camgrid%num_camgrid_lat, &
               camgrid%num_camgrid_levels, camgrid%num_camgrid_soil_levels

    read(31,*)nest
    if(nest.gt.ndomain)then
     print *,'the maxinum domains are set as five(5),please revise it'
     stop 
    endif

    read(31,*)ntt
    do i=1,nest
     read(31,*)ntbeg(i),ntend(i)
    enddo
    read(31,*) (dtstep(ne),ne=1,nest) !added by chenhs, time step in seconds
    read(31,*) (dtmet(ne),ne=1,nest) !added by chenhs,meteorology input frequency
    read(31,*) (dtout(ne),ne=1,nest) !added by chenhs,write out frequency in hours
    read(31,*)idifvert      ! to select the diffusion shceme 1: local 2: acm2
    read(31,*)ichemgas      ! to select gas chem
    read(31,*)idry      ! to select drydepostion scheme 1: constant 2: calculated
    read(31,*)iglobal       ! to select global conditions
    read(31,*)ifseacom,ifdustcom ! sea-salt composition,1:yes,0:no; dust composition,1:yes,0:no 
    read(31,*)imodis        ! to set modis landuse(1);0 is USGS landuse
    read(31,*)ifprocess        ! to decide if do the process analysis, by chenhs
    read(31,*)ifglobal        ! if use online global domain, by chenhs
    read(31,*)ifHg          ! if conduct Hg simulation, by chenhs
    read(31,*)press         ! to read the height(hpa) of NAQPMS top conditions from global
    read(31,*)hh            ! the top height of naqpms
    read(31,*)nzz
    do ne=1,nest
      read(31,*)nx(ne),ny(ne),nz(ne),nxlo(ne),nylo(ne)
    enddo
    read(31,*) (sratio(ne),ne=1,nest)  ! to determine the dx and dy resolution ratio between mother and child domian
    read(31,*) (imother(ne),ne=1,nest) !added by chenhs, to read info of the mother domain 
    read(31,*) isize        ! to determine the size of aerosols
    read(31,*) iaer         ! to determine the type of aerosols

    allocate(gravel(isize,iaer))

    do iiaer=1,iaer
      read(31,*)(gravel(iisize,iiaer),iisize=1,isize)
    enddo

    ALLOCATE (MSIZDIS(ISIZE+1),MSIZDID(ISIZE+1)) ! SEA SALT AND DUST PARTICLES SIZE DISTRIBUTION
    READ(31,*) (MSIZDID(IISIZE),IISIZE = 1, ISIZE+1 ) ! DUST
    READ(31,*) (MSIZDIS(IISIZE),IISIZE = 1, ISIZE+1 ) ! SEA
 
    KDUSTTOP = 6 ! DUST CAN TRANSPORTED TO 1KM TOP

    ! For Source Mark
    read(31,*)(ifsm(ne),ne=1,nest)    ! if we need Source Mark ?
    ifsmt=0
    do ne=1,nest
      ifsmt=ifsm(ne)+ifsmt
    enddo
    if(ifsmt>0)then
    read(31,*)idmSet,iSrcDefined,ismMax,iHgtLMax
              ! idmSet to define how many species
              ! iSrcDefined to define how many locations
              ! how many sourceis to mark; it should be ismMax=3+IsrcDefined*iHgtLMax
              ! iHgtLMax to define how many types of emissions
    ismMaxHere=3+IsrcDefined*iHgtLMax
    if(ismMax.ne.ismMaxHere)then
        print *, 'ismMax should be equel to 3+IsrcDefined*iHgtLMax'
        print *, ismMax,ismMaxHere
        stop
    endif

    allocate(igMark(idmSet))
    allocate(iaMarkAer(idmSet))
    allocate(iaMarkSiz(idmSet)) 
    ! to define how many types of emissions
    allocate(TmpSM(ismMax))
    allocate(contem0(idmSet))
    do idm=1,idmSet
       read(31,*)igMark(idm),iaMarkAer(idm),iaMarkSiz(idm) 
    enddo
    ! to define how many locations
    ! how many sourceis to mark
     endif
     close(31)
     
     if(ifHg.eq.1) then 
        igas = 105
     else
        igas = 102
     endif
     igasCBM = 74 
     iopen = 1
     iprecise = 1  ! 1 best precise, 2 the second
     ikosaline = 1 ! 1 on-line       2 off-line
     imasskeep = 0 ! 1 to keep mass-keeping, chenhs test
     iprocess=24   ! number of process analysis 
     time = 0
     NSOA = 6      ! TOTAL NUMBERS OF SOA
     ndustcom =  11 ! dust aerosols compositions numbers
     nseacom = 8  ! sea salt compositions numbers
     ifbalance = 1 ! 1 to keep mass-keeping, chenhs set
     itotalspe=igas+iaer*isize+2  !! plus 2 for u,v, by chenhs
     ! for Source Mark
     if(ifsmt>0)then
      itotalspe=itotalspe+ismMax*idmSet
      print *,'itotalspe=',itotalspe
     endif
  
    ! for sea salt and dust composition
    if(ifseacom.eq.1) then
      itotalspe = itotalspe + nseacom*isize 
    endif
    if(ifdustcom.eq.1) then
      itotalspe = itotalspe + ndustcom*isize
    endif

!--------------------------------------------------------------------
!    parallel
!--------------------------------------------------------------------   
  do ne= 1, nest
 
       dims(1,ne) = 0
       dims(2,ne) = 0
       call mpi_dims_create( numprocs, 2, dims(1,ne), ierr )
       call mpi_cart_create( mpicom_atm, 2, dims(1,ne),  &
			periods(1,ne), .true.,  &
			comm2d(ne), ierr )
       call mpi_comm_rank( comm2d(ne), myid, ierr )
       call mpi_cart_shift( comm2d(ne), 0, 1, nbrleft(ne), nbrright(ne), ierr )
       call mpi_cart_shift( comm2d(ne), 1, 1, nbrbottom(ne), nbrtop(ne), ierr )

       local_com = comm2d(ne)
       CALL mpi_cart_coords( comm2d(ne), myid, 2, coords(1,ne), ierr )
       myid_x = coords(1,ne)   ! col task (x)
       myid_y = coords(2,ne)   ! row task (y)       
       allocate(procs(0:dims(2,ne)-1,0:dims(1,ne)-1)) ! dims(2,ne) corresponds to ntasks_y
       procs(myid_y,myid_x)=myid       
       if (myid.eq.0) then
         do i=1,numprocs-1
            call MPI_Recv(buffer, 2, MPI_INTEGER, i, MPI_ANY_TAG, comm2d(ne), status, mpi_ierr)
            procs(buffer(1), buffer(2)) = status(MPI_SOURCE)
         end do
       else
         buffer(1) = myid_y
         buffer(2) = myid_x
         call MPI_Send(buffer, 2, MPI_INTEGER, 0, myid, comm2d(ne), mpi_ierr)
       end if

       allocate(cpuid(1:numprocs))
       if(myid.eq.0) then
       ig=0
       do i=0,dims(1,ne)-1
       do j=0,dims(2,ne)-1
         ig=ig+1
         cpuid(ig)=procs(j,i)
       end do
       end do
       endif
       call MPI_Bcast(cpuid, numprocs, MPI_INTEGER, 0, comm2d(ne), mpi_ierr)
       if(myid.ne.0) then
       ig=0
       do i=0,dims(1,ne)-1
       do j=0,dims(2,ne)-1
         ig=ig+1
         procs(j,i)=cpuid(ig)
       end do
       end do
       endif
       deallocate(cpuid)

!--------------------------------------------------------------------
!    Initial GEATM VARIABLES
!-------------------------------------------------------------------- 
       call initial_geatm_var
!--------------------------------------------------------------------
!    Initial WRFGRID and CAMGRID
!--------------------------------------------------------------------    
       call initial_wrfgrid(wrfgrid)
       call initial_camgrid(camgrid)
       call init_gea_pio

  !++++++++++++++++ chenhs,from geatm,east-west closure +++++++++++++++++++++++++++!
  if(ifglobal==1.and.ne==1.and.numprocs>1)then
  call exinfo(myid,numprocs,nbrleft(1),nbrright(1),nbrbottom(1),nbrtop(1),&
              stamark,memark)
  endif
 
  if(numprocs>1)call mpi_barrier( local_com,ierr )
  !++++++++++++++++ by chenhs,from geatm +++++++++++++++++++++++++++!
  call mpi_cart_get( comm2d(ne), 2, dims(1,ne), periods(1,ne),  &
                     coords(1,ne), ierr )
  call mpe_decomp1d( nx(ne), dims(1,ne), coords(1,ne), sx(ne), ex(ne) )
  call mpe_decomp1d( ny(ne), dims(2,ne), coords(2,ne), sy(ne), ey(ne) )

  !-----------------------for cam, added by juanxiong he---------------------------
  cnx(ne)=camgrid%num_camgrid_lon
  cny(ne)=camgrid%num_camgrid_lat
  if(ne.eq.1) then
  call mpe_decomp1d( cnx(ne), dims(1,ne), coords(1,ne), csx(ne), cex(ne) )
  call mpe_decomp1d( cny(ne), dims(2,ne), coords(2,ne), csy(ne), cey(ne) )
  end if
  !-----------------------for CAM, added by Juanxiong He---------------------------

  ! write the domain decomposition infomation
  cdnum='c000'
  write(cdnum(2:4),'(i3.3)')myid
  open(99,file='grid.dat'//cdnum,form='formatted')

!--------------------------------------------------------------------
!   Initial AVECT
!--------------------------------------------------------------------
       ! Initialize MCT gsMap, domain and attribute vectors
       ! initial processors map
       call geatm_SetgsMap_cam(camgrid, mpicom_atm, GEATMID, gsMap_ge ) 
       lsize = mct_gsMap_lsize(gsMap_ge, mpicom_atm)
       ! initial domain       
       call geatm_domain_cam( camgrid, lsize, gsMap_ge, dom_ge ) 
       
       !
       ! Initialize MCT attribute vectors
       !             
        call mct_aVect_init(x2c_c1, rList=seq_flds_x2ge_fields, lsize=lsize) 
        call mct_aVect_zero(x2c_c1)

        call mct_aVect_init(x2c_c2, rList=seq_flds_x2ge_fields, lsize=lsize) 
        call mct_aVect_zero(x2c_c2)

        call mct_aVect_init(c2x_c, rList=seq_flds_ge2x_fields, lsize=lsize)
        call mct_aVect_zero(c2x_c)     
      
!--------------------------------------------------------------------
!    Initial CB-IV Chemistry Schemes
!--------------------------------------------------------------------
	  
	if(myid == 0 )then
	cdnum2='d0'
	write(cdnum2(2:2),'(i1)')ne
	open(98,file='head.dat'//cdnum2,form='formatted')
	write(98,11) numprocs, nx(ne), ny(ne), nz(ne), &
	              ntbeg(ne), ntend(ne), isize, iaer
	iwritegas=0
	do ig=1,igas
	if(PrintGas(ig)==1)iwritegas=iwritegas+1
	enddo
	write(98,12)iwritegas
	12 format(1x,I6,3x,A,3x,A)
	iwritegas=0
	do ig=1,igas
	if(PrintGas(ig)==1)then
	iwritegas=iwritegas+1
	write(98,12)iwritegas, GC_NAME(ig), GC_Unit(ig)
	endif
	enddo
	
	! for Source Mark
	write(98,*)ifsmt
	if(ifsmt>0)then
	write(98,*)idmSet,iSrcDefined,ismMax,iHgtLMax
	do idm=1,idmSet
	write(98,*)igMark(idm),iaMarkAer(idm),iaMarkSiz(idm)
	enddo
	endif
	close(98)
	11 format(1x,7I6)
	endif
	
	write(99,10)myid, numprocs, sx(ne),ex(ne),ex(ne)-sx(ne)+1, &
	      sy(ne),ey(ne),ey(ne)-sy(ne)+1, ne
	10 format( "Process ", i3, " of ", i3, " sx-ex-nx: ",  &
	     3I4," sy-ey,ny: ",3I4, ' nest=',i4)
	enddo   ! end of ne
	close(99)
	
	if(numprocs .gt. 1) &
	call mpi_barrier( local_com, ierr )
	do i=0,numprocs-1
	cdnum='c000'
	write(cdnum(2:4),'(i3.3)')i
	open(99,file='grid.dat'//cdnum,form='formatted')
	do ne=1,nest
	read(99,10)ii, ii, sxc(i,ne),exc(i,ne),ii, &
	     syc(i,ne),eyc(i,ne),ii, ii
	enddo
	close(99)
	enddo
	
	if(myid==0)then
	open(99,file='grid.info')
	do ne=1,nest
	do i=0,numprocs-1
	write(99,*) i, 'process',  sxc(i,ne),exc(i,ne), syc(i,ne),eyc(i,ne)
	enddo
	enddo
	close(99)
	endif ! numprocs
  
!++++++++++++++++++++++++++ chenhs, polar transport ++++++++++++++++++++++++++++!
	  if(ifglobal.eq.-1) then
	   if(numprocs .gt. 1) call mpi_barrier( local_com,ierr )
	   if(myid==0)then
	   open(99,file='grid.polar')
	   do ne=1,nest
	   if(ne==1) then
	   do i=0,numprocs-1
	    !!!for south boundary
	    if(syc(i,ne)==1) then
	     do ixx=sxc(i,ne),exc(i,ne)
	       if(ixx.le.180) then
		  ixc=ixx+180
	       else
		  ixc=ixx-180
	       endif
	
	       do j=0,numprocs-1
		  if(syc(j,ne)==1) then
		    if(ixc.ge.sxc(j,ne).and.ixc.le.exc(j,ne)) then
		      write(99,*) i,j,ixx,ixc,syc(i,ne)
		    endif
		  endif
	       enddo
	     enddo
	    endif
	    !for north boundary
	    if(eyc(i,ne)==ny(ne)) then
	     do ixx=sxc(i,ne),exc(i,ne)
	       if(ixx.le.180) then
		  ixc=ixx+180
	       else
		  ixc=ixx-180
	       endif
	
	       do j=0,numprocs-1
		  if(eyc(j,ne)==ny(ne)) then
		    if(ixc.ge.sxc(j,ne).and.ixc.le.exc(j,ne)) then
		      write(99,*) i,j,ixx,ixc,eyc(i,ne)
		    endif
		  endif
	       enddo
	     enddo
	    endif
	
	    enddo
	   endif
	  enddo
	 close(99)
	 endif
	if(numprocs .gt. 1) call mpi_barrier( local_com,ierr )
	open(99,file='grid.polar',form='formatted')
	do i=1,ipolarnum
	   !ipolarmrk(5,720)
	   !j--> 1 (receive cpu)
	   !j--> 2 (send cpu)
	   !j--> 3 (reveive grid)
	   !j--> 4 (send grid)
	   !j--> 5 (boundary num,1 for south,180 for north)
	   read(99,*) (ipolarmrk(j,i),j=1,5)
	enddo  
	close(99)
	if(numprocs .gt. 1) call mpi_barrier( local_com,ierr )
	endif 
!++++++++++++++++++++++++++ chenhs, polar transport ++++++++++++++++++++++++++++!
!  to prepare my own send-rev
	do ne=2,nest   ! modified by chenhs for more children domain
	
	       lx0=nxlo(ne)
	       ly0=nylo(ne)
	       lx1=nxlo(ne)+nx(ne)/sratio(ne)
	       ly1=nylo(ne)+ny(ne)/sratio(ne)
	
	  
	  do iii=1,4    ! 1,south 2,north 3,west 4,east
	  bdysx(iii,ne) = 0
	  bdyex(iii,ne) = -1
	  bdysy(iii,ne) = 0
	  bdyey(iii,ne) = -1
	  enddo
	
	! south boundary
	   if(ly0 .ge. sy(imother(ne)) .and. ly0 .le. ey(imother(ne)))then
	    ifindit=1
	    ifindnum=0
	     do ikl=sx(imother(ne)),ex(imother(ne))
	     do ikk = lx0,lx1
	     if(ikk==ikl) then
		 ifindnum=ifindnum+1
		 if(ifindit==1)then
		    bdysx(1,ne)=ikl
		    ifindit=0
		   endif
	      endif
	     enddo
	     enddo
	     if(ifindnum>0) then
		bdyex(1,ne)=bdysx(1,ne)+ifindnum-1
		bdysy(1,ne)=ly0
		bdyey(1,ne)=ly0
	     endif
	  endif
	
	! north boundary
	    if(ly1 .ge. sy(imother(ne)) .and. ly1 .le. ey(imother(ne)))then
	    ifindit=1
	    ifindnum=0
	     do ikl=sx(imother(ne)),ex(imother(ne))
	     do ikk = lx0,lx1
	     if(ikk==ikl) then
		 ifindnum=ifindnum+1
		 if(ifindit==1)then
		    bdysx(2,ne)=ikl
		    ifindit=0
		   endif
	      endif
	     enddo
	     enddo
	     if(ifindnum>0) then
		bdyex(2,ne)=bdysx(2,ne)+ifindnum-1
		bdysy(2,ne)=ly1
		bdyey(2,ne)=ly1
	     endif
	  endif
	
	! West boundary
	   if(lx0 .ge. sx(imother(ne)) .and. lx0 .le. ex(imother(ne)))then
	    ifindit=1
	    ifindnum=0
	     do ikl=sy(imother(ne)),ey(imother(ne))
	     do ikk = ly0,ly1
	     if(ikk==ikl) then
		 ifindnum=ifindnum+1
		 if(ifindit==1)then
		    bdysy(3,ne)=ikl
		    ifindit=0
		   endif
	      endif
	     enddo
	     enddo
	     if(ifindnum>0) then
		bdyey(3,ne)=bdysy(3,ne)+ifindnum-1
		bdysx(3,ne)=lx0
		bdyex(3,ne)=lx0
	     endif
	  endif
	
	! East boundary
	   if(lx1 .ge. sx(imother(ne)) .and. lx1 .le. ex(imother(ne)))then
	    ifindit=1
	    ifindnum=0
	     do ikl=sy(imother(ne)),ey(imother(ne))
	     do ikk = ly0,ly1
	     if(ikk==ikl) then
		 ifindnum=ifindnum+1
		 if(ifindit==1)then
		    bdysy(4,ne)=ikl
		    ifindit=0
		   endif
	      endif
	     enddo
	     enddo
	     if(ifindnum>0) then
		bdyey(4,ne)=bdysy(4,ne)+ifindnum-1
		bdysx(4,ne)=lx1
		bdyex(4,ne)=lx1
	     endif
	  endif
	
	enddo ! ne
!---------------------------------------------------------------------
! to check the position(CPU) of the line in  nest domain
!---------------------------------------------------------------------
  ! now my id is myid
  if(numprocs .gt. 1) call mpi_barrier( local_com, ierr )

  cdnum='c000'
  write(cdnum(2:4),'(i3.3)')myid
  open(99,file='nest.bdy'//cdnum,form='formatted')

	do ne=2,nest  ! modified by chenhs for more children domain
	
	! for south and north boundary
	do iii=1,2
	   if( bdysy(iii,ne) .gt. 0 )then
	      if( bdysx(iii,ne).le. bdyex(iii,ne))then
	       write(99, 201) bdysy(iii,ne),bdysx(iii,ne),bdyex(iii,ne), &
			      bdyex(iii,ne)-bdysx(iii,ne)+1,&
			      iii,imother(ne),myid,ne  ! modified  by  chenhs
	      endif
	   endif
	enddo
	
	do iii=3,4
	
	   if( bdysx(iii,ne) .gt. 0 )then
	      if( bdysy(iii,ne).le. bdyey(iii,ne))then
	       write(99, 201) bdysx(iii,ne),bdysy(iii,ne),bdyey(iii,ne), &
			      bdyey(iii,ne)-bdysy(iii,ne)+1,&
			      iii,imother(ne),myid,ne
	       endif
	    endif
	 enddo
	
	enddo
	
	close(99)
	
	201 format(1x,12i5)
	
	  if(numprocs .gt. 1) &
	  call mpi_barrier( local_com, ierr )
	
	  nsndmrk=1
	  do i=0,numprocs-1
	  cdnum='c000'
	  write(cdnum(2:4),'(i3.3)')i
	  open(99,file='nest.bdy'//cdnum,form='formatted')
	  do ne=1,500
	  read(99,201,err=199,end=199)(isndmrk(j,nsndmrk),j=1,8) ! modified by chenhs, 8 for child domain
	  nsndmrk=nsndmrk+1
	  if(nsndmrk.gt.200)then
		   print *, '200 so small'
		   stop 200
		 endif
	  enddo
	  199  close(99)
	  enddo
	
	  nsndmrk=nsndmrk-1
	
	!---------- to send data from domain ne --> ne + 1
	  cdnum='c000'
	  write(cdnum(2:4),'(i3.3)')myid
	  open(99,file='locate.bdy'//cdnum,form='formatted')
	
	  do ir=1,nsndmrk
	     ! isndmrk(8,ir)
	     !         x--> 8 (ne)  child domain 
	     !         x--> 7 (cpu) 
	     !         x--> 6 mother domain
	     !         x--> 5 (boundary number mark: 1 south, 2 north, 3 west, 4 east)
	     !         x--> 4 ( grid numbers in this cpu need to send)
	     !         x--> 3 ( end point )
	     !         x--> 2 ( start point )
	     !         x--> 1 ( location, x or y depends on boundary numbert mark
	 
	    if( isndmrk(7,ir) == myid ) then    
		      !! this cpu has boundary condition to next domain
	      ne = isndmrk(6,ir)
	      iichild = isndmrk(8,ir)
	      !south boundary---------------------------------------
	      if(isndmrk(5,ir)==1)then   
		 do ipoint=isndmrk(2,ir),isndmrk(3,ir)   ! these big grids need to send
		  ! to put each point to other cpu
		  ipoint1b = (ipoint - nxlo(iichild))* sratio(iichild) + 1
		  ipoint1e = (ipoint - nxlo(iichild))* sratio(iichild) + sratio(iichild)
		  !    
		  do ipoint1=ipoint1b,ipoint1e
		    do icpu=0,numprocs-1   ! to check each cpu to receive
		      ! sxc(icpu,ne+1)-->exc(icpu,ne+1)
		      ! syc(icpu,ne+1)-->eyc(icpu,ne+1)
		      if( syc(icpu,iichild) .eq. 1 )then
			if(ipoint1.ge.sxc(icpu,iichild) .and.  &
			   ipoint1.le. exc(icpu,iichild))then 
				write(99,204)myid,icpu,ipoint,isndmrk(1,ir),&
					     isndmrk(5,ir),isndmrk(6,ir),isndmrk(8,ir), &
					     ipoint1, 0
					     !ipoint1, syc(icpu,ne+1)
			endif 
		      endif
		    enddo                  ! to check each cpu to receive
		  enddo 
	
		 enddo                                  ! these big grids need to send
	      endif
	     !north boundary---------------------------------------
	      if(isndmrk(5,ir)==2)then
		 do ipoint=isndmrk(2,ir),isndmrk(3,ir)
		  ! to put each point to other cpu
		  ipoint1b = (ipoint - nxlo(iichild))* sratio(iichild) + 1
		  ipoint1e = (ipoint - nxlo(iichild))* sratio(iichild) + sratio(iichild)
		  !    
		  do ipoint1=ipoint1b,ipoint1e
		    do icpu=0,numprocs-1
		      ! sxc(icpu,ne+1)-->exc(icpu,ne+1)
		      ! syc(icpu,ne+1)-->eyc(icpu,ne+1)
		      if( eyc(icpu,iichild) .eq. ny(iichild) )then
			if(ipoint1.ge.sxc(icpu,iichild) .and.  &
			   ipoint1.le. exc(icpu,iichild))then
				write(99,204)myid,icpu,ipoint,isndmrk(1,ir),&
					     isndmrk(5,ir),isndmrk(6,ir),isndmrk(8,ir), &
					     ipoint1, ny(iichild)+1
					     !ipoint1, syc(icpu,ne+1)
			endif
		      endif
		    enddo           
		   enddo           
		 enddo 
	      endif
	
	      if(isndmrk(5,ir)==3 )then   ! west boundary 
		 do jpoint=isndmrk(2,ir),isndmrk(3,ir)
		    ! to put each point to other cpu
		  jpoint1b = (jpoint - nylo(iichild))* sratio(iichild) + 1
		  jpoint1e = (jpoint - nylo(iichild))* sratio(iichild) + sratio(iichild)
		  do jpoint1 = jpoint1b,  jpoint1e
		    do icpu=0,numprocs-1
		      ! sxc(icpu,ne+1)-->exc(icpu,ne+1)
		      ! syc(icpu,ne+1)-->eyc(icpu,ne+1)
		      if( sxc(icpu,iichild) .eq. 1) then
			if(jpoint1.ge.syc(icpu,iichild) .and.  &
			   jpoint1.le. eyc(icpu,iichild))then
				write(99,204)myid,icpu,isndmrk(1,ir),jpoint,&
					     isndmrk(5,ir),isndmrk(6,ir),isndmrk(8,ir), &
					     0,jpoint1
					     !sxc(icpu,ne+1),jpoint1
			endif
		      endif
		    enddo
		    enddo
		 enddo
	      endif
	
	      if(isndmrk(5,ir)==4 )then   ! east boundary 
		 do jpoint=isndmrk(2,ir),isndmrk(3,ir)
		    ! to put each point to other cpu
		  jpoint1b = (jpoint - nylo(iichild))* sratio(iichild) + 1
		  jpoint1e = (jpoint - nylo(iichild))* sratio(iichild) + sratio(iichild)
		  do jpoint1 = jpoint1b,  jpoint1e
		    do icpu=0,numprocs-1
		      ! sxc(icpu,ne+1)-->exc(icpu,ne+1)
		      ! syc(icpu,ne+1)-->eyc(icpu,ne+1)
		      if( exc(icpu,iichild) .eq. nx(iichild)) then
			if(jpoint1.ge.syc(icpu,iichild) .and.  &
			   jpoint1.le. eyc(icpu,iichild))then
				write(99,204)myid,icpu,isndmrk(1,ir),jpoint,&
					     isndmrk(5,ir),isndmrk(6,ir),isndmrk(8,ir), &
					     nx(iichild)+1,jpoint1
					     !sxc(icpu,ne+1),jpoint1
			endif
		      endif
		    enddo
		   enddo
		 enddo
	      endif
	    endif
	  enddo  ! ne      
    204 format(1x,9I6)
    close(99)

	  if(numprocs .gt. 1) call mpi_barrier( local_com, ierr )
	
	  npsr=1   !numbers need to send and receive
	  do i=0,numprocs-1
	  cdnum='c000'
	  write(cdnum(2:4),'(i3.3)')i
	  open(99,file='locate.bdy'//cdnum,form='formatted')
	  do np=1,2000
	  read(99,204,err=399,end=399)iiscpu(npsr),iircpu(npsr),&
	     islocx(npsr),islocy(npsr),issnwe(npsr),isnest(npsr),irnest(npsr), &
	     irlocx(npsr),irlocy(npsr)
	     ! isndmrk(8,ir)
	     !         x--> 8 (ne)  child domain   
	     !         x--> 7 (cpu) 
	     !         x--> 6 mother domain  
	     !         x--> 5 (boundary number mark: 1 south, 2 north, 3 west, 4 east)
	     !         x--> 4 ( grid numbers in this cpu need to send)
	     !         x--> 3 ( end point )
	     !         x--> 2 ( start point )
	     !         x--> 1 ( location, x or y depends on boundary numbert mark 
	  npsr = npsr + 1
	  if(npsr.gt.4000)then
		   print *, '4000 so small'
		   stop 4000
		 endif
	  enddo
	  399  close(99)
	  enddo

	  npsr = npsr - 1
	  
	  if(myid == 0 )then
	      open(99,file='needsend.dat')
	      write(99,*)npsr
	      do ne=1,nest-1
	      do ib=1,4
	      do np=1,npsr
	      if(isnest(np)==ne .and. issnwe(np) == ib )then
	      write(99,204)iiscpu(np),iircpu(np),islocx(np),islocy(np), &
			   issnwe(np),isnest(np),irnest(np),irlocx(np),irlocy(np) 
	      endif
	      enddo
	      enddo
	      enddo
	      close(99)
	  endif
	       
           first_time=.false.    

else
!---------------------------------------------------------------------------    
!  geat or read the initial data, phase = 2
!---------------------------------------------------------------------------

    !----------------------------------------------------
    ! added by juanxiong he
    !----------------------------------------------------
    open(1000,file='wrfd01.dat',form='unformatted',&
         access='direct',recl=nx(1)*ny(1),status='old')
    irechgt=1
    ! to read SOIL TYPE
    i0=ip2mem(1)
    call read2d(myid,FSOIL(i0),sx(1),ex(1),sy(1),ey(1), &
              nx(1),ny(1),irechgt,1000)

    ! to read VEGETATION FRACTION (%)
    i0=ip2mem(1)
    call read2d(myid,FVEG(i0),sx(1),ex(1),sy(1),ey(1), &
              nx(1),ny(1),irechgt,1000)

    ! read terrain
    i0=ip2mem(1)
    call read2d(myid,HGT1(i0),sx(1),ex(1),sy(1),ey(1), &
              nx(1),ny(1),irechgt,1000)
       do j= sy(1),ey(1)
       do i= sx(1),ex(1)
           km=i0+(ex(1)-sx(1)+3)*(j-sy(1)+1)+i-sx(1)+1
           wrfgrid%ht(i,j)=HGT1(km)*1.0_8   
        end do
       end do
    ! the cubic spline interpolation needs the value varys from the large to the small
    ! Firstly the program reverses the wrfgrid from surface-to-upper to upper-to-surface
    ! for the sake of cubic spline interpolation.
    ! At subroutine wrfgrid_to_geatm, the program will reverse the export of wrfgrid to geatm from
    ! upper-to-surfaceto surface-to-upper.
    do j=1,ny(1)
    do k=1,wrfgrid%num_wrfgrid_levels
    do i=1,nx(1)
       wrfgrid%z3d(i,k,j) = wrfgrid%sigma(wrfgrid%num_wrfgrid_levels-k+1)*(20000.0-wrfgrid%ht(i,j))+wrfgrid%ht(i,j)
    enddo
    enddo
    enddo

    ! read landmask, 1 for land, 0 for sea
    i0=ip2mem(1)
    call read2d(myid,landmask(i0),sx(1),ex(1),sy(1),ey(1), &
              nx(1),ny(1),irechgt,1000)
    close(1000)

    call geatm_import_mct(x2c_c1, camgrid, loop)
    !----------------------------------------------------
    ! added by juanxiong he
    !----------------------------------------------------

           do ne=1,nest
                  call openfil1(myid,nx(ne),ny(ne),irec80(ne),irec60(ne),ne, &
                                iyear1,imonth1,idate1,ihour1)
                  if(numprocs .gt. 1)then
                   call mpi_type_vector( ey(ne)-sy(ne)+3, 1, ex(ne)-sx(ne)+3, &
                                mpi_real, stride(ne), ierr )
                   call mpi_type_commit( stride(ne), ierr )
                  endif
           enddo
	
	do ne=1,nest
	if(irec80(ne) .gt. 0 ) then
	  ! gas 	 
	  do ig=1,igas
	    if(InitPrintGas(ig)==1)then  ! by chenhs
	     do k=1,nzz
	     i0=ip4mem(k,ig,ne)
	     call read2d(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		      nx(ne),ny(ne),irec80(ne),80+ne)
	     enddo
	    endif
	  enddo
	  ! o1
	  if(ne.eq.0) then  ! by chenhs
	  do k=1,nzz
	  i0=ip3mem(k,ne)
	  call read2d(myid,jo1d(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		nx(ne),ny(ne),irec80(ne),80+ne)
	  enddo
	  ! no2
	  do k=1,nzz
	  i0=ip3mem(k,ne)
	  call read2d(myid,jno2(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		nx(ne),ny(ne),irec80(ne),80+ne)
	  enddo
	  endif
	  ! EXT
	  do k=1,nzz
	  i0=ip3mem(k,ne)
	  call read2d(myid,EXT(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		nx(ne),ny(ne),irec80(ne),80+ne)
	  enddo              
	  ! VISIB 
	  do k=1,nzz
	  i0=ip3mem(k,ne)
	  call read2d(myid,VISIB(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		nx(ne),ny(ne),irec80(ne),80+ne)
	  enddo
	
	  if(ne.eq.0) then
	  ! UVB	  
	  do k=1,nzz
	  i0=ip3mem(k,ne)
	  call read2d(myid,UVB(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		nx(ne),ny(ne),irec80(ne),80+ne)
	  enddo
	  ! UVBS
	  do k=1,nzz
	  i0=ip3mem(k,ne)
	  call read2d(myid,UVBS(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		nx(ne),ny(ne),irec80(ne),80+ne)
	  enddo
	  ! UVA
	  do k=1,nzz
	  i0=ip3mem(k,ne)
	  call read2d(myid,UVA(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		nx(ne),ny(ne),irec80(ne),80+ne)
	  enddo
	  ! VIS
	  do k=1,nzz
	  i0=ip3mem(k,ne)
	  call read2d(myid,VIS(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		nx(ne),ny(ne),irec80(ne),80+ne)
	  enddo  
	  endif !ne, chenhs 
	  !SSA
	  do k=1,nzz 
	  i0=ip3mem(k,ne)
	  call read2d(myid,SSA(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		nx(ne),ny(ne),irec80(ne),80+ne)
	  enddo
	  ! AOD
	  i0=ip2mem(ne)
	  call read2d(myid,AOD(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
	       nx(ne),ny(ne),irec80(ne),80+ne)
	  ! CLDOPD 
	  i0=ip2mem(ne)
	  call read2d(myid,CLDOPD(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
	       nx(ne),ny(ne),irec80(ne),80+ne)
	  ! DUSO2     
	  i0=ip2mem(ne)
	  call read2d(myid,DUSO2(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
	       nx(ne),ny(ne),irec80(ne),80+ne)
	  ! DUO3
	  i0=ip2mem(ne)
	  call read2d(myid,DUO3(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
	       nx(ne),ny(ne),irec80(ne),80+ne)
	  ! DUNO2      
	  i0=ip2mem(ne)
	  call read2d(myid,DUNO2(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
	       nx(ne),ny(ne),irec80(ne),80+ne)
	  ! ANA
	  do k=1,nzz
	  i0=ip3mem(k,ne)
	  call read2d(myid,ANA(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		nx(ne),ny(ne),irec80(ne),80+ne)
	  enddo
	  ! ASO4
	  do k=1,nzz
	  i0=ip3mem(k,ne)
	  call read2d(myid,ASO4(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		nx(ne),ny(ne),irec80(ne),80+ne)
	  enddo
	  ! ANH4
	  do k=1,nzz
	  i0=ip3mem(k,ne)
	  call read2d(myid,ANH4(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		nx(ne),ny(ne),irec80(ne),80+ne)
	  enddo
	  ! ANO3
	  do k=1,nzz
	  i0=ip3mem(k,ne)
	  call read2d(myid,ANO3(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		nx(ne),ny(ne),irec80(ne),80+ne)
	  enddo
	  ! ACL
	  do k=1,nzz
	  i0=ip3mem(k,ne)
	  call read2d(myid,ACL(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		nx(ne),ny(ne),irec80(ne),80+ne)
	  enddo
	  ! CPH             
	  do K = 1, NZZ
	  i0=ip3mem(k,ne)
	  call read2d(myid,CPH(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		nx(ne),ny(ne),irec80(ne),80+ne)
	  enddo
	  ! OPE
	  do K = 1, NZZ
	  i0=ip3mem(k,ne)
	  call read2d(myid,OPE(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		nx(ne),ny(ne),irec80(ne),80+ne)
	  enddo
	  ! AER       
	  do ia=1,iaer
	  do is=1,isize
	  do k=1,nzz
	  i0=ip5mem(k,is,ia,ne)
	     call read2d(myid,aer(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		      nx(ne),ny(ne),irec80(ne),80+ne)
	  enddo
	  enddo
	  enddo
	
	  ! ---- read dust init data
	  if(ne.eq.0) then
	  do ia=2,2
	  do is=1,isize
	  do k=1,nzz
	  i0=ip5mem(k,is,ia,ne)
	     call read2d(myid,aer(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		      nx(ne),ny(ne),irec60(ne),60+ne)
	  enddo
	  enddo
	  enddo
	  endif !! ne,chenhs
	  if(ifdustcom.eq.1) then
	  do iduc = 1, ndustcom
	  do is = 1, isize
	  do k =1 ,nzz
	     i0=ip5memc(k,is,iduc,ne)
	     call read2d(myid,dustcomp(i0),sx(ne),ex(ne),sy(ne),ey(ne),&
		   nx(ne),ny(ne),irec60(ne),60+ne)
	  enddo
	  enddo
	  enddo
	  endif
	 else
	  !the first choise is to use MOZART field for initial conditions, by chenhs!!
	  call openfilMOZART(myid,nx(ne),ny(ne),irecMOZART(ne),ne, &
				iyear1,imonth1,idate1,ihour1)
	  IF(irecMOZART(ne).GT.0) THEN
	    ! to read global concentration of NO2
	    do k=1,nzz
	       i0=ip4mem(k,6,ne)
	       call read2d(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne),&
			   nx(ne),ny(ne),irecMOZART(ne),350+ne )
	    enddo
	    ! to read global concentration of O3
	    do k=1,nzz
	       i0=ip4mem(k,11,ne)
	       call read2d(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne),&
			   nx(ne),ny(ne),irecMOZART(ne),350+ne )
	    enddo
	    ! to read global concentration of CO
	    do k=1,nzz
	       i0=ip4mem(k,17,ne)
	       call read2d(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne),&
			   nx(ne),ny(ne),irecMOZART(ne),350+ne )
	    enddo
	  ELSE
	   ! to ozone as 60ppb vertically
	   do ig=11,11
	     do k=1,nzz
	     i0=ip4mem(k,ig,ne)
	      if(k==1.or.k==2.or.k==3) &
		call puto3(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
			   nx(ne),ny(ne),20.)
	      if(k==4.or.k==5.or.k==6) &
		call puto3(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
			   nx(ne),ny(ne),10.)
	      if(k==7.or.k==8.or.k==9) &
		call puto3(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		      nx(ne),ny(ne),20.)
	      if(k>=10)  &
		call puto3(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		      nx(ne),ny(ne),30.)
	
	     enddo
	   enddo
	
	   do ig=17,17  !to CO as 200ppbv
	     do k=1,nzz
	     i0=ip4mem(k,ig,ne)
	     call puto3(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		      nx(ne),ny(ne),300.)
	     enddo
	   enddo
	  ENDIF !! if use MOZART initial 
	 
	 if(ne.eq.0) then !! by chenhs 
	 DO IA = 1, IAER  ! TO DUST AND SEA SALT 5 UG/M3
	 DO IS = 1, ISIZE
	 DO K =1 ,NZZ
	    I0 = IP5MEM(K,IS,IA,NE)
	  CALL PUTO3 (MYID, AER(I0), SX(NE), EX(NE), SY(NE), EY(NE),&
		       NX(NE),NY(NE),0.1)
	 ENDDO ! K
	 ENDDO ! IS
	 ENDDO ! IA
	 endif
	     !!!!!!!! put HG0 initial conditions, North Hemispere 1.6 ng/m3 !!!!!!
	     !!!!!!!! South Hemispere 1.2 ng/m3, decrease with altitude     !!!!!!
	     do ig=103,103
		do k=1,nzz
		   i0=ip4mem(k,ig,ne)
		   call putHG0(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
			       nx(ne),ny(ne),ne,k)
		enddo
	     enddo  

            close(350+ne)  ! for MOZART, juanxiong he
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	endif   ! if no exist of the init file, no input, irec80(ne).
	enddo
	
	 !!!!!!!!!!!!!!!!!!!!!
	 ! For Source Mark
	 if(ifsmt>0)then 
	 do ne=1,nest
	 if(ifsm(ne)==1)then !mmmmm
	  call openfilSM(myid,nx(ne),ny(ne),irecSM(ne),ne, &
			 iyear1,imonth1,idate1,ihour1)
	  if(irecSM(ne) .gt. 0 ) then !xxxxxxxxxxxxxxx
	     do idm=1,idmSet
	     do ism=1,ismMax  !ism
	      do k=1,nzz
	      i0=ipSMmem(k,ism,idm,ne)
	     call read2d(myid,SourceMark(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		       nx(ne),ny(ne),irecSM(ne),20+ne)
	      enddo
	     enddo            !ism
	      do k=1,nzz
	      ig=igMark(idm)
	      i0=ip4mem(k,ig,ne)
	      call read2d(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		       nx(ne),ny(ne),irecSM(ne),20+ne)
	      enddo
	     enddo
	   else
	     do idm=1,idmSet
	     do ism=1,ismMax
		do k=1,nzz-1
		i0=ipSMmem(k,ism,idm,ne)
		   if(ism==3)then   ! initial 100% (1 boundary,2 strato, 3 init)
	     call puto3(myid,SourceMark(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		      nx(ne),ny(ne),1.)
		   else
	     call puto3(myid,SourceMark(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		      nx(ne),ny(ne),0.)
		   endif
		 enddo !k
	! set top conditions 100%
	       do k=nzz,nzz
		   i0=ipSMmem(k,ism,idm,ne)
		  if(ism==2) then
	     call puto3(myid,SourceMark(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
			    nx(ne),ny(ne),1.) 
		  else
	     call puto3(myid,SourceMark(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		      nx(ne),ny(ne),0.)     
		  endif
	       enddo !k
	      enddo !ism
	     enddo !idm
	    endif                ! irecsm
           close(20+ne) ! juanxiong he
	  endif                  ! ifsm==1
	  enddo 
	  endif                  ! ifsm>0

          close(60+ne) ! juanxiong he
          close(80+ne) ! juanxiong he
!-------------------------------------------------------------
!============  Zifa add write initial condition =====
!===========   2006-07-07 
!-------------------------------------------------------------	
	do ne=1,nest
	
	iitime = (ntbeg(ne)-1)*3600    ! in seconds, need to modify, chenhs
	call getnewdate(iyear1,imonth1,idate1,ihour1,iitime, & 
			iyear2,imonth2,iday2, ihour2,iminute2)

	 !!!!!!!!!!!!!!!!!!!!!!!!!!
	 ! for Source Mark
	  if(ifsm(ne)==1)then 
	  ! to write source file <only save 3 levels> 1,3,5
           call output_mark(iyear2,imonth2,iday2,ihour2,ne,.false.)
	  IF(ihour2==0)THEN  !! to write down reinit source file
           call output_mark(iyear2,imonth2,iday2,ihour2,ne,.true.)
	  ENDIF
	  endif ! ifsm==1
	  
	enddo   ! ne, write for different domain

	end if ! first time

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

end subroutine geatm_init_mct

!================================================================================
	
subroutine geatm_run_mct(EClock_aa, EClock, cdata_a, x2c_c1, x2c_c2, c2x_c)
    use geatm_vartype
    use shr_sys_mod, only: shr_sys_flush
    
    ! 
    ! Arguments
    !
    type(ESMF_Clock)            ,intent(in)    :: EClock_aa ! cam
    type(ESMF_Clock)            ,intent(in)    :: EClock  ! geatm
    type(seq_cdata)             ,intent(inout) :: cdata_a ! geatm
    type(mct_aVect)             ,intent(inout) :: x2c_c1   ! cam -> cpl -> geatm  
    type(mct_aVect)             ,intent(inout) :: x2c_c2   ! cam -> cpl -> geatm
    type(mct_aVect)             ,intent(inout) :: c2x_c   ! geatm -> cpl -> cam
    
    include 'mpif.h'
    include 'chm1.inc'
    include 'gas1.inc'
    include 'params1'
    include 'tuv.inc'
    integer status(mpi_status_size)
    integer :: ymd, tod, curr_ymd, curr_tod, dtime
    integer,save ::it1 =1
    logical :: dosend

    ! Map input from mct to geatm data structure
    call geatm_import_mct( x2c_c1, camgrid, loop )
    ! get the next coupling time of CAM 
    call seq_timemgr_EClockGetData(EClock_aa, curr_ymd=curr_ymd, curr_tod=curr_tod)
    ! get the next coupling time of GEATM
    call seq_timemgr_EClockGetData(EClock, curr_ymd=ymd, curr_tod=tod, dtime=dtime)
    print *,'GEATM=',ymd,tod,'CAM=',curr_ymd,curr_tod
    ! get current time of GEATM
    call seq_timemgr_EClockGetData(EClock, prev_ymd=ymd, prev_tod=tod)
    iyear3=ymd/10000
    imonth3=(ymd-ymd/10000*10000)/100
    iday3=ymd-ymd/100*100
    ihour3=tod/3600

    call output_atm(iyear3,imonth3,iday3,ihour3,1) ! atm

    ! to integration for each time step         
    dosend = .false.
    do while (.not. dosend)

    iitime = (it1-1)*3600    ! in seconds
    call getnewdate(iyear1,imonth1,idate1,ihour1,iitime, & 
                    iyear2,imonth2,iday2, ihour2,iminute2)
    ymd = iyear2*10000 + imonth2*100 + iday2
    tod = ihour2*3600+iminute2*60
    dosend = (seq_timemgr_EClockDateInSync( EClock, ymd, tod))

    ! to read data
    do ne=1,nest ! nest

    if(mod(ihour2,dtmet(ne))==0) then  ! added by chenhs,to adjust met and emit input frequency 
    call openfil2(myid,nx(ne),ny(ne),irec(ne),irec_as(ne),&
                  irecglobal(ne),irechgt,ne, &
                  iyear2,imonth2,iday2,ihour2,ikosaline)  !added by lijie0500602

    if(it1.ge.ntbeg(ne) .and. it1.le.ntend(ne))then
    !---------------------lijie modify to global condition------
    ! to read global concentration of no2
    if(iglobal==1.and.mod(ihour2,1)==0) then !every 3hrs(0,3,6,9,12,15,18,21)
            do k=1,nzz
               i0=ip3mem(k,ne)
                   call read2d(myid,globalno2(i0),sx(ne),ex(ne),sy(ne),ey(ne),&
                        nx(ne),ny(ne),irecglobal(ne),300+ne )
            enddo
    ! to read global concentration of o3
            do k=1,nzz
               i0=ip3mem(k,ne)
                   call read2d(myid,globalo3(i0),sx(ne),ex(ne),sy(ne),ey(ne),&
                        nx(ne),ny(ne),irecglobal(ne),300+ne)
             enddo
   ! to read the  concentration of co
             do k=1,nzz
               i0=ip3mem(k,ne)
                   call read2d(myid,globalco(i0),sx(ne),ex(ne),sy(ne),ey(ne),&
                         nx(ne),ny(ne),irecglobal(ne),300+ne)
             enddo
    endif  !global
    close(300+ne) ! juanxionge he
    !-------------------------------finish----------------------
   	
   !zifa 2006/07/13/B
   iitime = (it1-1)*3600    ! in seconds
   call getnewdate(iyear1,imonth1,idate1,ihour1,iitime, &
                   iyear2,imonth2,iday2, ihour2,iminute2)
 
   imonthEmit = imonth1 ! added juanxiong he

   if(it1==ntbeg(ne) .or. imonthEmit.ne.imonth2) then  !Chenhs added
  
    !!!!!!  to get the information emissions type and species 
    close(999)
    open(999,file='emitinput.dat')
    read(999,*)  nemit    ! the number of emission species
    read(999,*)  nlay_em  ! the type of emission 
    allocate(ipig(nemit))
    read(999,*) (ipig(i),i=1,nemit) ! the ig ofemissions spece in model
 
    read(999,*)  ifnox        ! if add soil and lighting nox
    if(ifnox.eq.1) then
    read(999,*)  nemit_nox    ! the number of emission species 
    read(999,*)  nlay_nox     ! the type of nox emission
    allocate(ipig_nox(nemit_nox))
    read(999,*) (ipig_nox(i),i=1,nemit_nox) ! the ig of emissions spece in model
    endif

    read(999,*)  ifairc        ! if add aircraft emission
    if(ifairc.eq.1) then
    read(999,*)  nemit_airc    ! the number of emission species 
    read(999,*)  nlay_airc     ! the layers of aircraft emission
    allocate(ipig_airc(nemit_airc))
    read(999,*) (ipig_airc(i),i=1,nemit_airc) ! the ig of emissions spece in model
    endif

    if(ifhg.eq.1) then 
    read(999,*)  ifhga         ! if add hg anthro emission
    endif
    if(ifhga.eq.1) then
    read(999,*)  nemit_hga    ! the number of hg anthro emission species 
    read(999,*)  nlay_hga     ! the layers of hg anthro emission
    allocate(ipig_hga(nemit_hga))
    read(999,*) (ipig_hga(i),i=1,nemit_hga) ! the ig of hg anthro emissions spece in model 
    endif

    if(ifhg.eq.1) then
    read(999,*)  ifhgn         ! if add hg natural emission
    endif
    if(ifhgn.eq.1) then
    read(999,*)  nemit_hgn    ! the number of hg natural emission species 
    read(999,*)  nlay_hgn     ! the type of hg natural emission
    allocate(ipig_hgn(nemit_hgn))
    read(999,*) (ipig_hgn(i),i=1,nemit_hgn) ! the ig of hg natural emissions spece in model 
    endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! For Source Mark
  if(it1==ntbeg(ne) .and. ifsm(ne)==1)then  ! to get source/map
      call openfileSrcMap(myid,nx(ne),ny(ne),ne)
  i0=ip2mem(ne)
  call read2dmap(myid,MapSource(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),110+ne)
      close(110+ne) ! juanxiong he
  endif 
  !!!!!!!!!!!!!!!!!!!!!!!!!!

      if(ne==nest) imonthEmit = imonth2
      call openfileEmit(myid,nx(ne),ny(ne),irecg(ne),irecnox(ne),irecairc(ne),irecHgA(ne),irecHgN(ne),ne, &
                        iyear2,imonth2,iday2,ihour2,ifnox,ifairc,ifHgA,ifHgN)
  !zifa 2006/07/13/E
           ! to read dx
            do k=1,nzz
               i0=ip3mem(k,ne)
               call read2d(myid,dx(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
                           nx(ne),ny(ne),irecg(ne),50+ne)
            enddo
           ! to read dy
            do k=1,nzz
               i0=ip3mem(k,ne)
               call read2d(myid,dy(i0),sx(ne),ex(ne),sy(ne),ey(ne),&
                           nx(ne),ny(ne),irecg(ne),50+ne)
            enddo
           ! to read dz
            do k=1,nzz
               i0=ip3mem(k,ne)
               call read2d(myid,dz(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne),irecg(ne),50+ne)
            enddo
           ! to read heiz
            do k=1,nzz
               i0=ip3mem(k,ne)
               call read2d(myid,heiz(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne),irecg(ne),50+ne)
            enddo
	  i0=ip2mem(ne)
	  call read2d(myid,TERRAIN(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		      nx(ne),ny(ne),irecg(ne),50+ne)
	  call read2d(myid,LATITCRS(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		      nx(ne),ny(ne),irecg(ne),50+ne)
	  call read2d(myid,LONGICRS(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		      nx(ne),ny(ne),irecg(ne),50+ne)
	  call read2d(myid,LAND_USE(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
		      nx(ne),ny(ne),irecg(ne),50+ne)

     do img=1,NEMIT!   main emissions
       ig= IPIG(img)
        if(ig>0) then
       i0=ip2memGas(ig,ne)
       call read2d(myid,EmtAntGas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
               nx(ne),ny(ne),irecg(ne),50+ne)
         IF(NLAY_EM>=2) THEN        
       call read2d(myid,EmtShpGas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
               nx(ne),ny(ne),irecg(ne),50+ne)
         ENDIF
         IF(NLAY_EM>=3) THEN
       call read2d(myid,EmtBBGas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
               nx(ne),ny(ne),irecg(ne),50+ne)
         ENDIF     
         IF(NLAY_EM>=4) THEN
       call read2d(myid,EmtBioGas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
               nx(ne),ny(ne),irecg(ne),50+ne)
         ENDIF
         IF(NLAY_EM>=5) THEN
       call read2d(myid,EmtOceGas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
               nx(ne),ny(ne),irecg(ne),50+ne)
         ENDIF 
        endif
        !!!!!!!!!!!    chenhs to adjust emission    !!!!!!!!!!!
       !!!!!!!!!!! anthro and biomass NOx cut off 50%!!!!!!!!!
        IF(ig.eq.5.or.ig.eq.6) THEN
        do i=sx(ne),ex(ne)
        do j=sy(ne),ey(ne)
           ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
           i0e=ip2memGas(ig,ne)
           EmtAntGas(i0e+ixy)=EmtAntGas(i0e+ixy)*1.1  !! ant NOx +10%
        enddo
        enddo
        ENDIF
        IF( (img.ge.10.and.img.le.23).or.img.eq.26 ) THEN
        do i=sx(ne),ex(ne)
        do j=sy(ne),ey(ne)
           ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
           i0e=ip2memGas(ig,ne)
           EmtAntGas(i0e+ixy)=EmtAntGas(i0e+ixy)*0.9  !! ant NMVOC -10%
           EmtBioGas(i0e+ixy)=EmtBioGas(i0e+ixy)*0.9  !! BVOC -10%
        enddo
        enddo
        ENDIF
        IF( img.eq.22.or.img.eq.23 ) THEN  !! ISOP and TERP
        do i=sx(ne),ex(ne)
        do j=sy(ne),ey(ne)
           ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
           i0e=ip2memGas(ig,ne)
           EmtAntGas(i0e+ixy)=EmtAntGas(i0e+ixy)/12.
           EmtShpGas(i0e+ixy)=EmtShpGas(i0e+ixy)/12.
           EmtBBGas(i0e+ixy)=EmtBBGas(i0e+ixy)/12.
           EmtBioGas(i0e+ixy)=EmtBioGas(i0e+ixy)/12.
           EmtOceGas(i0e+ixy)=EmtOceGas(i0e+ixy)/12.
        enddo
        enddo
        ENDIF
        IF( img.eq.13 ) THEN  !! ETH*0.6
        do i=sx(ne),ex(ne)
        do j=sy(ne),ey(ne)
           ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
           i0e=ip2memGas(ig,ne)
           EmtAntGas(i0e+ixy)=EmtAntGas(i0e+ixy)*0.6
           EmtShpGas(i0e+ixy)=EmtShpGas(i0e+ixy)*0.6
           EmtBBGas(i0e+ixy)=EmtBBGas(i0e+ixy)*0.6
           EmtBioGas(i0e+ixy)=EmtBioGas(i0e+ixy)*0.6
           EmtOceGas(i0e+ixy)=EmtOceGas(i0e+ixy)*0.6
        enddo
        enddo
        ENDIF
        IF( img.eq.19 ) THEN  !! TOL*0.7
        do i=sx(ne),ex(ne)
        do j=sy(ne),ey(ne)
           ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
           i0e=ip2memGas(ig,ne)
           EmtAntGas(i0e+ixy)=EmtAntGas(i0e+ixy)*0.7
           EmtShpGas(i0e+ixy)=EmtShpGas(i0e+ixy)*0.7
           EmtBBGas(i0e+ixy)=EmtBBGas(i0e+ixy)*0.7
           EmtBioGas(i0e+ixy)=EmtBioGas(i0e+ixy)*0.7
           EmtOceGas(i0e+ixy)=EmtOceGas(i0e+ixy)*0.7
        enddo
        enddo
        ENDIF
        IF( img.eq.20 ) THEN  !! XYL*0.05
        do i=sx(ne),ex(ne)
        do j=sy(ne),ey(ne)
           ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
           i0e=ip2memGas(ig,ne)
           EmtAntGas(i0e+ixy)=EmtAntGas(i0e+ixy)*0.05
           EmtShpGas(i0e+ixy)=EmtShpGas(i0e+ixy)*0.05
           EmtBBGas(i0e+ixy)=EmtBBGas(i0e+ixy)*0.05
           EmtBioGas(i0e+ixy)=EmtBioGas(i0e+ixy)*0.05
           EmtOceGas(i0e+ixy)=EmtOceGas(i0e+ixy)*0.05
        enddo
        enddo
        ENDIF        
    enddo !img
       close(50+ne)  ! juanxiong he

    IF(ifnox.eq.1) then  !!SOIL and Lignting NOx emissions
      do img=1,NEMIT_NOx  
         ig=IPIG_NOx(img)
         if(ig>0) then
           i0=ip2memGas(ig,ne)
           call read2d(myid,EmtSoiGas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
                       nx(ne),ny(ne),irecnox(ne),100+ne)
           IF(NLAY_NOx>=2) THEN
           call read2d(myid,EmtLigGas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
                       nx(ne),ny(ne),irecnox(ne),100+ne)
           ENDIF
         endif
      enddo
    ENDIF
             close(100+ne) ! juanxiong he

    IF(ifairc.eq.1) then  !! aircraft emissions
      do img=1,NEMIT_AIRC
         ig=IPIG_AIRC(img)
         if(ig>0) then
           do k=1,NLAY_AIRC 
             i04=ip4mem(k,ig,ne)
             call read2d(myid,EmtAircGas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne),irecairc(ne),150+ne)
           enddo
         endif
      enddo
    ENDIF
             close(150+ne) ! juanxiong he

    IF(ifHg.eq.1) THEN
    IF(ifHgA.eq.1) then  !! Hg anthropogenic emissions
      do img=1,NEMIT_HgA
         ig=IPIG_HgA(img)
         if(ig>0) then
           do k=1,NLAY_HgA
             i04=ip4mem(k,ig,ne)
             call read2d(myid,EmtHgAGas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne),irecHgA(ne),200+ne)
           enddo
         endif
         !!!!!!!!!!! chenhs to adjust Asia Hg0 anthropogenic emission !!!!!!!!!!!
         if(ne.eq.0) then
         IF( img.eq.1 ) THEN  !! Hg0
          do i=sx(ne),ex(ne)
          do j=sy(ne),ey(ne)
             ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             i0f=ip2mem(ne)
             do k=1,NLAY_HgA
               i04e=ip4mem(k,ig,ne)
               if( LATITCRS(i0f+ixy).ge.17..and.LATITCRS(i0f+ixy).le.50. ) then
                 if( LONGICRS(i0f+ixy).ge.75..and.LONGICRS(i0f+ixy).le.135. ) then 
                   EmtHgAGas(i04e+ixy)=EmtHgAGas(i04e+ixy)*1.5
                 endif
               endif
             enddo
           enddo
           enddo
         ENDIF
         endif         
      enddo
    ENDIF !! ifHgA
             close(200+ne) ! juanxiong he

    IF(ifHgN.eq.1) then  !! Hg natural emissions
      do img=1,NEMIT_HgN
         ig=IPIG_HgN(img)
         if(ig>0) then
           i0=ip2memGas(ig,ne)
           call read2d(myid,EmtBBHg(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
               nx(ne),ny(ne),irecHgN(ne),250+ne)
           IF(NLAY_HgN>=2) THEN
           call read2d(myid,EmtOceHg(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
               nx(ne),ny(ne),irecHgN(ne),250+ne)
           ENDIF
           IF(NLAY_HgN>=3) THEN
           call read2d(myid,EmtGeoHg(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
               nx(ne),ny(ne),irecHgN(ne),250+ne)
           ENDIF 
           IF(NLAY_HgN>=4) THEN
           call read2d(myid,EmtReeHg(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
               nx(ne),ny(ne),irecHgN(ne),250+ne)
           ENDIF
         endif
         !!!!!!!!!!! chenhs to adjust Hg Ocean emission !!!!!!!!!!! 
         if(ne.eq.0) then
         IF( img.eq.1 ) THEN  !! Hg0
          do i=sx(ne),ex(ne)
          do j=sy(ne),ey(ne)
             ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
             i0e=ip2memGas(ig,ne)
             i0f=ip2mem(ne)
             EmtOceHg(i0e+ixy)=EmtOceHg(i0e+ixy)*1.667  !! 3000 to 5000
             if( LATITCRS(i0f+ixy).lt.-30. ) then
               EmtOceHg(i0e+ixy)=EmtOceHg(i0e+ixy)*1.6
             endif
             if( LATITCRS(i0f+ixy).gt.-30..and.LATITCRS(i0f+ixy).lt.-10. ) then
               EmtOceHg(i0e+ixy)=EmtOceHg(i0e+ixy)*1.5
             endif
             if( LATITCRS(i0f+ixy).gt.30. ) then
               if( LONGICRS(i0f+ixy).ge.-88..and.LONGICRS(i0f+ixy).le.33. ) then
                 EmtOceHg(i0e+ixy)=EmtOceHg(i0e+ixy)*0.8*1.2
               else 
                 EmtOceHg(i0e+ixy)=EmtOceHg(i0e+ixy)*0.8*0.8
               endif
             endif
             if( LATITCRS(i0f+ixy).gt.10..and.LATITCRS(i0f+ixy).lt.30. ) then
               EmtOceHg(i0e+ixy)=EmtOceHg(i0e+ixy)*0.667
             endif
             if( LATITCRS(i0f+ixy).gt.-10..and.LATITCRS(i0f+ixy).lt.10. ) then
               EmtOceHg(i0e+ixy)=EmtOceHg(i0e+ixy)*0.6
             endif
           enddo
           enddo
         ENDIF
         endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo 
    ENDIF !! ifHgN
           close(250+ne) ! juanxiong he
    ENDIF !! ifHg
    deallocate(IPIG,IPIG_NOx,IPIG_AIRC,IPIG_HgA,IPIG_HgN)

    ! CCCCCCC  TO CONVERT MODIS LANDUSE TO USGS CATAGORY
    IF(IMODIS==1) THEN
       i0= ip2mem(ne)
       CALL MODISLAND(myid,land_USE(i0),sx(ne),ex(ne),sy(ne),ey(ne),ne)
    ENDIF

  endif !! imonthEmit

  !!!!!!!!!! chenhs to adjust meteo !!!!!!!!!!!!!!!!!!!!!!
  do i=sx(ne),ex(ne)
    do j=sy(ne),ey(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        i0=ip2mem(ne)
        if( LATITCRS(i0+ixy).ge.-15..and.LATITCRS(i0+ixy).le.10. ) then
          if(landmask(i0+ixy).eq.0) then
            RAINCON(i0+ixy)=RAINCON(i0+ixy)*0.7
            RAINNON(i0+ixy)=RAINNON(i0+ixy)*0.7
          endif 
        endif 
        if( LATITCRS(i0+ixy).ge.10..and.LATITCRS(i0+ixy).le.70. ) then
           RAINCON(i0+ixy)=RAINCON(i0+ixy)*1.3
           RAINNON(i0+ixy)=RAINNON(i0+ixy)*1.3
        endif
        if( LATITCRS(i0+ixy).ge.-65..and.LATITCRS(i0+ixy).le.-30. ) then
           RAINCON(i0+ixy)=RAINCON(i0+ixy)*1.2
           RAINNON(i0+ixy)=RAINNON(i0+ixy)*1.2
        endif
    enddo
  enddo
  
   !-----------to find the layer of top conditons from global model and of pbl height----
    iwb=sx(ne)-1;ieb=ex(ne)+1
    jsb=sy(ne)-1;jeb=ey(ne)+1
    allocate(ppp(iwb:ieb,jsb:jeb,nzz),ttn(iwb:ieb,jsb:jeb,nzz),tp(iwb:ieb,jsb:jeb))

   do j = sy(ne),ey(ne)
    do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz
       i03=ip3mem(k,ne)
       ppp(i,j,k)=Plev(i03+ixy)
       ttn(i,j,k)=t(i03+ixy)
      enddo  !K
    enddo  !i
   enddo !j   

     plimu    = 45000.
     pliml    = 7500.
     plimlex  = 10000.
     dofill   = .TRUE.
    CALL tropo(myid,ttn, sx(ne),ex(ne),sy(ne),ey(ne),nzz,ppp, plimu, pliml, plimlex,dofill, tp, tperr)
    
   do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
       i0 = ip2mem(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       tropp(i0+ixy)=tp(i,j)/100. !PA-->HPA
       NPBL (i0+ixy) = 1
      do k=1,nzz
         i03 = ip3mem(k,ne)
         i02=ip2mem(ne)
         if(k==nzz) then
          i03_t = ip3mem(k,ne)
         else 
          i03_t = ip3mem(k+1,ne)  
         endif     
      if(Plev(i03+ixy)>=tropp(i0+ixy).and.Plev(i03_t+ixy)<tropp(i0+ixy)) then
          ktop(i02+ixy)=float(k)   
      endif 

        IF( (heiz(i03+ixy)-HGT1(i0+ixy)).LE. PBL_HGT(i0+ixy)) NPBL (i0+ixy) = K ! get the  layer of pbl eight, chenhs change 

      enddo !k

        
      if(i==53.and.j==30.and.ne==1) &
      print*,'ktop=',ktop(i02+ixy),tropp(i0+ixy)
      enddo !i
    enddo !j
  deallocate(ttn,ppp,tp)

!-----------------------------   
	if(numprocs .gt. 1 )then        
	call mpi_barrier( local_com, ierr )
	do k=1,nzz
	i03=ip3mem(k,ne)
	call exchng2( myid, u(i03), sx(ne), ex(ne), sy(ne), ey(ne),  &
		comm2d(ne), stride(ne),  nbrleft(ne), nbrright(ne), &
		nbrtop(ne), nbrbottom(ne) )
	call exchng2( myid, v(i03), sx(ne), ex(ne), sy(ne), ey(ne),  &
		comm2d(ne), stride(ne),  nbrleft(ne), nbrright(ne), &
		nbrtop(ne), nbrbottom(ne) )
	if(imasskeep==1)then
	call exchng2( myid, kpmass_m1(i03), sx(ne), ex(ne), sy(ne), ey(ne),  &
	      comm2d(ne), stride(ne),  nbrleft(ne), nbrright(ne), &
	      nbrtop(ne), nbrbottom(ne) )
	endif
	enddo
	i02=ip2mem(ne)
	call exchng2( myid, TERRAIN(i02), sx(ne), ex(ne), sy(ne), ey(ne),  &
		comm2d(ne), stride(ne),  nbrleft(ne), nbrright(ne), &
		nbrtop(ne), nbrbottom(ne) )
	call exchng2( myid, HGT1(i02), sx(ne), ex(ne), sy(ne), ey(ne),  &
		comm2d(ne), stride(ne),  nbrleft(ne), nbrright(ne), &
		nbrtop(ne), nbrbottom(ne) )  !! added by chenhs, for get_w
	call exchng2( myid, LAND_USE(i02), sx(ne), ex(ne), sy(ne), ey(ne),  &
		comm2d(ne), stride(ne),  nbrleft(ne), nbrright(ne), &
		nbrtop(ne), nbrbottom(ne) )  !! added by chenhs, for PUTSALT
	endif

    ! set boundary winds
    do k=1,nzz   ! get boundary

    i03=ip3mem(k,ne)
    if(ne.eq.1.and.ifglobal.eq.1) then 
    call setwindbound(myid,v(i03),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne))
    else
    call setwindbound(myid,u(i03),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne))
    call setwindbound(myid,v(i03),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne))
    if(k.eq.1) then  !! by chenhs
    i02=ip2mem(ne)
    call setwindbound(myid,HGT1(i02),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne))
    call setwindbound(myid,LAND_USE(i02),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne))
    endif
    endif
   enddo
 
!!++++++++++++++++++ chenhs,polar transport ++++++++++++!!
	if(ifglobal.eq.-1) then  !! when j=1 and j=180, v=0, so no need to do polar transport
	if(ne==1) then
	do ips=1,ipolarnum
	if(ipolarmrk(2,ips)==myid) then  !!! need to send
	do k=1,nzz
	   itsp=1
	   i0=ip3mem(k,ne)
	  call getvalue(myid,u(i0),sx(ne),ex(ne),sy(ne),ey(ne),ipolarmrk(4,ips),ipolarmrk(5,ips),atestsu(k,itsp))
	   itsp=itsp+1
	  call getvalue(myid,v(i0),sx(ne),ex(ne),sy(ne),ey(ne),ipolarmrk(4,ips),ipolarmrk(5,ips),atestsu(k,itsp))
	enddo
	do k=1,nzz
	   ipsmark = ips + 720 + (k-1)*ipolarnum
	   do ibeibei=1,2
	   atestsu0(ibeibei)=atestsu(k,ibeibei)
	  enddo
	call mpi_send(atestsu0,2,mpi_real,ipolarmrk(1,ips),ipsmark,comm2d(ne),ierr)
	enddo
	endif
	
	if(ipolarmrk(1,ips)==myid) then  !!! need to receive
	do k=1,nzz
	  ipsmark = ips + 720 + (k-1)*ipolarnum
	   call mpi_recv(atestru0,2,mpi_real,ipolarmrk(2,ips),ipsmark,  &
		comm2d(ne),status,ierr)
	  do ibeibei=1,2
	   atestru(k,ibeibei)=atestru0(ibeibei)
	  enddo
	enddo
	do k=1,nzz
	  if(ipolarmrk(5,ips)==1.or.ipolarmrk(5,ips)==ny(ne)) then
	     if(ipolarmrk(5,ips)==1) then
	      ipp=ipolarmrk(5,ips)-1
	     else
	      ipp=ipolarmrk(5,ips)+1
	     endif
	     i0=ip3mem(k,ne)
	     call putvalue(myid,u(i0),sx(ne),ex(ne),sy(ne),ey(ne),ipolarmrk(3,ips),ipp,atestru(k,1))
	     call putvalue(myid,v(i0),sx(ne),ex(ne),sy(ne),ey(ne),ipolarmrk(3,ips),ipp,atestru(k,2))
	  endif
	enddo
	endif
	enddo
	endif 
	endif
!!++++++++++++++++++ chenhs,polar transport ++++++++++++!! 
    
  ! to calculate w
  i02=ip2mem(ne)
  do  k=1,nzz
  i0=ip3mem(k,ne)
  if(k.gt.1)then
  i00=ip3mem(k-1,ne)
  else
  i00=ip3mem(1,ne)
  endif
  call get_w(myid,w(i0),w(i00),u(i0),v(i0),dx(i0),dy(i0),dz(i0), &
            HGT1(i02),landmask(i02),latitcrs(i02),hh,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),k,ne,ktop(i02)) !! modify by chenhs
  enddo

  endif   !! to check if this time need to calculate ?
endif !! added by chenhs,to adjust met and emit input frequency
enddo ! nest

if(ifprocess.eq.1) then
  GasTermBal = 0.
  gasOLD = gas
endif

! FOR DUST  emissions , dry and wet deposition rate kg/hr/m2
do ne=1,nest  !! by chenhs 
 if(mod(ihour2,dtout(ne))==0) then  !! added by chenhs,to adjust output frequency
   I02 = IP2MEM(NE)
   call set_default1( myid, DUSTEMISS(I02),DUSTDRY(I02),DUSTWET(I02),DUSTGRAV(I02), sx(ne), ex(ne), sy(ne), ey(ne))
  if(ifdustcom.eq.1) then 
   call set_default1( myid, DUSTDRYSO4(I02),DUSTDRYNO3(I02),DUSTDRYFeII(I02),DUSTDRYFeIII(I02), sx(ne), ex(ne), sy(ne), ey(ne))
   call set_default1( myid, DUSTWETSO4(I02),DUSTWETNO3(I02),DUSTWETFeII(I02),DUSTWETFeIII(I02), sx(ne), ex(ne), sy(ne), ey(ne))
   call set_default1( myid, DUSTGRAVSO4(I02),DUSTGRAVNO3(I02),DUSTGRAVFEII(I02),DUSTGRAVFEIII(I02), sx(ne), ex(ne), sy(ne), ey(ne))
  endif
  ! FOR SEA SALT EMISSIONS kg/m2/hr
  call set_default( myid, SEAEMISS(I02), sx(ne), ex(ne), sy(ne), ey(ne))   
  !!!!!  For wet deposition !!!!
  do ig=1,igas
    I0 = ip2memGas(ig,ne)
    call set_default2( myid, WETDEP(I0), WETDEP2(I0), DRYDEP2(I0), sx(ne), ex(ne), sy(ne), ey(ne))
  enddo
  !!!!!!!!!!!!!!!!!!!!!!! make balance !!!!!!!!!!!!!!!!!!!!!!!!
  if(ifbalance.eq.1.and.ne.eq.1) then
    if(myid.eq.0) then
    acdry=0. ; acwet=0.
    endif
  endif

 endif
enddo  ! nest

!-------------------------------------------------------------------
!         process analysis
!-------------------------------------------------------------------

dt=600.          ! by chenhs, to save time
ndt=3600./dt     ! chenhs, timestep 
do itt=1,ndt       ! 6*600=3600sec, the base unit is 600s

   do ne=1,nest    ! for nest

   !!!!!!!!!!!!!! make balance !!!!!!!!!!!!!!!!!!!!!!
   if(ifbalance.eq.1.and.ne.eq.1) then
   if(myid.eq.0) then
     close(11)
     open(11,file='real_balance.txt',position='append')
     write(11,*) "time step is:",it1,itt
     write(11,*) "+++++++++++++++++++++++++++++++++++++++++++++"
     close(11)
   endif
   endif
 
 if(it1.ge.ntbeg(ne) .and. it1.le.ntend(ne))then
      
 time(ne) = (it1-1)*ndt*dt +itt*dt     ! caculated time
 if(mod(time(ne),dtstep(ne)).eq.0) then  !! added by chenhs, to adjust integration time step

!------------for gas--------ppb-->ug/m3-------------------------------
  do j=sy(ne),ey(ne)
  do i=sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
  do k=1,nzz
      i03 = ip3mem(k,ne)
      do ig=1,igasCBM
        i04 = ip4mem(k,ig,ne)
        gas(i04+ixy) = gas(i04+ixy)*                  &
         GC_MOLWT(ig)/( (0.08206*T(i03+ixy))/(Plev(i03+ixy)/1000.) )
      enddo
  enddo
  enddo
  enddo
!-----------------------emission--------------------
!!!!!!!!!!!!!!!!!!!!!!! make balance !!!!!!!!!!!!!!!!!!!!!!!!
if(ifbalance.eq.1.and.ne.eq.1) then
  do ig=103,105
     amount_bf(myid+1,ig)=0.0
     do k=1,nzz
        I03=IP3MEM(K,NE)
        I04=IP4MEM(K,IG,NE)
        call N_Balance(MYID,gas(I04),amount_bf(myid+1,ig),DX(I03),DY(I03),DZ(I03),SX(NE),EX(NE),SY(NE),EY(NE))
     enddo
  enddo

  CDNUM='c000'
  WRITE(CDNUM(2:4),'(I3.3)')MYID
  OPEN(10,FILE='Amount.bf.dat'//CDNUM,FORM='FORMATTED')
  WRITE(10,*) MYID,amount_bf(myid+1,103),amount_bf(myid+1,104),amount_bf(myid+1,105)
  close(10)

  IF(NUMPROCS>1) CALL mpi_barrier( local_com,IERR )
   total_bf=total_af
   total_af=0.
   DO I=0,NUMPROCS-1
      CDNUM='c000'
      total_tmp=0.
      WRITE(CDNUM(2:4),'(I3.3)')I
      OPEN(10,FILE='Amount.bf.dat'//CDNUM,FORM='FORMATTED')
      read(10,*) iiii,total_tmp(1),total_tmp(2),total_tmp(3)
      CLOSE(10)
      total_af=total_af+sum(total_tmp(1:3))
   ENDDO
  if(myid.eq.0) then
   open(11,file='real_balance.txt',position='append')
   write(11,'(a60,f16.2,2x,f16.6)') "BeforeEmit,Global=",total_af,(total_af-total_bf)/total_bf*100
   close(11)
  endif
endif
 
 ! to calculate
 i02=ip2mem(ne)
 allocate(ratioemitAnt(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz,igas))
 allocate(ratioemitShp(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz,igas))
 allocate(ratioemitBB(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz,igas))
 allocate(ratioemitBio(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz,igas))
 allocate(ratioemitOce(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz,igas))
 allocate(ratioemitSoi(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz,igas))
 allocate(ratioemitLig(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz,igas))
 allocate(ratioemitAirc(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz,igas))

    ratioemitAnt = 0.0
    ratioemitShp = 0.0
    ratioemitBB  = 0.0
    ratioemitBio = 0.0
    ratioemitOce = 0.0
    ratioemitSoi = 0.0
    ratioemitLig = 0.0
    ratioemitAirc = 1.0
               
do iemittype=1,NLAY_EM+NLAY_NOx+2+NLAY_HgN ! 2 for Aircraft and Hg anthropogenic emission  

 do j = sy(ne),ey(ne)
   do i = sx(ne),ex(ne)    
     i02=ip2mem(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do ig =1 ,igas         
      call get_ratio_emit(myid,nzz,ig,LAND_USE(i02+ixy),latitcrs(i02+ixy),ratioemitAnt(i,j,:,ig),&
                          ratioemitBB(i,j,:,ig),ratioemitLig(i,j,:,ig)) 
      ratioemitShp(i,j,1,ig)= 1.0
      ratioemitBio(i,j,1,ig)= 1.0
      ratioemitOce(i,j,1,ig)= 1.0
      ratioemitSoi(i,j,1,ig)= 1.0  !! chenhs, close soil NOx emission 
     enddo ! igas
   enddo !i 
 enddo   !j 

 do ig=1,igas   ! for gas
       i02Gas= ip2memGas(ig,ne)
    do k=1,nzz !! modified by chenhs     
       i03=ip3mem(k,ne)
       i04=ip4mem(k,ig,ne)
    !!!!!!!!!!!!!!!!
    ! for Source Mark
    if(ifsm(ne)==1) then
    letdoit=0
    do idm=1,idmSet
       if(igMark(idm)==ig)letdoit=idm
    enddo
    if(letdoit>0)then
      i02=ip2mem(ne)
     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1 
        tmpMarkCon(i02+ixy)=gas(i04+ixy)
     enddo
     enddo
    endif
    endif !ifsm
    !!!!!!!!!!!!!!!!
     if(iemittype==1)then
       call putemit(myid, gas(i04),EmtAntGas(i02Gas),dz(i03), &
            sx(ne),ex(ne),sy(ne),ey(ne),dtstep(ne),ratioemitAnt,k,nzz,ne,ig,igas)  ! change dt,chenhs
     endif
    if(iemittype==2)then
       call putemit(myid, gas(i04),EmtShpGas(i02Gas),dz(i03), &
            sx(ne),ex(ne),sy(ne),ey(ne),dtstep(ne),ratioemitShp,k,nzz,ne,ig,igas)  ! change dt,chenhs
    endif
    if(iemittype==3)then
       call putemit(myid, gas(i04),EmtBBGas(i02Gas),dz(i03), &
            sx(ne),ex(ne),sy(ne),ey(ne),dtstep(ne),ratioemitBB,k,nzz,ne,ig,igas)  ! change dt,chenhs
    endif
    if(iemittype==4)then
       call putemit(myid, gas(i04),EmtBioGas(i02Gas),dz(i03), &
            sx(ne),ex(ne),sy(ne),ey(ne),dtstep(ne),ratioemitBio,k,nzz,ne,ig,igas) ! change dt,chenhs
    endif
    if(iemittype==5)then
       call putemit(myid, gas(i04),EmtOceGas(i02Gas),dz(i03), &
            sx(ne),ex(ne),sy(ne),ey(ne),dtstep(ne),ratioemitOce,k,nzz,ne,ig,igas) ! change dt,chenhs
    endif
    if(iemittype==6)then
       call putemit(myid, gas(i04),EmtSoiGas(i02Gas),dz(i03), &
            sx(ne),ex(ne),sy(ne),ey(ne),dtstep(ne),ratioemitSoi,k,nzz,ne,ig,igas) ! change dt,chenhs
    endif
    if(iemittype==7)then
       call putemit(myid, gas(i04),EmtLigGas(i02Gas),dz(i03), &
            sx(ne),ex(ne),sy(ne),ey(ne),dtstep(ne),ratioemitLig,k,nzz,ne,ig,igas) ! change dt,chenhs
    endif
    if(iemittype==8)then
       call putemit(myid, gas(i04),EmtAircGas(i04),dz(i03), &
            sx(ne),ex(ne),sy(ne),ey(ne),dtstep(ne),ratioemitAirc,k,nzz,ne,ig,igas) ! change dt,chenhs
    endif
    if(iemittype==9)then
       call putemit(myid, gas(i04),EmtHgAGas(i04),dz(i03), &
            sx(ne),ex(ne),sy(ne),ey(ne),dtstep(ne),ratioemitAirc,k,nzz,ne,ig,igas) ! change dt,chenhs
    endif
    if(iemittype==10)then
       call putemit(myid, gas(i04),EmtBBHg(i02Gas),dz(i03), &
            sx(ne),ex(ne),sy(ne),ey(ne),dtstep(ne),ratioemitBB,k,nzz,ne,ig,igas)  ! change dt,chenhs
    endif
    if(iemittype==11)then
       call putemit(myid, gas(i04),EmtOceHg(i02Gas),dz(i03), &
            sx(ne),ex(ne),sy(ne),ey(ne),dtstep(ne),ratioemitOce,k,nzz,ne,ig,igas)  ! change dt,chenhs
    endif
    if(iemittype==12)then
       call putemit(myid, gas(i04),EmtGeoHg(i02Gas),dz(i03), &
            sx(ne),ex(ne),sy(ne),ey(ne),dtstep(ne),ratioemitOce,k,nzz,ne,ig,igas)  ! change dt,chenhs
    endif
    if(iemittype==13)then
       call putemit(myid, gas(i04),EmtReeHg(i02Gas),dz(i03), &
            sx(ne),ex(ne),sy(ne),ey(ne),dtstep(ne),ratioemitOce,k,nzz,ne,ig,igas)  ! change dt,chenhs
    endif

    ! for Source Mark
    if(ifsm(ne)==1)then 
    if(letdoit>0)then
      i02=ip2mem(ne)
      i04=ip4mem(k,ig,ne)
     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1 
        OrgCon = tmpMarkCon(i02+ixy)
        DeltSpeMark=gas(i04+ixy)-OrgCon
        if(DeltSpeMark > 0 )then 
             do ism=1,ismMax
               i04sm=ipSMmem(k,ism,letdoit,ne)
               TmpSM(ism)=SourceMark(i04sm+ixy)
             enddo

        MapS=int(MapSource(i02+ixy))

         IHgtLev=iemittype-1
        call GetSMarkChange(DeltSpeMark,OrgCon,TmpSM,MapS, &
                  ISrcDefined,IHgtLev,ismMax)

             do ism=1,ismMax
               i04sm=ipSMmem(k,ism,letdoit,ne)
               SourceMark(i04sm+ixy)=TmpSM(ism)
             enddo

        endif
     enddo
     enddo
    endif
    endif !ifsm

    enddo
 enddo

enddo  ! iemittype Zifa 2007/03/03

! ****   PUT DUST AND SEA SALT EMISSIONS 
! ** SEA SALT EMISSIONS ***
  DO IS = 1, ISIZE 
    I03AER = IP3MEMAER(IS,1,NE)
    I03=IP3MEM(1,NE)
    I05=IP5MEM(1,IS,1,NE)    
    I02    = IP2MEM(NE)
    if(ifseacom.eq.1) then
    I05C1=ip5memcs(1,IS,1,NE)
    I05C2=ip5memcs(1,IS,2,NE)
    I05C3=ip5memcs(1,IS,3,NE)
    I05C4=ip5memcs(1,IS,4,NE)
    I05C5=ip5memcs(1,IS,5,NE)
    I05C6=ip5memcs(1,IS,6,NE)
    I05C7=ip5memcs(1,IS,7,NE)
    I05C8=ip5memcs(1,IS,8,NE)
   
   CALL PUTSALTCOM (MYID, AER(I05), SEACOMP(I05C1),SEACOMP(I05C2),SEACOMP(I05C3),SEACOMP(I05C4),&
                SEACOMP(I05C5),SEACOMP(I05C6),SEACOMP(I05C7),SEACOMP(I05C8),MSIZDIS(IS), &
                MSIZDIS(IS+1),LAND_USE(I02),FICE(I02), &
                U10(I02), V10(I02), RHSFC(I02),DZ(I03),DX(I03),SEAEMISS(I02),SX(NE),EX(NE), SY(NE), EY(NE), NE, dtstep(ne)) ! change dt,chenhs
    else
   CALL PUTSALT (MYID,AER(I05),MSIZDIS(IS),MSIZDIS(IS+1),LAND_USE(I02),FICE(I02), &
                U10(I02), V10(I02), RHSFC(I02),DZ(I03),DX(I03),SEAEMISS(I02),SX(NE),EX(NE), SY(NE), EY(NE), NE, dtstep(ne)) ! change dt,chenhs
    endif
  ENDDO ! IS

 ! ----- TO CALCULATE THE HEIGHT FACTOR FOR DUST EMNISSIONS
 ! TO CALCULATE THE TOTAL EMISSIONS FACTORS(TOTALDUST) BY SCICHINA D, 2011(4):234-242  
  DO K = 1, KDUSTTOP
    I02    = IP2MEM(NE)
    I03    = IP3MEM(K,NE)
   CALL DUSTHEIGHT  (MYID,HEIZ(I03), HGT1(I02),LAND_USE(I02),DUSTK(I03),TOTALDUST(I02),U(I03),V(I03),SX(NE), EX(NE), SY(NE), EY(NE),K,NE ) 
  ENDDO

 ! TO CALCULATE THE WEIGHTFACT BY DUSTK/TOTALDUST 
  DO K = 1, KDUSTTOP
   I02    = IP2MEM(NE)
   I03    = IP3MEM(K,NE)
   CALL DUSTHGTFACT (MYID,TOTALDUST(I02),DUSTK(I03),DUSTHGTF(I03),SX(NE),EX(NE), SY(NE), EY(NE),NE )
  ENDDO

  DO IS = 1, ISIZE
   I03AER = IP3MEMAER (IS, 2, NE)
   I03    = IP3MEM(1,NE)
   I02    = IP2MEM(NE) 
  CALL GETZ0(MYID,LAND_USE(I02),LATITCRS(I02),Z0(I02),IMONTH2,SX(NE), EX(NE), SY(NE), EY(NE),NE)
  CALL GETUST0(MYID,LAND_USE(I02), FSOIL(I02), FVEG(I02),UST0(I02),SOILT(I03),SOILRH(I03), RHSFC(I02),LONGICRS(I02),SX(NE),  EX(NE), SY(NE), EY(NE), NE)

   DO K = 1, KDUSTTOP
    I03_2 =  IP3MEM(K,NE)
    I05    = IP5MEM(K,IS,2,NE)
    if(ifdustcom.eq.1) then
    I05C1=ip5memc(K,IS,1,NE)
    I05C2=ip5memc(K,IS,2,NE)
    I05C3=ip5memc(K,IS,3,NE)
    I05C4=ip5memc(K,IS,4,NE)
    I05C5=ip5memc(K,IS,5,NE)
    I05C6=ip5memc(K,IS,6,NE)
    I05C7=ip5memc(K,IS,7,NE)
    I05C8=ip5memc(K,IS,8,NE)
    I05C9=ip5memc(K,IS,9,NE)
    I05C10=ip5memc(K,IS,10,NE)

    CALL PUTDUSTCOM (MYID,AER(I05),DUSTCOMP(I05C1),DUSTCOMP(I05C2),DUSTCOMP(I05C3),DUSTCOMP(I05C4),&
                  DUSTCOMP(I05C5),DUSTCOMP(I05C6),DUSTCOMP(I05C7),DUSTCOMP(I05C8),DUSTCOMP(I05C9),&
                  DUSTCOMP(I05C10), DUSTHGTF(I03_2), FICE(I02), FSNOW(I02), LAND_USE(I02), FSOIL(I02), &
                  FVEG(I02), Z0(I02),UST(I02), UST0(I02), SOILT(I03),T2(I02),SOILRH(I03), RHSFC(I02),&
                  U10(I02),V10(I02), EMITFACT(I02),DZ(I03_2),DUSTEMISS(I02),SX(NE),&
                 EX(NE), SY(NE), EY(NE), dtstep(ne), SFT(IS),K,NE ,IS, ITT)  !! change dt,chenhs
    else
    CALL PUTDUST (MYID,AER(I05),DUSTHGTF(I03_2), FICE(I02), FSNOW(I02), LAND_USE(I02), FSOIL(I02), &
                  FVEG(I02), Z0(I02),UST(I02), UST0(I02), SOILT(I03),T2(I02),SOILRH(I03), RHSFC(I02),&
                  U10(I02),V10(I02), EMITFACT(I02),DZ(I03_2),DUSTEMISS(I02),SX(NE),&
                 EX(NE), SY(NE), EY(NE), dtstep(ne), SFT(IS),K,NE ,IS, ITT)  !! change dt,chenhs   
    endif
   ENDDO ! K

  ENDDO ! IS

!Revised by Zifa
!---------for gas---------ug/m3-->ppb-------------------------------
 do j=sy(ne),ey(ne)
 do i=sx(ne),ex(ne)
    ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
 do k=1,nzz
      i03 = ip3mem(k,ne)
      do ig=1,igasCBM
        i04 = ip4mem(k,ig,ne)
        gas(i04+ixy) = gas(i04+ixy)*                          &
        (0.08206*T(i03+ixy))/(Plev(i03+ixy)/1000.)/GC_MOLWT(ig)
                              ! ug/m3 --> ppb
      enddo
  enddo
  enddo
  enddo

  deallocate(ratioemitAnt,ratioemitShp,ratioemitBB,ratioemitBio,&
              ratioemitOce,ratioemitSoi,ratioemitLig,ratioemitAirc)

if(ifglobal.ne.1) then  !  modified by chenhs, to use online global domain
  if(ne == 1 ) then    
!-----------------for global ---------------------------------------
   do k=1,nzz   ! get boundary

    i03=ip3mem(k,ne)
    call setwindbound(myid,globalo3(i03),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne))
    call setwindbound(myid,globalco(i03),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne))
    call setwindbound(myid,globalno2(i03),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne))

   enddo

  IF(iglobal==1) THEN  !from global model
    do k=1,nzz   ! get boundary
      do ig=11,11 ! ozone
       i04=ip4mem(k,ig,ne)
       i0=ip3mem(k,ne)
       call getboundnorth1(myid,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne),ig,globalo3(i0))
       call getboundsouth(myid,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne),ig,globalo3(i0))
       call getboundeast(myid,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                        nx(ne),ny(ne), ig,globalo3(i0))
       call getboundwest(myid,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne),ig,globalo3(i0))
      enddo    !ig
      do ig=17,17 ! CO
       i04=ip4mem(k,ig,ne)
       i0=ip3mem(k,ne)
       call getboundnorth1(myid,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne),ig,globalco(i0))
       call getboundsouth(myid,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne),ig,globalco(i0))
       call getboundeast(myid,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                        nx(ne),ny(ne), ig,globalco(i0))
       call getboundwest(myid,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne),ig,globalco(i0))
      enddo    !ig

      do ig=6,6 ! NO2
       i04=ip4mem(k,ig,ne)
       i0=ip3mem(k,ne)
       call getboundnorth1(myid,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne),ig,globalno2(i0))
       call getboundsouth(myid,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne),ig,globalno2(i0))
       call getboundeast(myid,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                        nx(ne),ny(ne), ig,globalno2(i0))
       call getboundwest(myid,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne),ig,globalno2(i0))
      enddo    !ig

   enddo!k
  ELSE  IF(iglobal==2) then  !fix
!-------------------------------------------------------------------

   do k=1,nzz   ! get boundary 
     do ig=11,11 ! ozone
       i04=ip4mem(k,ig,ne)
       call setboundnorth(myid,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne),ig,45.+float(k)*1.5)
       call setboundsouth(myid,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne),ig,20.+float(k)*0.8)
       call setboundeast(myid,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                        nx(ne),ny(ne),ig,25.+float(k)*0.8)
       call setboundwest(myid,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne),ig,40.+float(k)*1.5,k)
      enddo    !ig
      do ig=17,17 ! co
       i04=ip4mem(k,ig,ne)
       call setboundnorth(myid,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne),ig,150.-float(k)*5.)
       call setboundsouth(myid,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne),ig,100.-float(k)*5.)
       call setboundeast(myid,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                        nx(ne),ny(ne),ig,200.-float(k)*10.)
       call setboundwest(myid,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                         nx(ne),ny(ne),ig,100.-float(k)*5.,k)
      enddo !ig co

   enddo   !k
  ENDIF ! iglobal
 endif ! ne 
endif ! ifglobal
!============================set top boundary by li 05-04-21 ================== need to modify,chenhs
!--------------------------------global ------------------------------------
  IF(iglobal==1) THEN
   if(ne>=1) then ! chenhs 
      do k=1,nzz  ! chenhs
       do ig=11,11 !Ozone
         i04=ip4mem(k,ig,ne)
         i03=ip3mem(k,ne)
         i0=ip3mem(k,ne)
         i02=ip2mem(ne) 
         call setboundtop(myid,ig,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                           nx(ne),ny(ne),globalo3(i0),ne,k,ktop(i02))
       enddo      
      enddo

      do k=1,nzz
       do ig=17,17 !CO
         i04=ip4mem(k,ig,ne)
         i03=ip3mem(k,ne)
         i0=ip3mem(k,ne)
         i02=ip2mem(ne)
         call setboundtop(myid,ig,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                           nx(ne),ny(ne),globalco(i0),ne,k,ktop(i02))
       enddo
      enddo
     do k=1,nzz
       do ig=6,6 !no2
         i04=ip4mem(k,ig,ne)
         i03=ip3mem(k,ne)
         i0=ip3mem(k,ne)
         i02=ip2mem(ne)
         call setboundtop(myid,ig,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
                           nx(ne),ny(ne),globalno2(i0),ne,k,ktop(i02))
       enddo
      enddo

    endif
  ELSE
!-------------------------------finish--------------------------------------
    topo3=0. !! by chenhs 
    if(ne==1.or.ne==2.or.ne==3.or.ne==4) then
      do k=nzz,nzz
       do ig=11,11
         i04=ip4mem(k,ig,ne)
         i03=ip3mem(k,ne)
         i0=ip2mem(ne)
         i02=ip2mem(ne)
         call setboundtop1(myid,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &   !!! topo3 does not have value, need to modify, chenhs
                           nx(ne),ny(ne),topo3(i0),ne,k,ktop(i02))

       enddo
      enddo
    endif
ENDIF
!==============================================================================
 if(imasskeep==1)then
   do k=1,nzz   ! get boundary 
    i03=ip3mem(k,ne)
    call setboundnorthkpm(myid,kpmass_m1(i03),sx(ne),ex(ne), &
                         sy(ne),ey(ne), nx(ne),ny(ne),1000.)
    call setboundsouthkpm(myid,kpmass_m1(i03),sx(ne),ex(ne), &
                         sy(ne),ey(ne), nx(ne),ny(ne),1000.)
    call setboundeastkpm(myid,kpmass_m1(i03),sx(ne),ex(ne), &
                         sy(ne),ey(ne), nx(ne),ny(ne),1000.)
    call setboundwestkpm(myid,kpmass_m1(i03),sx(ne),ex(ne), &
                         sy(ne),ey(ne), nx(ne),ny(ne),1000.)
    enddo
  endif

!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
if(ifprocess.eq.1) then
   ! iprocess
   do k=1,nzz
       do ig=1,iPrintTermGas    ! for gas phase
           igg=IGGPOS(ig)
           i04=ip4mem(k,igg,ne)
           i05=ipGasTermBal(k,1,ig,ne) ! 1, emit
          call termbal(myid,gasOLD(i04),gas(i04),GasTermBal(i05), &
                k,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dtstep(ne))  !! change dt,chenhs
       enddo
    enddo
endif
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
!=====================================================================
!====    to write boundary conditions into next Domain           =====
!=====================================================================

time(ne) = (it1-1)*ndt*dt +itt*dt     ! change by chenhs
iitime = time(ne)
IF(NE.LE.NEST-1 .and. mod(iitime,int(dtstep(1))) == 0)THEN  !! change by chenhs,when all domain in the same itt then write boundary conditions
                                                       !! may need to modify when dtstep(ne) is changed

  ! following points need to be sent to other CPU
  !                  need to be got from other CPU
  DO IPS=1,NPSR

  if(IISCPU(IPS)==myid .and. ne==IsNest(IPS))then   ! need to send
        i=IsLocX(IPS)
        j=IsLocY(IPS)
     do k=1,nzz
        itsp=0
        ! for gas speceis
        do ig=1,igas
          itsp=itsp+1
          i04 = ip4mem(k,ig,ne)
     call getvalue(myid,gas(i04),sx(ne),ex(ne),sy(ne),ey(ne), &
              i,j,atestS(k,itsp))
        enddo
        ! for aerosols
        do ia=1,iaer
        do is=1,isize
           i05= ip5mem(k,is,ia,ne)
           itsp=itsp+1
     call getvalue(myid,aer(i05),sx(ne),ex(ne),sy(ne),ey(ne), &
              i,j,atestS(k,itsp))
        enddo
        enddo
        ! for u,v,by chenhs
        i03= ip3mem(k,ne)
        itsp=itsp+1
     call getvalue(myid,u(i03),sx(ne),ex(ne),sy(ne),ey(ne), &
              i,j,atestS(k,itsp))
        itsp=itsp+1
     call getvalue(myid,v(i03),sx(ne),ex(ne),sy(ne),ey(ne), &
              i,j,atestS(k,itsp))
         ! for Source Mark
       if(ifsmt>0)then    ! checking 
        do idm=1,idmSet
        do ism=1,ismMax
            i0=ipSMmem(k,ism,idm,ne)
            itsp=itsp+1
     call getvalue(myid,SourceMark(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
                   i,j,atestS(k,itsp))
        enddo
        enddo
       endif

        ! for dust and sea salt composition
       if(ifseacom.eq.1) then 
        do iduc = 1, nseacom
        do is =1 , isize
          i05c=ip5memcs(k,is,iduc,ne)
          itsp = itsp + 1
     call getvalue(myid, seacomp(i05c),sx(ne),ex(ne),sy(ne),ey(ne), &
             i,j,atestS(k,itsp))
        enddo
        enddo
       endif

       if(ifdustcom.eq.1) then
        do iduc = 1, ndustcom
        do is =1 , isize
          i05c=ip5memc(k,is,iduc,ne)
          itsp = itsp + 1
     call getvalue(myid, dustcomp(i05c) ,sx(ne),ex(ne),sy(ne),ey(ne), &
             i,j,atestS(k,itsp))
        enddo
        enddo
       endif
     enddo !k,1,nzz

  do k=1,nzz
     IPSMARK = IPS + 4000 + (k-1) * NPSR
 
  do ibeibei=1,isrnum
     atestS0(ibeibei)=atestS(k,ibeibei)
  enddo
     call MPI_Send(atestS0,isrnum,MPI_REAL,IIRCPU(IPS),IPSMARK,comm2d(ne),ierr)
  enddo 
  endif

  if(IIRCPU(IPS)==myid .and. ne == IsNest(IPS) )then   ! need to receive
    iichild=IrNest(IPS)  ! added by chenhs
    do k=1,nzz
    IPSMARK = IPS + 4000 + (k-1) * NPSR
     call MPI_Recv(atestR0,isrnum,MPI_REAL,IISCPU(IPS),IPSMARK,  &
                comm2d(ne),status,ierr)
      do ibeibei=1,isrnum
         atestR(k,ibeibei)=atestR0(ibeibei)
      enddo
    enddo
     i=IrLocX(IPS)
     j=IrLocY(IPS)
     do k=1,nzz

        ips1 = 0

        if(j==0 .or. j==(ny(iichild)+1))then
        do ip=i,i
           do ig=1,igas
              ips1 = ips1 + 1
              i04 = ip4mem(k,ig,iichild)
     call putvalue(myid,gas(i04),sx(iichild),ex(iichild),sy(iichild),ey(iichild), &
              ip,j,atestR(k,ips1))
           enddo
           do ia=1,iaer
           do is=1,isize
              ips1 = ips1 + 1 
              i05= ip5mem(k,is,ia,iichild)
     call putvalue(myid,aer(i05),sx(iichild),ex(iichild),sy(iichild),ey(iichild), &
              ip,j,atestR(k,ips1))
              enddo
              enddo
           !!! for u,v,by chenhs
           ips1 = ips1 + 1
           i03 = ip3mem(k,iichild)
     call putvalue(myid,u(i03),sx(iichild),ex(iichild),sy(iichild),ey(iichild), &
              ip,j,atestR(k,ips1))
           ips1 = ips1 + 1
     call putvalue(myid,v(i03),sx(iichild),ex(iichild),sy(iichild),ey(iichild), &
              ip,j,atestR(k,ips1))

           ! for Source Mark
        if(ifsmt>0 )then   !checking
         do idm=1,idmSet
         do ism=1,ismMax
            ips1 = ips1 + 1 
            i0=ipSMmem(k,ism,idm,iichild)
     call putvalue(myid,SourceMark(i0),sx(iichild),ex(iichild),sy(iichild),ey(iichild), &
              ip,j,atestR(k,ips1))
         enddo
         enddo
        endif

! for sea salt and dust compositions
         if(ifseacom.eq.1) then
          do iduc = 1, nseacom
          do is = 1, isize
            ips1 = ips1 + 1 
            i05c  = ip5memcs(k,is,iduc,iichild)
     call putvalue(myid,seacomp(i05c),sx(iichild),ex(iichild),sy(iichild),ey(iichild), &
              ip,j,atestR(k,ips1))        
          enddo
          enddo
         endif

         if(ifdustcom.eq.1) then
          do iduc = 1, ndustcom
          do is = 1, isize
            ips1 = ips1 + 1
            i05c  = ip5memc(k,is,iduc,iichild)
     call putvalue(myid,dustcomp(i05c),sx(iichild),ex(iichild),sy(iichild),ey(iichild), &
              ip,j,atestR(k,ips1))
          enddo
          enddo
         endif

        enddo
        endif

        ips1 = 0 
        if(i==0 .or. i==(nx(iichild)+1))then
        do jp=j,j
           do ig=1,igas
              ips1 = ips1 + 1
              i04 = ip4mem(k,ig,iichild)
     call putvalue(myid,gas(i04),sx(iichild),ex(iichild),sy(iichild),ey(iichild), &
              i,jp,atestR(k,ips1))
           enddo

           do ia=1,iaer
           do is=1,isize
              ips1 = ips1 + 1
              i05= ip5mem(k,is,ia,iichild)
     call putvalue(myid,aer(i05),sx(iichild),ex(iichild),sy(iichild),ey(iichild), &
              i,jp,atestR(k,ips1))
           enddo
           enddo
        !!!! for u,v by chenhs
        ips1 = ips1 + 1
        i03 = ip3mem(k,iichild)
     call putvalue(myid,u(i03),sx(iichild),ex(iichild),sy(iichild),ey(iichild), &
              i,jp,atestR(k,ips1))
        ips1 = ips1 + 1
     call putvalue(myid,v(i03),sx(iichild),ex(iichild),sy(iichild),ey(iichild), &
              i,jp,atestR(k,ips1))

        ! for Source Mark
        if(ifsm(ne)==1)then
         do idm=1,idmSet
         do ism=1,ismMax
            ips1 = ips1 + 1
            i0=ipSMmem(k,ism,idm,iichild)
     call putvalue(myid,SourceMark(i0),sx(iichild),ex(iichild),sy(iichild),ey(iichild), &
              i,jp,atestR(k,ips1))
         enddo
         enddo
        endif

! for sea salt and dust compositions
        if(ifseacom.eq.1) then
         do iduc = 1, nseacom
         do is = 1, isize
           ips1 =  ips1 + 1
            i05c= ip5memcs(k,is,iduc,iichild) 
     call putvalue(myid,seacomp(i05c),sx(iichild),ex(iichild),sy(iichild),ey(iichild), &
              i,jp,atestR(k,ips1))
         enddo
         enddo
        endif

        if(ifdustcom.eq.1) then
         do iduc = 1, ndustcom
         do is = 1, isize
           ips1 =  ips1 + 1
            i05c= ip5memc(k,is,iduc,iichild)
     call putvalue(myid,dustcomp(i05c),sx(iichild),ex(iichild),sy(iichild),ey(iichild), &
              i,jp,atestR(k,ips1))
         enddo
         enddo
        endif

        enddo
        endif

    enddo

  endif

  ENDDO

ENDIF

!!!!!!!!!!!!!!!!!!!!!!! make balance !!!!!!!!!!!!!!!!!!!!!!!!
if(ifbalance.eq.1.and.ne.eq.1) then
  do ig=103,105
     amount_bf(myid+1,ig)=0.0
     do k=1,nzz
        I03=IP3MEM(K,NE)
        I04=IP4MEM(K,IG,NE)
        call N_Balance(MYID,gas(I04),amount_bf(myid+1,ig),DX(I03),DY(I03),DZ(I03),SX(NE),EX(NE),SY(NE),EY(NE))
     enddo
  enddo

  CDNUM='c000'
  WRITE(CDNUM(2:4),'(I3.3)')MYID
  OPEN(10,FILE='Amount.bf.dat'//CDNUM,FORM='FORMATTED')
  WRITE(10,*) MYID,amount_bf(myid+1,103),amount_bf(myid+1,104),amount_bf(myid+1,105)
  close(10)

  CALL mpi_barrier( local_com,ierr )
   total_bf=total_af
   total_af=0.
   DO I=0,NUMPROCS-1
      CDNUM='c000'
      total_tmp=0.
      WRITE(CDNUM(2:4),'(I3.3)')I
      OPEN(10,FILE='Amount.bf.dat'//CDNUM,FORM='FORMATTED')
      read(10,*) iiii,total_tmp(1),total_tmp(2),total_tmp(3)
      CLOSE(10)
      total_af=total_af+sum(total_tmp(1:3))
   ENDDO
  if(myid.eq.0) then
   open(11,file='real_balance.txt',position='append')
   write(11,'(a60,f16.2,2x,f16.6)') "Before ADV_HORI,Global=",total_af,(total_af-total_bf)/total_bf*100
   close(11)
  endif
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

time(ne) = (it1-1)*ndt*dt +itt*dt     ! caculated time

!--------------------------advection & Diffution-------------
  kh=450.
!=========================================================================
 if(imasskeep==1)then
  kpmass_m1= 1000.  ! Zifa 2004/09/02
  i02=ip2mem(ne)
  do k=1,nzz-1   ! advection
    i03=ip3mem(k,ne)
    ! for mass keeping Zifa 2004/09/02
    call  adv_hori(myid,kpmass_m1(i03),u(i03),v(i03),dx(i03),dy(i03),&
            sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dtstep(ne))  !! change dt,chenhs
    ! to get mass conservation ratio/error
    call  GetMassRatio(myid,kpmass_m1(i03),RatioMass(i03), &
                 sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
  enddo 
  endif

  ! horizontal advection
    do ig=1,igas    ! for gas phase

      if(imasskeep==1)then
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)
         ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         do k=1,nzz-1
            i03=ip3mem(k,ne)
            i04=ip4mem(k,ig,ne)
            kpmass_m2(i03+ixy)=gas(i04+ixy)
         enddo
      enddo
      enddo
      endif
                                         
    do k=1,nzz-1
       i03=ip3mem(k,ne)
       i04=ip4mem(k,ig,ne)
       i02=ip2mem(ne)
    call  adv_hori(myid,gas(i04),u(i03),v(i03),dx(i03),dy(i03),&
          sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dtstep(ne))  !! change dt,chenhs
    enddo
   
!!!!!!!!!!!!!!!!!!!!!!!!
! for Source Mark
   if(ifsm(ne).eq.1) then
    letdoit=0
    do idm=1,idmSet
       if(igMark(idm)==ig)letdoit=idm
    enddo
    if(letdoit>0)then
    allocate(sm(ismMax,sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1))
    do k=1,nzz-1
       i03=ip3mem(k,ne)
       i04=ip4mem(k,ig,ne)
       i02=ip2mem(ne)
       do ism=1,ismMax
          i04sm=ipSMmem(k,ism,letdoit,ne)
              do i=sx(ne)-1,ex(ne)+1
              do j=sy(ne)-1,ey(ne)+1
               ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1     
               sm(ism,i,j)= SourceMark(i04sm+ixy)
              enddo
              enddo
       enddo
     call adv_hori_mark(myid,gas(i04),u(i03),v(i03),&
            dx(i03),dy(i03),sx(ne), ex(ne), sy(ne),&
            ey(ne),ne,nx(ne),ny(ne),dtstep(ne),k,ktop(i02),ISMMAX,SM) !! change dt,chenhs
       do ism=1,ismMax
          i04sm=ipSMmem(k,ism,letdoit,ne)
              do i=sx(ne)-1,ex(ne)+1
              do j=sy(ne)-1,ey(ne)+1
               ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1     
               SourceMark(i04sm+ixy)=sm(ism,i,j)
              enddo
              enddo
       enddo
    enddo

    deallocate(sm)
    endif 
   endif !ifsm
!!!!!!!!!!!!!!!!!!!!!!!!
   
    if(imasskeep==1)then
    do k=1,nzz-1   !mass conservation 
       i03=ip3mem(k,ne)
       i04=ip4mem(k,ig,ne)        
       call  balance(myid,gas(i04),kpmass_m2(i03),RatioMass(i03), &
                 sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
    enddo
    endif

    enddo    !ig

   ! for DUST AND SEA SALT  aerosols
    do ia=1,iaer
    do is=1,isize
     
      if(imasskeep==1)then
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)
         ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         do k=1,nzz-1
            i03=ip3mem(k,ne)
            i05=ip5mem(k,is,ia,ne)
            kpmass_m2(i03+ixy)=aer(i05+ixy)
         enddo
      enddo
      enddo
      endif

      do k=1,nzz-1
         i05=ip5mem(k,is,ia,ne)
         i03=ip3mem(k,ne)
         CALL ADV_HORI(MYID, AER(I05),u(i03),v(i03), dx(i03),dy(i03),sx(ne), ex(ne),&
                       sy(ne),ey(ne),nx(ne),ny(ne),dtstep(ne)) !! change dt,chenhs
      enddo ! k
       
      do k = 1, nzz - 1

        IF(ia==1) THEN ! sea salt
         if(ifseacom.eq.1) then
         do iduc = 1, nseacom
            i05c   = ip5memcs (k,is,iduc,ne)
            i03=ip3mem(k,ne)
          call ADV_HORI(MYID, SEACOMP(I05C),u(i03),v(i03), dx(i03),dy(i03),sx(ne),ex(ne),&
                       sy(ne),ey(ne),nx(ne),ny(ne),dtstep(ne)) !! change dt,chenhs
         enddo
         endif
        ELSE  IF(ia==2) THEN ! for dust
         if(ifdustcom.eq.1) then
         do iduc = 1, ndustcom
          i05c = ip5memc (k,is,iduc,ne)
          i03=ip3mem(k,ne)
          call ADV_HORI(MYID, DUSTCOMP(I05C),u(i03),v(i03), dx(i03),dy(i03),sx(ne),ex(ne),&
                       sy(ne),ey(ne),nx(ne),ny(ne),dtstep(ne)) !! change dt,chenhs
         enddo ! iduc   
         endif
        ENDIF ! ia

      enddo ! k

       if(imasskeep==1)then
       do k=1,nzz-1   !mass conservation 
          i03=ip3mem(k,ne)
          i05=ip5mem(k,is,ia,ne)
        call  balance(myid,aer(i05),kpmass_m2(i03),RatioMass(i03), &
                 sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
       enddo
       endif

     enddo        !isize
    enddo        !iaer

! advection in horizontal
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
if(ifprocess.eq.1) then
   ! iprocess
     do k=1,nzz
         do ig=1,iPrintTermGas    ! for gas phase
    igg=IGGPOS(ig)
    i03=ip3mem(k,ne)
    i04=ip4mem(k,igg,ne)
    i05=ipGasTermBal(k,2,ig,ne)    !2--> adv hor
    i059= ipGasTermBal(k,9,ig,ne)  !10-->north input
    i0510=ipGasTermBal(k,10,ig,ne) !11-->south input
    i0511=ipGasTermBal(k,11,ig,ne) !12-->west input
    i0512=ipGasTermBal(k,12,ig,ne) !13-->east input
    i0513=ipGasTermBal(k,13,ig,ne) !14-->west output
    i0514=ipGasTermBal(k,14,ig,ne) !15-->east output
    i0515=ipGasTermBal(k,15,ig,ne) !16-->south output
    i0516=ipGasTermBal(k,16,ig,ne) !17-->north output
    call termballi(myid,gasOLD(i04),gas(i04),GasTermBal(i059),GasTermBal(i0510),&
                  GasTermBal(i0511),GasTermBal(i0512),GasTermBal(i0513),&
                  GasTermBal(i0514),GasTermBal(i0515),&
                  GasTermBal(i0516),u(i03),v(i03),k,dx(i03),sx(ne),ex(ne),sy(ne),ey(ne),&
                  nx(ne),ny(ne),dtstep(ne))  !! change dt,chenhs
    call termbal(myid,gasOLD(i04),gas(i04),GasTermBal(i05),&
                  k,sx(ne),ex(ne),sy(ne), ey(ne),nx(ne),ny(ne),dtstep(ne)) !! change dt,chenhs
           enddo
     enddo
endif
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++

!-------------------------------------------------------------------
!!!!!!!!!!!!!!!! vertical advection !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-------------------------------------------------------------------
 iwb=sx(ne)-1;ieb=ex(ne)+1
 jsb=sy(ne)-1;jeb=ey(ne)+1
 allocate(ww(iwb:ieb,jsb:jeb,nzz),conc(iwb:ieb,jsb:jeb,nzz),&
          uu(iwb:ieb,jsb:jeb,nzz),vv(iwb:ieb,jsb:jeb,nzz),&
          ddx(iwb:ieb,jsb:jeb,nzz),ddy(iwb:ieb,jsb:jeb,nzz),ddz(iwb:ieb,jsb:jeb,nzz),&
          kktop(iwb:ieb,jsb:jeb),concmark(iwb:ieb,jsb:jeb,nzz),kp_tmp(iwb:ieb,jsb:jeb,nzz))
! for Gas
!--to get the step for stable-----
 nstep=1  !! by chenhs
 do j = sy(ne),ey(ne)
   do i = sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
     i02=ip3mem(k,ne)
      nstep0(i,j)=1.+dtstep(ne)/min(4000., dz(i02+ixy)/abs(1.E-09+w(i02+ixy))) !! change by chenhs to make larger time step
                                                                              !! add abs by chenhs
      if(nstep<=nstep0(i,j)) nstep=nstep0(i,j)  
    enddo
   enddo
 enddo  
!---------------------
! if(myid>=0) then 
do ig=1,igas
     !!!!!!!!!!!!!! vertical mass conservarion by chenhs !!!!!!!!!!!!!!!
     if ( ig.eq.1 ) then
        if(imasskeep==1)then
           kpmass_m1= 1000. 
        endif
     endif
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if(imasskeep==1)then
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)
         ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         do k=1,nzz
            i03=ip3mem(k,ne)
            i04=ip4mem(k,ig,ne)
            kpmass_m2(i03+ixy)=gas(i04+ixy)
         enddo
      enddo
      enddo
     endif
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 do ii=1,nstep
     dtt0=dtstep(ne)/float(nstep)   !! by chenhs
    do j=sy(ne)-1,ey(ne)+1    
     do i=sx(ne)-1,ex(ne)+1
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz
       i03=ip3mem(k,ne)  
       i04=ip4mem(k,ig,ne)    
       uu(i,j,k) = u (i03+ixy)
       vv(i,j,k) = v (i03+ixy)
       ww(i,j,k) = w (i03+ixy)
       ddx(i,j,k)= dx(i03+ixy)             
       ddy(i,j,k)= dy(i03+ixy)             
       ddz(i,j,k)= dz(i03+ixy)
       conc(i,j,k)=gas(i04+ixy)
       concmark(i,j,k)=gas(i04+ixy)
       if ( ig.eq.1 ) then
         if(imasskeep==1)then
           kp_tmp(i,j,k)= kpmass_m1(i03+ixy)     
         endif
       endif
      enddo  !k
       i0=ip2mem(ne)
       kktop(i,j) = ktop(i0+ixy)
     enddo   !j
    enddo    !i
       
     CALL ADV_VERT(MYID,conc,ww,uu,vv,ddz,ddx,ddy,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dtt0,ig)
     !!!!!!!!!!!!!! vertical mass conservarion by chenhs !!!!!!!!!!!!!!!
     if ( ig.eq.1 ) then
        if(imasskeep==1)then           
     CALL ADV_VERT(MYID,kp_tmp,ww,uu,vv,ddz,ddx,ddy,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dtt0,ig)
        endif
     endif
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
       i03=ip3mem(k,ne)     
       i04=ip4mem(k,ig,ne)
       gas(i04+ixy) = conc(i,j,k)
       if ( ig.eq.1 ) then
        if(imasskeep==1)then
           kpmass_m1(i03+ixy)=kp_tmp(i,j,k)
        endif
       endif
       enddo !k
      enddo  !i
     enddo    !j      
               
  !!!!!!!!!!!!!!!!!!!!!!!!!
  ! for Source Mark
  ! to close vertical source mark
  ! 
  if(ifsm(ne)==1)then
    letdoit=0
    do idm=1,idmSet
       if(igMark(idm)==ig)letdoit=idm
    enddo
    if(letdoit>0)then
    iwb=sx(ne)-1;ieb=ex(ne)+1
    jsb=sy(ne)-1;jeb=ey(ne)+1
    allocate(smconv(ismMax,iwb:ieb,jsb:jeb,nzz))
           
    do j=sy(ne)-1,ey(ne)+1
      do i=sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
        do ism=1,ismMax
          i04sm=ipSMmem(k,ism,letdoit,ne)
          smconv(ism,i,j,k)=SourceMark(i04sm+ixy)
        enddo !is
       enddo !k
      enddo!i
    enddo !j
        
      CALL ADV_VERT_MARK(MYID,concmark,ww,uu,vv,ddz,ddx,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dtt0,ig,ismMax,smconv,kktop)
      
      do j=sy(ne),ey(ne)
        do i=sx(ne),ex(ne)
          ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         do k=1,nzz
          do ism=1,ismMax
           i04sm=ipSMmem(k,ism,letdoit,ne)
           SourceMark(i04sm+ixy)=smconv(ism,i,j,k)
          enddo !is
         enddo !k
        enddo!i
      enddo !j
    deallocate(smconv)
  endif !letdoit
  endif !ifsm(ne)
!!!!!!!!!!!!!!!!!!!!!!!!
  enddo !ii  
  !!!!!!!!!!!!!! vertical mass conservarion by chenhs !!!!!!!!!!!!!!!
    if ( ig.eq.1 ) then
        if(imasskeep==1)then
          do k=1,nzz   ! advection
             i02=ip2mem(ne)
             i03=ip3mem(k,ne)
             ! to get mass conservation ratio/error
             call  GetMassRatio(myid,kpmass_m1(i03),RatioMass(i03), &
                                sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
           enddo  
        endif
    endif
    if(imasskeep==1)then
    do k=1,nzz   !mass conservation 
       i03=ip3mem(k,ne)
       i04=ip4mem(k,ig,ne)
        
     call  balance(myid,gas(i04),kpmass_m2(i03),RatioMass(i03), &
                 sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
    enddo
    endif 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 enddo  !ig    

!--- FOR DUST AND SEA SALT

  DO IA=1,IAER
  DO IS=1,ISIZE

      if(imasskeep==1)then
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)
         ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         do k=1,nzz-1
            i03=ip3mem(k,ne)
!            i04=ip4mem(k,ig,ne) ! added by juanxiong he
            I05=IP5MEM(K,IS,IA,NE)
            kpmass_m2(i03+ixy)=AER(i05+ixy)            
         enddo
      enddo
      enddo
      endif

    do ii=1,nstep
     dtt0=dtstep(ne)/float(nstep)  !! by chenhs
    do j=sy(ne)-1,ey(ne)+1

     do i=sx(ne)-1,ex(ne)+1
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz
       i03=ip3mem(k,ne)
!       i04=ip4mem(k,ig,ne)  ! added by juanxiong he
       I05=IP5MEM(K,IS,IA,NE)
       uu(i,j,k) = u (i03+ixy)
       vv(i,j,k) = v (i03+ixy)
       ww(i,j,k) = w (i03+ixy)
       ddx(i,j,k)= dx(i03+ixy)
       ddy(i,j,k)= dy(i03+ixy)
       ddz(i,j,k)= dz(i03+ixy)
       conc(i,j,k)=AER(i05+ixy)
       concmark(i,j,k)=AER(i05+ixy)
      enddo  !k
       i0=ip2mem(ne)
       kktop(i,j) = ktop(i0+ixy)
     enddo   !j
    enddo    !i

       CALL ADV_VERT(MYID,conc, ww,uu,vv,ddz,ddx,ddy,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dtt0,ig)

     do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
!       i04=ip4mem(k,ig,ne)  ! added by juanxiong he
       I05=IP5MEM(K,IS,IA,NE)
       AER(i05+ixy) = conc(i,j,k)
       enddo !k
      enddo  !i
     enddo    !j      

    IF(ia==1) THEN ! for Sea salt
     if(ifseacom.eq.1) then
     do iduc = 1, nseacom
      do j=sy(ne)-1,ey(ne)+1
       do i=sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz
         i03=ip3mem(k,ne)
         i05c = ip5memcs (k,is,iduc,ne)
         uu(i,j,k) = u (i03+ixy)
         vv(i,j,k) = v (i03+ixy)
         ww(i,j,k) = w (i03+ixy)
         ddx(i,j,k)= dx(i03+ixy)
         ddy(i,j,k)= dy(i03+ixy)
         ddz(i,j,k)= dz(i03+ixy)
         conc(i,j,k)= SEACOMP(i05c+ixy)
        enddo 
       i0=ip2mem(ne)
       kktop(i,j) = ktop(i0+ixy)
     enddo   !j
    enddo    !i

       CALL ADV_VERT (MYID,conc,ww,uu,vv,ddz,ddx,ddy,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dtt0,ig)

     do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
       i05c = ip5memcs (k,is,iduc,ne)
       SEACOMP(i05c+ixy) = AMAX1(conc(i,j,k), 1.E-20)
       enddo !k
      enddo  !i
     enddo    !j   

     enddo ! iduc
     endif ! ifseacom
    ELSE IF (ia==2) then
     if(ifdustcom.eq.1) then
     do iduc = 1, ndustcom
      do j=sy(ne)-1,ey(ne)+1
       do i=sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz
         i03=ip3mem(k,ne)
         i05c = ip5memc (k,is,iduc,ne)
         uu(i,j,k) = u (i03+ixy)
         vv(i,j,k) = v (i03+ixy)
         ww(i,j,k) = w (i03+ixy)
         ddx(i,j,k)= dx(i03+ixy)
         ddy(i,j,k)= dy(i03+ixy)
         ddz(i,j,k)= dz(i03+ixy)
         conc(i,j,k)= DUSTCOMP(i05c+ixy)
        enddo
       i0=ip2mem(ne)
       kktop(i,j) = ktop(i0+ixy)
     enddo   !j
    enddo    !i

        CALL ADV_VERT (MYID,conc,ww,uu,vv,ddz,ddx,ddy,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dtt0,ig)

     do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
       i05c = ip5memc (k,is,iduc,ne)
       DUSTCOMP(i05c+ixy) = AMAX1(conc(i,j,k), 1.E-20)
       enddo !k
      enddo  !i
     enddo    !j   

     enddo ! iduc
     endif ! ifdustcom

    ENDIF    ! IA IF

    enddo !ii nstep
    
  ENDDO !IS
  ENDDO ! IA

 !--
 deallocate(uu,vv,ww,ddz,ddx,ddy,conc,concmark,kktop,kp_tmp)
! endif !myid

!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
if(ifprocess.eq.1) then
  ! iprocess
    do k=1,nzz
        do ig=1,iPrintTermGas    ! for gas phase
           igg=IGGPOS(ig)
           i04=ip4mem(k,igg,ne)
           i05=ipGasTermBal(k,3,ig,ne) ! 3--> vert adv
    call termbal(myid,gasOLD(i04),gas(i04),GasTermBal(i05), &
            k,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dtstep(ne))  !! change dt,chenhs
         enddo
     enddo
endif
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++

!-------------------------------------------------------------------
!----------------------cloud and convection-------------------------
!-------------------------------------------------------------------

    allocate(T1(nzz),Q1(nzz),QS1(nzz),U1(nzz),V1(nzz),TRA1(nzz,1),P1(nzz),&
             PH1(nzz+1),FTRA1(nzz,1),FTRA1D(nzz,1),FTRA1U(nzz,1),FTRA1O(nzz,1),FTRA1E(nzz,1))  
           !FTRA1 is net, FTRA1D is downdraft, FTRA1U is updraft, 
           !FTRA1O is the export form the cell, FTRA1E is enchange of environment 
    i0=ip2mem(ne)
   
   do j=sy(ne)-1,ey(ne)+1
   do i=sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
   do ig=1,igas
     dt0=dtstep(ne)   !600 seconds,cannot alter unless yo have enough reasons
     nstep=1
     ttime=0  
     CBMF1(i0+ixy)=1.0E-20
    do k=1,nzz
      i04=ip4mem(k,ig,ne)
      i03=ip3mem(k,ne)
      if(k/=1) i03_1=ip3mem(k-1,ne)
      U1(k)=u(i03+ixy)
      V1(k)=v(i03+ixy)
      T1(k)=t(i03+ixy)
      Q1(k)=QVAPOR(i03+ixy)/(QVAPOR(i03+ixy)+1.) !specific huminity
      QS1(k)=Q1(k)/rh1(i03+ixy)*100.  !sturation specific huminity
      P1(k)=Plev(i03+ixy)
      TRA1(k,1)= AMAX1(gas(i04+ixy) , 1.E-20)
      FTRA1(k,1)=1.0E-20
      if(k==1) PH1(k)=PSFC(i0+ixy)/100.
      if(k>1)  PH1(k)=Plev(i03_1+ixy)+Plev(i03_1+ixy)-PH1(k-1)
    enddo !k=nzz
      i03_1=ip3mem(nzz,ne)
      PH1(nzz+1)=Plev(i03_1+ixy)+Plev(i03_1+ixy)-PH1(nzz)
     
   DO iconv=1,1000 !1000 no meaning ,just for safty
990  CONTINUE
    dtt=dt0/nstep
    ttime=ttime+dtt
    FCBMF=CBMF1(i0+ixy)
    IFLAG=0    
    CALL  CONVECT(T1,Q1,QS1,U1,V1,TRA1,P1,PH1,nzz,nzz-1,1,dtt,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E,FCBMF,IFLAG) 
    IF (IFLAG .ne. 1 .and.IFLAG .ne. 4) goto 992
    IF(IFLAG==4.and.nstep<=16) THEN
      nstep=nstep*2       
      ttime=ttime-dtt
      GOTO 990
    ENDIF  
        CBMF1(i0+ixy)=FCBMF 
    
    IF(ttime<=7.*dt0.and.ttime>=dt0*6.) THEN !after initial (2hr)time

  !!!!!!!!!!!!!!!!!!!!!!!!!
  ! for Source Mark
  ! to close vertical source mark
  !
    if(ifsm(ne)==1)then
     allocate(smconv1(ismMax,nzz), c00(nzz), smthis(ismMax),smother(ismMax))

    letdoit=0
    do idm=1,idmSet
       if(igMark(idm)==ig)letdoit=idm
    enddo
    if(letdoit>0)then
      
      do ism=1,ismMax
       do k=1,nzz
        i04sm=ipSMmem(k,ism,letdoit,ne)
        i04=ip4mem(k,ig,ne)
        smconv1 (ism,k)= SourceMark(i04sm +ixy)
        c00(k)=gas(i04+ixy)
        c00(k)=amax1(c00(k), 1.E-20)
       enddo !k
      enddo !ism
 
    DO nn=1,40 !!no meaning, just reduce the errors due to the order of up and down
      dttnn=dtt/40. 
        i02=ip2mem(ne)     
        kkbb=ktop(i02+ixy) 
      do k=1,kkbb
          deltc1=dttnn*FTRA1D(k,1)
          deltc2=dttnn*FTRA1U(k,1)
          deltc3=dttnn*FTRA1O(k,1)
          deltc4=dttnn*FTRA1E(k,1)
    !-----------downdraft
     IF (deltc1>0.0) THEN
         do is=1,ismMax
          smthis(is)=smconv1(is,k)
          if(k==nzz)  smother(is)=smconv1(is,k)
          if(k<nzz)   smother(is)=smconv1(is,k+1)
         enddo
         if(k==ktop(i02+ixy)) then
          call GetBounChange(deltc1,c00(k),smthis,2,ismMax)
         else if(k>ktop(i02+ixy)) then
          smthis(2)=1.
         else if(k<ktop(io2+ixy))then
          call SMmixing(c00(k),smthis,deltc1,smother,ismMax)
         endif !k
         
         do is=1,ismMax
           smconv1(is,k)=smthis(is)
         enddo
      ENDIF
         c00(k)=c00(k)+deltc1
         c00(k)=amax1(c00(k),1.E-20)
    !---------updraft
       IF(deltc2>0.0)THEN
         do is=1,ismMax
           smthis(is)=smconv1(is,k)
           if(k==1) smother(is)=smconv1(is,1)
           if(k>1) smother(is)=smconv1(is,k-1)
         enddo
       
         if(k>ktop(i02+ixy)) then
          smthis(2)=1.
         else if(k<=ktop(io2+ixy)) then
          call SMmixing(c00(k),smthis,deltc2,smother,ismMax)
         endif!k
  
         do is=1,ismMax
           smconv1(is,k)=smthis(is)
         enddo
        ENDIF
          c00(k)=c00(k)+deltc2
          c00(k)=c00(k)+deltc3
          c00(k)=c00(k)+deltc4
         c00(k)=amax1(c00(k),1.E-20)
      enddo !k=ktop
    ENDDO !nn

      do ism=1,ismMax
       do k=1,nzz
        i04sm=ipSMmem(k,ism,letdoit,ne)
        SourceMark(i04sm +ixy)=smconv1(ism,k)
       enddo !k
      enddo !ism
    
   endif !letdoit
   deallocate(smconv1,c00,smthis,smother)
 endif !ifsm
             
     do k=1,nzz
      i04=ip4mem(k,ig,ne)
      gas(i04+ixy)=TRA1(k,1)+FTRA1(k,1)*dtt  
     enddo!k
    ENDIF !ttime
    
     IF(nint(ttime)>dt0*7.) GOTO 991
     
   ENDDO !iconv
   
 991  CONTINUE
    enddo !igas   
 992 CONTINUE

  enddo ! j
  enddo ! i

  deallocate(T1,Q1,QS1,U1,V1,TRA1,P1,PH1,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E)

!------------------------------------------------
! ---- gas and sec aerosols
! ---- FOR SEA SALT AND DUST EMISSIONS ----
!------------------------------------------------
 allocate(T1(nzz),Q1(nzz),QS1(nzz),U1(nzz),V1(nzz),TRA1(nzz,1),P1(nzz),&
           PH1(nzz+1),FTRA1(nzz,1),FTRA1D(nzz,1),FTRA1U(nzz,1),FTRA1O(nzz,1),FTRA1E(nzz,1))
           !FTRA1 is net, FTRA1D is downdraft, FTRA1U is updraft, 
           !FTRA1O is the export form the cell, FTRA1E is enchange of
           !environment 

    i0=ip2mem(ne)

   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

  DO IA=1,IAER !! why chenhs
  DO IS=1,ISIZE
 
     dt0=dtstep(ne)   !600 seconds,cannot alter unless yo have enough reasons
     nstep=1
     ttime=0
     CBMF1(i0+ixy)=0.0
   do k=1,nzz
      i04=ip4mem(k,ig,ne)
      i03=ip3mem(k,ne)
      I05=IP5MEM(K,IS,IA,NE)
      if(k/=1) i03_1=ip3mem(k-1,ne)
      U1(k)=u(i03+ixy)
      V1(k)=v(i03+ixy)
      T1(k)=t(i03+ixy)
      Q1(k)=QVAPOR(i03+ixy)/(QVAPOR(i03+ixy)+1.) !specific huminity
      QS1(k)=Q1(k)/rh1(i03+ixy)*100.  !sturation specific huminity
      P1(k)=Plev(i03+ixy)
      TRA1(k,1)=AMAX1(AER(i05+ixy), 1.E-20)
      FTRA1(k,1)=0.0
      if(k==1) PH1(k)=PSFC(i0+ixy)/100.
      if(k>1)  PH1(k)=Plev(i03_1+ixy)+Plev(i03_1+ixy)-PH1(k-1)
    enddo !k=nzz
      i03_1=ip3mem(nzz,ne)
      PH1(nzz+1)=Plev(i03_1+ixy)+Plev(i03_1+ixy)-PH1(nzz)

   DO iconv=1,1000 !1000 no meaning ,just for safty
993  CONTINUE
    dtt=dt0/nstep
    ttime=ttime+dtt
    FCBMF=CBMF1(i0+ixy)
    IFLAG=0
    CALL CONVECT(T1,Q1,QS1,U1,V1,TRA1,P1,PH1,nzz,nzz-1,1,dtt,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E,FCBMF,IFLAG)
    IF (IFLAG .ne. 1 .and.IFLAG .ne. 4) goto 995
    IF(IFLAG==4.and.nstep<=16) THEN
      nstep=nstep*2
      ttime=ttime-dtt
      GOTO 993
    ENDIF
        CBMF1(i0+ixy)=FCBMF

    IF(ttime<=7.*dt0.and.ttime>=dt0*6.) THEN !after initial (2hr)time

     do k=1,nzz
      i04=ip4mem(k,ig,ne)
      I05=IP5MEM(K,IS,IA,NE)
      AER(i05+ixy)=TRA1(k,1)+FTRA1(k,1)*dtt
     enddo!k
    ENDIF !ttime

     IF(nint(ttime)>dt0*7.) GOTO 994

   ENDDO !iconv
994  CONTINUE

    IF(IA==1) THEN ! FOR SEA SALT
     if(ifseacom.eq.1) then
     DO iduc = 1, nseacom
     dt0=dtstep(ne)   !600 seconds,cannot alter unless yo have enough reasons
     nstep=1
     ttime=0
     CBMF1(i0+ixy)=0.0
   do k=1,nzz
      i04=ip4mem(k,ig,ne)
      i03=ip3mem(k,ne)
      i05c = ip5memcs(k,is,iduc,ne)
      if(k/=1) i03_1=ip3mem(k-1,ne)
      U1(k)=u(i03+ixy)
      V1(k)=v(i03+ixy)
      T1(k)=t(i03+ixy)
      Q1(k)=QVAPOR(i03+ixy)/(QVAPOR(i03+ixy)+1.) !specific huminity
      QS1(k)=Q1(k)/rh1(i03+ixy)*100.  !sturation specific huminity
      P1(k)=Plev(i03+ixy)
      TRA1(k,1)=AMAX1(SEACOMP(i05c+ixy), 1.E-20)
      FTRA1(k,1)=0.0
      if(k==1) PH1(k)=PSFC(i0+ixy)/100.
      if(k>1)  PH1(k)=Plev(i03_1+ixy)+Plev(i03_1+ixy)-PH1(k-1)
    enddo !k=nzz
      i03_1=ip3mem(nzz,ne)
      PH1(nzz+1)=Plev(i03_1+ixy)+Plev(i03_1+ixy)-PH1(nzz)
   DO iconv=1,1000 !1000 no meaning ,just for safty
996  CONTINUE
    dtt=dt0/nstep
    ttime=ttime+dtt
    FCBMF=CBMF1(i0+ixy)
    IFLAG=0
    CALL CONVECT(T1,Q1,QS1,U1,V1,TRA1,P1,PH1,nzz,nzz-1,1,dtt,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E,FCBMF,IFLAG)
    IF (IFLAG .ne. 1 .and.IFLAG .ne. 4) goto 995
    IF(IFLAG==4.and.nstep<=16) THEN
      nstep=nstep*2
      ttime=ttime-dtt
      GOTO 996
    ENDIF
        CBMF1(i0+ixy)=FCBMF

    IF(ttime<=7.*dt0.and.ttime>=dt0*6.) THEN !after initial (2hr)time

     do k=1,nzz
      i04=ip4mem(k,ig,ne)
      i05c = ip5memcs(k,is,iduc,ne)
      SEACOMP(i05c+ixy)=TRA1(k,1)+FTRA1(k,1)*dtt
     enddo!k
    ENDIF !ttime

     IF(nint(ttime)>dt0*7.) GOTO 997

   ENDDO !iconv
 997  CONTINUE

     ENDDO ! iduc
     endif ! ifseacom

    ELSE IF (IA==2) THEN ! FOR DUST 
     if(ifdustcom.eq.1) then
     DO iduc = 1, ndustcom
     dt0=dtstep(ne)   !600 seconds,cannot alter unless yo have enough reasons
     nstep=1
     ttime=0
     CBMF1(i0+ixy)=0.0
   do k=1,nzz
      i04=ip4mem(k,ig,ne)
      i03=ip3mem(k,ne)
      i05c = ip5memc (k,is,iduc,ne)
      if(k/=1) i03_1=ip3mem(k-1,ne)
      U1(k)=u(i03+ixy)
      V1(k)=v(i03+ixy)
      T1(k)=t(i03+ixy)
      Q1(k)=QVAPOR(i03+ixy)/(QVAPOR(i03+ixy)+1.) !specific huminity
      QS1(k)=Q1(k)/rh1(i03+ixy)*100.  !sturation specific huminity
      P1(k)=Plev(i03+ixy)
      TRA1(k,1)=AMAX1(DUSTCOMP(i05c+ixy), 1.E-20)
      FTRA1(k,1)=0.0
      if(k==1) PH1(k)=PSFC(i0+ixy)/100.
      if(k>1)  PH1(k)=Plev(i03_1+ixy)+Plev(i03_1+ixy)-PH1(k-1)
    enddo !k=nzz
      i03_1=ip3mem(nzz,ne)
      PH1(nzz+1)=Plev(i03_1+ixy)+Plev(i03_1+ixy)-PH1(nzz)
   DO iconv=1,1000 !1000 no meaning ,just for safty
892  CONTINUE
    dtt=dt0/nstep
    ttime=ttime+dtt
    FCBMF=CBMF1(i0+ixy)
    IFLAG=0
    CALL CONVECT(T1,Q1,QS1,U1,V1,TRA1,P1,PH1,nzz,nzz-1,1,dtt,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E,FCBMF,IFLAG)
    IF (IFLAG .ne. 1 .and.IFLAG .ne. 4) goto 995
    IF(IFLAG==4.and.nstep<=16) THEN
      nstep=nstep*2
      ttime=ttime-dtt
      GOTO 892
    ENDIF
        CBMF1(i0+ixy)=FCBMF

    IF(ttime<=7.*dt0.and.ttime>=dt0*6.) THEN !after initial (2hr)time

     do k=1,nzz
      i04=ip4mem(k,ig,ne)
      i05c = ip5memc (k,is,iduc,ne)
      DUSTCOMP(i05c+ixy)=TRA1(k,1)+FTRA1(k,1)*dtt
     enddo!k
    ENDIF !ttime

    IF(nint(ttime)>dt0*7.) GOTO 998

   ENDDO !iconv
 998  CONTINUE
     ENDDO ! iduc
     endif ! ifdustcom
    ENDIF ! IAIF

    ENDDO ! IS
    ENDDO ! IA  

 995 CONTINUE

   enddo !i
  enddo !j
    deallocate(T1,Q1,QS1,U1,V1,TRA1,P1,PH1,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E)

!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
if(ifprocess.eq.1) then
! iprocess
 do k=1,nzz
   do ig=1,iPrintTermGas    ! for gas phase
    i04=ip4mem(k,igg,ne)
    i05=ipGasTermBal(k,21,ig,ne) ! 21--> moist convect
    call termbal(myid,gasOLD(i04),gas(i04),GasTermBal(i05), &
      k,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dtstep(ne))  !! change dt,chenhs
   enddo
  enddo
endif
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++

!------------------------------------------------
!-----------------------diffusion
!------------------------------------------------
 do ig=1,igas    ! for gas phase
  do k=1,nzz-1   
    i03=ip3mem(k,ne)
    i04=ip4mem(k,ig,ne)
    i02=ip2mem(ne)
    call  dif_hori( myid, gas(i04), kh(i03), dx(i03), dy(i03), &
            sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dtstep(ne))  !! change dt,chenhs
  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!
  ! for Source Mark
  if(ifsm(ne)==1)then
    letdoit=0
    do idm=1,idmSet
       if(igMark(idm)==ig)letdoit=idm
    enddo
    if(letdoit>0)then
    allocate(sm(ismMax,sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1))
    do k=1,nzz-1
       i03=ip3mem(k,ne)
       i04=ip4mem(k,ig,ne)
       do ism=1,ismMax
          i04sm=ipSMmem(k,ism,letdoit,ne)
              do i=sx(ne)-1,ex(ne)+1
              do j=sy(ne)-1,ey(ne)+1
               ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
               sm(ism,i,j)= SourceMark(i04sm+ixy)
              enddo
              enddo
        enddo
     i02=ip2mem(ne)   
     call dif_hori_mark( myid, gas(i04), kh(i03), dx(i03),dy(i03), &
                  sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dtstep(ne),&  !! change dt,chenhs
                  ismMax,sm, ne,k,ktop(i02))
         do ism=1,ismMax
          i04sm=ipSMmem(k,ism,letdoit,ne)
              do i=sx(ne)-1,ex(ne)+1
              do j=sy(ne)-1,ey(ne)+1
               ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
               SourceMark(i04sm+ixy)=sm(ism,i,j)
              enddo
              enddo
        enddo
    enddo

    deallocate(sm)
    endif
  endif

 enddo    !ig

!-----------------------------------
!--- FOR DUST AND SEA SALT  -------- 
!-----------------------------------
  do ia=1,iaer
  do is=1,isize
      do k=1,nzz-1   
     i03=ip3mem(k,ne)
     i05=ip5mem(k,is,ia,ne)
     call  dif_hori( myid, aer(i05), kh(i03), dx(i03), dy(i03), &
                    sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dtstep(ne))  !! change dt,chenhs
       enddo   

   IF(IA==1) THEN ! For sea salt
    if(ifseacom.eq.1) then
    do iduc = 1, nseacom
     do k = 1, nzz-1
     i03=ip3mem(k,ne)
     i05c = ip5memcs (k,is,iduc,ne)
     call  dif_hori( myid, SEACOMP(i05c), kh(i03), dx(i03), dy(i03), &
            sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dtstep(ne))   !! change dt,chenhs
     enddo
    enddo ! iduc
    endif ! ifseacom
   ELSE IF(IA==2) THEN ! for dust
    if(ifdustcom.eq.1) then
    do iduc = 1, ndustcom
     do k = 1, nzz-1
     i03=ip3mem(k,ne)
     i05c = ip5memc (k,is,iduc,ne)
     call  dif_hori( myid, DUSTCOMP(i05c), kh(i03), dx(i03), dy(i03), &
            sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dtstep(ne))     !! change dt,chenhs
     enddo
    enddo ! iduc
    endif ! ifdustcom
   ENDIF ! IA IF
 
  enddo        !isize
  enddo        !iaer

!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
if(ifprocess.eq.1) then
! iprocess
  do k=1,nzz
      do ig=1,iPrintTermGas    ! for gas phase
          igg=IGGPOS(ig)
          i04=ip4mem(k,igg,ne)
          i05=ipGasTermBal(k,4,ig,ne) ! 4--> hor diff
          call termbal(myid,gasOLD(i04),gas(i04),GasTermBal(i05), &
                k,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dtstep(ne))  !! change dt,chenhs
      enddo
 enddo
endif
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++

!------------------------------------------------
!---------------------vertical diffusion---------
!------------------------------------------------
!======================to get the kv profile======================
     do k=1,nzz
      i02=ip3mem(k,ne)
      i0=ip2mem(ne)
      CALL eddyz(myid,UST(i0),PBL_HGT(i0),RMOL(i0),&
                 heiz(i02),HGT1(i0),kv(i02),sx(ne),ex(ne),sy(ne),ey(ne),k) !! h change to heiz, by chenhs
     enddo

!==================================kv======================== 
IF(idifvert==1) THEN
  iitime = time(ne)
  iwb=sx(ne)-1;ieb=ex(ne)+1
  jsb=sy(ne)-1;jeb=ey(ne)+1
 allocate(ppp(iwb:ieb,jsb:jeb,nzz),ttn(iwb:ieb,jsb:jeb,nzz),ffn(iwb:ieb,jsb:jeb,nzz),&
         conc(iwb:ieb,jsb:jeb,nzz),concmark(iwb:ieb,jsb:jeb,nzz),rkv(iwb:ieb,jsb:jeb,nzz),&
         dzz(iwb:ieb,jsb:jeb,nzz),atm(iwb:ieb,jsb:jeb,nzz),kktop(iwb:ieb,jsb:jeb))
 if(ifseacom.eq.1)  allocate(SEA(iwb:ieb,jsb:jeb,nzz,nseacom))
 if(ifdustcom.eq.1) allocate(DUST(iwb:ieb,jsb:jeb,nzz,ndustcom))
  do j=sy(ne)-1,ey(ne)+1
     do i=sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz
        i03=ip3mem(k,ne)
        i0_1=ip2mem2dconv(ne)
        i0=ip2mem(ne)
        IF(i>(sx(ne)-1).and.i<(ex(ne)+1).and.j>(sy(ne)-1).and.j<(ey(ne)+1))   THEN
         ppp(i,j,k)=Plev(i03+ixy)
         ttn(i,j,k)=t(i03+ixy)
         dzz(i,j,k)=dz(i03+ixy)
         rkv(i,j,k)=kv(i03+ixy)
         ffn(i,j,k)=rh1(i03+ixy)
         atm(i,j,k) =  PA2ATM*Plev(i03+ixy)*100.
        ENDIF
        enddo !k
         kktop(i,j) = ktop(i0+ixy)
      enddo !i
  enddo !j        
  
  DO ig=1,igas
  do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
    do k=1,nzz
     i04=ip4mem(k,ig,ne)  
     conc(i,j,k)=gas(i04+ixy)
     concmark(i,j,k)=gas(i04+ixy)
    enddo !k
   enddo!i
  enddo !j
   jsb1=sy(ne);jeb1=ey(ne)
   iwb1=sx(ne);ieb1=ex(ne)
   call diffus(myid,ig,iwb1,ieb1,jsb1,jeb1,nzz,dtstep(ne),rkv,dzz,ttn,ppp,conc,atm,int(kktop),GC_MOLWT(ig)) !! add GC_MOLWT,chenhs
   do j=sy(ne),ey(ne)
     do i=sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
       i04=ip4mem(k,ig,ne)
       i03=ip3mem(k,ne)
       gas(i04+ixy)=conc(i,j,k)
       enddo !k
     enddo!i
   enddo !j
  !!!!!!!!!!!!!!!!!!!!!!!!!
  ! for Source Mark
   if(ifsm(ne)==1)then
     letdoit=0
     do idm=1,idmSet
      if(igMark(idm)==ig) letdoit=idm
     enddo
        
     if(letdoit>0)then
        iwb=sx(ne)-1;ieb=ex(ne)+1
        jsb=sy(ne)-1;jeb=ey(ne)+1
       allocate(smconv(ismMax,iwb:ieb,jsb:jeb,nzz))
       do j=sy(ne),ey(ne)
        do i=sx(ne),ex(ne)
           ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz
         do ism=1,ismMax
           i04sm=ipSMmem(k,ism,letdoit,ne)
           smconv(ism,i,j,k)=SourceMark(i04sm+ixy) 
          enddo !is
        enddo !k
        enddo!i
        enddo !j
    jsb1=sy(ne);jeb1=ey(ne)
    iwb1=sx(ne);ieb1=ex(ne)
    call diffus_mark(myid,ig,iwb1,ieb1,jsb1,jeb1,nzz,dtstep(ne),rkv,dzz,ttn,ppp,concmark,atm,ismMax,smconv,int(kktop)) !! change dt,chenhs
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)    
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz
       do ism=1,ismMax
         i04sm=ipSMmem(k,ism,letdoit,ne)
         SourceMark(i04sm+ixy)=smconv(ism,i,j,k)
       enddo !is
      enddo !k
    enddo!i
    enddo !j
    
  deallocate(smconv)
 endif !letdoit
  endif !ifsm(ne)
  !!!!!!!!!!!!!!  
  ENDDO !ig

!----- FOR DUST AND SEA SALT ----
  DO IA=1,IAER
  DO IS=1,ISIZE
    do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1    
       do k=1,nzz       
        i05=ip5mem(k,is,ia,ne)
        conc(i,j,k)=aer(i05+ixy)    
        if(ifseacom.eq.1) then     
         do iduc = 1,nseacom
         i05c = ip5memcs (k,is,iduc,ne)
         sea(i,j,k,iduc) = SEACOMP(i05c+ixy)
         enddo
        endif
        if(ifdustcom.eq.1) then
         do iduc = 1, ndustcom
         i05c = ip5memc (k,is,iduc,ne)
         dust(i,j,k,iduc) = DUSTCOMP(i05c+ixy)        
         enddo ! iduc 
        endif
       enddo !k
      enddo!i
    enddo !j
    jsb1=sy(ne);jeb1=ey(ne)
    iwb1=sx(ne);ieb1=ex(ne)

    if(ia==1) then
     if(ifseacom.eq.1) then
     call diffus_ds_com(myid,ig,iwb1,ieb1,jsb1,jeb1,nzz,dtstep(ne),rkv,dzz,ttn,ppp,conc,sea,atm,int(kktop),nseacom,is,ia)  !! change dt,chenhs
     else
     call diffus_ds(myid,ig,iwb1,ieb1,jsb1,jeb1,nzz,dtstep(ne),rkv,dzz,ttn,ppp,conc,atm,int(kktop))  !! change dt,chenhs
     endif
    else if(ia==2)then 
     if(ifdustcom.eq.1) then
     call diffus_ds_com(myid,ig,iwb1,ieb1,jsb1,jeb1,nzz,dtstep(ne),rkv,dzz,ttn,ppp,conc,dust,atm,int(kktop),ndustcom,is,ia) !! change dt,chenhs
     else
     call diffus_ds(myid,ig,iwb1,ieb1,jsb1,jeb1,nzz,dtstep(ne),rkv,dzz,ttn,ppp,conc,atm,int(kktop)) !! change dt,chenhs
     endif
    endif

   do j=sy(ne),ey(ne)
     do i=sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
       i05=ip5mem(k,is,ia,ne)       
       aer(i05+ixy)=conc(i,j,k)
        
        if(ifseacom.eq.1) then
        do iduc = 1,nseacom
        i05c = ip5memcs (k,is,iduc,ne)
        SEACOMP(i05c+ixy) = sea(i,j,k,iduc)
        enddo
        endif
        if(ifdustcom.eq.1) then
        do iduc = 1, ndustcom
        i05c = ip5memc (k,is,iduc,ne)
        DUSTCOMP(i05c+ixy) = dust(i,j,k,iduc)
        enddo ! iduc
        endif

       enddo !k
     enddo!i
   enddo !j

  ENDDO ! IS
  ENDDO ! IA 

  deallocate(ppp,ttn,ffn,rkv,dzz,atm,conc,concmark,kktop)
  if(ifseacom.eq.1)  deallocate(sea)
  if(ifdustcom.eq.1) deallocate(dust)
ENDIF !idifvert  

!--------------------------------------------------------------
!---------------------- make balance -------------------------!
!--------------------------------------------------------------
if(ifbalance.eq.1.and.ne.eq.1) then
  do ig=103,105
     amount_bf(myid+1,ig)=0.0
     do k=1,nzz
        I03=IP3MEM(K,NE)
        I04=IP4MEM(K,IG,NE)
        call N_Balance(MYID,gas(I04),amount_bf(myid+1,ig),DX(I03),DY(I03),DZ(I03),SX(NE),EX(NE),SY(NE),EY(NE))
     enddo        
  enddo

  CDNUM='c000'
  WRITE(CDNUM(2:4),'(I3.3)')MYID
  OPEN(10,FILE='out/Amount.bf.dat'//CDNUM,FORM='FORMATTED')
  WRITE(10,*) MYID,amount_bf(myid+1,103),amount_bf(myid+1,104),amount_bf(myid+1,105)
  close(10)

  CALL mpi_barrier( local_com,IERR )
   total_bf=total_af
   total_af=0.
   DO I=0,NUMPROCS-1
      CDNUM='c000'
      total_tmp=0.
      WRITE(CDNUM(2:4),'(I3.3)')I
      OPEN(10,FILE='out/Amount.bf.dat'//CDNUM,FORM='FORMATTED')
      read(10,*) iiii,total_tmp(1),total_tmp(2),total_tmp(3)
      CLOSE(10)
      total_af=total_af+sum(total_tmp(1:3))
   ENDDO

   do ig=103,105
     amount_bf(myid+1,ig)=0.0
     do k=1,nzz
        I03=IP3MEM(K,NE)
        I04=IP4MEM(K,IG,NE)
        call make_balance(MYID,gas(I04),total_bf,total_af,DX(I03),DY(I03),DZ(I03),SX(NE),EX(NE),SY(NE),EY(NE))
        call N_Balance(MYID,gas(I04),amount_bf(myid+1,ig),DX(I03),DY(I03),DZ(I03),SX(NE),EX(NE),SY(NE),EY(NE))
     enddo
   enddo
   CDNUM='c000'
   WRITE(CDNUM(2:4),'(I3.3)')MYID
   OPEN(10,FILE='out/Amount.bf.dat'//CDNUM,FORM='FORMATTED')
   WRITE(10,*) MYID,amount_bf(myid+1,103),amount_bf(myid+1,104),amount_bf(myid+1,105)
   close(10)
    
   IF(NUMPROCS>1) CALL mpi_barrier( local_com,IERR )
   total_af=0.
   DO I=0,NUMPROCS-1
      CDNUM='c000'
      total_tmp=0.
      WRITE(CDNUM(2:4),'(I3.3)')I
      OPEN(10,FILE='out/Amount.bf.dat'//CDNUM,FORM='FORMATTED')
      read(10,*) iiii,total_tmp(1),total_tmp(2),total_tmp(3)
      CLOSE(10)
      total_af=total_af+sum(total_tmp(1:3))
   ENDDO

  if(myid.eq.0) then
   open(11,file='out/real_balance.txt',position='append')
   write(11,'(a60,f16.2,2x,f16.6)') "After DIF_VERT,Global=",total_af,(total_af-total_bf)/total_bf*100
   close(11)
  endif
endif ! ifbalance

!-------------------------------------------------------------------
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
!-------------------------------------------------------------------
if(ifprocess.eq.1) then
   ! iprocess
 do k=1,nzz
   do ig=1,iPrintTermGas    ! for gas phase
     igg=IGGPOS(ig)
     i04=ip4mem(k,igg,ne)
     i05=ipGasTermBal(k,5,ig,ne) ! 5--> ver diff
     call termbal(myid,gasOLD(i04),gas(i04),GasTermBal(i05), &
                   k,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dtstep(ne))  !! change dt,chenhs
   enddo
 enddo
endif ! ifprocess

!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
!       for dust and seasalt gravity deposition 
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
  dttmp=dtstep(ne)/5
  DO ITGRV=1,5
   DO K = 1 , NZZ 
    IF(K.GT.1)THEN
     I03_1=IP3MEM(K-1,NE)
    ELSE
     I03_1=IP3MEM(1,NE)
    ENDIF

     I03=IP3MEM(K,NE)

    IF(K.LT.NZZ)THEN
     I03P1=IP3MEM(K+1,NE)
    ELSE
     I03P1=IP3MEM(K,NE)
    ENDIF

  ! FOR DUST AND SEA SALT AEROSOLS
  DO IA=1,IAER
  DO IS=1,ISIZE
    I02 = IP2MEM(NE)

     IF(K.GT.1)THEN
       I05_1=IP5MEM(K-1,IS,IA,NE)
     ELSE
       I05_1=IP5MEM(1,IS,IA,NE)
     ENDIF

     I05=IP5MEM(K,IS,IA,NE)

     IF(K.LT.NZZ)THEN
       I05P1 = IP5MEM(K+1,IS,IA,NE)
     ELSE
       I05P1 = IP5MEM(K,IS,IA,NE)
     ENDIF

  ! TO ADD GRAVITY SETTLING ON JANUARY 09,2005
     GVEL=GRAVEL(IS,IA)  
     if(IA==2) then  !! added by chenhs
     CALL DUSTGRADEP(MYID,AER(I05),AER(I05P1),DUSTGRAV(I02),GVEL,K,NZZ,&
                      DZ(I03),DZ(I03P1),SX(NE),EX(NE),SY(NE),EY(NE),dttmp)  !! change dt,chenhs
     endif   
     CALL GRA_DEP_AER(MYID,AER(I05),AER(I05P1),GVEL,K,NZZ,&
                      DZ(I03),DZ(I03P1),SX(NE),EX(NE),SY(NE),EY(NE),dttmp) !! change dt,chenhs

    IF(IA==1)  THEN ! FOR SEA SALT 
      if(ifseacom.eq.1) then
      do iduc = 1, nseacom
       IF(K.GT.1)THEN
        I05_1=IP5MEMcs(K-1,IS,iduc,NE)
       ELSE
        I05_1=IP5MEMcs(1,IS,iduc,NE)
       ENDIF

        I05=IP5MEMcs(K,IS,iduc,NE)

       IF(K.LT.NZZ)THEN
        I05P1 = IP5MEMcs(K+1,IS,iduc,NE)
       ELSE
        I05P1 = IP5MEMcs(K,IS,iduc,NE)
       ENDIF

        GVEL=GRAVEL(IS,IA)
        CALL GRA_DEP_AER(MYID,SEACOMP(I05),SEACOMP(I05P1),GVEL,K,NZZ,&
                      DZ(I03),DZ(I03P1),SX(NE),EX(NE),SY(NE),EY(NE),dttmp)   !! change dt,chenhs 
       
       enddo ! iduc
       endif ! ifseacom
    ELSE IF (IA==2) THEN ! DUST
      if(ifdustcom.eq.1) then
      do iduc = 1, ndustcom
       IF(K.GT.1)THEN
        I05_1=IP5MEMc(K-1,IS,iduc,NE)
       ELSE
        I05_1=IP5MEMc(1,IS,iduc,NE)
       ENDIF

        I05=IP5MEMc(K,IS,iduc,NE)

       IF(K.LT.NZZ)THEN
        I05P1 = IP5MEMc(K+1,IS,iduc,NE)
       ELSE
        I05P1 = IP5MEMc(K,IS,iduc,NE)
       ENDIF

        GVEL=GRAVEL(IS,IA)
       IF(iduc==7) THEN
        CALL DUSTGRADEP(MYID, DUSTCOMP(I05), DUSTCOMP(I05P1),DUSTGRAVSO4(I02),GVEL,K,NZZ,&
                      DZ(I03),DZ(I03P1),SX(NE),EX(NE),SY(NE),EY(NE),dttmp)  !! change dt,chenhs
       ELSE IF (iduc==8) THEN
        CALL DUSTGRADEP(MYID, DUSTCOMP(I05), DUSTCOMP(I05P1),DUSTGRAVNO3(I02),GVEL,K,NZZ,&
                      DZ(I03),DZ(I03P1),SX(NE),EX(NE),SY(NE),EY(NE),dttmp)  !! change dt,chenhs
       ELSE IF(iduc==9) THEN
        CALL DUSTGRADEP(MYID, DUSTCOMP(I05),DUSTCOMP(I05P1),DUSTGRAVFEII(I02),GVEL,K,NZZ,&
                      DZ(I03),DZ(I03P1),SX(NE),EX(NE),SY(NE),EY(NE),dttmp)  !! change dt,chenhs
       ELSE IF(iduc==10) THEN
        CALL DUSTGRADEP(MYID, DUSTCOMP(I05), DUSTCOMP(I05P1),DUSTGRAVFEIII(I02),GVEL,K,NZZ,&
                      DZ(I03),DZ(I03P1),SX(NE),EX(NE),SY(NE),EY(NE),dttmp)  !! change dt,chenhs
       ENDIF

        CALL GRA_DEP_AER(MYID,DUSTCOMP(I05),DUSTCOMP(I05P1),GVEL,K,NZZ,&
                      DZ(I03),DZ(I03P1),SX(NE),EX(NE),SY(NE),EY(NE),dttmp)  !! change dt,chenhs
       
       enddo ! iduc
       endif ! ifdustcom

    ENDIF

  ENDDO ! IS
  ENDDO ! IA

  ENDDO !K
   
  ENDDO ! ITGRV

!ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
!ddddddddd             Dry    Deposition                   ddddddddddd
!ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
 if(idry==2) then
! This is to read drydeposition velocity from calculated value by lijie
!!!! for gas (CBM-Z)-----------------------------------------+
!                 1        2      3      4      5      6       !
!                 H2SO4   HNO3   HCL    NH3    NO     NO2
!                 7        8      9      10     11     12      +
!                 NO3     N2O5   HONO   HNO4   O3     O1D
!                 13       14     15     16     17     18      +
!                 O3P     OH     HO2    H2O2   CO     SO2
!                 19       20     21     22     23     24      +
!                 CH4     C2H6    CH3O2  ETHP   HCHO  CH3OH
!                 25       26     27     28     29     30
!                 ANOL    CH3OOH  ETHOOH ALD2   HCOOH RCOOH
!                 31       32   , 33   , 34     35     36
!                 C2O3    PAN     PAR   AONE   MGLY   ETH
!                 37       38     39     40     41     42
!                 OLET  , OLEI ,  TOL   XYL    CRES   TO2    
!                 43       44     45     46     47     48
!                 CRO     OPEN    ONIT  ROOH   RO2    ANO2
!                 49       50     51    52      53     54
!                 NAP     XO2     XPAR  ISOP   ISOPRD ISOPP
!                 55       56     57    58      59     60 
!                 ISOPN   ISOPO2  DMS  MSA     DMSO   DMSO2
!                 61       62         63     64     65     66
!                 CH3SO2H CH3SCH2OO CH3SO2  CH3SO3 CH3SO2OO CH3SO2CH2OO
!                  67     68      69     70    71    72
!                 SULFH   TERP    SV1    SV2   SV3   SV4
!                  73     74
!                  SV5    SV6
!                  75     76     77     78   79      80
!                   PM25  PM10   BC     OC   H+(AQ) NA+(AQ)
!                   81         82       83   84          85
!                 NH4+(AQ) CL-(AQ) SO4--(AQ) HSO4-(AQ)  NO3-(AQ)
!                  86      87        88         89        90
!                 NACL(S)  NA2SO4(S) NANO3(S)  NH42SO4(S) NH4NO3(S)
!                  91       92        93          94         95
!                  NH4CL(S) H2SO4(AQ) NH4HSO4(S)  NAHSO4(S) (NH4)4H(SO4)2(S)
!                  96      97        98      99    100   101, 102 
!                  SOA1    SOA2     SOA3    SOA4   SOA5  SOA6,AH2O

!!!-----------------------------------------------------------+    
  iwb=sx(ne)-1;ieb=ex(ne)+1
  jsb=sy(ne)-1;jeb=ey(ne)+1
  allocate(ppp(iwb:ieb,jsb:jeb,nzz),ttn(iwb:ieb,jsb:jeb,nzz),land(iwb:ieb,jsb:jeb),&
           tsurf(iwb:ieb,jsb:jeb),xlat(iwb:ieb,jsb:jeb),xlon(iwb:ieb,jsb:jeb),&
           QQ(iwb:ieb,jsb:jeb,nzz),height1(iwb:ieb,jsb:jeb,nzz),uu(iwb:ieb,jsb:jeb,nzz),&
           vv(iwb:ieb,jsb:jeb,nzz),water(iwb:ieb,jsb:jeb,nzz),SWDOWN1(iwb:ieb,jsb:jeb),&
           vdep(iwb:ieb,jsb:jeb),pblhgt(iwb:ieb,jsb:jeb),icddrt(iwb:ieb,jsb:jeb),icdsno(iwb:ieb,jsb:jeb))
  do j=sy(ne)-1,ey(ne)+1
    do i=sx(ne)-1,ex(ne)+1
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1      
     do k=1,nzz
      i03=ip3mem(k,ne)
      i0=ip2mem(ne)
!      i02Gas = ip2memGas(ig,ne)  ! added by juanxiong he
       ppp(i,j,k)=Plev(i03+ixy)
       ttn(i,j,k)=t(i03+ixy)
       QQ(i,j,k)=QVAPOR(i03+ixy)
       uu(i,j,k)=u(i03+ixy)
       vv(i,j,k)=v(i03+ixy)
       water(i,j,k)=QVAPOR(i03+ixy)*1.E06*29./18.!ppm
       if(k==1) height1(i,j,k) = 2*(heiz(i03+ixy)-HGT1(i0+ixy))     !! layer interface height, change h to heiz, by chenhs
       if(k/=1) height1(i,j,k) = 2*(heiz(i03+ixy)-HGT1(i0+ixy))-height1(i,j,k-1)
      enddo    !k     
      
        CALL getland(LAND_USE(i0+ixy),land(i,j))      
        
       tsurf(i,j)= T2(i0+ixy)
       xlat(i,j) = LATITCRS(i0+ixy)
       xlon(i,j) = LONGICRS(i0+ixy)
       SWDOWN1(i,j)=SWDOWN(i0+ixy)
       pblhgt(i,j) = PBL_HGT(i0+ixy)
       if(FSNOW(i0+ixy).ge.0.001) then !! m
         icdsno(i,j)=1  !! snow
       else  
         icdsno(i,j)=0
       endif
   enddo !j
  enddo !i
  
    lrddrt=.false.;icddrt=0  !if droughtness index
    lrdsno=.true.            !if snow index, by chenhs
     
!!!! ***    GAS SPECIES ****
  do ig=1,igasCBM+2  !! 2 for HG0 and HG2, by chenhs 
   CALL drydep_gas(myid,imonth2,tsurf,xlat,xlon,QQ,QQ,height1,ppp,uu,vv,&
               SWDOWN1,land,water,ttn,pblhgt,lrddrt,icddrt,lrdsno,icdsno,vdep,&
               sx(ne),ex(ne),sy(ne),ey(ne),nzz,ig)
   if(ig.eq.(igasCBM+1)) then
     i02Gas = ip2memGas(103,ne)
   else if(ig.eq.(igasCBM+2)) then
     i02Gas = ip2memGas(104,ne)
   else
     i02Gas = ip2memGas(ig,ne)
   endif                    
   do j=sy(ne)-1,ey(ne)+1
    do  i=sx(ne)-1,ex(ne)+1
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      !if(ig.eq.(igasCBM+1)) then
      !  DryVelGas(i02Gas+ixy) = vdep(i,j)*0.8 !! if no natural emission and re-emission, set Vd-HG0 to zero
      !else                         !! this is the typical option for regional modeling with only anthropogenic emissions
        if(ig.eq.11) then  !! o3 on snow cover
          i0=ip2mem(ne)
          if(FSNOW(i0+ixy).ge.0.001) then !! m
           DryVelGas(i02Gas+ixy) = vdep(i,j)/3.0
           if(imonth2.ge.6.and.imonth2.le.10.and.xlat(i,j).le.-60) then
           DryVelGas(i02Gas+ixy) = vdep(i,j)/2.0
           endif
          else
           DryVelGas(i02Gas+ixy) = vdep(i,j)*1.1
          endif
        else
          DryVelGas(i02Gas+ixy) = vdep(i,j)
        endif 
      !endif
    enddo
   enddo    
  enddo !ig
!!! ***    AEROSOL SPECIES ***
  do ig = igasCBM + 1, igas-2  !! igas-2=HGP, as PM2.5, by chenhs

   IF(ig==igasCBM+2) THEN ! PM10
      diam = 1.0E-6*sqrt(2.5*10.0)
   ELSE
      diam = 1.0E-6*sqrt(0.1*2.5)  ! THINK OTHER AEROSOLS AS PM25
   ENDIF
   
   CALL drydep_aer(myid,imonth2,tsurf,xlat,xlon,height1,ppp,uu,vv,&
                    land,ttn,pblhgt,diam,vdep,sx(ne),ex(ne),sy(ne),ey(ne),nzz,ig)
   
   if(ig.eq.(igas-2)) then
     i02Gas = ip2memGas(105,ne)
   else                 
     i02Gas = ip2memGas(ig,ne)
   endif
    
   do j=sy(ne)-1,ey(ne)+1
    do  i=sx(ne)-1,ex(ne)+1
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     DryVelGas(i02Gas+ixy) = vdep(i,j)
    enddo    
   enddo 

  enddo !ig

! *** FOR DUST AND SEA SALT ****
  DO IA=1,IAER
  DO IS=1,ISIZE
    IF(IA==1) THEN
      DIAM = 1.0E-6*SQRT(MSIZDIS(IS+1)*MSIZDIS(IS))
    ELSE IF(IA==2) THEN 
      DIAM = 1.0E-6*SQRT(MSIZDID(IS+1)*MSIZDID(IS))
    ENDIF
    CALL drydep_dust_salt(myid,imonth2,tsurf,xlat,xlon,height1,ppp,uu,vv,&
                    land,ttn,pblhgt,diam,vdep,sx(ne),ex(ne),sy(ne),ey(ne),nzz,IA)
    I03AER = IP3MEMAER(IS,IA,NE)
    do j=sy(ne)-1,ey(ne)+1
    do  i=sx(ne)-1,ex(ne)+1
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     DryVeldust(I03AER+ixy) = vdep(i,j)*0.45
    enddo
    enddo
  ENDDO ! IS
  ENDDO ! IA

! ****   
  deallocate(ppp,ttn,land,tsurf,xlat,xlon,QQ,height1,uu,vv,water,SWDOWN1,vdep,pblhgt,icddrt,icdsno)
 endif    !idry    

   do ig=1,igas !Gas
      i04 = ip4mem(1,ig,ne)
      i02Gas = ip2memGas(ig,ne)
      i03=ip3mem(1,ne)
      call dry_dep( myid, gas(i04),DryVelGas(i02gas),DRYDEP2(i02gas),dz(i03), &
                    t(i03),Plev(i03), sx(ne), ex(ne), sy(ne), ey(ne), dtstep(ne), ig, igasCBM, igas) !! change dt,chenhs
   enddo !ig

! *** FOR DUST AND SEA SALT ****
    DO IA=1,IAER
    DO IS=1,ISIZE   
      I03AER = IP3MEMAER(IS,IA,NE)       
      I05=IP5MEM(1,IS,IA,NE)
      i03=ip3mem(1,ne)
      I02 = IP2MEM(NE)


      IF(IA==2) THEN
         CALL DUSTDRYDEP (MYID, DUSTDRY(I02), AER(I05),DryVeldust(I03AER), dz(i03),sx(ne), ex(ne),&
                        sy(ne), ey(ne), dtstep(ne))  ! to calculate the DUST DRY DEPOSITION  UG/M2/HOUR
        if(ifdustcom.eq.1) then
        do iduc = 1, ndustcom
         i05c = ip5memc (1,is,iduc,ne)

         IF(iduc==7) then
          CALL DUSTDRYDEP(MYID, DUSTDRYSO4(I02),DUSTCOMP(I05C),DryVeldust(I03AER), dz(i03),sx(ne), ex(ne),&
                          sy(ne), ey(ne), dtstep(ne))  !  to calculate the DUST DSO4 DRY DEPOSITION  UG/M2/HOUR
         ELSE IF(iduc == 8) then                       ! change dt,chenhs
          CALL DUSTDRYDEP(MYID, DUSTDRYNO3(I02),DUSTCOMP(I05C),DryVeldust(I03AER), dz(i03),sx(ne), ex(ne),&   
                          sy(ne), ey(ne), dtstep(ne))  !  to calculate the DUST DNO3 DRY DEPOSITION  UG/M2/HOUR
         ELSE IF (iduc==9) then                        ! change dt,chenhs
          CALL DUSTDRYDEP(MYID, DUSTDRYFeII(I02),DUSTCOMP(I05C),DryVeldust(I03AER), dz(i03),sx(ne), ex(ne),&
                          sy(ne), ey(ne), dtstep(ne))  !  to calculate the DUST DFeII DRY DEPOSITION  UG/M2/HOUR
         ELSE IF(iduc==10) then                        ! change dt,chenhs
          CALL DUSTDRYDEP(MYID, DUSTDRYFeIII(I02),DUSTCOMP(I05C),DryVeldust(I03AER), dz(i03),sx(ne), ex(ne),&
                          sy(ne), ey(ne), dtstep(ne))  !  to calculate the DUST DFeIII DRY DEPOSITION  UG/M2/HOUR
         ENDIF                                         ! change dt,chenhs
        enddo
        endif

      ENDIF

      call dry_dep_gas_zhu( myid, AER(i05),DryVeldust(I03AER),dz(i03), &
            sx(ne), ex(ne), sy(ne), ey(ne), dtstep(ne))  !! change dt,chenhs

      IF(IA==1)  THEN !SEA SALT
       if(ifseacom.eq.1) then 
       do iduc = 1,nseacom
        i05c = ip5memcs (1,is,iduc,ne)
        CALL  dry_dep_gas_zhu( myid, SEACOMP(i05c),DryVeldust(I03AER),dz(i03), &
            sx(ne), ex(ne), sy(ne), ey(ne), dtstep(ne))   !! change dt,chenhs
       enddo ! iduc
       endif ! ifseacom
      ELSE IF(IA==2) THEN ! DUST
       if(ifdustcom.eq.1) then
       do iduc = 1,ndustcom
        i05c = ip5memc (1,is,iduc,ne)
        CALL  dry_dep_gas_zhu( myid, DUSTCOMP(i05c),DryVeldust(I03AER),dz(i03), &
            sx(ne), ex(ne), sy(ne), ey(ne), dtstep(ne))   !! change dt,chenhs
       enddo ! iduc       
       endif ! ifdustcom
      ENDIF

    ENDDO ! IS
    ENDDO ! IA  

!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
if(ifprocess.eq.1) then
   ! iprocess
     do k=1,nzz
         do ig=1,iPrintTermGas    ! for gas phase
             igg=IGGPOS(ig)
             i04=ip4mem(k,igg,ne)
             i05=ipGasTermBal(k,7,ig,ne) ! 7--> dry dep
     call termbal(myid,gasOLD(i04),gas(i04),GasTermBal(i05), &
              k,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dtstep(ne))   !! change dt,chenhs
         enddo
     enddo
endif
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++


 735  continue
!if(myid.eq.0) print *, "TEST08-CBMZ"
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccc     Gas  Chemistry with CBM-Z
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
iitime = time(ne)
call getnewdate(iyear1,imonth1,idate1,ihour1,iitime, &
                iyear2,imonth2,iday2,ihour2,iminute2 )
call JDAY(iday2,imonth2,iyear2, juday, iseason)
i02 = ip2mem(ne)
!!!!!!!!!!!!!!!!!!!!!!! make balance !!!!!!!!!!!!!!!!!!!!!!!!
if(ifbalance.eq.1.and.ne.eq.1) then
  do ig=103,105
     amount_bf(myid+1,ig)=0.0
     do k=1,nzz
        I03=IP3MEM(K,NE)
        I04=IP4MEM(K,IG,NE)
        call N_Balance(MYID,gas(I04),amount_bf(myid+1,ig),DX(I03),DY(I03),DZ(I03),SX(NE),EX(NE),SY(NE),EY(NE))
     enddo
  enddo

  CDNUM='c000'
  WRITE(CDNUM(2:4),'(I3.3)')MYID
  OPEN(10,FILE='out/Amount.bf.dat'//CDNUM,FORM='FORMATTED')
  WRITE(10,*) MYID,amount_bf(myid+1,103),amount_bf(myid+1,104),amount_bf(myid+1,105)
  close(10)
  
  IF(NUMPROCS>1) CALL mpi_barrier( local_com,IERR )
   total_bf=total_af
   total_af=0.
   DO I=0,NUMPROCS-1
      CDNUM='c000'
      total_tmp=0.
      WRITE(CDNUM(2:4),'(I3.3)')I
      OPEN(10,FILE='Amount.bf.dat'//CDNUM,FORM='FORMATTED')
      read(10,*) iiii,total_tmp(1),total_tmp(2),total_tmp(3)
      CLOSE(10)
      total_af=total_af+sum(total_tmp(1:3))
   ENDDO

  if(myid.eq.0) then
   acdry=acdry+(total_bf-total_af)
   open(11,file='real_balance.txt',position='append')
   write(11,'(a60,f16.2,2x,f16.6)') "After DRYDEP,Global=",total_af,(total_af-total_bf)/total_bf*100
   if(mod(ihour2,dtout(ne))==0.and.itt.eq.12) then
     write(11,'(a60,f16.4)') "acumulated drydep=",acdry
   endif
   close(11)
  endif
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------calculate the begin time from noon March 21 (UTC)---------------
!------------ for 2004 is 81,but for other is 80,because the eb is 29 in 2004
    IF(juday>=81.) THEN
      IF(ihour2>= 12) then
        tbeg_dd = juday - 81              !day
        tbeg_hh = ihour2- 12              !hour
        tbeg_mm = iminute2                !min   
        tbeg_ss = 0.0                     !second
      ELSE
        tbeg_dd = juday - 81 -1
        tbeg_hh = ihour2- 12 +24
        tbeg_mm = iminute2
        tbeg_ss = 0.0
      ENDIF
    ELSE
       IF(ihour2> 12) then
        tbeg_dd = juday - 81 +1
        tbeg_hh = ihour2- 12-24
        tbeg_mm = iminute2
        tbeg_ss = 0.0
      ELSE
        tbeg_dd = juday - 81
        tbeg_hh = ihour2- 12
        tbeg_mm = iminute2
        tbeg_ss = 0.0
      ENDIF
    ENDIF
    !!!! ++++++++++++++ for Hg to judge day and night, by chenhs ++++++++++++ !!!!
    if(ifHg.eq.1) then
      tmid_sec00  = ((tbeg_dd*24.+tbeg_hh)*60.+tbeg_mm)*60.+tbeg_ss
    endif
    !!!! ++++++++++++++        end Hg   +++++++++++++++++++++++++++++++++++++ !!!!
!------------------------------------------------------------------------
if(ichemgas == 1)then  !! if do CBM-Z gas chemistry
  !!!!!! chenhs to save time !!!!!!
  if(dtstep(ne).ge.1800.) then
     dttmp=3600
     dttmp0=60 
  else
     dttmp=1200
     dttmp0=20
  endif 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(mod(iitime,int(dttmp)) == 0 )then  !! need to change if the timestep is changed, chenhs
! if(myid.eq.0) print *,"chenhs-iitime=",iitime       
!cccccccccccccccccccccccccccc   set the run time   cccccccccccccccc
      trun_dd = 0.0       ! days
      trun_hh = 0.0       ! hours
      trun_ss = 0.0       !seconds
      trun_mm = dttmp0      ! minutes 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       dt_min = dttmp0      ! the dt of chemcial 
       msolar = 1         !(flag) 1 = diurnally varying photolysis; 2 = fixed phot
       mphoto = 2         !(flag) 1 = Old parameterization; 2 = New parameterization
       iprint = 0         !integer = freq of output. 
                          !Every iprint*dt_min mins.as the original,
                          !but for calculate ,it's not meaning
  !!!  emit_ratio !!!!
  allocate(ratioemitAnt(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz,igas))
  allocate(ratioemitShp(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz,igas))
  allocate(ratioemitBB(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz,igas))
  allocate(ratioemitBio(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz,igas))
  allocate(ratioemitOce(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz,igas))
  allocate(ratioemitSoi(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz,igas))
  allocate(ratioemitLig(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz,igas))
  allocate(ratioemitAirc(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz,igas))
    ratioemitAnt = 0.0
    ratioemitShp = 0.0
    ratioemitBB  = 0.0
    ratioemitBio = 0.0
    ratioemitOce = 0.0
    ratioemitSoi = 0.0
    ratioemitLig = 0.0
    ratioemitAirc = 1.0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(ifseacom.eq.1)  allocate(FSO4_SSA(ISIZE),FNO3_SSA(ISIZE))
  if(ifdustcom.eq.1) allocate(FSO4_DUST(ISIZE),FNO3_DUST(ISIZE))

! Zifa 2006/10/21
Gas_Chemistry:   do j = sy(ne),ey(ne)
             do i = sx(ne),ex(ne)
                  !i02=ip2mem(ne)
                  ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

       do k=1,nzz 
            do ig=1,igasCBM
               i04 = ip4mem(k,ig,ne)
               i02Gas= ip2memGas(ig,ne)
               i03 = ip3mem(k,ne)
               cnn(ig)= MAX(gas(i04+ixy) , 1.E-20)
               cnnzifa(ig)= gas(i04+ixy)
               !ppbv the initial concentration as the box chemical model
               species(ig)=GC_NAME(ig)
               !!!!  to get emission ratio !!!!!!
               if(k.eq.1) then
                  call get_ratio_emit(myid,nzz,ig,LAND_USE(i02+ixy),latitcrs(i02+ixy),ratioemitAnt(i,j,:,ig),&
                                      ratioemitBB(i,j,:,ig),ratioemitLig(i,j,:,ig))
                  ratioemitShp(i,j,1,ig)= 1.0
                  ratioemitBio(i,j,1,ig)= 1.0
                  ratioemitOce(i,j,1,ig)= 1.0
                  ratioemitSoi(i,j,1,ig)= 1.0
               endif
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          emission(ig)=((EmtAntGas(i02Gas+ixy)*ratioemitAnt(i,j,k,ig)+ &
                         EmtShpGas(i02Gas+ixy)*ratioemitShp(i,j,k,ig)+&
                         EmtBBGas(i02Gas+ixy)*ratioemitBB(i,j,k,ig)+&
                         EmtBioGas(i02Gas+ixy)*ratioemitBio(i,j,k,ig)+&
                         EmtOceGas(i02Gas+ixy)*ratioemitOce(i,j,k,ig)+&
                         EmtSoiGas(i02Gas+ixy)*ratioemitSoi(i,j,k,ig)+&
                         EmtLigGas(i02Gas+ixy)*ratioemitLig(i,j,k,ig)+&
                         EmtAircGas(i04+ixy)*ratioemitAirc(i,j,k,ig))*3600./ &
          dz(i03+ixy))*(0.08206*t(i03+ixy))/(Plev(i03+ixy)/1000.)/GC_MOLWT(ig)

                           !ug/s/m2-->ppbv/h
    !!!!!!!!!!!!!!!
    ! for Source Mark
    if(ifsm(ne)==1)then
    letdoit=0
    do idm=1,idmSet
       if(igMark(idm)==ig)letdoit=idm
    enddo
    if(letdoit>0)then
       contem0(letdoit)=gas(i04+ixy)
    endif
    endif
    !!!!!!!!!!!!!!!

            enddo  !ig
            
   !zifa for real-time forecast 
   ! to save time
   if((cnn(5)+cnn(6)).ge.0.0 ) then !NO+NO2 > 0.0
          rlon   =  longicrs(i02+ixy)   !rlat the box lat
          rlat   =  latitcrs(i02+ixy)   ! as the above ,but lon 
          zalt_m =  heiz(i03+ixy)       ! the altitude(asl) of box(m)
          RH     =  rh1(i03+ixy)        ! the RH
          te     =  t(i03+ixy)          ! the temp 
          pr_atm =  PA2ATM*Plev(i03+ixy)*100.  ! the pressure but as the atm
          TER    =  HGT1(i02+ixy)     
      ! for heterogeneous chemistry
          if(ifseacom.eq.1.AND.ifdustcom.eq.1) then 
          i04_bc = ip4mem(k,77,ne)
          I05_1=IP5MEM(K,1,2,NE)
          I05_2=IP5MEM(K,2,2,NE)
          I05_3=IP5MEM(K,3,2,NE)
          I05_4=IP5MEM(K,4,2,NE)

          DUST01= AER(I05_1+IXY)
          DUST02= AER(I05_2+IXY)
          DUST03= AER(I05_3+IXY)
          DUST04= AER(I05_4+IXY)

          I05_1=IP5MEM(K,1,1,NE)
          I05_2=IP5MEM(K,2,1,NE)
          I05_3=IP5MEM(K,3,1,NE)
          I05_4=IP5MEM(K,4,1,NE)

          SEA01 = AER (I05_1+IXY)
          SEA02 = AER (I05_2+IXY)
          SEA03 = AER (I05_3+IXY)
          SEA04 = AER (I05_4+IXY)

          PSO4 = AMAX1(ASO4(i02+ixy), 1.E-20)          
          PBC  = AMAX1(gas(i04_bc+ixy), 1.E-20)      
           SSA_SO4=0.;SSA_NO3=0.;DUST_SO4 = 0.;DUST_NO3=0.
          DO IS = 1, ISIZE
            i05c  =  ip5memcs (k,is,7,ne) ! SO4 on Ses Salt
            i05c1 =  ip5memcs (k,is,8,ne) ! NO3 on Sea Salt
            SSA_SO4 =  SSA_SO4+  SEACOMP(i05c+ixy)
            SSA_NO3 =  SSA_NO3+ SEACOMP(i05c1+ixy)
            i05c2 =  ip5memc  (k,is,7,ne) ! SO4 on DUST
            i05c3 =  ip5memc  (k,is,8,ne) ! NO3 on DUST
            DUST_SO4 = DUST_SO4 + DUSTCOMP(i05c2+ixy)
            DUST_NO3 = DUST_NO3 + DUSTCOMP(i05c3+ixy)
          ENDDO

          CALL HETERO_REAC( PSO4, PBC, DUST01,DUST02,DUST03,DUST04,&
                            SEA01,SEA02,SEA03,SEA04, &
                            DUST_SO4, DUST_NO3, SSA_SO4,SSA_NO3,&
                            FSO4_DUST,FNO3_DUST,FSO4_SSA,FNO3_SSA) ! CALCULATE THE RK_HET FOR HETEROGENEOUS CHEMISTRY 
          endif ! ifseacom and ifdustcom

         FCLD = 1.0
         
         if(ifdustcom.eq.1) then
          cnn(kdso4)=0.0  ! set intial value of dso4(heterogeneous) to 0.0 
          cnn(kdno3)=0.0
          RK_HETSO2_DUST  ( i03+ixy) = RK_HET(19) ! HETEROGENEOUS SO2 on dust
          RK_HETHNO3_DUST ( i03+ixy) = RK_HET(12) ! HETEROGENEOUS HNO3 on dust
         endif
         
         call cbmz(cppb,FCLD)  !! need to modify
        
         jo1d(i03+ixy) = rk_photo(jphoto_o3b) 
         jno2(i03+ixy) = rk_photo(jphoto_no2)
         
         
         call chemope(myid, OPE(i03+ixy), trun_mm,i,j,k,ne)

! allocate the so4 and no3 from heterogeneous on dust and ssa 4 bins
         if(ifseacom.eq.1.AND.ifdustcom.eq.1) then
         do is = 1, isize
          i05c  =  ip5memcs (k,is,7,ne) ! SO4 on Ses Salt
          i05c1 =  ip5memcs (k,is,8,ne) ! NO3 on Sea Salt
          i05c2 =  ip5memc  (k,is,7,ne) ! SO4 on DUST
          i05c3 =  ip5memc  (k,is,8,ne) ! NO3 on DUST
          

         CALL ALLD  (MYID, cppb(kdso4),cppb(kdno3), SEACOMP(i05c+ixy),SEACOMP(i05c1+ixy),&
                     DUSTCOMP(i05c2+ixy),DUSTCOMP(i05c3+ixy),&
                     FSO4_SSA(is),FNO3_SSA(is),FSO4_DUST(is),FNO3_DUST(is) )
         enddo ! is  
         endif ! ifseacom,ifdustcom
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
if(ifprocess.eq.1) then
        do ig=1,iPrintTermGas    ! for gas phase
          igg=IGGPOS(ig)
          if(igg==11) then  !for ozone
          i0517=ipGasTermBal(k,17,ig,ne) ! 17: production 11: ozone
          i0518=ipGasTermBal(k,18,ig,ne) ! 18: loss       11: ozone 
          i0519=ipGasTermBal(k,19,ig,ne) ! 19: radicals LOSS due to NOx 11:ozone
          i0520=ipGasTermBal(k,20,ig,ne) ! 20: radicals LOSS due to VOCs 11:ozone        

           call chemprodloss(myid,delta1,delta2,delta3,delta4,trun_mm,i,j,k,ne)
            GasTermBal(i0517+ixy)=delta1+GasTermBal(i0517+ixy)
            GasTermBal(i0518+ixy)=delta2+GasTermBal(i0518+ixy)
            GasTermBal(i0519+ixy)=delta3+GasTermBal(i0519+ixy)
            GasTermBal(i0520+ixy)=delta4+GasTermBal(i0520+ixy)

           endif !igg
        enddo
endif
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
                         
          !  to judge the change ratio
        iwrongchem=0
            ib=11  ! o3 
            deltcppb=abs(cppb(ib)-cnnzifa(ib))
            if(deltcppb.ge.100)iwrongchem=1 ! > 50ppb 
            ib=18  ! so2 
            deltcppb=abs(cppb(ib)-cnnzifa(ib))/(cnnzifa(ib)+0.00001)
            if(deltcppb.ge.0.10)iwrongchem=1 ! >  10% decrease
        if(iwrongchem==1)goto 4190  ! not to use the output of CBZ

          do ig=1,igasCBM
            i04 = ip4mem(k,ig,ne)
              gas(i04+ixy) = cppb(ig)   !cppb the output ppbv 

              call chemprod(ig,i,j,k,delta,trun_mm,ne) !lijie to get the
!                                                   ozone  production 
              
    !!!!!!!!!!!!!!!
    ! for Source Mark
    if(ifsm(ne)==1)then
    letdoit=0
    do idm=1,idmSet
       if(igMark(idm)==ig)letdoit=idm
    enddo
    if(letdoit>0)then
       i02 = ip2mem(ne)
       MapS=int(MapSource(i02+ixy))
       IHgtLev= iHgtLMax-1  ! for chemistry reaction part Zifa/2007/03/03
                   ! this is need to be defined outside for future   
       i04 = ip4mem(k,ig,ne)
       OrgCon     = contem0(letdoit)
       if(ig==11) then
         DeltSpeMark=delta
       else  
         DeltSpeMark= gas(i04+ixy)-OrgCon
       endif
         if(DeltSpeMark > 0. )then
            do ism=1,ismMax
                i04sm=ipSMmem(k,ism,letdoit,ne)
                TmpSM(ism)=SourceMark(i04sm+ixy)
            enddo
            call GetSMarkChange(DeltSpeMark,OrgCon,TmpSM,MapS, &
                 ISrcDefined,IHgtLev,ismMax)
            do ism=1,ismMax
               i04sm=ipSMmem(k,ism,letdoit,ne)
               SourceMark(i04sm+ixy)=TmpSM(ism)
            enddo
         endif
    endif
    endif
    !!!!!!!!!!!!!!!
          enddo !ig                    
 4190     continue
        !!!!++++++++++ Hg gas chemistry, by chenhs ++++++++++++++++++++++!!!!
        if(ifHg.eq.1) then 
         i04_HG0 = ip4mem(k,103,ne)
         i04_HG2 = ip4mem(k,104,ne)
         i04_O3  = ip4mem(k,11,ne)
         i04_H2O2= ip4mem(k,16,ne)
         i04_OH  = ip4mem(k,14,ne)
         CHG0   = gas(i04_HG0+ixy) 
         CHG2   = gas(i04_HG2+ixy) 
         CO3    = gas(i04_O3+ixy) 
         CH2O2  = gas(i04_H2O2+ixy)
         COH    = gas(i04_OH+ixy)
         dttmp=trun_mm*60 
         call HGGASCHEM ( CHG0, CHG2, CO3, CH2O2, COH, t(i03+ixy) , Plev(i03+ixy), latitcrs(i02+ixy),&
                          longicrs(i02+ixy), heiz(i03+ixy), landmask(i02+ixy), tmid_sec00, dttmp )  
        
         gas(i04_HG0+ixy) = CHG0
         gas(i04_HG2+ixy) = CHG2 
        endif
        !!!!++++++++++    end Hg gas chemistry     ++++++++++++++++++++++!!!!   
    endif  !NO+NO2 >= 0.0
        enddo !k
             enddo !i
                 enddo Gas_Chemistry  !j

     deallocate(ratioemitAnt,ratioemitShp,ratioemitBB,ratioemitBio,&
                ratioemitOce,ratioemitSoi,ratioemitLig,ratioemitAirc) 
     if(ifseacom.eq.1)  deallocate(FSO4_SSA,FNO3_SSA)
     if(ifdustcom.eq.1) deallocate(FSO4_DUST,FNO3_DUST)
endif !! iitime
endif !! ichemgas
!if(myid.eq.0) print *, "TEST0801-END_CBMZ"
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
if(ifprocess.eq.1) then
  ! iprocess
    do k=1,nzz
        do ig=1,iPrintTermGas    ! for gas phase
            igg=IGGPOS(ig)
            i04=ip4mem(k,igg,ne)
            i05=ipGasTermBal(k,6,ig,ne)  ! 6 for chemistry
    call termbal(myid,gasOLD(i04),gas(i04),GasTermBal(i05), &
            k,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dtstep(ne))   !! change dt,chenhs
        enddo
    enddo
endif
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++

200 continue

! ccccc   Fe evolution in atmosphere by lijie  ccccc
   if(ifdustcom.eq.1) then 
     do k = 1, nzz-1
       i03=ip3mem(k,ne)
       i04_so2=ip4mem(k,18,ne) ! so2
       i04_hno3=ip4mem(k,2,ne) ! hno3
       i02=ip2mem(ne)
      do is = 1, isize
       i05   = ip5memc (k,is,9, ne) ! FeII
       i05_1 = ip5memc (k,is,10,ne) ! Total FeIII
       i05_2 = ip5memc (k,is,11,ne) ! coated FeIII
         
      call  FEEVOLUTION(MYID,DUSTCOMP(I05_2),DUSTCOMP(I05_1), DUSTCOMP(I05), gas(i04_so2), gas(i04_hno3),&
                RK_HETSO2_DUST(i03),RK_HETHNO3_DUST(i03), clflo(i02), clfmi(i02), clfhi(i02),&
                SWDOWN(i02),RH1(i03),dtstep(ne), sx(ne), ex(ne), sy(ne), ey(ne),k,is)  !! change dt,chenhs
      enddo

     enddo
   endif ! ifdustcom
! ccc
!CCCCCCCCCCCC    AQUEOUS CHEMISTRY BY LI JIE FROM CAMX CCCCCCC
  
 i02 = ip2mem(ne)
 DO K = 1, nz(ne)
    i03 = ip3mem(k,ne)

    i04_SO2=ip4mem(k,18,ne)
    i04_HNO3=ip4mem(k,2,ne)
    i04_NXOY=ip4mem(k,8,ne) ! N2O5
    i04_NH3=ip4mem(k,4,ne) 
    i04_H2O2=ip4mem(k,16,ne)
    i04_O3=ip4mem(k,11,ne)
    i04_FOA=ip4mem(k,23,ne) ! formic acid      
    i04_H2SO4=ip4mem(k,1,ne)
    !!+++ for Hg ++++!!
    i04_HG0 = ip4mem(k,103,ne)
    i04_HG2 = ip4mem(k,104,ne)
    i04_OH  = ip4mem(k,14,ne)
    i04_HO2 = ip4mem(k,15,ne)
    !!+++++++++++++++!!
  
    DO J = SY(NE), EY(NE)
     DO I = SX(NE),EX(NE)

      ixy = (ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1

      CGAS ( 1) = gas(i04_SO2+IXY)      ! SO2
      CGAS ( 2) = gas(i04_HNO3+IXY)     ! hno3
      CGAS ( 3) = gas(i04_NXOY+IXY)*2.  ! NXOY
      CGAS ( 4) = 330.*1.e03            ! CO2 330 PPM
      CGAS ( 5) = gas(i04_NH3+IXY)      ! NH3
      CGAS ( 6) = gas(i04_H2O2+IXY)     ! H2O2
      CGAS ( 7) = gas(i04_O3+IXY)       ! O3
      CGAS ( 8) = gas(i04_FOA+IXY)      ! HCOOH
      CGAS ( 9) = 1.E-3                 ! MNP
      CGAS (10)= 1.E-3                  ! PAA
      CGAS (11)= gas(i04_H2SO4+IXY)     ! H2SO4(G)


      CAER ( 1 ) = ASO4(I03+IXY) ! ASO4
      CAER ( 2 ) = ANH4(I03+IXY) ! ANH4
      CAER ( 3 ) = ANO3(I03+IXY) ! ANO3
      CAER ( 4 ) = ASO4(I03+IXY)/1.5 ! ACACO3
      CAER ( 5 ) = ASO4(I03+IXY)/12. ! AMGCO3
      CAER ( 6 ) = ACL(I03+IXY)      ! ANACL
      CAER ( 7 ) = 0.01              ! FE+++
      CAER ( 8 ) = 0.005             ! Mn++
      CAER ( 9 ) = 0.000             ! potcl

      CLW_TMP = CLW(i03+ixy)*29.*44.9*(273./t(i03+ixy))*(Plev(i03+ixy)/1013.)  !! FROM KG/KG TO G/M3
      if(t(i03+ixy).LT.273) then
         CLW_TMP = AMAX1 (0., CLW_TMP*(t(i03+ixy)-243.)/(273.-243.))
      endif
       
      IF(CLW_TMP.GE.0.05.AND.t(i03+ixy).GE.243.) THEN  !! limit to call aqueous

      CALL AQUEOUS(t(i03+ixy), Plev(i03+ixy), QVAPOR(i03+ixy), CLW_TMP, CGAS, CAER, CPH(i03+ixy),&
                   dtstep(ne), sx(ne), ex(ne), sy(ne),ey(ne),I,J,K )          !! change dt,chenhs

      gas(i04_SO2+IXY)   = CGAS (1)    ! SO2
      gas(i04_HNO3+IXY)  = CGAS (2)    ! HNO3
      gas(i04_NXOY+IXY)  = CGAS (3)/2. ! N2O5
      gas(i04_NH3+IXY)   = CGAS (5)    ! NH3
      gas(i04_H2O2+IXY)  = CGAS (6)    ! H2O2
      gas(i04_O3+IXY)    = CGAS (7)    ! O3
      gas(i04_FOA+IXY)   = CGAS (8)    ! HCOOH
      gas(i04_H2SO4+IXY) = CGAS (11)   ! H2SO4 
      !!!!+++++++++++++++++++++++++ Hg Aqueous Chemistry, by chenhs ++++++++++++++++++++++++!!!!
      if(ifHg.eq.1) then
        CHG0   = gas(i04_HG0+ixy)
        CHG2   = gas(i04_HG2+ixy)
        CO3    = gas(i04_O3+ixy)
        CSO2   = gas(i04_SO2+IXY)    
        COH    = gas(i04_OH+ixy)
        CHO2   = gas(i04_HO2+ixy)
        CPM10  = 0.
        do ig=igasCBM+1,igas-4
           i04 = ip4mem(k,ig,ne)
           if(ig.ne.79.and.gas(i04+ixy).gt.0.) then
              CPM10 = CPM10 + gas(i04+ixy)  !! all PM conc
           endif
        enddo
 
        call hgaqschem( CHG0, CHG2, CO3, CSO2, COH, CHO2, CPM10,CPH(i03+ixy), CLW_TMP, t(i03+ixy), Plev(i03+ixy),& 
                        latitcrs(i02+ixy), longicrs(i02+ixy), heiz(i03+ixy), landmask(i02+ixy), tmid_sec00, dtstep(ne))

        gas(i04_HG0+ixy) = CHG0
        gas(i04_HG2+ixy) = CHG2 
      endif
      !!!!+++++++++++++++++++++++++    end Hg Aqueous Chemistry     ++++++++++++++++++++++++!!!!
      ELSE
      CPH(I03+IXY) = -1.E20   
      ENDIF

     ENDDO !I
    ENDDO !J


ENDDO !K


! CCCCCCCCCCCCCCCCCCCCC   END AQUEOUS   CCCCCCCCCCCCCCCCCCCCCCCCCCCC
!if(myid.eq.0) print *, "TEST0802-END_AQUEOUS"
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
if(ifprocess.eq.1) then
 ! iprocess
   do k=1,nzz
     do ig=1,iPrintTermGas    ! for AQUEOUS  process
       igg=IGGPOS(ig)
       i04=ip4mem(k,igg,ne)
       i05=ipGasTermBal(k,24,ig,ne)  ! 24 for AQUEOUS CHEMISTRY
    call termbal(myid,gasOLD(i04),gas(i04),GasTermBal(i05), &
            k,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dtstep(ne))  !! change dt,chenhs
     enddo
   enddo
endif
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++

!************** ISORROPIA AEROSOL THERMODUNAMICS MODEL ***************
!***  TO PARTITION THE SEMIVOLATILE INORGANIC AEROSOL COMPONENTS *****
!***    SODIUM-AMMONIUM-SULFATE-NITRATE-CHLORIDE AEROSOL         *****
!
!
    do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1

!!!
!   WI(1) : TOTAL SODIUM   AS EQUIVALENT NA
!   WI(2) : TOTAL SULFATE  AS EQUIVALENT H2SO4
!   WI(3) : TOTAL AMMONIUM AS EQUIVALENT NH3
!   WI(4) : TOTAL NITRATE  AS EQUIVALENT HNO3
!   WI(5) : TOTAL CHLORIDE AS EQUIVALENT HCL
!   WO :  17 AEROSOLS SPECIES
!   WOG:  4  GASEOUS  SPECIES

        i03 = ip3mem(k,ne)
        TEMPI     =  DBLE(t(i03+ixy ))  ! THE TEMPERATURE IN K
        RHI       =  DBLE(rh1(i03+ixy)) ! THE RH IN 0-100%   

  !!!!!!!!!!!!!!!
    ! for SourceMark
   if(ifsm(ne)==1)then
    do ig=1, igas
      i04 = ip4mem(k,ig,ne)
      
        letdoit=0
       do idm=1,idmSet
        if(igMark(idm)==ig)letdoit=idm
       enddo
       if(letdoit>0)then
        contem0(letdoit)=gas(i04+ixy)
       endif
       
    enddo  !ig
  endif! ifsm     
  !!!!!!!!!!!!!!!

!       
!  !!! CONVERT GAS IN PPBV TO UG/M3
    
        do ig=1,4
         i04 = ip4mem(k,ig,ne) 
          gas(i04+ixy) = gas(i04+ixy)*                  &
            GC_MOLWT(ig)/( (0.08206*T(i03+ixy))/(Plev(i03+ixy)/1000.) )
        enddo !ig 
!
    !!! GET WI(5)
        i04_1 = ip4mem(k,1,ne)   ! H2SO4(G)
        i04_2 = ip4mem(k,2,ne)   ! HNO3(G)
        i04_3 = ip4mem(k,3,ne)   ! HCL(G)
        i04_4 = ip4mem(k,4,ne)   ! NH3(G)
        i04_5 = ip4mem(k,80,ne)  ! Na+(AE)
        i04_6 = ip4mem(k,81,ne)  ! NH4+(AE)
        i04_7 = ip4mem(k,82,ne)  ! CL-(AE)
        i04_8 = ip4mem(k,83,ne)  ! SO4--(AE)
        i04_9 = ip4mem(k,84,ne)  ! HSO4-(AE)
        i04_10= ip4mem(k,85,ne)  ! NO3-(AE)
        i04_11= ip4mem(k,86,ne)  ! NACL(S)
        i04_12= ip4mem(k,87,ne)  ! NA2SO4(S)
        i04_13= ip4mem(k,88,ne)  ! NANO3(S)
        i04_14= ip4mem(k,89,ne)  ! (NH4)2SO4(S)
        i04_15= ip4mem(k,90,ne)  ! NH4NO3(S)
        i04_16= ip4mem(k,91,ne)  ! NH4CL(S)
        i04_17= ip4mem(k,92,ne)  ! H2SO4(AQ)
        i04_18= ip4mem(k,93,ne)  ! NH4HSO4(S)
        i04_19= ip4mem(k,94,ne)  ! NAHSO4(S)
        i04_20= ip4mem(k,95,ne)  ! (NH4)3H(SO4)2(S)
        i04_21= ip4mem(k,79,ne)  !  H+(AE)
        i04_22= ip4mem(k,102,ne) ! AH2O 
        
        WI    = 0.0
        WO    = 0.0        
        WOG   = 0.0
        
        WI(1) = DBLE(gas(i04_5+ixy) * 23.0/GC_MOLWT(80)&
              + gas(i04_11+ixy)* 23.0/GC_MOLWT(86)&
              + gas(i04_12+ixy)* 23.0/GC_MOLWT(87)*2.0 &
              + gas(i04_13+ixy)* 23.0/GC_MOLWT(88)&
              + gas(i04_19+ixy)* 23.0/GC_MOLWT(94))
! 
        WI(2) = DBLE(gas(i04_1+ixy) * 98.0/GC_MOLWT(1) &
              + gas(i04_8+ixy) * 98.0/GC_MOLWT(83)&
              + gas(i04_9+ixy) * 98.0/GC_MOLWT(84)&
              + gas(i04_12+ixy)* 98.0/GC_MOLWT(87)&
              + gas(i04_14+ixy)* 98.0/GC_MOLWT(89)&
              + gas(i04_17+ixy)* 98.0/GC_MOLWT(92)&
              + gas(i04_18+ixy)* 98.0/GC_MOLWT(93)&
              + gas(i04_19+ixy)* 98.0/GC_MOLWT(94)&
              + gas(i04_20+ixy)* 98.0/GC_MOLWT(95)*2.0)
!
        WI(3) = DBLE(gas(i04_4+ixy) * 17.0/GC_MOLWT(4) &
              + gas(i04_6+ixy) * 17.0/GC_MOLWT(81)&
              + gas(i04_14+ixy)* 17.0/GC_MOLWT(89)*2.0 &
              + gas(i04_15+ixy)* 17.0/GC_MOLWT(90)&
              + gas(i04_16+ixy)* 17.0/GC_MOLWT(91)&
              + gas(i04_18+ixy)* 17.0/GC_MOLWT(93)&
              + gas(i04_20+ixy)* 17.0/GC_MOLWT(95)*3.0)
!        
        WI(4) = DBLE(gas(i04_2+ixy) * 63.0/GC_MOLWT(2) &
              + gas(i04_10+ixy)* 63.0/GC_MOLWT(85)&
              + gas(i04_13+ixy)* 63.0/GC_MOLWT(88)&
              + gas(i04_15+ixy)* 63.0/GC_MOLWT(90))
!
        WI(5) = DBLE(gas(i04_3+ixy) * 36.5/GC_MOLWT(3) &
              + gas(i04_11+ixy)* 36.5/GC_MOLWT(86)&
              + gas(i04_16+ixy)* 36.5/GC_MOLWT(91))

        IF(WI(2).LE.0.1) GOTO 2880
        IF(WI(3).LE.0.1) GOTO 2880
        IF(WI(4).LE.0.1) GOTO 2880
                
!*****   TO CALL ISRPINTR ******
        CALL ISRPINTR(MYID,WI,RHI,TEMPI,WO,WOG(1),WOG(2),WOG(3),WOG(4),I,J,K)        

         gas(i04_21+ixy) = SNGL(WO(1))
         gas(i04_5+ixy)  = SNGL(WO(2))
         gas(i04_6+ixy)  = SNGL(WO(3))
         gas(i04_7+ixy)  = SNGL(WO(4))
         gas(i04_8+ixy)  = SNGL(WO(5))
         gas(i04_9+ixy)  = SNGL(WO(6))
         gas(i04_10+ixy) = SNGL(WO(7))
         gas(i04_11+ixy) = SNGL(WO(8))
         gas(i04_12+ixy) = SNGL(WO(9))
         gas(i04_13+ixy) = SNGL(WO(10))
         gas(i04_14+ixy) = SNGL(WO(11))
         gas(i04_15+ixy) = SNGL(WO(12))
         gas(i04_16+ixy) = SNGL(WO(13))
         gas(i04_17+ixy) = SNGL(WO(14))
         gas(i04_18+ixy) = SNGL(WO(15))
         gas(i04_19+ixy) = SNGL(WO(16))
         gas(i04_20+ixy) = SNGL(WO(17))
         gas(i04_22+IXY) = SNGL(AWATER)

         gas(i04_4+ixy)  = SNGL(WOG(1))
         gas(i04_2+ixy)  = SNGL(WOG(2))
         gas(i04_1+ixy)  = SNGL(WOG(3))
         gas(i04_3+ixy)  = SNGL(WOG(4))
!
                                                                                                                        
         ! Jian-Bin Wu change the calculation of these PM
        ANA(I03+IXY) = DBLE(gas(i04_5+ixy) * 23.0/GC_MOLWT(80)&
              + gas(i04_11+ixy)* 23.0/GC_MOLWT(86)&
              + gas(i04_12+ixy)* 23.0/GC_MOLWT(87)*2.0 &
              + gas(i04_13+ixy)* 23.0/GC_MOLWT(88)&
              + gas(i04_19+ixy)* 23.0/GC_MOLWT(94))
        ASO4(I03+IXY)= DBLE(gas(i04_8+ixy) * 96.0/GC_MOLWT(83)&
              + gas(i04_9+ixy) * 96.0/GC_MOLWT(84)&
              + gas(i04_12+ixy)* 96.0/GC_MOLWT(87)&
              + gas(i04_14+ixy)* 96.0/GC_MOLWT(89)&
              + gas(i04_17+ixy)* 96.0/GC_MOLWT(92)&
              + gas(i04_18+ixy)* 96.0/GC_MOLWT(93)&
              + gas(i04_19+ixy)* 96.0/GC_MOLWT(94)&
              + gas(i04_20+ixy)* 96.0/GC_MOLWT(95)*2.0)
        ANH4(I03+IXY)= DBLE(gas(i04_6+ixy) * 18.0/GC_MOLWT(81)&
              + gas(i04_14+ixy)* 18.0/GC_MOLWT(89)*2.0 &
              + gas(i04_15+ixy)* 18.0/GC_MOLWT(90)&
              + gas(i04_16+ixy)* 18.0/GC_MOLWT(91)&
              + gas(i04_18+ixy)* 18.0/GC_MOLWT(93)&
              + gas(i04_20+ixy)* 18.0/GC_MOLWT(95)*3.0)
        ANO3(I03+IXY)= DBLE(gas(i04_10+ixy)* 62.0/GC_MOLWT(85)&
              + gas(i04_13+ixy)* 62.0/GC_MOLWT(88)&
              + gas(i04_15+ixy)* 62.0/GC_MOLWT(90))
        ACL(I03+IXY) = DBLE(gas(i04_7+ixy)* 35.5/GC_MOLWT(82)&
              + gas(i04_11+ixy)* 35.5/GC_MOLWT(86)&
              + gas(i04_16+ixy)* 35.5/GC_MOLWT(91))
        
2880  CONTINUE
   
!!!  CONVERT GAS IN UG/M3 TO PPBV
        do ig=1,4
           i04 = ip4mem(k,ig,ne)
            gas(i04+ixy) = gas(i04+ixy)*                          &
             (0.08206*T(i03+ixy))/(Plev(i03+ixy)/1000.)/GC_MOLWT(ig)
        enddo    !ig   

!!!!!!!!!!!!!!!
  ! for Source Mark
 if(ifsm(ne)==1)then
 DO ig=1,igas
      letdoit=0
    do idm=1,idmSet
       if(igMark(idm)==ig)letdoit=idm
    enddo
    if(letdoit>0)then
       i02 = ip2mem(ne)    
       MapS=int(MapSource(i02+ixy))
       IHgtLev= iHgtLMax-1  ! for chemistry reaction part Zifa/2007/03/03
                   ! this is need to be defined outside for future
       i04 = ip4mem(k,ig,ne)
       OrgCon     = contem0(letdoit)
       DeltSpeMark= gas(i04+ixy)-OrgCon
       if(DeltSpeMark > 0. )then
            do ism=1,ismMax
                i04sm=ipSMmem(k,ism,letdoit,ne)
                TmpSM(ism)=SourceMark(i04sm+ixy)
            enddo
            call GetSMarkChange(DeltSpeMark,OrgCon,TmpSM,MapS, &
                 ISrcDefined,IHgtLev,ismMax)
            do ism=1,ismMax
               i04sm=ipSMmem(k,ism,letdoit,ne)
               SourceMark(i04sm+ixy)=TmpSM(ism)
            enddo
        endif
     endif ! letdoit
  ENDDO !IG
  endif  ! ifsm
!!!!!!

       enddo !k
      enddo  !i
    enddo    !j

!
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
if(ifprocess.eq.1) then
 ! iprocess
   do k=1,nzz
     do ig=1,iPrintTermGas    ! for ISORROPIA 
       igg=IGGPOS(ig)
       i04=ip4mem(k,igg,ne)
       i05=ipGasTermBal(k,22,ig,ne)  ! 22 for inorganic aerosol partioning
    call termbal(myid,gasOLD(i04),gas(i04),GasTermBal(i05), &
            k,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dtstep(ne))  !! change dt,chenhs
     enddo
   enddo
endif
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
!if(myid.eq.0) print *, "TEST0803-END_ISORROPIA"
! CCCCCCCCCCCCCCCCCCCCCCCCCCC   TO WET DEPOSITION  CCCCCCCCCCCCCCCCCC
  !!!!!! chenhs to save time !!!!!!
  if(dtstep(ne).ge.1800.) then
     dttmp=1800
  else
     dttmp=600
  endif 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(mod(iitime,int(dttmp)) == 0 )then  !! need to change if the timestep is changed, chenhs to save time
! if(myid.eq.0) print *,"chenhswet-iitime=",iitime       

allocate( clwc(nzz), rnwc(nzz),twet(nzz),pwet(nzz) , RR(NZZ), VOLRAT(NZZ)) 
allocate (pwr_c(nzz),cwc_c(nzz),con(nzz))
! CCC   WET DEPOSITION
iwb=sx(ne)-1;ieb=ex(ne)+1
jsb=sy(ne)-1;jeb=ey(ne)+1
ALLOCATE(DEPFLD(iwb:ieb,jsb:jeb,IGAS),DEPFLD2(iwb:ieb,jsb:jeb,IGAS))
!CCC  END

    do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
        i02 = ip2mem(ne)
        DO k = 1, nzz
         i03 = ip3mem(k,ne)
         clwc(k) = CLW(i03+ixy)
         rnwc(k) = RNW(i03+ixy)        
         TWET(k) = t(i03+ixy)
         PWET(k) = Plev(i03+ixy)                    
        ENDDO ! k 

! CCC  GET THE LAYER containin precipitation bottom/top
         CALL GETCLOUDDEPTH ( MYID,TWET,PWET,clwc,rnwc,CWC_C, PWR_C,KBOTC,KTOPC,NZZ,&
                 RR,VOLRAT,RAINCON(i02+ixy),RAINNON(i02+ixy),LPREC)
        
         IF(LPREC) THEN  !! if this grid has precipitation
!CCCCCCC   GAS
           DO IG = 1, IGASCBM+2  !! 2 for HG0 and HG2, by chenhs
            if(ig.eq.(igasCBM+1)) then
              i02Gas = ip2memGas(103,ne)
            else if(ig.eq.(igasCBM+2)) then
              i02Gas = ip2memGas(104,ne)
            else
              i02Gas = ip2memGas(ig,ne)
            endif 
             DEPFLD(I,J,IG) =  0.0  !! change by chenhs 
             DEPFLD2(I,J,IG) = 0.0  !! change by chenhs
           DO K = KTOPC, KBOTC,-1
            i03 = ip3mem(k,ne)  
            if(ig.eq.(igasCBM+1)) then      
              i04 = ip4mem(k,103,ne)
            else if(ig.eq.(igasCBM+2)) then
              i04 = ip4mem(k,104,ne)
            else
              i04 = ip4mem(k,ig,ne)
            endif

              con(k) = GAS(I04+IXY)

         
              CALL WETDEP_GAS (MYID, KBOTC, KTOPC, dttmp, DX(I03+IXY),DY(I03+IXY),&  !! change dt,chenhs
                 DZ(I03+IXY),TWET(K), PWET(K),CWC_C(K), PWR_C(K), 0., 0., VOLRAT(K), &
                 RR(K),CPH(I03+IXY),con(k), TMASS, DEPFLD(I,J,IG), DEPFLD2(I,J,IG), &
                 I, J, K, IG ,dtout(ne))
        
     
              GAS(I04+IXY) = con(k)

           ENDDO ! K
            WETDEP(i02gas+ixy)  =  WETDEP(i02gas+ixy)+DEPFLD(I,J,IG)   !! by chenhs
            WETDEP2(i02gas+ixy) =  WETDEP2(i02gas+ixy)+DEPFLD2(I,J,IG) !! by chenhs
           ENDDO ! IG

!CCCCCCC   AEROSOL 

           DO IG = IGASCBM + 1, IGAS-2 !! igas-2=HGP, as PM2.5, by chenhs
               if(ig.eq.(igas-2)) then
                 i02gas = ip2memGas(105,ne)                 
               else
                 i02gas = ip2memGas(ig,ne)
               endif
               DEPFLD(I,J,IG) =  0.0  !! change by chenhs for dt
               DEPFLD2(I,J,IG) = 0.0  !! change by chenhs for dt

            DO K = KTOPC, KBOTC,-1
             i03 = ip3mem(k,ne)
             if(ig.eq.(igas-2)) then
               i04 = ip4mem(k,105,ne)
             else
               i04 = ip4mem(k,ig,ne)
             endif
             con(k) = GAS(I04+IXY)
             IF(ig==igasCBM+2) THEN ! PM10
                diam = 1.0E-6*sqrt(2.5*10.0)
             ELSE
                diam = 1.0E-6*sqrt(0.1*2.5)  ! THINK OTHER AEROSOLS AS PM25
             ENDIF
                         
              CALL WETDEP_AER (MYID, KBOTC, KTOPC, dttmp, DX(I03+IXY), DY(I03+IXY),&  !! change dt,chenhs
                 DZ(I03+IXY),TWET(K), PWET(K), CWC_C(K), PWR_C(K), 0., 0., VOLRAT(K), &  ! juanxiong he
                 RR(K),con(k),diam, TMASS, DEPFLD(I,J,IG), DEPFLD2(I,J,IG), &
                 I, J, K, IG, dtout(ne))
          
               GAS(I04+IXY) = con(k)

             ENDDO ! K

              WETDEP  (i02gas+ixy) =  WETDEP(i02gas+ixy)+DEPFLD(I,J,IG)    !! by chenhs
              WETDEP2 (i02gas+ixy) =  WETDEP2 (i02gas+ixy)+DEPFLD2(I,J,IG) !! by chenhs

            ENDDO ! IG

         !ELSE  !! this grid has no precipitation 
         !  DO IG = 1, IGAS
         !    i02gas = ip2memGas(ig,ne)
         !    WETDEP  (i02gas+ixy) =  WETDEP  (i02gas+ixy)+1.E-20  !! by chenhs
         !    WETDEP2 (i02gas+ixy) =  WETDEP2 (i02gas+ixy)+1.E-20  !! by chenhs
         !  ENDDO
         ENDIF ! LPREC

      enddo 
    enddo

  deallocate(pwet,twet,clwc,rnwc,RR,VOLRAT,pwr_c,cwc_c,con)
  DEALLOCATE(DEPFLD,DEPFLD2)

endif !! iitime, chenhs

! FOR DUST AND SEA SALT ******
 DO K=1,11   ! WET DEP
      i03 = ip3mem(k,ne ) 
     DO IA=1,IAER
     DO IS=1,ISIZE
        I05=IP5MEM(K,IS,IA,NE)
        I02 = IP2MEM(NE)        
        IF( IA == 2 ) THEN
          !!! TO CALCULATE THE DUST WET DEPOSITION UG/M2/HR
          CALL DUSTWETDEP(MYID, DUSTWET(I02),AER(I05),RAINCON(I02),RAINNON(I02),&
                          SX(NE),EX(NE),SY(NE),EY(NE),dtstep(ne),K,dz(i03))  !! change dt,chenhs
          if(ifdustcom.eq.1) then
          do iduc = 1 , ndustcom
            i05c = ip5memc (k,is,iduc,ne)
            IF (iduc ==7 ) then
              CALL DUSTWETDEP(MYID,DUSTWETSO4(I02),DUSTCOMP(I05C),RAINCON(I02),RAINNON(I02),&
                          SX(NE),EX(NE),SY(NE),EY(NE),dtstep(ne),K,dz(i03))  !! change dt,chenhs
            ELSE IF(iduc == 8) then
              CALL DUSTWETDEP(MYID,DUSTWETNO3(I02),DUSTCOMP(I05C),RAINCON(I02),RAINNON(I02),&
                          SX(NE),EX(NE),SY(NE),EY(NE),dtstep(ne),K,dz(i03))          !! change dt,chenhs
            ELSE  IF(iduc==9) then
              CALL DUSTWETDEP(MYID,DUSTWETFeII(I02),DUSTCOMP(I05C),RAINCON(I02),RAINNON(I02),&
                          SX(NE),EX(NE),SY(NE),EY(NE),dtstep(ne),K,dz(i03))          !! change dt,chenhs
            ELSE IF(iduc==10) then
              CALL DUSTWETDEP(MYID,DUSTWETFeIII(I02),DUSTCOMP(I05C),RAINCON(I02),RAINNON(I02),&
                          SX(NE),EX(NE),SY(NE),EY(NE),dtstep(ne),K,dz(i03))          !! change dt,chenhs
            ENDIF
          enddo ! iduc
          endif ! ifdustcom
        ENDIF ! IA IF 
        CALL  WET_DEP_DUST(MYID,AER(I05),RAINCON(I02),RAINNON(I02),&
                          SX(NE),EX(NE),SY(NE),EY(NE),dtstep(ne))  !! change dt,chenhs
        IF(IA==1)  THEN ! SEA SALT
         if(ifseacom.eq.1) then
         do iduc = 1 ,nseacom
          i05c = ip5memcs (k,is,iduc,ne)
          CALL WET_DEP_DUST(MYID,SEACOMP(I05C),RAINCON(I02),RAINNON(I02),&
                          SX(NE),EX(NE),SY(NE),EY(NE),dtstep(ne))  !! change dt,chenhs
         enddo
         endif
        ELSE IF(IA==2) THEN
         if(ifdustcom.eq.1) then
         do iduc = 1 ,ndustcom
          i05c = ip5memc (k,is,iduc,ne)
          CALL WET_DEP_DUST(MYID,DUSTCOMP(I05C),RAINCON(I02),RAINNON(I02),&
                          SX(NE),EX(NE),SY(NE),EY(NE),dtstep(ne))  !! change dt,chenhs
         enddo
         endif
        ENDIF

     ENDDO        !ISIZE
     ENDDO        !IAER

  ENDDO !K

!CCCCCCCCCCCCCCCCC  END CCCC
 !!!!!!!!!!!!!!!!!!!!!!! make balance !!!!!!!!!!!!!!!!!!!!!!!!
if(ifbalance.eq.1.and.ne.eq.1) then
  do ig=103,105
     amount_bf(myid+1,ig)=0.0
     do k=1,nzz
        I03=IP3MEM(K,NE)
        I04=IP4MEM(K,IG,NE)
        call N_Balance(MYID,gas(I04),amount_bf(myid+1,ig),DX(I03),DY(I03),DZ(I03),SX(NE),EX(NE),SY(NE),EY(NE))
     enddo
  enddo

  CDNUM='c000'
  WRITE(CDNUM(2:4),'(I3.3)')MYID
  OPEN(10,FILE='out/Amount.bf.dat'//CDNUM,FORM='FORMATTED')
  WRITE(10,*) MYID,amount_bf(myid+1,103),amount_bf(myid+1,104),amount_bf(myid+1,105)
  close(10)

  IF(NUMPROCS>1) CALL mpi_barrier( local_com,IERR )
   total_bf=total_af
   total_af=0.
   DO I=0,NUMPROCS-1
      CDNUM='c000'
      total_tmp=0.
      WRITE(CDNUM(2:4),'(I3.3)')I
      OPEN(10,FILE='out/Amount.bf.dat'//CDNUM,FORM='FORMATTED')
      read(10,*) iiii,total_tmp(1),total_tmp(2),total_tmp(3)
      CLOSE(10)
      total_af=total_af+sum(total_tmp(1:3))
   ENDDO
  if(myid.eq.0) then
   acwet=acwet+(total_bf-total_af)
   open(11,file='out/real_balance.txt',position='append')
   write(11,'(a60,f16.2,2x,f16.6)') "After WETDEP,Global=",total_af,(total_af-total_bf)/total_bf*100
   if(mod(ihour2,dtout(ne))==0.and.itt.eq.12) then
     write(11,'(a60,f16.4)') "acumulated wetdep=",acwet
   endif
   close(11)
  endif
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
if(ifprocess.eq.1) then
 ! iprocess
   do k=1,nzz
     do ig=1,iPrintTermGas    ! for WETDEP process
       igg=IGGPOS(ig)
       i04=ip4mem(k,igg,ne)
       i05=ipGasTermBal(k,8,ig,ne)  ! 8 for WETDEP
    call termbal(myid,gasOLD(i04),gas(i04),GasTermBal(i05), &
            k,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dtstep(ne))   !! change dt,chenhs
     enddo
   enddo
endif
!if(myid.eq.0) print *, "TEST0804-END_WETDEP"
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
!!************** SOAP SECSONDARY ORGANIC EROSOL  MODEL ***************
!***  TO PARTITION THE SEMIVOLATILE ORGANIC AEROSOL COMPONENTS *****
!***    6 PRECIES INCLUDING ANTHROPOGENIC TOL AND XYL          *****
!***            AND BIOGENIC ISOP AND TERP                     *****

    do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1

        i03      =  ip3mem(k,ne)

  !!!!!!!!!!!!!!!
    ! for SourceMark
    do ig=1, igas
      i04 = ip4mem(k,ig,ne)
     if(ifsm(ne)==1)then
        letdoit=0
      do idm=1,idmSet
        if(igMark(idm)==ig)letdoit=idm
      enddo
       if(letdoit>0)then
        contem0(letdoit)=gas(i04+ixy)
       endif
     endif! ifsm
    enddo  !ig
  !!!!!!!!!!!!!!!    

        TEMP     =  t(i03+ixy )    ! THE TEMPERATURE IN K
        PRESS0   =  Plev(i03+ixy ) ! THE PRESSURE IN HPA
       
        i04_1  = ip4mem(k,69, ne)    ! SV1
        i04_2  = ip4mem(k,70, ne)    ! SV2
        i04_3  = ip4mem(k,71, ne)    ! SV3
        i04_4  = ip4mem(k,72, ne)    ! SV4
        i04_5  = ip4mem(k,73, ne)    ! SV5
        i04_6  = ip4mem(k,74, ne)    ! SV6
        i04_7  = ip4mem(k,96, ne)    ! SOA1
        i04_8  = ip4mem(k,97, ne)    ! SOA2
        i04_9  = ip4mem(k,98, ne)    ! SOA3
        i04_10 = ip4mem(k,99, ne)    ! SOA4
        i04_11 = ip4mem(k,100,ne)    ! SOA5
        i04_12 = ip4mem(k,101,ne)    ! SOA6
        i04_13 = ip4mem(k,78,ne)     ! POA
        
        POA    = gas(i04_13+ixy)     ! POA in ug/m3
        SVG(1) = gas(i04_1 +ixy)     ! SV1
        SVG(2) = gas(i04_2 +ixy)     ! SV2
        SVG(3) = gas(i04_3 +ixy)     ! SV3
        SVG(4) = gas(i04_4 +ixy)     ! SV4
        SVG(5) = gas(i04_5 +ixy)     ! SV5
        SVG(6) = gas(i04_6 +ixy)     ! SV6
        SOA(1) = gas(i04_7 +ixy)     ! SOA1
        SOA(2) = gas(i04_8 +ixy)     ! SOA2
        SOA(3) = gas(i04_9 +ixy)     ! SOA3
        SOA(4) = gas(i04_10+ixy)     ! SOA4
        SOA(5) = gas(i04_11+ixy)     ! SOA5
        SOA(6) = gas(i04_12+ixy)     ! SOA6
        
      CALL SOAP(NSOA, SOA, SVG, TEMP, PRESS0, POA, I, J, K)
      
       gas(i04_13+ixy) = POA
       gas(i04_1 +ixy) = SVG(1)
       gas(i04_2 +ixy) = SVG(2)
       gas(i04_3 +ixy) = SVG(3)
       gas(i04_4 +ixy) = SVG(4)
       gas(i04_5 +ixy) = SVG(5)
       gas(i04_6 +ixy) = SVG(6)
       gas(i04_7 +ixy) = SOA(1)
       gas(i04_8 +ixy) = SOA(2)
       gas(i04_9 +ixy) = SOA(3)
       gas(i04_10+ixy) = SOA(4)
       gas(i04_11+ixy) = SOA(5)
       gas(i04_12+ixy) = SOA(6)

!!!!!!!!!!!!!!!
  ! for Source Mark
 DO ig=1,igas
  if(ifsm(ne)==1)then
      letdoit=0
    do idm=1,idmSet
       if(igMark(idm)==ig)letdoit=idm
    enddo
    if(letdoit>0)then
       i02 = ip2mem(ne)
       MapS=int(MapSource(i02+ixy))
       IHgtLev= iHgtLMax-1  ! for chemistry reaction part Zifa/2007/03/03
                   ! this is need to be defined outside for future
       i04 = ip4mem(k,ig,ne)
       OrgCon     = contem0(letdoit)
       DeltSpeMark= gas(i04+ixy)-OrgCon
       if(DeltSpeMark > 0. )then
            do ism=1,ismMax
                i04sm=ipSMmem(k,ism,letdoit,ne)
                TmpSM(ism)=SourceMark(i04sm+ixy)
            enddo
            call GetSMarkChange(DeltSpeMark,OrgCon,TmpSM,MapS, &
                 ISrcDefined,IHgtLev,ismMax)
            do ism=1,ismMax
               i04sm=ipSMmem(k,ism,letdoit,ne)
               SourceMark(i04sm+ixy)=TmpSM(ism)
            enddo
        endif
     endif ! letdoit
    endif  ! ifsm
  ENDDO !IG
!!!!!!
       
       enddo !k
      enddo  !i
    enddo    !j
!!!!!!!!!

!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
if(ifprocess.eq.1) then
 ! iprocess
   do k=1,nzz
     do ig=1,iPrintTermGas    ! for SOA process
       igg=IGGPOS(ig)
       i04=ip4mem(k,igg,ne)
       i05=ipGasTermBal(k,23,ig,ne)  ! 23 for organic aerosol partioning
    call termbal(myid,gasOLD(i04),gas(i04),GasTermBal(i05), &
            k,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dtstep(ne))  !! change dt,chenhs
     enddo
   enddo
endif
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
!if(myid.eq.0) print *, "TEST0805-END_SOAP"
!!************** Exitinction and AOD  MODEL **************************
!!***         TO USE THE RRECONSTRUCTED METHOD          **************

   do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1

       ALLOCATE(EXT0(NZZ),EXTS(NZZ),DUSTEXT0(NZZ))    
       
       i02          =  ip2mem(ne)
       AOD(i02+ixy) =     0.0
       DUSTAOD(I02+IXY) = 0.0
       AODS         =     0.0
       PBLAOD(i02+ixy) =  0.0 
             
      do k=1,nzz-1

       i03      =  ip3mem(k,ne)
                  
       RH  =  rh1(i03+ixy)         ! the RH IN %

        i04_1  = ip4mem(k, 81, ne) ! NH4+(AE)
        i04_2  = ip4mem(k, 83, ne) ! SO42-(AE)
        i04_3  = ip4mem(k, 84, ne) ! HSO4-(AE)
        i04_4  = ip4mem(k, 85, ne) ! NO3-(AE)
        i04_5  = ip4mem(k, 87, ne) ! NA2SO4(S)
        i04_6  = ip4mem(k, 88, ne) ! NANO3(S)
        i04_7  = ip4mem(k, 89, ne) ! (NH4)2SO4(S)
        i04_8  = ip4mem(k, 90, ne) ! NH4NO3(S)
        i04_9  = ip4mem(k, 91, ne) ! NH4CL(S)
        i04_10 = ip4mem(k, 92, ne) ! H2SO4(AQ)
        i04_11 = ip4mem(k, 93, ne) ! NH4HSO4(S)
        i04_12 = ip4mem(k, 94, ne) ! NAHSO4(S)
        i04_13 = ip4mem(k, 95, ne) ! (NH4)3H(SO4)2(S)
        i04_14 = ip4mem(k, 96, ne) ! SOA1
        i04_15 = ip4mem(k, 97, ne) ! SOA2
        i04_16 = ip4mem(k, 98, ne) ! SOA3
        i04_17 = ip4mem(k, 99, ne) ! SOA4
        i04_18 = ip4mem(k,100, ne) ! SOA5
        i04_19 = ip4mem(k,101, ne) ! SOA6
        i04_20 = ip4mem(k, 78, ne) ! POA
        i04_21 = ip4mem(k, 77, ne) ! BC
        i04_22 = ip4mem(k, 76, ne) ! PM10
        i04_23 = ip4mem(k, 75, ne) ! PM25
        i04_24 = ip4mem(k, 80, ne) ! NA
                                   
        NH4 =  gas(i04_1+ixy) * 17.0/GC_MOLWT(81)&
             + gas(i04_7+ixy) * 17.0/GC_MOLWT(89)*2.0 &
             + gas(i04_8+ixy) * 17.0/GC_MOLWT(90)&
             + gas(i04_9+ixy) * 17.0/GC_MOLWT(91)&
             + gas(i04_11+ixy)* 17.0/GC_MOLWT(93)&
             + gas(i04_13+ixy)* 17.0/GC_MOLWT(95)*3.0

        SO4 =  gas(i04_2+ixy) * 98.0/GC_MOLWT(83)&
             + gas(i04_3+ixy) * 98.0/GC_MOLWT(84)&
             + gas(i04_5+ixy) * 98.0/GC_MOLWT(87)&
             + gas(i04_7+ixy) * 98.0/GC_MOLWT(89)&
             + gas(i04_10+ixy)* 98.0/GC_MOLWT(92)&
             + gas(i04_11+ixy)* 98.0/GC_MOLWT(93)&
             + gas(i04_12+ixy)* 98.0/GC_MOLWT(94)&     
             + gas(i04_13+ixy)* 98.0/GC_MOLWT(95)*2.0
             
        NO3 =  gas(i04_4+ixy) * 63.0/GC_MOLWT(85)&
             + gas(i04_6+ixy) * 63.0/GC_MOLWT(88)&
             + gas(i04_8+ixy) * 63.0/GC_MOLWT(90)
             
        NA  =  gas(i04_23+ixy) * 23.0/GC_MOLWT(80)&
             + gas(i04_5+ixy)  * 23.0/GC_MOLWT(87)*2.0&
             + gas(i04_6+ixy)  * 23.0/GC_MOLWT(88)&
             + gas(i04_12+ixy) * 23.0/GC_MOLWT(94)
             
        BC  =  gas(i04_21+ixy)   

        OC  =  gas(i04_14+ixy) + gas(i04_15+ixy) &
             + gas(i04_16+ixy) + gas(i04_17+ixy) &
             + gas(i04_18+ixy) + gas(i04_19+ixy) &
             + gas(i04_20+ixy)
         
        PM10 = gas(i04_22+ixy)        

        PM25 = gas(i04_24+ixy)     
  
          I05_1=IP5MEM(K,1,2,NE)
          I05_2=IP5MEM(K,2,2,NE)
          I05_3=IP5MEM(K,3,2,NE)
          I05_4=IP5MEM(K,4,2,NE)

          DUST01= AER(I05_1+IXY)
          DUST02= AER(I05_2+IXY)
          DUST03= AER(I05_3+IXY)
          DUST04= AER(I05_4+IXY)

          I05_1=IP5MEM(K,1,1,NE)
          I05_2=IP5MEM(K,2,1,NE)
          I05_3=IP5MEM(K,3,1,NE)
          I05_4=IP5MEM(K,4,1,NE)

          SEA01 = AER (I05_1+IXY)
          SEA02 = AER (I05_2+IXY)
          SEA03 = AER (I05_3+IXY)
          SEA04 = AER (I05_4+IXY)


        CALL GETEXT ( MYID, RH, SEA01, SEA02, SEA03, SEA04,DUST01,DUST02,DUST03,DUST04,&
                     NH4, SO4, NO3, NA, BC, OC, PM10, PM25,&
                     DUSTEXT0(K), EXT0(K),EXTS(K),2, I, J, K) 


        EXT(i03+ixy) = EXT0(K)
        DUSTEXT(I03+IXY) = DUSTEXT0(K)
        
        CALL GETVIS (MYID, RNW(i03+ixy),CLW(i03+ixy) , EXT(i03+ixy), VISIB(i03+ixy), I, J, K)
        
        AOD(i02+ixy) = AOD(i02+ixy) + EXT0(K)*dz(i03+ixy)/1.E03
        DUSTAOD(i02+ixy) = DUSTAOD(i02+ixy) + DUSTEXT0(K)*dz(i03+ixy)/1.E03
      

        AODS = AODS + EXTS(K)*dz(i03+ixy)/1.E03
        
        SSA(i03+ixy) = EXTS(K)/MAX(EXT0(K), 1.E-20)
        
        IF(K.LE.NPBL(i02+ixy)) &
        PBLAOD(i02+ixy) = PBLAOD(i02+ixy) + EXT0(K)*dz(i03+ixy)/1.E03
     enddo ! k
     
    DEALLOCATE(EXT0,EXTS,DUSTEXT0) 
   enddo !i
  enddo  !j  
!
!
!*********     TO GET THE Dobson (DU) FOR SO2, NO2, O3  *************
!
    DO j = sy(ne),ey(ne)
     DO  i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
      
      ALLOCATE(DSO2(NZZ),DNO2(NZZ),DO3(NZZ))
       i02          =  ip2mem(ne)
       DUSO2 (i02+ixy) = 0.0
       DUNO2 (i02+ixy) = 0.0
       DUO3  (i02+ixy) = 0.0      
     
       kkbb = ktop(i02+ixy) 
      DO k=1,kkbb-1
        i03      =  ip3mem(k,ne)
        temp     =  t(i03+ixy)          ! the temp
        PRESS0   =  Plev(i03+ixy)       ! pressuer in hpa

        i04_1  = ip4mem(k, 6,  ne) ! NO2
        i04_2  = ip4mem(k, 11, ne) ! O3
        i04_3  = ip4mem(k, 18, ne) ! SO2

        CALL GETDU(MYID, gas(i04_1+ixy),DNO2(K), PRESS0, temp )
        CALL GETDU(MYID, gas(i04_2+ixy),DO3 (K), PRESS0, temp )
        CALL GETDU(MYID, gas(i04_3+ixy),DSO2(K), PRESS0, temp )

        DUSO2 (i02+ixy) = DUSO2 (i02+ixy) + DSO2(K)*dz(i03+ixy)/2.69E+20
        DUO3  (i02+ixy) = DUO3  (i02+ixy) + DO3 (K)*dz(i03+ixy)/2.69E+20
        DUNO2 (i02+ixy) = DUNO2 (i02+ixy) + DNO2(K)*dz(i03+ixy)/2.69E+20

      ENDDO  !K

      DEALLOCATE(DSO2,DO3,DNO2)
    ENDDO  !I
    ENDDO  !J 

! if(myid.eq.0) print *, "TEST09--Chem End", itt
!=====================================================================
!====    to change boundary conditions                           =====
!=====================================================================

if(numprocs .gt. 1 .and. mod(it1,iprecise)==0 )then

call mpi_barrier( local_com, ierr )

do k=1,nzz
  i03=ip3mem(k,ne)
  call exchng2( myid, u(i03), sx(ne), ex(ne), sy(ne), ey(ne),  &
                comm2d(ne), stride(ne),  nbrleft(ne), nbrright(ne), &
                nbrtop(ne), nbrbottom(ne) )
  call exchng2( myid, v(i03), sx(ne), ex(ne), sy(ne), ey(ne),  &
                comm2d(ne), stride(ne),  nbrleft(ne), nbrright(ne), &
                nbrtop(ne), nbrbottom(ne) )
  do ia=1,iaer
  do is=1,isize
  i05=ip5mem(k,is,ia,ne)
  call exchng2( myid, aer(i05), sx(ne), ex(ne), sy(ne), ey(ne),  &
                comm2d(ne), stride(ne),  nbrleft(ne), nbrright(ne), &
                nbrtop(ne), nbrbottom(ne) )
  enddo
  enddo

  do ig=1,igas
  i04=ip4mem(k,ig,ne)
  call exchng2( myid, gas(i04), sx(ne), ex(ne), sy(ne), ey(ne),  &
                comm2d(ne), stride(ne),  nbrleft(ne), nbrright(ne), &
                nbrtop(ne), nbrbottom(ne) )
  enddo

  ! For Source Mark
  if(ifsm(ne)==1)then 
  do idm=1,idmSet
  do ism=1,ismMax 
  i04=ipSMmem(k,ism,idm,ne)
  call exchng2( myid, SourceMark(i04), sx(ne), ex(ne), sy(ne), ey(ne),  &
                comm2d(ne), stride(ne),  nbrleft(ne), nbrright(ne), &
                nbrtop(ne), nbrbottom(ne) )
  enddo 
  enddo 
  endif 

! for sea salt and dust
  if(ifseacom.eq.1) then
  do iduc = 1, nseacom
  do is = 1, isize
    i05c = ip5memcs (k, is, iduc, ne) 
    call exchng2 ( myid, seacomp(i05c), sx(ne), ex(ne), sy(ne), ey(ne),  &
                comm2d(ne), stride(ne),  nbrleft(ne), nbrright(ne), &
                nbrtop(ne), nbrbottom(ne))
  enddo
  enddo
  endif

  if(ifdustcom.eq.1) then
  do iduc = 1, ndustcom
  do is = 1, isize
    i05c = ip5memc (k, is, iduc, ne)
    call exchng2 ( myid, dustcomp(i05c), sx(ne), ex(ne), sy(ne), ey(ne),  &
                comm2d(ne), stride(ne),  nbrleft(ne), nbrright(ne), &
                nbrtop(ne), nbrbottom(ne))
  enddo
  enddo
  endif
  
enddo
!!++++++++++++++++++++++ chenhs,polar transport +++++++++++++++++++++++++!!
if(ifglobal.eq.-1) then
  IF(NE==1) THEN
      DO IPS=1,IPOLARNUM
       IF(IPOLARMRK(2,IPS)==MYID) THEN  !!! need to send
        DO K=1,NZZ
           ITSP=0
          ! For gas speceis
          DO IG=1,igas
           ITSP=ITSP+1
           I04 = IP4MEM(K,IG,NE)
          CALL GETVALUE(MYID,gas(I04),SX(NE),EX(NE),SY(NE),EY(NE),&
                       IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),ATESTS(K,ITSP))
          ENDDO
          !For aerosols
          DO IA=1,IAER
          DO IS=1,ISIZE
           I05= IP5MEM(K,IS,IA,NE)
           ITSP=ITSP+1
          CALL GETVALUE(MYID,AER(I05),SX(NE),EX(NE),SY(NE),EY(NE),&
              IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),ATESTS(K,ITSP))
          ENDDO
          ENDDO
           ! FOR U,V 
           I03= IP3MEM(K,NE)
           ITSP=ITSP+1
           CALL GETVALUE(MYID,U(I03),SX(NE),EX(NE),SY(NE),EY(NE),&
              IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),ATESTS(K,ITSP))
           ITSP=ITSP+1
           CALL GETVALUE(MYID,V(I03),SX(NE),EX(NE),SY(NE),EY(NE),&
              IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),ATESTS(K,ITSP))
        !!!!!!!!!!!!!!!!!!!!
        ! for Source Mark
        if(ifsmt>0)then    ! checking 
         do idm=1,idmSet
         do ism=1,ismMax
            i0=ipSMmem(k,ism,idm,ne)
            itsp=itsp+1
         call getvalue(myid,SourceMark(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),atestS(k,itsp))
         enddo
         enddo
        endif
        !!!!!!!!!!!!!!!!!!!!
        ! for dust and sea salt composition
        if(ifseacom.eq.1) then
        do iduc = 1, nseacom
        do is =1 , isize
          i05c=ip5memcs(k,is,iduc,ne)
          itsp = itsp + 1
     call getvalue(myid,seacomp(i05c),sx(ne),ex(ne),sy(ne),ey(ne), &
             IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),atestS(k,itsp))
        enddo
        enddo
        endif
       
        if(ifdustcom.eq.1) then
        do iduc = 1, ndustcom
        do is =1 , isize
          i05c=ip5memc(k,is,iduc,ne)
          itsp = itsp + 1
     call getvalue(myid, dustcomp(i05c) ,sx(ne),ex(ne),sy(ne),ey(ne), &
             IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),atestS(k,itsp))
        enddo
        enddo
        endif
        ENDDO !! k,1,nzz
        DO K=1,NZZ
           IPSMARK = IPS + 720 + (k-1)*IPOLARNUM
           do ibeibei=1,ISRNUM
           atestS0(ibeibei)=atestS(k,ibeibei)
          enddo
        call MPI_Send(atestS0,ISRNUM,MPI_REAL,IPOLARMRK(1,IPS),IPSMARK,comm2d(ne),ierr)
        ENDDO
      ENDIF

      IF(IPOLARMRK(1,IPS)==MYID) THEN  !!! need to receive
         DO K=1,NZZ
          IPSMARK = IPS + 720 + (k-1)*IPOLARNUM
           call MPI_Recv(atestR0,ISRNUM,MPI_REAL,IPOLARMRK(2,IPS),IPSMARK,  &
                comm2d(ne),status,ierr)
          do ibeibei=1,ISRNUM
           atestR(k,ibeibei)=atestR0(ibeibei)
          enddo
         ENDDO
         DO K=1,NZZ
           IF(IPOLARMRK(5,IPS)==1.OR.IPOLARMRK(5,IPS)==NY(NE)) THEN
             IF(IPOLARMRK(5,IPS)==1) THEN
               IPP=IPOLARMRK(5,IPS)-1
             ELSE
               IPP=IPOLARMRK(5,IPS)+1
             ENDIF

             ips1=0
             ! For gas speceis
             DO IG=1,igas
              ips1 = ips1 + 1
              I04 = IP4MEM(K,IG,NE)
              CALL PUTVALUE(MYID,gas(i04),SX(NE),EX(NE),SY(NE),EY(NE),&
              IPOLARMRK(3,IPS),IPP,ATESTR(K,ips1))
             ENDDO
              ! for aerosols
              DO IA=1,IAER
              DO IS=1,ISIZE
              ips1 = ips1 + 1
              I05= IP5MEM(K,IS,IA,NE)
              CALL PUTVALUE(MYID,AER(I05),SX(NE),EX(NE),SY(NE),EY(NE),&
              IPOLARMRK(3,IPS),IPP,ATESTR(K,IPS1))
              ENDDO
              ENDDO
              ! for uv
              ips1 = ips1 + 1
              I03 = IP3MEM(K,NE)
             CALL PUTVALUE(MYID,U(I03),SX(NE),EX(NE),SY(NE),EY(NE),&
              IPOLARMRK(3,IPS),IPP,ATESTR(K,IPS1))
              ips1 = ips1 + 1
              I03 = IP3MEM(K,NE)
             CALL PUTVALUE(MYID,V(I03),SX(NE),EX(NE),SY(NE),EY(NE),&
              IPOLARMRK(3,IPS),IPP,ATESTR(K,IPS1))
           !!!!!!!!!!!!!!!!!
          ! for Source Mark
          if(ifsmt>0 )then   !checking
            do idm=1,idmSet
            do ism=1,ismMax
             !ips1 = igas + iaer*isize +2+(idm-1)*ismMax + ism
             ips1 = ips1 + 1
             i0=ipSMmem(k,ism,idm,ne)
            call putvalue(myid,SourceMark(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              IPOLARMRK(3,IPS),IPP,atestR(k,ips1))
            enddo
            enddo
           endif
           !!!!!!!!!!!!!!!!!!!!
           ! for sea salt and dust compositions
          if(ifseacom.eq.1) then
          do iduc = 1, nseacom
          do is = 1, isize
            ips1 = ips1 + 1
            i05c  = ip5memcs(k,is,iduc,ne)
     call putvalue(myid,seacomp(i05c),sx(ne),ex(ne),sy(ne),ey(ne), &
              IPOLARMRK(3,IPS),IPP,atestR(k,ips1))
          enddo
          enddo
          endif
          
          if(ifdustcom.eq.1) then 
          do iduc = 1, ndustcom
          do is = 1, isize
            ips1 = ips1 + 1
            i05c  = ip5memc(k,is,iduc,ne)
     call putvalue(myid,dustcomp(i05c),sx(ne),ex(ne),sy(ne),ey(ne), &
              IPOLARMRK(3,IPS),IPP,atestR(k,ips1))
          enddo
          enddo
          endif
            
           ENDIF
         ENDDO

      ENDIF
     ENDDO
   ENDIF
 
endif
!!++++++++++++++++++++++ chenhs,polar transport +++++++++++++++++++++++++!!
endif ! change boundary conditions
!=====================================================================
  endif ! added by chenhs, to adjust integration time step 
endif   !
!----------------------------------------------------------------------
enddo  ! nest end calculation
!----------------------------------------------------------------------

enddo   ! end 1 hour itt<=12

!--------------------------------------------------------------------
! to write data
!--------------------------------------------------------------------	
	do ne=1,nest
	
	if(it1.ge.ntbeg(ne) .and. it1.le.ntend(ne))then
	if(mod(ihour2,dtout(ne))==0) then  !! added by chenhs,to adjust output frequency
	call output3(iyear2,imonth2,iday2,ihour2,it1,ne) ! aerosol
	call output_gas(iyear2,imonth2,iday2,ihour2,it1,ne) ! gas
	call output_dry(iyear2,imonth2,iday2,ihour2,ne) ! dry
	call output_wet(iyear2,imonth2,iday2,ihour2,ne) ! wet 
	call output_dust(iyear2,imonth2,iday2,ihour2,ne) ! dust
	call output_sea(iyear2,imonth2,iday2,ihour2,ne) ! sea
	
	!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
	if(ifprocess.eq.1) then
	call output_term(iyear2,imonth2,iday2,ihour2,ne)
	endif
	!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
	
	! for Source Mark
	if(ifsm(ne)==1)then 
        call output_mark(iyear2,imonth2,iday2,ihour2,ne,.false.)
	IF(ihour2==0)THEN  ! to write down reinit source file
        call output_mark(iyear2,imonth2,iday2,ihour2,ne,.true.)	
	endif 
	endif

       endif
       endif
       enddo
1188 continue

        it1=it1+1 ! advance time
  end do   ! end of the while integration

end subroutine geatm_run_mct
!================================================================================
subroutine geatm_final_mct()
   use geatm_vartype
   if(allocated(procs)) deallocate (procs)
   if(allocated(cam_local_points))    deallocate(cam_local_points)

   call final_camgrid(camgrid)
   call final_wrfgrid(wrfgrid)
   call final_geatm_var
   do ne=1,nest
    call MPI_TYPE_FREE( stride(ne), ierr )
    call MPI_COMM_FREE( comm2d(ne), ierr )
   enddo

end subroutine geatm_final_mct

!===============================================================================
	
  subroutine geatm_SetgsMap_cam( camgrid, local_com, GEAID, GSMap_gege )
!-------------------------------------------------------------------
!
! Arguments
!
      use geatm_vartype, only:csx,cex,csy,cey
      TYPE(camgrid_c) :: camgrid
      integer        , intent(in)  :: local_com
      integer        , intent(in)  :: GEAID
      type(mct_gsMap), intent(out) :: GSMap_gege
!
! Local variables
!
      integer, allocatable :: gindex(:)
      integer :: i, j, k, n,lsize,gsize
      integer :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,ips,ipe,jps,jpe,kps,kpe
      integer :: ier            ! error status
!-------------------------------------------------------------------
! Build the atmosphere grid numbering for MCT
! NOTE:  Numbering scheme is: West to East, bottom to level, and South to North
! starting at south pole.  Should be the same as what's used in SCRIP
! Determine global seg map
	
      ! prepare the decomposition	
      ips=csx(1)
      ipe=cex(1)
      jps=csy(1)
      jpe=cey(1)
      ids=1
      ide=camgrid%num_camgrid_lon
      jds=1
      jde=camgrid%num_camgrid_lat
      
      lsize=0
      do j=jps, jpe 
       do i=ips, ipe
             lsize = lsize+1  !local index
       end do
      end do
   
      gsize=(ide-ids+1)*(jde-jds+1)
      allocate(gindex(lsize))
      n=0
      do j=jps, jpe  
       do i=ips, ipe
          n=n+1
          gindex(n) =(j-1)*(ide-ids+1)+i  ! global index
       end do
      end do

      call mct_gsMap_init( gsMap_gege, gindex, local_com, GEAID, lsize, gsize)

      deallocate(gindex)

  end subroutine geatm_SetgsMap_cam
  
!===============================================================================

  subroutine geatm_domain_cam( camgrid, lsize, gsMap_cc, dom_c )

      use geatm_vartype, only:myid
!-------------------------------------------------------------------
! Arguments
!
      TYPE(camgrid_c) :: camgrid
      integer     , intent(in)   :: lsize
      TYpe(mct_gsMap), intent(in)   :: gsMap_cc
      type(mct_ggrid), intent(inout):: dom_c    
!
! Local Variables
!
      integer  :: i,j,k,mm,nn           ! indices	
      integer :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,ips,ipe,jps,jpe,kps,kpe
      real(r8), pointer  :: data(:)     ! temporary
      integer , pointer  :: idata(:)    ! temporary
                
! Initialize mct atm domain
         call mct_gGrid_init( GGrid=dom_c, CoordChars=trim(seq_flds_dom_coord), OtherChars=trim(seq_flds_dom_other), lsize=lsize )

! Allocate memory
         allocate(data(lsize))
 
! Initialize attribute vector with special value
         call mct_gsMap_orderedPoints(gsMap_cc, myid, idata)
         call mct_gGrid_importIAttr(dom_c,'GlobGridNum',idata,lsize)

! Determine domain (numbering scheme is: West to East and South to North to South pole)
! Initialize attribute vector with special value

        data(:) = -9999.0_R8 
        call mct_gGrid_importRAttr(dom_c,"lat"  ,data,lsize) 
        call mct_gGrid_importRAttr(dom_c,"lon"  ,data,lsize) 
        call mct_gGrid_importRAttr(dom_c,"area" ,data,lsize) 
        call mct_gGrid_importRAttr(dom_c,"aream",data,lsize) 
        data(:) = 0.0_R8     
        call mct_gGrid_importRAttr(dom_c,"mask" ,data,lsize) 
        data(:) = 1.0_R8
        call mct_gGrid_importRAttr(dom_c,"frac" ,data,lsize)
    
        if(associated(idata)) deallocate(idata)    
	deallocate(data)

	end subroutine geatm_domain_cam
	
!===============================================================================
  subroutine initial_camgrid( camgrid)
     TYPE(camgrid_c) :: camgrid
     
     ! Local variables	
     integer  :: i,j,k  ! indices    
     integer :: mids,mide,mjds,mjde,mkds,mkde,&
                 mims,mime,mjms,mjme,mkms,mkme,&
                 mips,mipe,mjps,mjpe,mkps,mkpe        
   
      mids=1
      mide=camgrid%num_camgrid_lon
      mjds=1
      mjde=camgrid%num_camgrid_lat
      mkps=1
      mkpe=camgrid%num_camgrid_levels
      if(.not.allocated(camgrid%z3d)) allocate(camgrid%z3d(mids:mide,mkps:mkpe,mjds:mjde))
      if(.not.allocated(camgrid%u3d)) allocate(camgrid%u3d(mids:mide,mkps:mkpe,mjds:mjde))
      if(.not.allocated(camgrid%v3d)) allocate(camgrid%v3d(mids:mide,mkps:mkpe,mjds:mjde))
      if(.not.allocated(camgrid%t3d)) allocate(camgrid%t3d(mids:mide,mkps:mkpe,mjds:mjde))
      if(.not.allocated(camgrid%p3d)) allocate(camgrid%p3d(mids:mide,mkps:mkpe,mjds:mjde))
      if(.not.allocated(camgrid%qv3d)) allocate(camgrid%qv3d(mids:mide,mkps:mkpe,mjds:mjde))
      if(.not.allocated(camgrid%qc3d)) allocate(camgrid%qc3d(mids:mide,mkps:mkpe,mjds:mjde))
      if(.not.allocated(camgrid%qi3d)) allocate(camgrid%qi3d(mids:mide,mkps:mkpe,mjds:mjde))      
      if(.not.allocated(camgrid%rh3d)) allocate(camgrid%rh3d(mids:mide,mkps:mkpe,mjds:mjde))
      if(.not.allocated(camgrid%taucldv3d)) allocate(camgrid%taucldv3d(mids:mide,mkps:mkpe,mjds:mjde))
      if(.not.allocated(camgrid%taucldi3d)) allocate(camgrid%taucldi3d(mids:mide,mkps:mkpe,mjds:mjde))      
      if(.not.allocated(camgrid%xlat)) allocate(camgrid%xlat(mids:mide,mjds:mjde))
      if(.not.allocated(camgrid%xlon)) allocate(camgrid%xlon(mids:mide,mjds:mjde))
      if(.not.allocated(camgrid%ps)) allocate(camgrid%ps(mids:mide,mjds:mjde))
      if(.not.allocated(camgrid%ht)) allocate(camgrid%ht(mids:mide,mjds:mjde))
      if(.not.allocated(camgrid%pblh)) allocate(camgrid%pblh(mids:mide,mjds:mjde))
      if(.not.allocated(camgrid%u10)) allocate(camgrid%u10(mids:mide,mjds:mjde))
      if(.not.allocated(camgrid%v10)) allocate(camgrid%v10(mids:mide,mjds:mjde))
      if(.not.allocated(camgrid%t2)) allocate(camgrid%t2(mids:mide,mjds:mjde))
      if(.not.allocated(camgrid%q2)) allocate(camgrid%q2(mids:mide,mjds:mjde))
      if(.not.allocated(camgrid%rh2)) allocate(camgrid%rh2(mids:mide,mjds:mjde))
      if(.not.allocated(camgrid%xland)) allocate(camgrid%xland(mids:mide,mjds:mjde))
      if(.not.allocated(camgrid%sst)) allocate(camgrid%sst(mids:mide,mjds:mjde))
      if(.not.allocated(camgrid%xice)) allocate(camgrid%xice(mids:mide,mjds:mjde))
      if(.not.allocated(camgrid%tsk)) allocate(camgrid%tsk(mids:mide,mjds:mjde))
      if(.not.allocated(camgrid%snowh)) allocate(camgrid%snowh(mids:mide,mjds:mjde))      
      if(.not.allocated(camgrid%ust)) allocate(camgrid%ust(mids:mide,mjds:mjde))
      if(.not.allocated(camgrid%rmol)) allocate(camgrid%rmol(mids:mide,mjds:mjde))
      if(.not.allocated(camgrid%raincv)) allocate(camgrid%raincv(mids:mide,mjds:mjde))
      if(.not.allocated(camgrid%rainncv)) allocate(camgrid%rainncv(mids:mide,mjds:mjde))
      if(.not.allocated(camgrid%swdown)) allocate(camgrid%swdown(mids:mide,mjds:mjde))
      if(.not.allocated(camgrid%clflo)) allocate(camgrid%clflo(mids:mide,mjds:mjde))
      if(.not.allocated(camgrid%clfmi)) allocate(camgrid%clfmi(mids:mide,mjds:mjde))
      if(.not.allocated(camgrid%clfhi)) allocate(camgrid%clfhi(mids:mide,mjds:mjde))
      
      if(.not.allocated(camgrid%soildepth)) allocate(camgrid%soildepth(mids:mide,1:camgrid%num_camgrid_soil_levels,mjds:mjde))
      if(.not.allocated(camgrid%soilthick)) allocate(camgrid%soilthick(mids:mide,1:camgrid%num_camgrid_soil_levels,mjds:mjde))
      if(.not.allocated(camgrid%soilt)) allocate(camgrid%soilt(mids:mide,1:camgrid%num_camgrid_soil_levels,mjds:mjde))
      if(.not.allocated(camgrid%soilm)) allocate(camgrid%soilm(mids:mide,1:camgrid%num_camgrid_soil_levels,mjds:mjde))

  end subroutine initial_camgrid
!-----------------------------------------------------------------------
  subroutine final_camgrid(camgrid)
    type(camgrid_c) :: camgrid
    if(allocated(camgrid%z3d))     deallocate(camgrid%z3d)
    if(allocated(camgrid%u3d))     deallocate(camgrid%u3d)
    if(allocated(camgrid%v3d))     deallocate(camgrid%v3d)
    if(allocated(camgrid%qv3d))   deallocate(camgrid%qv3d)
    if(allocated(camgrid%qc3d))   deallocate(camgrid%qc3d)
    if(allocated(camgrid%qi3d))   deallocate(camgrid%qi3d)       
    if(allocated(camgrid%t3d))   deallocate(camgrid%t3d)
    if(allocated(camgrid%rh3d))   deallocate(camgrid%rh3d)
    if(allocated(camgrid%p3d))   deallocate(camgrid%p3d)
    if(allocated(camgrid%taucldv3d))   deallocate(camgrid%taucldv3d)
    if(allocated(camgrid%taucldi3d))   deallocate(camgrid%taucldi3d)       
    if(allocated(camgrid%xlat))   deallocate(camgrid%xlat)
    if(allocated(camgrid%xlon))   deallocate(camgrid%xlon)
    if(allocated(camgrid%ps))   deallocate(camgrid%ps)
    if(allocated(camgrid%q2))   deallocate(camgrid%q2)
    if(allocated(camgrid%t2))   deallocate(camgrid%t2)
    if(allocated(camgrid%rh2))   deallocate(camgrid%rh2)
    if(allocated(camgrid%v10))   deallocate(camgrid%v10)
    if(allocated(camgrid%u10))   deallocate(camgrid%u10)
    if(allocated(camgrid%ht))   deallocate(camgrid%ht)
    if(allocated(camgrid%rmol))   deallocate(camgrid%rmol)
    if(allocated(camgrid%pblh))   deallocate(camgrid%pblh)
    if(allocated(camgrid%raincv))   deallocate(camgrid%raincv)
    if(allocated(camgrid%rainncv))   deallocate(camgrid%rainncv)
    if(allocated(camgrid%snowh))   deallocate(camgrid%snowh)       
    if(allocated(camgrid%sst))   deallocate(camgrid%sst)
    if(allocated(camgrid%xice))   deallocate(camgrid%xice)
    if(allocated(camgrid%tsk))   deallocate(camgrid%tsk)       
    if(allocated(camgrid%xland))   deallocate(camgrid%xland)
    if(allocated(camgrid%swdown))   deallocate(camgrid%swdown)
    if(allocated(camgrid%clflo))   deallocate(camgrid%clflo)
    if(allocated(camgrid%clfmi))   deallocate(camgrid%clfmi)
    if(allocated(camgrid%clfhi))   deallocate(camgrid%clfhi)

    if(allocated(camgrid%soildepth)) deallocate(camgrid%soildepth)
    if(allocated(camgrid%soilthick)) deallocate(camgrid%soilthick)
    if(allocated(camgrid%soilt)) deallocate(camgrid%soilt)
    if(allocated(camgrid%soilm)) deallocate(camgrid%soilm)

  end subroutine final_camgrid
!-----------------------------------------------------------------------
  subroutine initial_wrfgrid( wrfgrid)
     TYPE(wrfgrid_c) :: wrfgrid
     
     ! Local variables	
     integer  :: i,j,k  ! indices    
     integer :: mids,mide,mjds,mjde,mkds,mkde,&
                 mims,mime,mjms,mjme,mkms,mkme,&
                 mips,mipe,mjps,mjpe,mkps,mkpe        

      wrfgrid%num_wrfgrid_lon = 360
      wrfgrid%num_wrfgrid_lat = 180
      wrfgrid%num_wrfgrid_levels = 20
      wrfgrid%num_wrfgrid_soil_levels =4
     
      mids=1
      mide=wrfgrid%num_wrfgrid_lon
      mjds=1
      mjde=wrfgrid%num_wrfgrid_lat
      mkps=1
      mkpe=wrfgrid%num_wrfgrid_levels    
      if(.not.allocated(wrfgrid%sigma)) allocate(wrfgrid%sigma(mkps:mkpe))
      wrfgrid%sigma(1) =   0.00250_8
      wrfgrid%sigma(2) =   0.00800_8
      wrfgrid%sigma(3) =   0.01500_8
      wrfgrid%sigma(4) =   0.02350_8
      wrfgrid%sigma(5) =   0.03400_8
      wrfgrid%sigma(6) =   0.04650_8
      wrfgrid%sigma(7) =   0.06100_8
      wrfgrid%sigma(8) =   0.07850_8
      wrfgrid%sigma(9) =   0.09950_8
      wrfgrid%sigma(10) =   0.12500_8
      wrfgrid%sigma(11) =   0.15550_8
      wrfgrid%sigma(12) =   0.19200_8
      wrfgrid%sigma(13) =   0.23550_8
      wrfgrid%sigma(14) =   0.28800_8
      wrfgrid%sigma(15) =   0.35150_8
      wrfgrid%sigma(16) =   0.42700_8
      wrfgrid%sigma(17) =   0.51750_8
      wrfgrid%sigma(18) =   0.62650_8
      wrfgrid%sigma(19) =   0.75750_8
      wrfgrid%sigma(20) =   0.91450_8
 
      if(.not.allocated(wrfgrid%z3d)) allocate(wrfgrid%z3d(mids:mide,mkps:mkpe,mjds:mjde))
      if(.not.allocated(wrfgrid%u3d)) allocate(wrfgrid%u3d(mids:mide,mkps:mkpe,mjds:mjde))
      if(.not.allocated(wrfgrid%v3d)) allocate(wrfgrid%v3d(mids:mide,mkps:mkpe,mjds:mjde))
      if(.not.allocated(wrfgrid%t3d)) allocate(wrfgrid%t3d(mids:mide,mkps:mkpe,mjds:mjde))
      if(.not.allocated(wrfgrid%p3d)) allocate(wrfgrid%p3d(mids:mide,mkps:mkpe,mjds:mjde))
      if(.not.allocated(wrfgrid%qv3d)) allocate(wrfgrid%qv3d(mids:mide,mkps:mkpe,mjds:mjde))
      if(.not.allocated(wrfgrid%qc3d)) allocate(wrfgrid%qc3d(mids:mide,mkps:mkpe,mjds:mjde))
      if(.not.allocated(wrfgrid%qi3d)) allocate(wrfgrid%qi3d(mids:mide,mkps:mkpe,mjds:mjde))      
      if(.not.allocated(wrfgrid%rh3d)) allocate(wrfgrid%rh3d(mids:mide,mkps:mkpe,mjds:mjde))
      if(.not.allocated(wrfgrid%taucldv3d)) allocate(wrfgrid%taucldv3d(mids:mide,mkps:mkpe,mjds:mjde))
      if(.not.allocated(wrfgrid%taucldi3d)) allocate(wrfgrid%taucldi3d(mids:mide,mkps:mkpe,mjds:mjde))      
      ! for geatm variables output
      if(.not.allocated(wrfgrid%ps)) allocate(wrfgrid%ps(mids:mide,mjds:mjde))
      if(.not.allocated(wrfgrid%ht)) allocate(wrfgrid%ht(mids:mide,mjds:mjde))
      if(.not.allocated(wrfgrid%pblh)) allocate(wrfgrid%pblh(mids:mide,mjds:mjde))
      if(.not.allocated(wrfgrid%u10)) allocate(wrfgrid%u10(mids:mide,mjds:mjde))
      if(.not.allocated(wrfgrid%v10)) allocate(wrfgrid%v10(mids:mide,mjds:mjde))
      if(.not.allocated(wrfgrid%t2)) allocate(wrfgrid%t2(mids:mide,mjds:mjde))
      if(.not.allocated(wrfgrid%q2)) allocate(wrfgrid%q2(mids:mide,mjds:mjde))
      if(.not.allocated(wrfgrid%rh2)) allocate(wrfgrid%rh2(mids:mide,mjds:mjde))
      if(.not.allocated(wrfgrid%xland)) allocate(wrfgrid%xland(mids:mide,mjds:mjde))
      if(.not.allocated(wrfgrid%sst)) allocate(wrfgrid%sst(mids:mide,mjds:mjde))
      if(.not.allocated(wrfgrid%xice)) allocate(wrfgrid%xice(mids:mide,mjds:mjde))
      if(.not.allocated(wrfgrid%tsk)) allocate(wrfgrid%tsk(mids:mide,mjds:mjde))
      if(.not.allocated(wrfgrid%snowh)) allocate(wrfgrid%snowh(mids:mide,mjds:mjde))      
      if(.not.allocated(wrfgrid%ust)) allocate(wrfgrid%ust(mids:mide,mjds:mjde))
      if(.not.allocated(wrfgrid%rmol)) allocate(wrfgrid%rmol(mids:mide,mjds:mjde))
      if(.not.allocated(wrfgrid%raincv)) allocate(wrfgrid%raincv(mids:mide,mjds:mjde))
      if(.not.allocated(wrfgrid%rainncv)) allocate(wrfgrid%rainncv(mids:mide,mjds:mjde))
      if(.not.allocated(wrfgrid%swdown)) allocate(wrfgrid%swdown(mids:mide,mjds:mjde))
      if(.not.allocated(wrfgrid%clflo)) allocate(wrfgrid%clflo(mids:mide,mjds:mjde))
      if(.not.allocated(wrfgrid%clfmi)) allocate(wrfgrid%clfmi(mids:mide,mjds:mjde))
      if(.not.allocated(wrfgrid%clfhi)) allocate(wrfgrid%clfhi(mids:mide,mjds:mjde))  

      if(.not.allocated(wrfgrid%soildepth)) allocate(wrfgrid%soildepth(mids:mide,1:wrfgrid%num_wrfgrid_soil_levels,mjds:mjde))
      if(.not.allocated(wrfgrid%soilthick)) allocate(wrfgrid%soilthick(mids:mide,1:wrfgrid%num_wrfgrid_soil_levels,mjds:mjde))
      if(.not.allocated(wrfgrid%soilt)) allocate(wrfgrid%soilt(mids:mide,1:wrfgrid%num_wrfgrid_soil_levels,mjds:mjde))
      if(.not.allocated(wrfgrid%soilm)) allocate(wrfgrid%soilm(mids:mide,1:wrfgrid%num_wrfgrid_soil_levels,mjds:mjde))

      do j=mjds,mjde
      do i=mids,mide
      wrfgrid%soildepth(i,1,j)=0.05  
      wrfgrid%soildepth(i,2,j)=0.25
      wrfgrid%soildepth(i,3,j)=0.7
      wrfgrid%soildepth(i,4,j)=1.5 
      wrfgrid%soilthick(i,1,j)=0.1
      wrfgrid%soilthick(i,2,j)=0.3
      wrfgrid%soilthick(i,3,j)=0.6
      wrfgrid%soilthick(i,4,j)=1.0 
      end do
      end do

  end subroutine initial_wrfgrid
!-----------------------------------------------------------------------
	
  subroutine final_wrfgrid(wrfgrid)
    type(wrfgrid_c) :: wrfgrid

    if(allocated(wrfgrid%z3d))     deallocate(wrfgrid%z3d)
    if(allocated(wrfgrid%u3d))     deallocate(wrfgrid%u3d)
    if(allocated(wrfgrid%v3d))     deallocate(wrfgrid%v3d)
    if(allocated(wrfgrid%qv3d))   deallocate(wrfgrid%qv3d)
    if(allocated(wrfgrid%qc3d))   deallocate(wrfgrid%qc3d)
    if(allocated(wrfgrid%qi3d))   deallocate(wrfgrid%qi3d)       
    if(allocated(wrfgrid%t3d))   deallocate(wrfgrid%t3d)
    if(allocated(wrfgrid%rh3d))   deallocate(wrfgrid%rh3d)
    if(allocated(wrfgrid%p3d))   deallocate(wrfgrid%p3d)
    if(allocated(wrfgrid%taucldv3d))   deallocate(wrfgrid%taucldv3d)
    if(allocated(wrfgrid%taucldi3d))   deallocate(wrfgrid%taucldi3d)       
    if(allocated(wrfgrid%ps))   deallocate(wrfgrid%ps)
    if(allocated(wrfgrid%q2))   deallocate(wrfgrid%q2)
    if(allocated(wrfgrid%rh2))   deallocate(wrfgrid%rh2)
    if(allocated(wrfgrid%t2))   deallocate(wrfgrid%t2)
    if(allocated(wrfgrid%v10))   deallocate(wrfgrid%v10)
    if(allocated(wrfgrid%u10))   deallocate(wrfgrid%u10)
    if(allocated(wrfgrid%ht))   deallocate(wrfgrid%ht)
    if(allocated(wrfgrid%rmol))   deallocate(wrfgrid%rmol)
    if(allocated(wrfgrid%pblh))   deallocate(wrfgrid%pblh)
    if(allocated(wrfgrid%raincv))   deallocate(wrfgrid%raincv)
    if(allocated(wrfgrid%rainncv))   deallocate(wrfgrid%rainncv)
    if(allocated(wrfgrid%snowh))   deallocate(wrfgrid%snowh)       
    if(allocated(wrfgrid%sst))   deallocate(wrfgrid%sst)
    if(allocated(wrfgrid%xice))   deallocate(wrfgrid%xice)
    if(allocated(wrfgrid%tsk))   deallocate(wrfgrid%tsk)       
    if(allocated(wrfgrid%xland))   deallocate(wrfgrid%xland)
    if(allocated(wrfgrid%swdown))   deallocate(wrfgrid%swdown)
    if(allocated(wrfgrid%clflo))   deallocate(wrfgrid%clflo)
    if(allocated(wrfgrid%clfmi))   deallocate(wrfgrid%clfmi)
    if(allocated(wrfgrid%clfhi))   deallocate(wrfgrid%clfhi)

    if(allocated(wrfgrid%soildepth)) deallocate(wrfgrid%soildepth)
    if(allocated(wrfgrid%soilthick)) deallocate(wrfgrid%soilthick)
    if(allocated(wrfgrid%soilt)) deallocate(wrfgrid%soilt)
    if(allocated(wrfgrid%soilm)) deallocate(wrfgrid%soilm)

  end subroutine final_wrfgrid
!-----------------------------------------------------------------------
	
  subroutine geatm_import_mct( x2c_c, camgrid, loop )     
!-----------------------------------------------------------------------
! Arguments
!
     use geatm_vartype, only:csx,cex,csy,cey,myid	
     TYPE(camgrid_c) :: camgrid
     type(mct_aVect)   , intent(inout) :: x2c_c
     integer :: loop

     ! Local variables	
     integer  :: i,j,k,ig  ! indices    
     integer :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,ips,ipe,jps,jpe,kps,kpe   
     real*8,dimension(:,:),allocatable::xlat,xlon
     real*8,dimension(:),allocatable::templon,templat
     logical,save :: first_sst=.true.     
 
       ips=csx(1)
       ipe=cex(1)
       jps=csy(1)
       jpe=cey(1)
       
      ig=1
      do j=jps,jpe                                                        
       do i =ips, ipe
         
          camgrid%xlat(i,j)       = x2c_c%rAttr(index_x2ge_Sx_lat     ,ig)*180.0/3.1415926
          camgrid%xlon(i,j)       = x2c_c%rAttr(index_x2ge_Sx_lon     ,ig)*180.0/3.1415926 
	  if(camgrid%xlon(i,j).lt.0) camgrid%xlon(i,j) = camgrid%xlon(i,j)+360  
          camgrid%ps(i,j) = x2c_c%rAttr(index_x2ge_Sx_ps      ,ig)
          camgrid%ht(i,j) =  x2c_c%rAttr(index_x2ge_Sx_phis    ,ig)/9.81
          camgrid%t2(i,j) =  x2c_c%rAttr(index_x2ge_Sx_t2    ,ig) 
          camgrid%q2(i,j) =  x2c_c%rAttr(index_x2ge_Sx_q2    ,ig)
          camgrid%pblh(i,j) =  x2c_c%rAttr(index_x2ge_Sx_pblh    ,ig)
          camgrid%u10(i,j) =  x2c_c%rAttr(index_x2ge_Sx_u10    ,ig)
          camgrid%v10(i,j) =  x2c_c%rAttr(index_x2ge_Sx_v10    ,ig)
          camgrid%ust(i,j) =  x2c_c%rAttr(index_x2ge_Sx_ust    ,ig)
          camgrid%rmol(i,j) =  x2c_c%rAttr(index_x2ge_Sx_rmol    ,ig)
          camgrid%raincv(i,j) =  x2c_c%rAttr(index_x2ge_Sx_raincv    ,ig)
          camgrid%rainncv(i,j) =  x2c_c%rAttr(index_x2ge_Sx_rainncv    ,ig)
          camgrid%swdown(i,j) =  x2c_c%rAttr(index_x2ge_Sx_swdown    ,ig)
          camgrid%clflo(i,j) =  x2c_c%rAttr(index_x2ge_Sx_clflo    ,ig)
          camgrid%clfmi(i,j) =  x2c_c%rAttr(index_x2ge_Sx_clfmi    ,ig)
          camgrid%clfhi(i,j) =  x2c_c%rAttr(index_x2ge_Sx_clfhi    ,ig)
          camgrid%rh2(i,j) =  x2c_c%rAttr(index_x2ge_Sx_rh2    ,ig)
          camgrid%sst(i,j) = x2c_c%rAttr(index_x2ge_Sx_sst  ,ig)    
          camgrid%xice(i,j) = x2c_c%rAttr(index_x2ge_Sx_seaice  ,ig)
          camgrid%tsk(i,j) = x2c_c%rAttr(index_x2ge_Sx_ts  ,ig)      
          camgrid%xland(i,j)= 1.0-x2c_c%rAttr(index_x2ge_Sx_ocnfrac  ,ig) ! used as land
          camgrid%snowh(i,j) = (camgrid%xland(i,j) *x2c_c%rAttr(index_x2ge_Sx_snowhland  ,ig)  +camgrid%xice(i,j) * x2c_c%rAttr(index_x2ge_Sx_snowhice   ,ig))
           if(camgrid%xland(i,j).lt.0.5) then
            camgrid%xland(i,j)=0
           else
            camgrid%xland(i,j)=1
           endif
         
          ! from surface to underground 
          do k=1, camgrid%num_camgrid_soil_levels
            camgrid%soildepth(i,camgrid%num_camgrid_soil_levels-k+1,j)       =  x2c_c%rAttr(index_x2ge_Sx_soildepth(k),ig)  ! The unit is supposed as cm
            camgrid%soilthick(i,camgrid%num_camgrid_soil_levels-k+1,j)       = x2c_c%rAttr(index_x2ge_Sx_soilthick(k),ig)   !The unit is supposed as cm
            camgrid%soilt(i,camgrid%num_camgrid_soil_levels-k+1,j)   = x2c_c%rAttr(index_x2ge_Sx_soilt(k),ig)
            camgrid%soilm(i,camgrid%num_camgrid_soil_levels-k+1,j)   = x2c_c%rAttr(index_x2ge_Sx_soilm(k),ig) 
          enddo
           
	  do k = 1, camgrid%num_camgrid_levels
            camgrid%u3d(i,k,j) = x2c_c%rAttr(index_x2ge_Sx_u3d(k)     ,ig)
            camgrid%v3d(i,k,j) = x2c_c%rAttr(index_x2ge_Sx_v3d(k)     ,ig)
            camgrid%t3d(i,k,j) = x2c_c%rAttr(index_x2ge_Sx_t3d(k)     ,ig)
            camgrid%z3d(i,k,j) = x2c_c%rAttr(index_x2ge_Sx_z3d(k)     ,ig)            
            camgrid%p3d(i,k,j) = x2c_c%rAttr(index_x2ge_Sx_p3d(k)     ,ig)
	    camgrid%rh3d(i,k,j) = x2c_c%rAttr(index_x2ge_Sx_rh3d(k)     ,ig)  
            camgrid%qv3d(i,k,j)= x2c_c%rAttr(index_x2ge_Sx_qv3d(k),ig) 
            camgrid%qc3d(i,k,j)= x2c_c%rAttr(index_x2ge_Sx_qc3d(k),ig)
            camgrid%qi3d(i,k,j)= x2c_c%rAttr(index_x2ge_Sx_qi3d(k),ig)
            camgrid%taucldv3d(i,k,j) = x2c_c%rAttr(index_x2ge_Sx_taucldv3d(k)     ,ig)
            camgrid%taucldi3d(i,k,j) = x2c_c%rAttr(index_x2ge_Sx_taucldi3d(k)     ,ig)            
          end do  

           if( camgrid%ht(i,j) .lt.0)   camgrid%ht(i,j) = 0 ! it may be littler than zero in the ocean for the terrain
           do  k = 1, camgrid%num_camgrid_levels
              camgrid%z3d(i,k,j)  = camgrid%z3d(i,k,j)+camgrid%ht(i,j) ! the original height is the height above the surface not the sea level
           enddo

          ig=ig+1
         end do
       end do                                   

       ips=csx(1)
       ipe=cex(1)
       jps=csy(1)
       jpe=cey(1)
       kps=1
       kpe=camgrid%num_camgrid_levels       
       ids=1
       ide=camgrid%num_camgrid_lon
       jds=1
       jde=camgrid%num_camgrid_lat
       kds=1
       kde=camgrid%num_camgrid_levels
       call garther_camgrid(camgrid, ids, ide, jds, jde, kds, kde, &
                            camgrid%num_camgrid_levels, camgrid%num_camgrid_soil_levels, &
                            ips, ipe, jps, jpe, kps, kpe)

       if(first_sst) then
       call garther_camgrid_once(camgrid, ids, ide, jds, jde, kds, kde, &
                            camgrid%num_camgrid_levels,camgrid%num_camgrid_soil_levels, &
                            ips, ipe, jps, jpe, kps, kpe)
       allocate(xlon(wrfgrid%num_wrfgrid_lon,wrfgrid%num_wrfgrid_lat))
       allocate(xlat(wrfgrid%num_wrfgrid_lon,wrfgrid%num_wrfgrid_lat))
       allocate(templat(camgrid%num_camgrid_lat))
       allocate(templon(camgrid%num_camgrid_lon))
       do j=1,wrfgrid%num_wrfgrid_lat
       do i=1,wrfgrid%num_wrfgrid_lon
       xlon(i,j)=(i-1)*1.0_8+0.5_8
       end do
       end do
       do j=1,wrfgrid%num_wrfgrid_lat
       do i=1,wrfgrid%num_wrfgrid_lon
       xlat(i,j)=(j-1)*1.0_8-89.5_8
       end do       
       end do       
       templon(:)=1.0_8*camgrid%xlon(:,1)
       templat(:)=1.0_8*camgrid%xlat(1,:)
       call geatm_to_cam_mapping(templon, templat, camgrid%num_camgrid_lon, camgrid%num_camgrid_lat, &
            xlon, xlat,wrfgrid%num_wrfgrid_lon,wrfgrid%num_wrfgrid_lat)
       deallocate(xlat)
       deallocate(xlon)
       deallocate(templat)
       deallocate(templon)
        first_sst=.false.
       endif

       call cam_to_geatm(camgrid, wrfgrid)
       call wrfgrid_to_geatm(wrfgrid)

  end subroutine geatm_import_mct
!=============================================================================================================

  subroutine geatm_export_mct( x2ge_ge, camgrid, loop )
!-----------------------------------------------------------------------
! Arguments
!
     use geatm_vartype, only:csx,cex,csy,cey    
     TYPE(camgrid_c) :: camgrid
     type(mct_aVect)   , intent(inout) :: x2ge_ge
     integer :: loop
 end subroutine geatm_export_mct
!=============================================================================================================

  subroutine cam_to_geatm(camgrid, wrfgrid )
     use geatm_vartype, only:csx,cex,csy,cey,sx,ex,sy,ey

     TYPE(camgrid_c) :: camgrid
     TYPE(wrfgrid_c) :: wrfgrid

     ! Local variables  
     integer  :: i,j,k,ix,jy,km,local_index,global_index1,global_index2
     integer :: mids, mide, mjds, mjde, mkds, mkde, &
                mims, mime, mjms, mjme, mkms, mkme, &
                mips, mipe, mjps, mjpe, mkps, mkpe
     integer ::   ids, ide, jds, jde, kds, kde, &
                  ims, ime, jms, jme, kms, kme, &
                  ips, ipe, jps, jpe, kps, kpe 
     real(8),dimension(:,:,:),allocatable :: t3d_old
     real(8),dimension(:,:,:),allocatable :: t3d_new

       ! p3d 
       call cam_interpolate_to_geatm(camgrid, wrfgrid, 'p3d' )
       
       ! u3d              
       call cam_interpolate_to_geatm(camgrid, wrfgrid, 'u3d' )
       
       ! v3d              
       call cam_interpolate_to_geatm(camgrid, wrfgrid, 'v3d' )

       ! rh3d
       call cam_interpolate_to_geatm(camgrid, wrfgrid, 'rh3d' )

       ! qv3d
       call cam_interpolate_to_geatm(camgrid, wrfgrid, 'qv3d' )

       ! qc3d
       call cam_interpolate_to_geatm(camgrid, wrfgrid, 'qc3d' )

       ! qi3d
       call cam_interpolate_to_geatm(camgrid, wrfgrid, 'qi3d' )

       ! t3d
       call cam_interpolate_to_geatm(camgrid, wrfgrid, 't3d' )

       ! taucldv3d
       call cam_interpolate_to_geatm(camgrid, wrfgrid, 'taucldv3d' ) 

       ! taucldi3d
       call cam_interpolate_to_geatm(camgrid, wrfgrid, 'taucldi3d' )

       ! soilt
       call cam_interpolate_to_geatm_soil(camgrid, wrfgrid, 'soilt' )

       ! soilm
       call cam_interpolate_to_geatm_soil(camgrid, wrfgrid, 'soilm' )

       ips=sx(1)
       ipe=ex(1)
       jps=sy(1)
       jpe=ey(1)
       ids=1
       ide=wrfgrid%num_wrfgrid_lon
       jds=1
       jde=wrfgrid%num_wrfgrid_lat
       mips=csx(1)
       mipe=cex(1)
       mjps=csy(1)
       mjpe=cey(1)
       mids=1
       mide=camgrid%num_camgrid_lon
       mjds=1
       mjde=camgrid%num_camgrid_lat

       allocate(t3d_old(mids:mide,1:1,mjds:mjde))
       allocate(t3d_new(ips:ipe,1:1,jps:jpe))

       ! ps
       do j=mjds,mjde
       do i=mids,mide
       t3d_old(i,1,j)=camgrid%ps(i,j)
       end do
       end do
       call cam_horizontal_interpolate_to_geatm(t3d_old,t3d_new, &
       mips, mipe, mjps, mjpe, 1, 1, &
       mids, mide, mjds, mjde, 1, 1, &
       ips, ipe, jps, jpe, 1, 1, &
       ids, ide, jds, jde, 1, 1)
       do j=jps,jpe
       do i=ips,ipe
       wrfgrid%ps(i,j)=t3d_new(i,1,j)
       end do
       end do

       ! snowh
       do j=mjds,mjde
       do i=mids,mide
       t3d_old(i,1,j)=camgrid%snowh(i,j)
       end do
       end do
       call cam_horizontal_interpolate_to_geatm(t3d_old,t3d_new, &
       mips, mipe, mjps, mjpe, 1, 1, &
       mids, mide, mjds, mjde, 1, 1, &
       ips, ipe, jps, jpe, 1, 1, &
       ids, ide, jds, jde, 1, 1)
       do j=jps,jpe
       do i=ips,ipe
       wrfgrid%snowh(i,j)=t3d_new(i,1,j)
       end do
       end do

       ! xice
       do j=mjds,mjde
       do i=mids,mide
       t3d_old(i,1,j)=camgrid%xice(i,j)*1.0_8
       end do
       end do
       call cam_horizontal_interpolate_to_geatm(t3d_old,t3d_new, &
       mips, mipe, mjps, mjpe, 1, 1, &
       mids, mide, mjds, mjde, 1, 1, &
       ips, ipe, jps, jpe, 1, 1, &
       ids, ide, jds, jde, 1, 1)
       do j=jps,jpe
       do i=ips,ipe
       wrfgrid%xice(i,j)=t3d_new(i,1,j)
       end do
       end do

       ! rh2
       do j=mjds,mjde
       do i=mids,mide
       t3d_old(i,1,j)=camgrid%rh2(i,j)
       end do
       end do
       call cam_horizontal_interpolate_to_geatm(t3d_old,t3d_new, &
       mips, mipe, mjps, mjpe, 1, 1, &
       mids, mide, mjds, mjde, 1, 1, &
       ips, ipe, jps, jpe, 1, 1, &
       ids, ide, jds, jde, 1, 1)
       do j=jps,jpe
       do i=ips,ipe
       wrfgrid%rh2(i,j)=t3d_new(i,1,j)
       end do
       end do

       ! q2
       do j=mjds,mjde
       do i=mids,mide
       t3d_old(i,1,j)=camgrid%q2(i,j)
       end do
       end do
       call cam_horizontal_interpolate_to_geatm(t3d_old,t3d_new, &
       mips, mipe, mjps, mjpe, 1, 1, &   
       mids, mide, mjds, mjde, 1, 1, &
       ips, ipe, jps, jpe, 1, 1, &
       ids, ide, jds, jde, 1, 1)
       do j=jps,jpe
       do i=ips,ipe
       wrfgrid%q2(i,j)=t3d_new(i,1,j)
       end do
       end do

       ! t2
       do j=mjds,mjde
       do i=mids,mide
       t3d_old(i,1,j)=camgrid%t2(i,j)
       end do
       end do
       call cam_horizontal_interpolate_to_geatm(t3d_old,t3d_new, &
       mips, mipe, mjps, mjpe, 1, 1, &   
       mids, mide, mjds, mjde, 1, 1, &
       ips, ipe, jps, jpe, 1, 1, &
       ids, ide, jds, jde, 1, 1)
       do j=jps,jpe
       do i=ips,ipe
       wrfgrid%t2(i,j)=t3d_new(i,1,j)
       end do
       end do

       ! u10
       do j=mjds,mjde
       do i=mids,mide
       t3d_old(i,1,j)=camgrid%u10(i,j)
       end do
       end do
       call cam_horizontal_interpolate_to_geatm(t3d_old,t3d_new, &
       mips, mipe, mjps, mjpe, 1, 1, &   
       mids, mide, mjds, mjde, 1, 1, &
       ips, ipe, jps, jpe, 1, 1, &
       ids, ide, jds, jde, 1, 1)
       do j=jps,jpe
       do i=ips,ipe
       wrfgrid%u10(i,j)=t3d_new(i,1,j)
       end do
       end do

       ! v10
       do j=mjds,mjde
       do i=mids,mide
       t3d_old(i,1,j)=camgrid%v10(i,j)
       end do
       end do
       call cam_horizontal_interpolate_to_geatm(t3d_old,t3d_new, &
       mips, mipe, mjps, mjpe, 1, 1, &   
       mids, mide, mjds, mjde, 1, 1, &
       ips, ipe, jps, jpe, 1, 1, &
       ids, ide, jds, jde, 1, 1)
       do j=jps,jpe
       do i=ips,ipe
       wrfgrid%v10(i,j)=t3d_new(i,1,j)
       end do
       end do

       ! pblh
       do j=mjds,mjde
       do i=mids,mide
       t3d_old(i,1,j)=camgrid%pblh(i,j)
       end do
       end do
       call cam_horizontal_interpolate_to_geatm(t3d_old,t3d_new, &
       mips, mipe, mjps, mjpe, 1, 1, &   
       mids, mide, mjds, mjde, 1, 1, &
       ips, ipe, jps, jpe, 1, 1, &
       ids, ide, jds, jde, 1, 1)
       do j=jps,jpe
       do i=ips,ipe
       wrfgrid%pblh(i,j)=t3d_new(i,1,j)
       end do
       end do
       
       ! rmol
       do j=mjds,mjde
       do i=mids,mide
       t3d_old(i,1,j)=camgrid%rmol(i,j)
       end do
       end do
       call cam_horizontal_interpolate_to_geatm(t3d_old,t3d_new, &
       mips, mipe, mjps, mjpe, 1, 1, &   
       mids, mide, mjds, mjde, 1, 1, &
       ips, ipe, jps, jpe, 1, 1, &
       ids, ide, jds, jde, 1, 1)
       do j=jps,jpe
       do i=ips,ipe
       wrfgrid%rmol(i,j)=t3d_new(i,1,j)
       end do
       end do

       ! ust
       do j=mjds,mjde
       do i=mids,mide
       t3d_old(i,1,j)=camgrid%ust(i,j)
       end do
       end do
       call cam_horizontal_interpolate_to_geatm(t3d_old,t3d_new, &
       mips, mipe, mjps, mjpe, 1, 1, &   
       mids, mide, mjds, mjde, 1, 1, &
       ips, ipe, jps, jpe, 1, 1, &
       ids, ide, jds, jde, 1, 1)
       do j=jps,jpe
       do i=ips,ipe
       wrfgrid%ust(i,j)=t3d_new(i,1,j)
       end do
       end do
      
       if(1.eq.0) then 
       ! no use of sst and tsk during the integration
       ! sst
       do j=mjds,mjde
       do i=mids,mide
       t3d_old(i,1,j)=camgrid%sst(i,j)
       end do
       end do
       call cam_horizontal_interpolate_to_geatm(t3d_old,t3d_new, &
       mips, mipe, mjps, mjpe, 1, 1, &   
       mids, mide, mjds, mjde, 1, 1, &
       ips, ipe, jps, jpe, 1, 1, &
       ids, ide, jds, jde, 1, 1)
       do j=jps,jpe
       do i=ips,ipe
       wrfgrid%sst(i,j)=t3d_new(i,1,j)
       end do
       end do
      
       ! tsk
       do j=mjds,mjde
       do i=mids,mide
       t3d_old(i,1,j)=camgrid%tsk(i,j)
       end do
       end do
       call cam_horizontal_interpolate_to_geatm(t3d_old,t3d_new, &
       mips, mipe, mjps, mjpe, 1, 1, &   
       mids, mide, mjds, mjde, 1, 1, &
       ips, ipe, jps, jpe, 1, 1, &
       ids, ide, jds, jde, 1, 1)
       do j=jps,jpe
       do i=ips,ipe
       wrfgrid%tsk(i,j)=t3d_new(i,1,j)
       end do
       end do
       end if

       ! rainncv
       do j=mjds,mjde
       do i=mids,mide
       t3d_old(i,1,j)=camgrid%rainncv(i,j)
       end do
       end do
       call cam_horizontal_interpolate_to_geatm(t3d_old,t3d_new, &
       mips, mipe, mjps, mjpe, 1, 1, &   
       mids, mide, mjds, mjde, 1, 1, &
       ips, ipe, jps, jpe, 1, 1, &
       ids, ide, jds, jde, 1, 1)
       do j=jps,jpe
       do i=ips,ipe
       wrfgrid%rainncv(i,j)=t3d_new(i,1,j)
       end do
       end do

       ! raincv
       do j=mjds,mjde
       do i=mids,mide
       t3d_old(i,1,j)=camgrid%raincv(i,j)
       end do
       end do
       call cam_horizontal_interpolate_to_geatm(t3d_old,t3d_new, &
       mips, mipe, mjps, mjpe, 1, 1, &   
       mids, mide, mjds, mjde, 1, 1, &
       ips, ipe, jps, jpe, 1, 1, &
       ids, ide, jds, jde, 1, 1)
       do j=jps,jpe
       do i=ips,ipe
       wrfgrid%raincv(i,j)=t3d_new(i,1,j)
       end do
       end do
       
       ! swdown
       do j=mjds,mjde
       do i=mids,mide
       t3d_old(i,1,j)=camgrid%swdown(i,j)
       end do
       end do
       call cam_horizontal_interpolate_to_geatm(t3d_old,t3d_new, &
       mips, mipe, mjps, mjpe, 1, 1, &   
       mids, mide, mjds, mjde, 1, 1, &
       ips, ipe, jps, jpe, 1, 1, &
       ids, ide, jds, jde, 1, 1)
       do j=jps,jpe
       do i=ips,ipe
       wrfgrid%swdown(i,j)=t3d_new(i,1,j)
       end do
       end do
       
       ! clflo
       do j=mjds,mjde
       do i=mids,mide
       t3d_old(i,1,j)=camgrid%clflo(i,j)
       end do
       end do
       call cam_horizontal_interpolate_to_geatm(t3d_old,t3d_new, &
       mips, mipe, mjps, mjpe, 1, 1, &   
       mids, mide, mjds, mjde, 1, 1, &
       ips, ipe, jps, jpe, 1, 1, &
       ids, ide, jds, jde, 1, 1)
       do j=jps,jpe
       do i=ips,ipe
       wrfgrid%clflo(i,j)=t3d_new(i,1,j)
       end do
       end do

       ! clfmi
       do j=mjds,mjde
       do i=mids,mide
       t3d_old(i,1,j)=camgrid%clfmi(i,j)
       end do
       end do
       call cam_horizontal_interpolate_to_geatm(t3d_old,t3d_new, &
       mips, mipe, mjps, mjpe, 1, 1, &   
       mids, mide, mjds, mjde, 1, 1, &
       ips, ipe, jps, jpe, 1, 1, &
       ids, ide, jds, jde, 1, 1)
       do j=jps,jpe
       do i=ips,ipe
       wrfgrid%clfmi(i,j)=t3d_new(i,1,j)
       end do
       end do

       ! clfhi
       do j=mjds,mjde
       do i=mids,mide
       t3d_old(i,1,j)=camgrid%clfhi(i,j)
       end do
       end do
       call cam_horizontal_interpolate_to_geatm(t3d_old,t3d_new, &
       mips, mipe, mjps, mjpe, 1, 1, &   
       mids, mide, mjds, mjde, 1, 1, &
       ips, ipe, jps, jpe, 1, 1, &
       ids, ide, jds, jde, 1, 1)
       do j=jps,jpe
       do i=ips,ipe
       wrfgrid%clfhi(i,j)=t3d_new(i,1,j)
       end do
       end do

      deallocate(t3d_old)
      deallocate(t3d_new)

    end subroutine cam_to_geatm
    
!=============================================================================================================

    subroutine cam_interpolate_to_geatm(camgrid, wrfgrid, symbol )
     use geatm_vartype, only:csx,cex,csy,cey,sx,ex,sy,ey

     TYPE(camgrid_c) :: camgrid
     TYPE(wrfgrid_c) :: wrfgrid
     character(len=*) :: symbol
     
     ! Local variables  
     integer  :: i,j,k,ix,jy,kz,km,local_index,global_index1,global_index2
     real(8) :: s, wweight
     real(8),dimension(:,:,:),allocatable:: p3d, t3d, w3d, t3d_new, t3d_int 
     integer ::   ids, ide, jds, jde, kds, kde, &
                  ips, ipe, jps, jpe, kps, kpe 
                   
      ips=sx(1)
      ipe=ex(1)
      jps=sy(1)
      jpe=ey(1)

      allocate(p3d(1,1:camgrid%num_camgrid_levels,1))
      allocate(t3d(1,1:camgrid%num_camgrid_levels,1))      
      allocate(w3d(1,1:wrfgrid%num_wrfgrid_levels,1))
      allocate(t3d_new(1,1:wrfgrid%num_wrfgrid_levels,1))
      allocate(t3d_int(ips:ipe,1:wrfgrid%num_wrfgrid_levels,jps:jpe))

       global_index1 = 0
       t3d_int = 0.0_8
       do jy=jps,jpe
       do ix=ips,ipe

        do kz=1,wrfgrid%num_wrfgrid_levels
        w3d(1,kz,1) = wrfgrid%z3d(ix,kz,jy)
        end do

        local_index=0
        local_index=cam_local_points(ix,jy)
        if(local_index.gt.0) then

        do km=1,local_index
           global_index1=global_index1+1
           i=remap_cami(global_index1)
           j=remap_camj(global_index1)
           wweight=cam_weight(global_index1)
         
          select case(trim(adjustl(symbol)))
          ! p3d
          case('p3d')
            do kz=1,camgrid%num_camgrid_levels
             t3d(1,kz,1) = camgrid%p3d(i,kz,j)
            end do
          ! u3d
          case('u3d')
            do kz=1,camgrid%num_camgrid_levels
             t3d(1,kz,1) = camgrid%u3d(i,kz,j)
            end do
          ! v3d
          case('v3d')
            do kz=1,camgrid%num_camgrid_levels
             t3d(1,kz,1) = camgrid%v3d(i,kz,j)
            end do
          ! t3d
          case('t3d')
            do kz=1,camgrid%num_camgrid_levels
             t3d(1,kz,1) = camgrid%t3d(i,kz,j)
            end do
          ! qv3d
          case('qv3d')
            do kz=1,camgrid%num_camgrid_levels
             t3d(1,kz,1) = camgrid%qv3d(i,kz,j)
            end do
          ! qc3d
          case('qc3d')
            do kz=1,camgrid%num_camgrid_levels
             t3d(1,kz,1) = camgrid%qc3d(i,kz,j)
            end do
          ! qi3d
          case('qi3d')
            do kz=1,camgrid%num_camgrid_levels
             t3d(1,kz,1) = camgrid%qi3d(i,kz,j)
            end do 
          ! rh3d
          case('rh3d')
            do kz=1,camgrid%num_camgrid_levels
             t3d(1,kz,1) = camgrid%rh3d(i,kz,j)
            end do
          ! taucldv3d
          case('taucldv3d')
            do kz=1,camgrid%num_camgrid_levels
             t3d(1,kz,1)=camgrid%taucldv3d(i,kz,j)
            end do
          ! taucldi3d
          case('taucldi') 
            do kz=1,camgrid%num_camgrid_levels
             t3d(1,kz,1)=camgrid%taucldi3d(i,kz,j)
            end do
         end select       

         do kz=1,camgrid%num_camgrid_levels
         p3d(1,kz,1) = camgrid%z3d(i,kz,j)
         end do

         call interpolate_to_pressure_levels(1,camgrid%num_camgrid_levels &
               ,1,1,wrfgrid%num_wrfgrid_levels &
               ,p3d,w3d,t3d,t3d_new)
          
           do kz=1,wrfgrid%num_wrfgrid_levels
             t3d_int(ix,kz,jy) = t3d_int(ix,kz,jy)+t3d_new(1,kz,1)*wweight
           end do
      
       enddo
       end if ! local index
    
          select case(trim(adjustl(symbol)))
          ! p3d
          case('p3d')
           do kz=1,wrfgrid%num_wrfgrid_levels
            wrfgrid%p3d(ix,kz,jy) = t3d_int(ix,kz,jy)
           end do
          ! u3d
          case('u3d')
           do kz=1,wrfgrid%num_wrfgrid_levels
             wrfgrid%u3d(ix,kz,jy) = t3d_int(ix,kz,jy)
           end do 
          ! v3d
          case('v3d')
           do kz=1,wrfgrid%num_wrfgrid_levels
            wrfgrid%v3d(ix,kz,jy) = t3d_int(ix,kz,jy)
           end do
          ! t3d
          case('t3d')
           do kz=1,wrfgrid%num_wrfgrid_levels
            wrfgrid%t3d(ix,kz,jy) = t3d_int(ix,kz,jy)
           end do
          ! qv3d
          case('qv3d')
           do kz=1,wrfgrid%num_wrfgrid_levels
             wrfgrid%qv3d(ix,kz,jy) = t3d_int(ix,kz,jy)
           end do
          ! qc3d
          case('qc3d')
           do kz=1,wrfgrid%num_wrfgrid_levels
            wrfgrid%qc3d(ix,kz,jy) = t3d_int(ix,kz,jy)
           end do
          ! qi3d
          case('qi3d')
           do kz=1,wrfgrid%num_wrfgrid_levels
            wrfgrid%qi3d(ix,kz,jy) = t3d_int(ix,kz,jy)
           end do
          ! rh3d
          case('rh3d')
           do kz=1,wrfgrid%num_wrfgrid_levels
            wrfgrid%rh3d(ix,kz,jy) = t3d_int(ix,kz,jy)
           end do
          ! taucldv3d
          case('taucldv3d')
           do kz=1,wrfgrid%num_wrfgrid_levels
            wrfgrid%taucldv3d(ix,kz,jy) = t3d_int(ix,kz,jy)
           end do
          ! taucldi3d
          case('taucldi') 
           do kz=1,wrfgrid%num_wrfgrid_levels
             wrfgrid%taucldi3d(ix,kz,jy) = t3d_int(ix,kz,jy)
           end do
         end select       
       
      enddo 
      enddo 

      deallocate(p3d)
      deallocate(t3d)      
      deallocate(w3d)
      deallocate(t3d_new)
      deallocate(t3d_int)

      end subroutine cam_interpolate_to_geatm

!=============================================================================================================
    subroutine cam_interpolate_to_geatm_soil(camgrid, wrfgrid, symbol )
     use geatm_vartype, only:csx,cex,csy,cey,sx,ex,sy,ey

     TYPE(camgrid_c) :: camgrid
     TYPE(wrfgrid_c) :: wrfgrid
     character(len=*) :: symbol
     
     ! Local variables  
     integer  :: i,j,k,ix,jy,km,local_index,global_index1,global_index2
     real(8) :: s, wweight
     real(8),dimension(:,:,:),allocatable:: p3d, t3d, w3d, t3d_new,t3d_int 
     integer ::   ids, ide, jds, jde, kds, kde, &
                  ips, ipe, jps, jpe, kps, kpe 

      ips=sx(1)
      ipe=ex(1)
      jps=sy(1)
      jpe=ey(1)

      allocate(p3d(1,1:camgrid%num_camgrid_soil_levels,1))
      allocate(t3d(1,1:camgrid%num_camgrid_soil_levels,1))      
      allocate(w3d(1,1:wrfgrid%num_wrfgrid_soil_levels,1))
      allocate(t3d_new(1,1:wrfgrid%num_wrfgrid_soil_levels,1))
      allocate(t3d_int(ips:ipe,1:wrfgrid%num_wrfgrid_soil_levels,jps:jpe))

       global_index1 = 0
       t3d_int = 0
       do jy=jps,jpe
       do ix=ips,ipe

         do k=1,wrfgrid%num_wrfgrid_soil_levels
           w3d(1,k,1) = wrfgrid%soildepth(ix,wrfgrid%num_wrfgrid_soil_levels-k+1,jy)  ! use m
         enddo

        local_index=0
        local_index=cam_local_points(ix,jy)
        if(local_index.gt.0) then
        do km=1,local_index
           global_index1=global_index1+1
           i=remap_cami(global_index1)
           j=remap_camj(global_index1)
           wweight=cam_weight(global_index1)

         select case(trim(adjustl(symbol)))
         ! soilt
         case('soilt')
            do k=1,camgrid%num_camgrid_soil_levels
             t3d(1,k,1) = camgrid%soilt(i,k,j)
            end do
         ! soilm
         case('soilm')
            do k=1,camgrid%num_camgrid_soil_levels
             t3d(1,:,1) = camgrid%soilm(i,:,j)
            end do
         end select

         do k=1,camgrid%num_camgrid_soil_levels
           p3d(1,k,1) = camgrid%soildepth(i,k,j)
         end do
         call interpolate_to_pressure_levels(1,camgrid%num_camgrid_soil_levels &
               ,1,1,wrfgrid%num_wrfgrid_soil_levels &
               ,p3d,w3d,t3d,t3d_new)

         do k=1,wrfgrid%num_wrfgrid_soil_levels      
          t3d_int(ix,k,jy) = t3d_int(ix,k,jy)+t3d_new(1,k,1)*wweight
         end do 

       enddo
       end if ! local index
    
          ! sort to from the surfce to the underground
          do k=1,wrfgrid%num_wrfgrid_soil_levels/2 
           s = t3d_int(ix,k,jy)
           t3d_int(ix,k,jy) = t3d_int(ix,wrfgrid%num_wrfgrid_soil_levels-k+1,jy)
           t3d_int(ix,wrfgrid%num_wrfgrid_soil_levels-k+1,jy) = t3d_int(ix,k,jy)
          end do

         select case(trim(adjustl(symbol)))
         ! soilt
         case('soilt')
           do k=1,wrfgrid%num_wrfgrid_soil_levels
            wrfgrid%soilt(ix,k,jy) = t3d_int(ix,k,jy)
           end do
         ! soilm
         case('soilm')
           do k=1,wrfgrid%num_wrfgrid_soil_levels
            wrfgrid%soilm(ix,k,jy) = t3d_int(ix,k,jy)
           end do
         end select
       
      enddo 
      enddo 

      deallocate(p3d)
      deallocate(t3d)      
      deallocate(w3d)
      deallocate(t3d_new)
      deallocate(t3d_int)

      end subroutine cam_interpolate_to_geatm_soil
!=============================================================================================================

    subroutine cam_horizontal_interpolate_to_geatm(t3d_old,t3d_new, &
       mips, mipe, mjps, mjpe, mkps, mkpe, &   
       mids, mide, mjds, mjde, mkds, mkde, &
       ips, ipe, jps, jpe, kps, kpe, &
       ids, ide, jds, jde, kds, kde)

     integer ::  mids, mide, mjds, mjde, mkds, mkde, &
                  mips, mipe, mjps, mjpe, mkps, mkpe, & 
                  ids, ide, jds, jde, kds, kde, &
                  ips, ipe, jps, jpe, kps, kpe 
     real(8),dimension(mids:mide,1:1,mjds:mjde) :: t3d_old                   
     real(8),dimension(ips:ipe,1:1,jps:jpe) :: t3d_new
     
     ! Local variables  
     integer  :: i,j,k,ix,jy,km,local_index,global_index1,global_index2
     real(8) :: s, wweight 


       t3d_new = 0
       global_index1 = 0
       do jy=jps,jpe
       do ix=ips,ipe

        local_index=0
        local_index=cam_local_points(ix,jy)
        if(local_index.gt.0) then
        do km=1,local_index
           global_index1=global_index1+1
           i=remap_cami(global_index1)
           j=remap_camj(global_index1)
           wweight=cam_weight(global_index1)
           t3d_new(ix,1,jy) = t3d_new(ix,1,jy)+t3d_old(i,1,j)*wweight
        enddo
        end if ! local index
       
      enddo 
      enddo 

      end subroutine cam_horizontal_interpolate_to_geatm
      
!=============================================================================================================      
     subroutine wrfgrid_to_geatm(wrfgrid )
     use geatm_vartype, only: nx,ny,sx,ex,sy,ey, &
                              ip3mem, ip2mem, plev, heiz, &
                              u, v, t, qvapor, rh1, taucldc, &
                              taucldi, q2, t2, u10, h,clw, rnw, &
                              v10, rmol, psfc, rhsfc, &
                              rainnon, raincon, swdown, ust, &
                              clflo, clfmi, clfhi, soilt, soilrh, &
                              fsnow, fice, pbl_hgt

     TYPE(wrfgrid_c) :: wrfgrid

     ! Local variables  
     integer  :: i,j,k,ix,jy,km,kl,i0
       
       do k=1,wrfgrid%num_wrfgrid_levels
        i0=ip3mem(k,1)
        kl=wrfgrid%num_wrfgrid_levels-k+1 ! from the surface to upper troposphere
          do j= sy(1),ey(1)
          do i= sx(1),ex(1)
            km=i0+(ex(1)-sx(1)+3)*(j-sy(1)+1)+i-sx(1)+1
            plev(km)=real(wrfgrid%p3d(i,kl,j)/100.0) ! from Pa to hPa
            h(km)=real(wrfgrid%z3d(i,kl,j)/1000.0)  ! from m to km
            t(km)=real(wrfgrid%t3d(i,kl,j))
            u(km)=real(wrfgrid%u3d(i,kl,j))
            v(km)=real(wrfgrid%v3d(i,kl,j)) 
            qvapor(km)=real(wrfgrid%qv3d(i,kl,j)) 
            clw(km)=real(wrfgrid%qc3d(i,kl,j)) 
            rnw(km)=real(wrfgrid%qi3d(i,kl,j))
            rh1(km)=real(wrfgrid%rh3d(i,kl,j))
            if(rh1(km).gt.100.0) then
            rh1(km)=100.0
            elseif(rh1(km).lt.0.0) then
            rh1(km)=0.0
            endif
            taucldc(km)=real(wrfgrid%taucldv3d(i,kl,j))
            taucldi(km)=real(wrfgrid%taucldi3d(i,kl,j))
          end do
          end do
       end do  

       i0=1
       do k=1,wrfgrid%num_wrfgrid_soil_levels
          do j= sy(1),ey(1)
          do i= sx(1),ex(1)
            km=i0+(ex(1)-sx(1)+3)*(j-sy(1)+1)+i-sx(1)+1
            soilt(km)=real(wrfgrid%soilt(i,k,j))
            soilrh(km)=real(wrfgrid%soilm(i,k,j))
          end do
          end do
        i0=i0+(nx(1)+2)*(ny(1)+2)
       end do
       
       i0=ip2mem(1)
       do j= sy(1),ey(1)
       do i= sx(1),ex(1)
           km=i0+(ex(1)-sx(1)+3)*(j-sy(1)+1)+i-sx(1)+1
           q2(km)=real(wrfgrid%q2(i,j))
           t2(km)=real(wrfgrid%t2(i,j))
           u10(km)=real(wrfgrid%u10(i,j))
           v10(km)=real(wrfgrid%v10(i,j))
           pbl_hgt(km)=real(wrfgrid%pblh(i,j))
           rmol(km)=real(wrfgrid%rmol(i,j))
           ust(km)=real(wrfgrid%ust(i,j))
           psfc(km)=real(wrfgrid%ps(i,j))
           fsnow(km)=real(wrfgrid%snowh(i,j))
           fice(km)=real(wrfgrid%xice(i,j))
           rhsfc(km)=real(wrfgrid%rh2(i,j))
           if(rhsfc(km).gt.100) then
            rhsfc(km)=100.0
           elseif(rhsfc(km).lt.0) then
            rhsfc(km)=0.0
           end if
           rainnon(km)=real(wrfgrid%rainncv(i,j)*3600*100)  ! from m/s to cm/h
           raincon(km)=real(wrfgrid%raincv(i,j)*3600*100)   ! from m/s to cm/h
           swdown(km)=real(wrfgrid%swdown(i,j))
           clflo(km)=real(wrfgrid%clflo(i,j))
           clfmi(km)=real(wrfgrid%clfmi(i,j))
           clfhi(km)=real(wrfgrid%clfhi(i,j))
        end do
       end do
          
    end subroutine wrfgrid_to_geatm
!=============================================================================================================

   subroutine geatm_to_cam_mapping(templon, templat, lonlen, latlen, xlon, xlat, mlon, mlat)
     use geatm_vartype, only : sx,ex,sy,ey,csx,cex,csy,cey

     integer,intent(in) :: latlen, lonlen, mlon, mlat
     real*8,dimension(1:latlen),intent(inout) :: templat  ! source lat
     real*8,dimension(1:lonlen),intent(inout) :: templon  ! source lon
     real*8,dimension(1:mlon,1:mlat),intent(inout) :: xlat  ! target lat
     real*8,dimension(1:mlon,1:mlat),intent(inout) :: xlon  ! target lon

     integer :: lchnk, ncols, mm, n
     integer :: i,j,k,ix,jy,km,local_index,global_index
     real*8 :: s,ss,sn,se,sw,sa
     real*8 :: dlondeg,dlatdeg
     real*8 :: pi=3.1415926_8

     real*8,dimension(:,:),allocatable:: ws,we,wn,ww,camarea
     real*8,dimension(:),allocatable:: cmipx,cmipy,cmipa

     integer :: mids, mide, mjds, mjde, mkds, mkde, &
                mims, mime, mjms, mjme, mkms, mkme, &
                mips, mipe, mjps, mjpe, mkps, mkpe
     integer ::   ids, ide, jds, jde, kds, kde, &
                  ims, ime, jms, jme, kms, kme, &
                  ips, ipe, jps, jpe, kps, kpe

       mips=sx(1)
       mipe=ex(1)
       mjps=sy(1)
       mjpe=ey(1)
       mids=1
       mide=wrfgrid%num_wrfgrid_lon
       mjds=1
       mjde=wrfgrid%num_wrfgrid_lat
       ips=csx(1)
       ipe=cex(1)
       jps=csy(1)
       jpe=cey(1)
       ids=1
       ide=camgrid%num_camgrid_lon
       jds=1
       jde=camgrid%num_camgrid_lat
     
       ! source
       allocate(ww(lonlen,latlen))
       allocate(we(lonlen,latlen))
       allocate(ws(lonlen,latlen))
       allocate(wn(lonlen,latlen))
       allocate(camarea(mlon,mlat))
       allocate(cam_local_points(mlon,mlat))

       allocate(cmipx(cam_to_geatm_points))
       allocate(cmipy(cam_to_geatm_points))
       allocate(cmipa(cam_to_geatm_points))

       dlondeg=templon(2)-templon(1)
       if (templat(2)-templat(1).gt.0.0) then
        dlatdeg=templat(2)-templat(1)
       else
        dlatdeg=templat(1)-templat(2)
       end if

       do j=1,latlen
        do i=1,lonlen                 
          if(templon(i).lt.0) templon(i)=templon(i)+360.0_8
             we(i,j)=templon(i)+dlondeg/2
             ww(i,j)=templon(i)-dlondeg/2
          if(ww(i,j).lt.0)  then  ! cross the 0 longitude
             ww(i,j)=ww(i,j)+360.0_8
             we(i,j)=we(i,j)+360.0_8
          endif
          if (templat(2)-templat(1).gt.0.0) then
          if(j.gt.1.and.j.lt.latlen) then
           wn(i,j)=(templat(j)+templat(j+1))/2
           ws(i,j)=(templat(j)+templat(j-1))/2
          end if
          if(j.eq.1) then
           wn(i,j)=(templat(j)+templat(2))/2
           ws(i,j)=-90.0_8
          end if
          if(j.eq.latlen) then
           wn(i,j)=90.0_8
           ws(i,j)=(templat(j)+templat(j-1))/2
          end if
          else
          if(j.gt.1.and.j.lt.latlen) then
           wn(i,j)=(templat(j)+templat(j-1))/2
           ws(i,j)=(templat(j)+templat(j+1))/2
          end if
          if(j.eq.1) then
           ws(i,j)=(templat(j)+templat(2))/2
           wn(i,j)=90.0_8
          end if
          if(j.eq.latlen) then
           ws(i,j)=-90.0_8
           wn(i,j)=(templat(j)+templat(j-1))/2
          end if
          endif
          if(ws(i,j).lt.-90.0_8)  ws(i,j)=-90.0_8
          if(wn(i,j).gt.90.0_8)  wn(i,j)=90.0_8
          
        end do
       end do  
      
       global_index=0
         do jy=mjps,mjpe
           do ix=mips,mipe
              ax=ix+1 
              if(ax.gt.mide) ax=1  
              bx=ix-1 
              if(bx.lt.1) bx=mide  

              ay=jy+1 
              if(ay.gt.mjde) ay=1   
              by=jy-1  
              if(by.lt.1) by=mjde

              ! the target
              if(ix.eq.1) then
                 se=(xlon(ax,jy)+xlon(ix,jy))/2+360
                 sw=(xlon(mide,jy)+360)/2
              elseif(ix.eq.mide) then
                 se=(xlon(ix,jy)+360)/2
                 sw=(xlon(bx,jy)+xlon(ix,jy))/2
              else
                 se=(xlon(ax,jy)+xlon(ix,jy))/2
                 sw=(xlon(bx,jy)+xlon(ix,jy))/2    
              endif

              if(jy.eq.1) then
                 ss=-90
                 sn=(xlat(ix,ay)+xlat(ix,jy))/2
              elseif(jy.eq.mjde) then
                 ss=(xlat(ix,by)+xlat(ix,jy))/2
                 sn=90                    
              else
                 ss=(xlat(ix,by)+xlat(ix,jy))/2
                 sn=(xlat(ix,ay)+xlat(ix,jy))/2     
              endif
              s=sin(sn*pi/180)-sin(ss*pi/180)  
              camarea(ix,jy)=s*(se-sw)

           ! source to target
           local_index=0
           sa=0
           do j=1,latlen
           do i=1,lonlen

            ! case 1 -case 9 and case 10-1,2,3 and case 10-5,6,7
            ! are the templat and templon finer than xlat and xlon

            ! case 1
            if(wn(i,j).le.sn.and.ws(i,j).ge.ss.and.ww(i,j).ge.sw.and.we(i,j).le.se) then
              local_index=local_index+1
              cmipx(local_index)=i
              cmipy(local_index)=j
              cmipa(local_index)=(sin(wn(i,j)*pi/180)-sin(ws(i,j)*pi/180))*(we(i,j)-ww(i,j))       
              sa=sa+cmipa(local_index)
              goto 999
            endif

            ! case 2
            if(wn(i,j).le.sn.and.ws(i,j).ge.ss.and.ww(i,j).le.sw.and.&
               we(i,j).le.se.and.we(i,j).ge.sw) then
               local_index=local_index+1
               cmipx(local_index)=i
               cmipy(local_index)=j
               cmipa(local_index)=(sin(wn(i,j)*pi/180)-sin(ws(i,j)*pi/180))*(we(i,j)-sw)
               sa=sa+cmipa(local_index)
               goto 999
            endif

            ! case 3
            if(wn(i,j).le.sn.and.ws(i,j).ge.ss.and.ww(i,j).ge.sw.and.&
               ww(i,j).le.se.and.we(i,j).ge.se) then
               local_index=local_index+1
               cmipx(local_index)=i
               cmipy(local_index)=j
               cmipa(local_index)=(sin(wn(i,j)*pi/180)-sin(ws(i,j)*pi/180))*(se-ww(i,j))
               sa=sa+cmipa(local_index)
               goto 999 
            endif

            ! case4
            if(wn(i,j).ge.sn.and.ws(i,j).ge.ss.and.ws(i,j).le.sn.and.&
               ww(i,j).ge.sw.and.we(i,j).le.se) then 
               local_index=local_index+1 
               cmipx(local_index)=i 
               cmipy(local_index)=j 
               cmipa(local_index)=(sin(sn*pi/180)-sin(ws(i,j)*pi/180))*(we(i,j)-ww(i,j))
               sa=sa+cmipa(local_index)
               goto 999
            endif

            ! case5
            if(wn(i,j).le.sn.and.wn(i,j).ge.ss.and.ws(i,j).le.ss.and.&
               ww(i,j).ge.sw.and.we(i,j).le.se) then
               local_index=local_index+1
               cmipx(local_index)=i
               cmipy(local_index)=j
               cmipa(local_index)=(sin(wn(i,j)*pi/180)-sin(ss*pi/180))*(we(i,j)-ww(i,j))
               sa=sa+cmipa(local_index)
               goto 999
            endif 

            ! case6
            if(wn(i,j).ge.sn.and.ws(i,j).ge.ss.and.ws(i,j).le.sn.and.&
               we(i,j).ge.se.and.ww(i,j).ge.sw.and.ww(i,j).le.se) then
               local_index=local_index+1
               cmipx(local_index)=i
               cmipy(local_index)=j
               cmipa(local_index)=(sin(sn*pi/180)-sin(ws(i,j)*pi/180))*(se-ww(i,j))
               sa=sa+cmipa(local_index)
               goto 999
            endif

            ! case 7
            if(wn(i,j).ge.sn.and.ws(i,j).ge.ss.and.ws(i,j).le.sn.and.&
               we(i,j).le.se.and.we(i,j).ge.sw.and.ww(i,j).le.sw) then
               local_index=local_index+1
               cmipx(local_index)=i
               cmipy(local_index)=j
               cmipa(local_index)=(sin(sn*pi/180)-sin(ws(i,j)*pi/180))*(we(i,j)-sw)
               sa=sa+cmipa(local_index)
               goto 999
            endif

            ! case 8
            if(wn(i,j).le.sn.and.wn(i,j).ge.ss.and.ws(i,j).le.ss.and.&
               we(i,j).le.se.and.we(i,j).ge.sw.and.ww(i,j).le.sw) then
               local_index=local_index+1 
               cmipx(local_index)=i
               cmipy(local_index)=j 
               cmipa(local_index)=(sin(wn(i,j)*pi/180)-sin(ss*pi/180))*(we(i,j)-sw)
               sa=sa+cmipa(local_index)
               goto 999
            endif

            ! case 9
            if(wn(i,j).le.sn.and.wn(i,j).ge.ss.and.ws(i,j).le.ss.and.&
               we(i,j).ge.se.and.ww(i,j).ge.sw.and.ww(i,j).le.se) then
               local_index=local_index+1  
               cmipx(local_index)=i 
               cmipy(local_index)=j  
               cmipa(local_index)=(sin(wn(i,j)*pi/180)-sin(ss*pi/180))*(se-ww(i,j))
               sa=sa+cmipa(local_index)
               goto 999  
            endif
                   
            ! case 10,  templon doesn't across the 0 degree but cam's grid across the 0 longitude
            if((ww(i,j)+360-sw).ge.0.and.(ww(i,j)+360-se).le.0) then

             if((we(i,j)+360-se).le.0) then
              ! case 10-1
              if(wn(i,j).le.sn.and.wn(i,j).ge.ss.and.ws(i,j).le.ss) then
                 local_index=local_index+1
                 cmipx(local_index)=i
                 cmipy(local_index)=j
                 cmipa(local_index)=(sin(wn(i,j)*pi/180)-sin(ss*pi/180))*(we(i,j)-ww(i,j))
                 sa=sa+cmipa(local_index)
                 goto 999
              endif
              ! case 10-2
              if(wn(i,j).le.sn.and.ws(i,j).gt.ss) then
                 local_index=local_index+1
                 cmipx(local_index)=i 
                 cmipy(local_index)=j 
                 cmipa(local_index)=(sin(wn(i,j)*pi/180)-sin(ws(i,j)*pi/180))*(we(i,j)-ww(i,j))
                 sa=sa+cmipa(local_index)
                 goto 999
              endif
              ! case 10-3
              if(wn(i,j).gt.sn.and.ws(i,j).gt.ss.and.ws(i,j).le.sn) then
                 local_index=local_index+1
                 cmipx(local_index)=i
                 cmipy(local_index)=j
                 cmipa(local_index)=(sin(sn*pi/180)-sin(ws(i,j)*pi/180))*(we(i,j)-ww(i,j))
                 sa=sa+cmipa(local_index)
                 goto 999
              endif
              ! case 10-4
              if(wn(i,j).gt.sn.and.ws(i,j).lt.ss) then
                 local_index=local_index+1
                 cmipx(local_index)=i
                 cmipy(local_index)=j
                 cmipa(local_index)=(sin(sn*pi/180)-sin(ss*pi/180))*(we(i,j)-ww(i,j))
                 sa=sa+cmipa(local_index)
                 goto 999
              endif
             end if

             if((we(i,j)+360-se).gt.0) then
              ! case 10-5
              if(wn(i,j).le.sn.and.wn(i,j).ge.ss.and.ws(i,j).le.ss) then
                 local_index=local_index+1
                 cmipx(local_index)=i
                 cmipy(local_index)=j
                 cmipa(local_index)=(sin(wn(i,j)*pi/180)-sin(ss*pi/180))*(se-ww(i,j)-360)
                 sa=sa+cmipa(local_index)
                 goto 999
              endif
              ! case 10-6
              if(wn(i,j).le.sn.and.ws(i,j).gt.ss) then
                 local_index=local_index+1
                 cmipx(local_index)=i
                 cmipy(local_index)=j
                 cmipa(local_index)=(sin(wn(i,j)*pi/180)-sin(ws(i,j)*pi/180))*(se-ww(i,j)-360)
                 sa=sa+cmipa(local_index)
                 goto 999
              endif
              ! case 10-7
              if(wn(i,j).ge.sn.and.ws(i,j).gt.ss.and.ws(i,j).le.sn) then
                 local_index=local_index+1
                 cmipx(local_index)=i
                 cmipy(local_index)=j
                 cmipa(local_index)=(sin(sn*pi/180)-sin(ws(i,j)*pi/180))*(se-ww(i,j)-360)
                 sa=sa+cmipa(local_index)
                 goto 999
              endif
              ! case 10-8
              if(wn(i,j).ge.sn.and.ws(i,j).le.ss) then
                 local_index=local_index+1
                 cmipx(local_index)=i
                 cmipy(local_index)=j
                 cmipa(local_index)=(sin(sn*pi/180)-sin(ss*pi/180))*(se-ww(i,j)-360)
                 sa=sa+cmipa(local_index)
                 goto 999
              ENDIF
             end if

            end if

             ! case 11
             if(wn(i,j).ge.sn.and.ws(i,j).le.ss.and.ww(i,j).le.sw.and.we(i,j).ge.se) then
                 local_index=local_index+1
                 cmipx(local_index)=i
                 cmipy(local_index)=j
                 cmipa(local_index)=(sin(sn*pi/180)-sin(ss*pi/180))*(se-sw)
                 sa=sa+cmipa(local_index)
                 goto 999
             ENDIF

             ! case 12
             if(wn(i,j).ge.sn.and.ws(i,j).le.ss.and.ww(i,j).le.se.and.ww(i,j).ge.sw.and.&
                we(i,j).ge.se) then
                 local_index=local_index+1
                 cmipx(local_index)=i
                 cmipy(local_index)=j
                 cmipa(local_index)=(sin(sn*pi/180)-sin(ss*pi/180))*(se-ww(i,j))
                 sa=sa+cmipa(local_index)
                 goto 999
             ENDIF

             ! case 13
             if(wn(i,j).ge.sn.and.ws(i,j).le.ss.and.ww(i,j).le.sw.and.we(i,j).ge.sw.and.&
                we(i,j).le.se) then
                 local_index=local_index+1
                 cmipx(local_index)=i
                 cmipy(local_index)=j
                 cmipa(local_index)=(sin(sn*pi/180)-sin(ss*pi/180))*(we(i,j)-sw)
                 sa=sa+cmipa(local_index)
                 goto 999
             ENDIF

             ! case 14
             if(wn(i,j).gt.sn.and.ws(i,j).ge.ss.and.ws(i,j).le.sn.and. &
                ww(i,j).le.sw.and.we(i,j).ge.se) then
                 local_index=local_index+1
                 cmipx(local_index)=i
                 cmipy(local_index)=j
                 cmipa(local_index)=(sin(sn*pi/180)-sin(ws(i,j)*pi/180))*(se-sw)
                 sa=sa+cmipa(local_index)
                 goto 999
             ENDIF

             ! case 15
             if(wn(i,j).le.sn.and.wn(i,j).ge.ss.and.ws(i,j).le.ss.and. &
                ww(i,j).le.sw.and.we(i,j).ge.se) then
                 local_index=local_index+1
                 cmipx(local_index)=i
                 cmipy(local_index)=j
                 cmipa(local_index)=(sin(wn(i,j)*pi/180)-sin(ss*pi/180))*(se-sw)
                 sa=sa+cmipa(local_index)
                 goto 999
             ENDIF

             ! case 16
             if(wn(i,j).ge.sn.and.ws(i,j).le.ss.and. &
                ww(i,j).ge.sw.and.we(i,j).le.se) then
                 local_index=local_index+1
                 cmipx(local_index)=i
                 cmipy(local_index)=j
                 cmipa(local_index)=(sin(sn*pi/180)-sin(ss*pi/180))*(we(i,j)-ww(i,j))
                 sa=sa+cmipa(local_index)
                 goto 999
             ENDIF

             ! case 17
             if(wn(i,j).le.sn.and.ws(i,j).ge.ss.and. &
                ww(i,j).le.sw.and.we(i,j).ge.se) then
                 local_index=local_index+1
                 cmipx(local_index)=i
                 cmipy(local_index)=j
                 cmipa(local_index)=(sin(wn(i,j)*pi/180)-sin(ws(i,j)*pi/180))*(se-sw)
                 sa=sa+cmipa(local_index)
                 goto 999
             ENDIF

             ! case 18, the original grid across 0, the target doesn't across 0
             if(we(i,j).gt.360) then
                if(we(i,j)-360.ge.sw.and.(we(i,j)-360).le.se) then               
                 ! case 18-1
                 if(wn(i,j).le.sn.and.wn(i,j).ge.ss.and.ws(i,j).le.ss) then                 
                 local_index=local_index+1
                 cmipx(local_index)=i
                 cmipy(local_index)=j
                 cmipa(local_index)=(sin(wn(i,j)*pi/180)-sin(ss*pi/180))*(we(i,j)-360-sw)
                 sa=sa+cmipa(local_index)
                 goto 999
                 end if
                 ! case 18-2
                 if(wn(i,j).ge.sn.and.ws(i,j).le.ss) then
                 local_index=local_index+1
                 cmipx(local_index)=i
                 cmipy(local_index)=j
                 cmipa(local_index)=(sin(sn*pi/180)-sin(ss*pi/180))*(we(i,j)-360-sw)
                 sa=sa+cmipa(local_index)
                 goto 999
                 end if
                 ! case 18-3
                 if(wn(i,j).ge.sn.and.ws(i,j).ge.ss.and.ws(i,j).le.sn) then
                 local_index=local_index+1
                 cmipx(local_index)=i
                 cmipy(local_index)=j
                 cmipa(local_index)=(sin(sn*pi/180)-sin(ws(i,j)*pi/180))*(we(i,j)-360-sw)
                 sa=sa+cmipa(local_index)
                 goto 999
                 end if
                endif 
                if(we(i,j)-360.ge.se) then
                 ! case 18-4
                 if(wn(i,j).le.sn.and.wn(i,j).ge.ss.and.ws(i,j).le.ss) then
                 local_index=local_index+1
                 cmipx(local_index)=i
                 cmipy(local_index)=j
                 cmipa(local_index)=(sin(wn(i,j)*pi/180)-sin(ss*pi/180))*(se-sw)
                 sa=sa+cmipa(local_index)
                 goto 999
                 end if
                 ! case 18-5
                 if(wn(i,j).ge.sn.and.ws(i,j).le.ss) then
                 local_index=local_index+1
                 cmipx(local_index)=i
                 cmipy(local_index)=j
                 cmipa(local_index)=(sin(sn*pi/180)-sin(ss*pi/180))*(se-sw)
                 sa=sa+cmipa(local_index)
                 goto 999
                 end if
                 ! case 18-6
                 if(wn(i,j).ge.sn.and.ws(i,j).ge.ss.and.ws(i,j).le.sn) then
                 local_index=local_index+1
                 cmipx(local_index)=i
                 cmipy(local_index)=j
                 cmipa(local_index)=(sin(sn*pi/180)-sin(ws(i,j)*pi/180))*(se-sw)
                 sa=sa+cmipa(local_index)
                 goto 999
                 end if
                endif
             ENDIF

999               continue
                enddo
              enddo

              s=abs(sa-camarea(ix,jy))/camarea(ix,jy)
              if(s.gt.0.0001_8) then
               print *, s,local_index
               print *, ss,sn,sw,se,camarea(ix,jy)
              end if

              sa=0.0
              cam_local_points(ix,jy)=local_index
              do i=1,local_index
              if(s.gt.0.0001_8) then
               print *, ws(cmipx(i),cmipy(i)),wn(cmipx(i),cmipy(i)),ww(cmipx(i),cmipy(i)),we(cmipx(i),cmipy(i)),cmipa(i)
              end if 
              global_index=global_index+1
              remap_cami(global_index)=cmipx(i)
              remap_camj(global_index)=cmipy(i)
              cam_weight(global_index)=cmipa(i)/camarea(ix,jy)
              sa=sa+cam_weight(global_index)
              enddo
        enddo
        enddo

         deallocate(ws)
         deallocate(wn)
         deallocate(we)
         deallocate(ww)
         deallocate(cmipx)
         deallocate(cmipy)
         deallocate(cmipa)
         deallocate(camarea)

  end subroutine geatm_to_cam_mapping
!===============================================================================
   subroutine garther_camgrid_once(grid, ids, ide, jds, jde, kds, kde, &
                              num_metgrid_levels, num_metgrid_soil_levels, &
                              ips, ipe, jps, jpe, kps, kpe)
   use geatm_vartype, only : local_com
   include "mpif.h"

   type(camgrid_c)::grid
   integer,intent(in)::ids, ide, jds, jde, kds, kde,  &   
                       ips, ipe, jps, jpe, kps, kpe,  &
                       num_metgrid_levels, num_metgrid_soil_levels 
                               
   real(8), dimension(:,:,:), allocatable :: patch_array
   real(8), dimension(:,:,:), allocatable :: domain_array
   integer :: i, j, k, mpi_ierr
   
     allocate(patch_array(ips:ipe,1:num_metgrid_soil_levels,jps:jpe))
     allocate(domain_array(ids:ide,1:num_metgrid_soil_levels,jds:jde))

     ! soilthick
      do j=jps, jpe
       do k=1,num_metgrid_soil_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%soilthick(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, num_metgrid_soil_levels, jps, jpe, &
                              domain_array, ids, ide, 1, num_metgrid_soil_levels, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_soil_levels, MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      grid%soilthick=domain_array      

      ! soildepth
      do j=jps, jpe
       do k=1,num_metgrid_soil_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%soildepth(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, num_metgrid_soil_levels, jps, jpe, & 
                              domain_array, ids, ide, 1, num_metgrid_soil_levels, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_soil_levels, MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      grid%soildepth=domain_array
      
      deallocate(patch_array)
      deallocate(domain_array)
      
      allocate(patch_array(ips:ipe,1:1,jps:jpe))
      allocate(domain_array(ids:ide,1:1,jds:jde))
      
      ! xlat
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%xlat(i,j)
        enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      do j=jds, jde
        do i=ids, ide   
         grid%xlat(i,j)=domain_array(i,1,j)
        enddo
      enddo
      
      ! xlon
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%xlon(i,j)
        enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      do j=jds, jde
        do i=ids, ide                           
         grid%xlon(i,j)=domain_array(i,1,j)
        enddo
      enddo

      deallocate(patch_array)
      deallocate(domain_array)

   end subroutine garther_camgrid_once
!===============================================================================

   subroutine garther_camgrid(grid, ids, ide, jds, jde, kds, kde, &
                              num_metgrid_levels, num_metgrid_soil_levels, &
                              ips, ipe, jps, jpe, kps, kpe)
   use geatm_vartype, only : local_com
   include "mpif.h"

   type(camgrid_c)::grid
   integer,intent(in)::ids, ide, jds, jde, kds, kde,  &   
                       ips, ipe, jps, jpe, kps, kpe,  &
                       num_metgrid_levels, num_metgrid_soil_levels 
                               
   real(8), dimension(:,:,:), allocatable :: patch_array
   real(8), dimension(:,:,:), allocatable :: domain_array
   integer :: i, j, k, mpi_ierr
   
      allocate(patch_array(ips:ipe,1:num_metgrid_levels,jps:jpe))
      allocate(domain_array(ids:ide,1:num_metgrid_levels,jds:jde))

      ! p3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%p3d(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field(patch_array, ips, ipe, 1,num_metgrid_levels,jps, jpe, &
                              domain_array, ids, ide, 1,num_metgrid_levels,jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_levels, MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      grid%p3d=domain_array
      call mpi_barrier(local_com, mpi_ierr)

      ! u3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%u3d(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field(patch_array, ips, ipe, 1, num_metgrid_levels,jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels,jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_levels,MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      grid%u3d=domain_array
      call mpi_barrier(local_com, mpi_ierr)

      ! v3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%v3d(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field(patch_array, ips, ipe, 1, num_metgrid_levels,jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels,jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_levels,MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      grid%v3d=domain_array
      call mpi_barrier(local_com, mpi_ierr)

      ! t3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%t3d(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field(patch_array, ips, ipe, 1, num_metgrid_levels,jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels,jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_levels,MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      grid%t3d=domain_array
      call mpi_barrier(local_com, mpi_ierr)

      ! qv3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%qv3d(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field(patch_array, ips, ipe, 1, num_metgrid_levels,jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels,jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_levels,MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      grid%qv3d=domain_array
      call mpi_barrier(local_com, mpi_ierr)
      
      ! qc3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%qc3d(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field(patch_array, ips, ipe, 1, num_metgrid_levels,jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels,jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_levels,MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      grid%qc3d=domain_array
      call mpi_barrier(local_com, mpi_ierr)

      ! qi3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%qi3d(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field(patch_array, ips, ipe, 1, num_metgrid_levels,jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels,jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_levels,MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      grid%qi3d=domain_array
      call mpi_barrier(local_com, mpi_ierr)      

      ! rh3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%rh3d(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field(patch_array, ips, ipe, 1,num_metgrid_levels,jps, jpe, &
                                 domain_array, ids, ide, 1,num_metgrid_levels,jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_levels,MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      grid%rh3d=domain_array
      call mpi_barrier(local_com, mpi_ierr)

      ! z3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%z3d(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field(patch_array, ips, ipe, 1, num_metgrid_levels,jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels,jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_levels,MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      grid%z3d=domain_array
      call mpi_barrier(local_com, mpi_ierr)

      ! taucldv3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%taucldv3d(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field(patch_array, ips, ipe, 1, num_metgrid_levels,jps, jpe, &
                              domain_array, ids, ide, 1, num_metgrid_levels,jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_levels,MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      grid%taucldv3d=domain_array
      call mpi_barrier(local_com, mpi_ierr)      

      ! taucldi3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%taucldi3d(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field(patch_array, ips, ipe, 1, num_metgrid_levels,jps, jpe, &
                              domain_array, ids, ide, 1, num_metgrid_levels,jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_levels,MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      grid%taucldi3d=domain_array
      call mpi_barrier(local_com, mpi_ierr) 
      
      deallocate(patch_array)
      deallocate(domain_array)

     allocate(patch_array(ips:ipe,1:num_metgrid_soil_levels,jps:jpe))
     allocate(domain_array(ids:ide,1:num_metgrid_soil_levels,jds:jde))

      ! soilt
      do j=jps, jpe
       do k=1,num_metgrid_soil_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%soilt(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, num_metgrid_soil_levels, jps, jpe, &
                              domain_array, ids, ide, 1, num_metgrid_soil_levels, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_soil_levels, MPI_DOUBLE_PRECISION, 0, local_com, ierror)                                    
      grid%soilt=domain_array
      
      ! soilm
      do j=jps, jpe
       do k=1,num_metgrid_soil_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%soilm(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, num_metgrid_soil_levels, jps, jpe, &
                              domain_array, ids, ide, 1, num_metgrid_soil_levels, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_soil_levels, MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      grid%soilm=domain_array      
      
      deallocate(patch_array)
      deallocate(domain_array)
      
      allocate(patch_array(ips:ipe,1:1,jps:jpe))
      allocate(domain_array(ids:ide,1:1,jds:jde))
      
      ! snow
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%snowh(i,j)
        enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      do j=jds, jde
        do i=ids, ide                           
         grid%snowh(i,j)=domain_array(i,1,j)
        enddo
      enddo
                
     if(1.eq.0) then
     ! no use of xland, tsk, sst and ht during the integration
     ! xland
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%xland(i,j)
        enddo
      enddo
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      do j=jds, jde
        do i=ids, ide
         grid%xland(i,j)=domain_array(i,1,j)
        enddo
      enddo

      ! ht
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%ht(i,j)
        enddo
      enddo
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      do j=jds, jde
        do i=ids, ide
         grid%ht(i,j)=domain_array(i,1,j)
        enddo
      enddo

      ! tsk
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%tsk(i,j)
        enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0,local_com, ierror)
      do j=jds, jde
        do i=ids, ide
         grid%tsk(i,j)=domain_array(i,1,j)
        enddo
      enddo

      ! sst
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%sst(i,j)
        enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      do j=jds, jde
        do i=ids, ide                           
         grid%sst(i,j)=domain_array(i,1,j)
        enddo
      enddo
      
      end if

      ! xice
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%xice(i,j)
        enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      do j=jds, jde
        do i=ids, ide                           
         grid%xice(i,j)=domain_array(i,1,j)
        enddo
      enddo
      
      ! psfc
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%ps(i,j)
        enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      do j=jds, jde
        do i=ids, ide                           
         grid%ps(i,j)=domain_array(i,1,j)
        enddo
      enddo
      
      ! raincv
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%raincv(i,j)
        enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      do j=jds, jde
        do i=ids, ide                         
         grid%raincv(i,j)=domain_array(i,1,j)
        enddo
      enddo   
      
      ! rainncv
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%rainncv(i,j)
        enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      do j=jds, jde
        do i=ids, ide                         
         grid%rainncv(i,j)=domain_array(i,1,j)
        enddo
      enddo
    
      ! u10
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%u10(i,j)
        enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      do j=jds, jde
        do i=ids, ide                         
         grid%u10(i,j)=domain_array(i,1,j)
        enddo
      enddo
  
      ! v10
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%v10(i,j)
        enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      do j=jds, jde
        do i=ids, ide                         
         grid%v10(i,j)=domain_array(i,1,j)
        enddo
      enddo      
 
      ! t2
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%t2(i,j)
        enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      do j=jds, jde
        do i=ids, ide                         
         grid%t2(i,j)=domain_array(i,1,j)
        enddo
      enddo

      ! q2
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%q2(i,j)
        enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      do j=jds, jde
        do i=ids, ide                         
         grid%q2(i,j)=domain_array(i,1,j)
        enddo
      enddo
 
      ! rh2
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%rh2(i,j)
        enddo
      enddo
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0,local_com, ierror)
      do j=jds, jde
        do i=ids, ide
         grid%rh2(i,j)=domain_array(i,1,j)
        enddo
      enddo

      ! pblh
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%pblh(i,j)
        enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      do j=jds, jde
        do i=ids, ide                         
         grid%pblh(i,j)=domain_array(i,1,j)
        enddo
      enddo
   
      ! ust
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%ust(i,j)
        enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      do j=jds, jde
        do i=ids, ide                         
         grid%ust(i,j)=domain_array(i,1,j)
        enddo
      enddo

      ! rmol
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%rmol(i,j)
        enddo
      enddo
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0,local_com, ierror)
      do j=jds, jde
        do i=ids, ide
         grid%rmol(i,j)=domain_array(i,1,j)
        enddo
      enddo

      ! swdown
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%swdown(i,j)
        enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      do j=jds, jde
        do i=ids, ide                         
         grid%swdown(i,j)=domain_array(i,1,j)
        enddo
      enddo

      ! clflo
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%clflo(i,j)
        enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      do j=jds, jde
        do i=ids, ide                         
         grid%clflo(i,j)=domain_array(i,1,j)
        enddo
      enddo      
 
      ! clfmi
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%clfmi(i,j)
        enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      do j=jds, jde
        do i=ids, ide                         
         grid%clfmi(i,j)=domain_array(i,1,j)
        enddo
      enddo  
 
      ! clfhi
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%clfhi(i,j)
        enddo
      enddo 
      call gather_whole_field(patch_array, ips, ipe, 1, 1, jps, jpe, &
                              domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_DOUBLE_PRECISION, 0, local_com, ierror)
      do j=jds, jde
        do i=ids, ide                         
         grid%clfhi(i,j)=domain_array(i,1,j)
        enddo
      enddo  

      deallocate(patch_array)
      deallocate(domain_array)

   end subroutine garther_camgrid
!===============================================================================

      Subroutine interpolate_to_pressure_levels(kps,kpe &
                   ,nx,ny,nz &
                   ,P,P_int,field_data,field_data_int)

       implicit none
       integer,intent(in) :: kps,kpe,nx,ny,nz
       real*8, dimension(nx,kps:kpe,ny), intent(in) :: field_data
       real*8, dimension(nx,kps:kpe,ny), intent(in) :: p
       real*8, dimension(nx,nz,ny), intent(in) :: P_int
       real*8, dimension(nx,nz,ny), intent(out) :: field_data_int          
        
       real*8 :: X(kpe-kps+1),Y(kpe-kps+1),Y2(kpe-kps+1),XINT,YINT,yp1,ypn
       integer :: i,j,k

       yp1=2e30
       ypn=2e30
       ! interpolate from sigma to pressure coordinates:
       do j=1,ny
       do i=1,nx 
        do k=kps,kpe
         X(k-kps+1)=P(i,k,j)*1.0_8
         Y(k-kps+1)=field_data(i,k,j)*1.0_8
        enddo
        call spline(X,Y,yp1,ypn,kpe-kps+1,Y2)
        do k=1,nz
         XINT=P_int(i,k,j)*1.0_8
         call splint(X,Y,Y2,kpe-kps+1,XINT,YINT)
         field_data_int(i,k,j)=YINT
        end do
       end do
       end do
       return
      end subroutine interpolate_to_pressure_levels

!===============================================================================

       SUBROUTINE spline(x,y,yp1,ypn,n,y2) 
             implicit none
             INTEGER :: n
             INTEGER, parameter :: NMAX=500
             REAL*8 :: x(n),y(n),y2(n) 
             INTEGER :: i,k
             REAL*8 :: p,qn,sig,un,u(NMAX),yp1,ypn
  
           if(yp1.gt..99e30) then
             y2(1)=0
             u(1)=0
           else
             y2(1)=-0.5
             u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
           endif
  
             do i=2,n-1
             sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
             p=sig*y2(i-1)+2.
             y2(i)=(sig-1.)/p
             u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
             /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
             end do

           if(ypn.gt..99e30) then
             qn=0
             un=0
           else              
             qn=0.5
             un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
            endif
  
             y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
             
             do k=n-1,1,-1
             y2(k)=y2(k)*y2(k+1)+u(k)
             end do
             
             return
             
       end SUBROUTINE spline

!===============================================================================
        
       SUBROUTINE splint(xa,ya,y2a,n,x,y)
         implicit none
         INTEGER :: n
         REAL*8 :: x,y,xa(n),y2a(n),ya(n)
         INTEGER :: k,khi,klo
         REAL*8 :: a,b,h

         ! Eli: avoid actual extrapolation by using the end values:
         if (x<xa(n)) then
                 y=ya(n)
         elseif (x>xa(1)) then
                 y=ya(1)
         else
         ! Eli: end of my addition here.
                 klo=1
                 khi=n
   1    if (khi-klo.gt.1) then
             k=(khi+klo)/2
         if(xa(k).lt.x)then
                 khi=k
         else
                 klo=k
         end if
         goto 1
        end if
        h=xa(khi)-xa(klo)
        if (h.eq.0.) pause 'bad xa input in splint'
        a=(xa(khi)-x)/h
        b=(x-xa(klo))/h
        y=a*ya(klo)+b*ya(khi)+ &
          ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
        end if
        return
       end SUBROUTINE splint

end module geatm_comp_mct
