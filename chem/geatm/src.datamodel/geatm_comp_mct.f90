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
 
  public :: geatm_init_mct
  public :: geatm_run_mct
  public :: geatm_final_mct

  private :: geatm_SetgsMap_cam
  private :: geatm_domain_cam
 
  integer, parameter  :: nlen = 256     ! Length of character strings

!
! Time averaged counter for flux fields
!

  type camgrid_c
  integer :: num_camgrid_soil_levels, num_camgrid_levels, num_camgrid_lat, num_camgrid_lon
  end type camgrid_c
  type(camgrid_c):: camgrid

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
       
!--------------------------------------------------------------------
!   parameter reading
!--------------------------------------------------------------------       
    
    OPEN(31,FILE='input.dat',form='formatted')
      
    read(31,*)iyear1,imonth1,idate1,ihour1
    read(31,*) camgrid%num_camgrid_lon, camgrid%num_camgrid_lat, &
               camgrid%num_camgrid_levels, camgrid%num_camgrid_soil_levels
    read(31,*) ne
    close(31)
    
!--------------------------------------------------------------------
!    parallel
!--------------------------------------------------------------------   
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

  call mpi_cart_get( comm2d(ne), 2, dims(1,ne), periods(1,ne),  &
                     coords(1,ne), ierr )
  cnx(ne)=camgrid%num_camgrid_lon
  cny(ne)=camgrid%num_camgrid_lat
  call mpe_decomp1d( cnx(ne), dims(1,ne), coords(1,ne), csx(ne), cex(ne) )
  call mpe_decomp1d( cny(ne), dims(2,ne), coords(2,ne), csy(ne), cey(ne) )

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
      
        first_time=.false.    
        print *,camgrid%num_camgrid_lon, camgrid%num_camgrid_lat, &
               camgrid%num_camgrid_levels, camgrid%num_camgrid_soil_levels
else
end if ! first time

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

end subroutine geatm_init_mct

!================================================================================
	
subroutine geatm_run_mct(EClock_aa, EClock, cdata_a, x2c_c1, x2c_c2, c2x_c)
    ! 
    ! Arguments
    !
    type(ESMF_Clock)            ,intent(in)    :: EClock_aa ! cam
    type(ESMF_Clock)            ,intent(in)    :: EClock  ! geatm
    type(seq_cdata)             ,intent(inout) :: cdata_a ! geatm
    type(mct_aVect)             ,intent(inout) :: x2c_c1   ! cam -> cpl -> geatm  
    type(mct_aVect)             ,intent(inout) :: x2c_c2   ! cam -> cpl -> geatm
    type(mct_aVect)             ,intent(inout) :: c2x_c   ! geatm -> cpl -> cam
    
        
end subroutine geatm_run_mct
!================================================================================
subroutine geatm_final_mct()
   use geatm_vartype
   if(allocated(procs)) deallocate (procs)

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
	
end module geatm_comp_mct
