module cam_nudging
!----------------------------------------------------------------------- 
! module to handle reading of the master nudging file.
!----------------------------------------------------------------------- 
   use shr_kind_mod, only: r8 => shr_kind_r8
   use spmd_utils,   only: masterproc
   use ppgrid          
   use physics_types,only: physics_state, physics_tend
   use camsrfexch_types,only: cam_in_t
   use abortutils,   only: endrun
   use cam_logfile,  only: iulog
   use cam_pio_utils,only: pio_subsystem, fillvalue, cam_pio_openfile
   use pio,          only: file_desc_t, var_desc_t, io_desc_t, &
                           iotype_pnetcdf, iotype_netcdf, &
                           pio_noerr, pio_int, pio_real, pio_double, pio_char, &
                           pio_inq_varid, pio_read_darray, pio_openfile, pio_closefile, &
                           pio_initdecomp, pio_freedecomp
   use time_manager,   only: get_nstep,get_step_size,is_last_step, &
                             get_curr_date, get_curr_calday
   use nudging_driver, only: fdda, fdda_sfc, analysis_interval, guv, gt ,gq, &  
                             if_zfac_uv, if_zfac_t, if_zfac_q, & 
                             k_zfac_uv, k_zfac_t, k_zfac_q,& 
                             k_zbot_uv, k_zbot_t, k_zbot_q,& 
                             vert_int,data_src, data_sfcsrc
   implicit none

   ! CMIP
   integer :: latlen, lonlen, levlen,latlensfc,lonlensfc
   integer,parameter :: cmip_to_cam_points = 1000000
   real*8, dimension(1:cmip_to_cam_points) :: cmip_weight,cmip_weightsfc 
   integer, dimension(1:cmip_to_cam_points) :: cmip_local_points, cmip_local_pointssfc, & 
                                               remap_cmipi, remap_cmipj, remap_cmipisfc, remap_cmipjsfc
   real :: fill_value

   ! nudging data stored
   real*8,dimension(:,:),allocatable :: div_ndg_old, vort_ndg_old, u_ndg_old, &
                                       v_ndg_old, t_ndg_old, q_ndg_old, &
                                       div_ndg_new, vort_ndg_new, u_ndg_new, &
                                       v_ndg_new, t_ndg_new, q_ndg_new, &
                                       ps_ndg_old, ps_ndg_new, p_ndg_old, p_ndg_new,&
                                       ts_ndg_old,ts_ndg_new,wsx_ndg_old,wsx_ndg_new,&
                                       wsy_ndg_old,wsy_ndg_new,lhf_ndg_old,lhf_ndg_new,&
                                       shf_ndg_old,shf_ndg_new,tref_ndg_old,tref_ndg_new,&
                                       qref_ndg_old,qref_ndg_new,u10_ndg_old,u10_ndg_new,&
                                       ulwrf_ndg_old,ulwrf_ndg_new
   real*8,dimension(:),allocatable :: cmiplev,hyam,hybm

   ! FNL
   integer,parameter :: ndglev =26
   real*8,dimension(1:ndglev) :: fnllev
   data fnllev/1000.0D0,2000.0D0,3000.0D0,5000.0D0,7000.0D0,10000.0D0,15000.0D0,&
               20000.0D0,25000.0D0,30000.0D0,35000.0D0,40000.0D0,45000.0D0,50000.0D0,&
               55000.0D0,60000.0D0,65000.0D0,70000.0D0,75000.0D0,80000.0D0,85000.0D0,&
               90000.0D0,92500.0D0,95000.0D0,97500.0D0,100000.0D0/

!=========================================================================================
CONTAINS
!=========================================================================================
 subroutine cam_read_gfs( phys_state, nstep )

!----------------------------------------------------------------------- 
! Purpose: 
! Acquire and position the restart, master, primary and secondary
! datasets for a continuation run
! Author:Juanxiong He, 2012-06-01
!-----------------------------------------------------------------------
    use phys_grid, only: get_ncols_p, get_gcol_p, ngcols
    use pmgrid, only: plon, plev, plat
    implicit none

    type(physics_state), pointer :: phys_state(:)
    integer, intent(in) :: nstep

    integer :: lchnk, i, j, k, ii, m, n, csize, ncols, mm
    integer :: ntime, nnstep, fdda_nhtfrq, days, seconds
    real(r8) :: dt
    real, dimension(:,:,:), allocatable :: tmpfield1_3d,tmpfield2_3d
    real, dimension(:,:), allocatable :: xlat,xlon
    real*8, dimension(:), allocatable :: camlev
    real*8, dimension(:), allocatable :: fdatai,fdatao
    real*8, dimension(:), allocatable :: tmplat,tmplon

    fill_value=1.0e30
    dt = get_step_size()
    fdda_nhtfrq=int(analysis_interval/dt)
    ii=mod(nstep,fdda_nhtfrq )

    nnstep=nstep*dt/analysis_interval

    mm=0
    do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols
     n=get_gcol_p(lchnk,i)
     mm=mm+1
    end do
    end do

    if(.not.allocated(u_ndg_old)) then
    allocate(u_ndg_old(1:mm,1:ndglev))
    if(.not.allocated(v_ndg_old)) allocate(v_ndg_old(1:mm,1:ndglev))
    if(.not.allocated(t_ndg_old)) allocate(t_ndg_old(1:mm,1:ndglev))
    if(.not.allocated(q_ndg_old)) allocate(q_ndg_old(1:mm,1:ndglev))
    if(.not.allocated(u_ndg_new)) allocate(u_ndg_new(1:mm,1:ndglev))
    if(.not.allocated(v_ndg_new)) allocate(v_ndg_new(1:mm,1:ndglev))
    if(.not.allocated(t_ndg_new)) allocate(t_ndg_new(1:mm,1:ndglev))
    if(.not.allocated(q_ndg_new)) allocate(q_ndg_new(1:mm,1:ndglev))
    allocate(tmplat(181))
    allocate(tmplon(360))

    do j=1,181
    tmplat(j)=(j-90.0-1.0)*1.0_8
    enddo
    do i=1,360
    tmplon(i)=(i-1)*1.0_8
    enddo

    call cmip5_to_cam_mapping(tmplon,tmplat,360,181,phys_state)

     deallocate(tmplat)
     deallocate(tmplon)

    end if

    allocate(tmpfield1_3d(360,181,ndglev))
    allocate(tmpfield2_3d(360,181,ndglev))
    allocate(camlev(pver))
    allocate(fdatao(pver))
    allocate(fdatai(ndglev))

    if(ii.eq.0.and.(.not.is_last_step())) then ! read data. At the end of the simulation, don't read data

    ntime=nnstep+1

    open(299, file='cam.nudging.bin', form='unformatted', access='direct', &
         recl=4*360*181*ndglev,action='read')
    ! u
    m=(ntime-1)*4+1
    read(299,rec=m) (((tmpfield1_3d(i,j,k),i=1,360),j=1,181),k=1,ndglev)
    m=ntime*4+1
    read(299,rec=m) (((tmpfield2_3d(i,j,k),i=1,360),j=1,181),k=1,ndglev)

    call area_cmip5_3d(tmpfield1_3d,u_ndg_old,181,360,ndglev,mm,fill_value)
    call area_cmip5_3d(tmpfield2_3d,u_ndg_new,181,360,ndglev,mm,fill_value)

    ! v
    m=(ntime-1)*4+1
    read(299,rec=m) (((tmpfield1_3d(i,j,k),i=1,plon),j=1,plat),k=1,ndglev)
    m=ntime*4+1
    read(299,rec=m) (((tmpfield2_3d(i,j,k),i=1,plon),j=1,plat),k=1,ndglev)

    call area_cmip5_3d(tmpfield1_3d,v_ndg_old,181,360,ndglev,mm,fill_value)
    call area_cmip5_3d(tmpfield2_3d,v_ndg_new,181,360,ndglev,mm,fill_value)

    ! t
    m=(ntime-1)*4+1
    read(299,rec=m) (((tmpfield1_3d(i,j,k),i=1,plon),j=1,plat),k=1,ndglev)
    m=ntime*4+1
    read(299,rec=m) (((tmpfield2_3d(i,j,k),i=1,plon),j=1,plat),k=1,ndglev)

    call area_cmip5_3d(tmpfield1_3d,t_ndg_old,181,360,ndglev,mm,fill_value)
    call area_cmip5_3d(tmpfield2_3d,t_ndg_new,181,360,ndglev,mm,fill_value)

    ! q
    m=(ntime-1)*4+1
    read(299,rec=m) (((tmpfield1_3d(i,j,k),i=1,plon),j=1,plat),k=1,ndglev)
    m=ntime*4+1
    read(299,rec=m) (((tmpfield2_3d(i,j,k),i=1,plon),j=1,plat),k=1,ndglev)

    call area_cmip5_3d(tmpfield1_3d,q_ndg_old,181,360,ndglev,mm,fill_value)
    call area_cmip5_3d(tmpfield2_3d,q_ndg_new,181,360,ndglev,mm,fill_value)

   close(299)

   endif ! end of read data

    mm=0
    do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols

     camlev(:)=phys_state(lchnk)%pmid(i,:)
     mm=mm+1

     fdatai=u_ndg_old(mm,:)
     call vert_interpolation(fdatai,fnllev,fdatao,camlev,ndglev,pver,1,vert_int,fill_value)
     phys_state(lchnk)%u_ndg_old(i,:) = fdatao(:)
     fdatai=u_ndg_new(mm,:)
     call vert_interpolation(fdatai,fnllev,fdatao,camlev,ndglev,pver,1,vert_int,fill_value)
     phys_state(lchnk)%u_ndg_new(i,:) = fdatao(:)

     fdatai=v_ndg_old(mm,:)
     call vert_interpolation(fdatai,fnllev,fdatao,camlev,ndglev,pver,1,vert_int,fill_value)
     phys_state(lchnk)%v_ndg_old(i,:) = fdatao(:)
     fdatai=v_ndg_new(mm,:)
     call vert_interpolation(fdatai,fnllev,fdatao,camlev,ndglev,pver,1,vert_int,fill_value)
     phys_state(lchnk)%v_ndg_new(i,:) = fdatao(:)

     fdatai=t_ndg_old(mm,:)
     call vert_interpolation(fdatai,fnllev,fdatao,camlev,ndglev,pver,1,vert_int,fill_value)
     phys_state(lchnk)%t_ndg_old(i,:) = fdatao(:)
     fdatai=t_ndg_new(mm,:)
     call vert_interpolation(fdatai,fnllev,fdatao,camlev,ndglev,pver,1,vert_int,fill_value)
     phys_state(lchnk)%t_ndg_new(i,:) = fdatao(:)

     fdatai=q_ndg_old(mm,:)
     call vert_interpolation(fdatai,fnllev,fdatao,camlev,ndglev,pver,4,vert_int,fill_value)
     phys_state(lchnk)%q_ndg_old(i,:) = fdatao(:)
     fdatai=q_ndg_new(mm,:)
     call vert_interpolation(fdatai,fnllev,fdatao,camlev,ndglev,pver,4,vert_int,fill_value)
     phys_state(lchnk)%q_ndg_new(i,:) = fdatao(:)

    enddo
    enddo

    if(is_last_step()) then
      if(allocated(u_ndg_old)) deallocate(u_ndg_old)
      if(allocated(v_ndg_old)) deallocate(v_ndg_old)
      if(allocated(t_ndg_old)) deallocate(t_ndg_old)
      if(allocated(q_ndg_old)) deallocate(q_ndg_old)
      if(allocated(u_ndg_new)) deallocate(u_ndg_new)
      if(allocated(v_ndg_new)) deallocate(v_ndg_new)
      if(allocated(t_ndg_new)) deallocate(t_ndg_new)
      if(allocated(q_ndg_new)) deallocate(q_ndg_new)
   end if
   deallocate(camlev)
   deallocate(fdatao)
   deallocate(fdatai)
   deallocate(tmpfield1_3d)
   deallocate(tmpfield2_3d)

   end subroutine cam_read_gfs

!----------------------------------------------------------------------
 subroutine cam_read_ccsmrcp85_6hourly( phys_state, nstep )

!----------------------------------------------------------------------- 
! Purpose: 
! Acquire and position the restart, master, primary and secondary
! datasets for a continuation run
! Author:Juanxiong He, 2012-06-01
!-----------------------------------------------------------------------
    use phys_grid, only: get_ncols_p, get_gcol_p, ngcols
    implicit none

    include 'netcdf.inc'

    type(physics_state), pointer :: phys_state(:)
    integer, intent(in) :: nstep

    integer :: lchnk, i, j, k, ii, m, n, csize, ncols, mm
    integer :: ntime, nnstep, fdda_nhtfrq
    real*8 :: dt,p0mb
    real, dimension(:,:,:,:), allocatable :: tmpfield1_3d,tmpfield2_3d
    real, dimension(:,:,:), allocatable :: tmpfield1_2d,tmpfield2_2d
    real*8, dimension(:), allocatable :: camlev,fdatao
    real*8, dimension(:), allocatable :: templat,templon
    
    integer, dimension(1:1) :: starts, scounts
    integer, dimension(1:3) :: startp, scountp
    integer, dimension(1:4) :: start, scount
    integer :: latid, lonid, levid, rhid, ncid, status
    character*19 :: yearname, yearname1
 
    integer :: yr, mon, day, tod, day1
    integer :: byr,bmon,bday,nday,mday,mmday
    integer, dimension(1:12) :: nmonth

    data nmonth/31,28,31,30,31,30,31,31,30,31,30,31/
    fill_value=1.0e30

    ! get the exact day

    call get_curr_date(yr, mon, day, tod)

    if(yr.lt.2095) then
    day=(int(get_curr_calday())-1)*4+tod/21600+1   
    nday=(int(get_curr_calday())+(yr-2005)*365)/30 ! the number of 30
    day=mod((yr-2005)*1460+day,120)  ! the order in the file
    if(day.eq.0) day=120
    if (nday*30.ge.(yr-2005)*365.and.nday*30.lt.(yr-2004)*365) then
       byr=yr
    else
       byr=yr-1
    endif
    nday=nday*30+1-(byr-2005)*365
    mday=0
    do i=1,12
    if(nday.ge.mday.and.nday.le.mday+nmonth(i)) then
       bmon=i
       bday=nday-mday       
    end if
    mday=mday+nmonth(i)
    end do
    yearname=char(byr/1000+48)//char(byr/100-byr/1000*10+48)// &
             char(byr/10-byr/100*10+48)//char(byr-byr/10*10+48)//'-'// &
             char(bmon/10-bmon/100*10+48)//char(bmon-bmon/10*10+48)//'-'//  &
             char(bday/10-bday/100*10+48)//char(bday-bday/10*10+48)//'-00000.nc'

    day1=(int(get_curr_calday())-1)*4+tod/21600+1+1
    if(day.eq.1460) then
       day1=1
       yr=yr+1
    end if
    nday=((day1-1)/4+1+(yr-2005)*365)/30  ! the number of 30
    day1=mod((yr-2005)*1460+day1,120) ! the order in the file
    if(day1.eq.0) day1=120
    if (nday*30.ge.(yr-2005)*365.and.nday*30.lt.(yr-2004)*365) then
       byr=yr
    else
       byr=yr-1
    endif
    nday=nday*30+1-(byr-2005)*365
    mday=0
    do i=1,12
    if(nday.ge.mday.and.nday.le.mday+nmonth(i)) then
       bmon=i
       bday=nday-mday
    end if
    mday=mday+nmonth(i)
    end do
    yearname1=char(byr/1000+48)//char(byr/100-byr/1000*10+48)// &
             char(byr/10-byr/100*10+48)//char(byr-byr/10*10+48)//'-'// &
             char(bmon/10-bmon/100*10+48)//char(bmon-bmon/10*10+48)//'-'//&
             char(bday/10-bday/100*10+48)//char(bday-bday/10*10+48)//'-00000.nc'
    
    else

    day=(int(get_curr_calday())-1)*4+tod/21600+1-81
    nday=(day+(yr-2095)*1460)/120 ! the number of 120
    day=mod((yr-2095)*1460+day,120)  ! the order in the file
    if(day.eq.0) day=120
    if (nday*120+81.ge.(yr-2095)*1460.and.nday*120+81.lt.(yr-2094)*1460+2) then
       byr=yr
    else
       byr=yr-1
    endif
    nday=nday*120+1-(byr-2095)*1460
    mday=0
    do i=1,12
    if(nday+81.ge.mday.and.nday+81.le.mday+nmonth(i)*4) then
       bmon=i
       bday=(nday-mday+81-1)/4+1
    end if
    mday=mday+nmonth(i)*4
    end do
    yearname=char(byr/1000+48)//char(byr/100-byr/1000*10+48)// &
             char(byr/10-byr/100*10+48)//char(byr-byr/10*10+48)//'-'// &
             char(bmon/10-bmon/100*10+48)//char(bmon-bmon/10*10+48)//'-'//  &
             char(bday/10-bday/100*10+48)//char(bday-bday/10*10+48)//'-21600.nc'

    day1=(int(get_curr_calday())-1)*4+tod/21600+1+1-81
    if(day+81.eq.1460) then
       day1=1
       yr=yr+1
    end if
    nday=(day1+(yr-2095)*1460)/120  ! the number of 30
    day1=mod((yr-2095)*1460+day1,120) ! the order in the file
    if(day1.eq.0) day1=120
    if (nday*120+81.ge.(yr-2095)*1460.and.nday*120+81.lt.(yr-2094)*1460+2) then
       byr=yr
    else
       byr=yr-1
    endif
    nday=nday*120+1-(byr-2095)*1460
    mday=0
    do i=1,12
    if(nday+81.ge.mday.and.nday+81.le.mday+nmonth(i)*4) then
       bmon=i
       bday=(nday-mday+81-1)/4+1
    end if
    mday=mday+nmonth(i)*4
    end do
    yearname1=char(byr/1000+48)//char(byr/100-byr/1000*10+48)// &
             char(byr/10-byr/100*10+48)//char(byr-byr/10*10+48)//'-'// &
             char(bmon/10-bmon/100*10+48)//char(bmon-bmon/10*10+48)//'-'//&
             char(bday/10-bday/100*10+48)//char(bday-bday/10*10+48)//'-21600.nc'

    end if
    print *,yearname1,day1

    ! get the timestep
    dt = get_step_size()
    fdda_nhtfrq=int(analysis_interval/dt)
    ii=mod(nstep,fdda_nhtfrq )

    nnstep=nstep*dt/analysis_interval
   
    allocate(camlev(pver))
    allocate(fdatao(pver))
   
    starts(1)=1

    startp(1)=1
    startp(2)=1
    scountp(3)=1

    start(1)=1
    start(2)=1
    start(3)=1
    scount(4)=1

    ! read data. At the end of the simulation, don't read data
    if(ii.eq.0.and.(.not.is_last_step())) then 
 
    ! ps
    STATUS=NF_OPEN('/R3/CCSM-rcp8.5/atm/b40.rcp8_5.1deg.007.cam2.h3.'//yearname,NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'lat', LATID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LATID, LATLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'lon', LONID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LONID, LONLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'lev', LEVID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LEVID, LEVLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    allocate(tmpfield1_3d(lonlen,latlen,levlen,1))
    allocate(tmpfield2_3d(lonlen,latlen,levlen,1))
    allocate(tmpfield1_2d(lonlen,latlen,1))
    allocate(tmpfield2_2d(lonlen,latlen,1))

    scount(1)=lonlen
    scount(2)=latlen
    scount(3)=levlen

    ! allocate and remapping at the first time
    
    ! horizontal interpolation   
    mm=0
    do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols
     n=get_gcol_p(lchnk,i)
     mm=mm+1
    end do
    end do

    if(.not.allocated(u_ndg_old)) then
    
    allocate(u_ndg_old(1:mm,1:levlen))
    if(.not.allocated(ps_ndg_old)) allocate(ps_ndg_old(mm,1:1))
    if(.not.allocated(ps_ndg_new)) allocate(ps_ndg_new(mm,1:1))
    if(.not.allocated(p_ndg_old)) allocate(p_ndg_old(mm,1:levlen))
    if(.not.allocated(p_ndg_new)) allocate(p_ndg_new(mm,1:levlen))
    if(.not.allocated(v_ndg_old)) allocate(v_ndg_old(mm,1:levlen))
    if(.not.allocated(t_ndg_old)) allocate(t_ndg_old(mm,1:levlen))
    if(.not.allocated(q_ndg_old)) allocate(q_ndg_old(mm,1:levlen))
    if(.not.allocated(u_ndg_new)) allocate(u_ndg_new(mm,1:levlen))
    if(.not.allocated(v_ndg_new)) allocate(v_ndg_new(mm,1:levlen))
    if(.not.allocated(t_ndg_new)) allocate(t_ndg_new(mm,1:levlen))
    if(.not.allocated(q_ndg_new)) allocate(q_ndg_new(mm,1:levlen))
    
    if(.not.allocated(hyam)) allocate(hyam(levlen))
    if(.not.allocated(hybm)) allocate(hybm(levlen))
    if(.not.allocated(cmiplev)) allocate(cmiplev(levlen))
    if(.not.allocated(templat)) allocate(templat(latlen))
    if(.not.allocated(templon)) allocate(templon(lonlen))

    ! lat
    scounts(1)=latlen
    STATUS=NF_INQ_VARID(NCID,'lat',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,templat)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    ! lon
    scounts(1)=lonlen
    STATUS=NF_INQ_VARID(NCID,'lon',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,templon)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    call cmip5_to_cam_mapping(templon,templat,lonlen,latlen,phys_state)
    deallocate(templon)
    deallocate(templat)
    
    end if

    ! hyam
    scounts(1)=levlen
    STATUS=NF_INQ_VARID(NCID,'hyam',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,hyam)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    ! hybm
    scounts(1)=levlen
    STATUS=NF_INQ_VARID(NCID,'hybm',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,hybm)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    ! P0
    scounts(1)=1
    STATUS=NF_INQ_VARID(NCID,'P0',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,p0mb)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    scountp(1)=lonlen
    scountp(2)=latlen
    startp(3)=day
    STATUS=NF_INQ_VARID(NCID,'PS',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTP,SCOUNTP,tmpfield1_2d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS=NF_OPEN('/R3/CCSM-rcp8.5/atm/b40.rcp8_5.1deg.007.cam2.h3.'//yearname1,NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'PS',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    startp(3)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTP,SCOUNTP,tmpfield2_2d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    call area_cmip5_2d(tmpfield1_2d,ps_ndg_old,latlen,lonlen,1,mm,fill_value)
    call area_cmip5_2d(tmpfield2_2d,ps_ndg_new,latlen,lonlen,1,mm,fill_value)

    mm=0
    do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols
     n=get_gcol_p(lchnk,i)
     mm=mm+1
     do k=1,levlen
       p_ndg_old(mm,k) = ps_ndg_old(mm,1)*hyam(k)+p0mb*hybm(k)
       p_ndg_new(mm,k) = ps_ndg_new(mm,1)*hyam(k)+p0mb*hybm(k)
     end do
    end do
    end do
  
    ! t
    STATUS=NF_OPEN('/R3/CCSM-rcp8.5/atm/b40.rcp8_5.1deg.007.cam2.h3.'//yearname,NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'T',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS=NF_OPEN('/R3/CCSM-rcp8.5/atm/b40.rcp8_5.1deg.007.cam2.h3.'//yearname1,NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'T',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield2_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    call area_cmip5_3d(tmpfield1_3d,t_ndg_old,latlen,lonlen,levlen,mm,fill_value)
    call area_cmip5_3d(tmpfield2_3d,t_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    ! u
    STATUS=NF_OPEN('/R3/CCSM-rcp8.5/atm/b40.rcp8_5.1deg.007.cam2.h3.'//yearname,NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'U',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS=NF_OPEN('/R3/CCSM-rcp8.5/atm/b40.rcp8_5.1deg.007.cam2.h3.'//yearname1,NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'U',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield2_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    call area_cmip5_3d(tmpfield1_3d,u_ndg_old,latlen,lonlen,levlen,mm,fill_value)
    call area_cmip5_3d(tmpfield2_3d,u_ndg_new,latlen,lonlen,levlen,mm,fill_value)
 
    ! v
    STATUS=NF_OPEN('/R3/CCSM-rcp8.5/atm/b40.rcp8_5.1deg.007.cam2.h3.'//yearname,NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'V',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS=NF_OPEN('/R3/CCSM-rcp8.5/atm/b40.rcp8_5.1deg.007.cam2.h3.'//yearname1,NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'V',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield2_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    call area_cmip5_3d(tmpfield1_3d,v_ndg_old,latlen,lonlen,levlen,mm,fill_value)
    call area_cmip5_3d(tmpfield2_3d,v_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    ! q
    STATUS=NF_OPEN('/R3/CCSM-rcp8.5/atm/b40.rcp8_5.1deg.007.cam2.h3.'//yearname,NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'Q',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS=NF_OPEN('/R3/CCSM-rcp8.5/atm/b40.rcp8_5.1deg.007.cam2.h3.'//yearname1,NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'Q',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield2_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    call area_cmip5_3d(tmpfield1_3d,q_ndg_old,latlen,lonlen,levlen,mm,fill_value)
    call area_cmip5_3d(tmpfield2_3d,q_ndg_new,latlen,lonlen,levlen,mm,fill_value)

   deallocate(tmpfield1_3d)
   deallocate(tmpfield2_3d)
   deallocate(tmpfield1_2d)
   deallocate(tmpfield2_2d)

   endif ! end of read data

   ! vertical interpolation
   mm=0
   do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols
     n=get_gcol_p(lchnk,i)
     mm=mm+1

     camlev(:)=phys_state(lchnk)%pmid(i,:)

     ! old state
     cmiplev(:)=p_ndg_old(mm,:)
     call vert_interpolation(u_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%u_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(v_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%v_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(t_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%t_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(q_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,4,vert_int,fill_value)
     phys_state(lchnk)%q_ndg_old(i,:) = fdatao(:)

     ! new state
     cmiplev(:)=p_ndg_new(mm,:)
     call vert_interpolation(u_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%u_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(v_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%v_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(t_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%t_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(q_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,4,vert_int,fill_value)
     phys_state(lchnk)%q_ndg_new(i,:) = fdatao(:)

    end do
   end do

   if(is_last_step()) then
      if(allocated(hyam)) deallocate(hyam)
      if(allocated(hybm)) deallocate(hybm)
      if(allocated(cmiplev)) deallocate(cmiplev)
      if(allocated(u_ndg_old)) deallocate(u_ndg_old)
      if(allocated(v_ndg_old)) deallocate(v_ndg_old)
      if(allocated(t_ndg_old)) deallocate(t_ndg_old)
      if(allocated(q_ndg_old)) deallocate(q_ndg_old)
      if(allocated(u_ndg_new)) deallocate(u_ndg_new)
      if(allocated(v_ndg_new)) deallocate(v_ndg_new)
      if(allocated(t_ndg_new)) deallocate(t_ndg_new)
      if(allocated(q_ndg_new)) deallocate(q_ndg_new)
      if(allocated(ps_ndg_old)) deallocate(ps_ndg_old)
      if(allocated(ps_ndg_new)) deallocate(ps_ndg_new)
      if(allocated(p_ndg_old)) deallocate(p_ndg_old)
      if(allocated(p_ndg_new)) deallocate(p_ndg_new)
   end if

   deallocate(camlev)
   deallocate(fdatao)

 end subroutine  cam_read_ccsmrcp85_6hourly
	
!----------------------------------------------------------------------
 subroutine cam_read_cmip5_ccsm4_6hourly_eul( phys_state, nstep )
!----------------------------------------------------------------------- 
! Purpose: 
! Acquire and position the restart, master, primary and secondary
! datasets for a continuation run
! Author:Juanxiong He, 2012-06-01
!-----------------------------------------------------------------------
    use phys_grid, only: get_ncols_p, get_gcol_p, ngcols
    implicit none

    include 'netcdf.inc'

    type(physics_state), pointer :: phys_state(:)
    integer, intent(in) :: nstep

    integer :: lchnk, i, j, k, ii, m, n, csize, ncols, mm
    integer :: ntime, nnstep, fdda_nhtfrq
    real*8 :: dt,p0mb
    real, dimension(:,:,:,:), allocatable :: tmpfield1_3d,tmpfield2_3d
    real, dimension(:,:,:), allocatable :: tmpfield1_2d,tmpfield2_2d
    real*8, dimension(:), allocatable :: camlev,fdatao
    real*8, dimension(:), allocatable :: templat,templon
    
    integer, dimension(1:1) :: starts, scounts
    integer, dimension(1:3) :: startp, scountp
    integer, dimension(1:4) :: start, scount
    integer :: latid, lonid, levid, rhid, ncid, status
    character*4 :: yearname, yearname1
 
    integer :: yr, mon, day, tod, day1

    fill_value=1.0e30

    ! get the exact day
    call get_curr_date(yr, mon, day, tod)
    day=(int(get_curr_calday())-1)*4+tod/21600+1
    yearname=char(yr/1000+48)//char(yr/100-yr/1000*10+48)//char(yr/10-yr/100*10+48)// &
             char(yr-yr/10*10+48)
    day1=day+1
    yearname1=yearname
    if(day.eq.1460) then
       day1=1
       yr=yr+1
       yearname1=char(yr/1000+48)//char(yr/100-yr/1000*10+48)//char(yr/10-yr/100*10+48)// &
                 char(yr-yr/10*10+48)
    end if
     
    ! get the timestep
    dt = get_step_size()
    fdda_nhtfrq=int(analysis_interval/dt)
    ii=mod(nstep,fdda_nhtfrq )

    nnstep=nstep*dt/analysis_interval
   
    allocate(camlev(pver))
    allocate(fdatao(pver))
   
    starts(1)=1

    startp(1)=1
    startp(2)=1
    scountp(3)=1

    start(1)=1
    start(2)=1
    start(3)=1
    scount(4)=1

    ! read data. At the end of the simulation, don't read data
    if(ii.eq.0.and.(.not.is_last_step())) then 
 
    ! ps
    STATUS=NF_OPEN('/B0/CMIP5/6hourly/CCSM4/'//yearname//'.20thC.PS.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'lat', LATID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LATID, LATLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'lon', LONID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LONID, LONLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'lev', LEVID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LEVID, LEVLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    allocate(tmpfield1_3d(lonlen,latlen,levlen,1))
    allocate(tmpfield2_3d(lonlen,latlen,levlen,1))
    allocate(tmpfield1_2d(lonlen,latlen,1))
    allocate(tmpfield2_2d(lonlen,latlen,1))

    scount(1)=lonlen
    scount(2)=latlen
    scount(3)=levlen

    ! allocate and remapping at the first time
    
    ! horizontal interpolation   
    mm=0
    do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols
     n=get_gcol_p(lchnk,i)
     mm=mm+1
    end do
    end do

    if(.not.allocated(u_ndg_old)) then
    
    allocate(u_ndg_old(1:mm,1:levlen))
    if(.not.allocated(ps_ndg_old)) allocate(ps_ndg_old(mm,1:1))
    if(.not.allocated(ps_ndg_new)) allocate(ps_ndg_new(mm,1:1))
    if(.not.allocated(p_ndg_old)) allocate(p_ndg_old(mm,1:levlen))
    if(.not.allocated(p_ndg_new)) allocate(p_ndg_new(mm,1:levlen))
    if(.not.allocated(v_ndg_old)) allocate(v_ndg_old(mm,1:levlen))
    if(.not.allocated(t_ndg_old)) allocate(t_ndg_old(mm,1:levlen))
    if(.not.allocated(q_ndg_old)) allocate(q_ndg_old(mm,1:levlen))
    if(.not.allocated(u_ndg_new)) allocate(u_ndg_new(mm,1:levlen))
    if(.not.allocated(v_ndg_new)) allocate(v_ndg_new(mm,1:levlen))
    if(.not.allocated(t_ndg_new)) allocate(t_ndg_new(mm,1:levlen))
    if(.not.allocated(q_ndg_new)) allocate(q_ndg_new(mm,1:levlen))
    
    if(.not.allocated(hyam)) allocate(hyam(levlen))
    if(.not.allocated(hybm)) allocate(hybm(levlen))
    if(.not.allocated(cmiplev)) allocate(cmiplev(levlen))
    if(.not.allocated(templat)) allocate(templat(latlen))
    if(.not.allocated(templon)) allocate(templon(lonlen))

    ! lat
    scounts(1)=latlen
    STATUS=NF_INQ_VARID(NCID,'lat',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,templat)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    ! lon
    scounts(1)=lonlen
    STATUS=NF_INQ_VARID(NCID,'lon',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,templon)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    call cmip5_to_cam_mapping(templon,templat,lonlen,latlen,phys_state)
    deallocate(templon)
    deallocate(templat)
    
    end if

    ! hyam
    scounts(1)=levlen
    STATUS=NF_INQ_VARID(NCID,'hyam',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,hyam)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    ! hybm
    scounts(1)=levlen
    STATUS=NF_INQ_VARID(NCID,'hybm',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,hybm)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    ! P0
    scounts(1)=1
    STATUS=NF_INQ_VARID(NCID,'P0',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,p0mb)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    scountp(1)=lonlen
    scountp(2)=latlen
    startp(3)=day
    STATUS=NF_INQ_VARID(NCID,'PS',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTP,SCOUNTP,tmpfield1_2d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS=NF_OPEN('/B0/CMIP5/6hourly/CCSM4/'//yearname1//'.20thC.PS.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'PS',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    startp(3)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTP,SCOUNTP,tmpfield2_2d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    call area_cmip5_2d(tmpfield1_2d,ps_ndg_old,latlen,lonlen,1,mm,fill_value)
    call area_cmip5_2d(tmpfield2_2d,ps_ndg_new,latlen,lonlen,1,mm,fill_value)

    mm=0
    do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols
     n=get_gcol_p(lchnk,i)
     mm=mm+1
     do k=1,levlen
       p_ndg_old(mm,k) = ps_ndg_old(mm,1)*hyam(k)+p0mb*hybm(k)
       p_ndg_new(mm,k) = ps_ndg_new(mm,1)*hyam(k)+p0mb*hybm(k)
     end do
    end do
    end do

    ! t
    STATUS=NF_OPEN('/B0/CMIP5/6hourly/CCSM4/'//yearname//'.20thC.T.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'T',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS=NF_OPEN('/B0/CMIP5/6hourly/CCSM4/'//yearname1//'.20thC.T.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'T',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield2_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    call area_cmip5_3d(tmpfield1_3d,t_ndg_old,latlen,lonlen,levlen,mm,fill_value)
    call area_cmip5_3d(tmpfield2_3d,t_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    ! u
    STATUS=NF_OPEN('/B0/CMIP5/6hourly/CCSM4/'//yearname//'.20thC.U.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'U',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS=NF_OPEN('/B0/CMIP5/6hourly/CCSM4/'//yearname1//'.20thC.U.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'U',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield2_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    call area_cmip5_3d(tmpfield1_3d,u_ndg_old,latlen,lonlen,levlen,mm,fill_value)
    call area_cmip5_3d(tmpfield2_3d,u_ndg_new,latlen,lonlen,levlen,mm,fill_value)
 
    ! v
    STATUS=NF_OPEN('/B0/CMIP5/6hourly/CCSM4/'//yearname//'.20thC.V.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'V',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS=NF_OPEN('/B0/CMIP5/6hourly/CCSM4/'//yearname1//'.20thC.V.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'V',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield2_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    call area_cmip5_3d(tmpfield1_3d,v_ndg_old,latlen,lonlen,levlen,mm,fill_value)
    call area_cmip5_3d(tmpfield2_3d,v_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    ! q
    STATUS=NF_OPEN('/B0/CMIP5/6hourly/CCSM4/'//yearname//'.20thC.Q.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'Q',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS=NF_OPEN('/B0/CMIP5/6hourly/CCSM4/'//yearname1//'.20thC.Q.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'Q',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield2_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    call area_cmip5_3d(tmpfield1_3d,q_ndg_old,latlen,lonlen,levlen,mm,fill_value)
    call area_cmip5_3d(tmpfield2_3d,q_ndg_new,latlen,lonlen,levlen,mm,fill_value)

   deallocate(tmpfield1_3d)
   deallocate(tmpfield2_3d)
   deallocate(tmpfield1_2d)
   deallocate(tmpfield2_2d)

   endif ! end of read data

   ! vertical interpolation
   mm=0
   do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols
     n=get_gcol_p(lchnk,i)
     mm=mm+1

     camlev(:)=phys_state(lchnk)%pmid(i,:)

     ! old state
     cmiplev(:)=p_ndg_old(mm,:)
     call vert_interpolation(u_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%u_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(v_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%v_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(t_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%t_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(q_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,4,vert_int,fill_value)
     phys_state(lchnk)%q_ndg_old(i,:) = fdatao(:)

     ! new state
     cmiplev(:)=p_ndg_new(mm,:)
     call vert_interpolation(u_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%u_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(v_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%v_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(t_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%t_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(q_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,4,vert_int,fill_value)
     phys_state(lchnk)%q_ndg_new(i,:) = fdatao(:)


    end do
   end do

   if(is_last_step()) then
      if(allocated(hyam)) deallocate(hyam)
      if(allocated(hybm)) deallocate(hybm)
      if(allocated(cmiplev)) deallocate(cmiplev)
      if(allocated(u_ndg_old)) deallocate(u_ndg_old)
      if(allocated(v_ndg_old)) deallocate(v_ndg_old)
      if(allocated(t_ndg_old)) deallocate(t_ndg_old)
      if(allocated(q_ndg_old)) deallocate(q_ndg_old)
      if(allocated(u_ndg_new)) deallocate(u_ndg_new)
      if(allocated(v_ndg_new)) deallocate(v_ndg_new)
      if(allocated(t_ndg_new)) deallocate(t_ndg_new)
      if(allocated(q_ndg_new)) deallocate(q_ndg_new)
      if(allocated(ps_ndg_old)) deallocate(ps_ndg_old)
      if(allocated(ps_ndg_new)) deallocate(ps_ndg_new)
      if(allocated(p_ndg_old)) deallocate(p_ndg_old)
      if(allocated(p_ndg_new)) deallocate(p_ndg_new)
   end if

   deallocate(camlev)
   deallocate(fdatao)

   end subroutine cam_read_cmip5_ccsm4_6hourly_eul
!----------------------------------------------------------------------
 
   subroutine cam_read_cmip5_ccsm4_daily_eul( phys_state, nstep )
!----------------------------------------------------------------------- 
! Purpose: 
! Acquire and position the restart, master, primary and secondary
! datasets for a continuation run
! Author:Juanxiong He, 2012-06-01
!-----------------------------------------------------------------------
    use phys_grid, only: get_ncols_p, get_gcol_p, ngcols
    implicit none

    include 'netcdf.inc'

    type(physics_state), pointer :: phys_state(:)
    integer, intent(in) :: nstep

    integer :: lchnk, i, j, k, ii, m, n, csize, ncols, mm
    integer :: ntime, nnstep, fdda_nhtfrq
    real*8 :: dt
    real, dimension(:,:,:,:), allocatable :: tmpfield1_3d,tmpfield2_3d,tmpfield3_3d,tmpfield4_3d    
    real*8, dimension(:), allocatable :: camlev,fdatao
    real*8, dimension(:), allocatable :: templat,templon
    
    integer, dimension(1:1) :: starts, scounts
    integer, dimension(1:3) :: startp, scountp
    integer, dimension(1:4) :: start, scount
    integer :: latid, lonid, levid, rhid, ncid, status
    
    character*17, dimension(1:6) :: yearname
    integer :: yr, mon, day, tod, myyr, day1, myyr1

    data yearname/'19750101-19791231','19800101-19841231','19850101-19891231',&
                  '19900101-19941231','19950101-19991231','20000101-20051231'/
    
    fill_value=1.0e20
      
    ! get the exact day
    call get_curr_date(yr, mon, day, tod)
    day=int(get_curr_calday())
    if(yr.ge.1975.and.yr.le.1979) then
    	day=(yr-1975)*365+day
    	myyr=1 
    elseif(yr.ge.1980.and.yr.le.1984) then
    	day=(yr-1980)*365+day
    	myyr=2	
    elseif(yr.ge.1985.and.yr.le.1989) then
    	day=(yr-1985)*365+day
    	myyr=3
    elseif(yr.ge.1990.and.yr.le.1994) then
    	day=(yr-1990)*365+day
    	myyr=4
    elseif(yr.ge.1995.and.yr.le.1999) then
    	day=(yr-1995)*365+day
    	myyr=5
    elseif(yr.ge.2000.and.yr.le.2005) then
    	day=(yr-2000)*365+day
    	myyr=6	
    endif	
    
    day1=int(get_curr_calday())+1
    if(day1.lt.366) then
        day1=day+1
        myyr1=myyr
    else
     if(yr.eq.1979.or.yr.eq.1984.or.yr.eq.1989.or.yr.eq.1994.or. &
        yr.eq.1999.or.yr.eq.2005) then
        day1=1
        myyr1=myyr+1
     else
        day1=day+1
        myyr1=myyr
     endif
    endif

    ! get the timestep
    dt = get_step_size()
    fdda_nhtfrq=int(analysis_interval/dt)
    ii=mod(nstep,fdda_nhtfrq )

    nnstep=nstep*dt/analysis_interval
   
    allocate(camlev(pver))
    allocate(fdatao(pver))
   
    starts(1)=1

    startp(1)=1
    startp(2)=1
    scountp(3)=1

    start(1)=1
    start(2)=1
    start(3)=1
    scount(4)=1

    ! read data. At the end of the simulation, don't read data
    if(ii.eq.0.and.(.not.is_last_step())) then 
 
    ! t
    STATUS=NF_OPEN('/B3/cmip5/daily/CCSM4/ta_day_CCSM4_historical_r6i1p1_'//yearname(myyr)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'lat', LATID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LATID, LATLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'lon', LONID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LONID, LONLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'plev', LEVID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LEVID, LEVLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    if(.not.allocated(tmpfield1_3d)) allocate(tmpfield1_3d(lonlen,latlen,levlen,1))
    if(.not.allocated(tmpfield2_3d)) allocate(tmpfield2_3d(lonlen,latlen,levlen,1))
    if(.not.allocated(tmpfield3_3d)) allocate(tmpfield3_3d(lonlen,latlen,levlen,1))
    if(.not.allocated(tmpfield4_3d)) allocate(tmpfield4_3d(lonlen,latlen,levlen,1))

    scount(1)=lonlen
    scount(2)=latlen
    scount(3)=levlen

    ! allocate and remapping at the first time
    ! horizontal interpolation   
    mm=0
    do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols
     n=get_gcol_p(lchnk,i)
     mm=mm+1
    end do
    end do

    if(.not.allocated(u_ndg_old)) then
    allocate(u_ndg_old(1:mm,1:levlen))
    if(.not.allocated(v_ndg_old)) allocate(v_ndg_old(mm,1:levlen))
    if(.not.allocated(t_ndg_old)) allocate(t_ndg_old(mm,1:levlen))
    if(.not.allocated(q_ndg_old)) allocate(q_ndg_old(mm,1:levlen))
    if(.not.allocated(u_ndg_new)) allocate(u_ndg_new(mm,1:levlen))
    if(.not.allocated(v_ndg_new)) allocate(v_ndg_new(mm,1:levlen))
    if(.not.allocated(t_ndg_new)) allocate(t_ndg_new(mm,1:levlen))
    if(.not.allocated(q_ndg_new)) allocate(q_ndg_new(mm,1:levlen))
    if(.not.allocated(cmiplev)) allocate(cmiplev(levlen))
    if(.not.allocated(templat)) allocate(templat(latlen))
    if(.not.allocated(templon)) allocate(templon(lonlen))

    ! lat
    scounts(1)=latlen
    STATUS=NF_INQ_VARID(NCID,'lat',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,templat)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    ! lon
    scounts(1)=lonlen
    STATUS=NF_INQ_VARID(NCID,'lon',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,templon)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    ! lev
    scounts(1)=levlen
    STATUS=NF_INQ_VARID(NCID,'plev',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,cmiplev)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    cmiplev=cmiplev(levlen:1:-1)  ! reverse to the order from high to low
    
    call cmip5_to_cam_mapping(templon,templat,lonlen,latlen,phys_state)
    
    end if

    ! t
    STATUS=NF_INQ_VARID(NCID,'ta',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:) 
    call area_cmip5_3d(tmpfield1_3d,t_ndg_old,latlen,lonlen,levlen,mm,fill_value)

    STATUS=NF_OPEN('/B3/cmip5/daily/CCSM4/ta_day_CCSM4_historical_r6i1p1_'//yearname(myyr1)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'ta',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)
    call area_cmip5_3d(tmpfield1_3d,t_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    ! u
    STATUS=NF_OPEN('/B3/cmip5/daily/CCSM4/ua_day_CCSM4_historical_r6i1p1_'//yearname(myyr)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'ua',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:) 
    call area_cmip5_3d(tmpfield1_3d,u_ndg_old,latlen,lonlen,levlen,mm,fill_value)

    STATUS=NF_OPEN('/B3/cmip5/daily/CCSM4/ua_day_CCSM4_historical_r6i1p1_'//yearname(myyr1)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'ua',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:) 
    call area_cmip5_3d(tmpfield1_3d,u_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    ! v
    STATUS=NF_OPEN('/B3/cmip5/daily/CCSM4/va_day_CCSM4_historical_r6i1p1_'//yearname(myyr)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'va',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield2_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield2_3d=tmpfield2_3d(:,:,levlen:1:-1,:) 
    call area_cmip5_3d(tmpfield2_3d,v_ndg_old,latlen,lonlen,levlen,mm,fill_value)

    STATUS=NF_OPEN('/B3/cmip5/daily/CCSM4/va_day_CCSM4_historical_r6i1p1_'//yearname(myyr1)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'va',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield2_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield2_3d=tmpfield2_3d(:,:,levlen:1:-1,:)
    call area_cmip5_3d(tmpfield2_3d,v_ndg_new,latlen,lonlen,levlen,mm,fill_value)

   deallocate(tmpfield1_3d)
   deallocate(tmpfield2_3d)
   deallocate(tmpfield3_3d)
   deallocate(tmpfield4_3d)

   endif ! end of read data

   ! vertical interpolation
   mm=0
   do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols
     n=get_gcol_p(lchnk,i)
     mm=mm+1

     camlev(:)=phys_state(lchnk)%pmid(i,:)

     ! old state
     call vert_interpolation(u_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%u_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(u_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%u_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(v_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%v_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(v_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%v_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(t_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%t_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(t_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%t_ndg_new(i,:) = fdatao(:)

     phys_state(lchnk)%q_ndg_old(i,:) = -9999.0
     phys_state(lchnk)%q_ndg_new(i,:) = -9999.0

    end do
   end do

   if(is_last_step()) then
      if(allocated(cmiplev)) deallocate(cmiplev)
      if(allocated(u_ndg_old)) deallocate(u_ndg_old)
      if(allocated(v_ndg_old)) deallocate(v_ndg_old)
      if(allocated(t_ndg_old)) deallocate(t_ndg_old)
      if(allocated(q_ndg_old)) deallocate(q_ndg_old)
      if(allocated(u_ndg_new)) deallocate(u_ndg_new)
      if(allocated(v_ndg_new)) deallocate(v_ndg_new)
      if(allocated(t_ndg_new)) deallocate(t_ndg_new)
      if(allocated(q_ndg_new)) deallocate(q_ndg_new)
      if(allocated(templat)) deallocate(templat)
      if(allocated(templon)) deallocate(templon)
   end if

   deallocate(camlev)
   deallocate(fdatao)
   
   end subroutine cam_read_cmip5_ccsm4_daily_eul
   
!----------------------------------------------------------------------
 
   subroutine cam_read_cmip5_gfdlesm2m_daily_eul( phys_state, nstep )
!----------------------------------------------------------------------- 
! Purpose: 
! Acquire and position the restart, master, primary and secondary
! datasets for a continuation run
! Author:Juanxiong He, 2012-06-01
!-----------------------------------------------------------------------
    use phys_grid, only: get_ncols_p, get_gcol_p, ngcols
    implicit none

    include 'netcdf.inc'

    type(physics_state), pointer :: phys_state(:)
    integer, intent(in) :: nstep

    integer :: lchnk, i, j, k, ii, m, n, csize, ncols, mm
    integer :: ntime, nnstep, fdda_nhtfrq
    real*8 :: dt
    real, dimension(:,:,:,:), allocatable :: tmpfield1_3d,tmpfield2_3d,tmpfield3_3d,tmpfield4_3d
    real*8, dimension(:), allocatable :: camlev,fdatao
    real*8, dimension(:), allocatable :: templat,templon
    
    integer, dimension(1:1) :: starts, scounts
    integer, dimension(1:3) :: startp, scountp
    integer, dimension(1:4) :: start, scount
    integer :: latid, lonid, levid, rhid, ncid, status
    
    character*17, dimension(1:6) :: yearname
    integer :: yr, mon, day, tod, myyr, day1, myyr1

    data yearname/'19760101-19801231','19810101-19851231','19860101-19901231',&
                  '19910101-19951231','19960101-20001231','20010101-20051231'/
          
    fill_value=1.0e20
    ! get the exact day
    call get_curr_date(yr, mon, day, tod)
    day=int(get_curr_calday())
    if(yr.ge.1976.and.yr.le.1980) then
    	day=(yr-1976)*365+day
    	myyr=1 
    elseif(yr.ge.1981.and.yr.le.1985) then
    	day=(yr-1981)*365+day
    	myyr=2	
    elseif(yr.ge.1986.and.yr.le.1990) then
    	day=(yr-1986)*365+day
    	myyr=3
    elseif(yr.ge.1991.and.yr.le.1995) then
    	day=(yr-1991)*365+day
    	myyr=4
    elseif(yr.ge.1996.and.yr.le.2000) then
    	day=(yr-1996)*365+day
    	myyr=5
    elseif(yr.ge.2001.and.yr.le.2005) then
    	day=(yr-2001)*365+day
    	myyr=6	
    endif	
    
    day1=int(get_curr_calday())+1
    if(day1.lt.366) then
        day1=day+1
        myyr1=myyr
    else
     if(yr.eq.1980.or.yr.eq.1985.or.yr.eq.1990.or.yr.eq.1995.or. &
        yr.eq.2000.or.yr.eq.2005) then
        day1=1
        myyr1=myyr+1
     else
        day1=day+1
        myyr1=myyr
     endif
    endif

    ! get the timestep
    dt = get_step_size()
    fdda_nhtfrq=int(analysis_interval/dt)
    ii=mod(nstep,fdda_nhtfrq )

    nnstep=nstep*dt/analysis_interval

    allocate(camlev(pver))
    allocate(fdatao(pver))
   
    starts(1)=1

    startp(1)=1
    startp(2)=1
    scountp(3)=1

    start(1)=1
    start(2)=1
    start(3)=1
    scount(4)=1

    ! read data. At the end of the simulation, don't read data
    if(ii.eq.0.and.(.not.is_last_step())) then 
 
    ! t
    STATUS=NF_OPEN('/B3/cmip5/daily/GFDL-ESM2M/ta_day_GFDL-ESM2M_historical_r1i1p1_'//yearname(myyr)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'lat', LATID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LATID, LATLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'lon', LONID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LONID, LONLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'plev', LEVID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LEVID, LEVLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    if(.not.allocated(tmpfield1_3d)) allocate(tmpfield1_3d(lonlen,latlen,levlen,1))
    if(.not.allocated(tmpfield2_3d)) allocate(tmpfield2_3d(lonlen,latlen,levlen,1))
    if(.not.allocated(tmpfield3_3d)) allocate(tmpfield3_3d(lonlen,latlen,levlen,1))
    if(.not.allocated(tmpfield4_3d)) allocate(tmpfield4_3d(lonlen,latlen,levlen,1))

    scount(1)=lonlen
    scount(2)=latlen
    scount(3)=levlen

    ! allocate and remapping at the first time
    ! horizontal interpolation   
    mm=0
    do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols
     n=get_gcol_p(lchnk,i)
     mm=mm+1
    end do
    end do

    if(.not.allocated(u_ndg_old)) then

    allocate(u_ndg_old(1:mm,1:levlen))
    if(.not.allocated(v_ndg_old)) allocate(v_ndg_old(mm,1:levlen))
    if(.not.allocated(t_ndg_old)) allocate(t_ndg_old(mm,1:levlen))
    if(.not.allocated(q_ndg_old)) allocate(q_ndg_old(mm,1:levlen))
    if(.not.allocated(u_ndg_new)) allocate(u_ndg_new(mm,1:levlen))
    if(.not.allocated(v_ndg_new)) allocate(v_ndg_new(mm,1:levlen))
    if(.not.allocated(t_ndg_new)) allocate(t_ndg_new(mm,1:levlen))
    if(.not.allocated(q_ndg_new)) allocate(q_ndg_new(mm,1:levlen))

    if(.not.allocated(cmiplev)) allocate(cmiplev(levlen))
    if(.not.allocated(templat)) allocate(templat(latlen))
    if(.not.allocated(templon)) allocate(templon(lonlen))

    ! lat
    scounts(1)=latlen
    STATUS=NF_INQ_VARID(NCID,'lat',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,templat)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    ! lon
    scounts(1)=lonlen
    STATUS=NF_INQ_VARID(NCID,'lon',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,templon)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    ! lev
    scounts(1)=levlen
    STATUS=NF_INQ_VARID(NCID,'plev',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,cmiplev)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    cmiplev=cmiplev(levlen:1:-1)  ! reverse to the order from high to low
    
    call cmip5_to_cam_mapping(templon,templat,lonlen,latlen,phys_state)
    
    end if

    ! t
    STATUS=NF_INQ_VARID(NCID,'ta',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)
    call area_cmip5_3d(tmpfield1_3d,t_ndg_old,latlen,lonlen,levlen,mm,fill_value)
 
    STATUS=NF_OPEN('/B3/cmip5/daily/GFDL-ESM2M/ta_day_GFDL-ESM2M_historical_r1i1p1_'//yearname(myyr1)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'ta',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:) 
    call area_cmip5_3d(tmpfield1_3d,t_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    ! u
    STATUS=NF_OPEN('/B3/cmip5/daily/GFDL-ESM2M/ua_day_GFDL-ESM2M_historical_r1i1p1_'//yearname(myyr)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'ua',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:) 
    call area_cmip5_3d(tmpfield1_3d,u_ndg_old,latlen,lonlen,levlen,mm,fill_value)

    STATUS=NF_OPEN('/B3/cmip5/daily/GFDL-ESM2M/ua_day_GFDL-ESM2M_historical_r1i1p1_'//yearname(myyr1)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'ua',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)
    call area_cmip5_3d(tmpfield1_3d,u_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    ! v
    STATUS=NF_OPEN('/B3/cmip5/daily/GFDL-ESM2M/va_day_GFDL-ESM2M_historical_r1i1p1_'//yearname(myyr)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'va',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield2_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield2_3d=tmpfield2_3d(:,:,levlen:1:-1,:) 
    call area_cmip5_3d(tmpfield2_3d,v_ndg_old,latlen,lonlen,levlen,mm,fill_value)

    STATUS=NF_OPEN('/B3/cmip5/daily/GFDL-ESM2M/va_day_GFDL-ESM2M_historical_r1i1p1_'//yearname(myyr1)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'va',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield2_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield2_3d=tmpfield2_3d(:,:,levlen:1:-1,:)
    call area_cmip5_3d(tmpfield2_3d,v_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    ! q
    STATUS=NF_OPEN('/B3/cmip5/daily/GFDL-ESM2M/hur_day_GFDL-ESM2M_historical_r1i1p1_'//yearname(myyr)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'hur',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:) 
    call area_cmip5_3d(tmpfield1_3d,q_ndg_old,latlen,lonlen,levlen,mm,fill_value)

    STATUS=NF_OPEN('/B3/cmip5/daily/GFDL-ESM2M/hur_day_GFDL-ESM2M_historical_r1i1p1_'//yearname(myyr1)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'hur',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)
    call area_cmip5_3d(tmpfield1_3d,q_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    do j=1,levlen
    do i=1,mm
    if(q_ndg_old(i,j).ne.fill_value.and.t_ndg_old(i,j).ne.fill_value) then
    q_ndg_old(i,j)=exp((t_ndg_old(i,j)-273.15)*17.67/(t_ndg_old(i,j)-29.65))* &
                   611.2*q_ndg_old(i,j)/100
    q_ndg_old(i,j)=0.622*q_ndg_old(i,j)/(cmiplev(j)-0.378*q_ndg_old(i,j))
    if(q_ndg_old(i,j).le.0) q_ndg_old(i,j)=1d-12
    else
    q_ndg_old(i,j)=fill_value
    endif
    
    if(q_ndg_new(i,j).ne.fill_value.and.t_ndg_new(i,j).ne.fill_value) then
    q_ndg_new(i,j)=exp((t_ndg_new(i,j)-273.15)*17.67/(t_ndg_new(i,j)-29.65))* &
                   611.2*q_ndg_new(i,j)/100
    q_ndg_new(i,j)=0.622*q_ndg_new(i,j)/(cmiplev(j)-0.378*q_ndg_new(i,j))
    if(q_ndg_new(i,j).le.0) q_ndg_new(i,j)=1d-12
    else
    q_ndg_new(i,j)=fill_value
    endif
    end do
    end do

   deallocate(tmpfield1_3d)
   deallocate(tmpfield2_3d)
   deallocate(tmpfield3_3d)
   deallocate(tmpfield4_3d)

   endif ! end of read data

   ! vertical interpolation
   mm=0
   do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols
     n=get_gcol_p(lchnk,i)
     mm=mm+1

     camlev(:)=phys_state(lchnk)%pmid(i,:)

     ! old state
     call vert_interpolation(u_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%u_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(u_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%u_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(v_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%v_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(v_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%v_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(t_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%t_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(t_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%t_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(q_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,4,vert_int,fill_value)
     phys_state(lchnk)%q_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(q_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,4,vert_int,fill_value)
     phys_state(lchnk)%q_ndg_new(i,:) = fdatao(:)

    end do
   end do

   if(is_last_step()) then
      if(allocated(cmiplev)) deallocate(cmiplev)
      if(allocated(u_ndg_old)) deallocate(u_ndg_old)
      if(allocated(v_ndg_old)) deallocate(v_ndg_old)
      if(allocated(t_ndg_old)) deallocate(t_ndg_old)
      if(allocated(q_ndg_old)) deallocate(q_ndg_old)
      if(allocated(u_ndg_new)) deallocate(u_ndg_new)
      if(allocated(v_ndg_new)) deallocate(v_ndg_new)
      if(allocated(t_ndg_new)) deallocate(t_ndg_new)
      if(allocated(q_ndg_new)) deallocate(q_ndg_new)
      if(allocated(templat)) deallocate(templat)
      if(allocated(templon)) deallocate(templon)
   end if

   deallocate(camlev)
   deallocate(fdatao)
   
   end subroutine cam_read_cmip5_gfdlesm2m_daily_eul
   
!----------------------------------------------------------------------
 
   subroutine cam_read_cmip5_gfdlesm2m_6hourly_eul( phys_state, nstep )
!----------------------------------------------------------------------- 
! Purpose: 
! Acquire and position the restart, master, primary and secondary
! datasets for a continuation run
! Author:Juanxiong He, 2012-06-01
!-----------------------------------------------------------------------
    use phys_grid, only: get_ncols_p, get_gcol_p, ngcols
    implicit none

    include 'netcdf.inc'

    type(physics_state), pointer :: phys_state(:)
    integer, intent(in) :: nstep

    integer :: lchnk, i, j, k, ii, m, n, csize, ncols, mm
    integer :: ntime, nnstep, fdda_nhtfrq
    real*8 :: dt
    real, dimension(:,:,:,:), allocatable :: tmpfield1_3d,tmpfield2_3d,tmpfield3_3d,tmpfield4_3d
    real, dimension(:,:,:), allocatable :: tmpfield1_2d,tmpfield2_2d
    real*8, dimension(:), allocatable :: camlev,fdatao,a,b
    real*8, dimension(:), allocatable :: templat,templon
    
    integer, dimension(1:1) :: starts, scounts
    integer, dimension(1:3) :: startp, scountp
    integer, dimension(1:4) :: start, scount
    integer :: latid, lonid, levid, rhid, ncid, status
    
    real :: p0mb
    character*21, dimension(1:6) :: yearname
    integer :: yr, mon, day, tod, myyr, day1, myyr1

    data yearname/'1976010100-1980123123','1981010100-1985123123','1986010100-1990123123',&
                  '1991010100-1995123123','1996010100-2000123123','2001010100-2005123123'/
          
    fill_value=1.0e20
    ! get the exact day
    call get_curr_date(yr, mon, day, tod)

    day=(int(get_curr_calday())-1)*4+tod/21600+1
    if(yr.ge.1976.and.yr.le.1980) then
    	day=(yr-1976)*365*4+day
    	myyr=1 
    elseif(yr.ge.1981.and.yr.le.1985) then
    	day=(yr-1981)*365*4+day
    	myyr=2	
    elseif(yr.ge.1986.and.yr.le.1990) then
    	day=(yr-1986)*365*4+day
    	myyr=3
    elseif(yr.ge.1991.and.yr.le.1995) then
    	day=(yr-1991)*365*4+day
    	myyr=4
    elseif(yr.ge.1996.and.yr.le.2000) then
    	day=(yr-1996)*365*4+day
    	myyr=5
    elseif(yr.ge.2001.and.yr.le.2005) then
    	day=(yr-2001)*365*4+day
    	myyr=6	
    endif	

    day1=(int(get_curr_calday())-1)*4+tod/21600+1+1    
    if(day1.le.1460) then
        myyr1=myyr
    else
     if(yr.eq.1980.or.yr.eq.1985.or.yr.eq.1990.or.yr.eq.1995.or. &
        yr.eq.2000.or.yr.eq.2005) then
        day1=1
        myyr1=myyr+1
     else
        day1=day+1
        myyr1=myyr
     endif
    endif

    ! get the timestep
    dt = get_step_size()
    fdda_nhtfrq=int(analysis_interval/dt)
    ii=mod(nstep,fdda_nhtfrq )

    nnstep=nstep*dt/analysis_interval

    allocate(camlev(pver))
    allocate(fdatao(pver))
   
    starts(1)=1

    startp(1)=1
    startp(2)=1
    scountp(3)=1

    start(1)=1
    start(2)=1
    start(3)=1
    scount(4)=1

    ! read data. At the end of the simulation, don't read data
    if(ii.eq.0.and.(.not.is_last_step())) then 
 
    ! t
    STATUS=NF_OPEN('/R0/jhe/aaa/ta_6hrLev_GFDL-ESM2M_historical_r1i1p1_'//yearname(myyr)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'lat', LATID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LATID, LATLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'lon', LONID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LONID, LONLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'lev', LEVID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LEVID, LEVLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    if(.not.allocated(tmpfield1_2d)) allocate(tmpfield1_2d(lonlen,latlen,1))
    if(.not.allocated(tmpfield2_2d)) allocate(tmpfield2_2d(lonlen,latlen,1))
    if(.not.allocated(tmpfield1_3d)) allocate(tmpfield1_3d(lonlen,latlen,levlen,1))
    if(.not.allocated(tmpfield2_3d)) allocate(tmpfield2_3d(lonlen,latlen,levlen,1))
    if(.not.allocated(tmpfield3_3d)) allocate(tmpfield3_3d(lonlen,latlen,levlen,1))
    if(.not.allocated(tmpfield4_3d)) allocate(tmpfield4_3d(lonlen,latlen,levlen,1))

    scount(1)=lonlen
    scount(2)=latlen
    scount(3)=levlen

    scountp(1)=lonlen
    scountp(2)=latlen

    ! allocate and remapping at the first time
    ! horizontal interpolation   
    mm=0
    do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols
     n=get_gcol_p(lchnk,i)
     mm=mm+1
    end do
    end do

    if(.not.allocated(u_ndg_old)) then

    allocate(u_ndg_old(1:mm,1:levlen))    
    if(.not.allocated(v_ndg_old)) allocate(v_ndg_old(mm,1:levlen))
    if(.not.allocated(t_ndg_old)) allocate(t_ndg_old(mm,1:levlen))
    if(.not.allocated(q_ndg_old)) allocate(q_ndg_old(mm,1:levlen))
    if(.not.allocated(u_ndg_new)) allocate(u_ndg_new(mm,1:levlen))
    if(.not.allocated(v_ndg_new)) allocate(v_ndg_new(mm,1:levlen))
    if(.not.allocated(t_ndg_new)) allocate(t_ndg_new(mm,1:levlen))
    if(.not.allocated(q_ndg_new)) allocate(q_ndg_new(mm,1:levlen))
    if(.not.allocated(p_ndg_old)) allocate(p_ndg_old(mm,1:levlen))
    if(.not.allocated(p_ndg_new)) allocate(p_ndg_new(mm,1:levlen))
    if(.not.allocated(ps_ndg_old)) allocate(ps_ndg_old(mm,1:1))
    if(.not.allocated(ps_ndg_new)) allocate(ps_ndg_new(mm,1:1))

    if(.not.allocated(cmiplev)) allocate(cmiplev(levlen))
    if(.not.allocated(templat)) allocate(templat(latlen))
    if(.not.allocated(templon)) allocate(templon(lonlen))
    if(.not.allocated(hyam)) allocate(hyam(levlen))
    if(.not.allocated(hybm)) allocate(hybm(levlen))

    ! lat
    scounts(1)=latlen
    STATUS=NF_INQ_VARID(NCID,'lat',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,templat)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    ! lon
    scounts(1)=lonlen
    STATUS=NF_INQ_VARID(NCID,'lon',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,templon)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    call cmip5_to_cam_mapping(templon,templat,lonlen,latlen,phys_state)
    deallocate(templon)
    deallocate(templat)

    end if

    ! lev
    scounts(1)=levlen
    STATUS=NF_INQ_VARID(NCID,'a',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,hyam)
    STATUS=NF_INQ_VARID(NCID,'b',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,hybm)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    ! p0mb
    scounts(1)=1
    STATUS=NF_INQ_VARID(NCID,'p0',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTS,SCOUNTS,p0mb)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    ! t
    STATUS=NF_INQ_VARID(NCID,'ta',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)
    call area_cmip5_3d(tmpfield1_3d,t_ndg_old,latlen,lonlen,levlen,mm,fill_value)
 
    STATUS=NF_OPEN('/R0/jhe/aaa/ta_6hrLev_GFDL-ESM2M_historical_r1i1p1_'//yearname(myyr1)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'ta',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)
    call area_cmip5_3d(tmpfield1_3d,t_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    ! u
    STATUS=NF_OPEN('/R0/jhe/aaa/ua_6hrLev_GFDL-ESM2M_historical_r1i1p1_'//yearname(myyr)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'ua',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)
    call area_cmip5_3d(tmpfield1_3d,u_ndg_old,latlen,lonlen,levlen,mm,fill_value)

    STATUS=NF_OPEN('/R0/jhe/aaa/ua_6hrLev_GFDL-ESM2M_historical_r1i1p1_'//yearname(myyr1)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'ua',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)
    call area_cmip5_3d(tmpfield1_3d,u_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    ! v
    STATUS=NF_OPEN('/R0/jhe/aaa/va_6hrLev_GFDL-ESM2M_historical_r1i1p1_'//yearname(myyr)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'va',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)
    call area_cmip5_3d(tmpfield2_3d,v_ndg_old,latlen,lonlen,levlen,mm,fill_value)

    STATUS=NF_OPEN('/R0/jhe/aaa/va_6hrLev_GFDL-ESM2M_historical_r1i1p1_'//yearname(myyr1)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'va',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)
    call area_cmip5_3d(tmpfield2_3d,v_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    ! ps
    STATUS=NF_OPEN('/R0/jhe/aaa/ps_6hrLev_GFDL-ESM2M_historical_r1i1p1_'//yearname(myyr)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'ps',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    startp(3)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTP,SCOUNTP,tmpfield1_2d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    call area_cmip5_2d(tmpfield1_2d,ps_ndg_old,latlen,lonlen,1,mm,fill_value)

    STATUS=NF_OPEN('/R0/jhe/aaa/ps_6hrLev_GFDL-ESM2M_historical_r1i1p1_'//yearname(myyr1)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'ps',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    startp(3)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTP,SCOUNTP,tmpfield2_2d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    call area_cmip5_2d(tmpfield2_2d,ps_ndg_new,latlen,lonlen,1,mm,fill_value)

    mm=0
    do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols
     n=get_gcol_p(lchnk,i)
     mm=mm+1
     do k=1,levlen
       p_ndg_old(mm,k) = ps_ndg_old(mm,1)*hybm(k)+p0mb*hyam(k)
       p_ndg_new(mm,k) = ps_ndg_new(mm,1)*hybm(k)+p0mb*hyam(k)
     end do
    end do
    end do
    p_ndg_old = p_ndg_old(:,levlen:1:-1)
    p_ndg_new = p_ndg_new(:,levlen:1:-1)
    
    ! q
    STATUS=NF_OPEN('/R0/jhe/aaa/hus_6hrLev_GFDL-ESM2M_historical_r1i1p1_'//yearname(myyr)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'hus',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)
    call area_cmip5_3d(tmpfield1_3d,q_ndg_old,latlen,lonlen,levlen,mm,fill_value)

    STATUS=NF_OPEN('/R0/jhe/aaa/hus_6hrLev_GFDL-ESM2M_historical_r1i1p1_'//yearname(myyr1)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'hus',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)
    call area_cmip5_3d(tmpfield1_3d,q_ndg_new,latlen,lonlen,levlen,mm,fill_value)

   deallocate(tmpfield1_2d)
   deallocate(tmpfield2_2d)
   deallocate(tmpfield1_3d)
   deallocate(tmpfield2_3d)
   deallocate(tmpfield3_3d)
   deallocate(tmpfield4_3d)

   endif ! end of read data

   ! vertical interpolation
   mm=0
   do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols
     n=get_gcol_p(lchnk,i)
     mm=mm+1

     camlev(:)=phys_state(lchnk)%pmid(i,:)

     ! old state
     cmiplev(:)=p_ndg_old(mm,:)
     call vert_interpolation(u_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%u_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(v_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%v_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(t_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%t_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(q_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,4,vert_int,fill_value)
     phys_state(lchnk)%q_ndg_old(i,:) = fdatao(:)

     ! new state
     cmiplev(:)=p_ndg_new(mm,:)
     call vert_interpolation(u_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%u_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(v_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%v_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(t_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%t_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(q_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,4,vert_int,fill_value)
     phys_state(lchnk)%q_ndg_new(i,:) = fdatao(:)

    end do
   end do

   if(is_last_step()) then
      if(allocated(cmiplev)) deallocate(cmiplev)
      if(allocated(u_ndg_old)) deallocate(u_ndg_old)
      if(allocated(v_ndg_old)) deallocate(v_ndg_old)
      if(allocated(t_ndg_old)) deallocate(t_ndg_old)
      if(allocated(q_ndg_old)) deallocate(q_ndg_old)
      if(allocated(p_ndg_old)) deallocate(p_ndg_old)
      if(allocated(ps_ndg_old)) deallocate(ps_ndg_old)
      if(allocated(u_ndg_new)) deallocate(u_ndg_new)
      if(allocated(v_ndg_new)) deallocate(v_ndg_new)
      if(allocated(t_ndg_new)) deallocate(t_ndg_new)
      if(allocated(q_ndg_new)) deallocate(q_ndg_new)
      if(allocated(p_ndg_new)) deallocate(p_ndg_new)
      if(allocated(ps_ndg_new)) deallocate(ps_ndg_new)
      if(allocated(hyam)) deallocate(hyam)
      if(allocated(hybm)) deallocate(hybm)
   end if

   deallocate(camlev)
   deallocate(fdatao)
   
   end subroutine cam_read_cmip5_gfdlesm2m_6hourly_eul
   
!----------------------------------------------------------------------
 
   subroutine cam_read_cmip5_mpiesm2m_daily_eul( phys_state, nstep )
!----------------------------------------------------------------------- 
! Purpose: 
! Acquire and position the restart, master, primary and secondary
! datasets for a continuation run
! Author:Juanxiong He, 2012-06-01
!-----------------------------------------------------------------------
    use phys_grid, only: get_ncols_p, get_gcol_p, ngcols
    implicit none

    include 'netcdf.inc'

    type(physics_state), pointer :: phys_state(:)
    integer, intent(in) :: nstep

    integer :: lchnk, i, j, k, ii, m, n, csize, ncols, mm
    integer :: ntime, nnstep, fdda_nhtfrq
    real*8 :: dt
    real, dimension(:,:,:,:), allocatable :: tmpfield1_3d,tmpfield2_3d,tmpfield3_3d,tmpfield4_3d
    real*8, dimension(:), allocatable :: camlev,fdatao
    real*8, dimension(:), allocatable :: templat,templon
    
    integer, dimension(1:1) :: starts, scounts
    integer, dimension(1:3) :: startp, scountp
    integer, dimension(1:4) :: start, scount
    integer :: latid, lonid, levid, rhid, ncid, status
    
    character*17, dimension(1:71) :: yearname
    integer :: yr, mon, day, tod, myyr, day1, myyr1

    data yearname/'19500101-19501231','19510101-19511231','19520101-19521231',&
    '19530101-19531231','19540101-19541231','19550101-19551231','19560101-19561231',&
    '19570101-19571231','19580101-19581231','19590101-19591231','19600101-19601231',&
    '19610101-19611231','19620101-19621231','19630101-19631231','19640101-19641231',&
    '19650101-19651231','19660101-19661231','19670101-19671231','19680101-19681231',&
    '19690101-19691231','19700101-19701231','19710101-19711231','19720101-19721231',&
    '19730101-19731231','19740101-19741231','19750101-19751231','19760101-19761231',&
    '19770101-19771231','19780101-19781231','19790101-19791231','19800101-19801231',&
    '18610101-19811231','19820101-19821231','19830101-19831231','19840101-19841231',&
    '19850101-19851231','19860101-19861231','19870101-19871231','19880101-19881231',&
    '19890101-19891231','19900101-19901231','19910101-19911231','19920101-19921231',&
    '19930101-19931231','19940101-19941231','19950101-19951231','19960101-19961231',&
    '19970101-19971231','19980101-19981231','19990101-19991231','20000101-20001231',&
    '20010101-20011231','20020101-20021231','20030101-20031231','20040101-20041231',&
    '20050101-20051231','20060101-20061231','20070101-20071231','20080101-20081231',&
    '20090101-20091231','20100101-20101231','20110101-20111231','20120101-20121231',&
    '20130101-20131231','20140101-20141231','20150101-20151231','20160101-20161231',&
    '20170101-20171231','20180101-20181231','20190101-20191231','20100101-20101231'/

    fill_value=1.0e20
    ! get the exact day
    call get_curr_date(yr, mon, day, tod)
    day=int(get_curr_calday())
    myyr=yr-1949 
    
    day1=int(get_curr_calday())+1
    if(day1.lt.366) then
        day1=day+1
        myyr1=myyr
    else
        day1=1
        myyr1=myyr+1
    endif

    ! get the timestep
    dt = get_step_size()
    fdda_nhtfrq=int(analysis_interval/dt)
    ii=mod(nstep,fdda_nhtfrq )

    nnstep=nstep*dt/analysis_interval

    allocate(camlev(pver))
    allocate(fdatao(pver))
   
    starts(1)=1

    startp(1)=1
    startp(2)=1
    scountp(3)=1

    start(1)=1
    start(2)=1
    start(3)=1
    scount(4)=1

    ! read data. At the end of the simulation, don't read data
    if(ii.eq.0.and.(.not.is_last_step())) then 
 
    ! t
    STATUS=NF_OPEN('/B3/cmip5/daily/MPI-ESM-LR/ta_day_MPI-ESM-LR_historical_r1i1p1_'//yearname(myyr)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'lat', LATID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LATID, LATLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'lon', LONID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LONID, LONLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'plev', LEVID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LEVID, LEVLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    if(.not.allocated(tmpfield1_3d)) allocate(tmpfield1_3d(lonlen,latlen,levlen,1))
    if(.not.allocated(tmpfield2_3d)) allocate(tmpfield2_3d(lonlen,latlen,levlen,1))
    if(.not.allocated(tmpfield3_3d)) allocate(tmpfield3_3d(lonlen,latlen,levlen,1))
    if(.not.allocated(tmpfield4_3d)) allocate(tmpfield4_3d(lonlen,latlen,levlen,1))

    scount(1)=lonlen
    scount(2)=latlen
    scount(3)=levlen

    ! allocate and remapping at the first time
    ! horizontal interpolation   
    mm=0
    do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols
     n=get_gcol_p(lchnk,i)
     mm=mm+1
    end do
    end do

    if(.not.allocated(u_ndg_old)) then

    allocate(u_ndg_old(1:mm,1:levlen))
    if(.not.allocated(v_ndg_old)) allocate(v_ndg_old(mm,1:levlen))
    if(.not.allocated(t_ndg_old)) allocate(t_ndg_old(mm,1:levlen))
    if(.not.allocated(q_ndg_old)) allocate(q_ndg_old(mm,1:levlen))
    if(.not.allocated(u_ndg_new)) allocate(u_ndg_new(mm,1:levlen))
    if(.not.allocated(v_ndg_new)) allocate(v_ndg_new(mm,1:levlen))
    if(.not.allocated(t_ndg_new)) allocate(t_ndg_new(mm,1:levlen))
    if(.not.allocated(q_ndg_new)) allocate(q_ndg_new(mm,1:levlen))

    if(.not.allocated(cmiplev)) allocate(cmiplev(levlen))
    if(.not.allocated(templat)) allocate(templat(latlen))
    if(.not.allocated(templon)) allocate(templon(lonlen))

    ! lat
    scounts(1)=latlen
    STATUS=NF_INQ_VARID(NCID,'lat',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,templat)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    ! lon
    scounts(1)=lonlen
    STATUS=NF_INQ_VARID(NCID,'lon',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,templon)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    ! lev
    scounts(1)=levlen
    STATUS=NF_INQ_VARID(NCID,'plev',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_DOUBLE(NCID,RHID,STARTS,SCOUNTS,cmiplev)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    cmiplev=cmiplev(levlen:1:-1)  ! reverse to the order from high to low
    
    call cmip5_to_cam_mapping(templon,templat,lonlen,latlen,phys_state)
    
    end if

    ! t
    STATUS=NF_INQ_VARID(NCID,'ta',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)
    call area_cmip5_3d(tmpfield1_3d,t_ndg_old,latlen,lonlen,levlen,mm,fill_value)
 
    STATUS=NF_OPEN('/B3/cmip5/daily//MPI-ESM-LR/ta_day_MPI-ESM-LR_historical_r1i1p1_'//yearname(myyr1)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'ta',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:) 
    call area_cmip5_3d(tmpfield1_3d,t_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    ! u
    STATUS=NF_OPEN('/B3/cmip5/daily//MPI-ESM-LR/ua_day_MPI-ESM-LR_historical_r1i1p1_'//yearname(myyr)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'ua',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:) 
    call area_cmip5_3d(tmpfield1_3d,u_ndg_old,latlen,lonlen,levlen,mm,fill_value)

    STATUS=NF_OPEN('/B3/cmip5/daily/MPI-ESM-LR/ua_day_MPI-ESM-LR_historical_r1i1p1_'//yearname(myyr1)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'ua',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)
    call area_cmip5_3d(tmpfield1_3d,u_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    ! v
    STATUS=NF_OPEN('/B3/cmip5/daily/MPI-ESM-LR/va_day_MPI-ESM-LR_historical_r1i1p1_'//yearname(myyr)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'va',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield2_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield2_3d=tmpfield2_3d(:,:,levlen:1:-1,:) 
    call area_cmip5_3d(tmpfield2_3d,v_ndg_old,latlen,lonlen,levlen,mm,fill_value)

    STATUS=NF_OPEN('/B3/cmip5/daily/MPI-ESM-LR/va_day_MPI-ESM-LR_historical_r1i1p1_'//yearname(myyr1)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'va',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield2_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield2_3d=tmpfield2_3d(:,:,levlen:1:-1,:)
    call area_cmip5_3d(tmpfield2_3d,v_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    ! q
    STATUS=NF_OPEN('/B3/cmip5/daily/MPI-ESM-LR/hur_day_MPI-ESM-LR_historical_r1i1p1_'//yearname(myyr)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'hur',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:) 
    call area_cmip5_3d(tmpfield1_3d,q_ndg_old,latlen,lonlen,levlen,mm,fill_value)

    STATUS=NF_OPEN('/B3/cmip5/daily/MPI-ESM-LR/hur_day_MPI-ESM-LR_historical_r1i1p1_'//yearname(myyr1)//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'hur',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,tmpfield1_3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)
    call area_cmip5_3d(tmpfield1_3d,q_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    do j=1,levlen
    do i=1,mm
    if(q_ndg_old(i,j).ne.fill_value.and.t_ndg_old(i,j).ne.fill_value) then
    q_ndg_old(i,j)=exp((t_ndg_old(i,j)-273.15)*17.67/(t_ndg_old(i,j)-29.65))* &
                   611.2*q_ndg_old(i,j)/100
    q_ndg_old(i,j)=0.622*q_ndg_old(i,j)/(cmiplev(j)-0.378*q_ndg_old(i,j))
    if(q_ndg_old(i,j).le.0) q_ndg_old(i,j)=1d-12
    else
    q_ndg_old(i,j)=fill_value
    endif
    
    if(q_ndg_new(i,j).ne.fill_value.and.t_ndg_new(i,j).ne.fill_value) then
    q_ndg_new(i,j)=exp((t_ndg_new(i,j)-273.15)*17.67/(t_ndg_new(i,j)-29.65))* &
                   611.2*q_ndg_new(i,j)/100
    q_ndg_new(i,j)=0.622*q_ndg_new(i,j)/(cmiplev(j)-0.378*q_ndg_new(i,j))
    if(q_ndg_new(i,j).le.0) q_ndg_new(i,j)=1d-12
    else
    q_ndg_new(i,j)=fill_value
    endif
    end do
    end do

   deallocate(tmpfield1_3d)
   deallocate(tmpfield2_3d)
   deallocate(tmpfield3_3d)
   deallocate(tmpfield4_3d)

   endif ! end of read data

   ! vertical interpolation
   mm=0
   do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols
     n=get_gcol_p(lchnk,i)
     mm=mm+1

     camlev(:)=phys_state(lchnk)%pmid(i,:)

     ! old state
     call vert_interpolation(u_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%u_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(u_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%u_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(v_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%v_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(v_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%v_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(t_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%t_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(t_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%t_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(q_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,4,vert_int,fill_value)
     phys_state(lchnk)%q_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(q_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,4,vert_int,fill_value)
     phys_state(lchnk)%q_ndg_new(i,:) = fdatao(:)

    end do
   end do

   if(is_last_step()) then
      if(allocated(cmiplev)) deallocate(cmiplev)
      if(allocated(u_ndg_old)) deallocate(u_ndg_old)
      if(allocated(v_ndg_old)) deallocate(v_ndg_old)
      if(allocated(t_ndg_old)) deallocate(t_ndg_old)
      if(allocated(q_ndg_old)) deallocate(q_ndg_old)
      if(allocated(u_ndg_new)) deallocate(u_ndg_new)
      if(allocated(v_ndg_new)) deallocate(v_ndg_new)
      if(allocated(t_ndg_new)) deallocate(t_ndg_new)
      if(allocated(q_ndg_new)) deallocate(q_ndg_new)
      if(allocated(templat)) deallocate(templat)
      if(allocated(templon)) deallocate(templon)
   end if

   deallocate(camlev)
   deallocate(fdatao)

   end subroutine cam_read_cmip5_mpiesm2m_daily_eul

!----------------------------------------------------------------------- 

 subroutine cam_read_ncep2( phys_state, nstep )
!----------------------------------------------------------------------- 
! Purpose: 
! Acquire and position the restart, master, primary and secondary
! datasets for a continuation run
! Author:Juanxiong He, 2012-06-01
!-----------------------------------------------------------------------
    use phys_grid, only: get_ncols_p, get_gcol_p, ngcols
    implicit none

    include 'netcdf.inc'

    type(physics_state), pointer :: phys_state(:)
    integer, intent(in) :: nstep

    integer :: lchnk, i, j, k, ii, m, n, csize, ncols, mm
    integer :: ntime, nnstep, fdda_nhtfrq
    real*8 :: dt
    real, dimension(:,:,:,:), allocatable :: tmpfield1_3d,tmpfield2_3d
    real, dimension(:,:,:), allocatable :: tmpfield1_2d,tmpfield2_2d
    real*8, dimension(:), allocatable :: camlev,fdatao
    real*8, dimension(:), allocatable :: templat,templon
    real, dimension(:), allocatable :: xfield1d
    real, dimension(:,:,:), allocatable :: xfield2d
    real, dimension(:,:,:,:), allocatable :: xfield3d
    
    integer, dimension(1:1) :: starts, scounts
    integer, dimension(1:3) :: startp, scountp
    integer, dimension(1:4) :: start, scount
    integer :: latid, lonid, levid, rhid, ncid, status
    real :: scale_factor, add_offset
    character*4 :: yearname, yearname1
 
    integer :: yr, mon, day, tod, day1

    fill_value=1.0e30

    ! get the exact day
    call get_curr_date(yr, mon, day, tod)
    day=(int(get_curr_calday())-1)*4+tod/21600+1
    yearname=char(yr/1000+48)//char(yr/100-yr/1000*10+48)//char(yr/10-yr/100*10+48)// &
             char(yr-yr/10*10+48)
    day1=day+1
    yearname1=yearname
    if(day.eq.1460) then
       day1=1
       yr=yr+1
       yearname1=char(yr/1000+48)//char(yr/100-yr/1000*10+48)//char(yr/10-yr/100*10+48)// &
                 char(yr-yr/10*10+48)
    end if
     
    ! get the timestep
    dt = get_step_size()
    fdda_nhtfrq=int(analysis_interval/dt)
    ii=mod(nstep,fdda_nhtfrq )

    nnstep=nstep*dt/analysis_interval
   
    allocate(camlev(pver))
    allocate(fdatao(pver))
   
    starts(1)=1

    startp(1)=1
    startp(2)=1
    scountp(3)=1

    start(1)=1
    start(2)=1
    start(3)=1
    scount(4)=1

    ! read data. At the end of the simulation, don't read data
    if(ii.eq.0.and.(.not.is_last_step())) then 
 
    ! air
    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/air.'//yearname//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'lat', LATID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LATID, LATLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'lon', LONID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LONID, LONLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'level', LEVID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LEVID, LEVLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    allocate(tmpfield1_3d(lonlen,latlen,levlen,1))
    allocate(tmpfield2_3d(lonlen,latlen,levlen,1))
    allocate(xfield3d(lonlen,latlen,levlen,1))

    scount(1)=lonlen
    scount(2)=latlen
    scount(3)=levlen

    ! allocate and remapping at the first time
    
    ! horizontal interpolation   
    mm=0
    do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols
     n=get_gcol_p(lchnk,i)
     mm=mm+1
    end do
    end do

    if(.not.allocated(u_ndg_old)) then
    
    allocate(u_ndg_old(1:mm,1:levlen))
    if(.not.allocated(v_ndg_old)) allocate(v_ndg_old(mm,1:levlen))
    if(.not.allocated(t_ndg_old)) allocate(t_ndg_old(mm,1:levlen))
    if(.not.allocated(q_ndg_old)) allocate(q_ndg_old(mm,1:levlen))
    if(.not.allocated(u_ndg_new)) allocate(u_ndg_new(mm,1:levlen))
    if(.not.allocated(v_ndg_new)) allocate(v_ndg_new(mm,1:levlen))
    if(.not.allocated(t_ndg_new)) allocate(t_ndg_new(mm,1:levlen))
    if(.not.allocated(q_ndg_new)) allocate(q_ndg_new(mm,1:levlen))
    
    if(.not.allocated(cmiplev)) allocate(cmiplev(levlen))
    if(.not.allocated(templat)) allocate(templat(latlen))
    if(.not.allocated(templon)) allocate(templon(lonlen))

    ! lat
    allocate(xfield1d(latlen))
    scounts(1)=latlen
    STATUS=NF_INQ_VARID(NCID,'lat',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTS,SCOUNTS,xfield1d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    templat=xfield1d*1.0_8
    deallocate(xfield1d)

    ! lon
    scounts(1)=lonlen
    allocate(xfield1d(lonlen))
    STATUS=NF_INQ_VARID(NCID,'lon',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTS,SCOUNTS,xfield1d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    templon=xfield1d*1.0_8
    deallocate(xfield1d)
   
    call cmip5_to_cam_mapping(templon,templat,lonlen,latlen,phys_state)
    deallocate(templon)
    deallocate(templat)
    
    end if

    ! lev
    scounts(1)=levlen
    allocate(xfield1d(levlen))
    STATUS=NF_INQ_VARID(NCID,'level',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTS,SCOUNTS,xfield1d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    cmiplev=xfield1d*100.0_8
    cmiplev=cmiplev(levlen:1:-1)  ! reverse to the order from small value to large value

    STATUS=NF_INQ_VARID(NCID,'air',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=xfield3d*scale_factor+add_offset
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)

    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/air.'//yearname1//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'air',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield2_3d=xfield3d*scale_factor+add_offset
    tmpfield2_3d=tmpfield2_3d(:,:,levlen:1:-1,:)

    call area_cmip5_3d(tmpfield1_3d,t_ndg_old,latlen,lonlen,levlen,mm,fill_value)
    call area_cmip5_3d(tmpfield2_3d,t_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    ! u
    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/uwnd.'//yearname//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'uwnd',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=xfield3d*scale_factor+add_offset
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)

    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/uwnd.'//yearname1//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'uwnd',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield2_3d=xfield3d*scale_factor+add_offset
    tmpfield2_3d=tmpfield2_3d(:,:,levlen:1:-1,:)

    call area_cmip5_3d(tmpfield1_3d,u_ndg_old,latlen,lonlen,levlen,mm,fill_value)
    call area_cmip5_3d(tmpfield2_3d,u_ndg_new,latlen,lonlen,levlen,mm,fill_value)
 
    ! v
    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/vwnd.'//yearname//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'vwnd',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=xfield3d*scale_factor+add_offset
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)

    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/vwnd.'//yearname1//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'vwnd',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield2_3d=xfield3d*scale_factor+add_offset
    tmpfield2_3d=tmpfield2_3d(:,:,levlen:1:-1,:)

    call area_cmip5_3d(tmpfield1_3d,v_ndg_old,latlen,lonlen,levlen,mm,fill_value)
    call area_cmip5_3d(tmpfield2_3d,v_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    ! q
    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/rhum.'//yearname//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'rhum',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=xfield3d*scale_factor+add_offset
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)

    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/rhum.'//yearname1//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'rhum',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield2_3d=xfield3d*scale_factor+add_offset
    tmpfield2_3d=tmpfield2_3d(:,:,levlen:1:-1,:)

    call area_cmip5_3d(tmpfield1_3d,q_ndg_old,latlen,lonlen,levlen,mm,fill_value)
    call area_cmip5_3d(tmpfield2_3d,q_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    do j=1,levlen
    do i=1,mm
    if(q_ndg_old(i,j).ne.fill_value.and.t_ndg_old(i,j).ne.fill_value) then
    q_ndg_old(i,j)=exp((t_ndg_old(i,j)-273.15)*17.67/(t_ndg_old(i,j)-29.65))* &
                   611.2*q_ndg_old(i,j)/100
    q_ndg_old(i,j)=0.622*q_ndg_old(i,j)/(cmiplev(j)-0.378*q_ndg_old(i,j))
    if(q_ndg_old(i,j).le.0) q_ndg_old(i,j)=1d-12
    else
    q_ndg_old(i,j)=fill_value
    endif

    if(q_ndg_new(i,j).ne.fill_value.and.t_ndg_new(i,j).ne.fill_value) then
    q_ndg_new(i,j)=exp((t_ndg_new(i,j)-273.15)*17.67/(t_ndg_new(i,j)-29.65))* &
                   611.2*q_ndg_new(i,j)/100
    q_ndg_new(i,j)=0.622*q_ndg_new(i,j)/(cmiplev(j)-0.378*q_ndg_new(i,j))
    if(q_ndg_new(i,j).le.0) q_ndg_new(i,j)=1d-12
    else
    q_ndg_new(i,j)=fill_value
    endif
    end do
    end do


   deallocate(tmpfield1_3d)
   deallocate(tmpfield2_3d)
   deallocate(xfield1d)
   deallocate(xfield3d)

   endif ! end of read data

   ! vertical interpolation
   mm=0
   do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols
     n=get_gcol_p(lchnk,i)
     mm=mm+1

     camlev(:)=phys_state(lchnk)%pmid(i,:)

     ! old state
     call vert_interpolation(u_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%u_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(v_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%v_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(t_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%t_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(q_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,4,vert_int,fill_value)
     phys_state(lchnk)%q_ndg_old(i,:) = fdatao(:)

     ! new state
     call vert_interpolation(u_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%u_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(v_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%v_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(t_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%t_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(q_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,4,vert_int,fill_value)
     phys_state(lchnk)%q_ndg_new(i,:) = fdatao(:)

    end do
   end do

   if(is_last_step()) then
      if(allocated(cmiplev)) deallocate(cmiplev)
      if(allocated(u_ndg_old)) deallocate(u_ndg_old)
      if(allocated(v_ndg_old)) deallocate(v_ndg_old)
      if(allocated(t_ndg_old)) deallocate(t_ndg_old)
      if(allocated(q_ndg_old)) deallocate(q_ndg_old)
      if(allocated(u_ndg_new)) deallocate(u_ndg_new)
      if(allocated(v_ndg_new)) deallocate(v_ndg_new)
      if(allocated(t_ndg_new)) deallocate(t_ndg_new)
      if(allocated(q_ndg_new)) deallocate(q_ndg_new)
   end if

   deallocate(camlev)
   deallocate(fdatao)

   end subroutine cam_read_ncep2

!----------------------------------------------------------------------
 subroutine cam_read_ncep2_daily( phys_state, nstep )
!----------------------------------------------------------------------- 
! Purpose: 
! Acquire and position the restart, master, primary and secondary
! datasets for a continuation run
! Author:Juanxiong He, 2012-06-01
!-----------------------------------------------------------------------
    use phys_grid, only: get_ncols_p, get_gcol_p, ngcols
    implicit none

    include 'netcdf.inc'

    type(physics_state), pointer :: phys_state(:)
    integer, intent(in) :: nstep

    integer :: lchnk, i, j, k, ii, m, n, csize, ncols, mm
    integer :: ntime, nnstep, fdda_nhtfrq
    real*8 :: dt
    real, dimension(:,:,:,:), allocatable :: tmpfield1_3d,tmpfield2_3d
    real, dimension(:,:,:), allocatable :: tmpfield1_2d,tmpfield2_2d
    real*8, dimension(:), allocatable :: camlev,fdatao
    real*8, dimension(:), allocatable :: templat,templon
    real, dimension(:), allocatable :: xfield1d
    real, dimension(:,:,:), allocatable :: xfield2d
    real, dimension(:,:,:,:), allocatable :: xfield3d

    integer, dimension(1:1) :: starts, scounts
    integer, dimension(1:3) :: startp, scountp
    integer, dimension(1:4) :: start, scount
    integer :: latid, lonid, levid, rhid, ncid, status
    real :: scale_factor, add_offset
    character*4 :: yearname, yearname1

    integer :: yr, mon, day, tod, day1

    fill_value=1.0e30


    ! get the exact day
    call get_curr_date(yr, mon, day, tod)
    day=int(get_curr_calday())
    yearname=char(yr/1000+48)//char(yr/100-yr/1000*10+48)//char(yr/10-yr/100*10+48)// &
             char(yr-yr/10*10+48)
    day1=day+1
    yearname1=yearname
    if(day.eq.365) then
       day1=1
       yr=yr+1
       yearname1=char(yr/1000+48)//char(yr/100-yr/1000*10+48)//char(yr/10-yr/100*10+48)// &
                 char(yr-yr/10*10+48)
    end if

    ! get the timestep
    dt = get_step_size()
    fdda_nhtfrq=int(analysis_interval/dt)
    ii=mod(nstep,fdda_nhtfrq )

    nnstep=nstep*dt/analysis_interval

    allocate(camlev(pver))
    allocate(fdatao(pver))

    starts(1)=1

    startp(1)=1
    startp(2)=1
    scountp(3)=1

    start(1)=1
    start(2)=1
    start(3)=1
    scount(4)=1

    ! read data. At the end of the simulation, don't read data
    if(ii.eq.0.and.(.not.is_last_step())) then

    ! air
    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/air.'//yearname//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'lat', LATID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LATID, LATLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'lon', LONID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LONID, LONLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'level', LEVID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LEVID, LEVLEN)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    allocate(tmpfield1_3d(lonlen,latlen,levlen,1))
    allocate(tmpfield2_3d(lonlen,latlen,levlen,1))
    allocate(xfield3d(lonlen,latlen,levlen,1))

    scount(1)=lonlen
    scount(2)=latlen
    scount(3)=levlen

    ! allocate and remapping at the first time

    ! horizontal interpolation   
    mm=0
    do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols
     n=get_gcol_p(lchnk,i)
     mm=mm+1
    end do
    end do

    if(.not.allocated(u_ndg_old)) then

    allocate(u_ndg_old(1:mm,1:levlen))
    if(.not.allocated(v_ndg_old)) allocate(v_ndg_old(mm,1:levlen))
    if(.not.allocated(t_ndg_old)) allocate(t_ndg_old(mm,1:levlen))
    if(.not.allocated(q_ndg_old)) allocate(q_ndg_old(mm,1:levlen))
    if(.not.allocated(u_ndg_new)) allocate(u_ndg_new(mm,1:levlen))
    if(.not.allocated(v_ndg_new)) allocate(v_ndg_new(mm,1:levlen))
    if(.not.allocated(t_ndg_new)) allocate(t_ndg_new(mm,1:levlen))
    if(.not.allocated(q_ndg_new)) allocate(q_ndg_new(mm,1:levlen))

    if(.not.allocated(cmiplev)) allocate(cmiplev(levlen))
    if(.not.allocated(templat)) allocate(templat(latlen))
    if(.not.allocated(templon)) allocate(templon(lonlen))

    ! lat
    allocate(xfield1d(latlen))
    scounts(1)=latlen
    STATUS=NF_INQ_VARID(NCID,'lat',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTS,SCOUNTS,xfield1d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    templat=xfield1d*1.0_8
    deallocate(xfield1d)

    ! lon
    scounts(1)=lonlen
    allocate(xfield1d(lonlen))
    STATUS=NF_INQ_VARID(NCID,'lon',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTS,SCOUNTS,xfield1d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    templon=xfield1d*1.0_8
    deallocate(xfield1d)

    call cmip5_to_cam_mapping(templon,templat,lonlen,latlen,phys_state)
    deallocate(templon)
    deallocate(templat)

    end if

    ! lev
    scounts(1)=levlen
    allocate(xfield1d(levlen))
    STATUS=NF_INQ_VARID(NCID,'level',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTS,SCOUNTS,xfield1d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    cmiplev=xfield1d*100.0_8
    cmiplev=cmiplev(levlen:1:-1)  ! reverse to the order from small value to large value

    STATUS=NF_INQ_VARID(NCID,'air',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=xfield3d*scale_factor+add_offset
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)

    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/air.'//yearname1//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'air',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield2_3d=xfield3d*scale_factor+add_offset
    tmpfield2_3d=tmpfield2_3d(:,:,levlen:1:-1,:)

    call area_cmip5_3d(tmpfield1_3d,t_ndg_old,latlen,lonlen,levlen,mm,fill_value)
    call area_cmip5_3d(tmpfield2_3d,t_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    ! u
    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/uwnd.'//yearname//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'uwnd',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=xfield3d*scale_factor+add_offset
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)

    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/uwnd.'//yearname1//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'uwnd',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield2_3d=xfield3d*scale_factor+add_offset
    tmpfield2_3d=tmpfield2_3d(:,:,levlen:1:-1,:)

    call area_cmip5_3d(tmpfield1_3d,u_ndg_old,latlen,lonlen,levlen,mm,fill_value)
    call area_cmip5_3d(tmpfield2_3d,u_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    ! v
    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/vwnd.'//yearname//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'vwnd',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=xfield3d*scale_factor+add_offset
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)

    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/vwnd.'//yearname1//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'vwnd',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield2_3d=xfield3d*scale_factor+add_offset
    tmpfield2_3d=tmpfield2_3d(:,:,levlen:1:-1,:)

    call area_cmip5_3d(tmpfield1_3d,v_ndg_old,latlen,lonlen,levlen,mm,fill_value)
    call area_cmip5_3d(tmpfield2_3d,v_ndg_new,latlen,lonlen,levlen,mm,fill_value)

    ! q
    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/rhum.'//yearname//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'rhum',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=xfield3d*scale_factor+add_offset
    tmpfield1_3d=tmpfield1_3d(:,:,levlen:1:-1,:)

    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/rhum.'//yearname1//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'rhum',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield2_3d=xfield3d*scale_factor+add_offset
    tmpfield2_3d=tmpfield2_3d(:,:,levlen:1:-1,:)

    call area_cmip5_3d(tmpfield1_3d,q_ndg_old,latlen,lonlen,levlen,mm,fill_value)
    call area_cmip5_3d(tmpfield2_3d,q_ndg_new,latlen,lonlen,levlen,mm,fill_value)


    do j=1,levlen
    do i=1,mm
    if(q_ndg_old(i,j).ne.fill_value.and.t_ndg_old(i,j).ne.fill_value) then
    q_ndg_old(i,j)=exp((t_ndg_old(i,j)-273.15)*17.67/(t_ndg_old(i,j)-29.65))* &
                   611.2*q_ndg_old(i,j)/100
    q_ndg_old(i,j)=0.622*q_ndg_old(i,j)/(cmiplev(j)-0.378*q_ndg_old(i,j))
    if(q_ndg_old(i,j).le.0) q_ndg_old(i,j)=1d-12
    else
    q_ndg_old(i,j)=fill_value
    endif

    if(q_ndg_new(i,j).ne.fill_value.and.t_ndg_new(i,j).ne.fill_value) then
    q_ndg_new(i,j)=exp((t_ndg_new(i,j)-273.15)*17.67/(t_ndg_new(i,j)-29.65))* &
                   611.2*q_ndg_new(i,j)/100
    q_ndg_new(i,j)=0.622*q_ndg_new(i,j)/(cmiplev(j)-0.378*q_ndg_new(i,j))
    if(q_ndg_new(i,j).le.0) q_ndg_new(i,j)=1d-12
    else
    q_ndg_new(i,j)=fill_value
    endif
    end do
    end do

   deallocate(tmpfield1_3d)
   deallocate(tmpfield2_3d)
   deallocate(xfield1d)
   deallocate(xfield3d)

   endif ! end of read data

   ! vertical interpolation
   mm=0
   do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols
     n=get_gcol_p(lchnk,i)
     mm=mm+1

     camlev(:)=phys_state(lchnk)%pmid(i,:)

     ! old state
     call vert_interpolation(u_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%u_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(v_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%v_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(t_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%t_ndg_old(i,:) = fdatao(:)
     call vert_interpolation(q_ndg_old(mm,:),cmiplev,fdatao,camlev,levlen,pver,4,vert_int,fill_value)
     phys_state(lchnk)%q_ndg_old(i,:) = fdatao(:)

     ! new state
     call vert_interpolation(u_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%u_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(v_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%v_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(t_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,1,vert_int,fill_value)
     phys_state(lchnk)%t_ndg_new(i,:) = fdatao(:)
     call vert_interpolation(q_ndg_new(mm,:),cmiplev,fdatao,camlev,levlen,pver,4,vert_int,fill_value)
     phys_state(lchnk)%q_ndg_new(i,:) = fdatao(:)

    end do
   end do

   if(is_last_step()) then
      if(allocated(cmiplev)) deallocate(cmiplev)
      if(allocated(u_ndg_old)) deallocate(u_ndg_old)
      if(allocated(v_ndg_old)) deallocate(v_ndg_old)
      if(allocated(t_ndg_old)) deallocate(t_ndg_old)
      if(allocated(q_ndg_old)) deallocate(q_ndg_old)
      if(allocated(u_ndg_new)) deallocate(u_ndg_new)
      if(allocated(v_ndg_new)) deallocate(v_ndg_new)
      if(allocated(t_ndg_new)) deallocate(t_ndg_new)
      if(allocated(q_ndg_new)) deallocate(q_ndg_new)
   end if

   deallocate(camlev)
   deallocate(fdatao)

   end subroutine cam_read_ncep2_daily
!----------------------------------------------------------------------
 
 subroutine cam_read_ncep2_sfc( cam_in, phys_state, nstep )
!----------------------------------------------------------------------- 
! Purpose: 
! Acquire and position the restart, master, primary and secondary
! datasets for a continuation run
! Author:Juanxiong He, 2012-06-01
!-----------------------------------------------------------------------
    use phys_grid, only: get_ncols_p, get_gcol_p, ngcols
    implicit none

    include 'netcdf.inc'

    type(cam_in_t) :: cam_in(begchunk:endchunk)
    type(physics_state), pointer :: phys_state(:)
    integer, intent(in) :: nstep

    integer :: lchnk, i, j, k, ii, m, n, csize, ncols, mm
    integer :: ntime, nnstep, fdda_nhtfrq
    real*8 :: dt,p0mb
    real, dimension(:,:,:,:), allocatable :: tmpfield1_3d,tmpfield2_3d
    real, dimension(:,:,:), allocatable :: tmpfield1_2d,tmpfield2_2d
    real*8, dimension(:), allocatable :: camlev,fdatao
    real*8, dimension(:), allocatable :: templat,templon
    real, dimension(:), allocatable :: xlat,xlon
    real, dimension(:,:,:), allocatable :: xfield2d
    real, dimension(:,:,:,:), allocatable :: xfield3d
    
    integer, dimension(1:1) :: starts, scounts
    integer, dimension(1:3) :: startp, scountp
    integer, dimension(1:4) :: start, scount
    integer :: latid, lonid, levid, rhid, ncid, status
    real  :: add_offset, scale_factor
    character*4 :: yearname, yearname1
 
    integer :: yr, mon, day, tod, day1
    real*8  :: coef, tfac
    real*8  :: val_analysis
    real*8  :: xtime, xtime_old, xtime_new

    fill_value=1.0e30

    ! get the exact day
    call get_curr_date(yr, mon, day, tod)
    day=(int(get_curr_calday())-1)*4+tod/21600+1
    yearname=char(yr/1000+48)//char(yr/100-yr/1000*10+48)//char(yr/10-yr/100*10+48)// &
             char(yr-yr/10*10+48)
    day1=day+1
    yearname1=yearname
    if(day.eq.1460) then
       day1=1
       yr=yr+1
       yearname1=char(yr/1000+48)//char(yr/100-yr/1000*10+48)//char(yr/10-yr/100*10+48)// &
                 char(yr-yr/10*10+48)
    end if
     

    ! get the timestep
    dt = get_step_size()
    fdda_nhtfrq=int(analysis_interval/dt)
    ii=mod(nstep,fdda_nhtfrq )

    nnstep=nstep*dt/analysis_interval
   
    xtime = dt*(nstep+1)
    xtime_old = FLOOR(xtime/analysis_interval) * analysis_interval * 1.0_r8
    xtime_new = xtime_old + analysis_interval * 1.0_r8
    coef = (xtime-xtime_old)/(xtime_new-xtime_old)

    starts(1)=1

    startp(1)=1
    startp(2)=1
    scountp(3)=1

    start(1)=1
    start(2)=1
    start(3)=1
    scount(4)=1

    ! read data. At the end of the simulation, don't read data
    if(ii.eq.0.and.(.not.is_last_step())) then 
 
    ! tsk
    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/skt.sfc.gauss.'//yearname//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'lat', LATID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LATID, LATLENsfc)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    STATUS = NF_INQ_DIMID(NCID, 'lon', LONID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS = NF_INQ_DIMLEN(NCID, LONID, LONLENsfc)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    allocate(tmpfield1_3d(lonlensfc,latlensfc,1,1))
    allocate(tmpfield2_3d(lonlensfc,latlensfc,1,1))
    allocate(tmpfield1_2d(lonlensfc,latlensfc,1))
    allocate(tmpfield2_2d(lonlensfc,latlensfc,1))
    allocate(xfield2d(lonlensfc,latlensfc,1))
    allocate(xfield3d(lonlensfc,latlensfc,1,1))

    ! allocate and remapping at the first time
    ! horizontal interpolation   
    mm=0
    do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols
     n=get_gcol_p(lchnk,i)
     mm=mm+1
    end do
    end do

    if(.not.allocated(u10_ndg_old)) then
    
    allocate(u10_ndg_old(1:mm,1:1))
    if(.not.allocated(u10_ndg_new)) allocate(u10_ndg_new(mm,1:1))
    if(.not.allocated(ts_ndg_old)) allocate(ts_ndg_old(mm,1:1))
    if(.not.allocated(ts_ndg_new)) allocate(ts_ndg_new(mm,1:1))
    if(.not.allocated(wsx_ndg_old)) allocate(wsx_ndg_old(mm,1:1))
    if(.not.allocated(wsx_ndg_new)) allocate(wsx_ndg_new(mm,1:1))
    if(.not.allocated(wsy_ndg_old)) allocate(wsy_ndg_old(mm,1:1))
    if(.not.allocated(wsy_ndg_new)) allocate(wsy_ndg_new(mm,1:1))
    if(.not.allocated(lhf_ndg_old)) allocate(lhf_ndg_old(mm,1:1))
    if(.not.allocated(lhf_ndg_new)) allocate(lhf_ndg_new(mm,1:1))
    if(.not.allocated(shf_ndg_old)) allocate(shf_ndg_old(mm,1:1))
    if(.not.allocated(shf_ndg_new)) allocate(shf_ndg_new(mm,1:1))
    if(.not.allocated(tref_ndg_old)) allocate(tref_ndg_old(mm,1:1))
    if(.not.allocated(tref_ndg_new)) allocate(tref_ndg_new(mm,1:1))
    if(.not.allocated(qref_ndg_old)) allocate(qref_ndg_old(mm,1:1))
    if(.not.allocated(qref_ndg_new)) allocate(qref_ndg_new(mm,1:1))
    if(.not.allocated(ulwrf_ndg_old)) allocate(ulwrf_ndg_old(mm,1:1))
    if(.not.allocated(ulwrf_ndg_new)) allocate(ulwrf_ndg_new(mm,1:1))

    if(.not.allocated(templat)) allocate(templat(latlensfc))
    if(.not.allocated(templon)) allocate(templon(lonlensfc))
    if(.not.allocated(xlat)) allocate(xlat(latlensfc))
    if(.not.allocated(xlon)) allocate(xlon(lonlensfc))

    ! lat
    scounts(1)=latlensfc
    STATUS=NF_INQ_VARID(NCID,'lat',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTS,SCOUNTS,xlat)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    ! lon
    scounts(1)=lonlensfc
    STATUS=NF_INQ_VARID(NCID,'lon',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTS,SCOUNTS,xlon)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    templon=xlon*1.0_8
    templat=xlat*1.0_8
    call cmip5_to_cam_mapping(templon,templat,lonlensfc,latlensfc,phys_state)
    deallocate(templon)
    deallocate(templat)
    deallocate(xlat)
    deallocate(xlon)

    cmip_weightsfc =  cmip_weight
    cmip_local_pointssfc =  cmip_local_points
    remap_cmipisfc =  remap_cmipi
    remap_cmipjsfc =  remap_cmipj

    end if

    scount(1)=lonlensfc
    scount(2)=latlensfc
    scount(3)=1

    scountp(1)=lonlensfc
    scountp(2)=latlensfc
    startp(3)=day
    STATUS=NF_INQ_VARID(NCID,'skt',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTP,SCOUNTP,xfield2d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_2d=xfield2d*scale_factor+add_offset

    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/skt.sfc.gauss.'//yearname1//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'skt',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    startp(3)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTP,SCOUNTP,xfield2d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield2_2d=xfield2d*scale_factor+add_offset

    call area_cmip5_2dsfc(tmpfield1_2d,ts_ndg_old,latlensfc,lonlensfc,1,mm,fill_value)
    call area_cmip5_2dsfc(tmpfield2_2d,ts_ndg_new,latlensfc,lonlensfc,1,mm,fill_value)

    ! wsx
    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/uflx.sfc.gauss.'//yearname//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'uflx',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    startp(3)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTp,SCOUNTp,xfield2d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_2d=xfield2d*scale_factor+add_offset

    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/uflx.sfc.gauss.'//yearname1//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'uflx',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    startp(3)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTp,SCOUNTp,xfield2d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield2_2d=xfield2d*scale_factor+add_offset

    call area_cmip5_2dsfc(tmpfield1_2d,wsx_ndg_old,latlensfc,lonlensfc,1,mm,fill_value)
    call area_cmip5_2dsfc(tmpfield2_2d,wsx_ndg_new,latlensfc,lonlensfc,1,mm,fill_value)

    ! wsy
    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/vflx.sfc.gauss.'//yearname//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'vflx',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    startp(3)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTp,SCOUNTp,xfield2d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_2d=xfield2d*scale_factor+add_offset

    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/vflx.sfc.gauss.'//yearname1//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'vflx',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    startp(3)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTp,SCOUNTp,xfield2d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    tmpfield2_2d=xfield2d*scale_factor+add_offset
    call area_cmip5_2dsfc(tmpfield1_2d,wsy_ndg_old,latlensfc,lonlensfc,1,mm,fill_value)
    call area_cmip5_2dsfc(tmpfield2_2d,wsy_ndg_new,latlensfc,lonlensfc,1,mm,fill_value)

    ! tref
    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/air.2m.gauss.'//yearname//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'air',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=xfield3d*scale_factor+add_offset

    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/air.2m.gauss.'//yearname1//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'air',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    tmpfield2_3d=xfield3d*scale_factor+add_offset
    call area_cmip5_3dsfc(tmpfield1_2d,tref_ndg_old,latlensfc,lonlensfc,1,mm,fill_value)
    call area_cmip5_3dsfc(tmpfield2_2d,tref_ndg_new,latlensfc,lonlensfc,1,mm,fill_value)

    ! qref
    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/shum.2m.gauss.'//yearname//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'shum',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=xfield3d*scale_factor+add_offset

    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/shum.2m.gauss.'//yearname1//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'shum',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    tmpfield2_3d=xfield3d*scale_factor+add_offset
    call area_cmip5_3dsfc(tmpfield1_2d,qref_ndg_old,latlensfc,lonlensfc,1,mm,fill_value)
    call area_cmip5_3dsfc(tmpfield2_2d,qref_ndg_new,latlensfc,lonlensfc,1,mm,fill_value)

    ! ulwrf
        STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/ulwrf.sfc.gauss.'//yearname//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'ulwrf',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    startp(3)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTp,SCOUNTp,xfield2d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_2d=xfield2d*scale_factor+add_offset

    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/ulwrf.sfc.gauss.'//yearname1//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'ulwrf',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    startp(3)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTp,SCOUNTp,xfield2d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    tmpfield2_2d=xfield2d*scale_factor+add_offset
    call area_cmip5_2dsfc(tmpfield1_2d,ulwrf_ndg_old,latlensfc,lonlensfc,1,mm,fill_value)
    call area_cmip5_2dsfc(tmpfield2_2d,ulwrf_ndg_new,latlensfc,lonlensfc,1,mm,fill_value)

    ! shf
    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/shtfl.sfc.gauss.'//yearname//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'shtfl',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    startp(3)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTp,SCOUNTp,xfield2d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_2d=xfield2d*scale_factor+add_offset

    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/shtfl.sfc.gauss.'//yearname1//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'shtfl',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    startp(3)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTp,SCOUNTp,xfield2d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    tmpfield2_2d=xfield2d*scale_factor+add_offset
    call area_cmip5_2dsfc(tmpfield1_2d,shf_ndg_old,latlensfc,lonlensfc,1,mm,fill_value)
    call area_cmip5_2dsfc(tmpfield2_2d,shf_ndg_new,latlensfc,lonlensfc,1,mm,fill_value)

    ! lhf
    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/lhtfl.sfc.gauss.'//yearname//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'lhtfl',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    startp(3)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTp,SCOUNTp,xfield2d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_2d=xfield2d*scale_factor+add_offset

    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/lhtfl.sfc.gauss.'//yearname1//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'lhtfl',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    startp(3)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,STARTp,SCOUNTp,xfield2d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    tmpfield2_2d=xfield2d*scale_factor+add_offset
    call area_cmip5_2dsfc(tmpfield1_2d,lhf_ndg_old,latlensfc,lonlensfc,1,mm,fill_value)
    call area_cmip5_2dsfc(tmpfield2_2d,lhf_ndg_new,latlensfc,lonlensfc,1,mm,fill_value)

    ! u10
    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/uwnd.10m.gauss.'//yearname//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'uwnd',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=xfield3d*scale_factor+add_offset

    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/vwnd.10m.gauss.'//yearname//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'vwnd',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    tmpfield2_3d=xfield3d*scale_factor+add_offset
    tmpfield1_3d=sqrt(tmpfield1_3d*tmpfield2_3d+tmpfield1_3d*tmpfield2_3d)
    call area_cmip5_3dsfc(tmpfield1_3d,u10_ndg_old,latlensfc,lonlensfc,1,mm,fill_value)

    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/uwnd.10m.gauss.'//yearname1//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'uwnd',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    tmpfield1_3d=xfield3d*scale_factor+add_offset

    STATUS=NF_OPEN('/R0/jhe/data/ncep2-data/vwnd.10m.gauss.'//yearname1//'.nc',NF_NOWRITE,NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_INQ_VARID(NCID,'vwnd',RHID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    start(4)=day1
    STATUS=NF_GET_VARA_REAL(NCID,RHID,START,SCOUNT,xfield3d)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'add_offset',add_offset)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_GET_ATT_REAL(NCID,RHID,'scale_factor',scale_factor)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)
    STATUS=NF_CLOSE(NCID)
    IF (STATUS .NE. NF_NOERR) write(*,*) NF_STRERROR(STATUS)

    tmpfield2_3d=xfield3d*scale_factor+add_offset
    tmpfield1_3d=sqrt(tmpfield1_3d*tmpfield2_3d+tmpfield1_3d*tmpfield2_3d)
    call area_cmip5_3dsfc(tmpfield1_3d,u10_ndg_new,latlensfc,lonlensfc,1,mm,fill_value)

   deallocate(tmpfield1_2d)
   deallocate(tmpfield2_2d)
   deallocate(tmpfield1_3d)
   deallocate(tmpfield2_3d)
   deallocate(xfield2d)
   deallocate(xfield3d)

   endif ! end of read data

   mm=0
   do lchnk=begchunk,endchunk
     ncols=get_ncols_p(lchnk)
    do i=1,ncols
     n=get_gcol_p(lchnk,i)
     mm=mm+1

     ! old state
     cam_in(lchnk)%u10_ndg_old(i) = u10_ndg_old(mm,1)
     cam_in(lchnk)%wsx_ndg_old(i) = wsx_ndg_old(mm,1)
     cam_in(lchnk)%wsy_ndg_old(i) = wsy_ndg_old(mm,1)
     cam_in(lchnk)%lhf_ndg_old(i) = lhf_ndg_old(mm,1)
     cam_in(lchnk)%shf_ndg_old(i) = shf_ndg_old(mm,1)
     cam_in(lchnk)%ts_ndg_old(i) = ts_ndg_old(mm,1)
     cam_in(lchnk)%qref_ndg_old(i) = qref_ndg_old(mm,1)
     cam_in(lchnk)%tref_ndg_old(i) = tref_ndg_old(mm,1)
     cam_in(lchnk)%ulwrf_ndg_old(i) = ulwrf_ndg_old(mm,1)

     ! new state
     cam_in(lchnk)%u10_ndg_new(i) = u10_ndg_new(mm,1)
     cam_in(lchnk)%wsx_ndg_new(i) = wsx_ndg_new(mm,1)
     cam_in(lchnk)%wsy_ndg_new(i) = wsy_ndg_new(mm,1)
     cam_in(lchnk)%lhf_ndg_new(i) = lhf_ndg_new(mm,1)
     cam_in(lchnk)%shf_ndg_new(i) = shf_ndg_new(mm,1)
     cam_in(lchnk)%ts_ndg_new(i) = ts_ndg_new(mm,1)
     cam_in(lchnk)%qref_ndg_new(i) = qref_ndg_new(mm,1)
     cam_in(lchnk)%tref_ndg_new(i) = tref_ndg_new(mm,1)
     cam_in(lchnk)%ulwrf_ndg_new(i) = ulwrf_ndg_new(mm,1)

     cam_in(lchnk)%u10(i) =  cam_in(lchnk)%u10_ndg_old(i) *( 1.0_r8 -coef ) +   &
                          cam_in(lchnk)%u10_ndg_new(i) * coef
     cam_in(lchnk)%wsx(i) =  cam_in(lchnk)%wsx_ndg_old(i) *( 1.0_r8 -coef ) +   &
                         cam_in(lchnk)%wsx_ndg_new(i) * coef
     cam_in(lchnk)%wsy(i) =  cam_in(lchnk)%wsy_ndg_old(i) *( 1.0_r8 -coef ) +   &
                         cam_in(lchnk)%wsy_ndg_new(i) * coef
     cam_in(lchnk)%lhf(i) =  cam_in(lchnk)%lhf_ndg_old(i) *( 1.0_r8 -coef ) +   &
                         cam_in(lchnk)%lhf_ndg_new(i) * coef
     cam_in(lchnk)%shf(i) =  cam_in(lchnk)%shf_ndg_old(i) *( 1.0_r8 -coef ) +   &
                         cam_in(lchnk)%shf_ndg_new(i) * coef
     cam_in(lchnk)%ts(i) =  cam_in(lchnk)%ts_ndg_old(i) *( 1.0_r8 -coef ) +   &
                         cam_in(lchnk)%ts_ndg_new(i) * coef
     cam_in(lchnk)%qref(i) =  cam_in(lchnk)%qref_ndg_old(i) *( 1.0_r8 -coef ) + &
                         cam_in(lchnk)%qref_ndg_new(i) * coef
     cam_in(lchnk)%tref(i) =  cam_in(lchnk)%tref_ndg_old(i) *( 1.0_r8 -coef ) + &
                         cam_in(lchnk)%tref_ndg_new(i) * coef
     cam_in(lchnk)%lwup(i) =  cam_in(lchnk)%ulwrf_ndg_old(i) *( 1.0_r8 -coef ) + &
                         cam_in(lchnk)%ulwrf_ndg_new(i) * coef
   end do
   end do

   if(is_last_step()) then
      if(allocated(u10_ndg_old)) deallocate(u10_ndg_old)
      if(allocated(u10_ndg_new)) deallocate(u10_ndg_new)
      if(allocated(tref_ndg_old)) deallocate(tref_ndg_old)
      if(allocated(tref_ndg_new)) deallocate(tref_ndg_new)
      if(allocated(qref_ndg_old)) deallocate(qref_ndg_old)
      if(allocated(qref_ndg_new)) deallocate(qref_ndg_new)
      if(allocated(ts_ndg_old)) deallocate(ts_ndg_old)
      if(allocated(ts_ndg_new)) deallocate(ts_ndg_new)
      if(allocated(lhf_ndg_old)) deallocate(lhf_ndg_old)
      if(allocated(lhf_ndg_new)) deallocate(lhf_ndg_new)
      if(allocated(shf_ndg_old)) deallocate(shf_ndg_old)
      if(allocated(shf_ndg_new)) deallocate(shf_ndg_new)
      if(allocated(wsx_ndg_old)) deallocate(wsx_ndg_old)
      if(allocated(wsx_ndg_new)) deallocate(wsx_ndg_new)
      if(allocated(wsy_ndg_old)) deallocate(wsy_ndg_old)
      if(allocated(wsy_ndg_new)) deallocate(wsy_ndg_new)
      if(allocated(ulwrf_ndg_old)) deallocate(ulwrf_ndg_old)
      if(allocated(ulwrf_ndg_new)) deallocate(ulwrf_ndg_new)
   end if

   end subroutine cam_read_ncep2_sfc
 
!----------------------------------------------------------------------
   subroutine cmip5_to_cam_mapping(templon,templat,lonlen,latlen,phys_state)
     use phys_grid, only: get_ncols_p, get_gcol_p, ngcols
     use pmgrid, only: plon, plev, plat
     use commap,    only: londeg, latdeg
     implicit none

     integer,intent(in) :: latlen, lonlen
     real*8,dimension(1:latlen),intent(inout) :: templat
     real*8,dimension(1:lonlen),intent(inout) :: templon
     type(physics_state), pointer :: phys_state(:)

     integer :: lchnk, ncols, mm, n
     integer :: i,j,k,ix,jy,km,local_index,global_index
     real*8 :: s,ss,sn,se,sw,sa
     real*8 :: dlondeg,dlatdeg
     real*8 :: pi=3.1415926_8

     real*8,dimension(:),allocatable::xlat,xlon
     real*8,dimension(:,:),allocatable:: ws,we,wn,ww
     real*8,dimension(:),allocatable:: cmipx,cmipy,cmipa,camarea

       ! CAM lat and lon on the chunk
       allocate(xlat(cmip_to_cam_points))
       allocate(xlon(cmip_to_cam_points))
       mm=0
       do lchnk=begchunk,endchunk
       ncols=get_ncols_p(lchnk)
       do i=1,ncols
       mm=mm+1
       n=get_gcol_p(lchnk,i)
       xlat(mm)=phys_state(lchnk)%lat(i)*180.0/pi
       xlon(mm)=phys_state(lchnk)%lon(i)*180.0/pi
       end do
       end do

       ! cmip5
       allocate(ww(lonlen,latlen))
       allocate(we(lonlen,latlen))
       allocate(ws(lonlen,latlen))
       allocate(wn(lonlen,latlen))
       allocate(cmipx(cmip_to_cam_points))
       allocate(cmipy(cmip_to_cam_points))
       allocate(cmipa(cmip_to_cam_points))
       allocate(camarea(cmip_to_cam_points))
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

       ! the CMIP5 in CAM domain      
       cmip_local_points=0
       global_index=0
       do jy=1,mm

         ! CAM latitude vertices, please note the parallell can be symestry
         ! xlat(jy-1) doesn't means it is the neighbor of xlat(jy)
         ! it need use latdeg
         if ((xlat(jy)-latdeg(1)).lt.0.1_8.and. &
              (xlat(jy)-latdeg(1)).gt.-0.1_8 ) then
            ss = -90.0_8
            sn = (xlat(jy)+ latdeg(2))/2
            km=1     
         elseif((xlat(jy)-latdeg(plat)).lt.0.1_8.and. &
                (xlat(jy)-latdeg(plat)).gt.-0.1_8 ) then
            ss =  (xlat(jy)+ latdeg(plat-1))/2
            sn = 90.0_8
            km=plat
         else
            do ix=1,plat
             if((xlat(jy)-latdeg(ix)).lt.0.1_8.and. &
                (xlat(jy)-latdeg(ix)).gt.-0.1_8) then
                 ss=(xlat(jy)+latdeg(ix))/2
                 sn=(xlat(jy)+latdeg(ix+1))/2
                 km=ix
             end if
            end do
         endif

         ! CAM longitude vertice
         if((xlon(jy)-londeg(1,km)).lt.0.1_8.and. &
            (xlon(jy)-londeg(1,km)).gt.-0.1_8) then
             se=(xlon(jy)+londeg(2,km))/2+360
             sw=(londeg(plon,km)+360.0_8)/2

         elseif((xlon(jy)-londeg(plon,km)).lt.0.1_8.and. &
                (xlon(jy)-londeg(plon,km)).gt.-0.1_8) then
             se=(xlon(jy)+360.0_8)/2
             sw=(xlon(jy)+londeg(plon-1,km))/2

         else
            do ix=1,plon
             if((xlon(jy)-londeg(ix,km)).lt.0.1_8.and. &
                (xlon(jy)-londeg(ix,km)).gt.-0.1_8) then
             se=(xlon(jy)+londeg(ix+1,km))/2
             sw=(xlon(jy)+londeg(ix-1,km))/2    
             end if
            end do
         endif

          ! CAM area
          s=sin(sn*pi/180)-sin(ss*pi/180)  
          camarea(jy)=s*(se-sw)

           ! cmip5 
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

              s=abs(sa-camarea(jy))/camarea(jy)
              if(s.gt.0.001_8) then
               print *, s,local_index
               print *, ss,sn,sw,se,camarea(jy)
              end if
              cmip_local_points(jy)=local_index
              do i=1,local_index
              if(s.gt.0.001_8) then
               print *, ws(cmipx(i),cmipy(i)),wn(cmipx(i),cmipy(i)),ww(cmipx(i),cmipy(i)),we(cmipx(i),cmipy(i)),cmipa(i)
              end if 
              global_index=global_index+1
              remap_cmipi(global_index)=cmipx(i)
              remap_cmipj(global_index)=cmipy(i)
              cmip_weight(global_index)=cmipa(i)/camarea(jy)
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
         deallocate(xlat)
         deallocate(xlon)
  end subroutine cmip5_to_cam_mapping

!-------------------------------------------------------------------------------
   subroutine area_cmip5_2dsfc(input,output,latlen,lonlen,levlen,mm,fill_value)
   implicit none
   integer,intent(in) :: latlen,lonlen,levlen,mm
   real, dimension(1:lonlen,1:latlen,1:levlen),intent(in) :: input
   real*8, dimension(1:mm,1:levlen),intent(inout) :: output
   real, intent(in) :: fill_value

   integer :: global_index, local_index, km, jy, ix, i, j
   real :: s
  
    output=0.0_8
    global_index=0
     do jy=1,mm
      local_index =  cmip_local_pointssfc(jy)
       s=0
      do ix=1,local_index
       global_index=global_index+1
       i = remap_cmipisfc(global_index)
       j = remap_cmipjsfc(global_index)
       if(input(i,j,1).ne.fill_value) then
        output(jy,1)=output(jy,1)+cmip_weightsfc(global_index)*input(i,j,1)
        s=s+cmip_weightsfc(global_index)
       endif
      end do
       if(s.gt.0) then
         output(jy,1)=output(jy,1)/s
       else
         output(jy,1)=fill_value
       endif
     end do
   
   end subroutine area_cmip5_2dsfc
!-------------------------------------------------------------------------------
   subroutine area_cmip5_3dsfc(input,output,latlen,lonlen,levlen,mm,fill_value)
   implicit none
   integer,intent(in) :: latlen,lonlen,levlen,mm
   real, dimension(1:lonlen,1:latlen,1:levlen,1:1),intent(in) :: input
   real*8, dimension(1:mm,1:levlen),intent(inout) :: output
   real,intent(in)::fill_value

   integer :: global_index, local_index, km, jy, ix, i, j
   real,dimension(1:levlen) :: s

    output=0.0_8
    global_index=0
     do jy=1,mm
      local_index =  cmip_local_pointssfc(jy)
       s=0.0
      do ix=1,local_index
       global_index=global_index+1
       i = remap_cmipisfc(global_index)
       j = remap_cmipjsfc(global_index) 
      do km=1,levlen
       if(input(i,j,km,1).ne.fill_value) then
       output(jy,km)=output(jy,km)+cmip_weightsfc(global_index)*input(i,j,km,1)
       s(km)=s(km)+cmip_weightsfc(global_index)
       endif
      end do
      end do
      do km=1,levlen
       if(s(km).ne.0) then
         output(jy,km)=output(jy,km)/s(km)
       else
         output(jy,km)=fill_value
       endif
      end do
     end do

   end subroutine area_cmip5_3dsfc
!-------------------------------------------------------------------------------
   subroutine area_cmip5_2d(input,output,latlen,lonlen,levlen,mm,fill_value)
   implicit none
   integer,intent(in) :: latlen,lonlen,levlen,mm
   real, dimension(1:lonlen,1:latlen,1:levlen),intent(in) :: input
   real*8, dimension(1:mm,1:levlen),intent(inout) :: output
   real, intent(in) :: fill_value

   integer :: global_index, local_index, km, jy, ix, i, j
   real :: s
  
    output=0.0_8
    global_index=0
     do jy=1,mm
      local_index =  cmip_local_points(jy)
       s=0
      do ix=1,local_index
       global_index=global_index+1
       i = remap_cmipi(global_index)
       j = remap_cmipj(global_index)
       if(input(i,j,1).ne.fill_value) then
        output(jy,1)=output(jy,1)+cmip_weight(global_index)*input(i,j,1)
        s=s+cmip_weight(global_index)
       endif
      end do
       if(s.gt.0) then
         output(jy,1)=output(jy,1)/s
       else
         output(jy,1)=fill_value
       endif
     end do
   
   end subroutine area_cmip5_2d
!-------------------------------------------------------------------------------
   subroutine area_cmip5_3d(input,output,latlen,lonlen,levlen,mm,fill_value)
   implicit none
   integer,intent(in) :: latlen,lonlen,levlen,mm
   real, dimension(1:lonlen,1:latlen,1:levlen,1:1),intent(in) :: input
   real*8, dimension(1:mm,1:levlen),intent(inout) :: output
   real,intent(in)::fill_value

   integer :: global_index, local_index, km, jy, ix, i, j
   real,dimension(1:levlen) :: s

    output=0.0_8
    global_index=0
     do jy=1,mm
      local_index =  cmip_local_points(jy)
       s=0.0
      do ix=1,local_index
       global_index=global_index+1
       i = remap_cmipi(global_index)
       j = remap_cmipj(global_index) 
      do km=1,levlen
       if(input(i,j,km,1).ne.fill_value) then
       output(jy,km)=output(jy,km)+cmip_weight(global_index)*input(i,j,km,1)
       s(km)=s(km)+cmip_weight(global_index)
       endif
      end do
      end do
      do km=1,levlen
       if(s(km).ne.0) then
         output(jy,km)=output(jy,km)/s(km)
       else
         output(jy,km)=fill_value
       endif
      end do
     end do

   end subroutine area_cmip5_3d
!----------------------------------------------------------------------
   subroutine vert_interpolation(fdatai,pi,fdatao,po,kde,kme,ch,vh,fill_value)
       integer kde, kme, ch, vh
       real*8 fdatai(kde),fdatao(kme)
       real*8 pi(kde),po(kme)
       real fill_value

       real*8 X(kde),Y(kde),Y2(kde),XINT,YINT,yp1,ypn
       integer i,j,k

       if(vh.eq.1)  call cinterpolation(fdatai,pi,fdatao,po,kde,kme,ch,fill_value)
       if(vh.eq.2)  call liinterpolation(fdatai,pi,fdatao,po,kde,kme,ch,fill_value)
       if(vh.eq.3)  call lainterpolation(fdatai,pi,fdatao,po,kde,kme,ch,fill_value)

       return
   end subroutine vert_interpolation
!----------------------------------------------------------------------
   subroutine cinterpolation(fdatai,pi,fdatao,po,kde,kme,ch,fill_value)
       integer kde, kme, ch
       real*8 fdatai(kde),fdatao(kme)
       real*8 pi(kde),po(kme)
       real fill_value

       real*8 X(kde),Y(kde),Y2(kde),XINT,YINT,yp1,ypn
       integer i,j,k,m
    
       yp1=2d30
       ypn=2d30
       m=0
        do k=1,kde
         if(fdatai(k).ne.fill_value) then
         m=m+1
         X(m)=dlog(pi(k))
         Y(m)=fdatai(k)
         endif
        enddo
        call cspline(X,Y,yp1,ypn,m,Y2)
        do k=1,kme
         XINT=dlog(po(k))
         call csplint(X,Y,Y2,m,XINT,YINT)
         fdatao(k)=YINT
         if(ch.eq.4) then
           if(YINT.le.0) then
             fdatao(k)=(1d-12)*1.0D0
           endif
         endif
        end do
    
       return
   end subroutine cinterpolation

   SUBROUTINE cspline(x,y,yp1,ypn,n,y2)
       INTEGER n
       real*8 x(n),y(n),y2(n),yp1,ypn

       INTEGER i,k
       real*8  p,qn,sig,un,u(500)

       if(yp1.gt..99d30) then
         y2(1)=0.0D0
         u(1)=0.0D0
       else
         y2(1)=-0.5D0
         u(1)=(3.0D0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
       endif

       do i=2,n-1
       sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
       p=sig*y2(i-1)+2.0D0
       y2(i)=(sig-1.)/p
       u(i)=(6.0D0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
         /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
       end do

       if(ypn.gt..99d30) then
          qn=0.0D0
          un=0.0D0
       else
          qn=0.5D0
          un=(3.0D0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
       endif

       y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0D0)
       do k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
       end do

       return
   end SUBROUTINE cspline

   SUBROUTINE csplint(xa,ya,y2a,n,x,y)
       INTEGER n
       real*8 xa(n),y2a(n),ya(n),x,y

       INTEGER  k,khi,klo
       real*8 a,b,h

       if (x.lt.xa(1)) then
          y=ya(1)
       elseif(x.gt.xa(n)) then
          y=ya(n)
       else
        klo=1
        khi=n
   1    if(khi-klo.gt.1) then
         k=(khi+klo)/2
         if(xa(k).gt.x)then
                 khi=k
         else
                 klo=k
         end if
         goto 1
        end if
        h=xa(khi)-xa(klo)
        if (h.eq.0.) then
         write(*,*) 'bad xa input in splint'
         stop
        end if
        a=(xa(khi)-x)/h
        b=(x-xa(klo))/h
        y=a*ya(klo)+b*ya(khi)+ &
        ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6D0
        end if

        return
   end SUBROUTINE csplint
!----------------------------------------------------------------------
   subroutine lainterpolation(fdatai,pi,fdatao,po,kde,kme,ch,fill_value)
       integer kde, kme, ch
       real*8 fdatai(kde),fdatao(kme)
       real*8 pi(kde),po(kme)
       real fill_value

       real*8 xi(kde),yi(kde),y(kde)
       integer i,j,k,m

      fdatao=0.0_8
      do k=1,kme
        y=1.0d0
        m=0
        do i=1,kde
         if(fdatai(i).ne.fill_value) then
            m=m+1
            xi(m)=pi(m)
            yi(m)=fdatai(i)
         endif
        end do
        do i=1,m
         do j=1,m
          if(i.ne.j) then
           y(i)=y(i)*(dlog(po(k))-dlog(xi(j)))/(dlog(xi(i))-dlog(xi(j)))
          end if
         end do
         fdatao(k)=fdatao(k)+y(i)*yi(i)
         if(ch.eq.4) then
           if(fdatao(k).le.0) then
             fdatao(k)=1d-12
           endif
         endif
        end do
      end do
      return

   end subroutine lainterpolation
!----------------------------------------------------------------------
   subroutine liinterpolation(fdatai,pi,fdatao,po,kde,kme,ch,fill_value)
       integer kde, kme, ch
       real*8 fdatai(kde),fdatao(kme)
       real*8 pi(kde),po(kme)
       real fill_value

       integer i,j,k

      do k=1,kme
         do j=1,kde-1
          if(po(k).le.pi(j+1).and.po(k).ge.pi(j)) then
           if(fdatai(j).ne.fill_value.and.fdatai(j+1).ne.fill_value) then
           fdatao(k)=fdatai(j+1)*(dlog(po(k))-dlog(pi(j)))/(dlog(pi(j+1))-dlog(pi(j)))+&
                     fdatai(j)*(dlog(pi(j+1))-dlog(po(k)))/(dlog(pi(j+1))-dlog(pi(j)))
           else
           fdatao(k)=-9999.0
           endif
          end if
         end do

         if(po(k).gt.pi(kde)) then
           if(fdatai(kde).ne.fill_value) then
           fdatao(k)=(fdatai(kde)-fdatai(kde-1))/(dlog(pi(kde))-dlog(pi(kde-1)))* &
                     (dlog(po(k))-dlog(pi(kde)))+fdatai(kde)
           else
           fdatao(k)=-9999.0
           endif
         end if

         if(po(k).lt.pi(1)) then
           if(fdatai(1).ne.fill_value) then
           fdatao(k)=(fdatai(2)-fdatai(1))/(dlog(pi(2))-dlog(pi(1)))*&
                     (dlog(po(k))-dlog(pi(1)))+fdatai(1)
           else
           fdatao(k)=-9999.0
           endif
         end if

         if(ch.eq.4) then
           if(fdatao(k).le.0) then
             fdatao(k)=1d-12
           endif
         endif
      end do
      return

   end subroutine liinterpolation   
end module cam_nudging
