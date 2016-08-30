
!*
!* all datasets reached here has the coordinate system as following:
!* -179.5W -> 179.5E     89.5N -> -89.5S
!* also we could inquire its longitude & latitude information from its nlon & nlat variables
!* others informatino all in var & land variables.
!* NetCDF files will be closed according to its time coordinate records.
!*

#include <define.h>

MODULE metdata

#if(defined GSWP2 || defined PRINCETON)

   use netcdf
   use spmd
   use precision

   implicit none

   save
   private

   real(r8), parameter :: spv = -9999.        ! special value
   integer,  parameter :: max_ntpr = 750      ! maximum number of time records to pre-read

 ! grids mapping method: 
 ! mapping_mode=0 for 1to1 extact mapping, mapping_maxg = 1
 ! mapping_mode=1 for coarse lnd grids, mapping_maxg = ?
 ! mapping_mode=2 for fine lnd grids, mapping_maxg = ?
 ! integer,  parameter :: mapping_mode = 1 
 ! integer,  parameter :: mapping_maxg = 16   ! maximum number of mapping grids
   integer,  parameter :: mapping_mode = 0 
   integer,  parameter :: mapping_maxg = 1    ! maximum number of mapping grids

   character(len=256) :: fname_qair  = "null"
   character(len=256) :: fname_tair  = "null"
   character(len=256) :: fname_dlwrf = "null"
   character(len=256) :: fname_dswrf = "null"
   character(len=256) :: fname_pres  = "null"
   character(len=256) :: fname_wind  = "null"
   character(len=256) :: fname_snow  = "null"
   character(len=256) :: fname_rain  = "null"
   character(len=256) :: fname_rainc = "null"

   character(len=256) :: vname_qair  = "null"
   character(len=256) :: vname_tair  = "null"
   character(len=256) :: vname_dlwrf = "null"
   character(len=256) :: vname_dswrf = "null"
   character(len=256) :: vname_pres  = "null"
   character(len=256) :: vname_wind  = "null"
   character(len=256) :: vname_snow  = "null"
   character(len=256) :: vname_rain  = "null"
   character(len=256) :: vname_rainc = "null"

   character(len=256) :: datapath

   integer :: fid_qair  = -1
   integer :: fid_tair  = -1
   integer :: fid_dlwrf = -1
   integer :: fid_dswrf = -1
   integer :: fid_pres  = -1
   integer :: fid_wind  = -1
   integer :: fid_snow  = -1
   integer :: fid_rain  = -1
   integer :: fid_rainc = -1

   integer :: vid_qair  = -1
   integer :: vid_tair  = -1
   integer :: vid_dlwrf = -1
   integer :: vid_dswrf = -1
   integer :: vid_pres  = -1
   integer :: vid_wind  = -1
   integer :: vid_snow  = -1
   integer :: vid_rain  = -1
   integer :: vid_rainc = -1

   integer itpr   ! index of current record in pre-read data
   integer ntpr   ! number of records currently pre-read

!* Pointers to current index of preread data buffer

   real(r4), pointer :: qair (:)
   real(r4), pointer :: tair (:)
   real(r4), pointer :: dlwrf(:)
   real(r4), pointer :: dswrf(:)
   real(r4), pointer :: pres (:)
   real(r4), pointer :: wind (:)
   real(r4), pointer :: snow (:)
   real(r4), pointer :: rain (:)
   real(r4), pointer :: rainc(:)

!* Pointers to maximum records of pre-read data buffer

   real(r4), pointer :: pr_qair (:,:)
   real(r4), pointer :: pr_tair (:,:)
   real(r4), pointer :: pr_dlwrf(:,:)
   real(r4), pointer :: pr_dswrf(:,:)
   real(r4), pointer :: pr_pres (:,:)
   real(r4), pointer :: pr_wind (:,:)
   real(r4), pointer :: pr_snow (:,:)
   real(r4), pointer :: pr_rain (:,:)
   real(r4), pointer :: pr_rainc(:,:)

   real(r4), pointer :: lon  (:)
   real(r4), pointer :: lat  (:)
   integer,  pointer :: land (:)

   integer,  pointer :: gridnum(:,:)
   integer,  pointer :: gridmap(:,:,:)
   real(r4), pointer :: gridwt(:,:,:)

   integer :: origin_year, origin_month, origin_day, origin_second
   integer :: nlon, nlat, nland, ntime, itime, deltim

   interface metdata_init
      module procedure metdata_init
   end interface

   interface metdata_read
      module procedure metdata_read
   end interface

   interface metdata_close
      module procedure metdata_close
   end interface

   public metdata_init, metdata_read, metdata_close

CONTAINS

   SUBROUTINE gen_fnames(year,month,day)
      implicit none

      integer, intent(in) :: year
      integer, intent(in) :: month
      integer, intent(in) :: day

      character(len=256)  :: date

#if defined (GSWP2)

      character(len=256), parameter :: qair_hint  = 'Qair_cru/Qair_cru'
      character(len=256), parameter :: tair_hint  = 'Tair_cru/Tair_cru'
      character(len=256), parameter :: dlwrf_hint = 'LWdown_srb/LWdown_srb'
      character(len=256), parameter :: dswrf_hint = 'SWdown_srb/SWdown_srb'
    ! character(len=256), parameter :: dlwrf_hint = 'LWdown_ncep/LWdown_ncep'
    ! character(len=256), parameter :: dswrf_hint = 'SWdown_ncep/SWdown_ncep'
    ! character(len=256), parameter :: dlwrf_hint = 'LWdown_era/LWdown_era'
    ! character(len=256), parameter :: dswrf_hint = 'SWdown_era/SWdown_era'
    ! character(len=256), parameter :: dlwrf_hint = 'LWdown_isccp/LWdown_isccp'
    ! character(len=256), parameter :: dswrf_hint = 'SWdown_isccp/SWdown_isccp'
      character(len=256), parameter :: pres_hint  = 'PSurf_ecor/PSurf_ecor'
      character(len=256), parameter :: wind_hint  = 'Wind_ncep/Wind_ncep'
      character(len=256), parameter :: snow_hint  = 'Snowf_gswp/Snowf_gswp'
      character(len=256), parameter :: rain_hint  = 'Rainf_gswp/Rainf_gswp'
      character(len=256), parameter :: rainc_hint = 'Rainf_C_gswp/Rainf_C_gswp'

#elif defined (PRINCETON)

      character(len=256), parameter :: qair_hint  = 'shum/shum_30min_'
      character(len=256), parameter :: tair_hint  = 'tas/tas_30min_'
      character(len=256), parameter :: dlwrf_hint = 'dlwrf/dlwrf_30min_'
      character(len=256), parameter :: dswrf_hint = 'dswrf/dswrf_30min_'
      character(len=256), parameter :: pres_hint  = 'pres/pres_30min_'
      character(len=256), parameter :: wind_hint  = 'wind/wind_30min_'
      character(len=256), parameter :: snow_hint  = 'null'
      character(len=256), parameter :: rain_hint  = 'prcp/prcp_30min_'
      character(len=256), parameter :: rainc_hint = 'null'

#endif

   !* ---------------------------------

#if defined (GSWP2)

      write(date,'(I4.4,I2.2)') year,month

      fname_qair  = trim(datapath)//trim(qair_hint) //trim(date)//'_30min.nc'
      fname_tair  = trim(datapath)//trim(tair_hint) //trim(date)//'_30min.nc'
      fname_dlwrf = trim(datapath)//trim(dlwrf_hint)//trim(date)//'_30min.nc'
      fname_dswrf = trim(datapath)//trim(dswrf_hint)//trim(date)//'_30min.nc'
      fname_pres  = trim(datapath)//trim(pres_hint) //trim(date)//'_30min.nc'
      fname_wind  = trim(datapath)//trim(wind_hint) //trim(date)//'_30min.nc'
      fname_snow  = trim(datapath)//trim(snow_hint) //trim(date)//'_30min.nc'
      fname_rain  = trim(datapath)//trim(rain_hint) //trim(date)//'_30min.nc'
      fname_rainc = trim(datapath)//trim(rainc_hint)//trim(date)//'_30min.nc'

      vname_qair  = 'Qair'
      vname_tair  = 'Tair'
      vname_dlwrf = 'LWdown'
      vname_dswrf = 'SWdown'
      vname_pres  = 'PSurf'
      vname_wind  = 'Wind'
      vname_snow  = 'Snowf'
      vname_rain  = 'Rainf'
      vname_rainc = 'Rainf_C'

#elif defined (PRINCETON)

      write(date,'(I4.4)') year

      fname_qair  = trim(datapath)//trim(qair_hint) //trim(date)//'-'//trim(date)//'.nc'
      fname_tair  = trim(datapath)//trim(tair_hint) //trim(date)//'-'//trim(date)//'.nc'
      fname_dlwrf = trim(datapath)//trim(dlwrf_hint)//trim(date)//'-'//trim(date)//'.nc'
      fname_dswrf = trim(datapath)//trim(dswrf_hint)//trim(date)//'-'//trim(date)//'.nc'
      fname_pres  = trim(datapath)//trim(pres_hint) //trim(date)//'-'//trim(date)//'.nc'
      fname_wind  = trim(datapath)//trim(wind_hint) //trim(date)//'-'//trim(date)//'.nc'
      fname_snow  = 'null'
      fname_rain  = trim(datapath)//trim(rain_hint) //trim(date)//'-'//trim(date)//'.nc'
      fname_rainc = 'null'

      vname_qair  = 'shum'
      vname_tair  = 'tas'
      vname_dlwrf = 'dlwrf'
      vname_dswrf = 'dswrf'
      vname_pres  = 'pres'
      vname_wind  = 'wind'
      vname_snow  = 'null'
      vname_rain  = 'prcp'
      vname_rainc = 'null'

#endif

   END SUBROUTINE gen_fnames

   SUBROUTINE metdata_init(lon_points,lat_points,dtime,path)

      implicit none

      integer, intent(in) :: lon_points  ! longitude points of land model
      integer, intent(in) :: lat_points  ! latitude points of land model
      real(r8),intent(in) :: dtime
      character(len=255), intent(in) :: path

      itime = 0
      ntime = 0

      itpr  = 1
      ntpr  = 0

      nland = -1
      deltim = int(dtime)

      datapath = adjustl(path)
      datapath = trim(datapath)//'/'

      allocate(gridnum(lon_points,lat_points))
      allocate(gridmap(mapping_maxg,lon_points,lat_points))
      allocate(gridwt(mapping_maxg,lon_points,lat_points))

   END SUBROUTINE metdata_init

   SUBROUTINE metdata_read(idate,lon_points,lat_points,gmask,lonw,latn,longxy,latixy,&
                          tair2d,qair2d,pres2d,rainc2d,rainl2d,&
                          windu2d,windv2d,dswrf2d,dlwrf2d,tair_z,qair_z,wind_z)

      implicit none

      integer, intent(in)  :: idate(3)
      integer, intent(in)  :: lon_points
      integer, intent(in)  :: lat_points
      integer, pointer :: gmask(:,:)
      real(r8),pointer :: lonw (:)
      real(r8),pointer :: latn (:)
      real(r8),pointer :: longxy (:,:)
      real(r8),pointer :: latixy (:,:)
      real(r8),pointer :: tair2d (:,:)
      real(r8),pointer :: qair2d (:,:)
      real(r8),pointer :: pres2d (:,:)
      real(r8),pointer :: rainc2d(:,:)
      real(r8),pointer :: rainl2d(:,:)
      real(r8),pointer :: windu2d(:,:)
      real(r8),pointer :: windv2d(:,:)
      real(r8),pointer :: dswrf2d(:,:)
      real(r8),pointer :: dlwrf2d(:,:)
      real(r8),pointer :: tair_z (:,:)
      real(r8),pointer :: qair_z (:,:)
      real(r8),pointer :: wind_z (:,:)

      real(r8) tmpbuf1(mapping_maxg), tmpbuf2(mapping_maxg)

      integer  i, j, k

   !* -----------------------------------------------

      if(itime.lt.1) then
         CALL open_files(idate)
         CALL grid_mapping(lon_points,lat_points,gmask,lonw,latn,longxy,latixy)
      endif

      CALL read_files

      tair2d (1:lon_points,1:lat_points) = spv
      qair2d (1:lon_points,1:lat_points) = spv
      pres2d (1:lon_points,1:lat_points) = spv
      rainc2d(1:lon_points,1:lat_points) = spv
      rainl2d(1:lon_points,1:lat_points) = spv
      windu2d(1:lon_points,1:lat_points) = spv
      windv2d(1:lon_points,1:lat_points) = spv
      dswrf2d(1:lon_points,1:lat_points) = spv
      dlwrf2d(1:lon_points,1:lat_points) = spv
      tair_z (1:lon_points,1:lat_points) = spv
      qair_z (1:lon_points,1:lat_points) = spv
      wind_z (1:lon_points,1:lat_points) = spv

      do j = 1, lat_points
      do i = 1, lon_points
         k = gridnum(i,j)

         if (k.eq.0) cycle

         if (gmask(i,j).ne.p_iam) stop 'gmask conflicts with gridmap'

         if (fid_rainc.gt.0) then
            tmpbuf1(1:k) = rain(gridmap(1:k,i,j))*gridwt(1:k,i,j)
            tmpbuf2(1:k) = rainc(gridmap(1:k,i,j))*gridwt(1:k,i,j)

            rainc2d(i,j) = sum(tmpbuf2(1:k))
            rainl2d(i,j) = max(0.,sum(tmpbuf1(1:k)-tmpbuf2(1:k)))
         else
            tmpbuf1(1:k) = rain(gridmap(1:k,i,j))*gridwt(1:k,i,j)
            rainc2d(i,j) = 0.
            rainl2d(i,j) = sum(tmpbuf1(1:k))
         endif

         if(rainl2d(i,j).lt.0) write(6,*), 'rainl2d(i,j).lt.0', gridmap(1:k,i,j),gridwt(1:k,i,j),rain(gridmap(1:k,i,j))

         if (fid_snow.gt.0) then
            tmpbuf1(1:k) = snow(gridmap(1:k,i,j))*gridwt(1:k,i,j)
            rainl2d(i,j) = rainl2d(i,j) + sum(tmpbuf1(1:k))
         endif

         tair2d (i,j) = sum(tair(gridmap(1:k,i,j))*gridwt(1:k,i,j))
         qair2d (i,j) = sum(qair(gridmap(1:k,i,j))*gridwt(1:k,i,j))
         pres2d (i,j) = sum(pres(gridmap(1:k,i,j))*gridwt(1:k,i,j))
         windu2d(i,j) = sum(wind(gridmap(1:k,i,j))*gridwt(1:k,i,j))
         windv2d(i,j) = 0.
         dswrf2d(i,j) = sum(dswrf(gridmap(1:k,i,j))*gridwt(1:k,i,j))
         dlwrf2d(i,j) = sum(dlwrf(gridmap(1:k,i,j))*gridwt(1:k,i,j))
         tair_z (i,j) = 50.
         qair_z (i,j) = 50.
         wind_z (i,j) = 50.
      end do 
      end do 

      if (itime.gt.ntime) then
         if(itpr.le.ntpr) then
            write(6,*), 'Error in pre-read data indexing', itpr, ntpr, itime, ntime
            call abort
         end if

         itime = 0
         ntpr  = 0

         CALL close_files
      end if

   END SUBROUTINE metdata_read

   SUBROUTINE metdata_close

      call close_files

      deallocate (gridnum)
      deallocate (gridmap)
      deallocate (gridwt)

   END SUBROUTINE metdata_close

   SUBROUTINE open_files(idate)
      implicit none

      integer, intent(in) :: idate(3)     ! current model time

      integer year, month, day, second, k, vid_t
      integer months(0:12)
      logical leapyear

   !* ---------------------------------

      year   = idate(1)
      day    = idate(2)
      second = idate(3)

      leapyear = (mod(year,4)==0.and.mod(year,100)/=0) .or. mod(year,400)==0

      if (leapyear) then
         months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
      else
         months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
      end if

      do k=1, 12
         if (day.le.months(k)) then
            month = k; exit
         end if
      end do

      day = day - months(k-1)

      if(p_master) then
         call gen_fnames(year,month,day)

         call nc_open(fname_qair ,vname_qair ,fid_qair ,vid_qair )
         call nc_open(fname_tair ,vname_tair ,fid_tair ,vid_tair )
         call nc_open(fname_dlwrf,vname_dlwrf,fid_dlwrf,vid_dlwrf)
         call nc_open(fname_dswrf,vname_dswrf,fid_dswrf,vid_dswrf)
         call nc_open(fname_pres ,vname_pres ,fid_pres ,vid_pres )
         call nc_open(fname_wind ,vname_wind ,fid_wind ,vid_wind )
         call nc_open(fname_snow ,vname_snow ,fid_snow ,vid_snow )
         call nc_open(fname_rain ,vname_rain ,fid_rain ,vid_rain )
         call nc_open(fname_rainc,vname_rainc,fid_rainc,vid_rainc)
      end if

#if(defined SPMD)
      call mpi_bcast(fid_qair ,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(fid_tair ,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(fid_dlwrf,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(fid_dswrf,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(fid_pres ,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(fid_wind ,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(fid_snow ,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(fid_rain ,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(fid_rainc,1,mpi_integer,0,p_comm,p_err)

      call mpi_bcast(nlon ,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(nlat ,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(nland,1,mpi_integer,0,p_comm,p_err)
#endif

   !* allocate memory for pre-read and forcing data grids lon&lat
      call var_alloc

   !* assume all files have the same time stamp
      if(p_master) then
         call sanity(nf90_inq_varid(fid_dswrf,'origin_year',vid_t))
         call sanity(nf90_get_var(fid_dswrf,vid_t,origin_year))

         call sanity(nf90_inq_varid(fid_dswrf,'origin_month',vid_t))
         call sanity(nf90_get_var(fid_dswrf,vid_t,origin_month))

         call sanity(nf90_inq_varid(fid_dswrf,'origin_day',vid_t))
         call sanity(nf90_get_var(fid_dswrf,vid_t,origin_day))

         call sanity(nf90_inq_varid(fid_dswrf,'origin_second',vid_t))
         call sanity(nf90_get_var(fid_dswrf,vid_t,origin_second))

         call sanity(nf90_inq_varid(fid_dswrf,'lon',vid_t))
         call sanity(nf90_get_var(fid_dswrf,vid_t,lon(:)))

         call sanity(nf90_inq_varid(fid_dswrf,'lat',vid_t))
         call sanity(nf90_get_var(fid_dswrf,vid_t,lat(:)))

         call sanity(nf90_inq_varid(fid_dswrf,'land',vid_t))
         call sanity(nf90_get_var(fid_dswrf,vid_t,land(:)))

         call sanity(nf90_inq_dimid(fid_dswrf,'time',vid_t))
         call sanity(nf90_inquire_dimension(fid_dswrf,vid_t,len=ntime))

      !* itime value....
      !* forcing data has some fields to indicate its original time, 
      !* so here set the right "itime" according to the model time and forcing data file time.

         write (*,"(A, X, I4.4 ':' I2.2 ':' I2.2 ':' I6.6)"), &
                  'data  yy:mm:dd:ss', origin_year, origin_month, origin_day, origin_second
         write (*,"(A, X, I4.4 ':' I2.2 ':' I2.2 ':' I6.6)"), &
                  'model yy:mm:dd:ss', year, month, day, second

         call date2sec(origin_year,origin_month,origin_day,origin_second, &
                       year,month,day,second,itime)

         itime = itime/deltim
      end if

#if(defined SPMD)
      call mpi_bcast(itime,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(ntime,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(lon,size(lon),mpi_real,0,p_comm,p_err)
      call mpi_bcast(lat,size(lat),mpi_real,0,p_comm,p_err)
      call mpi_bcast(land,size(land),mpi_integer,0,p_comm,p_err)
#endif

   END SUBROUTINE open_files

   SUBROUTINE read_files
      implicit none

      integer tidx1, tidx2

    ! write(*,*), 'itime', itime

      if (itpr .gt. ntpr) then
         itpr = 1
         ntpr = min(ntime-itime+1,max_ntpr)

         tidx1 = itime
         tidx2 = tidx1+ntpr-1

         if(p_master) then
            call nc_read(fid_qair ,vid_qair ,tidx1,tidx2,pr_qair )
            call nc_read(fid_tair ,vid_tair ,tidx1,tidx2,pr_tair )
            call nc_read(fid_dlwrf,vid_dlwrf,tidx1,tidx2,pr_dlwrf)
            call nc_read(fid_dswrf,vid_dswrf,tidx1,tidx2,pr_dswrf)
            call nc_read(fid_pres ,vid_pres ,tidx1,tidx2,pr_pres )
            call nc_read(fid_wind ,vid_wind ,tidx1,tidx2,pr_wind )
            call nc_read(fid_snow ,vid_snow ,tidx1,tidx2,pr_snow )
            call nc_read(fid_rain ,vid_rain ,tidx1,tidx2,pr_rain )
            call nc_read(fid_rainc,vid_rainc,tidx1,tidx2,pr_rainc)
         end if

#if(defined SPMD)
         if(fid_qair .gt.0) &
            call mpi_bcast(pr_qair ,size(pr_qair) ,mpi_real,0,p_comm,p_err)
         if(fid_tair .gt.0) &
            call mpi_bcast(pr_tair ,size(pr_tair) ,mpi_real,0,p_comm,p_err)
         if(fid_dlwrf.gt.0) &
            call mpi_bcast(pr_dlwrf,size(pr_dlwrf),mpi_real,0,p_comm,p_err)
         if(fid_dswrf.gt.0) &
            call mpi_bcast(pr_dswrf,size(pr_dswrf),mpi_real,0,p_comm,p_err)
         if(fid_pres .gt.0) &
            call mpi_bcast(pr_pres ,size(pr_pres) ,mpi_real,0,p_comm,p_err)
         if(fid_wind .gt.0) &
            call mpi_bcast(pr_wind ,size(pr_wind) ,mpi_real,0,p_comm,p_err)
         if(fid_snow .gt.0) &
            call mpi_bcast(pr_snow ,size(pr_snow) ,mpi_real,0,p_comm,p_err)
         if(fid_rain .gt.0) &
            call mpi_bcast(pr_rain ,size(pr_rain) ,mpi_real,0,p_comm,p_err)
         if(fid_rainc.gt.0) &
            call mpi_bcast(pr_rainc,size(pr_rainc),mpi_real,0,p_comm,p_err)
#endif
      end if

      if(fid_qair .gt.0) qair  => pr_qair (:,itpr)
      if(fid_tair .gt.0) tair  => pr_tair (:,itpr)
      if(fid_dlwrf.gt.0) dlwrf => pr_dlwrf(:,itpr)
      if(fid_dswrf.gt.0) dswrf => pr_dswrf(:,itpr)
      if(fid_pres .gt.0) pres  => pr_pres (:,itpr)
      if(fid_wind .gt.0) wind  => pr_wind (:,itpr)
      if(fid_snow .gt.0) snow  => pr_snow (:,itpr)
      if(fid_rain .gt.0) rain  => pr_rain (:,itpr)
      if(fid_rainc.gt.0) rainc => pr_rainc(:,itpr)

      itpr = itpr+1
      itime = itime+1

   END SUBROUTINE read_files

   SUBROUTINE close_files
      implicit none

      call var_dealloc
   
      if(p_master) then
         call nc_close(fid_qair, vid_qair )
         call nc_close(fid_tair, vid_tair )
         call nc_close(fid_dlwrf,vid_dlwrf)
         call nc_close(fid_dswrf,vid_dswrf)
         call nc_close(fid_pres, vid_pres )
         call nc_close(fid_wind, vid_wind )
         call nc_close(fid_snow, vid_snow )
         call nc_close(fid_rain, vid_rain )
         call nc_close(fid_rainc,vid_rainc)
      end if

   END SUBROUTINE close_files

   SUBROUTINE nc_open(fname,vname,fid,vid)
      implicit none

      character(len=256), intent(in) :: fname
      character(len=256), intent(in) :: vname
      integer, intent(out) :: fid
      integer, intent(out) :: vid

      integer vid_t, nlon_t, nlat_t, nland_t

   !* ---------------------------------

      if (fname .ne. "null") then
         print *, 'open file: '//trim(fname)

         call sanity(nf90_open(path=trim(fname),mode=nf90_nowrite,ncid=fid))
         call sanity(nf90_inq_varid(fid,trim(vname),vid))

      !* call sanity(nf90_inq_varid(fid,'nlon',vid_t))
      !* call sanity(nf90_get_var(fid,vid_t,nlon_t))
      !* call sanity(nf90_inq_varid(fid,'nlat',vid_t))
      !* call sanity(nf90_get_var(fid,vid_t,nlat_t))

         call sanity(nf90_inq_dimid(fid,'lon',vid_t))
         call sanity(nf90_inquire_dimension(fid,vid_t,len=nlon))
         call sanity(nf90_inq_dimid(fid,'lat',vid_t))
         call sanity(nf90_inquire_dimension(fid,vid_t,len=nlat))
         call sanity(nf90_inq_dimid(fid,'land',vid_t))
         call sanity(nf90_inquire_dimension(fid,vid_t,len=nland))

         print *, 'nlat of data', nlat
         print *, 'nlon of data', nlon
         print *, 'nland of data', nland
      end if

   END SUBROUTINE nc_open

   SUBROUTINE nc_close(fid,vid)
      implicit none

      integer, intent(inout) :: fid
      integer, intent(inout) :: vid

   !* ---------------------------------

      if (fid .gt. 0) then
         call sanity(nf90_close(fid))
         fid = -1
         vid = -1
      end if

   END SUBROUTINE nc_close

   SUBROUTINE nc_read(fid,vid,tidx1,tidx2,vbuf)
      implicit none

      integer,  intent(in) :: fid
      integer,  intent(in) :: vid
      integer,  intent(in) :: tidx1      ! begin time record to read
      integer,  intent(in) :: tidx2      ! end time record to read
      real(r4), pointer    :: vbuf(:,:)  ! pre-read buffer for time record between <tidx1:tidx2>

   !* ---------------------------------

      if (fid .gt. 0) then
         call sanity(nf90_get_var(fid,vid,vbuf,start=(/1,tidx1/),count=(/nland,tidx2-tidx1+1/)))
      end if
   END SUBROUTINE nc_read

   SUBROUTINE date2sec(year1,month1,day1,second1,year2,month2,day2,second2,nsec)
      implicit none

      integer, intent(in) :: year1
      integer, intent(in) :: month1     ! month of the year1
      integer, intent(in) :: day1       ! day of the month1
      integer, intent(in) :: second1    ! second of the day1
      integer, intent(in) :: year2
      integer, intent(in) :: month2     ! month of the year2
      integer, intent(in) :: day2       ! day of the month2
      integer, intent(in) :: second2    ! second of the day2
      integer, intent(out):: nsec       

   !* nsec = [year1:month1:day1:second1]=>[year2:month2:day2-1:86400] + second2

      integer months(0:12), nday, k
      logical leapyear

   !* -----------------------

      nday = 0
      nsec = 0

      do k=year1, year2
         leapyear = (mod(k,4)==0.and.mod(k,100)/=0) .or. mod(k,400)==0

         if (leapyear) then
            months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
         else
            months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
         end if

         if (year1.eq.year2) then
            nday = (months(month2-1)+(day2-1)) - (months(month1-1)+day1) + 1
         else if (k.eq.year1) then
            if (leapyear) then
               nday = nday + (366-(months(month1-1)+day1)+1)
            else
               nday = nday + (365-(months(month1-1)+day1)+1)
            end if
         else if (k.eq.year2) then
            nday = nday + (months(month2-1)+day2-1)
         else
            if (leapyear) then
               nday = nday + 366
            else
               nday = nday + 365
            end if
         end if
      end do

      nsec = nday*86400 - second1 + second2

   END SUBROUTINE date2sec

   SUBROUTINE grid_mapping(lon_points,lat_points,gmask,lonw,latn,longxy,latixy)
   !
   !--A simple grid mapping method --
   !
      implicit none

      integer,  intent(in) :: lon_points
      integer,  intent(in) :: lat_points
      integer,  pointer :: gmask(:,:)
      real(r8), pointer :: lonw(:)
      real(r8), pointer :: latn(:)
      real(r8), pointer :: longxy(:,:)
      real(r8), pointer :: latixy(:,:)

      integer, allocatable :: xy2land(:,:)

      integer i, j, k, i1, j1, i2, j2

      integer x1, x2, y1, y2, ngrd

      real(r8) dx1, dx2, dy1, dy2, a, minlon, maxlon

      allocate(xy2land(nlon,nlat))

      xy2land(:,:) = 0

      gridnum(:,:) = 0
      gridmap(:,:,:) = 0
      gridwt(:,:,:) = 0.

      do k = 1, nland
         i = mod(land(k),nlon)
         j = land(k)/nlon+1

         if(i.eq.0) then
            i = nlon
            j = j-1
         end if

         xy2land(i,j) = k  ! native 1d->2d mapping in forcing data
      end do

      IF(mapping_mode .eq. 0) THEN

         do j = 1, lat_points
         do i = 1, lon_points
            if(gmask(i,j).ne.p_iam) cycle

            gridmap(1,i,j) = xy2land(i,j)
          ! gridmap(1,i,j) = xy2land(i,nlat-j+1)  ! forcing data is in reverse order 
            gridwt(1,i,j) = 1.
            gridnum(i,j) = 1
         end do
         end do

      ELSE IF(mapping_mode .eq. 1) THEN

       ! Original longitude range of meteorology data is -180W~180E
       ! Convert them according to the grids of land model
         minlon = minval(lonw)
         maxlon = maxval(lonw)

         if(maxlon.gt.180.) then
            do i = 1, nlon
               if(lon(i).lt.minlon) lon(i) = lon(i)+360.
               if(lon(i).gt.maxlon) stop 'error in longitude convertion'
            end do
         end if

         do j1 = 1, lat_points
         do i1 = 1, lon_points
            if(gmask(i1,j1).ne.p_iam) cycle

            ngrd = 0

            do j2 = 1, nlat
               if(latn(j1).lt.lat(j2) .and. lat(j2).lt.latn(j1+1)) then  !* land grids in S->N ascending order
             ! if(latn(j1).gt.lat(j2) .and. lat(j2).gt.latn(j1+1)) then  !* land grids in N->S ascending order
                  do i2 = 1, nlon
                     if(lonw(i1).lt.lon(i2) .and. lon(i2).lt.lonw(i1+1)) then
                        if(xy2land(i2,j2).gt.0) then
                           ngrd = ngrd+1
                           gridmap(ngrd,i1,j1) = xy2land(i2,j2)
                        end if
                     end if
                  end do
               end if
            end do

            if(ngrd.eq.0)then
               write(6,*) "fatal error in mapping(mode=1, ngrd=0):", lonw(i1),lonw(i1+1),latn(j1),latn(j1+1)
               call abort
            else if(ngrd.gt.mapping_maxg) then
               write(6,*) "fatal error in mapping(mode=1, ngrd>mapping_maxg)",ngrd,mapping_maxg
               call abort
            end if

            if(ngrd.gt.0) then
               gridnum(i1,j1) = ngrd
               gridwt(1:ngrd,i1,j1) = 1.0/ngrd
            end if

         end do
         end do

      ELSE IF(mapping_mode .eq. 2) THEN

       ! Original longitude range of meteorology data is -180W~180E
       ! Convert grids of land model here...

       ! begin of convertion

       ! end of convertion

         do j1 = 1, lat_points
         do i1 = 1, lon_points
            if(gmask(i1,j1).ne.p_iam) cycle

            x1 = 0
            x2 = 0
            y1 = 0
            y2 = 0
   
            do i2 = 1, nlon-1
               if(longxy(i1,j1).ge.lon(i2) .and. longxy(i1,j1).lt.lon(i2+1)) then
                  x1 = i2
                  x2 = i2+1
                  exit
               end if
            end do
   
            do j2 = 1, nlat-1
               if(latixy(i1,j1).gt.lat(j2+1) .and. latixy(i1,j1).le.lat(j2)) then
                  y1 = j2
                  y2 = j2+1
                  exit
               end if
            end do
   
            if(x1.gt.0 .and. y1.gt.0) then
               dx1 = longxy(i1,j1)-lon(x1)
               dx2 = lon(x2)-longxy(i1,j1)
               dy1 = lat(y1)-latixy(i1,j1)
               dy2 = latixy(i1,j1)-lat(y2)
   
               ngrd = 0
               a = (lon(x2)-lon(x1))*(lat(y1)-lat(y2))

               if(xy2land(x1,y1).gt.0) then
                  ngrd = ngrd+1
                  gridmap(ngrd,i1,j1) = xy2land(x1,y1)
                  gridwt(ngrd,i1,j1) =  dx2*dy2/a
               end if
   
               if(xy2land(x2,y1).gt.0) then
                  ngrd = ngrd+1
                  gridmap(ngrd,i1,j1) = xy2land(x2,y1)
                  gridwt(ngrd,i1,j1) = dx1*dy2/a
               end if
   
               if(xy2land(x2,y2).gt.0) then
                  ngrd = ngrd+1
                  gridmap(ngrd,i1,j1) = xy2land(x2,y2)
                  gridwt(ngrd,i1,j1) = dx1*dy1/a
               end if
   
               if(xy2land(x1,y2).gt.0) then
                  ngrd = ngrd+1
                  gridmap(ngrd,i1,j1) = xy2land(x1,y2)
                  gridwt(ngrd,i1,j1) = dx2*dy1/a
               end if

               a = sum(gridwt(1:ngrd,i1,j1))
   
               if(a.lt.0 .or. a.gt.1) then
                  write(6,*) 'Fatal Error: no valid data points near the lnd grid:', longxy(i1,j1), latixy(i1,j1)
                  call abort
               end if
   
               gridnum(i1,j1) = ngrd
               gridwt(1:ngrd,i1,j1) = gridwt(1:ngrd,i1,j1)/a
            else
               write(6,*) 'Fatal Error: no data points near the lnd grid:', longxy(i1,j1), latixy(i1,j1)
               call abort
            end if
         end do
         end do

      END IF

      deallocate (xy2land)

   END SUBROUTINE grid_mapping

   SUBROUTINE var_alloc
      implicit none

      allocate (lon(nlon))
      allocate (lat(nlat))
      allocate (land(nland))

   !* Allocate memory for pre-read
      if(fid_qair .gt.0) allocate(pr_qair (nland,max_ntpr))
      if(fid_tair .gt.0) allocate(pr_tair (nland,max_ntpr))
      if(fid_dlwrf.gt.0) allocate(pr_dlwrf(nland,max_ntpr))
      if(fid_dswrf.gt.0) allocate(pr_dswrf(nland,max_ntpr))
      if(fid_pres .gt.0) allocate(pr_pres (nland,max_ntpr))
      if(fid_wind .gt.0) allocate(pr_wind (nland,max_ntpr))
      if(fid_snow .gt.0) allocate(pr_snow (nland,max_ntpr))
      if(fid_rain .gt.0) allocate(pr_rain (nland,max_ntpr))
      if(fid_rainc.gt.0) allocate(pr_rainc(nland,max_ntpr))

   ENDSUBROUTINE var_alloc

   SUBROUTINE var_dealloc
      implicit none

      if (associated(lon  )) deallocate(lon  )
      if (associated(lat  )) deallocate(lat  )
      if (associated(land )) deallocate(land )

      if (associated(pr_qair )) deallocate(pr_qair )
      if (associated(pr_tair )) deallocate(pr_tair )
      if (associated(pr_dlwrf)) deallocate(pr_dlwrf)
      if (associated(pr_dswrf)) deallocate(pr_dswrf)
      if (associated(pr_pres )) deallocate(pr_pres )
      if (associated(pr_wind )) deallocate(pr_wind )
      if (associated(pr_snow )) deallocate(pr_snow )
      if (associated(pr_rain )) deallocate(pr_rain )
      if (associated(pr_rainc)) deallocate(pr_rainc)

   END SUBROUTINE var_dealloc

   SUBROUTINE sanity(ret)
      implicit none
      integer, intent(in) :: ret

      if (ret .ne. nf90_noerr) then
         write(6,*) trim(nf90_strerror(ret)); stop
      end if
   END SUBROUTINE sanity

#endif

END MODULE metdata
