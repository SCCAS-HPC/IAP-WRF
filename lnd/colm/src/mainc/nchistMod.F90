#include <define.h>

module nchistMod

 use precision
 use paramodel
 use colm_varctl
 use netcdf
#ifdef RTM
 use RunoffMod
#endif

 implicit none

 private

 type history_info
    integer :: nac                               ! number of accumulation
    integer :: y_interval                        ! year interval
    integer :: m_interval                        ! month interval
    integer :: d_interval                        ! day interval
    integer :: h_interval                        ! hour interval
    logical :: is_valid                          ! valid history tape
    logical :: is_ready                          ! ready to write out
    logical :: is_newyear                        ! another model year coming
    logical :: fluxmask(nflux)                   ! logical flag to write LSM fluxes
    logical :: fldvmask(nfldv)                   ! logical flag to write LSM fluxes
    real(r8), pointer :: fldv(:,:)               ! output fluxes on LSM grid average
#ifdef RTM
  !*Currently ONLY support writing RTM fluxes to ONE archive
    logical :: with_rtm                          ! with RTM coordinates & fluxes to write
    logical :: fluxmask_rtm(nflux_rtm)           ! logical flag to write RTM fluxes
    logical :: fldvmask_rtm(nfldv_rtm)           ! logical flag to write RTM fluxes
    integer :: lnd_nac                           ! number of accumulation
    integer :: ocn_nac                           ! number of accumulation
    real(r8), pointer :: lnd_ave(:)              ! output fluxes 1
    real(r8), pointer :: ocn_ave(:)              ! output fluxes 2
#endif
#ifdef DGVM
    logical :: with_dgvm                         ! with DGVM coordinates & fluxes to write
    logical :: fluxmask_dgvm(nflux_dgvm)         ! logical flag to wirte DGVM fluxes
    logical :: fldvmask_dgvm(nfldv_dgvm)         ! logical flag to wirte DGVM fluxes
    integer :: nac_dgvm                          ! number of accumulation
    real(r8), pointer :: fldv_dgvm(:,:)          ! output fluxes on LSM grid average
#endif
    character(len=255) :: fhistory               ! file name of history
 end type history_info

 type(history_info) :: histArray(nMaxHist)

 integer            :: flxmap(nfldv)             ! 2D form fluxes -> 2/3D fluxes mapping
 integer            :: flxbeg(nflux)
 integer            :: flxdims(nflux)            ! dimensions of fluxes

 integer            :: flxvid(nflux)             ! netCDF variable id
 character(LEN=31)  :: flxname(nflux)            ! short name of flux variables
 character(LEN=31)  :: flxunit(nflux)            ! unit of flux variables
 character(LEN=255) :: flxinfo(nflux)            ! long name of flux variables

#ifdef RTM
 integer            :: flxdims_rtm(nflux_rtm)    ! dimensions of fluxes
 integer            :: flxvid_rtm(nflux_rtm)     ! netCDF variable id
 character(LEN=31)  :: flxname_rtm(nflux_rtm)    ! short name of flux variables
 character(LEN=31)  :: flxunit_rtm(nflux_rtm)    ! unit of flux variables
 character(LEN=255) :: flxinfo_rtm(nflux_rtm)    ! long name of flux variables
#endif

#ifdef DGVM
 integer            :: flxmap_dgvm(nfldv_dgvm)
 integer            :: flxbeg_dgvm(nflux_dgvm)
 integer            :: flxdims_dgvm(nflux_dgvm)  ! dimensions of fluxes

 integer            :: flxvid_dgvm(nflux_dgvm)   ! netCDF variable id
 character(LEN=31)  :: flxname_dgvm(nflux_dgvm)  ! short name of flux variables
 character(LEN=31)  :: flxunit_dgvm(nflux_dgvm)  ! unit of flux variables
 character(LEN=255) :: flxinfo_dgvm(nflux_dgvm)  ! long name of flux variables
#endif

 integer            :: fhistid

 interface nchist_init
    module procedure nchist_init
 end interface

 interface nchist_create
    module procedure nchist_create
 end interface

 interface nchist_addfield
    module procedure nchist_addfield
 end interface

 interface nchist_close
    module procedure nchist_close
 end interface

 interface nchist_exit
    module procedure nchist_exit
 end interface

 public nchist_init, nchist_create, nchist_addfield, nchist_close, nchist_exit
 public history_info, histArray

contains

   subroutine nchist_init

      use colm_varctl
      use paramodel, only: nfldv
      use colm_varMod, only: numgrid
      use colm_varMod, only: idx2_mrsos, idx2_fsno, idx2_tg, idx2_scv, idx2_rnof
      implicit none

      integer n, num_lnd, num_ocn

#ifdef RTM
      call get_proc_rof_global (num_lnd, num_ocn)
#endif

      do n = 1, nMaxHist
         histArray(n)%is_valid = .false.
         histArray(n)%is_ready = .false.
         histArray(n)%is_newyear = .false.
      end do

      n = 1

      if(lhist_yearly) then
         histArray(n)%nac = 0
         histArray(n)%y_interval = 1
         histArray(n)%m_interval = 0
         histArray(n)%d_interval = 0
         histArray(n)%h_interval = 0
         histArray(n)%fluxmask(:) = .true.

         allocate(histArray(n)%fldv(nfldv,numgrid))
         histArray(n)%fldv(:,:) = 0._r8

#ifdef RTM
         histArray(n)%with_rtm = .true.
         histArray(n)%fluxmask_rtm(:) = .true.

         histArray(n)%lnd_nac = 0
         histArray(n)%ocn_nac = 0
         allocate(histArray(n)%lnd_ave(num_lnd))
         allocate(histArray(n)%ocn_ave(num_ocn))
         histArray(n)%lnd_ave(:) = 0._r8
         histArray(n)%ocn_ave(:) = 0._r8
#endif
#ifdef DGVM
         histArray(n)%with_dgvm = .true.
         histArray(n)%fluxmask_dgvm(:) = .true.
         histArray(n)%nac_dgvm = 0
         allocate(histArray(n)%fldv_dgvm(nfldv_dgvm,numgrid))
         histArray(n)%fldv_dgvm(:,:) = 0._r8
#endif

         histArray(n)%is_valid = .true.
      end if

      n = n+1 

      if(lhist_monthly) then
         histArray(n)%nac = 0
         histArray(n)%y_interval = 0
         histArray(n)%m_interval = 1
         histArray(n)%d_interval = 0
         histArray(n)%h_interval = 0
         histArray(n)%fluxmask(:) = .true.

         allocate(histArray(n)%fldv(nfldv,numgrid))
         histArray(n)%fldv(:,:) = 0.

#ifdef RTM
         histArray(n)%with_rtm = .true.
         histArray(n)%fluxmask_rtm(:) = .true.

         histArray(n)%lnd_nac = 0
         histArray(n)%ocn_nac = 0
         allocate(histArray(n)%lnd_ave(num_lnd))
         allocate(histArray(n)%ocn_ave(num_ocn))
         histArray(n)%lnd_ave(:) = 0._r8
         histArray(n)%ocn_ave(:) = 0._r8
#endif
#ifdef DGVM
         histArray(n)%with_dgvm = .true.
         histArray(n)%fluxmask_dgvm(:) = .false.
         histArray(n)%fluxmask_dgvm(1:nflux_dgvm_accum) = .true.  ! enable accumulate
         histArray(n)%nac_dgvm = 0
         allocate(histArray(n)%fldv_dgvm(nfldv_dgvm,numgrid))
         histArray(n)%fldv_dgvm(:,:) = 0._r8
#endif

         histArray(n)%is_valid = .true.
      end if

      n = n+1 

      if(lhist_daily) then
         histArray(n)%nac = 0
         histArray(n)%y_interval = 0
         histArray(n)%m_interval = 0
         histArray(n)%d_interval = 1
         histArray(n)%h_interval = 0
         histArray(n)%fluxmask(:) = .false.
#ifdef CMIP
         histArray(n)%fluxmask(idx2_mrsos) = .true.   ! mrsos
         histArray(n)%fluxmask(idx2_fsno)  = .true.   ! fsno
         histArray(n)%fluxmask(idx2_tg)    = .true.   ! tg
         histArray(n)%fluxmask(idx2_scv)   = .true.   ! scv
         histArray(n)%fluxmask(idx2_rnof)  = .true.   ! rnof
#endif
     
         allocate(histArray(n)%fldv(nfldv,numgrid))
         histArray(n)%fldv(:,:) = 0.

#ifdef RTM
         histArray(n)%with_rtm = .false.
         histArray(n)%fluxmask_rtm(:) = .false.
#endif
#ifdef DGVM
         histArray(n)%with_dgvm = .false.
         histArray(n)%fluxmask_dgvm(:) = .false.
#endif

         histArray(n)%is_valid = .true.
      end if

      n = n+1

      if(lhist_3hourly) then
         histArray(n)%nac = 0
         histArray(n)%y_interval = 0
         histArray(n)%m_interval = 0
         histArray(n)%d_interval = 0
         histArray(n)%h_interval = 3
         histArray(n)%fluxmask(:) = .false.
#ifdef CMIP
         histArray(n)%fluxmask(idx2_mrsos) = .true.   ! mrsos
         histArray(n)%fluxmask(idx2_tg)    = .true.   ! tg
         histArray(n)%fluxmask(idx2_rnof)  = .true.   ! rnof
#endif

         allocate(histArray(n)%fldv(nfldv,numgrid))
         histArray(n)%fldv(:,:) = 0.

#ifdef RTM
         histArray(n)%with_rtm = .false.
         histArray(n)%fluxmask_rtm(:) = .false.
#endif
#ifdef DGVM
         histArray(n)%with_dgvm = .false.
         histArray(n)%fluxmask_dgvm(:) = .false.
#endif

         histArray(n)%is_valid = .true.
      end if

      call build_flx_info

   end subroutine nchist_init

   subroutine nchist_exit

      implicit none

      integer n

      do n = 1, nMaxHist
         if(histArray(n)%is_valid) then
            histArray(n)%is_valid = .false.
            deallocate(histArray(n)%fldv)
#ifdef RTM
            if(histArray(n)%with_rtm) then
               deallocate(histArray(n)%lnd_ave)
               deallocate(histArray(n)%ocn_ave)
            end if
#endif
#ifdef DGVM
            if(histArray(n)%with_dgvm) then
               deallocate(histArray(n)%fldv_dgvm)
            end if
#endif
         end if
      end do

   end subroutine nchist_exit

   subroutine nchist_create(idate,cdate,histinfo)

      use colm_varMod
#ifdef RTM
      use colm_rtmVar, only: rtmlon, rtmlat
      use RtmMod, only: longxy_r, latixy_r
#endif
      implicit none

      real(r4), parameter :: missing_value = -9999.

      integer, intent(in) :: idate(3)
      character(len=255), intent(in) :: cdate
      type(history_info), intent(in) :: histinfo

      real(r4), allocatable :: lon(:)
      real(r4), allocatable :: lat(:)
      real(r4), allocatable :: soilz(:)
      real(r4), allocatable :: time(:)

      integer dim_lon, dim_lat, dim_soilz, dim_time
      integer vid_lon, vid_lat, vid_soilz, vid_time
      integer vid_year, vid_day, vid_second
      integer vid_longxy, vid_latixy, vid_area, vid_landfrac, vid_landmask

#ifdef DGVM
      integer dim_pft, vid_pft
      integer,allocatable :: pft(:)
#endif

#ifdef RTM
      real(r4), allocatable :: lonrof(:)
      real(r4), allocatable :: latrof(:)
      integer dim_lonrof, dim_latrof
      integer vid_lonrof, vid_latrof
#endif

      integer ntime, L

      ntime = 1

      allocate (lon(lon_points))
      allocate (lat(lat_points))
      allocate (soilz(nl_soil))
      allocate (time(ntime))

      lon = real(longxy(:,1))
      lat = real(latixy(1,:))

      soilz = (/0.0071006, 0.0279250, 0.062258, 0.118865, 0.212193, &
                0.3660658, 0.6197585, 1.038027, 1.727635, 2.864607/)

      time = (/0/)

#ifdef RTM
      if(histinfo%with_rtm) then
         allocate (lonrof(rtmlon))
         allocate (latrof(rtmlat))

         lonrof = real(longxy_r(:,1))
         latrof = real(latixy_r(1,:))
      end if
#endif

#ifdef DGVM
      if(histinfo%with_dgvm) then
         allocate (pft(numpft_nat))
         pft(:) = (/(L,L=1,numpft_nat)/)
      end if
#endif

      call sanity(nf90_create(path=trim(histinfo%fhistory)//".nc",cmode=nf90_clobber,ncid=fhistid))

      call sanity(nf90_def_dim(fhistid,'lon',lon_points,dim_lon))
      call sanity(nf90_def_dim(fhistid,'lat',lat_points,dim_lat))
      call sanity(nf90_def_dim(fhistid,'soilz',nl_soil,dim_soilz))
      call sanity(nf90_def_dim(fhistid,'time',NF90_UNLIMITED,dim_time))
#ifdef RTM
      if(histinfo%with_rtm) then
         call sanity(nf90_def_dim(fhistid,'lonrof',rtmlon,dim_lonrof))
         call sanity(nf90_def_dim(fhistid,'latrof',rtmlat,dim_latrof))
      end if
#endif
#ifdef DGVM
      if(histinfo%with_dgvm) then 
         call sanity(nf90_def_dim(fhistid,'pft',numpft_nat,dim_pft))
      end if
#endif

      call sanity(nf90_def_var(fhistid,'origin_year',NF90_INT,varid=vid_year))
      call sanity(nf90_def_var(fhistid,'origin_day',NF90_INT,varid=vid_day))
      call sanity(nf90_def_var(fhistid,'origin_second',NF90_INT,varid=vid_second))

      call sanity(nf90_def_var(fhistid,'lon',NF90_FLOAT, &
                               dimids=(/dim_lon/),varid=vid_lon))

      call sanity(nf90_def_var(fhistid,'lat',NF90_FLOAT, &
                               dimids=(/dim_lat/),varid=vid_lat))

#ifdef RTM
      if(histinfo%with_rtm) then
         call sanity(nf90_def_var(fhistid,'lonrof',NF90_FLOAT, &
                                  dimids=(/dim_lonrof/),varid=vid_lonrof))

         call sanity(nf90_def_var(fhistid,'latrof',NF90_FLOAT, &
                                  dimids=(/dim_latrof/),varid=vid_latrof))
      end if
#endif

#ifdef DGVM
      if(histinfo%with_dgvm) then 
         call sanity(nf90_def_var(fhistid,'pft',NF90_INT, &
                                  dimids=(/dim_pft/),varid=vid_pft))
      end if
#endif

      call sanity(nf90_def_var(fhistid,'soilz',NF90_FLOAT,  &
                               dimids=(/dim_soilz/),varid=vid_soilz))

      call sanity(nf90_def_var(fhistid,'time',NF90_FLOAT, &
                               dimids=(/dim_time/),varid=vid_time))

      call sanity(nf90_def_var(fhistid,"longxy",NF90_FLOAT, &
                               dimids=(/dim_lon,dim_lat/),varid=vid_longxy))

      call sanity(nf90_def_var(fhistid,"latixy",NF90_FLOAT, &
                               dimids=(/dim_lon,dim_lat/),varid=vid_latixy))

      call sanity(nf90_def_var(fhistid,"area",NF90_FLOAT, &
                               dimids=(/dim_lon,dim_lat/),varid=vid_area))

      call sanity(nf90_def_var(fhistid,"landfrac",NF90_FLOAT, &
                               dimids=(/dim_lon,dim_lat/),varid=vid_landfrac))

      call sanity(nf90_def_var(fhistid,"landmask",NF90_INT, &
                               dimids=(/dim_lon,dim_lat/),varid=vid_landmask))

      do L = 1, nflux
         if(histinfo%fluxmask(L)) then
            if(flxdims(L).eq.2) then
               call sanity(nf90_def_var(fhistid,trim(flxname(L)),NF90_FLOAT,&
                                        dimids=(/dim_lon,dim_lat,dim_time/),varid=flxvid(L)))
            else if(flxdims(L).eq.3) then
               call sanity(nf90_def_var(fhistid,trim(flxname(L)),NF90_FLOAT,&
                                        dimids=(/dim_lon,dim_lat,dim_soilz,dim_time/),varid=flxvid(L)))
            end if
         end if
      end do

#ifdef RTM
      if(histinfo%with_rtm) then
         do L = 1, nflux_rtm
            if(histinfo%fluxmask_rtm(L)) then
               call sanity(nf90_def_var(fhistid,trim(flxname_rtm(L)),NF90_FLOAT,&
                                        dimids=(/dim_lonrof,dim_latrof,dim_time/),varid=flxvid_rtm(L)))
            end if
         end do
      end if
#endif

#ifdef DGVM
      if(histinfo%with_dgvm) then 
         do L = 1, nflux_dgvm
            if(histinfo%fluxmask_dgvm(L)) then
               if(flxdims_dgvm(L).eq.2) then
                  call sanity(nf90_def_var(fhistid,trim(flxname_dgvm(L)),NF90_FLOAT,&
                                           dimids=(/dim_lon,dim_lat,dim_time/),varid=flxvid_dgvm(L)))
               else if(flxdims_dgvm(L).eq.3) then
                  call sanity(nf90_def_var(fhistid,trim(flxname_dgvm(L)),NF90_FLOAT,&
                                           dimids=(/dim_lon,dim_lat,dim_pft,dim_time/),varid=flxvid_dgvm(L)))
               end if
            end if
         end do
      end if
#endif

      call sanity(nf90_put_att(fhistid,vid_lon,"long_name","coordinate longitude"))
      call sanity(nf90_put_att(fhistid,vid_lon,"units","degrees_east"))

      call sanity(nf90_put_att(fhistid,vid_lat,"long_name","coordinate latitude"))
      call sanity(nf90_put_att(fhistid,vid_lat,"units","degrees_north"))

#ifdef RTM
      if(histinfo%with_rtm) then
         call sanity(nf90_put_att(fhistid,vid_lonrof,"long_name","coordinate longitude of RTM"))
         call sanity(nf90_put_att(fhistid,vid_lonrof,"units","degrees_east"))

         call sanity(nf90_put_att(fhistid,vid_latrof,"long_name","coordinate latitude of RTM"))
         call sanity(nf90_put_att(fhistid,vid_latrof,"units","degrees_north"))
      end if
#endif

#ifdef DGVM
      if(histinfo%with_dgvm) then 
         call sanity(nf90_put_att(fhistid,vid_pft,"long_name","PFT categories excluding 2 crops and 1 soil types"))
      end if
#endif

      call sanity(nf90_put_att(fhistid,vid_soilz,"long_name","coordinate soil levels"))
      call sanity(nf90_put_att(fhistid,vid_soilz,"units","m"))

      call sanity(nf90_put_att(fhistid,vid_time,"long_name","time"))
      call sanity(nf90_put_att(fhistid,vid_time,"units","days since "//trim(cdate)))

      call sanity(nf90_put_att(fhistid,vid_longxy,"long_name","longitude"))
      call sanity(nf90_put_att(fhistid,vid_longxy,"units","degrees_east"))

      call sanity(nf90_put_att(fhistid,vid_latixy,"long_name","latitude"))
      call sanity(nf90_put_att(fhistid,vid_latixy,"units","degrees_north"))

      call sanity(nf90_put_att(fhistid,vid_area,"long_name","grid cell areas"))
      call sanity(nf90_put_att(fhistid,vid_area,"units","km^2"))

      call sanity(nf90_put_att(fhistid,vid_landfrac,"long_name","land fraction"))
      call sanity(nf90_put_att(fhistid,vid_landfrac,"units","%"))

      call sanity(nf90_put_att(fhistid,vid_landmask,"long_name","land/ocean mask"))

      do L = 1, nflux
         if(histinfo%fluxmask(L)) then
            call sanity(nf90_put_att(fhistid,flxvid(L),"long_name",trim(flxinfo(L))))
            call sanity(nf90_put_att(fhistid,flxvid(L),"units",trim(flxunit(L))))
            call sanity(nf90_put_att(fhistid,flxvid(L),"_FillValue",missing_value))
         end if
      end do

#ifdef RTM
      if(histinfo%with_rtm) then
         do L = 1, nflux_rtm
            if(histinfo%fluxmask_rtm(L)) then
               call sanity(nf90_put_att(fhistid,flxvid_rtm(L),"long_name",trim(flxinfo_rtm(L))))
               call sanity(nf90_put_att(fhistid,flxvid_rtm(L),"units",trim(flxunit_rtm(L))))
               call sanity(nf90_put_att(fhistid,flxvid_rtm(L),"_FillValue",missing_value))
            end if
         end do
      end if
#endif

#ifdef DGVM
      if(histinfo%with_dgvm) then
         do L = 1, nflux_dgvm
            if(histinfo%fluxmask_dgvm(L)) then
               call sanity(nf90_put_att(fhistid,flxvid_dgvm(L),"long_name",trim(flxinfo_dgvm(L))))
               call sanity(nf90_put_att(fhistid,flxvid_dgvm(L),"units",trim(flxunit_dgvm(L))))
               call sanity(nf90_put_att(fhistid,flxvid_dgvm(L),"_FillValue",missing_value))
            end if
         end do
      end if
#endif

      call sanity(nf90_enddef(fhistid))

      call sanity(nf90_put_var(fhistid,vid_lon,lon(:)))
      call sanity(nf90_put_var(fhistid,vid_lat,lat(:)))
#ifdef RTM
      if(histinfo%with_rtm) then
         call sanity(nf90_put_var(fhistid,vid_lonrof,lonrof(:)))
         call sanity(nf90_put_var(fhistid,vid_latrof,latrof(:)))
      end if
#endif
#ifdef DGVM
      if(histinfo%with_dgvm) then
         call sanity(nf90_put_var(fhistid,vid_pft,pft(:)))
      end if
#endif
      call sanity(nf90_put_var(fhistid,vid_time,time(:)))
      call sanity(nf90_put_var(fhistid,vid_soilz,soilz(:)))
      call sanity(nf90_put_var(fhistid,vid_year,(/idate(1)/)))
      call sanity(nf90_put_var(fhistid,vid_day,(/idate(2)/)))
      call sanity(nf90_put_var(fhistid,vid_second,(/idate(3)/)))
      call sanity(nf90_put_var(fhistid,vid_longxy,longxy(:,:)))
      call sanity(nf90_put_var(fhistid,vid_latixy,latixy(:,:)))
      call sanity(nf90_put_var(fhistid,vid_area,area(:,:)))
      call sanity(nf90_put_var(fhistid,vid_landfrac,landfrac(:,:)))
      call sanity(nf90_put_var(fhistid,vid_landmask,landmask(:,:)))

      deallocate (lon)
      deallocate (lat)
      deallocate (soilz)
      deallocate (time)

#ifdef RTM
      if(histinfo%with_rtm) then
         deallocate (lonrof)
         deallocate (latrof)
      end if
#endif

#ifdef DGVM
      if(histinfo%with_dgvm) then
         deallocate (pft)
      end if
#endif
      
   end subroutine nchist_create

   subroutine nchist_addfield(L,nlon,nlat,fldxy,gridtype)

      implicit none

      integer,  intent(in) :: L
      integer,  intent(in) :: nlon
      integer,  intent(in) :: nlat
      real(r4), intent(in) :: fldxy(nlon,nlat)
      character(len=*), intent(in) :: gridtype

      integer k, beg, ndim, lev, vid

      if(trim(gridtype).eq.'lsm') then
         k = flxmap(L)
         beg = flxbeg(k)
         ndim = flxdims(k)
         vid = flxvid(k)

         if(ndim.eq.2) then
            call sanity(nf90_put_var(fhistid,vid,fldxy(:,:), &
                                     start=(/1,1,1/),count=(/nlon,nlat,1/)))
         else if(ndim.eq.3) then
            lev = L-beg+1
            call sanity(nf90_put_var(fhistid,vid,fldxy(:,:), &
                                     start=(/1,1,lev,1/),count=(/nlon,nlat,1,1/)))
         else
            write(6,*) 'fatal error, wrong LSM dimension number', ndim
            call abort
         end if
      else if(trim(gridtype).eq.'dgvm') then
#ifdef DGVM
         k = flxmap_dgvm(L)
         beg = flxbeg_dgvm(k)
         ndim = flxdims_dgvm(k)
         vid = flxvid_dgvm(k)

       ! write(6,*), 'DGVM nchist_addfield', L, NLON, NLAT, vid, trim(flxname_dgvm(k)), minval(fldxy), maxval(fldxy)

         if(ndim.eq.2) then
            call sanity(nf90_put_var(fhistid,vid,fldxy(:,:), &
                                     start=(/1,1,1/),count=(/nlon,nlat,1/)))
         else if(ndim.eq.3) then
            lev = L-beg+1
            call sanity(nf90_put_var(fhistid,vid,fldxy(:,:), &
                                     start=(/1,1,lev,1/),count=(/nlon,nlat,1,1/)))
         else
            write(6,*) 'fatal error, wrong DGVM dimension number', ndim
            call abort
         end if
#endif
      else if(trim(gridtype).eq.'rtm') then
#ifdef RTM
         vid = flxvid_rtm(L)

         call sanity(nf90_put_var(fhistid,vid,fldxy(:,:), &
                                  start=(/1,1,1/),count=(/nlon,nlat,1/)))
#endif
      end if

   end subroutine nchist_addfield

   subroutine nchist_close

      implicit none

      call sanity(nf90_close(fhistid))

   end subroutine nchist_close

   subroutine sanity(ret)

      implicit none

      integer, intent(in) :: ret

      if(ret /= nf90_noerr) then
         write(6,*) trim(nf90_strerror(ret))
         stop
      endif

   end subroutine sanity

   subroutine build_flx_info

      implicit none

      integer L, j, n

      flxvid(:) = -1
#ifdef RTM
      flxvid_rtm(:) = -1
#endif
#ifdef DGVM
      flxvid_dgvm(:) = -1
#endif

    ! PFT level

      L = 1
      flxname(L) = 'taux'
      flxunit(L) = 'kg/m/s2'
      flxinfo(L) = 'wind stress: E-W' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'tauy'
      flxunit(L) = 'kg/m/s2'
      flxinfo(L) = 'wind stress: N-S' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'fsena'
      flxunit(L) = 'W/m2'
      flxinfo(L) = 'sensible heat from canopy height to atmosphere' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'lfevpa'
      flxunit(L) = 'W/m2'
      flxinfo(L) = 'latent heat flux from canopy height to atmosphere' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'fevpa'
      flxunit(L) = 'mm/s'
      flxinfo(L) = 'evapotranspiration from canopy to atmosphere' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'fsenl'
      flxunit(L) = 'W/m2'
      flxinfo(L) = 'sensible heat from leaves' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'fevpl'
      flxunit(L) = 'mm/s'
      flxinfo(L) = 'evaporation+transpiration from leaves' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'etr'
      flxunit(L) = 'mm/s'
      flxinfo(L) = 'transpiration rate' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'fseng'
      flxunit(L) = 'W/m2'
      flxinfo(L) = 'sensible heat flux from ground' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'fevpg'
      flxunit(L) = 'mm/s'
      flxinfo(L) = 'evaporation heat flux from ground' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'fgrnd'
      flxunit(L) = 'W/m2'
      flxinfo(L) = 'ground heat flux' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'sabvsun'
      flxunit(L) = 'W/m2'
      flxinfo(L) = 'solar absorbed by sunlit canopy' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'sabvsha'
      flxunit(L) = 'W/m2'
      flxinfo(L) = 'solar absorbed by shaded' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'sabg'
      flxunit(L) = 'W/m2'
      flxinfo(L) = 'solar absorbed by ground ' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'olrg'
      flxunit(L) = 'W/m2'
      flxinfo(L) = 'outgoing long-wave radiation from ground+canopy' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'rnet'
      flxunit(L) = 'W/m2'
      flxinfo(L) = 'net radiation' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'zerr'
      flxunit(L) = 'W/m2'
      flxinfo(L) = 'the error of energy balance' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'assim'
      flxunit(L) = 'mol/m2/s'
      flxinfo(L) = 'canopy assimilation rate' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'respc'
      flxunit(L) = 'mol/m2/s'
      flxinfo(L) = 'respiration (plant+soil)' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'fmicr'
      flxunit(L) = 'mol/m2/s'
      flxinfo(L) = 'microbial respiration'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'tlsun'
      flxunit(L) = 'K'
      flxinfo(L) = 'sunlit leaf temperature' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'tlsha'
      flxunit(L) = 'K'
      flxinfo(L) = 'shaded leaf temperature' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'ldew'
      flxunit(L) = 'mm'
      flxinfo(L) = 'depth of water on foliage' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'sigf'
      flxunit(L) = '-'
      flxinfo(L) = 'fraction of veg cover, excluding snow-covered veg' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'green'
      flxunit(L) = '-'
      flxinfo(L) = 'leaf greenness'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'lai'
      flxunit(L) = 'm^2/m^2'
      flxinfo(L) = 'leaf area index'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'sai'
      flxunit(l) = 'm^2/m^2'
      flxinfo(L) = 'stem area index'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'avsdr'
      flxunit(L) = '-'
      flxinfo(L) = 'visible, direct averaged albedo' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'avsdf'
      flxunit(L) = '-'
      flxinfo(L) = 'visible, diffuse averaged albedo' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'anidr'
      flxunit(L) = '-'
      flxinfo(L) = 'near-infrared, direct averaged albedo' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'anidf'
      flxunit(L) = '-'
      flxinfo(L) = 'near-infrared,diffuse averaged albedo' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'emis'
      flxunit(L) = '-'
      flxinfo(L) = 'averaged bulk surface emissivity'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'z0ma'
      flxunit(L) = 'm'
      flxinfo(L) = 'effective roughness' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'trad'
      flxunit(L) = 'K'
      flxinfo(L) = 'radiative temperature of surface' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'ustar'
      flxunit(L) = 'm/s'
      flxinfo(L) = 'u* in similarity theory' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'tstar'
      flxunit(L) = 'kg/kg'
      flxinfo(L) = 't* in similarity theory' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'qstar'
      flxunit(L) = 'kg/kg'
      flxinfo(L) = 'q* in similarity theory' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'zol'
      flxunit(L) = '-'
      flxinfo(L) = 'dimensionless height (z/L) used in Monin-Obukhov theory'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'rib'
      flxunit(L) = '-'
      flxinfo(L) = 'bulk Richardson number in surface layer'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'fm'
      flxunit(L) = '-'
      flxinfo(L) = 'integral of profile function for momentum'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'fh'
      flxunit(L) = '-'
      flxinfo(L) = 'integral of profile function for heat'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'fq'
      flxunit(L) = '-'
      flxinfo(L) = 'integral of profile function for moisture'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'tref'
      flxunit(L) = 'kelvin'
      flxinfo(L) = '2 m height air temperature' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'qref'
      flxunit(L) = 'kg/kg'
      flxinfo(L) = '2 m height air specific humidity' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'u10m'
      flxunit(L) = 'm/s'
      flxinfo(L) = '10m u-velocity' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'v10m'
      flxunit(L) = 'm/s'
      flxinfo(L) = '10m v-velocity' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'f10m'
      flxunit(L) = '-'
      flxinfo(L) = 'integral of profile function for momentum at 10m' 
      flxdims(L) = 2

    ! COLUMN level

      L = L+1
      flxname(L) = 'xerr'
      flxunit(L) = 'mm/s'
      flxinfo(L) = 'the error of water banace' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'rsur'
      flxunit(L) = 'mm/s'
      flxinfo(L) = 'surface runoff' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'rnof'
      flxunit(L) = 'mm/s'
      flxinfo(L) = 'total runoff' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'tss'
      flxunit(L) = 'K'
      flxinfo(L) = 'soil temperature' 
      flxdims(L) = 3

      L = L+1
      flxname(L) = 'wliq'
      flxunit(L) = 'kg/m2'
      flxinfo(L) = 'liquid water in soil layers' 
      flxdims(L) = 3

      L = L+1
      flxname(L) = 'wice'
      flxunit(L) = 'kg/m2'
      flxinfo(L) = 'ice lens in soil layers' 
      flxdims(L) = 3

      L = L+1
      flxname(L) = 'tg'
      flxunit(L) = 'K'
      flxinfo(L) = 'ground surface temperature' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'scv'
      flxunit(L) = 'kg/m2'
      flxinfo(L) = 'snow cover, water equivalent' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'snowdp'
      flxunit(L) = 'meter'
      flxinfo(L) = 'snow depth' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'fsno'
      flxunit(L) = '-'
      flxinfo(L) = 'fraction of snow cover on ground'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'us'
      flxunit(L) = 'm/s'
      flxinfo(L) = 'wind in eastward direction' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'vs'
      flxunit(L) = 'm/s'
      flxinfo(L) = 'wind in northward direction' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'tm'
      flxunit(L) = 'kelvin'
      flxinfo(L) = 'temperature at reference height' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'qm'
      flxunit(L) = 'kg/kg'
      flxinfo(L) = 'specific humidity at reference height' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'prc'
      flxunit(L) = 'mm/s'
      flxinfo(L) = 'convective precipitation' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'prl'
      flxunit(L) = 'mm/s'
      flxinfo(L) = 'large scale precipitation'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'pbot'
      flxunit(L) = 'pa'
      flxinfo(L) = 'atmospheric pressure at the surface'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'frl'
      flxunit(L) = 'W/m2'
      flxinfo(L) = 'atmospheric infrared (longwave) radiation' 
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'solar'
      flxunit(L) = 'W/m2'
      flxinfo(L) = 'downward solar radiation at surface' 
      flxdims(L) = 2

#ifdef CMIP
      L = L+1
      flxname(L) = 'qsulb'
      flxunit(L) = 'kg/m2/s'
      flxinfo(L) = 'sublimation rate from snow pack'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'mrlsl'
      flxunit(L) = 'kg/m2'
      flxinfo(L) = 'mass of water of all phases in each soil layer'
      flxdims(L) = 3

      L = L+1
      flxname(L) = 'mrsos'
      flxunit(L) = 'kg/m2'
      flxinfo(L) = 'mass of water of all phases in the upper 0.1 meters of soil'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'mrso'
      flxunit(L) = 'kg/m2'
      flxinfo(L) = 'mass of water of all phases over all soil layers'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'mrfso'
      flxunit(L) = 'kg/m2'
      flxinfo(L) = 'mass of frozen water over all soil layers'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'lwsnl'
      flxunit(L) = 'kg/m2'
      flxinfo(L) = 'mass of liquid water of snow layers'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'snm'
      flxunit(L) = 'kg/m2/s'
      flxinfo(L) = 'surface snow melt'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'tsn'
      flxunit(L) = 'K'
      flxinfo(L) = 'snow internal temperature'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'nsnow'
      flxunit(L) = '-'
      flxinfo(L) = 'number of snow events'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'treeFrac'
      flxunit(L) = '%'
      flxinfo(L) = 'tree fraction'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'shrubFrac'
      flxunit(L) = '%'
      flxinfo(L) = 'shrub fraction'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'grassFrac'
      flxunit(L) = '%'
      flxinfo(L) = 'grass fraction'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'baresoilFrac'
      flxunit(L) = '%'
      flxinfo(L) = 'bare soil fraction'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'residualFrac'
      flxunit(L) = '%'
      flxinfo(L) = 'residual land fraction'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'soilfrac'
      flxunit(L) = '-'
      flxinfo(L) = 'soil areal fraction on gridlevel'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'urbanfrac'
      flxunit(L) = '-'
      flxinfo(L) = 'urban areal fraction on gridlevel'
      flxdims(L) = 2

      L = L+1
      flxname(L) = 'wetlandfrac'
      flxunit(L) = '-'
      flxinfo(L) = 'wetland areal fraction on gridlevel'
      flxdims(L) = 2
      L = L+1
      flxname(L) = 'icefrac'
      flxunit(L) = '-'
      flxinfo(L) = 'ice areal fraction on gridlevel'
      flxdims(L) = 2
      L = L+1
      flxname(L) = 'lakefrac'
      flxunit(L) = '-'
      flxinfo(L) = 'lake areal fraction on gridlevel'
      flxdims(L) = 2
#endif

      if(L.ne.nflux) stop 'build_flx_info: inconsistent nflux'

      j = 1

      do L = 1, nflux
         flxbeg(L) = j

         if(flxdims(L).eq.2) then
            flxmap(j) = L
            j = j+1
         else if(flxdims(L).eq.3) then
            flxmap(j:j+nl_soil-1) = L
            j = j+nl_soil
         end if
      end do

      if(j.ne.nfldv+1) stop'J.ne.nfldv+1 in build_flx_info'

      do n = 1, nMaxHist
         if(histArray(n)%is_valid) then
            histArray(n)%fldvmask(:) = .false.
   
            j = 1
            do L = 1, nflux
               if(flxdims(L).eq.2) then
                  if(histArray(n)%fluxmask(L)) &
                     histArray(n)%fldvmask(j) = .true.
                  j = j+1
               else if(flxdims(L).eq.3) then
                  if(histArray(n)%fluxmask(L)) &
                     histArray(n)%fldvmask(j:j+nl_soil-1) = .true.
                  j = j+nl_soil
               end if
            end do
         end if
      end do

#ifdef RTM
      L = 1
      flxname_rtm(L) = 'lndrof'
      flxunit_rtm(L) = 'm3/s'
      flxinfo_rtm(L) = 'river flow of land'
      flxdims_rtm(L) = 2

      L = L+1
      flxname_rtm(L) = 'ocnrof'
      flxunit_rtm(L) = 'm3/s'
      flxinfo_rtm(L) = 'river discharge into ocean'
      flxdims_rtm(L) = 2

      do n = 1, nMaxHist
         if(histArray(n)%is_valid) then
            histArray(n)%fldvmask_rtm(:) = .false.

            j = 1
            do L = 1, nflux_rtm
               if(flxdims_rtm(L).eq.2) then
                  if(histArray(n)%fluxmask_rtm(L)) &
                     histArray(n)%fldvmask_rtm(j) = .true.
                  j = j+1
               end if
            end do
         end if
      end do
#endif

#ifdef DGVM

    ! monthly variables

      L = 1
      flxname_dgvm(L) = 'leafc'
      flxunit_dgvm(L) = 'kg/m2'
      flxinfo_dgvm(L) = 'carbon mass in leaves'
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'woodc'
      flxunit_dgvm(L) = 'kg/m2'
      flxinfo_dgvm(L) = 'carbon mass in wood'
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'rootc'
      flxunit_dgvm(L) = 'kg/m2'
      flxinfo_dgvm(L) = 'carbon mass in roots'
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'vegc'
      flxunit_dgvm(L) = 'kg/m2'
      flxinfo_dgvm(L) = 'carbon mass in vegetation'
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'litc_ag'
      flxunit_dgvm(L) = 'kg/m2'
      flxinfo_dgvm(L) = 'carbon mass in above-ground litter'
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'litc_bg'
      flxunit_dgvm(L) = 'kg/m2'
      flxinfo_dgvm(L) = 'carbon mass in below-ground litter'
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'litc'
      flxunit_dgvm(L) = 'kg/m2'
      flxinfo_dgvm(L) = 'carbon mass in litter pool'
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'soic_fast'
      flxunit_dgvm(L) = 'kg/m2'
      flxinfo_dgvm(L) = 'carbon mass in fast soil pool'
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'soic_slow'
      flxunit_dgvm(L) = 'kg/m2'
      flxinfo_dgvm(L) = 'carbon mass in slow soil pool'
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'soic'
      flxunit_dgvm(L) = 'kg/m2'
      flxinfo_dgvm(L) = 'carbon mass in soil pool'
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'fveg2litter'
      flxunit_dgvm(L) = 'kg/m2/s'
      flxinfo_dgvm(L) = 'total carbon mass flux from vegetation to litter'
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'flitter2soil'
      flxunit_dgvm(L) = 'kg/m2/s'
      flxinfo_dgvm(L) = 'total carbon mass flux from litter to soil'
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'flitter2atmos'
      flxunit_dgvm(L) = 'kg/m2/s'
      flxinfo_dgvm(L) = 'total carbon mass flux from litter to atmosphere'
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'gpp'
      flxunit_dgvm(L) = 'kg/m2/s'
      flxinfo_dgvm(L) = 'carbon mass flux out of atmosphere due to GPP on land'
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'npp'
      flxunit_dgvm(L) = 'kg/m2/s'
      flxinfo_dgvm(L) = 'carbon mass flux out of atmosphere due to NPP on land'
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'nep'
      flxunit_dgvm(L) = 'kg/m2/s'
      flxinfo_dgvm(L) = 'net carbon mass flux out of atmosphere due to NEP on land'
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'nbp'
      flxunit_dgvm(L) = 'kg/m2/s'
      flxinfo_dgvm(L) = 'carbon mass flux out of atmosphere due to NBP on land'
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'ra'
      flxunit_dgvm(L) = 'kg/m2/s'
      flxinfo_dgvm(L) = 'carbon mass flux into atmosphere due to autotrophic(plant) respiration on land'
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'rh'
      flxunit_dgvm(L) = 'kg/m2/s'
      flxinfo_dgvm(L) = 'carbon mass flux into atmosphere due to heterotrophic respiration on land'
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'ffirec'
      flxunit_dgvm(L) = 'kg/m2/s'
      flxinfo_dgvm(L) = 'carbon mass flux into atmosphere due to CO2 emission from fire'
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'wt_pft'
      flxunit_dgvm(L) = '%'
      flxinfo_dgvm(L) = 'faction of each PFT on gridlevel'
      flxdims_dgvm(L) = 3

    ! yearly variables

      L = L+1
      flxname_dgvm(L) = 'fpc'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 3

      L = L+1
      flxname_dgvm(L) = 'bare'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'afirec'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'afiref'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'avegc'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'aestabc'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'anpp'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'amrh'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'alitcag'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'alitcbg'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'asoicf'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'asoics'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'npp_ind'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 3

      L = L+1
      flxname_dgvm(L) = 'lm'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 3

      L = L+1
      flxname_dgvm(L) = 'sm'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 3

      L = L+1
      flxname_dgvm(L) = 'hm'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 3

      L = L+1
      flxname_dgvm(L) = 'rm'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 3

      L = L+1
      flxname_dgvm(L) = 'ca'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 3

      L = L+1
      flxname_dgvm(L) = 'htop'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 3

      L = L+1
      flxname_dgvm(L) = 'nind'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 3

      L = L+1
      flxname_dgvm(L) = 'laimx'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 3

      L = L+1
      flxname_dgvm(L) = 'anngpp'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 3

      L = L+1
      flxname_dgvm(L) = 'annfrmf'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 3

      L = L+1
      flxname_dgvm(L) = 'annfrms'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 3

      L = L+1
      flxname_dgvm(L) = 'annfrmr'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 3

      L = L+1
      flxname_dgvm(L) = 'annfrg'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 3

#ifdef DyN
      L = L+1
      flxname_dgvm(L) = 'cnleaf'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 3

      L = L+1
      flxname_dgvm(L) = 'cnsap'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 3

      L = L+1
      flxname_dgvm(L) = 'cnroot'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 3

      L = L+1
      flxname_dgvm(L) = 'an_up'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'stress'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'avegn'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'alitnag'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'alitnbg'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'asoin'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'no3'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 2

      L = L+1
      flxname_dgvm(L) = 'nh4'
      flxunit_dgvm(L) = '-'
      flxinfo_dgvm(L) = ''
      flxdims_dgvm(L) = 2
#endif

      if(L.ne.nflux_dgvm) stop 'build_flx_info: inconsistent nflux_dgvm'

      j = 1

      do L = 1, nflux_dgvm
         flxbeg_dgvm(L) = j

         if(flxdims_dgvm(L).eq.2) then
            flxmap_dgvm(j) = L
            j = j+1
         else if(flxdims_dgvm(L).eq.3) then
            flxmap_dgvm(j:j+numpft_nat-1) = L
            j = j+numpft_nat
         end if
      end do

      if(j.ne.nfldv_dgvm+1) stop 'J.ne.nfldv_dgvm+1 in build_flx_info'

      do n = 1, nMaxHist
         if(histArray(n)%is_valid) then
            histArray(n)%fldvmask_dgvm(:) = .false.
   
            j = 1
            do L = 1, nflux_dgvm
               if(flxdims_dgvm(L).eq.2) then
                  if(histArray(n)%fluxmask_dgvm(L)) &
                     histArray(n)%fldvmask_dgvm(j) = .true.
                  j = j+1
               else if(flxdims_dgvm(L).eq.3) then
                  if(histArray(n)%fluxmask_dgvm(L)) &
                     histArray(n)%fldvmask_dgvm(j:j+numpft_nat-1) = .true.
                  j = j+numpft_nat
               end if
            end do
         end if
      end do
#endif

   end subroutine build_flx_info

end module nchistMod

