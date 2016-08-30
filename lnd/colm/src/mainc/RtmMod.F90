#include <define.h>

module RtmMod

#if (defined RTM)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: RtmMod
!
! !DESCRIPTION:
! River Routing Model (U. of Texas River Transport
! Model)~\cite{Branstetter:2001}
!
! !USES:
  use precision   , only : r8
  use colm_varMod , only : lsmlon => lon_points, lsmlat => lat_points
  use colm_rtmVar
  use shr_sys_mod , only : shr_sys_flush
  use abortutils  , only : endrun
  use spmd, only : p_iam
!
! !PUBLIC TYPES:
  implicit none
  save
  integer , parameter, public :: rtmloni = 1        ! RTM grid - per-proc beginning lon index
  integer , parameter, public :: rtmlonf = rtmlon   ! RTM grid - per-proc ending lon index
  integer , parameter, public :: rtmlati = 1        ! RTM grid - per-proc beginning lat index
  integer , parameter, public :: rtmlatf = rtmlat   ! RTM grid - per-proc ending lat index
  real(r8), public, pointer   :: latixy_r(:,:)      ! RTM grid - latitudes  of grid cells (degrees)
  real(r8), public, pointer   :: longxy_r(:,:)      ! RTM grid - longitudes of grid cells (degrees)
  real(r8), public, pointer   :: area_r(:,:)        ! RTM grid - gridcell area (km^2)
  integer , public, pointer   :: mask_r(:,:)        ! RTM grid - landmask (land=1,ocean=0)
  real(r8), allocatable, public :: volr(:,:)        ! RTM fluxes - water volume in cell (m^3)
!
! !PUBLIC MEMBER FUNCTIONS:
  public Rtmgridini    ! Initialize RTM grid and land mask
  public Rtmlandini    ! Initialize RTM-land interpolation weights
  public Rtmfluxini    ! Initialize RTM fluxout
  public Rtmriverflux  ! Interface with RTM river routing model
  public restart_rtm   ! Read/write RTM restart data
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!
! !PRIVATE MEMBER FUNCTIONS:
!
  private UpdateInput   ! Update rtm inputs
  private Rtm           ! River routing model (based on U. Texas code)
  private UpdateGlobal  ! Update global quantities
!
! !PRIVATE TYPES:
!
! RTM inputs at land model grid resolution
! rtmin_ave(1,:), rtmin_glob(1,:) and rtmin_loc(1,:) correspond to totrunin
! rtmin_ave(2,:), rtmin_glob(2,:) and rtmin_loc(2,:) correspond to prec
! rtmin_ave(3,:), rtmin_glob(3,:) and rtmin_loc(3,:) correspond to evap
!
  private
  real(r8), pointer :: rtmin_ave(:,:)        ! RTM local averaging buffer for runoff, prec, evap
  real(r8), pointer :: rtmin_glob(:,:)       ! RTM global input
  real(r8), pointer :: rtmin(:,:)            ! RTM local input
  integer  :: ncount_rtm                     ! RTM time averaging = number of time samples to average over
  real(r8) :: delt_rtm                       ! RTM time step
!
! RTM 1/2 degree resolution variables
!
  real(r8), parameter :: effvel = 0.35       ! RTM effective velocity (m/s)
  integer :: numlon_r(rtmlat)                ! RTM grid number of lon points at each lat
  integer :: mxovr_s2r                       ! RTM mapping - max number of overlapping cells
  integer :: novr_s2r(rtmlon,rtmlat)         ! RTM mapping - number    of overlapping lsm cells
  integer , allocatable :: iovr_s2r(:,:,:)   ! RTM mapping - lon index of overlapping land model cells
  integer , allocatable :: jovr_s2r(:,:,:)   ! RTM mapping - lat index of overlapping land model cells
  real(r8), allocatable :: wovr_s2r(:,:,:)   ! RTM mapping - weight    of overlapping land model cells
  real(r8), allocatable :: latsh(:)          ! RTM grid - southern edge of cells at rtm grid
  real(r8), allocatable :: lonwh(:,:)        ! RTM grid - western  edge of cells at rtm grid
  integer , allocatable :: rdirc(:,:)        ! RTM input - rtm river flow direction (0-8)
  real(r8), allocatable :: ddist(:,:)        ! RTM input - downstream distance (m)
  real(r8), allocatable :: rivarea(:,:)      ! RTM input - cell area (m^2)
  real(r8), allocatable :: totrunin_r(:,:)   ! RTM input - surface runoff (mm/s)
  real(r8), allocatable :: sfluxin(:,:)      ! RTM input - water flux into cell (m3/s)
  real(r8), allocatable :: fluxout(:,:)      ! RTM input/output - water flux out of cell (m^3/s)
  real(r8), allocatable :: runrtm(:,:)       ! RTM input  - input runoff on rtm grid (m**3/s)
  real(r8), allocatable :: flxlnd_r(:,:)     ! RTM output - river flux (m**3/s)
  real(r8), allocatable :: flxocn_r(:,:)     ! RTM output - river flux to the ocean (m**3/s)
  real(r8), allocatable :: dvolrdt_r(:,:)    ! RTM output - change in storage (mm/s)
  real(r8), allocatable :: volrtm(:,:)       ! RTM output - change in storage (m**3/s)
  real(r8) :: prec_global                    ! RTM global averaging - total precipitation (m^3/sec)
  real(r8) :: evap_global                    ! RTM global averaging - total evaporation (m^3/sec)
  real(r8) :: runlnd_global                  ! RTM global averaging - total input runoff on land grid (m^3/sec)
  real(r8) :: runrtm_global                  ! RTM global averaging - total input runoff on rtm grid (m^3/sec)
  real(r8) :: ocnrtm_global                  ! RTM global averaging - total ocean runoff on rtm grid (m^3/sec)
  real(r8) :: volrtm_global                  ! RTM global averaging - total change in storage on rtm (m^3/sec)
  integer  :: ncount_global                  ! RTM global averaging - global counter
  integer  :: yrold                          ! RTM global averaging - old year
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtmgridini
!
! !INTERFACE:
  subroutine Rtmgridini
!
! !DESCRIPTION:
! Initialize RTM grid and land mask.
!
! !USES:
    use precision    , only : r8
#if (defined SPMD)
    use spmd         , only : mpicom => p_comm, masterproc => p_master, MPI_REAL8, MPI_INTEGER
#else
    use spmd         , only : masterproc => p_master
#endif
    use areaMod      , only : celledge, cellarea
    use colm_rtmVar  , only : frivinp_rtm
    use phycon_module, only : re
    use shr_const_mod, only : SHR_CONST_PI
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!
! !LOCAL VARIABLES:
    real(r8), dimension(4) :: rtmedge = (/ 90., 180., -90., -180. /)  !N,E,S,W edges of rtm grid
    integer  :: ioff(0:8) = (/0,0,1,1,1,0,-1,-1,-1/) !calc dist as in hydra
    integer  :: joff(0:8) = (/0,1,1,0,-1,-1,-1,0,1/) !of grid cell down stream
    integer  :: i,j,k,n                       !loop indices
    integer  :: i2,j2                         !downstream i and j
    real(r8) :: deg2rad                       !pi/180
    real(r8) :: dx                            !lon dist. between grid cells (m)
    real(r8) :: dy                            !lat dist. between grid cells (m)
    real(r8) :: tempg(rtmlon,rtmlat)          !temporary buffer
    integer  :: tempgp(0:rtmlon+1,0:rtmlat+1) !temporary buffer
    integer  :: ier                           !error code
    character(len=16), dimension(50) :: river_name
    character(len=30), dimension(50) :: rivstat_name
    real(r8)         , dimension(50) :: rivstat_lon
    real(r8)         , dimension(50) :: rivstat_lat
!-----------------------------------------------------------------------

    ! Allocate rtm grid variables

    allocate (latixy_r(rtmlon,rtmlat), longxy_r(rtmlon,rtmlat), &
              latsh(rtmlat+1), lonwh(rtmlon+1,rtmlat), &
              area_r(rtmlon,rtmlat), mask_r(rtmlon,rtmlat), stat=ier)
    if (ier /= 0) then
       write(6,*)'Rtmgridini: Allocation error for ',&
            'latixy_r, longxy_r, latsy, lonwh, area_r, mask_r'
       call endrun
    end if

    ! Allocate rtm flux variables

    allocate (volr(rtmloni:rtmlonf,rtmlati:rtmlatf),&
              rdirc(rtmloni-1:rtmlonf+1,rtmlati-1:rtmlatf+1), &
              fluxout(rtmloni-1:rtmlonf+1,rtmlati-1:rtmlatf+1), &
              ddist(rtmloni:rtmlonf,rtmlati:rtmlatf), &
              rivarea(rtmloni:rtmlonf,rtmlati:rtmlatf), stat=ier)
    if (ier /= 0) then
       write(6,*)'Rtmgridini: Allocation error for ',&
            'rdirec, fluxout, ddist, rivarea'
       call endrun
    end if

    ! Allocate inputs and outputs to rtm at 1/2 degree resolution

    allocate (totrunin_r(rtmloni:rtmlonf,rtmlati:rtmlatf), &
              flxlnd_r(rtmloni:rtmlonf,rtmlati:rtmlatf), &
              flxocn_r(rtmloni:rtmlonf,rtmlati:rtmlatf), &
              dvolrdt_r(rtmloni:rtmlonf,rtmlati:rtmlatf), &
              sfluxin(rtmloni:rtmlonf,rtmlati:rtmlatf),  stat=ier)
    if (ier /= 0) then
       write(6,*)'Rtmgridini: Allocation error for ',&
            'totrunin_r, flxlnd_r, flxocn_r, dvolrdt_r, sfluxin'
       call endrun
    end if

    ! Useful constants and initial values

    deg2rad = SHR_CONST_PI / 180.
    volr(:,:) = 0.
    fluxout(:,:) = 0.
    flxocn_r(:,:) = 0.
    flxlnd_r(:,:) = 0.

    ! Open and read input data (river direction file)
    ! rtm operates from south to north and from the dateline
    ! River station data is currently not used by the model -
    ! it is only used by the diagnostic package
    ! If the river direction file is modified - the river station
    ! part must also be modified

    if (masterproc) then
       write(6,*)'Columns in RTM = ',rtmlon
       write(6,*)'Rows in RTM    = ',rtmlat

       open (10,file=frivinp_rtm)
       write(6,*)'opened river direction data'
       do j = 1,rtmlat
          numlon_r(j) = 0
          do i = 1,rtmlon
             read(10,*) latixy_r(i,j),longxy_r(i,j),tempg(i,j)
             if (longxy_r(i,j) /= 1.e36) numlon_r(j) = numlon_r(j) + 1
             tempgp(i,j) = nint(tempg(i,j))
          enddo
       enddo
       do n = 1,50
          read(10,10,iostat=ier) river_name(n), rivstat_lon(n), rivstat_lat(n), rivstat_name(n)
          if (ier /= 0) exit
10        format(1x,a16,f7.2,1x,f7.2,a30)
       end do
       close(10)
       write(6,*)'closed river direction data'
       write(6,*)
    endif

#if (defined SPMD)
    call mpi_bcast(numlon_r, size(numlon_r), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast(latixy_r, size(latixy_r), MPI_REAL8,   0, mpicom, ier)
    call mpi_bcast(longxy_r, size(longxy_r), MPI_REAL8,   0, mpicom, ier)
    call mpi_bcast(tempgp  , size(tempgp)  , MPI_INTEGER, 0, mpicom, ier)
#endif

    ! Determine RTM celledges, areas and interpolation masks

    call celledge (rtmlat    , rtmlon    , numlon_r  , longxy_r  , &
                   latixy_r  , rtmedge(1), rtmedge(2), rtmedge(3), &
                   rtmedge(4), latsh     , lonwh     )

    call cellarea (rtmlat    , rtmlon    , numlon_r  , latsh     , lonwh , &
                   rtmedge(1), rtmedge(2), rtmedge(3), rtmedge(4), area_r)

    ! Determine rtm mask, downstream distance and area

    do i=1,rtmlon
       tempgp(i,0)        = tempgp(mod(i+rtmlon/2-1,rtmlon)+1,1)
       tempgp(i,rtmlat+1) = tempgp(mod(i+rtmlon/2-1,rtmlon)+1,rtmlat)
       if (tempgp(i,0)        /= 0) tempgp(i,0)        = mod(tempgp(i,0)       +4-1,8)+1
       if (tempgp(i,rtmlat+1) /= 0) tempgp(i,rtmlat+1) = mod(tempgp(i,rtmlat+1)+4-1,8)+1
    enddo
    do j=0,rtmlat+1
       tempgp(0,j) =tempgp(rtmlon,j)
       tempgp(rtmlon+1,j)=tempgp(1,j)
    enddo

    ! Determine rtm river flow direction (0-8)

    do j=rtmlati-1,rtmlatf+1
       do i=rtmloni-1,rtmlonf+1
          rdirc(i,j)=tempgp(i,j)
       enddo
    enddo

    ! Determine rtm ocn/land mask

    do j=rtmlati,rtmlatf
       do i=rtmloni,rtmlonf
          if (rdirc(i,j) == 0) then
             mask_r(i,j) = 0
          else
             mask_r(i,j) = 1
          end if
       enddo
    enddo

    ! Determine downstream distance - instead of reading a distance file
    ! calculate the downstream distance

    do j=rtmlati,rtmlatf
       do i=rtmloni,rtmlonf
          i2 = i + ioff(tempgp(i,j))
          j2 = j + joff(tempgp(i,j))
          if (i2 == 0) i2 = 2                 !avoids i2 out of bounds in the following
          if (i2 == rtmlon+1) i2 = rtmlon-1   !avoids i2 out of bounds in the following
          dy = deg2rad * abs(latixy_r(i,j)-latixy_r(i2,j2)) * re*1000.
          dx = deg2rad * abs(longxy_r(i,j)-longxy_r(i2,j2)) * re*1000. &
               *0.5*(cos(latixy_r(i,j)*deg2rad)+cos(latixy_r(i2,j2)*deg2rad))
          ddist(i,j) = sqrt(dx*dx + dy*dy)
          rivarea(i,j)=1.e6 * area_r(i,j)     !convert into m**2
       enddo
    enddo

  end subroutine Rtmgridini

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtmlandini
!
! !INTERFACE:
  subroutine Rtmlandini
!
! !DESCRIPTION:
! Initialize RTM-land interpolation weights
! and variables related to runoff time averaging.
!
! !USES:
    use precision   , only : r8
    use colm_varMod , only : numlon, area, lats, lonw, landmask, numgrid, numgrid_glob
    use spmd        , only : masterproc => p_master
    use areaMod     , only : areaini_point, mkmxovr
    use timemgr     , only : get_curr_date
    use RunoffMod   , only : set_proc_rof_bounds, set_roflnd, set_rofocn
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
!@jidy
!   integer , pointer :: ixy(:)              ! gricell xy lon index
!   integer , pointer :: jxy(:)              ! gricell xy lat index
!@jidy
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: g,i,j,k,n                    ! indices
    integer  :: is,js                        ! land model grid indices
    integer  :: ir,jr                        ! rtm grid indices
    integer  :: pid                          ! processor id
    integer  :: ier                          ! error code
    real(r8) :: maskone_s(lsmlon,lsmlat)     ! dummy field: see below
    real(r8) :: maskone_r(rtmlon,rtmlat)     ! dummy field: see below
    integer  :: masktmp_r(rtmlon,rtmlat)     ! dummy mask
    integer  :: cplrof_mask(rtmlon,rtmlat)   ! rtm mask for ocean points with possible nonzero runoff
    integer  :: mon                          ! month (1, ..., 12)
    integer  :: day                          ! day of month (1, ..., 31)
    integer  :: ncsec                        ! seconds of current date
    real(r8) :: offset                       ! offset for interpolation from model->rtm grid
    real(r8) :: lonw_offset(lsmlon+1,lsmlat) ! longitudinal offset for interpolation from model->rtm grid
    integer  :: novr_i2o                     ! number of overlapping land model cells in given rtm cell
  !@jidy
  ! integer  :: begp, endp                   ! per-proc beginning and ending pft indices
  ! integer  :: begc, endc                   ! per-proc beginning and ending column indices
  ! integer  :: begl, endl                   ! per-proc beginning and ending landunit indices
  ! integer  :: begg, endg                   ! per-proc gridcell ending gridcell indices
  ! integer  :: numg                         ! total number of gridcells across all processors
  ! integer  :: numl                         ! total number of landunits across all processors
  ! integer  :: numc                         ! total number of columns across all processors
  ! integer  :: nump                         ! total number of pfts across all processors
  !@jidy
    integer  :: nroflnd
    integer  :: nrofocn
    integer , pointer :: iovr_i2o(:)         ! lon index of overlap input cell
    integer , pointer :: jovr_i2o(:)         ! lat index of overlap input cell
    real(r8), pointer :: wovr_i2o(:)         ! weight    of overlap input cell
!-----------------------------------------------------------------------

   ! Assign local pointers to derived type members (gridcell-level)

!@jidy
!   ixy => clm3%g%ixy
!   jxy => clm3%g%jxy
!@jidy

    ! --------------------------------------------------------------------
    ! The following section allows RTM and land model to coexist at different
    ! horizontal resolutions
    ! --------------------------------------------------------------------

    if (masterproc) then
       write(6,*)
       write(6,*) 'Initializing area-averaging interpolation for RTM.....'
    endif

    ! To find fraction of each land model grid cell that is land based on rtm grid.
    ! For this purpose, want all rtm grid cells to contribute to grid cell
    ! average on land model grid, i.e., all cells used regardless of whether land
    ! or ocean. Do this by setting [maskone_s] = 1

    ! [maskone_s]=1 means all grid cells on land model grid, regardless of whether
    ! land or ocean, will contribute to rtm grid.

    do j = 1,lsmlat
       do i = 1,numlon(j)
          maskone_s(i,j) = 1.
       end do
    end do

    ! [maskone_r] = 1 means all the rtm grid is land. Used as dummy
    ! variable so code will not abort with false, non-valid error check

    do j = rtmlati,rtmlatf
       do i = rtmloni,rtmlonf
          maskone_r(i,j) = 1.
       end do
    end do

    ! --------------------------------------------------------------------
    ! Map weights from land model grid to rtm grid
    ! --------------------------------------------------------------------

    if (masterproc) then
       write(6,*) 'Initializing land model -> rtm interpolation .....'
    endif

    ! For each rtm grid cell: get lat [jovr_s2r] and lon [iovr_s2r] indices
    ! and weights [wovr_s2r] of overlapping atm grid cells

    call mkmxovr (lsmlon, lsmlat, numlon  , lonw , lats , &
                  rtmlon, rtmlat, numlon_r, lonwh, latsh, &
                  mxovr_s2r     , novr_s2r)

    allocate(iovr_s2r(rtmloni:rtmlonf,rtmlati:rtmlatf,mxovr_s2r), &
             jovr_s2r(rtmloni:rtmlonf,rtmlati:rtmlatf,mxovr_s2r), &
             wovr_s2r(rtmloni:rtmlonf,rtmlati:rtmlatf,mxovr_s2r), &
             iovr_i2o(mxovr_s2r), jovr_i2o(mxovr_s2r), wovr_i2o(mxovr_s2r), &
             stat=ier)
    if (ier /= 0) then
       write(6,*)'Rtmlndini: Allocation error for ',&
            'iovr_s2r, jovr_s2r, wovr_s2r, iovr_i2o, jovr_i2o, wovr_i2o'
       call endrun
    end if

    ! Shift x-grid to locate periodic grid intersections. This
    ! assumes that all lonw(1,j) have the same value for all
    ! latitudes j and that the same holds for lonwh(1,j)

    if (lonw(1,1) < lonwh(1,1)) then
       offset = 360.0
    else
       offset = -360.0
    end if
    do js = 1, lsmlat
       do is = 1, numlon(js) + 1
          lonw_offset(is,js) = lonw(is,js) + offset
       end do
    end do

    ! Determine overlap indices and weights

    do jr = rtmlati,rtmlatf
       do ir = rtmloni,rtmlonf

          call areaini_point (ir           , jr         , lsmlon  , lsmlat  , numlon   , &
                              lonw         , lonw_offset, lats    , area    , maskone_s, &
                              rtmlon       , rtmlat     , numlon_r, lonwh   , latsh    , &
                              area_r(ir,jr), maskone_r(ir,jr), novr_i2o, iovr_i2o, jovr_i2o , &
                              wovr_i2o     , mxovr_s2r)

          if (novr_i2o /= novr_s2r(ir,jr)) then
             write(6,*)'Rtmlandini error: novr_i2o= ',novr_i2o,&
                  ' not equal to  novr_s2r ',novr_s2r(ir,jr),&
                  ' at ir,jr=',ir,jr
             call endrun
          else if (novr_i2o > mxovr_s2r) then
             write(6,*)'Rtmlandini error: novr_s2r= ',novr_s2r,&
                  ' greater than mxovr_s2r= ',mxovr_s2r
             call endrun
          endif

          do n = 1,novr_i2o
             iovr_s2r(ir,jr,n) = iovr_i2o(n)
             jovr_s2r(ir,jr,n) = jovr_i2o(n)
             wovr_s2r(ir,jr,n) = wovr_i2o(n)
          end do

       end do
    end do

    if (masterproc) then
       write(6,*) 'Successfully made land model -> rtm interpolation'
       write(6,*)
    endif

    ! Determine which ocean cells have runoff values
    ! First loop over all ocean points and determine which are at the
    ! end of rivers by examining if any neighboring points are land and
    ! if that land neighbor points into this ocean point. Next loop over all
    ! ocean points and determine which overlap with at least one land cell.
    ! Allocate ocean runoff vector and indices and determine indices
    ! need to reset cpl runoff size to 0 and do the counting again because need to first
    ! first count to allocate vector and must now count to actually determine indices

    nrofocn = 0
    masktmp_r(:,:) = 0
    do j=rtmlati,rtmlatf
       do i=rtmloni,rtmlonf
          if (mask_r(i,j) == 0) then
             if (rdirc(i  ,j-1)==1) masktmp_r(i,j) = 1
             if (rdirc(i-1,j-1)==2) masktmp_r(i,j) = 1
             if (rdirc(i-1,j  )==3) masktmp_r(i,j) = 1
             if (rdirc(i-1,j+1)==4) masktmp_r(i,j) = 1
             if (rdirc(i  ,j+1)==5) masktmp_r(i,j) = 1
             if (rdirc(i+1,j+1)==6) masktmp_r(i,j) = 1
             if (rdirc(i+1,j  )==7) masktmp_r(i,j) = 1
             if (rdirc(i+1,j-1)==8) masktmp_r(i,j) = 1
             if (masktmp_r(i,j) == 0) then
                do n=1,novr_s2r(i,j)
                   is = iovr_s2r(i,j,n)
                   js = jovr_s2r(i,j,n)
                   if (landmask(is,js)==1 .and. wovr_s2r(i,j,n)>0.) then
                      masktmp_r(i,j) = 1
                   end if
                end do
             endif
          endif
          if (masktmp_r(i,j) == 1) nrofocn = nrofocn +1
       enddo
    enddo

    call set_rofocn(rtmlon, rtmlat, masktmp_r, area_r, nrofocn)

    ! Determine which land cells have runoff values

    nroflnd = 0
    masktmp_r(:,:) = 0
    do j=rtmlati,rtmlatf
       do i=rtmloni,rtmlonf
          if (mask_r(i,j) == 1) then
             masktmp_r(i,j) = 1
             nroflnd = nroflnd +1
          end if
       end do
    end do

    call set_roflnd(rtmlon, rtmlat, masktmp_r, area_r, nroflnd)

    ! Deallocate memory for rtm grid  - needed to be done here because
    ! rtm grid information had to be sent to coupler between calls to
    ! Rtmgridini and Rtmlandini

    deallocate(latsh, lonwh, iovr_i2o, jovr_i2o, wovr_i2o)

    ! Determine per-processor runoff bounds

    call set_proc_rof_bounds()

    ! Allocate and initialize dynamic memory for rtm inputs

!@jidy
!   call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
!   call get_proc_global(numg, numl, numc, nump)
!@jidy

!@jidy - for consistance with CoLM's structure.
!   allocate (rtmin(3,numg), rtmin_glob(3,numg), rtmin_ave(3,numg), stat=ier)
    allocate (rtmin(3,numgrid), rtmin_glob(3,numgrid_glob), rtmin_ave(3,numgrid), stat=ier)
!@jidy
    if (ier /= 0) then
       write(6,*)'Rtmlandini: Allocation error for rtmin, rtmin_glob, rtmin_ave'
       call endrun
    end if
    rtmin(:,:)      = 0.
    rtmin_glob(:,:) = 0.
    rtmin_ave(:,:)  = 0.

    ! Determine current date info

    call get_curr_date(yrold, mon, day, ncsec)

    ! Initialize rtm time averaging variables.  Upon restart,
    ! the following variables will get new values  from the restart file.

    ncount_rtm    = 0
    ncount_global = 0
    prec_global   = 0.
    evap_global   = 0.
    runlnd_global = 0.
    runrtm_global = 0.
    volrtm_global = 0.
    ocnrtm_global = 0.

  end subroutine Rtmlandini

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtmfluxini()
!
! !INTERFACE:
  subroutine Rtmfluxini()
!
! !DESCRIPTION:
! Initialize RTM fluxout for case of initial run when initial data is
! read in. For restart run, RTM fluxout is read from restart dataset.
!
! !USES:
    use precision   , only : r8
!@jidy
!   use time_manager, only : get_step_size
!   use clm_varctl  , only : rtm_nsteps
    use timemgr     , only : get_step_size
    use colm_rtmVar , only : rtm_nsteps
!@jidy
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j         !indices
    integer :: delt        !delt for rtm
!-----------------------------------------------------------------------

    delt = rtm_nsteps*get_step_size()
    do j = rtmlati,rtmlatf
       do i = rtmloni,rtmlonf
          if (mask_r(i,j)==1) then
             fluxout(i,j) = volr(i,j) * effvel/ddist(i,j)
             fluxout(i,j) = min(fluxout(i,j), volr(i,j) / delt)
          else
             fluxout(i,j) = 0.
          endif
       enddo
    enddo

  end subroutine Rtmfluxini

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtmriverflux
!
! !INTERFACE:
  subroutine Rtmriverflux()
!
! !DESCRIPTION:
! Interface with RTM river routing model.
!
! !USES:
    use precision   , only : r8
    use spmd        , only : masterproc => p_master
    use colm_varMod , only : landfrac
    use RunoffMod   , only : UpdateRunoff
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!
    logical  :: do_rtm                       ! true => perform rtm calculation
    integer  :: i,j,k,n,g,l,c,p              ! indices
    integer  :: io,jo,ir,jr,is,js            ! mapping indices
    real(r8) :: wt                           ! weight
    real(r8) :: precxy(lsmlon,lsmlat)        ! precipitation (mm H2O /s)
    real(r8) :: evapxy(lsmlon,lsmlat)        ! evaporation (mm H2O /s)
    real(r8) :: totruninxy(lsmlon,lsmlat)    ! surface runoff (mm H2O /s)
!-----------------------------------------------------------------------

    ! Determine RTM inputs on land model grid

    call UpdateInput(do_rtm, precxy, evapxy, totruninxy)

    if (do_rtm) then

       ! Map RTM inputs from land model grid to RTM grid (1/2 degree resolution)

!$OMP PARALLEL DO PRIVATE (jr,ir,n,is,js,wt)
!CSD$ PARALLEL DO PRIVATE (jr,ir,n,is,js,wt)
       do jr = rtmlati,rtmlatf
          do ir = rtmloni,rtmlonf
             totrunin_r(ir,jr) = 0.
          end do
          do n = 1, mxovr_s2r
             do ir = rtmloni,rtmlonf
                if (n > novr_s2r(ir,jr)) cycle
                if (wovr_s2r(ir,jr,n) > 0.) then
                   is = iovr_s2r(ir,jr,n)
                   js = jovr_s2r(ir,jr,n)
                   wt = wovr_s2r(ir,jr,n)
                   totrunin_r(ir,jr) = totrunin_r(ir,jr) + wt*totruninxy(is,js)*landfrac(is,js)
                end if
             end do
          end do
       end do
!CSD$ END PARALLEL DO
!$OMP END PARALLEL DO

       ! Determine RTM runoff fluxes

     ! call t_startf('rtm_calc')
       call Rtm()
     ! call t_stopf('rtm_calc')

       ! Update coupler and history file info
     ! call t_startf('rtm_update')
       call UpdateRunoff(rtmloni, rtmlonf, rtmlati, rtmlatf, flxocn_r, flxlnd_r)
     ! call t_stopf('rtm_update')

       ! Determine global quantities and increment global counter

       if (masterproc) then
        ! call t_startf('rtm_global')
          call UpdateGlobal(totruninxy, precxy, evapxy)
        ! call t_stopf('rtm_global')
       end if

    end if

  end subroutine Rtmriverflux

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: UpdateInput
!
! !INTERFACE:
  subroutine UpdateInput(do_rtm, precxy, evapxy, totruninxy)
!
! !DESCRIPTION:
! Update RTM inputs.
!
! !USES:
    use precision      , only : r8
  ! use clmtype
  ! use decompMod      , only : get_proc_bounds, get_proc_global
  ! use clm_varpar     , only : lsmlon, lsmlat
  ! use clm_varsur     , only : area
  ! use clm_varctl     , only : rtm_nsteps
  ! use time_manager   , only : get_step_size, get_nstep
    use spmd_decomp    , only : gxmap_glob, gymap_glob
    use colm_varMod    , only : area, numgrid, numgrid_glob
    use colm_rtmVar    , only : rtm_nsteps, rnof, prc, prl, fevpa
    use timemgr        , only : get_step_size, get_nstep
#if (defined SPMD)
    use spmd           , only : masterproc => p_master, mpicom => p_comm
    use spmdGathScatMod, only : scatter_data_from_master, gather_data_to_master, allgather_data
#else
    use spmd           , only : masterproc => p_master
#endif
!
! !ARGUMENTS:
    implicit none
    logical , intent(out) :: do_rtm
    real(r8), intent(out) :: precxy(lsmlon,lsmlat)     ! precipitation (mm H2O /s)
    real(r8), intent(out) :: evapxy(lsmlon,lsmlat)     ! evaporation (mm H2O /s)
    real(r8), intent(out) :: totruninxy(lsmlon,lsmlat) ! surface runoff (mm H2O /s)
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Author: Sam Levis
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
!   integer , pointer :: ixy(:)               ! gricell xy lon index
!   integer , pointer :: jxy(:)               ! gricell xy lat index
!   integer , pointer :: cgridcell(:)         ! corresponding gridcell index for each column
!   real(r8), pointer :: wtgcell(:)           ! weight (relative to gridcell) for each column (0-1)
!   real(r8), pointer :: forc_rain(:)         ! rain rate [mm/s]
!   real(r8), pointer :: forc_snow(:)         ! snow rate [mm/s]
!   real(r8), pointer :: qflx_qrgwl(:)        ! qflx_surf at glaciers, wetlands, lakes
!   real(r8), pointer :: qflx_drain(:)        ! sub-surface runoff (mm H2O /s)
!   real(r8), pointer :: qflx_evap_tot(:)     ! qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
!   real(r8), pointer :: qflx_surf(:)         ! surface runoff (mm H2O /s)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer :: i,j,k,n,g,l,c,p                ! indices
    integer :: io,jo,ir,jr,is,js              ! mapping indices
!   integer :: begp, endp                     ! per-proc beginning and ending pft indices
!   integer :: begc, endc                     ! per-proc beginning and ending column indices
!   integer :: begl, endl                     ! per-proc beginning and ending landunit indices
!   integer :: begg, endg                     ! per-proc gridcell ending gridcell indices
!   integer :: numg                           ! total number of gridcells across all processors
!   integer :: numl                           ! total number of landunits across all processors
!   integer :: numc                           ! total number of columns across all processors
!   integer :: nump                           ! total number of pfts across all processors
    integer :: ier                            ! error status
    integer :: nstep                          ! time step index
!-----------------------------------------------------------------------

#if (defined TIMING_BARRIERS)
  ! call t_startf ('sync_clmrtm')
  ! call mpi_barrier (mpicom, ier)
  ! call t_stopf ('sync_clmrtm')
#endif

   ! Assign local pointers to derived type members (gridcell-level)

  ! ixy           => clm3%g%ixy
  ! jxy           => clm3%g%jxy

  ! forc_rain     => clm3%g%a2lf%forc_rain
  ! forc_snow     => clm3%g%a2lf%forc_snow

   ! Assign local pointers to derived type members (column-level)

  ! cgridcell     => clm3%g%l%c%gridcell
  ! wtgcell       => clm3%g%l%c%wtgcell
  ! qflx_qrgwl    => clm3%g%l%c%cwf%qflx_qrgwl
  ! qflx_drain    => clm3%g%l%c%cwf%qflx_drain
  ! qflx_evap_tot => clm3%g%l%c%cwf%pwf_a%qflx_evap_tot
  ! qflx_surf     => clm3%g%l%c%cwf%qflx_surf

    ! Determine subgrid bounds for this processor

  ! call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
  ! call get_proc_global(numg, numl, numc, nump)

    ! Make gridded representation of runoff
    ! total surface runoff = surface runoff on soils + runoff on glaciers,  wetlands, lakes (P-E)

    rtmin(:,:) = 0.

  ! do g = begg,endg
  !    rtmin(2,g) = forc_rain(g) + forc_snow(g)
  ! end do
  ! do c = begc, endc
  !    g = cgridcell(c)
  !    rtmin(1,g) = rtmin(1,g) + (qflx_surf(c) + qflx_qrgwl(c) + qflx_drain(c)) * wtgcell(c)
  !    rtmin(3,g) = rtmin(3,g) + qflx_evap_tot(c) * wtgcell(c)
  ! end do

    do g = 1, numgrid
       rtmin(1,g) = rnof(g)
       rtmin(2,g) = prc(g) + prl(g)
       rtmin(3,g) = fevpa(g)
    end do

    ! Average fluxes for RTM calculation if appropriate
    ! Note: in spmd mode it is assumed that each mpi process has
    ! all the values of totruninxy, precxy and evapxy

    ! RTM input averaging is not done

    if (rtm_nsteps <= 1) then

#if (defined SPMD)
       call allgather_data(rtmin, rtmin_glob, clmlevel='gridcell')
#else
       rtmin_glob => rtmin
#endif
       totruninxy(:,:) = 0.
       precxy(:,:) = 0.
       evapxy(:,:) = 0.
     ! do g = 1,numg
       do g = 1,numgrid_glob
        ! i = ixy(g)
        ! j = jxy(g)
          i = gxmap_glob(g)
          j = gymap_glob(g)
          totruninxy(i,j) = rtmin_glob(1,g)
          precxy(i,j)     = rtmin_glob(2,g)
          evapxy(i,j)     = rtmin_glob(3,g)
       end do
       delt_rtm = get_step_size()
       do_rtm = .true.

    end if

    ! RTM input averaging is done

    if (rtm_nsteps > 1) then

     ! do g = begg,endg
       do g = 1,numgrid
          rtmin_ave(1,g) = rtmin_ave(1,g) + rtmin(1,g)
          rtmin_ave(2,g) = rtmin_ave(2,g) + rtmin(2,g)
          rtmin_ave(3,g) = rtmin_ave(3,g) + rtmin(3,g)
       end do
       ncount_rtm = ncount_rtm + 1
       nstep = get_nstep()

       if ((mod(nstep,rtm_nsteps)==0) .and. (nstep>1)) then
#if (defined SPMD)
          call allgather_data(rtmin_ave, rtmin_glob, clmlevel='gridcell')
#else
          rtmin_glob = rtmin_ave
#endif
          totruninxy(:,:) = 0.
          precxy(:,:)     = 0.
          evapxy(:,:)     = 0.
        ! do g = 1,numg
          do g = 1,numgrid_glob
           ! i = ixy(g)
           ! j = jxy(g)
             i = gxmap_glob(g)
             j = gymap_glob(g)
             totruninxy(i,j) = rtmin_glob(1,g)/ncount_rtm
             precxy(i,j)     = rtmin_glob(2,g)/ncount_rtm
             evapxy(i,j)     = rtmin_glob(3,g)/ncount_rtm
          end do
          delt_rtm = ncount_rtm*get_step_size()   !compute delt for rtm
          ncount_rtm = 0                          !reset counter to 0
        ! do g = begg,endg
          do g = 1,numgrid
             rtmin_ave(1,g) = 0.                  !reset averager
             rtmin_ave(2,g) = 0.                  !reset averager
             rtmin_ave(3,g) = 0.                  !reset averager
          end do
          do_rtm = .true.
       else
          do_rtm = .false.
       endif

    endif

  end subroutine UpdateInput

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtm
!
! !INTERFACE:
  subroutine Rtm
!
! !DESCRIPTION:
! River routing model (based on U. Texas code).
! Input is totrunin\_r.
! Input/output is fluxout, volr.
! Outputs are dvolrdt\_r, flxocn\_r, flxlnd\_r.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine Rtmriverflux in this module
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i, j                        !loop indices
    real(r8) :: buffern(rtmlon)             !temp buffer
    real(r8) :: buffers(rtmlon)             !temp buffer
    real(r8) :: dvolrdt                     !change in storage (m3/s)
    real(r8) :: sumdvolr(rtmlat)            !global sum (m3/s)
    real(r8) :: sumrunof(rtmlat)            !global sum (m3/s)
    real(r8) :: sumdvolr_tot                !global sum (m3/s)
    real(r8) :: sumrunof_tot                !global sum (m3/s)
!-----------------------------------------------------------------------

    ! Determine fluxout at extended points and at southern and northern outer lats

    fluxout(rtmlon+1,rtmlat+1) = fluxout(1,1)
    fluxout(rtmlon+1,0)        = fluxout(1,rtmlat)
    fluxout(0,0)               = fluxout(rtmlon,rtmlat)
    fluxout(0,rtmlat+1)        = fluxout(rtmlon,1)

    do i=1,rtmlon
       fluxout(i,0)        = fluxout(i,rtmlat)
       fluxout(i,rtmlat+1) = fluxout(i,1)
       buffern(i)          = fluxout(i,rtmlat)
       buffers(i)          = fluxout(i,1)
    enddo
    do j=1,rtmlat
       fluxout(0,j)        = fluxout(rtmlon,j)
       fluxout(rtmlon+1,j) = fluxout(1,j)
    enddo
    do i=0,rtmlon+1
       fluxout(i,0)        = buffern(mod(i+rtmlon/2-1,rtmlon)+1)
       fluxout(i,rtmlat+1) = buffers(mod(i+rtmlon/2-1,rtmlon)+1)
    enddo

    ! Determine cell-to-cell transport - calculate sfluxin

!$OMP PARALLEL DO PRIVATE (i,j)
!CSD$ PARALLEL DO PRIVATE (i,j)
    do j=rtmlati,rtmlatf
       do i=rtmloni,rtmlonf
          sfluxin(i,j) = 0.
          if (rdirc(i  ,j-1)==1) sfluxin(i,j) = sfluxin(i,j) + fluxout(i  ,j-1)
          if (rdirc(i-1,j-1)==2) sfluxin(i,j) = sfluxin(i,j) + fluxout(i-1,j-1)
          if (rdirc(i-1,j  )==3) sfluxin(i,j) = sfluxin(i,j) + fluxout(i-1,j  )
          if (rdirc(i-1,j+1)==4) sfluxin(i,j) = sfluxin(i,j) + fluxout(i-1,j+1)
          if (rdirc(i  ,j+1)==5) sfluxin(i,j) = sfluxin(i,j) + fluxout(i  ,j+1)
          if (rdirc(i+1,j+1)==6) sfluxin(i,j) = sfluxin(i,j) + fluxout(i+1,j+1)
          if (rdirc(i+1,j  )==7) sfluxin(i,j) = sfluxin(i,j) + fluxout(i+1,j  )
          if (rdirc(i+1,j-1)==8) sfluxin(i,j) = sfluxin(i,j) + fluxout(i+1,j-1)
       enddo
    enddo
!CSD$ END PARALLEL DO
!$OMP END PARALLEL DO

    ! Loops above and below must remain separate because fluxout is updated below

    sumdvolr(:) = 0.
    sumrunof(:) = 0.
!$OMP PARALLEL DO PRIVATE (i,j,dvolrdt)
!CSD$ PARALLEL DO PRIVATE (i,j,dvolrdt)
    do j = rtmlati,rtmlatf
       do i = rtmloni,rtmlonf

          ! calculate change in cell storage volume change units for
          ! totrunin from kg/m2s==mm/s -> m3/s

          dvolrdt = sfluxin(i,j) - fluxout(i,j) + 0.001*totrunin_r(i,j)*rivarea(i,j)

          ! calculate flux out of a cell:
          ! land: do not permit change in cell storage volume greater than volume present
          ! make up for the difference with an extraction term (eg from aquifers)
          ! ocean: do not permit negative change in cell storage volume,
          ! because at ocean points cell storage volume equals zero
          ! water balance check (in mm/s), convert runinxy from mm/s to m/s (* 1.e-3)
          ! and land model area from km**2 to m**2 (* 1.e6)

          if (mask_r(i,j) == 1) then         ! land points
             volr(i,j)     = volr(i,j) + dvolrdt*delt_rtm
             fluxout(i,j)  = volr(i,j) * effvel/ddist(i,j)
             fluxout(i,j)  = min(fluxout(i,j), volr(i,j) / delt_rtm)
             flxlnd_r(i,j) = fluxout(i,j)
             flxocn_r(i,j) = 0.
          else                               ! ocean points
             flxlnd_r(i,j) = 0.
             flxocn_r(i,j) = dvolrdt
          endif
          sumdvolr(j) = sumdvolr(j) + dvolrdt
          sumrunof(j) = sumrunof(j) + totrunin_r(i,j)*1000.*area_r(i,j)
          dvolrdt_r(i,j) = 1000.*dvolrdt/rivarea(i,j)

       enddo
    enddo
!CSD$ END PARALLEL DO
!$OMP END PARALLEL DO

    ! Global water balance calculation and error check

    sumdvolr_tot = 0.
    sumrunof_tot = 0.
    do j = 1,rtmlat
       sumdvolr_tot = sumdvolr_tot + sumdvolr(j)
       sumrunof_tot = sumrunof_tot + sumrunof(j)
    end do
    if (abs((sumdvolr_tot-sumrunof_tot)/sumrunof_tot) > 0.01) then
       write(6,*) 'RTM Error: sumdvolr= ',sumdvolr_tot,&
            ' not equal to sumrunof= ',sumrunof_tot
       call endrun
    end if

  end subroutine Rtm

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: UpdateGlobal
!
! !INTERFACE:
  subroutine UpdateGlobal(totruninxy, precxy, evapxy)
!
! !DESCRIPTION:
! Update Global quantitities.
! Input is totrunin\_r.
! Input/output is fluxout, volr.
! Outputs are dvolrdt\_r, flxocn\_r, flxlnd\_r.
!
! !USES:
  ! use clm_varsur  , only : numlon, area, landfrac
  ! use time_manager, only : get_nstep, get_curr_date
    use colm_varMod , only : numlon, area, landfrac
    use timemgr     , only : get_nstep, get_curr_date
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(inout) :: totruninxy(lsmlon,lsmlat)  !surface runoff (mm H2O /s)
    real(r8), intent(inout) :: precxy(lsmlon,lsmlat)      !precipitation (mm H2O /s)
    real(r8), intent(inout) :: evapxy(lsmlon,lsmlat)      !evaporation (mm H2O /s)
!
! !CALLED FROM:
! subroutine Rtmriverflux in this module
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
!
! global balance
!
    integer  :: is,js,ir,jr   !indices
    integer  :: yrnew         !year (0, ...)
    integer  :: mon           !month (1, ..., 12)
    integer  :: day           !day of month (1, ..., 31)
    integer  :: ncsec         !seconds of current date
    integer  :: ncdate        !current date
    real(r8) :: prec_sum      !total precipitation (m^3/sec)
    real(r8) :: evap_sum      !total evaporation (m^3/sec)
    real(r8) :: runlnd_sum    !total input runoff on land grid (m^3/sec)
    real(r8) :: runrtm_sum    !total input runoff on rtm grid (m^3/sec)
    real(r8) :: ocnrtm_sum    !total ocean runoff on rtm grid (m^3/sec)
    real(r8) :: volrtm_sum    !total change in storage on rtm (m^3/sec)
    real(r8) :: runrtm(rtmlon,rtmlat) ! input runoff on rtm grid (m**3/s)
    real(r8) :: volrtm(rtmlon,rtmlat) ! change in storage (m**3/s)

    character(len=*),parameter :: F40="('(diag) ',a17,'    date  ', &
       & '   prec        evap        runoff(lnd)   runoff(rtm) dvoldt(rtm) runoff-ocn(rtm)  (m^3/sec)')"
    character(len=*),parameter :: F41="('(diag) ',a17,'   nstep  ', &
       & '   prec        evap        runoff(lnd)   runoff(rtm) dvoldt(rtm) runoff-ocn(rtm)  (m^3/sec)')"
    character(len=*),parameter :: F21="('(diag) ',a17,' ----------------------', &
       & 7('----------'))"
    character(len=*),parameter :: F22="('(diag) ',a17,i8,6(d13.4))"
!-----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE (jr,ir)
!CSD$ PARALLEL DO PRIVATE (jr,ir)
    do jr = rtmlati,rtmlatf
       do ir = rtmloni,rtmlonf
          runrtm(ir,jr) = totrunin_r(ir,jr)*1000.*area_r(ir,jr)
          volrtm(ir,jr) = dvolrdt_r(ir,jr)*1000.*area_r(ir,jr)
       end do
    end do
!CSD$ END PARALLEL DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE (js,is)
!CSD$ PARALLEL DO PRIVATE (js,is)
    do js = 1,lsmlat
       do is = 1,numlon(js)
          precxy(is,js) = precxy(is,js)*area(is,js)*1000.*landfrac(is,js)
          evapxy(is,js) = evapxy(is,js)*area(is,js)*1000.*landfrac(is,js)
          totruninxy(is,js) = totruninxy(is,js)*area(is,js)*1000.*landfrac(is,js)
       end do
    end do
!CSD$ END PARALLEL DO
!$OMP END PARALLEL DO

    prec_sum   = sum(precxy)
    evap_sum   = sum(evapxy)
    runlnd_sum = sum(totruninxy)
    runrtm_sum = sum(runrtm)
    volrtm_sum = sum(volrtm)
    ocnrtm_sum = sum(flxocn_r)

    prec_global   = prec_global   + prec_sum
    evap_global   = evap_global   + evap_sum
    runlnd_global = runlnd_global + runlnd_sum
    runrtm_global = runrtm_global + runrtm_sum
    volrtm_global = volrtm_global + volrtm_sum
    ocnrtm_global = ocnrtm_global + ocnrtm_sum

    ncount_global = ncount_global + 1

    ! Print out diagnostics if appropriate

#ifdef MYBUG
    write(6,*)
    write(6,F41)'water inst   '
    write(6,F21)'water inst   '
    write(6,F22)'water inst   ',get_nstep(), prec_sum, evap_sum, &
         runlnd_sum, runrtm_sum, volrtm_sum, ocnrtm_sum
    write(6,*)
#endif

    call get_curr_date(yrnew, mon, day, ncsec)
    ncdate = yrnew*10000 + mon*100 + day
    if (yrnew /= yrold) then
       prec_global   = prec_global/ncount_global
       evap_global   = evap_global/ncount_global
       runlnd_global = runlnd_global/ncount_global
       runrtm_global = runrtm_global/ncount_global
       volrtm_global = volrtm_global/ncount_global
       ocnrtm_global = ocnrtm_global/ncount_global
       ncount_global = 0
       write(6,*)
       write(6,F40)'water tavg   '
       write(6,F21)'water tavg   '
       write(6,F22)'water tavg   ',ncdate, prec_global, evap_global,&
            runlnd_global, runrtm_global, volrtm_global, ocnrtm_global
       write(6,*)
    endif
    yrold = yrnew

  end subroutine UpdateGlobal

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restart_rtm
!
! !INTERFACE:
  subroutine restart_rtm (nio, flag)
!
! !DESCRIPTION:
! Read/write RTM restart data.
!
! !USES:
    use precision      , only : r8
#if (defined SPMD)
    use spmd           , only : masterproc => p_master, mpicom => p_comm, MPI_INTEGER, MPI_REAL8
    use spmdGathScatMod, only : scatter_data_from_master, gather_data_to_master
#else
    use spmd           , only : masterproc => p_master
#endif
    use runoffMod      , only : runoff
    use colm_varMod    , only : numgrid_glob
    use spmd_decomp    , only : ggmap_glob
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: nio             ! restart unit
    character(len=*), intent(in) :: flag   ! 'read' or 'write'
!
! !CALLED FROM:
! subroutine restart in module restFileMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!   integer :: numg             ! total number of gridcells across all processors
!   integer :: numl             ! total number of landunits across all processors
!   integer :: numc             ! total number of columns across all processors
!   integer :: nump             ! total number of pfts across all processors
    integer :: ier              ! error status
    real(r8), pointer :: rtmin_ave_glob_dc(:,:)  ! temporary
    real(r8), pointer :: rtmin_ave_glob_sn(:,:)  ! temporary
    integer :: g
!-----------------------------------------------------------------------

    ! Allocate dynamic memory

  ! call get_proc_global(numg, numl, numc, nump)

    allocate(rtmin_ave_glob_sn(3,numgrid_glob), rtmin_ave_glob_dc(3,numgrid_glob), stat=ier)
    if (ier /= 0) then
       write (6,*) 'restart_rtm error : allocation error'
       call endrun()
    end if

    ! Read RTM restart - must put in logic to support restart_id 6 for backwards compatibility

    if (flag == 'read') then
       if (masterproc) then
          read(nio) volr
          read(nio) fluxout
          read(nio) ncount_rtm
          read(nio) rtmin_ave_glob_sn
          read(nio) ncount_global, yrold, prec_global, evap_global, &
                    runlnd_global, runrtm_global, volrtm_global, ocnrtm_global
          read(nio) runoff%ocn
          read(nio) runoff%lnd

        ! call map_sn2dc(rtmin_ave_glob_sn, rtmin_ave_glob_dc, lb1=1, ub1=3, type1d=nameg)
          do g = 1, numgrid_glob
             rtmin_ave_glob_dc(:,g) = rtmin_ave_glob_sn(:,ggmap_glob(g))
          end do
       endif

#if (defined SPMD)
       call mpi_bcast(ncount_rtm, 1               , MPI_INTEGER, 0, mpicom, ier)
       call mpi_bcast(volr      , size(volr)      , MPI_REAL8  , 0, mpicom, ier)
       call mpi_bcast(fluxout   , size(fluxout)   , MPI_REAL8  , 0, mpicom, ier)
       call mpi_bcast(runoff%ocn, size(runoff%ocn), MPI_REAL8  , 0, mpicom, ier)
       call mpi_bcast(runoff%lnd, size(runoff%lnd), MPI_REAL8  , 0, mpicom, ier)
       call scatter_data_from_master (rtmin_ave, rtmin_ave_glob_dc, clmlevel='gridcell')
#else
       rtmin_ave(:,:) = rtmin_ave_glob_dc(:,:)
#endif
    endif

    ! Write RTM restart

    if (flag == 'write') then
#if (defined SPMD)
       call gather_data_to_master (rtmin_ave, rtmin_ave_glob_dc, clmlevel='gridcell')
#else
       rtmin_ave_glob_dc(:,:) = rtmin_ave(:,:)
#endif

       if (masterproc) then
        ! call map_dc2sn(rtmin_ave_glob_dc, rtmin_ave_glob_sn, lb1=1, ub1=3, type1d=nameg)
          do g = 1, numgrid_glob
             rtmin_ave_glob_sn(:,ggmap_glob(g)) = rtmin_ave_glob_dc(:,g)
          end do

          write(nio) volr
          write(nio) fluxout
          write(nio) ncount_rtm
          write(nio) rtmin_ave_glob_sn
          write(nio) ncount_global, yrold, prec_global, evap_global, &
                     runlnd_global, runrtm_global, volrtm_global, ocnrtm_global
          write(nio) runoff%ocn
          write(nio) runoff%lnd
       endif
    endif

    ! Deallocate dynamic memory

    deallocate(rtmin_ave_glob_sn, rtmin_ave_glob_dc)

  end subroutine restart_rtm

#endif

end module RtmMod
