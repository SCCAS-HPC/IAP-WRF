#include <define.h>

#ifdef CPL7

module colm_cpl7

   use precision
   use spmd
   use spmd_decomp
   use forcedata
   use colm_cplMod
   use colm_rtmMod
   use colm_ioMod
   use colm_varMod
   use colm_varctl
   use landuse
   use nchistMod
   use shr_msg_mod        ! CSM message passing routines and variables 
   use shr_sys_mod        ! CSM system module

   implicit none

   real(r8)            :: nextsw_cday 

!  Unit Numbers
   integer, public     :: iulog = 6        ! "stdout" log file unit number, default is 6

   logical, save       :: noland    = .false.
   logical, save       :: downscale = .false.

   type domain_type
      integer          :: ns         ! global size of domain
      integer          :: ni,nj      ! global axis if 2d (nj=1 if unstructured)
      integer          :: nbeg,nend  ! local beg/end indices
   !  character(len=8) :: clmlevel   ! grid type
   !  logical          :: decomped   ! decomposed locally or global copy
   !  logical          :: regional   ! regional or global grid
   !  logical          :: areaset    ! has area been set 
      integer ,pointer :: mask(:)    ! land mask: 1 = land, 0 = ocean
      real(r8),pointer :: frac(:)    ! fractional land
   !  real(r8),pointer :: topo(:)    ! topography
      real(r8),pointer :: latc(:)    ! latitude of grid cell (deg)
      real(r8),pointer :: lonc(:)    ! longitude of grid cell (deg)
      real(r8),pointer :: area(:)    ! grid cell area (km**2)
      real(r8),pointer :: asca(:)    ! area scaling from CESM driver
   !  character*16     :: set        ! flag to check if domain is set 
   !  integer ,pointer :: glcmask(:) ! glc mask: 1=sfc mass balance required by GLC component
                                     !           0=SMB not required (default)
   end type domain_type

   type(domain_type),public :: adomain, ldomain

   !---global information on each pe
   type processor_type
   !  integer :: nclumps          ! number of clumps for processor_type iam 
   !  integer,pointer :: cid(:)   ! clump indices
   !  integer :: ncells           ! number of gridcells in proc
   !  integer :: nlunits          ! number of landunits in proc
   !  integer :: ncols            ! number of columns in proc
   !  integer :: npfts            ! number of pfts in proc
      integer :: begg, endg       ! beginning and ending gridcell index
   !  integer :: begl, endl       ! beginning and ending landunit index
   !  integer :: begc, endc       ! beginning and ending column index
   !  integer :: begp, endp       ! beginning and ending pft index
      integer :: abegg,aendg      ! beginning and ending atm gridcell index
   end type processor_type

   type(processor_type),public :: procinfo

   type decomp_type
      integer,pointer :: glo2gdc(:)    ! 1d glo to 1d gdc
      integer,pointer :: gdc2glo(:)    ! 1d gdc to 1d glo
   end type decomp_type

   type(decomp_type),public,target :: ldecomp
   type(decomp_type),public,target :: adecomp

!----------------------------------------------------
! atmosphere -> land variables structure
!----------------------------------------------------
   type atm2lnd_type
      real(r8), pointer :: forc_t(:)       !atmospheric temperature (Kelvin)
      real(r8), pointer :: forc_u(:)       !atm wind speed, east direction (m/s)
      real(r8), pointer :: forc_v(:)       !atm wind speed, north direction (m/s)
      real(r8), pointer :: forc_wind(:)    !atmospheric wind speed   
      real(r8), pointer :: forc_q(:)       !atmospheric specific humidity (kg/kg)
      real(r8), pointer :: forc_hgt(:)     !atmospheric reference height (m)
      real(r8), pointer :: forc_hgt_u(:)   !obs height of wind [m] (new)
      real(r8), pointer :: forc_hgt_t(:)   !obs height of temperature [m] (new)
      real(r8), pointer :: forc_hgt_q(:)   !obs height of humidity [m] (new)
      real(r8), pointer :: forc_pbot(:)    !atmospheric pressure (Pa)
      real(r8), pointer :: forc_th(:)      !atm potential temperature (Kelvin)
      real(r8), pointer :: forc_vp(:)      !atmospheric vapor pressure (Pa) 
      real(r8), pointer :: forc_rho(:)     !density (kg/m**3)
      real(r8), pointer :: forc_rh(:)      !atmospheric relative humidity (%)
      real(r8), pointer :: forc_psrf(:)    !surface pressure (Pa)
      real(r8), pointer :: forc_pco2(:)    !CO2 partial pressure (Pa)
      real(r8), pointer :: forc_lwrad(:)   !downwrd IR longwave radiation (W/m**2)
      real(r8), pointer :: forc_solad(:,:) !direct beam radiation (numrad) (vis=forc_sols , nir=forc_soll )
      real(r8), pointer :: forc_solai(:,:) !diffuse radiation (numrad) (vis=forc_solsd, nir=forc_solld)
      real(r8), pointer :: forc_solar(:)   !incident solar radiation
      real(r8), pointer :: forc_rain(:)    !rain rate [mm/s]
      real(r8), pointer :: forc_snow(:)    !snow rate [mm/s]
      real(r8), pointer :: forc_ndep(:)    !nitrogen deposition rate (gN/m2/s)
      real(r8), pointer :: rainf(:)        !ALMA rain+snow [mm/s]
#ifdef C13
      real(r8), pointer :: forc_pc13o2(:)  !C13O2 partial pressure (Pa)
#endif
      real(r8), pointer :: forc_po2(:)     !O2 partial pressure (Pa)
      real(r8), pointer :: forc_aer(:,:)   ! aerosol deposition array
   end type atm2lnd_type

   real(r8), pointer :: clm_a2l_rainc(:)
   real(r8), pointer :: clm_a2l_rainl(:)
   real(r8), pointer :: clm_a2l_snowc(:)
   real(r8), pointer :: clm_a2l_snowl(:)

!----------------------------------------------------
! land -> atmosphere variables structure
!----------------------------------------------------
   type lnd2atm_type
      real(r8), pointer :: t_rad(:)        !radiative temperature (Kelvin)
      real(r8), pointer :: t_ref2m(:)      !2m surface air temperature (Kelvin)
      real(r8), pointer :: q_ref2m(:)      !2m surface specific humidity (kg/kg)
      real(r8), pointer :: u_ref10m(:)     !10m surface wind speed (m/sec)
      real(r8), pointer :: h2osno(:)       !snow water (mm H2O)
      real(r8), pointer :: albd(:,:)       !(numrad) surface albedo (direct)
      real(r8), pointer :: albi(:,:)       !(numrad) surface albedo (diffuse)
      real(r8), pointer :: taux(:)         !wind stress: e-w (kg/m/s**2)
      real(r8), pointer :: tauy(:)         !wind stress: n-s (kg/m/s**2)
      real(r8), pointer :: eflx_lh_tot(:)  !total latent HF (W/m**2)  [+ to atm]
      real(r8), pointer :: eflx_sh_tot(:)  !total sensible HF (W/m**2) [+ to atm]
      real(r8), pointer :: eflx_lwrad_out(:) !IR (longwave) radiation (W/m**2)
      real(r8), pointer :: qflx_evap_tot(:)!qflx_evap_soi + qflx_evap_can + qflx_tran_veg
      real(r8), pointer :: fsa(:)          !solar rad absorbed (total) (W/m**2)
      real(r8), pointer :: nee(:)          !net CO2 flux (kg CO2/m**2/s) [+ to atm]
      real(r8), pointer :: ram1(:)         !aerodynamical resistance (s/m)
      real(r8), pointer :: fv(:)           !friction velocity (m/s) (for dust model)
      real(r8), pointer :: flxdst(:,:)       !dust flux (size bins)
      real(r8), pointer :: ddvel(:,:)        !dry deposition velocities
      real(r8), pointer :: flxvoc(:,:)       ! VOC flux (size bins)
   end type lnd2atm_type
   
!  type(atm2lnd_type),public,target :: atm_a2l      ! a2l fields on atm grid
!  type(lnd2atm_type),public,target :: atm_l2a      ! l2a fields on atm grid

   type(atm2lnd_type),public,target :: clm_a2l      ! a2l fields on clm grid
   type(lnd2atm_type),public,target :: clm_l2a      ! l2a fields on clm grid

   public :: init_atm2lnd_type
   public :: init_lnd2atm_type

   integer, parameter :: numrad = 2                 ! numrad wavebands (1=vis, 2=nir)

contains

   subroutine control_setNL (NLFile)

      character(len=*), intent(IN) :: NLFile ! Namelist filename

      NLFilename = NLFile

   end subroutine control_setNL

   subroutine build_procinfo

      integer i 

      procinfo%endg = 0

      do i = 0, p_iam
         procinfo%endg = procinfo%endg + numgrid_proc(i)
      end do

      procinfo%begg = procinfo%endg - numgrid_proc(p_iam) + 1

      procinfo%abegg = procinfo%begg
      procinfo%aendg = procinfo%endg

   end subroutine build_procinfo

   subroutine build_domaininfo

      integer g, i, x, y, begg, endg

      call get_proc_bounds_atm(begg,endg)

      adomain%ni = lon_points
      adomain%nj = lat_points
      adomain%ns = lon_points*lat_points

      ldomain%ni = lon_points
      ldomain%nj = lat_points
      ldomain%ns = lon_points*lat_points

      allocate(adomain%lonc(begg:endg))
      allocate(adomain%latc(begg:endg))
      allocate(adomain%area(begg:endg))
      allocate(adomain%mask(begg:endg))
      allocate(adomain%frac(begg:endg))
      allocate(adomain%asca(begg:endg))

      allocate(ldomain%lonc(begg:endg))
      allocate(ldomain%latc(begg:endg))
      allocate(ldomain%area(begg:endg))
      allocate(ldomain%mask(begg:endg))
      allocate(ldomain%frac(begg:endg))
      allocate(ldomain%asca(begg:endg))

      do g = begg, endg
         i = g-begg+1
         x = gxmap(i)
         y = gymap(i)

         adomain%lonc(g) = longxy(x,y)
         adomain%latc(g) = latixy(x,y)
         adomain%area(g) = area(x,y)
         adomain%mask(g) = landmask(x,y)
         adomain%frac(g) = landfrac(x,y)
         adomain%asca(g) = 1.0

         ldomain%lonc(g) = longxy(x,y)
         ldomain%latc(g) = latixy(x,y)
         ldomain%area(g) = area(x,y)
         ldomain%mask(g) = landmask(x,y)
         ldomain%frac(g) = landfrac(x,y)
         ldomain%asca(g) = 1.0
      end do

   end subroutine build_domaininfo

   subroutine get_proc_bounds_atm(begg, endg)

      integer, intent(out) :: begg
      integer, intent(out) :: endg

      begg = procinfo%abegg
      endg = procinfo%aendg

   end subroutine get_proc_bounds_atm

   subroutine get_proc_bounds(begg, endg)

      integer, intent(out) :: begg
      integer, intent(out) :: endg

      begg = procinfo%begg
      endg = procinfo%endg

   end subroutine get_proc_bounds

   subroutine build_decompinfo

      integer i, j, n, ag, an, numg, ni, nj
      integer lsize, gsize, beg, end
      integer pid

      ni = lon_points
      nj = lat_points

      numg = sum(numgrid_proc)

      allocate(adecomp%gdc2glo(numg), adecomp%glo2gdc(ni*nj))

      adecomp%gdc2glo(:)  = 0
      adecomp%glo2gdc(:)  = 0

      allocate(ldecomp%gdc2glo(numg), ldecomp%glo2gdc(ni*nj))

      ldecomp%gdc2glo(:)  = 0
      ldecomp%glo2gdc(:)  = 0

      ag = 0
      do pid = 0, p_nprocs-1
         do j = 1, nj
         do i = 1, ni
            an = (j-1)*ni + i
            if(gmask(i,j).eq.pid) then
               ag = ag + 1
               adecomp%gdc2glo(ag) = an
               adecomp%glo2gdc(an) = ag
               ldecomp%gdc2glo(ag) = an
               ldecomp%glo2gdc(an) = ag
            end if
         end do
         end do
      end do

!     ! set gsMap_atm_gdc2glo
!     call get_proc_bounds_atm(beg, end)
!     allocate(gindex(beg:end))
!     do n = beg,end
!        gindex(n) = adecomp%gdc2glo(n)
!     enddo
!     lsize = end-beg+1
!     gsize = ni * nj
!     call mct_gsMap_init(gsMap_atm_gdc2glo, gindex, mpicom, comp_id, lsize, gsize )
!     deallocate(gindex)

!     ! set gsMap_lnd_gdc2glo
!     call get_proc_bounds(beg, end)
!     allocate(gindex(beg:end))
!     do n = beg,end
!        gindex(n) = ldecomp%gdc2glo(n)
!     enddo
!     lsize = end-beg+1
!     gsize = ni * nj
!     call mct_gsMap_init(gsMap_lnd_gdc2glo, gindex, mpicom, comp_id, lsize, gsize )
!     deallocate(gindex)

   end subroutine build_decompinfo

   subroutine init_atm2lnd_type(beg, end, a2l)

      use nanMod, only: nan
      implicit none

      integer, intent(in) :: beg, end
      type (atm2lnd_type), intent(inout):: a2l
      real(r8) :: ival   ! initial value

      allocate(a2l%forc_t(beg:end))
      allocate(a2l%forc_u(beg:end))
      allocate(a2l%forc_v(beg:end))
      allocate(a2l%forc_wind(beg:end))
      allocate(a2l%forc_q(beg:end))
      allocate(a2l%forc_rh(beg:end))
      allocate(a2l%forc_hgt(beg:end))
      allocate(a2l%forc_hgt_u(beg:end))
      allocate(a2l%forc_hgt_t(beg:end))
      allocate(a2l%forc_hgt_q(beg:end))
      allocate(a2l%forc_pbot(beg:end))
      allocate(a2l%forc_th(beg:end))
      allocate(a2l%forc_vp(beg:end))
      allocate(a2l%forc_rho(beg:end))
      allocate(a2l%forc_psrf(beg:end))
      allocate(a2l%forc_pco2(beg:end))
      allocate(a2l%forc_lwrad(beg:end))
      allocate(a2l%forc_solad(beg:end,numrad))
      allocate(a2l%forc_solai(beg:end,numrad))
      allocate(a2l%forc_solar(beg:end))
      allocate(a2l%forc_rain(beg:end))
      allocate(a2l%forc_snow(beg:end))
      allocate(a2l%forc_ndep(beg:end))
      allocate(a2l%rainf(beg:end))
#if (defined C13)
      allocate(a2l%forc_pc13o2(beg:end))
#endif
      allocate(a2l%forc_po2(beg:end))
      allocate(a2l%forc_aer(beg:end,14))

      allocate(clm_a2l_rainc(beg:end))
      allocate(clm_a2l_rainl(beg:end))
      allocate(clm_a2l_snowc(beg:end))
      allocate(clm_a2l_snowl(beg:end))

      ! ival = nan      ! causes core dump in map_maparray, tcx fix
      ival = 0.0_r8

      a2l%forc_t(beg:end) = ival
      a2l%forc_u(beg:end) = ival
      a2l%forc_v(beg:end) = ival
      a2l%forc_wind(beg:end) = ival
      a2l%forc_q(beg:end) = ival
      a2l%forc_rh(beg:end) = ival
      a2l%forc_hgt(beg:end) = ival
      a2l%forc_hgt_u(beg:end) = ival
      a2l%forc_hgt_t(beg:end) = ival
      a2l%forc_hgt_q(beg:end) = ival
      a2l%forc_pbot(beg:end) = ival
      a2l%forc_th(beg:end) = ival
      a2l%forc_vp(beg:end) = ival
      a2l%forc_rho(beg:end) = ival
      a2l%forc_psrf(beg:end) = ival
      a2l%forc_pco2(beg:end) = ival
      a2l%forc_lwrad(beg:end) = ival
      a2l%forc_solad(beg:end,1:numrad) = ival
      a2l%forc_solai(beg:end,1:numrad) = ival
      a2l%forc_solar(beg:end) = ival
      a2l%forc_rain(beg:end) = ival
      a2l%forc_snow(beg:end) = ival
      a2l%forc_ndep(beg:end) = ival
      a2l%rainf(beg:end) = nan
#ifdef C13
      a2l%forc_pc13o2(beg:end) = ival
#endif
      a2l%forc_po2(beg:end) = ival
      a2l%forc_aer(beg:end,:) = ival

      clm_a2l_rainc(beg:end) = ival
      clm_a2l_rainl(beg:end) = ival
      clm_a2l_snowc(beg:end) = ival
      clm_a2l_snowl(beg:end) = ival

   end subroutine init_atm2lnd_type

   subroutine init_lnd2atm_type(beg, end, l2a)

      implicit none
      integer, intent(in) :: beg, end
      type (lnd2atm_type), intent(inout):: l2a
      real(r8) :: ival   ! initial value

      allocate(l2a%t_rad(beg:end))
      allocate(l2a%t_ref2m(beg:end))
      allocate(l2a%q_ref2m(beg:end))
      allocate(l2a%u_ref10m(beg:end))
      allocate(l2a%h2osno(beg:end))
      allocate(l2a%albd(beg:end,1:numrad))
      allocate(l2a%albi(beg:end,1:numrad))
      allocate(l2a%taux(beg:end))
      allocate(l2a%tauy(beg:end))
      allocate(l2a%eflx_lwrad_out(beg:end))
      allocate(l2a%eflx_sh_tot(beg:end))
      allocate(l2a%eflx_lh_tot(beg:end))
      allocate(l2a%qflx_evap_tot(beg:end))
      allocate(l2a%fsa(beg:end))
      allocate(l2a%nee(beg:end))
!     allocate(l2a%ram1(beg:end))
!     allocate(l2a%fv(beg:end))
!     allocate(l2a%flxdst(beg:end,1:ndst))
!     allocate(l2a%flxvoc(beg:end,1:nvoc))
!     if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
!        allocate(l2a%ddvel(beg:end,1:n_drydep))
!     end if

      ! ival = nan   ! causes core dump in map_maparray, tcx fix
      ival = 0.0_r8

      l2a%t_rad(beg:end) = ival
      l2a%t_ref2m(beg:end) = ival
      l2a%q_ref2m(beg:end) = ival
      l2a%u_ref10m(beg:end) = ival
      l2a%h2osno(beg:end) = ival
      l2a%albd(beg:end,1:numrad) = ival
      l2a%albi(beg:end,1:numrad) = ival
      l2a%taux(beg:end) = ival
      l2a%tauy(beg:end) = ival
      l2a%eflx_lwrad_out(beg:end) = ival
      l2a%eflx_sh_tot(beg:end) = ival
      l2a%eflx_lh_tot(beg:end) = ival
      l2a%qflx_evap_tot(beg:end) = ival
      l2a%fsa(beg:end) = ival
      l2a%nee(beg:end) = ival
!     l2a%ram1(beg:end) = ival
!     l2a%fv(beg:end) = ival
!     l2a%flxdst(beg:end,1:ndst) = ival
!     l2a%flxvoc(beg:end,1:nvoc) = ival
!     if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
!        l2a%ddvel(beg:end, : ) = ival
!     end if

   end subroutine init_lnd2atm_type

   subroutine build_clm_l2a(init)

      logical, optional, intent(in) :: init   ! at initial time(init==.true.), only set a subset of arguments

      integer begg, endg, g, i

      call get_proc_bounds(begg,endg)

      if (present(init)) then
         do g = begg, endg
            i = g-begg+1

            clm_l2a%h2osno(g)         = lnd_scv(i)
            clm_l2a%t_rad(g)          = lnd_trad(i)

          ! numrad wavebands (1=vis, 2=nir)
            clm_l2a%albd(g,1)         = lnd_avsdr(i)
            clm_l2a%albd(g,2)         = lnd_anidr(i)
            clm_l2a%albi(g,1)         = lnd_avsdf(i)
            clm_l2a%albi(g,2)         = lnd_anidf(i)
         end do
      else
         do g = begg, endg
            i = g-begg+1

            clm_l2a%h2osno(g)         = lnd_scv(i)
            clm_l2a%t_rad(g)          = lnd_trad(i)

          ! numrad wavebands (1=vis, 2=nir)
            clm_l2a%albd(g,1)         = lnd_avsdr(i)
            clm_l2a%albd(g,2)         = lnd_anidr(i)
            clm_l2a%albi(g,1)         = lnd_avsdf(i)
            clm_l2a%albi(g,2)         = lnd_anidf(i)

            clm_l2a%t_ref2m(g)        = lnd_tref(i)
            clm_l2a%q_ref2m(g)        = lnd_qref(i)
            clm_l2a%u_ref10m(g)       = sqrt(lnd_u10m(i)*lnd_u10m(i)+lnd_v10m(i)*lnd_v10m(i))
            clm_l2a%taux(g)           = lnd_taux(i)
            clm_l2a%tauy(g)           = lnd_tauy(i)
            clm_l2a%eflx_lh_tot(g)    = lnd_lhflx(i)
            clm_l2a%eflx_sh_tot(g)    = lnd_shflx(i)
            clm_l2a%eflx_lwrad_out(g) = lnd_lwup(i)
            clm_l2a%qflx_evap_tot(g)  = lnd_qflx(i)
            clm_l2a%fsa(g)            = lnd_swabs(i)
            clm_l2a%nee(g)            = lnd_nee(i)
         !  clm_l2a%ram1(g)           = ???
         !  clm_l2a%fv(g)             = ???
         !  clm_l2a%flxdst(g)         = ???
         !  clm_l2a%ddvel(g)          = ???
         !  clm_l2a%flxvoc(g)         = ???
         end do
      end if
      
   end subroutine build_clm_l2a

   subroutine build_clm_a2l()

      integer begg, endg, g, i
      real(r8) forc_rainc, forc_rainl, forc_snowc, forc_snowl

      call get_proc_bounds(begg,endg)

      do g = begg, endg
         i = g-begg+1

         forcg( 1,i) = clm_a2l%forc_pco2(g)

         forcg( 2,i) = clm_a2l%forc_po2(g)
         forcg( 3,i) = clm_a2l%forc_u(g)
         forcg( 4,i) = clm_a2l%forc_v(g)
         forcg( 5,i) = clm_a2l%forc_t(g)
         forcg( 6,i) = clm_a2l%forc_q(g)

         forc_rainc  = clm_a2l_rainc(g)
         forc_rainl  = clm_a2l_rainl(g)
         forc_snowc  = clm_a2l_snowc(g)
         forc_snowl  = clm_a2l_snowl(g)

         forcg( 7,i) = forc_rainc + forc_snowc
         forcg( 8,i) = forc_rainl + forc_snowl

         forcg( 9,i) = clm_a2l%forc_pbot(g)
         forcg(10,i) = clm_a2l%forc_pbot(g)

         forcg(11,i) = clm_a2l%forc_solad(g,1) ! bufR(g,index_c2l_Faxa_swvdr)  !forc_solsxy  Atm flux  W/m^2
         forcg(12,i) = clm_a2l%forc_solad(g,2) ! bufR(g,index_c2l_Faxa_swndr)  !forc_sollxy  Atm flux  W/m^2
         forcg(13,i) = clm_a2l%forc_solai(g,1) ! bufR(g,index_c2l_Faxa_swvdf)  !forc_solsdxy Atm flux  W/m^2
         forcg(14,i) = clm_a2l%forc_solai(g,2) ! bufR(g,index_c2l_Faxa_swndf)  !forc_solldxy Atm flux  W/m^2

         forcg(15,i) = clm_a2l%forc_lwrad(g)
         forcg(16,i) = clm_a2l%forc_hgt_u(g)
         forcg(17,i) = clm_a2l%forc_hgt_t(g)
         forcg(18,i) = clm_a2l%forc_hgt_q(g)
      end do

   end subroutine build_clm_a2l

   subroutine set_timemgr_init( calendar_in,      start_ymd_in,     start_tod_in, ref_ymd_in,        &
                                ref_tod_in,       stop_ymd_in,      stop_tod_in,  perpetual_run_in,  &
                                perpetual_ymd_in, nelapse_in,       dtime_in )
      use timemgr, only: idate
      use abortutils, only: endrun

      implicit none

      !---------------------------------------------------------------------------------
      ! set time manager startup values
      ! 
      ! Arguments
      character(len=*), optional, intent(IN) :: calendar_in       ! Calendar type
      integer         , optional, intent(IN) :: nelapse_in        ! Number of step (or days) to advance
      integer         , optional, intent(IN) :: start_ymd_in      ! Start date       (YYYYMMDD)
      integer         , optional, intent(IN) :: start_tod_in      ! Start time of day (sec)
      integer         , optional, intent(IN) :: ref_ymd_in        ! Reference date   (YYYYMMDD)
      integer         , optional, intent(IN) :: ref_tod_in        ! Reference time of day (sec)
      integer         , optional, intent(IN) :: stop_ymd_in       ! Stop date        (YYYYMMDD)
      integer         , optional, intent(IN) :: stop_tod_in       ! Stop time of day (sec)
      logical         , optional, intent(IN) :: perpetual_run_in  ! If in perpetual mode or not
      integer         , optional, intent(IN) :: perpetual_ymd_in  ! Perpetual date   (YYYYMMDD)
      integer         , optional, intent(IN) :: dtime_in          ! Time-step (sec)
      !
      character(len=*), parameter :: sub = 'clm::set_timemgr_init'

      integer months(13), yy, mm, dd

!     if ( timemgr_set ) call endrun( sub//":: timemgr_init or timemgr_restart already called" )
!     if (present(calendar_in)      ) calendar         = trim(calendar_in)
!     if (present(start_ymd_in)     ) start_ymd        = start_ymd_in
!     if (present(start_tod_in)     ) start_tod        = start_tod_in
!     if (present(ref_ymd_in)       ) ref_ymd          = ref_ymd_in
!     if (present(ref_tod_in)       ) ref_tod          = ref_tod_in
!     if (present(stop_ymd_in)      ) stop_ymd         = stop_ymd_in
!     if (present(stop_tod_in)      ) stop_tod         = stop_tod_in
!     if (present(perpetual_run_in) )then
!         tm_perp_calendar = perpetual_run_in
!         if ( tm_perp_calendar ) then
!            if ( .not. present(perpetual_ymd_in) .or. perpetual_ymd == uninit_int) &
!                call endrun( sub//":: perpetual_run set but NOT perpetual_ymd" )
!            perpetual_ymd    = perpetual_ymd_in
!         end if
!     end if
!     if (present(nelapse_in)       ) nelapse          = nelapse_in
!     if (present(dtime_in)         ) dtime            = dtime_in

      if (present(calendar_in) .and. present(start_ymd_in) .and. present(start_tod_in)) then
         if (trim(calendar_in)=='NO_LEAP') then
            yy = start_ymd_in/10000
            mm = (start_ymd_in-10000*yy)/100
            dd = (start_ymd_in-10000*yy-100*mm)

            months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)

            idate(1) = yy
            idate(2) = months(mm)+dd
            idate(3) = start_tod_in
         else
            call endrun("Wrong calendar in set_timemgr_init")
         end if
      else
         call endrun("Lack parameters in set_timemgr_init")
      end if

   end subroutine set_timemgr_init

   subroutine set_nextsw_cday( nextsw_cday_in )
   
    ! Set the next radiation calendar day, so that radiation step can be calculated
    ! Arguments
      real(r8), intent(IN) :: nextsw_cday_in ! input calday of next radiation computation
     
      character(len=*), parameter :: sub = 'clm::set_nextsw_cday'
   
      nextsw_cday = nextsw_cday_in
     
   end subroutine set_nextsw_cday

   subroutine update_rad_dtime(doalb)

      logical,intent(in) ::  doalb

    ! CoLM calculate surface radiation every step, so this subroutine is only a placeholder
      return

   end subroutine update_rad_dtime

   subroutine advance_timestep()

      use timemgr, only : TICKTIME, istep

    ! Calendar for NEXT time step
      CALL TICKTIME

      istep = istep + 1

   end subroutine advance_timestep

   subroutine initialize1

    ! Read namelist
      CALL readnml

      CALL colm_var_alloc1

    ! Read surface data, to construct grid info
      CALL readgridat

    ! MPI decomposition
      CALL task_decomp

      call build_procinfo

      call build_domaininfo

      call build_decompinfo

   end subroutine initialize1

   subroutine initialize2

      use timemgr, only : TICKTIME, istep

      integer begg, endg

#ifdef RTM
    ! Initialize RTM model
      CALL colm_rtm_init
#endif

    ! land-atmos flux initialize
      CALL colm_cpl_init

    ! Initialize time-constant and time-varying variables
      CALL colm_var_alloc2

      CALL readinidat

#ifdef VEGDATA
      CALL readvegdat
#endif

      CALL colm_var_dealloc1

    ! Initialize land use module
      CALL landuse_init

    ! advance one step to continue run
      if(nsrest.gt.0) then
         CALL TICKTIME
         istep = istep + 1
      endif

    ! Initialize netCDF history output
      CALL nchist_init

    ! call further initialization subroutines 

      call get_proc_bounds    (begg    , endg)
      call init_atm2lnd_type  (begg    , endg    , clm_a2l)
      call init_lnd2atm_type  (begg    , endg    , clm_l2a)

!     call get_proc_bounds_atm(begg_atm, endg_atm)
!     call init_atm2lnd_type  (begg_atm, endg_atm, atm_a2l)
!     call init_lnd2atm_type  (begg_atm, endg_atm, atm_l2a)

      call build_clm_l2a(init=.true.)

   end subroutine initialize2

   subroutine clm_drv(doalb, nextsw_cday, declinp1, declin, rstwr, nlend, rdate)

      implicit none
      logical,         intent(in) :: doalb       ! true if time for surface albedo calc
      real(r8),        intent(in) :: nextsw_cday ! calendar day for nstep+1
      real(r8),        intent(in) :: declinp1    ! declination angle for next time step
      real(r8),        intent(in) :: declin      ! declination angle for current time step
      logical,         intent(in) :: rstwr       ! true => write restart file this step
      logical,         intent(in) :: nlend       ! true => end of run on this step
      character(len=*),intent(in) :: rdate       ! restart file time stamp for name

!local variables:

      logical :: lwrite                          ! true: write out frequency
      logical :: doalb_lnd                       ! true if time for surface albedo calculation
      logical :: dolai_lnd                       ! true if time for time-varying vegetation paramter
      logical :: dosst_lnd                       ! true if time for update sst/ice/snow

      integer :: i,j,k,l,m                       ! looping indices

!ROUTINE:


! ======================================================================
! begin time stepping loop
! ======================================================================

       ! doalb_lnd is true when the next time step is a radiation time step
         doalb_lnd = .true.
         dolai_lnd = .true.
         dosst_lnd = .false.

         CALL build_clm_a2l
         CALL colm_cpl_a2l

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

       ! Average fluxes over interval if appropriate
       ! Surface states sent to the flux coupler states are not time averaged
         CALL colm_cpl_l2a
         CALL build_clm_l2a

#ifdef MYBUG
         if (p_master) then
            write(6,'(a, a10, i10, a8, i4.4, a8, i3.3, a8, i5.5)' )  &
                    '--------------------CoLM Main Count : Steps & Time : ', &
                    'Nstep=', istep, 'Year=', idate(1), 'Days=', idate(2), 'Secs=', idate(3) 

            call shr_sys_flush(6)
         endif
#endif
 
         CALL writehistdat
 
! ======================================================================
! end of time stepping loop
! ======================================================================

   end subroutine clm_drv

   subroutine finalize()

      CALL colm_cpl_exit
     
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

   end subroutine finalize

end module colm_cpl7

#endif
