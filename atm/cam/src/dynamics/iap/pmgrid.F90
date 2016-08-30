
module pmgrid

!----------------------------------------------------------------------- 
! 
! Purpose: Parameters and variables related to the dynamics grid
! 
! Author: 
! 
! Modified: ZhangHe, Hejuanxiong, Wujiangping
! Reviewed: ZhangHe, 2011-11-19
! Modified: Jiang Jinrong, October 2012, for 2D parallel
!           ZhangHe, 2013-01-30, added numbnd2, mod_gatscat
!-----------------------------------------------------------------------
   use decompmodule, only : decomptype

   implicit none

   public 

   integer, parameter :: plon   = PLON                  ! number of longitudes
   integer, parameter :: plev   = PLEV                  ! number of vertical levels
   integer, parameter :: plat   = PLAT                  ! number of latitudes
   integer, parameter :: plevp  = plev + 1              ! plev + 1
   integer, parameter :: plnlv  = plon*plev             ! Length of multilevel field slice
   integer, parameter :: numbnd2 = 0                    ! no.of latitudes passed N and S of forecast lat
!
!============================ zhh and jjr ==================================
   integer :: beglatdyn  ! beglat for dynamical framework
   integer :: endlatdyn  ! endlat for dynamical framework
   integer :: beglatxydyn  ! beglat for dynamical framework
   integer :: endlatxydyn  ! endlat for dynamical framework
   integer :: beglatxydynex  ! extended beglat for dynamical framework
   integer :: beglatdynex  ! extended beglat for dynamical framework
   integer :: endlatdynex  ! extended endlat for dynamical framework
   integer :: endlatxydynex  ! extended endlat for dynamical framework
   integer :: loc_JB     ! [loc_JB, loc_JE] of a processor related to [JB, JE].
   integer :: loc_JE     ! [loc_JB, loc_JE] of a processor related to [JB, JE].
!==============================================================================
   integer :: beglat     ! beg. index for latitudes owned by a given proc
   integer :: endlat     ! end. index for latitudes owned by a given proc
   integer :: numlats    ! number of latitudes owned by a given proc
   logical :: dyndecomp_set = .false. ! flag indicates dynamics grid has been set for history
!
!jjr added
   integer beglev     ! beg. index for levels owned by a given task
   integer endlev     ! end. index for levels owned by a given task
   integer endlevp1   ! end. index + 1 for levels owned by a given task
   integer endlevp    ! equals endlev, except in last subdomain where equals endlevp1

   integer myid_y     ! subdomain index (0-based) in latitude (y)
   integer myid_z     ! subdomain index (0 based) in level (z)
   integer npr_y      ! number of subdomains in y
   integer npr_z      ! number of subdomains in z
   integer omptotal   ! total number of omp threads in cd_core
   integer ompouter   ! number of outer omp threads in cd_core
   integer ompinner   ! number of inner omp threads in cd_core
   integer :: twod_decomp = 0  ! 1 for multi-2D decomposition with transposes, 0 otherwise
   integer :: mod_transpose = 0  ! method for computing mod_comm transposes
   integer :: mod_geopk = 0  ! method for computing mod_comm geopk
   integer :: mod_gatscat = 0  ! method for computing mod_comm gather/scatters

!  secondary xy decomposition used for remapping

   integer myidxy_x     ! subdomain index (0-based) in longitude (x) (second. decomp.)
   integer myidxy_y     ! subdomain index (0 based) in latitude (y) (second. decomp.)
   integer nprxy_x      ! number of subdomains in x (second. decomp.)
   integer nprxy_y      ! number of subdomains in y (second. decomp.)
   integer beglonxy     ! beg. index for longitudes (second. decomp.)
   integer endlonxy     ! end. index for longitudes (second. decomp.)
   integer beglatxy     ! beg. index for latitudes (second. decomp.)
   integer endlatxy     ! end. index for latitudes (second. decomp.)
!
   integer :: spmd_on = 0 ! 1 for Spmd, 0 for non-Spmd

   type(decomptype), save :: strip2d, strip3dxyz, strip3dxzy,              &
                       strip3dxyzp, strip3dxzyp, strip3zaty,               &
                       strip3zatypt, strip3yatz, strip3yatzp,              &
                       strip3kxyz, strip3kxzy, strip3kxyzp, strip3kxzyp
!

#if ( ! defined SPMD )
   parameter (beglat   = 1)
   parameter (endlat   = plat)
   parameter (numlats  = plat)
!============================ zhh(2007.5.24) ==================================
   parameter (beglatdyn = 1)
   parameter (endlatdyn = plat)
   parameter (beglatdynex = 1)             !zhh
   parameter (endlatdynex = plat)          !zhh
   parameter (loc_JB = 2)
   parameter (loc_JE = plat-1)
!==============================================================================
!
! jjr added 
   parameter (beglev = 1)
   parameter (endlev = plev)
   parameter (endlevp1 = plev+1)
   parameter (endlevp = plev+1)
   parameter (myid_y = 0)
   parameter (myid_z = 0)
   parameter (npr_y = 1)
   parameter (npr_z = 1)
!
! These are needed to pass strict run-time error checking
!
   parameter (myidxy_x=0)     ! subdomain index (0-based) in longitude (x) (second. decomp.)
   parameter (myidxy_y=0)     ! subdomain index (0 based) in latitude (y) (second. decomp.)
   parameter (nprxy_x=1)      ! number of subdomains in x (second. decomp.)
   parameter (nprxy_y=1)      ! number of subdomains in y (second. decomp.)
   parameter (beglonxy=1)     ! beg. index for longitudes (second. decomp.)
   parameter (endlonxy=plon)     ! end. index for longitudes (second. decomp.)
   parameter (beglatxy=1)     ! beg. index for latitudes (second. decomp.)
   parameter (endlatxy=plat)     ! end. index for latitudes (second. decomp.)
#endif

! Staggered grid parameters
! splon and splat may eventually need to become new parameters
! in params.h to define the size of the staggered grid arrays. - gg

   integer, parameter :: splon = plon     ! Number of longitudes on the staggered grid
   integer, parameter :: splat = plat     ! Number of latitudes on the staggered grid
!
end module pmgrid

