#define PFT
#define DGVM
#undef  DyN

module paramodel

!----------------------------------------------------------------------
! Define the dimension of model array
!----------------------------------------------------------------------

      integer nl_soil       ! number of soil layers
      integer maxsnl        ! max number of snow layers
      integer nfcon_col     ! number of time constant variables
      integer nfcon_pft     ! number of time constant variables
      integer nftune        ! number of clm tunable constants
      integer nfvar_col     ! number of time varying variables
      integer nfvar_pft     ! number of time varying variables
      integer nforc         ! number of forcing variables
      integer nfldv_col     ! number of output fluxes
      integer nfldv_pft     ! number of output fluxes
      integer nflai         ! number of leaf time varying variables
      integer nfldv         ! number of output fluxes (2D form)
      integer nflux         ! number of LSM flux variables
      integer maxpatch      ! maximum number of patches in model grid
      integer nlandcateg    ! number of land cover categories
      integer nsoilcateg    ! number of soil texture categories
      integer oceancateg    ! land cover category of ocean type
      integer grasscateg    ! land cover category of grass type

#ifdef RTM
      integer nfldv_rtm     ! number of RTM  fldv variables (2D form)
      integer nflux_rtm     ! number of RTM  flux variables
#endif

#ifdef DGVM
      integer numpft           ! number of PFTs (excluding soil)
      integer numpft_nat       ! number of natural PFTs (excluding soil & crop)
      integer nfldv_dgvm       ! number of output fluxes from DGVM
      integer nflux_dgvm       ! number of DGVM fluxes variables (2D & 3D)
      integer nfldv_dgvm_accum ! number of output fluxes from DGVM (accumulating)
      integer nflux_dgvm_accum ! number of DGVM fluxes variables (2D & 3D) (accumulating)
#endif

      parameter(nl_soil    = 10)
      parameter(maxsnl     = -5)
      parameter(nfcon_col  = 84)

#ifndef DGVM

      parameter(nfcon_pft  = 35)
      parameter(nfvar_col  = 81)
      parameter(nfvar_pft  = 45)

#else

      parameter(numpft     = 16)
      parameter(numpft_nat = 14)

      parameter(nfldv_dgvm_accum = 34)
      parameter(nflux_dgvm_accum = 21)

#ifndef DyN
      parameter(nfcon_pft  = 34+40)
      parameter(nfvar_col  = 81+3)
      parameter(nfvar_pft  = 45+68)
      parameter(nfldv_dgvm = 255)
      parameter(nflux_dgvm = 47)
#else 
      parameter(nfcon_pft  = 34+40+2)
      parameter(nfvar_col  = 81+3+6)
      parameter(nfvar_pft  = 45+68+17)
      parameter(nfldv_dgvm = 255+50)
      parameter(nflux_dgvm = 58)
#endif

#endif

      parameter(nftune     =  14)
      parameter(nforc      =  18)
      parameter(nflai      =   4)

#ifndef CMIP
      parameter(nfldv_col  =  46)
      parameter(nfldv_pft  =  47)
      parameter(nfldv      =  93)
      parameter(nflux      =  66)
#else
      parameter(nfldv_col  =  63)
      parameter(nfldv_pft  =  48)
      parameter(nfldv      = 121)
      parameter(nflux      =  85)
#endif

#ifdef RTM
      parameter(nfldv_rtm  =  2)
      parameter(nflux_rtm  =  2)
#endif

#if(defined PFT)
      parameter(grasscateg = 13)
      parameter(oceancateg = 22)
      parameter(nlandcateg = 22)
#elif(defined USGS)
      parameter(grasscateg = 7)
      parameter(oceancateg = 25)
      parameter(nlandcateg = 25)
#elif(defined OGE)
      parameter(grasscateg = 2)
      parameter(oceancateg = 15)
      parameter(nlandcateg = 96)
#elif(defined BATS)
      parameter(grasscateg = 2)
      parameter(oceancateg = 15)
      parameter(nlandcateg = 19)
#elif(defined IGBP)
      parameter(grasscateg = 10)
      parameter(oceancateg = 18)
      parameter(nlandcateg = 18)
#elif(defined SIB2)
      parameter(grasscateg = 9)
      parameter(oceancateg = 12)
      parameter(nlandcateg = 12)
#endif

      parameter(maxpatch   = nlandcateg)
      parameter(nsoilcateg = 17)

end module paramodel
