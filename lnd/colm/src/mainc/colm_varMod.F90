#include <define.h>

module colm_varMod

   use precision
   use paramodel

   implicit none

!
! constant flux table (fldv), consistant with fluxave.F90
! all idx_*** used with fldv*** variables
!

   integer, parameter :: idx_taux    = 1
   integer, parameter :: idx_tauy    = 2
   integer, parameter :: idx_shflx   = 3
   integer, parameter :: idx_lhflx   = 4
   integer, parameter :: idx_qflx    = 5
   integer, parameter :: idx_sabvsun = 12
   integer, parameter :: idx_sabvsha = 13
   integer, parameter :: idx_sabg    = 14
   integer, parameter :: idx_lwup    = 15
   integer, parameter :: idx_assim   = 18
   integer, parameter :: idx_respc   = 19
   integer, parameter :: idx_fmicr   = 20
   integer, parameter :: idx_avsdr   = 28
   integer, parameter :: idx_avsdf   = 29
   integer, parameter :: idx_anidr   = 30
   integer, parameter :: idx_anidf   = 31
   integer, parameter :: idx_trad    = 34
   integer, parameter :: idx_tref    = 43
   integer, parameter :: idx_qref    = 44
   integer, parameter :: idx_tg      = 81
   integer, parameter :: idx_scv     = 82
   integer, parameter :: idx_fsno    = 84

   integer, parameter :: idx_rnof    = 50
   integer, parameter :: idx_prc     = 89
   integer, parameter :: idx_prl     = 90

   integer, parameter :: idx_u10m    = 45
   integer, parameter :: idx_v10m    = 46
   integer, parameter :: idx_sublim  = 94

   integer, parameter :: idx_mrsos   = 105
   integer, parameter :: idx_tsn     = 110
   integer, parameter :: idx_nsnow   = 111

! 
! constant flux table (flux), consistant with fluxave.F90
! all idx2_*** used with fluxmask*** variables
!
   integer, parameter :: idx2_rnof   = 50
   integer, parameter :: idx2_tg     = 54
   integer, parameter :: idx2_scv    = 55
   integer, parameter :: idx2_fsno   = 57
   integer, parameter :: idx2_mrsos  = 69

! Orbital information needed as input to orbit_parms

  real(r8) :: eccen   ! Earth's eccentricity factor (unitless) (typically 0 to 0.1)

! Orbital information after processed by orbit_params

  real(r8) :: obliqr  ! Earth's obliquity in radians
  real(r8) :: lambm0  ! Mean longitude of perihelion at the vernal equinox (radians)
  real(r8) :: mvelpp  ! Earth's moving vernal equinox longitude of perihelion plus pi (radians)

!
! basic model grid info
!
   integer :: lon_points                      ! number of longitude points on model grid
   integer :: lat_points                      ! number of latitude points on model grid

   integer :: numgrid                         ! local grids number
   integer :: numpatch                        ! local patches number

   integer :: numgrid_glob                    ! global grids number
   integer :: numpatch_glob                   ! global patches number

#if(defined DGVM)
   integer :: numcolumn
   integer :: numcolumn_glob
#endif

! 
! land model grid location info
!
   integer , pointer :: numlon(:)             ! longitude points for each latitude strip
   real(r8), pointer :: latixy(:,:)           ! latitude of grid cell (degrees)
   real(r8), pointer :: longxy(:,:)           ! longitude of grid cell (degrees)
   real(r8), pointer :: area(:,:)             ! grid cell area (km**2)
   real(r8), pointer :: landarea              ! total land area for all gridcells (km^2)
   real(r8), pointer :: lats(:)               ! grid cell latitude, southern edge (degrees)
   real(r8), pointer :: lonw(:,:)             ! grid cell longitude, western edge (degrees)
!
! fractional land and mask
!
   real(r8), pointer :: landfrac(:,:)         ! fractional land
   integer , pointer :: landmask(:,:)         ! land mask: 1 = land. 0 = ocean
!
! patch and grid info
!
   real(r8), pointer :: wxy_patch(:)          ! patch weight

   integer , pointer :: ixy_patch_glob(:)     ! longitude index of patch
   integer , pointer :: jxy_patch_glob(:)     ! latitude index of patch
   real(r8), pointer :: wxy_patch_glob(:)     ! weight of patch

   integer,  pointer :: itypwat_glob(:)

#if(defined DGVM)
   real(r8), pointer :: wxy_column(:)         ! patch weight

   integer , pointer :: ixy_column_glob(:)    ! longitude index of patch
   integer , pointer :: jxy_column_glob(:)    ! latitude index of patch
   real(r8), pointer :: wxy_column_glob(:)    ! weight of patch
#endif
!
! model variables
!
   real(r8), pointer :: forc(:,:)             ! forcing variables
!  real(r8), pointer :: fcon(:,:)             ! time constant variables
!  real(r8), pointer :: fvar(:,:)             ! time varying variables

   real(r8), pointer :: oro(:)                ! ocean(0)/seaice(2)/ flag
   real(r8), pointer :: fldv(:,:)             ! output fluxes in grid average
   integer , pointer :: itypwat(:)            ! land water type

#if(defined DGVM)
   integer,  pointer :: numcolumn_lat(:)      ! number of columnes of grids at lon. strip
   real(r8), pointer :: fcon_col(:,:)         ! time constant variables
   real(r8), pointer :: fvar_col(:,:)         ! time varying variables
   real(r8), pointer :: fldv_col(:,:)         ! output fluxes
   real(r8), pointer :: fcon_pft(:,:)         ! time constant variables
   real(r8), pointer :: fvar_pft(:,:)         ! time varying variables
   real(r8), pointer :: fldv_pft(:,:)         ! output fluxes
   integer,  pointer :: numpatch_lat(:)       ! number of patches of grids at lon. strip
   real(r8), pointer :: fldv_dgvm(:,:)        ! fluxes from LPJ at column level
   real(r8), pointer :: nep_residual(:)       ! annual residual NEP = acfire + aestabc
   logical           :: lnep_adjust = .false.
   real(r8), pointer :: fLitterSoil(:)        ! cflux_litter_soil
   real(r8), pointer :: fLitterAtmos(:)       ! cflux_litter_atmos
#endif

   real(r8)  ftune(nftune)                    ! clm tunable constants

!
! forcing variables
!
   real(r8), pointer :: tair   (:,:)
   real(r8), pointer :: qair   (:,:)
   real(r8), pointer :: pres   (:,:)
   real(r8), pointer :: rainc  (:,:)
   real(r8), pointer :: rainl  (:,:)
   real(r8), pointer :: windu  (:,:)
   real(r8), pointer :: windv  (:,:)
   real(r8), pointer :: dswrf  (:,:)
   real(r8), pointer :: dlwrf  (:,:)
   real(r8), pointer :: tair_z (:,:)
   real(r8), pointer :: qair_z (:,:)
   real(r8), pointer :: wind_z (:,:)

   real(r8), pointer :: lai(:,:)
   real(r8), pointer :: sai(:,:)
   real(r8), pointer :: fveg(:,:)
   real(r8), pointer :: green(:,:)

#if(defined VEGDATA)
   real(r8), pointer :: mlai(:,:)            ! monthly LAI (/12,numpatch/)
   real(r8), pointer :: msai(:,:)            ! monthly SAI (/12,numpatch/)
#endif

   interface colm_var_alloc1
      module procedure colm_var_alloc1
   end interface 

   interface colm_var_alloc2
      module procedure colm_var_alloc2
   end interface 

   interface colm_var_dealloc1
      module procedure colm_var_dealloc1
   end interface 

   interface colm_var_dealloc2
      module procedure colm_var_dealloc2
   end interface 

CONTAINS

   subroutine colm_var_alloc1

      implicit none

! routine:

      allocate (numlon                (lat_points))
      allocate (lats                (lat_points+1))
      allocate (lonw     (lon_points+1,lat_points))
      allocate (area       (lon_points,lat_points))
      allocate (latixy     (lon_points,lat_points))
      allocate (longxy     (lon_points,lat_points))
      allocate (landfrac   (lon_points,lat_points))
      allocate (landmask   (lon_points,lat_points))

      allocate (ixy_patch_glob     (numpatch_glob))
      allocate (jxy_patch_glob     (numpatch_glob))
      allocate (wxy_patch_glob     (numpatch_glob))

#if(defined DGVM)
      allocate (ixy_column_glob   (numcolumn_glob))
      allocate (jxy_column_glob   (numcolumn_glob))
      allocate (wxy_column_glob   (numcolumn_glob))
#endif

      allocate (itypwat_glob      (numcolumn_glob))

   end subroutine colm_var_alloc1

   subroutine colm_var_alloc2

      use paramodel, only: nflai, nforc, nfvar_col, nfvar_pft, nfldv_dgvm, &
                           nfcon_col, nfcon_pft, nfldv_col, nfldv_pft, maxpatch

      implicit none

      allocate (oro                (numcolumn))
      allocate (itypwat            (numcolumn))
      allocate (forc         (nforc,numcolumn))

      allocate (wxy_column         (numcolumn))
      allocate (fcon_col (nfcon_col,numcolumn))
      allocate (fvar_col (nfvar_col,numcolumn))
      allocate (fldv_col (nfldv_col,numcolumn))

      allocate (wxy_patch           (numpatch))
      allocate (fcon_pft  (nfcon_pft,numpatch))
      allocate (fvar_pft  (nfvar_pft,numpatch))
      allocate (fldv_pft  (nfldv_pft,numpatch))

#if(defined VEGDATA)
      allocate (mlai             (12,numpatch))
      allocate (msai             (12,numpatch))
#endif

      allocate (fldv           (nfldv,numgrid))
      allocate (fldv_dgvm (nfldv_dgvm,numgrid))

#ifdef DGVM
      allocate (nep_residual         (numgrid))
      nep_residual(:) = 0.

      allocate (fLitterSoil         (numpatch))
      allocate (fLitterAtmos        (numpatch))
#endif

#ifndef COUP_CSM 
      allocate (tair   (lon_points,lat_points))
      allocate (qair   (lon_points,lat_points))
      allocate (pres   (lon_points,lat_points))
      allocate (rainc  (lon_points,lat_points))
      allocate (rainl  (lon_points,lat_points))
      allocate (windu  (lon_points,lat_points))
      allocate (windv  (lon_points,lat_points))
      allocate (dswrf  (lon_points,lat_points))
      allocate (dlwrf  (lon_points,lat_points))
      allocate (tair_z (lon_points,lat_points))
      allocate (qair_z (lon_points,lat_points))
      allocate (wind_z (lon_points,lat_points))
#endif

      allocate (lai    (lon_points,lat_points))
      allocate (sai    (lon_points,lat_points))
      allocate (fveg   (lon_points,lat_points))
      allocate (green  (lon_points,lat_points))

   end subroutine colm_var_alloc2

   subroutine colm_var_dealloc1

      implicit none

! routine:

      deallocate (ixy_patch_glob)
      deallocate (jxy_patch_glob)
      deallocate (wxy_patch_glob)

#ifdef DGVM
      deallocate (ixy_column_glob)
      deallocate (jxy_column_glob)
      deallocate (wxy_column_glob)
#endif

      deallocate (itypwat_glob)

   end subroutine colm_var_dealloc1

   subroutine colm_var_dealloc2

      implicit none

! routine:

      deallocate (numlon      )
      deallocate (lats        )
      deallocate (lonw        )
      deallocate (latixy      )
      deallocate (longxy      )
      deallocate (area        )
      deallocate (landmask    )
      deallocate (landfrac    )

      deallocate (oro         )
      deallocate (itypwat     )
      deallocate (wxy_patch   )
      deallocate (wxy_column  )
      deallocate (forc        )
      deallocate (fcon_col    )
      deallocate (fvar_col    )
      deallocate (fldv_col    )
      deallocate (fcon_pft    )
      deallocate (fvar_pft    )
      deallocate (fldv_pft    )

#if(defined VEGDATA)
      deallocate (mlai        )
      deallocate (msai        )
#endif

      deallocate (fldv        )
      deallocate (fldv_dgvm   )

#ifdef DGVM
      deallocate (nep_residual)
      deallocate (fLitterSoil )
      deallocate (fLitterAtmos)
#endif

#ifndef COUP_CSM 
      deallocate (tair        )
      deallocate (qair        )
      deallocate (pres        )
      deallocate (rainc       )
      deallocate (rainl       )
      deallocate (windu       )
      deallocate (windv       )
      deallocate (dswrf       )
      deallocate (dlwrf       )
      deallocate (tair_z      )
      deallocate (qair_z      )
      deallocate (wind_z      )
#endif

      deallocate (lai         )
      deallocate (sai         )
      deallocate (fveg        )
      deallocate (green       )

   end subroutine colm_var_dealloc2

end module colm_varMod
