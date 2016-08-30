 subroutine DGVMphenology(ivt,lai_ind,agdd0,tref10,assimn10,tmomin20,&
                          vegclass,raingreen,summergreen,t10min,& 
                          leafon,leafof,dphen,htop,lai,sai,rstfac)
!Called once per day.

    use precision
    implicit none
                          !INTENT IN VARIABLES
    integer , intent(in) :: ivt         ! the land cover type
    real(r8), intent(in) :: assimn10    ! 10-day running mean net photosynthesis
    real(r8), intent(in) :: vegclass    ! ecophys constant
    real(r8), intent(in) :: raingreen   ! ecophys constant
    real(r8), intent(in) :: summergreen ! ecophys constant 
    real(r8), intent(in) :: t10min      ! annual minimum of 10-day running mean (K)
    real(r8), intent(in) :: agdd0       ! accumulated growing degree days above 0 deg C
    real(r8), intent(in) :: tref10      ! 10-day averaged temperature at 2m
    real(r8), intent(in) :: tmomin20    ! 20 year running mean of monthly minimum
    real(r8), intent(in) :: rstfac      ! 

                          !INTENT INOUT VARIABLES
    real(r8), intent(inout) :: lai_ind     ! lai in this patch
    real(r8), intent(inout) :: leafon      ! leafon days
    real(r8), intent(inout) :: leafof      ! leafoff days 
    real(r8), intent(inout) :: dphen       ! phenology [0 to 1]
    real(r8), intent(inout) :: htop        ! canopy top (m)

                          !INTENT OUT VARIABLES
    real(r8), intent(out) :: lai        !leaf area index
    real(r8), intent(out) :: sai        !Stem area index
    
                          !LOCAL VARIABLES:
    real(r8), parameter :: ddfacu = 1.0 / 15.0 ! 'drop day factor of trees'
    real(r8), parameter :: ddfacl = 1.0 /  5.0 ! 'drop day factor of grass'
    real(r8), parameter :: l_long = 1.0 /  2.0 ! leaf longevity of raingreen trees (not good enough)
    real(r8), parameter :: T0 = 273.16         ! water freezing point
    real(r8) :: tthreshold

    if (vegclass .le. 2.) then   

       !---------------------------------------------------------------------
       !* * * Evergreen Tree and Shrub Phenology * * *
       !---------------------------------------------------------------------

       if ((summergreen .lt. 0.) .and. (raingreen .lt. 0.)) then
          dphen = 1.0
        ! lai = max (min(8.0,lai_ind) * dphen, 0.05) 
          lai = max (lai_ind * dphen, 0.05) 
          sai = lai * 0.25 
       endif

       ! ---------------------------------------------------------------------
       ! * * * Summergreen Tree & Shrub phenology * * *
       ! ---------------------------------------------------------------------
       ! temperature threshold for budburst and senescence
       ! temperature threshold is assumed to be 273.16K
       ! or 5 degrees warmer than the coldest monthly temperature
       ! tmomin20 is initialized to tfrz-5.0

       if (summergreen .gt. 0.) then

          tthreshold = max (T0, tmomin20 + 5.0)

       ! determine if growing degree days are initiated
       ! slevis:*tref10 = 10-day running mean air temperature (K)
       ! determine leaf display

          if (tref10 < tthreshold) then
             dphen = max (0.0, dphen - ddfacu)
          else
             dphen = min (1.0, max (0.0, agdd0 - 100.0) / 50.0)
          end if

        ! lai = max (min(7.0,lai_ind) * dphen, 0.05) 
          lai = max (lai_ind * dphen, 0.05) 
          sai = lai*0.25

       endif

       ! ---------------------------------------------------------------------
       ! * * * Raingreen Tree & Shrub phenology * * *
       ! ---------------------------------------------------------------------

       if (raingreen .gt. 0.) then

          if (vegclass .le. 1.) then        !tree zhq. 07/22/2010

             if (assimn10 <  0.0) dphen = max (0.1, dphen - ddfacu)
             if (assimn10 >= 0.0) dphen = min (1.0, dphen + ddfacu)

           ! trying out enforced drought phenology
             if (dphen > 0.95) leafon = leafon + 1.0
             if (leafon >= 365.0*l_long) then
                dphen = 0.1
                leafof = leafof + 1.0
                if (leafof >= 365.0*l_long) then
                   dphen  = 1.0
                   leafof = 0.0
                   leafon = 1.0
                endif
             endif 

          else                              !shrub zhq. 07/22/2010

           ! different with tree raingreen, to avoid weild and unrealistic behavior
             if (rstfac > 0.25) then
                dphen = min(1.0, dphen + 0.25)
           ! else if (0. < rstfac < 0.15) then
             else if (0. < rstfac .and. rstfac < 0.15) then
                dphen = max(0.0, dphen - 0.1)
             endif
           
          endif

        ! lai = max (min(7.0,lai_ind)*dphen, 0.05)
          lai = max (lai_ind * dphen, 0.05)
          sai = lai*0.25

       endif

       ! ---------------------------------------------------------------------
       ! * * * Grass phenology * * *
       ! ---------------------------------------------------------------------
       ! temperature threshold for budburst and senescence
       ! temperature threshold is assumed to be 273.16K

    else                

       tthreshold = T0

     ! determine leaf display

       if (tref10 < tthreshold) then                 ! cold phenology for grasses
          dphen = max (0.0, dphen - ddfacl)   
       else if (assimn10 < 0.0) then                 ! drought phenology for grasses
          dphen = max (0.1, dphen - ddfacl)   
       else
          dphen = min (1.0, dphen + ddfacl)          ! from CLM-DGVM 
       end if

     ! lai = max (min(7.0,lai_ind) * dphen, 0.05)    ! attention!! 
       lai = max (lai_ind * dphen, 0.05)             ! attention!! 
       sai = lai*0.05
       htop = max(0.25, lai * 0.25)

    endif
 
 end subroutine DGVMphenology

