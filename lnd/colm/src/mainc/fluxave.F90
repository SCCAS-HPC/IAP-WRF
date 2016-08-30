
#include <define.h>

! ----------------------------------------------------------------------
!                             FLUX TABLE                               !
! ----------------------------------------------------------------------
! perfrom grid-average from subgrid 1d vector
! subgrid to grid average mapping: average a subgrid input vector [fldv] 
! of length numpatch to a output array [fldv] of length numgrid
!
! Created by Yongjiu Dai
!--------------!-------------------------------------------------------
! pft level fluxes
!--------------!-------------------------------------------------------
! 01: taux     ! wind stress: E-W [kg/m/s2]
! 02: tauy     ! wind stress: N-S [kg/m/s2]
! 03: fsena    ! sensible heat from canopy height to atmosphere [W/m2]
! 04: lfevpa   ! latent heat flux from canopy height to atmosphere [W/m2]
! 05: fevpa    ! evapotranspiration from canopy to atmosphere [mm/s]
! 06: fsenl    ! sensible heat from leaves [W/m2]
! 07: fevpl    ! evaporation+transpiration from leaves [mm/s]
! 08: etr      ! transpiration rate [mm/s]
! 09: fseng    ! sensible heat flux from ground [W/m2]
! 10: fevpg    ! evaporation heat flux from ground [mm/s]
! 11: fgrnd    ! ground heat flux [W/m2]
! 12: sabvsun  ! solar absorbed by sunlit canopy [W/m2]
! 13: sabvsha  ! solar absorbed by shaded [W/m2]
! 14: sabg     ! solar absorbed by ground [W/m2]
! 15: olrg     ! outgoing long-wave radiation from ground+canopy [W/m2]
! 16: rnet     ! net radiation [W/m2]
! 17: zerr     ! the error of energy balance [W/m2]
! 18: assim    ! canopy assimilation rate [mol m-2 s-1]
! 19: respc    ! respiration (plant+soil) [mol m-2 s-1]
! 20: fmicr    ! microbial respiration    [mol m-2 s-1]
! 21: tlsun    ! sunlit leaf temperature [K]
! 22: tlsha    ! shaded leaf temperature [K]
! 23: ldew     ! depth of water on foliage [mm]
! 24: sigf     ! fraction of veg cover, excluding snow-covered veg [-]
! 25: green    ! leaf greenness
! 26: lai      ! leaf area index
! 27: sai      ! stem area index
! 28: alb(1,1) ! averaged albedo [visible, direct]
! 29: alb(1,2) ! averaged albedo [visible, diffuse]
! 30: alb(2,1) ! averaged albedo [near-infrared, direct]
! 31: alb(2,2) ! averaged albedo [near-infrared,diffuse]
! 32: emis     ! averaged bulk surface emissivity
! 33: z0ma     ! effective roughness [m]
!--------------!-------------------------------------------------------
! 34: trad     ! radiative temperature of surface [K]
! 35: ustar    ! u* in similarity theory [m/s]
! 36: tstar    ! t* in similarity theory [kg/kg]
! 37: qstar    ! q* in similarity theory [kg/kg]
! 38: zol      ! dimensionless height (z/L) used in Monin-Obukhov theory
! 39: rib      ! bulk Richardson number in surface layer
! 40: fm       ! integral of profile function for momentum
! 41: fh       ! integral of profile function for heat
! 42: fq       ! integral of profile function for moisture
!--------------!-------------------------------------------------------
! 43: tref     ! 2 m height air temperature [kelvin]
! 44: qref     ! 2 m height air specific humidity [kg/kg]
! 45: u10m     ! 10m u-velocity [m/s]
! 46: v10m     ! 10m v-velocity [m/s]
! 47: f10m     ! integral of profile function for momentum at 10m [-]
!--------------!-------------------------------------------------------
! column level fluxes
!--------------!-------------------------------------------------------
! 48: xerr     ! the error of water banace [mm/s]
! 49: rsur     ! surface runoff [mm/s]
! 50: rnof     ! total runoff [mm/s]
!--------------!-------------------------------------------------------
! 51:60: tss   ! soil temperature [K]
! 61:70: wliq  ! liquid water in soil layers [kg/m2]
! 71:80: wice  ! ice lens in soil layers [kg/m2]
!--------------!-------------------------------------------------------
! 81: tg       ! ground surface temperature [K]
! 82: scv      ! snow cover, water equivalent [mm]
! 83: snowdp   ! snow depth [meter]
! 84: fsno     ! fraction of snow cover on ground
!----------------------------------------------------------------------
! 85: us       ! wind in eastward direction [m/s]
! 86: vs       ! wind in northward direction [m/s]
! 87: tm       ! temperature at reference height [kelvin]
! 88: qm       ! specific humidity at reference height [kg/kg]
! 89: prc      ! convective precipitation [mm/s]
! 90: prl      ! large scale precipitation [mm/s]
! 91: pbot     ! atmospheric pressure at the surface [pa]
! 92: frl      ! atmospheric infrared (longwave) radiation [W/m2]
! 93: solar    ! downward solar radiation at surface [W/m2]
!----------------------------------------------------------------------
! pft level fluxes
!----------------------------------------------------------------------
! 94:  qsubl        ! sublimation rate from snow pack [kg/m2/s]
!----------------------------------------------------------------------
! column level fluxes
!----------------------------------------------------------------------
! 95:104: mrlsl     ! mass of water of all phases in each soil layer [kg/m2]
! 105: mrsos        ! mass of water of all phases in the upper 0.1 meters of soil [kg/m2]
! 106: mrso         ! mass of water of all phases over all soil layers [kg/m2]
! 107: mrfso        ! mass of frozen water over all soil layers [kg/m2]
! 108: lwsnl        ! mass of liquid water of snow layers [kg/m2]
! 109: sm           ! surface snow melt [kg/m2/s]
! 110: tsn          ! snow internal temperature [K]
! 111: nsnow        ! number of snow events [-]
!----------------------------------------------------------------------
! 112: treeFrac     ! tree fraction [-]
! 113: shrubFrac    ! shrub fraction [-] 
! 114: grassFrac    ! grass fraction [-]
! 115: baresoilFrac ! bair soil fraction [-]
! 116: residualFrac ! residual fraction [-]
! 117: soilFrac     ! soil fraction [-]
! 118: urbanFrac    ! urban fraction [-]
! 119: wetlandFrac  ! wetland fraction [-]
! 120: iceFrac      ! ice fraction [-]
! 121: lakeFrac     ! lake & river fraction [-]
!----------------------------------------------------------------------

 subroutine fluxave

      use precision
      use phycon_module, only: vonkar, stefnc, cpair, rgas, grav
      use paramodel, only: grasscateg, oceancateg, nforc, nfldv, nfcon_pft, nfvar_pft, nfldv_col, nfldv_pft
      use spmd_decomp, only: pgmap, pcmap, cgmap, gxmap, gymap
      use colm_varMod, only: numpatch, numcolumn, numgrid, &
                             wt_patch=>wxy_patch, wt_column=>wxy_column, &
                             itypwat, fldv_col, fldv_pft, &
#ifdef DGVM
                             fvar_pft,   &
#else
                             fcon_pft,   &
#endif
                             forc, fldv, &
                             idx_fmicr, idx_scv

      implicit none

! local variables

      integer  :: i,j,L,f,c,p,g                       ! indices
      real(r8) :: a                                   ! 
      real(r8) :: sumwt(nfldv,numgrid)                ! sum of wt
      logical  :: iswater(numgrid)                    ! full grid is lake or ocean
      real(r8) rhoair,thm,th,thv,ur,displa,zldis,hu,ht,hq
      real(r8) z0m,z0h,z0q,us,vs,tm,qm,pbot,psrf
      real(r8) obu,temp1,temp2,temp12m,temp22m
      real(r8) um,thvstar,beta,zii,wc,wc2
      integer  idx_treeFrac, idx_shrubFrac, idx_grassFrac, idx_basesoilFrac, idx_residualFrac
      integer  idx_soilFrac, idx_urbanFrac, idx_wetlandFrac, idx_iceFrac, idx_lakeFrac
      integer  ivt

      sumwt(:,:) = 0.0
      fldv(:,:)  = 0.0

      iswater(:) = .true.

      do c = 1, numcolumn
         g = cgmap(c)
         if(itypwat(c).le.3) iswater(g) = .false.
      end do

![1-33] Grid averages by area-weight over grid patches

    ! Mapping the patch [numpatch] to grid [numgrid]
      do p = 1, numpatch
         j = pgmap(p)
         do L = 1, 33
            fldv(L,j) = fldv(L,j) + wt_patch(p)*fldv_pft(L,p)
            sumwt(L,j) = sumwt(L,j) + wt_patch(p)
         enddo
      enddo

      do j = 1, numgrid
         do L = 1, 33
!           if(sumwt(L,j).gt.1.0 .or. sumwt(L,j).lt.0.0)then
!              write(6,*) 'summation of fraction patches = ', sumwt(L,j),L,j
!              call abort
!           endif
            if(sumwt(L,j).gt.0.)then
               fldv(L,j) = fldv(L,j)/sumwt(L,j)
            else
               fldv(L,j) = -9999.
               write(6,*) 'impossible grid 1', sumwt(L,j), L
               call abort
            endif
         enddo
      enddo

![34-42] Retrieve through averaged fluxes
!     do k = 1, numpatch
!        j = pgmap(k)
!        do L = 70, 78
!           if(k.eq.maxfp(j))then  ! take values as that at the largest patches
!              fldv(L,j) = fldv(L,k)
!              sumwt(L,j) = 1.0
!           endif
!        enddo
!     enddo

      do p = 1, numpatch
         j = pgmap(p)
         do L = 34, 47
            sumwt(L,j) = sumwt(L,j) + wt_patch(p)
         enddo
      enddo

      c = 1               ! column index
      do j = 1, numgrid   ! grid index
         do while(j.ne.cgmap(c))
            c = c+1
         enddo

         if(sumwt(34,j).gt.0.)then         !For land only defined
            z0m = fldv(33,j)
            z0h = fldv(33,j)
            z0q = fldv(33,j)
            displa = 2./3.*z0m/0.07
       
            hu = max(forc(16,c),5.+displa)
            ht = max(forc(17,c),5.+displa)
            hq = max(forc(18,c),5.+displa)
            zldis = hu-displa
       
            us = forc(3,c)
            vs = forc(4,c)
            tm = forc(5,c)
            qm = forc(6,c)
            pbot = forc(9,c)
            psrf = forc(10,c)
     
            rhoair = (pbot-0.378*qm*pbot/(0.622+0.378*qm))/(rgas*tm)

            fldv(34,j) = (fldv(15,j)/stefnc)**0.25 
            fldv(35,j) = sqrt(max(1.e-6,sqrt(fldv(1,j)**2+fldv(2,j)**2))/rhoair) 
            fldv(36,j) = -fldv(3,j)/(rhoair*fldv(35,j))/cpair
            fldv(37,j) = -fldv(5,j)/(rhoair*fldv(35,j))
 
            thm = tm + 0.0098*ht
            th = tm*(100000./psrf)**(rgas/cpair)
            thv = th*(1.+0.61*qm)       
 
            fldv(38,j) = zldis*vonkar*grav&
                * (fldv(36,j)+0.61*th*fldv(37,j))&
                / (fldv(35,j)**2*thv)
 
            if(fldv(38,j) .ge. 0.)then   !stable
               fldv(38,j) = min(2.,max(fldv(38,j),1.e-6))
            else                           !unstable
               fldv(38,j) = max(-100.,min(fldv(38,j),-1.e-6))
            endif

            beta = 1.
            zii = 1000.
            thvstar=fldv(36,j)+0.61*th*fldv(37,j)
            ur = sqrt(us*us+vs*vs)
            if(fldv(38,j) .ge. 0.)then
               um = max(ur,0.1)
            else
               wc = (-grav*fldv(35,j)*thvstar*zii/thv)**(1./3.)
              wc2 = beta*beta*(wc*wc)
               um = max(0.1,sqrt(ur*ur+wc2))
            endif

            obu = zldis/fldv(38,j)
            call moninobuk(hu,ht,hq,displa,z0m,z0h,z0q,&
                 obu,um,fldv(35,j),temp1,temp2,temp12m,temp22m,&
                 fldv(47,j),fldv(40,j),fldv(41,j),fldv(42,j))
 
            fldv(39,j) = fldv(38,j)*vonkar**3*fldv(35,j)**2/(temp1*um**2)
            fldv(39,j) = min(5.,fldv(39,j)) 
         else
            fldv(34:42,j) = -9999.
            write(6,*) 'impossible grid 2', sumwt(34,j)
            call abort
         endif
    
      enddo

![42-46] Modified by zhq 06/12/2009/ for matching the routine meteorological obs. If there is grass pft in a grid, fluxes are output as that calculated on grass pft; else fluxes are output as average on a grid.
      do p = 1, numpatch
         j = pgmap(p)
         do L = 43, 47
            fldv(L,j) = fldv(L,j) + wt_patch(p)*fldv_pft(L,p)
          ! sumwt(L,j) = sumwt(L,j) + wt_patch(p)
         enddo
      enddo

      do p = 1, numpatch
         j = pgmap(p)
         do L = 43, 47
#if(defined DGVM)
            ivt = nint(fvar_pft(104,p))      
#else
            ivt = nint(fcon_pft(1,p))
#endif
            if((ivt.eq.grasscateg .or. ivt.eq.oceancateg).and. wt_patch(p).gt.1.0E-2)then
               fldv(L,j) = fldv_pft(L,p)
               sumwt(L,j) = 1.0
            endif 
         enddo
      enddo

    ! for incredible tref output, qian's test
    ! print*,'tref',fldv(43:47,:)

    ! Mapping the column [numcolumn] to grid [numgrid]
      do c = 1, numcolumn
         j = cgmap(c)
         do L = 48, 93
            f = L-47       !index of fldv_col
            if(L.le.50)then             ! fluxes
![48-50] Grid averages by area-weight over grid patches
               fldv(L,j) = fldv(L,j) + wt_column(c)*fldv_col(f,c)
               sumwt(L,j) = sumwt(L,j) + wt_column(c)
            else if(L.le.80)then        ! soil temperature and water
![51-60 61-70 71-80] Area-weight over grid patches but excluding lake and ocean patches
               if(itypwat(c).le.3)then  ! lake and ocean excluded
                  fldv(L,j) = fldv(L,j) + wt_column(c)*fldv_col(f,c)
                  sumwt(L,j) = sumwt(L,j) + wt_column(c)
               else if(iswater(j)) then
                  fldv(L,j) = fldv_col(f,c)
                  sumwt(L,j) = 1.0_r8
               endif
            else if(L.le.84)then                       ! clm state variables
![81-84] Grid averages by area-weight over grid patches
               fldv(L,j) = fldv(L,j) + wt_column(c)*fldv_col(f,c)
               sumwt(L,j) = sumwt(L,j) + wt_column(c)
            else                                       ! forcing variables
![85-93] 
               sumwt(L,j) = sumwt(L,j) + wt_column(c)
            endif
         enddo
      enddo

      do j = 1, numgrid
         do L = 43, 84
!           if(sumwt(L,j).gt.1.0 .or. sumwt(L,j).lt.0.0)then
!              write(6,*) 'summation of fraction patches = ', sumwt(L,j),L,j
!              call abort
!           endif
            if(sumwt(L,j).gt.0.)then
               fldv(L,j) = fldv(L,j)/sumwt(L,j)
            else
               fldv(L,j) = -9999.
               write(6,*) 'impossible grid 3', sumwt(L,j), L
               call abort
            endif
         enddo
      enddo

![84-92] Meteorological forcing
      c = 1
      do j = 1, numgrid
         do while(j.ne.cgmap(c))
            c = c+1
         enddo

         if(sumwt(85,j).gt.0.)then      !For land only defined
            fldv(85,j) = forc(3,c)
            fldv(86,j) = forc(4,c)
            fldv(87,j) = forc(5,c)
            fldv(88,j) = forc(6,c)
            fldv(89,j) = forc(7,c)
            fldv(90,j) = forc(8,c)
            fldv(91,j) = forc(9,c)
            fldv(92,j) = forc(15,c)
            fldv(93,j) = forc(11,c)+forc(12,c)+forc(13,c)+forc(14,c)
         else
            fldv(85:93,j) = -9999.
            write(6,*) 'impossible grid 4', sumwt(85,j)
            call abort
         endif
      enddo

![20] Special treatment for fmicr

      fldv(idx_fmicr,:) = 0.

      do p = 1, numpatch
         c = pcmap(p)
         g = pgmap(p)
         if(itypwat(c).eq.0) then
            fldv(idx_fmicr,g) = fldv(idx_fmicr,g) + fldv_pft(idx_fmicr,p)
         end if
      end do

      do c = 1, numcolumn
         g = cgmap(c)
         if(itypwat(c).eq.0) then
            fldv(idx_fmicr,g) = fldv(idx_fmicr,g)*wt_column(c)
         end if
      enddo

#ifdef CMIP

!PFT LEVEL

    ! index of fldv: 94
    ! index of fldv_pft: 48

      L = 94

      do p = 1, numpatch
         g = pgmap(p)
         fldv(L,g) = fldv(L,g) + wt_patch(p)*fldv_pft(48,p)
         sumwt(L,g) = sumwt(L,g) + wt_patch(p)
      enddo

      do g = 1, numgrid
         if(sumwt(L,g).gt.0._r8)then
            fldv(L,g) = fldv(L,g)/sumwt(L,g)
         else
            write(6,*) 'summation of fraction patches = ', sumwt(L,g), L, g
            call abort
         endif
      enddo

!COLUMN LEVEL

    ! index of fldv: 95-111
    ! index of fldv_col: 47-63

      do c = 1, numcolumn
         g = cgmap(c)
         do L = 95, 109
            f = L-48       !index of fldv_col
            fldv(L,g) = fldv(L,g) + wt_column(c)*fldv_col(f,c)
            sumwt(L,g) = sumwt(L,g) + wt_column(c)
         enddo
      enddo

      do g = 1, numgrid
         do L = 95, 109
            if(sumwt(L,g).gt.0._r8)then
               fldv(L,g) = fldv(L,g)/sumwt(L,g)
            else
               write(6,*) 'summation of fraction patches = ', sumwt(L,g), L, g
               call abort
            endif
         enddo
      enddo

    ! special treatment for tsn(idx=110) & nsnow(idx=111)

      do c = 1, numcolumn
         g = cgmap(c)
         do L = 110, 110
            f = L-48       !index of fldv_col
            if(fldv_col(f,c).gt.1.0E-6) then
               fldv(L,g) = fldv(L,g) + wt_column(c)*fldv_col(idx_scv-47,c)*fldv_col(f,c)
               sumwt(L,g) = sumwt(L,g) + wt_column(c)*fldv_col(idx_scv-47,c)
            end if
         enddo
      enddo

      do g = 1, numgrid
         do L = 110, 110
            if(sumwt(L,g).gt.0.)then
               fldv(L,g) = fldv(L,g)/sumwt(L,g)
               fldv(L+1,g) = 1._r8
            else
               fldv(L,g)   = 0._r8
               fldv(L+1,g) = 0._r8
            endif
         enddo
      enddo

    ! index of fldv: 112-121

      idx_treeFrac     = 112
      idx_shrubFrac    = 113
      idx_grassFrac    = 114
      idx_basesoilFrac = 115
      idx_residualFrac = 116
      idx_soilFrac     = 117
      idx_urbanFrac    = 118
      idx_wetlandFrac  = 119
      idx_iceFrac      = 120
      idx_lakeFrac     = 121

      fldv(idx_treeFrac,    :) = 0.
      fldv(idx_shrubFrac,   :) = 0.
      fldv(idx_grassFrac,   :) = 0.
      fldv(idx_basesoilFrac,:) = 0.
      fldv(idx_residualFrac,:) = 0.
      fldv(idx_soilFrac,    :) = 0.
      fldv(idx_urbanFrac,   :) = 0.
      fldv(idx_wetlandFrac, :) = 0.
      fldv(idx_iceFrac,     :) = 0.
      fldv(idx_lakeFrac,    :) = 0.

      do p = 1, numpatch
         c = pcmap(p)
         g = pgmap(p)

         ivt = nint(fvar_pft(104,p)) 

         if(ivt.ge.1 .and. ivt.le.8) then
            fldv(idx_treeFrac,g) = fldv(idx_treeFrac,g) + wt_patch(p)
         else if(ivt.ge.9 .and. ivt.le.11) then
            fldv(idx_shrubFrac,g) = fldv(idx_shrubFrac,g) + wt_patch(p)
         else if(ivt.ge.12 .and. ivt.le.14) then
            fldv(idx_grassFrac,g) = fldv(idx_grassFrac,g) + wt_patch(p)
         else if(ivt.eq.17) then
            fldv(idx_basesoilFrac,g) = fldv(idx_basesoilFrac,g) + wt_patch(p)
         else
            fldv(idx_residualFrac,g) = fldv(idx_residualFrac,g) + wt_patch(p)
         end if
      end do

      do c = 1, numcolumn
         g = cgmap(c)

         if(itypwat(c).eq.0) then
            fldv(idx_soilFrac,g) = wt_column(c)
         else if(itypwat(c).eq.1) then
            fldv(idx_urbanFrac,g) = wt_column(c)
         else if(itypwat(c).eq.2) then
            fldv(idx_wetlandFrac,g) = wt_column(c)
         else if(itypwat(c).eq.3) then
            fldv(idx_iceFrac,g) = wt_column(c)
         else if(itypwat(c).eq.4) then
            fldv(idx_lakeFrac,g) = wt_column(c)
         end if
      end do

      do g = 1, numgrid
         a = fldv(idx_treeFrac,g) + fldv(idx_shrubFrac,g) + fldv(idx_grassFrac,g) &
           + fldv(idx_basesoilFrac,g) + fldv(idx_residualFrac,g)

         if(abs(a-1._r8).gt.1.e-6_r8) then
            write(6,*), "fluxave: fractions error"         ,&
                                  fldv(idx_treeFrac,g)     ,&
                                  fldv(idx_shrubFrac,g)    ,&
                                  fldv(idx_grassFrac,g)    ,&
                                  fldv(idx_basesoilFrac,g) ,&
                                  fldv(idx_residualFrac,g)
            call abort
         end if
      end do

#endif

 end subroutine fluxave
