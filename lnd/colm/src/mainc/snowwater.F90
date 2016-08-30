
 subroutine snowwater (lb, dtime, ssi, wimp, &
                       pg_rain, qseva, qsdew, qsubl, qfros, &
                       dz, wice, wliq, qout_snowb)

!-----------------------------------------------------------------------
! Original author : Yongjiu Dai, September 15, 1999
!
! Water flow wihtin snow is computed by an explicit and non-physical based scheme,
! which permits a part of liquid water over the holding capacity (a tentative value
! is used, i.e., equal to 0.033*porosity) to percolate into the underlying layer,
! except the case of that the porosity of one of the two neighboring layers is
! less than 0.05, the zero flow is assumed. The water flow out of the bottom
! snow pack will participate as the input of the soil water and runoff.
!-----------------------------------------------------------------------

  use precision
  use phycon_module, only : denice, denh2o  ! physical constant
  implicit none

!----------------------- dummy argument --------------------------------
  integer, INTENT(in) :: &
        lb          ! lower bound of array

  real(r8), INTENT(in) :: &
        dtime,     &! time step (s)
        ssi,       &! irreducible water saturation of snow
        wimp,      &! water impremeable if porosity less than wimp
        dz(lb:0),  &! layer thickness (m)

        pg_rain,   &! rainfall after removal of interception (mm h2o/s)
        qseva,     &! ground surface evaporation rate (mm h2o/s)
        qsdew,     &! ground surface dew formation (mm h2o /s) [+]
        qsubl,     &! sublimation rate from snow pack (mm h2o /s) [+]
        qfros       ! surface dew added to snow pack (mm h2o /s) [+]

  real(r8), INTENT(inout) :: &
        wice(lb:0),&! ice lens (kg/m2)
        wliq(lb:0)  ! liquid water (kg/m2)
  
  real(r8), INTENT(out) :: &
        qout_snowb  ! rate of water out of snow bottom (mm/s)

!----------------------- local variables --------------------------------
  integer j         ! k do loop/array indices

  real(r8) :: &
       qin,        &! water flow into the elmement (mm/s)
       qout,       &! water flow out of the elmement (mm/s)
       zwice,      &! the sum of ice mass of snow cover (kg/m2)
       wgdif,      &! ice mass after minus sublimation
    vol_liq(lb:0), &! partitial volume of liquid water in layer
    vol_ice(lb:0), &! partitial volume of ice lens in layer
 eff_porosity(lb:0) ! effective porosity = porosity - vol_ice

!=======================================================================
! renew the mass of ice lens (wice) and liquid (wliq) in the surface snow layer, 
! resulted by sublimation (frost) / evaporation (condense)

      wgdif = wice(lb) + (qfros - qsubl)*dtime
      wice(lb) = wgdif
      if(wgdif < 0.)then
         wice(lb) = 0.
         wliq(lb) = wliq(lb) + wgdif
      endif
      wliq(lb) = wliq(lb) + (pg_rain + qsdew - qseva)*dtime
      wliq(lb) = max(0., wliq(lb))

! Porosity and partitial volume
      do j = lb, 0
         vol_ice(j) = min(1., wice(j)/(dz(j)*denice))
         eff_porosity(j) = 1. - vol_ice(j)
         vol_liq(j) = min(eff_porosity(j), wliq(j)/(dz(j)*denh2o))
      enddo

! Capillary forces within snow are usually two or more orders of magnitude
! less than those of gravity. Only gravity term are considered. 
! the genernal expression for water flow is "K * ss**3", however, 
! no effective paramterization for "K". Thus, a very simple consideration 
! (not physical based) is introduced: 
! when the liquid water of layer exceeds the layer's holding 
! capacity, the excess meltwater adds to the underlying neighbor layer. 

      qin = 0.
      do j= lb, 0
         wliq(j) = wliq(j) + qin

         if(j <= -1)then
! no runoff over snow surface, just ponding on surface
           if(eff_porosity(j)<wimp .OR. eff_porosity(j+1)<wimp)then
             qout = 0.
           else
             qout = max(0.,(vol_liq(j)-ssi*eff_porosity(j))*dz(j))
             qout = min(qout,(1.-vol_ice(j+1)-vol_liq(j+1))*dz(j+1))
           endif
         else
           qout = max(0.,(vol_liq(j)-ssi*eff_porosity(j))*dz(j))
         endif

         qout    = qout*1000.
         wliq(j) = wliq(j) - qout
         qin     = qout

      enddo

      qout_snowb = qout/dtime

 end subroutine snowwater
