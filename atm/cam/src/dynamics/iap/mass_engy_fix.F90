subroutine mass_engy_fix (ztodt,  lat,    nlon,   u3m1,      &
                          u3,     v3m1,   v3,                &
                          t3m1,   t3,     q3m1,              &
                          q3,     psm1,   ps,                &
                          alpha,  etamid, qfcst,             &
                          qminus, beta,   hadv,              &
                          fu,     fv,     t2,                &
                          dqphy,  omga,   pdeldry )
!----------------------------------------------------------------------- 
! 
! Purpose: correct ps, q & T to conserve global mass and energy
! 
! Original: tfilt_massfixrun.F90 (cam3.1)
! 
! Modified by: Zhang He, 2008.4.18
! Update: Zhang He, 2008.09.11, added fu, fv, t2, and dqphy
!         Juanxiong He, 2010.08
!         2011.03.25 by Zhang He, do not use the energy fixer while running with ideal physics (H-S forcing)
! Reviewd: ZhangHe, 2011-11-18
! Modified: ZhangHe, 2012-01-15, changed pdeldry(:,:) to pdeldry(plon,plev)
!         : Jiang Jinrong, October 2012, for 2D parellel
! Reviewd: ZhangHe, 2012-10-31
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use cam_history, only: outfld
   use pmgrid,  only: plon, plev, plevp, plat, beglev, endlev, npr_z
   use pspect
   use commap
   use constituents, only: pcnst, qmin, cnst_cam_outfld, &
                           tottnam, tendnam, cnst_get_type_byind, fixcnam, &
                           hadvnam, vadvnam
   use time_manager, only: get_nstep
   use physconst,    only: cpair, gravit
   use eul_control_mod, only: fixmas
   use cam_control_mod, only: ideal_phys
#if ( defined SPMD )
   use spmd_dyn          , only:  comm_z
   use parutilitiesmodule, only: parcollective3d,sumop
#endif

   implicit none

!
! Input arguments
!
   real(r8), intent(in) :: ztodt                  ! timestep
   integer,  intent(in) :: lat
   integer,  intent(in) :: nlon

   real(r8), intent(in) :: qfcst(plon,beglev:endlev,pcnst) ! slt moisture forecast
   real(r8), intent(in) :: qminus(plon,beglev:endlev,pcnst)
   real(r8), intent(in) :: beta                   ! energy fixer coefficient
   real(r8), intent(in) :: hadv(plon,beglev:endlev,pcnst)  ! horizonal q advection tendency
   real(r8), intent(in) :: alpha(pcnst)
   real(r8), intent(in) :: etamid(plev)           ! vertical coords at midpoints 
   real(r8), intent(in) :: u3(plon,beglev:endlev)
   real(r8), intent(in) :: v3(plon,beglev:endlev)
   real(r8), intent(in) :: psm1(plon)
   real(r8), intent(in) :: u3m1(plon,beglev:endlev)
   real(r8), intent(in) :: v3m1(plon,beglev:endlev)
   real(r8), intent(in) :: t3m1(plon,beglev:endlev)
   real(r8), intent(in) :: q3m1(plon,beglev:endlev,pcnst)
   real(r8), intent(in) :: fu(plon,beglev:endlev)      ! du/dt due to physics
   real(r8), intent(in) :: fv(plon,beglev:endlev)      ! dv/dt due to physics
   real(r8), intent(in) :: t2(plon,beglev:endlev)      ! dT/dt due to physics
   real(r8), intent(in) :: dqphy(plon,beglev:endlev)   ! dq/dt due to physics
!
! Input/Output arguments
!
   real(r8), intent(inout) :: omga(plon,beglev:endlev)
   real(r8), intent(inout) :: t3(plon,beglev:endlev)
   real(r8), intent(inout) :: pdeldry(plon,beglev:endlev)      ! dry pressure difference at time n3
   real(r8), intent(inout) :: q3(plon,beglev:endlev,pcnst)
   real(r8), intent(inout) :: ps(plon)
!
!---------------------------Local workspace-----------------------------
!
   integer :: ifcnt                  ! Counter
   integer :: nstep                  ! current timestep number

   real(r8) :: q3tmp(nlon,plev,pcnst)
   real(r8) :: tfix    (plon)        ! T correction
   real(r8) :: engycorr(plon,beglev:endlev)   ! energy equivalent to T correction
   real(r8) :: pdel(plon,beglev:endlev)       ! pdel(k) = pint(k+1) - pint(k)
   real(r8) :: pint(plon,beglev:endlev+1)      ! pressure at model interfaces (n  )
   real(r8) :: pmid(plon,beglev:endlev)       ! pressure at model levels (time n)
   real(r8) :: utend(plon,beglev:endlev)      ! du/dt
   real(r8) :: vtend(plon,beglev:endlev)      ! dv/dt
   real(r8) :: ttend(plon,beglev:endlev)      ! dT/dt
   real(r8) :: qtend(plon,beglev:endlev,pcnst)! dq/dt
   real(r8) :: pstend(plon)          ! d(ps)/dt
   real(r8) :: vadv(plon,beglev:endlev,pcnst) ! vertical q advection tendency
   real(r8) :: pintm1(plon,plevp)    ! pressure at model interfaces (n-1)
   real(r8) :: pmidm1(plon,plev)     ! pressure at model levels (time n-1)
   real(r8) :: pdelm1(plon,plev)     ! pdelm1(k) = pintm1(k+1)-pintm1(k)
   real(r8) :: corm
   real(r8) :: wm
   real(r8) :: absf
   real(r8) :: worst
   logical  :: lfixlim               ! flag to turn on fixer limiter

   real(r8) :: ta(plon,beglev:endlev,pcnst)    ! total advection of constituents
   real(r8) :: dqfx3(plon,beglev:endlev,pcnst) ! q tendency due to mass adjustment
   integer  :: i, k, m,j,ixcldliq,ixcldice  ! loop indices

!-----------------------------------------------------------------------
   nstep   = get_nstep()
   lfixlim = .true.

!
! Set average dry mass to specified constant preserving horizontal
! gradients of ln(ps). Proportionality factor was calculated in STEPON
! for nstep=0 or SCAN2 otherwise from integrals calculated in INIDAT
! and SCAN2 respectively.
! Set p*.
!
   do i=1,nlon
      ps(i) = ps(i)*fixmas
   end do
!
! Set current time pressure arrays for model levels etc.
!
   call plevs00(nlon    ,plon   ,plev    ,ps      ,pint    ,pmid    ,pdel)
!
! Add temperature correction for energy conservation
!
   do k=beglev,endlev
      do i=1,nlon
         engycorr(i,k) = (cpair/gravit)*beta*pdel(i,k)/ztodt
         if (.not. ideal_phys) then     !zhh 2011.03.25  
            t3(i,k) = t3(i,k) + beta
         end if 
      end do
   end do
   do i=1,nlon
      tfix(i) = beta/ztodt
   end do
!
! Output Energy correction term
!
   call outfld ('ENGYCORR',engycorr ,plon   ,lat     )
   call outfld ('TFIX    ',tfix     ,plon   ,lat     )
!
! Compute q tendency due to mass adjustment
! If LFIXLIM = .T., then:
! Check to see if fixer is exceeding a desired fractional limit of the
! constituent mixing ratio ("corm").  If so, then limit the fixer to
! that specified limit.
!
   do m=1,pcnst
      if (cnst_get_type_byind(m).eq.'dry' ) then
         corm    = 1.e36
      else
         corm    = 0.1
      end if

      do k=beglev,endlev
         do i=1,nlon
!! see cam3 description P44, (3.238)
            dqfx3(i,k,m) = alpha(m)*etamid(k)*abs(qfcst(i,k,m) - qminus(i,k,m))
         end do
         if (lfixlim) then
            ifcnt = 0
            worst = 0.
            wm    = 0.
            do i = 1,nlon
               absf = abs(dqfx3(i,k,m))
               if (absf.gt.corm) then
                  ifcnt = ifcnt + 1
                  worst = max(absf,worst)
                  wm = wm + absf
                  dqfx3(i,k,m) = sign(corm,dqfx3(i,k,m))
!! if dqfx3 >= 0, dqfx3 = abs(corm); dqfx3 < 0, dqfx3 = -abs(corm)
               endif
            end do
            if (ifcnt.gt.0) then
               wm = wm/float(ifcnt)
            endif
         endif
         do i=1,nlon
            dqfx3(i,k,m) = qfcst(i,k,m)*dqfx3(i,k,m)/ztodt
            q3   (i,k,m) = qfcst(i,k,m) + ztodt*dqfx3(i,k,m)
            ta   (i,k,m) = (q3    (i,k,m) - qminus(i,k,m))/ztodt
            vadv (i,k,m) = (qfcst(i,k,m) - qminus(i,k,m))/ztodt - hadv(i,k,m)
         end do
      end do
   end do
   do k=beglev,endlev
      do i=1,nlon
         pdeldry(i,k) = pdel(i,k)*(1.-q3(i,k,1))
      end do ! i
   end do ! k

!
! Check for and correct invalid constituents
!
#if (defined SPMD)
   if (npr_z.gt.1) then
      q3tmp(:,:,:)=0.
   end if
#endif
   q3tmp(1:plon,beglev:endlev,1:pcnst)=q3(1:plon,beglev:endlev,1:pcnst)

#if (defined SPMD)
   if (npr_z.gt.1) then
      call parcollective3d( comm_z, sumop,plon ,plev,pcnst ,q3tmp)
   endif
#endif
!
   call qneg3 ('MASS_ENGY_FIX',lat   ,nlon    ,plon   ,plev    , &
               1, pcnst,qmin ,q3tmp(1,1,1))
!
   q3(1:plon,beglev:endlev,1:pcnst)=q3tmp(1:plon,beglev:endlev,1:pcnst)
!
! Send slt tendencies to the history tape
!
   do m=1,pcnst
      if ( cnst_cam_outfld(m) ) then
         call outfld(tottnam(m),ta(1,1,m),plon   ,lat     )
      end if
   end do
!
! correct vertical motion field according to ps
!
   do k=beglev,endlev
      do i=1,nlon
         omga(i,k) = omga(i,k)*fixmas
      end do
   end do
!
! Compute time tendencies:comment out since currently not on h-t
!
   do k=beglev,endlev
      do i=1,nlon
         ttend(i,k) = (t3(i,k)-t3m1(i,k))/ztodt
         utend(i,k) = (u3(i,k)-u3m1(i,k))/ztodt
         vtend(i,k) = (v3(i,k)-v3m1(i,k))/ztodt
      end do
   end do

   do m=1,pcnst
      do k=beglev,endlev
         do i=1,nlon
            qtend(i,k,m) = (q3(i,k,m) - q3m1(i,k,m))/ztodt
         end do
      end do
   end do

   do i=1,nlon
      pstend(i) = (ps(i) - psm1(i))/ztodt
   end do

   do m=1,pcnst
      if ( cnst_cam_outfld(m) ) then
         call outfld (tendnam(m),qtend(:,:,m),plon,lat)     !jjr
         call outfld (fixcnam(m),dqfx3(:,:,m),plon,lat)
         call outfld (hadvnam(m),hadv (:,:,m),plon,lat)
         call outfld (vadvnam(m),vadv (:,:,m),plon,lat)
      end if
   end do

   call outfld ('UTEND   ',utend,plon,lat)
   call outfld ('VTEND   ',vtend,plon,lat)
   call outfld ('TTEND   ',ttend,plon,lat)
   call outfld ('LPSTEN  ',pstend,plon,lat)
! ========================= zhh 2008.9.11 ===============================
!!   call outfld ('DUPHY   ',fu,plon,lat)
!!   call outfld ('DVPHY   ',fv,plon,lat)
!!   call outfld ('DTPHY   ',t2,plon,lat)
!!   call outfld ('DQPHY   ',dqphy,plon,lat)
! ========================= zhh 2008.9.11 ===============================

   return
1000 format(' mass_engy_fix: WARNING: fixer for tracer ',i3,' exceeded ', &
      f8.5,' for ',i5,' points at k,lat = ',2i4, &
      ' Avg/Worst = ',1p2e10.2)

end subroutine  mass_engy_fix

