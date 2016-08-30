subroutine dynpkg( adv_state, t2, fu, fv, qminus, flx_net, t3, u3, v3, q3, ps,     &
                    omga, phis, etamid, cwava, detam, dtime )
!-------------------------------------------------------------------------------------------------------
! Purpose: Driving routines for dynamics and horizontal diffusion.
! Author : ZhangHe
! Completed: 2007.4.28
! Update:  ZhangHe, 2007.12.17, 'q3(plon,plev,plat,npt)' ==> q3(plon,plev,ppcnst,plat,npt)
!          ZhangHe, 2008.4.20, added calling sub. mass_engy
!          ZhangHe, 2008.4.27, added calling sub. trans_af_mass
!          ZhangHe, 2010.8.10, added NSLT
!          juanxiong He, 2010.08
!          ZhangHe, 2011.05.04, do not call HDIFUS and mass_engy if adiabatic run
!          2011.7.10 Juanxiong He, adapted to CESM, 'q3(plon,plev,ppcnst,plat,npt)' ==> q3(plon,plev,pcnst,plat,npt)
!          ZhangHe, 2011-11-15
!          Jiang Jinrong, OCT 2012, 2D parellel
!          ZhangHe, 2012-11-08, added calling sub. init_trans_pd and sub. init_trans_IAP
!          ZhangHe, 2013-02-05, reviewed, removed redundant variables in use only statement
! Modified: ZhangHe, 2013-03-29, use dynamic arrays dqphy
!           zhangHe, 2013-04-17, removed calling sub. init_trans_pd and sub. init_trans_IAP
!--------------------------------------------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,     only: plon, plev, plat, beglat, endlat,beglev,endlev, beglatdyn, endlatdyn
   use flexib,     only: DTadv, NENG, NFST, DTIMEQ
   use IAP_prog,   only: U, V, T, PLY, PT, UT, VT, TT, GHS, P, Pstar1, PS2, psa, Q, Qliq
   use tendency,   only: NADDSS, DT, ST  
   use hdif,       only: HDIFUS
   use prognostics, only: n3, n3m1,ptimelevels  !juanxiong he, 2010.08  !zhh 2011-11-15
   use constituents, only: pcnst                     !juanxionghe 2010.08
   use scanslt,      only: advection_state, qfcst           !zhh 2008.4.21
   use Trans_coef,   only: trans_af_mass, PTU, PTV, PTT       !zhh 2008.4.27
   use stdatm,   only : TB, PSB
   use perf_mod, only : t_startf, t_stopf ! juanxiong he
   use cam_control_mod, only: adiabatic        !zhh, 2011-12-15
   use spmd_utils, only: masterproc     !zhh, 2013-02-18, for debug

   implicit none
!------------------------------Arguments--------------------------------
   real(r8), intent(in) :: fu (plon,beglev:endlev,beglat:endlat)   ! tendency of u-wind
   real(r8), intent(in) :: fv (plon,beglev:endlev,beglat:endlat)   ! tendency of v-wind
   real(r8), intent(in) :: t2 (plon,beglev:endlev,beglat:endlat)   ! tendency of temperature
   real(r8), intent(inout) :: qminus(plon,beglev:endlev,pcnst,beglat:endlat) ! constituents
   real(r8), intent(in) :: flx_net(plon,beglat:endlat)  ! net flux from physics
   real(r8), intent(in) :: dtime                        ! timestep of stepon
!
   real(r8), intent(inout) :: t3(plon,beglev:endlev,beglat:endlat,ptimelevels)  ! temperature
   real(r8), intent(inout) :: u3(plon,beglev:endlev,beglat:endlat,ptimelevels)  ! u-wind component
   real(r8), intent(inout) :: v3(plon,beglev:endlev,beglat:endlat,ptimelevels)  ! v-wind component
   real(r8), intent(inout) :: q3(plon,beglev:endlev,pcnst,beglat:endlat,ptimelevels)  ! specific humidity, juanxiong he, 201008
   real(r8), intent(inout) :: omga(plon,beglev:endlev,beglat:endlat)    ! p-surface vertical velocity
   real(r8), intent(inout) :: ps(plon,beglat:endlat,ptimelevels)       ! surface pressure
   real(r8), intent(inout) :: phis(plon,beglat:endlat)         ! surface geopotential
   type(advection_state), intent(inout) :: adv_state        ! Advection state data
   real(r8), intent(in) :: etamid(plev)     ! vertical coords at midpoints
   real(r8), intent(inout) :: cwava(plat)   ! weight applied to global integrals
   real(r8), intent(inout) :: detam(plev)   ! intervals between vert full levs.

!---------------------------Local workspace-----------------------------
   integer :: NSEQ     ! dynamical computation times of every dtime
   integer :: NSLT     ! slt computation times of every DTadv
   real(r8), allocatable :: dqphy(:,:,:)     ! q tendency due to physics, zhh 2008.9.10
   integer :: j,k,i,jd
!------------------------------------------------------------------------------------------------
! allocate arrays
   allocate ( dqphy(plon,beglev:endlev,beglat:endlat) )

! ======================= added by zhh 2008.9.10 =========================
! compute q tendency due to physics
   do j=beglat, endlat
      jd = plat + 1 -j
      do k=beglev,endlev
         do i=1,plon
            dqphy(i,k,j) = (qminus(i,k,1,j) - q3(i,k,1,j,n3m1))/dtime
         end do
      end do
   end do
! ======================= added by zhh 2008.9.10 =========================
!
   NSEQ = dtime / DTadv + 0.01D0
   NSLT = DTadv / DTIMEQ + 0.01D0
! ============================================================================
   call trans_grid_pd( fu, fv, t2, qminus, beglat, endlat, pcnst )
!=============================================================================
!
   IF ( NADDSS == 0 ) THEN
! ============ do first half-step horizontal diffusion ===================
      if(.not. adiabatic) CALL HDIFUS( dtime, PLY, U, V, T )
! ========================================================================
   ENDIF
!
   call t_startf('DYFRAM')
! ============== perform the dynamic integration cycle ==================
   call DYFRAM( NENG,NFST,NSEQ,NSLT, adv_state, detam, etamid, cwava,   &
                 t3, u3, v3, q3, ps, omga, phis, n3, ptimelevels )   !juanxiong he, 201008
! =======================================================================
   call t_stopf('DYFRAM')
!
! =======================================================================
!  do other half-step horizontal diffusion ( NADDSS == 0 ) .or.
!  do first half-step horizontal diffusion ( NADDSS /= 0 )
   if(.not. adiabatic) CALL HDIFUS( dtime, PLY, U, V, T)
! =======================================================================
!
   IF ( NADDSS /= 0 ) THEN
! ============ do other half-step horizontal diffusion ==================
      if(.not. adiabatic) CALL HDIFUS( dtime, PLY, U, V, T )
! =======================================================================
   ENDIF
!=============================== 2007.12.20 ===============================
   call copy_grid_dp( 1, t3, u3, v3, q3, omga, ps, phis, n3, ptimelevels ) !juanxiong he, 201008
!========================================================================
! mass and energy correction
   if(.not. adiabatic) call mass_engy (dtime, cwava, etamid, flx_net, fu, fv, t2, dqphy)  !zhh 2008.9.11
!
! Convert ps, t3, q3 to P, T, Q
   if(.not. adiabatic) call trans_af_mass( ps, t3, q3, pcnst, n3, ptimelevels )  
!
! deallocate arrays
   deallocate(dqphy)

! ============================ for test =================================
   do j = beglat, endlat
      do k = beglev, endlev
         if ( j==64 .and. k==26) then
            write(6,*) '--------------- At the end of sub. dynpkg --------------------'
            write(6,*) 'UT(4,26,65) =', UT(4,26,65), ' T(4,26,65) =', T(4,26,65)
            write(6,*) 'V(4,26,65) =', V(4,26,65), ' ps(4,64,n3) =', ps(4,64,n3)
            write(6,*) 'v3(4,26,64,n3) =', v3(4,26,64,n3)
            write(6,*) 'q3(4,26,1,64,n3) =', q3(4,26,1,64,n3)
            write(6,*) 'q3(4,26,4,64,n3) =', q3(4,26,4,64,n3)
            write(6,*) '--------------------------------------------------------------'
         end if
!!         if ( j==127 .and. k==1) then
!!            write(6,*) 'U(126,1,2) =', U(126,1,2)
!!         end if
      end do
   end do
!
   do j = beglatdyn, endlatdyn
      do k = beglev, endlev
         do i = 1, plon
            if ( abs(U(i,k,j))>300.0 .or. abs(V(i,k,j))>300.0 .or. T(i,k,j)<100.0 ) then
               print*, 'U(',i,k,j,') =', U(i,k,j)
               print*, 'V(',i,k,j,') =', V(i,k,j)
               print*, 'T(',i,k,j,') =', T(i,k,j)
!!               stop
           end if
         end do
!!         if ( j==115 .and. k==1) then
!!            write(6,*) 'U(87,1,115) =', U(87,1,115)
!!         end if
      end do
   end do
! =============================== zhh ===================================
!
   return

end subroutine dynpkg

