!=================================================================================
subroutine mass_engy (ztodt, cwava, etamid, flx_net, fu, fv, t2, dqphy)  !zhh 2008.9.11
!---------------------------------------------------------------------------------
! Purpose: Driving routines for mass and energy correction
! Original version: scan2.F90 (cam3.1)
! Modified : ZhangHe
! Completed : 2008.4.18
! Update:  Zhanghe, 2008.6.8, added sub. realloc5 for parallel version
! update: Juanxiong He, 2010.08	
! Reviewed: ZhangHe, 2011-11-18
!           ZhangHe, 2012-01-15
! Modified: Jiang Jinrong, OCT 2012, for 2D parallel 
!           ZhangHe, 2013-02-05, removed sub. realloc5 
!---------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plat, plev, plevp, beglat, endlat,  &
                           beglev, endlev, beglatdyn, npr_z
   use prognostics,  only: n3, n3m1, ps, u3, v3, q3, t3, qminus,   &
                           phis, omga, hadv, pdeld
   use rgrid,        only: nlon
   use scanslt,      only: hw1lat, engy1lat, qfcst
   use constituents, only: pcnst
   use physconst,    only: cpair, rga
   use massfix,      only: hw1,hw2,hw3,alpha
   use IAP_grid,     only: ie1,ex
#ifdef SPMD
   use mpishorthand, only: mpicom, mpir8
   use spmd_utils,   only: iam, masterproc    ! debug
   use parutilitiesmodule, only : sumop, parcollective  ! jjr
#endif
   use eul_control_mod, only: tmass, tmass0, tmassf, qmassf ,fixmas
   use cam_control_mod, only: adiabatic, ideal_phys 
   use perf_mod, only : t_startf, t_stopf, t_barrierf !juanxiong he, ZhangHe   

   implicit none
!
!--------------------------Input arguments------------------------------
!
   real(r8), intent(in) :: ztodt                ! timestep
   real(r8), intent(in) :: cwava(plat)          ! weight applied to global integrals
   real(r8), intent(in) :: etamid(plev)         ! vertical coords at midpoints 
   real(r8), intent(in) :: flx_net(plon,beglat:endlat) ! net flux from physics
! ========================= zhh 2008.9.11 ===============================
   real(r8), intent(in) :: fu(plon,beglev:endlev,beglat:endlat)  ! u-wind tendency due to physics 
   real(r8), intent(in) :: fv(plon,beglev:endlev,beglat:endlat)  ! v-wind tendency due to physics 
   real(r8), intent(in) :: t2(plon,beglev:endlev,beglat:endlat)  ! temperature tendency due to physics 
   real(r8), intent(in) :: dqphy(plon,beglev:endlev,beglat:endlat)  ! q tendency due to physics 
! ========================= zhh 2008.9.11 ===============================
!
!---------------------------Local workspace-----------------------------
!
   real(r8) :: engy1         ! component of global energy integral (for time step n)
   real(r8) :: engy2         ! component of global energy integral (for time step n+1)
   real(r8) :: difft         ! component of global delta-temp integral ( (n+1) - n )
   real(r8) :: beta          ! 
   real(r8) :: hwx(pcnst,4)  ! component of constituent global mass integral
   real(r8) :: engy2lat(plat)     ! lat contribution to total energy integral
   real(r8) :: difftlat(plat)     ! lat contribution to delta-temperature integral
   real(r8) :: hw2l(pcnst,plat)   !  | latitudinal contributions to the
   real(r8) :: hw3l(pcnst,plat)   ! <  components of global mass integrals
   real(r8) :: hwxl(pcnst,4,plat) !  |
   real(r8) :: temp(5*pcnst+2)
!
   real(r8) :: residual          ! residual energy integral
   real(r8) :: tmp               ! accumulator
   integer  :: lat, m, n, jdry, j, tt, i, k         ! indices
!
!-----------------------------------------------------------------------
!
! Set coefficient used for mass and energy fixer
!
   engy1  = 0.
   engy2  = 0.
   difft  = 0.
   do m=1,pcnst
      hw1(m)=0.
      hw2(m) = 0.
      hw3(m) = 0.
      do n=1,4
         hwx(m,n) = 0.
      end do
   end do
   tmassf = 0.
   do m=1,5*pcnst+2
      temp(m)=0.
   end do
!
   do lat=beglat,endlat
      call set_mass_engy (ztodt, lat, nlon(lat), cwava(lat),  &
                   qfcst(:,:,:,lat), qminus(:,:,:,lat), etamid, &
                   ps(:,lat,n3), u3(:,:,lat,n3), v3(:,:,lat,n3), &
                   t3(:,:,lat,n3), flx_net(1,lat), phis(1,lat), & 
                   ps(:,lat,n3m1), u3(:,:,lat,n3m1), v3(:,:,lat,n3m1), &
                   t3(:,:,lat,n3m1), hw2l(1,lat), hw3l(1,lat), &
                   hwxl(:,:,lat), engy1lat(lat), engy2lat(lat), difftlat(lat) )
!jjr
      temp(1)=temp(1)+engy1lat(lat)
      temp(2)=temp(2)+ engy2lat(lat)
      temp(3)=temp(3)+ difftlat(lat)
      temp(4)=temp(4)+hw1lat(1,lat)
      temp(5)=temp(5)+hw2l(1,lat)
      temp(6)=temp(6)+hw3l(1,lat)
      temp(7)=temp(7)+ tmass(lat)
!
      tt=7
      do m = 2,pcnst
         tt=tt+1
         temp(tt)=temp(tt)+ hw1lat(m,lat)

         do n = 1,4
            tt=tt+1
            temp(tt)=temp(tt)+hwxl(m,n,lat)
         end do
      end do
!
   end do
!
#ifdef SPMD
   call parcollective(mpicom,sumop,5*pcnst+2,temp)
#endif
!
   engy1  = temp(1) 
   engy2  = temp(2)
   difft  = temp(3)
   hw1(1) = temp(4)
   hw2(1) = temp(5)
   hw3(1) = temp(6)
   tmassf = temp(7)/npr_z
   tt=7
   do m = 2,pcnst
      tt=tt+1
      hw1(m) =temp(tt)

      do n = 1,4
         tt=tt+1
         hwx(m,n) = temp(tt)
      end do
   end do
!
! Compute atmospheric mass fixer coefficient
!
   tmassf = tmassf * 0.5
   qmassf = hw1(1)
!
   if (adiabatic .or. ideal_phys) then
      fixmas = tmass0/tmassf
   else
      fixmas = (tmass0 + qmassf)/tmassf
   end if
!
! Compute alpha for water ONLY
!
   hw2(1)    = fixmas*hw2(1)
   hw3(1)    = fixmas*hw3(1)
   if(hw3(1) .ne. 0.) then
      alpha(1)  = ( hw1(1) - hw2(1) )/hw3(1)
   else
      alpha(1)  = 1.
   endif
!
! Compute alpha for non-water constituents
!
   do m = 2,pcnst
      hw2(m) = hwx(m,1) - alpha(1)*hwx(m,2)
      hw3(m) = hwx(m,3) - alpha(1)*hwx(m,4)
      hw2(m) = fixmas*hw2(m)
      hw3(m) = fixmas*hw3(m)
      if(hw3(m) .ne. 0.) then
         alpha(m)  = ( hw1(m) - hw2(m) )/hw3(m)
      else
         alpha(m)  = 1.
      end if
   end do
!
! Compute beta for energy
!
   engy2    = fixmas*engy2
   difft    = fixmas*difft
   residual = (engy2 - engy1)/ztodt
   if(difft .ne. 0.) then
     beta = -residual*ztodt/(cpair*difft)
   else
     beta = 0.
   endif
!
   call t_startf ('mass_engy_fix')
!
!$OMP PARALLEL DO PRIVATE (LAT)
!
   do lat=beglat,endlat
       jdry = lat  !juanxiong he, 201008
      call mass_engy_fix (ztodt,     lat,   nlon(lat),        u3(:,:,lat,n3m1),    &
                          u3(:,:,lat,n3),   v3(:,:,lat,n3m1), v3(:,:,lat,n3),      &
                          t3(:,:,lat,n3m1), t3(:,:,lat,n3),   q3(:,:,:,lat,n3m1),  &
                          q3(:,:,:,lat,n3), ps(1,lat,n3m1),   ps(1,lat,n3),        &
                          alpha,            etamid,           qfcst(:,:,:,lat),    &
                          qminus(:,:,:,lat), beta,            hadv(:,:,:,lat) ,    &
                          fu(:,:,lat),      fv(:,:,lat),      t2(:,:,lat),         &
                          dqphy(:,:,lat),   omga(:,:,lat),    pdeld(:,:,jdry,n3) )
   end do
   call t_stopf ('mass_engy_fix')
!
   return
end subroutine mass_engy
!