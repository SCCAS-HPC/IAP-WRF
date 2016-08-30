module Dyn_const
!------------------------------------------------------------------------------------------------
! Purpose: constants mainly used in dynamical calculations depend on model frame of resolution
! Original version: STFRAM.f (IAP 9L)
! Reconstructed to module : ZhangHe
! Completed : 2005.8.19
! Update    : 2006.3.20,  ZhangHe, modify the coefficient for flexible leaping-point zonal difference
!           : 2006.10.12, ZhangHe, add the constants for solving Q equation
!           : 2006.12.12, ZhangHe, change the dumb parameter of sub. STFRAM
!           : 2007.11.11, ZhangHe, delete the dumb parameter of sub. STFRAM
!                                  add calling sub. setfle
!           : 2008.04.15, ZhangHe, add wlat(NLAT), modify calculation of DLAT & DLON 
!           : 2008.5, WuJianping, parallel version
!           : 2008.6.9,  ZhangHe, available for both serial & parallel
!           : 2008.9.11,  ZhangHe, redefine hyam, hyai, hybm, hybi
!           : 2010.10,  Juanxiong He
!           : October 2012, Jiang Jinrong, for 2D parellel
!-----------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,  only: NX, NY, NL, NZ, JB, JE, IM
   use physconst, only: SDAY, RAD
   use flexib,    only: PHALFC, setfle           ! zhh  2007.11.11
   use mathconst, only: ZERO, ONE, TWO, FOUR, HALF, PI, FOURTH
   use pmgrid,    only: beglatdyn, endlatdyn, loc_JB, loc_JE,  &
                        beglatdynex,endlatdynex !!jjr
   use commap,    only: w      !zhh 2008.4.15
   use spmd_utils,       only: masterproc,iam    !zhh test

   implicit none

   save
   public

   real(r8) :: PMTOP      ! pressure at the top level of the model (unit: hPa)
   real(r8) :: SIG(NZ)    ! sigma value on model layer
   real(r8) :: SIGL(NL)   ! sigma value on interface layer
   real(r8) :: DSIG(NL)   ! vertical stepsizes (integer layer)     (unit: 1)
   real(r8) :: DSIGL(NZ)  ! vertical stepsizes (half    layer)     (unit: 1)
   real(r8) :: DSGHL(NL)  ! DSGHL(K) = DSIG(K-1)/DSIG(K)
   real(r8) :: DLAT       ! step length of latitude
   real(r8) :: DLON       ! step length of longtitude
   real(r8) :: FIM        ! FIM = DBLE( NLON )
   real(r8) :: RDLATI     ! RDLATI  = 2 / [IM*a*d(sita)]       ! a = RAD
   real(r8) :: RDLATI2    ! RDLATI2 = RDLATI * SINU(2) / SINV(1)
   real(r8) :: RDLAT      ! RDLAT   = 1 / (a*DLAT)
   real(r8) :: RDLON      ! RDLON   = 1 / (a*DLON)
   real(r8) :: RDLATH     ! RDLATH  = RDLAT / 2
   real(r8) :: RDLATQ     ! RDLATQ  = RDLAT / 4
   real(r8) :: OUX(NY)    ! OUX(j)  = 1 / [a*DLON*sin(sita(j))]
   real(r8) :: OUXH(NY)   ! OUXH(j) = OUX(j) / 2
   real(r8) :: OUXQ(NY)   ! OUXQ(j) = OUX(j) / 4
   real(r8) :: CUR(NY)    ! CUR(j)  = ctg(sita(j)) / a
   real(r8) :: CURQ(NY)   ! CURQ(j) = CUR(j) / 4
   real(r8) :: OVXQ(NY)   ! OVXQ(j) = 1 / [4*a*DLAT*sin(sita(j+1/2))]
   real(r8) :: RUP(NY)    ! RUP(j)  = sin(sita(j+1/2)) / [a*DLAT*sin(sita(j)]
   real(r8) :: RUM(NY)    ! RUM(j)  = sin(sita(j-1/2)) / [a*DLAT*sin(sita(j)]
   real(r8) :: RUPH(NY)   ! RUPH(j) = RUP(j) / 2
   real(r8) :: RUMH(NY)   ! RUMH(j) = RUM(j) / 2
   real(r8) :: RUPPH(NY)  ! RUPPH(j)= sin(sita(j+1)) / [2*a*sin(sita(j))]
   real(r8) :: RUMMH(NY)  ! RUMMH(j)= sin(sita(j-1)) / [2*a*sin(sita(j))]
   real(r8) :: RUPQ(NY)   ! RUPQ(j) = RUP(j) / 4
   real(r8) :: RUMQ(NY)   ! RUMQ(j) = RUM(j) / 4
   real(r8) :: RVP(NY)    ! RVP(j)  = sin(sita(j+3/2)) / [4*a*DLAT*sin(sita(j+1/2))]
   real(r8) :: RVM(NY)    ! RVM(j)  = sin(sita(j-1/2)) / [4*a*DLAT*sin(sita(j+1/2))]
   real(r8) :: RUPD(NY)   ! RUPD(j) = sin(sita(j+1/2)) / [4*sin(sita(j))]
   real(r8) :: RUMD(NY)   ! RUMD(j) = sin(sita(j-1/2)) / [4*sin(sita(j))]
   real(r8) :: FF(NY) 
   real(r8) :: FFQ(NY)
   real(r8) :: SGH(NL)
   real(r8) :: SGQ(NL)
   real(r8) :: RSGD(NL)
   real(r8) :: PSL000
   real(r8) :: SCLPR(NL)
   real(r8) :: RINTQ(NL)
   real(r8) :: DXYP(NY)
   real(r8) :: SINL(NY)
   real(r8) :: SINV(NY)    ! SINV(J) = sin[sita(j+1/2)]    
   real(r8) :: COSL(NY)    ! COSL(J) = sin[sita(j)]
   real(r8) :: COSLN(IM)
   real(r8) :: SINLN(IM)
   real(r8) :: COLAT(NY)
   real(r8) :: DXVPN(NY)
   real(r8) :: DXVPS(NY)
   real(r8) :: DXPVN(NY)
   real(r8) :: DXPVS(NY)
   real(r8) :: DXYV(NY)
   real(r8) :: DSNP
   real(r8) :: DSSP
   real(r8) :: DTDLN(NY)   !zhh 2007.8.20
   real(r8) :: DTDLT(NY)   !zhh 2007.8.20
   real(r8) :: wlat(NY)    ! wlat(j)=d(sita) * sin[sita(J)] ; weight of latitude
!=======================================================================
   real(r8) :: Ru          ! Ru = dt / (a*DLON)
   real(r8) :: Rv          ! Rv = dt / (a*DLAT)
   real(r8) :: RuF(NY)     ! RuF = dt / [a*DLON*sin(sita(j))]
   real(r8) :: RvF(NY)     ! RvF = dt / [a*DLAT*sin(sita(j))]
   real(r8) :: RuvF(NY)    ! RuvF = RuF*Rv
   real(r8) :: RuH(NY)     ! RuH = dt / [a*DLON*sin(sita(j+1/2))]
   real(r8) :: RvH(NY)     ! RvH = dt / [a*DLAT*sin(sita(j+1/2))]
   real(r8) :: RuvH(NY)    ! RuvH = RuH*Rv
   real(r8) :: Rw(NL)      ! Rw = dt / DSIGL
   real(r8) :: Rw2(NL)     ! Rw = dt / DSIG
!=======================================================================
   real(r8) :: GC(NY)      ! GC(J) = COSL(J)
   real(r8) :: DTDSG(NL)
   real(r8) :: alpha(NY)   ! coefficient for flexible leaping-point zonal difference (at J)
   real(r8) :: alpha2(NY)  ! coefficient for flexible leaping-point zonal difference (at J+1/2)
   real(r8) :: betaw(NY)   ! coefficient for flexible leaping-point zonal difference (at J)
   real(r8) :: betaw2(NY)  ! coefficient for flexible leaping-point zonal difference (at J+1/2)
   integer  :: IP(NX) 
   real(r8),allocatable :: ksa(:,:)   !jjr

!================================================================================================
CONTAINS
!================================================================================================

!================================================================================================
   SUBROUTINE STFRAM          ! zhh 2007.11.11
!-----------------------------------------------------------------------------------------------
!  Set constants mainly used in dynamical calculations depend on model frame of resolution
!-----------------------------------------------------------------------------------------------
      use pmgrid, only: plev, plevp
      use cam_control_mod, only: adiabatic        !zhh, 2012-10-24
!!      use hycoef

      IMPLICIT NONE
!------------------------------------Arguments--------------------------------------
!!      real(r8), intent(in) :: PHALFC(NZ)  ! delete by zhh 2006.12.11
!!      real(r8), intent(in) :: ZPS0  ! tandard pressure at sea level by zhh 2006.12.11
!----------------------------------Local workspace----------------------------------
!-----------------------------------------------------------------------
      real(r8) :: ZPS0     ! tandard pressure at sea level 
      real(r8) :: SINU(NY) ! , SINV(NY)                        
      real(r8) :: PIH,PI18,FRDAY,FJM,DEGLON,DEGLAT,PESCC,YV,FF0,YU         &  ! delete PS0
                 ,SNU,CSU,SNV,RSNV,RSNU,OUYJ,OVYJ,SNVJ,SNVI,AA,OUXJ,CURJ   &   
                 ,RUPJ,RUMJ,FFJ,RUPQJ,RUMQJ,DSGK,PLAYC,RLATJ,DAN,DAS,DAV   &
                 ,DXYN,DXYS,DXYPN,DXYPS,APV,DAP,ALON,FIRST,DTIMEQ,Y1
      real(r8) :: bn, en, bs, es, inc, wt
      INTEGER  :: I, J, K, KP
!---------------------------------------------------------------------------------------------
      PIH       = PI      * HALF 
      PI18      = PI      / 180.0D0
      FRDAY     = TWO*PI  / SDAY         ! FRDAY: earth's angular velocity
      FIM       = DBLE(  IM  )
      FJM       = DBLE(  JE  )   
      DEGLON    = 360.0E0 / FIM          ! DEGLON = d(lamda)
      DEGLAT    = 180.0E0 / FJM          ! DEGLAT = d(sita)
!!      DLAT      = DEGLAT  * PI18
!!      DLON      = DEGLON  * PI18
      DLAT      = PI / dble(JE)
      DLON      = TWO * PI / dble(IM)
!! ======================== zhh 2008.4.15 =========================
!!      do j = loc_JB, loc_JE
      do j = JB, JE
         wlat(j) = PI / dble(JE) * sin( dble(j-1)*PI/dble(JE) )
      end do
      wlat(1)  = FOURTH * sin( HALF*DLAT ) * DLAT
      wlat(NY) = FOURTH * sin( HALF*DLAT ) * DLAT
!
      wt       = 0.0
!!      do j = beglatdyn, endlatdyn
      do j = 1, NY
         wt = wt + wlat(j)
      end do
!!      DO j = beglatdyn, endlatdyn
      DO J = 1, NY
         wlat(j) = wlat(j) * TWO / wt 
         w(j)    = wlat(j)  
      END DO
!! ======================== zhh 2008.4.15 =========================
!
      call setfle           !zhh 2007.11.11
! ====================== zhh ===========================
!!      PHALFC(NZ) = ZPS0          ! zhh 2006.12.12
      ZPS0  = PHALFC(NZ)          ! zhh 2007.11.11
!!zhh 2008.4.22      PMTOP = PHALFC(1)
      PMTOP = ZERO      !zhh 2008.4.22
      PESCC = ZPS0 - PMTOP   ! zhh 2007.11.7
! ===================== 2007.11.7 ===========================
!!zhh 2008.4.22      SIG (1)    = ZERO
!!      do K  = 2 ,NZ
      do K  = 1 ,NZ
         SIG (K)   = (PHALFC(K) - PMTOP) / PESCC
      end do
      do K  = 1 ,NL
         KP        = K + 1
         DSIG(K)   = SIG(KP) - SIG(K)
         SIGL(K)   = HALF * (SIG(KP) + SIG(K))
      end do
      do K  = 2 ,NL                       
         DSGHL(K)  = DSIG(K-1) / DSIG(K)       ! for sub. AVNEGQ
         DSIGL(K)  = SIGL(K) - SIGL(K-1)       ! add by zhh 2006.10.7
      end do
      DSIGL(1)  = SIGL(1) - SIG(1)             ! add by zhh 2006.10.7
      DSIGL(NZ) = SIG(NZ) - SIGL(NL)           ! add by zhh 2006.10.7
! ========================= zhh 2008.9.11 ===============================
!!      do K = 1, NZ
!!         hyai(K) = ZERO
!!         hybi(K) = SIG(K)
!!      end do
!!      do K = 1, NL
!!         hyam(K) = ZERO
!!         hybm(K) = SIGL(K)
!!      end do
! ========================= zhh 2008.9.11 ===============================
!
      RDLAT     = ONE   / (RAD * DLAT)  
      RDLON     = ONE   / (RAD * DLON)
      RDLATI    = RDLAT * TWO / FIM   
      RDLATH    = HALF  * RDLAT             ! add by zhh
      RDLATQ    = RDLAT / FOUR
!
!jjr      do J  = beglatdyn ,loc_JE
          do j=1,JE           !jjr
!***************************************************************************************
         YV        = DLAT  * (DBLE( J ) - HALF)
         YU        = DLAT  *  DBLE( J-1 )
         SINV(J)   = SIN( YV )             !  SINV(J) = sin¦È(j+1/2)
         SINU(J)   = SIN( YU )       !zhh 2005.8.18  SINU(J) = sin¦È(j)
!***************************************************************************************
      end do
      SINV(NY)  = ZERO
      SINU(NY)  = ZERO          !zhh 2007.8.23

      FF0       = FRDAY + FRDAY
!
#if (defined SPMD)
! Latitude distribution for plat=16 on 3 processors.             !wjp 2007.05
!  procId:      iam=0     |      iam=1     |      iam=2          !wjp 2007.05
!  CAM:    01 02 03 04 05 | 06 07 08 09 10 | 11 12 13 14 15 16   !wjp 2007.05
!  IAP:    16 15 14 13 12 | 11 10 09 08 07 | 06 05 04 03 02 01   !wjp 2007.05
!                left     |     current    |      right          !wjp 2007.05
!JJR      call mpi_move_right( SINV(beglatdyn), SINV(endlatdyn+1), 1)
!jjr      call mpi_move_left(SINV(endlatdyn), SINV(beglatdyn-1), 1)
!jjr      call mpi_move_right( SINU(beglatdyn), SINU(endlatdyn+1), 1)
!jjr      call mpi_move_left(SINU(endlatdyn), SINU(beglatdyn-1), 1)
#endif
!
!jjr      do J  = beglatdyn, endlatdyn
          do J=1,NY
         IF (JB .le. J .and. J .le. JE) then
!***************************************************************************************
            YU        = DLAT  *  DBLE( J-1 )
            SNU       = SIN( YU )
            CSU       = COS( YU )
            SNV       = SINV( J  )
!***************************************************************************************
            FF  (J)   = FF0   *  CSU
            RSNV      = ONE   / (RAD  * SNV)
            RSNU      = ONE   / (RAD  * SNU)
            CUR (J)   = CSU   *  RSNU
            OUX (J)   = RSNU  /  DLON
            OVXQ(J)   = RSNV  / (DLON * FOUR)
            OUYJ      = RSNU  /  DLAT
            OVYJ      = RSNV  /  DLAT
            SNVJ      = SINV(J-1)
            SNVI      = SINV(J+1)
            RUP (J)   = SNV   *  OUYJ
            RUPPH(J)  = SINU(J+1) *  OUYJ * HALF         !zhh 2005.8.18 
            RUM (J)   = SNVJ  *  OUYJ
            RUMMH(J)  = SINU(J-1) *  OUYJ * HALF         !zhh 2005.8.18
            RVP (J)   = SNVI  *  OVYJ / FOUR
            RVM (J)   = SNVJ  *  OVYJ / FOUR
         ELSE
            CUR (J)   = ZERO
            OUX (J)   = ZERO
            OVXQ(J)   = ZERO
            RUP (J)   = ZERO
            RUM (J)   = ZERO
            RVP (J)   = ZERO
            RVM (J)   = ZERO
       	    RUPPH(J)  = ZERO                            !zhh 2005.8.18
            RUMMH(J)  = ZERO                            !zhh 2005.8.18
         ENDIF
      end do
!wjp 2007.06      RDLATI2 = RDLATI * SINU(2) / SINV(1)      !zhh 2005.8.20
       RDLATI2 = RDLATI * SIN(DLAT)/SIN(DLAT/2)  ! SINU(2)=SIN(DLAT), SINV(1)=SIN(DLAT/2)
      
      FF  (1 )  = + FF0
      FF  (NY)  = - FF0
      AA        = ONE   / (RAD  * SINV(1))
      OVXQ(1 )  = AA    / (DLON * FOUR )
      OVYJ      = AA    /  DLAT
      RVP (1 )  = OVYJ  *  SINV(JB) / FOUR

!jjr      do J  = beglatdyn, endlatdyn
      do J=1,NY
         OUXJ      = OUX(J)
         CURJ      = CUR(J)
         RUPJ      = RUP(J)
         RUMJ      = RUM(J)
         FFJ       = FF (J)
         OUXH(J)   = OUXJ  * HALF
         RUPH(J)   = RUPJ  * HALF
         RUMH(J)   = RUMJ  * HALF
         FFQ (J)   = FFJ   / FOUR
         OUXQ(J)   = OUXJ  / FOUR
         CURQ(J)   = CURJ  / FOUR
         RUPQJ     = RUPJ  / FOUR
         RUMQJ     = RUMJ  / FOUR
         RUPQ(J)   = RUPQJ
         RUMQ(J)   = RUMQJ
         RUPD(J)   = RUPQJ / RDLAT
         RUMD(J)   = RUMQJ / RDLAT
      end do
!
      do K  = 1 ,NL
         DSGK      = ONE   / DSIG(K)      
         SGH (K)   = DSGK  * HALF
         SGQ (K)   = DSGK  / FOUR
         RSGD(K)   = RDLAT * SIGL(K)
      end do

! set coefficient for flexible leaping-point zonal difference (zhh 2005.8. 26, modify: 2006.3.20)
      bn  = 18.0
      en  = 58.0
      bs  = 180.0 - bn  ! 162.0
      es  = 180.0 - en  ! 122.0
      inc = 1.0 / (en - bn)

!jjr      do J = beglatdyn, endlatdyn
      do j=1,NY
         Y1 = DBLE(J-1) * DEGLAT
         if (Y1 <= bn .or. Y1 >= bs) then
            alpha(J) = 0.0            
         else if (bn < Y1 .and. Y1 < en)  then
            alpha(J) = (Y1 - bn) * inc
         else if (en <= Y1 .and. Y1 <= es) then
            alpha(J) = 1.0
         else if (es < Y1 .and. Y1 < bs) then
            alpha(J) = (bs - Y1) * inc
         end if
      end do
      
!jjr      do J = beglatdyn, endlatdyn
      do J=1,NY
         Y1 = DBLE(J-0.5) * DEGLAT
         if (Y1 <= bn .or. Y1 >= bs) then
            alpha2(J) = 0.0            
         else if (bn < Y1 .and. Y1 < en)  then
            alpha2(J) = (Y1 - bn) * inc
         else if (en <= Y1 .and. Y1 <= es) then
            alpha2(J) = 1.0
         else if (es < Y1 .and. Y1 < bs) then
            alpha2(J) = (bs - Y1) * inc
         end if
      end do

!jjr      do J = beglatdyn, endlatdyn
      do J=1,NY
         betaw(J)  = 1.0 - alpha(J)
         betaw2(J) = 1.0 - alpha2(J) 
      end do
!      
!!==========================================================
!!	  do J = beglatdyn, endlatdyn        ! for a test
!!	     alpha(J)  = 1.0 
!!         betaw(J)  = 0.0
!!         alpha2(J) = 1.0  
!!	     betaw2(J) = 0.0
!!	  end do
!!==========================================================

!     THE FOLLOWINGS ARE SPECIAL FOR PHYSICS' ROUTINES
!
      do K  = 2 ,NL
         RINTQ(K)  = (SIGL(K) - SIG(K)) / (SIGL(K) - SIGL(K-1))
      end do
      PSL000    = ZPS0
      do K  = 1 ,NL
         PLAYC     = PMTOP + SIGL(K)*PESCC
         SCLPR(K)  = PLAYC / ZPS0
      end do
      do J  = beglatdyn, endlatdyn
         RLATJ     = PIH   - DBLE( J-1 )*DLAT
         COLAT(J)  = PIH - RLATJ
         SINL (J)  = SIN( RLATJ )
         COSL (J)  = COS( RLATJ )    ! COSL(J) = sin(sita(j))
      end do
      COSL( 1)  =  ZERO
      COSL(NY)  =  ZERO
      SINL( 1)  =  ONE
      SINL(NY)  = -ONE
!
#if (defined SPMD)
! Latitude distribution for plat=16 on 3 processors.             !wjp 2007.05
!  procId:      iam=0     |      iam=1     |      iam=2          !wjp 2007.05
!  CAM:    01 02 03 04 05 | 06 07 08 09 10 | 11 12 13 14 15 16   !wjp 2007.05
!  IAP:    16 15 14 13 12 | 11 10 09 08 07 | 06 05 04 03 02 01   !wjp 2007.05
!                left     |     current    |      right          !wjp 2007.05
!      call mpi_move_right( COSL(beglatdyn), COSL(endlatdyn+1), 1)
!      call mpi_move_left(COSL(endlatdyn), COSL(beglatdyn-1), 1)
#endif
!
 
!jjr      do J  = loc_JB,loc_JE
      do J=JB,JE
         DXYP(J) = HALF*(HALF*(RAD*DLON)*(COSL(J-1) + COSL(J+1)+COSL(J)*2.0))*(RAD*DLAT) !(*1)
      end do
      DXYP(1 )   = (HALF*(RAD*DLON)*(COSL(2)+COSL(1)))*(RAD*DLAT) / FOUR
      DXYP(NY)   = (HALF*(RAD*DLON)*(COSL(JE)+COSL(NY)))*(RAD*DLAT)/FOUR

!jjr      do J  = beglatdyn ,loc_JE
       do J=1,JE
         DXYV(J)   = HALF*(RAD*DLON*(COSL(J+1)+COSL(J)))*(RAD*DLAT) !(*2)
      end do
!
#if (defined SPMD)
! Latitude distribution for plat=16 on 3 processors.             !wjp 2007.05
!  procId:      iam=0     |      iam=1     |      iam=2          !wjp 2007.05
!  CAM:    01 02 03 04 05 | 06 07 08 09 10 | 11 12 13 14 15 16   !wjp 2007.05
!  IAP:    16 15 14 13 12 | 11 10 09 08 07 | 06 05 04 03 02 01   !wjp 2007.05
!                left     |     current    |      right          !wjp 2007.05
!      call mpi_move_right( DXYP(beglatdyn), DXYP(endlatdyn+1), 1)
#endif
!
!jjr      do J  = loc_JB,min(JE-1,endlatdyn)
!zhh          do J=1,JE
      do J=JB,JE-1
         DAN       = DXYP(J  )
         DAS       = DXYP(J+1)
         DAV       = DAN + DAS
         DXYN      = DAN / DAV
         DXYS      = DAS / DAV
         DXPVN(J)  = DXYN
         DXPVS(J)  = DXYS
      end do
      DXYPN     = DXYP(1 )
      DXYPS     = DXYP(JB)
      APV       = DXYPN+DXYPN + DXYPS
      DXPVN(1)  = DXYPN/APV
      DXPVS(1)  = DXYPS/APV
      DXYPN     = DXYP(JE)
      DXYPS     = DXYP(NY)
      APV       = DXYPS+DXYPS + DXYPN
      DXPVN(JE) = DXYPN/APV
      DXPVS(JE) = DXYPS/APV
      DXPVN(NY) = ZERO
      DXPVS(NY) = ZERO
!
#if (defined SPMD)
! Latitude distribution for plat=16 on 3 processors.             !wjp 2007.05
!  procId:      iam=0     |      iam=1     |      iam=2          !wjp 2007.05
!  CAM:    01 02 03 04 05 | 06 07 08 09 10 | 11 12 13 14 15 16   !wjp 2007.05
!  IAP:    16 15 14 13 12 | 11 10 09 08 07 | 06 05 04 03 02 01   !wjp 2007.05
!                left     |     current    |      right          !wjp 2007.05
!      call mpi_move_left(DXYV(endlatdyn), DXYV(beglatdyn-1), 1)
#endif
!
!      do J  = loc_JB,loc_JE
       do J=JB,JE
         DAN       = DXYV(J-1)
         DAS       = DXYV(J  )
         DAP       = DAN + DAS
         DXYN      = DAN / DAP
         DXYS      = DAS / DAP
         DXVPN(J)  = DXYN  
         DXVPS(J)  = DXYS
      end do
      DXVPN(1)  = ZERO
      DXVPS(1)  = ZERO
      DXVPN(NY) = ZERO
      DXVPS(NY) = ZERO
 
      do I  = 1 ,IM
         ALON      = -PI + DBLE( I-1 )*DLON
         COSLN(I)  = COS( ALON )
         SINLN(I)  = SIN( ALON )
      end do
!jjr
     allocate(ksa(NX,beglatdynex:endlatdynex))
      do j = beglatdynex,endlatdynex
         do i = 1, NX
            if (adiabatic) then
               ksa(i,j) = 0.
            else
               ksa(i,j) = 0.1D0
            end if
         end do
      end do

      RETURN
   end subroutine STFRAM

!===============================================================================================
   subroutine CONPDA(FIRST)   ! modified by zhh  2006.10.10
!-----------------------------------------------------------------------------------------------
!     FOR THE USAGE OF SUB.MPDATA & VPDATA
!    DATA NONOS, IORD, ISOR, EP/0,3,3,1.E-10/
!-----------------------------------------------------------------------------------------------
      use flexib, only: DTIMEQ
!      use vapor,  only: ISOR, IORD
!------------------------------------Arguments--------------------------------------------------
      real(r8), intent(in) :: FIRST   ! initial control parameter , see dyfram.f
!----------------------------------Local workspace----------------------------------------------
      integer :: I, J, K              ! loop index
!----------------------------------------------------------------------------------------------- 
      IF ( FIRST.LE.ZERO ) THEN
!         IF ( ISOR.EQ.3 ) IORD = MAX0(IORD,3)
!     ISOR = 3: third order modification of water vapor scheme ; 
!     IORD: the times of iterations in water vapor scheme
         do I = 1 ,NX
            IP(I)  = MOD( I+IM/2-1,IM ) + 1
         end do
!=============================== zhh ====================================
! for new scheme
         DSNP      = TWO * DTIMEQ * RDLATI      !at the north polar
!*        DSNP      = 4*dt / [IM*a*d(sita)] 
         DSSP      = DSNP                       !at the south polar
! -----------------------------------------------------------------------
! for old scheme
!!         DSNP      = DXYP(JB) / (FIM*DXYP(01))  !at the north polar
!!         DSSP      = DXYP(JE) / (FIM*DXYP(NY))  !at the south polar
!!*        DSNP = DSSP =(about) 8.0 
!============================ 2007.12.5 =================================
!jjr         do J= beglatdyn, endlatdyn
        do j=1,NY
            GC(J)   = COSL(J)
         end do
      ENDIF

! -------------------- add by zhh 2006.9.5 ----------------------
      Ru = DTIMEQ / (RAD*DLON)
      Rv = DTIMEQ / (RAD*DLAT)

!jjr      do J = loc_JB, loc_JE
     do j=JB,JE
         RuF(J)  = Ru / COSL(J) 
         RvF(J)  = Rv / COSL(J) 
         RuvF(J) = Rv * RuF(J) 
      end do
!jjr      do J = beglatdyn, loc_JE
      do J=1,JE
         RuH(J)  = Ru / SINV(J)
         RvH(J)  = Rv / SINV(J) 
         RuvH(J) = Rv * RuH(J) 
      end do
! -------------------- delete by zhh 2006.9.5 ----------------------
!jjr      do J  = beglatdyn ,loc_JE
      do J=1,JE
         DTDLT(J)  = DTIMEQ / (RAD*DLAT) 
      end do
!      do J  = loc_JB,loc_JE
      do J=JB,JE
         DTDLN(J)  = DTIMEQ / (RAD*COSL(J)*DLON)
      end do
      DTDLT(NY) = ZERO     !2007.12.6
      DTDLN(1)  = ZERO
      DTDLN(NY) = ZERO
!-------------------------------------------------------------------
      do K  = 1 ,NL
         DTDSG(K)  = DTIMEQ / DSIG(K)
         Rw(K)  = DTIMEQ / DSIGL(K)       ! zhh 2006.10.11
         Rw2(K) = DTIMEQ / DSIG(K)        ! zhh 2006.10.12
      end do
!
      RETURN
   end subroutine CONPDA
end module Dyn_const
