module Trans_coef
!------------------------------------------------------------------------------------------------
! Purpose: Set IAP transformation coefficients
! Original version : MDINIT.f (IAP 9L)
! Reconstructed & revised : ZhangHe
! Completed : 2005.9.8
! Update : 2007.4.20, Zhanghe, 
!          2008.4.7, ZhangHe, added 'CALL DIAGHI' at sub. MDINIT
!          2008.4.27, ZhangHe, added sub. trans_af_mass
!          2008.5, WuJianping, parallel version
!          2008.6.11, ZhangHe, available for both serial & parallel
! Reviewed: ZhangHe, 2011-11-21
! Modified: Jiang Jinrong, 2012 October, for 2D parallel
! Reviewed: Zhang He, 2012-11-13
! Modified: Zhang He, 2013-01-25, new mp_send3d
!           Zhang He, 2013-03-12, removed transform of Q, Qliq, and Qice
!           Zhang He, 2013-02-21, moved calling of DIAGHI to stepon_init
!           Zhang He, 2013-03-28, removed UVW0 in MDINIT
!           Zhang He, 2013-04-01, use dynamic arrays WU et al.
!           Zhang He, 2013-04-07, moved calling of trans_IAP after DIAGHI in trans_af_mass
!           Zhang He, 2013-04-17, moved calling of DIAGPP and DIAGBB to sub. init_trans_IAP
!------------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,  only : NX, NY, NL, NZ, IB, IE, JB, JE, EX, period
   use mathconst, only : ZERO, HALF
   use pmgrid,    only : beglatdyn, endlatdyn, loc_JB, loc_JE,   &
                         beglatdynex, endlatdynex, beglev, endlev
   use spmd_utils, only: iam

   implicit none
   save
   public

   real(r8) :: PTU(NX,NY)     ! transformation coefficient from U to UT, UT = U * PTU
   real(r8) :: PTV(NX,NY)     ! transformation coefficient from V to VT, VT = V * PTV
   real(r8) :: PTT(NX,NY)     ! transformation coefficient from T'to TT, TT = T'* PTT
   real(r8) :: WW(NX,NY)     ! wjp 2008.5

!================================================================================================
CONTAINS
!================================================================================================

!================================================================================================
   SUBROUTINE MDINIT
!------------------------------------------------------------------------------------------------
! Purpose: INITIALIZE THE NORMAL DYNAMIC INTEGRATING CYCLE
!------------------------------------------------------------------------------------------------
      use stdatm,   only: DIAGBB
      use IAP_prog, only: WS, U, V, UVW0  !zhh 2007.11.7

	  IMPLICIT NONE
!----------------------------------Local workspace-----------------------------------------------
	  integer  :: I, J, K     ! loop index
!------------------------------------------------------------------------------------------------
      
!     SET BOUNDARY CONDITIONS AT THE SURFACE & MODEL TOP
      DO J = beglatdyn, endlatdyn
         DO I = 1 ,NX
            if(beglev.eq.1)     WS(I,1,J)  = ZERO !jjr
            if(endlev.eq.NL)    WS(I,NZ,J) = ZERO !jjr
         END DO
	  END DO
!!         DO J = beglatdyn, endlatdyn
!!            DO K = beglev ,endlev
!!               DO I = 1 ,NX
!!                  UVW0(I,K,J,1) = U (I,K,J)
!!                  UVW0(I,K,J,2) = V (I,K,J)
!!                  UVW0(I,K,J,3) = WS(I,K,J)
!!               END DO
!!            END DO
!!         END DO
!
!------------------ COMPUTE PRESSURE AT MODEL & INTERFACE SIGMA LAYERS ----------------------
!!      CALL DIAGPP
!--------------- SET MSA TEMPERATURE, GEOPOTENTIAL & CHARACTERISTIC PHASE-SPEED -------------
!!      CALL DIAGBB( 3 )
!--------------- COMPUTE GEOPOTENTIAL HEIGHT DEPARTURE ---------------------------
!!      CALL DIAGHI( 0 )            !zhh 2008.3.19
!!zhh 2013-02-21      CALL DIAGHI( 0 )            !zhh 2008.3.19
!---------------------------------------------------------------------------------
	  return
   end subroutine

!================================================================================================
   SUBROUTINE TRANSC
!------------------------------------------------------------------------------------------------
! Purpose: CALCULATE TRANSFORMATION COEFFECIENTS
!------------------------------------------------------------------------------------------------
      use physconst, only : RD, b0
      use IAP_prog,  only : PT
      use dynamics_vars,   only: T_FVDYCORE_STATE         !zhh 2013-01-25
      use dyn_internal_state,   only : get_dyn_state      !zhh 2013-01-25
#if (defined SPMD)
      use mod_comm, only: mp_send3d, mp_recv3d
      use pmgrid,     only : npr_y
#endif

	  IMPLICIT NONE
!----------------------------------Local workspace-----------------------------------------------
      type (T_FVDYCORE_STATE), pointer :: dyn_state   !zhh 2013-01-25
      real(r8) :: PT0, RDPT
	  integer  :: I, J, src, dest            ! loop index 
      integer  :: commyz     !zhh 2013-01-25
!------------------------------------------------------------------------------------------------
!
      dyn_state => get_dyn_state()
      commyz = dyn_state%grid%commyz
!
#if (defined SPMD)
      src = iam-1
      dest  = iam+1
      if ( mod(iam,npr_y) == 0 ) src = -1
      if ( mod(iam+1,npr_y) == 0 ) dest = -1
      call mp_send3d( commyz, dest, src, NX,  NY,1,                  &
                      1, NX, beglatdynex, endlatdynex,1,1,           &
                      1, NX, beglatdyn, beglatdyn, 1, 1, pt )
      call mp_recv3d( commyz, src, NX,  NY,1,                        &
                      1, NX, beglatdynex, endlatdynex,1,1,           &
                      1, NX, endlatdynex, endlatdynex, 1,1,pt )
#endif
      DO J = beglatdyn, endlatdyn
         if ( j >= JB .and. j <= JE ) then 
            DO I  = IB, IE
               PT0      = PT(I,J)
               RDPT     = RD   *   PT0
               PTU(I,J) = HALF * ( PT0 + PT(I-1,J) )
               PTV(I,J) = HALF * ( PT0 + PT(I,J+1) )
               PTT(I,J) = RDPT / b0
            END DO
!
         else if ( j == 1 ) then    ! at north pole
	        PT0 = PT(IB,1)
            DO I = IB, IE
               PTV(I,1) = HALF * ( PT0 + PT(I,JB) )  ! at the north polar
            END DO
            RDPT = RD * PT(IB,1)
            DO I = IB, IE
		       PTT(I,1) = RDPT / b0
	        END DO
            DO I = IB, IE
               PTU(I,1 )   = ZERO
            END DO
!
         else if ( j == NY ) then    ! at south pole
            RDPT = RD * PT(IB,NY)
            DO I = IB, IE
		       PTT(I,NY) = RDPT / b0
	        END DO
            DO I = IB, IE
               PTU(I,NY)   = ZERO
               PTV(I,NY)   = ZERO
            END DO
         end if
!
!   Set the spherical cyclicity condition
         call period ( PTU(1,J) )
         call period ( PTV(1,J) )
         call period ( PTT(1,J) )
	  END DO
!
	  RETURN
   END SUBROUTINE
!
!================================================================================================
   SUBROUTINE trans_IAP_tend(IC, DSU, DSV, DST)    
!------------------------------------------------------------------------------------------------
!  Purpose: BACK FROM PHYSICS ROUTINE TO DYNAMICS ROUTINE TO INITIALIZE THE DYNAMIC MODEL
!  Update : ZhangHe, 2007.12.15
!           ZhangHe, 2008.4.21
!------------------------------------------------------------------------------------------------
      use IAP_prog, only : P
      implicit none
!----------------------------------------Arguments-----------------------------------------------
	  integer,  intent(in)    :: IC
      real(r8), intent(inout) :: DSU(NX,beglev:endlev,beglatdynex:endlatdynex)    ! input sink & source tendency of U 
	  real(r8), intent(inout) :: DSV(NX,beglev:endlev,beglatdynex:endlatdynex)    ! input sink & source tendency of V
	  real(r8), intent(inout) :: DST(NX,beglev:endlev,beglatdynex:endlatdynex)    ! input sink & source tendency of T
!-------------------------------------Local workspace--------------------------------------------
      integer  :: I, J, K                         ! loop index
!------------------------------------------------------------------------------------------------
!
!! if Physics package can change P, then 
!------------------------------------------------------------------------------------------------
      CALL TRANSC   
!------------------------------------------------------------------------------------------------
! for DSU, DSV & DST
      DO J = beglatdyn, endlatdyn
         if ( j >= JB .and. j <= JE ) then 
            DO K = beglev ,endlev
               DO I = 1 ,NX       !! maybe need spherical cyclicity condition !!
                  DSU(I,K,J) = DSU(I,K,J) * PTU(I,J) 
                  DSV(I,K,J) = DSV(I,K,J) * PTV(I,J)
               END DO
            END DO
         else if ( j == 1 ) then    ! at north pole
            DO K = beglev ,endlev
               DO I = 1 ,NX
                  DSV(I,K,1) = DSV(IB,K,1) * PTV(IB,1)
                  DSU(I,K,1) = ZERO  
               END DO
            END DO
         else if ( j == NY ) then    ! at south pole
            DO K = beglev ,endlev
               DO I = 1 ,NX
                  DSU(I,K,NY)  = ZERO
                  DSV(I,K,NY)  = ZERO
               END DO
            END DO
         end if
      END DO
!
      DO J = beglatdyn, endlatdyn
         DO K = beglev ,endlev
            DO I = 1 ,NX
               DST(I,K,J) = DST(I,K,J) * PTT(I,J)
            END DO
         END DO
      END DO
!
      RETURN    
   END SUBROUTINE
!
!================================================================================================
   SUBROUTINE PHYSDM(IC, U, V, T)    !zhh 2008.4.21
!------------------------------------------------------------------------------------------------
!  Purpose: BACK FROM PHYSICS ROUTINE TO DYNAMICS ROUTINE TO INITIALIZE THE DYNAMIC MODEL
!------------------------------------------------------------------------------------------------
      use flexib,   only : NPHF
      use stdatm,   only : TB
      use IAP_prog, only : P, UT, VT, TT
!
	  implicit none
!----------------------------------------Arguments-----------------------------------------------
	  integer,  intent(in) :: IC
	  real(r8), intent(in) :: U(NX,beglev:endlev,beglatdynex:endlatdynex)
	  real(r8), intent(in) :: V(NX,beglev:endlev,beglatdynex:endlatdynex)
	  real(r8), intent(in) :: T(NX,beglev:endlev,beglatdynex:endlatdynex)
!-------------------------------------Local workspace--------------------------------------------
      real(r8) :: WU(NX,beglev:endlev,beglatdyn:endlatdyn)    
	  real(r8) :: WV(NX,beglev:endlev,beglatdyn:endlatdyn)
	  real(r8) :: WT(NX,beglev:endlev,beglatdyn:endlatdyn)
      integer  :: I, J, K, KK          ! loop index
!------------------------------------------------------------------------------------------------
!! if Physics package can change P, then 
!------------------------------------------------------------------------------------------------
      CALL TRANSC   
!------------------------------------------------------------------------------------------------
!
!     COMPUTE THE LATEST UT & VT & TT & QT
!
      DO J = beglatdyn, endlatdyn
         if ( j >= JB .and. j <= JE ) then 
            DO K = beglev ,endlev
               DO I = 1 ,NX
                  WU(I,K,J) = U(I,K,J) * PTU(I,J)
                  WV(I,K,J) = V(I,K,J) * PTV(I,J)
                  WT(I,K,J) = ( T(I,K,J)-TB(I,K,J) ) * PTT(I,J) 
               END DO
	        END DO
         else if ( j == 1 ) then    ! at north pole
            DO K = beglev ,endlev
               DO I = 1 ,NX
                  WU(I,K,J) = ZERO
                  WV(I,K,J) = V(I,K,J) * PTV(I,J)
                  WT(I,K,J) = ( T(I,K,J)-TB(I,K,J) ) * PTT(I,J) 
               END DO
	        END DO
         else if ( j == NY ) then    ! at south pole
            DO K = beglev ,endlev
               DO I = 1 ,NX
                  WU(I,K,J) = ZERO
                  WV(I,K,J) = ZERO
                  WT(I,K,J) = ( T(I,K,J)-TB(I,K,J) ) * PTT(I,J) 
               END DO
	        END DO
         end if
	  END DO
!
      IF ( NPHF.EQ.1 ) THEN    ! FILTER IS NEEDED DUE TO PHYSICS

! -------------------  FILTER DU & DV & DTT -----------------
         DO J = beglatdyn, endlatdyn
            DO K = beglev ,endlev
               DO I = 1 ,NX
                  WU(I,K,J)   = WU(I,K,J) - UT(I,K,J)   ! DU
                  WV(I,K,J)   = WV(I,K,J) - VT(I,K,J)   ! DV 
                  WT(I,K,J)   = WT(I,K,J) - TT(I,K,J)   ! DT
               END DO
		    END DO
	     END DO
!$DOACROSS LOCAL(KK,K)
         DO K = beglev ,endlev
            WW(:,beglatdyn:endlatdyn)=WT(:,K,beglatdyn:endlatdyn)
            CALL FILT2D( WW,0,1,IC )
            WT(:,K,beglatdyn:endlatdyn)=WW(:,beglatdyn:endlatdyn)
!
            WW(:,beglatdyn:endlatdyn)=WU(:,K,beglatdyn:endlatdyn)
            CALL FILT2D( WW,1,1,IC )
            WU(:,K,beglatdyn:endlatdyn)=WW(:,beglatdyn:endlatdyn)
!
            WW(:,beglatdyn:endlatdyn)=WV(:,K,beglatdyn:endlatdyn)
            CALL FILT2D( WW,1,2,IC )
            WV(:,K,beglatdyn:endlatdyn)=WW(:,beglatdyn:endlatdyn)
         ENDDO
!
!       RENEW UT & VT & TT
         DO J = beglatdyn, endlatdyn
            DO K = beglev ,endlev
               DO I = 1 ,NX
                  UT(I,K,J) = WU(I,K,J) + UT(I,K,J)
                  VT(I,K,J) = WV(I,K,J) + VT(I,K,J)
                  TT(I,K,J) = WT(I,K,J) + TT(I,K,J)
               END DO
		    END DO
	     END DO
!
      ELSE       ! FILTER IS NOT NEEDED  
         DO J = beglatdyn, endlatdyn
            DO K = beglev ,endlev
               DO I = 1 ,NX
                  UT(I,K,J)   = WU(I,K,J)
                  VT(I,K,J)   = WV(I,K,J)
                  TT(I,K,J)   = WT(I,K,J)
               END DO
		    END DO
	     END DO
      ENDIF
!
      RETURN
   END SUBROUTINE

!================================================================================================
   SUBROUTINE trans_antiIAP
!------------------------------------------------------------------------------------------------
!  Purpose: Convert dynamics variables to physics variables
!------------------------------------------------------------------------------------------------
      use stdatm,   only : TB
      use IAP_prog, only : P, UT, VT, TT, U, V, T

	  implicit none
!-------------------------------------Local workspace--------------------------------------------
      integer  :: I, J, K, KK                 ! loop index
!------------------------------------------------------------------------------------------------
!
!     SET FIELDS NEEDED  IN BACK & FOREWORD TRANSFORMATION
! ------------------------------------------------------------
      CALL TRANSC   
! ------------------------------------------------------------
!
!     DEDUCE   ( U,V,T,Q )  FROM   ( UT,VT,TT,QT ) 
      DO J = loc_JB, loc_JE
         DO K = beglev ,endlev
            DO I = 1 ,NX
               U(I,K,J) = UT(I,K,J) / PTU(I,J)  ! u=U/P
            END DO
		 END DO
	  END DO
      DO J = beglatdyn ,loc_JE
         DO K = beglev ,endlev
            DO I = 1 ,NX
               V(I,K,J) = VT(I,K,J) / PTV(I,J)  ! v=V/P
            END DO
		 END DO
	  END DO
      DO J = beglatdyn, endlatdyn
         DO K = beglev ,endlev
            DO I = 1 ,NX
               T(I,K,J) = TT(I,K,J) / PTT(I,J)  + TB(I,K,J) ! PTT=P*Rd/C0 
            END DO
		 END DO
	  END DO
!
      RETURN
   END SUBROUTINE
!================================================================================================
   SUBROUTINE trans_IAP
!------------------------------------------------------------------------------------------------
!  Purpose: Convert physics variables to dynamics variables
!  Completed: 2007.12.15
!------------------------------------------------------------------------------------------------
      use stdatm,   only : TB
      use IAP_prog, only : P, UT, VT, TT, U, V, T

	  implicit none
!-------------------------------------Local workspace--------------------------------------------
      integer  :: I, J, K                ! loop index
!------------------------------------------------------------------------------------------------
!
!     SET FIELDS NEEDED  IN BACK & FOREWORD TRANSFORMATION
! ------------------------------------------------------------
      CALL TRANSC   
! ------------------------------------------------------------
!
!     DEDUCE   ( U,V,T,Q )  FROM   ( UT,VT,TT,QT ) 
      DO J = beglatdyn, endlatdyn
         DO K = beglev ,endlev
            DO I = 1 ,NX
               UT(I,K,J) = U(I,K,J) * PTU(I,J)  
               VT(I,K,J) = V(I,K,J) * PTV(I,J)  
               TT(I,K,J) = ( T(I,K,J)-TB(I,K,J) ) * PTT(I,J)  
            END DO
		 END DO
	  END DO

      RETURN
   END SUBROUTINE

!================================================================================================
   SUBROUTINE trans_af_mass( ps, t3, q3, pcnst, n3, npt ) 
!------------------------------------------------------------------------------------------------
!  Purpose: Convert ps, t3, q3 to P, T, Q
!  Author: ZhangHe
!  Completed: 2008.4.27
!  Update: 2010.8, juanxiong He
!------------------------------------------------------------------------------------------------
      use Dyn_const, only: PMTOP        
      use IAP_prog,  only: P, PS2, Pstar1, Psa, PT, T, Q, Qliq, Qice
      use pmgrid,    only: plon, plev, plat, beglat, endlat
      use flexib,    only: Istar
      use stdatm,    only: PSB, P00, DIAGBB   
      use eul_control_mod, only: fixmas
	  implicit none

!----------------------------------------Arguments-----------------------------------------------
      integer,  intent(in) :: pcnst   ! number of advected constituents (including water vapor), juanxiong he
      integer,  intent(in) :: n3       ! current time level
      integer,  intent(in) :: npt      ! number of time levels in the dycore
      real(r8), intent(inout) :: t3(plon,beglev:endlev,beglat:endlat,npt)    ! temperature
      real(r8), intent(inout) :: q3(plon,beglev:endlev,pcnst,beglat:endlat,npt)   ! specific humidity, juanxiong he, 201008
      real(r8), intent(inout) :: ps(plon,beglat:endlat,npt)         ! surface pressure
!-------------------------------------Local workspace--------------------------------------------
      integer  :: I, J, Jp, Jd, K                ! loop index
!------------------------------------------------------------------------------------------------
!
! for ps
      DO Jd = beglatdyn, endlatdyn
         Jp = plat + 1 - Jd
         DO I = 1, plon
            PS2(I+EX,Jd) = ps(I,Jp,n3) / 100.0E0       
         END DO
         call period ( PS2(1,Jd) )
      END DO       
!
      DO J = beglatdyn, endlatdyn
         DO I = 1, NX
            P(I,J)   = PS2(I,J) - PMTOP
            Psa(I,J) = PS2(I,J) - PSB(I,J)
         END DO
      END DO       
!
      DO J = beglatdyn, endlatdyn
         DO I = 1, NX
            if (Istar == 1) then
               Pstar1(I,J) = P(I,J)
            else
               Pstar1(I,J) = Pstar1(I,J) * fixmas
            end if
            PT(I,J)  = sqrt( Pstar1(I,J)/P00 )
         END DO
      END DO       
!
! for T & Q
      DO Jd = beglatdyn, endlatdyn
         Jp = plat + 1 - Jd
         DO K = beglev, endlev
            DO I = 1, plon
               T   (I+EX,K,Jd) = t3(I,K,Jp,n3)       
!!               Q   (I+EX,K,Jd) = q3(I,K,1,Jp,n3) 
!!               Qliq(I+EX,K,Jd) = q3(I,K,2,Jp,n3)
!!               Qice(I+EX,K,Jd) = q3(I,K,3,Jp,n3)
            END DO
            call period( T   (1,K,Jd) )
!!            call period( Q   (1,K,Jd) )
!!            call period( Qliq(1,K,Jd) )
!!            call period( Qice(1,K,Jd) )
         END DO
      END DO 
!
! for UT, VT, TT      
!!      call trans_IAP
!------------------ COMPUTE PRESSURE AT MODEL & INTERFACE SIGMA LAYERS ----------------------
      CALL DIAGPP
!--------------- SET MSA TEMPERATURE, GEOPOTENTIAL & CHARACTERISTIC PHASE-SPEED -------------
      CALL DIAGBB( 3 )
!--------------- COMPUTE GEOPOTENTIAL HEIGHT DEPARTURE ---------------------------
      CALL DIAGHI( 0 )            

! for UT, VT, TT      
      call trans_IAP

      RETURN
   END SUBROUTINE

end module
