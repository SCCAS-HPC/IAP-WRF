SUBROUTINE tend_lin
!------------------------------------------------------------------------------------------------
! Purpose: Compute the tendence of P, T, U, V ( Linear(adaption) part )
! Author : ZhangHe
! Completed : 2005.9.1
! Update    : 2006.2.22, Zhanghe, revise to linear part
!             2007.4.23, ZhangHe, 1) removed dumb parameter 'NCTDCB'
!                                 2) changed sub. name from 'ptuvtend1' to 'tend_lin'
!                                 3) changed 'LADDSS' to 'NADDSS'
!             2007.10.11, ZhangHe, modified the calculation of TOO0
!             2008.4.7, ZhangHe, 1) moved statement of Kstd to module flexib
!                                2) modified use of kappa
!             2008.5, WuJianping, parallel version
!             2008.6.11, ZhangHe, available for both serial & parallel
!             2011.05.03, ZhangHe, let deltac = 0 if adiabatic run
! Modified: Jiang Jinrong, 2012 October, for 2D parallel
! Reviewed: Zhang He, 2012-11-13
! Modified: Zhang He, 2013-01-28, new mp_send3d, dyn_state was added
!           Zhang He, 2013-02-06, removed redundant variables in use only statement
! Reviewed: Zhang He, 2013-03-21
! Modified: Zhang He, 2013-04-01, use dynamic arrays DGH et al.
!------------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,  only: NX, NY, NM, NL, NZ, JB, JE, IB, IE, period
   use Dyn_const 
   use stdatm,    only: DPALIB, NA, CB, CBB, PSB, deltac, H0B, TB     !zhh 2007.4.23
   use IAP_prog,  only: PT, TT, UT, VT, WST, Psa, U, WS, WPV, PLY, GHI,     &
                        Pstar1, Pstar2, GZsm, deltap, P, V
   use tendency,  only: DPsa, DT, DU, DV, SU, SV, ST, NADDSS
   use physconst, only: b0, CAPA
   use mathconst, only: ZERO, HALF, ONE, TWO, FOURTH
   use flexib,    only:  b1, b21, b22, c10, kappa, kappa2, IBCFFT, Kstd
   use pmgrid,    only: beglatdyn, endlatdyn, loc_JB, loc_JE,beglev,endlev,plon,myid_y, &
                    endlevp,npr_z,npr_y,beglatdynex,endlatdynex,plat,plev,myid_z
   use stdatm,    only: P00
   use spmd_utils, only: masterproc, iam
   use cam_control_mod, only: adiabatic        !zhh, 2011-11-19
   use dynamics_vars,   only: T_FVDYCORE_STATE         !zhh 2013-01-28
   use dyn_internal_state,   only : get_dyn_state      !zhh 2013-01-28
   use perf_mod, only : t_startf, t_stopf              !zhh 2013-02-06
#if ( defined SPMD )
   use mod_comm, only: mp_send3d, mp_recv3d,mp_send3d_2, mp_recv3d_2
   use spmd_dyn          , only:  comm_z,zdist,lreq,rreq
   use mpishorthand, only: mpir8
   use parutilitiesmodule, only: parcollective,sumop,pargatherreal
#endif


   IMPLICIT NONE
!----------------------------------Local workspace-----------------------------------------------
   real(r8), allocatable :: WW(:,:)                  ! wjp  2008.5
   real(r8), allocatable :: DGH(:,:,:)               ! d(GHI) 
   real(r8), allocatable :: GHI1(:,:,:)              ! d(GHI) 
   real(r8), allocatable :: TT1(:,:,:,:)             ! jjr 2012-07-16
   real(r8), allocatable :: TT0(:,:,:,:)             ! jjr 2012-07-16 
   real(r8), allocatable :: PXW(:,:), PXW2(:,:)      ! auxiliary variables
   real(r8), allocatable :: PYP(:,:), PYM(:,:)       ! auxiliary variables
   real(r8), allocatable :: rhos(:,:)                ! surface density of atmosphere
   real(r8), allocatable :: pDrho(:,:)               ! pDrho = P'sa / rhos
   real(r8), allocatable :: DpDrX(:,:), DpDrY(:,:)   ! zonal & meridional gradient of pDrho
   real(r8), allocatable :: DIV(:,:,:)               ! divergence of (PT*VT)
   real(r8), allocatable :: Dsa(:,:)                 ! Dsa = DsaX + DsaY , diffusion of P'sa
   real(r8), allocatable :: PY(:,:,:)                ! meridional component of pressure-gradient force
   real(r8), allocatable :: PX(:,:,:)                ! zonal component of pressure-gradient force
!
   real(r8) :: DIVX, DIVY                   ! divergence of (PT*VT)
   real(r8) :: D1(NL), DWK(NL) , DPS1       ! auxiliary variables
   real(r8) :: D2(NL), D3(NL) !jjr
   real(r8) :: DsaX, DsaY                   ! Dsa = DsaX + DsaY , diffusion of P'sa
   real(r8) :: DsaY1(NX), DPsa1
   real(r8) :: WK5                          ! WK5 = 1 / PT 
   real(r8) :: DXP0, DXP1, DXP2, DXP3       ! zonal difference of pes**
   real(r8) :: DYP0, DYP2                   ! meridional difference of pes**
   real(r8) :: OO1, OO2X, OO2Y, TOO0        ! terms of DTT/Dt
   real(r8) :: O1P, O2P, OOP, TOOP          ! terms of DTT/Dt at polar
   real(r8) :: DTTP                         ! tendency of TT at polar
   real(r8) :: PV1(NX)                      ! divergence of (PT¡¤VT) at polar
   real(r8) :: W1(NL+1), WPV1(NL),W11       ! auxiliary variables
   real(r8) :: WSK(NL+1), D1K, DSK,WSK1     ! auxiliary variables
   real(r8) :: OPCT0, OPCT2, OPCT3, OPCT4, OPCT5   ! auxiliary variables
   real(r8) :: OUXAXP, OUXAXP2, OYAYP, OYVS ! auxiliary variables
   real(r8) :: GH0, GH1, GH2, GH3, GH4      ! the average value of GHI
   real(r8) :: PY1, PY2                     ! meridional component of pressure-gradient force
   real(r8) :: PX1, PX2                     ! zonal component of pressure-gradient force
   real(r8) :: FS0, FS1, FS2, FS3           ! Coriolis parameter 
   real(r8) :: fstar0, fstar1               ! Coriolis parameter at polar
   real(r8) :: FSU, FSV                     ! work terms of Coriolis
   real(r8) :: WPK, WPJ, WPI                ! weight as pressure
   real(r8) :: time_begin, time_end,ttt
   integer  :: KPK, KPI,src,dest,src1,dest1
   integer  :: KK, I, J, K,ierr                  ! loop index
   integer  :: tmp1,num1,recs(npr_z),disps(npr_z)  !test
!jjr
   real(r8), allocatable :: sbuf(:,:)
   real(r8), allocatable :: sbuf1(:,:)
   real(r8), allocatable :: rbuf(:,:)
   real(r8), allocatable :: rbuf1(:,:)
   real(r8):: temp1(IB:IE),temp2(IB:IE) !for vec
!jjr
   type (T_FVDYCORE_STATE), pointer :: dyn_state   !zhh 2013-01-28
   integer  :: commyz     !zhh 2013-01-28

!----------------------------------------------------------------------------------------------
! allocate arrays
   allocate ( WW(NX,NY) )
   allocate ( DGH(NX,NL,beglatdyn:endlatdyn) )
!!   allocate ( DGH(NX,beglev:endlev,beglatdyn:endlatdyn) )
   allocate ( GHI1(NX,NZ,beglatdyn:endlatdyn) )
   allocate ( TT1(IB:IE,beglatdyn:endlatdyn,2,NL) )
   allocate ( TT0(IB:IE,beglatdyn:endlatdyn,2,beglev:endlev) )
   allocate ( PXW(NX,beglatdyn:endlatdyn) )
   allocate ( PXW2(NX,beglatdyn:endlatdyn) )
   allocate ( PYP(NX,beglatdyn:endlatdyn) )
   allocate ( PYM(NX,beglatdyn:endlatdyn) )
   allocate ( rhos(NX,beglatdynex:endlatdynex) )
   allocate ( pDrho(NX,beglatdynex:endlatdynex) )
   allocate ( DpDrX(NX,beglatdynex:endlatdynex) )
   allocate ( DpDrY(NX,beglatdynex:endlatdynex) )
   allocate ( DIV(IB:IE,beglev:endlev,beglatdyn:endlatdyn) )
!!   allocate ( Dsa(NX,NY) )
   allocate ( Dsa(NX,beglatdyn:endlatdyn) )
!!   allocate ( PY(NX,NL,NY) )
   allocate ( PY(NX,beglev:endlev,beglatdyn:endlatdyn) )
!!   allocate ( PX(NX,NL,NY) )
   allocate ( PX(NX,beglev:endlev,beglatdyn:endlatdyn) )
   allocate ( sbuf(NX,2*(endlev-beglev+1)+4) )
   allocate ( sbuf1(NX,4*(endlev-beglev+1)+4) )
   allocate ( rbuf(NX,2*(endlev-beglev+1)+4) )
   allocate ( rbuf1(NX,4*(endlev-beglev+1)+4) )
!
   dyn_state => get_dyn_state()
   commyz = dyn_state%grid%commyz

!  *******************************************************************
!  *******************************************************************
!  **************        COMPUTE DPsa & DT        ********************
!  *******************************************************************
!  *******************************************************************

!------------------- COMPUTE FACTORS RELATED TO (I,J) ONLY ----------------------
!
   call t_startf('tend_lin')
#if (defined SPMD)
!jjr src dest,right

   num1=endlev-beglev+1
   sbuf(:,1)=pt(:,endlatdyn)
   sbuf(:,2)=ksa(:,endlatdyn)
   sbuf(:,3)=p(:,endlatdyn)
   sbuf(:,4)=Pstar2(:,endlatdyn)
   do kk=beglev,endlev
      sbuf(:,5+kk-beglev)=vt(:,kk,endlatdyn)
      sbuf(:,5+num1+kk-beglev)=v(:,kk,endlatdyn)
   end do
   call mpi_move_left(sbuf,rbuf,NX*(2*num1+4))

!        sbuf1(1:NX,1)=pt(1:NX,beglatdyn)
   sbuf1(1:NX,1)=ksa(1:NX,beglatdyn)
   sbuf1(1:NX,2)=p(1:NX,beglatdyn)
   sbuf1(1:NX,3)=Pstar2(1:NX,beglatdyn)
   sbuf1(1:NX,4)=Pstar1(1:NX,beglatdyn)

   do kk=beglev,endlev
      sbuf1(1:NX,5+kk-beglev)=tt(1:NX,kk,beglatdyn)
      sbuf1(1:NX,5+num1+kk-beglev)=deltap(1:NX,kk,beglatdyn)
   end do
   call mpi_w_a(lreq)
!       tmp1 = 0
!         do while(tmp1.eq.0)
!          call sleep(2)
!       enddo

   if (myid_y.ne.(npr_y-1)) then
      pt(:,beglatdynex)=rbuf(:,1)
   end if
!right
   src = iam-1
   dest  = iam+1
   if ( mod(iam,npr_y) == 0 ) src = -1
   if ( mod(iam+1,npr_y) == 0 ) dest = -1
!jjr src1,dest1,left
   src1 = iam+1
   dest1  = iam-1
   if ( mod(iam,npr_y) == 0 ) dest1 = -1
   if ( mod(iam+1,npr_y) == 0 ) src1 = -1

#endif

   DO J = loc_Jb,loc_JE 
      DO I = 1 ,NX
         PYP(I,J)   = RUPH(J) * (PT(I,J+1)+PT(I,J)) !RUPH(J) = sin¦È(j+1/2) / (2*a*¦¤¦È*sin¦È(j))
         PYM(I,J)   = RUMH(J) * (PT(I,J)+PT(I,J-1)) !RUMH(J) = sin¦È(j-1/2) / (2*a*¦¤¦È*sin¦È(j))
         rhos(i,j)  = 1 / H0B(i,j)
         pDrho(i,j) = Psa(i,j) / rhos(i,j)
      ENDDO
      DO I = IB,IE
         PXW(I,J)   = OUXH(J) * (PT(I,J)   + PT(I-1,J)) !OUXH(J) = 1 / (2*a*¦¤¦Ë*sin¦È(j))
         PXW2(I,J)  = OUXH(J) * (PT(I+1,J) + PT(I-2,J))            
         DpDrX(I,J) = OUXH(J) * (pDrho(I+1,J) - pDrho(I-1,J))                         !(*5)			
      ENDDO
      call period( PXW  (1,J) )                 
      call period( PXW2 (1,J) )
      call period( DpDrX(1,J) )
   end do
   if (myid_y.eq.0) then
      DO I = 1 ,NX
         PYP(I,NY)   = ZERO
         PYM(I,NY)   = ZERO
         rhos(i,NY)  = ZERO
         pDrho(i,NY) = ZERO
         PXW(I,NY)   = ZERO
         PXW2(I,NY)  = ZERO         
         DpDrX(I,NY) = ZERO		
      ENDDO
   else if(myid_y.eq.(npr_y-1)) then
      DO I = 1 ,NX
         PYP(I,1)   = ZERO
         PYM(I,1)   = ZERO
         rhos(i,1)  = ZERO
         pDrho(i,1) = ZERO
         PXW(I,1)   = ZERO
         PXW2(I,1)  = ZERO
         DpDrX(I,1) = ZERO
      ENDDO
   endif

#if (defined SPMD)

   call mp_send3d( commyz, dest, src, NX,NY,1,                     &
                   1, NX, beglatdynex, endlatdynex,1,1,       &
                   1, NX, beglatdyn, beglatdyn,1,1,  pDrho )
   do kk=beglev,endlev
      sbuf1(1:NX,5+2*num1+kk-beglev)=ut(1:NX,kk,beglatdyn)
      sbuf1(1:NX,5+3*num1+kk-beglev)=u(1:NX,kk,beglatdyn)
   end do

   call mp_recv3d( commyz, src, NX,  NY,1,                            &
                   1, NX, beglatdynex, endlatdynex,1,1,       &
                   1, NX, endlatdynex, endlatdynex,1,1,pDrho )
   call mpi_move_right(sbuf1,rbuf1,NX*(4*num1+4))

   if(myid_y.ne.(npr_y-1)) then
      do kk=beglev,endlev
         vt(:,kk,beglatdynex)=rbuf(:,5+kk-beglev)
         v(:,kk,beglatdynex)=rbuf(:,5+num1+kk-beglev)
      end do
   end if

   call mpi_w_a(rreq)

   call mp_send3d( commyz, dest1, src1, NX, NY,1,             &
                   1, NX, beglatdynex, endlatdynex,1,1,       &
                   1, NX, endlatdyn, endlatdyn, 1, 1, pDrho )
   if(myid_y.ne.0) then
      ksa(1:NX,endlatdynex)=rbuf1(1:NX,1)
      p(1:NX,endlatdynex)=rbuf1(1:NX,2)
      Pstar2(1:NX,endlatdynex)=rbuf1(1:NX,3)
      Pstar1(1:NX,endlatdynex)=rbuf1(1:NX,4)
   endif
!
   call mp_recv3d( commyz, src1, NX,  NY,1,                   &
                   1, NX, beglatdynex, endlatdynex,1,1,       &
                   1, NX, beglatdynex, beglatdynex,1,1, pDrho )
#endif

   DpDrY=0.
   DO J = loc_JB, loc_JE
      DO I = 1, NX
         DpDrY(i,j) = RDLATH*(pDrho(i,j+1) - pDrho(i,j-1))
      END DO
   END DO

!     COMPUTE    DTT/DT

#if (defined SPMD)

   call mp_send3d_2( commyz, dest, src, NX,NY,1,                     &
                     1, NX, beglatdynex,endlatdynex,1,1,       &
                     1, NX, beglatdyn, beglatdyn,1,1,  rhos,dpdry )
   if (myid_y.ne.0) then
      do kk=beglev,endlev
         tt(1:NX,kk,endlatdynex)=rbuf1(1:NX,5+kk-beglev)
         deltap(1:NX,kk,endlatdynex)=rbuf1(1:NX,5+num1+kk-beglev)
      end do
   endif

   call mp_recv3d_2( commyz, src, NX,  NY,1,                    &
                     1, NX, beglatdynex, endlatdynex,1,1,       &
                     1, NX, endlatdynex, endlatdynex,1,1,rhos,dpdry )
 
   call mp_send3d_2( commyz, dest1, src1, NX, NY,1,             &
                     1, NX, beglatdynex, endlatdynex,1,1,       &
                     1, NX, endlatdyn, endlatdyn, 1, 1, rhos,dpdry )
   if (myid_y.ne.(npr_y-1)) then
      ksa(:,beglatdynex)=rbuf(:,2)
      p(:,beglatdynex)=rbuf(:,3)
      Pstar2(:,beglatdynex)=rbuf(:,4)
   endif

   call mp_recv3d_2( commyz, src1, NX,  NY,1,                   &
                     1, NX, beglatdynex, endlatdynex,1,1,       &
                     1, NX, beglatdynex, beglatdynex,1,1, rhos,dpdry )

   if (npr_z.gt.1) then
      DIV=0.
      DO K=1,NL
          D2(K)=0.
          D3(K)=0.
      ENDDO
   endif

#endif

   DO J = loc_JB,loc_JE
! -------------------------------- COMPUTE DIVERGENCES & THEIR SUMS --------------------------
      DO K = beglev ,endlev
         DO I = IB, IE
            DIV(I,K,J) = alpha(J)*(PXW(I+1,J) *UT(I+1,K,J) - PXW(I,J)*UT(I  ,K,J))   &
                       + betaw(J)*(PXW2(I+2,J)*UT(I+2,K,J) - PXW2(I-1,J)*UT(I-1,K,J))/3.0 & !(4.21)
                       + (PYP(I,J)*VT(I,K,J)-PYM(I,J)*VT(I,K,J-1))                          !(4.22)
         ENDDO
      ENDDO
   end do !jjr j
   if (myid_y.eq.(npr_y-1)) then   ! at north polar 
      DO K = beglev ,endlev
         DO I = IB,IE
            PV1(I) = (PT(IB,1)+PT(I,JB))*VT(I ,K,1)    !PT(IB,1)=PT(I,1)       !(4.20)-(2)
         ENDDO
         D1K  = ZERO
         DO I = IB,IE
            D1K = D1K  + PV1(I)
         ENDDO
         D2(K)  = D1K  *  RDLATI     ! RDLATI = 2 / (IM*a*d(sita))
      ENDDO
      if(npr_z.gt.1) then
#if (defined SPMD)
         call parcollective(comm_z,sumop,NL,D2)
#endif
      end if
   endif
   if (myid_y.eq.0) then   ! at north polar
      DO K = beglev ,endlev
         DO I = IB,IE
            PV1(I) = (PT(IB,NY)+PT(I,JE))*VT(I,K,JE)                               !(4.20)-(2)
         ENDDO
         D1K = ZERO
         DO I = IB,IE
            D1K = D1K + PV1(I)
         ENDDO
         D3(K)  = -D1K  *  RDLATI
      ENDDO
      if (npr_z.gt.1) then
#if (defined SPMD)
         call parcollective(comm_z,sumop,NL,D3)
#endif
      end if
   ENDIF
      
#if (defined SPMD)
   if(myid_z.eq.0) then
#endif
      do k=beglev,endlev
         do j=beglatdyn,endlatdyn
            do i=IB,IE
               TT1(i,j,1,k)=DIV(i,k,j)
               TT1(i,j,2,k)=TT(i,k,j)
            end do
         enddo
      end do
#if (defined SPMD)
   endif
#endif

   if(npr_z.gt.1) then
#if (defined SPMD)
      do k=beglev,endlev
         do j=beglatdyn,endlatdyn
            do i=IB,IE
               TT0(i,j,1,k)=DIV(i,k,j)
               TT0(i,j,2,k)=TT(i,k,j)
            end do
         enddo
      end do
!jjr  TT0(IB:IE,beglatdyn:endlatdyn,beglev:endlev)=DIV(IB:IE,beglev:endlev,beglatdyn:endlatdyn)
      disps(1)=0
      do i=1,npr_z
         recs(i)=zdist(i)*(IE-IB+1)*(endlatdyn-beglatdyn+1)*2
         if(i.ne.npr_z) disps(i+1)=disps(i)+recs(i)
      end do

      num1=zdist(myid_z+1)*(endlatdyn-beglatdyn+1)*(IE-IB+1)*2
      call mpi_allgatherv(TT0,num1,mpir8,TT1,recs,disps,mpir8,comm_z,ierr)
#endif
   endif

     DO J = loc_JB,loc_JE 
         DO I = IB, IE
! -------------------------------- COMPUTE DIVERGENCES & THEIR SUMS --------------------------
            DO K = 1 ,NL
               D1(K)  = TT1(I,j,1,k)
               DWK(K) = D1(K) * DSIG(K)
            ENDDO
            DPS1 = -DWK(1)                 ! for Psa'
            DO K = 2 ,NL
               DPS1 = DPS1 - DWK(K)        ! sums for K
            ENDDO
!
! ******************************  compute Dsa *******************************
            DsaX       = OUXH(j)*(rhos(i+1,j)*ksa(i+1,j)*DpDrX(i+1,j)             &
                       -          rhos(i-1,j)*ksa(i-1,j)*DpDrX(i-1,j))                  !(4.23)
            DsaY       = RUPPH(j)*(rhos(i,j+1)*ksa(i,j+1)*DpDrY(i,j+1))           &
                       - RUMMH(j)*(rhos(i,j-1)*ksa(i,j-1)*DpDrY(i,j-1))                 !(4.23)
            Dsa(i,j)   = (DsaX + DsaY) 	   
			  
!-------------------------------- COMPUTE DP'sa/DT & D(SIGMA)/DT -------------------------------
            DPsa(i,j)  = P00*(DPS1) + kappa*Dsa(i,j)       
!!debug  if(beglev.eq.7.and.j.eq.180) print*,'test DPSA',DSAY,rhos(i,j+1),rhos(i,j-1),&
!! i,iam,testnum,ksa(i,j+1),ksa(i,j-1),Dpdry(i,j+1),DpDry(i,j-1)
!---------------------------------------------------------
            IF ( abs(DPS1).GT.0.001 ) THEN
               PRINT*, 'At sub. tend_lin, IN THE INNER AREA: DP=',DPS1, 'i=',I, 'j=', J
               PRINT*, 'DIV(i,26,j)=', D1(26), 'DIVX=', DIVX, 'DIVY=', DIVY
            ENDIF
!---------------------------------------------------------
            WK5 = ONE / PT(I,J)
            DPS1 = DPsa(i,j)
            DO K = 1 ,NM      ! NM = NL - 1
               DWK(K)  = WK5*((DPS1-kappa*Dsa(i,j))/P00 + D1(K))*DSIG(K)                 !(4.31)
            ENDDO
            WSK1  = ZERO
            DO K = 2 ,NL
               WSK1     = WSK1 - DWK(K-1)
               WSK(k)=WSK1
            end do
             WSK(NZ)=0.
             WSK(1)=0.
            do k=beglev,endlev+1
               WST(I,K,J)    = WSK(k)      
               WS(I,K,J)     = WST(I,K,J)/PT(I,J)  !WS = d(¦Ò)/dt; WS(1)=WS(NZ)=0
            ENDDO

         ENDDO       ! end i = IB,IE
!!!!

         call period( DPsa   (1,J) )
         DO K  = beglev ,endlev+1
            call period( WS (1,K,J) )
            call period( WST(1,K,J) )
         ENDDO
!------------------------ Compute P-surface vertical velocity (WPV) ---------------------------
         DO K = beglev ,endlev
            DO I = IB,IE
               DIVX = OUXQ(J) * (U(I,K,J) + U(I+1,K,J)) * (P(I+1,J) - P(I-1,J))     
               DIVY = RDLATQ  * (V(I,K,J) + V(I,K,J-1)) * (P(I,J+1) - P(I,J-1))
               WPV(I,K,J) = HALF * (WS(I,K+1,J) + WS(I,K,J)) * P(I,J)               &
                          + SIGL(K) * (DPsa(i,j) + DIVX + DIVY)                     ! P124 (8)
            ENDDO
            call period( WPV(1,K,J) )               		 
         ENDDO
!
!------------------ COMPUTE DTT/DT FOR J [JB,JE]  &  I [IB,IE]  &  K [1,NL] ---------------------
         DO K = beglev ,endlev
            DO I = IB,IE
               DXP0 = Pstar2(I  ,J) - Pstar2(I-1,J)
               DXP1 = Pstar2(I+1,J) - Pstar2(I,J)
               DXP2 = Pstar2(I  ,J) - Pstar2(I-3,J)
               DXP3 = Pstar2(I+3,J) - Pstar2(I,J)
               DYP0 = Pstar2(I,J+1) - Pstar2(I,J)
               DYP2 = Pstar2(I,  J) - Pstar2(I,J-1)
!
               OO1  = HALF*(WST(I,K+1,J) + WST(I,K,J)) / SIGL(K)                       &
                    - (DIV(I,K,J)/PT(I,J) + TWO*SGH(K)*(WST(I,K+1,J) - WST(I,K,J)))         !(4.19)
               OO2X = OUXH(J)/Pstar1(I,J)*(alpha(J)*(DXP1*UT(I+1,K,J)+DXP0*UT(I,K,J))  &
                                      +betaw(J)/3.0*(DXP3*UT(I+2,K,J)+DXP2*UT(I-1,K,J)))  !(4.18)
               OO2Y = (RUPH(J)*DYP0*VT(I,K,J) + RUMH(J)*DYP2*VT(I,K,J-1)) / Pstar1(I,J)   !(4.17)  
               TOO0 = (1. - deltap(I,K,J)) * ( b0 * (1. + deltac(I,K,J))               & 
                    + Kstd*CAPA*TT(I,K,J)/PT(I,J) ) * (b1*OO1 + b22*OO2X + b21*OO2Y)  
               DT(I,K,J)  = TOO0                       ! (3.3)
!!debug if(k.eq.7.and.j.eq.179) print*,'test DT00',DT(i,7,179),VT(i,7,178),WST(i,7,179),&
!! DIV(i,k,j),i,iam,testnum

            ENDDO
            call period( DT(1,K,J) )
         ENDDO
      enddo
!------------------------------- FOR THE POLAR POINTS -----------------------------------
      if(myid_y.eq.(npr_y-1)) then
         W11 = ZERO
! -------------------------COMPUTE DIVERGENCE AT NORTH POLE------------------------------
         DPS1 = ZERO
         DO K = 1 ,NL
            DPS1   = DPS1 - D2(K)*DSIG(K)                 ! for Psa'
         ENDDO
!
!ZRT---------------------------------------------------------
         IF (DPS1.GT.0.1) THEN
            PRINT*,'IN NORTH POLAR:    DP=',DPS1,I,J
            PRINT*,'DIV=',D2
         ENDIF
!ZRT---------------------------------------------------------
         DO I = IB,IE
            DsaY1(I) = rhos(I,2) * ksa(I,2) * DpDrY(I,2)                               !(4.23')
         ENDDO
         DSK = ZERO
         DO I = IB,IE
            DSK = DSK  + DsaY1(I)
         ENDDO
         DSK = DSK * RDLATI2
!----------------------COMPUTE DP/DT & D(SIGMA)/DT AT NORTH POLE------------------------
         DPsa1 = P00*(DPS1) + kappa*DSK
			
         DO K = 1 ,NM               ! NM = NL - 1
            DWK(K)    = ((DPsa1-kappa*DSK)/P00 + D2(K))*DSIG(K)/PT(IB,1)   !jjr
         ENDDO
         DO K = 2 ,NL
            W11     = W11 - DWK(K-1)  !W1(1) = ZERO
            W1(k)=W11
         ENDDO
         DO I = 1 ,NX
            DPsa(I,1) = DPsa1
         ENDDO
            W1(NZ)=0.
            W1(1)=0.

         DO K = beglev ,endlev+1
            DO I = 1 ,NX
               WST(I,K,1)  = W1(K)
               WS (I,K,1)  = WST(I,K,1)/PT(I,1)
            ENDDO
         ENDDO
         DO K = beglev ,endlev
            WPV1 (K) = HALF * (WS(IB,K+1,1) + WS(IB,K,1)) * P(IB,1)                  &
                     + SIGL(K) * DPsa(IB,1)          ! DIV(Pes) = 0
            DO I = 1 ,NX
               WPV(I,K,1) = WPV1 (K)
            ENDDO
         ENDDO

!---------------------------------COMPUTE DTT/DT AT NORTH POLE------------------------------------
         DO K = beglev ,endlev
!           COMPUTE HORIZONTAL & VERTICAL ADVECTION OF TT
            O2P    = ZERO
            DO I   = IB,IE
               O2P = O2P + (Pstar2(I,JB)-Pstar2(IB,1))*VT(I ,K,1) / Pstar1(IB,1)       !(4.17)-(2)
            ENDDO
            O2P    = RDLATI * O2P
            O1P    = HALF*(WST(IB,K+1,1) + WST(IB,K,1)) / SIGL(K)         &
                   - (D2(K)/PT(IB,1) + TWO*SGH(K)*(WST(IB,K+1,1) - WST(IB,K,1)))         !(4.19)
            TOOP = (1. - deltap(IB,K,1)) * ( b0 * (1. + deltac(IB,K,1))               & 
                 + Kstd*CAPA*TT(IB,K,1)/PT(IB,1) ) * (b1*O1P + b21*O2P)                                                         
            DTTP   = TOOP
            DO I = 1 ,NX
               DT(I,K,1) = DTTP
            ENDDO
         ENDDO
     endif !jjr
!  -------------------------------- FOR THE SOUTH POLAR POINTS --------------------------------------
       if(myid_y.eq.0) then
         W1(1) = ZERO
         W1(NZ)= ZERO
!  COMPUTE DIVERGENCE AT SOUTH POLE
         DPS1  = ZERO
         DO K = 1 ,NL
            DPS1   = DPS1 - D3(K)*DSIG(K)
         ENDDO

!ZRT---------------------------------------------------------
         IF (DPS1.GT.0.1) THEN
            PRINT*,'IN THE SOUTHEN POLAR:    DP=',DPS1,I,J
            PRINT*,'DIV=',D3
         ENDIF
!ZRT---------------------------------------------------------
         DO I = IB,IE
            DsaY1(I) = rhos(I,JE) * ksa(I,JE) * DpDrY(I,JE)                             !(4.23')
         ENDDO
         DSK = ZERO
         DO I = IB,IE
            DSK = DSK  + DsaY1(I)
         ENDDO
         DSK = -DSK * RDLATI2
!-------------------- COMPUTE DP/DT & D(SIGMA)/DT & WPV AT SOUTH POLE ------------------------
         DPsa1 = P00*(DPS1) + kappa*DSK
!
         DO K = 1 ,NM               ! NM = NL - 1
!            DWK(K)    = (kappa*(DPsa1-DSK)/P00 + D1(K))*DSIG(K)/PT(IB,NY)
            DWK(K)    = ((DPsa1-kappa*DSK)/P00 + D3(K))*DSIG(K)/PT(IB,NY)    !jjr
         ENDDO
         DO K = 2 ,NL
            W1(K)     = W1(K-1) - DWK(K-1)  !W1(1) = ZERO
         ENDDO
         DO I = 1 ,NX
            DPsa(I,NY)    = DPsa1
         ENDDO
         DO K = beglev ,endlev+1
            DO I = 1 ,NX
               WST(I,K,NY)  = W1(K)
               WS (I,K,NY)  = WST(I,K,NY)/PT(I,NY)
            ENDDO
         ENDDO
         DO K = beglev ,endlev
            WPV1 (K) = HALF * (WS(IB,K+1,NY) + WS(IB,K,NY)) * P(IB,NY)                 &
                     + SIGL(K) * DPsa(IB,NY)
            DO I = 1 ,NX
               WPV(I,K,NY) = WPV1 (K)
            ENDDO
         ENDDO
!
!-------------------------------- COMPUTE DTT/DT AT south POLES -------------------------------------
         DO K = beglev ,endlev
!           COMPUTE HORIZONTAL & VERTICAL ADVECTION OF TT
            O2P    = ZERO
            DO I   = IB,IE
               O2P = O2P + (Pstar2(IB,NY)-Pstar2(IB,JE))*VT(I,K,JE) / Pstar1(IB,NY)    !(4.17)-(3)
            ENDDO
            O2P    = RDLATI * O2P
            O1P    = HALF*(WST(IB,K+1,NY) + WST(IB,K,NY)) / SIGL(K)         &
                   - (D3(K)/PT(IB,NY) + TWO*SGH(K)*(WST(IB,K+1,NY) - WST(IB,K,NY)))      !(4.19)
            TOOP = (1. - deltap(IB,K,NY)) * ( b0 * (1. + deltac(IB,K,NY))               & 
                 + Kstd*CAPA*TT(IB,K,NY)/PT(IB,NY) ) * (b1*O1P + b21*O2P)                                                         
            DTTP   = TOOP

            DO I = 1 ,NX
               DT(I,K,NY) = DTTP
            ENDDO
         ENDDO
      ENDIF
!
      IF ( NADDSS.EQ.+1 ) THEN     !zhh 2007.4.23
!         ADD SOURCE & SINK
         do j=beglatdyn,endlatdyn
            DO K = beglev,endlev 
               DO I = 1 ,NX
                  DT(I,K,J) = DT(I,K,J) + ST(I,K,J)
!!debug if(k.eq.7.and.j.eq.179) print*,'test DT0',DT(i,7,179),ST(i,7,179),&
!! i,iam,testnum
               ENDDO
            ENDDO
          end do
      ENDIF
!!
   call t_startf('fft3')
   WW(:,beglatdyn:endlatdyn)=DPsa(:,beglatdyn:endlatdyn)   !zhh 2008.6.16
   CALL FILT2D( WW,0,1,IBCFFT )		   !zhh 2008.6.16
   DPsa(:,beglatdyn:endlatdyn)=WW(:,beglatdyn:endlatdyn)

   DO KK = beglev,endlev
      WW(:,beglatdyn:endlatdyn)=DT(:,KK,beglatdyn:endlatdyn)
      CALL FILT2D( WW,0,1,IBCFFT )   !OK
      DT(:,KK,beglatdyn:endlatdyn)=WW(:,beglatdyn:endlatdyn)
   ENDDO
   call t_stopf('fft3')
!
!     *******************************************************************
!     *******************************************************************
!     ****************         COMPUTE DU & DV        *******************
!     *******************************************************************
!     *******************************************************************
!
!******************* update PLY & deltac for next time step **********************

   DO J = beglatdyn,endlatdyn
      DO K = 1 ,NL
         DO I = 1 ,NX
            PLY(I,K,J) = P(I,J)*SIGL(K) + PMTOP
!jjr            deltap(I,K,J) = PMTOP / PLY (I,K,J) !jjr
         ENDDO
      ENDDO
      DO I = 1 ,NX
         PLY(I,NZ,J)  = P(I,J) + PMTOP
      ENDDO
      DO K = beglev ,endlev
         DO I = 1 ,NX
            WPK  = PLY(I,K,J) / DPALIB  ! DPALIB = 0.5
            KPK  = WPK + 1.0E-4
            WPJ  = WPK - KPK
            WPI  = ONE - WPJ
            KPI  = KPK + 1
            IF (KPI.GT.NA) THEN   ! NA = PEALIB/DPALIB = 2320
               PRINT*,'************************THERE IS AN ERROR IN CBB***************'     
               PRINT*,'PLY(I,K,J)=',PLY(I,K,J),'I,K,J=',I,K,J
               PRINT*,'P(I,J)=',P(I,J)
               PRINT*, 'stop: tend_lin'
               stop 'tend_lin'
            ENDIF
            CB (I,K,J)    = WPI*CBB(KPK) + WPJ*CBB(KPI)
            if (adiabatic) then
               deltac(I,K,J) = 0.0
            else
               deltac(I,K,J) = CB(I,K,J) * CB(I,K,J) / (b0 * b0) - 1.0
            end if
         ENDDO
      ENDDO
!********************************** Compute GHI *******************************************
      DO I = 1 ,NX
         GHI1(I,NZ,J) = H0B(I,J) * (PLY(I,NZ,J) - PSB(I,J)) + kappa2 * GZsm(I,J)    !P113 (4)    
!          H0B £½ Rd * TMSA(PSB)/PSB,   TMSA: temperature of standard atmosphere 
      ENDDO
      DO K = 1 ,NL
         DO I = IB ,IE			                                                                       
            DGH(I,K,J) = b0 * TT1(I,J,2,K) * (1.-deltap(I,K,J)) / (PT(I,J)*SIGL(K))*DSIG(K)
			                                                                       !  (3.5)
         ENDDO
      ENDDO
      DO K = NL,1 ,-1
         DO I = IB ,IE
            GHI1(I,K,J) = GHI1(I,K+1,J) + DGH(I,K,J)
         ENDDO
      ENDDO
      DO K=beglev,endlev+1
         DO I=IB,IE
           GHI(I,K,J)=GHI1(I,K,J)
         ENDDO
      call period(GHI(1,k,j))
     ENDDO

   ENDDO

!
#if (defined SPMD)
   call mp_send3d( commyz, dest, src, NX,NZ, NY,                     &
                   1, NX,beglev,endlev+1, beglatdynex, endlatdynex,       &
                   1, NX, beglev,endlev+1,beglatdyn, beglatdyn,  GHI )
   if(myid_y.ne.0) then
      num1=endlev-beglev+1
      do kk=beglev,endlev
         ut(1:NX,kk,endlatdynex)=rbuf1(1:NX,5+2*num1+kk-beglev)
         u(1:NX,kk,endlatdynex)=rbuf1(1:NX,5+3*num1+kk-beglev)
      end do
   endif

   call mp_recv3d( commyz, src, NX,NZ,  NY,                             &
                   1, NX,beglev,endlev+1, beglatdynex, endlatdynex,       &
                   1, NX,beglev,endlev+1, endlatdynex, endlatdynex, GHI )
#endif
!
   DO J  = beglatdyn, endlatdyn
! --------------------- COMPUTE PY AT THE NORTH POLAR (J=1) -------------------------------
!
      IF (J.EQ.1) THEN
         DO K = beglev ,endlev
            DO I = IB,IE
               DYP0  = Pstar2(I,JB) - Pstar2(IB,1)
               OYAYP = RDLATQ*(PT(IB,1)+PT(I,JB))                       
               OPCT0 = b0 * TT(IB,K,1) * (1. - deltap(IB,K,1)) / Pstar1(IB,1)
               OPCT3 = b0 * TT(I,K,JB) * (1. - deltap(I,K,JB)) / Pstar1(I,JB)
!
               PY1   = OYAYP * ( (GHI(I,K+1,JB) + GHI(I,K,JB))    &       
                     -           (GHI(IB,K+1,1) + GHI(IB,K,1)) )          ! (4.4)
               PY2   = RDLATH * ( OPCT0 + OPCT3 ) * DYP0                  ! (4.5)
               PY(I,K,1) = b1*PY1 + b21*PY2
               PX(I,K,1) = ZERO
            ENDDO
            call period( PX(1,K,J) )
            call period( PY(1,K,J) )               
         ENDDO
      ELSE IF(J.LE.JE) THEN
!
!     COMPUTE PX & PY FROM J=JB TO J=JE
!
         DO K = beglev ,endlev
            DO I = IB,IE
               DXP0  = Pstar2(I  ,J) - Pstar2(I-1,J)
               DXP1  = Pstar2(I+1,J) - Pstar2(I-2,J)
               OUXAXP = OUXQ(J) * ( PT(I  ,J) + PT(I-1,J) )
               OUXAXP2= OUXQ(J) * ( PT(I+1,J) + PT(I-2,J) )
               DYP0  = Pstar2(I,J+1) - Pstar2(I,J)
               OYAYP = RDLATQ  * ( PT(I,J) + PT(I,J+1) )
!
               OPCT0 = b0 * TT(I  ,K,J) * (1. - deltap(I  ,K,J)) / Pstar1(I  ,J)
               OPCT2 = b0 * TT(I-1,K,J) * (1. - deltap(I-1,K,J)) / Pstar1(I-1,J)
               OPCT3 = b0 * TT(I,K,J+1) * (1. - deltap(I,K,J+1)) / Pstar1(I,J+1)
               OPCT4 = b0 * TT(I+1,K,J) * (1. - deltap(I+1,K,J)) / Pstar1(I+1,J)
               OPCT5 = b0 * TT(I-2,K,J) * (1. - deltap(I-2,K,J)) / Pstar1(I-2,J)
               GH0   = GHI(I  ,K+1,J) + GHI(I  ,K,J)
               GH1   = GHI(I-1,K+1,J) + GHI(I-1,K,J)
               GH2   = GHI(I+1,K+1,J) + GHI(I+1,K,J)
               GH3   = GHI(I-2,K+1,J) + GHI(I-2,K,J)
               GH4   = GHI(I,K+1,J+1) + GHI(I,K,J+1)

               PX1   = alpha(J) * OUXAXP  * ( GH0 - GH1 )                      &
                     + betaw(J) * OUXAXP2 * ( GH2 - GH3 ) / 3.0                     ! (4.10)
               PX2   = alpha(J) * OUXH(J) * ( OPCT0 + OPCT2 ) * DXP0           &
                     + betaw(J) * OUXH(J) * ( OPCT4 + OPCT5 ) * DXP1 / 3.0          ! (4.11)
               PY1   = OYAYP  * ( GH4 - GH0)                                        ! (4.4 )
               PY2   = RDLATH * ( OPCT0 + OPCT3 ) * DYP0                            ! (4.5 )
               PX(I,K,J)  = b1*PX1 + b22*PX2
               PY(I,K,J)  = b1*PY1 + b21*PY2
            ENDDO
            call period( PX(1,K,J) )
            call period( PY(1,K,J) )
         ENDDO
      ELSE        ! AT THE SOUTH POLAR (J=NY)
         DO K = beglev ,endlev
            DO I = 1 ,NX
               PX(I,K,J)  = ZERO
               PY(I,K,J)  = ZERO
            ENDDO
         ENDDO
      ENDIF
!
!         COMPUTE DV/DT AT J=1
!
      IF (J.EQ.1) THEN    !AT THE NORTH POLAR  
         DO K = beglev ,endlev
            DO I = IB,IE
               fstar0 = FF(JB)+CUR(JB)  * U(I  ,K,JB)
               fstar1 = FF(JB)+CUR(JB)  * U(I+1,K,JB)
               FSU = FOURTH * ( fstar0 * UT(I,K,JB) + fstar1 * UT(I+1,K,JB) )
               DV(I,K,1) = - PY(I,K,1) + c10*FSU       ! (3.1)	
               DU(I,K,1) = ZERO
            ENDDO
            call period( DU(1,K,J) )
            call period( DV(1,K,J) )
         ENDDO
!
      ELSE IF(J.LE.JE) THEN
! ----------------- COMPUTE DU/DT & DV/DT FOR J [JB,JE] & I [IB,IE] & K [1,NL] ------------------
!         
         DO K = beglev ,endlev
            DO I = IB,IE
               FS0  = FF(J  ) + CUR(J  ) * U(I  ,K,J)                             ! P101 (18)
               FS1  = FF(J  ) + CUR(J  ) * U(I+1,K,J)
               FS2  = FF(J+1) + CUR(J+1) * U(I  ,K,J+1)
               FS3  = FF(J+1) + CUR(J+1) * U(I+1,K,J+1)
! 
               FSV  = FS0  * ( RUPD(J) * (VT(I,K,  J) + VT(I-1,K,  J))                      &
                    +          RUMD(J) * (VT(I,K,J-1) + VT(I-1,K,J-1)) )            ! (4.12)
               FSU  = FOURTH * ( FS0 * UT(I,K,J)   + FS1 * UT(I+1,K,J)                      &
                    +            FS2 * UT(I,K,J+1) + FS3 * UT(I+1,K,J+1) )          ! (4.6)
               DU(I,K,J) = - PX(I,K,J) - c10*FSV      ! (3.2)    
               DV(I,K,J) = - PY(I,K,J) + c10*FSU      ! (3.1)   
            ENDDO
            call period( DU(1,K,J) )
            call period( DV(1,K,J) )
         ENDDO
!
      ELSE   !  J = NY  AT THE SOUTH POLAR
         DO K = beglev ,endlev
            DO I = 1 ,NX
               DU(I,K,J)= ZERO
               DV(I,K,J)= ZERO
            ENDDO
         ENDDO
      ENDIF
!
!!      print*, 'NADDSS =', NADDSS
!  ADD SOURCE & SINK if necessary
      IF ( NADDSS.EQ.+1 ) THEN
         DO K = beglev ,endlev
            DO I = 1 ,NX
!=============================== zhh =====================================
               IF (abs(DU(I,K,J)) > 5E-1 .OR. abs(DV(I,K,J)) > 5E-1 .or.      &
                   abs(SU(I,K,J)) > 5E-1 .OR. abs(SV(I,K,J)) > 5E-1) THEN
                  print*, 'DU(',I,K,J,') =', DU(I,K,J)
                  print*, 'DV(',I,K,J,') =', DV(I,K,J)
                  print*, 'SU(',I,K,J,') =', SU(I,K,J)
                  print*, 'SV(',I,K,J,') =', SV(I,K,J)
                  print*, 'tend_lin--2'
 !!                 stop
               ENDIF
!============================ 2007.8.5 ====================================
               DU(I,K,J) = DU(I,K,J) + SU(I,K,J)
               DV(I,K,J) = DV(I,K,J) + SV(I,K,J)
            ENDDO
         ENDDO
      ENDIF
   ENDDO        !end J = 1,NY
!
!     HIGH-MID LAT FILTER  FOR DU/DT & DV/DT

   DO KK = beglev,endlev
      WW(:,beglatdyn:endlatdyn)=DU(:,KK,beglatdyn:endlatdyn)
      CALL FILT2D( WW,1,1,IBCFFT )	
      DU(:,kk,beglatdyn:endlatdyn)=WW(:,beglatdyn:endlatdyn)
!         
      WW(:,beglatdyn:endlatdyn)=DV(:,KK,beglatdyn:endlatdyn)
      CALL FILT2D( WW,1,2,IBCFFT )	
      DV(:,kk,beglatdyn:endlatdyn)=WW(:,beglatdyn:endlatdyn)
   end do
   call t_stopf('tend_lin')
!
! deallocate arrays
   deallocate(WW)
   deallocate(DGH)
   deallocate(GHI1)
   deallocate(TT1)
   deallocate(TT0)
   deallocate(PXW)
   deallocate(PXW2)
   deallocate(PYP)
   deallocate(PYM)
   deallocate(rhos)
   deallocate(DpDrX)
   deallocate(DpDrY)
   deallocate(DIV)
   deallocate(Dsa)
   deallocate(PX)
   deallocate(PY)
   deallocate(sbuf)
   deallocate(sbuf1)
   deallocate(rbuf)
   deallocate(rbuf1)
!
   RETURN
END
