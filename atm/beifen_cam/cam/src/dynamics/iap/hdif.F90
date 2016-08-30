module hdif
!-----------------------------------------------------------------------------------------------
! Purpose:  Compute horizontal diffusion of U, V, T & Q; 
!           non-linear diffusion of Smagorinsky (1963)
! Reference: [1]. W.M.Washington & D.L.Wiliamson, 1977: A description of the NCAR global 
!                 circulation models, General circulation models of the atmosphere. 
!                 (A78-10662 01-47) New York, Academic Press, Inc , pp. 111-172.
!            [2]. X.Z.Liang, 1986: The Design of IAP GCM and  the Simulation of Climate
!                 and Its Interseasonal Variability.  Ph.D Thesis, 250pp
! Original version : HDIFUS.f (IAP 9L)
! Recoded to module & revised by: ZhangHe
! Completed : 2006.4.30 
! Update : 2007.4.16, ZhangHe, added the dumb parameter 'deltat' to sub. HDIFUS
!          2007.4.25, ZhangHe, added the dumb parameter 'PLY, U, V, T, Q' to sub. HDIFUS 
!          2008.4.21, ZhangHe, removed diffusion of q
!          2008.5, Wujianping, parallel version
!          2008.6.10, ZhangHe, available for both serial & parallel
!          2012 OCT, Jiang Jinrong, 2D parellel
! Modified: Zhang He, 2013-03-08, new mp_send3d, dyn_state was added
!-----------------------------------------------------------------------------------------------
!
!                       i - 1/2        i          i + 1/2
!
!                         |            |            |
!                         |            |            |
!           j - 1/2  ---- Ds -------- D,v --------- Ds ----  j' - 1
!                         |            |            |
!                         |            |            |
!                         |            |            |
!                         |            |            |
!                         |            |            |
!             j      ---- u ---------- DT --------- u -----    j
!                         |            |            |
!                         |            |            |
!                         |            |            |
!                         |            |            |
!                         |            |            |
!           j + 1/2  ---- Ds -------- D,v --------- Ds ----    j'
!                         |            |            |
!                         |            |            |
!
!                         i'           i          i'+ 1
!
!       Fig. H.1 Distribution of variables related to horizontal diffusion on C-grid    
!-----------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,  only: NX, NY, NL, NZ, IB, IE, JB, JE, period
   use mathconst, only: ZERO, HALF, ONE, FOUR
   use pmgrid,    only: beglatdyn, endlatdyn, loc_JB, loc_JE, beglatdynex, endlatdynex,beglev,endlev,endlevp
   use dynamics_vars,   only: T_FVDYCORE_STATE         !zhh 2013-03-08
   use dyn_internal_state,   only : get_dyn_state      !zhh 2013-03-08
#if (defined SPMD)
   use mod_comm, only: mp_send3d, mp_recv3d
   use pmgrid,   only: npr_y
   use spmd_utils, only: masterproc, iam     !zhh 2013-03-08
#endif
 
   implicit none
   
   real(r8) :: FRDT(NY,3) 
   real(r8) :: FRDS(NY,3) 
   real(r8) :: FRDU(NY,3) 
   real(r8) :: FRDV(NY,3) 
   real(r8) :: FRDP(NY,3) 
   real(r8), parameter :: Kdif = 0.1

!================================================================================================
CONTAINS
!================================================================================================
 
!================================================================================================
!!   SUBROUTINE HDIFUS( deltat, PLY, U, V, T, Q )
   SUBROUTINE HDIFUS( deltat, PLY, U, V, T )
!-----------------------------------------------------------------------------------------------
!  Compute horizontal diffusion of U, V, T & Q
!-----------------------------------------------------------------------------------------------
      use flexib, only: DTHDFS
	  use physconst, only: CP, RD
      use Dyn_const, only: DXVPN, DXVPS
      use stdatm, only: TB  

	  
	  implicit none
!----------------------------------------Arguments-----------------------------------------------
      real(r8), intent(in) :: deltat            ! model timestep
      real(r8), intent(in) :: PLY(NX,NZ,beglatdynex:endlatdynex)     ! pressure at model sigma layer (unit: hPa)
      real(r8), intent(inout) :: U(NX,beglev:endlev,beglatdynex:endlatdynex)    ! latitudinal wind
      real(r8), intent(inout) :: V(NX,beglev:endlev,beglatdynex:endlatdynex)    ! meridional wind
      real(r8), intent(inout) :: T(NX,beglev:endlev,beglatdynex:endlatdynex)    ! air temperature
!!      real(r8), intent(inout) :: Q(NX,NL,beglatdynex:endlatdynex)    ! specific humidity
!----------------------------------Local workspace----------------------------------------------      
      real(r8) :: DU, DV, DdT, DQ   ! tendency of U, V, T & Q due to horizontal diffusion
      real(r8) :: Kd2
      real(r8) :: D(NX,beglatdynex:endlatdynex), DT(NX,beglatdynex:endlatdynex), DS(NX,beglatdynex:endlatdynex), &
                 DA(NX,beglatdynex:endlatdynex), DB(NX,beglatdynex:endlatdynex)
	  real(r8) :: ROT(NX,beglatdynex:endlatdynex) 
	  real(r8) :: VR(NX,beglatdynex:endlatdynex)
!!	  real(r8) :: QK(NX,NY), TK(NX,NY), VK(NX,NY), UK(NX,NY), TW(NX,NY)
      real(r8) :: TK(NX,beglatdynex:endlatdynex), VK(NX,beglatdynex:endlatdynex),&
             UK(NX,beglatdynex:endlatdynex), TW(NX,beglatdyn:endlatdyn)   !!zhh
	  real(r8) :: RLNT(NX,beglatdynex:endlatdynex)      ! RLNT(I,J) = 4 * rho(i', j')
	  real(r8) :: RDLN(NX,beglatdynex:endlatdynex)      ! RDLN(I,J) = 2 * rho(i', j ) * D(i', j )
	  real(r8) :: RDLT(NX,beglatdynex:endlatdynex)      ! RDLT(I,J) = 2 * rho(i , j') * D(i , j')
	  real(r8) :: FRDTN, FRDTS, FRDSI, FRDSJ
	  real(r8) :: TI, RIJ, RI, VRI, TDI
	  real(r8) :: FT1, FT2, FT3, FS1, FS2, FS3 
	  real(r8) :: U0, V0
	  real(r8) :: DTN, DTS, DVN
	  real(r8) :: DTJ, DSI, DIJ
	  real(r8) :: R0, R1
	  real(r8) :: D0, D1
	  real(r8) :: RT, RN
	  real(r8) :: DIJ1, DIJ2
	  real(r8) :: FU1, FU2, FU3
	  real(r8) :: FV1, FV2, FV3
	  real(r8) :: FA1, FA2, FA3
	  real(r8) :: VR0
	  real(r8) :: DS0, DT0, DA0, DB0
	  real(r8) :: ROT0
	  real(r8) :: RLNT0, RTA0, RSB0
	  real(r8) :: VRU, VRV
	  real(r8) :: RTAI, RSBI, RTAJ, RSBJ
	  real(r8) :: FB1, FB2, FB3, FB4
	  real(r8) :: TD0, Q0
!!	  real(r8) :: VRN, VRS, TDN, TDS, QKN, QKS, DQN, DQS, FAN, FAS, EAN, EAS, T0N, T0S, Q0N, Q0S
      real(r8) :: VRN, VRS, TDN, TDS, FAN, FAS, EAN, EAS, T0N, T0S
	  real(r8) :: FVI, FVJ
	  real(r8) :: DXVPNJ, DXVPSJ
      integer  :: I, J, K, J1, JJ, I1, II    ! loop index
      integer :: src,dest
      type (T_FVDYCORE_STATE), pointer :: dyn_state   !zhh 2013-03-08
      integer  :: commyz     !zhh 2013-03-08

!-----------------------------------------------------------------------------------------------
!
      dyn_state => get_dyn_state()
      commyz = dyn_state%grid%commyz

      DTHDFS = deltat / 2
      FRDTN = FRDT(1 ,1)  
      FRDTS = FRDT(NY,1)
      FRDSI = FRDS(1 ,1)
      FRDSJ = FRDS(1 ,2)
	  Kd2   = Kdif * Kdif
      
	  DO K = beglev, endlev
!
!     CALCULATE DENSITY & SET 2-D FIELDS TO BE DIFFUSED
         DO J = beglatdyn, endlatdyn
            DO I = 1, NX
               UK(I,J) = U(I,K,J)
               VK(I,J) = V(I,K,J)
!!               QK(I,J) = Q(I,K,J)
               TW(I,J) = T(I,K,J)
            END DO
         END DO
!      
!!	     DO J = loc_JB, loc_JE
         DO J = beglatdyn, endlatdyn
            if ( j >= JB .and. j <= JE ) then 
               DO I = IB, IE
                  TI       = TW(I,J)
                  RIJ      = PLY(I,K,J) / (TI*RD)   ! rho = p/(R*T)
                  ROT(I,J) = RIJ  
                  VR (I,J) = ONE / RIJ              ! VR = 1/rho
                  TK (I,J) = TI - TB(I,K,J)  
               END DO
            else   ! at north and south pole         
               TI  = TW (IB,J)
               RI  = PLY(IB,K,J) / (TI*RD)
               VRI = ONE / RI
               TDI = TI - TB(IB,K,J)
               DO I = IB, IE
                  ROT(I,J) = RI
                  VR (I,J) = VRI
                  TK (I,J) = TDI
               END DO
            end if
!--------------------- Set the spherical cyclicity condition ------------------
            call period( ROT(1,J) )
            call period( VR (1,J) )     
            call period( TK (1,J) )
         END DO
!------------------------------------------------------------------------------		 	 
!
#if (defined SPMD)
! Latitude distribution for plat=16 on 3 processors.
!  procId:      iam=0     |      iam=1     |      iam=2              !wjp 2007.05
!  CAM:    01 02 03 04 05 | 06 07 08 09 10 | 11 12 13 14 15 16       !wjp 2007.05
!  IAP:    16 15 14 13 12 | 11 10 09 08 07 | 06 05 04 03 02 01       !wjp 2007.05
!                left     |     current    |      right              !wjp 2007.05
!JJR right
      src = iam-1
      dest  = iam+1
      if ( mod(iam,npr_y) == 0 ) src = -1
      if ( mod(iam+1,npr_y) == 0 ) dest = -1
      call mp_send3d( commyz, dest, src, NX,  NY,1,              &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, beglatdyn, beglatdyn,1,1,  uk )
      call mp_recv3d( commyz, src, NX,  NY,1,                    &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, endlatdynex, endlatdynex,1,1,uk )
      call mp_send3d( commyz, dest, src, NX,  NY,1,              &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, beglatdyn, beglatdyn,1,1,  tk )
      call mp_recv3d( commyz, src, NX,  NY,1,                    &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, endlatdynex, endlatdynex,1,1,tk )
      call mp_send3d( commyz, dest, src, NX,  NY,1,              &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, beglatdyn, beglatdyn,1,1,  rot )
      call mp_recv3d( commyz, src, NX,  NY,1,                    &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, endlatdynex, endlatdynex,1,1,rot )
      call mp_send3d( commyz, dest, src, NX,  NY,1,              &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, beglatdyn, beglatdyn,1,1,  vr )
      call mp_recv3d( commyz, src, NX,  NY,1,                    &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, endlatdynex, endlatdynex,1,1,vr )

!JJR left
      src = iam+1
      dest  = iam-1
      if ( mod(iam,npr_y) == 0 ) dest = -1
      if ( mod(iam+1,npr_y) == 0 ) src = -1
      call mp_send3d( commyz, dest, src, NX, NY,1,               &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, endlatdyn, endlatdyn, 1, 1, vk )
      call mp_recv3d( commyz, src, NX,  NY,1,                    &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, beglatdynex, beglatdynex,1,1, vk )
      call mp_send3d( commyz, dest, src, NX, NY,1,               &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, endlatdyn, endlatdyn, 1, 1, tk )
      call mp_recv3d( commyz, src, NX,  NY,1,                    &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, beglatdynex, beglatdynex,1,1, tk )
#endif
!
!     CALCULATE DEFORMATION FIELDS DT, DS & D
!!         DO J = loc_JB, loc_JE
         DO J = beglatdyn, endlatdyn
            if ( j >= JB .and. j <= JE ) then 
               J1  = J - 1
               JJ  = J + 1
               FT1 = FRDT(J,1)
               FT2 = FRDT(J,2)
               FT3 = FRDT(J,3)
               FS1 = FRDS(J,1)
               FS2 = FRDS(J,2)
               FS3 = FRDS(J,3)
               DO I = IB, IE
                  U0      = UK(I,J)
                  V0      = VK(I,J)
                  DT(I,J) = FT1*(UK(I+1,J)-U0) - (FT2*V0 - FT3*VK(I,J1 ))   !(5.9)
                  DS(I,J) = FS1*(V0-VK(I-1,J)) + (FS2*UK(I,JJ) - FS3*U0)    !(5.10)
               END DO
            else if ( j == 1 ) then  ! at north pole         
	           DTN = ZERO
               DO I = IB, IE
                  DTN = DTN + VK(I,1 )  !at the north polar
               END DO
               DTN = DTN * FRDTN  !(5.9)-(2)
               DO I = IB, IE
                  DT(I,1 ) = DTN
                  DS(I,1 ) = FRDSI*(VK(I,1) - VK(I-1,1)) + FRDSJ*UK(I,JB) !at the north polar
               END DO
            else if ( j == NY ) then ! at south pole         
               DTS = ZERO
               DO I = IB, IE
                  DTS = DTS + VK(I,JE)  !at the south polar
               END DO
               DTS = DTS * FRDTS  !(5.9)-(3)
               DO I = IB, IE
                  DT(I,NY) = DTS
                  DS(I,NY) = ZERO
               END DO
            end if
!--------------------- Set the spherical cyclicity condition ------------------      
            call period( DT(1,J) )
            call period( DS(1,J) )     
!
         END DO
!
#if (defined SPMD)
! Latitude distribution for plat=16 on 3 processors.                 !wjp 2007.05
!  procId:      iam=0     |      iam=1     |      iam=2              !wjp 2007.05
!  CAM:    01 02 03 04 05 | 06 07 08 09 10 | 11 12 13 14 15 16       !wjp 2007.05
!  IAP:    16 15 14 13 12 | 11 10 09 08 07 | 06 05 04 03 02 01       !wjp 2007.05
!                left     |     current    |      right              !wjp 2007.05
!JJR 
     src = iam-1
      dest  = iam+1
      if ( mod(iam,npr_y) == 0 ) src = -1
      if ( mod(iam+1,npr_y) == 0 ) dest = -1
      call mp_send3d( commyz, dest, src, NX,  NY,1,              &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, beglatdyn, beglatdyn,1,1,  dt )
      call mp_recv3d( commyz, src, NX,  NY,1,                    &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, endlatdynex, endlatdynex,1,1,dt )
      src = iam+1
      dest  = iam-1
      if ( mod(iam,npr_y) == 0 ) dest = -1
      if ( mod(iam+1,npr_y) == 0 ) src = -1
      call mp_send3d( commyz, dest, src, NX, NY,1,               &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, endlatdyn, endlatdyn, 1, 1, ds )
      call mp_recv3d( commyz, src, NX,  NY,1,                    &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, beglatdynex, beglatdynex,1,1, ds )
#endif
!
! Deduce D from DT & DS
!
!!         DO J = beglatdyn, loc_JE
         DO J = beglatdyn, endlatdyn
            if ( j <= JE ) then 
               JJ = J + 1
               DO I = IB,IE
                  DTJ    = DT(I,JJ) + DT(I,J)
                  DSI    = DS(I,J ) + DS(I+1,J)   ! revised by zhh
                  DIJ    = DTJ*DTJ  + DSI*DSI
                  D(I,J) = HALF * SQRT( DIJ )     !(5.8)
               END DO
               call period( D(1,J) )
            else   ! at south pole
               DO I = 1 ,NX
                  D(I,NY) = ZERO
               ENDDO
            end if
         END DO
!
#if (defined SPMD)
! Latitude distribution for plat=16 on 3 processors.               !wjp 2007.05
!  procId:      iam=0     |      iam=1     |      iam=2            !wjp 2007.05
!  CAM:    01 02 03 04 05 | 06 07 08 09 10 | 11 12 13 14 15 16     !wjp 2007.05
!  IAP:    16 15 14 13 12 | 11 10 09 08 07 | 06 05 04 03 02 01     !wjp 2007.05
!                left     |     current    |      right            !wjp 2007.05
!JJR left     
      src = iam+1
      dest  = iam-1
      if ( mod(iam,npr_y) == 0 ) dest = -1
      if ( mod(iam+1,npr_y) == 0 ) src = -1
      call mp_send3d( commyz, dest, src, NX, NY,1,               &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, endlatdyn, endlatdyn, 1, 1, d )
      call mp_recv3d( commyz, src, NX,  NY,1,                    &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, beglatdynex, beglatdynex,1,1, d )
#endif
!
!         DO J = loc_JB, loc_JE
         DO J = beglatdyn, endlatdyn
            if ( j >= JB .and. j <= JE ) then 
               JJ = J + 1
               J1 = J - 1
               DO I = IB, IE
                  I1        = I - 1
                  R0        = ROT(I,J)
                  R1        = ROT(I,JJ)
                  D0        = D(I,J)                  ! D0 = D(i, j')
                  D1        = D(I,J1)                 ! D1 = D(i, j'-1)
                  RT        = R0 + R1                 ! RT = 2 * rho(i , j')
                  RN        = R0 + ROT(I1,J)          ! RN = 2 * rho(i', j )
                  RLNT(I,J) = RN + R1 + ROT(I1,JJ)    ! RLNT(I,J) = 4 * rho(i', j')
                  DIJ1      = HALF * (D0 + D(I1,J))                   ! DIJ1 = D(i', j') 
                  DIJ       = HALF * (D0 + D1)                        ! DIJ  = D(i , j )
                  DIJ2      = HALF * (DIJ1 + HALF*(D1 + D(I1,J1)) )   ! DIJ2 = D(i', j )
                  RDLN(I,J) = RN   *  DIJ2      ! RDLN(I,J) = 2 * rho(i', j ) * D(i', j  )
                  RDLT(I,J) = RT   *  D0        ! RDLT(I,J) = 2 * rho(i , j') * D(i , j' )
                  DA  (I,J) = DIJ
                  DB  (I,J) = DIJ1
               END DO
            else if ( j == 1 ) then ! at north pole         
               DO I = 1 ,NX
                  RDLN(I,J) = ZERO
                  DA(I,J)   = ZERO
               ENDDO
               R0 = ROT (IB,1)
               DO I = IB,IE     ! at the north polar
                  D0        = D(I,1)
                  DB  (I,1) = (D0 + D(I-1,1 )) * HALF
                  RDLT(I,1) = (R0 + ROT(I,JB)) * D0
                  RLNT(I,1) =  R0 + R0 + ROT(I-1,JB) + ROT(I,JB)
               END DO
            else if ( j == NY ) then ! at south pole         
               DO I = 1 ,NX
                  RDLN(I,J) = ZERO
                  DA(I,J)   = ZERO
                  RDLT(I,NY) = ZERO
                  RLNT(I,NY) = ZERO
                  DB(I,NY)   = ZERO
               ENDDO
            end if
!--------------------- Set the spherical cyclicity condition ------------------      
            call period( RLNT(1,J) )
            call period( RDLN(1,J) )     
            call period( RDLT(1,J) )
            call period( DA  (1,J) )     
            call period( DB  (1,J) )
         END DO
!
#if (defined SPMD)
! Latitude distribution for plat=16 on 3 processors.                   !wjp 2007.05
!  procId:      iam=0     |      iam=1     |      iam=2                !wjp 2007.05
!  CAM:    01 02 03 04 05 | 06 07 08 09 10 | 11 12 13 14 15 16         !wjp 2007.05
!  IAP:    16 15 14 13 12 | 11 10 09 08 07 | 06 05 04 03 02 01         !wjp 2007.05
!                left     |     current    |      right                !wjp 2007.05
!jjr right
     src = iam-1
      dest  = iam+1
      if ( mod(iam,npr_y) == 0 ) src = -1
      if ( mod(iam+1,npr_y) == 0 ) dest = -1
      call mp_send3d( commyz, dest, src, NX,  NY,1,              &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, beglatdyn, beglatdyn,1,1,  da )
      call mp_recv3d( commyz, src, NX,  NY,1,                    &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, endlatdynex, endlatdynex,1,1,da )
!jjr left     
      src = iam+1
      dest  = iam-1
      if ( mod(iam,npr_y) == 0 ) dest = -1
      if ( mod(iam+1,npr_y) == 0 ) src = -1
      call mp_send3d( commyz, dest, src, NX, NY,1,               &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, endlatdyn, endlatdyn, 1, 1, rlnt )
      call mp_recv3d( commyz, src, NX,  NY,1,                    &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, beglatdynex, beglatdynex,1,1, rlnt )
      call mp_send3d( commyz, dest, src, NX, NY,1,               &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, endlatdyn, endlatdyn, 1, 1, rdlt )
      call mp_recv3d( commyz, src, NX,  NY,1,                    &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, beglatdynex, beglatdynex,1,1, rdlt )
      call mp_send3d( commyz, dest, src, NX, NY,1,               &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, endlatdyn, endlatdyn, 1, 1, db )
      call mp_recv3d( commyz, src, NX,  NY,1,                    &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, beglatdynex, beglatdynex,1,1, db )

#endif
!------------------------------------------------------------------------------		 	 
!
!     UPDATE T & U & V & Q DUE TO THE HORIZONTAL DIFFUSION
!
!!         DO J = loc_JB, loc_JE
         DO J = beglatdyn, endlatdyn
            if ( j >= JB .and. j <= JE ) then 
               FU1 = FRDU(J,1)
               FU2 = FRDU(J,2)
               FU3 = FRDU(J,3)
               FV1 = FRDV(J,1)
               FV2 = FRDV(J,2)
               FV3 = FRDV(J,3)
               FA1 = FRDP(J,1)
               FA2 = FRDP(J,2)
               FA3 = FRDP(J,3)
               JJ  = J + 1
               J1  = J - 1
               DO I = IB, IE
                  II    = I + 1
                  I1    = I - 1
                  VR0   = VR  (I,J)          ! VR  = 1/rho
                  DS0   = DS  (I,J)          ! DS0 = DS(i', j')
                  DT0   = DT  (I,J)          ! DT0 = DT(i , j )
                  DA0   = DA  (I,J)          ! DA0 = D (i , j )
                  DB0   = DB  (I,J)          ! DB0 = D (i', j')
                  RLNT0 = RLNT(I,J)          ! RLNT0 = 4 * rho(i', j')
                  ROT0  = ROT (I,J)          ! ROT = rho(i,j)
                  RTA0  = DT0 * DA0 * ROT0   ! RTA0 = (rho*DT*D)(i ,j )
                  RSB0  = DS0 * DB0 * RLNT0  ! RSB0 = (rho*DS*D)(i',j') * 4
                  VRU   = VR0 + VR(I1,J)
                  VRV   = VR0 + VR(I ,JJ)
                  RTAI  = ROT (I1,J)*DT(I1,J)*DA(I1,J)  ! RTAI = (rho*DT*D)(i-1,j )
                  RSBI  = RLNT(I,J1)*DS(I,J1)*DB(I,J1)  ! RSBI = (rho*DS*D)(i' ,j'-1) * 4
                  RTAJ  = ROT (I,JJ)*DT(I,JJ)*DA(I,JJ)  ! RTAJ = (rho*DT*D)(i  ,j+1 )
                  RSBJ  = RLNT(II,J)*DS(II,J)*DB(II,J)  ! RSBJ = (rho*DS*D)(i'+1,j'-1) * 4
                  FB1   = FA1 * RDLN(II,J)    ! FB1 = 2 * rho(i'+1, j ) * D(i'+1, j  )
                  FB2   = FA1 * RDLN(I ,J)    ! FB2 = 2 * rho(i'  , j ) * D(i'  , j  )
                  FB3   = FA2 * RDLT(I ,J)    ! FB3 = FA2 * 2 * rho(i, j') * D(i, j' )
                  FB4   = FA3 * RDLT(I,J1)    ! FB4 = FA3 * 2 * rho(i, j'-1) * D(i, j'-1 )
                  TD0   = TK(I,J)
!!                  Q0    = QK(I,J)
!
                  DU    = VRU * ( FU1*(RTA0-RTAI) + (FU2*RSB0 - FU3*RSBI) ) * Kd2    !(5.11)
                  DV    = VRV * ( FV1*(RSBJ-RSB0) - (FV2*RTAJ - FV3*RTA0) ) * Kd2    !(5.12)
                  DdT   = VR0 * ( FB1*(TK(II,J)-TD0) - FB2*(TD0-TK(I1,J))          &
                        +         FB3*(TK(I,JJ)-TD0) - FB4*(TD0-TK(I,J1)) ) * Kd2    !(5.13)
!!                DQ    = VR0 * ( FB1*(QK(II,J)-Q0 ) - FB2*(Q0 -QK(I1,J))          &
!!                      +         FB3*(QK(I,JJ)-Q0 ) - FB4*(Q0 -QK(I,J1)) ) * Kd2    !(5.13)
!
                  U(I,K,J) = UK(I,J) + DU*DTHDFS  
                  V(I,K,J) = VK(I,J) + DV*DTHDFS
                  T(I,K,J) = TW(I,J) + DdT*DTHDFS
!!                Q(I,K,J) = QK(I,J) + DQ*DTHDFS
!   if(j.eq.91.and.k.eq.1.and.i.eq.50) then
!    print*,'test VK',DA0
!    end if
               END DO
            else if ( j == 1 ) then ! at north pole         
               VRN = VR(IB,1 )   
               TDN = TK(IB,1 )
!!             QKN = QK(IB,1 )
!!             DQN = ZERO
               DTN = ZERO
               FAN = FRDP(1 ,1) * VRN
               DO I = IB, IE      
                  EAN = RDLT(I,1 )  ! RDLT(I,1)  = (rho(I,1 ) + rho(I,JB)) * D(I,1)
                  DTN = DTN + EAN*(TK(I,JB) - TDN)  
!!                DQN = DQN + EAN*(QK(I,JB) - QKN)
               END DO
!!      
               T0N = TW(IB,1 ) + DTN * FAN * Kd2 * DTHDFS      !(5.13)-(2)
!!             Q0N = QK(IB,1 ) + DQN * FAN * Kd2 * DTHDFS
               FVI = FRDV(1,1)
               FVJ = FRDV(1,2)
               DO I = IB, IE
                  II        = I + 1
                  T(I ,K,1) = T0N
!!                Q(I ,K,1) = Q0N
                  U(I ,K,1) = ZERO
                  VRV       = VRN + VR(I,JB)
                  RSBI      = RLNT(II,1)*DS(II,1)*DB(II,1) !DB(II,1)=0.5*(D(I,1) + D(I+1,1)) 
                  RSB0      = RLNT(I ,1)*DS(I ,1)*DB(I ,1)
                  RTA0      = ROT (I,JB)*DT(I,JB)*DA(I,JB) !DA(I,JB)=0.5*(D(I,1) + D(I,JB))
                  DVN       = VRV * ( FVI*(RSBI - RSB0) - FVJ*RTA0 ) * Kd2    !(5.12)
                  V(I ,K,1) = VK(I,1) + DVN * DTHDFS
               END DO
!
            else if ( j == NY ) then ! at south pole         
               VRS = VR(IB,NY)
               TDS = TK(IB,NY)
!!             QKS = QK(IB,NY)
!!             DQS = ZERO
               DTS = ZERO
               FAS = FRDP(NY,1) * VRS
               DO I = IB, IE      
                  EAS = RDLT(I,JE)  ! RDLT(I,JE) = (rho(I,JE) + rho(I,NY)) * D(I,JE)
                  DTS = DTS + EAS*(TDS - TK(I,JE))
!!                DQS = DQS + EAS*(QKS - QK(I,JE))
               END DO
!      
               T0S = TW(IB,NY) + DTS * FAS * Kd2 * DTHDFS      !(5.13)-(3)
!!             Q0S = QK(IB,NY) + DQS * FAS * Kd2 * DTHDFS
               DO I = IB, IE
                  II        = I + 1
                  T(I,K,NY) = T0S
!!                Q(I,K,NY) = Q0S
                  U(I,K,NY) = ZERO
                  V(I,K,NY) = ZERO
               END DO
            end if
         END DO
!
!     ENERGY CONSERVATION DUE TO FRICTIONS          2006.4.30
!jjr
!!         DO J =beglatdyn,endlatdyn 
	     DO J = loc_JB, loc_JE
            DO I = IB, IE
               DA(I,J) = UK(I,J) * (U(I,K,J) - UK(I,J)) !DA(I,J)=UK(I,J)* DU=Fu*U
            END DO
         END DO
!--------------------- Set the spherical cyclicity conditions ------------------      
	     DO J = loc_JB, loc_JE
            call period( DA(1,J) )
         END DO
!------------------------------------------------------------------------------		 	 
! 
!   print*,'test VK',VK(50,91),V(50,1,91)
	     DO J = beglatdyn, loc_JE
            DO I = IB, IE
               DB(I,J) = VK(I,J) * (V(I,K,J) - VK(I,J))
            END DO
         END DO
!
#if (defined SPMD)
! Latitude distribution for plat=16 on 3 processors.                 !wjp 2007.05
!  procId:      iam=0     |      iam=1     |      iam=2              !wjp 2007.05
!  CAM:    01 02 03 04 05 | 06 07 08 09 10 | 11 12 13 14 15 16       !wjp 2007.05
!  IAP:    16 15 14 13 12 | 11 10 09 08 07 | 06 05 04 03 02 01       !wjp 2007.05
!                left     |     current    |      right              !wjp 2007.05
!jjr left
       src = iam+1
      dest  = iam-1
      if ( mod(iam,npr_y) == 0 ) dest = -1
      if ( mod(iam+1,npr_y) == 0 ) src = -1
      call mp_send3d( commyz, dest, src, NX, NY,1,               &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, endlatdyn, endlatdyn, 1, 1, db )
      call mp_recv3d( commyz, src, NX,  NY,1,                    &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, beglatdynex, beglatdynex,1,1, db )

!jjr         call mpi_move_left( DB(1,endlatdyn), DB(1,beglatdyn-1),NX)
#endif
 !  print*,'test TB',DB(50,92),DB(50,91),DA(50,92),DA(51,92),DXVPN(92),DXVPS(92),beglatdyn,endlatdyn
	     DO J = loc_JB, loc_JE
            J1     = J - 1 
            DXVPNJ = DXVPN(J)  ! sub. 'stfram' at module 'Dyn_const'
            DXVPSJ = DXVPS(J)
            DO I = IB, IE
               D(I,J)  = HALF * (DA(I,J)+DA(I+1,J)) + DXVPNJ*DB(I,J1) + DXVPSJ*DB(I,J)   !see note P70
               TK(I,J) = T(I,K,J)
            END DO
         END DO
      
	     DO J = loc_JB, loc_JE
            DO I = IB, IE
               T(I,K,J) = TK(I,J) - D(I,J)/CP
            END DO
         END DO
!
      END DO   ! K = 1, NL 
!
!--------------------- Set the spherical cyclicity conditions ------------------      
      DO K = beglev,endlev 
         DO J = beglatdyn, endlatdyn
            call period( U(1,K,J) )
            call period( V(1,K,J) )
            call period( T(1,K,J) )
!!            call period( Q(1,K,J) )
         END DO
      END DO
!------------------------------------------------------------------------------     
	  RETURN
   end subroutine

!================================================================================================
   SUBROUTINE STDFSC
!-----------------------------------------------------------------------------------------------
!  Set constants for computation of the horizontal diffusion
!-----------------------------------------------------------------------------------------------
	  use physconst, only: RAD
      use Dyn_const, only: DLAT, DLON, FIM

      implicit none
!----------------------------------Local workspace----------------------------------
      real(r8) :: SINU(NY), SINV(NY)      
      real(r8) :: YU, YV	  
      real(r8) :: DLONR, DLATR
      real(r8) :: ALN, ALT	  
      real(r8) :: ANT, AVT, ANTS
      real(r8) :: AUT, AVO, AUO	  
      real(r8) :: DNSI, FRTP
      real(r8) :: SNU, SNV 
      real(r8) :: VST	  
      real(r8) :: SVI, SVJ
      real(r8) :: SUI, SUJ
      real(r8) :: FPO	  
	  integer  :: J, K      ! loop index
      integer :: src,dest
!-----------------------------------------------------------------------------------------------

      DO J = 1, NY, JE
         SINU(J) = ZERO
         SINV(J) = ZERO
         DO K = 1, 3
            FRDT(J,K) = ZERO
            FRDS(J,K) = ZERO
            FRDU(J,K) = ZERO
            FRDV(J,K) = ZERO
            FRDP(J,K) = ZERO
         END DO
      END DO
!
!jjr      DO J = loc_JB, loc_JE
      DO J = 2,NY-1 
         YU      = DLAT *  DBLE(J-1)
         YV      = DLAT * (DBLE(J) - HALF)
         SINU(J) = SIN( YU )      ! SINU(J) = sin(sita(j))
         SINV(J) = SIN( YV )      ! SINV(J) = sin(sita(j+1/2))
      END DO
!      
	  SINV(1) = SIN( HALF*DLAT )
      DLONR   = DLON  * RAD
      DLATR   = DLAT  * RAD
      ALN     = ONE   / DLONR
      ALT     = ONE   / DLATR
      ANT     = DLON  / DLAT
      AVT     = DLONR * ANT
      AUT     = AVT   / FOUR
      AVO     = DLONR / FOUR  
      AUO     = DLONR            ! AUO  = a*d(lamda)
      ANTS    = ANT   * ANT      ! ANTS = [d(lamda)/d(sita)]^2
!zhh      DNSI    = DLON*DLON / FIM
      DNSI    = FOUR  * ANTS * SINV(1) * SINV(1) / FIM          !zhh
      FRTP    = FOUR  * ALT  / FIM   ! FRTP = 4/(a*d(sita)*IM)
!
#if (defined SPMD)
! Latitude distribution for plat=16 on 3 processors.             !wjp 2007.05
!  procId:      iam=0     |      iam=1     |      iam=2          !wjp 2007.05
!  CAM:    01 02 03 04 05 | 06 07 08 09 10 | 11 12 13 14 15 16   !wjp 2007.05
!  IAP:    16 15 14 13 12 | 11 10 09 08 07 | 06 05 04 03 02 01   !wjp 2007.05
!                left     |     current    |      right          !wjp 2007.05
!jjr left
!       src = iam+1
!      dest  = iam-1
!      if ( mod(iam,npr_y) == 0 ) dest = -1
!      if ( mod(iam+1,npr_y) == 0 ) src = -1
!      call mp_send3d( dest, src, NX, NY,1,                     &
!                      1, NX, beglatdynex, endlatdynex,1,1,       &
!                      1, NX, endlatdyn, endlatdyn, 1, 1, sinv )
!      call mp_recv3d( src, NX,  NY,1,                             &
!                      1, NX, beglatdynex, endlatdynex,1,1,       &
!                      1, NX, beglatdynex, beglatdynex,1,1, sinv )

!jjr      call mpi_move_left(SINV(endlatdyn), SINV(beglatdyn-1), 1)
#endif
!
!jjr      DO J = loc_JB, loc_JE
      DO J = 2,NY-1

         SNU       = SINU(J)
         FRDT(J,1) = ALN / SNU
         VST       = ALT / SNU
         FRDT(J,2) = VST * SINV(J)
         FRDT(J,3) = VST * SINV(J-1)
      END DO
      FRDT(1 ,1) = -FRTP
      FRDT(NY,1) =  FRTP
!
#if (defined SPMD)
!JJR
! Latitude distribution for plat=16 on 3 processors.             !wjp 2007.05
!  procId:      iam=0     |      iam=1     |      iam=2          !wjp 2007.05
!  CAM:    01 02 03 04 05 | 06 07 08 09 10 | 11 12 13 14 15 16   !wjp 2007.05
!  IAP:    16 15 14 13 12 | 11 10 09 08 07 | 06 05 04 03 02 01   !wjp 2007.05
!                left     |     current    |      right          !wjp 2007.05
!jjr right
!     src = iam-1
!      dest  = iam+1
!      if ( mod(iam,npr_y) == 0 ) src = -1
!      if ( mod(iam+1,npr_y) == 0 ) dest = -1
!      call mp_send3d( dest, src, NX,  NY,1,                     &
!                      1, NX, beglatdynex, endlatdynex,1,1,       &
!                      1, NX, beglatdyn, beglatdyn,1,1,  sinu )
!      call mp_recv3d( src, NX,  NY,1                             &
!                      1, NX, beglatdynex, endlatdynex,1,1,       &
!                      1, NX, endlatdynex, endlatdynex,1,1,sinu )

!jjr      call mpi_move_right(SINU(beglatdyn), SINU(endlatdyn+1), 1)
#endif
!
!jjr      DO J = beglatdyn, loc_JE
      DO J=1,NY-1
         SNV       = SINV(J)
         VST       = ALT / SNV
         FRDS(J,1) = ALN / SNV
         FRDS(J,2) = VST * SINU(J+1)
         FRDS(J,3) = VST * SINU(J)
      END DO
!
!      DO J = loc_JB, loc_JE
     DO J=2,NY-1
         SVI       = SINV(J)
         SVJ       = SINV(J-1)
         FRDU(J,1) = AUO * SINU(J)   ! FRDU(J,1) = a * d(lamda) * sin(sita(j))
         FRDU(J,2) = AUT * SVI*SVI   
!        FRDU(J,2) = a * [d(lamda)]^2 / [4*d(sita)] * [sin(sita(j'))]^2
         FRDU(J,3) = AUT * SVJ*SVJ
!        FRDU(J,3) = a * [d(lamda)]^2 / [4*d(sita)] * [sin(sita(j'-1))]^2
      END DO
!
!jjr      DO J = beglatdyn, loc_JE
      DO J=1,NY-1
         SUI       = SINU(J)
         SUJ       = SINU(J+1)          
         FRDV(J,1) = AVO * SINV(J)    ! FRDV(J,1) = a/4 * d(lamda) * sin(sita(j'))
         FRDV(J,2) = AVT * SUJ*SUJ
!        FRDV(J,2) = a * [d(lamda)]^2 / d(sita) * [sin(sita(j+1))]^2
         FRDV(J,3) = AVT * SUI*SUI
!        FRDV(J,3) = a * [d(lamda)]^2 / d(sita) * [sin(sita(j))]^2
      END DO
!
!jjr      DO J = loc_JB, loc_JE
      DO J=2,NY-1
         SVI       = SINV(J)
         SVJ       = SINV(J-1)
         FPO       = ANTS / SINU(J)
         FRDP(J,1) = ONE
         FRDP(J,2) = FPO  * SVI*SVI*SVI 
!        FRDP(J,2) = [d(lamda)/d(sita)]^2 * [sin(sita(j'))]^3 / [sin(sita(j))]  
         FRDP(J,3) = FPO  * SVJ*SVJ*SVJ
!        FRDP(J,3) = [d(lamda)/d(sita)]^2 * [sin(sita(j'-1))]^3 / [sin(sita(j))]  
      END DO
      FRDP(1 ,1) = + DNSI  
      FRDP(NY,1) = - DNSI
!
!!zhh      CDFS  = DTHDFS * DFS0*DFS0   ! CDFS = deltat*(Kdif)^2
!zhh      DTHDFS = 1800s ;   DFS0 = Kdif = 0.1
!!zhh      DO I = 1, NF
!!zhh         FC(I) = FC(I) * CDFS    ! EQUIVALENCE ( FC(1),FRDU(1,1) )  
!zhh   FRDU(J,1) = AUO * SINU(J) * CDFS                          
!!zhh      END DO
      
!!zhh	  DO K = 1, 3
!!         DO J = 1, NY
!!            FRDP(J,K) = 3.0E0 * FRDP(J,K)  !???????????
!!         END DO
!!zhh      END DO      
!	  
	  RETURN
   end subroutine

end module
