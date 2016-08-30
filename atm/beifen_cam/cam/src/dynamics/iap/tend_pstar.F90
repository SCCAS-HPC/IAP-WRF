SUBROUTINE tend_pstar( Istar )       
!------------------------------------------------------------------------------------------------
! Purpose: Compute the flexible substitute U*, V*, WS*, Pes* & Pes**
! Author : ZhangHe
! Completed : 2005.9.3
! Update    : 2006.3.25, ZhangHe, added two methods to define flexible substitute (Istar = 3 & 4)
!             2007.4.23, ZhangHe, removed dumb parameter 'SETUV'
!             2007.5.07, ZhangHe, changed subroutine's name from pstartend to tend_pstar
!             2008.4.07, ZhangHe, do not calculating DPes*/Dt when Istar=1    
!             2008.5, WuJianping, parallel version
!             2008.6.11, ZhangHe, available for both serial & parallel
! Modified: Jiang Jinrong, 2012 October, for 2D parallel
! Reviewed: Zhang He, 2012-11-13
! Modified: Zhang He, 2013-01-28, new mp_send3d, dyn_state was added
!           Zhang He, 2013-04-01, use dynamic arrays PXW3 et al.
!------------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,   only: NX, NY, NM, NL, NZ, JB, JE, IB, IE, period
   use Dyn_const 
   use IAP_prog,   only: PT, UT, VT, P, U, V, WS, PLY, Ustar, Vstar,    &
                         Pstar1, Pstar2, WSstar, deltap 
   use tendency,   only: DPstar1
   use mathconst,  only: ZERO, HALF, FOURTH
   use Trans_coef, only: PTU, PTV
   use flexib,     only: IBCFFT
   use sm9h,       only: SM9HAS
   use pmgrid,     only: beglatdyn, endlatdyn, loc_JB, loc_JE,npr_z,&
               beglatdynex,endlatdynex,beglev,endlev,endlevp,npr_y
   use spmd_utils, only: iam
   use dynamics_vars,   only: T_FVDYCORE_STATE         !zhh 2013-01-28
   use dyn_internal_state,   only : get_dyn_state      !zhh 2013-01-28
#if ( defined SPMD )
  use mod_comm, only: mp_send3d, mp_recv3d,mp_send3d_2,mp_recv3d_2
  use spmd_dyn          , only:  comm_z
  use parutilitiesmodule, only: parcollective3d,parcollective,sumop 
#endif

   IMPLICIT NONE
!------------------------------------Arguments---------------------------------------------------
   integer , intent(in) :: Istar     ! index of the methods to defeine flexible substitute
!----------------------------------Local workspace-----------------------------------------------
   real(r8), allocatable :: PXW3(:,:), PXW4(:,:)      ! auxiliary variables
   real(r8), allocatable :: PYP2(:,:), PYM2(:,:)      ! auxiliary variables
   real(r8), allocatable :: DIV2(:,:,:)               ! divergence of [(pes*)¡¤(v*)]
   real(r8), allocatable :: WW(:,:)                   ! wjp, 2008.5
   real(r8) :: DIVX2, DIVY2                  ! divergence of [(pes*)¡¤(v*)]
   real(r8) :: D2(NL), DWK2(NL), DPS2        ! auxiliary variables
   real(r8) :: D1(NL)        ! auxiliary variables
   real(r8) :: D3(NL)        ! auxiliary variables
   real(r8) :: WK6                           ! WK6 = 1 / Pes*
   real(r8) :: PV2(NX)                       ! divergence of [(pes*)¡¤(v*)] at polar
   real(r8) :: W2(NL),W22                    ! auxiliary variables
   real(r8) :: WSK2(NZ), D2K,WSK22           ! auxiliary variables
   real(r8) :: time_begin, time_end
   integer  :: I, J, K,src,dest,src1,dest1                       ! loop index
   type (T_FVDYCORE_STATE), pointer :: dyn_state   !zhh 2013-01-28
   integer  :: commyz     !zhh 2013-01-28

!------------------------------------------------------------------------------------------------
! allocate arrays
   allocate ( WW(NX,NY) )
   allocate ( PXW3(NX,beglatdyn:endlatdyn) )
   allocate ( PXW4(NX,beglatdyn:endlatdyn) )
   allocate ( PYP2(NX,beglatdyn:endlatdyn) )
   allocate ( PYM2(NX,beglatdyn:endlatdyn) )
   allocate ( DIV2(IB:IE,NL,beglatdyn:endlatdyn) )
!
   dyn_state => get_dyn_state()
   commyz = dyn_state%grid%commyz

!---------------------- Renew U & V for new timestep's space difference ---------------------
#if (defined SPMD)
!jjr   call mpi_move_right(PT(1,beglatdyn), PT(1,endlatdyn+1),NX)
   src = iam-1
   dest  = iam+1
   if ( mod(iam,npr_y) == 0 ) src = -1
   if ( mod(iam+1,npr_y) == 0 ) dest = -1
   call mp_send3d( commyz, dest, src, NX,NY,1,                &
                   1, NX, beglatdynex, endlatdynex,1,1,       &
                   1, NX, beglatdyn, beglatdyn,1,1,  pt )
   call mp_recv3d( commyz, src, NX,  NY,1,                    &
                   1, NX, beglatdynex, endlatdynex,1,1,       &
                   1, NX, endlatdynex, endlatdynex,1,1,pt )
   src1 = iam+1
   dest1  = iam-1
   if ( mod(iam,npr_y) == 0 ) dest1 = -1
   if ( mod(iam+1,npr_y) == 0 ) src1 = -1
!
#endif
   DO J = beglatdyn, endlatdyn
      IF (J.GE.JB.AND.J.LE.JE) THEN
         DO I = IB,IE
            PTU(I,J) = HALF * ( PT(I,J) + PT(I-1,J) )
            PTV(I,J) = HALF * ( PT(I,J) + PT(I,J+1) )
         ENDDO
         DO K = beglev ,endlev
            DO I = IB,IE
               U(I,K,J) = UT(I,K,J) / PTU(I,J)
               V(I,K,J) = VT(I,K,J) / PTV(I,J)
            ENDDO
            call period( V(1,K,J) )
            call period( U(1,K,J) )                
         ENDDO
      ELSE IF(J.EQ.1) THEN  ! at north polar 
         DO I = IB,IE
            PTV(I,J) = HALF * ( PT(I,J) + PT(I,J+1) )
         ENDDO
         DO K = beglev ,endlev
            DO I = IB,IE
               U(I,K,1) = ZERO
               V(I,K,J) = VT(I,K,J) / PTV(I,J)
            ENDDO
            call period( V(1,K,J) )
            call period( U(1,K,J) )
         ENDDO
      ELSE                 ! at south polar
         DO K = beglev ,endlev
            DO I = 1 ,NX
               U(I,K,J) = ZERO
               V(I,K,J) = ZERO
            ENDDO
         ENDDO
      ENDIF
   ENDDO
!
!----------------------------- Compute Ustar, Vstar, Pstar2 & deltap --------------------------
   DO J = beglatdyn, endlatdyn
      DO I = 1, NX
         Pstar2(I,J)    = P(I,J)
         DO K = beglev ,endlev
            Vstar(I,K,J)= V(I,K,J)
            Ustar(I,K,J)= U(I,K,J)
         END DO 
      END DO 
   END DO
    
!  IF ( Istar == 1 ) THEN       do not use flexible substitute   

   IF ( Istar == 2 ) THEN   
! --------------- define flexible substitute by 9-POINT HORIZONTAL AREAL SMOOTHING ------------
      call SM9HAS( Pstar2,1 )
      DO K = beglev, endlev
         WW(:,beglatdyn:endlatdyn)=Ustar(:,K,beglatdyn:endlatdyn)
         call SM9HAS( WW, 1 )
         Ustar(:,K,beglatdyn:endlatdyn)=WW(:,beglatdyn:endlatdyn)
!
         WW(:,beglatdyn:endlatdyn)=Vstar(:,K,beglatdyn:endlatdyn)
         call SM9HAS( WW, 2 )
         Vstar(:,K,beglatdyn:endlatdyn)=WW(:,beglatdyn:endlatdyn)
      END DO
   
   ELSE IF ( Istar == 3 ) THEN   
! ------------------------ define flexible substitute mixture with Istar = 1 & 2 --------------
      call SM9HAS( Pstar2,1 )
      DO K = beglev, endlev
         WW(:,beglatdyn:endlatdyn)=Ustar(:,K,beglatdyn:endlatdyn)
         call SM9HAS( WW, 1 )
         Ustar(:,K,beglatdyn:endlatdyn)=WW(:,beglatdyn:endlatdyn)
!
         WW(:,beglatdyn:endlatdyn)=Vstar(:,K,beglatdyn:endlatdyn)
         call SM9HAS( WW, 2 )
         Vstar(:,K,beglatdyn:endlatdyn)=WW(:,beglatdyn:endlatdyn)
      END DO

      DO J = beglatdyn, endlatdyn
         DO I = 1, NX
            Pstar2(I,J)    = alpha(J)*Pstar2(I,J) + betaw(J)*P(I,J)
            DO K = beglev ,endlev
               Vstar(I,K,J)= alpha2(J)*Vstar(I,K,J) + betaw2(J)*V(I,K,J)
               Ustar(I,K,J)= alpha(J)*Ustar(I,K,J) + betaw(J)*U(I,K,J)
            END DO 
         END DO 
      END DO
      
   ELSE IF ( Istar == 4 ) THEN   
! --------------- define flexible substitute by 5-POINT HORIZONTAL arithmetical mean SMOOTHING ------------
#if (defined SPMD)
      call mp_send3d( commyz, dest, src, NX,NY,1,                &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, beglatdyn, beglatdyn,1,1,  p )

      call mp_recv3d( commyz, src, NX,  NY,1,                    &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, endlatdynex, endlatdynex,1,1,p )
      call mp_send3d( commyz, dest1, src1, NX, NY,1,             &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, endlatdyn, endlatdyn, 1, 1, p )

      call mp_recv3d( commyz, src1, NX,  NY,1,                   &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, beglatdynex, beglatdynex,1,1, p )
!
#endif
      DO J = beglatdyn, endlatdyn
         IF (J.GE.JB.AND.J.LE.JE) THEN
            DO I = IB, IE
               Pstar2(I,J) = HALF*P(I,J) + ( P(I-1,J)+P(I+1,J)+P(I,J+1)+P(I,J-1) ) / 8.0
            END DO
         ELSE      ! at the polar
            DO I = IB, IE
               Pstar2(I,J) = P(IB,J)
            END DO 
         END IF
         call period( Pstar2(1,J) )
      END DO
#if (defined SPMD)
!jjr      call mpi_move_left(U(1,1,endlatdyn), U(1,1,beglatdyn-1),NX*NL)
!      call mpi_move_right( U(1,1,beglatdyn), U(1,1,endlatdyn+1),NX*NL)
!      call mpi_move_left(V(1,1,endlatdyn), V(1,1,beglatdyn-1),NX*NL)
!      call mpi_move_right( V(1,1,beglatdyn), V(1,1,endlatdyn+1),NX*NL)
      call mp_send3d_2( commyz, dest, src, NX,  NL,NY,                       &
                        1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                        1, NX,beglev,endlev, beglatdyn, beglatdyn,  u,v )

      call mp_recv3d_2( commyz, src, NX,  NL,NY,                             &
                        1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                        1, NX,beglev,endlev, endlatdynex, endlatdynex,u,v )

      call mp_send3d_2( commyz, dest1, src1, NX,NL, NY,                      &
                        1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                        1, NX, beglev,endlev,endlatdyn, endlatdyn,u,v  )

      call mp_recv3d_2( commyz, src1, NX, NL, NY,                            &
                        1, NX, beglev,endlev,beglatdynex, endlatdynex,       &
                        1, NX,beglev,endlev, beglatdynex, beglatdynex, u,v )
!
#endif
      DO K = beglev ,endlev
         DO J = beglatdyn, endlatdyn
            IF (J.EQ.1) THEN
               DO I = IB, IE
                  Vstar(I,K,1) = HALF * V(I,K,1) + FOURTH * ( V(I-1,K,1) + V(I+1,K,1) )
                  Ustar(I,K,1) = ZERO
               END DO
            ELSE IF (J.EQ.JB) THEN
               DO I = IB, IE
                  Vstar(I,K,JB) = HALF*V(I,K,J)+(V(I-1,K,J)+V(I+1,K,J)+V(I,K,J+1)+V(I,K,J-1)) / 8.0
                  Ustar(I,K,JB) = HALF * U(I,K,JB) + FOURTH * ( U(I-1,K,JB) + U(I+1,K,JB) )
               END DO
            ELSE IF (J.GT.JB.AND.J.LT.JE) THEN
               DO I = IB, IE
                  Vstar(I,K,J) = HALF*V(I,K,J)+(V(I-1,K,J)+V(I+1,K,J)+V(I,K,J+1)+V(I,K,J-1)) / 8.0
                  Ustar(I,K,J) = HALF*U(I,K,J)+(U(I-1,K,J)+U(I+1,K,J)+U(I,K,J+1)+U(I,K,J-1)) / 8.0
		    END DO
            ELSE IF (J.EQ.JE) THEN
               DO I = IB, IE
                  Vstar(I,K,JE) = HALF * V(I,K,JE) + FOURTH * ( V(I-1,K,JE) + V(I+1,K,JE) )
                  Ustar(I,K,JE) = HALF * U(I,K,JE) + FOURTH * ( U(I-1,K,JE) + U(I+1,K,JE) )
		       END DO
            ELSE    ! J = NY
               DO I = 1, NX
                  Vstar(I,K,NY)= ZERO
                  Ustar(I,K,NY)= ZERO
               END DO
            ENDIF
            call period( Vstar(1,K,J) )
            call period( Ustar(1,K,J) )      
         END DO
      END DO

   END IF

! ---------------------------- define deltap ----------------------------  
   DO J = beglatdyn, endlatdyn
      DO K = 1 ,NL
         DO I = 1, NX
            PLY(I,K,J) = P(I,J)*SIGL(K) + PMTOP
            deltap(I,K,J) = PMTOP / PLY (I,K,J)
         END DO
      END DO
   END DO

!************ compute Pstar1 & WSstar ****************
!------------------- COMPUTE FACTORS RELATED TO (I,J) ONLY ----------------------
!======================== zhh 2008.3.22 ===========================
   IF ( Istar == 2 .or. Istar == 3 .or. Istar == 4 ) THEN   
#if (defined SPMD)
!jjr      call mpi_move_left(Pstar1(1,endlatdyn), Pstar1(1,beglatdyn-1),NX)
!jjr      call mpi_move_right( Pstar1(1,beglatdyn), Pstar1(1,endlatdyn+1),NX)
      call mp_send3d( commyz, dest, src, NX,NY,1 ,               &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, beglatdyn, beglatdyn,1,1,  Pstar1 )


      call mp_recv3d( commyz, src, NX,  NY,1 ,                   &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, endlatdynex, endlatdynex,1,1,Pstar1 )
      call mp_send3d( commyz, dest1, src1, NX, NY,1,             &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, endlatdyn, endlatdyn, 1, 1, Pstar1 )
      call mp_recv3d( commyz, src1, NX,  NY,1,                   &
                      1, NX, beglatdynex, endlatdynex,1,1,       &
                      1, NX, beglatdynex, beglatdynex,1,1, Pstar1 )
      call mp_send3d( commyz, dest, src, NX,NL, NY,                        &
                      1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                      1, NX,beglev,endlev, endlatdyn, endlatdyn,  vstar )

#endif
      DO J = beglatdyn, endlatdyn
         IF (J.GE.JB.AND.J.LE.JE) THEN
            DO I = 1 ,NX
               PYP2(I,J)  = RUPH(J) * (Pstar1(I,J+1)+Pstar1(I,J))   
               PYM2(I,J)  = RUMH(J) * (Pstar1(I,J)  +Pstar1(I,J-1))
            END DO
            DO I = IB,IE
               PXW3(I,J)  = OUXH(J) * (Pstar1(I,J)  +Pstar1(I-1,J))
               PXW4(I,J)  = OUXH(J) * (Pstar1(I+1,J)+Pstar1(I-2,J))						 
            END DO
            call period( PXW3 (1,J) )
            call period( PXW4 (1,J) )
         ELSE
            DO I = 1 ,NX
               PYP2(I,J) = ZERO
               PYM2(I,J) = ZERO
               PXW3(I,J) = ZERO
               PXW4(I,J) = ZERO         
            ENDDO
         END IF
      END DO

!----------------------------------- Compute DPes*/Dt & WSstar ------------------------------
!jjr#if (defined SPMD)
!jjr      call mpi_move_left(Vstar(1,1,endlatdyn), Vstar(1,1,beglatdyn-1),NX*NL)
!#endif
#if (defined SPMD)
!jjr      call mpi_move_left(Vstar(1,1,endlatdyn), Vstar(1,1,beglatdyn-1),NX*NL)
      call mp_recv3d( commyz, src1, NX,NL, NY,                             &
                      1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                      1, NX, beglev,endlev,beglatdynex, beglatdynex, vstar)

    if(npr_z.gt.1) then
      DO J = beglatdyn, endlatdyn
            do k=1,NL
                 do I=IB,IE
                 DIV2(I,K,J)=0.
                 end do
            end do
      ENDDO
      DO K=1,NL
          D2(K)=0.
          D3(K)=0.
      ENDDO
    endif
#endif
      DO J = beglatdyn, endlatdyn
         IF (J.GE.JB.AND.J.LE.JE) THEN

          DO K = beglev ,endlev
            DO I = IB,IE
! -------------------------------- COMPUTE DIVERGENCES & THEIR SUMS --------------------------
!jjr               DO K = 1 ,NL
                  DIVX2 = alpha(J)*( PXW3(I+1,J)*Ustar(I+1,K,J)- PXW3(I,J)*Ustar(I,K,J) )      &
                        + betaw(J)*( PXW4(I+2,J)*Ustar(I+2,K,J)- PXW4(I-1,J)*Ustar(I-1,K,J) )/3.0 
					                                                                !P123 (3)
                  DIVY2 = (PYP2(I,J)*Vstar(I,K,J) - PYM2(I,J)*Vstar(I,K,J-1))          !P123 (4)
                  DIV2(I,K,J)  = DIVX2 + DIVY2
            ENDDO !jjr ie
         ENDDO
        ELSE  IF(J.EQ.1) THEN 
             DO K = beglev ,endlev
               DO I = IB,IE
                  PV2(I) = (Pstar1(IB,1)+Pstar1(I,JB)) * Vstar(I,K,1)
               END DO
               D2K  = ZERO
               DO I = IB,IE
                  D2K = D2K  + PV2(I)
               END DO
               D2(K)  = D2K  *  RDLATI     
            ENDDO
       ELSE
          DO K = beglev ,endlev
               DO I = IB,IE
                  PV2(I) = (Pstar1(IB,NY) + Pstar1(I,JE))*Vstar(I,K,JE)
               END DO
                  D2K = ZERO
               DO I = IB,IE
                  D2K = D2K + PV2(I)
               END DO
               D3(K)  = -D2K  *  RDLATI
         ENDDO
     ENDIF
    END DO
    if (npr_z .gt. 1) then
#if defined (SPMD)
       call parcollective(comm_z,sumop,NL,D2)
       call parcollective(comm_z,sumop,NL,D3)
       call parcollective3d( comm_z, sumop, IE-IB+1,NL, endlatdyn-beglatdyn+1,DIV2)
#endif
    endif

      DO J = beglatdyn, endlatdyn
         IF (J.GE.JB.AND.J.LE.JE) THEN
            DO I = IB,IE
! -------------------------------- COMPUTE DIVERGENCES & THEIR SUMS --------------------------
               DO K = 1 ,NL

                  D1(K)   = DIV2(I,K,J)
                  DWK2(K) = DIV2(I,K,J) * DSIG(K)
               ENDDO
               DPS2 = -DWK2(1)                ! for Pes*
               DO K = 2 ,NL
                  DPS2 = DPS2 - DWK2(K)       ! sums for K
               ENDDO
               DPstar1(I,J)  = DPS2                       ! P124  (5)
!			
               WK6 = ONE / Pstar1(I,J)
               DO K = 1 ,NM      ! NM = NL - 1
                  DWK2(K) = WK6*(DPS2 + D1(K))*DSIG(K)
               ENDDO
               WSK22 = ZERO
               DO K = 2 ,NL
                  WSK22    = WSK22 - DWK2(K-1)
                  WSK2(k)=WSK22 !jjr
               end do
               WSK2(1)=0.
               WSK2(NZ)=0.
               do k=beglev,endlev+1
                  WSstar(I,K,J) = WSK2(k)
               ENDDO
            ENDDO       ! end i = IB,IE
!jjr            call period( DPstar1(1,J) )               
!            DO K  = 1 ,NZ
!               call period( WSstar(1,K,J) )               
!            ENDDO

!------------------------------- FOR THE POLAR POINTS -----------------------------------
         ELSE IF(J.EQ.1) THEN         ! at north polar
            W2(1) = ZERO
! -------------------------COMPUTE DIVERGENCE AT NORTH POLE------------------------------
            DPS2 = ZERO
            DO K = 1 ,NL
               DPS2   = DPS2 - D2(K)*DSIG(K)      ! for Pes*
            END DO
!
            DO K = 1 ,NM               ! NM = NL - 1
               DWK2(K)   = (DPS2 + D2(K))*DSIG(K) / Pstar1(IB,1)
            END DO
            DO K = 2 ,NL
               W2(K)     = W2(K-1) - DWK2(K-1)
            END DO
            DO I = IB , IE             !jjr change nx to IE
               DPstar1(I,1) = DPS2
            END DO
              WSstar(:,endlev+1,1)=0.
            DO K = beglev ,endlev+1
               if(k.gt.1.and.k.le.NL) then
               DO I = IB ,IE  !jjr NX -->IE
                  WSstar(I,K,1)  = W2(K)
               END DO
               endif
            END DO
!
!  -------------------------------- FOR THE POLAR POINTS --------------------------------------
         ELSE                        ! J = NY  South polar
            W22 = ZERO
! ------------------------- COMPUTE DIVERGENCE AT SOUTH POLE ---------------------------
            DPS2  = ZERO
            DO K = 1 ,NL
               DPS2   = DPS2 - D3(K)*DSIG(K)
            END DO
            DO K = 1 ,NM               ! NM = NL - 1
               DWK2(K)   = (DPS2 + D3(K))*DSIG(K) / Pstar1(IB,NY)
            END DO
            DO K = 2 ,NL
               W22     = W22 - DWK2(K-1)
               W2(k)=W22 !jjr
            END DO
            DO I = IB , IE   !jjr NX -->IE
               DPstar1(I,NY) = DPS2
            END DO
              WSstar(:,endlev+1,NY)=0.
            DO K = beglev ,endlev+1
               if(k.gt.1.and.k.le.NL) then
               DO I = IB ,IE   !jjr NX-->IE
                  WSstar(I,K,NY)  = W2(K)
               END DO
               endif
            END DO
         END IF
      END DO 
           DO J=beglatdyn,endlatdyn
            call period( DPstar1(1,J) )
            DO K  = beglev ,endlev+1
               call period( WSstar(1,K,J) )
            ENDDO
          ENDDO



!     HIGH-MID LAT FILTER  FOR DPstar1/DT 
!!   write (*,*) 'before filter : DPstar1(206,2) =', DPstar1(206,2)
!!   CALL CPU_TIME ( time_begin )
! ===================== zhh 2008.6.16 =====================
      WW(:,beglatdyn:endlatdyn)=DPstar1(:,beglatdyn:endlatdyn)   !zhh 
!!      CALL FILT2D( DPstar1,0,1,IBCFFT ) 
      CALL FILT2D( WW,0,1,IBCFFT )		   !zhh 2008.6.16
      DPstar1(:,beglatdyn:endlatdyn)=WW(:,beglatdyn:endlatdyn)
! ===================== zhh 2008.6.16 =====================
!!   CALL CPU_TIME ( time_end )
!!   PRINT *, 'run time of FILT2D-Pstar was ', time_end - time_begin, 'sec'
!!   pause ' pstartend '
!!   write (*,*) 'after filter  : DPstar1(206,2) =', DPstar1(206,2)
   ELSE
      DO J = beglatdyn, endlatdyn
         DO I = 1, NX
            Pstar1(I,J) = P(I,J)
            DO K = beglev ,endlev+1
               WSstar(I,K,J) = WS(I,K,J)
            END DO 
         END DO 
      END DO
   END IF
!======================== zhh 2008.3.22 ===========================

! deallocate arrays
   deallocate(WW)
   deallocate(PXW3)
   deallocate(PXW4)
   deallocate(PYP2)
   deallocate(PYM2)
   deallocate(DIV2)
   
   return
END 
