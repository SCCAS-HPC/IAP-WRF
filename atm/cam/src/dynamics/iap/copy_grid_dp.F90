subroutine copy_grid_dp( Itrans, t3, u3, v3, q3, omga, ps, phis, n3, npt )
!---------------------------------------------------------------------
! Purpose: copy variables from dynamical structure to prognostics structure.
! Author : ZhangHe
! Completed: 2007.4.12
! Update : ZhangHe, 2007.10.11
!          ZhangHe, 2007.12.21
!          ZhangHe, 2008.04.20, add dumb parameter 'Itrans'
!          WuJianping, 2008.5, parallel version
!          ZhangHe, 2008.6.11, available for both serial & parallel
!          He Juanxiong, 2010.08, renamed from trans_grid_dp
!          Jiang Jinrong, 2012 OCT, for 2D parallel version
!          Zhang He, 2012-11-08, intent(out) :: t3,u3,... --> intent(inout) :: t3,u3,... 
!          Zhang He, 2013-03-12, removed transform of Q, Qliq, and Qice
!..................................................................................
!                                                                                 :
!             ==========================================                          :
!              The horizontal distribution of variables                           :
!                      C-grid system (dynamics)                                   :
!             ==========================================                          :
!                                                                                 :
!                       i - 1               i               i + 1    
!                                                                          
!                |        |        |        |        |        |
!                |        |        |        |        |        |
!    j - 1  ---- U ------ P ------ U ------ P ------ U ------ P ---  
!                |        |        |        |        |        |
!                |        |        |        |        |        |
!                |        |        |        |        |        |
!           ---- * ------ V ------ * ------ V ------ * ------ V --- j' - 1
!                |        |        |        |        |        |
!                |        |        |        |        |        |
!                |        |        |        |        |        |
!      j    ---- U ------ P ------ U ------ P ------ U ------ P ---   
!                |        |        |        |        |        |
!                |        |        |        |        |        |
!                |        |        |        |        |        |
!           ---- * ------ V ------ * ------ V ------ * ------ V ---   j' 
!                |        |        |        |        |        |
!                |        |        |        |        |        |
!                |        |        |        |        |        |
!    j + 1  ---- U ------ P ------ U ------ P ------ U ------ P ---    
!                |        |        |        |        |        |
!                |        |        |        |        |        |
!                |        |        |        |        |        |
!           ---- * ------ V ------ * ------ V ------ * ------ V --- j' + 1
!                |        |        |        |        |        |
!                |        |        |        |        |        |
!
!              i'- 1               i'              i'+ 1                     N
!                                                                            |
!                                                                            |
!                                                                      W <--------> E
!                                                                            |
!                                                                            |
!                                                                            S
!
!                                                                   
!       U CARRYED ON "U" GRID   V ON "V"   F(Q,P,T...) ON "P"       
!                                                                   
!      NOTE: V DOES NOT EXIST AT J'   = NY                          
!            FOR IAP   U(J = 1 OR NY) = 0.0E0  V(J' = NY) = 0.0E0   
!                                                                   
!......................................................................................

   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,   only: EX, JB, JE, NZ, NX, NY, NL
   use pmgrid,     only: plon, plat, plev, beglat, endlat, beglatdyn, endlatdyn,  &
                         beglev, endlev, beglatdynex, endlatdynex, npr_y, npr_z
   use IAP_prog,   only: U, V, T, P, WPV, Zg, Q, Qliq, Qice    !zhh 2013-03-12
   use Dyn_const,  only: DXVPN, DXVPS, PMTOP   
   use physconst,  only: GRAV
   use mathconst,  only: ZERO, HALF
   use constituents,  only: pcnst      !juanxiong he, 201008
   use spmd_utils, only: iam   !zhh 2013-03-13
   use dynamics_vars,   only: T_FVDYCORE_STATE         !zhh 2013-03-13
   use dyn_internal_state,   only : get_dyn_state      !zhh 2013-03-13
#if ( defined SPMD )
   use mod_comm, only: mp_send3d, mp_recv3d
   use spmd_dyn          , only:  comm_z
   use parutilitiesmodule, only: parcollective2d,parcollective,sumop
#endif

   implicit none

!----------------------------------------Arguments-----------------------------------------------
   integer,  intent(in)  :: Itrans                   ! Itrans = 1, transform all input vars;
!                                                      Itrans = 2, transform q and ps only
   integer,  intent(in)  :: n3                       ! current time level
   integer,  intent(in)  :: npt                      ! number of time levels in the dycore
   real(r8), intent(inout) :: t3(plon,beglev:endlev,beglat:endlat,npt)   ! temperature
   real(r8), intent(inout) :: u3(plon,beglev:endlev,beglat:endlat,npt)   ! u-wind component
   real(r8), intent(inout) :: v3(plon,beglev:endlev,beglat:endlat,npt)   ! v-wind component
   real(r8), intent(inout) :: q3(plon,beglev:endlev,pcnst,beglat:endlat,npt)   ! specific humidity, juanxiong he, 201008
   real(r8), intent(inout) :: omga(plon,beglev:endlev,beglat:endlat)     ! p-surface vertical velocity
   real(r8), intent(inout) :: ps(plon,beglat:endlat,npt)        ! surface pressure
   real(r8), intent(inout) :: phis(plon,beglat:endlat)          ! surface geopotential
!-------------------------------------Local workspace--------------------------------------------
   integer  :: I, Jd, Jp, K, src1, dest1
   type (T_FVDYCORE_STATE), pointer :: dyn_state   !zhh 2013-03-13
   integer  :: commyz     !zhh 2013-03-13
!---------------------------------------------------------------------------------------
!
   dyn_state => get_dyn_state()
   commyz = dyn_state%grid%commyz
!
#if (defined SPMD)
! Latitude distribution for plat=16 on 3 processors.
!  procId:      iam=0     |      iam=1     |      iam=2          !wjp 2007.05
!  CAM:    01 02 03 04 05 | 06 07 08 09 10 | 11 12 13 14 15 16   !wjp 2007.05
!  IAP:    16 15 14 13 12 | 11 10 09 08 07 | 06 05 04 03 02 01   !wjp 2007.05
!                left     |     current    |      right          !wjp 2007.05
!!   call mpi_move_left(V(1,1,endlatdyn), V(1,1,beglatdyn-1),plon*plev)
   src1 = iam+1
   dest1  = iam-1
   if ( mod(iam,npr_y) == 0 ) dest1 = -1
   if ( mod(iam+1,npr_y) == 0 ) src1 = -1
   call mp_send3d( commyz, dest1, src1,NX,NL ,NY,                       &
                   1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                   1, NX,beglev,endlev, endlatdyn, endlatdyn,  V )
   call mp_recv3d( commyz, src1, NX,NL,  NY,                            &
                   1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                   1, NX, beglev,endlev,beglatdynex, beglatdynex, V )
#endif
	   
   DO Jd = beglatdyn, endlatdyn    ! Scan latitude from north pole to south pole
!
!   transform latitude loop index of dynamics 'Jd' to that of physics 'Jp'
      Jp = plat + 1 - Jd
!
      DO I = 1, plon
         ps(I, Jp, n3) = ( P(I+EX,Jd)+PMTOP ) * 100.0E0
      END DO
!
!      SET UP MOISTURE WORKSHOPS
!
      DO K = beglev, endlev
         DO I = 1, plon
!            q3  (I, K,1,Jp, n3) = Q(I+EX,K,Jd)
!            q3  (I, K,2,Jp, n3) = Qliq(I+EX,K,Jd)
!            q3  (I, K,3,Jp, n3) = Qice(I+EX,K,Jd)
         END DO
      END DO
!
      if ( Itrans == 1 ) then 
!
         DO I = 1, plon
            phis(I,JP) = 0.
            if(endlev.eq.NL)   phis(I, Jp) = Zg(I+EX,NZ,Jd) * GRAV    !zhh 2007.12.5
         END DO
!
!      SET UP TEMPERATURE/vertical motion WORKSHOPS
!
         DO K = beglev, endlev
            DO I = 1, plon
               t3  (I, K,  Jp, n3) = T(I+EX,K,Jd)
               omga(I, K,  Jp    ) = WPV(I+EX,K,Jd) * 100.0E0
            END DO
         END DO
!
!      SET UP WIND WORKSHOPS ON physics GRIDS
!
         IF ( Jd.GE.JB .AND. Jd.LE.JE ) THEN
	        DO K = beglev, endlev
               DO I = 1, plon
                  u3(I, K, Jp, n3) = HALF * ( U(I+EX,K,Jd) + U(I+EX+1,K,Jd) )
                  v3(I, K, Jp, n3) = - (DXVPN(Jd)*V(I+EX,K,Jd-1) + DXVPS(Jd)*V(I+EX,K,Jd))
               END DO
            END DO
	     ELSE    ! at north and south pole
	        DO K = beglev, endlev
               DO I = 1, plon
                  u3(I, K, Jp, n3) = ZERO
                  v3(I, K, Jp, n3) = ZERO
               END DO
            END DO      
	     ENDIF
!
      end if
!
   END DO

   if (npr_z.gt.1) then
#if (defined SPMD)
      if ( Itrans == 1 ) then
         call parcollective(comm_z,sumop,plon,endlat-beglat+1,phis) !jjr
      endif
#endif
   endif
      
   return
end subroutine copy_grid_dp
