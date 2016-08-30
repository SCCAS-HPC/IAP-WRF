subroutine trans_grid_pd( fu, fv, t2, qminus, beglat, endlat, pcnst )
!----------------------------------------------------------------------------------
! Purpose: Transform tendencies from physics grid to dynamical grid.
! Author : ZhangHe
! Completed: 2007.4.11
! Update : ZhangHe, 2007.10.11
!          ZhangHe, 2007.12.21
!          WuJianping, 2008.5, parallel version
!          ZhangHe, 2008.6.11, available for both serial & parallel
!          He Juanxiong, 2010.08, revised	
! Modified: Jiang Jinrong, 2012 October, for 2D parallel
! Reviewed: Zhang He, 2012-11-13, removed the revision from He juanxiong in AUG 2010
! Modified: Zhang He, 2013-01-28, new mp_send3d, dyn_state is added
!           Zhang He, 2013-02-06, removed redundant variables in use only statement
!           Zhang He, 2013-03-12, removed transform of Q, Qliq, and Qice
!----------------------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,  only: NX, NY, NL, EX, IB, IE, period
   use pmgrid,    only: plon, plat, plev, beglatdyn, endlatdyn, beglev,endlev,beglatdynex,endlatdynex
   use IAP_prog,  only: Q, Qliq, Qice
   use prognostics, only: n3, q3        !zhh 2013-04-15
   use tendency,  only: SU, SV, ST
   use Dyn_const, only: DXPVN, DXPVS, PMTOP   
   use mathconst, only: ZERO, HALF
   use spmd_utils, only: iam
   use dynamics_vars,   only: T_FVDYCORE_STATE         !zhh 2013-01-28
   use dyn_internal_state,   only : get_dyn_state      !zhh 2013-01-28
#if (defined SPMD)
      use mod_comm, only: mp_send3d, mp_recv3d
      use pmgrid,     only : npr_y
#endif

   implicit none
!----------------------------------------Arguments-----------------------------------------------
   integer,  intent(in) :: beglat
   integer,  intent(in) :: endlat
   integer,  intent(in) :: pcnst
   real(r8), intent(in) :: fu (plon,beglev:endlev,beglat:endlat)
   real(r8), intent(in) :: fv (plon,beglev:endlev,beglat:endlat)
   real(r8), intent(in) :: t2 (plon,beglev:endlev,beglat:endlat)
   real(r8), intent(in) :: qminus(plon,beglev:endlev,pcnst,beglat:endlat)
!-------------------------------------Local workspace--------------------------------------------
   real(r8), allocatable :: futmp(:,:,:)    ! u wind tendency
   real(r8), allocatable :: fvtmp(:,:,:)    ! v wind tendency
   integer  :: I, Jp, Jd, K, J,m,dest,src
   type (T_FVDYCORE_STATE), pointer :: dyn_state   !zhh 2013-01-28
   integer  :: commyz     !zhh 2013-01-28
   
!-----------------------------------------------------------------------------
!
   dyn_state => get_dyn_state()
   commyz = dyn_state%grid%commyz
!
   allocate(futmp(NX,beglev:endlev,beglatdynex:endlatdynex))
   allocate(fvtmp(NX,beglev:endlev,beglatdynex:endlatdynex))
!
   DO Jp = beglat, endlat
      Jd = plat + 1 - Jp
      DO K = beglev,endlev 
         DO I = 1, plon
            futmp(I+EX,K,Jd) = fu(I,K,Jp)
            fvtmp(I+EX,K,Jd) =-fv(I,K,Jp)
            ST   (I+EX,K,Jd) = t2(I,K,Jp)
            do m = 1, pcnst
!!            do m = 1, 3
               q3(i,k,m,jp,n3) = qminus(i,k,m,jp)
            end do
!!            Q    (I+EX,K,Jd) = qminus(I,K,1,Jp)
!!            Qliq (I+EX,K,Jd) = qminus(I,K,2,Jp)
!!            Qice (I+EX,K,Jd) = qminus(I,K,3,Jp)
         END DO
         call period( futmp(1,K,Jd) )
         call period( fvtmp(1,K,Jd) )
         call period( ST   (1,K,Jd) )
!!         call period( Q    (1,K,Jd) )
!!         call period( Qliq (1,K,Jd) )
!!         call period( Qice (1,K,Jd) )
      END DO
   END DO       
!
!      TRANSFORM (DU,DV) FROM P-GRID TO (U,V)-GRID
!     
#if (defined SPMD)
!JJR
      src = iam-1
      dest  = iam+1
      if ( mod(iam,npr_y) == 0 ) src = -1
      if ( mod(iam+1,npr_y) == 0 ) dest = -1
      call mp_send3d( commyz, dest, src, NX,  NL,NY,                       &
                      1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                      1, NX,beglev,endlev, beglatdyn, beglatdyn,  fvtmp )
      call mp_recv3d( commyz, src, NX,  NL,NY,                             &
                      1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                      1, NX,beglev,endlev, endlatdynex, endlatdynex,fvtmp )

#endif
!
   DO J = beglatdyn, endlatdyn
      if ( J>1 .and. J<plat ) then
         DO K = beglev,endlev 
            DO I = IB, IE
               SU(I,K,J) = ( futmp(I,K,J) + futmp(I-1,K,J) ) * HALF
               SV(I,K,J) = fvtmp(I,K,J)*DXPVN(J) + fvtmp(I,K,J+1)*DXPVS(J)
            END DO
            call period( SU(1,K,J) )
            call period( SV(1,K,J) )
         END DO

      else if ( J==1 ) then
         DO K = beglev,endlev 
            DO I = IB, IE
               SU(I,K,J) = ZERO
               SV(I,K,J) = fvtmp(I,K,J)*DXPVN(J) + fvtmp(I,K,J+1)*DXPVS(J)
            END DO
            call period( SU(1,K,J) )
            call period( SV(1,K,J) )
         END DO

      else  ! J = plat       
         DO K = beglev,endlev 
            DO I = 1, NX
               SU(I,K,J) = ZERO
               SV(I,K,J) = ZERO
            END DO
         END DO
      end if
   END DO
!
   deallocate(futmp)
   deallocate(fvtmp)
   
   return
end subroutine trans_grid_pd
