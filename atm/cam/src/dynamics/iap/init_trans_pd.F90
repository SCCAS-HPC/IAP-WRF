subroutine init_trans_pd
!----------------------------------------------------------------------------------
! Purpose: Transform initial data from physics grid to dynamical grid.
! Author : ZhangHe
! Completed: 2007.5.7
! Update : ZhangHe, 2007.10.11, transform direction of V-wind from north(v3) to 
!                               south(V)
!          ZhangHe, 2007.12.21, add Qliq & Qice
!          ZhangHe, 2008.4.8, add calculation of P
!          WuJianping, 2008.05, parallel version
!          Juanxiong He, 201008
!          Jiang Jinrong and Zhang He, 2012-10-26, for 2D parallel
!          Zhang He, 2012-11-08, n3m2 --> n3m1
!          Zhang He, 2012-11-09, added is_first_step 
!          Zhang He, 2013-01-25, new mp_send3d, dyn_state is added
!          Zhang He, 2013-03-12, removed transform of Q, Qliq, and Qice
! Reviewed: Zhang He, 2013-03-27
! Modified: Zhang He, 2013-04-16, separated initial run and restart run 
!                     2013-04-17, removed calculation of P
!----------------------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,    only: NX, NY, NL, IB, IE, EX, period  !!
   use pmgrid,      only: plon, plat, plev, beglat, endlat, beglatdyn, endlatdyn,  &
                          beglatdynex, endlatdynex, beglev, endlev  !!
   use IAP_prog,    only: P, U, V, T, PS2, GHS, Q, Qliq, Qice
   use prognostics, only: u3, v3, t3, ps, phis, n3m1, qminus, q3
   use Dyn_const,   only: DXPVN, DXPVS, PMTOP, STFRAM   
   use mathconst,   only: ZERO, HALF
   use constituents, only: pcnst
   use flexib,      only: PHALFC
   use time_manager, only: is_first_step   !!zhh 2012-11-09
   use dynamics_vars,   only: T_FVDYCORE_STATE         !zhh 2013-01-25
   use dyn_internal_state,   only : get_dyn_state      !zhh 2013-01-25
#if (defined SPMD || defined COUP_CSM)
   use mpishorthand, only: mpir8, mpicom
   use spmd_utils,   only: iam, masterproc  !juanxiong he
   use mod_comm,     only: mp_send3d, mp_recv3d   ! jjr
   use pmgrid,       only: npr_y
#endif

   implicit none
!-------------------------------------Local workspace--------------------------------------------
   real(r8), allocatable :: utmp(:,:,:)    ! u wind 
   real(r8), allocatable :: vtmp(:,:,:)    ! v wind 
   type (T_FVDYCORE_STATE), pointer :: dyn_state   !zhh 2013-01-25
   integer  :: I, J, Jp, Jd, K, m, src, dest
   integer  :: commyz     !zhh 2013-01-25
!-----------------------------------------------------------------------------
!
   if (is_first_step()) then    ! for initial run only

      dyn_state => get_dyn_state()
      commyz = dyn_state%grid%commyz
!
      allocate(utmp(NX,beglev:endlev,beglatdynex:endlatdynex))
      allocate(vtmp(NX,beglev:endlev,beglatdynex:endlatdynex))
!
! for 2D variables
      DO Jp = beglat, endlat
         Jd = plat + 1 - Jp
         DO I = 1, plon
            GHS(I+EX,Jd) = phis(I,Jp)
         END DO
         call period( GHS(1,Jd) )
      END DO  
!
!   for 3D variables
      DO Jp = beglat, endlat
         Jd = plat + 1 - Jp
         DO K = beglev,endlev    ! jjr
            DO I = 1, plon
               utmp(I+EX,K,Jd) = u3(I,K,Jp,n3m1)
               vtmp(I+EX,K,Jd) = -v3(I,K,Jp,n3m1)
            END DO
            call period( utmp(1,K,Jd) )
            call period( vtmp(1,K,Jd) )
         END DO
      END DO       
!
!      TRANSFORM (U,V) FROM P-GRID TO (U,V)-GRID
!     
#if (defined SPMD)
      src = iam-1
      dest  = iam+1
      if ( mod(iam,npr_y) == 0 ) src = -1
      if ( mod(iam+1,npr_y) == 0 ) dest = -1
      call mp_send3d( commyz, dest, src, NX, NL, NY,                       &
                      1, NX,beglev,endlev, beglatdynex, endlatdynex,       &
                      1, NX,beglev,endlev, beglatdyn, beglatdyn,  vtmp )
      call mp_recv3d( commyz, src, NX, NL, NY,                             &
                      1, NX, beglev,endlev,beglatdynex, endlatdynex,       &
                      1, NX, beglev,endlev,endlatdynex, endlatdynex, vtmp )
#endif
!
      DO J = beglatdyn, endlatdyn
         if ( J>1 .and. J<plat ) then
            DO K = beglev, endlev
               DO I = IB, IE
                  U(I,K,J) = ( utmp(I,K,J) + utmp(I-1,K,J) ) * HALF
                  V(I,K,J) = vtmp(I,K,J)*DXPVN(J) + vtmp(I,K,J+1)*DXPVS(J)
               END DO
               call period( U(1,K,J) )
               call period( V(1,K,J) )
            END DO

         else if ( J==1 ) then
            DO K = beglev, endlev
               DO I = IB, IE
                  U(I,K,J) = ZERO
                  V(I,K,J) = vtmp(I,K,J)*DXPVN(J) + vtmp(I,K,J+1)*DXPVS(J)
               END DO
               call period( U(1,K,J) )
               call period( V(1,K,J) )
            END DO

         else  ! J = plat       
            DO K = beglev, endlev
               DO I = 1, NX
                  U(I,K,J) = ZERO
                  V(I,K,J) = ZERO
               END DO
            END DO
         end if
      END DO
!
      deallocate(utmp)
      deallocate(vtmp)
!
   end if
!
!  for both initial run and restart run

! for 2D variables
   DO Jp = beglat, endlat
      Jd = plat + 1 - Jp
      DO I = 1, plon
         PS2(I+EX,Jd) = ps(I,Jp,n3m1) / 100.0E0   !zhh 2007.7.19
      END DO
      call period( PS2(1,Jd) )
   END DO  
!!   DO J = beglatdyn, endlatdyn    
!!      DO I = 1, NX
!!         P(I,J) = PS2(I,J) - PMTOP         !zhh 2008.4.8
!!      END DO
!!   END DO       
!
!   for 3D variables
   DO Jp = beglat, endlat
      Jd = plat + 1 - Jp
      DO K = beglev,endlev    ! jjr
         DO I = 1, plon
            T   (I+EX,K,Jd) = t3(I,K,Jp,n3m1)
         END DO
         call period( T   (1,K,Jd) )
      END DO
   END DO       
!
   return
end subroutine init_trans_pd
