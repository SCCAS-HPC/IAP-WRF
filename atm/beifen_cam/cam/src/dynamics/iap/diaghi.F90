SUBROUTINE DIAGHI(IGZ)
!------------------------------------------------------------------------------------------------
! Purpose: Compute the diagnostic variables GHI, Zg & GZ
! Original version : DIAGHI.f (IAP 9L)
! Reconstructed & supplemented : ZhangHe, add the calculation of Zg
! Completed : 2005.9.1
! Update : 2006.6.14, ZhangHe, add the calculation of GZ
!          2008.4.7, ZhangHe, add the calculation of deltap
!          2008.4.30, ZhangHe, add the dumb parameter 'IGZ'
!          October,2012, Jiang Jinrong, 2D parellel
!------------------------------------------------------------------------------------------------
!!   use precision
   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,   only : NX, NY, NL, NZ, IB, IE, period
   use Dyn_const,  only : DSIG, SIGL, PMTOP
   use stdatm,     only : H0B, GHI0
   use IAP_prog,   only : PT, Psa, TT, GHI, GHS, deltap, Zg, GZsm, GZ, PLY
   use physconst,  only : b0, GRAV
   use flexib,     only : kappa2   
   use pmgrid,     only : npr_y, npr_z, myid_z, beglatdyn, endlatdyn,beglev,endlev,endlevp
#if (defined SPMD)
   use mod_comm, only: mp_send3d, mp_recv3d
   use spmd_dyn          , only:  comm_z
   use parutilitiesmodule, only: parcollective3d,sumop
#endif

   implicit none
!----------------------------------------Arguments-----------------------------------------------
   integer, intent(in) :: IGZ  ! index for calculation of GZ
!----------------------------------Local workspace-----------------------------------------------
   real(r8) :: DGH(NX,NL,beglatdyn:endlatdyn)     ! difference of GHI 
   real(r8) :: TT1(NX,NL,beglatdyn:endlatdyn)     ! difference of GHI 
   real(r8) :: GHI1(NX,NZ,beglatdyn:endlatdyn)     ! difference of GHI 
   integer  :: I,J,K,src,dest             ! loop index
!------------------------------------------------------------------------------------------------
   if(npr_z.gt.1) then
#if (defined SPMD)
      do j=beglatdyn,endlatdyn
         do k=1,nl
            do i=1,nx
               TT1(i,k,j)=0.
            end do
        end do
     end do
#endif
   endif
   do j=beglatdyn,endlatdyn
      do k=beglev,endlev
         do i=1,nx
            TT1(i,k,j)=TT(i,k,j)
         end do
      enddo
   end do
   if(npr_z.gt.1) then
#if (defined SPMD)
       call parcollective3d( comm_z, sumop, NX,NL, endlatdyn-beglatdyn+1,TT1)
#endif
   endif
!
   DO J = beglatdyn, endlatdyn
      DO I = 1 ,NX
       GHI1(I,NZ,J) = H0B(I,J) * Psa(I,J) + kappa2 * GZsm(I,J)    !P113 (4)    
!          H0B £½ Rd * TMSA(PSB)/PSB,   TMSA: temperature of standard atmosphere 
      END DO
   END DO
   DO J = beglatdyn, endlatdyn
      DO K = 1 ,NL
         DO I = 1 ,NX
            deltap(I,K,J) = PMTOP / PLY (I,K,J )
            DGH(I,K,J) = b0 * TT1(I,K,J) * (1.-deltap(I,K,J)) / (PT(I,J)*SIGL(K))*DSIG(K)
			                                                                       !  (3.5)
         END DO
      END DO
   END DO
!JJR  need wait
   DO J = beglatdyn, endlatdyn
      DO K = NL,1 ,-1
         DO I = 1 ,NX
            GHI1(I,K,J) = GHI1(I,K+1,J) + DGH(I,K,J)
         END DO
      END DO
   END DO
   do j = beglatdyn, endlatdyn
      do k = beglev , endlev+1
         do i = 1 , NX
        GHI(i,k,j)=GHI1(i,k,j)
         end do
      end do 
   end do

!    Compute geopotential height Zg & GZ
!!   ID = IB - 1
   do j = beglatdyn, endlatdyn
      do k = beglev , endlev+1  
         do i = 1 , NX
            Zg(i,k,j) = ( GHI(i,k,j) + GHI0(i,k,j) ) / GRAV
         end do
      end do
   end do
   if ( IGZ == 1 ) then 
      do j = beglatdyn, endlatdyn
         do k = beglev , endlev+1  
            do i = IB , IE
               GZ(i,k,j) = GHI(i,k,j) + GHI0(i,k,j) - GHS(i,j)    ! 2007.5.14
            end do
            call period(GZ(1,k,j))
         end do
      end do
   end if

   RETURN
END

