!-----------------------------------------------------------------------------
!   Processing some variables from flux coupler.
!-----------------------------------------------------------------------------
!
!        By Yongqiang Yu , 16 Apr., 1999
!
!
!
      SUBROUTINE post_cpl

#include <def-undef.h>


use param_mod
use pconst_mod
use tracer_mod
use forc_mod
use buf_mod
use control_mod
use shr_sys_mod
!use work_mod,only : wkj
use output_mod,only:spval
!
      implicit none
!
        tsf=0.0D0
        ssf=0.0D0
        swv=0.0D0
        su=0.0D0
        sv=0.0D0

!$OMP PARALLEL DO PRIVATE (i,j)
        do j=1,jmt
!        do i= 2,imt-1 ! 1,imt !LPF 20120822
        do i= 1,imt
           ! net heat flux
           TSF(i,j) = vit(i,j,1)*(lat1(i,j)+sen(i,j)+lwup(i,j)+lwdn(i,j )+netsw(i,j)+melth(i,j)) *OD0CP  
           ! net solar radiation
           SWV(i,j) = vit(i,j,1)*netsw(i,j) 
           ! none solar radiation !for BUOY
           NSWV(i,j)= vit(i,j,1)*(lat1(i,j)+sen(i,j)+lwup(i,j)+lwdn(i,j) +melth(i,j))
           SSF(i,j) =-vit(i,j,1)*(prec(i,j)+evap(i,j)+ meltw(i,j)+roff(i,j))  *34.7*1.0e-3/DZP(1)*OD0                                     ! P+E+melting !linpf 25->DZP(1)
        end do
        end do
        
        call exchange_2d(tsf,1,1)
        call exchange_2d(swv,1,1)
        call exchange_2d(ssf,1,1)
        call exchange_2d(roff,1,1)

!for ocn output
        lthf  = lat1        !latent flux
        sshf  = sen         !sensible flux
        lwv   = lwup+lwdn   !long wave flux
        fresh = ssf
        runoff= roff

        where(vit(:,:,1)<0.5) tsf=spval
        where(vit(:,:,1)<0.5) swv=spval
        where(vit(:,:,1)<0.5) nswv=spval
        where(vit(:,:,1)<0.5) ssf=spval

        where(vit(:,:,1)<0.5) lthf=spval
        where(vit(:,:,1)<0.5) sshf=spval
        where(vit(:,:,1)<0.5) lwv=spval
        where(vit(:,:,1)<0.5) fresh=spval
        where(vit(:,:,1)<0.5) runoff=spval
!
#ifdef USE_OCN_CARBON
        call exchange_2d(pco2)
#endif         
!
! surface stress
!$OMP PARALLEL DO PRIVATE (i,j)
        do j=1,jmt-1
        do i=2,imt
           SU(i,j) = (taux(i-1,j)+taux(i,j)+taux(i-1,j+1)+taux(i,j+1))*0.25
        end do
        end do
!$OMP PARALLEL DO PRIVATE (i,j)
        do j=2,jmt-1
        do i=2,imt
            SV(i,j) = (tauy(i-1,j)+tauy(i,j)+tauy(i-1,j+1)+tauy(i,j+1))*0.25 !linpf
        end do
        end do
!
        call exchange_2d(SU,1,1)
        call exchange_2d(SV,1,1)

        where(viv(:,:,1)<0.5) su=spval
        where(viv(:,:,1)<0.5) sv=spval

!calculate USTAR
      DO J = 1,jmt
         DO I = 1,imt
          USTAR(I,J)=sqrt(sqrt(taux(i,j)*taux(i,j)+tauy(i,j)*tauy(i,j))*OD0)*vit(i,j,1) 
         END DO
      END DO        
      
!$OMP PARALLEL DO PRIVATE (i,j)
       do i=1,imt
       do j=1,jmt
          licomqice (i,j)= 0.0D0
       end do
       end do

!      call chk_var2d(taux,"us",1)
!      call chk_var2d(tauy,"vs",1)
!      call chk_var2d(SU,"SU",0)
!      call chk_var2d(SV,"SV",0)
!      call chk_var2d(SWV,"sw",1)
!      call chk_var2d(NSWV,"sw",1)
!      call chk_var2d(USTAR,"US",1)

        return
        end
