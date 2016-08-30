!
!-----------------------------------------------------------------------------
!   Preparing oceanic variables for communication between OGCM and coupler.
!-----------------------------------------------------------------------------
!
!        By Yongqiang Yu , 16 Apr., 1999
!

module fluxcpl

#include <def-undef.h>   

      use param_mod
      use pconst_mod
      use buf_mod
      use tracer_mod
      use dyn_mod
      use cdf_mod
      use control_mod
      use msg_mod,only : nproc
      use shr_const_mod,only:SHR_CONST_SPVAL !linpf 20120816
#ifdef USE_OCN_CARBON      
      use coutput_mod, only : uptake
#endif      
!
!
contains

  SUBROUTINE flux_cpl
!
        implicit none
!
!$OMP PARALLEL DO PRIVATE (J,I)
        do j=1,jmt
        do i=1,imt
           if (licomqice(i,j) .gt. 0.0) then
               q(i,j)= licomqice(i,j)*D0*CP*DZP(1)/86400.
           else
               q(i,j)= (tbice-at(i,j,1,1))*D0*CP*DZP(1)/86400.  
           endif
        end do
        end do
!
!$OMP PARALLEL DO PRIVATE (J,I)
        T_CPL = AT(:,:,1,1) + 273.15
        S_CPL = AT(:,:,1,2)*1000. + 35.
        U_CPL = 0.
        V_CPL = 0.
        call exchange_2d(t_cpl,1,1)
        call exchange_2d(s_cpl,1,1)
!linpf 2012Jul26
        do j=1,jmt
         do i=1,imt
           if(vit(i,j,1)<0.5) t_cpl(i,j)=SHR_CONST_SPVAL  
           if(vit(i,j,1)<0.5) s_cpl(i,j)=SHR_CONST_SPVAL
         enddo
        enddo
 
       do j=2,jmt !ny
        do i=1,imt-1 !nx
           U_CPL(i,j)  =   (U(i,j,1)+U(i+1,j,1)+U(i,j-1,1)+U(i+1,j-1,1))*.25*vit(i,j,1)
           V_CPL(i,j)  =  -(V(i,j,1)+V(i+1,j,1)+V(i,j-1,1)+V(i+1,j-1,1))*.25*vit(i,j,1)
        end do
        end do
        if(mytid == 0) then
          do i = 1, imt - 1
             U_CPL(i, 1) = (U(i,1,1)+U(i+1,1,1))*.25 * vit(i,1,1)
          enddo
        endif
!
!$OMP PARALLEL DO PRIVATE (J,I)
        do j=2,jmt-1
        do i=2,imt-1
           dhdx (i,j)  =   (h0(i+1,j)-h0(i-1,j)) * 0.5*OTX(j)
           dhdy (i,j)  =   (h0(i,j+1)-h0(i,j-1)) * 0.5*OUY(j)
        end do
        end do

        ! TODO, consider here
        if(mytid == 0) then
          do i = 1, imt
            dhdy (i,1)  =   0.
          enddo
        else if (mytid == nproc-1) then
          do i=1,imt
           dhdy (i,jmt)  =   0.0
          end do
        end if
!linpf 2012Sep04
        do j=1,jmt
         do i=1,imt
           if(vit(i,j,1)<0.5) dhdx(i,j)=0.0D0 !SHR_CONST_SPVAL  
           if(vit(i,j,1)<0.5) dhdy(i,j)=0.0D0 !SHR_CONST_SPVAL
         enddo
        enddo
!
!
!$OMP PARALLEL DO PRIVATE (J,I)
        do i=1,imt
        do j=1,jmt
           licomqice(i,j)=0.0 
        end do
        end do
        
#ifdef USE_OCN_CARBON
      co2_cpl(:,:) = uptake(:,:)
#endif          
        return

  END subroutine flux_cpl


end module fluxcpl

