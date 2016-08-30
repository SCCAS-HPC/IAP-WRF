subroutine extxx (pkcnst, pkdim, fb, kloop)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Copy data to the longitude extensions of the extended array
! 
! Method: 
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!
! $Id: extx.F90,v 1.1.2.5 2004/09/23 17:20:08 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use scanslt,      only: plond, beglatex, endlatex, nxpt, nlonex
  use pmgrid
  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in) :: pkcnst    ! dimension construct for 3-D arrays
  integer , intent(in) :: pkdim     ! vertical dimension
  real(r8), intent(inout) :: fb(plond,beglev:endlev,pkcnst,beglatex:endlatex) ! constituents
  integer,  intent(in) :: kloop ! Limit extent of loop of pkcnst
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i                 ! longitude index
  integer j                 ! latitude  index
  integer k,l                ! vertical  index
  integer nlond             ! extended longitude dim
  integer i2pi              ! start of eastern long. extension
  integer pk                ! k extent to loop over
!-----------------------------------------------------------------------
!
! Fill west edge points.
!
!  pk = pkdim*kloop
  if(nxpt >= 1) then
     do j=beglatex,endlatex
        do i=1,nxpt
         do l=1,kloop
           do k=beglev,endlev
              fb(i,k,l,j) = fb(i+nlonex(j),k,l,j)
           end do
        end do
     end do
     end do
  end if
!
! Fill east edge points
!
  do j=beglatex,endlatex
     i2pi = nxpt + nlonex(j) + 1
     nlond = nlonex(j) + 1 + 2*nxpt
     do i=i2pi,nlond
        do l=1,kloop
        do k=beglev,endlev
           fb(i,k,l,j) = fb(i-nlonex(j),k,l,j)
        end do
     end do
  end do
 end do

  return
end subroutine extxx
