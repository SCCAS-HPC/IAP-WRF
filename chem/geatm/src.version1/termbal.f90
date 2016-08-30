module termbal_model
contains
  subroutine termbal( myid, c0, c ,d, k,sx, ex, sy, ey ,nx, ny, dt)
  integer myid, sx, ex, sy, ey
  real dt
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c0,c,d
  integer i, j


  do j = sy,ey
  do i = sx,ex

   d(i,j) = d(i,j) + c(i,j)-c0(i,j)
   c0(i,j)= c(i,j)
  enddo
  enddo

  return
  end subroutine

  subroutine termballi( myid,c0,c,d1,d2,d3,d4,O1,O2,O3,O4,u,v,k,dx,sx, ex, sy, ey ,nx, ny, dt)
  integer myid, sx, ex, sy, ey
  real dt
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c0,c,d,d1,d2,d3,d4,dx,u,v
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: O1,O2,O3,O4
  real, dimension(sx-1:ex+1,sy-1:ey+1) ::  u_tmp,v_tmp
  integer i, j

!  print*,c0(47,33),c(47,33),k,'+++++'
  do j=sy,ey
    do i=sx,ex
     u_tmp(i,j)=abs(u(i,j))
     v_tmp(i,j)=abs(v(i,j))
    enddo
  enddo



  do j = sy,ey
  do i = sx,ex
    
!   d(i,j) = d(i,j) + c(i,j)-c0(i,j)
!-------------------------------------------------
   if(v(i,j+1)<0.0) then
    d1(i,j)=dt*c0(i,j+1)*((dx(i,j)-u_tmp(i,j+1))*v_tmp(i,j+1))/(dx(i,j)*dx(i,j))+d1(i,j) 
                                                                               ! for north input
   else
     d1(i,j)=d1(i,j)
   endif
!-------------------------------------------------

!-------------------------------------------------
   if(v(i,j-1)>0.0) then
     d2(i,j)=dt*c0(i,j-1)*((dx(i,j)-u_tmp(i,j-1))*v_tmp(i,j-1))/(dx(i,j)*dx(i,j))+d2(i,j) 
                                                                               ! for south input
   else
     d2(i,j)=d2(i,j)
   endif
!-------------------------------------------------

!-------------------------------------------------
   if(u(i-1,j)>0.0) then
     d3(i,j)=dt*c0(i-1,j)*((dx(i,j)-v_tmp(i-1,j))*u_tmp(i-1,j))/(dx(i,j)*dx(i,j))+d3(i,j) 
                                                                                ! for west input
   else
     d3(i,j)=d3(i,j)
   endif
!-------------------------------------------------

!-------------------------------------------------
   if(u(i+1,j)<0.0) then
     d4(i,j)=dt*c0(i+1,j)*((dx(i,j)-v_tmp(i+1,j))*u_tmp(i+1,j))/(dx(i,j)*dx(i,j))+d4(i,j)
                                                                             ! for east input
   else
     d4(i,j)=d4(i,j)
   endif
!-------------------------------------------------

!-------------------------------------------------
  if(u(i,j)<0.0) then
     O4(i,j)=dt*c0(i,j)*((dx(i,j)-v_tmp(i,j))*u_tmp(i,j))/(dx(i,j)*dx(i,j))+O4(i,j)
                                                                           ! for west output
   else
     O3(i,j)=dt*c0(i,j)*((dx(i,j)-v_tmp(i,j))*u_tmp(i,j))/(dx(i,j)*dx(i,j))+O3(i,j)
                                                                            ! for east output
   endif
!-------------------------------------------------
!-------------------------------------------------
   if(v(i,j)<0.0) then
     O2(i,j)=dt*c0(i,j)*((dx(i,j)-u_tmp(i,j))*v_tmp(i,j))/(dx(i,j)*dx(i,j))+O2(i,j)
                                                                        ! for south output
   else
     O1(i,j)=dt*c0(i,j)*((dx(i,j)-u_tmp(i,j))*v_tmp(i,j))/(dx(i,j)*dx(i,j))+O1(i,j)
                                                                         ! for north output
   endif
!-------------------------------------------------



  enddo
  enddo

  return
  end subroutine
  
end module termbal_model
