module setboundary
contains
 subroutine getboundnorth(myid,c,sx,ex,sy,ey,nx,ny, &
                          sx1,ex1,sy1,ey1,nx1,ny1,nxlo1,nylo1,d)

 integer myid, sx, ex, sy, ey, nx, ny
 integer sx1,ex1,sy1,ey1,nx1,ny1,nxlo1,nylo1
 real, dimension(sx-1:ex+1,sy-1:ey+1) :: c
 real, dimension(sx1-1:ex1+1) :: d
  integer i, j, ii
  if(ey1==ny1)then         ! north boundary
    do i=sx1,ex1
      ii = nxlo1 + i/3
      if(ii.ge.1 .and. ii.le.nx)then
      d(i) = c(ii,nylo1)
      endif
    enddo
   endif
 return
 end subroutine
 
 subroutine setboundnorth(myid, c, sx, ex, sy, ey ,nx,ny,ig,d)
  integer myid, sx, ex, sy, ey,ig
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c
!  real, dimension(sx-1:ex+1) :: d
  real :: d
  integer i, j
  if(ey==ny)then         ! north boundary
    do i=sx,ex
      c(i,ey+1) = d
!      c(i,ey+1) = 30.
    enddo
   endif 
  return
  end subroutine

 subroutine setboundsouth(myid, c, sx, ex, sy, ey ,nx,ny,ig,d)
  integer myid, sx, ex, sy, ey,ig
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c
!  real, dimension(sx-1:ex+1) :: d
  real :: d
  integer i, j
  if(sy==1)then         ! south boundary
    do i=sx,ex     
      c(i,sy-1) = d
    enddo
   endif 
  return
  end subroutine

 subroutine setboundeast(myid, c, sx, ex, sy, ey ,nx,ny,ig,d)
  integer myid, sx, ex, sy, ey,ig
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c
!  real, dimension(sy-1:ey+1) :: d
   real :: d
  integer i, j
  if(ex==nx)then         ! east boundary
    do j=sy,ey
      c(ex+1,j) = d
    enddo
   endif 
  return
  end subroutine
  
 subroutine setboundwest(myid, c, sx, ex, sy, ey ,nx,ny,ig,d,k)
  integer myid, sx, ex, sy, ey,ig
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c
!  real, dimension(sy-1:ey+1) :: d
  real :: d
  integer i, j
  if(sx==1)then         ! west boundary
    do j=sy,ey
      IF(IG==11) THEN
        c(sx-1,j) = d*0.85
      ELSE  
        c(sx-1,j) = d
      ENDIF       
    enddo
   endif 
  return
  end subroutine

 subroutine setboundnorthkpm(myid, c, sx, ex, sy, ey ,nx,ny,d)
  integer myid, sx, ex, sy, ey
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c
  real d
  integer i, j
  if(ey==ny)then         ! north boundary
    do i=sx-1,ex+1
      c(i,ey+1) = d
    enddo
   endif
  return
  end subroutine

 subroutine setboundsouthkpm(myid, c, sx, ex, sy, ey ,nx,ny,d)
  integer myid, sx, ex, sy, ey
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c
  real d
  integer i, j
  if(sy==1)then         ! south boundary
    do i=sx-1,ex+1
      c(i,sy-1) = d
!      c(i,sy-1) = 30.
    enddo
   endif 
  return
 end subroutine

 subroutine setboundeastkpm(myid, c, sx, ex, sy, ey ,nx,ny,d)
  integer myid, sx, ex, sy, ey
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c
  real d
  integer i, j
  if(ex==nx)then         ! east boundary
    do j=sy-1,ey+1
      c(ex+1,j) = d
    enddo
   endif
  return
 end subroutine

 subroutine setboundwestkpm(myid, c, sx, ex, sy, ey ,nx,ny,d)
  integer myid, sx, ex, sy, ey
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c
  real d
  integer i, j
  if(sx==1)then         ! west boundary
    do j=sy-1,ey+1
      c(sx-1,j) = d
!      c(sx-1,j) = 30.
    enddo
   endif 
  return
 end subroutine
  
end module setboundary
