module getboundary
contains
 subroutine getboundnorth1(myid, c, sx, ex, sy, ey ,nx,ny,ig,d)
  integer myid, sx, ex, sy, ey
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c,d
  integer i, j
  if(ey==ny)then         ! north boundary
    do i=sx,ex
      c(i,ey+1) = d(i,ey+1) 
    enddo
   endif 
  return
  end subroutine

 subroutine getboundsouth(myid, c, sx, ex, sy, ey ,nx,ny,ig,d)
  integer myid, sx, ex, sy, ey
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c,d
  integer i, j
  if(sy==1)then         ! south boundary
    do i=sx,ex
      c(i,sy-1) = d(i,sy-1)
    enddo
   endif 
  return
  end subroutine

 subroutine getboundeast(myid, c, sx, ex, sy, ey ,nx,ny,ig,d)
  integer myid, sx, ex, sy, ey
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c,d
  integer i, j
  if(ex==nx)then         ! east boundary
    do j=sy,ey
      c(ex+1,j) = d(ex+1,j)
    enddo
   endif 
  return
  end subroutine
  
 subroutine getboundwest(myid, c, sx, ex, sy, ey ,nx,ny,ig,d)
  integer myid, sx, ex, sy, ey
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c,d
  integer i, j
  if(sx==1)then         ! west boundary
    do j=sy,ey
      IF(IG==11) THEN
         c(sx-1,j) = d(sx-1,j)*0.85
      ELSE    
         c(sx-1,j) = d(sx-1,j)*1.3
      ENDIF   
    enddo
   endif 
  return
  end subroutine

end module getboundary
