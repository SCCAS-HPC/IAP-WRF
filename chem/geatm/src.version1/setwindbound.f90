module setwindboundary
contains
  subroutine setwindbound(myid, c, sx, ex, sy, ey ,nx,ny)
  integer myid, sx, ex, sy, ey, it
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c
  integer i, j

  if(sx==1)then         ! west boundary
    do j=sy,ey
      c(sx-1,j) = c(sx,j)
    enddo
   endif 

  if(ex==nx)then         ! east boundary
    do j=sy,ey
      c(ex+1,j) = c(ex,j)
    enddo
   endif 

  if(sy==1)then         ! south boundary
    do i=sx,ex
      c(i,sy-1) = c(i,sy)
    enddo
   endif 

  if(ey==ny)then         ! north boundary
    do i=sx,ex
      c(i,ey+1) = c(i,ey)
    enddo
   endif 

  return
  end subroutine

end module setwindboundary
