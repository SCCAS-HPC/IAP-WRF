module setboundtop_model
contains
      subroutine setboundtop1(myid, c, sx, ex, sy, ey ,nx,ny,topo3,ne,k,ktop)
        integer myid, sx, ex, sy, ey
        real, dimension(sx-1:ex+1,sy-1:ey+1) :: c,topo3
        real, dimension(sx-1:ex+1,sy-1:ey+1) :: pv,ktop
        integer i, j,ne
        do j=sy,ey
         do i=sx,ex
            c(i,j)=topo3(i,j)
         enddo
        enddo
       return
      end subroutine

      subroutine setboundtop(myid, ig,c, sx, ex, sy, ey,nx,ny,globalo3,ne,k,ktop)
        integer myid, sx, ex, sy, ey,ig
        real, dimension(sx-1:ex+1,sy-1:ey+1) :: c,globalo3
        real, dimension(sx-1:ex+1,sy-1:ey+1) :: pv,ktop
        integer i, j, k ,ne
           
        do j=sy,ey
           do i=sx,ex
              !if(k.eq.(ktop(i,j)+1)) then !! the layer above tropopause, by chenhs 
              if(k.ge.ktop(i,j)) then !! the layer above tropopause, by chenhs 
                c(i,j)=globalo3(i,j)
              endif
           enddo
        enddo
           
       return
       end subroutine
end module setboundtop_model                                                                         
      

