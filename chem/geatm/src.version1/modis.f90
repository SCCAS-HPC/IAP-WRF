module modis
contains
      SUBROUTINE MODISLAND(myid,LANDUSE,sx,ex,sy,ey,ne)
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccc    TO CONVERT MODIS LANDUSE TO USGS LANDUSE    ccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer :: sx,ex,sy,ey,ne,myid
      REAL,dimension(sx-1:ex+1,sy-1:ey+1) :: LANDUSE
        

      DO I = sx, ex
       DO J = SY, EY 
      select case(int(landuse(I,J)))
      case(0)
        landuse(I,J) = 0 ! no data
      case(1)
        landuse(I,J) = 14. 
      case(2)
        landuse(I,J) = 13.
      case(3)
        landuse(I,J) = 12.
      case(4)
        landuse(I,J) = 11.
      case(5)
        landuse(I,J) = 15.
      case(6)
        landuse(I,J) = 8.
      case(7)
        landuse(I,J) = 8.
      case(8)
        landuse(I,J) = 10.
      case(9)
        landuse(I,J) = 10.
      case(10)
        landuse(I,J) = 10.
      case(11)
        landuse(I,J) = 7.
      case(12)
        landuse(I,J) = 2.
      case(13)
        landuse(I,J) = 1.
      case(14)
        landuse(I,J) = 5.
      case(15)
        landuse(I,J) = 24.
      case(16)
        landuse(I,J) = 19.      
      case(17)
        landuse(I,J) = 16.
      case(18)
        landuse(I,J) = 21.
      case(19)
        landuse(I,J) = 22.
      case(20)
        landuse(I,J) = 23.
      endselect
      ENDDO
      ENDDO
 

      RETURN
      END SUBROUTINE
end module modis
