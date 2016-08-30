module eddyz_model
contains
                 subroutine eddyz(myid,UST,PBLH,RMOL,HEIGHT,TERR,kz,sx,ex,sy,ey,k)
                     integer :: myid,sx,ex,sy,ey
                     real :: karman,gama, beta
                     real,dimension(sx-1:ex+1,sy-1:ey+1) :: UST,PBLH,RMOL,HEIGHT,TERR
                     real :: Z,ZSL,PBL
                     real :: kzmin,kzmax
                     real,dimension(sx-1:ex+1,sy-1:ey+1) :: kz
                      karman= 0.4
                      gama  = 19.3
                      beta = 6.0
                      kzmin = 0.5
                      kzmax = 500.
                      
                     do j=sy,ey
                     do i=sx,ex
                     
                      PBL = PBLH(i,j) 
                     
                     IF(HEIGHT(i,j)<100000.) THEN  !! by chnehs
                       Z=HEIGHT(i,j)-TERR(i,j)
                      IF(Z<=PBL) THEN
                        if(RMOL(i,j)<0.0) then  !unstable
                          ZSL=min(Z,0.1*PBL)
                          FM=1./(1.-gama*ZSL*RMOL(i,j))**0.25
                        else if(Z*RMOL(i,j)<=1) then !medium stable
                          ZSL=Z
                          FM=1.+beta*ZSL*RMOL(i,j)
                        else  !very stable
                          ZSL=Z
                          FM=beta+zsl*RMOL(i,j)
                        endif
                         kz(i,j)=karman*UST(i,j)*Z*((1-Z/PBL)**2.0)/FM
                         kz(i,j)=max(kzmin,kz(i,j))
                         kz(i,j)=min(kzmax,kz(i,j))
                      ELSE
                         kz(i,j)=kzmin
                      ENDIF !PBL
!                         if(i==99.and.j==38) print*,kz(i,j),PBLH(i,j),HEIGHT(i,j),TERR(i,j)
                     ENDIF !HEIGHT
                   enddo
                   enddo      
                     return
                     end subroutine
                     
end module eddyz_model
