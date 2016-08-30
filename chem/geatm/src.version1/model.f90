module model
contains
  subroutine sp2(tmpii1,tmpii2,heimii,hii,nzz,nz_mm5)

 real tmpii1(nz_mm5),tmpii2(nzz),heimii(nzz),hii(nz_mm5)

 tmpii2=0.

 do k=1,nzz
    hei=heimii(k)
    ! to find hei in which level of hii(:)
    
    kkk=2
    if(hei.le.hii(2))then
      kkk=2
      goto 200
      endif
    do kk=2,nz_mm5-1
      if(hei.le.hii(kk+1) .and. hei.gt.hii(kk))then
        kkk=kk+1
        goto 200
      endif
    enddo
 200 continue
    if(kkk==2)then
      tmpii2(k) = tmpii1(kkk)
    else
      ratio=(hei-hii(kkk-1))/(hii(kkK)-hii(kkk-1))
      tmpii2(k) = tmpii1(kkk-1)+(tmpii1(kkk)-tmpii1(kkk-1))*ratio
    endif
  enddo

 return
 end subroutine

 subroutine  get_w(myid,w,w0,u,v,dx,dy,dz,topo,landmask,lat,hh,sx,ex,sy,ey,nx,ny,k,ne,ktop)
  integer myid, sx, ex, sy, ey, k, ne, nx, ny
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: w,u,v,w0,dx,dy,dz,topo,landmask,lat,ktop
  real hh

if(k==1)then
  do j= sy,ey
  do i= sx,ex
    w(i,j)=0.
  enddo
  enddo
  return
endif

  do j= sy,ey
  do i= sx,ex
  w(i,j)=w0(i,j)-dz(i,j)*( (u(i+1,j)-u(i,j))/dx(i,j) &
                          +(v(i,j+1)-v(i,j))/dy(i,j)) &
                +dz(i,j)*( u(i,j)*(topo(i+1,j)-topo(i,j))/dx(i,j)  &
                          +v(i,j)*(topo(i,j+1)-topo(i,j))/dy(i,j)  &
                )/(hh-topo(i,j))
      !!!!!!!!!!! for nest domain boundary !!!!!!!!!!!!!!!!!!!
      !if(ne.gt.1) then
      !  if(i.eq.nx) w(i,j)=w(i-1,j)
      !  if(j.eq.ny) w(i,j)=w(i,j-1)
      !endif
      !!!!!!!!!!! test !!!!!!!!!!!!!
      !if(abs(w(i,j)).gt.1) then
      !  print *,"-----------------------------------------------"
      !  print *,"ne=",ne,"i=",i,"j=",j,"k=",k
      !  print *,'dx=',dx(i,j),'dy=',dy(i,j),"dz=",dz(i,j)
      !  print *,"u=",u(i,j),"v=",v(i,j),"w=",w(i,j)
      !  print *,"u(i+1)=",u(i+1,j),"v(j+1)=",v(i,j+1),"w0=",w0(i,j)
      !  print *,"topo=",topo(i,j),"topo(i+1)=",topo(i+1,j),"topo(j+1)=",topo(i,j+1)
      !  print *,"hh-topo(i,j)=",hh-topo(i,j)
      !  print *,"u(i+1,j)-u(i,j)=",u(i+1,j)-u(i,j)
      !  print *,"v(i,j+1)-v(i,j)=",v(i,j+1)-v(i,j)
      !  print *,"topo(i+1,j)-topo(i,j)=",topo(i+1,j)-topo(i,j)
      !  print *,"topo(i,j+1)-topo(i,j)=",topo(i,j+1)-topo(i,j)
      !  print *,"(u(i+1,j)-u(i,j))/dx(i,j)=",(u(i+1,j)-u(i,j))/dx(i,j)
      !  print *,"(v(i,j+1)-v(i,j))/dy(i,j))=",(v(i,j+1)-v(i,j))/dy(i,j)
      !  print *,"u(i,j)*(topo(i+1,j)-topo(i,j))/dx(i,j)=",u(i,j)*(topo(i+1,j)-topo(i,j))/dx(i,j)
      !  print *,"v(i,j)*(topo(i,j+1)-topo(i,j))/dy(i,j)=",v(i,j)*(topo(i,j+1)-topo(i,j))/dy(i,j)
      !  print *,"-----------------------------------------------"
      !endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!! for all grid !!!!!!!!!!!
      if(w(i,j).gt.1.0) w(i,j)=1.0
      if(w(i,j).lt.-1.0) w(i,j)=-1.0
      !!!! for polar !!!!!!!!!!!!!
      IF(ne.eq.1) THEN
        if(j.eq.1.or.j.eq.ny) then !! for 1x1 degree
           w(i,j)=0.0
        endif
      ENDIF
       !!! lat < 60S  !!! 
      if(lat(i,j).le.-60) then
        if(landmask(i,j).eq.1) then !! for land 
          if(w(i,j).gt.0.01) w(i,j)=0.01
          if(w(i,j).lt.-0.01) w(i,j)=-0.01
        endif
      endif
       !!! lat > 60N !!!
      if(lat(i,j).gt.60) then
        if(landmask(i,j).eq.0) then !! for sea 
         if(w(i,j).gt.0.02) w(i,j)=0.02
         if(w(i,j).lt.-0.02) w(i,j)=-0.02
        endif
      endif
      !!!!!!for layers above tropopause !!!!!!!!!!!!
      IF( k.gt.(ktop(i,j)+1) ) THEN
        if(w(i,j).gt.0.03) w(i,j)=0.03
        if(w(i,j).lt.-0.03) w(i,j)=-0.03
      ENDIF
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  enddo
  enddo
  return
  end subroutine

  subroutine  putkosa2(myid, c, emit, deltz, sx, ex, sy, ey ,dt)
  integer myid, sx, ex, sy, ey
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: emit,deltz,c
  real dt
  do j = sy,ey
     do i = sx,ex
      c(i,j)=c(i,j)+emit(i,j)*dt/deltz(i,j)
     enddo
  enddo
  return
  end subroutine

  subroutine putkosa( myid, c, U10,V10, UST, RHSFC, &
                      KOSASFC,emitratio,deltz, sx, ex, sy, ey ,dt)
  integer myid, sx, ex, sy, ey
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: U10,V10,UST, &
                                          RHSFC,KOSASFC,deltz,c
  real emitratio

  do j = sy,ey
      do i = sx,ex
        windspd = sqrt(U10(i,j)*U10(i,j)+V10(i,j)*V10(i,j))

      emitkosa = 0.
        rhratio = rhsfc(i,j)/100.
!      if(windspd.gt.5 .or. UST(i,j).ge.0.4)then
!************modified by zhao**********************
      if(UST(i,j).ge.0.4 .and. rhratio .le. 0.40)then
!     if(UST(i,j).ge.0.4 )then  
emitkosa=KOSASFC(i,j)*2.9E-2*ust(i,j)**2*(1-rhratio)
!    emitkosa=KOSASFC(i,j)*2.9E-2*ust(i,j)**2
!**************************************************
!      emitkosa=KOSASFC(i,j)*2.9E-2*ust(i,j)**2
      endif

      c(i,j)=c(i,j)+emitkosa*dt/deltz(i,j)*emitratio
    enddo
  enddo

  return
  end subroutine 

SUBROUTINE GETUST0(MYID, LANDUSE,SOIL,VEG,UST0,SOILT1,SOILRH,RHSFC,LONGI,SX,EX,SY,EY,NE)
INTEGER :: MYID, SX,EX,SY,EY,NE,IMODIS
REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: LANDUSE,SOIL,VEG,UST0,SOILT1,SOILRH,RHSFC,LONGI

! THE USTAR0 IS FROM THE WORK BY ZHU AND ZHAO (ACTA METEROLOGICAL SINICA,
! 2010,977-984)

DO J = SY, EY
 DO I = SX, EX

  IF(SOIL(I,J)< 5 .AND. LONGI(I,J) < 110 )   THEN
       UST0(I,J) = 0.35    ! DESERT EXCEPT HUNSHANDAKE
  ELSE IF (SOIL(I,J)< 5 .AND. LONGI(I,J) >= 110 )   THEN
       UST0(I,J) = 0.6     ! HUNSHANDAKE desert
  ELSE IF ( SOIL(I,J)==5 ) THEN 
       UST0(I,J) = 0.45    ! GOBI
  ELSE IF ( SOIL(I,J)==6 ) THEN
       UST0(I,J) = 0.7     ! HUANGTU PLATEAU
  ELSE 
       UST0 (I,J) = 0.4
  ENDIF
  !!!! chenhs set the threshold Ustar to 0.4 for all grid !!!!!!
  !!!!      based on (Yue and Wang et al., 2009, JGR)     !!!!!!
  UST0 (I,J) = 0.4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 ENDDO
ENDDO

RETURN
END SUBROUTINE

SUBROUTINE GETZ0(MYID,LANDUSE,LATITUDE,Z0,IMONTH0,SX,EX,SY,EY,NE)
INTEGER :: MYID, SX,EX,SY,EY,NE,IMODIS,IMONTH0,IMONTH
REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: Z0,LANDUSE,LATITUDE
DO J = SY, EY
 DO I = SX, EX
   !!!!!!!!!!!! by chenhs !!!!!!!!!!!!
   IF(LATITUDE(I,J).LT.0.) THEN
     IMONTH=IMONTH0+6
     if(IMONTH.GT.12) then
        IMONTH=IMONTH-12
     endif
   ELSE
     IMONTH=IMONTH0
   ENDIF
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   IF(IMONTH.GE.6.AND.IMONTH.LE.11) THEN ! SUMMER AND FALL
         select case(int(LANDUSE(I,J)))
             case(0)
                  Z0(I,J)=0 ! no data
             case(1)
                  Z0(I,J)=80. ! cm
             case(2)
                  Z0(I,J)=15.
             case(3)
                  Z0(I,J)=10.
             case(4)
                  Z0(I,J)=15.
             case(5)
                  Z0(I,J)=14.
             case(6)
                  Z0(I,J)=20.
             case(7)
                  Z0(I,J)=12.
             case(8)
                  Z0(I,J)=5.
             case(9)
                  Z0(I,J)=6.
             case(10)
                  Z0(I,J)=15.
             case(11)
                  Z0(I,J)=50.
             case(12)
                  Z0(I,J)=50.
             case(13)
                  Z0(I,J)=50.
             case(14)
                  Z0(I,J)=50.
             case(15)
                  Z0(I,J)=50.
             case(16)
                  Z0(I,J)=0.01
             case(17)
                  Z0(I,J)=20.
             case(18)
                  Z0(I,J)=40.
             case(19)
                  Z0(I,J)=1.
             case(20)
                  Z0(I,J)=10.
             case(21)
                  Z0(I,J)=30.
             case(22)
                  Z0(I,J)=15.
             case(23)
                  Z0(I,J)=10.
             case(24)
                  Z0(I,J)=0.1
             case(25)
                  Z0(I,J)=1.
             case(26)
                  Z0(I,J)=80.
             case(27)
                  Z0(I,J)=80.
             case(28)
                  Z0(I,J)=80.
             case(29)
                  Z0(I,J)=80.
             case(30)
                  Z0(I,J)=80.
             case(31)
                  Z0(I,J)=80.
             case(32)
                  Z0(I,J)=80.
             case(33)
                  Z0(I,J)=80.
          end select
   ELSE ! WINTER AND SPRING
         select case(int(LANDUSE(I,J)))
             case(0)
                  Z0(I,J)=0 ! no data
             case(1)
                  Z0(I,J)=80. ! cm
             case(2)
                  Z0(I,J)=5.
             case(3)
                  Z0(I,J)=2.
             case(4)
                  Z0(I,J)=5.
             case(5)
                  Z0(I,J)=5.
             case(6)
                  Z0(I,J)=20.
             case(7)
                  Z0(I,J)=10.
             case(8)
                  Z0(I,J)=1.
             case(9)
                  Z0(I,J)=1.
             case(10)
                  Z0(I,J)=15.
             case(11)
                  Z0(I,J)=50.
             case(12)
                  Z0(I,J)=50.
             case(13)
                  Z0(I,J)=50.
             case(14)
                  Z0(I,J)=50.
             case(15)
                  Z0(I,J)=20.
             case(16)
                  Z0(I,J)=0.01
             case(17)
                  Z0(I,J)=20.
             case(18)
                  Z0(I,J)=40.
             case(19)
                  Z0(I,J)=1.
             case(20)
                  Z0(I,J)=10.
             case(21)
                  Z0(I,J)=30.
             case(22)
                  Z0(I,J)=15.
             case(23)
                  Z0(I,J)=5.
             case(24)
                  Z0(I,J)=0.1
             case(25)
                  Z0(I,J)=1.
             case(26)
                  Z0(I,J)=15.
             case(27)
                  Z0(I,J)=1.
             case(28)
                  Z0(I,J)=80.
             case(29)
                  Z0(I,J)=80.
             case(30)
                  Z0(I,J)=80.
             case(31)
                  Z0(I,J)=80.
             case(32)
                  Z0(I,J)=80.
             case(33)
                  Z0(I,J)=80.
          end select

   ENDIF ! IMONTH
      Z0(I,J) = Z0(I,J)/100. ! cm-->m

 ENDDO ! I
ENDDO ! J

RETURN
END SUBROUTINE


SUBROUTINE DUSTHEIGHT(MYID,HEIZ,TERRAIN,LANDUSE,DUSTK,TOTALDUST, U,V,SX,EX,SY,EY,K,NE)
INTEGER MYID,SX,EX,SY,EY,NE
REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: HEIZ,TERRAIN,LANDUSE,DUSTK,TOTALDUST,U,V
REAL :: WIND
! TO CALCULATE THE ROUTINE EMISSIONS(Q) HEIGHT FACTOR FROM OBSERVATIONS BY SCICHIA D 2011(41):234-242
! IN DESERT C(=Q/V)(g3/M3) == 58.44*H**(-0.67)
! IN OTHER C = 30.17*H**(-0.48)
DO J=SY,EY
DO I=SX,EX
  IF(K ==1 )   TOTALDUST(I,J) = 0.0
  WIND = SQRT (U(I,J)*U(I,J)+ V(I,J)*V(I,J))
  IF (LANDUSE(I,J)==19) THEN
   DUSTK(I,J) = WIND*58.44 * MAX(HEIZ(I,J)-TERRAIN(I,J),10.)**(-0.67) 
  ELSE IF( LANDUSE(I,J)== 7 .OR. LANDUSE(I,J)==8.OR.LANDUSE(I,J)==9.OR.LANDUSE(I,J)==10) THEN
   DUSTK(I,J) = WIND*30.17 * MAX(HEIZ(I,J)-TERRAIN(I,J),10.)**(-0.48)
  ELSE
   DUSTK(I,J) = 1.E-20
  ENDIF

  TOTALDUST(I,J) = TOTALDUST(I,J) + DUSTK(I,J)

ENDDO
ENDDO

RETURN
END SUBROUTINE

SUBROUTINE DUSTHGTFACT(MYID,TOTALDUST,DUSTK,DUSTHGTF, SX,EX,SY,EY,NE)
INTEGER MYID,SX,EX,SY,EY,NE
REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: TOTALDUST,DUSTK,DUSTHGTF
! TO CALCULATE RELATIVE DUST EMISSION WEIGHT FACTOR DUST EMISSIONS 
!  Q = C*V V:M/S
DO J=SY,EY
DO I=SX,EX
 DUSTHGTF(I,J) =  DUSTK(I,J) / MAX( TOTALDUST(I,J), 1.E-20)
ENDDO
ENDDO

RETURN
END SUBROUTINE

SUBROUTINE PUTDUST(MYID,C,DUSTHGTF,ICE,SNOW,LANDUSE,&
        SOIL,VEG,Z0,USTWRF,UST0, SOILT1,T2,SOILRH,RHSFC,U10,V10,EMITFACT,&
        DELTZ,DUSTEMISS, SX,EX,SY,EY,DT,SFT,K,NE,IS,ITT)
INTEGER MYID,SX,EX,SY,EY,NE
REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: C,LANDUSE,ICE,SNOW,U10,V10,DELTZ,SOIL,VEG
!REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: DUSTCOMP1,DUSTCOMP2,DUSTCOMP3,DUSTCOMP4,&
!                                       DUSTCOMP5,DUSTCOMP6,DUSTCOMP7,DUSTCOMP8,&
!                                       DUSTCOMP9,DUSTCOMP10
REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: Z0,UST0,SOILT1,T2,RHSFC,SOILRH,USTWRF
REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: DUSTHGTF,DUSTEMISS ! DUSTEMISS is the emissions ug/m2/hr
REAL DT,EMITF,WIND,DTV,Gz1oz0,RB,BB,CONS1,CONS2,UST,SFT
REAL :: VEGF,SOILF ! CORRECTION FACTOR FROM VEG AND SOIL TYPE
REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: EMITFACT


DO J=SY,EY
DO I=SX,EX

EMITF = 0.0

 IF(LANDUSE(I,J).EQ.19) EMITF = 1.0
 IF(LANDUSE(I,J).EQ.7.OR.LANDUSE(I,J).EQ.10) EMITF = 0.3
 IF(LANDUSE(I,J).EQ.8.OR.LANDUSE(I,J).EQ.9) EMITF = 0.1

IF(VEG(I,J)<10.) THEN
 VEGF = 1.0
ELSE IF(VEG(I,J)<30.) THEN
 VEGF = 0.7
ELSE IF(VEG(I,J)<50) THEN
 VEGF = 0.5
ELSE
 VEGF = 0.0
ENDIF

IF(SOIL(I,J).LE.5 ) THEN
 SOILF = 1.0
ELSE IF(SOIL(I,J).EQ.6) THEN
 SOILF = 0.7
ELSE IF(SOIL(I,J).EQ.7) THEN
 SOILF = 0.85
ELSE
 SOILF = 0.0
ENDIF

EMITF =EMITF*VEGF*SOILF

EMITFACT(I,J)=EMITF

IF(SOILT1(I,J)>273.15) THEN
IF(ICE(I,J)<=0.0)THEN
IF(SNOW(I,J)<=0.0)THEN
IF(EMITF>0.0)THEN

WIND=U10(I,J)*U10(I,J)+V10(I,J)*V10(I,J)
WIND=MAX(WIND,0.1)
DTV=T2(I,J)-SOILT1(I,J)
Gz1oz0=ALOG(10.0/Z0(I,J))
!RB=9.8*DTV*10.0/(T2(I,J)*WIND)
RB=9.8*DTV*(10.0-Z0(I,J))/(T2(I,J)*WIND)  !chenhs
BB=70.0*0.4*0.4*SQRT(ABS(RB)*10.0/Z0(I,J))/(Gz1oz0*Gz1oz0)

CONS1=0.4*SQRT(WIND)/Gz1oz0

IF(RB>=0.0)THEN
  CONS2=1.0/(1.0+4.7*RB)
ELSE
  CONS2=SQRT( 1.0-9.4*RB/(1.0+BB) )
ENDIF

UST=CONS1*CONS2

! USE WRF OUTPUT UST
 UST = USTWRF(I,J)

IF(UST>UST0(I,J).AND.RHSFC(I,J)<40.)THEN
!  C(I,J)=C(I,J)+1.e-8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*&
!         (1.0-UST0(I,J)/UST)*(1.0-RHSFC(I,J)/40.)*1e09 ! kg-->ug
!   C(I,J)=C(I,J)+1.28e4*2.9e-11*DT*EMITF*UST*UST*SFT/DELTZ(I,J)*&
!          (1.0-UST0(I,J)/UST)*(1.0-RHSFC(I,J)/40.)*1e9  !chenhs atmospheric science ! kg-->ug

   C(I,J)=C(I,J)+1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*DUSTHGTF(I,J)*&
           (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)*1.E09 !kg-->ug  !chenhs doctor

   !DUSTCOMP1(I,J) = DUSTCOMP1(I,J) + 1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*DUSTHGTF(I,J)*&
   !        (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)*1.E09 *  0.07  ! CACO3

   !DUSTCOMP2(I,J) = DUSTCOMP2(I,J) + 1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*DUSTHGTF(I,J)*&
   !        (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)*1.E09 *  0.055 ! MGCO3

   !DUSTCOMP3(I,J) = DUSTCOMP3(I,J) + 1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*DUSTHGTF(I,J)*&
   !        (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)*1.E09 *  0.033 ! K2CO3

   !DUSTCOMP4(I,J) = DUSTCOMP4(I,J) + 1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*DUSTHGTF(I,J)*&
   !        (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)*1.E09 *  0.026 ! NA2CO3
   
   !DUSTCOMP5(I,J) = DUSTCOMP5(I,J) + 1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*DUSTHGTF(I,J)*&
   !        (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)*1.E09 *  0.60  ! SIO2

   !DUSTCOMP6(I,J) = DUSTCOMP6(I,J) + 1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*DUSTHGTF(I,J)*&
   !        (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)*1.E09 *  0.141 ! AL2O3

   !DUSTCOMP7(I,J) = DUSTCOMP7(I,J) + 1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*DUSTHGTF(I,J)*&
   !        (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)*1.E09 *  0.00  ! DSO4

   !DUSTCOMP8(I,J) = DUSTCOMP8(I,J) + 1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*DUSTHGTF(I,J)*&
   !        (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)*1.E09 *  0.00  ! DNO3

   !DUSTCOMP9(I,J) = DUSTCOMP9(I,J) + 1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*DUSTHGTF(I,J)*&
   !        (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)*1.E09 *  1.2E-4  ! FEII

   !DUSTCOMP10(I,J) = DUSTCOMP10(I,J) + 1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*DUSTHGTF(I,J)*&
   !        (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)*1.E09 *  0.028  ! FEIII 

   !DUSTCOMP1(I,J) = AMAX1( DUSTCOMP1(I,J), 1.E-20)
   !DUSTCOMP2(I,J) = AMAX1( DUSTCOMP2(I,J), 1.E-20)
   !DUSTCOMP3(I,J) = AMAX1( DUSTCOMP3(I,J), 1.E-20)
   !DUSTCOMP4(I,J) = AMAX1( DUSTCOMP4(I,J), 1.E-20)
   !DUSTCOMP5(I,J) = AMAX1( DUSTCOMP5(I,J), 1.E-20)
   !DUSTCOMP6(I,J) = AMAX1( DUSTCOMP6(I,J), 1.E-20)
   !DUSTCOMP7(I,J) = AMAX1( DUSTCOMP7(I,J), 1.E-20)
   !DUSTCOMP8(I,J) = AMAX1( DUSTCOMP8(I,J), 1.E-20)
   !DUSTCOMP9(I,J) = AMAX1( DUSTCOMP9(I,J), 1.E-20)
   !DUSTCOMP10(I,J) = AMAX1( DUSTCOMP10(I,J), 1.E-20)

!  THIS SUBROUTINE IS also TO ALLOCATE THE DUST COMPOSITIONS FROM THE
!  PUBLICATION Feng Yan  et al.,Global Modeling of Nitrate and Ammonium:
!  Interaction of Aerosols and Tropospheric Chemistry.
!  MASS WEIGHT CACO3: 7% MGCO3:5.5%; K2CO3 3.3%; NA2CO3 2.6% SIO2 60%; !  AL2O3
!  14.1%; Fe(III): 2.8 Fe(II) 1.2-04 ! Fe is from Zhao Phd thesis(2006)


   IF(K==1) DUSTEMISS(I,J)=  DUSTEMISS(I,J) + 1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT*&
           (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)   ! kg/m2/hr

ENDIF

ENDIF
ENDIF
ENDIF
ENDIF ! SOILT1

   DUSTEMISS(I,J) = AMAX1( DUSTEMISS(I,J), 1.E-20 )

ENDDO ! I
ENDDO ! J

RETURN

END SUBROUTINE


SUBROUTINE PUTDUSTCOM (MYID,C,DUSTCOMP1,DUSTCOMP2,DUSTCOMP3,DUSTCOMP4,DUSTCOMP5,DUSTCOMP6,&
        DUSTCOMP7,DUSTCOMP8,DUSTCOMP9,DUSTCOMP10,DUSTHGTF,ICE,SNOW,LANDUSE,&
        SOIL,VEG,Z0,USTWRF,UST0, SOILT1,T2,SOILRH,RHSFC,U10,V10,EMITFACT,&
        DELTZ,DUSTEMISS, SX,EX,SY,EY,DT,SFT,K,NE,IS,ITT)
INTEGER MYID,SX,EX,SY,EY,NE
REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: C,LANDUSE,ICE,SNOW,U10,V10,DELTZ,SOIL,VEG
REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: DUSTCOMP1,DUSTCOMP2,DUSTCOMP3,DUSTCOMP4,&
                                       DUSTCOMP5,DUSTCOMP6,DUSTCOMP7,DUSTCOMP8,&
                                       DUSTCOMP9,DUSTCOMP10
REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: Z0,UST0,SOILT1,T2,RHSFC,SOILRH,USTWRF
REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: DUSTHGTF,DUSTEMISS ! DUSTEMISS is the emissions ug/m2/hr
REAL DT,EMITF,WIND,DTV,Gz1oz0,RB,BB,CONS1,CONS2,UST,SFT
REAL :: VEGF,SOILF ! CORRECTION FACTOR FROM VEG AND SOIL TYPE
REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: EMITFACT


DO J=SY,EY
DO I=SX,EX

EMITF = 0.0

 IF(LANDUSE(I,J).EQ.19) EMITF = 1.0
 IF(LANDUSE(I,J).EQ.7.OR.LANDUSE(I,J).EQ.10) EMITF = 0.3
 IF(LANDUSE(I,J).EQ.8.OR.LANDUSE(I,J).EQ.9) EMITF = 0.1 


IF(VEG(I,J)<10.) THEN
 VEGF = 1.0
ELSE IF(VEG(I,J)<30.) THEN
 VEGF = 0.7
ELSE IF(VEG(I,J)<50) THEN
 VEGF = 0.5
ELSE
 VEGF = 0.0
ENDIF

IF(SOIL(I,J).LE.5 ) THEN
 SOILF = 1.0
ELSE IF(SOIL(I,J).EQ.6) THEN
 SOILF = 0.7
ELSE IF(SOIL(I,J).EQ.7) THEN
 SOILF = 0.85
ELSE
 SOILF = 0.0
ENDIF

EMITF =EMITF*VEGF*SOILF

EMITFACT(I,J)=EMITF

IF(SOILT1(I,J)>273.15) THEN
IF(ICE(I,J)<=0.0)THEN
IF(SNOW(I,J)<=0.0)THEN
IF(EMITF>0.0)THEN

WIND=U10(I,J)*U10(I,J)+V10(I,J)*V10(I,J)
WIND=MAX(WIND,0.1)
DTV=T2(I,J)-SOILT1(I,J)
Gz1oz0=ALOG(10.0/Z0(I,J))
!RB=9.8*DTV*10.0/(T2(I,J)*WIND)
RB=9.8*DTV*(10.0-Z0(I,J))/(T2(I,J)*WIND)  !chenhs
BB=70.0*0.4*0.4*SQRT(ABS(RB)*10.0/Z0(I,J))/(Gz1oz0*Gz1oz0)

CONS1=0.4*SQRT(WIND)/Gz1oz0

IF(RB>=0.0)THEN
  CONS2=1.0/(1.0+4.7*RB)
ELSE
  CONS2=SQRT( 1.0-9.4*RB/(1.0+BB) )
ENDIF

UST=CONS1*CONS2

! USE WRF OUTPUT UST
 UST = USTWRF(I,J)

IF(UST>UST0(I,J).AND.RHSFC(I,J)<40.)THEN
!  C(I,J)=C(I,J)+1.e-8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*&
!         (1.0-UST0(I,J)/UST)*(1.0-RHSFC(I,J)/40.)*1e09 ! kg-->ug
!   C(I,J)=C(I,J)+1.28e4*2.9e-11*DT*EMITF*UST*UST*SFT/DELTZ(I,J)*&
!          (1.0-UST0(I,J)/UST)*(1.0-RHSFC(I,J)/40.)*1e9  !chenhs atmospheric science ! kg-->ug

   C(I,J)=C(I,J)+1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*DUSTHGTF(I,J)*&
           (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)*1.E09 !kg-->ug  !chenhs doctor
   
   DUSTCOMP1(I,J) = DUSTCOMP1(I,J) + 1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*DUSTHGTF(I,J)*&
           (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)*1.E09 *  0.07  ! CACO3

   DUSTCOMP2(I,J) = DUSTCOMP2(I,J) + 1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*DUSTHGTF(I,J)*&
           (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)*1.E09 *  0.055 ! MGCO3

   DUSTCOMP3(I,J) = DUSTCOMP3(I,J) + 1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*DUSTHGTF(I,J)*&
           (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)*1.E09 *  0.033 ! K2CO3

   DUSTCOMP4(I,J) = DUSTCOMP4(I,J) + 1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*DUSTHGTF(I,J)*&
           (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)*1.E09 *  0.026 ! NA2CO3

   DUSTCOMP5(I,J) = DUSTCOMP5(I,J) + 1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*DUSTHGTF(I,J)*&
           (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)*1.E09 *  0.60  ! SIO2

   DUSTCOMP6(I,J) = DUSTCOMP6(I,J) + 1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*DUSTHGTF(I,J)*&
           (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)*1.E09 *  0.141 ! AL2O3

   DUSTCOMP7(I,J) = DUSTCOMP7(I,J) + 1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*DUSTHGTF(I,J)*&
           (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)*1.E09 *  0.00  ! DSO4

   DUSTCOMP8(I,J) = DUSTCOMP8(I,J) + 1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*DUSTHGTF(I,J)*&
           (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)*1.E09 *  0.00  ! DNO3

   DUSTCOMP9(I,J) = DUSTCOMP9(I,J) + 1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*DUSTHGTF(I,J)*&
           (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)*1.E09 *  1.2E-4  ! FEII

   DUSTCOMP10(I,J) = DUSTCOMP10(I,J) + 1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT/DELTZ(I,J)*DUSTHGTF(I,J)*&
           (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)*1.E09 *  0.028  ! FEIII 

   DUSTCOMP1(I,J) = AMAX1( DUSTCOMP1(I,J), 1.E-20)
   DUSTCOMP2(I,J) = AMAX1( DUSTCOMP2(I,J), 1.E-20)
   DUSTCOMP3(I,J) = AMAX1( DUSTCOMP3(I,J), 1.E-20)
   DUSTCOMP4(I,J) = AMAX1( DUSTCOMP4(I,J), 1.E-20)
   DUSTCOMP5(I,J) = AMAX1( DUSTCOMP5(I,J), 1.E-20)
   DUSTCOMP6(I,J) = AMAX1( DUSTCOMP6(I,J), 1.E-20)
   DUSTCOMP7(I,J) = AMAX1( DUSTCOMP7(I,J), 1.E-20)
   DUSTCOMP8(I,J) = AMAX1( DUSTCOMP8(I,J), 1.E-20)
   DUSTCOMP9(I,J) = AMAX1( DUSTCOMP9(I,J), 1.E-20)
   DUSTCOMP10(I,J) = AMAX1( DUSTCOMP10(I,J), 1.E-20)

  
!  THIS SUBROUTINE IS also TO ALLOCATE THE DUST COMPOSITIONS FROM THE
!  PUBLICATION Feng Yan  et al.,Global Modeling of Nitrate and Ammonium:
!  Interaction of Aerosols and Tropospheric Chemistry.
!  MASS WEIGHT CACO3: 7% MGCO3:5.5%; K2CO3 3.3%; NA2CO3 2.6% SIO2 60%; !  AL2O3
!  14.1%; Fe(III): 2.8 Fe(II) 1.2-04 ! Fe is from Zhao Phd thesis(2006)


   IF(K==1) DUSTEMISS(I,J)=  DUSTEMISS(I,J) + 1.0e-5*1.29/9.8*DT*EMITF*UST*UST*UST*SFT*&
           (1.0+UST0(I,J)/UST)*(1.0-(UST0(I,J)*UST0(I,J))/(UST*UST))*(1.0-RHSFC(I,J)/40.)   ! kg/m2/hr

ENDIF

ENDIF
ENDIF
ENDIF
ENDIF ! SOILT1

   DUSTEMISS(I,J) = AMAX1( DUSTEMISS(I,J), 1.E-20 )

ENDDO ! I
ENDDO ! J

RETURN

END SUBROUTINE

SUBROUTINE PUTSALT (MYID, C, R0,R1,LAND,ICE,U10,V10,RH,DELTZ,DELTX,&
                  SEAEMISS, SX,EX,SY,EY,NE,DT)
INTEGER MYID,SX,EX,SY,EY,NE
REAL :: DT
REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: C,LAND,ICE,U10,V10,DELTZ,RH,DELTX,SEAEMISS
!REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: SEACOMP1, SEACOMP2,SEACOMP3,SEACOMP4,SEACOMP5,&
!                                        SEACOMP6,SEACOMP7,SEACOMP8
REAL ::  RH_TMP,A1,A2,A3,A4, B, DENSITY, NM,X,C80,C0
REAL DR,CONST
REAL,DIMENSION(5000) :: RMID
REAL,DIMENSION(0:5000) :: REDGE
REAL :: R0,R1      ! start and end radius 
REAL,DIMENSION(5000) :: FM_CL06, FM_M03, FM, FM_DEL00
REAL :: TMP1,TMP2,TMP3,TMP4
LOGICAL :: ISURF ! 1: shore 2: open sea
INTEGER :: NR,R
!!! THE SSA EMISSIONS ARE FROM THE PAPER BY Athanasopoulou et al., acp, 2008
!  THIS SUBROUTINE IS ALSO TO ALLOCATE THE SEA SALT COMPOSITIONS FROM THE
!  PUBLICATION Athanasopoulou et al., ACP, 2008; The role of sea salt
!  emissions and heterogeneous chemistry in the air quality of polluted
!  coastal areas.
!  MASS WEIGHT CL : 55%; Na: 31%; ss-SO4 : 8%; Mg 4%; Ca:1% , K 1%

 DR=0.05
 NR=INT(((R1-R0)/DR)+0.05)-1

DO J=SY,EY !! by chenhs
DO I=SX,EX
RH_TMP = RH(I,J) /100.
RH_TMP = MIN( RH_TMP, 0.95 )  !%-->0-1
W10M=SQRT(U10(I,J)*U10(I,J)+V10(I,J)*V10(I,J))


IF(LAND(I,J)== 16)THEN

IF(ICE(I,J)<=0.0)THEN

 ISURF = .FALSE.

!IF(I.NE.SX.AND.I.NE.EX.AND.J.NE.SY.AND.J.NE.EY) THEN !! by chenhs
IF (LAND(I+1,J).NE.16 .OR. LAND(I-1,J).NE.16 .OR. LAND(I,J-1).NE.16 .OR. LAND(I,J+1).NE.16 ) THEN
  ISURF = .TRUE.
ENDIF ! ISURF
!ELSE
!  ISURF = .FALSE.
!ENDIF


 C80 = 1.82 * ( (1-RH_TMP)/ (2-RH_TMP))**0.33

 REDGE(0) = R0

IF (.NOT.ISURF) THEN

 DO R = 1, NR
   RMID(R)=REDGE(R-1)+(DR/2.0)
   REDGE(R)=REDGE(R-1)+DR

   A1 = -5.001E03+ 0.808E06*RMID(R)*C80 -1.980E07*(RMID(R)*C80)**2 + &
        2.188E8*(RMID(R)*C80)**3 - 1.144E09*(RMID(R)*C80)**4&
       + 2.290E09*(RMID(R)*C80)**5


   A2 = 3.854E03 + 1.168E04*(RMID(R)*C80)-6.572E04*(RMID(R)*C80)**2&
      + 1.003E05*(RMID(R)*C80)**3 - 6.407E04*(RMID(R)*C80)**4 + &
      1.493E04*(RMID(R)*C80)**5

   A3 = 4.498E02 + 0.839E03*(RMID(R)*C80) - 5.394E02*(RMID(R)*C80)**2 &
       + 1.218E02*(RMID(R)*C80)**3 - 1.213E01*(RMID(R)*C80)**4&
       + 4.514E-1*(RMID(R)*C80)**5

   A4 = 4.7 * (1 + 30.*RMID(R))** (-0.017*RMID(R)**(-1.44))

   B = ( 0.433 - LOG10(RMID(R)) ) / 0.433

   DENSITY =  1.E03 * ( 3.8033-16.248*RH_TMP + 46.085 * RH_TMP **2 - 68.317 * RH_TMP ** 3 &
                + 50.932*RH_TMP**4 - 15.261 * RH_TMP**5   ) / 1.E03 ! DENSITY IN G/M3

   X = 3.1657 - 19.079 * RH_TMP + 55.72 * RH_TMP**2 + 83.998 * RH_TMP**3 + &
            63.436*RH_TMP**4 - 19.248 * RH_TMP**5

   NM = 1.E-15 * 3.14 * 8. * RMID(R)**3 * DENSITY * X / 6.

   IF ( REDGE(R).LE.0.065 ) THEN  ! R<0.065UM
     FM_CL06(R) = DR/C80 * ( 1.E04*3.84*1E-06*W10M**3.41*A1*NM/ LOG(10.) / RMID (R) )
     FM_M03(R) = 0.0
     FM_DEL00(R) = 0.0
!   IF(I==65.AND.J==43) PRINT*, R,W10M,  FM_CL06(R), '111'

   ELSE IF ( REDGE(R).LE.0.6  ) THEN ! R<0.65UM
     TMP1 = (W10M **3.41)/( LOG(10.* RMID(R)) )
     TMP2 = A2*NM
     TMP3 = 1.E04*3.84*1.E-06
     FM_CL06(R) = TMP1 * TMP2 * TMP3 * (DR/C80)
!     FM_CL06(R) = DR/C80 * ( 1.E04*3.84*1E-06*W10M**3.41*A2*NM/ LOG(10.) / RMID (R) )
     FM_M03(R) = 0.0
     FM_DEL00(R) = 0.0
!   IF(I==65.AND.J==43) PRINT*, R, TMP1, TMP2, TMP3, FM_CL06(R), '222'
   ELSE IF ( REDGE(R).LE.4  ) THEN ! R < 4UM
     FM_CL06(R) = DR/C80 * ( 1.E04*3.84*1E-06*W10M**3.41*A3*NM/ LOG(10.) / RMID (R) )
     FM_M03(R) = 0.0
     FM_DEL00(R) = 0.0
!     IF(I==65.AND.J==43) PRINT*, R, W10M, FM_CL06(R), '333'
   ELSE IF ( REDGE(R).LE.10.)  THEN
     TMP1 = 1 + 0.057 * RMID(R)**3.45
     TMP2 = 1.607* EXP(-B**2)
     TMP3 = 10.** TMP2
     TMP4 = RMID(R) ** ( -A4)
     FM_CL06(R) = 0.0
     FM_M03(R) =  DR/C80 * C80 * NM * 1.373 * W10M ** 3.41 * TMP4 * TMP1 * TMP3
     FM_DEL00(R) = 0.0
!     IF(I==65.AND.J==43) PRINT*, R,W10M,FM_M03(R), '444'
   ENDIF

   FM(R) = FM_CL06(R) + FM_M03(R) + FM_DEL00(R)  ! G/M2/S

 ENDDO ! R


ELSE IF (ISURF) THEN ! OFFSHORE

 DO R = 1, NR
   RMID(R)=REDGE(R-1)+(DR/2.0)
   REDGE(R)=REDGE(R-1)+DR

   TMP1 = ((1.-RH_TMP)/(2.-RH_TMP))**0.33
   TMP2 = (35./38.5)**0.33
   TMP3 = 1 + 2.5 * 10. - 4.*(38.5-35)

   C0 = 3.7 * TMP1 * TMP2 * TMP3


   DENSITY =  1.E03 * ( 3.8033-16.248*RH_TMP + 46.085 * RH_TMP **2 - 68.317 * RH_TMP ** 3 &
                + 50.932*RH_TMP**4 - 15.261 * RH_TMP**5   ) /1.E03 ! DENSITY IN G/M3

   X = 3.1657 - 19.079 * RH_TMP + 55.72 * RH_TMP**2 + 83.998 * RH_TMP**3 + &
            63.436*RH_TMP**4 - 19.248 * RH_TMP**5

   NM = 2. * 1.E-15 * 3.14 * 8. * RMID(R)**3 * DENSITY * X /6.

   !FM_DEL00(R) = DR/C0 * NM * 1.1E7 * EXP (0.23 * W10M) * ( 2.* RMID(R) * C0)**(-1.65)*50./DELTX(I,J) ! G/M2/S 50m surf zone
   FM_DEL00(R) = DR/C0 * NM * 1.1E7 * EXP (0.23 * W10M) * ( 2.* RMID(R) * C0)**(-1.65)*0.0025 ! assume 0.25% surf-zone, by chenhs

   A1 = -5.001E03+ 0.808E06*RMID(R)*C80 -1.980E07*(RMID(R)*C80)**2 + &
        2.188E8*(RMID(R)*C80)**3 - 1.144E09*(RMID(R)*C80)**4&
       + 2.290E09*(RMID(R)*C80)**5


   A2 = 3.854E03 + 1.168E04*(RMID(R)*C80)-6.572E04*(RMID(R)*C80)**2&
      + 1.003E05*(RMID(R)*C80)**3 - 6.407E04*(RMID(R)*C80)**4 + &
      1.493E04*(RMID(R)*C80)**5

   A3 = 4.498E02 + 0.839E03*(RMID(R)*C80) - 5.394E02*(RMID(R)*C80)**2 &
       + 1.218E02*(RMID(R)*C80)**3 - 1.213E01*(RMID(R)*C80)**4&
       + 4.514E-1*(RMID(R)*C80)**5

   A4 = 4.7 * (1 + 30.*RMID(R))** (-0.017*RMID(R)**(-1.44))
  
   B = ( 0.433 - LOG10(RMID(R)) ) / 0.433 ! added by juanxiong he

   DENSITY =  1.E03 * ( 3.8033-16.248*RH_TMP + 46.085 * RH_TMP **2 - 68.317 * RH_TMP ** 3 &
                + 50.932*RH_TMP**4 - 15.261 * RH_TMP**5   ) / 1.E03 ! DENSITY IN G/M3

   X = 3.1657 - 19.079 * RH_TMP + 55.72 * RH_TMP**2 + 83.998 * RH_TMP**3 + &
            63.436*RH_TMP**4 - 19.248 * RH_TMP**5

   NM = 1.E-15 * 3.14 * 8. * RMID(R)**3 * DENSITY * X / 6.
   
   IF ( REDGE(R).LE.0.065 ) THEN  ! R<0.065UM
     !FM_CL06(R) = DR/C80 * ( 1.E04*3.84*1E-06*W10M**3.41*A1*NM/ LOG(10.) / RMID(R) ) * ( DELTX(I,J) - 50.) /DELTX(I,J) 
     FM_CL06(R) = DR/C80 * ( 1.E04*3.84*1E-06*W10M**3.41*A1*NM/ LOG(10.) / RMID(R) ) * (1-0.0025) ! assume 0.25% surf-zone, by chenhs
     FM_M03(R) = 0.0
!   IF(I==65.AND.J==43) PRINT*, R,W10M,  FM_CL06(R), '111'

   ELSE IF ( REDGE(R).LE.0.6  ) THEN ! R<0.65UM
     TMP1 = (W10M **3.41)/( LOG(10.* RMID(R)) )
     TMP2 = A2*NM
     TMP3 = 1.E04*3.84*1.E-06
     !FM_CL06(R) = TMP1 * TMP2 * TMP3 * (DR/C80)*( DELTX(I,J) - 50.) /DELTX(I,J)
     FM_CL06(R) = TMP1 * TMP2 * TMP3 * (DR/C80)*(1-0.0025) ! assume 0.25% surf-zone, by chenhs
     FM_M03(R) = 0.0
!   IF(I==65.AND.J==43) PRINT*, R, TMP1, TMP2, TMP3, FM_CL06(R), '222'
   ELSE IF ( REDGE(R).LE.4  ) THEN ! R < 4UM
     !FM_CL06(R) = DR/C80 * ( 1.E04*3.84*1E-06*W10M**3.41*A3*NM/ LOG(10.) / RMID(R) ) * ( DELTX(I,J) - 50.) /DELTX(I,J)
     FM_CL06(R) = DR/C80 * ( 1.E04*3.84*1E-06*W10M**3.41*A3*NM/ LOG(10.) / RMID(R) ) * (1-0.0025) ! assume 0.25% surf-zone, by chenhs
     FM_M03(R) = 0.0
!     IF(I==65.AND.J==43) PRINT*, R, W10M, FM_CL06(R), '333'
   ELSE IF ( REDGE(R).LE.10.)  THEN
     TMP1 = 1 + 0.057 * RMID(R)**3.45
     TMP2 = 1.607* EXP(-B**2)
     TMP3 = 10.** TMP2
     TMP4 = RMID(R) ** ( -A4)
     FM_CL06(R) = 0.0
     !FM_M03(R) =  DR/C80 * C80 * NM * 1.373 * W10M ** 3.41 * TMP4 * TMP1 * TMP3 * ( DELTX(I,J) - 50.) /DELTX(I,J)
     FM_M03(R) =  DR/C80 * C80 * NM * 1.373 * W10M ** 3.41 * TMP4 * TMP1 * TMP3 * (1-0.0025) ! assume 0.25% surf-zone, by chenhs
!     IF(I==65.AND.J==43) PRINT*, R,W10M,FM_M03(R), '444'
   ENDIF


   FM(R) = FM_CL06(R) + FM_M03(R) + FM_DEL00(R)  ! G/M2/S


 ENDDO ! R
ENDIF ! ISURF



  DO R = 1, NR
   C (I,J) = C(I,J) + DT * FM(R) / DELTZ(I,J) * 1E06 ! ug/m3 
   !SEACOMP1(I,J) = SEACOMP1(I,J ) +  DT * FM(R) / DELTZ(I,J) * 1E06 *  0.55 ! CL
   !SEACOMP2(I,J) = SEACOMP2(I,J ) +  DT * FM(R) / DELTZ(I,J) * 1E06 *  0.31 ! NA
   !SEACOMP3(I,J) = SEACOMP3(I,J ) +  DT * FM(R) / DELTZ(I,J) * 1E06 *  0.08 ! SS-SO4
   !SEACOMP4(I,J) = SEACOMP4(I,J ) +  DT * FM(R) / DELTZ(I,J) * 1E06 *  0.04 ! Mg
   !SEACOMP5(I,J) = SEACOMP5(I,J ) +  DT * FM(R) / DELTZ(I,J) * 1E06 *  0.01 ! Ca
   !SEACOMP6(I,J) = SEACOMP6(I,J ) +  DT * FM(R) / DELTZ(I,J) * 1E06 *  0.01 ! K
   !SEACOMP7(I,J) = SEACOMP7(I,J ) +  DT * FM(R) / DELTZ(I,J) * 1E06 *  0.0  ! DSO4
   !SEACOMP8(I,J) = SEACOMP8(I,J ) +  DT * FM(R) / DELTZ(I,J) * 1E06 *  0.0  ! DNO3

   SEAEMISS (I, J) = SEAEMISS(I,J) + FM(R)*DT/1.E03  ! kg/m2/hr
  ENDDO ! R

    !SEACOMP1(I,J) = AMAX1 (SEACOMP1(I,J), 1.E-20 )
    !SEACOMP2(I,J) = AMAX1 (SEACOMP2(I,J), 1.E-20 )
    !SEACOMP3(I,J) = AMAX1 (SEACOMP3(I,J), 1.E-20 )
    !SEACOMP4(I,J) = AMAX1 (SEACOMP4(I,J), 1.E-20 )
    !SEACOMP5(I,J) = AMAX1 (SEACOMP5(I,J), 1.E-20 )
    !SEACOMP6(I,J) = AMAX1 (SEACOMP6(I,J), 1.E-20 )
    !SEACOMP7(I,J) = AMAX1 (SEACOMP7(I,J), 1.E-20 )
    !SEACOMP8(I,J) = AMAX1 (SEACOMP8(I,J), 1.E-20 )



ENDIF ! ICE  

ENDIF ! LAND

   C (I,J) = AMAX1( C(I,J), 1.E-20)
ENDDO ! I
ENDDO ! J


RETURN
END SUBROUTINE


SUBROUTINE PUTSALTCOM (MYID, C,SEACOMP1, SEACOMP2,SEACOMP3,SEACOMP4,SEACOMP5,&
                    SEACOMP6,SEACOMP7,SEACOMP8, R0,R1,LAND,ICE,U10,V10,RH,DELTZ,DELTX,&
                  SEAEMISS, SX,EX,SY,EY,NE,DT)
INTEGER MYID,SX,EX,SY,EY,NE
REAL :: DT
REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: C,LAND,ICE,U10,V10,DELTZ,RH,DELTX,SEAEMISS
REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: SEACOMP1, SEACOMP2,SEACOMP3,SEACOMP4,SEACOMP5,&
                                        SEACOMP6,SEACOMP7,SEACOMP8
REAL ::  RH_TMP,A1,A2,A3,A4, B, DENSITY, NM,X,C80,C0
REAL DR,CONST
REAL,DIMENSION(5000) :: RMID
REAL,DIMENSION(0:5000) :: REDGE
REAL :: R0,R1      ! start and end radius 
REAL,DIMENSION(5000) :: FM_CL06, FM_M03, FM, FM_DEL00
REAL :: TMP1,TMP2,TMP3,TMP4
LOGICAL :: ISURF ! 1: shore 2: open sea
INTEGER :: NR,R
!!! THE SSA EMISSIONS ARE FROM THE PAPER BY Athanasopoulou et al., acp, 2008
!  THIS SUBROUTINE IS ALSO TO ALLOCATE THE SEA SALT COMPOSITIONS FROM THE
!  PUBLICATION Athanasopoulou et al., ACP, 2008; The role of sea salt
!  emissions and heterogeneous chemistry in the air quality of polluted
!  coastal areas.
!  MASS WEIGHT CL : 55%; Na: 31%; ss-SO4 : 8%; Mg 4%; Ca:1% , K 1%

 DR=0.05
 NR=INT(((R1-R0)/DR)+0.05)-1
  
DO J=SY,EY !! by chenhs
DO I=SX,EX
RH_TMP = RH(I,J) /100.
RH_TMP = MIN( RH_TMP, 0.95 )  !%-->0-1
W10M=SQRT(U10(I,J)*U10(I,J)+V10(I,J)*V10(I,J))


IF(LAND(I,J)== 16)THEN

IF(ICE(I,J)<=0.0)THEN

 ISURF = .FALSE.

!IF(I.NE.SX.AND.I.NE.EX.AND.J.NE.SY.AND.J.NE.EX) THEN  !! by chenhs
IF (LAND(I+1,J).NE.16 .OR. LAND(I-1,J).NE.16 .OR. LAND(I,J-1).NE.16 .OR. LAND(I,J+1).NE.16 ) THEN
  ISURF = .TRUE.
ENDIF ! ISURF
!ELSE 
!  ISURF = .FALSE. 
!ENDIF 


 C80 = 1.82 * ( (1-RH_TMP)/ (2-RH_TMP))**0.33

 REDGE(0) = R0

IF (.NOT.ISURF) THEN

 DO R = 1, NR
   RMID(R)=REDGE(R-1)+(DR/2.0)
   REDGE(R)=REDGE(R-1)+DR

   A1 = -5.001E03+ 0.808E06*RMID(R)*C80 -1.980E07*(RMID(R)*C80)**2 + &
        2.188E8*(RMID(R)*C80)**3 - 1.144E09*(RMID(R)*C80)**4&
       + 2.290E09*(RMID(R)*C80)**5


   A2 = 3.854E03 + 1.168E04*(RMID(R)*C80)-6.572E04*(RMID(R)*C80)**2&
      + 1.003E05*(RMID(R)*C80)**3 - 6.407E04*(RMID(R)*C80)**4 + &
      1.493E04*(RMID(R)*C80)**5   
    
   A3 = 4.498E02 + 0.839E03*(RMID(R)*C80) - 5.394E02*(RMID(R)*C80)**2 &
       + 1.218E02*(RMID(R)*C80)**3 - 1.213E01*(RMID(R)*C80)**4&
       + 4.514E-1*(RMID(R)*C80)**5

   A4 = 4.7 * (1 + 30.*RMID(R))** (-0.017*RMID(R)**(-1.44))
 
   B = ( 0.433 - LOG10(RMID(R)) ) / 0.433
 
   DENSITY =  1.E03 * ( 3.8033-16.248*RH_TMP + 46.085 * RH_TMP **2 - 68.317 * RH_TMP ** 3 &
                + 50.932*RH_TMP**4 - 15.261 * RH_TMP**5   ) / 1.E03 ! DENSITY IN G/M3

   X = 3.1657 - 19.079 * RH_TMP + 55.72 * RH_TMP**2 + 83.998 * RH_TMP**3 + &
            63.436*RH_TMP**4 - 19.248 * RH_TMP**5

   NM = 1.E-15 * 3.14 * 8. * RMID(R)**3 * DENSITY * X / 6.   
  

   IF ( REDGE(R).LE.0.065 ) THEN  ! R<0.065UM
     FM_CL06(R) = DR/C80 * ( 1.E04*3.84*1E-06*W10M**3.41*A1*NM/ LOG(10.) / RMID (R) ) 
     FM_M03(R) = 0.0
     FM_DEL00(R) = 0.0
!   IF(I==65.AND.J==43) PRINT*, R,W10M,  FM_CL06(R), '111'

   ELSE IF ( REDGE(R).LE.0.6  ) THEN ! R<0.65UM
     TMP1 = (W10M **3.41)/( LOG(10.* RMID(R)) )
     TMP2 = A2*NM
     TMP3 = 1.E04*3.84*1.E-06
     FM_CL06(R) = TMP1 * TMP2 * TMP3 * (DR/C80)
!     FM_CL06(R) = DR/C80 * ( 1.E04*3.84*1E-06*W10M**3.41*A2*NM/ LOG(10.) / RMID (R) )
     FM_M03(R) = 0.0
     FM_DEL00(R) = 0.0
!   IF(I==65.AND.J==43) PRINT*, R, TMP1, TMP2, TMP3, FM_CL06(R), '222'
   ELSE IF ( REDGE(R).LE.4  ) THEN ! R < 4UM
     FM_CL06(R) = DR/C80 * ( 1.E04*3.84*1E-06*W10M**3.41*A3*NM/ LOG(10.) / RMID (R) )
     FM_M03(R) = 0.0
     FM_DEL00(R) = 0.0
!     IF(I==65.AND.J==43) PRINT*, R, W10M, FM_CL06(R), '333'
   ELSE IF ( REDGE(R).LE.10.)  THEN
     TMP1 = 1 + 0.057 * RMID(R)**3.45
     TMP2 = 1.607* EXP(-B**2)
     TMP3 = 10.** TMP2
     TMP4 = RMID(R) ** ( -A4)
     FM_CL06(R) = 0.0 
     FM_M03(R) =  DR/C80 * C80 * NM * 1.373 * W10M ** 3.41 * TMP4 * TMP1 * TMP3
     FM_DEL00(R) = 0.0
!     IF(I==65.AND.J==43) PRINT*, R,W10M,FM_M03(R), '444'
   ENDIF      
 
      
   FM(R) = FM_CL06(R) + FM_M03(R) + FM_DEL00(R)  ! G/M2/S

 ENDDO ! R


ELSE IF (ISURF) THEN ! OFFSHORE

 DO R = 1, NR
   RMID(R)=REDGE(R-1)+(DR/2.0)
   REDGE(R)=REDGE(R-1)+DR

   TMP1 = ((1.-RH_TMP)/(2.-RH_TMP))**0.33
   TMP2 = (35./38.5)**0.33
   TMP3 = 1 + 2.5 * 10. - 4.*(38.5-35) 

   C0 = 3.7 * TMP1 * TMP2 * TMP3


   DENSITY =  1.E03 * ( 3.8033-16.248*RH_TMP + 46.085 * RH_TMP **2 - 68.317 * RH_TMP ** 3 &
                + 50.932*RH_TMP**4 - 15.261 * RH_TMP**5   ) /1.E03 ! DENSITY IN G/M3

   X = 3.1657 - 19.079 * RH_TMP + 55.72 * RH_TMP**2 + 83.998 * RH_TMP**3 + &
            63.436*RH_TMP**4 - 19.248 * RH_TMP**5

   NM = 2. * 1.E-15 * 3.14 * 8. * RMID(R)**3 * DENSITY * X /6.

   !FM_DEL00(R) = DR/C0 * NM * 1.1E7 * EXP (0.23 * W10M) * ( 2.* RMID(R) * C0)**(-1.65)*50./DELTX(I,J) ! G/M2/S 50m surf zone
   FM_DEL00(R) = DR/C0 * NM * 1.1E7 * EXP (0.23 * W10M) * ( 2.* RMID(R) * C0)**(-1.65)*0.0025 ! assume 0.25% surf-zone, by chenhs


   A1 = -5.001E03+ 0.808E06*RMID(R)*C80 -1.980E07*(RMID(R)*C80)**2 + &
        2.188E8*(RMID(R)*C80)**3 - 1.144E09*(RMID(R)*C80)**4&
       + 2.290E09*(RMID(R)*C80)**5


   A2 = 3.854E03 + 1.168E04*(RMID(R)*C80)-6.572E04*(RMID(R)*C80)**2&
      + 1.003E05*(RMID(R)*C80)**3 - 6.407E04*(RMID(R)*C80)**4 + &
      1.493E04*(RMID(R)*C80)**5

   A3 = 4.498E02 + 0.839E03*(RMID(R)*C80) - 5.394E02*(RMID(R)*C80)**2 &
       + 1.218E02*(RMID(R)*C80)**3 - 1.213E01*(RMID(R)*C80)**4&
       + 4.514E-1*(RMID(R)*C80)**5

   A4 = 4.7 * (1 + 30.*RMID(R))** (-0.017*RMID(R)**(-1.44))

   DENSITY =  1.E03 * ( 3.8033-16.248*RH_TMP + 46.085 * RH_TMP **2 - 68.317 * RH_TMP ** 3 &
                + 50.932*RH_TMP**4 - 15.261 * RH_TMP**5   ) / 1.E03 ! DENSITY IN G/M3

   X = 3.1657 - 19.079 * RH_TMP + 55.72 * RH_TMP**2 + 83.998 * RH_TMP**3 + &
            63.436*RH_TMP**4 - 19.248 * RH_TMP**5

   NM = 1.E-15 * 3.14 * 8. * RMID(R)**3 * DENSITY * X / 6.



   IF ( REDGE(R).LE.0.065 ) THEN  ! R<0.065UM
     !FM_CL06(R) = DR/C80 * ( 1.E04*3.84*1E-06*W10M**3.41*A1*NM/ LOG(10.) / RMID(R) ) * ( DELTX(I,J) - 50.) /DELTX(I,J) 
     FM_CL06(R) = DR/C80 * ( 1.E04*3.84*1E-06*W10M**3.41*A1*NM/ LOG(10.) / RMID(R) ) * (1-0.0025) ! assume 0.25% surf-zone, by chenhs
     FM_M03(R) = 0.0
!   IF(I==65.AND.J==43) PRINT*, R,W10M,  FM_CL06(R), '111'

   ELSE IF ( REDGE(R).LE.0.6  ) THEN ! R<0.65UM
     TMP1 = (W10M **3.41)/( LOG(10.* RMID(R)) )
     TMP2 = A2*NM
     TMP3 = 1.E04*3.84*1.E-06
     !FM_CL06(R) = TMP1 * TMP2 * TMP3 * (DR/C80)*( DELTX(I,J) - 50.) /DELTX(I,J)
     FM_CL06(R) = TMP1 * TMP2 * TMP3 * (DR/C80)*(1-0.0025) ! assume 0.25% surf-zone, by chenhs
     FM_M03(R) = 0.0
!   IF(I==65.AND.J==43) PRINT*, R, TMP1, TMP2, TMP3, FM_CL06(R), '222'
   ELSE IF ( REDGE(R).LE.4  ) THEN ! R < 4UM
     !FM_CL06(R) = DR/C80 * ( 1.E04*3.84*1E-06*W10M**3.41*A3*NM/ LOG(10.) / RMID(R) ) * ( DELTX(I,J) - 50.) /DELTX(I,J)
     FM_CL06(R) = DR/C80 * ( 1.E04*3.84*1E-06*W10M**3.41*A3*NM/ LOG(10.) / RMID(R) ) * (1-0.0025) ! assume 0.25% surf-zone, by chenhs
     FM_M03(R) = 0.0
!     IF(I==65.AND.J==43) PRINT*, R, W10M, FM_CL06(R), '333'
   ELSE IF ( REDGE(R).LE.10.)  THEN
     TMP1 = 1 + 0.057 * RMID(R)**3.45
     TMP2 = 1.607* EXP(-B**2)
     TMP3 = 10.** TMP2
     TMP4 = RMID(R) ** ( -A4)
     FM_CL06(R) = 0.0
     !FM_M03(R) =  DR/C80 * C80 * NM * 1.373 * W10M ** 3.41 * TMP4 * TMP1 * TMP3 * ( DELTX(I,J) - 50.) /DELTX(I,J)
     FM_M03(R) =  DR/C80 * C80 * NM * 1.373 * W10M ** 3.41 * TMP4 * TMP1 * TMP3 * (1-0.0025) ! assume 0.25% surf-zone, by chenhs
!     IF(I==65.AND.J==43) PRINT*, R,W10M,FM_M03(R), '444'
   ENDIF


   FM(R) = FM_CL06(R) + FM_M03(R) + FM_DEL00(R)  ! G/M2/S


 ENDDO ! R

ENDIF ! ISURF

   
        
  DO R = 1, NR
   C (I,J) = C(I,J) + DT * FM(R) / DELTZ(I,J) * 1E06 ! ug/m3 
   SEACOMP1(I,J) = SEACOMP1(I,J ) +  DT * FM(R) / DELTZ(I,J) * 1E06 *  0.55 ! CL
   SEACOMP2(I,J) = SEACOMP2(I,J ) +  DT * FM(R) / DELTZ(I,J) * 1E06 *  0.31 ! NA
   SEACOMP3(I,J) = SEACOMP3(I,J ) +  DT * FM(R) / DELTZ(I,J) * 1E06 *  0.08 ! SS-SO4
   SEACOMP4(I,J) = SEACOMP4(I,J ) +  DT * FM(R) / DELTZ(I,J) * 1E06 *  0.04 ! Mg
   SEACOMP5(I,J) = SEACOMP5(I,J ) +  DT * FM(R) / DELTZ(I,J) * 1E06 *  0.01 ! Ca
   SEACOMP6(I,J) = SEACOMP6(I,J ) +  DT * FM(R) / DELTZ(I,J) * 1E06 *  0.01 ! K
   SEACOMP7(I,J) = SEACOMP7(I,J ) +  DT * FM(R) / DELTZ(I,J) * 1E06 *  0.0  ! DSO4
   SEACOMP8(I,J) = SEACOMP8(I,J ) +  DT * FM(R) / DELTZ(I,J) * 1E06 *  0.0  ! DNO3
   
   SEAEMISS (I, J) = SEAEMISS(I,J) + FM(R)*DT/1.E03  ! kg/m2/hr
  ENDDO ! R

    SEACOMP1(I,J) = AMAX1 (SEACOMP1(I,J), 1.E-20 )
    SEACOMP2(I,J) = AMAX1 (SEACOMP2(I,J), 1.E-20 )
    SEACOMP3(I,J) = AMAX1 (SEACOMP3(I,J), 1.E-20 )
    SEACOMP4(I,J) = AMAX1 (SEACOMP4(I,J), 1.E-20 )
    SEACOMP5(I,J) = AMAX1 (SEACOMP5(I,J), 1.E-20 )
    SEACOMP6(I,J) = AMAX1 (SEACOMP6(I,J), 1.E-20 )
    SEACOMP7(I,J) = AMAX1 (SEACOMP7(I,J), 1.E-20 )
    SEACOMP8(I,J) = AMAX1 (SEACOMP8(I,J), 1.E-20 )
 

!      IF(I==65.AND.J==43) PRINT*,R0,R1,FFM!,SEAEMISS(I,J),C(I,J)

ENDIF ! ICE  

ENDIF ! LAND
   

  C (I,J) = AMAX1( C(I,J), 1.E-20)
ENDDO ! I
ENDDO ! J


RETURN
END SUBROUTINE

 subroutine putemit( myid, c, emit, deltz, sx, ex, sy, ey, dt, ratio,k,nzz,ne,ig,igas)
 integer myid, sx, ex, sy, ey,ne
 real, dimension(sx-1:ex+1,sy-1:ey+1)     :: c,emit,deltz,emit_tmp
 real, dimension(sx-1:ex+1,sy-1:ey+1,nzz,igas) :: ratio

   do j = sy,ey
      do i = sx,ex
        c(i,j)=c(i,j)+emit(i,j)*dt/deltz(i,j)*ratio(i,j,k,ig)
       enddo
   enddo
 return
 end subroutine
 
 subroutine get_ratio_emit( myid, nzz, ig, land_use, latitcrs, ratioemitAnt,ratioemitBB,ratioemitLig)
 integer myid, nzz, ig
 real land_use, latitcrs
 real, dimension(nzz) :: ratioemitAnt,ratioemitBB,ratioemitLig

       if(ig==18) THEN ! SO2
           ratioemitAnt (1)= 0.7
           ratioemitAnt (2)= 0.2
           ratioemitAnt (3)= 0.1
       else
           ratioemitAnt (1)= 0.8
           ratioemitAnt (2)= 0.15
           ratioemitAnt (3)= 0.05
       endif
        
       ratioemitBB(1)= 1.0
       if(LAND_USE <=15..and.LAND_USE >=11.) then
       ratioemitBB(1)= 0.05       !forest 
       ratioemitBB(2)= 0.025
       ratioemitBB(3)= 0.025
       ratioemitBB(4)= 0.025
       ratioemitBB(5)= 0.025
       ratioemitBB(6)= 0.025
       ratioemitBB(7)= 0.025
       ratioemitBB(8)= 0.1
       ratioemitBB(9)= 0.1
       ratioemitBB(10)= 0.15
       ratioemitBB(11)= 0.1
       ratioemitBB(12)= 0.25
       ratioemitBB(13)= 0.1
       ratioemitBB(14)= 0.0
       ratioemitBB(15)= 0.0
       ratioemitBB(16)= 0.0
       ratioemitBB(17)= 0.0
       ratioemitBB(18)= 0.0
       ratioemitBB(19)= 0.0
       ratioemitBB(20)= 0.0
     else
       ratioemitBB(1)= 0.3
       ratioemitBB(2)= 0.25
       ratioemitBB(3)= 0.15
       ratioemitBB(4)= 0.1
       ratioemitBB(5)= 0.05
       ratioemitBB(6)= 0.05
       ratioemitBB(7)= 0.05
       ratioemitBB(8)= 0.05
       ratioemitBB(9)= 0.0
       ratioemitBB(10)= 0.0
       ratioemitBB(11)= 0.0
       ratioemitBB(12)= 0.0
       ratioemitBB(13)= 0.0
       ratioemitBB(14)= 0.0
       ratioemitBB(15)= 0.0
       ratioemitBB(16)= 0.0
       ratioemitBB(17)= 0.0
       ratioemitBB(18)= 0.0
       ratioemitBB(19)= 0.0
       ratioemitBB(20)= 0.0
     endif
     if(abs(latitcrs)<30. ) then !lighting
       ratioemitLig(1)= 0.012             !TROPICAL 
       ratioemitLig(2)= 0.012
       ratioemitLig(3)= 0.012
       ratioemitLig(4)= 0.012
       ratioemitLig(5)= 0.012
       ratioemitLig(6)= 0.012
       ratioemitLig(7)= 0.012
       ratioemitLig(8)= 0.0095
       ratioemitLig(9)= 0.0095
       ratioemitLig(10)= 0.011
       ratioemitLig(11)= 0.011
       ratioemitLig(12)= 0.016
       ratioemitLig(13)= 0.011
       ratioemitLig(14)= 0.016
       ratioemitLig(15)= 0.046
       ratioemitLig(16)= 0.088
       ratioemitLig(17)= 0.172
       ratioemitLig(18)= 0.321
       ratioemitLig(19)= 0.205
     else
       ratioemitLig(1)= 0.029             !MIDLALTITUDE
       ratioemitLig(2)= 0.029
       ratioemitLig(3)= 0.029
       ratioemitLig(4)= 0.029
       ratioemitLig(5)= 0.029
       ratioemitLig(6)= 0.029
       ratioemitLig(7)= 0.029
       ratioemitLig(8)= 0.015
       ratioemitLig(9)= 0.015
       ratioemitLig(10)= 0.004
       ratioemitLig(11)= 0.004
       ratioemitLig(12)= 0.015
       ratioemitLig(13)= 0.034
       ratioemitLig(14)= 0.053
       ratioemitLig(15)= 0.074
       ratioemitLig(16)= 0.12
       ratioemitLig(17)= 0.300
       ratioemitLig(18)= 0.163
      endif
      !!!! chenhs to cut off 50% LIGNOx !!!!! 
      ratioemitLig=ratioemitLig*0.5
  return
 end subroutine
 
 subroutine wet_dep_gas( myid, c, raincon, rainnon, sx, ex, sy, ey ,igas,ig,dt)
 integer myid, sx, ex, sy, ey
 real, dimension(sx-1:ex+1,sy-1:ey+1) :: c, raincon, rainnon
 REAL, DIMENSION(88)                  :: A,B
 REAL                                 :: wetvel
 DATA A/-9.9E-09,  8.0E-04,-9.9E-09,-9.9E-09, 8.0E-06, 1.0E-05,&
        -9.9E-09, -9.9E-09,-9.9E-09,-9.9E-09,-9.9E-09,-9.9E-09,&
        -9.9E-09, -9.9E-09,-9.9E-09, 1.0E-04,-9.9E-09,-9.9E-09,&
        -9.9E-09, -9.9E-09,-9.9E-09,-9.9E-09,-9.9E-09,-9.9E-09,&
        -9.9E-09, -9.9E-09,-9.9E-09,-9.9E-09,-9.9E-09,-9.9E-09,&
        -9.9E-09, -9.9E-09,-9.9E-09,-9.9E-09,-9.9E-09,-9.9E-09,&
        -9.9E-09, -9.9E-09,-9.9E-09,-9.9E-09,-9.9E-09,-9.9E-09,&
        -9.9E-09, -9.9E-09,-9.9E-09,-9.9E-09,-9.9E-09,-9.9E-09,&
        -9.9E-09, -9.9E-09,-9.9E-09,-9.9E-09,-9.9E-09,-9.9E-09,&
        -9.9E-09, -9.9E-09,-9.9E-09,-9.9E-09,-9.9E-09,-9.9E-09,&
        -9.9E-09, -9.9E-09,-9.9E-09,-9.9E-09,-9.9E-09,-9.9E-09,&
        -9.9E-09, -9.9E-09,-9.9E-09,-9.9E-09,-9.9E-09,-9.9E-09,&
        -9.9E-09, -9.9E-09,-9.9E-09,-9.9E-09,-9.9E-09,-9.9E-09,&
        -9.9E-09,  1.0E-04, 1.0E-04, 1.0E-04, 1.0E-04,-9.9E-09,&
        -9.9E-09, -9.9E-09,-9.9E-09,-9.9E-09/
 DATA B/ 0.00   ,  0.62   , 0.00   , 0.00   , 0.62   , 0.62   ,&
         0.00   ,  0.00   , 0.00   , 0.00   , 0.00   , 0.00   ,&
         0.00   ,  0.00   , 0.00   , 0.62   , 0.00   , 0.62   ,&
         0.00   ,  0.00   , 0.00   , 0.00   , 0.00   , 0.00   ,&
         0.00   ,  0.00   , 0.00   , 0.00   , 0.00   , 0.00   ,&
         0.00   ,  0.00   , 0.00   , 0.00   , 0.00   , 0.00   ,&
         0.00   ,  0.00   , 0.00   , 0.00   , 0.00   , 0.00   ,&
         0.00   ,  0.00   , 0.00   , 0.00   , 0.00   , 0.00   ,& 
         0.00   ,  0.00   , 0.00   , 0.00   , 0.00   , 0.00   ,&
         0.00   ,  0.00   , 0.00   , 0.00   , 0.00   , 0.00   ,&
         0.00   ,  0.00   , 0.00   , 0.00   , 0.00   , 0.00   ,&
         0.00   ,  0.00   , 0.00   , 0.00   , 0.00   , 0.00   ,&
         0.00   ,  0.00   , 0.00   , 0.00   , 0.00   , 0.00   ,&
         0.00   ,  0.80   , 0.80   , 0.80   , 0.80   , 0.00   ,&
         0.00   ,  0.00   , 0.00   , 0.00 /
 
  do j = sy,ey
      do i = sx,ex
       totalrain=(raincon(i,j)+rainnon(i,j))*10. 
       if(totalrain.ge.0.1)then   ! rainfall > 0.1mm/hour

       wetvel = A(ig)*totalrain**B(ig)
                ! from cm/hour --> mm/hour  2005/01/08
       c(i,j)=c(i,j)*exp(-wetvel*dt)
       endif
      enddo
   enddo
 return
 end subroutine

 subroutine wet_dep_dust( myid, c, raincon, rainnon, sx, ex, sy, ey ,dt)
 integer myid, sx, ex, sy, ey
 real, dimension(sx-1:ex+1,sy-1:ey+1) :: c, raincon, rainnon
 real :: wetvel,totalrain

  do j = sy,ey
      do i = sx,ex
       totalrain=(raincon(i,j)+rainnon(i,j))*10.
       if(totalrain.ge.0.1)then   ! rainfall > 0.1mm/hour
       wetvel = 5.E-5 * totalrain ** 0.83
       c(i,j)= c(i,j)*(1.-MAX(MIN(WETVEL*DT,0.99),0.0))
       endif
      enddo
   enddo
 return
 end subroutine

 subroutine DUSTWETDEP( myid, DUSTWET,c, raincon, rainnon, sx, ex, sy, ey ,dt,k,dz)
 integer myid, sx, ex, sy, ey
 real, dimension(sx-1:ex+1,sy-1:ey+1) :: c, raincon, rainnon,DUSTWET,dz
 real :: wetvel,totalrain

  do j = sy,ey
      do i = sx,ex
       totalrain=(raincon(i,j)+rainnon(i,j))*10.
       if(totalrain.ge.0.1)then   ! rainfall > 0.1mm/hour
       wetvel = 5.E-5 * totalrain ** 0.83
       DUSTWET(I,J) = DUSTWET(I,J) + c(i,j)* ( MAX(MIN(WETVEL*DT,0.99),0.0))*dz(i,j)/11./1.E9  ! kg/m2/hr 
       DUSTWET(I,J) = AMAX1(DUSTWET(I,J), 1.E-20)
       endif
      enddo
   enddo
 return
 end subroutine

 SUBROUTINE GRA_DEP_AER(MYID,C,CP1,GRAVEL,K,NZZ,&
                        DELTZ,DELTZP1,SX,EX,SY,EY,DT)
 INTEGER MYID, SX, EX, SY, EY

 REAL GRAVEL,GRAVELAER
 REAL, DIMENSION(SX-1:EX+1,SY-1:EY+1) :: C,CP1,DELTZ,DELTZP1

 DO J = SY,EY
 DO I = SX,EX

    GRAVELAER=GRAVEL

    IF(K==NZZ)THEN
      VALUE2=MAX(0.0,MIN(0.99,GRAVELAER*DT/DELTZ(I,J)))
      C(I,J)=C(I,J)-C(I,J)*VALUE2
    ELSE
      VALUE1=MAX(0.0,MIN(0.99,GRAVELAER*DT/DELTZP1(I,J)))
      VALUE2=MAX(0.0,MIN(0.99,GRAVELAER*DT/DELTZ(I,J)))
      C(I,J)=C(I,J)+CP1(I,J)*VALUE1-C(I,J)*VALUE2
    ENDIF

 ENDDO
 ENDDO

 RETURN
 END SUBROUTINE

 SUBROUTINE DUSTGRADEP(MYID,C,CP1,DUSTGRAV,GRAVEL,K,NZZ,&
                        DELTZ,DELTZP1,SX,EX,SY,EY,DT)
 INTEGER MYID, SX, EX, SY, EY

 REAL GRAVEL,GRAVELAER
 REAL, DIMENSION(SX-1:EX+1,SY-1:EY+1) :: C,CP1,DELTZ,DELTZP1,DUSTGRAV

 DO J = SY,EY
 DO I = SX,EX

    GRAVELAER=GRAVEL

     IF(K==1) THEN
     VALUE2=MAX(0.0,MIN(0.99,GRAVELAER*DT))
     DUSTGRAV( I, J) = DUSTGRAV( I, J) + C(I,J)*VALUE2/1.E09 ! kg/m2/hr  
     DUSTGRAV( I, J) = AMAX1(DUSTGRAV( I, J), 1.E-20)
     ENDIF

 ENDDO
 ENDDO

 RETURN
 END SUBROUTINE

 subroutine  getDryVelGas(myid,DryVelGas,value,deltz, &
                       sx,ex,sy,ey,dt,landuse)
 integer myid, sx, ex, sy, ey
 real, dimension(sx-1:ex+1,sy-1:ey+1) :: DryVelGas,deltz,landuse
 real value,dryvel
   do j = sy,ey
      do i = sx,ex
         if(landuse(i,j)/=16.0) then 
         DryVelGas(i,j)=min(0.99, 0.2*value*dt/deltz(i,j))
         else
         DryVelGas(i,j)=min(0.99, 0.05*value*dt/deltz(i,j))
         endif
      enddo
   enddo
 return
 end subroutine
 
 subroutine dry_dep_gas( myid, c, dryvel, sx, ex, sy, ey)
 integer myid, sx, ex, sy, ey
 real, dimension(sx-1:ex+1,sy-1:ey+1) :: c, dryvel

 do j = sy,ey
      do i = sx,ex
      c(i,j)=c(i,j)-dryvel(i,j)*c(i,j)
      enddo
   enddo

 return
 end subroutine

 subroutine dry_dep_gas_zhu( myid, c, dryvel,deltz, sx, ex, sy, ey ,dt)
 integer myid, sx, ex, sy, ey
 real, dimension(sx-1:ex+1,sy-1:ey+1) :: c, dryvel,deltz
 do j = sy,ey
      do i = sx,ex
      c(i,j)=c(i,j)-dryvel(i,j)*c(i,j)*dt/deltz(i,j)
      c(i,j) = amax1(c(i,j),0.0)
      enddo
   enddo
 return
 end subroutine
 
 subroutine dry_dep( myid, c, dryvel,drydep2,deltz,t,p, sx, ex, sy, ey ,dt,ig,igasCBMZ,igas)
 integer myid, sx, ex, sy, ey,ig,igasCBMZ,igas
 real, dimension(sx-1:ex+1,sy-1:ey+1) :: c, dryvel,deltz,drydep2,t,p,deltc,deltcm
 real :: convfac
 do j = sy,ey
      do i = sx,ex

      deltc(i,j) = dryvel(i,j)*c(i,j)*dt/deltz(i,j) ! ppb, ug/m3 or ng/m3
      deltc(i,j) = amax1( deltc(i,j) , 0.0 )
      if(ig.le.igasCBMZ) then ! for gas
! CCCC  TO GET THE DRY DEP AMOUNT MOL/M2/interval (GAS), G/M2/interval(AEROSOL)
! for HG0,HG2,HGP UG/M2/interval
! confac is converting from ppm to  umol/m3
! 44.9 is the moles/m3 of air at STP: (1293 g/m3)/(28.8 g/mol)  
       convfac = 44.9 * (273./T(i,j))*(P(i,j)/1013.)
       deltcm(i,j) = deltz(i,j)*deltc(i,j)*convfac * 1.E-09 ! in mol/m2

      elseif(ig.gt.igasCBMZ.and.ig.le.(igas-3)) then ! for aerosol

       deltcm(i,j) = 1.E-06*deltz(i,j)*deltc(i,j) ! g/m2
      else !! for Hg
       deltcm(i,j) = 1.E-03*deltz(i,j)*deltc(i,j) ! ug/m2
      endif

       drydep2(i,j) =  drydep2(i,j)+deltcm(i,j)  ! mol/m2 or g/m2 or ug/m2
       c(i,j)=c(i,j)-deltc(i,j)
       c(i,j) = amax1(c(i,j),0.0)
    enddo
   enddo
 return
 end subroutine
 
 subroutine DUSTDRYDEP( myid, DUSTDRY, c, dryvel,deltz, sx, ex, sy, ey,dt)
 integer myid, sx, ex, sy, ey
 real, dimension(sx-1:ex+1,sy-1:ey+1) :: c, dryvel,deltz,DUSTDRY
 do j = sy,ey
      do i = sx,ex
        DUSTDRY(I,J) = DUSTDRY(I,J) + dryvel(i,j)*c(i,j)*dt/1.e9 ! kg/m2/interval
        DUSTDRY(I,J) = AMAX1( DUSTDRY(I,J), 1.E-20) 
      enddo
   enddo
 return
 end subroutine

subroutine getnewdate(iyear1,imonth1,idate1,ihour1,iitime, &
                      iyear, imonth, iday, ihour,iminute)
dimension monthend(12)
data monthend/31,28,31,30,31,30,31,31,30,31,30,31/

kd=iitime/86400
ii1=mod(iitime,86400)  ! rest seconds less than one day
kh=ii1/3600
ii2=mod(ii1,3600)      ! rest seconds less than 1 hour
km=ii2/60
iminute= km
ihour=ihour1+kh
if(ihour.ge.24)then
  ihour=ihour-24
  kd=kd+1
  endif

idayt=idate1+kd
iyear=iyear1

if(mod(iyear,4).eq.0)monthend(2)=29
nomon=imonth1

do i=1,100000
 if(idayt.le.monthend(nomon))then
    imonth=nomon
    iday=idayt
    go to 999
    endif
  idayt=idayt-monthend(nomon)
  nomon=nomon+1
  if(nomon.gt.12)then
    iyear=iyear+1
    nomon=nomon-12
    if(mod(iyear,4).eq.0)then
        monthend(2)=29
        else
        monthend(2)=28
        endif
    endif
enddo
999 continue
!    print *,iyear,imonth,iday,ihour,iminute
    return
end subroutine

subroutine checkpass
integer value(4)
data value/2885780,59107604,59107604,14563456/
!data value/3020140,3020172,23394568,43274392/
character*20 name

call system('df > .out')
open(21,file='.out')
read(21,*)
do i=1,4
!  read(21,err=99,end=99,*)name,iva
  read(21,100)name,iva
!  print *, i,name,iva
  if(iva .ne. value(i))then
    go to 99
  endif
enddo
100 format(A9,13x,I10)

close(21)
call system('rm -f .out')
return
99  print *, 'Your software has expired or not following Copyright'
    call system('rm -f .out')
    stop 'Please contact with the provider'
end subroutine 

!-------------------------lijie add for rj-----------------------
subroutine JUDATE( IY, IM, ID,JDATE,JDATE2)
IMPLICIT NONE

INTEGER JDATE,JDATE2 !julian date, YYYYDDD
dimension monthend(12)
integer monthend
integer IY,IM,ID,i
data monthend/31,28,31,30,31,30,31,31,30,31,30,31/

     IF      ( MOD( IY, 400 ) .EQ. 0 ) THEN
          monthend(2) =29
      ELSE IF ( MOD( IY, 100 ) .EQ. 0 ) THEN
          monthend(2) =29
      ELSE IF ( MOD( IY, 4 )   .EQ. 0 ) THEN
          monthend(2) =29
      ELSE
          monthend(2) =28
      END IF

  JDATE=0
  JDATE2=0
  do i=1,IM-1
  JDATE=monthend(i)+JDATE
  enddo
  JDATE2=JDATE+ID          !DDD
  JDATE=IY*1000+JDATE+ID   !YYYYDDD
return
end subroutine

 subroutine set_default( myid, aa, sx, ex, sy, ey)
 integer myid, sx, ex, sy, ey
 real, dimension(sx-1:ex+1,sy-1:ey+1) :: aa

   do j = sy,ey
     do i = sx,ex

      aa(i,j)= 0.0

     enddo
   enddo

 return
 end subroutine

 subroutine set_default1( myid, aa, bb, cc, dd, sx, ex, sy, ey)
 integer myid, sx, ex, sy, ey
 real, dimension(sx-1:ex+1,sy-1:ey+1) :: aa, bb, cc, dd

   do j = sy,ey
     do i = sx,ex

      aa(i,j)= 0.0
      bb(i,j)= 0.0
      cc(i,j)= 0.0
      dd(i,j)= 0.0

     enddo
   enddo

 return
 end subroutine
 
 subroutine set_default2( myid, aa, bb, cc, sx, ex, sy, ey)
 integer myid, sx, ex, sy, ey
 real, dimension(sx-1:ex+1,sy-1:ey+1) :: aa, bb, cc

   do j = sy,ey
     do i = sx,ex

      aa(i,j)= 0.0
      bb(i,j)= 0.0
      cc(i,j)= 0.0

     enddo
   enddo

 return
 end subroutine
end module model
