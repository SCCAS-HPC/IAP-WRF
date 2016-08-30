module adv_mark
use getsmark
contains
  subroutine adv_hori_mark( myid, c, u, v, deltx, delty, &
                       sx, ex, sy, ey ,ne,nx,ny,dt,k,ktop,&
                       ISMMAX,SM)
  integer myid, sx, ex, sy, ey, it,nx,ISMMAX
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c,c1, u, v, deltx, delty
  real, dimension(sx-1:ex+1) :: Q0,QN,u1D,DXX,DD0,DEN0,DEN1 
  real, dimension(sy-1:ey+1) :: Q0y,QNy,u1Dy,DXXy,DD0y,DEN0y,DEN1y 
  integer i, j, k
  real,dimension(sx-1:ex+1,sy-1:ey+1) :: ktop
  REAL,ALLOCATABLE,DIMENSION(:,:)     ::  TMPsmconv
  REAL,DIMENSION(ISMMAX,SX-1:EX+1,SY-1:EY+1) :: SM
  INTEGER                   :: ITYPE
  
  data d0/1./

  dt0 = dt+1  !! by chenhs
  do j = sy,ey
  do i = sx,ex
   dt0 = min(dt0, deltx(i,j)/abs(u(i,j)),delty(i,j)/abs(v(i,j)))
  enddo
  enddo

  istep=dt/dt0+1
  dtt=dt/istep

  do j=sy-1,ey+1
   do i=sx-1,ex+1
    c1(i,j)= c(i,j)
   enddo
  enddo   
  
  
do 10 itt=1,istep
 
  do j = sy,ey        ! do I advection
      do i = sx-1,ex+1
           Q0(i)=c1(i,j)
           u1D(i)=u(i,j) 
           DXX(i)=deltx(i,j)
           DEN0(i)=d0
           DEN1(i)=1.-dtt/DXX(i)*(d0*u(i,j)-d0*u(i-1,j))
           DD0(i)=DEN0(i)
      enddo
       KKKTOP   =  KTOP(I,J)

       ITYPE = 1

      ALLOCATE(TMPsmconv(ISMMAX,SX-1:EX+1)) 
      DO  is=1,ismMax
       DO I=SX-1,EX+1
        TMPsmconv(is,I) = SM(IS,I,J)
!         IF(I==84.AND.J==43.AND.K==1.AND.IS==3) PRINT*, SM(3,I,J) , 'inner1'
       ENDDO 
      ENDDO 
      
      
      call advec1d_mark2(sx,ex,Q0,QN,u1D,DEN0,DEN1,DXX,DD0,dtt,ISMMAX,&
                 TMPsmconv,KKKTOP,NE,NX,NY,ITYPE,K)

      DO  is=1,ismMax
       DO I=SX-1,EX+1
          SM(IS,I,J) = TMPsmconv(is,I)
!          IF(I==84.AND.J==43.AND.K==1.AND.IS==3) PRINT*, SM(3,I,J) , 'inner2'
       ENDDO   
      ENDDO 
                 
      do i = sx,ex
          c1(i,j)=amax1(QN(i),0.)
      enddo
      DEALLOCATE(TMPsmconv)
  enddo
  
   do i = sx,ex   ! do j direction
        do j = sy-1,ey+1
           Q0y(j)=c1(i,j)
           u1Dy(j)=v(i,j)
           DXXy(j)=delty(i,j)
           DEN0y(j)=1.-dtt/deltx(i,j)*(d0*u(i,j)-d0*u(i-1,j))
           DD0y(j)=d0
         enddo
      if(sy==1)then
           DEN1y(0)=DEN0y(0)
        do j = sy,ey+1
           DEN1y(j)=DEN0y(j)-dtt/delty(i,j)*(d0*v(i,j)-d0*v(i,j-1))
        enddo
       else
        do j = sy-1,ey+1
           DEN1y(j)=DEN0y(j)-dtt/delty(i,j)*(d0*v(i,j)-d0*v(i,j-1))
        enddo
       endif

       ITYPE = 2
       ALLOCATE(TMPsmconv(ISMMAX,SY-1:EY+1))
         
      DO  is=1,ismMax
        DO J=SY-1,EY+1
          TMPsmconv(is,J) = SM(IS,I,J)
        ENDDO  
      ENDDO  
      call advec1d_mark2(sy,ey,Q0y,QNy,u1Dy,DEN0y,DEN1y,DXXy,DD0y,dtt,ISMMAX,&
              TMPsmconv,KKKTOP,NE,NX,NY,ITYPE,K)
      
      DO  is=1,ismMax        
       DO J=SY-1,EY+1
        SM(IS,I,J) = TMPsmconv(is,J)
       ENDDO 
      ENDDO 

      DEALLOCATE(TMPsmconv)
      do j = sy,ey
          c1(i,j)=amax1(QNy(j),0.)
      enddo
  enddo

 10 continue
 return
  end subroutine

subroutine advec1d_mark2(sx,ex,Q0,QN,U,DEN0,DEN1,DXX,DD0,DT,ISMMAX,&
              TMPsmconv,KKKTOP,NE,NX,NY,ITYPE,K)
integer sx,ex
real, dimension(sx-1:ex+1)   :: Q0,QN,U,DEN0,DEN1,DXX,DD0
real, dimension(sx-1:ex+1)   :: FLUX,VCMAX,VCMIN
integer,dimension(sx-1:ex+1) :: IMXMN
INTEGER                      :: ITYPE,IMONTHERBOUNDATY
INTEGER                      :: ISMMAX,KKKTOP,NE,NX,NY
REAL,DIMENSION(ISMMAX,SX-1:EX+1) :: TMPsmconv
REAL,DIMENSION(ISMMAX):: smthis,smother
DATA ZR0,LSHP/0.,0./

! identify local max and min, specigy mxing ratio limit at new time

DO 5 I=sx,ex
   IMXMN(i)=0
   if(Q0(I).GE.AMAX1(Q0(I-1),Q0(I+1)).or.  &
      Q0(I).LE.AMIN1(Q0(I-1),Q0(I+1))    ) IMXMN(I)=1
   CK1=Q0(I)
   IF(U(I ).lt.ZR0)CK1=Q0(I+1)
   IF(U(I-1).ge.ZR0)CK1=Q0(I-1)
   DEN1(I) = amax1(DEN1(I) , 1.E-20)
   DXX(I) = amax1(DXX(I), 1.E-20)
   VCMAX(I)=aMAX1(Q0(I),CK1)
 5 VCMIN(I)=amin1(q0(i),ck1)


! update mixing ratios and limit Flux going up where u>0
 IF(U(sx-1).GE.ZR0) FLUX(sx-1) = Q0(sx-1)*U(sx-1)*DT*DD0(sx-1)
 DO 10  I=sx,ex
   IF(U(I).lt.ZR0) GOTO 10
   IF(U(I-1).lt.ZR0) THEN
      FLUX(I)=Q0(I)*U(I)*DT*DD0(I)  ! outflow only cell
      ELSE
        X1=DT*U(I)/DXX(I)
        CF1=CFCALC_mark2(X1,Q0(I),Q0(I-1),Q0(I+1),IMXMN(I-1),IMXMN(I+1),LSHP)
        FLUXI=DD0(I)*U(I)*DT*CF1/DXX(I)
        QN(I)=AMAX1( VCMIN(I) , AMIN1( VCMAX(I),  &
             (Q0(I)*DEN0(I) - FLUXI + FLUX(I-1)/DXX(I))/DEN1(I) ))
!!!!!!!!!!!!!!!!!!!!!!!!!SOURCE MARK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF(FLUX(I-1)>=0.0) THEN
          DO is=1,ismMax
            smthis(is) = TMPsmconv(is,I)
            smother(is)= TMPsmconv(is,I-1)
!            IF(I==84.AND.IJ==43.AND.K==1) PRINT*,IS,smthis(IS),'inner inner 01'
          ENDDO

                    
!         deltc1 = (Q0(I)*DEN0(I)/2.+FLUX(I-1)/DXX(I))/DEN1(I)-Q0(I)/2.
          deltc1 = FLUX(I-1)/DXX(I)/DEN1(I)
          deltc1 = MAX(deltc1,0.0)
          Q0(I) = MAX (Q0(I), 1.E-20)

!         IF(I==84.AND.IJ==43.AND.K==1) PRINT*,Q0(I),deltc1,smother(3)
      

         IMONTHERBOUNDARY = 0.0 
         IF(NE==1) THEN
          IF(ITYPE==1) THEN
            IF(I==1) IMONTHERBOUNDARY=1
          ELSE           
            IF(I==1) IMONTHERBOUNDARY=1
          ENDIF      
         ENDIF
          IF(IMONTHERBOUNDARY==1) THEN
             CALL GetBounChange(deltc1,Q0(I),smthis,1,ismMax)
          ELSE
             CALL SMmixing(Q0(I),smthis,deltc1,smother,ismMax)     
          ENDIF  ! IMONTHERBOUNDARY 

          DO is=1,ismMax
            TMPsmconv(is,I)=smthis(is)
          ENDDO 
!          IF(I==84.AND.IJ==43.AND.K==1) PRINT*,smthis(3),'inner inner 1'                  
        ENDIF !FLUX(I-1) 

        QQ0 = Q0(I) + FLUX(I-1)/DXX(I)/DEN1(I) 
        QQ0 = MAX(QQ0, 1.E-20)

        IF(FLUXI<=0.) THEN
          DO is=1,ismMax    
           smthis(is) = TMPsmconv(is,I)
           smother(is)= TMPsmconv(is,I+1) 
          ENDDO 

!          deltc1 = (Q0(I)*DEN0(I)/2. - FLUXI)/DEN1(I)-Q0(I)/2.
          deltc1 = - FLUXI/DEN1(I)
          
         IMONTHERBOUNDARY = 0.0
         IF(NE==1) THEN
           IF(ITYPE==1) THEN       
             IF(I==NX) IMONTHERBOUNDARY=1      
           ELSE
             IF(I==NY) IMONTHERBOUNDARY=1      
           ENDIF  
         ENDIF  

         IF(IMONTHERBOUNDARY==1) THEN
           CALL GetBounChange(deltc1,QQ0,smthis,1,ismMax)
         ELSE
           CALL SMmixing(QQ0,smthis,deltc1,smother,ismMax)
         ENDIF

         DO is=1,ismMax
           TMPsmconv(is,I)=smthis(is)
         ENDDO
  
       ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!         IF(I==84.AND.IJ==43.AND.K==1) PRINT*,smthis(3),'inner inner 2'
        FLUX(I)= DXX(I)*(Q0(I)*DEN0(I)-QN(I)*DEN1(I)) + FLUX(I-1)
   ENDIF
  
 10 CONTINUE 

! update mixing ratios and limit Flux going down where u<0

  IF(U(ex).lt.ZR0)FLUX(ex)=Q0(ex+1)*U(ex)*DT*DD0(ex)
  DO 20 I=ex,sx,-1
  IF(U(I-1).GE.ZR0)THEN   
         if(U(I).lt.ZR0) QN(I) =    &     ! inflow-only cell
         (Q0(I)*DEN0(I)-FLUX(I)/DXX(I)+FLUX(I-1)/DXX(I))/DEN1(I)  
!!!!!!!!!!!!!!!!!!!!!!!!!SOURCE MARK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(FLUX(I)<=0) THEN

        DO is=1,ismMax
         smthis(is) = TMPsmconv(is,I)
         smother(is)= TMPsmconv(is,I+1)
        ENDDO 

!        deltc1=(Q0(I)*DEN0(I)/2.-FLUX(I)/DXX(I))/DEN1(I)-Q0(I)/2.
         deltc1 = -FLUX(I)/DXX(I)/DEN1(I)
        deltc1 = MAX(deltc1,0.0)
        Q0(I) = MAX(Q0(I), 1.E-20)
              
        IMONTHERBOUNDARY = 0.0
        IF(NE==1) THEN
          IF(ITYPE==1) THEN       
           IF(I==NX) IMONTHERBOUNDARY=1
          ELSE 
           IF(I==NY) IMONTHERBOUNDARY=1   
          ENDIF
        ENDIF

        IF(IMONTHERBOUNDARY==1) THEN   
          CALL GetBounChange(deltc1,Q0(I),smthis,1,ismMax)       
        ELSE  
          CALL SMmixing(Q0(I),smthis,deltc1,smother,ismMax)
        ENDIF

        DO is=1,ismMax
         TMPsmconv(is,I)=smthis(is)
        ENDDO
         IF(I==84.AND.IJ==43.AND.K==1) PRINT*,smthis(3),'inner inner 3'

      ENDIF  ! FLUX(I)
        
        QQ0 = (Q0(I)*DEN0(I)- FLUX(I)/DXX(I))/DEN1(I) 
        QQ0 =MAX(1.E-20, QQ0)       
 
      IF(FLUX(I-1)>=0) THEN
        DO is=1,ismMax      
          smthis(is) = TMPsmconv(is,I)
          smother(is)= TMPsmconv(is,I-1)
        ENDDO  

!       deltc1=(Q0(I)*DEN0(I)/2.+FLUX(I-1)/DXX(I))/DEN1(I)-Q0(I)/2.
        deltc1 = FLUX(I-1)/DXX(I)/DEN1(I)
       deltc1 = MAX(deltc1,0.0)

       
       IMONTHERBOUNDARY=0.0
       IF(NE==1) THEN
        IF(ITYPE==1) THEN
         IF(I==1) IMONTHERBOUNDARY=1
        ELSE
         IF(I==1) IMONTHERBOUNDARY=1
        ENDIF 
       ENDIF 

       IF(IMONTHERBOUNDARY==1) THEN
         CALL GetBounChange(deltc1,QQ0,smthis,1,ismMax)
       ELSE
         CALL SMmixing(QQ0,smthis,deltc1,smother,ismMax)
       ENDIF

       DO is=1,ismMax
        TMPsmconv(is,I)=smthis(is)   
       ENDDO
!        IF(I==84.AND.IJ==43.AND.K==1) PRINT*,smthis(3),'inner inner 4'

     ENDIF ! FLUX(I-1)     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
   ELSE
     X1=DT*ABS(U(I-1))/DXX(I)
     CF1=CFCALC_mark2(X1,Q0(I),Q0(I+1),Q0(I-1),IMXMN(I+1),IMXMN(I-1),LSHP)
     IF(U(I).ge.ZR0)CF1=Q0(I)   ! outflow only cell
        FLUXex1=DD0(I-1)*U(I-1)*DT*CF1/DXX(I)
        QN(I)=AMAX1( VCMIN(I) , AMIN1( VCMAX(I),  &
             (Q0(I)*DEN0(I)-FLUX(I)/DXX(I)+FLUXex1)/DEN1(I)))
!!!!!!!!!!!!!!!!!!!!!!!!SOURCE MARK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF(FLUXex1>=0) THEN
      DO is=1,ismMax      
        smthis(is) = TMPsmconv(is,I)
        smother(is)= TMPsmconv(is,I-1)
      ENDDO

!      deltc1 = (Q0(I)*DEN0(I)/2.+FLUXex1)/DEN1(I) - Q0(I)/2.
       deltc1 = FLUXex1/DEN1(I)
       deltc1 = MAX(deltc1,0.0)
       Q0(I) = MAX( Q0(I), 1.E-20)

      
      IMONTHERBOUNDARY=0.0
      IF(NE==1) THEN
       IF(ITYPE==1) THEN       
        IF(I==1) IMONTHERBOUNDARY=1      
       ELSE 
        IF(I==1) IMONTHERBOUNDARY=1       
       ENDIF
      ENDIF

      IF(IMONTHERBOUNDARY==1) THEN
        CALL GetBounChange(deltc1,Q0(I),smthis,1,ismMax)
      ELSE
        CALL SMmixing(Q0(I),smthis,deltc1,smother,ismMax)
      ENDIF

      DO is=1,ismMax
        TMPsmconv(is,I)=smthis(is)
      ENDDO
       
     ENDIF  !FLUX(I-1)

     QQ0 = Q0(I)+FLUXex1/DEN1(I) 
     QQ0 = MAX( QQ0, 1.E-20)
     
     IF(FLUX(I)<=0) THEN
       DO is=1,ismMax       
        smthis(is) = TMPsmconv(is,I)
        smother(is)= TMPsmconv(is,I+1)
       ENDDO 
       
!       deltc1 = (Q0(I)*DEN0(I)/2.-FLUX(I)/DXX(I))/DEN1(I)- Q0(I)/2.
        deltc1 = -FLUX(I)/DXX(I)/DEN1(I) 
        deltc1 = MAX(deltc1,0.0)       
       IMONTHERBOUNDARY=0.0
       IF(NE==1) THEN
         IF(ITYPE==1) THEN
          IF(I==NX) IMONTHERBOUNDARY=1       
         ELSE 
          IF(I==NY) IMONTHERBOUNDARY=1  
         ENDIF
       ENDIF

       IF(IMONTHERBOUNDARY==1) THEN
         CALL GetBounChange(deltc1,QQ0,smthis,1,ismMax)
       ELSE
         CALL SMmixing(QQ0,smthis,deltc1,smother,ismMax)
       ENDIF
      
       DO is=1,ismMax
         TMPsmconv(is,I)=smthis(is)
       ENDDO
       
     ENDIF ! FLUX(I)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             
        FLUX(I-1)= DXX(I)*(QN(I)*DEN1(I)-Q0(I)*DEN0(I)) + FLUX(I)
   ENDIF 
 20 continue
   RETURN
   END subroutine
     
  FUNCTION CFCALC_mark2(X1,VCUP1,VCUP2,VCDW1,IMXUP2,IMXDN1,LSHP)
  ! function to calculate mixing ratio in the fuild moving into neighnor
  ! cell durint DT; X1=Courant No. set LSHP=1 for full sharpening
  ALFA=1.
  IF(IMXDN1.gt.0)ALFA=1.7763-0.5*X1
  IF(IMXUP2.gt.0)ALFA=1.2578+0.58502*X1
  IF(LSHP==1)ALFA=5.
  IF(X1.lt.0.5)THEN
     ALFA=MIN(ALFA,.98*6./(1.-2.*X1))
     CF=VCUP1*(1.+ALFA/6.*(2.*X1-1.)) + &
        VCDW1*ALFA/12.*(4-5.*X1) + VCUP2*ALFA/12.*(X1-2.)
    ELSE
     X1=MAX(X1,1.E-20)
     CF=VCUP1*(1.-ALFA*(1./X1 + 2.*X1 -3.)/6.) + &
        VCDW1*ALFA*(1./X1-X1)/12. + VCUP2*ALFA*(1./X1+5.*X1-6.)/12.
   ENDIF
! limit outflow mixing ratio to reasonable mixing ratio
    CFCALC_mark2=MIN(MAX(CF,MIN(VCUP1,VCDW1)),MAX(VCUP1,VCDW1))
    RETURN
    END FUNCTION
  
          SUBROUTINE ADV_VERT_MARK(MYID,C,W,U2D,V2D,DELTZ,&
                DELTX,SX,EX,SY,EY,NZZ,DT,IG,ISMMAX,SMCONV,KKTOP)

          INTEGER :: MYID,SX,EX,SY,EY,ISMMAX,NZZ
          REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1,NZZ) :: C,W,U2D,V2D
          REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1,NZZ) :: DELTX,DELTZ
          REAL    :: DT,D0 
          REAL,DIMENSION(NZZ)                     :: Q0,QN,DXX0,DXX, &
                                                     DEN0,DEN1,U,DD0
          REAL,DIMENSION(ISMMAX,SX-1:EX+1,SY-1:EY+1,NZZ) ::  SMCONV
          REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1)     ::  KKTOP
          REAL,DIMENSION(ISMMAX,NZZ)    ::  TMPsmconv
          DATA D0 /1.0/                                           
          
          DO I=SX,EX
           DO J=SY,EY
            DO  K=2,NZZ-1
             DXX0(K)  =  DELTX(I,J,K)
             DXX (K)  =  DELTZ(I,J,K)
             DEN0(K)  =  1.- DT/DXX0(K)*(D0*U2D(I,J,K)-D0*U2D(I-1,J,K))&
                           - DT/DXX0(K)*(D0*V2D(I,J,K)-D0*V2D(I,J-1,K))
             DEN1(K)  =  DEN0(K)-DT/DXX(K)*(D0*W(I,J,K)-D0*W(I,J,K-1)) 
             
             DD0 (K)  =  DEN0(K)
             Q0  (K)  =  C(I,J,K)
             U   (K)  =  W(I,J,K)
            ENDDO  !K

             DD0 (1)  =  D0
             U   (1)  =  W(I,J,1)
             Q0  (1)  =  C(I,J,1)
             Q0  (NZZ)=  C(I,J,NZZ)
             KKKTOP   =  KKTOP(I,J)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!source appointent
              DO  is=1,ismMax
              DO K=1,NZZ  
               TMPsmconv(is,k)=smconv(is,i,j,k)
              ENDDO !K
              ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             
            CALL ADVEC1D_MARK(Q0,QN,U,DEN0,DEN1,DT,DXX,DD0,NZZ,ISMMAX,&
                        TMPsmconv,KKKTOP,I,J,IG)

            DO K=2,NZZ-1

             C(I,J,K) = QN(K)
             DO  is=1,ismMax
              smconv(is,i,j,k)=TMPsmconv(is,k)
             ENDDO !IS 
             
            ENDDO  !K

            
           ENDDO   !J
          ENDDO !I

          RETURN
          END subroutine

         SUBROUTINE ADVEC1D_MARK(Q0,QN,U,DEN0,DEN1,DT,DXX,DD0,NZZ,&
                          ISMMAX,TMPsmconv,KKKTOP,II,JJ,IG)
         INTEGER              :: NZZ
         REAL,DIMENSION(NZZ)  :: Q0,QN,U,DEN0,DEN1,DXX,DD0
         REAL                 :: DT,CK1,X1,FLUXI,FLUXIM1,CF1,QQ0
         REAL,DIMENSION(NZZ)  :: FLUX,VCMAX,VCMIN,IMXMN
         INTEGER              :: ISMMAX,KKKTOP 
         REAL,DIMENSION(ISMMAX,NZZ) :: TMPsmconv
         REAL,DIMENSION(ISMMAX):: smthis,smother
         DATA ZR0,LSHP / 0.0,0.0/

         
         IMXMN = 0.0 
         
          
!!!!  INDENTIFY LOCAL MAX AND MIN,SPECIFY MIXING RATIO LIMITS AT NEW TIME
 
        DO 5 I=2,NZZ-1

         IMXMN(I) = 0
         IF(Q0(I).GE.MAX(Q0(I-1),Q0(I+1)).OR. &
            Q0(I).LE.MIN(Q0(I-1),Q0(I+1)) ) IMXMN(I)=1
         CK1 = Q0(I)
         IF(U(I).LT.ZR0)   CK1 = Q0(I+1)
         IF(U(I-1).GE.ZR0) CK1 = Q0(I-1)
         VCMAX(I) = MAX(Q0(I),CK1)
      5  VCMIN(I) = MIN(Q0(I),CK1)

!!!!    UPDATE MIXING RATIOS AND LIMITS FLUXES GOING UP WHERE U2D >0
             
        IF(U(1).GE.ZR0) FLUX(1) = Q0(1)*U(1)*DT*DD0(1)
        DO 10 I=2,NZZ-1
         IF(U(I).LT.ZR0) GOTO 10
         IF(U(I-1).LT.ZR0) THEN
            FLUX(I)=Q0(I)*U(I)*DT*DD0(I)              ! OUTFLOW-ONLY CELL    
         ELSE
            X1   = DT*U(I)/DXX(I)
            CF1  = CFCALC_MARK(X1,Q0(I),Q0(I-1),Q0(I+1),int(IMXMN(I-1)),int(IMXMN(I+1)),LSHP) ! juanxiong he
            FLUXI= DD0(I)*U(I)*DT*CF1/DXX(I)
            QN(I)= MAX(VCMIN(I) ,MIN(VCMAX(I), &
                  (Q0(I)*DEN0(I)-FLUXI+FLUX(I-1)/DXX(I))/DEN1(I) ) )
!!!!!!!!!!!!!!!!!!!!!!!!!SOURCE MARK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            IF(II==53.AND.JJ==30.AND.IG==11) &
!             PRINT*,I,Q0(I)*DEN0(I)/2./DEN1(I)+Q0(I)/2.-Q0(I),FLUX(I-1)/DXX(I)/DEN1(I),(-FLUXI)/DEN1(I)
                  
            IF(FLUX(I-1)>=0.0) THEN
              DO is=1,ismMax
               smthis(is) = TMPsmconv(is,I)
               smother(is)= TMPsmconv(is,I-1) 
              ENDDO      

!              deltc1 = (Q0(I)*DEN0(I)/2.+FLUX(I-1)/DXX(I))/DEN1(I)-Q0(I)/2.
              deltc1 =  FLUX(I-1)/DXX(I)/DEN1(I)
              deltc1 = MAX(deltc1,1.E-20)
              Q0(I)  = MAX(Q0(I),1.E-20)

              IF(I>KKKTOP) THEN
               do is =1 , ismMax
               if(is==2)then
                 smthis(is)=1.0
               else
                 smthis(is) =0.0
               endif
               enddo 
              ELSE
               call SMmixing(Q0(I),smthis,deltc1,smother,ismMax)        
              ENDIF !I

              DO is=1,ismMax
               TMPsmconv(is,I)=smthis(is)
              ENDDO
       
            ENDIF     ! FLUX 
            
             QQ0 = (Q0(I)*DEN0(I)+FLUX(I-1)/DXX(I))/DEN1(I)

            IF(FLUXI<= 0.) THEN
             DO is=1,ismMax       
              smthis(is) = TMPsmconv(is,I)
              smother(is)= TMPsmconv(is,I+1)
             ENDDO 
            
            deltc1 = (-FLUXI)/DEN1(I)
            deltc1 = max(deltc1,1.E-20)
            QQ0    = max(QQ0, 1.E-20)            

            IF(I==KKKTOP) THEN
              call GetBounChange(deltc1,QQ0,smthis,2,ismMax)       
            ELSE IF(I>KKKTOP)THEN
             do is = 1 ,ismMax
               if(is==2) then
                 smthis(is)=1.0
               else
                 smthis(is) = 0.0
               endif
             enddo
            ELSE
              call SMmixing(QQ0,smthis,deltc1,smother,ismMax)
            ENDIF !I             
            
            DO is=1,ismMax
              TMPsmconv(is,I)=smthis(is)
            ENDDO
 
           ENDIF !FLUX(I)            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                  
            FLUX(I) = DXX(I)*(Q0(I)*DEN0(I)-QN(I)*DEN1(I))+FLUX(I-1)   
         ENDIF        
      10 CONTINUE

!!!!    UPDATE MIXING RATIOS AND LIMIT FLUXES GOING DOWN WHERE U2D < 0
       
       IF(U(NZZ-1).LT.ZR0) FLUX(NZZ-1)= Q0(NZZ)*U(NZZ-1)*DT*DD0(NZZ-1)
       DO 20 I=NZZ-1,2,-1
        IF(U(I-1).GE.ZR0) THEN
          IF(U(I).LT.ZR0) QN(I) =  &                   !INFLOW-ONLY CELL
           (Q0(I)*DEN0(I)-FLUX(I)/DXX(I)+FLUX(I-1)/DXX(I))/DEN1(I)

!!!!!!!!!!!!!!!!!!!!!!!!!!  SOURCE MARK !!!!!!!!!!!!!!!!!!!!!!!!!!!
!        IF(II==53.AND.JJ==30.AND.IG==11) &
!          PRINT*,I,Q0(I)*DEN0(I)/2./DEN1(I)+Q0(I)/2.-Q0(I),-FLUX(I)/DXX(I)/DEN1(I),FLUX(I-1)/DXX(I)/DEN1(I)

           IF(FLUX(I)<=0) THEN
              DO is=1,ismMax
               smthis(is) = TMPsmconv(is,I)
               smother(is)= TMPsmconv(is,I+1)
              ENDDO
              
!              deltc1=(Q0(I)*DEN0(I)/2.-FLUX(I)/DXX(I))/DEN1(I)-Q0(I)/2.
              deltc1 = -FLUX(I)/DXX(I)/DEN1(I)
              deltc1 = MAX(deltc1,1.E-20)
              Q0(I)  = max(Q0(I),1.E-20)              

              IF(I==KKKTOP) THEN
                call GetBounChange(deltc1,Q0(I),smthis,2,ismMax)    
              ELSE IF(I>KKKTOP)THEN
               do is = 1 , ismMax
                if(is==2) then
                  smthis(is)=1.0
                else
                  smthis(is) = 0.0
                endif
               enddo
              ELSE
                call SMmixing(Q0(I),smthis,deltc1,smother,ismMax)      
              ENDIF !I

              DO is=1,ismMax
                TMPsmconv(is,I)=smthis(is)
              ENDDO   
             
            ENDIF  !FLUX(I)
             
             QQ0 = Q0(I) - FLUX(I)/DXX(I)/DEN1(I)
             
            IF(FLUX(I-1)>=0) THEN
              DO is=1,ismMax            
               smthis(is) = TMPsmconv(is,I)
               smother(is)= TMPsmconv(is,I-1)
              ENDDO 

!             deltc1=(Q0(I)*DEN0(I)/2.+FLUX(I-1)/DXX(I))/DEN1(I)-Q0(I)/2.
              deltc1 = FLUX(I-1)/DXX(I)/DEN1(I)
              deltc1 = MAX(deltc1,1.E-20)
              QQ0    = MAX(QQ0,1.E-20)
             
             IF(I>KKKTOP) THEN
               do is = 1 ,ismMax
                 if(is==2) then
                   smthis(is)=1.0      
                 else 
                   smthis(is) = 0.0
                 endif
               enddo
             ELSE
               call SMmixing(QQ0,smthis,deltc1,smother,ismMax)
             ENDIF  !I

             DO is=1,ismMax
               TMPsmconv(is,I)=smthis(is)
             ENDDO  
            
            ENDIF  ! FLUX(I-1) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
        ELSE
          X1  = DT*ABS(U(I-1))/DXX(I)
          CF1 = CFCALC_MARK(X1,Q0(I),Q0(I+1),Q0(I-1),int(IMXMN(I+1)),int(IMXMN(I-1)),LSHP)  ! junxiong he
          IF(U(I).GE.ZR0) CF1 = Q0(I)                  !OUTFLOW-ONLY CELL
          FLUXIM1 = DD0(I-1)*U(I-1)*DT*CF1/DXX(I)
          QN(I)= MAX( VCMIN(I), MIN( VCMAX(I), &
                 (Q0(I)*DEN0(I)-FLUX(I)/DXX(I)+FLUXIM1)/DEN1(I) ))
!!!!!!!!!!!!!!!!!!!!!!!!SOURCE MARK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         IF(II==53.AND.JJ==30.AND.IG==11) &
!           PRINT*,I,Q0(I)*DEN0(I)/2./DEN1(I)+Q0(I)/2.-Q0(I),FLUXIM1/DEN1(I),-FLUX(I)/DXX(I)/DEN1(I)
                 
          IF(FLUXIM1>=0) THEN               
            DO is=1,ismMax      
             smthis(is) = TMPsmconv(is,I)
             smother(is)= TMPsmconv(is,I-1)
            ENDDO     

!           deltc1 = (Q0(I)*DEN0(I)/2.+FLUXIM1)/DEN1(I) - Q0(I)/2.
            deltc1 = FLUXIM1/DEN1(I)
            deltc1 = MAX(deltc1,1.E-20)
            Q0(I) = MAX(Q0(I),1.E-20)              

            IF(I>KKKTOP) THEN
              do is = 1,ismMax
               if(is==2) then
                 smthis(is)=1.0
               else
                 smthis(is) = 0.0
               endif
              enddo
            ELSE
              call SMmixing(Q0(I),smthis,deltc1,smother,ismMax)      
            ENDIF  !I  

            DO is=1,ismMax
             TMPsmconv(is,I)=smthis(is)
            ENDDO 
            
          ENDIF  ! FLUX(I-1)

           QQ0 = (Q0(I)*DEN0(I)+FLUXIM1)/DEN1(I) 

          IF(FLUX(I)<=0) THEN
            DO is=1,ismMax
             smthis(is) = TMPsmconv(is,I)
             smother(is)= TMPsmconv(is,I+1)
            ENDDO
           
           
!           deltc1 = (Q0(I)*DEN0(I)/2.-FLUX(I)/DXX(I))/DEN1(I)-Q0(I)/2.  
           deltc1 = -FLUX(I)/DXX(I)/DEN1(I)
           deltc1 = MAX(deltc1,1.E-20)
           QQ0 = MAX(QQ0,1.E-20)
             
            IF(I==KKKTOP) THEN             
             call GetBounChange(deltc1,QQ0,smthis,2,ismMax)
            ELSE IF(I>KKKTOP) THEN
             do is = 1, ismMax
              if(is==2) then
               smthis(is)=1.0
              else
               smthis(is) = 0.0
              endif
             enddo 
            ELSE 
             call SMmixing(QQ0,smthis,deltc1,smother,ismMax)
            ENDIF      !I
            
           DO is=1,ismMax          
             TMPsmconv(is,I)=smthis(is)
           ENDDO
                              
          ENDIF ! FLUX(I)        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                 
          FLUX(I-1) = DXX(I)*(QN(I)*DEN1(I)-Q0(I)*DEN0(I)) + FLUX(I)
        ENDIF
      20 CONTINUE
      
        RETURN
       END   subroutine 

      FUNCTION CFCALC_MARK(X1,VCUP1,VCUP2,VCDW1,IMXUP2,IMXDN1,LSHP)
! function to calculate mixing ratio in the fuild moving into
! neighnor
! cell durint DT; X1=Courant No. set LSHP=1 for full
! sharpening
        REAL :: CFCALC_MARK   
        ALFA=1.
        IF(IMXDN1.gt.0)ALFA=1.7763-0.5*X1
        IF(IMXUP2.gt.0)ALFA=1.2578+0.58502*X1
        IF(LSHP==1)ALFA=5.
        IF(X1.lt.0.5)THEN
          ALFA=MIN(ALFA,.98*6./(1.-2.*X1))
          CF=VCUP1*(1.+ALFA/6.*(2.*X1-1.)) + &
                   VCDW1*ALFA/12.*(4-5.*X1) + VCUP2*ALFA/12.*(X1-2.)
        ELSE
          X1=MAX(X1,1.E-20)
          CF=VCUP1*(1.-ALFA*(1./X1 + 2.*X1 -3.)/6.) + &
             VCDW1*ALFA*(1./X1-X1)/12. + VCUP2*ALFA*(1./X1+5.*X1-6.)/12.
        ENDIF
! limit outflow mixing ratio to reasonable mixing ratio
        CFCALC_MARK=MIN(MAX(CF,MIN(VCUP1,VCDW1)),MAX(VCUP1,VCDW1))
        RETURN
        END  FUNCTION      

end module adv_mark
