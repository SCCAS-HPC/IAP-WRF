module diffusion_mark
use getsmark
contains

! Zifa 2006/10/05 revised
! For Source Mark
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine dif_hori_mark( myid, c0, kh, deltx, delty, &
                       sx, ex, sy, ey ,nx,ny,dt,   &
                      ismMax,sm, ne,k,ktop)
  integer,parameter :: ismMaxSet=500
  integer myid, sx, ex, sy, ey, it
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c0, kh, deltx, delty,ktop
  real, dimension(ismMax,sx-1:ex+1,sy-1:ey+1) :: sm

  integer i, j,k
  integer sx1,ex1,sy1,ey1
  real, allocatable, dimension(:,:) :: c
  real smthis(ismMaxSet)
  real smother(ismMaxSet)

  allocate(c(sx-1:ex+1,sy-1:ey+1))

  c=c0

  sx1=sx
  ex1=ex
  sy1=sy
  ey1=ey

    do j = sy1,ey1
      do i = sx1,ex1
         IF(k>ktop(i,j)) then
          do is =1 , ismMax
          if(is.eq.2) then
               sm(is,i,j)=1.0
          else
               sm(is,i,j)=0.0
          endif
          enddo
         ENDIF
      enddo
    enddo


  
  do j = sy1,ey1
    do i = sx1,ex1

      deltc1=dt*kh(i,j)*(c(i+1,j)-c(i,j))/(deltx(i,j)*deltx(i,j))
       
  
      if(deltc1 > 0 ) then !!!!!!!!!!!!!!!!!!
      c00=max(c(i,j),1.E-20)
         do is=1,ismMax
           smthis(is)=sm(is,i,j)
           smother(is)=sm(is,i+1,j)
         enddo
          imotherboundary =0  
          ! to judge if it is mother domain
            if(ne==1)then
              if(i==1)imotherboundary =1
              if(i==nx)imotherboundary =1
              if(j==1)imotherboundary =1
              if(j==ny)imotherboundary =1
            endif

          if(imotherboundary==1)then
!--------------------------lijie modify for sm----------------------
            IF(i==ex1) then 
            call GetBounChange(deltc1,c00,smthis,1,ismMax)
            ENDIF
!-------------------------finish------------------------------------
           else
            call SMmixing(c00,smthis,deltc1,smother,ismMax)
          endif
          do is=1,ismMax
           sm(is,i,j)=smthis(is)
          enddo
      endif  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       

      deltc2=dt*kh(i,j)*(c(i-1,j)-c(i,j))/(deltx(i,j)*deltx(i,j))
     if(deltc2 > 0 ) then !!!!!!!!!!!!!!!!!!
      c00=max(c(i,j),1.E-20)
         do is=1,ismMax    
           smthis(is)=sm(is,i,j)
           smother(is)=sm(is,i-1,j)
         enddo 
          imotherboundary =0
          ! to judge if it is mother domain
            if(ne==1)then
              if(i==1)imotherboundary =1 
              if(i==nx)imotherboundary =1
              if(j==1)imotherboundary =1
              if(j==ny)imotherboundary =1
            endif

          if(imotherboundary==1)then
!--------------------------lijie modify for sm----------------------
            IF(i==sx1) then
            call GetBounChange(deltc2,c00,smthis,1,ismMax)
            ENDIF
!-------------------------finish------------------------------------
           else
            call SMmixing(c00,smthis,deltc2,smother,ismMax)
          endif
          do is=1,ismMax
           sm(is,i,j)=smthis(is)
          enddo
      endif  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      c(i,j)=c(i,j)+deltc1+deltc2
    enddo
  enddo

  do j = sy1,ey1
    do i = sx1,ex1

     deltc1=dt*kh(i,j)*(c(i,j+1)-c(i,j))/(delty(i,j)*delty(i,j))

     if(deltc1 > 0 ) then !!!!!!!!!!!!!!!!!!
      c00=max(c(i,j),1.E-20)
         do is=1,ismMax
           smthis(is)=sm(is,i,j)
           smother(is)=sm(is,i,j+1)
         enddo
          imotherboundary =0
          ! to judge if it is mother domain
            if(ne==1)then
              if(i==1)imotherboundary =1
              if(i==nx)imotherboundary =1
              if(j==1)imotherboundary =1
              if(j==ny)imotherboundary =1
            endif

          if(imotherboundary==1)then
!--------------------------lijie modify for sm----------------------
            IF(j==ey1) then
            call GetBounChange(deltc1,c00,smthis,1,ismMax)
            ENDIF
!-------------------------finish------------------------------------
           else
            call SMmixing(c00,smthis,deltc1,smother,ismMax)
          endif
          do is=1,ismMax
           sm(is,i,j)=smthis(is)
          enddo
      endif  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     deltc2=dt*kh(i,j)*(c(i,j-1)-c(i,j))/(delty(i,j)*delty(i,j))

     if(deltc2 > 0 ) then !!!!!!!!!!!!!!!!!!
      c00=max(c(i,j),1.E-20)
         do is=1,ismMax
           smthis(is)=sm(is,i,j)
           smother(is)=sm(is,i,j-1)
         enddo
          imotherboundary =0
          ! to judge if it is mother domain
            if(ne==1)then
              if(i==1)imotherboundary =1
              if(i==nx)imotherboundary =1
              if(j==1)imotherboundary =1
              if(j==ny)imotherboundary =1
            endif

          if(imotherboundary==1)then
!--------------------------lijie modify for sm----------------------
            IF(j==sy1) then
            call GetBounChange(deltc2,c00,smthis,1,ismMax)
            ENDIF
!-------------------------finish------------------------------------
           else
            call SMmixing(c00,smthis,deltc2,smother,ismMax)
          endif
          do is=1,ismMax
           sm(is,i,j)=smthis(is)
          enddo
      endif  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      c(i,j)=c(i,j)+deltc1+deltc2
    enddo
  enddo

 deallocate(c)

  return
  end subroutine
!
  subroutine dif_vert_mark( myid, c10,c0,cp10,kz,deltz,  &
                       sx, ex, sy, ey ,k, nzz,dt,ismMax, &
             sm,smb,smp)
  integer,parameter :: ismMaxSet=100
  integer myid, sx, ex, sy, ey, it
  integer i, j, k
  real smthis(ismMaxSet)
  real smother(ismMaxSet)

  real, dimension(sx-1:ex+1,sy-1:ey+1) :: kz
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c10,c0,cp10
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: deltz
  real, dimension(ismMax,sx-1:ex+1,sy-1:ey+1):: sm,smb,smp

  real, allocatable, dimension(:,:) ::c1,c,cp1

  allocate(c1(sx-1:ex+1,sy-1:ey+1))
  allocate(c(sx-1:ex+1,sy-1:ey+1))
  allocate(cp1(sx-1:ex+1,sy-1:ey+1))


  c1=c10
  c=c0
  cp1=cp10

 ! goto 333

  do j = sy,ey
    do i = sx,ex
      deltc1=dt*kz(i,j)*(cp1(i,j)-c(i,j))/(deltz(i,j)*deltz(i,j))  !Zifa wang original
      c00=c(i,j)

      if(deltc1>0)then
        do is=1,ismMax
           smthis(is)=sm(is,i,j)
           smother(is)=smp(is,i,j)
          enddo
        if(k==nzz-1)then
          call GetBounChange(deltc1,c00,smthis,2,ismMax)
        else if(k>nzz-1) then
          smthis(2)=1.0
        else 
          call SMmixing(c00,smthis,deltc1,smother,ismMax)
        endif

      do is=1,ismMax
           sm(is,i,j)=smthis(is)
       enddo 

      endif

      deltc2=dt*kz(i,j)*(-c(i,j)+c1(i,j))/(deltz(i,j)*deltz(i,j))   !Zifa wang original
    if(deltc2>0)then
        do is=1,ismMax
           smthis(is)=sm(is,i,j)
           smother(is)=smb(is,i,j)
          enddo
        
        if(k>nzz-1) then
          smthis(2)=1.0
        else
          call SMmixing(c00,smthis,deltc2,smother,ismMax)
        endif   
       do is=1,ismMax
           sm(is,i,j)=smthis(is)
       enddo 
     endif

       c(i,j)=c(i,j)+deltc1+deltc2
    enddo
  enddo

 333 continue
  deallocate(c1,c,cp1)
  return
  end subroutine

       subroutine diffu_vert_mark(myid,&
              c1tmp,c_1tmp,ctmp,cp1tmp,&              
              kv1,kv_1,kv,kvp1,h1,h_1,h,hp1,terrain,&
              dep1,dep_1,dep,depp1,&
              fcon,pbl,kpbl,Gc_MOLWT,ig,igas,&
              ismMax,sm1, sm,smb,smp,ktop,&
              sx,ex,sy,ey,dt,k)
              
            integer, parameter :: ismMaxSet=500
            real smthis(ismMaxSet)
            real smother(ismMaxSet)
            integer  :: myid, sx, ex, sy, ey,it
            real, dimension(ismMax,sx-1:ex+1,sy-1:ey+1) :: sm1,sm,smb,smp
            
            
            real, dimension(sx-1:ex+1,sy-1:ey+1) :: c1tmp,c_1tmp,ctmp,cp1tmp
            real, dimension(sx-1:ex+1,sy-1:ey+1) :: kv1,kv_1,kv,kvp1
            real, dimension(sx-1:ex+1,sy-1:ey+1) :: h1,h_1,h,hp1
            real, dimension(sx-1:ex+1,sy-1:ey+1) :: dep1,dep_1,dep,depp1
            real, dimension(sx-1:ex+1,sy-1:ey+1) :: fcon,pbl,terrain
            integer, dimension(sx-1:ex+1,sy-1:ey+1) :: kpbl,ktop
            real :: M2u,M2di,M2dii,cm1tmp,cm_1tmp,cmtmp,cmp1tmp
            real :: dt,dtt,dt0,delta1,delta2,delta3,delta4,delta5
            real :: GC_MOLWT(igas)
           do j = sy,ey
           do i = sx,ex
           
            M2u= fcon(i,j)*kv1(i,j)/dep1(i,j)/(pbl(i,j)-h1(i,j)+terrain(i,j)-dep1(i,j)/2.)
            M2di=M2u*(pbl(i,j)-h(i,j)+terrain(i,j)+dep(i,j)/2.)/dep(i,j)
            M2dii=M2u*(pbl(i,j)-hp1(i,j)+terrain(i,j)+depp1(i,j)/2.)/depp1(i,j)
           

            dtt=dt
!------convert ppbv --> mass mixing ratio------

             cm1tmp=c1tmp(i,j)*GC_MOLWT(ig)/29.
             cm_1tmp=c_1tmp(i,j)*GC_MOLWT(ig)/29.
             cmtmp=ctmp(i,j)*GC_MOLWT(ig)/29.
             cmp1tmp=cp1tmp(i,j)*GC_MOLWT(ig)/29.
!-----------------------convert  over---------------------                  
!              if(i==46.and.j==32.and.k==2.and.ig==17.and.ne==1)print*,cm_1,cm,cmp1,'222' 
           delta1=0.0
           delta2=0.0
           delta3=0.0
           delta4=0.0
           delta5=0.0
           IF(kpbl(i,j)==1) goto 999
           if(k==1) then
            delta1=-cm1tmp*M2u*(pbl(i,j)-h1(i,j)+terrain(i,j)-dep1(i,j)/2.)/dep1(i,j)
            delta2=0.0
            delta3= M2dii*cmp1tmp*depp1(i,j)/dep(i,j)
           else if(k==kpbl(i,j)) then
            delta1= M2u*cm1tmp
            delta2=-M2di*cmtmp
            delta3= 0.0      
           else if(k<kpbl(i,j)) then 
            delta1= M2u*cm1tmp
            delta2=-M2di*cmtmp
            delta3= M2dii*cmp1tmp*depp1(i,j)/dep(i,j)
           endif 
999         continue
           if(k<=kpbl(i,j)) then
           else
            delta1=0.0;delta2=0.0;delta3=0.0 
           endif
             delta4=0.0; delta5=0.0

             c00=cmtmp
           if(delta3>0.0) then
             do is=1,ismMax
               smthis(is)=sm(is,i,j)
               smother(is)=smp(is,i,j)
               IF(i==46.and.j==55.and.k==1.and.is==3) PRINT*, smthis(is),smother(is),'delta3'
             enddo
                 
             if(k==ktop(i,j))then                   
               call GetBounChange(delta3*dtt,c00,smthis,2,ismMax)
             else if(k>ktop(i,j)) then
              do is = 1 ,ismMax
               if(is==2) then
                 smthis(is)=1.0
               else
                 smthis(is) =0.0 
               endif
              enddo
             else
               call SMmixing(c00,smthis,delta3*dtt,smother,ismMax)
             endif

             do is=1,ismMax
               sm(is,i,j)=smthis(is)
               IF(i==46.and.j==55.and.k==1.and.is==3) PRINT*, smthis(is),'delta3'    
             enddo 
           endif !delta3

            c00=c00+delta3*dtt 

            if(delta2>0.0) then
            do is=1,ismMax
              smthis(is)=sm(is,i,j)
              smother(is)=smb(is,i,j)
               IF(i==46.and.j==55.and.k==1.and.is==3) PRINT*, smthis(is),smother(is),'delta2'
            enddo       
  
            if(k>ktop(i,j)) then              
               do is = 1 ,ismMax
               if(is==2) then
                 smthis(is)=1.0
               else
                 smthis(is) =0.0
               endif
               enddo
            else
              call SMmixing(c00,smthis,delta2*dtt,smother,ismMax)      
            endif    

            do is=1,ismMax                               
              sm(is,i,j)=smthis(is)                           
               IF(i==46.and.j==55.and.k==1.and.is==3) PRINT*, smthis(is),'delta2'
            enddo  
           endif !delta2 
             c00=c00+delta2*dtt
  
           
           if(delta1>0.0) then     
            do is=1,ismMax                                              
             smthis(is)=sm(is,i,j)                              
             smother(is)=sm1(is,i,j)                                  
               IF(i==46.and.j==55.and.k==1.and.is==3) PRINT*, smthis(is),smother(is),'delta1'
            enddo      

            if(k>ktop(i,j)) then                                 
               do is = 1 ,ismMax
               if(is==2) then
                 smthis(is)=1.0
               else
                 smthis(is) =0.0
               endif
               enddo
            else  
             call SMmixing(c00,smthis,delta1*dtt,smother,ismMax)    
            endif          
                    
            do is=1,ismMax
              sm(is,i,j)=smthis(is)
               IF(i==46.and.j==55.and.k==1.and.is==3) PRINT*, smthis(is),smother(is),'delta1'
            enddo          

           endif !delta1
             c00=c00+delta1*dtt

                
           enddo !i
           enddo !j
            return 
            end subroutine

                   subroutine dif_vert_conv_mark( myid,c1m,cm,cp1m,c0m,kz,kz1,kzb,deltz,deltz1,&
                                deltzt,deltzb,pbl,&
                                xmol,heiz,heiz1,heizb,heizt,TERRAIN,&
                              sx, ex, sy, ey,dt,ig,k,nzz,ismMax,sm,smb,smp,sm0)
      integer,parameter :: ismMaxSet=100
      integer myid, sx, ex, sy, ey, it,ig
      real, dimension(sx-1:ex+1,sy-1:ey+1) :: kz,pbl,xmol,heiz,kz1,heiz1,heizb,heizt,kzb,TERRAIN
      !kz1 and heiz1 is the kz and height of l layers,kzb and heizb is the i-1/2 layers,heizt i+1/2 layer
      real, dimension(sx-1:ex+1,sy-1:ey+1) :: c1m,cm,cp1m,c0m !c0m is the mixing ratios in the 1st layer
      real, dimension(sx-1:ex+1,sy-1:ey+1) :: deltz,deltz1,deltzt,deltzb !delta1 is depth of 1st layer      
      real :: M2u,M2di,M2dit,fconv,kzz,kzzb 
      real :: const1,const2,karman,alpha1
      integer i, j, k     
      real :: smthis(ismMaxSet) 
      real :: smother(ismMaxSet) 
      real, dimension(ismMax,sx-1:ex+1,sy-1:ey+1):: sm,smb,smp,sm0

      real, allocatable, dimension(:,:) ::c1,c,cp1,c0
      allocate(c1(sx-1:ex+1,sy-1:ey+1))
      allocate(c(sx-1:ex+1,sy-1:ey+1))
      allocate(cp1(sx-1:ex+1,sy-1:ey+1))
      allocate(c0(sx-1:ex+1,sy-1:ey+1))

      c1=c1m
      c=cm
      cp1=cpm
      c0=c0m

! this scheme is using ACM2 scheme(non-local scheme)
      do j = sy,ey
          do i = sx,ex
             karman=0.4
             alpha1=7.2
             const1=(karman**(-2./3.))/(0.1*alpha1)
              if(xmol(i,j)<0) then
             const2=(-pbl(i,j)/(xmol(i,j)+1.e-6))**(-1./3.)
             fconv=1./(1.+const1*const2)
              else 
             fconv=0.
              endif 
              if((heiz(i,j)-TERRAIN(i,j))>pbl(i,j)) fconv=0.0 !confirm the highest layer in pbl
              fconv=min(fconv,0.3)
               fconv =0.      
             M2u=fconv*kz1(i,j)/(deltz1(i,j)*(pbl(i,j)-heiz1(i,j)+TERRAIN(i,j)))
             M2di=M2u*(pbl(i,j)-heizb(i,j)+TERRAIN(i,j))/deltz(i,j)
             M2dit=M2u*(pbl(i,j)-heiz(i,j)+TERRAIN(i,j))/deltzt(i,j)
             kzz=kz(i,j)*(1.-fconv)
             kzzb=kzb(i,j)*(1.-fconv)
             dt0=deltz(i,j)*deltz(i,j)/(10.*kz(i,j))
             iii= dt/dt0+1
             dtt=dt/iii
             do it=1,iii
                 c00=c(i,j)
                 deltam01=M2u*c0(i,j)
                 deltam02=-M2di*c(i,j)
                 deltam03= M2dit*cp1(i,j)*deltzt(i,j)/deltz(i,j)
                 deltam04a= (kzz*(cp1(i,j)-c(i,j))/deltz(i,j))/deltz(i,j)
                 deltam05a= -(kzzb*(c(i,j)-c1(i,j))/deltzb(i,j))/deltz(i,j)
                 deltam04b = (kzz*(cp1(i,j))/deltz(i,j))/deltz(i,j)
                 deltam05b = -(kzzb*(-c1(i,j))/deltz(i,j))/deltz(i,j)
               if(k==1) then                 
                  c(i,j)=c(i,j)+dtt*(deltam04+deltam05+M2dit*(cp1(i,j)-c(i,j)))
               else if((heizt(i,j)-TERRAIN(i,j))>pbl(i,j).and.(heiz(i,j)-TERRAIN(i,j))<pbl(i,j)) then
                  c(i,j)=c(i,j)+dtt*(deltam04+deltam05+deltam01+deltam02)
               else      
                  c(i,j)=c(i,j)+dtt*(deltam01+deltam02+deltam03+deltam04+deltam05)
               endif
                 
!----------------to estimate the contribution from 1st layers--------------
              if(k==1) then                 
                deltc1=0.0
              else if((heizt(i,j)-TERRAIN(i,j))>pbl(i,j).and.(heiz(i,j)-TERRAIN(i,j))<pbl(i,j)) then
                deltc1=dtt*(deltam01+deltam02)
!                 deltc1=0.0
              else      
                deltc1=dtt*(deltam01+deltam02)
              endif
               IF (deltc1>0.0) then
                do is=1,ismMax
                 smthis(is)=sm(is,i,j)
                 smother(is)=sm0(is,i,j)
                enddo
                if(k>nzz-1) then
                 smthis(2)=1.0
                else
                  call SMmixing(c00,smthis,deltc1,smother,ismMax)
                endif
                                                      
                 do  is=1,ismMax
                   sm(is,i,j)=smthis(is)
                   sm(is,i,j)=max(smthis(is),0.)
                   sm(is,i,j)=min(smthis(is),1.)
                 enddo
                 
                ENDIF
             c00=c00+deltc1
!!----------------to estimate the diffusion contribution from next lower layers--------------
             if(k==1) then
                 deltc2a=0.0
                 deltc2b=0.0
             else if((heizt(i,j)-TERRAIN(i,j))>pbl(i,j).and.(heiz(i,j)-TERRAIN(i,j))<pbl(i,j))  then
                 deltc2a=dtt*deltam05a
                 deltc2b=dtt*deltam05b
             else
                 deltc2a=dtt*deltam05a
                 deltc2b=dtt*deltam05b
             endif    

               IF(deltc2a>0.0) THEN
                 do is=1,ismMax
                   smthis(is)=sm(is,i,j)
                   smother(is)=smb(is,i,j)
                 enddo
                 if(k>nzz-1) then
                   smthis(2)=1.0
                 else
                   call SMmixing(c00,smthis,deltc2a,smother,ismMax)
                 endif

                 do  is=1,ismMax
                    sm(is,i,j)=smthis(is)
                    sm(is,i,j)=max(smthis(is),0.)
                    sm(is,i,j)=min(smthis(is),1.)
                 enddo
               ENDIF        
              c00=c00+deltc2a
!!!----------------to estimate the convective contribution from next higher layers--------------

             if(k==1) then
                deltc3=dtt*deltam03
             else if((heizt(i,j)-TERRAIN(i,j))>pbl(i,j).and.(heiz(i,j)-TERRAIN(i,j))<pbl(i,j))   then
                deltc3=0.0
             else
                deltc3=dtt*deltam03
             endif

             IF(deltc3>0.0) THEN
              do is=1,ismMax
                 smthis(is)=sm(is,i,j)
                 smother(is)=smp(is,i,j)
              enddo
              if(k==nzz-1)then
                 call GetBounChange(deltc3,c00,smthis,2,ismMax)
              else if(k>nzz-1) then
                 smthis(2)=1.0
              else
                 call SMmixing(c00,smthis,deltc3,smother,ismMax)
              endif 

              do  is=1,ismMax
                sm(is,i,j)=smthis(is)
                sm(is,i,j)=max(smthis(is),0.)
                sm(is,i,j)=min(smthis(is),1.)
              enddo
             ENDIF               
              c00=c00+deltc3
!!!----------------to estimate the diffusion  contribution from next higher layers--------------             

             if(k==1) then
                deltc4a=dtt*deltam04a
                deltc4b=dtt*deltam04b
             else if((heizt(i,j)-TERRAIN(i,j))>pbl(i,j).and.(heiz(i,j)-TERRAIN(i,j))<pbl(i,j))    then
                deltc4a=dtt*deltam04a
                deltc4b=dtt*deltam04b                
             else
                deltc4a=dtt*deltam04a
                deltc4b=dtt*deltam04b                
             endif
                                                       
             IF(deltc4a>0.0) THEN
              do is=1,ismMax
                  smthis(is)=sm(is,i,j)
                  smother(is)=smp(is,i,j)
              enddo
              if(k==nzz-1)then
                call GetBounChange(deltc4a,c00,smthis,2,ismMax)      
              else if(k>nzz-1) then
                smthis(2)=1.0      
              else
                call SMmixing(c00,smthis,deltc4a,smother,ismMax)
              endif
               
              do  is=1,ismMax
                sm(is,i,j)=smthis(is)
                sm(is,i,j)=max(smthis(is),0.)
                sm(is,i,j)=min(smthis(is),1.)
              enddo
             ENDIF
             c00=c00+deltc4a          

             c(i,j)=max(0.,c(i,j))
 
          enddo !it 
         enddo !i
        enddo  !j
          deallocate(c1,c,cp1,c0)
        return
       end    subroutine                   
 
      SUBROUTINE diffus_mark (myid, ig, sx, ex, sy, ey, nlay, deltat,   &
      rkv, depth, tempk, press, concmark, atm, ismMax, smconv, kktop)   
!                                                                       
!----CAMx v4.40 061025                                                  
!                                                                       
!     DIFFUS drives 3-D diffusion of concentrations.  This version also 
!     performs diffusion of sensitivities if DDM is enabled.            
!                                                                       
!     Copyright 1996-2006                                               
!     ENVIRON International Corporation                                 
!                                                                       
!     Modifications:                                                    
!        4/17/00   Revised diffusion equations to weight fluxes by densi
!       12/07/01   added instructions for OMP                           
!        1/13/03   added deposited mass array                           
!       10/12/04   Multiple substeps as f(Kv) applied for vertical diffu
!                                                                       
!     Input arguments:                                                  
!        igrd              grid index                                   
!        ncol              number of columns                            
!        nrow              number of rows                               
!        nlay              number of layers                             
!        nspc              number of species                            
!        nsen              number of species times number of DDM paramet
!                          or 1, whichever is larger                    
!        deltat            time step (s)                                
!        dx                cell size in x-direction (m)                 
!        dy                cell size in y-direction (m)                 
!        idfin             map of nested grids in this grid             
!        vdep              deposition velocity (m/s)                    
!        rkx/y             horizontal diffusion coefficient (m2/s)      
!        rkv               vertical diffusion coefficient (m2/s)        
!        depth             layer depth (m)                              
!        tempk             temperature (K)                              
!        press             pressure (mb)                                
!        mapscl            map-scale factor at cell centroid            
!        conc              species concentrations (umol/m3)             
!        sens              DDM sensitivities (umol/m3/parameter unit)   
!        tarray2           CPU timing arguments (s)                     
!        strz              string for labeling the z diffusion process  
!        strxy             string for labeling the x/y diffusion process
!        ipa_xy            2-D gridded array to identify if cell is     
!                          in a IPRM sub-domain                         
!        ipa_lay           3-D gridded array to identify which IPRM sub-
!                          each layer is in                             
!c        nparddm           number of parameters for which sensitivities
!                          are calculated.  Should be zero if DDM is not
!                          enabled.                                     
!       smconv              the fraction of every tracer                
!                                                                       
!     Output arguments:                                                 
!        concmark              species concentrations (umol/m3)         
!        sens              DDM sensitivities (umol/m3/parameter unit)   
!        fluxes            fluxes across the boundaries (umol)          
!        depfld            2-D array of dry deposited mass (mol/ha, g/ha
!                                                                       
!     Routines Called:                                                  
!        VDIFFIMP                                                       
!                                                                       
!     Called by:                                                        
!        EMISTRNS                                                       
!                                                                       
!                                                                       
!                                                                       
                                                                        
!======================== Process Analysis Begin =======================
!                                                                       
!                                                                       
      REAL fcup (nlay), fcdn (nlay) 
      LOGICAL ldoipts 
      INTEGER ::nddmsp 
!                                                                       
!========================= Process Analysis End ========================
!                                                                       
      INTEGER ::sx, ex, sy, ey 
      INTEGER ::MXTRSP 
                          ! the sensitity here no use                   
      PARAMETER (MXTRSP = 0) 
      DIMENSION concmark (sx - 1:ex + 1, sy - 1:ey + 1, nlay), vdep (sx &
      - 1:ex + 1, sy - 1:ey + 1), kktop (sx - 1:ex + 1, sy - 1:ey + 1)  
      REAL rkx (sx - 1:ex + 1, sy - 1:ey + 1, nlay), rky (sx - 1:ex + 1,&
      sy - 1:ey + 1, nlay), rkv (sx - 1:ex + 1, sy - 1:ey + 1, nlay),   &
      depth (sx - 1:ex + 1, sy - 1:ey + 1, nlay), tempk (sx - 1:ex + 1, &
      sy - 1:ey + 1, nlay), press (sx - 1:ex + 1, sy - 1:ey + 1, nlay), &
      mapscl (sx - 1:ex + 1, sy - 1:ey + 1), dx (sx - 1:ex + 1, sy - 1: &
      ey + 1, nlay), atm (sx - 1:ex + 1, sy - 1:ey + 1, nlay), dair (sx &
      - 1:ex + 1, sy - 1:ey + 1, nlay)                                  
      DIMENSION c1d (nlay + nlay * MXTRSP), d1d (nlay), rk1d (nlay),    &
      ro1d (nlay)                                                       
      INTEGER nstepv (sx:ex, sy:ey) 
      REAL dtv, dtmp 
      REAL convfac 
!=======================source appointment                              
      INTEGER ::ismMax, myid 
      REAL ::smconv (ismMax, sx - 1:ex + 1, sy - 1:ey + 1, nlay) 
      REAL ::TMPsmconv (ismMax, nlay) 
!============================                                           
                !lijie                                                  
      mapscl = 1. 
                ! here no dry depositions                               
      vdep = 0.0 
      nddmsp = 0 
!-----------convert ppbv to umol/m3---                                  
      DO i = sx, ex 
      DO j = sy, ey 
      DO k = 1, nlay 
!          if(ig==11.and.i==57.and.j==37.and.k==8)                      
!     &   print*,concmark(57,37,8),smconv(2,57,37,8),'mark begin'       
                                                                        
      convfac = 44.9 * (273. / tempk (i, j, k) ) * (press (i, j, k)     &
      / 1013.)                                                          
      concmark (i, j, k) = concmark (i, j, k) * convfac / 1e+03 
      enddo 
      enddo 
      enddo 
!-----Entry point                                                       
!                                                                       
!-----Vertical diffusion                                                
!                                                                       
!                                                                       
!-----Determine max vertical diffusion time step and                    
!     apply within a sub-step loop                                      
!                                                                       
      DO 601 j = sy, ey 
         i1 = sx 
         i2 = ex 
         DO 501 i = i1, i2 
            dtv = deltat 
            IF (deltat.lt.600.) goto 401 
            DO k = 1, nlay - 1 
            dtmp = 0.75 * depth (i, j, k) * depth (i, j, k + 1) / rkv ( &
            i, j, k)                                                    
            dtv = amin1 (dtv, dtmp) 
            enddo 
            dtv = amax1 (dtv, 600.) 
  401       nstepv (i, j) = INT (0.999 * deltat / dtv) + 1 
  501    END DO 
  601 END DO 
!                                                                       
!$omp parallel default(shared)                                          
!$omp&  private(ispc,i1,i2,i,j,k,ro1d,d1d,c1d,rk1d,lk,                  
!$omp&          ioff,fluxbot,isen,rokp1,ldoipts,fcup,                   
!$omp&          fcdn,ipa_idx,n,dtv)                                     
!                                                                       
!$omp do schedule(dynamic)                                              
!                                                                       
      DO 60 j = sy, ey 
         i1 = sx 
         i2 = ex 
         DO 50 i = i1, i2 
!                                                                       
!-----Skip cells occupied by child grids; load 1-D arrays               
!                                                                       
            dtv = deltat / FLOAT (nstepv (i, j) ) 
            DO n = 1, nstepv (i, j) 
                                                                        
!=======================source appointent                               
            DO is = 1, ismMax 
            DO k = 1, nlay 
            TMPsmconv (is, k) = smconv (is, i, j, k) 
            ENDDO 
            ENDDO 
!==========================================                             
            DO k = 1, nlay 
            ro1d (k) = press (i, j, k) / tempk (i, j, k) 
            d1d (k) = depth (i, j, k) 
            c1d (k) = concmark (i, j, k) 
            enddo 
            DO k = 1, nlay - 1 
            rokp1 = press (i, j, k + 1) / tempk (i, j, k + 1) 
            rk1d (k) = rkv (i, j, k) * (ro1d (k) + rokp1) / 2. 
            enddo 
            rk1d (nlay) = 0. 
!                                                                       
!                                                                       
!======================== Process Analysis Begin =======================
!                                                                       
            ldoipts = .TRUE. 
!                                                                       
!========================= Process Analysis End ========================
!                                                                       
            CALL vdiffimp_mark (nlay, dtv, vdep (i, j), d1d, ro1d, rk1d,&
            c1d, nddmsp, fcup, fcdn, ldoipts, ismMax, TMPsmconv, kktop (&
            i, j), i, j)                                                
!                                                                       
!                                                                       
!                                                                       
            DO k = 1, nlay 
            concmark (i, j, k) = c1d (k) 
            enddo 
                                                                        
            DO is = 1, ismMax 
            DO k = 1, nlay 
            smconv (is, i, j, k) = TMPsmconv (is, k) 
            ENDDO 
            ENDDO 
                  ! nstepv                                              
            enddo 
   50    END DO 
   60 END DO 
!                                                                       
!-----------convert umol/m3 to ppb---                                   
      DO i = sx, ex 
      DO j = sy, ey 
      DO k = 1, nlay 
      convfac = 44.9 * (273. / tempk (i, j, k) ) * (press (i, j, k)     &
      / 1013.)                                                          
      concmark (i, j, k) = 1e03 * concmark (i, j, k) / convfac 
!        if(ig==11.and.i==57.and.j==37.and.k==8)                        
!     & print*,concmark(57,37,8),smconv(2,57,37,8),'mark after'         
                                                                        
      enddo 
      enddo 
      enddo 
      RETURN 
      END SUBROUTINE diffus_mark                    
      
      SUBROUTINE vdiffimp_mark (nn, dt, vdep, depth, rho, rkv, rr,      &
      nparddm, fcup, fcdn, ldoipts, ismMax, TMPsmconv, kktop, i, j)     
!                                                                       
!----CAMx v4.40 061025                                                  
!                                                                       
!     VDIFFIMP performs vertical diffusion of concentrations using      
!     an implicit method, where a tri-diagonal matrix is solved.        
!     This version also performs vertical diffusion of sensitivities    
!     if DDM is enabled.                                                
!                                                                       
!     Copyright 1996-2006                                               
!     ENVIRON International Corporation                                 
!                                                                       
!     Modifications:                                                    
!        4/17/00   Revised diffusion equations to weight fluxes by densi
!                                                                       
!     Input arguments:                                                  
!        nn                number of layers                             
!        dt                time step (s)                                
!        vdep              deposition velocity (m/s)                    
!        depth             layer depth (m)                              
!        rho               atmospheric density (mb/K)                   
!        rkv               vertical diffusion coefficient (mb m2/s/K)   
!        rr                species concentrations (umol/m3) followed by 
!                          sensitivities (umol/m3/parameter unit).      
!        nparddm           number of parameters for which sensitivities 
!                          are calculated.  Should be zero if DDM is not
!                          enabled.                                     
!        ldoipts           flag to calculate data for Process Analysis  
!                                                                       
!     Output arguments:                                                 
!        rr                species concentrations (umol/m3) followed by 
!                          sensitivities (umol/m3/parameter unit)       
!        fcup              change in layer concentration due to flux acr
!                          upper interface (umol/m3) -- FOR Process Anal
!        fcdn              change in layer concentration due to flux acr
!                          lower interface (umol/m3) -- FOR Process Anal
!                                                                       
!     Routines Called:                                                  
!        TRDIAG                                                         
!                                                                       
!     Called by:                                                        
!        DIFFUS                                                         
!                                                                       
!                                                                       
!======================== Process Analysis Begin =======================
!                                                                       
!                                                                       
      REAL aa_old (nn), cc_old (nn) 
      REAL fcup (nn), fcdn (nn) 
      LOGICAL ldoipts 
!                                                                       
!========================= Process Analysis End ========================
!                                                                       
      INTEGER ::nparddm 
      DIMENSION rr (nn + nn * nparddm), depth (nn), rkv (nn), rho (nn) 
      DIMENSION aa (nn), bb (nn), cc (nn), rrtmp (nn + nn * nparddm) 
!                                                                       
!==============source appointment                                       
      integer, parameter::ismMaxSet = 100 
      REAL smthis (ismMaxSet) 
      REAL smother (ismMaxSet) 
      INTEGER ::ismMax, kktop 
      REAL ::TMPsmconv (ismMax, nn) 
!===============================                                        
!-----Entry point                                                       
!                                                                       
!-----Lower boundary condition                                          
!                                                                       
      aa (1) = 0. 
      bb (1) = 1. + dt / depth (1) * (vdep + 2. * rkv (1) / (depth (2)  &
      + depth (1) ) / rho (1) )                                         
!                                                                       
!-----Upper boundary condition                                          
!                                                                       
      cc (nn) = 0. 
      bb (nn) = 1. + dt / depth (nn) * 2. * rkv (nn - 1) / (depth (nn - &
      1) + depth (nn) ) / rho (nn)                                      
!                                                                       
      DO k = 2, nn 
      aa (k) = - dt / depth (k) * 2. * rkv (k - 1) / (depth (k - 1)     &
      + depth (k) ) / rho (k - 1)                                       
      enddo 
      DO k = 1, nn - 1 
      cc (k) = - dt / depth (k) * 2. * rkv (k) / (depth (k + 1) + depth &
      (k) ) / rho (k + 1)                                               
      enddo 
      DO k = 2, nn - 1 
      bb (k) = 1. + dt / depth (k) * 2. * rkv (k - 1) / (depth (k - 1)  &
      + depth (k) ) / rho (k) + dt / depth (k) * 2. * rkv (k) / (depth (&
      k + 1) + depth (k) ) / rho (k)                                    
      enddo 
!                                                                       
!======================== Process Analysis Begin =======================
!                                                                       
      IF (ldoipts) then 
         DO k = 1, nn 
         aa_old (k) = 0. 
         cc_old (k) = 0. 
         IF (k.gt.1) aa_old (k) = aa (k) * rho (k - 1) 
         IF (k.lt.nn) cc_old (k) = cc (k) * rho (k + 1) 
         enddo 
      ENDIF 
!                                                                       
!========================= Process Analysis End ========================
!=======================source appoint=======                           
      DO k = 1, nn 
      rrtmp (k) = rr (k) 
      enddo 
!=========================                                              
!                                                                       
!                                                                       
!-----Solve the equations                                               
!                                                                       
      CALL trdiag_mark (aa, bb, cc, rr, nn, 1 + nparddm) 
!                                                                       
!                                                                       
!======================== Process Analysis Begin =======================
!                                                                       
      IF (ldoipts) then 
         DO k = 2, nn - 1 
         fcup (k) = (rr (k) / rho (k) - rr (k + 1) / rho (k + 1) )      &
         * cc_old (k)                                                   
         fcdn (k) = (rr (k) / rho (k) - rr (k - 1) / rho (k - 1) )      &
         * aa_old (k)                                                   
         enddo 
!                                                                       
!-----Lower boundary                                                    
!                                                                       
         fcup (1) = (rr (1) / rho (1) - rr (2) / rho (2) ) * cc_old (1) 
         fcdn (1) = - rr (1) * dt / depth (1) * vdep 
!                                                                       
!-----Upper boundary                                                    
!                                                                       
         fcup (nn) = 0.0 
         fcdn (nn) = (rr (nn) / rho (nn) - rr (nn - 1) / rho (nn - 1) ) &
         * aa_old (nn)                                                  
      ENDIF 
!                                                                       
!========================= Process Analysis End ========================
!                                                                       
!========================source appointment========                     
      DO k = 1, nn 
      IF (fcup (k).ge. 0) then 
         DO is = 1, ismMax 
         smthis (is) = TMPsmconv (is, k) 
         IF (k .eq. nn) then 
            smother (is) = TMPsmconv (is, k) 
         ELSE 
            smother (is) = TMPsmconv (is, k + 1) 
         ENDIF 
         enddo 
         deltc1 = fcup (k) 
         deltc1 = max (fcup (k), 1.e-20) 
         rrtmp (k) = max (rrtmp (k), 1.e-20) 
                                                                        
                                                                        
         IF (k .eq. kktop) then 
            CALL GetBounChange (deltc1, rrtmp (k), smthis, 2, ismMax) 
         ELSEIF (k>kktop) then 
            DO is = 1, ismMax 
            IF (is .eq. 2) then 
               smthis (is) = 1.0 
            ELSE 
               smthis (is) = 0.0 
            ENDIF 
            enddo 
         ELSE 
            CALL SMmixing (rrtmp (k), smthis, deltc1, smother, ismMax) 
         ENDIF 
                                                                        
         DO is = 1, ismMax 
         TMPsmconv (is, k) = smthis (is) 
         enddo 
                                                                        
             !fcup                                                      
      ENDIF 
                                                                        
      IF (fcdn (k) .ge. 0) then 
         DO is = 1, ismMax 
         smthis (is) = TMPsmconv (is, k) 
         IF (k .eq. 1) then 
            smother (is) = TMPsmconv (is, k) 
         ELSE 
            smother (is) = TMPsmconv (is, k - 1) 
         ENDIF 
         enddo 
         deltc2 = fcdn (k) 
         deltc2 = max (fcdn (k), 1.e-20) 
         rrtmp (k) = max (rrtmp (k), 1.e-20) 
                                                                        
                                                                        
         IF (k>kktop) then 
            DO is = 1, ismMax 
            IF (is .eq. 2) then 
               smthis (2) = 1.0 
            ELSE 
               smthis (2) = 0.0 
            ENDIF 
            enddo 
         ELSE 
            CALL SMmixing (rrtmp (k) + fcup (k), smthis, deltc2,        &
            smother, ismMax)                                            
         ENDIF 
                                                                        
         DO is = 1, ismMax 
         TMPsmconv (is, k) = smthis (is) 
         enddo 
            !fcdn                                                       
      ENDIF 
            !k                                                          
      enddo 
      RETURN 
      END SUBROUTINE vdiffimp_mark                  

      SUBROUTINE trdiag_mark (a, b, c, d, n, m) 
!                                                                       
!----CAMx v4.40 061025                                                  
!                                                                       
!     TRDIAG solves a tridiagnal matrix equation with multiple          
!     right-hand side arrays.                                           
!                                                                       
!     Copyright 1996-2006                                               
!     ENVIRON International Corporation                                 
!                                                                       
!     Modifications:                                                    
!        none                                                           
!                                                                       
!     Input arguments:                                                  
!        n                  size of the coefficient matrix              
!        m                  number of right-hand side arrays            
!        a                  first coefficient array                     
!        b                  second coefficient array                    
!        c                  third coefficient array                     
!        d                  right-hand side arrays                      
!                                                                       
!     Output arguments:                                                 
!        d                  solution arrays                             
!                                                                       
!     Routines called:                                                  
!        none                                                           
!                                                                       
!     Called by:                                                        
!        VRTSLV                                                         
!        VDIFFIMP                                                       
!                                                                       
      DIMENSION a (n), b (n), c (n), d (n, m) 
!                                                                       
!-----Entry point                                                       
!                                                                       
      nm1 = n - 1 
      DO 10 i = 1, nm1 
         c (i) = - c (i) / b (i) 
         i1 = i + 1 
         b (i1) = c (i) * a (i1) + b (i1) 
   10 END DO 
!                                                                       
      DO 40 k = 1, m 
         d (1, k) = d (1, k) / b (1) 
         DO 20 i = 1, nm1 
            i1 = i + 1 
            d (i1, k) = (d (i1, k) - d (i, k) * a (i1) ) / b (i1) 
   20    END DO 
!                                                                       
         DO 30 j = 1, nm1 
            i = n - j 
            i1 = i + 1 
            d (i, k) = d (i1, k) * c (i) + d (i, k) 
   30    END DO 
   40 END DO 
!                                                                       
      RETURN 
      END SUBROUTINE trdiag_mark                    
      
end module diffusion_mark      
