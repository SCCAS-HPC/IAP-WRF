module diffusion
contains
  subroutine dif_hori( myid, c, kh, deltx, delty, &
                       sx, ex, sy, ey ,nx,ny,dt)
  integer myid, sx, ex, sy, ey, it
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c, kh, deltx, delty
  integer i, j
  integer sx1,ex1,sy1,ey1

  sx1=sx
  ex1=ex
  sy1=sy
  ey1=ey
!  if(sx1.eq.1)sx1=2
!  if(ex1.eq.nx)ex1=nx-1
!  if(sy1.eq.1)sy1=2
!  if(ey1.eq.ny)ey1=ny-1

  do j = sy1,ey1
    do i = sx1,ex1
!      dt0=deltx(i,j)*deltx(i,j)/(4.*kh(i,j))
!      iii= dt/dt0+1
!      dtt=dt/iii
!      do it=1,iii
!      c(i,j)=c(i,j)+dtt*kh(i,j)*(c(i+1,j)-2*c(i,j)+c(i-1,j))/ &
      c(i,j)=c(i,j)+dt*kh(i,j)*(c(i+1,j)-2*c(i,j)+c(i-1,j))/ &
                    (deltx(i,j)*deltx(i,j))
!      enddo
    enddo
  enddo

  do j = sy1,ey1
    do i = sx1,ex1
!      dt0=delty(i,j)*delty(i,j)/(4.*kh(i,j))
!      iii= dt/dt0+1
!      dtt=dt/iii
!      do it=1,iii
!      c(i,j)=c(i,j)+dtt*kh(i,j)*(c(i,j+1)-2*c(i,j)+c(i,j-1))/ &
      c(i,j)=c(i,j)+dt*kh(i,j)*(c(i,j+1)-2*c(i,j)+c(i,j-1))/ &
                    (delty(i,j)*delty(i,j))
!      enddo
    enddo
  enddo

  return
  end subroutine
  
  subroutine dif_vert( myid, c1,c,cp1,kz,deltz,  &
                       sx, ex, sy, ey , dt)
  integer myid, sx, ex, sy, ey, it,ig
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: kz
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c1,c,cp1
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: deltz

  integer i, j, k
!print*,'55555555555',ig,k
  do j = sy,ey
    do i = sx,ex
      dt0=deltz(i,j)*deltz(i,j)/(4.*kz(i,j))
      iii= dt/dt0+1
      dtt=dt/iii
      do it=1,iii
      c(i,j)=c(i,j)+dtt*kz(i,j)*(cp1(i,j)-2*c(i,j)+c1(i,j))/ &
                (deltz(i,j)*deltz(i,j))
      enddo

    enddo
  enddo

  return
  end  subroutine

     subroutine dif_vert_conv( myid,c1,c,cp1,c0,kz,kz1,kzb,deltz,deltz1,&
                                deltzt,deltzb,pbl,&
                                xmol,heiz,heiz1,heizb,heizt,TERRAIN,&
                              sx, ex, sy, ey ,dt,ig,k)
      integer myid, sx, ex, sy, ey, it,ig
      real, dimension(sx-1:ex+1,sy-1:ey+1) :: kz,pbl,xmol,heiz,kz1,heiz1,heizb,heizt,kzb,TERRAIN
      !kz1 and heiz1 is the kz and height of l layers,kzb and heizb is the i-1/2 layers,heizt i+1/2 layer
      real, dimension(sx-1:ex+1,sy-1:ey+1) :: c1,c,cp1,c0 !c0 is the mixing ratios in the 1st layer
      real, dimension(sx-1:ex+1,sy-1:ey+1) :: deltz,deltz1,deltzt,deltzb !delta1 is depth of 1st layer      
      real :: M2u,M2di,M2dit,fconv,kzz,kzzb 
      real :: const1,const2,karman,alpha1
      integer i, j, k     

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
              fconv=0.0
              
             M2u=fconv*kz1(i,j)/(deltz1(i,j)*(pbl(i,j)-heiz1(i,j)+TERRAIN(i,j)))
             M2di=M2u*(pbl(i,j)-heizb(i,j)+TERRAIN(i,j))/deltz(i,j)
             M2dit=M2u*(pbl(i,j)-heiz(i,j)+TERRAIN(i,j))/deltzt(i,j)
             kzz=kz(i,j)*(1.-fconv)
             kzzb=kzb(i,j)*(1.-fconv)
             dt0=deltz(i,j)*deltz(i,j)/(4.*kz(i,j))
             iii= dt/dt0+1
             dtt=dt/iii
             do it=1,iii
                 delta01=M2u*c0(i,j)
                 delta02=-M2di*c(i,j)
                 delta03=M2dit*cp1(i,j)*deltzt(i,j)/deltz(i,j)
                 delta04= (kzz*(cp1(i,j)-c(i,j))/deltz(i,j))/deltz(i,j)
                 delta05= -(kzzb*(c(i,j)-c1(i,j))/deltz(i,j))/deltz(i,j)
              if(k==1) then                 
                c(i,j)=c(i,j)+dtt*(delta04+delta05+M2dit*(cp1(i,j)-c(i,j)))
              else if((heizt(i,j)-TERRAIN(i,j))>pbl(i,j).and.(heiz(i,j)-TERRAIN(i,j))<pbl(i,j)) then
                c(i,j)=c(i,j)+dtt*(delta04+delta05+delta01+delta02)
              else      
                c(i,j)=c(i,j)+dtt*(delta01+delta02+delta03+delta04+delta05)
              endif
             enddo  !it
              c(i,j)=max(0.,c(i,j))
          enddo !i
        enddo  !j
       return
       end   subroutine        
       
          subroutine diffu_vert(myid,c1,c_1,c,cp1,&
              c1tmp,c_1tmp,ctmp,cp1tmp,&              
              kv1,kv_1,kv,kvp1,h1,h_1,h,hp1,terrain,&
                      dep1,dep_1,dep,depp1,&
                      fcon,pbl,kpbl,Gc_MOLWT,ig,igas,sx,ex,sy,ey,dt,k,ne,ii)
            integer  :: myid, sx, ex, sy, ey,it 
            real, dimension(sx-1:ex+1,sy-1:ey+1) :: c1,c_1,c,cp1
            real, dimension(sx-1:ex+1,sy-1:ey+1) :: c1tmp,c_1tmp,ctmp,cp1tmp
            real, dimension(sx-1:ex+1,sy-1:ey+1) :: kv1,kv_1,kv,kvp1
            real, dimension(sx-1:ex+1,sy-1:ey+1) :: h1,h_1,h,hp1
            real, dimension(sx-1:ex+1,sy-1:ey+1) :: dep1,dep_1,dep,depp1
            real, dimension(sx-1:ex+1,sy-1:ey+1) :: fcon,pbl,terrain
            integer, dimension(sx-1:ex+1,sy-1:ey+1) :: kpbl
            real :: M2u,M2di,M2dii,cm1,cm_1,cm,cmp1,cm1tmp,cm_1tmp,cmtmp,cmp1tmp
            real :: dt,dtt,dt0,delta1,delta2,delta3,delta4,delta5
            real :: GC_MOLWT(igas)
           do j = sy,ey
           do i = sx,ex
           
            M2u= fcon(i,j)*kv1(i,j)/dep1(i,j)/(pbl(i,j)-h1(i,j)+terrain(i,j)-dep1(i,j)/2.)
            M2di=M2u*(pbl(i,j)-h(i,j)+terrain(i,j)+dep(i,j)/2.)/dep(i,j)
            M2dii=M2u*(pbl(i,j)-hp1(i,j)+terrain(i,j)+depp1(i,j)/2.)/depp1(i,j)
           

            dtt=dt
!------convert ppbv --> mass mixing ratio------
            cm1=c1(i,j)*GC_MOLWT(ig)/29. !29 is mole mass of air
            cm_1=c_1(i,j)*GC_MOLWT(ig)/29.
            cm=c(i,j)*GC_MOLWT(ig)/29.
            cmp1=cp1(i,j)*GC_MOLWT(ig)/29.

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
            delta4= 2.*(1-fcon(i,j))*kvp1(i,j)*(cmp1tmp-cmtmp)&
                   /(dep(i,j)+depp1(i,j))/dep(i,j)
            delta5= -2.*(1-fcon(i,j))*kv_1(i,j)*(cmtmp-cm_1tmp)& 
                  /(dep(i,j)+dep_1(i,j))/dep(i,j)
           else
            delta4= 2.*kvp1(i,j)*(cmp1tmp-cmtmp)&
                        /(dep(i,j)+dep(i,j))/dep(i,j)
            delta5= -2.*kv_1(i,j)*(cmtmp-cm_1tmp)&
                        /(dep(i,j)+dep(i,j))/dep(i,j)       
           endif
              cm= cm+(delta1+delta2+delta3+delta4+delta5)*dtt
!------------convert mass mixing ratio --> ppbv------
            c(i,j)= cm*29./GC_MOLWT(ig)        
                
           enddo !i
           enddo !j
            return 
            end subroutine

      SUBROUTINE diffus (myid, ig, sx, ex, sy, ey, nlay, deltat, rkv,   &
      depth, tempk, press, conc, atm, kktop, gc_molwt)                  
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
!                                                                       
!     Output arguments:                                                 
!        conc              species concentrations (umol/m3)             
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
      DIMENSION conc (sx - 1:ex + 1, sy - 1:ey + 1, nlay), vdep (sx - 1:&
      ex + 1, sy - 1:ey + 1), kktop (sx - 1:ex + 1, sy - 1:ey + 1)      
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
      REAL convfac, gc_molwt 
                !lijie                                                  
      mapscl = 1. 
                ! here no dry depositions                               
      vdep = 0.0 
      nddmsp = 0 
!-----------convert conc unit, by chenhs ---                            
      DO i = sx, ex 
      DO j = sy, ey 
      DO k = 1, nlay 
                            !! ppb to umol/m3                           
      IF (ig.le.74) then 
         convfac = 44.9 * (273. / tempk (i, j, k) ) * (press (i, j, k)  &
         / 1013.)                                                       
         conc (i, j, k) = conc (i, j, k) * convfac / 1e+03 
                                              !! ug/m3 to umol/m3       
      ELSEIF (ig.gt.74.and.ig.le.102) then 
         conc (i, j, k) = conc (i, j, k) / gc_molwt 
               !! ng/m3 to umol/m3                                      
      ELSE 
         conc (i, j, k) = conc (i, j, k) * 1e-03 / gc_molwt 
      ENDIF 
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
                                                                        
            DO k = 1, nlay 
            ro1d (k) = press (i, j, k) / tempk (i, j, k) 
            d1d (k) = depth (i, j, k) 
            c1d (k) = conc (i, j, k) 
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
            CALL vdiffimp (nlay, dtv, vdep (i, j), d1d, ro1d, rk1d, c1d,&
            nddmsp, fcup, fcdn, ldoipts)                                
!                                                                       
!                                                                       
!                                                                       
            DO k = 1, nlay 
            conc (i, j, k) = c1d (k) 
            enddo 
                                                                        
                  ! nstepv                                              
            enddo 
   50    END DO 
   60 END DO 
!                                                                       
!-----------convert unit, by chenhs---                                  
      DO i = sx, ex 
      DO j = sy, ey 
      DO k = 1, nlay 
                              !! umol/m3 to ppb                         
      IF (ig.le.74) then 
         convfac = 44.9 * (273. / tempk (i, j, k) ) * (press (i, j, k)  &
         / 1013.)                                                       
         conc (i, j, k) = 1e+03 * conc (i, j, k) / convfac 
                                                !! umol/m3 to ug/m3     
      ELSEIF (ig.gt.74.and.ig.le.102) then 
         conc (i, j, k) = conc (i, j, k) * gc_molwt 
                 !! umol/m3 to ng/m3                                    
      ELSE 
         conc (i, j, k) = conc (i, j, k) * gc_molwt * 1e+03 
      ENDIF 
      conc (i, j, k) = amax1 (conc (i, j, k), 0.0) 
      enddo 
      enddo 
      enddo 
      RETURN 
      END SUBROUTINE diffus                         

      SUBROUTINE diffus_ds_com (myid, ig, sx, ex, sy, ey, nlay, deltat, &
      rkv, depth, tempk, press, conc, concom, atm, kktop, nspecom, IS,  &
      IA)                                                               
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
!                                                                       
!     Output arguments:                                                 
!        conc              species concentrations (umol/m3)             
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
      INTEGER ::sx, ex, sy, ey, nspecom 
      INTEGER ::MXTRSP 
                          ! the sensitity here no use                   
      PARAMETER (MXTRSP = 0) 
      DIMENSION conc (sx - 1:ex + 1, sy - 1:ey + 1, nlay), vdep (sx - 1:&
      ex + 1, sy - 1:ey + 1), kktop (sx - 1:ex + 1, sy - 1:ey + 1),     &
      concom (sx - 1:ex + 1, sy - 1:ey + 1, nlay, nspecom)              
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
      REAL ratio 
                !lijie                                                  
      mapscl = 1. 
                ! here no dry depositions                               
      vdep = 0.0 
      nddmsp = 0 
!-----------convert ppbv to umol/m3---                                  
!      do i=sx,ex                                                       
!       do j=sy,ey                                                      
!        do k=1,nlay                                                    
                                                                        
!          convfac=44.9*(273./tempk(i,j,k))*(press(i,j,k)/1013.)        
!          conc(i,j,k)=conc(i,j,k)*convfac/1e+03                        
                                                                        
!         do iduc = 1, nspecom                                          
!          concom(i,j,k,iduc) = concom(i,j,k,iduc) *convfac/1e+03       
!         enddo                                                         
                                                                        
!        enddo                                                          
!       enddo                                                           
!      enddo                                                            
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
  401       nstepv (i, j) = INT (0.999 * deltat / dtv) + 10 
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
                                                                        
            DO k = 1, nlay 
            ro1d (k) = press (i, j, k) / tempk (i, j, k) 
            d1d (k) = depth (i, j, k) 
            c1d (k) = conc (i, j, k) 
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
            CALL vdiffimp_ds (nlay, dtv, vdep (i, j), d1d, ro1d, rk1d,  &
            c1d, nddmsp, fcup, fcdn, ldoipts)                           
!                                                                       
!                                                                       
!                                                                       
            DO k = 1, nlay - 1 
                                                                        
            IF (c1d (k) - conc (i, j, k) .ne.0) then 
               ratio = (fcup (k) + fcdn (k) ) / (c1d (k) - conc (i, j,  &
               k) )                                                     
            ELSE 
               ratio = 1.0 
            ENDIF 
                                                                        
            IF (ratio.eq. 0.and. (c1d (k) - conc (i, j, k) ) .ne.0.)    &
            THEN                                                        
               fcup (k) = c1d (k) - conc (i, j, k) 
               fcdn (k) = 0.0 
               ratio = 1. 
            ENDIF 
                                                                        
            fcup (k) = fcup (k) / ratio 
            fcdn (k) = fcdn (k) / ratio 
                                                                        
            conc_tmp0 = conc (i, j, k) 
            conc_tmp1 = conc (i, j, k) + fcup (k) 
            conc_tmp2 = conc (i, j, k) + fcup (k) + fcdn (k) 
                                                                        
            conc (i, j, k) = amax1 (conc (i, j, k), 1.E-20) 
                                                                        
                                                                        
            DO iduc = 1, nspecom 
            concom (i, j, k, iduc) = amax1 (concom (i, j, k, iduc),     &
            1.E-20)                                                     
                                                                        
                                                                        
            IF (fcup (k) .GT.0.) THEN 
               IF (fcup (k) >conc (i, j, k + 1) ) GOTO 100 
               conc (i, j, k + 1) = amax1 (conc (i, j, k + 1), 1.E-20) 
               concom (i, j, k + 1, iduc) = amax1 (concom (i, j, k + 1, &
               iduc), 1.E-20)                                           
                                                                        
               concom (i, j, k, iduc) = concom (i, j, k, iduc) + fcup ( &
               k) * concom (i, j, k + 1, iduc) / conc (i, j, k + 1)     
                                                                        
            ELSE 
               IF (ABS (fcup (k) ) >conc_tmp0) GOTO 100 
               concom (i, j, k, iduc) = concom (i, j, k, iduc) + fcup ( &
               k) * concom (i, j, k, iduc) / conc_tmp0                  
                                                                        
                                                                        
            ENDIF 
                                                                        
  100       CONTINUE 
            IF (fcdn (k) .GT.0.) THEN 
               IF (K .eq. 1) THEN 
                  conc_1 = conc (i, j, 1) 
                  concom_1 = concom (i, j, 1, iduc) 
               ELSE 
                  conc_1 = conc (i, j, k - 1) 
                  concom_1 = concom (i, j, k - 1, iduc) 
               ENDIF 
                                                                        
               IF (fcdn (k) >conc_1) GOTO 101 
               conc_1 = amax1 (conc_1, 1.E-20) 
               concom_1 = amax1 (concom_1, 1.E-20) 
               concom (i, j, k, iduc) = concom (i, j, k, iduc) + fcdn ( &
               k) * concom_1 / conc_1                                   
            ELSE 
                                                                        
               IF (ABS (fcdn (k) ) >conc_tmp0) GOTO 101 
               concom (i, j, k, iduc) = concom (i, j, k, iduc) + fcdn ( &
               k) * concom (i, j, k, iduc) / conc_tmp1                  
            ENDIF 
            concom (i, j, k, iduc) = amax1 (concom (i, j, k, iduc),     &
            1.E-21)                                                     
  101       CONTINUE 
                                                                        
                    ! iduc                                              
            enddo 
                                                                        
                  ! k                                                   
            enddo 
                                                                        
            DO k = 1, nlay 
            conc (i, j, k) = c1d (k) 
            conc (i, j, k) = amax1 (conc (i, j, k), 1.E-20) 
                    ! k                                                 
            enddo 
                                                                        
                  ! nstepv                                              
            enddo 
   50    END DO 
   60 END DO 
!                                                                       
!-----------convert umol/m3 to ppb---                                   
!      do i=sx,ex                                                       
!       do j=sy,ey                                                      
!         do k=1,nlay                                                   
!          convfac=44.9*(273./tempk(i,j,k))*(press(i,j,k)/1013.)        
!          conc(i,j,k)=1e03*conc(i,j,k)/convfac                         
!          conc(i,j,k)=amax1(conc(i,j,k),0.0)                           
                                                                        
!          do iduc = 1, nspecom                                         
!           concom (i,j,k,iduc) = 1e03*concom (i,j,k,iduc)/convfac      
!          enddo                                                        
                                                                        
!         enddo                                                         
!        enddo                                                          
!      enddo                                                            
      RETURN 
      END SUBROUTINE diffus_ds_com                  
                                                                        
      SUBROUTINE diffus_ds (myid, ig, sx, ex, sy, ey, nlay, deltat, rkv,&
      depth, tempk, press, conc, atm, kktop)                            
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
!                                                                       
!     Output arguments:                                                 
!        conc              species concentrations (umol/m3)             
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
                            !,nspecom                                   
      INTEGER ::sx, ex, sy, ey 
      INTEGER ::MXTRSP 
                          ! the sensitity here no use                   
      PARAMETER (MXTRSP = 0) 
      DIMENSION conc (sx - 1:ex + 1, sy - 1:ey + 1, nlay), vdep (sx - 1:&
      ex + 1, sy - 1:ey + 1), kktop (sx - 1:ex + 1, sy - 1:ey + 1)      
!     & ,concom(sx-1:ex+1,sy-1:ey+1,nlay,nspecom)                       
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
      REAL ratio 
                !lijie                                                  
      mapscl = 1. 
                ! here no dry depositions                               
      vdep = 0.0 
      nddmsp = 0 
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
  401       nstepv (i, j) = INT (0.999 * deltat / dtv) + 10 
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
                                                                        
            DO k = 1, nlay 
            ro1d (k) = press (i, j, k) / tempk (i, j, k) 
            d1d (k) = depth (i, j, k) 
            c1d (k) = conc (i, j, k) 
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
            CALL vdiffimp_ds (nlay, dtv, vdep (i, j), d1d, ro1d, rk1d,  &
            c1d, nddmsp, fcup, fcdn, ldoipts)                           
!                                                                       
!                                                                       
!                                                                       
           !do k = 1,nlay-1                                             
                                                                        
           !   if(c1d(k)-conc(i,j,k).ne.0) then                         
           !    ratio = (fcup(k)+fcdn(k))/(c1d(k)-conc(i,j,k))          
           !   else                                                     
           !    ratio = 1.0                                             
           !   endif                                                    
                                                                        
           !    IF(ratio==0.and. (c1d(k)-conc(i,j,k)).ne.0.) THEN       
           !     fcup(k) = c1d(k)-conc(i,j,k)                           
           !     fcdn(k) = 0.0                                          
           !     ratio = 1.                                             
           !    ENDIF                                                   
                                                                        
           !   fcup(k)= fcup(k) /ratio                                  
           !   fcdn(k)= fcdn(k) / ratio                                 
                                                                        
           !   conc_tmp0 = conc(i,j,k)                                  
           !   conc_tmp1 = conc(i,j,k) + fcup(k)                        
           !   conc_tmp2 = conc(i,j,k) + fcup(k) + fcdn(k)              
                                                                        
           !   conc(i,j,k) = amax1(conc(i,j,k),1.E-20)                  
                                                                        
                                                                        
              !do iduc = 1,nspecom                                      
              !  concom(i,j,k,iduc) = amax1(concom(i,j,k,iduc),1.E-20)  
                                                                        
                                                                        
              ! IF(fcup(k).GT.0.) THEN                                  
              !   IF(fcup(k)>conc(i,j,k+1)) GOTO 100                    
              !   conc(i,j,k+1) = amax1(conc(i,j,k+1),1.E-20)           
              !   concom(i,j,k+1,iduc)=amax1(concom(i,j,k+1,iduc),1.E-20
                                                                        
              !  concom(i,j,k,iduc) = concom(i,j,k,iduc) +              
              !&               fcup(k)*concom(i,j,k+1,iduc)/conc(i,j,k+1
                                                                        
              !  ELSE                                                   
              !    IF(ABS(fcup(k))>conc_tmp0 ) GOTO 100                 
              !  concom(i,j,k,iduc) = concom(i,j,k,iduc) +              
              !&                 fcup(k)*concom(i,j,k,iduc)/conc_tmp0   
                                                                        
                                                                        
              ! ENDIF                                                   
                                                                        
!100   CONTINUE                                                         
              ! IF(fcdn(k).GT.0.) THEN                                  
              !  if(K==1) THEN                                          
              !    conc_1=conc(i,j,1)                                   
              !    concom_1 = concom(i,j,1,iduc)                        
              !  else                                                   
              !    conc_1 = conc(i,j,k-1)                               
              !    concom_1 = concom(i,j,k-1,iduc)                      
              !  endif                                                  
                                                                        
              !  IF(fcdn(k)>conc_1 ) GOTO 101                           
              !   conc_1 = amax1(conc_1,1.E-20)                         
              !   concom_1 = amax1(concom_1,1.E-20)                     
              !   concom(i,j,k,iduc) = concom(i,j,k,iduc) +             
              !&         fcdn(k)*concom_1/conc_1                        
              ! ELSE                                                    
                                                                        
              !   IF(ABS(fcdn(k))>conc_tmp0) GOTO 101                   
              !   concom(i,j,k,iduc) = concom(i,j,k,iduc) +             
              !&                 fcdn(k)*concom(i,j,k,iduc)/conc_tmp1   
              ! ENDIF                                                   
              !   concom(i,j,k,iduc) = amax1(concom(i,j,k,iduc),1.E-21) 
!101    CONTINUE                                                        
                                                                        
              !enddo ! iduc                                             
            !enddo ! k                                                  
                                                                        
            DO k = 1, nlay 
            conc (i, j, k) = c1d (k) 
            conc (i, j, k) = amax1 (conc (i, j, k), 1.E-20) 
                    ! k                                                 
            enddo 
                                                                        
                  ! nstepv                                              
            enddo 
   50    END DO 
   60 END DO 
!                                                                       
      RETURN 
      END SUBROUTINE diffus_ds
      
      SUBROUTINE vdiffimp (nn, dt, vdep, depth, rho, rkv, rr, nparddm,  &
      fcup, fcdn, ldoipts)                                              
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
      DIMENSION aa (nn), bb (nn), cc (nn) 
!                                                                       
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
!                                                                       
!                                                                       
!-----Solve the equations                                               
!                                                                       
      CALL trdiag (aa, bb, cc, rr, nn, 1 + nparddm) 
                                                                        
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
      RETURN 
      END SUBROUTINE vdiffimp                       

      SUBROUTINE vdiffimp_ds (nn, dt, vdep, depth, rho, rkv, rr,        &
      nparddm, fcup, fcdn, ldoipts)                                     
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
      DIMENSION aa (nn), bb (nn), cc (nn) 
!                                                                       
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
!                                                                       
!                                                                       
!-----Solve the equations                                               
!                                                                       
      CALL trdiag_ds (aa, bb, cc, rr, nn, 1 + nparddm) 
                                                                        
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
      RETURN 
      END SUBROUTINE vdiffimp_ds

        SUBROUTINE trdiag (a, b, c, d, n, m) 
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
      END SUBROUTINE trdiag                         
      
      SUBROUTINE trdiag_ds (a, b, c, d, n, m) 
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
      END SUBROUTINE trdiag_ds                      
      
end module diffusion
