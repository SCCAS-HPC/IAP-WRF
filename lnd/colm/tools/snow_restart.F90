!*
!* reduce column level snow mass: wliq, wice, scv, snowdp, dz, z
!*

PROGRAM MAIN

   use precision
   use paramodel, only: nfcon_col, nfvar_col, nfcon_pft, nfvar_pft, nftune, nl_soil, maxsnl
   implicit none

   real(r8),parameter :: factor = 0.01
   real(r8),parameter :: scv_threshold = 1.E5

   integer, parameter :: lat_points =  64
   integer, parameter :: lon_points = 128
   integer, parameter :: numcolumn = 5136
   integer, parameter :: numpatch  = 44400

   character(len=255) :: fgrid    = "/mnt/swgfs/jidy/piCtlMPI/lndPIC_1950/lnd/CoLM-T42-gridata-c"
   character(len=255) :: fconst   = "/mnt/swgfs/jidy/piCtlMPI/lndPIC_1950/lnd/CoLM-T42-const-c"
   character(len=255) :: frestart = "/mnt/swgfs/jidy/piCtlMPI/lndPIC_1950/lnd/CoLM-T42-restart-1949"
   character(len=255) :: fsbc     = "/mnt/swgfs/jidy/piCtlMPI/lndPIC_1950/lnd/CoLM-T42-restart-1949-sbc"

! local variables:

   integer  numlon(lat_points)
   real(r8) lats(lat_points+1)
   real(r8) lonw(lon_points+1,lat_points)
   real(r8) area(lon_points,lat_points)
   real(r8) latixy(lon_points,lat_points)
   real(r8) longxy(lon_points,lat_points)
   real(r8) landfrac(lon_points,lat_points)
   integer  landmask(lon_points,lat_points)
   real(r8) ftune(nftune)

   real(r8) avsdr(lon_points,lat_points)
   real(r8) avsdf(lon_points,lat_points)
   real(r8) anidr(lon_points,lat_points)
   real(r8) anidf(lon_points,lat_points)
   real(r8) trad (lon_points,lat_points)
   real(r8) xyscv(lon_points,lat_points)

   integer  itypwat(numcolumn)
   integer  ixy_column(numcolumn)
   integer  jxy_column(numcolumn)
   real(r8) wxy_column(numcolumn)
   integer  ixy_patch(numpatch)
   integer  jxy_patch(numpatch)
   real(r8) wxy_patch(numpatch)

   real(r8), target :: fcon_col(nfcon_col,numcolumn)
   real(r8), target :: fvar_col(nfvar_col,numcolumn)
   real(r8), target :: fcon_pft(nfcon_pft,numpatch)
   real(r8), target :: fvar_pft(nfvar_pft,numpatch)

   integer istep, idate_p(3), idate(3)

   real(r8), pointer :: rbuf(:)          ! buffer for reading fvar/fcon

   real(r8), pointer :: z   (:)          ! node depth [m] 
   real(r8), pointer :: dz  (:)          ! interface depth [m] 
   real(r8), pointer :: tss (:)          ! soil temperature [K] 
   real(r8), pointer :: wliq(:)          ! liquid water in layers [kg/m2]
   real(r8), pointer :: wice(:)          ! ice lens in layers [kg/m2]
   real(r8), pointer :: tg               ! ground surface temperature [K] 
   real(r8), pointer :: sag              ! non dimensional snow age [-] 
   real(r8), pointer :: scv              ! snow cover, water equivalent [mm]
   real(r8), pointer :: snowdp           ! snow depth [meter]
   real(r8), pointer :: coszen           ! cosine of solar zenith angle
   real(r8), pointer :: fsno             ! fraction of snow cover on ground

   integer i, j, p, c, g, jm, lc, uc, k, L
   real(r8) zi

! routine:

   call rdgrid(fgrid,lon_points,lat_points,numcolumn,numpatch, &
               numlon,lats,lonw,area,latixy,longxy,landfrac,landmask, &
               itypwat,ixy_column,jxy_column,wxy_column,ixy_patch,jxy_patch,wxy_patch)

   call rdinit(fconst,frestart,numcolumn,numpatch,istep,idate,idate_p,ftune, &
               fcon_col,fvar_col,fcon_pft,fvar_pft)

   jm = nl_soil + abs(maxsnl)

   do c = 1, numcolumn 
      i = ixy_column(c)
      j = jxy_column(c)

      lc = 1
      uc = jm

      z         => fvar_col(lc:uc,c)                  ; lc = uc+1; uc = uc+jm !1_ 
      dz        => fvar_col(lc:uc,c)                  ; lc = uc+1; uc = uc+jm !2_ 
      tss       => fvar_col(lc:uc,c)                  ; lc = uc+1; uc = uc+jm !3_ 
      wliq      => fvar_col(lc:uc,c)                  ; lc = uc+1; uc = uc+jm !4_ 
      wice      => fvar_col(lc:uc,c)                  ; uc = uc+1             !5_ 
      tg        => fvar_col(uc,c)                     ; uc = uc+1             !1  
      sag       => fvar_col(uc,c)                     ; uc = uc+1             !2  
      scv       => fvar_col(uc,c)                     ; uc = uc+1             !3  
      snowdp    => fvar_col(uc,c)                     ; uc = uc+1             !4  
      fsno      => fvar_col(uc,c)                     ; uc = uc+1             !5  
      coszen    => fvar_col(uc,c)                                             !6  

      if(scv.gt.scv_threshold) then
#ifdef MYBUG
         print *, 'scv', scv
         print *, 'snowdp', snowdp
         print *, 'dz', dz
         print *, 'z', z
         print *, 'wliq', wliq
         print *, 'wice', wice
#endif

         L = abs(maxsnl)

         dz(L) = dz(L)*factor
         snowdp = sum(dz(1:L))

         zi = 0.
         do k = L, 1, -1
            z(k) = zi-dz(k)*0.5
            zi = zi-dz(k)
         end do

         wliq(L) = wliq(L)*factor
         wice(L) = wice(L)*factor
         scv = sum(wice(1:L))

#ifdef MYBUG
         print *, 'scv', scv
         print *, 'snowdp', snowdp
         print *, 'dz', dz
         print *, 'z', z
         print *, 'wliq', wliq
         print *, 'wice', wice
         exit
#endif
      end if
   end do

   call rdsbc(fsbc,lon_points,lat_points,avsdr,avsdf,anidr,anidf,trad,xyscv)

   do j = 1, lat_points
   do i = 1, lon_points
      if (xyscv(i,j).gt.scv_threshold) then
          xyscv(i,j) = xyscv(i,j)*factor
      end if
   end do
   end do

   fconst = trim(fconst)//"x"
   frestart = trim(frestart)//"x"
   fsbc = trim(fsbc)//"x"

   call wrinit(fconst,frestart,numcolumn,numpatch,istep,idate,idate_p,ftune,fcon_col,fvar_col,fcon_pft,fvar_pft)
   
   call wtsbc(fsbc,lon_points,lat_points,avsdr,avsdf,anidr,anidf,trad,xyscv)

END PROGRAM MAIN

subroutine rdsbc(fsbc,lon_points,lat_points,avsdr,avsdf,anidr,anidf,trad,scv)

   use precision
   implicit none

   character(len=255), intent(in) :: fsbc
   integer, intent(in)  :: lon_points
   integer, intent(in)  :: lat_points
   real(r8),intent(out) :: avsdr(lon_points,lat_points)
   real(r8),intent(out) :: avsdf(lon_points,lat_points)
   real(r8),intent(out) :: anidr(lon_points,lat_points)
   real(r8),intent(out) :: anidf(lon_points,lat_points)
   real(r8),intent(out) :: trad (lon_points,lat_points)
   real(r8),intent(out) :: scv  (lon_points,lat_points)

   integer :: lusbc = 240

   OPEN(unit=lusbc,file=trim(fsbc),form='unformatted',status='old',action='read')

   read(lusbc) avsdr
   read(lusbc) avsdf
   read(lusbc) anidr
   read(lusbc) anidf
   read(lusbc) trad
   read(lusbc) scv

   CLOSE(lusbc)

end subroutine rdsbc

subroutine wtsbc(fsbc,lon_points,lat_points,avsdr,avsdf,anidr,anidf,trad,scv)

   use precision
   implicit none

   character(len=255), intent(in) :: fsbc
   integer, intent(in)  :: lon_points
   integer, intent(in)  :: lat_points
   real(r8),intent(in) :: avsdr(lon_points,lat_points)
   real(r8),intent(in) :: avsdf(lon_points,lat_points)
   real(r8),intent(in) :: anidr(lon_points,lat_points)
   real(r8),intent(in) :: anidf(lon_points,lat_points)
   real(r8),intent(in) :: trad (lon_points,lat_points)
   real(r8),intent(in) :: scv  (lon_points,lat_points)

   integer :: lusbc = 240

   OPEN(unit=lusbc,file=trim(fsbc),form='unformatted',status='unknown',action='write')

   write(lusbc) avsdr
   write(lusbc) avsdf
   write(lusbc) anidr
   write(lusbc) anidf
   write(lusbc) trad
   write(lusbc) scv

   CLOSE(lusbc)

end subroutine wtsbc

subroutine rdinit(fconst,frestart,numcolumn,numpatch,istep,idate,idate_p,ftune,fcon_col,fvar_col,fcon_pft,fvar_pft)

   use precision
   use paramodel, only: nfcon_col, nfvar_col, nfcon_pft, nfvar_pft, nftune
   implicit none

   character(len=255), intent(in) :: fconst
   character(len=255), intent(in) :: frestart
   integer,  intent(in)  :: numcolumn
   integer,  intent(in)  :: numpatch
   real(r8), intent(out) :: fcon_col(nfcon_col,numcolumn)
   real(r8), intent(out) :: fvar_col(nfvar_col,numcolumn)
   real(r8), intent(out) :: fcon_pft(nfcon_pft,numpatch)
   real(r8), intent(out) :: fvar_pft(nfvar_pft,numpatch)
   real(r8), intent(out) :: ftune(nftune)
   integer,  intent(out) :: istep
   integer,  intent(out) :: idate(3)
   integer,  intent(out) :: idate_p(3)

   real(r8), pointer :: rbuf(:)     ! buffer for reading fvar/fcon

   integer :: luconst = 200
   integer :: lurestart = 220

   integer i, k

    ! Open for model time invariant constant data
   OPEN(unit=luconst,file=trim(fconst),form='unformatted',status='old',action='read')

   read(luconst) ftune       !clm tunable constants

   allocate (rbuf(numcolumn))

   do k = 1, nfcon_col
      read(luconst) rbuf
      fcon_col(k,:) = rbuf
   end do

   deallocate (rbuf)

   allocate (rbuf(numpatch))

   do k = 1, nfcon_pft
      read(luconst) rbuf
      fcon_pft(k,:) = rbuf
   end do

   deallocate (rbuf)

   CLOSE(luconst)

   OPEN(unit=lurestart,file=trim(frestart),form='unformatted',status='old',action='read')

   read(lurestart) istep, idate, idate_p

   allocate (rbuf(numcolumn))

   do k = 1, nfvar_col
      read(lurestart) rbuf
      fvar_col(k,:) = rbuf
   end do

   deallocate (rbuf)

   allocate (rbuf(numpatch))

   do k = 1, nfvar_pft
      read(lurestart) rbuf
      fvar_pft(k,:) = rbuf
   end do

   deallocate (rbuf)

   CLOSE(lurestart)

end subroutine rdinit

subroutine rdgrid(fgrid,lon_points,lat_points,numcolumn,numpatch, &
               numlon,lats,lonw,area,latixy,longxy,landfrac,landmask, &
               itypwat,ixy_column,jxy_column,wxy_column,ixy_patch,jxy_patch,wxy_patch)

   use precision
   use paramodel, only: nftune
   implicit none

   character(len=255), intent(in) :: fgrid
   integer , intent(in)  :: lon_points
   integer , intent(in)  :: lat_points
   integer , intent(in)  :: numcolumn
   integer , intent(in)  :: numpatch
   integer , intent(out) :: numlon(lat_points)
   real(r8), intent(out) :: lats(lat_points+1)
   real(r8), intent(out) :: lonw(lon_points+1,lat_points)
   real(r8), intent(out) :: area(lon_points,lat_points)
   real(r8), intent(out) :: latixy(lon_points,lat_points)
   real(r8), intent(out) :: longxy(lon_points,lat_points)
   real(r8), intent(out) :: landfrac(lon_points,lat_points)
   integer , intent(out) :: landmask(lon_points,lat_points)

   integer , intent(out) :: itypwat(numcolumn)
   integer , intent(out) :: ixy_column(numcolumn)
   integer , intent(out) :: jxy_column(numcolumn)
   real(r8), intent(out) :: wxy_column(numcolumn)
   integer , intent(out) :: ixy_patch(numpatch)
   integer , intent(out) :: jxy_patch(numpatch)
   real(r8), intent(out) :: wxy_patch(numpatch)

   integer :: lugrid = 100

   OPEN(lugrid,file=trim(fgrid),form='unformatted',status='old')
   read(lugrid) numlon(:)
   read(lugrid) lats(:)
   read(lugrid) lonw(:,:)
   read(lugrid) area(:,:)
   read(lugrid) latixy(:,:)
   read(lugrid) longxy(:,:)
   read(lugrid) landfrac(:,:)
   read(lugrid) landmask(:,:)
   read(lugrid) itypwat(:)
   read(lugrid) ixy_column(:) !longitude index for each patch point
   read(lugrid) jxy_column(:) !latitude index for each patch point
   read(lugrid) wxy_column(:) !subgrid weight for each patch point
   read(lugrid) ixy_patch(:)  !longitude index for each patch point
   read(lugrid) jxy_patch(:)  !latitude index for each patch point
   read(lugrid) wxy_patch(:)  !subgrid weight for each patch point
   CLOSE(lugrid)

end subroutine rdgrid

subroutine wrinit(fconst,frestart,numcolumn,numpatch,istep,idate,idate_p,ftune,fcon_col,fvar_col,fcon_pft,fvar_pft)

   use precision
   use paramodel, only: nfcon_col, nfvar_col, nfcon_pft, nfvar_pft, nftune
   implicit none

   character(len=255), intent(in) :: fconst
   character(len=255), intent(in) :: frestart
   integer,  intent(in) :: numcolumn
   integer,  intent(in) :: numpatch
   real(r8), intent(in) :: fcon_col(nfcon_col,numcolumn)
   real(r8), intent(in) :: fvar_col(nfvar_col,numcolumn)
   real(r8), intent(in) :: fcon_pft(nfcon_pft,numpatch)
   real(r8), intent(in) :: fvar_pft(nfvar_pft,numpatch)
   real(r8), intent(in) :: ftune(nftune)
   integer,  intent(in) :: istep
   integer,  intent(in) :: idate(3)
   integer,  intent(in) :: idate_p(3)

   real(r8), pointer :: rbuf(:)     ! buffer for reading fvar/fcon

   integer :: luconst = 200
   integer :: lurestart = 220

   integer i, k

 ! Open for model time invariant constant data
   OPEN(unit=luconst,file=trim(fconst),form='unformatted',status='unknown',action='write')

   write(luconst) ftune       !clm tunable constants

   allocate (rbuf(numcolumn))

   do k = 1, nfcon_col
      rbuf = fcon_col(k,:)
      write(luconst) rbuf
   end do

   deallocate (rbuf)

   allocate (rbuf(numpatch))

   do k = 1, nfcon_pft
      rbuf = fcon_pft(k,:)
      write(luconst) rbuf
   end do

   deallocate (rbuf)

   CLOSE(luconst)

   OPEN(unit=lurestart,file=trim(frestart),form='unformatted',status='unknown',action='write')

   write(lurestart) istep, idate, idate_p

   allocate (rbuf(numcolumn))

   do k = 1, nfvar_col
      rbuf = fvar_col(k,:)
      write(lurestart) rbuf
   end do

   deallocate (rbuf)

   allocate (rbuf(numpatch))

   do k = 1, nfvar_pft
      rbuf = fvar_pft(k,:)
      write(lurestart) rbuf
   end do

   deallocate (rbuf)

   CLOSE(lurestart)

end subroutine wrinit
