!*
!* do column level mapping between two restart files
!* if column index same, copy all its column level and pft level fields to new file
!* one non-vegetation column has one patch to copy
!* one vegetation column has numPFT(16pfts+1soil=17) patches to copy
!*
!* Weights of column in new file are recalculated.
!*
!* Assume the source restart file has more columns than the destination file
!*

PROGRAM MAIN

   use precision
   use paramodel, only: nfcon_col, nfvar_col, nfcon_pft, nfvar_pft, nftune
   implicit none

   integer, parameter :: lat_points =  64
   integer, parameter :: lon_points = 128
   integer, parameter :: numcolumn1 = 5253
   integer, parameter :: numpatch1  = 45221
   integer, parameter :: numcolumn2 = 5136
   integer, parameter :: numpatch2  = 44400

   character(len=255) :: fgrid1    = "../data/CoLM_T42_MOM4p1_grid_pft"
   character(len=255) :: fgrid2    = "../data/CoLM_T42_MOM4p1_grid_pft_glc"
   character(len=255) :: fconst1   = "../data/CoLM_T42_MOM4p1_const_pft"
   character(len=255) :: fconst2   = "../data/CoLM_T42_MOM4p1_const_pft_glc"
   character(len=255) :: frestart1 = "../data/T42-restart-1382"
   character(len=255) :: frestart2 = "../data/CoLM_T42_MOM4p1_restart_pft_glc"

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

   integer  itypwat1(numcolumn1)
   integer  ixy_column1(numcolumn1)
   integer  jxy_column1(numcolumn1)
   real(r8) wxy_column1(numcolumn1)
   integer  ixy_patch1(numpatch1)
   integer  jxy_patch1(numpatch1)
   real(r8) wxy_patch1(numpatch1)

   integer  itypwat2(numcolumn2)
   integer  ixy_column2(numcolumn2)
   integer  jxy_column2(numcolumn2)
   real(r8) wxy_column2(numcolumn2)
   integer  ixy_patch2(numpatch2)
   integer  jxy_patch2(numpatch2)
   real(r8) wxy_patch2(numpatch2)

   real(r8) fcon_col1(nfcon_col,numcolumn1)
   real(r8) fvar_col1(nfvar_col,numcolumn1)
   real(r8) fcon_pft1(nfcon_pft,numpatch1)
   real(r8) fvar_pft1(nfvar_pft,numpatch1)

   real(r8) fcon_col2(nfcon_col,numcolumn2)
   real(r8) fvar_col2(nfvar_col,numcolumn2)
   real(r8) fcon_pft2(nfcon_pft,numpatch2)
   real(r8) fvar_pft2(nfvar_pft,numpatch2)

   logical g1mask(lon_points,lat_points)
   logical g2mask(lon_points,lat_points)

   integer g1begc(lon_points,lat_points), g1endc(lon_points,lat_points)
   integer g2begc(lon_points,lat_points), g2endc(lon_points,lat_points)

   integer c1begp(numcolumn1), c1endp(numcolumn1)
   integer c2begp(numcolumn2), c2endp(numcolumn2)

   integer ccmap(numcolumn2)
   integer ppmap(numpatch2)

   integer istep, idate_p(3), idate(3)

   real(r8), pointer :: rbuf(:)     ! buffer for reading fvar/fcon

   integer i, j, p, c, g, c1, c2, p1, p2, nump

   integer c1b, c1e, c2b, c2e

! routine:

   call rdgrid(fgrid1,lon_points,lat_points,numcolumn1,numpatch1, &
               numlon,lats,lonw,area,latixy,longxy,landfrac,landmask, &
               itypwat1,ixy_column1,jxy_column1,wxy_column1,ixy_patch1,jxy_patch1,wxy_patch1)

   call rdgrid(fgrid2,lon_points,lat_points,numcolumn2,numpatch2, &
               numlon,lats,lonw,area,latixy,longxy,landfrac,landmask, &
               itypwat2,ixy_column2,jxy_column2,wxy_column2,ixy_patch2,jxy_patch2,wxy_patch2)

   call rdinit(fconst1,frestart1,numcolumn1,numpatch1,istep,idate,idate_p,ftune, &
               fcon_col1,fvar_col1,fcon_pft1,fvar_pft1)

   call rdinit(fconst2,frestart2,numcolumn2,numpatch2,istep,idate,idate_p,ftune, &
               fcon_col2,fvar_col2,fcon_pft2,fvar_pft2)

   g1begc = -1
   g1endc = -1
   g1mask = .false.

   p = 1

   do c = 1, numcolumn1 
      i = ixy_column1(c)
      j = jxy_column1(c)
      g1mask(i,j) = .true.
      if(g1begc(i,j).lt.0) g1begc(i,j) = c
      if(g1begc(i,j).gt.0) g1endc(i,j) = c

      if(itypwat1(c).eq.0) then
         nump = 17
      else if(itypwat1(c).ge.1 .and. itypwat1(c).le.4)  then
         nump = 1
      else
         write(6,*) 'impossible itypwat1', itypwat1(c)
         call abort
      end if

      c1begp(c) = p
      c1endp(c) = p + nump - 1

      p = p + nump
   end do

   g2begc = -1
   g2endc = -1
   g2mask = .false.

   p = 1

   do c = 1, numcolumn2
      i = ixy_column2(c)
      j = jxy_column2(c)
      g2mask(i,j) = .true.
      if(g2begc(i,j).lt.0) g2begc(i,j) = c
      if(g2begc(i,j).gt.0) g2endc(i,j) = c

      if(itypwat2(c).eq.0) then
         nump = 17
      else if(itypwat2(c).ge.1 .and. itypwat2(c).le.4)  then
         nump = 1
      else
         write(6,*) 'impossible itypwat2', itypwat2(c)
         call abort
      end if

      c2begp(c) = p
      c2endp(c) = p + nump - 1

      p = p + nump
   end do

!* building ccmap firstly
!* following copies will base on ccmap and itypwat

   ccmap = -1

   do j = 1, lat_points
   do i = 1, lon_points
      if(g1mask(i,j) .and. g2mask(i,j)) then
         do c2 = g2begc(i,j), g2endc(i,j)
            do c1 = g1begc(i,j), g1endc(i,j)
               if(itypwat2(c2).eq.itypwat1(c1)) then
                  ccmap(c2) = c1
                  exit
               end if
            end do

            if(ccmap(c2).lt.0) then
               write(6,*) "+++++++++++ unmatch column +++++++++++"
               write(6,*) "columns in g2", g2begc(i,j), g2endc(i,j)
               write(6,*) "columns in g1", g1begc(i,j), g1endc(i,j)
               write(6,*) "itypwats in g2", itypwat2(g2begc(i,j):g2endc(i,j))
               write(6,*) "itypwats in g1", itypwat1(g1begc(i,j):g1endc(i,j))
            end if
         end do
      else if(g1mask(i,j) .and. (.not. g2mask(i,j))) then
         write(6,*) 'pass grid', i, j
      else if(g2mask(i,j) .and. (.not. g1mask(i,j))) then
         write(6,*) 'wrong grid', i, j
         call abort
      end if
   end do
   end do

!* building ppmap secondly

   p = 1

   ppmap = -1

   do c2 = 1, numcolumn2
      if(itypwat2(c2).eq.0) then
         nump = 17
      else if(itypwat2(c2).ge.1 .and. itypwat2(c2).le.4) then
         nump = 1
      else 
         write(6,*) 'impossible itypwat2', itypwat2(c2)
         call abort
      end if

      c1 = ccmap(c2)

      if(c1.gt.0) then
         if(itypwat1(c1).ne.itypwat2(c2)) then
            write(6,*) 'wrong ccmap on itypwat', c1, c2, itypwat1(c1), itypwat2(c2)
            call abort
         end if
   
         if((c1endp(c1)-c1begp(c1)+1).ne.nump) then
            write(6,*) 'wrong c1begp, c1endp, nump relationship', c1begp(c1), c1endp(c1), nump
            call abort
         end if
   
         if((c2endp(c2)-c2begp(c2)+1).ne.nump) then
            write(6,*) 'wrong c2begp, c2endp, nump relationship', c2begp(c2), c2endp(c2), nump
            call abort
         end if
   
       ! ppmap(p:p+nump-1) = (/(i,i=c1begp(c1),c1endp(c1))/)
         ppmap(c2begp(c2):c2endp(c2)) = (/(i,i=c1begp(c1),c1endp(c1))/)
      end if

      p = p + nump

      if(p.ne.c2endp(c2)+1) then
         write(6,*), 'error in building ppmap', p, c2endp(c2)
         call abort
      end if
   end do

!* do column level copying

   do c2 = 1, numcolumn2
      c1 = ccmap(c2)

      if(c1.gt.0) then
         fcon_col2(:,c2) = fcon_col1(:,c1)
         fvar_col2(:,c2) = fvar_col1(:,c1)
      end if
   end do

!* do patch level copying

   do p2 = 1, numpatch2
      p1 = ppmap(p2)

      if(p1.gt.0) then
         fcon_pft2(:,p2) = fcon_pft1(:,p1)
         fvar_pft2(:,p2) = fvar_pft1(:,p1)
      end if
   end do

   fconst2 = trim(fconst2)//"_updated"
   frestart2 = trim(frestart2)//"_updated"

   call wrinit(fconst2,frestart2,numcolumn2,numpatch2,istep,idate,idate_p,ftune,fcon_col2,fvar_col2,fcon_pft2,fvar_pft2)

END PROGRAM MAIN

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
   OPEN(unit=luconst,file=fconst,form='unformatted',status='old',action='read')

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

   OPEN(unit=lurestart,file=frestart,form='unformatted',status='old',action='read')

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

   OPEN(lugrid,file=fgrid,form='unformatted',status='old')
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
   OPEN(unit=luconst,file=fconst,form='unformatted',status='unknown',action='write')

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

   OPEN(unit=lurestart,file=frestart,form='unformatted',status='unknown',action='write')

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
