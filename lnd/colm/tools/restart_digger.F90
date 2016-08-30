!*
!* dig specified grid restart files from global dataset
!*

PROGRAM MAIN

   use precision
   use paramodel, only: nfcon_col, nfvar_col, nfcon_pft, nfvar_pft, nftune
   implicit none

   integer, parameter :: lat_points =  64
   integer, parameter :: lon_points = 128
   integer, parameter :: numcolumn  = 5136
   integer, parameter :: numpatch   = 44400

 ! grid point to extract
   integer, parameter :: jlat = 36 
   integer, parameter :: ilon = 98

   character(len=255) :: fgrid       = "/Users/jidy/Data/run9o/CoLM-T42-gridata-c"
   character(len=255) :: fconst      = "/Users/jidy/Data/run9o/CoLM-T42-const-c"
   character(len=255) :: frestart    = "/Users/jidy/Data/run9o/CoLM-T42-restart-1349"

   character(len=255) :: fgrid_sg    = "/Users/jidy/Cases/run10/CoLM-T42-gridata-c-sg"
   character(len=255) :: fconst_sg   = "/Users/jidy/Cases/run10/CoLM-T42-const-c-sg"
   character(len=255) :: frestart_sg = "/Users/jidy/Cases/run10/CoLM-T42-restart-c-sg"

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

   integer  itypwat(numcolumn)
   integer  ixy_column(numcolumn)
   integer  jxy_column(numcolumn)
   real(r8) wxy_column(numcolumn)
   integer  ixy_patch(numpatch)
   integer  jxy_patch(numpatch)
   real(r8) wxy_patch(numpatch)

   real(r8) fcon_col(nfcon_col,numcolumn)
   real(r8) fvar_col(nfvar_col,numcolumn)
   real(r8) fcon_pft(nfcon_pft,numpatch)
   real(r8) fvar_pft(nfvar_pft,numpatch)

   logical gmask(lon_points,lat_points)

   integer gbegc(lon_points,lat_points), gendc(lon_points,lat_points)

   integer cbegp(numcolumn), cendp(numcolumn)

   integer istep, idate_p(3), idate(3)

   integer i, j, p, c, g, c1, c2, p1, p2, numc, nump

! routine:

   call rdgrid(fgrid,lon_points,lat_points,numcolumn,numpatch, &
               numlon,lats,lonw,area,latixy,longxy,landfrac,landmask, &
               itypwat,ixy_column,jxy_column,wxy_column,ixy_patch,jxy_patch,wxy_patch)

   call rdinit(fconst,frestart,numcolumn,numpatch,istep,idate,idate_p,ftune, &
               fcon_col,fvar_col,fcon_pft,fvar_pft)

   gbegc = -1
   gendc = -1
   gmask = .false.

   p = 1

   do c = 1, numcolumn
      i = ixy_column(c)
      j = jxy_column(c)
      gmask(i,j) = .true.
      if(gbegc(i,j).lt.0) gbegc(i,j) = c
      if(gbegc(i,j).gt.0) gendc(i,j) = c

      if(itypwat(c).eq.0) then
         nump = 17
      else if(itypwat(c).ge.1 .and. itypwat(c).le.4)  then
         nump = 1
      else
         write(6,*) 'impossible itypwat', itypwat(c)
         call abort
      end if

      cbegp(c) = p
      cendp(c) = p + nump - 1

      p = p + nump
   end do

   c1 = gbegc(ilon,jlat)
   c2 = gendc(ilon,jlat)

   p1 = cbegp(c1)
   p2 = cendp(c2)

   nump = 0
   numc = c2-c1+1

   do c = c1, c2
      nump = nump + (cendp(c) - cbegp(c) + 1)
   end do

!*
!* Modify ixy_column, jxy_column, ixy_patch, jxy_patch directly *!
!*
   ixy_column(c1:c2) = 1
   jxy_column(c1:c2) = 1
   ixy_patch(p1:p2)  = 1
   jxy_patch(p1:p2)  = 1
!*
!*

   print *, 'numc', numc
   print *, 'nump', nump
   print *, 'lats', lats(jlat:jlat+1)
   print *, 'lonw', lonw(ilon:ilon+1,jlat)

   call wrgrid(fgrid_sg,1,1,numc,nump, &
               numlon(jlat),lats(jlat:jlat+1),lonw(ilon:ilon+1,jlat),area(ilon,jlat), &
               latixy(ilon,jlat),longxy(ilon,jlat),landfrac(ilon,jlat),landmask(ilon,jlat), &
               itypwat(c1:c2),ixy_column(c1:c2),jxy_column(c1:c2),wxy_column(c1:c2), &
               ixy_patch(p1:p2),jxy_patch(p1:p2),wxy_patch(p1:p2))

   call wrinit(fconst_sg,frestart_sg,numc,nump,istep,idate,idate_p,ftune, &
               fcon_col(:,c1:c2),fvar_col(:,c1:c2),fcon_pft(:,p1:p2),fvar_pft(:,p1:p2))

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

subroutine wrgrid(fgrid,lon_points,lat_points,numcolumn,numpatch, &
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
   integer , intent(in)  :: numlon(lat_points)
   real(r8), intent(in)  :: lats(lat_points+1)
   real(r8), intent(in)  :: lonw(lon_points+1,lat_points)
   real(r8), intent(in)  :: area(lon_points,lat_points)
   real(r8), intent(in)  :: latixy(lon_points,lat_points)
   real(r8), intent(in)  :: longxy(lon_points,lat_points)
   real(r8), intent(in)  :: landfrac(lon_points,lat_points)
   integer , intent(in)  :: landmask(lon_points,lat_points)

   integer , intent(in)  :: itypwat(numcolumn)
   integer , intent(in)  :: ixy_column(numcolumn)
   integer , intent(in)  :: jxy_column(numcolumn)
   real(r8), intent(in)  :: wxy_column(numcolumn)
   integer , intent(in)  :: ixy_patch(numpatch)
   integer , intent(in)  :: jxy_patch(numpatch)
   real(r8), intent(in)  :: wxy_patch(numpatch)

   integer :: lugrid = 100

   OPEN(lugrid,file=fgrid,form='unformatted',status='unknown')
   write(lugrid) numlon(:)
   write(lugrid) lats(:)
   write(lugrid) lonw(:,:)
   write(lugrid) area(:,:)
   write(lugrid) latixy(:,:)
   write(lugrid) longxy(:,:)
   write(lugrid) landfrac(:,:)
   write(lugrid) landmask(:,:)
   write(lugrid) itypwat(:)
   write(lugrid) ixy_column(:) !longitude index for each patch point
   write(lugrid) jxy_column(:) !latitude index for each patch point
   write(lugrid) wxy_column(:) !subgrid weight for each patch point
   write(lugrid) ixy_patch(:)  !longitude index for each patch point
   write(lugrid) jxy_patch(:)  !latitude index for each patch point
   write(lugrid) wxy_patch(:)  !subgrid weight for each patch point
   CLOSE(lugrid)

end subroutine wrgrid
