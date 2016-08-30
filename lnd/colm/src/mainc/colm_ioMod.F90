#include <define.h>

module colm_ioMod

 use precision
 use paramodel, only : nfldv, nfldv_dgvm, nfldv_dgvm_accum
 use nchistMod
 use colm_varctl

 implicit none

 character(len=255) :: site

 character(len=255) :: NLFilename

 integer            :: lugrid               ! logical unit of surface grid data
 integer            :: lusbcini             ! logical unit of initial land surface boundary condition

 integer            :: lurestart            ! logical unit number of restart time-varying file
 integer            :: luconst              ! logical unit number of restart time-invariant file

 integer            :: lulai                ! logical unit number of vegetation data forcing

 character(LEN=255) :: flai                 ! file name of time-varying vegetation data

 character(len=255) :: fgrid                ! file name of surface grids info
 character(len=255) :: fsbcini              ! file name of initial land surface boundary condition 

 character(len=255) :: frestart             ! file name of restart time-varying file
 character(len=255) :: fconst               ! file name of restart time-invariant file

 character(len=255) :: fout                 ! file name of output file (template)

 integer            :: flxctl(nfldv)        ! flag to control fluxes export

#ifdef COUP_CSM
 integer csm_date(3)
#endif

 interface readnml
    module procedure readnml
 end interface

 interface readgridat
    module procedure readgridat
 end interface

 interface readinidat
    module procedure readinidat
 end interface

#ifdef VEGDATA
 interface readvegdat
    module procedure readvegdat
 end interface
#endif

 interface writehistdat
    module procedure writehistdat
 end interface

CONTAINS

   subroutine readnml

      use spmd
      use timemgr
      use forcedata
      use colm_varMod
#ifdef RTM
      use colm_rtmVar
      use colm_rtmMod
#endif
#ifdef COUP_CSM
      use colm_cplMod, only: irad, nsrest
#ifdef CPL6
      use colm_csmMod, only: csm_dtime
#endif
#endif
      use landuse, only: fcarbon_emission
      use shr_sys_mod

      implicit none

! local variables:

      logical lflux  ! existence of flxctl.stdin 
      integer i, j

      namelist /clmexp/ site,             &
                        fgrid,            &
                        flai,             &
                        fmet,             &
                        fout,             &
                        fconst,           &
                        frestart,         &
                        lugrid,           &
                        lulai,            &
                        lumet,            &
                        lurestart,        &
                        luconst,          &
#ifdef RTM
                        frivinp_rtm,      &
                        rtm_nsteps,       &
#endif
                        fcarbon_emission, &
#ifdef COUP_CSM
                        csm_date,         &
                        fsbcini,          &
                        lusbcini,         &
                        irad,             &
                        nsrest,           &
                        csm_doflxave,     &
                        lnd_cflux_year,   &
#ifdef CPL6
                        co2_option,       &
#endif
#ifdef CPL7
                        co2_type,         &
                        co2_ppmv,         &
#endif
#endif
                        lhist_yearly,     &
                        lhist_monthly,    &
                        lhist_daily,      &
                        lhist_3hourly,    &
                        lon_points,       &
                        lat_points,       &
                        numpatch,         &
#ifdef DGVM
                        numcolumn,        &
#endif
                        dtime,            &
                        mstep

      namelist /flxctl_nml/ flxctl

! routine:

      istep = 0    ! COUP_CSM : istep = 0; OFFLINE : istep will be advanced to istep=1.

#ifdef RTM
      rtm_nsteps = -9999
#endif

#ifdef COUP_CSM
      nsrest = 0   ! initial run
      irad = -1
      csm_doflxave = .true.
#endif

      lugrid    = 130
      lusbcini  = 140
      lulai     = 150
      lumet     = 160
      luconst   = 170
      lurestart = 180

      if (p_master) then
         if (len_trim(NLFilename)>0) then
            open(11,file=trim(NLFilename),form='formatted',status='old')
            read(11,nml=clmexp)
            close(11)
         else
            read(5,nml=clmexp)
         end if
      end if

      if (p_master) then 
          write(6,clmexp)
          call shr_sys_flush(6)
          write(6,*) '-------CoLM namelist is ok------'
          call shr_sys_flush(6)
          write(6,*) 'lon =', lon_points
          call shr_sys_flush(6)
          write(6,*) 'lat =', lat_points
          call shr_sys_flush(6)
          write(6,*) 'numpatch =', numpatch
          call shr_sys_flush(6)
          write(6,*) 'numcolumn =', numcolumn
          call shr_sys_flush(6)
          write(6,*) 'dtime =', dtime
          call shr_sys_flush(6)
          write(6,*) 'mstep =', mstep
          call shr_sys_flush(6)
      end if 

      flxctl(:) = (/(i,i=1,nfldv)/)

      if (p_master) then
         inquire(file="./flxctl.stdin", exist=lflux)
         if (lflux) then
            OPEN(7,file="./flxctl.stdin",form='formatted',action='read')
            read(7,nml=flxctl_nml)
            CLOSE(7)
         end if
      end if

      if (trim(fcarbon_emission).ne."null") then
         lcarbon_emission = .TRUE.
      end if

#ifdef SPMD
      call mpi_bcast (dtime     ,1,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (mstep     ,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (numpatch  ,1,mpi_integer,0,p_comm,p_err)
#ifdef DGVM
      call mpi_bcast (numcolumn ,1,mpi_integer,0,p_comm,p_err)
#endif
      call mpi_bcast (lon_points,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (lat_points,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (flxctl,size(flxctl),mpi_integer,0,p_comm,p_err)
#ifdef RTM
      call mpi_bcast (rtm_nsteps,1,mpi_integer,0,p_comm,p_err)
#endif

#ifdef COUP_CSM
      call mpi_bcast (irad,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (nsrest,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (csm_doflxave,1,mpi_logical,0,p_comm,p_err)
      call mpi_bcast (csm_date,size(csm_date),mpi_integer,0,p_comm,p_err)
      call mpi_bcast (lnd_cflux_year,1,mpi_integer,0,p_comm,p_err)
#ifdef CPL6
      call mpi_bcast (co2_option,len(co2_option),mpi_character,0,p_comm,p_err)
#endif
#ifdef CPL7
      call mpi_bcast (co2_type,len(co2_type),mpi_character,0,p_comm,p_err)
      call mpi_bcast (co2_ppmv,1,mpi_real8,0,p_comm,p_err)
#endif
#endif

      call mpi_bcast (lhist_yearly,1,mpi_logical,0,p_comm,p_err)
      call mpi_bcast (lhist_monthly,1,mpi_logical,0,p_comm,p_err)
      call mpi_bcast (lhist_daily,1,mpi_logical,0,p_comm,p_err)
      call mpi_bcast (lhist_3hourly,1,mpi_logical,0,p_comm,p_err)
      call mpi_bcast (lhist_hourly,1,mpi_logical,0,p_comm,p_err)
      call mpi_bcast (lhist_steply,1,mpi_logical,0,p_comm,p_err)

      call mpi_bcast (lcarbon_emission,1,mpi_logical,0,p_comm,p_err)
#endif

      numpatch_glob = numpatch
#ifdef DGVM
      numcolumn_glob = numcolumn
#endif

#ifdef COUP_CSM
#ifdef CPL6
      csm_dtime = dtime
#endif

      if (irad < 0) irad = nint(-irad*3600./dtime)

      if (csm_doflxave .and. irad==1) then
         if (p_master) then
            write(6,*)'error: irad must be greater that one if', &
              ' flux averaging option is enabled'
         endif
         call abort
      endif
#endif

      if (p_master) then
#ifdef RTM
         write(6,*) 'RTM river data = ',trim(frivinp_rtm)
#endif
      end if

   end subroutine readnml

   subroutine readgridat

      use spmd
      use colm_varMod
!*#ifdef COUP_CSM
!*      use colm_csmMod, only : csm_recvgrid
!*#endif

      implicit none

!*#ifdef COUP_CSM
!*      integer , pointer :: cam_numlon(:)      !cam number of longitudes
!*      real(r8), pointer :: cam_longxy(:,:)    !cam lon values
!*      real(r8), pointer :: cam_latixy(:,:)    !cam lat values
!*      real(r8), pointer :: cam_landfrac(:,:)  !cam fractional land
!*      integer , pointer :: cam_landmask(:,:)  !cam land mask
!*#endif

      integer i, j 
 
! routine:

    ! Read in surface grid info.
      if (p_master) then
         OPEN(lugrid,file=fgrid,form='unformatted',status='old')
         read(lugrid) numlon(:)
         read(lugrid) lats(:)
         read(lugrid) lonw(:,:)
         read(lugrid) area(:,:)
         read(lugrid) latixy(:,:)
         read(lugrid) longxy(:,:)
         read(lugrid) landfrac(:,:)
         read(lugrid) landmask(:,:)
         read(lugrid) itypwat_glob(:)
         read(lugrid) ixy_column_glob(:) !longitude index for each patch point
         read(lugrid) jxy_column_glob(:) !latitude index for each patch point
         read(lugrid) wxy_column_glob(:) !subgrid weight for each patch point
         read(lugrid) ixy_patch_glob(:)  !longitude index for each patch point
         read(lugrid) jxy_patch_glob(:)  !latitude index for each patch point
         read(lugrid) wxy_patch_glob(:)  !subgrid weight for each patch point
         CLOSE(lugrid)
      end if

#ifdef SPMD
      call mpi_bcast (numlon         ,size(numlon)         ,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (lats           ,size(lats)           ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (lonw           ,size(lonw)           ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (area           ,size(area)           ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (latixy         ,size(latixy)         ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (longxy         ,size(longxy)         ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (landfrac       ,size(landfrac)       ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (landmask       ,size(landmask)       ,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (itypwat_glob   ,size(itypwat_glob)   ,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (ixy_column_glob,size(ixy_column_glob),mpi_integer,0,p_comm,p_err)
      call mpi_bcast (jxy_column_glob,size(jxy_column_glob),mpi_integer,0,p_comm,p_err)
      call mpi_bcast (wxy_column_glob,size(wxy_column_glob),mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (ixy_patch_glob ,size(ixy_patch_glob) ,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (jxy_patch_glob ,size(jxy_patch_glob) ,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (wxy_patch_glob ,size(wxy_patch_glob) ,mpi_real8  ,0,p_comm,p_err)
#endif

!*#ifdef COUP_CSM
!*      allocate (cam_numlon(lat_points))               !cam number of longitudes
!*      allocate (cam_longxy(lon_points,lat_points))    !cam lon values
!*      allocate (cam_latixy(lon_points,lat_points))    !cam lat values
!*      allocate (cam_landfrac(lon_points,lat_points))  !cam fractional land
!*      allocate (cam_landmask(lon_points,lat_points))  !cam land mask
!*
!*    ! Get grid and land mask back from flux coupler
!*
!*      call csm_recvgrid (cam_longxy, cam_latixy, cam_numlon, cam_landfrac, cam_landmask)
!*
!*    ! Determine land grid, land mask and land fraction
!*
!*    ! For CoLM making surface data by a offline fashion, consistance checking between 
!*    ! land grid info & cpl grid info is required here.
!*
!*      do j = 1, lat_points
!*         if(numlon(j).ne.cam_numlon(j)) then
!*            write(6,*) 'numlon(j).ne.cam_numlon(j), j=', j, numlon(j), cam_numlon(j)
!*            call abort
!*         end if
!*      end do
!*
!*      do j = 1,lat_points
!*         do i = 1,numlon(j)
!*            if(abs(longxy(i,j)-cam_longxy(i,j)).gt.1.0E-8) then
!*               write(6,*) 'longxy(i,j).ne.cam_longxy(i,j), i,j=', i,j, longxy(i,j), cam_longxy(i,j)
!*               call abort
!*            end if
!*            if(abs(latixy(i,j)-cam_latixy(i,j)).gt.1.0E-8) then
!*               write(6,*) 'latixy(i,j).ne.cam_latixy(i,j), i,j=', i,j, latixy(i,j), cam_latixy(i,j)
!*               call abort
!*            end if
!*            if(abs(landfrac(i,j)-cam_landfrac(i,j)).gt.1.0E-8) then
!*               write(6,*) 'landfrac(i,j).ne.cam_landfrac(i,j), i,j=', i,j, landfrac(i,j), cam_landfrac(i,j)
!*             ! call abort
!*            end if
!*            if(abs(landmask(i,j).ne.cam_landmask(i,j))) then
!*               write(6,*) 'landmask(i,j).ne.cam_landmask(i,j), i,j=', i,j, landmask(i,j), cam_landmask(i,j)
!*             ! call abort
!*            end if
!*
!*            longxy(i,j)   = cam_longxy(i,j)
!*            latixy(i,j)   = cam_latixy(i,j)
!*            landmask(i,j) = cam_landmask(i,j)
!*            landfrac(i,j) = cam_landfrac(i,j)
!*         end do
!*      end do
!*
!*      deallocate (cam_numlon  )
!*      deallocate (cam_longxy  )
!*      deallocate (cam_latixy  )
!*      deallocate (cam_landfrac)
!*      deallocate (cam_landmask)
!*#endif

   end subroutine readgridat

#ifdef VEGDATA
   subroutine readvegdat

      use spmd
      use spmd_decomp, only : pcmap, pgmap, gxmap, gymap
      use paramodel, only : numpft
      use colm_varMod, only : lon_points, lat_points, numpatch, mlai, msai, itypwat

      real(r8), allocatable :: laixy(:,:,:,:)
      real(r8), allocatable :: saixy(:,:,:,:)

      integer c, p, g, x, y, c0, p1, p2, p3

      allocate(laixy(lon_points,lat_points,numpft+1,12))   ! including bare soil
      allocate(saixy(lon_points,lat_points,numpft+1,12))   ! including bare soil

      if (p_master) then
         open(lulai,file=trim(flai),form='unformatted',status='old')
         read(lulai) laixy
         read(lulai) saixy
         close(lulai)
      end if

#ifdef SPMD
      call mpi_bcast(laixy,size(laixy),mpi_real8,0,p_comm,p_err)
      call mpi_bcast(saixy,size(laixy),mpi_real8,0,p_comm,p_err)
#endif

      mlai = 0.
      msai = 0.

      c0 = 0

      do p = 1, numpatch
         c = pcmap(p)
         g = pgmap(p)
         x = gxmap(g)
         y = gymap(g)

         if(c.ne.c0 .and. itypwat(c).eq.0) then
            p1 = p
            p2 = p+numpft
            do p3 = p1, p2
               mlai(:,p3) = laixy(x,y,p3-p1+1,:)
               msai(:,p3) = saixy(x,y,p3-p1+1,:)
            end do
         end if

         c0 = c
      end do

      deallocate(laixy)
      deallocate(saixy)

   end subroutine readvegdat
#endif

   subroutine readinidat

      use spmd
      use spmd_decomp, only : ppmap, ccmap, cgmap, gxmap, gymap
      use timemgr
      use colm_varMod
#ifdef RTM
      use RtmMod, only: restart_rtm
#endif
#ifdef COUP_CSM
      use colm_cplMod
#ifdef CPL6
      use colm_csmMod, only: csm_restart
#endif
#endif

      implicit none

! local variables:

      real(r8), pointer :: rbuf_glob(:)     ! buffer for reading fvar/fcon
      real(r8), pointer :: lndxy_avsdr(:,:) ! avsdr
      real(r8), pointer :: lndxy_avsdf(:,:) ! avsdf
      real(r8), pointer :: lndxy_anidr(:,:) ! anidr
      real(r8), pointer :: lndxy_anidf(:,:) ! anidf
      real(r8), pointer :: lndxy_scv(:,:)   ! snow water [mm]
      real(r8), pointer :: lndxy_trad(:,:)  ! radiative temp. [K]

#ifdef MYBUG
      real(r4), pointer :: frac(:,:,:)
#endif

      integer i, j, k, g

! routine:

      itypwat(:)    = itypwat_glob(ccmap)
      wxy_column(:) = wxy_column_glob(ccmap)
      wxy_patch(:)  = wxy_patch_glob(ppmap)

#ifdef MYBUG
      allocate(frac(lon_points,lat_points,5))
      frac(:,:,:) = -9999.

      do i = 1, numcolumn
         if(itypwat(i)<5) then
            k = itypwat(i)+1
            frac(gxmap(cgmap(i)),gymap(cgmap(i)),k) = wxy_column(i)
         end if
      end do

      do i = 1, 5
         write(100+p_iam), frac(:,:,i)
      end do

      deallocate(frac)
#endif

      write(6,*), p_iam, 'before const read'

      if (p_master) then
       ! Open for model time invariant constant data
         OPEN(unit=luconst,file=fconst,form='unformatted',&
                           status='old',action='read')

         read(luconst) ftune       !clm tunable constants
      end if

#ifdef SPMD
      call mpi_bcast (ftune,size(ftune),mpi_real8,0,p_comm,p_err)
#endif

      allocate (rbuf_glob(numcolumn_glob))

      do k = 1, nfcon_col
         if (p_master) read(luconst) rbuf_glob
#ifdef SPMD
         call mpi_bcast (rbuf_glob,size(rbuf_glob),mpi_real8,0,p_comm,p_err)
#endif
         fcon_col(k,:) = rbuf_glob(ccmap)
      end do

      deallocate (rbuf_glob)

      do i = 1, numcolumn
         if(itypwat(i) .ne. nint(fcon_col(3,i))) then
            write(6,*) 'fatal error in itypwat checking in readinidat', itypwat(i), nint(fcon_col(4,i))
            call abort
         end if
      end do

      allocate (rbuf_glob(numpatch_glob))

      do k = 1, nfcon_pft
         if (p_master) read(luconst) rbuf_glob
#ifdef SPMD
         call mpi_bcast (rbuf_glob,size(rbuf_glob),mpi_real8,0,p_comm,p_err)
#endif
         fcon_pft(k,:) = rbuf_glob(ppmap)
      end do

      deallocate (rbuf_glob)

      if (p_master) CLOSE(luconst)

      write(6,*), p_iam, 'after const read'

      if (p_master) then
       ! Open for model time varying data (model state variables)
         OPEN(unit=lurestart,file=frestart,form='unformatted',&
                             status='old',action='read')

         read(lurestart) istep, idate, idate_p
      end if

#ifdef SPMD
      call mpi_bcast (istep,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (idate,size(idate),mpi_integer,0,p_comm,p_err)
      call mpi_bcast (idate_p,size(idate_p),mpi_integer,0,p_comm,p_err)
#endif

#ifdef COUP_CSM
      if(nsrest.eq.0) then
         istep = 0  ! clear random istep in restart file
         idate_p = csm_date
         call set_curr_date(csm_date)
      end if
#else
      call set_curr_date(idate)
#endif

      allocate (rbuf_glob(numcolumn_glob))

      do k = 1, nfvar_col
         if (p_master) read(lurestart) rbuf_glob
#ifdef SPMD
         call mpi_bcast (rbuf_glob,size(rbuf_glob),mpi_real8,0,p_comm,p_err)
#endif
         fvar_col(k,:) = rbuf_glob(ccmap)
      end do

      deallocate (rbuf_glob)

      allocate (rbuf_glob(numpatch_glob))

      do k = 1, nfvar_pft
#ifdef RESTART_FMT1
         if ((k+5).le.nfvar_pft) then
            if (p_master) read(lurestart) rbuf_glob
         else
            if (p_master) rbuf_glob = 0._r8
         end if
#else
         if (p_master) read(lurestart) rbuf_glob
#endif
#ifdef SPMD
         call mpi_bcast (rbuf_glob,size(rbuf_glob),mpi_real8,0,p_comm,p_err)
#endif
         fvar_pft(k,:) = rbuf_glob(ppmap)
      end do

      deallocate (rbuf_glob)

#ifdef RTM
      if (nsrest.gt.0) call restart_rtm(lurestart,'read')
#endif

#ifdef COUP_CSM
#ifdef CPL6
      if (nsrest.gt.0) call csm_restart(lurestart,'read') 
#endif
#endif

      if (p_master) CLOSE(lurestart)

      write(6,*), p_iam, 'after restart read'

#ifdef COUP_CSM
      if (nsrest.eq.0) then   ! initial run

         allocate (lndxy_avsdr(lon_points,lat_points))
         allocate (lndxy_avsdf(lon_points,lat_points))
         allocate (lndxy_anidr(lon_points,lat_points))
         allocate (lndxy_anidf(lon_points,lat_points))
         allocate (lndxy_trad (lon_points,lat_points))
         allocate (lndxy_scv  (lon_points,lat_points))

         if (p_master) then
            OPEN(unit=lusbcini,file=fsbcini,form='unformatted',status='old',action='read')
            read(lusbcini) lndxy_avsdr(:,:) ! avsdr
            read(lusbcini) lndxy_avsdf(:,:) ! avsdf
            read(lusbcini) lndxy_anidr(:,:) ! anidr
            read(lusbcini) lndxy_anidf(:,:) ! anidf
            read(lusbcini) lndxy_trad
            read(lusbcini) lndxy_scv
            CLOSE(lusbcini)
         end if

         write(6,*), p_iam, 'after sbcini read'

#ifdef SPMD
         call mpi_bcast (lndxy_avsdr,size(lndxy_avsdr),mpi_real8,0,p_comm,p_err)
         call mpi_bcast (lndxy_avsdf,size(lndxy_avsdf),mpi_real8,0,p_comm,p_err)
         call mpi_bcast (lndxy_anidr,size(lndxy_anidr),mpi_real8,0,p_comm,p_err)
         call mpi_bcast (lndxy_anidf,size(lndxy_anidf),mpi_real8,0,p_comm,p_err)
         call mpi_bcast (lndxy_trad ,size(lndxy_trad) ,mpi_real8,0,p_comm,p_err)
         call mpi_bcast (lndxy_scv  ,size(lndxy_scv)  ,mpi_real8,0,p_comm,p_err)
#endif

         write(6,*), p_iam, 'after sbcini bcast'

         do g = 1,numgrid
            i = gxmap(g)
            j = gymap(g)

            lnd_avsdr(g) = lndxy_avsdr(i,j)
            lnd_avsdf(g) = lndxy_avsdf(i,j)
            lnd_anidr(g) = lndxy_anidr(i,j)
            lnd_anidf(g) = lndxy_anidf(i,j)
            lnd_trad(g)  = lndxy_trad(i,j)
            lnd_scv(g)   = lndxy_scv(i,j)/1000._r8
         end do

#ifdef MYBUG
         write(6,*), 'min & max avsdr', minval(lnd_avsdr), maxval(lnd_avsdr)
         write(6,*), 'min & max avsdf', minval(lnd_avsdf), maxval(lnd_avsdf)
         write(6,*), 'min & max anidr', minval(lnd_anidr), maxval(lnd_anidr)
         write(6,*), 'min & max anidf', minval(lnd_anidf), maxval(lnd_anidf)
         write(6,*), 'min & max trad', minval(lnd_trad), maxval(lnd_trad)
         write(6,*), 'min & max scv', minval(lnd_scv), maxval(lnd_scv)
#endif

         deallocate (lndxy_avsdr)
         deallocate (lndxy_avsdf)
         deallocate (lndxy_anidr)
         deallocate (lndxy_anidf)
         deallocate (lndxy_trad)
         deallocate (lndxy_scv)
      end if
#endif

   end subroutine readinidat

   subroutine writehistdat

      use spmd, only: p_master
      use timemgr, only: idate, idate_p
      use colm_varMod, only: fldv
#ifdef DGVM
      use colm_varMod, only: fldv_dgvm
#endif
#ifdef RTM
      use RtmMod, only: restart_rtm
      use colm_rtmMod
#endif
#ifdef COUP_CSM
#ifdef CPL6
      use colm_csmMod, only: csm_restart
#endif
#endif
      implicit none

! arguments:

! local variables:

      integer n, k
      character(len=255) :: cdate

! routine:

    ! Accumulating

      do n = 1, nMaxHist
         if(histArray(n)%is_valid) then
            histArray(n)%nac = histArray(n)%nac + 1
            do k = 1, nfldv
               if(histArray(n)%fldvmask(k)) then
                  histArray(n)%fldv(k,:) = histArray(n)%fldv(k,:) + fldv(k,:)
               end if
            end do
#ifdef RTM
            if(histArray(n)%with_rtm) then
               call colm_rtm_accum(histArray(n)%lnd_nac,histArray(n)%ocn_nac,&
                                   histArray(n)%lnd_ave,histArray(n)%ocn_ave)
            end if
#endif
#ifdef DGVM
            if(histArray(n)%with_dgvm) then
               histArray(n)%nac_dgvm = histArray(n)%nac_dgvm + 1
               do k = 1, nfldv_dgvm_accum
                  if(histArray(n)%fldvmask_dgvm(k)) then
                     histArray(n)%fldv_dgvm(k,:) = histArray(n)%fldv_dgvm(k,:) + fldv_dgvm(k,:)
                  end if
               end do
            end if
#endif
         end if
      end do

      frestart = "null"

      call lpwrite(idate_p,idate,fout,frestart,cdate)

      do n = 1, nMaxHist

         if (histArray(n)%is_valid .and. histArray(n)%is_ready) then

            if (p_master) then
               call nchist_create(idate,cdate,histArray(n))
            end if

            call write_colm_history(histArray(n))

#ifdef RTM
            if(histArray(n)%with_rtm) call write_rtm_history(histArray(n)%lnd_nac, &
                                                             histArray(n)%ocn_nac, &
                                                             histArray(n)%lnd_ave, &
                                                             histArray(n)%ocn_ave, &
                                                             histArray(n)%fldvmask_rtm)
#endif

#ifdef DGVM
            if(histArray(n)%with_dgvm) call write_dgvm_history(histArray(n))
#endif

            if (p_master) then
                call nchist_close
            end if

         end if

      end do

      if (trim(frestart).ne."null") then
         if (p_master) OPEN(unit=lurestart,file=frestart,access='sequential', &
                            form='unformatted',status='unknown')

         call write_colm_restart

#ifdef RTM
         call restart_rtm(lurestart,'write')
#endif

#ifdef COUP_CSM
#ifdef CPL6
         if (p_master) call csm_restart(lurestart,'write')
#endif
#endif

         if (p_master) CLOSE(lurestart)

         call write_colm_sbc

         if (p_master) call write_restart_nml(frestart)
      end if

    ! if((idate(1)*1000+idate(2)).gt.7360) lcycdump = .true.

      if (lcycdump) then
         call cycdump
      end if

   end subroutine writehistdat

   subroutine write_colm_history(histinfo)

      use spmd
      use spmd_decomp
      use colm_varMod

      implicit none

      type(history_info), intent(inout) :: histinfo

      real(r8), pointer :: fldv_loc(:)
      real(r8), pointer :: fldv_glob(:,:)
      real(r4), pointer :: fldxy(:,:)
      real(r8), pointer :: tsn(:)
      real(r8), pointer :: nsnow(:)

      integer i, j, g, L

    ! calculate average

      do L = 1, nfldv
         if (histinfo%fldvmask(L)) then
            histinfo%fldv(L,:) = histinfo%fldv(L,:)/float(histinfo%nac)
         end if
      end do

#ifdef CMIP
    ! Special treatment for TSN
      if(histinfo%fldvmask(idx_tsn)) then
         if(histinfo%fldvmask(idx_nsnow)) then
            tsn => histinfo%fldv(idx_tsn,:)
            nsnow => histinfo%fldv(idx_nsnow,:)

            do g = 1, numgrid
               if(nsnow(g).gt.0._r8) then
                  tsn(g) = tsn(g)/nsnow(g)
               else
                  tsn(g) = -9999.
               end if
            end do
         else
            write(6,*) 'fatal error on calculating tsn'
            call abort
         end if
      end if
#endif

      allocate (fldv_loc(numgrid))
      allocate (fldv_glob(numgrid_glob,nfldv))
      allocate (fldxy(lon_points,lat_points))

#ifdef SPMD
      p_counts(:) = numgrid_proc(:)
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numgrid_proc(0:i-1))
      end do
#endif

      do L = 1, nfldv
         if (histinfo%fldvmask(L)) then
#ifdef SPMD
            fldv_loc(:) = histinfo%fldv(L,:)
            call mpi_gatherv (fldv_loc,size(fldv_loc),mpi_real8,&
                              fldv_glob(:,L),p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
            fldv_glob(:,L) = histinfo%fldv(L,:)
#endif
         end if
      end do

      if (p_master) then
         do L = 1, nfldv
            if (histinfo%fldvmask(L)) then
               fldxy(:,:) = -9999.

               do g = 1, numgrid_glob
                  i = gxmap_glob(g)
                  j = gymap_glob(g)
                  fldxy(i,j) = fldv_glob(g,L)
               end do

               call nchist_addfield(L,lon_points,lat_points,fldxy,'lsm')
            end if
         end do
      end if

      deallocate (fldxy)
      deallocate (fldv_loc)
      deallocate (fldv_glob)

    ! reinitialize

      histinfo%nac = 0
      do L = 1, nfldv
         if(histinfo%fldvmask(L)) histinfo%fldv(L,:) = 0.
      end do

   end subroutine write_colm_history

#ifdef DGVM
   subroutine write_dgvm_history(histinfo)

      use spmd
      use spmd_decomp
      use paramodel, only: nfldv_dgvm
      use colm_varMod

      implicit none

      type(history_info), intent(inout) :: histinfo

      real(r8), pointer :: fldv_loc(:)
      real(r8), pointer :: fldv_glob(:,:)
      real(r4), pointer :: fldxy(:,:)

      integer i, j, c, g, L

    ! calculate average

      do L = 1, nfldv_dgvm
         if(histinfo%fldvmask_dgvm(L)) then
            if(L.le.nfldv_dgvm_accum) then
               histinfo%fldv_dgvm(L,:) = histinfo%fldv_dgvm(L,:)/float(histinfo%nac_dgvm)
            else
               histinfo%fldv_dgvm(L,:) = fldv_dgvm(L,:)
            end if
         end if
      end do

      allocate (fldv_loc(numgrid))
      allocate (fldv_glob(numgrid_glob,nfldv_dgvm))
      allocate (fldxy(lon_points,lat_points))

#ifdef SPMD
      p_counts(:) = numgrid_proc(:)
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numgrid_proc(0:i-1))
      end do
#endif

      do L = 1, nfldv_dgvm
         if(histinfo%fldvmask_dgvm(L)) then
#ifdef SPMD
            fldv_loc(:) = histinfo%fldv_dgvm(L,:)
            call mpi_gatherv(fldv_loc,size(fldv_loc),mpi_real8,&
                             fldv_glob(:,L),p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
            fldv_glob(:,L) = histinfo%fldv_dgvm(L,:)
#endif
         end if
      end do

      if(p_master) then
         do L = 1, nfldv_dgvm
            if(histinfo%fldvmask_dgvm(L)) then
               fldxy(:,:) = -9999.
   
               do g = 1, numgrid_glob
                  i = gxmap_glob(g)
                  j = gymap_glob(g)
                  fldxy(i,j) = fldv_glob(g,L)
               end do

               call nchist_addfield(L,lon_points,lat_points,fldxy,'dgvm')
            end if
         end do
      end if

      deallocate (fldxy)
      deallocate (fldv_loc)
      deallocate (fldv_glob)

    ! reinitialize

      histinfo%nac_dgvm = 0
      do L = 1, nfldv_dgvm
         if(histinfo%fldvmask_dgvm(L)) histinfo%fldv_dgvm(L,:) = 0.
      end do

   end subroutine write_dgvm_history
#endif

 ! Write out the model variables for restart run 

   subroutine write_colm_restart

      use spmd
      use spmd_decomp
      use timemgr, only: istep, idate, idate_p
      use colm_varMod

      implicit none

      real(r8), pointer :: fvar_loc(:)
      real(r8), pointer :: fvar_glob(:,:)
      real(r8), pointer :: fvar_buf(:)

      integer i, L

    ! column level

      allocate (fvar_buf(numcolumn_glob))
      allocate (fvar_loc(numcolumn))
      allocate (fvar_glob(numcolumn_glob,nfvar_col))

#ifdef SPMD
      p_counts(:) = numcolumn_proc(:)
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numcolumn_proc(0:i-1))
      end do
#endif

      do L = 1, nfvar_col
#ifdef SPMD
         fvar_loc(:) = fvar_col(L,:)
         call mpi_gatherv (fvar_loc,size(fvar_loc),mpi_real8,&
                           fvar_glob(:,L),p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
         fvar_glob(:,L) = fvar_col(L,:)
#endif
      end do

      if (p_master) then
         write(lurestart) istep, idate, idate_p

         do L = 1, nfvar_col
            fvar_buf(ccmap_glob) = fvar_glob(:,L)
            write(lurestart) fvar_buf
         end do
      end if

      deallocate (fvar_buf)
      deallocate (fvar_loc)
      deallocate (fvar_glob)

    ! patch level

      allocate (fvar_buf(numpatch_glob))
      allocate (fvar_loc(numpatch))
      allocate (fvar_glob(numpatch_glob,nfvar_pft))

#ifdef SPMD
      p_counts(:) = numpatch_proc(:)
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numpatch_proc(0:i-1))
      end do
#endif

      do L = 1, nfvar_pft
#ifdef SPMD
         fvar_loc(:) = fvar_pft(L,:)
         call mpi_gatherv (fvar_loc,size(fvar_loc),mpi_real8,&
                           fvar_glob(:,L),p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
         fvar_glob(:,L) = fvar_pft(L,:)
#endif
      end do

      if (p_master) then
         do L = 1, nfvar_pft
            fvar_buf(ppmap_glob) = fvar_glob(:,L)
            write(lurestart) fvar_buf
         end do
      end if

      deallocate (fvar_buf)
      deallocate (fvar_loc)
      deallocate (fvar_glob)

   end subroutine write_colm_restart

   subroutine write_colm_sbc

      use spmd
      use spmd_decomp
      use colm_varMod

      implicit none

      real(r8), pointer :: fldv_loc(:)
      real(r8), pointer :: fldv_glob(:,:)
      real(r8), pointer :: fldxy(:,:)

      logical :: fldvmask(nfldv)

      integer i, j, g, L

      fldvmask(:) = .false.

      fldvmask(idx_avsdr) = .true.
      fldvmask(idx_avsdf) = .true.
      fldvmask(idx_anidr) = .true.
      fldvmask(idx_anidf) = .true.
      fldvmask(idx_trad ) = .true.
      fldvmask(idx_scv  ) = .true.

      allocate (fldv_loc(numgrid))
      allocate (fldv_glob(numgrid_glob,nfldv))
      allocate (fldxy(lon_points,lat_points))

#ifdef SPMD
      p_counts(:) = numgrid_proc(:)
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numgrid_proc(0:i-1))
      end do
#endif

      do L = 1, nfldv
         if (fldvmask(L)) then
#ifdef SPMD
            fldv_loc(:) = fldv(L,:)
            call mpi_gatherv (fldv_loc,size(fldv_loc),mpi_real8,&
                              fldv_glob(:,L),p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
            fldv_glob(:,L) = fldv(L,:)
#endif
         end if
      end do

      if (p_master) then

         open(lusbcini,file=trim(frestart)//"-sbc",form="unformatted",status="unknown")

         do L = 1, nfldv
            if (fldvmask(L)) then
               fldxy(:,:) = -9999.

               do g = 1, numgrid_glob
                  i = gxmap_glob(g)
                  j = gymap_glob(g)
                  fldxy(i,j) = fldv_glob(g,L)
               end do

               write(lusbcini) fldxy
            end if
         end do

         close(lusbcini)

      end if

      deallocate (fldxy)
      deallocate (fldv_loc)
      deallocate (fldv_glob)

   end subroutine write_colm_sbc

   subroutine write_restart_nml(frestart)

      use timemgr
      use forcedata
      use colm_varMod
#ifdef RTM
      use colm_rtmVar
      use colm_rtmMod
#endif
#ifdef COUP_CSM
      use colm_cplMod, only: irad, nsrest
#endif

      implicit none

      character(len=255), intent(in) :: frestart

      integer :: lurestnml = 111

      open(unit=lurestnml,file="lnd.stdin.rst",form="formatted",status="unknown")

      write(lurestnml,10) "&clmexp"
      write(lurestnml,20) "site = ", "'"//trim(site)//"'"
      write(lurestnml,20) "fgrid = ", "'"//trim(fgrid)//"'"
      write(lurestnml,30) "lugrid = ", lugrid

#if(defined VEGDATA)
      write(lurestnml,20) "flai = ", "'"//trim(flai)//"'"
      write(lurestnml,30) "lulai = ", lulai
#endif
#ifdef OFFLINE
      write(lurestnml,20) "fmet = ", "'"//trim(fmet)//"'"
      write(lurestnml,30) "lumet = ", lumet
#endif

      write(lurestnml,20) "fconst = ", "'"//trim(fconst)//"'"
      write(lurestnml,30) "luconst = ", luconst

      write(lurestnml,20) "frestart = ", "'"//trim(frestart)//"'"
      write(lurestnml,30) "lurestart = ", lurestart

      write(lurestnml,20) "fout = ", "'"//trim(fout)//"'"

#ifdef RTM
      write(lurestnml,20) "frivinp_rtm = ", "'"//trim(frivinp_rtm)//"'"
#endif
#ifdef COUP_CSM
      write(lurestnml,30) "fsbcini = "//"''"
      write(lurestnml,30) "lusbcini = ", lusbcini
      write(lurestnml,30) "irad = ", irad
      write(lurestnml,30) "nsrest = ", 1
      write(lurestnml,60) "csm_doflxave = ", csm_doflxave
      write(lurestnml,40) "csm_date = ", idate(1), idate(2), idate(3)
      write(lurestnml,30) "lnd_cflux_year = ", lnd_cflux_year
#ifdef CPL6
      write(lurestnml,20) "co2_option = ", "'"//trim(co2_option)//"'"
#endif
#ifdef CPL7
      write(lurestnml,20) "co2_type = ", "'"//trim(co2_type)//"'"
      write(lurestnml,50) "co2_ppmv = ", co2_ppmv
#endif
#endif
      write(lurestnml,60) "lhist_yearly = ", lhist_yearly
      write(lurestnml,60) "lhist_monthly = ", lhist_monthly
      write(lurestnml,60) "lhist_daily = ", lhist_daily
      write(lurestnml,60) "lhist_3hourly = ", lhist_3hourly
#ifdef DGVM
      write(lurestnml,30) "numcolumn = ", numcolumn_glob
#endif
      write(lurestnml,30) "numpatch = ", numpatch_glob
      write(lurestnml,30) "lon_points = ", lon_points
      write(lurestnml,30) "lat_points = ", lat_points
      write(lurestnml,50) "dtime = ", dtime
      write(lurestnml,30) "mstep = ", mstep
      write(lurestnml,10), "/"

10 format(A)
20 format(A20,A)
30 format(A20,I10)
40 format(A20,I4,',',I3,',',I5)
50 format(A20,F10.2)
60 format(A20,L)

      close(lurestnml)

   end subroutine write_restart_nml

   subroutine cycdump

      use spmd
      use spmd_decomp
      use timemgr, only: istep, idate, idate_p
      use colm_varMod

      implicit none

      real(r8), pointer :: forc_loc(:)
      real(r8), pointer :: forc_glob(:,:)
      real(r8), pointer :: forc_buf(:)

      real(r8), pointer :: fvar_loc(:)
      real(r8), pointer :: fvar_glob(:,:)
      real(r8), pointer :: fvar_buf(:)

      integer i, L

      character(len=255) :: fcycdumpI

    ! forcing variables

#ifdef SPMD
      p_counts(:) = numcolumn_proc(:)
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numcolumn_proc(0:i-1))
      end do
#endif

      allocate (forc_buf(numcolumn_glob))
      allocate (forc_loc(numcolumn))
      allocate (forc_glob(numcolumn_glob,nforc))

      do L = 1, nforc
#ifdef SPMD
         forc_loc(:) = forc(L,:)
         call mpi_gatherv (forc_loc,size(forc_loc),mpi_real8,&
                           forc_glob(:,L),p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
         forc_glob(:,L) = forc(L,:)
#endif
      end do

      if (p_master) then
         write(fcycdumpI,'(I3.3)') mod(istep,ncycdump)
         fcycdumpI = trim(fcycdump)//"-"//trim(fcycdumpI)

         OPEN(unit=lucycdump,file=fcycdumpI,access='sequential', &
                             form='unformatted',status='unknown')
         write(lucycdump) istep, idate, idate_p

         do L = 1, nforc
            forc_buf(ccmap_glob) = forc_glob(:,L)
            write(lucycdump) forc_buf
         end do
      end if

      deallocate (forc_buf)
      deallocate (forc_loc)
      deallocate (forc_glob)

    ! column level variables

      allocate (fvar_buf(numcolumn_glob))
      allocate (fvar_loc(numcolumn))
      allocate (fvar_glob(numcolumn_glob,nfvar_col))

#ifdef SPMD
      p_counts(:) = numcolumn_proc(:)
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numcolumn_proc(0:i-1))
      end do
#endif

      do L = 1, nfvar_col
#ifdef SPMD
         fvar_loc(:) = fvar_col(L,:)
         call mpi_gatherv (fvar_loc,size(fvar_loc),mpi_real8,&
                           fvar_glob(:,L),p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
         fvar_glob(:,L) = fvar_col(L,:)
#endif
      end do

      if (p_master) then
         do L = 1, nfvar_col
            fvar_buf(ccmap_glob) = fvar_glob(:,L)
            write(lucycdump) fvar_buf
         end do
      end if

      deallocate (fvar_buf)
      deallocate (fvar_loc)
      deallocate (fvar_glob)

    ! patch level variables

      allocate (fvar_buf(numpatch_glob))
      allocate (fvar_loc(numpatch))
      allocate (fvar_glob(numpatch_glob,nfvar_pft))

#ifdef SPMD
      p_counts(:) = numpatch_proc(:)
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numpatch_proc(0:i-1))
      end do
#endif

      do L = 1, nfvar_pft
#ifdef SPMD
         fvar_loc(:) = fvar_pft(L,:)
         call mpi_gatherv (fvar_loc,size(fvar_loc),mpi_real8,&
                           fvar_glob(:,L),p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
         fvar_glob(:,L) = fvar_pft(L,:)
#endif
      end do

      if (p_master) then
         do L = 1, nfvar_pft
            fvar_buf(ppmap_glob) = fvar_glob(:,L)
            write(lucycdump) fvar_buf
         end do

         CLOSE(lucycdump)
      end if

      deallocate (fvar_buf)
      deallocate (fvar_loc)
      deallocate (fvar_glob)

   end subroutine cycdump

END MODULE colm_ioMod
