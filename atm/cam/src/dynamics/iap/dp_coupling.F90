
!-------------------------------------------------------------------------------
! dynamics - physics coupling module
! Modified: Jiang Jinrong and Zhang He, 2012-10-24
!           Zhang He, 2013-02-01, new mp_sendirr
!           Zhang He, 2013-03-05 
!           He Juanxiong, 2013-09-25, for nudging and two-way coupling
!-------------------------------------------------------------------------------
module dp_coupling

   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,        only: pcols, pver
   use rgrid,         only: nlon
   use pmgrid,        only: plev, plat, plon, beglat, endlat, &
                            twod_decomp, beglev, endlev,  &
                            beglonxy,endlonxy,beglatxy,endlatxy
   use phys_buffer,   only: pbuf_fld, pbuf_size_max
   use phys_grid
   use physics_types, only: physics_state, physics_tend
   use constituents,  only: pcnst
   use physconst,     only: cpair, gravit, rair, zvir
   use geopotential,  only: geopotential_t
   use wv_saturation, only: aqsat ! juanxiong he
   use check_energy,  only: check_energy_timestep_init
   use dynamics_vars, only: T_FVDYCORE_GRID, T_FVDYCORE_STATE   !zhh 2013-02-01
   use dyn_internal_state,   only : get_dyn_state      !zhh 2013-02-01
#if (defined SPMD)
   use spmd_dyn,      only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs  !zhh
   use mpishorthand, only : mpicom
   use parutilitiesmodule, only : sumop, parcollective
   use mod_comm, only: mp_sendirr, mp_recvirr
#endif
   use abortutils,    only: endrun
   use perf_mod
!======================== zhh debug =======================
   use spmd_utils,   only: masterproc,iam
!======================== zhh debug =======================

   implicit none

!===============================================================================
CONTAINS
!===============================================================================

!===============================================================================
  subroutine d_p_coupling(ps, t3, u3, v3, q3, &
                          omga, phis, phys_state, phys_tend, pbuf, pdeld)
!------------------------------------------------------------------------------
! Coupler for converting dynamics output variables into physics input variables
! also writes dynamics variables (on physics grid) to history file
!------------------------------------------------------------------------------
    use physconst,     only: cappa
    use constituents,  only: cnst_get_type_byind
    use physics_types, only: set_state_pdry

!------------------------------Arguments--------------------------------
    real(r8), intent(in) :: ps  (beglonxy:endlonxy, beglatxy:endlatxy)            ! surface pressure
    real(r8), intent(in) :: t3  (beglonxy:endlonxy, plev, beglatxy:endlatxy)  ! temperature
    real(r8), intent(in) :: u3  (beglonxy:endlonxy, plev, beglatxy:endlatxy)  ! u-wind component
    real(r8), intent(in) :: v3  (beglonxy:endlonxy, plev, beglatxy:endlatxy)  ! v-wind component
    real(r8), intent(in) :: q3  (beglonxy:endlonxy, plev, pcnst, beglatxy:endlatxy) ! constituents
    real(r8), intent(in) :: omga(beglonxy:endlonxy, plev, beglatxy:endlatxy)      ! vertical velocity
    real(r8), intent(in) :: phis(beglonxy:endlonxy, beglatxy:endlatxy)            ! Surface geopotential
    real(r8), intent(in) :: pdeld (beglonxy:endlonxy,plev,beglatxy:endlatxy)  

    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(out),   dimension(begchunk:endchunk) :: phys_tend
    type(pbuf_fld),      intent(inout), dimension(pbuf_size_max)     :: pbuf
!
!---------------------------Local workspace-----------------------------
#if (! defined SPMD)
    real(r8) :: buf1(1), buf2(1)     ! transpose buffers
    integer  :: buf1win, buf2win     ! MPI-2 window ids
    integer  :: spmdbuf_siz = 0
    integer  :: block_buf_nrecs = 0
    integer  :: chunk_buf_nrecs = 0
    integer  :: mpicom = 0
    logical  :: local_dp_map=.true. 
    integer :: test=10
#else
    integer :: test=-10
#endif

    integer :: i,k,j,m,lchnk,ib      ! indices
    integer :: ncol                  ! number of columns in current chunk
    integer :: lats(pcols)           ! array of latitude indices
    integer :: lons(pcols)           ! array of longitude indices
    integer :: tsize                 ! amount of data per grid point passed to physics
    integer, allocatable, dimension(:,:) :: bpter !jjr
    real(r8), allocatable, dimension(:) :: buf1,buf2 
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data
    logical :: wetq(pcnst)           ! 'moist-type' constituent flag
    real(r8) :: rlat(pcols)          ! array of latitudes (radians)
    real(r8) :: rlon(pcols)          ! array of longitudes (radians)
    integer:: blksiz
!
!-----------------------------------------------------------------------

! Determine which constituents are wet and which are dry
    do m=2,pcnst
       if (cnst_get_type_byind(m).eq.'wet') then  
          wetq(m) = .true.
       else
          wetq(m) = .false.
       endif
    enddo

!-----------------------------------------------------------------------
! copy data from dynamics data structure to physics data structure
!-----------------------------------------------------------------------
    if (local_dp_map) then

!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, I, K, M, LONS, LATS)
       do lchnk = begchunk,endchunk
          ncol = phys_state(lchnk)%ncol
          call get_lon_all_p(lchnk, ncol, lons)
          call get_lat_all_p(lchnk, ncol, lats)
          do i=1,ncol
             phys_state(lchnk)%ps   (i)     = ps  (lons(i),lats(i))
             phys_state(lchnk)%phis (i)     = phis(lons(i),lats(i))
          end do
  
          do k=1,plev
             do i=1,ncol
                phys_state(lchnk)%t    (i,k)   = t3  (lons(i),k,lats(i))
                phys_state(lchnk)%u    (i,k)   = u3  (lons(i),k,lats(i))
                phys_state(lchnk)%v    (i,k)   = v3  (lons(i),k,lats(i))
                phys_state(lchnk)%omega(i,k)   = omga(lons(i),k,lats(i))
                phys_state(lchnk)%q(i,k,1)     = q3  (lons(i),k,1,lats(i))
             end do
          end do

          do k=1,plev
             do i=1,ncol
                phys_state(lchnk)%pdeldry(i,k) = pdeld(lons(i),k,lats(i))
             end do
          end do

          ! convert moist-type constituents from dry to moist mixing ratio

          do m=2,pcnst
             if (wetq(m)) then  
                do k=1,plev
                   do i=1,ncol
                      phys_state(lchnk)%q(i,k,m) = q3(lons(i),k,m,lats(i))*(1._r8 - q3(lons(i),k,1,lats(i)))
                   end do
                end do
             else
                do k=1,plev
                   do i=1,ncol
                      phys_state(lchnk)%q(i,k,m) = q3(lons(i),k,m,lats(i))
                   end do
                end do
             endif
          end do

       end do

    else

       tsize = 5 + pcnst  

!jjr       if (tsize*max(block_buf_nrecs,chunk_buf_nrecs) > spmdbuf_siz) then
!jjr          call endrun ('p_d_coupling: communication buffers (spmdbuf_siz) too small')
!jjr       endif
!
       blksiz = (endlatxy-beglatxy+1)*(endlonxy-beglonxy+1)
       allocate(bpter(blksiz,0:plev))
       allocate(buf1(tsize*block_buf_nrecs))
       allocate(buf2(tsize*chunk_buf_nrecs))

#ifdef OUTER_OMP
!$OMP PARALLEL DO PRIVATE (J, BPTER, I, K, M)
#endif

       call block_to_chunk_send_pters(iam+1,blksiz,plev+1,tsize,bpter)

!DIR$ CONCURRENT
         ib = 0
         do j=beglatxy,endlatxy
          do i=beglonxy,endlonxy
            ib = ib + 1
             buf1(bpter(ib,0))   = ps  (i,j)
             buf1(bpter(ib,0)+1) = phis(i,j)

!DIR$ CONCURRENT
             do k=1,plev

                buf1(bpter(ib,k))   = t3  (i,k,j)
                buf1(bpter(ib,k)+1) = u3  (i,k,j)
                buf1(bpter(ib,k)+2) = v3  (i,k,j)
                buf1(bpter(ib,k)+3) = omga(i,k,j)
                buf1(bpter(ib,k)+4) = q3  (i,k,1,j)

                ! convert moist-type constituents from dry to moist mixing ratio

                do m=2,pcnst
                   if (wetq(m)) then  
                      buf1(bpter(ib,k)+3+m) = q3(i,k,m,j)*(1. - q3(i,k,1,j))
                   else
                      buf1(bpter(ib,k)+3+m) = q3(i,k,m,j)
                   endif
                end do

                buf1(bpter(ib,k)+4+pcnst) = pdeld(i,k,j)

             end do

          end do

       end do
!

       call t_barrierf ('sync_blk_to_chk', mpicom)
       call t_startf ('block_to_chunk')
!jjr       call transpose_block_to_chunk(tsize, buf1, buf2, buf2win)
       call transpose_block_to_chunk(tsize, buf1, buf2)
       call t_stopf  ('block_to_chunk')

!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, CPTER, I, K, M)
       do lchnk = begchunk,endchunk
          ncol = phys_state(lchnk)%ncol

          call block_to_chunk_recv_pters(lchnk,pcols,pver+1,tsize,cpter)

          do i=1,ncol
             phys_state(lchnk)%ps   (i)     = buf2(cpter(i,0))
             phys_state(lchnk)%phis (i)     = buf2(cpter(i,0)+1)
          end do

          do k=1,plev

             do i=1,ncol

                phys_state(lchnk)%t    (i,k)   = buf2(cpter(i,k))
                phys_state(lchnk)%u    (i,k)   = buf2(cpter(i,k)+1)
                phys_state(lchnk)%v    (i,k)   = buf2(cpter(i,k)+2)
                phys_state(lchnk)%omega (i,k)   = buf2(cpter(i,k)+3)

                do m=1,pcnst
                   phys_state(lchnk)%q (i,k,m) = buf2(cpter(i,k)+3+m)
                end do

                phys_state(lchnk)%pdeldry(i,k) = buf2(cpter(i,k)+4+pcnst)

             end do

          end do

       end do
       deallocate(bpter)
       deallocate(buf1)
       deallocate(buf2)

    endif

!-----------------------------------------------------------------------
! Fill auxilliary arrays in physics data structure
!-----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, I, K, M, LONS, LATS)

    do lchnk = begchunk,endchunk
       ncol = phys_state(lchnk)%ncol

! pressure arrays
       call plevs0(ncol, pcols, pver, &
                   phys_state(lchnk)%ps,   phys_state(lchnk)%pint,    &
                   phys_state(lchnk)%pmid, phys_state(lchnk)%pdel)

! log(pressure) arrays and Exner function
       do k=1,pver+1
          do i=1,ncol
             phys_state(lchnk)%lnpint(i,k) = log(phys_state(lchnk)%pint(i,k))
          end do
       end do
       do k=1,pver
          do i=1,ncol
             phys_state(lchnk)%rpdel(i,k)  = 1._r8/phys_state(lchnk)%pdel(i,k)
             phys_state(lchnk)%lnpmid(i,k) = log(phys_state(lchnk)%pmid(i,k))
             phys_state(lchnk)%exner (i,k) = (phys_state(lchnk)%pint(i,pver+1) &
                                             / phys_state(lchnk)%pmid(i,k))**cappa
          end do
       end do

! Compute initial geopotential heights
       call geopotential_t (phys_state(lchnk)%lnpint, phys_state(lchnk)%lnpmid  , phys_state(lchnk)%pint  , &
                            phys_state(lchnk)%pmid  , phys_state(lchnk)%pdel    , phys_state(lchnk)%rpdel , &
                            phys_state(lchnk)%t     , phys_state(lchnk)%q(1,1,1), rair,  gravit,  zvir    , &
                            phys_state(lchnk)%zi    , phys_state(lchnk)%zm      , ncol                )

! Compute initial dry static energy, include surface geopotential
       do k = 1, pver
          do i=1,ncol
             phys_state(lchnk)%s(i,k) = cpair*phys_state(lchnk)%t(i,k) &
                                      + gravit*phys_state(lchnk)%zm(i,k) + phys_state(lchnk)%phis(i)
          end do
       end do

! Compute other dry fields in phys_state, using pdeld copied from dynamics above
       call set_state_pdry(phys_state(lchnk),pdeld_calc=.false.)

! Compute energy and water integrals of input state
       call check_energy_timestep_init(phys_state(lchnk), phys_tend(lchnk), pbuf)

    end do

    return
  end subroutine d_p_coupling

!===============================================================================
  subroutine p_d_coupling(phys_state, phys_tend, t2, fu, fv, flx_net, qminus)
!------------------------------------------------------------------------------
! Coupler for converting physics output variables into dynamics input variables
!------------------------------Arguments--------------------------------
    use constituents, only: cnst_get_type_byind

    type(physics_state),intent(in), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend), intent(in), dimension(begchunk:endchunk) :: phys_tend

    real(r8), intent(out) :: t2(plon, beglev:endlev, beglat:endlat)        ! temp tendency
    real(r8), intent(out) :: fu(plon, beglev:endlev, beglat:endlat)        ! u wind tendency
    real(r8), intent(out) :: fv(plon, beglev:endlev, beglat:endlat)        ! v wind tendency
    real(r8), intent(out) :: flx_net(plon,beglat:endlat)          ! net flux
    real(r8), intent(out) :: qminus(plon, beglev:endlev, pcnst, beglat:endlat) ! constituents

!---------------------------Local workspace-----------------------------
    real(r8) :: dummy3(plon,beglev:endlev,beglat:endlat)
    real(r8) :: dummy3xy(beglonxy:endlonxy,plev,beglatxy:endlatxy)
    real(r8) :: t2xy(beglonxy:endlonxy, plev, beglatxy:endlatxy)        ! temp tendency
    real(r8) :: fuxy(beglonxy:endlonxy, plev, beglatxy:endlatxy)        ! u wind tendency
    real(r8) :: fvxy(beglonxy:endlonxy, plev, beglatxy:endlatxy)        ! v ind tendency
    real(r8) :: flx_netxy(beglonxy:endlonxy,beglatxy:endlatxy)          ! net flux
    real(r8) :: qminusxy(beglonxy:endlonxy, plev, pcnst, beglatxy:endlatxy) ! constituents
!
#if (! defined SPMD)
    real(r8) :: buf1(1), buf2(1)     ! transpose buffers
    integer  :: buf1win, buf2win     ! MPI-2 window ids
    integer  :: spmdbuf_siz = 0
    integer  :: block_buf_nrecs = 0
    integer  :: chunk_buf_nrecs = 0
    integer  :: mpicom = 0
    logical  :: local_dp_map=.true. 
#endif

    integer :: i,j,k,m,lchnk,ib         ! indices
    integer :: ncol                  ! number of columns in current chunk
    integer :: lats(pcols)           ! array of latitude indices
    integer :: lons(pcols)           ! array of longitude indices
    integer :: tsize                 ! amount of data per grid point passed to physics
    integer, allocatable, dimension(:,:) :: bpter !jjr
    real(r8), allocatable, dimension(:) :: buf1,buf2

!!    integer :: bpter(plon,0:plev)    ! offsets into block buffer for packing data
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data
    logical :: wetq(pcnst)           ! 'wet' constituent flag
    integer:: blksiz
    type (T_FVDYCORE_STATE), pointer :: dyn_state   !zhh 2013-02-01
    type (T_FVDYCORE_GRID) , pointer :: grid        ! For convenience
!-----------------------------------------------------------------------
!
   dyn_state => get_dyn_state()
   grid => dyn_state%grid     ! For convenience

! Determine which constituents are wet and which are dry
    do m=2,pcnst
       if (cnst_get_type_byind(m).eq.'wet') then  
          wetq(m) = .true.
       else
          wetq(m) = .false.
       endif
    enddo
!-----------------------------------------------------------------------
! copy data from physics data structure to dynamics data structure
!-----------------------------------------------------------------------
    if (local_dp_map) then

!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, I, K, M, LONS, LATS)

       do lchnk = begchunk,endchunk
          ncol = get_ncols_p(lchnk)
          call get_lon_all_p(lchnk, ncol, lons)
          call get_lat_all_p(lchnk, ncol, lats)

          do k=1,plev
             do i=1,ncol
                t2xy(lons(i),k,lats(i))   = phys_tend(lchnk)%dTdt (i,k)
                fuxy(lons(i),k,lats(i))   = phys_tend(lchnk)%dudt (i,k)
                fvxy(lons(i),k,lats(i))   = phys_tend(lchnk)%dvdt (i,k)
                qminusxy(lons(i),k,1,lats(i)) = phys_state(lchnk)%q(i,k,1)
             end do
          end do

          do i=1,ncol
             flx_netxy(lons(i),lats(i))   = phys_tend(lchnk)%flx_net(i)
          end do
          
          ! convert moist-type constituents from moist to dry mixing ratio

          do m=2,pcnst
             if (wetq(m)) then  
                do k=1,plev
                   do i=1,ncol
                      qminusxy(lons(i),k,m,lats(i)) = phys_state(lchnk)%q(i,k,m) /     &
                           (1._r8 - phys_state(lchnk)%q(i,k,1))
                   end do
                end do
             else
                do k=1,plev
                   do i=1,ncol
                      qminusxy(lons(i),k,m,lats(i)) = phys_state(lchnk)%q(i,k,m)
                   end do
                end do
             endif
          end do

       end do

    else

       tsize = 3 + pcnst
       blksiz = (endlatxy-beglatxy+1)*(endlonxy-beglonxy+1)
       allocate(bpter(blksiz,0:plev))
       allocate(buf1(tsize*block_buf_nrecs))
       allocate(buf2(tsize*chunk_buf_nrecs))

!       if (tsize*max(block_buf_nrecs,chunk_buf_nrecs) > spmdbuf_siz) then
!          call endrun ('d_p_coupling: communication buffers (spmdbuf_siz) too small')
!       endif

!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, CPTER, I, K, M)
       do lchnk = begchunk,endchunk
          ncol = get_ncols_p(lchnk)

          call chunk_to_block_send_pters(lchnk,pcols,pver+1,tsize,cpter)

          do i=1,ncol
             buf2(cpter(i,0)) = phys_tend(lchnk)%flx_net(i)
          end do

          do k=1,plev

             do i=1,ncol

                buf2(cpter(i,k))   = phys_tend(lchnk)%dTdt (i,k)
                buf2(cpter(i,k)+1) = phys_tend(lchnk)%dudt (i,k)
                buf2(cpter(i,k)+2) = phys_tend(lchnk)%dvdt (i,k)
                buf2(cpter(i,k)+3) = phys_state(lchnk)%q(i,k,1)

                ! convert moist-type constituents from moist to dry mixing ratio

                do m=2,pcnst
                   if (wetq(m)) then  
                      buf2(cpter(i,k)+2+m) = phys_state(lchnk)%q(i,k,m) /     &
                           (1._r8 - phys_state(lchnk)%q(i,k,1))
                   else
                      buf2(cpter(i,k)+2+m) = phys_state(lchnk)%q(i,k,m)
                   endif
                end do

             end do

          end do

       end do

       call t_barrierf ('sync_chk_to_blk', mpicom)
       call t_startf ('chunk_to_block')
!jjr       call transpose_chunk_to_block(tsize, buf2, buf1, buf1win)
       call transpose_chunk_to_block(tsize, buf2, buf1)
       call t_stopf ('chunk_to_block')

#ifdef OUTER_OMP
!$OMP PARALLEL DO PRIVATE (J, BPTER, I, K, M)
#endif
!!       do j=beglat,endlat

!!          call chunk_to_block_recv_pters(j,plon,plev+1,tsize,bpter)
       call chunk_to_block_recv_pters(iam+1,plon,plev+1,tsize,bpter)
       ib=0
       do j=beglatxy,endlatxy
          do i=beglonxy,endlonxy
             ib=ib+1
             flx_netxy(i,j)=buf1(bpter(ib,0))
             do k=1,plev

                t2xy(i,k,j) = buf1(bpter(ib,k))
                fuxy(i,k,j) = buf1(bpter(ib,k)+1)
                fvxy(i,k,j) = buf1(bpter(ib,k)+2)

                do m=1,pcnst
                   qminusxy(i,k,m,j) = buf1(bpter(ib,k)+2+m)
                end do

             end do

          end do

       end do
       deallocate(bpter)
       deallocate(buf1)
       deallocate(buf2)

    endif

    if( grid%twod_decomp.eq.1) then
#if defined (SPMD)
! Transpose from xy to yz decomposition
      do m=1,pcnst
         do k = 1,plev
            do j = beglatxy,endlatxy
               do i = beglonxy,endlonxy
                  dummy3xy(i,k,j) = qminusxy(i,k,m,j)
               enddo
            enddo
         enddo
         call mp_sendirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,         &
                         grid%ikj_xy_to_yz%RecvDesc, dummy3xy, dummy3,    &
                         modc=grid%modc_dynrun )
         call mp_recvirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,         &
                         grid%ikj_xy_to_yz%RecvDesc, dummy3xy, dummy3,    &
                         modc=grid%modc_dynrun )
         do k=beglev,endlev
            do j=beglat,endlat
               do i=1,plon
                  qminus(i,k,m,j)=dummy3(i,k,j)
               end do
            end do
         end do
      end do
!!
      do k = 1,plev
         do j = beglatxy,endlatxy
            do i = beglonxy,endlonxy
               dummy3xy(i,k,j) =flx_netxy(i,j)
            enddo
         enddo
      enddo
      call mp_sendirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,         &
                      grid%ikj_xy_to_yz%RecvDesc, dummy3xy, dummy3,    &
                      modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,         &
                      grid%ikj_xy_to_yz%RecvDesc, dummy3xy, dummy3,    &
                      modc=grid%modc_dynrun )
      do j=beglat,endlat
         do i=1,plon
            flx_net(i,j)=dummy3(i,beglev,j)    
         end do
      end do
!!
      call mp_sendirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,      &
                      grid%ikj_xy_to_yz%RecvDesc, t2xy, t2,         &
                      modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,      &
                      grid%ikj_xy_to_yz%RecvDesc, t2xy, t2,         &
                      modc=grid%modc_dynrun )
      call mp_sendirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,      &
                      grid%ikj_xy_to_yz%RecvDesc, fuxy, fu,         &
                      modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,      &
                      grid%ikj_xy_to_yz%RecvDesc, fuxy, fu,         &
                      modc=grid%modc_dynrun )
      call mp_sendirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,      &
                      grid%ikj_xy_to_yz%RecvDesc, fvxy, fv,         &
                      modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,      &
                      grid%ikj_xy_to_yz%RecvDesc, fvxy, fv,         &
                      modc=grid%modc_dynrun )

#endif
    else  !jjr need define _yz
        qminus=qminusxy
        flx_net=flx_netxy
        t2=t2xy
        fu=fuxy
        fv=fvxy
    endif


    return
  end subroutine p_d_coupling

!-----------------------------------------------------------------------
! juanxiong he
!-----------------------------------------------------------------------

  subroutine wrf_to_dynamics_tend(phys_state, phys_tend, t2, fu, fv)
!------------------------------------------------------------------------------
! Coupler for converting physics output variables into dynamics input variables
!------------------------------Arguments--------------------------------
    use constituents, only: cnst_get_type_byind

    type(physics_state),intent(in), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend), intent(in), dimension(begchunk:endchunk) :: phys_tend

    real(r8), intent(inout) :: t2(plon, beglev:endlev, beglat:endlat)        ! temp tendency
    real(r8), intent(inout) :: fu(plon, beglev:endlev, beglat:endlat)        ! u wind tendency
    real(r8), intent(inout) :: fv(plon, beglev:endlev, beglat:endlat)        ! v wind tendency

!---------------------------Local workspace-----------------------------
    real(r8) :: t2xy(beglonxy:endlonxy, plev, beglatxy:endlatxy)        ! temp tendency
    real(r8) :: fuxy(beglonxy:endlonxy, plev, beglatxy:endlatxy)        ! u wind tendency
    real(r8) :: fvxy(beglonxy:endlonxy, plev, beglatxy:endlatxy)        ! v ind tendency
    real(r8) :: t2temp(plon, beglev:endlev, beglat:endlat)        ! temp tendency
    real(r8) :: futemp(plon, beglev:endlev, beglat:endlat)        ! u wind tendency
    real(r8) :: fvtemp(plon, beglev:endlev, beglat:endlat)        ! v wind tendency

!
#if (! defined SPMD)
    real(r8) :: buf1(1), buf2(1)     ! transpose buffers
    integer  :: buf1win, buf2win     ! MPI-2 window ids
    integer  :: spmdbuf_siz = 0
    integer  :: block_buf_nrecs = 0
    integer  :: chunk_buf_nrecs = 0
    integer  :: mpicom = 0
    logical  :: local_dp_map=.true. 
#endif

    integer :: i,j,k,m,lchnk,ib         ! indices
    integer :: ncol                  ! number of columns in current chunk
    integer :: lats(pcols)           ! array of latitude indices
    integer :: lons(pcols)           ! array of longitude indices
    integer :: tsize                 ! amount of data per grid point passed to physics
    integer, allocatable, dimension(:,:) :: bpter !jjr
    real(r8), allocatable, dimension(:) :: buf1,buf2

    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data
    logical :: wetq(pcnst)           ! 'wet' constituent flag
    integer:: blksiz
    type (T_FVDYCORE_STATE), pointer :: dyn_state   !zhh 2013-02-01
    type (T_FVDYCORE_GRID) , pointer :: grid        ! For convenience
!-----------------------------------------------------------------------
!
   dyn_state => get_dyn_state()
   grid => dyn_state%grid     ! For convenience

!-----------------------------------------------------------------------
! copy data from physics data structure to dynamics data structure
!-----------------------------------------------------------------------
    if (local_dp_map) then

!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, I, K, M, LONS, LATS)

       do lchnk = begchunk,endchunk
          ncol = get_ncols_p(lchnk)
          call get_lon_all_p(lchnk, ncol, lons)
          call get_lat_all_p(lchnk, ncol, lats)

          do k=1,plev
             do i=1,ncol
                t2xy(lons(i),k,lats(i))   = phys_tend(lchnk)%dTdt (i,k)
                fuxy(lons(i),k,lats(i))   = phys_tend(lchnk)%dudt (i,k)
                fvxy(lons(i),k,lats(i))   = phys_tend(lchnk)%dvdt (i,k)
             end do
          end do

       end do

    else

       tsize = 3 
       blksiz = (endlatxy-beglatxy+1)*(endlonxy-beglonxy+1)
       allocate(bpter(blksiz,0:plev))
       allocate(buf1(tsize*block_buf_nrecs))
       allocate(buf2(tsize*chunk_buf_nrecs))

!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, CPTER, I, K, M)
       do lchnk = begchunk,endchunk
          ncol = get_ncols_p(lchnk)

          call chunk_to_block_send_pters(lchnk,pcols,pver+1,tsize,cpter)

          do k=1,plev

             do i=1,ncol

                buf2(cpter(i,k))   = phys_tend(lchnk)%dTdt (i,k)
                buf2(cpter(i,k)+1) = phys_tend(lchnk)%dudt (i,k)
                buf2(cpter(i,k)+2) = phys_tend(lchnk)%dvdt (i,k)

             end do

          end do

       end do

       call t_barrierf ('sync_chk_to_blk', mpicom)
       call t_startf ('chunk_to_block')
       call transpose_chunk_to_block(tsize, buf2, buf1)
       call t_stopf ('chunk_to_block')

#ifdef OUTER_OMP
!$OMP PARALLEL DO PRIVATE (J, BPTER, I, K, M)
#endif

       call chunk_to_block_recv_pters(iam+1,plon,plev+1,tsize,bpter)
       ib=0
       do j=beglatxy,endlatxy
          do i=beglonxy,endlonxy
             ib=ib+1
             do k=1,plev

                t2xy(i,k,j) = buf1(bpter(ib,k))
                fuxy(i,k,j) = buf1(bpter(ib,k)+1)
                fvxy(i,k,j) = buf1(bpter(ib,k)+2)

             end do

          end do

       end do
       deallocate(bpter)
       deallocate(buf1)
       deallocate(buf2)

    endif

    if( grid%twod_decomp.eq.1) then
#if defined (SPMD)
! Transpose from xy to yz decomposition
      call mp_sendirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,      &
                      grid%ikj_xy_to_yz%RecvDesc, t2xy, t2temp,         &
                      modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,      &
                      grid%ikj_xy_to_yz%RecvDesc, t2xy, t2temp,         &
                      modc=grid%modc_dynrun )
      call mp_sendirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,      &
                      grid%ikj_xy_to_yz%RecvDesc, fuxy, futemp,         &
                      modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,      &
                      grid%ikj_xy_to_yz%RecvDesc, fuxy, futemp,         &
                      modc=grid%modc_dynrun )
      call mp_sendirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,      &
                      grid%ikj_xy_to_yz%RecvDesc, fvxy, fvtemp,         &
                      modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,      &
                      grid%ikj_xy_to_yz%RecvDesc, fvxy, fvtemp,         &
                      modc=grid%modc_dynrun )
      t2=t2+t2temp
      fu=fu+futemp
      fv=fv+fvtemp
#endif
    else  !jjr need define _yz
        t2=t2+t2xy
        fu=fu+fuxy
        fv=fv+fvxy
    endif


    return
  end subroutine wrf_to_dynamics_tend

  subroutine wrf_to_dynamics(phys_state, ps, t3, u3, v3, q3)  ! juanxiong
!------------------------------------------------------------------------------
! Coupler for converting wrf input variables into dynamics input variables
!------------------------------Arguments--------------------------------
    use constituents, only: cnst_get_type_byind

    type(physics_state),intent(in), dimension(begchunk:endchunk) :: phys_state
    real(r8), intent(inout) :: ps(plon, beglat:endlat) ! surface pressure
    real(r8), intent(inout) :: t3(plon, beglev:endlev, beglat:endlat)  ! temperature
    real(r8), intent(inout) :: u3(plon, beglev:endlev, beglat:endlat)  ! u-wind component 
    real(r8), intent(inout) :: v3(plon, beglev:endlev, beglat:endlat)  ! v-wind component
    real(r8), intent(inout) :: q3(plon, beglev:endlev, 1, beglat:endlat) ! constituents

!---------------------------Local workspace-----------------------------
    real(r8) :: dummy31(plon,beglev:endlev,beglat:endlat)
    real(r8) :: dummy32(plon,beglev:endlev,beglat:endlat)
    real(r8) :: dummy33(plon,beglev:endlev,beglat:endlat)
    real(r8) :: dummy34(plon,beglev:endlev,beglat:endlat)
    real(r8) :: dummy31xy(beglonxy:endlonxy,plev,beglatxy:endlatxy)
    real(r8) :: dummy32xy(beglonxy:endlonxy,plev,beglatxy:endlatxy)
    real(r8) :: dummy33xy(beglonxy:endlonxy,plev,beglatxy:endlatxy)
    real(r8) :: dummy34xy(beglonxy:endlonxy,plev,beglatxy:endlatxy)
    real(r8) :: t3xy(beglonxy:endlonxy, plev, beglatxy:endlatxy)        ! temp tendency
    real(r8) :: u3xy(beglonxy:endlonxy, plev, beglatxy:endlatxy)        ! u wind tendency
    real(r8) :: v3xy(beglonxy:endlonxy, plev, beglatxy:endlatxy)        ! v wind tendency
    real(r8) :: q3xy(beglonxy:endlonxy, plev, 1, beglatxy:endlatxy) ! constituents
    real(r8) :: psxy(beglonxy:endlonxy,beglatxy:endlatxy)          ! net flux

!---------------------------Local workspace-----------------------------
#if (! defined SPMD)
    real(r8) :: buf1(1), buf2(1)     ! transpose buffers
    integer  :: buf1win, buf2win     ! MPI-2 window ids
    integer  :: spmdbuf_siz = 0
    integer  :: block_buf_nrecs = 0
    integer  :: chunk_buf_nrecs = 0
    integer  :: mpicom = 0
    logical  :: local_dp_map=.true.
#endif

    integer :: i,j,k,m,lchnk,ib         ! indices
    integer :: ncol                  ! number of columns in current chunk
    integer :: lats(pcols)           ! array of latitude indices
    integer :: lons(pcols)           ! array of longitude indices
    integer :: tsize                 ! amount of data per grid point passed to physics
    integer, allocatable, dimension(:,:) :: bpter 
    real(r8), allocatable, dimension(:) :: buf1,buf2
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data

!---------------------------Local workspace-----------------------------
    integer:: blksiz
    type (T_FVDYCORE_STATE), pointer :: dyn_state   
    type (T_FVDYCORE_GRID) , pointer :: grid       
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! copy data from physics data structure to dynamics data structure
!-----------------------------------------------------------------------

   dyn_state => get_dyn_state()
   grid => dyn_state%grid     ! For convenience

    if (local_dp_map) then

!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, I, K, M, LONS, LATS)

       do lchnk = begchunk,endchunk
          ncol = get_ncols_p(lchnk)
          call get_lon_all_p(lchnk, ncol, lons)
          call get_lat_all_p(lchnk, ncol, lats)

          do k=1,plev
             do i=1,ncol
               if(phys_state(lchnk)%t(i,k).ne.0.0) then
                t3xy(lons(i),k,lats(i))   = phys_state(lchnk)%t (i,k)
                u3xy(lons(i),k,lats(i))   = phys_state(lchnk)%u (i,k)
                v3xy(lons(i),k,lats(i))   = phys_state(lchnk)%v (i,k)
                q3xy(lons(i),k,1,lats(i)) = phys_state(lchnk)%q(i,k,1)
               else
                t3xy(lons(i),k,lats(i))   = 0.0
                u3xy(lons(i),k,lats(i))   = 0.0
                v3xy(lons(i),k,lats(i))   = 0.0
                q3xy(lons(i),k,1,lats(i)) = 0.0
               endif
             end do
          end do

          do i=1,ncol
            if(phys_state(lchnk)%ps(i).ne.0.0) then
              psxy(lons(i),lats(i))   = phys_state(lchnk)%ps(i)
            else
              psxy(lons(i),lats(i))   = 0.0
            end if
          end do

       end do

    else

       tsize = 4
       blksiz = (endlatxy-beglatxy+1)*(endlonxy-beglonxy+1)
       allocate(bpter(blksiz,0:plev))
       allocate(buf1(tsize*block_buf_nrecs))
       allocate(buf2(tsize*chunk_buf_nrecs))

!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, CPTER, I, K, M)
       do lchnk = begchunk,endchunk
          ncol = get_ncols_p(lchnk)

          call chunk_to_block_send_pters(lchnk,pcols,pver+1,tsize,cpter)

          do i=1,ncol
           if(phys_state(lchnk)%ps(i).ne.0)   buf2(cpter(i,0)) = phys_state(lchnk)%ps(i)
          end do


          do k=1,plev

             do i=1,ncol

                buf2(cpter(i,k))   = phys_state(lchnk)%t (i,k)
                buf2(cpter(i,k)+1) = phys_state(lchnk)%u (i,k)
                buf2(cpter(i,k)+2) = phys_state(lchnk)%v (i,k)
                buf2(cpter(i,k)+3) = phys_state(lchnk)%q(i,k,1)

             end do

          end do

       end do

       call t_barrierf ('sync_chk_to_blk', mpicom)
       call t_startf ('chunk_to_block')
       call transpose_chunk_to_block(tsize, buf2, buf1)
       call t_stopf ('chunk_to_block')

#ifdef OUTER_OMP
!$OMP PARALLEL DO PRIVATE (J, BPTER, I, K, M)
#endif
       call chunk_to_block_recv_pters(iam+1,plon,plev+1,tsize,bpter)
       ib=0
       do j=beglatxy,endlatxy
          do i=beglonxy,endlonxy
             ib=ib+1
              if( buf1(bpter(ib,0)) .ne. 0.0 ) then
                psxy(i,j)=buf1(bpter(ib,0))
              else
                psxy(i,j)=0.0
              end if
             do k=1,plev
               if( buf1(bpter(ib,k)) .ne. 0.0 ) then
                t3xy(i,k,j) = buf1(bpter(ib,k))
                u3xy(i,k,j) = buf1(bpter(ib,k)+1)
                v3xy(i,k,j) = buf1(bpter(ib,k)+2)
                q3xy(i,k,1,j) = buf1(bpter(ib,k)+3)
               else
                t3xy(i,k,j) = 0.0
                u3xy(i,k,j) = 0.0
                v3xy(i,k,j) = 0.0
                q3xy(i,k,1,j) = 0.0
               end if
             end do
          end do
       end do
       deallocate(bpter)
       deallocate(buf1)
       deallocate(buf2)
    
 endif

    if( grid%twod_decomp.eq.1) then
#if defined (SPMD)
      ! Transpose from xy to yz decomposition
      ! ps
      do k = 1,plev
         do j = beglatxy,endlatxy
            do i = beglonxy,endlonxy
               dummy31xy(i,k,j) =psxy(i,j)
            enddo
         enddo
      enddo
      call mp_sendirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,         &
                      grid%ikj_xy_to_yz%RecvDesc, dummy31xy, dummy31,    &
                      modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,         &
                      grid%ikj_xy_to_yz%RecvDesc, dummy31xy, dummy31,    &
                      modc=grid%modc_dynrun )
      do j=beglat,endlat
         do i=1,plon
           if(dummy31(i,beglev,j).ne.0) then
            ps(i,j)=dummy31(i,beglev,j)    
           endif
         end do
      end do

      ! t3, u3, v3, q3
         do k = 1,plev
            do j = beglatxy,endlatxy
               do i = beglonxy,endlonxy
                  dummy31xy(i,k,j) = t3xy(i,k,j)
                  dummy32xy(i,k,j) = u3xy(i,k,j)
                  dummy33xy(i,k,j) = v3xy(i,k,j)
                  dummy34xy(i,k,j) = q3xy(i,k,1,j)
               enddo
            enddo
         enddo
         call mp_sendirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,         &
                         grid%ikj_xy_to_yz%RecvDesc, dummy31xy, dummy31,    &
                         modc=grid%modc_dynrun )
         call mp_recvirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,         &
                         grid%ikj_xy_to_yz%RecvDesc, dummy31xy, dummy31,    &
                         modc=grid%modc_dynrun )
         call mp_sendirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,         &
                         grid%ikj_xy_to_yz%RecvDesc, dummy32xy, dummy32,    &
                         modc=grid%modc_dynrun )
         call mp_recvirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,         &
                         grid%ikj_xy_to_yz%RecvDesc, dummy32xy, dummy32,    &
                         modc=grid%modc_dynrun )
         call mp_sendirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,         &
                         grid%ikj_xy_to_yz%RecvDesc, dummy33xy, dummy33,    &
                         modc=grid%modc_dynrun )
         call mp_recvirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,         &
                         grid%ikj_xy_to_yz%RecvDesc, dummy33xy, dummy33,    &
                         modc=grid%modc_dynrun )
         call mp_sendirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,         &
                         grid%ikj_xy_to_yz%RecvDesc, dummy34xy, dummy34,    &
                         modc=grid%modc_dynrun )
         call mp_recvirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,         &
                         grid%ikj_xy_to_yz%RecvDesc, dummy34xy, dummy34,    &
                         modc=grid%modc_dynrun )
         do k=beglev,endlev
            do j=beglat,endlat
               do i=1,plon
                  if(dummy31(i,k,j).ne.0.0) then
                  t3(i,k,j)=dummy31(i,k,j)
                  u3(i,k,j)=dummy32(i,k,j)
                  v3(i,k,j)=dummy33(i,k,j)
                  q3(i,k,1,j)=dummy34(i,k,j)
                  end if
               end do
            end do
         end do 

#endif
    else  

       do lchnk = beglat,endlat
          do k = beglev,endlev
             do i = 1,plon
               if(t3xy(i,k,lchnk).ne.0.0) then
                t3(i,k,lchnk)   = t3xy (i,k,lchnk)
                u3(i,k,lchnk)   = u3xy (i,k,lchnk)
                v3(i,k,lchnk)   = v3xy (i,k,lchnk)
                q3(i,k,1,lchnk)   = q3xy (i,k,1,lchnk)
               endif
             end do
          end do

          do i=1,plon
            if(psxy(i,lchnk).ne.0.0)   then
              ps(i,lchnk)  =  psxy(i,lchnk)
            end if
          end do

       end do

    endif
    return
  end subroutine wrf_to_dynamics

!===============================================================================
  subroutine dynamics_to_wrf(phys_state, ps, t3, u3, v3, q3, pdeld)
!------------------------------------------------------------------------------
! Coupler for converting physics output variables into dynamics input variables
!------------------------------Arguments--------------------------------
    use constituents, only: cnst_get_type_byind
    use physconst,     only: cappa

    type(physics_state),intent(inout), dimension(begchunk:endchunk) :: phys_state
    real(r8), intent(in) :: t3(plon, beglev:endlev, beglat:endlat)        ! temp tendency
    real(r8), intent(in) :: u3(plon, beglev:endlev, beglat:endlat)        ! u wind tendency
    real(r8), intent(in) :: v3(plon, beglev:endlev, beglat:endlat)        ! v wind tendency
    real(r8), intent(in) :: pdeld(plon, beglev:endlev, beglat:endlat)     ! pressure thickness
    real(r8), intent(in) :: ps(plon, beglat:endlat)          ! net flux
    real(r8), intent(in) :: q3(plon, beglev:endlev, 1, beglat:endlat) ! constituents

!---------------------------Local workspace-----------------------------
    real(r8) :: dummy3(plon,beglev:endlev,beglat:endlat)
    real(r8) :: dummy3xy(beglonxy:endlonxy,plev,beglatxy:endlatxy)
    real(r8) :: t3xy(beglonxy:endlonxy, plev, beglatxy:endlatxy)        ! temp tendency
    real(r8) :: u3xy(beglonxy:endlonxy, plev, beglatxy:endlatxy)        ! u wind tendency
    real(r8) :: v3xy(beglonxy:endlonxy, plev, beglatxy:endlatxy)        ! v ind tendency
    real(r8) :: psxy(beglonxy:endlonxy,beglatxy:endlatxy)          ! net flux
    real(r8) :: q3xy(beglonxy:endlonxy, plev, 1, beglatxy:endlatxy) ! constituents
    real(r8) :: pdeldxy (beglonxy:endlonxy,plev,beglatxy:endlatxy)
    real(r8) :: tem2(pcols,pver) ! temporary workspace
    real(r8) :: ftem(pcols,pver) ! temporary workspace
!
#if (! defined SPMD)
    real(r8) :: buf1(1), buf2(1)     ! transpose buffers
    integer  :: spmdbuf_siz = 0
    integer  :: block_buf_nrecs = 0
    integer  :: chunk_buf_nrecs = 0
    integer  :: mpicom = 0
    logical  :: local_dp_map=.true. 
#endif

    integer :: i,j,k,m,lchnk,ib         ! indices
    integer :: ncol                  ! number of columns in current chunk
    integer :: lats(pcols)           ! array of latitude indices
    integer :: lons(pcols)           ! array of longitude indices
    integer :: tsize                 ! amount of data per grid point passed to physics
    integer, allocatable, dimension(:,:) :: bpter !jjr
    real(r8), allocatable, dimension(:) :: buf1,buf2

    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data
    logical :: wetq(pcnst)           ! 'wet' constituent flag
    integer:: blksiz
    type (T_FVDYCORE_STATE), pointer :: dyn_state   !zhh 2013-02-01
    type (T_FVDYCORE_GRID) , pointer :: grid        ! For convenience
!-----------------------------------------------------------------------
!
   dyn_state => get_dyn_state()
   grid => dyn_state%grid     ! For convenience
   
   if (grid%twod_decomp .eq. 1) then
#if defined( SPMD )
      do k = beglev,endlev
         do j = beglat,endlat
            do i = 1,plon
               dummy3(i,k,j) = ps(i,j)
            enddo
         enddo
      enddo
      call mp_sendirr(grid%commxy, grid%ijk_yz_to_xy%SendDesc,         &
                      grid%ijk_yz_to_xy%RecvDesc, dummy3, dummy3xy,    &
                      modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%ijk_yz_to_xy%SendDesc,         &
                      grid%ijk_yz_to_xy%RecvDesc, dummy3, dummy3xy,    &
                      modc=grid%modc_dynrun )
      do j = beglatxy,endlatxy
         do i = beglonxy,endlonxy
            psxy(i,j) = dummy3xy(i,1,j)
         enddo
      enddo

      call mp_sendirr(grid%commxy, grid%ikj_yz_to_xy%SendDesc,         &
                      grid%ikj_yz_to_xy%RecvDesc, u3(:,:,:),      &
                      u3xy(:,:,:), modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%ikj_yz_to_xy%SendDesc,         &
                      grid%ikj_yz_to_xy%RecvDesc, u3(:,:,:),      &
                      u3xy(:,:,:), modc=grid%modc_dynrun )
!
      call mp_sendirr(grid%commxy, grid%ikj_yz_to_xy%SendDesc,         &
                      grid%ikj_yz_to_xy%RecvDesc, v3(:,:,:),      &
                      v3xy(:,:,:), modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%ikj_yz_to_xy%SendDesc,         &
                      grid%ikj_yz_to_xy%RecvDesc, v3(:,:,:),      &
                      v3xy(:,:,:), modc=grid%modc_dynrun )
!
      call mp_sendirr(grid%commxy, grid%ikj_yz_to_xy%SendDesc,         &
                      grid%ikj_yz_to_xy%RecvDesc, t3(:,:,:),      &
                      t3xy(:,:,:), modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%ikj_yz_to_xy%SendDesc,         &
                      grid%ikj_yz_to_xy%RecvDesc, t3(:,:,:),      &
                      t3xy(:,:,:), modc=grid%modc_dynrun )

      call mp_sendirr(grid%commxy, grid%ikj_yz_to_xy%SendDesc,         &
                      grid%ikj_yz_to_xy%RecvDesc, pdeld(:,:,:),   &
                      pdeldxy(:,:,:), modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%ikj_yz_to_xy%SendDesc,         &
                      grid%ikj_yz_to_xy%RecvDesc, pdeld(:,:,:),   &
                      pdeldxy(:,:,:), modc=grid%modc_dynrun )

         do k = beglev,endlev
            do j = beglat,endlat
               do i =1,plon
                  dummy3(i,k,j) = q3(i,k,1,j)
               enddo
            enddo
         enddo
         call mp_sendirr(grid%commxy, grid%ijk_yz_to_xy%SendDesc,         &
                         grid%ijk_yz_to_xy%RecvDesc, dummy3, dummy3xy,    &
                         modc=grid%modc_dynrun )
         call mp_recvirr(grid%commxy, grid%ijk_yz_to_xy%SendDesc,         &
                         grid%ijk_yz_to_xy%RecvDesc, dummy3, dummy3xy,    &
                         modc=grid%modc_dynrun )
         do k=1,plev
            do j = beglatxy,endlatxy
               do i = beglonxy,endlonxy
                  q3xy(i,k,1,j) = dummy3xy(i,k,j)
               enddo
            end do
         enddo

#endif
   else !jjr need define _xy
      psxy(:,:)=ps(:,:)
      t3xy(:,:,:)=t3(:,:,:)
      u3xy(:,:,:)=u3(:,:,:)
      v3xy(:,:,:)= v3(:,:,:)
      pdeldxy(:,:,:)= pdeld(:,:,:)
      q3xy(:,:,:,:)=q3(:,:,:,:)
    end if

!-----------------------------------------------------------------------
! copy data from dynamics data structure to physics data structure
!-----------------------------------------------------------------------
    if (local_dp_map) then

!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, I, K, M, LONS, LATS)
       do lchnk = begchunk,endchunk
          ncol = phys_state(lchnk)%ncol
          call get_lon_all_p(lchnk, ncol, lons)
          call get_lat_all_p(lchnk, ncol, lats)
          do i=1,ncol
             phys_state(lchnk)%ps   (i)     = psxy  (lons(i),lats(i))
          end do
  
          do k=1,plev
             do i=1,ncol
                phys_state(lchnk)%t    (i,k)   = t3xy(lons(i),k,lats(i))
                phys_state(lchnk)%u    (i,k)   = u3xy(lons(i),k,lats(i))
                phys_state(lchnk)%v    (i,k)   = v3xy(lons(i),k,lats(i))
                phys_state(lchnk)%q  (i,k,1)   = q3xy(lons(i),k,1,lats(i))
                phys_state(lchnk)%pdel (i,k)   = pdeldxy(lons(i),k,lats(i))
             end do
          end do

       end do

    else

       tsize = 5  

       blksiz = (endlatxy-beglatxy+1)*(endlonxy-beglonxy+1)
       allocate(bpter(blksiz,0:plev))
       allocate(buf1(tsize*block_buf_nrecs))
       allocate(buf2(tsize*chunk_buf_nrecs))

#ifdef OUTER_OMP
!$OMP PARALLEL DO PRIVATE (J, BPTER, I, K, M)
#endif

       call block_to_chunk_send_pters(iam+1,blksiz,plev+1,tsize,bpter)

!DIR$ CONCURRENT
         ib = 0
         do j=beglatxy,endlatxy
          do i=beglonxy,endlonxy
            ib = ib + 1
             buf1(bpter(ib,0))   = psxy(i,j)

!DIR$ CONCURRENT
             do k=1,plev
                buf1(bpter(ib,k))   = t3xy(i,k,j)
                buf1(bpter(ib,k)+1) = u3xy(i,k,j)
                buf1(bpter(ib,k)+2) = v3xy(i,k,j)
                buf1(bpter(ib,k)+3) = q3xy(i,k,1,j)
                buf1(bpter(ib,k)+4) = pdeldxy(i,k,j)
             end do
          end do
         end do

       call t_barrierf ('sync_blk_to_chk', mpicom)
       call t_startf ('block_to_chunk')
       call transpose_block_to_chunk(tsize, buf1, buf2)
       call t_stopf  ('block_to_chunk')

!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, CPTER, I, K, M)
       do lchnk = begchunk,endchunk
          ncol = phys_state(lchnk)%ncol

          call block_to_chunk_recv_pters(lchnk,pcols,pver+1,tsize,cpter)

          do i=1,ncol
             phys_state(lchnk)%ps   (i)     = buf2(cpter(i,0))
          end do

          do k=1,plev
             do i=1,ncol

                phys_state(lchnk)%t    (i,k)   = buf2(cpter(i,k))
                phys_state(lchnk)%u    (i,k)   = buf2(cpter(i,k)+1)
                phys_state(lchnk)%v    (i,k)   = buf2(cpter(i,k)+2)
                phys_state(lchnk)%pdel (i,k)   = buf2(cpter(i,k)+2)
                phys_state(lchnk)%q (i,k,1) = buf2(cpter(i,k)+3)

             end do
          end do

       end do
       deallocate(bpter)
       deallocate(buf1)
       deallocate(buf2)

    endif

!-----------------------------------------------------------------------
! Fill auxilliary arrays in physics data structure
!-----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, I, K, M, LONS, LATS)

    do lchnk = begchunk,endchunk
       ncol = phys_state(lchnk)%ncol

! pressure arrays
       call plevs0(ncol, pcols, pver, &
                   phys_state(lchnk)%ps,   phys_state(lchnk)%pint,    &
                   phys_state(lchnk)%pmid, phys_state(lchnk)%pdel)

! log(pressure) arrays and Exner function
       do k=1,pver+1
          do i=1,ncol
             phys_state(lchnk)%lnpint(i,k) = log(phys_state(lchnk)%pint(i,k))
          end do
       end do
       do k=1,pver
          do i=1,ncol
             phys_state(lchnk)%rpdel(i,k)  = 1._r8/phys_state(lchnk)%pdel(i,k)
             phys_state(lchnk)%lnpmid(i,k) = log(phys_state(lchnk)%pmid(i,k))
             phys_state(lchnk)%exner (i,k) = (phys_state(lchnk)%pint(i,pver+1) &
                                             / phys_state(lchnk)%pmid(i,k))**cappa
          end do
       end do

! Compute initial geopotential heights
       call geopotential_t (phys_state(lchnk)%lnpint, phys_state(lchnk)%lnpmid  , phys_state(lchnk)%pint  , &
                            phys_state(lchnk)%pmid  , phys_state(lchnk)%pdel    , phys_state(lchnk)%rpdel , &
                            phys_state(lchnk)%t     , phys_state(lchnk)%q(1,1,1), rair,  gravit,  zvir    , &
                            phys_state(lchnk)%zi    , phys_state(lchnk)%zm      , ncol                )

! Compute RH
       call aqsat (phys_state(lchnk)%t, phys_state(lchnk)%pmid, tem2, ftem, pcols, &
                   ncol, pver, 1, pver    )
       phys_state(lchnk)%rh(:ncol,:) = phys_state(lchnk)%q(:ncol,:,1)/ftem(:ncol,:)*100._r8

    end do

    return
  end subroutine dynamics_to_wrf

end module dp_coupling
