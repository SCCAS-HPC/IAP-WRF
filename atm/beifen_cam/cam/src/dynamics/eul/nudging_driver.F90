module nudging_driver

! writed by juanxiong he

   use shr_kind_mod, only: r8 => shr_kind_r8
   integer :: fdda, fdda_sfc,  vert_int, data_src, data_sfcsrc
   integer :: analysis_interval ! nudging time interval 
   real :: guv, gt ,gq ! nudging coefficients 
   integer :: if_zfac_uv, if_zfac_t, if_zfac_q ! variables  take part in nudging 
   integer :: k_zfac_uv, k_zfac_t, k_zfac_q ! the top vertical level take part in nudging 
   integer :: k_zbot_uv, k_zbot_t, k_zbot_q ! the bottom vertical level take part in nudging 

contains

  subroutine nudging(phys_state,phys_tend,nudging_state,nudging_tend,&
                     nnstep,nstep)
   use phys_grid,        only: get_ncols_p  !juanxiong he
   use physics_types,     only: physics_state, physics_tend !juanxiong he
   use ppgrid,           only: begchunk, endchunk, pver ! juanxiong he
   use dp_coupling,      only: uv_to_physics,physics_to_uv, &
                               tq_to_physics,physics_to_tq ! juanxiong he
   use pmgrid,           only: plev, plat, plevp, plon, beglat, endlat ! juanxiong he
   use prognostics,      only: ps, u3, v3, t3, q3, div, vort, n3m2, n3m1 ! juanxiong
   use time_manager,      only: get_step_size, is_first_step ! juanxiong
   use perf_mod

   type(physics_state), pointer, intent(inout) :: phys_state(:), nudging_state(:)
   type(physics_tend), pointer, intent(inout) :: phys_tend(:), nudging_tend(:)
   integer,intent(inout) :: nnstep,nstep

   integer :: i,k,ncols !juanxiong he
   real(r8) :: timestep,dt

   call t_startf ('cam_nudging')

   timestep = get_step_size()
   dt = 2.0_r8*timestep
   !
   ! If initial time step adjust dt
   !
   if (is_first_step()) dt = timestep

   ! phys_state the previous timestep takes part in the physical schemes
   ! it needs nudging 
   do lchnk = begchunk,endchunk  !juanxiong he
   call grid_nudging(nudging_state(lchnk),nudging_tend(lchnk), &
                     dt, analysis_interval, nnstep, &
                     if_zfac_uv, if_zfac_t, if_zfac_q, &
                     k_zfac_uv, k_zfac_t, k_zfac_q, &
                     guv, gt, gq)
   end do 

   dt = timestep 
   call tq_to_physics (ps(1,beglat,n3m1), t3(1,1,beglat,n3m1), &
                       q3(1,1,1,beglat,n3m1),nudging_state)
   call uv_to_physics (ps(1,beglat,n3m1), u3(1,1,beglat,n3m1), &
                       v3(1,1,beglat,n3m1), nudging_state)
   do lchnk = begchunk,endchunk  
   call grid_nudging(nudging_state(lchnk), nudging_tend(lchnk), &
                     dt, analysis_interval, nstep, &
                     if_zfac_uv, if_zfac_t, if_zfac_q, &
                     k_zfac_uv, k_zfac_t, k_zfac_q, &
                     guv, gt, gq)
   end do 
   call physics_to_uv (nudging_state, ps(1,beglat,n3m1), u3(1,1,beglat,n3m1), &
                       v3(1,1,beglat,n3m1))
   call physics_to_tq (nudging_state, ps(1,beglat,n3m1), t3(1,1,beglat,n3m1), &
                       q3(1,1,1,beglat,n3m1))
   call t_stopf  ('cam_nudging')
  end subroutine nudging

  subroutine grid_nudging(state,tend,dt,analysis_interval,nstep,&
                          if_zfac_uv, if_zfac_t, if_zfac_q,&
                          k_zfac_uv,  k_zfac_t,  k_zfac_q,&
                          guv, gt, gq)

    use ppgrid,       only: pcols, pver
    use physics_types, only: physics_state, physics_tend, physics_ptend
    use geopotential, only: geopotential_dse
    use physconst,    only: zvir, gravit, cpair, rair

    implicit none
    type(physics_state), intent(inout) :: state
    type(physics_tend), intent(inout) :: tend

    real(r8), intent(in) :: dt
    integer, intent(in) :: nstep
    integer, intent(in) :: analysis_interval
    integer, intent(in) :: if_zfac_uv, if_zfac_t, if_zfac_q
    integer, intent(in) :: k_zfac_uv,  k_zfac_t,  k_zfac_q
    real, intent(in)  :: guv, gt, gq

    ! local
    integer :: ierr, lchnk
    integer :: i,j,k
    type(physics_ptend)  :: ptend   ! Parameterization tendencies

    integer :: ncol                                ! number of columns
    real(r8), dimension( :, :, : ), allocatable :: wpbl  ! 1: u, 2: v, 3: t, 4: q
    real(r8), dimension( :, : ), allocatable :: wzfac ! 1: u, 2: v, 3: t, 4: q
    real(r8)    :: coef, tfac
    real(r8)    :: val_analysis
    real(r8)    :: xtime, xtime_old, xtime_new
    real(r8)    :: lats,latn,lonw,lone

   ! time
   xtime = dt*(nstep+1)
   xtime_old = FLOOR(xtime/analysis_interval) * analysis_interval * 1.0_r8
   xtime_new = xtime_old + analysis_interval * 1.0_r8
   coef = (xtime-xtime_old)/(xtime_new-xtime_old)

   tfac=1.0_r8

   lats=-27*3.1415926/180.0*1.0_8
   latn=27*3.1415926/180.0*1.0_8
   lonw=51*3.1415926/180.0*1.0_8
   lone=158*3.1415926/180.0*1.0_8

   ! vertical 
   allocate( wzfac(1:pver,4))
   
   wzfac = 0.0_r8
   do k=1,pver
      IF( if_zfac_uv == 1 ) THEN
       IF( k <= k_zfac_uv   ) wzfac(k, 1:2) = 0.0_r8
       IF( k == k_zfac_uv+1 ) wzfac(k, 1:2) = 0.1_r8
       IF( k >  k_zfac_uv+1 ) wzfac(k, 1:2) = 1.0_r8
      ENDIF

      IF( if_zfac_t == 1 ) THEN
       IF( k <= k_zfac_t   ) wzfac(k, 3) = 0.0_r8
       IF( k == k_zfac_t+1 ) wzfac(k, 3) = 0.1_r8
       IF( k >  k_zfac_t+1 ) wzfac(k, 3) = 1.0_r8
      ENDIF

      IF( if_zfac_q == 1 ) THEN
       IF( k <= k_zfac_q   ) wzfac(k, 4) = 0.0_r8
       IF( k == k_zfac_q+1 ) wzfac(k, 4) = 0.1_r8
       IF( k >  k_zfac_q+1 ) wzfac(k, 4) = 1.0_r8
      ENDIF
   enddo

    wzfac((pver-k_zbot_uv):pver,:) = 0.0  ! the lowest layer do not take part in the nudging

   ! Set chunk id, number of columns, and coordinates
     ncol = state%ncol

     allocate( wpbl(1:ncol,1:pver,4))
     wpbl(:,:,:) = 1.0_r8

   ! Compute 3-D nudging tendencies for u, v, t and q
     do k=1,pver
     do i=1,ncol

!     if(.not.(state%lat(i).le.latn.and.state%lat(i).ge.lats.and. &
!        state%lon(i).le.lone.and.state%lon(i).ge.lonw)) then

      val_analysis = state%u_ndg_old(i,k) *( 1.0_r8 - coef ) +  state%u_ndg_new(i,k) * coef
      ptend%rundgdten(i,k) = guv * wpbl(i,k,1) * wzfac(k,1) * tfac * &
                            ( val_analysis -  state%u(i,k) )
      state%u(i,k) = state%u(i,k) + ptend%rundgdten(i,k) * dt

      val_analysis =  state%v_ndg_old(i,k) *( 1.0_r8 - coef ) +  state%v_ndg_new(i,k) * coef
      ptend%rvndgdten(i,k) =  guv * wpbl(i,k,2) * wzfac(k,2) * tfac * &
                     ( val_analysis -  state%v(i,k) )
      state%v(i,k) = state%v(i,k) + ptend%rvndgdten(i,k) * dt

      ! cam doesn't renew state%t but state%s and tend%dtdt, 
      ! so need to deduce the temperature then get the nudging tendecny of t
      state%t(i,k) = (state%s(i,k)-gravit*state%zm(i,k)-state%phis(i))/cpair
      val_analysis =  state%t_ndg_old(i,k) *( 1.0_r8 - coef ) +  state%t_ndg_new(i,k) * coef
      ptend%rtndgdten(i,k) =   gt * wpbl(i,k,3) * wzfac(k,3) * tfac * &
                          ( val_analysis -  state%t(i,k) )
      state%t(i,k) = state%t(i,k) + ptend%rtndgdten(i,k) * dt

      if(state%q_ndg_old(i,k).gt.0.and. state%q_ndg_new(i,k).gt.0) then
      val_analysis =  state%q_ndg_old(i,k) *( 1.0_r8 - coef ) +  state%q_ndg_new(i,k) * coef
      ptend%rqndgdten(i,k) =  gq * wpbl(i,k,4) * wzfac(k,4) * tfac * &
                          ( val_analysis -  state%q(i,k,1) )
      state%q(i,k,1) = state%q(i,k,1) + ptend%rqndgdten(i,k) * dt
      end if

!      endif
     enddo
     enddo

   deallocate( wpbl)
   deallocate( wzfac)

1000 continue

 end subroutine grid_nudging

end module nudging_driver
