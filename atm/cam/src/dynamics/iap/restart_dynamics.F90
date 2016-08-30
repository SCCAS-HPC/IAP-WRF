module restart_dynamics
!--------------------------------------------------------------------
! 
! Purpose: Read and write restart file
!
! Modified by ZhangHe, 2007.5.25
! Update : ZhangHe, 2008.6.5
!          ZhangHe, 2008.6.21
!          ZhangHe, 2008.7.21, writen the IAP dynamical fields from 
!                              (beglatdyn --> endlatdyn) to (beglat --> endlat) 
!          Juanxiong He, 201008, adapted to cesm	
! Reviewed: Zhang He, 2011-12-16
! Update: ZhangHe, 2011-12-17, added calling sub. initialize_IAPprog
!         ZhangHe, 2011-12-18, dyn_init() --> dyn_init(filename)
!         ZhangHe, 2012-01-19
!         ZhangHe, 2012-01-20, v4=qminus --> v5=q3
!         ZhangHe, 2012-02-09
!         ZhangHe, 2012-11-08, removed sub. set_r_dynmvar, sub. set_r_dynivar
!         ZhangHe, 2012-11-09, added GHS
!         ZhangHe, 2013-01-30, write variables as xy decomposition 
!         ZhangHe, 2013-03-14, moved yz2D_xy2D from sub. init_restart_varlist to sub. init_restart_dynamics
!                              revised the error in read_restart_dynamics
!         ZhangHe, 2013-03-21, removed n3m2
!         ZhangHe, 2013-03-29, removed redundant variables in use only
!         ZhangHe, 2013-04-15, added U and V
!--------------------------------------------------------------------
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use pio, only : var_desc_t, file_desc_t, pio_double, pio_unlimited, pio_def_var, &
	          pio_def_dim, io_desc_t, pio_offset, pio_put_var, pio_write_darray, &
                  pio_setdebuglevel, pio_setframe, pio_initdecomp, pio_freedecomp, &
                  pio_read_darray, pio_inq_varid, pio_get_var
                  
  use constituents, only: pcnst
  use prognostics,  only: u3, v3, t3, q3, pdeld, ps, phis, omga, ptimelevels
  use prognostics,  only: phisxy, GHSxy, psxy, omgaxy, lammpxy, phimpxy, sigmpxy,    &
                          WSxy, u3xy, v3xy, t3xy, pdeldxy, qfcstxy, q3xy, Uxy, Vxy       !zhh 2013-01-30
  use scanslt,      only: lammp, phimp, sigmp, qfcst
  use IAP_prog,     only: WS, GHS, U, V, initialize_IAPprog
   
#if ( defined BFB_CAM_SCAM_IOP )
  use iop,             only: dqfx3sav,divq3dsav,divt3dsav,t2sav,betasav
#endif
  use cam_logfile,  only: iulog
  use spmd_utils,   only: masterproc, iam

  implicit none
  private
  public :: read_restart_dynamics, init_restart_dynamics, write_restart_dynamics

  integer, parameter :: namlen=16

  type restart_var_t
     real(r8), pointer :: v1d(:)
     real(r8), pointer :: v2d(:,:)
     real(r8), pointer :: v3d(:, :, :)
     real(r8), pointer :: v4d(:, :, :, :)
     real(r8), pointer :: v5d(:, :, :, :, :)

     type(var_desc_t), pointer  :: vdesc
     integer           :: ndims
     integer           :: timelevels
     character(len=namlen) :: name
  end type restart_var_t
#if ( defined BFB_CAM_SCAM_IOP )
!!  integer, parameter :: restartvarcnt = 22    
  integer, parameter :: restartvarcnt = 19+pcnst*2    ! zhh   
#else
!!  integer, parameter :: restartvarcnt = 17  
  integer, parameter :: restartvarcnt = 14+pcnst*2     ! zhh 
#endif
  
  type(var_desc_t) :: timedesc, tmass0desc, fixmasdesc, hw1desc, hw2desc, hw3desc, alphadesc

  type(restart_var_t) :: restartvars(restartvarcnt)
  logical ::  restart_varlist_initialized=.false.
!
  
CONTAINS

  subroutine set_r_var(name, timelevels, index, v1, v2, v3, v4, v5)
    use abortutils,      only: endrun

    character(len=*), intent(in) :: name
    integer, intent(in) :: timelevels, index
    real(r8), target, optional :: v1(:), v2(:,:), v3(:,:,:), v4(:,:,:,:), v5(:,:,:,:,:)

    restartvars(index)%name=name
    restartvars(index)%timelevels = timelevels
    if(present(v1)) then
       restartvars(index)%ndims = 1
       restartvars(index)%v1d => v1
    else if(present(v2)) then
       restartvars(index)%ndims = 2
       restartvars(index)%v2d => v2
    else if(present(v3)) then
       restartvars(index)%ndims = 3
       restartvars(index)%v3d => v3
    else if(present(v4)) then
       restartvars(index)%ndims = 4
       restartvars(index)%v4d => v4
    else if(present(v5)) then
       restartvars(index)%ndims = 5
       restartvars(index)%v5d => v5
    else
       call endrun('bad ndims in call to set_r_var')
    end if
    allocate(restartvars(index)%vdesc)

  end subroutine set_r_var

!=====================================================================================  
  subroutine init_restart_varlist()
    use abortutils,      only: endrun
    use pmgrid,       only: beglat, endlat
    use constituents, only: cnst_name   !zhh 2013-01-18

    character(len=21) :: cnst_name_fcst(pcnst)     ! constituent names
    integer :: vcnt=1
    integer :: i,j,m
!-------------------------------------------------------------------------------------

! Should only be called once
    if(restart_varlist_initialized) return
    restart_varlist_initialized=.true.
!
    do m = 1, pcnst
       cnst_name_fcst(m) = trim(cnst_name(m))//'_fcst'
    end do

! prognostics

    call set_r_var('PHIS', 1, vcnt, v2=phisxy )
    
    vcnt=vcnt+1
    call set_r_var('GHS', 1, vcnt, v2=GHSxy )

    vcnt=vcnt+1
    call set_r_var('OMGA', 1, vcnt, v3=omgaxy )

    vcnt=vcnt+1
    call set_r_var('WS', 1, vcnt, v3=WSxy )

    vcnt=vcnt+1
    call set_r_var('U', 1, vcnt, v3=Uxy )

    vcnt=vcnt+1
    call set_r_var('V', 1, vcnt, v3=Vxy )

    vcnt=vcnt+1
    call set_r_var('PS',  ptimelevels,  vcnt, v3=psxy )
    
    vcnt=vcnt+1
    call set_r_var('U3',  ptimelevels,  vcnt, v4=u3xy )

    vcnt=vcnt+1
    call set_r_var('V3',  ptimelevels,  vcnt, v4=v3xy )
    
    vcnt=vcnt+1
    call set_r_var('T3',  ptimelevels,  vcnt, v4=t3xy )
    
    do m = 1, pcnst
       vcnt=vcnt+1
       call set_r_var(cnst_name(m),  ptimelevels,  vcnt, v4=q3xy(:,:,:,:,m) )
    end do

    vcnt=vcnt+1
    call set_r_var('PDELD', ptimelevels, vcnt, v4=pdeldxy )
    
! slt
	
    vcnt=vcnt+1
    call set_r_var('LAMMP', 1, vcnt, v3=lammpxy )
    
    vcnt=vcnt+1
    call set_r_var('PHIMP', 1, vcnt, v3=phimpxy )
    
    vcnt=vcnt+1
    call set_r_var('SIGMP', 1, vcnt, v3=sigmpxy )

    do m = 1, pcnst
       vcnt=vcnt+1
       call set_r_var(cnst_name_fcst(m), 1, vcnt, v3=qfcstxy(:,:,:,m) )    
    end do

#if ( defined BFB_CAM_SCAM_IOP )
!
! Write scam values
!
    vcnt=vcnt+1
    call set_r_var('DQFX', 1, vcnt, v4=dqfx3sav )

    vcnt=vcnt+1
    call set_r_var('DIVQ', 1, vcnt, v4=divq3dsav )

    vcnt=vcnt+1
    call set_r_var('DIVT', 1, vcnt, v3=divt3dsav )

    vcnt=vcnt+1
    call set_r_var('T2', 1, vcnt, v3=t2sav )

    vcnt=vcnt+1
    call set_r_var('BETA', 1, vcnt, v1=betasav )

#endif    

    if(vcnt.ne.restartvarcnt) then
       write(iulog,*) 'vcnt= ',vcnt, ' restartvarcnt=',restartvarcnt
       call endrun('bad restartvarcnt')
    end if

  end subroutine init_restart_varlist

!================================================================================
  subroutine init_restart_dynamics(File, hdimids, vdimids, dyn_out)
    use dyn_comp, only  : dyn_export_t
    use dyn_grid, only : get_horiz_grid_dim_d
    use dyn_comp,     only: yz2D_xy2D   !zhh 2013-01-21
    use pmgrid,       only: beglat, endlat, plat    !debug
    use prognostics,     only: n3m1, n3     !debug
    !
    ! Input arguments
    !
    type(File_desc_t), intent(in) :: File     ! Unit number
    type(dyn_export_t), intent(in) :: dyn_out ! Not used in eul dycore
    integer, intent(in) :: vdimids(2)
    integer, pointer :: hdimids(:)
    character(len=namlen) :: name

    integer :: alldims(4), alldims2d(3), qdims(5)
    integer :: timelevels_dimid, i, hdim1, hdim2, ierr
    type(var_desc_t), pointer :: vdesc
    integer :: ndims, timelevels
    integer :: j, jd   !debug
!---------------------------------------------------------------------------------    
!
    call get_horiz_grid_dim_d(hdim1, hdim2)

    allocate(hdimids(2))
    ierr = PIO_Def_Dim(File, 'lon',hdim1, hdimids(1))

    ierr = PIO_Def_Dim(File, 'lat',hdim2, hdimids(2))

    ierr = PIO_Def_Dim(File,'timelevels',PIO_UNLIMITED,timelevels_dimid)

    ierr = PIO_Def_Dim(File,'pcnst',pcnst, qdims(4))

    ierr = PIO_Def_Var(File, 'time', pio_double, (/timelevels_dimid/), timedesc)

    ierr = PIO_Def_var(File, 'tmass0', pio_double, tmass0desc)
    ierr = PIO_Def_var(File, 'fixmas', pio_double, fixmasdesc)
    ierr = PIO_Def_var(File, 'hw1', pio_double, qdims(4:4), hw1desc)
    ierr = PIO_Def_var(File, 'hw2', pio_double, qdims(4:4), hw2desc)
    ierr = PIO_Def_var(File, 'hw3', pio_double, qdims(4:4), hw3desc)
    ierr = PIO_Def_var(File, 'alpha', pio_double, qdims(4:4), alphadesc)

! slt and prognostics
    alldims(1:2) = hdimids(1:2)
    alldims(3) = vdimids(1)
    alldims(4) = timelevels_dimid

    alldims2d(1:2) = hdimids(1:2)
    alldims2d(3) = timelevels_dimid

    qdims(1:2) = hdimids(1:2)
    qdims(3) = vdimids(1)
    qdims(5) = timelevels_dimid

! Transpose yz 2D decomposition to xy 2D decomposition
    call yz2D_xy2D

    call init_restart_varlist()
        
    do i=1,restartvarcnt
    
       call get_restart_var(i, name, timelevels, ndims, vdesc)
       if(timelevels>1) then
          if(ndims==3) then
             ierr = PIO_Def_Var(File, name, pio_double, alldims2d, vdesc)
          else if(ndims==4) then
             ierr = PIO_Def_Var(File, name, pio_double, alldims, vdesc)
          else if(ndims==5) then
             if(masterproc) print*, '---- error for ndims=5 in sub. init_restart_dynamics ----------'
          end if
       else
          if(ndims==1) then
! broken i think
             ierr = PIO_Def_Var(File, name, pio_double, hdimids(2:2), vdesc)
          else if(ndims==2) then
             ierr = PIO_Def_Var(File, name, pio_double, alldims2d(1:2), vdesc)
          else if(ndims==3) then
             ierr = PIO_Def_Var(File, name, pio_double, alldims(1:3), vdesc)
          else if(ndims==4) then
             if(masterproc) print*, '---- error for ndims=4 in sub. init_restart_dynamics ----------'
          end if
       end if
    end do
    
  end subroutine init_restart_dynamics


!====================================================================================
  subroutine write_restart_dynamics (File, dyn_out)
    use cam_pio_utils, only : pio_subsystem
    use dyn_comp, only  : dyn_export_t
    use dyn_grid, only : get_horiz_grid_dim_d
    use time_manager, only: get_curr_time, get_step_size
    use prognostics,     only:  ptimelevels, n3m1, n3
    use pmgrid,       only: plon, plat, plev, beglat, endlat, beglev, endlev
!!    use ppgrid,          only: pver
    use massfix,         only: alpha, hw1, hw2, hw3
    use eul_control_mod, only: fixmas, tmass0
    use abortutils,      only: endrun

    !
    ! Input arguments
    !
    type(File_desc_t), intent(inout) :: File     ! Unit number
    type(Dyn_export_t), intent(in) :: dyn_out ! Not used in eul dycore

    !
    ! Local workspace
    !
    integer :: ierr   ! error status
    integer :: ndcur, nscur
    real(r8) :: time, dtime, mold(1)
    integer :: i, s3d(1), s2d(1), ct
    integer(kind=pio_offset) :: nt
    type(io_desc_t) :: iodesc3d, iodesc2d
    integer, pointer :: ldof(:)
    integer :: ndims, timelevels
    integer :: hdim1, hdim2
    type(var_desc_t), pointer :: vdesc
    character(len=namlen) :: name
    integer :: j, jd, k
    !
!---------------------------------------------------------------------------------    
!
!!    use_transfer = .true.
!! #if ( defined SPMD )
!!    dyn_state => get_dyn_state()
!!    if (dyn_state%grid%iam >= dyn_state%grid%npes_xy) then
!!       use_transfer = .false.
!!    end if       
!! #endif
!
    call get_curr_time(ndcur, nscur)
    dtime = get_step_size()
    call get_horiz_grid_dim_d(hdim1, hdim2)   !zhh

    ldof => get_restart_decomp(hdim1, hdim2, plev)
    call pio_initdecomp(pio_subsystem, pio_double, (/hdim1, hdim2, plev/), ldof, iodesc3d)
    deallocate(ldof)

    ldof => get_restart_decomp(hdim1, hdim2, 1)
    call pio_initdecomp(pio_subsystem, pio_double, (/hdim1, hdim2/), ldof, iodesc2d)
    deallocate(ldof)

    ierr = pio_put_var(File, tmass0desc, (/tmass0/))
    ierr = pio_put_var(File, fixmasdesc, (/fixmas/))

    ierr = pio_put_var(File, hw1desc, hw1)
    ierr = pio_put_var(File, hw2desc, hw2)
    ierr = pio_put_var(File, hw3desc, hw3)
    ierr = pio_put_var(File, alphadesc, alpha)

    ! slt and prognostic
    do nt=1,ptimelevels
       time = ndcur+(real(nscur,kind=r8)+ (nt-2)*dtime)/86400_r8
       ierr = pio_put_var(File,timedesc%varid, (/int(nt)/), time)
    end do
    do i=1,restartvarcnt
       call get_restart_var(i, name, timelevels, ndims, vdesc)

       if(timelevels==1) then
          if(ndims==2) then
             call pio_write_darray(File, vdesc, iodesc2d, transfer(restartvars(i)%v2d(:,:), mold), ierr)
          else if(ndims==3) then
             call pio_write_darray(File, vdesc, iodesc3d, transfer(restartvars(i)%v3d(:,:,:), mold), ierr)
          else if(ndims==4) then
             if(masterproc) print*, '---- error for ndims=4 in sub. write_restart_dynamics ----------'
          end if
       else
          do nt=1,timelevels
             if(nt==1) ct=n3m1
             if(nt==2) ct=n3
!!             if(nt==3) ct=n3
!
             call pio_setframe(vdesc, nt)
             if(ndims==3) then
                call pio_write_darray(File, vdesc, iodesc2d, transfer(restartvars(i)%v3d(:,:,ct), mold), ierr)
             else if(ndims==4) then
                call pio_write_darray(File, vdesc, iodesc3d, transfer(restartvars(i)%v4d(:,:,:,ct), mold), ierr)
             else if(ndims==5) then
                if(masterproc) print*, '---- error for ndims=5 in sub. write_restart_dynamics ----------'
             end if
             
          end do
          
       end if
    end do

    call pio_freedecomp(File, iodesc2d)
    call pio_freedecomp(File, iodesc3d)
    
!======================== zhh debug =======================
    do j = beglat, endlat
       jd = plat - j + 1
       do k = beglev, endlev
	      if (j==64 .and. k==30) then 
             print*, '-------- At the end of sub. write_restart_dynamics ---------'
             print*, 'j =', j, ' jd =', jd
             print*, 'n3 =', n3, ' n3m1 =', n3m1
             print*, 'phis(20,j) =', phis(20,j)
             print*, 'GHS(9,jd) =', GHS(9,jd)
             print*, 'omga(59,30,j) =', omga(59,30,j)
             print*, 'WS(26,30,jd) =', WS(26,30,jd)
             print*, 'U(96,30,jd) =', U(96,30,jd)
             print*, 'V(96,30,jd) =', V(96,30,jd)
             print*, 'ps(256,j,n3) =', ps(256,j,n3)
             print*, 'ps(256,j,n3m1) =', ps(256,j,n3m1)
             print*, 'u3(40,30,j,n3m1) =', u3(40,30,j,n3m1)
             print*, 'v3(4,26,j,n3m1) =', v3(4,26,j,n3m1)
             print*, 't3(4,26,j,n3m1) =', t3(4,26,j,n3m1) 
             print*, 'q3(4,30,2,j,n3) =', q3(4,30,2,j,n3)
             print*, 'q3(4,30,20,j,n3m1) =', q3(4,30,20,j,n3m1)
             print*, 'pdeld(9,30,j,n3) =', pdeld(9,30,j,n3) 
             print*, 'lammp(1,30,j) =', lammp(1,30,j)
             print*, 'phimp(100,30,j) =', phimp(100,30,j)
             print*, 'sigmp(240,30,j) =', sigmp(240,30,j)
             print*, 'qfcst(15,30,1,j) =', qfcst(15,30,1,j)
             print*, 'qfcst(15,26,4,j) =', qfcst(15,26,4,j)
             print*, '==============================================================='
          end if
       end do
    end do
!======================== zhh debug =======================
    
    return
  end subroutine write_restart_dynamics

!====================================================================================
  subroutine get_restart_var(i,name, timelevels, ndims, vdesc)
    integer, intent(in) :: i
    character(len=namlen), intent(out) :: name
    integer, intent(out) :: ndims, timelevels
    type(var_desc_t), pointer :: vdesc
!------------------------------------------------------------------------------------    

    name = restartvars(i)%name
    timelevels = restartvars(i)%timelevels
    ndims = restartvars(i)%ndims
    if(.not.associated(restartvars(i)%vdesc)) then
       allocate(restartvars(i)%vdesc)
    end if
    vdesc => restartvars(i)%vdesc
    call pio_setframe(vdesc, int(-1,pio_offset))

  end subroutine get_restart_var


!====================================================================================
  subroutine read_restart_dynamics (File, dyn_in, dyn_out, NLFileName)
    use dyn_comp, only : dyn_init, dyn_import_t, dyn_export_t, xy2D_yz2D  !zhh 2013-01-21
    use cam_pio_utils, only : pio_subsystem

    use pmgrid,          only: plon, plat, plev, beglat, endlat, beglonxy, endlonxy,  &
                               beglatxy, endlatxy, beglev, endlev
    use scanslt,         only: scanslt_alloc
#if ( defined BFB_CAM_SCAM_IOP )
    use iop,             only: init_iop_fields
#endif
    use massfix,         only: alpha, hw1, hw2, hw3
    use prognostics,     only: ptimelevels, n3m1, n3, initialize_prognostics
    use ppgrid,          only: pver
    use eul_control_mod, only: fixmas, tmass0
    use dynamics_vars,   only: T_FVDYCORE_STATE         !zhh 2013-01-30
    use dyn_internal_state, only : get_dyn_state   !zhh 2013-01-16

    !
    ! Input arguments
    !
    type(file_desc_t), intent(inout) :: File     ! PIO file handle
    type(dyn_import_t) :: dyn_in    ! not used by this dycore, included for compatibility
    type(dyn_export_t) :: dyn_out ! not used by this dycore, included for compatibility
    character(len=*), intent(in) :: NLFileName
    !
    ! Local workspace
    !
    type(io_desc_t) :: iodesc3d, iodesc2d
    integer, pointer :: ldof(:)
    integer :: ioerr   ! error status
    real(r8), allocatable :: tmp(:)
    type (T_FVDYCORE_STATE), pointer :: dyn_state   !zhh 2013-01-16
    !
    integer :: dims3d(3), dims2d(2)
    integer :: ierr, ct
    integer(kind=pio_offset) :: nt
    character(len=namlen) :: name
    integer :: ndims, timelevels, i, s2d, s3d
    type(var_desc_t), pointer :: vdesc
    integer :: j, jd, k
!---------------------------------------------------------------------------------------

   ! Initialize dynamics
!
    dyn_state => get_dyn_state()
    call dyn_init(dyn_state, NLFileName )    !zhh 2013-01-16

    call initialize_prognostics
    call initialize_IAPprog     !zhh 2011-12-17
	
    dims3d(1)=(endlonxy-beglonxy+1)
    dims3d(2)=(endlatxy-beglatxy+1)
    dims3d(3)=plev
    s2d = dims3d(1)*dims3d(2)
    s3d = s2d*dims3d(3)

    allocate(tmp(s3d))    

    ldof => get_restart_decomp(plon, plat, pver)
    call pio_initdecomp(pio_subsystem, pio_double, (/plon,plat,pver/), ldof, iodesc3d)
    deallocate(ldof)
    ldof => get_restart_decomp(plon, plat, 1)
    call pio_initdecomp(pio_subsystem, pio_double, (/plon,plat/), ldof, iodesc2d)
    deallocate(ldof)

    call scanslt_alloc()
    
    ierr = PIO_Inq_varid(File, 'tmass0', tmass0desc)
    ierr = pio_get_var(File, tmass0desc, tmass0)
    ierr = PIO_Inq_varid(File,'fixmas', fixmasdesc)
    ierr = pio_get_var(File, fixmasdesc, fixmas)

    ierr = PIO_Inq_varid(File, 'hw1', hw1desc)
    ierr = pio_get_var(File, hw1desc, hw1)
    ierr = PIO_Inq_varid(File, 'hw2', hw2desc)
    ierr = pio_get_var(File, hw2desc, hw2)
    ierr = PIO_Inq_varid(File, 'hw3', hw3desc)
    ierr = pio_get_var(File, hw3desc, hw3)
    ierr = PIO_Inq_varid(File,'alpha', alphadesc)
    ierr = pio_get_var(File, alphadesc, alpha)

! Initialize reatsrt variables' list 
    call init_restart_varlist()
    
#if ( defined BFB_CAM_SCAM_IOP )
    call init_iop_fields()
#endif
   
    !slt and prognostics
    do i=1,restartvarcnt
       call get_restart_var(i, name, timelevels, ndims, vdesc)

       ierr = PIO_Inq_varid(File, name, vdesc)
       if(timelevels == 1) then
          if(ndims==2) then
             call pio_read_darray(File, vdesc, iodesc2d, tmp(1:s2d), ierr)
             restartvars(i)%v2d(:,:) = reshape(tmp(1:s2d), dims3d(1:2))
          else if(ndims==3) then
             call pio_read_darray(File, restartvars(i)%vdesc, iodesc3d, tmp(1:s3d), ierr)
             restartvars(i)%v3d(:,:,:) = reshape(tmp(1:s3d), dims3d)
          else if(ndims==4) then
             if(masterproc) print*, '---- error for ndims=4 in dub. read_restart_dynamics ----------'
          end if

       else
          do nt=1,timelevels
             if(nt==1) ct=n3m1
             if(nt==2) ct=n3
!!             if(nt==3) ct=n3
             call pio_setframe(vdesc, nt)
             if(ndims==3) then
                call pio_read_darray(File, vdesc, iodesc2d, tmp(1:s2d), ierr)
                restartvars(i)%v3d(:,:,ct) = reshape(tmp(1:s2d), dims3d(1:2))
             else if(ndims==4) then
                call pio_read_darray(File, vdesc, iodesc3d, tmp(1:s3d), ierr)
                restartvars(i)%v4d(:,:,:,ct) = reshape(tmp(1:s3d), dims3d)
             else if(ndims==5) then
                if(masterproc) print*, '---- error for ndims=5 in dub. read_restart_dynamics ----------'
             end if

          end do
       end if
    end do

! Transpose xy 2D decomposition to yz 2D decomposition
    call xy2D_yz2D
    
    deallocate(tmp)
    call pio_freedecomp(File, iodesc2d)
    call pio_freedecomp(File, iodesc3d)

!======================== zhh debug =======================
    do j = beglat, endlat
       jd = plat - j + 1
       do k = beglev, endlev
	      if (j==64 .and. k==30) then 
             print*, '-------- At the end of sub. read_restart_dynamics ---------'
             print*, 'j =', j, ' jd =', jd
             print*, 'n3 =', n3, ' n3m1 =', n3m1
             print*, 'phis(20,j) =', phis(20,j)
             print*, 'GHS(9,jd) =', GHS(9,jd)
             print*, 'omga(59,30,j) =', omga(59,30,j)
             print*, 'WS(26,30,jd) =', WS(26,30,jd)
             print*, 'U(96,30,jd) =', U(96,30,jd)
             print*, 'V(96,30,jd) =', V(96,30,jd)
             print*, 'ps(256,j,n3) =', ps(256,j,n3)
             print*, 'ps(256,j,n3m1) =', ps(256,j,n3m1)
             print*, 'u3(40,30,j,n3m1) =', u3(40,30,j,n3m1)
             print*, 'v3(4,26,j,n3m1) =', v3(4,26,j,n3m1)
             print*, 't3(4,26,j,n3m1) =', t3(4,26,j,n3m1) 
             print*, 'q3(4,30,2,j,n3) =', q3(4,30,2,j,n3)
             print*, 'q3(4,30,20,j,n3m1) =', q3(4,30,20,j,n3m1)
             print*, 'pdeld(9,30,j,n3) =', pdeld(9,30,j,n3) 
             print*, 'lammp(1,30,j) =', lammp(1,30,j)
             print*, 'phimp(100,30,j) =', phimp(100,30,j)
             print*, 'sigmp(240,30,j) =', sigmp(240,30,j)
             print*, 'qfcst(15,30,1,j) =', qfcst(15,30,1,j)
             print*, 'qfcst(15,26,4,j) =', qfcst(15,26,4,j)
             print*, '==============================================================='
          end if
       end do
    end do
!======================== zhh debug =======================

    return

  end subroutine read_restart_dynamics

!======================================================================================
  function get_restart_decomp(hdim1, hdim2, nlev) result(ldof)
    use dyn_grid, only : get_dyn_grid_parm

    integer, intent(in) :: hdim1, hdim2, nlev
    integer, pointer :: ldof(:)
    integer :: i, k, j
    integer :: lcnt
    integer, allocatable :: gcols(:)

    integer :: beglatxy, beglonxy, endlatxy, endlonxy, plat
!--------------------------------------------------------------------------------------

    beglonxy = get_dyn_grid_parm('beglonxy')
    endlonxy = get_dyn_grid_parm('endlonxy')
    beglatxy = get_dyn_grid_parm('beglatxy')
    endlatxy = get_dyn_grid_parm('endlatxy')

    plat = get_dyn_grid_parm('plat')
    
    
    lcnt=(endlatxy-beglatxy+1)*nlev*(endlonxy-beglonxy+1)

    allocate(ldof(lcnt))
    lcnt=0
    ldof(:)=0	
    do j=beglatxy,endlatxy
       do k=1,nlev
          do i=beglonxy, endlonxy
             lcnt=lcnt+1
             ldof(lcnt)=i+(j-(plat-hdim2+1))*hdim1+(k-1)*hdim1*hdim2
          end do
       end do
    end do

  end function get_restart_decomp

  
end module restart_dynamics
