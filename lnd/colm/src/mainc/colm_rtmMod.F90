#include <define.h>

module colm_rtmMod

#ifdef RTM

   use spmd
   use colm_rtmVar
   implicit none

   interface colm_rtm_init
      module procedure colm_rtm_init
   end interface

   interface colm_rtm_drv
      module procedure colm_rtm_drv
   end interface

   interface colm_rtm_accum
      module procedure colm_rtm_accum
   end interface

   interface write_rtm_history
      module procedure write_rtm_history
   end interface

   interface colm_rtm_exit
      module procedure colm_rtm_exit
   end interface

CONTAINS

   subroutine colm_rtm_init

      use colm_varctl, only: rtm_dtime
      use colm_varMod, only: numgrid
      use timemgr    , only: get_step_size
      use RtmMod
      use RunoffMod
      implicit none

    ! integer num_lnd, num_ocn
      integer dtime

    ! Allocate memory to store CoLM flux

      allocate (fevpa(numgrid))
      allocate (rnof (numgrid))
      allocate (prc  (numgrid))
      allocate (prl  (numgrid))

    ! If rtm_nsteps was not entered in the namelist, 
    ! give it the following default value

      if (rtm_nsteps .eq. -9999) then
         dtime = get_step_size()
         rtm_nsteps = rtm_dtime/dtime
      end if

      if (p_master) then
         if (rtm_nsteps > 1) then
            write(6,*)'river runoff calculation performed only every ',rtm_nsteps,' nsteps'
         else
            write(6,*)'river runoff calculation performed every time step'
         endif
      endif

    ! Initialize RTM river routing grid and mask

      call Rtmgridini

    ! Initialize river routing model

      call Rtmlandini

      call Rtmfluxini

    ! Allocate memory for accum

    ! call get_proc_rof_global (num_lnd, num_ocn)

    ! allocate (lnd_ave(num_lnd))
    ! allocate (ocn_ave(num_ocn))

    ! lnd_ave(:) = 0.
    ! ocn_ave(:) = 0.

    ! lnd_nac = 0
    ! ocn_nac = 0

   end subroutine colm_rtm_init

   subroutine colm_rtm_drv

      use colm_varMod
      use RtmMod

      implicit none

      integer g

      do g = 1, numgrid
         fevpa(g) = fldv(idx_qflx,g)  ! Unit [mm/s]
          rnof(g) = fldv(idx_rnof,g)  ! Unit [mm/s]
           prc(g) = fldv(idx_prc,g)   ! Unit [mm/s]
           prl(g) = fldv(idx_prl,g)   ! Unit [mm/s]
      end do

      CALL Rtmriverflux()

   end subroutine colm_rtm_drv

   subroutine colm_rtm_accum(lnd_nac,ocn_nac,lnd_ave,ocn_ave)

      use RunoffMod
      implicit none

      integer, intent(inout) :: lnd_nac
      integer, intent(inout) :: ocn_nac
      real(r8),pointer,intent(inout) :: lnd_ave(:)
      real(r8),pointer,intent(inout) :: ocn_ave(:)

      integer beg_lnd, end_lnd, beg_ocn, end_ocn, g

      call get_proc_rof_bounds (beg_lnd, end_lnd, beg_ocn, end_ocn)

      lnd_nac = lnd_nac + 1

      do g = beg_lnd, end_lnd
         lnd_ave(g) = lnd_ave(g) + runoff%lnd(g)
      end do

      ocn_nac = ocn_nac + 1

      do g = beg_ocn, end_ocn
         ocn_ave(g) = ocn_ave(g) + runoff%ocn(g)
      end do

   end subroutine colm_rtm_accum

   subroutine write_rtm_history(lnd_nac,ocn_nac,lnd_ave,ocn_ave,flxmask)

      use RunoffMod
#ifdef SPMD
      use spmdGathScatMod, only : gather_data_to_master
#endif
      use nchistMod, only: nchist_addfield

      implicit none

      logical, intent(in)    :: flxmask(:)
      integer, intent(inout) :: lnd_nac
      integer, intent(inout) :: ocn_nac
      real(r8), pointer, intent(inout) :: lnd_ave(:)
      real(r8), pointer, intent(inout) :: ocn_ave(:)

      real(r8), pointer :: lnd_globdc(:)
      real(r8), pointer :: ocn_globdc(:)
      real(r4), pointer :: lnd_globxy(:,:)  ! RTM river flow [m3/s]
      real(r4), pointer :: ocn_globxy(:,:)  ! RTM river discharge into ocean [m3/s]

      integer num_lnd, num_ocn
      integer i, j, g

      call get_proc_rof_global (num_lnd, num_ocn)

      allocate (lnd_globdc(num_lnd))
      allocate (ocn_globdc(num_ocn))
      allocate (lnd_globxy(rtmlon,rtmlat))
      allocate (ocn_globxy(rtmlon,rtmlat))

#ifdef SPMD
      call gather_data_to_master (lnd_ave, lnd_globdc, clmlevel='lndrof')
      call gather_data_to_master (ocn_ave, ocn_globdc, clmlevel='ocnrof')
#else
      lnd_globdc(:) = lnd_ave(:)
      ocn_globdc(:) = ocn_ave(:)
#endif

      lnd_globdc(:) = lnd_globdc(:)/lnd_nac
      ocn_globdc(:) = ocn_globdc(:)/ocn_nac

      lnd_globxy(:,:) = -9999.
      ocn_globxy(:,:) = -9999.

      do g = 1,num_lnd
         i = runoff%lnd_ixy(g)
         j = runoff%lnd_jxy(g)
         lnd_globxy(i,j) = lnd_globdc(g)
      end do

      do g = 1,num_ocn
         i = runoff%ocn_ixy(g)
         j = runoff%ocn_jxy(g)
         ocn_globxy(i,j) = ocn_globdc(g)
      end do

    ! Write history file

      if (p_master) then
         if(flxmask(1)) call nchist_addfield(1,rtmlon,rtmlat,lnd_globxy,'rtm')
         if(flxmask(2)) call nchist_addfield(2,rtmlon,rtmlat,ocn_globxy,'rtm')
      end if

    ! Reset accum variables

      lnd_nac    = 0
      ocn_nac    = 0
      lnd_ave(:) = 0.
      ocn_ave(:) = 0.

      deallocate (lnd_globdc)
      deallocate (ocn_globdc)
      deallocate (lnd_globxy)
      deallocate (ocn_globxy)

   end subroutine write_rtm_history

   subroutine colm_rtm_exit

      deallocate (fevpa)
      deallocate ( rnof)
      deallocate (  prc)
      deallocate (  prl)

    ! deallocate (lnd_ave)
    ! deallocate (ocn_ave)

   end subroutine colm_rtm_exit

#endif

end module colm_rtmMod
