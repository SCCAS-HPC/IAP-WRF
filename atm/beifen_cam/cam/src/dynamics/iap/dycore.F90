module dycore
!
! Data and utility routines related to the dycore
!
! Modified: Zhang He, Oct 2012

   implicit none

PRIVATE

   public :: dycore_is, get_resolution

CONTAINS

   logical function dycore_is (name)
!
! Input arguments
!
      character(len=*), intent(in) :: name
      
      if (name == 'IAP' .or. name== 'iap') then   !zhh 2011-11-16
         dycore_is = .true.
      else
         dycore_is = .false.
      end if
      
      return
   end function dycore_is

   character(len=7) function get_resolution()

     use pmgrid, only: plat

     select case ( plat )
     case ( 128 )
        get_resolution = '1.4x1.4'
     case ( 181 )
        get_resolution = '1x1'
     case ( 361 )
        get_resolution = '0.5x0.5'
     case default
        get_resolution = 'UNKNOWN'
     end select

     return
   end function get_resolution

end module dycore


