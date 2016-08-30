module mathconst
!------------------------------------------------------------------------------------------------
! Purpose:    Set math constants used usually 
! Author :    ZhangHe
! Completed : 2005.8.26
!------------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none
   save
   public

   real(r8), parameter :: ZERO   = 0.0E0
   real(r8), parameter :: ONE    = 1.0E0
   real(r8), parameter :: TWO    = 2.0E0
   real(r8), parameter :: THREE  = 3.0E0
   real(r8), parameter :: FOUR   = 4.0E0
   real(r8), parameter :: FIVE   = 5.0E0
   real(r8), parameter :: SIX    = 6.0E0
   real(r8), parameter :: HALF   = 0.5E0
   real(r8), parameter :: FOURTH = 0.25E0
   real(r8), parameter :: PI     = 3.14159265358979323846
  
end module