module Dyn_const
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plev

   implicit none

   save
   public

   real(r8) :: PMTOP      ! pressure at the top level of the model (unit: hPa)
   real(r8) :: SIG(plev)    ! sigma value on model layer
   real(r8) :: SIGL(plev+1)   ! sigma value on interface layer

end module Dyn_const

