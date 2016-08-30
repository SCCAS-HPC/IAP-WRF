module SubgridMod
!---------------------------------------------------------------------------------
!Description:
!Utilities to perfrom subgrid averaging
!All subgrid structure utilities are contained in the module  07/07/2009/ zhq
!---------------------------------------------------------------------------------

use precision
implicit none

save


contains

  subroutine p2c(npft,num_filterp,filterp,wtpft,wtcol,varpft,varcol)

!intent in variables
     integer, intent(in) :: npft
     integer, intent(in) :: num_filterp
     integer, intent(in) :: filterp(npft)
     real(r8),intent(in) :: wtpft(npft)
     real(r8),intent(in) :: wtcol
     real(r8),intent(in) :: varpft(npft)

!intent out variables
     real(r8),intent(out) :: varcol

!local variables
     integer p,fp

     varcol = 0.
     do fp = 1, num_filterp
        p = filterp(fp)
        varcol = varcol + varpft(p) * wtpft(p) / wtcol
     enddo
   
  end subroutine p2c


  subroutine BuildFilter(npft, itypwat, wt_patch, num_filterp, filterp)

!------------------------------------------------------------------------
! Filter patches caculated in CLM processes.
! CALLED FROM: subroutine CLMMAIN
!------------------------------------------------------------------------

    integer , INTENT(in) :: npft                     ! total patches
    integer , INTENT(in) :: itypwat                  ! pft vegetation (pft level)
    real(r8), INTENT(in) :: wt_patch(npft)           ! whether this pft present in patch
    integer , INTENT(out):: num_filterp              ! number of pfts in naturally-vegetated filter
    integer , INTENT(out):: filterp(npft)            ! pft filter for naturally-vegetated points
!
! LOCAL VARIABLES:
    integer :: p
!-----------------------------------------------------------------------

    num_filterp = 0
    do p = 1,npft
       if (wt_patch(p) > 0._r8 .and. wt_patch(p) <= 1._r8+1.e-6_r8) then
           num_filterp = num_filterp + 1
           filterp(num_filterp) = p
       elseif (wt_patch(p) < 0._r8 .or. wt_patch(p) > 1._r8+1.e-6_r8) then
           stop 'can not filter inreasonable patch weight'
       end if
    end do

    if (itypwat == 0) then
       if(num_filterp < 0 .or. num_filterp > 17 ) then
           stop 'filtered incredible pft patch number'
       endif 
    else 
       if(num_filterp < 0 .or. num_filterp > 1 ) then
           stop 'filtered incredible col patch number'
       endif
    endif

  end subroutine BuildFilter

end module SubgridMod
