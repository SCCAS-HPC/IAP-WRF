#include <define.h>

  subroutine DGVMRespiration(dtime, ivt, fpcgrid, nind, dphen, lm_ind, sm_ind,&
                             rm_ind, respcoeff, l_cton, s_cton, r_cton, tl,&
                             tsoi25, assim, respc, frmf, dnpp, &
                             dgpp, dfrmf, dfrms, dfrmr, dfrg)  
                             
!--------------------------------------------------------------
! Calculates surface biogeochemical fluxes
! CALLED FROM: CLMMAIN
!-------------------------------------------------------------
    use precision
    implicit none

! INTENT IN VARIABLES:
    integer , INTENT(in) :: ivt         ! vegetation type for current pft
    real(r8), INTENT(in) :: dtime       ! model time step [second]
    real(r8), INTENT(in) :: fpcgrid     ! foliar projective cover on gridcell (fraction)
    real(r8), INTENT(in) :: nind        ! number of individuals
    real(r8), INTENT(in) :: dphen       ! phenology [0 to 1]
    real(r8), INTENT(in) :: lm_ind      ! individual leaf mass
    real(r8), INTENT(in) :: sm_ind      ! individual stem mass
    real(r8), INTENT(in) :: rm_ind      ! individual root mass
    real(r8), INTENT(in) :: respcoeff   ! respiration coefficient (LPJ)
    real(r8), INTENT(in) :: l_cton      ! c/n for leaves (LPJ)
    real(r8), INTENT(in) :: s_cton      ! c/n for stems (LPJ)
    real(r8), INTENT(in) :: r_cton      ! c/n for roots (LPJ)
    real(r8), INTENT(in) :: tl          ! vegetation temperature (Kelvin)
    real(r8), INTENT(in) :: tsoi25      ! soil temperature to 0.25 m (Kelvin)
    real(r8), INTENT(in) :: assim       ! photosynthesis (mol CO2 /m**2 /s)

! INTENT INOUT VARIABLES:

    real(r8), INTENT(inout) :: respc    ! autotropic respiration (mol CO2 /m**2/s)

! INTENT OUT VARIABLES:
    real(r8), INTENT(out) :: frmf        ! leaf maintenance respiration  (mol CO2 /m**2 /s)
    real(r8), INTENT(out) :: dnpp        ! total dry matter production (gC /m**2 vegt'd area/tstep)
    real(r8), INTENT(out) :: dgpp
    real(r8), INTENT(out) :: dfrmf
    real(r8), INTENT(out) :: dfrms
    real(r8), INTENT(out) :: dfrmr
    real(r8), INTENT(out) :: dfrg

! OTHER LOCAL VARIABLES:
!
    real(r8), parameter :: k = 0.0548 / 86400.0 ! from [/day] to [/second]
    real(r8), parameter :: dmcf = 24.0          ! co2-to-biomass conversion (g biomass/mol CO2)
!   real(r8) :: frmf       ! leaf maintenance respiration  (mol CO2 /m**2 /s)
    real(r8) :: frms       ! stem maintenance respiration  (mol CO2 /m**2 /s)
    real(r8) :: frmr       ! root maintenance respiration  (mol CO2 /m**2 /s)
    real(r8) :: frm        ! total maintenance respiration (mol CO2 /m**2/s)
    real(r8) :: frg        ! growth respiration (mol CO2 /m**2 /s)
    real(r8) :: tf1        ! temperature factor
    real(r8) :: tf2        ! temperature factor
!   real(r8) :: fco2       ! net CO2 flux (umol CO2 /m**2 /s) [+ = to atm]


    ! Set co2-to-biomass conversion (ug biomass / umol CO2)

    ! maintenance respiration: LPJ equations w/ units of [gC m-2 gridcell s-1]
    ! converted to units of [umol CO2 m-2 patch s-1]

    if (ivt .le. 16 .and. fpcgrid > 0.0) then

       if (tl >= 273.16-40.) then
          tf1 = exp(308.56 * (1.0/56.02 - 1.0/(tl-227.13)))
       else
          tf1 = 0.0
       end if
       if (tsoi25 >= 273.16-40.) then
          tf2 = exp(308.56 * (1.0/56.02 - 1.0/(tsoi25-227.13)))
       else
          tf2 = 0.0
       end if

#if (defined DyN)
!      frmf = respc 
       frmf = respcoeff * k * lm_ind * nind * l_cton &
            * tf1 * dphen * 2.0 / dmcf / fpcgrid
       frms = respcoeff * k * sm_ind * nind * s_cton &
            * tf1 * 2.0 / dmcf / fpcgrid
       frmr = respcoeff * k * rm_ind * nind * r_cton &
            * tf2 * dphen * 2.0 / dmcf / fpcgrid
#else
       frmf = respcoeff * k * lm_ind * nind / l_cton &
            * tf1 * dphen * 2.0 / dmcf / fpcgrid
       frms = respcoeff * k * sm_ind * nind / s_cton &
            * tf1 * 2.0 / dmcf / fpcgrid
       frmr = respcoeff * k * rm_ind * nind / r_cton &
            * tf2 * dphen * 2.0 / dmcf / fpcgrid
#endif
!print *,'l',l_cton,'s',s_cton,'r',r_cton
!print *,ivt, 'lm_ind', lm_ind, 'frmf', frmf
!print *,ivt, 'sm_ind', sm_ind, 'frms', frms
!print *,ivt, 'rm_ind', rm_ind, 'frmr',r_cton, frmr,dphen,fpcgrid,respcoeff,nind
    else
       frmf = 0.0
       frms = 0.0
       frmr = 0.0
    end if

    frm  = frmf + frms + frmr

    ! growth respiration and production

    frg = 0.25 * max(assim - frm, 0.0)      !changed to match LPJ
    respc = frm + frg
    ! npp gC/m2 vegt'd area /step
    dnpp = (assim - respc) * dmcf * 0.5 * dtime * fpcgrid

    dgpp  = assim * dmcf * 0.5 * dtime
    dfrmf = frmf  * dmcf * 0.5 * dtime
    dfrms = frms  * dmcf * 0.5 * dtime
    dfrmr = frmr  * dmcf * 0.5 * dtime
    dfrg  = frg   * dmcf * 0.5 * dtime

    ! bm_inc=[gC/m2 patch area] from dmi=[ug dry matter/m2 patch area/s]

!    bm_inc = bm_inc + dmi * dtime * 0.5

    ! microbial respiration

    ! DGVM calculates clm%fmicr in LitterSOM; problem with units in relation
    ! to history grid averages; {fmicr}=[gC/m2 gridcell vegetated area] in
    ! LPJ calculation => sum(fmicr) over a gridcell would give the total for
    ! the gridcell; it would be wrong to convert to [/m2 patch area] because
    ! soil carbon continues to exist even when the plants above it die and
    ! fpcgrid goes to 0; how to reconcile with history calculation which
    ! will take the following values and weight them by the corresponding
    ! weights of their patches? Could chg. soil carbon and plant litter to
    ! gridcell level pools; this will affect fco2, as well; could treat this
    ! issue simultaneously with soil water, which we also want converted to
    ! the gridcell level for plant competition. (slevis)

!    afmicr = afmicr + fmicr ![gC/m2 gridcell vegetated area]

    ! net CO2 flux

!    fco2 = -assim + respc + fmicr

 end subroutine DGVMRespiration

