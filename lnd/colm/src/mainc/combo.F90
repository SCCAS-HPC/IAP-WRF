

 subroutine combo ( dz, wliq,  wice,  t, &
                   dz2, wliq2, wice2, t2 )

!=======================================================================
! Original author: Yongjiu Dai, September 15, 1999
!
! combines two elements and returns the following combined
! variabless: dz, t, wliq, wice.
! the combined temperature is based on the equation:
! the sum of the enthalpies of the two elements = that of the combined element.
!
!=======================================================================

  use precision
  use phycon_module, only : cpice, cpliq, hfus, tfrz
  implicit none

!-------------------------- Dummy argument -----------------------------

  real(r8), INTENT(in) :: dz2     ! nodal thickness of 2 elements being combined [m]
  real(r8), INTENT(in) :: wliq2   ! liquid water of element 2 [kg/m2]
  real(r8), INTENT(in) :: wice2   ! ice of element 2 [kg/m2]
  real(r8), INTENT(in) :: t2      ! nodal temperature of element 2 [K]

  real(r8), INTENT(inout) :: dz   ! nodal thickness of 1 elements being combined [m]
  real(r8), INTENT(inout) :: wliq ! liquid water of element 1
  real(r8), INTENT(inout) :: wice ! ice of element 1 [kg/m2]
  real(r8), INTENT(inout) :: t    ! nodel temperature of elment 1 [K]

!----------------------- Local variables ------------------------------

  real(r8) dzc    ! Total thickness of nodes 1 and 2 (dzc=dz+dz2).
  real(r8) wliqc  ! Combined liquid water [kg/m2]
  real(r8) wicec  ! Combined ice [kg/m2]
  real(r8) tc     ! Combined node temperature [K]
  real(r8) h      ! enthalpy of element 1 [J/m2]
  real(r8) h2     ! enthalpy of element 2 [J/m2]
  real(r8) hc     ! temporary

!-----------------------------------------------------------------------

      dzc = dz+dz2
      wicec = (wice+wice2)
      wliqc = (wliq+wliq2)
      h   = (cpice*wice+cpliq*wliq)*(t-tfrz)+hfus*wliq
      h2  = (cpice*wice2+cpliq*wliq2)*(t2-tfrz)+hfus*wliq2

      hc = h + h2
      if(hc < 0.)then
         tc = tfrz + hc/(cpice*wicec+cpliq*wliqc)
      else if(hc.le.hfus*wliqc)then
         tc = tfrz
      else
         tc = tfrz + (hc - hfus*wliqc)/(cpice*wicec+cpliq*wliqc)
      endif

      dz = dzc
      wice = wicec 
      wliq = wliqc
      t = tc

 end subroutine combo
