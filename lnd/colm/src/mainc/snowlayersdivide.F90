
 subroutine snowlayersdivide (lb,snl,z,dz,zi,wliq,wice,tss)

!=======================================================================
! Original author : Yongjiu Dai, September 15, 1999
!
! subdivides snow layer when its thickness exceed the prescribed maximum
!=======================================================================

  use precision
  implicit none

!-------------------------- Dummy argument -----------------------------

   integer, INTENT(in) :: lb              ! lower bound of array
   integer, INTENT(inout) :: snl          ! Number of snow
  real(r8), INTENT(inout) :: wice(lb:0)   ! ice lens [kg/m2]
  real(r8), INTENT(inout) :: wliq(lb:0)   ! liquid water [kg/m2]
  real(r8), INTENT(inout) :: tss (lb:0)   ! Nodel temperature [K]
  real(r8), INTENT(inout) :: dz  (lb:0)   ! Layer thickness [m]
  real(r8), INTENT(inout) :: z   (lb:0)   ! Node depth [m]
  real(r8), INTENT(inout) :: zi  (lb-1:0) ! Depth of layer interface [m]

!----------------------- Local variables ------------------------------

! numbering from 1 (surface) msno (bottom)
  real(r8) :: drr      ! thickness of the combined [m]
  real(r8) :: dzsno(5) ! Snow layer thickness [m] 
  real(r8) :: swice(5) ! Partial volume of ice [m3/m3]
  real(r8) :: swliq(5) ! Partial volume of liquid water [m3/m3]
  real(r8) :: tsno(5)  ! Nodel temperature [K]

  integer k            ! number of do looping
  integer msno         ! number of snow layer 1 (top) to msno (bottom)

  real(r8) zwice,zwliq,propor

!-----------------------------------------------------------------------

      msno = abs(snl)
      do k = 1, msno
         dzsno(k) = dz  (k + snl)
         swice(k) = wice(k + snl)
         swliq(k) = wliq(k + snl)
         tsno(k)  = tss (k + snl)
      enddo

      if(msno == 1)then
         if(dzsno(1) > 0.03)then
         msno = 2
! Specified a new snow layer
         dzsno(1) = dzsno(1)/2.
         swice(1) = swice(1)/2.
         swliq(1) = swliq(1)/2.

         dzsno(2) = dzsno(1)
         swice(2) = swice(1)
         swliq(2) = swliq(1)
         tsno(2)  = tsno(1)
!        write(6,*)'Subdivided Top Node into two layer (1/2)'
         endif
      endif

      if(msno > 1)then
         if(dzsno(1) > 0.02)then
         drr = dzsno(1) - 0.02
         propor = drr/dzsno(1)
         zwice = propor*swice(1)
         zwliq = propor*swliq(1)

         propor = 0.02/dzsno(1)
         swice(1) = propor*swice(1)
         swliq(1) = propor*swliq(1)
         dzsno(1) = 0.02

         call combo(dzsno(2),swliq(2),swice(2),tsno(2), &
                    drr,zwliq,zwice,tsno(1))

!        write(6,*) 'Subdivided Top Node &
!                    20 mm combined into underlying neighbor'

         if(msno <= 2 .AND. dzsno(2) > 0.07)then
! subdivided a new layer
            msno = 3
            dzsno(2) = dzsno(2)/2.
            swice(2) = swice(2)/2.
            swliq(2) = swliq(2)/2.

            dzsno(3) = dzsno(2)
            swice(3) = swice(2)
            swliq(3) = swliq(2)
            tsno(3)  = tsno(2)
         endif
         endif
      endif

      if(msno > 2)then
         if(dzsno(2) > 0.05)then
         drr = dzsno(2) - 0.05
         propor = drr/dzsno(2)
         zwice = propor*swice(2)
         zwliq = propor*swliq(2)

         propor = 0.05/dzsno(2)
         swice(2) = propor*swice(2)
         swliq(2) = propor*swliq(2)
         dzsno(2) = 0.05

         call combo(dzsno(3),swliq(3),swice(3),tsno(3), &
                    drr,     zwliq,   zwice,   tsno(2))

!        write(6,*)'Subdivided 50 mm from the subsface layer &
!                   &and combined into underlying neighbor'

         if(msno <= 3 .AND. dzsno(3) > 0.18)then
! subdivided a new layer
            msno =  4
            dzsno(3) = dzsno(3)/2.
            swice(3) = swice(3)/2.
            swliq(3) = swliq(3)/2.

            dzsno(4) = dzsno(3)
            swice(4) = swice(3)
            swliq(4) = swliq(3)
            tsno(4)  = tsno(3)
         endif
         endif
      endif

      if(msno > 3)then
         if(dzsno(3) > 0.11)then
         drr = dzsno(3) - 0.11
         propor = drr/dzsno(3)
         zwice = propor*swice(3)
         zwliq = propor*swliq(3)

         propor = 0.11/dzsno(3)
         swice(3) = propor*swice(3)
         swliq(3) = propor*swliq(3)
         dzsno(3) = 0.11

         call combo(dzsno(4),swliq(4),swice(4),tsno(4), &
                    drr,     zwliq,   zwice,   tsno(3))

!        write(6,*)'Subdivided 110 mm from the third Node &
!                   &and combined into underlying neighbor'

         if(msno <= 4 .AND. dzsno(4) > 0.41)then
! subdivided a new layer
            msno = 5
            dzsno(4) = dzsno(4)/2.
            swice(4) = swice(4)/2.
            swliq(4) = swliq(4)/2.

            dzsno(5) = dzsno(4)
            swice(5) = swice(4)
            swliq(5) = swliq(4)
            tsno(5)  = tsno(4)
         endif
         endif
      endif
         
      if(msno > 4)then
         if(dzsno(4) > 0.23)then
         drr = dzsno(4) - 0.23
         propor = drr/dzsno(4)
         zwice = propor*swice(4)
         zwliq = propor*swliq(4)

         propor = 0.23/dzsno(4)
         swice(4) = propor*swice(4)
         swliq(4) = propor*swliq(4)
         dzsno(4) = 0.23

         call combo(dzsno(5),swliq(5),swice(5),tsno(5), &
                    drr,     zwliq,   zwice,   tsno(4))

!        write(6,*)'Subdivided 230 mm from the fourth Node &
!                   'and combined into underlying neighbor'
         endif
      endif
         
      snl = - msno

      do k = snl+1, 0
         dz(k)   = dzsno(k - snl)
         wice(k) = swice(k - snl)
         wliq(k) = swliq(k - snl)
         tss(k)  = tsno (k - snl)
      enddo

      do k = 0, snl+1, -1
         z(k)    = zi(k) - 0.5*dz(k)
         zi(k-1) = zi(k) - dz(k)
      enddo

 end subroutine snowlayersdivide
