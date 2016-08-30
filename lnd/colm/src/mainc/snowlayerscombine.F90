
 subroutine snowlayerscombine ( lb, snl, &
                     z, dz, zi, wliq, wice, tss, scv, snowdp )

!=======================================================================
! Original author : Yongjiu Dai, September 15, 1999
!
! checks for elements which are below prescribed minimum for thickness or mass.
! If snow element thickness or mass is less than a prescribed minimum,
! it is combined with neighboring element to be best combine with,
! and executes the combination of mass and energy in clm_combo.f90
! 
!=======================================================================

  use precision
  implicit none

!-------------------------- Dummy argument -----------------------------
  integer, INTENT(in) :: lb               ! lower bound of array

! numbering from 1 (bottom) mss (surface)
  real(r8), INTENT(inout) :: wice(lb:1)   ! ice lens [kg/m2]
  real(r8), INTENT(inout) :: wliq(lb:1)   ! liquid water {kg/m2]
  real(r8), INTENT(inout) :: tss (lb:1)   ! nodel temperature [K]
  real(r8), INTENT(inout) :: dz  (lb:1)   ! layer thickness [m]
  real(r8), INTENT(inout) :: z   (lb:1)   ! node depth [m]
  real(r8), INTENT(inout) :: zi  (lb-1:1) ! depth of layer interface [m]
  real(r8), INTENT(inout) :: snowdp       ! snow depth [m]
  real(r8), INTENT(inout) :: scv          ! snow mass - water equivalent [kg/m2]
  integer, INTENT(inout) :: snl           ! Number of snow

!----------------------- Local variables ------------------------------
  real(r8) :: drr           ! thickness of the combined [m]
  real(r8) :: dzmin(5)      ! minimum of snow layer 1 (top) to msn0 (bottom)
  real(r8) :: zwice         ! total ice mass in snow
  real(r8) :: zwliq         ! total liquid water in snow

  integer :: i              ! number of do looping
  integer :: j              ! node index
  integer :: k              ! number of do looping
  integer :: l              ! node index
  integer :: msn_old        ! number of snow layer 1 (top) to msn0 (bottom)
  integer :: mssi           ! node index
  integer :: neibor         ! adjacent node selected for combination

  data dzmin /0.010, 0.015, 0.025, 0.055, 0.115/

!-----------------------------------------------------------------------
! check the mass of ice lens of snow, when the total less than a small value,
! combine it with the underlying neighbor
      msn_old = snl
      do j = msn_old+1, 0
         if(wice(j) <= .1)then
            wliq(j+1) = wliq(j+1) + wliq(j)
            wice(j+1) = wice(j+1) + wice(j)

! shift all elements above this down one.
            if(j > snl+1 .AND. snl < -1)then
               do i =  j, snl+2, -1
                  tss(i) = tss(i-1)
                  wliq(i) = wliq(i-1)
                  wice(i) = wice(i-1)
                  dz(i) = dz(i-1)
               enddo
            endif

            snl = snl + 1
!*          write(6,*) 'one snow layer is gone'

         endif

      enddo

      if(snl == 0)then
         scv = 0.
         snowdp = 0.
!*       write(6,*) 'all snow has gone'
         return
      else
         scv = 0.
         snowdp = 0.
         zwice = 0.
         zwliq = 0.
         do j = snl + 1, 0
            scv = scv + wice(j) + wliq(j)
            snowdp = snowdp + dz(j)
            zwice = zwice + wice(j)
            zwliq = zwliq + wliq(j)
         enddo
      endif
!-----------------------------------------------------------------------
! check the snow depth

      if(snowdp < 0.01)then       !!! all snow gone 
          
         snl = 0
         scv = zwice
         if(scv <= 0.) snowdp = 0.

! the liquid water assumed ponding on soil surface
         wliq(1) = wliq(1) + zwliq
!*       write(6,'(17h all snow is gone)')
         return

      else                        !!! snow layers combined

! two or more layers 

        if(snl < -1)then
           msn_old = snl
           mssi = 1
           do i = msn_old+1, 0

! If top node is removed, combine with bottom neighbor
              if(dz(i) < dzmin(mssi))then
                 if(i == snl+1)then
                    neibor = i + 1

! If the bottom neighbor is not snow, combine with the top neighbor
                 else if(i == 0)then
                    neibor = i - 1

! If none of the above special cases apply, combine with the thinnest neighbor
                 else
                   neibor = i + 1
                   if((dz(i-1)+dz(i)) < (dz(i+1)+dz(i))) neibor = i-1
                 endif

! Node l and j are combined and stored as node j.

                 if(neibor > i)then
                    j = neibor
                    l = i
                 else
                    j = i
                    l = neibor
                 endif
                 call combo ( dz(j), wliq(j), wice(j), tss(j),&
                              dz(l), wliq(l), wice(l), tss(l) )

! Now shift all elements above this down one.

                 if(j-1 > snl+1) then
                    do k = j-1, snl+2, -1
                       tss(k) = tss(k-1)
                       wice(k) = wice(k-1)
                       wliq(k) = wliq(k-1)
                       dz(k) = dz(k-1)
                    enddo
                 endif

                 snl = snl + 1

!*    write(6,'(7h Nodes ,i4,4h and,i4,14h combined into,i4)') l,j,j

                 if(snl >= -1) EXIT

! The layer thickness great than the prescibed minimum value

              else
                 mssi = mssi + 1 
              endif
           enddo

        end if

! Reset the node depth and the depth of layer interface

        do k = 0, snl+1, -1
           z(k) = zi(k) - 0.5*dz(k)
           zi(k-1) = zi(k) - dz(k)
        enddo

      endif                       !!! snow layers combined 

 end subroutine snowlayerscombine
