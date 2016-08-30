
  subroutine diag_dynvar_ic(phis, ps, t3, u3, v3, q3)
!
!----------------------------------------------------------------------- 
! 
! Purpose: record state variables to IC file
!
! Modified: Jiang Jinrong, October 2012, for 2D parellel 
!           Zhang He, 2013-03-21, written as 2D 
!           Zhang He, 2013-03-25, revised an error
!-----------------------------------------------------------------------
!
    use shr_kind_mod, only: r8 => shr_kind_r8
    use pmgrid, only: beglonxy, endlonxy, beglatxy, endlatxy, plon, plat, plev
    use cam_history , only: outfld, write_inithist
    use constituents, only: pcnst, cnst_name
    use commap, only:clat,clon
    use dyn_grid,     only : get_horiz_grid_d
    implicit none
!
!-----------------------------------------------------------------------
!
! Arguments
!
    real(r8), intent(in) :: phis(beglonxy:endlonxy, beglatxy:endlatxy) ! Surface geopotential
    real(r8), intent(in) :: ps  (beglonxy:endlonxy, beglatxy:endlatxy) ! surface pressure
    real(r8), intent(in) :: t3  (beglonxy:endlonxy, plev, beglatxy:endlatxy) ! temperature
    real(r8), intent(in) :: u3  (beglonxy:endlonxy, plev, beglatxy:endlatxy) ! u-wind component
    real(r8), intent(in) :: v3  (beglonxy:endlonxy, plev, beglatxy:endlatxy) ! v-wind component
    real(r8), intent(in) :: q3  (beglonxy:endlonxy, plev, pcnst, beglatxy:endlatxy) ! constituents
!
! Local workspace
    real(r8) :: clat_plon(plon) ! constituents
    real(r8) :: phi(plat) ! constituents
    real(r8) :: lam(plon) ! constituents
!
!---------------------------Local workspace-----------------------------
    integer :: i, j, k, m, idim   ! indices
    real(r8):: tmp(beglonxy:endlonxy, plev)
!
!-----------------------------------------------------------------------
    idim = endlonxy - beglonxy + 1
!
    if( write_inithist() ) then

!$OMP PARALLEL DO PRIVATE (LAT, M)
       do j = beglatxy, endlatxy
!
          call outfld('PS&IC      ' , ps  (:,j), idim, j)
!
          do k = 1, plev
             do i = beglonxy, endlonxy
                tmp(i,k) = t3(i,k,j)
             enddo
          enddo
          call outfld('T&IC       ' , tmp, idim, j)
!
          do k = 1, plev
             do i = beglonxy, endlonxy
                tmp(i,k) = u3(i,k,j)
             enddo
          enddo
          call outfld('U&IC       ' , tmp, idim, j)
!
          do k = 1, plev
             do i = beglonxy, endlonxy
                tmp(i,k) = v3(i,k,j)
             enddo
          enddo
          call outfld('V&IC       ' , tmp, idim, j)
!
          do m = 1, pcnst
		     do k = 1, plev
                do i = beglonxy, endlonxy
                   tmp(i,k) = q3(i,k,m,j)
                end do
	         enddo
             call outfld(trim(cnst_name(m))//'&IC', tmp, idim, j)
          enddo
!
#if (defined BFB_CAM_SCAM_IOP) 
          clat_plon(:)=clat(lat)
          call outfld('CLAT1&IC    ', clat_plon,      plon, lat)
          call outfld('CLON1&IC    ', clon,      plon, lat)
          call get_horiz_grid_d(plat, clat_d_out=phi)
          call get_horiz_grid_d(plon, clon_d_out=lam)
          clat_plon(:)=phi(lat)
          call outfld('LAM&IC    ', lam,      plon, lat)
          call outfld('PHI&IC    ', clat_plon,      plon, lat)
#endif

       end do

    end if

    return
  end subroutine diag_dynvar_ic
