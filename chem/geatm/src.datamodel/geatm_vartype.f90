module geatm_vartype
  parameter(ndomain=5) ! added by chenhs
  integer, dimension(ndomain) :: imother !,ichild  ! added by chenhs
  integer, dimension(ndomain) :: nx, ny ,nz, ntbeg,ntend,dtout,dtmet  !added by chenhs
  integer, dimension(ndomain) :: nxlo, nylo

  integer, dimension(ndomain) :: comm2d, nbrleft, nbrright, nbrtop, nbrbottom, &
                                 sx, ex, sy, ey, stride
  integer, dimension(2,ndomain) :: dims, coords
  logical, dimension(2,ndomain) :: periods
  integer,dimension(:,:), allocatable :: procs ! for cam, added by Juanxiong He
  integer, dimension(ndomain) :: csx, cex, csy, cey, cnx, cny ! for cam, added by Juanxiong He
  integer myid, newid, numprocs, myid_x, myid_y, local_com ! for cam, added by Juanxiong He
  integer iyear1,imonth1,idate1,ihour1
contains
   
   subroutine final_geatm_var
   end subroutine final_geatm_var
end module geatm_vartype
 
