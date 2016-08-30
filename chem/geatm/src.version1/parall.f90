module parall
contains
 subroutine exchng2( myid, v, sx, ex, sy, ey,  &
                      comm2d, stride,  &
                      nbrleft, nbrright, nbrtop, nbrbottom  )
  include "mpif.h"
  integer myid, sx, ex, sy, ey, stride
  real v(sx-1:ex+1,sy-1:ey+1)
  integer nbrleft, nbrright, nbrtop, nbrbottom, comm2d
  integer status(mpi_status_size), ierr, nx
!
  nx = ex - sx + 1
!  these are just like the 1-d versions, except for less data
  call mpi_sendrecv( v(sx,ey),  nx, mpi_real,  &
                    nbrtop, 0,  &
                    v(sx,sy-1), nx, mpi_real,  &
                    nbrbottom, 0, comm2d, status, ierr )
  call mpi_sendrecv( v(sx,sy),  nx, mpi_real,  &
                    nbrbottom, 1,  &
                    v(sx,ey+1), nx, mpi_real,  &
                    nbrtop, 1, comm2d, status, ierr )
! this uses the "strided" datatype
!       v(ex,sy-1) = -100 - myid
  call mpi_sendrecv( v(ex,sy-1),  1, stride, nbrright,  2,  &
                     v(sx-1,sy-1), 1, stride, nbrleft,  2,  &
                     comm2d, status, ierr )
!       v(sx,sy-1) = -200 - myid
  call mpi_sendrecv( v(sx,sy-1),  1, stride, nbrleft,   3,  &
                     v(ex+1,sy-1), 1, stride, nbrright, 3,  &
                     comm2d, status, ierr )
  return
  end subroutine
!---------------------------------------------------------------------------------------
 subroutine mpe_decomp1d( n, numprocs, myid, s, e )
  integer n, numprocs, myid, s, e
  integer nlocal
  integer deficit
!
  nlocal  = n / numprocs
  s       = myid * nlocal + 1
  deficit = mod(n,numprocs)
  s       = s + min(myid,deficit)
  if (myid .lt. deficit) then
      nlocal = nlocal + 1
  endif
  e = s + nlocal - 1
  if (e .gt. n .or. myid .eq. numprocs-1) e = n
  return
  end subroutine
!---------------------------------------------------------------------------------------
  subroutine exinfo(myid,numprocs,aa1,aa2,aa3,aa4,&
                   smark,stadom)
  use geatm_vartype, only : local_com
  include "mpif.h"
  integer myid,numprocs,stadom,istd
  integer aa1, aa2, aa4, aa3
  integer status(mpi_status_size),ierr

  integer :: smark(stadom)

  do istd=1,4*numprocs  
    smark(istd)=0
  enddo

  smark(4*myid+1)=aa1
  smark(4*myid+2)=aa2
  smark(4*myid+3)=aa3
  smark(4*myid+4)=aa4

  idstd=4*myid+1

  if(myid>0)then
  do idetc=0,(myid-1)
  call mpi_send(smark(idstd),4,mpi_integer,idetc,0,local_com,ierr)
  enddo
  endif

  if(myid<(numprocs-1))then
  do idetc=(myid+1),(numprocs-1)
  call mpi_send(smark(idstd),4,mpi_integer,idetc,0,local_com,ierr)
  enddo
  endif

  if(numprocs .gt. 1)call mpi_barrier( local_com, ierr )

  if(myid>0)then
  do idetc=0,(myid-1)
  idstd=idetc*4+1
  call mpi_recv(smark(idstd),4,mpi_integer,idetc,0,local_com,status,ierr)
  enddo
  endif
  if(myid<(numprocs-1))then
  do idetc=(myid+1),(numprocs-1)
  idstd=idetc*4+1
  call mpi_recv(smark(idstd),4,mpi_integer,idetc,0,local_com,status,ierr)
  enddo
  endif

  if(numprocs .gt. 1)call mpi_barrier( local_com, ierr )

  if(aa1<0)then
    isid=smark(myid*4+2)
    do
      if(smark(isid*4+2)<0)then
        aa1=isid
        goto 1200
      else
        isid=smark(isid*4+2)
      endif
    enddo
  endif
  1200 continue

  if(aa2<0)then
    isid=smark(myid*4+1)
    do
      if(smark(isid*4+1)<0)then
        aa2=isid
        goto 1201
      else
        isid=smark(isid*4+1)
      endif
    enddo
  endif
  1201 continue

  return
 end subroutine exinfo
!---------------------------------------------------------------------------------------
 subroutine gather_whole_field(patch_array, ps1, pe1, ps2, pe2, ps3, pe3, &
                                domain_array, ds1, de1, ds2, de2, ds3, de3)
      use geatm_vartype, only : procs, local_com, myid, dims
      implicit none
      include "mpif.h"
      ! arguments
      integer, intent(in) :: ps1, pe1, ps2, pe2, ps3, pe3, &
                             ds1, de1, ds2, de2, ds3, de3
      real(8), dimension(ps1:pe1,ps2:pe2,ps3:pe3), intent(in) :: patch_array
      real(8), dimension(ds1:de1,ds2:de2,ds3:de3), intent(inout) :: domain_array
  
      ! local variables
      integer :: i, ii, j, jj, k, kk, m
      integer, dimension(2) :: idims, jdims
      integer :: mpi_ierr
      integer, dimension(mpi_status_size) :: mpi_stat
      real(8), dimension(:),allocatable ::  av
      
      if (myid .eq. 0) then
  
          do i=0,dims(2,1)-1
            do j=0,dims(1,1)-1
               if (procs(i,j) .ne. 0) then

                  call mpi_recv(jdims, 2, mpi_integer, procs(i,j),mpi_any_tag, local_com, mpi_stat, mpi_ierr)
                  call mpi_recv(idims, 2, mpi_integer, procs(i,j),mpi_any_tag, local_com, mpi_stat, mpi_ierr)

                  allocate(av((idims(2)-idims(1)+1)*(jdims(2)-jdims(1)+1)*(de2-ds2+1))) 
                  call mpi_recv(av, (idims(2)-idims(1)+1)*(jdims(2)-jdims(1)+1)*(de2-ds2+1), &
                                MPI_DOUBLE_PRECISION, procs(i,j), mpi_any_tag, local_com, mpi_stat, mpi_ierr)

                  m=0
                  do jj=jdims(1), jdims(2)
                    do kk=ds2, de2
                      do ii=idims(1), idims(2)
                        m=m+1
                        domain_array(ii,kk,jj)=av(m)
                       enddo
                    enddo
                  enddo
                  deallocate(av)

               else
                  domain_array(ps1:pe1,ps2:pe2,ps3:pe3) = patch_array(ps1:pe1,ps2:pe2,ps3:pe3)
               end if
            end do
           end do
  
      else
  
         jdims(1) = ps3
         jdims(2) = pe3
         call mpi_send(jdims, 2, mpi_integer, 0, myid, local_com,mpi_ierr)
         idims(1) = ps1
         idims(2) = pe1
         call mpi_send(idims, 2, mpi_integer, 0, myid, local_com,mpi_ierr)

         allocate(av((pe1-ps1+1)*(pe2-ps2+1)*(pe3-ps3+1)))
         m=0
         do jj=ps3, pe3
          do kk=ps2, pe2
           do ii=ps1, pe1 
             m = m+1
             av(m)= patch_array(ii,kk,jj)
           enddo
          enddo
         enddo
         call mpi_send(av, (pe1-ps1+1)*(pe2-ps2+1)*(pe3-ps3+1), &
                       MPI_DOUBLE_PRECISION, 0, myid, local_com, mpi_ierr)
         deallocate(av)
      end if
 
   end subroutine gather_whole_field
 end module parall
