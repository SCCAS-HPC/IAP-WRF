module inout2
contains
  subroutine openfileSrcMap(myid,nx,ny,id)
  integer myid,id
  character*1 cdnum
  if ( myid .lt. 0 ) then
      close( 11 )
      return
  endif
  write(cdnum(1:1),'(i1)')id

!  close(110+id)
  open(110+id, file= &
       'emit/mapgrid.d'//cdnum, &
        form='unformatted', access='direct',recl=nx*ny)
  if(myid==0)print *,  &
     ' opened emit/mapgrid.d'//cdnum
  return
  end subroutine

  subroutine openfilSM( myid,nx,ny,irecsm,id, &
                        iyear,imonth,iday, ihour)
  integer myid,id, irecsm
  character*1 cdnum
  character*20 fname
  character*10 date
  write(date(1:4),'(i4)')iyear
  write(date(5:6),'(i2.2)')imonth
  write(date(7:8),'(i2.2)')iday
  write(date(9:10),'(i2.2)')ihour

  write(cdnum(1:1),'(i1)')id
  write (fname(1:3),'(i3.3)') myid

!  close(20+id)
  open( 20+id, err=20, file='initSM/resmd'//cdnum//'.'// &
        date(1:10)//'.grd', form='unformatted',  &
        access='direct',recl=nx*ny, status='old')
  print *, 20+id,'initSM/resmd'//cdnum//'.'//date//'.grd', ' opened'
  irecsm=1
  return

20  irecsm=-1
  return
  end subroutine

subroutine openfil5(myid,nx,ny,id,iyear,imonth,iday, ihour)
  integer myid,id
  character*1 cdnum
  character*20 fname
  character*10 date
  write(date(1:4),'(i4)')iyear
  write(date(5:6),'(i2.2)')imonth
  write(date(7:8),'(i2.2)')iday
  write(date(9:10),'(i2.2)')ihour

  if ( myid .lt. 0 ) then
      close( 11 )
      return
  endif
  write(cdnum(1:1),'(i1)')id
  write (fname(1:3),'(i3.3)') myid
  close(90+id)
  call system("mkdir -p sm")
  call system("mkdir -p sm/tmp")
  call system("mkdir -p sm/tmp/"//date(1:8))
  open( 90+id, file='sm/tmp/'//date(1:8)// &
        '/makd'//cdnum//'.'//fname(1:3)//'.'//date,&
          form='unformatted',status='UNKNOWN' )
  if(ihour==0)then     ! to write restart file
  close(190+id)
  open( 190+id, file='sm/tmp/'//date(1:8)// &
        '/resd'//cdnum//'.'//fname(1:3)//'.'//date,&
          form='unformatted',status='UNKNOWN' )
  endif
  return
  end subroutine

subroutine writesm( myid, a, sx, ex, sy, ey, k,ism,idm,id)
  integer myid, sx, ex, sy, ey, k,ne
  real a(sx-1:ex+1,sy-1:ey+1)

  write(id)myid,k,ism,idm,id,sx,ex,sy,ey
  write(id)((a(i,j), i = sx,ex),j=sy,ey)

  return
  end subroutine

  subroutine read2dmap( myid, a, sx, ex, sy, ey, &
                     nx,ny,id)
  integer myid, sx, ex, sy, ey
  real a(sx-1:ex+1,sy-1:ey+1), v0(nx,ny)
  read(id,rec=1)((v0(i,j),i=1,nx),j=1,ny)

  do i= sx,ex
  do j= sy,ey
   a(i,j)=v0(i,j)
  enddo
  enddo
  return
  end subroutine
end module inout2
