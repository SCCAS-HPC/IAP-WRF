module inout

contains
  subroutine write2d( myid, a, sx, ex, sy, ey, k, id)
  integer myid, sx, ex, sy, ey, k
  real a(sx-1:ex+1,sy-1:ey+1)

  write(10+id)myid,k,sx,ex,sy,ey
  write(10+id)((a(i,j), i = sx,ex),j=sy,ey)
  return
  end subroutine

 subroutine writeterm( myid, a, sx, ex, sy, ey, k, ip, ig, id)
  integer myid, sx, ex, sy, ey, k ,ig, ip
  real a(sx-1:ex+1,sy-1:ey+1)

  write(70+id)myid,k,ip,ig,sx,ex,sy,ey
  write(70+id)((a(i,j), i = sx,ex),j=sy,ey)
  return
  end subroutine

! CCC AQUEOUS DRY
  subroutine writedry( myid, a, sx, ex, sy, ey, k,id)
  integer myid, sx, ex, sy, ey ,k
  real a(sx-1:ex+1,sy-1:ey+1)

  write(60+id) myid,k,sx,ex,sy,ey
  write(60+id)((a(i,j), i = sx,ex),j=sy,ey)

  return
  end subroutine

! CCC AQUEOUS WET
  subroutine writewet( myid, a, sx, ex, sy, ey, k,id)
  integer myid, sx, ex, sy, ey ,k
  real a(sx-1:ex+1,sy-1:ey+1)

  write(290+id) myid,k,sx,ex,sy,ey
  write(290+id)((a(i,j), i = sx,ex),j=sy,ey)

  return
  end subroutine 

! CCC DUST
  subroutine writedust( myid, a, sx, ex, sy, ey, k,id)
  integer myid, sx, ex, sy, ey ,k
  real a(sx-1:ex+1,sy-1:ey+1)

  write(400+id) myid,k,sx,ex,sy,ey
  write(400+id)((a(i,j), i = sx,ex),j=sy,ey)

  return
  end  subroutine

! CCC SEA SALT
  subroutine writesea( myid, a, sx, ex, sy, ey, k,id)
  integer myid, sx, ex, sy, ey ,k
  real a(sx-1:ex+1,sy-1:ey+1)

  write(500+id) myid,k,sx,ex,sy,ey
  write(500+id)((a(i,j), i = sx,ex),j=sy,ey)

  return
  end subroutine

  subroutine read2d( myid, a, sx, ex, sy, ey, &
                     nx,ny,irec,id)
  integer myid, sx, ex, sy, ey
  real a(sx-1:ex+1,sy-1:ey+1), v0(nx,ny)
  read(id,rec=irec)((v0(i,j),i=1,nx),j=1,ny)
  irec=irec+1

  do i= sx,ex
  do j= sy,ey
   a(i,j)=v0(i,j)
  enddo
  enddo
  return
  end subroutine

  subroutine puto3( myid, a, sx, ex, sy, ey, &
                     nx,ny,value)
  integer myid, sx, ex, sy, ey
  real a(sx-1:ex+1,sy-1:ey+1), value

  do i= sx,ex
  do j= sy,ey
   a(i,j)=value
  enddo
  enddo

  return
  end subroutine
   
  subroutine putHG0( myid, a, sx, ex, sy, ey, &
                     nx,ny,ne,k)
  integer myid, sx, ex, sy, ey, ne, k 
  real a(sx-1:ex+1,sy-1:ey+1)
  !!!! this subroutine need to modify if modeling domain is changed !!!!!!!!
  do i= sx,ex
  do j= sy,ey
     if(ne.eq.1) then
       if(j.ge.91) then
       a(i,j) = 1.6-0.8*(k-1)/20.
       else
       a(i,j) = 1.2-0.6*(k-1)/20.
       endif
     else
       a(i,j) = 1.6-0.8*(k-1)/20.
     endif
  enddo
  enddo
  if(myid.eq.0) then
    if(k.eq.1) then
     print *,"Use Global Background HG0 Concentration"
    endif
  endif
  return
  end subroutine

  subroutine openfil1( myid,nx,ny,irec80,irec60,id, &
                        iyear,imonth,iday, ihour)
  integer myid,id, irec80,irec60
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
  
  open( 80+id, err=20, file='init/testd'//cdnum//'.'// &
        date(1:10)//'.grd', form='unformatted',  &
        access='direct',recl=nx*ny, status='old')
  if(myid.eq.0) print *, 80+id,'init/testd'//cdnum//'.'//date//'.grd', ' opened'
  irec80=1

   open (60+id,err=21,file='init/dustd'//cdnum//'.'//&
         date(1:10)//'.grd',form='unformatted',&
         access = 'direct',recl=nx*ny,status='old')
  if(myid.eq.0) print*,60+id, 'init/dustd'//cdnum//'.'//date//'.grd',' opened'
  irec60 = 1  
  return

20 irec80 = -1 
21 irec60 = -1
  return   
  end subroutine

 subroutine openfileEmit( myid,nx,ny,irecg,irecnox,irecairc,irecHgA,irecHgN,id, &
                        iyear,imonth,iday, ihour,ifnox,ifairc,ifHgA,ifHgN)

  integer myid,irecg,irecnox,irecairc,irecHgA,irecHgN,id,ifnox,ifairc,ifHgA,ifHgN
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
  !!!  For main emissions !!!!
!  close(50+id)
  open(50+id,file='/R0/jhe/inputdata/chem/emis_data/emit_main.'//date(5:6)//'.d'//cdnum,&
                    form='unformatted', access='direct',recl=nx*ny,status='old')

  if(myid==0)print *,  &
      'opening emis_data/emit_main.'//date(5:6)//'.d'//cdnum

  irecg=1
  !!!  For soil and lightning NOx !!!!
  IF(ifnox.eq.1) then
!    close(100+id)
    open(100+id,file='/R0/jhe/inputdata/chem/emis_data/emit_nox.'//date(5:6)//'.d'//cdnum,&
                    form='unformatted', access='direct',recl=nx*ny,status='old')
    if(myid==0)print *,  &
      'opening emis_data/emit_nox.'//date(5:6)//'.d'//cdnum

    irecnox=1
  ELSE
    irecnox=-1
  ENDIF 
  !!! For aircraft emissions !!!!
  IF(ifairc.eq.1) then
!    close(150+id)
    open(150+id,file='/R0/jhe/inputdata/chem/emis_data/emit_aircraft.'//date(5:6)//'.d'//cdnum,&
                    form='unformatted', access='direct',recl=nx*ny,status='old')
    if(myid==0)print *,  &
      'opening emis_data/emit_aircraft.'//date(5:6)//'.d'//cdnum

    irecairc=1
  ELSE
    irecairc=-1
  ENDIF
  !!! For Hg anthropogenic emissions !!!!
  IF(ifHgA.eq.1) then
!    close(200+id)
    open(200+id,file='/R0/jhe/inputdata/chem/emis_data/emit_anthro_Hg.'//date(5:6)//'.d'//cdnum,&
                    form='unformatted', access='direct',recl=nx*ny,status='old')
    if(myid==0)print *,  &
      'opening emis_data/emit_anthro_Hg.'//date(5:6)//'.d'//cdnum

    irecHgA=1
  ELSE
    irecHgA=-1
  ENDIF
  !!! For Hg natural emissions !!!!
  IF(ifHgN.eq.1) then
!    close(250+id)
    open(250+id,file='/R0/jhe/inputdata/chem/emis_data/emit_natural_Hg.'//date(5:6)//'.d'//cdnum,&
                    form='unformatted', access='direct',recl=nx*ny,status='old')
    if(myid==0)print *,  &
      'opening emis_data/emit_natural_Hg.'//date(5:6)//'.d'//cdnum

    irecHgN=1
  ELSE
    irecHgN=-1
  ENDIF  

  return
  end subroutine

  subroutine openfil2( myid,nx,ny,irec,irec_as,irecglobal,irechgt,id,&! irectop added by lijie 05-06-3
                       iyear,imonth,iday,ihour,ikosaline)          ! ireckv added by lijie 050712
                                                           ! irecpbl and  irecxmol added by lijie 071203 ACM2
  integer myid,irec,irec_as,id,irecgtop,ireckv,irecglobal,irecpbl
  character*1 cdnum
  character*10 date
  
  write(date(1:4),'(i4)')iyear
  write(date(5:6),'(i2.2)')imonth
  write(date(7:8),'(i2.2)')iday
  write(date(9:10),'(i2.2)')ihour

  write(cdnum(1:1),'(i1)')id
!  close(30+id)
!  open(30+id, file='date(1:6)//'&
!       '/wrfd0'//cdnum//'_'//date(1:4)//'-'//date(5:6)// &
!       '-'//date(7:8)// '_'//date(9:10)//'.dat',  &
!       form='unformatted',access='direct',recl=nx*ny,status='old')
!  if(myid==0) print*, '/MM5orWRFp'//cdnum//'_'//date(1:4)//'-'//date(5:6)// &
!          '-'//date(7:8)// '_'//date(9:10)//'.grd is opened' 
!  irec=1

!  if(ikosaline==2)then
!  close(40+id)
!  open(40+id, file='/name8/kosa/estimate/data/'// &
!       'emitkosa.'//date(1:10), &
!        form='unformatted',access='direct',recl=nx*ny,status='old')
!  endif
!  irec_as = 1

 !---------------------global condition----------------
!  close(300+id)
  open(300+id,file='/R0/jhe/inputdata/chem/global/g2m/g2md0'//cdnum//'_'//date(1:4)//'-'//date(5:6)//&
       '-'//date(7:8)// '_'//date(9:10)//'.DAT',form='unformatted',&
       access='direct',recl=nx*ny,status='old')
  irecglobal=1
  
 !-----------------------------------------------------              
! close(1000) 
! open(1000,file='wrfd01.dat',form='unformatted',&
!      access='direct',recl=nx*ny,status='old')
! irechgt=1 !! modified by chenhs               

  return
  end subroutine
  
subroutine openfil3(myid,nx,ny,irecg,id,&
                    iyear,imonth,iday, ihour)
  integer myid,irec,id, irec80
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
  close(10+id)
  call system("mkdir -p out/tmp")
  call system("mkdir -p out/tmp/"//date(1:8))
  open( 10+id, file='out/tmp/'//date(1:8)// &
        '/food'//cdnum//'.'//fname(1:3)//'.'//date,&
          form='unformatted',status='UNKNOWN' )
!lijie change tmp to tmp1 050622 for run 2 apops.new at the same time

  close(70+id)
  call system("mkdir -p term/tmp")
  call system("mkdir -p term/tmp/"//date(1:8))
  open( 70+id, file='term/tmp/'//date(1:8)// &
         '/term'//cdnum//'.'//fname(1:3)//'.'//date,&
          form='unformatted',status='UNKNOWN' )
  return
  end subroutine

subroutine openfil4(myid,nx,ny,irecg,id,&
                      iyear,imonth,iday, ihour)
  integer irec,id, irec80
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
  close(10+id)
  call system("mkdir -p out/tmp")
  call system("mkdir -p out/tmp/"//date(1:8))
  open( 10+id, file='out/tmp/'//date(1:8)// &
        '/food'//cdnum//'.'//fname(1:3)//'.'//date,&
        form='unformatted',access='direct',recl=nx*ny,status='UNKNOWN')
  irecg=1
  return
  end subroutine

! CCCC DUST CONCENTRATION AND EMISSIONS, DRY, WET and GRAVITY DEPOSITION
  subroutine openfildust(myid,nx,ny,irecdust,id,iyear,imonth,iday,ihour)
   integer myid,irecdust,id
   character*1 cdnum
   character*20 fname
   character*10 date
    write(date(1:4),'(i4)')iyear
    write(date(5:6),'(i2.2)')imonth
    write(date(7:8),'(i2.2)')iday
    write(date(9:10),'(i2.2)')ihour

    write(cdnum(1:1),'(i1)')id
    write (fname(1:3),'(i3.3)') myid
    close(400+id)
    call system("mkdir -p dust/tmp")
    call system("mkdir -p dust/tmp/"//date(1:8))
    open( 400+id, file='dust/tmp/'//date(1:8)// &
        '/dustd'//cdnum//'.'//fname(1:3)//'.'//date,&
          form='unformatted',status='UNKNOWN' )
    irecdust = 1
    return
  end subroutine

! CCCC SEA CONCENTRATION AND EMISSIONS, DRY, WET and GRAVITY DEPOSITION
  subroutine openfilsea(myid,nx,ny,irecsea,id,iyear,imonth,iday,ihour)
   integer myid,irecsea,id
   character*1 cdnum
   character*20 fname
   character*10 date
    write(date(1:4),'(i4)')iyear
    write(date(5:6),'(i2.2)')imonth
    write(date(7:8),'(i2.2)')iday
    write(date(9:10),'(i2.2)')ihour

    write(cdnum(1:1),'(i1)')id
    write (fname(1:3),'(i3.3)') myid
    close(500+id)
    call system("mkdir -p seasalt/tmp")
    call system("mkdir -p seasalt/tmp/"//date(1:8))
    open( 500+id, file='seasalt/tmp/'//date(1:8)// &
        '/seasaltd'//cdnum//'.'//fname(1:3)//'.'//date,&
          form='unformatted',status='UNKNOWN' )
    irecsea = 1
    return
  end subroutine

!CCCC   AQUEOUS DEPOSITION  CCCCCCCCCC
  subroutine openfilwet(myid,nx,ny,irecwet,id,iyear,imonth,iday,ihour)
   integer myid,irecwet,id
   character*1 cdnum
   character*20 fname
   character*10 date
    write(date(1:4),'(i4)')iyear
    write(date(5:6),'(i2.2)')imonth
    write(date(7:8),'(i2.2)')iday
    write(date(9:10),'(i2.2)')ihour

    write(cdnum(1:1),'(i1)')id
    write (fname(1:3),'(i3.3)') myid
    close(290+id)
    call system("mkdir -p wetdep/tmp")
    call system("mkdir -p wetdep/tmp/"//date(1:8))
    open( 290+id, file='wetdep/tmp/'//date(1:8)// &
        '/wetdepd'//cdnum//'.'//fname(1:3)//'.'//date,&
          form='unformatted',status='UNKNOWN' )
    irecwet = 1
    return
  end subroutine
 
  subroutine openfildry(myid,nx,ny,irecdry,id,iyear,imonth,iday,ihour)
   integer myid,irecdry,id
   character*1 cdnum
   character*20 fname
   character*10 date
    write(date(1:4),'(i4)')iyear
    write(date(5:6),'(i2.2)')imonth
    write(date(7:8),'(i2.2)')iday
    write(date(9:10),'(i2.2)')ihour

    write(cdnum(1:1),'(i1)')id
    write (fname(1:3),'(i3.3)') myid
    close(60+id)
    call system("mkdir -p drydep/tmp")
    call system("mkdir -p drydep/tmp/"//date(1:8))
    open( 60+id, file='drydep/tmp/'//date(1:8)// &
        '/drydepd'//cdnum//'.'//fname(1:3)//'.'//date,&
          form='unformatted',status='UNKNOWN' )
    irecdry = 1
    return
  end subroutine
  
  subroutine openfilMOZART( myid,nx,ny,irecMOZART,id,iyear,imonth,iday,ihour)    
  integer myid,irecMOZART
  character*1 cdnum
  character*10 date

  write(date(1:4),'(i4)')iyear
  write(date(5:6),'(i2.2)')imonth
  write(date(7:8),'(i2.2)')iday
  write(date(9:10),'(i2.2)')ihour

  write(cdnum(1:1),'(i1)')id
!---------------------global initial condition----------------
!  close(350+id)
  open(350+id,err=10,file='/R0/jhe/inputdata/chem/global/g2m/g2md0'//cdnum//'_'//date(1:4)//'-'//date(5:6)//&
       '-'//date(7:8)// '_'//date(9:10)//'.DAT',form='unformatted',&
       access='direct',recl=nx*ny,status='old')
  if(myid.eq.0) print *, 350+id,'MOZART:global/g2m/g2md0'//cdnum//'_'//date(1:4)//'-'//date(5:6)//&
       '-'//date(7:8)// '_'//date(9:10)//'.DAT'
   irecMOZART=1
!-----------------------------------------------------
   return
10 irecMOZART = -1 ! we don't read the MOZART
  return
  end subroutine
    
  subroutine needSNDwrite( a, sx, ex, sy, ey, k,b,ntotalnum,nzz)
  integer sx, ex, sy, ey, k,  ntotalnum
  real a(sx-1:ex+1,sy-1:ey+1),b(ntotalnum,nzz)

  if(ntotalnum .ne. (ex-sx+1)*(ey-sy+1))then
    print *,'number setting error'
   stop
  endif

  ii=1
  do i=sx,ex
  do j=sy,ey
    b(ii,k)=a(i,j)
    ii=ii+1
    enddo
    enddo
  return
  end subroutine

  subroutine needGETwrite(a,nx,ny,sx,ex,sy,ey,k,b,ntotalnum,nzz)
  integer sx, ex, sy, ey, k, ntotalnum
  real a(nx,ny,nzz),b(ntotalnum,nzz)

  if(ntotalnum .ne. (ex-sx+1)*(ey-sy+1))then
    print *,'number setting error'
   stop
  endif

  ii=1
  do i=sx,ex
  do j=sy,ey
    a(i,j,k)=b(ii,k)
    ii=ii+1
    enddo
    enddo
  return
  end subroutine

  subroutine putvalue( myid, a, sx, ex, sy, ey, &
                     ix,iy,value)
  integer myid, sx, ex, sy, ey,ix,iy
  real a(sx-1:ex+1,sy-1:ey+1), value
   a(ix,iy)=value
  return
  end subroutine  

  subroutine getvalue( myid, a, sx, ex, sy, ey, &
                     ix,iy,value)
  integer myid, sx, ex, sy, ey, ix, iy
  real a(sx-1:ex+1,sy-1:ey+1), value
  
  value=a(ix,iy)

  return
  end subroutine

  subroutine write2dconv( myid, a, sx, ex, sy, ey, k, id)
  integer myid, sx, ex, sy, ey, k
  integer a(sx-1:ex+1,sy-1:ey+1)
     
  write(10+id)myid,k,sx,ex,sy,ey
  write(10+id)((float(a(i,j)), i = sx,ex),j=sy,ey)
  return
  end subroutine

  subroutine read2dconv( myid, a, sx, ex, sy, ey, &
                       nx,ny,irec,id)
  integer myid, sx, ex, sy, ey
  integer a(sx-1:ex+1,sy-1:ey+1)
  real v0(nx,ny)
   read(id,rec=irec)((v0(i,j),i=1,nx),j=1,ny)
   irec=irec+1
   do j= sy,ey
   do i= sx,ex
     a(i,j)=int(v0(i,j))
   enddo
   enddo
     return
  end subroutine
end module inout
