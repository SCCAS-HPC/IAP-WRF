  module inout3
  use pio, only : pio_write_darray, io_desc_t, file_desc_t, var_desc_t, iosystem_desc_t, &
      pio_int, pio_real, pio_double, pio_noerr, &
      PIO_iotype_binary, pio_iotype_netcdf, &
      pio_char, pio_write, pio_clobber,pio_noclobber,&
      pio_initdecomp,  pio_openfile, pio_closefile, pio_createfile,pio_freedecomp, &
      pio_finalize,  PIO_enddef,  PIO_def_dim,  PIO_def_var, PIO_put_att, PIO_put_var
  type(iosystem_desc_t), pointer, private :: pio_subsystem
  contains        

  subroutine output_atm(iyear, imonth, iday, ihour, ne)
   use geatm_vartype, only:sx,ex,sy,ey,&
       nx,ny,nz,nzz,isize,igas,iaer,ntt,&
       ip2mem,ip3mem,ip4mem,ip5mem,initprintgas,printgas, &
       u,v,t,w,rh1,plev,heiz,rhsfc,u10,v10,ust,&
       dz, tropp, ssa, aod, cldopd, duso2, &
       duo3, duno2, ana, aso4, anh4, ano3, &
       acl, cph, ope, aer, gas, &
       jo1d, jno2, uvb, uvbs, uva, vis
  
   integer :: iyear, imonth, iday, ihour, ne
   
   type(file_desc_t) :: file
   type(io_desc_t) :: iodesc2d,iodesc3d,iodesc4d,iodesc5d
   type(var_desc_t), pointer :: varid_x2a(:)
   integer :: i, j, k, i0, ig, is, ia, km, dim2d(2), dim3d(3), dim4d(4),dim5d(5)
   real :: xlat(ny(ne)),xlon(nx(ne)),sigma(nzz),b(isize),c(iaer),fill_value
   real, dimension(:,:),allocatable :: chem2d
   real, dimension(:,:,:),allocatable :: chem3d
   integer, pointer :: ldof2d(:),ldof3d(:)
   integer         :: m,mlen, mm        ! numbr of variables
   integer         :: rcode        ! return error code
   character*10 date
   character*30 fname
  
   date(1:4)=char(iyear/1000+48)//char(iyear/100-iyear/1000*10+48)//char(iyear/10-iyear/100*10+48)//&
             char(iyear-iyear/10*10+48)
   date(5:6)=char(imonth/10+48)//char(imonth-imonth/10*10+48)
   date(7:8)=char(iday/10+48)//char(iday-iday/10*10+48)
   date(9:10)=char(ihour/10+48)//char(ihour-ihour/10*10+48) 
   fillvalue = -99999.9
   
   do i=1,ny(ne)
   xlat(i)=(i-1)*1.0-89.5
   enddo
   do i=1,nx(ne)
   xlon(i)=(i-1)*1.0+0.5
   enddo
   do i=1,nzz
   sigma(i)=i
   enddo

   allocate(varid_x2a(16))
   
   allocate(chem2d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1))
   allocate(chem3d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1,nz(ne)))
   
   fname='geatm.atm.'//date//'.nc'
   rcode = pio_createfile(pio_subsystem, file, pio_iotype_netcdf, trim(adjustl(fname)), PIO_CLOBBER)
  
   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)
   allocate(ldof2d(m)) 
   m=0
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof2d(m)=(j-1)*nx(ne)+i
   end do
   end do

   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)*nz(ne)
   allocate(ldof3d(m)) 
   m=0
   do k=1,nz(ne)
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof3d(m)=(j-1)*nx(ne)+i+(k-1)*nx(ne)*ny(ne)
   end do
   end do
   end do

   dim2d=(/nx(ne),ny(ne)/)
   dim3d=(/nx(ne),ny(ne),nz(ne)/)
   call pio_initdecomp(pio_subsystem, pio_real, dim2d, ldof2d, iodesc2d)
   call pio_initdecomp(pio_subsystem, pio_real, dim3d, ldof3d, iodesc3d)   
   
   rcode = pio_def_dim(File,'lon',nx(ne),dim3d(1))
   rcode = pio_def_dim(File,'lat',ny(ne),dim3d(2))
   dim2d(1:2)=dim3d(1:2)
   rcode = pio_def_dim(File,'lev',nzz,dim3d(3))
   
   rcode = pio_def_var(File,'lon',PIO_real,dim3d(1:1),varid_x2a(1))
   rcode = pio_put_att (File, varid_x2a(1), 'long_name', 'longitude')
   rcode = pio_put_att (File, varid_x2a(1), 'units', 'degrees_east') 
   rcode = pio_def_var(File,'lat',PIO_REAL,dim3d(2:2),varid_x2a(2))
   rcode = pio_put_att (File, varid_x2a(2), 'long_name', 'latitude')
   rcode = pio_put_att (File, varid_x2a(2), 'units', 'degrees_north')
   rcode = pio_def_var(File,'lev',PIO_REAL,dim3d(3:3),varid_x2a(3))
   rcode = pio_def_var(File,'dz',PIO_REAL,dim3d,varid_x2a(4))   
   rcode = pio_def_var(File,'u',PIO_REAL,dim3d,varid_x2a(5))   
   rcode = pio_def_var(File,'v',PIO_REAL,dim3d,varid_x2a(6))   
   rcode = pio_def_var(File,'w',PIO_REAL,dim3d,varid_x2a(7))
   rcode = pio_def_var(File,'t',PIO_REAL,dim3d,varid_x2a(8))   
   rcode = pio_def_var(File,'rh1',PIO_REAL,dim3d,varid_x2a(9))
   rcode = pio_def_var(File,'plev',PIO_REAL,dim3d,varid_x2a(10))   
   rcode = pio_def_var(File,'heiz',PIO_REAL,dim3d,varid_x2a(11))
   rcode = pio_def_var(File,'tropp',PIO_REAL,dim2d,varid_x2a(12))
   rcode = pio_def_var(File,'rhsfc',PIO_REAL,dim2d,varid_x2a(13))
   rcode = pio_def_var(File,'u10',PIO_REAL,dim2d,varid_x2a(14))
   rcode = pio_def_var(File,'v10',PIO_REAL,dim2d,varid_x2a(15))
   rcode = pio_def_var(File,'ust',PIO_REAL,dim2d,varid_x2a(16))

   mm=16
   do i=1,mm
   rcode = pio_put_att(File,varid_x2a(i),"_fillvalue",fillvalue)
   end do
   rcode = pio_enddef(File) 
   
   rcode = pio_put_var(File, varid_x2a(1), xlon)
   rcode = pio_put_var(File, varid_x2a(2), xlat)
   rcode = pio_put_var(File, varid_x2a(3), sigma)
   
   ! dz
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
          do j= sy(ne),ey(ne)
          do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=dz(km) 
          end do
          end do
   end do  
   call pio_write_darray(File, varid_x2a(4), iodesc3d,chem3d, rcode)
   ! u
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
          do j= sy(ne),ey(ne)
          do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=u(km)                        
          end do
          end do
   end do     
   call pio_write_darray(File, varid_x2a(5), iodesc3d, chem3d, rcode)
   ! v
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
          do j= sy(ne),ey(ne)
          do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=v(km)                        
          end do
          end do
   end do        
   call pio_write_darray(File, varid_x2a(6), iodesc3d, chem3d, rcode)
   ! w
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
          do j= sy(ne),ey(ne)
          do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=w(km)                        
          end do
          end do
   end do        
   call pio_write_darray(File, varid_x2a(7), iodesc3d, chem3d, rcode)
   ! t
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
          do j= sy(ne),ey(ne)
          do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=t(km)                        
          end do
          end do
   end do        
   call pio_write_darray(File, varid_x2a(8), iodesc3d, chem3d, rcode)
   ! rh
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
          do j= sy(ne),ey(ne)
          do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=rh1(km)                        
          end do
          end do
   end do        
   call pio_write_darray(File, varid_x2a(9), iodesc3d, chem3d, rcode)
   ! p
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
          do j= sy(ne),ey(ne)
          do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=plev(km)                        
          end do
          end do
   end do        
   call pio_write_darray(File, varid_x2a(10), iodesc3d, chem3d, rcode)
   ! z
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
          do j= sy(ne),ey(ne)
          do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=heiz(km) 
          end do
          end do
   end do  
   call pio_write_darray(File, varid_x2a(11), iodesc3d, chem3d, rcode)
   ! tropp   
        i0=ip2mem(ne)
          do j= sy(ne),ey(ne)
          do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=tropp(km)                        
          end do
          end do   
   call pio_write_darray(File, varid_x2a(12), iodesc2d, chem2d, rcode)
   ! rhsfc
        i0=ip2mem(ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=rhsfc(km)
        end do
        end do
   call pio_write_darray(File, varid_x2a(13), iodesc2d, chem2d, rcode)
   ! u10
        i0=ip2mem(ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=u10(km)
        end do
        end do
   call pio_write_darray(File, varid_x2a(14), iodesc2d, chem2d, rcode)
   ! v10
        i0=ip2mem(ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=v10(km)
        end do
        end do
   call pio_write_darray(File, varid_x2a(15), iodesc2d, chem2d, rcode)
   ! ust
        i0=ip2mem(ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=ust(km)
        end do
        end do
   call pio_write_darray(File, varid_x2a(16), iodesc2d, chem2d, rcode)


   call pio_freedecomp(File,iodesc2d)
   call pio_freedecomp(File,iodesc3d)
   call pio_closefile(file)     

   deallocate(ldof2d)
   deallocate(ldof3d)
   deallocate(chem2d)
   deallocate(chem3d)
   deallocate(varid_x2a)
   end subroutine output_atm
!--------------------------------------------------------------------------------
  subroutine output3(iyear, imonth, iday, ihour, it1, ne)
   use geatm_vartype, only:sx,ex,sy,ey,&
       nx,ny,nz,nzz,isize,igas,iaer,ntt,&
       ip2mem,ip3mem,ip4mem,ip5mem,initprintgas,printgas, &
       u,v,t,w,rh1,plev,heiz,rhsfc,u10,v10,ust,&
       dz, tropp, ssa, aod, cldopd, duso2, &
       duo3, duno2, ana, aso4, anh4, ano3, &
       acl, cph, ope, aer, gas, &
       jo1d, jno2, uvb, uvbs, uva, vis
  
   integer :: iyear, imonth, iday, ihour, it1,  ne
   
   type(file_desc_t) :: file
   type(io_desc_t) :: iodesc2d,iodesc3d,iodesc4d,iodesc5d
   type(var_desc_t), pointer :: varid_x2a(:)
   integer :: i, j, k, i0, ig, is, ia, km, dim2d(2), dim3d(3), dim4d(4),dim5d(5)
   real :: xlat(ny(ne)),xlon(nx(ne)),sigma(nzz),b(isize),c(iaer),fill_value
   real, dimension(:,:),allocatable :: chem2d
   real, dimension(:,:,:),allocatable :: chem3d
   real, dimension(:,:,:,:,:),allocatable :: chem5d
   integer, pointer :: ldof2d(:),ldof3d(:),ldof5d(:)
   integer         :: m,mlen, mm        ! numbr of variables
   integer         :: rcode        ! return error code
   character*10 date
   character*30 fname
  
   date(1:4)=char(iyear/1000+48)//char(iyear/100-iyear/1000*10+48)//char(iyear/10-iyear/100*10+48)//&
             char(iyear-iyear/10*10+48)
   date(5:6)=char(imonth/10+48)//char(imonth-imonth/10*10+48)
   date(7:8)=char(iday/10+48)//char(iday-iday/10*10+48)
   date(9:10)=char(ihour/10+48)//char(ihour-ihour/10*10+48) 
   fillvalue = -99999.9
   
   do i=1,ny(ne)
   xlat(i)=(i-1)*1.0-89.5
   enddo
   do i=1,nx(ne)
   xlon(i)=(i-1)*1.0+0.5
   enddo
   do i=1,nzz
   sigma(i)=i
   enddo

   do i=1,isize
   b(i)=i
   enddo
   do i=1,iaer
   c(i)=i
   enddo
      
   allocate(varid_x2a(19))
   
   allocate(chem2d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1))
   allocate(chem3d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1,nz(ne)))
   allocate(chem5d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1,nz(ne),isize,iaer))
   
   fname='geatm.'//date//'.nc'
   rcode = pio_createfile(pio_subsystem, file, pio_iotype_netcdf, trim(adjustl(fname)), PIO_CLOBBER)
  
   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)
   allocate(ldof2d(m)) 
   m=0
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof2d(m)=(j-1)*nx(ne)+i
   end do
   end do

   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)*nz(ne)
   allocate(ldof3d(m)) 
   m=0
   do k=1,nz(ne)
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof3d(m)=(j-1)*nx(ne)+i+(k-1)*nx(ne)*ny(ne)
   end do
   end do
   end do

   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)*isize*nz(ne)*iaer
   allocate(ldof5d(m))
   m=0
   do ia=1,iaer
   do is=1,isize
   do k=1,nz(ne)
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof5d(m)=i+(j-1)*nx(ne)+(k-1)*nx(ne)*ny(ne)+(is-1)*nx(ne)*ny(ne)*nz(ne)+&
             (ia-1)*nx(ne)*ny(ne)*nz(ne)*isize
   end do
   end do
   end do
   end do
   end do

   dim2d=(/nx(ne),ny(ne)/)
   dim3d=(/nx(ne),ny(ne),nz(ne)/)
   dim5d=(/nx(ne),ny(ne),nz(ne),isize,iaer/)
   call pio_initdecomp(pio_subsystem, pio_real, dim2d, ldof2d, iodesc2d)
   call pio_initdecomp(pio_subsystem, pio_real, dim3d, ldof3d, iodesc3d)   
   call pio_initdecomp(pio_subsystem, pio_real, dim5d, ldof5d, iodesc5d)
   
   rcode = pio_def_dim(File,'lon',nx(ne),dim3d(1))
   rcode = pio_def_dim(File,'lat',ny(ne),dim3d(2))
   dim2d(1:2)=dim3d(1:2)
   rcode = pio_def_dim(File,'lev',nzz,dim3d(3))
   dim5d(1:3)=dim3d(1:3)
   rcode = pio_def_dim(File,'isize',isize,dim5d(4))
   rcode = pio_def_dim(File,'iaer',iaer,dim5d(5))
   
   rcode = pio_def_var(File,'lon',PIO_real,dim3d(1:1),varid_x2a(1))
   rcode = pio_put_att (File, varid_x2a(1), 'long_name', 'longitude')
   rcode = pio_put_att (File, varid_x2a(1), 'units', 'degrees_east') 
   rcode = pio_def_var(File,'lat',PIO_REAL,dim3d(2:2),varid_x2a(2))
   rcode = pio_put_att (File, varid_x2a(2), 'long_name', 'latitude')
   rcode = pio_put_att (File, varid_x2a(2), 'units', 'degrees_north')
   rcode = pio_def_var(File,'lev',PIO_REAL,dim3d(3:3),varid_x2a(3))
   rcode = pio_def_var(File,'isize',PIO_REAL,dim5d(4:4),varid_x2a(4))
   rcode = pio_def_var(File,'iaer',PIO_REAL,dim5d(5:5),varid_x2a(5))
   rcode = pio_def_var(File,'ssa',PIO_REAL,dim3d,varid_x2a(6))
   rcode = pio_def_var(File,'aod',PIO_REAL,dim2d,varid_x2a(7))
   rcode = pio_def_var(File,'cldopd',PIO_REAL,dim2d,varid_x2a(8))
   rcode = pio_def_var(File,'duso2',PIO_REAL,dim2d,varid_x2a(9))
   rcode = pio_def_var(File,'duo3',PIO_REAL,dim2d,varid_x2a(10))
   rcode = pio_def_var(File,'duno2',PIO_REAL,dim2d,varid_x2a(11))
   rcode = pio_def_var(File,'ana',PIO_REAL,dim3d,varid_x2a(12))
   rcode = pio_def_var(File,'aso4',PIO_REAL,dim3d,varid_x2a(13))
   rcode = pio_def_var(File,'anh4',PIO_REAL,dim3d,varid_x2a(14))
   rcode = pio_def_var(File,'ano3',PIO_REAL,dim3d,varid_x2a(15))
   rcode = pio_def_var(File,'acl',PIO_REAL,dim3d,varid_x2a(16))
   rcode = pio_def_var(File,'cph',PIO_REAL,dim3d,varid_x2a(17))
   rcode = pio_def_var(File,'ope',PIO_REAL,dim3d,varid_x2a(18))
   rcode = pio_def_var(File,'aer',PIO_REAL,dim5d,varid_x2a(19))

   mm=19
   do i=1,mm
   rcode = pio_put_att(File,varid_x2a(i),"_fillvalue",fillvalue)
   end do
   rcode = pio_enddef(File) 
   
   rcode = pio_put_var(File, varid_x2a(1), xlon)
   rcode = pio_put_var(File, varid_x2a(2), xlat)
   rcode = pio_put_var(File, varid_x2a(3), sigma)
   rcode = pio_put_var(File, varid_x2a(4), b)
   rcode = pio_put_var(File, varid_x2a(5), c)
   
   ! ssa
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=ssa(km)                           
        end do
        end do
   end do
   call pio_write_darray(File, varid_x2a(6), iodesc3d, chem3d, rcode)
   ! aod
        i0=ip2mem(ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=aod(km)                           
        end do
        end do
   call pio_write_darray(File, varid_x2a(7), iodesc2d, chem2d, rcode)   
   ! cldopd
        i0=ip2mem(ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=cldopd(km)                           
        end do
        end do
   call pio_write_darray(File, varid_x2a(8), iodesc2d, chem2d, rcode)
   ! duso2   
        i0=ip2mem(ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=duso2(km)                           
        end do
        end do
   call pio_write_darray(File, varid_x2a(9), iodesc2d, chem2d, rcode)
   ! duo3   
        i0=ip2mem(ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=duo3(km)                           
        end do
        end do        
   call pio_write_darray(File, varid_x2a(10), iodesc2d, chem2d, rcode)
   ! duno2   
        i0=ip2mem(ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=duno2(km)                           
        end do
        end do
   call pio_write_darray(File, varid_x2a(11), iodesc2d, chem2d, rcode)
   ! ana
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=ana(km)                           
        end do
        end do
   end do
   call pio_write_darray(File, varid_x2a(12), iodesc3d, chem3d, rcode)
   ! aso4   
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=aso4(km)                           
        end do
        end do
   end do  
   call pio_write_darray(File, varid_x2a(13), iodesc3d, chem3d, rcode)
   ! anh4   
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=anh4(km)                           
        end do
        end do
   end do  
   call pio_write_darray(File, varid_x2a(14), iodesc3d, chem3d, rcode)
   ! sno3   
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=ano3(km)                           
        end do
        end do
   end do  
   call pio_write_darray(File, varid_x2a(15), iodesc3d, chem3d, rcode)
   ! acl
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=acl(km) 
        end do
        end do
   end do  
   call pio_write_darray(File, varid_x2a(16), iodesc3d, chem3d, rcode)
   ! cph   
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=cph(km)                           
        end do
        end do
   end do  
   call pio_write_darray(File, varid_x2a(17), iodesc3d, chem3d, rcode)
   ! ope   
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=ope(km)                           
        end do
        end do
   end do  
   call pio_write_darray(File, varid_x2a(18), iodesc3d, chem3d, rcode)
   ! aer
   do ia=1,iaer
   do is=1,isize   
   do k=1,nz(ne)        
        i0=ip5mem(k,is,ia,ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem5d(i-sx(ne)+1,j-sy(ne)+1,k,is,ia)=aer(km)                           
        end do
        end do
   end do
   end do
   end do
   call pio_write_darray(File, varid_x2a(19), iodesc5d, chem5d, rcode)   

   call pio_freedecomp(File,iodesc2d)
   call pio_freedecomp(File,iodesc3d)
   call pio_freedecomp(File,iodesc5d)
   call pio_closefile(file)     

   deallocate(ldof2d)
   deallocate(ldof3d)
   deallocate(ldof5d)
   deallocate(chem2d)
   deallocate(chem3d)
   deallocate(chem5d)   
   deallocate(varid_x2a)
   end subroutine output3
!--------------------------------------------------------------------------------
  subroutine output_gas(iyear, imonth, iday, ihour, it1, ne)
   use geatm_vartype, only:sx,ex,sy,ey,&
       nx,ny,nz,nzz,isize,igas,iaer,ntt,&
       ip2mem,ip3mem,ip4mem,ip5mem,initprintgas,printgas, &
       u,v,t,w,rh1,plev,heiz,rhsfc,u10,v10,ust,&
       dz, tropp, ssa, aod, cldopd, duso2, &
       duo3, duno2, ana, aso4, anh4, ano3, &
       acl, cph, ope, aer, gas, &
       jo1d, jno2, uvb, uvbs, uva, vis
  
   integer :: iyear, imonth, iday, ihour, it1,  ne
   
   type(file_desc_t) :: file
   type(io_desc_t) :: iodesc3d
   type(var_desc_t), pointer :: varid_x2a(:)
   integer :: i, j, k, i0, ig, is, ia, km, dim3d(3)
   real :: xlat(ny(ne)),xlon(nx(ne)),sigma(nzz),fill_value
   real, dimension(:,:,:),allocatable :: chem3d
   integer, pointer :: ldof3d(:)
   integer         :: m,mlen, mm        ! numbr of variables
   integer         :: rcode        ! return error code
   character*10 date
   character*30 fname
   character*3 cgas(igas)
  
   date(1:4)=char(iyear/1000+48)//char(iyear/100-iyear/1000*10+48)//char(iyear/10-iyear/100*10+48)//&
             char(iyear-iyear/10*10+48)
   date(5:6)=char(imonth/10+48)//char(imonth-imonth/10*10+48)
   date(7:8)=char(iday/10+48)//char(iday-iday/10*10+48)
   date(9:10)=char(ihour/10+48)//char(ihour-ihour/10*10+48) 
   do i=1,igas
   cgas(i)=char(i/100+48)//char((i-i/100*100)/10+48)//char(i-i/10*10+48)
   end do
   fillvalue = -99999.9
   
   do i=1,ny(ne)
   xlat(i)=(i-1)*1.0-89.5
   enddo
   do i=1,nx(ne)
   xlon(i)=(i-1)*1.0+0.5
   enddo
   do i=1,nz(ne)
   sigma(i)=i
   enddo

   mlen=0
   if(it1.eq.ntt) then
     do i=1,igas
     if(InitPrintGas(i)==1)then
     mlen=mlen+1
     endif
     enddo
   else
     do i=1,igas
     if(PrintGas(i)==1)then
     mlen=mlen+1
     endif
     enddo
   end if

   allocate(varid_x2a(3+mlen))
   allocate(chem3d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1,nz(ne)))
   
   fname='geatm.gas.'//date//'.nc'
   rcode = pio_createfile(pio_subsystem, file, pio_iotype_netcdf, trim(adjustl(fname)), PIO_CLOBBER)
  
   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)*nz(ne)
   allocate(ldof3d(m))
   m=0
   do k=1,nz(ne)
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof3d(m)=(j-1)*nx(ne)+i+(k-1)*nx(ne)*ny(ne)
   end do
   end do
   end do

   dim3d=(/nx(ne),ny(ne),nz(ne)/)
   call pio_initdecomp(pio_subsystem, pio_real, dim3d, ldof3d, iodesc3d)
   
   rcode = pio_def_dim(File,'lon',nx(ne),dim3d(1))
   rcode = pio_def_dim(File,'lat',ny(ne),dim3d(2))
   rcode = pio_def_dim(File,'lev',nz(ne),dim3d(3))
   
   rcode = pio_def_var(File,'lon',PIO_real,dim3d(1:1),varid_x2a(1))
   rcode = pio_put_att (File, varid_x2a(1), 'long_name', 'longitude')
   rcode = pio_put_att (File, varid_x2a(1), 'units', 'degrees_east') 
   rcode = pio_def_var(File,'lat',PIO_REAL,dim3d(2:2),varid_x2a(2))
   rcode = pio_put_att (File, varid_x2a(2), 'long_name', 'latitude')
   rcode = pio_put_att (File, varid_x2a(2), 'units', 'degrees_north')
   rcode = pio_def_var(File,'lev',PIO_REAL,dim3d(3:3),varid_x2a(3))


   mlen=0
   if(it1.eq.ntt) then
    do ig=1,igas
     if(InitPrintGas(ig)==1)then
     mlen=mlen+1
     rcode = pio_def_var(File,'gas'//cgas(ig),PIO_REAL,dim3d,varid_x2a(3+mlen))
     endif
    end do
   else
    do ig=1,igas
     if(PrintGas(ig)==1)then
     mlen=mlen+1
     rcode = pio_def_var(File,'gas'//cgas(ig),PIO_REAL,dim3d,varid_x2a(3+mlen))
     end if
    end do
   endif

  
   do i=1,3+mlen
   rcode = pio_put_att(File,varid_x2a(i),"_fillvalue",fillvalue)
   end do
   rcode = pio_enddef(File) 
   
   rcode = pio_put_var(File, varid_x2a(1), xlon)
   rcode = pio_put_var(File, varid_x2a(2), xlat)
   rcode = pio_put_var(File, varid_x2a(3), sigma)
   
   ! gas
   mlen=0
   if(it1.eq.ntt) then
    do ig=1,igas        
     if(InitPrintGas(ig)==1)then        
     mlen=mlen+1    
     do k=1,nz(ne)
        i0=ip4mem(k,ig,ne)   
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
           km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
           chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=gas(km)
        end do
        end do        
     end do
     call pio_write_darray(File, varid_x2a(3+mlen), iodesc3d, chem3d, rcode)
     endif
    end do
   else
    do ig=1,igas        
     if(PrintGas(ig)==1)then        
     mlen=mlen+1
     do k=1,nz(ne)
        i0=ip4mem(k,ig,ne)   
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
           km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
           chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=gas(km)
        end do
        end do        
     end do
     call pio_write_darray(File, varid_x2a(3+mlen), iodesc3d, chem3d, rcode)
     end if
    end do
   endif

   call pio_freedecomp(File,iodesc3d)
   call pio_closefile(file)     

   deallocate(ldof3d)
   deallocate(chem3d)
   deallocate(varid_x2a)
   end subroutine output_gas
!--------------------------------------------------------------------

  subroutine output_dust(iyear, imonth, iday, ihour, ne)
   use geatm_vartype, only:sx,ex,sy,ey,&
       nx,ny,nz,nzz,isize,ndustcom,&
       ip2mem,ip3mem,ifdustcom,&
       DUSTEXT,DUSTAOD,DUSTEMISS,DUSTDRY,DUSTGRAV,DUSTWET,EMITFACT,DUSTCOMP,&
       DUSTDRYSO4,DUSTDRYNO3,DUSTDRYFeII,DUSTDRYFeIII,DUSTWETSO4,DUSTWETNO3,&
       DUSTWETFeII,DUSTWETFeIII,DUSTGRAVSO4,DUSTGRAVNO3,DUSTGRAVFeII,DUSTGRAVFeIII
      
   integer :: iyear, imonth, iday, ihour, ne
   
   type(file_desc_t) :: file
   type(io_desc_t) :: iodesc2d,iodesc3d,iodesc4d,iodesc5d
   type(var_desc_t), pointer :: varid_x2a(:)
   integer :: i, j, k, i0, ig, is, ia, dim2d(2), dim3d(3), dim4d(4),dim5d(5)
   integer :: mem3d, mem2d
   real :: xlat(ny(ne)), xlon(nx(ne)), sigma(nzz), a(isize), b(ndustcom), fill_value
   real, dimension(:,:),allocatable :: chem2d
   real, dimension(:,:,:),allocatable :: chem3d
!   real, dimension(:,:,:,:,:),allocatable :: chem5d
   integer, pointer :: ldof2d(:),ldof3d(:)
!   integer, pointer :: ldof5d(:)
   integer         :: m, mm        ! numbr of variables
   integer         :: rcode        ! return error code
   character*10 date
   character*30 fname

   date(1:4)=char(iyear/1000+48)//char(iyear/100-iyear/1000*10+48)//char(iyear/10-iyear/100*10+48)//&
             char(iyear-iyear/10*10+48)
   date(5:6)=char(imonth/10+48)//char(imonth-imonth/10*10+48)
   date(7:8)=char(iday/10+48)//char(iday-iday/10*10+48)
   date(9:10)=char(ihour/10+48)//char(ihour-ihour/10*10+48)

   fillvalue = -99999.9
   do i=1,ny(ne)
   xlat(i)=(i-1)*1.0-89.5
   enddo
   do i=1,nx(ne)
   xlon(i)=(i-1)*1.0+0.5
   enddo
   do i=1,nzz
   sigma(i)=i
   enddo
   do i=1,isize
   a(i)=i
   enddo
   do i=1,ndustcom
   b(i)=i
   enddo
      
   if(ifdustcom.eq.1) then
   allocate(varid_x2a(25))
   mm=25
   else
   allocate(varid_x2a(12))
   mm=12  
   endif
   allocate(chem2d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1))
   allocate(chem3d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1,nz(ne)))   
!   allocate(chem5d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1,nz(ne),isize,ndustcom))
   
   fname='geatm.dust.'//date//'.nc'
   rcode = pio_createfile(pio_subsystem, file, pio_iotype_netcdf, trim(adjustl(fname)), PIO_CLOBBER)
   
   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)
   allocate(ldof2d(m)) 
   m=0
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof2d(m)=(j-1)*nx(ne)+i
   end do
   end do

   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)*nz(ne)
   allocate(ldof3d(m)) 
   m=0
   do k=1,nz(ne)
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof3d(m)=(j-1)*nx(ne)+i+(k-1)*nx(ne)*ny(ne)
   end do
   end do
   end do

!   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)*isize*nz(ne)*ndustcom
!   allocate(ldof5d(m))
!   m=0
!   do ia=1,ndustcom
!   do is=1,isize
!   do k=1,nz(ne)
!   do j=sy(ne),ey(ne)
!   do i=sx(ne),ex(ne)
!   m=m+1
!   ldof5d(m)=i+(j-1)*nx(ne)+(k-1)*nx(ne)*ny(ne)+(is-1)*nx(ne)*ny(ne)*nz(ne)+&
!             (ia-1)*nx(ne)*ny(ne)*nz(ne)*isize
!   end do
!   end do
!   end do
!   end do
!   end do

   dim2d=(/nx(ne),ny(ne)/)
   dim3d=(/nx(ne),ny(ne),nz(ne)/)
   dim5d=(/nx(ne),ny(ne),nz(ne),isize,ndustcom/)
   call pio_initdecomp(pio_subsystem, pio_real, dim2d, ldof2d, iodesc2d)
   call pio_initdecomp(pio_subsystem, pio_real, dim3d, ldof3d, iodesc3d)
!   call pio_initdecomp(pio_subsystem, pio_real, dim5d, ldof5d, iodesc5d)
   
   rcode = pio_def_dim(File,'lon',nx(ne),dim5d(1))
   rcode = pio_def_dim(File,'lat',ny(ne),dim5d(2))
   dim2d(1:2)=dim5d(1:2)
   rcode = pio_def_dim(File,'lev',nzz,dim5d(3))  
   dim3d(1:3)=dim5d(1:3)
   rcode = pio_def_dim(File,'isize',isize,dim5d(4))
   rcode = pio_def_dim(File,'ndustcom',ndustcom,dim5d(5))   

   rcode = pio_def_var(File,'lon',PIO_real,dim5d(1:1),varid_x2a(1))
   rcode = pio_put_att (File, varid_x2a(1), 'long_name', 'longitude')
   rcode = pio_put_att (File, varid_x2a(1), 'units', 'degrees_east')
   rcode = pio_def_var(File,'lat',PIO_REAL,dim5d(2:2),varid_x2a(2))
   rcode = pio_put_att (File, varid_x2a(2), 'long_name', 'latitude')
   rcode = pio_put_att (File, varid_x2a(2), 'units', 'degrees_north')
   rcode = pio_def_var(File,'lev',PIO_REAL,dim5d(3:3),varid_x2a(3))
   rcode = pio_def_var(File,'isize',PIO_REAL,dim5d(4:4),varid_x2a(4))
   rcode = pio_def_var(File,'ndustcom',PIO_REAL,dim5d(5:5),varid_x2a(5))   
   rcode = pio_def_var(File,'DUSTEXT',PIO_REAL,dim3d,varid_x2a(6))      
   rcode = pio_def_var(File,'DUSTAOD',PIO_REAL,dim2d,varid_x2a(7))      
   rcode = pio_def_var(File,'DUSTEMISS',PIO_REAL,dim2d,varid_x2a(8))    
   rcode = pio_def_var(File,'DUSTDRY',PIO_REAL,dim2d,varid_x2a(9))     
   rcode = pio_def_var(File,'DUSTGRAV',PIO_REAL,dim2d,varid_x2a(10))    
   rcode = pio_def_var(File,'DUSTWET',PIO_REAL,dim2d,varid_x2a(11))
   rcode = pio_def_var(File,'EMITFACT',PIO_REAL,dim2d,varid_x2a(12))    
   mm=12
   if(ifdustcom.eq.1) then      
 !  rcode = pio_def_var(File,'DUSTDRYSO4',PIO_REAL,dim2d,varid_x2a(13))
 !  rcode = pio_def_var(File,'DUSTDRYNO3',PIO_REAL,dim2d,varid_x2a(14))
 !  rcode = pio_def_var(File,'DUSTDRYFeII',PIO_REAL,dim2d,varid_x2a(15))
 !  rcode = pio_def_var(File,'DUSTDRYFeIII',PIO_REAL,dim2d,varid_x2a(16))
 !  rcode = pio_def_var(File,'DUSTWETSO4',PIO_REAL,dim2d,varid_x2a(17))
 !  rcode = pio_def_var(File,'DUSTWETNO3',PIO_REAL,dim2d,varid_x2a(18))
 !  rcode = pio_def_var(File,'DUSTWETFeII',PIO_REAL,dim2d,varid_x2a(19))
 !  rcode = pio_def_var(File,'DUSTWETFeIII',PIO_REAL,dim2d,varid_x2a(20))
 !  rcode = pio_def_var(File,'DUSTGRAVSO4',PIO_REAL,dim2d,varid_x2a(21))
 !  rcode = pio_def_var(File,'DUSTGRAVNO3',PIO_REAL,dim2d,varid_x2a(22))
 !  rcode = pio_def_var(File,'DUSTGRAVFeII',PIO_REAL,dim2d,varid_x2a(23))
 !  rcode = pio_def_var(File,'DUSTGRAVFeIII',PIO_REAL,dim2d,varid_x2a(24))
 !  rcode = pio_def_var(File,'dustcomp',PIO_REAL,dim5d,varid_x2a(25))
   mm=25
   end if       
   
   do i=1,mm
   rcode = pio_put_att(File,varid_x2a(i),"_fillvalue",fillvalue)
   end do
   rcode = pio_enddef(File) 
   
   rcode = pio_put_var(File, varid_x2a(1), xlon)
   rcode = pio_put_var(File, varid_x2a(2), xlat)
   rcode = pio_put_var(File, varid_x2a(3), sigma)
   rcode = pio_put_var(File, varid_x2a(4), a)
   rcode = pio_put_var(File, varid_x2a(5), b)   
          
   ! dustext
   mem3d=1         
   do k=1,nz(ne)        
          do i= sx(ne),ex(ne)
          do j= sy(ne),ey(ne)
            km=mem3d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=DUSTEXT(km)                           
          end do
          end do
          mem3d=mem3d+(nx(ne)+2)*(ny(ne)+2)
   end do  
   call pio_write_darray(File, varid_x2a(6), iodesc3d,chem3d, rcode)    
   ! dustaod       
          i0=ip2mem(ne)
          do i= sx(ne),ex(ne)
          do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=dustaod(km) 
          end do
          end do   
   call pio_write_darray(File, varid_x2a(7), iodesc2d, chem2d, rcode)
   ! dustemiss
          i0=ip2mem(ne)
          do i= sx(ne),ex(ne)
          do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=dustemiss(km)
          end do
          end do   
   call pio_write_darray(File, varid_x2a(8), iodesc2d, chem2d, rcode)
   ! dustdry
        i0=ip2mem(ne)
        do i= sx(ne),ex(ne)
        do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=dustdry(km)                           
        end do
        end do
   call pio_write_darray(File, varid_x2a(9), iodesc3d, chem3d, rcode)
   ! dustgrav
        i0=ip2mem(ne)
        do i= sx(ne),ex(ne)
        do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=dustgrav(km)                           
        end do
        end do
   call pio_write_darray(File, varid_x2a(10), iodesc2d, chem2d, rcode)   
   ! dustwet
        i0=ip2mem(ne)
        do i= sx(ne),ex(ne)
        do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=dustwet(km)                           
        end do
        end do
   call pio_write_darray(File, varid_x2a(11), iodesc2d, chem2d, rcode)
   ! emitfact   
        i0=ip2mem(ne)
        do i= sx(ne),ex(ne)
        do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=emitfact(km)                           
        end do
        end do        
   call pio_write_darray(File, varid_x2a(12), iodesc2d, chem2d, rcode)
   if(ifdustcom.eq.1) then
   ! dustdryso4   
        i0=ip2mem(ne)
        do i= sx(ne),ex(ne)
        do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
!            chem2d(i-sx(ne)+1,j-sy(ne)+1)=dustdryso4(km)                           
        end do
        end do        
!   call pio_write_darray(File, varid_x2a(13), iodesc2d, chem2d, rcode)
   ! dustdryno3   
        i0=ip2mem(ne)
        do i= sx(ne),ex(ne)
        do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
!            chem2d(i-sx(ne)+1,j-sy(ne)+1)=dustdryno3(km)                           
        end do
        end do        
!   call pio_write_darray(File, varid_x2a(14), iodesc2d, chem2d, rcode)
   ! dustdryfeii   
        i0=ip2mem(ne)
        do i= sx(ne),ex(ne)
        do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
!            chem2d(i-sx(ne)+1,j-sy(ne)+1)=dustdryfeii(km)                           
        end do
        end do        
!   call pio_write_darray(File, varid_x2a(15), iodesc2d, chem2d, rcode)
   ! dustdryfeiii
        i0=ip2mem(ne)
        do i= sx(ne),ex(ne)
        do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
!            chem2d(i-sx(ne)+1,j-sy(ne)+1)=dustdryfeiii(km)                           
        end do
        end do        
!   call pio_write_darray(File, varid_x2a(16), iodesc2d, chem2d, rcode)
   ! dustwetso4   
        i0=ip2mem(ne)
        do i= sx(ne),ex(ne)
        do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
!            chem2d(i-sx(ne)+1,j-sy(ne)+1)=dustdryso4(km)                           
        end do
        end do        
!   call pio_write_darray(File, varid_x2a(17), iodesc2d, chem2d, rcode)
   ! dustwetno3   
        i0=ip2mem(ne)
        do i= sx(ne),ex(ne)
        do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
!            chem2d(i-sx(ne)+1,j-sy(ne)+1)=dustwetno3(km)                           
        end do
        end do        
!   call pio_write_darray(File, varid_x2a(18), iodesc2d, chem2d, rcode)
   ! dustwetfeii   
        i0=ip2mem(ne)
        do i= sx(ne),ex(ne)
        do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
!            chem2d(i-sx(ne)+1,j-sy(ne)+1)=dustwetfeii(km)                           
        end do
        end do        
!   call pio_write_darray(File, varid_x2a(19), iodesc2d, chem2d, rcode)
   ! dustwetfeiii
        i0=ip2mem(ne)
        do i= sx(ne),ex(ne)
        do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
!            chem2d(i-sx(ne)+1,j-sy(ne)+1)=dustwetfeiii(km)                           
        end do
        end do        
!   call pio_write_darray(File, varid_x2a(20), iodesc2d, chem2d, rcode)
   ! dustgravso4   
        i0=ip2mem(ne)
        do i= sx(ne),ex(ne)
        do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
!            chem2d(i-sx(ne)+1,j-sy(ne)+1)=dustgravso4(km)                           
        end do
        end do        
!   call pio_write_darray(File, varid_x2a(21), iodesc2d, chem2d, rcode)
   ! dustgravno3   
        i0=ip2mem(ne)
        do i= sx(ne),ex(ne)
        do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
!            chem2d(i-sx(ne)+1,j-sy(ne)+1)=dustgravno3(km)                           
        end do
        end do        
!   call pio_write_darray(File, varid_x2a(22), iodesc2d, chem2d, rcode)
   ! dustgravfeii   
        i0=ip2mem(ne)
        do i= sx(ne),ex(ne)
        do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
!            chem2d(i-sx(ne)+1,j-sy(ne)+1)=dustgravfeii(km)                           
        end do
        end do        
!   call pio_write_darray(File, varid_x2a(23), iodesc2d, chem2d, rcode)
   ! dustgravfeiii
        i0=ip2mem(ne)
        do i= sx(ne),ex(ne)
        do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
!            chem2d(i-sx(ne)+1,j-sy(ne)+1)=dustgravfeiii(km)                           
        end do
        end do        
!   call pio_write_darray(File, varid_x2a(24), iodesc2d, chem2d, rcode)      
   ! dustcomp
   mem5d=1        ! number of memory for 5D variables  
   do ia=1,ndustcom
   do is=1,isize   
   do k=1,nzz        
        do i= sx(ne),ex(ne)
        do j= sy(ne),ey(ne)
            km=mem5d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
!            chem5d(i-sx(ne)+1,j-sy(ne)+1,k,is,ia)=dustcomp(km)                           
        end do
        end do
        mem5d=mem5d+(nx(ne)+2)*(ny(ne)+2)
   end do
   end do
   end do
!   call pio_write_darray(File, varid_x2a(25), iodesc5d, chem5d, rcode)
  end if   
   
   call pio_freedecomp(File,iodesc2d)
   call pio_freedecomp(File,iodesc3d)
!   call pio_freedecomp(File,iodesc5d)
   call pio_closefile(file)     

   deallocate(ldof2d)
   deallocate(ldof3d)
!   deallocate(ldof5d)
   deallocate(chem2d)
   deallocate(chem3d)   
!   deallocate(chem5d)   
   deallocate(varid_x2a)
   end subroutine output_dust
!--------------------------------------------------------------------------------

subroutine output_dry(iyear, imonth, iday, ihour, ne)
   use geatm_vartype, only:sx,ex,sy,ey,&
       nx,ny,nz,nzz,igas,&
       ip2mem,ip3mem,&
       DryVelGas, DRYDEP2
      
   integer :: iyear, imonth, iday, ihour, ne
   
   type(file_desc_t) :: file
   type(io_desc_t) :: iodesc2d,iodesc3d,iodesc4d,iodesc5d
   type(var_desc_t), pointer :: varid_x2a(:)
   integer :: i, j, k, ig, is, ia, dim2d(2), dim3d(3), dim4d(4),dim5d(5)
   integer :: mem3d, mem2d
   real :: xlat(ny(ne)), xlon(nx(ne)), sigma(nzz), a(igas), fill_value
   real, dimension(:,:,:),allocatable :: chem3d
   integer,pointer :: ldof3d(:)
   integer         :: m, mm        ! numbr of variables
   integer         :: rcode        ! return error code
   character*10 date
   character*30 fname
   
   date(1:4)=char(iyear/1000+48)//char(iyear/100-iyear/1000*10+48)//char(iyear/10-iyear/100*10+48)//&
             char(iyear-iyear/10*10+48)
   date(5:6)=char(imonth/10+48)//char(imonth-imonth/10*10+48)
   date(7:8)=char(iday/10+48)//char(iday-iday/10*10+48)
   date(9:10)=char(ihour/10+48)//char(ihour-ihour/10*10+48)
   fillvalue = -99999.9

   do i=1,ny(ne)
   xlat(i)=(i-1)*1.0-89.5
   enddo
   do i=1,nx(ne)
   xlon(i)=(i-1)*1.0+0.5
   enddo
   do i=1,igas
   a(i)=i
   enddo
      
   allocate(varid_x2a(5))  
   allocate(chem3d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1,igas))
      
   fname='geatm.dry.'//date//'.nc'
   rcode = pio_createfile(pio_subsystem, file, pio_iotype_netcdf, trim(adjustl(fname)), PIO_CLOBBER)
      
   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)*igas
   allocate(ldof3d(m)) 
   m=0
   do k=1,igas
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof3d(m)=(j-1)*nx(ne)+i+(k-1)*nx(ne)*ny(ne)
   end do
   end do
   end do

   dim3d=(/nx(ne),ny(ne),igas/)
   call pio_initdecomp(pio_subsystem, pio_real, dim3d, ldof3d, iodesc3d) 
 
   rcode = pio_def_dim(File,'lon',nx(ne),dim3d(1))
   rcode = pio_def_dim(File,'lat',ny(ne),dim3d(2))
   rcode = pio_def_dim(File,'igas',igas,dim3d(3))   

   rcode = pio_def_var(File,'lon',PIO_real,dim3d(1:1),varid_x2a(1))
   rcode = pio_put_att (File, varid_x2a(1), 'long_name', 'longitude')
   rcode = pio_put_att (File, varid_x2a(1), 'units', 'degrees_east')
   rcode = pio_def_var(File,'lat',PIO_REAL,dim3d(2:2),varid_x2a(2))
   rcode = pio_put_att (File, varid_x2a(2), 'long_name', 'latitude')
   rcode = pio_put_att (File, varid_x2a(2), 'units', 'degrees_north')
   rcode = pio_def_var(File,'igas',PIO_REAL,dim3d(3:3),varid_x2a(3))
   rcode = pio_def_var(File,'DryVelGas',PIO_REAL,dim3d,varid_x2a(4))    
   rcode = pio_def_var(File,'DRYDEP2',PIO_REAL,dim3d,varid_x2a(5))      
        
   do i=1,5
   rcode = pio_put_att(File,varid_x2a(i),"_fillvalue",fillvalue)
   end do
   rcode = pio_enddef(File) 
   
   rcode = pio_put_var(File, varid_x2a(1), xlon)
   rcode = pio_put_var(File, varid_x2a(2), xlat)
   rcode = pio_put_var(File, varid_x2a(3), a)   
          
   ! DryVelGas
   mem3d=1         
   do k=1,igas        
          do j= sy(ne),ey(ne)
          do i= sx(ne),ex(ne)
            km=mem3d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=DryVelGas(km)                           
          end do
          end do
          mem3d=mem3d+(nx(ne)+2)*(ny(ne)+2)
   end do  
   call pio_write_darray(File, varid_x2a(4), iodesc3d,chem3d, rcode)    
   ! Drydep2
   mem3d=1         
   do k=1,igas        
          do j= sy(ne),ey(ne)
          do i= sx(ne),ex(ne)
            km=mem3d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=Drydep2(km)                           
          end do
          end do
          mem3d=mem3d+(nx(ne)+2)*(ny(ne)+2)
   end do  
   call pio_write_darray(File, varid_x2a(5), iodesc3d,chem3d, rcode)
        
   call pio_freedecomp(File,iodesc3d)
   call pio_closefile(file)     
   
   deallocate(ldof3d)
   deallocate(chem3d)
   deallocate(varid_x2a)
   end subroutine output_dry
   
!-------------------------------------------------------------------------------        

subroutine output_wet(iyear, imonth, iday, ihour, ne)
   use geatm_vartype, only:sx,ex,sy,ey,&
       nx,ny,nz,nzz,igas,&
       ip2mem,ip3mem,&
       WETDEP, WETDEP2, RAINCON, RAINNON
      
   integer :: iyear, imonth, iday, ihour, ne
   
   type(file_desc_t) :: file
   type(io_desc_t) :: iodesc2d,iodesc3d,iodesc4d,iodesc5d
   type(var_desc_t), pointer :: varid_x2a(:)
   integer :: i, j, k, ig, is, ia, dim2d(2), dim3d(3), dim4d(4), dim5d(5)
   integer :: mem3d, mem2d
   real :: xlat(ny(ne)), xlon(nx(ne)), sigma(nzz), a(igas), fill_value
   real, dimension(:,:),allocatable :: chem2d
   real, dimension(:,:,:),allocatable :: chem3d
   integer,pointer :: ldof3d(:),ldof2d(:)
   integer         :: m, mm        ! numbr of variables
   integer         :: rcode        ! return error code
   character*10 date
   character*30 fname
  
   date(1:4)=char(iyear/1000+48)//char(iyear/100-iyear/1000*10+48)//char(iyear/10-iyear/100*10+48)//&
             char(iyear-iyear/10*10+48)
   date(5:6)=char(imonth/10+48)//char(imonth-imonth/10*10+48)
   date(7:8)=char(iday/10+48)//char(iday-iday/10*10+48)
   date(9:10)=char(ihour/10+48)//char(ihour-ihour/10*10+48)
   fillvalue = -99999.9

   do i=1,ny(ne)
   xlat(i)=(i-1)*1.0-89.5
   enddo
   do i=1,nx(ne)
   xlon(i)=(i-1)*1.0+0.5
   enddo
   do i=1,igas
   a(i)=i
   enddo
      
   allocate(varid_x2a(7))
   allocate(chem2d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1))
   allocate(chem3d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1,igas))
      
   fname='geatm.wet.'//date//'.nc'
   rcode = pio_createfile(pio_subsystem, file, pio_iotype_netcdf, trim(adjustl(fname)), PIO_CLOBBER)
   
   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)
   allocate(ldof2d(m)) 
   m=0
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof2d(m)=(j-1)*nx(ne)+i
   end do
   end do

   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)*igas
   allocate(ldof3d(m)) 
   m=0
   do k=1,igas
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof3d(m)=(j-1)*nx(ne)+i+(k-1)*nx(ne)*ny(ne)
   end do
   end do
   end do

   dim2d=(/nx(ne),ny(ne)/)
   dim3d=(/nx(ne),ny(ne),igas/)
   call pio_initdecomp(pio_subsystem, pio_real, dim2d, ldof2d, iodesc2d)   
   call pio_initdecomp(pio_subsystem, pio_real, dim3d, ldof3d, iodesc3d)   
 
   rcode = pio_def_dim(File,'lon',nx(ne),dim3d(1))
   rcode = pio_def_dim(File,'lat',ny(ne),dim3d(2))
   dim2d(1:2)=dim3d(1:2)
   rcode = pio_def_dim(File,'igas',igas,dim3d(3))   

   rcode = pio_def_var(File,'lon',PIO_real,dim3d(1:1),varid_x2a(1))
   rcode = pio_put_att (File, varid_x2a(1), 'long_name', 'longitude')
   rcode = pio_put_att (File, varid_x2a(1), 'units', 'degrees_east')
   rcode = pio_def_var(File,'lat',PIO_REAL,dim3d(2:2),varid_x2a(2))
   rcode = pio_put_att (File, varid_x2a(2), 'long_name', 'latitude')
   rcode = pio_put_att (File, varid_x2a(2), 'units', 'degrees_north')
   rcode = pio_def_var(File,'igas',PIO_REAL,dim3d(3:3),varid_x2a(3))  
   rcode = pio_def_var(File,'RAINNON',PIO_REAL,dim2d,varid_x2a(4))
   rcode = pio_def_var(File,'RAINCON',PIO_REAL,dim2d,varid_x2a(5))
   rcode = pio_def_var(File,'WETDEP',PIO_REAL,dim3d,varid_x2a(6))       
   rcode = pio_def_var(File,'WETDEP2',PIO_REAL,dim3d,varid_x2a(7))      
        
   do i=1,7
   rcode = pio_put_att(File,varid_x2a(i),"_fillvalue",fillvalue)
   end do
   rcode = pio_enddef(File) 
   
   rcode = pio_put_var(File, varid_x2a(1), xlon)
   rcode = pio_put_var(File, varid_x2a(2), xlat)
   rcode = pio_put_var(File, varid_x2a(3), a)   

   ! rainnon
          i0=ip2mem(ne)
          do i= sx(ne),ex(ne)
          do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=rainnon(km)                           
          end do
          end do
   call pio_write_darray(File, varid_x2a(4), iodesc2d, chem2d, rcode)
   ! raincon
          i0=ip2mem(ne)
          do i= sx(ne),ex(ne)
          do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=raincon(km)                           
          end do
          end do
   call pio_write_darray(File, varid_x2a(5), iodesc2d, chem2d, rcode)   
   ! WEPDEP
   mem3d=1         
   do k=1,igas        
          do i= sx(ne),ex(ne)
          do j= sy(ne),ey(ne)
            km=mem3d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=WETDEP(km)                           
          end do
          end do
          mem3d=mem3d+(nx(ne)+2)*(ny(ne)+2)
   end do  
   call pio_write_darray(File, varid_x2a(6), iodesc3d, chem3d, rcode)    
   ! WEPDEP2
   mem3d=1
   do k=1,igas        
          do i= sx(ne),ex(ne)
          do j= sy(ne),ey(ne)
            km=mem3d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=wetdep2(km)                           
          end do
          end do
          mem3d=mem3d+(nx(ne)+2)*(ny(ne)+2)
   end do  
   call pio_write_darray(File, varid_x2a(7), iodesc3d, chem3d, rcode)
        
   call pio_freedecomp(File,iodesc2d)
   call pio_freedecomp(File,iodesc3d)
   call pio_closefile(file)     
   
   deallocate(ldof2d)
   deallocate(ldof3d)
   deallocate(chem2d)
   deallocate(chem3d)
   deallocate(varid_x2a)
   end subroutine output_wet
   
!-------------------------------------------------------------------------------

subroutine output_sea(iyear, imonth, iday, ihour, ne)
   use geatm_vartype, only:sx,ex,sy,ey,&
       nx,ny,nz,nzz,isize,ifseacom,nseacom,&
       ip2mem,ip3mem,&
       seacomp, seaemiss

   integer :: iyear, imonth, iday, ihour, ne
   
   type(file_desc_t) :: file
   type(io_desc_t) :: iodesc2d,iodesc3d,iodesc4d,iodesc5d
   type(var_desc_t), pointer :: varid_x2a(:)
   integer :: i, j, k, ig, is, ia, dim2d(2), dim3d(3), dim4d(4),dim5d(5)
   integer :: mem5d
   real :: xlat(ny(ne)), xlon(nx(ne)), sigma(nzz), a(isize), b(nseacom), fill_value
   real, dimension(:,:),allocatable :: chem2d
!   real, dimension(:,:,:,:,:),allocatable :: chem5d
   integer         :: m, mm        ! numbr of variables
   integer         :: rcode        ! return error code
   integer,pointer :: ldof2d(:)
!   integer,pointer :: ldof5d(:)
   character*10 date
   character*30 fname
   
   date(1:4)=char(iyear/1000+48)//char(iyear/100-iyear/1000*10+48)//char(iyear/10-iyear/100*10+48)//&
             char(iyear-iyear/10*10+48)
   date(5:6)=char(imonth/10+48)//char(imonth-imonth/10*10+48)
   date(7:8)=char(iday/10+48)//char(iday-iday/10*10+48)
   date(9:10)=char(ihour/10+48)//char(ihour-ihour/10*10+48)
   fillvalue = -99999.9

   do i=1,ny(ne)
   xlat(i)=(i-1)*1.0-89.5
   enddo
   do i=1,nx(ne)
   xlon(i)=(i-1)*1.0+0.5
   enddo
   do i=1,nz(ne)
   sigma(i)=i
   enddo
   do i=1,isize
   a(i)=i
   enddo
   do i=1,nseacom
   b(i)=i
   enddo
   
   if(ifseacom.eq.1) then   
   allocate(varid_x2a(7))
   mm=7
   else
   allocate(varid_x2a(6))
   mm=6
   end if
   
   allocate(chem2d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1))      
!   allocate(chem5d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1,nz(ne),isize,nseacom))
   
   fname='geatm.sea.'//date//'.nc'
   rcode = pio_createfile(pio_subsystem, file, pio_iotype_netcdf, trim(adjustl(fname)), PIO_CLOBBER)
   
   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)
   allocate(ldof2d(m)) 
   m=0
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof2d(m)=(j-1)*nx(ne)+i
   end do
   end do

!   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)*isize*nz(ne)*nseacom
!   allocate(ldof5d(m))
!   m=0
!   do ia=1,nseacom
!   do is=1,isize
!   do k=1,nz(ne)
!   do j=sy(ne),ey(ne)
!   do i=sx(ne),ex(ne)
!   m=m+1
!   ldof5d(m)=i+(j-1)*nx(ne)+(k-1)*nx(ne)*ny(ne)+(is-1)*nx(ne)*ny(ne)*nz(ne)+&
!             (ia-1)*nx(ne)*ny(ne)*nz(ne)*isize
!   end do
!   end do
!   end do
!   end do
!   end do

   dim2d=(/nx(ne),ny(ne)/)
   dim5d=(/nx(ne),ny(ne),nz(ne),isize,nseacom/)
   call pio_initdecomp(pio_subsystem, pio_real, dim2d, ldof2d, iodesc2d)    
!   call pio_initdecomp(pio_subsystem, pio_real, dim5d, ldof5d, iodesc5d)
   
   rcode = pio_def_dim(File,'lon',nx(ne),dim5d(1))
   rcode = pio_def_dim(File,'lat',ny(ne),dim5d(2))
   dim2d(1:2)=dim5d(1:2)
   rcode = pio_def_dim(File,'lev',nzz,dim5d(3))  
   rcode = pio_def_dim(File,'isize',isize,dim5d(4))
   rcode = pio_def_dim(File,'nseacom',nseacom,dim5d(5))   

   rcode = pio_def_var(File,'lon',PIO_real,dim5d(1:1),varid_x2a(1))
   rcode = pio_put_att (File, varid_x2a(1), 'long_name', 'longitude')
   rcode = pio_put_att (File, varid_x2a(1), 'units', 'degrees_east')
   rcode = pio_def_var(File,'lat',PIO_REAL,dim5d(2:2),varid_x2a(2))
   rcode = pio_put_att (File, varid_x2a(2), 'long_name', 'latitude')
   rcode = pio_put_att (File, varid_x2a(2), 'units', 'degrees_north')
   rcode = pio_def_var(File,'lev',PIO_REAL,dim5d(3:3),varid_x2a(3))
   rcode = pio_def_var(File,'isize',PIO_REAL,dim5d(4:4),varid_x2a(4))
   rcode = pio_def_var(File,'nseacom',PIO_REAL,dim5d(5:5),varid_x2a(5))    
   rcode = pio_def_var(File,'seaemiss',PIO_REAL,dim2d,varid_x2a(6))
   mm=6
   if(ifseacom.eq.1) then      
!   mm=7
!   rcode = pio_def_var(File,'seacomp',PIO_REAL,dim5d,varid_x2a(7))     
   end if
        
   do i=1,mm
   rcode = pio_put_att(File,varid_x2a(i),"_fillvalue",fillvalue)
   end do
   rcode = pio_enddef(File) 
   
   rcode = pio_put_var(File, varid_x2a(1), xlon)
   rcode = pio_put_var(File, varid_x2a(2), xlat)
   rcode = pio_put_var(File, varid_x2a(3), sigma)
   rcode = pio_put_var(File, varid_x2a(4), a)
   rcode = pio_put_var(File, varid_x2a(5), b)   
          
   ! seaemiss              
          i0=ip2mem(ne)
          do i= sx(ne),ex(ne)
          do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=seaemiss(km)                        
          end do
          end do   
   call pio_write_darray(File, varid_x2a(6), iodesc2d, chem2d, rcode)   
   if(ifseacom.eq.1) then
   ! seacomp
   mem5d=1        ! number of memory for 5D variables  
   do ia=1,nseacom
   do is=1,isize   
   do k=1,nzz        
        do i= sx(ne),ex(ne)
        do j= sy(ne),ey(ne)
            km=mem5d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
!            chem5d(i-sx(ne)+1,j-sy(ne)+1,k,is,ia)=seacomp(km)                           
        end do
        end do
        mem5d=mem5d+(nx(ne)+2)*(ny(ne)+2)
   end do
   end do
   end do
!   call pio_write_darray(File, varid_x2a(7), iodesc5d, chem5d, rcode)
  end if
        
   call pio_freedecomp(File,iodesc2d)
!   call pio_freedecomp(File,iodesc5d)
   call pio_closefile(file)     
   
   deallocate(ldof2d)
!   deallocate(ldof5d)
   deallocate(chem2d)
!   deallocate(chem5d)
   deallocate(varid_x2a)
   end subroutine output_sea
   
!-------------------------------------------------------------------------------

subroutine output_mark(iyear, imonth, iday, ihour, ne, restart)
   use geatm_vartype, only:sx,ex,sy,ey,&
       nx,ny,nz,nzz,isize,idmset,ismmax,&
       ip2mem,ip3mem,&
       sourcemark
      
   integer :: iyear, imonth, iday, ihour, ne
   logical :: restart
   
   type(file_desc_t) :: file
   type(io_desc_t) :: iodesc2d,iodesc3d,iodesc4d,iodesc5d
   type(var_desc_t), pointer :: varid_x2a(:)
   integer :: i, j, k, ig, is, ia, dim2d(2), dim3d(3), dim4d(4),dim5d(5)
   integer :: mem5d
   real :: xlat(ny(ne)), xlon(nx(ne)), sigma(nzz), a(ismmax), b(idmset)
   real, dimension(:,:,:,:,:),allocatable :: chem5d
   integer         :: m, mm        ! numbr of variables
   integer         :: rcode        ! return error code
   integer,pointer :: ldof5d(:)
   character*10 date
   character*30 fname

   date(1:4)=char(iyear/1000+48)//char(iyear/100-iyear/1000*10+48)//char(iyear/10-iyear/100*10+48)//&
             char(iyear-iyear/10*10+48)
   date(5:6)=char(imonth/10+48)//char(imonth-imonth/10*10+48)
   date(7:8)=char(iday/10+48)//char(iday-iday/10*10+48)
   date(9:10)=char(ihour/10+48)//char(ihour-ihour/10*10+48)   
   fillvalue = -99999.9

   do i=1,ny(ne)
   xlat(i)=(i-1)*1.0-89.5
   enddo
   do i=1,nx(ne)
   xlon(i)=(i-1)*1.0+0.5
   enddo
   do i=1,nzz
   sigma(i)=i
   enddo
   do i=1,ismMax
   a(i)=i
   enddo
   do i=1,idmSet
   b(i)=i
   enddo
   
   allocate(varid_x2a(6))
      
   allocate(chem5d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1,nz(ne),ismmax,idmset))
   
   if(restart) then
   fname='geatm.restd.'//date//'.nc'
   else
   fname='geatm.mark.'//date//'.nc'
   endif
   rcode = pio_createfile(pio_subsystem, file, pio_iotype_netcdf, trim(adjustl(fname)), PIO_CLOBBER)   
   
   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)*ismmax*nz(ne)*idmset
   allocate(ldof5d(m))
   m=0
   do ia=1,idmset
   do is=1,ismmax
   do k=1,nz(ne)
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof5d(m)=(j-1)*nx(ne)+(k-1)*nx(ne)*ny(ne)+(is-1)*nx(ne)*ny(ne)*nz(ne)+&
             (ia-1)*nx(ne)*ny(ne)*nz(ne)*ismmax
   end do
   end do
   end do
   end do
   end do

   dim5d=(/nx(ne),ny(ne),nz(ne),ismmax,idmset/)
   call pio_initdecomp(pio_subsystem, pio_real, dim5d, ldof5d, iodesc5d)
   
   rcode = pio_def_dim(File,'lon',nx(ne),dim5d(1))
   rcode = pio_def_dim(File,'lat',ny(ne),dim5d(2))
   rcode = pio_def_dim(File,'lev',nzz,dim5d(3))  
   rcode = pio_def_dim(File,'ismmax',ismmax,dim5d(4))
   rcode = pio_def_dim(File,'idmset',idmset,dim5d(5))   

   rcode = pio_def_var(File,'lon',PIO_real,dim5d(1:1),varid_x2a(1))
   rcode = pio_put_att (File, varid_x2a(1), 'long_name', 'longitude')
   rcode = pio_put_att (File, varid_x2a(1), 'units', 'degrees_east')
   rcode = pio_def_var(File,'lat',PIO_REAL,dim5d(2:2),varid_x2a(2))
   rcode = pio_put_att (File, varid_x2a(2), 'long_name', 'latitude')
   rcode = pio_put_att (File, varid_x2a(2), 'units', 'degrees_north')
   rcode = pio_def_var(File,'lev',PIO_REAL,dim5d(3:3),varid_x2a(3))
   rcode = pio_def_var(File,'ismmax',PIO_REAL,dim5d(4:4),varid_x2a(4))
   rcode = pio_def_var(File,'idmset',PIO_REAL,dim5d(5:5),varid_x2a(5))   
   rcode = pio_def_var(File,'sourcemark',PIO_REAL,dim5d,varid_x2a(6))
        
   do i=1,6
   rcode = pio_put_att(File,varid_x2a(i),"_fillvalue",fillvalue)
   end do
   rcode = pio_enddef(File) 
   
   rcode = pio_put_var(File, varid_x2a(1), xlon)
   rcode = pio_put_var(File, varid_x2a(2), xlat)
   rcode = pio_put_var(File, varid_x2a(3), sigma)
   rcode = pio_put_var(File, varid_x2a(4), a)
   rcode = pio_put_var(File, varid_x2a(5), b)   
          
   ! seacomp
   mem5d=1        ! number of memory for 5D variables  
   do ia=1,idmset
   do is=1,ismmax  
   do k=1,nzz        
        do i= sx(ne),ex(ne)
        do j= sy(ne),ey(ne)
            km=mem5d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem5d(i-sx(ne)+1,j-sy(ne)+1,k,is,ia)=sourcemark(km)                           
        end do
        end do
        mem5d=mem5d+(nx(ne)+2)*(ny(ne)+2)
   end do
   end do
   end do
   call pio_write_darray(File, varid_x2a(6), iodesc5d, chem5d, rcode)
        
   call pio_freedecomp(File,iodesc5d)
   call pio_closefile(file)     
      
   deallocate(ldof5d)
   deallocate(chem5d)
   deallocate(varid_x2a)
   end subroutine output_mark
!-------------------------------------------------------------------------------

subroutine output_term(iyear, imonth, iday, ihour, ne)
   use geatm_vartype, only:sx,ex,sy,ey,&
       nx,ny,nz,nzz,isize,igas,iprocess,&
       ip2mem,ip3mem,ip5mem,IGOPos,PrintTermGas,&
       ipGasTermBal,GasTermBal
      
   integer :: iyear, imonth, iday, ihour, ne
   
   type(file_desc_t) :: file
   type(io_desc_t) :: iodesc2d,iodesc3d,iodesc4d,iodesc5d
   type(var_desc_t), pointer :: varid_x2a(:)
   integer :: i, j, k, i05, igo, is, ia, mlen, dim2d(2), dim3d(3), dim4d(4),dim5d(5)
   integer :: mem5d
   real :: xlat(ny(ne)), xlon(nx(ne)), sigma(nzz), a(iprocess), b(igas), fill_value
   real, dimension(:,:,:,:,:),allocatable :: chem5d
   integer         :: m,mm        ! numbr of variables
   integer         :: rcode        ! return error code
   integer,pointer :: ldof5d(:)
   character*10 date
   character*30 fname
  
   date(1:4)=char(iyear/1000+48)//char(iyear/100-iyear/1000*10+48)//char(iyear/10-iyear/100*10+48)//&
             char(iyear-iyear/10*10+48)
   date(5:6)=char(imonth/10+48)//char(imonth-imonth/10*10+48)
   date(7:8)=char(iday/10+48)//char(iday-iday/10*10+48)
   date(9:10)=char(ihour/10+48)//char(ihour-ihour/10*10+48)
   fillvalue = -99999.9

   do i=1,ny(ne)
   xlat(i)=(i-1)*1.0-89.5
   enddo
   do i=1,nx(ne)
   xlon(i)=(i-1)*1.0+0.5
   enddo
   do i=1,nzz
   sigma(i)=i
   enddo
   do i=1,iprocess
   a(i)=i
   enddo

   mlen=0
   do i=1,igas
   if(PrintTermGas(i).eq.1)then
   mlen=mlen+1
   b(i)=mlen
   end if
   enddo
   if (mlen.eq.0) then
   go to 9999
   end if
   
   allocate(varid_x2a(6))
   
   allocate(chem5d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1,nz(ne),iprocess,mlen))
   
   fname='geatm.term.'//date//'.nc'   
   rcode = pio_createfile(pio_subsystem, file, pio_iotype_netcdf, trim(adjustl(fname)), PIO_CLOBBER)   
       
   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)*iprocess*nz(ne)*m
   allocate(ldof5d(m))
   m=0
   do ia=1,mlen
   do is=1,iprocess
   do k=1,nzz
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof5d(m)=(j-1)*nx(ne)+(k-1)*nx(ne)*ny(ne)+(is-1)*nx(ne)*ny(ne)*nz(ne)+&
             (ia-1)*nx(ne)*ny(ne)*nz(ne)*isize
   end do
   end do
   end do
   end do
   end do

   dim5d=(/nx(ne),ny(ne),nz(ne),iprocess,mlen/)
   call pio_initdecomp(pio_subsystem, pio_real, dim5d, ldof5d, iodesc5d)
   
   rcode = pio_def_dim(File,'lon',nx(ne),dim5d(1))
   rcode = pio_def_dim(File,'lat',ny(ne),dim5d(2))
   rcode = pio_def_dim(File,'lev',nzz,dim5d(3))  
   rcode = pio_def_dim(File,'iprocess',iprocess,dim5d(4))
   rcode = pio_def_dim(File,'igas',mlen,dim5d(5))   

   rcode = pio_def_var(File,'lon',PIO_real,dim5d(1:1),varid_x2a(1))
   rcode = pio_put_att (File, varid_x2a(1), 'long_name', 'longitude')
   rcode = pio_put_att (File, varid_x2a(1), 'units', 'degrees_east')
   rcode = pio_def_var(File,'lat',PIO_REAL,dim5d(2:2),varid_x2a(2))
   rcode = pio_put_att (File, varid_x2a(2), 'long_name', 'latitude')
   rcode = pio_put_att (File, varid_x2a(2), 'units', 'degrees_north')
   rcode = pio_def_var(File,'lev',PIO_REAL,dim5d(3:3),varid_x2a(3))
   rcode = pio_def_var(File,'iprocess',PIO_REAL,dim5d(4:4),varid_x2a(4))
   rcode = pio_def_var(File,'igas',PIO_REAL,dim5d(5:5),varid_x2a(5))   
   rcode = pio_def_var(File,'GasTermBal',PIO_REAL,dim5d,varid_x2a(6))
        
   do i=1,6
   rcode = pio_put_att(File,varid_x2a(i),"_fillvalue",fillvalue)
   end do
   rcode = pio_enddef(File) 
   
   rcode = pio_put_var(File, varid_x2a(1), xlon)
   rcode = pio_put_var(File, varid_x2a(2), xlat)
   rcode = pio_put_var(File, varid_x2a(3), sigma)
   rcode = pio_put_var(File, varid_x2a(4), a)
   rcode = pio_put_var(File, varid_x2a(5), b(1:mlen))   
          
   ! GasTermBal
   do ia=1,igas
   if(PrintTermGas(ia)==1)then
   igo = IGOPos(ia) 
   do is=1,iprocess
   do k=1,nzz
        i05=ipGasTermBal(k,is,igo,ne)
        do i= sx(ne),ex(ne)
        do j= sy(ne),ey(ne)
            km=i05+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem5d(i-sx(ne)+1,j-sy(ne)+1,k,is,ia)=GasTermBal(km)                           
        end do
        end do   
   end do
   end do
   endif
   end do
   call pio_write_darray(File, varid_x2a(6), iodesc5d, chem5d, rcode)
        
   call pio_freedecomp(File,iodesc5d)
   call pio_closefile(file)     
  
   deallocate(ldof5d)    
   deallocate(chem5d)
   deallocate(varid_x2a)
   9999 continue
   end subroutine output_term

!--------------------------------------------------------------------------------------   
  subroutine init_gea_pio()
    use seq_io_mod,       only: seq_io_getiosys, seq_io_getiotype

    pio_subsystem => seq_io_getiosys('GEA')

  end subroutine init_gea_pio

  end module inout3
