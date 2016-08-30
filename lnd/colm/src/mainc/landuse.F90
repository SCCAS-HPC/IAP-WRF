
#include <define.h>

module landuse

   use precision
   use netcdf
   use spmd
   use colm_varctl, only: lcarbon_emission
   use colm_varMod, only: lon_points, lat_points

   implicit none

   save
   private

   character(len=255) :: fcarbon_emission = "null"

   integer :: nlon, nlat, ntime

   real(r8), pointer :: lon(:)
   real(r8), pointer :: lat(:)
   real(r8), pointer :: time(:)
   real(r8), pointer :: carbon_emission(:,:,:) ! g(C)/m2/s

   interface landuse_init
      module procedure landuse_init
   end interface

   interface landuse_exit
      module procedure landuse_exit
   end interface

   interface landuse_cflux
      module procedure landuse_cflux
   end interface

   public fcarbon_emission, carbon_emission
   public landuse_init, landuse_exit, landuse_cflux

contains

   subroutine landuse_init

      integer fid, vid

      if (.not. lcarbon_emission) return

      if (p_master) then
         print *, 'open land use carbon emission file: '//trim(fcarbon_emission)

         call sanity(nf90_open(path=trim(fcarbon_emission),mode=nf90_nowrite,ncid=fid))

         call sanity(nf90_inq_dimid(fid,'lon',vid))
         call sanity(nf90_inquire_dimension(fid,vid,len=nlon))
         call sanity(nf90_inq_dimid(fid,'lat',vid))
         call sanity(nf90_inquire_dimension(fid,vid,len=nlat))
         call sanity(nf90_inq_dimid(fid,'time',vid))
         call sanity(nf90_inquire_dimension(fid,vid,len=ntime))

         if (nlon.ne.lon_points .or. nlat.ne.lat_points) then
            write(6,*) 'fatal error on dimensions of land use emission file', nlon, nlat
            call sanity(nf90_close(fid))
            call abort
         end if
      end if

#ifdef SPMD
      call mpi_bcast (nlon ,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (nlat ,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (ntime,1,mpi_integer,0,p_comm,p_err)
#endif
      
      allocate (lon(nlon))
      allocate (lat(nlat))
      allocate (time(ntime))
      allocate (carbon_emission(nlon,nlat,ntime))

      if (p_master) then
         call sanity(nf90_inq_varid(fid,'lon',vid))
         call sanity(nf90_get_var(fid,vid,lon(:)))

         call sanity(nf90_inq_varid(fid,'lat',vid))
         call sanity(nf90_get_var(fid,vid,lat(:)))

         call sanity(nf90_inq_varid(fid,'time',vid))
         call sanity(nf90_get_var(fid,vid,time(:)))

         call sanity(nf90_inq_varid(fid,'carbon_emission',vid))
         call sanity(nf90_get_var(fid,vid,carbon_emission(:,:,:)))

         call sanity(nf90_close(fid))
      end if

#ifdef SPMD
      call mpi_bcast (lon,size(lon),mpi_real8,0,p_comm,p_err)
      call mpi_bcast (lat,size(lat),mpi_real8,0,p_comm,p_err)
      call mpi_bcast (time,size(time),mpi_real8,0,p_comm,p_err)
      call mpi_bcast (carbon_emission,size(carbon_emission),mpi_real8,0,p_comm,p_err)
#endif

   end subroutine

   subroutine landuse_cflux(year,cflux)

      integer, intent(in) :: year
      real(r8), intent(out) :: cflux(lon_points,lat_points)

      integer i, itime

      cflux(:,:) = 0._r8

      if (lcarbon_emission) then
         itime = -9999

         if (year.lt.time(1)) then
            itime = 1
         else if (year.ge.time(ntime)) then
            itime = ntime
         else
            do i = 1, ntime-1
               if(year.ge.time(i) .and. year.lt.time(i+1)) then
                  itime = i
                  exit
               end if 
            end do
         end if

         if (itime.lt.0) then
            write(6,*) 'no available land emission data'
            call abort
         end if

         cflux(:,:) = carbon_emission(:,:,itime)
      end if

   end subroutine

   subroutine landuse_exit

      if (.not. lcarbon_emission) return

      deallocate (lon)
      deallocate (lat)
      deallocate (time)
      deallocate (carbon_emission)

   end subroutine

   subroutine sanity(ret)

      integer, intent(in) :: ret

      if (ret .ne. nf90_noerr) then
         write(6,*) trim(nf90_strerror(ret)); stop
      end if
   end subroutine

end module landuse
