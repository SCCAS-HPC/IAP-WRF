#include <define.h>

module spmd_decomp

   implicit none

   integer,  pointer :: gmask(:,:)            ! grid => process index mapping

   integer,  pointer :: pcmap(:)              ! local patch => local column mapping
   integer,  pointer :: pgmap(:)              ! local patch => local grid mapping
   integer,  pointer :: cgmap(:)              ! local column => local grid mapping

   integer,  pointer :: ppmap(:)              ! local patch => global patch mapping
   integer,  pointer :: ccmap(:)              ! local column => global column mapping
   integer,  pointer :: ggmap(:)              ! local grid  => global grid  mapping
   integer,  pointer :: gxmap(:)              ! local grid  => global X mapping
   integer,  pointer :: gymap(:)              ! local grid  => global Y mapping

   integer,  pointer :: ppmap_glob(:)         ! global ppmap
   integer,  pointer :: ccmap_glob(:)         ! global ccmap
   integer,  pointer :: ggmap_glob(:)         ! global ggmap
   integer,  pointer :: gxmap_glob(:)         ! global gxmap
   integer,  pointer :: gymap_glob(:)         ! global gymap

   integer,  pointer :: numgrid_proc(:)       ! grid number per process
   integer,  pointer :: numcolumn_proc(:)
   integer,  pointer :: numpatch_proc(:)      ! patch number per process

   interface task_decomp
      module procedure task_decomp
   end interface task_decomp

   interface spmd_var_dealloc
      module procedure spmd_var_dealloc
   end interface spmd_var_dealloc

contains

   subroutine task_decomp

      use precision
      use colm_varMod 
      use spmd

!local variables:

      integer, allocatable :: pmask(:)         ! patch -> p_iam mask
      integer, allocatable :: cmask(:)         ! column -> p_iam mask

      integer  i, j, k, n1, n2, p, g, pid, gxy1, gxy2

!routine:

      allocate (pmask          (numpatch_glob))
      allocate (cmask         (numcolumn_glob))
      allocate (gmask  (lon_points,lat_points))

      allocate (numgrid_proc    (0:p_nprocs-1))
      allocate (numcolumn_proc  (0:p_nprocs-1))
      allocate (numpatch_proc   (0:p_nprocs-1))

      pmask(:) = -1
      cmask(:) = -1
      gmask(:,:) = -1

      numgrid_proc(:)   = 0
      numcolumn_proc(:) = 0
      numpatch_proc(:)  = 0

      pid = 0
      numgrid_glob = 1
      numgrid_proc(pid) = 1
      gmask(ixy_column_glob(1),jxy_column_glob(1)) = pid

      gxy1 = ixy_column_glob(1) + (jxy_column_glob(1)-1)*lon_points

      do i = 1, numcolumn_glob
         gxy2 = ixy_column_glob(i) + (jxy_column_glob(i)-1)*lon_points

         if(i.gt.1 .and. gxy1.ne.gxy2) then
            pid = mod(pid+1,p_nprocs)
            numgrid_glob = numgrid_glob + 1
            numgrid_proc(pid) = numgrid_proc(pid)+1
            gmask(ixy_column_glob(i),jxy_column_glob(i)) = pid
         end if

         gxy1 = gxy2

         cmask(i) = pid
         numcolumn_proc(pid) = numcolumn_proc(pid)+1
      end do

      pid = 0
      gxy1 = ixy_patch_glob(1) + (jxy_patch_glob(1)-1)*lon_points

      do i = 1, numpatch_glob
         gxy2 = ixy_patch_glob(i) + (jxy_patch_glob(i)-1)*lon_points

         if(i.gt.1 .and. gxy1.ne.gxy2) then
            pid = mod(pid+1,p_nprocs)
         end if

         gxy1 = gxy2

         pmask(i) = pid
         numpatch_proc(pid) = numpatch_proc(pid)+1
      end do

      if (numgrid_glob .ne. sum(numgrid_proc(:))) write(6,*) 'global grids number error'
      if (numcolumn_glob .ne. sum(numcolumn_proc(:))) write(6,*) 'global columns number error'
      if (numpatch_glob .ne. sum(numpatch_proc(:))) write(6,*) 'global patches number error'

      if (numgrid_glob .lt. p_nprocs) then
         write(6,*) 'total land grids less than MPI processes number'
         call abort
      end if

      if (p_master) then
         write(6,*) 'number of grids per process:'
         write(6,*)  numgrid_proc
         write(6,*) 'number of columns per process:'
         write(6,*)  numcolumn_proc
         write(6,*) 'number of patches per process:'
         write(6,*)  numpatch_proc
      endif

      numgrid = numgrid_proc(p_iam)
      numcolumn = numcolumn_proc(p_iam)
      numpatch = numpatch_proc(p_iam)

    ! Building local patch & global patch mapping
    ! Building local patch & local grid mapping

      allocate (ppmap       (numpatch))
      allocate (pgmap       (numpatch))
      allocate (ppmap_glob  (numpatch_glob))

      n1 = 0
      n2 = 0
      gxy1 = 0

      do i = 1, numpatch_glob
         gxy2 = ixy_patch_glob(i) + (jxy_patch_glob(i)-1)*lon_points

         if (pmask(i)==p_iam) then
            n1 = n1+1
            ppmap(n1) = i

            if(gxy1 .ne. gxy2) then
               n2 = n2+1
            end if

            pgmap(n1) = n2
         end if 

         gxy1 = gxy2
      end do

#ifdef SPMD
      p_counts(:) = numpatch_proc(:)
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numpatch_proc(0:i-1))
      end do

      call mpi_gatherv (ppmap,size(ppmap),mpi_integer,ppmap_glob,&
                        p_counts,p_displs,mpi_integer,0,p_comm,p_err)
#else
      ppmap_glob(:) = ppmap(:)
#endif

    ! Building local column & global column mapping
    ! Building local column & local grid mapping

      allocate (ccmap       (numcolumn))
      allocate (cgmap       (numcolumn))
      allocate (ccmap_glob  (numcolumn_glob))

      n1 = 0
      n2 = 0
      gxy1 = 0

      do i = 1, numcolumn_glob
         gxy2 = ixy_column_glob(i) + (jxy_column_glob(i)-1)*lon_points

         if (cmask(i)==p_iam) then
            n1 = n1+1
            ccmap(n1) = i

            if(gxy1 .ne. gxy2) then
               n2 = n2+1
            end if

            cgmap(n1) = n2
         end if 

         gxy1 = gxy2
      end do

#ifdef SPMD
      p_counts(:) = numcolumn_proc(:)
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numcolumn_proc(0:i-1))
      end do

      call mpi_gatherv (ccmap,size(ccmap),mpi_integer,ccmap_glob,&
                        p_counts,p_displs,mpi_integer,0,p_comm,p_err)
#else
      ccmap_glob(:) = ccmap(:)
#endif

    ! Building local grid & global G-X-Y mapping

      allocate (ggmap        (numgrid))
      allocate (gxmap        (numgrid))
      allocate (gymap        (numgrid))

      allocate (ggmap_glob   (numgrid_glob))
      allocate (gxmap_glob   (numgrid_glob))
      allocate (gymap_glob   (numgrid_glob))

      n1 = 0
      n2 = 0

      do j = 1, lat_points
      do i = 1, lon_points
         if (gmask(i,j).ge.0) then
            n1 = n1+1
         end if

         if (gmask(i,j)==p_iam) then
            n2 = n2+1
            ggmap(n2) = n1
            gxmap(n2) = i
            gymap(n2) = j
         end if 
      end do
      end do

      if (n1.ne.numgrid_glob .or. n2.ne.numgrid) then
         write(6,*) 'failed on numgrid checking', n1, numgrid_glob, n2, numgrid
         call abort
      end if

#ifdef SPMD
      p_counts(:) = numgrid_proc(:)
      p_displs(0) = 0

      do i = 1, p_nprocs-1
         p_displs(i) = sum(numgrid_proc(0:i-1))
      end do

      call mpi_allgatherv (ggmap,size(ggmap),mpi_integer,ggmap_glob,&
                           p_counts,p_displs,mpi_integer,p_comm,p_err)

      call mpi_allgatherv (gxmap,size(gxmap),mpi_integer,gxmap_glob,&
                           p_counts,p_displs,mpi_integer,p_comm,p_err)

      call mpi_allgatherv (gymap,size(gymap),mpi_integer,gymap_glob,&
                           p_counts,p_displs,mpi_integer,p_comm,p_err)
#else
      ggmap_glob(:) = ggmap(:)
      gxmap_glob(:) = gxmap(:)
      gymap_glob(:) = gymap(:)
#endif

    ! Building pcmap

      allocate (pcmap(numpatch))

      write(6,*) 'min & max of ccmap', minval(ccmap), maxval(ccmap)

      p = 1
      do i = 1, numcolumn
         if(itypwat_glob(ccmap(i)).eq.0) then
            pcmap(p:p+16) = i
            p = p+17
         else
            pcmap(p) = i
            p = p+1
         end if
      end do

    ! Consistency check
      do i = 1, numpatch
         if(cgmap(pcmap(i)).ne.pgmap(i)) then
            write(6,*) 'fault error on pcmap, pgmap, cgmap',pcmap(i), cgmap(pcmap(i)), pgmap(i)
            call abort
         end if
      end do

      deallocate (cmask)
      deallocate (pmask)

   end subroutine task_decomp

   subroutine spmd_var_dealloc

      implicit none

!routine:

      deallocate (numgrid_proc)
      deallocate (numcolumn_proc)
      deallocate (numpatch_proc)

      deallocate (pgmap)

      deallocate (ppmap)
      deallocate (ccmap)
      deallocate (ggmap)
      deallocate (gxmap)
      deallocate (gymap)

      deallocate (ppmap_glob)
      deallocate (ccmap_glob)
      deallocate (ggmap_glob)
      deallocate (gxmap_glob)
      deallocate (gymap_glob)

   end subroutine spmd_var_dealloc

end module spmd_decomp
