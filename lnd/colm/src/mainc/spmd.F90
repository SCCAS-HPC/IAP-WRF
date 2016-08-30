#include <define.h>

module spmd

#ifndef SPMD

  implicit none

  integer :: p_iam = 0
  integer :: p_nprocs = 1
  logical :: p_master = .true. 
  logical :: p_slave  = .false.

  integer :: p_comm

#else

! use mpi

  implicit none

  include 'mpif.h'

  integer :: p_iam         !proc number
  integer :: p_nprocs      !number of processors
  logical :: p_master      !proc 0 logical for printing msgs
  logical :: p_slave 

  integer :: p_err
  integer :: p_tag
  integer :: p_comm
  integer :: p_stat(mpi_status_size)

  integer :: p_req(2)

  integer, pointer :: p_displs(:)
  integer, pointer :: p_counts(:)
  character(len=MPI_MAX_PROCESSOR_NAME), pointer :: p_nodes(:)

contains

  subroutine spmd_init(mpicom_lnd)

     integer, optional, intent(in) :: mpicom_lnd

     integer i,j 

#ifdef COUP_CSM
#ifdef CPL6
    ! Initialize mpi and set communication group
    ! Done in CLM.F90 => csm_setup(p_comm)
#endif
#ifdef CPL7
     if(present(mpicom_lnd)) then
        p_comm = mpicom_lnd
     end if
#endif
#else
     call mpi_init(p_err)
     p_comm = MPI_COMM_WORLD
#endif

     p_tag  = 0
     p_req(:) = MPI_REQUEST_NULL

! Get my processor id  

     call mpi_comm_rank(p_comm, p_iam, p_err)  

     if(p_iam==0)then 
        p_master = .true.
        p_slave  = .false.
     else
        p_master = .false.
        p_slave  = .true.
     endif

! Get number of processors

     call mpi_comm_size(p_comm, p_nprocs, p_err) 

! Get my processor names

     allocate (p_counts(0:p_nprocs-1))
     allocate (p_displs(0:p_nprocs-1))
     allocate (p_nodes(0:p_nprocs-1))

     call mpi_get_processor_name (p_nodes(p_iam),p_counts(p_iam),p_err)
     call mpi_allgather (p_counts(p_iam),1,mpi_integer,p_counts,1,mpi_integer,p_comm,p_err)

     do i=0,p_nprocs-1
        p_displs(i)=i*MPI_MAX_PROCESSOR_NAME
     enddo

     call mpi_gatherv (p_nodes(p_iam),p_counts(p_iam),mpi_character, &
                       p_nodes,p_counts,p_displs,mpi_character,0,p_comm,p_err)

     if (p_master) then
        write(6,100)p_nprocs
        write(6,200)
        write(6,220)
        do i=0,p_nprocs-1
           write(6,250)i,(p_nodes((i))(j:j),j=1,p_counts(i))
        enddo
     endif

100  format(i3," pes participating in computation")
200  format(/,35('-'))
220  format(/,"NODE#",2x,"NAME")
250  format("(",i3,")",2x,100a1)

  end subroutine spmd_init

  subroutine spmd_exit

     deallocate (p_counts)
     deallocate (p_displs)
     deallocate (p_nodes)

#ifndef COUP_CSM
     call mpi_barrier(p_comm,p_err)
     call mpi_finalize(p_err)
#endif

  end subroutine spmd_exit

#endif

end module spmd

