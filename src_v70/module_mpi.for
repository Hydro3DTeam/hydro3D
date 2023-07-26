!##########################################################################
      module mpi
!##########################################################################
      implicit none
      include 'mpif.h'

      SAVE
      INTEGER           :: nprocs     ! nr. of processors
      INTEGER           :: myrank      ! nr of this processor
      INTEGER           :: ierr
      INTEGER           :: MPI_FLT  
      INTEGER           :: status(MPI_STATUS_SIZE)
      REAL ( kind =8 )  :: wtime
      REAL ( kind =8 )  :: wtime2

      INTEGER,parameter :: tag_offset1=10000
      INTEGER,parameter :: tag_offset2=500
      INTEGER,parameter :: tag_offset3=150000      ! Reduce this number if compiler fails

      contains
!#########################################################################
      subroutine init_parallelisation
!#########################################################################
      implicit none

!  initialize MPI.
        call MPI_INIT ( ierr )
!  Get the number of processes.
        call MPI_COMM_SIZE ( MPI_COMM_WORLD, nprocs, ierr )
!  Get the individual process iD.
        call MPI_COMM_RANK ( MPI_COMM_WORLD, myrank, ierr )

        MPI_FLT   = MPI_DOUBLE_PRECISION

        if ( myrank == 0 ) then
           wtime = MPI_WTIME ( )
           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) 'HELLO_MPI - Master process:'
           write ( *, '(a)' ) '  FORTRAN90/MPI version'
           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) '  A 3D finite Difference MPI program.'
           write ( *, '(a)' ) ' '
           write ( *, '(a,i8)' ) '  The number of processes is ', nprocs
           write ( *, '(a)' ) ' '
        end if

        end subroutine init_parallelisation
!##########################################################################
      subroutine end_parallelisation
!##########################################################################
      implicit none

        if ( myrank == 0 ) then
           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) '3D_FD_MPI - Master process:'
           write (*,'(a)')'Normal end of execution:"Goodbye, world!".'
           wtime = MPI_WTIME ( ) - wtime
           write ( *, '(a)' ) ' '
           write (*,'(a,g14.6,a)') 'Elapsed wall clock time = ',wtime,
     &'seconds.'
        end if

        call MPI_BARRIER (MPI_COMM_WORLD,ierr)

        call MPI_FINALIZE ( ierr )

      end subroutine end_parallelisation
!##########################################################################

      end module mpi
!##########################################################################
