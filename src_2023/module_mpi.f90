!#######################################################################
      module mpi
!#######################################################################
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
!#######################################################################
      SUBROUTINE init_parallelisation
!-----------------------------------------------------------------------
!     Initialise the Message Passage Interface (MPI)  
!#######################################################################
      implicit none

!  Initialize MPI.
      call MPI_INIT(ierr)

!  Get the number of processors/cpu thread.
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

!  Get the individual processor/cpu thread iD.
      call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

      MPI_FLT = MPI_DOUBLE_PRECISION

      if (myrank .eq. 0) then
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
      
      RETURN
      END SUBROUTINE init_parallelisation
!#######################################################################
      SUBROUTINE end_parallelisation
!-----------------------------------------------------------------------
!     Terminate the MPI process
!#######################################################################
      implicit none

      if (myrank .eq. 0) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '3D_FD_MPI - Master process:'
        write (*,'(a)')'Normal end of execution:"Goodbye, world!".'
        wtime = MPI_WTIME ( ) - wtime
        write ( *, '(a)' ) ' '
        write (*,'(a,g14.6,a)') 'Elapsed wall clock time = ',wtime,'seconds.'
      end if

      call MPI_BARRIER (MPI_COMM_WORLD,ierr)

! Terminate the MPI process: 
      call MPI_FINALIZE(ierr)
      
      RETURN
      END SUBROUTINE end_parallelisation
!#######################################################################
      end module mpi
!#######################################################################
