# ActuatorLine\_Initial

File: actuator.for

```
######################################################################
      SUBROUTINE ActuatorLine_Initial(M)
!######################################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
	integer::M
	if(itime.eq.itime_start .and. myrank.ne.0) then
        call MPI_BCAST(c_act,maxnodeIBS,MPI_DOUBLE_PRECISION,
     &			master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Pit_act,maxnodeIBS,MPI_DOUBLE_PRECISION,
     &			master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(r_act,maxnodeIBS,MPI_DOUBLE_PRECISION,
     &			master,MPI_COMM_WORLD,ierr)
	endif
	END SUBROUTINE
```
