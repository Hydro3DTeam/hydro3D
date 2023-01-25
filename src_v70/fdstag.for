!#######################################################################
      program fdstag
!-----------------------------------------------------------------------
!     It is here where every simulation starts. All the subroutine
!     prior to flosol set up the initial condition for the simulation.
!     The flosol (flow solver) is the subroutine where every step time
!     of the simulation are performed (loop).
!#######################################################################
      use mpi
      use vars
      implicit none

!Starts the parallelisation:
      call init_parallelisation
     
!Reads multi-domain map:
      call read_mdmap

!Reads control file:
      call read_control

!Reads the information of each subdomain:
      call read_infodom

!Creation of the subdomain allocation to their respective processor:
      call alloc_dom

!Mesh generation:
      call localparameters

!Allocation of the main flow variables:
      call initial

!LSM: initialisation of the level set field and allocation of variables:
      IF (L_LSM .or. L_LSMbase)  CALL initial_LSM_3D_channel

!Initialisation of the flow field:
      call initflowfield

!Rough bed:
      if (LROUGH) then								
        if (.not.LRESTART) then
          call init_rough
        else
          call rough_restart
        end if
      end if

!IBM: reads geom.cin, genetare the geometries, opens files, etc.
      if (LIMB) call imb_initial

!IBM: Assignation of the Lagrangian markers to sub-domains and CPUs:
      if (LIMB) call PartLoc_Initial							

!Calculation of the initial flux in the simulation:
        call iniflux

!SCALAR: initial sediment concentration variables:
	 if(LSCALAR)   call sediment_init

        if(.not.LRESTART) then
           if (time_averaging) then
              call update_mean
	      if (noise.gt.0.0) call add_noise(noise)
           end if
        end if

!LPT: initialisation of the particles distribution:
	  if (LPT) then						
		if (myrank.eq.0) open(unit=202,file='LPT_particles.dat')
		call init_particle	
	  endif

        if ((solver.eq.2).and.(.not.L_LSM)) call coeff

        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
           write (numfile,*) '============START ITERATIONS========='
           write (6,*) '============START ITERATIONS========='
        end if

!All fields are initialised, now let's iterate the N-S eqs in time:
        call flosol

!End MPI parallelisation:
        call end_parallelisation

      end program
!#######################################################################
