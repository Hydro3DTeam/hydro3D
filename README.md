# VERSION 2023 Hydro3D

## = IBM =
* Geometry:
    * Now it is possible to import solidwork or gmsh geometry into Hydro3D ! 
    * Input file: in_geom.cin	(altered)
    * Script: cadtoibps.for 	(new)
	      ibm.for 		(altered)
	      module_ibm.for	(altered)
	      initial.for	(altered)

* Speed & Scalability:
  * The MPI load and transfer between the processor has been reduced to a minimum.
  * Now for fix geometry each processor assess the number of IBM and no information is transfer through the MPI_BCAST.
  * Limitation: It is unclear if the actuator line is used.
  * Input file: in_geom.cin 	
  * Script: ibm.for 		(altered)
            DeltaMLS.for 	(altered)
	    module_ibm.for 	(altered)
	    actuator.for 	(altered)

## = EULERIAN FIELD =
### Celltype:
  * This new function allows Hydro3D to distinguish three type of cells
    * -1: Emptycell (the velocity field is forced to be null=0)
    * 0: IBcell it allows a better calculation of the mass flow deficit (the cell has a immersed boundary point) 
    * 1: Fluid cell (normal cell)
  * Input file: in_celltype.cin (new)
  * Script: initial.for 	(altered)
	    press.for		(altered)
	    module_vars.for 	(altered)

### Field Inlet/Outlet:
  * This part of the code has been developed to input flow where ever in the field.
  * Depending on being an inlet or outlet the flow is accounted in the mass calculation.
  * Specific IB geometry can be imported and used as IB I/O, these geometry will not be account in the IBM.
  * Limitation: Challenging to add inclined I/O, you need at least two row of I/O for it to function correctly.
  * Input: HRS.cin  (new)
  * Script: HRS.for 		(new)
            module_HRS.for 	(new)
	    initial.for 	(altered)

### Scalar:
  * The boundary conditions have been simplified. This parameter were remove from the control.cin and moved to the new scalar input file.
  * Multiple field of active and passive scalar can now being created. (EG. Temperature & Concentration)
  * Input: scalar.cin		(new)
  * Script: scalar.for   	(Merge of: sediment.for & energy.for)
	    sediment.for 	(removed)
	    energy.for 	 	(removed)
	    initial.for 	(altered)
	    module_scalar.for 	(new)
	    module_vars.for 	(altered)

## = LPT =
* Evaporation model:
  * Now the particle can be subjected to eveporation. Coupled with the scalar field.
* The density of the particle & surrounding fluid can be specified in the input file.
* Input file: in_LPT.for (altered) 
* Script: LPT.for	 	(altered)
	  initial.for 		(altered)
          module_LPT.for	(altered)
	  module_vars.for	(altered)

## = EXPORT =
* The data file are not exported as .plt but as .dat (.plt is supposed to be for the binary file)
* For any instantaneous files can be exported at a specific time in control.cin
* Probes does not need to be specified as dom,i,j,k and should be specified as x,y,z and hydro3D will determine their domain location automatically.
* Export VTK type for Paraview uses. (https://visit-sphinx-github-user-manual.readthedocs.io/en/task-allen-vtk9_master_ospray/data_into_visit/VTKFormat.html)

## = PROBE =
* Probes can also have a lagrangian location, weight will be interpolated using trilinearity. 
  Computer: H3D10
  Sim Name: Coarse_LPT_16
  FLAG: ! PROBE UPDATE
  Input file: control.cin (altered)
  Script: initial.for (altered)
	  timesig.for (altered)
	  module_vars (altered)
