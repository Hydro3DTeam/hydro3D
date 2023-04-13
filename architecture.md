# Architecture

### Creation and parallesation of computational domain:

The user input file that allows the creation and parallel segmentation of the overall computational domain into subdomains are:

* control.cin: define the refinement of the grid (dx,dy,dz)
* infodom.cin: define the number and the bounds of the subdomains and the LMR resolution.
* mdmap.cin: assign the subdomains to the given number of processors.

The files in the code ensuring the creation and connectivity between the computational subdomains:

* module\_mpi.for
* module\_multidata.for
* module\_vars.for
* alloc\_dom.for
* localparameters.for
* exchange.for
* exchange\_bc.for
* exchange\_phi.for
* exchangep.for
* exchangep



