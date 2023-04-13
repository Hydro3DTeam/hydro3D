# Architecture

<details>

<summary>CAVITY  &#x26; CHANNEL SIMULATION </summary>

The **user input file** that allows the creation and parallel segmentation of the overall computational domain into subdomains are:

* **control.cin:** define the refinement of the grid (dx,dy,dz)
* **infodom.cin:** define the number and the bounds of the subdomains and the LMR resolution.
* **mdmap.cin:** assign the subdomains to the given number of processors.

The **output files** to analyse the velocity field:

* tecbin\*.bin  -  BINARY which export the instantaneous, 1st and 2nd order time average flow field.
* tecturb\*.dat - ASCII Tecplot files which export the instantaneous, 1st and 2nd order time average flow field.
* tecinst\*\_\*.dat - ASCII Tecplot files which export the instantaneous flow field.
* tecplane\*\_\*.dat - ASCII Tecplot files which export the instantaneous flow field of a given plane.
* unst\_\*.dat - ASCII Tecplot Files which export the instantaneous flow field at a given probe location.

#### Module Files:

* module\_mpi.for
* module\_multidata.for
* module\_vars.for

**Domain Creation Files:**

* alloc\_dom.for
* initial.for
* localparameters.for

**Initial Flow Fuild and Boundary Files:**

* initial.for
* bounds.for

**Subgrid-Scale (SGS) Files:**

* eddyvis\_1eqn.for
* eddyvis\_keps.for
* eddyvis\_smag.for
* eddyvis\_wale.for

**Navier-Stoke Numerical Files:**

* convection.for
* diffusion.for
* rungek.for
* press.for
* weno.for
* newsolv\_mg.for
* mgsolver.for
* sipsol.for

**Output Files:**

* post.for
* timesig.for

</details>

### Creation and parallelisation of computational domain:

The **user input file** that allows the creation and parallel segmentation of the overall computational domain into subdomains are:

* **control.cin:** define the refinement of the grid (dx,dy,dz)
* **infodom.cin:** define the number and the bounds of the subdomains and the LMR resolution.
* **mdmap.cin:** assign the subdomains to the given number of processors.

The **files** in the code ensuring the creation and connectivity between the computational subdomains:

* module\_mpi.for
* module\_multidata.for
* module\_vars.for
* alloc\_dom.for
* localparameters.for
* exchange.for
* exchange\_bc.for
* exchange\_phi.for
* exchangep.for
* exchangepp.for
* exchangeu.for
* exchangev.for
* exchangew.for

### **Free-Surface Initialisation and development:**

The **user input file** used to run a simulation with a free-surface:

* control.cin
* infodom.cin
* mdmap.cin
* in\_lsm.cin

