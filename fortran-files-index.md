# Fortran Files Index

<details>

<summary>actuator.for</summary>

n_**Purpose:**_ Computational method to represent an array of turbines for a lower cost.

_**Difficulty:**_ Hard | _**The user change likelihood:**_ Unlikely

_**Subroutines:**_

* act\_line\_vatt\_geom
* act\_line\_vatt
* act\_line\_geom
* actuatorline\_initial
* actuatorline
* actuatorline\_FEM

</details>

<details>

<summary>alloc_dom.for</summary>

_**Purpose:**_ Create the computational domain based on infodom.cin and mdmap.cin.

_**Difficulty:**_ Medium | _**The user change likelihood:**_ Very Unlikely

_**Subroutines:**_

* alloc\_dom
* read\_mdmap
* read\_infodom
* datainfo

</details>

<details>

<summary>averaging.for</summary>

_**Purpose:**_ Average all the flow field variables.

_**Difficulty:**_ Easy | _**The user change likelihood:**_ Unlikely

_**Subroutines:**_

* update\_mean
* add\_noise

</details>

<details>

<summary>bounds.for</summary>

_**Purpose:**_ Specific boundary conditions are applied to the flow field.

_**Difficulty:**_ Easy | _**The user change likelihood:**_ Very Likely

_**Subroutines:**_

* boundu
* boundv
* boundw
* boundcoeff

</details>

<details>

<summary>bounds_keps.for</summary>

_**Purpose:**_ Specify boundary conditions for RANS simulations.

_**Difficulty:**_ Easy | _**The user change likelihood:**_ Unlikely

_**Subroutines:**_

* boundksgs
* boune

</details>

<details>

<summary>bounds_lsm.for</summary>

_**Purpose:**_ Specify the boundary conditions for LSM simulations.

_**Difficulty:**_ Medium | _**The user change likelihood:**_ Likely

</details>

<details>

<summary>cadtoibp.for</summary>

_**Purpose:**_ Transform GMsh (.msh) or CAD (.STEP) into Hydro3D geometric file.

_**Difficulty:**_ High | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>checkdt.for</summary>

_**Purpose:**_ Check the CFL and adapt the time step size if set up as variable in control.cin

_**Difficulty:**_ Medium | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>convection.for</summary>

_**Purpose:**_ Numerical discretization of the convection term N-S.

_**Difficulty:**_ Medium | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>DeltaF_MLS.for</summary>

_**Purpose:**_ Calculate the delta interpolation function to couple Lagrangian points and the Eulerian mesh.

_**Difficulty:**_ Medium | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>diffusion.for</summary>

_**Purpose:**_ Numerical discretization of the diffusion term N-S.

_**Difficulty:**_ Medium | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>eddyvis_1eqn.for</summary>

_**Purpose:**_ Calculate the SGS viscosity created by the turbulence lower than the filter. Using the one-equation.

_**Difficulty:**_ Medium | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>eddyvis_keps.for</summary>

_**Purpose:**_ Calculate the SGS viscosity created by the turbulence lower than the filter. For RANS.

_**Difficulty:**_ Medium | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>eddyvis_smag.for</summary>

_**Purpose:**_ Calculate the SGS viscosity created by the turbulence lower than the filter. Using the Smagorosky equation.

_**Difficulty:**_ Medium | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>eddyvis_wale.for</summary>

_**Purpose:**_ Calculate the SGS viscosity created by the turbulence lower than the filter. Using the WALE algorithm.

_**Difficulty:**_ Medium | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>exchange.for</summary>

_**Purpose:**_ Menu to select the variable to exchange between the ghost-cell using MPI.

_**Difficulty:**_ Medium | _**The user change likelihood:**_ Unlikely

</details>

<details>

<summary>exchange_bc.for</summary>

_**Purpose:**_ Used for periodic boundaries to exchange data between the inlet and outlet of the main domain.

_**Difficulty:**_ High | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>exchange_bcphi.for</summary>

_**Purpose:**_ Used for periodic boundaries to exchange free-surface data between the inlet and outlet of the main domain.

_**Difficulty:**_ High | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>exchange_phi.for</summary>

_**Purpose:**_ Exchange the phi variable (free-surface) between the neighbouring subdomain ghost-cells using the MPI.

_**Difficulty:**_ High | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>exchangep.for</summary>

_**Purpose:**_ Exchange the p variable (pressure) between the neighbouring subdomain ghost-cells using the MPI.

_**Difficulty:**_ High | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>exchangepp.for</summary>

_**Purpose:**_ Exchange the pp variable (pseudo-pressure) between the neighbouring subdomain ghost-cells using the MPI.

_**Difficulty:**_ High | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>exchangesca.for</summary>

_**Purpose:**_ Exchange the sca variable (scalar) between the neighbouring subdomain ghost-cells using the MPI.

_**Difficulty:**_ High | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>exchangeu.for</summary>

u_**Purpose:**_ Exchange the u variable (streamwise-velocity) between the neighbouring subdomain ghost-cells using the MPI.

_**Difficulty:**_ High | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>exchangev.for</summary>

_**Purpose:**_ Exchange the v variable (spanwise-velocity) between the neighbouring subdomain ghost-cells using the MPI.

_**Difficulty:**_ High | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>exchangew.for</summary>

_**Purpose:**_ Exchange the w variable (vertical-velocity) between the neighbouring subdomain ghost-cells using the MPI.

_**Difficulty:**_ High | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>fdstag.for</summary>

_**Purpose:**_ The main skeleton of the code to run the simulations.

_**Difficulty:**_ Easy | _**The user change likelihood:**_ Very Likely

</details>

<details>

<summary>flosol.for</summary>

_**Purpose:**_ Squeletton to run each time-step.

_**Difficulty:**_ Easy | _**The user change likelihood:**_ Very Likely

</details>

<details>

<summary>HRS.for</summary>

_**Purpose:**_ Prescribe mass inflow or outflow at a specific domain location.

_**Difficulty:**_ Medium | _**The user change likelihood:**_ Likely

</details>

<details>

<summary>IBM.for</summary>

_**Purpose:**_ Perform the IBM to enforce a no-slip condition at the Lagrangian boundary of a geometry.

_**Difficulty:**_ High | _**The user change likelihood:**_ Unlikely

</details>

<details>

<summary>initial.for</summary>

_**Purpose:**_ Initialise most of the variable and initial field conditions of the simulation.

_**Difficulty:**_ Medium | _**The user change likelihood:**_ Very Likely

</details>

<details>

<summary>localparameters.for</summary>

_**Purpose:**_ Evaluate the multigrid level at which the simulation can be run. Check the LMR mapping.

_**Difficulty:**_ Medium | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>log_law.for</summary>

_**Purpose:**_ Provide log\_law boundary condition at each time step.

_**Difficulty:**_ Medium | _**The user change likelihood:**_ Unlikely

</details>

<details>

<summary>LPT.for</summary>

_**Purpose:**_ Initialise and perform the calculation for Lagrangian particles.

_**Difficulty:**_ High | _**The user change likelihood:**_ Unlikely

</details>

<details>

<summary>LSM.for</summary>

_**Purpose:**_ Initialise and perform the free-surface calculation at each time step.

_**Difficulty:**_ High | _**The user change likelihood:**_ Unlikely

</details>

<details>

<summary>mgsolver.for</summary>

_**Purpose:**_ Calculate the pressure from the velocity field using the poisson-pressure solver.

_**Difficulty:**_ High | _**The user change likelihood:**_ Unlikely

</details>

<details>

<summary>module_HRS.for</summary>

_**Purpose:**_ Declare global variables for the HRS.

_**Difficulty:**_ Easy | _**The user change likelihood:**_ Very Likely

</details>

<details>

<summary>module_IBM.for</summary>

_**Purpose:**_ Declare global variables for the IBM.

_**Difficulty:**_ Easy | _**The user change likelihood:**_ Very Likely

</details>

<details>

<summary>module_MPI.for</summary>

_**Purpose:**_ Declare and initialise global variables for the MESSAGE PASSING INTERFACE.

_**Difficulty:**_ Easy | _**The user change likelihood:**_ Very Likely

</details>

<details>

<summary>module_multidata.for</summary>

_**Purpose:**_ Declare the eulerian structure dom(ib) variables.

_**Difficulty:**_ Easy | _**The user change likelihood:**_ Very Likely

</details>

<details>

<summary>module_LPT.for</summary>

_**Purpose:**_ Declare global variables for the LPT.

_**Difficulty:**_ Easy | _**The user change likelihood:**_ Very Likely

</details>

<details>

<summary>module_LSM.for</summary>

_**Purpose:**_ Declare global variables for the LSM.

_**Difficulty:**_ Easy | _**The user change likelihood:**_ Very Likely

</details>

<details>

<summary>module_vars.for</summary>

_**Purpose:**_ Declare global variables for the basic simulations.

_**Difficulty:**_ Easy | _**The user change likelihood:**_ Very Likely

</details>

<details>

<summary>newsolv_mg.for</summary>

_**Purpose:**_ Run iteratively the poisson-pressure solver, and export the step print.

_**Difficulty:**_ Very High | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>parmove.for</summary>

_**Purpose:**_ Calculation for the bed sedimentation.

_**Difficulty:**_ High | _**The user change likelihood:**_ Unlikely

</details>

<details>

<summary>post.for</summary>

_**Purpose:**_ Export all the data files of the simulation.

_**Difficulty:**_ Easy | _**The user change likelihood:**_ Very Likely

</details>

<details>

<summary>press.for</summary>

_**Purpose:**_ Calculate the fractional-step velocity after the pressure-solver. Performed SIP solver.

_**Difficulty:**_ Medium | _**The user change likelihood:**_ Unlikely

</details>

<details>

<summary>roughness_function.for</summary>

_**Purpose:**_ Calculate the roughness function for porous beds.

_**Difficulty:**_ Medium | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>rungek.for</summary>

_**Purpose:**_ Numerical method for convection and diffusion spatial terms.

_**Difficulty:**_ High | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>scalar.for</summary>

_**Purpose:**_ Initial and perform the calculation for the scalar fields.

_**Difficulty:**_ High | _**The user change likelihood:**_ Likely

</details>

<details>

<summary>SEM.for</summary>

_**Purpose:**_ Initialise and perform calculations for SEM.

_**Difficulty:**_ High | _**The user change likelihood:**_ Unlikely

</details>

<details>

<summary>shape.for</summary>

_**Purpose:**_ Create specific shapes geometry.

_**Difficulty:**_ Medium | _**The user change likelihood:**_ Unlikely

</details>

<details>

<summary>sipsol.for</summary>

_**Purpose:**_ Perform Stone Implicit Pressure solver.

_**Difficulty:**_ Very High | _**The user change likelihood:**_ Very Unlikely

</details>

<details>

<summary>timesig.for</summary>

_**Purpose:**_ Export the variables for each probe.

_**Difficulty:**_ Easy | _**The user change likelihood:**_ Likely

</details>

<details>

<summary>wall_function.for</summary>

_**Purpose:**_ Calculate wall-function boundary conditions.

_**Difficulty:**_ Easy | _**The user change likelihood:**_ Unlikely

</details>

<details>

<summary>weno.for</summary>

_**Purpose:**_ Perform the WENO differencing scheme.

_**Difficulty:**_ Very High | _**The user change likelihood:**_ Very Unlikely

</details>
