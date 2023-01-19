# Home

Hydro3D is a finite difference Navier Stokes solver that performs accurate and efficient Large Eddy Simulation (LES) of turbulent flows. This software was used for a multitude of applications including turbulence in riverbeds, assessment of VATT and HATT turbine efficiency, and impact of fully submerged geometries onto the water's free-surface. Hydro3D has been developed over the years by Prof. Thorsten Stoesser and more than 20 collaborators.

#### Features

* Sub-grid Scale (SGS) turbulence modelling using Smagorinsky, 1-equation or WALE methods.
* Choice of spatial discretisation schemes including 4th order central differences.
* Choice of temporal solution schemes including 5th order ENO/WENO.
* Efficient multigrid solver with MPI and OpenMP parallelisation.
* Local Mesh Refinement for increased efficiency and accuracy (LMR).
* Passive and active scalar transport.&#x20;
* Free-surface flows simulation with Level-Set Method (LSM).
* Bubble column simulation using Lagrangian particle tracking (LPT).
* Importation of complex CAD geometries (.STEP or .msh format)&#x20;
* Optimised octree method to transform CAD format into an organised cloud of immersed boundary points. &#x20;
* Efficient computation of flows in and over complex geometries using an Immersed Boundary Method (IBM).

#### Ongoing development

* Non-Newtonian fluid solver.
* Improve the coupling of IBM and LSM for partially submerged structures.
* Sediment and scalar transport.
