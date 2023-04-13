---
description: Developed by Arthur Hajaali
---

# Ahmed Body

This aerodynamic benchmark ensures the integrity of the Immersed Boundary Method (IBM) and Local Mesh Refinement (LMR).  For the new version user, it will teach them the workflow to integrate any complex CAD model into Hydro3D using GMsh.

<mark style="color:red;">**Prerequisite:**</mark> The user is within docker or the user installed GMsh on his computer.

### **Transform CAD \*.step file into a triangular mesh (\*.msh) -  Only new version**

The following command is the first step towards transforming any CAD into an organised cloud of points required by Hydro3D to represent geometry.

```
gmsh --version                 #check that GMsh is installed
```

Now let's convert our parametric CAD file Ahmed\_Body.STEP into a mesh:

```
gmsh -match -clmax 20 -1 -2 -3 Ahmed_Body.STEP -o Ahmed_Body_0.02m.msh | tee Ahmed_Body_0.02m.log
```

The command can be intimidating at first, so let's break it down:

<table data-header-hidden><thead><tr><th></th><th></th><th data-hidden></th></tr></thead><tbody><tr><td>-match</td><td>selects all the geometry features found in Ahmed_Body.STEP</td><td></td></tr><tr><td>-clmax 20</td><td>contrains the maximum length between the mesh points to be 20mm</td><td></td></tr><tr><td>-1</td><td>performs the meshing of the Ahmed_Body edges (1D)</td><td></td></tr><tr><td>-2 </td><td>performs the meshing of the Ahmed_Body surfaces (2D)</td><td></td></tr><tr><td>-3</td><td>performs the meshing of the Ahmed_Body volume (3D)</td><td></td></tr><tr><td>-o</td><td>inform the name of the msh file</td><td></td></tr><tr><td>| tee</td><td>save the terminal print into the Ahmed_Body_0.02m.log</td><td></td></tr></tbody></table>

### Run the simulation:

```
mpirun -np 8 ./3dFDM.exe     # Start the simulation

In a new window:
tail -f worktime.dat         # Shows you the time that each section of the code takes
tail -f rms.dat              # Shows you the mass deficit of the simulation
htop                         # Shows you the ram and cpu usage of the simulation
```

### Analyse the Simulation Tecplot:

```
tec360 LAYOUT_Ahmed_Body.lay
```

### Analyse the data with Paraview:

```
paraview STATE_Ahmed_Body.pvsm
```
