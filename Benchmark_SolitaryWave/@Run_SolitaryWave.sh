#!/bin/bash

cpu=8
subdomains=42
itime=100
etime=15000
dtime=100
nprocs=80
Smax=0
filename='input_file.txt'

# Run_simulation:
mpirun -np $cpu ./3dFDM_2023.exe

# Computational result:
cp rms.dat worktime.dat output.dat forcn*.dat *.cin LSM_*.dat Results 


