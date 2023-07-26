#!/bin/bash

cpu=8
subdomains=8
nprocs=8
Smax=0
filename='input_file.txt'

# Run_simulation:
mpirun -np $cpu ./3dFDM_2023.exe

# Computational result:
cp rms.dat worktime.dat output.dat forcn*.dat *.cin Results 

# Write input file for postprocess:
echo "$subdomains,$nprocs" > "$filename"

# Run_post_processing:
cp ../PostProcessing_Script/POSTPROCESS_OpenMP.f90 .
ifort -traceback -qopenmp -parallel -fpp POSTPROCESS_OpenMP.f90 -o POSTPROCESS_OpenMP.exe
./POSTPROCESS_OpenMP.exe 


# Delete input file:
rm $filename
rm POSTPROCESS*
rm tec*.dat

