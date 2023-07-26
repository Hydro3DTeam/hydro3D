#############################################################
F90 = mpiifort
OPTIONS    = -traceback -c -real-size 64 -O2 -qopenmp 
LOPTIONS   = -traceback -O2 -qopenmp 
##############################################################

objects = \
module_vars.o\
module_multidata.o\
module_mpi.o\
module_LPT.o\
module_LSM.o\
module_SEM.o\
module_IBM.o\
main.o\
IBM.o\
actuator.o\
shapes.o\
initial.o\
domain.o\
post.o\
checkdt.o\
bounds.o\
sipsolver.o\
conservation.o\
convection.o\
diffusion.o\
mgsolver.o\
wall_function.o\
deltafunction.o\
LPT.o\
timesig.o\
weno.o\
energy.o\
press.o\
roughness_function.o\
averaging.o\
eddyvis_smag.o\
eddyvis_wale.o\
eddyvis_1eqn.o\
eddyvis_keps.o\
exchangeBC.o\
exchangeP.o\
exchangePP.o\
exchangeS.o\
exchange.o\
exchangeU.o\
exchangeV.o\
exchangeW.o\
exchangePHI.o\
LSM.o\
SEM.o\
sediment.o\
parmove.o

.SUFFIXES: .f90

.f90.o:
	$(F90) $(OPTIONS) -o $@ $<

3dFDM.exe: $(objects) 
	$(F90) $(objects) $(LOPTIONS) -o 3dFDM_2023.exe \

clean:
	rm -rfv *.o *.mod 3dFDM_2023.exe


domain.o : domain.f90 module_mpi.o module_vars.o module_multidata.o 
averaging.o : averaging.f90 module_vars.o module_multidata.o 
bounds.o : bounds.f90 IBM.o module_mpi.o module_multidata.o module_vars.o 
checkdt.o : checkdt.f90 module_multidata.o module_mpi.o module_vars.o 
conservation.o: conservation.f90 module_multidata.o module_mpi.o module_vars.o
convection.o : convection.f90 module_multidata.o module_mpi.o module_vars.o 
deltafunction.o : deltafunction.f90 module_multidata.o module_vars.o module_mpi.o
diffusion.o : diffusion.f90 module_multidata.o module_mpi.o module_vars.o 
eddyvis_1eqn.o : eddyvis_1eqn.f90 module_multidata.o module_vars.o 
eddyvis_smag.o : eddyvis_smag.f90 module_multidata.o module_vars.o 
eddyvis_wale.o : eddyvis_wale.f90 module_multidata.o module_vars.o 
eddyvis_keps.o : eddyvis_keps.f90 module_multidata.o module_vars.o 
energy.o : energy.f90 module_multidata.o module_mpi.o module_vars.o 
exchangeBC.o : exchangeBC.f90 module_vars.o module_mpi.o module_multidata.o 
exchange.o : exchange.f90 exchangePHI.o module_mpi.o module_multidata.o module_vars.o 
exchangeP.o : exchangeP.f90 module_vars.o module_mpi.o module_multidata.o 
exchangePHI.o : exchangePHI.f90 module_vars.o module_mpi.o module_multidata.o 
exchangePP.o : exchangePP.f90 module_vars.o module_mpi.o module_multidata.o 
exchangeS.o : exchangeS.f90 module_vars.o module_mpi.o module_multidata.o 
exchangeU.o : exchangeU.f90 module_vars.o module_mpi.o module_multidata.o 
exchangeV.o : exchangeV.f90 module_vars.o module_mpi.o module_multidata.o 
exchangeW.o : exchangeW.f90 module_vars.o module_mpi.o module_multidata.o 
IBM.o : IBM.f90 module_IBM.o module_mpi.o module_multidata.o module_vars.o
actuator.o : actuator.f90 module_IBM.o module_mpi.o module_multidata.o module_vars.o  
initial.o : initial.f90 module_mpi.o module_multidata.o module_vars.o 
LPT.o : LPT.f90 module_vars.o module_mpi.o module_multidata.o module_LPT.o
LSM.o : LSM.f90 module_mpi.o module_multidata.o module_vars.o module_LSM.o
sediment.o: sediment.f90 module_multidata.o module_vars.o module_mpi.o
main.o: main.f90 module_multidata.o module_mpi.o module_vars.o module_LPT.o
mgsolver.o : mgsolver.f90 module_multidata.o module_vars.o module_mpi.o
module_mpi.o : module_mpi.f90 
module_multidata.o : module_multidata.f90 
module_LPT.o : module_LPT.f90
module_LSM.o : module_LSM.f90
module_vars.o : module_vars.f90 
module_SEM.o     : module_SEM.f90
module_IBM.o : module_IBM.f90 
newsolv_mg.o : newsolv_mg.f90 module_multidata.o module_mpi.o module_vars.o 
post.o : post.f90 module_vars.o module_multidata.o 
press.o : press.f90 module_multidata.o module_mpi.o module_vars.o 
roughness_function.o : roughness_function.f90 module_multidata.o module_mpi.o module_vars.o 
shapes.o : shapes.f90 module_mpi.o IBM.o module_multidata.o module_vars.o 
sipsolver.o : sipsolver.f90 module_mpi.o module_multidata.o module_vars.o 
timesig.o : timesig.f90 module_mpi.o module_vars.o module_multidata.o 
wall_function.o : wall_function.f90 module_multidata.o module_vars.o 
weno.o : weno.f90 module_multidata.o module_vars.o 
SEM.o : SEM.f90 module_multidata.o module_vars.o module_SEM.o module_mpi.o
parmove.o : parmove.f90 module_multidata.o module_vars.o module_mpi.o 
