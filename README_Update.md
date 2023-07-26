## = File Structure = 
* Renamed:
  * alloc_dom.f90 -> domain.f90
  * lsm.f90 -> LSM.f90
  * ibm.f90 -> IBM.f90
  * module_ibm -> module_IBM.f90
  * Deltah_MLS.f90 -> deltafunction.f90
  * exchange_bc.f90 -> exchangeBC.f90
  * exchangep.f90 -> exchangeP.f90
  * exchangepp.f90 -> exchangePP.f90
  * exchangesca -> exchangeS.f90
  * exchangeu - > exchangeU.f90
  * exchangev - > exchangeV.f90
  * exchangew - > exchangeW.f90
  
* Merged:
  * newsolv_mg.f90 -> mgsolver.f90
  * flosol.f90 + fdstag.f90 -> main.f90
  * log_law.f90 -> wallfunction.f90
  * bounds_lsm.f90 -> LSM.f90
  * bounds_keps.f90 -> eddyvis_keps.f90
  * localparameters.f90 -> domain.f90
  
* Split:
  * rungek.f90 -> convection.f90 
    * rungek_conv2nd
    * rungek_conv4th
    * rungek_convWENO
  * rungek.f90 -> diffusion.f90
    * rungek_diff2nd
    * rungek_diff4th
    
* New:
  * conservation.f90
    * iniflux (initial.f90)
    * correctoutflux (initial.f90)
    * pressure_forcing (press.f90)
  * module_LPT.f90
    * module_vars (module_vars.f90)
  * module_LSM.f90
    * module_lsm (LSM.f90)

## = Optimisation = 
* Optimisation hierachic between loop and conditions:
  * bounds.f90: boundu, boundv, boundw  (if: L_LSM)
  * checkdt.f90: checkdt (if: SENERGY)
  * conservation.f90: iniflux, correctoutflux, pressure_forcing (if: L_LSM)
  * convection.f90: convection (if: differencing, conv_sch)
  * diffusion.f90: diffusion, rungek_diff2nd, rungek_diff4th (L_LSM, pressureforce)
  * mgsolver.f90: coeff, mgkcyc (if: L_LSM)
  * press.f90: calmas, calvel (if: L_LSM)
  * averaging.f90: averaging, add_noise (if: L_LSM)

* Optimisation of array leakage:
  * The declaration of the array using pointer can lead to memory leakage, it should be used only for fi =>
  * module_multidata.f90: pointer -> allocatable
  * .f90: local array: pointer -> allocatable

* Change doubt:
  * mgsolver.f90: coeff
  * bounds.f90: boundcoeff [ dealloacte(fi) ]

* Issue:
  * bounds.f90: boundu, boundv, boundw if(2,21: LSM take the whole field)
  * diffusion.f90: diffusion, rungek_diff2nd, rungek_diff4th (the way the pressure force is prescribed ) 
    

