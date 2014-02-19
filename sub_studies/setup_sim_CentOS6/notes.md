
E1F Nominal Beam Energy: 5.479 GeV
(http://clasweb.jlab.org/rungroups/e1/e1f/lists/LOAF_all.txt)

E1F MomCorr Beam Energy=5.49918 GeV 
(Marco Mizarita)

To compare MMs from ER and SR, for various BEs and GPP parameters
==================================================================

1. To compare ER-MMs for different BEs
	i.  In `aexp-2pi2`, ana_e1f.exp.2pi
			+ Makes `d2pi.root` which needs to be copied over as `d2pi_beXXXX.root`(by default ER-BE=5.499 GeV)
	ii. make_d2pi_sim.2pi.study [BE]:
			+ Makes $SETUPSIMCENTOS6_DATADIR /d2pi_beXXXX.root 
	iii. `study_BE_GPP_affect_MMs.py::anaMM2pi_ER(mmtyp="mm2")` to compare ER-MMs for different ER-BEs

2. To compare ER-,SR-MMs for different BEs and GPP pars
	i.   In `aexp-2pi2`, ana_e1f.exp.2pi
		 	+ Makes `d2pi.root` which needs to be copied over as `d2pi_beXXXX.root`(by default ER-BE=5.499 GeV)
	ii.  `make_d2pi_sim [BE]`
	iii. `study_BE_GPP_affect_MMs.py::anaMM2pi_ERSR(be_exp,be_sim,dtyps=28)` to compare ER-,SR-MMs for a particular BEs and range of GPP pars 
			+ Saves plots in $SETUPSIMCENTOS6_DATADIR/../../be_expXXXX_be_simXXXX`
	


To study Simulated-Detector response
====================================
1. `make_d2pi_sim.2pi.study [BE]` 
   Makes .root file with SR and ST trees and puts it in $SETUPSIMCENTOS6_DATADIR
2. `nb-study_E1F_Simulation-beR-[BE]` to study Simulated Detector 