Study objectives
================
To understand CLAS & Simulated-CLAS detector response and therefore understand differences in some key histograms that
are obtained by putting together measured kinematics. 

This study started when I noticed that distribution of ER- & SR-MM distributions did not match.  The ER-MMs are closer to the Expected Distribution and SR-MMs are shifted and narrower. To fix the SR-MMs, I tried:
	1. Tuning GPP pars
	2. Playing with SR-BE
	3. Applying Energy Loss corrections

[02-19-14] 
As it stands, even after correcting for shifts in SR hadron momentum by way of applying eloss corrections, shift in SR e- momentum will remain, causing the z-MM to remain offset and therfoere, the SR-MMs to remain Offset from ER-MMs. The shift in SR e- z-momentum cannot simply be accomodated by making a corresponding shift in SR-BE; it may well appear to fix the z-MM component which changes linearly with pz_e', however the t-MM component suffers even more strongly (see hand written notes), which cannot be accomodated by the shift in SR-BE. Therefore, the options may be:
	1. Shift SR-BE by an amount equal to shift in SR pz_e'; atleast this will fix the 3-momentum component of MM and therefore the 3-	momentum of the Missing Particle. However, the MM cuts in ER and SR, will be different (but have similar effects)
	2. Shift SR-BE by the amount needed to make the SR-MMs agree with ER-MMs; in this case the MM cuts in ER and SR will be the same, 
	   but the 3-momentum of the missing particle will be incorrect (But according to Gothe, I may not be sensitive to this shift, which I should verify)


1.Tuning GPP pars
-----------------
For now, I am using GPP pars picked by Evan (1.37,1.37,1.37). I tried 27 combinations in GPP-par space spanned by:
(0,0,0)*(1.37,1.37,1.37)*(4,4,4)
Among all of the the parameters used by Evan, matched ER-MMs most closely; however as mentioned before, they were shifted and narrower
compared with ER-MMs. Items 2 and 3 address the shift, but the narrowness of the MM distributions may need to be addressed by picking different GPP pars, but I will revisit this once 3 is done. 


2. Playing with SR-BE [02-19-14]
--------------------------------
I tried to change the SR-BE as a means for applying an "Effective Correction". But the SR-BE I had to pick was 5.485 GeV, which was 14 MeV off from ER-BE(5.499 GeV). This made me wonder what the matter was and led me to study the Simulated Detector response (See below). There I found that z component of the measured momentum (pz) was shifted viz-a-viz ST momentum. In the Q2W range = q2w2, this shift was ~ 3 Mev for e-,~2 MeV and almost neglible for the pip/pim. Although this shift is seems little, it affects the time componenent of the 4-vecs more strongly, since there the momentum components enter in quadratically and therefore, the error due to this little offset manifests itself more strongly (See handwritten notes for details)

3. Applying Energy Loss corrections
-----------------------------------
This may not solve the offset SR-MM offset problem, but should reduce it, since it should correct for the offsets of the measured pz for the hadrons and therefore, reduce the offets on the MMs. I think this correction will have to be applied for ER and SR and it will be interesting to note its affects. 




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