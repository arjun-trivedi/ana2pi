notes-integxsec
===============
Here I keep track of how my extraction of the Integrated Cross-section is coming along.

[12-02-13] at-integ-xsec < evgeny-integ-xsec
============================================= 
The story so far is that, although the shape of my Cross-sections look acceptable, they are factor of ~4 smaller than Evgeny's. The factors that could be affecting this, that I have so far identified are (In no particular order):

[ ] 1. Fiducial Cuts for hadrons are not applied  
[ ] 2. Acceptance may currently be overestimated (!LARGEST POSSIBLE SOURCE of DISCREPANCY!)
		+ Since I am using `mctk` banks in the h10-Tree created after Reconstruction. At this stage, the  `mctk' bank retains only those events that have been Reconstructed. Hence, I am not accounting for Acceptance using all the Thrown events.
[X] 3. MM cuts are "not optimized"
		+ I am not sure if this cut over- or under-estimates the Acceptance
		+ [01-31-14] I started looking into the details of this issue in sub-study setup_sim_CentOS6. For details look into the notes for this substudy. For now it seems that, with the current gpp-parameters(1.37,1.37,1.37,0.962), the EC yields may be underestimated by 1%-2% (Top1,3,4 ~ 1%, Top2 ~ 1%)
[ ] 4. Event selection not mindful of Evan's Run Quality check
		+ Since I am using the integrated luminosity as calculated by Evan, I should make sure that I am using ALL events that were within "good Luminosity Blocks" and NOT using events that were not part of LBs used to calculate the integrated Luminosity  




[10-10-13]Analyzing MM cuts
============================
Observations for MM distribtions:
+ mu-exp < mu-sim (mu-exp is closer to the True value)
+ sigma-exp > sigma-sim
+ Radiative tail
	+ Both distributions have a Radiative tail on the right side
+ 3 pion threshold in exp-data:
	+ MM(p,pip/p,pim) shows a feature starting at ~0.25GeV
	  (MM(p,pip/p,pim)= 0.28 GeV if pim,pi0 or pip,pi0 produced at threshold and are undetected)
	+ MM(pip,pim) shows a feature starting at ~1.05 GeV
	  (MM(pip,pim)=1.07 GeV if p,pi0 are produced at threshold)

Factors affecting MM distributions
-----------------------------------
Momenta of e',p',pip,pim affect MM distributions. To check the Reconstruction + "Correction" of the momenta of the particles:
- e': Use Elastic Events 
- p': ?Use Elastic Events
- pim,pip: Use Alignment of MM distributions in Data and Sim.

Note on Momentum Corrections
-----------------------------
Momentum Corrections:
	1. Systematic: DC chambers, B-field mis.-align in EXP
		+ For e': Glen, Gohn
	2. Due to Energy loss
		+ Gleb: proton
		+ Gohn: hadrons
	3. "Kijun's method"
