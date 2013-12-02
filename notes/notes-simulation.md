Quick-Notes Simulation
======================

#### 11-28-13
* I think I may be overestimating my Acceptance, for I am using Thrown events from the Cooked .root file made after Reconstruction. __This file does not have all the Thrown events__; hence my Acceptance will be overestimated and therefore, Yields, underestimated. Is it because: 
	1. ? Difference in __new__ `gsim.ffread` vs. __old__ `gsim_e1f.ffread`: `c NOMCDATA 'ALL' vs. NOMCDATA 'ALL'`
	1. ? Difference in __new__ `user_ana.tcl` vs. __old__ `recsis_e1fsim.tcl` 

