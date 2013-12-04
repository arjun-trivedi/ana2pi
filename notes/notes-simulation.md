Quick-Notes Simulation
======================

#### [11-28-13] Potential Acceptance Estimation Anomaly 
I think I may be overestimating my Acceptance, for I am using Thrown events from the Cooked .root file made after Reconstruction. __This file does not have all the Thrown events__; hence my Acceptance will be overestimated and therefore, Yields, underestimated. 

However, after moving migrating to the new CentOS 6 systems and starting to use $EPSIM as a basis for my simulation,
I don't know how, but the above problem seems to have disappeared. In my .root file after Reconstruction, I seem to keep all the Thrown Events. I am not sure what happened. Following are my attempts to track down this ussue:

1. ? Difference in __new__ `gsim.ffread` vs. __old__ `gsim_e1f.ffread`: `c NOMCDATA 'ALL' vs. NOMCDATA 'ALL'`
	+ Does not seem to be the case. See subdir/test2_old_gsim.ffread
1. ? Difference in __new__ `user_ana.tcl` vs. __old__ `recsis_e1fsim.tcl` 
	+ Does not seem to be the case. See subdir/test_user_ana.tcl.match_oldrecsis

