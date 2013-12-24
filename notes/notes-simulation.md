Simulation Notes
==================

[10-10-13]Starting new Simulation run
=====================================
!BEFORE STARTING, Verify following!:
+ Verify why after adding `reconhbook` mcnentr = 499 (out of 500) after Reconstruction.
+ Make sure beam energy for e1f is 5.497 EVERYWHERE (genev.inp, particle_constants, analysisconstants and other scripts)


Potential Causes for:
1. DC mismatch with old-sim
2. EC SF incorrect
-------------------------------------------------------------
+ gpp-pars
+ user_ana.tcl
+ addition of splitbos
+ addition of reconhbook, then reconroot

[12-10-13] q2w2_mQ2W 
---------------------
Should match `ep_nosplitbos_noreconhbook`/q2w2 conditions:
	+ genev
	+ gsim: 
		+ prod
		+ q2w2/gsim.ffread
	+ gpp: 
		+ prod
		+ ep-pars
		+ `splitbos` not used
	+ user_ana:
		+ q2w2/user_ana.tcl (matches old recsis)
		+ `reconhbook` not made i.e. direct ana.hbook --> recon.root

[12-10-13] q2w2_mQ2W-spbs 
---------------------
Same as q2w2_mQ2W, but with
	+ gpp: 
		+ ADDED `splitbos`

[12-11-13]
-----------
So far it seems like that the "new-gpp":
+ is shifting the vertex positions and therefore, affecting DC pars.
+ is overwriting the `Runnum` from GSIM and therefore requires `splitbos`


[12-23-13]Observations when comparing hbook produced by `nt10maker` & `user_ana` (had to do this when setting up Elastic Simulation)
=============================================================

+ Results for this study in comphbook/dbg8_atsim,dbg6_atelas; code macros/comphbook.C
	+ NOTE, by mistake, I continued to use rdcls-GPP pars for these studies!
+ See handwritten notes for full details
+ Implemented for `nt10maker` (see $N10MAKER/REAME.AT):
	+ -t9 = SEB+PART
	+ -t10 = SEB only

For atsim (Or more generally, when Evtgen Banks = MCTK)
-------------------------------------------------------
user_ana(SEB+MCTK+PART)=user_ana(SEB+MCTK)=nt10maker(-t3)
(This is why for 2pi-sim, I could have used either `user_ana` or `nt10maker`)
	

For atelas (Or more generally, when Evtgen Banks = PART)
--------------------------------------------------------
The.hbook made by `user_ana`{SEB+MCTK+PART,SEB+MCTK,SEB+PART,SEB} seem to have fewer entries as compared to the .hbook made by `nt10user`{-t3,-t9,-t10}. The distributions may as well be different.

In conclusion, It seems that the hbook file made by `user_ana` and `nt10maker` differ if along with the booking of SEB banks, the booking of:
	+ MCTK Bank is NOT required (or that MCTK Banks are not present in ANA.OUT)
AND/OR
	+ PART Banks ARE required 

The above is not an issue for 2pi simulations since the Evtgen uses MCTK banks. However, it is a concern for Elastic Events, for which Evtgen uses PART banks. At this point, I am keeping .hbook files from `user_ana` and `nt10maker` and see what compares best with data.

Miscellaneous Notes
===================
+ WALLTIME for 5K 2pi-sim jobs ~ 120 (check log for exact time taken)
+ WALLTIME for 25K elastic-sim jobs ~15. Therefore, 100K jobs = 60 min ?