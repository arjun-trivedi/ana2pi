Simulation Notes
======================

[10-10-13]Starting new Simulation run
=====================================
+ Verify why after adding `reconhbook` mcnentr = 499 (out of 500) after Reconstruction.
+ Make sure beam energy for e1f is 5.497 EVERYWHERE (genev.inp, particle_constants, analysisconstants and other scripts)

Potential Issues for Q2,W distributions not matching old-sim:
+ gpp-pars
+ user_ana.tcl
+ addition of splitbos
+ addition of reconhbook, then reconroot

[10-10-13] q2w2_mQ2W 
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

[10-10-13] q2w2_mQ2W-spbs 
---------------------
Same as q2w2_mQ2W, but with
	+ gpp: 
		+ ADDED `splitbos`

Miscellaneous Notes
===================
+ WALLTIME for 5K sim-jobs ~ 120 (check log for exact time taken)