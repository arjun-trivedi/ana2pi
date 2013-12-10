#!/usr/bin/python
from anah10 import *
nentries=10000
#dfs = ''
fold = '/e1f.sim2pi.datadir/comp_old_new_sim/Q2W__1.9-2.5__1.3-1.9/cooked/1.root'
fnew = ''
def runep_nosplitbos_noreconhbook(): #ep gpp-pars
	#$EPSIM/gsim.ffread
	#q2w2/user_ana.tcl 
	global fnew
	fnew = '/e1f.sim2pi.datadir/sim_range_study/q2w2/cooked/1.root'
	plot("runep_nosplitbos_noreconhbook")
	
def runye(): #ye gpp-pars
	#$EPSIM/gsim.ffread(diff NOMCDATA)
	#q2w2/user_ana.tcl 
	global fnew
	fnew = '/e1f.sim2pi.datadir/comp_old_new_sim/q2w2_ye_gpp/cooked/1.root'
	plot("runye")

def runep_exact(): #ep gpp-pars
	#$EPSIM/gsim.ffread
	#$EPSIM/user_ana.tcl
	#splitbos;reconhbook
	global fnew
	fnew = '/e1f.sim2pi.datadir/comp_old_new_sim/q2w2-ep-exact/cooked/1.root'
	plot("runep_exact")

def runat(): #at gpp-pars
	#$EPSIM/gsim.ffread
	#$EPSIM/user_ana.tcl
	#splitbos;reconhbook
	global fnew
	fnew = '/e1f.sim2pi.datadir/comp_old_new_sim/q2w2_old-gpp/cooked/1.root'
	plot("runat")

def plot(runtitle):
	dfs = import_h10([fold,fnew],nentries)
	fig = plt.figure('q2w_%s'%runtitle,figsize=(8,8))
	fig.suptitle('Reconstructed Q2,W')
	plot_qsq_w(dfs)
	fig=plt.figure('det_%s'%runtitle,figsize=(10,8))
	fig.suptitle('Compare (electron)data from for each sub-detector')
	plot_evt_sub_cols(dfs)


