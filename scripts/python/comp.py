#!/usr/bin/python
from anah10 import *
nentries=10000
#dfs = ''
fold = '/e1f.sim2pi.datadir/comp_old_new_sim/Q2W__1.9-2.5__1.3-1.9/cooked/1.root'
fnew = ''
def runep_nosplitbos_noreconhbook(): #[DC-Match:EC-BAD]
	#ep gpp-pars
	#$EPSIM/gsim.ffread
	#q2w2/user_ana.tcl (matches old-recsis)
	global fnew
	fnew = '/e1f.sim2pi.datadir/sim_range_study/q2w2/cooked/1.root'
	plot("runep_nosplitbos_noreconhbook")
	
def runye(): #[DC-NoMatch:EC-OK]
	#ye gpp-pars
	#$EPSIM/gsim.ffread(diff NOMCDATA)
	#q2w2/user_ana.tcl (matches old-recsis)
	global fnew
	fnew = '/e1f.sim2pi.datadir/comp_old_new_sim/q2w2_ye_gpp/cooked/1.root'
	plot("runye")

def runep_exact(): #[DC-NoMatch:EC-OK]
	#ep gpp-pars
	#$EPSIM/gsim.ffread
	#$EPSIM/user_ana.tcl
	#splitbos;reconhbook
	global fnew
	fnew = '/e1f.sim2pi.datadir/comp_old_new_sim/q2w2-ep-exact/cooked/1.root'
	plot("runep_exact")

def runat(): #[DC-NoMatch:EC-OK]
	#at gpp-pars
	#$EPSIM/gsim.ffread
	#$EPSIM/user_ana.tcl
	#splitbos;reconhbook
	global fnew
	fnew = '/e1f.sim2pi.datadir/comp_old_new_sim/q2w2_old-gpp/cooked/1.root'
	plot("runat")

def run_q2w2_mQ2W(): #[DC-Match:EC-BAD]
	#same as runep_nosplitbos_noreconhbook()
	global fnew
	fnew = '/e1f.sim2pi.datadir/comp_old_new_sim/q2w2-mQ2W/cooked/1.root'
	plot("run_q2w2-mQ2W")

def run_q2w2_mQ2W_spbs(): #[DC-NoMatch:EC-OK]
	#run_q2w2_mQ2W + splitbos
	global fnew
	fnew = '/e1f.sim2pi.datadir/comp_old_new_sim/q2w2-mQ2W-spbs/cooked/1.root'
	plot("run_q2w2-mQ2W_spbs")

#[12-15-13]
def run_nogpp_121513(): #[DC-NoMatch:EC-OK]
	#run_q2w2_mQ2W + splitbos
	global fnew
	fnew = '/e1f.sim2pi.datadir/comp_old_new_sim/nogpp_121513/cooked/1.root'
	plot("run_nogpp_121513")

def run_test3_gpp_121513(): #[DC-NoMatch:EC-OK]
	#run_q2w2_mQ2W + splitbos
	global fnew
	fnew = '/e1f.sim2pi.datadir/comp_old_new_sim/test3-gpp_121513/cooked/1.root'
	plot("run_test3_gpp_121513")

def run_rdcls_gpp_121513(): #[DC-NoMatch:EC-OK]
	#run_q2w2_mQ2W + splitbos
	global fnew
	fnew = '/e1f.sim2pi.datadir/comp_old_new_sim/rdcls-gpp_121513/cooked/1.root'
	plot("run_rdcls_gpp_121513")

def run_test_2pi_VII(): 
	#test run post addition of elast_gen in `subsim` & `runsim`
	global fnew
	fnew = '/e1f.sim2pi.datadir/comp_old_new_sim/test_2pi_VII/cooked/1.root'
	plot("run_test_2pi_VII")

def plot(runtitle):
	dfs = import_h10([fold,fnew],nentries)
	
	fig = plt.figure('q2w_%s'%runtitle,figsize=(8,6))
	fig.suptitle('Reconstructed Q2,W')
	plot_qsq_w(dfs)

	# fig=plt.figure('%s_MAIN'%runtitle,figsize=(10,6))
	# fig.suptitle('Electron data from MAIN')
	# hists(MAIN,dfs)
	
	fig=plt.figure('%s_DC'%runtitle,figsize=(10,6))
	fig.suptitle('Electron data from DC')
	hists(EL_DC,dfs)

	fig=plt.figure('%s_EC'%runtitle,figsize=(10,6))
	fig.suptitle('Electron data from EC')
	hists(EL_EC,dfs)

	fig=plt.figure('%s_SC'%runtitle,figsize=(10,6))
	fig.suptitle('Electron data from SC')
	hists(EL_SC,dfs)