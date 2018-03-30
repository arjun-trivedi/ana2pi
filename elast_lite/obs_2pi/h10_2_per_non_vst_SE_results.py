#!/usr/bin/python
from __future__ import division

import os,sys,datetime,subprocess,time

from collections import OrderedDict
import itertools

#from disp_obs import DispObs #! for SS
import disp_SS as tool_disp_SS#! for plot_SS

#! --  Additional options. For full details see h10looper_e1f.h -- 
#! + Default values on [12-14-15], in brackets, are toggled when any of the options are selected
# 1:  ECin (OFF)
# 2:  ECfid (OFF)
# 3:  zvtx (OFF)
# 4:  etot (ON)
# 5:  hitSC (ON) 
# 6:  dc_stat (ON)
# 7:  ep_efid ( ON, only for E1F)
# 8:  ep_pfid (ON) ([06-12-17] On also for E16. Hitherto was only on for E1F)
# 9:  Q2_var_binw_bng (ON)
# 10: t2 (OFF)
# 11: MM2_cut_EI (OFF)
# 12: reconcile (OFF)
# 13: gpart_pid (OFF)
# 14: hitSC_pid (ON)
# 15: stat_pid (OFF)
# 16: thesis Q2W (ON)
# 17: eff_scpd_at_mod (OFF)
# 18: ECfid_at_mod (OFF) (option #2 i.e. ECfid has to be ON for this to take effect)
# 19: MM2_cut_SS (OFF)
# 20: CC_cut_lse (OFF)
# 21: CC_cut_tgt (OFF)
# 22: CC_cut_pmtC_av (OFF) (either of option #20 or #21 have to be ON for this to take effect)
# 23: CC_cut_pmtC_L (OFF)  (either of option #20 or #21 have to be ON for this to take effect)
# 24: CC_cut_pmtC_R (OFF)  (either of option #20 or #21 have to be ON for this to take effect) 
# 25: cut_zvtx_etgt_bg_sub (ON, only for E16) (option #3 has to be on for this to take effect) (added on [06-12-17])
# ----

'''
[09-16-15]
+ This script provides a solution for making 'per_non_vst_SE_result' for each SE
+ The solution is hackish and not elegant (a more elegant solution is in development). 
+ This script is basically a copy of 'h10_2_Obs' with the following modifications:
	1. BASIC_ADTNL_OPTS changed: 
		+ ':1:3:11:17:' -> ':1:3:17:2:18:20:'
		  (ECin-on:zvtx-on:MM2_cut_EI-on:eff_scpd_atmod ->ECin-on:zvtx-on:eff_scpd_atmod:ECfid-ON:ECfid_atmod_ON:)
		+ This change is to update the BASIC_ADTNL_OPTS to include 'ECfid-ON:ECfid_atmod_ON:'(:2:18:) and ':CC_cut_lse:' (:20:)
		+ Note that ':11:' (':MM2_cut_EI-on:') is removed. However, this is purely due to technical reasons
		  because in the MM SE variations I need to turn it off so that 'MM2_cut_AT' is used 
		  and this I can only do by removing this option from BASIC_ADTNL_OPTS and toggling it ADTNL_OPTS_CODE_LIST
	2. This script now takes in an additional argument called SYSTEMATIC_EFFECT, based on which ADTNL_OPTS_CODE_LIST and ADTNL_OPTS_TAG_LIST
	   are setup.
	3. As a result of 2., CUTNSCORS are now setup after getting all user arguments
	4. For Obs_R2, using both True and False options for  EXTRACT_D_E_FOR_NON_ALPHA
	5. As a result of above, clean up of historical archive of the evolution of CUTSNCORS and 
	   its replacement by relevant documentation
	6. Removed code for, the options of which were obtained from user, for:
		i.   PROC_H10_LITE_TYPE: always run in parallel (no need found to run serially)
		ii.  MAKE_SS_PLOTS_ONLY: after new SE procedure, this is obsolete
		iii. SKIP_SS_PLOTS: after new SE procedure, this is obsolete


+ The output of this script, which orchestrates h10->Observables, is organized by jobtag=MAINTAG/cutsncors-cutsncorsnum, where
	+ MAINTAG=IDTFR_DATE. This serves, as the name implies, as the main tag for the user to tag the results from this script. Under this all the various processing of h10->Observables, starting with the  CUTNSCORS applied in the process {ER,SR}h10->d2pi, are organized. 
	+ cutsncorsnum = the numeric key in the dictionary CUTNCORS which contains all the cutsncors applied to {ER,SR}h10->d2pi
+ In this sense, the jobtag is only applicable for {ER,SR}h10->d2pi; {ST}h10->d2pi is independent of 'jobtag' and is directly made on farm labelled by DATE_ST
	+ To use the appropriate d2piT_<date>.root, the DATE_ST(EXPT,sim) is used.
	
+ Usage
	>h10_2_per_non_vst_SE_results expt<=e1f/e16> Q2min Q2max idtfr cumsiml[=SIMS[EXPT]] systematic_effect date_of_existing_jobtag[=today's MMDDYY] skip_make_d2pi[=False] skip_dobs_R2[=False]
	+ idtfr: string used to describe the job; can be entered as ""
	+ cumsiml: output created for every cumsim in cumsiml under jobtag/cumsim
	+ date_of_existing_jobtag: date should be of existing jobtag since the program assumes that jobtag is going to be updated with, as set up currently in the code,for example, addition of a new simulation.
'''
USAGE="h10_2_per_non_vst_SE_results expt<=e1f/e16> Q2min Q2max idtfr cumsiml[=SIMS[EXPT]] systematic_effect date_of_existing_jobtag[=today's MMDDYY] skip_make_d2pi[=False] skip_dobs_R2[=False]"

#! *** Setup global variables ***
EXPTS=['e1f','e16']
SEQS=['ER','SR','ST']
SIMS={'e1f':['TBD'],'e16':['sim4','sim5','sim6','sim7','sim8','sim9','sim10','sim11','sim12','sim13','sim14','sim15']}
#print EXPTS,SEQS,SIMS

H10PATH='/home/trivedia/ongoing/h10lsts_local'
H10LSTS={
('e1f','ER'):'TBD',
('e16','ER'):'%s/e16/exp/h10_2_h10-skim-SS_060717_lse_071017.lst'%(H10PATH), #! h10_2_h10-skim-SS_060717_lse_061217.lst, h10_2_h10-skim-SS_081716.lst, h10-skim-SS_030116.lst,h10-skim-e_010816.lst, h10-skim-e_122115.lst w/o ER:SF fix
('e1f','SR','sim-TBD'):'TBD',
('e16','SR','sim4'):  '%s/e16/sim/sim4/h10_2_h10-skim-SS_060717.lst'%(H10PATH), #! h10_2_h10-skim-SS_081716.lst,h10_2_h10-skim-SS_030116.lst,h10-skim-e_122115.lst
('e16','SR','sim5'):  '%s/e16/sim/sim5/h10_2_h10-skim-SS_060717.lst'%(H10PATH), #! h10_2_h10-skim-SS_081716.lst,h10_2_h10-skim-SS_030116.lst,h10-skim-e_122115.lst
('e16','SR','sim6'):  '%s/e16/sim/sim6/h10_2_h10-skim-SS_060717.lst'%(H10PATH), #! h10_2_h10-skim-SS_081716.lst,h10_2_h10-skim-SS_030116.lst,h10-skim-e_122915.lst
('e16','SR','sim7'):  '%s/e16/sim/sim7/h10_2_h10-skim-SS_060717.lst'%(H10PATH), #! h10_2_h10-skim-SS_081716.lst,h10_2_h10-skim-SS_030116.lst,h10-skim-e_010816.lst 
('e16','SR','sim8'):  '%s/e16/sim/sim8/h10_2_h10-skim-SS_060717.lst'%(H10PATH), #! h10_2_h10-skim-SS_081716.lst,h10_2_h10-skim-SS_030116.lst,h10-skim-e_012416.lst
('e16','SR','sim9'):  '%s/e16/sim/sim9/h10_2_h10-skim-SS_060717.lst'%(H10PATH), #! h10_2_h10-skim-SS_081716.lst,h10_2_h10-skim-SS_030116.lst,h10-skim-e_012416.lst
('e16','SR','sim10'):'%s/e16/sim/sim10/h10_2_h10-skim-SS_060717.lst'%(H10PATH), #! h10_2_h10-skim-SS_081716.lst,h10_2_h10-skim-SS_031516.lst
('e16','SR','sim11'):'%s/e16/sim/sim11/h10_2_h10-skim-SS_060717.lst'%(H10PATH), #! h10_2_h10-skim-SS_081716.lst,h10_2_h10-skim-SS_031516.lst
('e16','SR','sim12'):'%s/e16/sim/sim12/h10_2_h10-skim-SS_060717.lst'%(H10PATH), #! h10_2_h10-skim-SS_081716.lst,h10_2_h10-skim-SS_031516.lst
('e16','SR','sim13'):'%s/e16/sim/sim13/h10_2_h10-skim-SS_060717.lst'%(H10PATH), #! h10_2_h10-skim-SS_081716.lst,h10_2_h10-skim-SS_031516.lst
('e16','SR','sim14'):'%s/e16/sim/sim14/h10_2_h10-skim-SS_032318.lst'%(H10PATH),
('e16','SR','sim15'):'%s/e16/sim/sim15/h10_2_h10-skim-SS_032318.lst'%(H10PATH)
}
#print H10LSTS
DATE_ST={
('e1f','sim-TBD'):'TBD',
('e16','sim4') :'122215',
('e16','sim5') :'122215',
('e16','sim6') :'122915',
('e16','sim7') :'010816',
('e16','sim8') :'012416',
('e16','sim9') :'012416',
('e16','sim10'):'032116',
('e16','sim11'):'032116',
('e16','sim12'):'032116',
('e16','sim13'):'032116',
('e16','sim14'):'032318',
('e16','sim15'):'032318'
}

#! *** Define WMIN,WMAX ***
#! While Q2MIN,Q2MAX are input from the user, the main reason being the need to select between the two 
#! sets of Simulations and therefore related two sets of Observables, aside from the need to able to selecting 
#! bins in Q2 for debugging purposes, the W range in that sense is fixed.
WMIN='1.400'
WMAX='2.125'
#! *** End of global variables ***

#! *** Get input arguments from user ***
if len(sys.argv)<2:
	sys.exit("Please enter expt as per usage: %s"%USAGE)
EXPT=sys.argv[1]
if EXPT!="e16" and EXPT!="e1f":
        sys.exit("Please enter EXPT as either e16 or e1f")

Q2MIN=""
if len(sys.argv)>2:#i.e. Q2MIN entered by user
        Q2MIN=sys.argv[2]
else:
	sys.exit("Please enter Q2MIN as per usage: %s"%USAGE)

Q2MAX=""
if len(sys.argv)>3:#i.e. Q2MAX entered by user
        Q2MAX=sys.argv[3]
else:
        sys.exit("Please enter Q2MAX as per usage: %s"%USAGE)

IDTFR=""
if len(sys.argv)>4:#i.e. IDTFR entered by user
	IDTFR=sys.argv[4]
else:
	sys.exit("Please enter IDTFR as per usage: %s"%USAGE)

CUMSIML=[]
if len(sys.argv)>5:#i.e. CUMSIML entered by user
	cumsims=sys.argv[5].split("_")
	if len(cumsims)>0:
		for cumsim in cumsims:
			cumsim=cumsim.split(':')
			#print cumsim
			#cumsim=map(int,cumsim)
			for sim_user in cumsim:
				sim_user_valid=False
				for sim in SIMS[EXPT]:
					#print "sim_user=",sim_user
					#print "sim=",sim
					if sim_user==sim:
						sim_user_valid=True
				if not sim_user_valid:
					sys.exit("Please enter valid simX. %s for %s is not valid."%(sim_user,EXPT))
			CUMSIML.append(cumsim)
	else:
		sys.exit("Please enter cumsiml as sim1_sim1:sim2_sim1:sim2:sim3")
else:
	CUMSIML.append(SIMS[EXPT])

if len(sys.argv)>6:#i.e. SYSTEMATIC_EFFECT entered by user
 	SYSTEMATIC_EFFECT=sys.argv[6]
else:
	sys.exit("Please enter SYSTEMATIC_EFFECT as per usage: %s"%USAGE)


PROC_NEW_JOBTAG=True
if len(sys.argv)>7: #! i.e. date of existing jobtag entered by user i.e. update existing jobtag
	DATE=sys.argv[7]
        PROC_NEW_JOBTAG=False
else:
        DATE=datetime.datetime.now().strftime('%m%d%y')

if len(sys.argv)>8: #! i.e. skip_make_d2pi information entered by user
	if sys.argv[8]=="False":
		SKIP_MAKE_D2PI=False
	elif sys.argv[8]=="True":
		SKIP_MAKE_D2PI=True
	else:
		print "Please enter skip_make_d2pi=True or False only."
		sys.exit()
else:
	SKIP_MAKE_D2PI=False

if len(sys.argv)>9: #! i.e. skip_dobs_R2 information entered by user
	if sys.argv[9]=="False":
		SKIP_DOBS_R2=False
	elif sys.argv[9]=="True":
		SKIP_DOBS_R2=True
	else:
		print "Please enter skip_dobs_R2=True or False only."
		sys.exit()
else:
	SKIP_DOBS_R2=False

#! ******

#! *** Set up some other global variables after getting user input *** 

#! *** Setup dictionary for CUTSNCORS *** 
#! + The following CUTSNCORS include all cuts and corrections that affect the final cross-sections in ways that I have not fully understood yet.
#!   (For a historical development of CUTSNCORS upto this point see $ELAST_LITE/h10_2_Obs)
#! + syntax: CUTSNCORS={num:[procorder,adtnl_opts,adtnl_opts_tag]}
CUTSNCORS=OrderedDict()

#! 1. First set up the Basic procorder and adtnl_opts that are common to all 16 CUTSNCORS
BASIC_PROCORDER='eid:efid:pid:pfid:pcorr:eff:scpd:evtsel_2pi:'
BASIC_ADTNL_OPTS=':1:3:17:2:18:20:' #! ECin-on:zvtx-on:eff_scpd_atmod:ECfid-ON:ECfid_atmod_ON:CC_cut_lse:

#! 2. Now as per SYSTEMATIC_EFFECT
if SYSTEMATIC_EFFECT=="SSBands":
	ADTNL_OPTS_CODE_LIST=[':11:',                        ':11:13:15:']
	ADTNL_OPTS_TAG_LIST =[':MM-EI:gpart-pid-OFF:stat-pid-OFF:',':MM-EI:gpart-pid-ON:stat-pid-ON:']
elif SYSTEMATIC_EFFECT=="MM":
	ADTNL_OPTS_CODE_LIST=['',                                 ':11:']
	ADTNL_OPTS_TAG_LIST =[':MM-AT:gpart-pid-OFF:stat-pid-OFF:',':MM-EI:gpart-pid-OFF:stat-pid-OFF:']
elif SYSTEMATIC_EFFECT=="SSBands-off-off": #! [10-26-17] Added for testing observables as a function of simstats: used by $SUBSTUDIES/study_compare_AT_EI_results/make_obs_as_function_of_simstats.sh
	ADTNL_OPTS_CODE_LIST=[':11:']
	ADTNL_OPTS_TAG_LIST =[':MM-EI:gpart-pid-OFF:stat-pid-OFF:']
else:
	sys.exit("Please enter SYSTEMATIC_EFFECT as either SSBands or MM")


#! 3. Now finally, using 1. and 2., put together the CUTSNCORS
for i,lst in enumerate(zip(ADTNL_OPTS_CODE_LIST,ADTNL_OPTS_TAG_LIST)):
    adtnl_opt_code=BASIC_ADTNL_OPTS+''.join(lst[0])
    adtnl_opt_tag=''.join(lst[1])
    CUTSNCORS[i+1]=[BASIC_PROCORDER,adtnl_opt_code,adtnl_opt_tag]
if len(CUTSNCORS)==0: sys.exit("len(CUTSNCORS)==0. Exiting")
#sys.exit()
#! *** End: Define the CUTSNCORS ***

MAINTAG="%s_%s"%(IDTFR,DATE)

#! Create LOGDIR for the h10_2_per_non_vst_SE_results job ***
LOGDIR="/home/trivedia/h10_2_per_non_vst_SE_results_logs"
if not os.path.exists(LOGDIR):
	os.makedirs(LOGDIR)

#! Redirect stdout of this orchestrating script to LOGDIR_MAINTAG(=LOGDIR/MAINTAG/h10_2_per_non_vst_SE_results.log)
#! + All logs from various jobtags under this MAINTAG will be placed in this folder
LOGDIR_MAINTAG="%s/%s"%(LOGDIR,MAINTAG)
if not os.path.exists(LOGDIR_MAINTAG):
        os.makedirs(LOGDIR_MAINTAG)
if PROC_NEW_JOBTAG:
	sys.stdout=open('%s/h10_2_per_non_vst_SE_results.log'%(LOGDIR_MAINTAG), 'w')
else:
	sys.stdout=open('%s/h10_2_per_non_vst_SE_results_update_%s.log'%(LOGDIR_MAINTAG,datetime.datetime.now().strftime('%m%d%y_%H%M%S')), 'w')

#! Put all relevant global variables information in the file
print "h10_2_per_non_vst_SE_results start time",datetime.datetime.now().strftime('%H:%M:%S')

print "Total number of CUTSNCORS=",len(CUTSNCORS)
print "CUTSNCORS pretty print:"
for k in CUTSNCORS:
        print k,CUTSNCORS[k]

print "EXPT=",EXPT
print "Q2MIN=",Q2MIN
print "Q2MAX=",Q2MAX
print "IDTFR=",IDTFR
print "CUMSIML=",CUMSIML
print "SYSTEMATIC_EFFECT=",SYSTEMATIC_EFFECT
print "DATE=",DATE
print "PROC_NEW_JOBTAG=",PROC_NEW_JOBTAG
print "SKIP_MAKE_D2PI=",SKIP_MAKE_D2PI
print "SKIP_DOBS_R2=",SKIP_DOBS_R2
#sys.exit()
#! ******

#! *** Now start processing CUTSNCORS ***
for k in CUTSNCORS.keys():
	#continue #! debug
	#if k>2: continue #! debug
    	procorder,adtnl_opts,adtnl_opts_tag= CUTSNCORS[k][0],CUTSNCORS[k][1],CUTSNCORS[k][2]
	#! For SR, remove ':pcorr' from procorder 
	procorder_SR=procorder.replace(':pcorr','')
	#! Prepare jobtag
	cutsncorsnum=k
	jobtag='%s/cutsncors%d'%(MAINTAG,cutsncorsnum)

	logdir=""
	if PROC_NEW_JOBTAG:
		logdir="%s/%s"%(LOGDIR,jobtag)
	else:
		logdir="%s/%s_update_%s"%(LOGDIR,jobtag,datetime.datetime.now().strftime('%m%d%y_%H%M%S'))
	if not os.path.exists(logdir):
		os.makedirs(logdir)
	#! + While I can re-direct stdout/err of all subprocesses to logdir/<log> using Popen,
	#!   I have to use a hack-ish method to redirect stdout (mainly from print messages)
	#!   from the main program to logdir/main.log
	#! + Using this method, the main.log is not updated in real-time but after the file is
	#!   closed at the end of the program! Therefore, 'print' continue to accompany 'mainlog.write'
	#!   to see real-time progress
	mainlog=open("%s/main.log"%(logdir),'w')
	print "All stdout from h10_2_per_non_vst_SE_results for %s are in %s/%s"%(jobtag,logdir,"main.log")

	print "h10->Obs for jobtag=%s"%jobtag
        print "h10->Obs for expt:procorder:adtnl_opts:adtnl_opts_tag=(%s:%s:%s:%s)"%(EXPT,procorder,adtnl_opts,adtnl_opts_tag)
	mainlog.write("h10->Obs for jobtag=%s"%jobtag)
	mainlog.write("h10->Obs for expt:procorder:adtnl_opts:adtnl_opts_tag=(%s:%s:%s:%s)\n"%(EXPT,procorder,adtnl_opts,adtnl_opts_tag))
		
	#! 1. {ER,SR,ST}:{h10-skim-e->d2pi,h10->d2pi}
	print "*** 1. {ER,SR,ST}:{h10-skim-e->d2pi,h10->d2pi} ***"
	mainlog.write("*** 1. {ER,SR,ST}:{h10-skim-e->d2pi,h10->d2pi} ***\n")
	#! ER
	print "*** ER *** "
	mainlog.write("*** ER ***\n")
	sfx=''
	if EXPT=='e16':sfx='_E16'
	d2pidir_exp=os.path.join(os.environ['D2PIDIR_EXP%s'%sfx],jobtag)
	print 'd2pidir_exp=',d2pidir_exp
	mainlog.write('d2pidir_exp=%s\n'%d2pidir_exp)
	if not os.path.exists(d2pidir_exp):
		os.makedirs(d2pidir_exp)
	#! Process d2piR.root only if it does not exist OR if processing a non-existing jobtag 
	if (not SKIP_MAKE_D2PI) and (not os.path.isfile(os.path.join(d2pidir_exp,"d2piR.root")) or PROC_NEW_JOBTAG):
		h10lst=H10LSTS[EXPT,'ER']
		print 'h10lst=',h10lst
		mainlog.write('h10lst=%s\n'%h10lst)
		cmd=["proc_h10_lite","-i","%s"%h10lst,"-t","%s:exp:2pi:recon"%EXPT,"-c","%s"%procorder, "-o","%s/d2piR.root"%d2pidir_exp, "%s"%adtnl_opts]
		#! For debug
		#cmd=["proc_h10_lite","-i","%s"%h10lst,"-t","%s:exp:2pi:recon"%EXPT,"-c","%s"%procorder, "-o","%s/d2piR.root"%d2pidir_exp, "-n","10000","%s"%adtnl_opts]
        	print">>>%s >& %s/h10_2_d2pi_ER.log.log\n"%(cmd,logdir)
		mainlog.write(">>>%s >& %s/h10_2_d2pi_ER.log.log\n"%(cmd,logdir))
		logfile=open('%s/h10_2_d2pi_ER.log'%logdir,'w')
		#! + Popen called below without keeping any references to the Popen objects
		#!   (http://stackoverflow.com/questions/2760652/how-to-kill-or-avoid-zombie-processes-with-subprocess-module: ibz answered Aug 15 '11 at 9:53)
		#! + This way no zombie child processes are left and this is important
		#!   for the manner in which I am tracking and noting the completion of
		#!   ER,SR{sims}}:h10->d2pi processes below
		subprocess.Popen(cmd,stdout=logfile,stderr=subprocess.STDOUT)
	else:
		if SKIP_MAKE_D2PI:
			print "SKIP_MAKE_D2PI=True. Therefore not making one."
			mainlog.write("SKIP_MAKE_D2PI=True. Therefore not making one.")
		else:
			print os.path.join(d2pidir_exp,"d2piR.root"),"exists. Not making new one"
			mainlog.write("%s exists. Not making new one\n"%os.path.join(d2pidir_exp,"d2piR.root"))	

	#! SR
	print "SR"
	mainlog.write("SR\n")
	for CUMSIM in CUMSIML:
		print "*** SR for CUMSIM=%s***"%CUMSIM
		mainlog.write("*** SR for CUMSIM=%s ****\n"%CUMSIM)
		for i,sim in enumerate(CUMSIM): #enumerate(SIMS[EXPT]):
			sfx=''
        		if EXPT=='e16':sfx='_E16'
        		d2pidir_sim=os.path.join(os.environ['D2PIDIR_SIM%s'%sfx],sim,jobtag)
			if not os.path.exists(d2pidir_sim):
                		os.makedirs(d2pidir_sim)
			#! Process d2piR.root only if it does not exist OR if processing a non-existing jobtag 
			if (not SKIP_MAKE_D2PI) and (not os.path.isfile(os.path.join(d2pidir_sim,"d2piR.root")) or PROC_NEW_JOBTAG):
				h10lst=H10LSTS[EXPT,'SR',sim] 
				print 'sim=',sim
				mainlog.write('sim=%s\n'%sim)
				print 'h10lst=',h10lst
				print 'd2pidir_sim=',d2pidir_sim
				mainlog.write('h10lst=%s\n'%h10lst)
				mainlog.write('d2pidir_sim=%s\n'%d2pidir_sim)
				cmd=["proc_h10_lite","-i", "%s"%h10lst,"-t","%s:sim:2pi:recon"%EXPT,"-c","%s"%procorder_SR, "-o","%s/d2piR.root"%d2pidir_sim, "%s"%adtnl_opts]
				#! For debug
				#cmd=["proc_h10_lite","-i", "%s"%h10lst,"-t","%s:sim:2pi:recon"%EXPT,"-c","%s"%procorder_SR, "-o","%s/d2piR.root"%d2pidir_sim, "-n","10000","%s"%adtnl_opts]
				print">>>%s >& %s/h10_2_d2pi_SR_%s.log.log\n"%(cmd,logdir,sim)
				mainlog.write(">>>%s >& %s/h10_2_d2pi_SR_%s.log.log\n"%(cmd,logdir,sim))
				logfile=open('%s/h10_2_d2pi_SR_%s.log'%(logdir,sim),'w')
				#! + Popen called below without keeping any references to the Popen objects
        			#!   (http://stackoverflow.com/questions/2760652/how-to-kill-or-avoid-zombie-processes-with-subprocess-module: ibz answered Aug 15 '11 at 9:53)
        			#! + This way no zombie child processes are left and this is important
        			#!   for the manner in which I am tracking and noting the completion of
        			#!   ER,SR{sims}}:h10->d2pi processes below
				subprocess.Popen(cmd,stdout=logfile,stderr=subprocess.STDOUT)
			else:
				if SKIP_MAKE_D2PI:
					print "SKIP_MAKE_D2PI=True. Therefore not making one."
					mainlog.write("SKIP_MAKE_D2PI=True. Therefore not making one.")
				else:
					print os.path.join(d2pidir_sim,"d2piR.root"),"exists. Not making new one"
					mainlog.write("%s exists. Not making new one\n"%os.path.join(d2pidir_sim,"d2piR.root"))

			#! ST
			#! + For each sim,should directly already have been made on the farm on DATE_ST

	#! Wait for h10->d2pi for ER,SR{sims} to finish before proceeding
	while True:
		#! Following command taken from $HOME/ongoing/bin/track_mem_dR2
		process=subprocess.Popen("ps aux | grep \"[p]roc_h10_lite\" | awk '{print $4}'",shell=True,stdout=subprocess.PIPE)
		out,err=process.communicate()
		#print out
		if out=="":
			break
		else:
			print "proc_h10_lite(s) still running"
			mainlog.write("proc_h10_lite(s) still running\n")
			crnt_time=datetime.datetime.now().strftime('%H:%M:%S')
			print "%MEM usage (checked every 5min) of each running instance at",crnt_time,"(HH:MM:SS):"
			mainlog.write("MEM percentage usage (checked every 15min) of each running instance at %s(HH:MM:SS)\n"%crnt_time)
			print out
			mainlog.write("%s\n"%out)
			time.sleep(5*60)
			#! For debug
			#time.sleep(60)

		
	#! Add d2piR/T from sims
	#! NOTE:
	#! + d2pidir_SR_total is dependent on jobgtag and therefore on CUTSNCORS applied to {ER,SR}h10->d2pi
	#! 	+ d2pidir_SR_total=$D2PIDIR_SIM(_E16)/cumsim/jobtag
	#! + d2pidir_ST_total is *independent* of jobgtag and {ST}h10->d2pi is directly made on farm and labelled by dates in DATE_ST
	#! 	+ d2pidir_ST_total=$D2PIDIR_SIM(_E16)/cumsim/T_cumsim1-<DATE_ST(EXPT,cumsim1)>_cumsim2-<DATE_ST(EXPT,cumsim2)>_cumsimX-<DATE_ST(EXPT,cumsimX)>
	for CUMSIM in CUMSIML:
		print "*** hadd d2piR.root, d2piT.root for cumsim=%s"%CUMSIM
		mainlog.write("*** hadd d2piR.root, d2piT.root for cumsim=%s ***\n"%CUMSIM)
		sim_total="_".join(CUMSIM) #"_".join(SIMS[EXPT])
		sfx=''
		if EXPT=='e16':sfx='_E16'
		#! d2pidir_SR_total
		d2pidir_SR_total=os.path.join(os.environ['D2PIDIR_SIM%s'%sfx],sim_total,jobtag)
		#! d2pidir_ST_total
		d2pidir_ST_total=os.path.join(os.environ['D2PIDIR_SIM%s'%sfx],sim_total,'T')
		d2pidir_ST_total_sfx=''
		for i,sim in enumerate(CUMSIM):#enumerate(SIMS[EXPT]):
    			#print sim,DATE_ST[EXPT,sim]
    			d2pidir_ST_total_sfx+='_%s-%s'%(sim,DATE_ST[EXPT,sim])	
		d2pidir_ST_total+=d2pidir_ST_total_sfx
		#! create d2pidir_SR(ST)_total
		for drctry in [d2pidir_SR_total,d2pidir_ST_total]:
			if not os.path.exists(drctry):
				os.makedirs(drctry)	
		tgtfl=['d2piR.root','d2piT.root']
		for f in tgtfl:
			#! First determine if adding d2piR or d2piT
			addSR,addST=False,False
			logfilename=""
			if   "R" in f: 
				addSR=True
				hadd_outdir=d2pidir_SR_total
				logfilename="hadd_SR_sim_%s.log"%(sim_total)
			elif "T" in f: 
				addST=True
				hadd_outdir=d2pidir_ST_total
				logfilename="hadd_ST_sim_%s.log"%(sim_total)
			#! Now hadd
			if ( (addSR and len(CUMSIM)>1 and not SKIP_MAKE_D2PI) and (not os.path.isfile(os.path.join(hadd_outdir,f)) or PROC_NEW_JOBTAG) ) or (addST and not  os.path.isfile(os.path.join(hadd_outdir,f))):
				logfile=open('%s/%s'%(logdir,logfilename),'w')
				tgtf=os.path.join(hadd_outdir,f)
 				srcfl=[]
    				for sim in (CUMSIM):#(SIMS[EXPT]):
					if addSR:
						srcfl.append(os.path.join(os.environ['D2PIDIR_SIM%s'%sfx],sim,jobtag,'d2piR.root'))
					elif addST:
						srcfl.append(os.path.join(os.environ['D2PIDIR_SIM%s'%sfx],sim,'d2piT_%s.root'%DATE_ST[EXPT,sim]))
    				srcfl=sorted(srcfl)
    				cmd=["hadd","-f",tgtf] #! hadd with -f since I often have to remake the same tgtf
    				cmd+=srcfl
                		print">>>%s >& %s/%s\n"%(cmd,logdir,logfilename)
				mainlog.write(">>>%s >& %s/%s\n"%(cmd,logdir,logfilename))
				tool=subprocess.Popen(cmd,stdout=logfile,stderr=subprocess.STDOUT)
				tool.wait()
			else:
				if addSR and len(CUMSIM)==1:
					print "Not going to hadd since len(CUMSIM)==1 and therefore already %s/%s exists"%(hadd_outdir,f)
					mainlog.write("Not going to hadd since len(CUMSIM)==1 and therefore %s/%s exists\n"%(hadd_outdir,f))
				elif addSR and SKIP_MAKE_D2PI:
					print "Not going to hadd since SKIP_MAKE_D2PI=True"          
				else:
					print "Not going to hadd since %s/%s exists"%(hadd_outdir,f)
					mainlog.write("Not going to hadd since %s/%s exists\n"%(hadd_outdir,f))

		#! 2. d2pi->Observables
		print "*** 2. d2pi->Observables for %s***"%CUMSIM
		mainlog.write("*** 2. d2pi->Observables for %s ***\n"%CUMSIM)
		#! Create and prepare obsdir
		#logfile=open("%s/prep_obsdir_%s.log"%(logdir,sim_total),'w')
		sfx=''
		if EXPT=='e16':sfx='_E16'
		obsdir=os.path.join(os.environ['OBSDIR%s'%sfx],jobtag)
		print "obsdir=",obsdir
		mainlog.write("obsdir=%s\n"%obsdir)
		if not os.path.exists(obsdir):
			os.makedirs(obsdir)
		if not os.path.exists(os.path.join(obsdir,'d2pi_exp')):
			os.makedirs(os.path.join(obsdir,'d2pi_exp'))
		if not os.path.exists(os.path.join(obsdir,'d2pi_sim',sim_total)):
			os.makedirs(os.path.join(obsdir,'d2pi_sim',sim_total))

		#! Link files into obsdir/d2pi_exp(_sim)
		if not os.path.isfile(os.path.join(obsdir,"d2pi_exp/d2piR.root")) or PROC_NEW_JOBTAG:
			logfile=open("%s/prep_obsdir_%s.log"%(logdir,sim_total),'a')
			#! First unlink any existing link
			cmd=["unlink","%s/d2pi_exp/d2piR.root"%(obsdir)]
			print">>>%s\n"%cmd
			mainlog.write(">>>%s\n"%cmd)
			tool=subprocess.Popen(cmd,stdout=logfile,stderr=subprocess.STDOUT)
			tool.wait()
			#! Recreate link
			cmd=["ln", "-s", "%s/d2piR.root"%d2pidir_exp,"%s/d2pi_exp/d2piR.root"%(obsdir)]
			print">>>%s\n"%cmd
			mainlog.write(">>>%s\n"%cmd)
			tool=subprocess.Popen(cmd,stdout=logfile,stderr=subprocess.STDOUT)
			tool.wait()
		else:
			print os.path.join(obsdir,"d2pi_exp/d2piR.root"),"exists. Not linking new one"
			mainlog.write("%s exists. Not making new one\n"%os.path.join(obsdir,"d2pi_exp/d2piR.root"))
	
		if not os.path.isfile(os.path.join(obsdir,"d2pi_sim",sim_total,"d2piR.root")) or PROC_NEW_JOBTAG:
			logfile=open("%s/prep_obsdir_%s.log"%(logdir,sim_total),'a')
			#! First unlink any existing link
			cmd=["unlink","%s/d2pi_sim/%s/d2piR.root"%(obsdir,sim_total)]
			print">>>%s\n"%cmd
			mainlog.write(">>>%s\n"%cmd)
			tool=subprocess.Popen(cmd,stdout=logfile,stderr=subprocess.STDOUT)
			tool.wait()
			#! Recreate link
			cmd=(["ln", "-s", "%s/d2piR.root"%d2pidir_SR_total,"%s/d2pi_sim/%s/d2piR.root"%(obsdir,sim_total)])
			print">>>%s\n"%cmd
			mainlog.write(">>>%s\n"%cmd)
			tool=subprocess.Popen(cmd,stdout=logfile,stderr=subprocess.STDOUT)
			tool.wait()
		else:
			print os.path.join(obsdir,"d2pi_sim",sim_total,"d2piR.root"),"exists. Not linking new one"
			mainlog.write("%s exists. Not making new one\n"%os.path.join(obsdir,"d2pi_sim",sim_total,"d2piR.root"))

		if not os.path.isfile(os.path.join(obsdir,"d2pi_sim",sim_total,"d2piT.root")) or PROC_NEW_JOBTAG:
			logfile=open("%s/prep_obsdir_%s.log"%(logdir,sim_total),'a')
			#! First unlink any existing link
			cmd=["unlink","%s/d2pi_sim/%s/d2piT.root"%(obsdir,sim_total)]
			print">>>%s\n"%cmd
			mainlog.write(">>>%s\n"%cmd)
			tool=subprocess.Popen(cmd,stdout=logfile,stderr=subprocess.STDOUT)
			tool.wait()
 			#! Recreate link
			cmd=(["ln", "-s", "%s/d2piT.root"%(d2pidir_ST_total),"%s/d2pi_sim/%s/d2piT.root"%(obsdir,sim_total)])
			print">>>%s\n"%cmd
			mainlog.write(">>>%s\n"%cmd)
			tool=subprocess.Popen(cmd,stdout=logfile,stderr=subprocess.STDOUT)
			tool.wait()
		else:
			print os.path.join(obsdir,"d2pi_sim",sim_total,"d2piT.root"),"exists. Not linking new one"
			mainlog.write("%s exists. Not making new one\n"%os.path.join(obsdir,"d2pi_sim",sim_total,"d2piT.root"))

		#-- Now start to make Obervables
		#echo ">ph8 $obsdir $sim 2.00 3.00 1.400 2.125 >& $logdir_obs/ph8.log "
		if not os.path.isfile(os.path.join(obsdir,sim_total,"yield.root")) or PROC_NEW_JOBTAG:
			cmd=["ph8",obsdir,sim_total,Q2MIN,Q2MAX,WMIN,WMAX]
			logfile=open("%s/ph8_%s.log"%(logdir,sim_total),'w')
			print">>>%s >& %s/ph8_%s.log\n"%(cmd,logdir,sim_total)
			mainlog.write(">>>%s >& %s/ph8_%s.log\n"%(cmd,logdir,sim_total))
			tool=subprocess.Popen(cmd,stdout=logfile,stderr=subprocess.STDOUT)
			tool.wait();
		else:
			print os.path.join(obsdir,sim_total,"yield.root"),"exists. Not making new one"
                	mainlog.write("%s exists. Not making new one\n"%os.path.join(obsdir,sim_total,"yield.root")) 

		#echo ">dobs_1D $obsdir $sim norm 2.00 3.00 1.400 2.125 e16 "
		if not os.path.exists(os.path.join(obsdir,sim_total,"Obs_1D_norm")) or PROC_NEW_JOBTAG:
			cmd=["dobs_1D",obsdir,sim_total,"norm",Q2MIN,Q2MAX,WMIN,WMAX,EXPT]
			logfile=open("%s/dobs_1D_%s.log"%(logdir,sim_total),'w')
			print">>>%s >& %s/dobs_1D_%s.log\n"%(cmd,logdir,sim_total)
			mainlog.write(">>>%s >& %s/dobs_1D_%s.log\n"%(cmd,logdir,sim_total))
			tool=subprocess.Popen(cmd,stdout=logfile,stderr=subprocess.STDOUT)
			tool.wait()
		else:
			print os.path.join(obsdir,sim_total,"Obs_1D_norm"),"exists. Not making new one"
			mainlog.write("%s exists. Not making new one\n"%os.path.join(obsdir,sim_total,"Obs_1D_norm"))

        #! dobs_R2
        #! echo "dobs_R2 mthd1:mthd2:True $EXTRACT_D_E_FOR_NON_ALPHA $obsdir $sim $view 2.00 3.00 1.400 2.125 e16 >& $logdir_obs/dobs_R2.log"
        #! echo "dobs_R2 mthd2:False      $EXTRACT_D_E_FOR_NON_ALPHA $obsdir $sim EC_EH_ST 2.00 3.00 1.400 2.125 e16 
        #!
        #! [04-04-16]
        #! + Use views = 'EC_EH' and "EC_ST"
        #!   + Former to track contributions from holes; the latter to contrast new observables with current knowledge and also
        #!     till holes fill out in simulation, to show the function of EC that may till then, get drowned by EH.
        #! [04-07-16]
        #! + After fixing hole filling problem, use view="EC_EF_ST"
        #! + 'mthd1:mthd2:True' -> 'mthd2:False'
		if not SKIP_DOBS_R2:
			OBS_R2_VIEWL=['EC_EF_ST']#!['EC_EF','EC_ST']
			for view in OBS_R2_VIEWL:
				if not os.path.exists(os.path.join(obsdir,sim_total,"Obs_R2_%s"%view)) or PROC_NEW_JOBTAG:
					for EXTRACT_D_E_FOR_NON_ALPHA in ["True","False"]:
						#cmd=["dobs_R2","mthd1:mthd2:True","%s"%EXTRACT_D_E_FOR_NON_ALPHA,obsdir,sim_total,view,Q2MIN,Q2MAX,WMIN,WMAX,EXPT]
						cmd=["dobs_R2","mthd2:False","%s"%EXTRACT_D_E_FOR_NON_ALPHA,obsdir,sim_total,view,Q2MIN,Q2MAX,WMIN,WMAX,EXPT]
						logfile=open("%s/dobs_R2_%s_%s_extrct-D-E-non-alpha-%s.log"%(logdir,sim_total,view,EXTRACT_D_E_FOR_NON_ALPHA),'w')
						print">>>%s >& %s\n"%(cmd,logfile.name)
						mainlog.write(">>>%s >& %s\n"%(cmd,logfile.name))
						tool=subprocess.Popen(cmd,stdout=logfile,stderr=subprocess.STDOUT)
						tool.wait()
				else:
					print os.path.join(obsdir,sim_total,"Obs_R2_%s"%view),"exists. Not making new one"
					mainlog.write("%s exists. Not making new one\n"%os.path.join(obsdir,sim_total,"Obs_R2_%s"%view))
		else:
			print "SKIP_DOBS_R2=True. Therefore not doing dobs_R2"
			mainlog.write("SKIP_DOBS_R2=True. Therefore not doing dobs_R2")


		#-- simstats
		#echo ">ds_lite $obsdir sim1 >& $logdir_obs/ds_lite.log"
		if not os.path.isfile(os.path.join(obsdir,sim_total,"simstats/simstats.root")) or PROC_NEW_JOBTAG:
			cmd=['ds_lite',obsdir,sim_total]
			logfile=open('%s/ds_lite_%s.log'%(logdir,sim_total),'w')
			print">>>%s >& %s/ds_lite_%s.log\n"%(cmd,logdir,sim_total)
			mainlog.write(">>>%s >& %s/ds_lite_%s.log\n"%(cmd,logdir,sim_total))
        		tool=subprocess.Popen(cmd,stdout=logfile,stderr=subprocess.STDOUT)
			tool.wait()
		else:
			print os.path.join(obsdir,sim_total,"simstats/simstats.root"),"exists. Not making new one"
                        mainlog.write("%s exists. Not making new one\n"%os.path.join(obsdir,sim_total,"simstats/simstats.root"))

		#echo "ds_q2wbin_lite $obsdir sim1 >& $logdir_obs/ds_q2wbin_lite.log"
		if not os.path.isfile(os.path.join(obsdir,sim_total,"simstats/simstats_q2wbin.root")) or PROC_NEW_JOBTAG:
			cmd=['ds_q2wbin_lite',obsdir,sim_total]
			logfile=open('%s/ds_q2wbin_lite_%s.log'%(logdir,sim_total),'w')
			print">>>%s >& %s/ds_q2wbin_lite_%s.log\n"%(cmd,logdir,sim_total)
			mainlog.write(">>>%s >& %s/ds_q2wbin_lite_%s.log\n"%(cmd,logdir,sim_total))
        		tool=subprocess.Popen(cmd,stdout=logfile,stderr=subprocess.STDOUT)
			tool.wait()
		else:
			print os.path.join(obsdir,sim_total,"simstats/simstats.root"),"exists. Not making new one"
                        mainlog.write("%s exists. Not making new one\n"%os.path.join(obsdir,sim_total,"simstats/simstats_q2wbin.root"))

	mainlog.close()
#! *** End : Now start processing CUTSNCORS *** 

print "h10_2_per_non_vst_SE_results end time",datetime.datetime.now().strftime('%H:%M:%S')
