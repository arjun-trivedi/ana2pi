#!/usr/bin/python
from __future__ import division

import os,sys
import datetime
import glob
import subprocess

DATE=datetime.datetime.now().strftime('%m%d%y')

#! input and output dirs
D2PIRBYRUNDIR=os.path.join(os.environ['D2PIDIR_EXP_E16'],'h10-etgt_2_d2piR-byRun_%s'%DATE)
if not os.path.exists(D2PIRBYRUNDIR):
	os.makedirs(D2PIRBYRUNDIR)

#! log files
LOGDIR=os.path.join(D2PIRBYRUNDIR,'log')
if not os.path.exists(LOGDIR):
	os.makedirs(LOGDIR)
MAINLOG=open("%s/main.log"%(LOGDIR),'w')

#! Following taken directly from h10_2_per_non_vst_SE_results.py
BASIC_PROCORDER='eid:efid:pid:pfid:pcorr:eff:scpd:evtsel_2pi:'
BASIC_ADTNL_OPTS=':1:3:17:2:18:20:' #! ECin-on:zvtx-on:eff_scpd_atmod:ECfid-ON:ECfid_atmod_ON:CC_cut_lse:

#! Make temporary directory to keep all runnum.lst files
h10lstdir='/home/trivedia/ongoing/h10lsts_local/e16/exp/etgt'
#for h10lst in os.path(h10lstdir):
for dirname, dirnames, filenames in os.walk(h10lstdir):
	for filename in filenames:
		#! Now process run
		h10lst="%s/%s"%(dirname,filename)
		#print h10lst
		runnum=int(filename.split('.')[0])
		#print runnum 
		#sys.exit()
		logfile=open('%s/%d.log'%(LOGDIR,runnum),'w')
		cmd=["proc_h10_lite","-i","%s"%h10lst,"-t","e16:exp:2pi:recon","-c","%s"%BASIC_PROCORDER, "-o","%s/%d.root"%(D2PIRBYRUNDIR,runnum), "%s"%BASIC_ADTNL_OPTS]
		MAINLOG.write(">>>%s >& %s\n"%(cmd,logfile.name))
		t=subprocess.Popen(cmd,stdout=logfile,stderr=subprocess.STDOUT)
		t.wait()
