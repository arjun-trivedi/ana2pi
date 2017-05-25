#!/usr/bin/python
from __future__ import division

import os,sys
import datetime
import glob
import subprocess

DATE=datetime.datetime.now().strftime('%m%d%y')

#! input and output dirs
LUMBYRUNDIR=os.path.join(os.environ['D2PIDIR_EXP_E16'],'h10-etgt_2_dlum-byRun_%s'%DATE)
if not os.path.exists(LUMBYRUNDIR):
	os.makedirs(LUMBYRUNDIR)

#! log files
LOGDIR=os.path.join(LUMBYRUNDIR,'log')
if not os.path.exists(LOGDIR):
	os.makedirs(LOGDIR)
MAINLOG=open("%s/main.log"%(LOGDIR),'w')

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
		cmd=["proc_h10","-i","%s"%h10lst,"-t","e16:exp:2pi","-p","lum", "-o","%s/%d.root"%(LUMBYRUNDIR,runnum)]
		MAINLOG.write(">>>%s >& %s\n"%(cmd,logfile.name))
		t=subprocess.Popen(cmd,stdout=logfile,stderr=subprocess.STDOUT)
		t.wait()
