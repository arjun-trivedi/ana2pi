#!/usr/bin/python
from __future__ import division

import os,sys,datetime,subprocess,time

from collections import OrderedDict
import itertools


#! Setup constants

#! *** Define WMIN,WMAX ***
WMIN='1.400'
WMAX='2.125'

#! EXPT
EXPT='e16'

#! Setup Q2 ranges
NQ2RANGES=2
LQ2,HQ2=range(NQ2RANGES)
Q2RANGE=[0 for iq2r in range(NQ2RANGES)]
Q2RANGE[LQ2]=['2.00','3.00']
Q2RANGE[HQ2]=['3.00','5.00']

#! SIMTOTAL
SIMTOT=[0 for iq2r in range(NQ2RANGES)]
SIMTOT[LQ2]='sim4_sim5_sim6_sim7_sim8_sim13'
SIMTOT[HQ2]='sim9_sim10_sim11_sim12'

#! Setup OBSDIRS as per Acceptance Cuts
NACC_CUT=2
ACC_CUT_IDVAL=['005','01'] #! 0.5%,1%
ACC_CUT=['0.005','0.01'] #! 0.5%,1%
OBSDIR=[[0 for i in range(NACC_CUT)] for iq2r in range(NQ2RANGES)]
for iacc in range(NACC_CUT):
	acc_cut=ACC_CUT_IDVAL[iacc]
	for iq2r in range(NQ2RANGES):
		if   iq2r==LQ2: name,date="low",'102617'
		elif iq2r==HQ2: name,date="high",'102717'
		OBSDIR[iq2r][iacc]='%s/%sQ2_SSBands_off_off_%s_%s_acc_cut_%s/cutsncors1'%(os.environ['OBSDIR_E16'],name,SIMTOT[iq2r],date,acc_cut)
#print OBSDIR
#sys.exit()

#-- Now start to make Obervables (Obs_1D only)
for iq2r in range(NQ2RANGES):
	q2min,q2max=Q2RANGE[iq2r][0],Q2RANGE[iq2r][1]
	for obsdir,acc_cut in zip(OBSDIR[iq2r],ACC_CUT):
		print "*** Processing obs_1D for",obsdir,SIMTOT[iq2r],q2min,q2max,acc_cut,"***\n"
		#! ph8
		cmd=["ph8",obsdir,SIMTOT[iq2r],q2min,q2max,WMIN,WMAX,acc_cut]
		print "ph8 cmd=",cmd
		logfile=open("/tmp/ph8_%s_%s.log"%(SIMTOT[iq2r],acc_cut),'w')
		tool=subprocess.Popen(cmd,stdout=logfile,stderr=subprocess.STDOUT)
		tool.wait()
		#! dobs1D
		cmd=["dobs_1D",obsdir,SIMTOT[iq2r],"norm",q2min,q2max,WMIN,WMAX,EXPT]
		print "dobs_1D cmd=",cmd
		logfile=open("/tmp/dobs1D_%s_%s.log"%(SIMTOT[iq2r],acc_cut),'w')
		tool=subprocess.Popen(cmd,stdout=logfile,stderr=subprocess.STDOUT)
		tool.wait()


       