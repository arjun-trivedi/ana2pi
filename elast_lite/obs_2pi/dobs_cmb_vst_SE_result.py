#!/usr/bin/python
import sys,os
from collections import OrderedDict
import datetime

from lib_SS import *

from disp_obs import DispObs

'''
dobs_cmb_vst_SE_result.py <date> <q2> dbg=[False]
+ <date>= date from dobs_cmb_non_vst_SE_result.py  
'''
#! Take in user arguments 
if len(sys.argv)<2:
	sys.exit("Please enter date for h10_2_per_non_vst_SE_results")
DATE=sys.argv[1]

if len(sys.argv)<3:
        sys.exit("Please enter q2 = lowQ2/highQ2")
Q2=sys.argv[2]

if len(sys.argv)>3: #! dbg entered by user
	if sys.argv[3]=="False":
		DBG=False
	elif sys.argv[3]=="True":
		DBG=True
	else:
		print "Please enter dbg=True or False only."
		sys.exit()
else:
	DBG=False

#! setup vars as per user input
Q2MIN,Q2MAX=Q2LIMITS[Q2][0],Q2LIMITS[Q2][1]

#! OBSD,SS_TAG,SS_TYPE
#! OBSD
OBSD=OrderedDict()
if DBG==True: inputdir='SS/dbg/%s_cmb_non_vst_SE_%s'%(Q2,DATE)
else:         inputdir='SS/%s_cmb_non_vst_SE_%s'%(Q2,DATE)
OBSD[1]=[inputdir,'cmb_non_vst_SE']
#! SS_TAG
SS_TAG='%s_cmb_vst_SE'%(Q2)
#! SS_TYPE
SS_TYPE="cmb_vst_SE"

#! print all vars
print "DATE=",DATE
print "Q2=",Q2
print "DBG=",DBG
print "Q2MIN=",Q2MIN
print "Q2MAX=",Q2MAX
print "WMIN=",WMIN
print "WMAX=",WMAX
print "EXPT=",EXPT
print "OBSD pretty print:"
for k in OBSD.keys():
	print k, OBSD[k]
print "SS_TAG=",SS_TAG
print "SS_TYPE=",SS_TYPE
print "DBG=",DBG
print "******"
#sys.exit()
#! Check integrity of OBSD because these directories have to be present since they are the basic input to the process
for k in OBSD:
        drctry="%s/%s"%(os.environ['OBSDIR_E16'],OBSD[k][0])
	if not os.path.exists(drctry):
        	sys.exit("%s does not exist. Please check input arguments"%drctry)
#sys.exit()

#! Run
do=DispObs('NA','NA','NA',Q2MIN,Q2MAX,WMIN,WMAX,EXPT,DBG,obsd=OBSD,SS_tag=SS_TAG,SS_type=SS_TYPE)
do.disp_cmb_vst_SE_results_1D() 
for EXTRACT_D_E_FOR_NON_ALPHA in [True, False]:
	do.disp_cmb_vst_SE_results_R2(EXTRACT_D_E_FOR_NON_ALPHA)

print "***  dobs_cmb_vst_SE_result.py Done ***"
