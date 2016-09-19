#!/usr/bin/python
import sys,os
from collections import OrderedDict
import datetime

from lib_SS import *

from disp_obs import DispObs

'''
dobs_cmb_non_vst_SE_result.py <date1> <date2> <q2> dbg=[False]
+ <date1>= date tag from h10_2_per_non_vst_SE_results.py
+ <date2>= tag from dobs_per_non_vst_SE_result.py  
'''
#! Take in user arguments 
if len(sys.argv)<2:
	sys.exit("Please enter date1 from h10_2_per_non_vst_SE_results")
DATE1=sys.argv[1]

if len(sys.argv)<3:
        sys.exit("Please enter date2 from dobs_per_non_vst_SE_result.py")
DATE2=sys.argv[2]

if len(sys.argv)<4:
        sys.exit("Please enter q2 = lowQ2/highQ2")
Q2=sys.argv[3]

if len(sys.argv)>4: #! dbg entered by user
	if sys.argv[4]=="False":
		DBG=False
	elif sys.argv[4]=="True":
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
for i,se in enumerate(SES):
	if DBG==True: inputdir_se='SS/dbg/%s_%s_%s_%s_%s'%(Q2,se,DATE1,CUMSIMS[Q2],DATE2)
	else:         inputdir_se='SS/%s_%s_%s_%s_%s'%(Q2,se,DATE1,CUMSIMS[Q2],DATE2)
	OBSD[i+1]=[inputdir_se,se]
#! SS_TAG
SS_TAG='%s_cmb_non_vst_SE'%(Q2)
#! SS_TYPE
SS_TYPE="cmb_non_vst_SE"

#! print all vars
print "DATE1=",DATE1
print "DATE2=",DATE2
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
do.disp_SE_results_1D() 
for EXTRACT_D_E_FOR_NON_ALPHA in [True, False]:
	do.disp_SE_results_R2(EXTRACT_D_E_FOR_NON_ALPHA)

print "***  dobs_cmb_non_vst_SE_result.py Done ***"
