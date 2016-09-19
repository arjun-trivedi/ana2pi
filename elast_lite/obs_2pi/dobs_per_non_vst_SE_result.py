#!/usr/bin/python
import sys,os
from collections import OrderedDict
import datetime

from lib_SS import *

from disp_obs import DispObs

'''
dobs_per_non_vst_SE_result.py <date> <q2> <se> dbg=[False]
+ <date>= date tag from h10_2_per_non_vst_SE_results.py
'''
#! Take in user arguments 
if len(sys.argv)<2:
	sys.exit("Please enter date for h10_2_per_non_vst_SE_results")
DATE=sys.argv[1]

if len(sys.argv)<3:
        sys.exit("Please enter q2 = lowQ2/highQ2")
Q2=sys.argv[2]

if len(sys.argv)<4:
        sys.exit("Please enter SE")
SE=sys.argv[3]

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
ncutsncors=len(SECODES[SE])
for i in range(ncutsncors):
	cutncor_inputdir='%s_%s_%s/cutsncors%s/%s'%(Q2,SE,DATE,i+1,CUMSIMS[Q2])
	cutncor_code=SECODES[SE][i]
	OBSD[i+1]=[cutncor_inputdir,cutncor_code]
#! SS_TAG
SS_TAG='%s_%s_%s_%s'%(Q2,SE,DATE,CUMSIMS[Q2])
#! SS_TYPE
SS_TYPE="per_non_vst_SE"

#! print all vars
print "DATE=",DATE
print "Q2=",Q2
print "SE=",SE
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

print "*** dobs_per_non_vst_SE_result.py Done ***"
