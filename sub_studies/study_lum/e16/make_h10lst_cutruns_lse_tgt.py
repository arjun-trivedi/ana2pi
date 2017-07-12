#!/usr/bin/python
import sys,fileinput
import datetime

'''
+ Started as a copy of $SUBSTUDIES/study_tgt_BG/e16/make_h10lst_cutruns_lse_tgt.py

+ Takes the full h10lst as the input argument and produces h10lsts as per cuts defined in the code
+ These new h10lsts are put in the same directory as the original

+ >$SUBSTUDIES/study_lum/e16/make_h10lst_cutruns_lse_tgt.py h10lst_full

'''
#! Obtain runs that need to be cut away by lse and tgt cuts
#! [05-25-17] Taken from $STUDY_LUM_E16_DATADIR/results_052417/lum_results.txt
#! lse,tgt are nrm_N2pi cuts of 26 and 27, respectively
#cutruns_lse=[30543, 30563, 30579, 30580, 30585, 30589, 30605, 30634, 30682, 30701, 30709, 30730, 30731, 30738, 30756, 30768, 30771, 30785, 30793, 30834, 30842, 30944, 30965, 30971, 30993, 31036, 31057, 31058, 31066, 31079, 31080, 31085, 31086, 31087, 31145, 31146, 31147, 31162, 31185, 31199, 31200, 31201, 31207, 31231, 31249, 31276, 31363, 31433, 31457, 31466]
#cutruns_tgt=[30541, 30543, 30563, 30579, 30580, 30585, 30589, 30593, 30605, 30634, 30682, 30686, 30701, 30709, 30730, 30731, 30738, 30740, 30756, 30768, 30771, 30781, 30785, 30793, 30834, 30842, 30944, 30965, 30971, 30993, 31036, 31037, 31049, 31057, 31058, 31060, 31061, 31063, 31066, 31079, 31080, 31085, 31086, 31087, 31145, 31146, 31147, 31162, 31180, 31185, 31199, 31200, 31201, 31207, 31231, 31249, 31268, 31276, 31363, 31433, 31457, 31466, 31467]

#! [07-10-17] Taken from $STUDY_LUM_E16_DATADIR/results_052417/lum_results.txt
#! + Only 1 cut with lower and upper bound=(26,34)
#! + To continue using same code, will call this 1 cut 'lse' and 'tgt' will have no cuts applied
cutruns_lse=[30543, 30563, 30579, 30580, 30585, 30589, 30605, 30634, 30682, 30701, 30709, 30730, 30731, 30738, 30748, 30756, 30768, 30771, 30785, 30793, 30834, 30842, 30944, 30965, 30971, 30993, 30995, 31036, 31057, 31058, 31066, 31079, 31080, 31085, 31086, 31087, 31145, 31146, 31147, 31162, 31185, 31199, 31200, 31201, 31207, 31231, 31241, 31249, 31276, 31363, 31417, 31433, 31457, 31466]
cutruns_tgt=[]

num_cutruns_lse=len(cutruns_lse)
num_cutruns_tgt=len(cutruns_tgt)

#! Get user argument
h10lst_full=sys.argv[1]

#! Set outdir= path from the original file
outdir=h10lst_full.replace(h10lst_full.split('/')[-1],'')
#print "outdir=",outdir
#sys.exit()

#! Set appropriate names for files for lse and tgt
DATE=datetime.datetime.now().strftime('%m%d%y')
fname_full=h10lst_full.split('/')[-1]
fname_lse=fname_full.replace('.lst','_lse_%s.lst'%DATE)
fname_tgt=fname_full.replace('.lst','_tgt_%s.lst'%DATE)
#print "fname_lse=",fname_lse
#print "fname_tgt=",fname_tgt
#sys.exit()

#! First open ffull and create empy flse, ftgt
ffull = open(h10lst_full, 'r')
flse = open('%s/%s'%(outdir,fname_lse), 'w')
ftgt = open('%s/%s'%(outdir,fname_tgt), 'w')
#! Add entries to flse based on cutruns_lse
for line in ffull:
	if any(["%d"%run in line for run in cutruns_lse]): continue
	flse.write(line)
#! Add entries to ftgt based on cutruns_tgt
#! First go to begining of ffull 
ffull.seek(0)
for line in ffull:
	print "line=",line
	if any(["%d"%run in line for run in cutruns_tgt]): continue
	ftgt.write(line)
ffull.close()	
flse.close()
ftgt.close()
#sys.exit()

#! Check integrity of the job
for item in zip([open(flse.name,'r'),open(ftgt.name,'r')],[num_cutruns_lse,num_cutruns_tgt]):
	f,num_cutruns=item[0],item[1]
	print "*** Checking integrity for %s *** "%f.name
	nlines=sum(1 for line in f)
	print "Total runs in %s = %d"%(f.name,nlines)
	ndel=606-nlines
	if ndel!=num_cutruns:   
		print("ERROR! Total runs deleted != %d"%num_cutruns)
	else:                   
		print("Passed integrity check because total runs deleted = number of cutruns(=%d)"%num_cutruns)
		print("%d+%d=%d(=total files)"%(nlines,ndel,nlines+ndel))
	print "\n"
	

