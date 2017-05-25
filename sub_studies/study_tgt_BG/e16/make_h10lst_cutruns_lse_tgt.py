#!/usr/bin/python
import sys,fileinput

#! Obtain runs that need to be cut away by lse and tgt cuts
#! [05-25-17] Taken from $STUDY_LUM_E16_DATADIR/results_052417/lum_results.txt
#! lse,tgt are nrm_N2pi cuts of 26 and 27, respectively
cutruns_lse=[30543, 30563, 30579, 30580, 30585, 30589, 30605, 30634, 30682, 30701, 30709, 30730, 30731, 30738, 30756, 30768, 30771, 30785, 30793, 30834, 30842, 30944, 30965, 30971, 30993, 31036, 31057, 31058, 31066, 31079, 31080, 31085, 31086, 31087, 31145, 31146, 31147, 31162, 31185, 31199, 31200, 31201, 31207, 31231, 31249, 31276, 31363, 31433, 31457, 31466]
num_cutruns_lse=len(cutruns_lse)
cutruns_tgt=[30541, 30543, 30563, 30579, 30580, 30585, 30589, 30593, 30605, 30634, 30682, 30686, 30701, 30709, 30730, 30731, 30738, 30740, 30756, 30768, 30771, 30781, 30785, 30793, 30834, 30842, 30944, 30965, 30971, 30993, 31036, 31037, 31049, 31057, 31058, 31060, 31061, 31063, 31066, 31079, 31080, 31085, 31086, 31087, 31145, 31146, 31147, 31162, 31180, 31185, 31199, 31200, 31201, 31207, 31231, 31249, 31268, 31276, 31363, 31433, 31457, 31466, 31467]
num_cutruns_tgt=len(cutruns_tgt)

#! First make copies of the original file
forg = open('h10_ptgt_h10_2_h10-skim-SS_081716.lst', 'r')
flse = open('h10_ptgt_cutruns_lse_h10_2_h10-skim-SS_081716.lst', 'w')
ftgt = open('h10_ptgt_cutruns_tgt_h10_2_h10-skim-SS_081716.lst', 'w')
for line in forg:
	flse.write(line)
	ftgt.write(line)
forg.close()	
flse.close()
ftgt.close()

#! Now in the copied versions of the original, comment out lines (with #) that contain cutruns
for item in zip([flse,ftgt],[cutruns_lse,cutruns_tgt]):	
	f,cutruns=item[0],item[1]
	#! As per unutbu's solution in  https://stackoverflow.com/questions/7633485/insert-string-at-the-beginning-of-each-line
	for line in fileinput.input([f.name], inplace=True):
		if any(["%d"%run in line for run in cutruns]):
    			sys.stdout.write('#{l}'.format(l=line))
		else:
			sys.stdout.write(line)

#! Check integrity of the job
for item in zip([open(flse.name,'r'),open(ftgt.name,'r')],[num_cutruns_lse,num_cutruns_tgt]):
	f,num_cutruns=item[0],item[1]
	print "*** Checking integrity for %s *** "%f.name
	ncmntd=0
	for line in f:
		if "#" in line: ncmntd+=1
	print "Total runs commented out in %s = %d"%(f.name,ncmntd)
	if ncmntd!=num_cutruns: print("ERROR! Total runs commented out != %d"%num_cutruns)
	else:                   print("Total runs commented out = num_cutruns(=%d)"%num_cutruns)
	print "\n"
	

