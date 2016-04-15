#!/usr/bin/python

import subprocess, shlex
import os

obsdir=os.path.join(os.environ['OBSDIR_E16'],'lowQ2_new_EH_calc_040616')

for cutsncors in ['cutsncors1','cutsncors2','cutsncors3','cutsncors4','cutsncors5']:
	#! dobs_1D
	cmd="dobs_1D %s/%s sim4_sim5_sim6_sim7_sim8_sim13 norm 2.00 3.00 1.400 2.125 e16 True"%(obsdir,cutsncors)
	#print shlex.split(cmd)
	p=subprocess.Popen(shlex.split(cmd))
	p.wait()

	#! dobs_R2
	cmd="dobs_R2 mthd2:False True %s/%s  sim4_sim5_sim6_sim7_sim8_sim13 EC_EF_ST 2.00 4.00 1.400 2.125 e16 True"%(obsdir,cutsncors)
	p=subprocess.Popen(shlex.split(cmd))
	p.wait()

