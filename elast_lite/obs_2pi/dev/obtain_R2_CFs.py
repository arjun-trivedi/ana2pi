#!/usr/bin/python
from __future__ import division

import os,sys
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict
import math
import numpy as np

import scipy.integrate as integrate

import tools_extract_R2 as tR2

'''
[03-27-16]
+ This code documents the Correction Factors (CFs) applied to the theoretical R2 Extraction Factors (R2EFs) for mthd1 and mthd2 (For details see handwritten notes)
'''
OUTPUT="%s/%s"%(os.environ['ELAST_LITE'],'obs_2pi/dev/R2_EF_CFs.txt')
sys.stdout=open(OUTPUT, 'w')

#! print PHI_BIN_LEL declared in tR2
print "PHI_BIN_LEL=",tR2.PHI_BIN_LEL

def print_CFs(R2):
	if R2!='A' and R2!='B' and R2!='C' and R2!='D':
		print "get_CFs() not implemented for R2=",R2

	#! Obtain integral for R2 related sinusoids over these phi bins
	intgrll=tR2.get_intgrll_R2f(R2)
	print "integral for sinusoid related to R2 over phi bins:",intgrll

	print "+++ Obtaining CF for mthd1 +++"
	#! + Based on the assumption that the orthogonal binned integrals are 0
	#! + This assumption may not always be true and contributes to the limitations of mthd1
  
	#! 1. Get non-orthogonal binned integral i.e. the integral that separates R2
        binned_EF=tR2.get_mthd1_binned_EF(R2,intgrll)
        print "binned EF=",binned_EF
	#! 2. Get true EF
	true_EF=tR2.get_mthd1_true_EF(R2)
	print "true EF=",true_EF
	#! 3. Now calculate CF=sum(intgrl_mbyl)/EF
        CF_m1=binned_EF/true_EF
	
	print "CF for EF mthd1=%.5f"%CF_m1
	print "++++++"

	
	print "+++ Obtain CF list for mthd2 +++"
	#! Get list of CFs in each bin, which are identical and therefore a single CF is finally used
	CFl_m2=tR2.get_mthd2_CFl(R2,intgrll)
	print "CFl for m2",CFl_m2
	print "++++++"
		
for R2 in ['A','B','C','D','E']:
	print "*** %s ***"%R2
	print_CFs(R2)
	print "******\n"
#print "fctr used so far=",1/math.radians(phi_binw)

#! If wanting to keep TCanvas open till program exits				
if not ROOT.gROOT.IsBatch():
	plt.show()
	# wait for you to close the ROOT canvas before exiting
	wait(True)

#if __name__ == "__main__":	
#	if len(sys.argv)==2:
#		plot_fid(sys.argv[1])
#	elif len(sys.argv)==3:
#		plot_fid(sys.argv[1],int(sys.argv[2])) 
	
