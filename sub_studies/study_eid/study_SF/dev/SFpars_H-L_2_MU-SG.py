#!/usr/bin/python
from __future__ import division

import os,sys,time
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict

import numpy as np

"""
Usage: SFpars_H-L_2_MU-SG.py  fSFpars_H-L fSFpars_MU-SG 

+ Takes file with SF pars in format of cut_h/l(obtained for e16:SF) and converts them to mean/sigma(obtained for e1f:SF):
	s1:cut_h(pol3:p0:p1:p2:p3)
	s1:cut_l(pol3:p0:p1:p2:p3)
	.
	.
	.
	s6:cut_h(pol3:p0:p1:p2:p3)
	s6:cut_l(pol3:p0:p1:p2:p3)
	->
	s1:mean(pol3:p0:p1:p2:p3)
	s1:sigma(pol3:p0:p1:p2:p3)
	.
	.
	.
	s6:mean(pol3:p0:p1:p2:p3)
        s6:sigma(pol3:p0:p1:p2:p3)
+ I need to do be able to do this for various technical reasons, for example, when using 'hvy' machinery for e16, it can read 
  SF-pars only in mean/sigma format and therefore I have to convert SF-pars that are directly obtained in cut_h/l.


"""
#! Get args from user
fin=open(sys.argv[1],'r')
fout=open(sys.argv[2],'w')

NSCTR=6	
NPAR=4 #! using pol3 (1 more for constant!)

H=[[0 for j in range(NPAR)] for i in range(NSCTR)]
L=[[0 for j in range(NPAR)] for i in range(NSCTR)]
MU=[[0 for j in range(NPAR)] for i in range(NSCTR)]
SG=[[0 for j in range(NPAR)] for i in range(NSCTR)]
data=fin.readlines()
for line in data:
	words=line.split()
	for isctr in range(NSCTR):
		if "s%d:cut_h"%(isctr+1) in words[0]:
			for ipar in range(NPAR):
				H[isctr][ipar]=float(words[ipar+1])
		elif "s%d:cut_l"%(isctr+1) in words[0]:
			for ipar in range(NPAR):
				L[isctr][ipar]=float(words[ipar+1])
			
		#! Now set calculate MU,SG as:
		#! MU=(H+L)/2 
		#! SG=(H-L)/6
		#! where the following equations are used:
		#! +H=M+3*SG
		#! +L=M-3*SG
		for ipar in range(NPAR):
			h=H[isctr][ipar]
			l=L[isctr][ipar]
			MU[isctr][ipar]=(h+l)/2
			SG[isctr][ipar]=(h-l)/6
#! Write
for isctr in range(NSCTR):
	fout.write("s%d:mean(pol3:p0:p1:p2:p3) %f %f %f %f\n"%(isctr+1,MU[isctr][0],MU[isctr][1],MU[isctr][2],MU[isctr][3]))
	fout.write("s%d:sgma(pol3:p0:p1:p2:p3) %f %f %f %f\n"%(isctr+1,SG[isctr][0],SG[isctr][1],SG[isctr][2],SG[isctr][3]))
