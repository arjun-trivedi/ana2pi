#!/usr/bin/python
from __future__ import division

import os,sys
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict

import numpy as np

from cuts_theta_vs_p_e16 import NCUTLNS,H,L,CUT
print "NCUTLNS,H,L=",NCUTLNS,H,L
#print "CUT=",CUT

NSCTR=6

#! Setup particle information
NPRTCL=3
E,P,PIP=range(NPRTCL)
PRTCL_NAME=["e","p","pip"]

h=[[]for i in range(NPRTCL)]
h[E]=ROOT.TH2F("he","he",      100,0,5,100,0,60)
h[P]=ROOT.TH2F("hp","hp",      100,0,4,100,0,60)
h[PIP]=ROOT.TH2F("hpip","hpip",100,0,3,100,0,120)
#hpip=ROOT.TH2F("hpip","hpip",100,0,120,100,0,5) #! since definition of pip are cuts are inverse: p(theta) instead of theta(p)

#! Draw electron cuts
c=[[]for i in range(NPRTCL)]
for iprtcl in range(NPRTCL):
	prtcl=PRTCL_NAME[iprtcl]
	#if prtcl=="e" or prtcl=="p":continue
	cname=PRTCL_NAME[iprtcl]
	c[iprtcl]=ROOT.TCanvas(cname,cname)
	c[iprtcl].Divide(3,2)
	for isctr in range(NSCTR):
		if CUT[H][iprtcl][isctr]==None: continue
		c[iprtcl].cd(isctr+1)
		h[iprtcl].Draw()
		#if prtcl!="pip":h.Draw()
		#else:		hpip.Draw()
		for cuth,cutl in zip(CUT[H][iprtcl][isctr],CUT[L][iprtcl][isctr]):
			cuth.Draw("same")
			cutl.Draw("same")
#ce=ROOT.TCanvas("ce","ce")
#ce.Divide(3,2)
#for isctr in range(NSCTR):
#	if CUT[H][E][isctr]==None: continue
#	ce.cd(isctr+1)
#	h.Draw()
#	for cuth,cutl in zip(CUT[H][E][isctr],CUT[L][E][isctr]):
#		cuth.Draw("same")
#		cutl.Draw("same")



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
	
