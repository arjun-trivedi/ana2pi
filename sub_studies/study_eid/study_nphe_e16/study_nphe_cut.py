#!/usr/bin/python
from __future__ import division

import os,sys
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict
import itertools

import numpy as np

import math

import atlib as attools

'''
+ 

+ Usage: 
'''

USAGE=''

DBG=False
if len(sys.argv)>1: #i.e. debug entered by user
	if    sys.argv[1]=="True": DBG=True
	elif  sys.argv[1]=="False": DBG=False
	else: sys.exit('DBG=%s is not valid. usage: %s'%(DBG,USAGE))

print "DBG=",DBG
#sys.exit()

#! *** Prepare constants 
NSCTRS=6

NSGMNTS=18

NPMTS=2
L,R=range(NPMTS)

#! *** Prepare output datadir***
OUTDIR=os.environ['STUDY_EID_NPHE_E16_DATADIR']
if DBG==True: OUTDIR+="/dbg"
print "OUTDIR=",OUTDIR
if not os.path.exists(OUTDIR):
	os.makedirs(OUTDIR)
#sys.exit()
#! ***

#! *** Prepare NENTRIES ***
NENTRIES=1000000000
if DBG==True: NENTRIES=1000

#! *** Now fill FIN/T[dtyp][cut] ***
INPUT=os.path.join(os.environ['D2PIDIR_EXP_E16'],"data_nphe_071916")#! debug /bk
#print "INPUT=",INPUT
FIN=ROOT.TFile("%s/dnphe_nnphe.root"%(INPUT))
print "FIN=",FIN.GetName()
T=FIN.Get("d2piR/tR")
print "T=",T.GetName()
#sys.exit()
#! ***

#! *** Setup objects that will be created
hnphe=[[[ROOT.TH1F("hnphe_%d_%d_%d"%(i+1,j+1,k+1),"hnphe_sct%d_sgm%d_pmt%d"%(i+1,j+1,k+1),100,0,300) for k in range(NPMTS)]for j in range(NSGMNTS)]for i in range(NSCTRS)]
print hnphe
#sys.exit()
#! *** Now make plots ***

#! Process of looping over TTree in Python taken from Evan's mail of [02-7-15] 'pyroot and bytes in trees'
for i, ev in enumerate(T, 1):
	if (i%10000==0): print "processing ientry# %d"%(i)
	if DBG and i>1000: break
	
	#print "***entry #%d ***"%(i+1)

	#! Cuts
	#! 1. top cut: Use only top2 events even though top2' is actually top2+top1
	#! + This is because in top1, as cacluated by 'hvy', mmppip=0 
	#print "ev top=",ev.top
	if ev.top!=2: continue # and ev.top!=1: continue

	#! 2. Q2, W cut: to keep Q2 & W withing analysis range :[2.0,5.0], [1.400,2.125)
	if ev.Q2<2.0 or ev.Q2>=5.0: continue
	if ev.W<1.400 or ev.W>=2.125: continue

	isct=ev.sector-1;
	isgm=ev.cc_segm-1;
	if   ev.pmt==-1: ipmt=0
	elif ev.pmt==1: ipmt=1
	else: continue #! if pmt==0, then eff=1 

	#! Fill histograms
	hnphe[isct][isgm][ipmt].Fill(ev.nphe)
	

#! ***
#sys.exit()

#! Fit histograms
def get_ffnphe():
	'''
	+ fit function as obtained from YT (verified to same as KP)
	+ fit function init-pars as per YT
	'''
	fnphe=ROOT.TF1("fnphe", "[0]*TMath::Power(([1]),(x/[2]))*(TMath::Exp(-([1])))/TMath::Gamma((x/[2])+1)",50,300)
	fnphe.SetParName(0,"A")
	fnphe.SetParName(1,"L")
	fnphe.SetParName(1,"p")
	fnphe.SetParameters(58294,2.28921,38.0199)
	return fnphe
#! Now fit
d=itertools.product(*[range(NSCTRS),range(NSGMNTS),range(NPMTS)])
for r in d:
	isct,isgm,ipmt=r[0],r[1],r[2]
	f=get_ffnphe()
	hnphe[isct][isgm][ipmt].Fit(f)

#! Save histograms with fits in following format:
#! sector/canvas_L/R(6*3)
FOUT=ROOT.TFile("%s/study_nphe_cut.root"%OUTDIR,"RECREATE")
for isct in range(NSCTRS):
	FOUT.mkdir("sector%d"%(isct+1)).cd()
	for ipmt in range(NPMTS):
		cname=("sct%d_pmt%d"%(isct+1,ipmt+1))
		c=ROOT.TCanvas(cname,cname,3000,1000)
		c.Divide(6,3)
		for isgm in range(NSGMNTS):
			c.cd(isgm+1)
			hnphe[isct][isgm][ipmt].Draw()
		c.Write("%s"%cname,ROOT.gROOT.ProcessLine("TObject::kOverwrite"))#! for .root; latest namecycle
	
FOUT.Close()	

#! Obtain eff[sct][sgm][pmt]
eff=[[[0 for k in range(NPMTS)]for j in range(NSGMNTS)]for i in range(NSCTRS)]
d=itertools.product(*[range(NSCTRS),range(NSGMNTS),range(NPMTS)])
for r in d:
        isct,isgm,ipmt=r[0],r[1],r[2]
	f=hnphe[isct][isgm][ipmt].GetFunction("fnphe")
	print "Obtaining eff for isct,isgm,ipmt=",isct,isgm,ipmt
	if f!=None:
		print "Fit function found"
		print "eff=",f.Integral(20,400)/f.Integral(0,400)
		eff[isct][isgm][ipmt]=f.Integral(20,400)/f.Integral(0,400)#! extend 300 -> 400 in the sense of approaching infinity 
	else:
		print "Fit function not found" 

#! Write eff to file in the format:
#! sct sgm pmt eff
EFF_DATA=open("%s/eid_e16_exp_nphe_eff.txt"%OUTDIR,"w")
d=itertools.product(*[range(NSCTRS),range(NSGMNTS),range(NPMTS)])
for r in d:
	isct,isgm,ipmt=r[0],r[1],r[2]
	EFF_DATA.write("%d %d %d %f\n"%(isct+1,isgm+1,ipmt+1,eff[isct][isgm][ipmt]))
EFF_DATA.close()
	
#! If wanting to keep TCanvas open till program exits				
#if not ROOT.gROOT.IsBatch():
#	plt.show()
#	# wait for you to close the ROOT canvas before exiting
#	wait(True)

#if __name__ == "__main__":	
#	if len(sys.argv)==2:
#		plot_fid(sys.argv[1])
#	elif len(sys.argv)==3:
#		plot_fid(sys.argv[1],int(sys.argv[2])) 
