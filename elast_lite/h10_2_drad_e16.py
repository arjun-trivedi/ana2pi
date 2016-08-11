#!/usr/bin/python
from __future__ import division

import os,sys
from array import *
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict
import itertools

import numpy as np

import math

from multiprocessing import Process, Queue

'''
+ 

+ Usage: h10_2_drad.py <simnum>
'''

USAGE='h10_2_drad.py <simnum> make_hq2w=[True] debug=[False]'

SIMNUM=''
if len(sys.argv)<2:
        print("Please enter simnum as simx, where x=1,2,3,4 etc")
        sys.exit()
SIMNUM=sys.argv[1]

MAKE_HQ2W=True
if len(sys.argv)>2: #i.e. debug entered by user
        if    sys.argv[2]=="True": MAKE_HQ2W=True
        elif  sys.argv[2]=="False": MAKE_HQ2W=False
        else: sys.exit('MAKE_HQ2W=%s is not valid. usage: %s'%(sys.argv[2],USAGE))

DBG=False
if len(sys.argv)>3: #i.e. debug entered by user
        if    sys.argv[3]=="True": DBG=True
        elif  sys.argv[3]=="False": DBG=False
        else: sys.exit('DBG=%s is not valid. usage: %s'%(sys.argv[3],USAGE))

print "SIMNUM=",SIMNUM
print "MAKE_HQ2W=",MAKE_HQ2W
print "DBG=",DBG
#sys.exit()

#! *** setup some constants ***
ELECTRON=11
MASS_E=0.000511
MASS_P=0.93827203
E16_E0_P=5.754

lvE0=ROOT.TLorentzVector(0,0,0,0)
lvE0.SetXYZM(0,0,E16_E0_P,MASS_E)

lvP0=ROOT.TLorentzVector(0,0,0,0)
lvP0.SetXYZM(0,0,0,MASS_P)

#! Q2,W for hq2w: taken from '$ELAST_LITE/obs_2pi/ds_lite' which follows the master key for q2w bng: elast_lite/h8_bng.h
#! W contains some extra bins beyond the lower and upper limits of the analysis range to see overflow
WBINW=0.025
NWBINS,WMIN,WMAX=36,1.300,2.200# 
NQ2BINS,Q2MIN,Q2MAX=8,1.25,5.25 # place-holder bng, for it will be modified by variable bng
#! variable Q2 binning is different from 'elast_lite/h8_bng.h' in that there is 1 less bin: Q2=[1.5,2.0] is removed here
NBINS_Q2_VAR_BINW_BNG=5;
Q2_BINS_VAR_BINW_BNG=array("d",[2.00,2.40,3.00,3.50,4.20,5.00])


#! *** Prepare input data: H10LST[rad] ***
NRAD=2
ON,OFF=range(NRAD)
RAD_NAME=["on","off"]

#! *** Prepare H10LSTDIR and H10LST[rad]
H10LSTDIR=os.path.join(os.environ['D2PIDIR_SIM_E16'],'rad-eff-corr-sim',SIMNUM)
H10LST=[0 for i in range(NRAD)]
for irad in range(NRAD):
        if   irad==ON:  H10LST[irad]="%s/h10_RE-on.lst"%(H10LSTDIR)
        elif irad==OFF: H10LST[irad]="%s/h10_RE-off.lst"%(H10LSTDIR)
print "H10LST=",H10LST

#! *** OUTDIR etc
OUTDIR=os.path.join(os.environ['D2PIDIR_SIM_E16'],'rad-eff-corr-sim',SIMNUM)
FOUT_NAME="drad.root"
if DBG: FOUT_NAME="drad_dbg.root"

def Loop(irad,que=None):
	print("*** processing h10 for %s *** "%RAD_NAME[irad])
	h10chain=ROOT.TChain("h10");
	fc=ROOT.TFileCollection("fileList", "", H10LST[irad]);
	h10chain.AddFileInfoList(fc.GetList());
	#print h10chain.GetEntries();

	#! Create hq2w
	hq2w=ROOT.TH2F("hq2w_rad_%s"%(RAD_NAME[irad]),"Radiative Effects %s(Q2,W)"%(RAD_NAME[irad]),NWBINS,WMIN,WMAX,NQ2BINS,Q2MIN,Q2MAX)
        hq2w.GetYaxis().Set(NBINS_Q2_VAR_BINW_BNG,Q2_BINS_VAR_BINW_BNG)
	hq2w.SetXTitle("W [GeV]")
	hq2w.SetYTitle("Q^{2} [GeV]")
	#! Set Q2 bin labels 
	#! (makes it direct to label Q2 projections histograms which will all derive from this)
	nq2bins=hq2w.GetNbinsY()
	for iq2bin in range(nq2bins):
		le=hq2w.GetYaxis().GetBinLowEdge(iq2bin+1)
		ue=hq2w.GetYaxis().GetBinLowEdge(iq2bin+1+1)
		hq2w.GetYaxis().SetBinLabel(iq2bin+1,"[%.2f,%.2f)"%(le,ue))

	for i, ev in enumerate(h10chain, 1):
		if DBG:
			if (i%10000==0): print "processing ientry# %d"%i
                	if  i>10001: break
		else:
			if (i%100000==0): print "processing ientry# %d"%i

		lvE1=ROOT.TLorentzVector(0,0,0,0)
		lvQ=ROOT.TLorentzVector(0,0,0,0)
		lvW=ROOT.TLorentzVector(0,0,0,0)
		for idx in range(ev.mcnentr):
			id=ev.mcid[idx];
			if id != ELECTRON: continue
			mom  =ev.mcp[idx];
			theta=math.radians(ev.mctheta[idx])
			phi  =math.radians(ev.mcphi[idx])
			px=mom*math.cos(phi)*math.sin(theta);
			py=mom*math.sin(phi)*math.sin(theta);
			pz=mom*math.cos(theta);
			lvE1.SetPxPyPzE(px,py,pz,math.sqrt(mom*mom+MASS_E*MASS_E));
		lvQ=lvE0-(lvE1)
		lvW=lvQ+lvP0
		Q2=-1*lvQ.Mag2()
		W=lvW.Mag()
		#print "Q2,W=",Q2,W
		#print "Filling",hq2w.GetName()
		hq2w.Fill(W,Q2)
	if que!=None:
		#que.put(0)
		que.put(hq2w)
	return hq2w

def make_corr_factor_and_dgnstc_plts():
	fout=ROOT.TFile("%s/%s"%(OUTDIR,FOUT_NAME),"UPDATE")	
        hq2w_on=fout.Get("hq2w_rad_on")
	hq2w_off=fout.Get("hq2w_rad_off")

	#! Create corr factor(1/R)= rad_off/rad_on
	#! First create 2D i.e in Q2,W
	hq2w_on.Sumw2()
	hq2w_off.Sumw2()
	h2cf=hq2w_off.Clone("h2cf")
	h2cf.SetTitle("Radiative Effect Correction factor(Q^{2},W)")
	h2cf.Divide(hq2w_on)
	h2cf.Write(h2cf.GetName(),ROOT.gROOT.ProcessLine("TObject::kOverwrite"))
	#! Now make Q2 bin projections in W
	nq2bins=h2cf.GetNbinsY()
	hcf=[0 for i in range(nq2bins)]
	for iq2bin in range(nq2bins):
		hname="hcf_qbin_%d"%(iq2bin+1)
		htitle=h2cf.GetYaxis().GetBinLabel(iq2bin+1)
		hcf[iq2bin]=h2cf.ProjectionX(hname,iq2bin+1,iq2bin+1,"e")
		hcf[iq2bin].SetTitle("Radiative Effects Correction Factor(1/R) for Q2=%s"%htitle)
		hcf[iq2bin].SetXTitle("W [GeV]")
		hcf[iq2bin].SetYTitle("1/R")
		hcf[iq2bin].Write(hcf[iq2bin].GetName(),ROOT.gROOT.ProcessLine("TObject::kOverwrite"))
	

	fout.Close()	

if MAKE_HQ2W: #! then first make hq2w_on/of and then make_corr_factor_and_dgnstc_plts()
	#! Run Loop to make hq2w_on/off for rad=on and off in parallel
	#! + Used technique from http://stackoverflow.com/questions/7207309/python-how-can-i-run-python-functions-in-parallel
	qon=Queue()
	qof=Queue()
	pon=Process(target=Loop, args=(ON,qon))
	pon.start()
	pof=Process(target=Loop, args=(OFF,qof))
	pof.start()
	pon.join() # this blocks until the process terminates
	pof.join()
	hq2w_on=qon.get()
	hq2w_of=qof.get()
	#print "hq2w_on,hq2w_of=",hq2w_on.GetName(),hq2w_of.GetName()

	#! Write hq2w[rad] and related diagnostic plots to file
	fout=ROOT.TFile("%s/%s"%(OUTDIR,FOUT_NAME),"RECREATE")
	hq2w_on.Write(hq2w_on.GetName(),ROOT.gROOT.ProcessLine("TObject::kOverwrite"))
	hq2w_of.Write(hq2w_of.GetName(),ROOT.gROOT.ProcessLine("TObject::kOverwrite"))
	fout.Close()

	#! Now do make_corr_factor_and_dgnstc_plts()
	make_corr_factor_and_dgnstc_plts()
else:#! directly do 'make_corr_factor_and_dgnstc_plts()'
	make_corr_factor_and_dgnstc_plts()
