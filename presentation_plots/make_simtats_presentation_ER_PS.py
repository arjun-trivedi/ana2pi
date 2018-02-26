#!/usr/bin/python

import os,sys
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict

import numpy as np

from presentation_plots_lib import *

'''
[02-26-18]
Demonstrate typical simstats in ER PS using SSBands:cutsncors1 (gpart-PID*stat*PID=off*off)
'''

FNAME=[0 for i in range(NSIMRNG)]
FNAME[L]="%s/lowQ2_SSBands_122717/cutsncors1/sim4_sim5_sim6_sim7_sim8_sim13/simstats_ER_PS/simstats_ER_PS.root"%(os.environ['OBSDIR_E16'])
FNAME[H]="%s/highQ2_SSBands_122717/cutsncors1/sim9_sim10_sim11_sim12/simstats_ER_PS/simstats_ER_PS.root"%(os.environ['OBSDIR_E16'])

#! Setup names of hists to get from simstats.root
SIMSTAT_HISTS=["h_mu_ST",                           "h_mu_SR",                           "h_mu_SA"]
NAME_HISTS   =["h_mu_ST_ER_PS",                     "h_mu_SR_ER_PS",                     "h_mu_SA_ER_PS"]
TITLE_HISTS  =["#LTST^{5}#GT(Q^{2},W) in ER^{5} PS","#LTSR^{5}#GT(Q^{2},W) in ER^{5} PS","#LTSA^{5}#GT(Q^{2},W) in ER^{5} PS"]
#nbins-SH^{5}/nbins-ST^{5}
#! OUTDIR
OUTDIR="%s/figures/Acceptance/"%os.environ['ANANOTE']

for ihist,hist in enumerate(SIMSTAT_HISTS):
	#! Get simstat hist from each file
	h=[0 for i in range(NSIMRNG)]
	f=[0 for i in range(NSIMRNG)]
	print h
	for i in range(NSIMRNG):
		f[i]=ROOT.TFile(FNAME[i])
		h[i]=f[i].Get(hist)
	#! Add simstat hist from each file
	hf=h[L].Clone("full")
	hf.Add(h[H])

	#! Setup aesthetics
	#! Set maximum
	mxm=hf.GetBinContent(hf.GetMaximumBin())
	hf.SetMaximum(mxm)
	#! Adjust titles
	hf.SetTitle(TITLE_HISTS[ihist])
	#! Fix y-axis title since it was not set up correctly in 'ds_lite'
	hf.SetYTitle("Q^{2} [GeV^{2}]")

	c=ROOT.TCanvas()
	ROOT.gStyle.SetOptStat(0)
	hf.Draw("colz")
	#! Draw Q2-W grid
	lw=[ROOT.TLine() for j in range(len(WBINL))]
	lq=[ROOT.TLine() for j in range(len(Q2BINL))]
	for j,wbin in enumerate(WBINL):
		lw[j].DrawLine(wbin,Q2MIN,wbin,Q2MAX)
	for j,qbin in enumerate(Q2BINL):
		lq[j].DrawLine(WMIN,qbin,WMAX,qbin)

	c.Draw()
	c.SaveAs("%s/%s.png"%(OUTDIR,NAME_HISTS[ihist]))
	c.SaveAs("%s/%s.pdf"%(OUTDIR,NAME_HISTS[ihist]))

print "make_simtats_presentation_ER_PS.py done"


#! If wanting to keep TCanvas open till program exits                           
if not ROOT.gROOT.IsBatch():
	plt.show()
	# wait for you to close the ROOT canvas before exiting
	wait(True)




