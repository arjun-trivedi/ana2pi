#!/usr/bin/python

import os,sys
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict

import numpy as np

from presentation_plots_lib import *

'''
[10-18-16]

'''
F=ROOT.TFile("%s/rad-eff-corr-sim/sim1/drad.root"%os.environ['D2PIDIR_SIM_E16'])

#! Set up Q2-W bin dict as per q2 bins in file 
Q2BINS=["2.00-2.40","2.40-3.00","3.00-3.50","3.50-4.20","4.20-5.00"]
NQ2BINS=len(Q2BINS)

#! Create object structure to be created
#! +h[q2]
h=[0 for i in range(NQ2BINS)]
#print h
#! OUTDIR
OUTDIR="%s/figures/Radiative/"%os.environ['THESIS']

#! Now get object from file and make output objects
for i in range(NQ2BINS):
	h[i]=F.Get("hcf_qbin_%d"%(i+1))

#! Draw all Q2 bins on one canvas
CLRS_IQ2={(0):ROOT.gROOT.ProcessLine("kRed"),
		  (1):ROOT.gROOT.ProcessLine("kGreen"),
		  (2):ROOT.gROOT.ProcessLine("kCyan"),
		  (3):ROOT.gROOT.ProcessLine("kBlue"),
		  (4):ROOT.gROOT.ProcessLine("kMagenta")
		 }
CWDTH=600
CHGHT=500
c=ROOT.TCanvas("c","c",CWDTH,CHGHT)
ROOT.gStyle.SetOptStat(0)
#! legend
l1=ROOT.TLegend(0.50,0.75,0.9,0.9)
l1.SetNColumns(3)
l1.SetFillStyle(0)
l1.SetBorderSize(1)
l1.SetTextSize(0.03)#0.02
l1.SetHeader("Q^{2} bins:")
for i in range(NQ2BINS):
	draw_opt="same"
	if i==0: draw_opt=""
	#! Set title
	h[i].SetTitle("Radiative Effects Correction Factor (1/R)")
	#! Marker and color
	h[i].SetMarkerStyle(ROOT.gROOT.ProcessLine("kFullDotLarge"))
	h[i].SetMarkerColor(CLRS_IQ2[i])
	h[i].SetLineColor(CLRS_IQ2[i])
	#! Draw
	h[i].Draw(draw_opt)
	#! Add to legend
	l1.AddEntry(h[i],Q2BINS[i],"l")
#! Draw legend 
l1.Draw()

c.SaveAs("%s/%s.png"%(OUTDIR,"rad_eff_corr_fctr"))
c.SaveAs("%s/%s.pdf"%(OUTDIR,"rad_eff_corr_fctr"))

print "make_quantify_hole_filling.py done"

#! If wanting to keep TCanvas open till program exits                           
# if not ROOT.gROOT.IsBatch():
# 	plt.show()
# 	# wait for you to close the ROOT canvas before exiting
# 	wait(True)