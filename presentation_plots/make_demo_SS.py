#!/usr/bin/python

import os,sys
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict

import numpy as np

from presentation_plots_lib import *

'''
[09-29-17]

+ This script makes plots to demonstrate "statistically driven systematic uncertainties"
  as noted in the analysis noted by example of Variable Dependent Integrated Cross-Sections
'''

#FIN=ROOT.TFile("%s/SS/lowQ2_cmb_non_vst_SE_080317/Obs_itg.root"%(os.environ['OBSDIR_E16']),"READ")
FIN=ROOT.TFile("%s/SS/lowQ2_cmb_non_vst_SE_122817/Obs_itg.root"%(os.environ['OBSDIR_E16']),"READ")

#! OUTDIR
OUTDIR="%s/figures/Results/"%os.environ['ANANOTE']

NVSTS=3
VSTNAME=["var-set-1","var-set-2","var-set-3"]

#! Now get objects from file and make output objects
h=[0 for i in range(NVSTS)]
h[0]=FIN.Get("2.00-2.40_1.400-2.125/hW_obs_EF_1_THETA")
h[1]=FIN.Get("2.00-2.40_1.400-2.125/hW_obs_EF_2_THETA")
h[2]=FIN.Get("2.00-2.40_1.400-2.125/hW_obs_EF_3_THETA")

#! hist aesthetics
CLRS=[ROOT.gROOT.ProcessLine("kRed"),ROOT.gROOT.ProcessLine("kMagenta"),ROOT.gROOT.ProcessLine("kGreen+3")]
for i in range(NVSTS):
	h[i].SetLineColor(CLRS[i])
	h[i].SetMarkerColor(CLRS[i])

#! Draw and save
c=ROOT.TCanvas("c","c")
ROOT.gStyle.SetOptStat(0)
#! Give canvas a title and plot pad
pad_t=ROOT.TPad("pad_l","Legend pad",0.01,0.935,0.99,1.00)
pad_p=ROOT.TPad("pad_p","Plots pad", 0.01,0.935,0.99,0.01)
pad_t.Draw()
pad_p.Draw()
#! Set canvas title
pad_t.cd()
pt=ROOT.TPaveText(.01,.01,.99,.99)
#! Note the deliberate spacing between superscripts for Q2 and GeV below
pt.AddText("Variable Set Dependent Integrated Cross-Section for Q   ^{2}=[2.00, 2.40) GeV ^{2}")
pt.SetTextSize(0.65)
pt.Draw()
#! Now draw plot
#! To give rotate indices space
pad_p.cd()
pad_p.SetBottomMargin(0.20)
#! legend
l=ROOT.TLegend(0.7,0.75,0.9,0.9)
l.SetFillStyle(0)
l.SetBorderSize(1)
l.SetTextSize(0.03)#0.02
for i in range(NVSTS):
	draw_opt=""
	if i>0:draw_opt="same"
        h[i].Draw(draw_opt)
	l.AddEntry(h[i],"%s"%VSTNAME[i],"lp")
l.Draw()
c.SaveAs("%s/%s.png"%(OUTDIR,"demo_SS"))
c.SaveAs("%s/%s.pdf"%(OUTDIR,"demo_SS"))


# #! If wanting to keep TCanvas open till program exits                           
# if not ROOT.gROOT.IsBatch():
# 	plt.show()
# 	# wait for you to close the ROOT canvas before exiting
# 	wait(True)
