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

FNAME=[0 for i in range(NSIMRNG)]
FNAME[L]="%s/SS/lowQ2_cmb_vst_SE_080317/Obs_itg.root"%(os.environ['OBSDIR_E16'])  #!092716
FNAME[H]="%s/SS/highQ2_cmb_vst_SE_080317/Obs_itg.root"%(os.environ['OBSDIR_E16']) #!092716

#! Set up Q2-W bins that are in the respective files
Q2BINS=[0 for i in range(NSIMRNG)]
Q2BINS[L]=["2.00-2.40_1.400-2.125","2.40-3.00_1.400-2.125"]
Q2BINS[H]=["3.00-3.50_1.400-2.125","3.50-4.20_1.400-2.125","4.20-5.00_1.400-2.125"]
NQ2BINS=len(Q2BINS[L])+len(Q2BINS[H])

#! Create object structure to be created
#! +h[q2][seq] where seq={EC,EF}
#! +hdiff[q2] where diff=((EF-EC)/EC)*100
NSEQS=2
EC,EF=range(2)
h=[[0 for j in range(NSEQS)] for i in range(NQ2BINS)]
hdiff=[0 for i in range(NQ2BINS)]
#print h
#! OUTDIR
OUTDIR="%s/figures/Holes/"%os.environ['ANANOTE']
# #! For debugging
# OUTDIR="/data1/trivedia/tmp"

#! Now get object from file and make output objects
#! Create a file structure for each .root file to be opened (else there are technical problems)
f=[0 for i in range(NSIMRNG)]
iq2=0
for isim in range(NSIMRNG):
	f[isim]=ROOT.TFile(FNAME[isim])
	#print f[isim].GetName()
	for q2 in Q2BINS[isim]:
		#print "%s/hW_obs_%s_THETA"%(q2,'EC')
		h[iq2][EC]=f[isim].Get("%s/hW_obs_%s_THETA"%(q2,'EC'))
		h[iq2][EF]=f[isim].Get("%s/hW_obs_%s_THETA"%(q2,'EF'))
		#! Calculate hdiff
		#! First call Sumw2()
		h[iq2][EC].Sumw2()
		h[iq2][EF].Sumw2()
		hdiff[iq2]=h[iq2][EF].Clone("hdiff")
		hdiff[iq2].Sumw2()
		hdiff[iq2].Add(h[iq2][EC],-1)
		hdiff[iq2].Scale(1.0/2.0)
		hdiff[iq2].Divide(h[iq2][EF])
		hdiff[iq2].Scale(100)
		#print h[iq2][EC].GetName(), h[iq2][EF].GetName()
		iq2+=1

#! Draw all Q2 bins on one canvas
CLRS_IQ2={(0):ROOT.gROOT.ProcessLine("kRed"),
		  (1):ROOT.gROOT.ProcessLine("kGreen"),
		  (2):ROOT.gROOT.ProcessLine("kCyan"),
		  (3):ROOT.gROOT.ProcessLine("kBlue"),
		  (4):ROOT.gROOT.ProcessLine("kMagenta")
		 }
MRKS={(EC):ROOT.gROOT.ProcessLine("kOpenTriangleDown"),#kFullDiamond
	  (EF):ROOT.gROOT.ProcessLine("kFullDotLarge")
	 }
		  
CWDTH=800
CHGHT=1000
c=ROOT.TCanvas("c","c",CWDTH,CHGHT)
c.Divide(1,2)
ROOT.gStyle.SetOptStat(0)
#! legend
l1=ROOT.TLegend(0.50,0.75,0.9,0.9)
l1.SetNColumns(2)
l1.SetFillStyle(0)
l1.SetBorderSize(1)
l1.SetTextSize(0.03)#0.02
#l1.SetHeader("#splitline{Q^{2} bins:}{#bullet: After hole-filling}{#nabla: Before hole-filling}")#,"C")
#l1.SetHeader("#splitline{#splitline{#bullet: After hole-filling}{#nabla: Before hole-filling}}{Q^{2} bins:}")
#l1.SetHeader("#splitline{#bullet: After hole-filling, #nabla: Before hole-filling}{Q^{2} bins:}")
#l1.SetHeader("#bullet: After hole-filling, #nabla: Before hole-filling. Q^{2} bins:")
#l1.SetHeader("#splitline{Integrated cross-sections before (#nabla) and after (#bullet) hole-filling}{Q^{2} bins:}")
l1.SetHeader("Q^{2} bins:")
for iq2 in range(NQ2BINS):
	draw_opt="same"
	if iq2==0: draw_opt=""
	pad=c.cd(1)
	#! To give rotate indices space
	pad.SetBottomMargin(0.20)
	#! Titles
	h[iq2][EC].SetTitle("Integrated Cross-Sections before (#nabla) and after (#bullet) Hole Filling")
	h[iq2][EF].SetTitle("Integrated Cross-Sections before (#nabla) and after (#bullet) Hole Filling")

	h[iq2][EC].SetMarkerStyle(MRKS[EC])
	h[iq2][EF].SetMarkerStyle(MRKS[EF])
	h[iq2][EC].SetLineColor(CLRS_IQ2[iq2])
	h[iq2][EC].SetMarkerColor(CLRS_IQ2[iq2])
	h[iq2][EF].SetLineColor(CLRS_IQ2[iq2])
	h[iq2][EF].SetMarkerColor(CLRS_IQ2[iq2])
	h[iq2][EC].Draw(draw_opt)
	h[iq2][EF].Draw("same")
	#! legend
	if iq2<=1: 	q2bin=Q2BINS[L][iq2].split("_")[0]
	elif iq2>1: q2bin=Q2BINS[H][iq2-2].split("_")[0]
	l1.AddEntry(h[iq2][EF],q2bin,"l")
	
	pad=c.cd(2)
	#! Titles
	hdiff[iq2].SetTitle("Relative Systematic Uncertainty (%) to Integrated Cross-Sections from Hole Filling")
	pad.SetBottomMargin(0.20)
	#! Y-axis title: Remove \microbarns
	hdiff[iq2].SetYTitle("")
	hdiff[iq2].SetLineColor(CLRS_IQ2[iq2])
	hdiff[iq2].SetMarkerColor(CLRS_IQ2[iq2])
	hdiff[iq2].Draw(draw_opt)
#! legend on 1st and 2nd pad	
c.cd(1)
l1.Draw()
c.cd(2)
l1.Draw()
c.SaveAs("%s/%s.png"%(OUTDIR,"HF_EC_EF_itg"))
c.SaveAs("%s/%s.pdf"%(OUTDIR,"HF_EC_EF_itg"))

#! Plot values of hdiff[iq2] in a single histogram to get average
hdiff_hist=ROOT.TH1F("hdiff_hist","hdiff_hist",20,0,20)
hdiff_hist.SetTitle("#splitline{Distribution of Relative Systematic Uncertainty (%) to}{       Cross-Sections from Hole Filling}")
#! Not sure why, but the following line is needed to get stats to show after ROOT.gStyle.SetOptStat(0)
#! + I tried ROOT.gStyle.Reset(), but that did not seem to work
#ROOT.gStyle.Reset()
hdiff_hist.SetStats(1) 
ROOT.gStyle.SetStatY(0.9)
ROOT.gStyle.SetStatX(0.9)
ROOT.gStyle.SetStatW(0.2)
ROOT.gStyle.SetStatH(0.3) 
ROOT.gStyle.SetOptStat("mr") #! ("nemruo")
hdiff_hist.SetXTitle("Relative Systematic Uncertainty (%)")
for iq2 in range(NQ2BINS):
	nbins=hdiff[iq2].GetNbinsX()
	for ibin in range(nbins):
		binc=hdiff[iq2].GetBinContent(ibin+1)
		hdiff_hist.Fill(binc)
c1=ROOT.TCanvas("c1","c1",CWDTH,CHGHT)
hdiff_hist.Draw("hist")
c1.Update()
c1.SaveAs("%s/%s.png"%(OUTDIR,"rel_syst_err_HF_itg"))
c1.SaveAs("%s/%s.pdf"%(OUTDIR,"rel_syst_err_HF_itg"))

print "make_quantify_hole_filling.py done"

# c1=ROOT.TCanvas("c1","c1")
# m=ROOT.TMarker(0.5,0.8,ROOT.gROOT.ProcessLine("kFullTriangleUp"))
# m.SetMarkerColor(ROOT.gROOT.ProcessLine("kRed"))
# m.Draw()
# pt=ROOT.TPaveText(.05,.1,.95,.8)
# pt.AddText("#bullet")
# pt.AddText("#diamond")
# pt.SetTextSize(0.32)
# pt.Draw()

# #! If wanting to keep TCanvas open till program exits                           
# if not ROOT.gROOT.IsBatch():
# 	plt.show()
# 	# wait for you to close the ROOT canvas before exiting
# 	wait(True)
