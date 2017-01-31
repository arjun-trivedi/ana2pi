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

import atlib

'''
+ This script makes plots to be used in the thesis:
1. Simply consolidate plots already produced in .root files in '$STUDY_EID_NPHE_E16_DATADIR/study_nphe_cut.root'
   and saves them directly as .pdf in '$STUDY_EID_NPHE_E16_DATADIR/thesis_plots' that can be copied over to the thesis folder.

2. Histogram cut and efficiency values for the loose cuts and tight cuts
'''
USAGE="make_plots_thesis_nphe_cut.py"

#! *** Prepare DATADIR and OUTDIR***
DATADIR=os.path.join(os.environ['STUDY_EID_NPHE_E16_DATADIR'])
print "DATADIR=",DATADIR

OUTDIR=os.path.join(DATADIR,"thesis_plots")
if not os.path.exists(OUTDIR):
		os.makedirs(OUTDIR)
print "OUTDIR=",OUTDIR
#! ***

#! *** Setup some constants
NSCTRS=6
NSGMNTS=18
NPMTS=3
IPMT_L,IPMT_C,IPMT_R=range(NPMTS)
#! ***

# #! 1. Consolidate plots already produced in .root files in '$STUDY_EID_NPHE_E16_DATADIR/study_nphe_cut.root'
# #! *** Set up files to get data 
# F=ROOT.TFile("%s/study_nphe_cut.root"%DATADIR)
# print "F=",F.GetName()
# #! ***

# #! Now get relevant theta vs p TCanvas from the files and save them
# for isctr in range(NSCTRS):
# 	sctr=isctr+1
# 	#! Get TCanvas for each PMT
# 	c=[0 for i in range(NPMTS)]
# 	for ipmt in range(NPMTS):
# 		pmt=ipmt+1
# 		c[ipmt]=F.Get("sector%d/sct%d_pmt%d"%(sctr,sctr,pmt))
# 		print c[ipmt].GetName()
# 	#! Now save TCanvas
# 	for ipmt in range(NPMTS):
# 		pmt=ipmt+1
# 		c[ipmt].Draw()
# 		c[ipmt].SaveAs("%s/nphe_cut_sctr%d_pmt%d.pdf"%(OUTDIR,sctr,ipmt+1))
# 		c[ipmt].SaveAs("%s/nphe_cut_sctr%d_pmt%d.png"%(OUTDIR,sctr,ipmt+1))

#! If wanting to keep TCanvas open till program exits                           
#if not ROOT.gROOT.IsBatch():
#        plt.show()
#        # wait for you to close the ROOT canvas before exiting
#        wait(True)
			
#! 2. Histogram cut and efficiency values for the loose cuts and tight cuts
#! Read cut, eff for lse,tgt
NCUTLVLS=2 #! 0,1=lse,tgt
CUTLVLS=['lse','tgt']
fcut,feff=[0 for i in range(NCUTLVLS)],[0 for i in range(NCUTLVLS)]
for icutlvl,cutlvl in enumerate(CUTLVLS):
	fcut[icutlvl]=open("%s/eid_e16_exp_nphe_cut_%s.txt"%(DATADIR,cutlvl),"r")
	feff[icutlvl]=open("%s/eid_e16_exp_nphe_eff_%s.txt"%(DATADIR,cutlvl),"r")

#! Create structure to hold file data
#! cut[cutlvl][sct][sgm]
cut=[[[[]for k in range(NSGMNTS)]for j in range(NSCTRS)]for i in range(NCUTLVLS)]
#! eff[cutlvl][sct][sgm][pmt]
#! not 2 indices for PMTS: 1,3=L,R
eff=[[[[[]for l in range(2)]for k in range(NSGMNTS)]for j in range(NSCTRS)]for i in range(NCUTLVLS)]

#! Read file data into structure
for icutlvl in range(NCUTLVLS):
	for line in fcut[icutlvl].readlines():
		#print line
		words=line.split(" ")
		sct,sgm,val=int(words[0]),int(words[1]),float(words[2])
		#print "sct,sgm,val=",sct,sgm,val
		cut[icutlvl][sct-1][sgm-1]=val
	for line in feff[icutlvl].readlines():
		#print line
		words=line.split(" ")
		sct,sgm,pmt,val=int(words[0]),int(words[1]),int(words[2]),float(words[3])
		if   pmt==1: ipmt=0
		elif pmt==3: ipmt=1
		eff[icutlvl][sct-1][sgm-1][ipmt]=val

#! Fill hcut, heff histograms
hcut=[[0 for j in range(NSCTRS)] for i in range(NCUTLVLS)]
heff=[[[0 for k in range(2)]for j in range(NSCTRS)]for i in range(NCUTLVLS)]
for icutlvl in range(NCUTLVLS):
	for isct in range(NSCTRS):
		#! create hcut
		hname_cut="hcut_l%s_s%d"%(CUTLVLS[icutlvl],isct+1)
		htitle=""
		hcut[icutlvl][isct]=ROOT.TH1F(hname_cut,htitle,NSGMNTS,1-0.5,18+0.5)
		hcut[icutlvl][isct].SetXTitle("segment")
		hcut[icutlvl][isct].SetYTitle("cut-value [nphe x 10]")
		for isgm in range(NSGMNTS):
			hcut[icutlvl][isct].SetBinContent(isgm+1,cut[icutlvl][isct][isgm])
			hcut[icutlvl][isct].SetBinError(isgm+1,0)
			hcut[icutlvl][isct].SetMinimum(0)
		#! create heff
		for ipmt in range(2):
			hname_eff="heff_l%s_s%d_p%d"%(CUTLVLS[icutlvl],isct+1,pmt)
			htitle=""
			heff[icutlvl][isct][ipmt]=ROOT.TH1F(hname_eff,htitle,NSGMNTS,1-0.5,18+0.5)
			heff[icutlvl][isct][ipmt].SetXTitle("segment")
			heff[icutlvl][isct][ipmt].SetYTitle("efficiency [%]")
			for isgm in range(NSGMNTS):
				heff[icutlvl][isct][ipmt].SetBinContent(isgm+1,eff[icutlvl][isct][isgm][ipmt])
				heff[icutlvl][isct][ipmt].SetBinError(isgm+1,0)
#! plot cut,eff
ccut=ROOT.TCanvas("ccut","ccut",2000,1500)
ceff=ROOT.TCanvas("ceff","ceff",2000,1500)
ccut.Divide(3,2)
ceff.Divide(3,2)
#! legends
lcut=ROOT.TLegend(0.5,0.8,0.6,0.9)#(0.8,0.8,0.9,0.9)
leff=ROOT.TLegend(0.5,0.4,0.7,0.6)#(0.8,0.8,0.9,0.9)
for isct in range(NSCTRS):
	#! draw hcut
	p=ccut.cd(isct+1)
	#! set pad aesthetics
	p.SetGridx()
	p.SetLeftMargin(0.20)
	#! set hcut aesthetics
	for icutlvl in range(NCUTLVLS):
		hcut[icutlvl][isct].GetYaxis().SetTitleOffset(1.5)
		hcut[icutlvl][isct].SetStats(0)
		hcut[icutlvl][isct].SetLineColor(icutlvl+1)
	#! draw
	for icutlvl in [1,0]: #! draw tight first to that upper limit is set range(NCUTLVLS):
		draw_opt=""
		if icutlvl==0:draw_opt="same"
		hcut[icutlvl][isct].Draw(draw_opt)
		#! legend
		if isct==0:
			lcut.AddEntry(hcut[icutlvl][isct],CUTLVLS[icutlvl])
	#! Draw legend
	if isct==0:
		lcut.Draw()

	#! draw heff: both pmts on same canvas
	p=ceff.cd(isct+1)
	#! set pad aesthetics
	p.SetGridx()
	p.SetLeftMargin(0.20)
	#! set heff aesthetics
	for icutlvl in range(NCUTLVLS):
		for ipmt in range(2):
			heff[icutlvl][isct][ipmt].SetStats(0)
			heff[icutlvl][isct][ipmt].GetYaxis().SetTitleOffset(1.5)
			heff[icutlvl][isct][ipmt].SetLineColor(icutlvl+1)
			heff[icutlvl][isct][ipmt].SetLineStyle(ipmt+1)
	#! draw
	for icutlvl in range(NCUTLVLS):
		for ipmt in range(2):
			draw_opt=""
			if icutlvl>0 or ipmt>0: draw_opt="same"
			heff[icutlvl][isct][ipmt].Draw(draw_opt)
			if isct==0:
				if   ipmt==0: leff.AddEntry(heff[icutlvl][isct][ipmt],'%s-pmt-L'%CUTLVLS[icutlvl])
				elif ipmt==1: leff.AddEntry(heff[icutlvl][isct][ipmt],'%s-pmt-R'%CUTLVLS[icutlvl])
	if isct==0: leff.Draw()
#! draw and save canvas
ccut.SaveAs("%s/cut.png"%OUTDIR)
ccut.SaveAs("%s/cut.pdf"%OUTDIR)
ccut.Close()
ceff.SaveAs("%s/cut_eff.png"%OUTDIR)
ceff.SaveAs("%s/cut_eff.pdf"%OUTDIR)
ceff.Close()


