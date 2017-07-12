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

import datetime

'''
+ {e1f,e16}*{ER,SR}*2pi*{wzvtxcorr,nzvtxcorr}

+ Usage: study_zvtx.py expt=<e1f/e16> plots=ER-wzvtxcorr[:ER-nzvtxcorr:SR-nzvtxcorr]
'''

#! *** User entered parameters ***
if len(sys.argv)<2:
		sys.exit('usage: study_pcorr.py expt=(e1f/e16)')
expt=sys.argv[1]
if expt!='e1f' and expt!='e16':
	sys.exit('expt=%s is not valid. usage: study_pcorr.py expt=(e1f/e16)'%expt)
print "expt=",expt
#! For now ignore e1f
if expt=="e1f":
	sys.exit("Not implemented for e1f yet")

if len(sys.argv)>2: #! i.e. plots entered by user
	PLOTS=[item for item in sys.argv[2].split(":")]
else:
	PLOTS=["ER-wzvtxcorr","ER-nzvtxcorr","SR-nzvtxcorr"]

print "expt=",expt
print "PLOTS=",PLOTS
#sys.exit()

#! *** Prepare input data structure: FIN/T[dtyp][corr] ***
NDTYP=2
ER,SR=range(NDTYP)
DTYP_NAME=["ER","SR"]

NCORR=2
WZVTXCORR,NZVTXCORR=range(NCORR)
CORR_NAME=['wzvtxcorr','nzvtxcorr']
CORR_LEGEND_TITLE={'e1f':'e pcorr','e16':'e,p,#pi^{+} pcorr'}

FIN=[[[] for j in range(NCORR)] for i in range(NDTYP)]
T=  [[[] for j in range(NCORR)] for i in range(NDTYP)]
#print "FIN=",FIN
#print "T=",T
#! ***

#! *** Prepare output datadir***
if expt=="e1f":
	DATADIR_OUTPUT=os.environ['STUDY_ZVTX_E1F_DATADIR']
elif expt=="e16":
	DATADIR_OUTPUT=os.environ['STUDY_ZVTX_E16_DATADIR']
print "DATADIR_OUTPUT=",DATADIR_OUTPUT
#! ***

#! *** Now fill T[dtyp][rctn][corr] ***
DATADIR_INPUT=[[] for i in range(NDTYP)]
sfx=''
if expt=='e16':sfx='_E16'
DATADIR_INPUT[ER]=os.path.join(os.environ['D2PIDIR_EXP%s'%sfx],"data_zvtx_an_012617")#! debug /bk
DATADIR_INPUT[SR]=os.path.join(os.environ['D2PIDIR_SIM%s'%sfx],"data_zvtx_an_012617")#! debug /bl

print "DATADIR_INPUT[ER]=",DATADIR_INPUT[ER]
print "DATADIR_INPUT[SR]=",DATADIR_INPUT[SR]

d=[[ER,SR],[WZVTXCORR,NZVTXCORR]]
D=list(itertools.product(*d))
for r in D:
	idtyp,icorr=r[0],r[1]
	if DTYP_NAME[idtyp]=="SR" and CORR_NAME[icorr]=="wzvtxcorr": continue #! No {SR}*{wzvtxcorr}
	
	print "Getting T for",DTYP_NAME[idtyp],CORR_NAME[icorr]
	FIN[idtyp][icorr]=ROOT.TFile("%s/d2piR_%s.root"%(DATADIR_INPUT[idtyp],CORR_NAME[icorr]))
	T[idtyp][icorr]=FIN[idtyp][icorr].Get("d2piR/tR")
	#print FIN[idtyp][icorr].GetName()
	#print T[idtyp][icorr].GetName()
#print "T=",T
#sys.exit()
#! ***

#! *** Now make plots ***

#! sector information
NSCTR=6

#! *** Now make plots plots ***
#! outdir
DATE=datetime.datetime.now().strftime('%m%d%y')
outdir=os.path.join(DATADIR_OUTPUT,"results_%s"%DATE)
if not os.path.exists(outdir):
	os.makedirs(outdir)
#! draw_cmd
draw_cmd="vz_e>>hdcmd(100,-12,6)"
#! Create h[dtyp][sctr][corr]
h=[[[[] for k in range(NCORR)]for j in range(NSCTR)]for i in range(NDTYP)]
#print "h=",h
#! fill h
for r in D:
	idtyp,icorr=r[0],r[1]
	if DTYP_NAME[idtyp]=="SR" and CORR_NAME[icorr]=="wzvtxcorr": continue #! No {SR}*{wzvtxcorr}
	for isctr in range(NSCTR):
		cut=ROOT.TCut("sector==%d"%(isctr+1))
		print "TTree::Draw() for",DTYP_NAME[idtyp],CORR_NAME[icorr]
		#print T[idtyp][icorr].GetName()
		T[idtyp][icorr].Draw(draw_cmd,cut)#[idtyp][isctr][icorr],cut)
		#! Store histogram
		#ROOT.gStyle.SetOptStat("n")#! This line is critically needed here!
		htmp=ROOT.gDirectory.Get("hdcmd")#draw_cmd_hst_name[idtyp][isctr][icorr])
		h[idtyp][isctr][icorr]=htmp.Clone()
		h[idtyp][isctr][icorr].SetName("h_%s_%s"%(DTYP_NAME[idtyp],CORR_NAME[icorr]))
		h[idtyp][isctr][icorr].SetTitle("%s"%cut.GetTitle())
		#print h[idtyp][isctr][icorr].GetName()

#! Draw and save plots

#! Set up Cut TLines
#! + one for each (sector) pad on the canvas
#! exp
#! zvtx cut determined after Empty Target BG subtraction (etgt-bg-sub)
#! + This cut is only supposed to cut away the foil
#! + Value should be as coded in in $SUBSTUDIES/study_tgt_BG/e16/obtain_and_validate_R.py
#! + Latest plots in $STUDY_TGT_BG_E16_DATADIR/results_R_<latest-date>
ZVTX_CUTL_EI_EXP=[ROOT.TLine(-8.00,0,-8.00,0) for i in range(NSCTR)]
ZVTX_CUTH_EI_EXP=[ROOT.TLine(-0.75,0,-0.75,0) for i in range(NSCTR)] #! EI: (-0.8,0,-0.8,0)
#! [07-09-17] sim: since zvtx distribution for exp neq sim, sim cuts have to adapted
#! + sim cut is studied here for the first time
#! + bascically from ER by +0.5 offset by 0.5
#! + Note only exp cuts are studied in $SUBSTUDIES/study_tgt_BG/e16/obtain_and_validate_R.py
ZVTX_CUTL_EI_SIM=[ROOT.TLine(-7.50,0,-7.50,0) for i in range(NSCTR)]
ZVTX_CUTH_EI_SIM=[ROOT.TLine(-0.25,0,-0.25,0) for i in range(NSCTR)] #! 
for isctr in range(NSCTR):
	#! exp
	ZVTX_CUTL_EI_EXP[isctr].SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
	ZVTX_CUTH_EI_EXP[isctr].SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
	#! sim
	ZVTX_CUTL_EI_SIM[isctr].SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
	ZVTX_CUTH_EI_SIM[isctr].SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
	#! following same for exp and sim
	for line in [ZVTX_CUTL_EI_EXP[isctr],ZVTX_CUTH_EI_EXP[isctr],ZVTX_CUTL_EI_SIM[isctr],ZVTX_CUTH_EI_SIM[isctr]]:
		line.SetLineWidth(3)
		line.SetLineWidth(3)
		line.SetLineStyle(2)

ROOT.gStyle.SetOptStat("n")
ROOT.gStyle.SetOptFit(1)#111)
cname="czvtx"
czvtx=ROOT.TCanvas(cname,cname,2000,1500)
czvtx.Divide(3,2)
l=ROOT.TLegend(0.7,0.8,0.9,0.9)#,"","NDC");
for isctr in range(NSCTR):
	#cname="sector%d"%(isctr+1)
	#czvtx=ROOT.TCanvas(cname,cname)
	pad=czvtx.cd(isctr+1)
	if expt=="e1f":#! both {ER,SR} exist for e1f
		sys.exit("Not implemented for e1f for now")
	elif expt=="e16": #! only {ER} exist for e16
		h[ER][isctr][WZVTXCORR].SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
		h[ER][isctr][NZVTXCORR].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
		h[SR][isctr][NZVTXCORR].SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
		h[ER][isctr][WZVTXCORR].SetLineWidth(2)
		h[ER][isctr][NZVTXCORR].SetLineWidth(2)
		h[SR][isctr][NZVTXCORR].SetLineWidth(2)
		h[ER][isctr][WZVTXCORR].Draw("HIST")
		#! Axes title
		h[ER][isctr][WZVTXCORR].SetXTitle("z-vertex position [cm]")
		pad.SetLeftMargin(0.20)
		h[ER][isctr][WZVTXCORR].GetYaxis().SetTitleOffset(2.0)
		h[ER][isctr][WZVTXCORR].SetYTitle("N_{entries}")
		#! [01-22-17] Remove main title and stats box
		for htmp in [h[ER][isctr][WZVTXCORR],h[ER][isctr][NZVTXCORR],h[SR][isctr][NZVTXCORR]]:
			htmp.SetTitle("")
			htmp.SetStats(0)
		#! Draw rest of hists
		if "ER-nzvtxcorr" in PLOTS:
			h[ER][isctr][NZVTXCORR].Draw("HIST sames")
		#! scale SR and then draw
		if "SR-nzvtxcorr" in PLOTS:
               		max_ER= h[ER][isctr][WZVTXCORR].GetMaximum()
               		if max_ER==0:max_ER=1
               		max_SR=h[SR][isctr][NZVTXCORR].GetMaximum()
               		if max_SR==0:max_SR=1
               		scl_fctr_SR=max_ER/max_SR
               		if scl_fctr_SR==0:scl_fctr_SR=1
               		h[SR][isctr][NZVTXCORR].Scale(scl_fctr_SR)
			h[SR][isctr][NZVTXCORR].Draw("HIST sames")
		#! Ajust statbox
		# ROOT.gPad.Update()
		# #czvtx.Update()
		# #! statbox for wzvtxcorr
		# pt_wzvtxcorr=h[ER][isctr][WZVTXCORR].GetListOfFunctions().FindObject("stats")
  #               pt_wzvtxcorr.SetTextColor(ROOT.gROOT.ProcessLine("kBlue"))
		# #! statbox for nzvtxcorr
		# pt_nzvtxcorr=h[ER][isctr][NZVTXCORR].GetListOfFunctions().FindObject("stats")
  #               pt_nzvtxcorr.SetTextColor(ROOT.gROOT.ProcessLine("kBlack"))
  #               diff=pt_wzvtxcorr.GetY2NDC()-pt_wzvtxcorr.GetY1NDC()
  #               y2=pt_wzvtxcorr.GetY1NDC()
  #               y1=y2-diff
  #               pt_nzvtxcorr.SetY1NDC(y1)
  #               pt_nzvtxcorr.SetY2NDC(y2)
		# #czvtx.Update()
		# #! statbox for SR:nzvtxcorr
		# if "SR-nzvtxcorr" in PLOTS:
  #               	pt_SR_nzvtxcorr=h[SR][isctr][NZVTXCORR].GetListOfFunctions().FindObject("stats")
  #               	pt_SR_nzvtxcorr.SetTextColor(ROOT.gROOT.ProcessLine("kRed"))
  #               	diff=pt_nzvtxcorr.GetY2NDC()-pt_nzvtxcorr.GetY1NDC()
  #               	y2=pt_nzvtxcorr.GetY1NDC()
  #               	y1=y2-diff
  #               	pt_SR_nzvtxcorr.SetY1NDC(y1)
  #               	pt_SR_nzvtxcorr.SetY2NDC(y2)
  #               #czvtx.Update()
		# #! Draw cut lines
		cutlines=[ZVTX_CUTL_EI_EXP[isctr],ZVTX_CUTH_EI_EXP[isctr]]
		if "SR-nzvtxcorr" in PLOTS: 
			cutlines+=[ZVTX_CUTL_EI_SIM[isctr],ZVTX_CUTH_EI_SIM[isctr]]
		for line in cutlines:#[ZVTX_CUTL_EI_EXP[isctr],ZVTX_CUTH_EI_EXP[isctr],ZVTX_CUTL_EI_SIM[isctr],ZVTX_CUTH_EI_SIM[isctr]]:
			line.SetY1(h[ER][isctr][WZVTXCORR].GetMinimum())
			line.SetY2(h[ER][isctr][WZVTXCORR].GetMaximum())
			line.Draw("same")
		# ZVTX_CUTL_EI[isctr].SetY1(h[ER][isctr][WZVTXCORR].GetMinimum())
		# ZVTX_CUTL_EI[isctr].SetY2(h[ER][isctr][WZVTXCORR].GetMaximum())
		# ZVTX_CUTH_EI[isctr].SetY1(h[ER][isctr][WZVTXCORR].GetMinimum())
		# ZVTX_CUTH_EI[isctr].SetY2(h[ER][isctr][WZVTXCORR].GetMaximum())
		# ZVTX_CUTL_EI[isctr].Draw("same")
		# ZVTX_CUTH_EI[isctr].Draw("same")
		#! Draw legend #! (not necessary because statbox contain relevant information)
		#! + Note legend entries named differently when demoing ER-SR comparison and zvtx correction
		if isctr+1==1:
			if "SR-nzvtxcorr" in PLOTS: #! ER-SR comparison 
				l.AddEntry(h[ER][isctr][WZVTXCORR],'exp-z-vtx')
				l.AddEntry(h[SR][isctr][NZVTXCORR],'sim-z-vtx')
			else: #! zvtx correction
				l.AddEntry(h[ER][isctr][WZVTXCORR],'pst-z-vtx-corr')
				if "ER-nzvtxcorr" in PLOTS:
					l.AddEntry(h[ER][isctr][NZVTXCORR],'pre-z-vtx-corr')
			l.Draw()
czvtx.Update()
czvtx.SaveAs("%s/%s_%s.png"%(outdir,cname,'_'.join(PLOTS)))
czvtx.SaveAs("%s/%s_%s.pdf"%(outdir,cname,'_'.join(PLOTS)))
#! ***
#! debug
#sys.exit()
	

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
	
