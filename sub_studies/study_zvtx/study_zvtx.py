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

'''
+ {e1f,e16}*{ER,SR}*2pi*{wzvtxcorr,nzvtxcorr}

+ Usage: study_zvtx.py expt=<e1f/e16> plots=[ER-wzvtxcorr:ER-nzvtxcorr:SR-nzvtxcorr]
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
DATADIR_INPUT[ER]=os.path.join(os.environ['D2PIDIR_EXP%s'%sfx],"data_zvtx_011516")#! debug /bk
DATADIR_INPUT[SR]=os.path.join(os.environ['D2PIDIR_SIM%s'%sfx],"data_zvtx_011516")#! debug /bl

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
outdir=os.path.join(DATADIR_OUTPUT,"results")
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
ZVTX_CUTL_EI=ROOT.TLine(-8.0,0,-8.0,0);
ZVTX_CUTH_EI=ROOT.TLine(-0.8,0,-0.8,0);
ZVTX_CUTL_EI.SetLineColor(ROOT.gROOT.ProcessLine("kGreen"))
ZVTX_CUTH_EI.SetLineColor(ROOT.gROOT.ProcessLine("kGreen"))
ZVTX_CUTL_EI.SetLineWidth(3)
ZVTX_CUTH_EI.SetLineWidth(3)

ROOT.gStyle.SetOptStat("n")
ROOT.gStyle.SetOptFit(1)#111)
cname="czvtx"
czvtx=ROOT.TCanvas(cname,cname,2000,1500)
czvtx.Divide(3,2)
for isctr in range(NSCTR):
	#cname="sector%d"%(isctr+1)
	#czvtx=ROOT.TCanvas(cname,cname)
	czvtx.cd(isctr+1)
	#l=ROOT.TLegend(0.1,0.3,0.3,0.4)#,"","NDC");
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
		ROOT.gPad.Update()
		#czvtx.Update()
		#! statbox for wzvtxcorr
		pt_wzvtxcorr=h[ER][isctr][WZVTXCORR].GetListOfFunctions().FindObject("stats")
                pt_wzvtxcorr.SetTextColor(ROOT.gROOT.ProcessLine("kBlue"))
		#! statbox for nzvtxcorr
		pt_nzvtxcorr=h[ER][isctr][NZVTXCORR].GetListOfFunctions().FindObject("stats")
                pt_nzvtxcorr.SetTextColor(ROOT.gROOT.ProcessLine("kBlack"))
                diff=pt_wzvtxcorr.GetY2NDC()-pt_wzvtxcorr.GetY1NDC()
                y2=pt_wzvtxcorr.GetY1NDC()
                y1=y2-diff
                pt_nzvtxcorr.SetY1NDC(y1)
                pt_nzvtxcorr.SetY2NDC(y2)
		#czvtx.Update()
		#! statbox for SR:nzvtxcorr
		if "SR-nzvtxcorr" in PLOTS:
                	pt_SR_nzvtxcorr=h[SR][isctr][NZVTXCORR].GetListOfFunctions().FindObject("stats")
                	pt_SR_nzvtxcorr.SetTextColor(ROOT.gROOT.ProcessLine("kRed"))
                	diff=pt_nzvtxcorr.GetY2NDC()-pt_nzvtxcorr.GetY1NDC()
                	y2=pt_nzvtxcorr.GetY1NDC()
                	y1=y2-diff
                	pt_SR_nzvtxcorr.SetY1NDC(y1)
                	pt_SR_nzvtxcorr.SetY2NDC(y2)
                #czvtx.Update()
		#! Draw cut lines
		ZVTX_CUTL_EI.SetY1(h[ER][isctr][WZVTXCORR].GetMinimum())
                ZVTX_CUTL_EI.SetY2(h[ER][isctr][WZVTXCORR].GetMaximum())
                ZVTX_CUTH_EI.SetY1(h[ER][isctr][WZVTXCORR].GetMinimum())
                ZVTX_CUTH_EI.SetY2(h[ER][isctr][WZVTXCORR].GetMaximum())
                ZVTX_CUTL_EI.Draw("same")
                ZVTX_CUTH_EI.Draw("same")
		#! Draw legend (not necessary because statbox contain relevant information)
		#l.AddEntry(h[ER][isctr][WZVTXCORR],'wzvtxcorr')
		#l.AddEntry(h[ER][isctr][NZVTXCORR],'nzvtxcorr')
		#l.Draw()
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
	
