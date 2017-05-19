#!/usr/bin/python
from __future__ import division

import os,sys
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict
import itertools
import datetime

import numpy as np

import math

'''
+ Usage: study_tgt.py
'''

#! *** User entered parameters ***
#sys.exit()

#! *** Prepare input data structure: FIN/T[dtyp] ***
NDTYP=2
PTGT,ETGT=range(NDTYP)
DTYP_NAME=["PTGT","ETGT"]

FIN=[0 for i in range(NDTYP)]
T=  [0 for i in range(NDTYP)]
#print "FIN=",FIN
#print "T=",T
#! ***

#! *** Prepare output datadir***
DATADIR_OUTPUT=os.environ['STUDY_TGT_BG_E16_DATADIR']
print "DATADIR_OUTPUT=",DATADIR_OUTPUT
#! ***

#! *** Now fill T[dtyp][rctn][corr] ***
DATADIR_INPUT=os.path.join(os.environ['D2PIDIR_EXP_E16'],"data_tgt_051817")
print "DATADIR_INPUT",DATADIR_INPUT

for idtyp in range(NDTYP):
	if       DTYP_NAME[idtyp]=="PTGT": FIN[idtyp]=ROOT.TFile("%s/dtgt_ptgt.root"%DATADIR_INPUT)
	elif     DTYP_NAME[idtyp]=="ETGT": FIN[idtyp]=ROOT.TFile("%s/dtgt_etgt.root"%DATADIR_INPUT)
	
	print "Getting T for",DTYP_NAME[idtyp]
	T[idtyp]=FIN[idtyp].Get("d2piR/tR")
	#print FIN[idtyp].GetName()
	#print T[idtyp].GetName()
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
#outdir="/tmp/results" #! test
if not os.path.exists(outdir):
	os.makedirs(outdir)
#! draw_cmd
draw_cmd="vz_e>>hdcmd(100,-12,6)"
#! Create h[dtyp][sctr]
h=[[0 for j in range(NSCTR)]for i in range(NDTYP)]
print "h=",h
#sys.exit()
#! fill h
for idtyp in range(NDTYP):
	for isctr in range(NSCTR):
		cut=ROOT.TCut("sector==%d"%(isctr+1))
		print "TTree::Draw() for",DTYP_NAME[idtyp],"sector=",(isctr+1)
		#print T[idtyp][icorr].GetName()
		T[idtyp].Draw(draw_cmd,cut)#[idtyp][isctr][icorr],cut)
		#! Store histogram
		#ROOT.gStyle.SetOptStat("n")#! This line is critically needed here!
		htmp=ROOT.gDirectory.Get("hdcmd")#draw_cmd_hst_name[idtyp][isctr][icorr])
		h[idtyp][isctr]=htmp.Clone()
		h[idtyp][isctr].SetName("h_%s_s%d"%(DTYP_NAME[idtyp],isctr+1))
		h[idtyp][isctr].SetTitle("%s"%cut.GetTitle())
		#print h[idtyp][isctr][icorr].GetName()

#! Draw and save plots
ROOT.gStyle.SetOptStat("n")
ROOT.gStyle.SetOptFit(1)#111)
cname="czvtx"
czvtx=ROOT.TCanvas(cname,cname,2000,1500)
czvtx.Divide(3,2)
l=ROOT.TLegend(0.7,0.8,0.9,0.9)#,"","NDC");
for isctr in range(NSCTR):
	#! hist aesthetics
	h[PTGT][isctr].SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
	h[ETGT][isctr].SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
	h[PTGT][isctr].SetLineWidth(2)
	h[ETGT][isctr].SetLineWidth(2)
	#! [01-22-17] Remove main title and stats box
	for htmp in [h[PTGT][isctr],h[ETGT][isctr]]:
		htmp.SetTitle("")
		htmp.SetStats(0)
	#! Draw PTGT
	pad=czvtx.cd(isctr+1)
	h[PTGT][isctr].Draw("HIST")
	#! Axes title
	h[PTGT][isctr].SetXTitle("z-vertex position [cm]")
	pad.SetLeftMargin(0.20)
	h[PTGT][isctr].GetYaxis().SetTitleOffset(2.0)
	h[PTGT][isctr].SetYTitle("N_{entries}")
	#! Now draw etgt after scaling to Qrat (=Qptgt/Qetgt)
	Qrat=10
	h[ETGT][isctr].Sumw2()
	h[ETGT][isctr].Scale(Qrat)
	h[ETGT][isctr].Draw("HIST sames")
	# #! Draw legend 
	if isctr==0:
		l.AddEntry(h[PTGT][isctr],'ptgt-zvtx')
		l.AddEntry(h[ETGT][isctr],'etgt-zvtx')
		l.Draw()
czvtx.Update()
czvtx.SaveAs("%s/%s.png"%(outdir,cname)) 
czvtx.SaveAs("%s/%s.pdf"%(outdir,cname)) 
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
	
