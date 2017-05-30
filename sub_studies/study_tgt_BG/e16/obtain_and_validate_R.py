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
NDTYP=3
PTGT_LSE,PTGT_TGT,ETGT=range(NDTYP)
DTYP_NAME=["ptgt_lse","ptgt_tgt","etgt"]

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
	if       DTYP_NAME[idtyp]=="ptgt_lse": FIN[idtyp]=ROOT.TFile("%s/dtgt_ptgt_cutruns_lse.root"%DATADIR_INPUT)
	elif     DTYP_NAME[idtyp]=="ptgt_tgt": FIN[idtyp]=ROOT.TFile("%s/dtgt_ptgt_cutruns_tgt.root"%DATADIR_INPUT)
	elif     DTYP_NAME[idtyp]=="etgt": FIN[idtyp]=ROOT.TFile("%s/dtgt_etgt_cutruns.root"%DATADIR_INPUT)
	
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
outdir=os.path.join(DATADIR_OUTPUT,"results_R_%s"%DATE)
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

#! Draw and save plots:
#! 1. ETGT by itself
#! 2. Two canvases, one each for zvtx-PTGT_LSE and ztx-PTGT_TGT
#! + zvtx-ETGT is drawn scaled by R(LSE/TGT)=Q_PTGT(LSR/TGT)/Q_ETGT

#! Set up R(LSE/TGT)
#! Set Q(PTGT_LSE,PTGT_TGT,ETGT)
#! + Q for PTGT_LSE/TGT from: $STUDY_LUM_E16_DATADIR/results_052417/lum_results.txt
#! + Q for ETGT from:         $STUDY_LUM_E16_DATADIR/results_etgt_052617/lum_results.txt
Q=[0 for i in range(NDTYP)]
Q[PTGT_LSE]=21.33 #! mC
Q[PTGT_TGT]=20.90 #! mC
Q[ETGT]= 1.56 #! mC
R=[0 for i in range(NDTYP)] #! Note that the last index is not used
R[PTGT_LSE]=Q[PTGT_LSE]/Q[ETGT]
R[PTGT_TGT]=Q[PTGT_TGT]/Q[ETGT]
#! Write out these R factors to a file
foutR=open('%s/R.txt'%outdir,'w')
foutR.write('Q_etgt=%.2f mC\n'%Q[ETGT])
foutR.write('***  for cutruns_lse ***\n')
foutR.write("Q_ptgt_lse=%.2f mC\n"%Q[PTGT_LSE])
foutR.write("R_ptgt_lse=%.2f\n"%R[PTGT_LSE])
foutR.write('*** for cutruns_tgt ***\n')
foutR.write("Q_ptgt_tgt=%.2f mC\n"%Q[PTGT_TGT])
foutR.write("R_ptgt_tgt=%.2f\n"%R[PTGT_TGT])
foutR.close()

#! General ROOT plot aesthetics
ROOT.gStyle.SetOptStat("n")
ROOT.gStyle.SetOptFit(1)#111)


#! 1. ETGT by itself
cname="czvtx_etgt"
czvtx=ROOT.TCanvas(cname,cname,2000,1500)
czvtx.Divide(3,2)
for isctr in range(NSCTR):
	pad=czvtx.cd(isctr+1)
	#! hist aesthetic
	h[ETGT][isctr].SetLineWidth(2)
	#! Remove main title and stats box
	h[ETGT][isctr].SetTitle("")
	h[ETGT][isctr].SetStats(0)
	#! Axes title
	h[ETGT][isctr].SetXTitle("z-vertex position [cm]")
	pad.SetLeftMargin(0.20)
	h[ETGT][isctr].GetYaxis().SetTitleOffset(2.0)
	h[ETGT][isctr].SetYTitle("N_{entries}")
	#! Draw
	h[ETGT][isctr].Draw("HIST e")
czvtx.SaveAs("%s/%s.png"%(outdir,cname)) 
czvtx.SaveAs("%s/%s.pdf"%(outdir,cname))

#! 2. Two canvases, one each for zvtx-PTGT_LSE and ztx-PTGT_TGT
#! + zvtx-ETGT is drawn scaled by R(LSE/TGT)=Q_PTGT(LSR/TGT)/Q_ETGT
for ptgt in [PTGT_LSE,PTGT_TGT]:
	cname="czvtx_rto_%s_etgt"%DTYP_NAME[ptgt]
	czvtx=ROOT.TCanvas(cname,cname,2000,1500)
	czvtx.Divide(3,2)
	l=ROOT.TLegend(0.7,0.8,0.9,0.9)#,"","NDC");
	#! Data structure for keeping cloned copies for h[ETGT] (since it is scaled differently for each PTGT)
	hetgt=[0 for i in range(NSCTR)]
	for isctr in range(NSCTR):
		#! hist aesthetics
		h[ptgt][isctr].SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
		h[ETGT][isctr].SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
		h[ptgt][isctr].SetLineWidth(2)
		h[ETGT][isctr].SetLineWidth(2)
		#! [01-22-17] Remove main title and stats box
		for htmp in [h[ptgt][isctr],h[ETGT][isctr]]:
			htmp.SetTitle("")
			htmp.SetStats(0)
		#! Draw ptgt
		pad=czvtx.cd(isctr+1)
		h[ptgt][isctr].Draw("HIST e")
		#! Axes title
		h[ptgt][isctr].SetXTitle("z-vertex position [cm]")
		pad.SetLeftMargin(0.20)
		h[ptgt][isctr].GetYaxis().SetTitleOffset(2.0)
		h[ptgt][isctr].SetYTitle("N_{entries}")
		#! Now draw etgt after scaling to Qrat (=Qptgt/Qetgt)
		#! Firt get clone of h[ETGT]
		hetgt[isctr]=h[ETGT][isctr].Clone("hetgt")
		print "R(%s)=%.2f"%(ptgt,R[ptgt])
		hetgt[isctr].Sumw2()
		hetgt[isctr].Scale(R[ptgt])
		hetgt[isctr].Draw("HIST e sames")
		# #! Draw legend 
		if isctr==0:
			l.AddEntry(h[ptgt][isctr],'%s-zvtx'%DTYP_NAME[ptgt])
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
	
