#!/usr/bin/python
from __future__ import division

import os,sys,glob
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict

import numpy as np

OUTDIR=os.path.join("%s/study_W"%(os.environ['STUDY_ELASTIC_DATADIR']))
if not os.path.exists(OUTDIR):
	os.makedirs(OUTDIR)

NSCTR=6

NDTYP=2
ER,SR=range(NDTYP)
DTYP_NAME=['ER','SR']

#! Create FIN[DTYP],T[DTYP]
FIN=[[] for i in range(NDTYP)]
T=[[] for i in range(NDTYP)]

#! Obtain FIN[DTYP],T[DTYP]
FIN[ER]=ROOT.TFile("%s/study_elast_xsec_081715/delastR_tree.root"%os.environ['DELASTDIR_EXP'])
FIN[SR]=ROOT.TFile("%s/study_elast_xsec_081715/recon/delastR_tree.root"%os.environ['DELASTDIR_SIM'])
T[ER]=FIN[ER].Get("delast/monitor/t")
T[SR]=FIN[SR].Get("delast/monitor/t")

#FIN[ER]=ROOT.TFile(("%s/defid.root"%os.environ['STUDY_EFID_2PI_DATADIR_EXP']))
#FIN[SR]=ROOT.TFile(("%s/defid.root"%os.environ['STUDY_EFID_2PI_DATADIR_SIM']))	
#T[ER]=FIN[ER].Get("d2piR/tR")
#T[SR]=FIN[SR].Get("d2piR/tR")

DRAW_CMD_HST_NAME=[[] for i in range(NDTYP)]
DRAW_CMD=[[] for i in range(NDTYP)]
for idtyp in range(NDTYP):
	DRAW_CMD_HST_NAME[idtyp]="h%s"%DTYP_NAME[idtyp]	
	#DRAW_CMD[idtyp]="W>>%s(100,1.3,3.0)"%DRAW_CMD_HST_NAME[idtyp]
	DRAW_CMD[idtyp]="W>>%s(100,0.7,1.2)"%DRAW_CMD_HST_NAME[idtyp]

THETA_MIN=np.arange(14,50,4)
THETA_MAX=np.arange(18,54,4)
#H=[[[]for i in range(len(THETA_MIN))],[[]for i in range(len(THETA_MIN))]]
#C=[[]for i in range(len(THETA_MIN))]
#! H[NDTYP][NSCTR][THETA]
H=[[[[]for k in range(len(THETA_MIN))]for j in range(NSCTR)]for i in range(NDTYP)]
#! C[NSCTR][THETA]
C=[[[]for k in range(len(THETA_MIN))]for j in range(NSCTR)]
#! Make plots
c=ROOT.TCanvas("c","c") #! default canvas
for idtyp in range(NDTYP):
	for isctr in range(NSCTR):
		for itheta in range(len(THETA_MIN)):
			#if idtyp==ER: continue
			#T[idtyp].Draw("W>>h%s(100,0.7,1.2)"%DTYP_NAME[idtyp])
			cut=ROOT.TCut("sector==%d && theta_e>=%d && theta_e<=%d"%(isctr+1,THETA_MIN[itheta],THETA_MAX[itheta]));
			print("%s"%DRAW_CMD[idtyp])
			T[idtyp].Draw(DRAW_CMD[idtyp],cut)
	
			#! Store histogram
			ROOT.gStyle.SetOptStat("ne")
			htmp=ROOT.gDirectory.Get(DRAW_CMD_HST_NAME[idtyp])
			H[idtyp][isctr][itheta]=htmp.Clone()
			H[idtyp][isctr][itheta].SetTitle("%s"%cut.GetTitle())

#! Draw and save plots
for isctr in range(NSCTR):
	outdir="%s/sector%d"%(OUTDIR,isctr+1)
	if not os.path.exists(outdir):
        	os.makedirs(outdir)
	for itheta in range(len(THETA_MIN)):
		cname="%d_%d"%(THETA_MIN[itheta],THETA_MAX[itheta])
		C[isctr][itheta]=ROOT.TCanvas(cname,cname)
		s=H[ER][isctr][itheta].GetMaximum()/H[SR][isctr][itheta].GetMaximum()
		H[SR][isctr][itheta].Scale(s)
		H[ER][isctr][itheta].SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
		H[SR][isctr][itheta].SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
		H[ER][isctr][itheta].Draw()
		H[SR][isctr][itheta].Draw("same")
		C[isctr][itheta].SaveAs("%s/%s.png"%(outdir,cname))

#H[ER].SetLineColor(kBlue)
#H[SR].SetLineColor(kRed)
#s=H[ER].GetMaximum()/H[SR].GetMaximum()
#H[SR].Scale(s)
#H[ER][0].Draw()
#H[SR][0].Draw("same")

#! If wanting to keep TCanvas open till program exits				
#if not ROOT.gROOT.IsBatch():
#	plt.show()
#	# wait for you to close the ROOT canvas before exiting
#	wait(True)
