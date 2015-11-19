#!/usr/bin/python
from __future__ import division

import os,sys
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict

import math

import numpy as np

Q2MIN=1.25
Q2MAX=5.25
Q2BINW=0.5
WMIN=1.400
WMAX=2.125
WBINW=0.025

#Q2BIN_LEL=np.arange(Q2MIN,2.25,Q2BINW)
#WBIN_LEL=np.arange(WMIN,1.450,WBINW)
Q2BIN_LEL=np.arange(Q2MIN,Q2MAX,Q2BINW)
WBIN_LEL=np.arange(WMIN,WMAX,WBINW)
print Q2BIN_LEL
print WBIN_LEL

NDTYP=2
ER,SR=range(NDTYP)
DTYP_NAME=["ER","SR"]
FIN=[[] for i in range(NDTYP)]
T=[[] for i in range(NDTYP)]

NTOP=4
PPIPPIM,PPIP,PPIM,PIPPIM=range(4)
TOP_NAME=['p_pip_pim','p_pip','p_pim','pip_pim']

DRAW_CMD_MM =[[]for i in range(NTOP)]
DRAW_CMD_MM2=[[]for i in range(NTOP)]
for itop in range(NTOP):
	top=itop+1
	if top==1:
		cmd_mm="mmppippim>>hmm(20000, -0.50,  0.4)"
		cmd_mm2="mm2ppippim>>hmm2(600,-0.1,0.1)"	
	elif top==2:
		cmd_mm="mmppip>>hmm(200,-0.5,1)"
                cmd_mm2="mm2ppip>>hmm2(200,-0.25,1)"
	elif top==3:
		cmd_mm="mmppim>>hmmppim(200,-0.5,1)"
                cmd_mm2="mm2ppim>>hmm2ppim(200,-0.25,1)"
	elif top==4:
		cmd_mm="mmpippim>>hmm(100,0.00,2.0)"
                cmd_mm2="mm2pippim>>hmm2(100,0.00,2.0)"
	
	DRAW_CMD_MM[itop] =cmd_mm
	DRAW_CMD_MM2[itop]=cmd_mm2

#! MM cuts lines
MM2_T2_CUTL_EI=ROOT.TLine(-0.04,0,-0.04,0);
MM2_T2_CUTH_EI=ROOT.TLine(0.06,0,0.06,0);
MM_T2_CUTL_EI=ROOT.TLine(-math.sqrt(0.04),0,-math.sqrt(0.04),0);
MM_T2_CUTH_EI=ROOT.TLine(math.sqrt(0.06),0,math.sqrt(0.06),0);
MM2_T2_CUTL_EI.SetLineColor(ROOT.gROOT.ProcessLine("kGreen"))
MM2_T2_CUTH_EI.SetLineColor(ROOT.gROOT.ProcessLine("kGreen"))
MM_T2_CUTL_EI.SetLineColor(ROOT.gROOT.ProcessLine("kGreen"))
MM_T2_CUTH_EI.SetLineColor(ROOT.gROOT.ProcessLine("kGreen"))
MM2_T2_CUTL_EI.SetLineWidth(3)
MM2_T2_CUTH_EI.SetLineWidth(3)
MM_T2_CUTL_EI.SetLineWidth(3)
MM_T2_CUTH_EI.SetLineWidth(3)

MM2_T2_CUTL_AT=ROOT.TLine(0.00,0,0.00,0);
MM2_T2_CUTH_AT=ROOT.TLine(0.04,0,0.04,0);
MM_T2_CUTL_AT=ROOT.TLine(-math.sqrt(0.00),0,-math.sqrt(0.00),0);
MM_T2_CUTH_AT=ROOT.TLine(math.sqrt(0.04),0,math.sqrt(0.04),0);
MM2_T2_CUTL_AT.SetLineColor(ROOT.gROOT.ProcessLine("kYellow"))
MM2_T2_CUTH_AT.SetLineColor(ROOT.gROOT.ProcessLine("kYellow"))
MM_T2_CUTL_AT.SetLineColor(ROOT.gROOT.ProcessLine("kYellow"))
MM_T2_CUTH_AT.SetLineColor(ROOT.gROOT.ProcessLine("kYellow"))
MM2_T2_CUTL_AT.SetLineWidth(3)
MM2_T2_CUTH_AT.SetLineWidth(3)
MM_T2_CUTL_AT.SetLineWidth(3)
MM_T2_CUTH_AT.SetLineWidth(3)

FIN[ER]=ROOT.TFile("%s/d2pi_exp/d2piR_hvy.root"%os.environ['STUDY_BG_DATADIR'])
FIN[SR]=ROOT.TFile("%s/d2pi_sim/d2piR_hvy.root"%os.environ['STUDY_BG_DATADIR'])
#FIN[ER]=ROOT.TFile("%s/d2pi_exp/d2piR_hvy.root.small"%os.environ['STUDY_BG_DATADIR'])
#FIN[SR]=ROOT.TFile("%s/d2pi_sim/d2piR_hvy.root.small"%os.environ['STUDY_BG_DATADIR'])
T[ER]=FIN[ER].Get("d2piR/tR")
T[SR]=FIN[SR].Get("d2piR/tR")
#print FIN[ER].GetName()
#print T[ER].GetName()

OUTDIR="%s/results"%os.environ['STUDY_BG_DATADIR']
if not os.path.exists(OUTDIR):
	os.makedirs(OUTDIR)
#! Plotting aesthetics
CLRL=[ROOT.gROOT.ProcessLine("kBlue"),ROOT.gROOT.ProcessLine("kRed")]

#! Obtain HMM[DTYP][Q2][W][NTOP], HMM2[DTYP][Q2][W][NTOP]
HMM =[[[[[]for l in range(NTOP)]for k in range(len(WBIN_LEL))]for j in range(len(Q2BIN_LEL))] for i in range(NDTYP)]
HMM2=[[[[[]for l in range(NTOP)]for k in range(len(WBIN_LEL))]for j in range(len(Q2BIN_LEL))] for i in range(NDTYP)]
for idtyp in range(NDTYP):
	for i,q2bin_le in enumerate(Q2BIN_LEL):
		for j,wbin_le in enumerate(WBIN_LEL):
			q2min=q2bin_le
			q2max=q2bin_le+Q2BINW
			wmin=wbin_le
                        wmax=wbin_le+WBINW
			cut_q2w=ROOT.TCut("Q2>%f && Q2<%f & &W>%f && W<%f"%(q2min,q2max,wmin,wmax))
			for itop in range(NTOP):
				top=itop+1
				cut_top=ROOT.TCut("top==%d"%top)
				cut=ROOT.TCut(cut_q2w)
				cut+=cut_top
				T[idtyp].Draw(DRAW_CMD_MM[itop],cut)
				T[idtyp].Draw(DRAW_CMD_MM2[itop],cut)
				#! Store histogram
				ROOT.gStyle.SetOptStat("ne")
				#! hmm
				htmp=ROOT.gDirectory.Get("hmm")
				HMM[idtyp][i][j][itop]=htmp.Clone()
				#! hmm2
				htmp=ROOT.gDirectory.Get("hmm2")
                                HMM2[idtyp][i][j][itop]=htmp.Clone()
				#! Set Namei,Title and other plotting aesthetics
				HMM[idtyp][i][j][itop].SetName("%s"%DTYP_NAME[idtyp])
				HMM2[idtyp][i][j][itop].SetName("%s"%DTYP_NAME[idtyp])
				HMM[idtyp][i][j][itop].SetTitle("MM_{%s} %.2f-%.2f_%.3f-%.3f"%(TOP_NAME[itop],q2min,q2max,wmin,wmax))
				HMM2[idtyp][i][j][itop].SetTitle("MM2_{%s} %.2f-%.2f_%.3f-%.3f"%(TOP_NAME[itop],q2min,q2max,wmin,wmax))
				HMM[idtyp][i][j][itop].SetMarkerColor(CLRL[idtyp])
				HMM[idtyp][i][j][itop].SetLineColor(CLRL[idtyp])
				HMM2[idtyp][i][j][itop].SetMarkerColor(CLRL[idtyp])
                                HMM2[idtyp][i][j][itop].SetLineColor(CLRL[idtyp])

#! Plot HMM/MM2[DTYP][Q2][W][NTOP]
for itop in range(NTOP):
	for i,q2bin_le in enumerate(Q2BIN_LEL):
		for j,wbin_le in enumerate(WBIN_LEL):
			q2min=q2bin_le
                	q2max=q2bin_le+Q2BINW
                	wmin=wbin_le
                	wmax=wbin_le+WBINW

			outdir_mm="%s/MM_%s/%.2f-%.2f"%(OUTDIR,TOP_NAME[itop],q2min,q2max)
			if not os.path.exists(outdir_mm):
        			os.makedirs(outdir_mm)
			outdir_mm2="%s/MM2_%s/%.2f-%.2f"%(OUTDIR,TOP_NAME[itop],q2min,q2max)
                        if not os.path.exists(outdir_mm2):
                                os.makedirs(outdir_mm2)

			cname="c_%.3f-%.3f"%(wmin,wmax)
			#! mm
			cmm=ROOT.TCanvas(cname,cname)
			#! Scale hists for display
			HMM[ER][i][j][itop].Sumw2()
			HMM[SR][i][j][itop].Sumw2()
			max_ER=HMM[ER][i][j][itop].GetMaximum()
			if max_ER==0:max_ER=1
			max_SR=HMM[SR][i][j][itop].GetMaximum()
			if max_SR==0:max_SR=1
			scl_fctr_SR=max_ER/max_SR
			if scl_fctr_SR==0:scl_fctr_SR=1
			HMM[SR][i][j][itop].Scale(scl_fctr_SR)
			#! Draw
			HMM[ER][i][j][itop].Draw()
			HMM[SR][i][j][itop].Draw("same")
			if (itop+1==2):
				MM_T2_CUTL_EI.SetY1(HMM[ER][i][j][itop].GetMinimum())
                                MM_T2_CUTL_EI.SetY2(HMM[ER][i][j][itop].GetMaximum())
                                MM_T2_CUTH_EI.SetY1(HMM[ER][i][j][itop].GetMinimum())
                                MM_T2_CUTH_EI.SetY2(HMM[ER][i][j][itop].GetMaximum())
				MM_T2_CUTL_AT.SetY1(HMM[ER][i][j][itop].GetMinimum())
				MM_T2_CUTL_AT.SetY2(HMM[ER][i][j][itop].GetMaximum())
				MM_T2_CUTH_AT.SetY1(HMM[ER][i][j][itop].GetMinimum())
                                MM_T2_CUTH_AT.SetY2(HMM[ER][i][j][itop].GetMaximum())
				MM_T2_CUTL_EI.Draw("same")
				MM_T2_CUTH_EI.Draw("same")				
				MM_T2_CUTL_AT.Draw("same")
				MM_T2_CUTH_AT.Draw("same")
				
			#! Save canvas
			cmm.SaveAs("%s/%s.png"%(outdir_mm,cname))

			#! mm2
                        cmm2=ROOT.TCanvas(cname,cname)
                        #! Scale hists for display
                        HMM2[ER][i][j][itop].Sumw2()
                        HMM2[SR][i][j][itop].Sumw2()
                        max_ER=HMM2[ER][i][j][itop].GetMaximum()
                        if max_ER==0:max_ER=1
                        max_SR=HMM2[SR][i][j][itop].GetMaximum()
                        if max_SR==0:max_SR=1
                        scl_fctr_SR=max_ER/max_SR
                        if scl_fctr_SR==0:scl_fctr_SR=1
                        HMM2[SR][i][j][itop].Scale(scl_fctr_SR)
                        #! Draw
                        HMM2[ER][i][j][itop].Draw()
                        HMM2[SR][i][j][itop].Draw("same")
			if (itop+1==2):
                                MM2_T2_CUTL_EI.SetY1(HMM2[ER][i][j][itop].GetMinimum())
                                MM2_T2_CUTL_EI.SetY2(HMM2[ER][i][j][itop].GetMaximum())
                                MM2_T2_CUTH_EI.SetY1(HMM2[ER][i][j][itop].GetMinimum())
                                MM2_T2_CUTH_EI.SetY2(HMM2[ER][i][j][itop].GetMaximum())
                                MM2_T2_CUTL_AT.SetY1(HMM2[ER][i][j][itop].GetMinimum())
                                MM2_T2_CUTL_AT.SetY2(HMM2[ER][i][j][itop].GetMaximum())
                                MM2_T2_CUTH_AT.SetY1(HMM2[ER][i][j][itop].GetMinimum())
                                MM2_T2_CUTH_AT.SetY2(HMM2[ER][i][j][itop].GetMaximum())
                                MM2_T2_CUTL_EI.Draw("same")
                                MM2_T2_CUTH_EI.Draw("same")
                                MM2_T2_CUTL_AT.Draw("same")
                                MM2_T2_CUTH_AT.Draw("same")
                        #! Save canvas
                        cmm2.SaveAs("%s/%s.png"%(outdir_mm2,cname))
		
			
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
	
