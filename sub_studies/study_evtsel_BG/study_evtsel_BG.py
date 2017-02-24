#!/usr/bin/python
from __future__ import division

import os,sys,datetime
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict
import itertools

import numpy as np

import math

'''
+ Usage: study_evtsel.py debug=False
'''

USAGE='study_evtsel.py debug=False'

#! Get input data from user
DBG=False
if len(sys.argv)>1: #i.e. debug entered by user
	if    sys.argv[1]=="True":  DBG=True
	elif  sys.argv[1]=="False": DBG=False
	else: sys.exit('DBG=%s is not valid. usage: %s'%(DBG,USAGE))

print "DBG=",DBG
#sys.exit()

#! *** Prepare input data structure: FIN/T[dtyp] ***
NDTYP=3
ER,SR,SR3PI=range(NDTYP)
DTYP_NAME=["ER","SR","SR3PI"]

FIN=[0 for i in range(NDTYP)]
T  =[0 for i in range(NDTYP)]
#print "FIN=",FIN
#print "T=",T
#! ***

#! *** Prepare OUTDIR***
DATE=datetime.datetime.now().strftime('%m%d%y')
OUTDIR=os.path.join(os.environ['STUDY_EVTSEL_BG_DATADIR'],"evtsel_BG_%s"%DATE)
if DBG==True: OUTDIR=os.path.join(os.environ['STUDY_EVTSEL_BG_DATADIR'],"evtsel_BG_dbg_%s"%DATE)
print "OUTDIR=",OUTDIR
if not os.path.exists(OUTDIR):
	os.makedirs(OUTDIR)
#sys.exit()
#! ***

#! *** Prepare NENTRIES ***
NENTRIES=1000000000
#if DBG==True: NENTRIES=1000

#! *** Now prepare INDIR and input data structures: FIN/T[dtyp] ***
INDIR=[[] for i in range(NDTYP)]
INDIR[ER]   =os.path.join(os.environ['D2PIDIR_EXP_E16'],"data_pcorr_011116")#! debug /bk
INDIR[SR]   =os.path.join(os.environ['D2PIDIR_SIM_E16'],"data_pcorr_011116")#! debug /bl
INDIR[SR3PI]=os.path.join(os.environ['D2PIDIR_SIM_E16'],"data_evtsel_BG_020117")

print "INDIR[ER]=",INDIR[ER]
print "INDIR[SR]=",INDIR[SR]
print "INDIR[SR3PI]=",INDIR[SR3PI]
#sys.exit()

FIN[ER]   =ROOT.TFile("%s/d2piR_npcorr.root"%INDIR[ER])
FIN[SR]   =ROOT.TFile("%s/d2piR_npcorr.root"%INDIR[SR])
FIN[SR3PI]=ROOT.TFile("%s/d3piPS.root"%INDIR[SR3PI])
for idtyp in range(NDTYP):
	T[idtyp]=FIN[idtyp].Get("d2piR/tR")
#sys.exit()

#! *** Now make plots ***

#! sector information
NSCTR=6

#! *** Now make MM plots ***
WMIN=1.400
WMAX=2.125
WBINW=0.025
Q2BIN_LEL=[2.00,2.40,3.00,3.50,4.20,2.00]#[1.50,2.00,2.40,3.00,3.50,4.20] #debug
Q2BIN_UEL=[2.40,3.00,3.50,4.20,5.00,5.00]#[2.00,2.40,3.00,3.50,4.20,5.00] #debug
# Q2BIN_LEL=[2.00]
# Q2BIN_UEL=[5.00]13,acfjlmpqstvw-
if DBG:
	Q2BIN_LEL=[2.00]
	Q2BIN_UEL=[5.00]
WBIN_LEL=np.arange(WMIN,WMAX,WBINW)
#print Q2BIN_LEL
#print WBIN_LEL

#! MM cuts lines: EI and AT
#! EI
EI_MM2L,EI_MM2H=-0.04,0.06
EI_MML, EI_MMH=-math.sqrt(math.fabs(EI_MM2L)),math.sqrt(EI_MM2H)
MM2_T2_CUTL_EI=ROOT.TLine(EI_MM2L,0,EI_MM2L,0);
MM2_T2_CUTH_EI=ROOT.TLine(EI_MM2H,0,EI_MM2H,0);
MM_T2_CUTL_EI=ROOT.TLine(EI_MML,0,EI_MML,0);
MM_T2_CUTH_EI=ROOT.TLine(EI_MMH,0,EI_MMH,0);
for l in [MM2_T2_CUTL_EI,MM2_T2_CUTH_EI,MM_T2_CUTL_EI,MM_T2_CUTH_EI]:
	l.SetLineColor(ROOT.gROOT.ProcessLine("kMagenta"))
	l.SetLineStyle(2)
	l.SetLineWidth(3)
#! AT
AT_MM2L,AT_MM2H=-0.00,0.04
AT_MML, AT_MMH=-math.sqrt(math.fabs(AT_MM2L)),math.sqrt(AT_MM2H)
MM2_T2_CUTL_AT=ROOT.TLine(AT_MM2L,0,AT_MM2L,0);
MM2_T2_CUTH_AT=ROOT.TLine(AT_MM2H,0,AT_MM2H,0);
MM_T2_CUTL_AT=ROOT.TLine(AT_MML,0,AT_MML,0);
MM_T2_CUTH_AT=ROOT.TLine(AT_MMH,0,AT_MMH,0);
for l in [MM2_T2_CUTL_AT,MM2_T2_CUTH_AT,MM_T2_CUTL_AT,MM_T2_CUTH_AT]:
	l.SetLineColor(ROOT.gROOT.ProcessLine("kOrange+8"))
	l.SetLineStyle(2)
	l.SetLineWidth(3)

#! TLine depicting pion mass debug
MASS_PION=0.13957018;
MASS_PION_LINE=ROOT.TLine(MASS_PION,0,MASS_PION,0)
MASS_PION_LINE.SetLineColor(ROOT.gROOT.ProcessLine("kGreen"))
MASS_PION_LINE.SetLineWidth(3)

#! DRAW_CMD
NPLTS=2
MM,MM2=range(NPLTS)
PLT_NAME=['MM','MM2']
PLT_TITLE=['MM','MM^{2}']
PLT_UNIT=['GeV','GeV^{2}']

DRAW_CMD=[0 for i in range(NPLTS)]
DRAW_CMD[MM] ="mmppip>>hmmcmd(200,-0.5,1)"
DRAW_CMD[MM2]="mm2ppip>>hmmcmd(100,-0.2,0.2)"

#! Create hmm[dtyp][q2][w][plt]
hmm=[[[[0 for l in range(NPLTS)] for k in range(len(WBIN_LEL))] for j in range(len(Q2BIN_LEL))] for i in range(NDTYP)]
#print "hmm=",hmm
#sys.exit()

#! Create domain for making [ER,SR,SR3PI]*[MM,MM2] plots within a q2-w bin
d=[[ER,SR,SR3PI],[MM,MM2]]
D=list(itertools.product(*d))

#! Now fill hmm
for iq,q2bin_le in enumerate(Q2BIN_LEL):
	if DBG==True and iq>0: continue #! debug
	q2min=q2bin_le
	q2max=Q2BIN_UEL[iq]
	for iw,wbin_le in enumerate(WBIN_LEL):
		if DBG==True and (iw!=0 and iw!=15 and iw!=28): continue #! debug
		wmin=wbin_le
		wmax=wbin_le+WBINW
		cut_q2w=ROOT.TCut("Q2>%f && Q2<%f && W>%f && W<%f"%(q2min,q2max,wmin,wmax))
		print "q2w=",cut_q2w.GetTitle()
		for r in D:
			idtyp,iplt=r[0],r[1]
			
			cut_top=ROOT.TCut("top==2")
			cut=ROOT.TCut(cut_q2w)
			cut+=cut_top
			print "TTree::Draw() for",DTYP_NAME[idtyp],PLT_NAME[iplt]
			#print T[idtyp].GetName()
			T[idtyp].Draw(DRAW_CMD[iplt],cut,"",NENTRIES)#[idtyp][isctr][icorr],cut)

			#! Store histogram
			#ROOT.gStyle.SetOptStat("ne")
			htmp=ROOT.gDirectory.Get("hmmcmd")#(draw_cmd_hst_name[idtyp][isctr][icorr])
			hmm[idtyp][iq][iw][iplt]=htmp.Clone()#"hmm_%s_%s"%(DTYP_NAME[idtyp],CORR_NAME[icorr]))
			hmm[idtyp][iq][iw][iplt].SetName("h%s_%s"%(PLT_NAME[iplt],DTYP_NAME[idtyp]))
			title="%s for Q^{2}=[%.2f,%.2f) GeV^{2}, W=[%.3f,%.3f) GeV"%(PLT_TITLE[iplt],q2min,q2max,wmin,wmax)
			hmm[idtyp][iq][iw][iplt].SetTitle(title)#"%s %.2f-%.2f_%.3f-%.3f"%(PLT_TITLE[iplt],q2min,q2max,wmin,wmax))#("%s"%cut.GetTitle())
			#print hmm[idtyp][i][j][icorr].GetName()

#! Plot and save hmm
ROOT.gStyle.SetOptStat(0) #"n"
ROOT.gStyle.SetOptFit(1)#111)
#! Aesthetic related objects
CLRS_DTYP=[ROOT.gROOT.ProcessLine("kBlue"),ROOT.gROOT.ProcessLine("kRed"),ROOT.gROOT.ProcessLine("kGreen")]

# c=ROOT.TCanvas()
# hER=hmm[ER][0][0][MM]
# hSR=hmm[SR][0][0][MM]
# hSR.SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
# hSR3PI=hmm[SR3PI][0][0][MM]
# hER.Draw()
# hSR.Draw("same")
# print hER.GetEntries(),hSR.GetEntries()

for iq,q2bin_le in enumerate(Q2BIN_LEL):
	if DBG==True and iq>0: continue #! debug
	q2min=q2bin_le
	q2max=Q2BIN_UEL[iq]
	outdir_q2w="%s/%.2f-%.2f"%(OUTDIR,q2min,q2max)
	if not os.path.exists(outdir_q2w):
		os.makedirs(outdir_q2w)
	for iw,wbin_le in enumerate(WBIN_LEL):
		if DBG==True and (iw!=0 and iw!=15 and iw!=28): continue #! debug
		wmin=wbin_le
		wmax=wbin_le+WBINW
		
		print "Going to plot hmm and hmm2 for %.2f-%.2f_%.3f-%.3f"%(q2min,q2max,wmin,wmax)
		#! Create cmm[plt]
		cmm=[0 for i in range(NPLTS)]
		#! Now plot
		for iplt in range(NPLTS):
			outdir_q2w_plt="%s/%s"%(outdir_q2w,PLT_NAME[iplt])
			if not os.path.exists(outdir_q2w_plt):
				os.makedirs(outdir_q2w_plt)
			cname="c%s_%.3f-%.3f"%(PLT_NAME[iplt],wmin,wmax)
			cmm[iplt]=ROOT.TCanvas(cname,cname)
			#l=ROOT.TLegend(0.1,0.3,0.3,0.4)#,"","NDC");
			#! hists aesthetics	
			for idtyp in range(NDTYP):
				hmm[idtyp][iq][iw][iplt].SetLineColor(CLRS_DTYP[idtyp])
				hmm[idtyp][iq][iw][iplt].SetMarkerColor(CLRS_DTYP[idtyp])
			#! Draw hists
			#! First directly get copy of hists so as to not use tedious indices
			hER=hmm[ER][iq][iw][iplt]
			hSR=hmm[SR][iq][iw][iplt]
			hSR3PI=hmm[SR3PI][iq][iw][iplt]
			#! First draw hER
			hER.SetXTitle("%s [%s]"%(PLT_TITLE[iplt], PLT_UNIT[iplt]))
			cmm[iplt].SetLeftMargin(0.20)
			hER.GetYaxis().SetTitleOffset(1.5)
			hER.SetYTitle("N_{entries}")
			hER.Draw()
			#! scale SR and then draw
			max_ER=hER.GetMaximum()
			if max_ER==0:max_ER=1
			max_SR=hSR.GetMaximum()
			if max_SR==0:max_SR=1
			scl_fctr_SR=max_ER/max_SR
			if scl_fctr_SR==0:scl_fctr_SR=1
			hSR.Scale(scl_fctr_SR)
			hSR.Draw("same")
			#! Now scale SR3PI to hER-hSR and draw
			#! First get hdiff=hER-hSR
			hdiff=hER.Clone("hdiff")
			hdiff.Add(hSR,-1)
			hdiff.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
			hdiff.SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
			hdiff.Draw("same")
			#! scale
			max_diff=hdiff.GetMaximum()
			if max_diff==0:max_diff=1
			max_SR3PI=hSR3PI.GetMaximum()
			if max_SR3PI==0:max_SR3PI=1
			scl_fctr_SR3PI=max_diff/max_SR3PI
			if scl_fctr_SR3PI==0:scl_fctr_SR3PI=1
			hSR3PI.Scale(scl_fctr_SR3PI)
			#! Draw
			hSR3PI.Draw("same")
			#! Create and add entries: hists and cut-lines
			#! legend
			l=ROOT.TLegend(0.66,0.72,0.9,0.9)#,"","NDC");	
			for hist,label in zip([hER,hSR,hdiff,hSR3PI],["exp","sim","exp-sim","sim-3pi-PS"]):
				l.AddEntry(hist,label,"lp")
			#! Draw cut lines and add them to legend
			if iplt==MM:
				#! Draw cut lines
				#! EI
				MM_T2_CUTL_EI.SetY1(hER.GetMinimum())
				MM_T2_CUTL_EI.SetY2(hER.GetMaximum())
				MM_T2_CUTH_EI.SetY1(hER.GetMinimum())
				MM_T2_CUTH_EI.SetY2(hER.GetMaximum())
				MM_T2_CUTL_EI.Draw("same")
				MM_T2_CUTH_EI.Draw("same")
				#! AT
				MM_T2_CUTL_AT.SetY1(hER.GetMinimum())
				MM_T2_CUTL_AT.SetY2(hER.GetMaximum())
				MM_T2_CUTH_AT.SetY1(hER.GetMinimum())
				MM_T2_CUTH_AT.SetY2(hER.GetMaximum())
				MM_T2_CUTL_AT.Draw("same")
				MM_T2_CUTH_AT.Draw("same")
				#! add to legend
				l.AddEntry(MM_T2_CUTL_EI,"%.2f GeV < MM < %.2f GeV"%(EI_MML,EI_MMH),"l")
				l.AddEntry(MM_T2_CUTL_AT,"%.2f GeV < MM < %.2f GeV"%(AT_MML,AT_MMH),"l")
			elif iplt==MM2:
				#! EI
				MM2_T2_CUTL_EI.SetY1(hER.GetMinimum())
				MM2_T2_CUTL_EI.SetY2(hER.GetMaximum())
				MM2_T2_CUTH_EI.SetY1(hER.GetMinimum())
				MM2_T2_CUTH_EI.SetY2(hER.GetMaximum())
				MM2_T2_CUTL_EI.Draw("same")
				MM2_T2_CUTH_EI.Draw("same")
				#! AT
				MM2_T2_CUTL_AT.SetY1(hER.GetMinimum())
				MM2_T2_CUTL_AT.SetY2(hER.GetMaximum())
				MM2_T2_CUTH_AT.SetY1(hER.GetMinimum())
				MM2_T2_CUTH_AT.SetY2(hER.GetMaximum())
				MM2_T2_CUTL_AT.Draw("same")
				MM2_T2_CUTH_AT.Draw("same")
				#! add to legend
				l.AddEntry(MM2_T2_CUTL_EI,"%.2f GeV^{2} < MM^{2} < %.2f GeV^{2}"%(EI_MM2L,EI_MM2H),"l")
				l.AddEntry(MM2_T2_CUTL_AT,"%.2f GeV^{2} < MM^{2} < %.2f GeV^{2}"%(AT_MM2L,AT_MM2H),"l")
			#! Draw legend
			l.Draw()

			cmm[iplt].SaveAs("%s/%s.png"%(outdir_q2w_plt,cname))
			cmm[iplt].SaveAs("%s/%s.pdf"%(outdir_q2w_plt,cname))
#! ***

# #! If wanting to keep TCanvas open till program exits				
# if not ROOT.gROOT.IsBatch():
# 	plt.show()
# 	# wait for you to close the ROOT canvas before exiting
# 	wait(True)

#if __name__ == "__main__":	
#	if len(sys.argv)==2:
#		plot_fid(sys.argv[1])
#	elif len(sys.argv)==3:
#		plot_fid(sys.argv[1],int(sys.argv[2])) 
	
