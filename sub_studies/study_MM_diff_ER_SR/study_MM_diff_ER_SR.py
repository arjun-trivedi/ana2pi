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

#atlib
import atlib as atlib

'''
[07-13-17] 
+ Estimate SE due to the fact that the MM distributions in ER and SR are not the same. 
+ Plots are generated for EI MM cuts and EI compared with with AT MM cuts in separate folders called EI and EI_AT, respectively.
SE = (1-(x/y))*100 where:
  + x,y = fraction of events under MM-cut for ER,SR, respectively
    + x,y are obtained using Fits to distributions
    + The Fit limits are optimized in W bins
  + This formula is derived using:
    + xsec  (x,y; x=y=1)       ~ ER/(SR/ST)
    + xsec' (x,y; x,y neq 1)   ~ xER/(ySR/ST)  
    --> SE = relative diff = (xsec-xsec'/xsec)*100 =  (1-(x/y))*100
  + If sign is negative it means that xsec'>xsec and therefore the fact that the cross-sections are overestimated
  + If sign is positive it means that xsec'<xsec and therefore the fact that the cross-sections are underestimated

+ Usage: study_MM_diff_ER_SR.py debug=False
'''

USAGE='study_MM_diff_ER_SR.py debug=False'

#! Get input data from user
DBG=False
if len(sys.argv)>1: #i.e. debug entered by user
	if    sys.argv[1]=="True":  DBG=True
	elif  sys.argv[1]=="False": DBG=False
	else: sys.exit('DBG=%s is not valid. usage: %s'%(DBG,USAGE))

print "DBG=",DBG
#sys.exit()

#! *** Prepare input data structure: FIN/T[dtyp] ***
NDTYP=2
ER,SR=range(NDTYP)
DTYP_NAME=["ER","SR"]

FIN=[0 for i in range(NDTYP)]
T  =[0 for i in range(NDTYP)]
#print "FIN=",FIN
#print "T=",T
#! ***

#! *** Prepare OUTDIR***
DATE=datetime.datetime.now().strftime('%m%d%y')
OUTDIR=os.path.join(os.environ['STUDY_MM_DIFF_ER_SR_DATADIR'],"results_%s"%DATE)
if DBG==True: OUTDIR=os.path.join(os.environ['STUDY_MM_DIFF_ER_SR_DATADIR'],"results_dbg_%s"%DATE)
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

print "INDIR[ER]=",INDIR[ER]
print "INDIR[SR]=",INDIR[SR]
#sys.exit()

FIN[ER]   =ROOT.TFile("%s/d2piR_wpcorr.root"%INDIR[ER])
FIN[SR]   =ROOT.TFile("%s/d2piR_npcorr.root"%INDIR[SR])
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
# Q2BIN_UEL=[5.00]
if DBG:
	Q2BIN_LEL=[2.00] #[3.50]#[2.00]
	Q2BIN_UEL=[5.00] #[4.20]#[5.00]
WBIN_LEL=np.arange(WMIN,WMAX,WBINW)
#! [05-10-17] Remove 2.125 if it occures in WBIN_LEL
#! + Sometimes because of using decimal numbers with np.arange, 2.125 is not omitted 
#!   and therefore it has to be done manually
del_val=False
idel_val=-9999
for i,x in enumerate(WBIN_LEL):
	if np.isclose(x,2.125): 
		del_val=True
		idel_val=i
#print len(WBIN_LEL)
if del_val: 
    print "In WBIN_LEL going to delete",WBIN_LEL[idel_val],"at index",idel_val
WBIN_LEL=np.delete(WBIN_LEL,[idel_val])
NWBINS=len(WBIN_LEL)
#print WBIN_LEL
#sys.exit()
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

#! MM Fit limits =f (dtyp,wbin)
print "*** Going to set up MM fit limits ***"
MM_FIT_LIMS  =[[0 for j in range(NWBINS)] for i in range(NDTYP)]
MM2_FIT_LIMS =[[0 for j in range(NWBINS)] for i in range(NDTYP)]
for iw in range(NWBINS):
	wbin_le=WBIN_LEL[iw]
	wmin=wbin_le
	wmax=wbin_le+WBINW
	print "wmin,wmax=",wmin,wmax
	#if np.isclose(wmin,1.400) and np.isclose(wmax,1.425):
	if (np.isclose(wmin,1.400) and np.isclose(wmax,1.425)):	
		print "1"
		#!mm
		mml_ER,mmh_ER=0.10,0.15
		mml_SR,mmh_SR=0.10,0.18
		#!mm2
		mm2l_ER,mm2h_ER=0.0025,0.0225#mml_ER**2,mmh_ER**2
		mm2l_SR,mm2h_SR=0.0025,0.03#mml_SR**2,mmh_SR**2
	elif ( (wmin>1.400 and wmax>1.425) and (wmin<1.500 and wmax<1.525) and (not np.isclose(wmin,1.500) and not np.isclose(wmax,1.525)) ):
		print "2"
		#!mm
		mml_ER,mmh_ER=0.10,0.15
		mml_SR,mmh_SR=0.10,0.18
		#!mm2
		mm2l_ER,mm2h_ER=0.0025,0.03#0.03 mml_ER**2,mmh_ER**2
		mm2l_SR,mm2h_SR=0.0025,0.04#0.04 mml_SR**2,mmh_SR**2
	elif (np.isclose(wmin,1.500) and np.isclose(wmax,1.525)):
		print "3"
		#!mm
		mml_ER,mmh_ER=0.10,0.17
		mml_SR,mmh_SR=0.10,0.20
		#!mm2
		mm2l_ER,mm2h_ER=0.0025,0.03#0.04 mml_ER**2,mmh_ER**2
		mm2l_SR,mm2h_SR=0.0025,0.04#0.05 mml_SR**2,mmh_SR**2
	elif ( (wmin>1.500 and wmax>1.525) and (wmin<1.775 and wmax<1.800) and (not np.isclose(wmin,1.775) and np.isclose(wmax,1.800)) ):
		print "4"
		#!mm
		mml_ER,mmh_ER=0.10,0.17
		mml_SR,mmh_SR=0.10,0.20
		#!mm2
		mm2l_ER,mm2h_ER=0.0025,0.03#0.04 mml_ER**2,mmh_ER**2
		mm2l_SR,mm2h_SR=0.0025,0.04#0.05 mml_SR**2,mmh_SR**2
	elif (np.isclose(wmin,1.775) and np.isclose(wmax,1.800)):
		print "5"
		#!mm
		mml_ER,mmh_ER=0.10,0.17
		mml_SR,mmh_SR=0.10,0.20
		#!mm2
		mm2l_ER,mm2h_ER=0.0025,0.03#0.04 mml_ER**2,mmh_ER**2
		mm2l_SR,mm2h_SR=0.0025,0.04#0.05 mml_SR**2,mmh_SR**2
	elif ( (wmin>1.775 and wmax>1.800) and (wmin<2.100 and wmax<2.125) and (not np.isclose(wmin,2.100) and np.isclose(wmax,2.125)) ):
		print "6"
		#!mm
		mml_ER,mmh_ER=0.10,0.17
		mml_SR,mmh_SR=0.10,0.20
		#!mm2
		mm2l_ER,mm2h_ER=0.0025,0.04#mml_ER**2,mmh_ER**2
		mm2l_SR,mm2h_SR=0.0025,0.05#mml_SR**2,mmh_SR**2
	else:
		print "7"
		#!mm
		mml_ER,mmh_ER=0.10,0.16
		mml_SR,mmh_SR=0.10,0.20
		#!mm2
		mm2l_ER,mm2h_ER=0.0025,0.03#mml_ER**2,mmh_ER**2
		mm2l_SR,mm2h_SR=0.0025,0.04#mml_SR**2,mmh_SR**2
	#!MM
	MM_FIT_LIMS[ER][iw]=[mml_ER,mmh_ER]#EI_MMH]# since 0-artifact in MM #[0.08,0.19]
	MM_FIT_LIMS[SR][iw]=[mml_SR,mmh_SR] #EI_MMH]# since 0-artifact in MM #[0.08,0.22]
	#!MM2
	MM2_FIT_LIMS[ER][iw]=[mm2l_ER,mm2h_ER]#[0.000,0.035]#[0.005,0.030]
	MM2_FIT_LIMS[SR][iw]=[mm2l_SR,mm2h_SR]#[0.000,0.035]#[0.005,0.035]
	for idtyp in range(NDTYP):
		print "MMl,MMh (%s)=%.3f,%.3f"%(DTYP_NAME[idtyp],MM_FIT_LIMS[idtyp][iw][0],MM_FIT_LIMS[idtyp][iw][1])
	for idtyp in range(NDTYP):
		print "MM2l,MM2h (%s)=%.3f,%.3f"%(DTYP_NAME[idtyp],MM2_FIT_LIMS[idtyp][iw][0],MM2_FIT_LIMS[idtyp][iw][1])
print "***"
#sys.exit()

#! DRAW_CMD
NPLTS=2
MM,MM2=range(NPLTS)
PLT_NAME=['MM','MM2']
PLT_TITLE=['MM','MM^{2}']
PLT_UNIT=['GeV','GeV^{2}']

DRAW_CMD=[0 for i in range(NPLTS)]
DRAW_CMD[MM] ="mmppip>>hmmcmd(300,-0.5,1)" #! 200
DRAW_CMD[MM2]="mm2ppip>>hmmcmd(200,-0.2,0.2)" #! 100

#! Create hmm[dtyp][q2][w][plt]
hmm=[[[[0 for l in range(NPLTS)] for k in range(len(WBIN_LEL))] for j in range(len(Q2BIN_LEL))] for i in range(NDTYP)]
#print "hmm=",hmm
#sys.exit()

#! Create domain for making [ER,SR,SR3PI]*[MM,MM2] plots within a q2-w bin
d=[[ER,SR],[MM,MM2]]
D=list(itertools.product(*d))

#! Now fill hmm
for iq,q2bin_le in enumerate(Q2BIN_LEL):
	if DBG==True and iq>0: continue #! debug
	q2min=q2bin_le
	q2max=Q2BIN_UEL[iq]
	for iw,wbin_le in enumerate(WBIN_LEL):
		#!DBGW
		if DBG==True and (iw!=0 and iw!=15 and iw!=28): continue #! debug #! (iw!=0 and iw!=15 and iw!=28)
		#if DBG==True and iw+1 > 10: continue
		#if DBG==True and (iw+1!=5 and iw+1!=6): continue
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

#! + Plot and save hmm
ROOT.gStyle.SetOptStat(0) #"n"
ROOT.gStyle.SetOptFit(0)#111)
#! Aesthetic related objects
CLRS_DTYP=[ROOT.gROOT.ProcessLine("kBlue"),ROOT.gROOT.ProcessLine("kRed")]

#! [07-13-17] Create and fill SE structure in this loop, where SE is calculated as:
#! SE = (1-(x/y))*100 where:
#!  + x,y = fraction of events under MM-cut for ER,SR, respectively
#!  + This formula is derived using:
#!    + xsec  (x,y; x=y=1)       ~ ER/(SR/ST)
#!    + xsec' (x,y; x,y neq 1)   ~ xER/(ySR/ST)  
#!    --> SE = relative diff = (xsec-xsec'/xsec)*100 =  (1-(x/y))*100
#!  + If sign is negative it means that xsec'>xsec and therefore the fact that the cross-sections are overestimated
#!  + If sign is positive it means that xsec'<xsec and therefore the fact that the cross-sections are underestimated

#! Create structure to hold estimating SE(%):SE[q2][w][plt][cut]=(bg,qbin,wbin)
NCUTS=2
EI,AT=range(NCUTS)
SE=[[[[0 for l in range(NCUTS)] for k in range(NPLTS)] for j in range(len(WBIN_LEL))] for i in range(len(Q2BIN_LEL))]


for iq,q2bin_le in enumerate(Q2BIN_LEL):
	if DBG==True and iq>0: continue #! debug
	q2min=q2bin_le
	q2max=Q2BIN_UEL[iq]
	outdir_q2w="%s/%.2f-%.2f"%(OUTDIR,q2min,q2max)
	if not os.path.exists(outdir_q2w):
		os.makedirs(outdir_q2w)
	for iw,wbin_le in enumerate(WBIN_LEL):
		#!DBGW
		if DBG==True and (iw!=0 and iw!=15 and iw!=28): continue #! debug #!(iw!=0 and iw!=15 and iw!=28)
		#if DBG==True and iw+1 > 10: continue
		#if DBG==True and (iw+1!=5 and iw+1!=6): continue
		wmin=wbin_le
		wmax=wbin_le+WBINW
				
		print "Going to plot hmm and hmm2 for %.2f-%.2f_%.3f-%.3f"%(q2min,q2max,wmin,wmax)
		#! Create plots purely for EI in one and comparison of AT-EI in another
		NCTYP=2
		CEI,CEIAT=range(NCTYP)
		#! Create cmm[plt][typ]
		cmm=[[0 for j in range(NCTYP)] for i in range(NPLTS)]
		#! Now plot
		for iplt in range(NPLTS):
			outdir_q2w_plt_EI="%s/EI/%s"%(outdir_q2w,PLT_NAME[iplt])
			outdir_q2w_plt_EIAT="%s/EI_AT/%s"%(outdir_q2w,PLT_NAME[iplt])
			for outdir in [outdir_q2w_plt_EI,outdir_q2w_plt_EIAT]:
				if not os.path.exists(outdir):
					os.makedirs(outdir)
			cname="c%s_%.3f-%.3f"%(PLT_NAME[iplt],wmin,wmax)

			#! First plot CEIAT
			cmm[iplt][CEIAT]=ROOT.TCanvas(cname,cname)
			#l=ROOT.TLegend(0.1,0.3,0.3,0.4)#,"","NDC");
			#! hists aesthetics	
			for idtyp in range(NDTYP):
				hmm[idtyp][iq][iw][iplt].SetLineColor(CLRS_DTYP[idtyp])
				hmm[idtyp][iq][iw][iplt].SetMarkerColor(CLRS_DTYP[idtyp])
			#! Draw hists
			#! First directly get copy of hists so as to not use tedious indices
			hER=hmm[ER][iq][iw][iplt]
			hSR=hmm[SR][iq][iw][iplt]
			hER.Sumw2()
			hSR.Sumw2()
			#! First draw hER
			hER.SetXTitle("%s [%s]"%(PLT_TITLE[iplt], PLT_UNIT[iplt]))
			cmm[iplt][CEIAT].SetLeftMargin(0.20)
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

			#! Create and add entries: hists and cut-lines
			#! legend
			lgnd=ROOT.TLegend(0.66,0.72,0.9,0.9)#,"","NDC");	
			for hist,label in zip([hER,hSR],["exp","sim"]):
				lgnd.AddEntry(hist,label,"lp")

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
				lgnd.AddEntry(MM_T2_CUTL_EI,"%.2f GeV < MM < %.2f GeV"%(EI_MML,EI_MMH),"l")
				lgnd.AddEntry(MM_T2_CUTL_AT,"%.2f GeV < MM < %.2f GeV"%(AT_MML,AT_MMH),"l")
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
				lgnd.AddEntry(MM2_T2_CUTL_EI,"%.2f GeV^{2} < MM^{2} < %.2f GeV^{2}"%(EI_MM2L,EI_MM2H),"l")
				lgnd.AddEntry(MM2_T2_CUTL_AT,"%.2f GeV^{2} < MM^{2} < %.2f GeV^{2}"%(AT_MM2L,AT_MM2H),"l")
			#! Draw legend
			lgnd.Draw("same")
			
			# #! [07-13-17] fill SE(%)
			#! Firt prepare to Fit ER and SR and obtain data from fit
			fit_lims=[0 for i in range(NDTYP)]
			if iplt==MM:
				#! fit limits
				for idtyp in range(NDTYP):
					fit_lims[idtyp]=[MM_FIT_LIMS[idtyp][iw][0],MM_FIT_LIMS[idtyp][iw][1]]
				#! MM EI limits
				xmin_EI=EI_MML
				xmax_EI=EI_MMH
				bin_xmin_EI=hER.FindBin(EI_MML)
				bin_xmax_EI=hER.FindBin(EI_MMH)	
				#! MM AT limits
				xmin_AT=AT_MML
				xmax_AT=AT_MMH
				bin_xmin_AT=hER.FindBin(AT_MML)
				bin_xmax_AT=hER.FindBin(AT_MMH)	
				#! xmin,xmax of histogram
				xmin=hER.GetXaxis().GetXmin()
				xmax=hER.GetXaxis().GetXmax()
				bin_xmin=hER.FindBin(xmin)
				bin_xmax=hER.FindBin(xmax)
			elif iplt==MM2:
				#! fit limits
				for idtyp in range(NDTYP):
					fit_lims[idtyp]=[MM2_FIT_LIMS[idtyp][iw][0],MM2_FIT_LIMS[idtyp][iw][1]]
				#! MM2 EI limits
				xmin_EI=EI_MM2L
				xmax_EI=EI_MM2H
				bin_xmin_EI=hER.FindBin(EI_MM2L)
				bin_xmax_EI=hER.FindBin(EI_MM2H)
				#! MM2 AT limits
				xmin_AT=AT_MM2L
				xmax_AT=AT_MM2H
				bin_xmin_AT=hER.FindBin(AT_MM2L)
				bin_xmax_AT=hER.FindBin(AT_MM2H)
				#! xmin,xmax of histogram
				xmin=hER.GetXaxis().GetXmin()
				xmax=hER.GetXaxis().GetXmax()
				bin_xmin=hER.FindBin(xmin)
				bin_xmax=hER.FindBin(xmax)
			#! Fit
			hER.Fit("gaus","+0","",fit_lims[ER][0],fit_lims[ER][1])#fit_min,fit_max)
			hSR.Fit("gaus","+0","",fit_lims[SR][0],fit_lims[SR][1])#fit_min,fit_max)
			#! Draw fit functions and fit limits
			fER=hER.GetFunction("gaus")
			fSR=hSR.GetFunction("gaus")
			if fER!=None:
				fER.SetLineColor(CLRS_DTYP[ER])
				fER.Draw("same")
				#! draw fit lims
				fitlimlER=ROOT.TLine(fit_lims[ER][0],hER.GetMinimum(),fit_lims[ER][0],hER.GetMaximum())
				fitlimhER=ROOT.TLine(fit_lims[ER][1],hER.GetMinimum(),fit_lims[ER][1],hER.GetMaximum())
				for ln in [fitlimlER,fitlimhER]:
					ln.SetLineColor(CLRS_DTYP[ER])
					ln.Draw("same")
				#! Draw fit functions over full range
				fERf=ROOT.TF1("gausf_ER","gaus",xmin,xmax)
				fERf.SetParameters(fER.GetParameter(0),fER.GetParameter(1),fER.GetParameter(2))
				fERf.SetLineColor(ROOT.gROOT.ProcessLine("kCyan"))
				#fERf.SetLineStyle(2)
				fERf.Draw("same")
			else:
				print "ERROR: fER=None for plt=%s, q2wbin="%(PLT_NAME[iplt], q2min,q2max,wmin,wmax)
			if fSR!=None:
				fSR.SetLineColor(CLRS_DTYP[SR])
				fSR.Draw("same")
				#! draw fit lims
				fitlimlSR=ROOT.TLine(fit_lims[SR][0],hER.GetMinimum(),fit_lims[SR][0],hER.GetMaximum())
				fitlimhSR=ROOT.TLine(fit_lims[SR][1],hER.GetMinimum(),fit_lims[SR][1],hER.GetMaximum())
				for l in [fitlimlSR,fitlimhSR]:
					l.SetLineColor(CLRS_DTYP[SR])
					l.Draw("same")
				#! Draw fit functions over full range
				fSRf=ROOT.TF1("gausf_SR","gaus",xmin,xmax)
				fSRf.SetParameters(fSR.GetParameter(0),fSR.GetParameter(1),fSR.GetParameter(2))
				fSRf.SetLineColor(ROOT.gROOT.ProcessLine("kMagenta"))
				#fSRf.SetLineStyle(2)
				fSRf.Draw("same")
			else:
				print "ERROR: fSR=None for plt=%s, q2wbin="%(PLT_NAME[iplt], q2min,q2max,wmin,wmax)
			
			#! Now calculate SE
			#! EI
			# print "Integral-%s (EI) for ER,ER-full,SR,SR-full for %.2f-%.2f_%.3f-%.3f"%(PLT_NAME[iplt], q2min,q2max,wmin,wmax)
			# print fER.Integral(xmin_EI,xmax_EI)
			# print fER.Integral(xmin,xmax)
			# print fSR.Integral(xmin_EI,xmax_EI)
			# print fSR.Integral(xmin,xmax)
			#sys.exit()
			if fER!=None and fSR!=None:
				x_EI=fER.Integral(xmin_EI,xmax_EI)/fER.Integral(xmin,xmax)
				y_EI=fSR.Integral(xmin_EI,xmax_EI)/fSR.Integral(xmin,xmax)
				se_EI=(1-(x_EI/y_EI))*100
			else:
				se_EI=-9999
			#! AT
			if fER!=None and fSR!=None:
				x_AT=fER.Integral(xmin_AT,xmax_AT)/fER.Integral(xmin,xmax)
				y_AT=fSR.Integral(xmin_AT,xmax_AT)/fSR.Integral(xmin,xmax)
				se_AT=(1-(x_AT/y_AT))*100
			else:
				se_AT=-9999
			#! Fill BG structure
			SE[iq][iw][iplt][EI]=se_EI
			SE[iq][iw][iplt][AT]=se_AT

			#! save canvas
			cmm[iplt][CEIAT].SaveAs("%s/%s.png"%(outdir_q2w_plt_EIAT,cname))
			cmm[iplt][CEIAT].SaveAs("%s/%s.pdf"%(outdir_q2w_plt_EIAT,cname))

			#! Now plot CEI
			cmm[iplt][CEI]=ROOT.TCanvas(cname,cname)
			#! First draw hER
			hER.Draw()
			cmm[iplt][CEI].SetLeftMargin(0.20)
			hSR.Draw("same")

			#! Create and add entries: hists and cut-lines
			#! legend
			lgnd=ROOT.TLegend(0.66,0.72,0.9,0.9)#,"","NDC");	
			for hist,label in zip([hER,hSR],["exp","sim"]):
				lgnd.AddEntry(hist,label,"lp")

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
				#! add to legend
				lgnd.AddEntry(MM_T2_CUTL_EI,"%.2f GeV < MM < %.2f GeV"%(EI_MML,EI_MMH),"l")
			elif iplt==MM2:
				#! EI
				MM2_T2_CUTL_EI.SetY1(hER.GetMinimum())
				MM2_T2_CUTL_EI.SetY2(hER.GetMaximum())
				MM2_T2_CUTH_EI.SetY1(hER.GetMinimum())
				MM2_T2_CUTH_EI.SetY2(hER.GetMaximum())
				MM2_T2_CUTL_EI.Draw("same")
				MM2_T2_CUTH_EI.Draw("same")
				#! add to legend
				lgnd.AddEntry(MM2_T2_CUTL_EI,"%.2f GeV^{2} < MM^{2} < %.2f GeV^{2}"%(EI_MM2L,EI_MM2H),"l")
			#! Draw legend
			lgnd.Draw("same")
			
			#! Draw fit functions and fit lims
			fER=hER.GetFunction("gaus")
			fSR=hSR.GetFunction("gaus")
			if fER!=None:
				fER.SetLineColor(CLRS_DTYP[ER])
				#fER.SetLineWidth(1)
				fER.Draw("same")
				#! draw fit lims
				fitlimlER=ROOT.TLine(fit_lims[ER][0],hER.GetMinimum(),fit_lims[ER][0],hER.GetMaximum())
				fitlimhER=ROOT.TLine(fit_lims[ER][1],hER.GetMinimum(),fit_lims[ER][1],hER.GetMaximum())
				for l in [fitlimlER,fitlimhER]:
					l.SetLineColor(CLRS_DTYP[ER])
					l.Draw("same")
				#! Draw fit functions over full range
				fERf=ROOT.TF1("gausf_ER","gaus",xmin,xmax)
				fERf.SetParameters(fER.GetParameter(0),fER.GetParameter(1),fER.GetParameter(2))
				fERf.SetLineColor(ROOT.gROOT.ProcessLine("kCyan"))#(CLRS_DTYP[ER])
				fERf.SetLineWidth(1)
				#fERf.SetLineStyle(2)
				fERf.Draw("same")
			else:
				print "ERROR: fER=None for plt=%s, q2wbin="%(PLT_NAME[iplt], q2min,q2max,wmin,wmax)
			if fSR!=None:
				fSR.SetLineColor(CLRS_DTYP[SR])
				#fSR.SetLineWidth(1)
				#! draw fit lims
				fitlimlSR=ROOT.TLine(fit_lims[SR][0],hER.GetMinimum(),fit_lims[SR][0],hER.GetMaximum())
				fitlimhSR=ROOT.TLine(fit_lims[SR][1],hER.GetMinimum(),fit_lims[SR][1],hER.GetMaximum())
				for ln in [fitlimlSR,fitlimhSR]:
					ln.SetLineColor(CLRS_DTYP[SR])
					ln.Draw("same")
				fSR.Draw("same")
				#! Draw fit functions over full range
				fSRf=ROOT.TF1("gausf_SR","gaus",xmin,xmax)
				fSRf.SetParameters(fSR.GetParameter(0),fSR.GetParameter(1),fSR.GetParameter(2))
				fSRf.SetLineColor(ROOT.gROOT.ProcessLine("kMagenta"))#!CLRS_DTYP[SR])
				fSRf.SetLineWidth(1)
				#fSRf.SetLineStyle(2)
				fSRf.Draw("same")
			else:
				print "ERROR: fSR=None for plt=%s, q2wbin="%(PLT_NAME[iplt], q2min,q2max,wmin,wmax)
			# if fER!=None:
			# 	fER.SetLineColor(CLRS_DTYP[ER])
			# 	fER.Draw("same")
			# if fSR!=None:
			# 	fSR.SetLineColor(CLRS_DTYP[SR])
			# 	fSR.Draw("same")

			#! save canvas
			cmm[iplt][CEI].SaveAs("%s/%s.png"%(outdir_q2w_plt_EI,cname))
			cmm[iplt][CEI].SaveAs("%s/%s.pdf"%(outdir_q2w_plt_EI,cname))
			#cmm[iplt][CEI].SaveAs("%s/%s.eps"%(outdir_q2w_plt_EI,cname))


#! ***
#sys.exit()

# Now consolidate data from SE for plotting
#! + Basically arrange data to plotted vs W in bins for Q2

#! Some related aesthetics
#! LABEL[plt,cut]
LABEL=[[0 for j in range(NCUTS)] for i in range(NPLTS)] 
LABEL[MM][EI]="%.2f GeV < MM < %.2f GeV"%(EI_MML,EI_MMH)
LABEL[MM][AT]="%.2f GeV < MM < %.2f GeV"%(AT_MML,AT_MMH)
LABEL[MM2][EI]="%.2f GeV < " r"$MM^{2}$" " < %.2f " r"$GeV^{2}$"%(EI_MM2L,EI_MM2H)
LABEL[MM2][AT]="%.2f GeV < " r"$MM^{2}$" " < %.2f " r"$GeV^{2}$"%(AT_MM2L,AT_MM2H)
for iq,q2bin_le in enumerate(Q2BIN_LEL):
	if DBG==True and iq>0: continue #! debug
	q2min=q2bin_le
	q2max=Q2BIN_UEL[iq]
	outdir_q2w="%s/%.2f-%.2f"%(OUTDIR,q2min,q2max)
	if not os.path.exists(outdir_q2w):
		os.makedirs(outdir_q2w)
	#! set q2bin label
	q2binl='[%.2f,%.2f)'%(q2min,q2max)
	for iplt in range(NPLTS):
		outdir_q2w_plt_EI="%s/EI/%s"%(outdir_q2w,PLT_NAME[iplt])
		outdir_q2w_plt_EIAT="%s/EI_AT/%s"%(outdir_q2w,PLT_NAME[iplt])
		for outdir in [outdir_q2w_plt_EI,outdir_q2w_plt_EIAT]:
			if not os.path.exists(outdir):
				os.makedirs(outdir)
		#! structure to consolidate W data for q2,plt
		se_EI,se_AT=[],[]
		wbinn=[] #! wbin number
		wbinl=[] #! wbin label
		for iw,wbin_le in enumerate(WBIN_LEL):
			#!DBGW
			if DBG==True and (iw!=0 and iw!=15 and iw!=28): continue #! debug #!(iw!=0 and iw!=15 and iw!=28)
			#if DBG==True and iw+1 > 10: continue
			#if DBG==True and (iw+1!=5 and iw+1!=6): continue
			wmin=wbin_le
			wmax=wbin_le+WBINW
			se_EI.append(SE[iq][iw][iplt][EI])
			se_AT.append(SE[iq][iw][iplt][AT])
			wbinn.append(iw+1)
			wbinl.append('[%.3f,%.3f)'%(wmin,wmax))

		#! plot
		#! SE for EI and AT together
		fig=plt.figure()
		ax=plt.subplot(111)
		#! To give space below x axis for rotated W bin labels
		fig.subplots_adjust(bottom=0.25)
		ax.scatter(wbinn,se_EI,c='m',      label=LABEL[iplt][EI])
		ax.scatter(wbinn,se_AT,c='orange', label=LABEL[iplt][AT])
		#! Fix x axis
		# shift half a step to the left
		xmin=(3*wbinn[0]-wbinn[1])/2.
		# shift half a step to the right
		xmax=(3*wbinn[-1]-wbinn[-2])/2.
		ax.set_xlim(xmin,xmax)
		ax.set_xticks(wbinn)
		ax.set_xticklabels(wbinl,rotation='vertical')
		if iplt==MM:
			ax.set_title("Effect on the Measured Cross-sections due to the different" "\n" "Missing Mass Distributions in Experiment and Simulation") #for " r" $Q^{2}$" "=[%.2f,%.2f) " r"$GeV^{2}$"%(q2min,q2max))
		elif iplt==MM2:
			ax.set_title("Effect on the Measured Cross-sections due to the different" "\n" "Missing Mass Squared Distributions in Experiment and Simulation")# for " r" $Q^{2}$" "=[%.2f,%.2f) " r"$GeV^{2}$"%(q2min,q2max))
		ax.set_xlabel("W (GeV)")
		ax.set_ylabel(r"$\frac{\delta \sigma}{\sigma} (\%)$",fontsize=25)
		#! draw legend
		ax.legend(loc='upper right',prop={'size':8})
		#! save figure
		fig.savefig('%s/SE.png'%(outdir_q2w_plt_EIAT))
		fig.savefig('%s/SE.pdf'%(outdir_q2w_plt_EIAT))

		#! output as text
		fout=open("%s/SE.txt"%outdir_q2w_plt_EIAT,'w')
		fout.write("***EI***\n")
		for i in range(len(wbinn)):
			fout.write("%d %s %.2f\n"%(i+1,wbinl[i],se_EI[i]))
		fout.write("\n")
		fout.write("***AT***\n")
		for i in range(len(wbinn)):
			fout.write("%d %s %.2f\n"%(i+1,wbinl[i],se_AT[i]))
		fout.close()


		#! SE for EI only
		fig=plt.figure()
		ax=plt.subplot(111)
		#! To give space below x axis for rotated W bin labels
		fig.subplots_adjust(bottom=0.25)
		ax.scatter(wbinn,se_EI,c='m',      label=LABEL[iplt][EI])
		#! Fix x axis
		# shift half a step to the left
		xmin=(3*wbinn[0]-wbinn[1])/2.
		# shift half a step to the right
		xmax=(3*wbinn[-1]-wbinn[-2])/2.
		ax.set_xlim(xmin,xmax)
		ax.set_xticks(wbinn)
		ax.set_xticklabels(wbinl,rotation='vertical')
		if iplt==MM:
			ax.set_title("Effect on the Measured Cross-sections due to the different" "\n" "Missing Mass Distributions in Experiment and Simulation")# for " r" $Q^{2}$" "=[%.2f,%.2f) " r"$GeV^{2}$"%(q2min,q2max))
		elif iplt==MM2:
			ax.set_title("Effect on the Measured Cross-sections due to the different" "\n" "Missing Mass Squared Distributions in Experiment and Simulation")# for " r" $Q^{2}$" "=[%.2f,%.2f) " r"$GeV^{2}$"%(q2min,q2max))
		ax.set_xlabel("W (GeV)")
		ax.set_ylabel(r"$\frac{\delta \sigma}{\sigma} (\%)$",fontsize=25)
		#! draw legend
		#! legend not needed for only EI
		#ax.legend(loc='upper right',prop={'size':8})
		#! save figure
		fig.savefig('%s/SE.png'%(outdir_q2w_plt_EI))
		fig.savefig('%s/SE.pdf'%(outdir_q2w_plt_EI))

		#! output as text
		fout=open("%s/SE.txt"%outdir_q2w_plt_EI,'w')
		fout.write("***EI***\n")
		for i in range(len(wbinn)):
			fout.write("%d %s %.2f\n"%(i+1,wbinl[i],se_EI[i]))
		fout.close()


		


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
	
