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
+ {e1f,e16}*{ER,SR}*{W-elastic,MM-t2}*{wpcorr,npcorr}

+ Usage: study_pcorr.py expt=<e1f/e16> plot=<pcorr/MM> debug=False
'''

USAGE='study_pcorr.py expt=<e1f/e16> plot=<pcorr/MM> debug=[False]'

#! Get input data from user
if len(sys.argv)<2:
		sys.exit('usage: %s'%USAGE)
expt=sys.argv[1]
if expt!='e1f' and expt!='e16':
	sys.exit('expt=%s is not valid. usage: %s'%(expt,USAGE))
#! For now ignore e1f
if expt=="e1f":
	sys.exit("Not implemented for e1f yet")

if len(sys.argv)<3:
	sys.exit('usage: %s'%USAGE)
PLOT=sys.argv[2]
if PLOT!='pcorr' and PLOT!="MM":
	sys.exit('PLOT=%s is not valid. usage: %s'%(PLOT,USAGE))

DBG=False
if len(sys.argv)>3: #i.e. debug entered by user
	if    sys.argv[3]=="True": DBG=True
	elif  sys.argv[3]=="False": DBG=False
	else: sys.exit('DBG=%s is not valid. usage: %s'%(DBG,USAGE))


print "expt=",expt
print "PLOT=",PLOT
print "DBG=",DBG
#sys.exit()


#! *** Prepare input data structure: FIN/T[dtyp][rctn][corr] ***
NDTYP=2
ER,SR=range(NDTYP)
DTYP_NAME=["ER","SR"]

NRCTN=2
TPI,ELAST=range(NRCTN)
RCTN_NAME=["2pi","elast"]

NCORR=2
WPCORR,NPCORR=range(NCORR)
CORR_NAME=['wpcorr','npcorr']
CORR_LEGEND_TITLE={'e1f':'e pcorr','e16':'e,p,#pi^{+} pcorr'}

FIN=[[[[] for k in range(NCORR)] for j in range(NRCTN)] for i in range(NDTYP)]
T=[[[[] for k in range(NCORR)] for j in range(NRCTN)] for i in range(NDTYP)]
#print "FIN=",FIN
#print "T=",T
#! ***

#! *** Prepare output datadir***
if   expt=="e1f": DATADIR_OUTPUT=os.path.join(os.environ['STUDY_PCORR_E1F_DATADIR'],PLOT)
elif expt=="e16": DATADIR_OUTPUT=os.path.join(os.environ['STUDY_PCORR_E16_DATADIR'],PLOT)
if DBG==True: DATADIR_OUTPUT+="/dbg"
print "DATADIR_OUTPUT=",DATADIR_OUTPUT
#sys.exit()
#! ***

#! *** Prepare NENTRIES ***
NENTRIES=1000000000
if DBG==True: NENTRIES=1000

#! *** Now fill T[dtyp][rctn][corr] ***
DATADIR_INPUT=[[] for i in range(NDTYP)]
sfx=''
if expt=='e16':sfx='_E16'
DATADIR_INPUT[ER]=os.path.join(os.environ['D2PIDIR_EXP%s'%sfx],"data_pcorr_011116")#! debug /bk
DATADIR_INPUT[SR]=os.path.join(os.environ['D2PIDIR_SIM%s'%sfx],"data_pcorr_011116")#! debug /bl

print "DATADIR_INPUT[ER]=",DATADIR_INPUT[ER]
print "DATADIR_INPUT[SR]=",DATADIR_INPUT[SR]

d=[[ER,SR],[TPI,ELAST],[WPCORR,NPCORR]]
D=list(itertools.product(*d))
for r in D:
	idtyp,irctn,icorr=r[0],r[1],r[2]
	if DTYP_NAME[idtyp]=="SR" and CORR_NAME[icorr]=="wpcorr": continue #! No {SR}*{wpcorr}
	if expt=="e16" and DTYP_NAME[idtyp]=="SR" and RCTN_NAME[irctn]=="elast": continue #! No e16*SR*delastR	
	print "Getting T for",DTYP_NAME[idtyp],RCTN_NAME[irctn],CORR_NAME[icorr]
	FIN[idtyp][irctn][icorr]=ROOT.TFile("%s/d%sR_%s.root"%(DATADIR_INPUT[idtyp],RCTN_NAME[irctn],CORR_NAME[icorr]))
	if RCTN_NAME[irctn]=="elast": T[idtyp][irctn][icorr]=FIN[idtyp][irctn][icorr].Get("delast/monitor/t")
	elif RCTN_NAME[irctn]=="2pi": T[idtyp][irctn][icorr]=FIN[idtyp][irctn][icorr].Get("d2piR/tR")
	#print FIN[idtyp][irctn][icorr].GetName()
	#print T[idtyp][irctn][icorr].GetName()
#print "T=",T
#sys.exit()
#! ***

#! *** Now make plots ***

#! sector information
NSCTR=6

#! *** First make W-elastic plots ***
study_W_elast=False#True #! debug
if PLOT=='pcorr': study_W_elast=True
if study_W_elast:
	#! outdir
	outdir=os.path.join(DATADIR_OUTPUT,"W-elast")
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	#! draw_cmd
	draw_cmd="W>>hdcmd(100,0.7,1.2)"
	#! Create h[dtyp][sctr][corr]
	h=[[[[] for k in range(NCORR)]for j in range(NSCTR)]for i in range(NDTYP)]
	print "h=",h
	#! fill h
	for r in D:
		idtyp,irctn,icorr=r[0],r[1],r[2]
		if RCTN_NAME[irctn]=="2pi": continue
		if DTYP_NAME[idtyp]=="SR" and CORR_NAME[icorr]=="wpcorr": continue #! No {SR}*{wpcorr}
		if expt=="e16" and DTYP_NAME[idtyp]=="SR" and RCTN_NAME[irctn]=="elast": continue #! No e16*SR*delastR

		for isctr in range(NSCTR):
			cut=ROOT.TCut("sector==%d"%(isctr+1))
			print "TTree::Draw() for",DTYP_NAME[idtyp],RCTN_NAME[irctn],CORR_NAME[icorr]
			#print T[idtyp][irctn][icorr].GetName()
			T[idtyp][irctn][icorr].Draw(draw_cmd,cut,"",NENTRIES)#[idtyp][isctr][icorr],cut)

			#! Store histogram
			#ROOT.gStyle.SetOptStat("n")#! This line is critically needed here!
			htmp=ROOT.gDirectory.Get("hdcmd")#draw_cmd_hst_name[idtyp][isctr][icorr])
			h[idtyp][isctr][icorr]=htmp.Clone()
			h[idtyp][isctr][icorr].SetName("h_%s"%CORR_NAME[icorr])
			h[idtyp][isctr][icorr].SetTitle("%s"%cut.GetTitle())
			#print h[idtyp][isctr][icorr].GetName()
	#! Draw and save plots
	ROOT.gStyle.SetOptStat("n")
	ROOT.gStyle.SetOptFit(1)#111)
	cname="W"
	celast=ROOT.TCanvas(cname,cname,1300,800)
	celast.Divide(3,2)
	for isctr in range(NSCTR):
		pad=celast.cd(isctr+1)
		#cname="sector%d"%(isctr+1)
		#celast=ROOT.TCanvas(cname,cname,1300,800)
		#l=ROOT.TLegend(0.1,0.3,0.3,0.4)#,"","NDC");
		if expt=="e1f":#! both {ER,SR} exist for e1f
			sys.exit("Not implemented for e1f for now")
		elif expt=="e16": #! only {ER} exist for e16
			h[ER][isctr][WPCORR].SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
			h[ER][isctr][WPCORR].Fit("gaus","+0","",0.87,0.98)
			h[ER][isctr][NPCORR].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
			if isctr+1==1 or isctr+1==4 or isctr+1==5 or isctr+1==6: #! s1,s4 have a larger shift (This can also be seen in KP's thesis p86)
				h[ER][isctr][NPCORR].Fit("gaus","+0","",0.87,1.00)
			else:
				h[ER][isctr][NPCORR].Fit("gaus","+0","",0.87,0.98)
			h[ER][isctr][WPCORR].Draw()
			#! Axes title
                        h[ER][isctr][WPCORR].SetXTitle("W [GeV]")
                        celast.SetLeftMargin(0.20)
                        h[ER][isctr][WPCORR].GetYaxis().SetTitleOffset(1.5)
                        h[ER][isctr][WPCORR].SetYTitle("N_{entries}")
			#! other hists
			h[ER][isctr][NPCORR].Draw("sames")
			#! Draw the fit funtions
			f_wpcorr=h[ER][isctr][WPCORR].GetFunction("gaus")
			f_npcorr=h[ER][isctr][NPCORR].GetFunction("gaus")
			if f_wpcorr!=None: 
				f_wpcorr.SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
				f_wpcorr.Draw("same")
			if f_npcorr!=None: 
				f_npcorr.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
				f_npcorr.Draw("same")
			#f_wpcorr.Draw("same")
			#f_npcorr.Draw("same")
			#! Ajust statbox
			#celast.Update()
			pad.Update()
			#! statbox for wpcorr
			pt_wpcorr=h[ER][isctr][WPCORR].GetListOfFunctions().FindObject("stats")
                        pt_wpcorr.SetTextColor(ROOT.gROOT.ProcessLine("kBlue"))
			#! statbox for npcorr
			pt_npcorr=h[ER][isctr][NPCORR].GetListOfFunctions().FindObject("stats")
                        pt_npcorr.SetTextColor(ROOT.gROOT.ProcessLine("kBlack"))
                        diff=pt_wpcorr.GetY2NDC()-pt_wpcorr.GetY1NDC()
                        y2=pt_wpcorr.GetY1NDC()
                        y1=y2-diff
                        pt_npcorr.SetY1NDC(y1)
                        pt_npcorr.SetY2NDC(y2)
			#celast.Update()
			pad.Update()
			#! Draw legend (not necessary because statbox contain relevant information)
			#l.AddEntry(h[ER][isctr][WPCORR],'wpcorr')
			#l.AddEntry(h[ER][isctr][NPCORR],'npcorr')
			#l.Draw()
		#celast.Update()
		pad.Update()
		#celast.SaveAs("%s/%s.png"%(outdir,cname))
	celast.SaveAs("%s/%s.png"%(outdir,cname))
	celast.SaveAs("%s/%s.pdf"%(outdir,cname))
#! ***
#! debug
#sys.exit()
	
#! *** Now make MM plots ***
WMIN=1.400
WMAX=2.125
WBINW=0.025
Q2BIN_LEL=[2.00,2.40,3.00,3.50,4.20,2.00]#[1.50,2.00,2.40,3.00,3.50,4.20] #debug
Q2BIN_UEL=[2.40,3.00,3.50,4.20,5.00,5.00]#[2.00,2.40,3.00,3.50,4.20,5.00] #debug
WBIN_LEL=np.arange(WMIN,WMAX,WBINW)
#print Q2BIN_LEL
#print WBIN_LEL

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
#! TLine depicting pion mass debug
MASS_PION=0.13957018;
MASS_PION_LINE=ROOT.TLine(MASS_PION,0,MASS_PION,0)
MASS_PION_LINE.SetLineColor(ROOT.gROOT.ProcessLine("kGreen"))
MASS_PION_LINE.SetLineWidth(3)

#! outdir
outdir=os.path.join(DATADIR_OUTPUT,"2pi-MM-top2")
if not os.path.exists(outdir):
	os.makedirs(outdir)
#! draw_cmd
if PLOT=='pcorr': draw_cmd="mmppip>>hmmcmd(100,0,0.4)" #! make hist x-range narrower to bring out pcorr
elif PLOT=='MM':  draw_cmd="mmppip>>hmmcmd(200,-0.5,1)" #! make hist x-range wider to show validity of MMcut
#draw_cmd="mmppip>>hmmcmd(100,0,0.4)"#hmmcmd(200,-0.5,1) #debug
#! Create h[dtyp][q2][w][corr]
hmm=[[[[[]for l in range(NCORR)]for k in range(len(WBIN_LEL))]for j in range(len(Q2BIN_LEL))] for i in range(NDTYP)]
#print "hmm=",hmm
#! fill h
for i,q2bin_le in enumerate(Q2BIN_LEL):
	if DBG==True and i>0: continue #! debug
	q2min=q2bin_le
        q2max=Q2BIN_UEL[i]
	for j,wbin_le in enumerate(WBIN_LEL):
		wmin=wbin_le
		wmax=wbin_le+WBINW
		cut_q2w=ROOT.TCut("Q2>%f && Q2<%f && W>%f && W<%f"%(q2min,q2max,wmin,wmax))
		print "q2w=",cut_q2w.GetTitle()
		for r in D:
			idtyp,irctn,icorr=r[0],r[1],r[2]
			if RCTN_NAME[irctn]=="elast": continue
			if DTYP_NAME[idtyp]=="SR" and CORR_NAME[icorr]=="wpcorr": continue #! No {SR}*{wpcorr}
			if expt=="e16" and DTYP_NAME[idtyp]=="SR" and RCTN_NAME[irctn]=="elast": continue #! No e16*SR*delastR

			cut_top=ROOT.TCut("top==2")
			cut=ROOT.TCut(cut_q2w)
			cut+=cut_top
			print "TTree::Draw() for",DTYP_NAME[idtyp],RCTN_NAME[irctn],CORR_NAME[icorr]
			print T[idtyp][irctn][icorr].GetName()
			T[idtyp][irctn][icorr].Draw(draw_cmd,cut,"",NENTRIES)#[idtyp][isctr][icorr],cut)

			#! Store histogram
			#ROOT.gStyle.SetOptStat("ne")
			htmp=ROOT.gDirectory.Get("hmmcmd")#(draw_cmd_hst_name[idtyp][isctr][icorr])
			hmm[idtyp][i][j][icorr]=htmp.Clone()#"hmm_%s_%s"%(DTYP_NAME[idtyp],CORR_NAME[icorr]))
			hmm[idtyp][i][j][icorr].SetName("hmm_%s_%s"%(DTYP_NAME[idtyp],CORR_NAME[icorr]))
			hmm[idtyp][i][j][icorr].SetTitle("MM-top2 %.2f-%.2f_%.3f-%.3f"%(q2min,q2max,wmin,wmax))#("%s"%cut.GetTitle())
			#print hmm[idtyp][i][j][icorr].GetName()
#! Plot and save hmm
ROOT.gStyle.SetOptStat("n")
ROOT.gStyle.SetOptFit(1)#111)
for i,q2bin_le in enumerate(Q2BIN_LEL):
	if DBG==True and i>0: continue #! debug
	q2min=q2bin_le
        q2max=Q2BIN_UEL[i]
	outdir_q2w="%s/%.2f-%.2f"%(outdir,q2min,q2max)
	if not os.path.exists(outdir_q2w):
		os.makedirs(outdir_q2w)
	for j,wbin_le in enumerate(WBIN_LEL):
		wmin=wbin_le
		wmax=wbin_le+WBINW
		
		print "Going to plot hmm for %.2f-%.2f_%.3f-%.3f"%(q2min,q2max,wmin,wmax)

		cname="c_%.3f-%.3f"%(wmin,wmax)
		cmm=ROOT.TCanvas(cname,cname)
		#l=ROOT.TLegend(0.1,0.3,0.3,0.4)#,"","NDC");
		#! hists aesthetics	
		hmm[ER][i][j][WPCORR].SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
		hmm[ER][i][j][NPCORR].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
		hmm[SR][i][j][NPCORR].SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
		#! Draw hists
		hER_MM=None #! Use this to easily select whether to use hmm_ER_WPCORR or hmm_ER_NPCORR for illustrating MM (NPCORR compares better with sim)
		if PLOT=='MM': #! then for now draw ER-NPCORR and SR-NPCORR since ER-NPCORR agrees better with SR-NPCORR and together they justify the MM-cut
			hER_MM=hmm[ER][i][j][NPCORR] #! hmm[ER][i][j][WPCORR]
			hER_MM.Draw()
			#! Axes title
			hER_MM.SetXTitle("MM [GeV]")
                        cmm.SetLeftMargin(0.20)
                        hER_MM.GetYaxis().SetTitleOffset(1.5)
                        hER_MM.SetYTitle("N_{entries}")	
			#! scale SR and then draw
			max_ER=hER_MM.GetMaximum()
        	        if max_ER==0:max_ER=1
                	max_SR=hmm[SR][i][j][NPCORR].GetMaximum()
                	if max_SR==0:max_SR=1
                	scl_fctr_SR=max_ER/max_SR
                	if scl_fctr_SR==0:scl_fctr_SR=1
                	hmm[SR][i][j][NPCORR].Scale(scl_fctr_SR)
                	hmm[SR][i][j][NPCORR].Draw("sames")
		elif PLOT=='pcorr': #!then draw ER-NPCORR and ER-WPCORR to show effect of pcorr
			hmm[ER][i][j][WPCORR].Draw()
			#! Axes title
                        hmm[ER][i][j][WPCORR].SetXTitle("MM [GeV]")
                        cmm.SetLeftMargin(0.20)
                        hmm[ER][i][j][WPCORR].GetYaxis().SetTitleOffset(1.5)
                        hmm[ER][i][j][WPCORR].SetYTitle("N_{entries}")
			#! Now overlay NPCORR
			hmm[ER][i][j][NPCORR].Draw("sames") #! Note only "sames" seems to work here (with "same" statbox with name=hxx is drawn?!)
		#! Adjust statbox
		cmm.Update()
		if PLOT=='pcorr': #! then adjust accordingly statbox for ER:wpcorr
			pt_ER_npcorr=hmm[ER][i][j][NPCORR].GetListOfFunctions().FindObject("stats")
	                pt_ER_npcorr.SetTextColor(ROOT.gROOT.ProcessLine("kBlack"))
                	pt_ER_wpcorr=hmm[ER][i][j][WPCORR].GetListOfFunctions().FindObject("stats")
               		pt_ER_wpcorr.SetTextColor(ROOT.gROOT.ProcessLine("kBlue"))
                	diff=pt_ER_npcorr.GetY2NDC()-pt_ER_npcorr.GetY1NDC()
                	y2=pt_ER_npcorr.GetY1NDC()
                	y1=y2-diff
                	pt_ER_wpcorr.SetY1NDC(y1)
                	pt_ER_wpcorr.SetY2NDC(y2)
		if PLOT=='MM': #! then adjust accordingly statbox for SR:npcorr
			pt_ER=hER_MM.GetListOfFunctions().FindObject("stats")
			pt_ER.SetTextColor(hER_MM.GetLineColor())
			#! statbox for SR:npcorr
			pt_SR_npcorr=hmm[SR][i][j][NPCORR].GetListOfFunctions().FindObject("stats")
                	pt_SR_npcorr.SetTextColor(ROOT.gROOT.ProcessLine("kRed"))
                	diff=pt_ER.GetY2NDC()-pt_ER.GetY1NDC()
                	y2=pt_ER.GetY1NDC()
                	y1=y2-diff
                	pt_SR_npcorr.SetY1NDC(y1)
                	pt_SR_npcorr.SetY2NDC(y2)
                cmm.Update()
		if PLOT=='MM':
			#! Draw cut lines
			MM_T2_CUTL_EI.SetY1(hmm[ER][i][j][WPCORR].GetMinimum())
                	MM_T2_CUTL_EI.SetY2(hmm[ER][i][j][WPCORR].GetMaximum())
                	MM_T2_CUTH_EI.SetY1(hmm[ER][i][j][WPCORR].GetMinimum())
                	MM_T2_CUTH_EI.SetY2(hmm[ER][i][j][WPCORR].GetMaximum())
                	MM_T2_CUTL_EI.Draw("same")
                	MM_T2_CUTH_EI.Draw("same")
			#! Draw line depicting pion mass
			#MASS_PION_LINE.SetY1(hmm[ER][i][j][WPCORR].GetMinimum())
			#MASS_PION_LINE.SetY2(hmm[ER][i][j][WPCORR].GetMaximum())
			#MASS_PION_LINE.Draw("same")
		#! Draw legend (not drawing since relevant information in statbox)
		#l.AddEntry(hmm[ER][i][j][WPCORR],CORR_LEGEND_TITLE['e16'])
                #l.AddEntry(hmm[ER][i][j][NPCORR],'npcorr')
		#l.Draw()
		#print "here",hmm[ER][i][j][WPCORR].GetName(),hmm[ER][i][j][NPCORR].GetName()
	
		if PLOT=='pcorr':	
			#! Fit debug
			#! Fit limits optimized at Q2,W=[2.0,2.5],[1.650,1.675]
			hmm[ER][i][j][WPCORR].Fit("gaus","+0","",0.09,0.18)
			hmm[ER][i][j][NPCORR].Fit("gaus","+0","",0.09,0.19)
			hmm[SR][i][j][NPCORR].Fit("gaus","+0","",0.08,0.19)
			#! Draw the fit funtions
                	f_ER_wpcorr=hmm[ER][i][j][WPCORR].GetFunction("gaus")
			f_ER_npcorr=hmm[ER][i][j][NPCORR].GetFunction("gaus")
			f_SR_npcorr=hmm[SR][i][j][NPCORR].GetFunction("gaus")
                	if f_ER_wpcorr!=None: 
				f_ER_wpcorr.SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
				f_ER_wpcorr.Draw("same")
                	if f_ER_npcorr!=None: 
				f_ER_npcorr.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
				f_ER_npcorr.Draw("same")
			#if f_SR_npcorr!=None: 
			#	f_SR_npcorr.SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
			#	f_SR_npcorr.Draw("same")
                	#f_ER_wpcorr.Draw("same")
                	#f_ER_npcorr.Draw("same")
			#f_SR_npcorr.Draw("same")

		cmm.SaveAs("%s/%s.png"%(outdir_q2w,cname))
		cmm.SaveAs("%s/%s.pdf"%(outdir_q2w,cname))
#! ***

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
	
