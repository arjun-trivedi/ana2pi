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
+ 

+ Usage: 
'''

USAGE=''

DBG=False
if len(sys.argv)>1: #i.e. debug entered by user
	if    sys.argv[1]=="True": DBG=True
	elif  sys.argv[1]=="False": DBG=False
	else: sys.exit('DBG=%s is not valid. usage: %s'%(DBG,USAGE))

print "DBG=",DBG
#sys.exit()


#! *** Prepare input data structure: FIN/T[dtyp][rctn][corr] ***
NDTYP=2
ER,SR=range(NDTYP)
DTYP_NAME=["ER","SR"]

FIN=[0 for i in range(NDTYP)]
T  =[0 for i in range(NDTYP)]
#print "FIN=",FIN
#print "T=",T
#! ***

#! *** Prepare output datadir***
OUTDIR=os.environ['STUDY_SS_BANDS_DATADIR']
if DBG==True: OUTDIR+="/dbg"
print "OUTDIR=",OUTDIR
if not os.path.exists(OUTDIR):
	os.makedirs(OUTDIR)
#! create output .root file
FOUT=ROOT.TFile("%s/study.root"%OUTDIR,"RECREATE")
#sys.exit()
#! ***

#! *** Prepare NENTRIES ***
NENTRIES=1000000000
if DBG==True: NENTRIES=1000

#! *** Now fill T[dtyp][rctn][corr] ***
INPUT=[0 for i in range(NDTYP)]
INPUT[ER]=os.path.join(os.environ['D2PIDIR_EXP_E16'],"data_SS_bands_062316")#! debug /bk
INPUT[SR]=os.path.join(os.environ['D2PIDIR_SIM_E16'],"data_SS_bands_062316")#! debug /bl

print "INPUT[ER]=",INPUT[ER]
print "INPUT[SR]=",INPUT[SR]

FIN[ER]=ROOT.TFile("%s/d2piR.root"%(INPUT[ER]))
FIN[SR]=ROOT.TFile("%s/d2piR.root"%(INPUT[SR]))
T[ER]=FIN[ER].Get("d2piR/tR")
T[SR]=FIN[SR].Get("d2piR/tR")
print "T[ER]=",T[ER]
print "T[SR]=",T[SR]
#sys.exit()
#! ***

#! *** Now make plots ***

#! Define some constants
#! Q2-W binning
Q2BIN_LEL=[2.00,2.40,3.00,3.50,4.20,2.00]#[1.50,2.00,2.40,3.00,3.50,4.20] #debug
Q2BIN_UEL=[2.40,3.00,3.50,4.20,5.00,5.00]#[2.00,2.40,3.00,3.50,4.20,5.00] #debug
#WMIN=1.400
#WMAX=2.125
#WBINW=0.025
#WBIN_LEL=np.arange(WMIN,WMAX,WBINW)
#! Make W binnin coarser for sake of stats for now
WBIN_LEL=[1.400,1.600,1.800,2.000]
WBIN_UEL=[1.600,1.800,2.000,2.125]
NWBINS=len(WBIN_LEL)
#print Q2BIN_LEL
#print WBIN_LEL
def get_wbin_idx(W):
	''' Used to get W bin's index from a given W value'''	
	found_bin=False
	for ibin,wbin_le in enumerate(WBIN_LEL):
		wmin=wbin_le
		wmax=WBIN_UEL[ibin]#wbin_le+WBINW
		if W>=wmin and W<wmax:
			found_bin=True
			return ibin
	if found_bin==False:
		print "wbin idx not found. Going to return an invalid index of -9999"
		return -9999

#! MM cuts lines
MM2_CUT_EI_MIN=-0.04
MM2_CUT_EI_MAX=0.06
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
#! TLine depicting pion mass
MASS_PION=0.13957018;
MASS_PION_LINE=ROOT.TLine(MASS_PION,0,MASS_PION,0)
MASS_PION_LINE.SetLineColor(ROOT.gROOT.ProcessLine("kGreen"))
MASS_PION_LINE.SetLineWidth(3)

#! PID
PROTON=2212;
PIP=211;

NSSBANDS=4 
SSB1,SSB2,SSB3,SSB4=range(NSSBANDS)
SSB_NAME=["SSB1","SSB2","SSB3","SSB4"] #! in order of decreasing xsec: {gpart-pid}*{stat-pid}=1. off*on, 2. off*off, 3. on*on, 4. on*off
SSB_TITLE=["gpart-pid,stat-pid=OFF,ON","gpart-pid,stat-pid=OFF,OFF","gpart-pid,stat-pid=ON,ON","gpart-pid,stat-pid=ON,OFF"]

#! Create histogram to be made
#! hmm[dtyp][ssb][w]
hmm=[[[ROOT.TH1F("hmm%s_%s_wbin%d"%(DTYP_NAME[i],SSB_NAME[j],k+1),"mmppip_%s_%s_wbin%d"%(DTYP_NAME[i],SSB_TITLE[j],k+1),200,-0.5,1) for k in range(NWBINS)] for j in range(NSSBANDS)] for i in range(NDTYP)]
#! The following are used to analyze why the SSband yields differ for ER and SR
NANAHISTS=5
GPART_PID_STAT_ON,GPART_PID_STAT_OFF,NTRK_PID,STAT_PID_P,STAT_PID_PIP=range(NANAHISTS) #! NTRK is subset of GPART containing only positive or negative GPART tracks
ANAHISTS_NAME =["gpart-pid-stat-on","gpart-pid-stat-off","ntrk-pid","stat_pid_p","stat_pid_pip"]
ANAHISTS_TITLE=["gpart-pid-stat-on","gpart-pid-stat-off","ntrk-pid","stat_pid_p","stat_pid_pip"]
ANAHISTS_NBINS=[20,20,20,20,20]
ANAHISTS_XMIN=[0, 0, 0, -10, -10]
ANAHISTS_XMAX=[20,20,20, 10,  10]
#! hana[anahist][dtyp][w]
hana=[[[ROOT.TH1F("h%s_%s_wbin%d"%(DTYP_NAME[i],ANAHISTS_NAME[j],k+1),"%s_%s_wbin%d"%(DTYP_NAME[i],ANAHISTS_TITLE[j],k+1),ANAHISTS_NBINS[j],ANAHISTS_XMIN[j],ANAHISTS_XMAX[j]) for k in range(NWBINS)] for j in range(NANAHISTS)] for i in range(NDTYP)]

#! Process of looping over TTree in Python taken from Evan's mail of [02-7-15] 'pyroot and bytes in trees'
for idtyp in [ER,SR]:
	#if idtyp==SR: continue
	for i, ev in enumerate(T[idtyp], 1):
		if (i%10000==0): print "processing ientry# %d"%i
		if DBG and i>1000: break
	
		#print "***entry #%d ***"%(i+1)

		#! Cuts
		#! 1. top cut: Use only top2 events even though top2' is actually top2+top1
		#! + This is because in top1, as cacluated by 'hvy', mmppip=0 
		#print "ev top=",ev.top
		if ev.top!=2: continue # and ev.top!=1: continue

		#! 2. W cut: to keep W between [1.400,2.125)
		if ev.W<1.400 or ev.W>=2.125: continue

		#! Get iwbin
		iwbin=get_wbin_idx(ev.W)	
		#print "ev.W, iwbin=",ev.W,iwbin	
		#print "passed top=",ev.top
		gpart_pid=ev.gpart
		ntrk_pid=ev.ntrk
		#print "gpart_pid,ntrk_pid=",gpart_pid,ntrk_pid
		for itrk in range(ev.ntrk):
			#print ev.pid[itrk],ev.stat[itrk]
			if ev.pid[itrk]==PROTON: stat_pid_p=ev.stat[itrk]
			if ev.pid[itrk]==PIP:    stat_pid_pip=ev.stat[itrk]
		#print "stat_pid_p,stat_pid_pip=",stat_pid_p,stat_pid_pip
	
		#! Identify event to be in particular SS band
		#! + Note how ssb2 is directly set to true since it requires no cut on gpart_pid or stat_pid,
		#!   in other words, events are always a part of ssb2
		ssb1,ssb2,ssb3,ssb4=False,True,False,False
		if stat_pid_p>0 and stat_pid_pip>0: ssb1=True
		if (gpart_pid==3 or gpart_pid==4) and (stat_pid_p>0 and stat_pid_pip>0): ssb3=True
		if (gpart_pid==3 or gpart_pid==4): ssb4=True

		if ssb1: hmm[idtyp][SSB1][iwbin].Fill(ev.mmppip)
		if ssb2: hmm[idtyp][SSB2][iwbin].Fill(ev.mmppip)#! Note that ssb2 is always True
		if ssb3: hmm[idtyp][SSB3][iwbin].Fill(ev.mmppip)
		if ssb4: hmm[idtyp][SSB4][iwbin].Fill(ev.mmppip)

		#! Fill ANAHISTS 
		#! Make additional cut on MM here since the following needs to be studied only in signal region
		if ev.mm2ppip>MM2_CUT_EI_MIN and ev.mm2ppip<MM2_CUT_EI_MAX:
			hana[idtyp][GPART_PID_STAT_OFF][iwbin].Fill(gpart_pid)
			if stat_pid_p>0 and stat_pid_pip>0:
				hana[idtyp][GPART_PID_STAT_ON][iwbin].Fill(gpart_pid)
			hana[idtyp][NTRK_PID][iwbin].Fill(ntrk_pid)
			hana[idtyp][STAT_PID_P][iwbin].Fill(stat_pid_p)
			hana[idtyp][STAT_PID_PIP][iwbin].Fill(stat_pid_pip)
		#print "******"

#! ***

#! Plot
#! First set up aesthetics
#! Setup aesthetics for hmm
MRKR=[ROOT.gROOT.ProcessLine("kPlus"),ROOT.gROOT.ProcessLine("kCircle")]
for i in [ER,SR]:
	for j in range(NSSBANDS):
		for k in range(NWBINS):
			hmm[i][j][k].SetLineColor(j+1)
			hmm[i][j][k].SetMarkerColor(j+1)
			hmm[i][j][k].SetMarkerStyle(MRKR[i])
#! Setup aesthetics for hana
for i in [ER,SR]:
        for j in range(NANAHISTS):
                for k in range(NWBINS):
                        hana[i][j][k].SetLineColor(i+1)
                        hana[i][j][k].SetMarkerColor(i+1)

#! Draw ER*{NSSBANDS} and  SR*{NSSBANDS} on separate canvases
#! to bring about difference in these cuts on ER and SR 
#! output dir in .root 
FOUT.mkdir("comp_ER-SR_SSbands-cuts").cd()
#! anvases
c=[[0 for j in range(NWBINS)] for i in range(NDTYP)]
l=[[0 for j in range(NWBINS)] for i in range(NDTYP)]
for i in [ER,SR]:
	for j in range(NWBINS):
		cname="c%s_wbin%d"%(DTYP_NAME[i],j+1)
		c[i][j]=ROOT.TCanvas(cname,cname)
		l[i][j]=ROOT.TLegend(0.6,0.3,1.0,0.4)#,"","NDC");
		for iplot,k in enumerate([SSB2,SSB1,SSB3,SSB4]):#! Note order of drawing changed to draw max yield(=off*off) first #range(NSSBANDS):
			draw_opt=""
			if iplot>0:draw_opt="same"
			hmm[i][k][j].SetBit(ROOT.gROOT.ProcessLine("TH1::kNoTitle")) #! Turn of drawing of printing of hist title to canvas since legend contains the titles
			hmm[i][k][j].Draw(draw_opt)
			l[i][j].AddEntry(hmm[i][k][j],hmm[i][k][j].GetTitle())
		l[i][j].Draw()
		#! Draw cut lines
		MM_T2_CUTL_EI.SetY1(hmm[i][SSB1][j].GetMinimum())
		MM_T2_CUTL_EI.SetY2(hmm[i][SSB1][j].GetMaximum())
		MM_T2_CUTH_EI.SetY1(hmm[i][SSB1][j].GetMinimum())
		MM_T2_CUTH_EI.SetY2(hmm[i][SSB1][j].GetMaximum())
		MM_T2_CUTL_EI.Draw("same")
		MM_T2_CUTH_EI.Draw("same")
		#! Save canvas
		c[i][j].Write(c[i][j].GetName(),ROOT.gROOT.ProcessLine("TObject::kOverwrite"))#! for .root; latest namecycle
		c[i][j].Close()

#! Draw SSB1*{ER,SR},...,SSB4*{ER,SR} on separate canvases
#! to show that for all cuts there is no BG in the signal.
#! output dir in .root
FOUT.mkdir("background_SSbands-cut").cd()
#! canvases
cc=[[0 for j in range(NWBINS)] for i in range(NSSBANDS)]
ll=[[0 for j in range(NWBINS)] for i in range(NSSBANDS)]
#l=[0 for i in range(NDTYP)]
#cc=ROOT.TCanvas("cc","cc")
for i in range(NWBINS):
       	for j in range(NSSBANDS):
		cname="c%s_wbin%d"%(SSB_NAME[j],i+1)
		cc[j][i]=ROOT.TCanvas(cname,cname)
		#! 1. Draw ER first
		hmm[ER][j][i].SetBit(ROOT.gROOT.ProcessLine("TH1::kNoTitle")) #! Turn of drawing of printing of hist title to canvas since legend contains the titles
               	hmm[ER][j][i].Draw()
		#! 2. Normalize SR's max to ER's max and draw on same
		max_ER=hmm[ER][j][i].GetMaximum()
		if max_ER==0:max_ER=1
		max_SR=hmm[SR][j][i].GetMaximum()
		if max_SR==0:max_SR=1
               	scl_fctr_SR=max_ER/max_SR
               	if scl_fctr_SR==0:scl_fctr_SR=1
               	hmm[SR][j][i].Scale(scl_fctr_SR)
		hmm[ER][j][i].SetBit(ROOT.gROOT.ProcessLine("TH1::kNoTitle")) #! Turn of drawing of printing of hist title to canvas since legend contains the titles
               	hmm[SR][j][i].Draw("same")
       		#! Draw cut lines
       		MM_T2_CUTL_EI.SetY1(hmm[ER][j][i].GetMinimum())
       		MM_T2_CUTL_EI.SetY2(hmm[ER][j][i].GetMaximum())
       		MM_T2_CUTH_EI.SetY1(hmm[ER][j][i].GetMinimum())
       		MM_T2_CUTH_EI.SetY2(hmm[ER][j][i].GetMaximum())
       		MM_T2_CUTL_EI.Draw("same")
       		MM_T2_CUTH_EI.Draw("same")
		#! setup legend and draw
		ll[j][i]=ROOT.TLegend(0.6,0.3,1.0,0.4)
		ll[j][i].AddEntry(hmm[ER][j][i],hmm[ER][j][i].GetTitle())
		ll[j][i].AddEntry(hmm[SR][j][i],hmm[SR][j][i].GetTitle())
		ll[j][i].Draw()
		#! Draw and save canvas
		cc[j][i].Write(cc[j][i].GetName(),ROOT.gROOT.ProcessLine("TObject::kOverwrite"))#! for .root; latest namecycle
		cc[j][i].Close()

#! Draw gpart*{ER,SR},stat_p/pip*{ER,SR} on separate canvases
#! to compare their distributions in ER and SR
#! output dir in .root
FOUT.mkdir("ana_gpart-pid_stat-pid").cd()
#! canvases
c3=[[0 for j in range(NWBINS)] for i in range(NANAHISTS)]
l3=[[0 for j in range(NWBINS)] for i in range(NANAHISTS)]
#l=[0 for i in range(NDTYP)]
#cc=ROOT.TCanvas("cc","cc")
for i in range(NWBINS):
	for j in range(NANAHISTS):
		cname="c%s_wbin%d"%(ANAHISTS_NAME[j],i+1)
        	c3[j][i]=ROOT.TCanvas(cname,cname)
        	#! 1. Draw ER first
        	hana[ER][j][i].SetBit(ROOT.gROOT.ProcessLine("TH1::kNoTitle")) #! Turn of drawing of printing of hist title to canvas since legend contains the titles
        	hana[ER][j][i].Draw("hist")
        	#! 2. Normalize SR's max to ER's max and draw on same
		if j==STAT_PID_P or j==STAT_PID_PIP: #! then use integral to normalize since there are two peaks
			#print "j==STAT_PID_P or j==STAT_PID_PIP => using integral to normalize"
			max_ER=hana[ER][j][i].Integral()
                        if max_ER==0:max_ER=1
                        max_SR=hana[SR][j][i].Integral()
                        if max_SR==0:max_SR=1
                        scl_fctr_SR=max_ER/max_SR
                        if scl_fctr_SR==0:scl_fctr_SR=1
		else:
        		max_ER=hana[ER][j][i].GetMaximum()
        		if max_ER==0:max_ER=1
        		max_SR=hana[SR][j][i].GetMaximum()
        		if max_SR==0:max_SR=1
        		scl_fctr_SR=max_ER/max_SR
        		if scl_fctr_SR==0:scl_fctr_SR=1
        	hana[SR][j][i].Scale(scl_fctr_SR)
        	hana[ER][j][i].SetBit(ROOT.gROOT.ProcessLine("TH1::kNoTitle")) #! Turn of drawing of printing of hist title to canvas since legend contains the titles
        	hana[SR][j][i].Draw("hist same")
        	#! setup legend and draw
        	l3[j][i]=ROOT.TLegend(0.6,0.3,1.0,0.4)
        	l3[j][i].AddEntry(hana[ER][j][i],hana[ER][j][i].GetTitle())
        	l3[j][i].AddEntry(hana[SR][j][i],hana[SR][j][i].GetTitle())
        	l3[j][i].Draw()
		#! Draw TPaveText to say that top=2 if anahist=gpart
		if j==GPART_PID_STAT_ON or j==GPART_PID_STAT_OFF:
			pt=ROOT.TPaveText(.55,.7,.95,.8,"NDC")
  			pt.AddText("top=2 ( exclusive e',p,#pi^{+} event i.e. no #pi^{-} detected)")
  			#pt.SetTextSize(0.42)
  			pt.Draw()
        	#! Draw and save canvas
        	c3[j][i].Write(c3[j][i].GetName(),ROOT.gROOT.ProcessLine("TObject::kOverwrite"))#! for .root; latest namecycle
        	c3[j][i].Close()

#! Close FOUT
FOUT.Close()

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
	
