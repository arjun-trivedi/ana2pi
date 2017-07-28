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
+ Usage: study_3piPSsim.py debug=False
'''

USAGE='study_3piPSsim.py debug=False test_fit_only=False'

#! Get input data from user
DBG=False
if len(sys.argv)>1: #i.e. debug entered by user
	if    sys.argv[1]=="True":  DBG=True
	elif  sys.argv[1]=="False": DBG=False
	else: sys.exit('DBG=%s is not valid. usage: %s'%(DBG,USAGE))

#! *** Prepare OUTDIR***
DATE=datetime.datetime.now().strftime('%m%d%y')
OUTDIR=os.path.join(os.environ['STUDY_EVTSEL_BG_DATADIR'],"evtsel_study_3piPSsim_%s"%DATE) # os.environ['STUDY_EVTSEL_BG_DATADIR']
if DBG==True: OUTDIR=os.path.join(os.environ['STUDY_EVTSEL_BG_DATADIR'],"evtsel_study_3piPSsim_dbg_%s"%DATE) #os.environ['STUDY_EVTSEL_BG_DATADIR']
print "OUTDIR=",OUTDIR
if not os.path.exists(OUTDIR):
	os.makedirs(OUTDIR)
#sys.exit()
#! ***

#! *** Prepare NENTRIES ***
NENTRIES=1000000000
#if DBG==True: NENTRIES=1000

#! *** Now prepare INDIR and input data structures: FIN/T[dtyp] ***
INDIR=os.path.join(os.environ['D2PIDIR_SIM_E16'],"data_evtsel_BG_020117")

print "INDIR=",INDIR
#sys.exit()

NDTYP=2
DTYP_NAME=['SR','ST']
SR,ST=range(NDTYP)
FIN,T=[0 for i in range(NDTYP)],[0 for i in range(NDTYP)]
FIN[SR]=ROOT.TFile("%s/d3piPS.root"%INDIR)
FIN[ST]=ROOT.TFile("%s/d3piPS_ST.root"%INDIR)
T[SR]=FIN[SR].Get("d2piR/tR")
T[ST]=FIN[ST].Get("d2piT/tT")
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
	Q2BIN_LEL=[2.00]
	Q2BIN_UEL=[5.00]
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
#print WBIN_LEL
#sys.exit()
#print Q2BIN_LEL
#print WBIN_LEL

#! DRAW_CMD
NPLTS=2
MM,MM2=range(NPLTS)
PLT_NAME=['MM','MM2']
PLT_TITLE=['MM','MM^{2}']
PLT_UNIT=['GeV','GeV^{2}']

DRAW_CMD=[0 for i in range(NPLTS)]
DRAW_CMD[MM] ="mmppip>>hmmcmd(200,-0.5,1)"
DRAW_CMD[MM2]="mm2ppip>>hmmcmd(100,-0.2,0.2)"

#! Create hmm[q2][w][plt][dtyp]
hmm=[[[[0 for l in range(NDTYP)] for k in range(NPLTS)] for j in range(len(WBIN_LEL))] for i in range(len(Q2BIN_LEL))]
#print "hmm=",hmm
#sys.exit()

#! *** Setup some analysis constants ***

#! Fill hmm
for iq,q2bin_le in enumerate(Q2BIN_LEL):
	if DBG==True and iq>0: continue #! debug
	q2min=q2bin_le
	q2max=Q2BIN_UEL[iq]
	for iw,wbin_le in enumerate(WBIN_LEL):
		#if iw==0 or iw==1
		#! DBGW (3)
		#if DBG==True and (iw!=0 and iw!=15 and iw!=28): continue #! debug #! (iw!=0 and iw!=15 and iw!=28)
		#if DBG==True and iw+1 > 10: continue
		#if DBG==True and (iw+1!=5 and iw+1!=6 and iw+1!=7): continue
		wmin=wbin_le
		wmax=wbin_le+WBINW
		cut_q2w=ROOT.TCut("Q2>%f && Q2<%f && W>%f && W<%f"%(q2min,q2max,wmin,wmax))
		print "q2w=",cut_q2w.GetTitle()
		for iplt in [MM,MM2]:
			for idtyp in range(NDTYP):
				cut=ROOT.TCut(cut_q2w)
				if idtyp==SR:
					cut_top=ROOT.TCut("top==2")
					cut+=cut_top
				print "TTree::Draw() for",PLT_NAME[iplt],DTYP_NAME[idtyp]
				#print T[idtyp].GetName()
				T[idtyp].Draw(DRAW_CMD[iplt],cut,"",NENTRIES)#[idtyp][isctr][icorr],cut)

				#! Store histogram
				#ROOT.gStyle.SetOptStat("ne")
				htmp=ROOT.gDirectory.Get("hmmcmd")#(draw_cmd_hst_name[idtyp][isctr][icorr])
				hmm[iq][iw][iplt][idtyp]=htmp.Clone()#"hmm_%s_%s"%(DTYP_NAME[idtyp],CORR_NAME[icorr]))
				hmm[iq][iw][iplt][idtyp].SetName("h%s_%s"%(PLT_NAME[iplt],DTYP_NAME[idtyp]))
				title="%s-%s for Q^{2}=[%.2f,%.2f) GeV^{2}, W=[%.3f,%.3f) GeV"%(PLT_TITLE[iplt],DTYP_NAME[iplt],q2min,q2max,wmin,wmax)
				hmm[iq][iw][iplt][idtyp].SetTitle(title)#"%s %.2f-%.2f_%.3f-%.3f"%(PLT_TITLE[iplt],q2min,q2max,wmin,wmax))#("%s"%cut.GetTitle())
				#print hmm[idtyp][i][j][icorr].GetName()
#sys.exit()
#! + Plots
ROOT.gStyle.SetOptStat(0) #"n"
ROOT.gStyle.SetOptFit(1)#111)
#! Aesthetic related objects
CLRS_DTYP=[ROOT.gROOT.ProcessLine("kBlue"),ROOT.gROOT.ProcessLine("kRed"),ROOT.gROOT.ProcessLine("kGreen")]

LE=[0 for i in range(NPLTS)]
#LE[MM]=[0.278,0.417]
LE[MM]=[0.300,0.400]
LE[MM2]=[LE[MM][0]**2,LE[MM][1]**2]

#! TLine showing fit limits
LE_LINE_MIN=[0 for i in range(NPLTS)]
LE_LINE_MAX=[0 for i in range(NPLTS)]
LE_LINE_MIN[MM]=ROOT.TLine( LE[MM][0], 0, LE[MM][0], 0)
LE_LINE_MAX[MM]=ROOT.TLine( LE[MM][1], 0, LE[MM][1], 0)
LE_LINE_MIN[MM2]=ROOT.TLine( LE[MM2][0], 0, LE[MM2][0], 0)
LE_LINE_MAX[MM2]=ROOT.TLine( LE[MM2][1], 0, LE[MM2][1], 0)
for line in [LE_LINE_MIN[MM],LE_LINE_MAX[MM],LE_LINE_MIN[MM2],LE_LINE_MAX[MM2]]:
	line.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
	line.SetLineWidth(3)

#! TLine depicting pion mass debug
MASS_PION=0.13957018
MASS_PROTON=0.938
MASS_PION_LINE=ROOT.TLine(MASS_PION,0,MASS_PION,0)
MASS_PION_LINE.SetLineColor(ROOT.gROOT.ProcessLine("kGreen"))
MASS_PION_LINE.SetLineWidth(3)

MM_THRSHLD_2PI=ROOT.TLine(2*MASS_PION,0,2*MASS_PION,0)
MM_THRSHLD_2PI.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
MM_THRSHLD_2PI.SetLineWidth(3)

MM2_THRSHLD_2PI=ROOT.TLine((2*MASS_PION)**2,0,(2*MASS_PION)**2,0)
MM2_THRSHLD_2PI.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
MM2_THRSHLD_2PI.SetLineWidth(3)

for iq,q2bin_le in enumerate(Q2BIN_LEL):
	if DBG==True and iq>0: continue #! debug
	q2min=q2bin_le
	q2max=Q2BIN_UEL[iq]
	outdir_q2w="%s/%.2f-%.2f"%(OUTDIR,q2min,q2max)
	if not os.path.exists(outdir_q2w):
		os.makedirs(outdir_q2w)

	for iplt in range(NPLTS):
		for idtyp in range(NDTYP):
			outdir_q2w_plt_dtyp="%s/%s/%s"%(outdir_q2w,PLT_NAME[iplt],DTYP_NAME[idtyp])
			if not os.path.exists(outdir_q2w_plt_dtyp):
				os.makedirs(outdir_q2w_plt_dtyp)

			for iw,wbin_le in enumerate(WBIN_LEL):
				#! DBGW (3)
				#if DBG==True and (iw!=0 and iw!=15 and iw!=28): continue #! debug #!(iw!=0 and iw!=15 and iw!=28)
				#if DBG==True and iw+1 > 10: continue
				#if DBG==True and (iw+1!=2 and iw+1!=3 and iw+1!=4 and iw+1!=5): continue
				wmin=wbin_le
				wmax=wbin_le+WBINW

				cname="c%s_%.3f-%.3f"%(PLT_NAME[iplt],wmin,wmax)
				cmm=ROOT.TCanvas(cname,cname)
				hmm[iq][iw][iplt][idtyp].Draw()
				#! Draw fit limits
				if iplt==MM:
					lmin,lmax=LE_LINE_MIN[MM],LE_LINE_MAX[MM]			
				elif iplt==MM2:
					lmin,lmax=LE_LINE_MIN[MM2],LE_LINE_MAX[MM2]
				for line in [lmin,lmax]:
					line.SetY1(hmm[iq][iw][iplt][idtyp].GetMinimum())
					line.SetY2(hmm[iq][iw][iplt][idtyp].GetMaximum())
					line.Draw("same")

				cmm.SaveAs("%s/%s.png"%(outdir_q2w_plt_dtyp,cname))

	
	
				

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
	
