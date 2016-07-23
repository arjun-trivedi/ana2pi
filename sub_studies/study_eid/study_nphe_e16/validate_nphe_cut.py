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

import atlib as attools

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

#! *** Prepare constants 
NDTYP=2
ER,SR=range(NDTYP)
DTYP_NAME=["ER","SR"]

NCUTS=2
NNPHE,WNPHE=range(NDTYP)
CUT_NAME=["nnphe","wnphe"]

NSCTRS=6

NSGMNTS=18

NPMTS=2
L,R=range(NPMTS)

#! *** Prepare input data structure: FIN/T[dtyp][cut] ***
FIN=[[0 for j in range(NCUTS)] for i in range(NDTYP)]
T  =[[0 for j in range(NCUTS)] for i in range(NDTYP)]
#print "FIN=",FIN
#print "T=",T
#! ***

#! *** Prepare output datadir***
OUTDIR=os.environ['STUDY_EID_NPHE_E16_DATADIR']
if DBG==True: OUTDIR+="/dbg"
print "OUTDIR=",OUTDIR
if not os.path.exists(OUTDIR):
	os.makedirs(OUTDIR)
#sys.exit()
#! ***

#! *** Prepare NENTRIES ***
NENTRIES=1000000000
if DBG==True: NENTRIES=1000

#! *** Now fill FIN/T[dtyp][cut] ***
INPUT=[0 for i in range(NDTYP)]
INPUT[ER]=os.path.join(os.environ['D2PIDIR_EXP_E16'],"data_nphe_071916")#! debug /bk
INPUT[SR]=os.path.join(os.environ['D2PIDIR_SIM_E16'],"data_nphe_071916")#! debug /bl
#print "INPUT[ER]=",INPUT[ER]
#print "INPUT[SR]=",INPUT[SR]

FIN[ER][NNPHE]=ROOT.TFile("%s/dnphe_nnphe.root"%(INPUT[ER]))
FIN[ER][WNPHE]=ROOT.TFile("%s/dnphe_wnphe.root"%(INPUT[ER]))
FIN[SR][NNPHE]=ROOT.TFile("%s/dnphe_nnphe.root"%(INPUT[SR]))
FIN[SR][WNPHE]=None #! does not exist since nphe cut not applied for SR
T[ER][NNPHE]=FIN[ER][NNPHE].Get("d2piR/tR")
T[ER][WNPHE]=FIN[ER][WNPHE].Get("d2piR/tR")
T[SR][NNPHE]=FIN[SR][NNPHE].Get("d2piR/tR")
#print "FIN[ER][NNPHE]=",FIN[ER][NNPHE].GetName()
#print "FIN[ER][WNPHE]=",FIN[ER][WNPHE].GetName()
#print "FIN[SR][NNPHE]=",FIN[SR][NNPHE].GetName()
#print "T[ER][NNPHE]=",T[ER][NNPHE]
#print "T[ER][WNPHE]=",T[ER][WNPHE]
#print "T[SR][NNPHE]=",T[SR][NNPHE]
#sys.exit()
#! ***

#! *** Now make plots ***

#! Define some constants
#! Q2-W binning
Q2BIN_LEL=[2.00,2.40,3.00,3.50,4.20,2.00]#[1.50,2.00,2.40,3.00,3.50,4.20] #debug
Q2BIN_UEL=[2.40,3.00,3.50,4.20,5.00,5.00]#[2.00,2.40,3.00,3.50,4.20,5.00] #debug
NQ2BINS=len(Q2BIN_LEL)
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
def get_bin_idx(kin,val):
        ''' Used to get bin's index from a given kin(=Q2 or W) value'''
	if kin=="Q2":
		nbins=NQ2BINS
		lel,uel=Q2BIN_LEL,Q2BIN_UEL
	elif kin=="W":
		nbins=NWBINS
                lel,uel=WBIN_LEL,WBIN_UEL
	else:
		sys.exit("kin = %s not recognized. kin = Q2 or W only. "%kin)

        found_bin=False
        for ibin in range(nbins):
                min=lel[ibin]
                max=uel[ibin]
                if val>=min and val<max:
                        found_bin=True
                        return ibin
        if found_bin==False:
                print "bin idx not found. Going to return an invalid index of -9999"
                return -9999

def get_wbin_idx(W):
	''' Used to get W bin's index from a given W value'''	
	return get_bin_idx("W",W)
def get_q2bin_idx(Q2):
        ''' Used to get Q2 bin's index from a given Q2 value'''   
        return get_bin_idx("Q2",Q2)

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

#! Create histograms to be made
#! + hmm/hnphe[dtyp][cut][q2][w]
#! + For now use hnphe[dtyp][cut][q2][w] to monitor cut even though direct verification is done with hnphe [ER][cut][sct][sgm][pmt]
hmm=[[[[ROOT.TH1F("hmm%s_%s_q2bin%d_wbin%d"%(DTYP_NAME[i],CUT_NAME[j],k+1,l+1),"mmppip_%s_%s_q2bin_%.2f-%.2f_wbin_%.3f-%.3f"%(DTYP_NAME[i],CUT_NAME[j],Q2BIN_LEL[k],Q2BIN_UEL[k],WBIN_LEL[l],WBIN_UEL[l]),200,-0.5,1) for l in range(NWBINS)]for k in range(NQ2BINS)]for j in range(NCUTS)] for i in range(NDTYP)]
hnphe=[[[[ROOT.TH1F("hnphe%s_%s_q2bin%d_wbin%d"%(DTYP_NAME[i],CUT_NAME[j],k+1,l+1),"nphe_%s_%s_q2bin_%.2f-%.2f_wbin_%.3f-%.3f"%(DTYP_NAME[i],CUT_NAME[j],Q2BIN_LEL[k],Q2BIN_UEL[k],WBIN_LEL[l],WBIN_UEL[l]),100,0,300) for l in range(NWBINS)]for k in range(NQ2BINS)]for j in range(NCUTS)] for i in range(NDTYP)]
#print hmm
#print hnphe
#sys.exit()

#! Process of looping over TTree in Python taken from Evan's mail of [02-7-15] 'pyroot and bytes in trees'
d=list(itertools.product(*[[ER,SR],[NNPHE,WNPHE]]))
for r in d:
	idtyp,icut=r[0],r[1]
	if idtyp==SR and icut==WNPHE: continue #! a TTree does not exist
	for i, ev in enumerate(T[idtyp][icut], 1):
		if (i%10000==0): print "processing TTree for %s,%s: ientry# %d"%(DTYP_NAME[idtyp],CUT_NAME[icut],i)
		if DBG and i>1000: break
	
		#print "***entry #%d ***"%(i+1)

		#! Cuts
		#! 1. top cut: Use only top2 events even though top2' is actually top2+top1
		#! + This is because in top1, as cacluated by 'hvy', mmppip=0 
		#print "ev top=",ev.top
		if ev.top!=2: continue # and ev.top!=1: continue

		#! 2. Q2, W cut: to keep Q2 & W withing analysis range :[2.0,5.0], [1.400,2.125)
		if ev.Q2<2.0 or ev.Q2>=5.0: continue
		if ev.W<1.400 or ev.W>=2.125: continue

		#! Get iq2bin,iwbin
		#print "ev.Q2,ev.W=",ev.Q2,ev.Q2
		iq2bin=get_q2bin_idx(ev.Q2)
		iwbin=get_wbin_idx(ev.W)	
		#print "iq2bin,iwbin=",iq2bin,iwbin
		#print "passed top=",ev.top

		#! Fill histograms
		hmm[idtyp][icut][iq2bin][iwbin].Fill(ev.mmppip)
		hnphe[idtyp][icut][iq2bin][iwbin].Fill(ev.nphe)
	

#! ***
#sys.exit()

#! Plot and save histograms
#! First set up aesthetics
#! Setup aesthetics for hmm
MRKR=[ROOT.gROOT.ProcessLine("kPlus"),ROOT.gROOT.ProcessLine("kCircle")]
CLR=[ROOT.gROOT.ProcessLine("kBlue"),ROOT.gROOT.ProcessLine("kRed")]
d=list(itertools.product(*[[ER,SR],[NNPHE,WNPHE],range(NQ2BINS),range(NWBINS)]))
for r in d:
	idtyp,icut,iq2,iw=r[0],r[1],r[2],r[3]
	hmm[idtyp][icut][iq2][iw].SetLineColor(CLR[idtyp])
	hmm[idtyp][icut][iq2][iw].SetMarkerColor(CLR[idtyp])
	hmm[idtyp][icut][iq2][iw].SetMarkerStyle(MRKR[idtyp])
#sys.exit()

#! create output .root file and save hists in folders:
#! + prec/q2/w
#! + pstc/q2/w
FOUT=ROOT.TFile("%s/validate_nphe_cut.root"%OUTDIR,"RECREATE")	
CWDTH=500#2000#500
CHGHT=300#1500#300
for icut in range(NCUTS):
	cutdir=FOUT.mkdir(CUT_NAME[icut])
	#cutdir.cd()
	for iq2 in range(NQ2BINS):
		q2min,q2max=Q2BIN_LEL[iq2],Q2BIN_UEL[iq2]
		q2bin="%.1f-%.1f"%(q2min,q2max)
		q2bindir=cutdir.mkdir(q2bin)
		#q2bindir.cd()
		for iw in range(NWBINS):
                	wmin,wmax=WBIN_LEL[iw],WBIN_UEL[iw]
                	wbin="%.3f-%.3f"%(wmin,wmax)
                	q2bindir.mkdir(wbin).cd()	
			#! Draw plos on canvas and save
			#! nphe
			cname="hnphe_%d_%d_%d"%(icut+1,iq2+1,iw+1)
                        cnphe=ROOT.TCanvas(cname,cname,CWDTH,CHGHT)
                        hnphe[ER][icut][iq2][iw].Draw()
			cnphe.Write("%s"%cname,ROOT.gROOT.ProcessLine("TObject::kOverwrite"))#! for .root; latest namecycle
			#! mm
			#! + Note WNPHE does not exist for SR!
			cname="hmm_%d_%d_%d"%(icut+1,iq2+1,iw+1)
			cmm=ROOT.TCanvas(cname,cname,CWDTH,CHGHT)
			attools.scale_hist(hmm[ER][icut][iq2][iw],hmm[SR][NNPHE][iq2][iw])
			hmm[ER][icut][iq2][iw].Draw()
			hmm[SR][NNPHE][iq2][iw].Draw("sames")
			#hmm[ER][j][i].SetBit(ROOT.gROOT.ProcessLine("TH1::kNoTitle")) #! Turn of drawing of printing of hist title to canvas since legend contains the titles
			#! Setup ylim of cut lines and draw them
			attools.setup_cut_lines_ylim(MM_T2_CUTL_EI,MM_T2_CUTH_EI,hmm[ER][icut][iq2][iw])
			MM_T2_CUTL_EI.Draw("same")
			MM_T2_CUTH_EI.Draw("same")
			cmm.Write("%s"%cname,ROOT.gROOT.ProcessLine("TObject::kOverwrite"))#! for .root; latest namecycle
#sys.exit()	
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
	
