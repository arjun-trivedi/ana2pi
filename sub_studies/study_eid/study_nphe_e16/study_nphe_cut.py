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
NSCTRS=6

NSGMNTS=18

#! setup PMT related constants
#! NOTE that the following structure is carried on to h10looper_e1f(2pi).h/.cpp
#          L  C  R
#pmt-h10  -1  0  1
#pmt       1  2  3
#ipmt      0  1  2
NPMTS=3
L,C,R=range(NPMTS)

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
INPUT=os.path.join(os.environ['D2PIDIR_EXP_E16'],"data_nphe_071916")#! debug /bk
#print "INPUT=",INPUT
FIN=ROOT.TFile("%s/dnphe_nnphe.root"%(INPUT))
print "FIN=",FIN.GetName()
T=FIN.Get("d2piR/tR")
print "T=",T.GetName()
#sys.exit()
#! ***

#! *** Setup objects that will be created
#! 1. hnphe[sct][sgm][pmt]

#! + xmax for hnphe: usually 300, but can be higher for some pmts as noted below 
#! Init to 300
XMAX=[[[300 for k in range(NPMTS)]for j in range(NSGMNTS)]for i in range(NSCTRS)]
#! Now make special adjustments
d=itertools.product(*[range(NSCTRS),range(NSGMNTS),range(NPMTS)])
for r in d:
    isct,isgm,ipmt=r[0],r[1],r[2]
    sct,sgm,pmt=isct+1,isgm+1,ipmt+1

    #! For sgms=12,13,14,15 xmax appears higher
    if (sct==5 or sct==6) and (sgm==12 or sgm==13 or sgm==14 or sgm==15):
        XMAX[isct][isgm][ipmt]=500

    #! For sct:sgm:pmt=6:10:2 and 6:11:2 xmax is higher
    if sct==6 and sgm==10 and pmt==2:
	XMAX[isct][isgm][ipmt]=500
    if sct==6 and sgm==11 and pmt==2:
        XMAX[isct][isgm][ipmt]=500

#! Now finally create hnphe
hnphe=[[[ROOT.TH1F("hnphe_%d_%d_%d"%(i+1,j+1,k+1),"hnphe_sct%d_sgm%d_pmt%d"%(i+1,j+1,k+1),100,0,XMAX[i][j][k]) for k in range(NPMTS)]for j in range(NSGMNTS)]for i in range(NSCTRS)]
#print hnphe
#sys.exit()
#! *** Now make plots ***

#! Process of looping over TTree in Python taken from Evan's mail of [02-7-15] 'pyroot and bytes in trees'
for i, ev in enumerate(T, 1):
	if (i%10000==0): print "processing ientry# %d"%(i)
	if DBG and i>1000: break
	
	isct=ev.sector-1
        isgm=ev.cc_segm-1
	if   ev.pmt==-1: ipmt=0
	elif ev.pmt==0: ipmt=1
        elif ev.pmt==1: ipmt=2
        #else: continue #! if pmt==0, then eff=1
	sct=isct+1
	sgm=isgm+1

	#! debug
	#if (sct!=1 and sgm!=2 and pmt!=2) or (sct!=5 and sgm!=15 and pmt!=2) or (sct!=6 and sgm!=17 and pmt!=1):continue

	#! Cuts
	#! 1. top cut: Use only top2 events even though top2' is actually top2+top1
	#! + This is because in top1, as cacluated by 'hvy', mmppip=0
	if ev.top!=2: continue

	#! 2. Q2, W cut: to keep Q2 & W withing analysis range :[2.0,5.0], [1.400,2.125)
	#! Removed from here and made dependent on cc_segm. See below.
	#if ev.Q2<2.0 or ev.Q2>=5.0: continue
	#if ev.W<1.400 or ev.W>=2.125: continue

	#! 3. MM cut => moved into segment loop for manipulating statistics
	#if ev.mm2ppip<-0.04 or ev.mm2ppip>0.06: continue

	#! + Make Q2-W  && MM cuts only if segment!=1,16,17 or 18 because:
	#!   + I observed that with ana-Q2-W cut, all other segments had data but sgms=1,16,17,18 and
	#!     therefore I tried to remove the Q2-W wondering if that was the cause.
	#!   + However, it could well be that these segments generally have very less experimentally data, especially 1 and 18,
	#!     which may be just out of the kinematic reach of E16
	if sgm==1 or sgm==15 or sgm==16 or sgm==17 or sgm==18: #! do not apply Q2,W cut because it is because of that there is no data for these?
		#! Fill histograms
		hnphe[isct][isgm][ipmt].Fill(ev.nphe)
	else:
		#! MM cut
		if ev.mm2ppip<-0.04 or ev.mm2ppip>0.06: continue
		#! Q2, W cut: to keep Q2 & W withing analysis range :[2.0,5.0], [1.400,2.125)
		if ev.Q2<2.0 or ev.Q2>=5.0: continue
		if ev.W<1.400 or ev.W>=2.125: continue
		#! Fill histograms
		hnphe[isct][isgm][ipmt].Fill(ev.nphe)
	

#! ***
#sys.exit()

#! Fit histograms
#! Define nphe cut related objects: loose (lse) and tight (tgt) cuts. Loose cuts are conservative (allow greater
#! possibility for background) and tight are more aggresive (less possibility for background)
#! + NPHE_CUT_LSE: 30 for most, except for:
#!   + sct1, sgm>=15:40
#!   + sct2,3,4,5,6, sgm>=15:50
#! + NPHE_CUT_TGT=NPHE_CUT_LSE+20
#! Init to 30
NPHE_CUT_LSE=[[30 for j in range(NSGMNTS)] for i in range(NSCTRS)]
#! Now make adjustments for exceptions
for r in itertools.product(*[range(NSCTRS),range(NSGMNTS)]):
    isct,isgm=r[0],r[1]
    sct,sgm=isct+1,isgm+1
    if sgm>=15:
        if sct==1: 
            NPHE_CUT_LSE[isct][isgm]=40
        else:
            NPHE_CUT_LSE[isct][isgm]=50
#! Now setup NPHE_CUT_TGT
#! Init to 0
NPHE_CUT_TGT=[[0 for j in range(NSGMNTS)] for i in range(NSCTRS)]
#! Now set to NPHE_CUT_LSE+20
for r in itertools.product(*[range(NSCTRS),range(NSGMNTS)]):
    isct,isgm=r[0],r[1]
    NPHE_CUT_TGT[isct][isgm]=NPHE_CUT_LSE[isct][isgm]+20
#! Define CUTL_LSE/TGT as per NPHE_CUT_LSE/TGT
CUTL_LSE=[[ROOT.TLine(NPHE_CUT_LSE[i][j],0,NPHE_CUT_LSE[i][j],0) for j in range(NSGMNTS)] for i in range(NSCTRS)]
CUTL_TGT=[[ROOT.TLine(NPHE_CUT_TGT[i][j],0,NPHE_CUT_TGT[i][j],0) for j in range(NSGMNTS)] for i in range(NSCTRS)]
for r in itertools.product(*[range(NSCTRS),range(NSGMNTS)]):
	isct,isgm=r[0],r[1]
	CUTL_LSE[isct][isgm].SetLineColor(ROOT.gROOT.ProcessLine("kGreen"))
	CUTL_TGT[isct][isgm].SetLineColor(ROOT.gROOT.ProcessLine("kRed"))

#! Define function to return fit function
def get_ffnphe():
	'''
	+ fit function as obtained from YT (verified to same as KP)
	+ fit function init-pars as per YT
	'''
	fnphe=ROOT.TF1("fnphe", "[0]*TMath::Power(([1]),(x/[2]))*(TMath::Exp(-([1])))/TMath::Gamma((x/[2])+1)")#,50,300)
	fnphe.SetParName(0,"A")
	fnphe.SetParName(1,"L")
	fnphe.SetParName(2,"p")
	fnphe.SetParameters(58294,2.28921,38.0199)
	return fnphe
#! Define fit limits (=FITL[sct][sgm][pmt]) which is usually [NPHE_CUT_LSE,XMAX], except for some cases defined below
#! (NPHE_CUT_TGT arguably would be a better low limit for the fit since it moves past event from the noise-peak more. However, 
#!  fits seems to perform better when the lower limit was to LSE cuts. Anyway, I do not think the LSE/TGT low fit limits
#!  will make much difference)
#! 1. Init to [NPHE_CUT_LSE,XMAX]
FITL=[[[[NPHE_CUT_LSE[i][j],XMAX[i][j][k]] for k in range(NPMTS)]for j in range(NSGMNTS)]for i in range(NSCTRS)]
#! 2. Special cases
d=itertools.product(*[range(NSCTRS),range(NSGMNTS),range(NPMTS)])
for r in d:
    isct,isgm,ipmt=r[0],r[1],r[2]
    sct,sgm,pmt=isct+1,isgm+1,ipmt+1

    #! 1. FITL tuning for lower segments (sgm<5):
    #! +  For lower segments, with negligible noise-peak, sometimes the fit
    #!    needs to be constrained lower so as to go to 0 as x tends to 
    #! +  Sometimes purely for getting any reasonable fits the following need to done
    if sct==1 and sgm==2 and pmt==1:
        FITL[isct][isgm][ipmt]=[10,XMAX[isct][isgm][ipmt]]
    if sct==1 and sgm==2 and pmt==3:
        FITL[isct][isgm][ipmt]=[50,XMAX[isct][isgm][ipmt]]

    if sct==4 and sgm==2:# and pmt=1/3:
        FITL[isct][isgm][ipmt]=[10,XMAX[isct][isgm][ipmt]]

    if sct==5 and sgm==2 and pmt==1:
        FITL[isct][isgm][ipmt]=[20,XMAX[isct][isgm][ipmt]]
    if sct==5 and sgm==2 and pmt==3:
        FITL[isct][isgm][ipmt]=[30,XMAX[isct][isgm][ipmt]]
    

    #! 2. FITL tuning for higher segments (sgm>=15)
    #! 1. Fits that needs fine tuning mainly at lower and upper points of the distribution to address issues like:
    #!    + statistical fluctuations towards end of distribution
    #!    + simply to get a reasonable fit

    #! Generally is sgm>=15, start fit at higher nphe count
    if sgm==15 or sgm==16 or sgm==17 or sgm==18:
        FITL[isct][isgm][ipmt]=[60,XMAX[isct][isgm][ipmt]]

    #! Exceptions to the above generality
    if sct==2 and sgm==15 and pmt==3:
        FITL[isct][isgm][ipmt]=[80,XMAX[isct][isgm][ipmt]]
    
    if sct==3 and sgm==16 and pmt==3:
        FITL[isct][isgm][ipmt]=[80,XMAX[isct][isgm][ipmt]]

    if sct==4 and sgm==15 and pmt==1:
        FITL[isct][isgm][ipmt]=[80,XMAX[isct][isgm][ipmt]]
    if sct==4 and sgm==17 and pmt==1:
        FITL[isct][isgm][ipmt]=[40,XMAX[isct][isgm][ipmt]]

    if sct==5 and sgm==16 and pmt==1:
        FITL[isct][isgm][ipmt]=[100,XMAX[isct][isgm][ipmt]]
    if sct==5 and sgm==17 and pmt==1:
        FITL[isct][isgm][ipmt]=[100,XMAX[isct][isgm][ipmt]]

    if sct==6 and sgm==16 and pmt==3:
        FITL[isct][isgm][ipmt]=[100,XMAX[isct][isgm][ipmt]]

    #! 3. Fits that needs to be adjusted due to issues that I have not understood
    if sct==2 and sgm==11 and pmt==1:
        FITL[isct][isgm][ipmt]=[40,XMAX[isct][isgm][ipmt]]

    if sct==5 and sgm==12 and pmt==1:
	FITL[isct][isgm][ipmt]=[100,XMAX[isct][isgm][ipmt]]
    
#! Now fit
#! NOTE, do not fit pmt==2 i.e. center pmt
d=itertools.product(*[range(NSCTRS),range(NSGMNTS),range(NPMTS)])
for r in d:
	isct,isgm,ipmt=r[0],r[1],r[2]
	
	#! Skip fitting central pmt
	if ipmt==C: continue

	f=get_ffnphe()
	#! Draw histogram and fit to given range
	range_min=FITL[isct][isgm][ipmt][0]
	range_max=FITL[isct][isgm][ipmt][1]
	hnphe[isct][isgm][ipmt].Fit(f,"","",range_min,range_max)#NPHE_CUT_LSE,300)

#! Save histograms with fits, fitted to FITL and the FULLFIT (see below), in following format:
#! + sector/canvas_L/R(6*3)
CWDTH=500
CHGHT=500
#! hnphe YMAX for sgm:15,16,17
YMAX={15:500,16:500,17:500}
#1 FOUT
FOUT=ROOT.TFile("%s/study_nphe_cut.root"%OUTDIR,"RECREATE")
#! Objects to store and draw the the FULLFIT based in fits to FITL
FULLFITF=[[[0 for k in range(NPMTS)]for j in range(NSGMNTS)]for i in range(NSCTRS)]
for isct in range(NSCTRS):
	FOUT.mkdir("sector%d"%(isct+1)).cd()
	for ipmt in range(NPMTS):
		pmt=ipmt+1
		cname=("sct%d_pmt%d"%(isct+1,ipmt+1))
		c=ROOT.TCanvas(cname,cname,CWDTH,CHGHT)
		c.Divide(6,3)
		for isgm in range(NSGMNTS):
			sgm=isgm+1
			c.cd(isgm+1)
			#! Draw hist and associated fit function that is fitted to range specified in FITL
			#! Set ymax for sgm=15,16,17
			if sgm==15 or sgm==16 or sgm==17:
				 hnphe[isct][isgm][ipmt].SetMaximum(YMAX[sgm])	
			hnphe[isct][isgm][ipmt].Draw()
			#! FULLFIT: Draw a superimposed fit function that covers the full range to check the fits's integrity
			#! if pmt==2 then there is not fit
			if pmt==2: fitf=None
			else:      fitf=hnphe[isct][isgm][ipmt].GetFunction("fnphe")
        		if fitf!=None:
                		FULLFITF[isct][isgm][ipmt]=fitf.Clone("full_fit_sct%d_sgm%d_pmt%d"%(isct+1,isgm+1,ipmt+1))
                		FULLFITF[isct][isgm][ipmt].SetRange(0,XMAX[isct][isgm][ipmt])
                		FULLFITF[isct][isgm][ipmt].SetLineColor(ROOT.gROOT.ProcessLine("kGreen"))
                		FULLFITF[isct][isgm][ipmt].SetLineStyle(3)
				FULLFITF[isct][isgm][ipmt].SetLineWidth(3)
                		FULLFITF[isct][isgm][ipmt].Draw("same")
			#! Setup and draw cut line
			attools.setup_cut_lines_ylim(CUTL_LSE[isct][isgm],None,hnphe[isct][isgm][ipmt])
			attools.setup_cut_lines_ylim(CUTL_TGT[isct][isgm],None,hnphe[isct][isgm][ipmt])
			CUTL_LSE[isct][isgm].Draw("same")
			CUTL_TGT[isct][isgm].Draw("same")
		c.Write("%s"%cname,ROOT.gROOT.ProcessLine("TObject::kOverwrite"))#! for .root; latest namecycle
	
FOUT.Close()	

#! Obtain EFF_LSE[sct][sgm][pmt]
EFF_LSE=[[[0 for k in range(NPMTS)]for j in range(NSGMNTS)]for i in range(NSCTRS)]
EFF_TGT=[[[0 for k in range(NPMTS)]for j in range(NSGMNTS)]for i in range(NSCTRS)]
d=itertools.product(*[range(NSCTRS),range(NSGMNTS),range(NPMTS)])
for r in d:
        isct,isgm,ipmt=r[0],r[1],r[2]

	#! Skip central pmt
	if ipmt==C: continue

	print "isct:isgm:ipmt=%d:%d:%d"%(isct,isgm,ipmt)
	f=hnphe[isct][isgm][ipmt].GetFunction("fnphe")
	print "Obtaining EFF_LSE for isct,isgm,ipmt=",isct,isgm,ipmt
	if f!=None:
		print "Fit function found"
		cut_lse_itgrl=f.Integral(NPHE_CUT_LSE[isct][isgm],1000) #! extend XMAX -> 1000 in the sense of approaching infinity
		cut_tgt_itgrl=f.Integral(NPHE_CUT_TGT[isct][isgm],1000) #! extend XMAX -> 1000 in the sense of approaching infinity
		full_itgrl=f.Integral(0,1000) #! extend XMAX -> 1000 in the sense of approaching infinity
		if full_itgrl!=0:
			EFF_LSE[isct][isgm][ipmt]=cut_lse_itgrl/full_itgrl
			EFF_TGT[isct][isgm][ipmt]=cut_tgt_itgrl/full_itgrl
		else:	
			print "Something wrong with fit function. Setting EFF_LSE=-9999"
			EFF_LSE[isct][isgm][ipmt]=-9999
			EFF_TGT[isct][isgm][ipmt]=-9999
		print "EFF_LSE=",EFF_LSE[isct][isgm][ipmt]
		print "EFF_TGT=",EFF_TGT[isct][isgm][ipmt]
	else:
		print "Fit function NOT found" 

#! Write NPHE_CUT_LSE/TGT data to file in the format:
#! sct sgm
NPHE_CUT_LSE_DATA=open("%s/eid_e16_exp_nphe_cut_lse.txt"%OUTDIR,"w")
NPHE_CUT_TGT_DATA=open("%s/eid_e16_exp_nphe_cut_tgt.txt"%OUTDIR,"w")
d=itertools.product(*[range(NSCTRS),range(NSGMNTS)])
for r in d:
        isct,isgm=r[0],r[1]
        NPHE_CUT_LSE_DATA.write("%d %d %d\n"%(isct+1,isgm+1,NPHE_CUT_LSE[isct][isgm]))
	NPHE_CUT_TGT_DATA.write("%d %d %d\n"%(isct+1,isgm+1,NPHE_CUT_TGT[isct][isgm]))
NPHE_CUT_LSE_DATA.close()
NPHE_CUT_TGT_DATA.close()
 
#! Write EFF_LSE to file in the format:
#! sct sgm pmt EFF_LSE
EFF_LSE_DATA=open("%s/eid_e16_exp_nphe_eff_lse.txt"%OUTDIR,"w")
EFF_TGT_DATA=open("%s/eid_e16_exp_nphe_eff_tgt.txt"%OUTDIR,"w")
d=itertools.product(*[range(NSCTRS),range(NSGMNTS),range(NPMTS)])
for r in d:
	isct,isgm,ipmt=r[0],r[1],r[2]
	
	#! Skip central pmt 
        if ipmt==C: continue

	EFF_LSE_DATA.write("%d %d %d %f\n"%(isct+1,isgm+1,ipmt+1,EFF_LSE[isct][isgm][ipmt]))
	EFF_TGT_DATA.write("%d %d %d %f\n"%(isct+1,isgm+1,ipmt+1,EFF_TGT[isct][isgm][ipmt]))
EFF_LSE_DATA.close()
EFF_TGT_DATA.close()
	
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
