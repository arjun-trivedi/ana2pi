#!/usr/bin/python
from __future__ import division

import os,sys,time
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict

import numpy as np

"""
Usage: verify_SFvp_cuts.py  expt=<e1f/e16>

Obtain, from file, parameters for function to make SF cut, for each DTYP and in each SCTR:
fpSFvp[NDTYP][NSCTR][NMTHD][NFNC][NPAR], where
+ NMTHD=number of different ways SF cuts pars in be obtained. Currently 3:
        1. Fit full range of SF in each pbin
        2. Fit around peak of SF in each pin
        3. "current" pars (already being used for analysis)
+ NFNC=Number of functions: SF-high,-low and -mean
+ NPAR=number of parameters of the fit function. Current 4 since fit function = pol3
"""
#! Get args from user
if len(sys.argv)<2:
        sys.exit("Enter expt=e1f/e16")
else:
        expt=sys.argv[1]
        if expt!="e1f"and expt!="e16":
                sys.exit("Valid type for expt=e1f/e16 only")

if expt=="e1f":
        INDIR_EXP=os.path.join(os.environ['D2PIDIR_EXP'],"data_SF_100115")
        INDIR_SIM=os.path.join(os.environ['D2PIDIR_SIM'],"data_SF_100115")
	DATADIR="%s/results_SFvp/cutpars"%os.environ['STUDY_EID_SF_DATADIR']
	OUTDIR="%s/results_SFvp/"%os.environ['STUDY_EID_SF_DATADIR']
elif expt=="e16":
        INDIR_EXP=os.path.join(os.environ['D2PIDIR_EXP_E16'],"data_SF_010516")
        INDIR_SIM=os.path.join(os.environ['D2PIDIR_SIM_E16'],"data_SF_010516")
	DATADIR="%s/results_SFvp/cutpars"%os.environ['STUDY_EID_SF_E16_DATADIR']
	OUTDIR="%s/results_SFvp/"%os.environ['STUDY_EID_SF_E16_DATADIR']
print "INDIR_EXP",INDIR_EXP
print "INDIR_SIM",INDIR_SIM
print "DATADIR=",DATADIR
print "OUTDIR=",OUTDIR
#sys.exit()

NDTYP=2
EXP,SIM=range(2)
DTYP_NAME=['exp','sim']

NSCTR=6

PBINW=0.100
PMIN=0.64
PMAX=5.00
PBIN_LE=np.arange(PMIN,PMAX,PBINW)
NPBIN=len(PBIN_LE)
#print("NPBIN=",NPBIN)

#! The trees are used when testing SF cut
FIN_2pi=[0,0]
T_2pi=[0,0]
FIN_2pi[EXP]=ROOT.TFile("%s/dSF_2pi.root"%INDIR_EXP)
FIN_2pi[SIM]=ROOT.TFile("%s/dSF_2pi.root"%INDIR_SIM)
T_2pi[EXP]=FIN_2pi[EXP].Get("d2piR/tR")
T_2pi[SIM]=FIN_2pi[SIM].Get("d2piR/tR")

if expt=="e1f":	
	FIN_elast=[0,0]
	T_elast=[0,0]
	FIN_elast[EXP]=ROOT.TFile("%s/dSF_elast.root"%INDIR_EXP)
	FIN_elast[SIM]=ROOT.TFile("%s/dSF_elast.root"%INDIR_SIM)
	T_elast[EXP]=FIN_elast[EXP].Get("delast/cut/t")
	T_elast[SIM]=FIN_elast[SIM].Get("delast/cut/t")
elif expt=="e16": #! no elast sim  and therefore no elast data for e16
	FIN_elast=[0]
	T_elast=[0]
	FIN_elast[EXP]=ROOT.TFile("%s/dSF_elast.root"%INDIR_EXP)#os.environ['D2PIDIR_EXP'])
	T_elast[EXP]=FIN_elast[EXP].Get("delast/cut/t")

#! Create data structure fpSFvp: fpSFvp[NDTYP][NSCTR][NMTHD][NFNC][NPAR]
NMTHD=3 #! fit (1) Full, (2) peak of SF in each pbin and (3) "current" pars
F,P,C=range(NMTHD)
MTHD_NAME=['fullSF',"peakSF","current"]
MTHD_CLR=['kBlack','kMagenta','kOrange']
MTHD_LINE_STYLE=[1,2] #1=regular, 2=dashed
MTHD_MRKR_STYLE=['kOpenCircle','kOpenSquare','kOpenTriangleUp'] #1=regular, 2=dashed
MTHD_MRKR_STYLE=['kFullDotSmall','kFullDotSmall','kFullDotSmall']

NFNC=3 #! h,m,l=(high,low,mean)
H,L,M=range(NFNC)
FNC_NAME=['high','low','mean']
	
NPAR=4 #! using pol3 (1 more for constant!)

fpSFvp=[[[[[0 for ipar in range(NPAR)] for icut in range(NFNC)] for imthd in range(NMTHD)] for isctr in range(NSCTR)] for idtyp in range(NDTYP)]

#! Read files and fill in data structure
#DATADIR="%s/results_SFvp/cutpars"%os.environ['STUDY_EID_SF']
#DATADIR="%s/results_SFvp/cutpars"%os.environ['STUDY_EID_SF_DATADIR']
for idtyp in range(NDTYP):
	for imthd in range(NMTHD):
		fname="%s/%s_%s.txt"%(DATADIR,DTYP_NAME[idtyp],MTHD_NAME[imthd])
		with open(fname, 'r') as f:
			data=f.readlines()

		if MTHD_NAME[imthd]!="current":
			for line in data:
				words=line.split()
				for isctr in range(NSCTR):
					if "s%d:cut_h"%(isctr+1) in words[0]:
						for ipar in range(NPAR):
							fpSFvp[idtyp][isctr][imthd][H][ipar]=float(words[ipar+1])
					elif "s%d:cut_l"%(isctr+1) in words[0]:
						for ipar in range(NPAR):
							fpSFvp[idtyp][isctr][imthd][L][ipar]=float(words[ipar+1])
			#! Now set [M] calculated as:
			#! M=(H+L)/2 where,
			#! +H=M+3*SG
			#! +L=M-3*SG
			for isctr in range(NSCTR):
				for ipar in range(NPAR):
					high=fpSFvp[idtyp][isctr][imthd][H][ipar]
					low=fpSFvp[idtyp][isctr][imthd][L][ipar]
					fpSFvp[idtyp][isctr][imthd][M][ipar]=(high+low)/2
		elif MTHD_NAME[imthd]=="current": #! then read pars differently
			#! First get parms. that define function for :
			#! + SF_mean vs. p
			#! + SF_sgma vs. p
			mean=[[0 for ipar in range(NPAR)] for isctr in range(NSCTR)]
			sgma=[[0 for ipar in range(NPAR)] for isctr in range(NSCTR)]
			for line in data:
				words=line.split()
				for isctr in range(NSCTR):
					if "s%d:mean"%(isctr+1) in words[0]:
						for ipar in range(NPAR):
							mean[isctr][ipar]=float(words[ipar+1])
					if "s%d:sgma"%(isctr+1) in words[0]:
						for ipar in range(NPAR):
							sgma[isctr][ipar]=float(words[ipar+1])
							
			#! Now fill in fpSFvp[][][H/L/M][ipar] as:
			#! +[H][isctr][ipar]= SF_mean[isctr][ipar]+3*SF_sigma[isctr][ipar]
			#! +[L][isctr][ipar]= SF_mean[isctr][ipar]-3*SF_sigma[isctr][ipar]
			#! +[M][isctr][ipar]= SF_mean[isctr][ipar]
			for isctr in range(NSCTR):
				for ipar in range(NPAR):
					fpSFvp[idtyp][isctr][imthd][H][ipar]=mean[isctr][ipar]+3*sgma[isctr][ipar]
					fpSFvp[idtyp][isctr][imthd][L][ipar]=mean[isctr][ipar]-3*sgma[isctr][ipar]
					fpSFvp[idtyp][isctr][imthd][M][ipar]=mean[isctr][ipar]
#! Print fpSFvp
for idtyp in range(NDTYP):
	for isctr in range(NSCTR):
		for imthd in range(NMTHD):
			print("%s:%s:s%d:cut_h %f %f %f %f\n"%(DTYP_NAME[idtyp],MTHD_NAME[imthd],isctr+1,
				fpSFvp[idtyp][isctr][imthd][H][0],fpSFvp[idtyp][isctr][imthd][H][1],
				fpSFvp[idtyp][isctr][imthd][H][2],fpSFvp[idtyp][isctr][imthd][H][3]))
			print("%s:%s:s%d:cut_l %f %f %f %f\n"%(DTYP_NAME[idtyp],MTHD_NAME[imthd],isctr+1,
				fpSFvp[idtyp][isctr][imthd][L][0],fpSFvp[idtyp][isctr][imthd][L][1],
				fpSFvp[idtyp][isctr][imthd][L][2],fpSFvp[idtyp][isctr][imthd][L][3]))
			print("%s:%s:s%d:cut_m %f %f %f %f\n"%(DTYP_NAME[idtyp],MTHD_NAME[imthd],isctr+1,
				fpSFvp[idtyp][isctr][imthd][M][0],fpSFvp[idtyp][isctr][imthd][M][1],
				fpSFvp[idtyp][isctr][imthd][M][2],fpSFvp[idtyp][isctr][imthd][M][3]))

#! + Now various plots will be made
#! + Create fout to store plots
#OUTDIR="%s/results_SFvp/"%os.environ['STUDY_EID_SF']
#OUTDIR="%s/results_SFvp/"%os.environ['STUDY_EID_SF_DATADIR']
fout=ROOT.TFile("%s/verify_SF_cuts.root"%OUTDIR,"RECREATE")
#! Now draw cut functions
c=[0,0]
l=[0,0]
f=[[[[0 for icut in range(NFNC)] for imthd in range(NMTHD)] for isctr in range(NSCTR)] for idtyp in range(NDTYP)]
for idtyp in range(NDTYP):
	cname="cSF_fits_%s"%DTYP_NAME[idtyp]
	c[idtyp]=ROOT.TCanvas(cname,cname)
	c[idtyp].Divide(3,2)
	for isctr in range(NSCTR):
		#! First fill f structure
		for imthd in range(NMTHD):
			for ifnc in range(NFNC):
				#if FNC_NAME[ifnc]=='mean': continue #! not drawing mean
				fname="%s_s%d_%s_%s"%(DTYP_NAME[idtyp],isctr+1,MTHD_NAME[imthd],FNC_NAME[ifnc])
				f[idtyp][isctr][imthd][ifnc]=ROOT.TF1(fname,"pol3",0,5)
				f[idtyp][isctr][imthd][ifnc].SetLineColor(ROOT.gROOT.ProcessLine(MTHD_CLR[imthd]))
				for ipar in range(NPAR):
					f[idtyp][isctr][imthd][ifnc].FixParameter(ipar,fpSFvp[idtyp][isctr][imthd][ifnc][ipar])
		#! Now draw f structure			
		c[idtyp].cd(isctr+1)
		ctr=0
		for imthd in range(NMTHD):
			#if MTHD_NAME[imthd]!='fullSF': continue
			for ifnc in range(NFNC):
				#if FNC_NAME[ifnc]!='mean': continue #! not drawing mean
				draw_opt=""
				if ctr>0: draw_opt="same"
				f[idtyp][isctr][imthd][ifnc].SetMinimum(0)
				f[idtyp][isctr][imthd][ifnc].SetMaximum(0.5)
				f[idtyp][isctr][imthd][ifnc].Draw(draw_opt)
				ctr+=1	
		if (isctr+1)==1:
			#! legend
			l[idtyp]=ROOT.TLegend(0.1,0.8,0.4,1.0)#,"","NDC");
			l[idtyp].AddEntry(f[idtyp][isctr][F][M],MTHD_NAME[F])
			l[idtyp].AddEntry(f[idtyp][isctr][P][M],MTHD_NAME[P])
			l[idtyp].AddEntry(f[idtyp][isctr][C][M],MTHD_NAME[C])
			l[idtyp].Draw("same")	
c[EXP].Write()
c[SIM].Write()

#! + To test cut functions: hSFvp_pass[NDTYP][NSCTR][NMTHD]
#! + These histograms will later be used to quantify the cuts (see code that follows)
cSFvp_pass=[[0 for imthd in range(NMTHD)] for idtyp in range(NDTYP)]
hSFvp_pass=[[[0 for imthd in range(NMTHD)] for isctr in range(NSCTR)] for idtyp in range(NDTYP)]
#! Create hSFvp_pass[NDTYP][NSCTR][NMTHD] 
for idtyp in range(NDTYP):
	for isctr in range(NSCTR):
			for imthd in range(NMTHD):
				hname="hSFvp_pass_%s_s%d_%s"%(DTYP_NAME[idtyp],isctr+1,MTHD_NAME[imthd])
				hSFvp_pass[idtyp][isctr][imthd]=ROOT.TH2F(hname,hname,100,0,5,100,0,0.5)
#! Now fill hSFvp_pass[NDTYP][NSCTR][NMTHD]
for idtyp in range(NDTYP):
	if expt=="e1f" or (expt=="e16" and idtyp==0):
		trees=[T_2pi[idtyp],T_elast[idtyp]]
	else:
		trees=[T_2pi[idtyp]]
	for tr in trees:
		print "Going to iterate over tree-name,dtyp=",tr.GetName(),",",DTYP_NAME[idtyp]
		for i,ev in enumerate(tr):
			#if i>100: continue
			mom=ev.p
			sf=ev.etot/mom
			isctr=ev.sector-1
			for imthd in range(NMTHD):
				sf_high=f[idtyp][isctr][imthd][H].Eval(mom)
				sf_low=f[idtyp][isctr][imthd][L].Eval(mom)
				if (sf > sf_low and sf < sf_high):
					hSFvp_pass[idtyp][isctr][imthd].Fill(mom,sf)
#! Now draw hSFvp_pass[NDTYP][NSCTR][NMTHD]
for idtyp in range(NDTYP):
	for imthd in range(NMTHD):
		cname="cSFvp_pass_%s_%s"%(DTYP_NAME[idtyp],MTHD_NAME[imthd])
		cSFvp_pass[idtyp][imthd]=ROOT.TCanvas(cname,cname)
		cSFvp_pass[idtyp][imthd].Divide(3,2)
		for isctr in range(NSCTR):
			cSFvp_pass[idtyp][imthd].cd(isctr+1)
			hSFvp_pass[idtyp][isctr][imthd].Draw("colz")
		cSFvp_pass[idtyp][imthd].Write()

#! To plot differences between fits: hdiff[NDTYP][NSCTR][NMTHD]
hdiff=[[[0 for imthd in range(NMTHD)] for isctr in range(NSCTR)] for idtyp in range(NDTYP)]
cdiff=[0 for idtyp in range(NDTYP)]
#! First create and fill hdiff structure
for idtyp in range(NDTYP):
	for isctr in range(NSCTR):
		for imthd in range(NMTHD):
			hname="hSF_fits_diff_%s_s%d_%s"%(DTYP_NAME[idtyp],isctr+1,MTHD_NAME[imthd])
			hdiff[idtyp][isctr][imthd]=ROOT.TH1F(hname,hname,NPBIN,PMIN,PMAX)
			hdiff[idtyp][isctr][imthd].SetMarkerStyle(ROOT.gROOT.ProcessLine(MTHD_MRKR_STYLE[imthd]))
			hdiff[idtyp][isctr][imthd].SetMarkerColor(ROOT.gROOT.ProcessLine(MTHD_CLR[imthd]))
			for ipbin in range(NPBIN):
				pbin_min=round(PBIN_LE[ipbin],2)
				pbin_max=round(PBIN_LE[ipbin]+PBINW,2)
				binx1=hSFvp_pass[idtyp][isctr][imthd].GetXaxis().FindBin(pbin_min)
				binx2=hSFvp_pass[idtyp][isctr][imthd].GetXaxis().FindBin(pbin_max)
				nbinsy=hSFvp_pass[idtyp][isctr][imthd].GetYaxis().GetNbins()
				n=hSFvp_pass[idtyp][isctr][imthd].Integral(binx1,binx2,1,nbinsy)
				hdiff[idtyp][isctr][imthd].SetBinContent(ipbin+1,n)
				hdiff[idtyp][isctr][imthd].Sumw2()
		#! Now compute diff(%) = ((mthd-C)/C)*100
		for imthd in range(NMTHD):
			hdiff[idtyp][isctr][imthd].Add(hdiff[idtyp][isctr][C],-1)
			hdiff[idtyp][isctr][imthd].Divide(hdiff[idtyp][isctr][C])
			hdiff[idtyp][isctr][imthd].Scale(100)
			#! Fix related aesthetics
			#hdiff[idtyp][isctr][imthd].GetYaxis().SetTitle("%")
			hdiff[idtyp][isctr][imthd].SetTitle("#frac{``mthd''-current}{current}(%)")
#! Now plot hdiff structure
for idtyp in range(NDTYP):
	cname="cSF_fits_diff_%s"%(DTYP_NAME[idtyp])
	cdiff[idtyp]=ROOT.TCanvas(cname,cname)
	cdiff[idtyp].Divide(3,2)
	for isctr in range(NSCTR):
		cdiff[idtyp].cd(isctr+1)
		ctr=0
		for imthd in range(NMTHD):
			draw_opt=""
			if ctr>0: draw_opt="same"
			ROOT.gStyle.SetOptStat(0)
			hdiff[idtyp][isctr][imthd].SetMinimum(-5)
			hdiff[idtyp][isctr][imthd].SetMaximum(5)
			hdiff[idtyp][isctr][imthd].Draw("P %s"%draw_opt)
			ctr+=1
		l=ROOT.TLegend(0.1,0.7,0.4,0.9,"Methods")#,"","NDC");
		l.AddEntry(hdiff[idtyp][isctr][F],MTHD_NAME[F],"p")
		l.AddEntry(hdiff[idtyp][isctr][P],MTHD_NAME[P],"p")
		l.AddEntry(hdiff[idtyp][isctr][M],MTHD_NAME[M],"p")
		l.Draw("same")
	cdiff[idtyp].Write()	

#! Done making plots; close output file
# fout.Close() #! Commented out because it messes with TCanvas display in batch mode(?)
		
if not ROOT.gROOT.IsBatch():
	plt.show()
	# wait for you to close the ROOT canvas before exiting
	wait(True)




