#!/usr/bin/python
from __future__ import division

import os,sys,time
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict

import numpy as np

"""
Usage: verify_SFvp_cuts.py  expt=<e1f/e16> date
+ date= date tag for used for results_SFvp_<date>

+ This script verifies the cut parameters in results_SFvp_<date>/cutpars by plotting the TF1s

Obtain, from file, parameters for function to make SF cut, for each DTYP and in each SCTR:
fpSFvp[NDTYP][NSCTR][NMTHD][NFNC][NPAR], where
+ NMTHD=number of different ways SF cuts pars in be obtained. Currently 3:
        1. fullSF-lte: "lte style" (pars relate to low/high of SF cut)
        2. peakSF-lte: "lte style" (pars relate to low/high of SF cut)
        3. fullSF-hvy: "hvy style" (should be same as fullSF-lte, but transformed to be written in mu/sg i.e. 'hvy' style so that they can be directly used to 'hvy' code)
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

if len(sys.argv)<3:
	sys.exit("Enter date=MMDDYYYY")
else:
	date=sys.argv[2]
	
#! Exit if expt=e1f
if expt=="e1f": sys.exit("Not implemented for e1f")

if expt=="e1f":
	INDIR_EXP=os.path.join(os.environ['D2PIDIR_EXP'],"data_SF_100115")
	INDIR_SIM=os.path.join(os.environ['D2PIDIR_SIM'],"data_SF_100115")
	DATADIR="%s/results_SFvp/cutpars"%os.environ['STUDY_EID_SF_DATADIR']
	OUTDIR="%s/results_SFvp/"%os.environ['STUDY_EID_SF_DATADIR']
elif expt=="e16":
	INDIR_EXP=os.path.join(os.environ['D2PIDIR_EXP_E16'],"data_SF_an_011517")#!"data_SF_010516"
	INDIR_SIM=os.path.join(os.environ['D2PIDIR_SIM_E16'],"data_SF_an_011517")#!"data_SF_010516"
	DATADIR="%s/results_SFvp_%s/cutpars"%(os.environ['STUDY_EID_SF_E16_DATADIR'],date)
	OUTDIR ="%s/results_SFvp_%s/"%(os.environ['STUDY_EID_SF_E16_DATADIR'],date)

if not os.path.isdir(INDIR_EXP):
	sys.exit("INDIR_EXP=%s does not exist"%INDIR_EXP)
if not os.path.isdir(INDIR_SIM):
	sys.exit("INDIR_SIM=%s does not exist"%INDIR_SIM)
if not os.path.isdir(DATADIR):
	sys.exit("DATADIR=%s does not exist"%DATADIR)
print "INDIR_EXP",INDIR_EXP
print "INDIR_SIM",INDIR_SIM
print "DATADIR=",DATADIR
print "OUTDIR=",OUTDIR
#sys.exit()

NDTYP=2
EXP,SIM=range(2)
DTYP_NAME=['exp','sim']

NSCTR=6

#! Momentum binning information
#! + Should be same as in obtain_SFvp_cuts.py (the following is directly copied from there)
#! Initial set up
#! setup momentum binning for EXP and SIM
PBINW=0.100
PMIN=[0 for i in range(NDTYP)]
PMAX=[0 for i in range(NDTYP)]
PBIN_LE=[0 for i in range(NDTYP)]
NPBIN=[0 for i in range(NDTYP)]
#! + [01-19-16] Setup SIM P limits as per ana2pi kinematics
PMIN[SIM]=1.14
PMAX[SIM]=4.24
PBIN_LE[SIM]=np.arange(PMIN[SIM],PMAX[SIM]-PBINW,PBINW)
NPBIN[SIM]=len(PBIN_LE[SIM])
print("NPBIN[SIM]=",NPBIN[SIM])
time.sleep(3)
#! + [01-19-16] Setup EXP P with no limits, mainly for the sake of statistics and related integrity of the overall fit
#!              since because of low statistics at extremes, the overall fit can be biased
PMIN[EXP]=0.64
PMAX[EXP]=4.24#5.00
PBIN_LE[EXP]=np.arange(PMIN[EXP],PMAX[EXP]-PBINW,PBINW)
NPBIN[EXP]=len(PBIN_LE[EXP])
print("NPBIN[EXP]=",NPBIN[EXP])
time.sleep(3)


NMTHD=3 #! fit (1) fullSF-lte, (2) peakSF-lte (3) fullSF-hvy
F,P,C=range(NMTHD)
MTHD_NAME=['fullSF-lte',"peakSF-lte","fullSF-hvy"]
MTHD_CLR=['kBlack','kMagenta','kOrange']
MTHD_LINE_STYLE=[1,2,2] #1=regular, 2=dashed
MTHD_MRKR_STYLE=['kOpenCircle','kOpenSquare','kOpenTriangleUp'] #1=regular, 2=dashed
MTHD_MRKR_STYLE=['kFullDotSmall','kFullDotSmall','kFullDotSmall']

#! Setup files names for each MTHD: FILE_NAME[MTHD][DTYP]
#! + [01-27-17] Note peakSF-lte is no longer being used and its entires are set to None
FILE_NAME=[['exp_fullSF.txt','sim_fullSF.txt'],[None,None],['eid.exp.e16.out_fullSF','eid.mc.e16.out_fullSF']]

NFNC=3 #! h,m,l=(high,low,mean)
H,L,M=range(NFNC)
FNC_NAME=['high','low','mean']
	
NPAR=4 #! using pol3 (1 more for constant!)

#! Create data structure fpSFvp: fpSFvp[NDTYP][NSCTR][NMTHD][NFNC][NPAR]
fpSFvp=[[[[[0 for ipar in range(NPAR)] for icut in range(NFNC)] for imthd in range(NMTHD)] for isctr in range(NSCTR)] for idtyp in range(NDTYP)]

#! Read files and fill in data structure
for idtyp in range(NDTYP):
	for imthd in range(NMTHD):
		#! + [01-27-17] Note peakSF-lte is no longer being used
		if "peakSF" in MTHD_NAME[imthd]: continue

		fname="%s/%s"%(DATADIR,FILE_NAME[imthd][idtyp])#"%s/%s_%s.txt"%(DATADIR,DTYP_NAME[idtyp],MTHD_NAME[imthd])
		with open(fname, 'r') as f:
			data=f.readlines()

		if 'lte' in (MTHD_NAME[imthd]):#!="current":
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
		elif 'hvy' in MTHD_NAME[imthd]:#=="current": #! then read pars differently
			#! First get parms. that define function for :
			#! + mu vs. p
			#! + sg vs. p
			mu=[[0 for ipar in range(NPAR)] for isctr in range(NSCTR)]
			sg=[[0 for ipar in range(NPAR)] for isctr in range(NSCTR)]
			for line in data:
				words=line.split()
				val=float(words[2])
				if 'SFmean' in words[0] or 'SFsigma' in words[0]:
					sctr,tag,par=words[0].split("_")
					#! determine ipar
					if   "0" in par: ipar=0
					elif "1" in par: ipar=1
					elif "2" in par: ipar=2
					elif "3" in par: ipar=3
					#! determines isctr
					isctr=int(sctr)-1
					if   'SFmean'  in words[0]: mu[isctr][ipar]=val
					elif 'SFsigma' in words[0]: sg[isctr][ipar]=val
						
			#! Now fill in fpSFvp[][][H/L/M][ipar] as:
			#! +[H][isctr][ipar]= SF_mean[isctr][ipar]+3*SF_sigma[isctr][ipar]
			#! +[L][isctr][ipar]= SF_mean[isctr][ipar]-3*SF_sigma[isctr][ipar]
			#! +[M][isctr][ipar]= SF_mean[isctr][ipar]
			for isctr in range(NSCTR):
				for ipar in range(NPAR):
					fpSFvp[idtyp][isctr][imthd][H][ipar]=mu[isctr][ipar]+3*sg[isctr][ipar]
					fpSFvp[idtyp][isctr][imthd][L][ipar]=mu[isctr][ipar]-3*sg[isctr][ipar]
					fpSFvp[idtyp][isctr][imthd][M][ipar]=mu[isctr][ipar]
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
outdir="%s/%s"%(OUTDIR,"verify_SF_cuts")
if not os.path.exists(outdir):
	os.makedirs(outdir)
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
			#! + [01-27-17] Note peakSF-lte is no longer being used
			if "peakSF" in MTHD_NAME[imthd]: continue

			for ifnc in range(NFNC):
				#if FNC_NAME[ifnc]=='mean': continue #! not drawing mean
				fname="%s_s%d_%s_%s"%(DTYP_NAME[idtyp],isctr+1,MTHD_NAME[imthd],FNC_NAME[ifnc])
				f[idtyp][isctr][imthd][ifnc]=ROOT.TF1(fname,"pol3",PMIN[idtyp],PMAX[idtyp])
				f[idtyp][isctr][imthd][ifnc].SetRange(PMIN[idtyp],PMAX[idtyp])
				f[idtyp][isctr][imthd][ifnc].SetLineColor(ROOT.gROOT.ProcessLine(MTHD_CLR[imthd]))
				f[idtyp][isctr][imthd][ifnc].SetLineStyle(MTHD_LINE_STYLE[imthd])
				for ipar in range(NPAR):
					f[idtyp][isctr][imthd][ifnc].FixParameter(ipar,fpSFvp[idtyp][isctr][imthd][ifnc][ipar])
		#! Now draw f structure			
		c[idtyp].cd(isctr+1)
		ctr=0
		for imthd in range(NMTHD):
			#! + [01-27-17] Note peakSF-lte is no longer being used
			if "peakSF" in MTHD_NAME[imthd]: continue
			
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
			#! + [01-27-17] Note peakSF-lte is no longer being used
			#l[idtyp].AddEntry(f[idtyp][isctr][P][M],MTHD_NAME[P])
			l[idtyp].AddEntry(f[idtyp][isctr][C][M],MTHD_NAME[C])
			l[idtyp].Draw("same")	
c[EXP].SaveAs("%s/%s.png"%(outdir,c[EXP].GetName()))
c[SIM].SaveAs("%s/%s.png"%(outdir,c[SIM].GetName()))

# #! Done making plots; close output file
# # fout.Close() #! Commented out because it messes with TCanvas display in batch mode(?)
		
# if not ROOT.gROOT.IsBatch():
# 	plt.show()
# 	# wait for you to close the ROOT canvas before exiting
# 	wait(True)




