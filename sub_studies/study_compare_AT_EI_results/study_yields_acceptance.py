#!/usr/bin/python
from __future__ import division
import os,sys,datetime
from collections import OrderedDict
import array

import ROOT

import matplotlib.pyplot as plt
from rootpy.plotting import root2matplotlib as rplt
from rootpy.plotting import Hist

from rootpy.io import root_open, DoesNotExist

import math

import numpy as np

import ana_h5_stats

'''
Read documentation for function plot() to see what this code does.
'''

USAGE='study_yields_acceptance.py dbg[=False] stats_show_rel_err_dist[=False] plot_h5_stats_vst_var[=False]'

#! user inputs
DBG=False
if len(sys.argv)>1: #! dbg entered by user
	if    sys.argv[1]=="True":  DBG=True
	elif  sys.argv[1]=="False": DBG=False
	else: sys.exit('DBG=%s is not valid. usage: %s'%(sys.argv[1],USAGE))

STATS_SHOW_REL_ERR_DIST=False
if len(sys.argv)>2: #!  show_rel_err_dist entered by user
	if    sys.argv[2]=="True":  STATS_SHOW_REL_ERR_DIST=True
	elif  sys.argv[2]=="False": STATS_SHOW_REL_ERR_DIST=False
	else: sys.exit('STATS_SHOW_REL_ERR_DIST=%s is not valid. usage: %s'%(sys.argv[2],USAGE))

PLOT_H5_STATS_VST_VAR=False
if len(sys.argv)>3: #!  show_rel_err_dist entered by user
	if    sys.argv[3]=="True":  PLOT_H5_STATS_VST_VAR=True
	elif  sys.argv[3]=="False": PLOT_H5_STATS_VST_VAR=False
	else: sys.exit('PLOT_H5_STATS_VST_VAR=%s is not valid. usage: %s'%(sys.argv[3],USAGE))

print "DBG=",DBG
print "STATS_SHOW_REL_ERR_DIST=",STATS_SHOW_REL_ERR_DIST
print "PLOT_H5_STATS_VST_VAR=",PLOT_H5_STATS_VST_VAR
#sys.exit()

#! imports from proc_h8.py
sys.path.insert(0, '%s/elast_lite/obs_2pi'%os.environ['ANA2PI'])
from proc_h8 import H5_DIM
H4_PROJDIM=array.array('i',[H5_DIM['M1'],H5_DIM['M2'],H5_DIM['PHI'],H5_DIM['ALPHA']])
from disp_obs import VAR_NAMES_PLAIN

#! Set up THnTool
ROOT.gROOT.ProcessLine(".L THnTool.C+")
from ROOT import THnTool
thntool=THnTool()

#! Set up general ROOT plotting aesthetics
ROOT.gStyle.SetOptStat("nmMrReiuo")
ROOT.gStyle.SetStatStyle(0) #! transparent stats box

#! OUTDIR
DATE=datetime.datetime.now().strftime('%m%d%y')
OUTDIRNAME='results'
#! append identifiers
if STATS_SHOW_REL_ERR_DIST==True:
	OUTDIRNAME+='_with_rel_err_dist'
else:
	OUTDIRNAME+='_without_rel_err_dist'
OUTDIRNAME+='_%s'%DATE
#! Finally create OUTDIR
if DBG==True:
	OUTDIR=os.path.join(os.environ['STUDY_YIELDS_ACCEPTANCE_DATADIR'],'dbg',OUTDIRNAME)
else:
	OUTDIR=os.path.join(os.environ['STUDY_YIELDS_ACCEPTANCE_DATADIR'],OUTDIRNAME)
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)
print "OUTDIR=",OUTDIR
#sys.exit()

#! Setup input root files
#! + As per SSBANDS: on-on=cutsncors1, off-off=cutsncors2
#! + Different file for Q2 ranges (lowQ2,highQ2)
NQ2RANGES=2
LQ2,HQ2=range(NQ2RANGES)
FIN=[0 for i in range(NQ2RANGES)]
# FIN[LQ2]=root_open('%s/lowQ2_SSBands_080217/cutsncors1/sim4_sim5_sim6_sim7_sim8_sim13/yield.root'%(os.environ['OBSDIR_E16']),'r')
# FIN[HQ2]=root_open('%s/highQ2_SSBands_080217/cutsncors1/sim9_sim10_sim11_sim12/yield.root'%(os.environ['OBSDIR_E16']),'r')
FIN[LQ2]=root_open('%s/lowQ2_SSBands_121517/cutsncors1/sim4_sim5_sim6_sim7_sim8_sim13/yield.root'%(os.environ['OBSDIR_E16']),'r')
FIN[HQ2]=root_open('%s/highQ2_SSBands_121517/cutsncors1/sim9_sim10_sim11_sim12/yield.root'%(os.environ['OBSDIR_E16']),'r')


def get_q2wbinlist(q2min=0.00,q2max=6.00,wmin=0.000,wmax=3.000,dbg=False,dbg_bins=4,
	               dbg_binl=['2.00-2.40_1.425-1.450','2.00-2.40_1.450-1.475','2.00-2.40_1.475-1.500','2.00-2.40_1.500-1.525']):
		"""
		+ Taken from disp_obs.py and modified accordingly
		+ Note in dbg mode, functions works as expected when 'dbg_bins'=number of bins in 'dbg_binl'
		    + 'dbg_binl' contains bins which need immediate analysis
		"""
		q2wbinl=[]
		
		# print "DispObs::get_q2wbinlist() Going to Q2-W bins from file=",f.GetName()
		# print "DispObs::get_q2wbinlist() q2min,q2max,wmin,wmax=",q2min,q2max,wmin,wmax
		# if dbg==True:
		# 	print "DispObs::get_q2wbinlist() dbg=True"

		brk=False #! technical tool to break out of two nested for loops. Set in second (nested) for loop
		for f in [FIN[LQ2],FIN[HQ2]]:
			if brk==True: break
			for path,dirs,files in f.walk():
				if path=="":continue #! Avoid root path
				path_arr=path.split("/")
				if len(path_arr)==1:
					if DBG==True:
						if path in dbg_binl:
							q2wbinl.append(path)
					else:
						q2wbinl.append(path)
				if DBG==True and len(q2wbinl)==dbg_bins:
					brk=True
					break #! Uncomment/comment -> Get limited q2w-bins/Get all q2w-bins



		# #! Remove q2wbins that are not within [q2min,q2max],[wmin,wmax] 
		# q2wbins_remove=[]
		# for q2wbin in q2wbinl:
		# 	q2bin_le=q2wbin.split("_")[0].split("-")[0]
		# 	q2bin_ue=q2wbin.split("_")[0].split("-")[1]
		# 	wbin_le =q2wbin.split("_")[1].split("-")[0]
		# 	wbin_ue =q2wbin.split("_")[1].split("-")[1]
		# 	if float(q2bin_ue)<=q2min or float(q2bin_le)>=q2max or float(wbin_ue)<=wmin or float(wbin_le)>=wmax:
		# 		q2wbins_remove.append(q2wbin)
		# for q2wbin in q2wbins_remove:
		# 	q2wbinl.remove(q2wbin)

		return q2wbinl

def norm_1D_theta(hTheta):
	#! 1. Create normalization factor histogram
	hDCosTheta=hTheta.Clone("hDCosTheta")
	hDCosTheta.SetTitle("hDCosTheta")
	hDCosTheta.Reset()
	nbins=hTheta.GetNbinsX()
	for ibin in range(nbins):
		theta_a=hTheta.GetBinLowEdge(ibin+1)
		theta_b=hTheta.GetBinLowEdge(ibin+2)# + hTheta.GetBinWidth(ibin+1)
		DCosTheta=math.fabs(math.cos(math.radians(theta_b))-math.cos(math.radians(theta_a)))
		hDCosTheta.SetBinContent(ibin+1,DCosTheta)
		hDCosTheta.SetBinError(ibin+1,0.)
	#! Now divide hTheta by hDCosTheta
	#! Do Sumw2() so that errors are correctly propagated
	hTheta.Sumw2();
	hTheta.Divide(hDCosTheta)
	return hDCosTheta

#! main
if DBG==True:
	q2wbinl=get_q2wbinlist(dbg=True,dbg_bins=1,dbg_binl=['2.00-2.40_1.500-1.525'])
else:
	q2wbinl=get_q2wbinlist()
print "q2wbinl=",q2wbinl

#! Create structure to hold frctn_SA_zero_in_ER_PS from all q2wbin
frctn_SA_zero_in_ER_PS=[]
for q2wbin in q2wbinl:
	print "processing q2wbin=",q2wbin

	#! Get all the relevant h5(seq)
	#! + First establish if to use FIN[LQ2] or FIN[HQ2] for this q2wbin
	if   FIN[LQ2].FindObjectAny(q2wbin)!=None: fin=FIN[LQ2]
	elif FIN[HQ2].FindObjectAny(q2wbin)!=None: fin=FIN[HQ2]
	else:
		sys.exit("study_yields_acceptance.py:plot():Could not establish if to use FIN[LQ2] or FIN[HQ2] for q2wbin=%s"(q2wbin))
	print "fin=",fin.GetName()

	h5=OrderedDict()
	for seq in ['ER','ST','SR','SA']:
		h5[seq]=fin.Get('%s/%s/VST1/h5'%(q2wbin,seq))

	#! plot stats
	#! + Note ret=rctn_SA_zero_in_ER_PS and is only meaningul if seq=SA, else it is 'nan'
	#! yields
	ret=ana_h5_stats.plot_h5_stats(h5,'ER',q2wbin,OUTDIR,show_rel_err_dist=STATS_SHOW_REL_ERR_DIST)
	ret=ana_h5_stats.plot_h5_stats(h5,'ST',q2wbin,OUTDIR,show_rel_err_dist=STATS_SHOW_REL_ERR_DIST)
	ret=ana_h5_stats.plot_h5_stats(h5,'SR',q2wbin,OUTDIR,show_rel_err_dist=STATS_SHOW_REL_ERR_DIST)
	#! acceptance
	ret=ana_h5_stats.plot_h5_stats(h5,'SA',q2wbin,OUTDIR,show_rel_err_dist=STATS_SHOW_REL_ERR_DIST)
	frctn_SA_zero_in_ER_PS.append(ret)
	if PLOT_H5_STATS_VST_VAR:
		#! yields
		ret=ana_h5_stats.plot_h5_stats_vst_var(h5,'ER',q2wbin,OUTDIR,show_rel_err_dist=STATS_SHOW_REL_ERR_DIST)
		ret=ana_h5_stats.plot_h5_stats_vst_var(h5,'ST',q2wbin,OUTDIR,show_rel_err_dist=STATS_SHOW_REL_ERR_DIST)
		ret=ana_h5_stats.plot_h5_stats_vst_var(h5,'SR',q2wbin,OUTDIR,show_rel_err_dist=STATS_SHOW_REL_ERR_DIST)
		#! acceptance
		ret=ana_h5_stats.plot_h5_stats_vst_var(h5,'SA',q2wbin,OUTDIR,show_rel_err_dist=STATS_SHOW_REL_ERR_DIST)

#! Plot frctn_SA_zero_in_ER_PS
#! + I use rplt's Hist to create histogram, because it is simpler to do so
#! + But for plotting, I convert it to ROOT's TH1F
print "frctn_SA_zero_in_ER_PS=", frctn_SA_zero_in_ER_PS
h=Hist(100,0,1)
h.fill_array(frctn_SA_zero_in_ER_PS)
h.SetName('frctn_SA_zero_in_ER_PS')
h.SetTitle("Distribution of Fraction of SA^{5}=0 in ER^{5}-PS(Q^{2},W)")
h1=ROOT.TH1F("n","n",100,0,1)
h1=h
c=ROOT.TCanvas()
h1.Draw("hist")
c.SaveAs("%s/dist_frctn_SA_zero_in_ER_PS.png"%(OUTDIR))
c.SaveAs("%s/dist_frctn_SA_zero_in_ER_PS.pdf"%(OUTDIR))


