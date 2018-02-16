#!/usr/bin/python
from __future__ import division
import ROOT

import atlib as atlib
import q2w_bng

import collections
from array import *
import os,sys
import glob

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, '%s/elast_lite/obs_2pi'%os.environ['ANA2PI'])
from disp_obs import DispObs

from presentation_plots_lib import *

"""
+ Plots ST to ER yields ration in a 2D Q2-W grid

+ >make_2D_ST_2_ER_yields_ratio_plot.py
"""

DTYPS=['ST','ER','rto']
DTYPS_TITLE={'ST':'Thrown Events','ER':'Exp. Events', 'rto':'Thrown/Exp Events'}

#! Set up OBSDIR and SIMNUM
OBSDIR=[0 for i in range(NSIMRNG)]
OBSDIR[L]="%s/lowQ2_SSBands_122717/cutsncors1"%(os.environ['OBSDIR_E16'])
OBSDIR[H]="%s/highQ2_SSBands_122717/cutsncors1"%(os.environ['OBSDIR_E16'])
SIMNUM=[0 for i in range(NSIMRNG)]
SIMNUM[L]='sim4_sim5_sim6_sim7_sim8_sim13'
SIMNUM[H]='sim9_sim10_sim11_sim12'

#! OUTDIR
OUTDIR=os.path.join(os.environ['ANANOTE'],'figures','Acceptance','ST_2_ER_yields_ratio')
#OUTDIR="/tmp/ST_2_ER_yields_ratio"
if not os.path.exists(OUTDIR):
	os.makedirs(OUTDIR)

#! + The following binning should match up with h8_bng.h/q2w_bng.py
#! + It is slightly modified here for "aesthetic" purposes
WBINW=0.025
NWBINS,WMIN,WMAX=36,1.300,2.200# post removing W>2.125 GeV from analysis #68,1.3,3.0 #25 MeV/bin
NQ2BINS_TMP,Q2MIN_TMP,Q2MAX_TMP=8,Q2MIN,Q2MAX # place-holder bng, for it will be modified by variable bng
NBINS_Q2_VAR_BINW_BNG=5
Q2_BINS_VAR_BINW_BNG=array("d",[2.00,2.40,3.00,3.50,4.20,5.00])

#! Create histograms to store simstats(Q2,W) for ST,ER, and ratio
h=collections.OrderedDict()
for dtyp in DTYPS:
	h[dtyp]={}
	for par in ['nbins','N','mu','sg']:
		h[dtyp][par]=ROOT.TH2F("h_%s_%s"%(par,dtyp),"%s_%s(Q2,W)"%(par,dtyp),NWBINS,WMIN,WMAX,NQ2BINS_TMP,Q2MIN_TMP,Q2MAX_TMP)
		#! Set variable binning in Q2
		h[dtyp][par].GetYaxis().Set(NBINS_Q2_VAR_BINW_BNG,Q2_BINS_VAR_BINW_BNG)
		h[dtyp][par].SetXTitle("W [GeV]")
		h[dtyp][par].SetYTitle("Q^{2} [GeV^{2}]")

# #! Create histograms to display simstats(Q2,W)
# h={}
# for par in ['nbins','N','mu','sg']:
# 	h[par]=ROOT.TH2F("h_%s_%s"%(par,'ER'),"%s_%s(Q2,W)"%(par,'ER'),NWBINS,WMIN,WMAX,NQ2BINS_TMP,Q2MIN_TMP,Q2MAX_TMP)
# 	#! Set variable binning in Q2
# 	h[par].GetYaxis().Set(NBINS_Q2_VAR_BINW_BNG,Q2_BINS_VAR_BINW_BNG)
# 	h[par].SetXTitle("W [GeV]")
# 	h[par].SetYTitle("Q^{2} [GeV^{2}]")
		
#! Fill histograms
#! Prepare structure to get stats per sim range
STATS=collections.OrderedDict()
for dtyp in DTYPS:
	if dtyp=='rto': continue #! That will be calculated later once stats for ST and ER are obtained
	STATS[dtyp]=[0 for i in range(NSIMRNG)]
	for isimrng in range(NSIMRNG):
		obsdir=OBSDIR[isimrng]
		simnum=SIMNUM[isimrng]
		dy=DispObs(obsdir,simnum=simnum)#,tops=topcmbn)
		if dtyp=='ST': #! Note that this function returns simstats for ST,SR,SA and SH
			ss=dy.get_sim_stats()
			STATS[dtyp][isimrng]=ss['ST']
		elif dtyp=='ER':
			STATS[dtyp][isimrng]=dy.get_ER_stats()
		else:
			sys.exit("make_2D_ST_2_ER_yields_ratio_plot.py: %s unrecognized"%dtyp)
# print STATS[L]
# sys.exit()
#! Now concatenate NSIMRNG to get STATS_TOT for ST and ER
STATS_TOT=collections.OrderedDict()
for dtyp in DTYPS:
	if dtyp=='rto': continue #! That will be calculated later once stats for ST and ER are obtained
	STATS_TOT[dtyp]=STATS[dtyp][L]+STATS[dtyp][H]

#! Now calculate 'rto' stats from ST and ER
STATS_TOT['rto']=[]
for dST,dER in zip(STATS_TOT['ST'],STATS_TOT['ER']):
	q2_ST,w_ST,nbins_ST,N_ST,mu_ST,sg_ST=dST[0],dST[1],dST[2],dST[3],dST[4],dST[5]
	q2_ER,w_ER,nbins_ER,N_ER,mu_ER,sg_ER=dER[0],dER[1],dER[2],dER[3],dER[4],dER[5]
	#! Calculate ratios for everyting, but q2,w
	nbins_rto=nbins_ST/nbins_ER
	N_rto=N_ST/N_ER
	mu_rto=mu_ST/mu_ER
	sg_rto=sg_ST/sg_ER
	#! Fill STATS_TOT(rto)
	STATS_TOT['rto'].append([q2_ST,w_ST,nbins_rto,N_rto,mu_rto,sg_rto])

#! plot for ST,ER and rto
for dtyp in DTYPS:
	print "Making plot for %s"%dtyp
	for d in STATS_TOT[dtyp]:
		q2,w,nbins,N,mu,sg=d[0],d[1],d[2],d[3],d[4],d[5]
		binx=h[dtyp]['nbins'].GetXaxis().FindBin(w+(0.025/2));
		biny=h[dtyp]['nbins'].GetYaxis().FindBin(q2);
		#print "w,q2,binw,binq2=",w,q2,binx,biny
		bin=h[dtyp]['nbins'].GetBin(binx,biny)
		h[dtyp]['nbins'].SetBinContent(bin,nbins)
		h[dtyp]['N'].SetBinContent(bin,N)
		h[dtyp]['mu'].SetBinContent(bin,mu)
		h[dtyp]['sg'].SetBinContent(bin,sg)			
			
	#! Display "directly interesting" histograms to monitor simstats
	ROOT.gStyle.SetOptStat(0)
	#! N
	c=ROOT.TCanvas("cN_%s"%dtyp,"cN_%s"%dtyp)#,1000,1000)
	h[dtyp]['N'].SetTitle("%s(Q^{2},W)"%(DTYPS_TITLE[dtyp]))
	h[dtyp]['N'].Draw("colz")
	#! Adjust right margin so that palette is not cut off
	c.SetRightMargin(0.13)
	#! Draw Q2-W grid
	lw=[ROOT.TLine() for j in range(len(WBINL))]
	lq=[ROOT.TLine() for j in range(len(Q2BINL))]
	for j,wbin in enumerate(WBINL):
		lw[j].DrawLine(wbin,Q2MIN,wbin,Q2MAX)
	for j,qbin in enumerate(Q2BINL):
		lq[j].DrawLine(WMIN,qbin,WMAX,qbin)
	c.SaveAs("%s/%s_N_2D.png"%(OUTDIR,dtyp))
	c.SaveAs("%s/%s_N_2D.pdf"%(OUTDIR,dtyp))