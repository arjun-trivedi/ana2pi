#!/usr/bin/python
from __future__ import division
import os,sys,datetime
from collections import OrderedDict
import array

import ROOT

from rootpy.io import root_open, DoesNotExist

import math

import numpy as np
import matplotlib.pyplot as plt

import itertools

#! imports from proc_h8.py
sys.path.insert(0, '%s/elast_lite/obs_2pi'%os.environ['ANA2PI'])
from proc_h8 import H5_DIM
H4_PROJDIM=array.array('i',[H5_DIM['M1'],H5_DIM['M2'],H5_DIM['PHI'],H5_DIM['ALPHA']])
from disp_obs import VAR_NAMES_PLAIN

#! Set up THnTool
ROOT.gROOT.ProcessLine(".L THnTool.C+")
from ROOT import THnTool
thntool=THnTool()

#! Marker colors
#! + Assuming the FYLD will be listed as per simulation stats: hot -> cold = low -> high
CLRS=[ROOT.gROOT.ProcessLine("kBlue"),
	ROOT.gROOT.ProcessLine("kCyan"),
	ROOT.gROOT.ProcessLine("kGreen"),
	ROOT.gROOT.ProcessLine("kYellow"),
	#ROOT.gROOT.ProcessLine("kPink+1"),
	ROOT.gROOT.ProcessLine("kOrange"),
	ROOT.gROOT.ProcessLine("kRed")]

#! Define binning for hN distributions
BNG=OrderedDict()
BNG['ER']=[6,    -0.5, 5.5]
BNG['SA']=[200, -0.02,   0.2]#[1000, -0.5,   0.5]
BNG['SR']=[6,    -0.5, 5.5]#[100,  0.5,   100.5] #500
BNG['ST']=[6,    -0.5, 5.5]#[0,    0,   0] #! default binning

def general_plotting_aesthetics():
	#! General aesthetics
	ROOT.gStyle.Reset()
	ROOT.gStyle.SetOptStat("nmMrReiuo")

def plot_simstats_5D_vst_var(h5,q2wbin,outdir,show_rel_err_dist=False):
	'''
	+ Currently for vst-var = 1-THETA

	+ h5=h5[isim][seq]
	+ If h5=h5[seq], i.e. function called for a particular sim, 
	  then convert to h5[isim][seq]

	+ show_rel_err_dist not yet implemented. 
	  + To see how I had done it once, see branch ananote's commit at or before de09d220d65258e3e18a2b64f556aaa330aaa4ce (Date:   Thu Nov 9 08:01:25 2017 -0500)
		+ $SUBSTUDIES/study_compare_AT_EI_results/ana_obs_as_function_of_simstats.py:plot_simstats(h5l,q2wbin,q2r)
		OR
		+ $SUBSTUDIES/study_compare_AT_EI_results/study_yields_acceptance.py:plot(q2wbin)
	'''
	outdir_root="%s/simstats_5D_vst_var/%s"%(outdir,q2wbin)

	#! general aesthetics
	general_plotting_aesthetics()

	if isinstance(h5,list): #! h5=h5[isim][seq]
		nsim=len(h5)
	else: #! h5=h5[seq]
		#! convert h5[seq]->h5[isim][seq]
		#! 1. backup original h5
		h5tmp=h5
		#! 2. Create new converted version structure
		nsim=1
		h5=[0 for isim in range(nsim)]
		h5[0]=h5tmp

	#! Plot 5D-PS-vst-var-dist
	for seq in ['SA']: 
		#! Get VST1-THETA binning information from isim=0 (can use any isim)
		nbins=h5[0][seq].GetAxis(H5_DIM['THETA']).GetNbins()
		#print nbins
		binw=h5[0][seq].GetAxis(H5_DIM['THETA']).GetBinWidth(1)
		for ibin in range(nbins):
			h1   =[OrderedDict() for isim in range(nsim)]
			stats=[OrderedDict() for isim in range(nsim)]
			for isim in range(nsim):
					#if isim!=5: continue
					#! Set range and make projections as per bin
					binle=h5[isim][seq].GetAxis(H5_DIM['THETA']).GetBinLowEdge(ibin+1)
					binue=binle+binw
					print "processing seq,sim#,bin#: %s,%d,%d,[%.3f,%.3f)"%(seq,isim+1,ibin+1,binle,binue)
					#! Set range for seq, 'ST' and 'ER'
					h5[isim]['ST'].GetAxis(H5_DIM['THETA']).SetRange(ibin+1,ibin+1)
					h5[isim]['ER'].GetAxis(H5_DIM['THETA']).SetRange(ibin+1,ibin+1)
					h5[isim][seq].GetAxis(H5_DIM['THETA']).SetRange(ibin+1,ibin+1)
					#! h5->h4 for selected range
					h4=OrderedDict()
					h4['ST']=h5[isim]['ST'].Projection(4,H4_PROJDIM,"E")
					h4['ER']=h5[isim]['ER'].Projection(4,H4_PROJDIM,"E")
					h4[seq] =h5[isim][seq].Projection(4,H4_PROJDIM,"E")
					#! get histograms and related stats
					nbins,xmin,xmax=BNG[seq][0],BNG[seq][1],BNG[seq][2]
					for pstyp in ['ST','ER']:
							#if pstyp=='ST': continue
							stats[isim][seq,pstyp]=np.zeros(2,'f')
							h1[isim][seq,pstyp]=thntool.GetBinContentDistCommonBins2(h4[seq],h4[pstyp],stats[isim][seq,pstyp],nbins,xmin,xmax)
							h1[isim][seq,pstyp].SetName("%s_seq%s_sim%s_bin%d"%(pstyp,seq,isim+1,ibin+1))
							h1[isim][seq,pstyp].SetTitle("%s:%s"%(seq,pstyp))
							h1[isim][seq,pstyp].SetLineColor(CLRS[isim])

			#! Draw 
			for pstyp in ['ST','ER']:
				#if pstyp=='ST': continue
				#! Draw h1
				outdir="/%s/%s/%sPS/dist"%(outdir_root,seq,pstyp)
				if not os.path.exists(outdir):
						os.makedirs(outdir)
				c=ROOT.TCanvas()
				for isim in range(nsim):
					draw_opt="hist"
					if isim>0:draw_opt="hist sames"
					h1[isim][seq,pstyp].Draw(draw_opt)
				c.SaveAs("/%s/c%02d_%04.1f-%04.1f.png"%(outdir,ibin+1,binle,binue))
				#! Draw stats
				outdir="/%s/%s/%sPS/stats"%(outdir_root,seq,pstyp)
				if not os.path.exists(outdir):
					os.makedirs(outdir)
				ratio=[stats[isim][seq,pstyp][1]/stats[isim][seq,pstyp][0] for isim in range(nsim)]
				x=range(nsim)
				xticks=[isim+1 for isim in range(nsim)]
				fig,ax = plt.subplots(figsize=(10,5),nrows=1, ncols=1)
				fig.suptitle('%s:#empty-bins/#filled-bins:bin%d'%(pstyp,ibin+1),fontsize='xx-large')
				ax.scatter(x,ratio,s=50)
				#! xaxis
				ax.set_xlabel('sim#',size='xx-large')
				ax.set_xticks(x)
				ax.set_xticklabels(xticks)
				#! Fix x axis
				if nsim>1:
					# shift half a step to the left
					xmin=(3*x[0]-x[1])/2.
					# shift half a step to the right
					xmax=(3*x[-1]-x[-2])/2.
					ax.set_xlim(xmin,xmax)
				#! yaxis
				ax.set_ylabel('#empty-bin/#filled-bin',size='xx-large')
				#! save
				fig.savefig("%s/stats_%02d_%04.1f-%04.1f.png"%(outdir,ibin+1,binle,binue))        

def plot_simstats_5D(h5,q2wbin,outdir,show_rel_err_dist=False):
	'''
	+ h5=h5[isim][seq]
	+ If h5=h5[seq], i.e. function called for a particular sim, 
	  then convert to h5[isim][seq]

	+ show_rel_err_dist not yet implemented. 
	  + To see how I had done it once, see branch ananote's commit at or before de09d220d65258e3e18a2b64f556aaa330aaa4ce (Date:   Thu Nov 9 08:01:25 2017 -0500)
		+ $SUBSTUDIES/study_compare_AT_EI_results/ana_obs_as_function_of_simstats.py:plot_simstats(h5l,q2wbin,q2r)
		OR
		+ $SUBSTUDIES/study_compare_AT_EI_results/study_yields_acceptance.py:plot(q2wbin)
	'''
	outdir_root="%s/simstats_5D/%s"%(outdir,q2wbin)

	#! general aesthetics
	general_plotting_aesthetics()

	if isinstance(h5,list): #! h5=h5[isim][seq]
		nsim=len(h5)
	else: #! h5=h5[seq]
		#! convert h5[seq]->h5[isim][seq]
		#! 1. backup original h5
		h5tmp=h5
		#! 2. Create new converted version structure
		nsim=1
		h5=[0 for isim in range(nsim)]
		h5[0]=h5tmp
	
	h1   =[OrderedDict() for isim in range(nsim)]
	stats=[OrderedDict() for isim in range(nsim)]
	for isim in range(nsim):
		#! isim dbg
		# if isim!=5: 
		# 	for seq in ['SA']:#['ST','SR','SA','ER']:
		# 		for pstyp in ['ST','ER']:
		# 			stats[isim][seq,pstyp]=np.zeros(2,'f')
		# 	continue
		for seq in ['SA']:#['ST','SR','SA','ER']:
			for pstyp in ['ST','ER']:
				nbins,xmin,xmax=BNG[seq][0],BNG[seq][1],BNG[seq][2]
				stats[isim][seq,pstyp]=np.zeros(2,'f')
				h1[isim][seq,pstyp]=thntool.GetBinContentDistCommonBins2(h5[isim][seq],h5[isim][pstyp],stats[isim][seq,pstyp],nbins,xmin,xmax)
				h1[isim][seq,pstyp].SetName("%s_seq%s_sim%s"%(pstyp,seq,isim+1))
				h1[isim][seq,pstyp].SetTitle("%s:%s"%(seq,pstyp))
				h1[isim][seq,pstyp].SetLineColor(CLRS[isim])
	#! Draw 
	for seq in ['SA']:
		for pstyp in ['ST','ER']:
			#! Draw h1
			outdir="/%s/%s/%sPS/dist"%(outdir_root,seq,pstyp)
			if not os.path.exists(outdir):
				os.makedirs(outdir)
			c=ROOT.TCanvas()
			for isim in range(nsim):
				draw_opt="hist"
				if isim>0:draw_opt="hist sames"
				h1[isim][seq,pstyp].Draw(draw_opt)
				#! isim dbg
				# if isim!=5: continue
				# h1[isim][seq,pstyp].Draw("hist")
			c.SaveAs("/%s/c.png"%(outdir))
			#! Draw stats
			outdir="/%s/%s/%sPS/stats"%(outdir_root,seq,pstyp)
			if not os.path.exists(outdir):
				os.makedirs(outdir)
			ratio=[stats[isim][seq,pstyp][1]/stats[isim][seq,pstyp][0] for isim in range(nsim)]
			x=range(nsim)
			xticks=[isim+1 for isim in range(nsim)]
			fig,ax = plt.subplots(figsize=(10,5),nrows=1, ncols=1)
			fig.suptitle('%s:#empty-bins/#filled-bins'%(pstyp),fontsize='xx-large')
			ax.scatter(x,ratio,s=50)
			#! xaxis
			ax.set_xlabel('sim#',size='xx-large')
			ax.set_xticks(x)
			ax.set_xticklabels(xticks)
			#! Fix x axis
			if nsim>1:
				# shift half a step to the left
				xmin=(3*x[0]-x[1])/2.
				# shift half a step to the right
				xmax=(3*x[-1]-x[-2])/2.
				ax.set_xlim(xmin,xmax)
			#! yaxis
			ax.set_ylabel('#empty-bin/#filled-bin',size='xx-large')
			#! save
			fig.savefig("%s/stats.png"%(outdir)) 