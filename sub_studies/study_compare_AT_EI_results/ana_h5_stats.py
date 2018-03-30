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
CLRS=[ROOT.gROOT.ProcessLine("kViolet-7"),
      ROOT.gROOT.ProcessLine("kBlue"),
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

#! Titles for SEQ
SEQ_TITLES={'ER':"ER^{5}",
		    'ST':"ST^{5}",
		    'SR':"SR^{5}",
		    'SA':"SA^{5}"
}

SEQ_TITLES_PYPLOT={'ER':"$ER^{5}$",
		           'ST':"$ST^{5}$",
		           'SR':"$SR^{5}$",
		           'SA':"$SA^{5}$"
}

#! Simulation information  (Taken from ana_obs_as_function_of_simstats.py)
#! Setup Q2 ranges
NQ2RANGES=2
LQ2,HQ2=range(NQ2RANGES)
#! Setup simulations
ISIM=[0 for i in range(NQ2RANGES)]
ISIM[LQ2] =[0,1,2,3,4,5,6]
ISIM[HQ2]=[0,1,2,3,4]
SIMNAME=[0 for i in range(NQ2RANGES)]
SIMNAME[LQ2]=['sim4','sim4_sim5', 'sim4_sim5_sim6',  'sim4_sim5_sim6_sim7','sim4_sim5_sim6_sim7_sim8','sim4_sim5_sim6_sim7_sim8_sim13','sim4_sim5_sim6_sim7_sim8_sim13_sim14']
SIMNAME[HQ2]=['sim9','sim9_sim10','sim9_sim10_sim11','sim9_sim10_sim11_sim12','sim9_sim10_sim11_sim12_sim15']
NSIM=[0 for i in range(NQ2RANGES)]
NSIM[LQ2]=len(ISIM[LQ2])
NSIM[HQ2]=len(ISIM[HQ2])
SIMPCNTG=[0 for i in range(NQ2RANGES)]
SIMPCNTG[LQ2]=[14.3,28.6,42.8,57.1,71.4,85.7,100.0] #![16.7,33.3,50.0,66.7,83.3,100.0]
SIMPCNTG[HQ2]=[20.0,40.0,60.0,80.0,100.0]#1[25,50,75,100]
SIMFRCTN=[0 for i in range(NQ2RANGES)]
SIMFRCTN[LQ2]=[0.14,0.28,0.43,0.57,0.71,0.86,1.0]#[0.17,0.33,0.50,0.67,0.83,1.00]
SIMFRCTN[HQ2]=[0.20,0.40,0.60,0.80,1.0]#[0.25,0.50,0.75,1.00]
SIMTOT=[0 for i in range(NQ2RANGES)]
SIMTOT[LQ2]=['~7B']#['~6B']
SIMTOT[HQ2]=['~5B']#['~4B']

def general_plotting_aesthetics():
	#! General aesthetics
	ROOT.gStyle.Reset()
	ROOT.gStyle.SetOptStat("nmMrReiuo")

def plot_h5_stats_vst_var(h5,seq,q2wbin,outdir,show_rel_err_dist=False):
	'''
	+ Currently for vst-var = 1-THETA

	+ h5=h5[isim][seq]
	+ If h5=h5[seq], i.e. function called for a particular sim, 
	  then convert to h5[isim][seq]

	+ show_rel_err_dist not yet implemented for this function
	  + See plot_h5D_stats on how to do it here too
	'''
	outdir_root="%s/h5_stats_vst_var/%s"%(outdir,q2wbin)

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
				print "plot_h5_stats_vst_var: processing q2wbin,seq,sim#,bin#: %s,%s,%d,%d,[%.3f,%.3f)"%(q2wbin,seq,isim+1,ibin+1,binle,binue)
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
					if nsim>1: h1[isim][seq,pstyp].SetName("%s-vst-var-dist_%s-PS_bin%d_sim%s"%(seq,pstyp,ibin+1,isim+1))#("%s_seq%s_sim%s_bin%d"%(pstyp,seq,isim+1,ibin+1))
					else:      h1[isim][seq,pstyp].SetName("%s-vst-var dist_%s-PS_bin%d"%(seq,pstyp,ibin+1))
					h1[isim][seq,pstyp].SetTitle("%s vst-var Distribution in %s PS (Q^{2}_W=%s)"%(SEQ_TITLES[seq],SEQ_TITLES[pstyp],q2wbin))#("%s:%s"%(seq,pstyp))
					h1[isim][seq,pstyp].SetLineColor(CLRS[isim])
					h1[isim][seq,pstyp].SetXTitle("%s",SEQ_TITLES[seq])

		#! Draw 
		#! ROOT plot aesthetics
		if nsim>1: #! Many stats-boxes overlapped -> no point showing them. stat-box shown when nsim=1
			ROOT.gStyle.Reset()
			ROOT.gStyle.SetOptStat(0)
		for pstyp in ['ST','ER']:
			#if pstyp=='ST': continue
			#! Draw h1 and h1 err
			outdir="/%s/%s/%sPS/dist"%(outdir_root,seq,pstyp)
			if not os.path.exists(outdir):
					os.makedirs(outdir)
			c=ROOT.TCanvas()
			for isim in range(nsim):
				draw_opt="hist"
				if isim>0:draw_opt="hist sames"
				h1[isim][seq,pstyp].Draw(draw_opt)
			c.SaveAs("/%s/c%02d_%04.1f-%04.1f.png"%(outdir,ibin+1,binle,binue))
			c.SaveAs("/%s/c%02d_%04.1f-%04.1f.pdf"%(outdir,ibin+1,binle,binue))

			#! Draw stats
			outdir="/%s/%s/%sPS/ratio_empty_2_filled_bins"%(outdir_root,seq,pstyp)
			if not os.path.exists(outdir):
				os.makedirs(outdir)
			ratio=[stats[isim][seq,pstyp][1]/stats[isim][seq,pstyp][0] for isim in range(nsim)]

			x=range(nsim)
			if   nsim==NSIM[LQ2]: xticks=[SIMFRCTN[LQ2][isim] for isim in range(nsim)]
			elif nsim==NSIM[HQ2]: xticks=[SIMFRCTN[HQ2][isim] for isim in range(nsim)]
			elif nsim==1:         xticks=[isim+1 for isim in range(nsim)]
			else:
				sys.exit("ana_h5_stats.py: plot_h5_stats_vst_var():plot ratio_empty_2_filled_bins:nsim=%d not valid"%nsim)

			fig,ax = plt.subplots(figsize=(10,5),nrows=1, ncols=1)
			#fig.suptitle('empty-bins/filled-bins for %s in %s PS (vst-var-bin=%d,$Q^{2}$_W=%s)'%(SEQ_TITLES_PYPLOT[seq],SEQ_TITLES_PYPLOT[pstyp],ibin+1,q2wbin),fontsize='large')
			#fig.suptitle('Fraction of %s=0 in %s PS versus sim-stats fraction' '\n' 'vst-var-bin=%d, $Q^{2}$_W=%s'%(SEQ_TITLES_PYPLOT[seq],SEQ_TITLES_PYPLOT[pstyp],ibin+1,q2wbin),fontsize='large')
			#ax.set_title('Fraction of %s=0 in %s-PS versus sim-stats fraction' '\n' 'vst-var-bin=%d, $Q^{2}$_W=%s'%(SEQ_TITLES_PYPLOT[seq],SEQ_TITLES_PYPLOT[pstyp],ibin+1,q2wbin),fontsize='large')
			ax.set_title('Fraction of %s=0 in %s-PS versus sim-stats fraction (vst-var-bin=%d, $Q^{2}$_W=%s)'%(SEQ_TITLES_PYPLOT[seq],SEQ_TITLES_PYPLOT[pstyp],ibin+1,q2wbin),fontsize='large')
			ax.scatter(x,ratio,s=50)
			#! xaxis
			ax.set_xlabel('sim-stats fraction',size='large')
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
			#ax.set_ylabel('empty-bins/filled-bins',size='xx-large')
			ax.set_ylabel('Fraction of %s=0 in %s-PS'%(SEQ_TITLES_PYPLOT[seq],SEQ_TITLES_PYPLOT[pstyp]), size='large')
			#! save
			fig.savefig("%s/ratio_%02d_%04.1f-%04.1f.png"%(outdir,ibin+1,binle,binue))
			fig.savefig("%s/ratio_%02d_%04.1f-%04.1f.pdf"%(outdir,ibin+1,binle,binue))        
	#! Close all figures to avoid "/python2.7/matplotlib/pyplot.py:412: RuntimeWarning: More than 20 figures have been opened."
	plt.close('all')

def plot_h5_stats(h5,seq,q2wbin,outdir,show_rel_err_dist=False):
	'''
	+ h5=h5[isim][seq]
	+ If h5=h5[seq], i.e. function called for a particular sim, 
	  then convert to h5[isim][seq]

	+ if nsim=1 (i.e. max sim) and when seq=SA and pstyp=ER, then this function also returns 'frctn_SA_zero_in_ER_PS'
	  from all Q2-W bins, where
	  	+ frctn_SA_zero_in_ER_PS is fraction of holes in SA in ER PS:
	  		+ frctn_SA_zero_in_ER_PS=empty-bins/filled-bins when seq=SA and pstyp=ER
	'''
	outdir_root="%s/h5_stats/%s"%(outdir,q2wbin)

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
	
	#! Declare frctn_SA_zero_in_ER_PS='nan' and give it a meaningful value when nsim=1 for seq=SA and pstyp=ER
	frctn_SA_zero_in_ER_PS='nan'

	h1   =[OrderedDict() for isim in range(nsim)]
	stats=[OrderedDict() for isim in range(nsim)]
	if show_rel_err_dist==True:
		h1err=[OrderedDict() for isim in range(nsim)]
	for isim in range(nsim):
		print "plot_h5_stats: processing q2wbin,seq,sim#: %s,%s,%d"%(q2wbin,seq,isim+1)
		#! isim dbg
		# if isim!=5: 
		# 	for seq in ['SA']:#['ST','SR','SA','ER']:
		# 		for pstyp in ['ST','ER']:
		# 			stats[isim][seq,pstyp]=np.zeros(2,'f')
		# 	continue
		for pstyp in ['ST','ER']:
			nbins,xmin,xmax=BNG[seq][0],BNG[seq][1],BNG[seq][2]
			stats[isim][seq,pstyp]=np.zeros(2,'f')
			h1[isim][seq,pstyp]=thntool.GetBinContentDistCommonBins2(h5[isim][seq],h5[isim][pstyp],stats[isim][seq,pstyp],nbins,xmin,xmax)
			if nsim>1: h1[isim][seq,pstyp].SetName("%s-dist_%s-PS_sim%s"%(seq,pstyp,isim+1))#("%s_seq%s_sim%s"%(pstyp,seq,isim+1))
			else:	   h1[isim][seq,pstyp].SetName("%s-dist_%s-PS"%(seq,pstyp))
			h1[isim][seq,pstyp].SetTitle("%s Distribution in %s PS (Q^{2}_W=%s)"%(SEQ_TITLES[seq],SEQ_TITLES[pstyp],q2wbin))#("%s:%s"%(seq,pstyp))
			h1[isim][seq,pstyp].SetLineColor(CLRS[isim])
			h1[isim][seq,pstyp].SetXTitle("%s"%SEQ_TITLES[seq])
			if show_rel_err_dist==True:
				h1err[isim][seq,pstyp]=thntool.GetBinRelErrorDistCommonBins(h5[isim][seq],h5[isim][pstyp],100,0,1.5)
				# h1err[isim][seq,pstyp].SetName("%s_seq%srelerr_sim%s"%(pstyp,seq,isim+1))#("%s_seq%srelerr_sim%s"%(pstyp,seq,isim+1))
				# h1err[isim][seq,pstyp].SetTitle("%srelerr:%s"%(seq,pstyp))
				if nsim>1: h1err[isim][seq,pstyp].SetName("%s-relerr-dist_%s-PS_sim%s"%(seq,pstyp,isim+1))#("%s_seq%s_sim%s"%(pstyp,seq,isim+1))
				else:	   h1err[isim][seq,pstyp].SetName("%s-relerr-dist_%s-PS"%(seq,pstyp))
				h1err[isim][seq,pstyp].SetTitle("%s Relative Error Distribution in %s PS (Q^{2}_W=%s)"%(SEQ_TITLES[seq],SEQ_TITLES[pstyp],q2wbin))#("%s:%s"%(seq,pstyp))
				h1err[isim][seq,pstyp].SetLineColor(CLRS[isim])
				h1err[isim][seq,pstyp].SetXTitle("%s-relerr"%SEQ_TITLES[seq])

	#! Draw 
	#! ROOT plot aesthetics
	if nsim>1: #! Many stats-boxes overlapped -> no point showing them. stat-box shown when nsim=1
		ROOT.gStyle.Reset()
		ROOT.gStyle.SetOptStat(0)
	for pstyp in ['ST','ER']:
		#! Draw h1 and h1err
		outdir="/%s/%s/%sPS/dist"%(outdir_root,seq,pstyp)
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		c=ROOT.TCanvas()
		if show_rel_err_dist==True:
			c.SetCanvasSize(700,700)
			c.Divide(1,2)
			c.cd(1)
			for isim in range(nsim):
				draw_opt="hist"
				if isim>0:draw_opt="hist sames"
				h1[isim][seq,pstyp].Draw(draw_opt)
				#! isim dbg
				# if isim!=5: continue
				# h1[isim][seq,pstyp].Draw("hist")
			c.cd(2)
			for isim in range(nsim):
				draw_opt="hist"
				if isim>0:draw_opt="hist sames"
				h1err[isim][seq,pstyp].Draw(draw_opt)
				#! isim dbg
				# if isim!=5: continue
				# h1[isim][seq,pstyp].Draw("hist")
		else:
			for isim in range(nsim):
				draw_opt="hist"
				if isim>0:draw_opt="hist sames"
				h1[isim][seq,pstyp].Draw(draw_opt)
				#! isim dbg
				# if isim!=5: continue
				# h1[isim][seq,pstyp].Draw("hist")
		c.SaveAs("/%s/c.png"%(outdir))
		c.SaveAs("/%s/c.pdf"%(outdir))

		#! Draw stats
		outdir="/%s/%s/%sPS/ratio_empty_2_filled_bins"%(outdir_root,seq,pstyp)
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		ratio=[stats[isim][seq,pstyp][1]/stats[isim][seq,pstyp][0] for isim in range(nsim)]
		if nsim==1 and seq=='SA' and pstyp=='ER':
			frctn_SA_zero_in_ER_PS=ratio[0]
		
		x=range(nsim)
		if   nsim==NSIM[LQ2]: xticks=[SIMFRCTN[LQ2][isim] for isim in range(nsim)]
		elif nsim==NSIM[HQ2]: xticks=[SIMFRCTN[HQ2][isim] for isim in range(nsim)]
		elif nsim==1:         xticks=[isim+1 for isim in range(nsim)]
		else:
			sys.exit("ana_h5_stats.py: plot_h5_stats_vst_var():plot ratio_empty_2_filled_bins:nsim=%d not valid"%nsim)

		fig,ax = plt.subplots(figsize=(10,5),nrows=1, ncols=1)
		#fig.suptitle('empty-bins/filled-bins for %s in %s PS ($Q^{2}$_W=%s)'%(SEQ_TITLES_PYPLOT[seq],SEQ_TITLES_PYPLOT[pstyp],q2wbin),fontsize='large')
		#fig.suptitle('Fraction of %s=0 in %s PS versus sim-stats fraction' '\n' '$Q^{2}$_W=%s'%(SEQ_TITLES_PYPLOT[seq],SEQ_TITLES_PYPLOT[pstyp],q2wbin),fontsize='large')
		#ax.set_title('Fraction of %s=0 in %s-PS versus sim-stats fraction' '\n' '$Q^{2}$_W=%s'%(SEQ_TITLES_PYPLOT[seq],SEQ_TITLES_PYPLOT[pstyp],q2wbin),fontsize='large')
		ax.set_title('Fraction of %s=0 in %s-PS versus sim-stats fraction ($Q^{2}$_W=%s)'%(SEQ_TITLES_PYPLOT[seq],SEQ_TITLES_PYPLOT[pstyp],q2wbin),fontsize='large')
		ax.scatter(x,ratio,s=50)
		#! xaxis
		ax.set_xlabel('sim-stats fraction',size='large')
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
		ax.set_ylabel('empty-bins/filled-bins',size='xx-large')
		ax.set_ylabel('Fraction of %s=0 in %s-PS'%(SEQ_TITLES_PYPLOT[seq],SEQ_TITLES_PYPLOT[pstyp]), size='large')
		#! save
		fig.savefig("%s/ratio.png"%(outdir))
		fig.savefig("%s/ratio.pdf"%(outdir))  
	#! Close all figures to avoid "/python2.7/matplotlib/pyplot.py:412: RuntimeWarning: More than 20 figures have been opened."
	plt.close('all') 
	#! return
	return frctn_SA_zero_in_ER_PS