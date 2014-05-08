from __future__ import division

#rootpy
from root_numpy import root2array, root2rec, tree2rec
from root_numpy.testdata import get_filepath
from rootpy.plotting import Hist2D, Hist
from rootpy.plotting import root2matplotlib as rplt
from rootpy.interactive import wait
#PyROOT
import ROOT

#scipy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from math import *

import os
import shutil

W_TH = 1.216

DATADIR=os.environ['STUDY_STQ2WRANGE_DATADIR']
ANADIR=os.environ['STUDY_STQ2WRANGE_ANADIR']
OUTDIR="" #OUTDIR = ANADIR/q2wdir

dT=""
dR=""
VARS=""
DRVD_COLS=['px_e','py_e','pz_e','px_p','py_p','pz_p','px_pip','py_pip','pz_pip','px_pim','py_pim','pz_pim']
HRANGE=""
FRANGE=""
SEL_TOPS=""

Q2MIN=""
Q2MAX=""
NQ2BINS=""
Q2BINW=""
Q2BINS_LE=""
Q2BINS_UE=""
WMIN=""
WMAX=""
NWBINS=""
WBINW=""
WBINS_LE=""
WBINS_UE=""

Q2BINS=""
WBINS=""
HOFT=""
HRES=""

#def init(q2wdirs,nq2bins,nwbins,variables,hrange,frange,tops):
def init(q2wdirs,q2binw,wbinw,variables,hrange,frange,tops):
	"""
	This function intializes the necessary information needed in studying the Simulated-Detector
    
	Input arguments:
	----------------
	q2wdirs = list of q2wdirs relative to $STUDY_STQ2WRANGE_DATADIR,list of str
	# nq2bins = Number of Q2 bins,int
	# nwbins = Number of W bins,int
	q2binw = Q2 bin width bins,float
	wbinw= W bin width,float
	variables = Variables for which offset and resolution is to be extracted,list of str
	hrange=for each variable, range for ST-SR histogram,list of [min,max]
	frange=for each variable, range over which ST-SR histogram is fitted,list of [min,max]
	tops = Topology to be used, list of int
	"""
	#print init.__doc__
    
	#-- Import data 
	global dT,dR,VARS
	cols=['top']
	VARS=variables
	if any(items in VARS for items in DRVD_COLS):
		#print "Don't know what to do"
		#return
		drct_cols=[item for item in VARS if item not in DRVD_COLS]
		drvd_cols=[item for item in VARS if item in DRVD_COLS]
		if not 'theta_e' in VARS: drct_cols.extend(['theta_e'])#needed for drvd_cols
		if not 'phi_e' in VARS: drct_cols.extend(['phi_e'])#needed for drvd_cols    
 		cols.extend(drct_cols)
		print "Branches to load into DataFrame=",cols
		#f = os.path.join(DATADIR,q2wdir,'recon/d2pi.root')
		f=[]
		for q2wdir in q2wdirs:
			f.append(os.path.join(DATADIR,q2wdir,'recon/d2pi.root'))
		arrT = root2array(f,'d2piTR/T/tT',branches=cols)
		arrR = root2array(f,'d2piTR/R/tR',branches=cols)
		dT = pd.DataFrame(arrT)
		dR = pd.DataFrame(arrR)
		#-- Add drvd_cols
		for d in (dT,dR):
			d['px_e']=d['p_e']*np.sin(np.deg2rad(d['theta_e']))*np.cos(np.deg2rad(d['phi_e']))
			d['py_e']=d['p_e']*np.sin(np.deg2rad(d['theta_e']))*np.sin(np.deg2rad(d['phi_e']))
			d['pz_e']=d['p_e']*np.cos(np.deg2rad(d['theta_e']))
	else:
		cols.extend(VARS)
		print "Branches to load into DataFrame=",cols
		#f = os.path.join(DATADIR,q2wdir,'recon/d2pi.root')
		f=[]
		for q2wdir in q2wdirs:
			f.append(os.path.join(DATADIR,q2wdir,'recon/d2pi.root'))
		arrT = root2array(f,'d2piTR/T/tT',branches=cols)
		arrR = root2array(f,'d2piTR/R/tR',branches=cols)
		dT = pd.DataFrame(arrT)
		dR = pd.DataFrame(arrR)
        
	#-- get HRANGE,FRANGE
	global HRANGE,FRANGE
	HRANGE=hrange
	FRANGE=frange

	#-- Determine Topologies to be used
	global SEL_TOPS
	SEL_TOPS=""
	itr=0
	for top in tops:
		if itr==0:
			SEL_TOPS+="(dR['top']==%d)"%top
		else:
			SEL_TOPS+="|(dR['top']==%d)"%top
		itr+=1
	print "SEL_TOPS = %s"%SEL_TOPS
    
	#-- Determine kinematic range of Simulated-Thrown (ST) Events
	global Q2MIN,Q2MAX,WMIN,WMAX
	Q2MIN=round(min(dT['Q2']),2)
	Q2MAX=round(max(dT['Q2']),2)
	WMIN=round(min(dT['W']),2)
	WMAX=round(max(dT['W']),2)
	print "Simulated-Thrown Q2,W = [%.2f,%.2f]GeV^2,[%.2f,%.2f]GeV"%(Q2MIN,Q2MAX,WMIN,WMAX)
    
	#-- Determine Q2,W binning
	global NQ2BINS,Q2BINW,Q2BINS_LE,Q2BINS_UE
	global NWBINS,WBINW,WBINS_LE,WBINS_UE
	# NQ2BINS=nq2bins
	# Q2BINW=round((Q2MAX-Q2MIN)/NQ2BINS,2)
	Q2BINW=q2binw
	NQ2BINS=int(round((Q2MAX-Q2MIN)/Q2BINW,0))
	if NQ2BINS==0:NQ2BINS=1
	Q2BINS_LE=[Q2MIN+(i*Q2BINW) for i in range(NQ2BINS)]
	Q2BINS_UE=[Q2BINS_LE[i]+Q2BINW for i in range(NQ2BINS)]
	print "*** Q2 binning ***"
	print "NQ2BINS=%d,Q2BINW=%0.2f GeV^2"%(NQ2BINS,Q2BINW)
	print ["%0.2f" % i for i in Q2BINS_LE]
	print ["%0.2f" % i for i in Q2BINS_UE]
	# NWBINS=nwbins
	# WBINW=round((WMAX-WMIN)/NWBINS,2)
	WBINW=wbinw
	NWBINS=int(round((WMAX-WMIN)/WBINW,0))
	if NWBINS==0:NWBINS=1
	WBINS_LE=[WMIN+(i*WBINW) for i in range(NWBINS)]
  	WBINS_UE=[WBINS_LE[i]+WBINW for i in range(NWBINS)]
	print "*** W binning ***"
	print "NWBINS=%d,WBINW=%0.2f GeV"%(NWBINS,WBINW)
	print ["%0.2f" % i for i in WBINS_LE]
	print ["%0.2f" % i for i in WBINS_UE]
    
	#-- For the variables, create hres,hofst
	global Q2BINS,WBINS
	global HOFT,HRES
	HOFT=[]
	HRES=[]
	Q2BINS=Q2BINS_LE
	Q2BINS.extend([Q2BINS_UE[-1]])
	WBINS=WBINS_LE
	WBINS.extend([WBINS_UE[-1]])
	for var in VARS:
		HOFT.append(Hist2D(WBINS,Q2BINS,name="%s_OFT"%var,title="%s Offset"%var))
		HRES.append(Hist2D(WBINS,Q2BINS,name="%s_RES"%var,title="%s Resolution"%var))
#         HOFT.append(Hist2D(NWBINS,WMIN,WMAX,NQ2BINS,Q2MIN,Q2MAX,name="%s_OFT"%var,title="%s Offset"%var))
#         HRES.append(Hist2D(NWBINS,WMIN,WMAX,NQ2BINS,Q2MIN,Q2MAX,name="%s_RES"%var,title="%s Resolution"%var))
	print 'HOFT=',HOFT
	print 'HRES=',HRES
        
	# -- Create OUTDIR
	global OUTDIR
	#OUTDIR=os.path.join(ANADIR,q2wdir)
	OUTDIR=os.path.join(ANADIR,"_".join(q2wdirs))

	if not os.path.exists(OUTDIR):
		os.makedirs(OUTDIR)
	print "OUTDIR=%s"%OUTDIR

def plot_q2w_boundary(axes):
	axes.hlines(Q2MIN,WMIN,WMAX,colors='red',linewidth=1)#5DFC0A
	axes.hlines(Q2MAX,WMIN,WMAX,colors='red',linewidth=1)
	axes.vlines(WMIN,Q2MIN,Q2MAX,colors='red',linewidth=1)
	axes.vlines(WMAX,Q2MIN,Q2MAX,colors='red',linewidth=1)
def plot_q2_boundary(axes):
	ylim=axes.get_ylim()
	axes.vlines(Q2MIN,ylim[0],ylim[1],colors='red',linewidth=1)#5DFC0A
	axes.vlines(Q2MAX,ylim[0],ylim[1],colors='red',linewidth=1)
def plot_w_boundary(axes):
	ylim=axes.get_ylim()
	axes.vlines(WMIN,ylim[0],ylim[1],colors='red',linewidth=1)
	axes.vlines(WMAX,ylim[0],ylim[1],colors='red',linewidth=1)
def plot_WTh(axes):
	ylim=axes.get_ylim()
	axes.vlines(W_TH,ylim[0],ylim[1],colors='green',linewidth=2,linestyles='dashed')
    
def plot_q2w():
	fig=plt.figure(figsize=(30,15))
 	q2_T=dT['Q2']
	w_T=dT['W']
	sel=eval(SEL_TOPS)
	q2_R=dR[sel]['Q2']
	w_R=dR[sel]['W']
    
	nq2bins=4400
	nwbins=100

	axQ2vW_T=plt.subplot(221)
	axQ2vW_T.set_xticks(WBINS)
	axQ2vW_T.set_yticks(Q2BINS)
	axQ2vW_T.grid(1,linewidth=2)
	axQ2vW_T.set_title('Q2 vs. W: ST',fontsize='xx-large')
	axQ2vW_T.tick_params(axis='both', which='major', labelsize=20)
	hist,xbins,ybins=np.histogram2d(w_T,q2_T,bins=[nwbins,nq2bins],range=[[WMIN-0.1,WMAX+0.1],[Q2MIN-0.1,Q2MAX+0.1]])
#     hist,xbins,ybins=np.histogram2d(w_T,q2_T,bins=[WBINS,Q2BINS])
	extent = [xbins.min(),xbins.max(),ybins.min(),ybins.max()]
	im = plt.imshow(hist.T, interpolation='none', origin='lower', aspect='auto',extent=extent)
	cb = fig.colorbar(im, ax=axQ2vW_T)
	#-- ST region
	plot_q2w_boundary(axQ2vW_T)
	#-- W Th
	#plot_WTh(axQ2vW_T)

	axQ2vW_R=plt.subplot(222)
	axQ2vW_R.set_xticks(WBINS)
	axQ2vW_R.set_yticks(Q2BINS)
	axQ2vW_R.grid(1,linewidth=2)
	axQ2vW_R.set_title('Q2 vs. W: SR',fontsize='xx-large')
	axQ2vW_R.tick_params(axis='both', which='major', labelsize=20)
	hist,xbins,ybins=np.histogram2d(w_R,q2_R,bins=[nwbins,nq2bins],range=[[WMIN-0.1,WMAX+0.1],[Q2MIN-0.1,Q2MAX+0.1]])
	#hist,xbins,ybins=np.histogram2d(w_R,q2_R,bins=[WBINS,Q2BINS])
	extent = [xbins.min(),xbins.max(),ybins.min(),ybins.max()]
	im = plt.imshow(hist.T, interpolation='none', origin='lower', aspect='auto',extent=extent)
	cb = fig.colorbar(im, ax=axQ2vW_R)
	#-- ST region
	plot_q2w_boundary(axQ2vW_R)
	#-- W Th
	#plot_WTh(axQ2vW_R)
    
	axQ2=plt.subplot(223)
	axQ2.set_xticks(Q2BINS)
	axQ2.grid(1,linewidth=2)
	axQ2.set_title('Q2[GeV^2]',fontsize='xx-large')
	axQ2.set_xlabel('Q2[GeV^2]',fontsize='xx-large')
	axQ2.tick_params(axis='both', which='major', labelsize=20)
	r=plt.hist(q2_T,nq2bins,range=(Q2MIN-0.1,Q2MAX+0.1),alpha=0.5,label='ST',histtype='stepfilled')
	r=plt.hist(q2_R,nq2bins,range=(Q2MIN-0.1,Q2MAX+0.1),alpha=0.5,label='SR',histtype='stepfilled')
#     r=plt.hist(q2_T,bins=Q2BINS,histtype='stepfilled',alpha=0.5,label='ST')
#     r=plt.hist(q2_R,bins=Q2BINS,histtype='stepfilled',alpha=0.5,label='SR')
    #-- ST region
	plot_q2_boundary(axQ2)
	axQ2.legend()

	axW=plt.subplot(224)
	axW.set_xticks(WBINS)
	axW.grid(1,linewidth=2)
	axW.set_title('W',fontsize='xx-large')
	axW.set_xlabel('W[GeV]',fontsize='xx-large')
	axW.tick_params(axis='both', which='major', labelsize=20)
	r=plt.hist(w_T,nwbins,range=(WMIN-0.1,WMAX+0.1),alpha=0.5,label='ST',histtype='stepfilled')
	r=plt.hist(w_R,nwbins,range=(WMIN-0.1,WMAX+0.1),alpha=0.5,label='SR',histtype='stepfilled')
#     r=plt.hist(w_T,bins=WBINS,histtype='stepfilled',alpha=0.5,label='ST')
#     r=plt.hist(w_R,bins=WBINS,histtype='stepfilled',alpha=0.5,label='SR')
	#-- ST region
	plot_w_boundary(axW)
	axW.legend()
	#-- W Th
	#plot_WTh(axW)

def plot_res(min_entries=-1,max_spreading=1):
	"""
	Input arguments:
	----------------
	The input arguments control the quality of ST-SR histograms i.e. by observation, I have decided that 
	the following basic critereon should be satisfied before I consider a ST-SR histogram to be "good"
	enough to extract OFT and RES
	+ min_entries=-1: The minimum number of entries(Integral),int
	+ max_spreading=1: The maximum permissible ratio of (under- + over-flow)/Entries,float

	"""
	#-- Clear existing plots from outdir
	for ivar,var in enumerate(VARS):
		outdir=os.path.join(OUTDIR,var)
		#print "outdir=",outdir
		if os.path.exists(outdir):
			shutil.rmtree(outdir)
		os.makedirs(outdir)
		#-- Following were added when running in interactive mode
		HOFT[ivar].Reset()
		HRES[ivar].Reset()
    
	#-- For each q2wbin and therein, for each VAR:
	#-- 1. Plot and save ST-SR distribiotn
	#-- 2. Obtain offset and resolution from ST-SR distribution and fill HOFST,HRES    
	for iq2bin in range(NQ2BINS):
		for iwbin in range(NWBINS):
			#print "Q2=[%0.2f,%0.2f]"%(Q2BINS_LE[iq2bin],Q2BINS_UE[iq2bin])
			#print "W=[%0.2f,%0.2f]"%(WBINS_LE[iwbin],WBINS_UE[iwbin])
			#sel_q2w=(dR['Q2']>=Q2BINS_LE[iq2bin])&(dR['Q2']<Q2BINS_UE[iq2bin])&(dR['W']>=WBINS_LE[iwbin])&(dR['W']<WBINS_UE[iwbin])
			sel_q2w=(dT['Q2']>=Q2BINS_LE[iq2bin])&(dT['Q2']<Q2BINS_UE[iq2bin])&(dT['W']>=WBINS_LE[iwbin])&(dT['W']<WBINS_UE[iwbin])
			sel=sel_q2w&eval(SEL_TOPS)
            
			for ivar,var in enumerate(VARS):
				if len(dT[sel][var])==0: continue #if there is no data, continue 
				diff=dT[sel][var]-dR[sel][var]
				h=Hist(100,HRANGE[ivar][0],HRANGE[ivar][1],name='%s_hdiff_%02d_%02d'%(var,iq2bin+1,iwbin+1),
								title='%s:ST-SR for Q2,W=[%.2f,%.2f],[%.2f,%.2f]'%
								(var,Q2BINS_LE[iq2bin],Q2BINS_UE[iq2bin],WBINS_LE[iwbin],WBINS_UE[iwbin]))
				h.fill_array(diff)
				#-- Require minimum number of entries
				if h.Integral()<min_entries:continue 
				#-- Require that "the spreading is not too much"    
				entries=h.GetEntries()
				uoflow=h.GetBinContent(0)+h.GetBinContent(h.GetNbinsX()+1)
				if float(uoflow)/float(entries)>max_spreading: continue
				c=ROOT.TCanvas(h.GetName(),h.GetName())
				ROOT.gStyle.SetOptStat("nemMrRuo")
				ROOT.gStyle.SetOptFit(1111)
				h.Fit("gaus")#,"","",FRANGE[ivar][0],FRANGE[ivar][1])
				h.Draw()
				fgaus=h.GetFunction("gaus")
				mu,smga='',''
				if not fgaus:
					mu=0
					sgma=0
				else:
					mu=fgaus.GetParameter(1)
					sgma=fgaus.GetParameter(2)
				HOFT[ivar].Fill(WBINS_LE[iwbin]+WBINW/2,Q2BINS_LE[iq2bin]+Q2BINW/2,mu)
				HRES[ivar].Fill(WBINS_LE[iwbin]+WBINW/2,Q2BINS_LE[iq2bin]+Q2BINW/2,sgma)
				#print "w,q2,RMS=%.2f,%.2f,%.2f"%(WBINS_LE[iwbin],Q2BINS_LE[iq2bin],h.GetRMS())
				#c.SaveAs("%s/%02d_%02d.png"%(outdir,iq2bin+1,iwbin+1)) 
				outdir=os.path.join(OUTDIR,var)
				c.SaveAs("%s/%s.png"%(outdir,c.GetName())) 
				c.Close()
    
	#-- For each VAR, plot and save HOFST,HRES
	fig=plt.figure(figsize=(25,25))
	for ivar,var in enumerate(VARS):
		for ihist,hist in enumerate([HOFT[ivar],HRES[ivar]]):
			#-- Using matplotlib 
			plt.subplot(len(VARS),2,(ivar+1)+ivar+ihist,title=hist.GetName());
			plt.xticks(WBINS)
			plt.yticks(Q2BINS)
			plt.grid(1,linewidth=2)
			z = np.array(hist.z()).T
			im = plt.imshow(z,extent=[hist.xedges(0), hist.xedges(-1),
						hist.yedges(0), hist.yedges(-1)],
						interpolation='nearest',
						aspect='auto',
						origin='lower')
			cb = plt.colorbar(im)
			#-- Using TCanvas
			c=ROOT.TCanvas(hist.GetName(),hist.GetName())
			c.SetGrid();
			ROOT.gStyle.SetOptStat("ne")
			hist.Draw("colz");
			outdir=os.path.join(OUTDIR,var)
			c.SaveAs("%s/%s.png"%(outdir,c.GetName()))
			c.Close()
	fig.savefig("%s/simdet_study.eps"%OUTDIR)       
plt.show()

def plot_q2w_diagnostics():
	THETA_MIN=12
	THETA_MAX=55
	fig=plt.figure(figsize=(15,6))
	for i,d in enumerate([dT,dR]):
		plt.subplot(2,4,4*i+1)
		plt.scatter(d['Q2'], d['p_e'])
		plt.xlabel('Q2')
		plt.ylabel('p_e')

		plt.subplot(2,4,4*i+2)
		plt.scatter(d['Q2'], d['theta_e'])
		plt.xlabel('Q2')
		plt.ylabel('theta_e')
		xlim=plt.xlim()
		plt.hlines(THETA_MIN,xlim[0],xlim[1],colors='r')
		plt.hlines(THETA_MAX,xlim[0],xlim[1],colors='r')
		
		plt.subplot(2,4,4*i+3)
		plt.scatter(d['W'], d['p_e'])
		plt.xlabel('W')
		plt.ylabel('p_e')

		plt.subplot(2,4,4*i+4)
		plt.scatter(d['W'], d['theta_e'])
		plt.xlabel('W')
		plt.ylabel('theta_e')
		xlim=plt.xlim()
		plt.hlines(THETA_MIN,xlim[0],xlim[1],colors='r')
		plt.hlines(THETA_MAX,xlim[0],xlim[1],colors='r')
plt.show()