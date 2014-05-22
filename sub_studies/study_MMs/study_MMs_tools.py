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

#atlib
import atlib

from math import *

import os
import shutil

W_TH = 1.216

DATADIR=os.environ['STUDY_MMS_DATADIR']
ANADIR=os.environ['STUDY_MMS_ANADIR']
OUTDIR=None #OUTDIR = ANADIR/q2wdir

D=None
VARS=None
HRANGE,FRANGE=None,None
#SEL_TOPS=None

PRCSN_Q2,PRCSN_W="%.4f","%.4f"
Q2MIN,Q2MAX,Q2BINW,NQ2BINS,Q2BINS_LE,Q2BINS_UE,Q2BINS=None,None,None,None,None,None,None
WMIN, WMAX, WBINW, NWBINS, WBINS_LE, WBINS_UE, WBINS =None,None,None,None,None,None,None
HOFT,HRES=None,None

def init(q2wdir,q2binw,wbinw,tops,vars,hrange,frange,r):
	"""
	This function ?
    
	Input arguments:
	----------------
	q2wdir: q2wdir relative to DATADIR,list of str
	q2binw: Q2 bin width bins,float
	wbinw: W bin width,float
	tops: Topology to be used, list of int
	vars: Variables for which offset and resolution is to be extracted,list of str
	hrange: for each variable, range for ST-SR histogram,list of [min,max]
	frange: for each variable, range over which ST-SR histogram is fitted,list of [min,max]
	r: vector in multi-dim space which for now is spanned by [dtyp,mcor]. Each variable is a function in this space: var = f(r). 
	   For now r=[dtyp,mcor]; plan to extend it to r=[dtyp,mcor,ecor,effcor]
	   dtyp=(exp,sim)
	   mcor=(ymcor,nmcor)
	"""
	#print init.__doc__
    
    #-- Import rootfile(r) -> DataFrame(d)
	print "-----\nImporting d2pi.root -> DataFrame\n-----"
	global D
	f='%s/%s/%s/recon/d2pi.root'%(DATADIR,q2wdir,r)
	print "Going to use d2pi(%s)"%r
	print "File = %s"%f
	#-- Read in vars that need to be loaded
	global VARS
	VARS=vars
	D=atlib.tree2df(f,'d2piR/tR',vars)
	
	#-- In the DFs, keep data only from the relevant topologies
	print "-----\nTopology selection: tops to be used = ",tops,"\n-----"
	print "dbg:D before"
	print D.head()
	D=atlib.sel_tops(D,tops)
	print "dbg:D after"
	print D.head()
		
	#-- Determine Q2,W binning 
	#-- For this study, reference = ST events
	global Q2MIN,Q2MAX,Q2BINW,NQ2BINS,Q2BINS_LE,Q2BINS_UE,Q2BINS
	global WMIN, WMAX, WBINW, NWBINS, WBINS_LE, WBINS_UE, WBINS
	q2bng,wbng=atlib.init_q2wbng(q2binw,wbinw,range=[(1.9,2.5),(1.3,1.9)])
	Q2MIN,Q2MAX,Q2BINW,NQ2BINS,Q2BINS_LE,Q2BINS_UE,Q2BINS=q2bng
	WMIN, WMAX, WBINW, NWBINS, WBINS_LE, WBINS_UE, WBINS=wbng
	print "-----\nQ2 binning \n-----"
	print "NQ2BINS=%d,Q2BINW=%.4f GeV^2"%(NQ2BINS,Q2BINW)
	print ["%.4f" % i for i in Q2BINS]
	print "-----\nW binning \n-----"
	print "NWBINS=%d,WBINW=%.4f GeV"%(NWBINS,WBINW)
	print ["%.4f" % i for i in WBINS]
	
	#-- get HRANGE,FRANGE
	global HRANGE,FRANGE
	HRANGE=hrange
	FRANGE=frange

	# -- Create OUTDIR
	global OUTDIR
	OUTDIR=os.path.join(ANADIR,q2wdir,r)
	print "-----\nOUTDIR=%s \n-----"%OUTDIR
    
	#-- For the vars, create hres,hofst
	global HOFT,HRES
	HOFT=[]
	HRES=[]
	for var in VARS:
		HOFT.append(Hist2D(WBINS,Q2BINS,name="%s_OFT"%var,title="%s Offset"%var))
		HRES.append(Hist2D(WBINS,Q2BINS,name="%s_RES"%var,title="%s Resolution"%var))
#         HOFT.append(Hist2D(NWBINS,WMIN,WMAX,NQ2BINS,Q2MIN,Q2MAX,name="%s_OFT"%var,title="%s Offset"%var))
#         HRES.append(Hist2D(NWBINS,WMIN,WMAX,NQ2BINS,Q2MIN,Q2MAX,name="%s_RES"%var,title="%s Resolution"%var))
	print 'HOFT=',HOFT
	print 'HRES=',HRES
        
	
def plot_q2w():
	fig=plt.figure(figsize=(30,15))
 	q2=D['Q2']
	w=D['W']
	    
	nq2bins=int((Q2MAX-Q2MIN)/(0.010))#Avg.. Res. Q2~0.010 GeV^2 as obtained for q2w2
	nwbins=int((WMAX-WMIN)/(0.015))#Avg. Res. W~0.015 GeV^2 as obtained for q2w2

	axQ2vW=plt.subplot(221)
	axQ2vW.set_xticks(WBINS)
	axQ2vW.set_yticks(Q2BINS)
	axQ2vW.grid(1,linewidth=2)
	axQ2vW.set_title('Q2 vs. W',fontsize='xx-large')
	axQ2vW.tick_params(axis='both', which='major', labelsize=20)
	atlib.hist2D(w,q2,bins=[nwbins,nq2bins],range=[[WMIN-0.1,WMAX+0.1],[Q2MIN-0.1,Q2MAX+0.1]])

def plot_var(min_entries=-1,max_spreading=1,use_frange=False):
	"""
	Input arguments:
	----------------
	The followng input arguments control the quality the histograms, i.e. by observation, I have decided that 
	the following basic critereon should be satisfied before I consider histogram to be "good"
	enough to "trust" its Fit parameters
	+ min_entries=-1: The minimum number of entries(Integral),int
	+ max_spreading=1: The maximum permissible ratio of (under- + over-flow)/Entries,float
	+ use_frange: If to use user defined Fit range when fitting histograms
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
    
	#-- In each q2wbin, and therein, for each var, plot,fit & save histogam
	#-- 
	for iq2bin in range(NQ2BINS):
		for iwbin in range(NWBINS):
			sel_q2w=(D['Q2']>=Q2BINS_LE[iq2bin])&(D['Q2']<Q2BINS_UE[iq2bin])&(D['W']>=WBINS_LE[iwbin])&(D['W']<WBINS_UE[iwbin])
			for ivar,var in enumerate(VARS):
				data=D[var][sel_q2w]
				#-- Plot, Fit and save hvar
				h=Hist(100,HRANGE[ivar][0],HRANGE[ivar][1],name='%s_%02d_%02d'%(var,iq2bin+1,iwbin+1),
								title='%s for Q2,W=[%.4f,%.4f],[%.4f,%.4f]'%
								(var,Q2BINS_LE[iq2bin],Q2BINS_UE[iq2bin],WBINS_LE[iwbin],WBINS_UE[iwbin]))
				h.fill_array(data)
				
				ROOT.gStyle.SetOptFit(1111)
				c=ROOT.TCanvas(h.GetName())
				h.Fit("gaus")#,"","",fmin,fmax)
				f=h.GetFunction("gaus")
				mu,sg=None,None
				if f:
					mu=f.GetParameter(1)
					sg=f.GetParameter(2)
				c.SaveAs("%s/%s.png"%(outdir,c.GetName())) 
				c.Close()

def plot_MMs_comp(min_entries=-1,max_spreading=1,use_frange=False):
	"""
	Input arguments:
	----------------
	The followng input arguments control the quality of ST-SR histograms i.e. by observation, I have decided that 
	the following basic critereon should be satisfied before I consider a ST-SR histogram to be "good"
	enough to extract OFT and RES
	+ min_entries=-1: The minimum number of entries(Integral),int
	+ max_spreading=1: The maximum permissible ratio of (under- + over-flow)/Entries,float

	+ use_frange: If to use user defined Fit range when fitting ST_SR histograms

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
    
	#-- For Exp and Sim, in each q2wbin, obtain hmm[Exp] and hmm[Sim]:
	#-- 1. Plot and save MM distribution
	#-- 2. Obtain dmu/dsgma =mu/sgma(hmm[Sim])- mu/sgma(hmm[Exp])
	for iq2bin in range(NQ2BINS):
		for iwbin in range(NWBINS):
			#-- For dER & dSR, zoom into q2wbin and get data for hmm
			dfdict={'dER':dER,'dSR':dSR}
			sel_q2w={}
			data={}
			for dfkey in dfdict:
				print "dbg: %s before"%dfkey
				print dfdict[dfkey].head()
				sel_q2w[dfkey]="(%s['Q2']>=%f)&(%s['Q2']<%f)&(%s['W']>=%f)&(%s['W']<%f)"%(
				dfkey,Q2BINS_LE[iq2bin],dfkey,Q2BINS_UE[iq2bin],dfkey,WBINS_LE[iwbin],dfkey,WBINS_UE[iwbin])
				print sel_q2w[dfkey]
				dfdict[dfkey]=dfdict[dfkey][eval(sel_q2w[dfkey])]
				print "dbg: %s after"%dfkey
				print dfdict[dfkey].head()
				for ivar,var in enumerate(VARS):
					data[dfkey]=dfdict[dfkey][var]
			#-- Plot, Fit and save hmm[Exp] and hmm[Sim]
			fmin=0.1
			fmax=0.17
			hER=Hist(100,HRANGE[ivar][0],HRANGE[ivar][1],name='%s_ER_%02d_%02d'%(var,iq2bin+1,iwbin+1),
								title='%s for Q2,W=[%.4f,%.4f],[%.4f,%.4f]'%
								(var,Q2BINS_LE[iq2bin],Q2BINS_UE[iq2bin],WBINS_LE[iwbin],WBINS_UE[iwbin]))
			hER.SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
			hER.SetMarkerColor(ROOT.gROOT.ProcessLine("kBlue"))
			hER.fill_array(data['dER'])
			hSR=Hist(100,HRANGE[ivar][0],HRANGE[ivar][1],name='%s_SR_%02d_%02d'%(var,iq2bin+1,iwbin+1),
								title='%s for Q2,W=[%.4f,%.4f],[%.4f,%.4f]'%
								(var,Q2BINS_LE[iq2bin],Q2BINS_UE[iq2bin],WBINS_LE[iwbin],WBINS_UE[iwbin]))
			hSR.SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
			hSR.SetMarkerColor(ROOT.gROOT.ProcessLine("kRed"))
			hSR.fill_array(data['dSR'])

			ROOT.gStyle.SetOptFit(1111)
			c=ROOT.TCanvas('%s_%02d_%02d'%(var,iq2bin+1,iwbin+1),'%s_%02d_%02d'%(var,iq2bin+1,iwbin+1))
			hSRn=hSR.DrawNormalized("",10000)
			hSRn.Fit("gaus","","sames",fmin,fmax)
			c.Update()
			stSR=hSRn.GetListOfFunctions().FindObject("stats")
			# stSR.SetX1NDC(0.15)
			# stSR.SetX2NDC(0.40)
			stSR.SetY1NDC(0.50)
			stSR.SetY2NDC(0.25)
			stSR.SetTextColor(ROOT.gROOT.ProcessLine("kRed"))
			c.Update()
			hERn=hER.DrawNormalized("sames",10000)
			hERn.Fit("gaus","","sames",fmin,fmax)
			c.Update()
			f=hERn.GetFunction("gaus")
			f.SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
			stER=hERn.GetListOfFunctions().FindObject("stats")
			# stER.SetX1NDC(0.15)
			# stER.SetX2NDC(0.40)
			stER.SetY1NDC(0.75)
			stER.SetY2NDC(0.50)
			stER.SetTextColor(ROOT.gROOT.ProcessLine("kBlue"))
			c.Draw()
			outdir=os.path.join(OUTDIR,var)
			c.SaveAs("%s/%s.png"%(outdir,c.GetName())) 
			c.Close()
			#c.SaveAs("test.eps")			


	

def plot_electron_diagnostics():
	sel=eval(SEL_TOPS)
	theta_min=12
	theta_max=55

	nrows=2
	ncols=5
	fig=plt.figure(figsize=(25,nrows*5))
		
	for i,d in enumerate([dR]):#[dT,dR]):
		mc=''
		if i==0:
			mc='r'
		else:
			mc='b'

		pltnum=0 #keep track of plot number in subplot()
		for ipart,part in enumerate(['e']):
			#pltnum=1
			for ip,p in enumerate([['theta','p'],['px','py']]):
				rng=''
				if ip==0:
					rng=[[5,50],[1,5]]
				else:
					rng=[[-3,3],[-3,3]]
				pltnum=1
				print "debug=",p,pltnum+(ip*ncols),'%s_%s'%(p[0],part),'%s_%s'%(p[1],part)
				#plt.subplot(nrows,ncols,pltnum+(ipart*ncols))#ncols*i+1)
				plt.subplot(nrows,ncols,pltnum+(ip*ncols))#ncols*i+1)
				# plt.xlabel('theta_%s'%part)
				# plt.ylabel('p_%s'%part)
				plt.xlabel('%s_%s'%(p[0],part))
				plt.ylabel('%s_%s'%(p[1],part))
				atlib.hist2D(d['%s_%s'%(p[0],part)][sel],d['%s_%s'%(p[1],part)][sel],bins=100,range=rng)
							
				for ikin,kin in enumerate(['Q2','W']):
					xmin,xmax=0,0
					if kin=='Q2':
						xmin=Q2MIN
						xmax=Q2MAX
					else:
						xmin=WMIN
						xmax=WMAX

					pltnum+=1
					print "debug=",p,pltnum+(ip*ncols),kin,'%s_%s'%(p[0],part)
					#plt.subplot(nrows,ncols,pltnum+(ipart*ncols))#ncols*i+1)
					plt.subplot(nrows,ncols,pltnum+(ip*ncols))#ncols*i+1)
					plt.xlabel(kin)
					#plt.ylabel('p_%s'%part)
					plt.ylabel('%s_%s'%(p[0],part))
					atlib.hist2D(d[kin][sel], d['%s_%s'%(p[0],part)][sel],bins=100,range=[[xmin-0.1,xmax+0.1],rng[0]])
					
					pltnum+=1
					print "debug=",p,pltnum+(ip*ncols),kin,'%s_%s'%(p[1],part)
					#plt.subplot(nrows,ncols,pltnum+(ipart*ncols))#ncols*i+1)
					plt.subplot(nrows,ncols,pltnum+(ip*ncols))#ncols*i+1)
					plt.xlabel(kin)
					#plt.ylabel('theta_%s'%part)
					plt.ylabel('%s_%s'%(p[1],part))
					atlib.hist2D(d[kin][sel], d['%s_%s'%(p[1],part)][sel],bins=100,range=[[xmin-0.1,xmax+0.1],rng[1]])
plt.show()

def plot_hadron_diagnostics(hadrons=['p','pip','pim']):
	sel=eval(SEL_TOPS)
	theta_min=12
	theta_max=55

	nrows=3
	ncols=5
	fig=plt.figure(figsize=(25,nrows*5))
	# nrows=2
	for i,d in enumerate([dR]):#[dT,dR]):
		mc=''
		if i==0:
			mc='r'
		else:
			mc='b'

		pltnum=0 #keep track of plot number in subplot()
		for ipart,part in enumerate(hadrons):
			pltnum=1
			plt.subplot(nrows,ncols,pltnum+(ipart*ncols))#ncols*i+1)
			# plt.scatter(d['theta_%s'%part], d['p_%s'%part],c=mc,alpha=0.5)
			plt.xlabel('theta_%s'%part)
			plt.ylabel('p_%s'%part)
			atlib.hist2D(d['theta_%s'%part][sel],d['p_%s'%part][sel],bins=100,range=[[0,150],[0,3]])
						
			for ikin,kin in enumerate(['Q2','W']):
				xmin,xmax=0,0
				if kin=='Q2':
					xmin=Q2MIN
					xmax=Q2MAX
				else:
					xmin=WMIN
					xmax=WMAX

				pltnum+=1
				plt.subplot(nrows,ncols,pltnum+(ipart*ncols))#ncols*i+1)
				# plt.scatter(d[kin], d['p_%s'%part],c=mc,alpha=0.5)
				plt.xlabel(kin)
				plt.ylabel('p_%s'%part)
				atlib.hist2D(d[kin][sel], d['p_%s'%part][sel],bins=100,range=[[xmin-0.1,xmax+0.1],[0,3]])
				
				pltnum+=1
				plt.subplot(nrows,ncols,pltnum+(ipart*ncols))#ncols*i+1)
				# plt.scatter(d[kin], d['theta_%s'%part],c=mc,alpha=0.5)
				plt.xlabel(kin)
				plt.ylabel('theta_%s'%part)
				atlib.hist2D(d[kin][sel], d['theta_%s'%part][sel],bins=100,range=[[xmin-0.1,xmax+0.1],[0,150]])
plt.show()

def plot_mm(mm):
	#sel=eval(SEL_TOPS)
	#x=dR[mm][sel]
	dfs={'dER':dER,'dSR':dSR}
	# dfs['dER']=(dER)#('dER',dER)
	# dfs['dSR']=(dSR)#'dSR',dSR)
	sel_q2w={}
	wmin=1.7#1.3#1.7
	wmax=2.0#1.5#2.0
	for df in dfs:
		sel_q2w[df]="(%s['W']>=%f)&(%s['W']<%f)"%(df,wmin,df,wmax)
		print sel_q2w[df]
	# plt.hist(dER[mm][eval(sel_q2w['dER'])],bins=100,alpha=0.5,label='ER',histtype='stepfilled')
	# plt.hist(dSR[mm][eval(sel_q2w['dSR'])],bins=100,alpha=0.5,label='SR',histtype='stepfilled')
	# plt.legend()
	fmin=0.1
	fmax=0.17
	hER=Hist(100,0.0,0.2)
	hSR=Hist(100,0.0,0.2)
	ROOT.gStyle.SetOptFit(1111)
	c=ROOT.TCanvas()
	hER.fill_array(dER[mm][eval(sel_q2w['dER'])])
	hSR.fill_array(dSR[mm][eval(sel_q2w['dSR'])])
	hSRn=hSR.DrawNormalized("",10000)
	hSRn.Fit("gaus","","sames",fmin,fmax)
	c.Update()
	stSR=hSRn.GetListOfFunctions().FindObject("stats")
	stSR.SetX1NDC(0.15)
	stSR.SetX2NDC(0.40)
	stSR.SetY1NDC(0.50)
	stSR.SetY2NDC(0.25)
	stSR.SetTextColor(ROOT.gROOT.ProcessLine("kRed"))
	c.Update()
	hERn=hER.DrawNormalized("sames",10000)
	hERn.Fit("gaus","","sames",fmin,fmax)
	c.Update()
	f=hERn.GetFunction("gaus")
	f.SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
	stER=hERn.GetListOfFunctions().FindObject("stats")
	stER.SetX1NDC(0.15)
	stER.SetX2NDC(0.40)
	stER.SetTextColor(ROOT.gROOT.ProcessLine("kBlue"))
	c.Draw()
	c.SaveAs("test.eps")

