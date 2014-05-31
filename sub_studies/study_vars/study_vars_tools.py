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
import scipy.stats as stats
from scipy.integrate import quad

#atlib
import atlib

from math import *
from collections import OrderedDict

import os
import shutil

W_TH = 1.216

DATADIR=os.environ['STUDY_VARS_DATADIR']
ANADIR=os.environ['STUDY_VARS_ANADIR']
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
	hvar=[]
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
				h.Draw()
				hvar.append(h)
				# h.Fit("gaus")#,"","",fmin,fmax)
				# f=h.GetFunction("gaus")
				# mu,sg=None,None
				# if f:
				# 	mu=f.GetParameter(1)
				# 	sg=f.GetParameter(2)
				c.SaveAs("%s/%s.png"%(outdir,c.GetName())) 
				c.Close()
	return hvar

def gauss_ppi_hack(v, par):
	arg = 0;
	if (par[2] != 0): arg = (v[0] - par[1])/par[2];
	binw=((0.40-0.00)/100)
	fitval = par[0]*(1/(sqrt(2*pi)*par[2]))*exp(-0.5*arg*arg)*binw;
	return fitval;	

def plot_comp_var(hX,XMU,XCUT,OUTDIR,frange):
	ROOT.gStyle.SetOptStat("n")
	ROOT.gStyle.SetOptFit(1111)

	fgauss = ROOT.TF1("fgauss",gauss_ppi_hack,frange[0],frange[1],3)#0.0,0.2,3);
	fgauss.SetParameters(1,0,1);
	fgauss.SetParName(0,"Entries")
	fgauss.SetParName(1,"Mean")
	fgauss.SetParName(2,"Sigma")

	# clrs=["kBlack","kRed","kGreen","kCyan","kBlue","kYellow","kMagenta",
	#   "kBlack+3","kRed+3","kGreen+3","kCyan+3","kBlue+3","kYellow+3","kMagenta+3",
	#   "kBlack+6","kRed+6","kGreen+6","kCyan+6","kBlue+6","kYellow+6","kMagenta+6"]

	clrs=["kRed","kGreen","kBlue","kCyan","kYellow","kMagenta",
		  "kRed+3","kGreen+3","kBlue+3","kCyan+3","kYellow+3","kMagenta+3",
		  "kRed+6","kGreen+6","kBlue+6","kCyan+6","kYellow+6","kMagenta+6",
		  "kRed+9","kGreen+9","kBlue+9","kCyan+9","kYellow+9","kMagenta+9",
		  "kRed+12","kGreen+12","kBlue+12","kCyan+12","kYellow+12","kMagenta+12"]

	#-- Get R
	R=hX.keys()	
	#-- Xmu,Xsg=f(R)
	Xmu=OrderedDict()#{}
	Xsg=OrderedDict()#{}
	for ir,r in enumerate(R):
		Xmu[r]=[]
		Xsg[r]=[]
	
	#-- nq2wbins 
	nq2wbins=len(hX[R[0]])# nq2wbins should be same for all R    

	#-- prepate TLine for mmppip and cut
	l_X=ROOT.TLine(XMU,-50000,0.139,50000)
	l_cut=ROOT.TLine(XCUT,-50000,0.2,50000)


	for iq2w in range(nq2wbins):
		c=ROOT.TCanvas('%02d'%(iq2w+1),'%02d'%(iq2w+1))
		n1=None
		for ir,r in enumerate(R):
			hX[r][iq2w].SetLineColor(ROOT.gROOT.ProcessLine(clrs[ir]))
			hX[r][iq2w].SetMarkerColor(ROOT.gROOT.ProcessLine(clrs[ir]))
			fgauss.SetParameters(1,0,1)
			f=None
			if ir==0: # Directly fit 1st histogram; all others are normalized to the 1st hist
				hX[r][iq2w].Fit("fgauss","","",frange[0],frange[1])#0.0,0.2,3);
				f=hX[r][iq2w].GetFunction("fgauss")
				s=hX[r][iq2w].GetListOfFunctions().FindObject("stats")
				s.SetX1NDC(0.60)
				s.SetX2NDC(1.00)
				s.SetY1NDC(0.9-(ir*0.3))
				s.SetY2NDC(0.6-(ir*0.3))
				s.SetFillStyle(4000)
				s.SetTextColor(ROOT.gROOT.ProcessLine(clrs[ir]));
				c.Update()
				if f:
					n1=f.GetParameter(0)
					f.SetLineColor(ROOT.gROOT.ProcessLine(clrs[ir]))
					Xmu[r].append(f.GetParameter(1)) 
					Xsg[r].append(f.GetParameter(2))
					c.Update()
				else:
					n1=hX[r][iq2w].GetEntries()
					Xmu[r].append(0) 
					Xsg[r].append(0)
			else: # the following histograms are normalized
				hX[r][iq2w].Fit("fgauss","0","",frange[0],frange[1])#0.0,0.2,3);
				f=hX[r][iq2w].GetFunction("fgauss")
				if f:
					n=f.GetParameter(0)
				else:
					n=hX[r][iq2w].GetEntries()           
				norm=n1*(hX[r][iq2w].GetEntries()/n)
				hXn=hX[r][iq2w].DrawNormalized("sames",norm)
				hXn.Fit("fgauss","","sames",frange[0],frange[1])#0.0,0.2,3);
				s=hXn.GetListOfFunctions().FindObject("stats")
				s.SetX1NDC(0.60)
				s.SetX2NDC(1.00)
				s.SetY1NDC(0.9-(ir*0.3))
				s.SetY2NDC(0.6-(ir*0.3))
				s.SetFillStyle(4000)
				s.SetTextColor(ROOT.gROOT.ProcessLine(clrs[ir]));
				c.Update()
				f=hXn.GetFunction("fgauss")
				if f:
					f.SetLineColor(ROOT.gROOT.ProcessLine(clrs[ir]))
					Xmu[r].append(f.GetParameter(1)) 
					Xsg[r].append(f.GetParameter(2))
					c.Update()
				else:
					Xmu[r].append(0) 
					Xsg[r].append(0)
			l_X.Draw("same")
			l_cut.Draw("same")
		c.SaveAs("%s/%s.png"%(OUTDIR,c.GetName()))
		c.Close()
	return (Xmu,Xsg)

def plot_varpar_vs_q2w(Xmu,Xsg,XMU):
	#-- Get R
	R=Xmu.keys()
	clrs=['red','green','blue','yellow']
	fig=plt.figure(figsize=(20,5))
	for ir,r in enumerate(R):
		#print len(Xmu[r])
		clr=np.random.rand(3,1)
		ax=plt.subplot(131)
		ax.scatter(np.arange(len(Xmu[r])),Xmu[r],label=r,color=clrs[ir],s=50)#color=clrs[id])
		ax.set_ylim(0.10,0.2)
		ax.set_xlabel("W-bin")
		ax.set_ylabel("mean")
		ax.hlines(XMU,0,len(Xmu[r])-1)#1,25)
		ax=plt.subplot(132)
		ax.scatter(np.arange(len(Xsg[r])),Xsg[r],label=r,color=clrs[ir],s=50)#color=clrs[id])
		ax.set_ylim(0,0.06)#[0]=0
		ax.set_xlabel("W-bin")
		ax.set_ylabel("sigma")
		ax=plt.subplot(133)
		ax.scatter([1],[1],label=r,color=clrs[ir],s=50)#color=clrs[id])
		ax.legend(loc="center",prop={'size':10})

def plot_varpar_vs_R(Xmu,Xsg,XMU):
	#-- Get R
	R=Xmu.keys()
	clrs=['red','green','blue','yellow']
	fig=plt.figure(figsize=(20,5))

	Xmu_avg=OrderedDict()
	Xsg_avg=OrderedDict()
	for ir,r in enumerate(R):
		Xmu_avg[r]=np.mean(Xmu[r])
		Xsg_avg[r]=np.mean(Xsg[r])

	#clr=np.random.rand(3,1)
	ax=plt.subplot(131)
	ax.scatter(np.arange(len(R)),Xmu_avg.values())#,label=r,color=clrs[ir],s=50)#color=clrs[id])
	#ax.set_ylim(0.10,0.2)
	ax.set_xlabel("R")
	ax.set_xticks(np.arange(len(R)))
	ax.get_xaxis().set_ticklabels(R)
	ax.set_ylabel("mean_mu")
	#ax.hlines(XMU,1,25)
	ax=plt.subplot(132)
	ax.scatter(np.arange(len(R)),Xsg_avg.values())#,label=r,color=clrs[ir],s=50)#color=clrs[id])
	#ax.set_ylim(0,0.06)#[0]=0
	ax.set_xlabel("R")
	ax.set_xticks(np.arange(len(R)))
	ax.get_xaxis().set_ticklabels(R)
	ax.set_ylabel("mean_sg")
	# ax=plt.subplot(133)
	# ax.scatter([1],[1],label=r,color=clrs[ir],s=50)#color=clrs[id])
	# ax.legend(loc="center",prop={'size':10})

# def plot_ana_var_cut(Xmu,Xsg,XMU):
# 	#-- Get R
# 	R=Xmu.keys()
# 	zipped={}
# 	rv={}#rv(r)
# 	eff={}
# 	for ir,r in enumerate(R):
# 		zipped[r] = zip(Xmu[r],Xsg[r])
# 		rv[r]=[stats.norm(i[0],i[1]) for i in zipped[r]]
# 		eff[r]=np.array([quad(i.pdf,0,0.2)[0] for i in rv[r]])
# 	err={}
# 	err['exp_nmcor']=(1-(eff['exp_nmcor']/eff['sim_nmcor']))*100
# 	err['exp_ymcor']=(1-(eff['exp_ymcor']/eff['sim_nmcor']))*100

# 	fig=plt.figure(figsize=(10,5))
# 	fig.suptitle('eff-exp and rel-error in yield(%) vs. W', fontsize=14, fontweight='bold')

# 	plt.subplot(121)
# 	plt.scatter(np.arange(len(eff['exp_nmcor'])),eff['exp_nmcor'],c='red',label='exp_nmcor',s=50)
# 	plt.scatter(np.arange(len(eff['exp_ymcor'])),eff['exp_ymcor'],c='blue',label='exp_ymcor',s=50)
# 	plt.ylim(0.8,1.1)
# 	plt.ylabel("eff-exp")
# 	plt.xlabel("W-bin")
# 	plt.legend()
# 	plt.subplot(122)
# 	plt.scatter(np.arange(len(err['exp_nmcor'])),err['exp_nmcor'],c='red',label='exp_nmcor',s=50)
# 	plt.scatter(np.arange(len(err['exp_ymcor'])),err['exp_ymcor'],c='blue',label='exp_ymcor',s=50)
# 	plt.ylabel("rel-error in yield(%)")
# 	plt.xlabel("W-bin")
# 	plt.ylim(-20,20)
# 	plt.legend()



