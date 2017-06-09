#!/usr/bin/python

#! Imports
from __future__ import division
from math import *

from array import *
from collections import OrderedDict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools
import datetime

from multiprocessing import Process, Queue
import subprocess

from rootpy.io import root_open, DoesNotExist

#PyROOT
import ROOT

import os,sys

import atlib as atlib
import q2w_bng

ROOT.gROOT.ProcessLine(".L THnTool.C+")
from ROOT import THnTool

'''
+ [05-29-17] Started as a copy of $ELAST_LITE/obs_2pi/proc_h8.py

+ This script obtains the correction factor (CF), which is 1 - CNTMNF
 (the contamination fraction= of events from the target (tgt) in the total (H+tgt=ptgt) events):
	+ CF = 1 - CNTMNF, where
	+ CNTMNF = N-tgt/N-ptgt, where 
		+ N-tgt = N-etgt*R
		+ R = Q-ptgt/Q-etgt (obtained in $SUBSTUDIES/obtain_and_validate_R.py)
	=> CF=(N-etgt*R)/N-ptgt

+ This CF is then used to correct the cross-sections as (inpired by PK, p97):
	(dsigma/dX_ij)_corr[q2wb] = (dsigma/dX_ij)[q2wb] * (CF_ij)[q2wb], where
	+ CF_ij: 1 - CNTMNF_ij
	+ ij: W, 9 1D xsecs
	+ j = bin index
	+ CF_ij:  (N-etgt_ij * R)/(N-ptgt_ij)
	+ q2wb: q2w bin. Note:
		+ Binning in q2w can be specified by the user to return results in each (option 'all' )of the ana2pi analysis'
		  q2w bins or integrated over all the q2w bins (option 'none').
		+ When option 'none' is selected i.e. C_ij is obtained using data integrated over all analysis q2w bins,
		  there are some technical caveats that need to kept in mind:
			+ Because h8[Crs-W] in d2piR are directly used in this code and because of this Crs-W binning,
			  it is not possible to directly project out integrated contents over the entire q2w range. 
			+ Therefore, CF_ij[Crw-W] are obtained per Crs-W with integration over all Q2.
			+ These CF_ij[Crw-W] are then later averaged to obtain CF_ij using data integrated over all Q2-W
			+ Note when ij=W, then CF_ij[Crs-W] is not obtained (because of the limited range in W), and 
			  only averages are obtained.

+ Note that CF_ij is returned as a histogram object

+ Output is stored in:
	+ .root format: So that during making observables, CF histograms can be directly obtained and applied in a histogram operation.
	+ .png: For easy viewing and also for saving analysis plots that can be used in the thesis
'''
# Tools 
thntool=THnTool()

#! Get data from user
USAGE='$SUBSTUDIES/study_tgt_BG/e16/obtain_tgt_BG_CF.py q2wbng_opt[=none]<none/all> dbg[=False]<True/False>'

Q2WBNG_OPT='none'
if len(sys.argv)>1: #i.e. q2wbng_opt entered by user
	Q2WBNG_OPT=sys.argv[1]
	if Q2WBNG_OPT!='none' and Q2WBNG_OPT!='all':
		sys.exit("q2wbng_opt can only be none or all. Usage=%s"%USAGE)
print "Q2WBNG_OPT=",Q2WBNG_OPT

DBG=False
if len(sys.argv)>2: #i.e. debug entered by user
	if    sys.argv[2]=="True":  DBG=True
	elif  sys.argv[2]=="False": DBG=False
	else: sys.exit("dbg can only be True or False. Usage=%s"%USAGE)
print "DBG=",DBG

#! Setup other data
Q2MIN,Q2MAX=2.0,5.0
WMIN,WMAX=1.400,2.125
WBINW=0.025

DBGQ2WBIN1='1.50-2.00_1.600-1.625'
DBGQ2WBIN2='2.00-2.40_1.600-1.625'
#! Input data
TGTS=['ptgt_lse','ptgt_tgt','etgt']
NTGTS=len(TGTS)
#TGTS_NAME=['ptgt_lse','ptgt_tgt','etgt']
INDIR=os.path.join(os.environ['D2PIDIR_EXP_E16'],'data_tgt_051817')
FIN=OrderedDict()
FIN['ptgt_lse']=ROOT.TFile(os.path.join(INDIR,'d2piR_ptgt_cutruns_lse.root'))
FIN['ptgt_tgt']=ROOT.TFile(os.path.join(INDIR,'d2piR_ptgt_cutruns_tgt.root'))
FIN['etgt']    =ROOT.TFile(os.path.join(INDIR,'d2piR_etgt_cutruns.root'))

#! Output data
DATE=datetime.datetime.now().strftime('%m%d%y')
OUTDIR=os.path.join(os.environ['STUDY_TGT_BG_E16_DATADIR'],"results_CF_%s"%DATE)
if DBG==True:
	OUTDIR=os.path.join(os.environ['STUDY_TGT_BG_E16_DATADIR'],"results_CF_dbg_%s"%DATE)
if not os.path.exists(OUTDIR):
	os.makedirs(OUTDIR)

#! Setup constants
#! R 
R=OrderedDict() #! Note last index, ETGT, is not used, but used for code-readability
#! + from $STUDY_TGT_BG_E16_DATADIR/results_R_060917/R.txt
R['ptgt_lse']=13.67
R['ptgt_tgt']=13.40

#! H8 and H5 binning
#! See http://stackoverflow.com/questions/16553506/python-ordereddict-iteration for following code
H8_DIM=OrderedDict([('HEL',0),('Q2',1),('W',2),('M1',3),('M2',4),('THETA',5),('PHI',6),('ALPHA',7)])
H5_PROJDIM=array('i',[H8_DIM['M1'],H8_DIM['M2'],H8_DIM['THETA'],H8_DIM['PHI'],
					  H8_DIM['ALPHA']])
	
H5_DIM=OrderedDict([('M1',0),('M2',1),('THETA',2),('PHI',3),('ALPHA',4)])
VSTS=[1,2,3]
VARS=['M1','M2','THETA','PHI','ALPHA']

#! Create dictionary for VAR_NAMES(VST,VAR): 
#! + VAR_NAMES(VST,VAR) should correspond with 'DataAna::MakeYields()':
#! 		+VST1: Delta++ => (p,pip)pim
#! 		+VST2: Rho     => (pip,pim)p
#! 		+VST3: Delta0  => (p,pim)pip
VAR_NAMES={(1,VARS[0]):"M_{p#pi^{+}}",    (1,VARS[1]):"M_{#pi^{+}#pi^{-}}",
		   (1,VARS[2]):"#theta_{#pi^{-}}",(1,VARS[3]):"#phi_{#pi^{-}}",
		   (1,VARS[4]):"#alpha_{[p_{f}#pi^{+}][p#pi^{-}]}",
		   (2,VARS[0]):"M_{p#pi^{+}}",(2,VARS[1]):"M_{#pi^{+}#pi^{-}}",
		   (2,VARS[2]):"#theta_{p}",  (2,VARS[3]): "#phi_{p}",
		   (2,VARS[4]):"#alpha_{[#pi^{+}#pi^{-}][pp_{f}]}",
		   (3,VARS[0]):"M_{p#pi^{+}}",    (3,VARS[1]):"M_{p#pi^{-}}",
		   (3,VARS[2]):"#theta_{#pi^{+}}",(3,VARS[3]): "#phi_{#pi^{+}}",
		   (3,VARS[4]):"#alpha_{[p_{f}#pi^{-}][p#pi^{+}]}"}
#! Create a non-coded, simpler version of VAR_NAMES that can be used with plain text
VAR_NAMES_PLAIN={(1,VARS[0]):"M_ppip",   (1,VARS[1]):"M_pippim",
		         (1,VARS[2]):"theta_pim",(1,VARS[3]):"phi_pim",
		         (1,VARS[4]):"alpha_pfpip_ppim",
		         (2,VARS[0]):"M_ppip", (2,VARS[1]):"M_pippim",
		         (2,VARS[2]):"theta_p",(2,VARS[3]): "phi_p",
		         (2,VARS[4]):"alpha_pippim_ppf",
		         (3,VARS[0]):"M_ppip",   (3,VARS[1]):"M_ppim",
		         (3,VARS[2]):"theta_pip",(3,VARS[3]): "phi_pip",
		         (3,VARS[4]):"alpha_pfpim_ppip"}

VAR_UNIT_NAMES={VARS[0]:"GeV",VARS[1]:"GeV",VARS[2]:"deg",VARS[3]:"deg",VARS[4]:"deg"}
VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC={VARS[0]:"GeV",VARS[1]:"GeV",VARS[2]:"",VARS[3]:"rad",VARS[4]:"rad"}

#! Enumerate measured Helicity values as per their corresponding bin number in h8's Helicity dimension (see DataAna::MakrYields())
#! + Only NEG,ZERO & POS are valid measured Helicity values
#! + UNPOL does not correspond to a measured Helicity value, but helps in writing code
UNPOL,NEG,ZERO,POS=-9999,1,2,3
HEL_NAME={UNPOL:'UNPOL',NEG:'NEG',ZERO:'ZERO',POS:'POS'}

#! define some helper functions needed by main code below
def get_h8_seq_drct(iw):
	"""
	get h8(TGT,VST) for each Crs-W
	"""
	h8=OrderedDict()
	#! Fetch h8(SEQ_DRCT,vst) directly from file 
	for vst in VSTS:
		h8['ptgt_lse',vst]=FIN['ptgt_lse'].Get('d2pi/h8_%d_%d'%(iw+1,vst))
		h8['ptgt_tgt',vst]=FIN['ptgt_tgt'].Get('d2pi/h8_%d_%d'%(iw+1,vst))
		h8['etgt',vst]    =FIN['etgt'].Get('d2pi/h8_%d_%d'%(iw+1,vst))
	return h8

def setM1M2axisrange(h,vst,var,wmax):
	if (vst==1 and var=='M1') or (vst==3 and var=='M2'):
		h.GetXaxis().SetRangeUser(1.1000,wmax-0.14)
	elif (vst==2 and var=='M2'):
		h.GetXaxis().SetRangeUser(0.2780,wmax-0.938)
	
def get_q2w_bng(h8,q2wbng_opt):
	'''
	Using an h8 (bins=H,Q2,W,M1,M2,THETA,PHI,ALPHA) and depending on q2wbng_opt (see below), return Q2-W binning information:
	+ q2wb_crds: List of 2D bin coordinates of all Q2-W bins formed by combining the h8's
				 Q2 and W bins. 
	+ q2wb_lmts: List of bin edges for each of the 2D Q-2 bin coordinates.

	+ q2wbng_opt='all': 
		+ Return the 2D binning information for all Q2-W bins in the h8
		+ Enables the user to make projections from each Q2-W bin
	+ q2wbng_opt='none': 
		+ Return the entire expanse of 2D ana2pi-Q2-W (Q2MIN,Q2MAX,WMIN,WMAX defined in code) space as one bin
		+ Enables the user to project out, as a single bin, the entire Q2-W space
	'''
	#! Obtain q2wb_crds
	nq2bins=h8.GetAxis(H8_DIM['Q2']).GetNbins()
	nwbins=h8.GetAxis(H8_DIM['W']).GetNbins()

	if q2wbng_opt=='all':
		q2b_crds=range(1,nq2bins+1)
		wb_crds=range(1,nwbins+1)
		q2wb_crds=list(itertools.product(q2b_crds,wb_crds))

		#! Now obtain q2wb_lmts
		q2wb_lmts=[]
		for q2wb_crd in q2wb_crds:
			q2b_crd=q2wb_crd[0]
			q2b_min=(float)("%.2f"%h8.GetAxis(H8_DIM['Q2']).GetBinLowEdge(q2b_crd))
			q2b_max=(float)("%.2f"%h8.GetAxis(H8_DIM['Q2']).GetBinLowEdge(q2b_crd+1))
   
			wb_crd=q2wb_crd[1]
			wb_min=(float)("%.3f"%h8.GetAxis(H8_DIM['W']).GetBinLowEdge(wb_crd))
			wb_max=(float)("%.3f"%h8.GetAxis(H8_DIM['W']).GetBinLowEdge(wb_crd+1))
	
			q2wb_lmts.append(((q2b_min,q2b_max),(wb_min,wb_max)))
	if q2wbng_opt=='none': #! note that the q2wb_crds has a different form here
		q2b_min=(float)("%.2f"%h8.GetAxis(H8_DIM['Q2']).GetBinLowEdge(1))
		q2b_max=(float)("%.2f"%h8.GetAxis(H8_DIM['Q2']).GetBinLowEdge(nq2bins+1))
		#! Make sure q2b_min/max are within Q2MIN/MAX
		if q2b_min<Q2MIN: q2b_min=Q2MIN
		if q2b_max>Q2MAX: q2b_max=Q2MAX
		#! Now get q2 crds
		q2b_crd_min=h8.GetAxis(H8_DIM['Q2']).FindBin(q2b_min)
		q2b_crd_max=h8.GetAxis(H8_DIM['Q2']).FindBin(q2b_max)

		wb_min=(float)("%.3f"%h8.GetAxis(H8_DIM['W']).GetBinLowEdge(1))
		wb_max=(float)("%.3f"%h8.GetAxis(H8_DIM['W']).GetBinLowEdge(nwbins+1))		
		#! Make sure wb_min/max are within WMIN/MAX
		if wb_min<WMIN: wb_min=WMIN
		if wb_max>WMAX: wb_max=WMAX
		#! Now get w crds
		wb_crd_min=h8.GetAxis(H8_DIM['W']).FindBin(wb_min+(WBINW/2))
		wb_crd_max=h8.GetAxis(H8_DIM['W']).FindBin(wb_max-(WBINW/2))

		#! Form objects to be returned
		q2wb_crds=[((q2b_crd_min,q2b_crd_max),(wb_crd_min,wb_crd_max))]
		q2wb_lmts=[((q2b_min,q2b_max),(wb_min,wb_max))]

	return q2wb_crds,q2wb_lmts

def get_q2bin(q2wbin):
	return q2wbin.split('_')[0]
def get_wbin(q2wbin):
	return q2wbin.split('_')[1]
def get_q2bin_le(q2wbin):
	return float(get_q2bin(q2wbin).split('-')[0])
def get_wbin_le(q2wbin):
	return float(get_wbin(q2wbin).split('-')[0])
def get_q2bin_ue(q2wbin):
	return float(get_q2bin(q2wbin).split('-')[1])
def get_wbin_ue(q2wbin):
	return float(get_wbin(q2wbin).split('-')[1])

def merge_hW(hW,hWm):
	'''
	+ hW is the histogram that will be merged into hWm
	'''
	#! etgt
	for ibin in range(hW.GetNbinsX()):
		#! get info from q2wb hist
		binle=hW.GetBinLowEdge(ibin+1)
		binc=hW.GetBinContent(ibin+1)
		bine=hW.GetBinError(ibin+1)
		#! now put it in total hist
		#! Note bins less than 1.400 will be underflow
		bin=hWm.FindBin(binle+(WBINW/2))
		hWm.SetBinContent(bin,binc)
		hWm.SetBinError(bin,bine)

def get_q2wbinl_from_file(f,dbg=False):
	q2wbinl=[]
	i=0 #! for dbg_bins
	for path,dirs,files in f.walk():
		if path=="":continue #! Avoid root path
		path_arr=path.split("/")
		if len(path_arr)==1:
			q2wbinl.append(path)
			i+=1
		if dbg==True:
			if i>=dbg_bins: break
	return q2wbinl

def plot_1D_athtcs():
		#ROOT.gStyle.Reset()
		#! Stats Box
		ROOT.gStyle.SetOptStat(0)

		# ROOT.gStyle.SetLabelSize(0.5,"t")
		# ROOT.gStyle.SetTitleSize(0.5,"t")
		#ROOT.gStyle.SetPaperSize(20,26);
		#ROOT.gStyle.SetPadTopMargin(0.15)#(0.05);
		#ROOT.gStyle.SetPadRightMargin(0.15)#(0.05);
		ROOT.gStyle.SetPadBottomMargin(0.20)#(0.16);
		#ROOT.gStyle.SetPadLeftMargin(0.20)#(0.12);

		ROOT.gStyle.SetTitleW(10)# //title width 
		ROOT.gStyle.SetTitleFontSize(20)# 
		ROOT.gStyle.SetTitleH(0.15)# //title height 
		
		#! + The following options do not seem to work from here
		#! + I have to set them in label_hist_obs1D()
		#ROOT.gStyle.SetTitleFont(42,"xyz")
		#ROOT.gStyle.SetTitleSize(.35,"xyz")
		#ROOT.gStyle.SetTitleOffset(0.5,"xyz");

		#!get rid of X error bars and y error bar caps
		ROOT.gStyle.SetErrorX(0.001)

def hist_1D_athtcs(h1,hW):
	'''
	+ h1(vst,var)           
	+ hW	
	'''
	#! First do aesthetics for h1CF
	for k in h1.keys():	
		vst,var=k[0],k[1]
		
		#! Setup the main, X and Y titles
		h1[k].SetTitle("")
		h1[k].SetXTitle( "%s[%s]"%(VAR_NAMES[(vst,var)],VAR_UNIT_NAMES[var]) )
		h1[k].SetYTitle("CF")
		#! X-axis title aesthetics
		h1[k].GetXaxis().SetTitleOffset(0.7)
		h1[k].GetXaxis().SetLabelSize(0.05)
		h1[k].GetXaxis().SetTitleSize(0.10)
		#! Y-axis title aesthetics
		h1[k].GetYaxis().SetTitleOffset(0.7)
		h1[k].GetYaxis().SetLabelSize(0.05)
		h1[k].GetYaxis().SetTitleSize(0.10)
	#! Now do hW
	if hW!=None:
		hW.SetTitle("")
		hW.SetXTitle("W (GeV)")
		hW.SetYTitle("CF")
		

def plot_CF(q2wbin,h1CF,hWCF,outdir):
	"""
	"""
	#print "In plot_CF"
	#! Transfrom q2wbin from a.aa-b.bb_x.xx-y.yy => [a,aa GeV^2,b.bb GeV^2) [x.xx GeV,y.yy GeV^2) GeV
	q2bin_le=get_q2bin_le(q2wbin)
	q2bin_ue=get_q2bin_ue(q2wbin)
	wbin_le=get_wbin_le(q2wbin)
	wbin_ue=get_wbin_ue(q2wbin)
	q2wbin_title="Q^{2},W bin = [%.2f GeV^{2}, %.2f GeV^{2}), [%.3f GeV, %.3f GeV)"%(q2bin_le,q2bin_ue,wbin_le,wbin_ue)

	#! Set up plotting aesthetics
	plot_1D_athtcs()
	hist_1D_athtcs(h1CF,hWCF)

	#! First plot h1CF 
	#! TCanvas's pad_map[pad,vst,var] defined as per Gleb's display
	pad_map=[(1,1,"M1"),   (2,3,'M2'),   (3,2,'M2'),
			 (4,1,"THETA"),(5,3,'THETA'),(6,2,'THETA'),
			 (7,1,"ALPHA"),(8,3,'ALPHA'),(9,2,'ALPHA')]
						  
	c1=ROOT.TCanvas("c","c",1000,1000)
	pad_t=ROOT.TPad("pad_l","Legend pad",0.05,0.935,0.95,1.00)
	pad_p=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
	pad_p.Draw()
	pad_t.Draw()
	pad_t.cd()
	pt=ROOT.TPaveText(.05,.1,.95,.8)
	pt.AddText(q2wbin_title)
	pt.SetTextSize(0.40)
	pt.Draw()
	pad_p.Divide(3,3)
	for item in pad_map:
		pad,vst,var=item[0],item[1],item[2]
		#print "pad,vst,var=",pad,vst,var
		gpad=pad_p.cd(pad)
		h1CF[vst,var].Draw()
	#! Save canvas 
	# q2bin=self.get_q2bin(q2wbin)
	# wbin=self.get_wbin(q2wbin)
	# outdir_q2bin=os.path.join(self.OUTDIR_1D,"q%s"%q2bin)
	# if not os.path.exists(outdir_q2bin):
	# 	os.makedirs(outdir_q2bin)
	c1.SaveAs("%s/c_h1CF.png"%(outdir))
	cmd=['convert',"%s/c_h1CF.png"%(outdir),"%s/c_h1CF.pdf"%(outdir)]
	t=subprocess.Popen(cmd,stderr=subprocess.STDOUT)
	t.wait()
	

	#! Now plot hWCF
	if hWCF!=None:
		cW=ROOT.TCanvas("c","c")
		pad_t=ROOT.TPad("pad_l","Legend pad",0.05,0.935,0.95,1.00)
		pad_p=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
		pad_p.Draw()
		pad_t.Draw()
		pad_t.cd()
		pt=ROOT.TPaveText(.05,.1,.95,.8)
		pt.AddText(q2wbin_title)
		pt.SetTextSize(0.40)
		pt.Draw()
		pad_p.cd()
		hWCF.Draw()
		cW.SaveAs("%s/c_hWCF.png"%(outdir))
		cmd=['convert',"%s/c_hWCF.png"%(outdir),"%s/c_hWCF.pdf"%(outdir)]
		t=subprocess.Popen(cmd,stderr=subprocess.STDOUT)
		t.wait()		

	#print "*** plot_1D() Done ***\n"
	return

#! main code
"""
For very Crs-W bin:
	1. Get h8(tgt,vst) per Crs-W from FIN[tgt] 
	2. For h8(tgt,vst), ensure H8_DIM['HEL'] is appropriately set
	3. For q2w bins (depending on Q2WBNG_OPT) in h8:
		#! First do h1
		3.i.   h8(tgt,vst)=>h5(tgt,vst)
		3.ii.  h5(tgt,vst)=>h1(q2w,tgt,vst,var) 
		3.iii. hCNTNF(q2w,tgt,vst,var)=( h1(q2w,etgt,vst,var)*R(ptgt) )/h1(q2w,ptgt,vst,var)
		3.iv.  hCF(q2w,tgt,vst,var) = 1-hCNTMNF(q2w,tgt,vst,var)
		#! Now do hW (note hCNTMNF,hCF for hW is obtained only in integrated results [see step 5] form)
		3.v.   h8(tgt,vst)=>hW(q2w,tgt,vst)
	4. Obtain CF (and related quantities) for h1 after integrating data over all q2w bins
	5. Obtain CF (and related quantities) for hW after integrating data over all q2w bins
		+ For hW, vst avg is also calculated

	6. Save all relevant plots in .png form
	   (Note that for full debugging information, see the .root file)
"""
#! Before beginning prepare fout
fout=ROOT.TFile(os.path.join(OUTDIR,"CF.root"),"RECREATE")

#! Create structures to store:
#!   + h1(q2wb,tgt,vst,var) 
#!   + hW(q2wb,tgt,vst,var)
#! + This needs to be done so that later I can get q2wb integrated values
h1=OrderedDict()
hW=OrderedDict()

#! Loop over Crw-W bins
for iw in range(q2w_bng.NBINS_WCRS):
	print "**** Processing Crs-W bin %s ***"%(iw+1)

	#! 1. Get h8(TGT,VST)
	h8=get_h8_seq_drct(iw)
	# #! Debug h8
	# for k in h8:
	# 	print k,h8[k].GetEntries()

	#! 2. Make sure H8_DIM['HEL'] is appropriately set for h8(tgt,vst)
	for vst in VSTS:
		#! Get HEL binning information from any of the TGTs
		nbins=h8['ptgt_lse',vst].GetAxis(H8_DIM['HEL']).GetNbins()
		for tgt in TGTS:
			h8[tgt,vst].GetAxis(H8_DIM['HEL']).SetRange(1,nbins)
								
	#! 3. Loop of Q2-W bins 
	#! + Get Q2-W binning information from any of the h8(SEQ_DRCT,VST)
	#! + Q2-W binning information (for more details see 'get_q2w_bng()'):
	#! 	+ q2wb_crd=(q2b_crd,wb_crd)
	#! 	+ q2wb_lmt=((q2b_min,q2b_max),(wb_min,wb_max))
	q2wb_crds,q2wb_lmts=get_q2w_bng(h8['ptgt_lse',1],Q2WBNG_OPT)
	#! Debug Q2-W binning
	print "q2wb_crd=(q2b_crd,wb_crd) q2wb_lmt=((q2b_min,q2b_max),(wb_min,wb_max)):"
	for q2wb_crd,q2wb_lmt in zip(q2wb_crds,q2wb_lmts):
		print q2wb_crd,q2wb_lmt	
	#sys.exit()


	for q2wb_crd,q2wb_lmt in zip(q2wb_crds,q2wb_lmts):
		q2b_min,q2b_max=q2wb_lmt[0][0],q2wb_lmt[0][1]
		wb_min, wb_max =q2wb_lmt[1][0],q2wb_lmt[1][1]
		q2wb_name="%0.2f-%0.2f_%0.3f-%0.3f"%(q2b_min,q2b_max,wb_min,wb_max)#(q2wb_lmt[0][0],q2wb_lmt[0][1],q2wb_lmt[1][0],q2wb_lmt[1][1])
		if DBG==True and (q2wb_name!=DBGQ2WBIN1 and q2wb_name!=DBGQ2WBIN2):
			continue
		if DBG==True: print "Debug mode. Processing q2wb_name=",q2wb_name
		else:         print "Processing q2wb_name=",q2wb_name
			
		#! 3.i. h8(TGT,VST)=>h5(TGT,VST)
		h5=OrderedDict()
		#! 3.i. h8(TGT,VST)=>h5(TGT,VST)
		for tgt in TGTS:
			for vst in VSTS:
				#! Make Q2,W projections
				if Q2WBNG_OPT=='all':
					h8[tgt,vst].GetAxis(H8_DIM['Q2']).SetRange(q2wb_crd[0],q2wb_crd[0])
					h8[tgt,vst].GetAxis(H8_DIM['W']).SetRange(q2wb_crd[1],q2wb_crd[1])
				elif Q2WBNG_OPT=='none':
					h8[tgt,vst].GetAxis(H8_DIM['Q2']).SetRange(q2wb_crd[0][0],q2wb_crd[0][1])
					h8[tgt,vst].GetAxis(H8_DIM['W']).SetRange(q2wb_crd[1][0],q2wb_crd[1][1])
				h5[tgt,vst]=h8[tgt,vst].Projection(5,H5_PROJDIM,"E")
				h5[tgt,vst].SetName('h5')
				h5[tgt,vst].SetTitle('%s_%s_%d'%(q2wb_name,tgt,vst))
				thntool.SetUnderOverFLowBinsToZero(h5[tgt,vst])

		#! 3.ii.  h5(TGT,VST)=>h1(TGT,VST,VAR) 
		#! + For technical reasons (see below), h1s are directly created 
		#!   under the appropriate file directory when saving the histograms 
		#! + Additionally, it so happens, that this also saves memory
		#!   since the h1(SEQ_ALL,VST,VARS) structure does not have to be created
		#! 
		#!   + Technical reason
		#!   ------------------
		#! 	+ Write() seems to be needed only for THnSparse Objects, where after 'cd' into appropriate
		#!    fout's directory, THnSparse::Write() needs to be called.
		#!  + TH1::Write()' seems to save all cycles of the object(name;cycle) making fout look messy.
		#!	  Therefore, instead of first creating h1(SEQ_ALL,VST,VARS) and then using 'Write()'
		#!    to save each of them, it is simpler to use the default behaviour of TH1s, where 
		#!    if they are created after 'cd' into the appropriate fout's directory, they are automatically
		#!    saved when fout is written.

		#3.iii.  h5(TGT,VST) => h1(TGT,VST,VARS) (note h1 includes hW)
		#h1=OrderedDict()			
		q2wb_dir=fout.GetDirectory(q2wb_name)
		if q2wb_dir==None: 
			q2wb_dir=fout.mkdir(q2wb_name)
		q2wb_dir.cd() 	
		for item in list(itertools.product(TGTS,VSTS)):
			tgt,vst=item[0],item[1]
			#! Ensure seq_dir exists and cd into it
			tgt_dir=q2wb_dir.GetDirectory(tgt)
			if tgt_dir==None:
				tgt_dir=q2wb_dir.mkdir(tgt)
			tgt_dir.cd()
			#! Ensure vst_dir exists and cd into it
			vst_name="VST%d"%vst
			vst_dir=tgt_dir.GetDirectory(vst_name)
			if vst_dir==None:
				vst_dir=tgt_dir.mkdir(vst_name)
			vst_dir.cd()
			#! Save h5[seq,vst]
			h5[tgt,vst].Write()
			#! Create and save h1 for every [seq,vst,var]
			for var in VARS:
				h1[q2wb_name,tgt,vst,var]=h5[tgt,vst].Projection(H5_DIM[var],"E")
				h1[q2wb_name,tgt,vst,var].SetName('h1_%s'%var)
				h1[q2wb_name,tgt,vst,var].SetTitle("%s_%s_%d_%s"%(q2wb_name,tgt,vst,var))
				if var=='M1' or var=='M2':
					setM1M2axisrange(h1[q2wb_name,tgt,vst,var],vst,var,wb_max)
				#h1[q2wb_name,tgt,vst,var].Write()

		#! The following steps accmplished in this next loop
		#!  3.iii. hCNTNF('ptgt_lse/tgt',VST,VAR)=( h1(ETGT,VST,VAR)*R(PTGT) )/h1(PTGT,VST,VAR)
		#! 	3.iv.  hCF(VST,VAR) = 1-hCNTMNF(VST,VAR)
		#! Create required structures
		hCNTMNF=OrderedDict()
		hCF=OrderedDict()
		hetgt,hptgt=OrderedDict(),OrderedDict()
		hist1=OrderedDict()
		q2wb_dir.cd() 
		#! Note that 'etgt' in the following loop is purposefully omitted	
		for item in list(itertools.product(['ptgt_lse','ptgt_tgt'],VSTS,VARS)):
			tgt,vst,var=item[0],item[1],item[2]
			#! Make sure CF_dir/vst_dir exists
			#! CF_dir
			CF_dir=q2wb_dir.GetDirectory("CF_%s"%tgt)
			if CF_dir==None:
				CF_dir=q2wb_dir.mkdir("CF_%s"%tgt)
			CF_dir.cd()
			#! CF_dir/vst_dir
			vst_name="VST%d"%vst
			vst_dir=CF_dir.GetDirectory(vst_name)
			if vst_dir==None:
				vst_dir=CF_dir.mkdir(vst_name)
			vst_dir.cd()
			#! First make clones of h1['etgt',vst,var] and h1[tgt,vst,var] so that the originals 
			#! are not affected by Scaling and other hist operations
			#! + Call their Sumw2()
			hetgt[q2wb_name,tgt,vst,var]=h1[q2wb_name,'etgt',vst,var].Clone('hetgt_scaled_%s'%var)
			hetgt[q2wb_name,tgt,vst,var].SetTitle("%s_%s_%d_%s"%(q2wb_name,'hetgt_scaled',vst,var))
			hetgt[q2wb_name,tgt,vst,var].Sumw2()
			hptgt[q2wb_name,tgt,vst,var]=h1[q2wb_name,tgt,vst,var].Clone("hptgt_%s"%(var))
			hptgt[q2wb_name,tgt,vst,var].Sumw2()
			#! Cloned histograms have to explicity Written? hptgt[q2wb_name,tgt,vst,var].SetDirectory(ROOT.gROOT.CurrentDirectory())
			hptgt[q2wb_name,tgt,vst,var].Write() #! Duplicate, but valuable for efficient debugging
			#! Now scale hetgt by R[ptgt_lse/tgt]
			hetgt[q2wb_name,tgt,vst,var].Scale(R[tgt])
			hetgt[q2wb_name,tgt,vst,var].Write()
			#hetgt.Write()
			#! Now obtain hCNTMNF
			hCNTMNF[q2wb_name,tgt,vst,var]=hetgt[q2wb_name,tgt,vst,var].Clone('hCNTMNF_%s'%var)
			hCNTMNF[q2wb_name,tgt,vst,var].SetTitle("%s_%s_%d_%s"%(q2wb_name,'hCNTMNF',vst,var))
			hCNTMNF[q2wb_name,tgt,vst,var].Sumw2()
			hCNTMNF[q2wb_name,tgt,vst,var].Divide(hptgt[q2wb_name,tgt,vst,var])
			hCNTMNF[q2wb_name,tgt,vst,var].SetMinimum(0)
			hCNTMNF[q2wb_name,tgt,vst,var].Write()
			#! Now obtain hCF
			#! First obtain hist-1
			hist1[q2wb_name,tgt,vst,var]=hptgt[q2wb_name,tgt,vst,var].Clone("hist1_%s"%var)
			hist1[q2wb_name,tgt,vst,var].Reset()
			for ibin in range(hptgt[q2wb_name,tgt,vst,var].GetNbinsX()):
				hist1[q2wb_name,tgt,vst,var].SetBinContent(ibin+1,1)
				hist1[q2wb_name,tgt,vst,var].SetBinError(ibin+1,0)
			hist1[q2wb_name,tgt,vst,var].Sumw2()
			hist1[q2wb_name,tgt,vst,var].Write()
			#! Now do hCF=1-hCNTMNF
			hCF[q2wb_name,tgt,vst,var]=hist1[q2wb_name,tgt,vst,var].Clone('hCF_%s'%var)
			hCF[q2wb_name,tgt,vst,var].SetTitle("%s_%s_%d_%s"%(q2wb_name,'hCF',vst,var))
			hCF[q2wb_name,tgt,vst,var].Sumw2()
			hCF[q2wb_name,tgt,vst,var].Add(hCNTMNF[q2wb_name,tgt,vst,var],-1)
			hCF[q2wb_name,tgt,vst,var].SetMinimum(0)
			hCF[q2wb_name,tgt,vst,var].Write()

		#! 3.v.   h8(TGT,VST)=>hW(TGT,VST)
		#! + hW not obtained for each VAR since it is already shown the the integral
		#!   for all VARs within a VST is the same
		q2wb_dir.cd() 	
		for item in list(itertools.product(TGTS,VSTS)):
			tgt,vst=item[0],item[1]
			#! Ensure seq_dir exists and cd into it
			tgt_dir=q2wb_dir.GetDirectory(tgt)
			if tgt_dir==None:
				tgt_dir=q2wb_dir.mkdir(tgt)
			tgt_dir.cd()
			#! Ensure vst_dir exists and cd into it
			vst_name="VST%d"%vst
			vst_dir=tgt_dir.GetDirectory(vst_name)
			if vst_dir==None:
				vst_dir=tgt_dir.mkdir(vst_name)
			vst_dir.cd()
			#! h8->hW
			hW[q2wb_name,tgt,vst]=h8[tgt,vst].Projection(H8_DIM['W'],"E")
			hW[q2wb_name,tgt,vst].SetName('hW')
			hW[q2wb_name,tgt,vst].SetTitle("%s_%s_%d"%(q2wb_name,tgt,vst))

	print "Done processing Crs-W bin %s\n"%(iw+1)		

#! 4. Obtain q2w integrated quantities for h1
#! First make sure that h1(q2wb,tgt,vst,var) has 6*3*3*5=270 objects in it.
print "*** Going to compute integrated q2wb CF for h1 ***" 
print "Total number of objects in h1(q2wb,tgt,vst,var)=%d (There should be 6*3*3*5=270)"%len(h1)
#! Get all q2wb in h1
q2wbl=[]
for k in h1:
	q2wbl.append(k[0])
#! Keep only unique values
q2wbl=list(set(q2wbl))
print "q2wbs in h1=%s (total=%d)"%(q2wbl,len(q2wbl))
#! Now start computing integrated values
q2wb_name="%0.2f-%0.2f_%0.3f-%0.3f_itg"%(Q2MIN,Q2MAX,WMIN,WMAX)
q2wb_dir=fout.GetDirectory(q2wb_name)
if q2wb_dir==None: 
	q2wb_dir=fout.mkdir(q2wb_name)
q2wb_dir.cd() 
#! Create required structures
h1CNTMNFav=OrderedDict()
h1CFav=OrderedDict()
h1etgt,h1ptgt=OrderedDict(),OrderedDict()
hist11=OrderedDict()
for item in list(itertools.product(['ptgt_lse','ptgt_tgt'],VSTS,VARS)):
	tgt,vst,var=item[0],item[1],item[2]
	#! Make sure CF_dir/vst_dir exists
	#! CF_dir
	CF_dir=q2wb_dir.GetDirectory("CF_%s"%tgt)
	if CF_dir==None:
		CF_dir=q2wb_dir.mkdir("CF_%s"%tgt)
	CF_dir.cd()
	#! CF_dir/vst_dir
	vst_name="VST%d"%vst
	vst_dir=CF_dir.GetDirectory(vst_name)
	if vst_dir==None:
		vst_dir=CF_dir.mkdir(vst_name)
	vst_dir.cd()
	#! First sum all h1['etgt',vst,var] and h1[tgt,vst,var] 
	#! + Call their Sumw2()
	for i,q2wb in enumerate(q2wbl):
		if i==0: #! create by cloning
			h1etgt[tgt,vst,var]=h1[q2wb,'etgt',vst,var].Clone('hetgt_scaled_%s'%var)
			h1etgt[tgt,vst,var].SetTitle("%s_%s_%d_%s"%(q2wb_name,'hetgt_scaled',vst,var))
			h1etgt[tgt,vst,var].Sumw2()
			h1ptgt[tgt,vst,var]=h1[q2wb,tgt,vst,var].Clone('hptgt_%s'%var)
			h1ptgt[tgt,vst,var].SetTitle("%s_%s_%d_%s"%(q2wb_name,'hptgt_scaled',vst,var))
			h1ptgt[tgt,vst,var].Sumw2()
		#else: #! add
			h1etgt[tgt,vst,var].Add(h1[q2wb,'etgt',vst,var])
			h1ptgt[tgt,vst,var].Add(h1[q2wb,tgt,vst,var])
	#! Now scale hetgt by R[ptgt_lse/tgt]
	h1etgt[tgt,vst,var].Scale(R[tgt])
	#! Write
	h1etgt[tgt,vst,var].Write()
	h1ptgt[tgt,vst,var].Write()
	#! Now obtain hCNTMNF
	h1CNTMNFav[tgt,vst,var]=h1etgt[tgt,vst,var].Clone('hCNTMNF_%s'%var)
	h1CNTMNFav[tgt,vst,var].SetTitle("%s_%s_%d_%s"%(q2wb_name,'hCNTMNF',vst,var))
	h1CNTMNFav[tgt,vst,var].Sumw2()
	h1CNTMNFav[tgt,vst,var].Divide(h1ptgt[tgt,vst,var])
	h1CNTMNFav[tgt,vst,var].SetMinimum(0)
	h1CNTMNFav[tgt,vst,var].Write()
	#! Now obtain hCF
	#! First obtain hist-1
	hist11[tgt,vst,var]=h1ptgt[tgt,vst,var].Clone('hist1_%s'%var)
	hist11[tgt,vst,var].Reset()
	for ibin in range(h1ptgt[tgt,vst,var].GetNbinsX()):
		hist11[tgt,vst,var].SetBinContent(ibin+1,1)
		hist11[tgt,vst,var].SetBinError(ibin+1,0)
	hist11[tgt,vst,var].Sumw2()
	hist11[tgt,vst,var].Write()
	#! Now do hCF=1-hCNTMNF
	h1CFav[tgt,vst,var]=hist11[tgt,vst,var].Clone('hCF_%s'%var)
	h1CFav[tgt,vst,var].SetTitle("%s_%s_%d_%s"%(q2wb_name,'hCF',vst,var))
	h1CFav[tgt,vst,var].Sumw2()
	h1CFav[tgt,vst,var].Add(h1CNTMNFav[tgt,vst,var],-1)
	h1CFav[tgt,vst,var].SetMinimum(0)
	h1CFav[tgt,vst,var].Write()

#! 5. Obtain q2w integrated quantities for hW
#! First make sure that hW(q2wb,tgt,vst) has 6*3*3=54 objects in it.
print "*** Going to compute integrated q2wb integrated CF for hW ***" 
print "Total number of objects in hW(q2wb,tgt,vst)=%d (There should be 6*3*3=54)"%len(hW)
#! Get all q2wb in h1
q2wbl=[]
for k in hW:
	q2wbl.append(k[0])
#! Keep only unique values
q2wbl=list(set(q2wbl))
print "q2wbs in hW=%s (total=%d)"%(q2wbl,len(q2wbl))

#! Now compute integrated values
q2wb_name="%0.2f-%0.2f_%0.3f-%0.3f_itg"%(Q2MIN,Q2MAX,WMIN,WMAX)
q2wb_dir=fout.GetDirectory(q2wb_name)
if q2wb_dir==None: 
	q2wb_dir=fout.mkdir(q2wb_name)
q2wb_dir.cd() 
#! Create required structures
hWCNTMNFav=OrderedDict()
hWCFav=OrderedDict()
hWetgt,hWptgt=OrderedDict(),OrderedDict()
histW1=OrderedDict()
for item in list(itertools.product(['ptgt_lse','ptgt_tgt'],VSTS)):
	tgt,vst=item[0],item[1]
	#! Make sure CF_dir/vst_dir exists
	#! CF_dir
	CF_dir=q2wb_dir.GetDirectory("CF_%s"%tgt)
	if CF_dir==None:
		CF_dir=q2wb_dir.mkdir("CF_%s"%tgt)
	CF_dir.cd()
	#! CF_dir/vst_dir
	vst_name="VST%d"%vst
	vst_dir=CF_dir.GetDirectory(vst_name)
	if vst_dir==None:
		vst_dir=CF_dir.mkdir(vst_name)
	vst_dir.cd()
	#! First "merge" all hW['etgt',vst] and hW[tgt,vst] 
	#! + Call their Sumw2()
	#! + Binning for this merged hW is per ana2pi's final W limits and binning:29 bins,W=1.400 to 2.125
	#! + Note that entry for bins less than 1.400 will not be a problem as they will go into the underflow
	hWetgt[tgt,vst]=ROOT.TH1F("hetgt_W","hetgt_W",29,1.400,2.125)
	hWetgt[tgt,vst].SetName('hetgt_scaled_W')
	hWetgt[tgt,vst].SetTitle("%s_%s_%d"%(q2wb_name,'hetgt_scaled_W',vst))
	hWptgt[tgt,vst]=ROOT.TH1F("hptgt_W","hptgt_W",29,1.400,2.125)
	hWptgt[tgt,vst].SetName('hptgt_W')
	hWptgt[tgt,vst].SetTitle("%s_%s_%d"%(q2wb_name,'hptgt_W',vst))
	for i,q2wb in enumerate(q2wbl):
		merge_hW(hW[q2wb,'etgt',vst],hWetgt[tgt,vst])
		merge_hW(hW[q2wb,tgt,vst],hWptgt[tgt,vst])
	hWetgt[tgt,vst].Sumw2()
	hWptgt[tgt,vst].Sumw2()
	#! Now scale hetgt by R[ptgt_lse/tgt]
	hWetgt[tgt,vst].Scale(R[tgt])
	#hWetgt[tgt,vst].Write()
	#! Now obtain hCNTMNF
	hWCNTMNFav[tgt,vst]=hWetgt[tgt,vst].Clone('hCNTMNF_W')
	hWCNTMNFav[tgt,vst].SetTitle("%s_%s_%d_%s"%(q2wb_name,'hCNTMNF',vst,'W'))
	hWCNTMNFav[tgt,vst].Sumw2()
	hWCNTMNFav[tgt,vst].Divide(hWptgt[tgt,vst])
	hWCNTMNFav[tgt,vst].SetMinimum(0)
	#hWCNTMNFav[tgt,vst].Write()
	#! Now obtain hCF
	#! First obtain hist-1
	histW1[tgt,vst]=hWptgt[tgt,vst].Clone("hist1_W")
	histW1[tgt,vst].Reset()
	for ibin in range(hWptgt[tgt,vst].GetNbinsX()):
		histW1[tgt,vst].SetBinContent(ibin+1,1)
		histW1[tgt,vst].SetBinError(ibin+1,0)
	histW1[tgt,vst].Sumw2()
	#! Now do hCF=1-hCNTMNF
	hWCFav[tgt,vst]=histW1[tgt,vst].Clone('hCF_W')
	hWCFav[tgt,vst].SetTitle("%s_%s_%d_%s"%(q2wb_name,'hCF',vst,'W'))
	hWCFav[tgt,vst].Sumw2()
	hWCFav[tgt,vst].Add(hWCNTMNFav[tgt,vst],-1)
	hWCFav[tgt,vst].SetMinimum(0)
	#hWCFav[tgt,vst].Write()
#! Compute vst average for hW for each ptgt
hWCFav_vst=OrderedDict()
for tgt in ['ptgt_lse','ptgt_tgt']:
	CF_dir=q2wb_dir.GetDirectory("CF_%s"%tgt)
	if CF_dir==None:
		CF_dir=q2wb_dir.mkdir("CF_%s"%tgt)
	CF_dir.cd()
	#! Compute vst average
	hWCFav_vst[tgt]=hWCFav[tgt,1].Clone('hCF_W_vst_avg')
	hWCFav_vst[tgt].Add(hWCFav[tgt,2])
	hWCFav_vst[tgt].Add(hWCFav[tgt,3])
	hWCFav_vst[tgt].Scale(1/3)
	hWCFav_vst[tgt].SetMinimum(0)
	hWCFav_vst[tgt].Write()

#! Write file
fout.Write()
fout.Close()

#! 6. Save all relevant plots in .png form
outdir_png=os.path.join(OUTDIR,"CF_PNG")
if not os.path.exists(outdir_png):
	os.makedirs(outdir_png)
#! Open .root file to get plots from
f=root_open(os.path.join(OUTDIR,'CF.root'))
#! For every CF*q2wbin, get h1CF and hWCF, and plot them
h1CF=OrderedDict()
#! Not needed since here is only 1 per CF*q2w hWCF=OrderedDict()
for cf in ['CF_ptgt_lse','CF_ptgt_tgt']:
	outdir_png_cf=os.path.join(outdir_png,cf)
	if not os.path.exists(outdir_png_cf):
		os.makedirs(outdir_png_cf)
	#! Get q2wbinl from file and loop over it
	q2wbinl=get_q2wbinl_from_file(f)
	for q2wbin in q2wbinl:
		outdir_png_cf_q2wbin=os.path.join(outdir_png_cf,q2wbin)
		if not os.path.exists(outdir_png_cf_q2wbin):
			os.makedirs(outdir_png_cf_q2wbin)
		for vst,var in itertools.product(VSTS,VARS):
			h1CF[vst,var]=f.Get('%s/%s/VST%d/hCF_%s'%(q2wbin,cf,vst,var))
		#! get hWCF only for integrated q2wbin
		hWCF=None
		if "itg" in q2wbin:
			hWCF=f.Get('%s/%s/hCF_W_vst_avg'%(q2wbin,cf))  
		#! plot
		plot_CF(q2wbin,h1CF,hWCF,outdir_png_cf_q2wbin)

print "If appearing to be stuck, then either fout is still being written or Python is probably doing \"garbage collection\"(?); Wait a while!"

	

