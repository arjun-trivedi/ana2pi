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
	(dsigma/dX_ij)_corr = (dsigma/dX_ij) * (CF_ij), where
	+ CF_ij = 1 - CNTMNF_ij
	+ i = W, 9 1D xsecs
	+ j = bin index
	+ CF_ij =  (N-etgt_ij * R)/(N-ptgt_ij)

+ Note that CF_ij is returned as a histogram object
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
WMIN,WMAX=1.425,2.125

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
#! + from $STUDY_TGT_BG_E16_DATADIR/results_R_052917/R.txt
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
		wb_crd_min=h8.GetAxis(H8_DIM['W']).FindBin(wb_min)
		wb_crd_max=h8.GetAxis(H8_DIM['W']).FindBin(wb_max)

		#! Form objects to be returned
		q2wb_crds=[((q2b_crd_min,q2b_crd_max),(wb_crd_min,wb_crd_min))]
		q2wb_lmts=[((q2b_min,q2b_max),(wb_min,wb_max))]

	return q2wb_crds,q2wb_lmts

#! main code
"""
For very Crs-W bin:
	1. Get h8(TGT,VST) per Crs-W from FIN[TGT] 
	2. For h8(TGT,VST), ensure H8_DIM['HEL'] is appropriately set
	3. For every Q2-W bin in h8:
		3.i.   h8(TGT,VST)=>h5(TGT,VST)
		3.ii.  h5(TGT,VST)=>h1(TGT,VST,VAR) (Note h1 includes hW)
		3.iii. hCNTNF(VST,VAR)=( h1(ETGT,VST,VAR)*R(PTGT) )/h1(PTGT,VST,VAR)
		3.iv.  hCF(VST,VAR) = 1-hCNTMNF(VST,VAR)
		3.v.   Save hCF(VST,VAR)
"""
#! Before beginning prepare fout
fout=ROOT.TFile(os.path.join(OUTDIR,"CF.root"),"RECREATE")

#! Create structure to store h1(q2wb,tgt,vst,var) 
#! + This needs to be done so that later I can get q2wb averaged values
h1=OrderedDict()

#! Loop over Crw-W bins
for iw in range(q2w_bng.NBINS_WCRS):
	print "Processing Crs-W bin %s"%(iw+1)

	#! 1. Get h8(TGT,VST)
	h8=get_h8_seq_drct(iw)
	# #! Debug h8
	# for k in h8:
	# 	print k,h8[k].GetEntries()

	#! 2. Make sure H8_DIM['HEL'] is appropriately set for h8(SEQ_DRCT,vst)
	for vst in VSTS:
		nbins=h8['ptgt_lse',vst].GetAxis(H8_DIM['HEL']).GetNbins()
		h8['ptgt_lse',vst].GetAxis(H8_DIM['HEL']).SetRange(1,nbins)
		h8['ptgt_tgt',vst].GetAxis(H8_DIM['HEL']).SetRange(1,nbins)
		h8['etgt',vst].GetAxis(H8_DIM['HEL']).SetRange(1,nbins)
							
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
		
		#! 3.ii.  h5(TGT,VST)=>h1(TGT,VST,VAR) (Note h1 includes hW)
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

		#! The following steps accmplished in this next loop
		#!  3.iii. hCNTNF('ptgt_lse/tgt',VST,VAR)=( h1(ETGT,VST,VAR)*R(PTGT) )/h1(PTGT,VST,VAR)
		#! 	3.iv.  hCF(VST,VAR) = 1-hCNTMNF(VST,VAR)
		#! 	3.v.   Save hCF(VST,VAR): Done simply by creating histograms in steps 3.iii and 3.iv
		#! Create required structures
		hCNTMNF=OrderedDict()
		hCF=OrderedDict()
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
			hetgt=h1[q2wb_name,'etgt',vst,var].Clone("hetgt")
			hetgt.Sumw2()
			hpgt =h1[q2wb_name,tgt,vst,var].Clone("%s"%tgt)
			hpgt.Sumw2()
			#! Now scale hetgt by R[ptgt_lse/tgt]
			hetgt.Scale(R[tgt])
			hetgt.SetName('hetgt_scaled_%s'%var)
			hetgt.SetTitle("%s_%s_%d_%s"%(q2wb_name,'hetgt_scaled',vst,var))
			hetgt.Write()
			#! Now obtain hCNTMNF
			hCNTMNF[tgt,vst,var]=hetgt.Clone("hCNTMNF")
			hCNTMNF[tgt,vst,var].Sumw2()
			hCNTMNF[tgt,vst,var].Divide(hpgt)
			hCNTMNF[tgt,vst,var].SetMinimum(0)
			hCNTMNF[tgt,vst,var].SetName('hCNTMNF_%s'%var)
			hCNTMNF[tgt,vst,var].SetTitle("%s_%s_%d_%s"%(q2wb_name,'hCNTMNF',vst,var))
			hCNTMNF[tgt,vst,var].Write()
			#! Now obtain hCF
			#! First obtain hist-1
			hist1=hpgt.Clone("hist1")
			hist1.Reset()
			for ibin in range(hpgt.GetNbinsX()):
				hist1.SetBinContent(ibin+1,1)
				hist1.SetBinError(ibin+1,0)
			hist1.Sumw2()
			#! Now do hCF=1-hCNTMNF
			hCF[tgt,vst,var]=hist1.Clone("hCF")
			hCF[tgt,vst,var].Sumw2()
			hCF[tgt,vst,var].Add(hCNTMNF[tgt,vst,var],-1)
			hCF[tgt,vst,var].SetMinimum(0)
			hCF[tgt,vst,var].SetName('hCF_%s'%var)
			hCF[tgt,vst,var].SetTitle("%s_%s_%d_%s"%(q2wb_name,'hCF',vst,var))
			hCF[tgt,vst,var].Write()
	print "Done processing Crs-W bin %s"%(iw+1)		

#! Now compute q2wb averages CF
#! First make sure that h1(q2wb,tgt,vst,var) has 6*3*3*5=270 objects in it.
print "Going to compute averaged q2wb averaged CF" 
print "Total number of objects in h1(q2wb,tgt,vst,var)=%d (There should be 6*3*3*5=270)"%len(h1)
#! Get all q2wb in h1
q2wbl=[]
for k in h1:
	q2wbl.append(k[0])
#! Keep only unique values
q2wbl=list(set(q2wbl))
print "q2wbs in h1=%s (total=%d)"%(q2wbl,len(q2wbl))

#! Now compute averaged values
q2wb_name="%0.2f-%0.2f_%0.3f-%0.3f_avg"%(Q2MIN,Q2MAX,WMIN,WMAX)
q2wb_dir=fout.GetDirectory(q2wb_name)
if q2wb_dir==None: 
	q2wb_dir=fout.mkdir(q2wb_name)
q2wb_dir.cd() 
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
			hetgt=h1[q2wb,'etgt',vst,var].Clone("hetgt")
			hetgt.Sumw2()
			hpgt=h1[q2wb,tgt,vst,var].Clone("%s"%tgt)
			hpgt.Sumw2()
		else: #! add
			hetgt.Add(h1[q2wb,'etgt',vst,var])
			hpgt.Add(h1[q2wb,tgt,vst,var])
	#! Now scale hetgt by R[ptgt_lse/tgt]
	hetgt.Scale(R[tgt])
	hetgt.SetName('hetgt_scaled_%s'%var)
	hetgt.SetTitle("%s_%s_%d_%s"%(q2wb_name,'hetgt_scaled',vst,var))
	hetgt.Write()
	#! Now obtain hCNTMNF
	hCNTMNF[tgt,vst,var]=hetgt.Clone("hCNTMNF")
	hCNTMNF[tgt,vst,var].Sumw2()
	hCNTMNF[tgt,vst,var].Divide(hpgt)
	hCNTMNF[tgt,vst,var].SetMinimum(0)
	hCNTMNF[tgt,vst,var].SetName('hCNTMNF_%s'%var)
	hCNTMNF[tgt,vst,var].SetTitle("%s_%s_%d_%s"%(q2wb_name,'hCNTMNF',vst,var))
	hCNTMNF[tgt,vst,var].Write()
	#! Now obtain hCF
	#! First obtain hist-1
	hist1=hpgt.Clone("hist1")
	hist1.Reset()
	for ibin in range(hpgt.GetNbinsX()):
		hist1.SetBinContent(ibin+1,1)
		hist1.SetBinError(ibin+1,0)
	hist1.Sumw2()
	#! Now do hCF=1-hCNTMNF
	hCF[tgt,vst,var]=hist1.Clone("hCF")
	hCF[tgt,vst,var].Sumw2()
	hCF[tgt,vst,var].Add(hCNTMNF[tgt,vst,var],-1)
	hCF[tgt,vst,var].SetMinimum(0)
	hCF[tgt,vst,var].SetName('hCF_%s'%var)
	hCF[tgt,vst,var].SetTitle("%s_%s_%d_%s"%(q2wb_name,'hCF',vst,var))
	hCF[tgt,vst,var].Write()
#! Write file
fout.Write()
fout.Close()
print "If appearing to be stuck, then either fout is still being written or Python is probably doing \"garbage collection\"(?); Wait a while!"

	

