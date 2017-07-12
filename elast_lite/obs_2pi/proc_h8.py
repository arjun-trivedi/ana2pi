#! Imports
from __future__ import division
from math import *

from array import *
from collections import OrderedDict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools

from multiprocessing import Process, Queue

#PyROOT
import ROOT

import os,sys

import atlib as atlib
import q2w_bng

ROOT.gROOT.ProcessLine(".L THnTool.C+")
from ROOT import THnTool

# Tools 
thntool=THnTool()

#Constants
MASS_P=0.93827203
MASS_PIP=0.13957018
MASS_PIM=0.13957018

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

SEQ_ALL=['ST','SR','SA','SC','SH','SF','ER','EC','EH','EF']
SEQ_DRCT=['ST','SR','ER']
SEQ_CALC=['SA','SC','SH','SF','EC','EH','EF']

#! Setup data to perform Empty Target Background subtraction (etgt BG sub)
#! [06-09-17] Currently, R['ptgt_lse']=13.67 is being used to scale etgt events
R=OrderedDict()
#! + from $STUDY_TGT_BG_E16_DATADIR/results_R_060917/R.txt
R['ptgt_lse']=13.67
R['ptgt_tgt']=13.40

class ProcH8:
	"""
	+ Accomplishes d2pi.root(ER,ST,SR) -> yield.root  
		+ if self.USEHEL=false: d2pi.root -> yield.root
		+ if self.USEHEL=true:  d2pi.root -> yield_hel.root 

	+ The direct interface for the user is the 'execute()' method that calls 'proc()' for h8(VST,SEQ) in each Crs-W bin. 
		+ The 'execute()' method uses Python's multiprocessing module to accomplish the task of processing each h8(VST,SEQ) in 
	a separate thread, thereby not having to load all h8(VST,SEQ) from every Crs-W bin at once and exhausting the system's memory (gothe14=12GB RAM)
		+ For more technical details, see documentation of 'execute()'

	+ For each Q2-W bin (Please see details of proc() for details):
		+ if self.USEHEL=false: h8(SEQ_DRCT,VST) => h5(SEQ_DRCT,VST) => h5(SEQ_CALC,VST) => h1(SEQ_ALL,VST)
		+ if self.USEHEL=true:  Not implemented yet!
	
	"""
	def __init__(self,obsdir,simnum='siml',q2min=1.25,q2max=5.25,wmin=1.400,wmax=2.125,
		         acc_rel_err_cut=-1,
		         usehel=False,do_etgt_BG_sub=True,dbg=False):
		self.SIMNUM=simnum

		self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX=q2min,q2max,wmin,wmax
		print "Q2MIN,Q2MAX,WMIN,WMAX=",self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX

		self.ACC_REL_ERR_CUT=acc_rel_err_cut
		print "ACC_REL_ERR_CUT=",self.ACC_REL_ERR_CUT

		self.USEHEL=usehel
		print "USEHEL=",self.USEHEL

		self.DO_ETGT_BG_SUB=do_etgt_BG_sub
		print "DO_ETGT_BG_SUB=",self.DO_ETGT_BG_SUB

		self.DBG=dbg
		print "DBG=",self.DBG

		self.DATADIR=obsdir

		#! Prepare FIN
		self.FIN=OrderedDict()
		self.FIN['ER']=ROOT.TFile(os.path.join(self.DATADIR,'d2pi_exp','d2piR.root'))
		self.FIN['ST']=ROOT.TFile(os.path.join(self.DATADIR,'d2pi_sim',self.SIMNUM,'d2piT.root'))
		self.FIN['SR']=ROOT.TFile(os.path.join(self.DATADIR,'d2pi_sim',self.SIMNUM,'d2piR.root'))
		if self.DO_ETGT_BG_SUB: #! then get file with etgt d2piR data
			self.FIN['ET']=ROOT.TFile(os.path.join(os.environ['D2PIDIR_EXP_E16'],'data_tgt_051817','d2piR_etgt_cutruns.root'))
		#! Prepare OUTDIR
		self.OUTDIR=os.path.join(self.DATADIR,self.SIMNUM)
		if not os.path.exists(self.OUTDIR): 
			os.makedirs(self.OUTDIR)
		#! + Prepare FOUTNAME 
		#! + Note only FOUTNAME is created here; FOUT creation is handled, with attention to 
		#!   using multiprocessing, by 'proc()'. See doc. of 'proc()' for more details
		if self.USEHEL: self.FOUTNAME="yield_hel.root"
		else:           self.FOUTNAME="yield.root"
		print "DATADIR=%s\nFIN[ER]=%s\nFIN[ST]=%s\nFIN[SR]=%s\nOUTDIR=%s\nFOUTNAME=%s"%(self.DATADIR, self.FIN['ER'].GetName(),self.FIN['ST'].GetName(),self.FIN['SR'].GetName(), self.OUTDIR, self.FOUTNAME)
		

	def execute(self):
		'''
		+ The 'execute()' method uses Python's multiprocessing module to accomplish the task of processing
		  each h8(VST,SEQ) in a separate thread, thereby not having to load all h8(VST,SEQ) from every Crs-W bin
		  at once and exhausting the system's memory (gothe14=12GB RAM)
		+ Every Crs-W bin is seqentially processed in its own independent thread. This is accomplished by:
			+ Using multiprocessing.Process to call 'proc(iw,que)' for each Crs-W bin
			+ Waiting for each thread to complete before going on to the next bin.
				+ 'proc(iw,que)' puts in the que:
					+  0: if it returns at the end of the function (=Sucess)
					+ -1: if it returns some place earlier (=Failure)
		'''
		#! Loop over Crw-W bins

		queue=Queue()
		for iw in range(q2w_bng.NBINS_WCRS):
			print "*** execute(): Going to call proc(iw=%d) ***"%iw
			if self.DBG:
				if iw>2: break
			p=Process(target=self.proc, args=(iw,queue))
			p.start()
			p.join() # this blocks until the process terminates
			res=queue.get()
			print "*** execute(): Result from proc(iw=%d)=%d ***"%(iw,res)
		# self.FOUT.Write()
		# self.FOUT.Close()
		print "execute():Done"
		print "If the progam is not terminating, then Python is probably doing \"garbage collection\"(?); Wait a while!"
	
	def proc(self,iw,que=None):
		"""
		1. Get h8(SEQ_DRCT,VST) per Crs-W from FIN[SEQ_DRCT] 
		2. For h8(SEQ_DRCT,VST), ensure H8_DIM['HEL'] is appropriately set
		3. Open fout (self.FOUTNAME) in "appropriate" mode 
			+ if 1st Crs-W bin, then "RECREATE", else "UPDATE"
		4. For every Q2-W bin in h8:
			3.i.   h8(SEQ_DRCT,VST)=>h5(SEQ_DRCT,VST)
			3.ii.  h5(SEQ_DRCT,VST)=>h5(SEQ_CALC,VST)
			3.iii. If self.USEHEL=false, then h5(SEQ_ALL,VST)=>h1(SEQ_ALL,VST,VARS) 
			3.iv.  Save h5(SEQ_ALL,VST) and h1(SEQ_ALL,VST,VARS)
				
		"""
		print "proc(): Processing Crs-W bin %s"%(iw+1)

		#! 1. Get h8(SEQ_DRCT,VST)
		h8=self.get_h8_seq_drct(iw)
		# #! Debug h8
		# for k in h8:
		# 	print k,h8[k].GetEntries()

		#! 2. Make sure H8_DIM['HEL'] is appropriately set for h8(SEQ_DRCT,vst)
		for vst in VSTS:
			if self.USEHEL==False:
				nbins=h8['ST',vst].GetAxis(H8_DIM['HEL']).GetNbins()
				h8['ST',vst].GetAxis(H8_DIM['HEL']).SetRange(1,nbins)
				h8['SR',vst].GetAxis(H8_DIM['HEL']).SetRange(1,nbins)
				h8['ER',vst].GetAxis(H8_DIM['HEL']).SetRange(1,nbins)
				if self.DO_ETGT_BG_SUB:
					h8['ET',vst].GetAxis(H8_DIM['HEL']).SetRange(1,nbins)
			elif self.USEHEL==True:
				if que!=None:
					que.put(-1)
				sys.exit("proc(): Not yet implemented for helicity related Observables")
					
		#! 3. Prepare fout 
		if iw==0:
			fout=ROOT.TFile(os.path.join(self.OUTDIR,self.FOUTNAME),"RECREATE")
		else:
			fout=ROOT.TFile(os.path.join(self.OUTDIR,self.FOUTNAME),"UPDATE")


		#! 4. Loop of Q2-W bins 

		#! + First get Q2-W binning information from any of the h8(SEQ_DRCT,VST)
		#! + Q2-W binning information (for more details see 'get_q2w_bng()'):
		#! 	+ q2wb_crd=(q2b_crd,wb_crd)
		#! 	+ q2wb_lmt=((q2b_min,q2b_max),(wb_min,wb_max))
		q2wb_crds,q2wb_lmts=self.get_q2w_bng(h8['ER',1])
		#! Debug Q2-W binning
		# print "q2wb_crd=(q2b_crd,wb_crd) q2wb_lmt=((q2b_min,q2b_max),(wb_min,wb_max)):"
		# for q2wb_crd,q2wb_lmt in zip(q2wb_crds,q2wb_lmts):
		# 	print q2wb_crd,q2wb_lmt

		for q2wb_crd,q2wb_lmt in zip(q2wb_crds,q2wb_lmts):
			q2b_min,q2b_max=q2wb_lmt[0][0],q2wb_lmt[0][1]
			wb_min, wb_max =q2wb_lmt[1][0],q2wb_lmt[1][1]
			q2wb_name="%0.2f-%0.2f_%0.3f-%0.3f"%(q2b_min,q2b_max,wb_min,wb_max)#(q2wb_lmt[0][0],q2wb_lmt[0][1],q2wb_lmt[1][0],q2wb_lmt[1][1])
			print "proc(): Processing q2wb_name=",q2wb_name
			if q2b_min<self.Q2MIN or q2b_max>self.Q2MAX or wb_min<self.WMIN or wb_max>self.WMAX:
				print "proc(): Skipping q2wb_name=",q2wb_name
				continue
			
			#! Create h5(SEQ,VST)
			h5=OrderedDict()
			#! 3.i. h8(SEQ_DRCT,VST)=>h5(SEQ_DRCT,VST)
			for seq in (SEQ_DRCT+['ET'] if self.DO_ETGT_BG_SUB else SEQ_DRCT):
				for vst in VSTS:
					#! Make Q2,W projections
					h8[seq,vst].GetAxis(H8_DIM['Q2']).SetRange(q2wb_crd[0],q2wb_crd[0])
					h8[seq,vst].GetAxis(H8_DIM['W']).SetRange(q2wb_crd[1],q2wb_crd[1])
					h5[seq,vst]=h8[seq,vst].Projection(5,H5_PROJDIM,"E")
					h5[seq,vst].SetName('h5')
					h5[seq,vst].SetTitle('%s_%s_%d'%(q2wb_name,seq,vst))
					thntool.SetUnderOverFLowBinsToZero(h5[seq,vst])
			#! 3.ii.  h5(SEQ_DRCT,VST)=>h5(SEQ_CALC,VST)
			for vst in VSTS:
				#! SA
				h5['SA',vst]=h5['SR',vst].Clone()
				h5['SA',vst].Divide(h5['ST',vst])
				thntool.SetUnderOverFLowBinsToZero(h5['SA',vst])
				#! for SA,SR SetBinContentsAboveRelBinErrThresholdToZero()
				if (self.ACC_REL_ERR_CUT>0):
					thntool.SetBinContentsAboveRelBinErrThresholdToZero(h5['SA',vst],h5['SR',vst],self.ACC_REL_ERR_CUT)
				#! SC
				h5['SC',vst]=h5['SR',vst].Clone()
				h5['SC',vst].Divide(h5['SA',vst])
				thntool.SetUnderOverFLowBinsToZero(h5['SC',vst])
				#! SH
				h5['SH',vst]=h5['ST',vst].Clone()
				h5['SH',vst].Add(h5['SC',vst],-1)
				thntool.SetUnderOverFLowBinsToZero(h5['SH',vst])
				#! SF
				h5['SF',vst]=h5['SC',vst].Clone()
				h5['SF',vst].Add(h5['SH',vst],1)
				thntool.SetUnderOverFLowBinsToZero(h5['SF',vst])
				#! if self.DO_ETGT_BG_SUB==True: ER - (ET*R['ptgt_lse'])
				if self.DO_ETGT_BG_SUB:
					h5['ER',vst].Add(h5['ET',vst],-R['ptgt_lse'])
				#! EC
				h5['EC',vst]=h5['ER',vst].Clone()
				h5['EC',vst].Divide(h5['SA',vst])
				thntool.SetUnderOverFLowBinsToZero(h5['EC',vst]);
				#! EH
				#! 1. Set EH=SH
				h5['EH',vst]=h5['SH',vst].Clone()
				#! 2. Now scale EH by scale_factor 
				scale_factor=self.calc_EH_scale_factor(h5['EC',vst],h5['SC',vst])
				h5['EH',vst].Scale(scale_factor)
				thntool.SetUnderOverFLowBinsToZero(h5['EH',vst]);
				#! EF
				h5['EF',vst]=h5['EC',vst].Clone()
				h5['EF',vst].Add(h5['EH',vst],1)
				thntool.SetUnderOverFLowBinsToZero(h5['EF',vst]);
			#! Set the names of titles of h5(SEQ_CALC,VST) created
			for item in list(itertools.product(SEQ_CALC,VSTS)):
				seq,vst=item[0],item[1]
				h5[seq,vst].SetName('h5')
				h5[seq,vst].SetTitle('%s_%s_%d'%(q2wb_name,seq,vst))

			#! 3.iii. If self.USEHEL=false, then h5(SEQ_ALL,VST)=>h1(SEQ_ALL,VST,VARS) 
			#! + For technical reasons (see below), h1s are directly created 
			#!   under the appropriate file directory when saving the histograms (3.iv)
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


			#3.iv.  Save h5(SEQ_ALL,VST) and h1(SEQ_ALL,VST,VARS)			
			q2wb_dir=fout.GetDirectory(q2wb_name)
			if q2wb_dir==None: 
				q2wb_dir=fout.mkdir(q2wb_name)
			q2wb_dir.cd() 	
			for item in list(itertools.product(SEQ_ALL,VSTS)):
				seq,vst=item[0],item[1]
				#! Ensure seq_dir exists and cd into it
				seq_dir=q2wb_dir.GetDirectory(seq)
				if seq_dir==None:
					seq_dir=q2wb_dir.mkdir(seq)
				seq_dir.cd()
				#! Ensure vst_dir exists and cd into it
				vst_name="VST%d"%vst
				vst_dir=seq_dir.GetDirectory(vst_name)
				if vst_dir==None:
					vst_dir=seq_dir.mkdir(vst_name)
				vst_dir.cd()
				#! Save h5[seq,vst]
				h5[seq,vst].Write()
				#! Create and save h1 for every [seq,vst,var]
				if self.USEHEL==False:
					for var in VARS:
						h1=h5[seq,vst].Projection(H5_DIM[var],"E")
						h1.SetName('h1_%s'%var)
						h1.SetTitle("%s_%s_%d_%s"%(q2wb_name,seq,vst,var))
						if var=='M1' or var=='M2':
							self.setM1M2axisrange(h1,vst,var,wb_max)
			
		print "proc(): If appearing to be stuck, then either fout is still being written or Python is probably doing \"garbage collection\"(?); Wait a while!"
		if que!=None:
			que.put(0)
		fout.Write()
		fout.Close()
		print "proc(): Done processing Crs-W bin %s"%(iw+1)
		return 0
	

	def get_h8_seq_drct(self,iw):
		"""
		get h8(SEQ_DRCT,VST) for each Crs-W
		"""
		h8=OrderedDict()
		#! Fetch h8(SEQ_DRCT,vst) directly from file 
		for vst in VSTS:
			h8['ER',vst]=self.FIN['ER'].Get('d2pi/h8_%d_%d'%(iw+1,vst))
			h8['ST',vst]=self.FIN['ST'].Get('d2pi/h8_%d_%d'%(iw+1,vst))
			h8['SR',vst]=self.FIN['SR'].Get('d2pi/h8_%d_%d'%(iw+1,vst))
			if self.DO_ETGT_BG_SUB:
				h8['ET',vst]=self.FIN['ET'].Get('d2pi/h8_%d_%d'%(iw+1,vst))
		return h8

	
	#! [06-23-17] 
	#! + Replaced setM1M2axisrange() with "updated" version of function already in 
	#!   use in disp_obs.py
	#! + This "update" is with respect to the fact after inclusion of Obs_R2,
	#!   vst*{M1,M2} possibilities are extended
	#! + I expect, thus far, no impact of this change in proc_h8.py because no processing
	#!   with respect to Obs_R2 is being done here. 
	#! + Therefore, this change is mainly for integrity.
	# def setM1M2axisrange(self,h,vst,var,wmax):
	# 	if (vst==1 and var=='M1') or (vst==3 and var=='M2'):
	# 		h.GetXaxis().SetRangeUser(1.1000,wmax-0.14)
	# 	elif (vst==2 and var=='M2'):
	# 		h.GetXaxis().SetRangeUser(0.2780,wmax-0.938)

	#! [06-23-17] 
	#! + Based on updated version from disp_obs.y
	#! + Only difference is that wmax is directly passed in this version
	#!	 (In disp_obs.py, q2wbin is passed from which wmax is obtained)
	def setM1M2axisrange(self,h,vst,var,wmax):
		'''
		1,M1=p,pip 1,M2=pip,pim
		2,M1=p,pip 2,M2=pip,pim
		3,M1=p,pip 3,M2=p,pim
		'''
		#! Determine xmin,xmax=f(vst,M1/M2)
		if (var=='M1' and (vst==1 or vst==2 or vst==3)):
			xmin=MASS_P+MASS_PIP
			xmax=wmax-MASS_PIM
		elif (var=='M2' and (vst==1 or vst==2)):
			xmin=MASS_PIP+MASS_PIM
			xmax=wmax-MASS_P
		elif (var=='M2' and (vst==3)):
			xmin=MASS_P+MASS_PIM
			xmax=wmax-MASS_PIP
		#! Now set limits on h
		h.GetXaxis().SetRangeUser(xmin,xmax)
	
	def calc_EH_scale_factor(self,hN_EC,hN_SC):
		# scale_factor,nExpEvts,nSimEvts=0,0,0
		# ndims=hN_EC.GetNdimensions();
		# expBinCoord=np.zeros(ndims,'i')
		# nExpBins=hN_EC.GetNbins()
		# for iExpBin in range(nExpBins):
		# 	nExpEvts+=hN_EC.GetBinContent(iExpBin,expBinCoord)
		# 	iSimBin=hN_SC.GetBin(expBinCoord);
		# 	nSimEvts+=hN_SC.GetBinContent(iSimBin)
		# if (nSimEvts!=0):
		# 	scale_factor=nExpEvts/nSimEvts;
		# else:
		# 	print "norm for Holes=0 because nSC=0!"
		# return scale_factor

		#! [04-06-16] 
		#! Fixed this method to obtain scale_factor using SR Phase Space
		#! (Earlier I was using ER Phase Space and this was overestimating this factor)
		scale_factor,nExpEvts,nSimEvts=0,0,0
		ndims=hN_SC.GetNdimensions();
		simBinCoord=np.zeros(ndims,'i')
		nSimBins=hN_SC.GetNbins()
		for iSimBin in range(nSimBins):
			nSimEvts+=hN_SC.GetBinContent(iSimBin,simBinCoord)
			iExpBin=hN_EC.GetBin(simBinCoord);
			nExpEvts+=hN_EC.GetBinContent(iExpBin)
		if (nSimEvts!=0):
			scale_factor=nExpEvts/nSimEvts;
		else:
			print "norm for Holes=0 because nSC=0!"
		return scale_factor

	def get_q2w_bng(self,h8):
		'''
		Using an h8 (bins=H,Q2,W,M1,M2,THETA,PHI,ALPHA), return Q2-W binning information:
		+ q2wb_crds: List of 2D bin coordinates of all Q2-W bins formed by combining the h8's
		             Q2 and W bins. 
		+ q2wb_lmts: List of bin edges for each of the 2D Q-2 bin coordinates.
		'''
		#! Obtain q2wb_crds
		nq2bins=h8.GetAxis(H8_DIM['Q2']).GetNbins()
		nwbins=h8.GetAxis(H8_DIM['W']).GetNbins()

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

		return q2wb_crds,q2wb_lmts