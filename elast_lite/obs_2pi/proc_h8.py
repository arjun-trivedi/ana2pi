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
#! See http://stackoverflow.com/questions/16553506/python-ordereddict-iteration for following code
H8_DIM=OrderedDict([('HEL',0),('Q2',1),('W',2),('M1',3),('M2',4),('THETA',5),('PHI',6),('ALPHA',7)])
H5_PROJDIM=array('i',[H8_DIM['M1'],H8_DIM['M2'],H8_DIM['THETA'],H8_DIM['PHI'],
					  H8_DIM['ALPHA']])
	
H5_DIM=OrderedDict([('M1',0),('M2',1),('THETA',2),('PHI',3),('ALPHA',4)])
VARS=['M1','M2','THETA','PHI','ALPHA']

#! Enumerate measured Helicity values as per their corresponding bin number in h8's Helicity dimension (see DataAna::MakrYields())
#! + Only NEG,ZERO & POS are valid measured Helicity values
#! + UNPOL does not correspond to a measured Helicity value, but helps in writing code
UNPOL,NEG,ZERO,POS=-9999,1,2,3
HEL_NAME={UNPOL:'UNPOL',NEG:'NEG',ZERO:'ZERO',POS:'POS'}

SEQ_ALL=['ST','SR','SA','SC','SH','SF','ER','EC','EH','EF']
SEQ_DRCT=['ST','SR','ER']
SEQ_CALC=['SA','SC','SH','SF','EC','EH','EF']

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
	def __init__(self,obsdir,simnum='siml',q2min=1.25,q2max=5.25,wmin=1.300,wmax=2.125,vsts=[1,2,3],usehel=False,dbg=False):
		self.SIMNUM=simnum

		self.VSTS=vsts 

		self.USEHEL=usehel
		print "USEHEL=",self.USEHEL

		self.DBG=dbg
		print "DBG=",self.DBG

		self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX=q2min,q2max,wmin,wmax
		print "Q2MIN,Q2MAX,WMIN,WMAX=",self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX

		self.DATADIR=obsdir

		#! Prepare FIN
		self.FIN=OrderedDict()
		self.FIN['ER']=ROOT.TFile(os.path.join(self.DATADIR,'d2pi_exp','d2piR.root'))
		self.FIN['ST']=ROOT.TFile(os.path.join(self.DATADIR,'d2pi_sim',self.SIMNUM,'d2piT.root'))
		self.FIN['SR']=ROOT.TFile(os.path.join(self.DATADIR,'d2pi_sim',self.SIMNUM,'d2piR.root'))
		#! Prepare OUTDIR
		self.OUTDIR=os.path.join(self.DATADIR,self.SIMNUM)
		if not os.path.exists(self.OUTDIR): 
			sys.exit("Path %s does not exist. Exiting."%self.OUTDIR)
		#! + Prepare FOUTNAME 
		#! + Note only FOUTNAME is created here; FOUT creation is handled, with attention to 
		#!   using multiprocessing, by 'proc()'. See doc. of 'proc()' for more details
		if self.USEHEL: self.FOUTNAME="yield_hel.root"
		else:           self.FOUTNAME="yield.root"
		print "DATADIR=%s\nFIN[ER]=%s\nFIN[ST]=%s\nFIN[SR]=%s\nOUTDIR=%s\nFOUTNAME=%s"%(self.DATADIR, self.FIN['ER'].GetName(),self.FIN['ST'].GetName(),self.FIN['SR'].GetName(), self.OUTDIR, self.FOUTNAME)
		

		self.wmax=None #Used for setM1M2axisrange()

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
			3.iii. If self.USEHEL=false, then h5(SEQ_ALL,VST)=>h1(SEQ_ALL,VST) 
			3.iv.  Save h5(SEQ_ALL,VST) and h1(SEQ_ALL,VST)
				
		"""
		print "proc(): Processing Crs-W bin %s"%(iw+1)

		#! 1. Get h8(SEQ_DRCT,VST)
		h8=self.get_h8_seq_drct(iw)
		# #! Debug h8
		# for k in h8:
		# 	print k,h8[k].GetEntries()

		#! 2. Make sure H8_DIM['HEL'] is appropriately set for h8(SEQ_DRCT,vst)
		for vst in self.VSTS:
			if self.USEHEL==False:
				nbins=h8['ST',vst].GetAxis(H8_DIM['HEL']).GetNbins()
				h8['ST',vst].GetAxis(H8_DIM['HEL']).SetRange(1,nbins)
				h8['SR',vst].GetAxis(H8_DIM['HEL']).SetRange(1,nbins)
				h8['ER',vst].GetAxis(H8_DIM['HEL']).SetRange(1,nbins)
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
			if q2b_min<self.Q2MIN or q2b_max>self.Q2MAX:
				print "proc(): Skipping q2wb_name=",q2wb_name
				continue
			
			#! Create h5(SEQ,VST)
			h5=OrderedDict()
			#! 3.i. h8(SEQ_DRCT,VST)=>h5(SEQ_DRCT,VST)
			for seq in SEQ_DRCT:
				for vst in self.VSTS:
					#! Make Q2,W projections
					h8[seq,vst].GetAxis(H8_DIM['Q2']).SetRange(q2wb_crd[0],q2wb_crd[0])
					h8[seq,vst].GetAxis(H8_DIM['W']).SetRange(q2wb_crd[1],q2wb_crd[1])
					h5[seq,vst]=h8[seq,vst].Projection(5,H5_PROJDIM,"E")
					thntool.SetUnderOverFLowBinsToZero(h5[seq,vst]);
			#! 3.ii.  h5(SEQ_DRCT,VST)=>h5(SEQ_CALC,VST)
			for vst in self.VSTS:
				#! SA
				h5['SA',vst]=h5['SR',vst].Clone()
				h5['SA',vst].Divide(h5['ST',vst])
				#! SC
				h5['SC',vst]=h5['SR',vst].Clone()
				h5['SC',vst].Divide(h5['SA',vst])
				#! SH
				h5['SH',vst]=h5['ST',vst].Clone()
				h5['SH',vst].Add(h5['SC',vst],-1)
				#! SF
				h5['SF',vst]=h5['SC',vst].Clone()
				h5['SF',vst].Add(h5['SH',vst],1)
				#! EC
				h5['EC',vst]=h5['ER',vst].Clone()
				h5['EC',vst].Divide(h5['SA',vst])
				#! EH
				#! 1. Set EH=SH
				h5['EH',vst]=h5['SH',vst].Clone()
				#! 2. Now scale EH by scale_factor 
				scale_factor=self.calc_EH_scale_factor(h5['EC',vst],h5['SC',vst])
				h5['EH',vst].Scale(scale_factor)
				#! EF
				h5['EF',vst]=h5['EC',vst].Clone()
				h5['EF',vst].Add(h5['EH',vst],1)

			#! 3.iii. If self.USEHEL=false, then h5(SEQ_ALL,VST)=>h1(SEQ_ALL,VST) 

			#3.iv.  Save h5(SEQ_ALL,VST) and h1(SEQ_ALL,VST)			
			q2wb_dir=fout.GetDirectory(q2wb_name)
			if q2wb_dir==None: 
				q2wb_dir=fout.mkdir(q2wb_name)
			q2wb_dir.cd() 	

			for seq in SEQ_ALL:
				seq_dir=q2wb_dir.GetDirectory(seq)
				if seq_dir==None:
					seq_dir=q2wb_dir.mkdir(seq)
				seq_dir.cd()
				for vst in self.VSTS:
					vst_name="VST%d"%vst
					vst_dir=seq_dir.GetDirectory(vst_name)
					if vst_dir==None:
						vst_dir=seq_dir.mkdir(vst_name)
					vst_dir.cd()
					h5[seq,vst].SetName('h5')#_%s_%s'%(seq,vst_name))
					h5[seq,vst].SetTitle('%s_%s_%s'%(q2wb_name,seq,vst))
					h5[seq,vst].Write()
			
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
		for vst in self.VSTS:
			h8['ER',vst]=self.FIN['ER'].Get('d2pi/h8_%d_%d'%(iw+1,vst))
			h8['ST',vst]=self.FIN['ST'].Get('d2pi/h8_%d_%d'%(iw+1,vst))
			h8['SR',vst]=self.FIN['SR'].Get('d2pi/h8_%d_%d'%(iw+1,vst))
		return h8

	
	def setM1M2axisrange(self,h,vst,var):
		if (vst=="VST1" and var=='M1') or (vst=="VST3" and var=='M2'):
			h.GetXaxis().SetRangeUser(1.1000,self.wmax-0.14)
		if (vst=="VST2" and var=='M2'):
			h.GetXaxis().SetRangeUser(0.2780,self.wmax-0.938)
	
	def calc_EH_scale_factor(self,hN_EC,hN_SC):
		scale_factor,nExpEvts,nSimEvts=0,0,0
		ndims=hN_EC.GetNdimensions();
		expBinCoord=np.zeros(ndims,'i')
		nExpBins=hN_EC.GetNbins()
		for iExpBin in range(nExpBins):
			nExpEvts+=hN_EC.GetBinContent(iExpBin,expBinCoord)
			iSimBin=hN_SC.GetBin(expBinCoord);
			nSimEvts+=hN_SC.GetBinContent(iSimBin)
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