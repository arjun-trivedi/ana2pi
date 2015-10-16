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

class ProcH8:
	"""
	+ Accomplishes d2pi.root(ER,ST,SR) -> yield.root  
		+ if self.USEHEL=false: d2pi.root -> yield.root
		+ if self.USEHEL=true:  d2pi.root -> yield_hel.root 

	+ The direct interface for the user is the 'execute()' method that calls 'proc()' for h8(VST,SEQ) in each Crs-W bin. 
		+ The 'execute()' method uses Python's multiprocessing module to accomplish the task of processing each h8(VST,SEQ) in 
	a separate thread, thereby not having to load all h8(VST,SEQ) from every Crs-W bin at once and exhausting the system's memory (gothe14=12GB RAM)
		+ For more technical details, see documentation of 'execute()'

	+ For each q2w-bin (Please see details of proc(),proc_h5() and proc_h1() for details):
		+ if self.USEHEL=false: h8(VST,SEQ) => h5(VST,SEQ) => h1(VST,SEQ)
		+ if self.USEHEL=true:  h8(VST,SEQ) => h5(VST,SEQ) 
	
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
		1. Get h8(ER/ST/SR,VST) per Crs-W from FIN[SEQ] 
		2. Calculate h8(SA/SC/SH/SF/EC/EH/EF,VST) 
		3. Make sure H8_DIM['HEL'] is appropriately set
		4. Open fout (self.FOUTNAME) in "appropriate" mode:
			+ if 1st Crs-W bin, then "RECREATE", else "UPDATE"
			
				
		"""
		print "proc(): Processing Crs-W bin %s"%(iw+1)

		#! 1. Get h8(ER/ST/SR,VST)
		h8=self.get_h8(iw)
		# #! Debug h8
		# for k in h8:
		# 	print k,h8[k].GetEntries()

		#! 2. Calculate h8(SA/SC/SH/SF/EC/EH/EF,VST) 
		for vst in self.VSTS:
			#! SA
			h8['SA',vst]=h8['SR',vst].Clone()
			h8['SA',vst].Divide(h8['ST',vst])
			#! SC
			h8['SC',vst]=h8['SR',vst].Clone()
			h8['SC',vst].Divide(h8['SA',vst])
			#! SH
			h8['SH',vst]=h8['ST',vst].Clone()
			h8['SH',vst].Add(h8['SC',vst],-1)
			#! SF
			h8['SF',vst]=h8['SC',vst].Clone()
			h8['SF',vst].Add(h8['SH',vst],1)
			h8['SF',vst].SetName("h8_SF_VST%d"%(vst))
			#! EC
			h8['EC',vst]=h8['ER',vst].Clone()
			h8['EC',vst].Divide(h8['SA',vst])
			#! EH
			#! 1. Set EH=SH
			h8['EH',vst]=h8['SH',vst].Clone()
			#! 2. Now scale EH by scale_factor 
			scale_factor=self.calc_EH_scale_factor(h8['EC',vst],h8['SC',vst])
			h8['EH',vst].Scale(scale_factor)
			#! EF
			h8['EF',vst]=h8['EC',vst].Clone()
			h8['EF',vst].Add(h8['EH',vst],1)
			
		#! Debug h8
		for k in h8:
			print k,h8[k].GetEntries()

		#! 3. Make sure H8_DIM['HEL'] is appropriately set
		for seq in SEQ_ALL:
				for vst in self.VSTS:
					if self.USEHEL==False:
						h8[seq,vst].GetAxis(H8_DIM['HEL']).SetRange()
					elif self.USEHEL==True:
						if que!=None:
							que.put(-1)
						sys.exit("proc(): Not yet implemented for helicity related Observables")
			
		#! 4. Prepare fout
		if iw==0:
			fout=ROOT.TFile(os.path.join(self.OUTDIR,self.FOUTNAME),"RECREATE")
		else:
			fout=ROOT.TFile(os.path.join(self.OUTDIR,self.FOUTNAME),"UPDATE")

		#! 5. For every Q2-W bin in h8, make projections and write them to fout:
		#!   5.i.  h8(SEQ,VST)->h5(SEQ,VST)
		#!   5.ii. h8(SEQ,VST)->h1(SEQ,VST) only if self.USEHEL=false

		#! First get Q2-W binning information from any of the h8(SEQ,VST)
		q2wb_crds,q2wb_lmts=self.get_q2w_bng(h8['ER',1])
		#! Debug Q2-W binning
		# print "q2wb_crd=(q2b_crd,wb_crd) q2wb_lmt=((q2b_min,q2b_max),(wb_min,wb_max)):"
		# for q2wb_crd,q2wb_lmt in zip(q2wb_crds,q2wb_lmts):
		# 	print q2wb_crd,q2wb_lmt

		
		#! Now loop of Q2-W bins and make projections
		for q2wb_crd,q2wb_lmt in zip(q2wb_crds,q2wb_lmts):
			q2b_min,q2b_max=q2wb_lmt[0][0],q2wb_lmt[0][1]
			wb_min,wb_max  =q2wb_lmt[1][0],q2wb_lmt[1][1]
			q2wb_name="%0.2f-%0.2f_%0.3f-%0.3f"%(q2b_min,q2b_max,wb_min,wb_max)#(q2wb_lmt[0][0],q2wb_lmt[0][1],q2wb_lmt[1][0],q2wb_lmt[1][1])
			print "proc(): Processing q2wb_name=",q2wb_name
			if q2b_min<self.Q2MIN or q2b_max>self.Q2MAX:
				print "proc(): Skipping q2wb_name=",q2wb_name
				continue
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
					
					print "Processing seq,vst=%s,%s"%(seq,vst_name)
					#! 5.i.  h8(SEQ,VST)->h5(SEQ,VST)
					#! Make Q2,W projections
					h8[seq,vst].GetAxis(H8_DIM['Q2']).SetRange(q2wb_crd[0],q2wb_crd[0])
					h8[seq,vst].GetAxis(H8_DIM['W']).SetRange(q2wb_crd[1],q2wb_crd[1])
					h5=h8[seq,vst].Projection(5,H5_PROJDIM,"E")
					thntool.SetUnderOverFLowBinsToZero(h5);
					h5.SetName('h5')#_%s_%s'%(seq,vst_name))
					h5.SetTitle('%s_%s_%s'%(q2wb_name,seq,vst))
					h5.Write()

		print "proc(): If appearing to be stuck, then either fout is still being written or Python is probably doing \"garbage collection\"(?); Wait a while!"
		if que!=None:
			que.put(0)
		fout.Write()
		fout.Close()
		print "proc(): Done processing Crs-W bin %s"%(iw+1)
		return 0
									

		


		# #! 3. Loop over Q2-W bins
		# #! Get nq2bins & nwbins from h8(vst=1,seq=R)
		# NQ2BINS=h8(1,'R').GetAxis(H8_DIM['Q2']).GetNbins()
		# NWBINS= h8(1,'R').GetAxis(H8_DIM['W']).GetNbins()
		# #! Prepare list of seq to process
		# if self.DTYP=='exp':seql=['R']
		# if self.DTYP=='sim':seql=['T','R']
		# # Now start looping
		# for iq2bin in range(nq2bins):
		# 	for iwbin in range(nwbins):
		# 		for vst in SELF.VST:
		# 			for seq in seql:
		# 				#! Following is used to set up fout:Q2Wdir
		# 				q2bin_min=h8(vst,seq).GetAxis(H8_DIM['Q2']).GetBinLowEdge(iq2bin+1)
		# 				q2bin_max=h8(vst,seq).GetAxis(H8_DIM['Q2']).GetBinLowEdge(iq2bin+2)
		# 				wbin_min=h8(vst,seq).GetAxis(H8_DIM['W']).GetBinLowEdge(iwbin+1)
		# 				wbin_max=h8(vst,seq).GetAxis(H8_DIM['W']).GetBinLowEdge(iwbin+2)
		# 				q2wbin="%0.2f-%0.2f_%0.3f-%0.3f"%(q2bin_min,q2bin_max,)

		
		# #! 3. Set up Q2,W bng
		# self.Q2BNG,self.WBNG=None,None
		# try:
		# 	self.Q2BNG,self.WBNG=atlib.init_q2wbng2(q2w_bng.Q2W_BNG[iw])
		# except IndexError:
		# 	print "IndexError while setting q2w bng. Exiting"
		# except:
		# 	print "The following exception occured. Exiting."
		# 	raise

		# #! Test Q2WBNG, WBNG
		# print "Q2 bins=",self.Q2BNG['BINS']
		# print "W bins=",self.WBNG['BINS']
		
		# #! 4. Loop over Q2W-bins
		# for i in range(self.Q2BNG['NBINS']):
		# 	if self.Q2BNG['BINS_LE'][i]<self.Q2MIN or self.Q2BNG['BINS_UE'][i]>self.Q2MAX:continue
		# 	if self.DBG:
		# 		if i>1: break
		# 	for j in range(self.WBNG['NBINS']):
		# 		if self.WBNG['BINS_LE'][j]<self.WMIN or self.WBNG['BINS_UE'][j]>self.WMAX:continue
		# 		if self.DBG:
		# 			if j>1: break
		# 		q2wbin="%0.2f-%0.2f_%0.3f-%0.3f"%(self.Q2BNG['BINS_LE'][i],self.Q2BNG['BINS_UE'][i],self.WBNG['BINS_LE'][j],self.WBNG['BINS_UE'][j])
		# 		q2wbindir=fout.mkdir(q2wbin)
		# 		q2wbintitle="[%0.2f,%0.2f)_[%0.3f,%0.3f)"%(self.Q2BNG['BINS_LE'][i],self.Q2BNG['BINS_UE'][i],self.WBNG['BINS_LE'][j],self.WBNG['BINS_UE'][j])
		# 		self.wmax=self.WBNG['BINS_UE'][j]
		# 		print "*** Processing Q2,W %s ***"%q2wbintitle

		# 		#! 4.i   Set appropriate Q2 & W projection range for h8(VST,SEQ);SEQ=ST/SR/ER
		# 		print "*** h8(VST,SEQ) range setting for Q2 & W dimensions ... ***"
		# 		for vst in self.VSTS:
		# 			vstdir=q2wbindir.mkdir(vst)
									
		# 			if self.DTYP=='exp':seq_h8=['R']
		# 			if self.DTYP=='sim':seq_h8=['T','R']
		# 			for seq in seq_h8:
		# 				#!-- Q2
		# 				q2bin_le=h8[vst,seq].GetAxis(H8_DIM['Q2']).FindBin(self.Q2BNG['BINS_LE'][i]+self.Q2BNG['BINW']/2)
		# 				q2bin_ue=h8[vst,seq].GetAxis(H8_DIM['Q2']).FindBin(self.Q2BNG['BINS_UE'][i]-self.Q2BNG['BINW']/2)
		# 				h8[vst,seq].GetAxis(H8_DIM['Q2']).SetRange(q2bin_le,q2bin_ue)
		# 				#!-- W
		# 				wbin_le=h8[vst,seq].GetAxis(H8_DIM['W']).FindBin(self.WBNG['BINS_LE'][j]+self.WBNG['BINW']/2)
		# 				wbin_ue=h8[vst,seq].GetAxis(H8_DIM['W']).FindBin(self.WBNG['BINS_UE'][j]-self.WBNG['BINW']/2)
		# 				h8[vst,seq].GetAxis(H8_DIM['W']).SetRange(wbin_le,wbin_ue)
		# 				print "h8(%s,%s):finished setting range for Q2-,W-bin = %s ***"%(vst,seq,q2wbintitle)
		# 				#! Project out hq2w & save in FOUT 
		# 				if vstdir.GetDirectory(seq)==None:vstdir.mkdir(seq).cd()
		# 				else:vstdir.cd(seq)
		# 				hq2w=h8[vst,seq].Projection(H8_DIM['Q2'],H8_DIM['W'],"E")
		# 				hq2w.SetName('h_q2Vw')
		# 				hq2w.SetTitle("%s_%s_%s_q2w"%(q2wbintitle,vst,seq))
		# 				#hq2w.Write() #! Am using Write() only for THnSparse Objects; for regular TH1s, Write() seems to save all cycles of an object(name;cycle)
		# 		print "*** Done h8 Q2,W range setting"

		# 		#! 4.ii. h8(VST,SEQ) => h5(VST,SEQ) => (if self.USEHEL=false) h1(VST,SEQ) 
		# 		if self.USEHEL:
		# 			print "*** Processing h8(VST,SEQ) => h5_UNPOL(VST,SEQ),h5_POS(VST,SEQ),h5_NEG(VST,SEQ)  ... ***"
		# 			self.proc_h5(h8,q2wbin,q2wbindir,q2wbintitle,hel=UNPOL)
		# 			self.proc_h5(h8,q2wbin,q2wbindir,q2wbintitle,hel=POS)
		# 			self.proc_h5(h8,q2wbin,q2wbindir,q2wbintitle,hel=NEG)
		# 			print "*** Done processing h8(VST,SEQ) => h5_UNPOL(VST,SEQ),h5_POS(VST,SEQ),h5_NEG(VST,SEQ)  ... ***"
		# 		else:
		# 			print "*** Processing h8(VST,SEQ) => h5_UNPOL(VST,SEQ) ==> h1(VST,SEQ) ... ***"
		# 			h5=self.proc_h5(h8,q2wbin,q2wbindir,q2wbintitle,hel=UNPOL)#,vst,vstdir)
		# 			self.proc_h1(h5,q2wbin,q2wbindir,q2wbintitle)#,vst,vstdir)
		# 			print "*** Done processing h8(VST,SEQ) => h5(VST,SEQ) ==> (if self.USEHEL=false) h1(VST,SEQ) ... ***"
		
	def get_h8(self,iw):
		"""
		get h8(VST,SEQ) for each Crs-W
		"""
		h8=OrderedDict()
		#! Fetch h8-SR/ST/SR directly from file 
		for vst in self.VSTS:
			h8['ER',vst]=self.FIN['ER'].Get('d2pi/h8_%d_%d'%(iw+1,vst))
			h8['ST',vst]=self.FIN['ST'].Get('d2pi/h8_%d_%d'%(iw+1,vst))
			h8['SR',vst]=self.FIN['SR'].Get('d2pi/h8_%d_%d'%(iw+1,vst))
			#! Fetch h8-T
			# if self.DTYP=='sim':
			# 	print "Going to get",'d2pi/h8_%d_%d'%(iw+1,vst),"..."
			# 	h8[vst,'T']=self.FIN_T.Get('d2pi/h8_%d_%d'%(iw+1,vst))
			# 	print "Done getting",'d2pi/h8_%d_%d'%(iw+1,vst)
			# #! Fetch h8-R 
			# print "Going to get",'d2pi/h8_%d_%d'%(iw+1,vst),"..."
			# h8[vst,'R']=self.FIN_R.Get('d2pi/h8_%d_%d'%(iw+1,vst))
			# print "Done getting",'d2pi/h8_%d_%d'%(iw+1,vst)
		return h8

	def proc_h5(self,h8,q2wbin,q2wbindir,q2wbintitle,hel):#,vst,vstdir):
		"""
		+ This method is called by proc(), per Q2W-bin, after setting the appropriate Q2,W range on h8(VST,SEQ). 
		+ In this method, the following is accomplished
			+ h8(VST,SEQ) => h5(VST,SEQ)

		+ @hel=UNPOL,POS,NEG,ZERO
		"""
		h5=OrderedDict()
		print "*** Processing h8(VST,SEQ)->h5(VST,SEQ) for %s ***"%q2wbintitle

		#! seq to directly project from h8
		if self.DTYP=='exp':seq_h5=['R']
		if self.DTYP=='sim':seq_h5=['T','R']

		for vst in self.VSTS:
			vstdir=q2wbindir.GetDirectory(vst)

			#! 1. Project out h5-T/R
			for seq in seq_h5:
				#! First set HEL
				if hel==UNPOL:
					h8[vst,seq].GetAxis(H8_DIM['HEL']).SetRange()
				elif hel==POS or hel==NEG:
					h8[vst,seq].GetAxis(H8_DIM['HEL']).SetRange(hel,hel)

				if vstdir.GetDirectory(seq)==None:vstdir.mkdir(seq).cd()
				else:vstdir.cd(seq)
				h5[vst,seq]=h8[vst,seq].Projection(5,H5_PROJDIM,"E")
				thntool.SetUnderOverFLowBinsToZero(h5[vst,seq]);
				h5[vst,seq].SetName('h5_%s'%HEL_NAME[hel])
				h5[vst,seq].SetTitle('%s_%s_%s_%s'%(q2wbintitle,vst,seq,HEL_NAME[hel]))
				h5[vst,seq].Write()

			#! 2. Calculate h5-A
			if vstdir.GetDirectory('A')==None:vstdir.mkdir('A').cd()
			else:vstdir.cd('A')
			if self.DTYP=='sim':
				h5[vst,'A']=h5[vst,'R'].Clone()
				h5[vst,'A'].Divide(h5[vst,'T'])
				thntool.SetUnderOverFLowBinsToZero(h5[vst,'A']);
				h5[vst,'A'].SetName('h5_%s'%HEL_NAME[hel])
				h5[vst,'A'].SetTitle('%s_%s_%s_%s'%(q2wbintitle,vst,'A',HEL_NAME[hel]))
			if self.DTYP=='exp':
				h5[vst,'A']=self.FIN_SIMYIELD.Get('%s/%s/A/h5_UNPOL'%(q2wbin,vst))
			h5[vst,'A'].Write()

			#! 3. Calculate h5-C
			if vstdir.GetDirectory('C')==None:vstdir.mkdir('C').cd()
			else:vstdir.cd('C')
			h5[vst,'C']=h5[vst,'R'].Clone()
			h5[vst,'C'].Divide(h5[vst,'A'])
			thntool.SetUnderOverFLowBinsToZero(h5[vst,'C']);
			h5[vst,'C'].SetName('h5_%s'%HEL_NAME[hel])
			h5[vst,'C'].SetTitle('%s_%s_%s_%s'%(q2wbintitle,vst,'C',HEL_NAME[hel]))
			h5[vst,'C'].Write()

			#! 4. Calculate h5-H
			if vstdir.GetDirectory('H')==None:vstdir.mkdir('H').cd()
			else:vstdir.cd('H')
			if self.DTYP=='sim':
				h5[vst,'H']=h5[vst,'T'].Clone()
				h5[vst,'H'].Add(h5[vst,'C'],-1)
				thntool.SetUnderOverFLowBinsToZero(h5[vst,'H'])
				h5[vst,'H'].SetName('h5_%s'%HEL_NAME[hel])
				h5[vst,'H'].SetTitle('%s_%s_%s_%s'%(q2wbintitle,vst,'H',HEL_NAME[hel]))
			if self.DTYP=='exp':
				#!first copy sim-HOLE into exp-HOLE
				h5[vst,'H']=self.FIN_SIMYIELD.Get('%s/%s/H/h5_UNPOL'%(q2wbin,vst))
				#!now get SC for obtaining normalization factor
				h5_SIM_C=self.FIN_SIMYIELD.Get('%s/%s/C/h5_UNPOL'%(q2wbin,vst))
				#!calculate normalization factor
				norm=self.calcNorm(h5[vst,'C'],h5_SIM_C);
				#!normalize exp-HOLE
				h5[vst,'H'].Scale(norm);
				thntool.SetUnderOverFLowBinsToZero(h5[vst,'H']);
			h5[vst,'H'].Write()

			#! Calculate h5-F
			if vstdir.GetDirectory('F')==None:vstdir.mkdir('F').cd()
			else:vstdir.cd('F')
			h5[vst,'F']=h5[vst,'C'].Clone()
			h5[vst,'F'].Add(h5[vst,'H'],1)
			h5[vst,'F'].SetName('h5_%s'%HEL_NAME[hel])
			h5[vst,'F'].SetTitle('%s_%s_%s_%s'%(q2wbintitle,vst,'F',HEL_NAME[hel]))
			h5[vst,'F'].Write()
		print "*** done processing h8(VST,SEQ)->h5(VST,SEQ) for %s ***"%q2wbintitle
		return h5

	def proc_h1(self,h5,q2wbin,q2wbindir,q2wbintitle):#,vst,vstdir):
		"""
		+ This method is called by proc() after setting the appropriate Q2,W range on h8(VST,SEQ) for a q2w-bin and calling
		  proc_h5(). 
		+ In this method, the following is accomplished
			+ h5(VST,SEQ) => h1(VST,SEQ)
			+ Note that that an OrderedDict for h1 is not created, for the objects are directly saved and there 
			  is no need to keep track of them by keys
		"""
		print "*** Processing h5(VST,SEQ)->h1(VST,SEQ) for %s ***"%q2wbintitle

		#! Define SEQ for which h5s are directly projected to h1
		if self.DTYP=='exp':
			seq_h1=['R','C','A','H','F']
		if self.DTYP=='sim':
			seq_h1=['T','R','C','A','H','F']

		for vst in self.VSTS:
			vstdir=q2wbindir.GetDirectory(vst)
			for seq in seq_h1:
				vstdir.cd(seq)
				#! Project out h1s
				for var in VARS:
					h1=h5[vst,seq].Projection(H5_DIM[var],"E")
					if var=='M1' or var=='M2':
						self.setM1M2axisrange(h1,vst,var)
					h1.SetName('h_%s'%var)
					h1.SetTitle("%s_%s_%s_%s"%(q2wbintitle,vst,seq,var))
					#h1.Write() #! Am using Write() only for THnSparse Objects; for regular TH1s, Write() seems to save all cycles of an object(name;cycle)
		print "*** Processing h5(VST,SEQ)->h1(VST,SEQ) for %s ***"%q2wbintitle
		return

	def setM1M2axisrange(self,h,vst,var):
		if (vst=="VST1" and var=='M1') or (vst=="VST3" and var=='M2'):
			h.GetXaxis().SetRangeUser(1.1000,self.wmax-0.14)
		if (vst=="VST2" and var=='M2'):
			h.GetXaxis().SetRangeUser(0.2780,self.wmax-0.938)
	
	# def calcNorm(self,h5D_EC,h5D_SC):
	# 	norm,nExpEvts,nSimEvts=0,0,0
	# 	expBinCoord=np.zeros(5,'i')
	# 	nExpBins = h5D_EC.GetNbins()
	# 	for iExpBin in range(nExpBins):
	# 		nExpEvts+=h5D_EC.GetBinContent(iExpBin,expBinCoord)
	# 		iSimBin=h5D_SC.GetBin(expBinCoord);
	# 		nSimEvts+=h5D_SC.GetBinContent(iSimBin)
	# 	if (nSimEvts!=0):
	# 		norm=nExpEvts/nSimEvts;
	# 	else:
	# 		print "norm for Holes=0 because nSC=0!"
	# 	return norm

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
		nq2bins=h8.GetAxis(H8_DIM['Q2']).GetNbins()
		nwbins=h8.GetAxis(H8_DIM['W']).GetNbins()

		q2b_crds=range(1,nq2bins+1)
		wb_crds=range(1,nwbins+1)
		q2wb_crds=list(itertools.product(q2b_crds,wb_crds))

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