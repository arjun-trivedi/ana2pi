#! Imports
from __future__ import division
from math import *

from array import *
from collections import OrderedDict

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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
H7_DIM=OrderedDict([('Q2',0),('W',1),('M1',2),('M2',3),('THETA',4),('PHI',5),('ALPHA',6)])
H5_PROJDIM=array('i',[H7_DIM['M1'],H7_DIM['M2'],H7_DIM['THETA'],H7_DIM['PHI'],
					  H7_DIM['ALPHA']])
	
H5_DIM=OrderedDict([('M1',0),('M2',1),('THETA',2),('PHI',3),('ALPHA',4)])
VARS=['M1','M2','THETA','PHI','ALPHA']

#! Enumerate measured Helicity values as per their corresponding bin number in h7's Helicity dimension (see DataAna::MakrYields())
#! + Only NEG,ZERO & POS are valid measured Helicity values
#! + UNPOL does not correspond to a measured Helicity value, but helps in writing code
#UNPOL,NEG,ZERO,POS=-9999,1,2,3
#HEL_NAME={UNPOL:'UNPOL',NEG:'NEG',ZERO:'ZERO',POS:'POS'}

class ProcYieldsLite:
	"""
	+ Accomplishes d2pi.root -> yield.root
		+ if self.USEHEL=false: d2pi.root -> yield_exp/sim.root
		+ if self.USEHEL=true:  d2pi.root -> yield_exp_hel.root 

	+ The direct interface for the user is the 'execute()' method that calls 'proc()' for h7(VST,SEQ) in each Crs-W bin. 
		+ The 'execute()' method uses Python's multiprocessing module to accomplish the task of processing each h7(VST,SEQ) in 
	a separate thread, thereby not having to load all h7(VST,SEQ) from every Crw-W bin at once and exhausting the system's memory (atrivedi-
	laptop; 12GB RAM)

	+ For each q2w-bin (Please see details of proc(),proc_h5() and proc_h1() for details):
		+ if self.USEHEL=false: h7(VST,SEQ) => h5(VST,SEQ) => h1(VST,SEQ)
		+ if self.USEHEL=true:  h7(VST,SEQ) => h5(VST,SEQ) 
	
	"""
	def __init__(self,obsdir,dtyp,simnum='siml',q2min=1.25,q2max=5.25,wmin=1.300,wmax=2.125,vsts=[1,2,3],usehel=False,dbg=False):
		self.DTYP=""
		if dtyp=='sim':self.DTYP='sim'
		if dtyp=='exp':self.DTYP='exp'
		if not(self.DTYP=='exp' or self.DTYP=='sim'):
			sys.exit("dtyp is neither EXP or SIM! Exiting")
		print "dtyp=%s"%dtyp

		self.SIMNUM=simnum

		self.VSTS=vsts 

		# self.USEHEL=False
		# if self.DTYP=='exp':
		# 	self.USEHEL=usehel
		# print "USEHEL=",self.USEHEL

		self.DBG=dbg
		print "DBG=",self.DBG

		self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX=q2min,q2max,wmin,wmax
		print "Q2MIN,Q2MAX,WMIN,WMAX=",self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX

		self.DATADIR=obsdir

		#! Note that self.EXPT is not really used yet
		self.EXPT=""
		if "e1f" in self.DATADIR:#.find("e1f"):
			self.EXPT="e1f"
		elif "e16" in self.DATADIR:#self.DATADIR.find("e16"):
			self.EXPT="e16"
		if not(self.EXPT=='e1f' or self.EXPT=='e16'):
               		print "expt is neither e1f or e16"
		else:
			print "EXPT=",self.EXPT

		if self.DTYP=='exp':
			self.FIN_R=ROOT.TFile(os.path.join(self.DATADIR,'d2pi_exp','d2piR.root'))
			self.OUTDIR=os.path.join(self.DATADIR,self.SIMNUM)
			if not os.path.exists(self.OUTDIR): #! This path should already exist when making yield_exp
				sys.exit("Path %s does not exist. Exiting."%self.OUTDIR)
			self.FIN_SIMYIELD=ROOT.TFile(os.path.join( self.OUTDIR,"yield_sim.root"))
			self.FOUTNAME="yield_exp.root"
			# if self.USEHEL: self.FOUTNAME="yield_exp_hel.root"
			# else:           self.FOUTNAME="yield_exp.root"
			print "DATADIR=%s\nOUTDIR=%s\nFIN_R=%s\nFIN_SIMYIELD=%s\nFOUTNAME=%s/%s"%(self.DATADIR,self.OUTDIR,self.FIN_R.GetName(),self.FIN_SIMYIELD.GetName(),self.OUTDIR,self.FOUTNAME)
		if self.DTYP=='sim':
			self.FIN_T=ROOT.TFile(os.path.join(self.DATADIR,'d2pi_sim',self.SIMNUM,'d2piT.root'))
			self.FIN_R=ROOT.TFile(os.path.join(self.DATADIR,'d2pi_sim',self.SIMNUM,'d2piR.root'))
			self.OUTDIR=os.path.join(self.DATADIR,self.SIMNUM)
			if not os.path.exists(self.OUTDIR):
				os.makedirs(self.OUTDIR)
			self.FOUTNAME="yield_sim.root"
			print "DATADIR=%s\nOUTDIR=%s\nFIN_T=%s\nFIN_R=%s\nFOUT=%s/%s"%(self.DATADIR,self.OUTDIR,self.FIN_T.GetName(),self.FIN_R.GetName(),self.OUTDIR,self.FOUTNAME)

		self.wmax=None #Used for setM1M2axisrange()

	def execute(self):
		#! Loop over Crw-W bins

		queue=Queue()
		for iw in range(q2w_bng.NBINS_WCRS):
			if self.DBG:
				if iw>1: break
			p=Process(target=self.proc, args=(iw,queue))
			p.start()
			p.join() # this blocks until the process terminates
			res=queue.get()
			print "Result from proc(%d)=%d"%(iw+1,res)
		# self.FOUT.Write()
		# self.FOUT.Close()

		# #queue=Queue()
		# #for iw in range(q2w_bng.NBINS_WCRS):
		# for iw in range(2):
		# 	p=Process(target=self.proc, args=(iw,self.FOUT,None))
		# 	p.start()
		# 	#p.join() # this blocks until the process terminates
		# 	# res=queue.get()
		# 	# print "Result from proc(%d)=%d"%(iw+1,res)
		# # self.FOUT.Write()
		# self.FOUT.Close()

		# #queue=Queue()
		# #for iw in range(q2w_bng.NBINS_WCRS):
		# for iw in range(2):
		# 	p=Process(target=self.proc, args=(iw,None))
		# 	p.start()
		# 	p.join() # this blocks until the process terminates
		# 	# res=queue.get()
		# 	# print "Result from proc(%d)=%d"%(iw+1,res)
		# # self.FOUT.Write()
		# # self.FOUT.Close()

		# for iw in range(2):
		# 	self.proc(iw,self.FOUT)
		# #self.FOUT.Write()
		# self.FOUT.Close()

		# for iw in range(2):
		# 	self.proc(iw)
		# #self.FOUT.Write()
		# #self.FOUT.Close()

		print "execute():Done"
	
	def proc(self,iw,que=None):
		"""
		1. Get h7(VST,SEQ) per Crs-W from FIN_T,FIN_R (=d2piT.root,d2piR.root)
		2. Open fout (self.FOUTNAME) in "appropriate" mode (if 1st Crs-W bin, then "RECREATE", else "UPDATE")
		3. Setup Q2W binning as per Crs-W bin
		4. Loop over Q2W-bins
			4.i   Set appropriate Q2 & W projection range for h7(VST,SEQ)
			4.ii. h7(VST,SEQ) => h5(VST,SEQ) => (if self.USEHEL=false) h1(VST,SEQ) 
				
		"""
		print "*** Processing Crs-W bin %s ***"%(iw+1)

		#! 1. Get h7s from FIN_T,FIN_R (d2piT.root,d2piR.root)
		h7=self.get_h7(iw)
		#! Debug h7
		for k in h7:
			print k,h7[k].GetEntries()

		#! 2. Prepare fout
		if iw==0:
			fout=ROOT.TFile(os.path.join(self.OUTDIR,self.FOUTNAME),"RECREATE")
		else:
			fout=ROOT.TFile(os.path.join(self.OUTDIR,self.FOUTNAME),"UPDATE")
		
		#! 3. Set up Q2,W bng
		self.Q2BNG,self.WBNG=None,None
		try:
			self.Q2BNG,self.WBNG=atlib.init_q2wbng2(q2w_bng.Q2W_BNG[iw])
		except IndexError:
			print "IndexError while setting q2w bng. Exiting"
		except:
			print "The following exception occured. Exiting."
			raise

		#! Test Q2WBNG, WBNG
		print "Q2 bins=",self.Q2BNG['BINS']
		print "W bins=",self.WBNG['BINS']
		
		#! 4. Loop over Q2W-bins
		for i in range(self.Q2BNG['NBINS']):
			if self.Q2BNG['BINS_LE'][i]<self.Q2MIN or self.Q2BNG['BINS_UE'][i]>self.Q2MAX:continue
			if self.DBG:
				if i>1: break
			for j in range(self.WBNG['NBINS']):
				if self.WBNG['BINS_LE'][j]<self.WMIN or self.WBNG['BINS_UE'][j]>self.WMAX:continue
				if self.DBG:
					if j>1: break
				q2wbin="%0.2f-%0.2f_%0.3f-%0.3f"%(self.Q2BNG['BINS_LE'][i],self.Q2BNG['BINS_UE'][i],self.WBNG['BINS_LE'][j],self.WBNG['BINS_UE'][j])
				q2wbindir=fout.mkdir(q2wbin)
				q2wbintitle="[%0.2f,%0.2f)_[%0.3f,%0.3f)"%(self.Q2BNG['BINS_LE'][i],self.Q2BNG['BINS_UE'][i],self.WBNG['BINS_LE'][j],self.WBNG['BINS_UE'][j])
				self.wmax=self.WBNG['BINS_UE'][j]
				print "*** Processing Q2,W %s ***"%q2wbintitle

				#! 4.i   Set appropriate Q2 & W projection range for h7(VST,SEQ);SEQ=ST/SR/ER
				print "*** h7(VST,SEQ) range setting for Q2 & W dimensions ... ***"
				for vst in self.VSTS:
					vst_name='VST%d'%vst
					vstdir=q2wbindir.mkdir(vst_name)
									
					if self.DTYP=='exp':seq_h7=['R']
					if self.DTYP=='sim':seq_h7=['T','R']
					for seq in seq_h7:
						#!-- Q2
						q2bin_le=h7[vst_name,seq].GetAxis(H7_DIM['Q2']).FindBin(self.Q2BNG['BINS_LE'][i]+self.Q2BNG['BINW']/2)
						q2bin_ue=h7[vst_name,seq].GetAxis(H7_DIM['Q2']).FindBin(self.Q2BNG['BINS_UE'][i]-self.Q2BNG['BINW']/2)
						h7[vst_name,seq].GetAxis(H7_DIM['Q2']).SetRange(q2bin_le,q2bin_ue)
						#!-- W
						wbin_le=h7[vst_name,seq].GetAxis(H7_DIM['W']).FindBin(self.WBNG['BINS_LE'][j]+self.WBNG['BINW']/2)
						wbin_ue=h7[vst_name,seq].GetAxis(H7_DIM['W']).FindBin(self.WBNG['BINS_UE'][j]-self.WBNG['BINW']/2)
						h7[vst_name,seq].GetAxis(H7_DIM['W']).SetRange(wbin_le,wbin_ue)
						print "h7(%s,%s):finished setting range for Q2-,W-bin = %s ***"%(vst_name,seq,q2wbintitle)
						#! Project out hq2w & save in FOUT 
						if vstdir.GetDirectory(seq)==None:vstdir.mkdir(seq).cd()
						else:vstdir.cd(seq)
						hq2w=h7[vst_name,seq].Projection(H7_DIM['Q2'],H7_DIM['W'],"E")
						hq2w.SetName('h_q2Vw')
						hq2w.SetTitle("%s_%s_%s_q2w"%(q2wbintitle,vst_name,seq))
						#hq2w.Write() #! Am using Write() only for THnSparse Objects; for regular TH1s, Write() seems to save all cycles of an object(name;cycle)
				print "*** Done h7 Q2,W range setting"

				#! 4.ii. h7(VST,SEQ) => h5(VST,SEQ) => (if self.USEHEL=false) h1(VST,SEQ) 
				print "*** Processing h7(VST,SEQ) => h5(VST,SEQ) ==> h1(VST,SEQ) ... ***"
				h5=self.proc_h5(h7,q2wbin,q2wbindir,q2wbintitle)#,vst_name,vstdir)
				self.proc_h1(h5,q2wbin,q2wbindir,q2wbintitle)#,vst_name,vstdir)
				print "*** Done processing h7(VST,SEQ) => h5(VST,SEQ) ==> h1(VST,SEQ) ... ***"
				# if self.USEHEL:
				# 	print "*** Processing h7(VST,SEQ) => h5_UNPOL(VST,SEQ),h5_POS(VST,SEQ),h5_NEG(VST,SEQ)  ... ***"
				# 	self.proc_h5(h7,q2wbin,q2wbindir,q2wbintitle,hel=UNPOL)
				# 	self.proc_h5(h7,q2wbin,q2wbindir,q2wbintitle,hel=POS)
				# 	self.proc_h5(h7,q2wbin,q2wbindir,q2wbintitle,hel=NEG)
				# 	print "*** Done processing h7(VST,SEQ) => h5_UNPOL(VST,SEQ),h5_POS(VST,SEQ),h5_NEG(VST,SEQ)  ... ***"
				# else:
				# 	print "*** Processing h7(VST,SEQ) => h5_UNPOL(VST,SEQ) ==> h1(VST,SEQ) ... ***"
				# 	h5=self.proc_h5(h7,q2wbin,q2wbindir,q2wbintitle,hel=UNPOL)#,vst_name,vstdir)
				# 	self.proc_h1(h5,q2wbin,q2wbindir,q2wbintitle)#,vst_name,vstdir)
				# 	print "*** Done processing h7(VST,SEQ) => h5(VST,SEQ) ==> (if self.USEHEL=false) h1(VST,SEQ) ... ***"

		if que!=None:
			que.put(0)
		fout.Write()
		fout.Close()
		print "*** Done processing Crs-W bin %s ***"%(iw+1)
		return 0


	def get_h7(self,iw):
		"""
		get h7(VST,SEQ) for each Crs-W
		"""
		h7=OrderedDict()
		#! Fetch h7-ST/SR/ER directly from file (Add all tops)
		for vst in self.VSTS:
			vst_name='VST%d'%vst
			#! Fetch h7-T
			if self.DTYP=='sim':
				print "Going to get",'d2piT/h7_%d_%d'%(iw+1,vst),"..."
				h7[vst_name,'T']=self.FIN_T.Get('d2pi/h7_%d_%d'%(iw+1,vst))
				print "Done getting",'d2piT/h7_%d_%d'%(iw+1,vst)
			#! Fetch h7-R
			h7[vst_name,'R']=self.FIN_R.Get('d2pi/h7_%d_%d'%(iw+1,vst))
			# #! Fetch h7-R (Add all tops i.e. h7 = h7(t1)+h7(t2)+h7(t3)+h7(t4))
			# for itop,top in enumerate(self.TOPS):
			# 	if itop==0:
			# 		print "Going to get",'d2piR/top%d/h7_%d_%d'%(top,iw+1,vst),"..."
			# 		h7[vst_name,'R']=self.FIN_R.Get('d2piR/top%d/h7_%d_%d'%(top,iw+1,vst))
			# 		print "Done getting",'d2piR/top%d/h7_%d_%d'%(top,iw+1,vst)
			# 	else:
			# 		print "Going to get and add to top1",'d2piR/top%d/h7_%d_%d'%(top,iw+1,vst),"..."
			# 		h7[vst_name,'R'].Add(self.FIN_R.Get('d2piR/top%d/h7_%d_%d'%(top,iw+1,vst)))
			# 		print "Done gettng and adding to top1",'d2piR/top%d/h7_%d_%d'%(top,iw+1,vst)
		return h7

	def proc_h5(self,h7,q2wbin,q2wbindir,q2wbintitle):#,hel):#,vst_name,vstdir):
		"""
		+ This method is called by proc(), per Q2W-bin, after setting the appropriate Q2,W range on h7(VST,SEQ). 
		+ In this method, the following is accomplished
			+ h7(VST,SEQ) => h5(VST,SEQ)

		+ @hel=UNPOL,POS,NEG,ZERO
		"""
		h5=OrderedDict()
		print "*** Processing h7(VST,SEQ)->h5(VST,SEQ) for %s ***"%q2wbintitle

		#! seq to directly project from h7
		if self.DTYP=='exp':seq_h5=['R']
		if self.DTYP=='sim':seq_h5=['T','R']

		for vst in self.VSTS:
			vst_name='VST%d'%vst
			vstdir=q2wbindir.GetDirectory(vst_name)

			#! 1. Project out h5-T/R
			for seq in seq_h5:
				#! First set HEL
				# if hel==UNPOL:
				# 	h7[vst_name,seq].GetAxis(H7_DIM['HEL']).SetRange()
				# elif hel==POS or hel==NEG:
				# 	h7[vst_name,seq].GetAxis(H7_DIM['HEL']).SetRange(hel,hel)

				if vstdir.GetDirectory(seq)==None:vstdir.mkdir(seq).cd()
				else:vstdir.cd(seq)
				h5[vst_name,seq]=h7[vst_name,seq].Projection(5,H5_PROJDIM,"E")
				thntool.SetUnderOverFLowBinsToZero(h5[vst_name,seq]);
				h5[vst_name,seq].SetName('h5')#_%s'%HEL_NAME[hel])
				h5[vst_name,seq].SetTitle('%s_%s_%s'%(q2wbintitle,vst_name,seq))#,HEL_NAME[hel]))
				h5[vst_name,seq].Write()

			#! 2. Calculate h5-A
			if vstdir.GetDirectory('A')==None:vstdir.mkdir('A').cd()
			else:vstdir.cd('A')
			if self.DTYP=='sim':
				h5[vst_name,'A']=h5[vst_name,'R'].Clone()
				h5[vst_name,'A'].Divide(h5[vst_name,'T'])
				thntool.SetUnderOverFLowBinsToZero(h5[vst_name,'A']);
				h5[vst_name,'A'].SetName('h5')#_%s'%HEL_NAME[hel])
				h5[vst_name,'A'].SetTitle('%s_%s_%s'%(q2wbintitle,vst_name,'A'))
			if self.DTYP=='exp':
				h5[vst_name,'A']=self.FIN_SIMYIELD.Get('%s/%s/A/h5'%(q2wbin,vst_name))
			h5[vst_name,'A'].Write()

			#! 3. Calculate h5-C
			if vstdir.GetDirectory('C')==None:vstdir.mkdir('C').cd()
			else:vstdir.cd('C')
			h5[vst_name,'C']=h5[vst_name,'R'].Clone()
			h5[vst_name,'C'].Divide(h5[vst_name,'A'])
			thntool.SetUnderOverFLowBinsToZero(h5[vst_name,'C']);
			h5[vst_name,'C'].SetName('h5')#_%s'%HEL_NAME[hel])
			h5[vst_name,'C'].SetTitle('%s_%s_%s'%(q2wbintitle,vst_name,'C'))
			h5[vst_name,'C'].Write()

			#! 4. Calculate h5-H
			if vstdir.GetDirectory('H')==None:vstdir.mkdir('H').cd()
			else:vstdir.cd('H')
			if self.DTYP=='sim':
				h5[vst_name,'H']=h5[vst_name,'T'].Clone()
				h5[vst_name,'H'].Add(h5[vst_name,'C'],-1)
				thntool.SetUnderOverFLowBinsToZero(h5[vst_name,'H'])
				h5[vst_name,'H'].SetName('h5')#_%s'%HEL_NAME[hel])
				h5[vst_name,'H'].SetTitle('%s_%s_%s'%(q2wbintitle,vst_name,'H'))
			if self.DTYP=='exp':
				#!first copy sim-HOLE into exp-HOLE
				h5[vst_name,'H']=self.FIN_SIMYIELD.Get('%s/%s/H/h5'%(q2wbin,vst_name))
				#!now get SC for obtaining normalization factor
				h5_SIM_C=self.FIN_SIMYIELD.Get('%s/%s/C/h5'%(q2wbin,vst_name))
				#!calculate normalization factor
				norm=self.calcNorm(h5[vst_name,'C'],h5_SIM_C);
				#!normalize exp-HOLE
				h5[vst_name,'H'].Scale(norm);
				thntool.SetUnderOverFLowBinsToZero(h5[vst_name,'H']);
			h5[vst_name,'H'].Write()

			#! Calculate h5-F
			if vstdir.GetDirectory('F')==None:vstdir.mkdir('F').cd()
			else:vstdir.cd('F')
			h5[vst_name,'F']=h5[vst_name,'C'].Clone()
			h5[vst_name,'F'].Add(h5[vst_name,'H'],1)
			h5[vst_name,'F'].SetName('h5')#_%s'%HEL_NAME[hel])
			h5[vst_name,'F'].SetTitle('%s_%s_%s'%(q2wbintitle,vst_name,'F'))
			h5[vst_name,'F'].Write()
		print "*** done processing h7(VST,SEQ)->h5(VST,SEQ) for %s ***"%q2wbintitle
		return h5

	def proc_h1(self,h5,q2wbin,q2wbindir,q2wbintitle):#,vst_name,vstdir):
		"""
		+ This method is called by proc() after setting the appropriate Q2,W range on h7(VST,SEQ) for a q2w-bin and calling
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
			vst_name='VST%d'%vst
			vstdir=q2wbindir.GetDirectory(vst_name)
			for seq in seq_h1:
				vstdir.cd(seq)
				#! Project out h1s
				for var in VARS:
					h1=h5[vst_name,seq].Projection(H5_DIM[var],"E")
					if var=='M1' or var=='M2':
						self.setM1M2axisrange(h1,vst_name,var)
					h1.SetName('h_%s'%var)
					h1.SetTitle("%s_%s_%s_%s"%(q2wbintitle,vst_name,seq,var))
					#h1.Write() #! Am using Write() only for THnSparse Objects; for regular TH1s, Write() seems to save all cycles of an object(name;cycle)
		print "*** Processing h5(VST,SEQ)->h1(VST,SEQ) for %s ***"%q2wbintitle
		return

	def setM1M2axisrange(self,h,vst_name,var):
		if (vst_name=="VST1" and var=='M1') or (vst_name=="VST3" and var=='M2'):
			h.GetXaxis().SetRangeUser(1.1000,self.wmax-0.14)
		if (vst_name=="VST2" and var=='M2'):
			h.GetXaxis().SetRangeUser(0.2780,self.wmax-0.938)
	
	def calcNorm(self,h5D_EC,h5D_SC):
		norm,nExpEvts,nSimEvts=0,0,0
		expBinCoord=np.zeros(5,'i')
		nExpBins = h5D_EC.GetNbins()
		for iExpBin in range(nExpBins):
			nExpEvts+=h5D_EC.GetBinContent(iExpBin,expBinCoord)
			iSimBin=h5D_SC.GetBin(expBinCoord);
			nSimEvts+=h5D_SC.GetBinContent(iSimBin)
		if (nSimEvts!=0):
			norm=nExpEvts/nSimEvts;
		else:
			print "norm for Holes=0 because nSC=0!"
		return norm