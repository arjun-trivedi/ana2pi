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
H8_DIM=OrderedDict([('HEL',0),('Q2',1),('W',2),('M1',3),('M2',4),('THETA',5),('PHI',6),('ALPHA',7)])
H5_PROJDIM=array('i',[H8_DIM['M1'],H8_DIM['M2'],H8_DIM['THETA'],H8_DIM['PHI'],
					  H8_DIM['ALPHA']])
	
H5_DIM=OrderedDict([('M1',0),('M2',1),('THETA',2),('PHI',3),('ALPHA',4)])
VARS=['M1','M2','THETA','PHI','ALPHA']

class ProcYields:
	"""
	+ Accomplishes
		+ d2pi.root => yield_exp/sim.root
		+ d2pi.root => yield_exp_hel.root 

	+ The direct interface for the user is the 'execute()' method that calls 'proc()' for each h8[W]. The 'execute()'
	method uses Python's multiprocessing module to accomplish the task of processing each h8[W] in 
	a separate thread, thereby not having to load all h8[W] at once and exhausting the system's memory (atrivedi-
	laptop; 12GB RAM)

	+ For each q2w-bin and Varset therein, 'proc()' calls 'proc_h5()' and in the case for processing
	h8s irrespective of the beam helicity, calls 'proc_h1()' too.

	"""
	def __init__(self,dtyp,simnum='siml',tops=[1,2,3,4],vsts=[1,2,3]):
		self.EXP,self.SIM=False,False
		if dtyp=='sim':self.SIM=True
		if dtyp=='exp':self.EXP=True
		if not(self.EXP or self.SIM):
			sys.exit("dtyp is neither EXP or SIM! Exiting")
		print "dtyp=%s"%dtyp
		self.SIMNUM=simnum

		self.TOPS=tops 
		self.VSTS=vsts 

		if self.EXP:
			self.DATADIR=os.environ['D2PIDIR_EXP']
			self.FIN=ROOT.TFile(os.path.join(self.DATADIR,'d2pi.root'))
			self.OUTDIR=os.path.join(os.environ['OBSDIR'],self.SIMNUM)
			if not os.path.exists(self.OUTDIR):
				#! This path should already exist when making yield_exp
				sys.exit("Path %s does not exist. Exiting."%self.OUTDIR)
			self.FIN_SIMYIELD=ROOT.TFile(os.path.join(self.OUTDIR,"yield_sim.root"))
			self.FOUT=ROOT.TFile(os.path.join(self.OUTDIR,"yield_exp.root"),"RECREATE")
			print "DATADIR=%s\nOUTDIR=%s\nFIN=%s\nFIN_SIMYIELD=%s\nFOUT=%s"%(self.DATADIR,self.OUTDIR,self.FIN.GetName(),self.FIN_SIMYIELD.GetName(),self.FOUT.GetName())
		if self.SIM:
			self.DATADIR=os.path.join(os.environ['D2PIDIR_SIM'],self.SIMNUM)
			self.FIN=ROOT.TFile(os.path.join(self.DATADIR,'d2pi.root'))
			self.OUTDIR=os.path.join(os.environ['OBSDIR'],self.SIMNUM)
			if not os.path.exists(self.OUTDIR):
				os.makedirs(self.OUTDIR)
			self.FOUT=ROOT.TFile(os.path.join(self.OUTDIR,"yield_sim.root"),"RECREATE")
			print "DATADIR=%s\nOUTDIR=%s\nFIN=%s\nFOUT=%s"%(self.DATADIR,self.OUTDIR,self.FIN.GetName(),self.FOUT.GetName())

		self.wmax=None #Used for setM1M2axisrange()

	def execute(self):
		#! Loop over Crw-W bins

		queue=Queue()
		for iw in range(q2w_bng.NBINS_WCRS):
		#for iw in range(2):
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
		1. Get h8(VST,SEQ) per Crs-W from FIN
		2. Open FOUT in "appropriate" mode
		3. Setup Q2W binning as per Crs-W bin
		4. In each Q2W-bin and for each VST therein:
			4.i   Set appropriate Q2 & W projection range for h8(VST,SEQ)
			4.ii. proc_h5() & if self.USEHEL=true, proc_h1() 
				
		"""
		print "*** Processing Crs-W bin %s ***"%(iw+1)

		#! 1. Get h8s from FIN (d2pi.root)
		h8=self.get_h8(iw)
		#! Debug h8
		for k in h8:
			print k,h8[k].GetEntries()

		#! 2. Prepare FOUT (yield_exp/sim.root)
		if self.EXP:
			fname='yield_exp.root'
		if self.SIM:
			fname='yield_sim.root'
		if iw==0:
			f=ROOT.TFile(os.path.join(self.OUTDIR,fname),"RECREATE")
		else:
			f=ROOT.TFile(os.path.join(self.OUTDIR,fname),"UPDATE")
		
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
		
		#! 4. Loop over Q2W-bins and VSTs there in
		for i in range(self.Q2BNG['NBINS']):
			for j in range(self.WBNG['NBINS']):
				#if j>4: break
				q2wbin="%0.2f-%0.2f_%0.3f-%0.3f"%(self.Q2BNG['BINS_LE'][i],self.Q2BNG['BINS_UE'][i],self.WBNG['BINS_LE'][j],self.WBNG['BINS_UE'][j])
				q2wbindir=f.mkdir(q2wbin)
				q2wbintitle="[%0.2f,%0.2f)_[%0.3f,%0.3f)"%(self.Q2BNG['BINS_LE'][i],self.Q2BNG['BINS_UE'][i],self.WBNG['BINS_LE'][j],self.WBNG['BINS_UE'][j])
				self.wmax=self.WBNG['BINS_UE'][j]
				print "*** Processing Q2,W %s ***"%q2wbin
				for vst in self.VSTS:
					vst_name='VST%d'%vst
					vstdir=q2wbindir.mkdir(vst_name)           
					print "*** Processing %s ***"%(vst_name)
					#! 4.i. Set appropriate Q2 & W projection range for h8(VST,SEQ=ST/SR/ER) 
					print "*** h8 range setting for HEL,Q2 & W dimensions ... ***"
					if self.EXP:seq_h8=['R']
					if self.SIM:seq_h8=['T','R']
					for seq in seq_h8:
						#!-- Q2
						q2bin_le=h8[vst_name,seq].GetAxis(H8_DIM['Q2']).FindBin(self.Q2BNG['BINS_LE'][i]+self.Q2BNG['BINW']/2)
						q2bin_ue=h8[vst_name,seq].GetAxis(H8_DIM['Q2']).FindBin(self.Q2BNG['BINS_UE'][i]-self.Q2BNG['BINW']/2)
						h8[vst_name,seq].GetAxis(H8_DIM['Q2']).SetRange(q2bin_le,q2bin_ue)
						#!-- W
						wbin_le=h8[vst_name,seq].GetAxis(H8_DIM['W']).FindBin(self.WBNG['BINS_LE'][j]+self.WBNG['BINW']/2)
						wbin_ue=h8[vst_name,seq].GetAxis(H8_DIM['W']).FindBin(self.WBNG['BINS_UE'][j]-self.WBNG['BINW']/2)
						h8[vst_name,seq].GetAxis(H8_DIM['W']).SetRange(wbin_le,wbin_ue)
						print "For h8(%s),finished setting range for Q2-,W-bin = %s ***"%(seq,q2wbintitle)
						#! Project out hq2w & save (to FOUT & OS)
						hq2w=h8[vst_name,seq].Projection(H8_DIM['Q2'],H8_DIM['W'],"E")
						hq2w.SetName('h_q2Vw')
						hq2w.SetTitle("%s_%s_%s_q2w"%(q2wbintitle,vst_name,seq))
						vstdir.mkdir(seq).cd()
						hq2w.Write()
					print "*** Done h8 range setting"
					#! 4.ii. proc_h5() & if self.USEHEL=true, proc_h1() 
					#! h8->h5 
					h5=self.proc_h5(h8,q2wbin,q2wbindir,q2wbintitle,vst_name,vstdir)
					#! h5->h1
					self.proc_h1(h5,q2wbin,q2wbindir,q2wbintitle,vst_name,vstdir)
		if que!=None:
			que.put(0)
		f.Write()
		f.Close()
		print "*** Done processing Crs-W bin %s ***"%(iw+1)
		return 0


	def get_h8(self,iw):
		"""
		get h8(VST,SEQ) for each Crs-W
		"""
		h8=OrderedDict()
		#! Fetch h8-ST/SR/ER directly from file (Add all tops)
		for vst in self.VSTS:
			vst_name='VST%d'%vst
			#! Fetch h8-T
			if self.SIM:
				print "Going to get",'d2piT/h8_%d_%d'%(iw+1,vst),"..."
				h8[vst_name,'T']=self.FIN.Get('d2piT/h8_%d_%d'%(iw+1,vst))
				print "Done getting",'d2piT/h8_%d_%d'%(iw+1,vst)
			#! Fetch h8-R (Add all tops i.e. h8 = h8(t1)+h8(t2)+h8(t3)+h8(t4))
			for itop,top in enumerate(self.TOPS):
				if itop==0:
					print "Going to get",'d2piR/top%d/h8_%d_%d'%(top,iw+1,vst),"..."
					h8[vst_name,'R']=self.FIN.Get('d2piR/top%d/h8_%d_%d'%(top,iw+1,vst))
					print "Done getting",'d2piR/top%d/h8_%d_%d'%(top,iw+1,vst)
				else:
					print "Going to get and add to top1",'d2piR/top%d/h8_%d_%d'%(top,iw+1,vst),"..."
					h8[vst_name,'R'].Add(self.FIN.Get('d2piR/top%d/h8_%d_%d'%(top,iw+1,vst)))
					print "Done gettng and adding to top1",'d2piR/top%d/h8_%d_%d'%(top,iw+1,vst)
		return h8

	def proc_h5(self,h8,q2wbin,q2wbindir,q2wbintitle,vst_name,vstdir):
		h5=OrderedDict()
		print "*** Processing h8->h5 ***"

		#! Project out h5-T/R
		if self.EXP:
			seq_h5=['R']
		if self.SIM:
			seq_h5=['T','R']
		for seq in seq_h5:
			#! First set HEL: include all helicities
			h8[vst_name,seq].GetAxis(H8_DIM['HEL']).SetRange()
			
			h5[seq]=h8[vst_name,seq].Projection(5,H5_PROJDIM,"E")
			thntool.SetUnderOverFLowBinsToZero(h5[seq]);
			h5[seq].SetName('h5')
			h5[seq].SetTitle('%s_%s_%s'%(q2wbintitle,vst_name,seq))
			if vstdir.GetDirectory(seq)==None:
				vstdir.mkdir(seq).cd()
			else:
				vstdir.cd(seq)
			h5[seq].Write()

		#! Calculate Acceptance
		if self.SIM:
			h5['A']=h5['R'].Clone()
			h5['A'].Divide(h5['T'])
			thntool.SetUnderOverFLowBinsToZero(h5['A']);
			h5['A'].SetName('h5')
			h5['A'].SetTitle('%s_%s_%s'%(q2wbintitle,vst_name,'A'))
			# vstdir.mkdir('A').cd()
			# h5['A'].Write()
		if self.EXP:
			h5['A']=self.FIN_SIMYIELD.Get('%s/%s/A/h5'%(q2wbin,vst_name))
		vstdir.mkdir('A').cd()
		h5['A'].Write()
		#! Calculate Corrected yields
		h5['C']=h5['R'].Clone()
		h5['C'].Divide(h5['A'])
		thntool.SetUnderOverFLowBinsToZero(h5['C']);
		h5['C'].SetName('h5')
		h5['C'].SetTitle('%s_%s_%s'%(q2wbintitle,vst_name,'C'))
		vstdir.mkdir('C').cd()
		h5['C'].Write()
		#! Calculate H
		if self.SIM:
			h5['H']=h5['T'].Clone()
			h5['H'].Add(h5['C'],-1)
			thntool.SetUnderOverFLowBinsToZero(h5['H'])
			h5['H'].SetName('h5')
			h5['H'].SetTitle('%s_%s_%s'%(q2wbintitle,vst_name,'H'))
			# vstdir.mkdir('H').cd()
			# h5['H'].Write()
		if self.EXP:
			#!first copy sim-HOLE into exp-HOLE
			h5['H']=self.FIN_SIMYIELD.Get('%s/%s/H/h5'%(q2wbin,vst_name))
			#!now get SC for obtaining normalization factor
			h5_SIM_C=self.FIN_SIMYIELD.Get('%s/%s/C/h5'%(q2wbin,vst_name))
			#!calculate normalization factor
			norm=self.calcNorm(h5['C'],h5_SIM_C);
			#!normalize exp-HOLE
			h5['H'].Scale(norm);
			thntool.SetUnderOverFLowBinsToZero(h5['H']);
		vstdir.mkdir('H').cd()
		h5['H'].Write()
		#! Calculate F
		h5['F']=h5['C'].Clone()
		h5['F'].Add(h5['H'],1)
		h5['F'].SetName('h5')
		h5['F'].SetTitle('%s_%s_%s'%(q2wbintitle,vst_name,'F'))
		vstdir.mkdir('F').cd()
		h5['F'].Write()
		print "*** Done processing h8->h5 ***"
		return h5

	def proc_h1(self,h5,q2wbin,q2wbindir,q2wbintitle,vst_name,vstdir):
		print "*** Processing h5->h1 ... ***"
		if self.EXP:
			seq_h1=['R','C','A','H','F']
		if self.SIM:
			seq_h1=['T','R','C','A','H','F']
		for seq in seq_h1:
			vstdir.cd(seq)
			#! Project out h1s
			for var in VARS:
				h1=h5[seq].Projection(H5_DIM[var],"E")
				if var=='M1' or var=='M2':
					self.setM1M2axisrange(h1,vst_name,var)
				h1.SetName('h_%s'%var)
				h1.SetTitle("%s_%s_%s_%s"%(q2wbintitle,vst_name,seq,var))
				h1.Write()
		print "*** Done processing h5->h1 ***"

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