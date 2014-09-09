#! Imports
from __future__ import division
from math import *

from array import *
from collections import OrderedDict

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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
	def __init__(self,dtyp,sim_num='siml',tops=[1,2,3,4],vsts=[1,2,3]):
		self.EXP,self.SIM=False,False
		if dtyp=='sim':self.SIM=True
		if dtyp=='exp':self.EXP=True
		if not(self.EXP or self.SIM):
			sys.exit("dtyp is neither EXP or SIM! Exiting")
		print "dtyp=%s"%dtyp
		self.SIM_NUM=sim_num

		self.TOPS=tops 
		self.VSTS=vsts 

		if self.EXP:
			self.DATADIR='$HOME/ongoing/mem_test/exp/new-h8-bng'
			#!self.DATADIR=os.environ['OBS_DATADIR_EXP']
			self.FIN=ROOT.TFile(os.path.join(self.DATADIR,'d2pi.root'))
			self.OUTDIR=self.DATADIR
			#!self.OUTDIR=os.path.join(os.environ['OBS_DIR'],self.SIM_NUM,self.Q2W)
			if not os.path.exists(self.OUTDIR):
				#! This path should already exist when making yield_sim
				sys.exit("Path %s does not exist. Exiting."%self.OUTDIR)
				#os.makedirs(self.OUTDIR)
			self.FIN_SIMYIELD=ROOT.TFile('$HOME/ongoing/mem_test/sim/new-h8-bng/yield_sim.root')
			#self.FIN_SIMYIELD=ROOT.TFile(os.path.join(self.OUTDIR,"yield_sim.root"))
			self.FOUT=ROOT.TFile(os.path.join(self.OUTDIR,"yield_exp.root"),"RECREATE")
			print "DATADIR=%s\nOUTDIR=%s\nFIN=%s\nFIN_SIMYIELD=%s\nFOUT=%s"%(self.DATADIR,self.OUTDIR,self.FIN.GetName(),self.FIN_SIMYIELD.GetName(),self.FOUT.GetName())
		if self.SIM:
			self.DATADIR='$HOME/ongoing/mem_test/sim/new-h8-bng'
			#self.DATADIR=os.path.join(os.environ['OBS_DATADIR_SIM'],self.SIM_NUM,self.Q2W)
			self.FIN=ROOT.TFile(os.path.join(self.DATADIR,'d2pi.root'))
			self.OUTDIR=self.DATADIR
			#!self.OUTDIR=os.path.join(os.environ['OBS_DIR'],self.SIM_NUM,self.Q2W)
			if not os.path.exists(self.OUTDIR):
				os.makedirs(self.OUTDIR)
			self.FOUT=ROOT.TFile(os.path.join(self.OUTDIR,"yield_sim.root"),"RECREATE")
			print "DATADIR=%s\nOUTDIR=%s\nFIN=%s\nFOUT=%s"%(self.DATADIR,self.OUTDIR,self.FIN.GetName(),self.FOUT.GetName())

		self.wmax=None #Used for setM1M2axisrange()
	
	def proc(self):
		#! Loop over Crw-W bins
		for iw in range(q2w_bng.NBINS_WCRS):
			print "*** Processing Crs-W bin %s ***"%(iw+1)
			#! Get h8
			h8=self.get_h8(iw)
			#! Debug h8
			for k in h8:
				print k,h8[k].GetEntries()
			#! Set up Q2,W bng
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
			# for i in range(self.Q2BNG['NBINS']):
			# 	print "Q2 bin:",i+1
			# 	print self.Q2BNG['BINS_LE'][i]
			# 	print self.Q2BNG['BINS_UE'][i]
			# 	print self.Q2BNG['BINW']
			# for i in range(self.WBNG['NBINS']):
			# 	print "W bin:",i+1
			# 	print self.WBNG['BINS_LE'][i]
			# 	print self.WBNG['BINS_UE'][i]
			# 	print self.WBNG['BINW']

			#! Loop over [Q2BNG,WBNG],VSTS,SEQ, and project: h8->h5->h1
			for i in range(self.Q2BNG['NBINS']):
				for j in range(self.WBNG['NBINS']):
					#if j>4: break
					q2wbin="%0.2f-%0.2f_%0.3f-%0.3f"%(self.Q2BNG['BINS_LE'][i],self.Q2BNG['BINS_UE'][i],self.WBNG['BINS_LE'][j],self.WBNG['BINS_UE'][j])
					q2wbindir=self.FOUT.mkdir(q2wbin)
					q2wbintitle="[%0.2f,%0.2f)_[%0.3f,%0.3f)"%(self.Q2BNG['BINS_LE'][i],self.Q2BNG['BINS_UE'][i],self.WBNG['BINS_LE'][j],self.WBNG['BINS_UE'][j])
					self.wmax=self.WBNG['BINS_UE'][j]
					#hq2w,h5,h1=OrderedDict(),OrderedDict(),OrderedDict()
					# self.hq2w.clear()
					# self.h5.clear()
					# self.h1.clear()
					print "*** Processing Q2,W %s ***"%q2wbin
					for vst in self.VSTS:
						vst_name='VST%d'%vst
						vstdir=q2wbindir.mkdir(vst_name)           
						print "*** Processing %s ***"%(vst_name)
						#! First, for h8-ST/SR/ER, set appropriate ranges for HEL,Q2&W
						print "*** h8 range setting for HEL,Q2 & W dimensions ... ***"
						if self.EXP:seq_h8=['R']
						if self.SIM:seq_h8=['T','R']
						for seq in seq_h8:
							#!-- HEL: include all helicities
							h8[vst_name,seq].GetAxis(H8_DIM['HEL']).SetRange()
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
							# outdir=os.path.join(self.OUTDIR,q2wbin,vst_name,seq)
							# if not os.path.exists(outdir):
							# 	os.makedirs(outdir)
							#! cd into FOUT.q2wbindir.vstdir.seqdir
							vstdir.mkdir(seq).cd()
							hq2w.Write()
							# cq2w=ROOT.TCanvas(hq2w.GetName(),hq2w.GetTitle())
							# hq2w.Draw("colz")
							# cq2w.SaveAs("%s/%s.png"%(outdir,cq2w.GetName()))
						print "*** Done h8 range setting"
						#! 1. h8->h5 
						h5=self.proc_h5(h8,q2wbin,q2wbindir,q2wbintitle,vst_name,vstdir)
						#! 2. h5->h1
						self.proc_h1(h5,q2wbin,q2wbindir,q2wbintitle,vst_name,vstdir)
		self.FOUT.Close()

	def get_h8(self,iw):
		"""
		get h8[iw]
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
		#global h5
		h5=OrderedDict()
		print "*** Processing h8->h5 ***"
		#! Project out h5-T/R
		if self.EXP:seq_h5=['R']
		if self.SIM:seq_h5=['T','R']
		for seq in seq_h5:
			h5[vst_name,seq]=h8[vst_name,seq].Projection(5,H5_PROJDIM,"E")
			thntool.SetUnderOverFLowBinsToZero(h5[vst_name,seq]);
			h5[vst_name,seq].SetName('h5')
			h5[vst_name,seq].SetTitle('%s_%s_%s'%(q2wbintitle,vst_name,seq))
			if vstdir.GetDirectory(seq)==None:
				vstdir.mkdir(seq).cd()
			else:
				vstdir.cd(seq)
			h5[vst_name,seq].Write()
		#! Calculate Acceptance
		if self.SIM:
			h5[vst_name,'A']=h5[vst_name,'R'].Clone()
			h5[vst_name,'A'].Divide(h5[vst_name,'T'])
			thntool.SetUnderOverFLowBinsToZero(h5[vst_name,'A']);
			h5[vst_name,'A'].SetName('h5')
			h5[vst_name,'A'].SetTitle('%s_%s_%s'%(q2wbintitle,vst_name,'A'))
			# vstdir.mkdir('A').cd()
			# h5[vst_name,'A'].Write()
		if self.EXP:
			h5[vst_name,'A']=self.FIN_SIMYIELD.Get('%s/%s/A/h5'%(q2wbin,vst_name))
		vstdir.mkdir('A').cd()
		h5[vst_name,'A'].Write()
		#! Calculate Corrected yields
		h5[vst_name,'C']=h5[vst_name,'R'].Clone()
		h5[vst_name,'C'].Divide(h5[vst_name,'A'])
		thntool.SetUnderOverFLowBinsToZero(h5[vst_name,'C']);
		h5[vst_name,'C'].SetName('h5')
		h5[vst_name,'C'].SetTitle('%s_%s_%s'%(q2wbintitle,vst_name,'C'))
		vstdir.mkdir('C').cd()
		h5[vst_name,'C'].Write()
		#! Calculate H
		if self.SIM:
			h5[vst_name,'H']=h5[vst_name,'T'].Clone()
			h5[vst_name,'H'].Add(h5[vst_name,'C'],-1)
			thntool.SetUnderOverFLowBinsToZero(h5[vst_name,'H'])
			h5[vst_name,'H'].SetName('h5')
			h5[vst_name,'H'].SetTitle('%s_%s_%s'%(q2wbintitle,vst_name,'H'))
			# vstdir.mkdir('H').cd()
			# h5[vst_name,'H'].Write()
		if self.EXP:
			#!first copy sim-HOLE into exp-HOLE
			h5[vst_name,'H']=self.FIN_SIMYIELD.Get('%s/%s/H/h5'%(q2wbin,vst_name))
			#!now get SC for obtaining normalization factor
			h5_SIM_C=self.FIN_SIMYIELD.Get('%s/%s/C/h5'%(q2wbin,vst_name))
			#!calculate normalization factor
			norm=self.calcNorm(h5[vst_name,'C'],h5_SIM_C);
			#!normalize exp-HOLE
			h5[vst_name,'H'].Scale(norm);
			thntool.SetUnderOverFLowBinsToZero(h5[vst_name,'H']);
		vstdir.mkdir('H').cd()
		h5[vst_name,'H'].Write()
		#! Calculate F
		h5[vst_name,'F']=h5[vst_name,'C'].Clone()
		h5[vst_name,'F'].Add(h5[vst_name,'H'],1)
		h5[vst_name,'F'].SetName('h5')
		h5[vst_name,'F'].SetTitle('%s_%s_%s'%(q2wbintitle,vst_name,'F'))
		vstdir.mkdir('F').cd()
		h5[vst_name,'F'].Write()
		print "*** Done processing h8->h5 ***"
		return h5

	def proc_h1(self,h5,q2wbin,q2wbindir,q2wbintitle,vst_name,vstdir):
		#h1=OrderedDict()
		print "*** Processing h5->h1 ... ***"
		if self.EXP:
			seq_h1=['R','C','A','H','F']
		if self.SIM:
			seq_h1=['T','R','C','A','H','F']
		for seq in seq_h1:
			#! Create outdir in OS (I am currently saving 1D hist plots to view with gqview)
			# outdir=os.path.join(self.OUTDIR,q2wbin,vst_name,seq)
			# if not os.path.exists(outdir):
			# 	os.makedirs(outdir)
			#! cd into FOUT.q2wbindir.vstdir.seqdir(note seqdir already create in proc_h5)
			vstdir.cd(seq)
			#! Project out h1s
			for var in VARS:
				h1=h5[vst_name,seq].Projection(H5_DIM[var],"E")
				if var=='M1' or var=='M2':
					self.setM1M2axisrange(h1,vst_name,var)
				h1.SetName('h_%s'%var)
				h1.SetTitle("%s_%s_%s_%s"%(q2wbintitle,vst_name,seq,var))
				h1.Write()
				# c1=ROOT.TCanvas(h1.GetName(),h1.GetTitle())
				# h1.Draw()
				# c1.SaveAs("%s/%s.png"%(outdir,c1.GetName()))
				# c1.Close()
		print "*** Done processing h5->h1 ***"

	def setM1M2axisrange(self,h,vst_name,var):
		if (vst_name=="VST1" and var=='M1') or (vst_name=="VST3" and var=='M2'):
			h.GetXaxis().SetRangeUser(1.1000,self.wmax-0.14)
		if (vst_name=="VST2" and var=='M2'):
			h.GetXaxis().SetRangeUser(0.2780,self.wmax-0.938)
	
	def calcNorm(self,h5D_EA,h5D_SA):
		norm,nExpEvts,nSimEvts=0,0,0
		expBinCoord=np.zeros(5,'i')
		nExpBins = h5D_EA.GetNbins()
		for iExpBin in range(nExpBins):
			nExpEvts+=h5D_EA.GetBinContent(iExpBin,expBinCoord)
			iSimBin=h5D_SA.GetBin(expBinCoord);
			nSimEvts+=h5D_SA.GetBinContent(iSimBin)
		norm=nExpEvts/nSimEvts;
		return norm