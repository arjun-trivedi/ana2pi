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

ROOT.gROOT.ProcessLine(".L THnTool.C+")
from ROOT import THnTool

# Tools
thntool=THnTool()

# Constants
#h8_dim={'HEL':0,'Q2':1,'W':2,'M1':3,'M2':4,'THETA':6,'PHI':7,'ALPHA':8}
#! See http://stackoverflow.com/questions/16553506/python-ordereddict-iteration for following code
H8_DIM=OrderedDict([('HEL',0),('Q2',1),('W',2),('M1',3),('M2',4),('THETA',5),('PHI',6),('ALPHA',7)])
H5_PROJDIM=array('i',[H8_DIM['M1'],H8_DIM['M2'],H8_DIM['THETA'],H8_DIM['PHI'],
					  H8_DIM['ALPHA']])
# print H8_DIM 
# print H8_DIM.keys()
# print H8_DIM.values()

H5_DIM=OrderedDict([('M1',0),('M2',1),('THETA',2),('PHI',3),('ALPHA',4)])
# print H5_DIM 
# print H5_DIM.keys()
# print H5_DIM.values()

# SEQ_SIM=['ST','SR','SA','SC','SH','SF']
# SEQ_EXP=['ER','EC','EH','EF']
#SEQ=['T','R','A','C','H','F']

VARS=['M1','M2','THETA','PHI','ALPHA']


EXP,SIM,SIM_NUM=False,False,None
TOPS,VSTS=None,None
Q2W,Q2BNG,WBNG=None,None,None
DATADIR,ANADIR=None,None
FIN,FOUT=None,None
FIN_SIMYIELD=None

h8=OrderedDict()
q2wbin,q2wbindir=None,None
vst_name,vstdir=None,None
hq2w,h5,h1=OrderedDict(),OrderedDict(),OrderedDict() #for each q2wbin 


def run(dtyp,sim_num='sim_1',tops=[1,2,3,4],vsts=[1,2,3],q2w='q2w2'):
	global EXP,SIM,SIM_NUM
	global TOPS,VSTS
	global Q2W,Q2BNG,WBNG
	global DATADIR,ANADIR
	global FIN,FOUT
	global FIN_SIMYIELD

	if dtyp=='sim':SIM=True
	if dtyp=='exp':EXP=True
	if not(EXP or SIM):
		sys.exit("dtyp is neither EXP or SIM! Exiting")
	print "dtyp=%s"%dtyp
	SIM_NUM='sim_1'
	TOPS=tops # as per UI
	VSTS=vsts # as per UI

	Q2W=q2w # as per UI's q2w range
	Q2BNG,WBNG=None,None
	if Q2W=='q2w2':
		q2w_bng=[[2.0,2.5,0.5],[1.3,1.9,0.025]] #min,max,width; as per UI
		Q2BNG,WBNG=atlib.init_q2wbng2(q2w_bng)
	# #! Test Q2WBNG, WBNG
	# print "Q2 bins=",Q2BNG['BINS']
	# print "W bins=",WBNG['BINS']
	# for i in range(Q2BNG['NBINS']):
	#     print "Q2 bin:",i+1
	#     print Q2BNG['BINS_LE'][i]
	#     print Q2BNG['BINS_UE'][i]
	# for i in range(WBNG['NBINS']):
	#     print "W bin:",i+1
	#     print WBNG['BINS_LE'][i]
	#     print WBNG['BINS_UE'][i]

	if EXP:
		#SEQ=SEQ_EXP
		DATADIR=os.environ['OBS_DATADIR_EXP']
		ANADIR=os.path.join(os.environ['ANA2PI_OBS'],Q2W)
		FIN=ROOT.TFile(os.path.join(DATADIR,'d2pi.root'))
		FIN_SIMYIELD=ROOT.TFile(os.path.join(ANADIR,"yield_sim.root"))
		FOUT=ROOT.TFile(os.path.join(ANADIR,"yield_exp.root"),"RECREATE")
		print "DATADIR=%s\nANADIR=%s\nFIN=%s\nFIN_SIMYIELD=%s\nFOUT=%s"%(DATADIR,ANADIR,FIN,FIN_SIMYIELD,FOUT)
	if SIM:
	#SEQ=SEQ_SIM
		DATADIR=os.environ['OBS_DATADIR_SIM']
		ANADIR=os.path.join(os.environ['ANA2PI_OBS'],Q2W)
		FIN=ROOT.TFile(os.path.join(DATADIR,Q2W,SIM_NUM,'d2pi.root'))
		FOUT=ROOT.TFile(os.path.join(ANADIR,"yield_sim.root"),"RECREATE")
		print "DATADIR=%s\nANADIR=%s\nFIN=%s\nFOUT=%s"%(DATADIR,ANADIR,FIN,FOUT)

	print "DATADIR=",DATADIR
	print "ANADIR=",ANADIR
	if not os.path.exists(ANADIR):
		ANADIR=os.makedirs(ANADIR)
	
	#! Fetch h8{}
	#h8=OrderedDict()
	get_h8()

	# #! Debug h8{}
	# for k in h8:
	#     print k,h8[k].GetEntries()   
  
	global q2wbin,q2wbindir,vst_name,vstdir
	#! Loop over [Q2BNG,WBNG],VSTS,SEQ, and project: h8->h5->h1
	for i in range(Q2BNG['NBINS']):
		for j in range(WBNG['NBINS']):
			if j>0: break
			q2wbin="%0.1f-%0.1f_%0.3f-%0.3f"%(Q2BNG['BINS_LE'][i],Q2BNG['BINS_UE'][i],WBNG['BINS_LE'][i],WBNG['BINS_UE'][i])
			q2wbindir=FOUT.mkdir(q2wbin)
			#hq2w,h5,h1=OrderedDict(),OrderedDict(),OrderedDict()
			hq2w.clear()
			h5.clear()
			h1.clear()
			for vst in VSTS:
				vst_name='VST%d'%vst
				vstdir=q2wbindir.mkdir(vst_name)           
				print "*** Processing %s,%s ***"%(q2wbin,vst_name)
				#! First, for h8-ST/SR/ER, set appropriate ranges for HEL,Q2&W
				if EXP:seq_h8=['R']
				if SIM:seq_h8=['T','R']
				for seq in seq_h8:
					#!-- HEL: include all helicities
					h8[vst_name,seq].GetAxis(H8_DIM['HEL']).SetRange()
					#!-- Q2
					q2bin_le=h8[vst_name,seq].GetAxis(H8_DIM['Q2']).FindBin(Q2BNG['BINS_LE'][i])
					q2bin_ue=h8[vst_name,seq].GetAxis(H8_DIM['Q2']).FindBin(Q2BNG['BINS_UE'][i])
					h8[vst_name,seq].GetAxis(H8_DIM['Q2']).SetRange(q2bin_le,q2bin_ue)
					#!-- W
					wbin_le=h8[vst_name,seq].GetAxis(H8_DIM['W']).FindBin(WBNG['BINS_LE'][j])
					wbin_ue=h8[vst_name,seq].GetAxis(H8_DIM['W']).FindBin(WBNG['BINS_UE'][j])
					h8[vst_name,seq].GetAxis(H8_DIM['W']).SetRange(wbin_le,wbin_ue)
					print "For h8(%s,%s),finished setting Q2-,W-bin range = [%d,%d],[%d,%d] ***"%(vst_name,seq,q2bin_le,q2bin_ue,wbin_le,wbin_ue)
					#! Project out hq2w & save (to FOUT & OS)
					hq2w[vst_name,seq]=h8[vst_name,seq].Projection(H8_DIM['Q2'],H8_DIM['W'],"E")
					hq2w[vst_name,seq].SetName('h_q2Vw')
					hq2w[vst_name,seq].SetTitle("%s_%s_%s_q2w"%(q2wbin,vst_name,seq))
					outdir=os.path.join(ANADIR,q2wbin,vst_name,seq)
					if not os.path.exists(outdir):
						os.makedirs(outdir)
					#! cd into FOUT.q2wbindir.vstdir.seqdir
					vstdir.mkdir(seq).cd()
					hq2w[vst_name,seq].Write()
					cq2w=ROOT.TCanvas(hq2w[vst_name,seq].GetName(),hq2w[vst_name,seq].GetTitle())
					hq2w[vst_name,seq].Draw("colz")
					cq2w.SaveAs("%s/%s.png"%(outdir,cq2w.GetName()))
				#! 1. h8->h5 
				proc_h5()
				#! 2. h5->h1
				proc_h1()
	FOUT.Close()
# Following defined by UI
# TOPS=[1,2,3,4] # as per UI
# VSTS=[1,2,3] # as per UI

# Q2W='q2w2' # as per UI's q2w range
# Q2BNG,WBNG=None,None
# if Q2W=='q2w2':
# 	q2w=[[2.0,2.5,0.5],[1.3,1.9,0.025]] #min,max,width; as per UI
# 	Q2BNG,WBNG=atlib.init_q2wbng2(q2w)
# # #! Test Q2WBNG, WBNG
# # print "Q2 bins=",Q2BNG['BINS']
# # print "W bins=",WBNG['BINS']
# # for i in range(Q2BNG['NBINS']):
# #     print "Q2 bin:",i+1
# #     print Q2BNG['BINS_LE'][i]
# #     print Q2BNG['BINS_UE'][i]
# # for i in range(WBNG['NBINS']):
# #     print "W bin:",i+1
# #     print WBNG['BINS_LE'][i]
# #     print WBNG['BINS_UE'][i]

# #EXP/SIM was per UI
# EXP=False
# SIM=True
# if not(EXP or SIM):
# 	sys.exit("dtyp is neither EXP or SIM! Exiting")
# SIM_NUM='sim_1'
# if EXP:
# 	#SEQ=SEQ_EXP
# 	DATADIR=os.environ['OBS_DATADIR_EXP']
# 	ANADIR=os.path.join(os.environ['ANA2PI_OBS'],Q2W)
# 	FIN=os.path.join(DATADIR,'d2pi.root')
# 	FIN_SIMYIELD=ROOT.TFile(os.path.join(ANADIR,"yield_sim.root"))
# 	FOUT=ROOT.TFile(os.path.join(ANADIR,"yield_exp.root"),"RECREATE")
# if SIM:
# 	#SEQ=SEQ_SIM
# 	DATADIR=os.environ['OBS_DATADIR_SIM']
# 	ANADIR=os.path.join(os.environ['ANA2PI_OBS'],Q2W)
# 	FIN=os.path.join(DATADIR,Q2W,SIM_NUM,'d2pi.root')
# 	FOUT=ROOT.TFile(os.path.join(ANADIR,"yield_sim.root"),"RECREATE")

# print "DATADIR=",DATADIR
# print "ANADIR=",ANADIR
# if not os.path.exists(ANADIR):
# 	ANADIR=os.makedirs(ANADIR)

def get_h8():
	#! Fetch h8-ST/SR/ER directly from file (Add all tops)
	for vst in VSTS:
		vst_name='VST%d'%vst
		#! Fetch h8-T
		if SIM:
			print "Going to get",'d2piT/yield_varset%d'%vst,"..."
			h8[vst_name,'T']=FIN.Get('d2piT/yield_varset%d'%vst)
			print "Done getting",'d2piT/yield_varset%d'%vst
		#! Fetch h8-R (Add all tops i.e. h8 = h8(t1)+h8(t2)+h8(t3)+h8(t4))
		for itop,top in enumerate(TOPS):
			if itop==0:
				print "Going to get",'d2piR/top%d/yield_varset%d'%(top,vst),"..."
				h8[vst_name,'R']=FIN.Get('d2piR/top%d/yield_varset%d'%(top,vst))
				print "Done getting",'d2piR/top%d/yield_varset%d'%(top,vst)
			else:
				print "Going to get and add to top1",'d2piR/top%d/yield_varset%d'%(top,vst),"..."
				h8[vst_name,'R'].Add(FIN.Get('d2piR/top%d/yield_varset%d'%(top,vst)))
				print "Done gettng and adding to top1",'d2piR/top%d/yield_varset%d'%(top,vst)

def proc_h5():
	#global h5
	print "*** 1. Processing h8->h5 ***"
	# #! First, for h8-ST/SR/ER, set appropriate ranges for HEL,Q2&W
	# if EXP:SEQ_h8=['R']
	# if SIM:SEQ_h8=['T','R']
	# for seq in SEQ_h8:
	# 	#!-- HEL: include all helicities
	# 	h8[vst_name,seq].GetAxis(H8_DIM['HEL']).SetRange()
	# 	#!-- Q2
	# 	q2bin_le=h8[vst_name,seq].GetAxis(H8_DIM['Q2']).FindBin(Q2BNG['BINS_LE'][i])
	# 	q2bin_ue=h8[vst_name,seq].GetAxis(H8_DIM['Q2']).FindBin(Q2BNG['BINS_UE'][i])
	# 	h8[vst_name,seq].GetAxis(H8_DIM['Q2']).SetRange(q2bin_le,q2bin_ue)
	# 	#!-- W
	# 	wbin_le=h8[vst_name,seq].GetAxis(H8_DIM['W']).FindBin(WBNG['BINS_LE'][j])
	# 	wbin_ue=h8[vst_name,seq].GetAxis(H8_DIM['W']).FindBin(WBNG['BINS_UE'][j])
	# 	h8[vst_name,seq].GetAxis(H8_DIM['W']).SetRange(wbin_le,wbin_ue)
	# 	print "For h8(%s,%s),finished setting Q2-,W-bin range = [%d,%d],[%d,%d] ***"%(vst_name,seq,q2bin_le,q2bin_ue,wbin_le,wbin_ue)
	# 	#! Project out hq2w & save (to FOUT & OS)
	# 	hq2w[vst_name,seq]=h8[vst_name,seq].Projection(H8_DIM['Q2'],H8_DIM['W'],"E")
	# 	hq2w[vst_name,seq].SetName('h_q2_v_w')
	# 	hq2w[vst_name,seq].SetTitle("%s_%s_%s_q2wbin"%(q2wbin,vst_name,seq))
	# 	outdir=os.path.join(ANADIR,q2wbin,vst_name,seq)
	# 	if not os.path.exists(outdir):
	# 		os.makedirs(outdir)
	# 	#! cd into FOUT.q2wbindir.vstdir.seqdir
	# 	vstdir.mkdir(seq).cd()
	# 	hq2w[vst_name,seq].Write()
	# 	cq2w=ROOT.TCanvas(hq2w[vst_name,seq].GetName(),hq2w[vst_name,seq].GetTitle())
	# 	hq2w[vst_name,seq].Draw("colz")
	# 	cq2w.SaveAs("%s/%s.png"%(outdir,cq2w.GetName()))
			
	#! 1.2 Finally project out h5 & save them to file
	#! Project out h5-T/R
	if EXP:SEQ_h5=['R']
	if SIM:SEQ_h5=['T','R']
	for seq in SEQ_h5:
		h5[vst_name,seq]=h8[vst_name,seq].Projection(5,H5_PROJDIM,"E")
		thntool.SetUnderOverFLowBinsToZero(h5[vst_name,seq]);
		h5[vst_name,seq].SetName('h5')
		h5[vst_name,seq].SetTitle('%s_%s_%s'%(q2wbin,vst_name,seq))
		if vstdir.GetDirectory(seq)==None:
			vstdir.mkdir(seq).cd()
		else:
			vstdir.cd(seq)
		h5[vst_name,seq].Write()
	#! Calculate Acceptance
	if SIM:
		h5[vst_name,'A']=h5[vst_name,'R'].Clone()
		h5[vst_name,'A'].Divide(h5[vst_name,'T'])
		thntool.SetUnderOverFLowBinsToZero(h5[vst_name,'A']);
		h5[vst_name,'A'].SetName('h5')
		h5[vst_name,'A'].SetTitle('%s_%s_%s'%(q2wbin,vst_name,'A'))
		# vstdir.mkdir('A').cd()
		# h5[vst_name,'A'].Write()
	if EXP:
		h5[vst_name,'A']=FIN_SIMYIELD.Get('%s/%s/A/h5'%(q2wbin,vst_name))
	vstdir.mkdir('A').cd()
	h5[vst_name,'A'].Write()
	#! Calculate Corrected yields
	h5[vst_name,'C']=h5[vst_name,'R'].Clone()
	h5[vst_name,'C'].Divide(h5[vst_name,'A'])
	thntool.SetUnderOverFLowBinsToZero(h5[vst_name,'C']);
	h5[vst_name,'C'].SetName('h5')
	h5[vst_name,'C'].SetTitle('%s_%s_%s'%(q2wbin,vst_name,'C'))
	vstdir.mkdir('C').cd()
	h5[vst_name,'C'].Write()
	#! Calculate H
	if SIM:
		h5[vst_name,'H']=h5[vst_name,'T'].Clone()
		h5[vst_name,'H'].Add(h5[vst_name,'C'],-1)
		thntool.SetUnderOverFLowBinsToZero(h5[vst_name,'H'])
		h5[vst_name,'H'].SetName('h5')
		h5[vst_name,'H'].SetTitle('%s_%s_%s'%(q2wbin,vst_name,'H'))
		vstdir.mkdir('H').cd()
		h5[vst_name,'H'].Write()
	if EXP:
		#!first copy sim-HOLE into exp-HOLE
		h5[vst_name,'H']=FIN_SIMYIELD.Get('%s/%s/H/h5'%(q2wbin,vst_name))
		#!now get SC for obtaining normalization factor
		h5_SIM_C=FIN_SIMYIELD.Get('%s/%s/C/h5'%(q2wbin,vst_name))
		#!calculate normalization factor
		norm=calcNorm(h5[vst_name,'C'],h5_SIM_C);
		#!normalize exp-HOLE
		h5[vst_name,'H'].Scale(norm);
		thntool.SetUnderOverFLowBinsToZero(h5[vst_name,'H']);
	vstdir.mkdir('H').cd()
	h5[vst_name,'H'].Write()
	#! Calculate F
	h5[vst_name,'F']=h5[vst_name,'C'].Clone()
	h5[vst_name,'F'].Add(h5[vst_name,'H'],1)
	h5[vst_name,'F'].SetName('h5')
	h5[vst_name,'F'].SetTitle('%s_%s_%s'%(q2wbin,vst_name,'F'))
	vstdir.mkdir('F').cd()
	h5[vst_name,'F'].Write()
	print "*** Done processing h8->h5 ***"

def proc_h1():
	#global h1
	print "*** Processing h5->h1 ... ***"
	if EXP:seq_h1=['R','C','A','H','F']
	if SIM:seq_h1=['T','R','C','A','H','F']
	for seq in seq_h1:
		#! Create outdir in OS (I am currently saving 1D hist plots to view with gqview)
		outdir=os.path.join(ANADIR,q2wbin,vst_name,seq)
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		#! cd into FOUT.q2wbindir.vstdir.seqdir(note seqdir already create in proc_h5)
		vstdir.cd(seq)
		#! Project out h1s
		for var in VARS:
			h1[vst_name,seq,var]=h5[vst_name,seq].Projection(H5_DIM[var],"E")
			h1[vst_name,seq,var].SetName('h_%s'%var)
			h1[vst_name,seq,var].SetTitle("%s_%s_%s_%s"%(q2wbin,vst_name,seq,var))
			h1[vst_name,seq,var].Write()
			c1=ROOT.TCanvas(h1[vst_name,seq,var].GetName(),h1[vst_name,seq,var].GetTitle())
			h1[vst_name,seq,var].Draw()
			c1.SaveAs("%s/%s.png"%(outdir,c1.GetName()))
			c1.Close()
	print "*** Done processing h5->h1 ***"

def calcNorm(h5D_EA,h5D_SA):
	norm,nExpEvts,nSimEvts=0,0,0
	expBinCoord=np.zeros(5,'i')
	nExpBins = h5D_EA.GetNbins()
	for iExpBin in range(nExpBins):
		nExpEvts+=h5D_EA.GetBinContent(iExpBin,expBinCoord)
		iSimBin=h5D_SA.GetBin(expBinCoord);
		nSimEvts+=h5D_SA.GetBinContent(iSimBin)
	norm=nExpEvts/nSimEvts;
	return norm

# <codecell>

# #! Fetch h8{}
# h8=OrderedDict()
# f=ROOT.TFile(FIN)
# get_h8()

# # #! Debug h8{}
# # for k in h8:
# #     print k,h8[k].GetEntries()   
  

# # <codecell>

# #! Loop over [Q2BNG,WBNG],VSTS,SEQ, and project: h8->h5->h1
# FOUT=ROOT.TFile(os.path.join(ANADIR,"yield_sim.root"),"RECREATE")
# for i in range(Q2BNG['NBINS']):
# 	for j in range(WBNG['NBINS']):
# 		if j>0: break
# 		q2wbin="%0.1f-%0.1f_%0.3f-%0.3f"%(Q2BNG['BINS_LE'][i],Q2BNG['BINS_UE'][i],WBNG['BINS_LE'][i],WBNG['BINS_UE'][i])
# 		q2wbindir=FOUT.mkdir(q2wbin)
# 		hq2w,h5,h1=OrderedDict(),OrderedDict(),OrderedDict()
# 		for vst in VSTS:
# 			vst_name='VST%d'%vst
# 			vstdir=q2wbindir.mkdir(vst_name)           
# 			print "*** Processing %s,%s ***"%(q2wbin,vst_name)
# 			#! 1. h8->h5 
# 			proc_h5()
# 			#! 2. h5->h1
# 			proc_h1()
# FOUT.Close()

# <codecell>


