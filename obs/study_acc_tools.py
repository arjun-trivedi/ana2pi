from __future__ import division
import ROOT
import os
from array import *
from collections import OrderedDict
import numpy as np
from rootpy.interactive import wait
import matplotlib.pyplot as plt
import math

def plot_intg_acc(obsdate,expt,sim,tops,q2binl,vst="VST1"):
	'''
	+ Plot integrated acceptance in Mpippim dimension in W-bins defined in WBINL for various Q2-bins (input by user).
	+ Starting point are h5T and h5R located in $OBSDIR/simX

	Input paramaters:
	+ obsdate="obs_<date>""
	+ expt="e1f" or "e16" (default="e1f")
	+ sim=simnum entered as string. Ex. "sim1","siml" (default="siml")
	+ tops=List of topologies. Ex. [1,2,3,4]
	+ vst=varset entered as string. Ex: "VST1"
	+ q2binl=List of q2bins entered as string. Ex: ["1.25-1.75","1.75-2.25"]
	'''
	if expt=='e1f':
		DATADIR=os.path.join(os.environ['OBSDIR'],obsdate)
	elif expt=='e16':
		DATADIR=os.path.join(os.environ['OBSDIR_E16'],obsdate)
	else:
		sys.exit("expt is neither E1F or E16! Exiting")

	FIN=ROOT.TFile(os.path.join(DATADIR,sim,'yield_sim_top%s.root'%''.join(str(t) for t in tops)))
	OUTDIR=os.path.join(DATADIR,sim,'intg-eff_top%s.root'%''.join(str(t) for t in tops))
	if not os.path.exists(OUTDIR):
    		os.makedirs(OUTDIR)
	WBINL=["1.575-1.600","1.675-1.700","1.775-1.800","1.875-1.900","1.975-2.000"]
	for wbin in WBINL:
		q2wbinl=[]
		for q2bin in q2binl:
			q2wbinl.append("%s_%s"%(q2bin,wbin))
    		c=ROOT.TCanvas("c","c",800,700)
    		c.Divide(2,2)
    		coll=[ROOT.gROOT.ProcessLine("kBlue"),ROOT.gROOT.ProcessLine("kRed")]
    		l=ROOT.TLegend(0.85,0.5,1.0,0.65)
		hT,hR,hA=[],[],[]
    		for iq2wbin,q2wbin in enumerate(q2wbinl):
        		print "iq2wbin=",iq2wbin
        		h5T=FIN.Get("%s/%s/T/h5_UNPOL"%(q2wbin,vst))
        		h5R=FIN.Get("%s/%s/R/h5_UNPOL"%(q2wbin,vst))
        		hT.append(h5T.Projection(1,"E"))
        		hR.append(h5R.Projection(1,"E"))
        		hA.append(hR[iq2wbin].Clone("hA%d"%iq2wbin))
        		hA[iq2wbin].Divide(hT[iq2wbin])
        		dopt=""
        		if iq2wbin==1:
				dopt="sames"
        		c.cd(1)
        		hT[iq2wbin].Draw(dopt)
        		hT[iq2wbin].SetLineColor(coll[iq2wbin])
        		l.AddEntry(hT[iq2wbin],q2wbin.split("_")[0])
        		c.cd(2)
        		hR[iq2wbin].Draw(dopt)
        		hR[iq2wbin].SetLineColor(coll[iq2wbin])
        		c.cd(3)
        		if iq2wbin==0:
            			hA[iq2wbin].SetMinimum(0.)
            			hA[iq2wbin].SetMaximum(0.25)
        		hA[iq2wbin].Draw(dopt)
        		hA[iq2wbin].SetLineColor(coll[iq2wbin])
    		c.cd(1)
    		l.Draw()
    		c.SaveAs("%s/ceff_%s_%s.png"%(OUTDIR,wbin,vst))


def get_div_error(N,D):
    '''
    + Caculates relative error in Q=N/D.
    '''
    errN=math.sqrt(N)
    errD=math.sqrt(D)
    
    relErrN=errN/N 
    relErrD=errD/D
    
    relErrQ=math.sqrt(math.pow(relErrN,2)+math.pow(relErrD,2))
    return relErrQ

def plot_acc_q2wbin(sim,top=2,vst=1,q2wbin="1.25-1.75_1.575-1.600",crs_w_bin=3):
	'''
	For a given simX, the idea is to extract 5D acceptance in a given q2-w bin,starting from h8 and to display it
        in a meaningful way.
	'''
	ROOT.gROOT.ProcessLine(".L THnTool.C+")
	from ROOT import THnTool
	# Tools 
	thntool=THnTool()

	H8_DIM=OrderedDict([('HEL',0),('Q2',1),('W',2),('M1',3),('M2',4),('THETA',5),('PHI',6),('ALPHA',7)])
	H5_PROJDIM=array('i',[H8_DIM['M1'],H8_DIM['M2'],H8_DIM['THETA'],H8_DIM['PHI'],
					  H8_DIM['ALPHA']])
	
	H5_DIM=OrderedDict([('M1',0),('M2',1),('THETA',2),('PHI',3),('ALPHA',4)])
	VARS=['M1','M2','THETA','PHI','ALPHA']

	FINR=ROOT.TFile(os.path.join(os.environ['D2PIDIR_SIM'],sim,'d2piR.root'))
	FINT=ROOT.TFile(os.path.join(os.environ['D2PIDIR_SIM'],sim,'d2piT.root'))
	print "FINR=",FINR.GetName()
	print "FINT=",FINT.GetName()

	OUTDIR=os.path.join(os.environ['OBSDIR'],sim,'acc_debug')
	if not os.path.exists(OUTDIR):
    		os.makedirs(OUTDIR)
	h8R=FINR.Get('d2piR/top%d/h8_%d_%d'%(top,crs_w_bin,vst))
	h8T=FINT.Get('d2piT/h8_%d_%d'%(crs_w_bin,vst))
	q2min=float(q2wbin.split("_")[0].split("-")[0])
	q2max=float(q2wbin.split("_")[0].split("-")[1])
	wmin=float(q2wbin.split("_")[1].split("-")[0])
        wmax=float(q2wbin.split("_")[1].split("-")[1])
	#! --Q2 bin range
	q2bin_le=h8R.GetAxis(H8_DIM['Q2']).FindBin(q2min+0.5/2)
	q2bin_ue=h8R.GetAxis(H8_DIM['Q2']).FindBin(q2max-0.5/2)
	h8R.GetAxis(H8_DIM['Q2']).SetRange(q2bin_le,q2bin_ue)
	q2bin_le=h8T.GetAxis(H8_DIM['Q2']).FindBin(q2min+0.5/2)
	q2bin_ue=h8T.GetAxis(H8_DIM['Q2']).FindBin(q2max-0.5/2)
	h8T.GetAxis(H8_DIM['Q2']).SetRange(q2bin_le,q2bin_ue)
	#!-- W bin range
	wbin_le=h8R.GetAxis(H8_DIM['W']).FindBin(wmin+0.025/2)
	wbin_ue=h8R.GetAxis(H8_DIM['W']).FindBin(wmax-0.025/2)
	h8R.GetAxis(H8_DIM['W']).SetRange(wbin_le,wbin_ue)
	wbin_le=h8T.GetAxis(H8_DIM['W']).FindBin(wmin+0.025/2)
	wbin_ue=h8T.GetAxis(H8_DIM['W']).FindBin(wmax-0.025/2)
	h8T.GetAxis(H8_DIM['W']).SetRange(wbin_le,wbin_ue)
	cq2w=ROOT.TCanvas("Q2,W")
	cq2w.Divide(2,1)
	cq2w.cd(1)
	hq2wT=h8T.Projection(H8_DIM['Q2'],H8_DIM['W'],"E")
	hq2wT.Draw("colz")
	cq2w.cd(2)
	hq2wR=h8R.Projection(H8_DIM['Q2'],H8_DIM['W'],"E")
	hq2wR.Draw("colz")
	#cq2w.SaveAs("%s/c_q2w.png"%OUTDIR)
	#cq2w.Close()

	#! Get 5D Acceptance
	h5T=h8T.Projection(5,H5_PROJDIM,"E")
	h5R=h8R.Projection(5,H5_PROJDIM,"E")
	h5A=h5R.Clone()
	h5A.Divide(h5T)
	#binc_stats=np.zeros(2,'f')
	#thntool.GetBinContentDistStats(h5A,binc_stats)
	#mu=binc_stats[0]
	#sg=binc_stats[1]
	#print "mu,sg=",mu,sg

	#h5T_dist=thntool.GetBinContentDist(h5T,400,0,400)
	#h5R_dist=thntool.GetBinContentDist(h5R,60,0,60)
	#h5A_dist=thntool.GetBinContentDist(h5A,100,0,1)
	#c_binc_dist=ROOT.TCanvas()

	#c_binc_dist.Divide(2,2)
	#c_binc_dist.cd(1)
	#h5T_dist.Draw()
	#c_binc_dist.cd(2)
	#h5R_dist.Draw()
	#c_binc_dist.cd(3)
	#h5A_dist.Draw()
	#c_binc_dist.SaveAs("%s/c_binc_dist.png"%OUTDIR)

	#f=open("acc.out","w")
	nbins=h5T.GetNbins()
	hT=ROOT.TH1F("hT","hT",400,0,400)
	hR=ROOT.TH1F("hR","hR",60,0,60)
	hA=ROOT.TH1F("hA","hA",100,0,1)
	h_relErrA=ROOT.TH1F("hA","hA",100,0,1)
	#print nbins
	bincoordT=np.zeros(5,'i')
	for ibin in range(nbins):
    		bincT=h5T.GetBinContent(ibin,bincoordT)
    
    		binR=h5R.GetBin(bincoordT)
    		bincR=h5R.GetBinContent(binR)
       
    		binA=h5A.GetBin(bincoordT)
    		bincA=h5A.GetBinContent(binA)
		if bincA>0:
			rel_err_A=h5A.GetBinError(binA)/bincA
			#rel_err_A=get_div_error(bincR,bincT)
    
    		if bincA>0 and bincA<1:
        		hT.Fill(bincT)
        		hR.Fill(bincR)
        		hA.Fill(bincA)
        		h_relErrA.Fill(rel_err_A)
	c=ROOT.TCanvas()
	c.Divide(2,2)
	c.cd(1)
	hT.Draw()
	c.cd(2)
	hR.Draw()
	c.cd(3)
	hA.Draw()
	c.cd(4)
	h_relErrA.Draw()

	if not ROOT.gROOT.IsBatch():
    		plt.show()
    		# wait for you to close the ROOT canvas before exiting
    		wait(True)
def plot_acc_q2wbin_5D(obsdate,expt='e1f',sim='siml',q2wbin="1.75-2.25_1.575-1.600",tops=[1,2,3,4],vst=1):
	'''
	For a given simX, display acceptance in a given q2-w bin

	Arguments:
	-obsdate="obs_<date>"
	-expt="e1f" or "e16" (default='e1f')
	-sim="simnum" (default='siml')
	-q2wbin="q2min-q2max_wmin-wmax" (default="1.75-2.25_1.575-1.600")
	-tops=list of tops (default=[1,2,3,4])
	-vst=vst number (default=1)
	'''
	ROOT.gROOT.ProcessLine(".L THnTool.C+")
	from ROOT import THnTool
	# Tools 
	thntool=THnTool()

	if expt=='e1f':
		DATADIR=os.path.join(os.environ['OBSDIR'],obsdate)
	elif expt=='e16':
		DATADIR=os.path.join(os.environ['OBSDIR_E16'],obsdate)
	else:
		sys.exit("expt is neither E1F or E16! Exiting")

	FIN=ROOT.TFile(os.path.join( DATADIR,sim,"yield_sim_top%s.root"%(''.join(str(t) for t in tops)) ))
	print "FIN=",FIN.GetName()

	#! Get 5D histograms
	h5T=FIN.Get('%s/VST%d/T/h5_UNPOL'%(q2wbin,vst))
	h5R=FIN.Get('%s/VST%d/R/h5_UNPOL'%(q2wbin,vst))
	h5A=FIN.Get('%s/VST%d/A/h5_UNPOL'%(q2wbin,vst))
	
	nbins=h5T.GetNbins()
	hT=ROOT.TH1F("hT","hT",400,0,400)
	hR=ROOT.TH1F("hR","hR",60,0,60)
	hA=ROOT.TH1F("hA","hA",100,0,1)
	h_relErrA=ROOT.TH1F("hA_rel_err","hA_rel_err",100,0,1)
	#print nbins
	bincoordT=np.zeros(5,'i')
	for ibin in range(nbins):
    		bincT=h5T.GetBinContent(ibin,bincoordT)
    
    		binR=h5R.GetBin(bincoordT)
    		bincR=h5R.GetBinContent(binR)
       
    		binA=h5A.GetBin(bincoordT)
    		bincA=h5A.GetBinContent(binA)
		if bincA>0:
			rel_err_A=h5A.GetBinError(binA)/bincA
			#rel_err_A=get_div_error(bincR,bincT)
    
    		if bincA>0 and bincA<1:
        		hT.Fill(bincT)
        		hR.Fill(bincR)
        		hA.Fill(bincA)
        		h_relErrA.Fill(rel_err_A)
	c=ROOT.TCanvas()
	c.Divide(2,2)
	c.cd(1)
	hT.Draw()
	c.cd(2)
	hR.Draw()
	c.cd(3)
	hA.Draw()
	c.cd(4)
	h_relErrA.Draw()

	if not ROOT.gROOT.IsBatch():
    		plt.show()
    		# wait for you to close the ROOT canvas before exiting
    		wait(True)
