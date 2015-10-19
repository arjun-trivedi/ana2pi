from __future__ import division
import ROOT
from rootpy.io import root_open, DoesNotExist

from collections import OrderedDict
import array
import itertools

import os,sys

import numpy as np
import matplotlib.pyplot as plt

import math

import atlib as atlib

from proc_h8 import H5_DIM,VSTS,VARS,SEQ_ALL 

ROOT.gROOT.ProcessLine(".L THnTool.C+")
from ROOT import THnTool

# Tools 
thntool=THnTool()

#! For Luminosity & vgflux
LUM=19.844 #fb^-1
LUM_INVFB_TO_INVMICROB=1000000000

PI=3.14159265358979312
E1F_E0=5.499
FSC=0.00729735253
A=FSC
NA=6.02214129E23
QE=1.60217646E-19
MP=0.93827203

#! Create dictionary for VAR_NAMES(VST,VAR): 
#! + VAR_NAMES(VST,VAR) should correspond with 'DataAna::MakeYields()':
#! 		+VST1: Delta++ => (p,pip)pim
#! 		+VST2: Rho     => (pip,pim)p
#! 		+VST3: Delta0  => (p,pim)pip
VAR_NAMES={(1,VARS[0]):"M_{p#pi^{+}}",(1,VARS[1]):"M_{#pi^{+}#pi^{-}}",
           (1,VARS[2]):"#theta_{#pi^{-}}",(1,VARS[3]):"#phi_{#pi^{-}}",
           (1,VARS[4]):"#alpha_{[p_{f}#pi^{+}][p#pi^{-}]}",
           (2,VARS[0]):"M_{p#pi^{+}}",(2,VARS[1]):"M_{#pi^{+}#pi^{-}}",
           (2,VARS[2]):"#theta_{p}",(2,VARS[3]): "#phi_{p}",
           (2,VARS[4]):"#alpha_{[#pi^{+}#pi^{-}][pp_{f}]}",
           (3,VARS[0]):"M_{p#pi^{+}}",(3,VARS[1]):"M_{p#pi^{-}}",
           (3,VARS[2]):"#theta_{#pi^{+}}",(3,VARS[3]): "#phi_{#pi^{+}}",
           (3,VARS[4]):"#alpha_{[p_{f}#pi^{-}][p#pi^{+}]}"}

VAR_UNIT_NAMES={VARS[0]:"[GeV]",VARS[1]:"[GeV]",VARS[2]:"",VARS[3]:"",VARS[4]:""}

#! Following dictionaries for extracting R2
R2_NAMED={'A':'R2_{T}+R2_{L}','B':'R2_{LT}','C':'R2_{TT}','D':'R2_{LT^{\'}}'}
H5_MFUNCD={'A':'1','B':'cphi','C':'c2phi','D':'sphi'}
MTHD_NAMED={'mthd1':'h5-mply-itg','mthd2':'phi-proj-fit','mthd3':'phi-prof-mply-itg'}
#! The following binning information taken from h8_bng.h
NBINS={'M1':14,'M2':14,'THETA':10,'PHI':10,'ALPHA':10}	
#! The following for phi-projection canvas
NXPADS={'M1':2,'M2':2,'THETA':2,'PHI':2,'ALPHA':2}
NYPADS={'M1':7,'M2':7,'THETA':5,'PHI':5,'ALPHA':5}
#! Fit function (needed for mthd2)
FPHI=ROOT.TF1("fphi", "([0] + [1]*cos(x*TMath::DegToRad()) + [2]*cos(2*x*TMath::DegToRad()) + [3]*sin(x*TMath::DegToRad()))",0,360)
FPHI.SetParameter(0,1)
FPHI.SetParameter(1,10)
FPHI.SetParameter(2,20)
FPHI.SetParameter(3,100)
FPHI.SetParName(0, "A")
FPHI.SetParName(1, "B")
FPHI.SetParName(2, "C")
FPHI.SetParName(3, "D")#hPD

#! Constants relating to plotting and its aesthetics
PLT_SEQ_ALL=SEQ_ALL#['ST','SR','SA','SC','SH','SF','ER','EC','EH','EF']
PLT_SEQ_ALL.remove('SA') #! SA is not plotted, at least not as an Observable

COLS_PLT_SEQ_ALL={('ST'):ROOT.gROOT.ProcessLine("kGreen"),
			      ('SR'):ROOT.gROOT.ProcessLine("kMagenta"),
			      ('SC'):ROOT.gROOT.ProcessLine("kOrange"),
			      ('SH'):ROOT.gROOT.ProcessLine("kPink+1"),
			      ('SF'):ROOT.gROOT.ProcessLine("kRed"),
                  ('ER'):ROOT.gROOT.ProcessLine("kYellow"),
			      ('EC'):ROOT.gROOT.ProcessLine("kCyan"),
			      ('EH'):ROOT.gROOT.ProcessLine("kBlack"),
			      ('EF'):ROOT.gROOT.ProcessLine("kBlue")}

MRKS_PLT_SEQ_ALL={('ST'):ROOT.gROOT.ProcessLine("kPlus"),
			      ('SR'):ROOT.gROOT.ProcessLine("kOpenStar"),
			      ('SC'):ROOT.gROOT.ProcessLine("kFullDotLarge"),
			      ('SH'):ROOT.gROOT.ProcessLine("kCircle"),
			      ('SF'):ROOT.gROOT.ProcessLine("kFullDotLarge"),
                  ('ER'):ROOT.gROOT.ProcessLine("kOpenStar"),
			      ('EC'):ROOT.gROOT.ProcessLine("kFullDotLarge"),
			      ('EH'):ROOT.gROOT.ProcessLine("kCircle"),
			      ('EF'):ROOT.gROOT.ProcessLine("kFullDotLarge")}

def nu(w,q2):
	return (w*w-MP*MP+q2)/(2*MP)

def epsilon(w,q2):
	n=nu(w,q2)
	e0=E1F_E0
	e1=e0-n
	epsInv=1+2*(q2+n*n)/(4*e0*e1-q2)
	return 1.0/epsInv

def getvgflux(w,q2,e0=E1F_E0):
	eps=epsilon(w,q2)
	return A*w*(w*w-MP*MP)/(4*PI*e0*e0*MP*MP*q2*(1-eps))

class DispObs:
	def __init__(self,obsdir,simnum='siml',norm=False,q2min=1.25,q2max=5.25,wmin=1.400,wmax=2.125,dbg=False):
		self.SIM_NUM=simnum
		self.NORM=norm
		self.Q2MIN,self.Q2MAX=q2min,q2max
		self.WMIN,self.WMAX=wmin,wmax
		self.DBG=dbg
				
		self.FIN=root_open(os.path.join(obsdir,self.SIM_NUM,'yield.root'))
		self.OUTDIR=os.path.join(obsdir,self.SIM_NUM)
		print "obsdir=%s\nFIN=%s\nOUTDIR=%s"%(obsdir,self.FIN,self.OUTDIR)

	def plot_1D(self,h1,q2wbin):
		"""
		"""
		#! Set up plotting aesthetics
		self.plot_1D_athtcs()
		self.hist_1D_athtcs(h1)
				
		#! TCanvas's pad_map[pad,vst,var] defined as per Gleb's display
		pad_map=[(1,1,"M1"),   (2,3,'M2'),   (3,2,'M2'),
		         (4,1,"THETA"),(5,3,'THETA'),(6,2,'THETA'),
		         (7,1,"ALPHA"),(8,3,'ALPHA'),(9,2,'ALPHA')]
				          
		c=ROOT.TCanvas("c","c",1000,1000)
		pad_t=ROOT.TPad("pad_t","Title pad",0.25,0.935,0.75,1.00)
		pad_p=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
  		if (self.NORM!=True):pad_p.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
  		pad_p.Draw()
  		pad_t.Draw()
  		pad_t.cd()
  		pt=ROOT.TPaveText(.05,.1,.95,.8)
  		pt.AddText("Q2_W bin=%s"%q2wbin)
  		pt.SetTextSize(0.40)
  		pt.Draw()
  		pad_p.Divide(3,3)
		for item in pad_map:
			pad,vst,var=item[0],item[1],item[2]
			print "pad,vst,var=",pad,vst,var
			gpad=pad_p.cd(pad)
			if self.NORM!=True:#! Draw all hists normed to same integral (since I want to compare their distributions)
				#! seq that can be directly drawn normalized
				seq_drct=['ST','SR','SF','ER','EF']
				#! seq that need to be manipulated
				#! + In this case the following seqs need to be scaled as per 
				#!   their counterpart h1n['F'] histograms such that after normalization
				#!   the following relationship is still preserved:
				#!   h1n[F]=h1n[C]+h1n[H]
				seq_mnpl_sim=['SC','SH']
				seq_mnpl_exp=['EC','EH']
				h1n=OrderedDict() 

				#! First draw 'seq_drct'
				for i,seq in enumerate(seq_drct):
					draw_opt=""
					if i>0:
						draw_opt="sames"
					h1n[seq]=h1[seq,vst,var].DrawNormalized(draw_opt,1000)
				#! Set the minimum and maximum of y coordinate of normed histograms
				maxl=[h1n[seq].GetMaximum() for seq in seq_drct]#[hEF_n.GetMaximum(),hSF_n.GetMaximum(),hER_n.GetMaximum(),hSR_n.GetMaximum(),hST_n.GetMaximum()]
				maximum=max(maxl)
				for htmp in [h1n[seq] for seq in seq_drct]:#[hEF_n,hSF_n,hER_n,hSR_n,hST_n]:
					htmp.SetMinimum(0.)
					htmp.SetMaximum(maximum+10)

				#! Now process and then draw seq_mnpl_sim/exp
				#! Scale EC, EH as per EF
				#! Obtain scale factor
				h_scl_exp=h1n['EF'].Clone("h_scl_exp")
				h_scl_exp.Divide(h1['EF',vst,var])
				#! Scale and draw
				for seq in seq_mnpl_exp:
					h1n[seq]=h1[seq,vst,var].Clone()
					h1n[seq].Multiply(h_scl_exp)
					h1n[seq].Draw("sames")
				#! Scale SC, SH as per SF
				#! Obtain scale factor
				h_scl_sim=h1n['SF'].Clone("h_scl_sim")
				h_scl_sim.Divide(h1['SF',vst,var])
				#! Scale and draw
				for seq in seq_mnpl_sim:
					h1n[seq]=h1[seq,vst,var].Clone()
					h1n[seq].Multiply(h_scl_sim)
					h1n[seq].Draw("sames")

				#! Re-draw h1n['ST'] for verification
				#! + This makes h1n['ST'] appear as a green '+' on top of SF;
				#!   othewise it is hidden under h1n['SF'] which is drawn earlier
				h1n['ST'].Draw("same")
				
				gpad.Update()

				#! Add TLegend if pad==1
				if pad==1:
					l=ROOT.TLegend(0.75,0.50,0.90,0.90)
					l.SetFillStyle(0)
					l.SetBorderSize(0)
					l.SetTextSize(0.06)
					for seq in seq_drct+seq_mnpl_sim+seq_mnpl_exp:
						l.AddEntry(h1n[seq],seq,"p")#EF
					l.Draw()
		#! Save canvas by q2bin (create folder tagged by q2bin)
		q2bin=self.get_q2bin(q2wbin)
		wbin=self.get_wbin(q2wbin)
		outdir_q2bin=os.path.join(self.OUTDIR_1D,"q%s"%q2bin)
		if not os.path.exists(outdir_q2bin):
			os.makedirs(outdir_q2bin)
		c.SaveAs("%s/c_w%s_q%s.png"%(outdir_q2bin,wbin,q2bin))
		c.SaveAs("%s/c_w%s_q%s.eps"%(outdir_q2bin,wbin,q2bin))
		return

	def extract_and_plot_R2(self,h5d,R2l,mthd,viewl,plotphiproj,dtypl,seql):#,hel,dtypl,seql):
		"""
		+ Given h5d, extract R2s using the specified method, plots R2s in specied view
		"""
		print "In DispObs::extract_obs_R2()"
		
		clrd={('EXP','C'):ROOT.gROOT.ProcessLine("kCyan"),('EXP','F'):ROOT.gROOT.ProcessLine("kBlue"),('SIM','F'):ROOT.gROOT.ProcessLine("kRed")}
		mrkr_styl=ROOT.gROOT.ProcessLine("kFullCircle")
		
		print "Going to extract and plot R2s=%s for mthd %s..."%(R2l,mthd)
		if mthd=='mthd1':#'h5-mply-itg'
			#! h5d=>hR2d, directly
			hR2d={'UNP':{},'POS':{},'NEG':{},'ASM':{}}
			for R2 in R2l:
				print "Doing h5d=>hR2d for mthd %s:R2=%s"%(mthd,R2)
				for hel in h5d:
					for k in h5d[hel]:
						q2bin_le,wbin_le,vst,dtyp,seq=k[0],k[1],k[2],k[3],k[4]
						h5=h5d[hel][k]
						if R2!='D':
							h5m=thntool.MultiplyBy(h5,self.H5_MFUNCD[R2],1)	
						elif R2=='D':
							if   hel=='POS' or hel=='UNP': h5m=thntool.MultiplyBy(h5,self.H5_MFUNCD[R2],1)
							elif hel=='NEG':               h5m=thntool.MultiplyBy(h5,self.H5_MFUNCD[R2],-1)	
												
						for var in VARS:
							if var=='PHI': continue
							if vst==1 and var=='M2': continue
							if vst==2 and var=='M1': continue
							if vst==3 and var=='M1': continue
							hR2d[hel][q2bin_le,wbin_le,vst,var,dtyp,seq]=h5m.Projection(H5_DIM[var],"E")
							hR2d[hel][q2bin_le,wbin_le,vst,var,dtyp,seq].SetName("%s_VST%d_%s_%s_%s"%(R2,vst,var,seq,hel))
							hR2d[hel][q2bin_le,wbin_le,vst,var,dtyp,seq].SetTitle("%s(%s:%s:%.2f,%.3f)"%(self.R2_NAMED[R2],self.VAR_NAMES[(vst,var)],hel,q2bin_le,wbin_le))
							hR2d[hel][q2bin_le,wbin_le,vst,var,dtyp,seq].Scale(1/math.pi)
							hR2d[hel][q2bin_le,wbin_le,vst,var,dtyp,seq].Scale(1/50000)
				print "Done h5d=>hR2d for mthd %s:R2=%s"%(mthd,R2)

				for view in viewl:
					print "Going to plot hR2d for mthd %s:R2=%s:view=%s"%(mthd,R2,view)
					if   view=="view1":self.plot_obs_R2(hR2d,R2,dtypl,seql)
					elif view=="view2":self.plot_obs_R2_view2(hR2d,R2,dtypl,seql)
					print "Done plot hR2d for mthd %s:R2=%s:view=%s"%(mthd,R2,view)	
		elif mthd=='mthd2' or mthd=='mthd3':#proj-phi-fit or proj-phi-mply-itg
			#! 1. h5d=>hphiprojd
			print "Doing h5d=>hphiprojd for method %s..."%mthd
			hphiprojd={'UNP':{},'POS':{},'NEG':{},'ASM':{}}
			
			for hel in h5d:
				for k in h5d[hel]:
					q2bin_le,wbin_le,vst,dtyp,seq=k[0],k[1],k[2],k[3],k[4]
					h5=h5d[hel][k]
					for var in VARS:
						if var=='PHI': continue
						if vst==1 and var=='M2': continue
						if vst==2 and var=='M1': continue
						if vst==3 and var=='M1': continue
						nbins=self.NBINS[var]
						#! Now make projections on to PHI per bin
						for ibin in range(nbins):
							h5.GetAxis(H5_DIM[var]).SetRange(ibin+1,ibin+1)
							h=h5.Projection(H5_DIM['PHI'],"E")
							#! Modify the title so also state projection range
							#! + original title(set in proc_yields.py)=[q2min-q2max]_[wmin-wmax]_VSTX_SEQ_HEL=(POS/NEG/UNPOL)
							#! Get q2wbintitle
							# orig_ttl=h5.GetTitle().split("_")
							# q2wbintitle="%s_%s"%(orig_ttl[0],orig_ttl[1])
							#! get projection bin edges
							xmin=h5.GetAxis(H5_DIM[var]).GetBinLowEdge(ibin+1)
							xmax=h5.GetAxis(H5_DIM[var]).GetBinUpEdge(ibin+1)
							new_ttl="%s:proj. bin=[%.3f,%.3f]"%(self.VAR_NAMES[(vst,var)],xmin,xmax)
							h.SetTitle(new_ttl)
							h.SetMarkerStyle(mrkr_styl)
							h.SetMarkerColor(clrd[(dtyp,seq)])
							#! Finally, add this histogram to hphiprojd
							hphiprojd[hel][q2bin_le,wbin_le,vst,var,ibin+1,dtyp,seq]=h
							#! Before leaving projection loop, reset projection range for h5!
							h5.GetAxis(H5_DIM[var]).SetRange()
			print "Done h5d=>hphiprojd for method %s..."%mthd
			# for hel in hphiprojd:				
			# 	print "keys in hphiprojd[%s]"%hel
			# 	print hphiprojd[hel].keys()

			#! 2. hphiprojd=>fpard
			print "Doing hphiprojd=>fpard"
			fpard= {'UNP':{},'POS':{},'NEG':{},'ASM':{}}
			#fsd={'UNP':{},'POS':{},'NEG':{},'ASM':{}} #for storing fit stat
			for hel in hphiprojd:
				for k in hphiprojd[hel]:
					q2bin_le,wbin_le,vst,var,bin,dtyp,seq=k[0],k[1],k[2],k[3],k[4],k[5],k[6]
					f=self.FPHI
					print "Going to fit phi proj for",hel,k
					fstat=int(hphiprojd[hel][k].Fit(f,"NQ"))
					print "fstat=",fstat
					if fstat==0:
						A=   f.GetParameter(0)
						Aerr=f.GetParError(0)
						B=   f.GetParameter(1)
						Berr=f.GetParError(1)
						C=   f.GetParameter(2)
						Cerr=f.GetParError(2)
						D=   f.GetParameter(3)
						Derr=f.GetParError(3)
						fpard[hel][k]={'A':A,'Aerr':Aerr,'B':B,'Berr':Berr,'C':C,'Cerr':Cerr,'D':D,'Derr':Derr}
					else:
						fpard[hel][k]={'A':0,'Aerr':0,'B':0,'Berr':0,'C':0,'Cerr':0,'D':0,'Derr':0}
			print "Done hphiprojd=>fpard"
			# for hel in fpard:
			# 	for k in fpard[hel]:
			# 		print "dbg: fit pars for fpard[%s][%s]=%f,%f,%f,%f"%(hel,k,fpard[hel][k]['A'],fpard[hel][k]['B'],
			# 			fpard[hel][k]['C'],fpard[hel][k]['D'])

			#! 2. Obtain hR2 for each R2 from fpard and plot
			for R2 in R2l:
				print "Going to extract R2 from phiproj fits for mthd %s:R2=%s"%(mthd,R2)
				hR2d=self.extract_R2_from_phiproj(h5d,fpard,R2)
				print "Done to extract R2 from phiproj fits for mthd %s:R2=%s"%(mthd,R2)

				for view in viewl:
					print "Going to plot hR2d for mthd %s:R2=%s:view=%s"%(mthd,R2,view)
					if   view=="view1":self.plot_obs_R2(hR2d,R2,dtypl,seql)
					elif view=="view2":self.plot_obs_R2_view2(hR2d,R2,dtypl,seql)
					print "Done plot hR2d for mthd %s:R2=%s:view=%s"%(mthd,R2,view)

			#! 3. Finally, for visual verification, plot phiprojs
			if plotphiproj:
				print "Going to plot phi-proj and extract R2 for method %s"%mthd
				self.plot_phiproj(hphiprojd,fpard,dtypl,seql)
				print "Done to plot phi-proj and extract R2 for method %s"%mthd
			
			print "Done h5d=>hR2d for method %s..."%mthd
		print "Done to extract and plot R2s=%s for mthd %s..."%(R2l,mthd)		
			

	def plot_obs_R2(self,hR2d,R2,dtypl,seql):
		print "In DispObs::plot_obs_R2()"
		
		#! Set up ROOT-display aesthetics
		self.plot_obs_R2_athtcs()

		#! Now display hR2 
		#! Get all q2- and w-bins from the keys of hR2['UNPOL']
		#! (Code that actually makes the plots is not yet Function-alized and therefore
		#! I have to go through the tedious process of Loop-ing, which requires me to get	
		#! the q2- and w-bin informtion from h5d)
		q2bins_lel,wbins_lel=self.get_q2bin_lel_wbin_lel(hR2d['UNP'])
		print "going to plot R2 obs for:"
		print "Q2:"
		print q2bins_lel
		print "W:"
		print wbins_lel

		#! Russian-normalize theta distributions
		for hel in hR2d:
			for q2bin in q2bins_lel:
				for wbin in wbins_lel:
					for vst in self.VSTS:
						for dtyp in dtypl:
								for seq in seql:
									if hR2d[hel].has_key((q2bin,wbin,vst,'THETA',dtyp,seq)):
										self.russ_norm_theta_dist(hR2d[hel][q2bin,wbin,vst,'THETA',dtyp,seq])

		clrd={('EXP','C'):ROOT.gROOT.ProcessLine("kCyan"),('EXP','F'):ROOT.gROOT.ProcessLine("kBlue"),('SIM','F'):ROOT.gROOT.ProcessLine("kRed")}
		for vst in self.VSTS:
			for var in VARS:
				if var=='PHI': continue
				if vst==1 and var=='M2': continue
				if vst==2 and var=='M1': continue
				if vst==3 and var=='M1': continue

				for hel in hR2d.keys():
					outdir=os.path.join(self.OUTDIR_OBS_R2,"view1",R2,"VST%d_%s"%(vst,var),"%s"%hel)
					if not os.path.exists(outdir):
						os.makedirs(outdir)
					for q2bin_le in q2bins_lel:
						c=ROOT.TCanvas("c","c",6000,6000)
						pad_t=ROOT.TPad("pad_t","Title pad",0.05,0.97,0.95,1.00)
						#pad_t.SetFillColor(11)
  						pad_p=ROOT.TPad("pad_p","Plots pad",0.01,0.01,0.99,0.95);
  						pad_p.SetFillColor(11)
  						pad_t.Draw()
  						pad_p.Draw()
  						pad_t.cd()
  						pt=ROOT.TPaveText(.05,.1,.95,.8)
  						pt.AddText("%s(%s:hel=%s:Q2=%.2f)"%(self.R2_NAMED[R2],self.VAR_NAMES[(vst,var)],hel,q2bin_le))
  						pt.Draw()
 						pad_p.cd()
 						pad_p.Divide(7,10) #! Make this dynamic depending on W limits
						#pad_p.Divide(7,10)
						#c.Divide(7,10)
						ipad=0
						l,t,ln=[],[],[]
						for wbin_le in wbins_lel:
							pad=pad_p.cd(ipad+1)
							#! Legend for this pad
							l.append(ROOT.TLegend(0.1,0.8,0.2,0.9))
							#! Get all hists to be drawn on this pad
							hl=[]
							i=0 #index for histograms in hl
							for dtyp in dtypl:
								for seq in seql:
									if dtyp=='SIM' and seq=='C': continue
									if hR2d[hel].has_key((q2bin_le,wbin_le,vst,var,dtyp,seq)):
										hl.append(hR2d[hel][q2bin_le,wbin_le,vst,var,dtyp,seq])
										hl[i].SetMarkerStyle(ROOT.gROOT.ProcessLine("kFullCircle"))
										hl[i].SetMarkerColor(clrd[(dtyp,seq)])
										hl[i].SetTitle("")
										hl[i].SetXTitle( "%s%s"%(self.VAR_NAMES[(vst,var)],self.VAR_UNIT_NAMES[var]) )
										hl[i].GetXaxis().SetLabelSize(.05)
										hl[i].GetXaxis().SetTitleSize(.10)
										hl[i].GetXaxis().SetTitleOffset(.7)
										l[ipad].AddEntry(hl[i],"%s-%s"%(dtyp,seq),"p")
										i+=1
							#! Normalized histograms to the signed integral of the 1st hist in hl 
							itg=self.get_signed_integral(hl[0])
							#! Scale the rest of the hists to itg 
							#! and
							#! make list of min and max of each hist (used later to set limit)
							minl=[]
							maxl=[]
							for i in range(len(hl)):
								#! First call Sumw2() so that errors are also scaled
								hl[i].Sumw2()
								#! Now get scale factor and scale bin contents
								scale_factor=itg/self.get_signed_integral(hl[i])
								hl[i].Scale(scale_factor)
								#! minl, maxl
								minl.append(hl[i].GetMinimum())
								maxl.append(hl[i].GetMaximum())
							#! Get minimum and maximum from minl and maxl
							minimum=min(minl)
							maximum=max(maxl)
							#! Now draw histograms on the same pad
							fctr=1.5
							for i in range(len(hl)):
								if i==0:
									if R2!="A":
										# hl[i].SetMinimum(minimum-math.fabs((50/100)*minimum))
										# hl[i].SetMaximum(maximum+math.fabs((50/100)*maximum))
										hl[i].SetMinimum(minimum-math.fabs(fctr*minimum))
										hl[i].SetMaximum(maximum+math.fabs(fctr*maximum))
									hl[i].Draw()
									#! TLine at y=0
									ln.append(ROOT.TLine(hl[i].GetXaxis().GetXmin(),0,hl[i].GetXaxis().GetXmax(),0))
									#! Draw TLine only if ymin<0
									if hl[i].GetMinimum()<0:
										ln[ipad].Draw("same")
								else:
									hl[i].Draw("sames")
							#! Draw legend
							l[ipad].Draw()
							#! Add TText on pad to label W bin
							#pad.cd()
							t.append(ROOT.TText(0.5,0.5,"W=%s"%wbin_le))
							#print "W label=",t[ipad].GetTitle()
							t[ipad].SetNDC()
							t[ipad].SetTextSize(0.1)
							t[ipad].SetTextColor(ROOT.gROOT.ProcessLine("kMagenta"))
							t[ipad].Draw()
							#pad.Modified()
							#! Increment ipad for next W-bin
							ipad+=1

							#! [11-15-14] Before drawing Exp. and Sim on same pad
							# if hel=='UNP' and len(dtypl)==2: 
							# 	pad.Divide(1,2)
							# l.append(ROOT.TLegend(0.1,0.8,0.2,0.9))
							# i=0
							# for dtyp in dtypl:
							# 	for seq in seql:
							# 		if dtyp=='SIM' and seq=='C': continue
							# 		if hR2d[hel].has_key((q2bin_le,wbin_le,vst,var,dtyp,seq)):
							# 			#print "hR2 key=",q2bin_le,wbin_le,hel,vst,var,dtyp,seq,"exists"
							# 			h=hR2d[hel][q2bin_le,wbin_le,vst,var,dtyp,seq]
							# 			h.SetMarkerStyle(ROOT.gROOT.ProcessLine("kFullCircle"))
							# 			h.SetMarkerColor(clrd[(dtyp,seq)])
							# 			h.SetXTitle( "%s%s"%(self.VAR_NAMES[(vst,var)],self.VAR_UNIT_NAMES[var]) )
							# 			h.GetXaxis().SetLabelSize(.05)
							# 			h.GetXaxis().SetTitleSize(.10)
							# 			h.GetXaxis().SetTitleOffset(.7)
							# 			l[ipad].AddEntry(h,"%s-%s"%(dtyp,seq),"p")
							# 			if dtyp=='EXP': 
							# 				pad.cd(1)
							# 				if i==0:
							# 					h.Draw()
							# 					i+=1;
							# 				else:
							# 					h.Draw("sames")
							# 			if dtyp=='SIM':
							# 				pad.cd(2)
							# 				h.Draw()
							# 		#else:
							# 			#print "hR2 key=",q2bin_le,wbin_le,hel,vst,var,dtyp,seq,"does NOT exist"
							# l[ipad].Draw()
							# #! Add TText on pad to label W bin
							# #pad.cd()
							# t.append(ROOT.TText(0.5,0.5,"W=%s"%wbin_le))
							# #print "W label=",t[ipad].GetTitle()
							# t[ipad].SetNDC()
							# t[ipad].SetTextSize(0.1)
							# t[ipad].SetTextColor(ROOT.gROOT.ProcessLine("kMagenta"))
							# t[ipad].Draw()
							# #pad.Modified()
							# #! Increment ipad for next W-bin
							# ipad+=1

						c.SaveAs("%s/c_q%0.2f.png"%(outdir,q2bin_le))
						c.SaveAs("%s/c_q%0.2f.eps"%(outdir,q2bin_le))
						c.Close()
		print "Done DispObs::plot_obs_R2()"

	def plot_obs_R2_view2(self,hR2d,R2,dtypl,seql):
		"""
		+ Plot R2s just like 1D observables are plotted
		"""

		print "In DispObs::plot_obs_R2_view2()"

		#! Get all q2- and w-bins from the keys of hR2d['UNP']
		q2bins_le,wbins_le=self.get_q2bin_lel_wbin_lel(hR2d['UNP'])
		print "going to plot 1D obs for:"
		print "Q2:"
		print q2bins_le
		print "W:"
		print wbins_le

		#! Russian-normalize theta distributions
		#! [11-16-14] This is commented out for view2, since view1 is
		#! called before view2 and in the same process. Therefore, the 
		#! following Russian-normalizes the distribution again!
		# for hel in hR2d:
		# 	for q2bin in q2bins_le:
		# 		for wbin in wbins_le:
		# 			for vst in self.VSTS:
		# 				for dtyp in dtypl:
		# 						for seq in seql:
		# 							if hR2d[hel].has_key((q2bin,wbin,vst,'THETA',dtyp,seq)):
		# 								self.russ_norm_theta_dist(hR2d[hel][q2bin,wbin,vst,'THETA',dtyp,seq])

		
		#! Set up some plotting related styles and aesthetics 
		self.plot_obs_R2_view2_athtcs()
		clrd={('EXP','C'):ROOT.gROOT.ProcessLine("kCyan"),('EXP','F'):ROOT.gROOT.ProcessLine("kBlue"),('SIM','F'):ROOT.gROOT.ProcessLine("kRed")}
				
		#ivars=[0,1,2]
		ivars=[[0,2,4],[1,2,4],[1,2,4]]
		#! 3x3 TCanvas-pads as per VSTs: VST1=1,2,3; VST2=3,6,9; VST3=2,5,8 (To adapt to Gleb's display)
		pads=[[1,4,7],[3,6,9],[2,5,8]]
		#map_padnums_vars=[zip(pads[0],ivars),zip(pads[1],ivars),zip(pads[2],ivars)]
		map_padnums_vars=[zip(pads[0],ivars[0]),zip(pads[1],ivars[1]),zip(pads[2],ivars[2])]
		for wbin in wbins_le:
			for iq2bin,q2bin in enumerate(q2bins_le):
				for hel in hR2d:
					print "Plotting hR2 for R2=%s,w=%0.3f,q2=%0.2f,hel=%s"%(R2,wbin,q2bin,hel)
					c=ROOT.TCanvas("c","c",1000,1000)
					#pad_t=ROOT.TPad("pad_t","Title pad",0.05,0.97,0.95,1.00)
					pad_t=ROOT.TPad("pad_t","Title pad",0.15,0.945,0.85,1.00)
					#pad_t.SetFillColor(11)
  					pad_p=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
  					#pad_p.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
  					pad_p.Draw()
  					pad_t.Draw()
  					pad_t.cd()
  					pt=ROOT.TPaveText(.05,.1,.95,.8)
  					pt.AddText("%s for Q2,W =(%s,%s):hel=%s"%(self.R2_NAMED[R2],q2bin,wbin,hel))
  					pt.SetTextSize(0.42)
  					pt.Draw()
  					#pad_p.Draw()
  					#pad_p.cd()
					pad_p.Divide(3,3)
					ipad=0
					l,t,ln=[],[],[]
					for ivst,vst in enumerate(self.VSTS):
						for m in map_padnums_vars[ivst]:
							padnum=m[0]
							ivar=m[1]
							var=VARS[ivar]
							pad=pad_p.cd(padnum)
							#! Legend for this pad
							#l.append(ROOT.TLegend(0.1,0.8,0.2,0.9))
							l.append(ROOT.TLegend(0.70,0.75,0.90,0.90))
							l[ipad].SetFillStyle(0)
							l[ipad].SetTextSize(0.05)
							#! Get all hists to be drawn on this pad
							hl=[]
							i=0 #index for histograms in hl
							for dtyp in dtypl:
								for seq in seql:
									if dtyp=='SIM' and seq=='C': continue
									if hR2d[hel].has_key((q2bin,wbin,vst,var,dtyp,seq)):
										hl.append(hR2d[hel][q2bin,wbin,vst,var,dtyp,seq])
										hl[i].SetMarkerStyle(ROOT.gROOT.ProcessLine("kFullCircle"))
										hl[i].SetMarkerColor(clrd[(dtyp,seq)])
										hl[i].SetTitle("")
										hl[i].SetXTitle( "%s%s"%(self.VAR_NAMES[(vst,var)],self.VAR_UNIT_NAMES[var]) )
										hl[i].GetXaxis().SetLabelSize(.05)
										hl[i].GetXaxis().SetTitleSize(.10)
										hl[i].GetXaxis().SetTitleOffset(.7)
										l[ipad].AddEntry(hl[i],"%s-%s"%(dtyp,seq),"p")
										i+=1
							#! Normalized histograms to the signed integral of the 1st hist in hl 
							itg=self.get_signed_integral(hl[0])
							#! Scale the rest of the hists to itg 
							#! and
							#! make list of min and max of each hist (used later to set limit)
							minl=[]
							maxl=[]
							for i in range(len(hl)):
								#! First call Sumw2() so that errors are also scaled
								hl[i].Sumw2()
								#! Now get scale factor and scale bin contents
								scale_factor=itg/self.get_signed_integral(hl[i])
								hl[i].Scale(scale_factor)
								#! minl, maxl
								minl.append(hl[i].GetMinimum())
								maxl.append(hl[i].GetMaximum())
							#! Get minimum and maximum from minl and maxl
							minimum=min(minl)
							maximum=max(maxl)
							#! Now draw histograms on the same pad
							fctr=1.5
							for i in range(len(hl)):
								if i==0:
									if R2!="A":
										# hl[i].SetMinimum(minimum-math.fabs((50/100)*minimum))
										# hl[i].SetMaximum(maximum+math.fabs((50/100)*maximum))
										hl[i].SetMinimum(minimum-math.fabs(fctr*minimum))
										hl[i].SetMaximum(maximum+math.fabs(fctr*maximum))
									hl[i].Draw()
									#! TLine at y=0
									ln.append(ROOT.TLine(hl[i].GetXaxis().GetXmin(),0,hl[i].GetXaxis().GetXmax(),0))
									#! Draw TLine only if ymin<0
									if hl[i].GetMinimum()<0:
										ln[ipad].Draw("same")
								else:
									hl[i].Draw("sames")
							#! Draw legend
							if padnum==1:
								l[ipad].Draw()
							ipad+=1

							#! [11-15-14] Before drawing Exp. and Sim on same pad	
							# l.append(ROOT.TLegend(0.1,0.8,0.2,0.9))
							# if hel=='UNP' and len(dtypl)==2: 
							# 	pad.Divide(1,2)
							# l.append(ROOT.TLegend(0.1,0.8,0.2,0.9))
							# i=0
							# for dtyp in dtypl:
							# 	for seq in seql:
							# 		if dtyp=='SIM' and seq=='C': continue
							# 		if hR2d[hel].has_key((q2bin,wbin,vst,var,dtyp,seq)):
							# 			h=hR2d[hel][q2bin,wbin,vst,var,dtyp,seq]
							# 			h.SetMarkerStyle(ROOT.gROOT.ProcessLine("kFullCircle"))
							# 			h.SetMarkerColor(clrd[(dtyp,seq)])
							# 			h.SetTitle("")
							# 			h.SetXTitle( "%s%s"%(self.VAR_NAMES[(vst,var)],self.VAR_UNIT_NAMES[var]) )
							# 			h.GetXaxis().SetLabelSize(.05)
							# 			h.GetXaxis().SetTitleSize(.10)
							# 			h.GetXaxis().SetTitleOffset(.7)
							# 			l[ipad].AddEntry(h,"%s-%s"%(dtyp,seq),"p")
							# 			if dtyp=='EXP': 
							# 				pad.cd(1)
							# 				if i==0:
							# 					h.Draw()
							# 					i+=1;
							# 				else:
							# 					h.Draw("sames")
							# 			if dtyp=='SIM':
							# 				pad.cd(2)
							# 				h.Draw()
							# l[ipad].Draw()
							# ipad+=1
					#outdir=os.path.join(self.OUTDIR_OBS_R2,"view2","q%.2f"%q2bin,hel)
					outdir=os.path.join(self.OUTDIR_OBS_R2,"view2","q%.2f_w%.3f"%(q2bin,wbin),hel)
					if not os.path.exists(outdir):
						os.makedirs(outdir)
					# c.SaveAs("%s/c%s_w%.3f_q%.2f.png"%(outdir,R2,wbin,q2bin))
					# c.SaveAs("%s/c%s_w%.3f_q%.2f.eps"%(outdir,R2,wbin,q2bin))
					c.SaveAs("%s/c%s.png"%(outdir,R2))
					c.SaveAs("%s/c%s.eps"%(outdir,R2))
					c.Close()
		print "Done DispObs::plot_obs_R2_view2()"
		return

	def extract_R2_from_phiproj(self,h5d,fpard,R2):
		"""
		+ From fpard, for a particular R2, extract hR2
		"""
		print "In DispObs::extract_R2_from_phiproj for R2=%s"%(R2)

		clrd={('EXP','C'):ROOT.gROOT.ProcessLine("kCyan"),('EXP','F'):ROOT.gROOT.ProcessLine("kBlue"),('SIM','F'):ROOT.gROOT.ProcessLine("kRed")}

		hR2d={'UNP':{},'POS':{},'NEG':{},'ASM':{}}
		#! the following is a very "hack-ish" way to create hR2d by
		#1 basically doing h5d=>hR2d directly as per mthd1 and then reseting the contents of the hR2!
		for hel in h5d:
			for k in h5d[hel]:
				q2bin_le,wbin_le,vst,dtyp,seq=k[0],k[1],k[2],k[3],k[4]
				h5=h5d[hel][k]
				# if R2!='D':
				# 	h5m=thntool.MultiplyBy(h5,self.H5_MFUNCD[R2],1)	
				# elif R2=='D':
				# 	if   hel=='POS' or hel=='UNP': h5m=thntool.MultiplyBy(h5,self.H5_MFUNCD[R2],1)
				# 	elif hel=='NEG':               h5m=thntool.MultiplyBy(h5,self.H5_MFUNCD[R2],-1)	
											
				for var in VARS:
					if var=='PHI': continue
					if vst==1 and var=='M2': continue
					if vst==2 and var=='M1': continue
					if vst==3 and var=='M1': continue
					hR2d[hel][q2bin_le,wbin_le,vst,var,dtyp,seq]=h5.Projection(H5_DIM[var],"E")
					hR2d[hel][q2bin_le,wbin_le,vst,var,dtyp,seq].SetName("%s_VST%d_%s_%s_%s"%(R2,vst,var,seq,hel))
					hR2d[hel][q2bin_le,wbin_le,vst,var,dtyp,seq].SetTitle("%s(%s:%s:%.2f,%.3f)"%(self.R2_NAMED[R2],self.VAR_NAMES[(vst,var)],hel,q2bin_le,wbin_le))
					hR2d[hel][q2bin_le,wbin_le,vst,var,dtyp,seq].Reset()
					
		print "Going to extract R2 from phi-proj fits"
		for hel in fpard:
			for k in fpard[hel]:
				q2bin_le,wbin_le,vst,var,bin,dtyp,seq=k[0],k[1],k[2],k[3],k[4],k[5],k[6]
				#fpard[hel][q2bin_le,wbin_le,vst,var,bin,dtyp,seq]=f
				#! Now fill hR2 from fit
				if R2=='A':
					r2=    fpard[hel][k]['A']
					r2_err=fpard[hel][k]['Aerr']
				elif R2=='B':
					r2=    fpard[hel][k]['B']
					r2_err=fpard[hel][k]['Berr']
				elif R2=='C':
					r2=    fpard[hel][k]['C']
					r2_err=fpard[hel][k]['Cerr']
				elif R2=='D':
					if   hel=='POS' or hel=='UNP': r2=   fpard[hel][k]['D']
					elif hel=='NEG':               r2=-1*fpard[hel][k]['D']
					r2_err=fpard[hel][k]['Derr']
				#! Hack-ish way to get norm to match with mthd1
				if R2=='A': nmlzn_fctr=2/50000
				else:		nmlzn_fctr=1/50000
				hR2d[hel][q2bin_le,wbin_le,vst,var,dtyp,seq].SetBinContent(bin,r2*nmlzn_fctr)
				hR2d[hel][q2bin_le,wbin_le,vst,var,dtyp,seq].SetBinError(bin,r2_err*nmlzn_fctr)
		print "Done extracting R2 from phi-proj fits"
		return hR2d

	def plot_phiproj(self,hphiprojd,fpard,dtypl,seql):
		"""
		+ plot phiprojs
		"""
		print "In DispObs::plot_phiproj"

		#! Set up ROOT-display aesthetics
		self.plot_phiproj_athtcs()

		#! Now display phi-proj 
		#! Get all q2- and w-bins from the keys of hphiprojd['UNPOL']
		#! (Code that actually makes the plots is not yet Function-alized and therefore
		#! I have to go through the tedious process of Loop-ing, which requires me to get	
		#! the q2- and w-bin informtion from h5d)
		q2bins_lel,wbins_lel=self.get_q2bin_lel_wbin_lel(hphiprojd['UNP'])
		print "going to plot phi-proj for:"
		print "Q2:"
		print q2bins_lel
		print "W:"
		print wbins_lel

		clrd={('EXP','C'):ROOT.gROOT.ProcessLine("kCyan"),('EXP','F'):ROOT.gROOT.ProcessLine("kBlue"),('SIM','F'):ROOT.gROOT.ProcessLine("kRed")}
		for vst in self.VSTS:
			for var in VARS:
				if var=='PHI': continue
				if vst==1 and var=='M2': continue
				if vst==2 and var=='M1': continue
				if vst==3 and var=='M1': continue

				for hel in hphiprojd.keys():
					for q2bin_le in q2bins_lel:
						for wbin_le in wbins_lel:

							outdir=os.path.join(self.OUTDIR_OBS_R2,"phiprojs","VST%d_%s"%(vst,var),"%s"%hel,"q%.2f_w%.3f"%(q2bin_le,wbin_le))
							if not os.path.exists(outdir):
								os.makedirs(outdir)

							c=ROOT.TCanvas("c","c",5000,7000)
							pad_t=ROOT.TPad("pad_t","Title pad",0.05,0.97,0.95,1.00)
							#pad_t.SetFillColor(11)
  							pad_p=ROOT.TPad("pad_p","Plots pad",0.01,0.01,0.99,0.95);
  							pad_p.SetFillColor(11)
  							pad_t.Draw()
  							pad_p.Draw()
  							pad_t.cd()
  							pt=ROOT.TPaveText(.05,.1,.95,.8)
  							pt.AddText("#phi-projections(%s:hel=%s:Q2,W=%.2f,%.3f)"%(self.VAR_NAMES[(vst,var)],hel,q2bin_le,wbin_le))
  							pt.Draw()
 							pad_p.cd()
 							nbins=self.NBINS[var]
 							pad_p.Divide(self.NXPADS[var],self.NYPADS[var])
 							l=[]
 							f_exp1,f_exp2,f_sim=[],[],[]# used to draw fit parameters
 							pt_exp1,pt_exp2,pt_sim=[],[],[]# used to draw fit parameters
 							for ibin in range(nbins):
								pad=pad_p.cd(ibin+1)
								if hel=='UNP' and len(dtypl)==2:
									pad.Divide(1,2)
								l.append(ROOT.TLegend(0.1,0.8,0.2,0.9))
								i=0
								for dtyp in dtypl:
									for seq in seql:
										#if dtyp=='SIM' and seq=='C': continue
										if hphiprojd[hel].has_key((q2bin_le,wbin_le,vst,var,ibin+1,dtyp,seq)):
											k=q2bin_le,wbin_le,vst,var,ibin+1,dtyp,seq
											#h=hphiprojd[hel][q2bin_le,wbin_le,vst,var,ibin+1,dtyp,seq]
											h=hphiprojd[hel][k]
											h.SetXTitle( "%s%s"%(self.VAR_NAMES[(vst,'PHI')],self.VAR_UNIT_NAMES['PHI']) )
											h.GetXaxis().SetLabelSize(.05)
											h.GetXaxis().SetTitleSize(.10)
											h.GetXaxis().SetTitleOffset(.7)
											l[ibin].AddEntry(h,"%s-%s"%(dtyp,seq),"p")
											if dtyp=='EXP': 
												pad.cd(1)
												if i==0:
													#! Draw hist
													h.Draw()
													#! Draw Fit
													f_exp1.append(ROOT.TF1("fphi", "([0] + [1]*cos(x*TMath::DegToRad()) + [2]*cos(2*x*TMath::DegToRad()) + [3]*sin(x*TMath::DegToRad()))",0,360))
													f_exp1[ibin].SetParameter(0,fpard[hel][k]['A'])
													f_exp1[ibin].SetParError (0,fpard[hel][k]['Aerr'])
  													f_exp1[ibin].SetParameter(1,fpard[hel][k]['B'])
  													f_exp1[ibin].SetParError (1,fpard[hel][k]['Berr'])
  													f_exp1[ibin].SetParameter(2,fpard[hel][k]['C'])
  													f_exp1[ibin].SetParError (2,fpard[hel][k]['Cerr'])
  													f_exp1[ibin].SetParameter(3,fpard[hel][k]['D'])
  													f_exp1[ibin].SetParError (3,fpard[hel][k]['Derr'])
 													f_exp1[ibin].SetParName(0, "A")
  													f_exp1[ibin].SetParName(1, "B")
  													f_exp1[ibin].SetParName(2, "C")
  													f_exp1[ibin].SetParName(3, "D")#hPD
  													f_exp1[ibin].SetLineColor(h.GetMarkerColor())
													f_exp1[ibin].Draw("sames")
													#! Draw Fit stats
													pt_exp1.append(ROOT.TPaveText(.20,.6,.30,1.0,"NDC"))
													pt_exp1[ibin].AddText( "A=%.2f+/-%.2f"%(f_exp1[ibin].GetParameter(0),f_exp1[ibin].GetParError(0)) )
													pt_exp1[ibin].AddText( "B=%.2f+/-%.2f"%(f_exp1[ibin].GetParameter(1),f_exp1[ibin].GetParError(1)) )
													pt_exp1[ibin].AddText( "C=%.2f+/-%.2f"%(f_exp1[ibin].GetParameter(2),f_exp1[ibin].GetParError(2)) )
													pt_exp1[ibin].AddText( "D=%.2f+/-%.2f"%(f_exp1[ibin].GetParameter(3),f_exp1[ibin].GetParError(3)) )
													pt_exp1[ibin].SetTextColor(h.GetMarkerColor())
													pt_exp1[ibin].SetFillColor(11)
													pt_exp1[ibin].Draw()
													i+=1;
												else:
													#! Draw hist
													h.Draw("sames")
													#! Draw Fit
													f_exp2.append(ROOT.TF1("fphi", "([0] + [1]*cos(x*TMath::DegToRad()) + [2]*cos(2*x*TMath::DegToRad()) + [3]*sin(x*TMath::DegToRad()))",0,360))
													f_exp2[ibin].SetParameter(0,fpard[hel][k]['A'])
													f_exp2[ibin].SetParError (0,fpard[hel][k]['Aerr'])
  													f_exp2[ibin].SetParameter(1,fpard[hel][k]['B'])
  													f_exp2[ibin].SetParError (1,fpard[hel][k]['Berr'])
  													f_exp2[ibin].SetParameter(2,fpard[hel][k]['C'])
  													f_exp2[ibin].SetParError (2,fpard[hel][k]['Cerr'])
  													f_exp2[ibin].SetParameter(3,fpard[hel][k]['D'])
  													f_exp2[ibin].SetParError (3,fpard[hel][k]['Derr'])
 													f_exp2[ibin].SetParName(0, "A")
  													f_exp2[ibin].SetParName(1, "B")
  													f_exp2[ibin].SetParName(2, "C")
  													f_exp2[ibin].SetParName(3, "D")#hPD
  													f_exp2[ibin].SetLineColor(h.GetMarkerColor())
													f_exp2[ibin].Draw("sames")
													#! Draw Fit stats
													pt_exp2.append(ROOT.TPaveText(.70,.6,.80,1.0,"NDC"))
													pt_exp2[ibin].AddText( "A=%.2f+/-%.2f"%(f_exp2[ibin].GetParameter(0),f_exp2[ibin].GetParError(0)) )
													pt_exp2[ibin].AddText( "B=%.2f+/-%.2f"%(f_exp2[ibin].GetParameter(1),f_exp2[ibin].GetParError(1)) )
													pt_exp2[ibin].AddText( "C=%.2f+/-%.2f"%(f_exp2[ibin].GetParameter(2),f_exp2[ibin].GetParError(2)) )
													pt_exp2[ibin].AddText( "D=%.2f+/-%.2f"%(f_exp2[ibin].GetParameter(3),f_exp2[ibin].GetParError(3)) )
													pt_exp2[ibin].SetTextColor(h.GetMarkerColor())
													pt_exp2[ibin].SetFillColor(11)
													pt_exp2[ibin].Draw()
											if dtyp=='SIM':
												pad.cd(2)
												#! Draw hist
												h.Draw()
												#! Draw Fit
												f_sim.append(ROOT.TF1("fphi", "([0] + [1]*cos(x*TMath::DegToRad()) + [2]*cos(2*x*TMath::DegToRad()) + [3]*sin(x*TMath::DegToRad()))",0,360))
												f_sim[ibin].SetParameter(0,fpard[hel][k]['A'])
												f_sim[ibin].SetParError (0,fpard[hel][k]['Aerr'])
  												f_sim[ibin].SetParameter(1,fpard[hel][k]['B'])
  												f_sim[ibin].SetParError (1,fpard[hel][k]['Berr'])
  												f_sim[ibin].SetParameter(2,fpard[hel][k]['C'])
  												f_sim[ibin].SetParError (2,fpard[hel][k]['Cerr'])
  												f_sim[ibin].SetParameter(3,fpard[hel][k]['D'])
  												f_sim[ibin].SetParError (3,fpard[hel][k]['Derr'])
 												f_sim[ibin].SetParName(0, "A")
  												f_sim[ibin].SetParName(1, "B")
  												f_sim[ibin].SetParName(2, "C")
  												f_sim[ibin].SetParName(3, "D")#hPD
  												f_sim[ibin].SetLineColor(h.GetMarkerColor())
												f_sim[ibin].Draw("sames")
												#! Draw Fit stats
												pt_sim.append(ROOT.TPaveText(.20,.6,.30,1.0,"NDC"))
												pt_sim[ibin].AddText( "A=%.2f+/-%.2f"%(f_sim[ibin].GetParameter(0),f_sim[ibin].GetParError(0)) )
												pt_sim[ibin].AddText( "B=%.2f+/-%.2f"%(f_sim[ibin].GetParameter(1),f_sim[ibin].GetParError(1)) )
												pt_sim[ibin].AddText( "C=%.2f+/-%.2f"%(f_sim[ibin].GetParameter(2),f_sim[ibin].GetParError(2)) )
												pt_sim[ibin].AddText( "D=%.2f+/-%.2f"%(f_sim[ibin].GetParameter(3),f_sim[ibin].GetParError(3)) )
												pt_sim[ibin].SetTextColor(h.GetMarkerColor())
												pt_sim[ibin].SetFillColor(11)
												pt_sim[ibin].Draw()
								l[ibin].Draw()
															
							c.SaveAs("%s/c_q%0.2f_w%0.3f.png"%(outdir,q2bin_le,wbin_le))
							c.SaveAs("%s/c_q%0.2f_w%0.3f.eps"%(outdir,q2bin_le,wbin_le))
							c.Close()
		
		print "Done DispObs::plot_phiproj()"	
	
	def disp_1D(self):#,norm=False,q2min=1.25,q2max=5.25,wmin=1.400,wmax=2.125,dbg=False):
		"""
		+ Walk the ROOT file, per Q2-W bin:
			+ Extract h1(seq,vst,var)(=Obs_1D) where
				+ seq:
					+ =EC,EH if norm=True, 
					+ =SEQ_ALL
				+ var:	
					+ if vst==1:           VARS=(M1,THETA,ALPHA)
					+ if vst==2 or vst==3: VARS=(M2,THETA,ALPHA)

			+ If norm=True, then Normalize h1(seq,vst,var)
			+ Fill relevant bin of h2 (=Intg_yld)
		"""
		print "In DispObs::disp_1D()"

		#! Make OUTDIR_1D
		self.OUTDIR_1D=os.path.join(self.OUTDIR,"Obs_1D%s"%("_norm" if self.NORM==True else ""))
		if not os.path.exists(self.OUTDIR_1D):
			os.makedirs(self.OUTDIR_1D)

		#! + Get Q2-W binning from file
		#! + Q2-W binning in the file is as per Q2-W binning of h8s
		if self.DBG==True:
			q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX,dbg=True,dbg_bins=2)
		else:
			q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX)
		print "DispObs::disp_1D() Q2-W bins got from yield.root=",q2wbinl

		#! + Get binning information for Q2 and W, individually, from Q2-W binning
		q2bng=self.get_q2bng(q2wbinl)
		wbng=self.get_wbng(q2wbinl)
		print "DispObs::disp_1D() nq2bins,nwbins=",q2bng['NBINS'],wbng['NBINS']

		#! Make histogram to keep integrated yields in each Q2-W bin
		h2=self.make_h2_itg_yld(q2bng,wbng)

		print "DispObs::disp_1D() Going to begin processing Q2-W bins from file"
		#! Set seql to process
		if self.NORM==True:
			seql=['EC','EH']
		else:
			seql=SEQ_ALL
		#! Now, per Q2-W bin:
		#! + obtain h1 (=obs_1D)
		#! + Fill h2(=Intg_yld)
		for q2wbin in q2wbinl:
			print"Processing q2wbin=",q2wbin
			#! The following will be needed when filling h2
			iq2bin=q2bng['BINS'].index(self.get_q2bin_le(q2wbin))
			iwbin=wbng['BINS'].index(self.get_wbin_le(q2wbin))
			print "q2bin#,wbin#=",iq2bin+1,iwbin+1

			#! Obtain h1(seq,vst,var)
			h1=OrderedDict()
			for item in list(itertools.product(seql,VSTS)):
				seq,vst=item[0],item[1]
				if   vst==1:           varl=['M1','THETA','ALPHA']
				elif vst==2 or vst==3: varl=['M2','THETA','ALPHA']
				for var in varl:
					h1[seq,vst,var]=self.FIN.Get("%s/%s/VST%d/h1_%s"%(q2wbin,seq,vst,var))
			#! Plot h1(seq,vst,var)
			self.plot_1D(h1,q2wbin)



		
	def make_h2_itg_yld(self,q2bng,wbng):
		'''
		+ Make h2_itg_yld(Q2,W) with xaxis,yaxis=nwbins,nq2bins
			+ This is consistent with philosophy that binning is defined in h8s
			  and all subsequent analysis use that binning.
		'''
		#! First get binning in Q2 and W from Q2-W binning
		print "In DispObs::make_h2_itg_yld()"
		nq2bins=q2bng['NBINS']
		nwbins=wbng['NBINS']
		q2bins=array.array('f',range(1,nq2bins+1+1))#!Last '+1' to ensure 'range' includes top edge of last bin 
		wbins=array.array('f',range(1,nwbins+1+1))#!Last '+1' to ensure 'range' includes top edge of last bin
		h2=ROOT.TH2F("h_itg_yld","itg_yld(Q^{2},W)",nwbins,wbins,nq2bins,q2bins) 
		h2.SetXTitle("W (GeV)")
		h2.SetYTitle("Q^{2} (GeV^2)")
		#! Label the bins
		for iwbin in range(nwbins):
			le=wbng['BINS'][iwbin]
			ue=wbng['BINS'][iwbin+1]
			h2.GetXaxis().SetBinLabel(iwbin+1,"[%.3f,%.3f)"%(le,ue))
		for iq2bin in range(nq2bins):
			le=q2bng['BINS'][iq2bin]
			ue=q2bng['BINS'][iq2bin+1]
			h2.GetYaxis().SetBinLabel(iq2bin+1,"[%.2f,%.2f)"%(le,ue))
			
		print "make_h2_itg_yld() nq2bins,nwbins=%d,%d"%(h2.GetNbinsY(),h2.GetNbinsX())
		print "make_h2_itg_yld() wbin labels:" 
		for iwbin in range(nwbins):
			print "wbin#=%d,label=%s"%(iwbin+1,h2.GetXaxis().GetBinLabel(iwbin+1))
		print "make_h2_itg_yld() q2bin labels:" 
		for iq2bin in range(nq2bins):	
			print "q2bin#=%d,label=%s"%(iq2bin+1,h2.GetYaxis().GetBinLabel(iq2bin+1))
		return h2


	def disp_obs_R2(self,R2l,mthd,viewl,plotphiproj,q2min,q2max,wmin,wmax,dtypl,seql):
		"""
		+ Extract user specified R2s using used specified method(='mthd1'/'mthd2'/'mthd3' ='h5-mply-itg'/'phi-proj-fit'/'phi-proj-mply-itg')
			+ The user has to also specify:
				1. List of R2s to extract
				1. The Q2 & W limits i.e [q2min,q2max],[wmin,wmax]
				2. dtypl and seql 
		+ The following is the outline of the process:
			1. From file, get q2wbinl within specified q2min and q2max
			2. From file, get h5d[hel][(q2,w,dtyp,seql)]
			3. Extract hR2[hel][(q2,w,dtyp,seql)]
			4. Plot hR2[hel][(q2,w,dtyp,seql)]
		"""
		print "In DispObs::disp_obs_R2()"

		#self.OUTDIR_OBS_R2=os.path.join(self.OUTDIR,"Obs_R2",R2,"mthd_%s"%self.MTHD_NAMED[mthd],"Q2_%.2f-%.2f"%(q2min,q2max))
		self.OUTDIR_OBS_R2=os.path.join(self.OUTDIR,"Obs_R2","Q2_%.2f-%.2f"%(q2min,q2max),"mthd_%s"%self.MTHD_NAMED[mthd])
		if not os.path.exists(self.OUTDIR_OBS_R2):
			os.makedirs(self.OUTDIR_OBS_R2)

		#! 1. First get all q2wbin directories from file
		print "Getting q2wbinl"
		#q2wbinl=self.get_q2wbinlist(q2min=q2min,q2max=q2max,wmin=wmin,wmax=wmax,dbg=True,dbg_bins=10)
		q2wbinl=self.get_q2wbinlist(q2min=q2min,q2max=q2max,wmin=wmin,wmax=wmax)
		#print q2wbinl
		#! 1.1. Make a dictionary for the "bad" q2wbins
		q2wbinl_bad={}

		#! 2. Now get h5d
		print "Getting h5d"
		h5d=self.get_h5d(q2wbinl,q2wbinl_bad,dtypl,seql)
		print "Done getting h5d"
		# for hel in h5d:				
		# 	print "keys in h5d[%s]"%hel
		# 	print h5d[hel].keys()
		# ct=ROOT.TCanvas()
		# ct.Divide(1,3)
		# ct.cd(1)
		# h5d['POS'][1.25, 1.925, 3, 'EXP', 'C'].Projection(H5_DIM['PHI']).Draw()
		# ct.cd(2)
		# h5d['NEG'][1.25, 1.925, 3, 'EXP', 'C'].Projection(H5_DIM['PHI']).Draw()
		# ct.cd(3)
		# h5d['ASM'][1.25, 1.925, 3, 'EXP', 'C'].Projection(H5_DIM['PHI']).Draw()
		# ct.SaveAs("ctest.png")

		#! Write out the list of "bad" q2w bins and the reason why it is bad
		print "Writing out q2wbinl_bad"
		fout=open("%s/bad-q2wbins.txt"%self.OUTDIR_OBS_R2,"w")
		fout.write("Following are the \"bad\" q2w bins:\n")
		for k in q2wbinl_bad:
			fout.write("%s:%s\n"%(k,q2wbinl_bad[k]))
		fout.close()

		#! 3. and 4. Extract and plot R2
		print "Going to extract and plot R2 for mthd %s"%mthd
		self.extract_and_plot_R2(h5d,R2l,mthd,viewl,plotphiproj,dtypl,seql)
				
		print "Done DispObs::disp_obs_R2()"
		print "If the progam is not terminating, then Python is probably doing \"garbage collection\"(?); Wait a while!"
		return

	def disp_integ_yield(self,seql=['C','F'],vst='VST1',norm=False):
		"""
		Walk the ROOT file and plot y(w;seq,q2bin). 
		"""
		if norm==False:
			outdir=os.path.join(self.OUTDIR,"Obs_IntgYld_top%s"%(''.join(str(t) for t in self.TOPS)))
		else:
			outdir=os.path.join(self.OUTDIR,"Obs_IntgYld_norm_top%s"%(''.join(str(t) for t in self.TOPS)))
		if not os.path.exists(outdir):
			os.makedirs(outdir)

		#! 1. Get all q2wbins
		q2wbinl=self.get_q2wbinlist(q2max=3.25,wmin=1.400,wmax=2.125) #wbins to agree with Isupov
		#print q2wbinl

		q2bng=self.get_q2bng(q2wbinl)
		#print q2bng['BINW']
		print "Yield(w) will be plotted for the following Q2 bins:"
		print q2bng['BINS']

		#! 3. Now put together dictionary for yield: y{seq,q2bin:{w:yield}}
		#! First create dict
		y={}
		oy={} #Will be used later to store yields ordered
		for seq in seql:
			for q2wbin in q2wbinl:
				q2bin=q2wbin.split('_')[0]
				y[seq,q2bin]={}
				oy[seq,q2bin]={}
		#! Now fill dict
		for seq in seql:
			for q2wbin in q2wbinl:
				print "seq,q2wbin,norm=",seq,q2wbin,norm
				q2bin=q2wbin.split('_')[0]
				#! Get w, yield
				w=float(q2wbin.split('_')[1].split('-')[0])
				dw=float(q2wbin.split('_')[1].split('-')[1])-w
				#wbin=q2wbin.split('_')[1]
				h5_UNPOL=self.FEXP.Get("%s/%s/%s/h5_UNPOL"%(q2wbin,vst,seq))
				y[seq,q2bin][w]=thntool.GetIntegral(h5_UNPOL)
				if norm==True:
					q2=float(q2wbin.split('_')[0].split('-')[0])
					dq2=float(q2wbin.split('_')[0].split('-')[1])-q2
					#normf=LUM*LUM_INVFB_TO_INVMICROB*getvgflux(w,q2)*dw*dq2#!mub^-1
					normf=LUM*LUM_INVFB_TO_INVMICROB*getvgflux(w+(dw/2),q2+(dq2/2))*dw*dq2#!mub^-1
					#normf=LUM*LUM_INVFB_TO_INVMICROB*getvgflux(w+dw,q2+dq2)*dw*dq2#!mub^-1
					print "yield=",y[seq,q2bin][w]
					print "dw,dq2=",dw,dq2
					print "norm=",normf
					y[seq,q2bin][w]=y[seq,q2bin][w]/normf
					print "yield after norm=",y[seq,q2bin][w]
		#! Make sure y[seq,q2bin:w] are sorted by w
		for k in y.keys():
			oy[k]=OrderedDict(sorted(y[k].items()))
			
		#! 4. Now plot
		fig=plt.figure()
		ax=plt.subplot(111)
		clrd={'1.25-1.75':'red','1.75-2.25':'brown','2.25-2.75':'magenta','2.75-3.25':'orange',
		      '3.25-3.75':'yellow','3.75-4.25':'green','4.25-4.75':'cyan','4.75-5.25':'blue'}
		mrkrd={'C':'o','F':'^'}
		for k in oy.keys():
			seq=k[0]
			q2wbin=k[1]
			lbl='%s:[%s)'%(seq,q2wbin)
			ax.scatter(oy[k].keys(),oy[k].values(),label=lbl,c=clrd[q2wbin],marker=mrkrd[seq],s=50)
			ax.set_xticks(oy[k].keys())
			ax.set_xticklabels(oy[k].keys(),rotation='vertical')
		if norm==False:
			ax.set_ylim(0,600000)
			ax.set_ylabel(r'Yield [A.U.]',fontsize='xx-large')
		else:
			ax.grid(1)
			ax.set_ylim(0,20)
			ax.set_ylabel(r'$\mu b$',fontsize='xx-large')
		ax.legend()
		fig.savefig('%s/integ_yield_%s.png'%(outdir,vst))
		fig.savefig('%s/integ_yield_%s.eps'%(outdir,vst))	
			

	def get_sim_stats(self):
		"""
		Walk the ROOT file and obtain simstats(ss) for a h5_UNPOL in a Q2-W bin:
		ss={'T':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]],
		    'R':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]],
		    'A':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]],
		    'H':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]]}

		where:
			+ nbins=number of filled bins in a h5_UNPOL
			+ N=sum({n_i}) were n_i=events per bin (number of events in a h5_UNPOL)
			+ mu=average({n_i} (average of number of events per bin)
			+ sg=(RMS({n_i}) (RMS of number of events per bin)
		"""
		#! First get all q2wbin directories from file
		q2ws=self.get_q2wbinlist()
		#print "Processing sim_stats for %s:"%self.Q2W
		print q2ws

		ss={'T':[],'R':[],'A':[],'H':[]}
		f=ROOT.TFile(self.FSIM.GetName())
		for seq in ['T','R','A','H']:
			for q2w in q2ws:
				print "Processing sim_stats for %s"%q2w
				#! Determine q2,w
				q2bin=q2w.split('_')[0]
				wbin=q2w.split('_')[1]
				#print q2bin,wbin
				q2=float(q2bin.split('-')[0])
				w=float(wbin.split('-')[0])
				#print q2,w
				#! Determine nbins,N,mu,sg for this q2,w
				h5_UNPOL=f.Get("%s/VST1/%s/h5_UNPOL"%(q2w,seq))
				nbins=thntool.GetNbinsNotEq0(h5_UNPOL)
				N=thntool.GetIntegral(h5_UNPOL)
				binc_stats=np.zeros(2,'f')
				thntool.GetBinContentDistStats(h5_UNPOL,binc_stats)
				mu=binc_stats[0]
				sg=binc_stats[1]
				ss[seq].append([q2,w,nbins,N,mu,sg])
			# #! Compute average
			# ss[seq].append(nevts/len(q2ws))
			# ss[seq].append(nbins/len(q2ws))
		return ss

	def get_sim_stats_commonbins(self):
		"""
		The same as get_sim_stats(), but with the difference that mu,sg
		are obtained only for R-PS bins. I am currently not using this since
		this process takes impractically long.

		Walk the ROOT file and obtain simstats(ss) for a h5_UNPOL in a Q2-W bin:
		ss={'T':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]],
		    'R':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]],
		    'A':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]],
		    'H':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]]}

		where:
			+ nbins=number of filled bins in a h5_UNPOL
			+ N=sum({n_i}) were n_i=events per bin (number of events in a h5_UNPOL)
			+ mu=average({n_i}) (average of number of events per bin)
				+ Note, average computed over only R-PS bins
			+ sg=(RMS({n_i}) (RMS of number of events per bin)
				+ Note, average computed over only R-PS bins
		"""
		#! First get all q2wbin directories from file
		q2ws=self.get_q2wbinlist()
		print "Processing sim_stats for %s:"%self.Q2W
		print q2ws

		ss={'T':[],'R':[],'A':[],'H':[]}
		f=ROOT.TFile(self.FSIM.GetName())
		for q2w in q2ws:
			print "Processing %s..."%q2w
			#! Determine q2,w
			q2bin=q2w.split('_')[0]
			wbin=q2w.split('_')[1]
			#print q2bin,wbin
			q2=float(q2bin.split('-')[0])
			w=float(wbin.split('-')[0])
			#print q2,w
			#! First get all h5_UNPOLs
			h5_UNPOL={}
			for seq in ['T','R','A','H']:
				h5_UNPOL[seq]=f.Get("%s/VST1/%s/h5_UNPOL"%(q2w,seq))
			#! Now get simstats	
			for seq in ['T','R','A','H']:
				#! Determine nbins,N,mu,sg for this q2,w
				nbins=thntool.GetNbinsNotEq0(h5_UNPOL[seq])
				N=thntool.GetIntegral(h5_UNPOL[seq])
				binc_stats=np.zeros(2,'f')
				print "here"
				thntool.GetBinContentDistStatsCommonBins(h5_UNPOL[seq],h5_UNPOL['R'],binc_stats)
				print "here1"
				mu=binc_stats[0]
				sg=binc_stats[1]
				ss[seq].append([q2,w,nbins,N,mu,sg])
		return ss


	def get_q2wbinlist(self,q2min=0.00,q2max=6.00,wmin=0.000,wmax=3.000,dbg=False,dbg_bins=10):
		"""
		The ROOT file is arranged in a Tree Structure. The
		Observable histograms (obs-hists) are located as files in the following directory-path(dirpath):
		q2wbin/vst/seq/hists
		"""
		q2wbinl=[]
		#! File to use to extract q2wbinl
		f=self.FIN
				
		i=0 #! for dbg_bins
		for path,dirs,files in f.walk():
			if path=="":continue #! Avoid root path
			path_arr=path.split("/")
			if len(path_arr)==1:
				q2wbinl.append(path)
				i+=1
			if dbg==True:
				if i>=dbg_bins: break #! Uncomment/comment -> Get limited q2w-bins/Get all q2w-bins

		#! Remove q2wbins that are not within [q2min,q2max],[wmin,wmax] 
		q2wbins_remove=[]
		for q2wbin in q2wbinl:
			q2bin_le=q2wbin.split("_")[0].split("-")[0]
			q2bin_ue=q2wbin.split("_")[0].split("-")[1]
			wbin_le =q2wbin.split("_")[1].split("-")[0]
			wbin_ue =q2wbin.split("_")[1].split("-")[1]
			if float(q2bin_ue)<=q2min or float(q2bin_le)>=q2max or float(wbin_ue)<=wmin or float(wbin_le)>=wmax:
				q2wbins_remove.append(q2wbin)
		for q2wbin in q2wbins_remove:
			q2wbinl.remove(q2wbin)

		return q2wbinl


	def get_q2bng(self,q2wbinl):
		return self.get_xbng('Q2',q2wbinl)
	def get_wbng(self,q2wbinl):
		return self.get_xbng('W',q2wbinl)

	def get_xbng(self,x,q2wbinl):
		"""
		Gets x=(Q2 or W) binning information from q2wbinl
		"""

		#! 1. From q2wbinl, identify xbng 
		xbins=[]
		for q2wbin in q2wbinl:
			#print q2wbin
			if(x=='Q2'):
				xbins.append(float(q2wbin.split('_')[0].split('-')[0]))
				xbins.append(float(q2wbin.split('_')[0].split('-')[1]))
			elif(x=='W'):
				xbins.append(float(q2wbin.split('_')[1].split('-')[0]))
				xbins.append(float(q2wbin.split('_')[1].split('-')[1]))
			else:
				sys.exit("DispObs::get_xbng() x=%s not recognized"%x)
		xbins=set(xbins) #! Keep only unique entries
		xbins=list(xbins) #! convert 'set' back to 'list'
		xbins=sorted(xbins) #! Order entries
		#print xbins
		xbins_le=xbins[:-1]
		xbins_ue=xbins[1:]
		#print xbins_le
		#print xbins_ue
		xmin=min(xbins)
		xmax=max(xbins)
		xbinw=xbins[1]-xbins[0]
		nxbins=len(xbins_le)
		xbng={'MIN':xmin,'MAX':xmax,'BINW':xbinw,'NBINS':nxbins,
		   'BINS_LE':xbins_le,'BINS_UE':xbins_ue,'BINS':xbins}
		return xbng

	def get_q2bin(self,q2wbin):
		return q2wbin.split('_')[0]
	def get_wbin(self,q2wbin):
		return q2wbin.split('_')[1]
	def get_q2bin_le(self,q2wbin):
		return float(self.get_q2bin(q2wbin).split('-')[0])
	def get_wbin_le(self,q2wbin):
		return float(self.get_wbin(q2wbin).split('-')[0])
	# def get_q2bin_le(self,q2wbin):
	# 	return float(q2wbin.split('_')[0].split('-')[0])
	# def get_wbin_le(self,q2wbin):
	# 	return float(q2wbin.split('_')[1].split('-')[0])

	def get_h5d(self,q2wbinl,q2wbinl_bad,dtypl,seql):
		"""
		+ For given q2wbinl, from file, get and return:
			h5d[hel][(q2bin,wbin,dtyp,vst,dtyp,seq)]
			+ NOTE, currently h5d not filled if dtyp=='SIM' and seq=='C'
		+ Note bad q2wbins in q2wbinl_bad
		"""
		h5d={'UNP':{},'POS':{},'NEG':{},'ASM':{}}
		for q2wbin in q2wbinl:
			print "Getting h5 for",q2wbin
			#! First make sure this bin is "good"
			is_badbin=self.is_bad_q2wbin(q2wbin,q2wbinl_bad)
			if is_badbin: continue
			q2bin_le=self.get_q2bin_le(q2wbin)
			wbin_le =self.get_wbin_le(q2wbin)
			for dtyp in dtypl:
				if   dtyp=='EXP':f=self.FEXP_HEL
				elif dtyp=='SIM':f=self.FSIM
				for vst in self.VSTS:
					for seq in seql:
						if dtyp=='SIM' and seq=='C': continue
						if dtyp=='EXP':
							h5d['UNP'][q2bin_le,wbin_le,vst,dtyp,seq]=f.Get("%s/VST%d/%s/h5_%s"%(q2wbin,vst,seq,'UNPOL'))
							h5d['POS'][q2bin_le,wbin_le,vst,dtyp,seq]=f.Get("%s/VST%d/%s/h5_%s"%(q2wbin,vst,seq,'POS'))
							h5d['NEG'][q2bin_le,wbin_le,vst,dtyp,seq]=f.Get("%s/VST%d/%s/h5_%s"%(q2wbin,vst,seq,'NEG'))
							h5d['ASM'][q2bin_le,wbin_le,vst,dtyp,seq]=h5d['POS'][q2bin_le,wbin_le,vst,dtyp,seq].Clone()
							h5d['ASM'][q2bin_le,wbin_le,vst,dtyp,seq].Add(h5d['NEG'][q2bin_le,wbin_le,vst,dtyp,seq],-1)
							#! Set appropriate name and title for h5d['ASM']
							#! + original name(set in proc_yields.py)=h5_HEL(HEL=POS/NEG/UNPOL)
							#! + original title(set in proc_yields.py)=[q2min-q2max]_[wmin-wmax]_VSTX_SEQ_HEL=(POS/NEG/UNPOL)
							#! 
							orig_nm=h5d['ASM'][q2bin_le,wbin_le,vst,dtyp,seq].GetName().split("_")
							new_nm=orig_nm
							orig_nm[1]='ASM'
							orig_ttl=h5d['ASM'][q2bin_le,wbin_le,vst,dtyp,seq].GetTitle().split("_")
							new_ttl=orig_ttl
							new_ttl[4]='ASM'
							h5d['ASM'][q2bin_le,wbin_le,vst,dtyp,seq].SetName("%s_%s"%(new_nm[0],new_nm[1]))
							h5d['ASM'][q2bin_le,wbin_le,vst,dtyp,seq].SetTitle("%s_%s_%s_%s_%s"%(new_ttl[0],new_ttl[1],new_ttl[2],new_ttl[3],new_ttl[4]))
						elif dtyp=='SIM':
							h5d['UNP'][q2bin_le,wbin_le,vst,dtyp,seq]=f.Get("%s/VST%d/%s/h5_%s"%(q2wbin,vst,seq,'UNPOL'))
		return h5d					

	def hist_1D_athtcs(self,h1):
		for k in h1.keys():	
			seq,vst,var=k[0],k[1],k[2]
			h1[k].SetTitle("")
			h1[k].SetXTitle( "%s%s"%(VAR_NAMES[(vst,var)],VAR_UNIT_NAMES[var]) )
			h1[k].GetXaxis().SetLabelSize(.05)
			h1[k].GetXaxis().SetTitleSize(.10)
			h1[k].GetXaxis().SetTitleOffset(.7)
			h1[k].GetYaxis().SetTitleOffset(1.5)
			if self.NORM!=True:
				h1[k].SetMarkerColor(COLS_PLT_SEQ_ALL[seq])
				h1[k].SetLineColor(COLS_PLT_SEQ_ALL[seq])
				h1[k].SetMarkerStyle(MRKS_PLT_SEQ_ALL[seq])	
	
	def russ_norm_theta_dist(self,hTheta):
		#! 1. Create Russian-normalization factor histogram
		hDCosTheta=hTheta.Clone()
		nbins=hTheta.GetNbinsX()
		for ibin in range(nbins):
			theta_a=hTheta.GetBinLowEdge(ibin+1)
			theta_b=hTheta.GetBinLowEdge(ibin+1) + hTheta.GetBinWidth(ibin+1)
			DCosTheta=math.fabs(math.cos(math.radians(theta_b))-math.cos(math.radians(theta_a)))
			hDCosTheta.SetBinContent(ibin+1,DCosTheta)
			hDCosTheta.SetBinError(ibin+1,0.)
		#! Now divide hTheta by hDCosTheta
		#! Do Sumw2() so that errors are correctly propagated
		hTheta.Sumw2();
		hTheta.Divide(hDCosTheta)

	
	#! First make sure this bin is "good"
	def is_bad_q2wbin(self,q2wbin,q2wbinl_bad):
		is_badbin=False
		reason=''
		for vst in self.VSTS:
			h5_UNPOL=self.FEXP_HEL.Get("%s/VST%d/R/h5_UNPOL"%(q2wbin,vst))
			if thntool.GetIntegral(h5_UNPOL)==0:
				reason+="ER=0 for VST%d;"%vst
				is_badbin=True
		if is_badbin:
			q2wbinl_bad[q2wbin]=reason
		return is_badbin

	def get_q2bin_lel_wbin_lel(self,hd):
		"""
		From a histogram dictionary, which has multi-dimensional keys of which the first two 
		components are the low edges of q2bins and wbins, this function returns ordered
		lists, one each for the of low edges of q2bins and wbins. The lists contain only
		unique values and the values are ordered.
		
		@hd = histogram dictionary which has multi-dimensional keys. 
			+ Note, that for this function to work, the first two components
		      of the key i.e. key[0] and key[1], HAVE to be q2bin_le and wbin_le, respectively.
		"""
		q2bins_lel,wbins_lel=[],[]
		for k in hd.keys():
			q2bins_lel.append(k[0])
			wbins_lel.append(k[1])
		#! Keep only unique entries
		q2bins_lel=set(q2bins_lel) #! Keep only unique entries
		q2bins_lel=list(q2bins_lel) #! convert 'set' back to 'list'
		q2bins_lel=sorted(q2bins_lel)
		wbins_lel=set(wbins_lel) #! Keep only unique entries
		wbins_lel=list(wbins_lel) #! convert 'set' back to 'list'
		wbins_lel=sorted(wbins_lel)
		return q2bins_lel,wbins_lel

	def plot_obs_R2_athtcs(self):
		#! Stats Box
		ROOT.gStyle.SetOptStat(0)

		# # ROOT.gStyle.SetLabelSize(0.5,"t")
		# # ROOT.gStyle.SetTitleSize(0.5,"t")
		# #ROOT.gStyle.SetPaperSize(20,26);
		# ROOT.gStyle.SetPadTopMargin(0.15)#(0.05);
		# #ROOT.gStyle.SetPadRightMargin(0.09)#(0.05);
		# #ROOT.gStyle.SetPadBottomMargin(0.20)#(0.16);
		# #ROOT.gStyle.SetPadLeftMargin(0.15)#(0.12);

		ROOT.gStyle.SetTitleW(8)# //title width 
		ROOT.gStyle.SetTitleFontSize(10)# //title width 
		ROOT.gStyle.SetTitleH(0.10)# //title height 
		# ROOT.gStyle.SetTitleY(1)# //title Y location 

		#!get rid of X error bars and y error bar caps
		ROOT.gStyle.SetErrorX(0.001)

	def plot_obs_R2_view2_athtcs(self):
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

	def plot_phiproj_athtcs(self):
		#! Stats Box
		ROOT.gStyle.SetOptStat(0)
		#ROOT.gStyle.SetOptFit(1111)

		# # ROOT.gStyle.SetLabelSize(0.5,"t")
		# # ROOT.gStyle.SetTitleSize(0.5,"t")
		# #ROOT.gStyle.SetPaperSize(20,26);
		# ROOT.gStyle.SetPadTopMargin(0.15)#(0.05);
		# #ROOT.gStyle.SetPadRightMargin(0.09)#(0.05);
		# #ROOT.gStyle.SetPadBottomMargin(0.20)#(0.16);
		# #ROOT.gStyle.SetPadLeftMargin(0.15)#(0.12);

		ROOT.gStyle.SetTitleW(5)# //title width 
		ROOT.gStyle.SetTitleFontSize(10)# //title width 
		ROOT.gStyle.SetTitleH(0.10)# //title height 
		# ROOT.gStyle.SetTitleY(1)# //title Y location 

		#!get rid of X error bars and y error bar caps
		ROOT.gStyle.SetErrorX(0.001);
		return

	def plot_1D_athtcs(self):
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

	def get_signed_integral(self,h):
		itg=0
		for i in range(h.GetNbinsX()):
			itg+=math.fabs(h.GetBinContent(i+1))
		return itg

