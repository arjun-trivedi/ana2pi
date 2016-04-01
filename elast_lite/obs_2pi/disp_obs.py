from __future__ import division
import ROOT
from rootpy.io import root_open, DoesNotExist

from collections import OrderedDict
import array
import itertools

import os,sys
import time,datetime

import numpy as np
import matplotlib.pyplot as plt

import math

import atlib as atlib

from rootpy.interactive import wait

from proc_h8 import H5_DIM,VSTS,VARS,SEQ_ALL,MASS_P,MASS_PIP,MASS_PIM

ROOT.gROOT.ProcessLine(".L THnTool.C+")
from ROOT import THnTool

# Tools 
thntool=THnTool()

#! For Luminosity & vgflux
LUM_E1F=19.844 #fb^-1
E1F_E0=5.499

LUM_E16=28
E16_E0=5.754

LUM_INVFB_TO_INVMICROB=1000000000

PI=3.14159265358979312
FSC=0.00729735253
A=FSC
NA=6.02214129E23
QE=1.60217646E-19
MP=0.93827203

Q2BINW=0.5 #Gev^2
WBINW=0.025 # GeV

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

VAR_UNIT_NAMES={VARS[0]:"GeV",VARS[1]:"GeV",VARS[2]:"deg",VARS[3]:"deg",VARS[4]:"deg"}
VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC={VARS[0]:"GeV",VARS[1]:"GeV",VARS[2]:"",VARS[3]:"rad",VARS[4]:"rad"}

#! Following dictionaries for extracting R2
#! post VM-fdbk
#R2_NAMED={'A':'R2_{T}+R2_{L}','B':'R2_{LT}','C':'R2_{TT}','D':'R2_{LT^{p}}'} #'R2_{LT^{\'}}'
#R2_NAMED={'A':'R2_{T}+R2_{L}','B':'R2_{LT}','C':'R2_{TT}','D':'D','E':'E'}
R2_NAMED={'A':'R2_{T}+R2_{L}','B':'R2^{c}_{LT}','C':'R2^{c}_{TT}','D':'R2^{s}_{LT}','E':'R2^{s}_{TT}'}

H5_MFUNCD={'A':'1','B':'cphi','C':'c2phi','D':'sphi','E':'s2phi'}
MTHD_NAMED={'mthd1':'h5-mply-itg','mthd2':'phi-proj-fit','mthd3':'phi-prof-mply-itg'}
#! The following binning information taken from h8_bng.h
NBINS={'M1':14,'M2':14,'THETA':10,'PHI':10,'ALPHA':10}	
#! The following for phi-projection canvas
NXPADS={'M1':2,'M2':2,'THETA':2,'PHI':2,'ALPHA':2}
NYPADS={'M1':7,'M2':7,'THETA':5,'PHI':5,'ALPHA':5}
#! EF and CF needed for each mthd to finally extract R2s (see $ELAST_LITE/obs_2pi/dev/R2_EF_CFs.txt)
EF={('mthd1','A'):2*math.pi,('mthd1','B'):math.pi,('mthd1','C'):math.pi,('mthd1','D'):math.pi,('mthd1','E'):math.pi,
	('mthd2','A'):1,        ('mthd2','B'):1,      ('mthd2','C'):1,      ('mthd2','D'):1,      ('mthd2','E'):1}
CF={('mthd1','A'):1,      ('mthd1','B'):0.93549,('mthd1','C'):0.75683,('mthd1','D'):0.93549,('mthd1','E'):0.75683,
	('mthd2','A'):1.59155,('mthd2','B'):1.61803,('mthd2','C'):1.70131,('mthd2','D'):1.61803,('mthd2','E'):1.70131}

def get_fphi(name_sfx=""):
	fphi=ROOT.TF1("fphi%s"%name_sfx, "([0] + [1]*cos(x*TMath::DegToRad()) + [2]*cos(2*x*TMath::DegToRad()))",0,360)

	#! Set parameter names
	fphi.SetParName(0, "A")
	fphi.SetParName(1, "B")
	fphi.SetParName(2, "C")
	
	#! Set initial state of parameters
	fphi.SetParameter(0,.1)#!
	fphi.SetParameter(1,.1)#10
	fphi.SetParameter(2,.1)#20
	
	return fphi

def get_fphi_alpha(name_sfx=""):
	'''
	Add sphi and s2phi moments for phi dependence in alpha projection
	'''
	fphi=ROOT.TF1("fphi%s"%name_sfx, "([0] + [1]*cos(x*TMath::DegToRad()) + [2]*cos(2*x*TMath::DegToRad()) + [3]*sin(x*TMath::DegToRad()) + [4]*sin(2*x*TMath::DegToRad()))",0,360)
	
	#! Set parameter names
	fphi.SetParName(0, "A")
	fphi.SetParName(1, "B")
	fphi.SetParName(2, "C")
	fphi.SetParName(3, "D")
	fphi.SetParName(4, "E")

	#! Set initial state of parameters
	fphi.SetParameter(0,.1)#!
	fphi.SetParameter(1,.1)#10
	fphi.SetParameter(2,.1)#20
	fphi.SetParameter(3,.1)#100
	fphi.SetParameter(4,.1)#100
	
	return fphi

#! Constants relating to plotting and its aesthetics
PLT_SEQ_ALL=SEQ_ALL#['ST','SR','SA','SC','SH','SF','ER','EC','EH','EF']
PLT_SEQ_ALL.remove('SA') #! SA is not plotted, at least not as an Observable

CLRS_PLT_SEQ_ALL={('ST'):ROOT.gROOT.ProcessLine("kGreen"),
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

#! N unique identifiers in the space of markers,colors
#! Currently N=42: 7-markers * 6-colors
#! + This is useful when plotting several histograms on the same TCanvas
#! + This method can be extended by adding dimensions to the space and more points along each dimension
MRKS_N=[ROOT.gROOT.ProcessLine("kFullCircle"),ROOT.gROOT.ProcessLine("kFullTriangleUp"),
       ROOT.gROOT.ProcessLine("kFullTriangleDown"),ROOT.gROOT.ProcessLine("kFullSquare"),
       ROOT.gROOT.ProcessLine("kFullCross"),ROOT.gROOT.ProcessLine("kFullStar"),
       ROOT.gROOT.ProcessLine("kFullDiamond")]
CLRS_N=[ROOT.gROOT.ProcessLine("kRed"),ROOT.gROOT.ProcessLine("kMagenta"),ROOT.gROOT.ProcessLine("kBlue"),
        ROOT.gROOT.ProcessLine("kCyan"),ROOT.gROOT.ProcessLine("kGreen"),ROOT.gROOT.ProcessLine("kYellow")]

#! Get all coordinates in this space, where coordinate def.eq. (mrkr,clr)
MRKS_CLRS_N_COORDS=list(itertools.product(MRKS_N,CLRS_N))
#! Now put them into a dictionary labelled by numbers
MRKS_CLRS_N_COORDS_DICT=OrderedDict()
for i,coord in enumerate(MRKS_CLRS_N_COORDS):
	MRKS_CLRS_N_COORDS_DICT[i+1]=[coord[0],coord[1]]




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
	def __init__(self,obsdir,simnum='siml',view="norm",q2min=1.25,q2max=5.25,wmin=1.400,wmax=2.125,expt='e16',dbg=False,**SSargs):
		print "*** In DispObs::_init_() ***"
		
		print "expt=",expt

		#! setup self.DO_SS_STUDY
		self.DO_SS_STUDY=False
		if len(SSargs)>0: 
			print "SSargs exist. Will perform SS. \n"
			#print SSargs
			self.DO_SS_STUDY=True

		#! init independent of either SS or Regular
		self.Q2MIN,self.Q2MAX=q2min,q2max
		self.WMIN,self.WMAX=wmin,wmax
		if expt=='e1f':
			self.LUM=LUM_E1F
			self.E0=E1F_E0
		elif expt=='e16':
			self.LUM=LUM_E16
			self.E0=E16_E0

		self.DBG=dbg

		print "Q2MIN,Q2MAX=",self.Q2MIN,self.Q2MAX
		print "WMIN,WMAX=",self.WMIN,self.WMAX
		print "LUM=",self.LUM
		print "E0=",self.E0
		print "DBG=",self.DBG

		#! Now do init specific to SS or Regular
		if self.DO_SS_STUDY: #! Do SS related init
			#! set self.FIN=None
			#! + This is for technical reasons since in Regular init there is a self.FIN
			self.FIN=None

			#! Setup OBSDIR
			sfx=''
			if expt=='e16':sfx='_E16'
			self.OBSDIR=os.environ['OBSDIR%s'%sfx]

			#! Create OBSD
			self.OBSD=SSargs['obsd']

			#! Set up SEQS
			self.SEQS=['EC','EF']

			#! OUTDIR
			SS_tag=SSargs['SS_tag']
			date=datetime.datetime.now().strftime('%m%d%y')
			self.OUTDIR=os.path.join(self.OBSDIR,'SS','%s_%s'%(SS_tag,date))
			if self.DBG==True:
				self.OUTDIR=os.path.join(self.OBSDIR,'SS','dbg','%s_%s'%(SS_tag,date))

			print "OBSDIR=",self.OBSDIR
			#! pretty-print OBSD
			print "OBSD (pretty print):"
			for k in self.OBSD:
				print k,":",self.OBSD[k][0],self.OBSD[k][1]
			print "SEQS=",self.SEQS
			print "OUTDIR=",self.OUTDIR
			
		else: #! do Regular init	
			self.OBSDIR=obsdir
			self.SIM_NUM=simnum
			self.VIEW=view
			# self.Q2MIN,self.Q2MAX=q2min,q2max
			# self.WMIN,self.WMAX=wmin,wmax
			# if expt=='e1f':
			# 	self.LUM=LUM_E1F
			# 	self.E0=E1F_E0
			# elif expt=='e16':
			# 	self.LUM=LUM_E16
			# 	self.E0=E16_E0

			# self.DBG=dbg
				
			self.FIN=root_open(os.path.join(self.OBSDIR,self.SIM_NUM,'yield.root'))
			self.OUTDIR=os.path.join(self.OBSDIR,self.SIM_NUM)
			if self.DBG==True:
				self.OUTDIR=os.path.join(self.OBSDIR,self.SIM_NUM,'dbg')

			#! Set SEQS to process
			if   self.VIEW=="norm":
				self.SEQS=['EC','EF']
			elif self.VIEW=="fullana":
				self.SEQS=SEQ_ALL
			elif self.VIEW=="ERyield":
				self.SEQS=['ER']
			elif self.VIEW=="EC_SF_ST": #! views from this point on tuned for 'dobs_R2'
				self.SEQS=['EC','SF','ST'] 
			elif self.VIEW=="EC_EF":
				self.SEQS=['EC','EF']
			else:
				sys.exit("view={norm,fullana,ERyield,EC_SF_ST,EC_EF} only")

			print "OBSDIR=%s\nFIN=%s\nOUTDIR=%s"%(self.OBSDIR,self.FIN,self.OUTDIR)
			print "self.VIEW=",self.VIEW
			print "SEQS to process=",self.SEQS
			print "expt=",expt
			time.sleep(3)

		#! Some general ROOT related setup
		#! + Print only Warning and Error messages
		#! + Not print Info messages, for examples, from TCanvas::SaveAs()
		#! source: https://root.cern.ch/phpBB3/viewtopic.php?t=8787
		ROOT.gErrorIgnoreLevel=ROOT.gROOT.ProcessLine('kWarning')

		print "*** _init_() Done *** \n"

	def plot_1D(self,h1,q2wbin):
		"""
		"""
		print "In DispObs::plot_1D()"

		#! Set up plotting aesthetics
		self.plot_1D_athtcs()
		self.hist_1D_athtcs(h1)
				
		#! TCanvas's pad_map[pad,vst,var] defined as per Gleb's display
		pad_map=[(1,1,"M1"),   (2,3,'M2'),   (3,2,'M2'),
				 (4,1,"THETA"),(5,3,'THETA'),(6,2,'THETA'),
				 (7,1,"ALPHA"),(8,3,'ALPHA'),(9,2,'ALPHA')]
						  
		c=ROOT.TCanvas("c","c",1000,1000)
		pad_t=ROOT.TPad("pad_l","Legend pad",0.25,0.935,0.75,1.00)
		pad_p=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
		if (self.VIEW=="fullana"):pad_p.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
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
			if self.VIEW=="norm" or self.VIEW=="ERyield":
				#! First set minimum(=0) and maximum of y-axis
				fctr=1.1
				maxl=[h1[seq,vst,var].GetMaximum() for seq in self.SEQS]
				maximum=max(maxl)
				for htmp in [h1[seq,vst,var] for seq in self.SEQS]:
					htmp.SetMinimum(0.)
					htmp.SetMaximum(maximum*fctr)
				#! Now draw
				for i,seq in enumerate(self.SEQS):
					draw_opt="" if i==0 else "same"
					h1[seq,vst,var].Draw(draw_opt)
			elif self.VIEW=="fullana":#! Draw all hists normed to same integral (since I want to compare their distributions)
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

		print "*** plot_1D() Done ***\n"
		return

	def plot_1D_SS(self,h1,h1_SS_err,q2wbin):
		'''
		[02-29-16]
		+ 'plot_1D' was used as a tempate for this method
		+ The idea is extend the method to incorporate plotting of a 
		  list of SS-Obs and the associated SS-err
		+ Therefore only details relevant to the above were changed and the core functionality left intact

		+ Process details:
			+ NOTE that 1D-Obs and SS_err_1D_Obs are plotted for each of [EC,EH]
		'''

		print "In DispObs::plot_1D_SS()"

		#! Set up plotting aesthetics
		self.plot_1D_athtcs()
		self.hist_1D_athtcs(h1)
		self.hist_1D_athtcs(h1_SS_err)
				
		#! TCanvas's pad_map[pad,vst,var] defined as per Gleb's display
		pad_map=[(1,1,"M1"),   (2,3,'M2'),   (3,2,'M2'),
				 (4,1,"THETA"),(5,3,'THETA'),(6,2,'THETA'),
				 (7,1,"ALPHA"),(8,3,'ALPHA'),(9,2,'ALPHA')]
			
		#! Make plots for each of [EC,EH]
		#! + Setup variable: nlegend_canvas=0
		#! + The idea is to draw one master legend that will be used for all plots, Obs-1D and intg_xsec
		#! + This is because this legend is too large to draw on the same canvas as the plots.
		#! + Once this master legend in drawn and saved, no other such legend will be drawn.
		#! nlegend canvas
		nlegend_canvas=0	
		for seq in self.SEQS:			  
			#! Canvas for h1
			c=ROOT.TCanvas("c","c",1000,1000)
			pad_t_c=ROOT.TPad("pad_l","Title pad",0.25,0.935,0.75,1.00)
			pad_p_c=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
			pad_p_c.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
			pad_p_c.Draw()
			pad_t_c.Draw()
			pad_t_c.cd()
			pt=ROOT.TPaveText(.05,.1,.95,.8)
			pt.AddText("Q2_W bin=%s"%q2wbin)
			pt.SetTextSize(0.40)
			pt.Draw()
			pad_p_c.Divide(3,3)

			#! Canvas for h1_SS_err
			cerr=ROOT.TCanvas("cerr","cerr",1000,1000)
			pad_t_cerr=ROOT.TPad("pad_l","Title pad",0.25,0.935,0.75,1.00)
			pad_p_cerr=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
			pad_p_cerr.Draw()
			pad_t_cerr.Draw()
			pad_t_cerr.cd()
			# pt=ROOT.TPaveText(.05,.1,.95,.8)
			# pt.AddText("Q2_W bin=%s"%q2wbin)
			# pt.SetTextSize(0.40)
			pt.Draw()
			pad_p_cerr.Divide(3,3)

			for item in pad_map:
				pad,vst,var=item[0],item[1],item[2]
				print "pad,vst,var=",pad,vst,var
						
				#! Draw h1
				#! First set minimum and maximum of hists
				fctr=1.1
				maxl=[h1[obsnum,seq,vst,var].GetMaximum() for obsnum in self.OBSD.keys()]
				maximum=max(maxl)
				for obsnum in self.OBSD.keys():
					h1[obsnum,seq,vst,var].SetMinimum(0.)
					h1[obsnum,seq,vst,var].SetMaximum(fctr*maximum)
				#! Now draw
				gpad=pad_p_c.cd(pad)
				for i,obsnum in enumerate(self.OBSD.keys()):
					draw_opt="" if i==0 else "same"
					h1[obsnum,seq,vst,var].Draw(draw_opt)
				#! Add TLegend if pad==1
				# if pad==1:
				# 	l=ROOT.TLegend(0.60,0.60,1.00,0.90)
				# 	l.SetFillStyle(0)
				# 	l.SetBorderSize(0)
				# 	l.SetTextSize(0.04)
				# 	for k in self.OBSD.keys():
				# 		obsnum,obs,obstag=k,self.OBSD[k][0],self.OBSD[k][1]
				# 		l.AddEntry(h1[obsnum,seq,vst,var],"%s_%s"%(obstag,seq),"p")#EF
				# 	l.Draw()
				#! legend canvas
				if nlegend_canvas==0:
					#! Setup legend
					lm=ROOT.TLegend(0,0,1,1)
					lm.SetFillStyle(0)
					lm.SetBorderSize(0)
					lm.SetTextSize(0.03)
					for k in self.OBSD.keys():
						obsnum,obs,obstag=k,self.OBSD[k][0],self.OBSD[k][1]
						lm.AddEntry(h1[obsnum,seq,vst,var],"%s_%s"%(obstag,seq),"p")#EF
					
					clm=ROOT.TCanvas()
					clm.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
					lm.Draw()
					clm.SaveAs("%s/c_master_legend.png"%self.OUTDIR)
					clm.SaveAs("%s/c_master_legend.eps"%self.OUTDIR)
					clm.Close()
					nlegend_canvas+=1

				#! Draw h1_SS_err
				gpad=pad_p_cerr.cd(pad)
				h1_SS_err[seq,vst,var].SetMinimum(0)
				h1_SS_err[seq,vst,var].Draw()

			#! Save canvases for h1 and h1_SS_err
			q2bin=self.get_q2bin(q2wbin)
			wbin=self.get_wbin(q2wbin)
			#! h1 in self.OUTDIR_1D/seq/q2bin
			outdir_seq=os.path.join(self.OUTDIR_1D,seq)
			outdir_q2bin=os.path.join(outdir_seq,"q%s"%q2bin)
			if not os.path.exists(outdir_q2bin):
				os.makedirs(outdir_q2bin)
			if not os.path.exists(outdir_seq):
				os.makedirs(outdir_seq)
			c.SaveAs("%s/c_w%s_q%s.png"%(outdir_q2bin,wbin,q2bin))
			c.SaveAs("%s/c_w%s_q%s.eps"%(outdir_q2bin,wbin,q2bin))
			#! h1_SS_err in self.OUTDIR_1D_SS_ERR/seq/q2bin
			outdir_seq_err=os.path.join(self.OUTDIR_1D_SS_ERR,seq)
			outdir_q2bin_SS_err=os.path.join(outdir_seq_err,"q%s"%q2bin)
			if not os.path.exists(outdir_q2bin_SS_err):
				os.makedirs(outdir_q2bin_SS_err)
			if not os.path.exists(outdir_seq_err):
				os.makedirs(outdir_seq_err)
			cerr.SaveAs("%s/cerr_w%s_q%s.png"%(outdir_q2bin_SS_err,wbin,q2bin))
			cerr.SaveAs("%s/cerr_w%s_q%s.eps"%(outdir_q2bin_SS_err,wbin,q2bin))

		print "*** plot_1D_SS() Done ***\n"
		return

	def disp_1D(self):
		"""
		1. Get Q2-W bin list from yield.root
		2. Create h2(seq,vst,var) to store integrated yields=f(Q2-bin,W-bin)
		3. Per Q2-W bin:
			+ Extract h1(seq,vst,var)(=Obs_1D) where
				+ seq:
					+ =EC,EF       if view="norm" 
					+ =PLT_SEQ_ALL if view="fullana"
					+ =ER          if view="ERyield"
				+ var:	
					+ if vst==1:           VARS=(M1,THETA,ALPHA)
					+ if vst==2 or vst==3: VARS=(M2,THETA,ALPHA)

			+ If self.VIEW="norm", then Normalize h1(seq,vst,var)
			+ Fill relevant bin of h2 (=Intg_yld)
			+ Plot h1(seq,vst,var)
		4. Plot h2(seq,vst,var) and its projection on W
		"""
		print "*** In DispObs::disp_1D() ***"

		#! Make OUTDIR_1D
		#self.OUTDIR_1D=os.path.join(self.OUTDIR,"Obs_1D%s"%("_norm" if self.VIEW==True else ""))
		self.OUTDIR_1D=os.path.join(self.OUTDIR,"Obs_1D_%s"%self.VIEW)
		if not os.path.exists(self.OUTDIR_1D):
			os.makedirs(self.OUTDIR_1D)

		#! Make OUTDIR_2D
		#self.OUTDIR_ITG_YLD=os.path.join(self.OUTDIR,"Obs_Itg_Yld%s"%("_norm" if self.VIEW==True else ""))
		self.OUTDIR_ITG_YLD=os.path.join(self.OUTDIR,"Obs_Itg_Yld_%s"%self.VIEW)
		if not os.path.exists(self.OUTDIR_ITG_YLD):
			os.makedirs(self.OUTDIR_ITG_YLD)

		#! + Get Q2-W binning from file
		#! + Q2-W binning in the file is as per Q2-W binning of h8s
		if self.DBG==True:
			q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX,dbg=True,dbg_bins=2)
		else:
			q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX)
		print "DispObs::disp_1D() Q2-W bins got from yield.root=",q2wbinl

		#! Place holder for h2[seq,vst,var]: histogram to keep integrated yields in each Q2-W bin
		#! + Implementation appears sloppy for now, but it is created inside Q2-W bin loop
		h2=None
		#! + Get binning information for Q2 and W, individually, from Q2-W binning
		#! + This is needed later to specify the Q2- and W-binning when creating h2
		q2bng=self.get_q2bng(q2wbinl)
		wbng=self.get_wbng(q2wbinl)
		print "DispObs::disp_1D() nq2bins,nwbins=",q2bng['NBINS'],wbng['NBINS']

		print "DispObs::disp_1D() Going to begin processing Q2-W bins from file"
		#! Now, per Q2-W bin:
		#! + Obtain h1 (=obs_1D)
		#! + Create, if h2==None, and fill h2(=Intg_yld)
		for q2wbin in q2wbinl:
			print"Processing q2wbin=",q2wbin
			#! The following will be needed when filling h2
			iq2bin=q2bng['BINS'].index(self.get_q2bin_le(q2wbin))
			iwbin=wbng['BINS'].index(self.get_wbin_le(q2wbin))
			print "q2bin#,wbin#=",iq2bin+1,iwbin+1

			#! Obtain h1(seq,vst,var)
			h1=OrderedDict()
			for item in list(itertools.product(self.SEQS,VSTS)):
				seq,vst=item[0],item[1]
				if   vst==1:           varl=['M1','THETA','ALPHA']
				elif vst==2 or vst==3: varl=['M2','THETA','ALPHA']
				for var in varl:
					h1[seq,vst,var]=self.FIN.Get("%s/%s/VST%d/h1_%s"%(q2wbin,seq,vst,var))
					#! Call Sumw2() since the errors as currently in the h1s
					#! need to be propagated in all subsequent calculations that involve them
					h1[seq,vst,var].Sumw2()
			#! End of [seq,vst,var] loop

			#! If norm, then normalize h1(SEQ,VSTS,VARS)
			if self.VIEW=="norm":
				self.norm_1D(h1,q2wbin)
			else: #! just DCosTheta normalize theta distributions
				for k in h1:
					if k[2]=='THETA':
						self.norm_1D_theta(h1[k])

			#! Fill h2(seq,vst,var)
			if h2==None:
				h2=self.make_h2_itg_yld(q2bng,wbng,h1) #! Note h1 passed as a template to provide (seq,vst,var) keys
			self.fill_h2_itg_yld(h2,h1,iq2bin,iwbin)

			#! Plot h1(seq,vst,var)
			self.plot_1D(h1,q2wbin)
		#! End of q2wbin loop

		#! Plot h2
		self.plot_h2_itg_yld(h2)

		print "*** disp_1D() Done ***\n"


	def disp_1D_SS(self):
		"""
		[02-29-16]
		+ 'disp_1D' was used as a tempate for this method
		+ The idea is extend the method to incorporate a list of Observables representing various SS
		+ Therefore only details relevant to the above were changed and the core functionality left intact

		0. Make all needed output directories
		1. Get Q2-W bin list from any of the files from self.OBSD
		2. Create h2(obs,seq,vst,var) to store integrated yields=f(Q2-bin,W-bin)
		3. For every Q2-W bin*obsnum (in self.OBSD):
			3.i. Extract h1(obs,seq,vst,var)(=Obs_1D) where
				+ seq=EC,EF        
				+ var:	
					+ if vst==1:           VARS=(M1,THETA,ALPHA)
					+ if vst==2 or vst==3: VARS=(M2,THETA,ALPHA)
			3.ii.  Normalize h1(obs,seq,vst,var)
			3.iii. Get h1_SS_err(seq,vst,var) from h1(obs,seq,vst,var)
			3.iv.  Fill relevant bin of h2 (=Intg_yld)
			3.v.   Plot h1 and h1_SS_err(obs,seq,vst,var)
		4. Plot h2(obsnum,seq,vst,var) and its projection on W, which is the itg_xsec(W;Q2). 
		   The function that plots h2 and its projection on W i.e. the itg_xsec, also obtains
		   the SS error on the itg_xsec and plots it too.
		"""
		print "*** In DispObs::disp_1D_SS() ***"

		#! 0. Make output directories
		#! Make OUTDIR_1D,OUTDIR_1D_SS_ERR
		self.OUTDIR_1D=os.path.join(self.OUTDIR,"Obs_1D")
		self.OUTDIR_1D_SS_ERR=os.path.join(self.OUTDIR,"Obs_1D_SS_err")
		if not os.path.exists(self.OUTDIR_1D):
			os.makedirs(self.OUTDIR_1D)
		if not os.path.exists(self.OUTDIR_1D_SS_ERR):
			os.makedirs(self.OUTDIR_1D_SS_ERR)

		#! Make OUTDIR_ITG_YLD,OUTDIR_ITG_YLD_SS_ERR
		self.OUTDIR_ITG_YLD=os.path.join(self.OUTDIR,"Obs_Itg_Yld")
		self.OUTDIR_ITG_YLD_SS_ERR=os.path.join(self.OUTDIR,"Obs_Itg_Yld_SS_err")
		if not os.path.exists(self.OUTDIR_ITG_YLD):
			os.makedirs(self.OUTDIR_ITG_YLD)
		if not os.path.exists(self.OUTDIR_ITG_YLD_SS_ERR):
			os.makedirs(self.OUTDIR_ITG_YLD_SS_ERR)

		#! 1. Get Q2-W binning from file
		#! + Q2-W binning in the file is as per Q2-W binning of h8s
		if self.DBG==True: q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX,dbg=True,dbg_bins=2)
		else:              q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX)
		print "DispObs::disp_1D_SS() Q2-W bins got from yield.root=",q2wbinl

		#! 2. Create h2[obsnum,seq,vst,var]: histogram to keep integrated yields in each Q2-W bin
		#! + Implementation appears sloppy for now, but it is created inside Q2-W bin loop
		h2=None

		#! 3. Prepare to loop over Q2-W bin*obsnum (in self.OBSD):
		#! + Get binning information for Q2 and W, individually, from Q2-W binning
		#! + This is needed later to specify the Q2- and W-binning when creating h2
		q2bng=self.get_q2bng(q2wbinl)
		wbng=self.get_wbng(q2wbinl)
		print "DispObs::disp_1D() nq2bins,nwbins=",q2bng['NBINS'],wbng['NBINS']

		#! Now, process substeps of 3. for every Q2-W bin*obsnum (in self.OBSD):
		print "DispObs::disp_1D() Going to begin processing Q2-W bins from file"
		for q2wbin in q2wbinl:
			print"Processing q2wbin=",q2wbin

			#! 3.i. Create h1(obsnum,seq,vst,var)
			h1=OrderedDict()
			for k in self.OBSD:
				obsnum,obs,obstag=k,self.OBSD[k][0],self.OBSD[k][1]
				print "Processing obsnum:obs:obstag=",obsnum,obs,obstag

				#! Open file for this obs
				obsdir=os.path.join(self.OBSDIR,obs)
				f=root_open(os.path.join(obsdir,'yield.root'))
			
				#! The following will be needed when filling h2
				iq2bin=q2bng['BINS'].index(self.get_q2bin_le(q2wbin))
				iwbin=wbng['BINS'].index(self.get_wbin_le(q2wbin))
				print "q2bin#,wbin#=",iq2bin+1,iwbin+1

				for item in list(itertools.product(self.SEQS,VSTS)):
					seq,vst=item[0],item[1]
					if   vst==1:           varl=['M1','THETA','ALPHA']
					elif vst==2 or vst==3: varl=['M2','THETA','ALPHA']
					for var in varl:
						h1[obsnum,seq,vst,var]=f.Get("%s/%s/VST%d/h1_%s"%(q2wbin,seq,vst,var))

						#! + The following line was added to decouple histograms from the open
						#! ROOT file which can then be closed.
						#! + This had to be done because otherwise I was getting an error from
						#! rootpy saying that "Too many files are open" (perhaps related to 'ulimit': See notes in ZIM under Technology/UNIX):
						#! rootpy.ROOTError: level=5000, loc='TFile::TFile', msg='file /data/trivedia/e16/2pi/d2pi/highQ2_030916/cutsncors16/sim9/yield.root can not be opened for reading (Too many open files)'
						h1[obsnum,seq,vst,var].SetDirectory(0)

						#! Call Sumw2() since the errors as currently in the h1s
						#! need to be propagated in all subsequent calculations that involve them
						h1[obsnum,seq,vst,var].Sumw2()
				#! End of [seq,vst,var] loop

				#! Close file
				#! + Was getting intermittent error from rootpy saying "Too many files are open": 
				#! rootpy.ROOTError: level=5000, loc='TFile::TFile', msg='file /data/trivedia/e16/2pi/d2pi/highQ2_030916/cutsncors16/sim9/yield.root can not be opened for reading (Too many open files)'
				f.Close() 
			#! End of OBSD loop

			#! 3.ii. Normalize h1(obsnum,seq,vst,var)
			self.norm_1D(h1,q2wbin)

			#! 3.iii. Get h1_SS_err(seq,vst,var) from h1(obs,seq,vst,var)
			h1_SS_err=self.get_SS_err_1D(h1)

			#! 3.iv. Fill h2(obsnum,seq,vst,var)
			if h2==None:
				h2=self.make_h2_itg_yld(q2bng,wbng,h1) #! Note h1 passed as a template to provide (obsnum,seq,vst,var) keys
			self.fill_h2_itg_yld(h2,h1,iq2bin,iwbin)

			#! 3.v. Plot h1(obsnum,seq,vst,var)
			self.plot_1D_SS(h1,h1_SS_err,q2wbin)
		#! End of q2wbin loop

		#! 4. Plot h2
		self.plot_h2_itg_yld_SS(h2)

		print "*** disp_1D_SS() Done ***\n"
		print "DispObs::disp_1D() If not exiting immediately then Python is probably doing Garbage Collection which can take a while"

	def get_SS_err_1D(self,h1):
		'''
		[03-01-16]
		+ General comments:
			+ NOTE that the method actually returns the *relative* SS error.
			+ NOTE that the calculation of the error on this *relative* SS error is not yet implemented.
			+ Currently, h1 is either 1D-Obs or itg_xsec(W) and have keys [osbnum,seq,vst,var]

		+ Details:
			+ Using h1[osbnum,seq,vst,var] calculate and return h1_SS_err[seq,vst,var]
			+ Currently, the relative SS error is calcualated as follows:
				1. Calcualate arithmetic mean for all entries in a bin (mu_bin)
				2. Calcualate arithmetic sigma for all entries in a bin (sg_bin)
				3. Relative error = (sg_bin/mu_bin)*100
		  	+ In that sense, I still have to find a way to obtain the error of this SS error.
		'''
		h1_SS_err=OrderedDict()
		for item in list(itertools.product(self.SEQS,VSTS)):
			seq,vst=item[0],item[1]
			if   vst==1:           varl=['M1','THETA','ALPHA']
			elif vst==2 or vst==3: varl=['M2','THETA','ALPHA']
			for var in varl:
				#! Create h1_SS_err[seq,vst,var] by Cloning and then resetting
				#! the contents of the h1 of 1st SS-Obs in self.OBSD i.e. h1[1,seq,vst,var]
				h1_SS_err[seq,vst,var]=h1[1,seq,vst,var].Clone()
				h1_SS_err[seq,vst,var].Reset()
				#! Now loop over all bins can calculate SS err
				nbins=h1_SS_err[seq,vst,var].GetNbinsX()
				for ibin in range(nbins):
					binc,binerr=[],[]
					for obsnum in self.OBSD:
						binc.append(h1[obsnum,seq,vst,var].GetBinContent(ibin+1))
						binerr.append(h1[obsnum,seq,vst,var].GetBinError(ibin+1))
					SS_mu_bin=np.mean(binc)
					SS_sg_bin=np.std(binc)
					if SS_mu_bin==0: SS_rel_err=0
					else:            SS_rel_err=(SS_sg_bin/SS_mu_bin)*100
					h1_SS_err[seq,vst,var].SetBinContent(ibin+1,SS_rel_err)
					#! Till error on SS err is calculated, set err=0 (actually to an insignificantly small value so that ROOT can draw hist with "e" option)
					h1_SS_err[seq,vst,var].SetBinError(ibin+1,0.00000001)
		return h1_SS_err

	def norm_1D(self,h1,q2wbin):
		print "*** In DispObs::norm_1D() ***"
		
		#! Get virtual photon flux for this q2wbin
		q2bin_le=self.get_q2bin_le(q2wbin)
		wbin_le=self.get_wbin_le(q2wbin)
		# q2=q2bin_le+(Q2BINW/2)
		# w=wbin_le+(WBINW/2)
		q2bin_ue=self.get_q2bin_ue(q2wbin)
		q2binw=q2bin_ue-q2bin_le
		q2=q2bin_le+(q2binw/2)
		w=wbin_le+(WBINW/2)
		print "norm_1D() Going to get vgflux for %s at q2=%.2f,w=%.3f"%(q2wbin,q2,w)
		vgflux=getvgflux(w,q2,e0=self.E0)

		#! Create h1n[SEQ,VSTS,VARS]
		#! + SEQ and VSTS in this case is redundant, since normalization depends, aside from Q2,W
		#!   and vgflux, only on VAR; however this redundancy helps in code organization and readability
		h1n=OrderedDict()
		for k in h1:
			if self.DO_SS_STUDY:
				obsnum,seq,vst,var=k[0],k[1],k[2],k[3]
			else:
				seq,vst,var=k[0],k[1],k[2]
			nbins=h1[k].GetNbinsX()
			xmin=h1[k].GetXaxis().GetXmin()
			xmax=h1[k].GetXaxis().GetXmax()
			hname="hnorm_%s_VST%d_%s"%(seq,vst,var)
			h1n[k]=ROOT.TH1F(hname,hname,nbins,xmin,xmax)

			for ibin in range(nbins):
				#! Calculate var specific normalization factor
				if var=='M1' or var=='M2':
					DM=h1n[k].GetBinWidth(ibin+1)
					var_norm_fctr=DM
				elif var=='ALPHA':
					alpha_a=h1n[k].GetBinLowEdge(ibin+1)
					alpha_b=h1n[k].GetBinLowEdge(ibin+2)
					DAlpha=math.fabs(math.radians(alpha_b)-math.radians(alpha_a))
					var_norm_fctr=DAlpha
				elif var=='THETA':
					theta_a=h1n[k].GetBinLowEdge(ibin+1)
					theta_b=h1n[k].GetBinLowEdge(ibin+2)
					DCosTheta=math.fabs(math.cos(math.radians(theta_b))-math.cos(math.radians(theta_a)))
					var_norm_fctr=DCosTheta
				normf=self.LUM*LUM_INVFB_TO_INVMICROB*vgflux*q2binw*WBINW*var_norm_fctr #!Q2BINW
				h1n[k].SetBinContent(ibin+1,normf)
				h1n[k].SetBinError(ibin+1,0)

		#! Now normalize h1s
		for k in h1:
			h1[k].Divide(h1n[k])

		print "*** norm_1D() Done\n ***"

	def norm_1D_theta(self,hTheta):
		#! 1. Create normalization factor histogram
		hDCosTheta=hTheta.Clone()
		nbins=hTheta.GetNbinsX()
		for ibin in range(nbins):
			theta_a=hTheta.GetBinLowEdge(ibin+1)
			theta_b=hTheta.GetBinLowEdge(ibin+2)# + hTheta.GetBinWidth(ibin+1)
			DCosTheta=math.fabs(math.cos(math.radians(theta_b))-math.cos(math.radians(theta_a)))
			hDCosTheta.SetBinContent(ibin+1,DCosTheta)
			hDCosTheta.SetBinError(ibin+1,0.)
		#! Now divide hTheta by hDCosTheta
		#! Do Sumw2() so that errors are correctly propagated
		hTheta.Sumw2();
		hTheta.Divide(hDCosTheta)

	def plot_h2_itg_yld(self,h2):
		print "*** In DispObs::plot_h2_itg_yld() ***"

		#! Set up some plotting related aesthetics particular to this function
		#! Marker colors and styles
		#! + The idea is to use a uniqute MarkerColor for [vst,var]
		#! + As a further distinguishing feature, since colors take time to distinguish,
		#!   each var is assigned its own MarkerStyle, which can be additionally be classified
		#!   as 'Full' or 'Open' type markers, when for example, more than one seq is being plotted on the same plot.  
		#! 
		clrd={(1,'M1'):'kRed-6',     (2,'M2'):'kRed',     (3,'M2'):'kMagenta',
			  (1,'THETA'):'kBlue-6', (2,'THETA'):'kBlue', (3,'THETA'):'kCyan',
			  (1,'ALPHA'):'kGreen-6',(2,'ALPHA'):'kGreen',(3,'ALPHA'):'kOrange'}	
		mrkr_fulld={'M1':'kFullDotLarge', 'M2':'kFullDotLarge',     
					'THETA':'kFullSquare', 'ALPHA':'kFullTriangleUp'}
		mrkr_opend={'M1':'kOpenCircle', 'M2':'kOpenCircle',     
					'THETA':'kOpenSquare', 'ALPHA':'kOpenTriangleUp'}

		#! 1. Directly plot each h2 in its own canvas
		outdir_h2=os.path.join(self.OUTDIR_ITG_YLD,"Q2_W")
		if not os.path.exists(outdir_h2):
			os.makedirs(outdir_h2)
		#! Plotting aesthetics
		ROOT.gStyle.Reset() #! Reset aesthetics which may have been set by some other plot function
		ROOT.gStyle.SetOptStat("neiuo")
		ROOT.gStyle.SetErrorX(0.001)
		for k in h2:
			c=ROOT.TCanvas()
			h2[k].Draw("colz")
			c.SaveAs("%s/c_%s.png"%(outdir_h2,h2[k].GetName()))

		#! 2. Plot projection h2 on W bins i.e. int_yld(W) for various Q2 bins
		outdir_w_proj=os.path.join(self.OUTDIR_ITG_YLD,"W")
		if not os.path.exists(outdir_w_proj):
			os.makedirs(outdir_w_proj)
		nq2bins=h2[h2.keys()[0]].GetNbinsY()
		for iq2bin in range(nq2bins):
			if self.VIEW=="norm" or self.VIEW=="ERyield":
				#! Plotting aesthetics
				#ROOT.gROOT.SetOptStat("ne")
				
				#! + First create all projections: h1p[seq,vst,var]
				#! + This has to be done so that, before plotting, min and max for the Y-axis
				#!   can be determined from all the h1p[seq,vst,var]
				h1p=OrderedDict()
				for k in h2:
					seq,vst,var=k[0],k[1],k[2]
					hname="qbin_%d_%s_VST%d_%s"%(iq2bin+1,seq,vst,var)
					htitle="itg_xsec Q2=%s"%h2[k].GetYaxis().GetBinLabel(iq2bin+1)
					h1p[k]=h2[k].ProjectionX(hname,iq2bin+1,iq2bin+1,"e")
					h1p[k].SetTitle(htitle)
					if   self.VIEW=="norm": 
						h1p[k].SetYTitle("#sigma[#mub]")
						h1p[k].SetMinimum(0)
						if htitle=="itg_xsec Q2=[3.00,3.50)":
							h1p[k].SetMaximum(5.0)
						elif htitle=="itg_xsec Q2=[3.50,4.20)":
							h1p[k].SetMaximum(3.5)
						elif htitle=="itg_xsec Q2=[4.20,5.00)":
							h1p[k].SetMaximum(2.5)
						else:
							h1p[k].SetMaximum(10)#15
					elif self.VIEW=="ERyield": h1p[k].SetYTitle("counts")
					
					clr=clrd[vst,var]
					if seq=='EC':
						mrkr=mrkr_fulld[var]
					elif seq=='EF':
						mrkr=mrkr_opend[var]
					elif seq=="ER":
						mrkr=mrkr_fulld[var]
					
					h1p[k].SetMarkerStyle(ROOT.gROOT.ProcessLine(mrkr))
					h1p[k].SetLineColor(ROOT.gROOT.ProcessLine(clr))
					h1p[k].SetMarkerColor(ROOT.gROOT.ProcessLine(clr))
					
				#! Determine max and set Y-axis min(=0),max for all h1p[seq,vst,var]
				maxl=[h1p[k].GetMaximum() for k in h2]
				maximum=max(maxl)
				for k in h1p:
					h1p[k].SetMinimum(0.)
					h1p[k].SetMaximum(maximum)
				
				#! Now draw
				ROOT.gStyle.SetOptStat(0)
				c=ROOT.TCanvas("c","c",1200,800)
				#c.Divide(2,1) #! One for the plot and the other for the large legend
				pad_l=ROOT.TPad("pad_l","Title pad",0.70,0.00,1.00,1.00)
				pad_p=ROOT.TPad("pad_p","Plots pad",0.00,0.00,0.70,1.00)
				pad_l.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
				pad_p.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
				pad_l.Draw()
				pad_p.Draw()
				#! legend
				l=ROOT.TLegend(0.0,0.00,1.00,1.00)
				l.SetFillStyle(0)
				l.SetBorderSize(0)
				l.SetTextSize(0.03)#0.02
				#gpad=c.cd(1)
				#gpad.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
				pad_p.cd()
				for i,k in enumerate(h1p):
					seq,vst,var=k[0],k[1],k[2]
					draw_opt=""
					if i>0: draw_opt="same"
					h1p[k].Draw(draw_opt)
					l.AddEntry(h1p[k],"%s_VST%d_%s"%(seq,vst,var),"p")
				#gpad=c.cd(2)
				#gpad.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
				pad_l.cd()
				l.Draw()
				# #! [01-19-16] !DispObs does not know about jobtag! Add additional label for jobtag
				# pt=ROOT.TPaveText(.05,.90,.95,.95,"NDC")
				# pt.AddText("%s"%jobtag)
				# pt.Draw()
				#! Save canvas
				c.SaveAs("%s/c_qbin%d.png"%(outdir_w_proj,iq2bin+1))

			elif self.VIEW=="fullana": 
				#! Plotting aesthetics
				#ROOT.gROOT.SetOptStat("ne")
				#! + First create all projections: h1p[seq,vst,var]
				h1p=OrderedDict()
				for k in h2:
					seq,vst,var=k[0],k[1],k[2]
					if seq=='SH' or seq=='EH': continue
					hname="qbin_%d_%s_VST%d_%s"%(iq2bin+1,seq,vst,var)
					htitle="itg_yld Q2=%s"%h2[k].GetYaxis().GetBinLabel(iq2bin+1)
					h1p[k]=h2[k].ProjectionX(hname,iq2bin+1,iq2bin+1,"e")
					h1p[k].SetTitle(htitle)
					
					clr=clrd[vst,var]
					#! Set appropriate MarkerStyles for SEQS plotted on same cavas
					#! ER,SR on same
					if seq=='ER':
						mrkr=mrkr_opend[var]
					elif seq=='SR':
						mrkr=mrkr_fulld[var]
					#! EC,EF on same; SC,SF,ST on same	
					if seq=='EC' or seq=='SC':
						mrkr=mrkr_fulld[var]
					elif seq=='SF' or seq=='EF':
						mrkr=mrkr_opend[var]
					elif seq=='ST':
						mrkr="kPlus"
					
					h1p[k].SetMarkerStyle(ROOT.gROOT.ProcessLine(mrkr))
					h1p[k].SetLineColor(ROOT.gROOT.ProcessLine(clr))
					h1p[k].SetMarkerColor(ROOT.gROOT.ProcessLine(clr))
				
				#! Now draw
				#! Plotting aesthetics
				ROOT.gStyle.SetOptStat(0)
				#! Draw ER,SR normalized to same integral
				c1=ROOT.TCanvas("c1","c1",1000,500)#,1500,2000)
				c1.Divide(2,1)
				h1n=OrderedDict() #! To keep normed hists
				gpad=c1.cd(1)
				gpad.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
				l1=ROOT.TLegend(0.0,0.00,1.0,1.00)
				l1.SetFillStyle(0)
				l1.SetBorderSize(0)
				l1.SetTextSize(0.02)
				i=0
				for k in h1p:
					seq,vst,var=k[0],k[1],k[2]
					if seq!='ER' and seq!='SR': continue
					draw_opt="" if i==0 else "same"
					h1n[k]=h1p[k].DrawNormalized(draw_opt,1000)
					l1.AddEntry(h1n[k],"%s_VST%d_%s"%(k[0],k[1],k[2]),"p")
					i+=1
				#! Add additional label for ER,SR to say that they are normed to same integral
				pt=ROOT.TPaveText(.05,.90,.95,.95,"NDC")
				pt.AddText("Integral of all hists. normalized to 1000")
				pt.Draw()
				gpad=c1.cd(2)
				gpad.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
				l1.Draw()
				c1.SaveAs("%s/c_qbin%d_ER_SR.png"%(outdir_w_proj,iq2bin+1))

				#! Draw ST,SC,SF
				#! First set max of Y-axis
				maxl=[h1p[k].GetMaximum() for k in h1p if k[0]=='ST' or k[0]=='SC' or k[0]=='SF']
				maximum=max(maxl)
				for k in h1p:
					if k[0]=='ST' or k[0]=='SC' or k[0]=='SF':
						h1p[k].SetMinimum(0)
						h1p[k].SetMaximum(maximum)

				c2=ROOT.TCanvas("c2","c2",1000,500)
				c2.Divide(2,1)
				gpad=c2.cd(1)
				gpad.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
				l2=ROOT.TLegend(0.0,0.00,1.0,1.00)
				l2.SetFillStyle(0)
				l2.SetBorderSize(0)
				l2.SetTextSize(0.02)
				i=0
				#! First draw SC,SF
				for k in h1p:
					seq,vst,var=k[0],k[1],k[2]
					if seq!='SC' and seq != 'SF': continue
					draw_opt="" if i==0 else "same"
					h1p[k].Draw(draw_opt)
					l2.AddEntry(h1p[k],"%s_VST%d_%s"%(k[0],k[1],k[2]),"p")
					i+=1
				#! Now draw ST so that it appears on top
				for k in h1p:
					seq,vst,var=k[0],k[1],k[2]
					if seq!='ST': continue
					draw_opt="" if i==0 else "same"
					h1p[k].Draw(draw_opt)
					l2.AddEntry(h1p[k],"%s_VST%d_%s"%(k[0],k[1],k[2]),"p")
					i+=1
				gpad=c2.cd(2)
				gpad.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
				l2.Draw()
				c2.SaveAs("%s/c_qbin%d_ST_SC_SF.png"%(outdir_w_proj,iq2bin+1))

				#! Draw EC,EF
				#! First set max of Y-axis
				maxl=[h1p[k].GetMaximum() for k in h1p if k[0]=='EC' or k[0]=='EF']
				maximum=max(maxl)
				for k in h1p:
					if k[0]=='EC' or k[0]=='EF':
						h1p[k].SetMinimum(0)
						h1p[k].SetMaximum(maximum)

				c3=ROOT.TCanvas("c3","c3",1000,500)
				c3.Divide(2,1)
				gpad=c3.cd(1)
				gpad.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
				l3=ROOT.TLegend(0.0,0.00,1.0,1.00)
				l3.SetFillStyle(0)
				l3.SetBorderSize(0)
				l3.SetTextSize(0.02)
				i=0
				for k in h1p:
					seq,vst,var=k[0],k[1],k[2]
					if seq!='EC' and seq != 'EF': continue
					draw_opt="" if i==0 else "same"
					h1p[k].Draw(draw_opt)
					l3.AddEntry(h1p[k],"%s_VST%d_%s"%(k[0],k[1],k[2]),"p")
					i+=1
				gpad=c3.cd(2)
				gpad.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
				l3.Draw()
				c3.SaveAs("%s/c_qbin%d_EC_EH.png"%(outdir_w_proj,iq2bin+1))

		print "*** plot_h2_itg_yld() Done ***\n"

	def plot_h2_itg_yld_SS(self,h2):
		'''
		[02-29-16]
		+ 'plot_h2_itg_yld' was used as a tempate for this method
		+ The idea is extend the method to incorporate plotting of a 
		  list of SS-Obs and the associated SS-err
		+ Therefore only details relevant to the above were changed and the core functionality left intact

		+ Process:
			1. First h2 is directly drawn and saved.
			2. Then for every Q2 bin in h2:
				2.i.   Projection of h2 on W ais made, which is the itg_xsec(W;Q2)
				2.ii.  SS err for itg_yld(W;Q2) is obtained (SS_err_itg_xsec(W;Q2))
				2.iii. For every for every [EC,EF]*[vst]*[var] plot and save both 
				       itg_yld(W;Q2) and SS_err_itg_yld(W;Q2).
		'''

		print "*** In DispObs::plot_h2_itg_yld_SS() ***"

		#! 1. Directly plot each h2[obsnum,seq,vst,var] in its own canvas
		outdir_h2=os.path.join(self.OUTDIR_ITG_YLD,"Q2_W")
		if not os.path.exists(outdir_h2):
			os.makedirs(outdir_h2)
		#! Plotting aesthetics
		ROOT.gStyle.Reset() #! Reset aesthetics which may have been set by some other plot function
		ROOT.gStyle.SetOptStat("neiuo")
		ROOT.gStyle.SetErrorX(0.001)
		for k in h2:
			c=ROOT.TCanvas()
			h2[k].Draw("colz")
			c.SaveAs("%s/c_%s.png"%(outdir_h2,h2[k].GetName()))

		#! 2. Plot projections of h2 on W bins and SS error on this projections for various Q2 bins
		#!    i.e. intg_xsec(W) and SS_err_intg_yld(W) for various Q2 bins
		# outdir_w_proj=os.path.join(self.OUTDIR_ITG_YLD,"W")
		# outdir_w_proj_SS_err=os.path.join(self.OUTDIR_ITG_YLD_SS_ERR,"W")
		# if not os.path.exists(outdir_w_proj):
		# 	os.makedirs(outdir_w_proj)
		# if not os.path.exists(outdir_w_proj_SS_err):
		# 	os.makedirs(outdir_w_proj_SS_err)
		nq2bins=h2[h2.keys()[0]].GetNbinsY()
		for iq2bin in range(nq2bins):
			#! 2.i. First create all projections: h1p (itg_xsec(W;Q2))
			#! + This has to be done so that, before plotting, min and max for the Y-axis
			#!   can be determined from all the h1p[seq,vst,var]
			h1p=OrderedDict()
			for k in h2:
				obsnum,seq,vst,var=k[0],k[1],k[2],k[3]
				hname="qbin_%d_%d_%s_VST%d_%s"%(iq2bin+1,obsnum,seq,vst,var)
				htitle="intg xsec for Q2=%s from %s_VST%d_%s"%(h2[k].GetYaxis().GetBinLabel(iq2bin+1),seq,vst,var)
				h1p[k]=h2[k].ProjectionX(hname,iq2bin+1,iq2bin+1,"e")
				h1p[k].SetTitle(htitle)
				h1p[k].SetYTitle("#sigma[#mub]")
				h1p[k].SetMinimum(0)
				if   "Q2=[3.00,3.50)" in htitle: h1p[k].SetMaximum(5.0)
				elif "Q2=[3.50,4.20)" in htitle:	h1p[k].SetMaximum(3.5)
				elif "Q2=[4.20,5.00)" in htitle:	h1p[k].SetMaximum(2.5)
				else:                                   h1p[k].SetMaximum(10)#15
									
				mrk=MRKS_CLRS_N_COORDS_DICT[obsnum][0]
				clr=MRKS_CLRS_N_COORDS_DICT[obsnum][1]
				h1p[k].SetMarkerStyle(mrk)
				h1p[k].SetMarkerColor(clr)
				h1p[k].SetLineColor(clr)
							
			#! Determine max and set Y-axis min(=0),max for all h1p[seq,vst,var]
			maxl=[h1p[k].GetMaximum() for k in h2]
			maximum=max(maxl)
			for k in h1p:
				h1p[k].SetMinimum(0.)
				h1p[k].SetMaximum(maximum)

			#! 2.ii. Get SS err for the projected histograms (SS_err_itg_xsec(W;Q2))
			h1p_SS_err=self.get_SS_err_1D(h1p)
			#! fix labels for this histogram
			for k in h1p_SS_err:
				htitle="Systematic error for %s"%h1p_SS_err[k].GetTitle()
				h1p_SS_err[k].SetTitle(htitle)
				h1p_SS_err[k].SetYTitle("Systematic error [%] for itg. xsec")

			#! 2.iii. For every [EC,EF]*[vst]*[var] plot itg_xsec(W;Q2) and SS_err_itg_xsec(W;Q2) 
			#! + Construct output domain(d)=EC,EF]*[vst]*[var]
			#! + First construct vst_var to choose the 9 1D-Obs
			vst_var=[[1,'M1'],[1,'THETA'],[1,'ALPHA'],
			         [2,'M2'],[2,'THETA'],[2,'ALPHA'],
			         [3,'M2'],[3,'THETA'],[3,'ALPHA']]
			#! + Now take the outer product with [EC,EF] 
			d=list(itertools.product(self.SEQS,vst_var))
			#! Now for each item in the output domain, save plots
			for r in d:
				seq,vst,var=r[0],r[1][0],r[1][1]
				#! prepare outdir for itg_xsec and SS_err_itg_xsec
				outdir_w_proj=       os.path.join(self.OUTDIR_ITG_YLD,       seq,"%d_%s"%(vst,var),"W")
				outdir_w_proj_SS_err=os.path.join(self.OUTDIR_ITG_YLD_SS_ERR,seq,"%d_%s"%(vst,var),"W")
				if not os.path.exists(outdir_w_proj):
					os.makedirs(outdir_w_proj)
				if not os.path.exists(outdir_w_proj_SS_err):
					os.makedirs(outdir_w_proj_SS_err)

				#! Firt plot and save itg_xsec 
				ROOT.gStyle.SetOptStat(0)
				c=ROOT.TCanvas("c","c",1200,800)
				c.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
				#! legend
				#! legend removed after adding legend canvas in plot_1D_SS()
				# l=ROOT.TLegend(0.60,0.60,1.00,0.90)
				# l.SetFillStyle(0)
				# l.SetBorderSize(0)
				# l.SetTextSize(0.03)#0.02
				i=0
				for obsnum in self.OBSD:
					obstag=self.OBSD[obsnum][1]
					draw_opt=""
					if i>0: draw_opt="same"
					h1p[obsnum,seq,vst,var].Draw(draw_opt)
					i+=1
					#l.AddEntry(h1p[obsnum,seq,vst,var],"%s_%s_VST%d_%s"%(obstag,seq,vst,var),"p")
				#l.Draw()
				c.SaveAs("%s/c_qbin%d.png"%(outdir_w_proj,iq2bin+1))

				#! Now plot and save SS_err_itg_xsec
				ROOT.gStyle.SetOptStat(0)
				cerr=ROOT.TCanvas("cerr","cerr",1200,800)
				h1p_SS_err[seq,vst,var].SetMinimum(0)
				h1p_SS_err[seq,vst,var].SetMaximum(20)
				h1p_SS_err[seq,vst,var].Draw()
				cerr.SaveAs("%s/c_qbin%d.png"%(outdir_w_proj_SS_err,iq2bin+1))

		print "*** plot_h2_itg_yld_SS() Done ***\n"
		
		
	def fill_h2_itg_yld(self,h2,h1,iq2bin,iwbin):
		print "*** In DispObs::fill_h2_itg_yld ***"
		
		for k in h1:
			if self.DO_SS_STUDY:
				obsnum,seq,vst,var=k[0],k[1],k[2],k[3]
			else:
				seq,vst,var=k[0],k[1],k[2]
			if var=='M1' or var=='M2': #! then directly use TH1::Integral(opt); opt="width" when self.VIEW="norm", "" when self.VIEW="fullana"/"ERyield"
				if self.DO_SS_STUDY: opt="width"
				else:                opt="%s"%("width" if self.VIEW=="norm" else "")
				nbins=h1[k].GetNbinsX()
				itg_err=ROOT.Double(0)#https://root.cern.ch/phpBB3/viewtopic.php?t=2499
				itg=h1[k].IntegralAndError(1,nbins,itg_err,opt)
				h2[k].SetBinContent(iwbin+1,iq2bin+1,itg)
				h2[k].SetBinError(iwbin+1,iq2bin+1,itg_err)
			elif var=='THETA': #! manually compute integral to adapt for DCosTheta irrespective of self.VIEW (NOTE, theta distributions are atleast DCosTheta normalized irrespective of self.VIEW)
				itg,itg_err=0,0
				nbins=h1[k].GetNbinsX()
				for ibin in range(nbins):
					binc=h1[k].GetBinContent(ibin+1)
					bincerr=h1[k].GetBinError(ibin+1)
					theta_a=h1[k].GetBinLowEdge(ibin+1)
					theta_b=h1[k].GetBinLowEdge(ibin+2)
					DCosTheta=math.fabs(math.cos(math.radians(theta_b))-math.cos(math.radians(theta_a)))
					itg+=binc*DCosTheta
					itg_err+=bincerr*DCosTheta
				h2[k].SetBinContent(iwbin+1,iq2bin+1,itg)
				h2[k].SetBinError(iwbin+1,iq2bin+1,itg_err)
			elif var=='ALPHA': #! manually compute integral to adapt for DAlpha; DAlpha is in radians when self.VIEW="norm", 1 when self.VIEW="fullana"/"ERyield"
				itg,itg_err=0,0
				nbins=h1[k].GetNbinsX()
				for ibin in range(nbins):
					binc=h1[k].GetBinContent(ibin+1)
					bincerr=h1[k].GetBinError(ibin+1)
					alpha_a=h1[k].GetBinLowEdge(ibin+1)
					alpha_b=h1[k].GetBinLowEdge(ibin+2)# + h1[k].GetBinWidth(ibin+1)
					if self.DO_SS_STUDY:
						DAlpha=math.fabs(math.radians(alpha_b)-math.radians(alpha_a))
					else:
						if self.VIEW=="norm":
							DAlpha=math.fabs(math.radians(alpha_b)-math.radians(alpha_a))
						else:
							DAlpha=1
					itg+=binc*DAlpha
					itg_err+=bincerr*DAlpha
				h2[k].SetBinContent(iwbin+1,iq2bin+1,itg)
				h2[k].SetBinError(iwbin+1,iq2bin+1,itg_err)


		print "*** fill_h2_itg_yld() Done ***\n"	

	def make_h2_itg_yld(self,q2bng,wbng,h1):
		'''
		+ h1 needed to provide keys: either [seq,vst,var] or if self.DO_SS_STUDY== true, [obsnum,seq,vst,var]
		+ Make h2_itg_yld[seq,vst,var] with xaxis,yaxis=nwbins,nq2bins
			+ 'var' is included for debugging purposes since integrated yield 
			  from all VSTS*VARS should be the same.
		'''
		#! First get binning in Q2 and W from Q2-W binning
		print "In DispObs::make_h2_itg_yld()"
		nq2bins=q2bng['NBINS']
		nwbins=wbng['NBINS']
		q2bins=array.array('f',range(1,nq2bins+1+1))#!Last '+1' to ensure 'range' includes top edge of last bin 
		wbins=array.array('f',range(1,nwbins+1+1))#!Last '+1' to ensure 'range' includes top edge of last bin
		h2=OrderedDict()
		for k in h1:
			if self.DO_SS_STUDY:
				obsnum,seq,vst,var=k[0],k[1],k[2],k[3]
				obstag=self.OBSD[obsnum][1]
				name="h2_%s_%s_VST%d_%s"%(obstag,seq,vst,var)
			else:
				seq,vst,var=k[0],k[1],k[2]
				name="h2_%s_VST%d_%s"%(seq,vst,var)
			
			title="itg_yld(Q^{2},W)"
			h2[k]=ROOT.TH2F(name,title,nwbins,wbins,nq2bins,q2bins) 
			h2[k].SetXTitle("W (GeV)")
			h2[k].SetYTitle("Q^{2} (GeV^{2})")
			#! Label the bins
			for iwbin in range(nwbins):
				le=wbng['BINS'][iwbin]
				ue=wbng['BINS'][iwbin+1]
				h2[k].GetXaxis().SetBinLabel(iwbin+1,"[%.3f,%.3f)"%(le,ue))
			for iq2bin in range(nq2bins):
				le=q2bng['BINS'][iq2bin]
				ue=q2bng['BINS'][iq2bin+1]
				h2[k].GetYaxis().SetBinLabel(iq2bin+1,"[%.2f,%.2f)"%(le,ue))
			
		#! For debuggin print out h2 information for first key in h2
		print "make_h2_itg_yld() nq2bins,nwbins=%d,%d"%(h2[h2.keys()[0]].GetNbinsY(),h2[h2.keys()[0]].GetNbinsX())
		print "make_h2_itg_yld() wbin labels:" 
		for iwbin in range(nwbins):
			print "wbin#=%d,label=%s"%(iwbin+1,h2[h2.keys()[0]].GetXaxis().GetBinLabel(iwbin+1))
		print "make_h2_itg_yld() q2bin labels:" 
		for iq2bin in range(nq2bins):	
			print "q2bin#=%d,label=%s"%(iq2bin+1,h2[h2.keys()[0]].GetYaxis().GetBinLabel(iq2bin+1))

		print "***make_h2_itg_yld() Done \n***"
		return h2

	def disp_itg_yld_drct_from_h5(self):
		"""
		Walk the ROOT file and plot y(w;seq,vst,q2bin). 
		"""
		print "*** In DispObs::disp_itg_yld_drct_from_h5() ***"
		
		# if self.VIEW==False:
		# 	outdir=os.path.join(self.OUTDIR,"Obs_Itg_Yld_drct")
		# else:
		# 	outdir=os.path.join(self.OUTDIR,"Obs_Itg_Yld_drct_norm")
		outdir=os.path.join(self.OUTDIR,"Obs_Itg_Yld_drct_%s"%self.VIEW)
		if not os.path.exists(outdir):
			os.makedirs(outdir)

		#! 1. Get all q2wbins
		if self.DBG==True:
			q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX,dbg=True,dbg_bins=2)
		else:
			q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX)
		#print q2wbinl

		q2bng=self.get_q2bng(q2wbinl)
		#print q2bng['BINW']
		print "Yield(w) will be plotted for the following Q2 bins:"
		print q2bng['BINS']

		#! 3. Now put together dictionary for yield: y{seq,vst,q2bin:{w:yield}}
		#! First create dict
		y={}
		oy={} #Will be used later to store yields ordered
		for seq in self.SEQS:
			for vst in VSTS:
				for q2wbin in q2wbinl:
					q2bin=q2wbin.split('_')[0]
					y[seq,vst,q2bin]={}
					oy[seq,vst,q2bin]={}
		#! Now fill dict
		for seq in self.SEQS:
			for vst in VSTS:
				for q2wbin in q2wbinl:
					print "seq,vst,q2wbin",seq,vst,q2wbin
					q2bin=q2wbin.split('_')[0]
					#! Get w, yield
					w=float(q2wbin.split('_')[1].split('-')[0])
					dw=float(q2wbin.split('_')[1].split('-')[1])-w
					#wbin=q2wbin.split('_')[1]
					h5=self.FIN.Get("%s/%s/VST%d/h5"%(q2wbin,seq,vst))
					y[seq,vst,q2bin][w]=thntool.GetIntegral(h5)
					if self.VIEW=="norm":
						q2=float(q2wbin.split('_')[0].split('-')[0])
						dq2=float(q2wbin.split('_')[0].split('-')[1])-q2
						#normf=self.LUM*LUM_INVFB_TO_INVMICROB*getvgflux(w,q2,e0=self.E0)*dw*dq2#!mub^-1
						normf=self.LUM*LUM_INVFB_TO_INVMICROB*getvgflux(w+(dw/2),q2+(dq2/2),e0=self.E0)*dw*dq2#!mub^-1
						#normf=self.LUM*LUM_INVFB_TO_INVMICROB*getvgflux(w+dw,q2+dq2,e0=self.E0)*dw*dq2#!mub^-1
						print "yield=",y[seq,vst,q2bin][w]
						print "dw,dq2=",dw,dq2
						print "norm=",normf
						y[seq,vst,q2bin][w]=y[seq,vst,q2bin][w]/normf
						print "yield after norm=",y[seq,vst,q2bin][w]
		#! Make sure y[seq,q2bin:w] are sorted by w
		for k in y.keys():
			oy[k]=OrderedDict(sorted(y[k].items()))
			
		#! 4. Now plot
		fig=plt.figure()
		ax=plt.subplot(111)
		# clrd={('1.25-1.75',1):'red',('1.25-1.75',2):'brown',('1.25-1.75',3):'magenta',
		# 	  ('1.75-2.25',1):'cyan',('1.75-2.25',2):'blue',('1.75-2.25',3):'green',
		# 	  ('2.25-2.75',1):'orange',('2.25-2.75',2):'yellow',('2.25-2.75',3):'black',
		# 	  ('2.75-3.25',1):'cyan',('2.75-3.25',2):'blue',('2.75-3.25',3):'green',
		# 	  ('3.25-3.75',1):'cyan',('3.25-3.75',2):'blue',('3.25-3.75',3):'green',
		# 	  ('3.75-4.25',1):'cyan',('3.75-4.25',2):'blue',('3.75-4.25',3):'green',
		# 	  ('4.25-4.75',1):'cyan',('4.25-4.75',2):'blue',('4.25-4.75',3):'green',
		# 	  ('4.75-5.25',1):'orange',('4.75-5.25',2):'yellow',('4.75-5.25',3):'black'}
		clrd={('1.50-2.00',1):'red',('1.50-2.00',2):'brown',('1.50-2.00',3):'magenta',
			  ('2.00-2.40',1):'cyan',('2.00-2.40',2):'blue',('2.00-2.40',3):'green',
			  ('2.40-3.00',1):'orange',('2.40-3.00',2):'yellow',('2.40-3.00',3):'black',
			  ('3.00-3.50',1):'cyan',('3.00-3.50',2):'blue',('3.00-3.50',3):'green',
			  ('3.50-4.20',1):'cyan',('3.50-4.20',2):'blue',('3.50-4.20',3):'green',
			  ('4.20-5.00',1):'cyan',('4.20-5.00',2):'blue',('4.20-5.00',3):'green'}
		mrkrd={'EC':'o','EF':'^'}
		for k in oy.keys():
			seq,vst,q2wbin=k[0],k[1],k[2]
			lbl='%s:VST%d[%s)'%(seq,vst,q2wbin)
			ax.scatter(oy[k].keys(),oy[k].values(),label=lbl,c=clrd[q2wbin,vst],marker=mrkrd[seq],s=50)
			ax.set_xticks(oy[k].keys())
			ax.set_xticklabels(oy[k].keys(),rotation='vertical')
		if self.VIEW!="norm":
			ax.set_ylim(0,600000)
			ax.set_ylabel(r'Yield [A.U.]',fontsize='xx-large')
		else:
			ax.grid(1)
			ax.set_ylim(0,10)
			ax.set_ylabel(r'$\mu b$',fontsize='xx-large')
		ax.legend(loc='upper left',prop={'size':8}) #loc='lower right'
		fig.savefig('%s/integ_yield.png'%(outdir))
		#fig.savefig('%s/integ_yield.eps'%(outdir))

		print "*** disp_itg_yld_old_mchnry() Done ***\n"

	
	def get_sim_stats(self):
		"""
		Walk the ROOT file and obtain simstats(ss) for a h5 in a Q2-W bin:
		ss={'ST':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]],
			'SR':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]],
			'SA':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]],
			'SH':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]]}

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

		ss={'ST':[],'SR':[],'SA':[],'SH':[]}
		for seq in ['ST','SR','SA','SH']:
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
				h5=self.FIN.Get("%s/%s/VST1/h5"%(q2w,seq))
				nbins=thntool.GetNbinsNotEq0(h5)
				N=thntool.GetIntegral(h5)
				binc_stats=np.zeros(2,'f')
				thntool.GetBinContentDistStats(h5,binc_stats)
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

	def plot_sim_stats_q2wbin_5D(self,q2wbin='',acc_rel_err_cut=-1,cut_type='lt',vst=1):
		"""
		[11-16-15] 
		+ Orignal in study_acc_tools_lite.py::plot_acc_q2wbin_5D()
		+ Brought over to disp_obs.py in order to consolidate

		Notes:
		+ Unless q2wbin is specified, this function will plot and save 5D simstats 
          for all q2wbins in yield.root
        + If q2wbin is specified, then this function will make an interactive plot for that
          particular q2wbin
		"""
		is_intrctv=False
		if q2wbin=='':#! Get all q2wbin directories from file
			q2wbinl=self.get_q2wbinlist()
			#! Create output directory
			outdir=os.path.join(self.OUTDIR,'simstats')
			if not os.path.exists(outdir):
				os.makedirs(outdir)
			#! Create output.root file
			sfx=''
			if acc_rel_err_cut>0:sfx="_acc_relerr_%s_%.2f"%(cut_type,acc_rel_err_cut)
			fout=ROOT.TFile("%s/simstats_q2wbin%s.root"%(outdir,sfx),"RECREATE")
		else:#! interactive mode for particular q2wbin
			q2wbinl=[q2wbin]
			is_intrctv=True;

		print "Going to plot 5D simstats for q2wbinl=",q2wbinl

		ROOT.gStyle.SetOptStat("nemMrRuo")

		for q2wbin in q2wbinl:
			#! Get 5D histograms
			h5T=self.FIN.Get('%s/ST/VST%d/h5'%(q2wbin,vst))
			h5R=self.FIN.Get('%s/SR/VST%d/h5'%(q2wbin,vst))
			h5A=self.FIN.Get('%s/SA/VST%d/h5'%(q2wbin,vst))
	
			nbins=h5T.GetNbins()
			hT=ROOT.TH1F("hT","ST:%s"%q2wbin,400,0,400)
			hR=ROOT.TH1F("hR","SR:%s"%q2wbin,60,0,60)
			hA=ROOT.TH1F("hA","SA:%s"%q2wbin,100,0,1.5)
			h_relErrA=ROOT.TH1F("hA_rel_err","SA_rel_err:%s"%q2wbin,100,0,1.5)
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
					if acc_rel_err_cut>0: #! i.e. make a cut on acc_rel_err_cut
						pass_cut=False
						if   cut_type=='lt': pass_cut=rel_err_A<acc_rel_err_cut
						elif cut_type=='gt': pass_cut=rel_err_A>acc_rel_err_cut
						if pass_cut: 
							hT.Fill(bincT)
							hR.Fill(bincR)
							hA.Fill(bincA)
							h_relErrA.Fill(rel_err_A)
					else:
						hT.Fill(bincT)
						hR.Fill(bincR)
						hA.Fill(bincA)
						h_relErrA.Fill(rel_err_A)
	
				# if bincA>0:# and bincA<1:
				# 	hT.Fill(bincT)
				# 	hR.Fill(bincR)
				# 	hA.Fill(bincA)
				# 	h_relErrA.Fill(rel_err_A)
		
			#! Now draw hists
			cname="c_%s"%q2wbin
			c=ROOT.TCanvas(cname,cname)
			c.Divide(2,2)
			c.cd(1)
			hT.Draw()
			c.cd(2)
			hR.Draw()
			c.cd(3)
			hA.Draw()
			c.cd(4)
			h_relErrA.Draw()

			if is_intrctv:
				plt.show()
				# wait for you to close the ROOT canvas before exiting
				wait(True)
			else:
				c.Write()


	def get_q2wbinlist(self,q2min=0.00,q2max=6.00,wmin=0.000,wmax=3.000,dbg=False,dbg_bins=10):
		"""
		The ROOT file is arranged in a Tree Structure. The
		Observable histograms (obs-hists) are located as files in the following directory-path(dirpath):
		q2wbin/vst/seq/hists
		"""
		q2wbinl=[]
		#! File to use to extract q2wbinl
		if self.DO_SS_STUDY: #! Use obsdir of the first item in OBSD
			obsdir=os.path.join(self.OBSDIR,self.OBSD[1][0])
			#print "obsdir=",obsdir
			f=root_open(os.path.join(obsdir,'yield.root'))
			#print f.GetName()
		else: 
			f=self.FIN

		print "DispObs::get_q2wbinlist() Going to Q2-W bins from file=",f.GetName()
		print "DispObs::get_q2wbinlist() q2min,q2max,wmin,wmax=",q2min,q2max,wmin,wmax
		if dbg==True:
			print "DispObs::get_q2wbinlist() dbg=True"

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
	def get_q2bin_ue(self,q2wbin):
		return float(self.get_q2bin(q2wbin).split('-')[1])
	def get_wbin_ue(self,q2wbin):
		return float(self.get_wbin(q2wbin).split('-')[1])
	
	def hist_1D_athtcs(self,h1):
		'''
		[03-02-16]
		After addition of self.DO_SS_STUDY, this method was modified so that it can understand 
		the type of h1 being passed, which could be either of the following depending on the condition (noted
		in parenthesis and programmed into the logic of this function):
		+ h1(seq,vst,var)        (nkeys=3)
		+ h1(obsnum,seq,vst,var) (nkeys=4)
		+ h1_SS_err(seq,vst,var) (nkeys=3 and self.DO_SS_STUDY=true)
		'''
		#! First determine if h1 has 3 or 4 keys(3 keys=seq,vst,var; 4 keys=obsnum,seq,vst,var)
		#! + This structure of h1 depends on self.DO_SS_STUDY
		nkeys=len(h1.keys()[0])
		print("DispObs::hist_1D_athtcs(): nkeys=%d"%nkeys)
		for k in h1.keys():	
			if   nkeys==4: obsnum,seq,vst,var=k[0],k[1],k[2],k[3]
			elif nkeys==3: seq,vst,var=k[0],k[1],k[2]
			else: 		   sys.exit("DispObs::hist_1D_athtcs(): nkeys undetermined. Exiting.")

			h1[k].SetTitle("")
			h1[k].SetXTitle( "%s[%s]"%(VAR_NAMES[(vst,var)],VAR_UNIT_NAMES[var]) )
			#! X-axis title aesthetics
			hR2d[k].GetXaxis().SetTitleOffset(0.7)
			hR2d[k].GetXaxis().SetLabelSize(0.05)
			hR2d[k].GetXaxis().SetTitleSize(0.10)
			#! Y-axis title aesthetics
			hR2d[k].GetYaxis().SetTitleOffset(0.7)
			hR2d[k].GetYaxis().SetLabelSize(0.05)
			hR2d[k].GetYaxis().SetTitleSize(0.10)

			if nkeys==4:#s i.e. self.DO_SS_STUDY=true:
				if var!='THETA':
					h1[k].SetYTitle("#Delta#sigma/#Delta%s [#mub/%s]"%(VAR_NAMES[(vst,var)],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
				elif var=='THETA': #! only difference is to add the 'Cos' of DCosTheta
					h1[k].SetYTitle("#Delta#sigma/#Deltacos(%s) [#mub/%s]"%(VAR_NAMES[(vst,var)],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))

				mrk=MRKS_CLRS_N_COORDS_DICT[obsnum][0]
				clr=MRKS_CLRS_N_COORDS_DICT[obsnum][1]
				h1[k].SetMarkerStyle(mrk)
				h1[k].SetMarkerColor(clr)
				h1[k].SetLineColor(clr)
			elif nkeys==3 and self.DO_SS_STUDY: #! i.e. self.DO_SS_STUDY=true and it is h1_SS_err
				if var!='THETA':
					h1[k].SetYTitle("Systematic error [%%] for #Delta#sigma/#Delta%s"%(VAR_NAMES[(vst,var)]))
				elif var=='THETA': #! only difference is to add the 'Cos' of DCosTheta
					h1[k].SetYTitle("Systematic error [%%] for #Delta#sigma/#Deltacos(%s)"%(VAR_NAMES[(vst,var)]))
				
				h1[k].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
				h1[k].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
				h1[k].SetMarkerStyle(ROOT.gROOT.ProcessLine("kFullDotLarge"))
			elif nkeys==3:
				if self.VIEW=="norm":
					if var!='THETA':
						h1[k].SetYTitle("#Delta#sigma/#Delta%s [#mub/%s]"%(VAR_NAMES[(vst,var)],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
					elif var=='THETA': #! only difference is to add the 'Cos' of DCosTheta
						h1[k].SetYTitle("#Delta#sigma/#Deltacos(%s) [#mub/%s]"%(VAR_NAMES[(vst,var)],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
				elif self.VIEW=="fullana":
					h1[k].SetYTitle("\"Normalized\" counts")
				elif self.VIEW=="ERyield":
					h1[k].SetYTitle("Counts")

				if self.VIEW!="ERyield":
					h1[k].SetMarkerColor(CLRS_PLT_SEQ_ALL[seq])
					h1[k].SetLineColor(CLRS_PLT_SEQ_ALL[seq])
					h1[k].SetMarkerStyle(MRKS_PLT_SEQ_ALL[seq])

	def hist_R2_athtcs(self,hR2d,R2):
		print("DispObs::hist_R2_athtcs()")
		for k in hR2d:
			seq,vst,var=k[0],k[1],k[2]
			hR2d[k].SetMarkerStyle(MRKS_PLT_SEQ_ALL[seq]) #! ROOT.gROOT.ProcessLine("kFullCircle")
			hR2d[k].SetMarkerColor(CLRS_PLT_SEQ_ALL[seq])
			hR2d[k].SetLineColor(CLRS_PLT_SEQ_ALL[seq])
			hR2d[k].SetTitle("")
			hR2d[k].SetXTitle( "%s%s"%(VAR_NAMES[(vst,var)],VAR_UNIT_NAMES[var]) )
			#! post VM-fdbk
			#hR2d[k].SetYTitle("%s [#mub/%s]"%(R2_NAMED[R2],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
			#! tchncl problem with following
			#hR2d[k].SetYTitle("%s^{%s}_{%s} [#mub/%s]"%(R2_NAMED[R2],VAR_NAMES[(vst,var)],VAR_NAMES[(vst,'PHI')],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
			#! Solution 2: _R2 ^{X_{ij}}_{phi_{i}}
			hR2d[k].SetYTitle("%s ^{%s}_{%s} [#mub/%s]"%(R2_NAMED[R2],VAR_NAMES[(vst,var)],VAR_NAMES[(vst,'PHI')],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
			#! Solution 2: _{phi_{i}}R2^{X_{ij}}
			#hR2d[seq,vst,var].SetYTitle("_{%s}%s^{%s}_{%s} [#mub/%s]"%(VAR_NAMES[(vst,'PHI')],R2_NAMED[R2],VAR_NAMES[(vst,var)],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))

			#! X-axis title aesthetics
			hR2d[k].GetXaxis().SetTitleOffset(0.7)
			hR2d[k].GetXaxis().SetLabelSize(0.05)
			hR2d[k].GetXaxis().SetTitleSize(0.10)
			#! Y-axis title aesthetics
			hR2d[k].GetYaxis().SetTitleOffset(0.7)
			hR2d[k].GetYaxis().SetLabelSize(0.05)
			hR2d[k].GetYaxis().SetTitleSize(0.10)
		print("DispObs::hist_R2_athtcs()")


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

	def disp_R2(self,R2l,mthd,extract_D_E_for_non_alpha=True,plotphiproj=False,fopt=""):
		"""
		+ Extract user specified R2(s) in R2l using user specified method(='mthd1'/'mthd2'/'mthd3' ='h5-mply-itg'/'phi-proj-fit'/'phi-proj-mply-itg')
									
		+ The following is the outline of the process:
			1. Get Q2-W bin list from yield.root
			2. Loop over Q2-W bin list and
				2.1. Get h5d[self.SEQS,VST]
				2.2. From h5d[self.SEQS,VST], extract hR2d[self.SEQS,vst,var] as per specified method
				2.3. Plot hR2d[self.SEQS,vst,var]
		"""
		print "In DispObs::disp_R2()"

		#! Setup Class using variables passed in through this function
		#! + [03-18-16] I discovered only today that Class variables could be initialized from a 
		#!   member function and therefore of all such initializations I could do here, I have only
		#!   setup self.FOPT for now, which is a new argument I added for disp_R2() today.           
		self.FOPT=fopt
		print "self.FOPT=",self.FOPT
		self.EXTRACT_D_E_FOR_NON_ALPHA=extract_D_E_for_non_alpha
		print "self.EXTRACT_D_E_FOR_NON_ALPHA=",self.EXTRACT_D_E_FOR_NON_ALPHA

		self.OUTDIR_OBS_R2=os.path.join(self.OUTDIR,"Obs_R2_%s"%self.VIEW,"mthd_%s"%MTHD_NAMED[mthd])
		#! Append EXTRACT_D_E_FOR_NON_ALPHA in output dir name
		if self.EXTRACT_D_E_FOR_NON_ALPHA==True:
			self.OUTDIR_OBS_R2+="_w_non-alpha_DE"
		#! Append FOPT output dir name
		if mthd=="mthd2": self.OUTDIR_OBS_R2+="_%s"%self.FOPT
		if not os.path.exists(self.OUTDIR_OBS_R2):
			os.makedirs(self.OUTDIR_OBS_R2)

		#! 1. Get Q2-W binning from file
		#! + Q2-W binning in the file is as per Q2-W binning of h8s
		if self.DBG==True: 	q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX,dbg=True,dbg_bins=2)
		else:    			q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX)
		print "DispObs::disp_R2() Q2-W bins got from yield.root=",q2wbinl

		print "DispObs::disp_R2() Going to begin processing Q2-W bins from file"
		#! 2. Now, loop over Q2-W bin and per bin:
		#! 2.1. Obtain h5 
		for q2wbin in q2wbinl:
			print"DispObs::disp_R2() Processing q2wbin=",q2wbin
			
			#! 2.1. Obtain h5d(seq,vst)
			h5d=OrderedDict()
			for item in list(itertools.product(self.SEQS,VSTS)):
				seq,vst=item[0],item[1]
				h5d[seq,vst]=self.FIN.Get("%s/%s/VST%d/h5"%(q2wbin,seq,vst))
			
			#! 2.2. and 2.3. Extract and plot R2
			print "DispObs::disp_R2() Going to extract and plot R2 for mthd %s"%mthd
			self.extract_and_plot_R2(q2wbin,h5d,R2l,mthd,plotphiproj)
				
		print "Done DispObs::disp_R2()"
		print "If the progam is not terminating, then Python is probably doing \"garbage collection\"(?); Wait a while!"
		return

	def extract_and_plot_R2(self,q2wbin,h5d,R2l,mthd,plotphiproj):
		"""
		+ From h5d extract and plot hR2d using the specified method
		"""
		print "In DispObs::extract_and_plot_R2()"
		
		print "DispObs::extract_and_plot_R2() Going to extract and plot R2s=%s for mthd %s..."%(R2l,mthd)
		if mthd=='mthd1':#'h5-mply-itg'
			#! h5d=>hR2d, directly
			for R2 in R2l:
				print "DispObs::extract_and_plot_R2() Doing h5=>hR2 for mthd %s:R2=%s"%(mthd,R2)
				hR2d=OrderedDict()
				for k in h5d:
					seq,vst=k[0],k[1]
					h5=h5d[k]
					h5m=thntool.MultiplyBy(h5,H5_MFUNCD[R2],1)
					#! The following was valid when I was helicit information in E1F data
					# if R2!='D':
					# 	h5m=thntool.MultiplyBy(h5,H5_MFUNCD[R2],1)	
					# elif R2=='D':
					# 	h5m=thntool.MultiplyBy(h5,H5_MFUNCD[R2],1)
					# 	if   hel=='POS' or hel=='UNP': h5m=thntool.MultiplyBy(h5,H5_MFUNCD[R2],1)
					# 	elif hel=='NEG':               h5m=thntool.MultiplyBy(h5,H5_MFUNCD[R2],-1)	
												
					for var in VARS:
						if var=='PHI': continue
						#! post VM-fdbk
						# if vst==1 and var=='M2': continue
						# if vst==2 and var=='M1': continue
						# if vst==3 and var=='M1': continue
						hR2d[seq,vst,var]=h5m.Projection(H5_DIM[var],"E")
						#! Call Sumw2() since these errors will need to be propagated
						hR2d[seq,vst,var].Sumw2()
						#! setM1M2axisrange
						if var=='M1' or var=='M2': self.setM1M2axisrange(q2wbin,hR2d[seq,vst,var],vst,var)
						#! some histogram aesthetics
						hR2d[seq,vst,var].SetName("%s_%s_VST%d_%s"%(R2,seq,vst,var))
						#! post VM-fdbk: Label as R2 as R2^{X_ij}_{phi_i}
						#hR2d[seq,vst,var].SetTitle("%s(%s)"%(R2_NAMED[R2],VAR_NAMES[(vst,var)]))
						hR2d[seq,vst,var].SetTitle("%s^{%s}_{%s}"%(R2_NAMED[R2],VAR_NAMES[(vst,var)],VAR_NAMES[(vst,'PHI')]))
						
				#! + hR2s extracted thus far need to be divided by mthd2-dependent EF*SF, where
				#! 		+ EF=Extraction Factor and CF=Correction Factor
				#! (See handwritten notes for details)
				for k in hR2d:
					ef=EF[('mthd1',R2)]*CF[('mthd1',R2)]
					hR2d[k].Scale(1/ef)
					# if R2=='A': hR2d[k].Scale(1/(2*math.pi)) #! for A
					# else: 		hR2d[k].Scale(1/math.pi)     #! for B,C,D
										
				#! Normalize hR2
				self.norm_1D(hR2d,q2wbin)
				# if "norm" in self.VIEW:
				# 	self.norm_1D(hR2d,q2wbin)
				# else: #! just DCosTheta normalize theta distributions
				# 	for k in hR2d:
				# 		if k[2]=='THETA':
				# 			self.norm_1D_theta(hR2d[k])

				print "DispObs::extract_and_plot_R2(): Done h5d=>hR2d for mthd %s:R2=%s"%(mthd,R2)

				print "DispObs::extract_and_plot_R2(): Going to plot hR2d for mthd %s:R2=%s"%(mthd,R2)
				self.plot_R2(q2wbin,hR2d,R2)
				print "DispObs::extract_and_plot_R2(): Done plot hR2d for mthd %s:R2=%s"%(mthd,R2)	
		elif mthd=='mthd2' or mthd=='mthd3':#proj-phi-fit or proj-phi-mply-itg
			#! 1. h5d=>hphiprojd
			print "DispObs::extract_and_plot_R2(): Doing h5d=>hphiprojd for method %s..."%mthd
			
			hphiprojd=OrderedDict()
			for k in h5d:
				seq,vst=k[0],k[1]
				h5=h5d[k]
				for var in VARS:
					if var=='PHI': continue
					#! post VM-fdbk
					# if vst==1 and var=='M2': continue
					# if vst==2 and var=='M1': continue
					# if vst==3 and var=='M1': continue
					nbins=NBINS[var]
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
						#! post VM-fdbk
						#new_ttl="%s:proj. bin=[%.3f,%.3f]"%(VAR_NAMES[(vst,var)],xmin,xmax)
						new_ttl="%s:#phi_{%s} proj. bin=[%.3f,%.3f]"%(VAR_NAMES[(vst,var)],VAR_NAMES[(vst,'PHI')],xmin,xmax)
						h.SetTitle(new_ttl)
						h.SetMarkerStyle(ROOT.gROOT.ProcessLine("kFullCircle"))
						h.SetMarkerColor(CLRS_PLT_SEQ_ALL[seq])
						h.SetLineColor(CLRS_PLT_SEQ_ALL[seq])
						#! Finally, add this histogram to hphiprojd
						hphiprojd[seq,vst,var,ibin+1]=h
						#! Before leaving projection loop, reset projection range for h5!
						h5.GetAxis(H5_DIM[var]).SetRange()
			print "DispObs::extract_and_plot_R2(): Done h5d=>hphiprojd for method %s..."%mthd
			
			#! 2. hphiprojd=>fpard
			print "DispObs::extract_and_plot_R2(): Doing hphiprojd=>fpard"
			fpard=OrderedDict()
			for k in hphiprojd:
				seq,vst,var,bin=k[0],k[1],k[2],k[3]
				
				#! Get fphi
				if self.EXTRACT_D_E_FOR_NON_ALPHA:
					f=get_fphi_alpha()
				else:
					if   var!='ALPHA':	f=get_fphi()
					elif var=='ALPHA':	f=get_fphi_alpha()
				

				print "DispObs::extract_and_plot_R2(): Going to fit phi proj for",k
				fstat=int(hphiprojd[k].Fit(f,self.FOPT)) #! "NQI"
				print "DispObs::extract_and_plot_R2(): fstat=",fstat
				A,Aerr=0,0
				B,Berr=0,0
				C,Cerr=0,0
				D,Derr=0,0
				E,Eerr=0,0
				chisq_per_DOF=0
				prob=0
				#A,Aerr,B,Berr,C,Cerr=0,0,0,0,0,0

				if fstat==0:
					A,Aerr=f.GetParameter(0),f.GetParError(0)
					B,Berr=f.GetParameter(1),f.GetParError(1)
					C,Cerr=f.GetParameter(2),f.GetParError(2)
					if var=='ALPHA' or self.EXTRACT_D_E_FOR_NON_ALPHA:
						D,Derr=f.GetParameter(3),f.GetParError(3)
						E,Eerr=f.GetParameter(4),f.GetParError(4)
					#print "here",A,B,C,D,f.GetChisquare(),f.GetNDF()
					if f.GetNDF()==0: chisq_per_DOF=1000 #! I have not looked into why this happens sometimes
					else:             chisq_per_DOF=f.GetChisquare()/f.GetNDF()
					prob=f.GetProb()
				else:
					print "DispObs::extract_obs_R2(): AT-WARNING: Fit did not succeed for q2wbin,vst,var,bin=",q2wbin,vst,var,bin
					#! For now to identify these data points, use very large error bars
					Aerr,Berr,Cerr=10000,10000,10000
					if var=='ALPHA' or self.EXTRACT_D_E_FOR_NON_ALPHA:
						Derr,Err=10000,10000
					chisq_per_DOF=10000
					prob=10000

				fpard[k]={'A':A,'Aerr':Aerr,'B':B,'Berr':Berr,'C':C,'Cerr':Cerr,
				          'D':D,'Derr':Derr,'E':E,'Eerr':Eerr,
				          'chisq_per_DOF':chisq_per_DOF,'prob':prob}
				print "DispObs::extract_and_plot_R2(): A,B,C,D,E,chisq_per_DOF,prop=",A,B,C,D,E,chisq_per_DOF,prob
				#print "DispObs::extract_and_plot_R2(): A,B,C=",A,B,C
			print "Done hphiprojd=>fpard"
			
			#! 2. Obtain hR2 for each R2 from fpard and plot
			for R2 in R2l:
				print "DispObs::extract_and_plot_R2(): Going to extract R2 from phiproj fits for mthd %s:R2=%s"%(mthd,R2)
				hR2d=self.extract_R2_from_phiproj(q2wbin,h5d,fpard,R2)
				print "DispObs::extract_and_plot_R2(): Done to extract R2 from phiproj fits for mthd %s:R2=%s"%(mthd,R2)

				print "DispObs::extract_and_plot_R2(): Going to plot hR2d for mthd %s:R2=%s"%(mthd,R2)
				self.plot_R2(q2wbin,hR2d,R2)
				print "DispObs::extract_and_plot_R2(): Done plot hR2d for mthd %s:R2=%s"%(mthd,R2)

			#! 3. Finally, for visual verification, plot phiprojs
			if plotphiproj:
				print "DispObs::extract_and_plot_R2(): Going to plot phi-proj and extract R2 for method %s"%mthd
				self.plot_phiproj(q2wbin,hphiprojd,fpard)
				print "DispObs::extract_and_plot_R2(): Done to plot phi-proj and extract R2 for method %s"%mthd
			
			print "DispObs::extract_and_plot_R2(): Done h5d=>hR2d for method %s..."%mthd
		print "DispObs::extract_and_plot_R2(): Done to extract and plot R2s=%s for mthd %s..."%(R2l,mthd)
					

	def plot_R2(self,q2wbin,hR2d,R2):
		"""
		+ Plot R2s just like 1D observables are plotted
		"""

		print "In DispObs::plot_R2()"

		#! Set up some plotting related styles and aesthetics 
		self.plot_R2_athtcs()
		self.hist_R2_athtcs(hR2d,R2)
				
		#! TCanvas's pad_map[pad,vst,var] defined as per Gleb's display
		#! post VM-fdbk
		# pad_map=[(1,1,"M1"),   (2,3,'M2'),   (3,2,'M2'),
		# 		 (4,1,"THETA"),(5,3,'THETA'),(6,2,'THETA'),
		# 		 (7,1,"ALPHA"),(8,3,'ALPHA'),(9,2,'ALPHA')]
		pad_map=[(1,1,"M1"),    (2,3,'M1'),    (3,2,'M1'), 
				 (4,1,"M2"),    (5,3,'M2'),    (6,2,'M2'),
				 (7,1,"THETA"), (8,3,'THETA'), (9,2,'THETA'),
				 (10,1,"ALPHA"),(11,3,'ALPHA'),(12,2,'ALPHA')]
		npadsx=3
		npadsy=4
		npads=npadsx*npadsy
		print "DispObs::plot_R2() Plotting hR2 for R2=%s,q2wbin=%s"%(R2,q2wbin)
		c=ROOT.TCanvas("c","c",1000,1000)
		pad_t=ROOT.TPad("pad_t","Title pad",0.15,0.945,0.85,1.00)
		pad_p=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
  		pad_p.Draw()
  		pad_t.Draw()
  		pad_t.cd()
  		pt=ROOT.TPaveText(.05,.1,.95,.8)
  		pt.AddText("%s for Q2,W = %s"%(R2_NAMED[R2],q2wbin))
  		pt.SetTextSize(0.42)
  		pt.Draw()
  		pad_p.Divide(npadsx,npadsy)
  		#! TLine, for each pad, for ln at y=0
  		ln=[0 for i in range(npads)]
		for item in pad_map:
			pad,vst,var=item[0],item[1],item[2]
			print "pad,vst,var=",pad,vst,var
			gpad=pad_p.cd(pad)

			if self.EXTRACT_D_E_FOR_NON_ALPHA==False:
				if (R2=='D' or R2=='E') and var!= 'ALPHA':continue

			#! Prepare to draw histograms
			#! 1. If SF and EC in self.SEQS, then normalize its integral to that of EC
			#! [03-28-16] Normalization changed: max-amplitude of SF normalized to max-amplitude of EC
			if "SF" in self.SEQS and "EC" in self.SEQS:
				#! Integral normalization
				# print "DispObs::plot_R2() Scaling SF-itg to EC-itg"
				# itg_EC=self.get_signed_integral(hR2d['EC',vst,var])
				# #! First call Sumw2() so that errors are also scaled
				# hR2d['SF',vst,var].Sumw2()
				# #! Now get scale factor and scale bin contents
				# itg_SF=self.get_signed_integral(hR2d['SF',vst,var])
				# if itg_EC==0 or itg_SF==0: 	scale_factor=1
				# else:    					scale_factor=itg_EC/itg_SF
				# hR2d['SF',vst,var].Scale(scale_factor)

				#![03-18-16] Max-amplitude normalization
				print "DispObs::plot_R2() Scaling ampltd-SF to ampltd-EC"
				max_EC=hR2d['EC',vst,var].GetMaximum()
				min_EC=hR2d['EC',vst,var].GetMinimum()
				ampltd_EC=max([max_EC,math.fabs(min_EC)])
				#! First call Sumw2() so that errors are also scaled
				hR2d['SF',vst,var].Sumw2()
				#! Now get scale factor and scale bin contents
				max_SF=hR2d['SF',vst,var].GetMaximum()
				min_SF=hR2d['SF',vst,var].GetMinimum()
				ampltd_SF=max([max_SF,math.fabs(min_SF)])
				if ampltd_EC==0 or ampltd_SF==0: scale_factor=1
				else:    						 scale_factor=ampltd_EC/ampltd_SF
				hR2d['SF',vst,var].Scale(scale_factor)

			#! 1. If ST and EC in self.SEQS, then normalize its integral to that of EC
			#! [03-28-16] Normalization changed: max-amplitude of SF normalized to max-amplitude of EC
			if "ST" in self.SEQS and "EC" in self.SEQS:
				#! Integral normalization
				# print "DispObs::plot_R2() Scaling ST-itg to EC-itg"
				# itg_EC=self.get_signed_integral(hR2d['EC',vst,var])
				# #! First call Sumw2() so that errors are also scaled
				# hR2d['ST',vst,var].Sumw2()
				# #! Now get scale factor and scale bin contents
				# itg_ST=self.get_signed_integral(hR2d['ST',vst,var])
				# if itg_EC==0 or itg_SF==0: 	scale_factor=1
				# else:    					scale_factor=itg_EC/itg_ST
				# hR2d['ST',vst,var].Scale(scale_factor)

				#![03-18-16] Max-amplitude normalization
				print "DispObs::plot_R2() Scaling ampltd-ST to ampltd-EC"
				max_EC=hR2d['EC',vst,var].GetMaximum()
				min_EC=hR2d['EC',vst,var].GetMinimum()
				ampltd_EC=max([max_EC,math.fabs(min_EC)])
				#! First call Sumw2() so that errors are also scaled
				hR2d['ST',vst,var].Sumw2()
				#! Now get scale factor and scale bin contents
				max_ST=hR2d['ST',vst,var].GetMaximum()
				min_ST=hR2d['ST',vst,var].GetMinimum()
				ampltd_ST=max([max_ST,math.fabs(min_ST)])
				if ampltd_EC==0 or ampltd_ST==0: scale_factor=1
				else:    						 scale_factor=ampltd_EC/ampltd_ST
				hR2d['ST',vst,var].Scale(scale_factor)

			#! 2. determine and set maximum
			minl=[hR2d[seq,vst,var].GetMinimum() for seq in self.SEQS]
			maxl=[hR2d[seq,vst,var].GetMaximum() for seq in self.SEQS]
			minimum=min(minl)
			maximum=max(maxl)
			fctr_not_A=1.5
			fctr_A=1.1
			for seq in self.SEQS:		
				if R2!="A": #! B,C,D can oscillate around 0
					hR2d[seq,vst,var].SetMinimum(minimum-math.fabs(fctr_not_A*minimum))
					hR2d[seq,vst,var].SetMaximum(maximum+math.fabs(fctr_not_A*maximum))
				elif R2=="A": #! These distributions are same as Obs-1D
					hR2d[seq,vst,var].SetMinimum(0)
					hR2d[seq,vst,var].SetMaximum(fctr_A*maximum)
			#! 3. Finally draw
			if self.VIEW=="EC_SF_ST":
				hR2d['SF',vst,var].Draw()
				hR2d['ST',vst,var].Draw("same")
				hR2d['EC',vst,var].Draw("same") #! Draw EC at the end so that its error bar is on top
				#hR2d['SF',vst,var].Draw("same")
			elif self.VIEW=="EC_EF":
				hR2d['EC',vst,var].Draw()
				hR2d['EF',vst,var].Draw("same")
				#hR2d['SF',vst,var].Draw("same")			

			#! Draw TLine only if ymin<0
			if hR2d['EC',vst,var].GetMinimum()<0:
				ln[pad-1]=ROOT.TLine(hR2d['EC',vst,var].GetXaxis().GetXmin(),0,hR2d['EC',vst,var].GetXaxis().GetXmax(),0)
				ln[pad-1].Draw("same")			
			

		q2bin=self.get_q2bin(q2wbin)
		wbin=self.get_wbin(q2wbin)				
		outdir=os.path.join(self.OUTDIR_OBS_R2,"q%s"%q2bin,R2)
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		c.SaveAs("%s/c_w%s_q%s.png"%(outdir,wbin,q2bin))
		c.SaveAs("%s/c_w%s_q%s.eps"%(outdir,wbin,q2bin))
		c.Close()
		print "Done DispObs::plot_R2()"
		return

	def extract_R2_from_phiproj(self,q2wbin,h5d,fpard,R2):
		"""
		+ From fpard, for a particular R2, extract hR2
		"""
		print "In DispObs::extract_R2_from_phiproj() for R2=%s"%(R2)

		hR2d=OrderedDict()
		#! the following is a very "hack-ish" way to create hR2d by
		#! using projections from h5d and then resetting the contents of the histogram
		for k in h5d:
			seq,vst=k[0],k[1]
			h5=h5d[k]
														
			for var in VARS:
				if var=='PHI': continue
				#! post VM-fdbk
				# if vst==1 and var=='M2': continue
				# if vst==2 and var=='M1': continue
				# if vst==3 and var=='M1': continue
				hR2d[seq,vst,var]=h5.Projection(H5_DIM[var],"E")
				if var=='M1' or var=='M2': self.setM1M2axisrange(q2wbin,hR2d[seq,vst,var],vst,var)
				hR2d[seq,vst,var].SetName("%s_%s_VST%d_%s"%(R2,seq,vst,var))
				#! post VM-fdbk: Label as R2 as R2^{X_ij}_{phi_i}
				#hR2d[seq,vst,var].SetTitle("%s(%s)"%(R2_NAMED[R2],VAR_NAMES[(vst,var)]))
				hR2d[seq,vst,var].SetTitle("%s^{%s}_{%s}"%(R2_NAMED[R2],VAR_NAMES[(vst,var)],VAR_NAMES[(vst,'PHI')]))
				hR2d[seq,vst,var].Reset()
					
		print "Going to extract R2 from phi-proj fits"
		for k in fpard:
			seq,vst,var,bin=k[0],k[1],k[2],k[3]
			#fpard[hel][q2bin_le,wbin_le,vst,var,bin,dtyp,seq]=f
			#! Now fill hR2 from fit
			if R2=='A':
				r2=    fpard[k]['A']
				r2_err=fpard[k]['Aerr']
			elif R2=='B':
				r2=    fpard[k]['B']
				r2_err=fpard[k]['Berr']
			elif R2=='C':
				r2=    fpard[k]['C']
				r2_err=fpard[k]['Cerr']
			elif R2=='D':
				r2=    fpard[k]['D']
				r2_err=fpard[k]['Derr']
			elif R2=='E':
				r2=    fpard[k]['E']
				r2_err=fpard[k]['Eerr']
				#! The following was valid when I was helicit information in E1F data
				#if   hel=='POS' or hel=='UNP': r2=   fpard[hel][k]['D']
				#elif hel=='NEG':               r2=-1*fpard[hel][k]['D']
				#r2_err=fpard[hel][k]['Derr']
			hR2d[seq,vst,var].SetBinContent(bin,r2)
			hR2d[seq,vst,var].SetBinError(bin,r2_err)
		
		#! + hR2s extracted thus far need to be divided by mthd2-dependent EF*SF, where
		#! 		+ EF=Extraction Factor and CF=Correction Factor
		#! (See handwritten notes for details)
		#! + This is not needed if the Integral fit option is used
		if self.FOPT=="NQ" or self.FOPT=="NQL":
			ef=EF[('mthd2',R2)]*CF[('mthd2',R2)]
			for k in hR2d:
				hR2d[k].Scale(ef)
			# if R2=='A':
			# 	for k in hR2d:
			# 		sf=1/math.radians(36)
			# 		hR2d[k].Scale(sf)
			# else: #! i.e.B,C,D
			# 	for k in hR2d:
			# 		sf=1/math.radians(36)
			# 		hR2d[k].Scale(sf)
		elif self.FOPT=="NQI":	#! If fopt="I" i.e. the Integral fit option
			for k in hR2d:
				sf=1/math.radians(36)
				hR2d[k].Scale(sf)
		
		#! Call Sumw2() for hR2s since these errors will need to be propagated
		for k in hR2d:
			hR2d[k].Sumw2()

		#! Normalize hR2
		self.norm_1D(hR2d,q2wbin)
		# if "norm" in self.VIEW:
		# 	self.norm_1D(hR2d,q2wbin)
		# else: #! just DCosTheta normalize theta distributions
		# 	for k in hR2d:
		# 		if k[2]=='THETA':
		# 			self.norm_1D_theta(hR2d[k])

		print "DispObs::extract_R2_from_phiproj() Done extracting R2 from phi-proj fits"
		return hR2d

		print "Done DispObs::extract_R2_from_phiproj for R2=%s"%(R2)

	def plot_phiproj(self,q2wbin,hphiprojd,fpard):
		"""
		+ plot phiprojs
		"""
		print "In DispObs::plot_phiproj()"

		#! Set up ROOT-display aesthetics
		self.plot_phiproj_athtcs()

		for vst in VSTS:
			if self.DBG==True and vst!=1: continue
			for var in VARS:
				if var=='PHI': continue
				#! post VM-fdbk
				# if vst==1 and var=='M2': continue
				# if vst==2 and var=='M1': continue
				# if vst==3 and var=='M1': continue

				if self.DBG==True and var!='ALPHA': continue

				q2bin=self.get_q2bin(q2wbin)
				wbin=self.get_wbin(q2wbin)
				outdir=os.path.join(self.OUTDIR_OBS_R2,"phiprojs","VST%d_%s"%(vst,var),"q%s_w%s"%(q2bin,wbin))
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
  				pt.AddText("#phi-projections(%s:Q2,W=%s,%s)"%(VAR_NAMES[(vst,var)],q2bin,wbin))
  				pt.Draw()
 				pad_p.cd()
 				nbins=NBINS[var]
 				pad_p.Divide(NXPADS[var],NYPADS[var])
 				#! The following are used to draw the Fit, TPaveText with fit pars and the TLegend
 				f= [[0 for j in range(len(self.SEQS))]for i in range(nbins)]
 				pt=[[0 for j in range(len(self.SEQS))]for i in range(nbins)]
 				l= [[0 for j in range(len(self.SEQS))]for i in range(nbins)]
 				
 				for ibin in range(nbins):
					pad=pad_p.cd(ibin+1)
					pad.Divide(1,len(self.SEQS))
					
					for iseq,seq in enumerate(self.SEQS):
						if hphiprojd.has_key((seq,vst,var,ibin+1)):
							k=seq,vst,var,ibin+1
							h=hphiprojd[k]
							h.SetXTitle( "%s%s"%(VAR_NAMES[(vst,'PHI')],VAR_UNIT_NAMES['PHI']) )
							h.GetXaxis().SetLabelSize(.05)
							h.GetXaxis().SetTitleSize(.10)
							h.GetXaxis().SetTitleOffset(.7)
					
							if self.VIEW=="EC_SF_ST":										
								if   seq=='EC': pad.cd(1)
								elif seq=='SF': pad.cd(2)
								elif seq=='ST': pad.cd(3) 
							elif self.VIEW=="EC_EF":
								if   seq=='EC': pad.cd(1)
								elif seq=='EF': pad.cd(2)
							
							#! Draw hist
							h.SetMinimum(0)
							h.Draw()
							#! Draw Fit
							if   var!='ALPHA': f[ibin][iseq]=get_fphi("_%d_%s"%(ibin+1,seq))
							elif var=='ALPHA': f[ibin][iseq]=get_fphi_alpha("_%d_%s"%(ibin+1,seq))
							f[ibin][iseq].SetParameter(0,fpard[k]['A'])
							f[ibin][iseq].SetParError (0,fpard[k]['Aerr'])
  							f[ibin][iseq].SetParameter(1,fpard[k]['B'])
  							f[ibin][iseq].SetParError (1,fpard[k]['Berr'])
  							f[ibin][iseq].SetParameter(2,fpard[k]['C'])
  							f[ibin][iseq].SetParError (2,fpard[k]['Cerr'])
  							if var=='ALPHA':
  								f[ibin][iseq].SetParameter(3,fpard[k]['D'])
  								f[ibin][iseq].SetParError (3,fpard[k]['Derr'])
  								f[ibin][iseq].SetParameter(4,fpard[k]['E'])
  								f[ibin][iseq].SetParError (4,fpard[k]['Eerr'])
 							f[ibin][iseq].SetParName(0, "A")
  							f[ibin][iseq].SetParName(1, "B")
  							f[ibin][iseq].SetParName(2, "C")
  							if var=='ALPHA':
  								f[ibin][iseq].SetParName(3, "D")
  								f[ibin][iseq].SetParName(4, "E")
  							f[ibin][iseq].SetLineColor(h.GetMarkerColor())
							f[ibin][iseq].Draw("same")
							#! Draw Fit stats
							#pt[ibin].append(ROOT.TPaveText(.20,.6,.30,1.0,"NDC"))
							pt[ibin][iseq]=ROOT.TPaveText(.20,.6,.30,1.0,"NDC")
							pt[ibin][iseq].AddText( "A=%.2f+/-%.2f"%(f[ibin][iseq].GetParameter(0),f[ibin][iseq].GetParError(0)) )
							pt[ibin][iseq].AddText( "B=%.2f+/-%.2f"%(f[ibin][iseq].GetParameter(1),f[ibin][iseq].GetParError(1)) )
							pt[ibin][iseq].AddText( "C=%.2f+/-%.2f"%(f[ibin][iseq].GetParameter(2),f[ibin][iseq].GetParError(2)) )
							if var=='ALPHA':
								pt[ibin][iseq].AddText( "D=%.2f+/-%.2f"%(f[ibin][iseq].GetParameter(3),f[ibin][iseq].GetParError(3)) )
								pt[ibin][iseq].AddText( "E=%.2f+/-%.2f"%(f[ibin][iseq].GetParameter(4),f[ibin][iseq].GetParError(4)) )
							pt[ibin][iseq].AddText( "#chi^{2}/NDOF=%.2f"%fpard[k]['chisq_per_DOF'] )
							pt[ibin][iseq].AddText( "prob=%.2f"%fpard[k]['prob'] )
							pt[ibin][iseq].SetTextColor(h.GetMarkerColor())
							pt[ibin][iseq].SetFillColor(11)
							pt[ibin][iseq].Draw()
							
							#! Draw legend
							l[ibin][iseq]=ROOT.TLegend(0.1,0.8,0.2,0.9)
							l[ibin][iseq].AddEntry(h,seq,"p")
							l[ibin][iseq].Draw()
				c.SaveAs("%s/c_q%s_w%s.png"%(outdir,q2bin,wbin))
				c.SaveAs("%s/c_q%s_w%s.eps"%(outdir,q2bin,wbin))
				c.Close()
			print "Done DispObs::plot_phiproj()"	

	def plot_R2_athtcs(self):
		#ROOT.gStyle.Reset()
		#! Stats Box
		ROOT.gStyle.SetOptStat(0)

		# ROOT.gStyle.SetLabelSize(0.5,"t")
		# ROOT.gStyle.SetTitleSize(0.5,"t")
		#ROOT.gStyle.SetPaperSize(20,26);
		#ROOT.gStyle.SetPadTopMargin(0.15)#(0.05);
		#ROOT.gStyle.SetPadRightMargin(0.15)#(0.05);
		ROOT.gStyle.SetPadBottomMargin(0.20)#(0.16);
		ROOT.gStyle.SetPadLeftMargin(0.20)#(0.12);

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

	def setM1M2axisrange(self,q2wbin,h,vst,var):
		'''
		1,M1=p,pip 1,M2=pip,pim
		2,M1=p,pip 2,M2=pip,pim
		3,M1=p,pip 3,M2=p,pim
		'''
		#! Determine wmax
		wmax=self.get_wbin_ue(q2wbin)
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