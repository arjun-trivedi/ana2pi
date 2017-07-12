from __future__ import division
import ROOT
from rootpy.io import root_open, DoesNotExist

from collections import OrderedDict
import array
import itertools

import os,sys,subprocess
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

LUM_E16=28.18 #![08-10-17] 28.18 from From $STUDY_LUM_E16_DATADIR/results_071017: loose cut #![06-12-17] 28.19 from From $STUDY_LUM_E16_DATADIR/results_052417: loose cut #!Before [06-12-17] 28
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
VAR_NAMES={(1,VARS[0]):"M_{p#pi^{+}}",    (1,VARS[1]):"M_{#pi^{+}#pi^{-}}",
		   (1,VARS[2]):"#theta_{#pi^{-}}",(1,VARS[3]):"#phi_{#pi^{-}}",
		   (1,VARS[4]):"#alpha_{[p_{f}#pi^{+}][p#pi^{-}]}",
		   (2,VARS[0]):"M_{p#pi^{+}}",(2,VARS[1]):"M_{#pi^{+}#pi^{-}}",
		   (2,VARS[2]):"#theta_{p}",  (2,VARS[3]): "#phi_{p}",
		   (2,VARS[4]):"#alpha_{[#pi^{+}#pi^{-}][pp_{f}]}",
		   (3,VARS[0]):"M_{p#pi^{+}}",    (3,VARS[1]):"M_{p#pi^{-}}",
		   (3,VARS[2]):"#theta_{#pi^{+}}",(3,VARS[3]): "#phi_{#pi^{+}}",
		   (3,VARS[4]):"#alpha_{[p_{f}#pi^{-}][p#pi^{+}]}"}
#! Create a non-coded, simpler version of VAR_NAMES that can be used with plain text
VAR_NAMES_PLAIN={(1,VARS[0]):"M_ppip",   (1,VARS[1]):"M_pippim",
		         (1,VARS[2]):"theta_pim",(1,VARS[3]):"phi_pim",
		         (1,VARS[4]):"alpha_pfpip_ppim",
		         (2,VARS[0]):"M_ppip", (2,VARS[1]):"M_pippim",
		         (2,VARS[2]):"theta_p",(2,VARS[3]): "phi_p",
		         (2,VARS[4]):"alpha_pippim_ppf",
		         (3,VARS[0]):"M_ppip",   (3,VARS[1]):"M_ppim",
		         (3,VARS[2]):"theta_pip",(3,VARS[3]): "phi_pip",
		         (3,VARS[4]):"alpha_pfpim_ppip"}

VAR_UNIT_NAMES={VARS[0]:"GeV",VARS[1]:"GeV",VARS[2]:"deg",VARS[3]:"deg",VARS[4]:"deg"}
VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC={VARS[0]:"GeV",VARS[1]:"GeV",VARS[2]:"",VARS[3]:"rad",VARS[4]:"rad"}

#! Following dictionaries for extracting R2
R2S=['A','B','C','D','E']
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
CLRS_N=[ROOT.gROOT.ProcessLine("kRed"),    ROOT.gROOT.ProcessLine("kMagenta"),ROOT.gROOT.ProcessLine("kGreen+3"),
        ROOT.gROOT.ProcessLine("kCyan"),   ROOT.gROOT.ProcessLine("kGreen"),  ROOT.gROOT.ProcessLine("kYellow"),
        ROOT.gROOT.ProcessLine("kPink+1"), ROOT.gROOT.ProcessLine("kBlue")]
        #! Get all coordinates in this space, where coordinate def.eq. (mrkr,clr)
MRKS_CLRS_N_COORDS=list(itertools.product(MRKS_N,CLRS_N))
#! Now put them into a dictionary labelled by numbers
MRKS_CLRS_N_COORDS_DICT=OrderedDict()
for i,coord in enumerate(MRKS_CLRS_N_COORDS):
	MRKS_CLRS_N_COORDS_DICT[i+1]=[coord[0],coord[1]]

#! Alternate method inpired by having to distinguish
#! 1. upto 16 different cutsncors for Obs_1D and Obs_R2 on the same canvas
#! 2. upto 16*3 different cutsncors*vars for Itg_xsec for a particular vst
#! + Therefore, 
#! + As required by 1. simply use N different markers
#! + To accomodate 1., associate with each marker 3 different colors, one for every var
#MRKS_N2=[20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,5] #! Taken from https://root.cern.ch/doc/master/classTAttMarker.html
MRKS_N2=[20,24,21,25,22,26,23,32,29,30,33,27,34,28,31,5]
CLRS_N2=[ROOT.gROOT.ProcessLine("kRed"),ROOT.gROOT.ProcessLine("kGreen"),ROOT.gROOT.ProcessLine("kBlue")]


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
	def __init__(self,obsdir,simnum='siml',view="norm",q2min=1.25,q2max=5.25,wmin=1.400,wmax=2.125,
		         expt='e16',dbg=False,do_etgt_BG_sub_corr=True,**SSargs):
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

		#! [07-09-17] Empty target background subtraction correction factor
		#! + Only for e16
		#!
		#! + Since this SE is negative, the correction is supposed to increase the finally obtained cross-sections (itg,1D and R2)
		#!   by the absolute value of this percentage:
		#!   1. corr_factor=(1+|SE|)
		#!   2. xsec_new=xsec_old*corr_factor
		#! + Implemented in the code by modifying the 'normf' in norm_1D() as (also in disp_itg_yld_drct_from_h5(), which is no longer used):
		#!   - normf=normf/corr_factor
		#!   + In this way the corr_factor enters the cross-section formule as multiplicate factor to increase the cross-section:
		#!     - xsec = N/(normf/corr_factor) = (N/normf)*corr_factor
		#!
		#! + [07-09-17] Value obtained from $STUDY_TGT_BG_E16_DATADIR/results_R_070917/R.txt:SE of etgt BG sub for cutruns_lse:sector-average = -1.83%
		self.DO_ETGT_BG_SUB_CORR=do_etgt_BG_sub_corr
		if self.DO_ETGT_BG_SUB_CORR==True:
			if expt=='e16':
				self.ETGT_BG_SUB_CORR_FCTR=(1+0.0183) #! used multiplicatively, increases value by 1.83%
			else: #! no correction
				self.ETGT_BG_SUB_CORR_FCTR=1

		print "Q2MIN,Q2MAX=",self.Q2MIN,self.Q2MAX
		print "WMIN,WMAX=",self.WMIN,self.WMAX
		print "LUM=",self.LUM
		print "E0=",self.E0
		print "DBG=",self.DBG
		print "DO_ETGT_BG_SUB_CORR=",self.DO_ETGT_BG_SUB_CORR
		if self.DO_ETGT_BG_SUB_CORR:
			print "ETGT_BG_SUB_CORR_FCTR=",self.ETGT_BG_SUB_CORR_FCTR

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

			#! Set up SS_TYPE
			self.SS_TYPE=SSargs['SS_type']

			#! OUTDIR
			SS_tag=SSargs['SS_tag']
			date=datetime.datetime.now().strftime('%m%d%y')
			self.OUTDIR=os.path.join(self.OBSDIR,'SS','%s_%s'%(SS_tag,date))
			if self.DBG==True:
				self.OUTDIR=os.path.join(self.OBSDIR,'SS','dbg','%s_%s'%(SS_tag,date))
			if not os.path.exists(self.OUTDIR):
				os.makedirs(self.OUTDIR) 
			#! The following moved into 'disp_SE_results_1D()'
			# self.FOUT_OBS_1D=ROOT.TFile("%s/Obs_1D.root"%self.OUTDIR,"RECREATE")
			# self.FOUT_OBS_ITG=ROOT.TFile("%s/Obs_itg.root"%self.OUTDIR,"RECREATE")
			#! The following moved into 'disp_SE_results_R2()' because name of .root file
			#! depends on if EXTRACT_D_E_FOR_NON_ALPHA=True/False, which is known only to that method
			# if self.SS_TYPE!="cmb_vst_SE": #! because if so, this result is simply linked in from "cmb_non_vst_SE"
			# 	self.FOUT_OBS_R2=ROOT.TFile("%s/Obs_R2.root"%self.OUTDIR,"RECREATE")

			print "OBSDIR=",self.OBSDIR
			#! pretty-print OBSD
			print "OBSD (pretty print):"
			for k in self.OBSD:
				print k,":",self.OBSD[k][0],self.OBSD[k][1]
			print "SEQS=",self.SEQS
			print "SS_TYPE=",self.SS_TYPE
			print "OUTDIR=",self.OUTDIR
			#! The following moved into 'disp_SE_results_1D()'
			# print "FOUT_OBS_1D=",self.FOUT_OBS_1D.GetName()
			# print "FOUT_OBS_ITG=",self.FOUT_OBS_ITG.GetName()
			#! The following moved into 'disp_SE_results_R2()' because name of .root file
			#! depends on if EXTRACT_D_E_FOR_NON_ALPHA=True/False, which is known only to that method
			# if self.SS_TYPE!="cmb_vst_SE": #! because if so, this result is simply linked in from "cmb_non_vst_SE"
			# 	print "FOUT_OBS_R2=",self.FOUT_OBS_R2.GetName()
			
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
			elif self.VIEW=="norm_EC_EF_SF":
				self.SEQS=['EC','EF','SF']
			elif self.VIEW=="EC_SF_ST": #! views from this point on tuned for 'dobs_R2'
				self.SEQS=['EC','SF','ST'] 
			elif self.VIEW=="EC_EF":
				self.SEQS=['EC','EF']
			elif self.VIEW=="EC_ST": 
				self.SEQS=['EC','ST']
			elif self.VIEW=="EC_EF_ST": 
				self.SEQS=['EC','EF','ST']
			else:
				sys.exit("view={norm,norm_EC_EF_SF,fullana,ERyield,EC_SF_ST,EC_EF,EC_ST,EC_EF_ST} only")

			#! Setup related to applying rad-eff-corr
			#! if expt=='e16': set self.APPLY_RAD_EFF_CORR to True and obtain rad-eff-corr factors
			self.APPLY_RAD_EFF_CORR=False
			if expt=='e16': self.APPLY_RAD_EFF_CORR=True
			#! override for testing/debug
			#self.APPLY_RAD_EFF_CORR=False
			if self.APPLY_RAD_EFF_CORR==True:
				#! Obtain correction factors in self.RAD_EFF_CORR_FCTR[q2wbin]=[cf,cferr]
				#! + NOTE that q2wbin syntax should be similar to that in yield.root produced by proc_h8: q.qq-q.qq_w.www-w.www
				self.RAD_EFF_CORR_FCTR=OrderedDict()
				#! Open file 
				f=ROOT.TFile("%s/rad-eff-corr-sim/sim1/drad.root"%os.environ['D2PIDIR_SIM_E16'])
				#! Now read in data 
				#! + The file contains hcf_qbin_X where X=1,2,3,4,5 i.e. 5 Q2 bins for now
				for q2binnum in [1,2,3,4,5,6,7,8]:
					hname="hcf_qbin_%d"%q2binnum
					h=f.Get(hname)
					if h!=None:
						#! Get q2 bin in correct syntax from htitle
						htitle=h.GetTitle()
						#print htitle
						start=htitle.find("Q2")
						q2bin=htitle[start:]
						#print q2bin
						q2bin=q2bin.replace("Q2=","")
						q2bin=q2bin.replace("[","")
						q2bin=q2bin.replace(")","")	
						q2bin=q2bin.replace(",","-")
						#print q2bin
						for iwbin in range(h.GetNbinsX()):
							#! Get wbin in correct syntax
							wle=format(h.GetBinLowEdge(iwbin+1),'.3f')
							wue=format(h.GetBinLowEdge(iwbin+1+1),'.3f')
							wbin="%s-%s"%(wle,wue)
							#! Now get cf, cferr
							cf   =h.GetBinContent(iwbin+1)
							cferr=h.GetBinError(iwbin+1)
							self.RAD_EFF_CORR_FCTR["%s_%s"%(q2bin,wbin)]=[cf,cferr]
					else:
						print "DispObs::__init__():rad-eff-corr setup: %s not found"%hname
			
			print "expt=",expt
			print "OBSDIR=%s\nFIN=%s\nOUTDIR=%s"%(self.OBSDIR,self.FIN,self.OUTDIR)
			print "self.VIEW=",self.VIEW
			print "SEQS to process=",self.SEQS
			print "APPLY_RAD_EFF_CORR=",self.APPLY_RAD_EFF_CORR
			if self.APPLY_RAD_EFF_CORR==True:
				print "*** rad-eff-corr factors[q2bin,wbin]=[cf,cferr] ***"
				for k in self.RAD_EFF_CORR_FCTR:
					print k,self.RAD_EFF_CORR_FCTR[k]
				print "******"
			time.sleep(3)
			#! End Regular init

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

		#! For saving output to .root file, create and cd into q2wbindir
		if self.FOUT_1D!=None:
			q2wbin_dir=self.FOUT_1D.GetDirectory(q2wbin)
			if q2wbin_dir==None: 
				q2wbin_dir=self.FOUT_1D.mkdir(q2wbin)
			q2wbin_dir.cd()
				
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
				#! First set minimum and maximum of y-axis
				histl=[h1[seq,vst,var] for seq in self.SEQS]
				self.set_h1l_maximum_minimum(histl)
				#! [06-28-17] hitherto way
				# fctr=1.1
				# maxl=[h1[seq,vst,var].GetMaximum() for seq in self.SEQS]
				# maximum=max(maxl)
				# for htmp in [h1[seq,vst,var] for seq in self.SEQS]:
				# 	htmp.SetMinimum(0.)
				# 	htmp.SetMaximum(maximum*fctr)
				#! Now draw
				for i,seq in enumerate(self.SEQS):
					draw_opt="" if i==0 else "same"
					h1[seq,vst,var].Draw(draw_opt)
					#! Also write to .root file
					if self.FOUT_1D!=None:
						name="%s_"
						h1[seq,vst,var].Write()
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
				histl=[h1n[seq] for seq in seq_drct]
				self.set_h1l_maximum_minimum(histl)
				#! [06-28-17] hitherto way
				# maxl=[h1n[seq].GetMaximum() for seq in seq_drct]#[hEF_n.GetMaximum(),hSF_n.GetMaximum(),hER_n.GetMaximum(),hSR_n.GetMaximum(),hST_n.GetMaximum()]
				# maximum=max(maxl)
				# for htmp in [h1n[seq] for seq in seq_drct]:#[hEF_n,hSF_n,hER_n,hSR_n,hST_n]:
				# 	htmp.SetMinimum(0.)
				# 	htmp.SetMaximum(maximum+10)

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
			elif self.VIEW=="norm_EC_EF_SF":#! Draw SF hists normed to EF integral (since I want to compare their distributions)
				print "DispObs::plot_1D() Scaling intgrl-ST to intgrl-EF"
				intgrl_EF=h1['EF',vst,var].Integral()
				intgrl_SF=h1['SF',vst,var].Integral()
				#! First call Sumw2() so that errors are also scaled
				#! [10-26-16] + Note that SF hist normalized here is a tmp one. 
				#!            + This is because this normalization is not accurate ("old-normalization")
				#!            + The original SF histogram I will save in in its original form in the .root file
				#!              using which I can implement a "new normalization" when I make plots for the thesis
				h1_SF_tmp=h1['SF',vst,var].Clone()
				h1_SF_tmp.Sumw2()
				#! Now get scale factor and scale bin contents
				if intgrl_EF==0 or intgrl_SF==0: scale_factor=1
				else:    						 scale_factor=intgrl_EF/intgrl_SF
				h1_SF_tmp.Scale(scale_factor)

				#! Before drawing set minimum(=0) and maximum of y-axis
				histl=[h1['EC',vst,var],h1['EF',vst,var],h1_SF_tmp]
				self.set_h1l_maximum_minimum(histl)
				#! [06-28-17] hitherto way
				# fctr=1.1
				# #! [10-26-16] code below had to be modified following "new-normalization" change of [10-26-16]
				# #maxl=[h1[seq,vst,var].GetMaximum() for seq in self.SEQS]
				# maxl=[h1['EC',vst,var].GetMaximum(),h1['EF',vst,var].GetMaximum(),h1_SF_tmp.GetMaximum()]
				# maximum=max(maxl)
				# #for htmp in [h1[seq,vst,var] for seq in self.SEQS]:
				# for htmp in [h1['EC',vst,var],h1['EF',vst,var],h1_SF_tmp]:
				# 	htmp.SetMinimum(0.)
				# 	htmp.SetMaximum(maximum*fctr)

				#! Now draw
				h1['EC',vst,var].Draw()
				h1['EF',vst,var].Draw("same")
				h1_SF_tmp.Draw("same")
				#! Also write to .root file
				if self.FOUT_1D!=None:
					h1['EC',vst,var].Write()
					h1['EF',vst,var].Write()
					h1['SF',vst,var].Write()

				gpad.Update()

				#! Add TLegend if pad==1
				if pad==1:
					l=ROOT.TLegend(0.75,0.50,0.90,0.90)
					#l.SetFillStyle(0)
					l.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
					l.SetBorderSize(0)
					l.SetTextSize(0.06)
					for seq in self.SEQS:
						l.AddEntry(h1[seq,vst,var],seq,"p")#EF
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

	def get_cmb_vst_SE_results_1D(self,h1,results_EC,results_EF,q2wbin,h1err):
		'''
		[09-08-16]
		+ get_cmb_non_vst_SE_results() was used as a template for this method and that may be reflect in the code

		+ Format for results_EC/EF in this function (results refer to those obtained for vst_SE, for passed q2wbin, from Obs_itg): 
		  results_EC/EF=[f_vst_1,sg_f_vst_1,
		                 f_vst_2,sg_f_vst_2,
		                 f_vst_2,sg_f_vst_2,
		                 sg_vst_1D,sg_sg_vst_1D]

		+ Process
		'''

		print "In DispObs::plot_cmb_vst_SE_results_1D()"

		#! Get q2wbin information
		q2bin=self.get_q2bin(q2wbin)
		wbin=self.get_wbin(q2wbin)

		#! Get vst_SE results, per seq, for this qwbin
		#! Create structure for, per seq, vst related normalization factors and its error:
		#! + f[seq,vst]
		#! + sg_f[seq,vst]
		#! + sg_vst_1D[seq],sg_sg_vst_1D[seq]
		f   =OrderedDict()
		sg_f=OrderedDict()
		sg_vst_1D=OrderedDict()
		sg_sg_vst_1D=OrderedDict()
		#! Now fill the structures
		for seq in self.SEQS:
			if   seq=='EC':results=results_EC
			elif seq=='EF':results=results_EF
			else: print "ERROR DispObs::plot_cmb_vst_SE_results_1D():seq=%s not recognized"

			f[seq,1],sg_f[seq,1]=results[4],results[5]
			f[seq,2],sg_f[seq,2]=results[6],results[7]
			f[seq,3],sg_f[seq,3]=results[8],results[9]
			sg_vst_1D[seq],sg_sg_vst_1D[seq]=results[10],results[11]

		#! debug structure
		# print "*** Results received for q2bin,wbin=%s,%s: ***"%(q2bin,wbin)
		# for seq in self.SEQS:
		# 	if   seq=='EC':results=results_EC
		# 	elif seq=='EF':results=results_EF
		# 	else: print "ERROR DispObs::plot_cmb_vst_SE_results_1D():seq=%s not recognized"

		# 	print "For seq=%s:"%seq
		# 	print "f_vst_1,sg_f_vst_1=",f[seq,1],sg_f[seq,1]
		# 	print "f_vst_2,sg_f_vst_2=",f[seq,2],sg_f[seq,2]
		# 	print "f_vst_3,sg_f_vst_3=",f[seq,3],sg_f[seq,3]
		# 	print "sg_vst_1D,sg_sg_vst_1D=",sg_vst_1D[seq],sg_sg_vst_1D[seq]
		# print "******"

		#! Apply vst_SE to Obs_1D to h1/h1err[obsnum,seq,vst,var] to obtain:
		#! (NOTE, obsnum in h1/h1err is always=1)
		#! + h1_obs[seq,vst,var]
		#! + h1_err[seq,vst,var]
		#! + h1_rel_err[seq,vst,var]

		#! Create output object dicts
		h1_obs,h1_err,h1_rel_err,rslts=OrderedDict(),OrderedDict(),OrderedDict(),OrderedDict()		

		for item in list(itertools.product(self.SEQS,VSTS)):
			seq,vst=item[0],item[1]
			if   vst==1:           varl=['M1','THETA','ALPHA']
			elif vst==2 or vst==3: varl=['M2','THETA','ALPHA']
			for var in varl:
				#! Create output-h1s[seq,vst,var] by Cloning and then resetting
				#! the contents of the h1 of 1st obsnum in self.OBSD i.e. h1[1,seq,vst,var]
				h1_obs[seq,vst,var]=h1[1,seq,vst,var].Clone()
				h1_err[seq,vst,var]=h1[1,seq,vst,var].Clone()
				h1_rel_err[seq,vst,var]=h1[1,seq,vst,var].Clone()
				for h in [h1_obs[seq,vst,var],h1_err[seq,vst,var],h1_rel_err[seq,vst,var]]:
					h.Reset()
					h.SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
					h.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
					h.SetMarkerStyle(ROOT.gROOT.ProcessLine("kFullDotLarge"))
				#! Create output-rslts
				rslts[seq,vst,var]=[]
				#! Now loop over all bins and calculate results using vst_SE results passed in
				nbins=h1_obs[seq,vst,var].GetNbinsX()
				for ibin in range(nbins):
					#! Now calculate results (see handwritten notes for formulae)
					#! mu,sg_mu: Renormalize using renorm factor f[seq,var]
					#! Get pre norm mu,sg_mu 
					mu_old   =h1[1,seq,vst,var].GetBinContent(ibin+1)
					sg_mu_old=h1[1,seq,vst,var].GetBinError(ibin+1)
					#! Now renormalize them
					mu=(mu_old*f[seq,vst]) 
					if mu==0: #! This can happen if cross-section=0 as can be for M1,M2, and perhaps even others
						sg_mu=0 #! this is even mathematically true
					else:
						sg_mu=mu*math.sqrt( (sg_mu_old/mu_old)**2 + (sg_f[seq,vst]/f[seq,vst])**2 )
					# if self.DBG:
					# 	if n==0: 
					# 		print "DispObs::get_per_non_vst_SE_results(): n=0 for seq,vst,var,bin=",seq,vst,var,ibin+1
					# 	print "**** DispObs::get_per_non_vst_SE_results():: binc,binerr for seq,vst,var,bin=",seq,vst,var,ibin+1,"***"
					# 	print "mu,sg_mu=",mu,sg_mu
					# 	print "******"

					#! sg,sg_sg: Add sg_vst_1D in quadrature to sg_old
					#! Get sg_old,sg_sg_old
					sg_old   =h1err[1,seq,vst,var].GetBinContent(ibin+1)
					sg_sg_old=h1err[1,seq,vst,var].GetBinError(ibin+1)
					#! Calculate new sg,sg_sg
					#! Get binw
					binw=self.get_binw(h1[1,seq,vst,var],var)
					#print "binw for var=%s = %f"%(var,binw)
					sg=math.sqrt( sg_old**2 + (sg_vst_1D[seq]/(math.sqrt(nbins))*binw)**2 )
					if sg==0: #! This is true even in the analytical formulae used for cases when sg!=0
						sg_sg=0
					else:
						sg_sg_1=( 1/sg )
						sg_sg_2=math.sqrt( (sg_old**2*sg_sg_old**2) + ( (sg_vst_1D[seq]/(math.sqrt(nbins))*binw)**2*sg_sg_vst_1D[seq]**2) ) 
						sg_sg=sg_sg_1*sg_sg_2 

					#! rel_err,sg_rel_err
					if mu==0 or sg==0: 
						rel_err=0
						sg_rel_err=0
					else:
						rel_err=(sg/mu)*100
						sg_rel_err=(rel_err)*math.sqrt((sg_sg/sg)**2+(sg_mu/mu)**2)
					#! Store bin information (will be used for passing for full results rslts dict)
					bin_le=h1_obs[seq,vst,var].GetBinLowEdge(ibin+1)
					bin_ue=h1_obs[seq,vst,var].GetBinLowEdge(ibin+1+1)

					#! Now fill results in output objects
					#! mu,sg_mu => h1_obs
					h1_obs[seq,vst,var].SetBinContent(ibin+1,mu)
					h1_obs[seq,vst,var].SetBinError(ibin+1,sg_mu)
					#! sg,sg_sg => h1_err
					h1_err[seq,vst,var].SetBinContent(ibin+1,sg)
					h1_err[seq,vst,var].SetBinError(ibin+1,sg_sg)
					#! rel_err,sg_rel_err => h1_rel_err
					h1_rel_err[seq,vst,var].SetBinContent(ibin+1,rel_err)
					h1_rel_err[seq,vst,var].SetBinError(ibin+1,sg_rel_err)
					#! [bin,bin_le,bin_ue,mu,sg_mu,sg,sg_sg,rel_err,sg_rel_err] => rslts
					rslts[seq,vst,var].append([ibin+1, bin_le, bin_ue, mu, sg_mu, sg, sg_sg, rel_err, sg_rel_err])
				
		print "*** get_cmb_vst_SE_results_1D() Done ***\n"

		return h1_obs,h1_err,h1_rel_err,rslts

	def plot_SE_results_1D(self,h1,results,q2wbin,h1err):
		'''
		[09-04-16]
		+ 'h1err' added to accomodate plotting of "cmb_non_vst_SE"

		[08-29-16]
		+ 'plot_1D' was used as a tempate for this method
		+ The idea is extend the method to incorporate plotting of a 
		  list of SE-Obs and the associated set of 'results'
		+ Therefore only details relevant to the above were changed and the core functionality left intact

		+ The following are the ways in which the results are displayed:
			1. EH: display all 'results' for 'EH'
			2. EC: write only text part of 'results'
			3. EC_EH comparison: Display only 'h1_obs'*['EC','EH']
			4. Results stored in FOUT_OBS_1D: [EC,EH]* [h1_obs,h1_err,h1_rel_err] 
		
		'''

		print "In DispObs::plot_SE_results_1D()"

		#! unpack results(=h1_obs,h1_err,h1_rel_err,rslts)
		h1_obs,h1_err,h1_rel_err,rslts=results[0],results[1],results[2],results[3]

		#! Get q2wbin information
		q2bin=self.get_q2bin(q2wbin)
		wbin=self.get_wbin(q2wbin)

		#! TCanvas's pad_map[pad,vst,var] defined as per Gleb's display for Obs_1D
		#! + Note that this pad_map is also used by non-canvas,related displays (for example when writing text
		#!   version of resutls, div,, or when writing results to FOUT, 4.)  simply due to its value of organizing Obs_1D
		pad_map=[(1,1,"M1"),   (2,3,'M2'),   (3,2,'M2'),
				 (4,1,"THETA"),(5,3,'THETA'),(6,2,'THETA'),
				 (7,1,"ALPHA"),(8,3,'ALPHA'),(9,2,'ALPHA')]

		#! Set up plotting aesthetics
		self.plot_1D_athtcs()
		self.hist_1D_athtcs(h1)
		if self.SS_TYPE=="cmb_non_vst_SE" or self.SS_TYPE=="cmb_vst_SE":
			self.hist_1D_athtcs(h1err,"err")
		self.hist_1D_athtcs(h1_obs)
		self.hist_1D_athtcs(h1_err,"err")
		self.hist_1D_athtcs(h1_rel_err,"rel_err")
				
		#! 1. Make all displays (di,dii,diii and div) for only for EH 
		#! +  Note that text results of the displays, i.e. div., which also contain EC results, will be created outside
		#!    of the following loop, where "EC" is skipped
		#! +  For the same reason,EC-EH comparison, see 2., is also done outside this loop
		#! 
		#! + Setup variable: nlegend_canvas=0
		#! + The idea is to draw one master legend that will be used for all plots, Obs-1D and intg_xsec
		#! + This is because this legend is too large to draw on the same canvas as the plots.
		#! + Once this master legend in drawn and saved, no other such legend will be drawn.
		#! nlegend canvas
		nlegend_canvas=0	
		for seq in self.SEQS: #! di,dii and diii in this loop. div outside of loop since "EC" is skipped here.			  
			if seq=="EC": continue
			#! di. Canvas for h1
			c=ROOT.TCanvas("c","c",1000,1000)
			pad_t_c=ROOT.TPad("pad_l","Title pad",0.25,0.935,0.75,1.00)
			pad_p_c=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
			#pad_p_c.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
			pad_p_c.Draw()
			pad_t_c.Draw()
			pad_t_c.cd()
			pt=ROOT.TPaveText(.05,.1,.95,.8)
			pt.AddText("Q2_W bin=%s"%q2wbin)
			pt.SetTextSize(0.40)
			pt.Draw()
			pad_p_c.Divide(3,3)

			#! dii. Canvas for h1_obs and h1_err
			c_obs_err=ROOT.TCanvas("c_obs_err","c_obs_err",1000,1000)
			pad_t_c_obs_err=ROOT.TPad("pad_l","Title pad",0.25,0.935,0.75,1.00)
			pad_p_c_obs_err=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
			pad_p_c_obs_err.Draw()
			pad_t_c_obs_err.Draw()
			pad_t_c_obs_err.cd()
			# pt=ROOT.TPaveText(.05,.1,.95,.8)
			# pt.AddText("Q2_W bin=%s"%q2wbin)
			# pt.SetTextSize(0.40)
			pt.Draw()
			pad_p_c_obs_err.Divide(3,3)

			#! diii. Canvas for h1_rel_err
			c_err=ROOT.TCanvas("c_err","c_err",1000,1000)
			pad_t_c_err=ROOT.TPad("pad_l","Title pad",0.25,0.935,0.75,1.00)
			pad_p_c_err=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
			pad_p_c_err.Draw()
			pad_t_c_err.Draw()
			pad_t_c_err.cd()
			# pt=ROOT.TPaveText(.05,.1,.95,.8)
			# pt.AddText("Q2_W bin=%s"%q2wbin)
			# pt.SetTextSize(0.40)
			pt.Draw()
			pad_p_c_err.Divide(3,3)


			for item in pad_map:
				pad,vst,var=item[0],item[1],item[2]
				print "pad,vst,var=",pad,vst,var
						
				#! Draw h1[obsnum,seq,vst,var] for visual verification
				#! + if self.SSTYPE=="cmb_non_vst_SE" Draw h1err[obsnum,seq,vst,var]
				#! + Draw h1_obs,h1_err also to directly verify mu,sg etc calculated from obsum
				#! First set minimum and maximum of hists
				histl=[h1[obsnum,seq,vst,var] for obsnum in self.OBSD.keys()]
				self.set_h1l_maximum_minimum(histl)
				# #! [06-28-17] hitherto way
				# fctr=1.1
				# maxl=[h1[obsnum,seq,vst,var].GetMaximum() for obsnum in self.OBSD.keys()]
				# maximum=max(maxl)
				# for obsnum in self.OBSD.keys():
				# 	h1[obsnum,seq,vst,var].SetMinimum(0.)
				# 	h1[obsnum,seq,vst,var].SetMaximum(fctr*maximum)
				#! Now draw
				gpad=pad_p_c.cd(pad)
				for i,obsnum in enumerate(self.OBSD.keys()):
					draw_opt="" if i==0 else "same"
					h1[obsnum,seq,vst,var].Draw(draw_opt)
					if self.SS_TYPE=="cmb_non_vst_SE" or self.SS_TYPE=="cmb_vst_SE":
						h1err[obsnum,seq,vst,var].Draw("hist e same")

				#! Draw h1_obs and h1_err also 	
				h1_obs[seq,vst,var].Draw("same")
				h1_err[seq,vst,var].Draw("hist e same")

				#! legend canvas
				if nlegend_canvas==0:
					#! Setup legend
					lm=ROOT.TLegend(0,0,1,1)
					lm.SetFillStyle(0)
					lm.SetBorderSize(0)
					lm.SetTextSize(0.03)
					for k in self.OBSD.keys():
						obsnum,obs,obstag=k,self.OBSD[k][0],self.OBSD[k][1]
						if self.SS_TYPE=="cmb_vst_SE":
							lm.AddEntry(h1[obsnum,seq,vst,var],"%s"%(obstag),"p")
							lm.AddEntry(h1_obs[seq,vst,var],   "%s"%(self.SS_TYPE),"p")
						else:
							lm.AddEntry(h1[obsnum,seq,vst,var],"%s"%(obstag),"p")
					
					clm=ROOT.TCanvas()
					clm.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
					lm.Draw()
					clm.SaveAs("%s/c_master_legend.png"%self.OUTDIR)
					clm.SaveAs("%s/c_master_legend.eps"%self.OUTDIR)
					clm.Close()
					nlegend_canvas+=1

				#! Draw h1_obs and h1_err
				gpad=pad_p_c_obs_err.cd(pad)
				#! Make sure hist limits are Ok
				histl=[h1_obs[seq,vst,var],h1_err[seq,vst,var]]
				self.set_h1l_maximum_minimum(histl)
				#! [06-28-17] hitherto way
				# self.set_h1_maximum_minimum(h1_obs[seq,vst,var],10/100,minimum=0)
				# self.set_h1_maximum_minimum(h1_err[seq,vst,var],10/100,minimum=0)
				#! Draw
				h1_obs[seq,vst,var].Draw()
				h1_err[seq,vst,var].Draw("hist e same")

				#! Draw h1_rel_err
				gpad=pad_p_c_err.cd(pad)
				#! Make sure hist limits are Ok
				self.set_h1_maximum_minimum(h1_rel_err[seq,vst,var],10/100)
				#! [06-28-17] hitherto way
				# self.set_h1_maximum_minimum(h1_rel_err[seq,vst,var],10/100,minimum=0)
				#h1_rel_err[seq,vst,var].SetMaximum(40)
				#! Draw
				h1_rel_err[seq,vst,var].Draw()
				

			#! Save canvases for h1,h1_obs & h1_err, h1_rel_err
			#! c_obs_err => self.OUTDIR_1D
			#! c_err     => self.OUTDIR_1D_ERR
			#! c         => self.OUTDIR_1D_ERR_VSL
			for cvs,outdir in zip([c_obs_err,c_err,c],[self.OUTDIR_1D,self.OUTDIR_1D_ERR,self.OUTDIR_1D_ERR_VSL]):
				outdir_q2bin=os.path.join(outdir,"q%s"%q2bin)
				if not os.path.exists(outdir_q2bin):
					os.makedirs(outdir_q2bin)
				cvs.SaveAs("%s/c_w%s_q%s.png"%(outdir_q2bin,wbin,q2bin))
				cvs.SaveAs("%s/c_w%s_q%s.eps"%(outdir_q2bin,wbin,q2bin))

		#! + div. Save results in text file for EC and EH
		#!   rslts        => self.OUTDIR_1D_TXT_RSLT/q2bin/wbin/
		#!   rslts-EC     => self.OUTDIR_1D_TXT_RSLT/EC/q2vin/wbin
		#! + Write out full results (9 pieces of information) in text file as:
		#!   bin,bin_le,bin_ue = mu,sg_mu,sg,sg_sg,rel_err,sg_rel_err
		for seq in self.SEQS:
			#! as per seq, choose outdir
			if seq=="EC":
				outdir_txt_q2bin_wbin=os.path.join(self.OUTDIR_1D_TXT_RSLT,"EC","q%s"%q2bin,"w%s"%wbin)
			elif seq=="EF":
				outdir_txt_q2bin_wbin=os.path.join(self.OUTDIR_1D_TXT_RSLT,"q%s"%q2bin,"w%s"%wbin)
			else:
				print "ERROR: plot_SE_results_1D(): seq=%s not recognized while printing txt results"%seq
			if not os.path.exists(outdir_txt_q2bin_wbin):
				os.makedirs(outdir_txt_q2bin_wbin)
			#! Print results
			for item in pad_map:
				pad,vst,var=item[0],item[1],item[2]
				#ftxt=open("%s/%s.txt"%(outdir_txt_q2bin_wbin,VAR_NAMES_PLAIN[vst,var]),"w")
				ftxt=open("%s/%d_%s.txt"%(outdir_txt_q2bin_wbin,vst,var),"w")
				#! Loop over rslts per bin in rslts and write to file
				for d in rslts[seq,vst,var]:
					bin,bin_le,bin_ue,mu,sg_mu,sg,sg_sg,rel_err,sg_rel_err=d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8]
					ftxt.write("%d,%f,%f = (%f +/- %f),(%f +/- %f),(%f +/- %f)\n"%(bin,bin_le,bin_ue,mu,sg_mu,sg,sg_sg,rel_err,sg_rel_err))
				ftxt.close()
		
		#! 2. EC-EH comparison: h1_obs only
		#! Canvas for h1
		c=ROOT.TCanvas("c","c",1000,1000)
		pad_t_c=ROOT.TPad("pad_l","Title pad",0.25,0.935,0.75,1.00)
		pad_p_c=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
		#pad_p_c.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
		pad_p_c.Draw()
		pad_t_c.Draw()
		pad_t_c.cd()
		pt=ROOT.TPaveText(.05,.1,.95,.8)
		pt.AddText("Q2_W bin=%s"%q2wbin)
		pt.SetTextSize(0.40)
		pt.Draw()
		pad_p_c.Divide(3,3)	

		for item in pad_map:
			pad,vst,var=item[0],item[1],item[2]
			print "pad,vst,var=",pad,vst,var
						
			#! Draw h1_obs
			gpad=pad_p_c.cd(pad)
			for i,seq in enumerate(['EF','EC']): #! self.SEQS purposely modified so that 'EF' is drawn first 
				#! Make sure hist limits are Ok
				self.set_h1_maximum_minimum(h1_obs[seq,vst,var],10/100)
				#! hitherto way
				# self.set_h1_maximum_minimum(h1_obs[seq,vst,var],10/100,minimum=0)
				#! Setup aesthetics as per EC-EH comparison
				h1_obs[seq,vst,var].SetMarkerColor(CLRS_PLT_SEQ_ALL[seq])
				h1_obs[seq,vst,var].SetLineColor(CLRS_PLT_SEQ_ALL[seq])
				draw_opt=""
				if i>0: draw_opt="same"
				#! Draw
				h1_obs[seq,vst,var].Draw(draw_opt)

			#! Draw TLegend if pad==1 (since this pad contain Mppip distribution that tapers off towards the right)
			if pad==1:
				l=ROOT.TLegend(0.7,0.75,0.9,0.90)
				for seq in self.SEQS:
					l.AddEntry(h1_obs[seq,vst,var],seq,"p")
				l.Draw()

		outdir_q2bin=os.path.join(self.OUTDIR_1D,"EC_EH","q%s"%q2bin)
		if not os.path.exists(outdir_q2bin):
			os.makedirs(outdir_q2bin)
		c.SaveAs("%s/c_w%s_q%s.png"%(outdir_q2bin,wbin,q2bin))
		c.SaveAs("%s/c_w%s_q%s.eps"%(outdir_q2bin,wbin,q2bin))					

		#! 4. Results stored in FOUT_OBS_1D: [EC,EH]* [h1_obs,h1_err,h1_rel_err] 
		#! output structure: q2wbin/h1_obs(err/rel_err)_seq_vst_var
		fout_q2wdir=self.FOUT_OBS_1D.GetDirectory(q2wbin)
		if fout_q2wdir==None:
			fout_q2wdir=self.FOUT_OBS_1D.mkdir(q2wbin)
		fout_q2wdir.cd()

		#! Now save [EC,EH]*h1_obs,h1_err,h1_rel_err
		for item in pad_map:
			pad,vst,var=item[0],item[1],item[2]
			print "pad,vst,var=",pad,vst,var
			#! per seq, obtain h1_obs and h1_err and save them appropriately
			for i,seq in enumerate(self.SEQS): 
				#! Save hists
				#! + Note the hists have to be named, again because:
				#! 		i. somewhere earlier their titles are reset to "" 
				#!		AND, more importantly,
				#!		ii. 'obs','err', and 'rel_err' need to be added to their names
				h1_obs[seq,vst,var].SetName("h1_obs_%s_%d_%s"%(seq,vst,var))
				h1_obs[seq,vst,var].Write()
				h1_err[seq,vst,var].SetName("h1_err_%s_%d_%s"%(seq,vst,var))
				h1_err[seq,vst,var].Write()
				if seq=="EC": #! then set maximum (note already done so for EF in canvas displays)
					self.set_h1_maximum_minimum(h1_rel_err[seq,vst,var],10/100)
					#! [06-2-17] hitherto way
					# self.set_h1_maximum_minimum(h1_rel_err[seq,vst,var],10/100,minimum=0)
				h1_rel_err[seq,vst,var].SetName("h1_rel_err_%s_%d_%s"%(seq,vst,var))		
				h1_rel_err[seq,vst,var].Write()

		print "*** plot_SE_results_1D() Done ***\n"
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
					+ =EC,EF,SF    if view="norm_EC_EF_SF"
				+ var:	
					+ if vst==1:           VARS=(M1,THETA,ALPHA)
					+ if vst==2 or vst==3: VARS=(M2,THETA,ALPHA)

			+ If self.VIEW="norm" or self.VIEW=="norm_EC_EF_SF", then first Normalize and then, i.e. secondly, apply rad-effect-corr h1(seq,vst,var)
			+ Fill relevant bin of h2 (=Intg_yld)
			+ Plot h1(seq,vst,var)
		4. Plot h2(seq,vst,var) and its projection on W
		"""
		print "*** In DispObs::disp_1D() ***"

		#! Make OUTDIR_1D and .root for output
		#self.OUTDIR_1D=os.path.join(self.OUTDIR,"Obs_1D%s"%("_norm" if self.VIEW==True else ""))
		self.OUTDIR_1D=os.path.join(self.OUTDIR,"Obs_1D_%s"%self.VIEW)
		if not os.path.exists(self.OUTDIR_1D):
			os.makedirs(self.OUTDIR_1D)
		self.FOUT_1D=None
		if self.VIEW=="norm" or self.VIEW=="norm_EC_EF_SF":
			self.FOUT_1D=ROOT.TFile(os.path.join(self.OUTDIR_1D,"obs_1D.root"),"RECREATE")

		#! Make OUTDIR_2D and .root and .txt file for output
		#self.OUTDIR_ITG_YLD=os.path.join(self.OUTDIR,"Obs_Itg_Yld%s"%("_norm" if self.VIEW==True else ""))
		self.OUTDIR_ITG_YLD=os.path.join(self.OUTDIR,"Obs_Itg_Yld_%s"%self.VIEW)
		if not os.path.exists(self.OUTDIR_ITG_YLD):
			os.makedirs(self.OUTDIR_ITG_YLD)
		self.FOUT_ITG_YLD=None
		if self.VIEW=="norm":
			self.FOUT_ITG_YLD=ROOT.TFile(os.path.join(self.OUTDIR_ITG_YLD,"obs_itg_yld.root"),"RECREATE")

		#! + Get Q2-W binning from file
		#! + Q2-W binning in the file is as per Q2-W binning of h8s
		if self.DBG==True:
			q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX,dbg=True,dbg_bins=3)
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

					#! Give h1 its unique name
					#! + This has to be done here because while in FIN h1(var) are organized in
					#!   their unique seq/vst folders, in this program that structure is lost
					#!   and therefore they have to be directly given unique names
					h1[seq,vst,var].SetName("h1_%s_%s_%s"%(seq,vst,var))

					#! Call Sumw2() since the errors as currently in the h1s
					#! need to be propagated in all subsequent calculations that involve them
					h1[seq,vst,var].Sumw2()
			#! End of [seq,vst,var] loop

			#! + If "norm", then (1.)normalize and then (2.) apply_rad_eff_corr to h1(SEQ,VSTS,VARS)
			#! + If  "norm_EC_EF_SF", then also do the same, even for SF, since neither 1 or 2
			#!   on SF should be adverse to comparing its distribution with EC & EF. In fact,
			#!   applying rad_eff_corr on SF will only make this comparison truer
			if self.VIEW=="norm":
				self.norm_1D(h1,q2wbin)
				if self.APPLY_RAD_EFF_CORR==True:
					self.apply_rad_eff_corr(h1,q2wbin)
			elif self.VIEW=="norm_EC_EF_SF": 
				self.norm_1D(h1,q2wbin)
				if self.APPLY_RAD_EFF_CORR==True:
					self.apply_rad_eff_corr(h1,q2wbin)
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
		if self.VIEW=="norm":
			self.plot_h2_itg_yld_all_Q2_same_canvas(h2)

		print "*** disp_1D() Done ***\n"
		print "If the progam is not terminating, then Python is probably doing \"garbage collection\"(?); Wait a while!"


	def disp_SE_results_1D(self):
		"""
		[08-24-16]
		+ This function displays results, which includes errors, for Obs:itg,1D based on all variations
		for "per_non_vst_SE" and "cmb_non_vst_SE"

		+ This function uses the following output .root file from observables in self.OBSD:
			+ if self.SS_TYPE=="per_non_vst_SE"
		  		+ For Obs_1D:   Obs_1D_norm/obs_1D.root
		  		+ For Itg-xsec: Obs_Itg_Yld_norm/obs_itg_yld.root
		  	+ if self.SS_TYPE=="cmb_non_vst_SE"
		  		+ For Obs_1D:   Obs_1D.root
		  		+ For Itg-xsec: Obs_itg.root

		-1. First make sure SS_TYPE is valid

		0. Make all needed output directories
		
		1. Make displays for Obs_1D
		1.i.   Get Q2-W binning from any of the obs_1D.root files in self.OBSD
		1.ii.  Looper over Q2-W bins and for every Q2-W bin: 
			1.iii  Create:
				+ if self.SS_TYPE=="per_non_vst_SE": 
					+ h1(obsnum,seq,vst,var)
					+ h1err=None ('err' refers to Systematic Errors, which at this stage are not calculated)
				where h1 contains mu,sg_mu which at this point are xsec and its statistical error in each bin.	
				+ if self.SS_TYPE=="cmb_non_vst_SE": 
					+ h1(obsnum,seq,vst,var)     [renamed here from more appropriate named saved in file as h1_obs(obsnum,seq,vst,var)]
					+ h1err(obsnum,seq,vst,var)  [renamed here from more appropriate named saved in file as h1_err(obsnum,seq,vst,var)]
				where h1 and h1err contain mu,sg_mu (<xsec> and its error) and sg,sg_sg (RMS_<xsec> and it error)
				from "per_non_vst_SE" procedure. 
			1.iv.  For a [seq,vst,var], using all obsnum, obtain 'results'(seq,vst,var) from 
			       + if self.SS_TYPE=="per_non_vst_SE": 
			       		+ 'get_per_non_vst_SE_results(h1)'
			       + if self.SS_TYPE=="cmb_non_vst_SE":
			       		+ 'get_cmb_non_vst_SE_results(h1,h1err)'

			        where 'results' are (see functions for details):
			       	+ h1_obs[seq,vst,var]
			       	+ h1_err[seq,vst,var]
			       	+ h1_rel_err[seq,vst,var]
			       	+ rslts[seq,vst,var]

			1.v.   Display: call plot_SE_results_1D(h1,rslts,q2wbin,h1_err) to plot:
					di.   Final result(mu,sg_mu,sg,sg_sg_sg) from SE variations: h1_obs[seq,vst,var] and h1_err[seq,vst,var] together
					dii.  To just see relative error (rel_err, sg_rel_err):      h1_rel_err[seq,vst,var]
					diii. To visually see SE variations:                         h1(obsnum,seq,vst,var) directly
					div.  Text file of full results:                             put rslts[seq,vst,var] in a text file
		
		2. Make displays for Obs_itg
		2.i.   Get Q2-W bin list from any of the obs_itg_yld.root files from self.OBSD
		2.ii.  Create:
			+ if self.SS_TYPE=="per_non_vst_SE": 
				+ hW(q2wbin,obsnum,seq,vst,var)
				+ hWerr=None (err refers to Systematic Errors, which at this stage are not calculated)
			where hW contains mu,sg_mu which at this point are xsec and its statistical error in each bin.
			+ if self.SS_TYPE=="cmb_non_vst_SE": 
				+ hW(q2wbin,obsnum,seq,vst,var)     [renamed here from more appropriate named saved in file as hW_obs(q2wbin,obsnum,seq,vst,var)]
				+ hWerr(q2wbin,obsnum,seq,vst,var)  [renamed here from more appropriate named saved in file as hW_err(q2wbin,obsnum,seq,vst,var)]
			where hW and hWerr contain mu,sg_mu and sg,sg_sg from "per_non_vst_SE"	
		2.iii. Make plots by calling 'plot_SE_results_W(hW)' which calculates 'results' by calling
				+ if self.SS_TYPE=="per_non_vst_SE": 
					+ 'get_per_non_vst_SE_results(h1)'
			    + if self.SS_TYPE=="cmb_non_vst_SE": 
			    	+ 'get_cmb_non_vst_SE_results(h1,h1err)'
		"""
		print "*** In DispObs::disp_SE_results_1D() ***"

		#! -1. Make sure SS_TYPE is valid
		if self.SS_TYPE!="per_non_vst_SE" or self.SS_TYPE!="cmb_non_vst_SE":
			print "DispObs::disp_SE_results_1D(): self.SS_TYPE=%s is not valid. Exiting."%self.SS_TYPE

		#! 0. Make output directories
		#! Make OUTDIR_1D,OUTDIR_1D_ERR,OUTDIR_1D_ERR_VSL,OUTDIR_1D_TXT_RSTLS
		self.OUTDIR_1D              =os.path.join(self.OUTDIR,"Obs_1D")          #! Display i.
		self.OUTDIR_1D_ERR          =os.path.join(self.OUTDIR,"Obs_1D_err")      #! Display ii.
		self.OUTDIR_1D_ERR_VSL      =os.path.join(self.OUTDIR,"Obs_1D_err_vsl")  #! Display iii.
		self.OUTDIR_1D_TXT_RSLT     =os.path.join(self.OUTDIR,"Obs_1D_txt_rslt") #! Display iv.
		for d in [self.OUTDIR_1D, self.OUTDIR_1D_ERR, self.OUTDIR_1D_ERR_VSL, self.OUTDIR_1D_TXT_RSLT]:
			if not os.path.exists(d):
				os.makedirs(d)
		#! Make output root file
		self.FOUT_OBS_1D=ROOT.TFile("%s/Obs_1D.root"%self.OUTDIR,"RECREATE")

		#! Make OUTDIR_ITG,OUTDIR_ITG_ERR,OUTDIR_ITG_ERR_VSL,OUTDIR_ITG_TXT_RSTLS
		self.OUTDIR_ITG              =os.path.join(self.OUTDIR,"Obs_itg")          #! Display i.
		self.OUTDIR_ITG_ERR          =os.path.join(self.OUTDIR,"Obs_itg_err")      #! Display ii.
		self.OUTDIR_ITG_ERR_VSL      =os.path.join(self.OUTDIR,"Obs_itg_err_vsl")  #! Display iii.
		self.OUTDIR_ITG_TXT_RSLT     =os.path.join(self.OUTDIR,"Obs_itg_txt_rslt") #! Display iv.
		for d in [self.OUTDIR_ITG, self.OUTDIR_ITG_ERR, self.OUTDIR_ITG_ERR_VSL, self.OUTDIR_ITG_TXT_RSLT]:
			if not os.path.exists(d):
				os.makedirs(d)
		#! Make output root file
		self.FOUT_OBS_ITG=ROOT.TFile("%s/Obs_itg.root"%self.OUTDIR,"RECREATE")
		
		#! 1. Make displays for Obs_1D
		#! 1.i. Get Q2-W binning from file
		if self.DBG==True: q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX,dbg=True,dbg_bins=3)
		else:              q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX)
		print "DispObs::disp_SE_results_1D():Obs_1D Q2-W bins got from obs_1D.root=",q2wbinl

		#! 1.ii. Looper over Q2-W bins and get h1(obsnum,seq,vst,var)
		print "DispObs::disp_SE_results_1D():Obs_1D Going to begin processing Q2-W bins from file"
		for q2wbin in q2wbinl:
			print"DispObs::disp_SE_results_1D(): Processing q2wbin=",q2wbin

			#! 1.iii. Depending on SS_TYPE, Create h1(obsnum,seq,vst,var) or h1_obs(err)(obsnum,seq,vst,var)
			if self.SS_TYPE=="per_non_vst_SE":
				h1=OrderedDict()
				h1err=None 
			elif self.SS_TYPE=="cmb_non_vst_SE":
				h1=OrderedDict() 
				h1err=OrderedDict()

			for k in self.OBSD:
				obsnum,obs,obstag=k,self.OBSD[k][0],self.OBSD[k][1]
				print "DispObs::disp_SE_results_1D():Obs_1D Processing obsnum:obs:obstag=",obsnum,obs,obstag

				#! Open file for this obs
				obsdir=os.path.join(self.OBSDIR,obs)
				if   self.SS_TYPE=="per_non_vst_SE": f=root_open(os.path.join(obsdir,'Obs_1D_norm/obs_1D.root'))
				elif self.SS_TYPE=="cmb_non_vst_SE": f=root_open(os.path.join(obsdir,'Obs_1D.root'))
			
				for item in list(itertools.product(self.SEQS,VSTS)):
					seq,vst=item[0],item[1]
					if   vst==1:           varl=['M1','THETA','ALPHA']
					elif vst==2 or vst==3: varl=['M2','THETA','ALPHA']
					for var in varl:
						if self.SS_TYPE=="per_non_vst_SE":
							h1[obsnum,seq,vst,var]=f.Get("%s/h1_%s_%d_%s"%(q2wbin,seq,vst,var))

							#! + The following line was added to decouple histograms from the open
							#! ROOT file which can then be closed.
							#! + This had to be done because otherwise I was getting an error from
							#! rootpy saying that "Too many files are open" (perhaps related to 'ulimit': See notes in ZIM under Technology/UNIX):
							#! rootpy.ROOTError: level=5000, loc='TFile::TFile', msg='file /data/trivedia/e16/2pi/d2pi/highQ2_030916/cutsncors16/sim9/yield.root can not be opened for reading (Too many open files)'
							h1[obsnum,seq,vst,var].SetDirectory(0)

							#! Call Sumw2() since the errors as currently in the h1s
							#! need to be propagated in all subsequent calculations that involve them
							h1[obsnum,seq,vst,var].Sumw2()
						elif self.SS_TYPE=="cmb_non_vst_SE":
							h1[obsnum,seq,vst,var]    =f.Get("%s/h1_obs_%s_%d_%s"%(q2wbin,seq,vst,var))
							h1err[obsnum,seq,vst,var]=f.Get("%s/h1_err_%s_%d_%s"%(q2wbin,seq,vst,var))

							#! + The following line was added to decouple histograms from the open
							#! ROOT file which can then be closed.
							#! + This had to be done because otherwise I was getting an error from
							#! rootpy saying that "Too many files are open" (perhaps related to 'ulimit': See notes in ZIM under Technology/UNIX):
							#! rootpy.ROOTError: level=5000, loc='TFile::TFile', msg='file /data/trivedia/e16/2pi/d2pi/highQ2_030916/cutsncors16/sim9/yield.root can not be opened for reading (Too many open files)'
							h1[obsnum,seq,vst,var].SetDirectory(0)
							h1err[obsnum,seq,vst,var].SetDirectory(0)

							#! Call Sumw2() since the errors as currently in the h1s
							#! need to be propagated in all subsequent calculations that involve them
							h1[obsnum,seq,vst,var].Sumw2()
							h1err[obsnum,seq,vst,var].Sumw2()
				#! End of [seq,vst,var] loop

				#! Close file
				#! + Was getting intermittent error from rootpy saying "Too many files are open": 
				#! rootpy.ROOTError: level=5000, loc='TFile::TFile', msg='file /data/trivedia/e16/2pi/d2pi/highQ2_030916/cutsncors16/sim9/yield.root can not be opened for reading (Too many open files)'
				#f.Close() 
			#! End of OBSD loop

			#! 1.iv. Get results(seq,vst,var) from h1(obs,seq,vst,var)
			if   self.SS_TYPE=="per_non_vst_SE": rslts=self.get_per_non_vst_SE_results(h1)
			elif self.SS_TYPE=="cmb_non_vst_SE": rslts=self.get_cmb_non_vst_SE_results(h1,h1err)

			#! 1.v. Display: call plot_SE_results_1D(h1,rslts,q2wbin)
			self.plot_SE_results_1D(h1,rslts,q2wbin,h1err)
		#! End of q2wbin loop

		#! 2. Make displays for Obs_itg
		#! 2.i. Get Q2-W binning from file
		#! + Q2-W binning in the file is as per Q2-W binning of h8s
		if self.DBG==True: q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,1.400,2.125,dbg=True,dbg_bins=3,from_obs_itg_yld=True)
		else:              q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,1.400,2.125,from_obs_itg_yld=True)
		print "DispObs::disp_SE_results_1D():Obs_Itg_Yld Q2-W bins got from obs_itg_yld.root=",q2wbinl

		#! 2.ii. Create hW(q2wbin,obsnum,seq,vst,var) to store integrated yields
		if self.SS_TYPE=="per_non_vst_SE":
			hW=OrderedDict()
			hWerr=None
		elif self.SS_TYPE=="cmb_non_vst_SE":
			hW=OrderedDict()
			hWerr=OrderedDict()
		print "DispObs::disp_SE_results_1D():Obs_Itg_Yld Going to create hW(q2wbin,obsnum,seq,vst,var)"
		for q2wbin in q2wbinl:
			print"Processing q2wbin=",q2wbin

			for k in self.OBSD:
				obsnum,obs,obstag=k,self.OBSD[k][0],self.OBSD[k][1]
				print "DispObs::disp_SE_results_1D():Obs_Itg_Yld Processing obsnum:obs:obstag=",obsnum,obs,obstag

				#! Open file for this obs
				obsdir=os.path.join(self.OBSDIR,obs)
				if   self.SS_TYPE=="per_non_vst_SE":f=root_open(os.path.join(obsdir,'Obs_Itg_Yld_norm/obs_itg_yld.root'))
				elif self.SS_TYPE=="cmb_non_vst_SE":f=root_open(os.path.join(obsdir,'Obs_itg.root'))
			
				for item in list(itertools.product(self.SEQS,VSTS)):
					seq,vst=item[0],item[1]
					if   vst==1:           varl=['M1','THETA','ALPHA']
					elif vst==2 or vst==3: varl=['M2','THETA','ALPHA']
					for var in varl:
						if self.SS_TYPE=="per_non_vst_SE":
							hW[q2wbin,obsnum,seq,vst,var]=f.Get("%s/hW_%s_%d_%s"%(q2wbin,seq,vst,var))

							#! + The following line was added to decouple histograms from the open
							#! ROOT file which can then be closed.
							#! + This had to be done because otherwise I was getting an error from
							#! rootpy saying that "Too many files are open" (perhaps related to 'ulimit': See notes in ZIM under Technology/UNIX):
							#! rootpy.ROOTError: level=5000, loc='TFile::TFile', msg='file /data/trivedia/e16/2pi/d2pi/highQ2_030916/cutsncors16/sim9/yield.root can not be opened for reading (Too many open files)'
							hW[q2wbin,obsnum,seq,vst,var].SetDirectory(0)

							#! Call Sumw2() since the errors as currently in the h1s
							#! need to be propagated in all subsequent calculations that involve them
							hW[q2wbin,obsnum,seq,vst,var].Sumw2()
						elif self.SS_TYPE=="cmb_non_vst_SE":
							if var!='THETA': continue
							hW[q2wbin,obsnum,seq,vst,var]   =f.Get("%s/hW_obs_%s_%d_%s"%(q2wbin,seq,vst,var))
							hWerr[q2wbin,obsnum,seq,vst,var]=f.Get("%s/hW_err_%s_%d_%s"%(q2wbin,seq,vst,var))

							#! + The following line was added to decouple histograms from the open
							#! ROOT file which can then be closed.
							#! + This had to be done because otherwise I was getting an error from
							#! rootpy saying that "Too many files are open" (perhaps related to 'ulimit': See notes in ZIM under Technology/UNIX):
							#! rootpy.ROOTError: level=5000, loc='TFile::TFile', msg='file /data/trivedia/e16/2pi/d2pi/highQ2_030916/cutsncors16/sim9/yield.root can not be opened for reading (Too many open files)'
							hW[q2wbin,obsnum,seq,vst,var].SetDirectory(0)
							hWerr[q2wbin,obsnum,seq,vst,var].SetDirectory(0)

							#! Call Sumw2() since the errors as currently in the h1s
							#! need to be propagated in all subsequent calculations that involve them
							hW[q2wbin,obsnum,seq,vst,var].Sumw2()
							hWerr[q2wbin,obsnum,seq,vst,var].Sumw2()
				#! End of [seq,vst,var] loop

				#! Close file
				#! + Was getting intermittent error from rootpy saying "Too many files are open": 
				#! rootpy.ROOTError: level=5000, loc='TFile::TFile', msg='file /data/trivedia/e16/2pi/d2pi/highQ2_030916/cutsncors16/sim9/yield.root can not be opened for reading (Too many open files)'
				f.Close() 
		#! end q2wbinl loop
				
		#! 2.iii. Plot hW(q2wbin,obsnum,seq,vst,var)
		self.plot_SE_results_W(hW,hWerr)

		print "*** disp_SE_results_1D() Done ***\n"
		print "DispObs::disp_SE_results_1D() If not exiting immediately then Python is probably doing Garbage Collection which can take a while"

	def disp_cmb_vst_SE_results_1D(self):
		'''
		[09-07-16]
		+ This function displays results, which includes errors, for Obs:itg,1D based on all variations
		for "cmb_vst_SE"

		+ This function uses the following output .root file from observables in self.OBSD:
			+ For Obs_1D:   Obs_1D.root
		  	+ For Itg-xsec: Obs_itg.root

		+ For 'cmb_vst_SE', vst-SE are obtain from Obs_itg and therefore they are processed first, from where
		  errors are obtained, which are then used applied Obs_1D

		-1. First make sure SS_TYPE is valid

		0. Make all needed output directories

		1. Process Obs_itg and get vst_SE from them
		1.i.   Get Q2-W binning from file
		1.ii.  Create hW(q2wbin,obsnum,seq,vst,var) to store integrated yields
		1.iii. Loop through file and fill hW(q2wbin,obsnum,seq,vst,var)
			+ Only var='THETA' applicabale here
			+ obsnum=1 always here 
		1.iv. Plot cmb_vst_SE_results for Obs_itg: self.plot_cmb_vst_SE_results_hW(hW,hWerr)
		      NOTE that this functions also returns 'rslts' for cmb_vst_SE that are applied to Obs_1D,

		2. Make displays for Obs_1D
		2.i. Get Q2-W binning from file	
		2.ii. Loop over Q2-W bins and get h1/h1err(obsnum,seq,vst,var)
			+ obsnum=1 here always
		2.iii. Display: plot_cmb_vst_SE_results_1D(h1,rslts,q2wbin,h1err)
		'''
		print "*** In DispObs::disp_cmb_vst_SE_results_1D ***"

		#! -1. Make sure SS_TYPE is valid
		if self.SS_TYPE!="cmb_vst_SE":
			print "DispObs::disp_cmb_vst_SE_results_1D(): self.SS_TYPE=%s is not valid. Exiting."%self.SS_TYPE

		#! 0. Make output directories
		#! Make OUTDIR_1D,OUTDIR_1D_ERR,OUTDIR_1D_ERR_VSL,OUTDIR_1D_TXT_RSTLS
		self.OUTDIR_1D              =os.path.join(self.OUTDIR,"Obs_1D")          #! Display i.
		self.OUTDIR_1D_ERR          =os.path.join(self.OUTDIR,"Obs_1D_err")      #! Display ii.
		self.OUTDIR_1D_ERR_VSL      =os.path.join(self.OUTDIR,"Obs_1D_err_vsl")  #! Display iii.
		self.OUTDIR_1D_TXT_RSLT     =os.path.join(self.OUTDIR,"Obs_1D_txt_rslt") #! Display iv.
		for d in [self.OUTDIR_1D, self.OUTDIR_1D_ERR, self.OUTDIR_1D_ERR_VSL, self.OUTDIR_1D_TXT_RSLT]:
			if not os.path.exists(d):
				os.makedirs(d)
		#! Make output root file
		self.FOUT_OBS_1D=ROOT.TFile("%s/Obs_1D.root"%self.OUTDIR,"RECREATE")

		#! Make OUTDIR_ITG,OUTDIR_ITG_ERR,OUTDIR_ITG_ERR_VSL,OUTDIR_ITG_TXT_RSTLS
		self.OUTDIR_ITG              =os.path.join(self.OUTDIR,"Obs_itg")          #! Display i.
		self.OUTDIR_ITG_ERR          =os.path.join(self.OUTDIR,"Obs_itg_err")      #! Display ii.
		self.OUTDIR_ITG_ERR_VSL      =os.path.join(self.OUTDIR,"Obs_itg_err_vsl")  #! Display iii.
		self.OUTDIR_ITG_TXT_RSLT     =os.path.join(self.OUTDIR,"Obs_itg_txt_rslt") #! Display iv.
		for d in [self.OUTDIR_ITG, self.OUTDIR_ITG_ERR, self.OUTDIR_ITG_ERR_VSL, self.OUTDIR_ITG_TXT_RSLT]:
			if not os.path.exists(d):
				os.makedirs(d)
		#! Make output root file
		self.FOUT_OBS_ITG=ROOT.TFile("%s/Obs_itg.root"%self.OUTDIR,"RECREATE")
		
		#! 1. Process Obs_itg and get vst_SE from them
		#! 1.i. Get Q2-W binning from file
		#! + Q2-W binning in the file is as per Q2-W binning of h8s
		if self.DBG==True: q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,1.400,2.125,dbg=True,dbg_bins=3,from_obs_itg_yld=True)
		else:              q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,1.400,2.125,from_obs_itg_yld=True)
		print "DispObs::disp_cmb_vst_SE_results_1D():Obs_Itg_Yld Q2-W bins got from obs_itg_yld.root=",q2wbinl

		#! 1.ii. Create hW(q2wbin,obsnum,seq,vst,var) to store integrated yields
		hW=OrderedDict()
		hWerr=OrderedDict()
		print "DispObs::disp_cmb_vst_SE_results_1D():Obs_Itg_Yld Going to create hW(q2wbin,obsnum,seq,vst,var)"
		#! 1.iii. Loop through file and fill hW(q2wbin,obsnum,seq,vst,var)
		for q2wbin in q2wbinl:
			print"Processing q2wbin=",q2wbin

			for k in self.OBSD:
				obsnum,obs,obstag=k,self.OBSD[k][0],self.OBSD[k][1]
				print "DispObs::disp_cmb_vst_SE_results_1D():Obs_Itg_Yld Processing obsnum:obs:obstag=",obsnum,obs,obstag

				#! Open file for this obs
				obsdir=os.path.join(self.OBSDIR,obs)
				f=root_open(os.path.join(obsdir,'Obs_itg.root'))
			
				for item in list(itertools.product(self.SEQS,VSTS)):
					seq,vst=item[0],item[1]
					if   vst==1:           varl=['M1','THETA','ALPHA']
					elif vst==2 or vst==3: varl=['M2','THETA','ALPHA']
					for var in varl:
						if var!='THETA': continue
						hW[q2wbin,obsnum,seq,vst,var]   =f.Get("%s/hW_obs_%s_%d_%s"%(q2wbin,seq,vst,var))
						hWerr[q2wbin,obsnum,seq,vst,var]=f.Get("%s/hW_err_%s_%d_%s"%(q2wbin,seq,vst,var))

						#! + The following line was added to decouple histograms from the open
						#! ROOT file which can then be closed.
						#! + This had to be done because otherwise I was getting an error from
						#! rootpy saying that "Too many files are open" (perhaps related to 'ulimit': See notes in ZIM under Technology/UNIX):
						#! rootpy.ROOTError: level=5000, loc='TFile::TFile', msg='file /data/trivedia/e16/2pi/d2pi/highQ2_030916/cutsncors16/sim9/yield.root can not be opened for reading (Too many open files)'
						hW[q2wbin,obsnum,seq,vst,var].SetDirectory(0)
						hWerr[q2wbin,obsnum,seq,vst,var].SetDirectory(0)

						#! Call Sumw2() since the errors as currently in the h1s
						#! need to be propagated in all subsequent calculations that involve them
						hW[q2wbin,obsnum,seq,vst,var].Sumw2()
						hWerr[q2wbin,obsnum,seq,vst,var].Sumw2()
				#! End of [seq,vst,var] loop

				#! Close file
				#! + Was getting intermittent error from rootpy saying "Too many files are open": 
				#! rootpy.ROOTError: level=5000, loc='TFile::TFile', msg='file /data/trivedia/e16/2pi/d2pi/highQ2_030916/cutsncors16/sim9/yield.root can not be opened for reading (Too many open files)'
				f.Close() 
		#! end q2wbinl loop
				
		#! 1.iv. Plot cmb_vst_SE_results for Obs_itg: self.plot_cmb_vst_SE_results_hW(hW,hWerr)
		#!        NOTE that this functions also returns 'rslts' for cmb_vst_SE that are applied to Obs_1D, 
		#!        where 'rslts[q2bin,wbin,seq]' are the cmb_vst_SE 
		#!        (for more details see self.plot_cmb_vst_SE_results_hW(hW,hWerr))
		#!         
		rslts_itg=self.plot_cmb_vst_SE_results_hW(hW,hWerr)
		#! To print results for debug
		# print "results_itg from Obs_itg:"
		# for k in rslts_itg:
		# 	q2bin,wbin,seq=k[0],k[1],k[2]
		# 	mu,sg_mu,sg,sg_sg =rslts_itg[k][0],rslts_itg[k][1],rslts_itg[k][2],rslts_itg[k][3]
		# 	f_vst_1,sg_f_vst_1=rslts_itg[k][4],rslts_itg[k][5]
		# 	f_vst_2,sg_f_vst_2=rslts_itg[k][6],rslts_itg[k][7]
		# 	f_vst_3,sg_f_vst_3=rslts_itg[k][8],rslts_itg[k][9]
		# 	sg_vst_1D,sg_sg_vst_1D=rslts_itg[k][10],rslts_itg[k][11]
		# 	print "%s,%s,%s=%f,%f,%f,%f\n,%f,%f\n%f,%f\n %f,%f\n%f,%f"%(q2bin,wbin,seq,mu,sg_mu,sg,sg_sg,
		# 		                          f_vst_1,sg_f_vst_1,
		# 		                          f_vst_2,sg_f_vst_2,
		# 		                          f_vst_3,sg_f_vst_3,
		# 		                          sg_vst_1D,sg_sg_vst_1D)
		

		#! 2. Make displays for Obs_1D
		#! 2.i. Get Q2-W binning from file
		if self.DBG==True: q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX,dbg=True,dbg_bins=3)
		else:              q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX)
		print "DispObs::disp_cmb_vst_SE_results_1D():Obs_1D Q2-W bins got from obs_1D.root=",q2wbinl

		#! 2.ii. Loop over Q2-W bins and get h1/h1err(obsnum,seq,vst,var)
		print "DispObs::disp_cmb_vst_SE_results_1D():Obs_1D Going to begin processing Q2-W bins from file"
		for q2wbin in q2wbinl:
			#! get q2wbin information
			q2bin=self.get_q2bin(q2wbin)
			wbin=self.get_wbin(q2wbin)
			print"DispObs::disp_cmb_vst_SE_results_1D(): Processing q2wbin=",q2wbin

			#! Create h1(obsnum,seq,vst,var) or h1_obs(err)(obsnum,seq,vst,var)
			h1=OrderedDict() 
			h1err=OrderedDict()

			for k in self.OBSD:
				obsnum,obs,obstag=k,self.OBSD[k][0],self.OBSD[k][1]
				print "DispObs::disp_cmb_vst_SE_results_1D():Obs_1D Processing obsnum:obs:obstag=",obsnum,obs,obstag

				#! Open file for this obs
				obsdir=os.path.join(self.OBSDIR,obs)
				f=root_open(os.path.join(obsdir,'Obs_1D.root'))
			
				for item in list(itertools.product(self.SEQS,VSTS)):
					seq,vst=item[0],item[1]
					if   vst==1:           varl=['M1','THETA','ALPHA']
					elif vst==2 or vst==3: varl=['M2','THETA','ALPHA']
					for var in varl:
						h1[obsnum,seq,vst,var]    =f.Get("%s/h1_obs_%s_%d_%s"%(q2wbin,seq,vst,var))
						h1err[obsnum,seq,vst,var]=f.Get("%s/h1_err_%s_%d_%s"%(q2wbin,seq,vst,var))

						#! + The following line was added to decouple histograms from the open
						#! ROOT file which can then be closed.
						#! + This had to be done because otherwise I was getting an error from
						#! rootpy saying that "Too many files are open" (perhaps related to 'ulimit': See notes in ZIM under Technology/UNIX):
						#! rootpy.ROOTError: level=5000, loc='TFile::TFile', msg='file /data/trivedia/e16/2pi/d2pi/highQ2_030916/cutsncors16/sim9/yield.root can not be opened for reading (Too many open files)'
						h1[obsnum,seq,vst,var].SetDirectory(0)
						h1err[obsnum,seq,vst,var].SetDirectory(0)

						#! Call Sumw2() since the errors as currently in the h1s
						#! need to be propagated in all subsequent calculations that involve them
						h1[obsnum,seq,vst,var].Sumw2()
						h1err[obsnum,seq,vst,var].Sumw2()
				#! End of [seq,vst,var] loop

				#! Close file
				#! + Was getting intermittent error from rootpy saying "Too many files are open": 
				#! rootpy.ROOTError: level=5000, loc='TFile::TFile', msg='file /data/trivedia/e16/2pi/d2pi/highQ2_030916/cutsncors16/sim9/yield.root can not be opened for reading (Too many open files)'
				#f.Close() 
			#! End of OBSD loop

			#! 2.iii. Using 'rslts_itg', get cmb_vst_SE_results_1D(h1,rslts_q2wbin_EC,rslts_q2wbin_EF,q2wbin,h1err)
			#!        NOTE how only results appropriate to this bin are passed
			rslts=self.get_cmb_vst_SE_results_1D(h1,rslts_itg[q2bin,wbin,'EC'],rslts_itg[q2bin,wbin,'EF'],q2wbin,h1err)

			#! 2.iv. Display results
			self.plot_SE_results_1D(h1,rslts,q2wbin,h1err)
		#! End of q2wbin loop

		print "*** disp_cmb_vst_SE_results_1D() Done ***\n"
		print "DispObs::disp_cmb_vst_SE_results_1D If not exiting immediately then Python is probably doing Garbage Collection which can take a while"

	def disp_cmb_vst_SE_results_R2(self,extract_D_E_for_non_alpha):
		'''
		[09-13-16] 
		Since I did not find the process of obtain cmb_vst_SE results applicable for R2, this method simply links
		'cmb_non_vst_SE' R2 results in the 'cmb_vst_SE' directory
		be
		'''
		print "*** In DispObs::disp_cmb_vst_SE_results_R2 ***"

		self.EXTRACT_D_E_FOR_NON_ALPHA=extract_D_E_for_non_alpha
		if self.EXTRACT_D_E_FOR_NON_ALPHA!=True and self.EXTRACT_D_E_FOR_NON_ALPHA!=False:
			print "DispObs::disp_SE_results_R2(): self.EXTRACT_D_E_FOR_NON_ALPHA=%s is not valid. Exiting."%self.EXTRACT_D_E_FOR_NON_ALPHA
			sys.exit()

		#! -1. Make sure SS_TYPE is valid
		if self.SS_TYPE!="cmb_vst_SE":
			print "DispObs::disp_cmb_vst_SE_results_R2(): self.SS_TYPE=%s is not valid. Exiting."%self.SS_TYPE

		for k in self.OBSD:
			obsnum,obs,obstag=k,self.OBSD[k][0],self.OBSD[k][1]
			print "DispObs::disp_cmb_vst_SE_results_R2():Processing obsnum:obs:obstag=",obsnum,obs,obstag

			#! Open file for this obs
			obsdir=os.path.join(self.OBSDIR,obs)
			#f=root_open(os.path.join(obsdir,'Obs_itg.root'))

			if self.EXTRACT_D_E_FOR_NON_ALPHA==False:
				for rslt in ["Obs_R2", "Obs_R2_err", "Obs_R2_err_vsl", "Obs_R2_txt_rslt", "Obs_R2.root"]:
					cmd=["ln", "-s", "%s/%s"%(obsdir,rslt),"%s"%self.OUTDIR] 
					print">>>%s\n"%cmd
					#mainlog.write(">>>%s\n"%cmd)
					tool=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)#! stdout=logfile
					tool.wait()
			elif self.EXTRACT_D_E_FOR_NON_ALPHA==True: 
				cmd=["ln", "-s", "%s/%s"%(obsdir,"Obs_R2_w_non-alpha_DE"),"%s"%self.OUTDIR]
				print">>>%s\n"%cmd
				#mainlog.write(">>>%s\n"%cmd)
				tool=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)#! stdout=logfile
				tool.wait()

		print "*** disp_cmb_vst_SE_results_R2() Done ***\n"
		print "DispObs::disp_cmb_vst_SE_results_R2 If not exiting immediately then Python is probably doing Garbage Collection which can take a while"

	def get_per_non_vst_SE_results(self,h1):
		'''
		[08-24-16]
		+ General comments:
			+ 'h1' here refers to any of the 1 dim. histograms for Obs: itg,1D,R2
				+ if nkeys==4, then h1=hW[obsnum,seq,vst,var] or h1D[obsnum,seq,vst,var]
				+ if nkeys==5, then h1=hR2[obsnum,R2,seq,vst,var]
		
		+ Details:
			+ Using h1[keys] calculate (obsnum is the iterated index here and relate to numbered non_vst_SEs):
			 (for details on formulae used, see handwritten notes)
				1. (mu         +/- sg_mu)[keys]
				2. (sg         +/  sg_sg)[keys] 
				3. (rel_err[%] +/- sg_rel_err[%])[keys]
					+ where rel_err=sg/fabs(mu) (Absolute value is required for Obs_R2 which can be negative)

			  and return various display of these results, including the full array of results:
				1. h1_obs[keys]:     histogrammed (mu +/- sg_mu)[keys]
				2. h1_err[keys]:     histogrammed (sg +/- sg_sg)[keys] 
				3. h1_rel_err[keys]: histogrammed (rel_err +/- sg_rel_err)[keys]
				4. rslts[keys]:      results array in the following order (9 pieces of information):
					+ bin, bin_le, bin_ue, mu, sg_mu, sg, sg_sg, rel_err, sg_rel_err
		'''
		#! Create output object dicts
		h1_obs,h1_err,h1_rel_err,rslts=OrderedDict(),OrderedDict(),OrderedDict(),OrderedDict()		

		#! get nkeys to determine if h1 is h1D/hW or hR2 (hR2 has an additional key for R2)
		nkeys=len(h1.keys()[0])

		if nkeys==4: #! i.e. h1 is h1D/hW(obsnum,seq,vst,var)		
			for item in list(itertools.product(self.SEQS,VSTS)):
				seq,vst=item[0],item[1]
				if   vst==1:           varl=['M1','THETA','ALPHA']
				elif vst==2 or vst==3: varl=['M2','THETA','ALPHA']
				for var in varl:
					#! Create output-h1s[seq,vst,var] by Cloning and then resetting
					#! the contents of the h1 of 1st obsnum in self.OBSD i.e. h1[1,seq,vst,var]
					h1_obs[seq,vst,var]=h1[1,seq,vst,var].Clone()
					h1_err[seq,vst,var]=h1[1,seq,vst,var].Clone()
					h1_rel_err[seq,vst,var]=h1[1,seq,vst,var].Clone()
					for h in [h1_obs[seq,vst,var],h1_err[seq,vst,var],h1_rel_err[seq,vst,var]]:
						h.Reset()
						h.SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
						h.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
						h.SetMarkerStyle(ROOT.gROOT.ProcessLine("kFullDotLarge"))
					#! Create output-rslts
					rslts[seq,vst,var]=[]
					#! Now loop over all bins can calculate results using all obsnum
					nbins=h1_obs[seq,vst,var].GetNbinsX()
					for ibin in range(nbins):
						binc,binerr=[],[]
						for obsnum in self.OBSD:
							binc.append(h1[obsnum,seq,vst,var].GetBinContent(ibin+1))
							binerr.append(h1[obsnum,seq,vst,var].GetBinError(ibin+1))
						
						# if self.DBG:
						# 	print "**** DispObs::get_per_non_vst_SE_results():: binc,binerr for seq,vst,var,bin=",seq,vst,var,ibin+1,"***"
						# 	print "binc=",binc
						# 	print "binerr=",binerr
						# 	print "******"

						#! Now calculate results (see handwritten notes for formulae)
						#! mu,sg_mu 
						n=len(binc)
						mu=np.mean(binc)
						sg_mu=math.sqrt(sum([x**2 for x in binerr]))/n
						# if self.DBG:
						# 	if n==0: 
						# 		print "DispObs::get_per_non_vst_SE_results(): n=0 for seq,vst,var,bin=",seq,vst,var,ibin+1
						# 	print "**** DispObs::get_per_non_vst_SE_results():: binc,binerr for seq,vst,var,bin=",seq,vst,var,ibin+1,"***"
						# 	print "mu,sg_mu=",mu,sg_mu
						# 	print "******"
						#! sg,sg_sg
						sg=np.std(binc)
						if sg==0: #! This is true even in the analytical formulae used for cases when sg!=0
							sg_sg=0
						else:
							sg_sg_1=( 1/(sg*n) )
							sg_sg_2=math.sqrt(sum( (x-mu)**2 * (y**2+sg_mu**2) for x,y in zip(binc,binerr) ))
							sg_sg=sg_sg_1*sg_sg_2 
						#! rel_err,sg_rel_err
						if mu==0 or sg==0: 
							rel_err=0
							sg_rel_err=0
						else:
							rel_err=(sg/mu)*100
							sg_rel_err=(rel_err)*math.sqrt((sg_sg/sg)**2+(sg_mu/mu)**2)
						#! Store bin information (will be used for passing for full results rslts dict)
						#! Check to make sure h1 is not hW because for those binning has labels of form:[wmin-wmax]
						xtitle=h1_obs[seq,vst,var].GetXaxis().GetTitle()
						if "W" in xtitle:
							bin_label=h1_obs[seq,vst,var].GetXaxis().GetBinLabel(ibin+1)
							bin_le=float(bin_label.split(",")[0].strip("["))
							bin_ue=float(bin_label.split(",")[1].strip(")"))
						else:
							bin_le=h1_obs[seq,vst,var].GetBinLowEdge(ibin+1)
							bin_ue=h1_obs[seq,vst,var].GetBinLowEdge(ibin+1+1)

						#! Now fill results in output objects
						#! mu,sg_mu => h1_obs
						h1_obs[seq,vst,var].SetBinContent(ibin+1,mu)
						h1_obs[seq,vst,var].SetBinError(ibin+1,sg_mu)
						#! sg,sg_sg => h1_err
						h1_err[seq,vst,var].SetBinContent(ibin+1,sg)
						h1_err[seq,vst,var].SetBinError(ibin+1,sg_sg)
						#! rel_err,sg_rel_err => h1_rel_err
						h1_rel_err[seq,vst,var].SetBinContent(ibin+1,rel_err)
						h1_rel_err[seq,vst,var].SetBinError(ibin+1,sg_rel_err)
						#! [bin,bin_le,bin_ue,mu,sg_mu,sg,sg_sg,rel_err,sg_rel_err] => rslts
						rslts[seq,vst,var].append([ibin+1, bin_le, bin_ue, mu, sg_mu, sg, sg_sg, rel_err, sg_rel_err])
						
		elif nkeys==5: #! i.e. h1 is hR2(obsnum,R2,seq,vst,var)
			for item in list(itertools.product(R2S,self.SEQS,VSTS)):
				R2,seq,vst=item[0],item[1],item[2]
				for var in VARS:
					if var=='PHI': continue
					if self.EXTRACT_D_E_FOR_NON_ALPHA==False:
						if (R2=='D' or R2=='E') and var!= 'ALPHA':continue
					#! Create hR2_SS_err[R2,seq,vst,var] by Cloning and then resetting
					#! the contents of the hR2 of 1st SS-Obs in self.OBSD i.e. hR2[1,R2,seq,vst,var]
					h1_obs[R2,seq,vst,var]=h1[1,R2,seq,vst,var].Clone()
					h1_err[R2,seq,vst,var]=h1[1,R2,seq,vst,var].Clone()
					h1_rel_err[R2,seq,vst,var]=h1[1,R2,seq,vst,var].Clone()
					for h in [h1_obs[R2,seq,vst,var],h1_err[R2,seq,vst,var],h1_rel_err[R2,seq,vst,var]]:
						h.Reset()
						h.SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
						h.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
						h.SetMarkerStyle(ROOT.gROOT.ProcessLine("kFullDotLarge"))
					#! Create output-rslts
					rslts[R2,seq,vst,var]=[]
					#! Now loop over all bins can calculate SS err using all obsnum
					nbins=h1_obs[R2,seq,vst,var].GetNbinsX()
					for ibin in range(nbins):
						binc,binerr=[],[]
						for obsnum in self.OBSD:
							binc.append(h1[obsnum,R2,seq,vst,var].GetBinContent(ibin+1))
							binerr.append(h1[obsnum,R2,seq,vst,var].GetBinError(ibin+1))

						#! Now calculate results (see handwritten notes for formulae)
						#! mu,sg_mu 
						n=len(binc)
						mu=np.mean(binc)
						sg_mu=math.sqrt(sum([x**2 for x in binerr]))/n
						#! sg,sg_sg
						sg=np.std(binc)
						if sg==0: #! This is true even in the analytical formulae used for cases when sg!=0
							sg_sg=0
						else:
							sg_sg_1=( 1/(sg*n) )
							sg_sg_2=math.sqrt(sum( (x-mu)**2 * (y**2+sg_mu**2) for x,y in zip(binc,binerr) ))
							sg_sg=sg_sg_1*sg_sg_2 
						#! rel_err,sg_rel_err
						if mu==0 or sg==0: 
							rel_err=0
							sg_rel_err=0
						else:
							rel_err=(sg/math.fabs(mu))*100
							sg_rel_err=(rel_err)*math.sqrt((sg_sg/sg)**2+(sg_mu/mu)**2)
						#! Store bin information (will be used for passing for full results rslts dict)
						bin_le=h1_obs[R2,seq,vst,var].GetBinLowEdge(ibin+1)
						bin_ue=h1_obs[R2,seq,vst,var].GetBinLowEdge(ibin+1+1)

						#! Now fill results in output objects
						#! mu,sg_mu => h1_obs
						h1_obs[R2,seq,vst,var].SetBinContent(ibin+1,mu)
						h1_obs[R2,seq,vst,var].SetBinError(ibin+1,sg_mu)
						#! sg,sg_sg => h1_err
						h1_err[R2,seq,vst,var].SetBinContent(ibin+1,sg)
						h1_err[R2,seq,vst,var].SetBinError(ibin+1,sg_sg)
						#! rel_err,sg_rel_err => h1_rel_err
						h1_rel_err[R2,seq,vst,var].SetBinContent(ibin+1,rel_err)
						h1_rel_err[R2,seq,vst,var].SetBinError(ibin+1,sg_rel_err)
						#! [bin,bin_le,bin_ue,mu,sg_mu,sg,sg_sg,rel_err,sg_rel_err] => rslts
						rslts[R2,seq,vst,var].append([ibin+1, bin_le, bin_ue, mu, sg_mu, sg, sg_sg, rel_err, sg_rel_err])
		
		return h1_obs,h1_err,h1_rel_err,rslts

	def get_cmb_non_vst_SE_results(self,h1,h1err):
		'''
		[08-24-16]
		+ General comments:
			+ 'h1' and 'h1err' here refers to 1 dim. histograms for Obs: itg,1D,R2 containing:
				+ u,sg_mu(<xsec> and its error)
			    + sg,sg_sg(RMS_<xsec> and its error) 
			   respectively, obtained during "per_non_vst_SE".

				+ if nkeys==4, then h1=hW[obsnum,seq,vst,var] or h1D[obsnum,seq,vst,var]
				+ if nkeys==5, then h1=hR2[obsnum,R2,seq,vst,var]
		
		+ Details:
			+ Using h1[keys] and h1err[keys] calculate (obsnum is the iterated index here and relate to numbered cmb_non_vst_SE):
			 (for details on formulae used, see handwritten notes)
				1. (mu         +/- sg_mu)[keys]
				2. (sg         +/  sg_sg)[keys] 
				3. (rel_err[%] +/- sg_rel_err[%])[keys]
					+ where rel_err=sg/fabs(mu) (Absolute value is required for Obs_R2 which can be negative)

			  and return various display of these results, including the full array of results:
				1. h1_obs[keys]:     histogrammed (mu +/- sg_mu)[keys]
				2. h1_err[keys]:     histogrammed (sg +/- sg_sg)[keys] 
				3. h1_rel_err[keys]: histogrammed (rel_err +/- sg_rel_err)[keys]
				4. rslts[keys]:      results array in the following order (9 pieces of information):
					+ bin, bin_le, bin_ue, mu, sg_mu, sg, sg_sg, rel_err, sg_rel_err
		'''
		#! Create output object dicts
		h1_obs,h1_err,h1_rel_err,rslts=OrderedDict(),OrderedDict(),OrderedDict(),OrderedDict()		

		#! get nkeys to determine if h1 is h1D/hW or hR2 (hR2 has an additional key for R2)
		nkeys=len(h1.keys()[0])

		if nkeys==4: #! i.e. h1 is h1D/hW(obsnum,seq,vst,var)	
			#! determine if h1 is hW
			#! + needed because of h1 is hW, then when self.SS_TYPE=="cmb_non_vst_SE", only var='THETA' is applicable
			h1_is_hW=False
			if "W" in h1[1,'EF',1,'THETA'].GetXaxis().GetTitle():
				h1_is_hW=True

			for item in list(itertools.product(self.SEQS,VSTS)):
				seq,vst=item[0],item[1]
				if   vst==1:           varl=['M1','THETA','ALPHA']
				elif vst==2 or vst==3: varl=['M2','THETA','ALPHA']
				for var in varl:
					if self.SS_TYPE=="cmb_non_vst_SE" and h1_is_hW and var!='THETA': continue
					#! Create output-h1s[seq,vst,var] by Cloning and then resetting
					#! the contents of the h1 of 1st obsnum in self.OBSD i.e. h1[1,seq,vst,var]
					h1_obs[seq,vst,var]=h1[1,seq,vst,var].Clone()
					h1_err[seq,vst,var]=h1[1,seq,vst,var].Clone()
					h1_rel_err[seq,vst,var]=h1[1,seq,vst,var].Clone()
					for h in [h1_obs[seq,vst,var],h1_err[seq,vst,var],h1_rel_err[seq,vst,var]]:
						h.Reset()
						h.SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
						h.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
						h.SetMarkerStyle(ROOT.gROOT.ProcessLine("kFullDotLarge"))
					#! Create output-rslts
					rslts[seq,vst,var]=[]
					#! Now loop over all bins can calculate results using all obsnum
					nbins=h1_obs[seq,vst,var].GetNbinsX()
					for ibin in range(nbins):
						mu_l,sg_mu_l,sg_l,sg_sg_l=[],[],[],[]
						for obsnum in self.OBSD:
							mu_l.append(h1[obsnum,seq,vst,var].GetBinContent(ibin+1))
							sg_mu_l.append(h1[obsnum,seq,vst,var].GetBinError(ibin+1))
							sg_l.append(h1err[obsnum,seq,vst,var].GetBinContent(ibin+1))
							sg_sg_l.append(h1err[obsnum,seq,vst,var].GetBinError(ibin+1))
						
						# if self.DBG:
						# 	print "**** DispObs::get_per_non_vst_SE_results():: binc,binerr for seq,vst,var,bin=",seq,vst,var,ibin+1,"***"
						# 	print "binc=",binc
						# 	print "binerr=",binerr
						# 	print "******"

						#! Now calculate results (see handwritten notes for formulae)
						#! mu,sg_mu 
						n=len(mu_l)
						mu=np.mean(mu_l)
						sg_mu=math.sqrt(sum([x**2 for x in sg_mu_l]))/n
						# if self.DBG:
						# 	if n==0: 
						# 		print "DispObs::get_per_non_vst_SE_results(): n=0 for seq,vst,var,bin=",seq,vst,var,ibin+1
						# 	print "**** DispObs::get_per_non_vst_SE_results():: binc,binerr for seq,vst,var,bin=",seq,vst,var,ibin+1,"***"
						# 	print "mu,sg_mu=",mu,sg_mu
						# 	print "******"

						#! sg,sg_sg
						sg=math.sqrt(sum([x**2 for x in sg_l]))
						if sg==0: #! This is true even in the analytical formulae used for cases when sg!=0
							sg_sg=0
						else:
							sg_sg_1=( 1/sg )
							sg_sg_2=math.sqrt(sum( (x**2) * (y**2) for x,y in zip(sg_l,sg_sg_l) ))
							sg_sg=sg_sg_1*sg_sg_2 
						#! rel_err,sg_rel_err
						if mu==0 or sg==0: 
							rel_err=0
							sg_rel_err=0
						else:
							rel_err=(sg/mu)*100
							sg_rel_err=(rel_err)*math.sqrt((sg_sg/sg)**2+(sg_mu/mu)**2)
						#! Store bin information (will be used for passing for full results rslts dict)
						#! Check to make sure h1 is not hW because for those binning has labels of form:[wmin-wmax]
						xtitle=h1_obs[seq,vst,var].GetXaxis().GetTitle()
						if "W" in xtitle:
							bin_label=h1_obs[seq,vst,var].GetXaxis().GetBinLabel(ibin+1)
							bin_le=float(bin_label.split(",")[0].strip("["))
							bin_ue=float(bin_label.split(",")[1].strip(")"))
						else:
							bin_le=h1_obs[seq,vst,var].GetBinLowEdge(ibin+1)
							bin_ue=h1_obs[seq,vst,var].GetBinLowEdge(ibin+1+1)

						#! Now fill results in output objects
						#! mu,sg_mu => h1_obs
						h1_obs[seq,vst,var].SetBinContent(ibin+1,mu)
						h1_obs[seq,vst,var].SetBinError(ibin+1,sg_mu)
						#! sg,sg_sg => h1_err
						h1_err[seq,vst,var].SetBinContent(ibin+1,sg)
						h1_err[seq,vst,var].SetBinError(ibin+1,sg_sg)
						#! rel_err,sg_rel_err => h1_rel_err
						h1_rel_err[seq,vst,var].SetBinContent(ibin+1,rel_err)
						h1_rel_err[seq,vst,var].SetBinError(ibin+1,sg_rel_err)
						#! [bin,bin_le,bin_ue,mu,sg_mu,sg,sg_sg,rel_err,sg_rel_err] => rslts
						rslts[seq,vst,var].append([ibin+1, bin_le, bin_ue, mu, sg_mu, sg, sg_sg, rel_err, sg_rel_err])
						
		elif nkeys==5: #! i.e. h1 is hR2(obsnum,R2,seq,vst,var)
			for item in list(itertools.product(R2S,self.SEQS,VSTS)):
				R2,seq,vst=item[0],item[1],item[2]
				for var in VARS:
					if var=='PHI': continue

					if self.EXTRACT_D_E_FOR_NON_ALPHA==False:
						if (R2=='D' or R2=='E') and var!= 'ALPHA':continue

					#! Create hR2_SS_err[R2,seq,vst,var] by Cloning and then resetting
					#! the contents of the hR2 of 1st SS-Obs in self.OBSD i.e. hR2[1,R2,seq,vst,var]
					h1_obs[R2,seq,vst,var]=h1[1,R2,seq,vst,var].Clone()
					h1_err[R2,seq,vst,var]=h1[1,R2,seq,vst,var].Clone()
					h1_rel_err[R2,seq,vst,var]=h1[1,R2,seq,vst,var].Clone()
					for h in [h1_obs[R2,seq,vst,var],h1_err[R2,seq,vst,var],h1_rel_err[R2,seq,vst,var]]:
						h.Reset()
						h.SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
						h.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
						h.SetMarkerStyle(ROOT.gROOT.ProcessLine("kFullDotLarge"))
					#! Create output-rslts
					rslts[R2,seq,vst,var]=[]
					#! Now loop over all bins can calculate SS err using all obsnum
					nbins=h1_obs[R2,seq,vst,var].GetNbinsX()
					for ibin in range(nbins):
						mu_l,sg_mu_l,sg_l,sg_sg_l=[],[],[],[]
						for obsnum in self.OBSD:
							mu_l.append(h1[obsnum,R2,seq,vst,var].GetBinContent(ibin+1))
							sg_mu_l.append(h1[obsnum,R2,seq,vst,var].GetBinError(ibin+1))
							sg_l.append(h1err[obsnum,R2,seq,vst,var].GetBinContent(ibin+1))
							sg_sg_l.append(h1err[obsnum,R2,seq,vst,var].GetBinError(ibin+1))

						#! Now calculate results (see handwritten notes for formulae)
						#! mu,sg_mu 
						n=len(mu_l)
						mu=np.mean(mu_l)
						sg_mu=math.sqrt(sum([x**2 for x in sg_mu_l]))/n

						#! sg,sg_sg
						sg=math.sqrt(sum([x**2 for x in sg_l]))
						if sg==0: #! This is true even in the analytical formulae used for cases when sg!=0
							sg_sg=0
						else:
							sg_sg_1=( 1/sg )
							sg_sg_2=math.sqrt(sum( (x**2) * (y**2) for x,y in zip(sg_l,sg_sg_l) ))
							sg_sg=sg_sg_1*sg_sg_2 
						#! rel_err,sg_rel_err
						if mu==0 or sg==0: 
							rel_err=0
							sg_rel_err=0
						else:
							rel_err=(sg/math.fabs(mu))*100
							sg_rel_err=(rel_err)*math.sqrt((sg_sg/sg)**2+(sg_mu/mu)**2)
						#! Store bin information (will be used for passing for full results rslts dict)
						bin_le=h1_obs[R2,seq,vst,var].GetBinLowEdge(ibin+1)
						bin_ue=h1_obs[R2,seq,vst,var].GetBinLowEdge(ibin+1+1)

						#! Now fill results in output objects
						#! mu,sg_mu => h1_obs
						h1_obs[R2,seq,vst,var].SetBinContent(ibin+1,mu)
						h1_obs[R2,seq,vst,var].SetBinError(ibin+1,sg_mu)
						#! sg,sg_sg => h1_err
						h1_err[R2,seq,vst,var].SetBinContent(ibin+1,sg)
						h1_err[R2,seq,vst,var].SetBinError(ibin+1,sg_sg)
						#! rel_err,sg_rel_err => h1_rel_err
						h1_rel_err[R2,seq,vst,var].SetBinContent(ibin+1,rel_err)
						h1_rel_err[R2,seq,vst,var].SetBinError(ibin+1,sg_rel_err)
						#! [bin,bin_le,bin_ue,mu,sg_mu,sg,sg_sg,rel_err,sg_rel_err] => rslts
						rslts[R2,seq,vst,var].append([ibin+1, bin_le, bin_ue, mu, sg_mu, sg, sg_sg, rel_err, sg_rel_err])
		
		return h1_obs,h1_err,h1_rel_err,rslts


	def get_cmb_vst_SE_results_hW(self,h1,h1err):
		'''
		[09-07-16]
		+ get_cmb_non_vst_SE_results used as template and that will be reflected below

		+ General comments:
			+ 'h1' and 'h1err' here refers to 1 dim. histograms for Obs: itg containing:
				+ u,sg_mu(<xsec> and its error)
			    + sg,sg_sg(RMS_<xsec> and its error) 
			   respectively, obtained during "cmb_non_vst_SE".
			+ keys in h1,h1err:
				+h1=h1[obsnum,seq,vst,var] or h1err=h1err[obsnum,seq,vst,var]
					+ obsnum=1 always
					+ var='THETA' always
				

		+ Details:
			+ Using hW[keys] and hWerr[keys] calculate (vst is the iterated index here):
			 (for details on formulae used, see handwritten notes)
				1. (mu         +/- sg_mu)[keys]
				2. (sg         +/  sg_sg)[keys] 
				3. (rel_err[%] +/- sg_rel_err[%])[keys]
					+ where rel_err=sg/fabs(mu) (Absolute value is required for Obs_R2 which can be negative)
				4. (f_vst_i        +/- sg_f_vst_i)
					+ where f_vst_i=mu/mu_vst_i  i.e. the renormalization factor for Obs_1D within a VST
						+ where i =VST index
				5. (sg_vst_1D    +/- sg_sg_vst_1D) 
					+ These errors are used for propagation vst_SE for 1D

			  and return various display of these results, including the full array of results:
				1. h1_obs[keys]:     histogrammed (mu +/- sg_mu)[keys]
				2. h1_err[keys]:     histogrammed (sg +/- sg_sg)[keys] 
				3. h1_rel_err[keys]: histogrammed (rel_err +/- sg_rel_err)[keys]
				4. rslts[keys]:      results array in the following order (15=(9+6) pieces of information):
					+ bin, bin_le, bin_ue, mu, sg_mu, sg, sg_sg, rel_err, sg_rel_err,
					  f_vst_1,sg_f_vst_1,
					  f_vst_2,sg_f_vst_2,
					  f_vst_3,sg_f_vst_3,
					  sg_vst_1D, sg_sg_vst_1D
		'''
		#! Create output object dicts
		h1_obs,h1_err,h1_rel_err,rslts=OrderedDict(),OrderedDict(),OrderedDict(),OrderedDict()		

		for seq in self.SEQS:
			#! Create output-h1s[seq,'THETA'] by Cloning and then resetting
			#! the contents of the h1[1,seq,1,'THETA']
			h1_obs[seq,'THETA']=h1[1,seq,1,'THETA'].Clone()
			h1_err[seq,'THETA']=h1[1,seq,1,'THETA'].Clone()
			h1_rel_err[seq,'THETA']=h1[1,seq,1,'THETA'].Clone()
			for h in [h1_obs[seq,'THETA'],h1_err[seq,'THETA'],h1_rel_err[seq,'THETA']]:
				h.Reset()
				h.SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
				h.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
				h.SetMarkerStyle(ROOT.gROOT.ProcessLine("kFullDotLarge"))
			#! Create output-rslts
			rslts[seq,'THETA']=[]
			#! Now loop over all bins can calculate results using all obsnum
			nbins=h1_obs[seq,'THETA'].GetNbinsX()
			for ibin in range(nbins):
				mu_l,sg_mu_l,sg_l,sg_sg_l=[],[],[],[]
				for vst in VSTS:
					mu_l.append(h1[1,seq,vst,'THETA'].GetBinContent(ibin+1))
					sg_mu_l.append(h1[1,seq,vst,'THETA'].GetBinError(ibin+1))
					sg_l.append(h1err[1,seq,vst,'THETA'].GetBinContent(ibin+1))
					sg_sg_l.append(h1err[1,seq,vst,'THETA'].GetBinError(ibin+1))
					
				# if self.DBG:
				# 	print "**** DispObs::get_per_non_vst_SE_results():: binc,binerr for seq,vst,var,bin=",seq,vst,var,ibin+1,"***"
				# 	print "binc=",binc
				# 	print "binerr=",binerr
				# 	print "******"

				#! Now calculate results (see handwritten notes for formulae)
				#! mu,sg_mu 
				n=len(mu_l)
				mu=np.mean(mu_l)
				sg_mu=math.sqrt(sum([x**2 for x in sg_mu_l]))/n
				# if self.DBG:
				# 	if n==0: 
				# 		print "DispObs::get_per_non_vst_SE_results(): n=0 for seq,vst,var,bin=",seq,vst,var,ibin+1
				# 	print "**** DispObs::get_per_non_vst_SE_results():: binc,binerr for seq,vst,var,bin=",seq,vst,var,ibin+1,"***"
				# 	print "mu,sg_mu=",mu,sg_mu
				# 	print "******"

				#! sg,sg_sg
				sg=math.sqrt(sum([x**2 for x in sg_l])) #! TEST #! np.std(mu_l)
				if sg==0: #! This is true even in the analytical formulae used for cases when sg!=0
					sg_sg=0
				else: #! TEST
					sg_sg_1=( 1/sg )
					sg_sg_2=math.sqrt(sum( (x**2) * (y**2) for x,y in zip(sg_l,sg_sg_l) ))
					sg_sg=sg_sg_1*sg_sg_2 

				#! rel_err,sg_rel_err
				if mu==0 or sg==0: 
					rel_err=0
					sg_rel_err=0
				else:
					rel_err=(sg/mu)*100
					sg_rel_err=(rel_err)*math.sqrt((sg_sg/sg)**2+(sg_mu/mu)**2)

				#! Calculate Obs_1D vst renormalization factors and their errors
				f_vst=[0,0,0]
				sg_f_vst=[0,0,0]
				for i,vst in enumerate(VSTS):
					f_vst[i]=mu/mu_l[i]
					sg_f_vst[i]=f_vst[i]*math.sqrt( (sg_mu/mu)**2 + (sg_mu_l[i]/mu_l[i])**2 )

				#! Calculate sg_vst_1D and sg_sg_vst_1D which is nothing but RMS of mu_l and its prop. error
				std_dev=np.std(mu_l)
				if std_dev==0: #! This is true even in the analytical formulae used for cases when sg!=0
					sg_std_dev=0
				else:
					sg_std_dev_1=( 1/(std_dev*n) )
					sg_std_dev_2=math.sqrt(sum( (x-mu)**2 * (y**2+sg_mu**2) for x,y in zip(mu_l,sg_mu_l) ))
					sg_std_dev=sg_std_dev_1*sg_std_dev_2 
				sg_vst_1D=std_dev
				sg_sg_vst_1D=sg_std_dev


				#! Store bin information (will be used for passing for full results rslts dict)
				bin_label=h1_obs[seq,'THETA'].GetXaxis().GetBinLabel(ibin+1)
				bin_le=float(bin_label.split(",")[0].strip("["))
				bin_ue=float(bin_label.split(",")[1].strip(")"))
				
				#! Now fill results in output objects
				#! mu,sg_mu => h1_obs
				h1_obs[seq,'THETA'].SetBinContent(ibin+1,mu)
				h1_obs[seq,'THETA'].SetBinError(ibin+1,sg_mu)
				#! sg,sg_sg => h1_err
				h1_err[seq,'THETA'].SetBinContent(ibin+1,sg)
				h1_err[seq,'THETA'].SetBinError(ibin+1,sg_sg)
				#! rel_err,sg_rel_err => h1_rel_err
				h1_rel_err[seq,'THETA'].SetBinContent(ibin+1,rel_err)
				h1_rel_err[seq,'THETA'].SetBinError(ibin+1,sg_rel_err)
				rslts[seq,'THETA'].append([ibin+1, bin_le, bin_ue, mu, sg_mu, sg, sg_sg, rel_err, sg_rel_err,
					                       f_vst[0],sg_f_vst[0],
					                       f_vst[1],sg_f_vst[1],
					                       f_vst[2],sg_f_vst[2],
					                       sg_vst_1D,sg_sg_vst_1D])
		return h1_obs,h1_err,h1_rel_err,rslts

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
				if self.DO_ETGT_BG_SUB_CORR:
					normf=normf/self.ETGT_BG_SUB_CORR_FCTR
				h1n[k].SetBinContent(ibin+1,normf)
				h1n[k].SetBinError(ibin+1,0)

		#! Now normalize h1s
		for k in h1:
			h1[k].Divide(h1n[k])

		print "*** norm_1D() Done\n ***"

	def apply_rad_eff_corr(self,h1,q2wbin):
		print "*** In DispObs::apply_rad_eff_corr() ***"
		
		#! Get rad-eff-corr for this q2wbin
		cf,cferr=None,None
		if q2wbin in self.RAD_EFF_CORR_FCTR.keys():
			cf   =self.RAD_EFF_CORR_FCTR[q2wbin][0]
			cferr=self.RAD_EFF_CORR_FCTR[q2wbin][1]
		else:
			print "DispObs::apply_rad_eff_corr():WARNING: %s not found in self.RAD_EFF_CORR_FCTR. Using cf,cfferr=1,1"
			cf=1
			cferr=1

		#! Now apply cf,cferr to h1s
		#! + In order to propagate errors, use TH1::Multiply()
		for k in h1:
			#! 1. First make hm(=Multiply by histogram) all of whose binc,bine=cf,cferr
			#! + Do this using TH1::Clone()
			hm=h1[k].Clone("hm")
			hm.SetTitle("hm")
			#! Reset bin contents i.e. set binc,bine=0
			hm.Reset()
			#! Now set all binc,bine=cf,cferr
			for ibin in range(hm.GetNbinsX()):
				hm.SetBinContent(ibin+1,cf)
				hm.SetBinError(ibin+1,cferr)
			#! Set Sumw2() for hm so that errors are propagated
			hm.Sumw2()

			#! 2. Now apply rad-eff-corr
			h1[k].Multiply(hm)

		print "*** DispObs::apply_rad_eff_corr() Done\n ***"

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
		#! + The idea is to use a unique MarkerColor for [vst,var]
		#! + As a further distinguishing feature, since colors take time to distinguish,
		#!   each var is assigned its own MarkerStyle, which can be additionally be classified
		#!   as 'Full' or 'Open' type markers, when for example, more than one seq is being plotted on the same plot.  
		#! + [04-11-16] These aesthetics are optimized with reference to self.VIEW=="norm", where
		#!   all for each of 'EC' and 'EH', xsecs from all vst,var are plotted on the same canvas
		clrd={(1,'M1'):'kBlue',     (2,'M2'):'kBlue+2',     (3,'M2'):'kCyan',
			  (1,'THETA'):'kRed',   (2,'THETA'):'kRed+3',   (3,'THETA'):'kMagenta',
			  (1,'ALPHA'):'kGreen', (2,'ALPHA'):'kGreen+4', (3,'ALPHA'):'kOrange'}	
		mrkrd=     {'M1':'kFullTriangleUp', 'M2':'kFullTriangleUp', 'THETA':'kFullCircle', 'ALPHA':'kOpenSquare'}
		mrkr_fulld={'M1':'kFullDotLarge',   'M2':'kFullDotLarge',   'THETA':'kFullSquare', 'ALPHA':'kFullTriangleUp'}
		mrkr_opend={'M1':'kOpenCircle',     'M2':'kOpenCircle',     'THETA':'kOpenSquare', 'ALPHA':'kOpenTriangleUp'}

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
			c.SaveAs("%s/c_%s.eps"%(outdir_h2,h2[k].GetName()))

		#! 2. Plot projection h2 on W bins i.e. int_yld(W) for various Q2 bins
		outdir_w_proj=os.path.join(self.OUTDIR_ITG_YLD,"W")
		if not os.path.exists(outdir_w_proj):
			os.makedirs(outdir_w_proj)
		nq2bins=h2[h2.keys()[0]].GetNbinsY()
		for iq2bin in range(nq2bins):
			q2bin=h2[h2.keys()[0]].GetYaxis().GetBinLabel(iq2bin+1)
			#! Tranform q2bin (obtained from h2) to q2wbin, which is in the format used to name folders in .root files
			#! + q2bin=[q2min,q2max)
			#! + q2wbin=q2min-q2max_wmin-wmax (wmin,wmax here is equal to 1.400, 2.125)
			q2wbin=q2bin
			q2wbin=q2wbin.replace("[","")
			q2wbin=q2wbin.replace(")","")
			q2wbin=q2wbin.replace(",","-")
			q2wbin+="_1.400-2.125"
			#! For saving output to .root file, create and cd into q2wbindir
			if self.FOUT_ITG_YLD!=None:
				q2wbin_dir=self.FOUT_ITG_YLD.GetDirectory(q2wbin)
				if q2wbin_dir==None: 
					q2wbin_dir=self.FOUT_ITG_YLD.mkdir(q2wbin)
				q2wbin_dir.cd()

			if self.VIEW=="norm":
				#! Plotting aesthetics
				#ROOT.gROOT.SetOptStat("ne")
				
				#! + First create all projections: h1p[seq,vst,var]
				#! + This has to be done so that, before plotting, min and max for the Y-axis
				#!   can be determined from all the h1p[seq,vst,var]
				h1p=OrderedDict()
				for k in h2:
					seq,vst,var=k[0],k[1],k[2]
					hname="hW_%s_%d_%s"%(seq,vst,var)
					htitle="Integrated cross-section for Q2=%s"%q2bin#!h2[k].GetYaxis().GetBinLabel(iq2bin+1)
					h1p[k]=h2[k].ProjectionX(hname,iq2bin+1,iq2bin+1,"e")
					h1p[k].SetTitle(htitle)
					h1p[k].SetYTitle("#sigma[#mub]")
					h1p[k].SetMinimum(0)

				#! Determine max and set Y-axis min(=0),max for all h1p['EC',vst,var] and h1p['EF',vst,var]
				for seq in ['EC','EF']:
					maxl=[h1p[k].GetMaximum() for k in h2 if k[0]==seq]
					maximum=max(maxl)
					for k in h1p:
						h1p[k].SetMinimum(0.)
						h1p[k].SetMaximum(maximum+(.10*maximum))
				
				#! Now draw
				ROOT.gStyle.SetOptStat(0)
				#! + separate canvases for EC and EF
				for seq in ['EC','EF']:
					c=ROOT.TCanvas("c","c",1200,800)
					#c.Divide(2,1) #! One for the plot and the other for the large legend
					pad_l=ROOT.TPad("pad_l","Title pad",0.80,0.00,1.00,1.00)
					pad_p=ROOT.TPad("pad_p","Plots pad",0.00,0.00,0.80,1.00)
					pad_l.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
					pad_p.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
					pad_l.Draw()
					pad_p.Draw()
					#! legend
					l=ROOT.TLegend(0.0,0.50,1.00,1.00)
					l.SetFillStyle(0)
					l.SetBorderSize(0)
					l.SetTextSize(0.06)#0.02
					pad_p.cd()
					i=0
					for k in h1p:
						seqq,vst,var=k[0],k[1],k[2]
						if seqq!=seq: continue
						#! draw options
						draw_opt=""
						if i>0: draw_opt="same"
						#! Set up some aesthetics particular to hW
						#! + Rotate bin labels by 90 deg. and give them more space
						pad_p.SetBottomMargin(0.20)
						h1p[k].GetXaxis().LabelsOption("v")
						#! Setup markers and their colors
						clr=clrd[vst,var]
						mrkr=mrkrd[var]

						h1p[k].SetMarkerStyle(ROOT.gROOT.ProcessLine(mrkr))
						h1p[k].SetLineColor(ROOT.gROOT.ProcessLine(clr))
						h1p[k].SetMarkerColor(ROOT.gROOT.ProcessLine(clr))
						#! Draw
						h1p[k].Draw(draw_opt)
						i+=1
						if self.FOUT_ITG_YLD!=None:
							h1p[k].Write()
						l.AddEntry(h1p[k],"%s_VST%d_%s"%(seq,vst,var),"p")
					pad_l.cd()
					l.Draw()
					#! Save canvas
					c.SaveAs("%s/c_qbin%d_%s.png"%(outdir_w_proj,iq2bin+1,seq))
					c.SaveAs("%s/c_qbin%d_%s.eps"%(outdir_w_proj,iq2bin+1,seq))
			if self.VIEW=="norm_EC_EF_SF":
				#! Have not implemented plotting plotting for this view because not only
				#! is this view used mainly for comparing Obs_1D, but since in there intgrl-SF
				#! is normalized to intgrl-EF, here SF=EF and therefore not show anything new
				print("DispObs::plot_h2_itg_yld(): self.VIEW=norm_EC_EF_SF. Plotting not implented.")
			elif self.VIEW=="fullana": 
				#! Plotting aesthetics
				#ROOT.gROOT.SetOptStat("ne")
				#! + First create all projections: h1p[seq,vst,var]
				h1p=OrderedDict()
				for k in h2:
					seq,vst,var=k[0],k[1],k[2]
					if seq=='SH' or seq=='EH': continue
					hname="qbin_%d_%s_VST%d_%s"%(iq2bin+1,seq,vst,var)
					htitle="Integrated yield for Q2=%s"%h2[k].GetYaxis().GetBinLabel(iq2bin+1)
					h1p[k]=h2[k].ProjectionX(hname,iq2bin+1,iq2bin+1,"e")
					h1p[k].SetTitle(htitle)
					
					clr=clrd[vst,var]
					#! Set appropriate MarkerStyles for SEQS plotted on same cavas
					#! ER,SR on same
					if   seq=='ER': mrkr=mrkr_opend[var]
					elif seq=='SR': mrkr=mrkr_fulld[var]
					#! EC,EF on same; SC,SF,ST on same	
					if   seq=='EC' or seq=='SC':mrkr=mrkr_fulld[var]
					elif seq=='SF' or seq=='EF':mrkr=mrkr_opend[var]
					elif seq=='ST':				mrkr="kPlus"
					
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
				c1.SaveAs("%s/c_qbin%d_ER_SR.eps"%(outdir_w_proj,iq2bin+1))

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
				c2.SaveAs("%s/c_qbin%d_ST_SC_SF.eps"%(outdir_w_proj,iq2bin+1))

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
				c3.SaveAs("%s/c_qbin%d_EC_EH.eps"%(outdir_w_proj,iq2bin+1))

			elif self.VIEW=="ERyield":
				#! Plotting aesthetics
				#ROOT.gROOT.SetOptStat("ne")
				
				#! + First create all projections: h1p[seq,vst,var]
				#! + This has to be done so that, before plotting, min and max for the Y-axis
				#!   can be determined from all the h1p[seq,vst,var]
				h1p=OrderedDict()
				for k in h2:
					seq,vst,var=k[0],k[1],k[2]
					hname="hW_%s_%d_%s"%(seq,vst,var)
					htitle="itg_yld Q2=%s"%q2bin#!h2[k].GetYaxis().GetBinLabel(iq2bin+1)
					h1p[k]=h2[k].ProjectionX(hname,iq2bin+1,iq2bin+1,"e")
					h1p[k].SetTitle(htitle)
					h1p[k].SetYTitle("counts")
				
				#! Determine max and set Y-axis min(=0),max for all h1p[seq,vst,var]
				maxl=[h1p[k].GetMaximum() for k in h2]
				maximum=max(maxl)
				for k in h1p:
					h1p[k].SetMinimum(0.)
					h1p[k].SetMaximum(maximum)
				
				#! Now draw
				#! + If view=="norm" then separate canvases for EC and EH
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
				#! Save canvas
				c.SaveAs("%s/c_qbin%d.png"%(outdir_w_proj,iq2bin+1))
				c.SaveAs("%s/c_qbin%d.eps"%(outdir_w_proj,iq2bin+1))

		print "*** plot_h2_itg_yld() Done ***\n"

	def plot_h2_itg_yld_all_Q2_same_canvas(self,h2):
		'''
		[08-09-16]
		+ This function "simplifies"(simplifications described below) plotting of itg xsec as done in plot_h2_itg_yld(), 
		  where all possible results from seq,vst,var are plotted separately, and that too in different canvases for each Q2

		+ Simplications:
			+ Level 1: Plot only for seq,var,<vst>='EF','THETA',<VST1,VST2,VST3> because:
				+ 'EF' is used in the final result
				+  From the results of plot_h2_itg_yld() it can be seen that all VARs within a VST give same result. The only
				   variability is from different VSTs, hence the vst=average of all VSTs
			  NOTE, <vst> uses the error from TH1::Add() that used to obtain <vst>=(vst1+vst2+vst3)/3 

			+ Level 2: Plot all Q2 on the same TCanvas because:
				+ The Q2 evolution can directly be seen
		
		'''
		print "*** In DispObs::plot_h2_itg_yld_all_Q2_same_canvas ***"

		#! Define colors to differentiate between hW for different Q2
		clrd={1:'kBlue', 2:'kCyan', 3:'kRed', 4:'kMagenta', 5:'kGreen', 6:'kGreen+4', 7:'kOrange'}

		#! Plot projection h2 on W bins i.e. int_yld(W) for various Q2 bins on same canvas using:
		#! + seq,var,vst='EF','THETA',<VST1,VST2,VST3>
		#! Define TCanvas and TLegend
		c=ROOT.TCanvas("c","c",1200,800)
		l=ROOT.TLegend(0.8,0.80,1.00,1.00)
		l.SetFillStyle(0)
		l.SetBorderSize(0)
		#l.SetTextSize(0.06)#0.02
		#! Create output dir
		outdir_w_proj=os.path.join(self.OUTDIR_ITG_YLD,"W")
		if not os.path.exists(outdir_w_proj):
			os.makedirs(outdir_w_proj)
		#! Obtain number of q2 bins in h2
		nq2bins=h2[h2.keys()[0]].GetNbinsY()
		#! Create structure to store vst averaged histograms for each Q2
		h1p_vst_av=[0 for i in range(nq2bins)]
		#! Now loop over q2bins in h2
		for iq2bin in range(nq2bins):
			q2bin=h2[h2.keys()[0]].GetYaxis().GetBinLabel(iq2bin+1)
						
			#! Plotting aesthetics
			#ROOT.gROOT.SetOptStat("ne")
				
			#! + First create all projections: h1p['EF',vst,'THETA']
			#! + NOTE, call Sumw2() so that errors can be stored for later propagation
			h1p=OrderedDict()
			for k in h2:
				seq,vst,var=k[0],k[1],k[2]
				if seq=='EC': continue
				if var!='THETA': continue
				hname="hW_%s_%d_%s"%(seq,vst,var)
				htitle=""#! Remove title since legend will contain it
				h1p[k]=h2[k].ProjectionX(hname,iq2bin+1,iq2bin+1,"e")
				h1p[k].SetTitle(htitle)
				h1p[k].SetYTitle("#sigma[#mub]")
				h1p[k].SetMinimum(0)
				h1p[k].Sumw2()

			#! Now obtain h1p_vst_av=(h1p['EF','1','THETA']+h1p['EF','1','THETA']+h1p['EF','1','THETA'])/3
			#! 1. Create h1p_av using by Cloning h1p['EF','1','THETA']
			#! NOTE, set Sumw2() so that errors can be stored and then propagated
			h1p_vst_av[iq2bin]=h1p['EF',1,'THETA'].Clone("h1p_vst_av")
			h1p_vst_av[iq2bin].Sumw2()
			#! 2. Now Add to it h1p['EF',2,'THETA'] and h1p['EF',3,'THETA']
			h1p_vst_av[iq2bin].Add(h1p['EF',2,'THETA'])
			h1p_vst_av[iq2bin].Add(h1p['EF',3,'THETA'])
			#! 3. Divide by 3
			h1p_vst_av[iq2bin].Scale(1/3)

			#! SetMinimum=0 (for aesthetic reasons)
			h1p_vst_av[iq2bin].SetMinimum(0)

			#! Now draw
			clr=clrd[iq2bin+1]
			h1p_vst_av[iq2bin].SetMarkerStyle(ROOT.gROOT.ProcessLine('kFullCircle'))
			h1p_vst_av[iq2bin].SetLineColor(ROOT.gROOT.ProcessLine(clr))
			h1p_vst_av[iq2bin].SetMarkerColor(ROOT.gROOT.ProcessLine(clr))
			draw_opt=""
			if iq2bin>0: draw_opt="same"
			h1p_vst_av[iq2bin].Draw(draw_opt)
			l.AddEntry(h1p_vst_av[iq2bin],"Q2=%s"%q2bin,"p")

		#! Finally draw legend and save canvas
		l.Draw()
		#! Save canvas
		c.SaveAs("%s/c_all-qbin_vst-av_var-THETA.png"%(outdir_w_proj))	
		c.SaveAs("%s/c_all-qbin_vst-av_var-THETA.eps"%(outdir_w_proj))		

		print "*** Done DispObs::plot_h2_itg_yld_all_Q2_same_canvas ***"

	def plot_SE_results_W(self,hW,hWerr):
		'''
		[09-04-16]
		+ 'hWerr' added to accomodate plotting of "cmb_non_vst_SE"

		The following are the ways in which the results are displayed:
			1. EH: display all 'results' for 'EH'
			2. EC: write only text part of 'results'
			3. EC_EH comparison: Display only 'h1_obs'*['EC','EH']
			4. Results stored in FOUT: [EC,EH]* [hW_obs,hW_err,hW_rel_err] 

		Process:
		1. Get q2wbinl from hW(q2wbin,obsnum,seq,vst,var)
		2. Loop over the q2wbins and for each bin:
			2.i.   + Collect all hW in _hW(obsnum,seq,vst,var)
			       + also collect, if self.SSTYPE=="cmb_non_vst_SE", all hWerr in _hWerr(obsnum,seq,vst,var)
			2.ii.  For a particular seq,vst,var, using all obsnum, obtain 'results'(seq,vst,var) using:
				+ if self.SS_TYPE=="per_non_vst_SE": get_per_non_vst_SE_results()
				+ if self.SS_TYPE=="cmb_non_vst_SE": get_cmb_non_vst_SE_results()
				 where 'results':
				+ hW_obs[seq,vst,var]
			    + hW_err[seq,vst,var]
			    + hW_rel_err[seq,vst,var]
			    + rslts[seq,vst,var]
			2.iii. Collect all hW_obs in q2-dict: This will be used later to plot all xsec(W) for all Q2 in one canvas (step 3.)	 
			2.iv. For every EF*[vst]*THETA display:
				   (+ Note that full display i.e. [seq]*[vst]*[var] should be done individual SE results: $OBSDIR_E16/SE/cutsncorsX
				   	+ From there I have already learned that there is only variations due to vst and not in vars within a vst)
				di.   Final result(mu,sg_mu,sg,sg_sg_sg) from SE variations: hW_obs[seq,vst,var] and hW_err[seq,vst,var] together
				dii.  To just see relative error (rel_err, sg_rel_err):      hW_rel_err[seq,vst,var]
				diii. To visually see SE variations:                         _hW(obsnum,seq,vst,var) directly
				For every [EC,EF]*[vst]*THETA display:
				div.  Text file of full results:                             put rslts[seq,vst,var] in a text file
			2.v. EC-EH comparison:  For every [EC,EH]*[vst]*THETA display:
				i.   Final result(mu,sg_mu,sg,sg_sg_sg) from SE variations: hW_obs[seq,vst,var]
		3. Plot, per vst, xsec(W;EF,THETA) for all Q2 in one canvas.
				
		'''
		print "*** In DispObs::plot_SE_results_W() ***"
		#! Plotting aesthetics
		ROOT.gStyle.Reset() #! Reset aesthetics which may have been set by some other plot function
		ROOT.gStyle.SetOptStat("neiuo")
		ROOT.gStyle.SetErrorX(0.001)

		#! 1. Get list of q2wbins in hW(q2wbin,obsnum,seq,vst,var)
		#! Get all q2wbins
		q2wbinl=[k[0] for k in hW.keys()]
		#! Keep only unique bins
		q2wbinl=list(set(q2wbinl))
		print "DispObs::plot_SE_results_W() q2wbinl in hW=",q2wbinl

		#! Create a dictionary to keep all hW_obs labelled by Q2 (see steps 2.iii. and 3.)
		hW_obs_Q2=OrderedDict()
		#! 2. Loop over the q2wbins
		for iq2wbin,q2wbin_sel in enumerate(q2wbinl):
			#! 2.i. Collect all hW for a particular q2wbin in _hW(obsnum,seq,vst,var)
			_hW=OrderedDict()
			if self.SS_TYPE=="cmb_non_vst_SE":
				_hWerr=OrderedDict()
			for k in hW:
				q2wbin,obsnum,seq,vst,var=k[0],k[1],k[2],k[3],k[4]
				if q2wbin==q2wbin_sel:
					_hW[obsnum,seq,vst,var]=hW[k]
					if self.SS_TYPE=="cmb_non_vst_SE":
						_hWerr[obsnum,seq,vst,var]=hWerr[k]
					
			#! Determine max and set Y-axis min(=0),max for all h1p[seq,vst,var]
			maxl=[_hW[k].GetMaximum() for k in _hW]
			maximum=max(maxl)
			for k in _hW:
				_hW[k].SetMinimum(0.)
				_hW[k].SetMaximum(maximum)

			#! 2.ii. Get 'results' from these _hW 
			if   self.SS_TYPE=="per_non_vst_SE": results=self.get_per_non_vst_SE_results(_hW)
			elif self.SS_TYPE=="cmb_non_vst_SE": results=self.get_cmb_non_vst_SE_results(_hW,_hWerr)
			#! unpack results(=h1_obs,h1_err,h1_rel_err,rslts)
			hW_obs,hW_err,hW_rel_err,rslts=results[0],results[1],results[2],results[3]

			#! 2.iii. Collect all hW_obs in q2-dict: This will be used later to plot all xsec(W) for all Q2 in one canvas (step 3.)
			hW_obs_Q2[q2wbin_sel]=hW_obs

			#! 2.iv. For every EF*[vst]*THETA plot 'results'
			#! + Note that text results are written after the following loop, which skips "EC"
			#!   since text results also include EC
			#!Now constructed EF*[vst] domain to loop over and make plots
			d1=list(itertools.product(self.SEQS,VSTS))
			for r1 in d1:
				seq,vst=r1[0],r1[1]
				if seq=="EC": continue
				#! prepare outdirs
				outdir_hW_obs_err =os.path.join(self.OUTDIR_ITG,"VST%d"%(vst))
				outdir_hW_err     =os.path.join(self.OUTDIR_ITG_ERR,"VST%d"%(vst))
				outdir_hW_err_vsl =os.path.join(self.OUTDIR_ITG_ERR_VSL,"VST%d"%(vst))
				for d in [outdir_hW_obs_err,outdir_hW_err,outdir_hW_err_vsl]:
					if not os.path.exists(d):
						os.makedirs(d)
				
				#! di: Final result(mu,sg_mu,sg,sg_sg_sg) from SE variations: hW_obs[seq,vst,var] and hW_err[seq,vst,var] together
				ROOT.gStyle.SetOptStat(0)
				c=ROOT.TCanvas("c","c",1200,800)
				pad_p=ROOT.TPad("pad_p","Plots pad",0.00,0.00,1.00,0.90)
				pad_t=ROOT.TPad("pad_p","Plots pad",0.00,0.90,1.00,1.00)
				pad_p.Draw()
				pad_t.Draw()

				pad_t.cd()
				pt=ROOT.TPaveText(.05,.1,.95,.8)
				q2bin=self.get_q2bin(q2wbin_sel)
				pt.AddText("Integrated cross-section for Q2=%s,VST%d"%(q2bin,vst))
  				pt.SetTextSize(0.32)
  				pt.Draw()

				pad_p.cd()
				#! Set up some aesthetics particular to hW
				#! + Give already rotated bin labels by 90 deg. more space
				#! + Move xaxis title lower to accomodate 90 deg. rotated titles
				pad_p.SetBottomMargin(0.20)
				hW_obs[seq,vst,'THETA'].GetXaxis().SetTitleOffset(2.80)
				hW_err[seq,vst,'THETA'].GetXaxis().SetTitleOffset(2.80)
				#! Setup aesthetics for hW_err
				#! Since transparency is not working (for details see hist_R2_athtcs_non_vst_SE) , the following is a good replacement
				hW_err[seq,vst,'THETA'].SetFillColor(ROOT.gROOT.ProcessLine("kBlack"))
				hW_err[seq,vst,'THETA'].SetFillStyle(3003)
				hW_err[seq,vst,'THETA'].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
				hW_err[seq,vst,'THETA'].SetMarkerStyle(ROOT.gROOT.ProcessLine("kDot"))
				hW_err[seq,vst,'THETA'].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
				#! Draw
				hW_obs[seq,vst,'THETA'].SetTitle("")
				hW_err[seq,vst,'THETA'].SetTitle("")
				#! Make sure hist limits are Ok
				#print "old for hW_obs", q2wbin_sel," =",hW_obs[seq,vst,'THETA'].GetMaximum()
				self.set_h1_maximum_minimum(hW_obs[seq,vst,'THETA'],10/100,minimum=0)
				#print "new for hW_obs", q2wbin_sel," =",hW_obs[seq,vst,'THETA'].GetMaximum()
				self.set_h1_maximum_minimum(hW_err[seq,vst,'THETA'],10/100,minimum=0)
				#! Draw
				hW_obs[seq,vst,'THETA'].Draw()
				hW_err[seq,vst,'THETA'].Draw("hist e same")
				c.SaveAs("%s/c_q%s.png"%(outdir_hW_obs_err,q2wbin_sel))
				c.SaveAs("%s/c_q%s.eps"%(outdir_hW_obs_err,q2wbin_sel))

				#! dii.  To just see relative error (rel_err, sg_rel_err):      hW_rel_err[seq,vst,var]
				ROOT.gStyle.SetOptStat(0)
				c=ROOT.TCanvas("c","c",1200,800)
				pad_p=ROOT.TPad("pad_p","Plots pad",0.00,0.00,1.00,0.90)
				pad_t=ROOT.TPad("pad_p","Plots pad",0.00,0.90,1.00,1.00)
				pad_p.Draw()
				pad_t.Draw()

				pad_t.cd()
				pt=ROOT.TPaveText(.05,.1,.95,.8)
				q2bin=self.get_q2bin(q2wbin_sel)
				pt.AddText("Relative error for integrated cross-section for Q2=%s,VST%d"%(q2bin,vst))
  				pt.SetTextSize(0.32)
  				pt.Draw()

				pad_p.cd()
				#! Set up some aesthetics particular to hW
				#! + Give already rotated bin labels by 90 deg. more space
				#! + Move xaxis title lower to accomodate 90 deg. rotated titles
				pad_p.SetBottomMargin(0.20)
				hW_rel_err[seq,vst,'THETA'].GetXaxis().SetTitleOffset(2.80)
				#! Draw
				hW_rel_err[seq,vst,'THETA'].SetTitle("")
				hW_rel_err[seq,vst,'THETA'].SetYTitle("error [%] in #sigma")
				hW_rel_err[seq,vst,'THETA'].SetMinimum(0)
				#! Make sure maximum is displayed
				#print "old for", q2wbin_sel," =",hW_rel_err[seq,vst,'THETA'].GetMaximum()
				self.set_h1_maximum_minimum(hW_rel_err[seq,vst,'THETA'],10/100,minimum=0)
				#hW_rel_err[seq,vst,'THETA'].SetMaximum(40)
				#print "new for", q2wbin_sel," =",hW_rel_err[seq,vst,'THETA'].GetMaximum()
				#! Draw
				hW_rel_err[seq,vst,'THETA'].Draw()
				c.SaveAs("%s/c_q%s.png"%(outdir_hW_err,q2wbin_sel))
				c.SaveAs("%s/c_q%s.eps"%(outdir_hW_err,q2wbin_sel))

				#! diii.: Plot _hW(obsnum,seq,vst,var) directly
				ROOT.gStyle.SetOptStat(0)
				c=ROOT.TCanvas("c","c",1200,800)
				pad_p=ROOT.TPad("pad_p","Plots pad",0.00,0.00,1.00,0.90)
				pad_t=ROOT.TPad("pad_p","Plots pad",0.00,0.90,1.00,1.00)
				pad_p.Draw()
				pad_t.Draw()

				pad_t.cd()
				pt=ROOT.TPaveText(.05,.1,.95,.8)
				q2bin=self.get_q2bin(q2wbin_sel)
				pt.AddText("Integrated cross-section for Q2=%s,VST%d"%(q2bin,vst))
  				pt.SetTextSize(0.32)
  				pt.Draw()

				pad_p.cd()
				#! Set up some aesthetics particular to hW
				#! + Give already rotated bin labels by 90 deg. more space
				#! + Move xaxis title lower to accomodate 90 deg. rotated titles
				pad_p.SetBottomMargin(0.20)
				#! _hW[obsnum,seq,vst,var].GetXaxis().SetTitleOffset(2.80) is done inside obsnum loop that follows	
								
				for i,obsnum in enumerate(self.OBSD.keys()):#obsnum in obsnuml: 
					#var='THETA'
					obstag=self.OBSD[obsnum][1]
					draw_opt=""
					if i>0: draw_opt="same"
					#! setup histogram aesthetics
					_hW[obsnum,seq,vst,'THETA'].SetTitle("")
					#! + Move xaxis title lower to accomodate 90 deg. rotated titles
					_hW[obsnum,seq,vst,'THETA'].GetXaxis().SetTitleOffset(2.80)
					mrk=ROOT.gROOT.ProcessLine("kFullDotLarge")
					clr=CLRS_N[obsnum-1]
					_hW[obsnum,seq,vst,'THETA'].SetMarkerStyle(mrk)
					_hW[obsnum,seq,vst,'THETA'].SetMarkerColor(clr)
					_hW[obsnum,seq,vst,'THETA'].SetLineColor(clr)
					#! aesthetics for _hWerr
					if self.SS_TYPE=="cmb_non_vst_SE":
						_hWerr[obsnum,seq,vst,'THETA'].SetFillColor(clr)
						_hWerr[obsnum,seq,vst,'THETA'].SetFillStyle(3003+i) #! +i to modify fill style
						_hWerr[obsnum,seq,vst,'THETA'].SetMarkerStyle(ROOT.gROOT.ProcessLine("kDot"))
						_hWerr[obsnum,seq,vst,'THETA'].SetMarkerColor(clr)
					#! Draw
					_hW[obsnum,seq,vst,'THETA'].Draw(draw_opt)
					if self.SS_TYPE=="cmb_non_vst_SE":
						_hWerr[obsnum,seq,vst,'THETA'].Draw("hist e same")
				#! Also draw hW_obs and hW_err here
				hW_obs[seq,vst,'THETA'].Draw("same")
				hW_err[seq,vst,'THETA'].Draw("hist e same")
				c.SaveAs("%s/c_q%s.png"%(outdir_hW_err_vsl,q2wbin_sel))
				c.SaveAs("%s/c_q%s.eps"%(outdir_hW_err_vsl,q2wbin_sel))

			#! div.  Text file of full results: put rslts:
			#! + EF,vst,THETA -> OUTDIR_ITG_TXT_RSLT/vst/
			#! + EC,vst,THETA -> OUTDIR_ITG_TXT_RSLT/EC/vst/
			#! Prepate main outdir
			d1=list(itertools.product(self.SEQS,VSTS))
			for r1 in d1:
				seq,vst=r1[0],r1[1]
				#! make outdir
				if seq=="EF":
					outdir_hW_txt_rslt_q2wbin=os.path.join(self.OUTDIR_ITG_TXT_RSLT,"VST%d"%(vst),q2wbin_sel)
				elif seq=="EC":
					outdir_hW_txt_rslt_q2wbin=os.path.join(self.OUTDIR_ITG_TXT_RSLT,"EC","VST%d"%(vst),q2wbin_sel)
				else:
					print "ERROR: plot_SE_results_W(): seq=%s not recognized while printing txt results"%seq		
				if not os.path.exists(outdir_hW_txt_rslt_q2wbin):
					os.makedirs(outdir_hW_txt_rslt_q2wbin)
				#! now write
				ftxt=open("%s/result.txt"%(outdir_hW_txt_rslt_q2wbin),"w")
					#! Loop over rslts per bin in rslts and write to file
				for d in rslts[seq,vst,'THETA']:
					bin,bin_le,bin_ue,mu,sg_mu,sg,sg_sg,rel_err,sg_rel_err=d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8]
					ftxt.write("%d,%.3f,%.3f = (%f +/- %f),(%f +/- %f),(%f +/- %f)\n"%(bin,bin_le,bin_ue,mu,sg_mu,sg,sg_sg,rel_err,sg_rel_err))
				ftxt.close()
			

			#! 2.v. For every [EF,EC]*[vst]*THETA plot only hW_obs
			for vst in VSTS:
				#! prepare outdirs
				outdir_hW_obs=os.path.join(self.OUTDIR_ITG,"EC-EH","VST%d"%(vst))
				if not os.path.exists(outdir_hW_obs):
						os.makedirs(outdir_hW_obs)

				ROOT.gStyle.SetOptStat(0)
				c=ROOT.TCanvas("c","c",1200,800)
				pad_p=ROOT.TPad("pad_p","Plots pad",0.00,0.00,1.00,0.90)
				pad_t=ROOT.TPad("pad_p","Plots pad",0.00,0.90,1.00,1.00)
				pad_p.Draw()
				pad_t.Draw()
				pad_t.cd()
				pt=ROOT.TPaveText(.05,.1,.95,.8)
				q2bin=self.get_q2bin(q2wbin_sel)
				pt.AddText("Integrated cross-section for Q2=%s,VST%d"%(q2bin,vst))
  				pt.SetTextSize(0.32)
  				pt.Draw()
				#! Display i: Final result(mu,sg_mu,sg,sg_sg_sg) from SE variations: hW_obs[seq,vst,var] and hW_err[seq,vst,var] together
				for i,seq in enumerate(['EF','EC']): #! self.SEQS purposefully modified in this manner so that 'EF' is drawn first
					pad_p.cd()
					#! Set up some aesthetics particular to hW
					#! + Give already rotated bin labels by 90 deg. more space
					#! + Move xaxis title lower to accomodate 90 deg. rotated titles
					pad_p.SetBottomMargin(0.20)
					hW_obs[seq,vst,'THETA'].GetXaxis().SetTitleOffset(2.80)
					#! Remove title
					hW_obs[seq,vst,'THETA'].SetTitle("")
					#! Make sure hist limits are Ok
					self.set_h1_maximum_minimum(hW_obs[seq,vst,'THETA'],10/100,minimum=0)
					#! Setup aesthetics as per EC-EH comparison
					hW_obs[seq,vst,'THETA'].SetMarkerColor(CLRS_PLT_SEQ_ALL[seq])
					hW_obs[seq,vst,'THETA'].SetLineColor(CLRS_PLT_SEQ_ALL[seq])
					#! Draw
					draw_opt=""
					if i>0: draw_opt="sames"
					hW_obs[seq,vst,'THETA'].Draw(draw_opt)

				#! Draw TLegend 
				l=ROOT.TLegend(0.7,0.80,0.9,0.90)
				for seq in self.SEQS:
					l.AddEntry(hW_obs[seq,vst,'THETA'],seq,"p")
				l.Draw()

				c.SaveAs("%s/c_q%s.png"%(outdir_hW_obs,q2wbin_sel))
				c.SaveAs("%s/c_q%s.eps"%(outdir_hW_obs,q2wbin_sel))

			#! Results stored in FOUT: [EC,EH]* [hW_obs,hW_err,hW_rel_err]
			#! output structure: q2wbin/hW_obs(err/rel_err)_seq_vst_var
			fout_q2wdir=self.FOUT_OBS_ITG.GetDirectory(q2wbin_sel)
			if fout_q2wdir==None:
				fout_q2wdir=self.FOUT_OBS_ITG.mkdir(q2wbin_sel)
			fout_q2wdir.cd()

			#! Now save [EC,EH]*hW_obs,hW_err,hW_rel_err
			d1=list(itertools.product(self.SEQS,VSTS))
			for r1 in d1:
				seq,vst=r1[0],r1[1]
				#! Save hists
				#! + Note the hists have to be named, again because:
				#! 		i. somewhere earlier their titles are reset to "" 
				#!		AND, more importantly,
				#!		ii. 'obs','err', and 'rel_err' need to be added to their names
				hW_obs[seq,vst,'THETA'].SetName("hW_obs_%s_%d_%s"%(seq,vst,'THETA'))
				hW_obs[seq,vst,'THETA'].Write()
				hW_err[seq,vst,'THETA'].SetName("hW_err_%s_%d_%s"%(seq,vst,'THETA'))
				hW_err[seq,vst,'THETA'].Write()
				if seq=="EC": #! Set maximum for hW_rel_err (for EF is already set in canvas displays) 
					self.set_h1_maximum_minimum(hW_rel_err[seq,vst,'THETA'],10/100,minimum=0)
				hW_rel_err[seq,vst,'THETA'].SetName("hW_rel_err_%s_%d_%s"%(seq,vst,'THETA'))
				hW_rel_err[seq,vst,'THETA'].Write()
					
		#! 3. Plot, per vst, xsec(W;EF,THETA) for all Q2 in one canvas.
		for vst in VSTS:
			#! create canvas for vst
			c=ROOT.TCanvas("c","c",1200,800)
			pad_p=ROOT.TPad("pad_p","Plots pad",0.00,0.00,1.00,0.90)
			pad_t=ROOT.TPad("pad_p","Plots pad",0.00,0.90,1.00,1.00)
			pad_p.Draw()
			pad_t.Draw()
			pad_t.cd()
			pt=ROOT.TPaveText(.05,.1,.95,.8)
			pt.AddText("Integrated cross-section for VST%d"%(vst))
  			pt.SetTextSize(0.32)
  			pt.Draw()
  			#! create legend for this canvas
  			l=ROOT.TLegend(0.7,0.8,0.9,1.0)
  			#! Setup maximum of hists that depends on q2
  			#maxl=[hW_obs['EF',vst,'THETA'].GetMaximum() for q2wbin in hW_obs_Q2.keys()]
  			maxl=[]
  			for q2wbin in hW_obs_Q2.keys():
				hW_obs=hW_obs_Q2[q2wbin]
				maxl.append(hW_obs['EF',vst,'THETA'].GetMaximum())
  			mxm=max(maxl)
  			mxm=mxm+((10/100)*mxm)
			#! Loop over q2wbin keys in hW_obs_Q2.keys()
			for i,q2wbin in enumerate(hW_obs_Q2.keys()):
				#! get q2wbin information to be used later for, mainly, aesthetic reasons
				q2bin=self.get_q2bin(q2wbin)
				wbin=self.get_wbin(q2wbin)
				#! get histogram for this q2wbin
				hW_obs=hW_obs_Q2[q2wbin]
				#! Setup aesthetics
				hW_obs['EF',vst,'THETA'].SetLineColor(CLRS_N[i+1])
				hW_obs['EF',vst,'THETA'].SetMarkerColor(CLRS_N[i+1])
				#! Set up maximum
				hW_obs['EF',vst,'THETA'].SetMaximum(mxm)
				#! Draw
				pad_p.cd()
				draw_opt=""
				if i>0: draw_opt="sames"				
				hW_obs['EF',vst,'THETA'].Draw(draw_opt)
				#! add to legend
				l.AddEntry(hW_obs['EF',vst,'THETA'],"Q2=%s"%q2bin,"p")
			#! Draw legend
			l.Draw()
			#! Save canvas
			outdir=os.path.join(self.OUTDIR_ITG,"VST%d"%(vst))
			#! outdir should already exists, but I see no harm in still doing following
			if not os.path.exists(outdir):
				os.makedirs(outdir)
			c.SaveAs("%s/c_qall_w%s.png"%(outdir,wbin))
			c.SaveAs("%s/c_qall_w%s.eps"%(outdir,wbin))

		print "*** plot_SE_results_W() Done ***\n"

	def plot_cmb_vst_SE_results_hW(self,hW,hWerr):
		'''
		[09-07-16]
		+ Note that plot_SE_results_hW(self,hW,hWerr) was used as a the basis for this and that will reflect in 
		  the code

		+ This function:
			1. Obtains cmb_vst_SE related 'results': mu,sg_mu,sg,sg_sg,f_vst_i,sg_f_vst_i=f(q2b,wb,seq) from hW and hW_err, which are still
		       a function of vst.
			2. New cmb_vst_SE Obs_itg are created and plotted using these 'results' (hW_obs, hW_err)
			3. These 'results' are returned to be used for cmb_vst_SE Obs_1D

		The following are the ways in which the results are displayed:
			1. EH: display all 'results' for 'EH'
			2. EC: write only text part of 'results'
			3. EC_EH comparison: Display only 'h1_obs'*['EC','EH']
			4. Results stored in FOUT: [EC,EH]* [hW_obs,hW_err,hW_rel_err] 

		Note that the following dictionary keys are fixed
			+ obsnum=1 always, 
			+ var='THETA' always

		Process:
		0. Create 'rslts' dictionary that will store cmb_vst_SE related 'results': mu,sg_mu,sg,sg_sg,f_vst_i,sg_f_vst_i=f(q2b,wb,seq)
		   and be returned
		1. Get q2wbinl from hW(q2wbin,obsnum,seq,vst,var)
		2. Loop over the q2wbins and for each bin:
			2.i.   Collect all hW in _hW(1,seq,vst,var) hWerr in _hWerr(1,seq,vst,var)
			2.ii.  For a particular (seq,var), using all vst, obtain '_results'(seq,var) using get_cmb_vst_SE_results_hW():
				where '_results':
				(NOTE '_' used here for '_results(seq,var)' and '_rslts(seq,var)' to keep it separate from 'rslts(q2b,wb,seq)'
				 where it is understood that var='THETA' always )
				+ hW_obs[seq,var]
			    + hW_err[seq,var]
			    + hW_rel_err[seq,var]
			    + _rslts[seq,var]
			2.ii.i. Put '_results:_rslts[seq,var]' in 'rslts[q2b,wb,seq]' (the object that is returned by this function) 
			2.iii.  Collect all hW_obs in q2-dict: This will be used later to plot all xsec(W) for all Q2 in one canvas (step 3.)	 
			2.iv. For every EF*THETA display:
				   (+ Note that full display i.e. [seq]*[vst]*[var] should be done individual SE results: $OBSDIR_E16/SE/cutsncorsX
				   	+ From there I have already learned that there is only variations due to vst and not in vars within a vst)
				di.   Final result(mu,sg_mu,sg,sg_sg_sg) from vst_SE variations: hW_obs[seq,vst,var] and hW_err[seq,vst,var] together
				dii.  To just see relative error (rel_err, sg_rel_err):          hW_rel_err[seq,vst,var]
				diii. To visually see SE variations:                             _hW(obsnum,seq,vst,var) directly
				For every [EC,EF]*THETA display:
				div.  Text file of full results:                                 put rslts[seq,vst,var] in a text file
			2.v. EC-EH comparison:  For every [EC,EH]*THETA display:
				i.   Final result(mu,sg_mu,sg,sg_sg_sg) from vst_SE variations: hW_obs[seq,'THETA']
		3. Plot, per vst, xsec(W;EF,THETA) for all Q2 in one canvas.
		4. Return 'rslts'
				
		'''
		print "*** In DispObs::plot_cmb_vst_SE_results_hW() ***"
		#! Plotting aesthetics
		ROOT.gStyle.Reset() #! Reset aesthetics which may have been set by some other plot function
		ROOT.gStyle.SetOptStat("neiuo")
		ROOT.gStyle.SetErrorX(0.001)

		#! 0. Create 'rslts' dictionary
		rslts=OrderedDict()

		#! 1. Get list of q2wbins in hW(q2wbin,obsnum,seq,vst,var)
		#! Get all q2wbins
		q2wbinl=[k[0] for k in hW.keys()]
		#! Keep only unique bins
		q2wbinl=list(set(q2wbinl))
		print "DispObs::plot_cmb_vst_SE_results_hW() q2wbinl in hW=",q2wbinl

		#! Create a dictionary to keep all hW_obs labelled by Q2 (see steps 2.iii. and 3.)
		hW_obs_Q2=OrderedDict()
		#! 2. Loop over the q2wbins
		for iq2wbin,q2wbin_sel in enumerate(q2wbinl):
			#! Store q2bin information from q2wbin_sel
			#! wbin here is the full W range: 1,400-2.125 and not needed to be stored in any meaningful way
			q2bin=self.get_q2bin(q2wbin_sel)
			#! 2.i. Collect all hW for a particular q2wbin in _hW(obsnum,seq,vst,var)
			_hW=OrderedDict()
			_hWerr=OrderedDict()
			for k in hW:
				q2wbin,obsnum,seq,vst,var=k[0],k[1],k[2],k[3],k[4]
				if q2wbin==q2wbin_sel:
					_hW[obsnum,seq,vst,var]=hW[k]
					_hWerr[obsnum,seq,vst,var]=hWerr[k]
					
			#! Determine max and set Y-axis min(=0),max for all h1p[seq,vst,var]
			maxl=[_hW[k].GetMaximum() for k in _hW]
			maximum=max(maxl)
			for k in _hW:
				_hW[k].SetMinimum(0.)
				_hW[k].SetMaximum(maximum)

			#! 2.ii. Get 'results_q2wb' from these _hW 
			_results=self.get_cmb_vst_SE_results_hW(_hW,_hWerr)
			#! unpack results(=h1_obs,h1_err,h1_rel_err,rslts)
			hW_obs,hW_err,hW_rel_err,_rslts=_results[0],_results[1],_results[2],_results[3]

			#! 2.ii.i. Put '_results:_rslts[seq,var]' in 'rslts[q2b,wb,seq]' (the object that is returned by this function)
			for seq in self.SEQS:
				for d in _rslts[seq,'THETA']:
					bin,bin_le,bin_ue,mu,sg_mu,sg,sg_sg,rel_err,sg_rel_err=d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8]
					f_vst_1,sg_f_vst_1=d[9],d[10]
					f_vst_2,sg_f_vst_2=d[11],d[12]
					f_vst_3,sg_f_vst_3=d[13],d[14]
					sg_vst_1D,sg_sg_vst_1D=d[15],d[16]

					wbin="%.3f-%.3f"%(bin_le,bin_ue)
					rslts[q2bin,wbin,seq]=[mu,sg_mu,sg,sg_sg,
					                       f_vst_1,sg_f_vst_1,
					                       f_vst_2,sg_f_vst_2,
					                       f_vst_3,sg_f_vst_3,
					                       sg_vst_1D,sg_sg_vst_1D]
			
			#! 2.iii. Collect all hW_obs in q2-dict: This will be used later to plot all xsec(W) for all Q2 in one canvas (step 3.)
			hW_obs_Q2[q2wbin_sel]=hW_obs

			#! 2.iv. For every EF*THETA plot 'results'
			#! + Note that text results are written after the following loop, which skips "EC"
			#!   since text results also include EC
			#!Now constructed EF*[vst] domain to loop over and make plots
			d1=list(itertools.product(self.SEQS,VSTS))
			for r1 in d1:
				seq,var=r1[0],r1[1]
				if seq=="EC": continue
				#! prepare outdirs
				outdir_hW_obs_err =os.path.join(self.OUTDIR_ITG)
				outdir_hW_err     =os.path.join(self.OUTDIR_ITG_ERR)
				outdir_hW_err_vsl =os.path.join(self.OUTDIR_ITG_ERR_VSL)
				for d in [outdir_hW_obs_err,outdir_hW_err,outdir_hW_err_vsl]:
					if not os.path.exists(d):
						os.makedirs(d)
				
				#! di: Final result(mu,sg_mu,sg,sg_sg_sg) from vst_SE variations: hW_obs[seq,var] and hW_err[seq,var] together
				ROOT.gStyle.SetOptStat(0)
				c=ROOT.TCanvas("c","c",1200,800)
				pad_p=ROOT.TPad("pad_p","Plots pad",0.00,0.00,1.00,0.90)
				pad_t=ROOT.TPad("pad_p","Plots pad",0.00,0.90,1.00,1.00)
				pad_p.Draw()
				pad_t.Draw()

				pad_t.cd()
				pt=ROOT.TPaveText(.05,.1,.95,.8)
				q2bin=self.get_q2bin(q2wbin_sel)
				pt.AddText("Integrated cross-section for Q2=%s"%(q2bin))
  				pt.SetTextSize(0.32)
  				pt.Draw()

				pad_p.cd()
				#! Set up some aesthetics particular to hW
				#! + Give already rotated bin labels by 90 deg. more space
				#! + Move xaxis title lower to accomodate 90 deg. rotated titles
				pad_p.SetBottomMargin(0.20)
				hW_obs[seq,'THETA'].GetXaxis().SetTitleOffset(2.80)
				hW_err[seq,'THETA'].GetXaxis().SetTitleOffset(2.80)
				#! Setup aesthetics for hW_err
				#! Since transparency is not working (for details see hist_R2_athtcs_non_vst_SE) , the following is a good replacement
				hW_err[seq,'THETA'].SetFillColor(ROOT.gROOT.ProcessLine("kBlack"))
				hW_err[seq,'THETA'].SetFillStyle(3003)
				hW_err[seq,'THETA'].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
				hW_err[seq,'THETA'].SetMarkerStyle(ROOT.gROOT.ProcessLine("kDot"))
				hW_err[seq,'THETA'].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
				#! Draw
				hW_obs[seq,'THETA'].SetTitle("")
				hW_err[seq,'THETA'].SetTitle("")
				#! Make sure hist limits are Ok
				self.set_h1_maximum_minimum(hW_obs[seq,'THETA'],10/100,minimum=0)
				self.set_h1_maximum_minimum(hW_err[seq,'THETA'],10/100,minimum=0)
				#! Draw
				hW_obs[seq,'THETA'].Draw()
				hW_err[seq,'THETA'].Draw("hist e same")
				c.SaveAs("%s/c_q%s.png"%(outdir_hW_obs_err,q2wbin_sel))
				c.SaveAs("%s/c_q%s.eps"%(outdir_hW_obs_err,q2wbin_sel))

				#! dii.  To just see relative error (rel_err, sg_rel_err):      hW_rel_err[seq,vst,var]
				ROOT.gStyle.SetOptStat(0)
				c=ROOT.TCanvas("c","c",1200,800)
				pad_p=ROOT.TPad("pad_p","Plots pad",0.00,0.00,1.00,0.90)
				pad_t=ROOT.TPad("pad_p","Plots pad",0.00,0.90,1.00,1.00)
				pad_p.Draw()
				pad_t.Draw()

				pad_t.cd()
				pt=ROOT.TPaveText(.05,.1,.95,.8)
				pt.AddText("Relative error for integrated cross-section for Q2=%s"%(q2bin))
  				pt.SetTextSize(0.32)
  				pt.Draw()

				pad_p.cd()
				#! Set up some aesthetics particular to hW
				#! + Give already rotated bin labels by 90 deg. more space
				#! + Move xaxis title lower to accomodate 90 deg. rotated titles
				pad_p.SetBottomMargin(0.20)
				hW_rel_err[seq,'THETA'].GetXaxis().SetTitleOffset(2.80)
				#! Draw
				hW_rel_err[seq,'THETA'].SetTitle("")
				hW_rel_err[seq,'THETA'].SetYTitle("error [%] in #sigma")
				hW_rel_err[seq,'THETA'].SetMinimum(0)
				#! Make sure maximum is displayed
				self.set_h1_maximum_minimum(hW_rel_err[seq,'THETA'],10/100,minimum=0)
				#! Draw
				hW_rel_err[seq,'THETA'].Draw()
				c.SaveAs("%s/c_q%s.png"%(outdir_hW_err,q2wbin_sel))
				c.SaveAs("%s/c_q%s.eps"%(outdir_hW_err,q2wbin_sel))

				#! diii.: Plot _hW(obsnum,seq,vst,var) directly
				ROOT.gStyle.SetOptStat(0)
				c=ROOT.TCanvas("c","c",1200,800)
				pad_p=ROOT.TPad("pad_p","Plots pad",0.00,0.00,1.00,0.90)
				pad_t=ROOT.TPad("pad_p","Plots pad",0.00,0.90,1.00,1.00)
				pad_p.Draw()
				pad_t.Draw()

				pad_t.cd()
				pt=ROOT.TPaveText(.05,.1,.95,.8)
				q2bin=self.get_q2bin(q2wbin_sel)
				pt.AddText("Integrated cross-section for Q2=%s"%(q2bin))
  				pt.SetTextSize(0.32)
  				pt.Draw()

				pad_p.cd()
				#! Set up some aesthetics particular to hW
				#! + Give already rotated bin labels by 90 deg. more space
				#! + Move xaxis title lower to accomodate 90 deg. rotated titles
				pad_p.SetBottomMargin(0.20)
				#! _hW[obsnum,seq,vst,var].GetXaxis().SetTitleOffset(2.80) is done inside obsnum loop that follows	
				
				#! Create legend 
				l=ROOT.TLegend(0.7,0.80,0.9,0.90)
				for i,vst in enumerate(VSTS):#obsnum in obsnuml: 
					var='THETA'
					draw_opt=""
					if i>0: draw_opt="same"
					#! setup histogram aesthetics
					_hW[obsnum,seq,vst,var].SetTitle("")
					#! + Move xaxis title lower to accomodate 90 deg. rotated titles
					_hW[obsnum,seq,vst,var].GetXaxis().SetTitleOffset(2.80)
					mrk=ROOT.gROOT.ProcessLine("kFullDotLarge")
					clr=CLRS_N[vst-1]
					_hW[obsnum,seq,vst,var].SetMarkerStyle(mrk)
					_hW[obsnum,seq,vst,var].SetMarkerColor(clr)
					_hW[obsnum,seq,vst,var].SetLineColor(clr)
					#! aesthetics for _hWerr
					_hWerr[obsnum,seq,vst,var].SetFillColor(clr)
					_hWerr[obsnum,seq,vst,var].SetFillStyle(3003+i) #! +i to modify fill style
					_hWerr[obsnum,seq,vst,var].SetMarkerStyle(ROOT.gROOT.ProcessLine("kDot"))
					_hWerr[obsnum,seq,vst,var].SetMarkerColor(clr)
					#! Draw
					_hW[obsnum,seq,vst,var].Draw(draw_opt)
					_hWerr[obsnum,seq,vst,var].Draw("hist e same")
					#! Add to legend
					l.AddEntry(_hW[obsnum,seq,vst,var],"VST%d"%vst,"p")
				#! Also draw hW_obs and hW_err here
				hW_obs[seq,'THETA'].Draw("same")
				hW_err[seq,'THETA'].Draw("hist e same")
				#! Add to legend
				l.AddEntry(hW_obs[seq,'THETA'],"<VST>","p")
				l.Draw()
				#! Save
				c.SaveAs("%s/c_q%s.png"%(outdir_hW_err_vsl,q2wbin_sel))
				c.SaveAs("%s/c_q%s.eps"%(outdir_hW_err_vsl,q2wbin_sel))

			#! div.  Text file of full results: put rslts:
			#! + EF,THETA -> OUTDIR_ITG_TXT_RSLT/
			#! + EC,THETA -> OUTDIR_ITG_TXT_RSLT/EC/
			#! Prepate main outdir
			for seq in self.SEQS:
				#! make outdir
				if seq=="EF":
					outdir_hW_txt_rslt_q2wbin=os.path.join(self.OUTDIR_ITG_TXT_RSLT,q2wbin_sel)
				elif seq=="EC":
					outdir_hW_txt_rslt_q2wbin=os.path.join(self.OUTDIR_ITG_TXT_RSLT,"EC",q2wbin_sel)
				else:
					print "ERROR: plot_cmb_vst_SE_results_hW(): seq=%s not recognized while printing txt results"%seq		
				if not os.path.exists(outdir_hW_txt_rslt_q2wbin):
					os.makedirs(outdir_hW_txt_rslt_q2wbin)
				#! now write
				ftxt=open("%s/result.txt"%(outdir_hW_txt_rslt_q2wbin),"w")
					#! Loop over rslts per bin in rslts and write to file
				for d in _rslts[seq,'THETA']:
					bin,bin_le,bin_ue,mu,sg_mu,sg,sg_sg,rel_err,sg_rel_err=d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8]
					ftxt.write("%d,%.3f,%.3f = (%f +/- %f),(%f +/- %f),(%f +/- %f)\n"%(bin,bin_le,bin_ue,mu,sg_mu,sg,sg_sg,rel_err,sg_rel_err))
				ftxt.close()
			

			#! 2.v. For every [EF,EC]*THETA plot only hW_obs
			#! prepare outdirs
			outdir_hW_obs=os.path.join(self.OUTDIR_ITG,"EC-EH")
			if not os.path.exists(outdir_hW_obs):
					os.makedirs(outdir_hW_obs)

			ROOT.gStyle.SetOptStat(0)
			c=ROOT.TCanvas("c","c",1200,800)
			pad_p=ROOT.TPad("pad_p","Plots pad",0.00,0.00,1.00,0.90)
			pad_t=ROOT.TPad("pad_p","Plots pad",0.00,0.90,1.00,1.00)
			pad_p.Draw()
			pad_t.Draw()
			pad_t.cd()
			pt=ROOT.TPaveText(.05,.1,.95,.8)
			q2bin=self.get_q2bin(q2wbin_sel)
			pt.AddText("Integrated cross-section for Q2=%s"%(q2bin))
  			pt.SetTextSize(0.32)
  			pt.Draw()
			#! Display i: Final result(mu,sg_mu,sg,sg_sg_sg) from SE variations: hW_obs[seq,vst,var] and hW_err[seq,vst,var] together
			for i,seq in enumerate(['EF','EC']): #! self.SEQS purposefully modified in this manner so that 'EF' is drawn first
				pad_p.cd()
				#! Set up some aesthetics particular to hW
				#! + Give already rotated bin labels by 90 deg. more space
				#! + Move xaxis title lower to accomodate 90 deg. rotated titles
				pad_p.SetBottomMargin(0.20)
				hW_obs[seq,'THETA'].GetXaxis().SetTitleOffset(2.80)
				#! Remove title
				hW_obs[seq,'THETA'].SetTitle("")
				#! Make sure hist limits are Ok
				self.set_h1_maximum_minimum(hW_obs[seq,'THETA'],10/100,minimum=0)
				#! Setup aesthetics as per EC-EH comparison
				hW_obs[seq,'THETA'].SetMarkerColor(CLRS_PLT_SEQ_ALL[seq])
				hW_obs[seq,'THETA'].SetLineColor(CLRS_PLT_SEQ_ALL[seq])
				#! Draw
				draw_opt=""
				if i>0: draw_opt="sames"
				hW_obs[seq,'THETA'].Draw(draw_opt)

			#! Draw TLegend 
			l=ROOT.TLegend(0.7,0.80,0.9,0.90)
			for seq in self.SEQS:
				l.AddEntry(hW_obs[seq,'THETA'],seq,"p")
			l.Draw()
			c.SaveAs("%s/c_q%s.png"%(outdir_hW_obs,q2wbin_sel))
			c.SaveAs("%s/c_q%s.eps"%(outdir_hW_obs,q2wbin_sel))

			#! Results stored in FOUT: [EC,EH]* [hW_obs,hW_err,hW_rel_err]
			#! output structure: q2wbin/hW_obs(err/rel_err)_seq_vst_var
			fout_q2wdir=self.FOUT_OBS_ITG.GetDirectory(q2wbin_sel)
			if fout_q2wdir==None:
				fout_q2wdir=self.FOUT_OBS_ITG.mkdir(q2wbin_sel)
			fout_q2wdir.cd()

			#! Now save [EC,EH]*hW_obs,hW_err,hW_rel_err
			for seq in self.SEQS:
				#! Save hists
				#! + Note the hists have to be named, again because:
				#! 		i. somewhere earlier their titles are reset to "" 
				#!		AND, more importantly,
				#!		ii. 'obs','err', and 'rel_err' need to be added to their names
				hW_obs[seq,'THETA'].SetName("hW_obs_%s_%s"%(seq,'THETA'))
				hW_obs[seq,'THETA'].Write()
				hW_err[seq,'THETA'].SetName("hW_err_%s_%s"%(seq,'THETA'))
				hW_err[seq,'THETA'].Write()
				if seq=="EC": #! Set maximum for hW_rel_err (for EF is already set in canvas displays) 
					self.set_h1_maximum_minimum(hW_rel_err[seq,'THETA'],10/100,minimum=0)
				hW_rel_err[seq,'THETA'].SetName("hW_rel_err_%s_%d_%s"%(seq,vst,'THETA'))
				hW_rel_err[seq,'THETA'].Write()
					
		#! 3. Plot, xsec(W;EF,THETA) for all Q2 in one canvas.
		c=ROOT.TCanvas("c","c",1200,800)
		pad_p=ROOT.TPad("pad_p","Plots pad",0.00,0.00,1.00,0.90)
		pad_t=ROOT.TPad("pad_p","Plots pad",0.00,0.90,1.00,1.00)
		pad_p.Draw()
		pad_t.Draw()
		pad_t.cd()
		pt=ROOT.TPaveText(.05,.1,.95,.8)
		pt.AddText("Integrated cross-section")
  		pt.SetTextSize(0.32)
  		pt.Draw()
  		#! create legend for this canvas
  		l=ROOT.TLegend(0.7,0.8,0.9,1.0)
  		#! Setup maximum of hists that depends on q2
  		#! maxl=[hW_obs['EF','THETA'].GetMaximum() for q2wbin in hW_obs_Q2.keys()]
  		maxl=[]
  		for q2wbin in hW_obs_Q2.keys():
			hW_obs=hW_obs_Q2[q2wbin]
			maxl.append(hW_obs['EF','THETA'].GetMaximum())
  		mxm=max(maxl)
  		mxm=mxm+((10/100)*mxm)
  		#! Loop over q2wbin keys in hW_obs_Q2.keys()
		for i,q2wbin in enumerate(hW_obs_Q2.keys()):
			#! get q2wbin information to be used later for, mainly, aesthetic reasons
			q2bin=self.get_q2bin(q2wbin)
			wbin=self.get_wbin(q2wbin)
			#! get histogram for this q2wbin
			hW_obs=hW_obs_Q2[q2wbin]
			#! Setup aesthetics
			hW_obs['EF','THETA'].SetLineColor(CLRS_N[i+1])
			hW_obs['EF','THETA'].SetMarkerColor(CLRS_N[i+1])
			#! Set up maximum
			hW_obs['EF','THETA'].SetMaximum(mxm)
			#! Draw
			pad_p.cd()
			draw_opt=""
			if i>0: draw_opt="sames"				
			hW_obs['EF','THETA'].Draw(draw_opt)
			#! add to legend
			l.AddEntry(hW_obs['EF','THETA'],"Q2=%s"%q2bin,"p")
		#! Draw legend
		l.Draw()
		#! Save canvas
		outdir=os.path.join(self.OUTDIR_ITG)
		#! outdir should already exists, but I see no harm in still doing following
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		c.SaveAs("%s/c_qall_w%s.png"%(outdir,wbin))
		c.SaveAs("%s/c_qall_w%s.eps"%(outdir,wbin))

		#! 4. Return rslts
		return rslts

		print "*** plot_cmb_vst_SE_results_hW() Done ***\n"

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
				# av=round(round((ue+le),3)/2,3)
				# h2[k].GetXaxis().SetBinLabel(iwbin+1,"%.3f"%av)
			for iq2bin in range(nq2bins):
				le=q2bng['BINS'][iq2bin]
				ue=q2bng['BINS'][iq2bin+1]
				h2[k].GetYaxis().SetBinLabel(iq2bin+1,"[%.2f,%.2f)"%(le,ue))
			#h2[k].GetXaxis().LabelsOption("v")
			
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
			q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX,dbg=True,dbg_bins=3)
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
						if self.DO_ETGT_BG_SUB_CORR:
							normf=normf/self.ETGT_BG_SUB_CORR_FCTR
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


	def get_q2wbinlist(self,q2min=0.00,q2max=6.00,wmin=0.000,wmax=3.000,dbg=False,dbg_bins=10,from_obs_itg_yld=False,from_obs_R2=False):
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
			if from_obs_itg_yld==True:
				if   self.SS_TYPE=="per_non_vst_SE": f=root_open(os.path.join(obsdir,'Obs_Itg_Yld_norm/obs_itg_yld.root'))
				elif self.SS_TYPE=="cmb_non_vst_SE" or self.SS_TYPE=="cmb_vst_SE": f=root_open(os.path.join(obsdir,'Obs_itg.root'))
			elif from_obs_R2==True:
				if   self.SS_TYPE=="per_non_vst_SE": 
					if   self.EXTRACT_D_E_FOR_NON_ALPHA==True:  f=root_open(os.path.join(obsdir,'Obs_R2_EC_EF_ST/mthd_phi-proj-fit_w_non-alpha_DE_NQ/obs_R2.root'))
					elif self.EXTRACT_D_E_FOR_NON_ALPHA==False: f=root_open(os.path.join(obsdir,'Obs_R2_EC_EF_ST/mthd_phi-proj-fit_NQ/obs_R2.root'))
				elif self.SS_TYPE=="cmb_non_vst_SE": 
					if   self.EXTRACT_D_E_FOR_NON_ALPHA==True:  f=root_open(os.path.join(obsdir,'Obs_R2_w_non-alpha_DE','Obs_R2.root'))
					elif self.EXTRACT_D_E_FOR_NON_ALPHA==False: f=root_open(os.path.join(obsdir,'Obs_R2.root'))
				elif self.SS_TYPE=="cmb_vst_SE":
					print "DispObs::get_q2wbinlist(): SS_TYPE=%s not applicable for Obs_R2"%self.SS_TYPE
			else:
				if   self.SS_TYPE=="per_non_vst_SE":f=root_open(os.path.join(obsdir,'Obs_1D_norm/obs_1D.root'))
				elif self.SS_TYPE=="cmb_non_vst_SE" or self.SS_TYPE=="cmb_vst_SE":f=root_open(os.path.join(obsdir,'Obs_1D.root'))
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
	
	def hist_1D_athtcs(self,h1,err_type=None):
		'''
		[03-02-16]
		After addition of self.DO_SS_STUDY, this method was modified so that it can understand 
		the type of h1 being passed, which could be either of the following depending on the condition (noted
		in parenthesis and programmed into the logic of this function):

			+ h1(seq,vst,var)           (nkeys=3 and self.DO_SS_STUDY=false)
			
			+ h1(obsnum,seq,vst,var)    (nkeys=4 and self.DO_SS_STUDY=true and err_type=None)
			+ h1err(obsnum,seq,vst,var) (nkeys=4 and self.DO_SS_STUDY=true and err_type="err" and (self.SS_TYPE=="cmb_non_vst_SE" or self.SS_TYPE=="cmb_vst_SE")
			
			+ h1_obs(seq,vst,var)       (nkeys=3 and self.DO_SS_STUDY=true and err_type=None)
			+ h1_err(seq,vst,var)       (nkeys=3 and self.DO_SS_STUDY=true and err_type="err")
			+ h1_rel_err(seq,vst,var)   (nkeys=3 and self.DO_SS_STUDY=true and err_type="rel_err")
		
		'''
		#! First determine if h1 has 3 or 4 keys(3 keys=seq,vst,var; 4 keys=obsnum,seq,vst,var)
		#! + This structure of h1 depends on self.DO_SS_STUDY
		nkeys=len(h1.keys()[0])
		#print("DispObs::hist_1D_athtcs(): nkeys=%d"%nkeys)
		for k in h1.keys():	
			if   nkeys==4: obsnum,seq,vst,var=k[0],k[1],k[2],k[3]
			elif nkeys==3: seq,vst,var=k[0],k[1],k[2]
			else: 		   sys.exit("DispObs::hist_1D_athtcs(): nkeys undetermined. Exiting.")

			h1[k].SetTitle("")
			h1[k].SetXTitle( "%s[%s]"%(VAR_NAMES[(vst,var)],VAR_UNIT_NAMES[var]) )
			#! X-axis title aesthetics
			h1[k].GetXaxis().SetTitleOffset(0.7)
			h1[k].GetXaxis().SetLabelSize(0.05)
			h1[k].GetXaxis().SetTitleSize(0.10)
			#! Y-axis title aesthetics
			h1[k].GetYaxis().SetTitleOffset(0.7)
			h1[k].GetYaxis().SetLabelSize(0.05)
			h1[k].GetYaxis().SetTitleSize(0.10)

			if nkeys==4 and self.DO_SS_STUDY==True and err_type==None:#s i.e. self.DO_SS_STUDY=true:
				if var!='THETA':
					h1[k].SetYTitle("#Delta#sigma/#Delta%s [#mub/%s]"%(VAR_NAMES[(vst,var)],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
				elif var=='THETA': #! only difference is to add the 'Cos' of DCosTheta
					h1[k].SetYTitle("#Delta#sigma/#Deltacos(%s) [#mub/%s]"%(VAR_NAMES[(vst,var)],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))

				clr=CLRS_N[obsnum-1]
				mrk=ROOT.gROOT.ProcessLine("kFullDotLarge")
				h1[k].SetMarkerStyle(mrk)
				h1[k].SetMarkerColor(clr)
				h1[k].SetLineColor(clr)
			elif nkeys==4 and self.DO_SS_STUDY==True and err_type=="err" and (self.SS_TYPE=="cmb_non_vst_SE" or self.SS_TYPE=="cmb_vst_SE"):
				if var!='THETA':
					h1[k].SetYTitle("Absolute error in #Delta#sigma/#Delta%s"%(VAR_NAMES[(vst,var)]))
				elif var=='THETA': #! only difference is to add the 'Cos' of DCosTheta
					h1[k].SetYTitle("Absolute error in #Delta#sigma/#Deltacos(%s)"%(VAR_NAMES[(vst,var)]))
				
				#! Since transparency is not working (for details see hist_R2_athtcs_non_vst_SE) , the following is a good replacement
				clr=CLRS_N[obsnum-1]
				h1[k].SetFillColor(clr)
				h1[k].SetFillStyle(3003+(obsnum-1)) #! +(obsnum-1) to modify fill style
				h1[k].SetLineColor(clr)
				h1[k].SetMarkerStyle(ROOT.gROOT.ProcessLine("kDot"))
				h1[k].SetMarkerColor(clr)
			elif nkeys==3 and self.DO_SS_STUDY==False:
				if self.VIEW=="norm":
					if var!='THETA':
						h1[k].SetYTitle("#Delta#sigma/#Delta%s [#mub/%s]"%(VAR_NAMES[(vst,var)],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
					elif var=='THETA': #! only difference is to add the 'Cos' of DCosTheta
						h1[k].SetYTitle("#Delta#sigma/#Deltacos(%s) [#mub/%s]"%(VAR_NAMES[(vst,var)],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
				elif self.VIEW=="fullana":
					h1[k].SetYTitle("\"Normalized\" counts")
				elif self.VIEW=="ERyield":
					h1[k].SetYTitle("Counts")
				elif self.VIEW=="norm_EC_EF_SF":
					if var!='THETA':
						h1[k].SetYTitle("#Delta#sigma/#Delta%s [#mub/%s]"%(VAR_NAMES[(vst,var)],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
					elif var=='THETA': #! only difference is to add the 'Cos' of DCosTheta
						h1[k].SetYTitle("#Delta#sigma/#Deltacos(%s) [#mub/%s]"%(VAR_NAMES[(vst,var)],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
					#! marker color and style
					h1[k].SetMarkerColor(CLRS_PLT_SEQ_ALL[seq])
					h1[k].SetLineColor(CLRS_PLT_SEQ_ALL[seq])
					h1[k].SetMarkerStyle(MRKS_PLT_SEQ_ALL[seq])
				

				if self.VIEW!="ERyield":
					h1[k].SetMarkerColor(CLRS_PLT_SEQ_ALL[seq])
					h1[k].SetLineColor(CLRS_PLT_SEQ_ALL[seq])
					h1[k].SetMarkerStyle(MRKS_PLT_SEQ_ALL[seq])
			elif nkeys==3 and self.DO_SS_STUDY==True and err_type=="err": #! i.e. self.DO_SS_STUDY=true and it is h1_SS_err
				if var!='THETA':
					h1[k].SetYTitle("Absolute error in #Delta#sigma/#Delta%s"%(VAR_NAMES[(vst,var)]))
				elif var=='THETA': #! only difference is to add the 'Cos' of DCosTheta
					h1[k].SetYTitle("Absolute error in #Delta#sigma/#Deltacos(%s)"%(VAR_NAMES[(vst,var)]))
				
				#! Since transparency is not working (for details see hist_R2_athtcs_non_vst_SE) , the following is a good replacement
				h1[k].SetFillColor(ROOT.gROOT.ProcessLine("kBlack"))
				h1[k].SetFillStyle(3003)
				h1[k].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
				h1[k].SetMarkerStyle(ROOT.gROOT.ProcessLine("kDot"))
				h1[k].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))

			elif nkeys==3 and self.DO_SS_STUDY==True and err_type=="rel_err": #! i.e. self.DO_SS_STUDY=true and it is h1_SS_err
				if var!='THETA':
					#! Remove "Relative from title to save space"
					h1[k].SetYTitle("Error [%%] in #Delta#sigma/#Delta%s"%(VAR_NAMES[(vst,var)]))
				elif var=='THETA': #! only difference is to add the 'Cos' of DCosTheta
					#! Remove "Relative from title to save space"
					h1[k].SetYTitle("Error [%%] in #Delta#sigma/#Deltacos(%s)"%(VAR_NAMES[(vst,var)]))
				
				h1[k].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
				h1[k].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
				h1[k].SetMarkerStyle(ROOT.gROOT.ProcessLine("kFullDotLarge"))
			elif nkeys==3 and self.DO_SS_STUDY==True and err_type==None:
				if var!='THETA':
					h1[k].SetYTitle("#Delta#sigma/#Delta%s [#mub/%s]"%(VAR_NAMES[(vst,var)],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
				elif var=='THETA': #! only difference is to add the 'Cos' of DCosTheta
					h1[k].SetYTitle("#Delta#sigma/#Deltacos(%s) [#mub/%s]"%(VAR_NAMES[(vst,var)],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
				h1[k].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
				h1[k].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
				h1[k].SetMarkerStyle(ROOT.gROOT.ProcessLine("kFullDotLarge"))
			
	
	def hist_R2_athtcs(self,hR2d,R2):
		print("DispObs::hist_R2_athtcs()")
		for k in hR2d:
			seq,vst,var=k[0],k[1],k[2]
			if seq=="ST" and (self.VIEW=="EC_EF_ST" or self.VIEW=="EC_ST"): #! modify ST's defined marker color and style to contrast better with EC,EH
				hR2d[k].SetMarkerStyle(MRKS_PLT_SEQ_ALL['SF']) 
				hR2d[k].SetMarkerColor(CLRS_PLT_SEQ_ALL['SF'])
				hR2d[k].SetLineColor(CLRS_PLT_SEQ_ALL['SF'])
			else:
				hR2d[k].SetMarkerStyle(MRKS_PLT_SEQ_ALL[seq]) 
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

	def hist_R2_athtcs_non_vst_SE(self,hR2,err_type=None):
		'''
		[04-14-16]
		+ The only difference between 'hist_R2_athtcs()' and 'hist_R2_athtcs_non_vst_SE()' is in
		  the arguments of the functions and this is related to the differences of the process
		  of 'disp_R2' and 'disp_SE_results_R2': while the former process a single R2 at a time, the latter 
		  not only processes all R2s, but also calculates SS err on them

		+ The aesthetic setup below is as per the hR2 which are reconized by:
			+ hR2(obsnum,R2,seq,vst,var)     (nkeys=5 and err_type=None)
			+ hR2err(obsnum,R2,seq,vst,var)  (nkeys=5 and err_type="err")
			+ hR2_obs(R2,seq,vst,var)        (nkeys=4 and err_type=None)
			+ hR2_err(R2,seq,vst,var)        (nkeys=4 and err_type="err")
			+ hR2_rel_err(seq,vst,var)       (nkeys=4 and err_type="rel_err")
			
		'''
		print("DispObs::hist_R2_athtcs_non_vst_SE()")

		#! First determine number of keys in hR2 on which the subsequent algorithm depends
		nkeys=len(hR2.keys()[0])

		if nkeys==5 and err_type==None:
			for k in hR2:
				obsnum,R2,seq,vst,var=k[0],k[1],k[2],k[3],k[4]

				mrk=ROOT.gROOT.ProcessLine("kFullDotLarge")
				clr=CLRS_N[obsnum-1]
				hR2[k].SetMarkerStyle(mrk) 
				hR2[k].SetMarkerColor(clr)
				hR2[k].SetLineColor(clr)

				hR2[k].SetTitle("")
				hR2[k].SetXTitle( "%s%s"%(VAR_NAMES[(vst,var)],VAR_UNIT_NAMES[var]) )
				#! post VM-fdbk
				#hR2[k].SetYTitle("%s [#mub/%s]"%(R2_NAMED[R2],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
				#! tchncl problem with following
				#hR2[k].SetYTitle("%s^{%s}_{%s} [#mub/%s]"%(R2_NAMED[R2],VAR_NAMES[(vst,var)],VAR_NAMES[(vst,'PHI')],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
				#! Solution 2: _R2 ^{X_{ij}}_{phi_{i}}
				hR2[k].SetYTitle("%s ^{%s}_{%s} [#mub/%s]"%(R2_NAMED[R2],VAR_NAMES[(vst,var)],VAR_NAMES[(vst,'PHI')],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
				#! Solution 2: _{phi_{i}}R2^{X_{ij}}
				#hR2[seq,vst,var].SetYTitle("_{%s}%s^{%s}_{%s} [#mub/%s]"%(VAR_NAMES[(vst,'PHI')],R2_NAMED[R2],VAR_NAMES[(vst,var)],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))

				#! X-axis title aesthetics
				hR2[k].GetXaxis().SetTitleOffset(0.7)
				hR2[k].GetXaxis().SetLabelSize(0.05)
				hR2[k].GetXaxis().SetTitleSize(0.10)
				#! Y-axis title aesthetics
				hR2[k].GetYaxis().SetTitleOffset(0.7)
				hR2[k].GetYaxis().SetLabelSize(0.05)
				hR2[k].GetYaxis().SetTitleSize(0.10)
		elif nkeys==5 and err_type=="err":
			for k in hR2:
				obsnum,R2,seq,vst,var=k[0],k[1],k[2],k[3],k[4]

				clr=CLRS_N[obsnum-1]

				#! Since transparency is not working, the following is a good replacement
				hR2[k].SetFillColor(clr)
				hR2[k].SetFillStyle(3003+(obsnum-1)) #! +(obsnum-1) to modify fill style
				hR2[k].SetLineColor(clr)
				hR2[k].SetMarkerStyle(ROOT.gROOT.ProcessLine("kDot"))
				hR2[k].SetMarkerColor(clr)

				hR2[k].SetTitle("")
				hR2[k].SetXTitle( "%s%s"%(VAR_NAMES[(vst,var)],VAR_UNIT_NAMES[var]) )
				#! post VM-fdbk
				#hR2[k].SetYTitle("%s [#mub/%s]"%(R2_NAMED[R2],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
				#! tchncl problem with following
				#hR2[k].SetYTitle("%s^{%s}_{%s} [#mub/%s]"%(R2_NAMED[R2],VAR_NAMES[(vst,var)],VAR_NAMES[(vst,'PHI')],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
				#! Solution 2: _R2 ^{X_{ij}}_{phi_{i}}
				hR2[k].SetYTitle(" Absolute error in %s ^{%s}_{%s}"%(R2_NAMED[R2],VAR_NAMES[(vst,var)],VAR_NAMES[(vst,'PHI')]))
				#! Solution 2: _{phi_{i}}R2^{X_{ij}}
				#hR2[seq,vst,var].SetYTitle("_{%s}%s^{%s}_{%s} [#mub/%s]"%(VAR_NAMES[(vst,'PHI')],R2_NAMED[R2],VAR_NAMES[(vst,var)],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))

				#! X-axis title aesthetics
				hR2[k].GetXaxis().SetTitleOffset(0.7)
				hR2[k].GetXaxis().SetLabelSize(0.05)
				hR2[k].GetXaxis().SetTitleSize(0.10)
				#! Y-axis title aesthetics
				hR2[k].GetYaxis().SetTitleOffset(0.7)
				hR2[k].GetYaxis().SetLabelSize(0.05)
				hR2[k].GetYaxis().SetTitleSize(0.10)
		elif nkeys==4 and err_type==None:
			for k in hR2:
				R2,seq,vst,var=k[0],k[1],k[2],k[3]
			
				hR2[k].SetMarkerStyle(ROOT.gROOT.ProcessLine("kFullDotLarge")) 
				hR2[k].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
				hR2[k].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))

				hR2[k].SetTitle("")
				hR2[k].SetXTitle( "%s%s"%(VAR_NAMES[(vst,var)],VAR_UNIT_NAMES[var]) )
				#! post VM-fdbk
				#hR2[k].SetYTitle("%s [#mub/%s]"%(R2_NAMED[R2],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
				#! tchncl problem with following
				#hR2[k].SetYTitle("%s^{%s}_{%s} [#mub/%s]"%(R2_NAMED[R2],VAR_NAMES[(vst,var)],VAR_NAMES[(vst,'PHI')],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
				#! Solution 2: _R2 ^{X_{ij}}_{phi_{i}}
				hR2[k].SetYTitle("%s ^{%s}_{%s} [#mub/%s]"%(R2_NAMED[R2],VAR_NAMES[(vst,var)],VAR_NAMES[(vst,'PHI')],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
				#! Solution 2: _{phi_{i}}R2^{X_{ij}}
				#hR2[seq,vst,var].SetYTitle("_{%s}%s^{%s}_{%s} [#mub/%s]"%(VAR_NAMES[(vst,'PHI')],R2_NAMED[R2],VAR_NAMES[(vst,var)],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))

				#! X-axis title aesthetics
				hR2[k].GetXaxis().SetTitleOffset(0.7)
				hR2[k].GetXaxis().SetLabelSize(0.05)
				hR2[k].GetXaxis().SetTitleSize(0.10)
				#! Y-axis title aesthetics
				hR2[k].GetYaxis().SetTitleOffset(0.7)
				hR2[k].GetYaxis().SetLabelSize(0.05)
				hR2[k].GetYaxis().SetTitleSize(0.10)
		elif nkeys==4 and err_type=="err":
			for k in hR2:
				R2,seq,vst,var=k[0],k[1],k[2],k[3]
			
				#! Transparency for SetFillColor, which worked using the same code
				#! in iPython, does not work here and also in a simpler test.py
				# colmg=ROOT.gROOT.GetColor(ROOT.gROOT.ProcessLine("kBlue"))
				# colmg.SetAlpha(0.05)
				# print "colmg=",colmg.GetNumber()
				# hR2[k].SetFillColor(colmg.GetNumber())
				
				#! Since transparency is not working, the following is a good replacement
				hR2[k].SetFillColor(ROOT.gROOT.ProcessLine("kBlack"))
				hR2[k].SetFillStyle(3003)
				hR2[k].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
				hR2[k].SetMarkerStyle(ROOT.gROOT.ProcessLine("kDot"))
				hR2[k].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))

				hR2[k].SetTitle("")
				hR2[k].SetXTitle( "%s%s"%(VAR_NAMES[(vst,var)],VAR_UNIT_NAMES[var]) )
				#! post VM-fdbk
				#hR2[k].SetYTitle("%s [#mub/%s]"%(R2_NAMED[R2],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
				#! tchncl problem with following
				#hR2[k].SetYTitle("%s^{%s}_{%s} [#mub/%s]"%(R2_NAMED[R2],VAR_NAMES[(vst,var)],VAR_NAMES[(vst,'PHI')],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
				#! Solution 2: _R2 ^{X_{ij}}_{phi_{i}}
				hR2[k].SetYTitle(" Absolute error in %s ^{%s}_{%s}"%(R2_NAMED[R2],VAR_NAMES[(vst,var)],VAR_NAMES[(vst,'PHI')]))
				#! Solution 2: _{phi_{i}}R2^{X_{ij}}
				#hR2[seq,vst,var].SetYTitle("_{%s}%s^{%s}_{%s} [#mub/%s]"%(VAR_NAMES[(vst,'PHI')],R2_NAMED[R2],VAR_NAMES[(vst,var)],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))

				#! X-axis title aesthetics
				hR2[k].GetXaxis().SetTitleOffset(0.7)
				hR2[k].GetXaxis().SetLabelSize(0.05)
				hR2[k].GetXaxis().SetTitleSize(0.10)
				#! Y-axis title aesthetics
				hR2[k].GetYaxis().SetTitleOffset(0.7)
				hR2[k].GetYaxis().SetLabelSize(0.05)
				hR2[k].GetYaxis().SetTitleSize(0.10)
		elif nkeys==4 and err_type=="rel_err":
			for k in hR2:
				R2,seq,vst,var=k[0],k[1],k[2],k[3]
			
				hR2[k].SetMarkerStyle(ROOT.gROOT.ProcessLine("kFullDotLarge")) 
				hR2[k].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
				hR2[k].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))

				hR2[k].SetTitle("")
				hR2[k].SetXTitle( "%s%s"%(VAR_NAMES[(vst,var)],VAR_UNIT_NAMES[var]) )
				#! post VM-fdbk
				#hR2[k].SetYTitle("%s [#mub/%s]"%(R2_NAMED[R2],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
				#! tchncl problem with following
				#hR2[k].SetYTitle("%s^{%s}_{%s} [#mub/%s]"%(R2_NAMED[R2],VAR_NAMES[(vst,var)],VAR_NAMES[(vst,'PHI')],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))
				#! Solution 2: _R2 ^{X_{ij}}_{phi_{i}}
				hR2[k].SetYTitle("Error [%%] in %s ^{%s}_{%s}"%(R2_NAMED[R2],VAR_NAMES[(vst,var)],VAR_NAMES[(vst,'PHI')]))
				#! Solution 2: _{phi_{i}}R2^{X_{ij}}
				#hR2[seq,vst,var].SetYTitle("_{%s}%s^{%s}_{%s} [#mub/%s]"%(VAR_NAMES[(vst,'PHI')],R2_NAMED[R2],VAR_NAMES[(vst,var)],VAR_UNIT_NAMES_AFTR_NORM_FCTR_CALC[var]))

				#! X-axis title aesthetics
				hR2[k].GetXaxis().SetTitleOffset(0.7)
				hR2[k].GetXaxis().SetLabelSize(0.05)
				hR2[k].GetXaxis().SetTitleSize(0.10)
				#! Y-axis title aesthetics
				hR2[k].GetYaxis().SetTitleOffset(0.7)
				hR2[k].GetYaxis().SetLabelSize(0.05)
				hR2[k].GetYaxis().SetTitleSize(0.10) #0.10->0.08 to accomodate "Error (...)"
		
		print("DispObs::hist_R2_athtcs_non_vst_SE()")


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
		#! Create output .root file
		self.FOUT_R2=ROOT.TFile(os.path.join(self.OUTDIR_OBS_R2,"obs_R2.root"),"RECREATE")

		#! 1. Get Q2-W binning from file
		#! + Q2-W binning in the file is as per Q2-W binning of h8s
		if self.DBG==True: 	q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX,dbg=True,dbg_bins=3)
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

	def disp_SE_results_R2(self,extract_D_E_for_non_alpha):
		"""
		[08-26-16]
		This function uses the following output .root file from observables in self.OBSD:
			+ if self.SS_TYPE=="per_non_vst_SE"
				+ if self.EXTRACT_D_E_FOR_NON_ALPHA=True:  Obs_R2_EC_EF_ST/mthd_phi-proj-fit_w_non-alpha_DE_NQ/obs_R2.root 
				+ if self.EXTRACT_D_E_FOR_NON_ALPHA=False: Obs_R2_EC_EF_ST/mthd_phi-proj-fit_NQ/obs_R2.root
			+ if self.SS_TYPE=="cmb_non_vst_SE"
				+ if self.EXTRACT_D_E_FOR_NON_ALPHA=True:  Obs_R2_w_non-alpha_DE/Obs_R2.root
				+ if self.EXTRACT_D_E_FOR_NON_ALPHA=False: Obs_R2.root

		-1. Make sure SS_TYPE is valid
				
		0. Make all needed output directories
		
		1. Make displays for Obs_R2
		1.i.   Get Q2-W binning from any of the obs_R2.root files in self.OBSD
		1.ii.  Looper over Q2-W bins and for every Q2-W bin: 
			1.iii  Create:
				+ if self.SS_TYPE=="per_non_vst_SE": 
					+ hR2[obsnum,R2,seq,vst,var]
					+ hR2err=None ('err' refers to Systematic Errors, which at this stage are not calculated)
				where hR2 contains mu,sg_mu which at this point are xsec and its statistical error in each bin.	
				+ if self.SS_TYPE=="cmb_non_vst_SE":
					+ hR2[obsnum,R2,seq,vst,var]    [renamed here from more appropriate named saved in file as hR2_obs(obsnum,seq,vst,var)]
					+ hR2err[obsnum,R2,seq,vst,var] [renamed here from more appropriate named saved in file as hR2_err(obsnum,seq,vst,var)]
				where hR2 and hR2err contain mu,sg_mu (<xsec> and its error) and sg,sg_sg (RMS_<xsec> and it error)
				from "per_non_vst_SE" procedure	
			1.iv.  For [R2,seq,vst,var] using all obsnum, obtain 'results'(R2,seq,vst,var) from:
				+ if self.SS_TYPE=="per_non_vst_SE": 
					+ 'get_per_non_vst_SE_results(hR2)'
			    + if self.SS_TYPE=="cmb_non_vst_SE": 
			    	+ 'get_cmb_non_vst_SE_results(hR2,hR2err)' 
	        where 'results' are (see function for details):
	       		+ hR2_obs[R2,seq,vst,var]
	       		+ hR2_err[R2,seq,vst,var]
	       		+ hR2_rel_err[R2,seq,vst,var]
	       		+ rslts[R2,seq,vst,var]
			1.v.   Display Obs_R2
				i.   Final result(mu,sg_mu,sg,sg_sg_sg) from variations: hR2_obs[R2,seq,vst,var] and hR2_err[seq,vst,var] together
				ii.  To just see relative error (rel_err, sg_rel_err):      hR2_rel_err[R2,seq,vst,var]
				iii. To visually see SE variations:                         hR2(R2,obsnum,seq,vst,var) directly
				iv.  Text file of full results: 							put rslts[R2,seq,vst,var] in a text file
		
		"""
		print "*** In DispObs::disp_SE_results_R2() ***"

		#! -1. Make sure SS_TYPE is valid
		if self.SS_TYPE!="per_non_vst_SE" and self.SS_TYPE!="cmb_non_vst_SE":
			print "DispObs::disp_SE_results_R2(): self.SS_TYPE=%s is not valid. Exiting."%self.SS_TYPE
			sys.exit()

		self.EXTRACT_D_E_FOR_NON_ALPHA=extract_D_E_for_non_alpha
		if self.EXTRACT_D_E_FOR_NON_ALPHA!=True and self.EXTRACT_D_E_FOR_NON_ALPHA!=False:
			print "DispObs::disp_SE_results_R2(): self.EXTRACT_D_E_FOR_NON_ALPHA=%s is not valid. Exiting."%self.EXTRACT_D_E_FOR_NON_ALPHA
			sys.exit()

		#! 0. Make output directories and .root file
		#! Use intermediate dir "Obs_R2_w_non-alpha_DE" only if self.EXTRACT_D_E_FOR_NON_ALPHA=True
		dir_intrmdt=""
		if self.EXTRACT_D_E_FOR_NON_ALPHA==True:  
			dir_intrmdt="Obs_R2_w_non-alpha_DE"
		self.OUTDIR_R2              =os.path.join(self.OUTDIR,dir_intrmdt,"Obs_R2")          #! Display i.
		self.OUTDIR_R2_ERR          =os.path.join(self.OUTDIR,dir_intrmdt,"Obs_R2_err")      #! Display ii.
		self.OUTDIR_R2_ERR_VSL      =os.path.join(self.OUTDIR,dir_intrmdt,"Obs_R2_err_vsl")  #! Display iii.
		self.OUTDIR_R2_TXT_RSLT     =os.path.join(self.OUTDIR,dir_intrmdt,"Obs_R2_txt_rslt") #! Display iv.
		for d in [self.OUTDIR_R2, self.OUTDIR_R2_ERR, self.OUTDIR_R2_ERR_VSL, self.OUTDIR_R2_TXT_RSLT]:
			if not os.path.exists(d):
				os.makedirs(d)
		#! Make output root file
		if self.SS_TYPE!="cmb_vst_SE": #! because if so, this result is simply linked in from "cmb_non_vst_SE"
			self.FOUT_OBS_R2=ROOT.TFile("%s/Obs_R2.root"%(os.path.join(self.OUTDIR,dir_intrmdt)),"RECREATE")
			

		#! 1. Make displays for Obs_R2
		#! 1.i. Get Q2-W binning from file
		if self.DBG==True: q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX,dbg=True,dbg_bins=3,from_obs_R2=True)
		else:              q2wbinl=self.get_q2wbinlist(self.Q2MIN,self.Q2MAX,self.WMIN,self.WMAX,from_obs_R2=True)
		print "DispObs::disp_SE_results_R2():Obs_1D Q2-W bins got from obs_R2.root=",q2wbinl

		#! 1.ii. Looper over Q2-W bins and get hR2(obsnum,R2,seq,vst,var)
		print "DispObs::disp_SE_results_R2():Obs_1D Going to begin processing Q2-W bins from file"
		for q2wbin in q2wbinl:
			print"Processing q2wbin=",q2wbin

			#! 1.iii. Depending on SS_TYPE create hR2(obsnum,R2,seq,vst,var) or hR2_obs(err)(obsnum,R2,seq,vst,var)
			if self.SS_TYPE=="per_non_vst_SE":
				hR2=OrderedDict()
				hR2err=None
			if self.SS_TYPE=="cmb_non_vst_SE":
				hR2=OrderedDict()
				hR2err=OrderedDict()
			for k in self.OBSD:
				obsnum,obs,obstag=k,self.OBSD[k][0],self.OBSD[k][1]
				print "DispObs::disp_SE_results_R2():Obs_R2 Processing obsnum:obs:obstag=",obsnum,obs,obstag

				#! Open file for this obs
				obsdir=os.path.join(self.OBSDIR,obs)
				if self.SS_TYPE=="per_non_vst_SE": 
					if   self.EXTRACT_D_E_FOR_NON_ALPHA==True:  f=root_open(os.path.join(obsdir,'Obs_R2_EC_EF_ST/mthd_phi-proj-fit_w_non-alpha_DE_NQ/obs_R2.root'))
					elif self.EXTRACT_D_E_FOR_NON_ALPHA==False: f=root_open(os.path.join(obsdir,'Obs_R2_EC_EF_ST/mthd_phi-proj-fit_NQ/obs_R2.root'))
				elif self.SS_TYPE=="cmb_non_vst_SE": 
					if   self.EXTRACT_D_E_FOR_NON_ALPHA==True:  f=root_open(os.path.join(obsdir,'Obs_R2_w_non-alpha_DE','Obs_R2.root'))
					elif self.EXTRACT_D_E_FOR_NON_ALPHA==False: f=root_open(os.path.join(obsdir,'Obs_R2.root'))
					
				for item in list(itertools.product(R2S,self.SEQS,VSTS)):
					R2,seq,vst=item[0],item[1],item[2]
					for var in VARS:
						if var=='PHI': continue
						if self.EXTRACT_D_E_FOR_NON_ALPHA==False:
							if (R2=='D' or R2=='E') and var!= 'ALPHA':continue
						if self.SS_TYPE=="per_non_vst_SE":
							hR2[obsnum,R2,seq,vst,var]=f.Get("%s/hR2_%s_%s_%d_%s"%(q2wbin,R2,seq,vst,var))

							#! + The following line was added to decouple histograms from the open
							#! ROOT file which can then be closed.
							#! + This had to be done because otherwise I was getting an error from
							#! rootpy saying that "Too many files are open" (perhaps related to 'ulimit': See notes in ZIM under Technology/UNIX):
							#! rootpy.ROOTError: level=5000, loc='TFile::TFile', msg='file /data/trivedia/e16/2pi/d2pi/highQ2_030916/cutsncors16/sim9/yield.root can not be opened for reading (Too many open files)'
							hR2[obsnum,R2,seq,vst,var].SetDirectory(0)

							#! Call Sumw2() since the errors as currently in the h1s
							#! need to be propagated in all subsequent calculations that involve them
							hR2[obsnum,R2,seq,vst,var].Sumw2()
						elif self.SS_TYPE=="cmb_non_vst_SE":
							hR2[obsnum,R2,seq,vst,var]=f.Get("%s/hR2_obs_%s_%s_%d_%s"%(q2wbin,R2,seq,vst,var))
							hR2err[obsnum,R2,seq,vst,var]=f.Get("%s/hR2_err_%s_%s_%d_%s"%(q2wbin,R2,seq,vst,var))

							#! + The following line was added to decouple histograms from the open
							#! ROOT file which can then be closed.
							#! + This had to be done because otherwise I was getting an error from
							#! rootpy saying that "Too many files are open" (perhaps related to 'ulimit': See notes in ZIM under Technology/UNIX):
							#! rootpy.ROOTError: level=5000, loc='TFile::TFile', msg='file /data/trivedia/e16/2pi/d2pi/highQ2_030916/cutsncors16/sim9/yield.root can not be opened for reading (Too many open files)'
							hR2[obsnum,R2,seq,vst,var].SetDirectory(0)
							hR2err[obsnum,R2,seq,vst,var].SetDirectory(0)

							#! Call Sumw2() since the errors as currently in the h1s
							#! need to be propagated in all subsequent calculations that involve them
							hR2[obsnum,R2,seq,vst,var].Sumw2()
							hR2err[obsnum,R2,seq,vst,var].Sumw2()
				#! End of [seq,vst,var] loop

				#! Close file
				#! + Was getting intermittent error from rootpy saying "Too many files are open": 
				#! rootpy.ROOTError: level=5000, loc='TFile::TFile', msg='file /data/trivedia/e16/2pi/d2pi/highQ2_030916/cutsncors16/sim9/yield.root can not be opened for reading (Too many open files)'
				#f.Close() 
			#! End of OBSD loop

			#! 1.iv. Get results(R2,seq,vst,var) from hR2(obsnum,R2,seq,vst,var)
			if   self.SS_TYPE=="per_non_vst_SE": rslts=self.get_per_non_vst_SE_results(hR2)
			elif self.SS_TYPE=="cmb_non_vst_SE": rslts=self.get_cmb_non_vst_SE_results(hR2,hR2err)
			#print "debug",hR2_SS_err.keys()

			#! 1.v. Plot hR2(obsnum,R2,seq,vst,var)
			self.plot_SE_results_R2(hR2,rslts,q2wbin,hR2err)
		#! End of q2wbin loop

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
						hR2d[seq,vst,var].SetName("hR2_%s_%s_%d_%s"%(R2,seq,vst,var))
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
										
				#! Normalize and apply_rad_eff_corr to hR2
				self.norm_1D(hR2d,q2wbin)
				if self.APPLY_RAD_EFF_CORR==True:
					self.apply_rad_eff_corr(hR2d,q2wbin)
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

		#! For saving output to .root file, create and cd into q2wbindir
		if self.FOUT_R2 !=None:
			q2wbin_dir=self.FOUT_R2.GetDirectory(q2wbin)
			if q2wbin_dir==None: 
				q2wbin_dir=self.FOUT_R2.mkdir(q2wbin)
			q2wbin_dir.cd()
				
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
  		#! [10-26-16] In "new normalization" fix, "old normalization" hists are kept at tmp
  		hR2d_tmp=OrderedDict()
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
				#![03-18-16] Max-amplitude normalization
				#! [10-26-16] + Note that SF hist normalized here is a tmp one. 
				#!            + This is because this normalization is not accurate ("old-normalization")
				#!            + The original SF histogram I will save in in its original form in the .root file
				#!              using which I can implement a "new normalization" when I make plots for the thesis
				print "DispObs::plot_R2() Scaling ampltd-SF to ampltd-EC"
				max_EC=hR2d['EC',vst,var].GetMaximum()
				min_EC=hR2d['EC',vst,var].GetMinimum()
				ampltd_EC=max([max_EC,math.fabs(min_EC)])
				#! First call Sumw2() so that errors are also scaled
				hR2d_tmp['SF',vst,var]=hR2d['SF',vst,var].Clone()
				hR2d_tmp['SF',vst,var].Sumw2()
				#! Now get scale factor and scale bin contents
				max_SF=hR2d_tmp['SF',vst,var].GetMaximum()
				min_SF=hR2d_tmp['SF',vst,var].GetMinimum()
				ampltd_SF=max([max_SF,math.fabs(min_SF)])
				if ampltd_EC==0 or ampltd_SF==0: scale_factor=1
				else:    						 scale_factor=ampltd_EC/ampltd_SF
				hR2d_tmp['SF',vst,var].Scale(scale_factor)

			#! 1. If ST and EC in self.SEQS, then normalize its integral to that of EC
			#! [03-28-16] Normalization changed: max-amplitude of SF normalized to max-amplitude of EC
			if "ST" in self.SEQS and "EC" in self.SEQS:
				#![03-18-16] Max-amplitude normalization
				#! [10-26-16] + Note that ST hist normalized here is a tmp one. 
				#!            + This is because this normalization is not accurate ("old-normalization")
				#!            + The original SF histogram I will save in in its original form in the .root file
				#!              using which I can implement a "new normalization" when I make plots for the thesis
				print "DispObs::plot_R2() Scaling ampltd-ST to ampltd-EC"
				max_EC=hR2d['EC',vst,var].GetMaximum()
				min_EC=hR2d['EC',vst,var].GetMinimum()
				ampltd_EC=max([max_EC,math.fabs(min_EC)])
				#! First call Sumw2() so that errors are also scaled
				hR2d_tmp['ST',vst,var]=hR2d['ST',vst,var].Clone()#"hR2d_ST_%s_%s"%(vst,var))
				hR2d_tmp['ST',vst,var].Sumw2()
				#! Now get scale factor and scale bin contents
				max_ST=hR2d_tmp['ST',vst,var].GetMaximum()
				min_ST=hR2d_tmp['ST',vst,var].GetMinimum()
				ampltd_ST=max([max_ST,math.fabs(min_ST)])
				if ampltd_EC==0 or ampltd_ST==0: scale_factor=1
				else:    						 scale_factor=ampltd_EC/ampltd_ST
				hR2d_tmp['ST',vst,var].Scale(scale_factor)
			
			#! 2. determine and set maximum
			#! #! [10-26-16] code below had to be modified following "new-normalization" change of [10-26-16]
			if "SF" in self.SEQS and "EC" in self.SEQS:
				#! Add hR2d_SF_tmp.GetMinimum() separately
				minl=[hR2d[seq,vst,var].GetMinimum() for seq in self.SEQS if seq != "SF"]
				minl.append(hR2d_tmp['SF',vst,var].GetMinimum())
				#! Add hR2d_SF_tmp.GetMaximum() separately
				maxl=[hR2d[seq,vst,var].GetMaximum() for seq in self.SEQS if seq != "SF"]
				maxl.append(hR2d_tmp['SF',vst,var].GetMaximum())
			elif "ST" in self.SEQS and "EC" in self.SEQS:
				#! Add hR2d_ST_tmp.GetMinimum() separately
				minl=[hR2d[seq,vst,var].GetMinimum() for seq in self.SEQS if seq != "ST"]
				minl.append(hR2d_tmp['ST',vst,var].GetMinimum())
				#! Add hR2d_ST_tmp.GetMaximum() separately
				maxl=[hR2d[seq,vst,var].GetMaximum() for seq in self.SEQS if seq != "ST"]
				maxl.append(hR2d_tmp['ST',vst,var].GetMaximum())
			else:
				minl=[hR2d[seq,vst,var].GetMinimum() for seq in self.SEQS]
				maxl=[hR2d[seq,vst,var].GetMaximum() for seq in self.SEQS]
			#! Set maximum and minimum
			minimum=min(minl)
			maximum=max(maxl)
			fctr_not_A=1.5
			fctr_A=1.1
			if "SF" in self.SEQS and "EC" in self.SEQS:
				for seq in self.SEQS:
					if seq=="SF":continue #! will be done separately		
					if R2!="A": #! B,C,D can oscillate around 0
						hR2d[seq,vst,var].SetMinimum(minimum-math.fabs(fctr_not_A*minimum))
						hR2d[seq,vst,var].SetMaximum(maximum+math.fabs(fctr_not_A*maximum))
					elif R2=="A": #! These distributions are same as Obs-1D
						hR2d[seq,vst,var].SetMinimum(0)
						hR2d[seq,vst,var].SetMaximum(fctr_A*maximum)
				#! For SF
				if R2!="A": #! B,C,D can oscillate around 0
					hR2d_tmp['SF',vst,var].SetMinimum(minimum-math.fabs(fctr_not_A*minimum))
					hR2d_tmp['SF',vst,var].SetMaximum(maximum+math.fabs(fctr_not_A*maximum))
				elif R2=="A": #! These distributions are same as Obs-1D
					hR2d_tmp['SF',vst,var].SetMinimum(0)
					hR2d_tmp['SF',vst,var].SetMaximum(fctr_A*maximum)
			elif "ST" in self.SEQS and "EC" in self.SEQS:
				for seq in self.SEQS:
					if seq=="ST":continue #! will be done separately		
					if R2!="A": #! B,C,D can oscillate around 0
						hR2d[seq,vst,var].SetMinimum(minimum-math.fabs(fctr_not_A*minimum))
						hR2d[seq,vst,var].SetMaximum(maximum+math.fabs(fctr_not_A*maximum))
					elif R2=="A": #! These distributions are same as Obs-1D
						hR2d[seq,vst,var].SetMinimum(0)
						hR2d[seq,vst,var].SetMaximum(fctr_A*maximum)
				#! For ST
				if R2!="A": #! B,C,D can oscillate around 0
					hR2d_tmp['ST',vst,var].SetMinimum(minimum-math.fabs(fctr_not_A*minimum))
					hR2d_tmp['ST',vst,var].SetMaximum(maximum+math.fabs(fctr_not_A*maximum))
				elif R2=="A": #! These distributions are same as Obs-1D
					hR2d_tmp['ST',vst,var].SetMinimum(0)
					hR2d_tmp['ST',vst,var].SetMaximum(fctr_A*maximum)
			else:
				for seq in self.SEQS:		
					if R2!="A": #! B,C,D can oscillate around 0
						hR2d[seq,vst,var].SetMinimum(minimum-math.fabs(fctr_not_A*minimum))
						hR2d[seq,vst,var].SetMaximum(maximum+math.fabs(fctr_not_A*maximum))
					elif R2=="A": #! These distributions are same as Obs-1D
						hR2d[seq,vst,var].SetMinimum(0)
						hR2d[seq,vst,var].SetMaximum(fctr_A*maximum)

			#! 3. Finally draw
			if self.VIEW=="EC_SF_ST":
				hR2d_tmp['SF',vst,var].Draw()
				hR2d['ST',vst,var].Draw("same")
				hR2d['EC',vst,var].Draw("same") #! Draw EC at the end so that its error bar is on top
			elif self.VIEW=="EC_EF":
				hR2d['EC',vst,var].Draw()
				hR2d['EF',vst,var].Draw("same")
			elif self.VIEW=="EC_ST":
				hR2d['ST',vst,var].Draw()
				hR2d['EC',vst,var].Draw("same") #! Draw EC at the end so that its error bar is on top			
			elif self.VIEW=="EC_EF_ST":
				hR2d_tmp['ST',vst,var].Draw()
				hR2d['EF',vst,var].Draw("same")
				hR2d['EC',vst,var].Draw("same")	#! Draw EC at the end so that its error bar is on top
			#! save to .root file
			for seq in self.SEQS:
				hR2d[seq,vst,var].Write()
			
			#! Draw TLine only if ymin<0
			if hR2d['EC',vst,var].GetMinimum()<0:
				ln[pad-1]=ROOT.TLine(hR2d['EC',vst,var].GetXaxis().GetXmin(),0,hR2d['EC',vst,var].GetXaxis().GetXmax(),0)
				ln[pad-1].Draw("same")			

			#! Draw TLegend if pad==3 (since this pad, in the upper right corner, is not shrouded by the TPaveText-Q2W label)
			if pad==3:
				l=ROOT.TLegend(0.7,0.8,0.9,1.0)
				if "SF" in self.SEQS and "EC" in self.SEQS:
					for seq in self.SEQS:
						if seq=="SF": continue #! done separately
						l.AddEntry(hR2d[seq,vst,var],seq,"p")
					l.AddEntry(hR2d_tmp['SF',vst,var],seq,"p")
				if "ST" in self.SEQS and "EC" in self.SEQS:
					for seq in self.SEQS:
						if seq=="ST": continue #! done separately
						l.AddEntry(hR2d[seq,vst,var],seq,"p")
					l.AddEntry(hR2d_tmp['ST',vst,var],seq,"p")
				for seq in self.SEQS:
					l.AddEntry(hR2d[seq,vst,var],seq,"p")
				l.Draw()
			

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

	def plot_SE_results_R2(self,hR2,results,q2wbin,hR2err):
		"""
		[09-04-16]
		+ 'hR2err' added to accomodate plotting of "cmb_non_vst_SE"

		[04-14-16]
		+ 'plot_R2' was used as a tempate for this method
		+ The idea is extend the method to incorporate plotting of a 
		  list of SE-Obs and the associated set of 'results'
		+ Therefore only details relevant to the above were changed and the core functionality left intact

		The following are the ways in which the results are displayed:
			1. EH: display all 'results' for 'EH'
			2. EC: write only text part of 'results'
			3. EC_EH comparison: Display only 'h1_obs'*['EC','EH']
			4. Results stored in FOUT: [EC,EH]* [hR2_obs,hR2_err,hR2_rel_err] 
		"""

		print "In DispObs::plot_SE_results_R2()"

		#! unpack results(=h1_obs,h1_err,h1_rel_err,rslts)
		hR2_obs,hR2_err,hR2_rel_err,rslts=results[0],results[1],results[2],results[3]

		#! Set up some plotting related styles and aesthetics 
		self.plot_R2_athtcs()
		self.hist_R2_athtcs_non_vst_SE(hR2)
		if self.SS_TYPE=="cmb_non_vst_SE": self.hist_R2_athtcs_non_vst_SE(hR2err,"err")
		self.hist_R2_athtcs_non_vst_SE(hR2_obs)
		self.hist_R2_athtcs_non_vst_SE(hR2_err,"err")
		self.hist_R2_athtcs_non_vst_SE(hR2_rel_err,"rel_err")

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

		#! Get q2wbin information that will be used below
		q2bin=self.get_q2bin(q2wbin)
		wbin=self.get_wbin(q2wbin)

		#! 1. EH: display all results for 'EH'
		#! + NOTE, text results are plotted outside of this loop because EC results are written alongwith EF
		#!   and EC is skipped in this loop
		for item in list(itertools.product(R2S,self.SEQS)): #! di,dii, and diii in this loop. div outside since it needs "EC" which this loop skips
			R2,seq=item[0],item[1]
			if seq=="EC": continue
			print "DispObs::plot_SE_results_R2() Plotting hR2 for q2wbin=%s,R2=%s,seq=%s"%(q2wbin,R2,seq)

			#! di. Draw hR2_obs,hR2_err to directly verify mu,sg etc calculated from obsum
			c1=ROOT.TCanvas("c1","c1",1000,1000)
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

				#! 3. Finally draw
				#! Make sure maximum and minimum are displayed
				#! Note that for the following the methods are commented out because default plotting looks OK
				# self.set_h1_maximum_minimum(hR2_obs[R2,seq,vst,var],10/100)
				# self.set_h1_maximum_minimum(hR2_err[R2,seq,vst,var],10/100,minimum=0)
				hR2_obs[R2,seq,vst,var].Draw()
				hR2_err[R2,seq,vst,var].Draw("hist e same")
							
				#! Draw TLine only if ymin<0 for any of hR2 for various obsnum
				if hR2_obs[R2,seq,vst,var].GetMinimum()<0:
					ln[pad-1]=ROOT.TLine(hR2_obs[R2,seq,vst,var].GetXaxis().GetXmin(),0,hR2_obs[R2,seq,vst,var].GetXaxis().GetXmax(),0)
					ln[pad-1].Draw("same")
			
			# outdir_R2=os.path.join(self.OUTDIR_R2,R2)
			# outdir_q2bin=os.path.join(outdir_R2,"q%s"%q2bin)
			outdir_q2bin=os.path.join(self.OUTDIR_R2,"q%s"%q2bin)
			outdir_R2=os.path.join(outdir_q2bin,R2)
			if not os.path.exists(outdir_R2):
				os.makedirs(outdir_R2)
			c1.SaveAs("%s/c_w%s_q%s.png"%(outdir_R2,wbin,q2bin))
			c1.SaveAs("%s/c_w%s_q%s.eps"%(outdir_R2,wbin,q2bin))
			c1.Close()

			#! dii. Draw hR2_rel_err 
			c2=ROOT.TCanvas("c2","c2",1000,1000)
			pad_t=ROOT.TPad("pad_t","Title pad",0.15,0.945,0.85,1.00)
			pad_p=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
  			pad_p.Draw()
  			pad_t.Draw()
  			pad_t.cd()
  			pt=ROOT.TPaveText(.05,.1,.95,.8)
  			pt.AddText("Error in %s for Q2,W = %s"%(R2_NAMED[R2],q2wbin))
  			pt.SetTextSize(0.42)
  			pt.Draw()
  			pad_p.Divide(npadsx,npadsy)
  			
			for item in pad_map:
				pad,vst,var=item[0],item[1],item[2]
				print "pad,vst,var=",pad,vst,var
				gpad=pad_p.cd(pad)

				if self.EXTRACT_D_E_FOR_NON_ALPHA==False:
					if (R2=='D' or R2=='E') and var!= 'ALPHA':continue

				#! 3. Finally draw
				hR2_rel_err[R2,seq,vst,var].SetMinimum(0)
				#! Make sure maximum is drawn
				#! + Note that minimum,maximum=0,200 work better for these
				#!   because some errors are far in excess of 200, at least in test case using all 4 SSBands, (up to >~ 1000) and pcnt method, therefore,
				#!   blows up the scale too much
				#self.set_h1_maximum_minimum(hR2_rel_err[R2,seq,vst,var],10/100,minimum=0, maximum=200)
				#! Keep "blown up" error-scale
				self.set_h1_maximum_minimum(hR2_rel_err[R2,seq,vst,var],10/100,minimum=0)
				#hR2_rel_err[R2,seq,vst,var].SetMaximum(200)
				#! Draw
				hR2_rel_err[R2,seq,vst,var].Draw()
							
			# outdir_R2=os.path.join(self.OUTDIR_R2_ERR,R2)
			# outdir_q2bin=os.path.join(outdir_R2,"q%s"%q2bin)
			# if not os.path.exists(outdir_q2bin):
			# 	os.makedirs(outdir_q2bin)
			# c2.SaveAs("%s/c_w%s_q%s.png"%(outdir_q2bin,wbin,q2bin))
			# c2.SaveAs("%s/c_w%s_q%s.eps"%(outdir_q2bin,wbin,q2bin))
			# c2.Close()
			outdir_q2bin=os.path.join(self.OUTDIR_R2_ERR,"q%s"%q2bin)
			outdir_R2=os.path.join(outdir_q2bin,R2)
			if not os.path.exists(outdir_R2):
				os.makedirs(outdir_R2)
			c2.SaveAs("%s/c_w%s_q%s.png"%(outdir_R2,wbin,q2bin))
			c2.SaveAs("%s/c_w%s_q%s.eps"%(outdir_R2,wbin,q2bin))
			c2.Close()

			#! diii. Draw hR2[obsnum,R2,seq,vst,var] for visual verification
			#! + Draw hR2_obs,hR2_err also to directly verify mu,sg etc calculated from obsum
			#!c=ROOT.TCanvas("c","c",1000,1000)
			c3=ROOT.TCanvas("c3","c3",1000,1000)
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
				#! [04-14-16] I am not sure why the following is not working as expected for plot_SE_results_R2:
				#! the scale gets very large here when I do this. 
				#!	+ This could be because these histograms already have max and min set via plot_R2
				#!    and therefore saved accordinly in the .root file from where I am getting them.
				#! 2. determine and set maximum
				# minl=[hR2[obsnum,R2,seq,vst,var].GetMinimum() for obsnum in self.OBSD.keys()]
				# maxl=[hR2[obsnum,R2,seq,vst,var].GetMaximum() for obsnum in self.OBSD.keys()]
				# minimum=min(minl)
				# maximum=max(maxl)
				# fctr_not_A=1.5
				# fctr_A=1.1
				# for obsnum in self.OBSD.keys():
				# 	if R2!="A": #! B,C,D can oscillate around 0
				# 		hR2[obsnum,R2,seq,vst,var].SetMinimum(minimum-math.fabs(fctr_not_A*minimum))
				# 		hR2[obsnum,R2,seq,vst,var].SetMaximum(maximum+math.fabs(fctr_not_A*maximum))
				# 	elif R2=="A": #! These distributions are same as Obs-1D
				# 		hR2[obsnum,R2,seq,vst,var].SetMinimum(0)
				# 		hR2[obsnum,R2,seq,vst,var].SetMaximum(fctr_A*maximum)
				#! 3. Finally draw
				i=0
				for obsnum in self.OBSD.keys():
					draw_opt=""
					if i>0: draw_opt="sames"
					hR2[obsnum,R2,seq,vst,var].Draw(draw_opt)
					i+=1
					if self.SS_TYPE=="cmb_non_vst_SE":
						hR2err[obsnum,R2,seq,vst,var].Draw("hist e same")
				#! Draw hR2_obs and hR2_err here also
				hR2_obs[R2,seq,vst,var].Draw("same")
				hR2_err[R2,seq,vst,var].Draw("hist e same")
							
				#! Draw TLine only if ymin<0 for any of hR2 for various obsnum
				for obsnum in self.OBSD.keys():
					if hR2[obsnum,R2,seq,vst,var].GetMinimum()<0:
						ln[pad-1]=ROOT.TLine(hR2[obsnum,R2,seq,vst,var].GetXaxis().GetXmin(),0,hR2[obsnum,R2,seq,vst,var].GetXaxis().GetXmax(),0)
						ln[pad-1].Draw("same")
						break #! once a line is drawn, no more needed; break out of the loop		

			outdir_q2bin=os.path.join(self.OUTDIR_R2_ERR_VSL,"q%s"%q2bin)
			outdir_R2=os.path.join(outdir_q2bin,R2)
			if not os.path.exists(outdir_R2):
				os.makedirs(outdir_R2)
			c3.SaveAs("%s/c_w%s_q%s.png"%(outdir_R2,wbin,q2bin))
			c3.SaveAs("%s/c_w%s_q%s.eps"%(outdir_R2,wbin,q2bin))
			c3.Close()

		#! div. Save full results in text file
		#! + EF,vst,var -> OUTDIR_R2_TXT_RSLT/q2bin/wbin/
		#! + EC,vst,var -> OUTDIR_R2_TXT_RSLT/EC/q2bin/wbin/
		for item in list(itertools.product(R2S,self.SEQS)):
			R2,seq=item[0],item[1]
			#! create outdir as per seq
			if seq=="EC":
				outdir_q2bin_wbin=os.path.join(self.OUTDIR_R2_TXT_RSLT,"EC","q%s"%q2bin,"w%s"%wbin)
			elif seq=="EF":
				outdir_q2bin_wbin=os.path.join(self.OUTDIR_R2_TXT_RSLT,"q%s"%q2bin,"w%s"%wbin)
			else:
				print "ERROR: plot_SE_results_R2: seq=%s not recognized while written txt results"%seq
			outdir_R2=os.path.join(outdir_q2bin_wbin,R2)
			if not os.path.exists(outdir_R2):
				os.makedirs(outdir_R2)
			#! now write
			for item in pad_map:
				pad,vst,var=item[0],item[1],item[2]
				if self.EXTRACT_D_E_FOR_NON_ALPHA==False:
					if (R2=='D' or R2=='E') and var!= 'ALPHA':continue
				#ftxt=open("%s/%s_%s.txt"%(outdir_R2,VAR_NAMES_PLAIN[vst,var],VAR_NAMES_PLAIN[(vst,'PHI')]),"w")
				ftxt=open("%s/%d_%s.txt"%(outdir_R2,vst,var),"w")
				#! Loop over rslts per bin in rslts and write to file
				for d in rslts[R2,seq,vst,var]:
					bin,bin_le,bin_ue,mu,sg_mu,sg,sg_sg,rel_err,sg_rel_err=d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8]
					ftxt.write("%d,%f,%f = (%f +/- %f),(%f +/- %f),(%f +/- %f)\n"%(bin,bin_le,bin_ue,mu,sg_mu,sg,sg_sg,rel_err,sg_rel_err))
				ftxt.close()

		#! 3. EC-EH comparison
		for R2 in R2S:
			c1=ROOT.TCanvas("c1","c1",1000,1000)
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

			print "DispObs::plot_SE_results_R2() Plotting hR2:EC-EH comparison for q2wbin=%s,R2=%s"%(q2wbin,R2)

			#! TLine, for each pad, for ln at y=0
  			ln=[0 for i in range(npads)]
			for item in pad_map:
				pad,vst,var=item[0],item[1],item[2]
				print "pad,vst,var=",pad,vst,var
				gpad=pad_p.cd(pad)

				if self.EXTRACT_D_E_FOR_NON_ALPHA==False:
					if (R2=='D' or R2=='E') and var!= 'ALPHA':continue
				for i,seq in enumerate(['EF','EC']): #! self.SEQS purposefully modified thus to plot 'EF' first
					#! Setup aesthetics for EC-EH comparison
					hR2_obs[R2,seq,vst,var].SetMarkerColor(CLRS_PLT_SEQ_ALL[seq])
					hR2_obs[R2,seq,vst,var].SetLineColor(CLRS_PLT_SEQ_ALL[seq])
					#! 3. Finally draw
					#! Make sure maximum and minimum are displayed
					#! Note that for the following the methods are commented out because default plotting looks OK
					# self.set_h1_maximum_minimum(hR2_obs[R2,seq,vst,var],10/100)
					# self.set_h1_maximum_minimum(hR2_err[R2,seq,vst,var],10/100,minimum=0)
					draw_opt="" if i==0 else "same"
					hR2_obs[R2,seq,vst,var].Draw(draw_opt)
					
							
					#! Draw TLine only if ymin<0 for any of hR2 for various obsnum
					if hR2_obs[R2,seq,vst,var].GetMinimum()<0:
						ln[pad-1]=ROOT.TLine(hR2_obs[R2,seq,vst,var].GetXaxis().GetXmin(),0,hR2_obs[R2,seq,vst,var].GetXaxis().GetXmax(),0)
						ln[pad-1].Draw("same")

				#! Draw TLegend if pad==1 (since this pad contain Mppip distribution that tapers off towards the right)
				if pad==1:
					l=ROOT.TLegend(0.7,0.75,0.9,0.90)
					for seq in self.SEQS:
						l.AddEntry(hR2_obs[R2,seq,vst,var],seq,"p")
					l.Draw()	
			
			outdir_q2bin=os.path.join(self.OUTDIR_R2,"EC-EH","q%s"%q2bin)
			outdir_R2=os.path.join(outdir_q2bin,R2)
			if not os.path.exists(outdir_R2):
				os.makedirs(outdir_R2)
			c1.SaveAs("%s/c_w%s_q%s.png"%(outdir_R2,wbin,q2bin))
			c1.SaveAs("%s/c_w%s_q%s.eps"%(outdir_R2,wbin,q2bin))
			c1.Close()

		#! 4. Results stored in FOUT: [EC,EH]* [hR2_obs,hR2_err,hR2_rel_err]
		#! output structure: q2wbin/hR2_obs(err/rel_err)_seq_vst_var
		fout_q2wdir=self.FOUT_OBS_R2.GetDirectory(q2wbin)
		if fout_q2wdir==None:
			fout_q2wdir=self.FOUT_OBS_R2.mkdir(q2wbin)
		fout_q2wdir.cd()
		#! Now save [EC,EH]*hR2_obs,hR2_err,hR2_rel_err	
		for R2 in R2S:
			for item in pad_map:
				pad,vst,var=item[0],item[1],item[2]
				print "pad,vst,var=",pad,vst,var
				
				if self.EXTRACT_D_E_FOR_NON_ALPHA==False:
					if (R2=='D' or R2=='E') and var!= 'ALPHA':continue
				for i,seq in enumerate(self.SEQS): 
					#! Save hists
					#! + Note the hists have to be named, again because:
					#! 		i. somewhere earlier their titles are reset to "" 
					#!		AND, more importantly,
					#!		ii. 'obs','err', and 'rel_err' need to be added to their names
					hR2_obs[R2,seq,vst,var].SetName("hR2_obs_%s_%s_%d_%s"%(R2,seq,vst,var))
					hR2_obs[R2,seq,vst,var].Write()
					hR2_err[R2,seq,vst,var].SetName("hR2_err_%s_%s_%d_%s"%(R2,seq,vst,var))
					hR2_err[R2,seq,vst,var].Write()
					if seq=="EC": #! then set limits on rel_err (note this is already done for EF for canvas displays, where the upper limit of 200 was established)
						#self.set_h1_maximum_minimum(hR2_rel_err[R2,seq,vst,var],10/100,minimum=0, maximum=200)
						#! Keep "blown-up" error scale
						self.set_h1_maximum_minimum(hR2_rel_err[R2,seq,vst,var],10/100,minimum=0)
					hR2_rel_err[R2,seq,vst,var].SetName("hR2_rel_err_%s_%s_%d_%s"%(R2,seq,vst,var))
					hR2_rel_err[R2,seq,vst,var].Write()
				

		print "Done DispObs::plot_SE_results_R2()"
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
				hR2d[seq,vst,var].SetName("hR2_%s_%s_%d_%s"%(R2,seq,vst,var))
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
		
		#! Call Sumw2() for hR2s since these errors will need to be propagated
		#! [09-25-16]:
		#! + Note that hitherto, this was done after applying EF*SF
		#!   (However this change should not make a difference: see below)
		#! + However, this shold not make a difference because calling TH1::SetBinError(),
		#!   calls TH1::Sumw2(), and therefore, at this point it should already be set.
		for k in hR2d:
			hR2d[k].Sumw2()

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
		
		#! Normalize and apply_rad_eff_corr to hR2
		self.norm_1D(hR2d,q2wbin)
		if self.APPLY_RAD_EFF_CORR==True:
			self.apply_rad_eff_corr(hR2d,q2wbin)
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

		#! For saving output to .root file, create and cd into q2wbindir
		#! Not saving phiprojs to .root file for now
		# if self.FOUT_R2 !=None:
		# 	q2wbin_dir=self.FOUT_R2.GetDirectory(q2wbin)
		# 	if q2wbin_dir==None: 
		# 		q2wbin_dir=self.FOUT_R2.mkdir(q2wbin)
		# phiproj_dir=q2wbin_dir.mkdir("phiprojs")
		# phiproj_dir.cd()

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
							elif self.VIEW=="EC_ST":
								if   seq=='EC': pad.cd(1)
								elif seq=='ST': pad.cd(2)
							elif self.VIEW=="EC_EF_ST":
								if   seq=='EC': pad.cd(1)
								elif seq=='EF': pad.cd(2)
								elif seq=='ST': pad.cd(3)
							
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

	# def set_h1_maximum(self,h1,pcnt=10/100):
	# 	#! Get current maximum
	# 	#! The following is susceptible to getting maximum being overridden due to various plotting related issues
	# 	#mxm=h1.GetMaximum()
	# 	#! To get the maximum value of bins in the histogram regardless of whether the value has been overridden
	# 	#! + from https://root.cern.ch/doc/master/classTH1.html#acb53c99a65ab29a045cbdc2789e55250
	# 	mxm=h1.GetBinContent(h1.GetMaximumBin())
	# 	#print "old=",mxm
	# 	#! Set new maximum=maximum+(pct*maximum)
	# 	mxm=mxm+(pcnt*mxm)
	# 	h1.SetMaximum(mxm)
	# 	#print "new=",h1.GetMaximum()

	def set_h1_maximum_minimum(self,h1,pcnt=10/100,**kwargs):
		'''
		+ [06-28-17] Note that this method was updated to use the most generalized
		  way of setting maximum and minimum by using the bin-content and bin-error
		  of each bin, as implemented in a similar method today which takes in a list of h1s: set_h1l_maximum_minimum()
		+ Set maximum and minimum for h1 so that both are displayed based:
			+ on pcnt
			OR
			+ 'minimum' and 'maximum' as kwargs, which overrides value based on pcnt

		+ Note that the direct methods GetMaximum/Minimum are not used because they are susceptible to 
		  obtaining overridden values from previous plotting related code and therefore the following is used:
		  	+ mxm/mnm=h1.GetBinContent(h1.GetMaximumBin()/GetMinimumBin())
		  	(from https://root.cern.ch/doc/master/classTH1.html#acb53c99a65ab29a045cbdc2789e55250)
		'''
		#! Get current maximum and minimum
		yvals=[]
		nbins=h1.GetNbinsX()
		for ibin in range(nbins):
			binc=h1.GetBinContent(ibin+1)
			bine=h1.GetBinError(ibin+1)
			yvals.append(binc+bine)
			yvals.append(binc-bine)
		mnm=min(yvals)
		if mnm>0: mnm=0
		mxm=max(yvals)
		#! 10% padding
		if mnm!=0:
			mnm=mnm-(pcnt*mnm)
		mxm=mxm+(pcnt*mxm)
		h1.SetMaximum(mxm)
		h1.SetMinimum(mnm)
		#! [06-28-17] hitherto way
		# #! Get current maximum and minimum
		# mxm=h1.GetBinContent(h1.GetMaximumBin())
		# mnm=h1.GetBinContent(h1.GetMinimumBin())
		# #! Set new maximum/minimum=maximum/minimum +/-(pct*maximum/minimum)
		# mxm=mxm+(pcnt*mxm)
		# mnm=mnm-(pcnt*mnm)
		# h1.SetMaximum(mxm)
		# h1.SetMinimum(mnm)

		#! Override minimum and maximum
		if 'maximum' in kwargs.keys():
			h1.SetMaximum(kwargs['maximum'])
		if 'minimum' in kwargs.keys():
			h1.SetMinimum(kwargs['minimum'])

	def get_binw(self,h1,var):
		'''
		Return binw of h1 

		+ Based on code in fill_h2_itg_yld()
		'''

		if var=='M1' or var=='M2': #! then binw=uniform and directly encoded as histogram's binw
			binw=h1.GetBinWidth(1)
		elif var=='THETA': #! binw=variable and not that can be directly obtained from TH1::GetBinWidth(). Therefore, for now, return <binw>
			binw=0
			nbins=h1.GetNbinsX()
			for ibin in range(nbins):
				theta_a=h1.GetBinLowEdge(ibin+1)
				theta_b=h1.GetBinLowEdge(ibin+2)
				DCosTheta=math.fabs(math.cos(math.radians(theta_b))-math.cos(math.radians(theta_a)))
				binw+=DCosTheta
			#! compute average
			binw=binw/nbins	
		elif var=='ALPHA': #! binw=variable and not that can be directly obtained from TH1::GetBinWidth(). Therefore, for now, return <binw>
			binw=0
			nbins=h1.GetNbinsX()
			for ibin in range(nbins):
				alpha_a=h1.GetBinLowEdge(ibin+1)
				alpha_b=h1.GetBinLowEdge(ibin+2)
				DAlpha=math.fabs(math.radians(alpha_b)-math.radians(alpha_a))
				binw+=DAlpha
			#! compute average
			binw=binw/nbins
		return binw

	def set_h1l_maximum_minimum(self,h1l):
		'''
		+ From a list of h1s, obtain and set the minimum and maximum limits of the y-axis
		  on which all h1s can be completely displayed.
		+ If minimum>0, then simply set minimum=0
		+ Note that this is done in the most general way by looping over each bin of every
		  h1 in the list and using the bin-content and bin-error to determine the y-axis
		  limits
		'''
		#! obtain minimum and maximum
		yvals=[]
		for h1 in h1l:
			nbins=h1.GetNbinsX()
			for ibin in range(nbins):
				binc=h1.GetBinContent(ibin+1)
				bine=h1.GetBinError(ibin+1)
				yvals.append(binc+bine)
				yvals.append(binc-bine)
		maximum=max(yvals)
		minimum=min(yvals)
		if minimum>0: minimum=0
		#! Now set minimum and maximum with 10% padding (as ROOT does too)
		#! + Note special treatment if minimum is set to 0
		pdng=10/100
		if minimum!=0:
			minimum=minimum-(pdng*minimum)
		maximum=maximum+(pdng*maximum)
		for h1 in h1l:
			h1.SetMinimum(minimum)
			h1.SetMaximum(maximum)
		return
		
