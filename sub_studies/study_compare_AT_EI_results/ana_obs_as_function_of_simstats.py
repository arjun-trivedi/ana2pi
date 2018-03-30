#!/usr/bin/python
from __future__ import division
import os,sys,datetime
from collections import OrderedDict
import array

import ROOT

from rootpy.io import root_open, DoesNotExist

import math

import numpy as np
import matplotlib.pyplot as plt

import itertools

import ana_h5_stats

'''

'''

USAGE='study_obs_as_function_of_simstats dbg[=False] simstats_show_rel_err_dist[=False] plot_h5_stats_vst_var[=False]'

#! user inputs
DBG=False
if len(sys.argv)>1: #! dbg entered by user
	if    sys.argv[1]=="True":  DBG=True
	elif  sys.argv[1]=="False": DBG=False
	else: sys.exit('DBG=%s is not valid. usage: %s'%(sys.argv[1],USAGE))

SIMSTATS_SHOW_REL_ERR_DIST=False
if len(sys.argv)>2: #!  show_rel_err_dist entered by user
	if    sys.argv[2]=="True":  SIMSTATS_SHOW_REL_ERR_DIST=True
	elif  sys.argv[2]=="False": SIMSTATS_SHOW_REL_ERR_DIST=False
	else: sys.exit('SIMSTATS_SHOW_REL_ERR_DIST=%s is not valid. usage: %s'%(sys.argv[3],USAGE))

PLOT_H5_STATS_VST_VAR=False
if len(sys.argv)>3: #!  show_rel_err_dist entered by user
	if    sys.argv[3]=="True":  PLOT_H5_STATS_VST_VAR=True
	elif  sys.argv[3]=="False": PLOT_H5_STATS_VST_VAR=False
	else: sys.exit('PLOT_H5_STATS_VST_VAR=%s is not valid. usage: %s'%(sys.argv[3],USAGE))

print "DBG=",DBG
print "SIMSTATS_SHOW_REL_ERR_DIST=",SIMSTATS_SHOW_REL_ERR_DIST
print "PLOT_H5_STATS_VST_VAR=",PLOT_H5_STATS_VST_VAR
#sys.exit()

#! imports from proc_h8.py
sys.path.insert(0, '%s/elast_lite/obs_2pi'%os.environ['ANA2PI'])
from proc_h8 import H5_DIM
H4_PROJDIM=array.array('i',[H5_DIM['M1'],H5_DIM['M2'],H5_DIM['PHI'],H5_DIM['ALPHA']])
from disp_obs import VAR_NAMES_PLAIN

#! Set up THnTool
ROOT.gROOT.ProcessLine(".L THnTool.C+")
from ROOT import THnTool
thntool=THnTool()

#! Set up general ROOT plotting aesthetics
#ROOT.gStyle.SetOptStat("nmMrReiuo")
#ROOT.gStyle.SetStatStyle(0) #! transparent stats box

#! OUTDIR
DATE=datetime.datetime.now().strftime('%m%d%y')
OUTDIRNAME='results'
if SIMSTATS_SHOW_REL_ERR_DIST==True:
	OUTDIRNAME+='_w_SA_rel_err'
OUTDIRNAME+='_%s'%DATE
#! Finally create OUTDIR
if DBG==True:
	OUTDIR=os.path.join(os.environ['STUDY_OBS_AS_FUNCTION_OF_SIMSTATS_DATADIR'],'dbg',OUTDIRNAME)
else:
	OUTDIR=os.path.join(os.environ['STUDY_OBS_AS_FUNCTION_OF_SIMSTATS_DATADIR'],OUTDIRNAME)
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)
print "OUTDIR=",OUTDIR
#sys.exit()

#! Setup input root files: FYLD[q2r][sim], FO1D[q2r][sim],FOIT[q2r][sim]
#! Setup Q2 ranges
NQ2RANGES=2
LQ2,HQ2=range(NQ2RANGES)
#! Setup simulations
ISIM=[0 for i in range(NQ2RANGES)]
ISIM[LQ2] =[0,1,2,3,4,5,6]
ISIM[HQ2]=[0,1,2,3,4]
SIMNAME=[0 for i in range(NQ2RANGES)]
SIMNAME[LQ2]=['sim4','sim4_sim5', 'sim4_sim5_sim6',  'sim4_sim5_sim6_sim7','sim4_sim5_sim6_sim7_sim8','sim4_sim5_sim6_sim7_sim8_sim13','sim4_sim5_sim6_sim7_sim8_sim13_sim14']
SIMNAME[HQ2]=['sim9','sim9_sim10','sim9_sim10_sim11','sim9_sim10_sim11_sim12','sim9_sim10_sim11_sim12_sim15']
NSIM=[0 for i in range(NQ2RANGES)]
NSIM[LQ2]=len(ISIM[LQ2])
NSIM[HQ2]=len(ISIM[HQ2])
SIMPCNTG=[0 for i in range(NQ2RANGES)]
SIMPCNTG[LQ2]=[14.3,28.6,42.8,57.1,71.4,85.7,100.0] #![16.7,33.3,50.0,66.7,83.3,100.0]
SIMPCNTG[HQ2]=[20.0,40.0,60.0,80.0,100.0]#1[25,50,75,100]
SIMFRCTN=[0 for i in range(NQ2RANGES)]
SIMFRCTN[LQ2]=[0.14,0.28,0.43,0.57,0.71,0.86,1.0]#[0.17,0.33,0.50,0.67,0.83,1.00]
SIMFRCTN[HQ2]=[0.20,0.40,0.60,0.80,1.0]#[0.25,0.50,0.75,1.00]
SIMTOT=[0 for i in range(NQ2RANGES)]
SIMTOT[LQ2]=['~7B']#['~6B']
SIMTOT[HQ2]=['~5B']#['~4B']
#! Finally setup and fill FYLD[q2r][sim],FO1D[q2r][sim],FOIT[q2r][sim]
FYLD=[[0 for i in range(NSIM[LQ2])],[0 for i in range(NSIM[HQ2])]]
FO1D=[[0 for i in range(NSIM[LQ2])],[0 for i in range(NSIM[HQ2])]]
FOIT=[[0 for i in range(NSIM[LQ2])],[0 for i in range(NSIM[HQ2])]]
#print FYLD,FYLD[LQ2],FYLD[HQ2]
#sys.exit()
#tmp_obsdir='%s/tmp/obsdir'%os.environ['GDRIVE'] #! os.environ['OBSDIR_E16'] ALSO -> !cutsncor-cutncors for lowQ2 only!
for q2 in range(NQ2RANGES):
	if q2==LQ2:
		for isim in ISIM[LQ2]:
			#if isim==0: continue #! data corrurpt. remaking. skip for now
			if   isim==0: date='032818'
			#elif isim==6: date='032918'
			else:         date='032918'
			FYLD[q2][isim]=root_open('%s/lowQ2_SSBands_off_off_%s_%s/cutsncors1/%s/yield.root'%(os.environ['OBSDIR_E16'],SIMNAME[q2][isim],date,SIMNAME[q2][isim]),'r')
			FO1D[q2][isim]=root_open('%s/lowQ2_SSBands_off_off_%s_%s/cutsncors1/%s/Obs_1D_norm/obs_1D.root'%(os.environ['OBSDIR_E16'],SIMNAME[q2][isim],date,SIMNAME[q2][isim]),'r')
			FOIT[q2][isim]=root_open('%s/lowQ2_SSBands_off_off_%s_%s/cutsncors1/%s/Obs_Itg_Yld_norm/obs_itg_yld.root'%(os.environ['OBSDIR_E16'],SIMNAME[q2][isim],date,SIMNAME[q2][isim]),'r')
	elif q2==HQ2:
		for isim in ISIM[HQ2]:
			if   isim==4: date='033018'
			#elif isim==4: date='032718'
			else:         date='032918'
			FYLD[q2][isim]=root_open('%s/highQ2_SSBands_off_off_%s_%s/cutsncors1/%s/yield.root'%(os.environ['OBSDIR_E16'],SIMNAME[q2][isim],date,SIMNAME[q2][isim]),'r')
			FO1D[q2][isim]=root_open('%s/highQ2_SSBands_off_off_%s_%s/cutsncors1/%s/Obs_1D_norm/obs_1D.root'%(os.environ['OBSDIR_E16'],SIMNAME[q2][isim],date,SIMNAME[q2][isim]),'r')
			FOIT[q2][isim]=root_open('%s/highQ2_SSBands_off_off_%s_%s/cutsncors1/%s/Obs_Itg_Yld_norm/obs_itg_yld.root'%(os.environ['OBSDIR_E16'],SIMNAME[q2][isim],date,SIMNAME[q2][isim]),'r')
if DBG==True:
	print FYLD
	print FO1D
	print FOIT
#sys.exit()

#! Set up other analysis constants
#! TCanvas's pad_map[pad,vst,var] defined as per Gleb's display
PAD_MAP_1D=[(1,1,"M1"),   (2,3,'M2'),   (3,2,'M2'),
		    (4,1,"THETA"),(5,3,'THETA'),(6,2,'THETA'),
		    (7,1,"ALPHA"),(8,3,'ALPHA'),(9,2,'ALPHA')]

#! Marker colors as per simulation stats: inspired by ROOT's color pallette:hot -> cold = low -> high
CLRS=[ROOT.gROOT.ProcessLine("kViolet-7"),
      ROOT.gROOT.ProcessLine("kBlue"),
      ROOT.gROOT.ProcessLine("kCyan"),
      ROOT.gROOT.ProcessLine("kGreen"),
      ROOT.gROOT.ProcessLine("kYellow"),
      #ROOT.gROOT.ProcessLine("kPink+1"),
      ROOT.gROOT.ProcessLine("kOrange"),
      ROOT.gROOT.ProcessLine("kRed")]

#! Marker colors as per simulation stats: to correlate with colors defined by ROOT's color pallette:hot -> cold = low -> high
CLRS_MPLT=['rebeccapurple'
           'b',
           'c',
           'g',
           'y',
           #ROOT.gROOT.ProcessLine("kPink+1"),
           'darkorange',
           'r']

#! Marker types as per sequence
MRKS_PLT_SEQ={('ST'):ROOT.gROOT.ProcessLine("kFullDotLarge"),
			  ('SR'):ROOT.gROOT.ProcessLine("kFullDotLarge"),
			  ('SA'):ROOT.gROOT.ProcessLine("kFullDotLarge"),  
			  ('ER'):ROOT.gROOT.ProcessLine("kFullDotLarge"),#!kCircle
			  ('EC'):ROOT.gROOT.ProcessLine("kFullDotLarge"),
			  ('EF'):ROOT.gROOT.ProcessLine("kFullDotLarge")}

#! PSTYPs
NPSTYP=2
AL5D,ER5D=range(NPSTYP)
PSTYPS=['AL5D','ER5D']

#! THETA binning
NTHETABINS=10
THETAMIN,THETAMAX=0,180
THETABINW=18
THETABINLE=np.arange(THETAMIN,THETAMAX,THETABINW)

def get_q2wbinlist(f,q2min=0.00,q2max=6.00,wmin=0.000,wmax=3.000,dbg=False,dbg_bins=2,
	               dbg_binl=['2.00-2.40_1.425-1.450','2.00-2.40_1.450-1.475','2.00-2.40_1.475-1.500','2.00-2.40_1.500-1.525']):
		"""
		+ Taken from disp_obs.py and modified accordingly
		+ Note in dbg mode, functions works as expected when 'dbg_bins'=number of bins in 'dbg_binl'
		    + 'dbg_binl' contains bins which need immediate analysis
		"""
		q2wbinl=[]
		
		# print "DispObs::get_q2wbinlist() Going to Q2-W bins from file=",f.GetName()
		# print "DispObs::get_q2wbinlist() q2min,q2max,wmin,wmax=",q2min,q2max,wmin,wmax
		# if dbg==True:
		# 	print "DispObs::get_q2wbinlist() dbg=True"

		#brk=False #! technical tool to break out of two nested for loops. Set in second (nested) for loop
		#if brk==True: break
		for path,dirs,files in f.walk():
			#print "path,dirs,files=",path,dirs,files
			if path=="":continue #! Avoid root path
			path_arr=path.split("/")
			if len(path_arr)==1:
				if dbg==True:
					if path in dbg_binl:
						q2wbinl.append(path)
				else:
					q2wbinl.append(path)
			if dbg==True and len(q2wbinl)==dbg_bins:
				#brk=True
				break #! Uncomment/comment -> Get limited q2w-bins/Get all q2w-bins

		return q2wbinl

def get_q2bin(q2wbin):
	return q2wbin.split('_')[0]

def norm_1D_theta(hTheta):
	#! 1. Create normalization factor histogram
	hDCosTheta=hTheta.Clone("hDCosTheta")
	hDCosTheta.SetTitle("hDCosTheta")
	hDCosTheta.Reset()
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
	return hDCosTheta

def plot_1D_athtcs():
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

def set_h1l_maximum_minimum(h1l):
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

def hist_1D_athtcs(hobs):
	'''
	'''
	for isim,h1 in enumerate(hobs):
		for k in h1.keys():
			seq,vst,var=k[0],k[1],k[2]	
			mrk=MRKS_PLT_SEQ[seq]
			clr=CLRS[isim]
			h1[k].SetMarkerStyle(mrk)
			h1[k].SetMarkerColor(clr)
			h1[k].SetLineColor(clr)

def hist_itg_athtcs(hobs):
	'''
	'''
	for isim,h1 in enumerate(hobs):
		for k in h1.keys():
			seq=k	
			mrk=MRKS_PLT_SEQ[seq]
			clr=CLRS[isim]
			h1[k].SetMarkerStyle(mrk)
			h1[k].SetMarkerColor(clr)
			h1[k].SetLineColor(clr)
			
def plot_1D(hobs,q2wbin,q2r):
	"""
	"""
	print "In plot_1D()"

	#! Set up plotting aesthetics
	plot_1D_athtcs()
	hist_1D_athtcs(hobs)

	# #! For saving output to .root file, create and cd into q2wbindir
	# if self.FOUT_1D!=None:
	# 	q2wbin_dir=self.FOUT_1D.GetDirectory(q2wbin)
	# 	if q2wbin_dir==None: 
	# 		q2wbin_dir=self.FOUT_1D.mkdir(q2wbin)
	# 	q2wbin_dir.cd()
				
	# #! TCanvas's pad_map[pad,vst,var] defined as per Gleb's display
	# pad_map=[(1,1,"M1"),   (2,3,'M2'),   (3,2,'M2'),
	# 		 (4,1,"THETA"),(5,3,'THETA'),(6,2,'THETA'),
	# 		 (7,1,"ALPHA"),(8,3,'ALPHA'),(9,2,'ALPHA')]
		
	#! A canvas each for each sequence
	for seq in ['EC','EF']:				  
		c=ROOT.TCanvas("c","c",1000,1000)
		pad_t=ROOT.TPad("pad_l","Legend pad",0.25,0.935,0.75,1.00)
		pad_p=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
		#if (self.VIEW=="fullana"):pad_p.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
		pad_p.Draw()
		pad_t.Draw()
		pad_t.cd()
		pt=ROOT.TPaveText(.05,.1,.95,.8)
		pt.AddText("%s Q2_W bin=%s"%(seq,q2wbin))
		pt.SetTextSize(0.40)
		pt.Draw()
		pad_p.Divide(3,3)
		#! TLegend
		l=ROOT.TLegend(0.53,0.60,0.90,0.90)
		l.SetFillStyle(0)
		#l.SetBorderSize(0)
		l.SetTextSize(0.05)
		l.SetHeader("sim-stats fraction")
		for item in PAD_MAP_1D:
			pad,vst,var=item[0],item[1],item[2]
			#print "pad,vst,var=",pad,vst,var
			gpad=pad_p.cd(pad)
			#! First set maximum of y axis as per all hists that go in this canvas
			nsims=len(hobs)
			histl=[]
			for isim,h1 in enumerate(hobs):
				histl.append(h1[seq,vst,var])
			set_h1l_maximum_minimum(histl)
			#! Now draw
			iplt=0
			for isim,h1 in enumerate(hobs):
				# print 'sim#=',isim+1
				# for k in h1.keys():
				# 	print k,h1[k].GetName()
				#if isim==0: continue #! "skip"
				draw_opt="" if iplt==0 else "same"
				h1[seq,vst,var].Draw(draw_opt)
				iplt+=1
			#! Draw legend
			if pad==1:
				for isim,h1 in enumerate(hobs):
					l.AddEntry(h1[seq,vst,var],"%s"%SIMFRCTN[q2r][isim],"p")
				l.Draw()
		#! Save canvas 
		outdir=os.path.join(OUTDIR,"Obs_1D",q2wbin)
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		c.SaveAs("%s/c_1D_%s.png"%(outdir,seq))
		
	print "*** plot_1D() Done ***\n"
	return

def plot_itg(hobs,q2wbin,q2r):
	"""
	"""
	print "In plot_itg()"

	#! Determine nsims
	nsims=len(hobs)
		
	#! General aesthetics
	ROOT.gStyle.SetOptStat(0)

	#! Histogram aesthetics
	hist_itg_athtcs(hobs)

	#! A canvas each for each sequence
	for seq in ['EC','EF']:				  
		c=ROOT.TCanvas("c","c",1200,800)
		pad_p=ROOT.TPad("pad_p","Plots pad",0.00,0.00,1.00,0.90)
		pad_t=ROOT.TPad("pad_p","Plots pad",0.00,0.90,1.00,1.00)
		pad_p.Draw()
		pad_t.Draw()
		pad_t.cd()
		pt=ROOT.TPaveText(.05,.1,.95,.8)
		q2bin=get_q2bin(q2wbin)
		pt.AddText("Integrated cross-section for %s,Q2=%s (vst,var=1,THETA)"%(seq,q2bin))
		pt.SetTextSize(0.32)
  		pt.Draw()
		pad_p.cd()
		#! TLegend
		l=ROOT.TLegend(0.60,0.70,0.90,0.90)
		l.SetFillStyle(0)
		#l.SetBorderSize(0)
		l.SetTextSize(0.04)
		l.SetHeader("sim-stats fraction")
		#! Set maximum of y axis as per all hists that go in this canvas
		nsims=len(hobs)
		histl=[]
		for isim,h1 in enumerate(hobs):
			histl.append(h1[seq])
		set_h1l_maximum_minimum(histl)
		#! Begin plotting
		for isim in range(nsims):
			draw_opt="" if isim==0 else "same"
			#! Set up some aesthetics particular to hW
			#! + Give already rotated bin labels by 90 deg. more space
			#! + Move xaxis title lower to accomodate 90 deg. rotated titles
			pad_p.SetBottomMargin(0.20)
			hobs[isim][seq].GetXaxis().SetTitleOffset(2.80)
			hobs[isim][seq].GetXaxis().SetTitleOffset(2.80)
			#! Draw
			hobs[isim][seq].SetTitle("")
			hobs[isim][seq].Draw(draw_opt)
		#! Legend
		for isim in range(nsims):
			l.AddEntry(hobs[isim][seq],"%s"%SIMFRCTN[q2r][isim],"p")
		l.Draw()
		#! Save canvas 
		outdir=os.path.join(OUTDIR,"Obs_itg",q2bin)
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		c.SaveAs("%s/c_itg_%s.png"%(outdir,seq))
		
	print "*** plot_itg() Done ***\n"
	return
#! main

#! debug
#print FYLD[LQ2][0]
#print FYLD[HQ2][0]
#print FO1D[LQ2][1]
#sys.exit()

#! First get Q2WBINL from ISIM=0 (any ISIM should get the same list)
Q2WBINL=[0 for i in range(NQ2RANGES)]
if DBG==True:
	#Q2WBINL[LQ2]=get_q2wbinlist(FYLD[LQ2][0],dbg=True,dbg_bins=2,dbg_binl=['2.00-2.40_1.500-1.525','2.40-3.00_1.500-1.525'])
	Q2WBINL[LQ2]=get_q2wbinlist(FYLD[LQ2][0],dbg=True,dbg_bins=1,dbg_binl=['2.00-2.40_1.500-1.525'])
	#Q2WBINL[HQ2]=get_q2wbinlist(FYLD[HQ2][0],dbg=True,dbg_bins=1,dbg_binl=['3.00-3.50_1.500-1.525'])
else:
	#Q2WBINL[LQ2]=get_q2wbinlist(FYLD[LQ2][0])
	#Q2WBINL[HQ2]=get_q2wbinlist(FYLD[HQ2][0])
	# Q2WBINL[LQ2]=get_q2wbinlist(FYLD[LQ2][0],dbg=True,dbg_bins=8,dbg_binl=['2.00-2.40_1.425-1.450','2.00-2.40_1.450-1.475','2.00-2.40_1.475-1.500','2.00-2.40_1.500-1.525',
	# 	                                                                   '2.40-3.00_1.425-1.450','2.40-3.00_1.450-1.475','2.40-3.00_1.475-1.500','2.40-3.00_1.500-1.525'])
	# Q2WBINL[HQ2]=get_q2wbinlist(FYLD[HQ2][0],dbg=True,dbg_bins=8,dbg_binl=['3.50-4.20_1.425-1.450','3.50-4.20_1.450-1.475','3.50-4.20_1.475-1.500','3.50-4.20_1.500-1.525',
	# 	                                                                   '4.20-5.00_1.425-1.450','4.20-5.00_1.450-1.475','4.20-5.00_1.475-1.500','4.20-5.00_1.500-1.525'])	                                                                	
	Q2WBINL[LQ2]=get_q2wbinlist(FYLD[LQ2][0])
	Q2WBINL[HQ2]=get_q2wbinlist(FYLD[HQ2][0])
if DBG==True:
	print "Q2WBINL[LQ2]=",Q2WBINL[LQ2]
	print "Q2WBINL[HQ2]=",Q2WBINL[HQ2]
#sys.exit()

#! 1. Loop over Q2WBINL and for each q2wbin make Obs_1D=f(sim)
for q2r in range(NQ2RANGES):
	if isinstance(Q2WBINL[q2r],list)==False: continue #! in debug mode Q2WBINL[q2r] may not be obtained
	for q2wbin in Q2WBINL[q2r]:
		#! Create structure to hold hobs[isim][seq,vst,var]
		hobs=[OrderedDict() for isim in range(NSIM[q2r])]
		for isim in ISIM[q2r]:
			#if q2r==LQ2 and isim==0: continue #! data corrurpt. remaking. skip for now
			#! Get hobs
			for item in PAD_MAP_1D:
				pad,vst,var=item[0],item[1],item[2]
				for seq in ['EC','EF']:
					hobs[isim][seq,vst,var]=FO1D[q2r][isim].Get("%s/h1_%s_%d_%s"%(q2wbin,seq,vst,var))
		if DBG==True:
			print "***hobs pretty print for q2wbin",q2wbin,"***"
			for isim,h1 in enumerate(hobs):
				print "** sim#=",isim+1,"**"
				for k in h1:
					print k,h1[k].GetName()
		#! plot
		plot_1D(hobs,q2wbin,q2r)

#! 2. Plot Obs_itg=f(sim) 
#! + use vst=1, var=THETA
#! First get Q2WBINL (different format compared to Obs_1D) from ISIM=0 (any ISIM should get the same list)
Q2WBINL_ITG=[0 for i in range(NQ2RANGES)]
#! Note that no debug mode needed since there only a few bins
Q2WBINL_ITG[LQ2]=get_q2wbinlist(FOIT[LQ2][0])
Q2WBINL_ITG[HQ2]=get_q2wbinlist(FOIT[HQ2][0])	
print "Q2WBINL_ITG[LQ2]=",Q2WBINL_ITG[LQ2]
print "Q2WBINL_ITG[HQ2]=",Q2WBINL_ITG[HQ2]
for q2r in range(NQ2RANGES):
	for q2wbin in Q2WBINL_ITG[q2r]:
		#! Create structure to hold hobs[isim][seq,vst,var]
		hobs=[OrderedDict() for isim in range(NSIM[q2r])]
		for isim in ISIM[q2r]:
			#! Get hobs
			for seq in ['EC','EF']:
				hobs[isim][seq]=FOIT[q2r][isim].Get("%s/hW_%s_1_THETA"%(q2wbin,seq))
		if DBG==True:
			print "***hobs pretty print for q2wbin",q2wbin,"***"
			for isim,h1 in enumerate(hobs):
				print "** sim#=",isim+1,"**"
				for k in h1:
					print k,h1[k].GetName()
		#! plot
		plot_itg(hobs,q2wbin,q2r)

#! 3. Loop over Q2WBINL and for each q2wbin make simstats=f(sim)
for q2r in range(NQ2RANGES):
	if isinstance(Q2WBINL[q2r],list)==False: continue #! in debug mode Q2WBINL[q2r] may not be obtained
	for q2wbin in Q2WBINL[q2r]:
		#! Create structure to hold h5[isim][seq] (note vst=VST1 for now)
		h5l=[OrderedDict() for isim in range(NSIM[q2r])]
		for isim in ISIM[q2r]:
			#! Get h5
			for seq in ['ST','SR','SA','ER']:
				h5l[isim][seq]=FYLD[q2r][isim].Get("%s/%s/VST%d/h5"%(q2wbin,seq,1))
		if DBG==True:
			print "***h5 pretty print for q2wbin",q2wbin,"***"
			for isim,h5 in enumerate(h5l):
				print "** sim#=",isim+1,"**"
				for k in h5:
					print k,h5[k].GetName()
		#! plot
		#! The following 2 steps take a while, especially plot_simstats_5D()
		ana_h5_stats.plot_h5_stats(h5l,'SA',q2wbin,OUTDIR,show_rel_err_dist=SIMSTATS_SHOW_REL_ERR_DIST)
		if PLOT_H5_STATS_VST_VAR:
			ana_h5_stats.plot_h5_stats_vst_var(h5l,'SA',q2wbin,OUTDIR,show_rel_err_dist=SIMSTATS_SHOW_REL_ERR_DIST)