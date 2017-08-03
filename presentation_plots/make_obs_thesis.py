#!/usr/bin/python
from __future__ import division

import os,sys,subprocess
import ROOT
from rootpy.io import root_open, DoesNotExist
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict
import itertools

import numpy as np

import math

import time,datetime

'''
[10-24-16]
+ Makes plots that will be finally used in the thesis. Idea is to do further processing of observables that will be
  displayed in the thesis:
  	1. Final results with systematic effects are processed in the following way (see next bullet point for why):
  		+ Only for 'Obs_itg', use results from the full process of systematic studies developed
  		+ For 'Obs_1D' and 'Obs_R2', use results from 'SSBands/cutsncors1' and propagate systematic errors estimated
  		  for 'Obs_itg' using full chain of results in 'cmb_vst_SE_<date>'
			+ Adjust error on systematic error: currently setting it to 0.
	+ Why plot systematic effects in this way?
		+ This is because only for 'Obs_itg' are there enough statistics to determine precisely the error due to systematic effects i.e. the variations due to different systematic effects are beyond statistical errors. Therefore the determined error on the systematic error is low.
		+ For 'Obs_1D' and 'Obs_R2', the variations are well within statistical errors because the statistics are much lower, which leads to unacceptable levels of precision on the determined systematic effects i.e. the error on the systematic error is too high, reaching levels of 100%. Therefore, it was decided to use the systematic effects as determined from 'Obs_itg' as a measure for systematic effects for these observables.

	2. Display ST, EC and EF together for 'Obs_1D' and 'Obs_R2':
		+ see which normalization scheme works best for displaying ST: currently normalizing in a way to compare shape with EF

+ Two different normalization for ST histograms can be used though the former is preferred: "shape" and "shape_and_amplitude"

+ Note that the main code follows a 3 step process, each within its own loop (inefficient, but practical for now)
	step 1. First process Obs_itg from 'cmb_vst_SE' and get SYST_ERR
	step 2. If NORM_ST=="shape_and_amplitude", get NORM_FACTOR(q2-sim) using Obs_1D from 'SSBands/cutsncors1'
	step 3. Process Obs_1D and Obs_R2 from 'SSBands/cutsncors1'
		+ In this step normalization for ST is done as per NORM_ST

+ >make_obs_thesis.py [dbg=False] [norm_ST=shape] [study_norm_factor_only=False]
'''

#! Get input from user
if len(sys.argv)>1: #! User input dbg, else use dbg=False
	if sys.argv[1]=="False":
		DBG=False
	elif sys.argv[1]=="True":
		DBG=True
	else:
		print "Please enter dbg=True or False!"
		sys.exit()
else:
	DBG=False

if len(sys.argv)>2: #! User input norm_ST, else use norm_ST=shape
	if sys.argv[2]=="shape":
		NORM_ST="shape"
	elif sys.argv[2]=="shape_and_amplitude":
		NORM_ST="shape_and_amplitude"
	else:
		print "Please enter norm_ST=shape or shape_and_amplitude!"
		sys.exit()
else:
	NORM_ST="shape"	


if len(sys.argv)>3: #! User input study_norm_factor_only, else use study_norm_factor_only=False
	if sys.argv[3]=="False":
		STUDY_NORM_FACTOR_ONLY=False
	elif sys.argv[3]=="True":
		STUDY_NORM_FACTOR_ONLY=True
	else:
		print "Please enter study_norm_factor_only=True or False!"
		sys.exit()
else:
	STUDY_NORM_FACTOR_ONLY=False

print "DBG=",DBG
print "NORM_ST=",NORM_ST
print "STUDY_NORM_FACTOR_ONLY=",STUDY_NORM_FACTOR_ONLY
#sys.exit()

OBSDIR_E16=os.environ['OBSDIR_E16']
OBSDIR_E162=os.environ['OBSDIR_E162']
OBS=["1D","itg","R2"]
DATE=date=datetime.datetime.now().strftime('%m%d%y')

#! Create OUTDIR
OUTDIR="%s/thesis_obs_norm_ST_%s_%s/"%(OBSDIR_E162,NORM_ST,DATE)
if DBG==True: OUTDIR="%s/thesis_obs_norm_ST_%s_dbg_%s/"%(OBSDIR_E162,NORM_ST,DATE)
print "OUTDIR=",OUTDIR

#! Setup files that will be used to obtain ST data
F_ST=OrderedDict()
F_ST['1D',"lowQ2"] ="%s/lowQ2_SSBands_080217/cutsncors1/sim4_sim5_sim6_sim7_sim8_sim13/Obs_1D_norm_EC_EF_SF/obs_1D.root"%OBSDIR_E16 #! lowQ2_SSBands_071017, lowQ2_SSBands_061217, lowQ2_SSBands_092516, lowQ2_SSBands_061417
F_ST['1D',"highQ2"]="%s/highQ2_SSBands_080217/cutsncors1/sim9_sim10_sim11_sim12/Obs_1D_norm_EC_EF_SF/obs_1D.root"%OBSDIR_E16        #! highQ2_SSBands_071017,highQ2_SSBands_061217, highQ2_SSBands_092516, highQ2_SSBands_061417 

F_ST['itg',"lowQ2"]=None
F_ST['itg',"highQ2"]=None

F_ST['R2',"lowQ2"] ="%s/lowQ2_SSBands_080217/cutsncors1/sim4_sim5_sim6_sim7_sim8_sim13/Obs_R2_EC_EF_ST/mthd_phi-proj-fit_NQ/obs_R2.root"%OBSDIR_E16 #! lowQ2_SSBands_071017, lowQ2_SSBands_061217, lowQ2_SSBands_092516, lowQ2_SSBands_061417
F_ST['R2',"highQ2"]="%s/highQ2_SSBands_080217/cutsncors1/sim9_sim10_sim11_sim12/Obs_R2_EC_EF_ST/mthd_phi-proj-fit_NQ/obs_R2.root"%OBSDIR_E16        #! highQ2_SSBands_071017,highQ2_SSBands_061217, highQ2_SSBands_092516, highQ2_SSBands_061417

#! Setup files that will be used to get results
F_RSLT=OrderedDict()
F_RSLT['1D',"lowQ2"] =F_ST['1D',"lowQ2"]
F_RSLT['1D',"highQ2"]=F_ST['1D',"highQ2"]

F_RSLT['itg',"lowQ2"] ="%s/SS/lowQ2_cmb_vst_SE_080317/Obs_itg.root"%OBSDIR_E16  #!lowQ2_cmb_vst_SE_071117, lowQ2_cmb_vst_SE_061317,lowQ2_cmb_vst_SE_092716,  lowQ2_cmb_vst_SE_061517
F_RSLT['itg',"highQ2"]="%s/SS/highQ2_cmb_vst_SE_080317/Obs_itg.root"%OBSDIR_E16 #!highQ2_cmb_vst_SE_071117,highQ2_cmb_vst_SE_061317,highQ2_cmb_vst_SE_092716, highQ2_cmb_vst_SE_061517

F_RSLT['R2',"lowQ2"] =F_ST['R2',"lowQ2"]
F_RSLT['R2',"highQ2"]=F_ST['R2',"highQ2"]

#! Setup Q2,W limits
Q2MIN,Q2MAX=OrderedDict(),OrderedDict()
Q2MIN['lowQ2']=2.00
Q2MAX['lowQ2']=3.00
Q2MIN['highQ2']=3.00
Q2MAX['highQ2']=5.00
WMIN=1.400
WMAX=2.125

#! Define some objects that will be needed in steps 1 and step 2
#! PAD_MAP[obs]=(pad,vst,var)
PAD_MAP=OrderedDict()
PAD_MAP['1D']=[(1,1,"M1"),   (2,3,'M2'),   (3,2,'M2'),
			   (4,1,"THETA"),(5,3,'THETA'),(6,2,'THETA'),
			   (7,1,"ALPHA"),(8,3,'ALPHA'),(9,2,'ALPHA')]

PAD_MAP['R2']=[(1,1,"M1"),    (2,3,'M1'),    (3,2,'M1'), 
			   (4,1,"M2"),    (5,3,'M2'),    (6,2,'M2'),
			   (7,1,"THETA"), (8,3,'THETA'), (9,2,'THETA'),
			   (10,1,"ALPHA"),(11,3,'ALPHA'),(12,2,'ALPHA')]

R2L=['A','B','C','D','E']
R2_NAMED={'A':'R2_{T}+R2_{L}','B':'R2^{c}_{LT}','C':'R2^{c}_{TT}','D':'R2^{s}_{LT}','E':'R2^{s}_{TT}'}

CLRS_SEQ={('ST'):ROOT.gROOT.ProcessLine("kRed"),
		  ('EC'):ROOT.gROOT.ProcessLine("kCyan"),
		  ('EF'):ROOT.gROOT.ProcessLine("kBlue")}

MRKS_SEQ={('ST'):ROOT.gROOT.ProcessLine("kFullDotLarge"),
	      ('EC'):ROOT.gROOT.ProcessLine("kFullDotLarge"),
		  ('EF'):ROOT.gROOT.ProcessLine("kFullDotLarge")}

#! setup functions that will be needed
def get_q2wbinlist(f,q2min=0.00,q2max=6.00,wmin=0.000,wmax=3.000,dbg=False,dbg_bins=10):
		"""
		The ROOT file is arranged in a Tree Structure. The
		Observable histograms (obs-hists) are located as files in the following directory-path(dirpath):
		q2wbin/vst/seq/hists
		"""
		q2wbinl=[]
		
		print "get_q2wbinlist() Going to Q2-W bins from file=",f.GetName()
		print "get_q2wbinlist() q2min,q2max,wmin,wmax=",q2min,q2max,wmin,wmax
		if dbg==True:
			print "get_q2wbinlist() dbg=True"

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

def get_q2bin(q2wbin):
	return q2wbin.split('_')[0]
def get_wbin(q2wbin):
	return q2wbin.split('_')[1]

def set_h1_maximum_minimum(h1,pcnt=10/100,**kwargs):
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
	#! Now set minimum and maximum with 10% padding 
	#! + Note special treatment if minimum is set to 0
	pdng=10/100
	if minimum!=0:
		minimum=minimum-(pdng*minimum)
	maximum=maximum+(pdng*maximum)
	#! now set limits
	for h1 in h1l:
		h1.SetMinimum(minimum)
		h1.SetMaximum(maximum)
	return

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

#def plot_and_write_Obs_itg(q2wbin,outdir,hobs,herr):
def plot_and_write_Obs_itg(q2wbin,hobs,herr,outdir,froot):
	#! Get q2bin information
	q2bin=get_q2bin(q2wbin)

	#! 1. Plot both seq together
	ROOT.gStyle.SetOptStat(0)
	c=ROOT.TCanvas("c","c",1200,800)
	pad_p=ROOT.TPad("pad_p","Plots pad",0.00,0.00,1.00,0.90)
	pad_t=ROOT.TPad("pad_p","Plots pad",0.00,0.90,1.00,1.00)
	pad_p.Draw()
	pad_t.Draw()

	pad_t.cd()
	pt=ROOT.TPaveText(.05,.1,.95,.8)
	pt.AddText("Integrated cross-section for Q2=%s"%(q2bin))
  	pt.SetTextSize(0.32)
  	pt.Draw()

	pad_p.cd()
	#!get rid of X error bars and y error bar caps
	ROOT.gStyle.SetErrorX(0.001)
	#! Set up some aesthetics particular to hW
	#! + Give already rotated bin labels by 90 deg. more space
	#! + Move xaxis title lower to accomodate 90 deg. rotated titles
	pad_p.SetBottomMargin(0.20)
	#! Set minimum and maximum
	set_h1_maximum_minimum(hobs['EF'],10/100,minimum=0)
	#! Marker aesthetics
	hobs['EC'].SetLineColor(ROOT.gROOT.ProcessLine("kCyan"))
	hobs['EC'].SetMarkerColor(ROOT.gROOT.ProcessLine("kCyan"))
	hobs['EF'].SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
	hobs['EF'].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlue"))
	#! Draw obs
	hobs['EF'].Draw()
	hobs['EC'].Draw("same")
	#! Draw err
	herr['EF'].SetFillColor(ROOT.gROOT.ProcessLine("kBlack"))
	herr['EF'].SetFillStyle(3003)
	herr['EF'].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
	herr['EF'].SetMarkerStyle(ROOT.gROOT.ProcessLine("kDot"))
	herr['EF'].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
	herr['EF'].Draw("hist e same")
	#! Draw TLegend 
	l=ROOT.TLegend(0.7,0.80,0.9,0.90)
	for seq in ['EC','EF']:
		l.AddEntry(hobs[seq],seq,"p")
	l.AddEntry(herr['EF'],"syst-err","l")
	l.Draw()
	#! Save canvas
	c.SaveAs("%s/c_q%s.png"%(outdir,q2wbin))
	#c.SaveAs("%s/c_q%s.eps"%(outdir,q2wbin))
	#! Convert .png->.eps->.pdf since in .eps and .pdf produced directly by ROOT
	#! the y-axis titles are cutoff
	#! 1. .png->.eps
	cmd=['convert',"%s/c_q%s.png"%(outdir,q2wbin),"%s/c_q%s.eps"%(outdir,q2wbin)]
	tool=subprocess.Popen(cmd,stderr=subprocess.STDOUT)
	tool.wait()
	#! 2. .png->.pdf
	cmd=['convert',"%s/c_q%s.png"%(outdir,q2wbin),"%s/c_q%s.pdf"%(outdir,q2wbin)]
	tool=subprocess.Popen(cmd,stderr=subprocess.STDOUT)
	tool.wait()
	
	#! 2. Save output in .root file
	froot_q2wdir=froot.GetDirectory(q2wbin)
	if froot_q2wdir==None:
		froot_q2wdir=froot.mkdir(q2wbin)
	froot_q2wdir.cd()
	for seq in ['EC','EF']:
		hobs[seq].Write()
		herr[seq].Write()

	#! 3. txt output
	outdir_txt=("%s/txt/%s"%(outdir,q2wbin))
	if not os.path.exists(outdir_txt):
		os.makedirs(outdir_txt)

	#! Loop over rslts, per seq, per bin in rslts and write to file
	for seq in ['EC','EF']:
		ftxt=open("%s/results_%s.txt"%(outdir_txt,seq),"w")
		nbins=herr[seq].GetNbinsX()
		for ibin in range(nbins):
			bin=ibin+1

			mu   =hobs[seq].GetBinContent(ibin+1)
			sg_mu=hobs[seq].GetBinError(ibin+1)
			sg   =herr[seq].GetBinContent(ibin+1)
			sg_sg=herr[seq].GetBinError(ibin+1)

			#! rel_err,sg_rel_err
			if mu==0 or sg==0: 
				rel_err=0
				sg_rel_err=0
			else:
				rel_err=(sg/mu)*100
				sg_rel_err=(rel_err)*math.sqrt((sg_sg/sg)**2+(sg_mu/mu)**2)

			#! Get wbin information that will be needed to label data
			bin_label=herr[seq].GetXaxis().GetBinLabel(ibin+1)
			bin_le=float(bin_label.split(",")[0].strip("["))
			bin_ue=float(bin_label.split(",")[1].strip(")"))
			wbin="%.3f-%.3f"%(bin_le,bin_ue)

			#SYST_ERR[q2bin,wbin]=rel_err
			#! Add "rest of systematic error estimates" to 'rel_err': see thesis's Systematic Error Table:
			#! (Assume that error on these systematic error estimates=0 i.e. 'sg_rel_err' for these errors = 0)
			#!
			#! [pre 08-03-17] Hitherto, i.e. before getting close to final ananote results, I did not re-check this values
			#! because I did not expect them to change and they did not in any significant way
			#! I Errors estimated in h10->per_non_vst_SE->cmb_non_vst_SE->cmb_vst_SE process ('rel_err' in this code):
			#!   1. Hadron Id (SS-Bands)=2.5
			#!   2. Event Selection (MM-cut)=2.9
			#!   3. Variable set dependent extraction of Obs=5.3
			#!   --> Total = sqrt(2.5**2+5.3**2+2.9**2)=6.5
			#! II Rest of the Errors:
			#!   1. Electron Id: 1
			#!   2. Electron fiducial boundary selection: 1
			#!   3. Hadron fiducial boundary selection: 1
			#!   4. Detector inefficiency identification: 1
			#!   5. Momentum and Energy Loss Corrections: 1
			#!   6. Acceptance Calculation: 1
			#!   7. Radiative Effects correction: 5
			#!   8. Estimation of Experimental Yields in Kinematical Holes: 5
			#!   9. Luminosity measurement: 5
			#!   --> Total= sqrt(1**2 + 1**2 + 1**2 + 1**2 + 1**2 + 1**2 + 5**2 + 5**2 + 5**2) = 9
			#! Total(I+II)=sqrt(6.5**2 + 9**2)=11.1 (NOTE THAT THIS IS Q2WBIN AVERAGED but final observables have fully binned version)
			#!
			#! [08-03-17] 're-obtain-obs-3'
			#! Note categories of SE:
			#! 1. ana : calculated using SE analysis (<q2wbin> average number listed; location listed next to err)
			#! 2. guess-est: guess estimated
			#! 3. ref : referenced from another source (ref listed in ananote)
			#! + Errors estimated in h10->per_non_vst_SE->cmb_non_vst_SE->cmb_vst_SE process ('rel_err' in this code):
			#!   1. Hadron Id (SS-Bands):                     3.12 (ana) (lQ2,hQ2=2.83,3.41 -> avg = 3.12) ($OBSDIR_E16/SS/lowQ2_SSBands_080217_sim4_sim5_sim6_sim7_sim8_sim13_080317/SE_plots/Obs_itg/itg.png; hQ2=corresponding folder)
			#!   2. Event Selection (MM-cut):                 NA (ana) removed in 're-obtain-obs-3'
			#!   3. Variable set dependent extraction of Obs: 5.55 (ana)(lQ2,hQ2=4.88,6.22 -> avg = 5.55) ($OBSDIR_E16/SS/lowQ2_cmb_vst_SE_080317/SE_plots/Obs_itg/itg.png; hQ2=corresponging folder)
			#!   --> Total = sqrt(3.12**2 + 5.55**2)=6.36
			#! + Rest of the Errors:
			#!   1. Electron Id:                                            1 (guess-est, conservative upper bound)
			#!   2. Electron fiducial boundary selection:                   1 (guess-est, conservative upper bound)
			#!   3. Hadron fiducial boundary selection                      1 (guess-est, conservative upper bound)
			#!   4. Detector inefficiency identification:                   1 (guess-est, conservative upper bound)
			#!   5. Momentum and Energy Loss Corrections:                   2 (ana) (Estimated from $STUDY_MM_DIFF_ER_SR_DATADIR/results_080117. For full details, see Tomboy:MM_diff_ER_SR:DateLog:080117 or handwritten notes)
			#!   6. Background:                                             1 (ana) (Estimated from curve of $STUDY_EVTSEL_BG_DATADIR/evtsel_BG_072717/2.00-5.00/MM2/bg.png) 
			#!   7. Acceptance Calculation:                                 1 (guess-est, conservative upper bound)
			#!   8. Radiative Effects correction:                           5 (ref)
			#!   9. Estimation of Experimental Yields in Kinematical Holes: 5 (ana) (For details see $ANANOTE/figures/Holes)
			#!  10. Luminosity measurement:                                 5 (ref)
			#!   --> Total = sqrt(1**2 + 1**2 + 1**2 + 1**2 + 2**2 + 1**2 + 1**2 + 5**2 + 5**2 + 5**2)=9.21
			#! Total(I+II)=sqrt(6.36**2 + 9.21**2)=11.19 (NOTE THAT THIS IS Q2WBIN AVERAGED but final observables have fully binned version)  

			#! [pre 08-03-17]
			#rel_err_full=math.sqrt(rel_err**2+9**2)
			#! [08-03-17] 're-obtain-obs-3'
			rel_err_full=math.sqrt(rel_err**2+9.21**2)
			SYST_ERR[q2bin,wbin]=rel_err_full
								
			ftxt.write("%d,%f,%f = (%f +/- %f),(%f +/- %f),(%f +/- %f)\n"%(bin,bin_le,bin_ue,mu,sg_mu,sg,sg_sg,rel_err,sg_rel_err))
		ftxt.close()

def plot_and_write_Obs_1D(q2wbin,hobs,herr,outdir,froot):
	#! Get q2bin, wbin information
	q2bin=get_q2bin(q2wbin)
	wbin=get_wbin(q2wbin)

	#! General plotting aesthetics
	plot_1D_athtcs()

	#! For saving output to .root file, create and cd into q2wbindir
	q2wbin_dir_1D=froot.GetDirectory(q2wbin)
	if q2wbin_dir_1D==None: 
		q2wbin_dir_1D=froot.mkdir(q2wbin)
	q2wbin_dir_1D.cd()

	#! 1. Plot all seqs together with ST appropriately normed to EF
	c=ROOT.TCanvas("c","c",1000,1000)
	pad_t=ROOT.TPad("pad_l","Legend pad",0.25,0.935,0.75,1.00)
	pad_p=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
	pad_p.Draw()
	pad_t.Draw()
	pad_t.cd()
	pt=ROOT.TPaveText(.05,.1,.95,.8)
	pt.AddText("Q2_W bin=%s"%q2wbin)
	pt.SetTextSize(0.40)
	pt.Draw()
	pad_p.Divide(3,3)
	for item in PAD_MAP['1D']:
		pad,vst,var=item[0],item[1],item[2]
		#print "pad,vst,var=",pad,vst,var
		gpad=pad_p.cd(pad)

		#! First set minimum and maximum of y-axis
		histl=[hobs[seq,vst,var] for seq in ['EC','EF','ST']]
		#! append herr also so that it can also be used in determining and setting y-axes limits
		histl.append(herr['EF',vst,var]) 
		set_h1l_maximum_minimum(histl)
		#! [06-28-17] hitherto way
		# fctr=1.1
		# maxl=[hobs[seq,vst,var].GetMaximum() for seq in ['EC','EF','ST']]
		# maximum=max(maxl)
		# for htmp in [hobs[seq,vst,var] for seq in ['EC','EF','ST']]:
		# 	htmp.SetMinimum(0.)
		# 	htmp.SetMaximum(maximum*fctr)
		#! Now draw hobs
		for i,seq in enumerate(['ST','EC','EF']):
			draw_opt="" if i==0 else "same"
			hobs[seq,vst,var].Draw(draw_opt)
			#! Also write to .root file
			hobs[seq,vst,var].Write()
		#! Draw herr
		herr['EF',vst,var].Draw("e hist same")
		herr['EF',vst,var].Write()			

		#! Add TLegend if pad==1
		if pad==1:
			l=ROOT.TLegend(0.60,0.70,0.90,0.90)
			#l=ROOT.TLegend(0.55,0.80,0.90,0.90)
			#l.SetNColumns(3)
			#l.SetFillStyle(0)
			#l.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
			l.SetBorderSize(1)
			l.SetTextSize(0.04)
			for seq in ['EC','EF','ST']:
				lbl=seq
				if seq=="ST": lbl="%s-norm"%(seq)
				l.AddEntry(hobs[seq,vst,var],lbl,"p")#EF
			l.AddEntry(herr['EF',vst,var],"syst-err","l")
			l.Draw()

	#! Save canvas
	outdir_1D_q2bin=("%s/%s"%(outdir,q2bin))
	if not os.path.exists(outdir_1D_q2bin):
		os.makedirs(outdir_1D_q2bin)
	c.SaveAs("%s/c_w%s_q%s.png"%(outdir_1D_q2bin,wbin,q2bin))
	#c.SaveAs("%s/c_w%s_q%s.eps"%(outdir_1D_q2bin,wbin,q2bin))
	#! Convert .png->.eps->.pdf since in .eps and .pdf produced directly by ROOT
	#! the y-axis titles are cutoff
	#! 1. .png->.eps
	cmd=['convert',"%s/c_w%s_q%s.png"%(outdir_1D_q2bin,wbin,q2bin),"%s/c_w%s_q%s.eps"%(outdir_1D_q2bin,wbin,q2bin)]
	tool=subprocess.Popen(cmd,stderr=subprocess.STDOUT)
	tool.wait()
	#! 2. .png->.pdf
	cmd=['convert',"%s/c_w%s_q%s.png"%(outdir_1D_q2bin,wbin,q2bin),"%s/c_w%s_q%s.pdf"%(outdir_1D_q2bin,wbin,q2bin)]
	tool=subprocess.Popen(cmd,stderr=subprocess.STDOUT)
	tool.wait()
	
	#! 3. txt output
	outdir_txt_1D=("%s/txt/q%s/w%s"%(outdir,q2bin,wbin))
	if not os.path.exists(outdir_txt_1D):
		os.makedirs(outdir_txt_1D)

	#! Loop over rslts, per seq, per bin in rslts and write to file
	for seq in ['EC','EF','ST']:
		for item in PAD_MAP['1D']:
			pad,vst,var=item[0],item[1],item[2]
			ftxt=open("%s/%d_%s_%s.txt"%(outdir_txt_1D,vst,var,seq),"w")
			nbins=hobs[seq,vst,var].GetNbinsX()
			for ibin in range(nbins):
				bin=ibin+1

				mu   =hobs[seq,vst,var].GetBinContent(ibin+1)
				sg_mu=hobs[seq,vst,var].GetBinError(ibin+1)
				sg   =herr[seq,vst,var].GetBinContent(ibin+1)
				sg_sg=herr[seq,vst,var].GetBinError(ibin+1)


				#! rel_err,sg_rel_err
				if mu==0 or sg==0: 
					rel_err=0
					sg_rel_err=0
				else:
					rel_err=(sg/mu)*100
					sg_rel_err=(rel_err)*math.sqrt((sg_sg/sg)**2+(sg_mu/mu)**2)

				#! Get wbin information that will be needed to label data
				bin_le=hobs[seq,vst,var].GetBinLowEdge(ibin+1)
				bin_ue=hobs[seq,vst,var].GetBinLowEdge(ibin+1+1)
												
				ftxt.write("%d,%f,%f = (%f +/- %f),(%f +/- %f),(%f +/- %f)\n"%(bin,bin_le,bin_ue,mu,sg_mu,sg,sg_sg,rel_err,sg_rel_err))
			ftxt.close()

def plot_and_write_Obs_R2(q2wbin,hobs,herr,outdir,froot):
	#! Get q2bin, wbin information
	q2bin=get_q2bin(q2wbin)
	wbin=get_wbin(q2wbin)

	#! General plotting aesthetics
	plot_1D_athtcs()

	#! For saving output to .root file, create and cd into q2wbindir
	q2wbin_dir_R2=froot.GetDirectory(q2wbin)
	if q2wbin_dir_R2==None: 
		q2wbin_dir_R2=froot.mkdir(q2wbin)
	q2wbin_dir_R2.cd()

	for R2 in R2L:
		#! 1. Plot all seqs together with ST appropriately normed to EF
		npadsx=3
		npadsy=4
		npads=npadsx*npadsy
		c=ROOT.TCanvas("c","c",1000,1000)
		pad_t=ROOT.TPad("pad_t","Title pad",0.15,0.945,0.85,1.00)
		pad_p=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
		pad_p.Draw()
		pad_t.Draw()
		pad_t.cd()
		pt=ROOT.TPaveText(.05,.1,.95,.8)
		pt.AddText("%s for Q2_W bin=%s"%(R2_NAMED[R2],q2wbin))
		pt.SetTextSize(0.42)
		pt.Draw()
		pad_p.Divide(npadsx,npadsy)
		#! TLine, for each pad, for ln at y=0
  		ln=[0 for i in range(npads)]
		for item in PAD_MAP['R2']:
			pad,vst,var=item[0],item[1],item[2]
			#print "pad,vst,var=",pad,vst,var
			gpad=pad_p.cd(pad)

			if (R2=='D' or R2=='E') and var!= 'ALPHA':continue

			#! First set minimum and maximum of y-axis
			# minl=[hobs[R2,seq,vst,var].GetMinimum() for seq in ['EC','EF','ST']]
			# maxl=[hobs[R2,seq,vst,var].GetMaximum() for seq in ['EC','EF','ST']]
			minl=[hobs[R2,seq,vst,var].GetBinContent(hobs[R2,seq,vst,var].GetMinimumBin()) for seq in ['EC','EF','ST']]
			maxl=[hobs[R2,seq,vst,var].GetBinContent(hobs[R2,seq,vst,var].GetMaximumBin()) for seq in ['EC','EF','ST']]
			minimum=min(minl)
			maximum=max(maxl)
			fctr_not_A=1.5#1.5
			fctr_A=1.1
			for seq in ['EC','EF','ST']:		
				if R2!="A": #! B,C,D can oscillate around 0
					hobs[R2,seq,vst,var].SetMinimum(minimum-math.fabs(fctr_not_A*minimum))
					hobs[R2,seq,vst,var].SetMaximum(maximum+math.fabs(fctr_not_A*maximum))
				elif R2=="A": #! These distributions are same as Obs-1D
					hobs[R2,seq,vst,var].SetMinimum(0)
					hobs[R2,seq,vst,var].SetMaximum(fctr_A*maximum)

			#! Now draw hobs
			for i,seq in enumerate(['ST','EC','EF']):
				draw_opt="" if i==0 else "same"
				hobs[R2,seq,vst,var].Draw(draw_opt)
				#! Also write to .root file
				hobs[R2,seq,vst,var].Write()
			#! Draw herr
			herr[R2,'EF',vst,var].Draw("e hist same")
			herr[R2,'EF',vst,var].Write()

			#! Draw TLine only if ymin<0
			if hobs[R2,'EC',vst,var].GetMinimum()<0:
				ln[pad-1]=ROOT.TLine(hobs[R2,'EC',vst,var].GetXaxis().GetXmin(),0,hobs[R2,'EC',vst,var].GetXaxis().GetXmax(),0)
				ln[pad-1].Draw("same")				

			#! Draw TLegend if pad==3 (since this pad, in the upper right corner, is not shrouded by the TPaveText-Q2W label)
			if pad==1:
				l=ROOT.TLegend(0.60,0.70,0.90,0.90)
				#l=ROOT.TLegend(0.55,0.80,0.90,0.90)
				#l.SetNColumns(3)
				#l.SetFillStyle(0)
				#l.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
				l.SetBorderSize(1)
				l.SetTextSize(0.04)
				#l=ROOT.TLegend(0.7,0.8,0.9,1.0)
				for seq in ['EC','EF','ST']:
					lbl=seq
					if seq=="ST": lbl="%s-norm"%(seq)
					l.AddEntry(hobs[R2,seq,vst,var],lbl,"p")
				l.AddEntry(herr[R2,seq,vst,var],"syst-err","l")
				l.Draw()

		#! Save canvas
		outdir_R2_q2bin=("%s/%s/%s"%(outdir,q2bin,R2))
		if not os.path.exists(outdir_R2_q2bin):
			os.makedirs(outdir_R2_q2bin)
		c.SaveAs("%s/c_w%s_q%s.png"%(outdir_R2_q2bin,wbin,q2bin))
		#c.SaveAs("%s/c_w%s_q%s.eps"%(outdir_R2_q2bin,wbin,q2bin))
		#! Convert .png->.eps->.pdf since in .eps and .pdf produced directly by ROOT
		#! the y-axis titles are cutoff
		#! 1. .png->.eps
		cmd=['convert',"%s/c_w%s_q%s.png"%(outdir_R2_q2bin,wbin,q2bin),"%s/c_w%s_q%s.eps"%(outdir_R2_q2bin,wbin,q2bin)]
		tool=subprocess.Popen(cmd,stderr=subprocess.STDOUT)
		tool.wait()
		#! 2. .png->.pdf
		cmd=['convert',"%s/c_w%s_q%s.png"%(outdir_R2_q2bin,wbin,q2bin),"%s/c_w%s_q%s.pdf"%(outdir_R2_q2bin,wbin,q2bin)]
		tool=subprocess.Popen(cmd,stderr=subprocess.STDOUT)
		
		#! 3. txt output
		outdir_txt_R2=("%s/txt/q%s/w%s/%s"%(outdir,q2bin,wbin,R2))
		if not os.path.exists(outdir_txt_R2):
			os.makedirs(outdir_txt_R2)

		#! Loop over rslts, per seq, per bin in rslts and write to file
		for seq in ['EC','EF','ST']:
			for item in PAD_MAP['R2']:
				pad,vst,var=item[0],item[1],item[2]
				
				if (R2=='D' or R2=='E') and var!= 'ALPHA':continue

				ftxt=open("%s/%d_%s_%s.txt"%(outdir_txt_R2,vst,var,seq),"w")
				nbins=hobs[R2,seq,vst,var].GetNbinsX()
				for ibin in range(nbins):
					bin=ibin+1

					mu   =hobs[R2,seq,vst,var].GetBinContent(ibin+1)
					sg_mu=hobs[R2,seq,vst,var].GetBinError(ibin+1)
					sg   =herr[R2,seq,vst,var].GetBinContent(ibin+1)
					sg_sg=herr[R2,seq,vst,var].GetBinError(ibin+1)


					#! rel_err,sg_rel_err
					if mu==0 or sg==0: 
						rel_err=0
						sg_rel_err=0
					else:
						rel_err=(sg/mu)*100
						sg_rel_err=(rel_err)*math.sqrt((sg_sg/sg)**2+(sg_mu/mu)**2)

					#! Get wbin information that will be needed to label data
					bin_le=hobs[R2,seq,vst,var].GetBinLowEdge(ibin+1)
					bin_ue=hobs[R2,seq,vst,var].GetBinLowEdge(ibin+1+1)
												
					ftxt.write("%d,%f,%f = (%f +/- %f),(%f +/- %f),(%f +/- %f)\n"%(bin,bin_le,bin_ue,mu,sg_mu,sg,sg_sg,rel_err,sg_rel_err))
				ftxt.close()

# def get_norm_factor_ST():
# 	#! Fill norm_factor
# 	for pm in PAD_MAP['1D']:
# 		vst,var=pm[1],pm[2]
# 		itgrl_EF=hobs['EF',vst,var].Integral()
# 		itgrl_ST=hobs['ST',vst,var].Integral()
# 		if itgrl_EF==0 or itgrl_ST==0: nf=1
# 		else:    					   nf=itgrl_EF/itgrl_ST
# 		norm_factor[q2bin,wbin,vst,var]=nf
	
def proc_Obs_itg():
	print "*** Processing Obs_itg from 'cmb_vst_SE' and storing SYST_ERR ***"
	for q2 in ["lowQ2","highQ2"]:
		#! Create outdir=<q2>_thesis_<date>/Obs_itg
		# outdir="%s/%s_thesis_%s/Obs_itg"%(OBSDIR_E162,q2,DATE)
		# if DBG==True:outdir="%s/%s_thesis_%s_dbg/Obs_itg"%(OBSDIR_E162,q2,DATE)
		outdir="%s/%s/Obs_itg"%(OUTDIR,q2)
		#if DBG==True: outdir="%s_dbg/%s/Obs_itg"%(OUTDIR,q2)
		#print outdir
		if not os.path.exists(outdir):
			os.makedirs(outdir)

		#! Setup q2min, q2max as per q2
		q2min,q2max=Q2MIN[q2],Q2MAX[q2]
		print "Processing q2,q2min,q2max=",q2,q2min,q2max

		#! 1. Open appropriate .root files
		f_RSLT=root_open(F_RSLT['itg',q2])
		q2wbinl=get_q2wbinlist(f_RSLT,q2min,q2max,WMIN,WMAX)
	
		#! 2. For every q2bin,
		#! i.   Get from f_RSLT, hobs[seq],herr[seq] where seq={EC,EF}
		#! ii.  Adjust sg_sg in herr[seq]
		#! iii. Plot and write to file(.txt and .root), for both seq, hobs[seq],herr[seq]
		#! iv.  Update SYST_ERR

		#! Create .root file
		froot=ROOT.TFile("%s/Obs_itg.root"%outdir,"RECREATE")
		#! Create .txt file
		#! Create inside loop for every q2wbin
		for q2wbin in q2wbinl:
			# #! Get q2bin information
			# q2bin=get_q2bin(q2wbin)

			#! i. Get from f_RSLT, hobs[seq],herr[seq] where seq={EC,EF}
			hobs,herr=OrderedDict(),OrderedDict()
			for seq in ['EC','EF']:
				hobs[seq]=f_RSLT.Get("%s/hW_obs_%s_THETA"%(q2wbin,seq))
				herr[seq]=f_RSLT.Get("%s/hW_err_%s_THETA"%(q2wbin,seq))

			#! ii.  Update herr[seq] after adjust sg_sg in herr[seq]
			#! + Current adjustment is to set sg_sg=0
			for seq in ['EC','EF']:
				nbins=herr['EF'].GetNbinsX()
				for ibin in range(nbins):
					herr[seq].SetBinError(ibin+1,0)
				
			#! iii. Plot and write to file, for both seq, hobs[seq],herr[seq]
			#! + Note that within this function 'iv.  Update SYST_ERR' is also done
			plot_and_write_Obs_itg(q2wbin,hobs,herr,outdir,froot)

def get_norm_factor_ST():
	print "*** Processing Obs_1D to get norm_factor from SSBands/cutsncors1 ***"
	for q2 in ["lowQ2","highQ2"]:
		#! Create outdir=<q2>_thesis_<date>/Obs_<obs>
		# outdir="%s/%s_thesis_%s/Obs_1D"%(OBSDIR_E162,q2,DATE)
		# if DBG==True:outdir="%s/%s_thesis_%s_dbg/Obs_1D"%(OBSDIR_E162,q2,DATE)
		outdir="%s/%s/Obs_1D"%(OUTDIR,q2)
		#if DBG==True: outdir="%s_dbg/%s/Obs_1D"%(OUTDIR,q2)
		#print outdir
		if not os.path.exists(outdir):
			os.makedirs(outdir)

		#! Setup q2min, q2max as per q2
		q2min,q2max=Q2MIN[q2],Q2MAX[q2]
		print "Processing q2,q2min,q2max=",q2,q2min,q2max

		#! 1. Open appropriate .root files
		f_RSLT=root_open(F_RSLT['1D',q2])

		if DBG==True: q2wbinl=get_q2wbinlist(f_RSLT,q2min,q2max,WMIN,WMAX,dbg=True,dbg_bins=3)
		else:         q2wbinl=get_q2wbinlist(f_RSLT,q2min,q2max,WMIN,WMAX)
	
		#! 2. For every q2bin
			#! i.   Get from f_RSLT, for every seq in {EC,SF,ST}, hobs[seq,vst,var] as per respective PAD_MAP[obs]=(pad,vst,var)
			#! ii. Get norm_factor(q2,w,vst,var)
		#! 3. Plot and study norm_factor(q2,w,vst,var)
				
		#! Create structure to keep norm_factor
		#! +norm_factor[q2bin,wbin,vst,var]
		norm_factor=OrderedDict()
		for q2wbin in q2wbinl:
			#! Get q2bin, wbin information
			q2bin=get_q2bin(q2wbin)
			wbin=get_wbin(q2wbin)

			#! i. Get from f_RSLT, hobs[seq,vst,var]
			hobs=OrderedDict()
			for seq in ['EC','EF','ST']:
				for pm in PAD_MAP['1D']:
					vst,var=pm[1],pm[2]
					#! for 1D, SF was used instead of ST (SF=ST because of acceptance corr. and hole-filling)
					if seq=='ST': hobs[seq,vst,var]=f_RSLT.Get("%s/h1_%s_%d_%s"%(q2wbin,'SF',vst,var))
					else:		  hobs[seq,vst,var]=f_RSLT.Get("%s/h1_%s_%d_%s"%(q2wbin,seq,vst,var))
				
					
			#! ii. Fill norm_factor(q2,w,vst,var)
			for pm in PAD_MAP['1D']:
				vst,var=pm[1],pm[2]
				itgrl_EF=hobs['EF',vst,var].Integral()
				itgrl_ST=hobs['ST',vst,var].Integral()
				if itgrl_EF==0 or itgrl_ST==0: nf=1
				else:    					   nf=itgrl_EF/itgrl_ST
				#! remove "bad bins"
				#if q2=="lowQ2" and (wbin=="1.400-1.425" or wbin=="1.425-1.450" or wbin=="1.450-1.475" or wbin=="1.475-1.500"): continue
				# if wbin=="1.400-1.425": continue
				# if q2=="lowQ2" and (wbin=="1.425-1.450" or wbin=="1.450-1.475"): continue
				#!fill
				norm_factor[q2bin,wbin,vst,var]=nf
			#! end q2wbinl loop
		#! 3. Plot and study norm_factor(q2,w,vst,var)
		#fig=plt.figure()
		#plt.hist(norm_factor.values())
		#plt.title(r"nf(q2,w,vst,var) \n  \mu,sg,rel_err(%%)=%.5f,%.5f,%.5f"%(mu,sg,rel_err))
		fig,ax=plt.subplots(figsize=(12,8))
		ax.hist(norm_factor.values())
		mu=np.mean(norm_factor.values())
		sg=np.std(norm_factor.values())
		rel_err=(sg/mu)*100
		fig.suptitle("sim=%s"%q2,fontsize='xx-large')
		fig.subplots_adjust(top=0.85)
		ax.set_title(r'nf(q2bin,wbin,vst,var)' "\n" r'$\mu=%.5f,\sigma=%.5f,\mu/\sigma=%.2f$(%%)'%(mu,sg,rel_err),fontsize='xx-large')
		fig.savefig("%s/nf_%s.png"%(outdir,q2))
		if q2=="lowQ2":
			for k in norm_factor:
				if norm_factor[k]>0.006: print k			
		#! Save norm_factor in NORM_FACTOR[q2]
		NORM_FACTOR[q2]=mu
#! end q2 loop
print "*** Done: Processing Obs_%s to get norm_factor from SSBands/cutsncors1 ***"

def norm_ST_shape_and_amplitude(hobs,q2):
	#! first determine obs from number of keys in hobs
	if   len(hobs.keys()[0])==3: obs='1D'
	elif len(hobs.keys()[0])==4: obs='R2'
	else: sys.exit("norm_ST: Could not determine obs. Exiting.")

	for pm in PAD_MAP[obs]:
		vst,var=pm[1],pm[2]
		if obs=='1D':
			hobs['ST',vst,var].Sumw2()
			hobs['ST',vst,var].Scale(NORM_FACTOR[q2])
		elif obs=='R2': 
			for R2 in R2L:
				if (R2=='D' or R2=='E') and var!= 'ALPHA':continue
				hobs[R2,'ST',vst,var].Sumw2()
				hobs[R2,'ST',vst,var].Scale(NORM_FACTOR[q2])

def norm_ST_shape(hobs):
	#! first determine obs from number of keys in hobs
	if   len(hobs.keys()[0])==3: obs='1D'
	elif len(hobs.keys()[0])==4: obs='R2'
	else: sys.exit("norm_ST: Could not determine obs. Exiting.")

	for pm in PAD_MAP[obs]:
		vst,var=pm[1],pm[2]
		if obs=='1D':
			#! First call Sumw2()
			hobs['EF',vst,var].Sumw2()
			hobs['ST',vst,var].Sumw2()
			#! Now obtain the norm factor
			intgrl_EF=hobs['EF',vst,var].Integral()
			intgrl_ST=hobs['ST',vst,var].Integral()
			#! 1. Get scale factor
			if intgrl_EF==0 or intgrl_ST==0: scale_factor=1
			else:    						 scale_factor=intgrl_EF/intgrl_ST
			#! 2. Scale bin contents of ST
			hobs['ST',vst,var].Scale(scale_factor)
		elif obs=='R2': 
			for R2 in R2L:
				if (R2=='D' or R2=='E') and var!= 'ALPHA':continue
				#! First call Sumw2()
				hobs[R2,'EF',vst,var].Sumw2()
				hobs[R2,'ST',vst,var].Sumw2()
				if R2=='A': #! then do intgrl norm like for 1D
					#! Now obtain the norm factor
					intgrl_EF=hobs[R2,'EF',vst,var].Integral()
					intgrl_ST=hobs[R2,'ST',vst,var].Integral()
					#! 1. Get scale factor
					if intgrl_EF==0 or intgrl_ST==0: scale_factor=1
					else:    						 scale_factor=intgrl_EF/intgrl_ST
					#! 2. Scale bin contents of ST
					hobs[R2,'ST',vst,var].Scale(scale_factor)
				else: #! Do ampltd normalization
					#! get EF ampltd
					max_EF=hobs[R2,'EF',vst,var].GetBinContent(hobs[R2,'EF',vst,var].GetMaximumBin())
					min_EF=hobs[R2,'EF',vst,var].GetBinContent(hobs[R2,'EF',vst,var].GetMinimumBin())#GetMinimum()
					ampltd_EF=max([max_EF,math.fabs(min_EF)])
					#! get ST ampltd
					max_ST=hobs[R2,'ST',vst,var].GetBinContent(hobs[R2,'ST',vst,var].GetMaximumBin())
					min_ST=hobs[R2,'ST',vst,var].GetBinContent(hobs[R2,'ST',vst,var].GetMinimumBin())
					ampltd_ST=max([max_ST,math.fabs(min_ST)])
					#! Now get scale factor and scale bin contents
					if ampltd_EF==0 or ampltd_ST==0: scale_factor=1
					else:    						 scale_factor=ampltd_EF/ampltd_ST
					hobs[R2,'ST',vst,var].Scale(scale_factor)

def proc_Obs_1D_and_R2():
	for obs in ['1D','R2']:
		print "*** Processing Obs_%s from SSBands/cutsncors1 ***"%obs
		for q2 in ["lowQ2","highQ2"]:
			#! Create outdir=<q2>_thesis_<date>/Obs_<obs>
			# outdir="%s/%s_thesis_%s/Obs_%s"%(OBSDIR_E162,q2,DATE,obs)
			# if DBG==True: outdir="%s/%s_thesis_%s_dbg/Obs_%s"%(OBSDIR_E162,q2,DATE,obs)
			outdir="%s/%s/Obs_%s"%(OUTDIR,q2,obs)
			#if DBG==True: outdir="%s_dbg/%s/Obs_%s"%(OUTDIR,q2,obs)
			#print outdir
			if not os.path.exists(outdir):
				os.makedirs(outdir)

			#! Setup q2min, q2max as per q2
			q2min,q2max=Q2MIN[q2],Q2MAX[q2]
			print "Processing q2,q2min,q2max=",q2,q2min,q2max

			#! 1. Open appropriate .root files
			f_RSLT=root_open(F_RSLT[obs,q2])

			if DBG==True: q2wbinl=get_q2wbinlist(f_RSLT,q2min,q2max,WMIN,WMAX,dbg=True,dbg_bins=3)
			else:         q2wbinl=get_q2wbinlist(f_RSLT,q2min,q2max,WMIN,WMAX)
			if DBG==True:
				print "q2wbinl=",q2wbinl
		
			#! 2. For every qw2bin
			#! i.   Get from f_RSLT, for every seq in {EC,SF,ST}, hobs[seq,vst,var] as per respective PAD_MAP[obs]=(pad,vst,var)
			#! ii.  Create herr[seq,vst,var] = SYST_ERR[q2bin,wbin]*hobs[seq,vst,var] (SYST_ERR contains relative error)
			#! iii. Plot and write to file(.txt and .root), for all seq, hobs[seq],herr[seq]
		
			#! Create .root file
			froot=ROOT.TFile("%s/Obs_%s.root"%(outdir,obs),"RECREATE")

			#! Create .txt file
			#! Create inside loop for every q2wbin

			for q2wbin in q2wbinl:
				#! Get q2bin, wbin information
				q2bin=get_q2bin(q2wbin)
				wbin=get_wbin(q2wbin)

				#! i. Get from f_RSLT, hobs[seq,vst,var]
				hobs=OrderedDict()
				for seq in ['EC','EF','ST']:
					for pm in PAD_MAP[obs]:
						vst,var=pm[1],pm[2]
						if obs=='1D': #! for 1D, SF was used instead of ST (SF=ST because of acceptance corr. and hole-filling)
							if seq=='ST': hobs[seq,vst,var]=f_RSLT.Get("%s/h1_%s_%d_%s"%(q2wbin,'SF',vst,var))
							else:		  hobs[seq,vst,var]=f_RSLT.Get("%s/h1_%s_%d_%s"%(q2wbin,seq,vst,var))
						#! Set aesthetics
							hobs[seq,vst,var].SetMarkerStyle(MRKS_SEQ[seq])
							hobs[seq,vst,var].SetMarkerColor(CLRS_SEQ[seq])
							hobs[seq,vst,var].SetLineColor(CLRS_SEQ[seq])
						elif obs=='R2':
							for R2 in R2L:
								if (R2=='D' or R2=='E') and var != 'ALPHA': continue
								hobs[R2,seq,vst,var]=f_RSLT.Get("%s/hR2_%s_%s_%d_%s"%(q2wbin,R2,seq,vst,var))
								#! Set aesthetics
								hobs[R2,seq,vst,var].SetMarkerStyle(MRKS_SEQ[seq])
								hobs[R2,seq,vst,var].SetMarkerColor(CLRS_SEQ[seq])
								hobs[R2,seq,vst,var].SetLineColor(CLRS_SEQ[seq])	
				if DBG==True:
					print "hobs:"
					print hobs		

				#! ii.  Create herr[seq,vst,var]:
				#! + if seq='EF' or seq=='EC': err=SYST_ERR[q2bin,wbin]/100)*hobs['EF',vst,var]
				#! + if seq='ST':              err=0
				herr=OrderedDict()
				for k in hobs:
					#! Obtain seq 
					if obs=='1D':
						seq,vst,var=k[0],k[1],k[2]
					elif obs=='R2':
						R2,seq,vst,var=k[0],k[1],k[2],k[3]
						if (R2=='D' or R2=='E') and var!= 'ALPHA':continue
									
					herr[k]=hobs[k].Clone("herr_%s"%hobs[k].GetName())
					#! Set bin contents
					nbins=herr[k].GetNbinsX()
					for ibin in range(nbins):
						#! obtain err as per seq
						if seq=='EF' or seq=='EC':
							binc=hobs[k].GetBinContent(ibin+1)
							#print "SYST_ERR[q2bin,wbin],binc",SYST_ERR[q2bin,wbin],binc
							err=(SYST_ERR[q2bin,wbin]/100)*binc
							err_err=0
						elif seq=='ST':
							err=0
							err_err=0
						#! Now set in histogram
						herr[k].SetBinContent(ibin+1,err)
						herr[k].SetBinError(ibin+1,err_err)
					#! Adjust y-axis title for herr
					herr[k].SetYTitle("Error [%%] in %s"%(hobs[k].GetYaxis().GetTitle()))
					#! Set general aesthetics
					herr[k].SetFillColor(ROOT.gROOT.ProcessLine("kBlack"))
					herr[k].SetFillStyle(3003)
					herr[k].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
					herr[k].SetMarkerStyle(ROOT.gROOT.ProcessLine("kDot"))
					herr[k].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
				if DBG==True:
					print "herr:"
					print herr

				#! iii. Normalize ST using NORM_FACTOR[q2]
				if   NORM_ST=="shape_and_amplitude": norm_ST_shape_and_amplitude(hobs,q2)
				elif NORM_ST=="shape":				 norm_ST_shape(hobs)
				#! iv. Plot and write to file(.txt and .root), for all seq, hobs[seq],herr[seq]
				if obs=='1D':
					plot_and_write_Obs_1D(q2wbin,hobs,herr,outdir,froot)
				elif obs=='R2': #! TEST TMP
					#print "not done"
					plot_and_write_Obs_R2(q2wbin,hobs,herr,outdir,froot)
			#! end q2wbinl loop
		#! end q2 loop
	print "*** Done: Processing Obs_%s from SSBands/cutsncors1 ***"%obs
#! end OBS loop

#! main
if NORM_ST=="shape_and_amplitude" and STUDY_NORM_FACTOR_ONLY: #! then only do the indepenent step 2.: study_norm_factor
	print" Begin: Step 2. \n"
	#! step 2. Get norm_factor using Obs_1D from SSBands/cutsncors1
	#! Create NORM_FACTOR(q2) to be obtained for each simulation
	NORM_FACTOR=OrderedDict()
	get_norm_factor_ST()
	print "NORM_FACTOR=",NORM_FACTOR
	print" Done: Step 2. \n"
	sys.exit()

print" Begin: Step 1. \n"
#! 1. First process Obs_itg from 'cmb_vst_SE'
#! + Before plotting and writing results update sg_sg to 0 in 'herr'
#!   (In a discussion with RG on [10-25-16], it was decided to "simply" set the error on systematic err=0, 
#!    for now atleast)
#! + Create and store SYST_ERR[q2bin,wbin] that will be used later for Obs_1D and Obs_R2

#! Create dictionary to store SYST_ERR[q2bin,wbin]
SYST_ERR=OrderedDict()
proc_Obs_itg()
print" Done: Step 1. \n"

if NORM_ST=="shape_and_amplitude":
	print" Begin: Step 2. \n"
	#! step 2. Get norm_factor using Obs_1D from SSBands/cutsncors1
	#! Create NORM_FACTOR(q2) to be obtained for each simulation
	NORM_FACTOR=OrderedDict()
	get_norm_factor_ST()
	print "NORM_FACTOR=",NORM_FACTOR
	print" Done: Step 2. \n"

print" Begin: Step 3. \n"
#! step 3. Now process Obs_1D and Obs_R2 from SSBands/cutsncors1
#! + Apply SYST_ERR(q2bin,wbin) obtained from Obs_itg
#! + Along with EC,EH, use ST in plotting -- use NORM_FACTOR(q2) from step 2. --  and writing results
proc_Obs_1D_and_R2()
print" Done: Step 3. \n"

	
	



		

