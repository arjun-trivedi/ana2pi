#!/usr/bin/python
from __future__ import division
import sys,os
from collections import OrderedDict

from rootpy.io import root_open, DoesNotExist
import itertools

import ROOT

import math

from proc_h8 import VARS
from disp_obs import VAR_NAMES_PLAIN, R2S, R2_NAMED

import numpy as np

'''
+ Idea is to use this script to plot the difference between two sets of observables extracted for the ananote
+ These observables are assumed to be in their final form

>plot_diff_obs.py obsdir1 obsdir2 expctd_diff expctd_diff_tlrnc dbg[=False]
+ obsdir1= previous result
+ obsdir2= new result
+ expctd_diff=expected level of relative difference (%)
+ expctd_diff_tlrnc=tolerance level for expctd_diff (%)
+ Difference(%) is plotted relative to obs1 i.e. diff=(obs2-obs1)/obs1
+ Plots are put in obsdir2
+ Assumptions:
	+ path of obsdir1 and obsdir2 is the same
'''
USAGE=">plot_diff_obs.py obsdir1 obsdir2 expctd_diff expctd_diff_tlrnc dbg[=False]"

NOBS=2
OBS1,OBS2=range(NOBS)
OBSDIR=[0 for i in range(NOBS)]

#! Get user args
#! First make sure all args are present
if len(sys.argv)<5:
	sys.exit("Not all arguments provided. Note usage=%s"%USAGE)
OBSDIR[OBS1]=sys.argv[1]
OBSDIR[OBS2]=sys.argv[2]
EXPCTD_DIFF=float(sys.argv[3])
EXPCTD_DIFF_TLRNC=float(sys.argv[4])
print "OBSDIR1=",OBSDIR[OBS1]
print "OBSDIR2=",OBSDIR[OBS2]
print "EXPCTD_DIFF=",EXPCTD_DIFF
print "EXPCTD_DIFF_TLRNC=",EXPCTD_DIFF_TLRNC

if len(sys.argv)>5: #! User input dbg, else use dbg=False
	if sys.argv[5]=="False":
		DBG=False
	elif sys.argv[5]=="True":
		DBG=True
	else:
		print "Please enter dbg=True or False!"
		sys.exit()
else:
	DBG=False

#! Create OUTDIR
obsname1=os.path.basename(os.path.normpath(OBSDIR[OBS1]))
obsname2=os.path.basename(os.path.normpath(OBSDIR[OBS2]))
#print "obsname1=",obsname1
#print "obsname2=",obsname2
obsnamediff='%s_%s'%(obsname1,obsname2)
#! use same path for OUTDIR as OBS2
path=OBSDIR[OBS2].replace(obsname2,'')
path=os.path.normpath(path)
#print "path=",path
OUTDIR="%s/diff_obs/%s"%(path,obsnamediff)
if DBG:
	OUTDIR="%s/diff_obs/%s_dbg"%(path,obsnamediff)
if not os.path.exists(OUTDIR):
	os.makedirs(OUTDIR)
print "OUTDIR=",OUTDIR
#sys.exit()

#! Create FIN(obsn,q2,obst)
#! + obsn=OBS1,OBS2
#! + q2=lowQ2,highQ2
#! + obst=1D,R2,itg
FIN=OrderedDict()
for item in itertools.product([OBS1,OBS2],['lowQ2','highQ2'],['Obs_1D','Obs_R2','Obs_itg']):
	obsn,q2,obst=item[0],item[1],item[2]
	#FIN[(obsn,q2,obst)]=root_open("%s/%s/%s/%s.root"%(OBSDIR[obsn],q2,obst,obst))
	FIN[(obsn,q2,obst)]="%s/%s/%s/%s.root"%(OBSDIR[obsn],q2,obst,obst)
print "Number of files in FIN=(should be 12)",len(FIN)
print "*** FIN pretty print: ***"
for k in FIN:
	obsn,q2,obst=k[0],k[1],k[2]
	print "FIN(%d,%s,%s)=%s"%(obsn,q2,obst,FIN[(obsn,q2,obst)])
print "******"
#sys.exit()

#! Define some objects that will be needed later
#! PAD_MAP[obs]=(pad,vst,var)
PAD_MAP=OrderedDict()
PAD_MAP['Obs_1D']=[(1,1,"M1"),   (2,3,'M2'),   (3,2,'M2'),
				   (4,1,"THETA"),(5,3,'THETA'),(6,2,'THETA'),
				   (7,1,"ALPHA"),(8,3,'ALPHA'),(9,2,'ALPHA')]

PAD_MAP['Obs_R2']=[(1,1,"M1"),    (2,3,'M1'),    (3,2,'M1'), 
				   (4,1,"M2"),    (5,3,'M2'),    (6,2,'M2'),
				   (7,1,"THETA"), (8,3,'THETA'), (9,2,'THETA'),
				   (10,1,"ALPHA"),(11,3,'ALPHA'),(12,2,'ALPHA')]

#! Define binning for hdiff hists that will be created
NBINS=200
XMIN,XMAX=-100,100

#! Define some additional level of diff(%) beyond which the user should be notified
DIFF_THRSHLD_VALS=[10,20,40,60,80,100]

#! CLRS and MRKS
CLRS_SEQ={('EC'):ROOT.gROOT.ProcessLine("kCyan"),
		  ('EF'):ROOT.gROOT.ProcessLine("kBlue")}
MRKS_OBS={(OBS1):ROOT.gROOT.ProcessLine("kFullTriangleDown"),
		  (OBS2):ROOT.gROOT.ProcessLine("kFullTriangleUp")}

CLRS_OBS={(OBS1):ROOT.gROOT.ProcessLine("kBlack"),
		  (OBS2):ROOT.gROOT.ProcessLine("kBlue")}
MRKS_SEQ={('EC'):ROOT.gROOT.ProcessLine("kOpenCircle"),
		  ('EF'):ROOT.gROOT.ProcessLine("kFullCircle")}
		  

#! setup functions that will be needed
def get_q2wbinlist(f,q2min=0.00,q2max=6.00,wmin=0.000,wmax=3.000,dbg=False,dbg_bins=3):
	"""
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
def get_q2bin_le(q2wbin):
	return float(get_q2bin(q2wbin).split('-')[0])
def get_wbin_le(q2wbin):
	return float(get_wbin(q2wbin).split('-')[0])
def get_q2bin_ue(q2wbin):
	return float(get_q2bin(q2wbin).split('-')[1])
def get_wbin_ue(q2wbin):
	return float(get_wbin(q2wbin).split('-')[1])

def calc_rel_diff(obs1,obs2):
	'''
	+ diff(%) is calculated relative to obs1: ((obs2-obs1)/abs(obs1))*100
	+ Note that the sign for the relative difference falls out in the numerator
	+ special cases are noted
	'''
	if np.isclose(obs1,0):
		if np.isclose(obs2,0): #! then diff =0
			diff=0
		else: #! if obs2!=0, then diff has to be calculated relative to obs2, which is always 100%
			diff=((obs2-obs1)/math.fabs(obs2))*100 #! should always work out to 100%
	else:
		diff=((obs2-obs1)/math.fabs(obs1))*100
	return diff

def get_abs_hist(h):
	'''
	Returns habs:
	+ habs_bin_i=abs(h_bin_i)
	'''
	habs=h.Clone()
	nbins=h.GetNbinsX()
	for ibin in range(nbins):
		binc=h.GetBinContent(ibin+1)
		bine=h.GetBinError(ibin+1)
		habs.SetBinContent(ibin+1,math.fabs(binc))
		habs.SetBinError(ibin+1,bine)
	return habs

def get_hlindiff(hobs1,hobs2):
	'''
	+ Histogram version of calc_rel_diff()
	+ returned hlindiff has error bars set to 0 because statistical errors when comparing 
	  difference between hobs1 and hobs2 are accounted for in another part of the code. See line:
		'abs_diff_gt_staterr= math.fabs(abs_diff)>3*stat_err1 and math.fabs(abs_diff)>3*stat_err2'
		 + where stat_err1 and stat_err2 are stat errs, for every bin, for hobs1 and hobs2, respectively

	'''
	#! obs2-obs1
	hlindiff=hobs2.Clone()
	hlindiff.Add(hobs1,-1)
	#! /abs(obs1)
	hobs1_abs=get_abs_hist(hobs1)
	hlindiff.Divide(hobs1_abs)
	#! *100
	hlindiff.Scale(100)
	#! Set error bars to zero
	nbins=hlindiff.GetNbinsX()
	for ibin in range(nbins):
		hlindiff.SetBinError(ibin+1,0)
	#! return
	return hlindiff

def return_unique_2D_list(l):
	'''
	+ l = 2D list
	+ Algorithm:
		+ Go through list elements
		+ If encountering an element the first time, fill it in unique list
		+ If encountering the same element again, ignore
	'''
	#! First set, for all elements, filled(element)=false
	filled=OrderedDict()
	for x in l:
		v1,v2=x[0],x[1]
		filled[(v1,v2)]=False
	#! Now make new unique list
	#! + Go through list
	#! + if filled(element)=false, fill it and set filled(element)=true
	l2=[]
	for x in l:
		v1,v2=x[0],x[1]
		if filled [(v1,v2)]==False:
			l2.append([v1,v2])
			filled [(v1,v2)]=True
		else:
			continue
	return l2

def proc_Obs_1D_and_R2():
	print "*** In proc_Obs_1D_and_R2() ***\n"

	#! create structure to contain:
	#! + hobs[obsn,obst]=f(q2wb,seq,vst,var)
	#! + hdiff[obst]    =f(q2wb,seq,vst,var)
	#! hobs
	hobs=OrderedDict()
	for obsn in range(NOBS): 
		hobs[(obsn,'Obs_1D')]=OrderedDict()
		hobs[(obsn,'Obs_R2')]=OrderedDict()
	#! hdiff
	#! + Note two sets of hists: one for all entries and only for those where diff>3*stat_err
	hdiff,hdiff_gt_staterr=OrderedDict(),OrderedDict()
	for obst in ['Obs_1D','Obs_R2']:
		hdiff[obst]=OrderedDict()
		hdiff_gt_staterr[obst]=OrderedDict()
	#! hlindiff
	#! + only for those where diff>3*stat_err
	hlindiff=OrderedDict()
	for obst in ['Obs_1D','Obs_R2']:
		hlindiff[obst]=OrderedDict()


	#! create file to store, in text form, any difference found to be greater than a threshold
	fout=OrderedDict()
	#! 1D
	OUTDIR_1D="%s/Obs_1D"%(OUTDIR)
	if not os.path.exists(OUTDIR_1D):
		os.makedirs(OUTDIR_1D)
	fout['Obs_1D']=open("%s/1D.txt"%(OUTDIR_1D),'w')
	#! R2
	OUTDIR_R2="%s/Obs_R2"%(OUTDIR)
	if not os.path.exists(OUTDIR_R2):
		os.makedirs(OUTDIR_R2)
	fout['Obs_R2']=open("%s/R2.txt"%(OUTDIR_R2),'w')
	#! Descriptive opening remark in fout
	for obst in ['Obs_1D','Obs_R2']:
		#fout[obst].write("*** This file will note bins where diff in %s is greater than %.2f%% ***\n"%(obst,EXPCTD_DIFF))
		fout[obst].write("*** This file will note bins where absolute diff in %s is greater than stat errs ***\n"%(obst))
		fout[obst].write("\n")

	#! Now begin to fill this structure
	for item in itertools.product(['Obs_1D','Obs_R2'],['lowQ2','highQ2']):
		obst,q2=item[0],item[1]
	
		#! open f[obsn]
		f=[0 for i in range(NOBS)]
		f[OBS1]=root_open(FIN[(OBS1,q2,obst)])
		f[OBS2]=root_open(FIN[(OBS2,q2,obst)])

		#! get q2wbinl from any of the files
		if DBG==True: q2wbinl=get_q2wbinlist(f[OBS1],dbg=True,dbg_bins=3)
		else:         q2wbinl=get_q2wbinlist(f[OBS1])
		if DBG==True:
			print "q2wbinl=",q2wbinl
		
		#! 2. For every qw2bin
		#! i. Get from f[obsn], for every seq in {EC,EF}, hobs[obsn,obst]=h[seq,vst,var] as per respective PAD_MAP[obs]=(pad,vst,var)
					
		for iq2wbin,q2wbin in enumerate(q2wbinl):
			#! Get q2bin, wbin information
			q2bin=get_q2bin(q2wbin)
			wbin=get_wbin(q2wbin)
			
			#! i. Get from f[obsn], hobs[obsn,obst][q2wb,seq,vst,var]
			for iseq,seq in enumerate(['EC','EF']):
				for pm in PAD_MAP[obst]:
					vst,var=pm[1],pm[2]
					if obst=='Obs_1D': #! for 1D, SF was used instead of ST (SF=ST because of acceptance corr. and hole-filling)
						hobs[OBS1,obst][q2wbin,seq,vst,var]=f[OBS1].Get("%s/h1_%s_%d_%s"%(q2wbin,seq,vst,var))
						hobs[OBS2,obst][q2wbin,seq,vst,var]=f[OBS2].Get("%s/h1_%s_%d_%s"%(q2wbin,seq,vst,var))
						#! Create and fill hdiff
						hname ="hdiff_%d_1D_%s_VST%d_%s"%(iq2wbin+1,seq,vst,var)
						htitle="hdiff_%s_1D_%s_VST%d_%s"%(q2wbin,seq,vst,var)
						hdiff[obst][q2wbin,seq,vst,var]=ROOT.TH1F(hname,htitle,NBINS,XMIN,XMAX)
						hdiff_gt_staterr[obst][q2wbin,seq,vst,var]=ROOT.TH1F("%s_gt_staterr_"%hname,"%s_gt_staterr"%htitle,NBINS,XMIN,XMAX)
						#! loop over bins of hobs1 and fill ((hobs1-hobs2)/hobs1)*100
						nbins=hobs[OBS1,obst][q2wbin,seq,vst,var].GetNbinsX()
						#! + Initialize hlindiff_drawn=False.
						#! + This will be used later for cases where hlindiff=hobs1-hobs2 needs to be drawn
						#! + This is required because diffs are checked bin by bin, and if there is a diff the whole hist needs to be drawn
						#! + However, if there more bins are different, the hist does not need to be re-drawn
						hlindiff_drawn=False
						for ibin in range(nbins):
							obs1=hobs[OBS1,obst][q2wbin,seq,vst,var].GetBinContent(ibin+1)
							obs2=hobs[OBS2,obst][q2wbin,seq,vst,var].GetBinContent(ibin+1)
							stat_err1=hobs[OBS1,obst][q2wbin,seq,vst,var].GetBinError(ibin+1)
							stat_err2=hobs[OBS2,obst][q2wbin,seq,vst,var].GetBinError(ibin+1)
							diff=calc_rel_diff(obs1,obs2)
							abs_diff=math.fabs(obs2-obs1)
							#! Fill hdiff
							hdiff[obst][q2wbin,seq,vst,var].Fill(diff)
							#! Check if abs_diff is greater than stat errs
							abs_diff_gt_staterr= math.fabs(abs_diff)>3*stat_err1 and math.fabs(abs_diff)>3*stat_err2
							#! If abs_diff is greater than stat errs:
							#!    1. get hlindiff (#! hlindiff is basically histogram version of calc_rel_diff())
							#!    2. Fill hdiff_gt_staterr 
							#!    3. Write diff info to file
							if abs_diff_gt_staterr:
								#! 1. get hlindiff,
								if hlindiff_drawn==False:
									hlindiff[obst][q2wbin,seq,vst,var]=get_hlindiff(hobs[OBS1,obst][q2wbin,seq,vst,var],hobs[OBS2,obst][q2wbin,seq,vst,var])
									# if var=='M1' or var=='M2':
									# 	wmax=get_wbin_ue(q2wbin)
									# 	setM1M2axisrange(hlindiff[obst][q2wbin,seq,vst,var],vst,var,wmax)
									hlindiff_drawn=True
								#! 2. Fill hdiff_gt_staterr
								hdiff_gt_staterr[obst][q2wbin,seq,vst,var].Fill(diff)
								#! Write info to file
								#! note binn and binle
								binn=ibin+1
								binle=hobs[OBS1,obst][q2wbin,seq,vst,var].GetBinLowEdge(ibin+1)
								#! 3. Now write to to file
								# fout[obst].write("diff=%.2f for q2wb,seq,vst,var,var_name,binn,binle=%s,%s,%d,%s,%s,%d,%.3f\n"
								# 	            %(diff,q2wbin,seq,vst,var,VAR_NAMES_PLAIN[(vst,var)],binn,binle))
								fout[obst].write("*** diff(%), q2wb, seq, vst, var, var_name, binn, binle, obs1+/3*stat_err1, obs2+/-3*stat_err2:***\n")
								fout[obst].write("%.2f, %s, %s, %d, %s, %s, %d, %.3f, %f+/-%f, %f+/-%f\n"%(diff,q2wbin,seq,vst,var,VAR_NAMES_PLAIN[(vst,var)],
												 binn,binle,obs1,3*stat_err1,obs2,3*stat_err2))
								#! check if diff greater than values in DIFF_THRSHLD_VALS
								if any([math.fabs(diff)>x for x in DIFF_THRSHLD_VALS]):
									thrshld_val=get_closest_diff_thrshld_val(diff)
									fout[obst].write("diff exceeded-threshold-%d\n"%thrshld_val)
								fout[obst].write("\n")
					elif obst=='Obs_R2':
						for R2 in R2S:
							if (R2=='D' or R2=='E') and var != 'ALPHA': continue
							hobs[OBS1,obst][q2wbin,R2,seq,vst,var]=f[OBS1].Get("%s/hR2_%s_%s_%d_%s"%(q2wbin,R2,seq,vst,var))
							hobs[OBS2,obst][q2wbin,R2,seq,vst,var]=f[OBS2].Get("%s/hR2_%s_%s_%d_%s"%(q2wbin,R2,seq,vst,var))
							#! Create and fill hdiff
							hname ="hdiff_%d_%s_%s_VST%d_%s"%(iq2wbin+1,R2,seq,vst,var)
							htitle="hdiff_%s_%s_%s_VST%d_%s"%(q2wbin,R2,seq,vst,var)
							hdiff[obst][q2wbin,R2,seq,vst,var]=ROOT.TH1F(hname,htitle,NBINS,XMIN,XMAX)
							hdiff_gt_staterr[obst][q2wbin,R2,seq,vst,var]=ROOT.TH1F("%s_gt_staterr"%hname,"%s_gt_staterr"%htitle,NBINS,XMIN,XMAX)
							#! loop over bins of hobs1 and fill ((hobs1-hobs2)/hobs1)*100
							nbins=hobs[OBS1,obst][q2wbin,R2,seq,vst,var].GetNbinsX()
							#!  Initialize hlindiff_drawn=False.
							#! + This will be used later for cases where hlindiff=hobs1-hobs2 needs to be drawn
							#! + This is required because diffs are checked bin by bin, and if there is a diff the whole hist needs to be drawn
							#! + However, if there more bins are different, the hist does not need to be re-drawn
							hlindiff_drawn=False
							for ibin in range(nbins):
								obs1=hobs[OBS1,obst][q2wbin,R2,seq,vst,var].GetBinContent(ibin+1)
								obs2=hobs[OBS2,obst][q2wbin,R2,seq,vst,var].GetBinContent(ibin+1)
								stat_err1=hobs[OBS1,obst][q2wbin,R2,seq,vst,var].GetBinError(ibin+1)
								stat_err2=hobs[OBS2,obst][q2wbin,R2,seq,vst,var].GetBinError(ibin+1)
								diff=calc_rel_diff(obs1,obs2)
								abs_diff=math.fabs(obs2-obs1)
								#! Fill hdiff
								hdiff[obst][q2wbin,R2,seq,vst,var].Fill(diff)
								#! Check if abs_diff is greater than stat errs
								abs_diff_gt_staterr= math.fabs(abs_diff)>3*stat_err1 and math.fabs(abs_diff)>3*stat_err2
								#! If abs_diff is greater than stat errs:
								#! 1. get hlindiff (#! hlindiff is basically histogram version of calc_rel_diff())
								#! 2. Fill hdiff_gt_staterr 
								#! 3. Write diff info to file
								if abs_diff_gt_staterr:
									#! 1. get hlindiff
									if hlindiff_drawn==False:
										hlindiff[obst][q2wbin,R2,seq,vst,var]=get_hlindiff(hobs[OBS1,obst][q2wbin,R2,seq,vst,var],hobs[OBS2,obst][q2wbin,R2,seq,vst,var])
										# if var=='M1' or var=='M2':
										# 	wmax=get_wbin_ue(q2wbin)
										# 	print "q2wbin,wamx=",q2wbin,wmax
										# 	setM1M2axisrange(hlindiff[obst][q2wbin,R2,seq,vst,var],vst,var,wmax)
										hlindiff_drawn=True
									#! 2. Fill hist_gt_staterr
									hdiff_gt_staterr[obst][q2wbin,R2,seq,vst,var].Fill(diff)
									#! 3. Write diff info to file
									#! note binn and binle
									binn=ibin+1
									binle=hobs[OBS1,obst][q2wbin,R2,seq,vst,var].GetBinLowEdge(ibin+1)
									fout[obst].write("*** diff(%), q2wb, R2(named), seq, vst, var, var_name, phi_name, binn, binle, obs1+/3*stat_err1, obs2+/-3*stat_err2:***\n")
									fout[obst].write("%.2f, %s, %s(%s), %s, %d, %s, %s, %s, %d, %.3f ,%f+/-%f , %f+/-%f\n"%
													(diff,q2wbin,R2,R2_NAMED[R2],seq,vst,var,VAR_NAMES_PLAIN[(vst,var)],VAR_NAMES_PLAIN[(vst,'PHI')],
													 binn,binle,obs1,3*stat_err1,obs2,3*stat_err2))
									if any([math.fabs(diff)>x for x in DIFF_THRSHLD_VALS]):
										thrshld_val=get_closest_diff_thrshld_val(diff)
										fout[obst].write("diff exceeded-threshold-%d\n"%thrshld_val)
									fout[obst].write("\n")
	#! close fout
	for obst in ['Obs_1D','Obs_R2']:
		fout[obst].close()

	#! Check integrity of hobs[obsn,obst]=f(q2wb,seq,vst,var)
	if DBG==True:
		nq2wbins=6 #! 3 each from lowQ2 and highQ2 ->3+3=6
	else: #! ana2pi nq2wbins 
		nq2wbins=5*29
	nhobs=len(hobs[OBS1,'Obs_1D'])+len(hobs[OBS2,'Obs_1D'])+len(hobs[OBS1,'Obs_R2'])+len(hobs[OBS2,'Obs_R2'])
	nhobs_true_1D=NOBS*nq2wbins*2*9  #!2,9=nseqs(EC,EF) + 9 1D xsecs (see thesis 12)
	nhobs_true_R2=NOBS*nq2wbins*2*42 #!2,42=nseqs(EC,EF)+ 42 R2 xsecs (see thesis 12)
	nhobs_true=nhobs_true_1D+nhobs_true_R2
	print "*** Integrity test for hobs ***"
	print "Number of objects in hobs=",nhobs
	print "True number of hobs=",nhobs_true
	if nhobs==nhobs_true:
		print "hobs PASSED integrity test"
	else:
		print "hobs FAILED integrity test"
	print "******"

	#! Plots

	#! Plots: hobs1,hobs2,hlindiff
	#! Plotting aesthetics for hlindiff
	ROOT.gStyle.Reset() #! Reset aesthetics which may have been set by some other plot function
	ROOT.gStyle.SetOptStat(0)

	#! hlindiff: Obs_1D
	#! Get list of q2wbins from hlindiff['Obs_1D']
	q2wbinl_1D=[k[0] for k in hlindiff['Obs_1D']]
	#! keep only unique values
	q2wbinl_1D=list(set(q2wbinl_1D))
	#print q2wbinl_1D
	#! Now loop over q2wbinl and plot
	#! + Note plot structure is that of Obs_1D
	#! + {hobs1,hobs2,hlindiff}*{EC,EF} are drawn together
	for q2wbin in q2wbinl_1D:
		#! prepare canvas for q2wbin
		c=ROOT.TCanvas("c","c",1000,1000)
		pad_t=ROOT.TPad("pad_l","Legend pad",0.25,0.935,0.75,1.00)
		pad_p=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
		#! create another transparent pad on which hlindiff will be drawn 
		#! using its own y-axis scale
		pad_p2=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
		pad_p2.SetFillStyle(4000)# transparent
		pad_p.Draw()
		pad_p2.Draw()
		pad_t.Draw()
		pad_t.cd()
		pt=ROOT.TPaveText(.05,.1,.95,.8)
		pt.AddText("Q2_W bin=%s"%q2wbin)
		pt.SetTextSize(0.40)
		pt.Draw()
		pad_p.Divide(3,3)
		pad_p2.Divide(3,3)
		#! Following for hist draw option with "same"
		pad_first_plot_drawn=OrderedDict()
		for pad in [1,2,3,4,5,6,7,8,9]:
			pad_first_plot_drawn[pad]=False
		#! Needed for separate yaxis for hlindiff
		axis=OrderedDict()
		#! legend
		l=OrderedDict()
		#! TLine at y=0 for hlindiff
		tl=OrderedDict()
		#! Now loop over hlindiff['Obs_1D'] and plot only where q2wbin1=q2wbin
		for k in hlindiff['Obs_1D']:
			q2wbin1,seq,vst,var=k[0],k[1],k[2],k[3]
			if q2wbin1!=q2wbin: continue
			for item in PAD_MAP['Obs_1D']:
				pad,vst_pm,var_pm=item[0],item[1],item[2]
				if vst_pm==vst and var_pm==var:
					crntpad=pad_p.cd(pad)
					#! first plot on pad
					draw_opt_first_plot=""
					if pad_first_plot_drawn[pad]==True: draw_opt_first_plot="same"
					#! aesthetics
					hobs[OBS1,'Obs_1D'][k].SetMarkerStyle(MRKS_SEQ[seq])
					hobs[OBS1,'Obs_1D'][k].SetMarkerColor(CLRS_OBS[OBS1])
					hobs[OBS1,'Obs_1D'][k].SetLineColor(CLRS_OBS[OBS1])
					#! draw
					hobs[OBS1,'Obs_1D'][k].Draw(draw_opt_first_plot)
					pad_first_plot_drawn[pad]=True

					#! rest of the plots: obs2
					#! aesthetics
					hobs[OBS2,'Obs_1D'][k].SetMarkerStyle(MRKS_SEQ[seq])
					hobs[OBS2,'Obs_1D'][k].SetMarkerColor(CLRS_OBS[OBS2])
					hobs[OBS2,'Obs_1D'][k].SetLineColor(CLRS_OBS[OBS2])
					#! draw
					hobs[OBS2,'Obs_1D'][k].Draw("same")

					#! rest of the plots: hlindiff
					#! + This needs a separate scale for the y-axis
					#! + Solution taken from https://root.cern.ch/root/html/tutorials/hist/transpad.C.html
					#! + Currently being plotted only when seq=='EF'
					if seq=='EF' or seq=='EC':
						ymin=-50 #hlindiff['Obs_1D'][k].GetMinimum()
						ymax=50 #hlindiff['Obs_1D'][k].GetMaximum()
						dy=(ymax-ymin)/0.8 #!10 per cent margins top and bottom
						if var=='M1' or var=='M2': #! because of setM1M2axisrange
							xmin=crntpad.GetUxmin()
							xmax=crntpad.GetUxmax()
						else:
							xmin=hlindiff['Obs_1D'][k].GetXaxis().GetXmin()
							xmax=hlindiff['Obs_1D'][k].GetXaxis().GetXmax()
						dx=(xmax-xmin)/0.8 #!10 per cent margins left and right
						crntpad2=pad_p2.cd(pad)
						if crntpad2==None:
							print "Null pointer at: proc_Obs_1D_and_R2(): Plot hlindiff 1D: crntpad2=pad_p2.cd(pad)\n"
							sys.exit("q2wbin1,seq,vst,var=\n",q2wbin1,seq,vst,var)
						crntpad2.SetFillStyle(4000)# transparent
						print "ymin,ymax,xmin,xmax=",ymin,ymax,xmin,xmax
						crntpad2.Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy)
						crntpad2.Draw()
						crntpad2.cd()
						hlindiff['Obs_1D'][k].SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
						hlindiff['Obs_1D'][k].SetMarkerColor(ROOT.gROOT.ProcessLine("kRed"));
						hlindiff['Obs_1D'][k].SetMarkerStyle(MRKS_SEQ[seq])
						hlindiff['Obs_1D'][k].Draw("p][same");
						#!draw axis on the right side of the pad
						axis[k]=ROOT.TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,50510,"+L");
						axis[k].SetLabelColor(ROOT.gROOT.ProcessLine("kRed"))
						axis[k].SetTitle("rel. diff. [%]")
						axis[k].SetTitleColor(ROOT.gROOT.ProcessLine("kRed"))
						axis[k].Draw()
						#! draw TLine at y=0
						tl[k]=ROOT.TLine(xmin,0,xmax,0)
						tl[k].SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
						tl[k].Draw()

					#! legend
					#! + Note the legend purposefully omits 'EC' entries to reduce crowd
					#! + It is implicit that kOpenCircle is for EC
					if seq=='EF':
						l[k]=ROOT.TLegend(0.7,0.80,0.85,0.90)
						l[k].SetFillStyle(0)
						l[k].SetBorderSize(0)
						l[k].SetTextSize(0.05)
						l[k].AddEntry(hobs[OBS1,'Obs_1D'][k],"obs1","p")
						l[k].AddEntry(hobs[OBS2,'Obs_1D'][k],"obs2","p")
						l[k].AddEntry(hlindiff['Obs_1D'][k],"rel-diff","p")
						l[k].Draw()
		#! save
		q2bin=get_q2bin(q2wbin)
		wbin=get_wbin(q2wbin)
		outdir="%s/%s"%(OUTDIR_1D,q2bin)
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		c.SaveAs("%s/c_%s.png"%(outdir,wbin))
		#print "outdir=",outdir
	

	#! hlindiff: Obs_R2
	# #print "Going to plot hlindiff:Obs_R2\n"
	#! Get list of unique combination of q2wbins,R2 from hlindiff['Obs_R2']
	q2wbin_R2l=[[k[0],k[1]] for k in hlindiff['Obs_R2']]
	#! keep only unique items
	q2wbin_R2l=return_unique_2D_list(q2wbin_R2l)
	# print "q2wbin_R2l"
	# print q2wbin_R2l
	# for q2wbin,R2 in q2wbin_R2l:
	# 	print q2wbin,R2
		# for k in hlindiff['Obs_R2']:
		# 	q2wbin1,R21,seq,vst,var=k[0],k[1],k[2],k[3],k[4]
		# 	if q2wbin1!=q2wbin or R21!=R2: continue
		# 	print "making canvas for q2wbin,R2,s"
	#sys.exit()
	#! Now loop over q2wbinl*R2 and plot
	#! + Note plot structure is that of Obs_R2
	#! + {hobs1,hobs2,hlindiff}*{EC,EF} are drawn together
	for q2wbin,R2 in q2wbin_R2l:
		#! prepare canvas for q2wbin/R2
		c=ROOT.TCanvas("c","c",1000,1000)
		pad_t=ROOT.TPad("pad_l","Legend pad",0.15,0.945,0.85,1.00)
		pad_p=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
		#! create another transparent pad on which hlindiff will be drawn 
		#! using its own y-axis scale
		pad_p2=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
		pad_p2.SetFillStyle(4000)# transparent
		pad_p.Draw()
		pad_p2.Draw()
		pad_t.Draw()
		pad_t.cd()
		pt=ROOT.TPaveText(.05,.1,.95,.8)
		pt.AddText("R2=%s, Q2_W bin=%s"%(R2,q2wbin))
		pt.SetTextSize(0.40)
		pt.Draw()
		pad_p.Divide(3,4)
		pad_p2.Divide(3,4)
		#! Following for hist draw option with "same"
		pad_first_plot_drawn=OrderedDict()
		for pad in [1,2,3,4,5,6,7,8,9,10,11,12]:
			pad_first_plot_drawn[pad]=False
		#! Needed for separate yaxis for hlindiff
		axis=OrderedDict()
		#! legend
		l=OrderedDict()
		#! draw TLine at y=0
		tl=OrderedDict()
		#! Now loop over hlindiff['Obs_R2'] and plot only where q2wbin1=q2wbin and R21=R2
		for k in hlindiff['Obs_R2']:
			q2wbin1,R21,seq,vst,var=k[0],k[1],k[2],k[3],k[4]
			#if q2wbin1!=q2wbin or R21!=R2: continue
			if not(q2wbin1==q2wbin and R21==R2): continue
			#if R21!=R2: continue
			for item in PAD_MAP['Obs_R2']:
				pad,vst_pm,var_pm=item[0],item[1],item[2]
				if vst_pm==vst and var_pm==var:
					crntpad=pad_p.cd(pad)
					#! first plot on pad
					draw_opt_first_plot=""
					if pad_first_plot_drawn[pad]==True: draw_opt_first_plot="same"
					#! aesthetics
					hobs[OBS1,'Obs_R2'][k].SetMarkerStyle(MRKS_SEQ[seq])
					hobs[OBS1,'Obs_R2'][k].SetMarkerColor(CLRS_OBS[OBS1])
					hobs[OBS1,'Obs_R2'][k].SetLineColor(CLRS_OBS[OBS1])
					#! draw
					hobs[OBS1,'Obs_R2'][k].Draw(draw_opt_first_plot)
					pad_first_plot_drawn[pad]=True

					#! rest of the plots: obs2
					#! aesthetics
					hobs[OBS2,'Obs_R2'][k].SetMarkerStyle(MRKS_SEQ[seq])
					hobs[OBS2,'Obs_R2'][k].SetMarkerColor(CLRS_OBS[OBS2])
					hobs[OBS2,'Obs_R2'][k].SetLineColor(CLRS_OBS[OBS2])
					#! draw
					hobs[OBS2,'Obs_R2'][k].Draw("same")

					#! rest of the plots: hlindiff
					#! + This needs a separate scale for the y-axis
					#! + Solution taken from https://root.cern.ch/root/html/tutorials/hist/transpad.C.html
					#! + Currently being plotted only when seq=='EF'
					if seq=='EF' or seq=='EC':
						ymin=-50 #hlindiff['Obs_1D'][k].GetMinimum()
						ymax=50 #hlindiff['Obs_1D'][k].GetMaximum()
						dy=(ymax-ymin)/0.8 #!10 per cent margins top and bottom
						if var=='M1' or var=='M2': #! because of setM1M2axisrange
							xmin=crntpad.GetUxmin()
							xmax=crntpad.GetUxmax()
						else:
							xmin=hlindiff['Obs_R2'][k].GetXaxis().GetXmin()
							xmax=hlindiff['Obs_R2'][k].GetXaxis().GetXmax()
						dx=(xmax-xmin)/0.8 #!10 per cent margins left and right
						crntpad2=pad_p2.cd(pad)
						if crntpad2==None:
							print "Null pointer at: proc_Obs_1D_and_R2(): Plot hlindiff R2: crntpad2=pad_p2.cd(pad)\n"
							sys.exit("q2wbin1,R21,seq,vst,var=\n",q2wbin1,R21,seq,vst,var)
						crntpad2.SetFillStyle(4000)# transparent
						print "ymin,ymax,xmin,xmax=",ymin,ymax,xmin,xmax
						crntpad2.Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy)
						crntpad2.Draw();
						crntpad2.cd();
						hlindiff['Obs_R2'][k].SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
						hlindiff['Obs_R2'][k].SetMarkerColor(ROOT.gROOT.ProcessLine("kRed"))
						hlindiff['Obs_R2'][k].SetMarkerStyle(MRKS_SEQ[seq])
						hlindiff['Obs_R2'][k].Draw("p][same");
						#!draw axis on the right side of the pad
						axis[k]=ROOT.TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,50510,"+L");
						axis[k].SetLabelColor(ROOT.gROOT.ProcessLine("kRed"))
						axis[k].SetTitle("rel. diff. [%]")
						axis[k].SetTitleColor(ROOT.gROOT.ProcessLine("kRed"))
						axis[k].Draw()
						#! draw TLine at y=0
						tl[k]=ROOT.TLine(xmin,0,xmax,0)
						tl[k].SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
						tl[k].Draw()

					#! legend
					#! + Note the legend purposefully omits 'EC' entries to reduce crowd
					#! + It is implicit that kOpenCircle is for EC
					if seq=='EF':
						l[k]=ROOT.TLegend(0.7,0.80,0.85,0.90)
						l[k].SetFillStyle(0)
						l[k].SetBorderSize(0)
						l[k].SetTextSize(0.05)
						l[k].AddEntry(hobs[OBS1,'Obs_R2'][k],"obs1","p")
						l[k].AddEntry(hobs[OBS2,'Obs_R2'][k],"obs2","p")
						l[k].AddEntry(hlindiff['Obs_R2'][k],"rel-diff","p")
						l[k].Draw()
		#! save
		q2bin=get_q2bin(q2wbin)
		wbin=get_wbin(q2wbin)
		outdir="%s/%s/%s"%(OUTDIR_R2,R2,q2bin)
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		c.SaveAs("%s/c_%s.png"%(outdir,wbin))
		print "outdir=",outdir
	#! end q2wbin*R2 loop
	#print "Done :Going to plot hlindiff:Obs_R2\n"
		
	#! Plots: hdiff
	#! Plotting aesthetics 
	ROOT.gStyle.Reset() #! Reset aesthetics which may have been set by some other plot function
	#! Keep track of under and overflow bins for cases when diff is in these limits, which it should not be!
	ROOT.gStyle.SetOptStat("neMRiuo")

	#! Define TLine that will be drawn at EXPCTD_DIFF
	#! + Should be able to use in all following plots
	expctd_diff_line=ROOT.TLine(EXPCTD_DIFF,0,EXPCTD_DIFF,0)
	expctd_diff_tlrnc_line_min=ROOT.TLine(EXPCTD_DIFF-EXPCTD_DIFF_TLRNC,0,EXPCTD_DIFF-EXPCTD_DIFF_TLRNC,0)
	expctd_diff_tlrnc_line_max=ROOT.TLine(EXPCTD_DIFF+EXPCTD_DIFF_TLRNC,0,EXPCTD_DIFF+EXPCTD_DIFF_TLRNC,0)
	for l in [expctd_diff_line,expctd_diff_tlrnc_line_min,expctd_diff_tlrnc_line_max]:
		l.SetLineColor(ROOT.gROOT.ProcessLine("kGreen"))
		l.SetLineStyle(4)
		l.SetLineWidth(1)

	#! For Obs_1D, make plot integrated over all q2wb,seq,vst,var
	#! Add hists for all q2wb,vst,seq,var
	i=0
	outdir=OUTDIR_1D
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	for k in hdiff['Obs_1D']:
		if i==0:
			hdiff_itg=hdiff['Obs_1D'][k].Clone("hdiff_itg")
			hdiff_itg.SetTitle("hdiff_itg")
			hdiff_itg_gt_staterr=hdiff_gt_staterr['Obs_1D'][k].Clone("hdiff_itg_gt_staterr")
			hdiff_itg_gt_staterr.SetTitle("hdiff_itg_gt_staterr")
		else: #! add
			hdiff_itg.Add(hdiff['Obs_1D'][k])
			hdiff_itg_gt_staterr.Add(hdiff_gt_staterr['Obs_1D'][k])
		i+=1
	#! hist aesthetics
	hdiff_itg.SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
	hdiff_itg_gt_staterr.SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
	hdiff_itg_gt_staterr.SetLineStyle(2)
	#! Now draw and save hdiff_itg
	c=ROOT.TCanvas(hdiff_itg.GetName(),hdiff_itg.GetName())
	hdiff_itg.Draw("hist")
	hdiff_itg_gt_staterr.Draw("hist sames")
	#! Reposition stat boxes
	#! for hdiff_itg
	c.Update()
	st_itg=hdiff_itg.GetListOfFunctions().FindObject("stats")
	st_itg.SetTextColor(ROOT.gROOT.ProcessLine("kBlue"))
	st_itg.SetX1NDC(0.60)
	st_itg.SetX2NDC(0.80)
	st_itg.SetY1NDC(0.70)
	st_itg.SetY2NDC(0.90)
	c.Update()
	#! for hdiff_itg_gt_staterr
	st_itg_gt_staterr=hdiff_itg_gt_staterr.GetListOfFunctions().FindObject("stats")
	st_itg_gt_staterr.SetTextColor(ROOT.gROOT.ProcessLine("kRed"))
	st_itg_gt_staterr.SetX1NDC(0.80)
	st_itg_gt_staterr.SetX2NDC(1.00)
	st_itg_gt_staterr.SetY1NDC(0.70)
	st_itg_gt_staterr.SetY2NDC(0.90)
	c.Update()

	#! Draw expctd_diff_line
	for l in [expctd_diff_line,expctd_diff_tlrnc_line_min,expctd_diff_tlrnc_line_max]:
		l.SetY1(hdiff_itg.GetMinimum())
		l.SetY2(hdiff_itg.GetMaximum())
		l.Draw("same")
	#! legend
	l=ROOT.TLegend(0.58,0.50,0.73,0.65)	
	l.SetFillStyle(0)
	l.SetBorderSize(0)
	l.SetTextSize(0.04)
	n_itg=hdiff_itg.GetEntries()
	n_itg_gt_staterr=hdiff_itg_gt_staterr.GetEntries()
	prcnt_gt_st_staterr=(n_itg_gt_staterr/n_itg)*100
	l.AddEntry(hdiff_itg,"all(%d bins)"%n_itg,"l")
	l.AddEntry(hdiff_itg_gt_staterr,"stat-sig(%d bins)=%.2f%%"%(n_itg_gt_staterr,prcnt_gt_st_staterr),"l")
	l.AddEntry(expctd_diff_line,"\"expected diff\"(%.2f%%)"%EXPCTD_DIFF,"l")
	l.Draw()
	#! save
	c.SaveAs("%s/%s.png"%(outdir,c.GetName()))

	#! For Obs_1D, for each seq, make plot integrated over all q2wb, vst, var
	for seq in ['EC','EF']:
		outdir=OUTDIR_1D
		# outdir="%s/%s"%(OUTDIR_1D,seq)
		# if not os.path.exists(outdir):
		# 	os.makedirs(outdir)
		#! Add hists for all q2wb,vst,var
		i=0
		for k in hdiff['Obs_1D']:
			q2wbin,seq1,vst,var=k[0],k[1],k[2],k[3]
			if seq1!=seq: continue
			if i==0: #! then make clone
				hdiff_itg=hdiff['Obs_1D'][k].Clone("hdiff_itg_%s"%seq)
				hdiff_itg.SetTitle("hdiff_itg_%s"%seq)
				hdiff_itg_gt_staterr=hdiff_gt_staterr['Obs_1D'][k].Clone("hdiff_itg_%s"%seq)
				hdiff_itg_gt_staterr.SetTitle("hdiff_itg_gt_staterr_%s"%seq)
			else: #! add
				hdiff_itg.Add(hdiff['Obs_1D'][k])
				hdiff_itg_gt_staterr.Add(hdiff_gt_staterr['Obs_1D'][k])
			i+=1
		#! hist aesthetics
		hdiff_itg.SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
		hdiff_itg_gt_staterr.SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
		hdiff_itg_gt_staterr.SetLineStyle(2)
		#! Now draw and save hdiff_itg
		c=ROOT.TCanvas(hdiff_itg.GetName(),hdiff_itg.GetName())
		hdiff_itg.Draw("hist")
		hdiff_itg_gt_staterr.Draw("hist sames")
		#! Reposition stat boxes
		#! for hdiff_itg
		c.Update()
		st_itg=hdiff_itg.GetListOfFunctions().FindObject("stats")
		st_itg.SetTextColor(ROOT.gROOT.ProcessLine("kBlue"))
		st_itg.SetX1NDC(0.60)
		st_itg.SetX2NDC(0.80)
		st_itg.SetY1NDC(0.70)
		st_itg.SetY2NDC(0.90)
		c.Update()
		#! for hdiff_itg_gt_staterr
		st_itg_gt_staterr=hdiff_itg_gt_staterr.GetListOfFunctions().FindObject("stats")
		st_itg_gt_staterr.SetTextColor(ROOT.gROOT.ProcessLine("kRed"))
		st_itg_gt_staterr.SetX1NDC(0.80)
		st_itg_gt_staterr.SetX2NDC(1.00)
		st_itg_gt_staterr.SetY1NDC(0.70)
		st_itg_gt_staterr.SetY2NDC(0.90)
		c.Update()
		#! draw lines for expected diff
		for l in [expctd_diff_line,expctd_diff_tlrnc_line_min,expctd_diff_tlrnc_line_max]:
			l.SetY1(hdiff_itg.GetMinimum())
			l.SetY2(hdiff_itg.GetMaximum())
			l.Draw("same")
		#! legend
		l=ROOT.TLegend(0.58,0.50,0.73,0.65)	
		l.SetFillStyle(0)
		l.SetBorderSize(0)
		l.SetTextSize(0.04)
		n_itg=hdiff_itg.GetEntries()
		n_itg_gt_staterr=hdiff_itg_gt_staterr.GetEntries()
		prcnt_gt_st_staterr=(n_itg_gt_staterr/n_itg)*100
		l.AddEntry(hdiff_itg,"all(%d bins)"%n_itg,"l")
		l.AddEntry(hdiff_itg_gt_staterr,"stat-sig(%d bins)=%.2f%%"%(n_itg_gt_staterr,prcnt_gt_st_staterr),"l")
		l.AddEntry(expctd_diff_line,"\"expected diff\"(%.2f%%)"%EXPCTD_DIFF,"l")
		l.Draw()
		#! save
		c.SaveAs("%s/%s.png"%(outdir,c.GetName()))

	#! For Obs_R2, make plot integrated over all q2wb,R2,seq,vst,var
	#! Add hists for all q2wb,vst,seq,var
	i=0
	outdir=OUTDIR_R2
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	for k in hdiff['Obs_R2']:
		if i==0:
			hdiff_itg=hdiff['Obs_R2'][k].Clone("hdiff_itg")
			hdiff_itg.SetTitle("hdiff_itg")
			hdiff_itg_gt_staterr=hdiff_gt_staterr['Obs_R2'][k].Clone("hdiff_itg_gt_staterr")
			hdiff_itg_gt_staterr.SetTitle("hdiff_itg_gt_staterr")
		else: #! add
			hdiff_itg.Add(hdiff['Obs_R2'][k])
			hdiff_itg_gt_staterr.Add(hdiff_gt_staterr['Obs_R2'][k])
		i+=1
	#! hist aesthetics
	hdiff_itg.SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
	hdiff_itg_gt_staterr.SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
	hdiff_itg_gt_staterr.SetLineStyle(2)
	#! Now draw and save hdiff_itg
	c=ROOT.TCanvas(hdiff_itg.GetName(),hdiff_itg.GetName())
	hdiff_itg.Draw("hist")
	hdiff_itg_gt_staterr.Draw("hist sames")
	#! Reposition stat boxes
	#! for hdiff_itg
	c.Update()
	st_itg=hdiff_itg.GetListOfFunctions().FindObject("stats")
	st_itg.SetTextColor(ROOT.gROOT.ProcessLine("kBlue"))
	st_itg.SetX1NDC(0.60)
	st_itg.SetX2NDC(0.80)
	st_itg.SetY1NDC(0.70)
	st_itg.SetY2NDC(0.90)
	c.Update()
	#! for hdiff_itg_gt_staterr
	st_itg_gt_staterr=hdiff_itg_gt_staterr.GetListOfFunctions().FindObject("stats")
	st_itg_gt_staterr.SetTextColor(ROOT.gROOT.ProcessLine("kRed"))
	st_itg_gt_staterr.SetX1NDC(0.80)
	st_itg_gt_staterr.SetX2NDC(1.00)
	st_itg_gt_staterr.SetY1NDC(0.70)
	st_itg_gt_staterr.SetY2NDC(0.90)
	c.Update()
	#! Draw expctd_diff_line
	for l in [expctd_diff_line,expctd_diff_tlrnc_line_min,expctd_diff_tlrnc_line_max]:
		l.SetY1(hdiff_itg.GetMinimum())
		l.SetY2(hdiff_itg.GetMaximum())
		l.Draw("same")
	#! legend
	l=ROOT.TLegend(0.58,0.50,0.73,0.65)	
	l.SetFillStyle(0)
	l.SetBorderSize(0)
	l.SetTextSize(0.04)
	n_itg=hdiff_itg.GetEntries()
	n_itg_gt_staterr=hdiff_itg_gt_staterr.GetEntries()
	prcnt_gt_st_staterr=(n_itg_gt_staterr/n_itg)*100
	l.AddEntry(hdiff_itg,"all(%d bins)"%n_itg,"l")
	l.AddEntry(hdiff_itg_gt_staterr,"stat-sig(%d bins)=%.2f%%"%(n_itg_gt_staterr,prcnt_gt_st_staterr),"l")
	l.AddEntry(expctd_diff_line,"\"expected diff\"(%.2f%%)"%EXPCTD_DIFF,"l")
	l.Draw()
	#! save
	c.SaveAs("%s/%s.png"%(outdir,c.GetName()))

	#! For Obs_R2, for each seq and R2, make plot integrated over all q2wb, vst, var
	for seq in ['EC','EF']:
		for R2 in R2S:
			outdir=OUTDIR_R2
			# outdir="%s/%s/%s"%(OUTDIR_R2,R2,seq)
			# if not os.path.exists(outdir):
			# 	os.makedirs(outdir)
			#! Add hists for all q2wb,vst,var
			i=0
			for k in hdiff['Obs_R2']:
				q2wbin,R21,seq1,vst,var=k[0],k[1],k[2],k[3],k[4]
				if R21!=R2 or seq1!=seq: continue
				if i==0: #! then make clone
					hdiff_itg=hdiff['Obs_R2'][k].Clone("hdiff_itg_%s_%s"%(R2,seq))
					hdiff_itg.SetTitle("hdiff_itg_%s_%s"%(R2,seq))
					hdiff_itg_gt_staterr=hdiff_gt_staterr['Obs_R2'][k].Clone("hdiff_itg_gt_staterr_%s_%s"%(R2,seq))
					hdiff_itg_gt_staterr.SetTitle("hdiff_itg_gt_staterr_%s_%s"%(R2,seq))
				else: #! add
					hdiff_itg.Add(hdiff['Obs_R2'][k])
					hdiff_itg_gt_staterr.Add(hdiff_gt_staterr['Obs_R2'][k])
				i+=1
			#! hist aesthetics
			hdiff_itg.SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
			hdiff_itg_gt_staterr.SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
			hdiff_itg_gt_staterr.SetLineStyle(2)
			#! Now draw and save hdiff_itg
			c=ROOT.TCanvas(hdiff_itg.GetName(),hdiff_itg.GetName())
			hdiff_itg.Draw("hist")
			hdiff_itg_gt_staterr.Draw("hist sames")
			#! Reposition stat boxes
			#! for hdiff_itg
			c.Update()
			st_itg=hdiff_itg.GetListOfFunctions().FindObject("stats")
			st_itg.SetTextColor(ROOT.gROOT.ProcessLine("kBlue"))
			st_itg.SetX1NDC(0.60)
			st_itg.SetX2NDC(0.80)
			st_itg.SetY1NDC(0.70)
			st_itg.SetY2NDC(0.90)
			c.Update()
			#! for hdiff_itg_gt_staterr
			st_itg_gt_staterr=hdiff_itg_gt_staterr.GetListOfFunctions().FindObject("stats")
			st_itg_gt_staterr.SetTextColor(ROOT.gROOT.ProcessLine("kRed"))
			st_itg_gt_staterr.SetX1NDC(0.80)
			st_itg_gt_staterr.SetX2NDC(1.00)
			st_itg_gt_staterr.SetY1NDC(0.70)
			st_itg_gt_staterr.SetY2NDC(0.90)
			c.Update()
			#! Draw expctd_diff_line
			for l in [expctd_diff_line,expctd_diff_tlrnc_line_min,expctd_diff_tlrnc_line_max]:
				expctd_diff_line.SetY1(hdiff_itg.GetMinimum())
				expctd_diff_line.SetY2(hdiff_itg.GetMaximum())
				expctd_diff_line.Draw("same")
			#! legend
			l=ROOT.TLegend(0.58,0.50,0.73,0.65)	
			l.SetFillStyle(0)
			l.SetBorderSize(0)
			l.SetTextSize(0.04)
			n_itg=hdiff_itg.GetEntries()
			n_itg_gt_staterr=hdiff_itg_gt_staterr.GetEntries()
			prcnt_gt_st_staterr=(n_itg_gt_staterr/n_itg)*100
			l.AddEntry(hdiff_itg,"all(%d bins)"%n_itg,"l")
			l.AddEntry(hdiff_itg_gt_staterr,"stat-sig(%d bins)=%.2f%%"%(n_itg_gt_staterr,prcnt_gt_st_staterr),"l")
			l.AddEntry(expctd_diff_line,"\"expected diff\"(%.2f%%)"%EXPCTD_DIFF,"l")
			l.Draw()
			#! save
			c.SaveAs("%s/%s.png"%(outdir,c.GetName()))
	print "*** Done proc_Obs_1D_and_R2() ***\n"
	print "If the function is not exiting, then Python is probably doing \"garbage collection\"(?); Wait a while!"

def proc_Obs_itg():
	print "*** In proc_Obs_itg() ***\n"

	#! create structure to contain:
	#! + hobs[obsn]=f(q2wb,seq) #! note vst is NA because it is averaged over and var=THETA
	#! + hdiff     =f(q2wb,seq)
	#! hobs
	hobs=OrderedDict()
	for obsn in range(NOBS): 
		hobs[obsn]=OrderedDict()
	#! hdiff
	#! + Note two sets of hists: one for all entries and only for those where diff>3*stat_err
	hdiff,hdiff_gt_staterr=OrderedDict(),OrderedDict()
	#! hlindiff
	#! + only for those where diff>3*stat_err
	hlindiff=OrderedDict()
			
	#! create file to store, in text form, any difference found to be greater than a threshold
	OUTDIR_ITG="%s/Obs_itg"%(OUTDIR)
	if not os.path.exists(OUTDIR_ITG):
		os.makedirs(OUTDIR_ITG)
	fout=open("%s/itg.txt"%(OUTDIR_ITG),'w')
	#! Descriptive opening remark in fout
	#fout.write("*** This file will note bins where diff in Obs_itg is greater than %.2f%% ***\n"%(EXPCTD_DIFF))
	fout.write("*** This file will note bins where absolute diff in Obs_itg is greater than stat errs ***\n")
	fout.write("\n")

	#! Now begin to fill this structure
	for q2 in ['lowQ2','highQ2']:
		
		#! open f[obsn]
		f=[0 for i in range(NOBS)]
		f[OBS1]=root_open(FIN[(OBS1,q2,'Obs_itg')])
		f[OBS2]=root_open(FIN[(OBS2,q2,'Obs_itg')])

		#! get q2wbinl from any of the files
		#! + DBG mode does not matter here because number of q2wbins are just 2 and 3 in lowQ2 and highQ2, respectively
		#! + This is because q2wbins for Obs_itg are seprated only by Q2 and all W bins are in a particular q2wbin
		q2wbinl=get_q2wbinlist(f[OBS1])
		if DBG==True:
			print "q2wbinl=",q2wbinl
		
		#! 2. For every qw2bin
		#! i. Get from f[obsn], for every seq in {EC,EF}, hobs[obsn]=f[q2wb,seq]
		for iq2wbin,q2wbin in enumerate(q2wbinl):
			#! Get q2bin, wbin information
			q2bin=get_q2bin(q2wbin)
			wbin=get_wbin(q2wbin)

			#! i. Get from f[obs], hobs[seq,vst,var]
			for seq in ['EC','EF']:
				hobs[OBS1][q2wbin,seq]=f[OBS1].Get("%s/hW_obs_%s_THETA"%(q2wbin,seq))
				hobs[OBS2][q2wbin,seq]=f[OBS2].Get("%s/hW_obs_%s_THETA"%(q2wbin,seq))
				#! Create and fill hdiff
				hname ="hdiff_%d_itg_%s"%(iq2wbin+1,seq)
				htitle="hdiff_%s_itg_%s"%(q2wbin,seq)
				hdiff[q2wbin,seq]=ROOT.TH1F(hname,htitle,NBINS,XMIN,XMAX)
				hdiff_gt_staterr[q2wbin,seq]=ROOT.TH1F("%s_gt_staterr"%hname,"%s_gt_staterr"%htitle,NBINS,XMIN,XMAX)
				#! loop over bins of hobs1 and fill ((hobs1-hobs2)/hobs1)*100
				nbins=hobs[OBS1][q2wbin,seq].GetNbinsX()
				#! + Initialize hlindiff_drawn=False.
				#! + This will be used later for cases where hlindiff=hobs1-hobs2 needs to be drawn
				#! + This is required because diffs are checked bin by bin, and if there is a diff the whole hist needs to be drawn
				#! + However, if there more bins are different, the hist does not need to be re-drawn
				hlindiff_drawn=False
				for ibin in range(nbins):
					obs1=hobs[OBS1][q2wbin,seq].GetBinContent(ibin+1)
					obs2=hobs[OBS2][q2wbin,seq].GetBinContent(ibin+1)
					stat_err1=hobs[OBS1][q2wbin,seq].GetBinError(ibin+1)
					stat_err2=hobs[OBS2][q2wbin,seq].GetBinError(ibin+1)
					diff=calc_rel_diff(obs1,obs2)
					abs_diff=math.fabs(obs2-obs1)
					#! Fill hdiff
					hdiff[q2wbin,seq].Fill(diff)
					#! Check if abs_diff is greater than stat errs
					abs_diff_gt_staterr= math.fabs(abs_diff)>3*stat_err1 and math.fabs(abs_diff)>3*stat_err2
					#! If abs_diff is greater than stat errs:
					#!    1. get hlindiff (#! hlindiff is basically histogram version of calc_rel_diff())
					#!    2. Fill hdiff_gt_staterr 
					#!    3. Write diff info to file
					if abs_diff_gt_staterr:
						#! 1. get hlindiff,
						if hlindiff_drawn==False:
							hlindiff[q2wbin,seq]=get_hlindiff(hobs[OBS1][q2wbin,seq],hobs[OBS2][q2wbin,seq])
						#! 2. Fill hdiff_gt_staterr
						hdiff_gt_staterr[q2wbin,seq].Fill(diff)
						#! note binn and binle
						binn=ibin+1
						binle=hobs[OBS1][q2wbin,seq].GetXaxis().GetLabels().At(ibin).GetName()
						#! 3. Write diff info to to file
						#fout.write("diff=%.2f for q2wbin,seq,binn,binle=%s,%s,%d,%s\n"%(diff,q2wbin,seq,binn,binle))
						fout.write("*** diff(%), q2wbin, seq, binn, binle, obs1+/-3*stat_err1, obs2+/-3*stat_err2: ***\n")
						fout.write("%.2f, %s, %s, %d, %s, %f+/-%f, %f+/-%f\n"%(diff,q2wbin,seq,binn,binle,obs1,3*stat_err1,obs2,3*stat_err2))
						# if math.fabs(diff)>3*stat_err1:
						# 	fout[obst].write("diff is greater than 3*stat_err1\n")
						# if math.fabs(diff)>3*stat_err2:
						# 	fout[obst].write("diff is greater than 3*stat_err2\n")
						#! check if diff greater than values in DIFF_THRSHLD_VALS
						if any([math.fabs(diff)>x for x in DIFF_THRSHLD_VALS]):
							thrshld_val=get_closest_diff_thrshld_val(diff)
							fout.write("diff exceeded-threshold-%d\n"%thrshld_val)
						fout.write("\n")
	#! close fout
	fout.close()

	#! Check integrity of hobs[obsn]=f(q2wb,seq,vst,var)
	nq2wbins=5
	nhobs=len(hobs[OBS1])+len(hobs[OBS2])
	nhobs_true=NOBS*nq2wbins*2*1  #!2,1=nseqs(EC,EF) + 1 itg xsec
	print "*** Integrity test for hobs ***"
	print "Number of objects in hobs=",nhobs
	print "True number of hobs=",nhobs_true
	if nhobs==nhobs_true:
		print "hobs PASSED integrity test"
	else:
		print "hobs FAILED integrity test"
	print "******"

	#! Plots
	#! Plots:hlindiff
	#! Plotting aesthetics for hlindiff
	ROOT.gStyle.Reset() #! Reset aesthetics which may have been set by some other plot function
	ROOT.gStyle.SetOptStat(0)
	#! Get list of q2wbins from hlindiff['Obs_1D']
	q2wbinl=[k[0] for k in hlindiff]
	#! keep only unique values
	q2wbinl=list(set(q2wbinl))
	#! Now loop through q2wbinl and plot
	for q2wbin in q2wbinl:
		c=ROOT.TCanvas()
		pad=ROOT.TPad("pad1","",0,0,1,1)
		#! pad2 for hlindiff which has a differnet yaxis scale 
		pad2=ROOT.TPad("pad1","",0,0,1,1)
		pad2.SetFillStyle(4000)#! transparent
		pad.Draw()
		pad2.Draw()
		#! Following for hist draw option with "same"
		pad_first_plot_drawn=False
		#! axis
		axis=OrderedDict()
		#! legend
		l=OrderedDict()
		#! TLine at y=0
		tl=OrderedDict()
		for k in hlindiff:
			q2wbin1,seq=k[0],k[1]
			if q2wbin1!=q2wbin: continue
			#! draw obs1
			pad.cd()
			#! aesthetics
			hobs[OBS1][k].SetTitle("%s"%q2wbin1)
			hobs[OBS1][k].SetMarkerStyle(MRKS_SEQ[seq])
			hobs[OBS1][k].SetMarkerColor(CLRS_OBS[OBS1])
			hobs[OBS1][k].SetLineColor(CLRS_OBS[OBS1])
			draw_opt_first_plot="same"
			if pad_first_plot_drawn==False:draw_opt_first_plot=""
			hobs[OBS1][k].Draw(draw_opt_first_plot)
			pad_first_plot_drawn=True

			#! draw obs2
			#! aesthetics
			hobs[OBS2][k].SetMarkerStyle(MRKS_SEQ[seq])
			hobs[OBS2][k].SetMarkerColor(CLRS_OBS[OBS2])
			hobs[OBS2][k].SetLineColor(CLRS_OBS[OBS2])
			hobs[OBS2][k].Draw("same")

			#! hlindiff
			pad2.cd()
			if seq=="EF" or seq=='EC':
				ymin=-50
				ymax=50
				dy=(ymax-ymin)/0.8 #!10 per cent margins top and bottom
				xmin=hlindiff[k].GetXaxis().GetXmin()
				xmax=hlindiff[k].GetXaxis().GetXmax()
				dx=(xmax-xmin)/0.8 #!10 per cent margins left and right
				pad2.SetFillStyle(4000)# transparent
				print "ymin,ymax,xmin,xmax=",ymin,ymax,xmin,xmax
				pad2.Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy)
				pad2.Draw();
				pad2.cd();
				hlindiff[k].SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
				hlindiff[k].SetMarkerColor(ROOT.gROOT.ProcessLine("kRed"))
				hlindiff[k].SetMarkerStyle(MRKS_SEQ[seq])
				hlindiff[k].Draw("p][same")
				#!draw axis on the right side of the pad
				axis[k]=ROOT.TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,50510,"+L")
				axis[k].SetLabelColor(ROOT.gROOT.ProcessLine("kRed"))
				axis[k].SetTitle("rel. diff. [%]")
				axis[k].SetTitleColor(ROOT.gROOT.ProcessLine("kRed"))
				axis[k].Draw()
				#! draw TLine at y=0
				tl[k]=ROOT.TLine(xmin,0,xmax,0)
				tl[k].SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
				tl[k].Draw()

			#! legend
			if seq=='EF':
				l[k]=ROOT.TLegend(0.7,0.80,0.85,0.90)
				l[k].SetFillStyle(0)
				l[k].SetBorderSize(0)
				l[k].SetTextSize(0.05)
				l[k].AddEntry(hobs[OBS1][k],"obs1","p")
				l[k].AddEntry(hobs[OBS2][k],"obs2","p")
				l[k].AddEntry(hlindiff[k],"rel-diff","p")
				l[k].Draw()
		#! save
		q2bin=get_q2bin(q2wbin)
		wbin=get_wbin(q2wbin)
		outdir="%s/%s"%(OUTDIR_ITG,q2bin)
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		c.SaveAs("%s/c_%s.png"%(outdir,wbin))


	#! Plots: hdiff
	#! Plotting aesthetics
	ROOT.gStyle.Reset() #! Reset aesthetics which may have been set by some other plot function
	#! Keep track of under and overflow bins for cases when diff is in these limits, which it should not be!
	ROOT.gStyle.SetOptStat("neiMRuo")

	#! Define TLine that will be drawn at EXPCTD_DIFF
	#! + Should be able to use in all following plots
	expctd_diff_line=ROOT.TLine(EXPCTD_DIFF,0,EXPCTD_DIFF,0)
	expctd_diff_tlrnc_line_min=ROOT.TLine(EXPCTD_DIFF-EXPCTD_DIFF_TLRNC,0,EXPCTD_DIFF-EXPCTD_DIFF_TLRNC,0)
	expctd_diff_tlrnc_line_max=ROOT.TLine(EXPCTD_DIFF+EXPCTD_DIFF_TLRNC,0,EXPCTD_DIFF+EXPCTD_DIFF_TLRNC,0)
	for l in [expctd_diff_line,expctd_diff_tlrnc_line_min,expctd_diff_tlrnc_line_max]:
		l.SetLineColor(ROOT.gROOT.ProcessLine("kGreen"))
		l.SetLineStyle(4)
		l.SetLineWidth(1)

	#! Make plot integrated over all q2wb,seq
	#! Add hists for all q2wb,seq
	i=0
	outdir=OUTDIR_ITG
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	for k in hdiff:
		if i==0:
			hdiff_itg=hdiff[k].Clone("hdiff_itg")
			hdiff_itg.SetTitle("hdiff_itg")
			hdiff_itg_gt_staterr=hdiff_gt_staterr[k].Clone("hdiff_itg_gt_staterr")
			hdiff_itg_gt_staterr.SetTitle("hdiff_itg_gt_staterr")
		else: #! add
			hdiff_itg.Add(hdiff[k])
			hdiff_itg_gt_staterr.Add(hdiff_gt_staterr[k])
		i+=1
	#! hist aesthetics
	hdiff_itg.SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
	hdiff_itg_gt_staterr.SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
	hdiff_itg_gt_staterr.SetLineStyle(2)
	#! Now draw and save hdiff_itg
	c=ROOT.TCanvas(hdiff_itg.GetName(),hdiff_itg.GetName())
	hdiff_itg.Draw("hist")
	hdiff_itg_gt_staterr.Draw("hist sames")
	#! Reposition stat boxes
	#! for hdiff_itg
	c.Update()
	st_itg=hdiff_itg.GetListOfFunctions().FindObject("stats")
	st_itg.SetTextColor(ROOT.gROOT.ProcessLine("kBlue"))
	st_itg.SetX1NDC(0.60)
	st_itg.SetX2NDC(0.80)
	st_itg.SetY1NDC(0.70)
	st_itg.SetY2NDC(0.90)
	c.Update()
	#! for hdiff_itg_gt_staterr
	st_itg_gt_staterr=hdiff_itg_gt_staterr.GetListOfFunctions().FindObject("stats")
	st_itg_gt_staterr.SetTextColor(ROOT.gROOT.ProcessLine("kRed"))
	st_itg_gt_staterr.SetX1NDC(0.80)
	st_itg_gt_staterr.SetX2NDC(1.00)
	st_itg_gt_staterr.SetY1NDC(0.70)
	st_itg_gt_staterr.SetY2NDC(0.90)
	c.Update()
	#! Draw expctd_diff_line
	for l in [expctd_diff_line,expctd_diff_tlrnc_line_min,expctd_diff_tlrnc_line_max]:
		l.SetY1(hdiff_itg.GetMinimum())
		l.SetY2(hdiff_itg.GetMaximum())
		l.Draw("same")
	#! legend
	l=ROOT.TLegend(0.58,0.50,0.73,0.65)	
	l.SetFillStyle(0)
	l.SetBorderSize(0)
	l.SetTextSize(0.04)
	n_itg=hdiff_itg.GetEntries()
	n_itg_gt_staterr=hdiff_itg_gt_staterr.GetEntries()
	prcnt_gt_st_staterr=(n_itg_gt_staterr/n_itg)*100
	l.AddEntry(hdiff_itg,"all(%d bins)"%n_itg,"l")
	l.AddEntry(hdiff_itg_gt_staterr,"stat-sig(%d bins)=%.2f%%"%(n_itg_gt_staterr,prcnt_gt_st_staterr),"l")
	l.AddEntry(expctd_diff_line,"\"expected diff\"(%.2f%%)"%EXPCTD_DIFF,"l")
	l.Draw()
	#! save
	c.SaveAs("%s/%s.png"%(outdir,c.GetName()))

	#! Now, for each seq, make plot integrated over all q2wb
	for seq in ['EC','EF']:
		outdir=OUTDIR_ITG
		# outdir="%s/%s"%(OUTDIR_ITG,seq)
		# if not os.path.exists(outdir):
		# 	os.makedirs(outdir)
		#! Add hists for all q2wb,vst,var
		i=0
		for k in hdiff:
			q2wbin,seq1=k[0],k[1]
			if seq1!=seq: continue
			if i==0: #! then make clone
				hdiff_itg=hdiff[k].Clone("hdiff_itg_%s"%seq)
				hdiff_itg.SetTitle("hdiff_itg_%s"%seq)
				hdiff_itg_gt_staterr=hdiff_gt_staterr[k].Clone("hdiff_itg_gt_staterr_%s"%seq)
				hdiff_itg_gt_staterr.SetTitle("hdiff_itg_gt_staterr_%s"%seq)
			else: #! add
				hdiff_itg.Add(hdiff[k])
				hdiff_itg_gt_staterr.Add(hdiff_gt_staterr[k])
			i+=1
		#! hist aesthetics
		hdiff_itg.SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
		hdiff_itg_gt_staterr.SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
		hdiff_itg_gt_staterr.SetLineStyle(2)
		#! Now draw and save hdiff_itg
		c=ROOT.TCanvas(hdiff_itg.GetName(),hdiff_itg.GetName())
		hdiff_itg.Draw("hist")
		hdiff_itg_gt_staterr.Draw("hist sames")
		#! Reposition stat boxes
		#! for hdiff_itg
		c.Update()
		st_itg=hdiff_itg.GetListOfFunctions().FindObject("stats")
		st_itg.SetTextColor(ROOT.gROOT.ProcessLine("kBlue"))
		st_itg.SetX1NDC(0.60)
		st_itg.SetX2NDC(0.80)
		st_itg.SetY1NDC(0.70)
		st_itg.SetY2NDC(0.90)
		c.Update()
		#! for hdiff_itg_gt_staterr
		st_itg_gt_staterr=hdiff_itg_gt_staterr.GetListOfFunctions().FindObject("stats")
		st_itg_gt_staterr.SetTextColor(ROOT.gROOT.ProcessLine("kRed"))
		st_itg_gt_staterr.SetX1NDC(0.80)
		st_itg_gt_staterr.SetX2NDC(1.00)
		st_itg_gt_staterr.SetY1NDC(0.70)
		st_itg_gt_staterr.SetY2NDC(0.90)
		c.Update()
		#! Draw expctd_diff_line
		for l in [expctd_diff_line,expctd_diff_tlrnc_line_min,expctd_diff_tlrnc_line_max]:
			l.SetY1(hdiff_itg.GetMinimum())
			l.SetY2(hdiff_itg.GetMaximum())
			l.Draw("same")
		#! legend
		l=ROOT.TLegend(0.58,0.50,0.73,0.65)	
		l.SetFillStyle(0)
		l.SetBorderSize(0)
		l.SetTextSize(0.04)
		n_itg=hdiff_itg.GetEntries()
		n_itg_gt_staterr=hdiff_itg_gt_staterr.GetEntries()
		prcnt_gt_st_staterr=(n_itg_gt_staterr/n_itg)*100
		l.AddEntry(hdiff_itg,"all(%d bins)"%n_itg,"l")
		l.AddEntry(hdiff_itg_gt_staterr,"stat-sig(%d bins)=%.2f%%"%(n_itg_gt_staterr,prcnt_gt_st_staterr),"l")
		l.AddEntry(expctd_diff_line,"\"expected diff\"(%.2f%%)"%EXPCTD_DIFF,"l")
		l.Draw()
		#! save
		c.SaveAs("%s/%s.png"%(outdir,c.GetName()))
	print "*** Done proc_Obs_itg() ***\n"	
	print "If the function is not exiting, then Python is probably doing \"garbage collection\"(?); Wait a while!"

def get_closest_diff_thrshld_val(diff):
	'''
	Returns value from DIFF_THRSHLD_VALS that diff is closest to
	'''
	return min(DIFF_THRSHLD_VALS, key=lambda x:math.fabs(x-math.fabs(diff)))


#! main
proc_Obs_1D_and_R2()
proc_Obs_itg()