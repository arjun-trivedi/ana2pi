#!/usr/bin/python

import sys,os
from collections import OrderedDict

from rootpy.io import root_open, DoesNotExist
import itertools

import ROOT

import math

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
obsname1=OBSDIR[OBS1].split("/")[-1]
obsname2=OBSDIR[OBS2].split("/")[-1]
obsnamediff='diff_%s_%s'%(obsname1,obsname2)
#! use same path for OUTDIR as OBS2
path=OBSDIR[OBS2].replace(obsname2,'')
#print path
OUTDIR="%s/%s"%(path,obsnamediff)
if DBG:
	OUTDIR="%s/%s_dbg"%(path,obsnamediff)
if not os.path.exists(OUTDIR):
	os.makedirs(OUTDIR)
print "OUTDIR=",OUTDIR

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

R2L=['A','B','C','D','E']
R2_NAMED={'A':'R2_{T}+R2_{L}','B':'R2^{c}_{LT}','C':'R2^{c}_{TT}','D':'R2^{s}_{LT}','E':'R2^{s}_{TT}'}

#! Define binning for hdiff hists that will be created
NBINS=200
XMIN,XMAX=-100,100

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

def calc_diff(obs1,obs2):
	'''
	+ diff(%) is calculated relative to obs1: ((obs2-obs1)/obs1)*100
	+ special cases are noted
	'''
	if obs1==0:
		if obs2==0: #! then diff =0
			diff=0
		else: #! if obs2!=0, then diff has to be calculated relative to obs2, which is always 100%
			diff=((obs2-obs1)/obs2)*100 #! should always work out to 100%
	else:
		diff=((obs2-obs1)/obs1)*100
	return diff

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
	hdiff=OrderedDict()
	hdiff['Obs_1D']=OrderedDict()
	hdiff['Obs_R2']=OrderedDict()

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
		fout[obst].write("*** This file will note bins where diff in %s is greater than %.2f%% ***\n"%(obst,EXPCTD_DIFF))

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

			#! i. Get from f_RSLT, hobs[seq,vst,var]
			for seq in ['EC','EF']:
				for pm in PAD_MAP[obst]:
					vst,var=pm[1],pm[2]
					if obst=='Obs_1D': #! for 1D, SF was used instead of ST (SF=ST because of acceptance corr. and hole-filling)
						hobs[OBS1,obst][q2wbin,seq,vst,var]=f[OBS1].Get("%s/h1_%s_%d_%s"%(q2wbin,seq,vst,var))
						hobs[OBS2,obst][q2wbin,seq,vst,var]=f[OBS2].Get("%s/h1_%s_%d_%s"%(q2wbin,seq,vst,var))
						#! Create and fill hdiff
						hname ="hdiff_%d_1D_%s_VST%d_%s"%(iq2wbin+1,seq,vst,var)
						htitle="hdiff_%s_1D_%s_VST%d_%s"%(q2wbin,seq,vst,var)
						hdiff[obst][q2wbin,seq,vst,var]=ROOT.TH1F(hname,htitle,NBINS,XMIN,XMAX)
						#! loop over bins of hobs1 and fill ((hobs1-hobs2)/hobs1)*100
						nbins=hobs[OBS1,obst][q2wbin,seq,vst,var].GetNbinsX()
						for ibin in range(nbins):
							obs1=hobs[OBS1,obst][q2wbin,seq,vst,var].GetBinContent(ibin+1)
							obs2=hobs[OBS2,obst][q2wbin,seq,vst,var].GetBinContent(ibin+1)
							diff=calc_diff(obs1,obs2)
							#! Fill hdiff
							hdiff[obst][q2wbin,seq,vst,var].Fill(diff)
							#! note if diff is too large
							#if not (math.fabs(EXPCTD_DIFF)-EXPCTD_DIFF_TLRNC) < math.fabs(diff) < (math.fabs(EXPCTD_DIFF)+EXPCTD_DIFF_TLRNC):
							if math.fabs(diff) > (math.fabs(EXPCTD_DIFF)+EXPCTD_DIFF_TLRNC): 
								#! note binn and binle
								binn=ibin+1
								binle=hobs[OBS1,obst][q2wbin,seq,vst,var].GetBinLowEdge(ibin+1)
								#! Now write to to file
								fout[obst].write("diff=%.2f for q2wb,seq,vst,var,binn,binle=%s,%s,%d,%s,%d,%.3f\n"%(diff,q2wbin,seq,vst,var,binn,binle))
					elif obst=='Obs_R2':
						for R2 in R2L:
							if (R2=='D' or R2=='E') and var != 'ALPHA': continue
							hobs[OBS1,obst][q2wbin,R2,seq,vst,var]=f[OBS1].Get("%s/hR2_%s_%s_%d_%s"%(q2wbin,R2,seq,vst,var))
							hobs[OBS2,obst][q2wbin,R2,seq,vst,var]=f[OBS2].Get("%s/hR2_%s_%s_%d_%s"%(q2wbin,R2,seq,vst,var))
							#! Create and fill hdiff
							hname ="hdiff_%d_%s_%s_VST%d_%s"%(iq2wbin+1,R2,seq,vst,var)
							htitle="hdiff_%s_%s_%s_VST%d_%s"%(q2wbin,R2,seq,vst,var)
							hdiff[obst][q2wbin,R2,seq,vst,var]=ROOT.TH1F(hname,htitle,NBINS,XMIN,XMAX)
							#! loop over bins of hobs1 and fill ((hobs1-hobs2)/hobs1)*100
							nbins=hobs[OBS1,obst][q2wbin,R2,seq,vst,var].GetNbinsX()
							for ibin in range(nbins):
								obs1=hobs[OBS1,obst][q2wbin,R2,seq,vst,var].GetBinContent(ibin+1)
								obs2=hobs[OBS2,obst][q2wbin,R2,seq,vst,var].GetBinContent(ibin+1)
								diff=calc_diff(obs1,obs2)
								#! Fill hdiff
								hdiff[obst][q2wbin,R2,seq,vst,var].Fill(diff)
								#! note if diff is too large
								#if not (math.fabs(EXPCTD_DIFF)-EXPCTD_DIFF_TLRNC) < math.fabs(diff) < (math.fabs(EXPCTD_DIFF)+EXPCTD_DIFF_TLRNC):
								if math.fabs(diff) > (math.fabs(EXPCTD_DIFF)+EXPCTD_DIFF_TLRNC):
									#! note binn and binle
									binn=ibin+1
									binle=hobs[OBS1,obst][q2wbin,R2,seq,vst,var].GetBinLowEdge(ibin+1)
									fout[obst].write("diff=%.2f for q2wb,seq,vst,var,binn,binle=%s,%s,%d,%s,%d,%.3f\n"%(diff,q2wbin,seq,vst,var,binn,binle))
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
		else: #! add
			hdiff_itg.Add(hdiff['Obs_1D'][k])
		i+=1
	#! Now draw and save hdiff_itg
	c=ROOT.TCanvas(hdiff_itg.GetName(),hdiff_itg.GetName())
	hdiff_itg.Draw("hist")
	#! Draw expctd_diff_line
	for l in [expctd_diff_line,expctd_diff_tlrnc_line_min,expctd_diff_tlrnc_line_max]:
		l.SetY1(hdiff_itg.GetMinimum())
		l.SetY2(hdiff_itg.GetMaximum())
		l.Draw("same")
	c.SaveAs("%s/%s.png"%(outdir,c.GetName()))

	#! For Obs_1D, for each seq, make plot integrated over all q2wb, vst, var
	for seq in ['EC','EF']:
		outdir="%s/%s"%(OUTDIR_1D,seq)
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		#! Add hists for all q2wb,vst,var
		i=0
		for k in hdiff['Obs_1D']:
			q2wbin,seq1,vst,var=k[0],k[1],k[2],k[3]
			if seq1!=seq: continue
			if i==0: #! then make clone
				hdiff_itg=hdiff['Obs_1D'][k].Clone("hdiff_itg_%s"%seq)
				hdiff_itg.SetTitle("hdiff_itg_%s"%seq)
			else: #! add
				hdiff_itg.Add(hdiff['Obs_1D'][k])
			i+=1
		#! Now draw and save hdiff_itg
		c=ROOT.TCanvas(hdiff_itg.GetName(),hdiff_itg.GetName())
		hdiff_itg.Draw("hist")
		for l in [expctd_diff_line,expctd_diff_tlrnc_line_min,expctd_diff_tlrnc_line_max]:
			l.SetY1(hdiff_itg.GetMinimum())
			l.SetY2(hdiff_itg.GetMaximum())
			l.Draw("same")
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
		else: #! add
			hdiff_itg.Add(hdiff['Obs_R2'][k])
		i+=1
	#! Now draw and save hdiff_itg
	c=ROOT.TCanvas(hdiff_itg.GetName(),hdiff_itg.GetName())
	hdiff_itg.Draw("hist")
	#! Draw expctd_diff_line
	for l in [expctd_diff_line,expctd_diff_tlrnc_line_min,expctd_diff_tlrnc_line_max]:
		l.SetY1(hdiff_itg.GetMinimum())
		l.SetY2(hdiff_itg.GetMaximum())
		l.Draw("same")
	c.SaveAs("%s/%s.png"%(outdir,c.GetName()))

	#! For Obs_R2, for each seq and R2, make plot integrated over all q2wb, vst, var
	for seq in ['EC','EF']:
		for R2 in R2L:
			outdir="%s/%s/%s"%(OUTDIR_R2,R2,seq)
			if not os.path.exists(outdir):
				os.makedirs(outdir)
			#! Add hists for all q2wb,vst,var
			i=0
			for k in hdiff['Obs_R2']:
				q2wbin,R21,seq1,vst,var=k[0],k[1],k[2],k[3],k[4]
				if R21!=R2 or seq1!=seq: continue
				if i==0: #! then make clone
					hdiff_itg=hdiff['Obs_R2'][k].Clone("hdiff_itg_%s_%s"%(R2,seq))
					hdiff_itg.SetTitle("hdiff_itg_%s_%s"%(R2,seq))
				else: #! add
					hdiff_itg.Add(hdiff['Obs_R2'][k])
				i+=1
			#! Now draw and save hdiff_itg
			c=ROOT.TCanvas(hdiff_itg.GetName(),hdiff_itg.GetName())
			hdiff_itg.Draw("hist")
			#! Draw expctd_diff_line
			expctd_diff_line.SetY1(hdiff_itg.GetMinimum())
			expctd_diff_line.SetY2(hdiff_itg.GetMaximum())
			expctd_diff_line.Draw("same")
			c.SaveAs("%s/%s.png"%(outdir,c.GetName()))
	print "*** Done proc_Obs_1D_and_R2() ***\n"
	print "If the function is not exiting, then Python is probably doing \"garbage collection\"(?); Wait a while!"

def proc_Obs_itg():
	print "*** In proc_Obs_itg ***\n"

	#! create structure to contain:
	#! + hobs[obsn]=f(q2wb,seq) #! note vst is NA because it is averaged over and var=THETA
	#! + hdiff     =f(q2wb,seq)
	#! hobs
	hobs=OrderedDict()
	for obsn in range(NOBS): 
		hobs[obsn]=OrderedDict()
	#! hdiff
	hdiff=OrderedDict()
	
	#! create file to store, in text form, any difference found to be greater than a threshold
	OUTDIR_ITG="%s/Obs_itg"%(OUTDIR)
	if not os.path.exists(OUTDIR_ITG):
		os.makedirs(OUTDIR_ITG)
	fout=open("%s/itg.txt"%(OUTDIR_ITG),'w')
	#! Descriptive opening remark in fout
	fout.write("*** This file will note bins where diff in Obs_itg is greater than %.2f%% ***\n"%(EXPCTD_DIFF))

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

			#! i. Get from f_RSLT, hobs[seq,vst,var]
			for seq in ['EC','EF']:
				hobs[OBS1][q2wbin,seq]=f[OBS1].Get("%s/hW_obs_%s_THETA"%(q2wbin,seq))
				hobs[OBS2][q2wbin,seq]=f[OBS2].Get("%s/hW_obs_%s_THETA"%(q2wbin,seq))
				#! Create and fill hdiff
				hname ="hdiff_%d_itg_%s"%(iq2wbin+1,seq)
				htitle="hdiff_%s_itg_%s"%(q2wbin,seq)
				hdiff[q2wbin,seq]=ROOT.TH1F(hname,htitle,NBINS,XMIN,XMAX)
				#! loop over bins of hobs1 and fill ((hobs1-hobs2)/hobs1)*100
				nbins=hobs[OBS1][q2wbin,seq].GetNbinsX()
				for ibin in range(nbins):
					obs1=hobs[OBS1][q2wbin,seq].GetBinContent(ibin+1)
					obs2=hobs[OBS2][q2wbin,seq].GetBinContent(ibin+1)
					diff=calc_diff(obs1,obs2)
					#! Fill hdiff
					hdiff[q2wbin,seq].Fill(diff)
					#! note if diff is too large
					#if not (math.fabs(EXPCTD_DIFF)-EXPCTD_DIFF_TLRNC) < math.fabs(diff) < (math.fabs(EXPCTD_DIFF)+EXPCTD_DIFF_TLRNC):
					if math.fabs(diff) > (math.fabs(EXPCTD_DIFF)+EXPCTD_DIFF_TLRNC):
						#! note binn and binle
						binn=ibin+1
						binle=hobs[OBS1][q2wbin,seq].GetXaxis().GetLabels().At(ibin).GetName()
						#! Now write to to file
						fout.write("diff=%.2f for q2wbin,seq,binn,binle=%s,%s,%d,%s\n"%(diff,q2wbin,seq,binn,binle))
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
		else: #! add
			hdiff_itg.Add(hdiff[k])
		i+=1
	#! Now draw and save hdiff_itg
	c=ROOT.TCanvas(hdiff_itg.GetName(),hdiff_itg.GetName())
	hdiff_itg.Draw("hist")
	#! Draw expctd_diff_line
	for l in [expctd_diff_line,expctd_diff_tlrnc_line_min,expctd_diff_tlrnc_line_max]:
		l.SetY1(hdiff_itg.GetMinimum())
		l.SetY2(hdiff_itg.GetMaximum())
		l.Draw("same")
	c.SaveAs("%s/%s.png"%(outdir,c.GetName()))

	#! Now, for each seq, make plot integrated over all q2wb
	for seq in ['EC','EF']:
		outdir="%s/%s"%(OUTDIR_ITG,seq)
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		#! Add hists for all q2wb,vst,var
		i=0
		for k in hdiff:
			q2wbin,seq1=k[0],k[1]
			if seq1!=seq: continue
			if i==0: #! then make clone
				hdiff_itg=hdiff[k].Clone("hdiff_itg_%s"%seq)
				hdiff_itg.SetTitle("hdiff_itg_%s"%seq)
			else: #! add
				hdiff_itg.Add(hdiff[k])
			i+=1
		#! Now draw and save hdiff_itg
		c=ROOT.TCanvas(hdiff_itg.GetName(),hdiff_itg.GetName())
		hdiff_itg.Draw("hist")
		#! Draw expctd_diff_line
		for l in [expctd_diff_line,expctd_diff_tlrnc_line_min,expctd_diff_tlrnc_line_max]:
			l.SetY1(hdiff_itg.GetMinimum())
			l.SetY2(hdiff_itg.GetMaximum())
			l.Draw("same")
		c.SaveAs("%s/%s.png"%(outdir,c.GetName()))
	print "*** Done proc_Obs_itg ***\n"	
	print "If the function is not exiting, then Python is probably doing \"garbage collection\"(?); Wait a while!"

#! main
proc_Obs_1D_and_R2()
proc_Obs_itg()