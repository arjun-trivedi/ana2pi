#!/usr/bin/python
from __future__ import division

import os,sys,datetime
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict
import itertools
from array import *

from rootpy.io import root_open, DoesNotExist

import numpy as np
import pandas as pd

import math

import re

from rootpy.plotting import Hist2D, Hist

USAGE='write_obs_VM_format dbg[=False]'

#! user inputs
DBG=False
if len(sys.argv)>1: #! dbg entered by user
	if    sys.argv[1]=="True":  DBG=True
	elif  sys.argv[1]=="False": DBG=False
	else: sys.exit('DBG=%s is not valid. usage: %s'%(sys.argv[1],USAGE))

#! NOBS
NOBS=3
O1D,OIT,OR2=range(NOBS)
OBSFILE=["Obs_1D/Obs_1D.root","Obs_itg/Obs_itg.root","Obs_R2/Obs_R2.root"]

#! INDIR
INDIR='%s/thesis_obs_norm_ST_shape_032818'%os.environ['OBSDIR_E162'] #! 122817,121817,111317,080317
#! Use the following when output of make_obs_thesis.py is in some debugging area
#INDIR='%s/tmp/thesis_obs_norm_ST/dbg'%os.environ['OBSDIR_E163']
#INDIR='%s/tmp/thesis_obs_norm_ST'%os.environ['OBSDIR_E163']

#! OUTDIR
DATE=datetime.datetime.now().strftime('%m%d%y')
OUTDIR=[0 for i in range(NOBS)]
#! First create the root directory of OUTDIR
#! + For now create in INDIR, from where the data is read and converted to VM format
OUTDIRROOT=os.path.join(INDIR,"vm_format")#!os.environ['GDRIVE'],'tmp','obs_VM_format_%s'%DATE)
if DBG:
	OUTDIRROOT=os.path.join("/tmp/vm_format_%s"%DATE)#os.environ['GDRIVE'],'tmp','obs_VM_format_dbg_%s'%DATE)
if not os.path.exists(OUTDIRROOT):
	os.makedirs(OUTDIRROOT)
OUTDIR[O1D]="%s/1D"%OUTDIRROOT
OUTDIR[OIT]="%s/Itg"%OUTDIRROOT
OUTDIR[OR2]="%s/R2"%OUTDIRROOT

#! NQ2RANGES
NQ2RANGES=2
LQ2,HQ2=range(NQ2RANGES)
Q2RANGEDIR=["lowQ2","highQ2"]

#! FIN[obs][q2r]
FIN=[[None for j in range(NQ2RANGES)] for i in range(NOBS)]
#print FIN
for iobs,iq2r in list(itertools.product(range(NOBS),range(NQ2RANGES))):
	FIN[iobs][iq2r]=root_open("%s/%s/%s"%(INDIR,Q2RANGEDIR[iq2r],OBSFILE[iobs]))
if DBG:
	print FIN

def get_q2wbinlist(f,q2min=0.00,q2max=6.00,wmin=0.000,wmax=3.000,dbg=False,dbg_bins=2,
				   dbg_binl=['2.00-2.40_1.425-1.450','2.00-2.40_1.450-1.475','2.00-2.40_1.475-1.500','2.00-2.40_1.500-1.525']):
	"""
	+ Taken from disp_obs.py and modified accordingly
	+ Note in dbg mode, functions works as expected when 'dbg_bins'=number of bins in 'dbg_binl'
		+ 'dbg_binl' contains bins which need immediate analysis
	"""
	q2wbinl=[]
		
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
def get_wbin(q2wbin):
	return q2wbin.split('_')[1]

def get_q2bin_le_vm_format(q2wbin):
	val="%.1f"%float(get_q2bin(q2wbin).split('-')[0])
	#! remove decimal point
	val=val.replace(".","")
	return val
def get_wbin_le_vm_format(q2wbin):
	val="%.3f"%float(get_wbin(q2wbin).split('-')[0])
	#! remove decimal point
	val=val.replace(".","")
	return val
def get_q2bin_ue_vm_format(q2wbin):
	val="%.1f"%float(get_q2bin(q2wbin).split('-')[1])
	#! remove decimal point
	val=val.replace(".","")
	return val
def get_wbin_ue_vm_format(q2wbin):
	val="%.3f"%float(get_wbin(q2wbin).split('-')[1])
	#! remove decimal point
	val=val.replace(".","")
	return val

def write_1D():

	#! Note PAD_MAP extended to include explicit variable name (as per EI's syntax)
	PAD_MAP=[(1,1,"M1","mpippr"),      (2,3,'M2',"mpimpr"),      (3,2,'M2','mpippim'),
			 (4,1,"THETA","theta_pim"),(5,3,'THETA','theta_pip'),(6,2,'THETA','theta_pr'),
			 (7,1,"ALPHA",'psi_pim'),  (8,3,'ALPHA','psi_pip'),  (9,2,'ALPHA','psi_pr')]

	#! Get q2wbinl from file
	q2wbinl=[0 for i in range(NQ2RANGES)]
	if DBG:
		q2wbinl[LQ2]=get_q2wbinlist(FIN[O1D][LQ2],dbg=True,dbg_bins=1,dbg_binl=['2.00-2.40_1.500-1.525'])
		q2wbinl[HQ2]=get_q2wbinlist(FIN[O1D][HQ2],dbg=True,dbg_bins=1,dbg_binl=['3.00-3.50_1.500-1.525'])
	else:
		q2wbinl[LQ2]=get_q2wbinlist(FIN[O1D][LQ2])
		q2wbinl[HQ2]=get_q2wbinlist(FIN[O1D][HQ2])
	if DBG:
		print "q2wbinl[LQ2]=",q2wbinl[LQ2]
		print "q2wbinl[HQ2]=",q2wbinl[HQ2]

	#sys.exit()
	#! Loop over q2wbinl and write out observables
	for q2r in range(NQ2RANGES):
		for q2wbin in q2wbinl[q2r]:
			#! Get hobs/herr(seq,vst,var)
			#! + Note that, by design, herr['EC']=herr['EF']
			hobs=OrderedDict()
			herr=OrderedDict()
			for seq in ['EC','EF']:
				for item in PAD_MAP:
					pad,vst,var=item[0],item[1],item[2]
					hobs[seq,vst,var]=FIN[O1D][q2r].Get("%s/h1_%s_%d_%s"%(q2wbin,seq,vst,var))
					#! Note how herr['EC']=herr['EF']
					herr[seq,vst,var]=FIN[O1D][q2r].Get("%s/herr_h1_%s_%d_%s"%(q2wbin,'EF',vst,var))
			#print hobs
			#print herr

			#! Write out hobs/herr(seq,vst,var) to file
			#! Prepare outdir
			q2le,q2ue=get_q2bin_le_vm_format(q2wbin),get_q2bin_ue_vm_format(q2wbin)
			wle,wue  =get_wbin_le_vm_format(q2wbin), get_wbin_ue_vm_format(q2wbin)
			#print q2le,q2ue
			#print wle,wue
			q2id="q2_%s_%s"%(q2le,q2ue)
			wid ="w_%s_%s"%(wle,wue)
			outdir="%s/%s/%s"%(OUTDIR[O1D],q2id,wid)
			if not os.path.exists(outdir):
				os.makedirs(outdir)
			#! Now write result
			#! + Currently on for EF
			for seq in ['EC','EF']:
				if seq=='EC':continue
				for item in PAD_MAP:
					pad,vst,var,varname=item[0],item[1],item[2],item[3]
					ftxt=open("%s/xsec1d_%s.dat"%(outdir,varname),"w")
					nbins=hobs[seq,vst,var].GetNbinsX()
					for ibin in range(nbins):
						binc=hobs[seq,vst,var].GetBinCenter(ibin+1)
						xsec=hobs[seq,vst,var].GetBinContent(ibin+1)
						staterr=hobs[seq,vst,var].GetBinError(ibin+1)
						systerr=herr[seq,vst,var].GetBinContent(ibin+1)
						ftxt.write("%f %f %f %f\n"%(binc,xsec,staterr,systerr))
						#! Update ERRS
						if xsec!=0: #! To avoid div by 0, which there are very few cases of
							rel_syst_err=(systerr/xsec)*100
							rel_stat_err=(staterr/xsec)*100
							rel_tot_err=math.sqrt(rel_syst_err**2+rel_stat_err**2)
							ERRS["syst",O1D].append(rel_syst_err)
							ERRS["stat",O1D].append(rel_stat_err)
							ERRS["total",O1D].append(rel_tot_err)
					ftxt.close()

def check_1D():
	NFILES_TOTAL=5*29*9 #!nQ2*nW*n1D=1305
	nfiles_total=0
	for path,dirs,files in os.walk(OUTDIR[O1D]):
		#! Consider only files in q2 and w directory structure
		if re.search("q2_{1}.*w_{1}",path)==None: continue
		#print path
		nfiles=len(files)
		nfiles_total+=nfiles
		if nfiles!=9: print "Files missingt in %s"%path
	if DBG==True: return #! then ignore following because not all data is written out
	if nfiles_total!=NFILES_TOTAL:
		print "Some files are missing. #written-files!=Expected #written-files (%d!=%d)"%(nfiles_total,NFILES_TOTAL)
	else:
		print "All 1D data written. #written-files=Expected #written-files (%d=%d)"%(nfiles_total,NFILES_TOTAL)

def get_bin_center_obs_itg(hobs,ibin):
	label=hobs.GetXaxis().GetBinLabel(ibin+1)
	le=float(label.split(',')[0].replace('[',''))
	ue=float(label.split(',')[1].replace(')',''))
	bc=(le+ue)/2
	return bc	

def write_itg():
	#! Store Q2[q2r] bin names for AT and EI in one structure
	Q2L=[0 for i in range(NQ2RANGES)]
	Q2L[LQ2]=[["2.00-2.40_1.400-2.125","cs20_24"],["2.40-3.00_1.400-2.125","cs24_30"]]
	Q2L[HQ2]=[["3.00-3.50_1.400-2.125","cs30_35"],["3.50-4.20_1.400-2.125","cs35_42"],["4.20-5.00_1.400-2.125","cs42_50"]]

	for q2r in range(NQ2RANGES):
		for iq2,q2l in enumerate(Q2L[q2r]):
			q2b=q2l[0]
			q2n=q2l[1]

			#! get histograms
			hobs=FIN[OIT][q2r].Get('%s/hW_obs_EF_THETA'%q2b)
			herr=FIN[OIT][q2r].Get('%s/hW_err_EF_THETA'%q2b)
			#! Write data from histograms
			outdir="%s"%OUTDIR[OIT]
			if not os.path.exists(outdir):
				os.makedirs(outdir)
			ftxt=open("%s/%s.txt"%(outdir,q2n),"w")
			nbins=hobs.GetNbinsX()
			for ibin in range(nbins):
				binc=get_bin_center_obs_itg(hobs,ibin)
				xsec=hobs.GetBinContent(ibin+1)
				staterr=hobs.GetBinError(ibin+1)
				systerr=herr.GetBinContent(ibin+1)
				ftxt.write("%f %f %f %f\n"%(binc,xsec,staterr,systerr))
				#! Update ERRS
				if xsec!=0: #! To avoid div by 0, which there are very few cases of
					rel_syst_err=(systerr/xsec)*100
					rel_stat_err=(staterr/xsec)*100
					rel_tot_err=math.sqrt(rel_syst_err**2+rel_stat_err**2)
					ERRS["syst",OIT].append(rel_syst_err)
					ERRS["stat",OIT].append(rel_stat_err)
					ERRS["total",OIT].append(rel_tot_err)
			ftxt.close()	

def check_itg():
	NFILES_TOTAL=5*1 #!nQ2*nItg=5
	nfiles_total=0
	for path,dirs,files in os.walk(OUTDIR[OIT]):
		nfiles=len(files)
		nfiles_total+=nfiles
		if nfiles!=5: print "Files missingt in %s"%path
	if DBG==True: return #! then ignore following because not all data is written out
	if nfiles_total!=NFILES_TOTAL:
		print "Some files are missing. #written-files!=Expected #written-files (%d!=%d)"%(nfiles_total,NFILES_TOTAL)
	else:
		print "All Itg data written. #written-files=Expected #written-files (%d=%d)"%(nfiles_total,NFILES_TOTAL)

def write_R2():
	#! Following dictionaries for extracting R2
	R2S=['A','B','C','D','E']
	#! Following modified from its original in disp_obs.py for easier reading in text file
	R2_NAMED={'A':'R2-T_plus_R2-L','B':'R2-c-LT','C':'R2-c-TT','D':'R2-s-LT','E':'R2-s-TT'}

	#! Note PAD_MAP extended to include explicit variable name (as per EI's syntax)
	PAD_MAP=[(1,1, 'M1',   'mpippr_phi_pim'),    (2,3, 'M1',   'mpippr_phi_pip'),    (3,2, 'M1',   'mpippr_phi_pr'), 
			 (4,1, 'M2',   'mpippim_phi_pim'),   (5,3, 'M2',   'mpimpr_phi_pip'),    (6,2, 'M2',   'mpippim_phi_pr'),
			 (7,1, 'THETA','theta_pim_phi_pim'), (8,3, 'THETA','theta_pip_phi_pip'), (9,2, 'THETA','theta_pr_phi_pr'),
			 (10,1,'ALPHA','psi_pim_phi_pim'),   (11,3,'ALPHA','psi_pip_phi_pip'),   (12,2,'ALPHA','psi_pr_phi_pr')]
	
	#! Get q2wbinl from file
	q2wbinl=[0 for i in range(NQ2RANGES)]
	if DBG:
		q2wbinl[LQ2]=get_q2wbinlist(FIN[OR2][LQ2],dbg=True,dbg_bins=1,dbg_binl=['2.00-2.40_1.500-1.525'])
		q2wbinl[HQ2]=get_q2wbinlist(FIN[OR2][HQ2],dbg=True,dbg_bins=1,dbg_binl=['3.00-3.50_1.500-1.525'])
	else:
		q2wbinl[LQ2]=get_q2wbinlist(FIN[OR2][LQ2])
		q2wbinl[HQ2]=get_q2wbinlist(FIN[OR2][HQ2])
	if DBG:
		print "q2wbinl[LQ2]=",q2wbinl[LQ2]
		print "q2wbinl[HQ2]=",q2wbinl[HQ2]

	#sys.exit()
	#! Loop over q2wbinl and write out observables
	for q2r in range(NQ2RANGES):
		for q2wbin in q2wbinl[q2r]:
			#! Get hobs/herr(R2,seq,vst,var)
			#! + Note that, by design, herr['EC']=herr['EF']
			hobs=OrderedDict()
			herr=OrderedDict()
			for R2 in R2S:
				for seq in ['EC','EF']:
					for item in PAD_MAP:
						pad,vst,var=item[0],item[1],item[2]
						if (R2=='D' or R2=='E') and var!= 'ALPHA':continue
						hobs[R2,seq,vst,var]=FIN[OR2][q2r].Get("%s/hR2_%s_%s_%d_%s"%(q2wbin,R2,seq,vst,var))
						#! Note how herr['EC']=herr['EF']
						herr[R2,seq,vst,var]=FIN[OR2][q2r].Get("%s/herr_hR2_%s_%s_%d_%s"%(q2wbin,R2,'EF',vst,var))
			#print hobs
			#print herr

			#! Write out hobs/herr(seq,vst,var) to file
			#! Prepare outdir
			q2le,q2ue=get_q2bin_le_vm_format(q2wbin),get_q2bin_ue_vm_format(q2wbin)
			wle,wue  =get_wbin_le_vm_format(q2wbin), get_wbin_ue_vm_format(q2wbin)
			#print q2le,q2ue
			#print wle,wue
			q2id="q2_%s_%s"%(q2le,q2ue)
			wid ="w_%s_%s"%(wle,wue)
			#! Now write result
			#! + Currently on for EF
			for R2 in R2S:
				outdir="%s/%s/%s/%s"%(OUTDIR[OR2],R2_NAMED[R2],q2id,wid)
				if not os.path.exists(outdir):
					os.makedirs(outdir)
				for seq in ['EC','EF']:
					if seq=='EC':continue
					for item in PAD_MAP:
						pad,vst,var,varname=item[0],item[1],item[2],item[3]
						if (R2=='D' or R2=='E') and var!= 'ALPHA':continue
						ftxt=open("%s/xsec%s_%s.dat"%(outdir,R2_NAMED[R2],varname),"w")
						nbins=hobs[R2,seq,vst,var].GetNbinsX()
						for ibin in range(nbins):
							binc=hobs[R2,seq,vst,var].GetBinCenter(ibin+1)
							xsec=hobs[R2,seq,vst,var].GetBinContent(ibin+1)
							staterr=hobs[R2,seq,vst,var].GetBinError(ibin+1)
							systerr=herr[R2,seq,vst,var].GetBinContent(ibin+1)
							ftxt.write("%f %f %f %f\n"%(binc,xsec,staterr,systerr))
							#! Update ERRS
							if xsec!=0: #! To avoid div by 0, which there are very few cases of
								rel_syst_err=(systerr/xsec)*100
								rel_stat_err=(staterr/xsec)*100
								rel_tot_err=math.sqrt(rel_syst_err**2+rel_stat_err**2)
								ERRS["syst",OR2].append(rel_syst_err)
								ERRS["stat",OR2].append(rel_stat_err)
								ERRS["total",OR2].append(rel_tot_err)
						ftxt.close()	
def check_R2():
	NFILES_TOTAL=5*29*42 #!nQ2*nW*nR2=6090
	nfiles_total=0
	for path,dirs,files in os.walk(OUTDIR[OR2]):
		#! Consider only files in q2 and w directory structure
		if re.search("q2_{1}.*w_{1}",path)==None: continue
		#print path
		nfiles=len(files)
		nfiles_total+=nfiles
		if "-c-" in path:
			if nfiles!=12: print "Files missingt in %s"%path
		elif "-s-" in path:#! applicable only for alpha
			if nfiles!=3: print "Files missingt in %s"%path
	if DBG==True: return #! then ignore following because not all data is written out
	if nfiles_total!=NFILES_TOTAL:
		print "Some files are missing. #written-files!=Expected #written-files (%d!=%d)"%(nfiles_total,NFILES_TOTAL)
	else:
		print "All R2 data written. #written-files=Expected #written-files (%d=%d)"%(nfiles_total,NFILES_TOTAL)

#! main

#! Create some structures to hold {syst-err,stat-err,tot-err}*{1D,itg,R2}
#! + These are simple structures in that they hold a list of errors over all data points: q2,w,vst,var,bin etc
#! + Since currently on EF data is written out, these structures hold data only for EF
#! + These structures will be filled in their respective write functions
#! + The values stored are relative errors (%)
ERRNAMES=["syst","stat","total"]
ERRS=OrderedDict()
for err in ERRNAMES:
	for obs in [O1D,OIT,OR2]:
		ERRS[err,obs]=[]
		

#! Now start writing data and filling error structures
write_1D()
check_1D()
write_itg()
check_itg()
write_R2()
check_R2()

#! Plot errors
#! .png/.pdf output
OUTDIR_SE="%s/syst-errs"%OUTDIRROOT
if not os.path.exists(OUTDIR_SE):
	os.makedirs(OUTDIR_SE)
#! .root output
fseroot=ROOT.TFile("%s/syst-errs.root"%OUTDIR_SE,"RECREATE")
#! ROOT plot aesthetics
ROOT.gStyle.SetOptStat("nemruo")
for err in ERRNAMES:
	for item in zip([O1D,OIT,OR2],["1D","itg","R2"]):
		obs=item[0]
		obsname=item[1]
		#! SYST_ERR
		h=Hist(200,0,200)
		name="%s-%s"%(err,obsname)
		h.SetNameTitle(name,name)
		h.fill_array(ERRS[err,obs])
		c=ROOT.TCanvas()
		h.Draw()
		c.SaveAs("%s/%s.png"%(OUTDIR_SE,name))
		c.SaveAs("%s/%s.pdf"%(OUTDIR_SE,name))
		h.Write()
fseroot.Close()

print "If the progam is not terminating, then Python is probably doing \"garbage collection\"(?); Wait a while!"