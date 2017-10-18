#!/usr/bin/python
from __future__ import division
import os,sys,datetime
from collections import OrderedDict
import array

import ROOT

from rootpy.io import root_open, DoesNotExist

'''
Read documentation for function plot() to see what this code does.
'''

USAGE='study_yields_acceptance dbg[=False] ana_ER_commonbins[=True] show_rel_err_dist[=True] SSBands=[off-off]'

#! user inputs
DBG=False
if len(sys.argv)>1: #! dbg entered by user
	if    sys.argv[1]=="True":  DBG=True
	elif  sys.argv[1]=="False": DBG=False
	else: sys.exit('DBG=%s is not valid. usage: %s'%(sys.argv[1],USAGE))

ANA_ER_COMMONBINS=True
if len(sys.argv)>2: #!  ana_ER_commonbins entered by user
	if    sys.argv[2]=="True":  ANA_ER_COMMONBINS=True
	elif  sys.argv[2]=="False": ANA_ER_COMMONBINS=False
	else: sys.exit('ANA_ER_COMMONBINS=%s is not valid. usage: %s'%(sys.argv[2],USAGE))

SHOW_REL_ERR_DIST=True
if len(sys.argv)>3: #!  show_rel_err_dist entered by user
	if    sys.argv[3]=="True":  SHOW_REL_ERR_DIST=True
	elif  sys.argv[3]=="False": SHOW_REL_ERR_DIST=False
	else: sys.exit('SHOW_REL_ERR_DIST=%s is not valid. usage: %s'%(sys.argv[3],USAGE))

SSBANDS='off-off'
if len(sys.argv)>4: #!  SSBands entered by user
	if sys.argv[4]=="on-on" or sys.argv[4]=="off-off":
	   SSBANDS=sys.argv[4]
	else:
		sys.exit('SSBANDS=%s is not valid. usage: %s'%(sys.argv[4],USAGE))

print "DBG=",DBG
print "ANA_ER_COMMONBINS=",ANA_ER_COMMONBINS
print "SHOW_REL_ERR_DIST=",SHOW_REL_ERR_DIST
print "SSBANDS=",SSBANDS
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
ROOT.gStyle.SetOptStat("nmMrReiuo")
ROOT.gStyle.SetStatStyle(0) #! transparent stats box

#! OUTDIR
DATE=datetime.datetime.now().strftime('%m%d%y')
OUTDIRNAME='results'
#! append identifiers
if ANA_ER_COMMONBINS==True:
	OUTDIRNAME+='_ER-5D-CommonBins'
else:
	OUTDIRNAME+='_All-5D-Bins'
if SHOW_REL_ERR_DIST==True:
	OUTDIRNAME+='_with_rel_err_dist'
else:
	OUTDIRNAME+='_without_rel_err_dist'
if SSBANDS=='on-on':
	OUTDIRNAME+='_SSBands-on-on'
elif SSBANDS=='off-off':
	OUTDIRNAME+='_SSBands-off-off'
OUTDIRNAME+='_%s'%DATE
#! Finally create OUTDIR
if DBG==True:
	OUTDIR=os.path.join(os.environ['STUDY_YIELDS_ACCEPTANCE_DATADIR'],'dbg',OUTDIRNAME)
else:
	OUTDIR=os.path.join(os.environ['STUDY_YIELDS_ACCEPTANCE_DATADIR'],OUTDIRNAME)
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)
print "OUTDIR=",OUTDIR
#sys.exit()

#! Setup input root files
#! + As per SSBANDS: on-on=cutsncors1, off-off=cutsncors2
if   SSBANDS=='off-off': CUTSNCORS='cutsncors1'
elif SSBANDS=='on-on':   CUTSNCORS='cutsncors2'
#! + Different file for Q2 ranges (lowQ2,highQ2)
NQ2RANGES=2
LQ2,HQ2=range(NQ2RANGES)
FIN=[0 for i in range(NQ2RANGES)]
FIN[LQ2]=root_open('%s/lowQ2_SSBands_080217/%s/sim4_sim5_sim6_sim7_sim8_sim13/yield.root'%(os.environ['OBSDIR_E16'],CUTSNCORS),'r')
FIN[HQ2]=root_open('%s/highQ2_SSBands_080217/%s/sim9_sim10_sim11_sim12/yield.root'%(os.environ['OBSDIR_E16'],CUTSNCORS),'r')


def get_q2wbinlist(q2min=0.00,q2max=6.00,wmin=0.000,wmax=3.000,dbg=False,dbg_bins=4,
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

		brk=False #! technical tool to break out of two nested for loops. Set in second (nested) for loop
		for f in [FIN[LQ2],FIN[HQ2]]:
			if brk==True: break
			for path,dirs,files in f.walk():
				if path=="":continue #! Avoid root path
				path_arr=path.split("/")
				if len(path_arr)==1:
					if DBG==True:
						if path in dbg_binl:
							q2wbinl.append(path)
					else:
						q2wbinl.append(path)
				if DBG==True and len(q2wbinl)==dbg_bins:
					brk=True
					break #! Uncomment/comment -> Get limited q2w-bins/Get all q2w-bins



		# #! Remove q2wbins that are not within [q2min,q2max],[wmin,wmax] 
		# q2wbins_remove=[]
		# for q2wbin in q2wbinl:
		# 	q2bin_le=q2wbin.split("_")[0].split("-")[0]
		# 	q2bin_ue=q2wbin.split("_")[0].split("-")[1]
		# 	wbin_le =q2wbin.split("_")[1].split("-")[0]
		# 	wbin_ue =q2wbin.split("_")[1].split("-")[1]
		# 	if float(q2bin_ue)<=q2min or float(q2bin_le)>=q2max or float(wbin_ue)<=wmin or float(wbin_le)>=wmax:
		# 		q2wbins_remove.append(q2wbin)
		# for q2wbin in q2wbins_remove:
		# 	q2wbinl.remove(q2wbin)

		return q2wbinl


def plot(q2wbin):	
	'''
	+ For {ER,ST,SR,SA} plot:
	    i.  5D-PS distributions of entries and their corresponding rel_err
	    ii. 1D projection (currently only for VST1,THETA i.e. theta_pim)
	+ Additionally for {ST,SR,SA}
	    iii.  5D-PS-vst-var distributions of entries and their corresponding rel_err
	    	  + 5D-PS-vst-var= 5D PS distributions but only for a particular subspace of 
	    	    the total 5D-PS, which corresponds to a particular VST-VAR's bin.
	          + Currently this is done only for VST1,THETA i.e. theta_pim
	'''
	print "processing q2wbin=",q2wbin

	var=VAR_NAMES_PLAIN[1,'THETA']

	SEQ=['ER','ST','SR','SA']

	#! Get all the h5(seq)
	#! + First establish if to use FIN[LQ2] or FIN[HQ2] for this q2wbin
	if   FIN[LQ2].FindObjectAny(q2wbin)!=None: fin=FIN[LQ2]
	elif FIN[HQ2].FindObjectAny(q2wbin)!=None: fin=FIN[HQ2]
	else:
		sys.exit("study_yields_acceptance.py:plot():Could not establish if to use FIN[LQ2] or FIN[HQ2] for q2wbin=%s"(q2wbin))
	print "fin=",fin.GetName()
	
	h5=OrderedDict()
	for seq in SEQ:
		h5[seq]=fin.Get('%s/%s/VST1/h5'%(q2wbin,seq))
	
	#! outdir_q2w
	outdir_q2w="%s/%s"%(OUTDIR,q2wbin)
	if not os.path.exists(outdir_q2w):
		os.makedirs(outdir_q2w)

	#! Now start processing
	for seq in SEQ:
		print "processing seq=",seq

		#! i. 5D-PS distributions
		#! BinContentDist
		outdir="%s/%s/%s"%(outdir_q2w,seq,'5D-PS')
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		if seq=='ER': #! Note ANA_ER_COMMONBINS not applicable for ER
			h1tot=thntool.GetBinContentDist(h5[seq],5,0.5,5.5)
		if seq=='SA':
			if ANA_ER_COMMONBINS: h1tot=thntool.GetBinContentDistCommonBins(h5[seq],h5['ER'],1000,0,0.5)
			else:       		  h1tot=thntool.GetBinContentDist(h5[seq],1000,0,0.5)
			
		elif seq=='SR':
			if ANA_ER_COMMONBINS: h1tot=thntool.GetBinContentDistCommonBins(h5[seq],h5['ER'],100,0,500)
			else:                 h1tot=thntool.GetBinContentDist(h5[seq],100,0,500)
			
		elif seq=='ST':
			if ANA_ER_COMMONBINS: h1tot=thntool.GetBinContentDistCommonBins(h5[seq],h5['ER'])
			else:                 h1tot=thntool.GetBinContentDist(h5[seq])
			
		#! BinRelErrDist
		if seq=='ER': #! Note ANA_ER_COMMONBINS not applicable for ER
			h1errtot=thntool.GetBinRelErrorDist(h5[seq],100,0,1.5)
		else:
			if ANA_ER_COMMONBINS: h1errtot=thntool.GetBinRelErrorDistCommonBins(h5[seq],h5['ER'],100,0,1.5)
			else:                 h1errtot=thntool.GetBinRelErrorDist(h5[seq],100,0,1.5)
		
		#! Some hist aesthetics
		h1tot.SetTitle("h5_BinContentDist_%s_%s"%(seq,q2wbin))
		h1errtot.SetTitle("h5_BinRelErrDist_%s_%s"%(seq,q2wbin))
		
		#! Draw and save
		c=ROOT.TCanvas()
		if SHOW_REL_ERR_DIST:
			c.Divide(1,2)
			c.cd(1)
			h1tot.Draw("e")
			c.cd(2)
			h1errtot.Draw("e")
		else:
			h1tot.Draw("e")
		c.SaveAs("%s/h1tot.png"%outdir)

		#! ii. 1D projection (currently only for VST1,THETA i.e. theta_pim)
		outdir="%s/%s/%s"%(outdir_q2w,seq,'1D-proj-%s'%var)
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		h1d=h5[seq].Projection(H5_DIM['THETA'],"E")
		h1d.SetTitle("1D_proj_%s_%s_%s"%(seq,q2wbin,var))
		c=ROOT.TCanvas()
		h1d.Draw("e")
		c.SaveAs("%s/h1d.png"%outdir)
		#sys.exit()

		#! iii. 5D-PS-vst-var distributions
		#! + only for {'ST','SR','SA'}
		if seq=='ER': continue
		outdir_var="%s/%s/%s"%(outdir_q2w,seq,'5D-PS-%s'%var)
		if not os.path.exists(outdir_var):
			os.makedirs(outdir_var)
		#! Get VST1-THETA binning information
		nbins=h5[seq].GetAxis(H5_DIM['THETA']).GetNbins()
		#print nbins
		binw=h5[seq].GetAxis(H5_DIM['THETA']).GetBinWidth(1)
		for ibin in range(nbins):
			#! Set range and make projections as per bin
			binle=h5[seq].GetAxis(H5_DIM['THETA']).GetBinLowEdge(ibin+1)
			binue=binle+binw
			print "processing bin number: %d [%.3f,%.3f)"%(ibin+1,binle,binue)
			#! Set range
			h5['ER'].GetAxis(H5_DIM['THETA']).SetRange(ibin+1,ibin+1)
			h5[seq].GetAxis(H5_DIM['THETA']).SetRange(ibin+1,ibin+1)
			#! h5->h4 for selected range
			h4ER=h5['ER'].Projection(4,H4_PROJDIM,"E")
			h4=h5[seq].Projection(4,H4_PROJDIM,"E")
			#! get h1
			if seq=='SA':
				if ANA_ER_COMMONBINS: h1=thntool.GetBinContentDistCommonBins(h4,h4ER,1000,0,0.5)
				else:                 h1=thntool.GetBinContentDist(h4,1000,0,0.5)
			elif seq=='SR':
				if ANA_ER_COMMONBINS: h1=thntool.GetBinContentDistCommonBins(h4,h4ER,100,0,200)
				else:                 h1=thntool.GetBinContentDist(h4,100,0,200)
			else:
				if ANA_ER_COMMONBINS: h1=thntool.GetBinContentDistCommonBins(h4,h4ER)
				else:                 h1=thntool.GetBinContentDist(h4)
			
			#! get h1err
			if ANA_ER_COMMONBINS: h1err=thntool.GetBinRelErrorDistCommonBins(h4,h4ER,100,0,1)
			else:                 h1err=thntool.GetBinRelErrorDist(h4)
			
			#! Some hist aesthetics
			h1.SetTitle("h5_BinContentDist_%s_%s_%s_bin%d"%(seq,q2wbin,var,ibin+1))
			h1err.SetTitle("h5_BinRelErrDist_%s_%s_%s_bin%d"%(seq,q2wbin,var,ibin+1))
			
			#! Draw and save
			c=ROOT.TCanvas()
			if SHOW_REL_ERR_DIST:
				c.Divide(1,2)
				c.cd(1)
				h1.Draw("e")
				c.cd(2)
				h1err.Draw("e")
			else:
				h1.Draw("e")
			c.SaveAs("%s/h1_%02d_%04.1f-%04.1f.png"%(outdir_var,ibin+1,binle,binue))

			#! Reset Range
			h5['ER'].GetAxis(H5_DIM['THETA']).SetRange()
			h5[seq].GetAxis(H5_DIM['THETA']).SetRange()

#! main
q2wbinl=get_q2wbinlist()
print "q2wbinl=",q2wbinl
for q2wbin in q2wbinl:
	plot(q2wbin)