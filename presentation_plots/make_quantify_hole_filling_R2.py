#!/usr/bin/python
from __future__ import division

import os,sys,subprocess
import ROOT
from collections import OrderedDict
import math

'''
+ Note that this script makes plots pretty for ananote using data already produced by
  'make_obs_thesis.py'
+ In this script the effect of hole-filling is demonstrated with the aid 
  of Obs_R2 for Q2WBINS=['2.40-3.00_1.725-1.750','4.20-5.00_2.100-2.125']

+ syst-err for hole-filling at the 1D level is also plotted
'''
#! Follwing method taken from make_obs_thesis.py
#! + Minor modifications made for this code, for eample argument 'obs' removed
def get_herr_HF_absolute(hEC,hEF,vst,var):
	'''
	In each bin, calculate absolute SE due to Hole-Filling:
		+ se=abs(binc_EF-bincEC)/2
		+ se_err=0 (for now at least)
	+ Note there are 8 unique cases that the above formula needs to be aware of: EC(0,<0,>0)* EF(0,<0,>0):
		+ EC=0 & :
			+ 1. EF=0: se=formula(=0)
			+ 2. EF>0: se=formula      (EXPECTED)
			+ 3. EF<0: se=0            (SHOULD NOT OCCUR, but if so set se=0 for now, and investigate why)
		+ EF=0 & :
			+    EC=0: se=formula      (same as 1.)
			+ 4. EC<0: se=0            (SHOULD NOT OCCUR, but if so set se=0 for now, and investigate why)
			+ 5. EC>0: se=0            (SHOULD NOT OCCUR, but if so set se=0 for now, and investigate why)
		+ EC>0 &:
			+ 6. EF>0: se=formula      (EXPECTED case, applicable for EF>EC and EF<EC, the latter should not occur, but even if so, obtain se from formula)
			+ 7. EF<0: se=0            (SHOULD NOT OCCUR, but if so set se=0 for now, and investigate why)
		+ EF>0:
			+    EC>0: se=formula      (same as 6.)
			+    EC<0: se=0            (SHOULD NOT OCCUR, but if so set se=0 for now, and investigate why)
	'''
	#! First make herr_HF as a clone of hEC and resetting its bin contents
	herr=hEC.Clone("herr_abs_%d_%s"%(vst,var))
	herr.Reset()
	#! Now obtain SE due to Hole-Filling in each bin
	nbins=hEC.GetNbinsX()
	for ibin in range(nbins):
		binc_EF=hEF.GetBinContent(ibin+1)
		binc_EC=hEC.GetBinContent(ibin+1)
		#! Obtain se, se_err as per 8 unique cases of formula (see documentation)
		#! + Only for cases 1, 2 and 6 se has non-zero value, else it is 0
		if binc_EC==0 and binc_EF==0: #! case 1.
			se=math.fabs(binc_EF-binc_EC)/2
			se_err=0
		elif binc_EC==0 and binc_EF>0: #! case 2.
			se=math.fabs(binc_EF-binc_EC)/2
			se_err=0
		elif binc_EC>0 and binc_EF>0: #! case 6. for EF>EC and EF<EC (latter should not occur, but even if so)
			se=math.fabs(binc_EF-binc_EC)/2
			se_err=0
		else: #! for all other cases
			se=0
			se_err=0
		#! Set se, se_err in herr's bin
		herr.SetBinContent(ibin+1,se)
		herr.SetBinError(ibin+1,se_err)
	return herr
#! OUTDIR
OUTDIR="%s/figures/Holes/"%os.environ['ANANOTE']
#OUTDIR="/tmp"

#! Define some objects that will be needed in steps 1 and step 2
#! PAD_MAP[obs]=(pad,vst,var)
PAD_MAP=OrderedDict()
PAD_MAP['R2']=[(1,1,"M1"),    (2,3,'M1'),    (3,2,'M1'), 
			   (4,1,"M2"),    (5,3,'M2'),    (6,2,'M2'),
			   (7,1,"THETA"), (8,3,'THETA'), (9,2,'THETA'),
			   (10,1,"ALPHA"),(11,3,'ALPHA'),(12,2,'ALPHA')]

R2L=['A','B','C','D','E']
R2_NAMED={'A':'R2_{T}+R2_{L}','B':'R2^{c}_{LT}','C':'R2^{c}_{TT}','D':'R2^{s}_{LT}','E':'R2^{s}_{TT}'}

fobs_lQ2    =ROOT.TFile("%s/thesis_obs_norm_ST_shape_122817/lowQ2/Obs_R2/Obs_R2.root"%os.environ['OBSDIR_E162'],"READ")
fobs_hQ2    =ROOT.TFile("%s/thesis_obs_norm_ST_shape_122817/highQ2/Obs_R2/Obs_R2.root"%os.environ['OBSDIR_E162'],"READ")
fsyserrs=ROOT.TFile("%s/thesis_obs_norm_ST_shape_122817/syst-errs/syst-errs.root"%os.environ['OBSDIR_E162'],"READ")

#! 1. First show a direct plot of hobs for EC,EF and herr
hobs=OrderedDict()
#! Following contains absolute syst-err due to hole filling
herr=OrderedDict()
SEQS=['EC','EF']
Q2WBINS=['2.40-3.00_1.725-1.750','4.20-5.00_2.100-2.125'] #'4.20-5.00_1.725-1.750'
for q2wbin in Q2WBINS:#'2.40-3.00_1.725-1.750'
	if   q2wbin=='2.40-3.00_1.725-1.750': fobs=fobs_lQ2
	elif q2wbin=='4.20-5.00_2.100-2.125': fobs=fobs_hQ2
	for R2 in R2L:
		for pm in PAD_MAP['R2']:
			vst,var=pm[1],pm[2]
			if (R2=='D' or R2=='E') and var != 'ALPHA': continue
			#! get hobs and herr from .root file
			hobs['EC',vst,var]=fobs.Get("%s/hR2_%s_%s_%d_%s"%(q2wbin,R2,'EC',vst,var))
			hobs['EF',vst,var]=fobs.Get("%s/hR2_%s_%s_%d_%s"%(q2wbin,R2,'EF',vst,var))
			#! get herr (unfortunately, it has to be recalculated because in make_obs_thesis.py, it is added
			#! to the rest of the syst-errs)
			herr[vst,var]=get_herr_HF_absolute(hobs['EC',vst,var],hobs['EF',vst,var],vst,var)
			herr[vst,var].SetFillColor(ROOT.gROOT.ProcessLine("kBlack"))
			herr[vst,var].SetFillStyle(3003)
			herr[vst,var].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
			herr[vst,var].SetMarkerStyle(ROOT.gROOT.ProcessLine("kDot"))
			herr[vst,var].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
		
		#! plot
		ROOT.gStyle.SetOptStat(0)
		#! legend
		l=ROOT.TLegend(0.50,0.70,0.90,0.90)
		l.SetFillStyle(0)
		#l.SetBorderSize(0)
		l.SetTextSize(0.06)
		#! Canvas
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
		pad_p.Divide(3,4)
		for item in PAD_MAP['R2']:
			pad,vst,var=item[0],item[1],item[2]
			if (R2=='D' or R2=='E') and var != 'ALPHA': continue
			#print "pad,vst,var=",pad,vst,var
			gpad=pad_p.cd(pad)
			#! Draw
			hobs['EC',vst,var].Draw()
			hobs['EF',vst,var].Draw("same")
			herr[vst,var].Draw("e hist same")
			#! legend
			if (R2=='D' or R2=='E') and pad==10:
				l.AddEntry(hobs['EC',vst,var],"EC","p")
				l.AddEntry(hobs['EF',vst,var],"EF","p")
				l.AddEntry(herr[vst,var],"syst-err-HF","p")
				l.Draw()
			elif pad==1:
				l.AddEntry(hobs['EC',vst,var],"EC","p")
				l.AddEntry(hobs['EF',vst,var],"EF","p")
				l.AddEntry(herr[vst,var],"syst-err-HF","p")
				l.Draw()
		c.SaveAs("%s/HF_EC_EF_R2_%s_%s.png"%(OUTDIR,q2wbin,R2))
		#! Convert .png->.eps->.pdf since in .eps and .pdf produced directly by ROOT
		#! the y-axis titles are cutoff
		# #! 1. .png->.eps
		# cmd=['convert',"%s/c_q%s.png"%(outdir,q2wbin),"%s/c_q%s.eps"%(outdir,q2wbin)]
		# tool=subprocess.Popen(cmd,stderr=subprocess.STDOUT)
		# tool.wait()
		#! 2. .png->.pdf
		cmd=['convert',"%s/HF_EC_EF_R2_%s_%s.png"%(OUTDIR,q2wbin,R2),"%s/HF_EC_EF_R2_%s_%s.pdf"%(OUTDIR,q2wbin,R2)]
		tool=subprocess.Popen(cmd,stderr=subprocess.STDOUT)
		tool.wait()

#! Now make plot of the relative systematic error
hrelse=fsyserrs.Get("SYST_ERR_HF_REL_R2")
hrelse.SetTitle("#splitline{      Distribution of Relative Systematic Uncertainty (%) to}{Photon Polarization Dependent Cross-Sections from Hole Filling}")
#! Not sure why, but the following line is needed to get stats to show after ROOT.gStyle.SetOptStat(0)
#! + I tried ROOT.gStyle.Reset() and hrelse.SetStats(1), but neither seem to work
hrelse.SetName("")
#ROOT.gStyle.Reset()
#hrelse.SetStats(1)
ROOT.gStyle.SetOptStat("mMrRuo")#"emMrRuo"
# ROOT.gStyle.SetTitleW(8.5)# //title width 
# ROOT.gStyle.SetTitleFontSize(20)# 
# ROOT.gStyle.SetTitleH(0.05)
ROOT.gStyle.SetTitleX(.44)
ROOT.gStyle.SetTitleY(1.0)
ROOT.gStyle.SetTitleSize(0.04,"other")
c1=ROOT.TCanvas("c1","c1",800,500)
# pad_p=ROOT.TPad("pad_p","Plots pad",0.00,0.00,1.00,1.00)
# pad_p.Draw()
# pad_p.cd()
hrelse.GetXaxis().SetRangeUser(0,60)
#c1.SetLogy()
hrelse.SetXTitle("Relative Systematic Uncertainty (%)")
hrelse.Draw("hist")
c1.Update()
c1.SaveAs("%s/rel_syst_err_HF_R2.png"%(OUTDIR))
c1.SaveAs("%s/rel_syst_err_HF_R2.pdf"%(OUTDIR))
# cmd=['convert',"%s/rel_syst_err_HF_1D.png"%(OUTDIR),"%s/rel_syst_err_HF_1D.pdf"%(OUTDIR)]
# tool=subprocess.Popen(cmd,stderr=subprocess.STDOUT)
# tool.wait()


