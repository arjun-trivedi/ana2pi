#!/usr/bin/python
from __future__ import division

import os,sys
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict
import itertools

import time,datetime

import numpy as np

import math

import atlib

'''
+ This script simply consolidates plots, listed below, already produced by 'study_det_ineff.py' in .root files in
 '$STUDY_DET_INEFF_E16/results_top2' and saves them directly as .pdf in '$STUDY_DET_INEFF_E16/results_top2/thesis_plots'
  that can be copied over to the thesis folder. Plots consolidated:
  + "pre-cut" plots in '$STUDY_DET_INEFF_E16/results_top2/kinematics.root'
  + "pst-cut" plots in '$STUDY_DET_INEFF_E16/results_top2/sc_pd_theta_vs_p/bad_good_sc_pd_theta_vs_p.root"
NOTE that "cut" here refers to kinematic theta_vs_p cuts and scpd cuts as summarized in handwritten notes of [02-24-16]
'''
USAGE="make_plots_thesis.py"

ROOT.gROOT.ProcessLine(".L misc-CINT-funcs.C")

#! *** Prepare DATADIR and OUTDIR***
DATADIR=os.path.join(os.environ['STUDY_DET_INEFF_E16'],'results_top2')
print "DATADIR=",DATADIR

DATE=datetime.datetime.now().strftime('%m%d%y')
OUTDIR=os.path.join(DATADIR,"thesis_plots_%s"%DATE)
#OUTDIR=os.path.join("/tmp/thesis_plots_%s"%DATE)
if not os.path.exists(OUTDIR):
		os.makedirs(OUTDIR)
print "OUTDIR=",OUTDIR
#! ***

#! *** Setup some constants
NCUTS=2
PRE,PST=range(NCUTS)
CUT_NAME=["pre_cut","pst_cut"]

NPRTCLS=3
E,P,PIP=range(NPRTCLS)
PRTCL_NAME=["e","p","pip"]

NDTYP=2
ER,SR=range(NDTYP)
DTYPE_NAME=["ER","SR"]

NSCTR=6

#! *** Define function to return hists from canvas
def get_hists_from_canvas(c,idtyp,iprtcl,cut):
	'''
	Each canvas has 6 pads for each sector, and within each pad are contained two types of hists:
		1. theta_vs_p: 1 only
		2. theta_vs_p_pd: Varies 
		3. cut lines: Varies

	Return ,per sector,4 objects:
	1. and normalized version of 1, 2, and 3
	'''
	#! declare list of objects 
	h  =[0 for isctr in range(NSCTR)]
	hn =[0 for isctr in range(NSCTR)]
	hpd=[[] for isctr in range(NSCTR)]
	ln =[[] for isctr in range(NSCTR)]

	#! Get canvas per sector
	for isctr in range(NSCTR):
			#print "Getting hists for canvas for sector",isctr+1
			print "c_%s_theta_vs_p_%s_%s_2.00-5.00_%d"%(PRTCL_NAME[iprtcl],'' if cut==PRE else 'pdall',DTYPE_NAME[idtyp],isctr+1)
			csctr=c.GetPrimitive("c_%s_theta_vs_p%s%s_2.00-5.00_%d"%(PRTCL_NAME[iprtcl],'_' if cut==PRE else '_pdall_',DTYPE_NAME[idtyp],isctr+1))
			#! Get list of primitives from the canvas
			lop=csctr.GetListOfPrimitives()
			#! Loop over these primitives
			for obj in lop:
					objname=obj.GetName()
					#print "Reading: ",obj.GetName()
					if obj.InheritsFrom("TH1"):
							#print "histo: ",obj.GetName()
							if "pd" in objname:
									#print "pd histogram"
									hpd[isctr].append(csctr.GetPrimitive(objname))
							else:
									h[isctr]=csctr.GetPrimitive(objname)
									hn[isctr]=ROOT.norm2D(h[isctr])
									#! Set  titles
									#! + Don't know where they got messed up
									for htmp in [h[isctr],hn[isctr]]:
										xt="p [GeV]"#htmp.GetXaxis().GetTitle()
										yt="#theta [deg]"#htmp.GetYaxis().GetTitle()
										htmp.SetTitle(";%s;%s"%(xt,yt))
										htmp.GetXaxis().SetTitleOffset(0.7)
										htmp.GetYaxis().SetTitleOffset(1.3)
					elif obj.InheritsFrom("TF1"):
						#! Have to rename each object because it appears that in .root file each line was not given a unique name
						#! which causes problems in storing them here
						objc=csctr.GetPrimitive(objname)
						objc.SetName("objname_%d"%(len(ln[isctr])+1))
						objcname=objc.GetName()
						# print "TF1 orig name:",objname
						# print "TF1 copy name:",objcname
						ln[isctr].append(objc)
	return h,hn,hpd,ln

#! ***

#! *** Set up files to get data 
FIN=[0 for i in range(NCUTS)]
FIN[PRE]=ROOT.TFile("%s/kinematics.root"%DATADIR)
FIN[PST]=ROOT.TFile("%s/sc_pd_theta_vs_p/bad_good_sc_pd_theta_vs_p.root"%DATADIR)
print "FIN[PRE]=",FIN[PRE].GetName()
print "FIN[PST=]",FIN[PST].GetName()
#! ***

#! *** General ROOT plotting aesthetics
ROOT.gStyle.SetOptStat(0)

#! Now get relevant theta vs p TCanvas from the files and save them
for cut in range(NCUTS):
	#! Get TCanvas
	#c=[[0 for j in range(NDTYP)] for i in range(NPRTCLS)]
	if cut==PRE:
		for idtyp in range(NDTYP):
			for iprtcl in range(NPRTCLS):
				c=FIN[cut].Get("2.00-5.00/c_%s_theta_vs_p_%s_2.00-5.00"%(PRTCL_NAME[iprtcl],DTYPE_NAME[idtyp]))
				h,hn,hpd,ln=get_hists_from_canvas(c,idtyp,iprtcl,cut)
				#! Save hists in a canvas each for theta_vs_p and theta_vs_p_norm
				cc=ROOT.TCanvas("cc","cc",2000,1500)
				cc.Divide(3,2)
				ccn=ROOT.TCanvas("ccn","ccn",2000,1500)
				ccn.Divide(3,2)
				for isctr in range(NSCTR):
					#! theta_vs_p
					cc.cd(isctr+1)
					h[isctr].Draw("colz")
					#! theta_vs_p_norm
					ccn.cd(isctr+1)
					hn[isctr].Draw("colz")
				cc.SaveAs("%s/theta_vs_p_%s_%s_%s.pdf"%(OUTDIR,PRTCL_NAME[iprtcl],DTYPE_NAME[idtyp],CUT_NAME[cut]))
				ccn.SaveAs("%s/theta_vs_p_norm_%s_%s_%s.pdf"%(OUTDIR,PRTCL_NAME[iprtcl],DTYPE_NAME[idtyp],CUT_NAME[cut]))
	elif cut==PST:
		for idtyp in range(NDTYP):
			for iprtcl in range(NPRTCLS):
				c=FIN[cut].Get("2.00-5.00/pdall/c_%s_theta_vs_p_pdall_%s_2.00-5.00"%(PRTCL_NAME[iprtcl],DTYPE_NAME[idtyp]))
				h,hn,hpd,ln=get_hists_from_canvas(c,idtyp,iprtcl,cut)
				#! Save hists in a canvas each for theta_vs_p and theta_vs_p_norm
				cc=ROOT.TCanvas("cc","cc",2000,1500)
				cc.Divide(3,2)
				ccn=ROOT.TCanvas("ccn","ccn",2000,1500)
				ccn.Divide(3,2)
				for isctr in range(NSCTR):
					#! theta_vs_p
					cc.cd(isctr+1)
					h[isctr].Draw("colz")
					for hp in hpd[isctr]:
						hp.Draw("same")
					for l in ln[isctr]:
						l.Draw("same")
					#! theta_vs_p_norm
					ccn.cd(isctr+1)
					hn[isctr].Draw("colz")
					for hp in hpd[isctr]:
						hp.Draw("same")
					for l in ln[isctr]:
						l.Draw("same")
				cc.SaveAs("%s/theta_vs_p_%s_%s_%s.pdf"%(OUTDIR,PRTCL_NAME[iprtcl],DTYPE_NAME[idtyp],CUT_NAME[cut]))
				ccn.SaveAs("%s/theta_vs_p_norm_%s_%s_%s.pdf"%(OUTDIR,PRTCL_NAME[iprtcl],DTYPE_NAME[idtyp],CUT_NAME[cut]))

			# c[E][idtyp]=FIN[cut].Get("2.00-5.00/pdall/c_e_theta_vs_p_pdall_%s_2.00-5.00"%DTYPE_NAME[idtyp])
			# mod_canvas(c[E][idtyp])
			# c[P][idtyp]=FIN[cut].Get("2.00-5.00/pdall/c_p_theta_vs_p_pdall_%s_2.00-5.00"%DTYPE_NAME[idtyp])
			# mod_canvas(c[P][idtyp])
			# c[PIP][idtyp]=FIN[cut].Get("2.00-5.00/pdall/c_pip_theta_vs_p_pdall_%s_2.00-5.00"%DTYPE_NAME[idtyp])
			# mod_canvas(c[PIP][idtyp])
	#! Now save TCanvas
	#outdir=("%s/%s"%(OUTDIR,CUT_NAME[cut]))
	#if not os.path.exists(outdir):
		#	os.makedirs(outdir)
	# for prtcl in range(NPRTCLS):
	# 	for idtyp in range(NDTYP):
	# 	#c[prtcl].SetCanvasSize(2000,1500)
	# 		c[prtcl][idtyp].Draw()
	# 		c[prtcl][idtyp].SaveAs("%s/theta_vs_p_%s_%s_%s.pdf"%(OUTDIR,PRTCL_NAME[prtcl],DTYPE_NAME[idtyp],CUT_NAME[cut]))

#! If wanting to keep TCanvas open till program exits                           
#if not ROOT.gROOT.IsBatch():
#        plt.show()
#        # wait for you to close the ROOT canvas before exiting
#        wait(True)
			
	


