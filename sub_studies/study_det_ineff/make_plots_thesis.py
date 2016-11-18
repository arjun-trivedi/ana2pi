#!/usr/bin/python
from __future__ import division

import os,sys
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict
import itertools

import numpy as np

import math

import atlib

'''
+ This script simply consolidates plots, listed below, already produced by 'study_det_ineff.py' in .root files in '$STUDY_DET_INEFF_E16/results_top2' and saves them directly as .pdf in '$STUDY_DET_INEFF_E16/results_top2/thesis_plots' that can be copied over to the thesis folder. Plots consolidated:
  + "pre-cut" plots in '$STUDY_DET_INEFF_E16/results_top2/kinematics.root'
  + "pst-cut" plots in '$STUDY_DET_INEFF_E16/results_top2/sc_pd_theta_vs_p/bad_good_sc_pd_theta_vs_p.root"
NOTE that "cut" here refers to kinematic theta_vs_p cuts and scpd cuts as summarized in handwritten notes of [02-24-16]
'''
USAGE="make_plots_thesis.py"

#! *** Prepare DATADIR and OUTDIR***
DATADIR=os.path.join(os.environ['STUDY_DET_INEFF_E16'],'results_top2')
print "DATADIR=",DATADIR

OUTDIR=os.path.join(DATADIR,"thesis_plots")
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
#! ***

#! *** Set up files to get data 
FIN=[0 for i in range(NCUTS)]
FIN[PRE]=ROOT.TFile("%s/kinematics.root"%DATADIR)
FIN[PST]=ROOT.TFile("%s/sc_pd_theta_vs_p/bad_good_sc_pd_theta_vs_p.root"%DATADIR)
print "FIN[PRE]=",FIN[PRE].GetName()
print "FIN[PST=]",FIN[PST].GetName()
#! ***

#! Define function to modify histograms in cavnas by adding axes titles
def mod_canvas(c):
	for pnum in [1,2,3,4,5,6]:
                p=c.GetPad(pnum)
                l=p.GetListOfPrimitives()
                for obj in l:
                        if (obj.InheritsFrom("TH1")):
                                print obj
                                #! Axes title
                                obj.SetXTitle("p [GeV]")
				#! The following line to SetLeftMargin is purposefully commented out
				#! because of a techincal issue
                                #p.SetLeftMargin(0.20)
                                obj.GetYaxis().SetTitleOffset(1.5)
                                obj.SetYTitle("#theta [deg]")
#! Now get relevant theta vs p TCanvas from the files and save them
for cut in range(NCUTS):
	#! Get TCanvas
	c=[0 for i in range(NPRTCLS)]
	if cut==PRE:
		c[E]=FIN[cut].Get("2.00-5.00/c_e_theta_vs_p_ER_2.00-5.00")
		mod_canvas(c[E])
		c[P]=FIN[cut].Get("2.00-5.00/c_p_theta_vs_p_ER_2.00-5.00")
		mod_canvas(c[P])
		c[PIP]=FIN[cut].Get("2.00-5.00/c_pip_theta_vs_p_ER_2.00-5.00")
		mod_canvas(c[PIP])
	elif cut==PST:
		c[E]=FIN[cut].Get("2.00-5.00/pdall/c_e_theta_vs_p_pdall_ER_2.00-5.00")
		mod_canvas(c[E])
                c[P]=FIN[cut].Get("2.00-5.00/pdall/c_p_theta_vs_p_pdall_ER_2.00-5.00")
		mod_canvas(c[P])
		c[PIP]=FIN[cut].Get("2.00-5.00/pdall/c_pip_theta_vs_p_pdall_ER_2.00-5.00")
		mod_canvas(c[PIP])
	#! Now save TCanvas
	#outdir=("%s/%s"%(OUTDIR,CUT_NAME[cut]))
	#if not os.path.exists(outdir):
        #	os.makedirs(outdir)
	for prtcl in range(NPRTCLS):
		#c[prtcl].SetCanvasSize(2000,1500)
		c[prtcl].Draw()
		c[prtcl].SaveAs("%s/theta_vs_p_%s_%s.pdf"%(OUTDIR,PRTCL_NAME[prtcl],CUT_NAME[cut]))

#! If wanting to keep TCanvas open till program exits                           
#if not ROOT.gROOT.IsBatch():
#        plt.show()
#        # wait for you to close the ROOT canvas before exiting
#        wait(True)
			
	


