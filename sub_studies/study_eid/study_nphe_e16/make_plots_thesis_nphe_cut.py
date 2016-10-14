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
+ This script simply consolidates plots, listed below, already produced in .root files in '$STUDY_EID_NPHE_E16_DATADIR/study_nphe_cut.root' and saves them directly as .pdf in '$STUDY_EID_NPHE_E16_DATADIR/thesis_plots' that can be copied over to the thesis folder.
'''
USAGE="make_plots_thesis_nphe_cut.py"

#! *** Prepare DATADIR and OUTDIR***
DATADIR=os.path.join(os.environ['STUDY_EID_NPHE_E16_DATADIR'])
print "DATADIR=",DATADIR

OUTDIR=os.path.join(DATADIR,"thesis_plots")
if not os.path.exists(OUTDIR):
        os.makedirs(OUTDIR)
print "OUTDIR=",OUTDIR
#! ***

#! *** Setup some constants
NSCTRS=6
NPMTS=3
IPMT_L,IPMT_C,IPMT_R=range(NPMTS)
#! ***
#! *** Set up files to get data 
F=ROOT.TFile("%s/study_nphe_cut.root"%DATADIR)
print "F=",F.GetName()
#! ***

#! Now get relevant theta vs p TCanvas from the files and save them
for isctr in range(NSCTRS):
	sctr=isctr+1
	#! Get TCanvas for each PMT
	c=[0 for i in range(NPMTS)]
	for ipmt in range(NPMTS):
		pmt=ipmt+1
		c[ipmt]=F.Get("sector%d/sct%d_pmt%d"%(sctr,sctr,pmt))
		print c[ipmt].GetName()
	#! Now save TCanvas
	for ipmt in range(NPMTS):
		pmt=ipmt+1
		c[ipmt].Draw()
		c[ipmt].SaveAs("%s/nphe_cut_sctr%d_pmt%d.pdf"%(OUTDIR,sctr,ipmt+1))

#! If wanting to keep TCanvas open till program exits                           
#if not ROOT.gROOT.IsBatch():
#        plt.show()
#        # wait for you to close the ROOT canvas before exiting
#        wait(True)
			
	


