#!/usr/bin/python
from __future__ import division

import os,sys
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict

import numpy as np
NENTRIES=1000000000
if len(sys.argv)>1:# i.e. nentries entered by user
	NENTRIES=int(sys.argv[1])
print "NENTRIES=",NENTRIES

FEXP=ROOT.TFile(os.path.join(os.environ['DELASTDIR_EXP'],'study_elast_xsec_081715/delastR_tree.root'))
T=FEXP.Get("delast/monitor/t")

NSCTR=6

FOUT=ROOT.TFile("out.root","RECREATE")
#! Study nphe for W[0,0.848] & W=[0.848,1.028]
NWCUT=2
w_cut=["W<0.848","W>0.848&&W<1.028"]
w_cut_name=['W_lt_0848','W_gt_0848_lt_1028']

ROOT.gStyle.SetOptStat("nemMrRuo")
#! Study nphe for W[0,0.848] & W=[0.848,1.028]
h_nphe=[[[]for j in range(NSCTR)]for i in range(NWCUT)]
for iw in range(NWCUT):
	for isctr in range(NSCTR):
		print w_cut[iw]
		cut=ROOT.TCut(w_cut[iw])
		cut+=ROOT.TCut("sector==%d"%(isctr+1))
		T.Draw("nphe>>h(100,0,300)",cut,"",NENTRIES)
		htmp=ROOT.gDirectory.Get("h")
		h_nphe[iw][isctr]=htmp.Clone()
c_nphe=[[]for i in range(NWCUT)]
for iw in range(NWCUT):
	cname=("nphe_%s"%w_cut_name[iw])
	c_nphe[iw]=ROOT.TCanvas(cname,cname)
	c_nphe[iw].Divide(3,2)
	for isctr in range(NSCTR):
		c_nphe[iw].cd(isctr+1)
		h_nphe[iw][isctr].Draw("HIST")
	FOUT.WriteTObject(c_nphe[iw])

#ROOT.gStyle.SetOptStat(0)
#! Study theta for W[0,0.848] & W=[0.848,1.028] X no cc-cut & cc cut
NCCCUT=3
cc_cut=["","cc>0","cc>0&nphe>30"]
cc_cut_name=["cc_none","cc_gt_0","cc_gt_0_nphe_gt_30"]
cc_cut_clr=["kBlue","kGreen","kRed"]
h_theta=[[[[]for k in range(NSCTR)]for j in range(NCCCUT)]for i in range(NWCUT)]
for iw in range(NWCUT):
        for isctr in range(NSCTR):
                print w_cut[iw]
		for icc in range(NCCCUT):
                	cut=ROOT.TCut(w_cut[iw])
			cut+=cc_cut[icc]
                	cut+=ROOT.TCut("sector==%d"%(isctr+1))
                	T.Draw("theta_e>>h(32,14,46)",cut,"",NENTRIES)
                	htmp=ROOT.gDirectory.Get("h")
                	h_theta[iw][icc][isctr]=htmp.Clone()
c_theta=[[]for i in range(NWCUT)]
c_theta_rto=[[]for i in range(NWCUT)]
h_rto_0_0=[[[]for j in range(NSCTR)]for i in range(NWCUT)]
h_rto_1_0=[[[]for j in range(NSCTR)]for i in range(NWCUT)]
h_rto_2_0=[[[]for j in range(NSCTR)]for i in range(NWCUT)]
l=[[]for i in range(NWCUT)]
l_rto=[[]for i in range(NWCUT)]
for iw in range(NWCUT):
        cname=("theta_%s"%w_cut_name[iw])
        c_theta[iw]=ROOT.TCanvas(cname,cname)
        c_theta[iw].Divide(3,2)
	cname_rto=("theta_rto_%s"%w_cut_name[iw])
        c_theta_rto[iw]=ROOT.TCanvas(cname_rto,cname_rto)
        c_theta_rto[iw].Divide(3,2)
	l[iw]=ROOT.TLegend(0.6,0.7,1.0,0.9)
	l_rto[iw]=ROOT.TLegend(0.6,0.7,1.0,0.9)
        for isctr in range(NSCTR):
                c_theta[iw].cd(isctr+1)
		for icc in range(NCCCUT):
			h_theta[iw][icc][isctr].SetLineColor(ROOT.gROOT.ProcessLine("%s"%cc_cut_clr[icc]))
	                h_theta[iw][icc][isctr].SetName("%s"%cc_cut_name[icc])
			if icc==0:
				 h_theta[iw][icc][isctr].Draw("HIST")
			else:
				 h_theta[iw][icc][isctr].Draw("HIST same")
		if (isctr==0):
                        for icc in range(NCCCUT):
                                l[iw].AddEntry(h_theta[iw][icc][isctr],h_theta[iw][icc][isctr].GetName(),"l")
                        l[iw].Draw()
		#! rto
		#! First set Sumw2() for correct error prop
		h_theta[iw][0][isctr].Sumw2()
		h_theta[iw][1][isctr].Sumw2()
		h_theta[iw][2][isctr].Sumw2()
		#! Now make rto
		c_theta_rto[iw].cd(isctr+1)
		h_rto_0_0[iw][isctr]=h_theta[iw][0][isctr].Clone()
                h_rto_0_0[iw][isctr].Divide(h_theta[iw][0][isctr])
		h_rto_1_0[iw][isctr]=h_theta[iw][1][isctr].Clone()
		h_rto_1_0[iw][isctr].Divide(h_theta[iw][0][isctr])
		h_rto_2_0[iw][isctr]=h_theta[iw][2][isctr].Clone()
                h_rto_2_0[iw][isctr].Divide(h_theta[iw][0][isctr])
		h_rto_0_0[iw][isctr].SetMinimum(0)
		h_rto_0_0[iw][isctr].SetMaximum(2)
		h_rto_1_0[iw][isctr].SetMinimum(0)
                h_rto_1_0[iw][isctr].SetMaximum(2)
		h_rto_2_0[iw][isctr].SetMinimum(0)
                h_rto_2_0[iw][isctr].SetMaximum(2)
		h_rto_0_0[iw][isctr].Draw("")
		h_rto_1_0[iw][isctr].Draw("same")
		h_rto_2_0[iw][isctr].Draw("same")
		if (isctr==0):
                        l_rto[iw].AddEntry(h_rto_0_0[iw][isctr],h_rto_0_0[iw][isctr].GetName(),"l")
			l_rto[iw].AddEntry(h_rto_1_0[iw][isctr],h_rto_1_0[iw][isctr].GetName(),"l")
			l_rto[iw].AddEntry(h_rto_2_0[iw][isctr],h_rto_2_0[iw][isctr].GetName(),"l")
                        l_rto[iw].Draw()
        FOUT.WriteTObject(c_theta[iw])
	FOUT.WriteTObject(c_theta_rto[iw])

#! Study SF for W[0,0.848] & W=[0.848,1.028]
h_sf=[[[]for j in range(NSCTR)]for i in range(NWCUT)]
for iw in range(NWCUT):
        for isctr in range(NSCTR):
                print w_cut[iw]
                cut=ROOT.TCut(w_cut[iw])
                cut+=ROOT.TCut("sector==%d"%(isctr+1))
                T.Draw("etot/p:p>>h(100,0,10,100,0,0.5)",cut,"colz",NENTRIES)
                htmp=ROOT.gDirectory.Get("h")
                h_sf[iw][isctr]=htmp.Clone()
c_sf=[[]for i in range(NWCUT)]
for iw in range(NWCUT):
        cname=("sf_%s"%w_cut_name[iw])
        c_sf[iw]=ROOT.TCanvas(cname,cname)
        c_sf[iw].Divide(3,2)
        for isctr in range(NSCTR):
                c_sf[iw].cd(isctr+1)
                h_sf[iw][isctr].Draw("colz")
        FOUT.WriteTObject(c_sf[iw])

#! If wanting to keep TCanvas open till program exits				
if not ROOT.gROOT.IsBatch():
	plt.show()
 	# wait for you to close the ROOT canvas before exiting
 	wait(True)

	
