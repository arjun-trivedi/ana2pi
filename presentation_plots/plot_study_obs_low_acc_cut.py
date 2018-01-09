#!/usr/bin/python
from __future__ import division

import os,sys,subprocess
import ROOT
from collections import OrderedDict
import math

'''
+ In this script the results for making cuts on low acceptance , <0.5% and < 1%,
on Obs_1D are shown for q2wbin='2.40-3.00_1.725-1.750'

+ NOTE that results from seq=='EC' are used because with 'EF' the sensitivity of this
study is obfuscated.
'''

#! OUTDIR
OUTDIR="%s/figures/Results/"%os.environ['ANANOTE']
#OUTDIR="/tmp"

#! Define some objects that will be needed in steps 1 and step 2
#! PAD_MAP[obs]=(pad,vst,var)
PAD_MAP=OrderedDict()
PAD_MAP['1D']=[(1,1,"M1"),   (2,3,'M2'),   (3,2,'M2'),
			   (4,1,"THETA"),(5,3,'THETA'),(6,2,'THETA'),
			   (7,1,"ALPHA"),(8,3,'ALPHA'),(9,2,'ALPHA')]

fobs=OrderedDict()
obsdir=OrderedDict()
obsdir['nocut']="%s/lowQ2_SSBands_off_off_sim4_sim5_sim6_sim7_sim8_sim13_102617/cutsncors1/sim4_sim5_sim6_sim7_sim8_sim13/Obs_1D_norm"%os.environ['OBSDIR_E16']
obsdir['005']  ="%s/lowQ2_SSBands_off_off_sim4_sim5_sim6_sim7_sim8_sim13_102617_acc_cut_005/cutsncors1/sim4_sim5_sim6_sim7_sim8_sim13/Obs_1D_norm"%os.environ['OBSDIR_E16']
obsdir['01']   ="%s/lowQ2_SSBands_off_off_sim4_sim5_sim6_sim7_sim8_sim13_102617_acc_cut_01/cutsncors1/sim4_sim5_sim6_sim7_sim8_sim13/Obs_1D_norm"%os.environ['OBSDIR_E16']
fobs['nocut']=ROOT.TFile("%s/obs_1D.root"%obsdir['nocut'],"READ")
fobs['005']=ROOT.TFile("%s/obs_1D.root"%obsdir['005'],"READ")
fobs['01']=ROOT.TFile("%s/obs_1D.root"%obsdir['01'],"READ")

CUTNAME={'nocut':"no cut", '005':"<0.5%", '01':"<1.0%"}

hobs=OrderedDict()
seq='EC'
q2wbin='2.40-3.00_1.725-1.750'
for pm in PAD_MAP['1D']:
	vst,var=pm[1],pm[2]
	#! 1D
	for cut in ['nocut','005','01']:
		hobs[cut,seq,vst,var]=fobs[cut].Get("%s/h1_%s_%d_%s"%(q2wbin,seq,vst,var))
	

#! plot
ROOT.gStyle.SetOptStat(0)
#! legend
l=ROOT.TLegend(0.60,0.70,0.90,0.90)
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
pad_p.Divide(3,3)
for item in PAD_MAP['1D']:
	pad,vst,var=item[0],item[1],item[2]
	#print "pad,vst,var=",pad,vst,var
	gpad=pad_p.cd(pad)
	#! Draw
	draw_opt=""
	CLRS=[ROOT.gROOT.ProcessLine("kBlue"),ROOT.gROOT.ProcessLine("kGreen"),ROOT.gROOT.ProcessLine("kRed")]
	for icut,cut in enumerate(['nocut','005','01']):
		hobs[cut,seq,vst,var].SetLineColor(CLRS[icut])
		hobs[cut,seq,vst,var].SetMarkerColor(CLRS[icut])
		if icut>0: draw_opt="same"
		hobs[cut,seq,vst,var].Draw(draw_opt)
	#! legend
	if pad==1:
		for icut,cut in enumerate(['nocut','005','01']):
			l.AddEntry(hobs[cut,seq,vst,var],"%s"%CUTNAME[cut],"p")
		l.Draw()
c.SaveAs("%s/comp_1D_low_acc_cut_%s.png"%(OUTDIR,q2wbin))
#! Convert .png->.eps->.pdf since in .eps and .pdf produced directly by ROOT
#! the y-axis titles are cutoff
# #! 1. .png->.eps
# cmd=['convert',"%s/c_q%s.png"%(outdir,q2wbin),"%s/c_q%s.eps"%(outdir,q2wbin)]
# tool=subprocess.Popen(cmd,stderr=subprocess.STDOUT)
# tool.wait()
#! 2. .png->.pdf
cmd=['convert',"%s/comp_1D_low_acc_cut_%s.png"%(OUTDIR,q2wbin),"%s/comp_1D_low_acc_cut_%s.pdf"%(OUTDIR,q2wbin)]
tool=subprocess.Popen(cmd,stderr=subprocess.STDOUT)
tool.wait()
# <codecell>


