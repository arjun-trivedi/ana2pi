#!/usr/bin/python
from __future__ import division

import os,sys,subprocess
import ROOT
from collections import OrderedDict
import math

'''
In this script the equivalence for Obs_1D and Obs_R2:A is shown 
for q2wbin='2.40-3.00_1.725-1.750'
'''

#! OUTDIR
OUTDIR="%s/figures/Results/"%os.environ['ANANOTE']

#! Define some objects that will be needed in steps 1 and step 2
#! PAD_MAP[obs]=(pad,vst,var)
PAD_MAP=OrderedDict()
PAD_MAP['1D']=[(1,1,"M1"),   (2,3,'M2'),   (3,2,'M2'),
			   (4,1,"THETA"),(5,3,'THETA'),(6,2,'THETA'),
			   (7,1,"ALPHA"),(8,3,'ALPHA'),(9,2,'ALPHA')]

PAD_MAP['R2']=[(1,1,"M1"),    (2,3,'M1'),    (3,2,'M1'), 
			   (4,1,"M2"),    (5,3,'M2'),    (6,2,'M2'),
			   (7,1,"THETA"), (8,3,'THETA'), (9,2,'THETA'),
			   (10,1,"ALPHA"),(11,3,'ALPHA'),(12,2,'ALPHA')]

R2_NAMED={'A':'R2_{T}+R2_{L}','B':'R2^{c}_{LT}','C':'R2^{c}_{TT}','D':'R2^{s}_{LT}','E':'R2^{s}_{TT}'}

fobs=OrderedDict()
#OBSDIR="%s/thesis_obs_norm_ST_shape_122817"%os.environ['OBSDIR_E162']
OBSDIR="%s/thesis_obs_norm_ST_shape_041818"%os.environ['OBSDIR_E162']
fobs['1D']=ROOT.TFile("%s/lowQ2/Obs_1D/Obs_1D.root"%OBSDIR,"READ")
fobs['R2']=ROOT.TFile("%s/lowQ2/Obs_R2/Obs_R2.root"%OBSDIR,"READ")

# <codecell>

hobs=OrderedDict()
seq='EF'
R2L=['A','B','C','D','E']
q2wbin='2.40-3.00_1.725-1.750'
for pm in PAD_MAP['1D']:
	vst,var=pm[1],pm[2]
	#! 1D
	hobs['1D',seq,vst,var]=fobs['1D'].Get("%s/h1_%s_%d_%s"%(q2wbin,seq,vst,var))
	print hobs['1D',seq,vst,var].GetName()
	#! R2:A only
	hobs['R2A',seq,vst,var]=fobs['R2'].Get("%s/hR2_A_%s_%d_%s"%(q2wbin,seq,vst,var))
	print hobs['R2A',seq,vst,var].GetName()



# <codecell>

#! plot
ROOT.gStyle.SetOptStat(0)
#! legend
l=ROOT.TLegend(0.30,0.75,0.90,0.90)
l.SetFillStyle(0)
#l.SetBorderSize(0)
#l.SetTextSize(0.06)
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
	#! Draw 1D
	print hobs['1D',seq,vst,var].GetName()
	hobs['1D',seq,vst,var].Draw()
	#! Draw R2A
	#! First scale by 2pi
	hobs['R2A',seq,vst,var].Scale(2*math.pi)
	hobs['R2A',seq,vst,var].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
	hobs['R2A',seq,vst,var].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
	hobs['R2A',seq,vst,var].Draw("same")
	#! legend
	if pad==1:
		l.AddEntry(hobs['1D',seq,vst,var],"Single-Differential","p")
		l.AddEntry(hobs['R2A',seq,vst,var],"(%s)#upoint2#pi"%R2_NAMED['A'],"p")
		l.Draw()
c.SaveAs("%s/comp_1D_R2A_%s.png"%(OUTDIR,q2wbin))
#! Convert .png->.eps->.pdf since in .eps and .pdf produced directly by ROOT
#! the y-axis titles are cutoff
# #! 1. .png->.eps
# cmd=['convert',"%s/c_q%s.png"%(outdir,q2wbin),"%s/c_q%s.eps"%(outdir,q2wbin)]
# tool=subprocess.Popen(cmd,stderr=subprocess.STDOUT)
# tool.wait()
#! 2. .png->.pdf
cmd=['convert',"%s/comp_1D_R2A_%s.png"%(OUTDIR,q2wbin),"%s/comp_1D_R2A_%s.pdf"%(OUTDIR,q2wbin)]
tool=subprocess.Popen(cmd,stderr=subprocess.STDOUT)
tool.wait()
# <codecell>


