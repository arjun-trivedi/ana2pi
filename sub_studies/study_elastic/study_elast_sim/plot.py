#!/usr/bin/python
from __future__ import division

import os,sys
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict

import numpy as np

NDTYP=5
ER,SR,SRE16,ST,STE16=range(NDTYP)
DTYP_NAME=['ER','SR','SRE16','ST','STE16']
clr=['kBlue','kRed','kBlack','kRed','kBlack']

F=[0 for i in range(NDTYP)]
F[ER]=ROOT.TFile("%s/test_elast_sim_073015/delastR.root"%os.environ['DELASTDIR_EXP'])
F[SR]=ROOT.TFile("%s/test_elast_sim_073015/delastR.root"%os.environ['DELASTDIR_SIM'])
F[SRE16]=ROOT.TFile("%s/test_elast_sim_073015/e16/delastR.root"%os.environ['DELASTDIR_SIM'])
F[ST]=ROOT.TFile("%s/test_elast_sim_073015/delastT.root"%os.environ['DELASTDIR_SIM'])
F[STE16]=ROOT.TFile("%s/test_elast_sim_073015/e16/delastT.root"%os.environ['DELASTDIR_SIM'])

helast=[0 for i in range(NDTYP)]
hW=[0 for i in range(NDTYP)]
for i in range(5):
	helast[i]=F[i].Get("delast/monitor/helastic")
	hW[i]=F[i].Get("delast/monitor/hW")
	helast[i].SetLineColor(ROOT.gROOT.ProcessLine("%s"%clr[i]))
	hW[i].SetLineColor(ROOT.gROOT.ProcessLine("%s"%clr[i]))
	helast[i].SetName("%s_%s"%(helast[i].GetName,DTYP_NAME[i]))

c=ROOT.TCanvas()
c.Divide(2,2)

l=ROOT.TLegend(0.6,0.7,0.8,0.8)
c.cd(1)
#helast[ER].Draw()
helast[SR].Draw("")
helast[SRE16].Scale(helast[SR].GetMaximum()/helast[SRE16].GetMaximum())
helast[SRE16].Draw("same")
if len(sys.argv)==2:
	helast[ER].Scale(helast[SR].GetMaximum()/helast[ER].GetMaximum())
	helast[ER].Draw("same")
l.AddEntry(helast[SR],"AT-SR-E1F")#DTYP_NAME[SR])
l.AddEntry(helast[SRE16],"EI-SR-E16")#DTYP_NAME[SRE16])
if len(sys.argv)==2:
	l.AddEntry(helast[ER],"ER")#DTYP_NAME[ER])
l.Draw()

#helast[SR].DrawNormalized("",1)
#helast[ER].DrawNormalized("sames",1)
c.cd(2)
hW[SR].Draw("")
hW[SRE16].Scale(hW[SR].GetMaximum()/hW[SRE16].GetMaximum())
hW[SRE16].Draw("same")

c.cd(3)
helast[ST].Draw("")
s=helast[ST].GetMaximum()/helast[STE16].GetMaximum()
helast[STE16].Scale(s)
helast[STE16].Draw("same")

c.cd(4)
hW[ST].Draw("")
#hW[STE16].Scale(hW[ST].GetMaximum()/hW[STE16].GetMaximum())
hW[STE16].Scale(s)
hW[STE16].Draw("same")		

#! If wanting to keep TCanvas open till program exits				
if not ROOT.gROOT.IsBatch():
	plt.show()
	# wait for you to close the ROOT canvas before exiting
	wait(True)

#if __name__ == "__main__":	
#	if len(sys.argv)==2:
#		plot_fid(sys.argv[1])
#	elif len(sys.argv)==3:
#		plot_fid(sys.argv[1],int(sys.argv[2])) 
	
