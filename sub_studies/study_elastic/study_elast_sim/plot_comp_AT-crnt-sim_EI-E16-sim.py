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

c.cd(1) #! W:ST
hW[ST].Draw("")
sST=helast[ST].GetMaximum()/helast[STE16].GetMaximum()
#print "sT=",sST
hW_STE16=hW[STE16].Clone()
if sys.argv[1]=="s":hW_STE16.Scale(sST)
hW_STE16.Draw("same")

c.cd(2) #! elast:ST
helast[ST].Draw("")
helast_STE16=helast[STE16].Clone()
if sys.argv[1]=="s":helast_STE16.Scale(sST)
helast_STE16.Draw("same")

c.cd(3) #! W:SR
hW[SR].Draw("")
sSR=helast[SR].GetMaximum()/helast[SRE16].GetMaximum()
#print "sR=",sSR
hW_SRE16=hW[SRE16].Clone()
if sys.argv[1]=="s":hW_SRE16.Scale(sSR)
hW_SRE16.Draw("same")

l=ROOT.TLegend(0.6,0.7,0.8,0.8)
c.cd(4) #!elast:SR and ER
helast[SR].Draw("")
helast_SRE16=helast[SRE16].Clone()
if sys.argv[1]=="s":helast_SRE16.Scale(sSR)
helast_SRE16.Draw("same")
l.AddEntry(helast[SR],"AT-SR-E1F")
l.AddEntry(helast_SRE16,"EI-SR-E16")
if len(sys.argv)==3:
	#! Scale ER
	sER=helast[SR].GetMaximum()/helast[ER].GetMaximum()
	helast_ER=helast[ER].Clone()
	if sys.argv[1]=="s":helast_ER.Scale(sER)
	helast_ER.Draw("same")
	l.AddEntry(helast_ER,"ER")
l.Draw()

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
	
