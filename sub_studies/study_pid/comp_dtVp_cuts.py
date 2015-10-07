#!/usr/bin/python
from __future__ import division

import os,sys
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict

import numpy as np

#! get .root files
fin=[]
for i in range(len(sys.argv)):
        if i==0: continue
	fin.append(ROOT.TFile(sys.argv[i]))

NDTYP=2
EXP,SIM=range(NDTYP)
DTYP_NAME=["exp","sim"]

NPRTCL=3
P,PIP,PIM=range(NPRTCL)
PRTCL_NAME=["p","pip","pim"]

#! get cut functions
cutl=[[[0 for iprtcl in range(NPRTCL)]for idtyp in range(NDTYP)] for ifin in range(len(fin))]
cuth=[[[0 for iprtcl in range(NPRTCL)]for idtyp in range(NDTYP)] for ifin in range(len(fin))]

for i in range(len(fin)):
	for j in range(NDTYP):
		for k in range(NPRTCL):
			canvas=fin[i].Get("%s/%s/c_cut_%s_%s"%(DTYP_NAME[j],PRTCL_NAME[k],DTYP_NAME[j],PRTCL_NAME[k]))
			cutl[i][j][k]=canvas.GetPrimitive("fl")		
			cuth[i][j][k]=canvas.GetPrimitive("fh")

#! Draw cut functions
c=[0 for i in range(NDTYP)]
l=[[0 for j in range(NPRTCL)] for i in range(NDTYP)]
for i in range(NDTYP):
	c[i]=ROOT.TCanvas(DTYP_NAME[i],DTYP_NAME[i])
	c[i].Divide(2,2)
	for j in range(NPRTCL):
		c[i].cd(j+1)
		l[i][j]=ROOT.TLegend(0.5,0.7,0.8,0.9)
		ctr=0
		for k in range(len(fin)):
			draw_opt=""
			if ctr>0: draw_opt="same"
			#! Draw cutl
			#! aesthetics
			cutl[k][i][j].SetMinimum(-5)
			cutl[k][i][j].SetMaximum(5)
			cutl[k][i][j].SetLineColor(k+1)
  			#! draw
			cutl[k][i][j].Draw(draw_opt)
			#! Draw cuth
			cuth[k][i][j].SetMinimum(-5)
                        cuth[k][i][j].SetMaximum(5)
			cuth[k][i][j].SetLineColor(k+1)
			#! draw
                        cuth[k][i][j].Draw("same")	
			
			#! legend
			#if (j==0):
			fname=fin[k].GetName()
			l[i][j].SetHeader("%s_%s"%(DTYP_NAME[i],PRTCL_NAME[j]))
			tag=fname[fname.find("dtVp_cuts"):fname.find('fout')]
			l[i][j].AddEntry(cutl[k][i][j],tag,"l")

			ctr+=1
		#if j==0:
		l[i][j].Draw("same")

#t[t.find("dtVp_cuts"):t.find('fout')]	
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
	
