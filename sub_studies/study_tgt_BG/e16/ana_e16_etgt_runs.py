#!/usr/bin/python
from __future__ import division

import os,sys
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict

import numpy as np
import itertools

from rootpy.plotting import Hist, HistStack, Legend, Canvas
import rootpy.plotting.root2matplotlib as rplt
from rootpy import asrootpy

import datetime

'''
+ This scripts analyzes plots the vz distribution from etgt runs given to me by EI (See email of [05-16-17])
+ Given these vz distributions, it can be verified if these indeed are etgt runs
+ Runs given to me by EI:
  30675,30676,30746,30825,30867,30874,30875,30886,30910,30962,31029,31104,31128,31252,31254,31299,31300,31344,31382,31385,31391
'''

USAGE='$SUBSTUDIES/study_tgt_BG/e16/ana_e16_etgt_runs.py >& /tmp/log&' #log messages show list of files in chain for each run
def get_vz_hist(lcn,runnum):
	chain=ROOT.TChain("h10")
	chain.Add("%s/*%d*"%(lcn,runnum))
	arr=chain.GetListOfFiles()
	narr=arr.GetEntries()
	if narr==0: #! empty chain. create empty hist
		htmp=Hist(100,-12,6)
		return htmp
	#! Print out list of files in chain
	for i in range(narr):
		print arr.At(i).GetTitle()
	#! Draw
	chain.Draw("vz[0]>>hvz(100,-12,6)","id[0]==11 && dc[0]>0 && cc[0]>0 && ec[0]>0 && sc[0]>0")
	htmp=ROOT.gDirectory.Get("hvz")
	return htmp

DATE=datetime.datetime.now().strftime('%m%d%y')

RUNL=[30675,30676,30746,30825,30867,30874,30875,30886,30910,30962,31029,
      31104,31128,31252,31254,31299,31300,31344,31382,31385,31391]

OUTDIR=os.environ['STUDY_TGT_BG_E16_DATADIR']

hvz=[0 for i in range(len(RUNL))]
print hvz
for irun,run in enumerate(RUNL):
	print "*** Getting hvz for run %d ***"%run
	if run==30675 or run==30676: lcn="%s/h10-sample"%os.environ['H10DIR_EXP_E16']
	else:                        lcn="%s/e16_etgt_run_files"%os.environ['H10DIR_EXP_E16_2']
	hvz[irun]=get_vz_hist(lcn,run)
	hvz[irun]=asrootpy(hvz[irun])
	hvz[irun].SetName("vz_r%s"%run)
	print "hvz name =",hvz[irun].GetName()
	print "******"
print hvz

NCOLS=6
NROWS=int(len(RUNL)/NCOLS)
if (len(RUNL)%NCOLS)>0: NROWS+=1
RIDX=range(NROWS)
CIDX=range(NCOLS)
GRID=list(itertools.product(RIDX,CIDX))
#print NCOLS,NROWS
#print GRID
fig,axs=plt.subplots(figsize=(12,8),nrows=NROWS,ncols=NCOLS)
fig.subplots_adjust(hspace=.5)
#for irun in range(len(RUNL)):
for irun,grid in enumerate(GRID):
	if irun+1>len(RUNL): break
	l,m=grid[0],grid[1]
	#print irun,l,m
	axs[l][m].set_title("%s"%hvz[irun].GetName())
	rplt.errorbar(hvz[irun], xerr=False, emptybins=False, axes=axs[l][m])
fig.savefig("%s/ana_EI_etgt_runs_%s.png"%(OUTDIR,DATE))

#! If wanting to keep TCanvas open till program exits				
# if not ROOT.gROOT.IsBatch():
# 	plt.show()
# 	# wait for you to close the ROOT canvas before exiting
# 	wait(True)

#if __name__ == "__main__":	
#	if len(sys.argv)==2:
#		plot_fid(sys.argv[1])
#	elif len(sys.argv)==3:
#		plot_fid(sys.argv[1],int(sys.argv[2])) 
	
