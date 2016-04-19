#!/usr/bin/python
from disp_obs import DispObs
import sys,os
import glob

def plot_SS(q2min,q2max,wmin,wmax,expt,dbg,obsd,SS_tag):
	print "disp_SS::plot_SS(): Going to display SS for q2min q2max wmin wmax expt dbg=",q2min,q2max,wmin,wmax,expt,dbg
	print "obsd (pretty print):"
	for k in obsd:
		obsnum,obsdir,obstag=k,obsd[k][0],obsd[k][1]
		print obsnum,obsdir,obstag
	print "SS_tag=",SS_tag

	#! Just to be safe, in the sense of resources (mainly memory), use a new instance of
	#! of DispObs for disp_1D_SS() and disp_R2_SS()
	#! 1. First do disp_1D_SS()
	do1=DispObs('NA','NA','NA',q2min,q2max,wmin,wmax,expt,dbg,obsd=obsd,SS_tag=SS_tag)
	do1.disp_1D_SS()
	#! 2. Now do disp_R2_SS()
	do2=DispObs('NA','NA','NA',q2min,q2max,wmin,wmax,expt,dbg,obsd=obsd,SS_tag=SS_tag)
        do2.disp_R2_SS()
