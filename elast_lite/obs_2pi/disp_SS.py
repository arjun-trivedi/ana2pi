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

	do=DispObs('NA','NA','NA',q2min,q2max,wmin,wmax,expt,dbg,obsd=obsd,SS_tag=SS_tag)
	do.disp_1D_SS()
