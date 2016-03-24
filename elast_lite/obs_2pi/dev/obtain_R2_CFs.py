#!/usr/bin/python
from __future__ import division

import os,sys
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict
import math
import numpy as np

import scipy.integrate as integrate

def print_integrand(R2):
	if R2!='A' and R2!='B' and R2!='C' and R2!='D':
		print "get_R2() not implemented for R2=",R2

	phi_bin_lel=range(0,360,36)
	print phi_bin_lel

	integrandl=[]
	for phi in phi_bin_lel:
		if R2=='A':
			val=integrate.quad(lambda x: 1,math.radians(phi),math.radians(phi+36))
		elif R2=='B':
			val=integrate.quad(lambda x: math.cos(x),math.radians(phi),math.radians(phi+36))
		elif R2=='C':
			val=integrate.quad(lambda x: math.cos(2*x),math.radians(phi),math.radians(phi+36))
		elif R2=='D':
			val=integrate.quad(lambda x: math.sin(x),math.radians(phi),math.radians(phi+36))
		integrandl.append(val[0])
	print integrandl

	fctrl=[]
	for ibin,phi in enumerate(phi_bin_lel):
		if R2=='A':
			fctrl.append("%.2f"%(1/integrandl[ibin]))
		elif R2=='B':
			fctrl.append("%.2f"%(math.cos(math.radians(phi+(36/2)))/integrandl[ibin]))
		elif R2=='C':
			fctrl.append("%.2f"%(math.cos(2*math.radians(phi+(36/2)))/integrandl[ibin]))
		elif R2=='D':
			fctrl.append("%.2f"%(math.sin(math.radians(phi+(36/2)))/integrandl[ibin]))
	print fctrl

	print "for mthd2"
	integrand_mbyl=[]
	for ibin,phi in enumerate(phi_bin_lel):
                if R2=='A':
                        integrand_mbyl.append(1*integrandl[ibin])
                elif R2=='B':
                        integrand_mbyl.append(math.cos(math.radians(phi))*integrandl[ibin])
                elif R2=='C':
                        integrand_mbyl.append(math.cos(2*math.radians(phi))*integrandl[ibin])
                elif R2=='D':
                        integrand_mbyl.append(math.sin(math.radians(phi))*integrandl[ibin])
	print integrand_mbyl
	if R2=='A':
		fctr=2*math.pi/sum(integrand_mbyl)	
	else:
		fctr=math.pi/sum(integrand_mbyl)
	print "fctr for mthd2=%.5f"%fctr
		
for R2 in ['A','B','C','D']:
	print "*** %s ***"%R2
	print_integrand(R2)
	print "******"
print "fctr used so far=",1/math.radians(36)

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
	
