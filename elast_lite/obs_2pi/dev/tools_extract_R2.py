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

#! phi binning as per analysis
PHI_MIN,PHI_MAX=0,360
PHI_BINW=36
PHI_BIN_LEL=range(PHI_MIN,PHI_MAX,PHI_BINW)

def get_intgrll_R2f(R2):
	intgrll=[]
	for phi in PHI_BIN_LEL:
		phi_min=math.radians(phi)
		phi_max=math.radians(phi+PHI_BINW)
		if R2=='A':
			val=integrate.quad(lambda x: 1,phi_min,phi_max)
		elif R2=='B':
			val=integrate.quad(lambda x: math.cos(x),phi_min,phi_max)
		elif R2=='C':
			val=integrate.quad(lambda x: math.cos(2*x),phi_min,phi_max)
		elif R2=='D':
			val=integrate.quad(lambda x: math.sin(x),phi_min,phi_max)
		elif R2=='E':
			val=integrate.quad(lambda x: math.sin(2*x),phi_min,phi_max)
		intgrll.append(val[0])
	return intgrll

def get_mthd1_binned_EF(R2,intgrll):
	intgrl_mbyl=[]
	for ibin,phi in enumerate(PHI_BIN_LEL):
		if R2=='A':
			val=1*intgrll[ibin]
		elif R2=='B':
			val=math.cos(math.radians(phi))*intgrll[ibin]
		elif R2=='C':
			val=math.cos(2*math.radians(phi))*intgrll[ibin]
		elif R2=='D':
			val=math.sin(math.radians(phi))*intgrll[ibin]
		elif R2=='E':
			val=math.sin(2*math.radians(phi))*intgrll[ibin]
		intgrl_mbyl.append(val)
	return sum(intgrl_mbyl)

def get_mthd1_true_EF(R2):
	if R2=='A': #! 2*math.pi
		EF=integrate.quad(lambda x: 1,math.radians(PHI_MIN),math.radians(PHI_MAX))[0]
	elif R2=='B': #! math.pi
		EF=integrate.quad(lambda x: math.cos(x)*math.cos(x),math.radians(PHI_MIN),math.radians(PHI_MAX))[0]
	elif R2=='C': #! math.pi
		EF=integrate.quad(lambda x: math.cos(2*x)*math.cos(2*x),math.radians(PHI_MIN),math.radians(PHI_MAX))[0]
	elif R2=='D': #! math.pi
		EF=integrate.quad(lambda x: math.sin(x)*math.sin(x),math.radians(PHI_MIN),math.radians(PHI_MAX))[0]
 	elif R2=='E': #! math.pi
 		EF=integrate.quad(lambda x: math.sin(2*x)*math.sin(2*x),math.radians(PHI_MIN),math.radians(PHI_MAX))[0]
	return EF

def get_mthd2_CFl(R2,intgrll):
	'''
	Based on the fact that ROOT fits, unless specified, fit to bin content at the center of the bin
	'''
	CFl_m2=[]
	for ibin,phi in enumerate(PHI_BIN_LEL):
		if R2=='A':
			CFl_m2.append("%.5f"%(1/intgrll[ibin]))
		elif R2=='B':
			CFl_m2.append("%.5f"%(math.cos(math.radians(phi+(PHI_BINW/2)))/intgrll[ibin]))
		elif R2=='C':
			CFl_m2.append("%.5f"%(math.cos(2*math.radians(phi+(PHI_BINW/2)))/intgrll[ibin]))
		elif R2=='D':
			CFl_m2.append("%.5f"%(math.sin(math.radians(phi+(PHI_BINW/2)))/intgrll[ibin]))
		elif R2=='E':
			CFl_m2.append("%.5f"%(math.sin(2*math.radians(phi+(PHI_BINW/2)))/intgrll[ibin]))
	return CFl_m2
