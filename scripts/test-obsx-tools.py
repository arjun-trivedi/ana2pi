from __future__ import division
#python imports
import os
from array import array
import math
#scipy imports
from matplotlib import pyplot as plt
import numpy as np
#rootpy imports
from rootpy.io import root_open, File
import rootpy.plotting.root2matplotlib as rplt
#import root2matplot as r2m
from rootpy.plotting import Hist

from ROOT import gSystem, gROOT, THnSparseF, TCanvas, TString, TLine, TPaveText, TText
gSystem.Load('myTHnTool_C')
from ROOT import myTHnTool
mythnt = myTHnTool()


"""Test Case for plotR2"""
var = {'M1':0, 'M2':1, 'THETA':2, 'PHI':3, 'ALPHA':4};
print var


#INPUT data
anadir =  os.environ['E1F_2PI_ANADIR2']
fname = os.path.join(anadir,'1:2:3:4__1-2.000-2.400__24-1.300-1.900__pol__exp.root')#test.root')
f = root_open(fname)

"""rootpy"""
# for dirs in f.walk('.').next()[1]:
# 	h5 = {}
# 	print dirs
# 	obj = os.path.join('.',dirs,'hY5D_POS','Varset1','hY5D_ACC_CORR')
# 	print 'obj=',obj
# 	h5[('EC','POS')] = obj
# 	mythnt.MultiplyBySinPhi(h5[('EC','POS')]);

"""PyRoot"""
keys = f.GetListOfKeys()
#OUTPUT data
outdir_root = os.path.join(anadir,'polobs.new')
if not os.path.isdir(outdir_root):
	os.makedirs(outdir_root)

for q2wdir in keys:
	outdir = os.path.join( ('%s/%s/')%(outdir_root,q2wdir.GetName()) )
	if not os.path.isdir(outdir):
		os.makedirs(outdir);

	h5 = {}	
	h5[('EC','POS')]=f.Get('%s/hY5D_POS/Varset1/hY5D_ACC_CORR'%
							q2wdir.GetName());
	h5[('EC','NEG')]=f.Get('%s/hY5D_NEG/Varset1/hY5D_ACC_CORR'%
							q2wdir.GetName());
	mythnt.MultiplyBySinPhi(h5[('EC','POS')]);
	mythnt.MultiplyBySinPhi(h5[('EC','NEG')],-1);
	
	norm = 50000
	hR2 = {}
	hR2[('THETA','POS')] = h5[('EC','POS')].Projection(var['THETA'])
	hR2[('THETA','POS')].Scale(1/math.pi)
	hR2[('THETA','POS')].Scale(1/norm)
	hR2[('THETA','NEG')] = h5[('EC','NEG')].Projection(var['THETA'])
	hR2[('THETA','NEG')].Scale(1/math.pi)
	hR2[('THETA','NEG')].Scale(1/norm)
	hR2[('THETA','AVG')] = hR2[('THETA','POS')].Clone("avg")
	hR2[('THETA','AVG')].Add(hR2[('THETA','NEG')])
	hR2[('THETA','AVG')].Scale(0.5)
	hR2[('THETA','AVG')].SetMinimum(-0.003)
	hR2[('THETA','AVG')].SetMaximum(0.003)
	hR2[('THETA','AVG')].SetLineColor(gROOT.ProcessLine("kMagenta"));
	hR2[('THETA','AVG')].SetMarkerStyle(gROOT.ProcessLine("kFullCircle"));
	#Make Titles nice
	hR2[('THETA','AVG')].SetTitle("")
	pt = TPaveText(0.3, 0.85, 0.7, 1.0, "NDC")
	q2wt = pt.AddText('[Q^{2}][W] = %s'%q2wdir.GetName())
	q2wt.SetTextColor(gROOT.ProcessLine("kBlue"))
	vart = pt.AddText(("D^{%s} vs. %s")%('test','test'));
	vart.SetTextSize(0.05);

	
	cR2 = {}
	l = TLine(0,0,180,0)
	cR2[('THETA','AVG')] = TCanvas("RvVar", "RvVar")
	hR2[('THETA','AVG')].Draw("ep")
	l.Draw("same")
	pt.Draw()
	cR2[('THETA','AVG')].SaveAs( ('%s/%s.png')%(outdir,cR2[('THETA','AVG')].GetName()))