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

from ROOT import gSystem, gROOT, gStyle, THnSparseF, TCanvas, TString, TLine, TPaveText, TText
gSystem.Load('myTHnTool_C')
from ROOT import myTHnTool
mythnt = myTHnTool()

"""Test Case for plotR2"""
nVARSETS=3

nVARS=5
M1,M2,THETA,PHI,ALPHA=range(nVARS)

nSEQ_SIM=6
ST,SR,SA,SC,SH,SF=range(nSEQ_SIM)

nSEQ_EXP=3
ER,EC,EF=range(nSEQ_EXP)

nPOLS=4
POS,NEG,UNP,AVG=range(nPOLS)

nPOBS=4
A,B,C,D=range(nPOBS) # [A,B,C,D]=[Rt+Rl,Rlt,Rtt,Rlt']

#print var
varname = ['M1','M2','theta','phi','alpha']
vartitle = [ 
				["M_{p#pi^{+}}", "M_{#pi^{+}#pi^{-}}","#theta_{#pi^{-}}", "#phi_{#pi^{-}}", "#alpha_{[p^{'}#pi^{+}][p#pi^{-}]}"],
			 	["M_{p#pi^{+}}", "M_{#pi^{+}#pi^{-}}","#theta_{p}", "#phi_{p}", "#alpha_{[#pi^{+}#pi^{-}][p^{'}p]"],
			 	["M_{p#pi^{+}}", "M_{p#pi^{-}}", "#theta_{#pi^{+}}", "#phi_{#pi^{+}}", "#alpha_{[p^{'}#pi^{-}][p#pi^{+}]"] 
		   ]
def plotR2_D(h5):
	norm = 50000
	hR2 = {} #indexed by (ij,alpha,pol)
	cR2 = {}
	for var in range(0,nVARS):
		if var==PHI or var==ALPHA:
			continue
		hR2[(var,POS,D)] = h5[(EC,POS,D)].Projection(var)
		hR2[(var,POS,D)].Scale(1/math.pi)
		hR2[(var,POS,D)].Scale(1/norm)
		hR2[(var,NEG,D)] = h5[(EC,NEG,D)].Projection(var)
		hR2[(var,NEG,D)].Scale(1/math.pi)
		hR2[(var,NEG,D)].Scale(1/norm)
		hR2[(var,AVG,D)] = hR2[(var,POS,D)].Clone("avg")
		hR2[(var,AVG,D)].Add(hR2[(var,NEG,D)])
		hR2[(var,AVG,D)].Scale(0.5)
		hR2[(var,AVG,D)].SetMinimum(-0.003)
		hR2[(var,AVG,D)].SetMaximum(0.003)
		hR2[(var,AVG,D)].SetLineColor(gROOT.ProcessLine("kMagenta"));
		hR2[(var,AVG,D)].SetMarkerStyle(gROOT.ProcessLine("kFullCircle"));
		#Make Titles nice
		hR2[(var,AVG,D)].SetTitle("")
		pt = TPaveText(0.3, 0.85, 0.7, 1.0, "NDC")
		q2wt = pt.AddText('[Q^{2}][W] = %s'%q2wdir.GetName())
		q2wt.SetTextColor(gROOT.ProcessLine("kBlue"))
		vart = pt.AddText(("D^{%s} vs. %s")%(vartitle[0][var],vartitle[0][var]));
		vart.SetTextSize(0.05);

	
		l = TLine(0,0,180,0)
		cname = ('R2%sV%s')%(varname[var],varname[var])
		cR2[(var,AVG,D)] = TCanvas(cname,cname)#"RvVar", "RvVar")
		hR2[(var,AVG,D)].Draw("ep")
		l.Draw("same")
		pt.Draw()
		csavename = ('%s/%s')%(outdir,cR2[(var,AVG, D)].GetName())
		cR2[(var,AVG, D)].SaveAs( ('%s.png')%(csavename))


#COSMETICS
gStyle.SetOptStat(0);
gStyle.SetOptFit(1111);
#INPUT data
anadir =  os.environ['E1F_2PI_ANADIR2']
fname = os.path.join(anadir,'1:2:3:4__1-2.000-2.400__24-1.300-1.900__pol__exp.root')#test.root')
f = root_open(fname)
#OUTPUT data
outdir_root = os.path.join(anadir,'polobs.new')
if not os.path.isdir(outdir_root):
	os.makedirs(outdir_root)

#MAIN part of program
# Using rootpy
# for dirs in f.walk('.').next()[1]:
# 	h5 = {}
# 	print dirs
# 	obj = os.path.join('.',dirs,'hY5D_POS','Varset1','hY5D_ACC_CORR')
# 	print 'obj=',obj
# 	h5[(EC,POS)] = obj
# 	mythnt.MultiplyBySinPhi(h5[(EC,POS)]);

#Using PyRoot
"""
First, for each q2wbin, make all R2^{ij}_{alpha} where:
	-ij    = Index over nVARSETS and nVAR respectively
	-alpha = Index over nPOBS (=[A,B,C,D])
"""
keys = f.GetListOfKeys()
for q2wdir in keys:
	outdir = os.path.join( ('%s/%s/')%(outdir_root,q2wdir.GetName()) )
	if not os.path.isdir(outdir):
		os.makedirs(outdir);

	h5 = {}	
	h5[(EC,POS)]=f.Get('%s/hY5D_POS/Varset1/hY5D_ACC_CORR'%
							q2wdir.GetName());
	h5[(EC,NEG)]=f.Get('%s/hY5D_NEG/Varset1/hY5D_ACC_CORR'%
							q2wdir.GetName());
	h5[(EC,POS,D)] = mythnt.MultiplyBySinPhi(h5[(EC,POS)]);
	h5[(EC,NEG,D)] = mythnt.MultiplyBySinPhi(h5[(EC,NEG)],-1);
	print 'h5=',h5

	plotR2_D(h5)
	
"""
Now put all q2wbin/R2^{ij}_{alpha} in a single pdf
"""