#!/usr/bin/python
from __future__ import division
import math
import os,sys
from collections import OrderedDict
import ROOT
from rootpy.plotting import Hist, HistStack, Legend, Canvas
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
import re

import elaslib

XMIN=14
XMAX=46
NBINS=32

BE=5.499
E1F_TARGET_RAD_LENGTH=0.00562
WCUT=1.028

fout=ROOT.TFile("TT.root","RECREATE")
hTTnorm=ROOT.TH1F("hTTnorm","hTTnorm",NBINS,XMIN,XMAX)
#hTTnorm.SetYTitle("d\sigma [#frac{\mu b}{sr}]")
hTTnorm_norad=ROOT.TH1F("hTTnorm_norad","hTTnorm_norad",NBINS,XMIN,XMAX)

for ibin in range(NBINS):
	theta=hTTnorm.GetBinLowEdge(ibin+1)
	binc=elaslib.elasrad(BE,theta,E1F_TARGET_RAD_LENGTH,WCUT)
	binc_norad=elaslib.elas(BE,theta)
	hTTnorm.SetBinContent(ibin+1,binc)
	hTTnorm_norad.SetBinContent(ibin+1,binc_norad)
	hTTnorm.SetBinError(ibin+1,0)
	hTTnorm_norad.SetBinError(ibin+1,0)

fout.Write()
