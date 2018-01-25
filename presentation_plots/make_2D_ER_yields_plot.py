#!/usr/bin/python
from __future__ import division
import ROOT

import atlib as atlib
import q2w_bng

import collections
from array import *
import os,sys
import glob

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, '%s/elast_lite/obs_2pi'%os.environ['ANA2PI'])
from disp_obs import DispObs

from presentation_plots_lib import *

"""
+ Plots ER yields in a 2D Q2-W grid

+ >make_2D_ER_yields_plot.py
"""

#! Set up OBSDIR and SIMNUM
OBSDIR=[0 for i in range(NSIMRNG)]
OBSDIR[L]="%s/lowQ2_SSBands_122717/cutsncors1"%(os.environ['OBSDIR_E16'])
OBSDIR[H]="%s/highQ2_SSBands_122717/cutsncors1"%(os.environ['OBSDIR_E16'])
SIMNUM=[0 for i in range(NSIMRNG)]
SIMNUM[L]='sim4_sim5_sim6_sim7_sim8_sim13'
SIMNUM[H]='sim9_sim10_sim11_sim12'

#! OUTDIR
OUTDIR=os.path.join(os.environ['ANANOTE'],'figures','Results')
#OUTDIR="/tmp/ER_2D_yields"
if not os.path.exists(OUTDIR):
	os.makedirs(OUTDIR)

#! + The following binning should match up with h8_bng.h/q2w_bng.py
#! + It is slightly modified here for "aesthetic" purposes
WBINW=0.025
NWBINS,WMIN,WMAX=36,1.300,2.200# post removing W>2.125 GeV from analysis #68,1.3,3.0 #25 MeV/bin
NQ2BINS_TMP,Q2MIN_TMP,Q2MAX_TMP=8,Q2MIN,Q2MAX # place-holder bng, for it will be modified by variable bng
NBINS_Q2_VAR_BINW_BNG=5
Q2_BINS_VAR_BINW_BNG=array("d",[2.00,2.40,3.00,3.50,4.20,5.00])

#! Create histograms to display simstats(Q2,W)
h={}
for par in ['nbins','N','mu','sg']:
	h[par]=ROOT.TH2F("h_%s_%s"%(par,'ER'),"%s_%s(Q2,W)"%(par,'ER'),NWBINS,WMIN,WMAX,NQ2BINS_TMP,Q2MIN_TMP,Q2MAX_TMP)
	#! Set variable binning in Q2
	h[par].GetYaxis().Set(NBINS_Q2_VAR_BINW_BNG,Q2_BINS_VAR_BINW_BNG)
	h[par].SetXTitle("W [GeV]")
	h[par].SetYTitle("Q^{2} [GeV^{2}]")
		
#! Fill histograms
#! Prepare structure to get ER stats per sim range
ES=[0 for i in range(NSIMRNG)]
for isimrng in range(NSIMRNG):
	obsdir=OBSDIR[isimrng]
	simnum=SIMNUM[isimrng]
	dy=DispObs(obsdir,simnum=simnum)#,tops=topcmbn)
	ES[isimrng]=dy.get_ER_stats()
# print ES[L]
# sys.exit()
#! Now concatenate to get ES_TOT
ES_TOT=ES[L]+ES[H]

for d in ES_TOT:
	q2,w,nbins,N,mu,sg=d[0],d[1],d[2],d[3],d[4],d[5]
	binx=h['nbins'].GetXaxis().FindBin(w+(0.025/2));
	biny=h['nbins'].GetYaxis().FindBin(q2);
	print "w,q2,binw,binq2=",w,q2,binx,biny
	bin=h['nbins'].GetBin(binx,biny)
	h['nbins'].SetBinContent(bin,nbins)
	h['N'].SetBinContent(bin,N)
	h['mu'].SetBinContent(bin,mu)
	h['sg'].SetBinContent(bin,sg)			
			
#! Display "directly interesting" histograms to monitor simstats
ROOT.gStyle.SetOptStat(0)
#! N
c=ROOT.TCanvas("cN","cN")#,1000,1000)
h['N'].SetTitle("Experimental Yield(Q^{2},W)")
h['N'].Draw("colz")
#! Adjust right margin so that palette is not cut off
c.SetRightMargin(0.13)
#! Draw Q2-W grid
lw=[ROOT.TLine() for j in range(len(WBINL))]
lq=[ROOT.TLine() for j in range(len(Q2BINL))]
for j,wbin in enumerate(WBINL):
	lw[j].DrawLine(wbin,Q2MIN,wbin,Q2MAX)
for j,qbin in enumerate(Q2BINL):
	lq[j].DrawLine(WMIN,qbin,WMAX,qbin)
c.SaveAs("%s/ER_N_2D.png"%OUTDIR)
c.SaveAs("%s/ER_N_2D.pdf"%OUTDIR)
# c_pars=ROOT.TCanvas("cpars","cpars",1000,1000)
# c_pars.Divide(2,2)
# c_pars.cd(1)
# h_rtio_nbinsT_nbins=divide(h['ST','nbins'],h_nbins_vs_Q2W)
# #h_rtio_nbinsT_nbins.SetTitle("nbins_T(Q2,W)/nbins(Q2,W))")
# h_rtio_nbinsT_nbins.SetTitle("Fraction of max. possible PS-bins that are ST filled(Q^{2},W)")# for top %s"%''.join(str(t) for t in topcmbn))
# #h_rtio_nbinsT_nbins.GetZaxis().SetRangeUser(0,1)
# h_rtio_nbinsT_nbins.Draw("colz")
# c_pars.cd(2)
# h_rtio_nbinsH_nbinsT=divide(h['SH','nbins'],h['ST','nbins'])
# h_rtio_nbinsH_nbinsT.GetZaxis().SetRangeUser(0,1)
# #h_rtio_nbinsH_nbinsT.SetTitle("nbins_H(Q2,W)/nbins_T(Q2,W))")
# h_rtio_nbinsH_nbinsT.SetTitle("Fraction of T PS-bins that are Holes (Q^{2},W)")# for top %s"%''.join(str(t) for t in topcmbn))
# h_rtio_nbinsH_nbinsT.Draw("colz")
# c_pars.cd(3)
# h['SA','mu'].SetTitle("#LTSA/bin#GT(Q^{2},W)")# for top=%s"%''.join(str(t) for t in topcmbn))
# h['SA','mu'].Draw("colz")
# h['SA','mu'].GetZaxis().SetRangeUser(0,0.3)
# # if ''.join(str(t) for t in topcmbn)=="1234":
# # 	h['A','mu'].GetZaxis().SetRangeUser(0,0.3)
# # elif ''.join(str(t) for t in topcmbn)=="1":
# # 	h['A','mu'].GetZaxis().SetRangeUser(0,0.06)
# # else:
# # 	h['A','mu'].GetZaxis().SetRangeUser(0,0.2)#0.1
# c_pars.cd(4)
# h['SA','sg'].SetTitle("RMS_{SA/bin}(Q^{2},W)")#" for top=%s"%''.join(str(t) for t in topcmbn))
# h['SA','sg'].Draw("colz")
# h['SA','sg'].GetZaxis().SetRangeUser(0,0.3)
# # if ''.join(str(t) for t in topcmbn)=="1234":
# # 	h['A','sg'].GetZaxis().SetRangeUser(0,0.3)
# # else:
# # 	h['A','sg'].GetZaxis().SetRangeUser(0,0.1)
# # h_rtio_sgA_muA=divide(h['A','sg'],h['A','mu'])
# # h_rtio_sgA_muA.GetZaxis().SetRangeUser(0,2)
# # h_rtio_sgA_muA.SetTitle("#frac{RMS_{A}}{#LTA#GT}(Q^{2},W)")
# # h_rtio_sgA_muA.Draw("colz")
# c_pars.SaveAs("%s/pars1.png"%outdir)

# c_pars2=ROOT.TCanvas("cpars2","cpars2",1000,1000)
# c_pars2.Divide(2,2)
# c_pars2.cd(1)
# h['ST','mu'].SetTitle("#LTST/bin#GT(Q^{2},W)")#" for top %s"%''.join(str(t) for t in topcmbn))
# h['ST','mu'].Draw("colz")
# c_pars2.cd(2)
# h_rtio_sgT_muT=divide(h['ST','sg'],h['ST','mu'])
# h_rtio_sgT_muT.GetZaxis().SetRangeUser(0,2)
# h_rtio_sgT_muT.SetTitle("#frac{RMS_{ST/bin}}{#LTST/bin#GT}(Q^{2},W)")# for top %s"%''.join(str(t) for t in topcmbn))
# h_rtio_sgT_muT.Draw("colz")
# c_pars2.cd(3)
# h['SR','mu'].SetTitle("#LTSR/bin#GT(Q^{2},W)")#" for top %s"%''.join(str(t) for t in topcmbn))
# h['SR','mu'].Draw("colz")
# c_pars2.cd(4)
# h_rtio_sgR_muR=divide(h['SR','sg'],h['SR','mu'])
# h_rtio_sgR_muR.GetZaxis().SetRangeUser(0,2)
# h_rtio_sgR_muR.SetTitle("#frac{RMS_{SR/bin}}{#LTSR/bin#GT}(Q^{2},W)")#" for top %s"%''.join(str(t) for t in topcmbn))
# h_rtio_sgR_muR.Draw("colz")
# c_pars2.SaveAs("%s/pars2.png"%outdir)

# c_pars3=ROOT.TCanvas("cpars3","cpars3",1000,1000)
# c_pars3.Divide(2,1)
# c_pars3.cd(1)
# h['ST','N'].SetTitle("N-ST(Q^{2},W)")#" for top %s"%''.join(str(t) for t in topcmbn))
# h['ST','N'].Draw("colz")
# c_pars3.cd(2)
# h['SR','N'].SetTitle("N-SR(Q^{2},W)")# for top %s"%''.join(str(t) for t in topcmbn))
# h['SR','N'].Draw("colz")
# c_pars3.SaveAs("%s/pars3.png"%outdir)


# #! Save all the histograms to file
# #f=ROOT.TFile("%s/simstats_top%s.root"%(outdir,''.join(str(t) for t in topcmbn)),"RECREATE")
# f=ROOT.TFile("%s/simstats.root"%outdir,"RECREATE")
# for seq in ['ST','SR','SA','SH']:
# 	for par in ['nbins','N','mu','sg']:
# 		h[seq,par].Write()
# h_nbins_vs_Q2W.Write()
# h_rtio_nbinsT_nbins.Write()
# h_rtio_nbinsH_nbinsT.Write()
# #h_rtio_sgA_muA.Write()
# h_rtio_sgT_muT.Write()
# h_rtio_sgR_muR.Write()
# f.Close()
