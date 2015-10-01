#!/usr/bin/python
from __future__ import division

import os,sys,time
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict

import numpy as np

#! Get args from user
PBINW=0.100 #! 100 MeV
if len(sys.argv)>=2:
	PBINW=float(sys.argv[1])

NENTRIES=10000000000
if len(sys.argv)>=3:
	NENTRIES=int(sys.argv[2])
print "PBINW,NENTRIES=",PBINW,"GeV",NENTRIES

#! Initial set up
#! setup momentum binning
PMIN=0.64
PMAX=5.00
PBIN_LE=np.arange(PMIN,PMAX,PBINW)
NPBIN=len(PBIN_LE)
print("NPBIN=",NPBIN)
time.sleep(3)

NDTYP=2
EXP,SIM=range(2)
DTYP_NAME=['exp','sim']

NSCTR=6

def make_fout_SFvp():
	FIN_2pi=[0,0]
        T_2pi=[0,0]
        FIN_2pi[EXP]=ROOT.TFile("%s/data_SF_100115/dSF_2pi.root"%os.environ['D2PIDIR_EXP'])
        FIN_2pi[SIM]=ROOT.TFile("%s/data_SF_100115/dSF_2pi.root"%os.environ['D2PIDIR_SIM'])
        T_2pi[EXP]=FIN_2pi[EXP].Get("d2piR/tR")
        T_2pi[SIM]=FIN_2pi[SIM].Get("d2piR/tR")
	
	FIN_elast=[0,0]
        T_elast=[0,0]
        FIN_elast[EXP]=ROOT.TFile("%s/data_SF_100115/dSF_elast.root"%os.environ['D2PIDIR_EXP'])
        FIN_elast[SIM]=ROOT.TFile("%s/data_SF_100115/dSF_elast.root"%os.environ['D2PIDIR_SIM'])
        T_elast[EXP]=FIN_elast[EXP].Get("delast/cut/t")
        T_elast[SIM]=FIN_elast[SIM].Get("delast/cut/t")

        #! Make hSFvp[NDTYP][NSCTR]
        draw_cmd="etot/p:p>>hdcmd(100,0,5,100,0,0.5)"
        hSFvp=[[[] for isctr in range(NSCTR)] for idtyp in range(NDTYP)]
        for idtyp in range(NDTYP):
                for isctr in range(NSCTR):
			print "Making SFvp for",DTYP_NAME[idtyp],isctr+1
                        cut_sector=ROOT.TCut("sector==%d"%(isctr+1))
			#! First draw hist from 2pi data
                        T_2pi[idtyp].Draw(draw_cmd,cut_sector,"colz",NENTRIES)
                        #! put hist in hSFvp[NDTYP][NSCTR]
                        htmp=ROOT.gDirectory.Get("hdcmd")
                        hname="h_%s_s%d"%(DTYP_NAME[idtyp],isctr+1)
                        hSFvp[idtyp][isctr]=htmp.Clone(hname)

			#! Now draw hist from elast data
			T_elast[idtyp].Draw(draw_cmd,cut_sector,"colz",NENTRIES)
			#! Add hist to already created hSFvp[NDTYP][NSCTR]
			htmp=ROOT.gDirectory.Get("hdcmd")
			hSFvp[idtyp][isctr].Add(htmp)
			

        #! Plot and write hSFvp[NDTYP][NSCTR] to file
        ROOT.gStyle.SetOptStat("ne")
        C=[0,0]
        FOUT_SFvp=ROOT.TFile(fout_sfvp_name,"RECREATE")
        for idtyp in range(NDTYP):
                C[idtyp]=ROOT.TCanvas("c_%s"%DTYP_NAME[idtyp],"c_%s"%DTYP_NAME[idtyp])
                C[idtyp].Divide(3,2)
                for isctr in range(NSCTR):
                        C[idtyp].cd(isctr+1)
                        hSFvp[idtyp][isctr].Draw("colz")
                        hSFvp[idtyp][isctr].Write()
        C[EXP].Write()
        C[SIM].Write()
        FOUT_SFvp.Close()

def plot_SF():
        #! Create, plot and save: hSF[NDTYP][NSCTR][NPIN]
	FOUT_SFvp=ROOT.TFile(fout_sfvp_name)
        hSF=[[[[] for ipbin in range(NPBIN)] for isctr in range(NSCTR)] for idtyp in range(NDTYP)]
        for idtyp in range(NDTYP):
                for isctr in range(NSCTR):
                        hSFvp=FOUT_SFvp.Get("h_%s_%d"%(DTYP_NAME[idtyp],isctr+1))
                        outdir="%s/SFvp/%s/s%d"%(os.environ['STUDY_EID_SF'],DTYP_NAME[idtyp],isctr+1)
                        if not os.path.exists(outdir):
                                os.makedirs(outdir)
                        for ipbin in range(NPBIN):
                                hname="%s_s%d_pbin%d"%(DTYP_NAME[idtyp],isctr+1,ipbin+1)
                                htitle="p=[%.2f,%.2f)"%(PBIN_LE[ipbin],PBIN_LE[ipbin]+PBINW)
                                print "Making SF projection for",hname
                                bin1=hSFvp.GetXaxis().FindBin(PBIN_LE[ipbin])
                                bin2=hSFvp.GetXaxis().FindBin(PBIN_LE[ipbin]+PBINW)
                                hSF[idtyp][isctr][ipbin]=hSFvp.ProjectionY(hname,bin1,bin2)
                                hSF[idtyp][isctr][ipbin].SetTitle(htitle)
                                c=ROOT.TCanvas()
                                hSF[idtyp][isctr][ipbin].Draw()
                                c.SaveAs("%s/c_pbin%02d.png"%(outdir,ipbin+1))

#! Start of main program
fout_sfvp_name="%s/fSFvp.root"%os.environ['STUDY_EID_SF']
if not os.path.exists(fout_sfvp_name):
	make_fout_SFvp()
	plot_SF()
else:
        print fout_sfvp_name,"exists"
	plot_SF()
	
if not ROOT.gROOT.IsBatch():
	plt.show()
	# wait for you to close the ROOT canvas before exiting
	wait(True)

#if __name__ == "__main__":	
#	if len(sys.argv)==2:
#		plot_fid(sys.argv[1])
#	elif len(sys.argv)==3:
#		plot_fid(sys.argv[1],int(sys.argv[2])) 
	
