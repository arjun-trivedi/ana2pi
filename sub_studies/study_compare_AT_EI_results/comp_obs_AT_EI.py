#!/usr/bin/python
from __future__ import division

import os,sys,datetime
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict
import itertools
from array import *

import numpy as np
import pandas as pd

import math

#General global vars
SYST_ERR_EI=0.14 #14% taken from his analysis note

DATE=datetime.datetime.now().strftime('%m%d%y')
OUTDIR=os.path.join(os.environ['STUDY_COMPARE_AT_EI_RESULTS_DATADIR'],'results_%s'%DATE)
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

NDTYP=2
AT,EI=range(NDTYP)
DTYP_NAME=["AT","EI"]
DTYP_TITLE=["This analysis","E. Isupov"]

INDIR=[0 for i in range(NDTYP)]
#INDIR[AT]='%s/thesis_obs_norm_ST_shape_080317'%os.environ['OBSDIR_E162']
INDIR[AT]='%s/thesis_obs_norm_ST_shape_122817'%os.environ['OBSDIR_E162']
INDIR[EI]='%s/EI_results'%os.environ['OBSDIR_E162']

CLRS=[ROOT.gROOT.ProcessLine("kBlue"),ROOT.gROOT.ProcessLine("kBlack")]
MRKS=[ROOT.gROOT.ProcessLine("kFullCircle"),ROOT.gROOT.ProcessLine("kFullCircle")]

#! The following is needed to make h1D for EI
MASS_P=0.93827203
MASS_PIP=0.13957018
MASS_PIM=0.13957018
#! + Coded in as per PAD_MAP
#! + Differences from AT binning:
#!    + M:diff: nbins AT/EI=14/10
#!    + ALPHA:diff: nbins AT/EI=10/6
BNG_EI={(1,"M1"):[10,MASS_P+MASS_PIP,'NA'],    #! 'mpippr'
        (3,'M2'):[10,MASS_P+MASS_PIM,'NA'],    #! 'mpimpr' 
        (2,'M2'):[10,MASS_PIP+MASS_PIM,'NA'], #! 'mpippim'
        (1,"THETA"):[10,0,180],
        (3,'THETA'):[10,0,180],
        (2,'THETA'):[10,0,180],
        (1,"ALPHA"):[6,0,360],
        (3,'ALPHA'):[6,0,360], 
        (2,'ALPHA'):[6,0,360]}
# #! + M (diff: nbins AT=14)
# MASS_P=0.93827203
# MASS_PIP=0.13957018
# MASS_PIM=0.13957018
# NBINS_M_EI=10 #! Note that this binning may change in EI's data. This is as obtained from q2_24_30w_1725_1750 and other bins may be different

# #! Note that XMAX limits are not set because depend on W and will be determined when processing q2wbin
# XMIN_1_M1=MASS_P+MASS_PIP   #! Mppip (EI:mpippr)
# XMIN_3_M2=MASS_P+MASS_PIM   #! Mppim (EI:mpimpr)
# XMIN_2_M2=MASS_PIP+MASS_PIM #! Mpippim (EI:mpippim)
# #! + THETA (diff: none)
# NBINS_THETA_EI=10
# XMIN_THETA_EI=0
# XMAX_THETA_EI=180
# #! + ALPHA (diff: nbins AT=10)
# NBINS_ALPHA_EI=6
# XMIN_ALPHA_EI=0
# XMAX_ALPHA_EI=360

def compare_itg():
    outdir="%s/itg"%OUTDIR
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    #! Store Q2 bin names for AT and EI in one structure
    Q2L=[["2.00-2.40_1.400-2.125","cs20_24"],["2.40-3.00_1.400-2.125","cs24_30"],["3.00-3.50_1.400-2.125","cs30_35"],
         ["3.50-4.20_1.400-2.125","cs35_42"],["4.20-5.00_1.400-2.125","cs42_50"]]

    #! Create structure to store data, in the 2nd and 3rd Res. region, to plot sgW vs Q2 to see Q2 evolution 
    #! R2 W=1.500-1.525, 5th wbin
    #! R3 W=1.725-1.750, 14th wbin
    R2WBIN=5
    R3WBIN=14
    R2WBIN_TITLE="[1.500, 1.525)"
    R3WBIN_TITLE="[1.725, 1.750)"
    sgR2=[[] for i in range(NDTYP)]
    erR2=[[] for i in range(NDTYP)]
    sgR3=[[] for i in range(NDTYP)]
    erR3=[[] for i in range(NDTYP)]
    #! Store Q2 information: center,bin-width
    #! + Calculated from Q2L
    Q2C=[2.20,2.70,3.25,3.85,4.60]
    Q2W=[0.20,0.30,0.25,0.35,0.40]

    for iq2,q2n in enumerate(Q2L):
        q2=[0 for i in range(NDTYP)]
        q2[AT]=q2n[0]
        q2[EI]=q2n[1]

        #! structure for storing W and W-err: hobs, herr
        hobs=[0 for i in range(NDTYP)]
        herr=[0 for i in range(NDTYP)]

        #! get AT data
        if iq2+1<=2:
            fin=ROOT.TFile("%s/lowQ2/Obs_itg/Obs_itg.root"%INDIR[AT],"READ")
        else:
            fin=ROOT.TFile("%s/highQ2/Obs_itg/Obs_itg.root"%INDIR[AT],"READ")
        print '%s/hW_obs_EF_THETA'%q2[AT]
        hobs[AT]=fin.Get('%s/hW_obs_EF_THETA'%q2[AT])
        herr[AT]=fin.Get('%s/hW_err_EF_THETA'%q2[AT])
        print fin.GetName()
        print hobs[AT].GetName()
        #! Add syst-err to stat-err
        nbins=hobs[AT].GetNbinsX()
        for ibin in range(nbins):
            stat_err=hobs[AT].GetBinError(ibin+1)
            syst_err=herr[AT].GetBinContent(ibin+1)
            #print stat_err, syst_err
            #! calculate total error and use that in hat
            tot_err=math.sqrt(stat_err**2+syst_err**2)
            hobs[AT].SetBinError(ibin+1,tot_err)
        #! get sgR2/W2 information
        sgR2[AT].append(hobs[AT].GetBinContent(R2WBIN))
        sgR3[AT].append(hobs[AT].GetBinContent(R3WBIN))
        erR2[AT].append(hobs[AT].GetBinError(R2WBIN))
        erR3[AT].append(hobs[AT].GetBinError(R3WBIN))
    
        #! get EI-data
        #! + Get in df version (tried reading line by line and then parsing, but there values are 
        #!   separated by uneven number of space-delmiter)
        #! + Code source: #! https://stackoverflow.com/questions/43058462/reading-a-variable-white-space-delimited-table-in-python
        df=pd.read_csv('%s/%s.txt'%(INDIR[EI],q2[EI]),sep='\s{1,}',header=None,engine='python',thousands=',')
        #! Get number of rows (corresponds to number of wbins)
        nrows=len(df.index)
        #print df.head()
        #! create hei from hat.Clone and Reset
        hobs[EI]=hobs[AT].Clone("hei")
        hobs[EI].Reset()
        hobs[EI].SetMarkerColor(1)
        hobs[EI].SetLineColor(1)
        #hobs[EI].SetMarkerStyle(2)
        for ibin in range(nrows):
            wbin,val,stat_err=df.loc[ibin,0],df.loc[ibin,1],df.loc[ibin,2]
            #print ibin+1,wbin,val,stat_err
            #! Calculate total error
            tot_err=math.sqrt(stat_err**2+(SYST_ERR_EI*val)**2)
            #! Now fill histogram
            hobs[EI].SetBinContent(ibin+1,val)
            hobs[EI].SetBinError(ibin+1,tot_err)
            #! get sgR2/W2 information
            if ibin+1==R2WBIN:  
                sgR2[EI].append(val)
                erR2[EI].append(tot_err)
            if ibin+1==R3WBIN:  
                sgR3[EI].append(val)
                erR3[EI].append(tot_err)

        #! fine tuning
        #! +ymax
        if   q2[AT]=="2.40-3.00_1.400-2.125":ymax=8.2
        elif q2[AT]=="3.00-3.50_1.400-2.125":ymax=5.8
        elif q2[AT]=="3.50-4.20_1.400-2.125":ymax=4.0
        if q2[AT]=="2.40-3.00_1.400-2.125" or q2[AT]=="3.00-3.50_1.400-2.125" or q2[AT]=="3.50-4.20_1.400-2.125":
            hobs[AT].SetMaximum(ymax)

        #! Draw
        #! Histogram aesthetics
        #! Reset any previous settings
        ROOT.gStyle.Reset()
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPadBottomMargin(0.20)
        #!get rid of X error bars
        ROOT.gStyle.SetErrorX(0.001)

        c=ROOT.TCanvas()
        hobs[AT].Draw()
        hobs[EI].Draw("same")
        l=ROOT.TLegend(0.60,0.80,0.90,0.90)
        l.SetBorderSize(1)
        l.SetTextSize(0.04)
        l.AddEntry(hobs[AT],"%s"%DTYP_TITLE[AT],"p")
        l.AddEntry(hobs[EI],"%s"%DTYP_TITLE[EI],"p")
        l.Draw()
        c.SaveAs("%s/c%s.png"%(outdir,q2[AT]))
        c.SaveAs("%s/c%s.pdf"%(outdir,q2[AT]))

    #! Compare W vs Q2 evolution
    outdir_itg_wvsq="%s/sgvsQ2"%outdir
    if not os.path.exists(outdir_itg_wvsq):
        os.makedirs(outdir_itg_wvsq)

    NQ2BINS=5
    Q2BINS=array("d",[2.00,2.40,3.00,3.50,4.20,5.00])

    #! Create histograms
    hR2=[0 for i in range(NDTYP)]
    hR3=[0 for i in range(NDTYP)]
    for i in range(NDTYP):
        #! R2
        hR2[i]=ROOT.TH1F("hR2_%s"%DTYP_NAME[i],"hR2_%s"%DTYP_NAME[i],NQ2BINS,Q2BINS)
        hR2[i].SetTitle("Q^{2} Evolution of the Integrated Cross-Section for W Bin=%s GeV"%R2WBIN_TITLE)
        #! R3
        hR3[i]=ROOT.TH1F("hR3_%s"%DTYP_NAME[i],"hR3_%s"%DTYP_NAME[i],NQ2BINS,Q2BINS)
        hR3[i].SetTitle("Q^{2} Evolution of the Integrated Cross-Section for W Bin=%s GeV"%R3WBIN_TITLE)
        for hist in [hR2[i],hR3[i]]:
            hist.SetXTitle("Q^{2} (GeV^{2})")
            hist.SetYTitle("#sigma (#mub)")
            hist.SetLineColor(CLRS[i])
            hist.SetMarkerColor(CLRS[i])
            hist.SetMarkerStyle(MRKS[i])

    #! Fill histogrms
    for i in range(NDTYP):
        #! R2
        for ibin,sg in enumerate(sgR2[i]):
            hR2[i].SetBinContent(ibin+1,sg)
            hR2[i].SetBinError(ibin+1,erR2[i][ibin])
        #! R3
        for ibin,sg in enumerate(sgR3[i]):
            hR3[i].SetBinContent(ibin+1,sg)
            hR3[i].SetBinError(ibin+1,erR3[i][ibin])
    
    #! Define Dipole fit function 
    ffR2=[0 for i in range(NDTYP)]
    ffR3=[0 for i in range(NDTYP)]
    for i in range(NDTYP):
        #! R2
        ffR2[i]=ROOT.TF1("ffR2_%s"%DTYP_NAME[i],"[0]*x**(-2)",2,5)
        ffR2[i].SetLineColor(CLRS[i])
        #! R3
        ffR3[i]=ROOT.TF1("ffR2_%s"%DTYP_NAME[i],"[0]*x**(-2)",2,5)
        ffR3[i].SetLineColor(CLRS[i])


    #! Do integral-fit
    for i in range(NDTYP):
        #! R2
        hR2[i].Fit(ffR2[i].GetName(),"+0I")
        #! R3
        hR3[i].Fit(ffR3[i].GetName(),"+0I")

    #! fine-tuning
    #! yaxis
    #hR2[AT].SetMaximum(10)
    hR2[AT].SetMinimum(0)
    hR3[AT].SetMaximum(10)
    hR3[AT].SetMinimum(0)

    #! Draw
    #! Histogram aesthetics
    #! Reset any previous settings
    ROOT.gStyle.Reset()
    ROOT.gStyle.SetOptStat(0)
    #ROOT.gStyle.SetTitleFontSize(1)
    #ROOT.gStyle.SetTitleW(5)
    #!get rid of X error bars
    ROOT.gStyle.SetErrorX(0.001)

    #! R2
    c=ROOT.TCanvas()
    hR2[AT].Draw()
    hR2[EI].Draw("same")
    fat=hR2[AT].GetFunction(ffR2[AT].GetName())
    fat.Draw("same")
    fei=hR2[EI].GetFunction(ffR2[EI].GetName())
    fei.Draw("same")
    l=ROOT.TLegend(0.60,0.80,0.90,0.90)
    l.SetBorderSize(1)
    l.SetTextSize(0.04)
    l.AddEntry(hR2[AT],"%s"%DTYP_TITLE[AT],"p")
    l.AddEntry(hR2[EI],"%s"%DTYP_TITLE[EI],"p")
    l.Draw()
    c.SaveAs("%s/R2.png"%outdir_itg_wvsq)
    c.SaveAs("%s/R2.pdf"%outdir_itg_wvsq)
    #! R3
    c=ROOT.TCanvas()
    hR3[AT].Draw()
    hR3[EI].Draw("same")
    fat=hR3[AT].GetFunction(ffR3[AT].GetName())
    fat.Draw("same")
    fei=hR3[EI].GetFunction(ffR3[EI].GetName())
    fei.Draw("same")
    l=ROOT.TLegend(0.60,0.80,0.90,0.90)
    l.SetBorderSize(1)
    l.SetTextSize(0.04)
    l.AddEntry(hR3[AT],"%s"%DTYP_TITLE[AT],"p")
    l.AddEntry(hR3[EI],"%s"%DTYP_TITLE[EI],"p")
    l.Draw()
    c.SaveAs("%s/R3.png"%outdir_itg_wvsq)
    c.SaveAs("%s/R3.pdf"%outdir_itg_wvsq)
    

def compare_1D():

    outdir="%s/1D"%OUTDIR
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    #! Note PAD_MAP extended to include explicit variable name (as per EI's syntax)
    PAD_MAP=[(1,1,"M1","mpippr"),      (2,3,'M2',"mpimpr"),      (3,2,'M2','mpippim'),
	        (4,1,"THETA","theta_pim"),(5,3,'THETA','theta_pip'),(6,2,'THETA','theta_pr'),
		   (7,1,"ALPHA",'psi_pim'),  (8,3,'ALPHA','psi_pip'),    (9,2,'ALPHA','psi_pr')]

    #! Set Q2WBINS in which comaparison is to be made
    #! + AT and EI have different naming convention
    Q2WBIN_DIRS=[[] for i in range(NDTYP)]
    Q2WBIN_DIRS[AT]=["2.00-2.40_1.500-1.525","2.00-2.40_1.725-1.750","2.40-3.00_1.500-1.525","2.40-3.00_1.725-1.750"]
    Q2WBIN_DIRS[EI]=["q2_20_24/w_1500_1525", "q2_20_24/w_1725_1750", "q2_24_30/w_1500_1525", "q2_24_30/w_1725_1750"]#["q2_24_30w_1725_1750"]

    #! Create structures to hold hobs,herr for Obs_1D
    hobs=[0 for i in range(NDTYP)]
    herr=[0 for i in range(NDTYP)]

    for item in zip(Q2WBIN_DIRS[AT],Q2WBIN_DIRS[EI]):
        q2wbin_dir=[0 for i in range(NDTYP)]
        q2wbin_dir[AT]=item[AT]
        q2wbin_dir[EI]=item[EI]
        #! get wmin,wmax from q2wbin (need for setting up M dimension for EI:1D)
        wmin=float(q2wbin_dir[AT].split("_")[1].split("-")[0])
        wmax=float(q2wbin_dir[AT].split("_")[1].split("-")[1])
        
        #! Get AT data
        f=ROOT.TFile("%s/lowQ2/Obs_1D/Obs_1D.root"%INDIR[AT],"READ")
        hobs[AT]=OrderedDict()
        herr[AT]=OrderedDict()
        for item in PAD_MAP:
            pad,vst,var=item[0],item[1],item[2]
            #print pad,vst,var
            hobs[AT][vst,var]=f.Get("%s/h1_EF_%d_%s"%(q2wbin_dir[AT],vst,var))
            herr[AT][vst,var]=f.Get("%s/herr_h1_EF_%d_%s"%(q2wbin_dir[AT],vst,var))
            hobs[AT][vst,var].SetMarkerColor(CLRS[AT])
            hobs[AT][vst,var].SetLineColor(CLRS[AT])
        #! Add syst-err to stat-err
        for k in hobs[AT]:
            nbins=hobs[AT][k].GetNbinsX()
            for ibin in range(nbins):
                stat_err=hobs[AT][k].GetBinError(ibin+1)
                syst_err=herr[AT][k].GetBinContent(ibin+1)
                #print stat_err, syst_err
                #! calculate total error and use that in hat
                tot_err=math.sqrt(stat_err**2+syst_err**2)
                hobs[AT][k].SetBinError(ibin+1,tot_err)
    
        #! Get EI data
        hobs[EI]=OrderedDict()
        for item in PAD_MAP:
            pad,vst,var,varname=item[0],item[1],item[2],item[3]
            df=pd.read_csv('%s/%s/xsec1d_%s.dat'%(INDIR[EI],q2wbin_dir[EI],varname),sep='\s{1,}',header=None,engine='python',thousands=',')
            #print df.head()
            #! Setup binning for hobs(vst,var) for EI
            hname="h_%d_%s"%(vst,var)
            nbins,xmin,xmax=BNG_EI[vst,var][0],BNG_EI[vst,var][1],BNG_EI[vst,var][2]
            if (var=='M1' and (vst==1 or vst==2 or vst==3)):
                xmax=wmax-MASS_PIM
            elif (var=='M2' and (vst==1 or vst==2)):
                xmax=wmax-MASS_P
            elif (var=='M2' and (vst==3)):
                xmax=wmax-MASS_PIP
            #print hname,nbins,xmin,xmax
            #! now create hobs(vst,var) for EI
            hobs[EI][vst,var]=ROOT.TH1F(hname,hname,nbins,xmin,xmax) #hobs[AT][vst,var].Clone()
            #hobs[EI][vst,var].Reset()
            hobs[EI][vst,var].SetMarkerColor(CLRS[EI])
            hobs[EI][vst,var].SetLineColor(CLRS[EI])
            hobs[EI][vst,var].SetMarkerStyle(MRKS[EI])
            nrows=len(df.index)
            for ibin in range(nrows):
                binv,val,stat_err=df.loc[ibin,0],df.loc[ibin,1],df.loc[ibin,2]
                #print ibin+1,binv,val,stat_err
                #! Calculate total error
                tot_err=math.sqrt(stat_err**2+(SYST_ERR_EI*val)**2)
                #! Now fill histogram
                hobs[EI][vst,var].SetBinContent(ibin+1,val)
                hobs[EI][vst,var].SetBinError(ibin+1,tot_err)
        
        #! draw together
        #! Histogram aesthetics
        #! Reset any previous settings
        ROOT.gStyle.Reset()
        ROOT.gStyle.SetOptStat(0)
        #!get rid of X error bars
        ROOT.gStyle.SetErrorX(0.001)

        c=ROOT.TCanvas("c","c",1000,1000)
        pad_t=ROOT.TPad("pad_l","Legend pad",0.25,0.935,0.75,1.00)
        pad_p=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
        pad_p.Draw()
        pad_t.Draw()
        pad_t.cd()
        pt=ROOT.TPaveText(.05,.1,.95,.8)
        pt.AddText("Q2_W bin=%s"%q2wbin_dir[AT])
        pt.SetTextSize(0.40)
        pt.Draw()
        pad_p.Divide(3,3)
        for item in PAD_MAP:
            pad,vst,var=item[0],item[1],item[2]
            #print "pad,vst,var=",pad,vst,var
            gpad=pad_p.cd(pad)

            #! fine tuning aesthetics
            #! + ymax
            if q2wbin_dir[AT]=="2.00-2.40_1.725-1.750":
                if var=="THETA": ymax=6.0
            if q2wbin_dir[AT]=="2.40-3.00_1.725-1.750":
                if   var=="M1":    ymax=32
                elif var=="M2":    ymax=25
                elif var=="THETA": ymax=5
                elif var=="ALPHA": ymax=1.5
            if q2wbin_dir[AT]=="2.00-2.40_1.725-1.750":
                if var=="THETA":
                    hobs[AT][vst,var].SetMaximum(ymax)
            if q2wbin_dir[AT]=="2.40-3.00_1.725-1.750":
                hobs[AT][vst,var].SetMaximum(ymax)
            #! Draw
            hobs[AT][vst,var].Draw()
            hobs[EI][vst,var].Draw("same")

            #! Add TLegend if pad==1
            if pad==1:
                l=ROOT.TLegend(0.60,0.80,0.90,0.90)
                l.SetBorderSize(1)
                l.SetTextSize(0.04)
                l.AddEntry(hobs[AT][vst,var],"%s"%DTYP_TITLE[AT],"p")
                l.AddEntry(hobs[EI][vst,var],"%s"%DTYP_TITLE[EI],"p")
                l.Draw()

        #! Save canvas
        c.SaveAs("%s/c%s.png"%(outdir,q2wbin_dir[AT]))
        c.SaveAs("%s/c%s.pdf"%(outdir,q2wbin_dir[AT]))

#! main
compare_itg()
compare_1D()


