#!/usr/bin/python

from ROOT import *
from array import array

#SET ROOT ENVIRONMENT
gSystem.Load('libEpOmega')
gSystem.AddIncludePath("-I/home/ephelps/w/epOmega")
gStyle.SetPalette(1);
gStyle.SetOptStat(0000)
gStyle.SetOptFit(0000)
gROOT.ProcessLine(".x styBABAR.C");
gROOT.ProcessLine(".x styPUB.C");
gStyle.SetCanvasDefH(300);
gStyle.SetCanvasDefW(300);
goodcolors = [ kRed+1, kGreen+1, kBlue, kYellow+1, kMagenta+1, kCyan+1, 9 ]
goodcolors.extend(reversed(range(7,49)))

fin = TFile("/home/ephelps/w/epOmega/sim.root")
fout = TFile("acc_lm.root","recreate")
fin.ls()
bdt_tmp = fin.Get('tops/mc/top1/yield_morand')
bdr_tmp = fin.Get('tops2/top1/yield_morand')
kindims = array('i',[0,1,2,3])
bdt = bdt_tmp.Projection(4,kindims)
bdr = bdr_tmp.Projection(4,kindims)
bcoords = array('i',[1,1,1,1])
acc = bdr.Clone("acc_top1")
acc.Sumw2()
bdt.Sumw2()
totalbins = acc.GetNbins()
for i in range(0,totalbins):
    bdt.SetBinError(i,sqrt(bdt.GetBinContent(i)))
    acc.SetBinError(i,sqrt(acc.GetBinContent(i)))
acc.Divide(bdt)

nbinsQ2 = acc.GetAxis(0).GetNbins()
nbinsXb = acc.GetAxis(1).GetNbins()
nbinsT = acc.GetAxis(2).GetNbins()
nbinsPhi = acc.GetAxis(3).GetNbins()
q2edges = []
xbedges = []
for q2bin in range(1,nbinsQ2+2):
    q2edges.append(acc.GetAxis(0).GetBinCenter(q2bin))
for xbbin in range(1,nbinsXb+2):
    xbedges.append(acc.GetAxis(1).GetBinCenter(xbbin))
xbmin = xbedges[0]
xbmax = xbedges[nbinsXb]
q2min = q2edges[0]
q2max = q2edges[nbinsQ2]
xbrang = xbmax-xbmin
q2rang = q2max-q2min

cqx = TCanvas("cqx_tVphi","t vs phi",900,600)
cx1 = xbmin-0.1*xbrang
cx2 = xbmax+0.05*xbrang;
cy1 = q2min-0.1*q2rang;
cy2 = q2max+0.1*q2rang;
cqx.Range(cx1,cy1,cx2,cy2);
px1 = (xbmin-cx1)/(cx2-cx1)+0.01
px2 = (xbmax-cx1)/(cx2-cx1)
py1 = (q2min-cy1)/(cy2-cy1)+0.01
py2 = (q2max-cy1)/(cy2-cy1)
qxpad = TPad("qxpad","qxpad",px1,py1,px2,py2)
qxpad.Draw()
qxpad.Divide(nbinsXb,nbinsQ2)
f1 = TF1("f1","-x",q2min,q2max)
xaxisg = TGaxis(xbmin,q2min,xbmax,q2min,xbmin,xbmax,nbinsXb,"")
yaxisg = TGaxis(xbmin,q2min,xbmin,q2max,"f1",nbinsQ2,"-") #q2min,q2max,nbinsQ2,"")
yaxisg.SetTitle("Q^{2} [Gev^{2}]")
xaxisg.SetTitle("x_{b}")
xaxisg.Draw()
yaxisg.Draw()


histnum = 0
hists = []
acchists = []
errhists = []
accmaster = TH1F("hacchist","acceptance distribution",35,0,0.35)
errmaster = TH1F("herrhist","relative error distribution",50,0,0.5)
for q2bin in range(1,nbinsQ2+1):
    q2val = acc.GetAxis(0).GetBinCenter(q2bin)
    acc.GetAxis(0).SetRange(q2bin,q2bin)
    fout.cd()
    q2dir = gDirectory.mkdir("q2_%.2f"%q2val)
    for xbbin in range(1,nbinsXb+1):
        histnum+=1
        xbval = acc.GetAxis(1).GetBinCenter(xbbin)
        acc.GetAxis(1).SetRange(xbbin,xbbin)
        fout.cd()
        htvphi = acc.Projection(2,3)
        htvphi.SetTitle("Q^{2} = %.3f, x_{b} = %.3f" % (q2val,xbval))
        hists.append(htvphi)
        qxpad.cd(histnum)
        deletethis = hists[histnum-1].GetEntries()<=0
        if deletethis:
            htvphi.SetName('trash%i'%histnum)
            gDirectory.Delete('%s;*'%htvphi.GetName())
        else:
            dir = q2dir.mkdir("xb_%.2f"%xbval)
            dir.cd()
            hists[histnum-1].SetDirectory(dir)
            hists[histnum-1].Draw("colz")
            nbins = nbinsT*nbinsPhi
            acchist = TH1F("hacchist","Q^{2} = %.3f, x_{b} = %.3f" % (q2val,xbval),20,0,0.2)
            errhist = TH1F("herrhist","Q^{2} = %.3f, x_{b} = %.3f" % (q2val,xbval),200,0,2)
            acchists.append(acchist)
            errhists.append(errhist)
            for b in range(1,nbins+1):
                bc = hists[histnum-1].GetBinContent(b)
                be = hists[histnum-1].GetBinError(b)
                if bc > 0:
                    accmaster.Fill(bc)
                    errmaster.Fill(be/bc)
                    acchist.Fill(bc)
                    errhist.Fill(be/bc)
fout.WriteObject(xaxisg,"xaxisg")
fout.WriteObject(yaxisg,"yaxisg")
fout.Append(f1)
fout.WriteObject(cqx,"cqx")
fout.Append(acc)
fout.Write()
