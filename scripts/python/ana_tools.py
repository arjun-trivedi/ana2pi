from __future__ import division
import ROOT as ROOT
#from ROOT import TFile, TCanvas, TF1, gROOT, gStyle,TMath, gDirectory
from rootpy.io import root_open
from rootpy.plotting import Hist, HistStack,Legend, Canvas
import rootpy.plotting.root2matplotlib as rplt
from rootpy.interactive import wait

import matplotlib.pyplot as plt
import numpy as np
from array import *
import os
from math import *


def gauss(v, par):
    arg = 0;
    if (par[2] != 0): arg = (v[0] - par[1])/par[2];
    fitval = par[0]*(1/(sqrt(2*pi)*par[2]))*exp(-0.5*arg*arg)
    return fitval
def gauss_ppi_hack(v, par):
    arg = 0;
    if (par[2] != 0): arg = (v[0] - par[1])/par[2];
    binw=((0.40-0.00)/100)
    fitval = par[0]*(1/(sqrt(2*pi)*par[2]))*exp(-0.5*arg*arg)*binw;
    return fitval;
def gauss_pippim_hack(v, par):
    arg = 0;
    if (par[2] != 0): arg = (v[0] - par[1])/par[2];
    binw=((1.20-0.80)/100)#((1.40-0.60)/100)
    fitval = par[0]*(1/(sqrt(2*pi)*par[2]))*exp(-0.5*arg*arg)*binw;
    return fitval;

def plot_ana2pi_MMs(be,dtyps=28):#be=beam energy,dtypes for user control
	gpppars_name=[]
	gpppars=[0.0,1.4,4.0]
	for i in range(len(gpppars)):
		for j in range(len(gpppars)):
			for k in range(len(gpppars)):
				gpppars_name.append("%4.2f_%4.2f_%4.2f"%(gpppars[i],gpppars[j],gpppars[k]))
	#print len(gpppars_name)

	DTYPS=dtyps
	ER=0# [1,27] = indices for Simulation runs with differnet parameters
	DTYPS_NAME=[]
	for i in range(DTYPS):
		if i==0:DTYPS_NAME.append("ER")
		else:	DTYPS_NAME.append(gpppars_name[i-1])

	NMM=4
	hmms=[[] for i in range(28)]
	hmm2s=[[] for i in range(28)]
	#print hmms

	f=[]
	f.append(ROOT.TFile("/datadir2/e1f/ana-2pi/exp/q2w2/d2pi.root"))
	for i in range(27):
		f.append(ROOT.TFile("/data/trivedia/e1f/simulation_2pi/ana_new-sim/q2w2_gpptest_%d_011914/recon/d2pi_be%d.root"%(i+1,be)))

	
	OUTDIR_ROOT=os.path.join(os.environ['E1F_SIM2PI_DATADIR'],'ana_new-sim')
	OUTDIR=os.path.join(OUTDIR_ROOT,'be%d'%be)
	if not os.path.isdir(OUTDIR):
		os.mkdir(OUTDIR)

	ROOT.gStyle.SetOptStat("n")
	ROOT.gStyle.SetOptFit(1011)
	# ROOT.gStyle.SetStatY(0.9);
	# ROOT.gStyle.SetStatX(0.9);	
	# ROOT.gStyle.SetStatW(0.4);
	# ROOT.gStyle.SetStatH(0.5);
	CWIDTH=1000
	CHEIGHT=800	
	MARKER_SIZE=2
	ERMAXS=[2000,10000,3000,4000]

	for idt in range(DTYPS):
		topdir=""
		beame=""
		if idt==ER:	
			topdir="top"
			beame=5497
		else:
		    topdir="top2"
		    beame=be
		hmms[idt].append(f[idt].Get("/%s/hmm2ppippimVw"%topdir).ProjectionY('hmm2ppippim_%s_be%d'%(DTYPS_NAME[idt],beame)))
		hmms[idt].append(f[idt].Get("/%s/hmmppipVw"%topdir).ProjectionY('hmmppip_%s_be%d'%(DTYPS_NAME[idt],beame)))
		hmms[idt].append(f[idt].Get("/%s/hmmppimVw"%topdir).ProjectionY('hmmppim_%s_be%d'%(DTYPS_NAME[idt],beame)))
		hmms[idt].append(f[idt].Get("/%s/hmmpippimVw"%topdir).ProjectionY('hmmpippim_%s_be%d'%(DTYPS_NAME[idt],beame)))
		# hmms[idt].append(f[idt].Get("/%s/top1/hmmppippimVw"%topdir).ProjectionY('hmmppippim_%s'%DTYPS_NAME[idt]))
		# hmms[idt].append(f[idt].Get("/%s/top2/hmmppipVw"%topdir).ProjectionY('hmmppip_%s'%DTYPS_NAME[idt]))
		# hmms[idt].append(f[idt].Get("/%s/top3/hmmppimVw"%topdir).ProjectionY('hmmppim_%s'%DTYPS_NAME[idt]))
		# hmms[idt].append(f[idt].Get("/%s/top4/hmmpippimVw"%topdir).ProjectionY('hmmpippim_%s'%DTYPS_NAME[idt]))

		# hmm2s[idt].append(f[idt].Get("/%s/hmm2ppippimVw"%topdir).ProjectionY('hmm2ppippim_%s'%DTYPS_NAME[idt]))
		# hmm2s[idt].append(f[idt].Get("/%s/hmm2ppipVw"%topdir).ProjectionY('hmm2ppip_%s'%DTYPS_NAME[idt]))
		# hmm2s[idt].append(f[idt].Get("/%s/hmm2ppimVw"%topdir).ProjectionY('hmm2ppim_%s'%DTYPS_NAME[idt]))
		# hmm2s[idt].append(f[idt].Get("/%s/hmm2pippimVw"%topdir).ProjectionY('hmm2pippim_%s'%DTYPS_NAME[idt]))

	NFITPARS=2
	MEAN,SGMA=range(NFITPARS)
	mm_fitpars_sim=np.zeros((NFITPARS,3,27),'d')
	mm_fitpars_exp=np.zeros((NFITPARS,3,27),'d')#3rd index is redundant for EXP, but kept for code readability

	yields_SR_mmcut=np.zeros((3,27),'d')
	yields_ER_mmcut=np.zeros((3,27),'d')#3rd index is redundant EXP, but kept for code readability

	lcut_t2t3=ROOT.TLine(0.2,0,0.2,50000)
	lcut_t4=ROOT.TLine(1,0,1,50000)
	for idt in range(1,DTYPS):
		cmm = ROOT.TCanvas("mm_%s"%(gpppars_name[idt-1]),"mm_%s"%(gpppars_name[idt-1]),CWIDTH,CHEIGHT)
		cmm.Divide(2,2)
		for imm in range(NMM):
			pad = cmm.cd(imm+1)
			hmms[ER][imm].SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
			hmms[ER][imm].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlue"))
			hmms[ER][imm].SetMaximum(ERMAXS[imm])
			hmms[ER][imm].Draw()
			pad.Update();
			st=hmms[ER][imm].GetListOfFunctions().FindObject("stats")
			st.SetX1NDC(0.60)
			st.SetX2NDC(1.00)
			st.SetY1NDC(0.625)
			st.SetY2NDC(0.350)
			st.SetFillStyle(4000)
			st.SetTextColor(ROOT.gROOT.ProcessLine("kBlue"));
			
			
			fexp=None
			if imm==1 or imm==2:
				fgauss = ROOT.TF1("fgauss",gauss_ppi_hack,0.0,0.2,3);
				fgauss.SetParameters(1,0,1);
				fgauss.SetParName(0,"Entries")
				fgauss.SetParName(1,"Mean")
				fgauss.SetParName(2,"Sigma")
				hmms[ER][imm].Fit("fgauss","","",0.1,0.17)
				fexp=hmms[ER][imm].GetFunction("fgauss")
				fexp.SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
			elif imm==3:
				fgauss = ROOT.TF1("fgauss",gauss_pippim_hack,0.6,2.0,3);
				fgauss.SetParameters(1,0,1);
				fgauss.SetParName(0,"Entries")
				fgauss.SetParName(1,"Mean")
				fgauss.SetParName(2,"Sigma")
				hmms[ER][imm].Fit("fgauss","","",0.9,0.96)
				fexp=hmms[ER][imm].GetFunction("fgauss")
				fexp.SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))

			nsignal_exp=None
			if fexp is None:nsignal_exp=1000
			else:			nsignal_exp=fexp.GetParameter(0)
			#print "imm,norm=",imm,norm
			
			hmms[idt][imm].SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
			hmms[idt][imm].SetMarkerColor(ROOT.gROOT.ProcessLine("kRed"))
			#estimate signal to background ratio, used to determine norm
			fsim=None
			if imm==1 or imm==2:
				fgauss = ROOT.TF1("fgauss",gauss_ppi_hack,0.0,0.2,3);
				fgauss.SetParameters(1,0,1);
				hmms[idt][imm].Fit("fgauss","0","",0.1,0.17)
				fsim=hmms[idt][imm].GetFunction("fgauss")
			elif imm==3:
				fgauss = ROOT.TF1("fgauss",gauss_pippim_hack,0.6,2.0,3);
				fgauss.SetParameters(1,0,1);
				hmms[idt][imm].Fit("fgauss","0","",0.9,0.96)
				fsim=hmms[idt][imm].GetFunction("fgauss")

			nsignal_sim=None
			if fsim is None:nsignal_sim=nsignal_exp
			else:			nsignal_sim=fsim.GetParameter(0)
				
			norm=nsignal_exp*(hmms[idt][imm].GetEntries()/nsignal_sim)
			hsim=hmms[idt][imm].DrawNormalized("sames",norm)
			pad.Update();
			if imm==1 or imm==2:
				# lcut_t2t3.SetY1(pad.GetUymin())
				# lcut_t2t3.SetY2(pad.GetUymax())
				lcut_t2t3.Draw("same")
			elif imm==3:
				# lcut_t4.SetY1(pad.GetUymin())
				# lcut_t4.SetY2(pad.GetUymax())
				lcut_t4.Draw("same")
			st=hsim.GetListOfFunctions().FindObject("stats")
			st.SetX1NDC(0.60)
			st.SetX2NDC(1.00)
			st.SetY1NDC(0.900)
			st.SetY2NDC(0.625)#(0.625)
			st.SetTextColor(ROOT.gROOT.ProcessLine("kRed"));
			
			fsim=0
			if imm==1 or imm==2:
				hsim.Fit("gaus","","",0.1,0.17)
				fsim = hsim.GetFunction("gaus")
				fsim.SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
			elif imm==3:
				hsim.Fit("gaus","","",0.9,0.96)
				fsim = hsim.GetFunction("gaus")
				fsim.SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
										
			if imm!=0:
				mm_fitpars_sim[MEAN][imm-1][idt-1]=fsim.GetParameter(1)
				mm_fitpars_sim[SGMA][imm-1][idt-1]=fsim.GetParameter(2)
				mm_fitpars_exp[MEAN][imm-1][idt-1]=fexp.GetParameter(1)
				mm_fitpars_exp[SGMA][imm-1][idt-1]=fexp.GetParameter(2)

				fsim_pdf = ROOT.TF1("fsim_pdf",gauss,0.0,0.2,3);
				fsim_pdf.SetParameters(1,fsim.GetParameter(1),fsim.GetParameter(2))
				fexp_pdf = ROOT.TF1("fexp_pdf",gauss,0.0,0.2,3);
				fexp_pdf.SetParameters(1,fexp.GetParameter(1),fexp.GetParameter(2))

				if imm==1 or imm==2:
					yields_SR_mmcut[imm-1][idt-1]=fsim_pdf.Integral(0,0.2)
					yields_ER_mmcut[imm-1][idt-1]=fexp_pdf.Integral(0,0.2)
				elif imm==3:
					yields_SR_mmcut[imm-1][idt-1]=fsim_pdf.Integral(sqrt(0.8),1)
					yields_ER_mmcut[imm-1][idt-1]=fexp_pdf.Integral(sqrt(0.8),1)
		cmm.SaveAs("%s/%s.png"%(OUTDIR,cmm.GetName()))
		cmm.Close()	
		#cmm.Delete()
	
	##--Plot [Mean&Sigma of SR-MM distributions]-[Mean&Sigma of ER-MM distributions] Vs. gpppars
	lgpp_sel=ROOT.TLine(13,-50000,13,50000)
	gpppars_cmbns=array('d',range(27))
	#delta_mm_fitpars=np.subtract(mm_fitpars_exp,mm_fitpars_sim)
	delta_mm_fitpars=np.subtract(mm_fitpars_sim,mm_fitpars_exp)

	
	cname="delta_mm_fitpars" #Note addition of '2', since top1 mm is not analyzed
	c=ROOT.TCanvas(cname,cname,2*CWIDTH,2*CHEIGHT)
	c.Divide(1,2)
	gdmeanVgpp=[]
	gdsgmaVgpp=[]
	leg_mm_fitpars=ROOT.TLegend(0.10,0.90,0.20,0.80);
	for imm in range(3):
		delta_mm_mean=array('d',delta_mm_fitpars[MEAN][imm])
		delta_mm_sgma=array('d',delta_mm_fitpars[SGMA][imm])
		pad=c.cd(1)
		pad.SetGridx()
		gdmeanVgpp.append(ROOT.TGraph(len(gpppars_cmbns),gpppars_cmbns,delta_mm_mean))
		gdmeanVgpp[imm].SetTitle("#mu_{SR}-#mu_{ER} Vs. gpp-pars (be=%d)"%(be))
		gdmeanVgpp[imm].SetMarkerStyle(20+imm)
		gdmeanVgpp[imm].SetMarkerSize(MARKER_SIZE)
		gdmeanVgpp[imm].SetMarkerColor(imm+1)
		x1=gdmeanVgpp[imm].GetHistogram().GetXaxis().GetXmin()#GetBinLowEdge(1)
		x2=gdmeanVgpp[imm].GetHistogram().GetXaxis().GetXmax()#GetBinUpEdge(gfpVgp.GetNbins())
		gdmeanVgpp[imm].GetHistogram().GetXaxis().Set(len(gpppars_cmbns),-0.5,26.5)
		for j in range(len(gpppars_cmbns)):
			gdmeanVgpp[imm].GetHistogram().GetXaxis().SetBinLabel(j+1,gpppars_name[j])
		gdmeanVgpp[imm].GetYaxis().SetTitle("#mu_{SR}-#mu_{ER}(GeV)")
		gdmeanVgpp[imm].SetMinimum(-0.006)
		gdmeanVgpp[imm].SetMaximum(0.015)
		gdmeanVgpp[imm].SetMarkerColor(imm+1)
		leg_mm_fitpars.AddEntry(gdmeanVgpp[imm],"Top%d"%(imm+2),"P")
		if (imm==0):gdmeanVgpp[imm].Draw("AP")
		else:gdmeanVgpp[imm].Draw("P")
		lgpp_sel.Draw("same")
		pad=c.cd(2)
		pad.SetGridx()
		gdsgmaVgpp.append(ROOT.TGraph(len(gpppars_cmbns),gpppars_cmbns,delta_mm_sgma))
		gdsgmaVgpp[imm].SetTitle("#sigma_{SR}-#sigma_{ER} Vs. gpp-pars (be=%d)"%(be))
		gdsgmaVgpp[imm].SetMarkerStyle(20+imm)
		gdsgmaVgpp[imm].SetMarkerSize(MARKER_SIZE)
		gdsgmaVgpp[imm].SetMarkerColor(imm+1)
		x1=gdsgmaVgpp[imm].GetHistogram().GetXaxis().GetXmin()#GetBinLowEdge(1)
		x2=gdsgmaVgpp[imm].GetHistogram().GetXaxis().GetXmax()#GetBinUpEdge(gfpVgp.GetNbins())
		gdsgmaVgpp[imm].GetHistogram().GetXaxis().Set(len(gpppars_cmbns),-0.5,26.5)
		for j in range(len(gpppars_cmbns)):
			gdsgmaVgpp[imm].GetHistogram().GetXaxis().SetBinLabel(j+1,gpppars_name[j])
		gdsgmaVgpp[imm].GetYaxis().SetTitle("#sigma_{SR}-#sigma_{ER}(GeV)")
		gdsgmaVgpp[imm].SetMinimum(-0.02)
		gdsgmaVgpp[imm].SetMaximum(0.03)
		if (imm==0):gdsgmaVgpp[imm].Draw("AP")
		else:gdsgmaVgpp[imm].Draw("P")
		lgpp_sel.Draw("same")
	c.cd(1)
	leg_mm_fitpars.Draw("same")	
	c.SaveAs("%s/%s.png"%(OUTDIR,c.GetName()))
	c.Close()	

	##-- Plot dyield(%) vs gpppars	
	##   dyield = [EC(from EA)-EC(from SA)]/EC(from EA) = 1-[ER(mmcut)/SR(mmcut)]
	dyields=np.divide(yields_ER_mmcut,yields_SR_mmcut)
	# dyields=np.multiply(dyields,-1)
	# dyields=np.add(dyields,1)
	dyields=np.add(dyields,-1)
	dyields=np.multiply(dyields,100)
	leg_dyields=ROOT.TLegend(0.10,0.90,0.20,0.80);
	cname="dyields_top" #Note addition of '2', since top1 mm is not analyzed
	c=ROOT.TCanvas(cname,cname,4*CWIDTH,2*CHEIGHT)
	c.SetGridx()
	g=[]
	for imm in range(3):
		dy=array('d',dyields[imm])
		g.append(ROOT.TGraph(len(gpppars_cmbns),gpppars_cmbns,dy))
		g[imm].SetTitle("#DeltaEC/EC Vs. gpp-pars(be=%d)"%(be))
		g[imm].SetMarkerStyle(20+imm)
		g[imm].SetMarkerSize(MARKER_SIZE+2)
		g[imm].SetMarkerColor(imm+1)
		leg_dyields.AddEntry(g[imm],"Top%d"%(imm+2),"P")
		x1=g[imm].GetHistogram().GetXaxis().GetXmin()#GetBinLowEdge(1)
		x2=g[imm].GetHistogram().GetXaxis().GetXmax()#GetBinUpEdge(gfpVgp.GetNbins())
		g[imm].GetHistogram().GetXaxis().Set(len(gpppars_cmbns),-0.5,26.5)#x1,x2);
		for j in range(len(gpppars_cmbns)):
			g[imm].GetHistogram().GetXaxis().SetBinLabel(j+1,gpppars_name[j])
		g[imm].GetYaxis().SetTitle("#DeltaEC/EC(%)")
		g[imm].SetMinimum(-5.)
		g[imm].SetMaximum(20.)
		if imm==0:g[imm].Draw("AP")
		else:g[imm].Draw("P")
		lgpp_sel.Draw("same")
	leg_dyields.Draw("same")
	c.SaveAs("%s/%s.png"%(OUTDIR,c.GetName()))
	c.Close()
	
	if not ROOT.gROOT.IsBatch():
		plt.show()
		# wait for you to close the ROOT canvas before exiting
		wait(True)

def plot_elastic_W(nentries=1000000000):
	ROOT.gStyle.SetOptFit(1111)
	
	DTYPS=5
	ER_NMCR,ER_YMCR,SR_NGPP,SR_YGPP,ST = range(DTYPS)
	DTYPS_NAME=['ER-nomcorr','ER-yesmcorr','SR-nogpp','SR-yesgpp','ST']

	f=[]
	f.append(ROOT.TFile("/datadir2/e1f/ana-elast/exp/delast.root"))
	f.append(ROOT.TFile("/datadir2/e1f/ana-elast/exp/delast.root"))
	f.append(ROOT.TFile("/datadir2/e1f/ana-elast/sim/elast_gpp-no_011014/delast.root"))
	f.append(ROOT.TFile("/datadir2/e1f/ana-elast/sim/elast_gpp-yes_011014/delast.root"))
	f.append(ROOT.TFile("/datadir2/e1f/ana-elast/sim/elast_gpp-yes_011014/delast_ST.root"))

	hW=[]
	hp=[[],[],[]]
	hv=[[],[],[]]
		
	t=""
	for i in range(DTYPS):
		print "iteration #:",i
		if i==ST: break
		f[i].cd()
		if i==ER_NMCR: 
			print "i=ER_NMCR"
			t=f[i].Get("/delast/t")
		elif i==ER_YMCR: 
			print "i=ER_YMCR"
			t=f[i].Get("/delast2/t")
		elif i==ST:
			print "i=ST"
			t=f[i].Get("/delast/t_ST");
		else:
			print "i=SR"
			t=f[i].Get("/delast/t");

		# print "Tree title:",t.GetTitle()
		# print "Tree entries:",t.GetEntries()
		if i==ST:
			t.Draw("W>>hW(600,0.6,1.2)","","",nentries)
		else:
			t.Draw("W>>hW(100,0.6,1.2)","","",nentries)
		cut=ROOT.TCut("W<1.1");# && lvE.Theta()*TMath::RadToDeg()>15")
		t.Draw("lvE.X()>>hpx(100,-2,2)",cut,"",nentries)
		t.Draw("lvE.Y()>>hpy(100,-2,2)",cut,"",nentries)
		t.Draw("lvE.Z()>>hpz(100,0,6)",cut,"",nentries)
		t.Draw("vxE.X()>>hvx(100,-1,1)",cut,"",nentries)
		t.Draw("vxE.Y()>>hvy(100,-1,1)",cut,"",nentries)
		t.Draw("vxE.Z()>>hvz(100,-50,0)",cut,"",nentries)
		hW.append(ROOT.gDirectory.Get("hW"))
		hp[0].append(ROOT.gDirectory.Get("hpx"))
		hp[1].append(ROOT.gDirectory.Get("hpy"))
		hp[2].append(ROOT.gDirectory.Get("hpz"))
		hv[0].append(ROOT.gDirectory.Get("hvx"))
		hv[1].append(ROOT.gDirectory.Get("hvy"))
		hv[2].append(ROOT.gDirectory.Get("hvz"))
		# print "pwd:"
		# ROOT.gDirectory.pwd()
		# print "---- gDirectory Content ---- :"
		# ROOT.gDirectory.ls()
		
	c_hW = Canvas(name="W",title="W")
	c_hW.Divide(3,2)
	c_DC=Canvas(name="DC",title="DC")
	c_DC.Divide(3,2);
	norm=1000;
	for i in range(DTYPS):
		print "iteration #:",i
		if i==ST: break
		ROOT.gStyle.SetOptStat("nemMrRiuo")
		c_hW.cd(i+1)
		hW[i].SetName(DTYPS_NAME[i]+"_"+hW[i].GetName())
		hW[i].Draw()
		if i!=ST:
			hW[i].Fit("gaus","","",0.90,0.98)

		ROOT.gStyle.SetOptStat("nemMrR")
		for j in range(3):
			#pad=c_hp[j].cd()
			pad=c_DC.cd(j+1)
			hp[j][i].SetLineColor(i+1)
			hp[j][i].SetName(DTYPS_NAME[i]+"_"+hp[j][i].GetName())
			if i==0:
				hn=hp[j][i].DrawNormalized("HIST",norm)
			else:
				hn=hp[j][i].DrawNormalized("HIST sames",norm)
			pad.Update();
			#st=hp[j][i].GetListOfFunctions().FindObject("stats");
			st=hn.GetListOfFunctions().FindObject("stats");
			st.SetTextColor(i+1);
			st.Draw()

		for j in range(3):
			#pad=c_hv[j].cd()
			pad=c_DC.cd(j+1+3)
			hv[j][i].SetLineColor(i+1)
			hv[j][i].SetOption("")
			hv[j][i].SetName(DTYPS_NAME[i]+"_"+hv[j][i].GetName())
			if i==0:
				hn=hv[j][i].DrawNormalized("HIST",norm)
			else:
				hn=hv[j][i].DrawNormalized("HIST sames",norm)
			pad.Update();
			#st=hv[j][i].GetListOfFunctions().FindObject("stats");
			st=hn.GetListOfFunctions().FindObject("stats");
			st.SetTextColor(i+1);
			st.Draw()

	if not ROOT.gROOT.IsBatch():
		plt.show()
		# wait for you to close the ROOT canvas before exiting
		wait(True)
