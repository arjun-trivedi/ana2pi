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

def plot_ana2pi_MMs(be):#be=beam energy
	gpppars_name=[]
	gpppars=[0.0,1.4,4.0]
	for i in range(len(gpppars)):
		for j in range(len(gpppars)):
			for k in range(len(gpppars)):
				gpppars_name.append("%4.2f_%4.2f_%4.2f"%(gpppars[i],gpppars[j],gpppars[k]))
	#print len(gpppars_name)

	DTYPS=28
	ER=0
	DTYPS_NAME=[]
	for i in range(DTYPS):
		if i==0:
			DTYPS_NAME.append("ER")
		else:
			#DTYPS_NAME.append("SR%d"%(i))
			DTYPS_NAME.append(gpppars_name[i-1])

	NMM=4
	# hmms=[[],[]]
	# hmm2s=[[],[]]
	hmms=[[] for i in range(28)]
	hmm2s=[[] for i in range(28)]
	#print hmms

	ROOT.gStyle.SetOptFit(1111)
	CWIDTH=1000
	CHEIGHT=800	

	OUTDIR_ROOT=os.path.join(os.environ['E1F_SIM2PI_DATADIR'],'ana_new-sim')
	OUTDIR=os.path.join(OUTDIR_ROOT,'be%d'%be)
	if not os.path.isdir(OUTDIR):
		os.mkdir(OUTDIR)

	f=[]
	f.append(ROOT.TFile("/datadir2/e1f/ana-2pi/exp/q2w2/d2pi.root"))
	for i in range(27):
		f.append(ROOT.TFile("/data/trivedia/e1f/simulation_2pi/ana_new-sim/q2w2_gpptest_%d_011914/recon/d2pi_be%d.root"%(i+1,be)))

	# f.append(ROOT.TFile("/datadir2/e1f/ana-2pi/exp/q2w2/d2pi.root"))
	# f.append(ROOT.TFile("/e1f.2pi.anadir2/simdir/yield.root"))
	

	for idt in range(DTYPS):
		topdir=""
		if idt==ER:
			topdir="top"
		# elif idt!=ER:
		# 	topdir="top2"
		else:
			topdir="top2"
		# print idt,topdir,f[idt].GetName()
		hmms[idt].append(f[idt].Get("/%s/hmmppippimVw"%topdir).ProjectionY('hmmppippim_%s'%DTYPS_NAME[idt]))
		hmms[idt].append(f[idt].Get("/%s/hmmppimVw"%topdir).ProjectionY('hmmppim_%s'%DTYPS_NAME[idt]))
		hmms[idt].append(f[idt].Get("/%s/hmmppipVw"%topdir).ProjectionY('hmmppip_%s'%DTYPS_NAME[idt]))
		hmms[idt].append(f[idt].Get("/%s/hmmpippimVw"%topdir).ProjectionY('hmmpippim_%s'%DTYPS_NAME[idt]))
		hmm2s[idt].append(f[idt].Get("/%s/hmm2ppippimVw"%topdir).ProjectionY('hmm2ppippim_%s'%DTYPS_NAME[idt]))
		hmm2s[idt].append(f[idt].Get("/%s/hmm2ppimVw"%topdir).ProjectionY('hmm2ppim_%s'%DTYPS_NAME[idt]))
		hmm2s[idt].append(f[idt].Get("/%s/hmm2ppipVw"%topdir).ProjectionY('hmm2ppip_%s'%DTYPS_NAME[idt]))
		hmm2s[idt].append(f[idt].Get("/%s/hmm2pippimVw"%topdir).ProjectionY('hmm2pippim_%s'%DTYPS_NAME[idt]))

	# for idt in range(DTYPS):
	# 	for imm in range(NMM):
	# 		print "hmm[%d][%d]"%(idt+1,imm+1),hmms[idt][imm].GetName()
	
	#mm_fitpars=[[[] for j in range(3)] for i in range(27)]
	#mm_fitpars=[[[[]for i in range(27)] for i in range(3)],[[[]for i in range(27)] for i in range(3)]]
	mm_fitpars=np.zeros((2,3,27),'d')
	mm_fitpars_exp=np.zeros((2,3),'d')
	for idt in range(1,28):
		cmm = ROOT.TCanvas("mm_%s"%(gpppars_name[idt-1]),"mm_%s"%(gpppars_name[idt-1]),CWIDTH,CHEIGHT)
		cmm.Divide(2,2)
		for imm in range(NMM):
			pad = cmm.cd(imm+1)
			hmms[idt][imm].SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
			hmms[idt][imm].SetMarkerColor(ROOT.gROOT.ProcessLine("kRed"))
			hsim = hmms[idt][imm].DrawNormalized("",10000)
			pad.Update();
			st=hsim.GetListOfFunctions().FindObject("stats")
			st.SetX1NDC(0.50)
			st.SetX2NDC(1.00)
			st.SetY1NDC(1.00)
			st.SetY2NDC(0.625)
			print "hsim[%d][%d]"%(idt+1,imm+1),hsim.GetName()
			hmms[ER][imm].SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
			hmms[ER][imm].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlue"))
			hexp=hmms[ER][imm].DrawNormalized("sames",10000)
			pad.Update();
			st=hexp.GetListOfFunctions().FindObject("stats")
			st.SetX1NDC(0.50)
			st.SetX2NDC(1.00)
			st.SetY1NDC(0.625)
			st.SetY2NDC(0.250)
			print "hexp=",hexp.GetName()
			
			
			if imm==1 or imm==2:
				hsim.Fit("gaus","","",0.1,0.17)
				hexp.Fit("gaus","","",0.1,0.17)
			elif imm==3:
				hsim.Fit("gaus","","",0.9,0.96)
				hexp.Fit("gaus","","",0.9,0.96)
			
			if imm!=0:
				fsim = hsim.GetFunction("gaus")
				fsim.SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
				fexp = hexp.GetFunction("gaus")
				fexp.SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
				pad.Update()
				# mm_fitpars[idt-1][imm-1].append(fsim.GetParameter(1))
				# mm_fitpars[idt-1][imm-1].append(fsim.GetParameter(2))
				mm_fitpars[0][imm-1][idt-1]=fsim.GetParameter(1)
				mm_fitpars[1][imm-1][idt-1]=fsim.GetParameter(2)
				mm_fitpars_exp[0][imm-1]=fexp.GetParameter(1)
				mm_fitpars_exp[1][imm-1]=fexp.GetParameter(2)
		cmm.SaveAs("%s/%s.png"%(OUTDIR,cmm.GetName()))
		cmm.Close()	

	cmm_fitpars = ROOT.TCanvas("fit_pars","fit_pars",2*CWIDTH,2*CHEIGHT)
	cmm_fitpars.Divide(1,2)
	# gpppars_cmbns=range(27)
	gpppars_cmbns=np.zeros((27),'d')
	gpppars_cmbns=np.arange(27)
	# print len(gpppars_cmbns),len(mm_fitpars[0][0])
	# print gpppars_cmbns
	# print mm_fitpars[0][0]
	tx=array('d',gpppars_cmbns)
	tmean=array('d',mm_fitpars[0][0])
	tmean_diff=np.subtract(tmean,mm_fitpars_exp[0][0])
	tsigma=array('d',mm_fitpars[1][0])
	tsigma_diff=np.subtract(tsigma,mm_fitpars_exp[1][0])
	#gfpVgp = ROOT.TGraph(len(gpppars_cmbns),gpppars_cmbns,mm_fitpars[0][0])
	print len(tx),len(tmean)
	print tx
	print tmean
	gfpVgp=[]
	gfpVgp.append(ROOT.TGraph(len(tx),tx,tmean_diff))
	gfpVgp.append(ROOT.TGraph(len(tx),tx,tsigma_diff))
	for i in range(2):
		x1=gfpVgp[i].GetHistogram().GetXaxis().GetXmin()#GetBinLowEdge(1)
		x2=gfpVgp[i].GetHistogram().GetXaxis().GetXmax()#GetBinUpEdge(gfpVgp.GetNbins())
		gfpVgp[i].GetHistogram().GetXaxis().Set(len(tx),-0.5,26.5)#x1,x2);
		for j in range(len(tx)):
			gfpVgp[i].GetHistogram().GetXaxis().SetBinLabel(j+1,gpppars_name[j])
		gfpVgp[i].GetXaxis().SetTitle("gpp-pars")
		if i ==0:
			gfpVgp[i].GetYaxis().SetTitle("delta-mean")
		elif i==1:
			gfpVgp[i].GetYaxis().SetTitle("delta-sigma")
		cmm_fitpars.cd(i+1)	
		gfpVgp[i].Draw("ALP")	
	cmm_fitpars.SaveAs("%s/%s.png"%(OUTDIR,cmm_fitpars.GetName()))
	cmm_fitpars.Close()

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
