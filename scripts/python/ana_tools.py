import ROOT as ROOT
#from ROOT import TFile, TCanvas, TF1, gROOT, gStyle,TMath, gDirectory
from rootpy.io import root_open
from rootpy.plotting import Hist, HistStack,Legend, Canvas
import rootpy.plotting.root2matplotlib as rplt
from rootpy.interactive import wait

import matplotlib.pyplot as plt

def plot_ana2pi_MMs():
	DTYPS=2
	EXP,SIM = range(DTYPS)
	DTYPS_NAME=['exp','sim']

	ROOT.gStyle.SetOptFit(1111)
	NMM=4
	hmms=[[],[]]
	hmm2s=[[],[]]

	f=[]
	f.append(ROOT.TFile("/datadir2/e1f/ana-2pi/exp/q2w1/d2pi.root"))
	#f.append(ROOT.TFile("/e1f.2pi.anadir1/simdir/yield.root"))
	f.append(ROOT.TFile("/data/trivedia/e1f/simulation_2pi/test_sim/q2wf_gpp-ep_011714/recon/d2pi.root"))
	
	# f.append(ROOT.TFile("/datadir2/e1f/ana-2pi/exp/q2w2/d2pi.root"))
	# f.append(ROOT.TFile("/e1f.2pi.anadir2/simdir/yield.root"))
	

	for idt in range(DTYPS):
		topdir=""
		if idt==EXP:
			topdir="top"
		elif idt==SIM:
			topdir="top2"
		hmms[idt].append(f[idt].Get("/%s/hmmppippimVw"%topdir).ProjectionY('hmmppippim_%s'%DTYPS_NAME[idt]))
		hmms[idt].append(f[idt].Get("/%s/hmmppimVw"%topdir).ProjectionY('hmmppim_%s'%DTYPS_NAME[idt]))
		hmms[idt].append(f[idt].Get("/%s/hmmppipVw"%topdir).ProjectionY('hmmppip_%s'%DTYPS_NAME[idt]))
		hmms[idt].append(f[idt].Get("/%s/hmmpippimVw"%topdir).ProjectionY('hmmpippim_%s'%DTYPS_NAME[idt]))
		hmm2s[idt].append(f[idt].Get("/%s/hmm2ppippimVw"%topdir).ProjectionY('hmm2ppippim_%s'%DTYPS_NAME[idt]))
		hmm2s[idt].append(f[idt].Get("/%s/hmm2ppimVw"%topdir).ProjectionY('hmm2ppim_%s'%DTYPS_NAME[idt]))
		hmm2s[idt].append(f[idt].Get("/%s/hmm2ppipVw"%topdir).ProjectionY('hmm2ppip_%s'%DTYPS_NAME[idt]))
		hmm2s[idt].append(f[idt].Get("/%s/hmm2pippimVw"%topdir).ProjectionY('hmm2pippim_%s'%DTYPS_NAME[idt]))

	cmm = ROOT.TCanvas("mm","mm")
	cmm.Divide(2,2)
	# fsimc=[]
	for imm in range(NMM):
		pad = cmm.cd(imm+1)
		hmms[SIM][imm].SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
		hmms[SIM][imm].SetMarkerColor(ROOT.gROOT.ProcessLine("kRed"))
		hsim = hmms[SIM][imm].DrawNormalized("",10000)
		hmms[EXP][imm].SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
		hmms[EXP][imm].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlue"))
		hexp=hmms[EXP][imm].DrawNormalized("sames",10000)

		if imm==1 or imm==2:
			hsim.Fit("gaus","","",0.1,0.17)
			hexp.Fit("gaus","","",0.1,0.17)
		elif imm==3:
			hsim.Fit("gaus","","",0.9,1.0)
			hexp.Fit("gaus","","",0.9,0.96)
		if imm!=0:
			fsim = hsim.GetFunction("gaus")
			fsim.SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
			fexp = hexp.GetFunction("gaus")
			fexp.SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
			pad.Update()
			# fsimc.append(ROOT.TF1(fsim))
			# fsimc[imm-1].SetLineColor(ROOT.gROOT.ProcessLine("kGreen"))
			# fsimc[imm-1].SetRange(hmms[SIM][imm].GetXaxis().GetXmin(),hmms[SIM][imm].GetXaxis().GetXmax())
			# fsimc[imm-1].Draw("same")
			#pad.Update()

	cmm2 = ROOT.TCanvas("mm2","mm2")
	cmm2.Divide(2,2)
	for imm in range(NMM):
		pad = cmm2.cd(imm+1)
		hmm2s[SIM][imm].SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
		hmm2s[SIM][imm].SetMarkerColor(ROOT.gROOT.ProcessLine("kRed"))
		hsim = hmm2s[SIM][imm].DrawNormalized("",10000)
		hmm2s[EXP][imm].SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
		hmm2s[EXP][imm].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlue"))
		hexp=hmm2s[EXP][imm].DrawNormalized("sames",10000)

		if imm==1 or imm==2:
			hsim.Fit("gaus","","",0.1,0.17)
			hexp.Fit("gaus","","",0.1,0.17)
		elif imm==3:
			hsim.Fit("gaus","","",0.9,1.0)
			hexp.Fit("gaus","","",0.9,0.96)
		if imm!=0:
			hsim.GetFunction("gaus").SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
			hexp.GetFunction("gaus").SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
			pad.Update()

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
