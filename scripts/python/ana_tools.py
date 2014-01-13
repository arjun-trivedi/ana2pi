import ROOT as ROOT
#from ROOT import TFile, TCanvas, TF1, gROOT, gStyle,TMath, gDirectory
from rootpy.io import root_open
from rootpy.plotting import Hist, HistStack,Legend, Canvas
import rootpy.plotting.root2matplotlib as rplt
from rootpy.interactive import wait

import matplotlib.pyplot as plt

DTYPS=2
EXP,SIM = range(DTYPS)
DTYPS_NAME=['exp','sim']

def plotMM():
	gStyle.SetOptFit(1111)
	NMM=4
	hmms=[[],[]]
	hmm2s=[[],[]]

	f=[]
	f.append(root_open('/e1f.2pi.anadir2/yield.root'))
	f.append(root_open('/e1f.2pi.anadir2/simdir/yield.root'))
	#f = TF1("gaus","gaus",)

	for idt in range(DTYPS):
		hmms[idt].append(f[idt].top.hmmppippimVw.ProjectionY('hmmppippim_%s'%DTYPS_NAME[idt]))
		hmms[idt].append(f[idt].top.hmmppipVw.ProjectionY('hmmppip_%s'%DTYPS_NAME[idt]))
		hmms[idt].append(f[idt].top.hmmppimVw.ProjectionY('hmmppim_%s'%DTYPS_NAME[idt]))
		hmms[idt].append(f[idt].top.hmmpippimVw.ProjectionY('hmmpippim_%s'%DTYPS_NAME[idt]))
		hmm2s[idt].append(f[idt].top.hmm2ppippimVw.ProjectionY('hmm2ppippim_%s'%DTYPS_NAME[idt]))
		hmm2s[idt].append(f[idt].top.hmm2ppipVw.ProjectionY('hmm2ppip_%s'%DTYPS_NAME[idt]))
		hmm2s[idt].append(f[idt].top.hmm2ppimVw.ProjectionY('hmm2ppim_%s'%DTYPS_NAME[idt]))
		hmm2s[idt].append(f[idt].top.hmm2pippimVw.ProjectionY('hmm2pippim_%s'%DTYPS_NAME[idt]))
		
		# hmms[idt].append(f[idt].Get("/top/hmmppippimVw").ProjectionY('hmmppippim_%s'%DTYPS_NAME[idt]))
		# hmms[idt].append(f[idt].Get("/top/hmmppimVw").ProjectionY('hmmppim_%s'%DTYPS_NAME[idt]))
		# hmms[idt].append(f[idt].Get("/top/hmmppipVw").ProjectionY('hmmppip_%s'%DTYPS_NAME[idt]))
		# hmms[idt].append(f[idt].Get("/top/hmmpippimVw").ProjectionY('hmmpippim_%s'%DTYPS_NAME[idt]))

	cmm = TCanvas("mm","mm")
	cmm.Divide(2,2)
	fsimc=[]
	for imm in range(NMM):
		pad = cmm.cd(imm+1)
		# hmms[EXP][imm].SetLineColor(gROOT.ProcessLine("kBlue"))
		# hmms[EXP][imm].DrawNormalized("",10000)
		hmms[SIM][imm].SetLineColor(gROOT.ProcessLine("kRed"))
		hmms[SIM][imm].SetMarkerColor(gROOT.ProcessLine("kRed"))
		hsim = hmms[SIM][imm].DrawNormalized("",10000)
		hmms[EXP][imm].SetLineColor(gROOT.ProcessLine("kBlue"))
		hmms[EXP][imm].SetMarkerColor(gROOT.ProcessLine("kBlue"))
		hexp=hmms[EXP][imm].DrawNormalized("sames",10000)

		if imm==1 or imm==2:
			hsim.Fit("gaus","","",0.1,0.17)
			hexp.Fit("gaus","","",0.1,0.17)
		elif imm==3:
			hsim.Fit("gaus","","",0.9,1.0)
			hexp.Fit("gaus","","",0.9,0.96)
		if imm!=0:
			fsim = hsim.GetFunction("gaus")
			fsim.SetLineColor(gROOT.ProcessLine("kBlack"))
			fexp = hexp.GetFunction("gaus")
			fexp.SetLineColor(gROOT.ProcessLine("kBlack"))
			pad.Update()
			fsimc.append(TF1(fsim))
			fsimc[imm-1].SetLineColor(gROOT.ProcessLine("kGreen"))
			fsimc[imm-1].SetRange(hmms[SIM][imm].GetXaxis().GetXmin(),hmms[SIM][imm].GetXaxis().GetXmax())
			fsimc[imm-1].Draw("same")
			#pad.Update()

	cmm2 = TCanvas("mm2","mm2")
	cmm2.Divide(2,2)
	for imm in range(NMM):
		pad = cmm2.cd(imm+1)
		# hmms[EXP][imm].SetLineColor(gROOT.ProcessLine("kBlue"))
		# hmms[EXP][imm].DrawNormalized("",10000)
		hmm2s[SIM][imm].SetLineColor(gROOT.ProcessLine("kRed"))
		hmm2s[SIM][imm].SetMarkerColor(gROOT.ProcessLine("kRed"))
		hsim = hmm2s[SIM][imm].DrawNormalized("",10000)
		hmm2s[EXP][imm].SetLineColor(gROOT.ProcessLine("kBlue"))
		hmm2s[EXP][imm].SetMarkerColor(gROOT.ProcessLine("kBlue"))
		hexp=hmm2s[EXP][imm].DrawNormalized("sames",10000)

		if imm==1 or imm==2:
			hsim.Fit("gaus","","",0.1,0.17)
			hexp.Fit("gaus","","",0.1,0.17)
		elif imm==3:
			hsim.Fit("gaus","","",0.9,1.0)
			hexp.Fit("gaus","","",0.9,0.96)
		if imm!=0:
			hsim.GetFunction("gaus").SetLineColor(gROOT.ProcessLine("kBlack"))
			hexp.GetFunction("gaus").SetLineColor(gROOT.ProcessLine("kBlack"))
			pad.Update()

	if not gROOT.IsBatch():
		plt.show()
		# wait for you to close the ROOT canvas before exiting
		wait(True)

def elastic_study(nentries=1000000000):
	ROOT.gStyle.SetOptFit(1111)
	
	DTYPS=4
	ER,SR_NGPP,SR_YGPP,ST = range(DTYPS)
	DTYPS_NAME=['exp-recon','sim-recon-nogpp','sim-recon-yesgpp','sim-thrown']

	f=[]
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
		f[i].cd()
		if i==ER: 
			print "i=ER"
			t=f[i].Get("/delast2/t")#;delast2.t#Get("/delast2/t");
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
		cut=ROOT.TCut("W<1.1 && lvE.Theta()*TMath::RadToDeg()>15")
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
	c_hW.Divide(2,2)
	c_DC=Canvas(name="DC",title="DC")
	c_DC.Divide(3,2);
	norm=1000;
	for i in range(DTYPS):
		print "iteration #:",i
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
