from ROOT import TFile, TCanvas, TF1, gROOT, gStyle,TMath
from rootpy.io import root_open
from rootpy.plotting import Hist, HistStack
import rootpy.plotting.root2matplotlib as rplt

import matplotlib.pyplot as plt

DTYPS=2
EXP,SIM = range(DTYPS)
DTYPS_NAME=['exp','sim']

f=[]
f.append(root_open('/e1f.2pi.anadir2/yield.root'))
f.append(root_open('/e1f.2pi.anadir2/simdir/yield.root'))
# f.append(TFile('/e1f.2pi.anadir2/yield.root'))
# f.append(TFile('/e1f.2pi.anadir2/simdir/yield.root'))

def plotMM():
	gStyle.SetOptFit(1111)
	NMM=4
	hmms=[[],[]]
	hmm2s=[[],[]]

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
