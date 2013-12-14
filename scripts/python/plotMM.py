from ROOT import TFile, TCanvas, gROOT
from rootpy.io import root_open
from rootpy.plotting import Hist, HistStack
import rootpy.plotting.root2matplotlib as rplt

import matplotlib.pyplot as plt

DTYPS=2
EXP,SIM = range(DTYPS)
DTYPS_NAME=['exp','sim']

f=[]
# f.append(root_open('/e1f.2pi.anadir2/yield.root'))
# f.append(root_open('/e1f.2pi.anadir2/simdir/yield.root'))
f.append(TFile('/e1f.2pi.anadir2/yield.root'))
f.append(TFile('/e1f.2pi.anadir2/simdir/yield.root'))

def plotMM():
	NMM=4
	hmms=[[],[]]

	for idt in range(DTYPS):
		# hmms[idt].append(f[idt].top.hmmppippimVw.ProjectionY("0",0,-1,"e"))
		# hmms[idt].append(f[idt].top.hmmppimVw.ProjectionY("1",0,-1,"e"))
		# hmms[idt].append(f[idt].top.hmmppipVw.ProjectionY("2",0,-1,"e"))
		# hmms[idt].append(f[idt].top.hmmpippimVw.ProjectionY("3",0,-1,"e"))
		# hmms[idt].append(f[idt].Get("/top/hmmppippimVw"))#.ProjectionY())
		# hmms[idt].append(f[idt].Get("/top/hmmppimVw"))#.ProjectionY())
		# hmms[idt].append(f[idt].Get("/top/hmmppipVw"))#.ProjectionY())
		# hmms[idt].append(f[idt].Get("/top/hmmpippimVw"))#.ProjectionY())
		hmms[idt].append(f[idt].Get("/top/hmmppippimVw").ProjectionY('hmmppippim_%s'%DTYPS_NAME[idt]))
		hmms[idt].append(f[idt].Get("/top/hmmppimVw").ProjectionY('hmmppim_%s'%DTYPS_NAME[idt]))
		hmms[idt].append(f[idt].Get("/top/hmmppipVw").ProjectionY('hmmppip_%s'%DTYPS_NAME[idt]))
		hmms[idt].append(f[idt].Get("/top/hmmpippimVw").ProjectionY('hmmpippim_%s'%DTYPS_NAME[idt]))

	cmm = TCanvas()
	cmm.Divide(2,2)
	for imm in range(NMM):
		cmm.cd(imm+1)
		# hmms[EXP][imm].SetLineColor(gROOT.ProcessLine("kBlue"))
		# hmms[EXP][imm].DrawNormalized("",10000)
		hmms[SIM][imm].SetLineColor(gROOT.ProcessLine("kRed"))
		hmms[SIM][imm].DrawNormalized("",10000)
		hmms[EXP][imm].SetLineColor(gROOT.ProcessLine("kBlue"))
		hmms[EXP][imm].DrawNormalized("sames",10000)
