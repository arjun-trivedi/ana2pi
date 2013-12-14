from ROOT import TFile, TCanvas
DTYPE=2
EXP,SIM = range(DTYPE)
print EXP,SIM
f=[]
f.append(TFile('/e1f.2pi.anadir2/yield.root'))
f.append(TFile('/e1f.2pi.anadir2/simdir/yield.root'))

NMM=4
hmm=[[],[]]
hmm=[[],[]]

print f[EXP].GetName()

for idt in range(DTYPE):
	#for imm in range(NMM)
	# print f[idt].ls("top/hmm")
	# h=f[idt].Get("/top/hmmppippimVw")
	# print h.GetName()
	hmm[idt].append(f[idt].Get("/top/hmmppippimVw"))
	hmm[idt].append(f[idt].Get("/top/hmmppimVw"))
	hmm[idt].append(f[idt].Get("/top/hmmppipVw"))
	hmm[idt].append(f[idt].Get("/top/hmmpippimVw"))
# hmm[EXP][0].GetName();
cmm = TCanvas()
cmm.Divide(2,2)
for imm in range(NMM):
	cmm.cd(imm+1)
	hmm[EXP][imm].Draw()


#fexp = TFile('/e1f.2pi.anadir2/yiel

