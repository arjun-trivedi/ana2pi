#!/usr/bin/python
import os,sys
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict

#! for getting function norm2D
ROOT.gROOT.ProcessLine(".L misc-CINT-funcs.C")

NDTYP=2
EXP,SIM=range(NDTYP)
DTYP_NAME=["exp","sim"]
FIN=[[] for i in range(NDTYP)]
T=[[] for i in range(NDTYP)]
FOUT=[[] for i in range(NDTYP)]

NPRT=4
E,P,PIP,PIM=range(NPRT)
PRT_NAME=["e","p","pip","pim"]

NSCTR=6
SCTR_PHI_MIN=[-30,30,90, 150,210,270];
SCTR_PHI_MAX=[ 30,90,150,210,270,330];

#! Number of different monitoring histograms
NHST=1

#! Number of different view for monitoring hists(=LIN and NRM)
NVW=2
LIN,NRM=range(NVW)

#! Create structure H[DTYP][NHST][NPRT][NSCTR][NVW]
H=[]
for i in range(NDTYP):
	H.append([[[[[] for m in range(NVW)] for l in range(NSCTR)] for k in range(NPRT)] for j in range(NHST)])

#! Create structures:
#! + HST_NAME[NHST][NPRT]
#! + DRW_CMD[NHST][NPRT][NSCTR]
DRW_CMD=[]
HST_NAME=[]
for i in range(NHST):
	HST_NAME.append([[] for j in range(NPRT)])
        DRW_CMD.append([[[] for k in range(NSCTR)] for j in range(NPRT)])
#! Now fill DRW_CMD
for iprt in range(NPRT):
	HST_NAME[0][iprt]=("phi_%sVtheta_%s"
                           %(PRT_NAME[iprt],PRT_NAME[iprt]))
	for isctr in range(NSCTR):
		DRW_CMD[0][iprt][isctr]=("phi_%s:theta_%s>>%s(100,0,60,100,%d,%d)"
					  %(PRT_NAME[iprt],PRT_NAME[iprt],HST_NAME[0][iprt],
					    SCTR_PHI_MIN[isctr],SCTR_PHI_MAX[isctr])
                                        )
#! Print DRW_CMD
#for ihst in range(NHST):
#	for iprt in range(NPRT):
#        	for isctr in range(NSCTR):
#			print DRW_CMD[ihst][iprt][isctr]

#! Canvases for plotting NHST=C[NDTYP][NHST][NPRT][NVW]
C=[]
for i in range(NDTYP):
	C.append([[[[] for l in range(NVW)] for k in range(NPRT)] for j in range(NHST)])

def plot_fid(rctn,nentries=1000000000):
	#! Get FIN, T and make FOUT	
	if (rctn=='2pi'):
		#! Input objects
		FIN[EXP]=ROOT.TFile(("%s/defid.root"%os.environ['STUDY_EFID_2PI_DATADIR_EXP']))
		FIN[SIM]=ROOT.TFile(("%s/defid.root"%os.environ['STUDY_EFID_2PI_DATADIR_SIM']))	
		T[EXP]=FIN[EXP].Get("d2piR/tR")
		T[SIM]=FIN[SIM].Get("d2piR/tR")
		#! Output objects
		OUTDIR=("%s/hists"%os.environ['STUDY_EFID_2PI_DATADIR'])
		if not os.path.exists(OUTDIR):
			os.makedirs(OUTDIR)
		FOUT[EXP]=ROOT.TFile("%s/fexp.root"%(OUTDIR),"RECREATE");
		FOUT[SIM]=ROOT.TFile("%s/fsim.root"%(OUTDIR),"RECREATE");
		
        elif (rctn=='elast'):
		#! Input objects
		FIN[EXP]=ROOT.TFile(("%s/defid.root"%os.environ['STUDY_EFID_ELAST_DATADIR_EXP']))
                FIN[SIM]=ROOT.TFile(("%s/defid.root"%os.environ['STUDY_EFID_ELAST_DATADIR_SIM']))
		T[EXP]=FIN[EXP].Get("delast/t");
 		T[SIM]=FIN[SIM].Get("delast/t");
		#! Output objects
                OUTDIR=("%s/hists"%os.environ['STUDY_EFID_ELAST_DATADIR'])
                if not os.path.exists(OUTDIR):
                        os.makedirs(OUTDIR)
                FOUT[EXP]=ROOT.TFile("%s/fexp.root"%(OUTDIR),"RECREATE");
                FOUT[SIM]=ROOT.TFile("%s/fsim.root"%(OUTDIR),"RECREATE");	
	else:
		sys.exit("rctn=%s not recognized"%rctn)

		
	#! Start making plots H[DTYP][NHST][NPRT][NSCTR]
	c=ROOT.TCanvas() #! default canvas
	for idtyp in range(NDTYP):
		#if idtyp==EXP: continue
		for ihst in range(NHST):
			for iprt in range(NPRT):
				for isctr in range(NSCTR):
					if (iprt==E): sctr_cut=ROOT.TCut("sector==%d"%(isctr+1))
					else:         sctr_cut=ROOT.TCut("sector_%s==%d"%(PRT_NAME[iprt],isctr+1))
					top_cut=ROOT.TCut("")
					if (iprt==P or iprt==PIP or iprt==PIM): top_cut=ROOT.TCut("top==1")
					#! The following two line combine sctr_cut and top_cut using
					#! '+=' instead of '&&' because the later cannot be used in PyROOT
					#! Ref: https://root.cern.ch/phpBB3/viewtopic.php?t=3939
					cut=ROOT.TCut(sctr_cut)
					cut+=top_cut
					T[idtyp].Draw(DRW_CMD[ihst][iprt][isctr],cut,"colz",nentries)

					#! Store histogram
					ROOT.gStyle.SetOptStat("ne")
					htmp=ROOT.gDirectory.Get(HST_NAME[ihst][iprt])	
					H[idtyp][ihst][iprt][isctr][LIN]=htmp.Clone()
					H[idtyp][ihst][iprt][isctr][LIN].SetName("%s_%s_s%d"%(DTYP_NAME[idtyp],HST_NAME[ihst][iprt],isctr+1))
					H[idtyp][ihst][iprt][isctr][LIN].SetTitle("%s %s"%(HST_NAME[ihst][iprt],cut.GetTitle()))
					H[idtyp][ihst][iprt][isctr][NRM]=ROOT.norm2D(H[idtyp][ihst][iprt][isctr][LIN])


	#! Draw and Save H[DTYP][NHST][NPRT][NSCTR][NVW]
	#! Note, that the following code is optimized to:
	#! + Draw plots intereactively
	#! + Save to .root file and .jpg
	for idtyp in range(NDTYP):
		#if idtyp==EXP: continue
		#! For .jpg output
		outdir="%s/%s"%(OUTDIR,DTYP_NAME[idtyp])
		if not os.path.exists(outdir):
                        os.makedirs(outdir)		
		for ihst in range(NHST):
			for iprt in range(NPRT):
				#C[idtyp][ihst][iprt][LIN]=ROOT.TCanvas("%s_%s_lin"%(DTYP_NAME[idtyp],HST_NAME[ihst][iprt]))
				cname="%s_%s_lin"%(DTYP_NAME[idtyp],HST_NAME[ihst][iprt])
				C[idtyp][ihst][iprt][LIN]=ROOT.TCanvas(cname,cname)
				C[idtyp][ihst][iprt][LIN].Divide(3,2);
				#C[idtyp][ihst][iprt][NRM]=ROOT.TCanvas("%s_%s_nrm"%(DTYP_NAME[idtyp],HST_NAME[ihst][iprt]))
				cname="%s_%s_nrm"%(DTYP_NAME[idtyp],HST_NAME[ihst][iprt])
				C[idtyp][ihst][iprt][NRM]=ROOT.TCanvas(cname,cname)
				C[idtyp][ihst][iprt][NRM].Divide(3,2);
				for isctr in range(NSCTR):
					C[idtyp][ihst][iprt][LIN].cd(isctr+1)
					H[idtyp][ihst][iprt][isctr][LIN].Draw("colz")
					C[idtyp][ihst][iprt][NRM].cd(isctr+1)
					H[idtyp][ihst][iprt][isctr][NRM].Draw("colz")
				#! .root file output
				FOUT[idtyp].WriteTObject(C[idtyp][ihst][iprt][LIN])
				FOUT[idtyp].WriteTObject(C[idtyp][ihst][iprt][NRM])
				#! .jpg output
				C[idtyp][ihst][iprt][LIN].SaveAs("%s/%s.jpg"%(outdir,C[idtyp][ihst][iprt][LIN].GetName()))
				C[idtyp][ihst][iprt][NRM].SaveAs("%s/%s.jpg"%(outdir,C[idtyp][ihst][iprt][NRM].GetName()))
	
	#! Write FOUT
	#! Commented out because ALL objects are written to file!
	#FOUT[EXP].Write();
	#FOUT[SIM].Write();
	#! Close FOUT
	#! Commented out because it causes interactive displays to be empty?!
	#FOUT[EXP].Close();
        #FOUT[SIM].Close();

	#! If wanting to keep TCanvas open till program exits				
	if not ROOT.gROOT.IsBatch():
    		plt.show()
    		# wait for you to close the ROOT canvas before exiting
    		wait(True)

if __name__ == "__main__":	
	if len(sys.argv)==2:
		plot_fid(sys.argv[1])
	elif len(sys.argv)==3:
		plot_fid(sys.argv[1],int(sys.argv[2])) 
	
