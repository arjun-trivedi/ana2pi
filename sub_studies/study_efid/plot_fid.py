#!/usr/bin/python
from __future__ import division

import os,sys
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict

import numpy as np

#! for getting function norm2D
ROOT.gROOT.ProcessLine(".L misc-CINT-funcs.C")
#! for Evan's e- fiducial cuts
ROOT.gROOT.ProcessLine(".L fidfuncs.C");
#1 for old e- fiducial cuts
ROOT.gROOT.ProcessLine(".L cuts_CINT.h");
#! For Evan's Hadron fiducial cuts
ROOT.gROOT.ProcessLine(".L fidfuncs_hadrons.C")
#! For old Hadron fiducial cuts
ROOT.gROOT.ProcessLine(".L fidfuncs_hadrons_e16.C")

def plot_fid(rctn,nentries=1000000000):
	NDTYP=2
	EXP,SIM=range(NDTYP)
	DTYP_NAME=["exp","sim"]
	FIN=[[] for i in range(NDTYP)]
	T=[[] for i in range(NDTYP)]
	FOUT=[[] for i in range(NDTYP)]

	if rctn=='2pi':
		NPRT=4
		E,P,PIP,PIM=range(NPRT)
		PRT_NAME=["e","p","pip","pim"]
		PRT_ID=[11,2212,211,-211]
	elif rctn=='elast':
		NPRT=1
		E=0
		PRT_ID[11]
		PRT_NAME=["e"]
	else:
		sys.exit("rctn=%s not recognized"%rctn)

	NSCTR=6
	SCTR_PHI_MIN=[-30,30,90, 150,210,270];
	SCTR_PHI_MAX=[ 30,90,150,210,270,330];

	#! Create p bins for NPRT
	PBIN_MIN=[[] for i in range(NPRT)]
	PBIN_MAX=[[] for i in range(NPRT)]
	NPBIN=[0 for i in range(NPRT)]
	#! binw=200MeV for e-
	PBIN_MIN[E]=np.arange(2.0,5.0,0.2)
	PBIN_MAX[E]=np.arange(2.2,5.2,0.2)
	NPBIN[E]=len(PBIN_MIN[E])
	#! 1 bin for p,pip,pim
	for iprt in range(NPRT):
		if (iprt==E): continue
		PBIN_MIN[iprt]=[0]
		PBIN_MAX[iprt]=[10]
		NPBIN[iprt]=len(PBIN_MIN[iprt])
	#! Debug print p-bin information
	# for iprt in range(NPRT):
	#     print "pbins for",PRT_NAME[iprt]
	#     print PBIN_MIN[iprt]
	#     print PBIN_MAX[iprt]
	#     print NPBIN[iprt]



	#! Number of different monitoring histograms
	NHST=1
	HST_NAME=['phiVtheta']

	#! Create structure H[DTYP][NHST][NPRT][NSCTR][NPBIN][NVW]
	H=[]
	#! First create H[DTYP][NHST][NPRT][NSCTR]
	#! Number of different monitoring histograms
	NHST=1

	#! Number of different view for monitoring hists(=LIN and NRM)
	NVW=2
	VW_NAME=['lin','nrm']
	LIN,NRM=range(NVW)

	#! Create structure H[DTYP][NHST][NPRT][NSCTR][NPBIN_<PRT>][NVW]
	#! First create H[DTYP][NHST][NPRT][NSCTR]
	H=[[[[[] for l in range(NSCTR)] for k in range(NPRT)] for j in range(NHST)] for i in range(NDTYP)]
	#! Now add [NPBIN_<NPRT>][NVW]
	for idtyp in range(NDTYP):
		for ihst in range(NHST):
			for iprt in range(NPRT):
				for isctr in range(NSCTR):
					H[idtyp][ihst][iprt][isctr]=[[[] for ivw in range(NVW)] for ipbin in range(NPBIN[iprt])]
	#! Use the following to debug H structure
	# for idtyp in range(NDTYP):
	# 	for ihst in range(NHST):
	# 		for iprt in range(NPRT):
	# 			for isctr in range(NSCTR):
	# 				for ipbin in range(NPBIN[iprt]):
	# 					for ivw in range(NVW):
	# 						H[idtyp][ihst][iprt][isctr][ipbin][ivw]=VW_NAME[ivw]
	# 						H[idtyp][ihst][iprt][isctr][ipbin][ivw]=VW_NAME[ivw]
	# print H
	#print H[0][0][0][0][0][0]
	
	
	#! Create structures:
	#! + DRW_CMD_HST_NAME[NHST][NPRT]
	#! + DRW_CMD[NHST][NPRT][NSCTR]
	DRW_CMD_HST_NAME=[[[] for j in range(NPRT)] for i in range(NHST)]
	DRW_CMD=[[[[] for k in range(NSCTR)] for j in range(NPRT)] for i in range(NHST)]
	# #! Now add [NPBIN_<PRT>]
	# for ihst in range(NHST):
	# 	for iprt in range(NPRT):
	# 		for isctr in range(NSCTR):
	# 				DRW_CMD[ihst][iprt][isctr]=[[] for ipbin in range(NPBIN[iprt])]
	#! Now fill DRW_CMD
	for iprt in range(NPRT):
		DRW_CMD_HST_NAME[0][iprt]=("phi_%sVtheta_%s"%(PRT_NAME[iprt],PRT_NAME[iprt]))
		for isctr in range(NSCTR):
			if rctn=="2pi":
				DRW_CMD[0][iprt][isctr]=("phi_%s:theta_%s>>%s(200,0,70,200,%d,%d)"
								     	 %(PRT_NAME[iprt],PRT_NAME[iprt],DRW_CMD_HST_NAME[0][iprt],
					    	             SCTR_PHI_MIN[isctr],SCTR_PHI_MAX[isctr])
                	                )
			elif rctn=='elast':
				DRW_CMD[0][iprt][isctr]=("phi1:theta1>>%s(200,0,70,200,%d,%d)"
								     	 %(DRW_CMD_HST_NAME[0][iprt],
					    	             SCTR_PHI_MIN[isctr],SCTR_PHI_MAX[isctr])
                	                )

	#! Debug: Print DRW_CMD
	# for ihst in range(NHST):
	# 	for iprt in range(NPRT):
	# 		for isctr in range(NSCTR):
	# 				print DRW_CMD[ihst][iprt][isctr]

	# #! Canvases for plotting: C[NDTYP][NHST][NPRT][NPBIN_<PRT>][NVW]
	C=[[[[] for k in range(NPRT)] for j in range(NHST)] for i in range(NDTYP)]
	#! Now add [NPBIN_<PRT>][NVW]
	for idtyp in range(NDTYP):
		for ihst in range(NHST):
			for iprt in range(NPRT):
				C[idtyp][ihst][iprt]=[[[] for ivw in range(NVW)] for ipbin in range(NPBIN[iprt])]
	#! Now add [NPBIN_<PRT>][NVW]
	# C=[]
	# for i in range(NDTYP):
	# 	C.append([[[[] for l in range(NVW)] for k in range(NPRT)] for j in range(NHST)])
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

		
	#! Start making plots H[DTYP][NHST][NPRT][NSCTR][NPBIN_<PRT>][NVW]
	c=ROOT.TCanvas() #! default canvas
	for idtyp in range(NDTYP):
		#if idtyp==EXP: continue
		for ihst in range(NHST):
			for iprt in range(NPRT):
				for isctr in range(NSCTR):
					# if (iprt==E): sctr_cut=ROOT.TCut("sector==%d"%(isctr+1))
					# else:         sctr_cut=ROOT.TCut("sector_%s==%d"%(PRT_NAME[iprt],isctr+1))
					if rctn=='2pi':
						sctr_cut=ROOT.TCut("sector_%s==%d"%(PRT_NAME[iprt],isctr+1))
					elif rctn=='elast':
						sctr_cut=ROOT.TCut("sector==%d"%(isctr+1))
					top_cut=ROOT.TCut("")
					if rctn=='2pi':
						if (iprt==P or iprt==PIP or iprt==PIM): top_cut=ROOT.TCut("top==1")
					for ipbin in range(NPBIN[iprt]):
						pmin=PBIN_MIN[iprt][ipbin]
						pmax=PBIN_MAX[iprt][ipbin]
						p_cut=ROOT.TCut("p_%s>=%.2f && p_%s<=%.2f"%(PRT_NAME[iprt],pmin,PRT_NAME[iprt],pmax))
						#! The following two line combine sctr_cut and top_cut using
						#! '+=' instead of '&&' because the later cannot be used in PyROOT
						#! Ref: https://root.cern.ch/phpBB3/viewtopic.php?t=3939
						cut=ROOT.TCut(sctr_cut)
						cut+=top_cut
						cut+=p_cut
						print DRW_CMD[ihst][iprt][isctr]
						T[idtyp].Draw(DRW_CMD[ihst][iprt][isctr],cut,"colz",nentries)

						#! Store histogram
						ROOT.gStyle.SetOptStat("ne")
						htmp=ROOT.gDirectory.Get(DRW_CMD_HST_NAME[ihst][iprt])	
						#print idtyp,ihst,iprt,isctr,ipbin,LIN
						#print htmp.GetName()
						H[idtyp][ihst][iprt][isctr][ipbin][LIN]=htmp.Clone()
						H[idtyp][ihst][iprt][isctr][ipbin][LIN].SetName("%s_s%d_pbin%d"%(DRW_CMD_HST_NAME[ihst][iprt],isctr+1,ipbin+1))
						H[idtyp][ihst][iprt][isctr][ipbin][LIN].SetTitle("%s"%(cut.GetTitle()))
						H[idtyp][ihst][iprt][isctr][ipbin][NRM]=ROOT.norm2D(H[idtyp][ihst][iprt][isctr][ipbin][LIN])


	#! Draw and Save H[DTYP][NHST][NPRT][NSCTR][NPBIN_<PRT>][NVW] 
	#! Note, that the following code is optimized to:
	#! + Draw plots intereactively
	#! + Save to .root file and .jpg
	for idtyp in range(NDTYP):
		#if idtyp==EXP: continue
		#! For .jpg output
		#outdir="%s/%s"%(OUTDIR,DTYP_NAME[idtyp])
		#if not os.path.exists(outdir):
                #        os.makedirs(outdir)		
		for ihst in range(NHST):
			for iprt in range(NPRT):
				#! outdir for .jpg
				outdir="%s/%s/%s/%s/"%(OUTDIR,DTYP_NAME[idtyp],HST_NAME[ihst],PRT_NAME[iprt])
				if not os.path.exists(outdir):
					os.makedirs(outdir)
				#! outdir for .root
				outdir_root="%s/%s/"%(HST_NAME[ihst],PRT_NAME[iprt])
				if FOUT[idtyp].GetDirectory(HST_NAME[ihst])==None:
					FOUT[idtyp].mkdir(HST_NAME[ihst])
				hst_dir=FOUT[idtyp].GetDirectory(HST_NAME[ihst])
				if hst_dir.GetDirectory(PRT_NAME[iprt])==None:
					hst_dir.mkdir(PRT_NAME[iprt])
				prt_dir=hst_dir.GetDirectory(PRT_NAME[iprt])

				for ipbin in range(NPBIN[iprt]):
					#C[idtyp][ihst][iprt][LIN]=ROOT.TCanvas("%s_%s_lin"%(DTYP_NAME[idtyp],DRW_CMD_HST_NAME[ihst][iprt]))
					#cname="%s_%s_pbin%d_lin"%(DTYP_NAME[idtyp],DRW_CMD_HST_NAME[ihst][iprt],ipbin+1)
					cname="%s_pbin%02d_lin"%(HST_NAME[ihst],ipbin+1)
					C[idtyp][ihst][iprt][ipbin][LIN]=ROOT.TCanvas(cname,cname,1300,800)
					C[idtyp][ihst][iprt][ipbin][LIN].Divide(3,2);
					#C[idtyp][ihst][iprt][NRM]=ROOT.TCanvas("%s_%s_nrm"%(DTYP_NAME[idtyp],DRW_CMD_HST_NAME[ihst][iprt]))
					#cname="%s_%s_pbin%d_nrm"%(DTYP_NAME[idtyp],DRW_CMD_HST_NAME[ihst][iprt],ipbin+1)
					cname="%s_pbin%02d_nrm"%(HST_NAME[ihst],ipbin+1)
					C[idtyp][ihst][iprt][ipbin][NRM]=ROOT.TCanvas(cname,cname,1300,800)
					C[idtyp][ihst][iprt][ipbin][NRM].Divide(3,2);
					for isctr in range(NSCTR):
						if PRT_NAME[iprt]=="e":
							p=(PBIN_MIN[iprt][ipbin]+PBIN_MAX[iprt][ipbin])/2
							f=ROOT.fPhiFid_mod(p,isctr+1);
							f_old=ROOT.GetFidF_mod(11,p,isctr+1);
							f_old.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
						elif PRT_NAME[iprt]=='p' or PRT_NAME[iprt]=='pip':
							f_h_exp=ROOT.fPhiFid_hdrn_h_mod(PRT_ID[iprt],"exp",isctr+1,1);
							f_l_exp=ROOT.fPhiFid_hdrn_l_mod(PRT_ID[iprt],"exp",isctr+1,1);
							f_h_exp.SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
							f_l_exp.SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
							f_h_sim=ROOT.fPhiFid_hdrn_h_mod(PRT_ID[iprt],"sim",isctr+1,1);
							f_l_sim=ROOT.fPhiFid_hdrn_l_mod(PRT_ID[iprt],"sim",isctr+1,1);
							f_h_sim.SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
							f_l_sim.SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
							#! t0
							lt0_exp=ROOT.ROOT.lt0(PRT_ID[iprt],"exp");
							lt0_sim=ROOT.ROOT.lt0(PRT_ID[iprt],"sim");
							lt0_exp.SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
							lt0_sim.SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
							lt0_exp.SetLineWidth(3)
							lt0_sim.SetLineWidth(3)
							#! old fid cuts
							f_old_h=ROOT.fPhiFid_e16_hdrn_h_mod(isctr+1)
							f_old_l=ROOT.fPhiFid_e16_hdrn_l_mod(isctr+1)
							f_old_h.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
							f_old_l.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))


						C[idtyp][ihst][iprt][ipbin][LIN].cd(isctr+1)
						H[idtyp][ihst][iprt][isctr][ipbin][LIN].Draw("colz")
						if PRT_NAME[iprt]=="e":
							f.Draw("same")
							f_old.Draw("same")
						elif PRT_NAME[iprt]=='p' or PRT_NAME[iprt]=='pip':
							f_h_exp.Draw("same")
							f_l_exp.Draw("same")
							f_h_sim.Draw("same")
							f_l_sim.Draw("same")
							#! Set Y1 and Y2 for lt0 and draw
							ymin=H[idtyp][ihst][iprt][isctr][ipbin][LIN].GetYaxis().GetXmin()
							ymax=H[idtyp][ihst][iprt][isctr][ipbin][LIN].GetYaxis().GetXmax()
							lt0_exp.SetY1(ymin)
							lt0_exp.SetY2(ymax)
							lt0_sim.SetY1(ymin)
							lt0_sim.SetY2(ymax)
							lt0_exp.Draw("same")
							lt0_sim.Draw("same")
							#! old cuts
							f_old_h.Draw("same")
							f_old_l.Draw("same")

						C[idtyp][ihst][iprt][ipbin][NRM].cd(isctr+1)
						H[idtyp][ihst][iprt][isctr][ipbin][NRM].Draw("colz")
						if PRT_NAME[iprt]=="e":
							f.Draw("same")
							f_old.Draw("same")
						elif PRT_NAME[iprt]=='p' or PRT_NAME[iprt]=='pip':
							f_h_exp.Draw("same")
							f_l_exp.Draw("same")
							f_h_sim.Draw("same")
							f_l_sim.Draw("same")
							#! Set Y1 and Y2 for lt0 and draw
							ymin=H[idtyp][ihst][iprt][isctr][ipbin][NRM].GetYaxis().GetXmin()
							ymax=H[idtyp][ihst][iprt][isctr][ipbin][NRM].GetYaxis().GetXmax()
							lt0_exp.SetY1(ymin)
							lt0_exp.SetY2(ymax)
							lt0_sim.SetY1(ymin)
							lt0_sim.SetY2(ymax)
							lt0_exp.Draw("same")
							lt0_sim.Draw("same")
							#! old cuts
							f_old_h.Draw("same")
							f_old_l.Draw("same")
					#! .root file output
					#FOUT[idtyp].WriteTObject(C[idtyp][ihst][iprt][ipbin][LIN])
					#FOUT[idtyp].WriteTObject(C[idtyp][ihst][iprt][ipbin][NRM])
					prt_dir.WriteTObject(C[idtyp][ihst][iprt][ipbin][LIN])
					prt_dir.WriteTObject(C[idtyp][ihst][iprt][ipbin][NRM])
					#! .jpg output
					C[idtyp][ihst][iprt][ipbin][LIN].SaveAs("%s/%s.jpg"%(outdir,C[idtyp][ihst][iprt][ipbin][LIN].GetName()))
					C[idtyp][ihst][iprt][ipbin][NRM].SaveAs("%s/%s.jpg"%(outdir,C[idtyp][ihst][iprt][ipbin][NRM].GetName()))
					C[idtyp][ihst][iprt][ipbin][LIN].Close()
					C[idtyp][ihst][iprt][ipbin][NRM].Close()
	
# 	#! Write FOUT
# 	#! Commented out because ALL objects are written to file!
# 	#FOUT[EXP].Write();
# 	#FOUT[SIM].Write();
# 	#! Close FOUT
# 	#! Commented out because it causes interactive displays to be empty?!
# 	#FOUT[EXP].Close();
#         #FOUT[SIM].Close();

# 	#! If wanting to keep TCanvas open till program exits				
# 	if not ROOT.gROOT.IsBatch():
# 			plt.show()
# 			# wait for you to close the ROOT canvas before exiting
# 			wait(True)

if __name__ == "__main__":	
	if len(sys.argv)==2:
		plot_fid(sys.argv[1])
	elif len(sys.argv)==3:
		plot_fid(sys.argv[1],int(sys.argv[2])) 
	
