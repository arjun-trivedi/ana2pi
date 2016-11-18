#!/usr/bin/python
from __future__ import division

import os,sys,time
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict

import numpy as np

'''
Usage: obtain_SFvp_cuts.py expt=<e1f/e16> [any gibberish(!) for debug] 
'''
#! Get args from user
if len(sys.argv)<2:
	sys.exit("Enter expt=e1f/e16")
else:
	expt=sys.argv[1]
	if expt!="e1f"and expt!="e16":
		sys.exit("Valid type for expt=e1f/e16 only")

DEBUG=False
if len(sys.argv)==3:
	DEBUG=True

print "Going to study SF for expt=",expt
if expt=="e1f":
	INDIR_EXP=os.path.join(os.environ['D2PIDIR_EXP'],"data_SF_100115")
	INDIR_SIM=os.path.join(os.environ['D2PIDIR_SIM'],"data_SF_100115")
	OUTDIR="%s/results_SFvp/"%os.environ['STUDY_EID_SF_DATADIR']
elif expt=="e16":
	INDIR_EXP=os.path.join(os.environ['D2PIDIR_EXP_E16'],"data_SF_010516")
        INDIR_SIM=os.path.join(os.environ['D2PIDIR_SIM_E16'],"data_SF_010516")
	OUTDIR="%s/results_SFvp/"%os.environ['STUDY_EID_SF_E16_DATADIR']
print "INDIR_EXP",INDIR_EXP
print "INDIR_SIM",INDIR_SIM
print "OUTDIR=",OUTDIR
#sys.exit()

PBINW=0.100 #! 100 MeV
# if len(sys.argv)>=2:
# 	PBINW=float(sys.argv[1])

NENTRIES=10000000000
# if len(sys.argv)>=3:
# 	NENTRIES=int(sys.argv[2])
print "PBINW,NENTRIES=",PBINW,"GeV",NENTRIES

#! Initial set up
#! setup momentum binning
PMIN=0.64
PMAX=5.00
PBIN_LE=np.arange(PMIN,PMAX,PBINW)
NPBIN=len(PBIN_LE)
print("NPBIN=",NPBIN)
time.sleep(3)

NDTYP=2
EXP,SIM=range(2)
DTYP_NAME=['exp','sim']

NSCTR=6

def make_fout_SFvp_pre_fits():
	FIN_2pi=[0,0]
	T_2pi=[0,0]
	FIN_2pi[EXP]=ROOT.TFile("%s/dSF_2pi.root"%INDIR_EXP)#os.environ['D2PIDIR_EXP'])
	FIN_2pi[SIM]=ROOT.TFile("%s/dSF_2pi.root"%INDIR_SIM)#os.environ['D2PIDIR_SIM'])
	T_2pi[EXP]=FIN_2pi[EXP].Get("d2piR/tR")
	T_2pi[SIM]=FIN_2pi[SIM].Get("d2piR/tR")

	if expt=="e1f":
		FIN_elast=[0,0]
		T_elast=[0,0]
		FIN_elast[EXP]=ROOT.TFile("%s/dSF_elast.root"%INDIR_EXP)#os.environ['D2PIDIR_EXP'])
		FIN_elast[SIM]=ROOT.TFile("%s/dSF_elast.root"%INDIR_SIM)#os.environ['D2PIDIR_SIM'])
		T_elast[EXP]=FIN_elast[EXP].Get("delast/cut/t")
		T_elast[EXP]=FIN_elast[EXP].Get("delast/cut/t")
	elif expt=="e16": #! no elast sim  and therefore no elast data for e16
		FIN_elast=[0]
                T_elast=[0]
		FIN_elast[EXP]=ROOT.TFile("%s/dSF_elast.root"%INDIR_EXP)#os.environ['D2PIDIR_EXP'])
		T_elast[EXP]=FIN_elast[EXP].Get("delast/cut/t")

	#! Make hSFvp[NDTYP][NSCTR]
	draw_cmd="etot/p:p>>hdcmd(100,0,5,100,0,0.5)"
	hSFvp=[[[] for isctr in range(NSCTR)] for idtyp in range(NDTYP)]
	for idtyp in range(NDTYP):
		for isctr in range(NSCTR):
			print "Making SFvp for",DTYP_NAME[idtyp],isctr+1
			cut_sector=ROOT.TCut("sector==%d"%(isctr+1))
			#! First draw hist from 2pi data
			T_2pi[idtyp].Draw(draw_cmd,cut_sector,"colz",NENTRIES)
			#! put hist in hSFvp[NDTYP][NSCTR]
			htmp=ROOT.gDirectory.Get("hdcmd")
			hname="h_%s_s%d"%(DTYP_NAME[idtyp],isctr+1)
			hSFvp[idtyp][isctr]=htmp.Clone(hname)

			if expt=="e1f" or (expt=="e16" and idtyp==0):
				#! Now draw hist from elast data
				T_elast[idtyp].Draw(draw_cmd,cut_sector,"colz",NENTRIES)
				#! Add hist to already created hSFvp[NDTYP][NSCTR]
				htmp=ROOT.gDirectory.Get("hdcmd")
				hSFvp[idtyp][isctr].Add(htmp)
			

	#! Plot and write hSFvp[NDTYP][NSCTR] to file
	ROOT.gStyle.SetOptStat("ne")
	C=[0,0]
	FOUT_SFvp=ROOT.TFile(fout_sfvp_pre_fits_name,"RECREATE")
	for idtyp in range(NDTYP):
		C[idtyp]=ROOT.TCanvas("c_%s"%DTYP_NAME[idtyp],"c_%s"%DTYP_NAME[idtyp])
		C[idtyp].Divide(3,2)
		for isctr in range(NSCTR):
			ROOT.gStyle.SetOptStat("ne")
			C[idtyp].cd(isctr+1)
			hSFvp[idtyp][isctr].Draw("colz")
			hSFvp[idtyp][isctr].Write()
	C[EXP].Write()
	C[SIM].Write()
	FOUT_SFvp.Close()

def obtain_SFcut_pars():
	"""
	Obtain parameters for function to make SF cut, for each DTYP and in each SCTR: fpSFvp[NDTYP][NSCTR][NMTHD][NFNC][NPAR], where
	+ NMTHD=number of different ways SF cuts pars in be obtained. Currently 2:
		1. Fit full range of SF in each pbin
		2. Fit around peak of SF in each pin
	+ NFNC=Number of functions: SF-high,-low and -mean
	+ NPAR=number of parameters of the fit function. Current 4 since fit function = pol3
	"""
	#! Open input root file
	FIN=ROOT.TFile(fout_sfvp_pre_fits_name)
	#! Create data structure fpSFvp: fpSFvp[NDTYP][NSCTR][NMTHD][NFNC][NPAR]
	NMTHD=2 #! fit (1) Full and (2) peak of SF in each pbin
	F,P=range(NMTHD)
	MTHD_NAME=['fullSF',"peakSF"]
	MTHD_CLR=['kBlack','kMagenta']
	MTHD_LINE_STYLE=[1,2] #1=regular, 2=dashed
	MTHD_MRKR_STYLE=['kOpenCircle','kOpenSquare'] #1=regular, 2=dashed

	#! Fit ranges for peakSF: peakSF_fit_range[NPRNG][NSCTR], where
	#! +NPRNG: Number of momentum ranges where peakSF fit range needs to be tuned
	PRNG=[[PMIN,0.94],[0.94,1.74],[1.74,4.04],[4.04,4.44],[4.44,PMAX]]
	NPRNG=len(PRNG)
	peakSF_fit_range=[[] for iprng in range(NPRNG)]
	peakSF_fit_range[0]=[[0.23,0.32],[0.26,0.34],[0.26,0.34],[0.26,0.34],[0.23,0.32],[0.23,0.32]]
	peakSF_fit_range[1]=[[0.26,0.32],[0.26,0.38],[0.26,0.36],[0.26,0.35],[0.26,0.32],[0.26,0.32]]
	peakSF_fit_range[2]=[[0.28,0.35],[0.32,0.42],[0.31,0.38],[0.28,0.38],[0.28,0.35],[0.28,0.35]]
	peakSF_fit_range[3]=[[0.30,0.36],[0.35,0.42],[0.32,0.39],[0.30,0.38],[0.28,0.36],[0.30,0.36]]
	peakSF_fit_range[4]=[[0.32,0.37],[0.34,0.40],[0.35,0.42],[0.32,0.38],[0.31,0.38],[0.30,0.35]]


	NFNC=3 #! h,m,l=(high,low,mean)
	H,L,M=range(NFNC)
	FNC_NAME=['high','low','mean']
	

	NPAR=4 #! using pol3 (1 more for constant!)

	fpSFvp=[[[[[0 for ipar in range(NPAR)] for icut in range(NFNC)] for imthd in range(NMTHD)] for isctr in range(NSCTR)] for idtyp in range(NDTYP)]
	
	#! + Define the 2 Gaus functions for the 2 NMTHD for fitting SF
	#! + This makes plotting, displaying and obtain the fits from the 2 Gaus functions simpler
	gf=[0,0]
	for imthd in range(NMTHD):
		gf[imthd]=ROOT.TF1(MTHD_NAME[imthd],"gaus")
		gf[imthd].SetLineColor(ROOT.gROOT.ProcessLine(MTHD_CLR[imthd]))
		#gf[imthd].SetLineStyle(MTHD_LINE_STYLE[imthd])
		
	for idtyp in range(NDTYP):
			if DEBUG and DTYP_NAME[idtyp]!="exp": continue

			fout_name="%s/fSFvp_post_fits_%s.root"%(OUTDIR,DTYP_NAME[idtyp])
			if DEBUG: 
				fout_name="%s/fSFvp_post_fits_%s_dbg.root"%(OUTDIR,DTYP_NAME[idtyp])
			FOUT=ROOT.TFile(fout_name,"RECREATE")

			#! [04-28-16] thesis-phase addition
			#! Create TCanvas to store in one canvas for all sectors "relevant data" for showing SF cut in thesis
			#! + "relevant data" created under each sector's loop iteration
			#! + Such work will typical in thesis-phase
			cname="cSFvp_thesis_%s"%(DTYP_NAME[idtyp])
			cSFvp_thesis=ROOT.TCanvas(cname,cname,2000,1500)	
			cSFvp_thesis.Divide(3,2)

			#! Create hcutSF[NMTHD][NFNC][NSCTR]
			#! [04-28-16] thesis-phase: The following is now created outside sctr loop so that all the relevant
                        #!                          objects to be stored in cSFvp_thesis have their own index
                        hcutSFvp=[[[[] for isctr in range(NSCTR)] for ifnc in range(NFNC)] for imthd in range(NMTHD)]
			for isctr in range(NSCTR):
				if DEBUG and (isctr+1)!=1: continue

				FOUT.mkdir("s%d"%(isctr+1)).cd()
				outdir_png="%s/%s_pngs/s%d"%(OUTDIR,DTYP_NAME[idtyp],isctr+1)
				if not os.path.exists(outdir_png):
					os.makedirs(outdir_png)

				#! Get hSFvp[NDTYP]
				hSFvp=FIN.Get("h_%s_s%d"%(DTYP_NAME[idtyp],isctr+1))
				hSFvp.SetTitle("sector=%d"%(isctr+1))

				#! Create hcutSF[NMTHD][NFNC]
				#! [04-28-16] thesis-phase: The following is now created outside sctr loop so that all the relevant
				#!                          objects to be stored in cSFvp_thesis have their own index
				#hcutSFvp=[[[] for ifnc in range(NFNC)] for imthd in range(NMTHD)]
				nbins=hSFvp.GetNbinsX()
				xmin=hSFvp.GetXaxis().GetXmin()
				xmax=hSFvp.GetXaxis().GetXmax()
				for imthd in range(NMTHD):
					for ifnc in range(NFNC):
						hname="%s_%s_s%d"%(MTHD_NAME[imthd],FNC_NAME[ifnc],isctr+1)
						hcutSFvp[imthd][ifnc][isctr]=ROOT.TH1F(hname,hname,nbins,xmin,xmax)
						#! aesthetics
						hcutSFvp[imthd][ifnc][isctr].SetMarkerStyle(ROOT.gROOT.ProcessLine(MTHD_MRKR_STYLE[imthd]))
						hcutSFvp[imthd][ifnc][isctr].SetMarkerColor(ROOT.gROOT.ProcessLine(MTHD_CLR[imthd]))
						
				for ipbin in range(NPBIN):
					pbin_min=round(PBIN_LE[ipbin],2)
					pbin_max=round(PBIN_LE[ipbin]+PBINW,2)
					if pbin_min>=4.84: continue #! These bins have no data in them

					hname="%s_s%d_pbin%d"%(DTYP_NAME[idtyp],isctr+1,ipbin+1)
					htitle="p=[%.2f GeV, %.2f GeV)"%(pbin_min,pbin_max)
					print "Making SF projection and fitting for",hname
					bin1=hSFvp.GetXaxis().FindBin(pbin_min)
					bin2=hSFvp.GetXaxis().FindBin(pbin_max)
					hSF=hSFvp.ProjectionY(hname,bin1,bin2)
					hSF.SetTitle(htitle)
					#hSF.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
					#! fit projection and get SF-cut-pars(pbin)
					if DTYP_NAME[idtyp]=="exp":
						#! fullSF
						if pbin_min>=0.84:
							hSF.Fit(gf[F].GetName(),"+0","",0.25,hSF.GetXaxis().GetXmax())
						else:
							hSF.Fit(gf[F].GetName(),"+0","")
						#! peakSF
						if pbin_min<0.94:
							sfmin=peakSF_fit_range[0][isctr][0]
							sfmax=peakSF_fit_range[0][isctr][1]
						if pbin_min>=0.94 and pbin_min<1.74:
							sfmin=peakSF_fit_range[1][isctr][0]
							sfmax=peakSF_fit_range[1][isctr][1]
							#hSF.Fit(gf[P].GetName(),"+0","",0.26,0.32)
						elif pbin_min>=1.74 and pbin_min<4.04:
							sfmin=peakSF_fit_range[2][isctr][0]
							sfmax=peakSF_fit_range[2][isctr][1]
							#hSF.Fit(gf[P].GetName(),"+0","",0.29,0.35)
						elif pbin_min>=4.04 and pbin_min<4.44:
							sfmin=peakSF_fit_range[3][isctr][0]
							sfmax=peakSF_fit_range[3][isctr][1]
							#hSF.Fit(gf[P].GetName(),"+0","",0.30,0.36)
						elif pbin_min>=4.44 and PMAX:
							sfmin=peakSF_fit_range[4][isctr][0]
							sfmax=peakSF_fit_range[4][isctr][1]
							#hSF.Fit(gf[P].GetName(),"+0","",0.32,0.37)
						# else:
						# 	hSF.Fit(gf[P].GetName(),"+0","",0.24,0.32)
						hSF.Fit(gf[P].GetName(),"+0","",sfmin,sfmax)
					else:#! i.e. DTYP="sim"
						#! fullSF
						hSF.Fit(gf[F].GetName(),"+0","")
						#! peakSF
						hSF.Fit(gf[P].GetName(),"+0","",0.2,0.25)

					#! Now get fit pars
					for imthd in range(NMTHD):
						ff=hSF.GetFunction(gf[imthd].GetName());
						if ff==None:
							mu,sg,mu_err,sg_err=0,0,0,0
						else:
							mu=ff.GetParameter(1);
							sg=ff.GetParameter(2);
							mu_err=ff.GetParError(1);
							sg_err=ff.GetParError(2);
						cut_h=mu+3*sg;
						cut_l=mu-3*sg;
						cut_h_err=np.sqrt(np.power(mu_err,2)+np.power(3*sg_err,2))
						cut_l_err=cut_h_err
						#! Fill SF-cut-pars(pbin) in hcutSF
						#! Note the bin used here is that contains pavg=(pbin_min+pbin_max)/2
						pavg=(pbin_min+pbin_max)/2
						bin=hcutSFvp[imthd][H][isctr].FindBin(pavg)
						hcutSFvp[imthd][H][isctr].SetBinContent(bin,cut_h)
						hcutSFvp[imthd][L][isctr].SetBinContent(bin,cut_l)
						hcutSFvp[imthd][M][isctr].SetBinContent(bin,mu)
						#! + PyROOT will not fit till hist bins have errors set!
						#! + If error bars are not set, it returns "Fit data is empty"
						hcutSFvp[imthd][H][isctr].SetBinError(bin,cut_h_err)
						hcutSFvp[imthd][L][isctr].SetBinError(bin,cut_l_err)
						hcutSFvp[imthd][M][isctr].SetBinError(bin,mu_err)
					#! Draw and save fitted hSF along with fits (as .png and in .root file) 
					#! plotting aesthetics
					ROOT.gStyle.SetOptStat("ne")
					ROOT.gStyle.SetOptFit(1111)
					ROOT.gStyle.SetStatX(0.9)
					ROOT.gStyle.SetStatW(0.1)
					ROOT.gStyle.SetStatH(0.1)
					c=ROOT.TCanvas("c_%s"%hSF.GetName(),"c_%s"%hSF.GetName())
					hSF.Draw()
					#! Axes title
                                	hSF.SetXTitle("Sampling Fraction (SF) [GeV]")
                                	c.SetLeftMargin(0.20)
                                	hSF.GetYaxis().SetTitleOffset(1.5)
                                	hSF.SetYTitle("N_{entries}")
					#! Get fit funcs and draw them
					gfF=hSF.GetFunction(gf[F].GetName())
					gfP=hSF.GetFunction(gf[P].GetName())
					if gfF!=None:
						gfF.Draw("same")
					if gfP!=None:
						gfP.Draw("same")
					#! legend
					l=ROOT.TLegend(0.8,0.5,0.9,0.6)#,"","NDC");
					l.AddEntry(gfF,MTHD_NAME[F])
					l.AddEntry(gfP,MTHD_NAME[P])
					l.Draw("same")
					#! .png and .pdf  save
					c.SaveAs("%s/c_pbin%02d.png"%(outdir_png,ipbin+1))
					c.SaveAs("%s/c_pbin%02d.pdf"%(outdir_png,ipbin+1))
					#! .root file save
					#hSF.Write()
					c.Write()
				#! Now fit hcutSFvp[NMTHD][NFNC][NSCTR], save fit pars in fpSFvp[NDTYP][NSCTR][NMTHD][NFNC][NPAR]
				for imthd in range(NMTHD):
					fitf=[0,0,0]
					for ifnc in range(NFNC):
						fitf[ifnc]=ROOT.TF1("fitf_%s"%FNC_NAME[ifnc],"pol3",PMIN,PMAX)
						fitf[ifnc].SetLineColor(ROOT.gROOT.ProcessLine(MTHD_CLR[imthd]))
						#fitf[ifnc].SetLineStyle(MTHD_LINE_STYLE[imthd])
					print "Fitting hcutSFvp for mthd:fnc:sctr=%s,%s,%d"%(MTHD_NAME[imthd],FNC_NAME[H],isctr+1)
					hcutSFvp[imthd][H][isctr].Fit(fitf[H].GetName(),"","",PMIN,PMAX)
					print "Fitting hcutSFvp for mthd:fnc:sctr=%s,%s,%d"%(MTHD_NAME[imthd],FNC_NAME[L],isctr+1)
					hcutSFvp[imthd][L][isctr].Fit(fitf[L].GetName(),"","",PMIN,PMAX)
					print "Fitting hcutSFvp for mthd:fnc:sctr=%s,%s,%d"%(MTHD_NAME[imthd],FNC_NAME[M],isctr+1)
					hcutSFvp[imthd][M][isctr].Fit(fitf[M].GetName(),"","",PMIN,PMAX)
					for ipar in range(NPAR):
						fpSFvp[idtyp][isctr][imthd][H][ipar]=fitf[H].GetParameter(ipar)
						fpSFvp[idtyp][isctr][imthd][L][ipar]=fitf[L].GetParameter(ipar)
						fpSFvp[idtyp][isctr][imthd][M][ipar]=fitf[M].GetParameter(ipar)
				#! Draw hSFvp, hcutSFvp (and their fits) on same canvas and save	
				c=ROOT.TCanvas("cSFvp","cSFvp")
				#! plotting aesthetics
				ROOT.gStyle.SetOptStat("ne")
				ROOT.gStyle.SetStatX(0.9);
				hSFvp.Draw("colz")
				#! Axes title
				hSFvp.SetXTitle("p [GeV]")
				c.SetLeftMargin(0.20)
				hSFvp.GetYaxis().SetTitleOffset(2.0)
				hSFvp.SetYTitle("Sampling Fraction (SF)")
				for imthd in range(NMTHD):
					hcutSFvp[imthd][H][isctr].Draw("P same")
					hcutSFvp[imthd][L][isctr].Draw("P same")
					hcutSFvp[imthd][M][isctr].Draw("P same")
				#! legend
				l=ROOT.TLegend(0.1,0.8,0.3,0.9)#,"","NDC");
				for imthd in range(NMTHD):
					l.AddEntry(hcutSFvp[imthd][H][isctr],MTHD_NAME[imthd],"p")
				l.Draw("same")
				#! save as .png
				c.SaveAs("%s/cSFvp.png"%(outdir_png))
				#! save in .root file
				c.Write()

				#! + [04-28-16] Also add this histograms to cSFvp_thesis
				print cSFvp_thesis.GetName()
				pad=cSFvp_thesis.cd(isctr+1)
				#! plotting aesthetics
                                #ROOT.gStyle.SetOptStat("ne")
				ROOT.gStyle.SetOptStat(0)
                                ROOT.gStyle.SetStatX(0.9);
                                hSFvp.Draw("colz")
				#! Axes title
                                hSFvp.SetXTitle("p [GeV]")
                                pad.SetLeftMargin(0.20)
                                hSFvp.GetYaxis().SetTitleOffset(2.0)
                                hSFvp.SetYTitle("Sampling Fraction (SF)")
                                for imthd in range(NMTHD):
                                        hcutSFvp[imthd][H][isctr].Draw("P same")
                                        hcutSFvp[imthd][L][isctr].Draw("P same")
                                        hcutSFvp[imthd][M][isctr].Draw("P same")
                                #! legend
				if (isctr+1==1):
                                	l_thesis=ROOT.TLegend(0.2,0.8,0.4,0.9)#,"","NDC");
                                	for imthd in range(NMTHD):
                                        	l_thesis.AddEntry(hcutSFvp[imthd][H][isctr],MTHD_NAME[imthd],"p")
                                	l_thesis.Draw("same")

			#! Write cSFvp_thesis in .root file and as .png and .pdf
			FOUT.cd()
			cSFvp_thesis.Write()
			outdir_png="%s/%s_pngs"%(OUTDIR,DTYP_NAME[idtyp])
                        if not os.path.exists(outdir_png):
                        	os.makedirs(outdir_png)
			cSFvp_thesis.SaveAs("%s/cSFvp.png"%(outdir_png))
			cSFvp_thesis.SaveAs("%s/cSFvp.pdf"%(outdir_png))
			
				
				
	#! Write fpSFvp[NDTYP][NSCTR][NMTHD][NFNC][NPAR] to file
	#print "fpSFvp[idtyp][isctr][imthd][H][0]=",fpSFvp[0][0][0][H][0]
	for idtyp in range(NDTYP):
		for imthd in range(NMTHD):
			outdir=("%s/cutpars"%OUTDIR)
			if not os.path.exists(outdir):
				os.makedirs(outdir)
			fcutpars=open("%s/%s_%s.txt"%(outdir,DTYP_NAME[idtyp],MTHD_NAME[imthd]), "w")
			for isctr in range(NSCTR):
				#! First write the high cut pars
				fcutpars.write("s%d:cut_h(pol3:p0:p1:p2:p3) %f %f %f %f\n"%(isctr+1,
					fpSFvp[idtyp][isctr][imthd][H][0],fpSFvp[idtyp][isctr][imthd][H][1],
					fpSFvp[idtyp][isctr][imthd][H][2],fpSFvp[idtyp][isctr][imthd][H][3]))
				#! Now write the low cut pars
				fcutpars.write("s%d:cut_l(pol3:p0:p1:p2:p3) %f %f %f %f\n"%(isctr+1,
					fpSFvp[idtyp][isctr][imthd][L][0],fpSFvp[idtyp][isctr][imthd][L][1],
					fpSFvp[idtyp][isctr][imthd][L][2],fpSFvp[idtyp][isctr][imthd][L][3]))
			fcutpars.close()
	#fcutpars.write("Purchase Amount: %s" % TotalAmount)
	

#! Start of main program
#OUTDIR="%s/results_SFvp/"%os.environ['STUDY_EID_SF']
#OUTDIR="%s/results_SFvp/"%os.environ['STUDY_EID_SF_DATADIR']
if not os.path.exists(OUTDIR):
	os.makedirs(OUTDIR)
#! 'touch' OUTDIR to show it has been worked upon
os.system('touch %s'%OUTDIR)

fout_sfvp_pre_fits_name="%s/fSFvp_pre_fits_exp_sim.root"%OUTDIR

if not os.path.exists(fout_sfvp_pre_fits_name):
	make_fout_SFvp_pre_fits()
	obtain_SFcut_pars()
else:
	print fout_sfvp_pre_fits_name,"exists"
	obtain_SFcut_pars()
	
if not ROOT.gROOT.IsBatch():
	plt.show()
	# wait for you to close the ROOT canvas before exiting
	wait(True)

#if __name__ == "__main__":	
#	if len(sys.argv)==2:
#		plot_fid(sys.argv[1])
#	elif len(sys.argv)==3:
#		plot_fid(sys.argv[1],int(sys.argv[2])) 
	
