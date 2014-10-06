from __future__ import division
import ROOT
from rootpy.io import root_open, DoesNotExist

from collections import OrderedDict

import os,sys

import numpy as np
import matplotlib.pyplot as plt

import math

import atlib as atlib

from proc_yields import H5_DIM,VARS #! for H5_DIM 

ROOT.gROOT.ProcessLine(".L THnTool.C+")
from ROOT import THnTool

# Tools 
thntool=THnTool()

class DispYields:
	def __init__(self,simnum='siml'):
		self.SIM_NUM=simnum
		#self.Q2W=q2w
		# self.FEXP=root_open('$HOME/ongoing/mem_test/exp/new-h8-bng/yield_exp.root')
		# self.FSIM=root_open('$HOME/ongoing/mem_test/sim/new-h8-bng/yield_sim.root')
		# self.OUTDIR='/home/trivedia/ongoing/mem_test/obs'
		self.FEXP=root_open(os.path.join(os.environ['OBSDIR'],self.SIM_NUM,'yield_exp.root'))
		self.FEXP_HEL=root_open(os.path.join(os.environ['OBSDIR'],self.SIM_NUM,'yield_exp_hel.root'))
		self.FSIM=root_open(os.path.join(os.environ['OBSDIR'],self.SIM_NUM,'yield_sim.root'))
		self.OUTDIR=os.path.join(os.environ['OBSDIR'],self.SIM_NUM)
		if not os.path.exists(self.OUTDIR):
			sys.exit("%s does not exist!"%self.OUTDIR)

		self.VSTS=[1,2,3]

		#! Create dictionary for variable_names(vst,var) where
		#! + vst in VSTS=[1,2,3]
		#! + var in VARS=['M1','M2','THETA','PHI','ALPHA']
		#! NOTE: Should be same as DataAna::MakeYields()
		self.VAR_NAMES={(1,VARS[0]):"M_{p#pi^{+}}",(1,VARS[1]):"M_{#pi^{+}#pi^{-}}",
                        (1,VARS[2]):"#theta_{#pi^{-}}",(1,VARS[3]):"#phi_{#pi^{-}}",
                        (1,VARS[4]):"#alpha_{[p_{f}#pi^{+}][p#pi^{-}]}",
                        (2,VARS[0]):"M_{p#pi^{+}}",(2,VARS[1]):"M_{#pi^{+}#pi^{-}}",
                        (2,VARS[2]):"#theta_{p}",(2,VARS[3]): "#phi_{p}",
                        (2,VARS[4]):"#alpha_{[#pi^{+}#pi^{-}][pp_{f}]}",
                        (3,VARS[0]):"M_{p#pi^{+}}",(3,VARS[1]):"M_{p#pi^{-}}",
                        (3,VARS[2]):"#theta_{#pi^{+}}",(3,VARS[3]): "#phi_{#pi^{+}}",
                        (3,VARS[4]):"#alpha_{[p_{f}#pi^{-}][p#pi^{+}]}"}

        #self.VAR_UNIT_NAMES={VARS[0]:"[GeV]",VARS[1]:"[GeV]",VARS[2]:"[#degree]",VARS[3]:"[#degree]",VARS[4]:"[#degree]"}

	def plot_obs_1D(self,hVST1,hVST2,hVST3,view="q2_evltn",dtyp='EXP',seq='F'):
		"""
		+ Plot 1D-Obs as per "view", where:
			+ view="full_ana" (dtyp,seq not required to be specified)
			+ view="q2_evltn" (dtyp and seq need to specific)
		"""

		if view=="q2_evltn":
			print "Going to display 1D-Obs: q2-evol,%s,%s"%(dtyp,seq)
			outdir=os.path.join(self.OUTDIR_OBS_1D,"Q2_Evolution_%s_%s"%(dtyp,seq))
			if not os.path.exists(outdir):
				os.makedirs(outdir)
		elif view=="full_ana":
			print "Going to display 1D-Obs: full_ana"
			outdir=os.path.join(self.OUTDIR_OBS_1D,"full_ana")
			if not os.path.exists(outdir):
				os.makedirs(outdir)
		else:
			sys.exit("view=%s not recognized. Exiting."%view)

		#! Get all q2- and w-bins from the keys of hVSTX
		q2bins_le,wbins_le=[],[]
		for k in hVST1.keys():
			q2bins_le.append(k[0])
			wbins_le.append(k[1])
		#! Keep only unique entries
		q2bins_le=set(q2bins_le) #! Keep only unique entries
		q2bins_le=list(q2bins_le) #! convert 'set' back to 'list'
		q2bins_le=sorted(q2bins_le)
		wbins_le=set(wbins_le) #! Keep only unique entries
		wbins_le=list(wbins_le) #! convert 'set' back to 'list'
		wbins_le=sorted(wbins_le)
		print "going to plot 1D obs for:"
		print "Q2:"
		print q2bins_le
		print "W:"
		print wbins_le

		#! Set up some plotting related styles and aesthetics 
		if view=="q2_evltn":
			colors=["kRed","kOrange","kYellow","kGreen+3","kGreen","kCyan","kBlue","kMagenta"]
			coll=[]
			for iq2bin in range(len(q2bins_le)):
				coll.append(ROOT.gROOT.ProcessLine(colors[iq2bin]))
		elif view=="full_ana":
			coll={('EXP','C'):ROOT.gROOT.ProcessLine("kCyan"),
			      ('EXP','F'):ROOT.gROOT.ProcessLine("kBlue"),
			      ('EXP','H'):ROOT.gROOT.ProcessLine("kBlack"),
			      ('SIM','F'):ROOT.gROOT.ProcessLine("kRed"),
			      ('SIM','T'):ROOT.gROOT.ProcessLine("kGreen")}
		#! Label histograms
		q2wbintitle,wbintitle=self.label_hist_obs1D(hVST1,hVST2,hVST3)
		for k in wbintitle:
			print k,wbintitle[k]
		#! Stats Box
		ROOT.gStyle.SetOptStat(0)

		# ROOT.gStyle.SetLabelSize(0.5,"t")
		# ROOT.gStyle.SetTitleSize(0.5,"t")
		#ROOT.gStyle.SetPaperSize(20,26);
		ROOT.gStyle.SetPadTopMargin(0.15)#(0.05);
		#ROOT.gStyle.SetPadRightMargin(0.09)#(0.05);
		#ROOT.gStyle.SetPadBottomMargin(0.20)#(0.16);
		#ROOT.gStyle.SetPadLeftMargin(0.15)#(0.12);

		ROOT.gStyle.SetTitleW(10)# //title width 
		ROOT.gStyle.SetTitleFontSize(20)# //title width 
		ROOT.gStyle.SetTitleH(0.15)# //title height 
		ROOT.gStyle.SetTitleY(1)# //title Y location 

		#!get rid of X error bars and y error bar caps
		ROOT.gStyle.SetErrorX(0.001);



		ivars=[0,1,2]
		#! 3x3 TCanvas-pads as per VSTs: VST1=1,2,3; VST2=2,5,6; VST3=3,6,9
		pads=[[1,4,7],[2,5,8],[3,6,9]]
		map_padnums_vars=[zip(pads[0],ivars),zip(pads[1],ivars),zip(pads[2],ivars)]
		for wbin in wbins_le:
			if view=="q2_evltn" and wbintitle.has_key((wbin)):
				c=ROOT.TCanvas("c","c",1000,1000)
				pad_t=ROOT.TPad("pad_t","Title pad",0.05,0.97,0.95,1.00)
				#pad_t.SetFillColor(11)
  				pad_p=ROOT.TPad("pad_p","Plots pad",0.01,0.01,0.99,0.95);
  				pad_t.Draw()
  				pad_p.Draw()
  				pad_t.cd()
  				pt=ROOT.TPaveText(.05,.1,.95,.8)
  				pt.AddText("W bin=%s"%wbintitle[wbin])
  				pt.Draw()
 				#pad_p.cd()
				pad_p.Divide(3,3)
			for iq2bin,q2bin in enumerate(q2bins_le):
				if view=="full_ana" and q2wbintitle.has_key((q2bin,wbin)):
					c=ROOT.TCanvas("c","c",1000,1000)
					pad_t=ROOT.TPad("pad_t","Title pad",0.05,0.97,0.95,1.00)
					#pad_t.SetFillColor(11)
  					pad_p=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
  					pad_t.Draw()
  					pad_p.Draw()
  					pad_t.cd()
  					pt=ROOT.TPaveText(.05,.1,.95,.8)
  					pt.AddText("Q2_W bin=%s"%q2wbintitle[q2bin,wbin])
  					pt.Draw()
  					#pad_p.Draw()
  					#pad_p.cd()
					pad_p.Divide(3,3)
				print "Plotting h1D for w=%0.3f,q2=%0.2f"%(wbin,q2bin)
				for ivst,vst in enumerate(self.VSTS):
					if   vst==1: h=hVST1
					elif vst==2: h=hVST2
					elif vst==3: h=hVST3
					if h.has_key((q2bin,wbin,dtyp,seq)):
						# pads=padsl[ivst]
						# map_var_pad=zip(pads,ivars)
						for m in map_padnums_vars[ivst]:
							padnum=m[0]
							ivar=m[1]
							pad=pad_p.cd(padnum)
							if view=="q2_evltn":
								drawopt="same"
								if iq2bin==0:
									drawopt=""
									h[q2bin,wbin,dtyp,seq][ivar].Sumw2()
									h[q2bin,wbin,dtyp,seq][ivar].Scale(0.50)
									h[q2bin,wbin,dtyp,seq][ivar].SetMinimum(0)
								h[q2bin,wbin,dtyp,seq][ivar].SetMarkerColor(coll[iq2bin])
								h[q2bin,wbin,dtyp,seq][ivar].SetLineColor(coll[iq2bin])
								h[q2bin,wbin,dtyp,seq][ivar].Draw(drawopt)
							elif view=="full_ana":
								h[q2bin,wbin,'EXP','C'][ivar].SetMarkerColor(coll['EXP','C'])
								h[q2bin,wbin,'EXP','F'][ivar].SetMarkerColor(coll['EXP','F'])
								h[q2bin,wbin,'EXP','H'][ivar].SetMarkerColor(coll['EXP','H'])
								h[q2bin,wbin,'SIM','F'][ivar].SetMarkerColor(coll['SIM','F'])
								h[q2bin,wbin,'SIM','T'][ivar].SetMarkerColor(coll['SIM','T'])
								h[q2bin,wbin,'EXP','C'][ivar].SetLineColor(coll['EXP','C'])
								h[q2bin,wbin,'EXP','F'][ivar].SetLineColor(coll['EXP','F'])
								h[q2bin,wbin,'EXP','H'][ivar].SetLineColor(coll['EXP','H'])
								h[q2bin,wbin,'SIM','F'][ivar].SetLineColor(coll['SIM','F'])
								h[q2bin,wbin,'SIM','T'][ivar].SetLineColor(coll['SIM','T'])
								h[q2bin,wbin,'SIM','T'][ivar].SetMarkerStyle(ROOT.gROOT.ProcessLine("kPlus"))
								
								hexp=h[q2bin,wbin,'EXP','F'][ivar].DrawNormalized("",1000)
								hsim=h[q2bin,wbin,'SIM','F'][ivar].DrawNormalized("sames",1000)
								hexp.SetMinimum(0.)
								hsim.SetMinimum(0.)
								maximum=hexp.GetMaximum()
								if hsim.GetMaximum()>hexp.GetMaximum():
									maximum=hsim.GetMaximum()
								hexp.SetMaximum(maximum+10)
								hsim.SetMaximum(maximum+10)
								hF=hexp.Clone()
								hF.Divide(h[q2bin,wbin,'EXP','F'][ivar])
								h[q2bin,wbin,'EXP','C'][ivar].Multiply(hF)
								h[q2bin,wbin,'EXP','C'][ivar].Draw("sames")
								h[q2bin,wbin,'EXP','H'][ivar].Multiply(hF)
								h[q2bin,wbin,'EXP','H'][ivar].Draw("sames")
								h[q2bin,wbin,'SIM','T'][ivar].DrawNormalized("sames",1000)
								pad.Update()
				if view=="full_ana" and q2wbintitle.has_key((q2bin,wbin)):
					outdir_q2bin=os.path.join(outdir,str(q2bin))
					if not os.path.exists(outdir_q2bin):
						os.makedirs(outdir_q2bin)
					c.SaveAs("%s/c%s_%s.png"%(outdir_q2bin,wbin,q2bin))
					c.Close()
			if view=="q2_evltn" and wbintitle.has_key((wbin)):
				c.SaveAs("%s/c%s.png"%(outdir,wbin))
				c.Close()
		return

	def plot_obs_R2(self,h5l,dtyp='EXP',seq='C'):
		"""
		+ Plot Obs_R2
		"""
		#! Define some objects needed by this method
		R2l=['T+L','LT','TT','LTp']

		# outdir=os.path.join(self.OUTDIR_OBS_R2)
		# if not os.path.exists(outdir):
		# 	os.makedirs(outdir)

		#! Get all q2- and w-bins from the keys of hVSTX
		q2bins_le,wbins_le=[],[]
		for k in h5l.keys():
			q2bins_le.append(k[0])
			wbins_le.append(k[1])
		#! Keep only unique entries
		q2bins_le=set(q2bins_le) #! Keep only unique entries
		q2bins_le=list(q2bins_le) #! convert 'set' back to 'list'
		q2bins_le=sorted(q2bins_le)
		wbins_le=set(wbins_le) #! Keep only unique entries
		wbins_le=list(wbins_le) #! convert 'set' back to 'list'
		wbins_le=sorted(wbins_le)
		print "going to plot R2 obs for:"
		print "Q2:"
		print q2bins_le
		print "W:"
		print wbins_le

		#! Extract hR2(EC,vst,var) in every q2w-bin
		for wbin in wbins_le:
			for q2bin in q2bins_le:
				for vst in self.VSTS:
					for seq in ['C']:
						for EC in R2l:
							if EC=='LTp':
								if h5l.has_key((q2bin,wbin,dtyp,vst,seq,'POS')) and h5l.has_key((q2bin,wbin,dtyp,vst,seq,'NEG')):
									h5_POS=h5l[q2bin,wbin,dtyp,vst,seq,'POS']
									h5_NEG=h5l[q2bin,wbin,dtyp,vst,seq,'NEG']

									h5_POSm=thntool.MultiplyBy(h5_POS,'sphi',1)
									h5_NEGm=thntool.MultiplyBy(h5_NEG,'sphi',-1)

									for var in VARS:
										if var=='PHI': continue
										if vst==1 and var=='M2': continue
										if vst==2 and var=='M1': continue
										if vst==3 and var=='M1': continue

										hR2_POS=h5_POSm.Projection(H5_DIM[var],"E")
										hR2_POS.SetName("D_%d_%s_%s_POS"%(vst,var,seq))
										hR2_POS.SetTitle("R2_{LT^{\'}}_VST%d_%s_%s_POS(%.2f,%.3f)"%(vst,var,seq,q2bin,wbin))
										hR2_POS.Scale(1/math.pi)
										hR2_POS.Scale(1/50000)
										hR2_NEG=h5_NEGm.Projection(H5_DIM[var],"E")
										hR2_NEG.SetName("D_%d_%s_%s_NEG"%(vst,var,seq))
										hR2_NEG.SetTitle("R2_{LT^{\'}}_VST%d_%s_%s_NEG(%.2f,%.3f)"%(vst,var,seq,q2bin,wbin))
										hR2_NEG.Scale(1/math.pi)
										hR2_NEG.Scale(1/50000)
										hR2_AVG=hR2_POS.Clone("avg")
										hR2_AVG.SetName("D_%d_%s_%s_AVG"%(vst,var,seq))
										hR2_AVG.SetTitle("R2_{LT^{\'}}_VST%d_%s_%s_AVG(%.2f,%.3f)"%(vst,var,seq,q2bin,wbin))
										hR2_AVG.Add(hR2_NEG,1)
										hR2_AVG.Scale(0.5)

										#ROOT.gStyle.SetOptStat("nemruo")

										outdir=os.path.join(self.OUTDIR_OBS_R2,"D_VST%d_%s_%s"%(vst,var,seq))
										if not os.path.exists(outdir):
											os.makedirs(outdir)
										c=ROOT.TCanvas()
										hR2_AVG.Draw()
										c.SaveAs("%s/cD_%.3f_%0.2f_VST%d_%s.png"%(outdir,wbin,q2bin,vst,var))
										c.Close()
							else: #! Currently, not implemented if EC neq LTp
								continue




		
	def disp_1D(self,view="q2_evltn",dtypl=['EXP','SIM'],seql=['T','C','H','F']):
		"""
		Walk the ROOT file and extract:
			+ hVST1(q2bin,wbin,dtyp,seq)=[VST1.M1, VST1.THETA, VST1.ALPHA]
			+ hVST2(q2bin,wbin,dtyp,seq)=[VST2.M2, VST2.THETA, VST2.ALPHA]
			+ hVST3(q2bin,wbin,dtyp,seq)=[VST3.M2, VST3.THETA, VST3.ALPHA]
		where (as per DataAna::makeYields())
			+ VST1=delta_pp,pim
			+ VST2=rho,p
			+ VST3=delta_0,pip
		"""
		self.OUTDIR_OBS_1D=os.path.join(self.OUTDIR,"Obs_1D")
		if not os.path.exists(self.OUTDIR_OBS_1D):
			os.makedirs(self.OUTDIR_OBS_1D)

		print "In DispYields::disp_1D()"
		#! 1. First get all q2wbin directories from file
		q2wbinl=self.get_q2wbinlist()
		print q2wbinl

		#! 2. Get q2bng and wbng from file
		q2bng=self.get_q2bng()
		wbng=self.get_wbng()
		print "1D Observables will be extracted in the following Q2-bin,W-bin 2D space:"
		print "Q2 bins:"
		print q2bng['BINS']
		print "W bins:"
		print wbng['BINS']

		#! 2. Now get relevant histograms
		hVST1,hVST2,hVST3={},{},{}
		q2wbinl_bad={} 
		for q2wbin in q2wbinl:
			print "Getting hVST1,hVST2,hVST3 for",q2wbin
			#! First make sure this bin is "good"
			badbin=False
			reason=''
			for vst in self.VSTS:
				ht_R=self.FEXP.Get("%s/VST%d/R/h_M1"%(q2wbin,vst))
				ht_C=self.FEXP.Get("%s/VST%d/C/h_M1"%(q2wbin,vst))
				if ht_R.Integral()==0:
					reason+="ER=0 for VST%d;"%vst
					badbin=True
				elif ht_C.Integral()==0:
					reason+="SA=0 in ER bins for VST%d;"%vst
					badbin=True
			if badbin:
				q2wbinl_bad[q2wbin]=reason
				print "\"Bad\" bin: %s"%q2wbinl_bad[q2wbin]
				continue
			#! Now that bin is good, get hists
			q2bin_le=self.get_q2bin_le(q2wbin)
			wbin_le=self.get_wbin_le(q2wbin)
			#h_dpp,h_rho,h_dzr=OrderedDict(),OrderedDict(),OrderedDict()
			for dtyp in dtypl:
				for seq in seql:
					if dtyp=='EXP' and seq=='T': continue
					if dtyp=='SIM' and (seq=='H' or seq=='C'): continue

					if dtyp=='EXP':
						f=self.FEXP
					if dtyp=='SIM':
						f=self.FSIM
					
					hVST1[q2bin_le,wbin_le,dtyp,seq]=[f.Get("%s/VST1/%s/h_M1"%(q2wbin,seq)),f.Get("%s/VST1/%s/h_THETA"%(q2wbin,seq)),
								 	f.Get("%s/VST1/%s/h_ALPHA"%(q2wbin,seq))]
					hVST2[q2bin_le,wbin_le,dtyp,seq]=[f.Get("%s/VST2/%s/h_M2"%(q2wbin,seq)),f.Get("%s/VST2/%s/h_THETA"%(q2wbin,seq)),
								 	f.Get("%s/VST2/%s/h_ALPHA"%(q2wbin,seq))]
					hVST3[q2bin_le,wbin_le,dtyp,seq]=[f.Get("%s/VST3/%s/h_M2"%(q2wbin,seq)),f.Get("%s/VST3/%s/h_THETA"%(q2wbin,seq)),
								 	f.Get("%s/VST3/%s/h_ALPHA"%(q2wbin,seq))]
		#print "keys in hVST1:"
		#print hVST1.keys()
		#! Write out the list of "bad" q2w bins and the reason why it is bad
		fout=open("%s/bad-q2wbins.txt"%self.OUTDIR_OBS_1D,"w")
		fout.write("Following are the \"bad\" q2w bins:\n")
		for k in q2wbinl_bad:
			fout.write("%s:%s\n"%(k,q2wbinl_bad[k]))
		fout.close()
		print "Finished getting hVST1,hVST2,hVST3. Now going to display yields"
		if view=="q2_evltn":
			self.plot_obs_1D(hVST1,hVST2,hVST3,view="q2_evltn",dtyp='EXP',seq='C')
			self.plot_obs_1D(hVST1,hVST2,hVST3,view="q2_evltn",dtyp='EXP',seq='F')
			self.plot_obs_1D(hVST1,hVST2,hVST3,view="q2_evltn",dtyp='SIM',seq='F')
		elif view=="full_ana":
			self.plot_obs_1D(hVST1,hVST2,hVST3,view="full_ana")
		else:
			sys.exit("view=%s not recognized. Exiting."%view)
		print "Done DispYields::disp_1D()"
		print "If the progam is not terminating, then Python is probably doing \"garbage collection\"(?); Wait a while!"
		return

	def disp_R2(self,dtypl=['EXP'],seql=['C']):
		"""
		Walk the ROOT file and extract:
			+ h5(q2bin,wbin,dtyp,vst,seq,hel)
		"""
		self.OUTDIR_OBS_R2=os.path.join(self.OUTDIR,"Obs_R2")
		if not os.path.exists(self.OUTDIR_OBS_R2):
			os.makedirs(self.OUTDIR_OBS_R2)

		print "In DispYields::disp_R2()"
		#! 1. First get all q2wbin directories from file
		q2wbinl=self.get_q2wbinlist()
		print q2wbinl

		#! 2. Now get relevant histograms
		h5={}
		q2wbinl_bad={} 
		for q2wbin in q2wbinl:
			print "Getting h5 for",q2wbin
			#! First make sure this bin is "good"
			badbin=False
			reason=''
			for vst in self.VSTS:
				h5_UNPOL=self.FEXP_HEL.Get("%s/VST%d/R/h5_UNPOL"%(q2wbin,vst))
				if thntool.GetIntegral(h5_UNPOL)==0:
					reason+="ER=0 for VST%d;"%vst
					badbin=True
			if badbin:
				q2wbinl_bad[q2wbin]=reason
				print "\"Bad\" bin: %s"%q2wbinl_bad[q2wbin]
				continue
			#! Now that bin is "good", get h5
			q2bin_le=self.get_q2bin_le(q2wbin)
			wbin_le=self.get_wbin_le(q2wbin)
			#h_dpp,h_rho,h_dzr=OrderedDict(),OrderedDict(),OrderedDict()
			for dtyp in dtypl:
				if   dtyp=='EXP':f=self.FEXP_HEL
				elif dtyp=='SIM':f=self.FSIM
				for vst in self.VSTS:
					for seq in seql:
						h5[q2bin_le,wbin_le,dtyp,vst,seq,'UNPOL']=f.Get("%s/VST%d/%s/h5_UNPOL"%(q2wbin,vst,seq))
						h5[q2bin_le,wbin_le,dtyp,vst,seq,'POS']=f.Get("%s/VST%d/%s/h5_POS"%(q2wbin,vst,seq))
						h5[q2bin_le,wbin_le,dtyp,vst,seq,'NEG']=f.Get("%s/VST%d/%s/h5_NEG"%(q2wbin,vst,seq))

		print "keys in h5:"
		print h5.keys()

		#! Write out the list of "bad" q2w bins and the reason why it is bad
		fout=open("%s/bad-q2wbins.txt"%self.OUTDIR_OBS_R2,"w")
		fout.write("Following are the \"bad\" q2w bins:\n")
		for k in q2wbinl_bad:
			fout.write("%s:%s\n"%(k,q2wbinl_bad[k]))
		fout.close()
		print "Finished getting h5s. Now going to display R2s"
		self.plot_obs_R2(h5)
		print "Done DispYields::disp_R2()"
		print "If the progam is not terminating, then Python is probably doing \"garbage collection\"(?); Wait a while!"
		return

	def disp_integ_yield(self,seq='F'):
		"""
		Walk the ROOT file and plot y(w;q2). 
		"""
		outdir=os.path.join(self.OUTDIR,"Obs_IntgYld")
		if not os.path.exists(outdir):
			os.makedirs(outdir)

		#! 1. First of all get q2bng information
		q2bng=self.get_q2bng()
		#print q2bng['BINW']
		print "Yield(w) will be plotted for the following Q2 bins:"
		print q2bng['BINS']

		#! 2. Get all q2wbins
		q2wbinl=self.get_q2wbinlist()
		#print q2wbinl

		#! 3. Now put together y[q2bin]=(w,yield)
		y=[{} for i in range(q2bng['NBINS'])]
		for q2wbin in q2wbinl:
			#! Get iq2bin
			q2bin_le=float(q2wbin.split('_')[0].split('-')[0])
			iq2bin=q2bng['BINS_LE'].index(q2bin_le)
			#! Get w, yield
			w=float(q2wbin.split('_')[1].split('-')[0])
			h5_UNPOL=self.FEXP.Get("%s/VST1/%s/h5_UNPOL"%(q2wbin,seq))
			y[iq2bin][w]=thntool.GetIntegral(h5_UNPOL)
		#! Make sure y[q2bin]=(w,yield) are sorted by w
		oy=[{} for i in range(len(y))]
		for iq2bin in range(len(y)):
			oy[iq2bin]=OrderedDict(sorted(y[iq2bin].items()))
			# print iq2bin
			# print y[iq2bin]
			# print oy[iq2bin] 

		#! 4. Now plot
		fig=plt.figure()
		ax=plt.subplot(111)
		clrs=['red','green','cyan','blue','black','yellow','brown','orange']
		for iq2bin in range(len(oy)):
			lbl='[%.2f-%.2f]'%(q2bng['BINS_LE'][iq2bin],q2bng['BINS_UE'][iq2bin])
			clr=clrs[iq2bin]
			ax.scatter(oy[iq2bin].keys(),oy[iq2bin].values(),label=lbl,c=clr)
		ax.set_ylim(0,600000)
		ax.legend()
		fig.savefig('%s/integ_yield.png'%(outdir))	
		#fig.savefig('test.png')		

	def get_sim_stats(self):
		"""
		Walk the ROOT file and obtain simstats(ss) for a h5_UNPOL in a Q2-W bin:
		ss={'T':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]],
		    'R':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]],
		    'A':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]],
		    'H':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]]}

		where:
			+ nbins=number of filled bins in a h5_UNPOL
			+ N=sum({n_i}) were n_i=events per bin (number of events in a h5_UNPOL)
			+ mu=average({n_i} (average of number of events per bin)
			+ sg=(RMS({n_i}) (RMS of number of events per bin)
		"""
		#! First get all q2wbin directories from file
		q2ws=self.get_q2wbinlist()
		#print "Processing sim_stats for %s:"%self.Q2W
		print q2ws

		ss={'T':[],'R':[],'A':[],'H':[]}
		f=ROOT.TFile(self.FSIM.GetName())
		for seq in ['T','R','A','H']:
			for q2w in q2ws:
				print "Processing sim_stats for %s"%q2w
				#! Determine q2,w
				q2bin=q2w.split('_')[0]
				wbin=q2w.split('_')[1]
				#print q2bin,wbin
				q2=float(q2bin.split('-')[0])
				w=float(wbin.split('-')[0])
				#print q2,w
				#! Determine nbins,N,mu,sg for this q2,w
				h5_UNPOL=f.Get("%s/VST1/%s/h5_UNPOL"%(q2w,seq))
				nbins=thntool.GetNbinsNotEq0(h5_UNPOL)
				N=thntool.GetIntegral(h5_UNPOL)
				binc_stats=np.zeros(2,'f')
				thntool.GetBinContentDistStats(h5_UNPOL,binc_stats)
				mu=binc_stats[0]
				sg=binc_stats[1]
				ss[seq].append([q2,w,nbins,N,mu,sg])
			# #! Compute average
			# ss[seq].append(nevts/len(q2ws))
			# ss[seq].append(nbins/len(q2ws))
		return ss

	def get_sim_stats_commonbins(self):
		"""
		The same as get_sim_stats(), but with the difference that mu,sg
		are obtained only for R-PS bins. I am currently not using this since
		this process takes impractically long.

		Walk the ROOT file and obtain simstats(ss) for a h5_UNPOL in a Q2-W bin:
		ss={'T':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]],
		    'R':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]],
		    'A':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]],
		    'H':[[q21,w1,nbins,N,mu,sg],...,[q2N,wN,nbins,N,mu,sg]]}

		where:
			+ nbins=number of filled bins in a h5_UNPOL
			+ N=sum({n_i}) were n_i=events per bin (number of events in a h5_UNPOL)
			+ mu=average({n_i}) (average of number of events per bin)
				+ Note, average computed over only R-PS bins
			+ sg=(RMS({n_i}) (RMS of number of events per bin)
				+ Note, average computed over only R-PS bins
		"""
		#! First get all q2wbin directories from file
		q2ws=self.get_q2wbinlist()
		print "Processing sim_stats for %s:"%self.Q2W
		print q2ws

		ss={'T':[],'R':[],'A':[],'H':[]}
		f=ROOT.TFile(self.FSIM.GetName())
		for q2w in q2ws:
			print "Processing %s..."%q2w
			#! Determine q2,w
			q2bin=q2w.split('_')[0]
			wbin=q2w.split('_')[1]
			#print q2bin,wbin
			q2=float(q2bin.split('-')[0])
			w=float(wbin.split('-')[0])
			#print q2,w
			#! First get all h5_UNPOLs
			h5_UNPOL={}
			for seq in ['T','R','A','H']:
				h5_UNPOL[seq]=f.Get("%s/VST1/%s/h5_UNPOL"%(q2w,seq))
			#! Now get simstats	
			for seq in ['T','R','A','H']:
				#! Determine nbins,N,mu,sg for this q2,w
				nbins=thntool.GetNbinsNotEq0(h5_UNPOL[seq])
				N=thntool.GetIntegral(h5_UNPOL[seq])
				binc_stats=np.zeros(2,'f')
				print "here"
				thntool.GetBinContentDistStatsCommonBins(h5_UNPOL[seq],h5_UNPOL['R'],binc_stats)
				print "here1"
				mu=binc_stats[0]
				sg=binc_stats[1]
				ss[seq].append([q2,w,nbins,N,mu,sg])
		return ss


	def get_q2wbinlist(self):
		"""
		The ROOT file is arranged in a Tree Structure. The
		Observable histograms (obs-hists) are located as files in the following directory-path(dirpath):
		q2wbin/vst/seq/hists
		"""
		q2ws=[]
		i=0 #! for getting limted q2wbins
		for path,dirs,files in self.FEXP.walk():
			if path=="":continue #! Avoid root path
			path_arr=path.split("/")
			if len(path_arr)==1:
				q2ws.append(path)
				i+=1
			#if i>10: break #! Uncomment/comment -> Get limited q2w-bins/Get all q2w-bins
		return q2ws

	def get_q2bng(self):
		return self.get_xbng(x='Q2')
	def get_wbng(self):
		return self.get_xbng(x='W')

	def get_xbng(self,x):
		"""
		Gets x=(Q2 or W) binning information from yield_exp.root
		"""

		#! First get all q2w-bins
		q2wbinl=self.get_q2wbinlist()
		
		#! 2. From q2wbinl, identify xbng 
		xbins=[]
		for q2wbin in q2wbinl:
			#print q2wbin
			if(x=='Q2'):
				xbins.append(float(q2wbin.split('_')[0].split('-')[0]))
				xbins.append(float(q2wbin.split('_')[0].split('-')[1]))
			elif(x=='W'):
				xbins.append(float(q2wbin.split('_')[1].split('-')[0]))
				xbins.append(float(q2wbin.split('_')[1].split('-')[1]))
			else:
				sys.exit("DispYields::get_xbng() x=%s not recognized"%x)
		xbins=set(xbins) #! Keep only unique entries
		xbins=list(xbins) #! convert 'set' back to 'list'
		xbins=sorted(xbins) #! Order entries
		#print xbins
		xbins_le=xbins[:-1]
		xbins_ue=xbins[1:]
		#print xbins_le
		#print xbins_ue
		xmin=min(xbins)
		xmax=max(xbins)
		xbinw=xbins[1]-xbins[0]
		nxbins=len(xbins_le)
		xbng={'MIN':xmin,'MAX':xmax,'BINW':xbinw,'NBINS':nxbins,
		   'BINS_LE':xbins_le,'BINS_UE':xbins_ue,'BINS':xbins}
		return xbng

	def get_q2bin_le(self,q2wbin):
		return float(q2wbin.split('_')[0].split('-')[0])
	def get_wbin_le(self,q2wbin):
		return float(q2wbin.split('_')[1].split('-')[0])

	def label_hist_obs1D(self,hVST1,hVST2,hVST3):
		q2wbintitle={}
		wbintitle={}
		for vst in self.VSTS:
			if   vst==1: 
				hl=hVST1
				m='M1'
			elif vst==2: 
				hl=hVST2
				m='M2'
			elif vst==3: 
				hl=hVST3
				m='M2'

			for k in hl.keys():
				#! extract 'q2wbintitle' that is already a part of title set in proc_yields.py:
				#! 	+ title=[q2min-q2max]_[wmin-wmax]_VSTX_SEQ_VAR
				title_orig=hl[k][0].GetTitle().split("_")
				#q2wbintitle="%s_%s"%(title_orig[0],title_orig[1])
				if not q2wbintitle.has_key((k[0],k[1])):
					q2wbintitle[k[0],k[1]]="%s_%s"%(title_orig[0],title_orig[1])
				if not wbintitle.has_key(k[1]):
					wbintitle[k[1]]="%s"%(title_orig[1])

			for k in hl.keys():	
				#! Now set title
				# hl[k][0].SetTitle("%s:%s"%(self.VAR_NAMES[(vst,m)],q2wbintitle))
				# hl[k][1].SetTitle("%s:%s"%(self.VAR_NAMES[(vst,'THETA')],q2wbintitle))
				# hl[k][2].SetTitle("%s:%s"%(self.VAR_NAMES[(vst,'ALPHA')],q2wbintitle))

				hl[k][0].SetTitle("%s"%(self.VAR_NAMES[(vst,m)]))
				hl[k][1].SetTitle("%s"%(self.VAR_NAMES[(vst,'THETA')]))
				hl[k][2].SetTitle("%s"%(self.VAR_NAMES[(vst,'ALPHA')]))
		return q2wbintitle,wbintitle

	# def set_q2wbintitle(self,hl):
	# 	q2wbintitle={}
	# 	for k in hl.keys():
	# 		if k[0]==q2bin and k[1]==wbin:
	# 			title_orig=hl[k][0].GetTitle().split("_")
	# 			q2wbintitle[]"%s_%s"%(title_orig[0],title_orig[1])
	# 			break
	# 	return q2wbintitle



