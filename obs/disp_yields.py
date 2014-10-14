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

#! For Luminosity & vgflux
LUM=19.844 #fb^-1
LUM_INVFB_TO_INVMICROB=1000000000

PI=3.14159265358979312
E1F_E0=5.499
FSC=0.00729735253
A=FSC
NA=6.02214129E23
QE=1.60217646E-19
MP=0.93827203

def nu(w,q2):
	return (w*w-MP*MP+q2)/(2*MP)

def epsilon(w,q2):
	n=nu(w,q2)
	e0=E1F_E0
	e1=e0-n
	epsInv=1+2*(q2+n*n)/(4*e0*e1-q2)
	return 1.0/epsInv

def getvgflux(w,q2,e0=E1F_E0):
	eps=epsilon(w,q2)
	return A*w*(w*w-MP*MP)/(4*PI*e0*e0*MP*MP*q2*(1-eps))

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

		#! Russian-normalize theta distributions
		print("Russian-normalizing theta distributions ...")
		for hd in [hVST1,hVST2,hVST3]:
			for k in hd:
				self.russ_norm_theta_dist(hd[k][1])

		#! Set up some plotting related styles and aesthetics 
		if view=="q2_evltn":
			colors=["kRed","kOrange","kYellow","kGreen+3","kGreen","kCyan","kBlue","kMagenta"]
			coll=[]
			for iq2bin in range(len(q2bins_le)):
				coll.append(ROOT.gROOT.ProcessLine(colors[iq2bin]))
		elif view=="full_ana":
			coll={('EXP','R'):ROOT.gROOT.ProcessLine("kYellow"),
			      ('EXP','C'):ROOT.gROOT.ProcessLine("kCyan"),
			      ('EXP','F'):ROOT.gROOT.ProcessLine("kBlue"),
			      ('EXP','H'):ROOT.gROOT.ProcessLine("kBlack"),
			      ('SIM','T'):ROOT.gROOT.ProcessLine("kGreen"),
			      ('SIM','R'):ROOT.gROOT.ProcessLine("kMagenta"),
			      ('SIM','F'):ROOT.gROOT.ProcessLine("kRed")}
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
		#! 3x3 TCanvas-pads as per VSTs: VST1=1,2,3; VST2=3,6,9; VST3=2,5,8 (To adapt to Gleb's display)
		pads=[[1,4,7],[3,6,9],[2,5,8]]
		map_padnums_vars=[zip(pads[0],ivars),zip(pads[1],ivars),zip(pads[2],ivars)]
		for wbin in wbins_le:
			if view=="q2_evltn" and wbintitle.has_key((wbin)):
				c=ROOT.TCanvas("c","c",1000,1000)
				pad_t=ROOT.TPad("pad_t","Title pad",0.05,0.97,0.95,1.00)
				#pad_t.SetFillColor(11)
  				pad_p=ROOT.TPad("pad_p","Plots pad",0.01,0.01,0.99,0.95);
  				pad_p.SetFillColor(11)
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
  					pad_p.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
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
								h[q2bin,wbin,'EXP','R'][ivar].SetMarkerColor(coll['EXP','R'])
								h[q2bin,wbin,'EXP','C'][ivar].SetMarkerColor(coll['EXP','C'])
								h[q2bin,wbin,'EXP','F'][ivar].SetMarkerColor(coll['EXP','F'])
								h[q2bin,wbin,'EXP','H'][ivar].SetMarkerColor(coll['EXP','H'])
								h[q2bin,wbin,'SIM','T'][ivar].SetMarkerColor(coll['SIM','T'])
								h[q2bin,wbin,'SIM','R'][ivar].SetMarkerColor(coll['SIM','R'])
								h[q2bin,wbin,'SIM','F'][ivar].SetMarkerColor(coll['SIM','F'])
								h[q2bin,wbin,'EXP','R'][ivar].SetLineColor(coll['EXP','R'])
								h[q2bin,wbin,'EXP','C'][ivar].SetLineColor(coll['EXP','C'])
								h[q2bin,wbin,'EXP','F'][ivar].SetLineColor(coll['EXP','F'])
								h[q2bin,wbin,'EXP','H'][ivar].SetLineColor(coll['EXP','H'])
								h[q2bin,wbin,'SIM','T'][ivar].SetLineColor(coll['SIM','T'])
								h[q2bin,wbin,'SIM','R'][ivar].SetLineColor(coll['SIM','R'])
								h[q2bin,wbin,'SIM','F'][ivar].SetLineColor(coll['SIM','F'])
								h[q2bin,wbin,'EXP','R'][ivar].SetMarkerStyle(ROOT.gROOT.ProcessLine("kOpenStar"))
								h[q2bin,wbin,'SIM','R'][ivar].SetMarkerStyle(ROOT.gROOT.ProcessLine("kOpenStar"))
								h[q2bin,wbin,'EXP','H'][ivar].SetMarkerStyle(ROOT.gROOT.ProcessLine("kCircle"))
								h[q2bin,wbin,'SIM','T'][ivar].SetMarkerStyle(ROOT.gROOT.ProcessLine("kPlus"))
								
								#! Draw noramlized histograms (since I want to compare their distributions)
								hEF_n=h[q2bin,wbin,'EXP','F'][ivar].DrawNormalized("",1000)
								hSF_n=h[q2bin,wbin,'SIM','F'][ivar].DrawNormalized("sames",1000)
								hER_n=h[q2bin,wbin,'EXP','R'][ivar].DrawNormalized("sames",1000)
								hSR_n=h[q2bin,wbin,'SIM','R'][ivar].DrawNormalized("sames",1000)
								hST_n=h[q2bin,wbin,'SIM','T'][ivar].DrawNormalized("sames",1000)
								#! Set the minimum and maximum of y coordinate of histograms
								maxl=[hEF_n.GetMaximum(),hSF_n.GetMaximum(),hER_n.GetMaximum(),hSR_n.GetMaximum(),hST_n.GetMaximum()]
								maximum=max(maxl)
								for htmp in [hEF_n,hSF_n,hER_n,hSR_n,hST_n]:
									htmp.SetMinimum(0.)
									htmp.SetMaximum(maximum+10)
								#! Scale EC, EH as per hEF_n
								hF=hEF_n.Clone("hF")
								hF.Divide(h[q2bin,wbin,'EXP','F'][ivar])
								h[q2bin,wbin,'EXP','C'][ivar].Multiply(hF)
								h[q2bin,wbin,'EXP','C'][ivar].Draw("sames")
								h[q2bin,wbin,'EXP','H'][ivar].Multiply(hF)
								h[q2bin,wbin,'EXP','H'][ivar].Draw("sames")
								#! Upate pad
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

	def plot_obs_R2(self,h5l,R2,hel,dtypl,seql):
		"""
		+ Plot Obs_R2
		"""
		#! Define some objects needed by this method
		R2_named={'A':'R2_{T+L}','B':'R2_{LT}','C':'R2_{TT}','D':'R2_{LT^{\'}}'}
		mfuncd={'A':'1','B':'cphi','C':'c2phi','D':'sphi'}

		#! Stats Box
		# ROOT.gStyle.SetOptStat(0)

		# # ROOT.gStyle.SetLabelSize(0.5,"t")
		# # ROOT.gStyle.SetTitleSize(0.5,"t")
		# #ROOT.gStyle.SetPaperSize(20,26);
		# ROOT.gStyle.SetPadTopMargin(0.15)#(0.05);
		# #ROOT.gStyle.SetPadRightMargin(0.09)#(0.05);
		# #ROOT.gStyle.SetPadBottomMargin(0.20)#(0.16);
		# #ROOT.gStyle.SetPadLeftMargin(0.15)#(0.12);

		# ROOT.gStyle.SetTitleW(10)# //title width 
		# ROOT.gStyle.SetTitleFontSize(20)# //title width 
		# ROOT.gStyle.SetTitleH(0.15)# //title height 
		# ROOT.gStyle.SetTitleY(1)# //title Y location 

		#!get rid of X error bars and y error bar caps
		ROOT.gStyle.SetErrorX(0.001);

		# outdir=os.path.join(self.OUTDIR_OBS_R2)
		# if not os.path.exists(outdir):
		# 	os.makedirs(outdir)

		#! Get all q2- and w-bins from the keys of hVSTX
		q2bins_lel,wbins_lel=[],[]
		for k in h5l.keys():
			q2bins_lel.append(k[0])
			wbins_lel.append(k[1])
		#! Keep only unique entries
		q2bins_lel=set(q2bins_lel) #! Keep only unique entries
		q2bins_lel=list(q2bins_lel) #! convert 'set' back to 'list'
		q2bins_lel=sorted(q2bins_lel)
		wbins_lel=set(wbins_lel) #! Keep only unique entries
		wbins_lel=list(wbins_lel) #! convert 'set' back to 'list'
		wbins_lel=sorted(wbins_lel)
		print "going to plot R2 obs for:"
		print "Q2:"
		print q2bins_lel
		print "W:"
		print wbins_lel

		#! Extract hR2(q2bin_le,wbin_le,vst,var,hel,dtyp,seq)
		hR2={}
		for wbin_le in wbins_lel:
			for q2bin_le in q2bins_lel:
				for vst in self.VSTS:
						for dtyp in dtypl:
							for seq in seql:
								if h5l.has_key((q2bin_le,wbin_le,dtyp,vst,seq,hel)):
									h5=h5l[q2bin_le,wbin_le,dtyp,vst,seq,hel]
									if   hel=='UNPOL': h5m=thntool.MultiplyBy(h5,mfuncd[R2],1)	
									elif hel=='POS':   h5m=thntool.MultiplyBy(h5,mfuncd[R2],1)	
									elif hel=='NEG':   h5m=thntool.MultiplyBy(h5,mfuncd[R2],-1)	
									else: print hel,"is not recognized. hel=POS/NEG/UNPOL"

									for var in VARS:
										if var=='PHI': continue
										if vst==1 and var=='M2': continue
										if vst==2 and var=='M1': continue
										if vst==3 and var=='M1': continue

										hR2[q2bin_le,wbin_le,hel,vst,var,dtyp,seq]=h5m.Projection(H5_DIM[var],"E")
										hR2[q2bin_le,wbin_le,hel,vst,var,dtyp,seq].SetName("%s_VST%d_%s_%s_%s"%(R2,vst,var,seq,hel))
										hR2[q2bin_le,wbin_le,hel,vst,var,dtyp,seq].SetTitle("%s(%s:%s:%.2f,%.3f)"%(R2_named[R2],self.VAR_NAMES[(vst,var)],hel,q2bin_le,wbin_le))
										hR2[q2bin_le,wbin_le,hel,vst,var,dtyp,seq].Scale(1/math.pi)
										hR2[q2bin_le,wbin_le,hel,vst,var,dtyp,seq].Scale(1/50000)
				print "keys in hR2:\n",hR2.keys()

				#! Now display hR2
				for vst in self.VSTS:
					for var in VARS:
						if var=='PHI': continue
						if vst==1 and var=='M2': continue
						if vst==2 and var=='M1': continue
						if vst==3 and var=='M1': continue

						outdir=os.path.join(self.OUTDIR_OBS_R2,"VST%d_%s"%(vst,var))
						if not os.path.exists(outdir):
								os.makedirs(outdir)
						for wbin_le in wbins_lel:
							for q2bin_le in q2bins_lel:
								if len(dtypl)==2:
									c=ROOT.TCanvas("","",800,1000)
									c.Divide(1,2)
								else:
									c=ROOT.TCanvas()
								# pad_exp=ROOT.TPad("pad_exp","",0,0,1,1)
								# pad_sim=ROOT.TPad("pad_sim","",0,0,1,1)
								#pad_sim.SetFillStyle(4000)# will be transparent
								#pad_exp.Draw()
								#pad_sim.Draw()
								l=ROOT.TLegend(0.1,0.8,0.2,0.9)
								mrkrd={('EXP','C'):ROOT.gROOT.ProcessLine("kCyan"),('EXP','F'):ROOT.gROOT.ProcessLine("kBlue"),('SIM','F'):ROOT.gROOT.ProcessLine("kRed")}
								i=0
								for dtyp in dtypl:
									for seq in seql:
										if dtyp=='SIM' and seq=='C': continue
										if hR2.has_key((q2bin_le,wbin_le,hel,vst,var,dtyp,seq)):
											h=hR2[q2bin_le,wbin_le,hel,vst,var,dtyp,seq]
											h.SetMarkerStyle(ROOT.gROOT.ProcessLine("kFullCircle"))
											h.SetMarkerColor(mrkrd[(dtyp,seq)])
											l.AddEntry(h,"%s-%s"%(dtyp,seq),"p")
											if dtyp=='EXP': 
												c.cd(1)
												if i==0:
													h.Draw()
													i+=1;
												else:
													h.Draw("sames")
												# pad_exp.Modified()
												# c.cd()
											if dtyp=='SIM':
												c.cd(2)
												h.Draw()
												# pad_sim.cd()
												# ymin=h.GetMinimum()
												# ymax=h.GetMaximum()
												# dy=(ymax-ymin)/0.8#10 per cent margins top and bottom
												# xmin=h.GetXaxis().GetXmin()
												# xmax=h.GetXaxis().GetXmax()
												# dx=(xmax-xmin)/0.8#10 per cent margins top and bottom
												# pad_sim.Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy)
												# pad_sim.Draw()
												# pad_sim.cd()
												# h.Draw("][sames")
												# pad_sim.Update()
												# ax=ROOT.TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,50510,"+L")
												# ax.SetLabelColor(ROOT.gROOT.ProcessLine("kRed"))
												# ax.Draw()
								c.cd(1)				
								l.Draw("same")
								c.SaveAs("%s/c_w%.3f_q%0.2f.png"%(outdir,wbin_le,q2bin_le))
								c.Close()
						


				# outdir=os.path.join(self.OUTDIR_OBS_R2,"w%s_q%s"%(wbin,q2bin))
				# if not os.path.exists(outdir):
				# 	os.makedirs(outdir)
				# for vst in self.VSTS:
				# 	for var in VARS:
				# 			if var=='PHI': continue
				# 			if vst==1 and var=='M2': continue
				# 			if vst==2 and var=='M1': continue
				# 			if vst==3 and var=='M1': continue
				# 			c=ROOT.TCanvas()
				# 			i=0
				# 			for dtyp in dtypl:
				# 				for seq in seql:
				# 					if i==0:
				# 						hR2[vst,var,dtyp,seq].Draw()
				# 					else:
				# 						hR2[vst,var,dtyp,seq].Draw("sames")
				# 			c.SaveAs("%s/c%s_%.3f_%0.2f_VST%d_%s.png"%(outdir,R2,wbin,q2bin,vst,var))
				# 			c.Close()




		
	def disp_1D(self,view="q2_evltn",dtypl=['EXP','SIM'],seql=['T','R','C','H','F']):
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
		q2bng=self.get_q2bng(q2wbinl)
		wbng=self.get_wbng(q2wbinl)
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

	def disp_R2(self,R2='D',q2min=1.25,q2max=1.75,hel='POS',dtypl=['EXP'],seql=['C','F']):
		"""
		Walk the ROOT file and extract:
			+ h5(q2bin,wbin,dtyp,vst,seq,hel)
		"""
		self.OUTDIR_OBS_R2=os.path.join(self.OUTDIR,"Obs_R2",R2,hel,"Q2_%.2f-%.2f"%(q2min,q2max))
		if not os.path.exists(self.OUTDIR_OBS_R2):
			os.makedirs(self.OUTDIR_OBS_R2)

		print "In DispYields::disp_R2()"
		#! 1. First get all q2wbin directories from file
		q2wbinl=self.get_q2wbinlist(q2min=q2min,q2max=q2max,dbg=True,dbg_bins=10)
		print q2wbinl

		#! 2. Now get relevant histograms
		h5={}
		q2wbinl_bad={} 
		for q2wbin in q2wbinl:
			print "Getting h5 for",q2wbin
			#! First make sure this bin is "good"
			is_badbin=self.is_bad_q2wbin(q2wbin,q2wbinl_bad)
			if is_badbin: continue
			#! Now that bin is "good", get h5
			q2bin_le=self.get_q2bin_le(q2wbin)
			wbin_le=self.get_wbin_le(q2wbin)
			for dtyp in dtypl:
				if   dtyp=='EXP':f=self.FEXP_HEL
				elif dtyp=='SIM':f=self.FSIM
				for vst in self.VSTS:
					for seq in seql:
						if dtyp=='SIM' and (hel=='POS' or hel=='NEG'): continue
						h5[q2bin_le,wbin_le,dtyp,vst,seq,hel]=f.Get("%s/VST%d/%s/h5_%s"%(q2wbin,vst,seq,hel))
							# if dtyp=='EXP':
							# 	h5[q2bin_le,wbin_le,dtyp,vst,seq,hel]=f.Get("%s/VST%d/%s/h5_%s"%(q2wbin,vst,seq,hel))
							# 	h5[q2bin_le,wbin_le,dtyp,vst,seq,hel]=f.Get("%s/VST%d/%s/h5_%s"%(q2wbin,vst,seq,hel))

		print "keys in h5:"
		print h5.keys()

		#! Write out the list of "bad" q2w bins and the reason why it is bad
		fout=open("%s/bad-q2wbins.txt"%self.OUTDIR_OBS_R2,"w")
		fout.write("Following are the \"bad\" q2w bins:\n")
		for k in q2wbinl_bad:
			fout.write("%s:%s\n"%(k,q2wbinl_bad[k]))
		fout.close()
		print "Finished getting h5s. Now going to display R2s"
		#self.plot_obs_R2(h5)
		self.plot_obs_R2(h5,R2,hel,dtypl,seql)
		#plot_obs_R2(self,h5l,R2,hel,dtypl,seql)
		print "Done DispYields::disp_R2()"
		print "If the progam is not terminating, then Python is probably doing \"garbage collection\"(?); Wait a while!"
		return

	def disp_integ_yield(self,seql=['C','F'],norm=False):
		"""
		Walk the ROOT file and plot y(w;seq,q2bin). 
		"""
		if norm==False:
			outdir=os.path.join(self.OUTDIR,"Obs_IntgYld")
		else:
			outdir=os.path.join(self.OUTDIR,"Obs_IntgYld_norm")
		if not os.path.exists(outdir):
			os.makedirs(outdir)

		#! 1. Get all q2wbins
		q2wbinl=self.get_q2wbinlist(q2max=2.25)
		#print q2wbinl

		q2bng=self.get_q2bng(q2wbinl)
		#print q2bng['BINW']
		print "Yield(w) will be plotted for the following Q2 bins:"
		print q2bng['BINS']

		#! 3. Now put together dictionary for yield: y{seq,q2bin:{w:yield}}
		#! First create dict
		y={}
		oy={} #Will be used later to store yields ordered
		for seq in seql:
			for q2wbin in q2wbinl:
				q2bin=q2wbin.split('_')[0]
				y[seq,q2bin]={}
				oy[seq,q2bin]={}
		#! Now fill dict
		for seq in seql:
			for q2wbin in q2wbinl:
				print "seq,q2wbin,norm=",seq,q2wbin,norm
				q2bin=q2wbin.split('_')[0]
				#! Get w, yield
				w=float(q2wbin.split('_')[1].split('-')[0])
				#wbin=q2wbin.split('_')[1]
				h5_UNPOL=self.FEXP.Get("%s/VST1/%s/h5_UNPOL"%(q2wbin,seq))
				y[seq,q2bin][w]=thntool.GetIntegral(h5_UNPOL)
				if norm==True:
					q2=float(q2wbin.split('_')[0].split('-')[0])
					normf=LUM*LUM_INVFB_TO_INVMICROB*getvgflux(w,q2)#!mub^-1
					print "yield=",y[seq,q2bin][w]
					print "norm=",normf
					y[seq,q2bin][w]=y[seq,q2bin][w]/normf
					print "yield after norm=",y[seq,q2bin][w]
		#! Make sure y[seq,q2bin:w] are sorted by w
		for k in y.keys():
			oy[k]=OrderedDict(sorted(y[k].items()))
			
		#! 4. Now plot
		fig=plt.figure()
		ax=plt.subplot(111)
		clrd={'1.25-1.75':'red','1.75-2.25':'brown','2.25-2.75':'magenta','2.75-3.25':'orange',
		      '3.25-3.75':'yellow','3.75-4.25':'green','4.25-4.75':'cyan','4.75-5.25':'blue'}
		mrkrd={'C':'o','F':'^'}
		for k in oy.keys():
			seq=k[0]
			q2wbin=k[1]
			lbl='%s:[%s)'%(seq,q2wbin)
			ax.scatter(oy[k].keys(),oy[k].values(),label=lbl,c=clrd[q2wbin],marker=mrkrd[seq],s=50)
		if norm==False:
			ax.set_ylim(0,600000)
			ax.set_ylabel(r'Yield [A.U.]',fontsize='xx-large')
		else:
			ax.set_ylim(0,0.05)
			ax.set_ylabel(r'$\mu b$',fontsize='xx-large')
		ax.legend()
		fig.savefig('%s/integ_yield.png'%(outdir))	
			

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


	def get_q2wbinlist(self,q2min=0.00,q2max=6.00,wmin=0.000,wmax=3.000,dbg=False,dbg_bins=10):
		"""
		The ROOT file is arranged in a Tree Structure. The
		Observable histograms (obs-hists) are located as files in the following directory-path(dirpath):
		q2wbin/vst/seq/hists
		"""
		q2wbinl=[]
		i=0 #! for dbg_bins
		for path,dirs,files in self.FEXP.walk():
			if path=="":continue #! Avoid root path
			path_arr=path.split("/")
			if len(path_arr)==1:
				q2wbinl.append(path)
				i+=1
			if dbg==True:
				if i>=dbg_bins: break #! Uncomment/comment -> Get limited q2w-bins/Get all q2w-bins

		#! Remove "particular" q2wbins as implemented below
		# q2max=2.25#1.75
		# wmin=0.00#1.600
		# wmax=3.00#2.200
		q2wbins_remove=[]
		for q2wbin in q2wbinl:
			q2bin_le=q2wbin.split("_")[0].split("-")[0]
			q2bin_ue=q2wbin.split("_")[0].split("-")[1]
			wbin_le =q2wbin.split("_")[1].split("-")[0]
			wbin_ue =q2wbin.split("_")[1].split("-")[1]
			if float(q2bin_ue)<=q2min or float(q2bin_le)>=q2max or float(wbin_ue)<=wmin or float(wbin_le)>=wmax:
				q2wbins_remove.append(q2wbin)
		for q2wbin in q2wbins_remove:
			q2wbinl.remove(q2wbin)

		return q2wbinl


	def get_q2bng(self,q2wbinl):
		return self.get_xbng('Q2',q2wbinl)
	def get_wbng(self,q2wbinl):
		return self.get_xbng('W',q2wbinl)

	def get_xbng(self,x,q2wbinl):
		"""
		Gets x=(Q2 or W) binning information from q2wbinl
		"""

		#! 1. From q2wbinl, identify xbng 
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

	def russ_norm_theta_dist(self,hTheta):
		#! 1. Create Russian-normalization factor histogram
		hDCosTheta=hTheta.Clone()
		nbins=hTheta.GetNbinsX()
		for ibin in range(nbins):
			theta_a=hTheta.GetBinLowEdge(ibin+1)
			theta_b=hTheta.GetBinLowEdge(ibin+1) + hTheta.GetBinWidth(ibin+1)
			DCosTheta=math.fabs(math.cos(math.radians(theta_b))-math.cos(math.radians(theta_a)))
			hDCosTheta.SetBinContent(ibin+1,DCosTheta)
			hDCosTheta.SetBinError(ibin+1,0.)
		#! Now divide hTheta by hDCosTheta
		#! Do Sumw2() so that errors are correctly propagated
		hTheta.Sumw2();
		hTheta.Divide(hDCosTheta)

	#! First make sure this bin is "good"
	def is_bad_q2wbin(self,q2wbin,q2wbinl_bad):
		is_badbin=False
		reason=''
		for vst in self.VSTS:
			h5_UNPOL=self.FEXP_HEL.Get("%s/VST%d/R/h5_UNPOL"%(q2wbin,vst))
			if thntool.GetIntegral(h5_UNPOL)==0:
				reason+="ER=0 for VST%d;"%vst
				is_badbin=True
		if is_badbin:
			q2wbinl_bad[q2wbin]=reason
		return is_badbin
		



