from __future__ import division
import math
import os
from collections import OrderedDict
import ROOT
from rootpy.plotting import Hist, HistStack, Legend, Canvas
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt

import elaslib

'''
The following code is to: delast->yields->display
'''
VARS=['THETA','PHI']
H2_DIM=OrderedDict([('THETA',0),('PHI',1)])

SCTRL=[1,2,3,4,5,6]

PHI_PROJ_BINS_THETA_MIN_14={
				1:[[354,356],[356,358],[358,360],[0,2],    [2,4],    [4,6]],
				2:[[54,56],  [56,58],  [58,60],  [60,62],  [62,64],  [64,66]],
				3:[[114,116],[116,118],[118,120],[120,122],[122,124],[124,126]],
				4:[[174,176],[176,178],[178,180],[180,182],[182,184],[184,186]],
				5:[[234,236],[236,238],[238,240],[240,242],[242,244],[244,246]],
				6:[[294,296],[296,298],[298,300],[300,302],[302,304],[304,306]]
}
DPHI=2
DTHETA=1

#! For Luminosity & vgflux
LUM=19.844 #fb^-1
LUM_INVFB_TO_INVMICROB=1000000000

#! For TCanvas
CANVAS_W=1000
CANVAS_H=800

class StudyElasticTools:
	def __init__(self,theta_min):
		self.THETA_MIN=theta_min
		#! Acoording to THETA_MIN, set up PHI_PROJ_BINS
		if   self.THETA_MIN==14: self.PHI_PROJ_BINS=PHI_PROJ_BINS_THETA_MIN_14
		elif self.THETA_MIN==25: self.PHI_PROJ_BINS=PHI_PROJ_BINS_THETA_MIN_25

	def proc_yield(self):
		#! Get all input delast files
		FIN={}
		FIN['ER']=ROOT.TFile(os.path.join(os.environ['DELASTDIR_EXP'],'delastR.root'))
		if self.THETA_MIN==14:
			FIN['SR']=ROOT.TFile(os.path.join(os.environ['DELASTDIR_SIM'],'siml','delastR.root'))
			FIN['ST']=ROOT.TFile(os.path.join(os.environ['DELASTDIR_SIM'],'siml','delastT.root'))
			FOUT=ROOT.TFile(os.path.join(os.environ['OBSDIR_ELASTIC'],'yield.root'),"RECREATE")
		elif self.THETA_MIN==25:
			FIN['SR']=ROOT.TFile(os.path.join(os.environ['DELASTDIR_SIM'],'siml_theta_min_25','delastR.root'))
			FIN['ST']=ROOT.TFile(os.path.join(os.environ['DELASTDIR_SIM'],'siml_theta_min_25','delastT.root'))
			FOUT=ROOT.TFile(os.path.join(os.environ['OBSDIR_ELASTIC_THETA_MIN_25'],'yield.root'),"RECREATE")
		#print FIN['ER'].GetName()
	
		
		#! Begin process of obtaining Raw and Acceptance corrected yields

		#! 1. Get all hN2s, calculate Acceptances and cast into h2 for easy viewing
		hN2={}
		for seq in ['ST','SR','SA','SC','ER','EC']:
			if seq=='ST' or seq=='SR' or seq=='ER':#then directly get hN2
				fin=FIN[seq]
				hN2[seq]=fin.Get("delast/yield")
				FOUT.mkdir(seq).cd()
				hN2[seq].Write()
			elif seq=='SA':#then calculate Acceptance from Simulation
				hN2[seq]=hN2['SR'].Clone()
				hN2[seq].Divide(hN2['ST'])
				FOUT.mkdir(seq).cd()
				hN2[seq].Write()
			elif seq=='SC' or seq=='EC':#then calculate SC and EC from SA
				if   seq=='SC':hN2[seq]=hN2['SR'].Clone()
				elif seq=='EC':hN2[seq]=hN2['ER'].Clone()
				hN2[seq].Divide(hN2['SA'])
				FOUT.mkdir(seq).cd()
				hN2[seq].Write()
			#! h2
			h2=hN2[seq].Projection(H2_DIM['PHI'],H2_DIM['THETA'])
			h2.SetName("h2")
			h2.Write()	

		#! 2. Now make projection on to THETA and PHI: Direct and on to THETA for PHI_PROJ_BINS
		for seq in ['ST','SR','SA','SC','ER','EC']:
			#! First project out THETA from the appropriate minimum THETA simulated
			bin_min=hN2[seq].GetAxis(H2_DIM['THETA']).FindBin(self.THETA_MIN+DTHETA/2)
			hN2[seq].GetAxis(H2_DIM['THETA']).SetRange(bin_min)
			#! Direct projection on to THETA and PHI 
			seqdir=FOUT.GetDirectory(seq)
			seqdir.cd()
			h={}
			for var in VARS:
				h[var]=hN2[seq].Projection(H2_DIM[var])
				h[var].SetName("h%s"%var)
				h[var].Write()
			#! Projection on to THETA for PHI_PROJ_BINS
			for sector in self.PHI_PROJ_BINS:
				sectordir=seqdir.mkdir("sector%d"%sector)
				sectordir.cd()
				print self.PHI_PROJ_BINS[sector]
				for iphibinnum,phibin in enumerate(self.PHI_PROJ_BINS[sector]):
					print phibin
					phi_min=phibin[0]
					phi_max=phibin[1]
					print "sector:phibinnum:phibin:phi_min:phi_max=",sector,iphibinnum+1,phibin,phi_min,phi_max
					#phibindir=sectordir.mkdir("%d-%d"%(phi_min,phi_max)).cd()
					phibindir=sectordir.mkdir("phibinnum%d"%(iphibinnum+1)).cd()
					bin_min=hN2[seq].GetAxis(H2_DIM['PHI']).FindBin(phi_min+DPHI/2)
					bin_max=hN2[seq].GetAxis(H2_DIM['PHI']).FindBin(phi_max-DPHI/2)
					hN2[seq].GetAxis(H2_DIM['PHI']).SetRange(bin_min,bin_max)
					hTHETA=hN2[seq].Projection(H2_DIM['THETA'])
					hTHETA.SetName("hTHETA")
					hTHETA.SetTitle("#theta projection for #phi=[%d,%d)"%(phi_min,phi_max))
					hTHETA.Write()
				
	def disp_yields(self,thry='wrad'):
		#!get rid of X error bars and y error bar caps
		ROOT.gStyle.SetErrorX(0.001)
		FIN=ROOT.TFile(os.path.join(os.environ['OBSDIR_ELASTIC'],'yield.root'))	
		if thry=='wrad':
			OUTDIR=os.path.join(os.environ['OBSDIR_ELASTIC'],'displays_thry_%s'%thry)
		elif thry=='nrad':
			OUTDIR=os.path.join(os.environ['OBSDIR_ELASTIC'],'displays_thry_%s'%thry)
		else:
			sys.exit('thry has to be wrad or nrad. Entered value= %s'%thry)
		if not os.path.exists(OUTDIR):
			os.makedirs(OUTDIR)
	
		#! Set up histogram aesthetics			
		cold={('ER'):ROOT.gROOT.ProcessLine("kYellow"),
			  ('EC'):ROOT.gROOT.ProcessLine("kCyan"),
			  #('EF'):ROOT.gROOT.ProcessLine("kBlue"),
			  #('EH'):ROOT.gROOT.ProcessLine("kBlack"),
			  ('ST'):ROOT.gROOT.ProcessLine("kGreen"),
			  ('SR'):ROOT.gROOT.ProcessLine("kMagenta"),
			  #('SF'):ROOT.gROOT.ProcessLine("kRed"),#Using SC instead till Hole Filling 
			  ('SC'):ROOT.gROOT.ProcessLine("kRed")}

		mrkrd={('ER'):ROOT.gROOT.ProcessLine("kOpenStar"),
			   ('EC'):ROOT.gROOT.ProcessLine("kFullCircle"),
			  #('EF'):ROOT.gROOT.ProcessLine("kBlue"),
			  #('EH'):ROOT.gROOT.ProcessLine("kBlack"),
			  ('ST'):ROOT.gROOT.ProcessLine("kPlus"),
			  ('SR'):ROOT.gROOT.ProcessLine("kOpenStar"),
			  #('SF'):ROOT.gROOT.ProcessLine("kRed"),#Using SC instead till Hole Filling 
			  ('SC'):ROOT.gROOT.ProcessLine("kFullCircle")}
   
		#! Get all hTHETA in PHI_PROJ_BINS and set up their aesthetics
		hTHETA={}
		for seq in ['ST','SR','SC','ER','EC']:
			for sector in self.PHI_PROJ_BINS:
				for iphibinnum,phibin in enumerate(self.PHI_PROJ_BINS[sector]):
					phi_min=phibin[0]
					phi_max=phibin[1]
					#print "seq:sector:phibinnum:phibin:phi_min:phi_max=",seq,sector,iphibinnum+1,phibin,phi_min,phi_max
					hTHETA[seq,"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)]=FIN.Get("%s/sector%d/phibinnum%d/hTHETA"%(seq,sector,(iphibinnum+1)))
					hTHETA[seq,"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetLineColor(cold[seq])
					hTHETA[seq,"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetMarkerColor(cold[seq])
					hTHETA[seq,"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetMarkerStyle(mrkrd[seq])
	
		#! Generate theoretical Cross Sections
		#! hTTnorm=dsigma/dOmega
		#! hTT=sigma; this histogram is helpful in direct comparison with Yields
		E1F_TARGET_RAD_LENGTH=0.00562 #! see nb_play_elaslib
		nbins=hTHETA['EC','sector1','phibinnum1'].GetNbinsX()
		xmin= hTHETA['EC','sector1','phibinnum1'].GetXaxis().GetXmin()
		xmax= hTHETA['EC','sector1','phibinnum1'].GetXaxis().GetXmax()
		hTTnorm=ROOT.TH1F("hTTnorm","hTTnorm",nbins,xmin,xmax)
		hTTnorm.SetYTitle("d\sigma [#frac{\mu b}{sr}]")
		hTT=ROOT.TH1F("hTT","hTT",nbins,xmin,xmax)
		hTT.SetYTitle("#d\sigma} [\mu b]")
		for ibin in range(nbins):
			theta=hTTnorm.GetBinLowEdge(ibin+1)
			if theta==0: continue
			if thry=='wrad':
				binc=elaslib.elasrad(5.499,theta,E1F_TARGET_RAD_LENGTH,1.1)
			elif thry=='nrad':
				binc=elaslib.elas(5.499,theta)
			#print "theta,xsec=",theta,binc
			hTTnorm.SetBinContent(ibin+1,binc)
			hTTnorm.SetBinError(ibin+1,0)
			#hTTnorm.GetXaxis().SetRangeUser(12,50)
			theta=hTTnorm.GetBinLowEdge(ibin+1)
			if theta==0: continue
			dtheta=hTTnorm.GetBinWidth(ibin+1)
			dOmega=( math.sin(math.radians(theta)) )*(dtheta*math.pi/180)*(DPHI*math.pi/180)
			hTT.SetBinContent(ibin+1,binc*dOmega)
			hTT.SetBinError(ibin+1,0)
	
		#! Now plot hTHETA in PHI_PROJ_BINS
		
		#! First create a TCanvas & TLegend
		#! NOTE, that this TCanvas and TLegend are reused
		c=ROOT.TCanvas("c","c",CANVAS_W,CANVAS_H)
		num_phibins=len(self.PHI_PROJ_BINS[1])
		nrows=2
		ncols=int(num_phibins/nrows)
		print "nrows,ncols",nrows,ncols
		c.Divide(ncols,nrows)
		l=ROOT.TLegend(0.45,0.6,1.0,0.85)
		l.SetFillStyle(0)
		l.SetBorderSize(0)
		l.SetTextSize(0.03)

		#! 1. Checkout ST comparison with TT
		hrto_ST_thry={}
		outdir=os.path.join(OUTDIR,"ST_TT")
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		print "outdir=",outdir
		for sector in self.PHI_PROJ_BINS:
			for iphibinnum,phibin in enumerate(self.PHI_PROJ_BINS[sector]):
				pad=c.cd(iphibinnum+1)
				pad.SetLogy()
				hSTn=hTHETA['ST',"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].DrawNormalized("",1000)
				hTTn=hTT.DrawNormalized("sames",1000)
				if (iphibinnum+1)==1: #then draw TLegend
					l.Clear()
					l.AddEntry(hSTn,"ST-norm")
					l.AddEntry(hTTn,"d\sigma_{thry-%s}-norm"%thry)
					l.Draw()
				#! Calculate ratio
				hSTn.Sumw2()
				hrto_ST_thry["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)]=hSTn.Clone()
				hrto_ST_thry["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].Divide(hTTn)
			c.SaveAs("%s/c_sector%d.png"%(outdir,sector))
		#! Plot ratio
		outdir=os.path.join(OUTDIR,"rto_ST_TT")
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		print "outdir=",outdir
		for sector in self.PHI_PROJ_BINS:
			for iphibinnum,phibin in enumerate(self.PHI_PROJ_BINS[sector]):
				pad=c.cd(iphibinnum+1)
				pad.SetLogy(0)
				hrto_ST_thry["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetMinimum(0.)
				hrto_ST_thry["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetMaximum(8.)
				hrto_ST_thry["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
				hrto_ST_thry["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
				#hrto_ST_thry["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetMaximum(6)
				#hrto["phibin%d"%binnum].GetXaxis().SetRangeUser(12,50)
				hrto_ST_thry["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetYTitle("")
				old_title=hrto_ST_thry["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].GetTitle()
				hrto_ST_thry["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetTitle("#frac{ST}{TT} for %s"%old_title)
				hrto_ST_thry["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].Draw()
				if (iphibinnum+1)==1: #then draw TLegend
					l.Clear()
					l.AddEntry(hrto_ST_thry["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)],"#frac{ST-norm}{d\sigma_{thry-%s}-norm}"%thry)
					l.Draw()
			c.SaveAs("%s/c_sector%d.png"%(outdir,sector)) 

		#! 2. Checkout Simulation: ST,SR and SC
		outdir=os.path.join(OUTDIR,"ST_SR_SC")
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		print "outdir=",outdir
		for sector in self.PHI_PROJ_BINS:
			for iphibinnum,phibin in enumerate(self.PHI_PROJ_BINS[sector]):
				pad=c.cd(iphibinnum+1)
				pad.SetLogy()
				hTHETA['ST',"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].Draw()#DrawNormalized("",1000)
				hTHETA['SC',"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].Draw("sames")#DrawNormalized("sames",1000)
				hTHETA['SR',"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].Draw("sames")
				if (iphibinnum+1)==1: #then draw TLegend
					l.Clear()
					l.AddEntry(hTHETA['ST',"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)],"ST")
					l.AddEntry(hTHETA['SC',"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)],"SC")
					l.AddEntry(hTHETA['SR',"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)],"SR")
					l.Draw()
			c.SaveAs("%s/c_sector%d.png"%(outdir,sector))
	
		#! 3. Checkout that the distributions for ER and SR match
		outdir=os.path.join(OUTDIR,"ER_SR")
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		print "outdir=",outdir
		for sector in self.PHI_PROJ_BINS:
			for iphibinnum,phibin in enumerate(self.PHI_PROJ_BINS[sector]):
				pad=c.cd(iphibinnum+1)
				pad.SetLogy()
				hERn=hTHETA['ER',"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].DrawNormalized("",1000)
				hSRn=hTHETA['SR',"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].DrawNormalized("sames",1000)
				#print hSRn.GetMaximum()
				#! Set the minimum and maximum of y coordinate of histograms
				maxl=[hERn.GetMaximum(),hSRn.GetMaximum()]
				maximum=max(maxl)
				#print maxl,maximum
				for htmp in [hERn,hSRn]:
						htmp.SetMinimum(0.01)
						htmp.SetMaximum(maximum+10)
				if (iphibinnum+1)==1: #then draw TLegend
					l.Clear()
					l.AddEntry(hERn,"ER-norm")
					l.AddEntry(hSRn,"SR-norm")
					l.Draw()
			c.SaveAs("%s/c_sector%d.png"%(outdir,sector))  
	
		#! 4. Checkout that the distributions for EC and SC match
		hrto_EC_SC={}
		outdir=os.path.join(OUTDIR,"EC_SC")
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		print "outdir=",outdir
		for sector in self.PHI_PROJ_BINS:
			for iphibinnum,phibin in enumerate(self.PHI_PROJ_BINS[sector]):
				pad=c.cd(iphibinnum+1)
				pad.SetLogy()
				hECn=hTHETA['EC',"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].DrawNormalized("",1000)
				hSCn=hTHETA['SC',"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].DrawNormalized("sames",1000)
				#print hSCn.GetMaximum()
				#! Set the minimum and maximum of y coordinate of histograms
				maxl=[hECn.GetMaximum(),hSCn.GetMaximum()]
				maximum=max(maxl)
				#print maxl,maximum
				for htmp in [hECn,hSCn]:
						htmp.SetMinimum(0.01)
						htmp.SetMaximum(maximum+10)
				if (iphibinnum+1)==1: #then draw TLegend
					l.Clear()
					l.AddEntry(hECn,"EC-norm")
					l.AddEntry(hSCn,"SC-norm")
					l.Draw()
				#! Calculate difference
				hECn.Sumw2()
				hSCn.Sumw2()
				hrto_EC_SC["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)]=hECn.Clone()
				hrto_EC_SC["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].Divide(hSCn)
			c.SaveAs("%s/c_sector%d.png"%(outdir,sector)) 
		#! Plot difference
		outdir=os.path.join(OUTDIR,"rto_EC_SC")
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		print "outdir=",outdir
		for sector in self.PHI_PROJ_BINS:
			for iphibinnum,phibin in enumerate(self.PHI_PROJ_BINS[sector]):
				pad=c.cd(iphibinnum+1)
				pad.SetLogy(0)
				hrto_EC_SC["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
				hrto_EC_SC["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
				hrto_EC_SC["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetMaximum(6)
				#hrto["phibin%d"%binnum].GetXaxis().SetRangeUser(12,50)
				hrto_EC_SC["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetYTitle("")
				old_title=hrto_EC_SC["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].GetTitle()
				hrto_EC_SC["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetTitle("#frac{EC}{SC} for %s"%old_title)
				if (iphibinnum+1)==1: #then draw TLegend
					l.Clear()
					l.AddEntry(hrto_EC_SC["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)],"#frac{EC-norm}{SC-norm}")
					l.Draw()
				hrto_EC_SC["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].Draw()
			c.SaveAs("%s/c_sector%d.png"%(outdir,sector)) 

		#! 5. Checkout EC comparison with TT
		hrto_EC_TT={}
		outdir=os.path.join(OUTDIR,"EC_TT")
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		print "outdir=",outdir
		for sector in self.PHI_PROJ_BINS:
			for iphibinnum,phibin in enumerate(self.PHI_PROJ_BINS[sector]):
				pad=c.cd(iphibinnum+1)
				pad.SetLogy()
				hECn=hTHETA['EC',"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].DrawNormalized("",1000)
				#hTTnormn=hTTnorm.DrawNormalized("sames",1000)
				hTTn=hTT.DrawNormalized("sames",1000)
				#! Calculate ratio
				hECn.Sumw2()
				hrto_EC_TT["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)]=hECn.Clone()
				hrto_EC_TT["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].Divide(hTTn)
				if (iphibinnum+1)==1: #then draw TLegend
					l.Clear()
					l.AddEntry(hECn,"EC-norm")
					l.AddEntry(hTTn,"d\sigma_{thry-%s}-norm"%thry)
					l.Draw()
			c.SaveAs("%s/c_sector%d.png"%(outdir,sector))
		#! Plot ratio
		outdir=os.path.join(OUTDIR,"rto_EC_TT")
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		print "outdir=",outdir
		for sector in self.PHI_PROJ_BINS:
			for iphibinnum,phibin in enumerate(self.PHI_PROJ_BINS[sector]):
				pad=c.cd(iphibinnum+1)
				pad.SetLogy(0)
				hrto_EC_TT["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetMinimum(0.)
				hrto_EC_TT["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetMaximum(8.)
				hrto_EC_TT["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
				hrto_EC_TT["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
				#hrto_EC_TT["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetMaximum(6)
				#hrto["phibin%d"%binnum].GetXaxis().SetRangeUser(12,50)
				hrto_EC_TT["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetYTitle("")
				old_title=hrto_EC_TT["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].GetTitle()
				hrto_EC_TT["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetTitle("#frac{EC}{TT} for %s"%old_title)
				hrto_EC_TT["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].Draw()
				if (iphibinnum+1)==1: #then draw TLegend
					l.Clear()
					l.AddEntry(hrto_EC_TT["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)],"#frac{EC-norm}{d\sigma_{thry-%s}-norm}"%thry)
					l.Draw()
			c.SaveAs("%s/c_sector%d.png"%(outdir,sector))
	
		#! 6. Normalize EC to Luminosity and calculate rto between EClumnorm and TTnorm
		hrto_EClumnorm_TTnorm={}
		outdir=os.path.join(OUTDIR,"EClumnorm_TTnorm")
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		print "outdir=",outdir
		for sector in self.PHI_PROJ_BINS:
			for iphibinnum,phibin in enumerate(self.PHI_PROJ_BINS[sector]):
				pad=c.cd(iphibinnum+1)
				pad.SetLogy(1)
				hEClumnorm=hTHETA['EC',"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)]#.Clone() #! Doing Clone() hear breaks things?!
				hEClumnorm.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
				hEClumnorm.SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
				hEClumnorm.SetYTitle("#frac{d\sigma}{d\Omega} [#frac{\mu b}{sr}]")
				old_title=hEClumnorm.GetTitle()
				hEClumnorm.SetTitle("Differential Cross Section for %s"%old_title)
				for ibin in range(hEClumnorm.GetNbinsX()):
					binc=hEClumnorm.GetBinContent(ibin+1)
					binerr=hEClumnorm.GetBinError(ibin+1)
					theta=hEClumnorm.GetBinLowEdge(ibin+1)
					if theta==0: continue
					dtheta=hEClumnorm.GetBinWidth(ibin+1)
					dOmega=( math.sin(math.radians(theta)) )*(dtheta*math.pi/180)*(DPHI*math.pi/180)
					#print "theta,dtheta,dphi,dOmega=",theta,dtheta*math.pi/180,DPHI*math.pi/180,dOmega
					norm=LUM*LUM_INVFB_TO_INVMICROB*dOmega
					#print "norm=",norm
					hEClumnorm.SetBinContent(ibin+1,binc/norm)
					hEClumnorm.SetBinError(ibin+1,binerr/norm)

				maxl=[hEClumnorm.GetMaximum(),hTTnorm.GetMaximum()]
				maximum=max(maxl)
				for htmp in [hEClumnorm,hTTnorm]:
					htmp.SetMinimum(0.00000001)
					htmp.SetMaximum(maximum)
				hEClumnorm.Draw()
				hTTnorm.Draw("sames e")
				if (iphibinnum+1)==1: #then draw TLegend
					l.Clear()
					l.AddEntry(hEClumnorm,"#frac{d\sigma}{d\Omega}^{EC}")
					l.AddEntry(hTTnorm,"#frac{d\sigma_{thry-%s}}{d\Omega}"%thry)
					l.Draw()
				#! Calculate difference append from theory
				hEClumnorm.Sumw2()
				hTTnorm.Sumw2()
				hrto_EClumnorm_TTnorm["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)]=hEClumnorm.Clone()
				hrto_EClumnorm_TTnorm["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].Divide(hTTnorm)
			c.SaveAs("%s/c_sector%d.png"%(outdir,sector)) 
		#! Plot ratio
		outdir=os.path.join(OUTDIR,"rto_EClumnorm_TTnorm")
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		print "outdir=",outdir
		for sector in self.PHI_PROJ_BINS:
			for iphibinnum,phibin in enumerate(self.PHI_PROJ_BINS[sector]):
				pad=c.cd(iphibinnum+1)
				pad.SetLogy(0)
				hrto_EClumnorm_TTnorm["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetMaximum(6)
				#hrto_EClumnorm_TTnorm["phibin%d"%binnum].GetXaxis().SetRangeUser(12,50)
				hrto_EClumnorm_TTnorm["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetYTitle("")
				old_title=hrto_EClumnorm_TTnorm["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].GetTitle()
				hrto_EClumnorm_TTnorm["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetTitle("#frac{#frac{d\sigma}{d\Omega}^{EC}}{#frac{d\sigma_{thry-%s}}{d\Omega}} for %s"%(thry,old_title))
				hrto_EClumnorm_TTnorm["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].Draw()
				if (iphibinnum+1)==1: #then draw TLegend
					l.Clear()
					l.AddEntry(hrto_EClumnorm_TTnorm["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)],"#frac{#frac{d\sigma}{d\Omega}^{EC}}{#frac{d\sigma_{thry-%s}}{d\Omega}}"%thry)
					l.Draw()
			c.SaveAs("%s/c_sector%d.png"%(outdir,sector))
'''
The following code is used to play with elaslib.f
'''
#! Following are parameters for obtain legnth of E1F target in radiation length.
#! Sources=Rad_length_materials_src1_CERN.pdf, Rad_length_materials_src2_LBNL.pdf
X0=63.047 #g/cm^2 for atomic hydrogrm
DENSITY=0.0708 #g/cm^3 for liquid H2
RAD_LENGTH=890.45 #cm = X0/DENSITY
E1F_TARGET_RAD_LENGTH=0.00562 # =5cm/RAD_LENGTH 
def play_elaslib_low_theta(be,nbins,xmin,xmax,wcut,rto_min=None,rto_max=None):
	h={}
	cold={'nrad':'kRed','wrad':'kBlue'}
	for radtyp in ['nrad','wrad']:
		h[radtyp]=Hist(nbins,xmin,xmax,title='%s'%radtyp)
		h[radtyp].SetLineColor(ROOT.gROOT.ProcessLine(cold[radtyp]))
		h[radtyp].SetMarkerColor(ROOT.gROOT.ProcessLine(cold[radtyp]))
	for ibin in range(nbins):
		theta=h['nrad'].GetBinLowEdge(ibin+1)
		if theta==0: continue
		binc_nrad=elaslib.elas(be,theta)
		binc_wrad=elaslib.elasrad(be,theta,E1F_TARGET_RAD_LENGTH,wcut)
		#print "theta,xsec=",theta,binc
		h['nrad'].SetBinContent(ibin+1,binc_nrad)
		h['nrad'].SetBinError(ibin+1,0)
		if math.isnan(binc_wrad):
			h['wrad'].SetBinContent(ibin+1,0)
			h['wrad'].SetBinError(ibin+1,0)
		else:
			h['wrad'].SetBinContent(ibin+1,binc_wrad)
			h['wrad'].SetBinError(ibin+1,0)

	hrto=h['wrad'].Clone()
	#hrto.SetTitle('wrad/nrad_%.2f'%wcut)
	hrto.Divide(h['nrad'])
	hrto.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
	hrto.SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
	return h['nrad'],h['wrad'],hrto

def play_elaslib(be,nbins=179,xmin=1,xmax=180,wcutl=[0.98,1.00,1.10,1.50,2.00,2.50,3.00,3.50,4.00,4.50,5.00],rto_min=None,rto_max=None):
	outdir=os.path.join(os.environ['STUDY_ELASTIC_DATADIR'],'elaslib_ana')
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	h={}
	hrto={}
	#! For theta<1
	h_lt={}
        hrto_lt={}

	cold={'nrad':'kRed','wrad':'kBlue'}
	for wcut in wcutl:
		for radtyp in ['nrad','wrad']:
			h[radtyp,wcut]=Hist(nbins,xmin,xmax,title='%s'%radtyp)
			h[radtyp,wcut].SetLineColor(ROOT.gROOT.ProcessLine(cold[radtyp]))
			h[radtyp,wcut].SetMarkerColor(ROOT.gROOT.ProcessLine(cold[radtyp]))
		for ibin in range(nbins):
			theta=h['nrad',wcut].GetBinLowEdge(ibin+1)
			if theta==0: continue
			binc_nrad=elaslib.elas(be,theta)
			binc_wrad=elaslib.elasrad(be,theta,E1F_TARGET_RAD_LENGTH,wcut)
			#print "theta,xsec=",theta,binc
			h['nrad',wcut].SetBinContent(ibin+1,binc_nrad)
			h['nrad',wcut].SetBinError(ibin+1,0)
			if math.isnan(binc_wrad):
        			h['wrad',wcut].SetBinContent(ibin+1,0)
        			h['wrad',wcut].SetBinError(ibin+1,0)
			else:
				h['wrad',wcut].SetBinContent(ibin+1,binc_wrad)
                		h['wrad',wcut].SetBinError(ibin+1,0)
        
		hrto[wcut]=h['wrad',wcut].Clone()
		hrto[wcut].SetTitle('wrad/nrad_%.2f'%wcut)
		hrto[wcut].Divide(h['nrad',wcut])
		hrto[wcut].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
		hrto[wcut].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
   
		#! Get low theta
		h_lt['nrad',wcut],h_lt['wrad',wcut],hrto_lt[wcut]=play_elaslib_low_theta(be,nbins=50,xmin=0.00001,xmax=1,wcut=wcut)
		hrto_lt[wcut].SetTitle('wrad/nrad_%.2f'%wcut)

	#! First plot all Cross Sections for different W cuts in separate figures
	for wcut in wcutl:
		fig=plt.figure(figsize=(10, 6), dpi=100)
                fig.suptitle("Beam energy=%.3f[GeV]. W cut=%.2f[GeV]"%(be,wcut),fontsize=20)
		#! First plot low theta
		ax=plt.subplot(211)
		stack = HistStack()
		stack.Add(h_lt['nrad',wcut])
		stack.Add(h_lt['wrad',wcut])
		rplt.errorbar(stack, xerr=False, emptybins=True, axes=ax)
		ax.set_title("ep Cross section")
		ax.set_ylabel(r"$\frac{d\sigma}{d\Omega}[\mu b]$",fontsize=20)
		ax.set_yscale('log')
		ymin=min(elaslib.elas(be,1),elaslib.elasrad(be,1,E1F_TARGET_RAD_LENGTH,wcut))
		ax.set_ylim(ymin)
		plt.legend()
		#! Now plot high theta
		ax=plt.subplot(212)
		stack = HistStack()
		stack.Add(h['nrad',wcut])
		stack.Add(h['wrad',wcut])
		rplt.errorbar(stack, xerr=False, emptybins=True, axes=ax)
		ax.set_title("ep Cross section")
		ax.set_ylabel(r"$\frac{d\sigma}{d\Omega}[\mu b]$",fontsize=20)
		ax.set_yscale('log')
		ymin=min(elaslib.elas(be,xmax),elaslib.elasrad(be,xmax,E1F_TARGET_RAD_LENGTH,wcut))
		ax.set_ylim(ymin)
		plt.legend()
		fig.savefig("%s/xsec_be%.3f_wcut%.2f.png"%(outdir,be,wcut))

	#! Now plot all ratios in a single figure
	fig=plt.figure(figsize=(10, 6), dpi=100)
        fig.suptitle("Beam energy=%.3f[GeV]"%(be),fontsize=20)
	#! First plot low theta
        ax=plt.subplot(211)
        stack=HistStack()
	#! Get appropriate range for histograms
	minl=[]
	maxl=[]
	for wcut in wcutl:
		minl.append( hrto_lt[wcut].GetMinimum())
		maxl.append( hrto_lt[wcut].GetMaximum())
	minimum=min(minl)
	maximum=max(maxl)		
	for iwcut,wcut in enumerate(wcutl):
		hrto_lt[wcut].SetLineColor(ROOT.gROOT.ProcessLine('kViolet+%d'%(iwcut+1)))
		hrto_lt[wcut].SetMarkerColor(ROOT.gROOT.ProcessLine('kViolet+%d'%(iwcut+1)))
		stack.Add(hrto_lt[wcut])
	rplt.errorbar(stack, xerr=False, emptybins=True, axes=ax)
	ax.set_ylim(minimum,maximum)
	plt.legend(loc='lower right')
	#! Now plot high theta
        ax=plt.subplot(212)
        stack=HistStack()
	#! Get appropriate range for histograms
        minl=[]
        maxl=[]
        for wcut in wcutl:
                minl.append( hrto[wcut].GetMinimum())
                maxl.append( hrto[wcut].GetMaximum())
	minimum=min(minl)
        maximum=max(maxl)
        for iwcut,wcut in enumerate(wcutl):
		hrto[wcut].SetLineColor(ROOT.gROOT.ProcessLine('kViolet+%d'%(iwcut+1)))
		hrto[wcut].SetMarkerColor(ROOT.gROOT.ProcessLine('kViolet+%d'%(iwcut+1)))
                stack.Add(hrto[wcut])
        rplt.errorbar(stack, xerr=False, emptybins=True, axes=ax)
	ax.set_ylim(minimum,maximum)
        plt.legend(loc='lower right')
        fig.savefig("%s/rto_be%.3f.png"%(outdir,be))
	