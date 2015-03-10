from __future__ import division
import math
import os
from collections import OrderedDict
import ROOT

import elaslib

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
				
	def disp_yields(self):
		#!get rid of X error bars and y error bar caps
		ROOT.gStyle.SetErrorX(0.001)
		if self.THETA_MIN==14:
			FIN=ROOT.TFile(os.path.join(os.environ['OBSDIR_ELASTIC'],'yield.root'))	
			OUTDIR=os.path.join(os.environ['OBSDIR_ELASTIC'],'displays')
		elif self.THETA_MIN==25:
			FIN=ROOT.TFile(os.path.join(os.environ['OBSDIR_ELASTIC_THETA_MIN_25'],'yield.root'))	
			OUTDIR=os.path.join(os.environ['OBSDIR_ELASTIC_THETA_MIN_25'],'displays')
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
		nbins=hTHETA['EC','sector1','phibinnum1'].GetNbinsX()
		xmin= hTHETA['EC','sector1','phibinnum1'].GetXaxis().GetXmin()
		xmax= hTHETA['EC','sector1','phibinnum1'].GetXaxis().GetXmax()
		hThrtcl=ROOT.TH1F("hThrtcl","hThrtcl",nbins,xmin,xmax)
		hThrtcl.SetYTitle("#frac{d\sigma}{d\Omega} [\mu b]")
		for ibin in range(nbins):
			theta=hThrtcl.GetBinLowEdge(ibin+1)
			if theta==0: continue
			#binc=elaslib.elas(5.499,theta)
			binc=elaslib.elasrad(5.499,theta,2.5*0.00577,1.1)
			#print "theta,xsec=",theta,binc
			hThrtcl.SetBinContent(ibin+1,binc)
			hThrtcl.SetBinError(ibin+1,0)
			#hThrtcl.GetXaxis().SetRangeUser(12,50)
	
		#! Now plot hTHETA in PHI_PROJ_BINS
		
		#! First create a TCanvas
		#! NOTE, that this TCanvas is reused
		c=ROOT.TCanvas("c","c",CANVAS_W,CANVAS_H)
		num_phibins=len(self.PHI_PROJ_BINS[1])
		nrows=2
		ncols=int(num_phibins/nrows)
		print "nrows,ncols",nrows,ncols
		c.Divide(ncols,nrows)

		#! Checkout ST comparison with Theory
		outdir=os.path.join(OUTDIR,"ST_Theory")
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		print "outdir=",outdir
		for sector in self.PHI_PROJ_BINS:
			for iphibinnum,phibin in enumerate(self.PHI_PROJ_BINS[sector]):
				pad=c.cd(iphibinnum+1)
				pad.SetLogy()
				hTHETA['ST',"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].DrawNormalized("",1000)
				hThrtcl.DrawNormalized("sames",1000)
			c.SaveAs("%s/c_sector%d.png"%(outdir,sector))

		#! Checkout Simulation: ST,SR and SC
		outdir=os.path.join(OUTDIR,"ST_SR_SC")
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		print "outdir=",outdir
		for sector in self.PHI_PROJ_BINS:
			for iphibinnum,phibin in enumerate(self.PHI_PROJ_BINS[sector]):
				pad=c.cd(iphibinnum+1)
				pad.SetLogy(0)
				hTHETA['ST',"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].Draw()#DrawNormalized("",1000)
				hTHETA['SC',"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].Draw("sames")#DrawNormalized("sames",1000)
				hTHETA['SR',"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].Draw("sames")
				#! Set the minimum and maximum of y coordinate of histograms
				#maxl=[hSTn.GetMaximum(),hSCn.GetMaximum()]
				#maxl=[h['ST','THETA'].GetMaximum(),h['SC','THETA'].GetMaximum(),h['SR','THETA'].GetMaximum()]
				#maximum=max(maxl)
				#for htmp in [hSTn,hSCn]:
				#for htmp in [h['ST','THETA'],h['SC','THETA'],h['SR','THETA']]:
					#htmp.SetMinimum(0.)
					#htmp.SetMaximum(maximum+10)
			c.SaveAs("%s/c_sector%d.png"%(outdir,sector))
	
		# #! Checkout that the distributions for ER and SR match
		outdir=os.path.join(OUTDIR,"ER_SR")
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		print "outdir=",outdir
		for sector in self.PHI_PROJ_BINS:
			for iphibinnum,phibin in enumerate(self.PHI_PROJ_BINS[sector]):
				c.cd(iphibinnum+1)
				hERn=hTHETA['ER',"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].DrawNormalized("",1000)
				hSRn=hTHETA['SR',"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].DrawNormalized("sames",1000)
				#print hSRn.GetMaximum()
				#! Set the minimum and maximum of y coordinate of histograms
				maxl=[hERn.GetMaximum(),hSRn.GetMaximum()]
				maximum=max(maxl)
				#print maxl,maximum
				for htmp in [hERn,hSRn]:
						htmp.SetMinimum(0.)
						htmp.SetMaximum(maximum+10)
			c.SaveAs("%s/c_sector%d.png"%(outdir,sector))  
	
		# #! Checkout that the distributions for EC and SC match
		hrto_EC_SC={}
		outdir=os.path.join(OUTDIR,"EC_SC")
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		print "outdir=",outdir
		for sector in self.PHI_PROJ_BINS:
			for iphibinnum,phibin in enumerate(self.PHI_PROJ_BINS[sector]):
				c.cd(iphibinnum+1)
				hECn=hTHETA['EC',"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].DrawNormalized("",1000)
				hSCn=hTHETA['SC',"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].DrawNormalized("sames",1000)
				#print hSCn.GetMaximum()
				#! Set the minimum and maximum of y coordinate of histograms
				maxl=[hECn.GetMaximum(),hSCn.GetMaximum()]
				maximum=max(maxl)
				#print maxl,maximum
				for htmp in [hECn,hSCn]:
						htmp.SetMinimum(0.)
						htmp.SetMaximum(maximum+10)
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
				c.cd(iphibinnum+1)
				hrto_EC_SC["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
				hrto_EC_SC["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
				hrto_EC_SC["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetMaximum(6)
				#hrto["phibin%d"%binnum].GetXaxis().SetRangeUser(12,50)
				hrto_EC_SC["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetYTitle("")
				old_title=hrto_EC_SC["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].GetTitle()
				hrto_EC_SC["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetTitle("#frac{EC}{SC} for %s"%old_title)
				hrto_EC_SC["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].Draw()
			c.SaveAs("%s/c_sector%d.png"%(outdir,sector)) 
	
		#! Normalize EC to Luminosity and calculate rto between EClumnorm and theoretical yield
		hrto_EClumnorm_thrtcl={}
		outdir=os.path.join(OUTDIR,"EClumnorm_thrtcl")
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
				hEClumnorm.SetYTitle("#frac{d\sigma}{d\Omega} [\mu b]")
				old_title=hEClumnorm.GetTitle()
				hEClumnorm.SetTitle("#frac{d\sigma}{d\Omega} [\mu b] for %s"%old_title)
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

				maxl=[hEClumnorm.GetMaximum(),hThrtcl.GetMaximum()]
				maximum=max(maxl)
				for htmp in [hEClumnorm,hThrtcl]:
					htmp.SetMinimum(0.00000001)
					htmp.SetMaximum(maximum)
				hEClumnorm.Draw()
				hThrtcl.Draw("sames e")
				#! Calculate difference append from theory
				hEClumnorm.Sumw2()
				hThrtcl.Sumw2()
				hrto_EClumnorm_thrtcl["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)]=hEClumnorm.Clone()
				hrto_EClumnorm_thrtcl["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].Divide(hThrtcl)
			c.SaveAs("%s/c_sector%d.png"%(outdir,sector)) 
		#! Plot ratio
		outdir=os.path.join(OUTDIR,"rto_EClumnorm_thrtcl")
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		print "outdir=",outdir
		for sector in self.PHI_PROJ_BINS:
			for iphibinnum,phibin in enumerate(self.PHI_PROJ_BINS[sector]):
				pad=c.cd(iphibinnum+1)
				pad.SetLogy(0)
				hrto_EClumnorm_thrtcl["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetMaximum(6)
				#hrto_EClumnorm_thrtcl["phibin%d"%binnum].GetXaxis().SetRangeUser(12,50)
				hrto_EClumnorm_thrtcl["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetYTitle("")
				old_title=hrto_EClumnorm_thrtcl["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].GetTitle()
				hrto_EClumnorm_thrtcl["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetTitle("#frac{EC}{Theoretical} for %s"%old_title)
				hrto_EClumnorm_thrtcl["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].Draw()
			c.SaveAs("%s/c_sector%d.png"%(outdir,sector))