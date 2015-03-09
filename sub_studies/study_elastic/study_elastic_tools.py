from __future__ import division
import math
import os
from collections import OrderedDict
import ROOT

import elaslib

VARS=['THETA','PHI']
H2_DIM=OrderedDict([('THETA',0),('PHI',1)])

#PHI_PROJ_BINS={1:[358,2],2:[58,62],3:[118,122],4:[178,182],5:[238,242],6:[298,302]}
#PHI_PROJ_BINS={1:[358,360],2:[58,60],3:[118,120],4:[178,180],5:[238,240],6:[298,300]}

PHI_PROJ_BINS={1:[54,56],2:[56,58],3:[58,60],4:[60,62],5:[62,64],6:[64,66]}
#PHI_PROJ_BINS={1:[234,236],2:[236,238],3:[238,240],4:[240,242],5:[242,244],6:[246,248]}
DPHI=4#degrees
#! For Luminosity & vgflux
LUM=19.844 #fb^-1
LUM_INVFB_TO_INVMICROB=1000000000

def proc_yield():
        #! Get all input delast files
	FIN={}
	FIN['ER']=ROOT.TFile(os.path.join(os.environ['DELASTDIR_EXP'],'delastR.root'))
	FIN['SR']=ROOT.TFile(os.path.join(os.environ['DELASTDIR_SIM'],'siml','delastR.root'))
	FIN['ST']=ROOT.TFile(os.path.join(os.environ['DELASTDIR_SIM'],'siml','delastT.root'))
	#print FIN['ER'].GetName()
	
	FOUT=ROOT.TFile(os.path.join(os.environ['OBSDIR_ELASTIC'],'yield.root'),"RECREATE")

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
		#! Direct projection on to THETA and PHI 
		seqdir=FOUT.GetDirectory(seq)
		seqdir.cd()
		h={}
		for var in VARS:
			h[var]=hN2[seq].Projection(H2_DIM[var])
			h[var].SetName("h%s"%var)
			h[var].Write()
		#! Projection on to THETA for PHI_PROJ_BINS
		for binnum in PHI_PROJ_BINS:
			seqdir.mkdir("phibin%d"%binnum).cd()	
			phi_min=PHI_PROJ_BINS[binnum][0]
			phi_max=PHI_PROJ_BINS[binnum][1]
			print "binnum:phi_min:phi_max=",binnum,phi_min,phi_max
			if binnum==1:#Have to seperately project out this bin in  2 parts because phibin1=[358,2]
				bin_min_1=hN2[seq].GetAxis(H2_DIM['PHI']).FindBin(phi_min+2/2)
				bin_max_1=bin_min_1
				hN2[seq].GetAxis(H2_DIM['PHI']).SetRange(bin_min_1,bin_max_1)
				hTHETA_1=hN2[seq].Projection(H2_DIM['THETA'])
                       		bin_min_2=hN2[seq].GetAxis(H2_DIM['PHI']).FindBin(phi_max-2/2)
				bin_max_2=bin_min_2
				hN2[seq].GetAxis(H2_DIM['PHI']).SetRange(bin_min_2,bin_max_2)
				hTHETA_2=hN2[seq].Projection(H2_DIM['THETA'])	
				hTHETA=hTHETA_1.Clone()
				hTHETA.Add(hTHETA_2)
				hTHETA.SetName("hTHETA")
				hTHETA.Write()
			else:
				bin_min=hN2[seq].GetAxis(H2_DIM['PHI']).FindBin(phi_min+2/2)
				bin_max=hN2[seq].GetAxis(H2_DIM['PHI']).FindBin(phi_max-2/2)
				hN2[seq].GetAxis(H2_DIM['PHI']).SetRange(bin_min,bin_max)
				hTHETA=hN2[seq].Projection(H2_DIM['THETA'])
				hTHETA.SetName("hTHETA")
				hTHETA.Write()
				
def disp_yields():
	#!get rid of X error bars and y error bar caps
	ROOT.gStyle.SetErrorX(0.001)

        FIN=ROOT.TFile(os.path.join(os.environ['OBSDIR_ELASTIC'],'yield.root'))	
	OUTDIR=os.path.join(os.environ['OBSDIR_ELASTIC'],'displays')
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
		for binnum in PHI_PROJ_BINS:
                	print "binnum=",binnum
                        hTHETA[seq,"phibin%d"%binnum]=FIN.Get("%s/phibin%d/hTHETA"%(seq,binnum))
			hTHETA[seq,"phibin%d"%binnum].SetLineColor(cold[seq])
			hTHETA[seq,"phibin%d"%binnum].SetMarkerColor(cold[seq])
			hTHETA[seq,"phibin%d"%binnum].SetMarkerStyle(mrkrd[seq])
			#hTHETA[seq,"phibin%d"%binnum].GetXaxis().SetRangeUser(12,50)
	#! Generate theoretical Cross Sections
        nbins=hTHETA['EC','phibin1'].GetNbinsX()
        xmin=hTHETA['EC','phibin1'].GetXaxis().GetXmin()
        xmax=hTHETA['EC','phibin1'].GetXaxis().GetXmax()
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
	
	#! Checkout Simulation
	c=ROOT.TCanvas()
	c.Divide(3,2)
	for i,binnum in enumerate(PHI_PROJ_BINS):
		c.cd(i+1)
		#hSTn=h['ST','THETA'].DrawNormalized("",1000)
		#hSCn=h['SC','THETA'].DrawNormalized("sames",1000)
		hTHETA['ST',"phibin%d"%binnum].Draw()#DrawNormalized("",1000)
        	hTHETA['SC',"phibin%d"%binnum].Draw("sames")#DrawNormalized("sames",1000)
        	hTHETA['SR',"phibin%d"%binnum].Draw("sames")
		#! Set the minimum and maximum of y coordinate of histograms
		#maxl=[hSTn.GetMaximum(),hSCn.GetMaximum()]
		#maxl=[h['ST','THETA'].GetMaximum(),h['SC','THETA'].GetMaximum(),h['SR','THETA'].GetMaximum()]
		#maximum=max(maxl)
		#for htmp in [hSTn,hSCn]:
		#for htmp in [h['ST','THETA'],h['SC','THETA'],h['SR','THETA']]:
			#htmp.SetMinimum(0.)
			#htmp.SetMaximum(maximum+10)
	c.SaveAs("%s/c_ST_SR_SC.png"%OUTDIR)	

	#! Checkout that the distributions for ER and SR match
	#c=ROOT.TCanvas()
	for i,binnum in enumerate(PHI_PROJ_BINS):
		c.cd(i+1)
        	hERn=hTHETA['ER',"phibin%d"%binnum].DrawNormalized("",1000)
        	hSRn=hTHETA['SR',"phibin%d"%binnum].DrawNormalized("sames",1000)
		#print hSRn.GetMaximum()
        	#! Set the minimum and maximum of y coordinate of histograms
        	maxl=[hERn.GetMaximum(),hSRn.GetMaximum()]
        	maximum=max(maxl)
		#print maxl,maximum
        	for htmp in [hERn,hSRn]:
                	htmp.SetMinimum(0.)
                	htmp.SetMaximum(maximum+10)
       	c.SaveAs("%s/c_ER_SR.png"%OUTDIR)  

	#! Checkout that the distributions for EC and SC match
	hrto_EC_SC={}
	for i,binnum in enumerate(PHI_PROJ_BINS):
		c.cd(i+1)
        	hECn=hTHETA['EC',"phibin%d"%binnum].DrawNormalized("",1000)
        	hSCn=hTHETA['SC',"phibin%d"%binnum].DrawNormalized("sames",1000)
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
                hrto_EC_SC["phibin%d"%binnum]=hECn.Clone()
                hrto_EC_SC["phibin%d"%binnum].Divide(hSCn)
       	c.SaveAs("%s/c_EC_SC.png"%OUTDIR)
        #! Plot difference
        for i,binnum in enumerate(PHI_PROJ_BINS):
                c.cd(i+1)
		hrto_EC_SC["phibin%d"%binnum].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
		hrto_EC_SC["phibin%d"%binnum].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
                hrto_EC_SC["phibin%d"%binnum].SetMaximum(6)
                #hrto["phibin%d"%binnum].GetXaxis().SetRangeUser(12,50)
                hrto_EC_SC["phibin%d"%binnum].SetYTitle("")
                hrto_EC_SC["phibin%d"%binnum].SetTitle("#frac{EC}{SC}")
                hrto_EC_SC["phibin%d"%binnum].Draw()
        c.SaveAs("%s/c_rto_EC_SC.png"%OUTDIR)

	#! Normalize EC to Luminosity and calculate rto between ECln and theoretical yield
	hrto={}
	for i,binnum in enumerate(PHI_PROJ_BINS):
                c.cd(i+1)
		hECln=hTHETA['EC',"phibin%d"%binnum]#.Clone() #! Doing Clone() hear breaks things?!
		hECln.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
		hECln.SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
		hECln.SetYTitle("#frac{d\sigma}{d\Omega} [\mu b]")
		for ibin in range(hECln.GetNbinsX()):
			binc=hECln.GetBinContent(ibin+1)
			binerr=hECln.GetBinError(ibin+1)
			theta=hECln.GetBinLowEdge(ibin+1)
			if theta==0: continue
			dtheta=hECln.GetBinWidth(ibin+1)
			dOmega=( math.sin(math.radians(theta)) )*(dtheta*math.pi/180)*(DPHI*math.pi/180)
			#print "theta,dtheta,dphi,dOmega=",theta,dtheta*math.pi/180,DPHI*math.pi/180,dOmega
			norm=LUM*LUM_INVFB_TO_INVMICROB*dOmega
			#print "norm=",norm
			hECln.SetBinContent(ibin+1,binc/norm)
			hECln.SetBinError(ibin+1,binerr/norm)
		maxl=[hECln.GetMaximum(),hThrtcl.GetMaximum()]
                maximum=max(maxl)
                for htmp in [hECln,hThrtcl]:
                        htmp.SetMinimum(0.)
                        htmp.SetMaximum(maximum)
		hECln.Draw()
                hThrtcl.Draw("sames e")
		#! Calculate difference append from theory
		hECln.Sumw2()
                hThrtcl.Sumw2()
		hrto["phibin%d"%binnum]=hECln.Clone()
		hrto["phibin%d"%binnum].Divide(hThrtcl)
	c.SaveAs("%s/c_EC_lumnorm.png"%OUTDIR)
	#! Plot difference
	for i,binnum in enumerate(PHI_PROJ_BINS):
                c.cd(i+1)
		hrto["phibin%d"%binnum].SetMaximum(6)
		#hrto["phibin%d"%binnum].GetXaxis().SetRangeUser(12,50)
		hrto["phibin%d"%binnum].SetYTitle("")
		hrto["phibin%d"%binnum].SetTitle("#frac{d\sigma}{d\Omega}(#frac{EC}{Theoretical})")
		hrto["phibin%d"%binnum].Draw()
	c.SaveAs("%s/c_rto_ECln_thrtcl.png"%OUTDIR)
