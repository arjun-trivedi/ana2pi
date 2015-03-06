from __future__ import division
import os
from collections import OrderedDict
import ROOT

VARS=['THETA','PHI']
H2_DIM=OrderedDict([('THETA',0),('PHI',1)])

PHI_PROJ_BINS={1:[358,2],2:[58,62],3:[118,122],4:[178,182],5:[238,242],6:[298,302]}

#! For Luminosity & vgflux
LUM=19.844 #fb^-1
LUM_INVFB_TO_INVMICROB=1000000000

def proc_yield():
        #! Get all input delast files
	FIN={}
	FIN['ER']=ROOT.TFile(os.path.join(os.environ['DELASTDIR_EXP'],'delastR.root'))
	FIN['SR']=ROOT.TFile(os.path.join(os.environ['DELASTDIR_SIM'],'siml','delastR.root'))
	FIN['ST']=ROOT.TFile(os.path.join(os.environ['DELASTDIR_SIM'],'siml','delastT.root'))
	print FIN['ER'].GetName()
	
	FOUT=ROOT.TFile(os.path.join(os.environ['OBSDIR_ELASTIC'],'yield.root'),"RECREATE")

	#! Begin process of obtaining Raw and Acceptance corrected yields

	#! 1. Get all hN2s, calculate Acceptances and cast into h2 for easy viewing
	hN2={}
	for seq in ['ST','SR','SA','SC','ER','EC']:
		print seq
		if seq=='ST' or seq=='SR' or seq=='ER':#then directly get hN2
			print "passed"
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
			print "binnum=",binnum
			seqdir.mkdir("phibin%d"%binnum).cd()	
			phi_min=PHI_PROJ_BINS[binnum][0]
			print "phi_min=",phi_min
			phi_max=PHI_PROJ_BINS[binnum][1]
			print "phi_max=",phi_max
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
        FIN=ROOT.TFile(os.path.join(os.environ['OBSDIR_ELASTIC'],'yield.root'))	
	OUTDIR=os.path.join(os.environ['OBSDIR_ELASTIC'],'displays')
	if not os.path.exists(OUTDIR):
		os.makedirs(OUTDIR)
	
	#! "Basic Checkout"
	#! Set up histogram aesthetics			
	coll={('ER'):ROOT.gROOT.ProcessLine("kYellow"),
	      ('EC'):ROOT.gROOT.ProcessLine("kCyan"),
	      #('EF'):ROOT.gROOT.ProcessLine("kBlue"),
	      #('EH'):ROOT.gROOT.ProcessLine("kBlack"),
	      ('ST'):ROOT.gROOT.ProcessLine("kGreen"),
	      ('SR'):ROOT.gROOT.ProcessLine("kMagenta"),
	      #('SF'):ROOT.gROOT.ProcessLine("kRed"),#Using SC instead till Hole Filling 
	      ('SC'):ROOT.gROOT.ProcessLine("kRed")}
   
	#! Get all hTHETA in PHI_PROJ_BINS and set up their aesthetics
	for binnum in PHI_PROJ_BINS:
		hTHETA={}
		for seq in ['ST','SR','SC','ER','EC']:	
                        print "binnum=",binnum
                        hTHETA[seq]=FIN.Get("%s/phibin%d/hTHETA"%(seq,binnum))
			hTHETA[seq].SetLineColor(coll[seq])
			hTHETA[seq].SetMarkerColor(coll[seq])
			hTHETA[seq].GetXaxis().SetRangeUser(10,50)
		outdir=os.path.join(OUTDIR,"phibin%d"%binnum)
        	if not os.path.exists(outdir):
                	os.makedirs(outdir)
		#! Checkout Simulation
		c_ST_SR_SC=ROOT.TCanvas()
		#hSTn=h['ST','THETA'].DrawNormalized("",1000)
		#hSCn=h['SC','THETA'].DrawNormalized("sames",1000)
		hTHETA['ST'].Draw()#DrawNormalized("",1000)
        	hTHETA['SC'].Draw("sames")#DrawNormalized("sames",1000)
        	hTHETA['SR'].Draw("sames")
		#! Set the minimum and maximum of y coordinate of histograms
		#maxl=[hSTn.GetMaximum(),hSCn.GetMaximum()]
		#maxl=[h['ST','THETA'].GetMaximum(),h['SC','THETA'].GetMaximum(),h['SR','THETA'].GetMaximum()]
		#maximum=max(maxl)
		#for htmp in [hSTn,hSCn]:
		#for htmp in [h['ST','THETA'],h['SC','THETA'],h['SR','THETA']]:
			#htmp.SetMinimum(0.)
			#htmp.SetMaximum(maximum+10)
		c_ST_SR_SC.SaveAs("%s/c_ST_SR_SC.png"%outdir)	

		#! Checkout that the distributions for ER and SR match
		c_ER_SR=ROOT.TCanvas()
        	hERn=hTHETA['ER'].DrawNormalized("",1000)
        	hSRn=hTHETA['SR'].DrawNormalized("sames",1000)
		#print hSRn.GetMaximum()
        	#! Set the minimum and maximum of y coordinate of histograms
        	maxl=[hERn.GetMaximum(),hSRn.GetMaximum()]
        	maximum=max(maxl)
		print maxl,maximum
        	for htmp in [hERn,hSRn]:
                	htmp.SetMinimum(0.)
                	htmp.SetMaximum(maximum+10)
        	c_ER_SR.SaveAs("%s/c_ER_SR.png"%outdir)  

		#! Checkout that the distributions for EC and SC match
        	c_EC_SC=ROOT.TCanvas()
        	hECn=hTHETA['EC'].DrawNormalized("",1000)
        	hSCn=hTHETA['SC'].DrawNormalized("sames",1000)
        	#print hSCn.GetMaximum()
        	#! Set the minimum and maximum of y coordinate of histograms
        	maxl=[hECn.GetMaximum(),hSCn.GetMaximum()]
        	maximum=max(maxl)
        	print maxl,maximum
        	for htmp in [hECn,hSCn]:
                	htmp.SetMinimum(0.)
                	htmp.SetMaximum(maximum+10)
        	c_EC_SC.SaveAs("%s/c_EC_SC.png"%outdir)

		#! Normalize EC to Luminosity
		c_EC_lumnorm=ROOT.TCanvas()
		hECln=hTHETA['EC'].Clone()
		hECln.Sumw2();
		hECln.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
		hECln.SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
		norm=LUM*LUM_INVFB_TO_INVMICROB
		hECln.Scale(1/norm)
		hECln.Draw()
		c_EC_lumnorm.SaveAs("%s/c_EC_lumnorm.png"%outdir)

	
