import os
from collections import OrderedDict
import ROOT

VARS=['THETA','PHI']
H2_DIM=OrderedDict([('THETA',0),('PHI',1)])
def run():
        #! Get all input delast files
	FIN={}
	FIN['ER']=ROOT.TFile(os.path.join(os.environ['DELASTDIR_EXP'],'delastR.root'))
	FIN['SR']=ROOT.TFile(os.path.join(os.environ['DELASTDIR_SIM'],'siml','delastR.root'))
	FIN['ST']=ROOT.TFile(os.path.join(os.environ['DELASTDIR_SIM'],'siml','delastT.root'))
	print FIN['ER'].GetName()
	
	FOUT=ROOT.TFile(os.path.join(os.environ['OBSDIR_ELASTIC'],'yield.root'),"RECREATE")

	#! Begin calculation
	hN2={}
	for seq in ['ST','SR','SA','SC','ER','EC']:
		print seq
		if seq=='ST' or seq=='SR' or seq=='ER':#then directly make hists
			print "passed"
			fin=FIN[seq]
			hN2[seq]=fin.Get("delast/yield")
			FOUT.mkdir(seq).cd()
			hN2[seq].Write()
			#! h2
			h2=hN2[seq].Projection(H2_DIM['PHI'],H2_DIM['THETA'])
			h2.SetName("h2")
			h2.Write()
			#! h1s
			h1={}
			for var in VARS:
				h1[var]=hN2[seq].Projection(H2_DIM[var])
				h1[var].SetName("h%s"%var)
				h1[var].Write()
		elif seq=='SA':#then calculate Acceptance from Simulation
			hN2[seq]=hN2['SR'].Clone()
			hN2[seq].Divide(hN2['ST'])
			FOUT.mkdir(seq).cd()
			hN2[seq].Write()
			#! h2
			h2=hN2[seq].Projection(H2_DIM['PHI'],H2_DIM['THETA'])
                        h2.SetName("h2")
                        h2.Write()
			#! h1s
                        h1={}
                        for var in VARS:
                                h1[var]=hN2[seq].Projection(H2_DIM[var])
                                h1[var].SetName("h%s"%var)
                                h1[var].Write()
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
			#! h1s
                        h1={}
                        for var in VARS:
                                h1[var]=hN2[seq].Projection(H2_DIM[var])
                                h1[var].SetName("h%s"%var)
                                h1[var].Write()
			
def disp():
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
	      #('SF'):ROOT.gROOT.ProcessLine("kRed"),#Using SC instead for the moment
	      ('SC'):ROOT.gROOT.ProcessLine("kRed")}
   
	#! Get all histograms and set up their aesthetics
	h={}
        for seq in ['ST','SR','SC','ER','EC']:
                for var in VARS:
                        h[seq,var]=FIN.Get("%s/h%s"%(seq,var))
			h[seq,var].SetLineColor(coll[seq])
			h[seq,var].SetMarkerColor(coll[seq])
			h[seq,var].GetXaxis().SetRangeUser(10,50)
	c=ROOT.TCanvas()
	#hSTn=h['ST','THETA'].DrawNormalized("",1000)
	#hSCn=h['SC','THETA'].DrawNormalized("sames",1000)
	h['ST','THETA'].Draw()#DrawNormalized("",1000)
        h['SC','THETA'].Draw("sames")#DrawNormalized("sames",1000)
	#! Set the minimum and maximum of y coordinate of histograms
	#maxl=[hSTn.GetMaximum(),hSCn.GetMaximum()]
	maxl=[h['ST','THETA'].GetMaximum(),h['SC','THETA'].GetMaximum()]
	maximum=max(maxl)
	#for htmp in [hSTn,hSCn]:
	for htmp in [h['ST','THETA'],h['SC','THETA']]:
		htmp.SetMinimum(0.)
		htmp.SetMaximum(maximum+10)
	c.SaveAs("%s/basic_checkout.png"%OUTDIR)	
