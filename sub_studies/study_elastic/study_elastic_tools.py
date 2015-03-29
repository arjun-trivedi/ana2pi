from __future__ import division
import math
import os
from collections import OrderedDict
import ROOT
from rootpy.plotting import Hist, HistStack, Legend, Canvas
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
import re

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
	def __init__(self,obsdate,simnum='siml'):
		self.THETA_MIN=14 #! Minimum theta in Simulation
		#! Acoording to THETA_MIN, set up PHI_PROJ_BINS
		self.PHI_PROJ_BINS=PHI_PROJ_BINS_THETA_MIN_14
		#! Setup DATADIR
		self.SIMNUM=simnum
		self.DATADIR=os.path.join(os.path.join(os.environ['OBSDIR_ELASTIC'],obsdate))
		
	def proc_yield(self):
		#! Get all input delast files
		FIN={}
		FIN['ER']=ROOT.TFile(os.path.join(self.DATADIR,'delast_exp','delastR.root'))
		if self.THETA_MIN==14:
			FIN['SR']=ROOT.TFile(os.path.join(self.DATADIR,'delast_sim',self.SIMNUM,'delastR.root'))
			FIN['ST']=ROOT.TFile(os.path.join(self.DATADIR,'delast_sim',self.SIMNUM,'delastT.root'))
			FOUT=ROOT.TFile(os.path.join(self.DATADIR,self.SIMNUM,'yield.root'),"RECREATE")
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
				#print self.PHI_PROJ_BINS[sector]
				for iphibinnum,phibin in enumerate(self.PHI_PROJ_BINS[sector]):
					#print phibin
					phi_min=phibin[0]
					phi_max=phibin[1]
					#print "sector:phibinnum:phibin:phi_min:phi_max=",sector,iphibinnum+1,phibin,phi_min,phi_max
					#phibindir=sectordir.mkdir("%d-%d"%(phi_min,phi_max)).cd()
					phibindir=sectordir.mkdir("phibinnum%d"%(iphibinnum+1)).cd()
					bin_min=hN2[seq].GetAxis(H2_DIM['PHI']).FindBin(phi_min+DPHI/2)
					bin_max=hN2[seq].GetAxis(H2_DIM['PHI']).FindBin(phi_max-DPHI/2)
					hN2[seq].GetAxis(H2_DIM['PHI']).SetRange(bin_min,bin_max)
					hTHETA=hN2[seq].Projection(H2_DIM['THETA'])
					hTHETA.SetName("hTHETA")
					hTHETA.SetTitle("#theta projection for #phi=[%d,%d)"%(phi_min,phi_max))
					hTHETA.Write()
				
	def disp_yields(self,thry='wrad',wcut=1.1):
		#!get rid of X error bars and y error bar caps
		ROOT.gStyle.SetErrorX(0.001)
		FIN=ROOT.TFile(os.path.join(self.DATADIR,self.SIMNUM,'yield.root'))	
		if thry=='wrad':
			self.OUTDIR=os.path.join(self.DATADIR,self.SIMNUM,'displays_thry_%s'%thry)
		elif thry=='nrad':
			self.OUTDIR=os.path.join(self.DATADIR,self.SIMNUM,'displays_thry_%s'%thry)
		else:
			sys.exit('thry has to be wrad or nrad. Entered value= %s'%thry)
		if not os.path.exists(self.OUTDIR):
			os.makedirs(self.OUTDIR)
	
		#! Set up histogram aesthetics			
		cold={('ER'):ROOT.gROOT.ProcessLine("kYellow"),
			  ('EC'):ROOT.gROOT.ProcessLine("kCyan"),
			  ('ST'):ROOT.gROOT.ProcessLine("kGreen"),
			  ('SR'):ROOT.gROOT.ProcessLine("kMagenta"),
			  ('SC'):ROOT.gROOT.ProcessLine("kRed")}

		mrkrd={('ER'):ROOT.gROOT.ProcessLine("kOpenStar"),
			   ('EC'):ROOT.gROOT.ProcessLine("kFullCircle"),
			   ('ST'):ROOT.gROOT.ProcessLine("kPlus"),
			   ('SR'):ROOT.gROOT.ProcessLine("kOpenStar"),
			   ('SC'):ROOT.gROOT.ProcessLine("kFullCircle")}
   
		#! Get all hTHETA in PHI_PROJ_BINS and set up their aesthetics
		self.hTHETA={}
		for seq in ['ST','SR','SC','ER','EC']:
			for sector in self.PHI_PROJ_BINS:
				for iphibinnum,phibin in enumerate(self.PHI_PROJ_BINS[sector]):
					phi_min=phibin[0]
					phi_max=phibin[1]
					#print "seq:sector:phibinnum:phibin:phi_min:phi_max=",seq,sector,iphibinnum+1,phibin,phi_min,phi_max
					self.hTHETA[seq,"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)]=FIN.Get("%s/sector%d/phibinnum%d/hTHETA"%(seq,sector,(iphibinnum+1)))
					self.hTHETA[seq,"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetLineColor(cold[seq])
					self.hTHETA[seq,"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetMarkerColor(cold[seq])
					self.hTHETA[seq,"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetMarkerStyle(mrkrd[seq])
	
		#! Generate theoretical Cross Sections
		E1F_TARGET_RAD_LENGTH=0.00562 #! see nb_play_elaslib
		nbins=self.hTHETA['EC','sector1','phibinnum1'].GetNbinsX()
		xmin= self.hTHETA['EC','sector1','phibinnum1'].GetXaxis().GetXmin()
		xmax= self.hTHETA['EC','sector1','phibinnum1'].GetXaxis().GetXmax()
		self.hTHETA['TTnorm'],self.hTHETA['TT']=self.get_thrtcl_xsec(thry,5.499,E1F_TARGET_RAD_LENGTH,wcut,nbins,xmin,xmax)
			
		#! 1. Checkout ST comparison with TT
		
		self.comp("STnorm","TTnorm",logy=True)
		self.comp("ERnorm","SRnorm",logy=True)
		self.comp("ECnorm","SCnorm",logy=True)
		self.comp("ECnorm","STnorm",logy=True)
		self.comp("ECnorm","TTnorm",logy=True)
		self.comp("EClumnorm","TTnorm",logy=True,draw_normalized=False)

	def get_thrtcl_xsec(self,thry,be,trgt_lgth,wcut,nbins,xmin,xmax):
		'''
		Generate theoretical Cross Sections
		+ hTTnorm=dsigma/DOmega
		+ hTT=dsigma=hTTnorm*dOmega; this histogram is helpful in direct comparison with Yields
		'''
		hTTnorm=ROOT.TH1F("hTTnorm","hTTnorm",nbins,xmin,xmax)
		hTTnorm.SetYTitle("d\sigma [#frac{\mu b}{sr}]")
		hTT=ROOT.TH1F("hTT","hTT",nbins,xmin,xmax)
		hTT.SetYTitle("#d\sigma} [\mu b]")
		
		for ibin in range(nbins):
			theta=hTTnorm.GetBinLowEdge(ibin+1)
			if theta==0: continue
			if thry=='wrad':
				binc=elaslib.elasrad(be,theta,trgt_lgth,wcut)
			elif thry=='nrad':
				binc=elaslib.elas(be,theta)
			#print "theta,xsec=",theta,binc
			hTTnorm.SetBinContent(ibin+1,binc)
			hTTnorm.SetBinError(ibin+1,0)
			theta_bin_low=hTTnorm.GetBinLowEdge(ibin+1)
			theta_bin_max=hTTnorm.GetBinLowEdge(ibin+1)+hTTnorm.GetBinWidth(ibin+1)
			if theta_bin_low==0: continue
			DOmega=self.get_DOmega(theta_bin_low,theta_bin_max,self.PHI_PROJ_BINS[1][0][0],self.PHI_PROJ_BINS[1][1][1])
			hTT.SetBinContent(ibin+1,binc*DOmega)
			hTT.SetBinError(ibin+1,0)

		return hTTnorm,hTT

	def get_DOmega(self,theta_bin_min,theta_bin_max,phi_bin_min,phi_bin_max):
			DCosTheta=math.fabs(math.cos(math.radians(theta_bin_max))-math.cos(math.radians(theta_bin_min)))
			DPhi=math.radians(phi_bin_max)-math.radians(phi_bin_min)
			DOmega=DCosTheta*DPhi
			return DOmega

	def norm_hist(self,h,phi_bin_low,phi_bin_max,use_lum=False):
		for ibin in range(h.GetNbinsX()):
			binc=h.GetBinContent(ibin+1)
			binerr=h.GetBinError(ibin+1)
			theta_bin_min=h.GetBinLowEdge(ibin+1)
			theta_bin_max=h.GetBinLowEdge(ibin+1) + h.GetBinWidth(ibin+1)
			if theta_bin_min==0: continue
			DOmega=self.get_DOmega(theta_bin_min,theta_bin_max,phi_bin_low,phi_bin_max)
			norm=DOmega
			if use_lum:
				norm=DOmega*LUM*LUM_INVFB_TO_INVMICROB
			#print "norm=",norm
			h.SetBinContent(ibin+1,binc/norm)
			h.SetBinError(ibin+1,binerr/norm)

	def comp(self,seqA,seqB,logy=False,draw_normalized=True):#,A_norm=False,A_lumnorm=False,B_norm=False,B_lumnorm=False,logy=False):
		print "Generating comparison plots for %s and %s"%(seqA,seqB)
		#! First create a TCanvas & TLegend
		c=ROOT.TCanvas("c","c",CANVAS_W,CANVAS_H)
		num_phibins=len(self.PHI_PROJ_BINS[1])
		nrows=2
		ncols=int(num_phibins/nrows)
		#print "nrows,ncols",nrows,ncols
		c.Divide(ncols,nrows)
		l=ROOT.TLegend(0.45,0.6,1.0,0.85)
		l.SetFillStyle(0)
		l.SetBorderSize(0)
		l.SetTextSize(0.03)

		outdir=os.path.join(self.OUTDIR,"%s_%s"%(seqA,seqB))
		outdir_rto=os.path.join(self.OUTDIR,"rto_%s_%s"%(seqA,seqB))
		
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		#print "outdir=",outdir

		h={}
		hrto={}
		for sector in self.PHI_PROJ_BINS:
			for iphibinnum,phibin in enumerate(self.PHI_PROJ_BINS[sector]):
				pad=c.cd(iphibinnum+1)
				if logy==True:
					pad.SetLogy()

				for seq in [seqA,seqB]:
					if seq.find("TT")<0:#! if not TT, TTnorm histograms
						if re.match(".{2}norm",seq):#if seq.find("norm")>0:
							h[seq,sector,iphibinnum]=self.hTHETA["%s"%seq.strip("norm"),"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].Clone("seq_%d"%sector) 
							self.norm_hist(h[seq,sector,iphibinnum],phibin[0],phibin[1],use_lum=False)
						if re.match(".{2}lumnorm",seq):#if seq.find("lumnorm")>0:
							h[seq,sector,iphibinnum]=self.hTHETA["%s"%seq.strip("lumnorm"),"sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].Clone("seq_%d"%sector) 
							self.norm_hist(h[seq,sector,iphibinnum],phibin[0],phibin[1],use_lum=True)
						if seq=="EClumnorm":
							h[seq,sector,iphibinnum].SetYTitle("#frac{d\sigma}{d\Omega} [#frac{\mu b}{sr}]")
						else:
							h[seq,sector,iphibinnum].SetYTitle("Yield [A.U.]")
					else:
						h[seq,sector,iphibinnum]=self.hTHETA[seq].Clone()

				if draw_normalized:
					hseqAn=h[seqA,sector,iphibinnum].DrawNormalized("",1000)
					hseqBn=h[seqB,sector,iphibinnum].DrawNormalized("sames",1000)
				else:
					h[seqA,sector,iphibinnum].Draw()
					h[seqB,sector,iphibinnum].Draw("sames")

				if (iphibinnum+1)==1: #then draw TLegend
					l.Clear()
					if draw_normalized:
						l.AddEntry(hseqAn,seqA)
						l.AddEntry(hseqBn,seqB)
					else:
						l.AddEntry(h[seqA,sector,iphibinnum],seqA)
						l.AddEntry(h[seqB,sector,iphibinnum],seqB)
					l.Draw()

				#! Calculate ratio
				if draw_normalized:
					hseqAn.Sumw2()
					hseqBn.Sumw2()
					hrto["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)]=hseqAn.Clone()
					hrto["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].Divide(hseqBn)
				else:
					h[seqA,sector,iphibinnum].Sumw2()
					h[seqB,sector,iphibinnum].Sumw2()
					hrto["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)]=h[seqA,sector,iphibinnum].Clone()
					hrto["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].Divide(h[seqB,sector,iphibinnum])
			c.SaveAs("%s/c_sector%d.png"%(outdir,sector))

		#! Plot ratio
		if not os.path.exists(outdir_rto):
			os.makedirs(outdir_rto)
		#print "outdir_rto=",outdir_rto
		for sector in self.PHI_PROJ_BINS:
			for iphibinnum,phibin in enumerate(self.PHI_PROJ_BINS[sector]):
				pad=c.cd(iphibinnum+1)
				pad.SetLogy(0)
				hrto["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetMinimum(0)
				hrto["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetMaximum(6)
				hrto["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
				hrto["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
				hrto["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetYTitle("")
				old_title=hrto["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].GetTitle()
				hrto["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].SetTitle("#frac{%s}{%s} for %s"%(seqA,seqB,old_title))
				hrto["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)].Draw()
				if (iphibinnum+1)==1: #then draw TLegend
					l.Clear()
					l.AddEntry(hrto["sector%d"%sector,"phibinnum%d"%(iphibinnum+1)],"#frac{%s}{%s}"%(seqA,seqB))
					l.Draw()
			c.SaveAs("%s/c_sector%d.png"%(outdir_rto,sector)) 

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
    cold={'norad':'kRed','wrad':'kBlue'}
    for radtyp in ['norad','wrad']:
        h[radtyp]=Hist(nbins,xmin,xmax,title='%s'%radtyp)
        h[radtyp].SetLineColor(ROOT.gROOT.ProcessLine(cold[radtyp]))
        h[radtyp].SetMarkerColor(ROOT.gROOT.ProcessLine(cold[radtyp]))
    for ibin in range(nbins):
        theta=h['norad'].GetBinLowEdge(ibin+1)
        if theta==0: continue
        binc_norad=elaslib.elas(be,theta)
        binc_wrad=elaslib.elasrad(be,theta,E1F_TARGET_RAD_LENGTH,wcut)
        #print "theta,xsec=",theta,binc
        h['norad'].SetBinContent(ibin+1,binc_norad)
        h['norad'].SetBinError(ibin+1,0)
        if math.isnan(binc_wrad):
                h['wrad'].SetBinContent(ibin+1,0)
                h['wrad'].SetBinError(ibin+1,0)
        else:
                h['wrad'].SetBinContent(ibin+1,binc_wrad)
                h['wrad'].SetBinError(ibin+1,0)

    hrto=h['wrad'].Clone()
    hrto.SetTitle('wrad/norad')
    hrto.Divide(h['norad'])
    hrto.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
    hrto.SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
    return h['norad'],h['wrad'],hrto

def play_elaslib(be,nbins=179,xmin=1,xmax=180,wcut=1.1,rto_min=None,rto_max=None):
    #outdir="/tmp/elastic_xsec_ana"
    outdir=os.path.join(os.environ['ANA2PI_STUDY_ELASTIC'],'elaslib_ana')
    if not os.path.exists(outdir):
	os.makedirs(outdir)
    h={}
    cold={'norad':'kRed','wrad':'kBlue'}
    for radtyp in ['norad','wrad']:
        h[radtyp]=Hist(nbins,xmin,xmax,title='%s'%radtyp)
        h[radtyp].SetLineColor(ROOT.gROOT.ProcessLine(cold[radtyp]))
        h[radtyp].SetMarkerColor(ROOT.gROOT.ProcessLine(cold[radtyp]))
    for ibin in range(nbins):
        theta=h['norad'].GetBinLowEdge(ibin+1)
        if theta==0: continue
        binc_norad=elaslib.elas(be,theta)
        binc_wrad=elaslib.elasrad(be,theta,E1F_TARGET_RAD_LENGTH,wcut)
        #print "theta,xsec=",theta,binc
        h['norad'].SetBinContent(ibin+1,binc_norad)
        h['norad'].SetBinError(ibin+1,0)
	if math.isnan(binc_wrad):
        	h['wrad'].SetBinContent(ibin+1,0)
        	h['wrad'].SetBinError(ibin+1,0)
	else:
		h['wrad'].SetBinContent(ibin+1,binc_wrad)
                h['wrad'].SetBinError(ibin+1,0)
        
    hrto=h['wrad'].Clone()
    hrto.SetTitle('wrad/norad')
    hrto.Divide(h['norad'])
    hrto.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
    hrto.SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
   
    #! Get low theta
    h['norad_lt'],h['wrad_lt'],hrto_lt=play_elaslib_low_theta(be,nbins=50,xmin=0.00001,xmax=1,wcut=wcut)

    #! In a 4x4 grid plot histograms
    #! First plot low theta
    fig=plt.figure(figsize=(10, 6), dpi=100)
    fig.suptitle("Beam energy=%.3f[GeV]. W cut=%.2f[GeV]"%(be,wcut),fontsize=20)
    ax=plt.subplot(221)
    stack = HistStack()
    stack.Add(h['norad_lt'])
    stack.Add(h['wrad_lt'])
    rplt.errorbar(stack, xerr=False, emptybins=True, axes=ax)
    #rplt.errorbar(hThrtcl['wrad'],xerr=False, emptybins=False, axes=ax)
    ax.set_title("ep Cross section")
    ax.set_ylabel(r"$\frac{d\sigma}{d\Omega}[\mu b]$",fontsize=20)
    ax.set_yscale('log')
    ymin=min(elaslib.elas(be,1),elaslib.elasrad(be,1,E1F_TARGET_RAD_LENGTH,wcut))
    ax.set_ylim(ymin)
    plt.legend()
    ax=plt.subplot(222)
    #print "Ratio min:max(low theta)=%.2f,%.2f"%(hrto_lt.GetMinimum(),hrto_lt.GetMaximum())
    rplt.errorbar(hrto_lt, xerr=False, emptybins=True, axes=ax)
    ax.set_ylim(hrto_lt.GetMinimum(),hrto_lt.GetMaximum())
    plt.legend(loc='lower right')

    #! Now plot high theta
    ax=plt.subplot(223)
    stack = HistStack()
    stack.Add(h['norad'])
    stack.Add(h['wrad'])
    rplt.errorbar(stack, xerr=False, emptybins=True, axes=ax)
    #rplt.errorbar(hThrtcl['wrad'],xerr=False, emptybins=False, axes=ax)
    ax.set_title("ep Cross section")
    ax.set_ylabel(r"$\frac{d\sigma}{d\Omega}[\mu b]$",fontsize=20)
    ax.set_yscale('log')
    ymin=min(elaslib.elas(be,xmax),elaslib.elasrad(be,xmax,E1F_TARGET_RAD_LENGTH,wcut))
    ax.set_ylim(ymin)
    plt.legend()
    ax=plt.subplot(224)
    #print "Ratio min:max=%.2f,%.2f"%(hrto.GetMinimum(),hrto.GetMaximum())
    rplt.errorbar(hrto, xerr=False, emptybins=True, axes=ax)
    if rto_min==None and rto_max==None:
        if hrto.GetMinimum()<1:
    		ax.set_ylim(hrto.GetMinimum(),hrto.GetMaximum())
	else:
        	ax.set_ylim(0.9,hrto.GetMaximum())
    else:
    	ax.set_ylim(rto_min,rto_max)
    plt.legend(loc='upper left')
    fig.savefig("%s/be%.3f_wcut%.2f.png"%(outdir,be,wcut))
    #return hThrtcl['norad'],hThrtcl['wrad'],hrto
