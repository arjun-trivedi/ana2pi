#!/usr/bin/python
from __future__ import division

import os,sys
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict
import math
import numpy as np

def get_fphi(extract_D,name_sfx=""):
	#! Approach where D*sin(phi) is added/removed depending on extract_D
	if extract_D: #! fphi with D*sin(phi)
		fphi=ROOT.TF1("fphi%s"%name_sfx, "([0] + [1]*cos(x*TMath::DegToRad()) + [2]*cos(2*x*TMath::DegToRad()) + [3]*sin(x*TMath::DegToRad()))",0,360)
	else: #! fphi without D*sin(theta)
		fphi=ROOT.TF1("fphi%s"%name_sfx, "([0] + [1]*cos(x*TMath::DegToRad()) + [2]*cos(2*x*TMath::DegToRad()))",0,360)

	#! Set parameter names
	fphi.SetParName(0, "A")
	fphi.SetParName(1, "B")
	fphi.SetParName(2, "C")
	if extract_D:
		fphi.SetParName(3, "D")#hPD
	#! Set initial state of parameters
	fphi.SetParameter(0,.1)#!
	fphi.SetParameter(1,.1)#10
	fphi.SetParameter(2,.1)#20
	if extract_D:
		fphi.SetParameter(3,.1)#100
	
	return fphi

def get_R2_mthd1(R2,h):
	if R2!='A' and R2!='B' and R2!='C' and R2!='D':
		print "get_R2_mthd1() not implemented for R2=",R2

	hc=h.Clone("hc")
	hc.Reset()

	nbins=h.GetNbinsX()
	for ibin in range(nbins):
		phi=math.radians(h.GetBinLowEdge(ibin+1))
		binc=h.GetBinContent(ibin+1)
		bine=h.GetBinError(ibin+1)
		if R2=='A':
			fctr=1
		elif R2=='B':
			fctr=math.cos(phi)
		elif R2=='C':
			fctr=math.cos(2*phi)
		elif R2=='D':
			fctr=math.sin(phi)
		hc.SetBinContent(ibin+1,binc*fctr)
                hc.SetBinError(ibin+1,bine*fctr)

	if R2=='A':
		EF=2*math.pi
		CF=1
	elif R2=='B':#! for B,C and D
		EF=math.pi
		CF=0.93549 #!(1/1.06896)
	elif R2=='C':
		EF=math.pi
		CF=0.75683 #!(1/1.32131)
		scl_fctr=math.pi*CF
	elif R2=='D':
		EF=math.pi
		CF=0.93549 #! (1/1.06896)
                scl_fctr=math.pi*CF

	itg_err=ROOT.Double(0)
	itg=hc.IntegralAndError(1,hc.GetNbinsX(),itg_err)
	return [itg/(EF*CF),itg_err]
	#return hc.Integral()/scl_fctr
		

#! Create histogram
h=ROOT.TH1F("h","h",10,0,360)

#! Fill with fphi distribution with no 'D'
#! First define fphi 
fphi=get_fphi(False)
fphi.SetParameter(0,400000)
fphi.SetParameter(1,20000)
fphi.SetParameter(2,30000)
#! Now fill hist
for i in range(5000000):
	h.Fill(fphi.GetRandom())

#! Test fitting h
#! 1. Fit with no 'D'
fit_fphi_nD=get_fphi(False,"_nD")
h.Fit(fit_fphi_nD,"N")
#! 1. Fit with 'D'
fit_fphi_wD=get_fphi(True,"_wD")
h.Fit(fit_fphi_wD,"N")

#! Draw
ROOT.gStyle.SetOptFit(1111)
c=ROOT.TCanvas("cfit","cfit")
c.Divide(2,1)

c.cd(1)
h.Draw()
fit_fphi_nD.Draw("same")
pt_nD=ROOT.TPaveText(.20,.6,.50,1.0,"NDC")
pt_nD.AddText("nD")
pt_nD.AddText( "A=%.2f+/-%.2f"%(fit_fphi_nD.GetParameter(0),fit_fphi_nD.GetParError(0)) )
pt_nD.AddText( "B=%.2f+/-%.2f"%(fit_fphi_nD.GetParameter(1),fit_fphi_nD.GetParError(1)) )
pt_nD.AddText( "C=%.2f+/-%.2f"%(fit_fphi_nD.GetParameter(2),fit_fphi_nD.GetParError(2)) )
pt_nD.AddText( "#chi^{2}/NDOF=%.2f"%(fit_fphi_nD.GetChisquare()/fit_fphi_nD.GetNDF()) )
pt_nD.Draw()

c.cd(2)
h.Draw()
fit_fphi_wD.Draw("same")
pt_wD=ROOT.TPaveText(.20,.6,.50,1.0,"NDC")
pt_wD.AddText("wD")
pt_wD.AddText( "A=%.2f+/-%.2f"%(fit_fphi_wD.GetParameter(0),fit_fphi_wD.GetParError(0)) )
pt_wD.AddText( "B=%.2f+/-%.2f"%(fit_fphi_wD.GetParameter(1),fit_fphi_wD.GetParError(1)) )
pt_wD.AddText( "C=%.2f+/-%.2f"%(fit_fphi_wD.GetParameter(2),fit_fphi_wD.GetParError(2)) )
pt_wD.AddText( "D=%.2f+/-%.2f"%(fit_fphi_wD.GetParameter(3),fit_fphi_wD.GetParError(3)) )
pt_wD.AddText( "#chi^{2}/NDOF=%.2f"%(fit_fphi_wD.GetChisquare()/fit_fphi_wD.GetNDF()) )
pt_wD.Draw()
	
#! from mthd1
print "*** mthd1 ***"
print "A(err)=%.2f(%.2f)"%(get_R2_mthd1('A',h)[0],get_R2_mthd1('A',h)[1])
print "B(err)=%.2f(%.2f)"%(get_R2_mthd1('B',h)[0],get_R2_mthd1('B',h)[1])
print "C(err)=%.2f(%.2f)"%(get_R2_mthd1('C',h)[0],get_R2_mthd1('C',h)[1])
print "D(err)=%.2f(%.2f)"%(get_R2_mthd1('D',h)[0],get_R2_mthd1('D',h)[1])
print "******"

#! from mthd2
CF_A=1.59155 #!=1/math.radians(36)
CF_B=1.61803
CF_C=1.70131
CF_D=1.61803
#! nD
print "*** mthd2 nD ***"
print "A(err)=%.2f(%.2f)"%(fit_fphi_nD.GetParameter(0)*CF_A,fit_fphi_nD.GetParError(0))
print "B(err)=%.2f(%.2f)"%(fit_fphi_nD.GetParameter(1)*CF_B,fit_fphi_nD.GetParError(1))
print "C(err)=%.2f(%.2f)"%(fit_fphi_nD.GetParameter(2)*CF_C,fit_fphi_nD.GetParError(2))
print "******"
#! from mthd2:wD
print "*** wD ***"
print "A(err)=%.2f(%.2f)"%(fit_fphi_wD.GetParameter(0)*CF_A,fit_fphi_wD.GetParError(0))
print "B(err)=%.2f(%.2f)"%(fit_fphi_wD.GetParameter(1)*CF_B,fit_fphi_wD.GetParError(1))
print "C(err)=%.2f(%.2f)"%(fit_fphi_wD.GetParameter(2)*CF_C,fit_fphi_wD.GetParError(2))
print "D(err)=%.2f(%.2f)"%(fit_fphi_wD.GetParameter(3)*CF_D,fit_fphi_wD.GetParError(3))
print "******"


#! If wanting to keep TCanvas open till program exits				
if not ROOT.gROOT.IsBatch():
	plt.show()
	# wait for you to close the ROOT canvas before exiting
	wait(True)

#if __name__ == "__main__":	
#	if len(sys.argv)==2:
#		plot_fid(sys.argv[1])
#	elif len(sys.argv)==3:
#		plot_fid(sys.argv[1],int(sys.argv[2])) 
	
