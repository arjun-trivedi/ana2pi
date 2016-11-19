#!/usr/bin/python
from __future__ import division

import os,sys
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict

import numpy as np

#! Get kin from user
if len(sys.argv)<2:
	sys.exit("Please enter kin=p,theta or phi")
kin=sys.argv[1]

#! Get simrgn from user
if len(sys.argv)<3:
	sys.exit("Please enter lQ2 or hQ2")
simrgn=sys.argv[2]
#! As per simrgn setup Q2
if simrgn=="lQ2":
	Q2MIN,Q2MAX=2.00,3.00
elif simrgn=="hQ2":
	Q2MIN,Q2MAX=3.00,5.00
else:
	sys.exit("simrgn not understood")
#! W independent of sim region
WMIN,WMAX=1.400,2.125


top=1

NDTYP=3;
f=[0 for i in range(NDTYP)]# //0=ER,1=SR,2=ST
ER,SR,ST=range(NDTYP);
print ER,SR,ST
dtyp_name=["ER","SR","ST"]
dtyp_name_englsh=["exp-reco","sim-reco","sim-thrn"]
clr=[ROOT.gROOT.ProcessLine("kBlue"),ROOT.gROOT.ProcessLine("kRed"),ROOT.gROOT.ProcessLine("kGreen")]
f[0]=ROOT.TFile("$D2PIDIR_EXP_E16/data_eid_w_eff_scpd_071516/deid.root")
f[1]=ROOT.TFile("$D2PIDIR_SIM_E16/data_eid_w_eff_scpd_071516/deid.root")
f[2]=ROOT.TFile("$D2PIDIR_SIM_E16/data_eid_w_eff_scpd_071516/deidT.root")
t=[0 for i in range(NDTYP)]
for idtyp in range(NDTYP):
	if   idtyp==ER or idtyp==SR: t[idtyp]=f[idtyp].Get("d2piR/tR")
	elif idtyp==ST:              t[idtyp]=f[idtyp].Get("d2piT/tT")

NPART=4;#//0=e,1=p,2=pip,3=pim
part_name=["e","p","pip","pim"]
part_name_ltx=["e","p","#pi^{+}","#pi^{-}"]
h_p=[[0 for j in range(NPART)] for i in range(NDTYP)]#[NDTYP][NPART]; 

#! Global q2w cut
cut_q2w=ROOT.TCut("Q2>%.2f && Q2<%.2f && W>%.3f && W<%.3f"%(Q2MIN,Q2MAX,WMIN,WMAX))

#! Misc global
NENTRIES=100000000
#NENTRIES=10000# debug

c=ROOT.TCanvas("c","c")
for idtyp in range(NDTYP):
	c.cd()
	#! establish top cut
	cut_top=ROOT.TCut("top==%d"%top)
	if idtyp==ST: cut_top="" #!no cut for ST
	#if idtyp==ER: continue #tmp
	#if idtyp==SR: continue #tmp
	#if idtyp==ST: continue #tmp
	#! total cut=top+q2w
	cut=ROOT.TCut(cut_top)
	cut+=cut_q2w
	print("Processing %s\n"%dtyp_name[idtyp]);
	print("cut="),cut.Print()
	#! Now draw
	print("Tree="),t[idtyp].GetName()
	if kin=="p":
		t[idtyp].Draw("p_e>>h_e(100,0,5)",cut,"",NENTRIES)
		t[idtyp].Draw("p_p>>h_p(100,0,5)",cut,"",NENTRIES)
		t[idtyp].Draw("p_pip>>h_pip(100,0,5)",cut,"",NENTRIES)
		t[idtyp].Draw("p_pim>>h_pim(100,0,5)",cut,"",NENTRIES)
	elif kin=="theta":
		t[idtyp].Draw("theta_e>>h_e(60,0,60)",cut,"",NENTRIES)
		t[idtyp].Draw("theta_p>>h_p(60,0,60)",cut,"",NENTRIES)
		t[idtyp].Draw("theta_pip>>h_pip(100,0,180)",cut,"",NENTRIES)
		t[idtyp].Draw("theta_pim>>h_pim(100,0,180)",cut,"",NENTRIES)
	elif kin=="phi":
		t[idtyp].Draw("phi_e>>h_e(100,-50,350)",cut,"",NENTRIES)
		t[idtyp].Draw("phi_p>>h_p(100,-50,350)",cut,"",NENTRIES)
		t[idtyp].Draw("phi_pip>>h_pip(100,-50,350)",cut,"",NENTRIES)
		t[idtyp].Draw("phi_pim>>h_pim(100,-50,350)",cut,"",NENTRIES)
	else:
		print kin,"not implemented"
		sys.exit()

	print("here")#tmp

	ROOT.gStyle.SetOptStat(0)#("ne") 
	h_p[idtyp][0]=ROOT.gDirectory.Get("h_e")
	h_p[idtyp][1]=ROOT.gDirectory.Get("h_p")
	h_p[idtyp][2]=ROOT.gDirectory.Get("h_pip")
	h_p[idtyp][3]=ROOT.gDirectory.Get("h_pim")
	for iprt in range(NPART):
		h_p[idtyp][iprt].SetName("h_%s_%s_%s"%(kin,dtyp_name[idtyp],part_name[iprt]))
		h_p[idtyp][iprt].SetLineColor(clr[idtyp])
		h_p[idtyp][iprt].SetMarkerColor(clr[idtyp])
	print("here2")#tmp

c.Close()

#sys.exit()#tmp

print("here3")#tmp
cname="c_%s_top%d"%(kin,top)
c_p=ROOT.TCanvas(cname,cname);
c_p.Divide(2,2)
#! Setup legend
leg=ROOT.TLegend(0.70,0.70,0.85,0.90)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.05)
for iprt in range(NPART):
	pad=c_p.cd(iprt+1);
	nhist=0
	#! to keep cloned hists
	htmp=[0 for l in range(NDTYP)] #! should not be more than 2 (NDTYP)
	for idtyp in range(NDTYP):
		#if idtyp==ER: continue #tmp
		#if idtyp==SR: continue #tmp
		#if idtyp==ST: continue #tmp
		#! Main title
		if kin=="p": title="%s_{%s}"%(kin,part_name_ltx[iprt])
		else:		 title="#%s_{%s}"%(kin,part_name_ltx[iprt])
		h_p[idtyp][iprt].SetTitle(title)
		#! Axes title
		if kin=="p": xtitle="%s [GeV]"%kin
		else:		 xtitle="#%s [deg]"%kin
		h_p[idtyp][iprt].SetXTitle(xtitle)
		c_p.SetLeftMargin(0.20)
		h_p[idtyp][iprt].GetYaxis().SetTitleOffset(1.5)
		h_p[idtyp][iprt].SetYTitle("N_{entries}")
		#! draw
		if nhist==0:  htmp[nhist]=h_p[idtyp][iprt].DrawNormalized("hist",1000);
		else:         htmp[nhist]=h_p[idtyp][iprt].DrawNormalized("hist sames",1000);
		#! increment nhist
		nhist+=1
	# #! set y axis maximum
	fctr=1.2
	maxl=[]
	print htmp
	for htmp2 in htmp:
		if htmp2==0: continue
		maxl.append(htmp2.GetMaximum())
	maximum=max(maxl)
	for htmp2 in htmp:
		if htmp2==0: continue
		htmp2.SetMinimum(0.)
		htmp2.SetMaximum(maximum*fctr)
	pad.Update()
	#! add legend
	if iprt+1==1:
		print "Drawing legend"
		for idtyp in range(NDTYP):
			#if idtyp==ST: continue
			#l.AddEntry(h_p[idtyp][iprt],dtyp_name_englsh[idtyp],"l")
			leg.AddEntry(htmp[idtyp],dtyp_name_englsh[idtyp],"l")
		leg.Draw("same")
	pad.Update()
#! Save canvas
print("here4")#tmp
outdir="/data/trivedia/e16/study_d2pi_kinematics/%s"%(simrgn)
if not os.path.exists(outdir):
	os.makedirs(outdir)
c_p.SaveAs("%s/%s.png"%(outdir,c_p.GetName()))
c_p.SaveAs("%s/%s.pdf"%(outdir,c_p.GetName()))
  
#! If wanting to keep TCanvas open till program exits                           
#if not ROOT.gROOT.IsBatch():
#        plt.show()
#        # wait for you to close the ROOT canvas before exiting
#        wait(True)

#if __name__ == "__main__":     
#       if len(sys.argv)==2:
#               plot_fid(sys.argv[1])
#       elif len(sys.argv)==3:
#               plot_fid(sys.argv[1],int(sys.argv[2]))
