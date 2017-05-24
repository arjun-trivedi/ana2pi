#!/usr/bin/python

from __future__ import division
import os,datetime,sys
import study_lum_e16_lib as lib
import matplotlib.pyplot as plt
import numpy as np
import ROOT

DATE=datetime.datetime.now().strftime('%m%d%y')
INDIR=os.path.join(os.environ['D2PIDIR_EXP_E16'],'h10_2_dlum-byRun-e16_051717')
OUTDIR=os.path.join(os.environ['STUDY_LUM_E16_DATADIR'],'results_%s'%DATE)
if not os.path.exists(OUTDIR):
	os.makedirs(OUTDIR)

#! dlum => dlum.txt	
dlumtxt=os.path.join(os.path.join(OUTDIR,'dlum.txt'))
lib.make_dlum_txt(INDIR,dlumtxt)

#! Obtain data from dlum.txt in DataFrame format
lumd=lib.make_lum_df(dlumtxt)
#print lumd.head() # debug

#! Make some diagnostic plots
fig,axs=plt.subplots(figsize=(20,25),nrows=4,ncols=2)
plt.subplots_adjust(wspace=0.15)
#! q vs run
axs[0][0].scatter(lumd['run'],lumd['q'])
axs[0][0].set_ylabel("q")
#! nrm_Ntrg vs run
axs[1][0].scatter(lumd['run'],lumd['nrm_Ntrg'])
axs[1][0].set_ylabel("nrm_Ntrg")
#! nrm_N2pi vs run
axs[2][0].scatter(lumd['run'],lumd['nrm_N2pi'])
axs[2][0].set_ylabel("nrm_N2pi")
#! histogram for nrm_Ntrg
ret=axs[1][1].hist(lumd['nrm_Ntrg'],bins=100,range=(0,100000))
axs[1][1].set_xlabel("nrm_Ntrg")
#! histogram for nrm_N2pi
ret=axs[2][1].hist(lumd['nrm_N2pi'],bins=100,range=(0,40))
axs[2][1].set_xlabel("nrm_N2pi")
#! zoomed version of nrm_N2pi vs run to depict cuts 
axs[3][0].scatter(lumd['run'],lumd['nrm_N2pi'])
axs[3][0].set_ylim(15,40)
axs[3][0].set_yticks(np.arange(15, 40+1, 1.0))
axs[3][0].grid()
axs[3][0].set_ylabel("nrm_N2pi_zoomed")
#! Draw cut values
CUTS=[26,27]
CLRS=['g','r']
xmin,xmax=axs[3][0].get_xlim()[0],axs[3][0].get_xlim()[1]
for i,cut in enumerate(CUTS):
	axs[3][0].hlines(cut,xmin,xmax,color=CLRS[i],label="nrm_N2pi cut = %d"%cut)
axs[3][0].legend()
fig.savefig("%s/goodruns_cut_ana_plts.png"%(OUTDIR))

#! Draw ROOT-fitted version of zoomed version of nrm_N2pi vs run to depict cuts 
h=ROOT.TH1F("h_nrm_N2pi_zoomed","h_nrm_N2pi_zoomed",30,15,40)
h.SetXTitle("nrm_N2pi")
for x in list(lumd['nrm_N2pi']):
	h.Fill(x)
#h.Sumw2()
h.Fit("gaus")
ROOT.gStyle.SetOptFit(1111)
c=ROOT.TCanvas()
h.Draw("e")
#! Prepare and draw cut lines
CLRS_ROOT=[ROOT.gROOT.ProcessLine("kGreen"),ROOT.gROOT.ProcessLine("kRed")]
TCUTL=[0,0]
for i,cut in enumerate(CUTS):
	TCUTL[i]=ROOT.TLine(CUTS[i],h.GetMinimum(),CUTS[i],h.GetMaximum())
	TCUTL[i].SetLineColor(CLRS_ROOT[i])
	TCUTL[i].SetLineStyle(2)
	TCUTL[i].SetLineWidth(3)
#! Draw
for i,cut in enumerate(CUTS):
	TCUTL[i].Draw("sames")
c.SaveAs("%s/goodruns_cut_ana_nrm_N2pi_zoomed_fitted.png"%OUTDIR)

#! + Obtain luminosity based on CUTS on nrm_N2pi
#! + Save other relevant data based on this selection

#! Make a copy of runs that are cut away
XRUNL=[None for i in range(len(CUTS))]
for i in range(len(CUTS)):
	XRUNL[i]=list(lumd['run'][lumd['nrm_N2pi']<CUTS[i]])
#print XRUNL[0]
#print XRUNL[1]

#! Now using CUTRUNL make XQL, XN2PIL (X=cut) and their complements, GQL, GN2PIL (G=good)
XQL,    GQL   =[0 for i in range(len(CUTS))], [0 for i in range(len(CUTS))]
XN2PIL, GN2PIL=[0 for i in range(len(CUTS))], [0 for i in range(len(CUTS))]
#! Total X/GQ,X/GN2PI 
TXQ,    TGQ   =[None for i in range(len(CUTS))], [0 for i in range(len(CUTS))]
TXN2PI, TGN2PI=[None for i in range(len(CUTS))], [0 for i in range(len(CUTS))]
#print XQL
#print GQL
#sys.exit()

for i in range(len(CUTS)):
	#! X
	XQL[i]   =list(lumd['q'][lumd['run'].isin(XRUNL[i])])
	XN2PIL[i]=list(lumd['N2pi'][lumd['run'].isin(XRUNL[i])])
	#! X-total
	TXQ[i]   =sum(XQL[i])
	TXN2PI[i]=sum(XN2PIL[i])

	#! G
	GQL[i]   =list(lumd['q'][~lumd['run'].isin(XRUNL[i])])
	GN2PIL[i]=list(lumd['N2pi'][~lumd['run'].isin(XRUNL[i])])
	#! G-Total
	TGQ[i]   =sum(GQL[i])
	TGN2PI[i]=sum(GN2PIL[i])

#! Write and plot all this information
#! Write
fout=open(os.path.join('%s/lum_results.txt'%OUTDIR), 'w')
for i in range(len(CUTS)):
	fout.write("Information for nrm_N2pi cut=%d\n"%CUTS[i])
	fout.write("================================\n")
	fout.write("\n")

	#! First calculate print out some integrity related sanity-check data 
	#! which ensure that the number of runs=606
	txruns,tgruns=len(XQL[i]),len(GQL[i])
	truns=txruns+tgruns
	pxruns=(txruns/606)*100
	fout.write("*** Sanity-check based on number of total runs in $E16D_EI=%d ***\n"%606)
	fout.write("#cutruns=%d, #goodruns=%d => #total runs=%d\n"%(txruns,tgruns,tgruns))
	fout.write("Percentage cut runs=%.2f %%\n"%pxruns)
	if (truns!=606):fout.write("ERROR! Something is wrong. #totalruns != 606")
	#fout.write("******\n")
	fout.write("\n")

	#! Print cut run list
	txruns=len(XRUNL[i])
	fout.write("*** Cut Runs ***\n")
	fout.write("cut runs=%s\n"%XRUNL[i])
	#fout.write("******\n")
	fout.write("\n")

	#! Print out charge information (in mC, converted directly in print statement)
	txq,tgq=TXQ[i],TGQ[i] #! in microC here
	pxq=(txq/(txq+tgq))*100
	#! calculate p_diff_q based on KP (p156) ad BG (p66). Note PK = KP
	KP_Q=21.287*1e3  #! mC => microC
	BG_Q=20.6*1e3    #! mC => microC 
	p_diff_q_KP=((tgq-KP_Q)/KP_Q)*100
	p_diff_q_BG=((tgq-BG_Q)/BG_Q)*100
	fout.write("*** Faraday Cup Charge information***\n")
	fout.write("cutq=%.2f mC, goodq=%.2f mC\n"%(txq/1000,tgq/1000))
	fout.write("Percent cutq=%.2f %%\n"%pxq)
	fout.write("Percentage difference relative to KP-Q (%.3f mC)=%.2f %%\n"%(KP_Q/1000,p_diff_q_KP))
	fout.write("Percentage difference relative to BG-Q (%.3f mC)=%.2f %%\n"%(BG_Q/1000,p_diff_q_BG))
	#fout.write("******\n")
	fout.write("\n")

	#! Pring out luminosity information
	lum=lib.calc_lum(tgq*1e-6)
	#! calculate p_diff_lum based on KP ad BG (p156) ad BG (p66). Note PK = KP
	KP_LUM=lib.calc_lum(KP_Q*1e-6) # should be 28.13
	BG_LUM=lib.calc_lum(BG_Q*1e-6) # shoule be 27.22
	p_diff_lum_KP=((lum-KP_LUM)/KP_LUM)*100
	p_diff_lum_BG=((lum-BG_LUM)/BG_LUM)*100
	#! calculate p_diff_xsec also based on KP abd BG
	#! + assume that everything else in xsec forumula stays the same
	p_diff_xsec_KP=(((1/KP_LUM)-(1/lum))/(1/KP_LUM))*100
	p_diff_xsec_BG=(((1/BG_LUM)-(1/lum))/(1/BG_LUM))*100
	fout.write("*** Luminosity information ***\n")
	fout.write("Luminosity=%.2f fb^-1\n"%lum)
	fout.write("Percentage difference relative to KP-lum (%.2f fb^-1)=%.2f %%\n"%(KP_LUM,p_diff_lum_KP))
	fout.write("Percentage difference relative to BG-lum (%.2f fb^-1)=%.2f %%\n"%(BG_LUM,p_diff_lum_BG))
	fout.write("Percentage difference in xsec relative to xsec from KP-lum=%.2f %%\n"%p_diff_xsec_KP)
	fout.write("Percentage difference in xsec relative to xsec from BG-lum=%.2f %%\n"%p_diff_xsec_BG)	
	#fout.write("******\n")
	fout.write("\n")

	#! Print out n2pi information
	txn2pi,tgn2pi=TXN2PI[i],TGN2PI[i]
	pxn2pi=(txn2pi/(txn2pi+tgn2pi))*100
	fout.write("*** n2pi information *** \n")
	fout.write("cutn2pi=%d, goodn2pi=%d \n"%(txn2pi,tgn2pi))
	fout.write("Percent cutq=%.2f %%\n"%pxn2pi)
	#fout.write("******\n")
	fout.write("\n")

#! plot
fig,axs=plt.subplots(figsize=(15,10),nrows=2,ncols=2)
plt.subplots_adjust(wspace=0.15)
#! n2pi hist for each cut
for i in range(len(CUTS)):
	#! Draw cut n2pi
	r=axs[0][i].hist(XN2PIL[i],bins=250,range=(10,4000))
	axs[0][i].set_title("Cut n2pi distribution for nrm_n2pi cut=%d"%CUTS[i])
	#! Draw good n2pi
	r=axs[1][i].hist(GN2PIL[i],bins=250,range=(10,4000))
	axs[1][i].set_title("Good n2pi distribution for nrm_n2pi cut=%d"%CUTS[i])
fig.savefig("%s/lum_results.png"%(OUTDIR))