#!/usr/bin/python

from __future__ import division
import os,datetime,sys
import study_lum_e16_lib as lib
import matplotlib.pyplot as plt
import numpy as np
import ROOT
import time

'''
> study_lum_e16 tgt[=ptgt]<ptgt/etgt>
'''

#! Get input from user
TGT='ptgt'
if len(sys.argv)==2: #i.e. tgt info entered by user
	TGT=sys.argv[1]
print "TGT=%s"%TGT
if TGT!='ptgt' and TGT!='etgt':
	sys.exit('TGT=%s is not valid. TGT=ptgt or etgt only'%TGT)
#sys.exit()

DATE=datetime.datetime.now().strftime('%m%d%y')
if TGT=='ptgt':
	DLUMDIR=os.path.join(os.environ['D2PIDIR_EXP_E16'],'h10_2_dlum-byRun-e16_051717')
	D2PIDIR=os.path.join(os.environ['D2PIDIR_EXP_E16'],'h10-skim-SS_2_d2piR-byRun_052217')
	OUTDIR=os.path.join(os.environ['STUDY_LUM_E16_DATADIR'],'results_ptgt_%s'%DATE)
if TGT=='etgt':
	DLUMDIR=os.path.join(os.environ['D2PIDIR_EXP_E16'],'h10-etgt_2_dlum-byRun_052517')
	D2PIDIR=os.path.join(os.environ['D2PIDIR_EXP_E16'],'h10-etgt_2_d2piR-byRun_052517')
	OUTDIR=os.path.join(os.environ['STUDY_LUM_E16_DATADIR'],'results_etgt_%s'%DATE)

print "Going to analyze and calculate luminosity for TGT=%s"%TGT
print "DLUMDIR=",DLUMDIR
print "D2PIDIR=",D2PIDIR
print "OUTDIR=",OUTDIR
time.sleep(5)

if not os.path.exists(OUTDIR):
	os.makedirs(OUTDIR)

#! dlum => dlum.txt	
dlumtxt=os.path.join(os.path.join(OUTDIR,'dlum.txt'))
lib.make_dlum_txt(DLUMDIR,D2PIDIR,dlumtxt)

#! Obtain data from dlum.txt in DataFrame format
lumd=lib.make_lum_df(dlumtxt)
#print lumd.head() # debug

#! begin etgt analysis
if TGT=='etgt': #! Then do simplified anaylysis and exit
	#! Make some diagnostic plots
	fig,axs=plt.subplots(figsize=(20,25),nrows=3,ncols=2)
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
	#! Draw cut values
	#! Note compared to ptgt, there is only cut value, but with a min and max limit
	CUTS=[[0.35,0.45]]
	CLRS=['r']
	xmin,xmax=axs[2][0].get_xlim()[0],axs[2][0].get_xlim()[1]
	for i in range(len(CUTS)):
		cutmin,cutmax=CUTS[0][0],CUTS[0][1]
		axs[2][0].hlines(cutmin,xmin,xmax,color=CLRS[i],label="nrm_N2pi cut min = %.3f"%cutmin)
		axs[2][0].hlines(cutmax,xmin,xmax,color=CLRS[i],label="nrm_N2pi cut max = %.3f"%cutmax)
	axs[2][0].legend()
	fig.savefig("%s/goodruns_cut_ana_plts.png"%(OUTDIR))

	#! Make a copy of runs that are within cuts i.e. good runs(=G)
	GRUNL=[None for i in range(len(CUTS))]
	for i in range(len(CUTS)):
		cutmin,cutmax=CUTS[0][0],CUTS[0][1]
		GRUNL[i]=list(lumd['run'][(lumd['nrm_N2pi']>cutmin) & (lumd['nrm_N2pi']<cutmax)])

	#! print out relevant information
	fout=open(os.path.join('%s/lum_results.txt'%OUTDIR), 'w')
	for i in range(len(CUTS)):
		cutmin,cutmax=CUTS[0][0],CUTS[0][1]
		fout.write("Information for nrm_N2pi cutmin,cutmax=%f,%f\n"%(cutmin,cutmax))
		fout.write("================================\n")
		fout.write("\n")
	
		#! First calculate print out some integrity related sanity-check data 
		#! which ensure that the number of runs=606
		truns=len(lumd)
		fout.write("*** Sanity-check based on total number of selected E16 etgt runs for analysis =%d ***\n"%8)
		fout.write("#total etgt runs selected for analysis=%d\n"%(truns))
		if (truns!=8):fout.write("ERROR! Something is wrong. #total selected etgt runs != 8")
		#fout.write("******\n")
		fout.write("\n")

		#! Print etgt run list
		runs=list(lumd['run'][lumd['run'].isin(GRUNL[i])])
		fout.write("*** etgt runs used for analysis ***\n")
		fout.write("selected etgt runs=%s\n"%runs)
		#fout.write("******\n")
		fout.write("\n")

		#! Print out charge information (in mC, converted directly in print statement)
		fout.write("*** Q from selected etgt runs ***\n")
		q=sum(lumd['q'][lumd['run'].isin(GRUNL[i])]) #! in microC here
		fout.write("Q=%.2f mC"%(q/1000))
		#fout.write("******\n")
		fout.write("\n")
	sys.exit()
#! begin etgt analysis

#! begin ptgt analysis 
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
#! Draw cuts
#! [07-10-17] 
#! + Updated two low bound cuts of 26 and 27 (lse and tgt respectively)
#!   to a single double sided cut of [26,34] as per fit to Gaussuan norm_2pi distribution (see below)
#!   and obtaining these as rounded values from (mu-3*sg,mu+3*sg): 
#!   [29.75-3*1.398,29.75+3*1.398] = [25.56,33.94] => (rounded) [26,34]
CUTS=[[26,34]] #! CUTS=[26,27]
CLRS=['g','r']
xmin,xmax=axs[3][0].get_xlim()[0],axs[3][0].get_xlim()[1]
for i,cut in enumerate(CUTS):
	cutmin,cutmax=CUTS[i][0],CUTS[i][1]
	axs[2][0].hlines(cutmin,xmin,xmax,color=CLRS[i],label="nrm_N2pi cut min = %.3f"%cutmin)
	axs[2][0].hlines(cutmax,xmin,xmax,color=CLRS[i],label="nrm_N2pi cut max = %.3f"%cutmax)
axs[3][0].legend()
fig.savefig("%s/goodruns_cut_ana_plts.png"%(OUTDIR))

#! Draw ROOT-fitted version of zoomed version of nrm_N2pi vs run to depict cuts 
h=ROOT.TH1F("h_nrm_N2pi_zoomed","h_nrm_N2pi_zoomed",30,15,40)
h.SetXTitle("nrm_N2pi")
for x in list(lumd['nrm_N2pi']):
	h.Fill(x)
#h.Sumw2()
h.Fit("gaus","","",27,33)
ROOT.gStyle.SetOptFit(1111)
c=ROOT.TCanvas()
h.Draw("e")
#! Prepare and draw cut lines
CLRS_ROOT=[ROOT.gROOT.ProcessLine("kGreen"),ROOT.gROOT.ProcessLine("kRed")]
TCUTL=[0 for i in range(len(CUTS))]
TCUTH=[0 for i in range(len(CUTS))]
for i,cut in enumerate(CUTS):
	cutmin,cutmax=CUTS[i][0],CUTS[i][1]
	TCUTL[i]=ROOT.TLine(cutmin,h.GetMinimum(),cutmin,h.GetMaximum())
	TCUTH[i]=ROOT.TLine(cutmax,h.GetMinimum(),cutmax,h.GetMaximum())
	for line in [TCUTL[i],TCUTH[i]]:
		line.SetLineColor(CLRS_ROOT[i])
		line.SetLineStyle(2)
		line.SetLineWidth(3)
#! Draw
for i,cut in enumerate(CUTS):
	for line in [TCUTL[i],TCUTH[i]]:
		line.Draw("sames")
c.SaveAs("%s/goodruns_cut_ana_nrm_N2pi_zoomed_fitted.png"%OUTDIR)

#! + Obtain luminosity based on CUTS on nrm_N2pi
#! + Save other relevant data based on this selection

#! Make a copy of runs that are cut away
XRUNL=[None for i in range(len(CUTS))]
for i in range(len(CUTS)):
	cutmin,cutmax=CUTS[0][0],CUTS[0][1]
	print cutmin,cutmax
	XRUNL[i]=list(lumd['run'][(lumd['nrm_N2pi']<cutmin) | (lumd['nrm_N2pi']>cutmax)])
	#XRUNL[i]=list(lumd['run'][lumd['nrm_N2pi']<CUTS[i]])
#print XRUNL[0]
#print XRUNL[1]
#sys.exit()

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
	cutmin,cutmax=CUTS[0][0],CUTS[0][1]
	fout.write("Information for nrm_N2pi cut=(%d,%d)\n"%(cutmin,cutmax))
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
	cutmin,cutmax=CUTS[0][0],CUTS[0][1]
	#! Draw cut n2pi
	r=axs[0][i].hist(XN2PIL[i],bins=250,range=(10,4000))
	axs[0][i].set_title("Cut n2pi distribution for nrm_n2pi cut=(%d,%d)"%(cutmin,cutmax))
	#! Draw good n2pi
	r=axs[1][i].hist(GN2PIL[i],bins=250,range=(10,4000))
	axs[1][i].set_title("Good n2pi distribution for nrm_n2pi cut=(%d,%d)"%(cutmin,cutmax))
fig.savefig("%s/lum_results.png"%(OUTDIR))