#!/usr/bin/python
from __future__ import division

import os,sys,datetime
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict
import itertools

import numpy as np

import math

#atlib
import atlib as atlib

'''
	+ Note on bg_display_modes: full and simple (default)
	The background obtained from sim-3pi-PS can be displayed (and justified )in two ways:
	(Note that this relates *only to the visual display* of the background i.e. its estimate
	is obtained using just one way (normalizing hSR3PI to hdiff), but is displayed in two ways, one of which (simple) is easier
	to interpret)
	1. full: This display is true to the exact process used to obtain the bg estimate in that hSR3PIn is drawn
	         directly in comparison with hdiff, which is the histogram that is used to normalized it. However, hdiff has a wiggle-
	         structure around 0 because hER is not equal to hSR, which is complicated to explain and it detracts
	         from the main point of this study.
	2. simple: To avoid potential detraction from the main result of the study, a simple display can be made
	           where hSR3PIn+hSR is plotted and directly compared with hER:
	             hSR3PIn+hSR=hER
	           Note that this is consistent with 'full' because:
	             hSR3PIn=hER-hSR
'''

USAGE='study_evtsel_BG.py debug=False test_fit_only=False bg_display_mode=simple <simple/full> norm_opt_SR_to_ER=integral <integral/maximum>'

#! Get input data from user
DBG=False
if len(sys.argv)>1: #i.e. debug entered by user
	if    sys.argv[1]=="True":  DBG=True
	elif  sys.argv[1]=="False": DBG=False
	else: sys.exit('DBG=%s is not valid. usage: %s'%(sys.argv[1],USAGE))

TEST_FIT_ONLY=False
if len(sys.argv)>2: #i.e. test_fit_only entered by user
	if    sys.argv[2]=="True":  TEST_FIT_ONLY=True
	elif  sys.argv[2]=="False": TEST_FIT_ONLY=False
	else: sys.exit('TEST_FIT_ONLY=%s is not valid. usage: %s'%(sys.argv[2],USAGE))

BG_DISPLAY_MODE='simple'
if len(sys.argv)>3: #i.e. bg_display_mode entered by user
	if    sys.argv[3]=="simple":  BG_DISPLAY_MODE=sys.argv[3]
	elif  sys.argv[3]=="full":    BG_DISPLAY_MODE=sys.argv[3]
	else: sys.exit('BG_DISPLAY_MODE=%s is not valid. usage: %s'%(sys.argv[3],USAGE))

NORM_OPT_SR_TO_ER='integral'
if len(sys.argv)>4: #i.e. norm_opt_SR_to_ER entered by user
	if    sys.argv[4]=="integral":   NORM_OPT_SR_TO_ER=sys.argv[4]
	elif  sys.argv[4]=="maximum":    NORM_OPT_SR_TO_ER=sys.argv[4]
	else: sys.exit('NORM_OPT_SR_TO_ER=%s is not valid. usage: %s'%(sys.argv[4],USAGE))

print "DBG=",DBG
print "TEST_FIT_ONLY=",TEST_FIT_ONLY
print "BG_DISPLAY_MODE=",BG_DISPLAY_MODE
print "NORM_OPT_SR_TO_ER=",NORM_OPT_SR_TO_ER
#sys.exit()

#! *** Prepare input data structure: FIN/T[dtyp] ***
NDTYP=3
ER,SR,SR3PI=range(NDTYP)
DTYP_NAME=["ER","SR","SR3PI"]

FIN=[0 for i in range(NDTYP)]
T  =[0 for i in range(NDTYP)]
#print "FIN=",FIN
#print "T=",T
#! ***

#! *** Prepare OUTDIR***
DATE=datetime.datetime.now().strftime('%m%d%y')
OUTDIR=os.path.join(os.environ['STUDY_EVTSEL_BG_DATADIR'],"evtsel_BG_bg-%s_norm-%s_%s"%(BG_DISPLAY_MODE,NORM_OPT_SR_TO_ER,DATE))
if DBG==True: OUTDIR=os.path.join(os.environ['STUDY_EVTSEL_BG_DATADIR'],"evtsel_BG_dbg_bg-%s_norm-%s_%s"%(BG_DISPLAY_MODE,NORM_OPT_SR_TO_ER,DATE))
print "OUTDIR=",OUTDIR
if not os.path.exists(OUTDIR):
	os.makedirs(OUTDIR)
#sys.exit()
#! ***

#! *** Prepare NENTRIES ***
NENTRIES=1000000000
#if DBG==True: NENTRIES=1000

#! *** Now prepare INDIR and input data structures: FIN/T[dtyp] ***
INDIR=[[] for i in range(NDTYP)]
INDIR[ER]   =os.path.join(os.environ['D2PIDIR_EXP_E16'],"data_pcorr_011116")#! debug /bk
INDIR[SR]   =os.path.join(os.environ['D2PIDIR_SIM_E16'],"data_pcorr_011116")#! debug /bl
INDIR[SR3PI]=os.path.join(os.environ['D2PIDIR_SIM_E16'],"data_evtsel_BG_020117")

print "INDIR[ER]=",INDIR[ER]
print "INDIR[SR]=",INDIR[SR]
print "INDIR[SR3PI]=",INDIR[SR3PI]
#sys.exit()

FIN[ER]   =ROOT.TFile("%s/d2piR_wpcorr.root"%INDIR[ER]) #!npcorr
FIN[SR]   =ROOT.TFile("%s/d2piR_npcorr.root"%INDIR[SR])
FIN[SR3PI]=ROOT.TFile("%s/d3piPS.root"%INDIR[SR3PI])
for idtyp in range(NDTYP):
	T[idtyp]=FIN[idtyp].Get("d2piR/tR")
#sys.exit()

#! *** Now make plots ***

#! sector information
NSCTR=6

#! *** Now make MM plots ***
WMIN=1.400
WMAX=2.125
WBINW=0.025
Q2BIN_LEL=[2.00,2.40,3.00,3.50,4.20,2.00]#[1.50,2.00,2.40,3.00,3.50,4.20] #debug
Q2BIN_UEL=[2.40,3.00,3.50,4.20,5.00,5.00]#[2.00,2.40,3.00,3.50,4.20,5.00] #debug
# Q2BIN_LEL=[2.00]
# Q2BIN_UEL=[5.00]
if DBG:
	Q2BIN_LEL=[2.00]
	Q2BIN_UEL=[5.00]
WBIN_LEL=np.arange(WMIN,WMAX,WBINW)
#! [05-10-17] Remove 2.125 if it occures in WBIN_LEL
#! + Sometimes because of using decimal numbers with np.arange, 2.125 is not omitted 
#!   and therefore it has to be done manually
del_val=False
idel_val=-9999
for i,x in enumerate(WBIN_LEL):
	if np.isclose(x,2.125): 
		del_val=True
		idel_val=i
#print len(WBIN_LEL)
if del_val: 
    print "In WBIN_LEL going to delete",WBIN_LEL[idel_val],"at index",idel_val
WBIN_LEL=np.delete(WBIN_LEL,[idel_val])
#print WBIN_LEL
#sys.exit()
#print Q2BIN_LEL
#print WBIN_LEL

#! MM cuts lines: EI and AT
#! EI
EI_MM2L,EI_MM2H=-0.04,0.06
EI_MML, EI_MMH=-math.sqrt(math.fabs(EI_MM2L)),math.sqrt(EI_MM2H)
MM2_T2_CUTL_EI=ROOT.TLine(EI_MM2L,0,EI_MM2L,0);
MM2_T2_CUTH_EI=ROOT.TLine(EI_MM2H,0,EI_MM2H,0);
MM_T2_CUTL_EI=ROOT.TLine(EI_MML,0,EI_MML,0);
MM_T2_CUTH_EI=ROOT.TLine(EI_MMH,0,EI_MMH,0);
for l in [MM2_T2_CUTL_EI,MM2_T2_CUTH_EI,MM_T2_CUTL_EI,MM_T2_CUTH_EI]:
	l.SetLineColor(ROOT.gROOT.ProcessLine("kMagenta"))
	l.SetLineStyle(2)
	l.SetLineWidth(3)

#! TLine depicting pion mass debug
MASS_PION=0.13957018;
MASS_PION_LINE=ROOT.TLine(MASS_PION,0,MASS_PION,0)
MASS_PION_LINE.SetLineColor(ROOT.gROOT.ProcessLine("kGreen"))
MASS_PION_LINE.SetLineWidth(3)

#! DRAW_CMD
NPLTS=2
MM,MM2=range(NPLTS)
PLT_NAME=['MM','MM2']
PLT_TITLE=['MM','MM^{2}']
PLT_UNIT=['GeV','GeV^{2}']

DRAW_CMD=[0 for i in range(NPLTS)]
DRAW_CMD[MM] ="mmppip>>hmmcmd(200,-0.5,1)"
DRAW_CMD[MM2]="mm2ppip>>hmmcmd(100,-0.2,0.2)"

#! Normalization-range limits in which hSR3PI will be normed to hdiff
NL=[0 for i in range(NPLTS)]
#NL[MM]=[0.278,0.417]
NL[MM]=[0.300,0.400]
NL[MM2]=[NL[MM][0]**2,NL[MM][1]**2]

NL1=[0 for i in range(NPLTS)]
#NL1[MM]=[0.278,0.417]
NL1[MM]=[0.300,0.400]
NL1[MM2]=[NL1[MM][0]**2,NL1[MM][1]**2]

#! TLine showing fit limits
NL_LINE_MIN=[0 for i in range(NPLTS)]
NL_LINE_MAX=[0 for i in range(NPLTS)]
NL_LINE_MIN[MM]=ROOT.TLine( NL[MM][0], 0, NL[MM][0], 0)
NL_LINE_MAX[MM]=ROOT.TLine( NL[MM][1], 0, NL[MM][1], 0)
NL_LINE_MIN[MM2]=ROOT.TLine( NL[MM2][0], 0, NL[MM2][0], 0)
NL_LINE_MAX[MM2]=ROOT.TLine( NL[MM2][1], 0, NL[MM2][1], 0)
for line in [NL_LINE_MIN[MM],NL_LINE_MAX[MM],NL_LINE_MIN[MM2],NL_LINE_MAX[MM2]]:
	line.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
	line.SetLineWidth(3)

NL1_LINE_MIN=[0 for i in range(NPLTS)]
NL1_LINE_MAX=[0 for i in range(NPLTS)]
NL1_LINE_MIN[MM]=ROOT.TLine( NL1[MM][0], 0, NL1[MM][0], 0)
NL1_LINE_MAX[MM]=ROOT.TLine( NL1[MM][1], 0, NL1[MM][1], 0)
NL1_LINE_MIN[MM2]=ROOT.TLine( NL1[MM2][0], 0, NL1[MM2][0], 0)
NL1_LINE_MAX[MM2]=ROOT.TLine( NL1[MM2][1], 0, NL1[MM2][1], 0)
for line in [NL1_LINE_MIN[MM],NL1_LINE_MAX[MM],NL1_LINE_MIN[MM2],NL1_LINE_MAX[MM2]]:
	line.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
	line.SetLineWidth(3)


#! TLine showing 0 of y-axis (after making y-axis go below to see fluctuations in hdiff=exp-sim)
ZERO_YAXIS=[0 for i in range(NPLTS)]
ZERO_YAXIS[MM]=ROOT.TLine( -0.5, 0, 1.0, 0)
ZERO_YAXIS[MM2]=ROOT.TLine(-0.2, 0, 0.2, 0)

#! Create hmm[dtyp][q2][w][plt]
hmm=[[[[0 for l in range(NPLTS)] for k in range(len(WBIN_LEL))] for j in range(len(Q2BIN_LEL))] for i in range(NDTYP)]
#print "hmm=",hmm
#sys.exit()

#! Create domain for making [ER,SR,SR3PI]*[MM,MM2] plots within a q2-w bin
d=[[ER,SR,SR3PI],[MM,MM2]]
D=list(itertools.product(*d))

#! *** Setup some analysis constants ***

#! Setup W bin up to which sim3pi cannot be normalized
WBINS_NORM3PI_NOT_POSSIBLE=[[1.400,1.425],[1.425,1.450],[1.450,1.475],[1.475,1.500],[1.500,1.525],
                            [1.525,1.550],[1.550,1.575],[1.575,1.600],[1.600,1.625],[1.625,1.650]]

#! Setup W bin number up to which fit to bg should be ignored
#! + only for W bin >=[1.650,1.675) (11th W bin, idx=10) 
FIT_IWBIN_GT_ET=10
#! + All bins
#FIT_IWBIN_GT_ET=0
#! ***

#! define some functions needed
def norm_hSR3PI_based_on_bg_fit(hER,hSR3PI,bg,cut_min,cut_max):
	'''
	+ For cases where hSR3PI cannot be normalized to hER in the leading-edge,
	  it is normalized based on fit to bg-data where the bg can be estimated because this normalization is possible
	+ From the fit the bg in the area where this normalization fails can be estimated and based on that the normalization
	  factor for hSR3PI can be extracted based on the fact that:

	    + (int_cut_rgn(hSR3PI)*sf)/int_cut_rgn(hER) * 100 = bg, where
	      + bg is the bg(%) estimated from the fit
	      + sf = normalization factor (in technical terms, a scale factor)

	    Based on the above, the formula to obtain the sf is:
	    + sf= (bg*int_cut_rgn(hER))/(int_cut_rgn(hSR3PI)*100)
	'''
	#! obtain sf
	bin_min=hER.FindBin(cut_min)
	bin_max=hER.FindBin(cut_max)
	if bg!=0:
		sf=(bg*hER.Integral(bin_min,bin_max))/(hSR3PI.Integral(bin_min,bin_max)*100)
	else:
		sf=0
	print "norm_hSR3PI_based_on_bg_fit():sf=%.3f"%sf
	#! Return normalized copy of hSR3PI
	hSR3PIn=hSR3PI.Clone("hSR3PIn")
	hSR3PIn.Sumw2()
	hSR3PIn.Scale(sf)
	return hSR3PIn

def plot_fit_and_write_bg(BG):
	#! Some related aesthetics
	#! LABEL[plt,cut]
	LABEL=[0 for i in range(NPLTS)] 
	LABEL[MM]="%.2f GeV < MM < %.2f GeV"%(EI_MML,EI_MMH)
	LABEL[MM2]="%.2f GeV < " r"$MM^{2}$" " < %.2f " r"$GeV^{2}$"%(EI_MM2L,EI_MM2H)
	for iq,q2bin_le in enumerate(Q2BIN_LEL):
		if DBG==True and iq>0: continue #! debug
		q2min=q2bin_le
		q2max=Q2BIN_UEL[iq]
		outdir_q2w="%s/with_BG_ana/%.2f-%.2f"%(OUTDIR,q2min,q2max)
		if not os.path.exists(outdir_q2w):
			os.makedirs(outdir_q2w)
		#! set q2bin label
		q2binl='[%.2f,%.2f)'%(q2min,q2max)
		for iplt in range(NPLTS):
			outdir_q2w_plt="%s/%s"%(outdir_q2w,PLT_NAME[iplt])
			if not os.path.exists(outdir_q2w_plt):
				os.makedirs(outdir_q2w_plt)
			#! structure to consolidate W data for q2,plt
			bg=[]
			wbinn=[] #! wbin number
			wbinl=[] #! wbin label
			for iw,wbin_le in enumerate(WBIN_LEL):
				#! DBGW (3)
				if DBG==True and (iw!=0 and iw!=15 and iw!=28): continue #! debug #!(iw!=0 and iw!=15 and iw!=28)
				#if DBG==True and iw+1 > 10: continue
				#if DBG==True and (iw+1!=2 and iw+1!=3 and iw+1!=4 and iw+1!=5): continue
				wmin=wbin_le
				wmax=wbin_le+WBINW
				bg.append(BG[iq][iw][iplt])
				wbinn.append(iw+1)
				wbinl.append('[%.3f,%.3f)'%(wmin,wmax))

			#! plot
			fig=plt.figure()
			ax=plt.subplot(111)
			#! To give space below x axis for rotated W bin labels
			fig.subplots_adjust(bottom=0.25)
			ax.scatter(wbinn,bg,c='m', label=LABEL[iplt])
			#! Fix x axis
			# shift half a step to the left
			xmin=(3*wbinn[0]-wbinn[1])/2.
			# shift half a step to the right
			xmax=(3*wbinn[-1]-wbinn[-2])/2.
			ax.set_xlim(xmin,xmax)
			ax.set_xticks(wbinn)
			ax.set_xticklabels(wbinl,rotation='vertical')
			ax.set_title(r"Background for" r" $Q^{2}$" "=[%.2f,%.2f) " r"$GeV^{2}$"%(q2min,q2max))
			ax.set_xlabel("W (GeV)")
			ax.set_ylabel("Background (%)")
			#ax.legend(loc='upper right',prop={'size':8})

			#! fit with 3 order poly
			#! + constrain to go to 0 by using setting BG for 1st W bin to 0
			pars=np.polyfit([wbinn[0]]+wbinn[FIT_IWBIN_GT_ET:],[0]+bg[FIT_IWBIN_GT_ET:],3)
			
			#! Get fitted bg
			#! note highest order power coefficient is returned first: 
			#!  + https://docs.scipy.org/doc/numpy/reference/generated/numpy.polyfit.html

			#bg_fitl=[pars[8] + pars[7]*x**1 + pars[6]*x**2 + pars[5]*x**3 + pars[4]*x**4+ pars[3]*x**5+ pars[2]*x**6 + pars[1]*x**7 + pars[0]*x**8 for x in wbinn]
			bg_fitl=[pars[3] + pars[2]*x**1 + pars[1]*x**2 + pars[0]*x**3 for x in wbinn]
			
			#! Refit if any value from bg_fitl<0, which should only happen in the extrapolated region 
			#! and therein most probably in the first 3 bins
			#! Note some refitting settings to "keep things simple":
			#! 1. The precision of comparison to 0 is at 0.01 level i.e. negative values of greater order are treated as 0 because:
			#!   + Precision for background at 0.01% level is enough
			#!   + Additionally, fitting can get complicated in trying to get values to >0 with such levels of precision
			#! 2. Note that refits are limited upto a maximum of 10 times to avoid infinite loops
			refit=True
			nitr_refit=0
			print "refitting for iplt==",PLT_NAME[iplt]
			while (refit==True and nitr_refit<10):
				vals_lt_0=[]
				for i,x in enumerate(wbinn):
					if bg_fitl[i]<0 and not np.isclose(math.fabs(bg_fitl[i]),0,atol=1e-2): 
						vals_lt_0.append(bg_fitl[i])
				if len(vals_lt_0)>0:
					#print len(vals_lt_0)
					print vals_lt_0
					val=math.fabs(min(vals_lt_0))
					cnstrnt_rft=val
					print "Adjusting contraint for fit to go to %f instead of 0 in the first W bin"%cnstrnt_rft
					pars=np.polyfit([wbinn[0]]+wbinn[FIT_IWBIN_GT_ET:],[cnstrnt_rft]+bg[FIT_IWBIN_GT_ET:],3)
					bg_fitl=[pars[3] + pars[2]*x**1 + pars[1]*x**2 + pars[0]*x**3 for x in wbinn]
					nitr_refit+=1
				else:
					refit=False
			print "Number of refit iterations for %s=%d"%(PLT_NAME[iplt],nitr_refit)
			print "bg_fitl for %s:"%(PLT_NAME[iplt])
			print bg_fitl
			
			#! plot
			ax.plot(wbinn,bg_fitl,'-',c='m')
			#! draw legend
			ax.legend(loc='upper right',prop={'size':8})
			#! save figure
			fig.savefig('%s/bg.png'%(outdir_q2w_plt))
			fig.savefig('%s/bg.pdf'%(outdir_q2w_plt))
			fig.savefig('%s/bg.eps'%(outdir_q2w_plt))

			#! Save fit data in text file
			#! Note that if bg_fit<0, it is set to 0
			#! + This is true after re-refitting where attempt is made, within "keep things simple", to not let fit
			#!   go below 0. Read refitting procedure comments for more details
			fout=open("%s/bg_from_fit.txt"%(outdir_q2w_plt),'w')
			for i,x in enumerate(wbinn):
				wbin=wbinl[i]
				if bg_fitl[i]<0: bgf=0
				else:            bgf=bg_fitl[i]
				fout.write("%d %s %.2f\n"%(x,wbin,bgf))
			fout.close()

			#! Save bg from  data points in text file
			fout=open("%s/bg_from_data.txt"%(outdir_q2w_plt),'w')
			for i,x in enumerate(wbinn):
				wbin=wbinl[i]
				fout.write("%d %s %.2f\n"%(x,wbin,bg[i]))
			fout.close()

def plot_without_BG_analysis():
	'''
	plot MM and MM2 without BG analysis
	'''

	for iq,q2bin_le in enumerate(Q2BIN_LEL):
		if DBG==True and iq>0: continue #! debug
		q2min=q2bin_le
		q2max=Q2BIN_UEL[iq]
		outdir_q2w="%s/no_BG_ana/%.2f-%.2f"%(OUTDIR,q2min,q2max)
		if not os.path.exists(outdir_q2w):
			os.makedirs(outdir_q2w)
		for iw,wbin_le in enumerate(WBIN_LEL):
			#! DBGW (3)
			if DBG==True and (iw!=0 and iw!=15 and iw!=28): continue #! debug #!(iw!=0 and iw!=15 and iw!=28)
			#if DBG==True and iw+1 > 10: continue
			#if DBG==True and (iw+1!=2 and iw+1!=3 and iw+1!=4 and iw+1!=5): continue
			wmin=wbin_le
			wmax=wbin_le+WBINW
		
			print "Going to plot hmm and hmm2 for %.2f-%.2f_%.3f-%.3f"%(q2min,q2max,wmin,wmax)

			#! Create cmm[plt]
			cmm=[0 for i in range(NPLTS)]
			#! Now plot
			for iplt in range(NPLTS):
				outdir_q2w_plt="%s/%s"%(outdir_q2w,PLT_NAME[iplt])
				if not os.path.exists(outdir_q2w_plt):
					os.makedirs(outdir_q2w_plt)
				cname="c%s_%.3f-%.3f"%(PLT_NAME[iplt],wmin,wmax)
				cmm[iplt]=ROOT.TCanvas(cname,cname)
				#l=ROOT.TLegend(0.1,0.3,0.3,0.4)#,"","NDC");
				#! hists aesthetics	
				for idtyp in range(NDTYP):
					hmm[idtyp][iq][iw][iplt].SetLineColor(CLRS_DTYP[idtyp])
					hmm[idtyp][iq][iw][iplt].SetMarkerColor(CLRS_DTYP[idtyp])
				#! Draw hists
				#! First directly get copy of hists so as to not use tedious indices
				hER=hmm[ER][iq][iw][iplt]
				hSR=hmm[SR][iq][iw][iplt]
				hSR3PI=hmm[SR3PI][iq][iw][iplt]
				hER.Sumw2()
				hSR.Sumw2()
				hSR3PI.Sumw2()
				#! First draw hER
				hER.SetXTitle("%s [%s]"%(PLT_TITLE[iplt], PLT_UNIT[iplt]))
				cmm[iplt].SetLeftMargin(0.20)
				hER.GetYaxis().SetTitleOffset(1.5)
				hER.SetYTitle("N_{entries}")
				hER.Draw()
				ZERO_YAXIS[iplt].Draw("same")
				#! scale SR and then draw
				max_ER=hER.GetMaximum()
				if max_ER==0:max_ER=1
				max_SR=hSR.GetMaximum()
				if max_SR==0:max_SR=1
				scl_fctr_SR=max_ER/max_SR
				if scl_fctr_SR==0:scl_fctr_SR=1
				hSR.Scale(scl_fctr_SR)
				hSR.Draw("same")

				#! Create and add entries: hists and cut-lines
				#! legend
				l=ROOT.TLegend(0.66,0.72,0.9,0.9)#,"","NDC");	
				for hist,label in zip([hER,hSR],["exp","sim"]):
					l.AddEntry(hist,label,"lp")
				#! Draw cut lines and add them to legend
				if iplt==MM:
					#! Draw cut lines
					MM_T2_CUTL_EI.SetY1(hER.GetMinimum())
					MM_T2_CUTL_EI.SetY2(hER.GetMaximum())
					MM_T2_CUTH_EI.SetY1(hER.GetMinimum())
					MM_T2_CUTH_EI.SetY2(hER.GetMaximum())
					MM_T2_CUTL_EI.Draw("same")
					MM_T2_CUTH_EI.Draw("same")
					#! add to legend
					l.AddEntry(MM_T2_CUTL_EI,"%.2f GeV < MM < %.2f GeV"%(EI_MML,EI_MMH),"l")
				elif iplt==MM2:
					#! EI
					MM2_T2_CUTL_EI.SetY1(hER.GetMinimum())
					MM2_T2_CUTL_EI.SetY2(hER.GetMaximum())
					MM2_T2_CUTH_EI.SetY1(hER.GetMinimum())
					MM2_T2_CUTH_EI.SetY2(hER.GetMaximum())
					MM2_T2_CUTL_EI.Draw("same")
					MM2_T2_CUTH_EI.Draw("same")
					#! add to legend
					l.AddEntry(MM2_T2_CUTL_EI,"%.2f GeV^{2} < MM^{2} < %.2f GeV^{2}"%(EI_MM2L,EI_MM2H),"l")
				#! Draw legend
				l.Draw()
				cmm[iplt].SaveAs("%s/%s.png"%(outdir_q2w_plt,cname))
				cmm[iplt].SaveAs("%s/%s.pdf"%(outdir_q2w_plt,cname))

def test_fit_bg():
	'''
	Use data bg_from_data.txt (not bg_from_fit.txt) to draw the data points and try fitting them
	'''

	#! Some related aesthetics
	#! LABEL[plt,cut]
	LABEL=[0 for i in range(NPLTS)] 
	LABEL[MM]="%.2f GeV < MM < %.2f GeV"%(EI_MML,EI_MMH)
	LABEL[MM2]="%.2f GeV < " r"$MM^{2}$" " < %.2f " r"$GeV^{2}$"%(EI_MM2L,EI_MM2H)
	for iq,q2bin_le in enumerate(Q2BIN_LEL):
		if DBG==True and iq>0: continue #! debug
		q2min=q2bin_le
		q2max=Q2BIN_UEL[iq]
		outdir_q2w="%s/%.2f-%.2f"%(OUTDIR,q2min,q2max)
		if not os.path.exists(outdir_q2w):
			continue
		#! set q2bin label
		q2binl='[%.2f,%.2f)'%(q2min,q2max)
		for iplt in range(NPLTS):
			outdir_q2w_plt="%s/%s"%(outdir_q2w,PLT_NAME[iplt])
			if not os.path.exists(outdir_q2w_plt):
				continue
			#! file to read data from
			f=open("%s/bg_from_data.txt"%(outdir_q2w_plt),'r')
			#! structure to consolidate W data for q2,plt
			bg=[]
			wbinn=[] #! wbin number
			wbinl=[] #! wbin label
			for line in f:
				arr=line.rstrip().split(" ")
				wbinn.append(int(arr[0]))
				wbinl.append(arr[1])
				bg.append(float(arr[2]))
			
			#! plot
			fig=plt.figure()
			ax=plt.subplot(111)
			#! To give space below x axis for rotated W bin labels
			fig.subplots_adjust(bottom=0.25)
			ax.scatter(wbinn,bg,c='m', label=LABEL[iplt])
			#! Fix x axis
			# shift half a step to the left
			xmin=(3*wbinn[0]-wbinn[1])/2.
			# shift half a step to the right
			xmax=(3*wbinn[-1]-wbinn[-2])/2.
			ax.set_xlim(xmin,xmax)
			ax.set_xticks(wbinn)
			ax.set_xticklabels(wbinl,rotation='vertical')
			ax.set_title(r"Background for" r" $Q^{2}$" "=[%.2f,%.2f) " r"$GeV^{2}$"%(q2min,q2max))
			ax.set_xlabel("W (GeV)")
			ax.set_ylabel("Background (%)")
			#ax.legend(loc='upper right',prop={'size':8})

			#! Fit condition 1. Fit under condition that there should be no increase in BG for low W
			#! + Therefore fit only in range where BG is decreasing
			#! + constrain to go to 0 by using setting BG for 1st W bin to 0
			pars=np.polyfit([wbinn[0]]+wbinn[FIT_IWBIN_GT_ET:],[0]+bg[FIT_IWBIN_GT_ET:],3)
			#pars=np.polyfit([wbinn[0]]+wbinn[FIT_IWBIN_GT_ET:],[0.03]+bg[FIT_IWBIN_GT_ET:],3)
			#pars=np.polyfit([wbinn[0]]+wbinn[FIT_IWBIN_GT_ET:27],[0]+bg[FIT_IWBIN_GT_ET:27],3)
			#pars=np.polyfit([wbinn[0]]+wbinn[FIT_IWBIN_GT_ET:],[0]+bg[FIT_IWBIN_GT_ET:28]+[bg[27]],3)

			#pars=np.polyfit([wbinn[0]]+wbinn[FIT_IWBIN_GT_ET:],[0]+bg[FIT_IWBIN_GT_ET:],2)

			#pars=np.polyfit([wbinn[0]]+wbinn[6:],[0]+bg[6:],4)
			#pars=np.polyfit([wbinn[0]]+wbinn[6:],[0]+bg[6:],5)
			#pars=np.polyfit([wbinn[0]]+wbinn[11:],[0]+bg[11:],6)
			#pars=np.polyfit([wbinn[0]]+wbinn[FIT_IWBIN_GT_ET:],[0]+bg[FIT_IWBIN_GT_ET:],3)
			
			#! Fit condition 2: Assume that there is an increase at low W, but there are 
			#!  bins in this region where sim3pi does not match exp, leading to increased estimation of BG
			#! + Remove these bins 
			#pars=np.polyfit(wbinn[0:2]+wbinn[6:],bg[0:2]+bg[6:],8)

			#! obtain bg from fit
			bg_fitl=[pars[3] + pars[2]*x**1 + pars[1]*x**2 + pars[0]*x**3 for x in wbinn]
			#bg_fitl=[pars[2] + pars[1]*x**1 + pars[0]*x**2 for x in wbinn]
			#ax.plot(wbinn,bg_fitl,'-',c='m')
			
			#! Refit if any value from bg_fitl<0, which should only happen in the extrapolated region 
			#! and therein most probably in the first 3 bins
			#! Note some refitting settings to "keep things simple":
			#! 1. The precision of comparison to 0 is at 0.01 level i.e. negative values of greater order are treated as 0 because:
			#!   + Precision for background at 0.01% level is enough
			#!   + Additionally, fitting can get complicated in trying to get values to >0 with such levels of precision
			#! 2. Note that refits are limited upto a maximum of 10 times to avoid infinite loops
			refit=True
			nitr_refit=0
			print "refitting for iplt==",PLT_NAME[iplt]
			while (refit==True and nitr_refit<10):
				vals_lt_0=[]
				for i,x in enumerate(wbinn):
					if bg_fitl[i]<0 and not np.isclose(math.fabs(bg_fitl[i]),0,atol=1e-2): 
						vals_lt_0.append(bg_fitl[i])
				if len(vals_lt_0)>0:
					#print len(vals_lt_0)
					print vals_lt_0
					val=math.fabs(min(vals_lt_0))
					cnstrnt_rft=val
					print "Adjusting contraint for fit to go to %f instead of 0 in the first W bin"%cnstrnt_rft
					pars=np.polyfit([wbinn[0]]+wbinn[FIT_IWBIN_GT_ET:],[cnstrnt_rft]+bg[FIT_IWBIN_GT_ET:],3)
					bg_fitl=[pars[3] + pars[2]*x**1 + pars[1]*x**2 + pars[0]*x**3 for x in wbinn]
					nitr_refit+=1
				else:
					refit=False
			print "Number of refit iterations for %s=%d"%(PLT_NAME[iplt],nitr_refit)
			print "bg_fitl for %s:"%(PLT_NAME[iplt])
			print bg_fitl
			#! plot
			ax.plot(wbinn,bg_fitl,'-',c='m')

			#! draw legend
			ax.legend(loc='upper right',prop={'size':8})
			#! save figure
			fig.savefig('%s/bg_test_fit.png'%(outdir_q2w_plt))
			fig.savefig('%s/bg_test_fit.pdf'%(outdir_q2w_plt))
			fig.savefig('%s/bg_test_fit.eps'%(outdir_q2w_plt))

def get_GP():
	'''
	+ Get parameters of gaussian fits to MM and MM2 distributions for ER and SR
	+ These parameters are used when integral-normalizing SR to ER
	+ Uses GP.txt output by $SUBSTUDIES/study_MM_diff_ER_SR/study_MM_diff_ER_SR.py
		+ GP.txt format: wbinnum,wbinlabel,mu_ER,sg_ER,mu_SR,sg_SR
	'''
	#! Create structure to fill and return: GP[w][plt][dtyp]=(mu,sg)
	#! + Note there is no q2 index as pars are obtained for integrated q2
	#! + First setup indexing vars
	NDTYP=2
	ER,SR=range(NDTYP)
	DTYP_NAME=["ER","SR"]

	NPLTS=2
	MM,MM2=range(NPLTS)
	PLT_NAME=['MM','MM2']

	NWBINS=29

	GP=[[[0 for k in range(NDTYP)] for j in range(NPLTS)] for i in range(len(WBIN_LEL))]

	#! Setup input data file
	fin=[0 for i in range(NPLTS)]
	datadir=os.path.join(os.environ['STUDY_MM_DIFF_ER_SR_DATADIR'],'results_092817/2.00-5.00/EI')
	fin[MM] =open(os.path.join(datadir, 'MM/GP.txt'),'r')
	fin[MM2]=open(os.path.join(datadir, 'MM2/GP.txt'),'r')

	#! Begin to fill structure
	for iplt in range(NPLTS):
		data=fin[iplt].readlines()
		for line in data:
			words=line.rstrip().split(' ')
			#print words
			iw,mu_ER,sg_ER,mu_SR,sg_SR=int(words[0])-1,float(words[2]),float(words[3]),float(words[4]),float(words[5])
			#print iw,mu_ER,sg_ER,mu_SR,sg_SR
			GP[iw][iplt][ER]=[mu_ER,sg_ER]
			GP[iw][iplt][SR]=[mu_SR,sg_SR]
	return GP
# GP=get_GP()
# print GP
# sys.exit()

#! begin main
if TEST_FIT_ONLY:
	test_fit_bg()
	sys.exit()

if NORM_OPT_SR_TO_ER=='integral':
	GP=get_GP()

#! Fill hmm
for iq,q2bin_le in enumerate(Q2BIN_LEL):
	if DBG==True and iq>0: continue #! debug
	q2min=q2bin_le
	q2max=Q2BIN_UEL[iq]
	for iw,wbin_le in enumerate(WBIN_LEL):
		#! DBGW (3)
		if DBG==True and (iw!=0 and iw!=15 and iw!=28): continue #! debug #! (iw!=0 and iw!=15 and iw!=28)
		#if DBG==True and iw+1 > 10: continue
		#if DBG==True and (iw+1!=2 and iw+1!=3 and iw+1!=4 and iw+1!=5): continue
		wmin=wbin_le
		wmax=wbin_le+WBINW
		cut_q2w=ROOT.TCut("Q2>%f && Q2<%f && W>%f && W<%f"%(q2min,q2max,wmin,wmax))
		print "q2w=",cut_q2w.GetTitle()
		for r in D:
			idtyp,iplt=r[0],r[1]
			
			cut_top=ROOT.TCut("top==2")
			cut=ROOT.TCut(cut_q2w)
			cut+=cut_top
			print "TTree::Draw() for",DTYP_NAME[idtyp],PLT_NAME[iplt]
			#print T[idtyp].GetName()
			T[idtyp].Draw(DRAW_CMD[iplt],cut,"",NENTRIES)#[idtyp][isctr][icorr],cut)

			#! Store histogram
			#ROOT.gStyle.SetOptStat("ne")
			htmp=ROOT.gDirectory.Get("hmmcmd")#(draw_cmd_hst_name[idtyp][isctr][icorr])
			hmm[idtyp][iq][iw][iplt]=htmp.Clone()#"hmm_%s_%s"%(DTYP_NAME[idtyp],CORR_NAME[icorr]))
			hmm[idtyp][iq][iw][iplt].SetName("h%s_%s"%(PLT_NAME[iplt],DTYP_NAME[idtyp]))
			title="%s for Q^{2}=[%.2f,%.2f) GeV^{2}, W=[%.3f,%.3f) GeV"%(PLT_TITLE[iplt],q2min,q2max,wmin,wmax)
			hmm[idtyp][iq][iw][iplt].SetTitle(title)#"%s %.2f-%.2f_%.3f-%.3f"%(PLT_TITLE[iplt],q2min,q2max,wmin,wmax))#("%s"%cut.GetTitle())
			#print hmm[idtyp][i][j][icorr].GetName()

#! + Plots
ROOT.gStyle.SetOptStat(0) #"n"
ROOT.gStyle.SetOptFit(1)#111)
#! Aesthetic related objects
CLRS_DTYP=[ROOT.gROOT.ProcessLine("kBlue"),ROOT.gROOT.ProcessLine("kRed"),ROOT.gROOT.ProcessLine("kGreen")]

#! First plot MM/MM2 without BG analysis
#! + i.e. simply ER and SR plots together as was done before BG analysis
plot_without_BG_analysis()

#! Now plot with BG analysis

#! [07-13-17] Create and fill bg structure in this loop
#! + Create structure to hold estimating BG(%):BG[q2][w][plt]=(bg,qbin,wbin)
#! + Initialize to 'nan' because then where BG is not calculated will not be plotted
BG=[[[float('nan') for k in range(NPLTS)] for j in range(len(WBIN_LEL))] for i in range(len(Q2BIN_LEL))]

#! 2 iterations:
#! 1. For W bins where hSR3PI *can* be normalized to hdiff
#!   + Using data from this iteration calculate bg and obtain fit
#! 2. For W bins where hSR3PI *cannot* be normalized to hdiff
#!   + Normalize hSR3PI using norm_hSR3PI_based_on_bg_fit(
#!   + No bg to be calculate here
for itr in range(2):
	for iq,q2bin_le in enumerate(Q2BIN_LEL):
		if DBG==True and iq>0: continue #! debug
		q2min=q2bin_le
		q2max=Q2BIN_UEL[iq]
		outdir_q2w="%s/with_BG_ana/%.2f-%.2f"%(OUTDIR,q2min,q2max)
		if not os.path.exists(outdir_q2w):
			os.makedirs(outdir_q2w)
		for iw,wbin_le in enumerate(WBIN_LEL):
			#! DBGW (3)
			if DBG==True and (iw!=0 and iw!=15 and iw!=28): continue #! debug #!(iw!=0 and iw!=15 and iw!=28)
			#if DBG==True and iw+1 > 10: continue
			#if DBG==True and (iw+1!=2 and iw+1!=3 and iw+1!=4 and iw+1!=5): continue
			wmin=wbin_le
			wmax=wbin_le+WBINW
		
			wbin_norm3pi_not_possible=False
			if True in [np.isclose(wmin,wbin[0]) and np.isclose(wmax,wbin[1]) for wbin in WBINS_NORM3PI_NOT_POSSIBLE]:
				wbin_norm3pi_not_possible=True
			
			#! If using 2-iterations mode where in the 1st iteration bg-fit is made where wbin_norm3pi_IS_possible and
			#! in the 2nd iteration normalization in wbin_norm3pi_not_possible is based on nf extrapolated using the bg-fit
			#! + Note that some code may need to be commented/un-commented manually in this mode
			# if   (itr+1==1) and wbin_norm3pi_not_possible==True: continue
			# elif (itr+1==2) and wbin_norm3pi_not_possible==False: continue

			#! If using 1-iteration mode where wbin_norm3pi_not_possible is simply not normalized and 
			#! bg-fit from wbin_norm3pi_IS_possible is used to extrapolate bg where wbin_norm3pi_not_possible
			if (itr+1==2): continue

			print "Going to plot hmm and hmm2 for %.2f-%.2f_%.3f-%.3f"%(q2min,q2max,wmin,wmax)

			#! Create cmm[plt]
			cmm=[0 for i in range(NPLTS)]
			#! Now plot
			for iplt in range(NPLTS):
				outdir_q2w_plt="%s/%s"%(outdir_q2w,PLT_NAME[iplt])
				if not os.path.exists(outdir_q2w_plt):
					os.makedirs(outdir_q2w_plt)
				cname="c%s_%.3f-%.3f"%(PLT_NAME[iplt],wmin,wmax)
				cmm[iplt]=ROOT.TCanvas(cname,cname)
				#l=ROOT.TLegend(0.1,0.3,0.3,0.4)#,"","NDC");
				#! hists aesthetics	
				for idtyp in range(NDTYP):
					hmm[idtyp][iq][iw][iplt].SetLineColor(CLRS_DTYP[idtyp])
					hmm[idtyp][iq][iw][iplt].SetMarkerColor(CLRS_DTYP[idtyp])
				#! Draw hists
				#! First directly get copy of hists so as to not use tedious indices
				hER=hmm[ER][iq][iw][iplt]
				hSR=hmm[SR][iq][iw][iplt]
				hSR3PI=hmm[SR3PI][iq][iw][iplt]
				hER.Sumw2()
				hSR.Sumw2()
				hSR3PI.Sumw2()
				#! setup hER
				hER.SetXTitle("%s [%s]"%(PLT_TITLE[iplt], PLT_UNIT[iplt]))
				cmm[iplt].SetLeftMargin(0.20)
				hER.GetYaxis().SetTitleOffset(1.5)
				hER.SetYTitle("N_{entries}")

				#! setup SR 
				if NORM_OPT_SR_TO_ER=='maximum':
					#! max-norm
					max_ER=hER.GetMaximum()
					if max_ER==0:max_ER=1
					max_SR=hSR.GetMaximum()
					if max_SR==0:max_SR=1
					scl_fctr_SR=max_ER/max_SR
					if scl_fctr_SR==0:scl_fctr_SR=1
					hSR.Scale(scl_fctr_SR)
				elif NORM_OPT_SR_TO_ER=='integral':
					#! itgrl-norm
					#! + setup integral limits
					mu_ER,sg_ER=GP[iw][iplt][ER][0],GP[iw][iplt][ER][1]
					mu_SR,sg_SR=GP[iw][iplt][SR][0],GP[iw][iplt][SR][1]

					min_ER,max_ER=mu_ER-3*sg_ER,mu_ER+3*sg_ER
					min_SR,max_SR=mu_SR-3*sg_SR,mu_SR+3*sg_SR

					bin_min_ER,bin_max_ER=hER.FindBin(min_ER),hER.FindBin(max_ER)
					bin_min_SR,bin_max_SR=hSR.FindBin(min_SR),hER.FindBin(max_SR)
					#! + itgrl for hER
					itgrl_ER=hER.Integral(bin_min_ER,bin_max_ER)
					#! + itgrl for hSR
					itgrl_SR=hSR.Integral(bin_min_SR,bin_max_SR)
					#! + scale hSR as per itgrl ratio
					scl_fctr_SR=itgrl_ER/itgrl_SR
					hSR.Scale(scl_fctr_SR)
	
				#! Now normalize SR3PI to hER-hSR (=hdiff) and draw
				#! 1. get hdiff=hER-hSR
				hdiff=hER.Clone("hdiff")
				hdiff.Sumw2()
				hdiff.Add(hSR,-1)
				hdiff.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
				hdiff.SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
				#! 2. Now do the normalization or obtain it from bg-fit (done in a previous iteration)
				if wbin_norm3pi_not_possible==True:
					#! + Neither normalize nor obtain nf from bg-fit which is done the previous iteration
					continue
				else:
					#! range norm (note that the obtained hist has Sumw2() set in the method)
					print "range norm using NL"
					hSR3PIn,nf=atlib.norm_hist_range(hdiff,hSR3PI,NL[iplt])
					#NF[iq][iw][iplt]=nf
				
				if BG_DISPLAY_MODE=='full':
					#! Adjust lims of y-axis to go below 0 to show fluctuations in hdiff 
					#! and draw line at y=0
					min_hdiff=hdiff.GetMinimum()
					min_hdiff=min_hdiff-(10/100)*math.fabs(min_hdiff)
					hER.SetMinimum(min_hdiff)

					#! Now draw hER,hSR, hdiff and hSR3PIn
					hER.Draw()
					hSR.Draw("same")
					hdiff.Draw("same")
					ZERO_YAXIS[iplt].Draw("same")
					hSR3PIn.Draw("same")
				elif BG_DISPLAY_MODE=='simple':
					#! Draw hER,hSR, hsum(=hSR3PIn+hSR)
					#! 1. Get hsum
					hsum=hSR3PIn.Clone("hsum")
					hsum.Sumw2()
					hsum.Add(hSR)
					hsum.SetLineColor(ROOT.gROOT.ProcessLine("kBlack"))
					hsum.SetMarkerColor(ROOT.gROOT.ProcessLine("kBlack"))
					#! 2. Draw hER,hSR, and hsum
					hER.Draw()
					hSR.Draw("same")
					hSR3PIn.Draw("same")
					hsum.Draw("same")

				#! Create and add entries: hists and cut-lines
				#! legend
				l=ROOT.TLegend(0.66,0.72,0.9,0.9)#,"","NDC");	
				if BG_DISPLAY_MODE=='full':
					for hist,label in zip([hER,hSR,hdiff,hSR3PIn],["exp","sim","exp - sim","sim-3pi-PS"]):
						l.AddEntry(hist,label,"lp")
				if BG_DISPLAY_MODE=='simple':
					for hist,label in zip([hER,hSR,hSR3PIn,hsum],["exp","sim","sim-3pi-PS","sim + sim-3pi-PS"]):
						l.AddEntry(hist,label,"lp")

				#! Draw cut lines and add them to legend
				if iplt==MM:
					#! Draw cut lines
					MM_T2_CUTL_EI.SetY1(hER.GetMinimum())
					MM_T2_CUTL_EI.SetY2(hER.GetMaximum())
					MM_T2_CUTH_EI.SetY1(hER.GetMinimum())
					MM_T2_CUTH_EI.SetY2(hER.GetMaximum())
					MM_T2_CUTL_EI.Draw("same")
					MM_T2_CUTH_EI.Draw("same")
					#! add to legend
					l.AddEntry(MM_T2_CUTL_EI,"%.2f GeV < MM < %.2f GeV"%(EI_MML,EI_MMH),"l")
					#! Draw NL line
					if BG_DISPLAY_MODE=='full':
						if wbin_norm3pi_not_possible==True:
							lmin,lmax=NL1_LINE_MIN[MM],NL1_LINE_MAX[MM]
						else:
							lmin,lmax=NL_LINE_MIN[MM],NL_LINE_MAX[MM]
						for line in [lmin,lmax]:
							line.SetY1(hER.GetMinimum())
							line.SetY2(hER.GetMaximum())
							line.Draw("same")
						#! add to legend
						l.AddEntry(lmin,"%.2f GeV < norm-range < %.2f GeV"%(lmin.GetX1(),lmax.GetX1()),"l")
					elif BG_DISPLAY_MODE=='simple':
						if wbin_norm3pi_not_possible==True:
							lmax=NL1_LINE_MAX[MM]
						else:
							lmax=NL_LINE_MAX[MM]
						for line in [lmax]:
							line.SetY1(hER.GetMinimum())
							line.SetY2(hER.GetMaximum())
							line.Draw("same")
						#! add to legend
						#! + Note that this normalization range is not technically true and it should be same as for 'full'
						#! + However it is arguably consistent with the idea of making a 'simple' display
						l.AddEntry(lmax,"norm-range < %.2f GeV"%(lmax.GetX1()),"l")
				elif iplt==MM2:
					#! EI
					MM2_T2_CUTL_EI.SetY1(hER.GetMinimum())
					MM2_T2_CUTL_EI.SetY2(hER.GetMaximum())
					MM2_T2_CUTH_EI.SetY1(hER.GetMinimum())
					MM2_T2_CUTH_EI.SetY2(hER.GetMaximum())
					MM2_T2_CUTL_EI.Draw("same")
					MM2_T2_CUTH_EI.Draw("same")
					#! add to legend
					l.AddEntry(MM2_T2_CUTL_EI,"%.2f GeV^{2} < MM^{2} < %.2f GeV^{2}"%(EI_MM2L,EI_MM2H),"l")
					#! Draw NL line
					if BG_DISPLAY_MODE=='full':
						if wbin_norm3pi_not_possible==True:
							lmin,lmax=NL1_LINE_MIN[MM2],NL1_LINE_MAX[MM2]
						else:
							lmin,lmax=NL_LINE_MIN[MM2],NL_LINE_MAX[MM2]
						for line in [lmin,lmax]:
							line.SetY1(hER.GetMinimum())
							line.SetY2(hER.GetMaximum())
							line.Draw("same")
						#! add to legend
						l.AddEntry(lmin,"%.2f GeV < norm-range < %.2f GeV"%(lmin.GetX1(),lmax.GetX1()),"l")
					elif BG_DISPLAY_MODE=='simple':
						if wbin_norm3pi_not_possible==True:
							lmax=NL1_LINE_MAX[MM2]
						else:
							lmax=NL_LINE_MAX[MM2]
						for line in [lmax]:
							line.SetY1(hER.GetMinimum())
							line.SetY2(hER.GetMaximum())
							line.Draw("same")
						#! add to legend
						#! + Note that this normalization range is not technically true and it should be same as for 'full'
						#! + However it is arguably consistent with the idea of making a 'simple' display
						l.AddEntry(lmax,"norm-range < %.2f GeV"%(lmax.GetX1()),"l")
				#! Draw legend
				l.Draw()
				cmm[iplt].SaveAs("%s/%s.png"%(outdir_q2w_plt,cname))
				cmm[iplt].SaveAs("%s/%s.pdf"%(outdir_q2w_plt,cname))
				if itr+1==1:
					#! [07-13-17] fill BG(%)
					if wbin_norm3pi_not_possible==True:
						BG[iq][iw][iplt]=float('nan')#float('nan')
					else:
						if   iplt==MM:  cut_min,cut_max=EI_MML,EI_MMH
						elif iplt==MM2: cut_min,cut_max=EI_MM2L,EI_MM2H
						#! Get integration limits in terms of bin number
						cut_bin_min=hER.FindBin(cut_min)
						cut_bin_max=hER.FindBin(cut_max)
						#! Now calculate BG
						num=hSR3PIn.Integral(cut_bin_min,cut_bin_max)
						den=hER.Integral(cut_bin_min,cut_bin_max)
						bg=(num/den)*100
						#! Fill BG structure
						BG[iq][iw][iplt]=bg
	if itr+1==1:
		#! Plot and write BG
		plot_fit_and_write_bg(BG)

# #! If wanting to keep TCanvas open till program exits				
# if not ROOT.gROOT.IsBatch():
# 	plt.show()
# 	# wait for you to close the ROOT canvas before exiting
# 	wait(True)

#if __name__ == "__main__":	
#	if len(sys.argv)==2:
#		plot_fid(sys.argv[1])
#	elif len(sys.argv)==3:
#		plot_fid(sys.argv[1],int(sys.argv[2])) 
	
