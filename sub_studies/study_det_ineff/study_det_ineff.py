#!/usr/bin/python
from __future__ import division

import os,sys
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict
import itertools

import numpy as np

import math

from cuts_theta_vs_p_e16 import NCUTLNS,H,L,CUT
print "NCUTLNS,H,L=",NCUTLNS,H,L

'''
+ {e1f,e16}*{"Q2,W"}*(sector}*{prtcl}*{"plots"}*{AT-cuts,EI-cuts}, where
	+ "Q2,W":
		+ full Q2-W range
		+ Q2 bins
		+ Q2,W bins (if there is enough statistics)
	+ "plots"
		+ theta
		+ p
		+ theta vs. p
		+ theta vs. phi
		+ SC paddles
		+ paddle vs. theta vs. p (how to visualize this?)

+ Usage: study_det_ineff.py expt=<e1f/e16> top=[2] make_kin_plots[=True] stdy_bad_sc_pd_theta_vs_p_crltn[=False] debug[=False]
'''
USAGE="study_det_ineff.py expt=<e1f/e16> top=[2] make_kin_plots[=True] stdy_bad_sc_pd_theta_vs_p_crltn[=False] debug[=False]"
#! *** Get arguments from user *** #!
if len(sys.argv)<2:
		sys.exit('usage: %s'%USAGE)
expt=sys.argv[1]
if expt!='e1f' and expt!='e16':
	sys.exit('expt=%s is not valid. usage: USAGE'%(expt,USAGE))
print "expt=",expt
#! For now ignore e1f
if expt=="e1f":
	sys.exit("Not implemented for e1f yet")

TOP=2
if len(sys.argv)>2: #! i.e. top entered by user
	TOP=int(sys.argv[2])
if TOP<1 or TOP>4:
	sys.exit("Valid tops=1,2,3 or 4")
print "TOP=",TOP

MAKE_KIN_PLOTS=True
if len(sys.argv)>3: #! i.e. make_kin_plots entered by user
        if sys.argv[3]=="True":
                MAKE_KIN_PLOTS=True
        elif sys.argv[3]=="False":
                MAKE_KIN_PLOTS=False
        else:
                sys.exit("Please enter make_kin_plots as True/False only.")
print "MAKE_KIN_PLOTS=",MAKE_KIN_PLOTS

STDY_BAD_SC_PD_THETA_VS_P_CRLTN=False
if len(sys.argv)>4: #! i.e. stdy_bad_sc_pd_theta_vs_p_crltn entered by user
	if sys.argv[4]=="True":
        	STDY_BAD_SC_PD_THETA_VS_P_CRLTN=True
	elif sys.argv[4]=="False":
		STDY_BAD_SC_PD_THETA_VS_P_CRLTN=False
	else:
		sys.exit("Please enter stdy_bad_sc_pd_theta_vs_p_crltn as True/False only.")
print "STDY_BAD_SC_PD_THETA_VS_P_CRLTN=",STDY_BAD_SC_PD_THETA_VS_P_CRLTN

DEBUG=False
if len(sys.argv)>5: #! i.e. debug entered by user
        if sys.argv[5]=="True":
                DEBUG=True
        elif sys.argv[5]=="False":
                DEBUG=False
        else:
                sys.exit("Please enter debug as True/False only.")
print "DEBUG=",DEBUG
#sys.exit()
#! ****

#! *** Prepare input data structure: FIN/T[dtyp] ***
NDTYP=3
ER,SR,ST=range(NDTYP)
DTYP_NAME=["ER","SR","ST"]

FIN=[[] for i in range(NDTYP)]
T=[[] for i in range(NDTYP)]
#print "FIN=",FIN
#print "T=",T
#! ***

#! *** Prepare output datadir***
if expt=="e1f":
	DATADIR_OUTPUT=os.environ['STUDY_DET_INEFF_E1F']
elif expt=="e16":
	DATADIR_OUTPUT=os.environ['STUDY_DET_INEFF_E16']
print "DATADIR_OUTPUT=",DATADIR_OUTPUT
#! ***

#! *** Now fill T[dtyp][rctn] ***
DATADIR_INPUT=[[] for i in range(NDTYP)]
sfx=''
if expt=='e16':sfx='_E16'
DATADIR_INPUT[ER]=os.path.join(os.environ['D2PIDIR_EXP%s'%sfx],"data_det_ineff_012416")
DATADIR_INPUT[SR]=os.path.join(os.environ['D2PIDIR_SIM%s'%sfx],"data_det_ineff_012416")
DATADIR_INPUT[ST]=os.path.join(os.environ['D2PIDIR_SIM%s'%sfx],"data_det_ineff_012416")

print "DATADIR_INPUT[ER]=",DATADIR_INPUT[ER]
print "DATADIR_INPUT[SR]=",DATADIR_INPUT[SR]
print "DATADIR_INPUT[ST]=",DATADIR_INPUT[ST]

for idtyp in range(NDTYP): 
	print "Getting T for",DTYP_NAME[idtyp]#,RCTN_NAME[irctn]
	FIN[idtyp]=ROOT.TFile("%s/d2pi%s_ddetineff.root"%(DATADIR_INPUT[idtyp],DTYP_NAME[idtyp]))

	if DTYP_NAME[idtyp]=="ER" or DTYP_NAME[idtyp]=="SR":
		T[idtyp]=FIN[idtyp].Get("d2piR/tR")
	elif DTYP_NAME[idtyp]=="ST":
		T[idtyp]=FIN[idtyp].Get("d2piT/tT")
#print "T=",T
#print FIN[ER].GetName()
#print FIN[SR].GetName()
#print FIN[ST].GetName()
#print T[ER].GetName()
#print T[ER].GetEntries()
#print T[SR].GetName()
#print T[SR].GetEntries()
#print T[ST].GetName()
#print T[ST].GetEntries()
#sys.exit()
#! ***

#! *** Now make plots ***

#! *** Setup Q2,W information ***
WMIN=1.400
WMAX=2.125
WBINW=0.025
#! Note that last the elements of the Q2BIN_LEL/UEL represent a "hackish" way to 
#! to plot variables over full Q2 range
Q2BIN_LEL=[2.00,2.40,3.00,3.50,4.20,2.00] 
Q2BIN_UEL=[2.40,3.00,3.50,4.20,5.00,5.00]
WBIN_LEL=np.arange(WMIN,WMAX,WBINW)
WBIN_UEL=np.arange(WMIN+WBINW,WMAX+WBINW,WBINW)
#print Q2BIN_LEL
#print WBIN_LEL

Q2Wdomain="Q2" #! should be made a user input
#! + Set up appropriate Q2,W domain list
if Q2Wdomain=="Q2-W":
        dle=[Q2BIN_LEL,WBIN_LEL]
        due=[Q2BIN_UEL,WBIN_UEL]
elif Q2Wdomain=="Q2":
        dle=[Q2BIN_LEL]
        due=[Q2BIN_UEL]
DLE=list(itertools.product(*dle))
DUE=list(itertools.product(*due))
#! ***

#! Setup sector information 
NSCTR=6
#! ***

#! Setup particle information
NPRTCL=4
E,P,PIP,PIM=range(NPRTCL)
PRTCL_NAME=["e","p","pip","pim"]
#! ***

#! Setup "plot" information"
NPLT=5
THETA,P,THETA_VS_P,THETA_VS_PHI,SC_PD=range(NPLT)
PLT_NAME=["theta","p","theta_vs_p","theta_vs_phi","sc_pd"]
if DEBUG and (not STDY_BAD_SC_PD_THETA_VS_P_CRLTN):
	NPLT=1
	#THETA=range(NPLT)
	#PLT_NAME=["theta"]
	SC_PD=range(NPLT)
	PLT_NAME=["sc_pd"]
elif DEBUG and STDY_BAD_SC_PD_THETA_VS_P_CRLTN:
	NPLT=1
	THETA_VS_P=range(NPLT)
	PLT_NAME=["theta_vs_p"]

#! setup binning information for plots
BNG={
("theta","e"):[100,0,60],
("theta","p"):[100,0,60],
("theta","pip"):[100,0,120],
("theta","pim"):[100,0,120],
("p","e"):[100,0,5],
("p","p"):[100,0,5],
("p","pip"):[100,0,5],
("p","pim"):[100,0,5],
("theta_vs_p","e"):[100,0,5,100,0,60],
("theta_vs_p","p"):[100,0,5,100,0,60],
("theta_vs_p","pip"):[100,0,5,100,0,120],
("theta_vs_p","pim"):[100,0,5,100,0,120],
("theta_vs_phi","e"):[100,0,360,100,0,60],
("theta_vs_phi","p"):[100,0,360,100,0,60],
("theta_vs_phi","pip"):[100,0,360,100,0,120],
("theta_vs_phi","pim"):[100,0,360,100,0,120],
("sc_pd","e"):[[70,100,170],[70,200,270],[70,300,370],[70,400,470],[70,500,570],[70,600,670]],
("sc_pd","p"):[[70,100,170],[70,200,270],[70,300,370],[70,400,470],[70,500,570],[70,600,670]],
("sc_pd","pip"):[[70,100,170],[70,200,270],[70,300,370],[70,400,470],[70,500,570],[70,600,670]],
("sc_pd","pim"):[[70,100,170],[70,200,270],[70,300,370],[70,400,470],[70,500,570],[70,600,670]]
}
#! *** 

#! Create and fill draw_cmd[sctr][prtcl][nplt]
draw_cmd=[[[[]for k in range(NPLT)]for j in range(NPRTCL)]for i in range(NSCTR)]
dl=[range(NSCTR),range(NPRTCL),range(NPLT)]
d=list(itertools.product(*dl))
for r in d:
	isctr,iprtcl,iplt=r[0],r[1],r[2]
	plt=PLT_NAME[iplt]
        prtcl=PRTCL_NAME[iprtcl]
        if plt=="theta" or plt=="p" or plt=="sc_pd":
        	if plt=="sc_pd": nbins,xmin,xmax=BNG[plt,prtcl][isctr][0],BNG[plt,prtcl][isctr][1],BNG[plt,prtcl][isctr][2]
                else:           nbins,xmin,xmax=BNG[plt,prtcl][0],BNG[plt,prtcl][1],BNG[plt,prtcl][2]
                draw_cmd[isctr][iprtcl][iplt]="%s_%s>>hcmd(%d,%d,%d)"%(plt,prtcl,nbins,xmin,xmax)
        elif plt=="theta_vs_p":
        	nbinsx,xmin,xmax=BNG[plt,prtcl][0],BNG[plt,prtcl][1],BNG[plt,prtcl][2]
                nbinsy,ymin,ymax=BNG[plt,prtcl][3],BNG[plt,prtcl][4],BNG[plt,prtcl][5]
                draw_cmd[isctr][iprtcl][iplt]="theta_%s:p_%s>>hcmd(%f,%f,%f,%f,%f,%f)"%(prtcl,prtcl,nbinsx,xmin,xmax,nbinsy,ymin,ymax)
        elif plt=="theta_vs_phi":
        	nbinsx,xmin,xmax=BNG[plt,prtcl][0],BNG[plt,prtcl][1],BNG[plt,prtcl][2]
                nbinsy,ymin,ymax=BNG[plt,prtcl][3],BNG[plt,prtcl][4],BNG[plt,prtcl][5]
                draw_cmd[isctr][iprtcl][iplt]="theta_%s:phi_%s>>hcmd(%f,%f,%f,%f,%f,%f)"%(prtcl,prtcl,nbinsx,xmin,xmax,nbinsy,ymin,ymax)
#print draw_cmd
#sys.exit()
#! ***

#! outdir
outdir=os.path.join(DATADIR_OUTPUT,"results_top%d"%TOP)
if not os.path.exists(outdir):
	os.makedirs(outdir)
#! ***

#! Define any global TCuts
#! top==2
#TOP=2
cut_top=ROOT.TCut("top==%d"%TOP)
#! ***

#! First create structure to store hkin
if Q2Wdomain=="Q2-W": #! Create h[dtyp][q2-w][sctr][prtcl][nplt]
	hkin=[[[[[[]for m in range(NPLT)]for l in range(NPRTCL)]for k in range(NSCTR)]for j in range(len(Q2BIN_LEL)*len(WBIN_LEL))]for i in range(NDTYP)]
elif Q2Wdomain=="Q2": #! Create h[dtyp][q2][sctr][prtcl][nplt]
	hkin=[[[[[[]for m in range(NPLT)]for l in range(NPRTCL)]for k in range(NSCTR)]for j in range(len(Q2BIN_LEL))]for i in range(NDTYP)]

if MAKE_KIN_PLOTS:
	#! Now fill hkin
	for iq2wb,q2wbin_le in enumerate(DLE):
		if DEBUG and iq2wb>0: continue #! debug
		q2min,q2max=q2wbin_le[0],DUE[iq2wb][0]
		if   Q2Wdomain=="Q2-W": wmin,wmax=q2wbin_le[1],DUE[iq2wb][1]
		elif Q2Wdomain=="Q2":   wmin,wmax=WMIN,WMAX
		dl=[range(NSCTR),range(NDTYP),range(NPRTCL),range(NPLT)]
		d=list(itertools.product(*dl))
		for r in d:
			isctr,idtyp,iprtcl,iplt=r[0],r[1],r[2],r[3]
			sctr=isctr+1		
			dtyp=DTYP_NAME[idtyp]
			prtcl=PRTCL_NAME[iprtcl]
			plt=PLT_NAME[iplt]
			#! Make TCut
			cut_q2w=ROOT.TCut("Q2>%f && Q2<%f && W>%f && W<%f"%(q2min,q2max,wmin,wmax))
			cut=ROOT.TCut(cut_q2w)
			if dtyp=="ER" or dtyp=="SR": #! not for ST
				cut+=cut_top #! cut_top defined earlier
				cut_sctr=ROOT.TCut("sector_%s==%d"%(prtcl,sctr))
				cut+=cut_sctr
			print "TTree::Draw() for q2min,q2max,wmin,wmax,sctr,dtyp,prtcl,plt=",q2min,q2max,wmin,wmax,sctr,dtyp,prtcl,plt
			if T[idtyp].GetEntries()>0:
				T[idtyp].Draw(draw_cmd[isctr][iprtcl][iplt],cut)
				#! Store histogram
				#ROOT.gStyle.SetOptStat("ne")
				htmp=ROOT.gDirectory.Get("hcmd")
				hkin[idtyp][iq2wb][isctr][iprtcl][iplt]=htmp.Clone()
				hkin[idtyp][iq2wb][isctr][iprtcl][iplt].SetName("h_%s_%s_%s_s%d"%(dtyp,prtcl,plt,sctr))
				hkin[idtyp][iq2wb][isctr][iprtcl][iplt].SetTitle("%s_%s_%s %.2f-%.2f_%.3f-%.3f"%(prtcl,plt,dtyp,q2min,q2max,wmin,wmax))
			else:
				print "TTree name=",T[idtyp].GetName()
				print "TTree entries=",T[idtyp].GetEntries()
				sys.exit("TTree has no entries. Because of what appears to be a bug in ROOT, this can lead, ina subtle manner, to deeper problems: TTree::Draw() should return a Null pointer in this case, however, it seems that ROOT returns the pointer to the previous histogram draw by TTree::Draw and therefore I end up referencing the wrong histogram. A workaround could be to use unique names for the histograms in the 'draw_cmd'. However, for now, this will have to do. ")
	
	#! Plot on TCanvas and save h:
	#! + Directly as .png
	#! + in .root file
	fname="kinematics"
	if DEBUG:
                fname+="_dbg"
	fout_root=ROOT.TFile("%s/%s.root"%(outdir,fname),"RECREATE")

	#! CWDTH,CHGHT defined as per 3,2 TCanvas
	CWDTH=2000#500
        CHGHT=1500#300
	for iq2wb,q2wbin_le in enumerate(DLE):
		if DEBUG and iq2wb>0: continue #! debug
		q2min,q2max=q2wbin_le[0],DUE[iq2wb][0]
		if Q2Wdomain=="Q2-W": wmin,wmax=q2wbin_le[1],DUE[iq2wb][1]
		elif Q2Wdomain=="Q2": wmin,wmax=WMIN,WMAX
		q2wbin="%.2f-%.2f"%(q2min,q2max)
		#outdir_q2w="%s/%s"%(outdir,q2wbin) #! for .png
		#if not os.path.exists(outdir_q2w):
                #        os.makedirs(outdir_q2w)
		q2wbindir_root=fout_root.mkdir(q2wbin)#! for .root
		q2wbindir_root.cd()
		dl=[range(NDTYP),range(NPRTCL),range(NPLT)]
		d=list(itertools.product(*dl))
		for r in d:
			idtyp,iprtcl,iplt=r[0],r[1],r[2]
			prtcl=PRTCL_NAME[iprtcl]
			plt=PLT_NAME[iplt]
			print "Going to plot for [%.2f-%.2f_%.3f-%.3f],%s,%s"%(q2min,q2max,wmin,wmax,prtcl,plt)
			if plt=="theta" or plt=="p" or plt=="sc_pd":#! then draw on same canvas
				draw_opt=""
				if plt=="sc_pd": draw_opt="hist"
				cname="c_%s_%s_%s"%(prtcl,plt,q2wbin)
				c=ROOT.TCanvas(cname,cname,CWDTH,CHGHT)
				c.Divide(3,2)
				for isctr in range(NSCTR):
					#l=ROOT.TLegend(0.1,0.3,0.3,0.4)#,"","NDC");
					#! Draw hists
					c.cd(isctr+1)
					hkin[ER][iq2wb][isctr][iprtcl][iplt].SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
					hkin[SR][iq2wb][isctr][iprtcl][iplt].SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
					hkin[ST][iq2wb][isctr][iprtcl][iplt].SetLineColor(ROOT.gROOT.ProcessLine("kGreen"))
					hkin[ER][iq2wb][isctr][iprtcl][iplt].Draw()
					hkin[ER][iq2wb][isctr][iprtcl][iplt].Draw("sames %s"%draw_opt) #! Note only "sames" seems to work here (with "same" statbox with name=hxx is drawn?!)
					#! obtain max_ER(=max height of distribution) used to scale SR,ST 
					max_ER=hkin[ER][iq2wb][isctr][iprtcl][iplt].GetMaximum()
					if max_ER==0:max_ER=1
					#! scale SR and draw
					max_SR=hkin[SR][iq2wb][isctr][iprtcl][iplt].GetMaximum()
					if max_SR==0:max_SR=1
					scl_fctr_SR=max_ER/max_SR
					if scl_fctr_SR==0:scl_fctr_SR=1
					hkin[SR][iq2wb][isctr][iprtcl][iplt].Scale(scl_fctr_SR)
					hkin[SR][iq2wb][isctr][iprtcl][iplt].Draw("sames %s"%draw_opt)
					#! scale ST and draw
					if plt!="sc_pd":
						max_ST=hkin[ST][iq2wb][isctr][iprtcl][iplt].GetMaximum()
						if max_ST==0:max_ST=1
						scl_fctr_ST=max_ER/max_ST
						if scl_fctr_ST==0:scl_fctr_ST=1
						hkin[ST][iq2wb][isctr][iprtcl][iplt].Scale(scl_fctr_ST)
						hkin[ST][iq2wb][isctr][iprtcl][iplt].Draw("sames %s"%draw_opt)
						
					#! Adjust statbox
					#c.Update()
					ROOT.gPad.Update()
					#! statbox for ER
					pt_ER=hkin[ER][iq2wb][isctr][iprtcl][iplt].GetListOfFunctions().FindObject("stats")
					pt_ER.SetTextColor(ROOT.gROOT.ProcessLine("kBlue"))
					#! statbox for SR
					pt_SR=hkin[SR][iq2wb][isctr][iprtcl][iplt].GetListOfFunctions().FindObject("stats")
					pt_SR.SetTextColor(ROOT.gROOT.ProcessLine("kRed"))
					diff=pt_ER.GetY2NDC()-pt_ER.GetY1NDC()
					y2=pt_ER.GetY1NDC()
					y1=y2-diff
					pt_SR.SetY1NDC(y1)
					pt_SR.SetY2NDC(y2)
					#c.Update()
					ROOT.gPad.Update()
					#! statbox for ST
					if plt!="sc_pd":
						pt_ST=hkin[ST][iq2wb][isctr][iprtcl][iplt].GetListOfFunctions().FindObject("stats")
						pt_ST.SetTextColor(ROOT.gROOT.ProcessLine("kGreen"))
						diff=pt_ER.GetY2NDC()-pt_ER.GetY1NDC()
						y2=pt_SR.GetY1NDC()
						y1=y2-diff
						pt_ST.SetY1NDC(y1)
						pt_ST.SetY2NDC(y2)
						#c.Update()
						ROOT.gPad.Update()
					#! Draw legend (not drawing since relevant information in statbox)
					#l.AddEntry(hmm[ER][i][iprtcl][iplt],'ER')
					#l.AddEntry(hmm[ER][i][iprtcl][iplt],'SR')
					#l.Draw()
					#c.SaveAs("%s/%s.png"%(outdir_q2w,cname))#! for .png
					c.Write("%s"%cname,ROOT.gROOT.ProcessLine("TObject::kOverwrite"))#! for .root; latest namecycle
			elif plt=="theta_vs_p" or plt=="theta_vs_phi": #! draw ER and ER on separate canvases
					cname="c_%s_%s_ER_%s"%(prtcl,plt,q2wbin)
					c=ROOT.TCanvas(cname,cname,CWDTH,CHGHT)
					c.Divide(3,2)
					for isctr in range(NSCTR):
						c.cd(isctr+1)
						hkin[ER][iq2wb][isctr][iprtcl][iplt].Draw("colz")
					#c.SaveAs("%s/%s.png"%(outdir_q2w,cname))
					c.Write("%s"%cname,ROOT.gROOT.ProcessLine("TObject::kOverwrite"))#! for .root; latest namecycle

					cname="c_%s_%s_SR_%s"%(prtcl,plt,q2wbin)
					c=ROOT.TCanvas(cname,cname,CWDTH,CHGHT)
					c.Divide(3,2)
					for isctr in range(NSCTR):
						c.cd(isctr+1)
						hkin[SR][iq2wb][isctr][iprtcl][iplt].Draw("colz")
					#c.SaveAs("%s/%s.png"%(outdir_q2w,cname))
					c.Write("%s"%cname,ROOT.gROOT.ProcessLine("TObject::kOverwrite"))#! for .root; latest namecycle

					cname="c_%s_%s_ST_%s"%(prtcl,plt,q2wbin)
					c=ROOT.TCanvas(cname,cname,CWDTH,CHGHT)
					c.Divide(3,2)
					for isctr in range(NSCTR):
						c.cd(isctr+1)
						hkin[ST][iq2wb][isctr][iprtcl][iplt].Draw("colz")
					#c.SaveAs("%s/%s.png"%(outdir_q2w,cname)) #! for .png
					c.Write("%ss"%cname,ROOT.gROOT.ProcessLine("TObject::kOverwrite"))#! for .root; latest namecycle
	fout_root.Close()					
#! ***
#! Now study bad_sc_pd and theta_vs_p dip correlation
if STDY_BAD_SC_PD_THETA_VS_P_CRLTN:

	outdir_sc_pd_theta_vs_p=os.path.join(outdir,"sc_pd_theta_vs_p")
	if not os.path.exists(outdir_sc_pd_theta_vs_p):
		os.makedirs(outdir_sc_pd_theta_vs_p)

	#! Define some functions to get sector and paddle from sc_pd_code
	def get_sector(sc_pd_code):
			return (sc_pd_code-(sc_pd_code%100))/100
	def get_paddle(sc_pd_code):
			return sc_pd_code%100 
	#! Define some other functions used below
	def found_hist_on_tpad(pad):
		ret=False
		next=ROOT.TIter(pad.GetListOfPrimitives())
                while(1):
                	obj=next()
                        if obj==None: break
                        print "Reading",obj.GetName()
                        if obj.InheritsFrom("TH1"):
                        	print "histo:",obj.GetName()
                                ret=True
				break
		return ret
	
	def draw_EI_theta_vs_p_cuts(iprtcl,isctr,pad):
        	if CUT[H][iprtcl][isctr]!=None:
			print "Going to draw EI's theta_vs_p cuts for prtcl,sctr=",PRTCL_NAME[iprtcl],isctr+1
			#! setup draw opt
			if found_hist_on_tpad(pad): draw_opt="same"
                	else:                       draw_opt=""
			#! ***
                	for cuth,cutl in zip(CUT[H][iprtcl][isctr],CUT[L][iprtcl][isctr]):
                        	cuth.Draw(draw_opt)
                                cutl.Draw("same")
		else:
			print "Skip (since not cut exists) Going to draw EI's theta_vs_p cuts for prtcl,sctr=",PRTCL_NAME[iprtcl],isctr+1

	#! Make list of bad_sc_pd_code[NPRTCL][NDTYP]; NPRTCL=(E,P,PIP),NDTYP=(ER,SR)
	#! + Truly, if a paddle is bad, it is bad for for all particles. 
	#! + However, given the kinematics of each particle, every paddle is not accessible to each particle.
	#! + Therefore, a comprehensive list of bad pad paddles is identified by looking at each particle separately 
	#bad_sc_pd_code_list=[
	#[[],[311,505]], #!E
	#[[],[311,324,520]], #!P
	#[[145,340,342,345,346,347,632,633,634,635,636,637,638,639],[145,311,324,336,337,340,342,434,520,532,540,542,636,637,644]] #! PIP
	#]

	#! [02-14-16] Paddles and added from above list after a closer look at sc_pd histograms prompted by the first observations of output from this script (there may be some other changes related to arraning the paddles in numerical order where they were not so ordered)
	#! Paddles removed:
	#! + pip
	#! 1. SR:340
	#! 2. SR:336
	#! Paddles added:
	#! + p
	#! ER:324
	#! + pip
	#! ER:245,337,338,324,434,542,644
	#! SR:338
	#bad_sc_pd_code_list=[
        #[[],[311,505]], #!E
        #[[324],[311,324,520]], #!P
        #[[145,245,324,337,338,340,342,345,346,347,434,542,632,633,634,635,636,637,638,639,644],[145,311,324,337,338,342,434,520,532,540,542,636,637,644]] #! PIP
        #]

	#! [02-16-16] The list on [02-14-16] makes note of bad paddles exlusively in exp and sim, respective. The list below, is made keepign in mind that when the bad pd cuts are finally applied, if a paddle is removed in exp, it will also have to be removed in sim and viceversa, i.e. the list below is compiled keeping in mind the final cut display that will be required in my thesis. Basically, in the following list, both exp and sim will have a union of the list of paddles that were identfied exclusively for each.
	bad_sc_pd_code_list=[
        [[311,205,505],[311,205,505]], #!E
        [[311,324,520],[311,324,520]], #!P
        [[145,245,311,324,337,338,340,342,345,346,347,434,520,532,540,542,632,633,634,635,636,637,638,639,644],[145,245,311,324,337,338,340,342,345,346,347,434,520,532,540,542,632,633,634,635,636,637,638,639,644]] #! PIP
        ]

	#! Now make bad_sc_pd_list[NPRTCL][DTYP]{pd:[sctrl]}
	bad_sc_pd_list=[[[]for j in range(NDTYP-1)]for i in range(NPRTCL-1)]
	for i in (range(NPRTCL-1)):
    		for j in (range(NDTYP-1)):
        		bad_sc_pd_list[i][j]={}
        		for code in bad_sc_pd_code_list[i][j]:
            			pd,sctr=get_paddle(code),get_sector(code)
            			if pd in bad_sc_pd_list[i][j].keys():
                			#print "here for",PRTCL_NAME[i],DTYP_NAME[j],pd
                			bad_sc_pd_list[i][j][pd].append(sctr)
            			else:
                			bad_sc_pd_list[i][j][pd]=[sctr]
	#! direct print
	#print bad_sc_pd_list
	#! pretty print
	for i in (range(NPRTCL-1)):
    		for j in (range(NDTYP-1)):
        		print "Going to list bad paddles for:",PRTCL_NAME[i],DTYP_NAME[j]
        		for pd in bad_sc_pd_list[i][j].keys():
            			print "pd:sctrs=",pd,":",bad_sc_pd_list[i][j][pd]

	#! Now make good_sc_pd_list[prtcl][dtyp]{pd:[sctrl]}
	#! + This structure contains for a particular bad paddle in a given sector, list of sectors where it is good
	#! + Therefore this list can be used to ascertain the dip in theta_vs_p from these bad paddles 
	good_sc_pd_list=[[[]for j in range(NDTYP-1)]for i in range(NPRTCL-1)]
	for i in (range(NPRTCL-1)):
    		for j in (range(NDTYP-1)):
        		good_sc_pd_list[i][j]={}
        		for pd in bad_sc_pd_list[i][j].keys():
            			#! Assume all sectors are good
            			sctrl=set([1,2,3,4,5,6])
            			#! Now get bad sectors for this paddle,remove them from sctrl and form sctr_goodl
            			sctr_badl=bad_sc_pd_list[i][j][pd]
            			sctr_goodl=sctrl-set(sctr_badl)
            			#! Now add this to good_sc_pd_list[i][j][pd]
            			good_sc_pd_list[i][j][pd]=list(sctr_goodl)

	#! direct print
	#print good_sc_pd_list
	#! pretty print
	for i in (range(NPRTCL-1)):
    		for j in (range(NDTYP-1)):
        		print "Going to list good paddles for:",PRTCL_NAME[i],DTYP_NAME[j]
        		for pd in good_sc_pd_list[i][j].keys():
            			print "pd:sctrs=",pd,":",good_sc_pd_list[i][j][pd]

	#! Now using good_sc_pd_list, make plots for correlated dips in theta_vs_p for bad_sc_pd_list
	#! First create structure to store h
	NPDL=70 #! actually less than 70, but just to be safe
	if Q2Wdomain=="Q2-W": #! Create h_bad_sc_pd[dtyp][q2-w][sctr][prtcl][plt][pd]
		h_bad_sc_pd=[[[[[[[]for n in range(NPDL)]for m in range(NPLT)]for l in range(NPRTCL)]for k in range(NSCTR)]for j in range(len(Q2BIN_LEL)*len(WBIN_LEL))]for i in range(NDTYP)]
	elif Q2Wdomain=="Q2": #! Create h_bad_sc_pd[dtyp][q2][sctr][prtcl][plt]
		h_bad_sc_pd=[[[[[[[]for n in range(NPDL)]for m in range(NPLT)]for l in range(NPRTCL)]for k in range(NSCTR)]for j in range(len(Q2BIN_LEL))]for i in range(NDTYP)]
	 
	#! Now fill h
	for iq2wb,q2wbin_le in enumerate(DLE):
		if DEBUG and iq2wb>0: continue #! debug
		q2min,q2max=q2wbin_le[0],DUE[iq2wb][0]
		if   Q2Wdomain=="Q2-W": wmin,wmax=q2wbin_le[1],DUE[iq2wb][1]
		elif Q2Wdomain=="Q2":   wmin,wmax=WMIN,WMAX
		dl=[range(NDTYP),range(NPRTCL),range(NPLT)]
		d=list(itertools.product(*dl))
		for r in d:
			idtyp,iprtcl,iplt=r[0],r[1],r[2]
			dtyp=DTYP_NAME[idtyp]
			prtcl=PRTCL_NAME[iprtcl]
			plt=PLT_NAME[iplt]

			if dtyp=="ST": continue #! this study not done for ST
			if prtcl=="pim": continue #! this study not done for PIM
			if plt!="theta_vs_p": continue #! this study only uses theta_vs_p

			#print PRTCL_NAME[iprtcl],DTYP_NAME[idtyp]
			#print "good_sc_pd_list[iprtcl][idtyp].keys()=",good_sc_pd_list[iprtcl][idtyp].keys()
			for pd in good_sc_pd_list[iprtcl][idtyp].keys():
				ipd=pd-1

				#! Make TCut
				#! + pd:sectors used here will be where the pd is good
				good_sctrl=good_sc_pd_list[iprtcl][idtyp][pd]
				cut_q2w=ROOT.TCut("Q2>%f && Q2<%f && W>%f && W<%f"%(q2min,q2max,wmin,wmax))
				cut=ROOT.TCut(cut_q2w)
				if dtyp=="ER" or dtyp=="SR": #! not for ST
					cut+=cut_top #! cut_top defined globally
					#! Construct cut_str_good_sctrs and add to cut
					for z,sctr in enumerate(good_sctrl):
						if z==0: cut_str_good_sctrs="sector_%s==%d"%(prtcl,sctr)
						else:    cut_str_good_sctrs+="||sector_%s==%d"%(prtcl,sctr)
					cut_good_sctrs=ROOT.TCut(cut_str_good_sctrs)
					cut+=cut_good_sctrs
					#! Re-code  good_pd_codel and add it to cut
					good_pd_codel=[]
					for sctr in good_sctrl:
						pdcode="%d%02d"%(sctr,pd)
						good_pd_codel.append(pdcode)
					for z,pd_code in enumerate(good_pd_codel):
						if z==0: cut_str_good_pd_codes="sc_pd_%s==%s"%(prtcl,pd_code)
						else:    cut_str_good_pd_codes+="||sc_pd_%s==%s"%(prtcl,pd_code)
					cut_good_pd_codes=ROOT.TCut(cut_str_good_pd_codes)
					cut+=cut_good_pd_codes
			
				#! Now indentify sectors where pd is bad and loop over them to fill h
				#! + isctr for filling h will be used from this list
				bad_sctrl=bad_sc_pd_list[iprtcl][idtyp][pd]
				for isctr in range(NSCTR):
					sctr=isctr+1
					if sctr in bad_sctrl:
						print "TTree::Draw() for q2min,q2max,wmin,wmax,dtyp,prtcl,plt,pd,sctr,good_sctrl,good_pdcodel=",q2min,q2max,wmin,wmax,dtyp,prtcl,plt,pd,sctr,good_sctrl,good_pd_codel
						#print "TCut=",cut.GetTitle()
						if T[idtyp].GetEntries()>0:
							#print "drawing histogram for idtyp,i,isctr,iprtcl,iplt,ipd=",idtyp,i,isctr,iprtcl,iplt,ipd
							T[idtyp].Draw(draw_cmd[isctr][iprtcl][iplt],cut)
							#! Store histogram
							#ROOT.gStyle.SetOptStat("ne")
							htmp=ROOT.gDirectory.Get("hcmd")#(draw_cmd_hst_name[idtyp][isctr][icorr])
							h_bad_sc_pd[idtyp][iq2wb][isctr][iprtcl][iplt][ipd]=htmp.Clone()#"hmm_%s_%s"%(DTYP_NAME[idtyp],CORR_NAME[icorr]))
							h_bad_sc_pd[idtyp][iq2wb][isctr][iprtcl][iplt][ipd].SetName("h_%s_%s_%s_pd%d_s%d"%(dtyp,prtcl,plt,pd,sctr))
							h_bad_sc_pd[idtyp][iq2wb][isctr][iprtcl][iplt][ipd].SetTitle("%s_%s %.2f-%.2f_%.3f-%.3f pd%d s_%d"%(prtcl,plt,q2min,q2max,wmin,wmax,pd,sctr))#("%s"%cut.GetTitle())
							print "histogram entries=",h_bad_sc_pd[idtyp][iq2wb][isctr][iprtcl][iplt][ipd].GetEntries()
						else:
							print "TTree name=",T[idtyp].GetName()
							print "TTree entries=",T[idtyp].GetEntries()
							sys.exit("TTree has no entries. Because of what appears to be a bug in ROOT, this can lead, ina subtle manner, to deeper problems: TTree::Draw() should return a Null pointer in this case, however, it seems that ROOT returns the pointer to the previous histogram draw by TTree::Draw and therefore I end up referencing the wrong histogram. A workaround could be to use unique names for the histograms in the 'draw_cmd'. However, for now, this will have to do. ")
					else:
						print "TTree::Draw() for q2min,q2max,wmin,wmax,dtyp,prtcl,plt,pd,sctr,good_sctrl,good_pdcodel=",q2min,q2max,wmin,wmax,dtyp,prtcl,plt,pd,sctr,good_sctrl,good_pd_codel,":Skip(because sctr in good_sctrl)"
						h_bad_sc_pd[idtyp][iq2wb][isctr][iprtcl][iplt][ipd]=None

	#print "test",  h_bad_sc_pd[ER][0][5][PIP][0][31].GetEntries()
	#print "test2", h_bad_sc_pd[0][0][5][2][0][31].GetEntries()
	#c2=ROOT.TCanvas("c2","c2")
	#h_bad_sc_pd[0][0][5][2][0][31].Draw("colz")
	#c2.SaveAs("/tmp/canvas/testw.png")

	#! Plot on TCanvas and save h:
	fname="bad_sc_pd_theta_vs_p"
	if MAKE_KIN_PLOTS:
		fname="bad_good_sc_pd_theta_vs_p"
	if DEBUG:
		fname+="_dbg"
	fout_root=ROOT.TFile("%s/%s.root"%(outdir_sc_pd_theta_vs_p,fname),"RECREATE")
	#! CWDTH,CHGHT defined as per 3,2 TCanvas
	CWDTH=2000#500
	CHGHT=1500#300
	for iq2wb,q2wbin_le in enumerate(DLE):
		if DEBUG and iq2wb>0: continue #! debug
		q2min,q2max=q2wbin_le[0],DUE[iq2wb][0]
		if   Q2Wdomain=="Q2-W": wmin,wmax=q2wbin_le[1],DUE[iq2wb][1]
		elif Q2Wdomain=="Q2":   wmin,wmax=WMIN,WMAX
		q2wbin="%.2f-%.2f"%(q2min,q2max)
		#outdir_q2w="%s/%s"%(outdir_sc_pd_theta_vs_p,q2wbin) #! for .png
		#if not os.path.exists(outdir_q2w):
                #       os.makedirs(outdir_q2w)
		q2wbindir_root=fout_root.mkdir(q2wbin)#! for .root
		q2wbindir_root.cd()
		dl=[range(NDTYP),range(NPRTCL),range(NPLT)]
		d=list(itertools.product(*dl))
		for r in d:
			idtyp,iprtcl,iplt=r[0],r[1],r[2]
			dtyp=DTYP_NAME[idtyp]
			prtcl=PRTCL_NAME[iprtcl]
			plt=PLT_NAME[iplt]

			if dtyp=="ST": continue #! this study not done for ST
			if prtcl=="pim": continue #! this study not done for pim
			if plt!="theta_vs_p": continue #! this study only uses theta_vs_p
	
			#! Make directory for storing plots from all pds
			pd_all_dir_name="pdall"
			if q2wbindir_root.GetDirectory(pd_all_dir_name)==None:
				pd_all_dir_root=q2wbindir_root.mkdir(pd_all_dir_name)
			else:
                        	pd_all_dir_root=q2wbindir_root.GetDirectory(pd_all_dir_name)
			#! Make canvas for plot for all pd
			cname_pd_all="c_%s_%s_pdall_%s_%s"%(prtcl,plt,dtyp,q2wbin)
                        c_pd_all=ROOT.TCanvas(cname_pd_all,cname_pd_all,CWDTH,CHGHT)
                        c_pd_all.Divide(3,2)
			htmp=[]
			for isctr in range(NSCTR):
				sctr=isctr+1
                                c_pd_all.cd(sctr)
				if MAKE_KIN_PLOTS:
					hkin[idtyp][iq2wb][isctr][iprtcl][iplt].Draw("colz")
				else: #! then make histogram first because then that will provide TF1::Draw the axes to draw on
					hname="%s_%s_s%d"%(plt,prtcl,sctr) 
                                	nbinsx,xmin,xmax=BNG[plt,prtcl][0],BNG[plt,prtcl][1],BNG[plt,prtcl][2]
                               		nbinsy,ymin,ymax=BNG[plt,prtcl][3],BNG[plt,prtcl][4],BNG[plt,prtcl][5]
                               		htmp.append(ROOT.TH2F(hname,hname,nbinsx,xmin,xmax,nbinsy,ymin,ymax))
                                  	htmp[isctr].Draw()

			#! Now loop over paddles
			for pd in good_sc_pd_list[iprtcl][idtyp].keys():
				ipd=pd-1
				pddir_name="pd%d"%pd
				if q2wbindir_root.GetDirectory(pddir_name)==None:
					pddir_root=q2wbindir_root.mkdir(pddir_name)
				else:
					pddir_root=q2wbindir_root.GetDirectory(pddir_name)
				pddir_root.cd()

				cname_pd="c_%s_%s_pd%d_%s_%s"%(prtcl,plt,pd,dtyp,q2wbin)
				c_pd=ROOT.TCanvas(cname_pd,cname_pd,CWDTH,CHGHT)
				c_pd.Divide(3,2)
				for isctr in range(NSCTR):
					sctr=isctr+1
					print "Going to plot for q2min,q2max,wmin,wmax,dtyp,prtcl,plt,pd,sctr",q2min,q2max,wmin,wmax,dtyp,prtcl,plt,pd,sctr
					#! First make plots on canvas for individual paddle
					c_pd.cd(sctr)
					#! First draw h
					if MAKE_KIN_PLOTS:
						hkin[idtyp][iq2wb][isctr][iprtcl][iplt].Draw("colz")
					#! Now superimpose h_bad_sc_pd
					#! This condition should not be met: if len(h_bad_sc_pd[idtyp][iq2wb][isctr][iprtcl][iplt][ipd])==0: continue
					if h_bad_sc_pd[idtyp][iq2wb][isctr][iprtcl][iplt][ipd]!=None: #!good_sc_pd:[isctr] not plotted
						print "h_bad_sc_pd!=None -> Plot"
						if MAKE_KIN_PLOTS: draw_opt="sames"
						else:          draw_opt="colz"
						h_bad_sc_pd[idtyp][iq2wb][isctr][iprtcl][iplt][ipd].Draw(draw_opt)
					else:
						print "h_bad_sc_pd==None -> Skip Plot"
					
					#! Now make plots on canvas for all paddles
					pad=c_pd_all.cd(sctr)
					if MAKE_KIN_PLOTS: #! hkin[idtyp][iq2wb][isctr][iprtcl][iplt] already plotted earlier!
						if h_bad_sc_pd[idtyp][iq2wb][isctr][iprtcl][iplt][ipd]!=None: #!good_sc_pd:[isctr] not plotted
							print "h_bad_sc_pd!=None -> Plot"
                                                	h_bad_sc_pd[idtyp][iq2wb][isctr][iprtcl][iplt][ipd].Draw("sames")
						else:
							print "h_bad_sc_pd==None -> Skip Plot"
                                        else:
						if h_bad_sc_pd[idtyp][iq2wb][isctr][iprtcl][iplt][ipd]!=None: #!good_sc_pd:[isctr] not plotted
							print "h_bad_sc_pd!=None -> Plot"
							if found_hist_on_tpad(pad): draw_opt="colz sames"
							else:                       draw_opt="colz"
							h_bad_sc_pd[idtyp][iq2wb][isctr][iprtcl][iplt][ipd].Draw(draw_opt)
						else:
							print "h_bad_sc_pd==None -> Skip Plot"
				#c.SaveAs("/tmp/canvas/%s.png"%(cname))
				c_pd.Write("%s"%cname_pd,ROOT.gROOT.ProcessLine("TObject::kOverwrite"))#! for .root; latest namecycl
			#! On canvas for all paddles, superimpose EI's theta_vs_p cuts
			for isctr in range(NSCTR):
                        	sctr=isctr+1
                                pad=c_pd_all.cd(sctr)
                                #! Superimpose EI's theta_vs_p cuts
                                draw_EI_theta_vs_p_cuts(iprtcl,isctr,pad)
			pd_all_dir_root.cd()
			c_pd_all.Write("%s"%cname_pd_all,ROOT.gROOT.ProcessLine("TObject::kOverwrite"))#! for .root; latest namecycl
	fout_root.Close()

#! If wanting to keep TCanvas open till program exits				
#if not ROOT.gROOT.IsBatch():
#	plt.show()
#	# wait for you to close the ROOT canvas before exiting
#	wait(True)

#if __name__ == "__main__":	
#	if len(sys.argv)==2:
#		plot_fid(sys.argv[1])
#	elif len(sys.argv)==3:
#		plot_fid(sys.argv[1],int(sys.argv[2])) 
	
