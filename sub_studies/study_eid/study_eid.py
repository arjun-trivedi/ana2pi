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

'''
+ {e1f,e16}*{Q2,W}*{dtyp}*{cutlvl}*{plots}, where
	+ Q2,W:
		+ full Q2-W range
		+ Q2 bins
		+ Q2,W bins (if there is enough statistics)
	+ cutlvel:
		+ mon (monitor: i.e. not eid cut applied)
		+ cut (eid cut)
		+ evt (after 2pi event selection)
	+ dtyp:
		+ ER 
		+ SR
	+ plots
		+ EC:
			+ U,V,W
			+ eout_vs_ein 
			+ etot_vs_p
			+ SF_vs_p
		+ CC
			+ theta_vs_seg 
			+ nphe

+ Usage: study_eid.py expt=<e1f/e16> top[=2] debug[=False]
'''
USAGE="study_eid.py expt=<e1f/e16> top[=2] debug[=False]"
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

DEBUG=False
if len(sys.argv)>3: #! i.e. debug entered by user
        if sys.argv[3]=="True":
                DEBUG=True
        elif sys.argv[3]=="False":
                DEBUG=False
        else:
                sys.exit("Please enter debug as True/False only.")
print "DEBUG=",DEBUG
#sys.exit()
#! ****

#! Set up constants
NDTYP=2
ER,SR=range(NDTYP)
DTYP_NAME=["ER","SR"]

NCUTLVL=3
MON,CUT,EVT=range(NCUTLVL)
CUTLVL_NAME=["mon","cut","evt"]

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

#! *** Prepare input data structure: FIN[dtyp],T[dtyp][cutlvl] ***
FIN=[[] for i in range(NDTYP)]
T=[[[]for j in range(NCUTLVL)] for i in range(NDTYP)]
#print "FIN=",FIN
#print "T=",T
#! ***

#! *** Prepare output datadir***
if expt=="e1f":
	DATADIR_OUTPUT=os.environ['STUDY_EID_E1F']
elif expt=="e16":
	DATADIR_OUTPUT=os.environ['STUDY_EID_E16']
print "DATADIR_OUTPUT=",DATADIR_OUTPUT
#! ***

#! *** Now fill T[dtyp][cutlvl] ***
DATADIR_INPUT=[[] for i in range(NDTYP)]
sfx=''
if expt=='e16':sfx='_E16'
DATADIR_INPUT[ER]=os.path.join(os.environ['D2PIDIR_EXP%s'%sfx],"data_eid_020816")
DATADIR_INPUT[SR]=os.path.join(os.environ['D2PIDIR_SIM%s'%sfx],"data_eid_020816")

print "DATADIR_INPUT[ER]=",DATADIR_INPUT[ER]
print "DATADIR_INPUT[SR]=",DATADIR_INPUT[SR]

for idtyp,icutlvl in zip(range(NDTYP),range(NCUTLVL)):
	print "Getting FIN,T for",DTYP_NAME[idtyp]
	FIN[idtyp]=ROOT.TFile("%s/deid.root"%DATADIR_INPUT[idtyp])

	T[idtyp][MON]=FIN[idtyp].Get("eid/monitor/t")
	T[idtyp][CUT]=FIN[idtyp].Get("eid/cut/t")
	T[idtyp][EVT]=FIN[idtyp].Get("eid2/cut/t")
#print FIN[ER].GetName()
#print FIN[SR].GetName()
#print T[SR][MON].GetName()
#print T[SR][CUT].GetName()
#print T[SR][EVT].GetName()
#sys.exit()
#! ***

#! *** Now make plots ***
#! Setup plot information
NPLT=8
ECE0_VS_ECEI,SFO_VS_SFI,SF_VS_P,ECU,ECV,EVW,NPHE,CC_THETA_VS_CC_SEGM=range(NPLT)
PLT_NAME=["ec_eo_vs_ec_ei","sfo_vs_sfi","sf_vs_p","ecU","ecV","ecW","nphe","cc_theta_vs_cc_segm"]
PLT_DIM=[2,2,2,1,1,1,1,2]
#if DEBUG:
	#NPLT=1
	#ECU=range(NPLT)
	#PLT_NAME=["ecU"]
	#NPLT=4
	#ECE0_VS_ECEI,ECU,ECV,EVW=range(NPLT)
	#PLT_NAME=["ec_eo_vs_ec_ei","ecU","ecV","ecW"]

#! Setup NENTRIES to be used in TTree::Draw()
NENTRIES=1000000000
if DEBUG:
	NENTRIES=100000

#! setup binning information for plots
BNG={
("ec_eo_vs_ec_ei"):[150,0,1,150,0,1],
("sfo_vs_sfi"):[150,0,0.5,150,0,0.5],
("sf_vs_p"):[160,0,5,150,0,0.5],
("ecU"):[100,0,450],
("ecV"):[100,0,450],
("ecW"):[100,0,450],
("nphe"):[100,0,300],
("cc_theta_vs_cc_segm"):[20,0,20,100,0,50]
}
#! *** 

#! Create and fill draw_cmd[sctr][nplt]
draw_cmd=[[[]for j in range(NPLT)]for i in range(NSCTR)]
dl=[range(NSCTR),range(NPLT)]
d=list(itertools.product(*dl))
for r in d:
	isctr,iplt=r[0],r[1]
	plt=PLT_NAME[iplt]
	plt_dim=PLT_DIM[iplt]
        #if plt=="ecU" or plt=="ecV" or plt=="ecW":
	if plt_dim==1:
        	nbins,xmin,xmax=BNG[plt][0],BNG[plt][1],BNG[plt][2]
                draw_cmd[isctr][iplt]="%s>>hcmd(%d,%d,%d)"%(plt,nbins,xmin,xmax)	
	else:#! if 2D plots
		nbinsx,xmin,xmax=BNG[plt][0],BNG[plt][1],BNG[plt][2]
                nbinsy,ymin,ymax=BNG[plt][3],BNG[plt][4],BNG[plt][5]
        	if   plt=="ec_eo_vs_ec_ei":      draw_cmd_tmp="ec_eo:ec_ei"
		elif plt=="sfo_vs_sfi":          draw_cmd_tmp="ec_eo/p:ec_ei/p"
		elif plt=="sf_vs_p":             draw_cmd_tmp="etot/p:p"
		elif plt=="cc_theta_vs_cc_segm": draw_cmd_tmp="cc_theta:cc_segm"
		draw_cmd[isctr][iplt]="%s>>hcmd(%f,%f,%f,%f,%f,%f)"%(draw_cmd_tmp,nbinsx,xmin,xmax,nbinsy,ymin,ymax)
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

#! First create structure to store heid
if Q2Wdomain=="Q2-W": #! Create heid[dtyp][q2-w][sctr][cutlvl][nplt]
	heid=[[[[[[]for m in range(NPLT)]for l in range(NCUTLVL)]for k in range(NSCTR)]for j in range(len(Q2BIN_LEL)*len(WBIN_LEL))]for i in range(NDTYP)]
elif Q2Wdomain=="Q2": #! Create heid[dtyp][q2][sctr][cutlvl][nplt]
	heid=[[[[[[]for m in range(NPLT)]for l in range(NCUTLVL)]for k in range(NSCTR)]for j in range(len(Q2BIN_LEL))]for i in range(NDTYP)]
#print heid
#sys.exit()

#! Now fill heid
for iq2wb,q2wbin_le in enumerate(DLE):
	if DEBUG and iq2wb>0: continue #! debug
	q2min,q2max=q2wbin_le[0],DUE[iq2wb][0]
	if   Q2Wdomain=="Q2-W": wmin,wmax=q2wbin_le[1],DUE[iq2wb][1]
	elif Q2Wdomain=="Q2":   wmin,wmax=WMIN,WMAX
	dl=[range(NSCTR),range(NDTYP),range(NCUTLVL),range(NPLT)]
	d=list(itertools.product(*dl))
	for r in d:
		isctr,idtyp,icutlvl,iplt=r[0],r[1],r[2],r[3]
		sctr=isctr+1		
		dtyp=DTYP_NAME[idtyp]
		cutlvl=CUTLVL_NAME[icutlvl]
		plt=PLT_NAME[iplt]
		#! Make TCut: begin by creating with Q2,W cut
		cut_q2w=ROOT.TCut("Q2>%f && Q2<%f && W>%f && W<%f"%(q2min,q2max,wmin,wmax))
		cut=ROOT.TCut(cut_q2w)
		#! Add cut_sector
		cut_sctr=ROOT.TCut("sector==%d"%sctr)
                cut+=cut_sctr
		#! Add cut_top defined globally
		if cutlvl=="evt":
			cut+=cut_top #! cut_top defined earlier
		#! Now plot
		print "TTree::Draw() for q2min,q2max,wmin,wmax,sctr,dtyp,cutlvl,plt=",q2min,q2max,wmin,wmax,sctr,dtyp,cutlvl,plt
		if T[idtyp][icutlvl].GetEntries()>0:
			print "draw_cmd=",draw_cmd[isctr][iplt]
			print "cut=",cut
			T[idtyp][icutlvl].Draw(draw_cmd[isctr][iplt],cut,"",NENTRIES)
			#! Store histogram
			#ROOT.gStyle.SetOptStat("ne")
			htmp=ROOT.gDirectory.Get("hcmd")
			heid[idtyp][iq2wb][isctr][icutlvl][iplt]=htmp.Clone()
			heid[idtyp][iq2wb][isctr][icutlvl][iplt].SetName("h_%s_%s_%s_s%d"%(dtyp,cutlvl,plt,sctr))
			heid[idtyp][iq2wb][isctr][icutlvl][iplt].SetTitle("%s_%s %.2f-%.2f_%.3f-%.3f"%(cutlvl,plt,q2min,q2max,wmin,wmax))
		else:
			print "TTree name=",T[idtyp][icutlvl].GetName()
			print "TTree entries=",T[idtyp][icutlvl].GetEntries()
			sys.exit("TTree has no entries. Because of what appears to be a bug in ROOT, this can lead, ina subtle manner, to deeper problems: TTree::Draw() should return a Null pointer in this case, however, it seems that ROOT returns the pointer to the previous histogram draw by TTree::Draw and therefore I end up referencing the wrong histogram. A workaround could be to use unique names for the histograms in the 'draw_cmd'. However, for now, this will have to do. ")

#sys.exit()

#! Plot on TCanvas and save h:
#! + Directly as .png
#! + in .root file
fname="eid"
if DEBUG:
	fname+="_dbg"
fout_root=ROOT.TFile("%s/%s.root"%(outdir,fname),"RECREATE")

#! CWDTH,CHGHT defined as per 3,2 TCanvas
CWDTH=500
CHGHT=300
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
	dl=[range(NDTYP),range(NCUTLVL),range(NPLT)]
	d=list(itertools.product(*dl))
	for r in d:
		idtyp,icutlvl,iplt=r[0],r[1],r[2]
		cutlvl=CUTLVL_NAME[icutlvl]
		plt=PLT_NAME[iplt]
		plt_dim=PLT_DIM[iplt]
		print "Going to plot for [%.2f-%.2f_%.3f-%.3f],%s,%s"%(q2min,q2max,wmin,wmax,cutlvl,plt)
		if plt_dim==1:#plt=="ecU" or plt=="ecV" or plt=="ecW":#! then draw on same canvas
			draw_opt=""
			if plt=="ecU" or plt=="ecV" or plt=="ecW": draw_opt="hist"
			cname="c_%s_%s_%s"%(cutlvl,plt,q2wbin)
			c=ROOT.TCanvas(cname,cname,CWDTH,CHGHT)
			c.Divide(3,2)
			for isctr in range(NSCTR):
				#l=ROOT.TLegend(0.1,0.3,0.3,0.4)#,"","NDC");
				#! Draw hists
				c.cd(isctr+1)
				heid[ER][iq2wb][isctr][icutlvl][iplt].SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
				heid[SR][iq2wb][isctr][icutlvl][iplt].SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
				heid[ER][iq2wb][isctr][icutlvl][iplt].Draw()
				heid[ER][iq2wb][isctr][icutlvl][iplt].Draw("sames %s"%draw_opt) #! Note only "sames" seems to work here (with "same" statbox with name=hxx is drawn?!)
				#! obtain max_ER(=max height of distribution) used to scale SR,ST 
				max_ER=heid[ER][iq2wb][isctr][icutlvl][iplt].GetMaximum()
				if max_ER==0:max_ER=1
				#! scale SR and draw
				max_SR=heid[SR][iq2wb][isctr][icutlvl][iplt].GetMaximum()
				if max_SR==0:max_SR=1
				scl_fctr_SR=max_ER/max_SR
				if scl_fctr_SR==0:scl_fctr_SR=1
				heid[SR][iq2wb][isctr][icutlvl][iplt].Scale(scl_fctr_SR)
				heid[SR][iq2wb][isctr][icutlvl][iplt].Draw("sames %s"%draw_opt)
						
				#! Adjust statbox
				#c.Update()
				ROOT.gPad.Update()
				#! statbox for ER
				pt_ER=heid[ER][iq2wb][isctr][icutlvl][iplt].GetListOfFunctions().FindObject("stats")
				pt_ER.SetTextColor(ROOT.gROOT.ProcessLine("kBlue"))
				#! statbox for SR
				pt_SR=heid[SR][iq2wb][isctr][icutlvl][iplt].GetListOfFunctions().FindObject("stats")
				pt_SR.SetTextColor(ROOT.gROOT.ProcessLine("kRed"))
				diff=pt_ER.GetY2NDC()-pt_ER.GetY1NDC()
				y2=pt_ER.GetY1NDC()
				y1=y2-diff
				pt_SR.SetY1NDC(y1)
				pt_SR.SetY2NDC(y2)
				#c.Update()
				ROOT.gPad.Update()
				#! Draw legend (not drawing since relevant information in statbox)
				#l.AddEntry(hmm[ER][i][icutlvl][iplt],'ER')
				#l.AddEntry(hmm[ER][i][icutlvl][iplt],'SR')
				#l.Draw()
				#c.SaveAs("%s/%s.png"%(outdir_q2w,cname))#! for .png
				c.Write("%s"%cname,ROOT.gROOT.ProcessLine("TObject::kOverwrite"))#! for .root; latest namecycle
		else: #! if 2D plots
		#elif plt=="ec_eo_vs_ec_ei": #! draw ER and ER on separate canvases
			cname="c_%s_%s_ER_%s"%(cutlvl,plt,q2wbin)
			c=ROOT.TCanvas(cname,cname,CWDTH,CHGHT)
			c.Divide(3,2)
			for isctr in range(NSCTR):
				c.cd(isctr+1)
				heid[ER][iq2wb][isctr][icutlvl][iplt].Draw("colz")
			#c.SaveAs("%s/%s.png"%(outdir_q2w,cname))
			c.Write("%s"%cname,ROOT.gROOT.ProcessLine("TObject::kOverwrite"))#! for .root; latest namecycle

			cname="c_%s_%s_SR_%s"%(cutlvl,plt,q2wbin)
			c=ROOT.TCanvas(cname,cname,CWDTH,CHGHT)
			c.Divide(3,2)
			for isctr in range(NSCTR):
				c.cd(isctr+1)
				heid[SR][iq2wb][isctr][icutlvl][iplt].Draw("colz")
			#c.SaveAs("%s/%s.png"%(outdir_q2w,cname))
			c.Write("%s"%cname,ROOT.gROOT.ProcessLine("TObject::kOverwrite"))#! for .root; latest namecycle

fout_root.Close()					
#! ***

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
	
