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

+ Usage: study_det_ineff.py expt=<e1f/e16> top=[2]
'''

if len(sys.argv)<2:
		sys.exit('usage: study_det_ineff.py expt=(e1f/e16) top[=2]')
expt=sys.argv[1]
if expt!='e1f' and expt!='e16':
	sys.exit('expt=%s is not valid. usage: study_det_ineff.py expt=(e1f/e16)'%expt)
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
#sys.exit()
#! *** Prepare input data structure: FIN/T[dtyp][rctn] ***
NDTYP=3
ER,SR,ST=range(NDTYP)
DTYP_NAME=["ER","SR","ST"]

NRCTN=2
TPI,ELAST=range(NRCTN)
RCTN_NAME=["2pi","elast"]

FIN=[[[] for j in range(NRCTN)] for i in range(NDTYP)]
T=[[[] for j in range(NRCTN)] for i in range(NDTYP)]
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

d=[[ER,SR,ST],[TPI,ELAST]]
D=list(itertools.product(*d))
for r in D:
	idtyp,irctn=r[0],r[1]
	if RCTN_NAME[irctn]=="elast": continue #! No data made from Elastic

	print "Getting T for",DTYP_NAME[idtyp],RCTN_NAME[irctn]
	FIN[idtyp][irctn]=ROOT.TFile("%s/d%s%s_ddetineff.root"%(DATADIR_INPUT[idtyp],RCTN_NAME[irctn],DTYP_NAME[idtyp]))

	if RCTN_NAME[irctn]=="elast": 
		T[idtyp][irctn]=FIN[idtyp][irctn].Get("delast/monitor/t")
	elif RCTN_NAME[irctn]=="2pi": 
		if DTYP_NAME[idtyp]=="ER" or DTYP_NAME[idtyp]=="SR":
			T[idtyp][irctn]=FIN[idtyp][irctn].Get("d2piR/tR")
		elif DTYP_NAME[idtyp]=="ST":
			T[idtyp][irctn]=FIN[idtyp][irctn].Get("d2piT/tT")
	#print FIN[idtyp][irctn][icorr].GetName()
	#print T[idtyp][irctn][icorr].GetName()
#print "T=",T
#print FIN[ER][TPI].GetName()
#print FIN[SR][TPI].GetName()
#print FIN[ST][TPI].GetName()
#print T[ER][TPI].GetName()
#print T[ER][TPI].GetEntries()
#print T[SR][TPI].GetName()
#print T[SR][TPI].GetEntries()
#print T[ST][TPI].GetName()
#print T[ST][TPI].GetEntries()
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

#! Setup sector information 
NSCTR=6

#! Setup particle information
NPRTCL=4
E,P,PIP,PIM=range(NPRTCL)
PRTCL_NAME=["e","p","pip","pim"]

#! Setup "plot" information"
NPLT=5
THETA,P,THETA_VS_P,THETA_VS_PHI,SC_PD=range(NPLT)
PLT_NAME=["theta","p","theta_vs_p","theta_vs_phi","sc_pd"]
#! debug
#NPLT=1
#THETA=range(NPLT)
#PLT_NAME=["theta"]
#SC_PD=range(NPLT)
#PLT_NAME=["sc_pd"]

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


Q2Wdomain="Q2" #! should be made a user input
#! + Set up appropriate Q2,W domain list
#! + Set up appropriate hh[dtyp]["q2-w"][sctr][prtcl][nplt]
if Q2Wdomain=="Q2-W":
        dle=[Q2BIN_LEL,WBIN_LEL]
        due=[Q2BIN_UEL,WBIN_UEL]
	#! Create h[dtyp][q2-w][sctr][prtcl][nplt]
	h=[[[[[[]for m in range(NPLT)]for l in range(NPRTCL)]for k in range(NSCTR)]for j in range(len(Q2BIN_LEL)*len(WBIN_LEL))]for i in range(NDTYP)]
elif Q2Wdomain=="Q2":
        dle=[Q2BIN_LEL]
        due=[Q2BIN_UEL]
	#! Create h[dtyp][q2][sctr][prtcl][nplt]
	h=[[[[[[]for m in range(NPLT)]for l in range(NPRTCL)]for k in range(NSCTR)]for j in range(len(Q2BIN_LEL))]for i in range(NDTYP)]
DLE=list(itertools.product(*dle))
DUE=list(itertools.product(*due))

#! Create and fill draw_cmd[sctr][prtcl][nplt]
draw_cmd=[[[[]for k in range(NPLT)]for j in range(NPRTCL)]for i in range(NSCTR)]
for isctr in range(NSCTR):
	for iprtcl in range(NPRTCL):
		for iplt in range(NPLT):
			plt=PLT_NAME[iplt]
			prtcl=PRTCL_NAME[iprtcl]
			if plt=="theta" or plt=="p" or plt=="sc_pd":
				if plt=="sc_pd":
					nbins,xmin,xmax=BNG[plt,prtcl][isctr][0],BNG[plt,prtcl][isctr][1],BNG[plt,prtcl][isctr][2]
				else:
					nbins,xmin,xmax=BNG[plt,prtcl][0],BNG[plt,prtcl][1],BNG[plt,prtcl][2]
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

#! Define any TCuts
#! top==2
#TOP=2
cut_top=ROOT.TCut("top==%d"%TOP)

#! outdir
outdir=os.path.join(DATADIR_OUTPUT,"results_top%d"%TOP)
if not os.path.exists(outdir):
	os.makedirs(outdir)

STDY_BAD_SC_PD_THETA_VS_P_CRLTN=True#! debug

if not STDY_BAD_SC_PD_THETA_VS_P_CRLTN:
	#! fill h
	for i,q2wbin_le in enumerate(DLE):
		#if i>0: continue #! debug
		q2min,q2max=q2wbin_le[0],DUE[i][0]
		if Q2Wdomain=="Q2-W":
			wmin,wmax=q2wbin_le[1],DUE[i][1]
		elif Q2Wdomain=="Q2":
			wmin,wmax=WMIN,WMAX
		for isctr in range(NSCTR):
			for r in D:
				idtyp,irctn=r[0],r[1]
				if RCTN_NAME[irctn]=="elast": continue

				for iprtcl in range(NPRTCL):
					for iplt in range(NPLT):
						prtcl=PRTCL_NAME[iprtcl]
						plt=PLT_NAME[iplt]
						#! Make TCut
						cut_q2w=ROOT.TCut("Q2>%f && Q2<%f && W>%f && W<%f"%(q2min,q2max,wmin,wmax))
						cut=ROOT.TCut(cut_q2w)
						if DTYP_NAME[idtyp]=="ER" or DTYP_NAME[idtyp]=="SR": #! not for ST
							cut+=cut_top #! cut_top defined earlier
							cut_sctr=ROOT.TCut("sector_%s==%d"%(prtcl,isctr+1))
							cut+=cut_sctr
						print "TTree::Draw() for q2min,q2max,wmin,wmax,sctr,dtyp,rctn,prtcl,plt=",q2min,q2max,wmin,wmax,isctr+1,DTYP_NAME[idtyp],RCTN_NAME[irctn],prtcl,plt
						if T[idtyp][irctn].GetEntries()>0:
							T[idtyp][irctn].Draw(draw_cmd[isctr][iprtcl][iplt],cut)#[idtyp][isctr][icorr],cut)

							#! Store histogram
							#ROOT.gStyle.SetOptStat("ne")
							htmp=ROOT.gDirectory.Get("hcmd")#(draw_cmd_hst_name[idtyp][isctr][icorr])
							h[idtyp][i][isctr][iprtcl][iplt]=htmp.Clone()#"hmm_%s_%s"%(DTYP_NAME[idtyp],CORR_NAME[icorr]))
							h[idtyp][i][isctr][iprtcl][iplt].SetName("h_%s_%s_%s"%(DTYP_NAME[idtyp],prtcl,plt))
							#h[idtyp][i][isctr][iprtcl][iplt].SetTitle("%s_%s %.2f-%.2f_%.3f-%.3f"%(prtcl,plt,q2min,q2max,wmin,wmax))#("%s"%cut.GetTitle())
						else:
							print "TTree name=",T[idtyp][irctn].GetName()
							print "TTree entries=",T[idtyp][irctn].GetEntries()
							sys.exit("TTree has no entries. Because of what appears to be a bug in ROOT, this can lead, ina subtle manner, to deeper problems: TTree::Draw() should return a Null pointer in this case, however, it seems that ROOT returns the pointer to the previous histogram draw by TTree::Draw and therefore I end up referencing the wrong histogram. A workaround could be to use unique names for the histograms in the 'draw_cmd'. However, for now, this will have to do. ")
				#print hmm[idtyp][i][j][icorr].GetName()
	#! Plot on TCanvas and save h:
	#! + Directly as .png
	#! + in .root file
	fout_root=ROOT.TFile("%s/results.root"%outdir,"RECREATE")

	#! CWDTH,CHGHT defined as per 3,2 TCanvas
	CWDTH=1000
	CHGHT=800
	for i,q2wbin_le in enumerate(DLE):
		#if i>0: continue #! debug
		q2min,q2max=q2wbin_le[0],DUE[i][0]
		if Q2Wdomain=="Q2-W":
			wmin,wmax=q2wbin_le[1],DUE[i][1]
		elif Q2Wdomain=="Q2":
			wmin,wmax=WMIN,WMAX
		q2wbin="%.2f-%.2f"%(q2min,q2max)
		outdir_q2w="%s/%s"%(outdir,q2wbin) #! for .png
		q2wbindir_root=fout_root.mkdir(q2wbin)#! for .root
		q2wbindir_root.cd()
		if not os.path.exists(outdir_q2w):
			os.makedirs(outdir_q2w)
		for r in D:
			idtyp,irctn=r[0],r[1]
			if RCTN_NAME[irctn]=="elast": continue
			for iprtcl in range(NPRTCL):
				for iplt in range(NPLT):
					prtcl=PRTCL_NAME[iprtcl]
					plt=PLT_NAME[iplt]
		
					print "Going to plot for [%.2f-%.2f_%.3f-%.3f],%s,%s"%(q2min,q2max,wmin,wmax,prtcl,plt)
				
					if plt=="theta" or plt=="p" or plt=="sc_pd":#! then draw on same canvas
						draw_opt=""
						if plt=="sc_pd": draw_opt="hist"
						cname="c_%s_%s"%(prtcl,plt)
						c=ROOT.TCanvas(cname,cname,CWDTH,CHGHT)
						c.Divide(3,2)
						for isctr in range(NSCTR):
							#l=ROOT.TLegend(0.1,0.3,0.3,0.4)#,"","NDC");
							#! Draw hists
							c.cd(isctr+1)
							h[ER][i][isctr][iprtcl][iplt].SetLineColor(ROOT.gROOT.ProcessLine("kBlue"))
							h[SR][i][isctr][iprtcl][iplt].SetLineColor(ROOT.gROOT.ProcessLine("kRed"))
							h[ST][i][isctr][iprtcl][iplt].SetLineColor(ROOT.gROOT.ProcessLine("kGreen"))
							h[ER][i][isctr][iprtcl][iplt].Draw()
							h[ER][i][isctr][iprtcl][iplt].Draw("sames %s"%draw_opt) #! Note only "sames" seems to work here (with "same" statbox with name=hxx is drawn?!)
							#! obtain max_ER(=max height of distribution) used to scale SR,ST 
							max_ER=h[ER][i][isctr][iprtcl][iplt].GetMaximum()
							if max_ER==0:max_ER=1
							#! scale SR and draw
							max_SR=h[SR][i][isctr][iprtcl][iplt].GetMaximum()
							if max_SR==0:max_SR=1
							scl_fctr_SR=max_ER/max_SR
							if scl_fctr_SR==0:scl_fctr_SR=1
							h[SR][i][isctr][iprtcl][iplt].Scale(scl_fctr_SR)
							h[SR][i][isctr][iprtcl][iplt].Draw("sames %s"%draw_opt)
							#! scale ST and draw
							if plt!="sc_pd":
								max_ST=h[ST][i][isctr][iprtcl][iplt].GetMaximum()
								if max_ST==0:max_ST=1
								scl_fctr_ST=max_ER/max_ST
								if scl_fctr_ST==0:scl_fctr_ST=1
								h[ST][i][isctr][iprtcl][iplt].Scale(scl_fctr_ST)
								h[ST][i][isctr][iprtcl][iplt].Draw("sames %s"%draw_opt)
							
							#! Adjust statbox
							#c.Update()
							ROOT.gPad.Update()
							#! statbox for ER
							pt_ER=h[ER][i][isctr][iprtcl][iplt].GetListOfFunctions().FindObject("stats")
							pt_ER.SetTextColor(ROOT.gROOT.ProcessLine("kBlue"))
							#! statbox for SR
							pt_SR=h[SR][i][isctr][iprtcl][iplt].GetListOfFunctions().FindObject("stats")
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
								pt_ST=h[ST][i][isctr][iprtcl][iplt].GetListOfFunctions().FindObject("stats")
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
						c.SaveAs("%s/%s.png"%(outdir_q2w,cname))#! for .png
						c.Write("%s_%s"%(cname,q2wbin),ROOT.gROOT.ProcessLine("TObject::kOverwrite"))#! for .root; latest namecycle
					elif plt=="theta_vs_p" or plt=="theta_vs_phi": #! draw ER and ER on separate canvases
						cname="c_%s_%s_ER"%(prtcl,plt)
						c=ROOT.TCanvas(cname,cname,CWDTH,CHGHT)
						c.Divide(3,2)
						for isctr in range(NSCTR):
							c.cd(isctr+1)
							h[ER][i][isctr][iprtcl][iplt].Draw("colz")
						c.SaveAs("%s/%s.png"%(outdir_q2w,cname))
						c.Write("%s_%s"%(cname,q2wbin),ROOT.gROOT.ProcessLine("TObject::kOverwrite"))#! for .root; latest namecycle

						cname="c_%s_%s_SR"%(prtcl,plt)
						c=ROOT.TCanvas(cname,cname,CWDTH,CHGHT)
						c.Divide(3,2)
						for isctr in range(NSCTR):
							c.cd(isctr+1)
							h[SR][i][isctr][iprtcl][iplt].Draw("colz")
						c.SaveAs("%s/%s.png"%(outdir_q2w,cname))
						c.Write("%s_%s"%(cname,q2wbin),ROOT.gROOT.ProcessLine("TObject::kOverwrite"))#! for .root; latest namecycle
						cname="c_%s_%s_ST"%(prtcl,plt)
						c=ROOT.TCanvas(cname,cname,CWDTH,CHGHT)
						c.Divide(3,2)
						for isctr in range(NSCTR):
							c.cd(isctr+1)
							h[ST][i][isctr][iprtcl][iplt].Draw("colz")
						c.SaveAs("%s/%s.png"%(outdir_q2w,cname)) #! for .png
						c.Write("%s_%s"%(cname,q2wbin),ROOT.gROOT.ProcessLine("TObject::kOverwrite"))#! for .root; latest namecycle
					
#! ***
#! Now study bad_sc_pd and theta_vs_p dip correlation
if STDY_BAD_SC_PD_THETA_VS_P_CRLTN:

	#! Define some functions to sector and paddle from sc_pd_code
	def get_sector(sc_pd_code):
    		return (sc_pd_code-(sc_pd_code%100))/100
	def get_paddle(sc_pd_code):
    		return sc_pd_code%100 

	#! Make list of bad_sc_pd_code[NPRTCL][NDTYP]; NPRTCL=(E,P,PIP),NDTYP=(ER,SR)
	#! + Truly, if a paddle is bad, it is bad for for all particles. 
	#! + However, given the kinematics of each particle, every paddle is not accessible to each particle.
	#! + Therefore, a comprehensive list of bad pad paddles is identified by looking at each particle separately 
	bad_sc_pd_code_list=[
	[[],[311,505]], #!E
	[[],[311,324,520]], #!P
	[[145,340,342,345,346,347,632,634,635,636,637,638,639],[145,311,324,336,337,340,342,434,520,532,540,542,636,637,644]] #! PIP
	]

	#! Now make bad_sc_pd_list[NPRTCL][DTYP]{pd:[sctrl]}
	bad_sc_pd_list=[[[]for j in range(NDTYP-1)]for i in range(NPRTCL-1)]
	for i in (E,P,PIP):
    		for j in (ER,SR):
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
	for i in (E,P,PIP):
    		for j in (ER,SR):
        		print "Going to list bad paddles for:",PRTCL_NAME[i],DTYP_NAME[j]
        		for pd in bad_sc_pd_list[i][j].keys():
            			print "pd:sctrs=",pd,":",bad_sc_pd_list[i][j][pd]

#! Transform bad_sc_pd_code_list[NPRTCL][NDTYP] -> bad_sc_pd_list[NPRTCL][NDTYP][NSCTR]

#! Now, using bad_sc_pd_code_list, make good_sc_pd_code_list, 
#! which will be used to get an idea for the theta_vs_p distributions corresponding to the bad_sc_pd

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
	
