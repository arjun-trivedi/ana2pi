#!/usr/bin/python
from __future__ import division

import os,sys,glob
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait
from collections import OrderedDict

import numpy as np

DATADIR=os.path.join(os.environ['DELASTDIR_SIM'],'study_elast_sim')
OUTDIR=os.path.join(os.environ['STUDY_ELAST_SIM'])

# #! Get simdirl
# simdirl=glob.glob("%s/*"%DATADIR)
# #! Sort simdirl as per f par
# simdirl.sort(key=lambda s: s.split("/")[8].split("_")[2])
# print simdirl

#! Get simdirl in a more orderly way i.e. as per pars
simdirl=[]
F=[0.2,0.5,0.8,1]
D=[0.005,0.01,0.02]
G=[2,3]
STRG=[0,1]
for f in F:
        for d in D:
                for g in G:
                        for strg in STRG:
                                if strg==0: #! since during simulation _strg suffux was not specified when strg=0       
                                        if f!=1 and d!=0.005:
                                                sim="elast_pars_f%.1f_d%.2f_g%d"%(f,d,g)
                                        elif f==1 and d!=0.005:
                                                sim="elast_pars_f%.0f_d%.2f_g%d"%(f,d,g)
                                        elif f!=1 and d==0.005:
                                                sim="elast_pars_f%.1f_d%.3f_g%d"%(f,d,g)
                                        elif f==1 and d==0.005:
                                                sim="elast_pars_f%.0f_d%.3f_g%d"%(f,d,g)
                                        else:
                                            sys.exit("Case unaccounted for f,d,g,strg=",f,d,g,strg)

                                else:
                                        if f!=1 and d!=0.005:
                                                sim="elast_pars_f%.1f_d%.2f_g%d_strg1"%(f,d,g)
                                        elif f==1 and d!=0.005:
                                                sim="elast_pars_f%.0f_d%.2f_g%d_strg1"%(f,d,g)
                                        elif f!=1 and d==0.005:
                                                sim="elast_pars_f%.1f_d%.3f_g%d_strg1"%(f,d,g)
                                        elif f==1 and d==0.005:
                                                sim="elast_pars_f%.0f_d%.3f_g%d_strg1"%(f,d,g)
                                        else:
											sys.exit("Case unaccounted for f,d,g,strg=",f,d,g,strg)
                                simdirl.append('%s/%s'%(DATADIR,sim))
#print len(simdirl)
#print simdirl

#! Remove simdirs that do not exist or do not have delastT,R
simdirdl=[]
for simdir in simdirl:
	#print "checking",simdir 
	#! First check if simdir exists
	if os.path.exists(simdir)==False:
		#print "removing here", simdir
		#simdirl.remove(simdir)
		simdirdl.append(simdir)
		continue
	#! Now check if delastT/R exists
	if os.path.exists("%s/delastT_080415.root"%simdir)==False or os.path.exists("%s/delastR_080415.root"%simdir)==False:
		#print "removing here2",simdir
		#simdirl.remove(simdir)
		simdirdl.append(simdir)
for simdird in simdirdl:
	#print "going to remove",simdird
	simdirl.remove(simdird)

nsim=len(simdirl)
#print nsim

NSEQ=2
T,R=range(2)

#! create fin[nsim][NSEQ]
fin=[[[]for j in range(NSEQ)] for i in range(nsim)]

#!create simname[nsim],clr[nsim]
simname=[[] for i in range(nsim)]
clr=[[] for i in range(nsim)]
#! create hW[nsim][NSEQ],helast[nsim][NSEQ],cnvs[nsim]
hW=[[[]for j in range(NSEQ)] for i in range(nsim)]
helast=[[[]for j in range(NSEQ)] for i in range(nsim)]
cnvs=[[] for i in range(nsim)]
#print hW
#print helast

#! Obtain ER and Isupov data
fin_ER=ROOT.TFile("%s/test_elast_sim_073015/delastR.root"%os.environ['DELASTDIR_EXP'])
fin_ER_W_PID=ROOT.TFile("%s/test_elast_sim_073015/delastR_w_pid.root"%os.environ['DELASTDIR_EXP'])
fin_SR_E16=ROOT.TFile("%s/test_elast_sim_073015/e16/delastR.root"%os.environ['DELASTDIR_SIM'])
fin_ST_E16=ROOT.TFile("%s/test_elast_sim_073015/e16/delastT.root"%os.environ['DELASTDIR_SIM'])

hW_ER=fin_ER.Get("delast/monitor/hW")
helast_ER=fin_ER.Get("delast/monitor/helastic")
hW_ER_W_PID=fin_ER_W_PID.Get("delast/monitor/hW")
helast_ER_W_PID=fin_ER_W_PID.Get("delast/monitor/helastic")
hW_ST_E16=fin_ST_E16.Get("delast/monitor/hW")
helast_ST_E16=fin_ST_E16.Get("delast/monitor/helastic")
hW_SR_E16=fin_SR_E16.Get("delast/monitor/hW")
helast_SR_E16=fin_SR_E16.Get("delast/monitor/helastic")

hW_ER.SetLineColor(ROOT.gROOT.ProcessLine('kBlue'))
helast_ER.SetLineColor(ROOT.gROOT.ProcessLine('kBlue'))

hW_ER_W_PID.SetLineColor(ROOT.gROOT.ProcessLine('kGreen+4'))
helast_ER_W_PID.SetLineColor(ROOT.gROOT.ProcessLine('kGreen+4'))

hW_ST_E16.SetLineColor(ROOT.gROOT.ProcessLine('kBlack'))
helast_ST_E16.SetLineColor(ROOT.gROOT.ProcessLine('kBlack'))

hW_SR_E16.SetLineColor(ROOT.gROOT.ProcessLine('kBlack'))
helast_SR_E16.SetLineColor(ROOT.gROOT.ProcessLine('kBlack'))

#! Now obtain fin, simname,clr,hW, helast
for isim,simdir in enumerate(simdirl):
	pos=simdir.find("_f")
        simname[isim]=simdir[pos+1:]

	fin[isim][T]=ROOT.TFile("%s/delastT_080415.root"%simdir)
	fin[isim][R]=ROOT.TFile("%s/delastR_080415.root"%simdir) 

	hW[isim][T]=fin[isim][T].Get("delast/monitor/hW")
	helast[isim][T]=fin[isim][T].Get("delast/monitor/helastic")

	hW[isim][R]=fin[isim][R].Get("delast/monitor/hW")
        helast[isim][R]=fin[isim][R].Get("delast/monitor/helastic")

	#! hist aesthetic
	hW[isim][T].SetLineColor(ROOT.gROOT.ProcessLine('kRed'))
        helast[isim][T].SetLineColor(ROOT.gROOT.ProcessLine('kRed'))

	hW[isim][R].SetLineColor(ROOT.gROOT.ProcessLine('kRed'))
        helast[isim][R].SetLineColor(ROOT.gROOT.ProcessLine('kRed'))

	#! Make plots
	ROOT.gStyle.SetOptStat("nemriuo")
	ROOT.gStyle.SetStatW(0.35)
	ROOT.gStyle.SetStatH(0.4)
	cnvs[isim]=ROOT.TCanvas("%s"%simname,"%s"%simname)
	cnvs[isim].Divide(2,2)

	cnvs[isim].cd(1)#! W:ST
	hW[isim][T].Draw()
	#! Determine scale factor to scale ST_E16
	#s=helast[isim][T].GetMaximum()/helast_ST_E16.GetMaximum()
	#hW_ST_E16_tmp=hW_ST_E16.Clone()
	#hW_ST_E16_tmp.Scale(s)
	#hW_ST_E16_tmp.Draw("same")

	cnvs[isim].cd(2)#! elast:ST
        helast[isim][T].Draw()
	#! Determine scale factor to scale ST_E16
        #s=helast[isim][T].GetMaximum()/helast_ST_E16.GetMaximum()
        #helast_ST_E16_tmp=helast_ST_E16.Clone()
        #helast_ST_E16_tmp.Scale(s)
        #helast_ST_E16_tmp.Draw("same")

	cnvs[isim].cd(3)#! W:SR
        hW[isim][R].Draw()
	#! Determine scale factor to scale SR_E16
        #s=helast[isim][R].GetMaximum()/helast_SR_E16.GetMaximum()
        #hW_SR_E16_tmp=hW_SR_E16.Clone()
        #hW_SR_E16_tmp.Scale(s)
        #hW_SR_E16_tmp.Draw("same")

	l=ROOT.TLegend(0.15,0.7,0.4,0.9)
        cnvs[isim].cd(4)#! elast:SR and ER!
        helast[isim][R].Draw()
	#! Determine scale factor to scale SR_E16 and ER
        #sSR=helast[isim][R].GetMaximum()/helast_SR_E16.GetMaximum()
	sER=helast[isim][R].GetMaximum()/helast_ER.GetMaximum()
	sER_W_PID=helast[isim][R].GetMaximum()/helast_ER_W_PID.GetMaximum()
        #helast_SR_E16_tmp=helast_SR_E16.Clone()
        #helast_SR_E16_tmp.Scale(sSR)
	helast_ER_tmp=helast_ER.Clone()
        helast_ER_tmp.Scale(sER)
	helast_ER_W_PID_tmp=helast_ER_W_PID.Clone()
        helast_ER_W_PID_tmp.Scale(sER_W_PID)
        #helast_SR_E16_tmp.Draw("same")
	helast_ER_tmp.Draw("same")
	helast_ER_W_PID_tmp.Draw("same")
	l.AddEntry(helast[isim][R],"%s"%simname[isim],"l")
	l.AddEntry(helast_ER_tmp,"ER-n-pid","l")
	l.AddEntry(helast_ER_W_PID_tmp,"ER-w-pid","l")
	l.Draw("same")

	#outdir="/tmp/plots/%s"%simname[isim]
	outdir="%s/direct_comp"%OUTDIR
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	#cnvs[isim].SaveAs("%s/c.png"%outdir)
	cnvs[isim].SaveAs("%s/c%s.png"%(outdir,simname[isim]))
		
	
	#! histogram aesthetics when plotting all on one canvas
        #if isim>=0 and isim<10:
        #        clr[isim]=ROOT.gROOT.ProcessLine("kPink+%d"%isim)
        #elif isim>=10 and isim<20:
        #        clr[isim]=ROOT.gROOT.ProcessLine("kViolet+%d"%(isim-10))
        #elif isim>=20 and isim<30:
        #        clr[isim]=ROOT.gROOT.ProcessLine("kTeal+%d"%(isim-10))
        #elif isim>=30 and isim<40:
        #        clr[isim]=ROOT.gROOT.ProcessLine("kSpring+%d"%(isim-10))
        #elif isim>=40 and isim<50:
        #        clr[isim]=ROOT.gROOT.ProcessLine("kOrange+%d"%(isim-10))

	#hW[isim][T].SetLineColor(clr[isim])
	#helast[isim][T].SetLineColor(clr[isim])

	#hW[isim][R].SetLineColor(clr[isim])
        #helast[isim][R].SetLineColor(clr[isim])
#print simname
#print hW
#print helast

#! Plots (all in one canvas)
#BGCLR=12
#c=ROOT.TCanvas()
#c.Divide(2,2)
##l=ROOT.TLegend(0.6,0.7,0.8,0.8)
#l=ROOT.TLegend(0.1,0.1,0.9,0.9)
#p=c.cd(1)
#for isim,h in enumerate(hW):
#	p.SetFillColor(BGCLR)
#	if isim==0:h[T].Draw()
#	else:h[T].Draw("same")
#	l.AddEntry(h[T],simname[isim],"l")
#
#p=c.cd(2)
#for isim,h in enumerate(helast):
#	p.SetFillColor(BGCLR)
#        if isim==0:h[T].Draw()
#        else:h[T].Draw("same")
#
#p=c.cd(3)
#for isim,h in enumerate(hW):
#	p.SetFillColor(BGCLR)
#        if isim==0:h[R].Draw()
#        else:h[R].Draw("same")
#
#p=c.cd(4)
#for isim,h in enumerate(helast):
#	p.SetFillColor(BGCLR)
#        if isim==0:h[R].Draw()
#        else:h[R].Draw("same")
#
##! Separate canvas to draw the long legend
#c1=ROOT.TCanvas()#("c1","c1",100,1000)
#l.SetFillColor(BGCLR)
#l.Draw()

#! If wanting to keep TCanvas open till program exits				
if not ROOT.gROOT.IsBatch():
	plt.show()
	# wait for you to close the ROOT canvas before exiting
	wait(True)
#
##if __name__ == "__main__":	
##	if len(sys.argv)==2:
##		plot_fid(sys.argv[1])
##	elif len(sys.argv)==3:
##		plot_fid(sys.argv[1],int(sys.argv[2])) 
#	
