#!/usr/bin/python
from __future__ import division
import os,sys,datetime
from collections import OrderedDict
import array

import ROOT

from rootpy.io import root_open, DoesNotExist

import math

import numpy as np
import matplotlib.pyplot as plt

import itertools

'''

'''

USAGE='study_obs_as_function_of_simstats dbg[=False] show_rel_err_dist[=False]'

#! user inputs
DBG=False
if len(sys.argv)>1: #! dbg entered by user
	if    sys.argv[1]=="True":  DBG=True
	elif  sys.argv[1]=="False": DBG=False
	else: sys.exit('DBG=%s is not valid. usage: %s'%(sys.argv[1],USAGE))

SHOW_REL_ERR_DIST=False
if len(sys.argv)>2: #!  show_rel_err_dist entered by user
	if    sys.argv[2]=="True":  SHOW_REL_ERR_DIST=True
	elif  sys.argv[2]=="False": SHOW_REL_ERR_DIST=False
	else: sys.exit('SHOW_REL_ERR_DIST=%s is not valid. usage: %s'%(sys.argv[3],USAGE))

print "DBG=",DBG
print "SHOW_REL_ERR_DIST=",SHOW_REL_ERR_DIST
#sys.exit()

#! imports from proc_h8.py
sys.path.insert(0, '%s/elast_lite/obs_2pi'%os.environ['ANA2PI'])
from proc_h8 import H5_DIM
H4_PROJDIM=array.array('i',[H5_DIM['M1'],H5_DIM['M2'],H5_DIM['PHI'],H5_DIM['ALPHA']])
from disp_obs import VAR_NAMES_PLAIN

#! Set up THnTool
ROOT.gROOT.ProcessLine(".L THnTool.C+")
from ROOT import THnTool
thntool=THnTool()

#! Set up general ROOT plotting aesthetics
#ROOT.gStyle.SetOptStat("nmMrReiuo")
#ROOT.gStyle.SetStatStyle(0) #! transparent stats box

#! OUTDIR
DATE=datetime.datetime.now().strftime('%m%d%y')
OUTDIRNAME='results'
#! append identifiers
if SHOW_REL_ERR_DIST==True:
	OUTDIRNAME+='_w_relerr_dist'
else:
	OUTDIRNAME+='_wo_relerr_dist'
OUTDIRNAME+='_%s'%DATE
#! Finally create OUTDIR
if DBG==True:
	OUTDIR=os.path.join(os.environ['STUDY_OBS_AS_FUNCTION_OF_SIMSTATS_DATADIR'],'dbg',OUTDIRNAME)
else:
	OUTDIR=os.path.join(os.environ['STUDY_OBS_AS_FUNCTION_OF_SIMSTATS_DATADIR'],OUTDIRNAME)
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)
print "OUTDIR=",OUTDIR
#sys.exit()

#! Setup input root files: FYLD[q2r][sim], FO1D[q2r][sim],FOIT[q2r][sim]
#! Setup Q2 ranges
NQ2RANGES=2
LQ2,HQ2=range(NQ2RANGES)
#! Setup simulations
ISIM=[0 for i in range(NQ2RANGES)]
ISIM[LQ2] =[0,1,2,3,4,5]
ISIM[HQ2]=[0,1,2,3]
SIMNAME=[0 for i in range(NQ2RANGES)]
SIMNAME[LQ2]=['sim4','sim4_sim5', 'sim4_sim5_sim6',  'sim4_sim5_sim6_sim7','sim4_sim5_sim6_sim7_sim8','sim4_sim5_sim6_sim7_sim8_sim13']
SIMNAME[HQ2]=['sim9','sim9_sim10','sim9_sim10_sim11','sim9_sim10_sim11_sim12']
NSIM=[0 for i in range(NQ2RANGES)]
NSIM[LQ2]=len(ISIM[LQ2])
NSIM[HQ2]=len(ISIM[HQ2])
SIMPCNTG=[0 for i in range(NQ2RANGES)]
SIMPCNTG[LQ2]=[16.7,33.3,50.0,66.7,83.3,100.0]
SIMPCNTG[HQ2]=[25,50,75,100]
SIMFRCTN=[0 for i in range(NQ2RANGES)]
SIMFRCTN[LQ2]=[0.17,0.33,0.50,0.67,0.83,1.00]
SIMFRCTN[HQ2]=[0.25,0.50,0.75,1.00]
SIMTOT=[0 for i in range(NQ2RANGES)]
SIMTOT[LQ2]=['~6B']
SIMTOT[HQ2]=['~4B']
#! Finally setup and fill FYLD[q2r][sim],FO1D[q2r][sim],FOIT[q2r][sim]
FYLD=[[0 for i in range(NSIM[LQ2])],[0 for i in range(NSIM[HQ2])]]
FO1D=[[0 for i in range(NSIM[LQ2])],[0 for i in range(NSIM[HQ2])]]
FOIT=[[0 for i in range(NSIM[LQ2])],[0 for i in range(NSIM[HQ2])]]
#print FYLD,FYLD[LQ2],FYLD[HQ2]
#sys.exit()
#tmp_obsdir='%s/tmp/obsdir'%os.environ['GDRIVE'] #! os.environ['OBSDIR_E16'] ALSO -> !cutsncor-cutncors for lowQ2 only!
for q2 in range(NQ2RANGES):
	if q2==LQ2:
		for isim in ISIM[LQ2]:
			#if isim==0: continue #! data corrurpt. remaking. skip for now
			if isim==0: date='102717'
			else:       date='102617'
			FYLD[q2][isim]=root_open('%s/lowQ2_SSBands_off_off_%s_%s/cutsncors1/%s/yield.root'%(os.environ['OBSDIR_E16'],SIMNAME[q2][isim],date,SIMNAME[q2][isim]),'r')
			FO1D[q2][isim]=root_open('%s/lowQ2_SSBands_off_off_%s_%s/cutsncors1/%s/Obs_1D_norm/obs_1D.root'%(os.environ['OBSDIR_E16'],SIMNAME[q2][isim],date,SIMNAME[q2][isim]),'r')
			FOIT[q2][isim]=root_open('%s/lowQ2_SSBands_off_off_%s_%s/cutsncors1/%s/Obs_Itg_Yld_norm/obs_itg_yld.root'%(os.environ['OBSDIR_E16'],SIMNAME[q2][isim],date,SIMNAME[q2][isim]),'r')
	elif q2==HQ2:
		for isim in ISIM[HQ2]:
			if isim==0: date='102617'
			else:       date='102717'
			FYLD[q2][isim]=root_open('%s/highQ2_SSBands_off_off_%s_%s/cutsncors1/%s/yield.root'%(os.environ['OBSDIR_E16'],SIMNAME[q2][isim],date,SIMNAME[q2][isim]),'r')
			FO1D[q2][isim]=root_open('%s/highQ2_SSBands_off_off_%s_%s/cutsncors1/%s/Obs_1D_norm/obs_1D.root'%(os.environ['OBSDIR_E16'],SIMNAME[q2][isim],date,SIMNAME[q2][isim]),'r')
			FOIT[q2][isim]=root_open('%s/highQ2_SSBands_off_off_%s_%s/cutsncors1/%s/Obs_Itg_Yld_norm/obs_itg_yld.root'%(os.environ['OBSDIR_E16'],SIMNAME[q2][isim],date,SIMNAME[q2][isim]),'r')
if DBG==True:
	print FYLD
	print FO1D
	print FOIT
#sys.exit()

#! Set up other analysis constants
#! TCanvas's pad_map[pad,vst,var] defined as per Gleb's display
PAD_MAP_1D=[(1,1,"M1"),   (2,3,'M2'),   (3,2,'M2'),
		    (4,1,"THETA"),(5,3,'THETA'),(6,2,'THETA'),
		    (7,1,"ALPHA"),(8,3,'ALPHA'),(9,2,'ALPHA')]

#! Marker colors as per simulation stats: hot -> cold = low -> high
CLRS=[ROOT.gROOT.ProcessLine("kBlue"),
      ROOT.gROOT.ProcessLine("kCyan"),
      ROOT.gROOT.ProcessLine("kGreen"),
      ROOT.gROOT.ProcessLine("kYellow"),
      #ROOT.gROOT.ProcessLine("kPink+1"),
      ROOT.gROOT.ProcessLine("kOrange"),
      ROOT.gROOT.ProcessLine("kRed")]

#! Marker colors as per simulation stats: hot -> cold = low -> high
CLRS_MPLT=['b',
           'c',
           'g',
           'y',
           #ROOT.gROOT.ProcessLine("kPink+1"),
           'darkorange',
           'r']

#! Marker types as per sequence
MRKS_PLT_SEQ={('ST'):ROOT.gROOT.ProcessLine("kFullDotLarge"),
			  ('SR'):ROOT.gROOT.ProcessLine("kFullDotLarge"),
			  ('SA'):ROOT.gROOT.ProcessLine("kFullDotLarge"),  
			  ('ER'):ROOT.gROOT.ProcessLine("kFullDotLarge"),#!kCircle
			  ('EC'):ROOT.gROOT.ProcessLine("kFullDotLarge"),
			  ('EF'):ROOT.gROOT.ProcessLine("kFullDotLarge")}

#! PSTYPs
NPSTYP=2
AL5D,ER5D=range(NPSTYP)
PSTYPS=['AL5D','ER5D']

#! THETA binning
NTHETABINS=10
THETAMIN,THETAMAX=0,180
THETABINW=18
THETABINLE=np.arange(THETAMIN,THETAMAX,THETABINW)

def get_q2wbinlist(f,q2min=0.00,q2max=6.00,wmin=0.000,wmax=3.000,dbg=False,dbg_bins=2,
	               dbg_binl=['2.00-2.40_1.425-1.450','2.00-2.40_1.450-1.475','2.00-2.40_1.475-1.500','2.00-2.40_1.500-1.525']):
		"""
		+ Taken from disp_obs.py and modified accordingly
		+ Note in dbg mode, functions works as expected when 'dbg_bins'=number of bins in 'dbg_binl'
		    + 'dbg_binl' contains bins which need immediate analysis
		"""
		q2wbinl=[]
		
		# print "DispObs::get_q2wbinlist() Going to Q2-W bins from file=",f.GetName()
		# print "DispObs::get_q2wbinlist() q2min,q2max,wmin,wmax=",q2min,q2max,wmin,wmax
		# if dbg==True:
		# 	print "DispObs::get_q2wbinlist() dbg=True"

		#brk=False #! technical tool to break out of two nested for loops. Set in second (nested) for loop
		#if brk==True: break
		for path,dirs,files in f.walk():
			#print "path,dirs,files=",path,dirs,files
			if path=="":continue #! Avoid root path
			path_arr=path.split("/")
			if len(path_arr)==1:
				if dbg==True:
					if path in dbg_binl:
						q2wbinl.append(path)
				else:
					q2wbinl.append(path)
			if dbg==True and len(q2wbinl)==dbg_bins:
				#brk=True
				break #! Uncomment/comment -> Get limited q2w-bins/Get all q2w-bins



		# #! Remove q2wbins that are not within [q2min,q2max],[wmin,wmax] 
		# q2wbins_remove=[]
		# fdor q2wbin in q2wbinl:
		# 	q2bin_le=q2wbin.split("_")[0].split("-")[0]
		# 	q2bin_ue=q2wbin.split("_")[0].split("-")[1]
		# 	wbin_le =q2wbin.split("_")[1].split("-")[0]
		# 	wbin_ue =q2wbin.split("_")[1].split("-")[1]
		# 	if float(q2bin_ue)<=q2min or float(q2bin_le)>=q2max or float(wbin_ue)<=wmin or float(wbin_le)>=wmax:
		# 		q2wbins_remove.append(q2wbin)
		# for q2wbin in q2wbins_remove:
		# 	q2wbinl.remove(q2wbin)

		return q2wbinl

def get_q2bin(q2wbin):
	return q2wbin.split('_')[0]

def norm_1D_theta(hTheta):
	#! 1. Create normalization factor histogram
	hDCosTheta=hTheta.Clone("hDCosTheta")
	hDCosTheta.SetTitle("hDCosTheta")
	hDCosTheta.Reset()
	nbins=hTheta.GetNbinsX()
	for ibin in range(nbins):
		theta_a=hTheta.GetBinLowEdge(ibin+1)
		theta_b=hTheta.GetBinLowEdge(ibin+2)# + hTheta.GetBinWidth(ibin+1)
		DCosTheta=math.fabs(math.cos(math.radians(theta_b))-math.cos(math.radians(theta_a)))
		hDCosTheta.SetBinContent(ibin+1,DCosTheta)
		hDCosTheta.SetBinError(ibin+1,0.)
	#! Now divide hTheta by hDCosTheta
	#! Do Sumw2() so that errors are correctly propagated
	hTheta.Sumw2();
	hTheta.Divide(hDCosTheta)
	return hDCosTheta

def plot_simstats(h5l,q2wbin,q2r):	
	'''
	Following plots are made:
	1. 5D-PS (For All-5D and ER-Common-5D)
		1.i. Total statistics for {ER,SR,ST} vs sim#              
		1.ii 5D-PS distribution for {ER,ST,SR,SA} for each sim#   
	
	2. 5D-PS-vst-var (For All-5D and ER-Common-5D) (current vst,var=1,THETA i.e. theta_pim)
		2.i.  Total statistics for {ER,SR,ST} vs vst-var, for each sim# 
		2.ii. 5D-PS-vst-var distributions for each 1D-bin for each sim# 
		2.iii 5D-PS-vst-var-dist-avg for {ST,SR,SA} vs vst-var, for each sim# 
	'''
	print "processing q2wbin=",q2wbin

	var=VAR_NAMES_PLAIN[1,'THETA']

	#SEQ=['ER','ST','SR','SA']

	#! outdir_q2w
	outdir_q2w="%s/simstats/%s"%(OUTDIR,q2wbin)
	if not os.path.exists(outdir_q2w):
		os.makedirs(outdir_q2w)

	#! Get total number of sims
	nsims=len(h5l)

	#! General aesthetics
	plot_simstats_athtcs()

	#! Now start processing

	#! Total statistics distributions: 1.i and 2.i
	#! 1.i. Total statistics for {ER,SR,ST} vs sim#:
	#!    stats[isim][seq,pstyp]
	#! + get data
	stats=[OrderedDict() for isim in range(nsims)]
	for isim in range(nsims):
		for seq in ['ER','SR','ST']:
			if seq=='ER': #! Note 'ER5D' not applicable for ER
				stats[isim][seq,'AL5D']=thntool.GetIntegral(h5l[isim][seq])
			else:
				stats[isim][seq,'AL5D']=thntool.GetIntegral(h5l[isim][seq])
				stats[isim][seq,'ER5D']=thntool.GetIntegralCommonBins(h5l[isim][seq],h5l[isim]['ER'])
	#! plot data {'ER','SR','ST'}*{pstyp} vs sim on one plot
	fig_size_x=20
	fig_size_y=15
	fig,axs = plt.subplots(figsize=(fig_size_x,fig_size_y),nrows=3, ncols=2,sharex=True)
	fig.suptitle('Total statistics',fontsize='xx-large')
	#! Set up x axis = sim#
	x=range(nsims)
	xticks=[SIMFRCTN[q2r][isim] for isim in range(nsims)]
	for seq in ['ST','SR','ER']:
		print "Processing",seq
		#! Set up axs for 'AL5D' ('ER5D'=ir+1,ic)
		if   seq=='ST': ir,ic=0,0
		elif seq=='SR':	ir,ic=0,1
		elif seq=='ER': ir,ic=2,0
		for pstyp in PSTYPS:
			print "Processing %s"%pstyp
			if seq=='ER' and pstyp=='ER5D': continue
			if   pstyp=='AL5D': ax=axs[ir][ic]
			elif pstyp=='ER5D': ax=axs[ir+1][ic]
			y=[stats[isim][seq,pstyp] for isim in range(nsims)]
			print x
			print y
			ax.set_title("%s-%s"%(seq,pstyp))
			ax.scatter(x,y,s=50)
			#! aesthetics
			#! yaxis
			ax.set_ylabel("N")
			ymin,ymax=0,( max(y)+(10/100)*max(y) )
			ax.set_ylim(ymin,ymax)
			#! xaxis
			ax.set_xlabel("Simulation Fraction")
			ax.set_xticks(x)
			ax.set_xticklabels(xticks)
			#! Fix x axis
			# shift half a step to the left
			xmin=(3*x[0]-x[1])/2.
			# shift half a step to the right
			xmax=(3*x[-1]-x[-2])/2.
			ax.set_xlim(xmin,xmax)
	#! Save
	outdir="%s/5D-PS-tot-stat/"%(outdir_q2w)
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	fig.savefig("/%s/totstats.png"%(outdir))
			

	#! 5D-PS distributions: 1.ii and 2.ii
	#! 1.ii 5D-PS distribution for {ER,ST,SR,SA} for each sim# 
	for seq in ['ST','SR','SA','ER']:
		print "5D-PS distributions: processing seq=",seq

		#! BinContentDist and BinRelErrDist: h1tot[isim]['AL5D','ER5D'] and h1errtot[isim]['AL5D','ER5D'] 
		h1tot   =[OrderedDict() for i in range(nsims)]
		h1errtot=[OrderedDict() for i in range(nsims)]
		#! setup binning information bng(seq)=(nbins,xmin,xmax)
		bng=OrderedDict()
		bng['ER']=[5,    0.5, 5.5]
		bng['SA']=[1000, 0,   0.5]
		bng['SR']=[100,  0.5,   100.5] #500
		bng['ST']=[0,    0,   0] #! default binning
		#! Now get histograms
		#! BinContentDist
		for isim in range(nsims):
			nbins,xmin,xmax=bng[seq][0],bng[seq][1],bng[seq][2]
			if seq=='ER': #! Note 'ER5D' not applicable for ER
				h1tot[isim]['AL5D']=thntool.GetBinContentDist(h5l[isim][seq],nbins,xmin,xmax)
				h1errtot[isim]['AL5D']=thntool.GetBinRelErrorDist(h5l[isim][seq],100,0,1.5)
			else:
				h1tot[isim]['AL5D']=thntool.GetBinContentDist(h5l[isim][seq],nbins,xmin,xmax)
				h1tot[isim]['ER5D']=thntool.GetBinContentDistCommonBins(h5l[isim][seq],h5l[isim]['ER'],nbins,xmin,xmax)
				h1errtot[isim]['AL5D']=thntool.GetBinRelErrorDist(h5l[isim][seq],100,0,1.5)
				h1errtot[isim]['ER5D']=thntool.GetBinRelErrorDistCommonBins(h5l[isim][seq],h5l[isim]['ER'],100,0,1.5)
		#! Some hist aesthetics
		for isim in range(nsims):
			for k in h1tot[isim].keys():
				#h1tot[isim][k].SetTitle("h5_BinContentDist_sim%d_%s_%s"%(isim+1,seq,q2wbin))
				h1tot[isim][k].SetTitle("h5_BinContentDist_%s_%s_%s"%(k,seq,q2wbin))
				clr=CLRS[isim]
				mrk=MRKS_PLT_SEQ[seq]
				# h1tot[isim][k].SetMarkerStyle(mrk)
				# h1tot[isim][k].SetMarkerColor(clr)
				h1tot[isim][k].SetLineColor(clr)
				# h1tot[isim][k].SetFillColor(clr)
				# h1tot[isim][k].SetFillStyle(3001)
			for k in h1errtot[isim].keys():
				#h1errtot[isim][k].SetTitle("h5_BinRelErrDist_sim%d_%s_%s"%(isim+1,seq,q2wbin))
				h1errtot[isim][k].SetTitle("h5_BinRelErrDist_%s_%s_%s"%(k,seq,q2wbin))
				clr=CLRS[isim]
				mrk=MRKS_PLT_SEQ[seq]
				# h1errtot[isim][k].SetMarkerStyle(mrk)
				# h1errtot[isim][k].SetMarkerColor(clr)
				h1errtot[isim][k].SetLineColor(clr)
				# h1errtot[isim][k].SetFillColor(clr)
				# h1errtot[isim][k].SetFillStyle(3001)
				
		#! Draw and save
		#! create canvas, one each for 'AL5D' and 'ER5D'
		for pstype in ['AL5D','ER5D']:
			if seq=='ER' and pstype=='ER5D': continue
			#! create outdir
			outdir="%s/%s/%s/%s"%(outdir_q2w,'5D-PS-dist',pstype,seq)
			if not os.path.exists(outdir):
				os.makedirs(outdir)
			c=ROOT.TCanvas()
			if SHOW_REL_ERR_DIST:
				c.SetCanvasSize(700,700)
				c.Divide(1,2)
				c.cd(1)
				#! before drawing set min and max
				h1l=[h1tot[isim][pstype] for isim in range(nsims)]
				set_h1l_maximum_minimum(h1l)
				#! finally draw
				for isim in range(nsims):
					draw_opt="hist"#"e"
					if isim>0:draw_opt="hist sames" #"e sames"
					h1tot[isim][pstype].Draw(draw_opt)
				c.cd(2)
				#! before drawing set min and max
				h1l=[h1errtot[isim][pstype] for isim in range(nsims)]
				set_h1l_maximum_minimum(h1l)
				#! finally draw
				for isim in range(nsims):
					draw_opt="hist"#"e"
					if isim>0:draw_opt="hist sames" #"e sames"
					h1errtot[isim][pstype].Draw(draw_opt)
			else:
				#! before drawing set min and max
				h1l=[h1tot[isim][pstype] for isim in range(nsims)]
				set_h1l_maximum_minimum(h1l)
				#! finally draw
				for isim in range(nsims):
					draw_opt="hist"#"e"
					if isim>0:draw_opt="hist sames" #"e sames"
					h1tot[isim][pstype].Draw(draw_opt)
			c.SaveAs("%s/h1.png"%(outdir))

	#! 2.ii. 5D-PS-vst-var-dist for each 1D-bin for each sim# 
	#! + Also, in this loop store data for 2.iii 5D-PS-vst-var-dist-avg for {ST,SR,SA} vs vst-var, for each sim#:
	#!    + avg[isim][seq,pstyp,ibin]
	avg=[OrderedDict() for isim in range(nsims)]
	for seq in ['ST','SR','SA']:
		print "5D-PS-vst-var distributions: processing seq=",seq
		#! Get VST1-THETA binning information from isim=0 (can use any sim#)
		nbins=h5l[0][seq].GetAxis(H5_DIM['THETA']).GetNbins()
		#print nbins
		binw=h5l[0][seq].GetAxis(H5_DIM['THETA']).GetBinWidth(1)
		for ibin in range(nbins):
			#! BinContentDist and BinRelErrDist: h1tot[isim]['AL5D','ER5D'] and h1errtot[isim]['AL5D','ER5D'] 
			h1   =[OrderedDict() for i in range(nsims)]
			h1err=[OrderedDict() for i in range(nsims)]
			#! setup binning information bng(seq)=(nbins,xmin,xmax)
			#! Use same bng as for 1.ii 5D-PS distribution for {ER,ST,SR,SA} for each sim# 
			# bng=OrderedDict()
			# bng['ER']=[5,    0.5, 5.5]
			# bng['SA']=[1000, 0,   0.5]
			# bng['SR']=[100,  0,   100] #500
			# bng['ST']=[0,    0,   0] #! default binning
			for isim in range(nsims):
				#! Set range and make projections as per bin
				binle=h5l[isim][seq].GetAxis(H5_DIM['THETA']).GetBinLowEdge(ibin+1)
				binue=binle+binw
				print "processing seq,sim#,bin#: %s,%d%d,[%.3f,%.3f)"%(seq,isim+1,ibin+1,binle,binue)
				#! Set range for seq and 'ER'
				h5l[isim]['ER'].GetAxis(H5_DIM['THETA']).SetRange(ibin+1,ibin+1)
				h5l[isim][seq].GetAxis(H5_DIM['THETA']).SetRange(ibin+1,ibin+1)
				#! h5->h4 for selected range
				h4ER=h5l[isim]['ER'].Projection(4,H4_PROJDIM,"E")
				h4  =h5l[isim][seq].Projection(4,H4_PROJDIM,"E")
				#! get h1 and h1err
				nbins,xmin,xmax=bng[seq][0],bng[seq][1],bng[seq][2]
				h1[isim]['AL5D']=thntool.GetBinContentDist(h4,nbins,xmin,xmax)
				h1[isim]['ER5D']=thntool.GetBinContentDistCommonBins(h4,h4ER,nbins,xmin,xmax)
				h1err[isim]['AL5D']=thntool.GetBinRelErrorDist(h4,100,0,1.5)
				h1err[isim]['ER5D']=thntool.GetBinRelErrorDistCommonBins(h4,h4ER,100,0,1.5)
				#! get data for avg[isim][seq,pstyp,ibin]
				avg[isim][seq,'AL5D',ibin]=np.zeros(2,'f')
				avg[isim][seq,'ER5D',ibin]=np.zeros(2,'f')
				if seq=='ER': #! Note 'ER5D' not applicable for ER
					thntool.GetBinContentDistStats(h4,avg[isim][seq,'AL5D',ibin])
				else:
					thntool.GetBinContentDistStats(h4,avg[isim][seq,'AL5D',ibin])
					thntool.GetBinContentDistStatsCommonBins(h4,h4ER,avg[isim][seq,'ER5D',ibin])
				#! Reset Range
				h5l[isim]['ER'].GetAxis(H5_DIM['THETA']).SetRange()
				h5l[isim][seq].GetAxis(H5_DIM['THETA']).SetRange()
							
			#! Some hist aesthetics
			for isim in range(nsims):
				for k in h1[isim].keys():
					#h1[isim][k].SetTitle("h5_BinContentDist_sim%d_%s_bin%d_%s"%(isim+1,seq,ibin+1,q2wbin))
					h1[isim][k].SetTitle("h5_BinContentDist_%s_%s_bin%d_%s"%(k,seq,ibin+1,q2wbin))
					clr=CLRS[isim]
					mrk=MRKS_PLT_SEQ[seq]
					# h1[isim][k].SetMarkerStyle(mrk)
					# h1[isim][k].SetMarkerColor(clr)
					h1[isim][k].SetLineColor(clr)
					# h1[isim][k].SetFillColor(clr)
					# h1[isim][k].SetFillStyle(3001)
				for k in h1err[isim].keys():
					h1err[isim][k].SetTitle("h5_BinRelErrDist_%s_%s_bin%d_%s"%(k,seq,ibin+1,q2wbin))
					clr=CLRS[isim]
					mrk=MRKS_PLT_SEQ[seq]
					# h1err[isim][k].SetMarkerStyle(mrk)
					# h1err[isim][k].SetMarkerColor(clr)
					h1err[isim][k].SetLineColor(clr)
					# h1err[isim][k].SetFillColor(clr)
					# h1err[isim][k].SetFillStyle(3001)

			#! Draw and save
			#! create canvas, one each for 'AL5D' and 'ER5D'
			for pstype in ['AL5D','ER5D']:
				#! create outdir
				outdir="%s/%s/%s/%s"%(outdir_q2w,'5D-PS-vst-var-dist',pstype,seq)
				if not os.path.exists(outdir):
					os.makedirs(outdir)
				c=ROOT.TCanvas()
				if SHOW_REL_ERR_DIST:
					c.SetCanvasSize(700,700)
					c.Divide(1,2)
					c.cd(1)
					#! before drawing set min and max
					h1l=[h1[isim][pstype] for isim in range(nsims)]
					set_h1l_maximum_minimum(h1l)
					#! finally draw
					for isim in range(nsims):
						draw_opt="hist"#"e"
						if isim>0:draw_opt="hist sames" #"e sames"
						h1[isim][pstype].Draw(draw_opt)
					c.cd(2)
					#! before drawing set min and max
					h1l=[h1err[isim][pstype] for isim in range(nsims)]
					set_h1l_maximum_minimum(h1l)
					#! finally draw
					for isim in range(nsims):
						draw_opt="hist"#"e"
						if isim>0:draw_opt="hist sames" #"e sames"
						h1err[isim][pstype].Draw(draw_opt)
				else:
					#! before drawing set min and max
					h1l=[h1[isim][pstype] for isim in range(nsims)]
					set_h1l_maximum_minimum(h1l)
					#! finally draw
					for isim in range(nsims):
						draw_opt="hist"#"e"
						if isim>0:draw_opt="hist sames" #"e sames"
						h1[isim][pstype].Draw(draw_opt)
				#c.SaveAs("%s/h1_%s.png"%(outdir,pstype))
				c.SaveAs("%s/h1_%02d_%04.1f-%04.1f.png"%(outdir,ibin+1,binle,binue))			

	#! Plot 2.iii data obtain earlier in avg[isim][seq,pstyp,ibin]
	fig_size_x=20
	fig_size_y=15
	fig,axs = plt.subplots(figsize=(fig_size_x,fig_size_y),nrows=4, ncols=2,sharex=True)
	fig.suptitle('Average',fontsize='xx-large')
	#! Set up x axis = theta bins
	x=range(NTHETABINS)
	xticks=[THETABINLE[ibin] for ibin in range(NTHETABINS)]
	for seq in ['ST','SR','SA']: #! ignore 'ST','SR','ER' for now
		print "Processing",seq
		#! Set up axs for 'AL5D' ('ER5D'=ir+1,ic)
		if   seq=='ST': ir,ic=0,0
		elif seq=='SR':	ir,ic=0,1
		elif seq=='SA': ir,ic=2,0
		for pstyp in PSTYPS:
			print "Processing %s"%pstyp
			if seq=='ER' and pstyp=='ER5D': continue
			if   pstyp=='AL5D': ax=axs[ir][ic]
			elif pstyp=='ER5D': ax=axs[ir+1][ic]
			#! get avg,avg_err as a function of sim
			y=[0 for isim in range(nsims)]
			yerr=[0 for isim in range(nsims)]
			yerr0=[0 for isim in range(nsims)]
			for isim in range(nsims):
				y[isim]=   [avg[isim][seq,pstyp,ibin][0] for ibin in range(NTHETABINS)]
				yerr[isim]=[avg[isim][seq,pstyp,ibin][1] for ibin in range(NTHETABINS)]
				yerr0[isim]=[0 for ibin in range(NTHETABINS)]
			if DBG==True:
				print x
				for isim in range(nsims):
					print "y[%d]:"%isim
					print y[isim]
					print "yerr[%d]:"%isim
					print yerr[isim]
			ax.set_title("%s-%s"%(seq,pstyp))
			for isim in range(nsims):
				# #! with error bars
				# ax.errorbar(x,y[isim],yerr=yerr[isim],fmt='o',label='sim-frc-%.1f'%SIMFRCTN[q2r][isim],markersize=8,c=CLRS_MPLT[isim])
				#! 0error bars
				ax.errorbar(x,y[isim],yerr=yerr0[isim],fmt='o',label='sim-frc-%.1f'%SIMFRCTN[q2r][isim],markersize=8,c=CLRS_MPLT[isim])
			#! aesthetics
			#! yaxis
			ax.set_ylabel("<%s>"%seq)
			lists=[y[isim] for isim in range(nsims)]
			ytot=[item for sublist in lists for item in sublist]
			ymin,ymax=0,( max(ytot)+(10/100)*max(ytot) )
			ax.set_ylim(ymin,ymax)
			#! xaxis
			ax.set_xlabel(r"$\theta_{\pi^{-}}$",size='xx-large')
			ax.set_xticks(x)
			ax.set_xticklabels(xticks)
			#! Fix x axis
			# shift half a step to the left
			xmin=(3*x[0]-x[1])/2.
			# shift half a step to the right
			xmax=(3*x[-1]-x[-2])/2.
			ax.set_xlim(xmin,xmax)
	#! Save
	outdir="%s/5D-PS-vst-var-dist-avg/"%(outdir_q2w)
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	fig.savefig("/%s/avg.png"%(outdir))

def plot_simstats_athtcs():
	#! General aesthetics
	#! + Reset any previous
	ROOT.gStyle.Reset()
	ROOT.gStyle.SetOptStat("nmMrReiuo")

def plot_1D_athtcs():
	#ROOT.gStyle.Reset()
	#! Stats Box
	ROOT.gStyle.SetOptStat(0)

	# ROOT.gStyle.SetLabelSize(0.5,"t")
	# ROOT.gStyle.SetTitleSize(0.5,"t")
	#ROOT.gStyle.SetPaperSize(20,26);
	#ROOT.gStyle.SetPadTopMargin(0.15)#(0.05);
	#ROOT.gStyle.SetPadRightMargin(0.15)#(0.05);
	ROOT.gStyle.SetPadBottomMargin(0.20)#(0.16);
	#ROOT.gStyle.SetPadLeftMargin(0.20)#(0.12);

	ROOT.gStyle.SetTitleW(10)# //title width 
	ROOT.gStyle.SetTitleFontSize(20)# 
	ROOT.gStyle.SetTitleH(0.15)# //title height 
		
	#! + The following options do not seem to work from here
	#! + I have to set them in label_hist_obs1D()
	#ROOT.gStyle.SetTitleFont(42,"xyz")
	#ROOT.gStyle.SetTitleSize(.35,"xyz")
	#ROOT.gStyle.SetTitleOffset(0.5,"xyz");

	#!get rid of X error bars and y error bar caps
	ROOT.gStyle.SetErrorX(0.001)

def set_h1l_maximum_minimum(h1l):
	'''
	+ From a list of h1s, obtain and set the minimum and maximum limits of the y-axis
	  on which all h1s can be completely displayed.
	+ If minimum>0, then simply set minimum=0
	+ Note that this is done in the most general way by looping over each bin of every
	  h1 in the list and using the bin-content and bin-error to determine the y-axis
	  limits
	'''
	#! obtain minimum and maximum
	yvals=[]
	for h1 in h1l:
		nbins=h1.GetNbinsX()
		for ibin in range(nbins):
			binc=h1.GetBinContent(ibin+1)
			bine=h1.GetBinError(ibin+1)
			yvals.append(binc+bine)
			yvals.append(binc-bine)
	maximum=max(yvals)
	minimum=min(yvals)
	if minimum>0: minimum=0
	#! Now set minimum and maximum with 10% padding (as ROOT does too)
	#! + Note special treatment if minimum is set to 0
	pdng=10/100
	if minimum!=0:
		minimum=minimum-(pdng*minimum)
	maximum=maximum+(pdng*maximum)
	for h1 in h1l:
		h1.SetMinimum(minimum)
		h1.SetMaximum(maximum)
	return

def hist_1D_athtcs(hobs):
	'''
	'''
	for isim,h1 in enumerate(hobs):
		for k in h1.keys():
			seq,vst,var=k[0],k[1],k[2]	
			mrk=MRKS_PLT_SEQ[seq]
			clr=CLRS[isim]
			h1[k].SetMarkerStyle(mrk)
			h1[k].SetMarkerColor(clr)
			h1[k].SetLineColor(clr)

def hist_itg_athtcs(hobs):
	'''
	'''
	for isim,h1 in enumerate(hobs):
		for k in h1.keys():
			seq=k	
			mrk=MRKS_PLT_SEQ[seq]
			clr=CLRS[isim]
			h1[k].SetMarkerStyle(mrk)
			h1[k].SetMarkerColor(clr)
			h1[k].SetLineColor(clr)
			
def plot_1D(hobs,q2wbin):
	"""
	"""
	print "In plot_1D()"

	#! Set up plotting aesthetics
	plot_1D_athtcs()
	hist_1D_athtcs(hobs)

	# #! For saving output to .root file, create and cd into q2wbindir
	# if self.FOUT_1D!=None:
	# 	q2wbin_dir=self.FOUT_1D.GetDirectory(q2wbin)
	# 	if q2wbin_dir==None: 
	# 		q2wbin_dir=self.FOUT_1D.mkdir(q2wbin)
	# 	q2wbin_dir.cd()
				
	# #! TCanvas's pad_map[pad,vst,var] defined as per Gleb's display
	# pad_map=[(1,1,"M1"),   (2,3,'M2'),   (3,2,'M2'),
	# 		 (4,1,"THETA"),(5,3,'THETA'),(6,2,'THETA'),
	# 		 (7,1,"ALPHA"),(8,3,'ALPHA'),(9,2,'ALPHA')]
		
	#! A canvas each for each sequence
	for seq in ['EC','EF']:				  
		c=ROOT.TCanvas("c","c",1000,1000)
		pad_t=ROOT.TPad("pad_l","Legend pad",0.25,0.935,0.75,1.00)
		pad_p=ROOT.TPad("pad_p","Plots pad",0.01,0.97,0.99,0.01)
		#if (self.VIEW=="fullana"):pad_p.SetFillColor(ROOT.gROOT.ProcessLine("kGray+2"))
		pad_p.Draw()
		pad_t.Draw()
		pad_t.cd()
		pt=ROOT.TPaveText(.05,.1,.95,.8)
		pt.AddText("%s Q2_W bin=%s"%(seq,q2wbin))
		pt.SetTextSize(0.40)
		pt.Draw()
		pad_p.Divide(3,3)
		for item in PAD_MAP_1D:
			pad,vst,var=item[0],item[1],item[2]
			#print "pad,vst,var=",pad,vst,var
			gpad=pad_p.cd(pad)
			#! First set maximum of y axis as per all hists that go in this canvas
			nsims=len(hobs)
			histl=[]
			for isim,h1 in enumerate(hobs):
				histl.append(h1[seq,vst,var])
			set_h1l_maximum_minimum(histl)
			#! Now draw
			iplt=0
			for isim,h1 in enumerate(hobs):
				# print 'sim#=',isim+1
				# for k in h1.keys():
				# 	print k,h1[k].GetName()
				#if isim==0: continue #! "skip"
				draw_opt="" if iplt==0 else "same"
				h1[seq,vst,var].Draw(draw_opt)
				iplt+=1
		#! Save canvas 
		outdir=os.path.join(OUTDIR,"Obs_1D",q2wbin)
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		c.SaveAs("%s/c_1D_%s.png"%(outdir,seq))
		
	print "*** plot_1D() Done ***\n"
	return

def plot_itg(hobs,q2wbin):
	"""
	"""
	print "In plot_itg()"

	#! Determine nsims
	nsims=len(hobs)
		
	#! General aesthetics
	ROOT.gStyle.SetOptStat(0)

	#! Histogram aesthetics
	hist_itg_athtcs(hobs)

	#! A canvas each for each sequence
	for seq in ['EC','EF']:				  
		c=ROOT.TCanvas("c","c",1200,800)
		pad_p=ROOT.TPad("pad_p","Plots pad",0.00,0.00,1.00,0.90)
		pad_t=ROOT.TPad("pad_p","Plots pad",0.00,0.90,1.00,1.00)
		pad_p.Draw()
		pad_t.Draw()
		pad_t.cd()
		pt=ROOT.TPaveText(.05,.1,.95,.8)
		q2bin=get_q2bin(q2wbin)
		pt.AddText("Integrated cross-section for %s,Q2=%s (vst,var=1,THETA)"%(seq,q2bin))
		pt.SetTextSize(0.32)
  		pt.Draw()
		pad_p.cd()
		#! Set maximum of y axis as per all hists that go in this canvas
		nsims=len(hobs)
		histl=[]
		for isim,h1 in enumerate(hobs):
			histl.append(h1[seq])
		set_h1l_maximum_minimum(histl)
		#! Begin plotting
		for isim in range(nsims):
			draw_opt="" if isim==0 else "same"
			#! Set up some aesthetics particular to hW
			#! + Give already rotated bin labels by 90 deg. more space
			#! + Move xaxis title lower to accomodate 90 deg. rotated titles
			pad_p.SetBottomMargin(0.20)
			hobs[isim][seq].GetXaxis().SetTitleOffset(2.80)
			hobs[isim][seq].GetXaxis().SetTitleOffset(2.80)
			#! Draw
			hobs[isim][seq].SetTitle("")
			hobs[isim][seq].Draw(draw_opt)
		#! Save canvas 
		outdir=os.path.join(OUTDIR,"Obs_itg",q2bin)
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		c.SaveAs("%s/c_itg_%s.png"%(outdir,seq))
		
	print "*** plot_itg() Done ***\n"
	return
#! main

#! debug
#print FYLD[LQ2][0]
#print FYLD[HQ2][0]
#print FO1D[LQ2][1]
#sys.exit()

#! First get Q2WBINL from ISIM=0 (any ISIM should get the same list)
Q2WBINL=[0 for i in range(NQ2RANGES)]
if DBG==True:
	#Q2WBINL[LQ2]=get_q2wbinlist(FYLD[LQ2][0],dbg=True,dbg_bins=2,dbg_binl=['2.00-2.40_1.500-1.525','2.40-3.00_1.500-1.525'])
	Q2WBINL[LQ2]=get_q2wbinlist(FYLD[LQ2][0],dbg=True,dbg_bins=1,dbg_binl=['2.00-2.40_1.500-1.525'])
	#Q2WBINL[HQ2]=get_q2wbinlist(FYLD[HQ2][0],dbg=True,dbg_bins=1,dbg_binl=['3.00-3.50_1.500-1.525'])
else:
	#Q2WBINL[LQ2]=get_q2wbinlist(FYLD[LQ2][0])
	#Q2WBINL[HQ2]=get_q2wbinlist(FYLD[HQ2][0])
	Q2WBINL[LQ2]=get_q2wbinlist(FYLD[LQ2][0],dbg=True,dbg_bins=8,dbg_binl=['2.00-2.40_1.425-1.450','2.00-2.40_1.450-1.475','2.00-2.40_1.475-1.500','2.00-2.40_1.500-1.525',
		                                                                   '2.40-3.00_1.425-1.450','2.40-3.00_1.450-1.475','2.40-3.00_1.475-1.500','2.40-3.00_1.500-1.525'])
	Q2WBINL[HQ2]=get_q2wbinlist(FYLD[HQ2][0],dbg=True,dbg_bins=8,dbg_binl=['3.50-4.20_1.425-1.450','3.50-4.20_1.450-1.475','3.50-4.20_1.475-1.500','3.50-4.20_1.500-1.525',
		                                                                   '4.20-5.00_1.425-1.450','4.20-5.00_1.450-1.475','4.20-5.00_1.475-1.500','4.20-5.00_1.500-1.525'])	                                                                	
print "Q2WBINL[LQ2]=",Q2WBINL[LQ2]
print "Q2WBINL[HQ2]=",Q2WBINL[HQ2]
#sys.exit()

# #! 1. Loop over Q2WBINL and for each q2wbin make Obs_1D=f(sim)
# for q2r in range(NQ2RANGES):
# 	if isinstance(Q2WBINL[q2r],list)==False: continue #! in debug mode Q2WBINL[q2r] may not be obtained
# 	for q2wbin in Q2WBINL[q2r]:
# 		#! Create structure to hold hobs[isim][seq,vst,var]
# 		hobs=[OrderedDict() for isim in range(NSIM[q2r])]
# 		for isim in ISIM[q2r]:
# 			#if q2r==LQ2 and isim==0: continue #! data corrurpt. remaking. skip for now
# 			#! Get hobs
# 			for item in PAD_MAP_1D:
# 				pad,vst,var=item[0],item[1],item[2]
# 				for seq in ['EC','EF']:
# 					hobs[isim][seq,vst,var]=FO1D[q2r][isim].Get("%s/h1_%s_%d_%s"%(q2wbin,seq,vst,var))
# 		if DBG==True:
# 			print "***hobs pretty print for q2wbin",q2wbin,"***"
# 			for isim,h1 in enumerate(hobs):
# 				print "** sim#=",isim+1,"**"
# 				for k in h1:
# 					print k,h1[k].GetName()
# 		#! plot
# 		plot_1D(hobs,q2wbin)

# #! 2. Plot Obs_itg=f(sim) 
# #! + use vst=1, var=THETA
# #! First get Q2WBINL (different format compared to Obs_1D) from ISIM=0 (any ISIM should get the same list)
# Q2WBINL_ITG=[0 for i in range(NQ2RANGES)]
# #! Note that no debug mode needed since there only a few bins
# Q2WBINL_ITG[LQ2]=get_q2wbinlist(FOIT[LQ2][0])
# Q2WBINL_ITG[HQ2]=get_q2wbinlist(FOIT[HQ2][0])	
# print "Q2WBINL_ITG[LQ2]=",Q2WBINL_ITG[LQ2]
# print "Q2WBINL_ITG[HQ2]=",Q2WBINL_ITG[HQ2]
# for q2r in range(NQ2RANGES):
# 	for q2wbin in Q2WBINL_ITG[q2r]:
# 		#! Create structure to hold hobs[isim][seq,vst,var]
# 		hobs=[OrderedDict() for isim in range(NSIM[q2r])]
# 		for isim in ISIM[q2r]:
# 			#! Get hobs
# 			for seq in ['EC','EF']:
# 				hobs[isim][seq]=FOIT[q2r][isim].Get("%s/hW_%s_1_THETA"%(q2wbin,seq))
# 		if DBG==True:
# 			print "***hobs pretty print for q2wbin",q2wbin,"***"
# 			for isim,h1 in enumerate(hobs):
# 				print "** sim#=",isim+1,"**"
# 				for k in h1:
# 					print k,h1[k].GetName()
# 		#! plot
# 		plot_itg(hobs,q2wbin)

#! 3. Loop over Q2WBINL and for each q2wbin make simstats=f(sim)
for q2r in range(NQ2RANGES):
	if isinstance(Q2WBINL[q2r],list)==False: continue #! in debug mode Q2WBINL[q2r] may not be obtained
	for q2wbin in Q2WBINL[q2r]:
		#! Create structure to hold h5[isim][seq] (note vst=VST1 for now)
		h5l=[OrderedDict() for isim in range(NSIM[q2r])]
		for isim in ISIM[q2r]:
			#! Get h5
			for seq in ['ST','SR','SA','ER']:
				h5l[isim][seq]=FYLD[q2r][isim].Get("%s/%s/VST%d/h5"%(q2wbin,seq,1))
		if DBG==True:
			print "***h5 pretty print for q2wbin",q2wbin,"***"
			for isim,h5 in enumerate(h5l):
				print "** sim#=",isim+1,"**"
				for k in h5:
					print k,h5[k].GetName()
		#! plot
		plot_simstats(h5l,q2wbin,q2r)