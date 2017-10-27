#!/usr/bin/python
from __future__ import division
import os,sys,datetime
from collections import OrderedDict
import array

import ROOT

from rootpy.io import root_open, DoesNotExist

import math

'''

'''

USAGE='study_obs_as_function_of_simstats dbg[=False] '

#! user inputs
DBG=False
if len(sys.argv)>1: #! dbg entered by user
	if    sys.argv[1]=="True":  DBG=True
	elif  sys.argv[1]=="False": DBG=False
	else: sys.exit('DBG=%s is not valid. usage: %s'%(sys.argv[1],USAGE))

print "DBG=",DBG
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
ROOT.gStyle.SetOptStat("nmMrReiuo")
#ROOT.gStyle.SetStatStyle(0) #! transparent stats box

#! OUTDIR
DATE=datetime.datetime.now().strftime('%m%d%y')
OUTDIRNAME='results'
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

#! Setup input root files: FYLD[q2r][sim], FOBS[q2r][sim]
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
#! Finally setup and fill FYLD[q2r][sim],FOBS[q2r][sim]
FYLD=[[0 for i in range(NSIM[LQ2])],[0 for i in range(NSIM[HQ2])]]
FOBS=[[0 for i in range(NSIM[LQ2])],[0 for i in range(NSIM[HQ2])]]
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
			FOBS[q2][isim]=root_open('%s/lowQ2_SSBands_off_off_%s_%s/cutsncors1/%s/Obs_1D_norm/obs_1D.root'%(os.environ['OBSDIR_E16'],SIMNAME[q2][isim],date,SIMNAME[q2][isim]),'r')
	elif q2==HQ2:
		for isim in ISIM[HQ2]:
			if isim==0: date='102617'
			else:       date='102717'
			FYLD[q2][isim]=root_open('%s/highQ2_SSBands_off_off_%s_%s/cutsncors1/%s/yield.root'%(os.environ['OBSDIR_E16'],SIMNAME[q2][isim],date,SIMNAME[q2][isim]),'r')
			FOBS[q2][isim]=root_open('%s/highQ2_SSBands_off_off_%s_%s/cutsncors1/%s/Obs_1D_norm/obs_1D.root'%(os.environ['OBSDIR_E16'],SIMNAME[q2][isim],date,SIMNAME[q2][isim]),'r')
if DBG==True:
	print FYLD
	print FOBS
#sys.exit()

#! Set up other analysis constants
#! TCanvas's pad_map[pad,vst,var] defined as per Gleb's display
PAD_MAP_1D=[(1,1,"M1"),   (2,3,'M2'),   (3,2,'M2'),
		    (4,1,"THETA"),(5,3,'THETA'),(6,2,'THETA'),
		    (7,1,"ALPHA"),(8,3,'ALPHA'),(9,2,'ALPHA')]

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
				if DBG==True:
					if path in dbg_binl:
						q2wbinl.append(path)
				else:
					q2wbinl.append(path)
			if DBG==True and len(q2wbinl)==dbg_bins:
				#brk=True
				break #! Uncomment/comment -> Get limited q2w-bins/Get all q2w-bins



		# #! Remove q2wbins that are not within [q2min,q2max],[wmin,wmax] 
		# q2wbins_remove=[]
		# for q2wbin in q2wbinl:
		# 	q2bin_le=q2wbin.split("_")[0].split("-")[0]
		# 	q2bin_ue=q2wbin.split("_")[0].split("-")[1]
		# 	wbin_le =q2wbin.split("_")[1].split("-")[0]
		# 	wbin_ue =q2wbin.split("_")[1].split("-")[1]
		# 	if float(q2bin_ue)<=q2min or float(q2bin_le)>=q2max or float(wbin_ue)<=wmin or float(wbin_le)>=wmax:
		# 		q2wbins_remove.append(q2wbin)
		# for q2wbin in q2wbins_remove:
		# 	q2wbinl.remove(q2wbin)

		return q2wbinl

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

def plot(q2wbin):	
	'''
	+ For {ER,ST,SR,SA} plot:
	    i.  5D-PS distributions of entries and their corresponding rel_err
	    ii. 1D projection (currently only for VST1,THETA i.e. theta_pim)
	+ Additionally for {ST,SR,SA}
	    iii.  5D-PS-vst-var distributions of entries and their corresponding rel_err
	    	  + 5D-PS-vst-var= 5D PS distributions but only for a particular subspace of 
	    	    the total 5D-PS, which corresponds to a particular VST-VAR's bin.
	          + Currently this is done only for VST1,THETA i.e. theta_pim
	'''
	print "processing q2wbin=",q2wbin

	var=VAR_NAMES_PLAIN[1,'THETA']

	SEQ=['ER','ST','SR','SA']

	#! Get all the h5(seq)
	#! + First establish if to use FIN[LQ2] or FIN[HQ2] for this q2wbin
	if   FIN[LQ2].FindObjectAny(q2wbin)!=None: fin=FIN[LQ2]
	elif FIN[HQ2].FindObjectAny(q2wbin)!=None: fin=FIN[HQ2]
	else:
		sys.exit("study_yields_acceptance.py:plot():Could not establish if to use FIN[LQ2] or FIN[HQ2] for q2wbin=%s"(q2wbin))
	print "fin=",fin.GetName()

	h5=OrderedDict()
	for seq in SEQ:
		h5[seq]=fin.Get('%s/%s/VST1/h5'%(q2wbin,seq))
	
	#! outdir_q2w
	outdir_q2w="%s/%s"%(OUTDIR,q2wbin)
	if not os.path.exists(outdir_q2w):
		os.makedirs(outdir_q2w)

	#! Now start processing
	for seq in SEQ:
		print "processing seq=",seq

		#! i. 5D-PS distributions
		#! BinContentDist
		outdir="%s/%s/%s"%(outdir_q2w,seq,'5D-PS')
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		if seq=='ER': #! Note ANA_ER_COMMONBINS not applicable for ER
			h1tot=thntool.GetBinContentDist(h5[seq],5,0.5,5.5)
		if seq=='SA':
			if ANA_ER_COMMONBINS: h1tot=thntool.GetBinContentDistCommonBins(h5[seq],h5['ER'],1000,0,0.5)
			else:       		  h1tot=thntool.GetBinContentDist(h5[seq],1000,0,0.5)
			
		elif seq=='SR':
			if ANA_ER_COMMONBINS: h1tot=thntool.GetBinContentDistCommonBins(h5[seq],h5['ER'],100,0,500)
			else:                 h1tot=thntool.GetBinContentDist(h5[seq],100,0,500)
			
		elif seq=='ST':
			if ANA_ER_COMMONBINS: h1tot=thntool.GetBinContentDistCommonBins(h5[seq],h5['ER'])
			else:                 h1tot=thntool.GetBinContentDist(h5[seq])
			
		#! BinRelErrDist
		if seq=='ER': #! Note ANA_ER_COMMONBINS not applicable for ER
			h1errtot=thntool.GetBinRelErrorDist(h5[seq],100,0,1.5)
		else:
			if ANA_ER_COMMONBINS: h1errtot=thntool.GetBinRelErrorDistCommonBins(h5[seq],h5['ER'],100,0,1.5)
			else:                 h1errtot=thntool.GetBinRelErrorDist(h5[seq],100,0,1.5)
		
		#! Some hist aesthetics
		h1tot.SetTitle("h5_BinContentDist_%s_%s"%(seq,q2wbin))
		h1errtot.SetTitle("h5_BinRelErrDist_%s_%s"%(seq,q2wbin))
		
		#! Draw and save
		c=ROOT.TCanvas()
		if SHOW_REL_ERR_DIST:
			c.SetCanvasSize(700,700)
			c.Divide(1,2)
			c.cd(1)
			h1tot.Draw("e")
			c.cd(2)
			h1errtot.Draw("e")
		else:
			h1tot.Draw("e")
		c.SaveAs("%s/h1tot.png"%outdir)

		#! ii. 1D projection (currently only for VST1,THETA i.e. theta_pim)
		outdir="%s/%s/%s"%(outdir_q2w,seq,'1D-proj-%s'%var)
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		h1d=h5[seq].Projection(H5_DIM['THETA'],"E")
		h1d.SetMinimum(0)
		h1d.SetTitle("1D_proj_%s_%s_%s"%(seq,q2wbin,var))
		c=ROOT.TCanvas("c","c",700,700)
		c.Divide(1,2)
		c.cd(1)
		h1d.Draw("e")
		#! Save DCosTheta normalized version too
		c.cd(2)
		h1d_theta_norm=h1d.Clone("h1d_theta_norm")
		h1d_theta_norm.SetTitle("1D_proj_theta_norm_%s_%s_%s"%(seq,q2wbin,var))
		hDCosTheta=norm_1D_theta(h1d_theta_norm)
		h1d_theta_norm.Draw("e")
		c.SaveAs("%s/h1d.png"%outdir)
		#! Save hDCosTheta too
		c1=ROOT.TCanvas("c1","c1")
		hDCosTheta.Draw()
		c1.SaveAs("%s/hDCosTheta.png"%outdir)

		#sys.exit()

		#! iii. 5D-PS-vst-var distributions
		#! + only for {'ST','SR','SA'}
		if seq=='ER': continue
		outdir_var="%s/%s/%s"%(outdir_q2w,seq,'5D-PS-%s'%var)
		if not os.path.exists(outdir_var):
			os.makedirs(outdir_var)
		#! Get VST1-THETA binning information
		nbins=h5[seq].GetAxis(H5_DIM['THETA']).GetNbins()
		#print nbins
		binw=h5[seq].GetAxis(H5_DIM['THETA']).GetBinWidth(1)
		for ibin in range(nbins):
			#! Set range and make projections as per bin
			binle=h5[seq].GetAxis(H5_DIM['THETA']).GetBinLowEdge(ibin+1)
			binue=binle+binw
			print "processing bin number: %d [%.3f,%.3f)"%(ibin+1,binle,binue)
			#! Set range
			h5['ER'].GetAxis(H5_DIM['THETA']).SetRange(ibin+1,ibin+1)
			h5[seq].GetAxis(H5_DIM['THETA']).SetRange(ibin+1,ibin+1)
			#! h5->h4 for selected range
			h4ER=h5['ER'].Projection(4,H4_PROJDIM,"E")
			h4=h5[seq].Projection(4,H4_PROJDIM,"E")
			#! get h1
			if seq=='SA':
				if ANA_ER_COMMONBINS: h1=thntool.GetBinContentDistCommonBins(h4,h4ER,1000,0,0.5)
				else:                 h1=thntool.GetBinContentDist(h4,1000,0,0.5)
			elif seq=='SR':
				if ANA_ER_COMMONBINS: h1=thntool.GetBinContentDistCommonBins(h4,h4ER,100,0,200)
				else:                 h1=thntool.GetBinContentDist(h4,100,0,200)
			else:
				if ANA_ER_COMMONBINS: h1=thntool.GetBinContentDistCommonBins(h4,h4ER)
				else:                 h1=thntool.GetBinContentDist(h4)
			
			#! get h1err
			if ANA_ER_COMMONBINS: h1err=thntool.GetBinRelErrorDistCommonBins(h4,h4ER,100,0,1)
			else:                 h1err=thntool.GetBinRelErrorDist(h4)
			
			#! Some hist aesthetics
			h1.SetTitle("h5_BinContentDist_%s_%s_%s_bin%d"%(seq,q2wbin,var,ibin+1))
			h1err.SetTitle("h5_BinRelErrDist_%s_%s_%s_bin%d"%(seq,q2wbin,var,ibin+1))
			
			#! Draw and save
			c=ROOT.TCanvas()
			if SHOW_REL_ERR_DIST:
				c.SetCanvasSize(700,700)
				c.Divide(1,2)
				c.cd(1)
				h1.Draw("e")
				c.cd(2)
				h1err.Draw("e")
			else:
				h1.Draw("e")
			c.SaveAs("%s/h1_%02d_%04.1f-%04.1f.png"%(outdir_var,ibin+1,binle,binue))

			#! Reset Range
			h5['ER'].GetAxis(H5_DIM['THETA']).SetRange()
			h5[seq].GetAxis(H5_DIM['THETA']).SetRange()

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
	#! Marker colors as per simulation stats: hot -> cold = low -> high
	CLRS=[ROOT.gROOT.ProcessLine("kBlue"),
	      ROOT.gROOT.ProcessLine("kCyan"),
	      ROOT.gROOT.ProcessLine("kGreen"),
	      ROOT.gROOT.ProcessLine("kYellow"),
	      #ROOT.gROOT.ProcessLine("kPink+1"),
	      ROOT.gROOT.ProcessLine("kOrange"),
	      ROOT.gROOT.ProcessLine("kRed")]

	#! Marker types as per sequence
	MRKS_PLT_SEQ={('EC'):ROOT.gROOT.ProcessLine("kFullDotLarge"),
				  ('EF'):ROOT.gROOT.ProcessLine("kFullDotLarge")}

	for isim,h1 in enumerate(hobs):
		for k in h1.keys():
			seq,vst,var=k[0],k[1],k[2]	
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
		pt.AddText("Q2_W bin=%s"%q2wbin)
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
		outdir=os.path.join(OUTDIR,q2wbin)
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		c.SaveAs("%s/c_1D_%s.png"%(outdir,seq))
		
	print "*** plot_1D() Done ***\n"
	return
#! main

#! debug
#print FYLD[LQ2][0]
#print FYLD[HQ2][0]
#print FOBS[LQ2][1]
#sys.exit()

#! First get Q2WBINL from ISIM=0 (any ISIM should get the same list)
Q2WBINL=[0 for i in range(NQ2RANGES)]
if DBG==True:
	Q2WBINL[LQ2]=get_q2wbinlist(FYLD[LQ2][0],dbg_bins=2,dbg_binl=['2.00-2.40_1.500-1.525','2.40-3.00_1.500-1.525'])
	Q2WBINL[HQ2]=get_q2wbinlist(FYLD[HQ2][0],dbg_bins=1,dbg_binl=['3.00-3.50_1.500-1.525'])
else:
	Q2WBINL[LQ2]=get_q2wbinlist(FYLD[LQ2][0])
	Q2WBINL[HQ2]=get_q2wbinlist(FYLD[LH2][0])	
print "Q2WBINL[LQ2]=",Q2WBINL[LQ2]
print "Q2WBINL[HQ2]=",Q2WBINL[HQ2]
#sys.exit()
#! Now loop over Q2WBINL and for each q2wbin make obs=f(sim)
for q2r in range(NQ2RANGES):
	for q2wbin in Q2WBINL[q2r]:
		#! Create structure to hold hobs[isim][seq,vst,var]
		hobs=[OrderedDict() for isim in range(NSIM[q2r])]
		for isim in ISIM[q2r]:
			#if q2r==LQ2 and isim==0: continue #! data corrurpt. remaking. skip for now
			#! Get hobs
			for item in PAD_MAP_1D:
				pad,vst,var=item[0],item[1],item[2]
				for seq in ['EC','EF']:
					hobs[isim][seq,vst,var]=FOBS[q2r][isim].Get("%s/h1_%s_%d_%s"%(q2wbin,seq,vst,var))
		if DBG==True:
			print "***hobs pretty print for q2wbin",q2wbin,"***"
			for isim,h1 in enumerate(hobs):
				print "** sim#=",isim+1,"**"
				for k in h1:
					print k,h1[k].GetName()
		#! plot
		plot_1D(hobs,q2wbin)


