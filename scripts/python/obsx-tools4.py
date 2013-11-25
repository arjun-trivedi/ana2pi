#added extraction of 1D
#made creating of DF use more of vectorized operations
from __future__ import division

import os
from array import array
import math
import subprocess

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

from rootpy.io import root_open, File
import rootpy.plotting.root2matplotlib as rplt
from rootpy.plotting import Hist

from ROOT import gSystem, gROOT, gStyle, THnSparseF, TCanvas, TString, TLine, TPaveText, TText, TH1D
gSystem.Load('myTHnTool_C')
from ROOT import myTHnTool
mythnt = myTHnTool()

"""Test Case for plotR2"""
out=open('/tmp/stdout.text','w')
err=open('/tmp/stderr.text','w')

NDTYPE=2
EXP,SIM=range(NDTYPE)
DTYPES_COLOR=['kGreen','kRed']

NVARSETS=3

NVARS=5
M1,M2,THETA,PHI,ALPHA=range(NVARS)
VARS_NAME = ['M1','M2','THETA','PHI','ALPHA']
VARS_TITLE = [ 
				["M_{p#pi^{+}}", "M_{#pi^{+}#pi^{-}}","#theta_{#pi^{-}}", "#phi_{#pi^{-}}", "#alpha_{[p^{'}#pi^{+}][p#pi^{-}]}"],
			 	["M_{p#pi^{+}}", "M_{#pi^{+}#pi^{-}}","#theta_{p}", "#phi_{p}", "#alpha_{[#pi^{+}#pi^{-}][p^{'}p]"],
			 	["M_{p#pi^{+}}", "M_{p#pi^{-}}", "#theta_{#pi^{+}}", "#phi_{#pi^{+}}", "#alpha_{[p^{'}#pi^{-}][p#pi^{+}]"] 
		   ]
VARS_UNIT = ['GeV','GeV','#circ','#circ','#circ'] 

# NSEQ_SIM=6
# ST,SR,SA,SC,SH,SF=range(NSEQ_SIM)

# NSEQ_EXP=3
# ER,EC,EF=range(NSEQ_EXP)

NSEQ=9
ST,SR,SA,SC,SH,SF,ER,EC,EF=range(NSEQ)
SEQ_NAME=['ST','SR','SA','SC','SH','SF','ER','EC','EF']

NPOLS=3
POS,NEG,UNP=range(NPOLS)
POLS_COLOR=['kRed','kBlack','kBlue']
POLS_NAME=['POS','NEG','UNP']

ANADIR = os.environ['E1F_2PI_ANADIR2']
TOP = "1:2:3:4"
Q2WBNG = '1-2.000-2.400__24-1.300-1.900'
SIMDIR = os.path.join(ANADIR,'simdir')

SEQ_POLS_H5FILE = np.zeros((NSEQ,NPOLS),object)
SEQ_POLS_H5FILE[ST:SF+1,POS:NEG+1]  ='%s/%s__%s__pol__sim.root'%(SIMDIR,TOP,Q2WBNG)
SEQ_POLS_H5FILE[ST:SF+1,UNP:]       ='%s/%s__%s__sim.root'%(SIMDIR,TOP,Q2WBNG)
SEQ_POLS_H5FILE[ER:,POS:NEG+1]      ='%s/%s__%s__pol__exp.root'%(ANADIR,TOP,Q2WBNG)
SEQ_POLS_H5FILE[ER:,UNP:]           ='%s/%s__%s__exp.root'%(ANADIR,TOP,Q2WBNG)
#print SEQ_POLS_H5FILE

SEQ_POLS_H5=[
	['hY5D/Varset1/hY5D_TH',      'hY5D/Varset1/hY5D_TH',      'hY5D/Varset1/hY5D_TH'],
	['hY5D/Varset1/hY5D_RECO',    'hY5D/Varset1/hY5D_RECO',    'hY5D/Varset1/hY5D_RECO'],
	['hY5D/Varset1/hY5D_ACC',     'hY5D/Varset1/hY5D_ACC',     'hY5D/Varset1/hY5D_ACC'],
	['hY5D/Varset1/hY5D_ACC_CORR','hY5D/Varset1/hY5D_ACC_CORR','hY5D/Varset1/hY5D_ACC_CORR'],
	['hY5D/Varset1/hY5D_HOLE',    'hY5D/Varset1/hY5D_HOLE',    'hY5D/Varset1/hY5D_HOLE'],
	['hY5D/Varset1/hY5D_FULL',    'hY5D/Varset1/hY5D_FULL',    'hY5D/Varset1/hY5D_FULL'],
	['hY5D_POS/Varset1/hY5D_RECO',    'hY5D_NEG/Varset1/hY5D_RECO',    'hY5D/Varset1/hY5D_RECO'],
	['hY5D_POS/Varset1/hY5D_ACC_CORR','hY5D_NEG/Varset1/hY5D_ACC_CORR','hY5D/Varset1/hY5D_ACC_CORR'],
	['hY5D_POS/Varset1/hY5D_FULL',    'hY5D_NEG/Varset1/hY5D_FULL',    'hY5D/Varset1/hY5D_FULL']
]

NPOBS=4
A,B,C,D=range(NPOBS) # [A,B,C,D]=[Rt+Rl,Rlt,Rtt,Rlt']
POBS_NAME=['A','B','C','D']
#print var

def plot1D():
	"""On a TCanvas, plot 1D xsec"""
	outdir = os.path.join(outdir_1D_root,q2wbin)
	if not os.path.isdir(outdir):
		os.makedirs(outdir)

	#cname='1D_q2wbin_%d'%d_q2w['q2wbinnum'].tolist()[0]
	#ctitle='1D:%s'%q2wbin
	cname=ctitle='1D'
	c1D = TCanvas(cname,ctitle)
	c1D.Divide(3,3)
	npads=9
	hname=[]
	hname=['h1_1M1',    'h1_1M2',   'h1_3M2',#TODO: h1_1M2->h1_2M2
	       'h1_1THETA', 'h1_2THETA','h1_3THETA',
	       'h1_1ALPHA', 'h1_2ALPHA','h1_3ALPHA']
	htitle=[]
	htitle=[VARS_TITLE[0][M1],   VARS_TITLE[0][M2],   VARS_TITLE[1][M2],
	        VARS_TITLE[0][THETA],VARS_TITLE[1][THETA],VARS_TITLE[2][THETA],
	        VARS_TITLE[0][ALPHA],VARS_TITLE[1][ALPHA],VARS_TITLE[2][ALPHA]]
	hxtitle=['GeV','GeV','GeV','#circ','#circ','#circ']
	
	activepads=[1,4,2] 


	EF_UNP_sel = (d_q2w['SEQ']==EF) & (d_q2w['POL']==UNP)
	SF_UNP_sel = (d_q2w['SEQ']==SF) & (d_q2w['POL']==UNP)
	d_q2w_EF_UNP = d_q2w[EF_UNP_sel]
	d_q2w_SF_UNP = d_q2w[SF_UNP_sel]
	for ipad in range(0,npads):
		if ipad+1 not in activepads: continue
		c1D.cd(ipad+1)
		h1_exp = d_q2w_EF_UNP.iloc[0][hname[ipad]]
		h1_sim = d_q2w_SF_UNP.iloc[0][hname[ipad]]

		h1_exp.SetTitle('%s:%s'%(htitle[ipad],q2wbin))
		h1_exp.SetXTitle('%s[%s]'%(htitle[ipad],hxtitle[ipad]))
		h1_sim.SetLineColor(1)
		
		h1_exp.DrawNormalized("ep",10000)
		h1_sim.DrawNormalized("ep same",10000)

	csavename=('%s/%s')%(outdir,c1D.GetName())
	c1D.SaveAs(('%s.png')%(csavename))
	print ('>>>convert %s.png %s.pdf')%(csavename,csavename)
	rc = subprocess.call(['convert', '%s.png'%csavename, '%s.pdf'%csavename])
	if rc!=0:print '.png to .pdf failed for %s'%csavename
	

def plotR2(seq_pol_sell):
	"""On a TCanvas, plot R2 for a list of particular SEQ and POL selections
	
	Keyword arguments:
	seq_pol_sell -- A List of DataFrame selections based on SEQ and POL

	"""
	l = TLine(0,0,180,0)# will be used for all hR2 objects

	for seq_pol_sel in seq_pol_sell:
		d_q2w_seq_pol = d_q2w[seq_pol_sel]
		#print dq2w_seq_pol
		seq=d_q2w_seq_pol.iloc[0]['SEQ']
		pol=d_q2w_seq_pol.iloc[0]['POL']
		outdir_seq_pol=os.path.join(outdir_q2w,SEQ_NAME[seq],POLS_NAME[pol])
		if not os.path.isdir(outdir_seq_pol):
			os.makedirs(outdir_seq_pol)

		for pob in range(0,NPOBS):
			if pob==A: continue
			for var in range(0,NVARS):
				if var==PHI or var==ALPHA: continue
				hR2_name = 'hR2_%s_1%s'%(POBS_NAME[pob],VARS_NAME[var])
				hR2 = d_q2w_seq_pol.iloc[0][hR2_name]
				cR2 = TCanvas(hR2_name,hR2_name)
				hR2.SetLineColor(gROOT.ProcessLine("%s"%POLS_COLOR[pol]))
				hR2.SetTitle("")# will be made "prettier" later
				hR2.Draw("ep")

				#make Title of hR2 "pretty"
				l.Draw("same")
				pt = TPaveText(0.3, 0.85, 0.7, 1.0, "NDC")
				q2wt = pt.AddText('[Q^{2}][W] = %s'%q2wbin)
				q2wt.SetTextColor(gROOT.ProcessLine("kBlue"))
				vart = pt.AddText(("%s,%s: %s^{%s} vs. %s")%(SEQ_NAME[seq],POLS_NAME[pol],POBS_NAME[pob],VARS_TITLE[0][var],VARS_TITLE[0][var]))
				vart.SetTextSize(0.05);
				pt.Draw()

				csavename=('%s/%s')%(outdir_seq_pol,cR2.GetName())
				cR2.SaveAs(('%s.png')%(csavename))
				print ('>>>convert %s.png %s.pdf')%(csavename,csavename)
				rc = subprocess.call(['convert', '%s.png'%csavename, '%s.pdf'%csavename])
				if rc!=0:print '.png to .pdf failed for %s'%csavename

def makepdf():
	"""For every SEQ and POL, combines R2 from every q2wbin into a single pdf.
	If R2 does not exist, it is simply ignored

	"""
	for seq in range(0,NSEQ):
		for pol in range(0,NPOLS):
			#First check to see if in the last q2wbin processed,
			#SEQ/POl exists and if so if there are any R2s plotted for SEQ/POL
			R2s_exist=False
			outdir_seq_pol=os.path.join(outdir_q2w,SEQ_NAME[seq],POLS_NAME[pol])
			if os.path.isdir(outdir_seq_pol):
				files = os.listdir(outdir_seq_pol)
  				if len(files) != 0:
					R2s_exist=True
						
			if not R2s_exist: continue
			#print 'R2s exist for', outdir_seq_pol
			outdir_pdf=os.path.join(outdir_root,SEQ_NAME[seq],POLS_NAME[pol])
			if not os.path.isdir(outdir_pdf):os.makedirs(outdir_pdf)
			for var in range(0,NVARS):
				if var==PHI or var==ALPHA: continue
				for pob in range(0,NPOBS):
					if pob==A: continue
					#Following are arguments for UNIX shell command
					pdfname=('hR2_%s_1%s.pdf')%(POBS_NAME[pob],VARS_NAME[var])
					q2w_pdfs=('%s/*/%s/%s/%s')%(outdir_root,SEQ_NAME[seq],POLS_NAME[pol],pdfname)
					out_pdf=('%s/%s')%(outdir_pdf,pdfname)

					print '>>>ls %s > /tmp/tmp'%q2w_pdfs
					rc=subprocess.call('ls %s > /tmp/tmp'%q2w_pdfs,shell=True)
					if rc!=0: print 'command failed!'

					print '>>>echo %s >> /tmp/tmp'%out_pdf
					rc=subprocess.call('echo %s >> /tmp/tmp'%out_pdf,shell=True)
					if rc!=0: print 'command failed!'

					print '>>>cat /tmp/tmp | xargs pdfunite'
					rc=subprocess.call('cat /tmp/tmp | xargs pdfunite',shell=True)
					if rc!=0: print 'command failed!'

def make1Dpdf():
	"""put all q2wbin/1D.pdf in a single pdf"""

	pdfname='1D.pdf'
	q2w_pdfs=('%s/*/%s')%(outdir_1D_root,pdfname)
	out_pdf=('%s/%s')%(outdir_1D_root,pdfname)

	print '>>>ls %s > /tmp/tmp'
	rc=subprocess.call('ls %s > /tmp/tmp'%q2w_pdfs,shell=True)
	if rc!=0: print 'command failed!'

	print '>>>echo %s >> /tmp/tmp'%out_pdf
	rc=subprocess.call('echo %s >> /tmp/tmp'%out_pdf,shell=True)
	if rc!=0: print 'command failed!'

	print '>>>cat /tmp/tmp | xargs pdfunite'
	rc=subprocess.call('cat /tmp/tmp | xargs pdfunite',shell=True)
	if rc!=0: print 'command failed!'

def makedf():
	"""For a particular Q2W analysis range, starting at the level of h5s, put all analysis
	objects into a DataFrame

	output -- DF5.h5
	contains DataFrame d with analysis data represented by columns:
	q2wbinnum -- q2wbin -- SEQ -- POL -- h5 -- h1{ij} -- h5{p} -- hR2_{p}^{ij}

	"""
	outfile = os.path.join(ANADIR,'DF5.h5')
	store = pd.HDFStore(outfile)

	norm = 50000*math.pi

	d = pd.DataFrame()
	ftemplate = root_open(SEQ_POLS_H5FILE[0][0])
	keys = ftemplate.GetListOfKeys()

	#1. First do the "looping part" of process of creating DF
	q2wbinnum=0
	dl_counter=0 #counter for number of dls (=Data-Lists, defined later) insereted into DF
	for q2wdir in keys:
		q2wbinnum+=1
		for seq in range(0,NSEQ):
			for pol in range(0,NPOLS):#poll:
				dl_counter+=1
				
				f=root_open(SEQ_POLS_H5FILE[seq][pol])
				h5=f.Get('%s/%s'%(q2wdir.GetName(),SEQ_POLS_H5[seq][pol]))
				f.Close()
									
				#Create Data-List (dl) to be added to the DataFrame			
				dl = [q2wbinnum,q2wdir.GetName(),seq,pol,h5]
				rindex=['q2wbinnum','q2wbin','SEQ','POL','h5']
				print 'len(dl)=',len(dl)
				print 'len(rindex)=',len(rindex)
				if not d:
					data = pd.DataFrame({'s1':dl},index=rindex) # Data for 1st. Column 
					d = d.append(data)
				else:
					d['s%d'%dl_counter]=dl
	
	dt = d.transpose()
	
	#2. Now use semi-vectorized operation to fill up rest of the DF5
	h5s = dt['h5']
	h1s=[]
	for i in range(len(h5s)):
		h1s.append(h5s[i].Projection(M1))
		#dt['h1_1M1']=dt['h5']
	dt['h1_1M1']=h1s

	store['d']=dt


# MAIN part of the program
#COSMETICS
gStyle.SetOptStat(0);
gStyle.SetOptFit(1111);
#INPUT data
infile = os.path.join(ANADIR,'DF5.h5')
#OUTPUT data
outdir_1D_root = os.path.join(ANADIR,'1D')
outdir_root = os.path.join(ANADIR,'polobs.new')
if not os.path.isdir(outdir_root):
	os.makedirs(outdir_root)

#First make DataFrame that will be used in the analysis
if not os.path.isfile(infile):
	makedf()
#See what d looks like
store = pd.HDFStore(infile)
d=store['d']
#print d.loc['s1':'s3','h5B':'h5C']
print d.loc['s1':'s3',:]
#Now use the DataFrame to access histograms"""
d_grpd_q2wbin=d.groupby('q2wbin')
nq2wbins = len(d_grpd_q2wbin)
print 'nq2wbins=',nq2wbins

# for q2wbin in d_grpd_q2wbin.groups:
# 		d_q2w=d_grpd_q2wbin.get_group(q2wbin)
# 		outdir_q2w=os.path.join(outdir_root,q2wbin)
# 		if not os.path.isdir(outdir_q2w):
# 			os.makedirs(outdir_q2w)
		
# 		seq_pol_sell =[
# 			(d_q2w['SEQ']==SF) & (d_q2w['POL']==UNP),
# 			(d_q2w['SEQ']==EF) & (d_q2w['POL']==UNP),
# 			(d_q2w['SEQ']==EF) & (d_q2w['POL']==POS),
# 			(d_q2w['SEQ']==EF) & (d_q2w['POL']==NEG),
# 			(d_q2w['SEQ']==EC) & (d_q2w['POL']==AVG) #dnp results
# 		]
		
		#plotR2(seq_pol_sell)
		#plot1D()

#Now put all q2wbin/seq/R2^{ij}_{alpha} in a single pdf
#makepdf()
#make1Dpdf()
