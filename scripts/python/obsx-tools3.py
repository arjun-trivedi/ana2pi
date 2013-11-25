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

# NSEQ_SIM=6
# ST,SR,SA,SC,SH,SF=range(NSEQ_SIM)

# NSEQ_EXP=3
# ER,EC,EF=range(NSEQ_EXP)

NSEQ=9
ST,SR,SA,SC,SH,SF,ER,EC,EF=range(NSEQ)
SEQ_NAME=['ST','SR','SA','SC','SH','SF','ER','EC','EF']

NPOLS=4
POS,NEG,UNP,AVG=range(NPOLS)
POLS_COLOR=['kRed','kBlack','kBlue','kMagenta']
POLS_NAME=['POS','NEG','UNP','AVG']

ANADIR = os.environ['E1F_2PI_ANADIR2']
TOP = "1:2:3:4"
Q2WBNG = '1-2.000-2.400__24-1.300-1.900'
SIMDIR = os.path.join(ANADIR,'simdir')
SEQ_POLS_H5FILE = np.zeros((NSEQ,NPOLS),object)
SEQ_POLS_H5FILE[ST:SF+1,POS:NEG+1]='%s/%s__%s__pol__sim.root'%(SIMDIR,TOP,Q2WBNG)
SEQ_POLS_H5FILE[ST:SF+1,UNP:AVG]  ='%s/%s__%s__sim.root'%(SIMDIR,TOP,Q2WBNG)
SEQ_POLS_H5FILE[ER:,POS:NEG+1]    ='%s/%s__%s__pol__exp.root'%(ANADIR,TOP,Q2WBNG)
SEQ_POLS_H5FILE[ER:,UNP:AVG]      ='%s/%s__%s__exp.root'%(ANADIR,TOP,Q2WBNG)
#print SEQ_POLS_H5FILE

SEQ_POLS_H5=[
	['hY5D/Varset1/hY5D_TH','hY5D/Varset1/hY5D_TH','hY5D/Varset1/hY5D_TH',''],
	['hY5D/Varset1/hY5D_RECO','hY5D/Varset1/hY5D_RECO','hY5D/Varset1/hY5D_RECO',''],
	['hY5D/Varset1/hY5D_ACC','hY5D/Varset1/hY5D_ACC','hY5D/Varset1/hY5D_ACC',''],
	['hY5D/Varset1/hY5D_ACC_CORR','hY5D/Varset1/hY5D_ACC_CORR','hY5D/Varset1/hY5D_ACC_CORR',''],
	['hY5D/Varset1/hY5D_HOLE','hY5D/Varset1/hY5D_HOLE','hY5D/Varset1/hY5D_HOLE',''],
	['hY5D/Varset1/hY5D_FULL','hY5D/Varset1/hY5D_FULL','hY5D/Varset1/hY5D_FULL',''],
	['hY5D_POS/Varset1/hY5D_RECO',    'hY5D_NEG/Varset1/hY5D_RECO',    'hY5D/Varset1/hY5D_RECO',''],
	['hY5D_POS/Varset1/hY5D_ACC_CORR','hY5D_NEG/Varset1/hY5D_ACC_CORR','hY5D/Varset1/hY5D_ACC_CORR',''],
	['hY5D_POS/Varset1/hY5D_FULL',    'hY5D_NEG/Varset1/hY5D_FULL',    'hY5D/Varset1/hY5D_FULL','']
]

NPOBS=4
A,B,C,D=range(NPOBS) # [A,B,C,D]=[Rt+Rl,Rlt,Rtt,Rlt']
POBS_NAME=['A','B','C','D']
#print var

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

def makedf():
	"""For a particular Q2W analysis range, starting at the level of h5s, put all analysis
	objects into a DataFrame

	output -- ana_store.h5
	contains DataFrame d with analysis data represented by columns:
	q2wbinnum * q2wbin * SEQ * POL * h5 * h1{ij} * h5{p} * hR2_{p}^{ij}

	"""
	outfile = os.path.join(ANADIR,'ana_store.h5')
	store = pd.HDFStore(outfile)

	# seql = [EC]
	# poll = [POS,NEG,AVG]

	d = pd.DataFrame()
	ftemplate = root_open(SEQ_POLS_H5FILE[0][0])
	keys = ftemplate.GetListOfKeys()

	q2wbinnum=0
	dl_counter=0 #counter for number of dls (=Data-Lists, defined later) insereted into DF
	for q2wdir in keys:
		q2wbinnum+=1
		for seq in range(0,NSEQ):
			# The following indexing is required for reason of
			# *practicality* and *ease-of-reading* the code
			# 	- *ease-of-reading* code requires POLS
			#	- *Practical* reason requires POBS,VARSETS,VARS since
			#	  after obtain hR2(POLS=POS,NEG,UNP), hR2(POL=AVG) needs to be 
			#	  calculated for all VARSETS,VARS and this is done in a seperate loop
			h5={}#indexed by POLS
			h1={}#indexed by POLS
			h5p={}#indexed by POLS,POBS
			hR2={}#indexed by POLS,POBS,VARSETS,VARS
			norm = 50000*math.pi
			for pol in range(0,NPOLS):#poll:
				hists=[]
				dl_counter+=1
				if pol!=AVG:
					f=root_open(SEQ_POLS_H5FILE[seq][pol])
					h5[pol]=f.Get('%s/%s'%(q2wdir.GetName(),SEQ_POLS_H5[seq][pol]))
					f.Close()
					hists.append(h5[pol])
					
					#1. do 1D obs
					for varset in range(0,NVARSETS):
							if varset==1 or varset==2: continue
							for var in range(0,NVARS):
								h1[pol]=h5[pol].Projection(var)
								hists.append(h1[pol])
					#2. do Pol. Obs.			
					for pob in range(0,NPOBS):
						if   pob==A: continue
						elif pob==B: h5p[(pol,pob)]=mythnt.MultiplyBy(h5[pol],'cphi')
						elif pob==C: h5p[(pol,pob)]=mythnt.MultiplyBy(h5[pol],'c2phi')
						elif pob==D:
							if   pol==POS or pol==UNP:
								h5p[(pol,pob)]=mythnt.MultiplyBy(h5[pol],'sphi')
							elif pol==NEG:
								h5p[(pol,pob)]=mythnt.MultiplyBy(h5[pol],'sphi',-1)
						hists.append(h5p[(pol,pob)])	
						for varset in range(0,NVARSETS):
							if varset==1 or varset==2: continue
							for var in range(0,NVARS):
								hR2[(pol,pob,varset,var)]=h5p[(pol,pob)].Projection(var)
								hR2[(pol,pob,varset,var)].Scale(1/norm)
								hists.append(hR2[(pol,pob,varset,var)])
				elif pol==AVG:
					hists.append('') #blank version of h5[pol]
					#1. do 1D obs
					for varset in range(0,NVARSETS):
							if varset==1 or varset==2: continue
							for var in range(0,NVARS):
								hists.append('') #blank version of h1_ij
					#2. do Pol.Obs
					for pob in range(0,NPOBS):
						if pob==A: continue
						hists.append('') #blank version of h5p[pol]
						for varset in range(0,NVARSETS):
							if varset==1 or varset==2: continue
							for var in range(0,NVARS):
								hR2[(AVG,pob,varset,var)]=hR2[(POS,pob,varset,var)].Clone("avg")
								hR2[(AVG,pob,varset,var)].Add(hR2[(NEG,pob,varset,var)])
								hR2[(AVG,pob,varset,var)].Scale(0.5)
								if pob==D:
									hR2[(AVG,pob,varset,var)].SetMinimum(-0.003)
									hR2[(AVG,pob,varset,var)].SetMaximum(0.003)
								hists.append(hR2[(AVG,pob,varset,var)])		
				"""Create Data-List (dl) to be added to the DataFrame"""			
				#dl = [q2wbinnum,q2wdir.GetName(),SEQ_NAME[seq],POLS_NAME[pol]]
				dl = [q2wbinnum,q2wdir.GetName(),seq,pol]
				dl+=hists
				rindex=['q2wbinnum','q2wbin','SEQ','POL',
				'h5', 'h1_1M1',   'h1_1M2',   'h1_1THETA',   'h1_1PHI',   'h1_1ALPHA',
				'h5B','hR2_B_1M1','hR2_B_1M2','hR2_B_1THETA','hR2_B_1PHI','hR2_B_1ALPHA',
				'h5C','hR2_C_1M1','hR2_C_1M2','hR2_C_1THETA','hR2_C_1PHI','hR2_C_1ALPHA',
				'h5D','hR2_D_1M1','hR2_D_1M2','hR2_D_1THETA','hR2_D_1PHI','hR2_D_1ALPHA']
				print 'len(dl)=',len(dl)
				print 'len(rindex)=',len(rindex)
				# print 'dl='
				# print dl
				# print 'rindex='
				# print rindex
				if not d:
					data = pd.DataFrame({'s1':dl},index=rindex) # Data for 1st. Column 
					d = d.append(data)
				else:
					d['s%d'%dl_counter]=dl
	
	dt = d.transpose()
	store['d']=dt
	#return dt					

# MAIN part of the program
#COSMETICS
gStyle.SetOptStat(0);
gStyle.SetOptFit(1111);
#INPUT data
infile = os.path.join(ANADIR,'ana_store.h5')
#OUTPUT data
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

for q2wbin in d_grpd_q2wbin.groups:
		d_q2w=d_grpd_q2wbin.get_group(q2wbin)
		outdir_q2w=os.path.join(outdir_root,q2wbin)
		if not os.path.isdir(outdir_q2w):
			os.makedirs(outdir_q2w)
		
		seq_pol_sell =[
			(d_q2w['SEQ']==SF) & (d_q2w['POL']==UNP),
			(d_q2w['SEQ']==EF) & (d_q2w['POL']==UNP),
			(d_q2w['SEQ']==EF) & (d_q2w['POL']==POS),
			(d_q2w['SEQ']==EF) & (d_q2w['POL']==NEG),
			(d_q2w['SEQ']==EC) & (d_q2w['POL']==AVG) #dnp results
		]
		
		#plotR2(seq_pol_sell)

#Now put all q2wbin/seq/R2^{ij}_{alpha} in a single pdf
#makepdf()
