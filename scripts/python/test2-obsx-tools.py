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
	norm = 50000
	for seq_pol_sel in seq_pol_sell:
		dq2w_seq_pol = dq2w[seq_pol_sel]
		print dq2w_seq_pol
		# seq = dq2w_seq_pol['SEQ']
		# print 'seq=',seq
		# hR2 = dq2w_seq_pol['hR2_B_1THETA']
		# print 'hR2=',hR2
		# hR2.Draw()

		# tl = dq2w_seq_pol.values
		# print 'tl='
		# print tl
		# print tl[0][6]
		# h = tl[0][6]
		# print h.GetName()
		outdir_seq_pol=os.path.join(outdir_q2w,dq2w_seq_pol.iloc[0]['SEQ'],dq2w_seq_pol.iloc[0]['POL'])
		if not os.path.isdir(outdir_seq_pol):
			os.makedirs(outdir_seq_pol)

		hR2 = dq2w_seq_pol.iloc[0]['hR2_D_1THETA']
		cR2 = TCanvas("R2_D_1THETA","R2_D_1THETA")
		hR2.Draw("ep")
		csavename=('%s/%s')%(outdir_seq_pol,cR2.GetName())
		cR2.SaveAs(('%s.png')%(csavename))
		print ('>>>convert %s.png %s.pdf')%(csavename,csavename)
		rc = subprocess.call(['convert', '%s.png'%csavename, '%s.pdf'%csavename])
		if rc!=0:print '.png to .pdf failed for %s'%csavename

	

def makepdf(seql,poll):
	for seq in seql:
		for pol in poll:
			outdir_pdf=os.path.join(outdir_root,SEQ_NAME[seq],POLS_NAME[pol])
			if not os.path.isdir(outdir_pdf):os.makedirs(outdir_pdf)
			for var in range(0,NVARS):
				if var==PHI or var==ALPHA:continue
				for pob in range(0,NPOBS):
					if pob==A: continue
					#Following are arguments for UNIX shell command
					pdfname=('R2_%s_1%s.pdf')%(POBS_NAME[pob],VARS_NAME[var])
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


#COSMETICS
gStyle.SetOptStat(0);
gStyle.SetOptFit(1111);
#OUTPUT data
outdir_root = os.path.join(ANADIR,'polobs.new')
if not os.path.isdir(outdir_root):
	os.makedirs(outdir_root)

#MAIN part of program
# Using rootpy
# for dirs in f.walk('.').next()[1]:
# 	h5 = {}
# 	print dirs
# 	obj = os.path.join('.',dirs,'hY5D_POS','Varset1','hY5D_ACC_CORR')
# 	print 'obj=',obj
# 	h5[(EC,POS)] = obj
# 	mythnt.MultiplyBySinPhi(h5[(EC,POS)]);

#Using PyRoot
"""
First, for each q2wbin/seql, make all R2^{ij}_{alpha} where:
	-ij    = Index over NVARSETS and nVAR respectively
	-alpha = Index over NPOBS (=[A,B,C,D])
"""
#1.
# seql = [EF,SF]
# poll = [UNP]

#2.
seql = [EC]
poll = [POS,NEG,AVG]

d = pd.DataFrame()
ftemplate = root_open(SEQ_POLS_H5FILE[0][0])
keys = ftemplate.GetListOfKeys()

"""First create the DataFrame that will be used in the analysis"""
q2wbinnum=0
counter=0
for q2wdir in keys:
	q2wbinnum+=1
	for seq in seql:
		# hR2_B_1THETA={} #indexed by pol
		# hR2_C_1THETA={}
		# hR2_D_1THETA={}
		#dhists=[]
		h5m={}#indexed by POLS,POBS
		hR2={}#indexed by POLS,POBS,VARSETS,VARS
		norm = 50000*math.pi
		for pol in poll:
			dhists=[]
			counter+=1
			if pol!=AVG:
				f=root_open(SEQ_POLS_H5FILE[seq][pol])
				h5=f.Get('%s/%s'%(q2wdir.GetName(),SEQ_POLS_H5[seq][pol]))
				f.Close()
				
				for pob in range(0,NPOBS):
					if   pob==A: continue
					elif pob==B: h5m[(pol,pob)]=mythnt.MultiplyBy(h5,'cphi')
					elif pob==C: h5m[(pol,pob)]=mythnt.MultiplyBy(h5,'c2phi')
					elif pob==D:
						if   pol==POS or pol==UNP:
							h5m[(pol,pob)]=mythnt.MultiplyBy(h5,'sphi')
						elif pol==NEG:
							h5m[(pol,pob)]=mythnt.MultiplyBy(h5,'sphi',-1)
					dhists.append(h5m[(pol,pob)])	
					for varset in range(0,NVARSETS):
						if varset==1 or varset==2: continue
						for var in range(0,NVARS):
							hR2[(pol,pob,varset,var)]=h5m[(pol,pob)].Projection(var)
							hR2[(pol,pob,varset,var)].Scale(1/norm)
							dhists.append(hR2[(pol,pob,varset,var)])
			elif pol==AVG:
				for pob in range(0,NPOBS):
					if pob==A: continue
					dhists.append('')
					for varset in range(0,NVARSETS):
						if varset==1 or varset==2: continue
						for var in range(0,NVARS):
							hR2[(AVG,pob,varset,var)]=hR2[(POS,pob,varset,var)].Clone("avg")
							hR2[(AVG,pob,varset,var)].Add(hR2[(NEG,pob,varset,var)])
							hR2[(AVG,pob,varset,var)].Scale(0.5)
							if pob==D:
								hR2[(AVG,pob,varset,var)].SetMinimum(-0.003)
								hR2[(AVG,pob,varset,var)].SetMaximum(0.003)
							dhists.append(hR2[(AVG,pob,varset,var)])		
				# hR2_B_1THETA[AVG] = hR2_B_1THETA[POS].Clone("avg")
				# hR2_B_1THETA[AVG].Add(hR2_B_1THETA[NEG])
				# hR2_B_1THETA[AVG].Scale(0.5)
				
				# hR2_C_1THETA[AVG] = hR2_C_1THETA[POS].Clone("avg")
				# hR2_C_1THETA[AVG].Add(hR2_C_1THETA[NEG])
				# hR2_C_1THETA[AVG].Scale(0.5)
				
				# hR2_D_1THETA[AVG] = hR2_D_1THETA[POS].Clone("avg")
				# hR2_D_1THETA[AVG].Add(hR2_D_1THETA[NEG])
				# hR2_D_1THETA[AVG].Scale(0.5)
				# hR2_D_1THETA[AVG].SetMinimum(-0.003)
				# hR2_D_1THETA[AVG].SetMaximum(0.003)
			
			# dl = [q2wbinnum,q2wdir.GetName(),SEQ_NAME[seq],POLS_NAME[pol],
			# h5B,hR2_B_1THETA[pol],h5C,hR2_C_1THETA[pol],h5D,hR2_D_1THETA[pol]]
			dl = [q2wbinnum,q2wdir.GetName(),SEQ_NAME[seq],POLS_NAME[pol]]
			dl+=dhists
			print 'dl='
			print dl
			print len(dl)
			rindex=['q2wbinnum','q2wbin','SEQ','POL',
			'h5B','hR2_B_1M1','hR2_B_1M2','hR2_B_1THETA','hR2_B_1PHI','hR2_B_1ALPHA',
			'h5C','hR2_C_1M1','hR2_C_1M2','hR2_C_1THETA','hR2_C_1PHI','hR2_C_1ALPHA',
			'h5D','hR2_D_1M1','hR2_D_1M2','hR2_D_1THETA','hR2_D_1PHI','hR2_D_1ALPHA']
			print len(rindex)
			if not d:
				data = pd.DataFrame({'s1':dl},index=rindex) # Data for 1st. Column 
				d = d.append(data)
			else:
				d['s%d'%counter]=dl


#print a few columns of d to see what it looks like
print 'dataframe = '
print d.loc[:,'s1':'s4']
# Test Drawing a histogram from d
# c = TCanvas("test","test")
# d['s2'][6].Draw("ep")
# c.SaveAs("test.png")

dt = d.transpose()
print dt.loc['s1':'s3','h5B':'h5C']
dt_grpd_q2wbin=dt.groupby('q2wbin')
nq2wbins = len(dt_grpd_q2wbin)
print 'nq2wbins=',nq2wbins


"""Now use the DataFrame to access histograms"""
for q2wbin in dt_grpd_q2wbin.groups:
		dq2w=dt_grpd_q2wbin.get_group(q2wbin)
		outdir_q2w=os.path.join(outdir_root,q2wbin)
		if not os.path.isdir(outdir_q2w):
			os.makedirs(outdir_q2w)
		# seq_pol_sell = [
		# 			(dq2w['SEQ']=='EC') & (dq2w['POL']=='POS'),
		# 			(dq2w['SEQ']=='EC') & (dq2w['POL']=='NEG'),
		# 	   ]
		seq_pol_sell = [(dq2w['SEQ']=='EC') & (dq2w['POL']=='AVG')]
		plotR2(seq_pol_sell)

	
"""
Now put all q2wbin/seq/R2^{ij}_{alpha} in a single pdf
"""
makepdf(seql,poll)
