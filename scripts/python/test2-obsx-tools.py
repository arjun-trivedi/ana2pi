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

from ROOT import gSystem, gROOT, gStyle, THnSparseF, TCanvas, TString, TLine, TPaveText, TText
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
VARS_NAME = ['M1','M2','theta','phi','alpha']
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

def plotR2(h5,seql,poll):
	norm = 50000
	for var in range(0,NVARS):
		if var==PHI or var==ALPHA:continue
		for pobs in range(0,NPOBS):
			if pobs==A: continue
			for seq in seql:
				outdir_seq=os.path.join(outdir_q2w,SEQ_NAME[seq])
				if not os.path.isdir(outdir_seq):os.makedirs(outdir_seq)

				cname = ('R2_%s_1%s')%(POBS_NAME[pobs],VARS_NAME[var])
				l = TLine(0,0,180,0)
				cR2 = TCanvas(cname,cname)#"RvVar", "RvVar")
				hR2 = {}
				#pltnum=0
				for pol in poll:
					outdir_pol=os.path.join(outdir_seq,POLS_NAME[pol])
					if not os.path.isdir(outdir_pol):os.makedirs(outdir_pol)
					if pol!=AVG:
						hR2[(pol)] = h5[(seq,pol,pobs)].Projection(var)
						hR2[(pol)].Scale(1/math.pi)
						hR2[(pol)].Scale(1/norm)
						hR2[(pol)].SetLineColor(gROOT.ProcessLine("%s"%POLS_COLOR[pol]));
						hR2[(pol)].SetMarkerStyle(gROOT.ProcessLine("kFullCircle"));
					elif pol==AVG:
						hR2[(AVG)] = hR2[(POS)].Clone("avg")
						hR2[(AVG)].Add(hR2[(NEG)])
						hR2[(AVG)].Scale(0.5)
						if pobs==D:
							hR2[(AVG)].SetMinimum(-0.003)
							hR2[(AVG)].SetMaximum(0.003)			
					#Make Titles nice
					hR2[(pol)].SetTitle("")
					#hR2[(NEG)].SetTitle("")
					#hR2[(AVG)].SetTitle("")
					

					#cR2 = {}
					#l = TLine(0,0,180,0)
					#cname = ('R2%sV%s')%(VARS_NAME[var],VARS_NAME[var])
					#if pltnum==0:
					hR2[(pol)].Draw("ep")
					#else:
					#	hR2[(pol)].Draw("ep same")
					# hR2[(NEG)].Draw("ep same")
					#hR2[(AVG)].Draw("ep same")
					l.Draw("same")
					pt = TPaveText(0.3, 0.85, 0.7, 1.0, "NDC")
					q2wt = pt.AddText('[Q^{2}][W] = %s'%q2wdir.GetName())
					q2wt.SetTextColor(gROOT.ProcessLine("kBlue"))
					vart = pt.AddText(("%s^{%s} vs. %s")%(POBS_NAME[pobs],VARS_TITLE[0][var],VARS_TITLE[0][var]));
					vart.SetTextSize(0.05);
					pt.Draw()
					csavename = ('%s/%s')%(outdir_pol,cR2.GetName())
					cR2.SaveAs( ('%s.png')%(csavename))
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
				for pobs in range(0,NPOBS):
					if pobs==A: continue
					#Following are arguments for UNIX shell command
					pdfname=('R2_%s_1%s.pdf')%(POBS_NAME[pobs],VARS_NAME[var])
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

q2wbinnum=0
counter=0
for q2wdir in keys:
	q2wbinnum+=1
	for seq in seql:
		for pol in poll:
			if pol==AVG: continue
			counter+=1	
			f=root_open(SEQ_POLS_H5FILE[seq][pol])
			#h5[(seq,pol)]=f.Get('%s/%s'%(q2wdir.GetName(),SEQ_POLS_H5[seq][pol]))
			h5=f.Get('%s/%s'%(q2wdir.GetName(),SEQ_POLS_H5[seq][pol]))
			f.Close()
			h5B=mythnt.MultiplyBy(h5,'cphi')
			hR2_B_1M1 = h5B.Projection(M1)

			l = [q2wbinnum,q2wdir.GetName(),SEQ_NAME[seq],POLS_NAME[pol],h5,h5B,hR2_B_1M1]
			rindex=['q2wbinnum','q2wbin','SEQ','POL','h5','h5B','hR2_B_1M1']
			if not d:
				data = pd.DataFrame({'s1':l},index=rindex) # Data for 1st. Column 
				d = d.append(data)
			else:
				d['s%d'%counter]=l
#print a few columns of d to see what it looks like
print 'dataframe = '
print d.loc[:,'s1':'s4']
# Test Drawing a histogram from d
# c = TCanvas("test","test")
# d['s2'][6].Draw("ep")
# c.SaveAs("test.png")

dt = d.transpose()
print dt.loc['s1':'s4',:]
dt_grpd_q2wbinnum=dt.groupby('q2wbinnum')
nq2wbins = len(dt_grpd_q2wbinnum)
print 'nq2wbins=',nq2wbins

for q2wbin in dt_grpd_q2wbinnum.groups:
	#print q2wbin
	df=dt_grpd_q2wbinnum.get_group(q2wbin)
	sel = (df['SEQ']=='EC') & (df['POL']=='POS')
	print df[sel]


			

	#plotR2(h5,seql,poll)

	
"""
Now put all q2wbin/seq/R2^{ij}_{alpha} in a single pdf
"""
#makepdf(seql,poll)
