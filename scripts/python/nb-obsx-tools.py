# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

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

# <codecell>

ftemplate = root_open(SEQ_POLS_H5FILE[0][0])
keys = ftemplate.GetListOfKeys()

q2wdirnum=0
for q2wdir in keys:
    q2wdirnum+=1;
    for seq in range(0,NSEQ):
        for pol in range(0,NPOLS):
            if pol==AVG: continue
            #print 'q2wdir:seq:pol=%s:%s:%s'%(q2wdir.GetName(),
             #                                SEQ_NAME[seq],
              #                               POLS_NAME[pol])
            f=root_open(SEQ_POLS_H5FILE[seq][pol])
            s=[q2wdirnum,q2wdir.GetName(),SEQ_NAME[seq],POLS_NAME[pol]]
            f.Close()

# <codecell>


