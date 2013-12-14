#!/usr/bin/python
import ROOT
from root_numpy import root2array, root2rec, tree2rec
from root_numpy.testdata import get_filepath

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.gridspec as gridspec

import math

MASS_E = 0.000511
MASS_P = 0.93827203
E1F_P = 5.499#5.479#5.499; 
lvE0 = ROOT.TLorentzVector(0,0,E1F_P,math.sqrt(E1F_P*E1F_P+MASS_E*MASS_E));
lvP0 = ROOT.TLorentzVector(0,0,0,MASS_P);

#MAIN_COLS = ['gpart','vx','vy','vz'] 
MAIN = {}
MAIN['COLS']=['gpart'] 


DC={}
DC['COLS']=['p','cx', 'cy',  'cz','vx','vy','vz']
EL_DC={}
EL_DC['COLS']=['el_p','el_cx', 'el_cy',  'el_cz','el_vx','el_vy','el_vz']
EL_DC['NBINS'] =      [150, 100,  100,  100, 100,100,100]
EL_DC['XMIN']=        [2.0, -1.0, -1.0, 0.9, -2, -4,  -40]
EL_DC['XMAX']=        [4.0, 1.0,  1.0,  1.0, 2,  4,  0]
# DC_COLS=['p','cx', 'cy',  'cz']
# NBINS_DC_COLS =      [150, 100,  100,  100]
# XMIN_DC_COLS=        [2.0, -1.0, -1.0, 0.9]
# XMAX_DC_COLS=        [4.0, 1.0,  1.0,  1.0]

EC={}
EC['COLS']=['etot','ec_ei','ec_eo']
EL_EC={}
EL_EC['COLS']=['el_etot','el_ec_ei','el_ec_eo']
EL_EC['NBINS'] =      [100, 100, 100]
EL_EC['XMIN']=        [0.0, 0.0, 0.0]
EL_EC['XMAX']=        [1.0, 1.0, 1.0]
# EC_COLS=['etot','ec_ei','ec_eo']
# NBINS_EC_COLS =      [100, 100, 100]
# XMIN_EC_COLS=        [0.0, 0.0, 0.0]
# XMAX_EC_COLS=        [1.0, 1.0, 1.0]

SC={}
SC['COLS']=['sc_t']
EL_SC={}
EL_SC['COLS']=['el_sc_t']
EL_SC['NBINS'] =      [100]
EL_SC['XMIN']=        [0.0]
EL_SC['XMAX']=        [30]
# SC_COLS=['sc_t']
# NBINS_SC_COLS =      [100]
# XMIN_SC_COLS=        [0.0]
# XMAX_SC_COLS=        [30]

#XMIN=        [0.0, -1.0, -1.0, 0.8,  0.0,   0.0,    0.0,    0.0]
#XMAX=        [5.0, 1.0,  1.0,  1.0,  1.0,   1.0,    1.0,    30]

MC={}
MC['COLS']=['mcnentr','mcid','mcp','mctheta','mcphi']
#MC_COLS=['mcnentr','mcid','mcp','mctheta','mcphi']


#COLS=MAIN_COLS + DC_COLS + EC_COLS + SC_COLS + MC_COLS
COLS=MAIN['COLS'] + DC['COLS'] + EC['COLS'] + SC['COLS'] + MC['COLS']


NENTRIES=2000
def import_h10(rfiles,nentries=NENTRIES,cols=COLS):
    """"From each of the rfiles:
            1. Extract and Convert h10-TTree to h10-DF and put
              in a List(h10-DFs).
            2. New Columns are added to the h10-DFs
            3. The h10-DFs are returned
    """
    
    h10dfs = []
    for rfile in rfiles:        
        ar = root2array(rfile,'h10',stop=nentries,branches=cols)
        df = pd.DataFrame(ar)
        h10dfs.append(df)

    add_cols(h10dfs)        
    return h10dfs;
    
def add_cols(dfs):
    for df in dfs:
        for col in DC['COLS']:
            add_elcol(col,df)
            add_elcol(col,df)
        for col in EC['COLS']:
            add_elcol(col,df)
            add_elcol(col,df)
    
        add_reco_elcols(df)
        add_reco_elcols(df)

def add_elcol(col,df):
    d=[]
    for ievt in np.arange(0,len(df[col])):
        if len(df[col][ievt])==0:        
            d.append(-9999)
        else:
            d.append(df[col][ievt][0])
    df['el_%s'%(col)]=d

def add_reco_elcols(d):
    qsq = []
    w = []
    for ievt in np.arange(0,len(d.el_p)):
        p = d.el_p[ievt]
        px = p * d.el_cx[ievt]
        py = p * d.el_cy[ievt]
        pz = p * d.el_cz[ievt]
        lve = ROOT.TLorentzVector(0,0,0,0)
        e = math.sqrt(p*p+MASS_E*MASS_E) 
        lve.SetPxPyPzE(px,py,pz,e)
        lvq = lvE0-lve
        lvw = lvq+lvP0
        qsq.append(-lvq.Mag2())
        w.append(lvw.Mag())
    d['qsq']=qsq 
    d['w']=w

    mc_qsq=[]
    mc_w=[]
    for ievt in np.arange(0,len(d.mcnentr)):
        for ipart in np.arange(0,d.mcnentr[ievt]):
            if d.mcid[ievt][ipart] != 11: continue
            p=d.mcp[ievt][ipart]
            theta=math.radians(d.mctheta[ievt][ipart])
            phi=math.radians(d.mcphi[ievt][ipart])
            px=p*math.sin(theta)*math.cos(phi)
            py=p*math.sin(theta)*math.sin(phi)
            pz=p*math.cos(theta)
            e=math.sqrt(p*p+MASS_E*MASS_E)
            lve = ROOT.TLorentzVector(px,py,pz,e)
            lvq = lvE0-lve
            lvw = lvq+lvP0
            mc_qsq.append(-lvq.Mag2())
            mc_w.append(lvw.Mag())
    d['mc_qsq']=mc_qsq 
    d['mc_w']=mc_w


def plot_qsq_w(h10dfs):
    #cosmetics
    colors = cm.rainbow(np.linspace(0,1,len(h10dfs)))
    labels = ['%d'%i for i in range(len(h10dfs))]
    #bng
    wlow=1.2
    whigh=2.2
    wbins=100
    q2low=1.5
    q2high=3.0
    q2bins=100
    
    gs = gridspec.GridSpec(2,2)
    
    ax_w=plt.subplot(gs[0])
    ax_w.set_title('rec-W')
    ax_w.set_xlabel('W(GeV)')
    for h10df,c,l in zip(h10dfs,colors,labels):
        plt.hist(h10df.w,wbins,(wlow,whigh),
            histtype='step',color=c,label='rec_%s'%l)
        plt.hist(h10df.mc_w,wbins,(wlow,whigh),
            histtype='step',color='green',label='mc_%s'%l,linestyle='dashed')
    
    ax_qsq=plt.subplot(gs[1])
    ax_qsq.set_title('rec-Q2')
    ax_qsq.set_xlabel('Q2(GeV^2)')
    for h10df,c,l in zip(h10dfs,colors,labels):
        plt.hist(h10df.qsq,q2bins,(q2low,q2high),
            histtype='step',color=c,label='rec_%s'%l)
        plt.hist(h10df.mc_qsq,q2bins,(q2low,q2high),
            histtype='step',color='green',label='mc_%s'%l,linestyle='dashed')
        plt.legend(loc=2,prop={'size':8})

    ax_mcw=plt.subplot(gs[2])
    ax_mcw.set_title('mc-W')
    ax_mcw.set_xlabel('W(GeV)')
    for h10df,c,l in zip(h10dfs,colors,labels):
        plt.hist(h10df.mc_w,wbins,(wlow,whigh),
            histtype='step',color=c,label='mc_%s'%l,linestyle='dashed')
        

    ax_mcqsq=plt.subplot(gs[3])
    ax_mcqsq.set_title('mc-Q2')
    ax_mcqsq.set_xlabel('Q2(GeV^2)')
    for h10df,c,l in zip(h10dfs,colors,labels):
        plt.hist(h10df.mc_qsq,q2bins,(q2low,q2high),
            histtype='step',color=c,label='mc_%s'%l,linestyle='dashed')
        plt.legend(loc=2,prop={'size':8})

def hists(dgrp,dfs):
    "dgrp = MAIN,DC,EC,SC,MC"
    #set up some cosmetics
    colors = cm.rainbow(np.linspace(0,1,len(dfs)))
    labels = ['h10_%d'%i for i in range(len(dfs))]
    
    ncols=len(dgrp['COLS'])
    gs =''
    if ncols > 4: 
        gs = gridspec.GridSpec(2,4)
    else: 
        gs = gridspec.GridSpec(1,ncols)

    for icol in np.arange(ncols):
        ax=plt.subplot(gs[icol])
        col=dgrp['COLS'][icol]
        nbins=dgrp['NBINS'][icol]
        xmin=dgrp['XMIN'][icol]
        xmax=dgrp['XMAX'][icol]
        #print 'col=%s:nbins=%d:xmin=%d:xmax%d'%(col,nbins,xmin,xmax)
        ax.set_title(col)
        ax.set_xlabel(col)
        for c,l,df in zip(colors,labels,dfs):
            #print df[col]
            plt.hist(df[col],nbins,(xmin,xmax), histtype='step',color=c,label=l)
    plt.legend(loc=2,prop={'size':8})

def plot_DC_COLS(dfs):
    #set up some cosmetics
    colors = cm.rainbow(np.linspace(0,1,len(dfs)))
    labels = ['h10_%d'%i for i in range(len(dfs))]
    
    gs = gridspec.GridSpec(len(DC_COLS)/4,len(DC_COLS)/2)
    for icol in np.arange(0,len(DC_COLS)):
        ax=plt.subplot(gs[icol])
        ax.set_title(DC_COLS[icol])
        ax.set_xlabel(DC_COLS[icol])
        for c,l,df in zip(colors,labels,dfs):
            plt.hist(df['el_%s'%DC_COLS[icol]],NBINS[icol],[XMIN[icol],XMAX[icol]],
                histtype='step',color=c,label=l)
    plt.legend()     