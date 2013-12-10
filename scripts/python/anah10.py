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
E1F_P = 5.499; 
lvE0 = ROOT.TLorentzVector(0,0,E1F_P,math.sqrt(E1F_P*E1F_P+MASS_E*MASS_E));
lvP0 = ROOT.TLorentzVector(0,0,0,MASS_P);

EVT_COLS = ['gpart'] 
EVT_SUB_COLS=['p','cx', 'cy',  'cz','etot','ec_ei','ec_eo','sc_t']
NBINS =      [150, 100,  100,  100,  100,   100,    100,    100]
XMIN=        [2.0, -1.0, -1.0, 0.9,  0.0,   0.0,    0.0,    0.0]
XMAX=        [4.0, 1.0,  1.0,  1.0,  1.0,   1.0,    1.0,    30]
#XMIN=        [0.0, -1.0, -1.0, 0.8,  0.0,   0.0,    0.0,    0.0]
#XMAX=        [5.0, 1.0,  1.0,  1.0,  1.0,   1.0,    1.0,    30]

MC_COLS=['mcnentr','mcid','mcp','mctheta','mcphi']

COLS=EVT_COLS + EVT_SUB_COLS + MC_COLS

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
        for col in EVT_SUB_COLS:
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
    
    gs = gridspec.GridSpec(2,2)
    
    ax_w=plt.subplot(gs[0])
    ax_w.set_title('rec-W')
    ax_w.set_xlabel('W(GeV)')
    for h10df,c,l in zip(h10dfs,colors,labels):
        plt.hist(h10df.w,100,(1.0,2.0),
            histtype='step',color=c,label='rec_%s'%l)
    
    ax_qsq=plt.subplot(gs[1])
    ax_qsq.set_title('rec-Q2')
    ax_qsq.set_xlabel('Q2(GeV^2)')
    for h10df,c,l in zip(h10dfs,colors,labels):
        plt.hist(h10df.qsq,100,(1.0,2.0),
            histtype='step',color=c,label='rec_%s'%l)
        plt.legend(loc=2)

    ax_mcw=plt.subplot(gs[2])
    ax_mcw.set_title('mc-W')
    ax_mcw.set_xlabel('W(GeV)')
    for h10df,c,l in zip(h10dfs,colors,labels):
        plt.hist(h10df.mc_w,100,(1.0,2.0),
            histtype='step',color=c,label='mc_%s'%l,linestyle='dashed')
        

    ax_mcqsq=plt.subplot(gs[3])
    ax_mcqsq.set_title('mc-Q2')
    ax_mcqsq.set_xlabel('Q2(GeV^2)')
    for h10df,c,l in zip(h10dfs,colors,labels):
        plt.hist(h10df.mc_qsq,100,(1.0,2.0),
            histtype='step',color=c,label='mc_%s'%l,linestyle='dashed')
        plt.legend(loc=2)

def plot_evt_sub_cols(dfs):
    #set up some cosmetics
    colors = cm.rainbow(np.linspace(0,1,len(dfs)))
    labels = ['h10_%d'%i for i in range(len(dfs))]
    
    gs = gridspec.GridSpec(len(EVT_SUB_COLS)/4,len(EVT_SUB_COLS)/2)
    for icol in np.arange(0,len(EVT_SUB_COLS)):
        ax=plt.subplot(gs[icol])
        ax.set_title(EVT_SUB_COLS[icol])
        ax.set_xlabel(EVT_SUB_COLS[icol])
        for c,l,df in zip(colors,labels,dfs):
            plt.hist(df['el_%s'%EVT_SUB_COLS[icol]],NBINS[icol],[XMIN[icol],XMAX[icol]],
                histtype='step',color=c,label=l)
    plt.legend()     