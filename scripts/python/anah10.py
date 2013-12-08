#!/usr/bin/python
import ROOT
from root_numpy import root2array, root2rec, tree2rec
from root_numpy.testdata import get_filepath

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#import matplotlib.pyplot as plt

import math

MASS_E = 0.000511
MASS_P = 0.93827203
E1F_P = 5.499; 
lvE0 = ROOT.TLorentzVector(0,0,E1F_P,math.sqrt(E1F_P*E1F_P+MASS_E*MASS_E));
lvP0 = ROOT.TLorentzVector(0,0,0,MASS_P);

EVT_COLS = ['gpart'] 
EVT_SUB_COLS=['p','cx','cy','cz','etot','ec_ei','ec_eo']
XMIN=[0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0]
XMAX=[5.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

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

    addcols(h10dfs)        
    return h10dfs;
    
def addcols(dfs):
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


def comp_qsq_w():
    fig = plt.figure(figsize=(8,8))
    fig.suptitle('old-new sim comparison: Reconstructed Q2,W')
    ax_w=fig.add_subplot(211)
    ax_w.set_title('W')
    ax_w.set_xlabel('W(GeV)')
    plt.hist(pnl_h10['OLD'].w,100,(1.0,2.0),
         histtype='step',color='black',label='old')
    plt.hist(pnl_h10['NEW'].w,100,(1.0,2.0),
         histtype='step',color='red',label='new')
    plt.hist(pnl_h10['OLD'].mc_w,100,(1.0,2.0),
         histtype='step',color='green',linestyle='solid',label='mc_old')
    plt.hist(pnl_h10['NEW'].mc_w,100,(1.0,2.0),
         histtype='step',color='green',linestyle='dashed',label='mc_new')

    ax_qsq=fig.add_subplot(212)
    ax_qsq.set_title('Q2')
    ax_qsq.set_xlabel('Q2(GeV)')
    plt.hist(pnl_h10['OLD'].qsq,100,(1.0,2.0),
         histtype='step',color='black',label='old')
    plt.hist(pnl_h10['NEW'].qsq,100,(1.0,2.0),
         histtype='step',color='red',label='new')
    plt.hist(pnl_h10['OLD'].mc_qsq,100,(1.0,2.0),
         histtype='step',color='green',linestyle='solid',label='mc_old')
    plt.hist(pnl_h10['NEW'].mc_qsq,100,(1.0,2.0),
         histtype='step',color='green',linestyle='dashed',label='mc_new')
    plt.legend(loc=2)

def comp_basic(dfs):
    #set up some cosmetics
    colors = cm.rainbow(np.linspace(0,1,len(dfs)))
    
    fig=plt.figure(figsize=(10,20))
    fig.suptitle('old-new sim comparison: Directly measured data')
    for icol in np.arange(0,len(EVT_SUB_COLS)):
        ax=fig.add_subplot(len(EVT_SUB_COLS),1,icol+1)
        ax.set_title(EVT_SUB_COLS[icol])
        ax.set_xlabel(EVT_SUB_COLS[icol])
        for c,df in zip(colors,dfs):
            plt.hist(df['el_%s'%EVT_SUB_COLS[icol]],100,(XMIN[icol],XMAX[icol]),
                            histtype='step',color=c,label='old')
            # plt.hist(df['el_%s'%EVT_SUB_COLS[icol]],100,(XMIN[icol],XMAX[icol]),
            #                 histtype='step',color='red',label='new')

    plt.legend()

    # fig2=plt.figure(figsize=(8,8))
    # fig2.suptitle('2D')
    # ax_old=fig2.add_subplot(211)
    # h2,x,y = np.histogram2d(pnl_h10['OLD'].el_p,pnl_h10['OLD'].el_cz,100,[[0,5],[0.5,1]])
    # plt.pcolormesh(x,y,h2)
    # ax_new=fig2.add_subplot(212)
    # h2,x,y = np.histogram2d(pnl_h10['NEW'].el_p,pnl_h10['NEW'].el_cz,100,[[0,5],[0.5,1]])
    # plt.pcolormesh(x,y,h2)

    
    