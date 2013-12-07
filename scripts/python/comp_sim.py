#!/usr/bin/python
import ROOT
from root_numpy import root2array, root2rec, tree2rec
from root_numpy.testdata import get_filepath

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.pyplot as plt

import math

OLD,NEW     = range(0,2)
MCTK,COOKED = range(0,2)
sar_h10 = {} 
pnl_h10 = '' 

MASS_E = 0.000511
MASS_P = 0.93827203
E1F_P = 5.499; 
lvE0 = ROOT.TLorentzVector(0,0,E1F_P,math.sqrt(E1F_P*E1F_P+MASS_E*MASS_E));
lvP0 = ROOT.TLorentzVector(0,0,0,MASS_P);

ELCOLS=['p','cx','cy','cz','etot','ec_ei','ec_eo']

NENTRIES=2000
def import_data(fold,fnew,nentries=NENTRIES):
    """"return pnl_h10: items,rows,cols = OLD/NEW,h10-rows,h10-cols"""
    
    global pnl_h10

    #arr={}
    sar_h10[OLD] = root2array(fold,'h10',stop=NENTRIES)
    sar_h10[NEW] = root2array(fnew,'h10',stop=NENTRIES)

    df={}
    df[OLD] = pd.DataFrame(sar_h10[OLD])
    df[NEW] = pd.DataFrame(sar_h10[NEW])

    dpnl = {'OLD':df[OLD], 'NEW':df[NEW]}
    pnl_h10 = pd.Panel(dpnl)
    print "Number of rows:cols in pnl_h10['OLD']=",pnl_h10['OLD'].shape
    print "Number of rows:cols in pnl_h10['NEW']=",pnl_h10['NEW'].shape
    
def prep_data():
    for col in ELCOLS:
        add_basic_elcol(col,pnl_h10['OLD'])
        add_basic_elcol(col,pnl_h10['NEW'])
    
    add_reco_elcols(pnl_h10['OLD'])
    add_reco_elcols(pnl_h10['NEW'])

def add_basic_elcol(col,df):
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

def comp_basic():
    fig=plt.figure(figsize=(10,10))
    fig.suptitle('old-new sim comparison: Directly measured data')
    for icol in np.arange(0,len(ELCOLS)):
        ax=fig.add_subplot(len(ELCOLS),1,icol+1)
        ax.set_title(ELCOLS[icol])
        ax.set_xlabel(ELCOLS[icol])
        plt.hist(pnl_h10['OLD']['el_%s'%ELCOLS[icol]],100,(0.0,1.0),
                            histtype='step',color='black',label='old')
        plt.hist(pnl_h10['NEW']['el_%s'%ELCOLS[icol]],100,(0.0,1.0),
                            histtype='step',color='red',label='new')
    
    plt.legend()