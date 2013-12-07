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

def at_hw():
    print 'hello world'

NENTRIES=1000
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
    for i in np.arange(0,len(df[col])):
        if len(df[col][i])==0:        
            d.append(-9999)
        else:
            d.append(df[col][i][0])
    df['el_%s'%(col)]=d

def add_reco_elcols(d):
    qsq = []
    w = []
    for i in np.arange(0,len(d.el_p)):
        p = d.el_p[i]
        px = p * d.el_cx[i]
        py = p * d.el_cy[i]
        pz = p * d.el_cz[i]
        lve = ROOT.TLorentzVector(0,0,0,0)
        e = math.sqrt(p*p+MASS_E*MASS_E) 
        lve.SetPxPyPzE(px,py,pz,e)
        lvq = lvE0-lve
        lvw = lvq+lvP0
        qsq.append(-lvq.Mag2())
        w.append(lvw.Mag())
            
    d['qsq']=qsq 
    d['w']=w

def comp_qsq_w():
    fig = plt.figure(figsize=(5,5))
    fig.suptitle('old-new sim comparison: Reconstructed Q2,W')
    ax_w=fig.add_subplot(211)
    ax_w.set_title('W')
    ax_w.set_xlabel('W(GeV)')
    n,bins,patches=plt.hist(pnl_h10['OLD'].w,100,(1.0,2.0),
                            histtype='step',color='black',label='old')
    n,bins,patches=plt.hist(pnl_h10['NEW'].w,100,(1.0,2.0),
                            histtype='step',color='red',label='new')
    ax_qsq=fig.add_subplot(212)
    ax_qsq.set_title('Q2')
    ax_qsq.set_xlabel('Q2(GeV)')
    n,bins,patches=plt.hist(pnl_h10['OLD'].qsq,100,(1.0,2.0),
                            histtype='step',color='black',label='old')
    n,bins,patches=plt.hist(pnl_h10['NEW'].qsq,100,(1.0,2.0),
                            histtype='step',color='red',label='new')
    plt.legend()

def comp_basic():
    fig=plt.figure()
    fig.suptitle('old-new sim comparison: Directly measured data')
    ax_etot=fig.add_subplot(111)
    ax_etot.set_title('etot')
    ax_etot.set_xlabel('etot(GeV)')
    n,bins,patches=plt.hist(pnl_h10['OLD'].el_etot,100,(0.0,1.0),
                            histtype='step',color='black',label='old')
    n,bins,patches=plt.hist(pnl_h10['NEW'].el_etot,100,(0.0,1.0),
                            histtype='step',color='red',label='new')
    
    plt.legend()