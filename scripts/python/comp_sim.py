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
ph10 = ''

MASS_E = 0.000511
MASS_P = 0.93827203
E1F_P = 5.499; 
lvE0 = ROOT.TLorentzVector(0,0,E1F_P,math.sqrt(E1F_P*E1F_P+MASS_E*MASS_E));
lvP0 = ROOT.TLorentzVector(0,0,0,MASS_P);


def at_hw():
    print 'hello world'

NENTRIES=1000
def import_data(fold,fnew,nentries=NENTRIES):
    """"return ph10: items,rows,cols = OLD/NEW,h10-rows,h10-cols"""
    
    global ph10

    arr={}
    arr[OLD,COOKED] = root2array(fold,'h10',stop=NENTRIES)
    arr[NEW,COOKED] = root2array(fnew,'h10',stop=NENTRIES)

    d={}
    d[OLD,COOKED] = pd.DataFrame(arr[OLD,COOKED])
    d[NEW,COOKED] = pd.DataFrame(arr[NEW,COOKED])

    data = {'OLD':d[OLD,COOKED], 'NEW':d[NEW,COOKED]}
    ph10 = pd.Panel(data)
    print "Number of rows:cols in ph10['OLD']=",ph10['OLD'].shape
    print "Number of rows:cols in ph10['NEW']=",ph10['NEW'].shape
    
def prep_data():
    #addcol_qsq_w(ph10['OLD'])
    #addcol_qsq_w(ph10['NEW'])
    addcols(ph10['OLD'])
    addcols(ph10['NEW'])

def addcol_qsq_w(d):
    qsq = []
    w = []
    for i in np.arange(0,len(d.p)):
        if len(d.p[i])==0:
            qsq.append(-9999)
            w.append(-9999)
        else:
            p = d.p[i][0]
            px = p * d.cx[i][0]
            py = p * d.cy[i][0]
            pz = p * d.cz[i][0]
            lve = ROOT.TLorentzVector(0,0,0,0)
            e = math.sqrt(p*p+MASS_E*MASS_E) 
            lve.SetPxPyPzE(px,py,pz,e)
            lvq = lvE0-lve
            lvw = lvq+lvP0
            qsq.append(-lvq.Mag2())
            w.append(lvw.Mag())
            
    d['qsq']=qsq 
    d['w']=w

def addcols(d):
    #new cols to be added to df
    qsq = []
    w = []
    for i in np.arange(0,len(d.p)):
        if len(d.p[i])==0:
            qsq.append(-9999)
            w.append(-9999)
        else:
            p = d.p[i][0]
            px = p * d.cx[i][0]
            py = p * d.cy[i][0]
            pz = p * d.cz[i][0]
            lve = ROOT.TLorentzVector(0,0,0,0)
            e = math.sqrt(p*p+MASS_E*MASS_E) 
            lve.SetPxPyPzE(px,py,pz,e)
            lvq = lvE0-lve
            lvw = lvq+lvP0
            qsq.append(-lvq.Mag2())
            w.append(lvw.Mag())
    d['qsq']=qsq 
    d['w']=w

    el_etot=[]
    for i in np.arange(0,len(d.etot)):
        if len(d.etot[i])==0:        
            el_etot.append(-9999)
        else:
            el_etot.append(d.etot[i][0])
    d['el_etot']=el_etot


def comp_qsq_w():
    fig = plt.figure(figsize=(5,5))
    fig.suptitle('old-new sim comparison: Reconstructed Q2,W')
    ax_w=fig.add_subplot(211)
    ax_w.set_title('W')
    ax_w.set_xlabel('W(GeV)')
    n,bins,patches=plt.hist(ph10['OLD'].w,100,(1.0,2.0),
                            histtype='step',color='black',label='old')
    n,bins,patches=plt.hist(ph10['NEW'].w,100,(1.0,2.0),
                            histtype='step',color='red',label='new')
    ax_qsq=fig.add_subplot(212)
    ax_qsq.set_title('Q2')
    ax_qsq.set_xlabel('Q2(GeV)')
    n,bins,patches=plt.hist(ph10['OLD'].qsq,100,(1.0,2.0),
                            histtype='step',color='black',label='old')
    n,bins,patches=plt.hist(ph10['NEW'].qsq,100,(1.0,2.0),
                            histtype='step',color='red',label='new')
    plt.legend()

def comp_basic():
    fig=plt.figure()
    fig.suptitle('old-new sim comparison: Directly measured data')
    ax_etot=fig.add_subplot(111)
    ax_etot.set_title('etot')
    ax_etot.set_xlabel('etot(GeV)')
    n,bins,patches=plt.hist(ph10['OLD'].el_etot,100,(0.0,1.0),
                            histtype='step',color='black',label='old')
    n,bins,patches=plt.hist(ph10['NEW'].el_etot,100,(0.0,1.0),
                            histtype='step',color='red',label='new')
    plt.legend()