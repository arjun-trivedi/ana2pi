#!/usr/bin/python
import ROOT
from root_numpy import root2array, root2rec, tree2rec
from root_numpy.testdata import get_filepath

import pandas as pd
import numpy as np

import math

OLD,NEW = range(0,2)
MCTK,COOKED = range(0,2)

def at_hw():
    print 'hello world'

def import_data():
    arr={}
    arr[OLD,COOKED] = root2array('/e1f.sim2pi.datadir/comp_old_new_sim/Q2W__1.9-2.5__1.3-1.9/cooked/1.root',
                 'h10')
    arr[NEW,COOKED] = root2array('/e1f.sim2pi.datadir/comp_old_new_sim/q2w2_ye_gpp/cooked/1.root',
                 'h10')

    d={}
    d[OLD,COOKED] = pd.DataFrame(arr[OLD,COOKED])
    d[NEW,COOKED] = pd.DataFrame(arr[NEW,COOKED])

    d[OLD,COOKED] = d[OLD,COOKED].ix[0:9999] #make size of dn=do

    return d[OLD,COOKED],d[NEW,COOKED]

    #print 'shape of old and new Tree (nrows,ncols) or (nevts,nbranches):'
    #print d[OLD,COOKED].shape
    #print d[NEW,COOKED].shape

MASS_E = 0.000511
MASS_P = 0.93827203
E1F_P = 5.499; 
lvE0 = ROOT.TLorentzVector(0,0,E1F_P,math.sqrt(E1F_P*E1F_P+MASS_E*MASS_E));
lvP0 = ROOT.TLorentzVector(0,0,0,MASS_P);

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

