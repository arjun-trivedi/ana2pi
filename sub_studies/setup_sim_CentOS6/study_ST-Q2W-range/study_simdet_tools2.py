from __future__ import division

#rootpy
from root_numpy import root2array, root2rec, tree2rec
from root_numpy.testdata import get_filepath
from rootpy.plotting import Hist2D, Hist
from rootpy.plotting import root2matplotlib as rplt
from rootpy.interactive import wait
#PyROOT
import ROOT

#scipy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from math import *

import os

W_TH = 1.216

DATADIR=os.environ['STUDY_STQ2WRANGE_DATADIR']
ANADIR=os.environ['STUDY_STQ2WRANGE_ANADIR']
OUTDIR="" #OUTDIR = ANADIR/q2wdir

dT=""
dR=""
VARS=""
DRVD_COLS=['px_e','py_e','pz_e','px_p','py_p','pz_p','px_pip','py_pip','pz_pip','px_pim','py_pim','pz_pim']
HRANGE=""
FRANGE=""
SEL_TOPS=""

Q2MIN=""
Q2MAX=""
NQ2BINS=""
Q2BINW=""
Q2BINS_LE=""
Q2BINS_UE=""
WMIN=""
WMAX=""
NWBINS=""
WBINW=""
WBINS_LE=""
WBINS_UE=""

Q2BINS=""
WBINS=""
HOFT=""
HRES=""

def init(q2wdirs,nq2bins,nwbins,variables,hrange,frange,tops):
    """
    This function intializes the necessary information needed in studying the Simulated-Detector
    
    Input arguments:
    ----------------
    q2wdirs = list of q2wdirs relative to $STUDY_STQ2WRANGE_DATADIR,list of str
    nq2bins = Number of Q2 bins,int
    nwbins = Number of W bins,int
    variables = Variables for which offset and resolution is to be extracted,list of str
    hrange=for each variable, range for ST-SR histogram,list of [min,max]
    frange=for each variable, range over which ST-SR histogram is fitted,list of [min,max]
    tops = Topology to be used, list of int
    """
    #print init.__doc__
    
    #-- Import data 
    global dT,dR,VARS
    cols=['top']
    VARS=variables
    if any(items in VARS for items in DRVD_COLS):
        #print "Don't know what to do"
        #return
        drct_cols=[item for item in VARS if item not in DRVD_COLS]
        drvd_cols=[item for item in VARS if item in DRVD_COLS]
        if not 'theta_e' in VARS: drct_cols.extend(['theta_e'])#needed for drvd_cols
        if not 'phi_e' in VARS: drct_cols.extend(['phi_e'])#needed for drvd_cols    
        cols.extend(drct_cols)
        print "Branches to load into DataFrame=",cols
        #f = os.path.join(DATADIR,q2wdir,'recon/d2pi.root')
        f=[]
        for q2wdir in q2wdirs:
            f.append(os.path.join(DATADIR,q2wdir,'recon/d2pi.root'))
        arrT = root2array(f,'d2piTR/T/tT',branches=cols)
        arrR = root2array(f,'d2piTR/R/tR',branches=cols)
        dT = pd.DataFrame(arrT)
        dR = pd.DataFrame(arrR)
        #-- Add drvd_cols
        for d in (dT,dR):
            d['px_e']=d['p_e']*np.sin(np.deg2rad(d['theta_e']))*np.cos(np.deg2rad(d['phi_e']))
            d['py_e']=d['p_e']*np.sin(np.deg2rad(d['theta_e']))*np.sin(np.deg2rad(d['phi_e']))
            d['pz_e']=d['p_e']*np.cos(np.deg2rad(d['theta_e']))
    else:
        cols.extend(VARS)
        print "Branches to load into DataFrame=",cols
        #f = os.path.join(DATADIR,q2wdir,'recon/d2pi.root')
        f=[]
        for q2wdir in q2wdirs:
            f.append(os.path.join(DATADIR,q2wdir,'recon/d2pi.root'))
        arrT = root2array(f,'d2piTR/T/tT',branches=cols)
        arrR = root2array(f,'d2piTR/R/tR',branches=cols)
        dT = pd.DataFrame(arrT)
        dR = pd.DataFrame(arrR)
        
    #-- get HRANGE,FRANGE
    global HRANGE,FRANGE
    HRANGE=hrange
    FRANGE=frange

    #-- Determine Topologies to be used
    global SEL_TOPS
    SEL_TOPS=""
    itr=0
    for top in tops:
        if itr==0:
            SEL_TOPS+="(dR['top']==%d)"%top
        else:
            SEL_TOPS+="|(dR['top']==%d)"%top
        itr+=1
    print "SEL_TOPS = %s"%SEL_TOPS
    
    #-- Determine kinematic range of Simulated-Thrown (ST) Events
    global Q2MIN,Q2MAX,WMIN,WMAX
    Q2MIN=round(min(dT['Q2']),2)
    Q2MAX=round(max(dT['Q2']),2)
    WMIN=round(min(dT['W']),2)
    WMAX=round(max(dT['W']),2)
    print "Simulated-Thrown Q2,W = [%.2f,%.2f]GeV^2,[%.2f,%.2f]GeV"%(Q2MIN,Q2MAX,WMIN,WMAX)
    
    #-- Determine Q2,W binning
    global NQ2BINS,Q2BINW,Q2BINS_LE,Q2BINS_UE
    global NWBINS,WBINW,WBINS_LE,WBINS_UE
    NQ2BINS=nq2bins
    Q2BINW=round((Q2MAX-Q2MIN)/NQ2BINS,2)
    Q2BINS_LE=[Q2MIN+(i*Q2BINW) for i in range(NQ2BINS)]
    Q2BINS_UE=[Q2BINS_LE[i]+Q2BINW for i in range(NQ2BINS)]
    print "*** Q2 binning ***"
    print "NQ2BINS=%d,Q2BINW=%0.2f GeV^2"%(NQ2BINS,Q2BINW)
    print ["%0.2f" % i for i in Q2BINS_LE]
    print ["%0.2f" % i for i in Q2BINS_UE]
    NWBINS=nwbins
    WBINW=round((WMAX-WMIN)/NWBINS,2)
    WBINS_LE=[WMIN+(i*WBINW) for i in range(NWBINS)]
    WBINS_UE=[WBINS_LE[i]+WBINW for i in range(NWBINS)]
    print "*** W binning ***"
    print "NWBINS=%d,WBINW=%0.2f GeV"%(NWBINS,WBINW)
    print ["%0.2f" % i for i in WBINS_LE]
    print ["%0.2f" % i for i in WBINS_UE]
    
    #-- For the variables, create hres,hofst
    global Q2BINS,WBINS
    global HOFT,HRES
    HOFT=[]
    HRES=[]
    Q2BINS=Q2BINS_LE
    Q2BINS.extend([Q2BINS_UE[-1]])
    WBINS=WBINS_LE
    WBINS.extend([WBINS_UE[-1]])
    for var in VARS:
        HOFT.append(Hist2D(WBINS,Q2BINS,name="%s_OFT"%var,title="%s Offset"%var))
        HRES.append(Hist2D(WBINS,Q2BINS,name="%s_RES"%var,title="%s Resolution"%var))
#         HOFT.append(Hist2D(NWBINS,WMIN,WMAX,NQ2BINS,Q2MIN,Q2MAX,name="%s_OFT"%var,title="%s Offset"%var))
#         HRES.append(Hist2D(NWBINS,WMIN,WMAX,NQ2BINS,Q2MIN,Q2MAX,name="%s_RES"%var,title="%s Resolution"%var))
    print 'HOFT=',HOFT
    print 'HRES=',HRES
        
    # -- Create OUTDIR
    global OUTDIR
    #OUTDIR=os.path.join(ANADIR,q2wdir)
    OUTDIR=os.path.join(ANADIR,"_".join(q2wdirs))

    if not os.path.exists(OUTDIR):
        os.makedirs(OUTDIR)
    print "OUTDIR=%s"%OUTDIR