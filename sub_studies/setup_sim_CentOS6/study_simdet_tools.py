from __future__ import division
import ROOT
from root_numpy import root2array, root2rec, tree2rec
from root_numpy.testdata import get_filepath
from rootpy.plotting import Hist2D, Hist
from rootpy.interactive import wait
import ROOT

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from math import *

import os


dfT=""
dfR=""
datadir=''
be=''

NPARTS=4
PARTS=['e','p','pip','pim']

NCOLS=13+12# 3- components of momentum for each of the 4 particles
# COLS=['p_e',  'p_p',  'p_pip','p_pim', 'Q2',    'W',   'M_ppip', 'M_ppim', 'M_pippim', 'mm2ppippim', 'mmppip','mmppim','mmpippim']
# XMIN=[-0.20, -0.08,   -0.08,  -0.08,   -0.06,  -0.06,  -0.06,    -0.06,    -0.06,       -0.001,       -0.2,-0.2,-0.2]
# XMAX=[ 0.20,  0.08,    0.08,   0.08,    0.06,   0.06,   0.06,     0.06,     0.06,        0.001,        0.2,0.2,0.2]
COLS=['p_e','p_p','p_pip','p_pim', 'px_e','py_e','pz_e', 'px_p','py_p','pz_p', 'px_pip','py_pip','pz_pip', 'px_pim','py_pim','pz_pim','Q2','W',    'M_ppip','M_ppim','M_pippim', 'mm2ppippim','mmppip','mmppim','mmpippim']
XMIN=[-0.20,-0.08,-0.08,-0.08,     -0.08,-0.08,-0.08,    -0.08,-0.08,-0.08,    -0.08,-0.08,-0.08,          -0.08,-0.08,-0.08,         -0.06,-0.06, -0.06,-0.06,-0.06,           -0.001,-0.2,-0.2,-0.2]
XMAX=[ 0.20, 0.08, 0.08, 0.08,      0.08, 0.08, 0.08,     0.08, 0.08, 0.08,     0.08, 0.08, 0.08,           0.08, 0.08, 0.08,          0.06, 0.06,  0.06, 0.06, 0.06,            0.001, 0.2, 0.2, 0.2]
NBINS=100
dP_E,dP_P,dP_PIP,dP_PIM, dPX_E,dPY_E,dPZ_E, dPX_P,dPY_P,dPZ_P, dPX_PIP,dPY_PIP,dPZ_PIP, dPX_PIM,dPY_PIM,dPZ_PIM, dQ2,dW,  dM_PPIP,dM_PPIM,dM_PIPPIM, dMM2PPIPPIM,dMMPPIP,dMMPPIM,dMMPIPPIM = range(NCOLS)
# means=[[] for i in range(NCOLS)]
# sgmas=[[] for i in range(NCOLS)]
    
NTOPS=4
TOPS=[1,2,3,4]

def load_data(beam_energy):
	"""
	Load data
	"""
	print load_data.__doc__

	global dfT
	global dfR
	global datadir
	global be

	be=beam_energy
	datadir=os.environ['SETUPSIMCENTOS6_DATADIR']
	f = os.path.join(datadir,'d2pi_be%d.root'%be)
	arrT = root2array(f,'d2piTR/T/tT')#,stop=10000)#,start=1,stop=5)#,start=1,stop=5)#start=1,stop=2)
	arrR = root2array(f,'d2piTR/R/tR')#,stop=10000)#,start=1,stop=5)#,start=1,stop=5)#start=1,stop=2)

	dfT = pd.DataFrame(arrT)
	dfR = pd.DataFrame(arrR)

	#Add some columns to resolve momentum into its x,y,z components
	for df in (dfT,dfR):
		df['px_e']=df['p_e']*np.sin(np.deg2rad(df['theta_e']))*np.cos(np.deg2rad(df['phi_e']))
		df['py_e']=df['p_e']*np.sin(np.deg2rad(df['theta_e']))*np.sin(np.deg2rad(df['phi_e']))
		df['pz_e']=df['p_e']*np.cos(np.deg2rad(df['theta_e']))
		df['px_p']=df['p_p']*np.sin(np.deg2rad(df['theta_p']))*np.cos(np.deg2rad(df['phi_p']))
		df['py_p']=df['p_p']*np.sin(np.deg2rad(df['theta_p']))*np.sin(np.deg2rad(df['phi_p']))
		df['pz_p']=df['p_p']*np.cos(np.deg2rad(df['theta_p']))
		df['px_pip']=df['p_pip']*np.sin(np.deg2rad(df['theta_pip']))*np.cos(np.deg2rad(df['phi_pip']))
		df['py_pip']=df['p_pip']*np.sin(np.deg2rad(df['theta_pip']))*np.sin(np.deg2rad(df['phi_pip']))
		df['pz_pip']=df['p_pip']*np.cos(np.deg2rad(df['theta_pip']))
		df['px_pim']=df['p_pim']*np.sin(np.deg2rad(df['theta_pim']))*np.cos(np.deg2rad(df['phi_pim']))
		df['py_pim']=df['p_pim']*np.sin(np.deg2rad(df['theta_pim']))*np.sin(np.deg2rad(df['phi_pim']))
		df['pz_pim']=df['p_pim']*np.cos(np.deg2rad(df['theta_pim']))

def plot_q2w():
	"""
	Look at Q2 Vs. W
	"""
	print plot_q2w.__doc__

	q2T=dfT['Q2']
	wT=dfT['W']
	q2wT=np.array((wT,q2T))
	#print q2T
	#print wT
	#print q2wT
	#print q2wT.transpose()

	q2R=dfR['Q2']
	wR=dfR['W']
	q2wR=np.array((wR,q2R))

	hq2wT=Hist2D(100,1.2,2.0,100,1.8,2.6)	
	hq2wR=Hist2D(100,1.2,2.0,100,1.8,2.6)
	hq2wT.fill_array(q2wT.transpose())
	hq2wR.fill_array(q2wR.transpose())
	zT = np.array(hq2wT.z()).T
	zR = np.array(hq2wR.z()).T
	fig, (ax1,ax2) = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))
	ax1.set_title('q2wT')
	im1 = ax1.imshow(zT,
    		extent=[hq2wT.xedges(0), hq2wT.xedges(-1),
            hq2wT.yedges(0), hq2wT.yedges(-1)],
    		interpolation='nearest',
    		aspect='auto',
    		origin='lower')
	ax2.set_title(q2wR)
	im2 = ax2.imshow(zR,
    		extent=[hq2wR.xedges(0), hq2wR.xedges(-1),
            hq2wR.yedges(0), hq2wR.yedges(-1)],
    		interpolation='nearest',
    		aspect='auto',
    		origin='lower')

def plot_kinematics():
    """
    Plot Reconstructed kinematics of e,p,pip,pim at the e vertex 
    Compare with Thrown kinematics
    """
    print plot_kinematics.__doc__
    NPART=4
    COLS=['p_e','p_p','p_pip','p_pim']
    XLOWS =[2.5,0.0,0.0,0.0]
    XHIGHS=[4.1,3.0,3.0,3.0]
    alpha=0.5
    figT = plt.figure(figsize=(20,8))
    for i in range(NPART):
        plt.subplot(1,NPART,i+1)
        ret=plt.hist(dfT[COLS[i]],100,range=(XLOWS[i],XHIGHS[i]),alpha=alpha,label='ST',color=['blue'])#,color=[red,blue])
        ret=plt.hist(dfR[COLS[i]],100,range=(XLOWS[i],XHIGHS[i]),alpha=alpha,label='SR',color=['red'])
        plt.legend(loc='best')
    plt.show()

def get_detector_resolution_and_offset():
    """
    Plot (Simulated)detector "Offset" & "Resolution"
    Offset:     Mean of ST-SR distribution
    Resolution: RMS of ST-SR distribution
    """
    print get_detector_resolution_and_offset.__doc__
        
    ROOT.gStyle.SetOptFit(1111)
    ROOT.gStyle.SetStatW(0.4); 
    ROOT.gStyle.SetStatH(0.2); 
    CWIDTH=1200
    CHEIGHT=1000
        
    deltas=[[] for i in range(NCOLS)]
    hs=[[] for i in range(NCOLS)]
    means=[[] for i in range(NCOLS)]
    sgmas=[[] for i in range(NCOLS)]
    cs=[]
    outdir='%s/be%d'%(datadir,be)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    for icol in range(NCOLS):
        #print "processing %s"%COLS[ipart]
        cs.append(ROOT.TCanvas('Delta_%s'%COLS[icol],'Delta_%s'%COLS[icol],CWIDTH,CHEIGHT))
        cs[icol].Divide(2,2)
        for itop in range(NTOPS):
            #print "processing top %d"%TOPS[itop]
         	sel=(dfR['top']==TOPS[itop])
        	deltas[icol].append(dfR[COLS[icol]][sel]-dfT[COLS[icol]][sel])
        	pad=cs[icol].cd(itop+1)
        	pad.SetGrid()
         	hs[icol].append(Hist(NBINS,XMIN[icol],XMAX[icol]))
        	hs[icol][itop].fill_array(deltas[icol][itop])
         	hs[icol][itop].Draw()
        	if COLS[icol]=='mmppip' or COLS[icol]=='mmppim':
        		hs[icol][itop].Fit("gaus","","",-0.03,0.03)
        	elif COLS[icol]=='mmpippim':
        		hs[icol][itop].Fit("gaus","","",-0.07,0.03)
        	else:
        		hs[icol][itop].Fit("gaus","","",-0.02,0.02)        
        	pad.Update()
        	f=hs[icol][itop].GetFunction("gaus")
        	if f!=None:
        		means[icol].append(f.GetParameter(1))
        		sgmas[icol].append(f.GetParameter(2))
        	else:
        		means[icol].append(0)
        		sgmas[icol].append(0)
        cs[icol].SaveAs("%s/%s.png"%(outdir,cs[icol].GetName()))
  #   if not ROOT.gROOT.IsBatch():
		# plt.show()
	 # 	# wait for you to close the ROOT canvas before exiting
	 # 	wait(True)
        
    return means,sgmas      
    #plt.show()

def plot_detector_resolution_and_offset(means,sgmas):
	fig_nrows=8
	fig_ncols=2
	fig,(axs)=plt.subplots(nrows=fig_nrows,ncols=fig_ncols,figsize=(20,40))
	subplot_datas=[(dP_E,dP_P,dP_PIP,dP_PIM),
    			   (dPX_E,dPY_E,dPZ_E),
    			   (dPX_P,dPY_P,dPZ_P),
    			   (dPX_PIP,dPY_PIP,dPZ_PIP),
    			   (dPX_PIM,dPY_PIM,dPZ_PIM),
    			   (dQ2,dW),
    			   (dM_PPIP,dM_PPIM,dM_PIPPIM),
    			   (dMM2PPIPPIM,dMMPPIP,dMMPPIM,dMMPIPPIM)]
	marker_styles=['b^','g>','rv','k<']
    
	for irow in range(fig_nrows):
		ax=axs[irow][0]
		lns=[]
		iln=0
		if irow!=7:
			for idata in subplot_datas[irow]:
				lns.append(ax.plot(TOPS,means[idata],
					       marker_styles[iln],markersize=10,label=COLS[idata]))
				iln+=1
		else:
			ax.plot(TOPS,[means[dMM2PPIPPIM][0],means[dMMPPIP][1],means[dMMPPIM][2],means[dMMPIPPIM][3]],
					       marker_styles[iln],markersize=10)
		ax.set_xlim(0,5)
		ax.set_ylim(ax.get_ylim()[0]-0.001,ax.get_ylim()[1]+0.001)
		ax.set_xlabel('Top',size='xx-large')
		ax.set_ylabel('Offset(SR-ST)',size='xx-large')
		ax.legend(loc='best',prop={'size':10})
		ax.grid(True)
		
		ax=axs[irow][1]
		lns=[]
		iln=0
		if irow!=7:
			for idata in subplot_datas[irow]:
				lns.append(ax.plot(TOPS,sgmas[idata],
					       marker_styles[iln],markersize=10,label=COLS[idata]))
				iln+=1
		else:
			ax.plot(TOPS,[sgmas[dMM2PPIPPIM][0],sgmas[dMMPPIP][1],sgmas[dMMPPIM][2],sgmas[dMMPIPPIM][3]],
					       marker_styles[iln],markersize=10)
		ax.set_xlim(0,5)
		ax.set_ylim(ax.get_ylim()[0]-0.001,ax.get_ylim()[1]+0.001)
		ax.set_xlabel('Top',size='xx-large')
		ax.set_ylabel('Resolution',size='xx-large')
		ax.legend(loc='best',prop={'size':10})
		ax.grid(True)	
   

# def plot_detector_resolution_and_offset(means,sgmas):
#     fig,(axs)=plt.subplots(nrows=8,ncols=2,figsize=(20,15))
    
#     ##-- momenta
#     ax=axs[0][0]
#     lns=ax.plot(TOPS,means[0],b^,
#     			TOPS,means[1],g>,
#                 TOPS,means[2],rv,
#                 TOPS,means[3],k<)
#     ax.set_xlim(0,5)
#     ax.set_ylim(ax.get_ylim()[0]-0.001,ax.get_ylim()[1]+0.001)
#     ax.set_xlabel(top)
#     ax.set_ylabel(offset:SR-ST)
#     ax.legend(lns,[COLS[0],COLS[1],COLS[2],COLS[3]],loc=best,prop={size:9})
#     for i in range(len(lns)):
#     	lns[i].set_markersize(10)
        
#     ax=axs[0][1]
#     lns=ax.plot(TOPS,sgmas[0],b^,
#                 TOPS,sgmas[1],g>,
#                 TOPS,sgmas[2],rv,
#                 TOPS,sgmas[3],k<)
#     ax.set_xlim(0,5)
#     ax.set_ylim(ax.get_ylim()[0]-0.001,ax.get_ylim()[1]+0.001)
#     ax.set_xlabel(top)
#     ax.set_ylabel(resolution)
#     ax.legend(lns,[COLS[0],COLS[1],COLS[2],COLS[3]],loc=best,prop={size:9})
#     for i in range(len(lns)):
#     	lns[i].set_markersize(10)

#     ##-- Q2,W
#     ax=axs[1][0]
#     lns=ax.plot(TOPS,means[4],b^,
#                 TOPS,means[5],g>,)                  
#     ax.set_xlim(0,5)
#     ax.set_ylim(ax.get_ylim()[0]-0.001,ax.get_ylim()[1]+0.001)
#     #ax.set_ylim(-0.003,0.003)
#     ax.set_xlabel(top)
#     ax.set_ylabel(offset)
#     ax.grid(True)
#     ax.legend(lns,[COLS[4],COLS[5]],loc=best,prop={size:9})
#     for i in range(len(lns)):
#     	lns[i].set_markersize(10)
    
#     ax=axs[1][1]
#     lns=ax.plot(TOPS,sgmas[4],b^,
#                 TOPS,sgmas[5],g>,)                  
#     ax.set_xlim(0,5)
#     ax.set_ylim(ax.get_ylim()[0]-0.001,ax.get_ylim()[1]+0.001)
#     ax.set_xlabel(top)
#     ax.set_ylabel(resolution)
#     ax.legend(lns,[COLS[4],COLS[5]],loc=best,prop={size:9})
#     for i in range(len(lns)):
#     	lns[i].set_markersize(10)
    
#     #-- {Mij}
#     ax=axs[2][0]
#     lns=ax.plot(TOPS,means[6],b^,
#                 TOPS,means[7],g>,
#                 TOPS,means[8],rv)
#     ax.set_xlim(0,5)
#     ax.set_ylim(ax.get_ylim()[0]-0.001,ax.get_ylim()[1]+0.001)
#     ax.set_xlabel(top)
#     ax.set_ylabel(offset:SR-ST)
#     ax.legend(lns,[COLS[6],COLS[7],COLS[8]],loc=best,prop={size:9})
#     for i in range(len(lns)):
#     	lns[i].set_markersize(10)
    
    
#     ax=axs[2][1]
#     lns=ax.plot(TOPS,sgmas[6],b^,
#                 TOPS,sgmas[7],g>,
#                 TOPS,sgmas[8],rv)
#     ax.set_xlim(0,5)
#     ax.set_ylim(ax.get_ylim()[0]-0.001,ax.get_ylim()[1]+0.001)
#     ax.set_xlabel(top)
#     ax.set_ylabel(resolution)
#     ax.legend(lns,[COLS[6],COLS[7],COLS[8]],loc=best,prop={size:9})
#     for i in range(len(lns)):
#     	lns[i].set_markersize(10)
    
#     #-- MMs
#     ax=axs[3][0]
#     lns=ax.plot(TOPS,[means[9][0],means[10][1],means[11][2],means[12][3]],b^)
#     ax.set_xlim(0,5)
#     #ax.set_ylim(-0.001,0.006)
#     ax.set_xlabel(top)
#     ax.set_ylabel(offset:SR-ST)
#     ax.grid(True)
#     #ax.legend(lns,[COLS[6],COLS[7],COLS[8]],loc=best,prop={size:9})
#     for i in range(len(lns)):
#     	lns[i].set_markersize(10)
    
#     ax=axs[3][1]
#     lns=ax.plot(TOPS,[sgmas[9][0],sgmas[10][1],sgmas[11][2],sgmas[12][3]],b^)
#     ax.set_xlim(0,5)
#     ax.set_ylim(-0.001,0.030)
#     ax.set_xlabel(top)
#     ax.set_ylabel(resolution)
#     ax.grid(True)
#     for i in range(len(lns)):
#     	lns[i].set_markersize(10)
