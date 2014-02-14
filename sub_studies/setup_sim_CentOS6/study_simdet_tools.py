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

import os

dfT=""
dfR=""
datadir=''
xx=None

NPARTS=4
PARTS=['e','p','pip','pim']

NCOLS=13
COLS=['p_e',  'p_p',  'p_pip','p_pim', 'Q2',    'W',   'M_ppip', 'M_ppim', 'M_pippim', 'mm2ppippim', 'mmppip','mmppim','mmpippim']
XMIN=[-0.20, -0.08,   -0.08,  -0.08,   -0.06,  -0.06,  -0.06,    -0.06,    -0.06,       -0.001,       -0.2,-0.2,-0.2]
XMAX=[ 0.20,  0.08,    0.08,   0.08,    0.06,   0.06,   0.06,     0.06,     0.06,        0.001,        0.2,0.2,0.2]
NBINS=100
means=[[] for i in range(NCOLS)]
sgmas=[[] for i in range(NCOLS)]
    
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
	ax2.set_title('q2wR')
	im2 = ax2.imshow(zR,
    		extent=[hq2wR.xedges(0), hq2wR.xedges(-1),
            hq2wR.yedges(0), hq2wR.yedges(-1)],
    		interpolation='nearest',
    		aspect='auto',
    		origin='lower')

def plot_kinematics():
    """
    Plot Reconstructed kinematics of e',p',pip,pim at the e' vertex 
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
        ret=plt.hist(dfT[COLS[i]],100,range=(XLOWS[i],XHIGHS[i]),alpha=alpha,label='T',color=['blue'])#,color=['red','blue'])
        ret=plt.hist(dfR[COLS[i]],100,range=(XLOWS[i],XHIGHS[i]),alpha=alpha,label='R',color=['red'])
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
    global means#=[[] for i in range(NCOLS)]
    global sgmas#=[[] for i in range(NCOLS)]
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
            hs[icol].append(Hist(NBINS,XMIN[icol],XMAX[icol]))
            hs[icol][itop].fill_array(deltas[icol][itop])
            hs[icol][itop].Draw()
            if COLS[icol]=='mmppip' or COLS[icol]=='mmppim':
                hs[icol][itop].Fit("gaus","","",-0.03,0.03)
            elif COLS[icol]=='mmpippim':
                hs[icol][itop].Fit("gaus","","",-0.07,0.03)
            else:
                hs[icol][itop].Fit("gaus")        
            pad.Update()
            f=hs[icol][itop].GetFunction("gaus")
            if f!=None:
                means[icol].append(f.GetParameter(1))
                sgmas[icol].append(f.GetParameter(2))
            else:
                means[icol].append(0)
                sgmas[icol].append(0)
        cs[icol].SaveAs('%s/%s.png'%(outdir,cs[icol].GetName()))
  #   if not ROOT.gROOT.IsBatch():
		# plt.show()
	 # 	# wait for you to close the ROOT canvas before exiting
	 # 	wait(True)
        
    #return means,sgmas      
    #plt.show()

def plot_detector_resolution_and_offset():
    fig,(ax)=plt.subplots(nrows=4,ncols=2,figsize=(15,10))
    
    ##-- momenta
    ax1=ax[0][0]
    lns1=ax1.plot(TOPS,means[0],'b^',
                  TOPS,means[1],'g>',
                  TOPS,means[2],'rv',
                  TOPS,means[3],'k<')
    ax1.set_xlim(0,5)
    ax1.set_ylim(ax1.get_ylim()[0]-0.001,ax1.get_ylim()[1]+0.001)
    ax1.set_xlabel('top')
    ax1.set_ylabel('offset:SR-ST')
    ax1.legend(lns1,[COLS[0],COLS[1],COLS[2],COLS[3]],loc='best',prop={'size':9})
    #ax1.legend(loc='best')
    
    ax2=ax[0][1]
    lns2=ax2.plot(TOPS,sgmas[0],'b^',
                  TOPS,sgmas[1],'g>',
                  TOPS,sgmas[2],'rv',
                  TOPS,sgmas[3],'k<')
    ax2.set_xlim(0,5)
    ax2.set_ylim(ax2.get_ylim()[0]-0.001,ax2.get_ylim()[1]+0.001)
    ax2.set_xlabel('top')
    ax2.set_ylabel('resolution')
    ax2.legend(lns2,[COLS[0],COLS[1],COLS[2],COLS[3]],loc='best',prop={'size':9})
    
    ##-- Q2,W
    ax3=ax[1][0]
    lns3=ax3.plot(TOPS,means[4],'b^',
                  TOPS,means[5],'g>',)                  
    ax3.set_xlim(0,5)
    ax3.set_ylim(ax3.get_ylim()[0]-0.001,ax3.get_ylim()[1]+0.001)
    #ax3.set_ylim(-0.003,0.003)
    ax3.set_xlabel('top')
    ax3.set_ylabel('offset')
    ax3.grid(True)
    ax3.legend(lns2,[COLS[4],COLS[5]],loc='best',prop={'size':9})
    
    ax4=ax[1][1]
    lns4=ax4.plot(TOPS,sgmas[4],'b^',
                  TOPS,sgmas[5],'g>',)                  
    ax4.set_xlim(0,5)
    ax4.set_ylim(ax4.get_ylim()[0]-0.001,ax4.get_ylim()[1]+0.001)
    ax4.set_xlabel('top')
    ax4.set_ylabel('resolution')
    ax4.legend(lns4,[COLS[4],COLS[5]],loc='best',prop={'size':9})
    
    #-- {Mij}
    ax5=ax[2][0]
    lns5=ax5.plot(TOPS,means[6],'b^',
                  TOPS,means[7],'g>',
                  TOPS,means[8],'rv')
    ax5.set_xlim(0,5)
    ax5.set_ylim(ax5.get_ylim()[0]-0.001,ax5.get_ylim()[1]+0.001)
    ax5.set_xlabel('top')
    ax5.set_ylabel('offset:SR-ST')
    ax5.legend(lns5,[COLS[6],COLS[7],COLS[8]],loc='best',prop={'size':9})
    
    
    ax6=ax[2][1]
    lns6=ax6.plot(TOPS,sgmas[6],'b^',
                  TOPS,sgmas[7],'g>',
                  TOPS,sgmas[8],'rv')
    ax6.set_xlim(0,5)
    ax6.set_ylim(ax6.get_ylim()[0]-0.001,ax6.get_ylim()[1]+0.001)
    ax6.set_xlabel('top')
    ax6.set_ylabel('resolution')
    ax6.legend(lns6,[COLS[6],COLS[7],COLS[8]],loc='best',prop={'size':9})
    
    #-- MMs
    ax7=ax[3][0]
    lns7=ax7.plot(TOPS,[means[9][0],means[10][1],means[11][2],means[12][3]],'b^')
    ax7.set_xlim(0,5)
    #ax7.set_ylim(-0.001,0.006)
    ax7.set_xlabel('top')
    ax7.set_ylabel('offset:SR-ST')
    ax7.grid(True)
    #ax7.legend(lns7,[COLS[6],COLS[7],COLS[8]],loc='best',prop={'size':9})
    
    ax8=ax[3][1]
    lns8=ax8.plot(TOPS,[sgmas[9][0],sgmas[10][1],sgmas[11][2],sgmas[12][3]],'b^')
    ax8.set_xlim(0,5)
    ax8.set_ylim(-0.001,0.030)
    ax8.set_xlabel('top')
    ax8.set_ylabel('resolution')
    ax8.grid(True)


# def plot_detector_resolution_and_offset():
#     fig,(ax)=plt.subplots(nrows=4,ncols=2,figsize=(15,10))
    
#     ##-- momenta
#     ax1=ax[0][0]
#     lns1=ax1.plot(TOPS,means[0],'b^',
#                   TOPS,means[1],'g>',
#                   TOPS,means[2],'rv',
#                   TOPS,means[3],'k<')
#     ax1.set_xlim(0,5)
#     ax1.set_ylim(ax1.get_ylim()[0]-0.001,ax1.get_ylim()[1]+0.001)
#     ax1.set_xlabel('top')
#     ax1.set_ylabel('offset:SR-ST')
#     ax1.legend(lns1,[COLS[0],COLS[1],COLS[2],COLS[3]],loc='best',prop={'size':9})
#     #ax1.legend(loc='best')
    
#     ax2=ax[0][1]
#     lns2=ax2.plot(TOPS,sgmas[0],'b^',
#                   TOPS,sgmas[1],'g>',
#                   TOPS,sgmas[2],'rv',
#                   TOPS,sgmas[3],'k<')
#     ax2.set_xlim(0,5)
#     ax2.set_ylim(ax2.get_ylim()[0]-0.001,ax2.get_ylim()[1]+0.001)
#     ax2.set_xlabel('top')
#     ax2.set_ylabel('resolution')
#     ax2.legend(lns2,[COLS[0],COLS[1],COLS[2],COLS[3]],loc='best',prop={'size':9})
    
#     ##-- Q2,W
#     ax3=ax[1][0]
#     lns3=ax3.plot(TOPS,means[4],'b^',
#                   TOPS,means[5],'g>',)                  
#     ax3.set_xlim(0,5)
#     ax3.set_ylim(ax3.get_ylim()[0]-0.001,ax3.get_ylim()[1]+0.001)
#     #ax3.set_ylim(-0.003,0.003)
#     ax3.set_xlabel('top')
#     ax3.set_ylabel('offset')
#     ax3.grid(True)
#     ax3.legend(lns2,[COLS[4],COLS[5]],loc='best',prop={'size':9})
    
#     ax4=ax[1][1]
#     lns4=ax4.plot(TOPS,sgmas[4],'b^',
#                   TOPS,sgmas[5],'g>',)                  
#     ax4.set_xlim(0,5)
#     ax4.set_ylim(ax4.get_ylim()[0]-0.001,ax4.get_ylim()[1]+0.001)
#     ax4.set_xlabel('top')
#     ax4.set_ylabel('resolution')
#     ax4.legend(lns4,[COLS[4],COLS[5]],loc='best',prop={'size':9})
    
#     #-- {Mij}
#     ax5=ax[2][0]
#     lns5=ax5.plot(TOPS,means[6],'b^',
#                   TOPS,means[7],'g>',
#                   TOPS,means[8],'rv')
#     ax5.set_xlim(0,5)
#     ax5.set_ylim(ax5.get_ylim()[0]-0.001,ax5.get_ylim()[1]+0.001)
#     ax5.set_xlabel('top')
#     ax5.set_ylabel('offset:SR-ST')
#     ax5.legend(lns5,[COLS[6],COLS[7],COLS[8]],loc='best',prop={'size':9})
    
    
#     ax6=ax[2][1]
#     lns6=ax6.plot(TOPS,sgmas[6],'b^',
#                   TOPS,sgmas[7],'g>',
#                   TOPS,sgmas[8],'rv')
#     ax6.set_xlim(0,5)
#     ax6.set_ylim(ax6.get_ylim()[0]-0.001,ax6.get_ylim()[1]+0.001)
#     ax6.set_xlabel('top')
#     ax6.set_ylabel('resolution')
#     ax6.legend(lns6,[COLS[6],COLS[7],COLS[8]],loc='best',prop={'size':9})
    
#     #-- MMs
#     ax7=ax[3][0]
#     lns7=ax7.plot(TOPS,[means[9][0],means[10][1],means[11][2],means[12][3]],'b^')
#     ax7.set_xlim(0,5)
#     #ax7.set_ylim(-0.001,0.006)
#     ax7.set_xlabel('top')
#     ax7.set_ylabel('offset:SR-ST')
#     ax7.grid(True)
#     #ax7.legend(lns7,[COLS[6],COLS[7],COLS[8]],loc='best',prop={'size':9})
    
#     ax8=ax[3][1]
#     lns8=ax8.plot(TOPS,[sgmas[9][0],sgmas[10][1],sgmas[11][2],sgmas[12][3]],'b^')
#     ax8.set_xlim(0,5)
#     ax8.set_ylim(-0.001,0.030)
#     ax8.set_xlabel('top')
#     ax8.set_ylabel('resolution')
#     ax8.grid(True)
    

# 	