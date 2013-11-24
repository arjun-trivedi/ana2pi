#!/usr/bin/python

from __future__ import division
import os,sys,getopt, shutil, datetime, time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
pd.set_printoptions(max_columns=20)

MAX_TOPS,MIN_TOPS = range(2)
NTOPS=5
TOPS = [1,2,3,4,5];
TOP_LABELS  = ['Top1','Top2','Top3','Top4','Top5']
TOP_MARKERS = ['s','<','8' ,'p','>']
TOP_COLORS  =  ['r','k','r','r','k']

NSTATS=4
nFB_ST,nEB_SR,nFB_SR,nEB_SA = range(0,NSTATS)
STATS_NAME = ['nFB_ST','nEB_SR','nFB_SR','nEB_SA']

def plot_track_stats():
    for q2wbinnum in df_grpd_q2wbinnum.groups:
        df_q2wbinnum = df_grpd_q2wbinnum.get_group(q2wbinnum)
        q2wbinname = df_q2wbinnum['q2wbin'].tolist()[0]
        #print 'q2wbinname=',q2wbinname
        for stat in range(0,NSTATS):
            #print STATS_NAME[stat]
            fig=plt.figure('%d:%s'%(q2wbinnum,STATS_NAME[stat]))
            axs = []
            axs.append(fig.add_subplot(1,1,1,title='%s:%s'%(STATS_NAME[stat],q2wbinname)))
            axs.append(axs[MAX_TOPS].twinx())
            axs[MAX_TOPS].set_ylabel('%s t2,5'%STATS_NAME[stat])
            axs[MIN_TOPS].set_ylabel('%s t1,3,4'%STATS_NAME[stat],color='r')
            for tl in axs[MIN_TOPS].get_yticklabels():
                tl.set_color('r')
        
            lns = []
            for top in TOPS:
                itop=top-1
                if top==2 or top==5:
                    #print df['Sim'][df['Top']==top]
                    sel= (df_q2wbinnum['Top']==top)
                    lns.append(axs[MAX_TOPS].scatter(df_q2wbinnum['Sim'][sel],df_q2wbinnum[STATS_NAME[stat]][sel],color='k',
                                     marker=TOP_MARKERS[itop],label=TOP_LABELS[itop],s=50))
                else:
                   sel= (df_q2wbinnum['Top']==top)
                   lns.append(axs[MIN_TOPS].scatter(df_q2wbinnum['Sim'][sel],df_q2wbinnum[STATS_NAME[stat]][sel],color='r',
                                    marker=TOP_MARKERS[itop],label=TOP_LABELS[itop],s=50))
            labs = [l.get_label() for l in lns]
            axs[MAX_TOPS].legend(lns, labs, loc=2)
                    
            #save plots
            outdir = os.path.join(anadir,outdirname,q2wbinname)
            if not os.path.isdir(outdir):
                os.makedirs(outdir)
            #plt.figure(stat)
            fig.savefig('%s/%s.jpg'%(outdir,STATS_NAME[stat]))

def plot_latest_stats():
    for stat in range(0,NSTATS):
        fig=plt.figure('%s'%(STATS_NAME[stat]))
        ax = []
        ax.append(plt.subplot(1,1,1,title=STATS_NAME[stat],xticks=np.arange(1,nq2wbins+1,1),xlabel='q2wbinnum'))
        ax.append(ax[MAX_TOPS].twinx())
        ax[MAX_TOPS].set_ylabel('%s t2,5'%STATS_NAME[stat])
        ax[MIN_TOPS].set_ylabel(('%s t1,3,4'%STATS_NAME[stat]),color='r')
        for tl in ax[MIN_TOPS].get_yticklabels():
            tl.set_color('r')

        #Fill nbins[top][q2wbin]
        nbins = np.zeros((NTOPS,nq2wbins));
        q2wbs = np.arange(1,nq2wbins+1,1);
        for q2wbinnum in df_grpd_q2wbinnum.groups:
            #print 'q2wbinnum=',q2wbinnum
            iq2wbin = q2wbinnum-1
            df_q2wbinnum = df_grpd_q2wbinnum.get_group(q2wbinnum)
            siml = df_q2wbinnum['Sim'].max()
            for top in TOPS:
                itop = top-1
                sel = (df_q2wbinnum['Sim']==siml) & (df_q2wbinnum['Top']==(top))
                nbins[itop][iq2wbin] = df_q2wbinnum[STATS_NAME[stat]][sel]
    
        lns = []
        for top in TOPS:
            itop = top-1
            nbins_top=nbins[itop]
            if top==2 or top==5:
                lns.append(ax[MAX_TOPS].scatter(q2wbs, nbins_top, color='k', marker=TOP_MARKERS[itop], s=50,
                       label=TOP_LABELS[itop]))
            else:
                lns.append(ax[MIN_TOPS].scatter(q2wbs, nbins_top, color='r', marker=TOP_MARKERS[itop], s=50,
                       label=TOP_LABELS[itop]))    
        labs = [l.get_label() for l in lns]
        ax[MAX_TOPS].legend(lns, labs, loc=2)

        #savefig
        outdir = os.path.join(anadir,outdirname,'cum')
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        fig.savefig('%s/%s.jpg'%(outdir,STATS_NAME[stat]))
 

#MAIN part of the program

#get options
# global anadir
# global outdirname
# global csvfname
# global csvf
csvfname = 'simstats_vm.csv'
outdirname = 'simstats'
opts, args = getopt.getopt(sys.argv[1:],"h",["help", "e1fs1", "e1fs2", "hel="])

for opt, arg in opts:
    if opt in ('-h', '--help'):
        print 'plot_simstats.py --<e1fs1/e1fs2> --hel=<UNPOL/POS/NEG>'
        sys.exit()
    elif opt == "--e1fs1":
        anadir = os.environ['E1F_2PI_ANADIR1']
    elif opt == "--e1fs2":
        anadir = os.environ['E1F_2PI_ANADIR2']
    elif opt == "--hel":       
        if arg=='POS': 
            csvfname = 'simstats_vm_POS.csv'
            outdirname = 'simstats_POS'
        elif arg=='NEG': 
            csvfname = 'simstats_vm_NEG.csv'
            outdirname = 'simstats_NEG'
        else:
            print arg,' is not recognized'
            sys.exit()
            
csvf = os.path.join(anadir,csvfname)
print 'anadir = ', anadir
print 'outdirname = ',outdirname
print 'csvf = ', csvf

#anadir = os.environ['E1F_2PI_ANADIR1']
#anadir = os.environ['E1F_2PI_ANADIR2']
filename = os.path.join(anadir, "simstats_vm.csv")
df = pd.read_csv(csvf, na_values=['\n'])
#print df

#group df by q2wbinnum column
df_grpd_q2wbinnum = df.groupby('q2wbinnum')
nq2wbins = len(df_grpd_q2wbinnum.groups)
print 'nq2wbins=',nq2wbins

plot_track_stats()
plot_latest_stats()