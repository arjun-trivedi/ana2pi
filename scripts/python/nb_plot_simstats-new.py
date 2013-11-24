# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ### This notebook has various tests for the DataFrame object made for Simulation Statistics

# <codecell>

from __future__ import division
import numpy as npy
import os
from pandas import *
from matplotlib.lines import Line2D
pandas.set_printoptions(max_columns=20)

#load .csv file and create DataFrame for simstats(=dfss)
anadir = os.environ['E1F_2PI_ANADIR1']
filename = os.path.join(anadir, "simstats_vm.csv")
dfss = pandas.read_csv(filename, na_values=['\n'])
#print dfss

#group dfss by q2wbinnum column
dfss_grpd_q2wbinnum = dfss.groupby('q2wbinnum')
nq2wbins = len(dfss_grpd_q2wbinnum.groups)
print 'nq2wbins=',nq2wbins

#stat = columns of df [nFB_ST,nEB_SR,nFB_SR,nEB_SA]
#tdraw(stat) plots, for each topology & latest simulation:
#  1. stat(q2wbin) Vs. q2wbin
#  2. frc(q2wbin) Vs. q2wbin, where frc(q2wbin)=stat(q2wbin)/nFB_ST(q2wbin)

MAX_TOPS,MIN_TOPS = range(2)
ntops=5
tops = [1,2,3,4,5];
top_labels = ['Top1','Top2','Top3','Top4','Top5']
top_markers = ['s','<','8' ,'p','>']
top_colors =  ['r','k','r','r','k']
def plot_latest_stats(stat):
    ax = []
    ax.append(plt.subplot(1,1,1,ylabel=stat))
    ax.append(ax[MAX_TOPS].twinx())
    ax[MAX_TOPS].set_ylabel('%s t2,5'%stat)
    ax[MIN_TOPS].set_ylabel(('%s top1,3,4'%stat),color='r')
    for tl in ax[MIN_TOPS].get_yticklabels():
        tl.set_color('r')
    #m_nbins = npy.zeros((nq2wbins,ntops));
    m_nbins = npy.zeros((ntops,nq2wbins));
    l_q2wb = npy.arange(1,nq2wbins+1,1);
    for q2wbinnum in dfss_grpd_q2wbinnum.groups:
        #print 'q2wbinnum=',q2wbinnum
        iq2wbin = q2wbinnum-1
        df = dfss_grpd_q2wbinnum.get_group(q2wbinnum)
        siml = df['Sim'].max()
        for top in tops:
            itop = top-1
            sel = (df['Sim']==siml) & (df['Top']==(top))
            m_nbins[itop][iq2wbin] = df[stat][sel]
    
    lns = []
    for top in tops:
        itop = top-1
        l_nbins=m_nbins[itop]
        if top==2 or top==5:
            lns.append(ax[MAX_TOPS].scatter(l_q2wb, l_nbins, color='k', marker=top_markers[itop], s=50,
                           label=top_labels[itop]))
        else:
            lns.append(ax[MIN_TOPS].scatter(l_q2wb, l_nbins, color='r', marker=top_markers[itop], s=50,
                           label=top_labels[itop]))    
    labs = [l.get_label() for l in lns]
    ax[MAX_TOPS].legend(lns, labs, loc=2)
    plt.show()

plot_latest_stats('nEB_SR')
plot_latest_stats('nFB_SR')
plot_latest_stats('nFB_ST')

# <codecell>

def plot_track_stats(stat):
    for q2wbinnum in dfss_grpd_q2wbinnum.groups:
        iq2wbin = q2wbinnum-1
        if q2wbinnum>1: continue
        df = dfss_grpd_q2wbinnum.get_group(q2wbinnum)
        q2wbinname=df.loc[0]['q2wbin']
        fig=plt.figure(stat)
        ax = []
        ax.append(fig.add_subplot(1,1,1,ylabel=stat,title=q2wbinname))
        ax.append(ax[MAX_TOPS].twinx())
        ax[MAX_TOPS].set_ylabel('%s t2,5'%stat)
        ax[MIN_TOPS].set_ylabel(('%s top1,3,4'%stat),color='r')
        for tl in ax[MIN_TOPS].get_yticklabels():
            tl.set_color('r')
        
        #nsims = len(df['Sim'])
        #m_nbins = npy.zeros((ntops,nsims));
        #l_sims = npy.arange(1,nsims+1,1);
        #sel_maxtops = (df['Top']==(2)) | (df['Top']==(5))
        #sel_mintops = (df['Top']==(1)) | (df['Top']==(3)) | (df['Top']==(4))
        lns = []
        for top in tops:
            itop=top-1
            if top==2 or top==5:
                
                #print df['Sim'][df['Top']==top]
                sel= (df['Top']==top)
                lns.append(ax[MAX_TOPS].scatter(df['Sim'][sel],df[stat][sel],color='k',
                                     marker=top_markers[itop],label=top_labels[itop]))
            else:
                sel= (df['Top']==top)
                lns.append(ax[MIN_TOPS].scatter(df['Sim'][sel],df[stat][sel],color='r',
                                     marker=top_markers[itop],label=top_labels[itop]))
        labs = [l.get_label() for l in lns]
        ax[MAX_TOPS].legend(lns, labs, loc=2)
        #plt.show()
        
        #save plots
        outdir = os.path.join(anadir,'simstats.new',q2wbinname)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        #plt.figure(stat)
        fig.savefig('%s/%s.jpg'%(outdir,stat))
        

plot_track_stats('nFB_SR')
plt.show()
    

# <codecell>

NSTATS=4
nFB_ST,nEB_SR,nFB_SR,nEB_SA = range(0,NSTATS)
STATS_NAME = ['nFB_ST','nEB_SR','nFB_SR','nEB_SA']
def plot_track_stats():
    for q2wbinnum in dfss_grpd_q2wbinnum.groups:
        iq2wbin = q2wbinnum-1
        #if q2wbinnum>1: continue
        df = dfss_grpd_q2wbinnum.get_group(q2wbinnum)
        q2wbinname = df['q2wbin'].tolist()[0]
        for stat in range(0,NSTATS):
            print STATS_NAME[stat]
            fig=plt.figure(STATS_NAME[stat])
            ax = []
            ax.append(fig.add_subplot(1,1,1,ylabel=stat,title='%s:%s'%(STATS_NAME[stat],q2wbinname)))
            ax.append(ax[MAX_TOPS].twinx())
            ax[MAX_TOPS].set_ylabel('%s t2,5'%STATS_NAME[stat])
            ax[MIN_TOPS].set_ylabel(('%s top1,3,4'%STATS_NAME[stat]),color='r')
            for tl in ax[MIN_TOPS].get_yticklabels():
                tl.set_color('r')
        
            lns = []
            for top in tops:
                itop=top-1
                if top==2 or top==5:
                    #print df['Sim'][df['Top']==top]
                    sel= (df['Top']==top)
                    lns.append(ax[MAX_TOPS].scatter(df['Sim'][sel],df[STATS_NAME[stat]][sel],color='k',
                                     marker=top_markers[itop],label=top_labels[itop]))
                else:
                    sel= (df['Top']==top)
                    lns.append(ax[MIN_TOPS].scatter(df['Sim'][sel],df[STATS_NAME[stat]][sel],color='r',
                                     marker=top_markers[itop],label=top_labels[itop]))
            labs = [l.get_label() for l in lns]
            ax[MAX_TOPS].legend(lns, labs, loc=2)
             #plt.show()
        
            #save plots
            outdir = os.path.join(anadir,'simstats.new',q2wbinname)
            if not os.path.isdir(outdir):
                os.makedirs(outdir)
            #plt.figure(stat)
            fig.savefig('%s/%s.jpg'%(outdir,STATS_NAME[stat]))
            
plot_track_stats()
        

# <codecell>


