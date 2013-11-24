from __future__ import division
import numpy as npy
import os
from pandas import *
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
pandas.set_printoptions(max_columns=20)

anadir = os.environ['E1F_2PI_ANADIR1']
filename = os.path.join(anadir, "simstats_vm.csv")
df = pandas.read_csv(filename, na_values=['\n'])
#print df

#group df by q2wbinnum column
df_q2wbinnum = df.groupby('q2wbinnum')
nq2wbins = len(df_q2wbinnum.groups)
print 'nq2wbins=',nq2wbins

#stat = columns of df [nFB_ST,nEB_SR,nFB_SR,nEB_SA]
#tdraw(stat) plots, for each topology & latest simulation:
#  1. stat(q2wbin) Vs. q2wbin
#  2. frc(q2wbin) Vs. q2wbin, where frc(q2wbin)=stat(q2wbin)/nFB_ST(q2wbin)

MAX_TOPS,MIN_TOPS = range(2)
NTOPS=5
tops = [1,2,3,4,5];
top_labels = ['Top1','Top2','Top3','Top4','Top5']
top_markers = ['s','<','8' ,'p','>']
top_colors =  ['r','k','r','r','k']

NSTATS=4
nFB_ST,nEB_SR,nFB_SR,nEB_SA = range(0,NSTATS)
STATS_NAME = ['nFB_ST','nEB_SR','nFB_SR','nEB_SA']
def plot_track_stats():
    for q2wbinnum in df_q2wbinnum.groups:
        #iq2wbin = q2wbinnum-1
        #if q2wbinnum>1: continue
        df = df_q2wbinnum.get_group(q2wbinnum)
        q2wbinname = df['q2wbin'].tolist()[0]
        print 'q2wbinname=',q2wbinname
        for stat in range(0,NSTATS):
            print STATS_NAME[stat]
            fig=plt.figure('%d:%s'%(q2wbinnum,STATS_NAME[stat]))
            ax = []
            ax.append(fig.add_subplot(1,1,1,title='%s:%s'%(STATS_NAME[stat],q2wbinname)))
            ax.append(ax[MAX_TOPS].twinx())
            #print len(ax)
            ax[MAX_TOPS].set_ylabel('%s t2,5'%STATS_NAME[stat])
            ax[MIN_TOPS].set_ylabel('%s t1,3,4'%STATS_NAME[stat],color='r')
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
                    
            #save plots
            outdir = os.path.join(anadir,'simstats.new',q2wbinname)
            if not os.path.isdir(outdir):
                os.makedirs(outdir)
            #plt.figure(stat)
            fig.savefig('%s/%s.jpg'%(outdir,STATS_NAME[stat]))
            del fig
            del ax
            
plot_track_stats()