# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ### This notebook has various tests for the DataFrame object made for Simulation Statistics
# - Simulation Statistics file = anae1fd2/simstats_vm.csv    

# <markdowncell>

# ### First load the .csv file and create the DataFrame for simstats(=dfss)

# <codecell>

import os
from pandas import *
pandas.set_printoptions(max_columns=20)

datadir = os.environ['E1F_2PI_ANADIR1']
filename = os.path.join(datadir, "simstats_vm.csv.bk.112213")
dfss = pandas.read_csv(filename, na_values=['\n'])
#print df

# <markdowncell>

# ### Group dfss by q2wbin column

# <codecell>

dfss_grpd_q2wbin = dfss.groupby('Q2Wbin')
nq2wbins = len(dfss_grpd_q2wbin.groups)
ngridx = 4
ngridy = nq2wbins/ngridx
print ngridx, ngridy

# <markdowncell>

# ### Test gr_nFthVq2wbin

# <codecell>

from __future__ import division
import numpy as npy

#create Figure
plt.figure(figsize=(8,8))
# create Axes that go in the Figure
ax_nFTH = plt.subplot(3,1,1,ylabel = 'nFth') 
ax_nESR = plt.subplot(3,1,2,ylabel='nEsr')
ax_frc = plt.subplot(3,1,3,xlabel='q2wbin',ylabel='frc')
sels = {}
ntops=5
tops = [1,2,3,4,5];
nFTH = npy.zeros((nq2wbins,ntops));
nESR = npy.zeros((nq2wbins,ntops));
frc = npy.zeros((nq2wbins,ntops));
q2wb = npy.arange(1,nq2wbins+1,1);
for q2wbin in dfss_grpd_q2wbin.groups:
    iq2wbin = q2wbin-1
    #print 'q2wbin', q2wbin
    df = dfss_grpd_q2wbin.get_group(q2wbin)
    #print 'sim:'
    #print df['Sim']
    
    siml = df['Sim'].max()
    for top in tops:
        itop = top-1
        sels[itop] = (df['Sim']==siml) & (df['Top']==(itop+1))
        nFTH[iq2wbin][itop] = df.nFth[sels[itop]]
        nESR[iq2wbin][itop] = df.nEsr[sels[itop]]
        frc[iq2wbin][itop] = nESR[iq2wbin][itop]/nFTH[iq2wbin][itop]

nFTH_trps = numpy.vsplit(npy.transpose(nFTH),ntops)
nESR_trps = numpy.vsplit(npy.transpose(nESR),ntops)
frc_trps = numpy.vsplit(npy.transpose(frc),ntops)
#print frc
#print frc[0]
#print q2wb

colors = ['r','g','b','k','y']
for top in tops:
    itop = top-1
    ax_nFTH.scatter(q2wb, nFTH_trps[itop], color=colors[itop])
    ax_nESR.scatter(q2wb, nESR_trps[itop], color=colors[itop])
    ax_frc.scatter(q2wb, frc_trps[itop], color=colors[itop])
    
#ax1.scatter(q2wb, tp[0])
#ax2.scatter(q2wb,ne)
#ax3.scatter(q2wb,frc)
#ax3.axes.xaxis.label = 'Q2WBIN'
#plt.show()

# <codecell>

#create Figure
from matplotlib.lines import Line2D
plt.figure("tfig",figsize=(8,8))
# create Axes that go in the Figure
#ax_nFTH = plt.subplot(3,1,1,ylabel = 'nFth')
ax_nESR = []
ax_nESR.append(plt.subplot(1,1,1,ylabel='nEsr'))
ax_nESR[0].set_ylabel('nESR top1,3&4')
#ax_nESR[0].legend(loc=0)
ax_nESR.append(ax_nESR[0].twinx())
ax_nESR[1].set_ylabel('nESR top2&5',color='r')
#ax_nESR[1].legend(loc=0)
for tl in ax_nESR[1].get_yticklabels():
    tl.set_color('r')

q2wb = npy.arange(1,nq2wbins+1,1);
for q2wbin in dfss_grpd_q2wbin.groups:
    iq2wbin = q2wbin-1
    #print 'q2wbin', q2wbin
    df = dfss_grpd_q2wbin.get_group(q2wbin)
    #print 'sim:'
    #print df['Sim']
    
    siml = df['Sim'].max()
    for top in tops:
        itop = top-1
        sels[itop] = (df['Sim']==siml) & (df['Top']==(itop+1))
        #nFTH[iq2wbin][itop] = df.nFth[sels[itop]]
        #nESR[iq2wbin][itop] = df.nEsr[sels[itop]]
        nESR[iq2wbin][itop] = df['nEsr'][sels[itop]]
        #frc[iq2wbin][itop] = nESR[iq2wbin][itop]/nFTH[iq2wbin][itop]

nESR_trps = numpy.vsplit(npy.transpose(nESR),ntops)
#frc_trps = numpy.vsplit(npy.transpose(frc),ntops)
#print frc
#print frc[0]
#print q2wb

labels = ['Top1','Top2','Top3','Top4','Top5']
markers = ['s','<','8' ,'p','>']
colors =  ['b','r','g','y','r']
linestyles = ['_', '-', '--', ':']
lns = []
for top in tops:
    itop = top-1
    if top==2 or top==5:
        #print markers[itop]
        lns.append(ax_nESR[0].scatter(q2wb, nESR_trps[itop], color='k', marker=markers[itop], s=50,
                           label=labels[itop]))
    else:
        lns.append(ax_nESR[1].scatter(q2wb, nESR_trps[itop], color='r', marker=markers[itop], s=50,
                           label=labels[itop]))    
    #ax_nFTH.scatter(q2wb, nFTH_tp[itop], color=colors[itop])
    #ax_nESR.scatter(q2wb, nESR_tp[itop], color=colors[itop])
    #ax_frc.scatter(q2wb, frc_tp[itop], color=colors[itop])
#p1, = plot([1,2,3], label="test1")
#p2, = plot([3,2,1], label="test2")

#l1 = legend([ax_nESR[1],ax_nESR[0]], ["Label 1"], loc=1)
#l2 = legend([ax_nESR[0]], ["Label 2"], loc=4) # this removes l1 from the axes.
#gca().add_artist(l1) # add l1 as a separate artist to the axes

#plt.legend()#, ["Label 1"], loc=1)
#plt.figure("tfig")
labs = [l.get_label() for l in lns]
ax_nESR[0].legend(lns, labs, loc=0)
plt.show()

# <codecell>

#seq = [ST,SR,SA,SC,SH,SF], [ER,EC,EF]
#nbins = nbins(seq)
#tdraw(seq) plots, for each topology:
#  1. nbins(seq) Vs. q2wbin
#  2. frc(seq) Vs. q2wbin, where frc(seq)=nbins(seq)/nbins(ST)
def tdraw(seq):
    ax = []
    ax.append(plt.subplot(1,1,1,ylabel=seq))
    ax.append(ax[0].twinx())
    ax[0].set_ylabel('%s t2,5'%seq)
    ax[1].set_ylabel(('%s top1,3,4'%seq),color='r')
    for tl in ax[0].get_yticklabels():
        tl.set_color('r')
    nbins = npy.zeros((nq2wbins,ntops));
    q2wb = npy.arange(1,nq2wbins+1,1);
    for q2wbin in dfss_grpd_q2wbin.groups:
        iq2wbin = q2wbin-1
        df = dfss_grpd_q2wbin.get_group(q2wbin)
        siml = df['Sim'].max()
        for top in tops:
            itop = top-1
            sels[itop] = (df['Sim']==siml) & (df['Top']==(itop+1))
            nbins[iq2wbin][itop] = df[seq][sels[itop]]
    
    nbins_trps = numpy.vsplit(npy.transpose(nbins),ntops)
    
    labels = ['Top1','Top2','Top3','Top4','Top5']
    markers = ['s','<','8' ,'p','>']
    colors =  ['b','r','g','y','r']
    linestyles = ['_', '-', '--', ':']
    lns = []
    for top in tops:
        itop = top-1
        if top==2 or top==5:
            lns.append(ax[1].scatter(q2wb, nbins_trps[itop], color='r', marker=markers[itop], s=50,
                           label=labels[itop]))
        else:
            lns.append(ax[0].scatter(q2wb, nbins_trps[itop], color='k', marker=markers[itop], s=50,
                           label=labels[itop]))    
    labs = [l.get_label() for l in lns]
    ax[0].legend(lns, labs, loc=0)
    plt.show()
    
tdraw('nEsr')
tdraw('nFsr')
tdraw('nFth')

# <markdowncell>

# ### Test Division Operation between two different columumns of D

# <codecell>

#from __future__ import division
#for q2wbin in dfss_grpd_q2wbin.groups:
#    df = dfss_grpd_q2wbin.get_group(q2wbin)
#    if (q2wbin==1):
#        #x= df.nEsr/df.nFth
#        #print x
#        dfsel = df[(df['Top']==1) & (df['Varset']==1)]
#        num = dfsel.nEsr.values
#        print num
#        #num.reindex
#        #print num
#        
#        den = df[(df['Top']==2) & (df['Varset']==1)].nEsr.values
#        print den
#        #den.reindex
#        #print den
#        
#        x = num/den
#        print x
#        
#        print dfsel
#        print df[(df['Top']==2) & (df['Varset']==1)]

# <codecell>


