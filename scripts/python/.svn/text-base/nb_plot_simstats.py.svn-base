# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os
from pandas import *
pandas.set_printoptions(max_columns=20)

datadir = os.environ['E1F_2PI_ANADIR1']
filename = os.path.join(datadir, "simstats_vm.csv")
dfss = pandas.read_csv(filename, na_values=['\n'])
#print df

# <codecell>

dfss_grpd_q2wbin = dfss.groupby('Q2Wbin')
nq2wbins = len(dfss_grpd_q2wbin.groups)
ngridx = 4
ngridy = nq2wbins/ngridx
print ngridx, ngridy

# <markdowncell>

# ### Test Division between columns of two different DFs

# <codecell>

from __future__ import division
for q2wbin in dfss_grpd_q2wbin.groups:
    df = dfss_grpd_q2wbin.get_group(q2wbin)
    if (q2wbin==1):
        #x= df.nEsr/df.nFth
        #print x
        dfsel = df[(df['Top']==1) & (df['Varset']==1)]
        num = dfsel.nEsr.values
        print num
        #num.reindex
        #print num
        
        den = df[(df['Top']==2) & (df['Varset']==1)].nEsr.values
        print den
        #den.reindex
        #print den
        
        x = num/den
        print x
        
        print dfsel
        print df[(df['Top']==2) & (df['Varset']==1)]

# <codecell>

itops = {0,1,2,3,4} #NOTE, that these are indices of tops, not tops; vm_tops used for now
colors = ['b','r','g','m','k']
markers = ['s','o','^','v','<']
labels = ['Top1','Top2','Top3','Top4','Top5']
sels  = {}
gr_Fth = {}
gr_Fsr = {}
gr_Esr = {}
gr_EsrVsEac = {}


for q2wbin in dfss_grpd_q2wbin.groups:
    df = dfss_grpd_q2wbin.get_group(q2wbin)
    
    fg_Fth_name  = str.format('Fth_%d'%q2wbin)
    fg_Fth_title = str.format('Fth[Q2W=%d]'%q2wbin)
    fg_Fth = plt.figure(fg_Fth_name,figsize=(6,5))
    ax_Fth = fg_Fth.add_subplot(1,1,1,title=fg_Fth_title,xticks=df.Sim, xlabel=df.Sim.name,ylabel=df.nFth.name)
    
    fg_Fsr_name  = str.format('Fsr_%d'%q2wbin)
    fg_Fsr_title = str.format('Fsr[Q2W=%d]'%q2wbin)
    fg_Fsr = plt.figure(fg_Fsr_name,figsize=(6,5))
    ax_Fsr = fg_Fsr.add_subplot(1,1,1,title=fg_Fsr_title,xticks=df.Sim, xlabel=df.Sim.name,ylabel=df.nFsr.name)
    
    fg_EsrVsEac_name  = str.format('EsrVsEac_%d'%q2wbin)
    fg_EsrVsEac_title = str.format('EsrVsEac[Q2W=%d]'%q2wbin)
    fg_EsrVsEac = plt.figure(fg_EsrVsEac_name,figsize=(6,5))
    ax_EsrVsEac = fg_EsrVsEac.add_subplot(1,1,1,title=fg_EsrVsEac_title,xlabel=df.nEac.name,ylabel=df.nEsr.name)
    
    
    for itop in itops:
        sels[itop] = (df['Top']==itop+1) & (df['Varset']==1)
        dfsel = df[sels[itop]]
        gr_Fth[itop] = ax_Fth.scatter(dfsel.Sim, dfsel.nFth,
                                      s=50,c=colors[itop],marker=markers[itop],label=labels[itop])
        gr_Fsr[itop] = ax_Fsr.scatter(dfsel.Sim, dfsel.nFsr,
                                      s=50,c=colors[itop],marker=markers[itop],label=labels[itop])
        gr_EsrVsEac[itop] = ax_EsrVsEac.scatter(dfsel.nEac,dfsel.nEsr,
                                                s=50,c=colors[itop],marker=markers[itop],label=labels[itop])
    
    plt.figure(fg_Fth_name)
    plt.legend()
    fg_Fth.savefig("%s.jpg"%fg_Fth_name)
    plt.figure(fg_Fsr_name)
    plt.legend()
    fg_Fsr.savefig("%s.jpg"%fg_Fsr_name)
    plt.figure(fg_EsrVsEac_name)
    plt.legend()
    fg_EsrVsEac.savefig("%s.jpg"%fg_EsrVsEac_name)

# <codecell>


