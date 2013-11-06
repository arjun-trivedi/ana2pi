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

# <markdowncell>

# ### Test gr_nFthVq2wbin

# <codecell>

from __future__ import division
dfss_grpd_q2wbin = dfss.groupby('Q2Wbin')
nq2wbins = len(dfss_grpd_q2wbin.groups)
ngridx = 4
ngridy = nq2wbins/ngridx
print ngridx, ngridy

#gr_FthVq2wbin

for q2wbin in dfss_grpd_q2wbin.groups:
    print 'q2wbin', q2wbin
    df = dfss_grpd_q2wbin.get_group(q2wbin)
    #print 'sim:'
    #print df['Sim']
    
    siml = df['Sim'].max()
    sel = (df['Sim']==siml) & (df['Top']==2)
    #nFth_siml = df.nFth(sel)
    nFth_siml = df.nFth[sel]
    print siml,nFth_siml

# <codecell>


