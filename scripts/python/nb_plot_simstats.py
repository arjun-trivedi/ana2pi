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

#create Figure
plt.figure(figsize=(8,8))
# create Axes that go in the Figure
ax1 = plt.subplot(3,1,1,ylabel = 'nFth') 
ax2 = plt.subplot(3,1,2,ylabel='nEsr')
ax3 = plt.subplot(3,1,3,xlabel='q2wbin',ylabel='frc')
#sels = {}
nf = [[],[],[],[]]
ne = [[],[],[],[]]
frc = [[],[],[],[]]
q2wb = []
for q2wbin in dfss_grpd_q2wbin.groups:
    #print 'q2wbin', q2wbin
    df = dfss_grpd_q2wbin.get_group(q2wbin)
    #print 'sim:'
    #print df['Sim']
    
    siml = df['Sim'].max()
    selt1 = (df['Sim']==siml) & (df['Top']==1)
    selt2 = (df['Sim']==siml) & (df['Top']==2)
    #nFth_siml = df.nFth(sel)
    #nFth = df.nFth[sels[2]]
    #nEsr = df.nEsr[sels[2]]
    f = nEsr/nFth
    #print siml,nFth_siml
    nf[0].append(df.nFth[sels[1]])
    ne[0].append(df.nEsr[sels[1]])
    frc.append(f)
    q2wb.append(q2wbin)
    
ax1.scatter(q2wb,nf)
ax2.scatter(q2wb,ne)
ax3.scatter(q2wb,frc)
#ax3.axes.xaxis.label = 'Q2WBIN'
plt.show()

# <codecell>

### Test Division Operation between two different columumns of DF

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

