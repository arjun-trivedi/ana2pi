# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ### ye's gpp parameters
# 
# * /group/clas/builds/bin/gpp -P0x1f -oGPP.OUT ./GSIM.OUT

# <codecell>

from anah10 import *

fold = '/e1f.sim2pi.datadir/comp_old_new_sim/Q2W__1.9-2.5__1.3-1.9/cooked/1.root'
fnew = '/e1f.sim2pi.datadir/comp_old_new_sim/q2w2_ye_gpp/cooked/1.root'

dfs = import_h10([fold,fnew])
plot_evt_sub_cols(dfs)

# <codecell>

dfn=ph10['NEW']
dfo=ph10['OLD']
no=plt.hist(dfo.el_cz,100,[0.9,1],histtype='step',color='black')
nn=plt.hist(dfn.el_cz,100,[0.9,1],histtype='step',color='red')
print mean(no[0])
print mean(nn[0])
#val = dfo.el_cz
#print val
#print mean(val)
#print mean(dfn.el_cz)

# <codecell>

n=plt.hist(dfo.el_cx,100,[0,1],histtype='step',color='black')
n=plt.hist(dfn.el_cx,100,[0,1],histtype='step',color='red')
print mean(dfo.el_cx)
print mean(dfn.el_cx)

# <codecell>


