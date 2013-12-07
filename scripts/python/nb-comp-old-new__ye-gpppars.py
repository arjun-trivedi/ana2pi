# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ### ye's gpp parameters
# 
# * /group/clas/builds/bin/gpp -P0x1f -oGPP.OUT ./GSIM.OUT

# <codecell>

from comp_sim import *

fold = '/e1f.sim2pi.datadir/comp_old_new_sim/Q2W__1.9-2.5__1.3-1.9/cooked/1.root'
fnew = '/e1f.sim2pi.datadir/comp_old_new_sim/q2w2_ye_gpp/cooked/1.root'

ph10 = import_data(fold,fnew,10000)
prep_data();

comp_qsq_w()
comp_basic()

# <codecell>

dfn=ph10['NEW']
dfo=ph10['OLD']
n=plt.hist(dfo.el_cz,100,[0.9,1],histtype='step',color='black')
n=plt.hist(dfn.el_cz,100,[0.9,1],histtype='step',color='red')
#plt.

#plt.pcolormesh(x,y,h2)

# <codecell>


