# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ### Evan's gpp parameters
# 
# * /group/clas/builds/bin/gpp -P0x1f -a1.357 -b1.357 -c1.357 -f0.962 -Y -R38121 -oGPP.OUT ./GSIM.OUT

# <codecell>

from anah10 import *

fold = '/e1f.sim2pi.datadir/comp_old_new_sim/Q2W__1.9-2.5__1.3-1.9/cooked/1.root'
fnew = '/e1f.sim2pi.datadir/sim_range_study/q2w2/cooked/1.root'

dfs = import_h10([fold,fnew])
plot_evt_sub_cols(dfs)

# <codecell>

x = randn(10000)

