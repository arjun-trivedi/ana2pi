# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from comp_sim import *
at_hw()
#d={}
OLD,NEW = range(0,2)
MCTK,COOKED = range(0,2)
#d[OLD,COOKED],d[NEW,COOKED]=import_data()
p = import_data('/e1f.sim2pi.datadir/comp_old_new_sim/Q2W__1.9-2.5__1.3-1.9/cooked/1.root',
                '/e1f.sim2pi.datadir/sim_range_study/q2w2/cooked/1.root')

#addcol_qsq_w(dn)

# <codecell>

print p['OLD'].shape
print p['NEW'].shape
#print d[NEW,COOKED].shape

# <codecell>

addcol_qsq_w(p['OLD'])
addcol_qsq_w(p['NEW'])

# <codecell>

fig = plt.figure(figsize(5,5))
fig.suptitle('old-new sim comparison')
ax_w=fig.add_subplot(211)
ax_w.set_title('W')
ax_w.set_xlabel('W(GeV)')
n,bins,patches=plt.hist(p['OLD'].w,100,(1.0,2.0),
                        histtype='step',color='black',label='old')
n,bins,patches=plt.hist(p['NEW'].w,100,(1.0,2.0),
                        histtype='step',color='red',label='new')
ax_qsq=fig.add_subplot(212)
ax_qsq.set_title('Q2')
ax_qsq.set_xlabel('Q2(GeV)')
n,bins,patches=plt.hist(p['OLD'].qsq,100,(1.0,2.0),
                        histtype='step',color='black',label='old')
n,bins,patches=plt.hist(p['NEW'].qsq,100,(1.0,2.0),
                        histtype='step',color='red',label='new')
plt.legend()

# <codecell>

comp_qsq_w()

# <codecell>


