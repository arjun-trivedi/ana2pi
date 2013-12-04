# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from comp_sim import *
at_hw()
d={}
OLD,NEW = range(0,2)
MCTK,COOKED = range(0,2)
d[OLD,COOKED],d[NEW,COOKED]=import_data()

#addcol_qsq_w(dn)

# <codecell>

print d[OLD,COOKED].shape
print d[NEW,COOKED].shape

# <codecell>

addcol_qsq_w(d[OLD,COOKED])
addcol_qsq_w(d[NEW,COOKED])

# <codecell>

fig = plt.figure(figsize(5,5))
fig.suptitle('old-new sim comparison')
ax_w=fig.add_subplot(211)
ax_w.set_title('W')
ax_w.set_xlabel('W(GeV)')
n,bins,patches=plt.hist(d[OLD,COOKED].w,100,(1.0,2.0),
                        histtype='step',color='black',label='old')
n,bins,patches=plt.hist(d[NEW,COOKED].w,100,(1.0,2.0),
                        histtype='step',color='red',label='new')
ax_qsq=fig.add_subplot(212)
ax_qsq.set_title('Q2')
ax_qsq.set_xlabel('Q2(GeV)')
n,bins,patches=plt.hist(d[OLD,COOKED].qsq,100,(1.0,2.0),
                        histtype='step',color='black',label='old')
n,bins,patches=plt.hist(d[NEW,COOKED].qsq,100,(1.0,2.0),
                        histtype='step',color='red',label='new')
plt.legend()

# <codecell>


