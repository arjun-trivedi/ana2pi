# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ### Test Structured Arrays

# <codecell>

from comp_sim import *
fold = '/e1f.sim2pi.datadir/comp_old_new_sim/Q2W__1.9-2.5__1.3-1.9/cooked/1.root'
arr = root2array(fold,'h10',start=1,stop=2)

# <codecell>

arr

# <codecell>

print arr['gpart']
print arr['p']
print arr['p'][0]
#print arr['p'][1]
print arr['p'][0][0]
print arr['p'][0][1]

# <codecell>

x = np.zeros((2,),dtype=('i4,f4,a10'))
x
x[:] = [(1,2.,'Hello'),(2,3.,"World")]
x

# <codecell>

y = np.zeros(3, dtype='3int8, float32, (2,3)float64')
y

# <codecell>

y[0] = ([1,2,3],4,[[1,2,3],[4,5,6]])

# <codecell>

y

# <codecell>

y['f0'][0][0]

# <codecell>


