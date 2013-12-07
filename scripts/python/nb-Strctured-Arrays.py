# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ### Test Structured Arrays

# <codecell>

from comp_sim import *
fold = '/e1f.sim2pi.datadir/comp_old_new_sim/Q2W__1.9-2.5__1.3-1.9/cooked/1.root'
arr = root2array(fold,'h10',start=1,stop=5)#start=1,stop=2)

#arr['gpart']
print ">>>arr['gpart']:"
print arr['gpart']
#print ">>>arr['etot']:"
#print arr['etot']

tarr = np.array(arr['etot'])#ndmin=2)
print tarr
print tarr.shape
print len(tarr)



#print tarr[0][0]


# <codecell>

df = pd.DataFrame(arr)
print df['mcp']

# <codecell>

x = np.zeros((2,),dtype=('i4,f4,a10'))
x[:] = [(1,2.,'Hello'),(2,3.,"World")]
x

# <codecell>

y = np.zeros(3, dtype='3int8, float32, (2,3)float64')
y

# <codecell>

y[0] = ([1,2,3],4,[[1,2,3],[4,5,6]])
y[1] = ([10,20,30],40,[[10,20,30],[40,50,60]])
y[2] = ([20,40,60],80,[[20,40,20],[80,90,90]])

# <codecell>

print y['f0']
print y.dtype.fields['f0']

# <codecell>

print y['f0']
t = np.array(y['f0'])
print t.shape
print t[:,0:1]
print t
n = plt.hist(t)

# <codecell>

m =np.zeros((2,2))
m[:,1]=1
print m

# <codecell>

a=b=[]

t = [a,b]
a.append(1)
b.append(2)
t[0][0]

# <codecell>

wp = pd.Panel(randn(2, 5, 4), items=['Item1', 'Item2'],
               major_axis=['1','2','3','4','5',],
               minor_axis=['A', 'B', 'C', 'D'])
wp
#for item in len(wp.items):
    #wp[item]

# <codecell>


