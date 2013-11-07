# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ## Lists and Arrays in Python 
# * ### Arrays are made of of Lists
# * ### Lists are 1D "horizontal" data structures

# <rawcell>

# Example of List & implementation of 2D data structure using List

# <codecell>

l = [] #Create an empty List
print "Empty list = ",l
l = [1,2,3] #Initialize List
print "Initialized list =",l

md = [[],[]] #Create a List of Lists (="2D Array")
print "Empty \"2D Array\" = ",md
md[0]=[1,2,3] #Initialize "2D Array"=(Lists of Lists)
md[1]=[4,5,6]
print "Initiliazed \"2D Array\" = ",md
#print md[1][1] #Access a particular "2D Array" element


# <rawcell>

# Direct use of Numpy to create Arrays

# <codecell>

#from numpy import * 
import numpy as np
a = np.arange(10).reshape(2, 5)
print "array A shape:size =",a.shape,":",a.ndim
print a

b = np.array([(1,2,1),(3,4,1),(3,4,1),(4,5,1)])
print "array B shape:size =",b.shape,":",b.ndim
print b

c = np.zeros((5,5,5))
print "array C shape:size =",c.shape,":",c.ndim
print c
c[0][0][0] = 1
print c[0][0][0]

# <codecell>

#print b
bt = np.transpose(b)
print bt 
np.vsplit(bt,3)
#vs = np.vsplit(bt,2)
#hs = np.hsplit(bt,(0,1,2))
#print vs
#print hs

#print vs[0][0]
#print hs[0][0]

# <codecell>

a = floor(10*random.random((12,1)))
print a

# <codecell>

np.hsplit(a,1)

# <codecell>


