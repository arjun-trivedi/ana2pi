# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ## Lists and Arrays in Python (Arrays are made up of Lists in Python)

# <markdowncell>

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


# <codecell>

import matplotlib.pyplot as plt

x = range(100)
y = range(100,200)
fig = plt.figure()

ax1 = fig.add_subplot(221,title="ax1")
opts = 'c=\'b\', marker=\"s\"'
print opts
#ax1.scatter(x[:4], y[:4], s=10, opts)
ax1.scatter(y[4:], x[4:], s=10, c='r', marker="o")


ax2 = fig.add_subplot(222, title="ax2")
ax2.scatter(x[:4], y[:4], s=1, c='b', marker="s")
ax2.scatter(y[4:], x[4:], s=1, c='r', marker="o")

#plt.show()
plt.savefig("test2.jpg")

# <codecell>


