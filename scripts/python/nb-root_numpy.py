# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <rawcell>

# Data in `h10 Tree` from .root files -> Structured Array -> Pandas DataFrame

# <codecell>

import ROOT
from root_numpy import root2array, root2rec, tree2rec
from root_numpy.testdata import get_filepath

#filename = get_filepath('test.root')
#filename='/e1f.2pi.datadir1/h10/qskim_38462.root'

# Get the TTree from the ROOT file
#rfile = ROOT.TFile(filename)
#intree = rfile.Get('tree')
#intree=rfile.Get('h10')
#print 'Number of Branches in Tree = ',len(intree.GetListOfBranches())
#print 'Number of entries in Tree = ',intree.GetEntries()

# Convert a TTree in a ROOT file into a NumPy structured array
#arr = root2array(filename, 'tree')
#arr = root2array(filename, 'h10')
arr = root2array('/e1f.2pi.datadir1/h10/qskim_3860*.root','h10')
print 'Number of Branches in Tree = ',len(arr[0])
print 'Number of entries in Tree = ',arr.size

#Convert structured array to dataframe
import pandas as pd
d = pd.DataFrame(arr)

# <codecell>

e_p=[]
for i in arange(0,len(d.p)):
    e_p.append(d.p[i][0])
#print e_p   
h=numpy.histogram(e_p)
bins=numpy.linspace(0,5,100)
n = plt.hist(e_p,100,(0,5))

#print 'bins = ',bins
#print 'num bins = ',len(bins)

#h = numpy.histogram(e_p,100,(0.0,5.0))
#ax = plt.hist(h)

#plt.show()
#d.p[0][0]
#d.hist('p[0][0]')

# <codecell>

#chainlist = '/e1f.2pi.anadir1/h10.lst'
#rchain= ROOT.TChain("h10")
#fc = ROOT.TFileCollection("fileList", "", chainlist);
#rchain.AddFileInfoList(fc.GetList())
#print rchain.GetEntries()

#arr1 = root2array('/e1f.2pi.datadir1/h10/qskim_3860*','h10')

#import rootpy
#from rootpy.tree import Tree, TreeModel, TreeChain, FloatCol, IntCol
#from rootpy.io import root_open

import pylab as P

#
# The hist() function now has a lot more options
#

#
# first create a single histogram
#
mu, sigma = 200, 25
x = mu + sigma*P.randn(10000)

# the histogram of the data with histtype='step'
n, bins, patches = plt.hist(x, 50, normed=1, histtype='stepfilled')
P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
y = P.normpdf( bins, mu, sigma)
l = P.plot(bins, y, 'k--', linewidth=1.5)
P.show()

# <codecell>


