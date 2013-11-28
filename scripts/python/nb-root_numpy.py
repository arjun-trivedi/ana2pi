# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <rawcell>

# imports

# <codecell>

import ROOT
from root_numpy import root2array, root2rec, tree2rec
from root_numpy.testdata import get_filepath



# <rawcell>

# Decide on example .root file to use

# <codecell>

#filename = get_filepath('test.root')
filename='/e1f.2pi.datadir1/h10/qskim_38462.root'

# <codecell>

# Get the TTree from the ROOT file
rfile = ROOT.TFile(filename)
#intree = rfile.Get('tree')
intree=rfile.Get('h10')
print intree.GetEntries()
print len(intree.GetListOfBranches())

# <codecell>

# Convert a TTree in a ROOT file into a NumPy structured array
#arr = root2array(filename, 'tree')
arr = root2array(filename, 'h10')
arr
print arr.size
print len(arr[0])

# <codecell>


