# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <rawcell>

# From a list of ROOT files, gather h10 Trees -> Structured Array -> Pandas DataFrame

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

MASS_E = 0.000511
MASS_P = 0.93827203

e_p=[]
qsq = []
w = []
for i in arange(0,len(d.p)):
    p = d.p[i][0]
    px = p * d.cx[i][0]
    py = p * d.cy[i][0]
    pz = p * d.cz[i][0]
    lve = ROOT.TLorentzVector(0,0,0,0)
    e = sqrt(p*p+MASS_E*MASS_E) 
    lve.SetPxPyPzE(px,py,pz,e)
    qsq.append(lve.Mag())
 
#Directly plot from qsq[]
#n = plt.hist(qsq,20,(-0.9,0.9))

#Or add qsq[] to DF and plot using DF
d['qsq']=qsq 
n=plt.hist(d.qsq,20,(-0.9,0.9))


#print 'bins = ',bins
#print 'num bins = ',len(bins)

#h = numpy.histogram(e_p,100,(0.0,5.0))
#ax = plt.hist(h)

#plt.show()
#d.p[0][0]
#d.hist('p[0][0]')

# <codecell>

    

