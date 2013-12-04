# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <rawcell>

# - From a list of ROOT files, gather h10 Trees -> Structured Array -> Pandas DataFrame (DF)
# - Then perform analysis using DF

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
#arr = root2array('/e1f.2pi.datadir1/h10/qskim_3860*.root','h10')
arr = root2array('/data/trivedia/e1f/simulation_2pi/sim_range_study/q2wA/cooked/*.root',
                 'h10')
print 'Number of Branches in Tree = ',len(arr[0])
print 'Number of entries in Tree = ',arr.size

#Convert structured array to dataframe
import pandas as pd
d = pd.DataFrame(arr)

# <codecell>

MASS_E = 0.000511
MASS_P = 0.93827203
E1F_P = 5.499; 
lvE0 = ROOT.TLorentzVector(0,0,E1F_P,sqrt(E1F_P*E1F_P+MASS_E*MASS_E));
lvP0 = ROOT.TLorentzVector(0,0,0,MASS_P);

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
    lvq = lvE0-lve
    lvw = lvq+lvP0
    
    qsq.append(-lvq.Mag2())
    w.append(lvw.Mag())
 
#Directly plot from qsq[]
#n = plt.hist(qsq,20,(-0.9,0.9))

#Or add qsq[] to DF and plot using DF
d['qsq']=qsq 
d['w']=w
#n=plt.hist(d.qsq,20,(1.0,2.0))
#n=plt.hist2d(w,qsq)


#print 'bins = ',bins
#print 'num bins = ',len(bins)

#h = numpy.histogram(e_p,100,(0.0,5.0))
#ax = plt.hist(h)

#plt.show()
#d.p[0][0]
#d.hist('p[0][0]')

# <codecell>

#n=plt.hist(d.qsq,20,(1.0,2.0))

h2,x,y = np.histogram2d(d.qsq,d.w,10,[[1,2],[1,2]])
#h2,x,y = np.histogram2d(d.w,d.w,bins=10)
#h2,x,y = np.histogram2d(d.qsq,d.qsq,20,[[1,2],[1,2]])
#h2,x,y = np.histogram2d([2,2,2,2],[1,2,3,4],3,[[2,5],[2,5]])
print x
print y
#print h2
#plt.imshow(h2)
plt.pcolormesh(x,y,h2)

# <codecell>


# <codecell>


