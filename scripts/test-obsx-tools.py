# import os
# #from rootpy.io import root_open
# anadir =  os.environ['E1F_SIM2PI_ANADIR1']
# print anadir
# #myfile = root_open(anadir/2:1__1-2.000-2.400__24-1.300-1.900__pol__exp.root)

# import rootpy
# from rootpy.plotting import Hist
# hist = Hist(100, 0, 10)
# hist.Fill(3) # Small bug: trying to plot an empty histogram gives an error

# # The created object is automatically drawn because
# # IPython matplotlib integration is already there
# import rootpy.plotting.root2matplotlib as rplt
# rplt.hist(hist);

from matplotlib import pyplot as plt
fg = plt.figure('test')
axes = plt.axes()

#plt.show()

import rootpy.plotting.root2matplotlib as rplt
from rootpy.plotting import Hist
hist = Hist(100, 0, 10, name='test hist')
hist.Fill(3) # Small bug: trying to plot an empty histogram gives an error
print hist.GetName()
#rplt.hist(hist)
rplt.hist(hist, axes=axes)
#plt.show()

import os
# print os.environ

from rootpy.io import root_open
anadir =  os.environ['E1F_2PI_ANADIR2']
print anadir
rfilename = os.path.join(anadir,'2:1__1-2.000-2.400__24-1.300-1.900__pol__exp.root')#test.root')
rfile = root_open(rfilename)
keys = rfile.GetListOfKeys()
#print keys
print keys.GetSize()
h = []
for x in keys:
	#print x.GetName()
	h.append(rfile.Get('%s/hY5D_POS/Varset1/hY5D_ACC_CORR'%
			x.GetName()));
	#print h.GetNbins()

from ROOT import gSystem
gSystem.Load('myTHnTool_C')
from ROOT import myTHnTool
_mythnt = myTHnTool()
print _mythnt.GetNbinsNotEq0(h[0])



