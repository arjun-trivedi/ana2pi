import os

from matplotlib import pyplot as plt

from rootpy.io import root_open
import rootpy.plotting.root2matplotlib as rplt
from rootpy.plotting import Hist

from ROOT import gSystem
gSystem.Load('myTHnTool_C')
from ROOT import myTHnTool
mythnt = myTHnTool()

#an example of plotting with rootpy and pyplot
fg = plt.figure('test')
axes = plt.axes()
hist = Hist(100, 0, 10, name='test hist')
hist.Fill(3) # Small bug: trying to plot an empty histogram gives an error
print hist.GetName()
#rplt.hist(hist)
rplt.hist(hist, axes=axes)
#plt.show()




anadir =  os.environ['E1F_2PI_ANADIR2']
print anadir
rfilename = os.path.join(anadir,'2:1__1-2.000-2.400__24-1.300-1.900__pol__exp.root')#test.root')
rfile = root_open(rfilename)
keys = rfile.GetListOfKeys()
print keys.GetSize()
h = []
for x in keys:
	h.append(rfile.Get('%s/hY5D_POS/Varset1/hY5D_ACC_CORR'%
			x.GetName()));

print mythnt.GetNbinsNotEq0(h[0])



