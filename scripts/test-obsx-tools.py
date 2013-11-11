#python imports
import os
from array import array
import math
#scipy imports
from matplotlib import pyplot as plt
import numpy as np
#rootpy imports
from rootpy.io import root_open
import rootpy.plotting.root2matplotlib as rplt
from rootpy.plotting import Hist

from ROOT import gSystem, THnSparseF
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


"""Test Case for plotR2"""
var = {'M1':0, 'M2':1, 'THETA':2, 'PHI':3, 'ALPHA':4};
print var


#INPUT data
anadir =  os.environ['E1F_2PI_ANADIR2']
rfilename = os.path.join(anadir,'1:2:3:4__1-2.000-2.400__24-1.300-1.900__pol__exp.root')#test.root')
rfile = root_open(rfilename)
keys = rfile.GetListOfKeys()
#OUTPUT data
outdir_root = os.path.join(anadir,'polobs.new')
if not os.path.isdir(outdir_root):
	os.makedirs(outdir_root)

h5 = {}
for q2wdir in keys:
	outdir = os.path.join( ('%s/%s/')%
		                   (outdir_root,q2wdir.GetName()) )
	if not os.path.isdir(outdir):
		os.makedirs(outdir);

	h5[('EC','POS')]=rfile.Get('%s/hY5D_POS/Varset1/hY5D_ACC_CORR'%
							q2wdir.GetName());
	h5[('EC','NEG')]=rfile.Get('%s/hY5D_NEG/Varset1/hY5D_ACC_CORR'%
							q2wdir.GetName());

	nbins=h5[('EC','POS')].GetNbins()
	axphi = h5[('EC','POS')].GetAxis(var['PHI'])
	for ibin in (0,nbins):
		bincoord = array('i',[0, 0, 0, 0, 0])
		binc = h5[('EC','POS')].GetBinContent(ibin,bincoord)
		binerr = h5[('EC','POS')].GetBinError(ibin)
		phi = math.radians( axphi.GetBinLowEdge(bincoord[var['PHI']]) )
		#print binc
		h5[('EC','POS')].SetBinContent(ibin,math.sin(phi)*binc)

	#h.append(rfile.Get('%s/hY5D_POS/Varset1/hY5D_ACC_CORR'%
			#q2wdir.GetName()));

#print mythnt.GetNbinsNotEq0(h[0])



