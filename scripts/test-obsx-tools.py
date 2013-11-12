#python imports
import os
from array import array
import math
#scipy imports
from matplotlib import pyplot as plt
import numpy as np
#rootpy imports
from rootpy.io import root_open, File
import rootpy.plotting.root2matplotlib as rplt
#import root2matplot as r2m
from rootpy.plotting import Hist

from ROOT import gSystem, THnSparseF
gSystem.Load('myTHnTool_C')
from ROOT import myTHnTool
mythnt = myTHnTool()

#an example of plotting with rootpy and pyplot
# fg = plt.figure('test')
# axes = plt.axes()
# hist = Hist(100, 0, 10, name='test hist')
# hist.Fill(3) # Small bug: trying to plot an empty histogram gives an error
# print hist.GetName()
# #rplt.hist(hist)
# rplt.hist(hist, axes=axes)
# #plt.show()


"""Test Case for plotR2"""
var = {'M1':0, 'M2':1, 'THETA':2, 'PHI':3, 'ALPHA':4};
print var


#INPUT data
anadir =  os.environ['E1F_2PI_ANADIR2']
fname = os.path.join(anadir,'1:2:3:4__1-2.000-2.400__24-1.300-1.900__pol__exp.root')#test.root')
f = root_open(fname)
#f = File(fname,'read')
# for path, dirs, objects in rfile.walk():
# 	print path,dirs#,objects
tl = list(f)
print tl
for root,subfolders,files in f.walk():
	print 'root=', root
	print 'subfolders=',subfolders
	print 'files=',files
	# for iroot in root:
	# 	obj = os.path.join(iroot,'hY5D_POS','Varset1','hY5D_ACC_CORR')
	# 	print 'obj =', obj
# for i in t:
# 	print i
#print f.keys()
# for dirs,paths,objects in f.walk('',maxdepth=0):
# 	for d in dirs:
# 		print d
	#h5[('EC','POS')] = rfile.q2wdir.hY5D_POS.Varset1.hY5D_ACC_CORR
	# print dirs
	# for objects in dirs:
	# 	print objects

# keys = rfile.GetListOfKeys()
# #OUTPUT data
# outdir_root = os.path.join(anadir,'polobs.new')
# if not os.path.isdir(outdir_root):
# 	os.makedirs(outdir_root)

# for q2wdir in keys:
# 	outdir = os.path.join( ('%s/%s/')%(outdir_root,q2wdir.GetName()) )
# 	if not os.path.isdir(outdir):
# 		os.makedirs(outdir);

# 	h5 = {}	
# 	h5[('EC','POS')]=rfile.Get('%s/hY5D_POS/Varset1/hY5D_ACC_CORR'%
# 							q2wdir.GetName());
# 	h5[('EC','NEG')]=rfile.Get('%s/hY5D_NEG/Varset1/hY5D_ACC_CORR'%
# 							q2wdir.GetName());
# 	mythnt.MultiplyBySinPhi(h5[('EC','POS')]);
# 	mythnt.MultiplyBySinPhi(h5[('EC','NEG')],-1);
# 	hR2 = {}
# 	fig = plt.figure('test')
# 	ax = fig.add_subplot(111)
# 	hR2[('THETA','POS')] = h5[('EC','POS')].Projection(var['THETA'])
# 	#rplt.hist(hR2[('THETA','POS')],axes=ax)
# 	print hR2[('THETA','POS')].GetName()

# 	fig.savefig("test.png")
# 	fig.show()
# 	#rplt.hist(hR2[('THETA','POS')],axes=ax)

	





