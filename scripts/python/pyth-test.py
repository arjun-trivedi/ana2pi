from __future__ import division
from math import *
"""Testing dict"""
# h = {}
# h['T'] = 'Thrown'
# h['R'] = 'Reco'
# print h['T']
# print h.keys()
# print h.values()

# h2 = {}
# h2[('R','P')] = 'Reco. Pos'
# h2[('A','N')] = 'Acc Corr, Neg'
# print h2

"""Test for loop"""
# for i in xrange(0,10):
# 	print i


"""Test passing array to ROOT"""
# from ROOT import gSystem, gROOT, gStyle, gPad, TCanvas, TF1, TH1F, THnSparseF
# from array import array

# #c1 = rootnotes.default_canvas()

# hdim = 8
# bins = array('i',[10, 10, 10, 10, 10, 10, 10, 10])
# xmin = array('d',[1,  1,  1,  1,  1,  1,  1,  1])
# xmax = array('d',[11, 11, 11, 11, 11, 11, 11, 11])

# h8a = THnSparseF("h8a", "8D THnSparse", hdim, bins, xmin, xmax)
# h8a.Sumw2();
# coord1  = array('d',[1, 1, 1, 1, 1, 1, 1, 1])
# coord10 = array('d',[10,10,10,10,10,10,10,10])
# h8a.Fill(coord1)
# print h8a.GetNbins()

# bincoord = array('i',[0, 0, 0, 0, 0])
# binc = h8a.GetBinContent(0,bincoord)
# print binc

# import math
# rad = math.radians(90)
# print rad
# print math.sin(rad)

"""Test 2 by 2 DS"""
#matrix = [[1, 2, 3, 4],[5, 6, 7, 8],[9, 10, 11, 12] ]
#var = {'M1':0, 'M2':1, 'THETA':2, 'PHI':3, 'ALPHA':4};
# nVAR=5
# M1,M2,THETA,PHI,ALPHA = range(nVAR)
# for var in range(0,nVAR):
# 	print var

# varTitle = [ 
# 				["M_{p#pi^{+}}", "M_{#pi^{+}#pi^{-}}","#theta_{#pi^{-}}", "#phi_{#pi^{-}}", "#alpha_{[p^{'}#pi^{+}][p#pi^{-}]}"],
# 			 	["M_{p#pi^{+}}", "M_{#pi^{+}#pi^{-}}","#theta_{p}", "#phi_{p}", "#alpha_{[#pi^{+}#pi^{-}][p^{'}p]"],
# 			 	["M_{p#pi^{+}}", "M_{p#pi^{-}}", "#theta_{#pi^{+}}", "#phi_{#pi^{+}}", "#alpha_{[p^{'}#pi^{-}][p#pi^{+}]"] 
# 		   ]
#print varTitle[0][M1]

"""Test walking directory structure"""
# import os
# anadir =  os.environ['E1F_2PI_ANADIR2']
# outdir = os.path.join(anadir,'polobs.new')

# for root, subfolders, files in os.walk(outdir):
# 	print 'root=',root
# 	print 'subfolders=', subfolders
# 	print 'files=',files

# for dirs in os.walk(outdir).next()[1]:
# 	print 'dirs=',dirs

# ldirs = [dirs[0] for dirs in os.walk(outdir)]
# print ldirs
	
"""Test Subprocess"""
# import subprocess
# out=open('/tmp/stdout.text','w')
# err=open('/tmp/stderr.text','w')
# #code = subprocess.call(["ls", "test"],stdout=out,stderr=err)
# opt = '-lrt'
# #code = subprocess.call(["ls",opt,'> /tmp/tmp'])#,stdout=out,stderr=err)
# subprocess.Popen('ls /home/trivedia/*.d',shell=True)
# #print 'code=',code
# # if code!=0:
# # 	print 'failed!'


""" Filling Matrix Strings in Blocks"""
# #http://stackoverflow.com/questions/14639496/python-numpy-array-of-arbitrary-length-strings
# import numpy as np
# s = np.zeros((5,5),object)
# s[0:4,0:1]="arjun"
# s[0:4,1:2]="trivedi"
# print s
# print s[0,0]

# z = np.zeros((5,5),int)
# z[0:4,0:1] = 3
# print z

""" Use Case for storing Pandas DataFrame with ROOT Histograms in 
a HDF5 file
Instructions from: http://pandas.pydata.org/pandas-docs/dev/io.html#hdf5-pytables

"""
# import numpy as np
# randn = np.random.randn
# import pandas as pd

# from ROOT import TH1D

# def make_test_hdf5():
# 	""" Create a DF with just one Column of N TH1Ds,
#     indexed by [0,...,N] and store it in a HDF5 file.
#     Instructions from: http://pandas.pydata.org/pandas-docs/dev/io.html#hdf5-pytables

#     """
# 	print 'Going to create store.h5:'
# 	#1. If store.hdf5 exists, delete it
# 	if os.path.isfile('store.h5'):
# 		os.remove('store.h5')
# 	#2. Create HDF5 store & see what it contains
# 	store = pd.HDFStore('store.h5')
# 	print 'store =',store

# 	#3. Create Histograms to store
# 	NHISTS=10
# 	hl = []
# 	for i in range(0,NHISTS):
# 		name = 'h%d'%(i+1)
# 		hl.append(TH1D(name,name,100,0,5))
		
# 	#4. Insert 1st columns into df
# 	index = np.arange(0,NHISTS)
# 	df = pd.DataFrame()
# 	if not df:
# 		data = pd.DataFrame({'H1':hl},index=index) # Data for 1st. Column 
# 		df = df.append(data)

# 	#5. Add df to store
# 	store['df']=df
# 	#6. See the stored data type 
# 	print 'Type of stored data = ',store.root.df._v_attrs.pandas_type
# 	#7. See what the file finally contains
# 	print 'store=',store

# def read_test_hdf5():
# 	print 'Going to read contents of store.h5:'
# 	if os.path.isfile('store.h5'):
# 		store = pd.HDFStore('store.h5')
# 		print 'store =',store
# 		print store['df']
# 	else:
# 		print 'store.h5 does not exist!'
	

# make_test_hdf5()
# read_test_hdf5()
#print store['df']

""" Use Case for creating multi-dimensional list
"""
#y = [[[] for i in range(1)] for i in range(1)]
# y = [[[] for j in range(4)] for i in range(10)]
# print y

# mm_fitpars=[[[[]for i in range(5)] for i in range(3)],[[[]for i in range(5)] for i in range(3)]]
# print mm_fitpars

# import numpy as np
# arrnd = np.zeros((2,3,5),'i')
# print arrnd

# from array import *
# arr=array('i',arrnd[0][0])
# arr = np.subtract(arr,10)
# print arr



""" Misc. Use Cases
"""
# import ROOT as ROOT
# def gausst(v, par):
#     arg = 0;
#     if (par[2] != 0): arg = (v[0] - par[1])/par[2];
#     fitval = par[0]*(1/(sqrt(2*pi)*par[2]))*exp(-0.5*arg*arg);
#     return fitval;

# fgausst = ROOT.TF1("fgauss_test",gausst,-3,3,3);
# fgausst.SetParameters(100,0,1);
# c3=ROOT.TCanvas("c3","c3",400,400)
# fgausst.Draw()
# print fgausst.Integral(-10,10)

# import numpy as np
# arrnd = np.zeros((2.),'i')
# print arrnd
# arrnd=np.add(arrnd,100)
# print arrnd
# arrnd=np.divide(arrnd,50)
# print arrnd
# arrnd=np.multiply(arrnd,-1)
# print arrnd
# arrnd=np.add(arrnd,2)
# print arrnd

from array import *
gpppars_cmbns=array('d',range(27))
print gpppars_cmbns