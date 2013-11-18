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
import os
#print os.__file__
anadir =  os.environ['E1F_2PI_ANADIR2']
outdir = os.path.join(anadir,'polobs.new')

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
a HDF5 file"""
import numpy as np
randn = np.random.randn
import pandas as pd

from ROOT import TH1D

store = pd.HDFStore('store.h5')
print 'store =',store


# h1 = TH1D("test1","test1",100,0,5)
# h2 = TH1D("test2","test2",100,0,5)
# h3 = TH1D("test3","test3",100,0,5)
# hl = [h1,h2,h3]
# hl = []
# for i in range(0,8):
# 	name = 'h%d'%(i+1)
# 	hl.append(TH1D(name,name,100,0,5))
# hl[0].Draw()
# print hl[0]


	
# #index = pd.date_range('1/1/2000', periods=8)
# index = np.arange(0,8)
# # s = pd.Series(randn(5), index=['a', 'b', 'c', 'd', 'e'])
# #df = pd.DataFrame(randn(8, 3), index=index,
# #                  columns=['A', 'B', 'C'])
# df = pd.DataFrame()
# if not df:
# 	data = pd.DataFrame({'s1':hl},index=index) # Data for 1st. Column 
# 	df = df.append(data)

# df['H'] = hl
# print 'df=',df

# store['df']=df
# #store.root.wp._v_attrs.pandas_type
# store

print store['df']