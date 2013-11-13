# # Testing dict
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

# # Test for loop
# for i in xrange(0,10):
# 	print i


# #Test passing array to ROOT
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

#Test 2 by 2 DS
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

#Test walking directory structure
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

ldirs = [dirs[0] for dirs in os.walk(outdir)]
print ldirs
	
#Test Subprocess
import subprocess
out=open('/tmp/stdout.text','w')
err=open('/tmp/stderr.text','w')
code = subprocess.call(["ls", "test"],stdout=out,stderr=err)
print 'code=',code
if code!=0:
	print 'failed!'
	
