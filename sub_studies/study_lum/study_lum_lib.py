from __future__ import division
import os,glob#,subprocess
import ROOT

import matplotlib.pyplot as plt

def make_dlum_txt(dlumdir,flumname,runl=None):
	'''
	+ dlumdir=directory where dlum is located
	+ flumname=name of output text file that contains lum data
	+ runl=optional; if provided then dlum_txt only made from runl provided
	'''
	lumdl=glob.glob(os.path.join(dlumdir,"*.root"))
	lumdl.sort()
        fout=open(flumname, 'w')
	#print lumdl
	for lumd in lumdl:
		runnum=int(lumd.split("/")[-1].split(".")[0])
		#print runnum
		if runl is not None:
			if runnum not in runl: continue
		f=ROOT.TFile(lumd)
		#print f.GetName()
		hevtsum=f.Get("lum/hevtsum")
		hQ=f.Get("lum/hQ")
		#print hevtsum.GetName()
		q=hQ.GetBinContent(1)
		Ntrigger=hevtsum.GetBinContent(1)
		Ngoode=hevtsum.GetBinContent(2)
		Ngoodebos=hevtsum.GetBinContent(3)
		fout.write("%d %.2f %d %d %d\n"%(runnum,q,Ntrigger,Ngoode,Ngoodebos))
	fout.close()

def disp_lum(flumname):
	'''
	+ flumname=name of output text file that contains lum data
	+ Returns Q(=list of Q per Run)
        '''
	ff=open(flumname)
	runnuml=[]
	Q=[]
	norm_Ntrigger=[]
	norm_Ngoode=[]
	norm_Ngoodebos=[]
	for l in ff.readlines():
		#print l
    		dl=l.split(" ")
    		runnum=int(dl[0].strip())
    		q=float(dl[1].strip())
    		Ntrigger=int(dl[2].strip())
    		Ngoode=int(dl[3].strip())
    		Ngoodebos=int(dl[4].strip())
    		#print runnum,q,Ntrigger,Ngoode
    		runnuml.append(runnum)
		Q.append(q)
    		if q==0.0:
        		norm_Ntrigger.append(0)
        		norm_Ngoode.append(0)
        		norm_Ngoodebos.append(0)
    		else:
        		norm_Ntrigger.append(Ntrigger/q)
        		norm_Ngoode.append(Ngoode/q)
        		norm_Ngoodebos.append(Ngoodebos/q)
	fig,axs=plt.subplots(figsize=(15,5),nrows=1,ncols=4)
	axs[0].scatter(runnuml,Q)
	axs[1].scatter(runnuml,norm_Ntrigger)
	axs[2].scatter(runnuml,norm_Ngoode)
	axs[3].scatter(runnuml,norm_Ngoodebos)
	return Q

KHETARPAL_Q=21.287*1e-3
def calc_lum(Q):
    	'''
	+ Q=Total Faraday Cup charge in Coulombs
    	+ Returns Luminosity in fb^-1 (inverse femtobarns)

        The formula for Luminosity for me is not an exact calculation. It makes sense to me if think of it
	as following:
	+ The goal is to calculate total number of possible interactions per unit area.
	Therefore, 
	1. Target: Even though it has a 5 cm length, the formula seems to be such that the incoming beam particles "see" all the target particles distributed in 2D space
	2. Beam: The beam is such that within the 2D extent of the target, all target particles present to the particles in the beam the same probability of interaction. 
	
    	'''
	#! + The formula for Luminosity calculation is as per Kijun Park's thesis, page 139
        #! + If I enter Q as obtained by Khetarpal for E16(=21.287 mC), L returned ~29 fb^-1 (=E16 Lum)
    	rho=0.073 #[g/cm^3] Liquid Hydrogen target density (Kijun=0.0708 at 20K)
    	Na=6.022e23 #[mol^-1] Avogadro's number
    	Ltg=5 #[cm] E1F target length
    	M_H=1.007 #[g/mol] Molar mass of Hydrogen 
    	e=1.602e-19 #[C] electron charge
    
    	L=(Q/e)*((Na*rho*Ltg)/M_H) #[cm^-2]
    	return L*1e-39 #[fb^-1]
