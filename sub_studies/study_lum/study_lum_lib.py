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

def calc_lum(Q):
    	'''
	+ Q=Total Faraday Cup charge in Coulombs
    	+ Returns Luminosity in fb^-1 (inverse femtobarns)
    	'''
	#! + The source for Luminosity calculation is Kijun Park's thesis, page 139
        #! + If I enter Q as obtained by Kijun(=2128.7181 microC), L returned ~ 28 fb^-1 (=E16 Lum)
    	rho=0.0708 #[g/cm^3] Liquid Hydrogen target density at 20K
    	Na=6.022136e23 #[mol^-1] Avogadro's number
    	Ltg=5 #[cm] E1F target length
    	amu=1.00795 #[g/mol] Atomic Mass Unit
    	e=1.602176e-19 #[C] electron charge
    
    	L=(Q/e)*((Na*rho*Ltg)/amu) #[cm^-2]
    	return L*10e-39 #[fb^-1]
