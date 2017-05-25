from __future__ import division
import os,glob#,subprocess
import ROOT

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def make_dlum_txt(dlumdir,d2pidir,flumname,runl=None,runl_ignore=None):
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
		if runl_ignore is not None:
			if runnum in runl_ignore:
				print "Ignoring run",runnum
				continue
		#! Get dlum and d2piR files for run
		flum=ROOT.TFile(lumd)
		fd2pi=ROOT.TFile('%s/%d.root'%(d2pidir,runnum))
		#print f.GetName()
		#! Get info from dlum
		hevtsum=flum.Get("lum/hevtsum")
		hQ=flum.Get("lum/hQ")
		#print hevtsum.GetName()
		q=hQ.GetBinContent(1)
		Ntrigger=hevtsum.GetBinContent(1)
		Ngoode=hevtsum.GetBinContent(2)
		Ngoodebos=hevtsum.GetBinContent(3)
		#! Get info from d2pi
		hevt=fd2pi.Get("hevt")
		N2pi=hevt.GetBinContent(9)
		fout.write("%d %.2f %d %d %d %d\n"%(runnum,q,Ntrigger,Ngoode,Ngoodebos,N2pi))
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
	norm_N2pi=[]
	#! "bad runs" and their excluded charge
	bad_runl=[]
	bad_run_ql=[]
	for l in ff.readlines():
		#print l
		dl=l.split(" ")
		runnum=int(dl[0].strip())
		q=float(dl[1].strip())
		Ntrigger=int(dl[2].strip())
		Ngoode=int(dl[3].strip())
		Ngoodebos=int(dl[4].strip())
		N2pi=int(dl[5].strip())
		#print runnum,q,Ntrigger,Ngoode
		runnuml.append(runnum)
		Q.append(q)
		if q==0.0:
			norm_Ntrigger.append(0)
			norm_Ngoode.append(0)
			norm_Ngoodebos.append(0)
			norm_N2pi.append(0)
		else:
			norm_Ntrigger.append(Ntrigger/q)
			norm_Ngoode.append(Ngoode/q)
			norm_Ngoodebos.append(Ngoodebos/q)
			norm_N2pi.append(N2pi/q)
	fig,axs=plt.subplots(figsize=(15,5),nrows=3,ncols=3)
	axs[0][0].scatter(runnuml,Q)
	axs[0][0].set_ylabel("Q")
	axs[0][1].scatter(runnuml,norm_Ntrigger)
	axs[0][1].set_ylabel("Ntrigger/Q")
	axs[0][2].scatter(runnuml,norm_Ngoode)
	axs[0][2].set_ylabel("Ntrigger/goode")
	axs[1][0].scatter(runnuml,norm_Ngoodebos)
	axs[1][0].set_ylabel("Ntrigger/goodebos")
	axs[1][1].scatter(Q,norm_Ntrigger)
	axs[1][1].set_xlabel("run#")
	axs[1][1].set_ylabel("Ntrigger/Q")
	axs[1][2].hist(norm_Ntrigger,bins=50,range=(60000,100000))
	axs[1][2].set_xlabel("Ntrigger/Q")
	axs[2][0].scatter(runnuml,norm_N2pi)
	axs[2][0].set_xlabel("N2pi/Q")
	axs[2][1].hist(norm_N2pi)#,bins=50,range=(60000,100000))
	axs[2][1].set_xlabel("N2pi/Q")

	return Q,bad_runl,bad_run_ql

def make_lum_df(flumname):
	'''
	+ flumname=name of output text file that contains lum data
	+ Returns a df containing runnum,q,Ntrigger,Ngoode,norm_Ngoodebos,N2pi
        '''
	ff=open(flumname)
	q,Ntrg,Nge,Ngeb,N2pi={},{},{},{},{}
	for l in ff.readlines():
		#print l
		dl=l.split(" ")
		runnum=int(dl[0].strip())
		q[runnum]=float(dl[1].strip())
		Ntrg[runnum]=int(dl[2].strip())
		Nge[runnum]=int(dl[3].strip())
		Ngeb[runnum]=int(dl[4].strip())
		N2pi[runnum]=int(dl[5].strip())
	#! Put all these dicts in a df
	#! + [05-227-17] 'Nge' not added currently because its values are = 0 (bug in ProcEid::ProcEid()?)
	df=pd.DataFrame({'q':q.values(),'Ntrg':Ntrg.values(),'Ngeb':Ngeb.values(),'N2pi':N2pi.values()},
		             index=q.keys())
	df=pd.DataFrame({'q':q.values(),'Ntrg':Ntrg.values(),'Ngeb':Ngeb.values(),'N2pi':N2pi.values(),
		             'run':q.keys()})#,index=q.keys())
	#! Add normalized columns to df
	df['nrm_Ntrg']=df['Ntrg']/df['q']
	df['nrm_Ngeb']=df['Ngeb']/df['q']
	df['nrm_N2pi']=df['N2pi']/df['q']

	#! Return df sorter by runnum
	#df.sort_index(inplace=True)
	df.sort('run',inplace=True)
	#! Reset index after this sorting
	df=df.reset_index(drop=True)

	#! Convert inf,Nan to 0
	df.replace([np.inf, -np.inf], np.nan,inplace=True)
	df.fillna(0,inplace=True)

	#! pd will reorder columns alphabetically. The following puts back original order
	df=df[['run','q', 'Ntrg','Ngeb','N2pi','nrm_Ntrg','nrm_Ngeb','nrm_N2pi']]
	return df

	KP_Q_E16=21.287*1e-3 #(p156) => L=28.13 fb-1 using rho=0.0708 g/cm3
	PK_Q_E16=21.287*1e-3 #(p96)  => L=28.13 fb-1 using rho=0.0708 g/cm3
	BG_Q_E16=20.600*1e-3 #(p69)  => L=27.34 fb-1 using rho=0.0708 g/cm3
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
	#! + Note various rho values:
	#!    + KP, PK (E16)=0.0708 g/cm3          : Seems to be most accurate for e16 based on http://clasweb.jlab.org/shift/e1-6/target.html, see p156 of his KP's thesis
	#!    + BG= 0.711 g/cm3 (E16)              : Seems to be on the higher side.
	#!    + EI (E16), GF (?)=0.073 g/cm3:      : Too high. May have been accurate for GF, which was then used by EI.
	rho=0.0708 #[g/cm^3] Liquid Hydrogen target density for E16 (from KP thesis p155, which refers to ref [74]: http://clasweb.jlab.org/shift/e1-6/target.html)
	Na=6.022e23 #[mol^-1] Avogadro's number
	Ltg=5 #[cm] E16 and E1F target length
	M_H=1.007 #[g/mol] Molar mass of Hydrogen 
	e=1.602e-19 #[C] electron charge
    
	L=(Q/e)*((Na*rho*Ltg)/M_H) #[cm^-2]
	return L*1e-39 #[fb^-1]
