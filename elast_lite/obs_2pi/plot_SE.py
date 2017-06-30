#!/usr/bin/python
from __future__ import division
from collections import OrderedDict
import math
import os,sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches


'''
[09-06-16]
+ Makes plots to visualize systematic errors
+ >plot_SE.py OBSDIR Q2BIN EXTRACT_D_E_FOR_NON_ALPHA
	+ Q2BIN is used only for aesthetic purposes, for examples, labelling plots

+ First the functions needed are defined
+ These functions are then used in the code that follows
'''

#! 1. Define functions needed
def get_hist_n(d,min,max):
	'''
	Returns total number of list of elements that will be plotted by matplotlib.hist for a 
	given list of data, d, and specified min, and max
	'''
	n=0
	dplotl=[]
	for x in d:
		if x>=min and x<=max: 
			n+=1
			dplotl.append(x)
	return n

def get_hist_plotl(d,min,max):
	'''
	Returns list of elements that will be plotted by matplotlib.hist for a 
	given list of data, d, and specified min, and max
	'''
	n=0
	dplotl=[]
	for x in d:
		if x>=min and x<=max: 
			n+=1
			dplotl.append(x)
	return dplotl

def get_err(line):
	'''
	+ For a given line in text results of Obs, 1D,itg or R2, obtain  rel_err,sg_rel_err (called err,err_err here)
	+ line format:
		bin,bin_le,bin_ue = (mu +/- sg_mu),(sg +/- sg_sg),(rel_err +/- sg_rel_err)
	'''
	rel_err=line.split("=")[1].split(',')[2]
	#! first get err
	b=rel_err.find("(")+1
	e=rel_err.find("+/-") #! Till [06-27-17] used to be: e=rel_err.find("+"). While this had no adverse effect, the fix is truer to the line format.
	err= float(rel_err[b:e]) 
	#! Now get err_err
	b=rel_err.find("+/-")+3 #!  Till [06-27-17] used to be: b=rel_err.find("-")+1. This was causing problems when number following +/- was negative.
	e=rel_err.find(")")
	err_err= float(rel_err[b:e])
	#print "err,err_err=",err,err_err,err-err_err
	return err,err_err

def plot_generic(SE,OUTDIR,obs_type,Q2BIN):
	R2=None
	OBSTYPEL=["1D","itg","R2_all","R2_A","R2_B","R2_C","R2_D","R2_E"]
	if obs_type not in OBSTYPEL: 
		print "ERROR! Please pass obs_type from",OBSTYPEL
		return
	if "R2" in obs_type:
		if "A" in obs_type: R2="A"
		if "B" in obs_type: R2="B"
		if "C" in obs_type: R2="C"
		if "D" in obs_type: R2="D"
		if "E" in obs_type: R2="E"
		if "all" in obs_type: R2="all"


	if R2==None or R2=="all":
		errl    =[SE[k][0] for k in SE.keys()]
		err_errl=[SE[k][1] for k in SE.keys()]
	else:
		errl    =[SE[k][0] for k in SE.keys() if k[0]==R2]
		err_errl=[SE[k][1] for k in SE.keys() if k[0]==R2]
	#! Calculate relative error of relative error
	#! + Note cases where err=0, err_err set to 0, which is true mathematically also
	rel_err_errl=[(x/y)*100 if y>0 else 0 for x,y in zip(err_errl,errl)]
	#! Remove 0s since rel_err=0 where err=0, and therefore is meaningless
	rel_err_errl=[x for x in rel_err_errl if x != 0]

	#! Get total number of data points
	nt=len(errl)

	#! Start making plots
	fig,axs=plt.subplots(figsize=(12,8),nrows=2,ncols=1)
	fig.subplots_adjust(hspace=.5)

	#! Plot err
	#! Determine mnm,mxm
	mnm=min(errl)
	if mnm<0: mnm=int(round(min(errl))-1)
	else: mnm=0
	mxm=int(round(max(errl))+1)
	if mxm>100:
		if obs_type=="1D" or obs_type=="itg":
			mxm=50
		else: #! for obs_type=R2
			mxm=1000 #! for obs_type=R2
	#! Get number of plotted data points
	plotl=get_hist_plotl(errl,mnm,mxm)
	n=len(plotl)
	f=(n/nt)*100
	mean=np.mean(plotl)
	axs[0].set_title("Estimated systematic error for Obs_%s(%s) \n average=%.2f %% (n=%d, Fraction of total points=%.2f %%)"%(obs_type,Q2BIN,mean,n,f))
	axs[0].hist(plotl,bins=(range(mnm,mxm+1,1)))

	#! Plot rel_err_err
	#! Determine mnm,mxm
	mnm=min(rel_err_errl)
	if mnm<0: mnm=int(round(min(rel_err_errl))-1)
	else: mnm=0
	mxm=int(round(max(rel_err_errl))+1)
	if mxm>100: mxm=500
	#! Get number of plotted data points
	plotl=get_hist_plotl(rel_err_errl,mnm,mxm)
	n=len(plotl)
	f=(n/nt)*100
	mean=np.mean(plotl)
	axs[1].set_title("Error on estimated systematic error for Obs_%s(%s) \n average=%.2f %% (n=%d, Fraction of total points=%.2f %%) "%(obs_type,Q2BIN,mean,n,f))
	axs[1].hist(plotl,bins=(range(mnm,mxm+1,10)))

	if R2==None:
		outdir="%s/Obs_%s"%(OUTDIR,obs_type)
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		fig.savefig("%s/%s.png"%(outdir,obs_type))
		fig.savefig("%s/%s.eps"%(outdir,obs_type))
	else:
		outdir="%s/Obs_R2"%(OUTDIR)
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		fig.savefig("%s/%s.png"%(outdir,R2))
		fig.savefig("%s/%s.eps"%(outdir,R2))


def plot_R2_zoom(SE_R2,OUTDIR,R2,Q2BIN):
	'''
	Uses split x-axis since range of SE for R2 is large
	+ Split x-axis solution taken from:
	  http://stackoverflow.com/questions/5656798/python-matplotlib-is-there-a-way-to-make-a-discontinuous-axis/5669301#5669301
	'''
	R2L=["all","A","B","C","D","E"]
	if R2 not in R2L: 
		print "Please enter R2 from:",R2L
		return
	if R2=="all":
		errl    =[SE_R2[k][0] for k in SE_R2.keys()]
		err_errl=[SE_R2[k][1] for k in SE_R2.keys()]
	else:
		errl    =[SE_R2[k][0] for k in SE_R2.keys() if k[0]==R2]
		err_errl=[SE_R2[k][1] for k in SE_R2.keys() if k[0]==R2]
	#! Calculate relative error of relative error
	#! + Note cases where err=0, err_err set to 0, which is true mathematically also
	rel_err_errl=[(x/y)*100 if y>0 else 0 for x,y in zip(err_errl,errl)]
	#! Remove 0s since rel_err=0 where err=0, and therefore is meaningless
	rel_err_errl=[x for x in rel_err_errl if x != 0]

	#! keep track of all number of data points
	nt_err=len(errl)
	nt_rel_err=len(rel_err_errl)

	#! set up axes: 2*3
	#! + for err:     ax[0][0],ax[0][1],ax[0][2] = xaxis split into lower,middle,upper 
	#! + for err_err: ax[1][0],ax[1][1],ax[1][2] (xaxis not yet spliced)
	#fig,(ax1,ax2,ax3,ax4,ax5,ax6) = plt.subplots(figsize=(20,8),nrows=2,ncols=3, sharey=True)
	fig,axs = plt.subplots(figsize=(20,15),nrows=2,ncols=3, sharey=True)
	if R2=="all":
		figtitle="R2_all"
	else:
		figtitle="R2_%s"%R2
	fig.suptitle("Estimated systematic error for Obs_%s(%s) (average=%.2f %%) (top row) \n Error on estimated systematic error (average=%.2E %%) (bottom row) \n (These averages are for all points some of which may beyond the upper limits of the plots below)"%(figtitle,Q2BIN,np.mean(errl),np.mean(rel_err_errl)) ,fontsize='xx-large')
	fig.subplots_adjust(top=0.82)

	LMIN,LMAX,LBINW=0,10,1
	MMIN,MMAX,MBINW=10,100,10
	UMIN,UMAX,UBINW=100,4000,50
	# print LMIN,LMAX,LBINW
	# print MMIN,MMAX,MBINW
	# print UMIN,UMAX,UBINW

	# plot the same data on all axes
	#! lower limit
	#! err
	ax=axs[0][0]
	plotl=get_hist_plotl(errl,LMIN,LMAX)
	n=len(plotl)
	f=(n/nt_err)*100
	mean=np.mean(plotl)
	ax.set_title("Average=%.2f %% \n n= %d (Fraction of data points %.2f %%)"%(mean,n,f))
	ax.hist(plotl,bins=(range(LMIN,LMAX+LBINW,LBINW)))
	#! rel_err_errl
	ax=axs[1][0]
	#! Get appropriate data from rel_err_errl
	d=[]
	for x in zip(errl,rel_err_errl):
		err,rel_err_err=x[0],x[1]
		if err>=LMIN and err<=LMAX:
			d.append(rel_err_err)
	mnm=min(d)
	if mnm<0: mnm=int(round(min(d))-1)
	else: mnm=0
	mxm=int(max(d)+1)
	if mxm>100:mxm=500
	plotl=get_hist_plotl(d,mnm,mxm)
	n=len(plotl)
	f=(n/nt_err)*100
	mean=np.mean(plotl)
	ax.set_title("Average=%.2f %% \n n= %d (Fraction of data points %.2f %%)"%(mean,n,f))
	ax.hist(plotl,bins=(range(mnm,mxm+1,10)))
	
	#! middle limit
	#! err
	ax=axs[0][1]
	plotl=get_hist_plotl(errl,MMIN,MMAX)
	n=len(plotl)
	f=(n/nt_err)*100
	mean=np.mean(plotl)
	ax.set_title("Average=%.2f %% \n n= %d (Fraction of data points %.2f %%)"%(mean,n,f))
	ax.hist(plotl,bins=(range(MMIN,MMAX+MBINW,MBINW)))
	#! err_errl
	ax=axs[1][1]
	#! Get appropriate data from rel_err_errl
	d=[]
	for x in zip(errl,rel_err_errl):
		err,rel_err_err=x[0],x[1]
		if err>MMIN and err<=MMAX:
			d.append(rel_err_err)
	mnm=min(d)
	if mnm<0: mnm=int(round(min(d))-1)
	else: mnm=0
	mxm=int(max(d)+1)
	if mxm>100:mxm=500
	plotl=get_hist_plotl(d,mnm,mxm)
	n=len(plotl)
	f=(n/nt_err)*100
	mean=np.mean(plotl)
	ax.set_title("Average=%.2f %% \n n= %d (Fraction of data points %.2f %%)"%(mean,n,f))
	ax.hist(plotl,bins=(range(mnm,mxm+1,10)))
		
	#! upper limits
	#! err
	ax=axs[0][2]
	plotl=get_hist_plotl(errl,UMIN,UMAX)
	n=len(plotl)
	f=(n/nt_err)*100
	mean=np.mean(plotl)
	ax.set_title("Average=%.2f %% \n n= %d (Fraction of data points %.2f %%)"%(mean,n,f))
	ax.hist(plotl,bins=(range(UMIN,UMAX+UBINW,UBINW)))
	#! err_errl
	ax=axs[1][2]
	#! Get appropriate data from rel_err_errl
	d=[]
	for x in zip(errl,rel_err_errl):
		err,rel_err_err=x[0],x[1]
		if err>=UMIN and err<=UMAX:
			d.append(rel_err_err)
	mnm=min(d)
	if mnm<0: mnm=int(round(min(d))-1)
	else: mnm=0
	mxm=int(max(d)+1)
	if mxm>100:mxm=500
	plotl=get_hist_plotl(d,mnm,mxm)
	n=len(plotl)
	f=(n/nt_err)*100
	mean=np.mean(plotl)
	ax.set_title("Aaverage=%.2f %% \n n= %d (Fraction of data points %.2f %%)"%(mean,n,f))
	ax.hist(plotl,bins=(range(mnm,mxm+1,10)))
			
	# zoom-in / limit the view to different portions of the data
	axs[0][0].set_xlim(LMIN,LMAX+LBINW) # most of the data
	axs[0][1].set_xlim(MMIN,MMAX+MBINW) # outliers only
	axs[0][2].set_xlim(UMIN,UMAX+UBINW)

	# hide the spines between ax and ax2
	axs[0][0].spines['right'].set_visible(False)
	axs[0][1].spines['left'].set_visible(False)
	axs[0][1].spines['right'].set_visible(False)
	axs[0][2].spines['left'].set_visible(False)
	axs[0][0].yaxis.tick_left()
	axs[0][0].tick_params(labeltop='off') # don't put tick labels at the top
	axs[0][2].yaxis.tick_right()

	axs[1][2].yaxis.tick_right()

	# Make the spacing between the two axes a bit smaller
	plt.subplots_adjust(wspace=0.15)

	# #! I think some entries are being lost at edges to comparison-precision
	# print nl,nm,nh,nt_err
	# print "diff=",(nl+nm+nl)-nt
	# print sum(1 for x in errl if x>UMAX)
	# print sum(1 for x in errl if x<LMIN)

	outdir="%s/Obs_R2_zoom"%OUTDIR
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	if R2=="all":
		fig.savefig("%s/all.png"%outdir)
		fig.savefig("%s/all.eps"%outdir)
	else:
		fig.savefig("%s/%s.png"%(outdir,R2))
		fig.savefig("%s/%s.eps"%(outdir,R2))

def fill_and_get_SE_1D(OBSDIR):
	print "*** In  fill_and_get_SE_1D ()***"
	
	SE_1D=OrderedDict()

	rootDir = '%s/Obs_1D_txt_rslt'%(OBSDIR)
	for dirName, subdirList, fileList in os.walk(rootDir):
		#! Avoid root path
		if dirName==rootDir:continue

		#! Avoid "EC"
		if "EC" in dirName: continue

		#! Use following to debug
		# print('Found directory: %s' % dirName)
		# for fname in fileList:
		# 	print('\t%s' % fname)

		if len(fileList)==9: 
			#! obtain q2,w info from dirName
			q2=dirName.split("/")[-2]#[10]
			w=dirName.split("/")[-1]#[11]
			q2=q2.replace("q","")
			w=w.replace("w","")
			#print "q2,w=",q2,w

			for fname in fileList:
				#! obtain vst,var information from file
				#print fname
				vst=int(fname.split(".txt")[0].split("_")[0])
				var=fname.split(".txt")[0].split("_")[1]
				#print "vst,var=",vst,var
				#! now read file
				f=open("%s/%s"%(dirName,fname),"r")
				print "Reading data from", f.name
				fdata=f.readlines()
				for ibin,line in enumerate(fdata):
					err,err_err=get_err(line)
					SE_1D[q2,w,vst,var,ibin+1]=[err,err_err]
	return SE_1D
	print "*** Done  fill_and_get_SE_1D ()***"


def fill_and_get_SE_ITG(OBSDIR,SE_TYPE):
	print "*** In  fill_and_get_SE_ITG() ***"

	SE_ITG=OrderedDict()

	if SE_TYPE=="cmb_vst_SE": VST=[""]
	else:                     VST=[1,2,3]
	for vst in [1,2,3]:
		if SE_TYPE=="cmb_vst_SE": rootDir = '%s/Obs_itg_txt_rslt'%(OBSDIR)
		else:                     rootDir = '%s/Obs_itg_txt_rslt/VST%d'%(OBSDIR,vst)
		for dirName, subdirList, fileList in os.walk(rootDir):
			#! Avoid root path
			if dirName==rootDir:continue

			#! Avoid "EC"
			if "EC" in dirName: continue
	
			#! Use following to debug
			print('Found directory: %s' % dirName)
			for fname in fileList:
				print('\t%s' % fname)

			if len(fileList)==1 and fileList[0]=="result.txt":
				f=open("%s/%s"%(dirName,fileList[0]),"r")
				print "%s/%s"%(dirName,fileList[0])
				#print f.name
				#sys.exit()
				q2w=dirName.split("/")[-1]#[11]
				#print "q2w=",q2w
				print "Reading data from", f.name
				fdata=f.readlines()
				for ibin,line in enumerate(fdata):
					err,err_err=get_err(line)
					SE_ITG[q2w,vst,ibin+1]=[err,err_err]
	return SE_ITG
	print "*** Dobe fill_and_get_SE_ITG() ***"

def fill_and_get_SE_R2(OBSDIR_R2,EXTRACT_D_E_FOR_NON_ALPHA):
	print "*** In  fill_and_get_SE_R2() ***"

	SE_R2=OrderedDict()

	rootDir = '%s/Obs_R2_txt_rslt'%(OBSDIR_R2)
	for dirName, subdirList, fileList in os.walk(rootDir):
		#! Avoid root path
		if dirName==rootDir:continue

		#! Avoid "EC"
		if "EC" in dirName: continue

		# print('Found directory: %s' % dirName)
		# for fname in fileList:
		# 	print('\t%s' % fname)

		if len(fileList)==12 or (EXTRACT_D_E_FOR_NON_ALPHA==False and len(fileList)==3):
			#! Obtain q2,w,R2 info from dirName 
			q2=dirName.split("/")[-3]#[10]
			w=dirName.split("/")[-2]#[11]
			R2=dirName.split("/")[-1]#[12]
			q2=q2.replace("q","")
			w=w.replace("w","")
			#print "q2,w,R2=",q2,w,R2

			for fname in fileList:
				#! obtain vst,var information from file
				#print fname
				vst=int(fname.split(".txt")[0].split("_")[0])
				var=fname.split(".txt")[0].split("_")[1]
				#print "vst,var=",vst,var
				#! now read file
				f=open("%s/%s"%(dirName,fname),"r")
				print "Reading data from", f.name
				fdata=f.readlines()
				for ibin,line in enumerate(fdata):
					err,err_err=get_err(line)
					SE_R2[R2,q2,w,vst,var,ibin+1]=[err,err_err]
	return SE_R2
	print "*** Done fill_and_get_SE_R2() ***"
#sys.exit()

#! 2. Use functions to make plots

#! Get inputs from user
if len(sys.argv)<2:
	sys.exit("Please enter OBSDIR")
else:
	OBSDIR=sys.argv[1]
print "OBSDIR=",OBSDIR

if len(sys.argv)<3:
	sys.exit("Please enter Q2BIN")
else:
	Q2BIN=sys.argv[2]
print "OBSDIR=",OBSDIR

if len(sys.argv)<4:
	sys.exit("Please enter EXTRACT_D_E_FOR_NON_ALPHA")
else:
	if   sys.argv[3]=="True":  EXTRACT_D_E_FOR_NON_ALPHA=True
	elif sys.argv[3]=="False": EXTRACT_D_E_FOR_NON_ALPHA=False
print "EXTRACT_D_E_FOR_NON_ALPHA=",EXTRACT_D_E_FOR_NON_ALPHA


#! Determine SE_TYPE
SE_TYPE=""
#print os.listdir('%s/Obs_itg_txt_rslt'%OBSDIR)
if not "VST" in os.listdir('%s/Obs_itg_txt_rslt'%OBSDIR):
	SE_TYPE="cmb_vst_SE"
print "SE_TYPE=",SE_TYPE
#sys.exit()

#! Define needed objects
#! OBSDIR,OUTDIR for 1D and ITG
OUTDIR="%s/SE_plots"%OBSDIR
if not os.path.exists(OUTDIR):
	os.makedirs(OUTDIR)
print "OBSDIR=",OBSDIR
print "OUTDIR=",OUTDIR
#! OBSDIR,OUTDIR for R2
if   EXTRACT_D_E_FOR_NON_ALPHA==True:    OBSDIR_R2='%s/%s'%(OBSDIR,"Obs_R2_w_non-alpha_DE")
elif EXTRACT_D_E_FOR_NON_ALPHA==False:   OBSDIR_R2=OBSDIR
OUTDIR_R2="%s/SE_plots"%(OBSDIR_R2)
print "OBSDIR_R2=",OBSDIR_R2
print "OUTDIR_R2=",OUTDIR_R2
#sys.exit()


SE_1D=fill_and_get_SE_1D(OBSDIR)
plot_generic(SE_1D,OUTDIR,"1D",Q2BIN)

SE_ITG=fill_and_get_SE_ITG(OBSDIR,SE_TYPE)
plot_generic(SE_ITG,OUTDIR,"itg",Q2BIN)

SE_R2=fill_and_get_SE_R2(OBSDIR_R2,EXTRACT_D_E_FOR_NON_ALPHA)
plot_generic(SE_R2,OUTDIR_R2,"R2_all",Q2BIN)
plot_R2_zoom(SE_R2,OUTDIR_R2,"all",Q2BIN)
for R2 in ['A','B','C','D','E']:
	plot_generic(SE_R2,OUTDIR_R2,"R2_%s"%(R2),Q2BIN)
	plot_R2_zoom(SE_R2,OUTDIR_R2,R2,Q2BIN)

print "plot_SE.py done"
