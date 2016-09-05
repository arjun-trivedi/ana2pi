#!/usr/bin/python
from __future__ import division
from collections import OrderedDict
import math
import os,sys
import matplotlib.pyplot as plt
import numpy as np


def get_err(line):
	'''
	+ For a given line in text results of Obs, 1D,itg or R2, obtain  rel_err,sg_rel_err (called err,err_err here)
	+ line format:
		bin,bin_le,bin_ue = (mu +/- sg_mu),(sg +/- sg_sg),(rel_err +/- sg_rel_err)
	'''
	#print line
	#print line.split("=")
	#print line.split("=")[1].split(',')
	#print line.split("=")[1].split(',')[2]
	rel_err=line.split("=")[1].split(',')[2]
	#! first get err
	b=rel_err.find("(")+1
	e=rel_err.find("+")
	err= float(rel_err[b:e])
	#! Now get err_err
	b=rel_err.find("-")+1
	e=rel_err.find(")")
	err_err= float(rel_err[b:e])
	#print "err,err_err=",err,err_err,err-err_err
	return err,err_err

def plot_1D(SE_1D,OUTDIR):	
	fig,ax=plt.subplots(figsize=(12,8),nrows=1,ncols=1)
	vals=[SE_1D[k][0] for k in SE_1D.keys()]
	mnm=min(vals)
	if mnm<0: mnm=int(round(min(vals))-1)
	else: mnm=0
	mxm=int(round(max(vals))+1)
	ax.set_title("1D %f"%np.mean(vals))
	ax.hist(vals,bins=(range(mnm,mxm+1,1)))
	
	outdir="%s/Obs_1D"%OUTDIR
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	fig.savefig("%s/1D.png"%(outdir))
	fig.savefig("%s/1D.eps"%(outdir))

def plot_itg(SE_ITG,OUTDIR):	
	fig,ax=plt.subplots(figsize=(12,8),nrows=1,ncols=1)
	vals=[SE_ITG[k][0] for k in SE_ITG.keys()]
	mnm=min(vals)
	if mnm<0: mnm=int(round(min(vals))-1)
	else: mnm=0
	mxm=int(round(max(vals))+1)
	ax.set_title("1D %f"%np.mean(vals))
	ax.hist(vals,bins=(range(mnm,mxm+1,1)))

	outdir="%s/Obs_itg"%OUTDIR
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	fig.savefig("%s/itg.png"%(outdir))
	fig.savefig("%s/itg.eps"%(outdir))
	
def plot_R2(SE_R2,OUTDIR,R2=None):
	'''
	Uses split x-axis since range of SE for R2 is large
	+ Split x-axis solution taken from:
	  http://stackoverflow.com/questions/5656798/python-matplotlib-is-there-a-way-to-make-a-discontinuous-axis/5669301#5669301
	'''

	#! set up axes: 2*3
	#! + for err:     ax[0][0],ax[0][1],ax[0][2] = xaxis split into lower,middle,upper 
	#! + for err_err: ax[1][0],ax[1][1],ax[1][2] (xaxis not yet spliced)
	#fig,(ax1,ax2,ax3,ax4,ax5,ax6) = plt.subplots(figsize=(20,8),nrows=2,ncols=3, sharey=True)
	fig,axs = plt.subplots(figsize=(20,15),nrows=2,ncols=3, sharey=True)
	if R2==None:
		errl    =[SE_R2[k][0] for k in SE_R2.keys()]
		err_errl=[SE_R2[k][1] for k in SE_R2.keys()]
	else:
		errl    =[SE_R2[k][0] for k in SE_R2.keys() if k[0]==R2]
		err_errl=[SE_R2[k][1] for k in SE_R2.keys() if k[0]==R2]
	#! keep track of all number of data points
	nt=len(errl)

	LMIN,LMAX,LBINW=0,10,1
	MMIN,MMAX,MBINW=10,100,10
	UMIN,UMAX,UBINW=100,4000,50
	# print LMIN,LMAX,LBINW
	# print MMIN,MMAX,MBINW
	# print UMIN,UMAX,UBINW

	# plot the same data on all axes
	#! lower limit
	#! err
	nl=sum(1 for x in errl if x>=LMIN and x<=LMAX)
	pl=(nl/nt)*100
	axs[0][0].set_title("%f %%"%(pl))
	axs[0][0].hist(errl,bins=(range(LMIN,LMAX+LBINW,LBINW)))
	#! err_errl
	#! Get appropriate data from err_errl
	d=[]
	for x in zip(errl,err_errl):
		err,err_err=x[0],x[1]
		if err>=LMIN and err<=LMAX:
			d.append(err_err)
	axs[1][0].hist(d)


	#! middle limit
	nm=sum(1 for x in errl if x>MMIN and x<=MMAX)
	pm=(nm/nt)*100
	axs[0][1].set_title("%f %%"%(pm))
	axs[0][1].hist(errl,bins=(range(MMIN,MMAX+MBINW,MBINW)))
	#! err_errl
	#! Get appropriate data from err_errl
	d=[]
	for x in zip(errl,err_errl):
		err,err_err=x[0],x[1]
		if err>MMIN and err<=MMAX:
			d.append(err_err)
	axs[1][1].hist(d)

	#! upper limits
	nh=sum(1 for x in errl if x>UMIN and x<=UMAX)
	ph=(nh/nt)*100
	axs[0][2].set_title("%f %%"%(ph))
	axs[0][2].hist(errl,bins=(range(UMIN,UMAX+UBINW,UBINW)))
	#! err_errl
	#! Get appropriate data from err_errl
	d=[]
	for x in zip(errl,err_errl):
		err,err_err=x[0],x[1]
		if err>=LMIN and err<=LMAX:
			d.append(err_err)
	axs[1][2].hist(errl)

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

	# Make the spacing between the two axes a bit smaller
	plt.subplots_adjust(wspace=0.15)

	# #! I think some entries are being lost at edges to comparison-precision
	# print nl,nm,nh,nt
	# print "diff=",(nl+nm+nl)-nt
	# print sum(1 for x in errl if x>UMAX)
	# print sum(1 for x in errl if x<LMIN)

	outdir="%s/Obs_R2"%OUTDIR
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	if R2==None:
		fig.savefig("%s/all.png"%outdir)
		fig.savefig("%s/all.eps"%outdir)
	else:
		fig.savefig("%s/%s.png"%(outdir,R2))
		fig.savefig("%s/%s.eps"%(outdir,R2))

SE_1D,SE_ITG,SE_R2=OrderedDict(),OrderedDict(),OrderedDict()
OBSDIR='/data/trivedia/e16/2pi/d2pi/SS/dbg/lowQ2_cmb_non_vst_SE_090516'

#! Fill SE_1D
print "*** Filling SE_1D ***"
rootDir = '%s/Obs_1D_txt_rslt'%(OBSDIR)
for dirName, subdirList, fileList in os.walk(rootDir):
	#! Avoid root path
	if dirName==rootDir:continue

	#! Avoid "EC"
	if "EC" in dirName: continue

	# print('Found directory: %s' % dirName)
	# for fname in fileList:
	# 	print('\t%s' % fname)

	if len(fileList)==9: 
		#! obtain q2,w info from dirName
		q2=dirName.split("/")[10]
		w=dirName.split("/")[11]
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
print "******"
#sys.exit()

#! Fill SE_ITG
print "*** Filling SE_ITG ***"
for vst in [1,2,3]:
	rootDir = '%s/Obs_itg_txt_rslt/VST%d'%(OBSDIR,vst)
	for dirName, subdirList, fileList in os.walk(rootDir):
		#! Avoid root path
		if dirName==rootDir:continue
	
		#print('Found directory: %s' % dirName)
		#for fname in fileList:
		#	print('\t%s' % fname)
		f=open("%s/%s"%(dirName,fileList[0]),"r")
		print "%s/%s"%(dirName,fileList[0])
		#print f.name
		#sys.exit()
		q2w=dirName.split("/")[11]
		#print "q2w=",q2w
		print "Reading data from", f.name
		fdata=f.readlines()
		for ibin,line in enumerate(fdata):
			err,err_err=get_err(line)
			SE_ITG[q2w,vst,ibin+1]=[err,err_err]
print "******"

#! Fill SE_R2
print "*** Filling SE_R2 ***"
rootDir = '%s/Obs_R2_txt_rslt'%(OBSDIR)
for dirName, subdirList, fileList in os.walk(rootDir):
	#! Avoid root path
	if dirName==rootDir:continue

	#! Avoid "EC"
	if "EC" in dirName: continue

	# print('Found directory: %s' % dirName)
	# for fname in fileList:
	# 	print('\t%s' % fname)

	if len(fileList)==12:
		#! Obtain q2,w,R2 info from dirName 
		q2=dirName.split("/")[10]
		w=dirName.split("/")[11]
		R2=dirName.split("/")[12]
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
print "******"
#sys.exit()


outdir="/tmp/figs"
if not os.path.exists(outdir):
	os.makedirs(outdir)

plot_R2(SE_R2,outdir)
for R2 in ['A','B','C','D','E']:
	plot_R2(SE_R2,outdir,R2)

plot_1D(SE_1D,outdir)
plot_itg(SE_ITG,outdir)

