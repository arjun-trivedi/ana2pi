#!/usr/bin/python
from __future__ import division
import sys,os
import math
import scipy.integrate as integrate
import itertools
import matplotlib.pyplot as plt

'''
[03-30-16]
+ This script tests the binned-integrals for orthogonal sinusoid functions that are used to extract R2s using mthd1
+ All the binned-integrals combinations of the expected moments, cos,sin,cos2,sin2, are tested. 
+ The default phi_binw=36 degrees, which the user can vary as an input to this script
+ The output of the script is in files labelled 'm1_integral_phibinw_<phi_binw>.pdf' where all of the binned-integrals per bin
  are plotted and cases where the sum of these binned-integrals deviates from the exact-integral highlighted in red.

+ Main observations is that while most of the binned-integrals results are in agreement with the exact-integral, 
  the binned integrals for cos*sin and cos2*sin2 are not equal to zero and this deviation from zero is a f(phi_binw)
+ Therefore this deviation of the binned-integral from the exact-integral causes errors in the extraction of R2_alpha
  because d2(sigma)/d(alpha)d(phi) constain sin and sin2 moments
'''
USAGE='./demo_limitation_m1_extract_R2.py [phi_binw=36]'

def get_binned_integral(f,fp):
	y=[]
	for i in phi_bin_lel:
		phi_le=math.radians(i)
		phi_ue=math.radians(i+phi_bin_w)
		#! compute binned integral
		if f=='c':
			intg=integrate.quad(lambda x: math.cos(x),phi_le,phi_ue)
		elif f=='s':
			intg=integrate.quad(lambda x: math.sin(x),phi_le,phi_ue)
		elif f=='c2':
			intg=integrate.quad(lambda x: math.cos(2*x),phi_le,phi_ue)
		elif f=='s2':
			intg=integrate.quad(lambda x: math.sin(2*x),phi_le,phi_ue)    
		#! multiply binned integral with fp
		if fp=='c':
			y.append(math.cos(phi_le)*intg[0])
		elif fp=='s':
			y.append(math.sin(phi_le)*intg[0]) 
		elif fp=='c2':
			y.append(math.cos(2*phi_le)*intg[0])
		elif fp=='s2':
			y.append(math.sin(2*phi_le)*intg[0])
		
	return y

def get_exact_integral(f,fp):
	y=[]
	for i in phi_bin_lel:
		phi_le=math.radians(i)
		phi_ue=math.radians(i+phi_bin_w)
		if f=='c' and fp=='c':
			intg=integrate.quad(lambda x: math.cos(x)*math.cos(x),phi_le,phi_ue)
		elif (f=='c' and fp=='s') or (f=='s' and fp=='c'):
			intg=integrate.quad(lambda x: math.cos(x)*math.sin(x),phi_le,phi_ue)
		elif (f=='c' and fp=='c2') or (f=='c2' and fp=='c'):
			intg=integrate.quad(lambda x: math.cos(x)*math.cos(2*x),phi_le,phi_ue)
		elif (f=='c' and fp=='s2') or (f=='s2' and fp=='c'):
			intg=integrate.quad(lambda x: math.cos(x)*math.sin(2*x),phi_le,phi_ue)
		elif f=='s' and fp=='s':
			intg=integrate.quad(lambda x: math.sin(x)*math.sin(x),phi_le,phi_ue)
		elif (f=='s' and fp=='c2') or (f=='c2' and fp=='s'):
			intg=integrate.quad(lambda x: math.sin(x)*math.cos(2*x),phi_le,phi_ue)
		elif (f=='s' and fp=='s2') or (f=='s2' and fp=='s'):
			intg=integrate.quad(lambda x: math.sin(x)*math.sin(2*x),phi_le,phi_ue)
		elif f=='c2' and fp=='c2':
			intg=integrate.quad(lambda x: math.cos(2*x)*math.cos(2*x),phi_le,phi_ue)
		elif (f=='c2' and fp=='s2') or (f=='s2' and fp=='c2'):
			intg=integrate.quad(lambda x: math.cos(2*x)*math.sin(2*x),phi_le,phi_ue)
		elif f=='s2' and fp=='s2':
			intg=integrate.quad(lambda x: math.sin(2*x)*math.sin(2*x),phi_le,phi_ue)    
		y.append(intg[0])
	return y
	

#! setup phi binning
phi_bin_w=36
if len(sys.argv)>1:
	phi_bin_w=int(sys.argv[1])#5#36
phi_bin_lel=range(0,360,phi_bin_w)

#! setup list of functions and all possible combination of integrals (f,fp) 
fl=['c','s','c2','s2']
ffpl=list(itertools.product(fl,fl))
print "All combination of integrals=",ffpl
print "Total number of combinations=",len(ffpl)

# <codecell>

NMOMENTS=len(fl) #!c,s,c2,s2
NROWS,NCOLS=NMOMENTS,NMOMENTS
plt.figure(figsize=(25,15))
for i,ffp in enumerate(ffpl):
	#! Get binned and true integral per bin
	f,fp=ffp[0],ffp[1]
	binned_itgrl=get_binned_integral(f,fp)
	true_itgrl=get_exact_integral(f,fp)
	
	#! Draw integral per bin and label them by their sum
	ax=plt.subplot(NROWS,NCOLS,i+1)
	ax.set_title('itgrl(%s,%s):true:binned=%.2f,%.2f'%(f,fp,sum(true_itgrl),sum(binned_itgrl)))
	ax.plot(phi_bin_lel,binned_itgrl,'bs',label='binned')
	ax.plot(phi_bin_lel,true_itgrl,'gs',label='true')

	#! flag plot if sum(binned_itgrl) does not match what is expected from orthogonality of sinusoid functions
	flag=False
	if f==fp and round(sum(binned_itgrl),2)==0:flag=True
	if f!=fp and round(sum(binned_itgrl),2)!=0:flag=True
	if flag:
		ax.spines['bottom'].set_color('red')
    		ax.spines['top'].set_color('red') 
    		ax.spines['right'].set_color('red')
    		ax.spines['left'].set_color('red')

	#! Legend
	if i==0:#! Draw legend for only one of the plots on the grid
		plt.legend()
pltname="%s/%s"%(os.environ['ELAST_LITE'],'obs_2pi/dev/m1_integral_phibinw_%d.pdf'%phi_bin_w)
plt.savefig(pltname)
#     #! Draw diff
#     ax2=plt.twinx()
#     diff=[x - y for x, y in zip(true_itgrl, binned_itgrl)]
#     #diff=binned_itgrl-true_itgrl
#     ax2.plot(phi_bin_lel,diff,'rs',label='diff')
#     ax2.set_ylabel('diff',color='r')
	
	   

# <codecell>


