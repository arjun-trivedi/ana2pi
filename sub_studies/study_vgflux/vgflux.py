from __future__ import division
import math

'''
To verify calculation of vgflux.

Lesson learnt:
1. Gleb and Ye have a slightly different value for A that affects the result
   at the 7th decimal place:
   + vglflux(A,MP precision = GF/YT)= 0.0006011789703446808
   + vglflux(A,MP precision = AT)=    0.0006010210385911745

2. I tried two different methods for getting theta' in the calculation of epsilon and t   they both give the same result for epsilon:
	i.  Original: I forgot the derivation, but it is true
	ii. from Q2=4EE'sin^2(theta'/2) (GF uses this method)
'''
PI=3.14159265358979312
E1F_E0=5.499
FSC=0.00729927007299270073# Ye,Gleb: 0.00729927007299270073 #AT:0.00729735253
A=FSC
MP=0.93827203 #Ye,Gleb: 0.938272# AT:0.93827203

def nu(w,q2):
	#print "nu=",(w*w-MP*MP+q2)/(2*MP)
	return (w*w-MP*MP+q2)/(2*MP)	

def epsilon(w,q2):
	n=nu(w,q2)
	e0=E1F_E0
	e1=e0-n

	#! 1. and 2. give same resuts

	#! 1. Original: tan^2(theta/2)=f(Q2,nu) 
	epsInv=1+2*(q2+n*n)/(4*e0*e1-q2)

	#! 2. theta/2=f(Q2,nu) [theta obtain from Q2=4EE'sin^2(theta/2)] 
	#theta=2*math.asin(math.sqrt(q2/(4*e0*(e0-n))))
	#epsInv=1+2*(1+((n*n)/q2))*math.tan(theta/2)**2

	return 1.0/epsInv

def getvgflux(w,q2,e0=E1F_E0):
	eps=epsilon(w,q2)
	#! Original
	return A*w*(w*w-MP*MP)/(4*PI*e0*e0*MP*MP*q2*(1-eps))

	#! test: Original re-written
	# return (A*w*(w**2-MP**2))/(4*PI*e0*e0*MP*MP*q2*(1-eps))

	#! Gothe
	#return ( A/(2*PI**2) )*( (MP**2+2*MP*e0-q2-w**2)/(2*MP*e0) )*( (w**2-MP**2+q2)/(2*MP) )*(1/(q2*(1-eps))) 

	#! Ye: "cleaner"
	# return A*w*(w-MP*MP)/(4*PI*e0*e0*MP*MP*q2*(1-eps))

	#! Ye: direct
	# t1=(A/(4*math.pi))
	# t2=(1/((e0*e0)*(MP*MP)))
	# t3=((w*(w-MP*MP))/((1-eps)*q2))
	# #print "t1,t2,t3=",t1,t2,t3
	# return t1*t2*t3
