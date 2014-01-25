# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from __future__ import division
MP=0.93827203
E1F_E0=5.497
FSC=0.00729735253
A=FSC
LUM=20#in fb^-1
FACTOR=1000000000#to obtain final xsec in nb

def nu(w,q2):
    return (w*w-MP*MP+q2)/(2*MP);

def epsilon(w,q2):
    n = nu(w,q2);
    e0 = E1F_E0;
    e1 = e0-n;
    epsInv = 1+2*(q2+n*n)/(4*e0*e1-q2);
    return 1.0/epsInv;

def getvgflux(w,q2,e0 = E1F_E0):
    eps = epsilon(w,q2);
    return A*w*(w*w-MP*MP)/(4*pi*e0*e0*MP*MP*q2*(1-eps));

def get_dsigma_dwdq2(w,q2,dw,dq2,n,dn_sys):
    norm=LUM*getvgflux(w,q2)*dq2*dw*FACTOR
    xsec=n/norm
    stt_abs_err = sqrt(n)/norm
    stt_rel_err=stt_abs_err/xsec*100
    sys_abs_err=((dn_sys/100)*n)/norm
    sys_rel_err=sys_abs_err/xsec*100
    res=[xsec,stt_abs_err,stt_rel_err,sys_abs_err,sys_rel_err]
    return res

# <codecell>

res = get_dsigma_dwdq2(1.625,2.0,0.025,0.4,201000,2)
print "dsigma_dwdq2=%5.3f(+/-%5.3f)(+/-%5.3f)"%(res[0],res[1],res[3])
print "Stat-Err(abs)=%5.3f"%res[1]
print "Stat-Err(rel)=%5.3f%%"%res[2]
print "Syst-Err(abs)=%5.3f"%res[3]
print "Syst-Err(rel)=%5.3f%%"%res[4]

# <codecell>


