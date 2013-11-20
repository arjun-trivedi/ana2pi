# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ## Toy Model of $\phi$ distributions and extraction of parameters
# 
# Assume following form for $\phi$ Cross Sections (see nb-beamasym-obs-mtg-091013)
# 
# $$ \left(\frac{d^2\sigma_{v}}{ {dX^{ij}d\phi^{j}} }\right)^{h}  = A^{ij} + B^{ij}\cos\phi^{j} + C^{ij}\cos2\phi^{j} + D^{ij}\sin\phi^{j} $$
# 
# >  where
# 
# * $A^{ij},B^{ij},C^{ij},D^{ij} = f(Q^{2},W,X^{ij})$

# <markdowncell>

# Assume a form for the $\phi$ projection of the Double Differential Cross Section for a particular $X^{ij}$ bin
# 
# * `diffxsec(x) = A + Bcos(x) + Ccos(2x) + Dsin(x)`
# 
# Given binning resolution in $\phi$, see 
# 
# 1. How well `xsec(x)` can be `Reconstructed` 
# 2. How well `Method 3.` can then be used to extract `R2`

# <codecell>

from __future__ import division
from scipy import integrate
import math

N = 10
xmax=360
xmin=0
dx=(xmax-xmin)/N
x=arange(xmin,xmax,dx)
#print 'x =',x
#print 'dx =',dx
#dx = 2*math.pi/N
#x  = arange(0,2*math.pi,dx)
A=15#16
B=2.3#3.1
C=2.7#3.6
D=0.0026#-0.021
diffxsec = A + B*cos(radians(x)) + C*cos(radians(2*x)) + D*sin(radians(x))
f_sinx = sin(radians(x))

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111, title='R2 Extraction Test (Method 3)')

#plt.scatter(x,fsinphi,marker='o', color='r',s=50,label='fsinphi')
#plt.scatter(x,f,marker='^', color='g',s=50,label='fphi')
#marker='^', color='g',s=50,label='fphi')
#ax.plot(

#Reconstructed phi distribution(r_fphi)
diffxsec_lambda  = lambda t: A + B*cos(radians(t)) + C*cos(radians(2*t)) + D*sin(radians(t))
def integrate_diffxsec(xmin,xmax):
    return integrate.quad(diffxsec_lambda,xmin,xmax)[0]

r_xsec=[]
#r_n=[]
for i in range(0,N):
    #add a statistical modelled error to integral!
    xsec = integrate_diffxsec(x[i],x[i]+dx)
    r_xsec.append(xsec)
    #r_nphi.append(integral/dx)

#ax.scatter(x,r_fphidx,marker='o',color='r',s=50, label='r_fphidx')
ax.scatter(x,r_xsec,marker='.',color='r',s=50,label='r_fphi')
ax.set_xticks(x)
#ax.set_yticks(arange(int(min(r_xsec)),int(max(r_xsec))+2,1))
ax.grid()
#plt.legend()

#integral=0
#for i in range(0,N):
    #integral+=r_nphidx[i]*sin(radians(x[i]))
    #integral+=r_fphi[i]*sin(x[i])*dx
#print integral
#print 'coefficient, D =', integral/math.pi

# <codecell>

#Idealphi distribution(fphi)
N = 3
dx = 2*math.pi/N
x  = arange(0,2*math.pi,dx)
A=16
B=3.1
C=3.6
D=-0.021

f_A = A
f_B = B*cos(x)
f_C = C*cos(2*x)
f_D = D*sin(x)
f = f_A + f_B + f_C + f_D #A + B*cos(x) + C*cos(2*x) + D*sin(x)

fsinphi = sin(x)

fig,axs = plt.subplots(4, sharex=True,figsize=(5,8))

axs[0].scatter(x,f_B,marker='o', color='k',s=50,label='f_B')
axs[0].set_ylabel('f_B')
#axs[0].scatter(x,fsinphi,marker='o', color='r',s=50,label='f_sinphi')

axs[1].scatter(x,f_C,marker='o', color='k',s=50,label='f_C')
axs[1].set_ylabel('f_C')
#axs[1].scatter(x,fsinphi,marker='o', color='r',s=50,label='f_sinphi')

axs[2].scatter(x,f_D,marker='o', color='k',s=50,label='f_D')
axs[2].set_ylabel('f_D')
#axs[2].scatter(x,fsinphi,marker='o', color='r',s=50,label='f_sinphi')

axs[3].scatter(x,f,marker='o', color='k',s=50,label='f')
axs[3].set_ylabel('f')
#axs[3].scatter(x,fsinphi,marker='o', color='r',s=50,label='f_sinphi')

integral = {}
integral['B']=integral['C']=integral['D']=integral['T']=0
for i in range(0,N):
   integral['B'] += f_B[i]*sin(x[i])*dx
   integral['C'] += f_C[i]*sin(x[i])*dx
   integral['D'] += f_D[i]*sin(x[i])*dx
   integral['T'] += f[i]*sin(x[i])*dx
   #integral+=r_fphi[i]*sin(x[i])*dx
#print integral
print 'B =',integral['B']/math.pi
print 'C =',integral['C']/math.pi
print 'D =',integral['D']/math.pi
print 'T =',integral['T']/math.pi

#integral_T_check = integral['B']+integral['C']+integral['D']
#print 'cross check by directly adding integrals = ',integral_T_check

# <codecell>

ft  = lambda t: sin(radians(t))#*sin(radians(t))
print 'using degrees = ',integrate.quad(ft,0,180)[0]

ft  = lambda t: sin(t)#*sin(t)
print 'using radians = ',integrate.quad(ft,0,math.pi)[0]

#def int_fphi(x,dx):
 #   return integrate.quad(fphi,x,x+dx)[0]

# <codecell>

sin(radians(90))

# <codecell>

print divide(360,100)

# <codecell>


