# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ## Toy Model of $\phi$ distributions and extraction of parameters
# 
# Assume following form of $\phi$ distribution (see nb-beamasym-obs-mtg-091013)
# 
# $$ \left(\frac{d^2\sigma_{v}}{ {dX^{ij}d\phi^{j}} }\right)^{h}  = A^{ij} + B^{ij}\cos\phi^{j} + C^{ij}\cos2\phi^{j} + D^{ij}\sin\phi^{j} $$
# 
# >  where
# 
# * $A,B,C,\mathcal{R}^{ij}_{k} = f(Q^{2},W,X^{ij}_{k})$

# <markdowncell>

# * Assume a form for fphi = 1 + 10*cos(x) + 20*cos(2*x) + 100*sin(x)
# * Given binning resolution in x, see how well fphi can be "reconstructed" (r_fphi)

# <codecell>

from __future__ import division
from scipy import integrate
import math

#Ideal phi distribution(fphi)
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
f = A + B*cos(radians(x)) + C*cos(radians(2*x)) + D*sin(radians(x))
fsinphi = sin(radians(x))
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111, title=('fphi%d'%N))

#plt.scatter(x,fsinphi,marker='o', color='r',s=50,label='fsinphi')
plt.scatter(x,f,marker='^', color='g',s=50,label='fphi')

#Reconstructed phi distribution(r_fphi)
fphi  = lambda t: A + B*cos(radians(t)) + C*cos(radians(2*t)) + D*sin(radians(t))
def int_fphi(x,dx):
    return integrate.quad(fphi,x,x+dx)[0]

r_fphidx=[]
r_fphi=[]
for i in range(0,N):
    #add a statistical modelled error to integral!
    integral = int_fphi(x[i],dx)
    r_fphidx.append(integral)
    r_fphi.append(integral/dx)

#ax.scatter(x,r_fphidx,marker='o',color='r',s=50, label='r_fphidx')
ax.scatter(x,r_fphi,marker='.',color='r',s=50,label='r_fphi')
ax.set_xticks(x)
ax.set_yticks(arange(int(min(r_fphi)),int(max(r_fphi))+2,1))
ax.grid()
#plt.legend()

integral=0
for i in range(0,N):
    integral+=r_fphidx[i]*sin(radians(x[i]))
    #integral+=r_fphi[i]*sin(x[i])*dx
#print integral
print 'coefficient, D =', integral/math.pi

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


