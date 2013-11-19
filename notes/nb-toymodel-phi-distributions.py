# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ## Toy Model of $\phi$ distributions and extraction of parameters
# 
# Assume following form of $\phi$ distribution (see nb-beamasym-obs-mtg-091013)
# 
# $$ \left(\frac{d^2\sigma_{v}}{ {dX^{ij}_{k}d\phi^{i}} }\right)^{h}  = A_{k} + B_{k}\cos\phi + C\cos2\phi + hP(\mathcal{R}^{ij}_{k})_{LT^{'}}\sin\phi $$
# 
# >  where
# 
# * $A,B,C,\mathcal{R}^{ij}_{k} = f(Q^{2},W,X^{ij}_{k})$

# <markdowncell>

# * Assume a form for fphi = 1 + 10*cos(x) + 20*cos(2*x) + 100*sin(x)
# * Given binning resolution in x, see how well fphi can be "reconstructed" (r_fphi)

# <codecell>

from scipy import integrate
import math

#Ideal phi distribution(fphi)
N = 10
dx = 2*math.pi/N
x  = arange(0,2*math.pi,dx)
A=16
B=3.1
C=3.6
D=-0.021
f = A + B*cos(x) + C*cos(2*x) + D*sin(x)
fsinphi = sin(x)
fig = plt.figure(figsize=(3,3))
ax = fig.add_subplot(111, title=('fphi%d'%N))
plt.scatter(x,fsinphi,marker='o', color='r',s=50,label='fsinphi')
plt.scatter(x,f,marker='^', color='g',s=50,label='fphi')

#Reconstructed phi distribution(r_fphi)
fphi  = lambda x: A + B*cos(x) + C*cos(2*x) + D*sin(x)
def int_fphi(x,dx):
    return integrate.quad(fphi,x,x+dx)[0]

r_fphidx=[]
r_fphi=[]
for i in x:
    #add a statistical modelled error to integral!
    integral = int_fphi(i,dx)
    r_fphidx.append(integral)
    r_fphi.append(divide(integral,dx))

#ax.scatter(x,r_fphidx,marker='v',color='r',s=50, label='r_fphidx')
ax.scatter(x,r_fphi,marker='<',color='b',s=50,label='r_fphi')
#plt.legend()

integral=0
for i in range(0,N):
    integral+=r_fphidx[i]*sin(x[i])
    #integral+=r_fphi[i]*sin(x[i])*dx
#print integral
print 'integral =', integral/math.pi

# <codecell>

#Ideal phi distribution(fphi)
N = 10
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
print 'integral_B =',integral['B']/math.pi
print 'integral_C =',integral['C']/math.pi
print 'integral_D =',integral['D']/math.pi
print 'integral_T =',integral['T']/math.pi

#integral_T_check = integral['B']+integral['C']+integral['D']
#print 'cross check by directly adding integrals = ',integral_T_check

# <codecell>

t = {}
t['A']= t['B']=1
print t

# <codecell>


