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
A=10
B=10
C=10
D=10
f = A + B*cos(x) + C*cos(2*x) + D*sin(x)
fig = plt.figure(figsize=(3,3))
ax = fig.add_subplot(111, title=('fphi%d'%N))
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
plt.legend()


# <markdowncell>

# ### Try to retreive constant multiplying $\sin(\phi)$ using Method 2

# <codecell>

integral=0
for i in range(0,N):
    integral+=r_fphidx[i]*sin(x[i])
    #integral+=r_fphi[i]*sin(x[i])*dx
#print integral
print integral/math.pi

# <codecell>

dl=[0.0003,-0.0009,0.002,0.0068,0.0026,-0.0021,0.0022,0.0049,0.0052,-0.0037]
dlerr = np.empty(10)
dlerr[:]=0.02
xl = np.arange(0,10)
#plt.scatter(xl,dl)
plt.errorbar(xl,dl,yerr=dlerr)

# <codecell>


