# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from __future__ import division

q2min=np.arange(2,5,0.5)
print 'q2min = ',q2min
q2max=np.arange(2.5,5.5,0.5)
print 'q2max = ',q2max

sigma = (-1/3)*((1/q2max**3)-(1/q2min**3))
print 'sigma = ',sigma
norm_sigma = sigma/sigma[0]
print 'norm_sigma = ',norm_sigma
print 'sum(norm_sigma)=',sum(norm_sigma)

# <codecell>

1500*1.9

# <codecell>

qsq = q2min
print 'q2 = ',qsq
sigma_naive = 1/qsq**4
print 'sigma_naive = ',sigma_naive
norm_sigma_naive = sigma_naive/sigma_naive[0]
print 'norm_sigma_naive = ',norm_sigma_naive
print 'sum(norm_sigma_naive)=',sum(norm_sigma_naive)

# <codecell>

q2l=np.linspace(2,5,5)
print 'q2l = ',q2min

sigma = (-1/3)*((1/(q2l+0.75)**3)-(1/q2l**3))
print 'sigma = ',sigma
norm_sigma = sigma/sigma[0]
print 'norm_sigma = ',norm_sigma
print 'sum(norm_sigma)=',sum(norm_sigma)

# <codecell>

20/.03

# <codecell>


