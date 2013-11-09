# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ### Tutorial to understand basic Data Structures in Python
# (http://pandas.pydata.org/pandas-docs/dev/dsintro.html)

# <codecell>

import numpy as np
randn = np.random.randn

import pandas as pd

# Example of Series (sub-class = NDFrame; earlier was ndarray)
s=pd.Series(np.sin(np.arange(0,np.pi,np.pi/4))); #series with default index [1,2,3...]
print s

# <codecell>


