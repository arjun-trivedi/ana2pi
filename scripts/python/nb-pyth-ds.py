# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ## Tutorial to understand basic Data Structures in Python
# (http://pandas.pydata.org/pandas-docs/dev/dsintro.html)

# <codecell>

import numpy as np
randn = np.random.randn

import pandas as pd

# Example of Series (sub-class = NDFrame; earlier was ndarray)
#series with default index [1,2,3...]
s=pd.Series(np.sin(np.arange(0,np.pi,np.pi/4))); 
print s

# <markdowncell>

# ### Example of Series (sub-class = 'NDFrame'; earlier was 'ndarray')
# * Series is __dict-like__ i.e an array of __values__ with __labels__(__indices__)

# <codecell>

# Series with default index
s=pd.Series(np.sin(np.arange(0,np.pi,np.pi/4))); 
print s

# <markdowncell>

# ### DataFrame
# * Can be thought of as a __dict of Series__ object or more generally, as a 2D labeled data structure with columns of potentially different types (as in a spreadsheet)
# * Since __index__ is reserved as labels for __values__ in a Series,  __columns__ are used to access a particular __Series__ inserted in a DF
#     * __index__ = row labels (or label for each element of a series)
#     * __columns__ = column label (or label for each Series in a DF)
#     
#     OR more abstractly (for example if imagining a flat 2D structure or horizontal rows and vertical columns)
#     
#     * __index__ = __horizontal__ label = __axis 0__
#     * __columns__ = __vertical__ label = __axis 1__

# <codecell>

# Example: DF from dict of Series
boy1 = ['arjun','trivedi',768]
boy2 = ['rohit','bagaria',741]
d = {'s1':boy1, 's2':boy2} # dict of Series
df = pd.DataFrame(d)
print df

# <markdowncell>

# ### Example for "basics of indexing are as follows (search this in reference for details):

# <codecell>

print "Select column df[col]: >>>df[\'s1\']"
print df['s1']
print "\n"
print "Select row by label df.loc[label]: >>>df.loc[0]"
print df.loc[0]
print "\n"
print "Select row by integer location df.iloc[loc]: >>>df.iloc[0]"
print df.iloc[0]
print "\n"
print "Slice rows df[rlabel1:rlabel2]: >>>df[0:2]"
print df[0:2]



# <markdowncell>

# ### More sophisticated indexing
# (http://pandas.pydata.org/pandas-docs/dev/indexing.html)

# <markdowncell>

# ### Boolean Indexing
# 
# * This is what I have been using for simstats DF to select data. 
# * The idea is to set up a __boolean vector__ that can be passed to the indexing operator __[]__

# <codecell>

d = pd.DataFrame(randn(3,5), index=['r1','r2','r3'], 
                 columns=['c1','c2','c3','c4','c5'])
print d

# <codecell>

#set up boolean vector =  selection
sel = d['c1']>0.1
print d[sel]

# <codecell>


