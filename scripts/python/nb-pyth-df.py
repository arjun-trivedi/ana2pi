# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ### Following is an example for adding Columns to an empty data frame
# * Take note of the _special_ manner in which the first column is appended

# <codecell>

import pandas as pd
df = pd.DataFrame()
if not df:
    print 'df is empty'

rindex=['name','surname','roll#']   
boy1 = ['saswat','sarda',772]
data = pd.DataFrame({"A": boy1},index=rindex) # Data for 1st. Column 
df = df.append(data) # Add Data for 1st column
if not df:
    print 'df is empty after 1st append'

boy2 = ['arjun','trivedi',768]
df['B']=boy2

boy3 = ['rohit','bagaria',741]
df['C']=boy3

print df

# <codecell>


