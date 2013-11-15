# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ### Following is an example for adding Columns to an empty data frame
# * Take note of the _special_ manner in which the first column is appended

# <codecell>

import pandas as pd
df = pd.DataFrame()
data = pd.DataFrame({"A": range(3)}) # Data for 1st. Column 
df = df.append(data) # Add Data for 1st column
df['B']=data2
df['C']=pd.Series(['the','champions','of'])
print df

# <codecell>


