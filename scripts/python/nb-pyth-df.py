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

df_grpd_A = df.groupby('A')
for grp in df_grpd_A.groups:
    print grp
    d=df_grpd_A.get_group(grp)
    print d
    #name = d.values[0][0]
    #print name
    
    if d['B']==768:
        print 'here'
        num1 = d['A']#,d['B'],d['C']
        num2 = d['B']
        print num1
        print num2
        print num1+num2


# <codecell>


