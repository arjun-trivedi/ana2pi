# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ### {Pandas DataFrame}{ROOT Hisotgrams}{HDF5 file format}

# <rawcell>

# Test Case for using Pandas DataFrame to store (to HDF5 files) and perform various operations in a "functional" manner on ROOT Histograms.  

# <codecell>

""" Use Case for storing Pandas DataFrame with ROOT Histograms in 
a HDF5 file
Instructions from: http://pandas.pydata.org/pandas-docs/dev/io.html#hdf5-pytables

"""
import os
import numpy as np
randn = np.random.randn
import pandas as pd

from ROOT import TH1D

def make_test_hdf5():
	print 'Going to create store.h5:'
	#1. If store.hdf5 exists, delete it
	if os.path.isfile('store.h5'):
		os.remove('store.h5')
	#2. Create HDF5 store & see what it contains
	store = pd.HDFStore('store.h5')
	print 'store =',store

	#3. Create Histograms to store
	hl = []
	for i in range(0,8):
		name = 'h%d'%(i+1)
		hl.append(TH1D(name,name,100,0,5))
		
	index = np.arange(0,8)
	df = pd.DataFrame()

	#4. Insert 1st columns into df
	if not df:
		data = pd.DataFrame({'s1':hl},index=index) # Data for 1st. Column 
		df = df.append(data)

	#5. Insert 2nd (duplicate) column into df
	df['H'] = hl
	#print 'df=',df

	#6. Add df to store
	store['df']=df
	#7. See the stored data type 
	print 'Type of stored data = ',store.root.df._v_attrs.pandas_type
	#8. See what the file finally contains
	print 'store.h5 contains', store

def read_test_hdf5():
	print 'Going to read contents of store.h5:'
	if os.path.isfile('store.h5'):
		store = pd.HDFStore('store.h5')
		print 'store =',store
		print store['df']
	else:
		print 'store.h5 does not exist!'
        
def get_df_hdf5():
    if os.path.isfile('store.h5'):
		store = pd.HDFStore('store.h5')
		df = store['df']
    else:
        print 'store.h5 does not exist!'
    return df
            
	

# <codecell>

make_test_hdf5()

# <codecell>

read_test_hdf5()

# <codecell>

df = get_df_hdf5()
print len(df.index)
dl=[]
for i in range(0,len(df.index)):
    dl.append(df.iloc[i]['H'])
df['test']=dl
print df
#df['test']=df['H'].Clone("test")

# <codecell>


