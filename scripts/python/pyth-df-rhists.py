# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ### {Pandas DataFrame}{ROOT Hisotgrams}{HDF5 file format}

# <rawcell>

# In this notebook, I want to demonstrate that I can use Pandas DataFrame (DF) to store and manipulate ROOT histograms. 
# 
# I cannot, however, take advantage of vectorized operations for a DF that holds ROOT objects. Nevertheless, once a DF is prepared with ROOT objects, the ease of accessibility of the objects, combined with "familiar" semantics of Python, makes for writing readable and terse code.  

# <markdowncell>

# ####Define 
# * make_test_hdf5(): Create a DF with just one Column(H1) of 10 TH1Ds and store it in store.h5
# * read_test_hdf5(): Read contents of store.h5
# * get_df_hdf5(): Get DF from store.h5

# <codecell>

import os
import numpy as np
randn = np.random.randn
import pandas as pd

from ROOT import TH1D

def make_test_hdf5():
	""" Create a DF with just one Column(H1) of N TH1Ds,
    indexed by [0,...,N] and store it in a HDF5 file.
    Instructions from: http://pandas.pydata.org/pandas-docs/dev/io.html#hdf5-pytables

    """
	print 'Going to create store.h5:'
	#1. If store.hdf5 exists, delete it
	if os.path.isfile('store.h5'):
		os.remove('store.h5')
	#2. Create HDF5 store & see what it contains
	store = pd.HDFStore('store.h5')
	#print 'store =',store

	#3. Create Histograms to store
	NHISTS=10
	hl = []
	for i in range(0,NHISTS):
		name = 'h%d'%(i+1)
		hl.append(TH1D(name,name,100,0,5))
		
	#4. Insert 1st columns into df
	index = np.arange(0,NHISTS)
	df = pd.DataFrame()
	if not df:
		data = pd.DataFrame({'H1':hl},index=index) # Data for 1st. Column 
		df = df.append(data)

	#5. Add df to store
	store['df']=df
	#6. See the stored data type 
	#print 'Type of stored data = ',store.root.df._v_attrs.pandas_type
	#7. See what the file finally contains
	#print 'store=',store

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
df = get_df_hdf5()

# <codecell>

def print_hname(h):
    print h.GetName()
def print_entries(h):
    print h.GetEntries()
def fill1(h):
    h.Fill(1)
def fill2(h):
    h.Fill(2)
    
hs=df['H1']
for i in range(0,len(hs)):
    print_hname(hs[i])
    print_entries(hs[i])
    #fill(hs[i])
    print_entries(hs[i])
    

# <codecell>


