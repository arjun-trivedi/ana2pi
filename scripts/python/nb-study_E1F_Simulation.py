# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import ROOT
from root_numpy import root2array, root2rec, tree2rec
from root_numpy.testdata import get_filepath

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

f = '/data/trivedia/e1f/simulation_2pi/setup_sim_CentOS6/gpppars/try1/q2w2_gpptest_14_011914/recon/d2pi.root'
arr = root2array(f,'d2piR/tR',start=1,stop=5)#,start=1,stop=5)#start=1,stop=2)
df = pd.DataFrame(arr)
print df

# <codecell>


