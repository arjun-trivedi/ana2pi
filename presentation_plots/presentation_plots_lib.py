#!/usr/bin/python
import numpy as np

#! Q2,W for thesis
Q2MIN,Q2MAX=2.0,5.0
Q2BINL=[2.0,2.4,3.0,3.5,4.2,5.0]
WMIN,WMAX,WBINW=1.400,2.125,0.025
WBINL=np.arange(WMIN,WMAX+WBINW,WBINW)

#! Q2 ranges for low and high sim
NSIMRNG=2
L,H=range(NSIMRNG)
SIMRNG_IDTFR=["lowQ2","highQ2"]
Q2MINSIM =[0 for i in range(NSIMRNG)]
Q2MAXSIM =[0 for i in range(NSIMRNG)]
Q2BINLSIM=[0 for i in range(NSIMRNG)]
#! Q2L
Q2MINSIM[L],Q2MAXSIM[L]=2.0,3.0
Q2BINLSIM[L]=[2.0,2.4,3.0]
#! Q2H
Q2MINSIM[H],Q2MAXSIM[H]=3.0,5.0
Q2BINLSIM[H]=[3.0,3.5,4.2,5.0]
