#!/usr/bin/python
from __future__ import division
import ROOT

from sympy.solvers import solve
from sympy import Symbol,sympify

'''
+ All cuts setup here are taken from EI's e16 analysis.
+ The original(in Fortrain) implementation of these cuts for h10->d2pi is in $ELAST_LITE/cut_theta_vs_p_e16.f. 
+ In this file the Fortran implementation of the cuts is converted to display their functional form using ROOT::TF1.
+ These functional forms will be used to justify the cuts and subsquently be used in my thesis.

'''

#! Setup sector information 
NSCTR=6
#! ***

#! Setup particle information
NPRTCL=3
E,P,PIP=range(NPRTCL)
PRTCL_NAME=["e","p","pip"]
#! ***

#! Setup cut
NCUTLNS=2
H,L=range(NCUTLNS)
#! CUT[h/l][prtcl][sctr]
CUT=[[[[]for k in range(NSCTR)]for j in range(NPRTCL)]for i in range(NCUTLNS)]

#! Set up electron cuts
#! + Cut exist for sectors: 2,5,6
#! + Sector 1,3,4 (no cuts)

#! First set all cuts=None
for sctr in [1,2,3,4,5,6]:
	isctr=sctr-1
	CUT[H][E][isctr]=CUT[L][E][isctr]=None

#! Sector 2
#! Remove: The following 3 cuts replaced as they were deemed unnecessary
#CUT[H][E][1]=[]
#CUT[L][E][1]=[]
#! cut1
#CUT[H][E][1].append(ROOT.TF1("H","35.5-4.0*sqrt(x-1.5)",2.00,3.05))
#CUT[L][E][1].append(ROOT.TF1("L","34.0-4.0*sqrt(x-1.5)",2.00,3.05))
#! cut2
#CUT[H][E][1].append(ROOT.TF1("H","27.3-4.8*sqrt(x-1.5)",2.60,4.00))
#CUT[L][E][1].append(ROOT.TF1("L","26.7-4.8*sqrt(x-1.5)",2.60,4.00))
#! cut3
#CUT[H][E][1].append(ROOT.TF1("H","28.5-4.8*sqrt(x-1.5)",2.60,4.00))
#CUT[L][E][1].append(ROOT.TF1("L","26.1-4.8*sqrt(x-1.5)",2.60,4.00))

#! Sector 5
CUT[H][E][4]=[]
CUT[L][E][4]=[]
#! cut1
CUT[H][E][4].append(ROOT.TF1("H","20.5",               2.90,4.20))
CUT[L][E][4].append(ROOT.TF1("L","21.3-2.4*sqrt(x-2.)",2.90,4.20))
#! cut2
CUT[H][E][4].append(ROOT.TF1("H","26.3",2.40,3.40))
CUT[L][E][4].append(ROOT.TF1("L","24.7",2.40,3.40))

#! Sector 6
#! Remove: The following 3 cuts replaced as they were deemed unnecessary
#CUT[H][E][5]=[]
#CUT[L][E][5]=[]
#! cut1
#CUT[H][E][5].append(ROOT.TF1("H","31.0-5.0*sqrt(x-1.5)",2.40,3.70))
#CUT[L][E][5].append(ROOT.TF1("L","29.7-5.0*sqrt(x-1.5)",2.40,3.70))

#! Set up proton cuts
#! + Cut exist for sectors:2,3,5,6 
#! Sector 1,4 (no cuts)

#! First set all cuts=None
for sctr in [1,2,3,4,5,6]:
        isctr=sctr-1
        CUT[H][P][isctr]=CUT[L][P][isctr]=None

#! Sector 2
CUT[H][P][1]=[]
CUT[L][P][1]=[]
#! cut1
CUT[H][P][1].append(ROOT.TF1("H","27.2",0.60,4.00))
CUT[L][P][1].append(ROOT.TF1("L","24.7",0.60,4.00))

#! Sector 3
#! Remove: The following EI-cut replaced by scpd cut 324
#CUT[H][P][2]=[]
#CUT[L][P][2]=[]
#! cut1
#CUT[H][P][2].append(ROOT.TF1("H","37.7+6.0*sqrt(x-0.9)",0.90,2.00))
#CUT[L][P][2].append(ROOT.TF1("L","31.5+9.0*sqrt(x-0.9)",0.90,2.00))

#! Sector 5
CUT[H][P][4]=[]
CUT[L][P][4]=[]
#! cut1
CUT[H][P][4].append(ROOT.TF1("H","25.2",1.00,3.80))
CUT[L][P][4].append(ROOT.TF1("L","24.0",1.00,3.80))
#! cut2
CUT[H][P][4].append(ROOT.TF1("H","11.2",1.00,4.40))
CUT[L][P][4].append(ROOT.TF1("L","10.1",1.00,4.40))
#! cut3
CUT[H][P][4].append(ROOT.TF1("H","11.0+3.50*sqrt(x-0.9)",0.90,4.00))
CUT[L][P][4].append(ROOT.TF1("L","13.5-0.38*(x-4)*(x-4)",0.90,4.00))

#! Sector 6
CUT[H][P][5]=[]
CUT[L][P][5]=[]
#! cut1
CUT[H][P][5].append(ROOT.TF1("H","29.2-1.9*(x-3.2)*(x-3.2)",0.90,2.60))
CUT[L][P][5].append(ROOT.TF1("L","26.0-1.3*(x-3.5)*(x-3.5)",0.90,2.60))

#! Set up pip cuts
#! + Note that for pions EI used p(theta) cuts. Therefore, in order to maintain consistency
#!   with the cuts for the electron and proton, I  have had to obtain theta(p) after inverting 
#!   EI's p(theta). Following are the two main functional forms used by EI and the inverted form
#!   that I obtained:
#! 	+ p=a+b(theta+c)^2 -> theta=sqrt((p-a)/b)-c
#!	+ p=a+b*theta      -> theta=(p-a)/b
#! + Additionally, as a result of the above, I had to obtain the corresponding p_limits for 
#!   the theta(p) cuts. This was done by requiring theta(p) to equal EI's determined theta_limits
#!   and then solving for p.
#! + Still additionally note that as a result of this inversion, EI's high and low cuts defined
#!   with respect to p(theta) get swapped when I invert them to obtain theta(p):
#!	+ p(theta)_H -> theta(p)_L
#!	+ p(theta)_L -> theta(p)_H
#! 
#! + Cut exist for sectors:2,3,5,6 
#! Sector 1,4 (no cuts)

#! First set all cuts=None
for sctr in [1,2,3,4,5,6]:
        isctr=sctr-1
        CUT[H][PIP][isctr]=CUT[L][PIP][isctr]=None

#! + Define x as a symbol that can therefore be used in algebraic equations
#!   and solved for using sympy.solvers.solve
#! + This is needed to obtain p_limits as listed previously
x=Symbol('x')

#! Define a function that encapsulates the process of setting up
#! cuts for pip
def setup_cut_pip(sctr,fh,fl,fh_str,fl_str,theta_min,theta_max):
	isctr=sctr-1
	#! Determine pmin,pmax for each of the cuts(cuth and cutl)
	cuth_pmin=solve(fh-theta_min,x)
	cuth_pmax=solve(fh-theta_max,x)
	cutl_pmin=solve(fl-theta_min,x)
	cutl_pmax=solve(fl-theta_max,x)
	#! Now define ROOT.TF1
	CUT[H][PIP][isctr].append(ROOT.TF1("H",fh_str,cuth_pmin[0],cuth_pmax[0]))
	CUT[L][PIP][isctr].append(ROOT.TF1("L",fl_str,cutl_pmin[0],cutl_pmax[0]))

#! Sector 2
CUT[H][PIP][1]=[]
CUT[L][PIP][1]=[]
#! cut1=cut1.1+cut1.2
#! cut1.1 (no dependence on p)
CUT[H][PIP][1].append(ROOT.TF1("H","26",0.72,0.64))
CUT[L][PIP][1].append(ROOT.TF1("L","6", 0.72,0.64))
#! cut1.2
fh=((x-0.68)/0.002070)**(1/2)-(-24.96)
fl=((x-0.75)/0.003116)**(1/2)-(-24.64)
fh_str="((x-0.68)/0.002070)**(1/2)-(-24.96)"
fl_str="((x-0.75)/0.003116)**(1/2)-(-24.64)"
theta_min,theta_max=26,46
setup_cut_pip(2,fh,fl,fh_str,fl_str,theta_min,theta_max)
#! cut2
fh=((x-(-1.93))/0.0235)
fl=((x-(-2.65))/0.0350)
fh_str="((x-(-1.93))/0.0235)"
fl_str="((x-(-2.65))/0.0350)"
theta_min,theta_max=80,98
setup_cut_pip(2,fh,fl,fh_str,fl_str,theta_min,theta_max)
#! cut3
#! Remove: The following EI-cut replaced by scpd cut 245
fh=(x-(-2.21))/0.0235
fl=(x-(-2.07))/0.0235
fh_str="(x-(-2.21))/0.0235"
fl_str="(x-(-2.07))/0.0235"
theta_min,theta_max=94,107
#setup_cut_pip(2,fh,fl,fh_str,fl_str,theta_min,theta_max)

#! Sector 3
#! Removed all cuts for sector 3: See below
#CUT[H][PIP][2]=[]
#CUT[L][PIP][2]=[]
#! cut1=cut1.1+cut1.2
#! Remove: The following EI-cut replaced by scpd cut 324
#! cut1.1 (no dependence on p)
#CUT[H][PIP][2].append(ROOT.TF1("H","26.3",0.60,0.69))
#CUT[L][PIP][2].append(ROOT.TF1("L","6"   ,0.60,0.69))
#! cut1.2
fh=((x-0.611)/0.002070)**(1/2)-(-24.96)
fl=((x-0.700)/0.003116)**(1/2)-(-24.64)
fh_str="((x-0.611)/0.002070)**(1/2)-(-24.96)"
fl_str="((x-0.700)/0.003116)**(1/2)-(-24.64)"
theta_min,theta_max=26,46
#setup_cut_pip(3,fh,fl,fh_str,fl_str,theta_min,theta_max)
#! cut2
#! Remove: The following EI-cut replaced by scpd cut 337,338
fh=((x-0.193)/0.000840)**(1/2)-(-54.04)
fl=((x-0.308)/0.001466)**(1/2)-(-54.77)
fh_str="((x-0.193)/0.000840)**(1/2)-(-54.04)"
fl_str="((x-0.308)/0.001466)**(1/2)-(-54.77)"
theta_min,theta_max=52,79
#setup_cut_pip(3,fh,fl,fh_str,fl_str,theta_min,theta_max)
#! cut3
#! Remove: The following EI-cut replaced by scpd cut 342
fh=(x-(-0.68))/0.011
fl=((x-0.190)/0.00091)**(1/2)-(-68.80)
fh_str="(x-(-0.68))/0.011"
fl_str="((x-0.190)/0.00091)**(1/2)-(-68.80)"
theta_min,theta_max=73,98
#setup_cut_pip(3,fh,fl,fh_str,fl_str,theta_min,theta_max)
#! cut4
#! Remove: The following EI-cut replaced by scpd cut 345,346,347
fh=130#!purposely created to be out of bounds so that the upper bound is not displayed. This is because EI has no upper bound defined i.e. cut is only on lower bound
fl=(x-(-2.02))/0.0235
fh_str="130"
fl_str="(x-(-2.02))/0.0235"
theta_min,theta_max=94,130
cuth_pmin=[0]#solve(fh-theta_min,x) NA for this cut
cuth_pmax=[5]#solve(fh-theta_max,x) NA for this cut
cutl_pmin=solve(fl-theta_min,x)
cutl_pmax=solve(fl-theta_max,x)
#CUT[H][PIP][2].append(ROOT.TF1("H",fh_str,cuth_pmin[0],cuth_pmax[0]))
#CUT[L][PIP][2].append(ROOT.TF1("L",fl_str,cutl_pmin[0],cutl_pmax[0]))

#! Sector 4
#! Remove: Cuts for this sector removed. See below
#CUT[H][PIP][3]=[]
#CUT[L][PIP][3]=[]
#! cut1
#! Remove: The following EI-cut replaced by scpd cut 434
fh=((x-0.266)/0.000766)**(1/2)-(-43.10)
fl=((x-0.317)/0.001050)**(1/2)-(-42.60)
fh_str="((x-0.266)/0.000766)**(1/2)-(-43.10)"
fl_str="((x-0.317)/0.001050)**(1/2)-(-42.60)"
theta_min,theta_max=42,69
#setup_cut_pip(4,fh,fl,fh_str,fl_str,theta_min,theta_max)

#! Sector 5
#! Remove: Cuts for this sector removed. See below
#CUT[H][PIP][4]=[]
#CUT[L][PIP][4]=[]
#! cut1
#! Remove: The following EI-cut replaced by scpd cut 542
fh=((x-0.101)/0.00020)**(1/2)-(-56.7)
fl=((x-0.135)/0.00064)**(1/2)-(-64.2)
fh_str="((x-0.101)/0.00020)**(1/2)-(-56.7)"
fl_str="((x-0.135)/0.00064)**(1/2)-(-64.2)"
theta_min,theta_max=76,95
#setup_cut_pip(5,fh,fl,fh_str,fl_str,theta_min,theta_max)

#! Sector 6
CUT[H][PIP][5]=[]
CUT[L][PIP][5]=[]
#! cut1
fh=(x-(-0.12))/0.06
fl=-10#!purposely created to be out of bounds so that the lower bound is not displayed. This is because EI has no lower bound defined i.e. cut is only on upper bound
fh_str="(x-(-0.12))/0.06"
fl_str="-10"
theta_min,theta_max=8,36
cuth_pmin=solve(fh-theta_min,x)
cuth_pmax=solve(fh-theta_max,x)#solve(fh-theta_max,x) NA for this cut
cutl_pmin=[0]#solve(fh-theta_min,x) NA for this cut
cutl_pmax=[5]#solve(fh-theta_min,x) NA for this cut
CUT[H][PIP][5].append(ROOT.TF1("H",fh_str,cuth_pmin[0],cuth_pmax[0]))
CUT[L][PIP][5].append(ROOT.TF1("L",fl_str,cutl_pmin[0],cutl_pmax[0]))
#! cut2
fh=((x-0.47)/0.0036)**(1/2)-(-22.7)
fl=((x-0.62)/0.0040)**(1/2)-(-22.7)
fh_str="((x-0.47)/0.0036)**(1/2)-(-22.7)"
fl_str="((x-0.62)/0.0040)**(1/2)-(-22.7)"
theta_min,theta_max=23,37
setup_cut_pip(6,fh,fl,fh_str,fl_str,theta_min,theta_max)
#! cut3
fh=((x-0.34)/0.0040)**(1/2)-(-32.8)
fl=((x-0.46)/0.0044)**(1/2)-(-32.8)
fh_str="((x-0.34)/0.0040)**(1/2)-(-32.8)"
fl_str="((x-0.46)/0.0044)**(1/2)-(-32.8)"
theta_min,theta_max=33,45
setup_cut_pip(6,fh,fl,fh_str,fl_str,theta_min,theta_max)
#! cut4
#! Remove: The following EI-cut replaced by scpd cut 644
fh=(x-(-1.250))/0.015
fl=(x-(-1.075))/0.015
fh_str="(x-(-1.250))/0.015"
fl_str="(x-(-1.075))/0.015"
theta_min,theta_max=82,106
#setup_cut_pip(6,fh,fl,fh_str,fl_str,theta_min,theta_max)
#! cut5
#! Remove: The following EI-cut replaced by scpd cut 632,633,634,635,636,637,638,639
fh=((x-0.25)/0.00054)**(1/2)-(-48.5)
fl=(x-(-0.64))/0.027
fh_str="((x-0.25)/0.00054)**(1/2)-(-48.5)"
fl_str="(x-(-0.64))/0.027"
theta_min,theta_max=35.5,77
#setup_cut_pip(6,fh,fl,fh_str,fl_str,theta_min,theta_max)
