#!/usr/bin/python
from __future__ import division
import sys,os,math
from collections import OrderedDict

sys.path.insert(0, "%s/obs_2pi/"%os.environ['ELAST_LITE'])
import disp_obs as disp_obs


'''
Following is optimized for numerical demonstration of polarization observables for
+ alpha =[0deg,36deg] w.r.t phi_pim (VST1-ALPHA, [pf,pip][p,pim], psi_pim_phi_pim) 
+ Q2,W=2.40-2.30,1.725-1.750

+ The A, B, C, D, and E fit pars are obtained from:
$OBSDIR_E163/lowQ2_SSBands_122717_w_phi_fit_curves_for_ananote/cutsncors1/sim4_sim5_sim6_sim7_sim8_sim13/dbg/Obs_R2_EC_EF_ST/mthd_phi-proj-fit_NQ/phiprojs/
VST1_ALPHA/q_2.40-2.30_w1.725-1.750/c_EF(EC)_q2.40-2.30_w1.725-1.750.pdf

+ Using these parameters and other relevant information (see code for details),
 this program calculates the value of the polarization observables.

+ These values can be compared with (and are found to be exactly the same):
A:$OBSDIR_E162/thesis_obs_norm_ST_shape_122817/vm_format/R2/R2-T_plus_R2-L/q2_24_30/w_1725_1750/xsecR2-T_plus_R2-L_psi_pim_phi_pim.dat
B:$OBSDIR_E162/thesis_obs_norm_ST_shape_122817/vm_format/R2/R2-c-LT/q2_24_30/w_1725_1750/xsecR2-c-LT_psi_pim_phi_pim.dat
C:$OBSDIR_E162/thesis_obs_norm_ST_shape_122817/vm_format/R2/R2-c-TT/q2_24_30/w_1725_1750/xsecR2-c-TT_psi_pim_phi_pim.dat
D:$OBSDIR_E162/thesis_obs_norm_ST_shape_122817/vm_format/R2/R2-s-LT/q2_24_30/w_1725_1750/xsecR2-s-LT_psi_pim_phi_pim.dat
E:$OBSDIR_E162/thesis_obs_norm_ST_shape_122817/vm_format/R2/R2-s-TT/q2_24_30/w_1725_1750/xsecR2-s-TT_psi_pim_phi_pim.dat
'''
LUM_E16=28.18
LUM_INVFB_TO_INVMICROB=1000000000
E16_E0=5.754
CF={"A":1.59155,'B':1.61803,'C':1.70131,'D':1.61803, 'E':1.70131}
#! rad-eff-corr for Q2,W=2.00-2.40,1.725-1.750
RADEFFCORR=1.1688429117202759
RADEFFCORR_ERR=0.0012284496823476696
#! etgt-bg-sub-corr
ETGT_BG_SUB_CORR_FCTR=(1+0.0183)
def get_R2(R2n,R2v,q2min,q2max,wmin,wmax,alphamin,alphamax):
    DQ2=q2max-q2min
    DW=wmax-wmin
    q2c=q2min+(DQ2/2)
    wc=wmin+(DW/2)
    DX=math.fabs(math.radians(alphamax)-math.radians(alphamin))
    vgflux=disp_obs.getvgflux(wc,q2c,E16_E0)
    v=(R2v*CF[R2n])/(LUM_E16*LUM_INVFB_TO_INVMICROB*vgflux*DQ2*DW*DX)
    print "q2c,2c=",q2c,wc
    print "vgflux=",vgflux
    print "DX=",DX
    #! apply rad-eff-corr
    v=v*RADEFFCORR
    #! apply etgt-bg-sub-corr
    v=v*ETGT_BG_SUB_CORR_FCTR
    return v

R2=OrderedDict()
for item in zip(["A","B","C","D","E"],[2622.22, -657.15, 54.45, 265.50, -220.61]):
    r2n,r2v=item[0],item[1]
    R2[r2n]=get_R2(r2n,r2v,2.4,3.0,1.725,1.750,0,36)
print R2


# <codecell>


