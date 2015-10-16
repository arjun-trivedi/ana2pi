"""
+This file is used by proc_yields.py to make appropriate
Q2,W projections from h8[W]

+ Q2,W binning should match that in ../h8_bng.h

[02-23-15]
Adjusted binning to not process W>2.125 GeV
"""
# NBINS_WCRS=13
# #! Define Q2,W binning in each Crs-W bin
# Q2W_BNG=[
# [[1.25,5.25,0.5],[1.300,1.425,0.025]],#q2w1
# [[1.25,5.25,0.5],[1.425,1.575,0.025]],#q2w2
# [[1.25,5.25,0.5],[1.575,1.725,0.025]],#q2w3
# [[1.25,5.25,0.5],[1.725,1.850,0.025]],#q2w4
# [[1.25,5.25,0.5],[1.850,2.000,0.025]],#q2w5
# [[1.25,5.25,0.5],[2.000,2.125,0.025]],#q2w6
# [[1.25,5.25,0.5],[2.125,2.275,0.025]],#q2w7
# [[1.25,5.25,0.5],[2.275,2.425,0.025]],#q2w8
# [[1.25,5.25,0.5],[2.425,2.550,0.025]],#q2w9
# [[1.25,5.25,0.5],[2.550,2.700,0.025]],#q2w10
# [[1.25,5.25,0.5],[2.700,2.825,0.025]],#q2w11
# [[1.25,5.25,0.5],[2.825,2.950,0.025]],#q2w12
# [[1.25,5.25,0.5],[2.950,3.000,0.025]],#q2w13
# ]

NBINS_WCRS=6
#! Define Q2,W binning in each Crs-W bin
Q2W_BNG=[
[[1.25,5.25,0.5],[1.300,1.425,0.025]],#q2w1
[[1.25,5.25,0.5],[1.425,1.575,0.025]],#q2w2
[[1.25,5.25,0.5],[1.575,1.725,0.025]],#q2w3
[[1.25,5.25,0.5],[1.725,1.850,0.025]],#q2w4
[[1.25,5.25,0.5],[1.850,2.000,0.025]],#q2w5
[[1.25,5.25,0.5],[2.000,2.125,0.025]] #q2w6
]