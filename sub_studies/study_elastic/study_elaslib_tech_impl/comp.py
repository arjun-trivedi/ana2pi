#!/usr/bin/python

from __future__ import division
import os
import ROOT
import matplotlib.pyplot as plt
from rootpy.interactive import wait

"""
Compare technical implementations, in Python and C++, of different versions, by date, of elaslib.f. 
+ Beam energy=5.499GeV,Wcut=1.028GeV,E1F_TARGET_RAD_LENGTH=0.00562
+ Result: No difference between technical implementation x version
+ For comparing measured Elastic xsec, I am using python-2012. 
+ The dates are as per the date Gleb first gave me the respective elaslib.f
+ Note the 2012 version elaslib.f, for both python and C++, differ similarly from their original + The 2015 version is the original version


	2012 | 2015
python | X   |  - |
C++    | X   |  X |

Result: python_2012=C++_2012 and python_2012=C++_2015 => C++_2012=C++_2015
"""
NCOMPS=3
P12,C12,C15=range(NCOMPS)

F=[]
F.append(ROOT.TFile("python/2012/TT.root"))
F.append(ROOT.TFile("CPP/2012/TT.root"))
F.append(ROOT.TFile("CPP/2015/TT.root"))
print "Files used:"
print [f.GetName() for f in F]

H=[]
H.append(F[P12].Get("hTTnorm"))
H.append(F[C12].Get("hTTnorm"))
H.append(F[C15].Get("hTTnorm"))
H[P12].SetTitle("hTTnorm-P12")
H[C12].SetTitle("hTTnorm-CPP12")
H[C15].SetTitle("hTTnorm-CPP15")
print "Histograms used:"
print [h.GetTitle() for h in H]

#! make dsigma/domega plots
c_xsec=ROOT.TCanvas()
c_xsec.Divide(2,2)
c_xsec.cd(1)
c_xsec.cd(1).SetLogy()
H[P12].Draw()
c_xsec.cd(3)
c_xsec.cd(3).SetLogy()
H[C12].Draw()
c_xsec.cd(4)
c_xsec.cd(4).SetLogy()
H[C15].Draw()


#! Make comp plots
h_rto_p12_c12=H[P12].Clone()
h_rto_p12_c12.Divide(H[C12])
h_rto_p12_c12.SetTitle("Python-12/C++-12")

h_rto_p12_c15=H[P12].Clone()
h_rto_p12_c15.Divide(H[C15])
h_rto_p12_c15.SetTitle("Python-12/C++-15")

c_comp=ROOT.TCanvas()
c_comp.Divide(2,1)
c_comp.cd(1)
h_rto_p12_c12.Draw()
c_comp.cd(2)
h_rto_p12_c15.Draw()

#! If wanting to keep TCanvas open till program exits                           
if not ROOT.gROOT.IsBatch():
       plt.show()
       # wait for you to close the ROOT canvas before exiting
       wait(True)

