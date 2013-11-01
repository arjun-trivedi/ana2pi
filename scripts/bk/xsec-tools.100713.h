#include <TString.h>
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>

#include "constants.h"

/* private-type variables */
const TString Q2W_BNG1 = "1-1.400-1.500__8-1.600-1.800";
const TString Q2W_BNG2 = "1-2.000-2.400__24-1.300-1.900";
const TString AT_TOPNAMES[5] = 
                        {
                          "p#pi^{+}#pi^{-}", 
                          "p#pi^{+}#slash{#pi}^{-}", 
                          "p#slash{#pi}^{+}#pi^{-}", 
                          "#slash{p}#pi^{+}#pi^{-}", 
                          "p#pi^{+}#pi^{-} + p#pi^{+}#slash{#pi}^{-} + p#slash{#pi}^{+}#pi^{-} + #slash{p}#pi^{+}#pi^{-}"
                        };
const TString VM_TOPNAMES[5] = 
                        {
                          "p#pi^{+}#pi^{-}", 
                          "p#pi^{+}#slash{#pi}^{-} + p#pi^{+}#pi^{-}", 
                          "p#slash{#pi}^{+}#pi^{-} + p#pi^{+}#pi^{-}", 
                          "#slash{p}#pi^{+}#pi^{-} + p#pi^{+}#pi^{-}", 
                          "p#pi^{+}#pi^{-} + p#pi^{+}#slash{#pi}^{-} + p#slash{#pi}^{+}#pi^{-} + #slash{p}#pi^{+}#pi^{-}"
                        };
const TString XSECTITLE = "#sigma^{ #gamma^{*}p #rightarrow p#pi^{+}#pi^{-} } [#mub]";

/* private-type variables */
//Depending on Q2W_BNG(1/2) & xsectype (AT/VM),
//setup following
TString _q2w_bng;
TString _dtype;
TString _topNames[5];
TFile* _fyexp[5];
TFile* _fysim[5];
//since there is only 1 Q2 bin, the following can be directly obtained
float _q2min;
float _q2max;
float _dq2;

/* public-type functions */
void plotxsec(TString xsectype = "vm", bool ploty=kFALSE, bool sim=kFALSE,bool comp=kFALSE);
void plotcompxsec(TString xsectype = "vm", bool ploty=kFALSE, bool comp=kFALSE);
void plotxsec_CommonBins(seq_t seq=ACC_CORR, bool sim=kFALSE, int Q2Wbin=0);

/* private-type functions */
Bool_t setup(TString xsectype = "vm", TString q2w_bng="", bool pol=false);
TH1F* normalizeYield(TH1F* hYW);

//! each plot function gets a main canvas "c" and tnail canvas "ct" via
//  c->AddExec("tnail", "tnail()")
TCanvas* c;
TCanvas* ct;
TPad* selold_tnail;
/************** for Luminosity & vgflux ********************/
double nu(double w, double q2) {
    return (w*w-MP*MP+q2)/(2*MP);
}

double epsilon(double w, double q2) {
  double n = nu(w,q2);
  double e0 = E1F_E0;
  double e1 = e0-n;
  double epsInv = 1+2*(q2+n*n)/(4*e0*e1-q2);
    return 1.0/epsInv;
}

double getvgflux(double w, double q2, double e0 = E1F_E0) {
  double eps = epsilon(w,q2);
  return A*w*(w*w-MP*MP)/(4*PI*e0*e0*MP*MP*q2*(1-eps));
}
/**************************/