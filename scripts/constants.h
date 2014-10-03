#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <TString.h>

const Int_t nTOP          = 4;
const Int_t nVARSET       = 1; //atrivedi:080213
const Int_t nSEQ          = 6; 
const Int_t nHEL          = 5; 
const Int_t nEVAR         = 3;
const Int_t nVAR          = 5;


//! hY8D bins !//
//bin index: {0, 1,  2, 3,  4,  5,     6,   7}
//           {h, Q2, W, M1, M2, theta, phi, alpha}

const Int_t PROJDIMS[5] = {3,4,5,6,7}; //to hY8D -> hY5D

enum seq_t{
	TH=0, RECO=1, ACC=2, ACC_CORR=3, HOLE=4, FULL=5
};
enum hel_t{ //enum helicities NEG, ZERO and POS as per THnSparse-->h variable axis. UNPOL = 0 && ASYM = 4  
	UNPOL=0, NEG=1, ZERO=2, POS=3, ASYM=4
};
enum evar_t{ 
	H=0, Q2=1, W=2
};
enum var_t{ 
	M1=0, M2=1, THETA=2, PHI=3, ALPHA=4
};


TString evarName[nEVAR]  = {"h", "Q2", "W"};
TString varName[nVAR]    = {"M1", "M2", "theta", "phi", "alpha"};  
TString evarTitle[nEVAR] = {"Helicity", "Q^{2}", "W"};
TString varTitle[nVARSET][nVAR] = { 
                                    "M_{p#pi^{+}}", "M_{#pi^{+}#pi^{-}}",
                                    "#theta_{#pi^{-}}", "#phi_{#pi^{-}}",
                                    "#alpha_{[p^{'}#pi^{+}][p#pi^{-}]}"};
                                    //atrivedi:080213
                                    /*"M_{p#pi^{+}}", "M_{#pi^{+}#pi^{-}}",
                                    "#theta_{p}", "#phi_{p}",
                                    "#alpha_{[#pi^{+}#pi^{-}][pp^{'}]}",
                                    "M_{p#pi^{+}}", "M_{p#pi^{-}}",
                                    "#theta_{#pi^{+}}", "#phi_{#pi^{+}}",
                                    "#alpha_{[p^{'}#pi^{-}][p#pi^{+}]}"
                                };*/
TString evarUnitName[nEVAR]   = {"", "(GeV^{2})", "(GeV)"};                                 
TString varUnitName[nVAR]     = {"[GeV]", "[GeV]", "[#degree]", "[#degree]", "[#degree]"}; 
TString helTitle[nHEL]        = {"UNPOL", "NEG", "0", "POS", "ASYM"};
TString evtSel2PiName[nTOP]         = {"p#pi^{+}#pi^{-}", "p#pi^{+}", "p#pi^{-}", "#pi^{+}#pi^{-}" };
TString asgnmtName[nVARSET]   = {"Varset1"};//, "Varset2", "Varset3"}; #atrivedi:080213
TString asgnmtTitle[nVARSET]  = {"#Delta^{++}"};//, "#rho", "#Delta^{0}"}; atrivedi:080213
TString seqTitle[nSEQ]      = {"TH", "RECO", "ACC", "ACC_CORR", "HOLE", "FULL"};

const float LUM = 19.844;

static const float PI = 3.14159265358979312;
static const float E1F_E0 = 5.499;
static const float FSC = 0.00729735253;
static const float A = FSC;
static const float NA = 6.02214129E23;
static const float QE = 1.60217646E-19;
static const float MP = 0.93827203;

/*double nu(double w, double q2) {
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
}*/

#endif /* CONSTANTS_H_ */
