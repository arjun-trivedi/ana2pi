#ifndef CUTEFF_H
#define CUTEFF_H

#include <TF1.h>
#include <TString.h>

using namespace std;
/*
/* [03-20-15]
/* + This Class implements cutting away e,p,pip and pim if they pass
/*   through the inefficient regions of the CLAS detector in the E1F run. 
*/
class CutEff {

public:
	CutEff();
	~CutEff();
	
	Bool_t InEfficientRegion(Int_t pid,Int_t sector,Float_t theta,Float_t p);

private:
	TString _prtcl[4];
	Int_t _theta_min[4];
	Int_t _theta_max[4];
	Int_t _p_min[4];
	Int_t _p_max[4];

	TF1** _cut_lw[4][6];
	TF1** _cut_hg[4][6];

	void setup_cuts();
};

#endif // CUTEFF_H
