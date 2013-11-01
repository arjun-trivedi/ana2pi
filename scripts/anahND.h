#include <THnSparse.h>
#include <TString.h>
#include <TH1.h>
#include <TCanvas.h>

#include "constants.h"

TString _type;

const Int_t nDtype = 2;
enum dtype_t{EXP, SIM};
TString dtypeTitle[nSeq] = {"EXP", "SIM"};

THnSparse* h5D[nDtype][nSeq]; 
TH1D*      h1D[nDtype][nSeq];
TH1D*      h1D0[nDtype][nSeq]; //!empty bin visualization only for SIM h1Ds 
	
TString h5Dpath[nSeq];
Int_t nBins_h5D[nDtype][nSeq] = {0, 0, 0, 0,
	                             0, 0, 0, 0};
Int_t nFBins; 

TString Q2Wdirname;	
	
void tnail_dcStudy(TCanvas* c_dcStudy, TCanvas* c_tnail_dcStudy);
void makeHists_h1D();

void FillLinHist(TH1* h1, THnSparse* hN, Int_t iFBin, Int_t* FBincoord);
Int_t GetNum0Bins(TH1* h1);
