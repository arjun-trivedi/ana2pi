#include "TF1.h"
using namespace TMath;

float a0p[]={24.,24.,23.,23.5,24.5,24.5};
float a0m[]={25.,26.,26.,25.5,27.,26.};
float a1p[]={0.22,0.23,0.2,0.2,0.22,0.22};
float a1m[]={0.22,0.22,0.22,0.22,0.16,0.16};
float a2p[]={8.,8.,8.,8.,8.,8.};
float a2m[]={8.,8.,8.,8.,8.,8.};
float a3p[]={1.,1.,1.,1.,1.,1.};
float a3m[]={1.,1.,1.,1.,1.,1.};


Double_t dPhiFid_e16_hdrn_h(Double_t *x, Double_t *parms) {
        Double_t theta = x[0];
	Int_t s=parms[0]-1;
        Double_t phib = a0p[s]*(1-Exp(-a1p[s]*(theta-a2p[s])))+a3p[s];
        return phib;
}

Double_t dPhiFid_e16_hdrn_l(Double_t *x, Double_t *parms) {
        Double_t theta = x[0];
        Int_t s=parms[0]-1;
        Double_t phib = -a0m[s]*(1-Exp(-a1m[s]*(theta-a2m[s])))+a3m[s];
        return phib;
}

TF1* fPhiFid_e16_hdrn_h(Int_t sector) {
        TF1* retFunc = new TF1(TString::Format("fphifid_e16_hdrn_l"),dPhiFid_e16_hdrn_h,0,70,1);
        retFunc->SetParameter(0,sector);
        return retFunc;
}

TF1* fPhiFid_e16_hdrn_l(Int_t sector) {
        TF1* retFunc = new TF1(TString::Format("fphifid_e16_hdrn_l"),dPhiFid_e16_hdrn_l,0,70,1);
        retFunc->SetParameter(0,sector);
        return retFunc;
}

//! mod functions
Double_t dPhiFid_e16_hdrn_h_mod(Double_t *x, Double_t *parms) {
        Double_t theta = x[0];
        Int_t s=parms[0]-1;
        Double_t phib = a0p[s]*(1-Exp(-a1p[s]*(theta-a2p[s])))+a3p[s];
	//! trivedia modified 
        float offst[]={0,60,120,180,240,300}; //offst to get -[30,30] to values appropriate per sector
        phib+=offst[s];
        //
        return phib;
}

Double_t dPhiFid_e16_hdrn_l_mod(Double_t *x, Double_t *parms) {
        Double_t theta = x[0];
        Int_t s=parms[0]-1;
        Double_t phib = -a0m[s]*(1-Exp(-a1m[s]*(theta-a2m[s])))+a3m[s];
	//! trivedia modified 
        float offst[]={0,60,120,180,240,300}; //offst to get -[30,30] to values appropriate per sector
        phib+=offst[s];
        //
        return phib;
}

TF1* fPhiFid_e16_hdrn_h_mod(Int_t sector) {
        TF1* retFunc = new TF1(TString::Format("fphifid_e16_hdrn_l_mod"),dPhiFid_e16_hdrn_h_mod,0,70,1);
        retFunc->SetParameter(0,sector);
        return retFunc;
}

TF1* fPhiFid_e16_hdrn_l_mod(Int_t sector) {
        TF1* retFunc = new TF1(TString::Format("fphifid_e16_hdrn_l_mod"),dPhiFid_e16_hdrn_l_mod,0,70,1);
        retFunc->SetParameter(0,sector);
        return retFunc;
}


