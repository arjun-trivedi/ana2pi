#include "TF1.h"
#include "TMath.h"
//using namespace TMath;

/*
[05-06-16]
+ Direct C++ implemention of function 'fidu_e_sub' from $ANA2PI/cut_fid_e16.f, with the following modifications:
	1. EI's 'del_phie'([0,30]) -> +/-'del_phie'([-30,30])
	2. +/-'del_phie'([-30,30]) -> limits appropriate as per sector
*/
Double_t pshift=0.14;
Double_t c1=12.0, c2=18.5, c3=0.25, c4=25.0;
Double_t factor=0.416667;

//! AT-MOD vars
Double_t offst[]={0,60,120,180,240,300}; //offst to get -[30,30] to values appropriate per sector
Double_t pi=TMath::Pi();

Double_t phi(Double_t *x, Double_t *parms,Int_t sign) {
        Double_t thetael=x[0];

        Double_t pel=parms[0];
        Int_t isctr=parms[1]-1;
        printf("inputs: theta=%.2f,p=%.2f,isctr=%d,offset=%.2f\n",thetael,pel,isctr,offst[isctr]);
        printf("inputs: pshift,c1,c2,c3,c4,factor=%.2f,%.2f,%.2f,%.2f,%.2f,%.6f\n",pshift,c1,c2,c3,c4,factor);

        Double_t thetacut=(c1+c2/((pel+pshift)));
        printf("thetacut=%.2f\n",thetacut);

        Double_t phi=offst[isctr];
        if (thetael>=thetacut){
          //Double_t expon=c3*(TMath::Power(pel,factor));
          //phih=c4*TMath::Power(TMath::Sin((thetael-thetacut)*pi/180.),expon);
          Double_t expon=c3*(pel**factor);
          phi=sign*c4*TMath::Sin((thetael-thetacut)*pi/180.)**expon;
          //printf("expon,phi=%.2f,%.2f\n",expon,phi);
          //! AT-MOD 
          phi+=offst[isctr];
          //printf("phioffset=%.2f\n",phi);
        }
        return phi;
}

Double_t phi_h(Double_t *x, Double_t *parms) {
	return phi(x,parms,1);
}
Double_t phi_l(Double_t *x, Double_t *parms) {
        return phi(x,parms,-1);
}

TF1* fphi_h(Double_t p,Int_t sector) {
        TF1* retFunc = new TF1(TString::Format("phi_h"),phi_h,0,40,2);
	retFunc->SetParameter(0,p);
        retFunc->SetParameter(1,sector);
        return retFunc;
}
TF1* fphi_l(Double_t p,Int_t sector) {
        TF1* retFunc = new TF1(TString::Format("phi_l"),phi_l,10,40,2);
        retFunc->SetParameter(0,p);
        retFunc->SetParameter(1,sector);
        return retFunc;
}




