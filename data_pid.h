#ifndef DATAPID_H
#define DATAPID_H

#include <TROOT.h>

class DataPid
{
public:
	DataPid();
	virtual ~DataPid();
	void Clear();
        Int_t h10IdxP;
        Int_t h10IdxPip;
        Int_t h10IdxPim;
        Int_t sectorP;
        Int_t sectorPip;
        Int_t sectorPim;
        Float_t betaP;
        Float_t betaPip;
        Float_t betaPim;
	Float_t pP;
	Float_t pPip;
	Float_t pPim;
        Float_t betaStrP;
	Float_t betaStrPip;
	Float_t betaStrPim;
        Float_t dtP;
	Float_t dtPip;
	Float_t dtPim;
        /*CC, EC signature*/
        //from EC
        Float_t P_ec_ei;
        Float_t P_ec_eo;
        Float_t P_etot;
        Float_t Pip_ec_ei;
        Float_t Pip_ec_eo;
        Float_t Pip_etot;
        Float_t Pim_ec_ei;
        Float_t Pim_ec_eo;
        Float_t Pim_etot;        
        //from CC
        Int_t P_nphe;
        Int_t P_cc_segm;
        Float_t P_cc_theta;  
        Int_t Pip_nphe;
        Int_t Pip_cc_segm;
        Float_t Pip_cc_theta;
        Int_t Pim_nphe;
        Int_t Pim_cc_segm;
        Float_t Pim_cc_theta;
};

#endif 
