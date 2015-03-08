#ifndef DATAPIDELAST_H
#define DATAPIDELAST_H

#include <TROOT.h>

class DataPidElast
{
public:
	DataPidElast();
	virtual ~DataPidElast();
	void Clear();
        Int_t h10IdxP;
        Int_t sectorP;
        Float_t betaP;
        Float_t pP;
	Float_t betaStrP;
	Float_t dtP;
	/*CC, EC signature*/
        //from EC
        Float_t P_ec_ei;
        Float_t P_ec_eo;
        Float_t P_etot;
        //from CC
        Int_t P_nphe;
        Int_t P_cc_segm;
        Float_t P_cc_theta;  
};

#endif 
