#ifndef DATAEID_H
#define DATAEID_H

#include <TROOT.h>

class DataEid
{
public:
	DataEid();
	virtual ~DataEid();
	void Clear();
	Int_t sector;
        Float_t beta;
        Float_t betaStrE;
	Float_t dtE;
	Float_t p;
	//from EC
	Float_t ec_ei;
	Float_t ec_eo;
	Float_t etot;
	//from CC
	Int_t nphe;
	Int_t cc_segm;
	Float_t cc_theta;
};

#endif // DATAEID_H
