#ifndef DATAEFID_H
#define DATAEFID_H

#include <TROOT.h>

class DataEFid
{

public:
	DataEFid();
	virtual ~DataEFid();
	void Clear();
	Bool_t fidE;
        Int_t sector;
    Float_t p;
	Float_t phi;
	Float_t theta;
        Float_t dc_xsc;
        Float_t dc_ysc;
        Float_t ech_x;
        Float_t ech_y;
};

#endif // DATAFID_H
