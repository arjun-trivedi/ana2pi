#ifndef DATAEEFF_H
#define DATAEEFF_H

#include <TROOT.h>

class DataEEff
{

public:
	DataEEff();
	virtual ~DataEEff();
	void Clear();
	Int_t sector;
	Float_t theta;
	Float_t p;
};

#endif // DATAFID_H
