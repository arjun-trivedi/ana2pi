#ifndef DATAPFIDELAST_H
#define DATAPFIDELAST_H

#include <TROOT.h>

class DataPFidElast
{

public:
	DataPFidElast();
	virtual ~DataPFidElast();
	void Clear();
	Int_t sector_p;
	Float_t theta_p;
	Float_t phi_p;
	Float_t p_p;
};

#endif // DATAFID_H
