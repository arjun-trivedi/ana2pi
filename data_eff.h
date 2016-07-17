#ifndef DATAEFF_H
#define DATAEFF_H

#include <TROOT.h>

class DataEff
{

public:
	DataEff();
	virtual ~DataEff();
	void Clear();
	Int_t sector_e;
	Int_t sector_p;
	Int_t sector_pip;
	Float_t theta_e;
	Float_t theta_p;
	Float_t theta_pip;
	Float_t p_e;
	Float_t p_p;
	Float_t p_pip;
};

#endif // DATAEFF_H
