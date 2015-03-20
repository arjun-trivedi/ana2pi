#ifndef DATAPEFF_H
#define DATAPEFF_H

#include <TROOT.h>

class DataPEff
{

public:
	DataPEff();
	virtual ~DataPEff();
	void Clear();
	Int_t sector_p;
	Int_t sector_pip;
	Int_t sector_pim;
	Float_t theta_p;
	Float_t theta_pip;
	Float_t theta_pim;
	Float_t p_p;
	Float_t p_pip;
	Float_t p_pim;
};

#endif // DATAFID_H
