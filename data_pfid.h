#ifndef DATAPFID_H
#define DATAPFID_H

#include <TROOT.h>

class DataPFid
{

public:
	DataPFid();
	virtual ~DataPFid();
	void Clear();
	Int_t sector_p;
	Int_t sector_pip;
	Int_t sector_pim;
	Float_t theta_p;
	Float_t theta_pip;
	Float_t theta_pim;
	Float_t phi_p;
	Float_t phi_pip;
	Float_t phi_pim;
	Float_t p_p;
	Float_t p_pip;
	Float_t p_pim;
};

#endif // DATAFID_H
