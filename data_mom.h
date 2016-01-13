#ifndef DATAMOM_H
#define DATAMOM_H

#include <TROOT.h>

/*
[01-10-16]
Added variables for pcorr:{pip,p}
*/
class DataMom
{
public:
	DataMom();
	virtual ~DataMom();
	void Clear();
	//!momentum corrections for electron
	Int_t sector;
	Float_t p;
	Float_t dcx, dcy, dcz, dp; 
	//!momentum corrections for pip
	Int_t sector_pip;
	Float_t p_pip;
	Float_t dcx_pip, dcy_pip, dcz_pip, dp_pip;
	//!momentum corrections for p
	Int_t sector_p;
	Float_t p_p;
	Float_t dcx_p, dcy_p, dcz_p, dp_p;
};

#endif // DATAMOM_H
