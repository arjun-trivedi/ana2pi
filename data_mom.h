#ifndef DATAMOM_H
#define DATAMOM_H

#include <TROOT.h>

class DataMom
{
public:
	DataMom();
	virtual ~DataMom();
	void Clear();
        Int_t sector;
        Float_t p;
	Float_t dcx, dcy, dcz, dp; //momentum corrections for electron
};

#endif // DATAMOM_H
