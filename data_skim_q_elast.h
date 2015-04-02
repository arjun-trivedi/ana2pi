#ifndef DATASKIMQELAST_H
#define DATASKIMQELAST_H

#include <TROOT.h>

class DataSkimQElast
{
public:
	DataSkimQElast();
	virtual ~DataSkimQElast();
	void Clear();
    Bool_t isEVT_EQGT_1POS;
};

#endif // DATASKIMQELAST_H
