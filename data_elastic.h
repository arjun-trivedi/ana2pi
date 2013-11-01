#ifndef DATAELASTIC_H
#define DATAELASTIC_H

#include <TROOT.h>

class DataElastic
{

public:
	DataElastic();
	virtual ~DataElastic();
	void Clear();
	
	Float_t mm2e;
	Float_t theta;
	Float_t phi;
};

#endif // DATAELASTIC_H
