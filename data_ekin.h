#ifndef DATAEKIN_H
#define DATAEKIN_H

#include <TROOT.h>

class DataEkin
{

public:
	DataEkin();
	virtual ~DataEkin();
	void Clear();
	Int_t sector;
    Float_t W;
	Float_t Q2;
	Float_t nu;
	Float_t xb;
	Float_t E1;
	Float_t theta1;
	Float_t phi1;
	Float_t theta;
	Float_t phi;
	Float_t vx;
	Float_t vy;
	Float_t vz;
};

#endif // DATAEKIN_H
