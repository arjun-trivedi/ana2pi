#ifndef DATA3PI_H
#define DATA3PI_H

#include <TROOT.h>

class Data3pi
{

public:
	Data3pi();
	virtual ~Data3pi();
	void Clear();
	
	Float_t mm2p;
	Float_t mm2ppip;
	Float_t mm2ppim;
	Float_t mm2ppippim;
	Float_t s;
	Float_t t;   //morand
	Float_t t0;
	Float_t u;
	Float_t W;
	Float_t Q2;  //morand
	Float_t nu;
	Float_t xb;  //morand
	Float_t phi; //morand
	Float_t theta;
};

#endif // DATA3PI_H
