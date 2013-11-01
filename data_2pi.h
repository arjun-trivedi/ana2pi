#ifndef DATA2PI_H
#define DATA2PI_H

#include <TROOT.h>

class Data2pi
{

public:
	Data2pi();
	virtual ~Data2pi();
	void Clear();
	
	Float_t mm2ppippim;
	Float_t mm2ppip;
	Float_t mm2ppim;
	Float_t mm2pippim;
	Float_t W;
	Float_t Q2;  //morand
	Char_t  h;
	struct vars{
		Float_t M1;
	    Float_t M2;
	    Float_t theta;
	    Float_t phi;
	    Float_t alpha;
    } asgnmt1, asgnmt2, asgnmt3;
};

#endif // DATA2PI_H
