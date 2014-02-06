#ifndef DATA2PI_H
#define DATA2PI_H

#include <TROOT.h>

class Data2pi
{

public:
	Data2pi();
	virtual ~Data2pi();
	void Clear();
	
	Float_t Q2; 
	Float_t W;
	Int_t top;
	Char_t  h;
	Float_t mm2ppippim;
    Float_t mmppippim;
	Float_t mm2ppip;
    Float_t mmppip;
	Float_t mm2ppim;
    Float_t mmppim;
	Float_t mm2pippim;
    Float_t mmpippim; 
    struct vars{
		Float_t M1;
	    Float_t M2;
	    Float_t theta;
	    Float_t phi;
	    Float_t alpha;
    } varset1, varset2, varset3;
};

#endif // DATA2PI_H
