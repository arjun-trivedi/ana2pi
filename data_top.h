#ifndef DATATOP_H
#define DATATOP_H

#include <TROOT.h>

class DataTop
{

public:
	DataTop();
	virtual ~DataTop();
	void Clear();
	
	Float_t mm2ppippim;
    Float_t mmppippim;
	Float_t mm2ppip;
    Float_t mmppip;
	Float_t mm2ppim;
    Float_t mmppim;
	Float_t mm2pippim;
    Float_t mmpippim; 
	Float_t W;
	Float_t Q2;  //morand
	Char_t  h;
	struct vars{
		Float_t M1;
	    Float_t M2;
	    Float_t theta;
	    Float_t phi;
	    Float_t alpha;
    } varset1, varset2, varset3;
};

#endif // DATATOP_H
