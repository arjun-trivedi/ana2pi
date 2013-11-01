#ifndef DATAEVT_H
#define DATAEVT_H

#include <TROOT.h>

class DataEvt
{

public:
	DataEvt();
	virtual ~DataEvt();
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
    } varset1, varset2, varset3;
};

#endif // DATAEVT_H
