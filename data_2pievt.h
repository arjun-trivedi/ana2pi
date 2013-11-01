#ifndef DATA2PIEVT_H
#define DATA2PIEVT_H

#include <TROOT.h>

class Data2piEvt
{

public:
	Data2piEvt();
	virtual ~Data2piEvt();
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

#endif // DATA2PIEVT_H
