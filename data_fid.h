#ifndef DATAFID_H
#define DATAFID_H

#include <TROOT.h>

class DataFid
{

public:
	DataFid();
	virtual ~DataFid();
	void Clear();
	Bool_t fidE;
	Bool_t fidP;
	Bool_t fidPip;
	Bool_t fidPim;
};

#endif // DATAFID_H
