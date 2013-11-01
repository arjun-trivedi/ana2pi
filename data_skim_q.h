#ifndef DATASKIMQ_H
#define DATASKIMQ_H

#include <TROOT.h>

class DataSkimQ
{
public:
	DataSkimQ();
	virtual ~DataSkimQ();
	void Clear();
        Bool_t isEVT_2POS_EX;
        Bool_t isEVT_1POS1NEG_EX;
        Bool_t isEVT_2POS1NEG_EX;
};

#endif // DATASKIMQ_H
