#ifndef DATASKIMQ_H
#define DATASKIMQ_H

#include <TROOT.h>

class DataSkimQ
{
public:
	DataSkimQ();
	virtual ~DataSkimQ();
	void Clear();
		Bool_t isEVT_ETGT_2POS_ETGT_1NEG;
        Bool_t isEVT_ETGT_2POS;
        Bool_t isEVT_ETGT_1POS_ETGT_1NEG;
};

#endif // DATASKIMQ_H
