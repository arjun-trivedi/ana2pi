#ifndef PROCMORAND_H
#define PROCMORAND_H

#include "ep_processor.h" // Base class: EpProcessor

class ProcMorand : public EpProcessor {

public:
	ProcMorand(TDirectory *td);
	virtual ~ProcMorand();
	virtual void handle(DataOmega* d);
};

#endif // PROCMORAND_H
